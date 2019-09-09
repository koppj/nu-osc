// ----------------------------------------------------------------------------
// Supplement to snu.c providing an interface to nuSQuIDS
// ----------------------------------------------------------------------------
// Author: Joachim Kopp (JGU Mainz and CERN, jkopp@uni-mainz.de)
// ----------------------------------------------------------------------------
#include <gsl/gsl_matrix.h>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/marray.h>
#include <nuSQuIDS/tools.h>
#include <nusquids_decay/nusquids_decay.h>
#include "osc_decay/osc_decay.h"
#include "snu.h"

#ifdef NU_USE_NUSQUIDS

//#define NUSQUIDS_DEBUG

// Constants
#define GLB_Ne_MANTLE       0.5        // Effective electron numbers for calculation
#define GLB_Ne_CORE         0.468      //   of MSW potentials

// Macros
#define SQR(x)      ((x)*(x))          // x^2

// Parameter structure for Ivan's code
regeneration::Param p_oscdecay; // Parameter vector for Ivan's code

using namespace nusquids;

// ----------------------------------------------------------------------------
//         N u S Q u I D S - B A S E D   P R O B A B I L I T I E S
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
int snu_probability_matrix_nusquids_internal(double P[][2][SNU_MAX_FLAVORS],
      unsigned n_E, double *E, double ini_state_nu[][3], double ini_state_nubar[][3],
      int psteps, const double *length, const double *density,
      unsigned n_flavors, double th[SNU_MAX_FLAVORS+1][SNU_MAX_FLAVORS+1],
      double delta[SNU_MAX_PHASES], double dmsq[SNU_MAX_FLAVORS-1],
      double M_A_prime, double g_prime)
// ----------------------------------------------------------------------------
// Calculates the neutrino oscillation probability matrix using nuSQuIDS. As
// nuSQuIDS simulates also processes like tau regeneration or neutrino decay,
// there can be migration between energy bins, hence the function computes
// the probability matrix for a whole range of energy bins at once.
// ----------------------------------------------------------------------------
// Parameters:
//   P:         Buffer for the storage of the final weights. The indices of P
//              correspond to energy, cp sign (nu/nubar), final flavor.
//   cp_sign:   +1 for neutrinos, -1 for antineutrinos
//   n_E:       number of energy bins
//   E:         list of neutrino energies (in GeV)
//   ini_state_nu/nubar: initial state for neutrinos/antineutrinos. Each entry
//              corresponds to a particular energy bin and gives the weights of
//              the different flavors at that energy
//   psteps:    Number of layers in the matter density profile
//   length:    Lengths of the layers in the matter density profile in km
//   density:   The matter densities in g/cm^3
//   n_flavors: number of neutrino flavors considered
//   th:        mixing angles
//   delta:     CP violating phases
//   dmsq:      mass squared differences
//   M_A_prime: mass of dark force carrier
//   g_prime:   coupling of A' to sterile neutrinos
// ----------------------------------------------------------------------------
{
  squids::Const units;

  // Create coupling matrix
  gsl_matrix *couplings = gsl_matrix_calloc(n_flavors, n_flavors);
  gsl_matrix_set(couplings, 3, 2, g_prime);  // g_{43} := 1
  gsl_matrix_set(couplings, 3, 1, g_prime);  // g_{42} := 1
  gsl_matrix_set(couplings, 3, 0, g_prime);  // g_{41} := 1
     // FIXME FIXME Need to reflect correct flavor structure of couplings here!

  unsigned dk = 1; //FIXME
  unsigned n_red = (n_E - 1) / dk + 1;

  // Create nuSQuIDS objects
  std::vector<double> m_nu(n_flavors);
  std::vector<double> x(psteps+1);
  std::vector<double> rho(psteps+1);
  std::vector<double> ye(psteps+1);
  marray<double,1>    bin_centers({n_red});
  marray<double,3>    ini_state({n_red, 2, n_flavors});
  m_nu[0] = 0.0;
  for (unsigned i=1; i < n_flavors; i++)
    m_nu[i] = sqrt(dmsq[i-1]);
  for (unsigned k=0; k < n_red; k++)
  {
    unsigned m = k*dk + dk/2;
    if (m > n_E - 1)
      m = n_E - 1;
    bin_centers[k] = E[m] * units.GeV;
    for (unsigned j=0; j < 3; j++)
    {
      ini_state[k][0][j] = ini_state_nu[m][j];
      ini_state[k][1][j] = ini_state_nubar[m][j];
    }
    for (unsigned j=3; j < n_flavors; j++)
      ini_state[k][0][j] = ini_state[k][1][j] = 0.0;
  }
  x[0] = 0.0;
  for (int i=0; i < psteps; i++)
  {
    x[i+1] = x[i] + length[i];
    rho[i] = density[i];
    ye[i]  = GLB_Ne_MANTLE;
  }
  rho[psteps] = density[psteps-1];
  ye[psteps]  = ye[psteps-1];
 
  nuSQUIDSDecay nusquids(bin_centers,  // Array of energy bins
                         n_flavors,    // number of flavor/mass eigenstates
                         both,         // simulate both nu and nubar
                         false,        // incoherent interactions
                         true,         // decay and regeneration
                         false,        // scalar (no pseudoscalar) couplings
                         m_nu,         // vector of neutrino masses
                         couplings);   // coupling matrix
  nusquids.Set_TauRegeneration(false);

#ifdef NUSQUIDS_DEBUG
  printf("# Initial state:\n");
  for (unsigned k=0; k < n_red; k++)
  {
    printf("#   E=%7.5g -- nu = (", bin_centers[k] / units.GeV);
    for (unsigned j=0; j < n_flavors; j++)
      printf(" %10.7g", ini_state[k][0][j]);
    printf(" ), nubar = (");
    for (unsigned j=0; j < n_flavors; j++)
      printf(" %10.7g", ini_state[k][1][j]);
    printf(" )\n");
  }
#endif

  // Set up mixing parameters
  for (unsigned i=0; i < n_flavors; i++)
    for (unsigned j=i+1; j < n_flavors; j++)
      nusquids.Set_MixingAngle(i, j, th[i+1][j+1]);
  for (unsigned i=1; i < n_flavors; i++)
    nusquids.Set_SquareMassDifference(i, dmsq[i-1]);

  // FIXME include CP phases!

  // Set up matter density profile
  if (psteps < 2)
  {
    nusquids.Set_Body(std::make_shared<ConstantDensity>(rho[0], ye[0]));
    nusquids.Set_Track(std::make_shared<ConstantDensity::Track>(x[psteps]*units.km));
  }
  else
  {
    nusquids.Set_Body(std::make_shared<VariableDensity>(x, rho, ye));
    nusquids.Set_Track(std::make_shared<VariableDensity::Track>(x[psteps]*units.km));
  }

  // Set up numerical integration
  nusquids.Set_GSL_step(gsl_odeiv2_step_rkf45);
  nusquids.Set_rel_error(1e-5);
  nusquids.Set_abs_error(1e-5);
  nusquids.Set_initial_state(ini_state, flavor); // initial state given in flavor basis

  // Run nuSQuIDS
  nusquids.EvolveState();

#ifdef NUSQUIDS_DEBUG
  printf("# Final state:\n");
  for (unsigned k=0; k < n_red; k++)
  {
    printf("#   E=%7.5g -- nu = (", bin_centers[k] / units.GeV);
    for (unsigned j=0; j < n_flavors; j++)
      printf(" %10.7g", nusquids.EvalFlavor(j, bin_centers[k], 0));
    printf(" ), nubar = (");
    for (unsigned j=0; j < n_flavors; j++)
      printf(" %10.7g", nusquids.EvalFlavor(j, bin_centers[k], 1));
    printf(" )\n");
  }
#endif

  // Store output
  for (unsigned k=0; k < n_E; k++)
  {
    if (k % dk == 0)
    {
      for (unsigned j=0; j < n_flavors; j++)
      {
        P[k][0][j] = nusquids.EvalFlavor(j, bin_centers[k/dk], 0);
        P[k][1][j] = nusquids.EvalFlavor(j, bin_centers[k/dk], 1);
      }
    }
    else
    {
      for (unsigned j=0; j < n_flavors; j++)
      {
        P[k][0][j] = P[k-1][0][j];
        P[k][1][j] = P[k-1][1][j];
      }
    }

    for (unsigned j=n_flavors; j < SNU_MAX_FLAVORS; j++)
    {
      P[k][0][j] = P[k][1][j] = 0.0;
   }
  }

  // Cleanup
  if (couplings)  { gsl_matrix_free(couplings);  couplings = NULL; }

  return 0;
}


// ----------------------------------------------------------------------------
int snu_set_oscillation_parameters_osc_decay_internal(int n_flavors, glb_params params)
// ----------------------------------------------------------------------------
// Stores the oscillation parameters in the data structure expected by
// Ivan's oscillation and decay code
// ----------------------------------------------------------------------------
// Parameters:
//   n_flavors: number of neutrino flavors in snu oscillation engine
//   params:    GLoBES parameter vector
// ----------------------------------------------------------------------------
{
  int status = 0;

  p_oscdecay.snq12 = SQR(sin(glbGetOscParamByName(params, "TH12")));
  p_oscdecay.snq13 = SQR(sin(glbGetOscParamByName(params, "TH13")));
  p_oscdecay.snq23 = SQR(sin(glbGetOscParamByName(params, "TH23")));
  p_oscdecay.dlt13 = glbGetOscParamByName(params, "DELTA_0");
  if (n_flavors <= 3)
  { 
    p_oscdecay.snq14 = p_oscdecay.snq24 = p_oscdecay.snq34 = p_oscdecay.dmq41 = 0.0;
    if (glbGetOscParamByName(params, "DM31") > 0.)  // Normal ordering
      p_oscdecay.m_nu = { 0.,
                 sqrt(glbGetOscParamByName(params, "DM21")),
                 sqrt(glbGetOscParamByName(params, "DM31")) };
    else                                            // Inverted ordering
      p_oscdecay.m_nu = { 0.,
                 sqrt(-glbGetOscParamByName(params, "DM31")),
                 sqrt(-glbGetOscParamByName(params, "DM31") + 
                       glbGetOscParamByName(params, "DM21")) };
  }
  else if (n_flavors == 4)
  {
    p_oscdecay.snq14 = SQR(sin(glbGetOscParamByName(params, "TH14")));
    p_oscdecay.snq24 = SQR(sin(glbGetOscParamByName(params, "TH24")));
    p_oscdecay.snq34 = SQR(sin(glbGetOscParamByName(params, "TH34")));
    p_oscdecay.dmq41 = glbGetOscParamByName(params, "DM41");
    if (glbGetOscParamByName(params, "DM31") > 0.)  // Normal ordering
      p_oscdecay.m_nu = { 0.,
                 sqrt(glbGetOscParamByName(params, "DM21")),
                 sqrt(glbGetOscParamByName(params, "DM31")),
                 sqrt(glbGetOscParamByName(params, "DM41")) };
    else                                            // Inverted ordering
      p_oscdecay.m_nu = { 0.,
                 sqrt(-glbGetOscParamByName(params, "DM31")),
                 sqrt(-glbGetOscParamByName(params, "DM31") + 
                       glbGetOscParamByName(params, "DM21")),
                 sqrt(-glbGetOscParamByName(params, "DM31") + 
                       glbGetOscParamByName(params, "DM41")) };
    // FIXME include extra complex phases
  }
  else   // > 4 flavors not supported yet
    status = -1000;

  p_oscdecay.m_A = glbGetOscParamByName(params, "M_A_PRIME");
  p_oscdecay.g   = glbGetOscParamByName(params, "G_PRIME");
  if (n_flavors > 3)
  {
    if (fabs(glbGetOscParamByName(params, "MA_OVER_M4")) > 1e-30)
    {
      if (fabs(p_oscdecay.m_A) > 1e-30)
      {
        fprintf(stderr, "snu_set_oscillation_parameters_osc_decay_internal: "
                "Warning: Both M_A_PRIME and MA_OVER_M4 != 0. Ignoring MA_OVER_M4.\n");
        status = -1001;
      }
      else
        p_oscdecay.m_A = p_oscdecay.m_nu[3] * glbGetOscParamByName(params, "MA_OVER_M4");
    }

    if (fabs(glbGetOscParamByName(params, "M4_GAMMA")) > 1e-12)
    {
      if (p_oscdecay.g > 1e-12)
      {
        fprintf(stderr, "snu_set_oscillation_parameters_osc_decay_internal: "
                "Warning: Both G_PRIME and M4_GAMMA given. Ignoring M4_GAMMA.\n");
        status = -1002;
      }
      else
        p_oscdecay.g = p_oscdecay.get_g(glbGetOscParamByName(params, "M4_GAMMA"));
    }
  } // if (n_flavors > 3)
 
  return status;
}


// ----------------------------------------------------------------------------
int snu_get_oscillation_parameters_osc_decay_internal(int n_flavors, glb_params params)
// ----------------------------------------------------------------------------
// Stores the parameters relevant for decay in the given GLoBES parameter
// vector. Note: all other entries of this vector are left unaffected!
// ----------------------------------------------------------------------------
// Parameters:
//   n_flavors: number of neutrino flavors in snu oscillation engine
//   params:    GLoBES parameter vector
// ----------------------------------------------------------------------------
{
  glbSetOscParamByName(params, p_oscdecay.m_A, "M_A_PRIME");
  glbSetOscParamByName(params, p_oscdecay.g,   "G_PRIME");
  if (n_flavors > 3  &&  p_oscdecay.m_nu.size() == (unsigned) n_flavors)
  {
    glbSetOscParamByName(params, p_oscdecay.m_A / p_oscdecay.m_nu[3], "MA_OVER_M4");
    glbSetOscParamByName(params, snu_get_m4Gamma_osc_decay(), "M4_GAMMA");
  }
  else
  {
    glbSetOscParamByName(params, NAN, "MA_OVER_M4");
    glbSetOscParamByName(params, NAN, "M4_GAMMA");
  }

  return 0;
}


// ----------------------------------------------------------------------------
double snu_get_m4Gamma_osc_decay()
// ----------------------------------------------------------------------------
// Returns the sterile neutrino decay width, multiplied by the sterile neutrino
// mass, from Ivan's oscillation+decay code
// ----------------------------------------------------------------------------
{
  if (p_oscdecay.m_nu.size() > 3)
  {
    if (p_oscdecay.g > 0.)
    {
      double Gamma_tot = 0.0;
      for (int i=0; i < 3; i++)
        Gamma_tot += p_oscdecay.getGamma(i);
      return Gamma_tot * p_oscdecay.m_nu[3];
    }
    else
      return 0.0;
  }
  else
    return NAN;
}


// ----------------------------------------------------------------------------
int snu_probability_matrix_osc_decay_internal(double P[][2][SNU_MAX_FLAVORS],
      unsigned n_E, double *E, double ini_state_nu[][3], double ini_state_nubar[][3],
      int psteps, const double *length, const double *density,
      unsigned n_flavors, const double filter_value)
// ----------------------------------------------------------------------------
// Calculates the neutrino oscillation probability matrix using analytical
// approximations for oscillations + decay (based on Ivan Esteban's code).
// As there can be migration between energy bins, the function computes
// the probability matrix for a whole range of energy bins at once.
// ----------------------------------------------------------------------------
// Parameters:
//   P:         Buffer for the storage of the final weights. The indices of P
//              correspond to energy, cp sign (nu/nubar), final flavor.
//   cp_sign:   +1 for neutrinos, -1 for antineutrinos
//   n_E:       number of energy bins
//   E:         list of neutrino energies (in GeV)
//   ini_state_nu/nubar: initial state for neutrinos/antineutrinos. Each entry
//              corresponds to a particular energy bin and gives the weights of
//              the different flavors at that energy
//   psteps:    Number of layers in the matter density profile
//   length:    Lengths of the layers in the matter density profile in km
//   density:   The matter densities in g/cm^3
//   n_flavors: number of neutrino flavors considered
//   filter_value: width of low-pass filter for smoothing fast oscillations
//              (implemented as in GLoBES, see GLoBES manual for details)
// ----------------------------------------------------------------------------
{
  using namespace regeneration;

  // Create vectors of energy bins and initial state fluxes
  marray<double,1> bin_centers({n_E});
  marray<double,3> ini_state({n_E, 2, n_flavors});
  for (unsigned k=0; k < n_E; k++)
  {
    bin_centers[k] = E[k] * units::GeV;
    for (unsigned j=0; j < 3; j++)
    {
      ini_state[k][0][j] = ini_state_nu[k][j];
      ini_state[k][1][j] = ini_state_nubar[k][j];
    }
    for (unsigned j=3; j < n_flavors; j++)
      ini_state[k][0][j] = ini_state[k][1][j] = 0.0;
  }
  double L = 0.0;
  for (int i=0; i < psteps; i++)
    L += length[i];

  // Initialize oscillation engine
  regProb prob(bin_centers, ini_state, OSC_DECAY_MAJORANA);
  prob.setParam(p_oscdecay);

#ifdef NUSQUIDS_DEBUG
  printf("# Initial state:\n");
  for (unsigned k=0; k < n_E; k++)
  {
    printf("#   E=%7.5g -- nu = (", bin_centers[k] / units::GeV);
    for (unsigned j=0; j < n_flavors; j++)
      printf(" %10.7g", ini_state[k][0][j]);
    printf(" ), nubar = (");
    for (unsigned j=0; j < n_flavors; j++)
      printf(" %10.7g", ini_state[k][1][j]);
    printf(" )\n");
  }
#endif

  // Compute and store output
  marray<double,2> final_state({2, n_flavors});
  for (unsigned k=0; k < n_E; k++)
  {
    final_state = prob.getFinalFlux(bin_centers[k], L * units::km, filter_value * units::GeV);
    for (unsigned c=0; c < 2; c++)
    {
      for (unsigned j=0; j < n_flavors; j++)
        P[k][c][j] = final_state[c][j];
      for (unsigned j=n_flavors; j < SNU_MAX_FLAVORS; j++)
        P[k][c][j] = 0.0;
    }
  }

#ifdef NUSQUIDS_DEBUG
  printf("# Final state:\n");
  for (unsigned k=0; k < n_E; k++)
  {
    printf("#   E=%7.5g -- nu = (", bin_centers[k] / units::GeV);
    for (unsigned j=0; j < n_flavors; j++)
      printf(" %10.7g", P[k][0][j]);
    printf(" ), nubar = (");
    for (unsigned j=0; j < n_flavors; j++)
      printf(" %10.7g", P[k][0][j]);
    printf(" )\n");
  }
#endif

  return 0;
}

#endif // NU_USE_NUSQUIDS
