/***************************************************************************
 * Functions for MiniBooNE fit                                             *
 ***************************************************************************
 * Author: Joachim Kopp                                                    *
 ***************************************************************************/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <globes/globes.h>
#include "osc_decay/osc_decay.h"
#include "nu.h"

// Flags that affect chi^2 calculation
//#define OSC_NORM       // Take into account impact of \nu_\mu disapp. on normalization
//#define OSC_BG         // Include oscillations of \nu_e backgrounds
//#define OSC_NUMU       // Allow also the \nu_\mu sample to oscillate

// Data structures from Ivan's oscillation+decay code
#ifdef NU_USE_NUSQUIDS
  using namespace regeneration;
  static regProb prb_nu(OSC_DECAY_MAJORANA);
  static regProb prb_nu_bar(OSC_DECAY_MAJORANA);
  static regProb prb_bg_nue_from_mu(OSC_DECAY_MAJORANA);
  static regProb prb_bg_nue_from_K(OSC_DECAY_MAJORANA);
#endif

#define MB_DATA_DIR        "glb/mb-jk-2018/"
#define MB_COV_MATRIX_FILE "glb/mb-jk-2018/cov_matrix.h"

static const double MB_baseline = 0.520;   // [km] - MiniBooNE baseline
static const double SB_baseline = 0.100;   // [km] - SciBooNE baseline

// Energ binning
static const int E_true_bins   =  680;       // binning in E_true
static const double E_true_min =  120.0;     // [MeV]
static const double E_true_max = 3520.0;     // [MeV]
//static const double E_true_max = 5480.0;     // [MeV]

#define NE (E_reco_bins_e)
static const int E_reco_bins_e = 11;
static const double E_reco_bin_edges_e[] = {
                             200.0, 300.0,  375.0,  475.0,  550.0,  675.0,
                             800.0, 950.0, 1100.0, 1250.0, 1500.0, 3000.0 }; // [MeV]
static const double E_reco_min_e = E_reco_bin_edges_e[0];
static const double E_reco_max_e = E_reco_bin_edges_e[E_reco_bins_e];

#define NMU (E_reco_bins_mu)
static const int E_reco_bins_mu = 8;
static const double E_reco_bin_edges_mu[] = {    0.,   500.0,  700.0,  900.0, 1100.0, 
                                              1300.0, 1500.0, 1700.0, 1900.0 };
                        // fig. 17 from https://arxiv.org/abs/1805.12028
static const double E_reco_min_mu = E_reco_bin_edges_mu[0];
static const double E_reco_max_mu = E_reco_bin_edges_mu[E_reco_bins_mu];

// data and predicted backgrounds
static const double data_e[]      = { 497.0, 283.0, 313.0, 167.0, 184.0, 163.0,
                                      140.0, 115.0,  97.0,  98.0, 130.0 };
static const double data_e_bar[]  = { 122.0,  70.0,  65.0,  43.0,  57.0,  39.0,
                                       37.0,  23.0,  22.0,  30.0,  43.0 };
static const double data_mu[]     = { 37676.0, 59515.0, 53126.0, 37050.0,
                                      22150.0, 11478.0,  5374.0,  2547.0 };
static const double data_mu_bar[] = {  9481.0, 13581.0, 11308.0,  7667.0,
                                       4682.0,  2371.0,   985.0,   380.0 };

static const double pred_mu[] =
  { 38564.217639, 59339.405335, 53069.519495, 37171.337542,
    23002.153188, 12423.361945,  6012.845025,  2801.295291 };
static const double pred_mu_bar[] =
  {  9998.957451, 13461.301884, 11298.240453,  7604.960449,
     4331.886940,  2125.537108,   891.222608,   336.987112 };

static const double bg_e_unosc[] =              // total neutrino mode BG
  { 361.002334, 216.002142, 239.436776, 127.517957, 179.035344, 133.901816,
    139.020389, 113.446978,  81.204519,  98.603919, 137.953204 };
static const double bg_e_frac_nue_from_mu[] =   // fractional contrib. of \nu_e from \mu decay
  { 0.073, 0.152, 0.241, 0.303, 0.358, 0.463, 0.443, 0.454, 0.452, 0.466, 0.397 };
static const double bg_e_frac_nue_from_K[] =    // fractional contrib. of \nu_e from K decay
  { 0.034, 0.06, 0.095, 0.178, 0.199, 0.267, 0.306, 0.325, 0.391, 0.462, 0.663 };
static const double bg_e_frac_other[] =         // fractional contrib. of non-\nu_e BG
  { 0.894, 0.788, 0.664, 0.519, 0.442, 0.269, 0.251, 0.221, 0.157, 0.072, 0.0 };

static const double bg_e_bar_unosc[] =          // total antineutrino mode BG
  {  90.289907,  53.077595,  57.098801,  32.937945,  43.159072,  34.174322,
     36.383542,  28.737807,  22.339305,  26.509072,  42.697791 };

// covariance matrix
#define NCOV6 (n_cov_6)
#define NCOV4 (n_cov_4)
static const int n_cov_6 = 2 * (2*E_reco_bins_e + E_reco_bins_mu);
static const int n_cov_4 = 2 * (  E_reco_bins_e + E_reco_bins_mu);
#include MB_COV_MATRIX_FILE
static gsl_matrix *M6        = NULL; // 6x6 block convariance matrix
static gsl_matrix *M4        = NULL; // 4x4 block convariance matrix
static gsl_matrix *M4inv     = NULL; // inverse of M4
static gsl_permutation *perm = NULL;

// Global data structures
static double R_nu[E_true_bins][E_reco_bins_e];   // E_true -> E_reco mapping based on MC
static double R_nu_bar[E_true_bins][E_reco_bins_e];


/***************************************************************************
 * Initialize data structures required for MiniBooNE chi^2 calculation *
 ***************************************************************************/
int chiMB_jk_init()
{
  printf("# Flags in Joachim's MiniBooNE code: ");
  #ifdef OSC_BG
    printf("OSC_BG ");
  #endif
  #ifdef OSC_NORM
    printf("OSC_NORM ");
  #endif
  #ifdef OSC_NUMU
    printf("OSC_NUMU ");
  #endif
  printf("\n");

  FILE *f;

  // Sort MC events into detector response matrices
  memset(R_nu, 0, E_true_bins*E_reco_bins_e*sizeof(R_nu[0][0]));
  f = fopen(MB_DATA_DIR "miniboone_numunuefullosc_ntuple.txt", "r");
  if (!f)
  {
    fprintf(stderr, "chiMB_jk_init: MC n-tuple file for neutrino mode not found.\n");
    return -1;
  }
  int n_nu = 0;    // Event counter
  while (!feof(f))
  {
    int status, i_true, i_reco;
    double E_true, E_reco, L, w;
    status = fscanf(f, "%lg %lg %lg %lg", &E_true, &E_reco, &L, &w);
    if (status == EOF)
      break;
    i_true = E_true_bins * (E_true - E_true_min)   / (E_true_max   - E_true_min);
    i_reco = E_reco_bins_e;
    while (E_reco < E_reco_bin_edges_e[i_reco])
        i_reco--;
    if (i_true >= 0  &&  i_true < E_true_bins  &&
        i_reco >= 0  &&  i_reco < E_reco_bins_e)
      R_nu[i_true][i_reco] += w;
    n_nu++;
  }
  if (f) fclose(f);

  memset(R_nu_bar, 0, E_true_bins*E_reco_bins_e*sizeof(R_nu_bar[0][0]));
  f = fopen(MB_DATA_DIR "miniboone_nubarfullosc_ntuple.txt", "r");
  if (!f)
  {
    fprintf(stderr, "chiMB_jk_init: MC n-tuple file for antineutrino mode not found.\n");
    return -2;
  }
  int n_nu_bar = 0;
  while (!feof(f))
  {
    int status, i_true, i_reco;
    double E_true, E_reco, L, w;
    status = fscanf(f, "%lg %lg %lg %lg", &E_true, &E_reco, &L, &w);
    if (status == EOF)
      break;
    i_true = E_true_bins * (E_true - E_true_min)   / (E_true_max   - E_true_min);
    i_reco = E_reco_bins_e;
    while (E_reco < E_reco_bin_edges_e[i_reco])
        i_reco--;
    if (i_true >= 0  &&  i_true < E_true_bins  &&
        i_reco >= 0  &&  i_reco < E_reco_bins_e)
      R_nu_bar[i_true][i_reco] += w;
    n_nu_bar++;
  }
  if (f) fclose(f);

  // Normalize detector response matrices (see MB instructions)
  for (int it=0; it < E_true_bins; it++)
    for (int ir=0; ir < E_reco_bins_e; ir++)
    {
      R_nu[it][ir]     /= n_nu;
      R_nu_bar[it][ir] /= n_nu_bar;
    }

#ifdef NU_USE_NUSQUIDS
  // Compute initial spectrum for signal prediction by projecting R
  marray<double,3> inistate_nu{E_true_bins, 2, 4};
  marray<double,3> inistate_nu_bar{E_true_bins, 2, 4};
  std::fill(inistate_nu.begin(), inistate_nu.end(), 0.0);
  std::fill(inistate_nu_bar.begin(), inistate_nu_bar.end(), 0.0);
  for (int it=0; it < E_true_bins; it++)
    for (int ir=0; ir < E_reco_bins_e; ir++)
    {
      inistate_nu[it][0][1]     += R_nu[it][ir];
      inistate_nu_bar[it][1][1] += R_nu_bar[it][ir];
    }
  marray<double,1> E_true_range = regeneration::linspace(
      (E_true_min + 0.5 * (E_true_max - E_true_min) / E_true_bins) * units::MeV,
      (E_true_max - 0.5 * (E_true_max - E_true_min) / E_true_bins) * units::MeV,
      E_true_bins);                                     // Array of E_true bin centers
  prb_nu.setIniState(E_true_range, inistate_nu);
  prb_nu_bar.setIniState(E_true_range, inistate_nu_bar);

  // Initial spectrum for background prediction is just the MiniBooNE BG predictions
  marray<double,3> inistate_bg_nue_from_mu{E_reco_bins_e, 2, 4};
  marray<double,3> inistate_bg_nue_from_K{E_reco_bins_e, 2, 4};
  marray<double,1> E_reco_range{E_reco_bins_e};
  std::fill(inistate_bg_nue_from_mu.begin(), inistate_bg_nue_from_mu.end(), 0.0);
  std::fill(inistate_bg_nue_from_K.begin(), inistate_bg_nue_from_K.end(), 0.0);
  std::fill(E_reco_range.begin(), E_reco_range.end(), 0.0);
  for (int ir=0; ir < E_reco_bins_e; ir++)
  {
    E_reco_range[ir] = 0.5 * (E_reco_bin_edges_e[ir] + E_reco_bin_edges_e[ir+1]) * units::MeV;
    inistate_bg_nue_from_mu[ir][0][0] = bg_e_unosc[ir] * bg_e_frac_nue_from_mu[ir];
    inistate_bg_nue_from_K[ir][0][0]  = bg_e_unosc[ir] * bg_e_frac_nue_from_K[ir];
  }
  prb_bg_nue_from_mu.setIniState(E_reco_range, inistate_bg_nue_from_mu);
  prb_bg_nue_from_K.setIniState(E_reco_range, inistate_bg_nue_from_K);
#endif

//  for (int ir=0; ir < E_reco_bins_e; ir++)
//  {
//    double E_reco = 0.5 * (E_reco_bin_edges_e[ir] + E_reco_bin_edges_e[ir+1]);
//    marray<double,2> inistate_nue_from_mu;
//    marray<double,2> inistate_nue_from_K;
//    marray<double,2> finalstate_nue_from_mu;
//    marray<double,2> finalstate_nue_from_K;
//    inistate_nue_from_mu   = prb_bg_nue_from_mu.getInitialFlux(E_reco * units::MeV);
//    inistate_nue_from_K    = prb_bg_nue_from_K.getInitialFlux(E_reco * units::MeV);
//    finalstate_nue_from_mu = prb_bg_nue_from_mu.getFinalFlux(E_reco * units::MeV,
//                                                             MB_baseline * units::km);
//    finalstate_nue_from_K  = prb_bg_nue_from_K.getFinalFlux(E_reco * units::MeV,
//                                                            MB_baseline * units::km);
//    printf("%3d  %10.5g %10.5g   %10.5g %10.5g\n", ir,
//           finalstate_nue_from_K[0][0], finalstate_nue_from_K[1][0],
//           inistate_nue_from_K[0][0], inistate_nue_from_K[1][0]); //FIXME
//  }
//  getchar();

  // Initialize data structures for covariance matrix
  M6    = gsl_matrix_alloc(NCOV6, NCOV6);
  M4    = gsl_matrix_alloc(NCOV4, NCOV4);
  M4inv = gsl_matrix_alloc(NCOV4, NCOV4);
  perm  = gsl_permutation_alloc(NCOV4);

  return 0;
}


/***************************************************************************
 * Cleanup GSL data structures required for MiniBooNE chi^2 calculation    *
 ***************************************************************************/
int chiMB_jk_clear()
{
  if (perm)  { gsl_permutation_free(perm); perm  = NULL; }
  if (M4inv) { gsl_matrix_free(M4inv);     M4inv = NULL; }
  if (M4)    { gsl_matrix_free(M4);        M4    = NULL; }
  if (M6)    { gsl_matrix_free(M6);        M6    = NULL; }
  return 0;
}


/***************************************************************************
 * Calculate chi^2 for the MiniBooNE analysis                              *
 * see https://www-boone.fnal.gov/for_physicists/data_release/nuebar2010/  *
 * for the instructions which we loosely follow                            *
 ***************************************************************************
 * Parameters:                                                             *
 *   print_spectrum: whether (1) or not (0) to output the event spectra    *
 ***************************************************************************/
double chiMB_jk(int print_spectrum)
{
  // Evolve initial spectrum though Ivan's code
  double rates_e[E_reco_bins_e];        // Oscillated flux: \nu_e signal (nu mode)
  double rates_e_bar[E_reco_bins_e];    //                  \nu_e signal (anti-nu mode)
  double rates_mu[E_reco_bins_e];       //                  \nu_\mu signal (nu mode)
  double rates_mu_bar[E_reco_bins_e];   //                  \nu_\mu signal (anti-nu mode)
  double rates_bg_e[E_reco_bins_e];     //                  \nu_e BG (nu mode)
  double rates_bg_e_bar[E_reco_bins_e]; //                  \nu_e BG (anti-nu mode)

  // FIXME FIXME
//  {
//    extern Param p_oscdecay;
//    prb_nu.setParam(p_oscdecay);
//    prb_nu_bar.setParam(p_oscdecay);
//    prb_bg_nue_from_mu.setParam(p_oscdecay);
//    prb_bg_nue_from_K.setParam(p_oscdecay);
//    marray<double,2> inistate_tmp;
//    marray<double,2> finalstate_tmp;
//    double rates_nusquids[E_reco_bins_e];
//    double rates_globes[E_reco_bins_e];
//    memset(rates_nusquids, 0, E_reco_bins_e * sizeof(rates_nusquids[0]));
//    memset(rates_globes,   0, E_reco_bins_e * sizeof(rates_globes[0]));
//    for (int it=0; it < E_true_bins; it++)
//    {
//      double P[3][3];
//      double rho = 0.0;
//      double E_true  = E_true_min + (it+0.5) * (E_true_max - E_true_min) / E_true_bins;
//      snu_probability_matrix(P, +1, 0.001 * E_true, 1, &MB_baseline, &rho, 0., NULL);
//      inistate_tmp   = prb_nu.getInitialFlux(E_true * units::MeV);
//      finalstate_tmp = prb_nu.getFinalFlux(E_true * units::MeV, MB_baseline * units::km);
////      printf("ini %3d %7.5g   %10.5g %10.5g     %10.5g %10.5g %10.5g\n",
////              it, E_true, P[1][0] / P[1][1], finalstate_tmp[0][0] / finalstate_tmp[0][1],
////              inistate_tmp[0][1], finalstate_tmp[0][1], finalstate_tmp[0][0]);
//      for (int ir=0; ir < E_reco_bins_e; ir++)
//      {
//        if (finalstate_tmp[0][1] + finalstate_tmp[1][1] > 1e-10)         // neutrino mode
//          rates_nusquids[ir] += R_nu[it][ir] * (finalstate_tmp[0][0] + finalstate_tmp[1][0])
//                                        / (finalstate_tmp[0][1] + finalstate_tmp[1][1]);
//        rates_globes[ir] += R_nu[it][ir] * P[1][0] / P[1][1];
////        printf("   %3d %3d %10.7g    %10.7g %10.7g\n", it, ir, R_nu[it][ir],
////               R_nu[it][ir] * (finalstate_tmp[0][0] + finalstate_tmp[1][0])
////                            / (finalstate_tmp[0][1] + finalstate_tmp[1][1]),
////               R_nu[it][ir] * P[1][0] / P[1][1]);
////        printf("     ** %10.7g %10.7g\n", rates_nusquids[ir], rates_globes[ir]);
//      }
//    }
//    for (int ir=0; ir < E_reco_bins_e; ir++)
//    {
//      printf("final %3d %10.7g %10.7g\n", ir, rates_nusquids[ir], rates_globes[ir]);
//    }
//    getchar(); //FIXME
//  }

  memset(rates_e,     0, E_reco_bins_e * sizeof(rates_e[0]));
  memset(rates_e_bar, 0, E_reco_bins_e * sizeof(rates_e_bar[0]));
  for (int ir=0; ir < E_reco_bins_e; ir++)
  {
    rates_bg_e[ir] = bg_e_unosc[ir];
    rates_bg_e_bar[ir] = bg_e_bar_unosc[ir];
  }
  for (int ir=0; ir < E_reco_bins_mu; ir++)
  {
    rates_mu[ir]     = pred_mu[ir];
    rates_mu_bar[ir] = pred_mu_bar[ir];
  }

  // FIXME Sub-dominant wrong-sign contributions may oscillate differently

  // Neutrino oscillations from Ivan's code
  // --------------------------------------
#ifdef NU_USE_NUSQUIDS
  extern Param p_oscdecay;
  prb_nu.setParam(p_oscdecay);
  prb_nu_bar.setParam(p_oscdecay);
  prb_bg_nue_from_mu.setParam(p_oscdecay);
  prb_bg_nue_from_K.setParam(p_oscdecay);

  marray<double,2> inistate_nu;
  marray<double,2> inistate_nu_bar;
  marray<double,2> finalstate_nu;
  marray<double,2> finalstate_nu_bar;
  for (int it=0; it < E_true_bins; it++)
  {
    double E_true     = E_true_min + (it+0.5) * (E_true_max - E_true_min) / E_true_bins;
    inistate_nu       = prb_nu.getInitialFlux(E_true * units::MeV);
    inistate_nu_bar   = prb_nu_bar.getInitialFlux(E_true * units::MeV);
    finalstate_nu     = prb_nu.getFinalFlux(E_true * units::MeV, MB_baseline * units::km);
    finalstate_nu_bar = prb_nu_bar.getFinalFlux(E_true * units::MeV, MB_baseline * units::km);

    #ifdef OSC_NORM   // Including adjustment of normalization due to oscillations
      for (int ir=0; ir < E_reco_bins_e; ir++)
      {
        if (finalstate_nu[0][NU_MU] + finalstate_nu[1][NU_MU] > 1e-10)         // nu mode
          rates_e[ir] += R_nu[it][ir] * (finalstate_nu[0][NU_E]  + finalstate_nu[1][NU_E])
                                      / (finalstate_nu[0][NU_MU] + finalstate_nu[1][NU_MU]);
        if (finalstate_nu_bar[0][NU_MU] + finalstate_nu_bar[1][NU_MU] > 1e-10) // nubar mode
          rates_e_bar[ir] += R_nu_bar[it][ir]
                               * (finalstate_nu_bar[0][NU_E]  + finalstate_nu_bar[1][NU_E])
                               / (finalstate_nu_bar[0][NU_MU] + finalstate_nu_bar[1][NU_MU]);
                  /* Pme/Pmm - see discussion with W Louis:
                   * Flux is normalized to \nu_\mu rate, i.e. if there is \nu_\mu
                   * disappearance, the flux is underestimated by 1 / Pmm */
      } // for(ir)

    #else    // Onlye \nu_\mu -> \nu_e oscillations as in the official MiniBooNE analysis
      for (int ir=0; ir < E_reco_bins_e; ir++)
      {
        if (inistate_nu[0][1] > 1e-10)                    // neutrino mode
          rates_e[ir]     += R_nu[it][ir]
                               * (finalstate_nu[0][NU_E] + finalstate_nu[1][NU_E])
                               / inistate_nu[0][NU_MU];
        if (inistate_nu_bar[1][1] > 1e-10)                // antineutrino mode
          rates_e_bar[ir] += R_nu_bar[it][ir]
                               * (finalstate_nu_bar[0][NU_E] + finalstate_nu_bar[1][NU_E])
                               / inistate_nu_bar[1][NU_MU];
      } // for(ir)
    #endif // ifdef OSC_NORM
  } // for(it)

  // Background oscillations for \nu_e sample
  //   - we compute oscillation probabilities only at the centers of the E_reco bins
  //   - \nu_e BG from muon decays is normalized to \nu_\mu rate at MiniBooNE
  //   - \nu_e BG from kaon decays is normalized to \nu_\mu rate at SciBooNE
  #ifdef OSC_BG
    marray<double,2> inistate_numu;
    marray<double,2> inistate_nue_from_mu;
    marray<double,2> inistate_nue_from_K;
    marray<double,2> finalstate_numu_MB;
    marray<double,2> finalstate_numu_SB;
    marray<double,2> finalstate_nue_from_mu;
    marray<double,2> finalstate_nue_from_K;
    for (int ir=0; ir < E_reco_bins_e; ir++)
    {
      double E_reco = 0.5 * (E_reco_bin_edges_e[ir] + E_reco_bin_edges_e[ir+1]);
      inistate_numu          = prb_nu.getInitialFlux(E_reco * units::MeV);
      inistate_nue_from_mu   = prb_bg_nue_from_mu.getInitialFlux(E_reco*units::MeV);
      inistate_nue_from_K    = prb_bg_nue_from_K.getInitialFlux(E_reco*units::MeV);
      finalstate_numu_MB     = prb_nu.getFinalFlux(E_reco*units::MeV, MB_baseline*units::km);
      finalstate_numu_SB     = prb_nu.getFinalFlux(E_reco*units::MeV, SB_baseline*units::km);
      finalstate_nue_from_mu = prb_bg_nue_from_mu.getFinalFlux(E_reco*units::MeV,
                                                               MB_baseline*units::km);
      finalstate_nue_from_K  = prb_bg_nue_from_K.getFinalFlux(E_reco*units::MeV,
                                                              MB_baseline*units::km);

      rates_bg_e[ir] = bg_e_unosc[ir]
          * ( bg_e_frac_other[ir]
            + bg_e_frac_nue_from_mu[ir] 
                * (finalstate_nue_from_mu[0][NU_E] + finalstate_nue_from_mu[1][NU_E])
                / (inistate_nue_from_mu[0][NU_E] + inistate_nue_from_mu[1][NU_E])
                / (finalstate_numu_MB[0][NU_MU] + finalstate_numu_MB[1][NU_MU])
                * (inistate_numu[0][NU_MU] + inistate_numu[1][NU_MU])
            + bg_e_frac_nue_from_K[ir]
                * (finalstate_nue_from_K[0][NU_E] + finalstate_nue_from_K[1][NU_E])
                / (inistate_nue_from_K[0][NU_E] + inistate_nue_from_K[1][NU_E])
                / (finalstate_numu_SB[0][NU_MU] + finalstate_numu_SB[1][NU_MU])
                * (inistate_numu[0][NU_MU] + inistate_numu[1][NU_MU])
            );
      rates_bg_e_bar[ir] = bg_e_bar_unosc[ir];
    }
  #endif // ifdef(OSC_BG)

  #ifdef OSC_NUMU  // Oscillate muon neutrino rates
    for (int ir=0; ir < E_reco_bins_mu; ir++)
    {
      double E_reco = 0.5 * (E_reco_bin_edges_mu[ir] + E_reco_bin_edges_mu[ir+1]);
      inistate_nu       = prb_nu.getInitialFlux(E_reco*units::MeV);
      inistate_nu_bar   = prb_nu_bar.getInitialFlux(E_reco*units::MeV);
      finalstate_nu     = prb_nu.getFinalFlux(E_reco*units::MeV, MB_baseline*units::km);
      finalstate_nu_bar = prb_nu_bar.getFinalFlux(E_reco*units::MeV, MB_baseline*units::km);
      rates_mu[ir]     *= (finalstate_nu[0][NU_MU] + finalstate_nu[1][NU_MU])
                        / (inistate_nu[0][NU_MU] + inistate_nu[1][NU_MU]);
      rates_mu_bar[ir] *= (finalstate_nu_bar[0][NU_MU] + finalstate_nu_bar[1][NU_MU])
                        / (inistate_nu_bar[0][NU_MU] + inistate_nu_bar[1][NU_MU]);
    }
  #endif // ifdef OSC_NUMU

  // Oscillation probabilities from GLoBES
  // -------------------------------------
#else
  for (int it=0; it < E_true_bins; it++)
  {
    double E_true = E_true_min + (it+0.5) * (E_true_max - E_true_min) / E_true_bins;
    double P[3][3], Pbar[3][3];
    double rho = 0.0;
    snu_probability_matrix(P,    +1, 0.001 * E_true, 1, &MB_baseline, &rho, 0., NULL);
    snu_probability_matrix(Pbar, -1, 0.001 * E_true, 1, &MB_baseline, &rho, 0., NULL);
    for (int ir=0; ir < E_reco_bins_e; ir++)
    {
      #ifdef OSC_NORM
        rates_e[ir]     += R_nu[it][ir] * P[NU_MU][NU_E] / P[NU_MU][NU_MU];
        rates_e_bar[ir] += R_nu_bar[it][ir] * Pbar[NU_MU][NU_E] / Pbar[NU_MU][NU_MU];
                  /* Pme/Pmm - see discussion with W Louis:
                   * Flux is normalized to \nu_\mu rate, i.e. if there is \nu_\mu
                   * disappearance, the flux is underestimated by 1 / Pmm */
      #else
        rates_e[ir]     += R_nu[it][ir] * P[NU_MU][NU_E];
        rates_e_bar[ir] += R_nu_bar[it][ir] * Pbar[NU_MU][NU_E];
      #endif
    } // for (ir)
  } // for (it)

  #ifdef OSC_BG
    for (int ir=0; ir < E_reco_bins_e; ir++)
    {
      double P_MB[3][3], P_SB[3][3];
      double rho = 0.0;
      double E_reco = 0.5 * (E_reco_bin_edges_e[ir] + E_reco_bin_edges_e[ir+1]);
      snu_probability_matrix(P_MB, +1, 0.001 * E_reco, 1, &MB_baseline, &rho, 0., NULL);
      snu_probability_matrix(P_SB, +1, 0.001 * E_reco, 1, &SB_baseline, &rho, 0., NULL);
      rates_bg_e[ir] = bg_e_unosc[ir]
          * ( bg_e_frac_other[ir]
            + bg_e_frac_nue_from_mu[ir] * P_MB[NU_E][NU_E] / P_MB[NU_MU][NU_MU]
            + bg_e_frac_nue_from_K[ir]  * P_MB[NU_E][NU_E] / P_SB[NU_MU][NU_MU] );
      rates_bg_e_bar[ir] = bg_e_bar_unosc[ir];
    }
  #endif // ifdef(FULL_OSC)

  #ifdef OSC_NUMU  // Oscillate muon neutrino rates
    for (int ir=0; ir < E_reco_bins_mu; ir++)
    {
      double E_reco = 0.5 * (E_reco_bin_edges_mu[ir] + E_reco_bin_edges_mu[ir+1]);
      double P[3][3], Pbar[3][3];
      double rho = 0.0;
      snu_probability_matrix(P,    +1, 0.001 * E_reco, 1, &MB_baseline, &rho, 0., NULL);
      snu_probability_matrix(Pbar, -1, 0.001 * E_reco, 1, &MB_baseline, &rho, 0., NULL);
      rates_mu[ir]     *= P[NU_MU][NU_MU];
      rates_mu_bar[ir] *= P[NU_MU][NU_MU];
    }
  #endif // ifdef OSC_NUMU
#endif // ifdef(NU_USE_NUSQUIDS)


  /* Construct covariance matrix in 3x3 block form */
  double (*_M6)[NCOV6]    = (double (*)[NCOV6]) gsl_matrix_ptr(M6, 0, 0);
  double (*_M4)[NCOV4]    = (double (*)[NCOV4]) gsl_matrix_ptr(M4, 0, 0);
  double (*_M4inv)[NCOV4] = (double (*)[NCOV4]) gsl_matrix_ptr(M4inv, 0, 0);
  double P[NCOV6];        /* Vector of predicted nu_e signal, nu_e bg, nu_mu signal,
                             nu_e_bar signal, nu_e_bar bg, nu_mu_bar signal */
  int k = 0;
  for (int ir=0; ir < NE; ir++)
    P[k++] = rates_e[ir];
  for (int ir=0; ir < NE; ir++)
    P[k++] = rates_bg_e[ir];
  for (int ir=0; ir < NMU; ir++)
    P[k++] = rates_mu[ir];
  for (int ir=0; ir < NE; ir++)
    P[k++] = rates_e_bar[ir];
  for (int ir=0; ir < NE; ir++)
    P[k++] = rates_bg_e_bar[ir];
  for (int ir=0; ir < NMU; ir++)
    P[k++] = rates_mu_bar[ir];

  for (int i=0; i < NCOV6; i++)
    for (int j=0; j < NCOV6; j++)
      _M6[i][j] = Mfrac[i][j] * P[i] * P[j];
  for (int ir=0; ir < NE; ir++)
  {
    _M6[ir][ir]              += rates_e[ir];
    _M6[NCOV6/2+ir][NCOV6/2] += rates_e_bar[ir];
  }

  /* Collapse covariance matrix to 4x4 block form */
  gsl_matrix_set_zero(M4);
  for (int i=0; i < NCOV6; i++)        /* project horizontally */
  {
    for (int j=0; j < NE; j++)
      _M6[i][j]          = _M6[i][j] + _M6[i][j+NE];
    for (int j=0; j < NMU; j++)
      _M6[i][j+NE]       = _M6[i][j+2*NE];
    for (int j=0; j < NE; j++)
      _M6[i][j+NE+NMU]   = _M6[i][j+2*NE+NMU] + _M6[i][j+3*NE+NMU];
    for (int j=0; j < NMU; j++)
      _M6[i][j+2*NE+NMU] = _M6[i][j+4*NE+NMU];
  }
  for (int j=0; j < NCOV4; j++)        /* project vertically */
  {
    for (int i=0; i < NE; i++)
      _M4[i][j]          = _M6[i][j] + _M6[i+NE][j];
    for (int i=0; i < NMU; i++)
      _M4[i+NE][j]       = _M6[i+2*NE][j];
    for (int i=0; i < NE; i++)
      _M4[i+NE+NMU][j]   = _M6[i+2*NE+NMU][j] + _M6[i+3*NE+NMU][j];
    for (int i=0; i < NMU; i++)
      _M4[i+2*NE+NMU][j] = _M6[i+4*NE+NMU][j];
  }

  //FIXME
//  for (int i=0; i < NE + NMU; i++)
//    for (int j=0; j < NE + NMU; j++)
//      printf("cov %3d %3d  %10.5g\n", i, j, _M4[i][j]);
//  getchar();

  /* Invert covariance matrix and compute log-likelihood */
  int signum;
  gsl_linalg_LU_decomp(M4, perm, &signum);
  gsl_linalg_LU_invert(M4, perm, M4inv);

  double P4[NCOV4];
  for (int i=0; i < NE; i++)
  {
    P4[i]          = data_e[i] - (rates_e[i] + rates_bg_e[i]);
//    P4[i] = P4[i] * P4[i] / data_e[i]; //FIXME
//    printf("e     %3d  %10.5g  %10.5g  %10.5g\n", i, data_e[i], rates_e[i] + 0*rates_bg_e[i],
//        P4[i]); //FIXME 
  }
  for (int i=0; i < NMU; i++)
  {
    P4[i+NE]       = data_mu[i] - rates_mu[i];
//    P4[i+NE] = P4[i+NE] * P4[i+NE] / data_mu[i]; //FIXME
//    printf("mu    %3d  %10.5g  %10.5g  %10.5g\n", i, data_mu[i], rates_mu[i],
//        P4[i+NE]); //FIXME 
  }
  for (int i=0; i < NE; i++)
  {
    P4[i+NE+NMU]   = data_e_bar[i] - (rates_e_bar[i] + rates_bg_e_bar[i]);
//    P4[i+NE+NMU] = P4[i+NE+NMU] * P4[i+NE+NMU] / data_e_bar[i]; //FIXME
//    printf("ebar  %3d  %10.5g  %10.5g  %10.5g\n", i, data_e_bar[i],
//        rates_e_bar[i] + rates_bg_e_bar[i], P4[i+NE+NMU]); //FIXME 
  }
  for (int i=0; i < NMU; i++)
  {
    P4[i+2*NE+NMU] = data_mu_bar[i] - rates_mu_bar[i];
//    P4[i+2*NE+NMU] = P4[i+2*NE+NMU] * P4[i+2*NE+NMU] / data_mu_bar[i]; //FIXME
//    printf("mubar %3d  %10.5g  %10.5g  %10.5g\n", i, data_mu_bar[i], rates_mu_bar[i],
//        P4[i+2*NE+NMU]); //FIXME 
  }
//  getchar();  //FIXME

  double chi2 = 0.0;
//  for (int i=0; i < NCOV4; i++)
//    chi2 += P4[i]; 
//  printf("chi2 = %g\n", chi2); //FIXME

  for (int i=0; i < NCOV4; i++)
    for (int j=0; j < NCOV4; j++)
    {
      chi2 += P4[i] * _M4inv[i][j] * P4[j];
//      if (fabs(P4[i] * _M4inv[i][j] * P4[j]) > 1)//FIXME
//        printf("%3d %3d   %10.7g  %10.7g\n", i, j, P4[i] * _M4inv[i][j] * P4[j], chi2);//FIXME
    }
  chi2 += gsl_linalg_LU_lndet(M4);
//  printf(" chi2 = %g\n", chi2); //FIXME
//  getchar();


  // Output event spectrum if requested
  // (format is [signal, bg, 0, 0] for compatibility with Pedro's code)
  if (print_spectrum)
  {
    for (int ir=0; ir < NE; ir++)
    {
      printf("# MBJKSPECT   %10.7g %10.7g %10.7g %10.7g\n",
             rates_e[ir], rates_bg_e[ir], 0., 0.);
    }
  }

  return chi2;
}


