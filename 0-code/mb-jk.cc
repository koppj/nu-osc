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
#include "glb_error.h"
#include "osc_decay/osc_decay.h"
#include "nu.h"

// Flags that affect chi^2 calculation
//   OSC_NORM     : Take into account impact of \nu_\mu disapp. on normalization
//   OSC_BG       : Include oscillations of \nu_e backgrounds
//   OSC_NUMU     : Allow also the \nu_\mu sample to oscillate
//   OSC_NO_NUBAR : Omit anti-neutrino data from the fit
//   OSC_NO_NUMU  : Omit muon neutrinos and anti-neutrinos from the fit
//                  (I'm not sure this is consistent - use for debugging only)
// TAG_DEFS - This signals ./compile-and-run where to insert new #define's
#define OSC_NORM     // TAG_DEF - tells ./compile-and-run which #define's to replace
#define OSC_BG       // TAG_DEF
#define OSC_NUMU     // TAG_DEF
#define OSC_NO_NUBAR // TAG_DEF // FIXME

// Data structures from Ivan's oscillation+decay code
#ifdef NU_USE_NUSQUIDS
  using namespace regeneration;
  static regProb prb_nu(OSC_DECAY_MAJORANA);
  static regProb prb_nu_bar(OSC_DECAY_MAJORANA);
  static regProb prb_bg_nue_from_mu(OSC_DECAY_MAJORANA);
  static regProb prb_bg_nue_from_K(OSC_DECAY_MAJORANA);
  static regProb prb_bg_nue_bar_from_mu(OSC_DECAY_MAJORANA);
  static regProb prb_bg_nue_bar_from_K(OSC_DECAY_MAJORANA);
#endif

#define MB_DATA_DIR        "glb/mb-jk-2018/"
#define MB_COV_MATRIX_FILE "glb/mb-jk-2018/cov_matrix.h"

static const double MB_baseline = 0.520;   // [km] - MiniBooNE baseline
static const double SB_baseline = 0.100;   // [km] - SciBooNE baseline

// MiniBooNE low-pass filter
double mb_lowpass_width = 5.0;             // [MeV]

// Energ binning
static const int E_true_bins   =  680;       // binning in E_true
static const double E_true_min =  120.0;     // [MeV]
static const double E_true_max = 3520.0;     // [MeV]

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

static double pred_mu[] =
  { 38564.217639, 59339.405335, 53069.519495, 37171.337542,
    23002.153188, 12423.361945,  6012.845025,  2801.295291 };
static double pred_mu_bar[] =
  {  9998.957451, 13461.301884, 11298.240453,  7604.960449,
     4331.886940,  2125.537108,   891.222608,   336.987112 };

static double bg_e_unosc[] =                // total neutrino mode BG
  { 361.002334, 216.002142, 239.436776, 127.517957, 179.035344, 133.901816,
    139.020389, 113.446978,  81.204519,  98.603919, 137.953204 };
static double bg_e_frac_nue_from_mu[] =     // fractional contrib. of \nu_e from \mu decay
  { 0.073, 0.152, 0.241, 0.303, 0.358, 0.463, 0.443, 0.454, 0.452, 0.466, 0.397 };
static double bg_e_frac_nue_from_K[] =      // fractional contrib. of \nu_e from K decay
  { 0.034, 0.06, 0.095, 0.178, 0.199, 0.267, 0.306, 0.325, 0.391, 0.462, 0.603 };
static double bg_e_frac_other[] =           // fractional contrib. of non-\nu_e BG
  { 0.894, 0.788, 0.664, 0.519, 0.442, 0.269, 0.251, 0.221, 0.157, 0.072, 0.0 };

static const double bg_e_cc[] =             // official prediction for \nu_e CC BG
  {  39.839,   47.4743, 83.291, 63.1215, 103.11,
     98.5888, 105.197, 88.962,  63.9345,  86.81, 70.185 };
static const double bg_e_delta[] =          // official prediction for Delta BG
  {  40.531,   42.669,  54.947, 17.1765, 13.09,
      4.91125,  1.9665,  1.9665, 1.968,   1.6425, 0. };
static const double bg_e_pi0[] =            // official prediction for pi^0 BG
  { 213.771,   90.2422, 70.648, 27.4823, 40.0912,
     21.2825,  25.5495, 21.6315,12.789,  11.485,  0.};

static double bg_e_bar_unosc[] =            // total antineutrino mode BG
  {  90.289907,  53.077595,  57.098801,  32.937945,  43.159072,  34.174322,
     36.383542,  28.737807,  22.339305,  26.509072,  42.697791 };
static double bg_e_bar_frac_nue_from_mu[] = // fractional contrib. of \nu_e from \mu decay
  { 0.092, 0.148, 0.2, 0.253, 0.323, 0.35, 0.345, 0.361, 0.367, 0.319, 0.003 };
static double bg_e_bar_frac_nue_from_K[] =  // fractional contrib. of \nu_e from K decay
  { 0.071, 0.135, 0.209, 0.295, 0.343, 0.406, 0.442, 0.452, 0.475, 0.553, 0.81 };
static double bg_e_bar_frac_other[] =       // fractional contrib. of non-\nu_e BG
  { 0.837, 0.717, 0.59, 0.451, 0.334, 0.244, 0.213, 0.187, 0.158, 0.128, 0.187 };

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
int chiMB_jk_init(const char *bg_tune)
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
  #ifdef OSC_NO_NUBAR
    printf("OSC_NO_NUBAR ");
  #endif
  #ifdef OSC_NO_NUMU
    printf("OSC_NO_NUMU ");
  #endif
  if (bg_tune && strlen(bg_tune) > 0)
    printf("bg_tune=%s", bg_tune);
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
    volatile int status, i_true, i_reco;
      // Here, I ran into a GCC bug: with -O2, the code sometimes claims
      // that i_reco>=0 even though i_reco=-1 (indicating an energy outside
      // the simulated range). I assume this is because i_reco is relocated
      // into a processor register and then treated as unsigned int.
      // Including "volatile" is a workaround for this.
    double E_true, E_reco, L, w;
    status = fscanf(f, "%lg %lg %lg %lg", &E_true, &E_reco, &L, &w);
    if (status == EOF)
      break;
    i_true = E_true_bins * (E_true - E_true_min) / (E_true_max   - E_true_min);
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
    volatile int status, i_true, i_reco;
    double E_true, E_reco, L, w;
    status = fscanf(f, "%lg %lg %lg %lg", &E_true, &E_reco, &L, &w);
    if (status == EOF)
      break;
    i_true = E_true_bins * (E_true - E_true_min) / (E_true_max   - E_true_min);
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
  marray<double,3> inistate_bg_nue_bar_from_mu{E_reco_bins_e, 2, 4};
  marray<double,3> inistate_bg_nue_bar_from_K{E_reco_bins_e, 2, 4};
  marray<double,1> E_reco_range{E_reco_bins_e};
  std::fill(inistate_bg_nue_from_mu.begin(), inistate_bg_nue_from_mu.end(), 0.0);
  std::fill(inistate_bg_nue_from_K.begin(),  inistate_bg_nue_from_K.end(), 0.0);
  std::fill(inistate_bg_nue_bar_from_mu.begin(), inistate_bg_nue_bar_from_mu.end(), 0.0);
  std::fill(inistate_bg_nue_bar_from_K.begin(),  inistate_bg_nue_bar_from_K.end(), 0.0);
  std::fill(E_reco_range.begin(), E_reco_range.end(), 0.0);
  for (int ir=0; ir < E_reco_bins_e; ir++)
  {
    E_reco_range[ir] = 0.5 * (E_reco_bin_edges_e[ir] + E_reco_bin_edges_e[ir+1]) * units::MeV;
    inistate_bg_nue_from_mu[ir][0][0] = bg_e_unosc[ir] * bg_e_frac_nue_from_mu[ir];
    inistate_bg_nue_from_K[ir][0][0]  = bg_e_unosc[ir] * bg_e_frac_nue_from_K[ir];
    inistate_bg_nue_bar_from_mu[ir][0][0] = bg_e_bar_unosc[ir] * bg_e_bar_frac_nue_from_mu[ir];
    inistate_bg_nue_bar_from_K[ir][0][0]  = bg_e_bar_unosc[ir] * bg_e_bar_frac_nue_from_K[ir];
  }
  prb_bg_nue_from_mu.setIniState(E_reco_range, inistate_bg_nue_from_mu);
  prb_bg_nue_from_K.setIniState(E_reco_range,  inistate_bg_nue_from_K);
  prb_bg_nue_bar_from_mu.setIniState(E_reco_range, inistate_bg_nue_bar_from_mu);
  prb_bg_nue_bar_from_K.setIniState(E_reco_range,  inistate_bg_nue_bar_from_K);
#endif

  // Initialize data structures for covariance matrix
  M6    = gsl_matrix_alloc(NCOV6, NCOV6);
  M4    = gsl_matrix_alloc(NCOV4, NCOV4);
  M4inv = gsl_matrix_alloc(NCOV4, NCOV4);
  perm  = gsl_permutation_alloc(NCOV4);

  // load backgrounds from file
  if (bg_tune && strlen(bg_tune) > 0)
  {
    const int n_bg = 4;
    const char *bg[n_bg] = { "pi0", "singlephoton", "beam-nue", "numu-control" };
      // when adding entries here adjust also correction bg_e_frac_XXX!
    double *buffer[n_bg][2];
    for (int j=0; j < n_bg; j++)
    {
      char bg_file[512];
      int n_rows=-1;
      int status;
      sprintf(bg_file, "%s/%s-%s.dat", MB_DATA_DIR "/joachim", bg[j], bg_tune);
      if ((status=LoadNdAlloc(bg_file, buffer[j], 2, &n_rows)) != GLB_SUCCESS)
        return status;

      // check that number of rows and energy binning are correct
      if (strstr(bg[j], "numu"))
      {
        if (n_rows != E_reco_bins_mu)
        {
          fprintf(stderr, "chiMB_jk_init: invalid number of rows in file %s.\n", bg_file);
          return GLBERR_INVALID_FILE_FORMAT;
        }
        for (int i=0; i < E_reco_bins_mu; i++)
          if (fabs(buffer[j][0][i] - 0.5*(E_reco_bin_edges_mu[i] + E_reco_bin_edges_mu[i+1]))
                     / buffer[j][0][i] > 1e-5)
          {
            fprintf(stderr, "chiMB_jk_init: invalid energy (%g MeV) in row %d of file %s.\n",
                    buffer[j][0][i], i+1, bg_file);
            return GLBERR_INVALID_FILE_FORMAT;
          }

        // convert events/MeV to absolute event numbers
        for (int i=0; i < E_reco_bins_mu; i++)
          buffer[j][1][i] *= E_reco_bin_edges_mu[i+1] - E_reco_bin_edges_mu[i];
      } // end (numu)
      else
      {
        if (n_rows != E_reco_bins_e)
        {
          fprintf(stderr, "chiMB_jk_init: invalid number of rows in file %s.\n", bg_file);
          return GLBERR_INVALID_FILE_FORMAT;
        }
        for (int i=0; i < E_reco_bins_e; i++)
          if (fabs(buffer[j][0][i] - 0.5*(E_reco_bin_edges_e[i] + E_reco_bin_edges_e[i+1]))
                     / buffer[j][0][i] > 1e-5)
          {
            fprintf(stderr, "chiMB_jk_init: invalid energy (%g MeV) in row %d of file %s.\n",
                    buffer[j][0][i], i+1, bg_file);
            return GLBERR_INVALID_FILE_FORMAT;
          }

        // convert events/MeV to absolute event numbers
        for (int i=0; i < E_reco_bins_e; i++)
          buffer[j][1][i] *= E_reco_bin_edges_e[i+1] - E_reco_bin_edges_e[i];
      }
    } // end if (numu)

    // adjust BG prediction
    // FIXME Do the same also for antineutrino mode
    for (int i=0; i < E_reco_bins_e; i++)
    {
      double t = bg_e_unosc[i] - bg_e_pi0[i]   + buffer[0][1][i]
                               - bg_e_delta[i] + buffer[1][1][i]
                               - bg_e_cc[i]    + buffer[2][1][i];
      bg_e_frac_nue_from_mu[i] *= bg_e_unosc[i] / t;
      bg_e_frac_nue_from_K[i]  *= bg_e_unosc[i] / t;
      bg_e_frac_other[i]        = 1 - bg_e_frac_nue_from_mu[i] - bg_e_frac_nue_from_K[i];
      bg_e_unosc[i]             = t;
    }

    // adjust predictions for nu_mu control sample
    for (int i=0; i < E_reco_bins_mu; i++)
      pred_mu[i] = buffer[3][1][i];
  } // end if (bg_tune)

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
 *   print_spectrum: 0: no extra output                                    *
 *                   1: output signal and total BG rates                   *
 *                   2: output signal + separate BG rates + osc. params    *
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
  prb_bg_nue_bar_from_mu.setParam(p_oscdecay);
  prb_bg_nue_bar_from_K.setParam(p_oscdecay);

  marray<double,2> inistate_nu;
  marray<double,2> inistate_nu_bar;
  marray<double,2> finalstate_nu;
  marray<double,2> finalstate_nu_bar;
  for (int it=0; it < E_true_bins; it++)
  {
    double E_true     = E_true_min + (it+0.5) * (E_true_max - E_true_min) / E_true_bins;
    inistate_nu       = prb_nu.getInitialFlux(E_true * units::MeV);
    inistate_nu_bar   = prb_nu_bar.getInitialFlux(E_true * units::MeV);
    finalstate_nu     = prb_nu.getFinalFlux(E_true * units::MeV, MB_baseline * units::km,
                                            mb_lowpass_width * units::MeV);
    finalstate_nu_bar = prb_nu_bar.getFinalFlux(E_true * units::MeV, MB_baseline * units::km,
                                                mb_lowpass_width * units::MeV);

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
    marray<double,2> inistate_numu_bar;
    marray<double,2> inistate_nue_bar_from_mu;
    marray<double,2> inistate_nue_bar_from_K;
    marray<double,2> finalstate_numu_MB;
    marray<double,2> finalstate_numu_SB;
    marray<double,2> finalstate_nue_from_mu;
    marray<double,2> finalstate_nue_from_K;
    marray<double,2> finalstate_numu_bar_MB;
    marray<double,2> finalstate_numu_bar_SB;
    marray<double,2> finalstate_nue_bar_from_mu;
    marray<double,2> finalstate_nue_bar_from_K;
    for (int ir=0; ir < E_reco_bins_e; ir++)
    {
      double E_reco = 0.5 * (E_reco_bin_edges_e[ir] + E_reco_bin_edges_e[ir+1]);
      inistate_numu              = prb_nu.getInitialFlux(E_reco * units::MeV);
      inistate_nue_from_mu       = prb_bg_nue_from_mu.getInitialFlux(E_reco*units::MeV);
      inistate_nue_from_K        = prb_bg_nue_from_K.getInitialFlux(E_reco*units::MeV);
      inistate_numu_bar          = prb_nu_bar.getInitialFlux(E_reco * units::MeV);
      inistate_nue_bar_from_mu   = prb_bg_nue_bar_from_mu.getInitialFlux(E_reco*units::MeV);
      inistate_nue_bar_from_K    = prb_bg_nue_bar_from_K.getInitialFlux(E_reco*units::MeV);
      finalstate_numu_MB         = prb_nu.getFinalFlux(E_reco*units::MeV,
                                     MB_baseline*units::km, mb_lowpass_width*units::MeV);
      finalstate_numu_SB         = prb_nu.getFinalFlux(E_reco*units::MeV,
                                     SB_baseline*units::km, mb_lowpass_width*units::MeV);
      finalstate_nue_from_mu     = prb_bg_nue_from_mu.getFinalFlux(E_reco*units::MeV,
                                     MB_baseline*units::km, mb_lowpass_width*units::MeV);
      finalstate_nue_from_K      = prb_bg_nue_from_K.getFinalFlux(E_reco*units::MeV,
                                     MB_baseline*units::km, mb_lowpass_width*units::MeV);
      finalstate_numu_bar_MB     = prb_nu_bar.getFinalFlux(E_reco*units::MeV,
                                     MB_baseline*units::km, mb_lowpass_width*units::MeV);
      finalstate_numu_bar_SB     = prb_nu_bar.getFinalFlux(E_reco*units::MeV,
                                     SB_baseline*units::km, mb_lowpass_width*units::MeV);
      finalstate_nue_bar_from_mu = prb_bg_nue_bar_from_mu.getFinalFlux(E_reco*units::MeV,
                                     MB_baseline*units::km, mb_lowpass_width*units::MeV);
      finalstate_nue_bar_from_K  = prb_bg_nue_bar_from_K.getFinalFlux(E_reco*units::MeV,
                                     MB_baseline*units::km, mb_lowpass_width*units::MeV);

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
      rates_bg_e_bar[ir] = bg_e_bar_unosc[ir]
          * ( bg_e_bar_frac_other[ir]
            + bg_e_bar_frac_nue_from_mu[ir] 
                * (finalstate_nue_bar_from_mu[0][NU_E] + finalstate_nue_bar_from_mu[1][NU_E])
                / (inistate_nue_bar_from_mu[0][NU_E] + inistate_nue_bar_from_mu[1][NU_E])
                / (finalstate_numu_bar_MB[0][NU_MU] + finalstate_numu_bar_MB[1][NU_MU])
                * (inistate_numu_bar[0][NU_MU] + inistate_numu_bar[1][NU_MU])
            + bg_e_bar_frac_nue_from_K[ir]
                * (finalstate_nue_bar_from_K[0][NU_E] + finalstate_nue_bar_from_K[1][NU_E])
                / (inistate_nue_bar_from_K[0][NU_E] + inistate_nue_bar_from_K[1][NU_E])
                / (finalstate_numu_bar_SB[0][NU_MU] + finalstate_numu_bar_SB[1][NU_MU])
                * (inistate_numu_bar[0][NU_MU] + inistate_numu_bar[1][NU_MU])
            );
    }
  #endif // ifdef(OSC_BG)

  #ifdef OSC_NUMU  // Oscillate muon neutrino rates
    for (int ir=0; ir < E_reco_bins_mu; ir++)
    {
      double E_reco = 0.5 * (E_reco_bin_edges_mu[ir] + E_reco_bin_edges_mu[ir+1]);
      inistate_nu       = prb_nu.getInitialFlux(E_reco*units::MeV);
      inistate_nu_bar   = prb_nu_bar.getInitialFlux(E_reco*units::MeV);
      finalstate_nu     = prb_nu.getFinalFlux(E_reco*units::MeV, MB_baseline*units::km,
                                              mb_lowpass_width*units::MeV);
      finalstate_nu_bar = prb_nu_bar.getFinalFlux(E_reco*units::MeV, MB_baseline*units::km,
                                                  mb_lowpass_width*units::MeV);
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
    double Pe[3][3], Pbare[3][3];
    double rho = 0.0;
    double filter_sigma = 0.0; //0.1 * (E_true_max - E_true_min) / E_true_bins;
    snu_probability_matrix(Pe,   +1,0.001*E_true,1,&MB_baseline,&rho,filter_sigma,NULL);
    snu_probability_matrix(Pbare,-1,0.001*E_true,1,&MB_baseline,&rho,filter_sigma,NULL);
    for (int ir=0; ir < E_reco_bins_e; ir++)
    {
      #ifdef OSC_NORM
        rates_e[ir]     += R_nu[it][ir] * Pe[NU_MU][NU_E] / Pe[NU_MU][NU_MU];
        rates_e_bar[ir] += R_nu_bar[it][ir] * Pbare[NU_MU][NU_E] / Pbare[NU_MU][NU_MU];
                  /* Pme/Pmm - see discussion with W Louis:
                   * Flux is normalized to \nu_\mu rate, i.e. if there is \nu_\mu
                   * disappearance, the flux is underestimated by 1 / Pmm */
      #else
        rates_e[ir]     += R_nu[it][ir] * Pe[NU_MU][NU_E];
        rates_e_bar[ir] += R_nu_bar[it][ir] * Pbare[NU_MU][NU_E];
      #endif
    } // for (ir)
  } // for (it)

  #ifdef OSC_BG
    for (int ir=0; ir < E_reco_bins_e; ir++)
    {
//      double P_MB[3][3], P_SB[3][3];
//      double P_bar_MB[3][3], P_bar_SB[3][3];
//      double rho = 0.0;
//      double E_reco = 0.5 * (E_reco_bin_edges_e[ir] + E_reco_bin_edges_e[ir+1]);
//      double filter_sigma = 0.0;//0.8 * (E_reco_max_e - E_reco_min_e) / E_reco_bins_e;
//      snu_probability_matrix(P_MB,+1,0.001*E_reco,1,&MB_baseline,&rho,filter_sigma,NULL);
//      snu_probability_matrix(P_SB,+1,0.001*E_reco,1,&SB_baseline,&rho,filter_sigma,NULL);
//      snu_probability_matrix(P_bar_MB,-1,0.001*E_reco,1,&MB_baseline,&rho,filter_sigma,NULL);
//      snu_probability_matrix(P_bar_SB,-1,0.001*E_reco,1,&SB_baseline,&rho,filter_sigma,NULL);
//      rates_bg_e[ir] = bg_e_unosc[ir]
//          * ( bg_e_frac_other[ir]
//            + bg_e_frac_nue_from_mu[ir] * P_MB[NU_E][NU_E] / P_MB[NU_MU][NU_MU]
//            + bg_e_frac_nue_from_K[ir]  * P_MB[NU_E][NU_E] / P_SB[NU_MU][NU_MU] );
//      rates_bg_e_bar[ir] = bg_e_bar_unosc[ir]
//          * ( bg_e_bar_frac_other[ir]
//            + bg_e_bar_frac_nue_from_mu[ir]*P_bar_MB[NU_E][NU_E]/P_bar_MB[NU_MU][NU_MU]
//            + bg_e_bar_frac_nue_from_K[ir] *P_bar_MB[NU_E][NU_E]/P_bar_SB[NU_MU][NU_MU] );

      const int n = 20;
      double dE = (E_reco_bin_edges_e[ir+1] - E_reco_bin_edges_e[ir]) / n;
      double P_MB_ee=0., P_MB_mumu=0., P_bar_MB_ee=0., P_bar_MB_mumu=0.;
      double P_SB_mumu=0., P_bar_SB_mumu=0.;
      double P_MB[3][3], P_SB[3][3];
      double P_bar_MB[3][3], P_bar_SB[3][3];
      double rho = 0.0;
      double filter_sigma = 0.2 * dE;
      for (int j=0; j < n; j++)
      {
        double E = E_reco_bin_edges_e[ir] + (j+0.5) * dE;
        snu_probability_matrix(P_MB,+1,0.001*E,1,&MB_baseline,&rho,filter_sigma,NULL);
        snu_probability_matrix(P_SB,+1,0.001*E,1,&SB_baseline,&rho,filter_sigma,NULL);
        snu_probability_matrix(P_bar_MB,-1,0.001*E,1,&MB_baseline,&rho,filter_sigma,NULL);
        snu_probability_matrix(P_bar_SB,-1,0.001*E,1,&SB_baseline,&rho,filter_sigma,NULL);
        P_MB_ee       += P_MB[NU_E][NU_E]       / n;
        P_MB_mumu     += P_MB[NU_MU][NU_MU]     / n;
        P_bar_MB_ee   += P_bar_MB[NU_E][NU_E]   / n;
        P_bar_MB_mumu += P_bar_MB[NU_MU][NU_MU] / n;
        P_SB_mumu     += P_SB[NU_MU][NU_MU]     / n;
        P_bar_SB_mumu += P_bar_SB[NU_MU][NU_MU] / n;
      }
      rates_bg_e[ir] = bg_e_unosc[ir]
          * ( bg_e_frac_other[ir]
            + bg_e_frac_nue_from_mu[ir] * P_MB_ee / P_MB_mumu
            + bg_e_frac_nue_from_K[ir]  * P_MB_ee / P_SB_mumu );
      rates_bg_e_bar[ir] = bg_e_bar_unosc[ir]
          * ( bg_e_bar_frac_other[ir]
            + bg_e_bar_frac_nue_from_mu[ir]*P_bar_MB_ee/P_bar_MB_mumu
            + bg_e_bar_frac_nue_from_K[ir] *P_bar_MB_ee/P_bar_SB_mumu );
    }
  #endif // ifdef(OSC_BG)

  #ifdef OSC_NUMU  // Oscillate muon neutrino rates
    for (int ir=0; ir < E_reco_bins_mu; ir++)
    {
//      double E_reco = 0.5 * (E_reco_bin_edges_mu[ir] + E_reco_bin_edges_mu[ir+1]);
//      double Pmu[3][3], Pbarmu[3][3];
//      double rho = 0.0;
//      double filter_sigma = 0.0;//0.8 * (E_reco_max_mu - E_reco_min_mu) / E_reco_bins_mu;
//      snu_probability_matrix(Pmu,   +1,0.001*E_reco,1,&MB_baseline,&rho,filter_sigma,NULL);
//      snu_probability_matrix(Pbarmu,-1,0.001*E_reco,1,&MB_baseline,&rho,filter_sigma,NULL);
//      rates_mu[ir]     *= Pmu[NU_MU][NU_MU];
//      rates_mu_bar[ir] *= Pbarmu[NU_MU][NU_MU];

      const int n = 20;
      double P=0., Pbar=0.;
      double dE = (E_reco_bin_edges_mu[ir+1] - E_reco_bin_edges_mu[ir]) / n;
      double Pmu[3][3], Pbarmu[3][3];
      double rho = 0.0;
      double filter_sigma = 0.2 * dE;// FIXME FIXME MeV -> GeV conversion!!
      for (int j=0; j < n; j++)
      {
        double E = E_reco_bin_edges_mu[ir] + (j+0.5) * dE;
        snu_probability_matrix(Pmu,   +1,0.001*E,1,&MB_baseline,&rho,filter_sigma,NULL);
        snu_probability_matrix(Pbarmu,-1,0.001*E,1,&MB_baseline,&rho,filter_sigma,NULL);
        P    += Pmu[NU_MU][NU_MU]    / n;
        Pbar += Pbarmu[NU_MU][NU_MU] / n;
      }
      rates_mu[ir]     *= P;
      rates_mu_bar[ir] *= Pbar;
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

  #ifdef OSC_NO_NUBAR
    for (int i=0; i < NE+NMU; i++)
      for (int j=NE+NMU; j < NCOV4; j++)
        _M4[i][j] = 0.;
    for (int i=NE+NMU; i < NCOV4; i++)
    {
      for (int j=0; j < NCOV4; j++)
        _M4[i][j] = 0.;
      _M4[i][i] = 1.;
    }
  #endif

  #ifdef OSC_NO_NUMU
    for (int i=0; i < NCOV4; i++)
    {
      for (int j=NE; j < NE+NMU; j++)
        _M4[i][j] = _M4[j][i] = 0.;
      for (int j=2*NE+NMU; j < NCOV4; j++)
        _M4[i][j] = _M4[j][i] = 0.;
    }
    for (int i=NE; i < NE+NMU; i++)
      _M4[i][i] = 1.;
    for (int i=2*NE+NMU; i < NCOV4; i++)
      _M4[i][i] = 1.;
  #endif

  /* Invert covariance matrix and compute log-likelihood */
  int signum;
  gsl_linalg_LU_decomp(M4, perm, &signum);
  gsl_linalg_LU_invert(M4, perm, M4inv);

  double P4[NCOV4];
  for (int i=0; i < NE; i++)
    P4[i]          = data_e[i] - (rates_e[i] + rates_bg_e[i]);
  for (int i=0; i < NMU; i++)
    P4[i+NE]       = data_mu[i] - rates_mu[i];
  for (int i=0; i < NE; i++)
    P4[i+NE+NMU]   = data_e_bar[i] - (rates_e_bar[i] + rates_bg_e_bar[i]);
  for (int i=0; i < NMU; i++)
    P4[i+2*NE+NMU] = data_mu_bar[i] - rates_mu_bar[i];

  #ifdef OSC_NO_NUBAR
    for (int i=NE+NMU; i < NCOV4; i++)
      P4[i] = 0.;
  #endif

  #ifdef OSC_NO_NUMU
    for (int i=NE; i < NE+NMU; i++)
      P4[i] = 0.;
    for (int i=2*NE+NMU; i < NCOV4; i++)
      P4[i] = 0.;
  #endif

  double chi2 = 0.0;
  for (int i=0; i < NCOV4; i++)
    for (int j=0; j < NCOV4; j++)
    {
      chi2 += P4[i] * _M4inv[i][j] * P4[j];
    }
  chi2 += gsl_linalg_LU_lndet(M4);

  // Debug code for printing the covariance matrix
//  printf("# P4 ");
//  for (int i=0; i < NCOV4; i++)
//    printf("%10.5g ", P4[i]);
//  printf("\n");
//  for (int i=0; i < NCOV4; i++)
//  {
//    printf("# McovInv ");
//    for (int j=0; j < NCOV4; j++)
//      printf("%10.5g ", _M4inv[i][j]);
//    printf("\n");
//  }

  // Output event spectrum if requested
  // (format is [signal, bg, 0, 0] for compatibility with Pedro's code)
  if (print_spectrum)
  {
    for (int ir=0; ir < NE; ir++)
    {
      printf("# MBJKSPECT   %10.7g %10.7g %10.7g %10.7g\n",
             rates_e[ir], rates_bg_e[ir], rates_e_bar[ir], rates_bg_e_bar[ir]);
    }

    // print also individual BG components
    if (print_spectrum >= 2)
    {
      for (int ir=0; ir < NMU; ir++)
        printf("# MBJKSPECT MU      %10.7g %10.7g\n", rates_mu[ir], rates_mu_bar[ir]);

    #ifndef NU_USE_NUSQUIDS
      for (int ir=0; ir < NE; ir++)
      {
        double P_MB[3][3], P_SB[3][3];
        double rho = 0.0;
        double E_reco = 0.5 * (E_reco_bin_edges_e[ir] + E_reco_bin_edges_e[ir+1]);
        snu_probability_matrix(P_MB,+1,0.001 * E_reco,1,&MB_baseline,&rho,0.,NULL);
        snu_probability_matrix(P_SB,+1,0.001 * E_reco,1,&SB_baseline,&rho,0.,NULL);
        #ifdef OSC_BG
          printf("# MBJKSPECT BG      %10.7g %10.7g %10.7g\n",
            bg_e_unosc[ir]*bg_e_frac_other[ir],
            bg_e_unosc[ir]*bg_e_frac_nue_from_mu[ir]*P_MB[NU_E][NU_E]/P_MB[NU_MU][NU_MU],
            bg_e_unosc[ir]*bg_e_frac_nue_from_K[ir] *P_MB[NU_E][NU_E]/P_SB[NU_MU][NU_MU]);
        #else
          printf("# MBJKSPECT BG      %10.7g %10.7g %10.7g\n",
            bg_e_unosc[ir]*bg_e_frac_other[ir],
            bg_e_unosc[ir]*bg_e_frac_nue_from_mu[ir],
            bg_e_unosc[ir]*bg_e_frac_nue_from_K[ir]);
        #endif
      }

      glb_params p = glbAllocParams();
      snu_get_oscillation_parameters(p, NULL);
      printf("# MBJKSPECT PARAMS ");
      my_print_params(p);
      glbFreeParams(p);
    }
    #endif

    printf("# CHI2   %10.7g\n", chi2);
  }

  return chi2;
}


