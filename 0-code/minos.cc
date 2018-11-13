/***************************************************************************
 * Functions for MINOS CC and NC analysis                                  *
 ***************************************************************************
 * Author: J. Kopp                                                         *
 ***************************************************************************/
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <globes/globes.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_erf.h>
#include <root/TROOT.h>
#include <root/TFile.h>
#include <root/TMatrixD.h>
#include <root/TH1.h>
#include <root/TH2.h>
#include "glb_smear.h"
#include "nu.h"
using namespace std;

#define MINOS_2017_USE_GLOBES_LOWPASS // Use GLoBES low-pass filter (or MINOS averaging)
//#define MINOS_2017_FREE_NORM          // Leave overall normalization free? FIXME
#define MINOS_2017_MINIMIZE      // Minimize over nuisance parametersin 2017 analysis?
#ifdef MINOS_2017_FREE_NORM
  #define MINOS_2017_N_NUISANCE  5
#else
  #define MINOS_2017_N_NUISANCE  4 // # of nuisance parameters in MINOS/MINOS+ 2017 analysi
#endif

// ROOT file containing MINOS data
static TFile *f = NULL;

// Data structures for MINOS-2016 analysis
// ---------------------------------------

// Near detector MC histograms (response functions)
static TH2D *NDNC_TrueNC=NULL, *NDNC_NuMu=NULL, *NDNC_BeamNue=NULL,
            *NDNC_AppNue=NULL, *NDNC_AppNuTau=NULL;
static TH2D *NDCC_TrueNC=NULL, *NDCC_NuMu=NULL, *NDCC_BeamNue=NULL,
            *NDCC_AppNue=NULL, *NDCC_AppNuTau;
static double *sm_nd_nc_truenc=NULL, *sm_nd_nc_numu=NULL, *sm_nd_nc_beamnue=NULL,
              *sm_nd_nc_appnue=NULL, *sm_nd_nc_appnutau=NULL;;
static double *sm_nd_cc_truenc=NULL, *sm_nd_cc_numu=NULL, *sm_nd_cc_beamnue=NULL,
              *sm_nd_cc_appnue=NULL, *sm_nd_cc_appnutau=NULL;;

// Far detector MC histograms (responcse functions)
static TH2D *FDNC_TrueNC=NULL, *FDNC_NuMu=NULL, *FDNC_BeamNue=NULL,
            *FDNC_AppNue=NULL, *FDNC_AppNuTau;
static TH2D *FDCC_TrueNC=NULL, *FDCC_NuMu=NULL, *FDCC_BeamNue=NULL,
            *FDCC_AppNue=NULL, *FDCC_AppNuTau;
static double *sm_fd_nc_truenc=NULL, *sm_fd_nc_numu=NULL, *sm_fd_nc_beamnue=NULL,
              *sm_fd_nc_appnue=NULL, *sm_fd_nc_appnutau=NULL;;
static double *sm_fd_cc_truenc=NULL, *sm_fd_cc_numu=NULL, *sm_fd_cc_beamnue=NULL,
              *sm_fd_cc_appnue=NULL, *sm_fd_cc_appnutau=NULL;;

// Data
static TH1D *FD_dataNC=NULL, *FD_dataCC=NULL, *ND_dataNC=NULL, *ND_dataCC;
static TH1D *RdataNC=NULL, *RdataCC=NULL;
static double *data_r_nc=NULL, *data_r_cc=NULL;

// Covariance matrices
TMatrixD *CoVarCCinvert;
TMatrixD *CoVarNCinvert;
Double_t *inv_cov_matrix_nc=NULL, *inv_cov_matrix_cc=NULL;
double sigma_nd_norm_nc, sigma_nd_norm_cc;  // ND normalization uncertainty
double data_total_nd_nc, data_total_nd_cc;


// Data structures for 2017 MINOS/MINOS+ analysis
// ----------------------------------------------

// Near detector MC histograms (response functions)
static TH2D *NDNC_TrueNC_minos=NULL, *NDNC_NuMu_minos=NULL, *NDNC_BeamNue_minos=NULL,
            *NDNC_AppNue_minos=NULL, *NDNC_AppNuTau_minos=NULL;
static TH2D *NDCC_TrueNC_minos=NULL, *NDCC_NuMu_minos=NULL, *NDCC_BeamNue_minos=NULL,
            *NDCC_AppNue_minos=NULL, *NDCC_AppNuTau_minos;
static TH2D *NDNC_TrueNC_mplus=NULL, *NDNC_NuMu_mplus=NULL, *NDNC_BeamNue_mplus=NULL,
            *NDNC_AppNue_mplus=NULL, *NDNC_AppNuTau_mplus=NULL;
static TH2D *NDCC_TrueNC_mplus=NULL, *NDCC_NuMu_mplus=NULL, *NDCC_BeamNue_mplus=NULL,
            *NDCC_AppNue_mplus=NULL, *NDCC_AppNuTau_mplus;
static double *sm_nd_nc_truenc_minos=NULL, *sm_nd_nc_numu_minos=NULL, *sm_nd_nc_beamnue_minos=NULL,
              *sm_nd_nc_appnue_minos=NULL, *sm_nd_nc_appnutau_minos=NULL;;
static double *sm_nd_cc_truenc_minos=NULL, *sm_nd_cc_numu_minos=NULL, *sm_nd_cc_beamnue_minos=NULL,
              *sm_nd_cc_appnue_minos=NULL, *sm_nd_cc_appnutau_minos=NULL;;
static double *sm_nd_nc_truenc_mplus=NULL, *sm_nd_nc_numu_mplus=NULL, *sm_nd_nc_beamnue_mplus=NULL,
              *sm_nd_nc_appnue_mplus=NULL, *sm_nd_nc_appnutau_mplus=NULL;;
static double *sm_nd_cc_truenc_mplus=NULL, *sm_nd_cc_numu_mplus=NULL, *sm_nd_cc_beamnue_mplus=NULL,
              *sm_nd_cc_appnue_mplus=NULL, *sm_nd_cc_appnutau_mplus=NULL;;

// Far detector MC histograms (responcse functions)
static TH2D *FDNC_TrueNC_minos=NULL, *FDNC_NuMu_minos=NULL, *FDNC_BeamNue_minos=NULL,
            *FDNC_AppNue_minos=NULL, *FDNC_AppNuTau_minos;
static TH2D *FDCC_TrueNC_minos=NULL, *FDCC_NuMu_minos=NULL, *FDCC_BeamNue_minos=NULL,
            *FDCC_AppNue_minos=NULL, *FDCC_AppNuTau_minos;
static TH2D *FDNC_TrueNC_mplus=NULL, *FDNC_NuMu_mplus=NULL, *FDNC_BeamNue_mplus=NULL,
            *FDNC_AppNue_mplus=NULL, *FDNC_AppNuTau_mplus;
static TH2D *FDCC_TrueNC_mplus=NULL, *FDCC_NuMu_mplus=NULL, *FDCC_BeamNue_mplus=NULL,
            *FDCC_AppNue_mplus=NULL, *FDCC_AppNuTau_mplus;
static double *sm_fd_nc_truenc_minos=NULL, *sm_fd_nc_numu_minos=NULL, *sm_fd_nc_beamnue_minos=NULL,
              *sm_fd_nc_appnue_minos=NULL, *sm_fd_nc_appnutau_minos=NULL;;
static double *sm_fd_cc_truenc_minos=NULL, *sm_fd_cc_numu_minos=NULL, *sm_fd_cc_beamnue_minos=NULL,
              *sm_fd_cc_appnue_minos=NULL, *sm_fd_cc_appnutau_minos=NULL;;
static double *sm_fd_nc_truenc_mplus=NULL, *sm_fd_nc_numu_mplus=NULL, *sm_fd_nc_beamnue_mplus=NULL,
              *sm_fd_nc_appnue_mplus=NULL, *sm_fd_nc_appnutau_mplus=NULL;;
static double *sm_fd_cc_truenc_mplus=NULL, *sm_fd_cc_numu_mplus=NULL, *sm_fd_cc_beamnue_mplus=NULL,
              *sm_fd_cc_appnue_mplus=NULL, *sm_fd_cc_appnutau_mplus=NULL;;

// Data
static TH1D *FD_dataNC_minos=NULL, *FD_dataCC_minos=NULL;
static TH1D *ND_dataNC_minos=NULL, *ND_dataCC_minos;
static TH1D *FD_dataNC_mplus=NULL, *FD_dataCC_mplus=NULL;
static TH1D *ND_dataNC_mplus=NULL, *ND_dataCC_mplus;
static double *data_nd_nc=NULL, *data_nd_cc=NULL;
static double *data_fd_nc=NULL, *data_fd_cc=NULL;

// Predicted rates
static double *pred_nd_nc_minos=NULL, *pred_nd_cc_minos=NULL;
static double *pred_fd_nc_minos=NULL, *pred_fd_cc_minos=NULL;
static double *pred_nd_nc_mplus=NULL, *pred_nd_cc_mplus=NULL;
static double *pred_fd_nc_mplus=NULL, *pred_fd_cc_mplus=NULL;
static gsl_vector *pred_nc=NULL, *pred_cc=NULL;
static gsl_vector *temp_nc=NULL, *temp_cc=NULL; // Temporary space needed during computation

// Weights for beam focusing nuisance parameters
static TH1D *CC_IMiscal_ND_minos=NULL, *CC_IDist_ND_minos=NULL;
static TH1D *CC_IMiscal_FD_minos=NULL, *CC_IDist_FD_minos=NULL;
static TH1D *NC_IMiscal_ND_minos=NULL, *NC_IDist_ND_minos=NULL;
static TH1D *NC_IMiscal_FD_minos=NULL, *NC_IDist_FD_minos=NULL;
static TH1D *CC_IMiscal_ND_mplus=NULL, *CC_IDist_ND_mplus=NULL;
static TH1D *CC_IMiscal_FD_mplus=NULL, *CC_IDist_FD_mplus=NULL;
static TH1D *NC_IMiscal_ND_mplus=NULL, *NC_IDist_ND_mplus=NULL;
static TH1D *NC_IMiscal_FD_mplus=NULL, *NC_IDist_FD_mplus=NULL;
static double *wgt_miscal_nd_cc_minos=NULL, *wgt_dist_nd_cc_minos=NULL;
static double *wgt_miscal_fd_cc_minos=NULL, *wgt_dist_fd_cc_minos=NULL;
static double *wgt_miscal_nd_nc_minos=NULL, *wgt_dist_nd_nc_minos=NULL;
static double *wgt_miscal_fd_nc_minos=NULL, *wgt_dist_fd_nc_minos=NULL;
static double *wgt_miscal_nd_cc_mplus=NULL, *wgt_dist_nd_cc_mplus=NULL;
static double *wgt_miscal_fd_cc_mplus=NULL, *wgt_dist_fd_cc_mplus=NULL;
static double *wgt_miscal_nd_nc_mplus=NULL, *wgt_dist_nd_nc_mplus=NULL;
static double *wgt_miscal_fd_nc_mplus=NULL, *wgt_dist_fd_nc_mplus=NULL;

// Covariance matrices
TMatrixD *CoVarCC_relative;
TMatrixD *CoVarNC_relative;
double *cov_matrix_nc_rel=NULL, *cov_matrix_cc_rel=NULL; // Cov. matrices from files
gsl_matrix *cov_matrix_nc_abs=NULL, *cov_matrix_cc_abs=NULL; // Rescaled cov. matrices
gsl_permutation *cov_matrix_perm_nc=NULL, *cov_matrix_perm_cc=NULL; // Permutation for LU decomp.
//double sigma_nd_norm_nc, sigma_nd_norm_cc;  // ND normalization uncertainty
//double data_total_nd_nc, data_total_nd_cc;


/* Square of real number */
inline static double square(double x)
{
  return x*x;
}

/* Poisson likelihood */
inline static double poisson_likelihood(double true_rate, double fit_rate)
{
  double res;
  res = fit_rate - true_rate;
  if (true_rate > 0)
  {
    if (fit_rate <= 0.0)
      res = 1e100;
    else
      res += true_rate * log(true_rate/fit_rate);
  }
  else
    res = fabs(res);

  return 2.0 * res;
}


#define TOL 1.0E-7

/***************************************************************************
 * Function minos_smear                                                    *
 ***************************************************************************
 * Modify the smearing matrix in such a way that the effect of the         *
 * finite length of the decay pipe is approximately accounted for.         *
 * This function is based on glb_filter_compensate.                        *
 * See also Mathematica worksheet devel.nb, in which the expression for    *
 * the additional smearing width sigma is developed.                       *
 ***************************************************************************
 * Parameters:                                                             *
 *   smear:  glb_smear structure containing metadata about the smearing    *
 *   matrix: The smearing matrix entries                                   *
 *   lower:  Array of indices of lowest nonzero entry in each row          *
 *   upper:  Array of indices of highest nonzero entry in each row         *
 ***************************************************************************/
int minos_smear(glb_smear *s, double **matrix, int *lower, int *upper)
{
  fprintf(stderr, "minos_smear: This function has not been debugged yet");
  exit(-1);

  if (!s || !lower || !upper || !matrix)
  {
    fprintf(stderr, "minos_smear: NULL pointer argument encountered.");
    return -1;
  }

  int factor = 16;
  int hires_bins = factor * s->simbins;

  gsl_matrix *Sin  = gsl_matrix_alloc(s->numofbins, s->simbins);
  gsl_matrix *Sout = gsl_matrix_alloc(s->numofbins, s->simbins);
  gsl_matrix *R    = gsl_matrix_alloc(s->simbins, hires_bins);
  gsl_matrix *SR   = gsl_matrix_alloc(s->numofbins, hires_bins);
  gsl_matrix *F    = gsl_matrix_alloc(s->simbins, s->simbins);
  gsl_matrix *V    = gsl_matrix_alloc(s->simbins, s->simbins);
  gsl_vector *b    = gsl_vector_alloc(hires_bins);
  gsl_vector *sv   = gsl_vector_alloc(s->simbins);
  gsl_vector *x    = gsl_vector_alloc(s->simbins);

  if (!Sin || !Sout | !R || !SR || !F || !V || !b || !sv || !x)
  {
    fprintf(stderr, "minos_smear: Unable to allocate temporary memory for smearing matrix.");
    return -2;
  }

  /* Copy existing smearing matrix into a non-sparse data structure */
  for (int i=0; i < s->numofbins; i++)
  {
    if (!matrix[i])
    {
      fprintf(stderr, "minos_smear: Incomplete smearing matrix encountered.");
      return -3;
    }
    for (int j=0; j < s->simbins; j++)
    {
      if (j >= lower[i] && j <= upper[i])
        gsl_matrix_set(Sin, i, j, matrix[i][j-lower[i]]);
      else
        gsl_matrix_set(Sin, i, j, 0.0);
    }
  }

  /* Generate matrix that projects a vector with highres_bins entries onto
   * a vector with simbins entries */
//  gsl_matrix_set_zero(R);
//  for (int i=0; i < s->simbins; i++)
//    for (int j=0; j < factor; j++)
//      gsl_matrix_set(R, i, factor*i + j, 1.0);
//  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Sin, R, 0.0, SR);

  /* Generate matrix describing effect of decay pipe length. For numerical stability, we
   * use erfc instead of erf if erf would be too close to unity. Maybe this
   * does not really help, but at least it doesn't hurt. */
  for (int i=0; i < s->simbins; i++)
  {
    for (int j=0; j < s->simbins; j++)
    {
      double E = 0.5 * (glb_upper_sbin_boundary(i,s) + glb_lower_sbin_boundary(i,s)) * GEV;
      double sigma = 0.5 * MIN(1, 7.8045*METER * 4.7*E/M_PI_PLUS) * E/GEV;

      double c = 0.5 * (glb_upper_sbin_boundary(j,s) + glb_lower_sbin_boundary(j,s));
      double x1 = (glb_upper_sbin_boundary(i,s) - c) / (sqrt(2.0)*sigma);
      double x2 = (glb_lower_sbin_boundary(i,s) - c) / (sqrt(2.0)*sigma);
      double erf1 = gsl_sf_erf(x1);
      double erf2 = gsl_sf_erf(x2);
      if (erf1 < 0.5 && erf2 < 0.5)
        gsl_matrix_set(F, i, j, 0.5 * (erf1 - erf2));
      else if (erf1 < 0.5 && erf2 >= 0.5)
        gsl_matrix_set(F, i, j, 0.5 * ((erf1 + gsl_sf_erfc(x2)) - 1.0));
      else if (erf1 >= 0.5 && erf2 < 0.5)
        gsl_matrix_set(F, i, j, 0.5 * ((-erf2 - gsl_sf_erfc(x1)) + 1.0));
      else
        gsl_matrix_set(F, i, j, 0.5 * (gsl_sf_erfc(x2) - gsl_sf_erfc(x1)));
    }
  }

  /* Multiply original smearing matrix S with new contribution F */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Sin, F, 0.0, Sout);

  /* Solve linear system F^T (S R F^-1)^T = (S R)^T to find (S R F^-1) */
//  gsl_matrix_transpose_memcpy(FT, F);
//  gsl_linalg_SV_decomp(FT, V, sv, x);
////  gsl_linalg_SV_decomp_jacobi(FT, V, sv);
//  if (gsl_vector_min(sv) < TOL  ||  gsl_vector_max(sv) / gsl_vector_min(sv) > 1/TOL)
//    glb_error("minos_smear: Unfolding of filter is numerically unstable. "
//              "Reduce filter width.");
//  for (int i=0; i < s->numofbins; i++)
//  {
//    gsl_matrix_get_row(b, SR, i);
//    gsl_linalg_SV_solve(FT, V, sv, b, x);
//    gsl_matrix_set_row(Sout, i, x);
//  }

  /* Copy smearing matrix back into GLoBES data structures */
  for (int i=0; i < s->numofbins; i++)
  {
    /* Comment about TOL: The value value TOL/2.0 makes this criterion
     * equivalent to the one in SmearMatrixC and SmearMatrixA. It is
     * responsible for tiny differences observed between the results for
     * no compensation and compensation with F = id. */
    lower[i] = 0;
    while (lower[i] < s->simbins-1  &&  fabs(gsl_matrix_get(Sout, i, lower[i])) < TOL/2.0)
      lower[i]++;

    upper[i] = s->simbins - 1;
    while (upper[i] > lower[i]  &&  fabs(gsl_matrix_get(Sout, i, upper[i])) < TOL/2.0)
      upper[i]--;

    if (matrix[i])
      free(matrix[i]);
    matrix[i] =  (double*) malloc(sizeof(double)*(upper[i]-lower[i]+1));
    if (!matrix[i])
    {
      fprintf(stderr, "minos_smear: Unable to allocate smearing matrix.");
      return -4;
    }

    for (int j=0; j < upper[i]-lower[i]+1; j++)
      matrix[i][j] = gsl_matrix_get(Sout, i, lower[i]+j);
  }

  gsl_vector_free(x);
  gsl_vector_free(sv);
  gsl_vector_free(b);
  gsl_matrix_free(V);
  gsl_matrix_free(F);
  gsl_matrix_free(SR);
  gsl_matrix_free(R);
  gsl_matrix_free(Sout);
  gsl_matrix_free(Sin);

  return 0;
}


// -------------------------------------------------------------------------
double chiMINOS_2016(int exp, int rule, int n_params, double *x, double *sys_errors,
                     void *user_data)
// -------------------------------------------------------------------------
// Calculate chi^2 for NC or CC events in the MINOS sterile neutrino search.
// The function produces sensible results for n_params=0 (systematics OFF)
// and for n_params=5 (systematics ON). In the latter case, the following
// nuisance parameters are included:
//   x[0]: Overall normalization of F/N ratio
//   x[1]: Error in BG normalization - near
//   x[2]: Error in BG normalization - far
//   x[3]: Signal energy calibration error - far
//   x[4]: BG energy calibration error - far
// user_data should be pointing to an integer 
// -------------------------------------------------------------------------
{
  // far/near ratio for CC selection (https://arxiv.org/abs/1607.01176, fig. 2)
  // and its error bars
  static const double fn_ratio_CC[] = {
    0.19758e-3, 0.08964e-3, 0.03586e-3, 0.05922e-3, 0.06566e-3, 0.06566e-3,
    0.08543e-3, 0.08975e-3, 0.11605e-3, 0.16580e-3, 0.16865e-3, 0.22289e-3,
    0.25490e-3, 0.20259e-3, 0.28451e-3, 0.33535e-3, 0.31520e-3, 0.35768e-3,
    0.35474e-3, 0.28873e-3, 0.33765e-3, 0.36386e-3, 0.38086e-3, 0.42904e-3,
    0.47014e-3, 0.52927e-3, 0.52550e-3, 0.46555e-3, 0.63141e-3, 0.55969e-3,
    0.58148e-3, 0.66194e-3, 0.54921e-3, 0.61302e-3, 0.67628e-3, 0.65567e-3,
    0.59343e-3, 0.45460e-3, 0.40780e-3, 0.48218e-3, 0.45037e-3, 0.50674e-3,
    0.39898 };
  static const double fn_errors_CC[] = {
    0.05843e-3, 0.03839e-3, 0.01981e-3, 0.02367e-3, 0.01787e-3, 0.02073e-3,
    0.01511e-3, 0.0152e-3, 0.02073e-3, 0.02237e-3, 0.02246e-3, 0.02532e-3,
    0.02935e-3, 0.02659e-3, 0.03439e-3, 0.03662e-3, 0.03635e-3, 0.04527e-3,
    0.05033e-3, 0.04517e-3, 0.05547e-3, 0.02973e-3, 0.03396e-3, 0.04086e-3,
    0.04223e-3, 0.0505e-3, 0.05445e-3, 0.05483e-3, 0.0665e-3, 0.06282e-3,
    0.06989e-3, 0.07596e-3, 0.07128e-3, 0.08204e-3, 0.08820e-3, 0.09033e-3,
    0.06294e-3, 0.06358e-3, 0.06349e-3, 0.0749e-3, 0.075820e-3, 0.054190e-3,
    0.0667 };

  // far/near ratio for NC selection (https://arxiv.org/abs/1607.01176e-3, fig. 2)
  // and its error bars
  static const double fn_ratio_NC[] = {
    0.33369e-3, 0.30551e-3, 0.19701e-3, 0.21231e-3, 0.20336e-3, 0.15623e-3,
    0.18407e-3, 0.17201e-3, 0.19943e-3, 0.19116e-3, 0.17649e-3, 0.17098e-3,
    0.19522e-3, 0.24732e-3, 0.18028e-3, 0.19483e-3, 0.20743e-3, 0.21997e-3,
    0.17644e-3, 0.14571e-3, 0.14895e-3, 0.24740e-3, 0.20875e-3, 0.25848e-3,
    0.26785e-3, 0.29789e-3, 0.24030e-3, 0.39489e-3, 0.29755 };
  static const double fn_errors_NC[] = {
    0.02814e-3, 0.04412e-3, 0.03302e-3, 0.03426e-3, 0.03287e-3, 0.0286e-3,
    0.03025e-3, 0.03012e-3, 0.03335e-3, 0.03177e-3, 0.03183e-3, 0.03086e-3,
    0.03336e-3, 0.03879e-3, 0.03536e-3, 0.03562e-3, 0.03928e-3, 0.04307e-3,
    0.04168e-3, 0.03858e-3, 0.04182e-3, 0.01723e-3, 0.02405e-3, 0.03267e-3,
    0.04294e-3, 0.05624e-3, 0.04734e-3, 0.07717e-3, 0.06491 };

  const double *data, *errors;
  if (user_data)
  {
    switch (*((int *) user_data))
    {
      case MINOS_NC:
        data   = fn_ratio_NC;
        errors = fn_errors_NC;
        break;
      case MINOS_CC:
        data   = fn_ratio_CC;
        errors = fn_errors_CC;
        break;
      default:
        return -2e-10;
    }
  }
  else
    return -1e10;

  int exp_near = exp + 1;
  int exp_far  = exp;
  int n_bins = glbGetNumberOfBins(exp);
  double *sig_N = glbGetSignalFitRatePtr(exp_near, rule);
  double *bg_N  = glbGetBGFitRatePtr(exp_near, rule);
  int ew_low, ew_high;
  double emin, emax;
  double chi2 = 0.0;
  int i;

  glbGetEminEmax(exp, &emin, &emax);
  glbGetEnergyWindowBins(exp_far, rule, &ew_low, &ew_high);
  if (n_params)                    // Systematics ON
  {
    double sig_F[n_bins], bg_F[n_bins];
    double norm_tot = 1.0 + x[0];
    double norm_bgN = 1.0 + x[1];
    double norm_bgF = 1.0 + x[2];
    // FIXME Energy scale uncertainty switched off because it leads to
    // convergence problems in the systematics minimization in some cases
    // (in particular at very large th34)
    glbShiftEnergyScale(0.0*x[3], glbGetSignalFitRatePtr(exp_far, rule), sig_F, n_bins, emin, emax);
    glbShiftEnergyScale(0.0*x[4], glbGetBGFitRatePtr(exp_far, rule), bg_F, n_bins, emin, emax);
    for (i=ew_low; i <= ew_high; i++)
    {
      double nf_fit_ratio = norm_tot * (sig_F[i] + norm_bgF*bg_F[i])
                                     / (sig_N[i] + norm_bgN*bg_N[i]);
      chi2 += SQR(data[i] - nf_fit_ratio) / SQR(errors[i]);
    }
  }
  else                             // Systematics OFF
  {
    double *sig_F = glbGetSignalFitRatePtr(exp_far, rule);
    double *bg_F  = glbGetBGFitRatePtr(exp_far, rule);
    for (i=ew_low; i <= ew_high; i++)
    {
      double nf_fit_ratio = (sig_F[i] + bg_F[i]) / (sig_N[i] + bg_N[i]);
      chi2 += SQR(data[i] - nf_fit_ratio) / SQR(errors[i]);
    }
  }

  // Systematical part of chi^2 (= priors)
  for (i=0; i < n_params; i++)
    chi2 += square(x[i] / sys_errors[i]);

  return chi2;
}


// -------------------------------------------------------------------------
double chiMINOS_2011(int exp, int rule, int n_params, double *x, double *errors,
                     void *user_data)
// -------------------------------------------------------------------------
// Calculate chi^2 for NC or CC events in the MINOS sterile neutrino search.
// The function produces sensible results for n_params=0 (systematics OFF)
// and for n_params=5 (systematics ON). In the latter case, the following
// nuisance parameters are included:
//   x[0]: Overall normalization of F/N ratio
//   x[1]: Error in BG normalization - near
//   x[2]: Error in BG normalization - far
//   x[3]: Signal energy calibration error - far
//   x[4]: BG energy calibration error - far
// user_data should be pointing to an integer 
// -------------------------------------------------------------------------
{
  // NC data shown at Neutrino 2010
  static const double data_NC_N[] =  { 23.07692e4, 45.44379e4, 47.04142e4, 41.36095e4,
      28.40237e4, 17.92899e4, 12.78106e4, 9.76331e4, 7.98817e4, 6.56805e4,
       5.50296e4,  4.61538e4,  4.08284e4, 3.55030e4, 3.01775e4, 2.66272e4,
       2.48521e4,  2.13018e4,  1.77515e4, 1.77515e4 };
  static const double data_NC_F[] = { 92, 129, 106, 100, 73, 47, 52, 27, 22, 21, 10, 16,
      12, 8, 9, 10, 2, 8, 8, 2 };

  // CC data set from 1103.0340 - rebinned to 1 GeV bins between 0 and 20 GeV
  static const double data_CC_N[] = { 79474.4, 1.08292e6, 2.75959e6, 2.9265e6,
    1.66333e6, 957837., 722416., 638357., 582196., 526039., 477828., 429627., 393376.,
     361109.,  320876., 280643., 248376., 220098., 199777., 175483.};
  static const double data_CC_F[] = { 8.59002, 28.308, 145.64, 237.007, 175.705,
    94.8805, 86.2907, 89.0238, 92.9283, 89.0239, 86.6811, 64.8156, 68.7202, 51.5401,
    56.2256, 63.2538, 46.8547, 42.9501, 40.6074, 32.7983, 163.991, 109.328 };

  const double *data_N, *data_F;
  if (user_data)
  {
    switch (*((int *) user_data))
    {
      case MINOS_NC:
        data_N = data_NC_N;
        data_F = data_NC_F;
        break;
      case MINOS_CC:
        data_N = data_CC_N;
        data_F = data_CC_F;
        break;
      default:
        return -2e-10;
    }
  }
  else
    return -1e10;

  int exp_near = exp + 1;
  int exp_far  = exp;
  int n_bins = glbGetNumberOfBins(exp);
  double *sig_N = glbGetSignalFitRatePtr(exp_near, rule);
  double *bg_N  = glbGetBGFitRatePtr(exp_near, rule);
  int ew_low, ew_high;
  double emin, emax;
  double fit_rate;
  double chi2   = 0.0;
  int i;

  glbGetEminEmax(exp, &emin, &emax);
  glbGetEnergyWindowBins(exp_far, rule, &ew_low, &ew_high);
  if (n_params)                    // Systematics ON
  {
    double sig_F[n_bins], bg_F[n_bins];
    double norm_tot = 1.0 + x[0];
    double norm_bgN = 1.0 + x[1];
    double norm_bgF = 1.0 + x[2];
    // FIXME Energy scale uncertainty switched off because it leads to
    // convergence problems in the systematics minimization in some cases
    // (in particular at very large th34)
    glbShiftEnergyScale(0.0*x[3], glbGetSignalFitRatePtr(exp_far, rule), sig_F, n_bins, emin, emax);
    glbShiftEnergyScale(0.0*x[4], glbGetBGFitRatePtr(exp_far, rule), bg_F, n_bins, emin, emax);
    for (i=ew_low; i <= ew_high; i++)
    {
      // Predicted rate at far detector is F_th/N_th * N_data
      fit_rate = norm_tot * (sig_F[i] + norm_bgF*bg_F[i])/(sig_N[i] + norm_bgN*bg_N[i]) * data_N[i];
      chi2 += poisson_likelihood(data_F[i], fit_rate);
    }
  }
  else                             // Systematics OFF
  {
    double *sig_F = glbGetSignalFitRatePtr(exp_far, rule);
    double *bg_F  = glbGetBGFitRatePtr(exp_far, rule);
    for (i=ew_low; i <= ew_high; i++)
    {
      // Predicted rate at far detector is F_th/N_th * N_data
      fit_rate = (sig_F[i] + bg_F[i])/(sig_N[i] + bg_N[i]) * data_N[i];
      chi2 += poisson_likelihood(data_F[i], fit_rate);
    }
  }

  // Systematical part of chi^2 (= priors)
  for (i=0; i < n_params; i++)
    chi2 += square(x[i] / errors[i]);

  return chi2;
}


// -------------------------------------------------------------------------
double chiMINOS_2010(int exp, int rule, int n_params, double *x, double *errors,
                     void *user_data)
// -------------------------------------------------------------------------
// Calculate chi^2 for NC or CC events in the MINOS sterile neutrino search.
// The function produces sensible results for n_params=0 (systematics OFF)
// and for n_params=5 (systematics ON). In the latter case, the following
// nuisance parameters are included:
//   x[0]: Overall normalization of F/N ratio
//   x[1]: Error in BG normalization - near
//   x[2]: Error in BG normalization - far
//   x[3]: Signal energy calibration error - far
//   x[4]: BG energy calibration error - far
// user_data should be pointing to an integer 
// -------------------------------------------------------------------------
{
  // NC data set from 1001.0336 //2010
//  double data_NC_N[] =  { 25.58376e4, 29.34010e4, 26.70051e4, 24.06091e4, //2010
//    17.36041e4, 11.06599e4, 7.51269e4, 5.58376e4, 4.46701e4, 3.75635e4, 3.24873e4, //2010
//     2.63959e4,  2.33503e4, 2.03046e4, 1.72589e4, 1.52284e4, 1.31980e4, 1.21827e4, //2010
//     1.01523e4,  0.91371e4 }; //2010
//  double data_NC_F[] = { 42, 58, 41, 52, 40, 23, 22, 19, 11, 12, 8, 4, 7, //2010
//     5, 5, 1, 6, 0, 2, 5 }; //2010

  // NC data shown at Neutrino 2010 //2010
  static const double data_NC_N[] =  { 23.07692e4, 45.44379e4, 47.04142e4, 41.36095e4, //2010
      28.40237e4, 17.92899e4, 12.78106e4, 9.76331e4, 7.98817e4, 6.56805e4, //2010
       5.50296e4,  4.61538e4,  4.08284e4, 3.55030e4, 3.01775e4, 2.66272e4, //2010
       2.48521e4,  2.13018e4,  1.77515e4, 1.77515e4 }; //2010
  static const double data_NC_F[] = { 92, 129, 106, 100, 73, 47, 52, 27, 22, 21, 10, 16, //2010
      12, 8, 9, 10, 2, 8, 8, 2 }; //2010

  // CC data set from 1001.0336 //2010
  static const double data_CC_N[] = { 1.96970e4, 20.77135e4, 50.67493e4, 56.04683e4, //2010
    34.91735e4, 20.59229e4, 14.86226e4, 12.35537e4, 11.10193e4, 10.20661e4, 9.49036e4, //2010
     8.59504e4,  7.87879e4,  7.16253e4,  6.62534e4,  5.73003e4,  5.19284e4, 4.47658e4, //2010
     3.93939e4,  3.58127e4 }; //2010
  static const double data_CC_F[] = { 2, 11, 51, 80, 64, 21, 31, 28, 22, 22, 12, 21, 16, //2010
    20, 19, 17, 20, 11, 14, 7 }; //2010

  const double *data_N, *data_F; //2010
//  const char *rule_name = glbValueToName(exp, "rule", rule); //2010
//  if (strcmp(rule_name, "#rule_NC") == 0) //2010
  if (user_data) //2010
  { //2010
    switch (*((int *) user_data)) //2010
    { //2010
      case MINOS_NC: //2010
        data_N = data_NC_N; //2010
        data_F = data_NC_F; //2010
        break; //2010
      case MINOS_CC: //2010
        data_N = data_CC_N; //2010
        data_F = data_CC_F; //2010
        break; //2010
      default: //2010
        return -2e-10; //2010
    } //2010
  } //2010
  else //2010
    return -1e10; //2010
 //2010
  int exp_near = exp + 1; //2010
  int exp_far  = exp; //2010
  int n_bins = glbGetNumberOfBins(exp); //2010
  double *sig_N = glbGetSignalFitRatePtr(exp_near, rule); //2010
  double *bg_N  = glbGetBGFitRatePtr(exp_near, rule); //2010
  int ew_low, ew_high; //2010
  double emin, emax; //2010
  double fit_rate; //2010
  double chi2   = 0.0; //2010
  int i; //2010
 //2010
  glbGetEminEmax(exp, &emin, &emax); //2010
  glbGetEnergyWindowBins(exp_far, rule, &ew_low, &ew_high); //2010
  if (n_params)                    // Systematics ON //2010
  { //2010
    double sig_F[n_bins], bg_F[n_bins]; //2010
    double norm_tot = 1.0 + x[0]; //2010
    double norm_bgN = 1.0 + x[1]; //2010
    double norm_bgF = 1.0 + x[2]; //2010
    glbShiftEnergyScale(x[3], glbGetSignalFitRatePtr(exp_far, rule), sig_F, n_bins, emin, emax); //2010
    glbShiftEnergyScale(x[4], glbGetBGFitRatePtr(exp_far, rule), bg_F, n_bins, emin, emax); //2010
    for (i=ew_low; i <= ew_high; i++) //2010
    { //2010
      // Predicted rate at far detector is F_th/N_th * N_data //2010
      fit_rate = norm_tot * (sig_F[i] + norm_bgF*bg_F[i])/(sig_N[i] + norm_bgN*bg_N[i]) * data_N[i]; //2010
 //2010
      // NC only: Rescale to pot from 1001.0336 to verify their sensitivities //2010
      if (rule == 0) // Beware: This is error-prone if the rules in the glb file change //2010
        chi2 += poisson_likelihood(3.18/7.2*data_F[i], 3.18/7.2*fit_rate); //2010
      else //2010
        chi2 += poisson_likelihood(data_F[i], fit_rate); //2010
    } //2010
  } //2010
  else                             // Systematics OFF //2010
  { //2010
    double *sig_F = glbGetSignalFitRatePtr(exp_far, rule); //2010
    double *bg_F  = glbGetBGFitRatePtr(exp_far, rule); //2010
    for (i=ew_low; i <= ew_high; i++) //2010
    { //2010
      // Predicted rate at far detector is F_th/N_th * N_data //2010
      fit_rate = (sig_F[i] + bg_F[i])/(sig_N[i] + bg_N[i]) * data_N[i]; //2010
      chi2 += poisson_likelihood(data_F[i], fit_rate); //2010
    } //2010
  } //2010
 //2010
  // Systematical part of chi^2 (= priors) //2010
  for (i=0; i < n_params; i++) //2010
    chi2 += square(x[i] / errors[i]); //2010
 //2010
  return chi2; //2010
} //2010


// -------------------------------------------------------------------------
//   A N A L Y S I S   B A S E D   O N   2 0 1 6   D A T A   R E L E A S E
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
int MINOS_2016_init()
// -------------------------------------------------------------------------
// A new implementation of the MINOS fit, based on the code provided in the
// data release (http://www-numi.fnal.gov/PublicInfo/forscientists.html).
// Bypasses GLoBES except for the computation of the oscillation
// probabilities.
// -------------------------------------------------------------------------
{
  f = new TFile("glb/minos-2016/data/dataRelease.root");
  CoVarCCinvert = (TMatrixD*)f->Get("TotalInvertCC"); assert(CoVarCCinvert);
  CoVarNCinvert = (TMatrixD*)f->Get("TotalInvertNC"); assert(CoVarNCinvert);
  
  f->GetObject("hRecoToTrueNDNCSelectedTrueNC",   NDNC_TrueNC);   assert(NDNC_TrueNC);
  f->GetObject("hRecoToTrueNDNCSelectedNuMu",     NDNC_NuMu);     assert(NDNC_NuMu);
  f->GetObject("hRecoToTrueNDNCSelectedBeamNue",  NDNC_BeamNue);  assert(NDNC_BeamNue);
  f->GetObject("hRecoToTrueNDNCSelectedAppNue",   NDNC_AppNue);   assert(NDNC_AppNue);
  f->GetObject("hRecoToTrueNDNCSelectedAppNuTau", NDNC_AppNuTau); assert(NDNC_AppNuTau);

  f->GetObject("hRecoToTrueNDCCSelectedTrueNC",   NDCC_TrueNC);   assert(NDCC_TrueNC);
  f->GetObject("hRecoToTrueNDCCSelectedNuMu",     NDCC_NuMu);     assert(NDCC_NuMu);
  f->GetObject("hRecoToTrueNDCCSelectedBeamNue",  NDCC_BeamNue);  assert(NDCC_BeamNue);
  f->GetObject("hRecoToTrueNDCCSelectedAppNue",   NDCC_AppNue);   assert(NDCC_AppNue);
  f->GetObject("hRecoToTrueNDCCSelectedAppNuTau", NDCC_AppNuTau); assert(NDCC_AppNuTau);

  f->GetObject("hRecoToTrueFDNCSelectedTrueNC",   FDNC_TrueNC);   assert(FDNC_TrueNC);
  f->GetObject("hRecoToTrueFDNCSelectedNuMu",     FDNC_NuMu);     assert(FDNC_NuMu);
  f->GetObject("hRecoToTrueFDNCSelectedBeamNue",  FDNC_BeamNue);  assert(FDNC_BeamNue);
  f->GetObject("hRecoToTrueFDNCSelectedAppNue",   FDNC_AppNue);   assert(FDNC_AppNue);
  f->GetObject("hRecoToTrueFDNCSelectedAppNuTau", FDNC_AppNuTau); assert(FDNC_AppNuTau);

  f->GetObject("hRecoToTrueFDCCSelectedTrueNC",   FDCC_TrueNC);   assert(FDCC_TrueNC);
  f->GetObject("hRecoToTrueFDCCSelectedNuMu",     FDCC_NuMu);     assert(FDCC_NuMu);
  f->GetObject("hRecoToTrueFDCCSelectedBeamNue",  FDCC_BeamNue);  assert(FDCC_BeamNue);
  f->GetObject("hRecoToTrueFDCCSelectedAppNue",   FDCC_AppNue);   assert(FDCC_AppNue);
  f->GetObject("hRecoToTrueFDCCSelectedAppNuTau", FDCC_AppNuTau); assert(FDCC_AppNuTau);

  f->GetObject("dataFDNC", FD_dataNC); assert(FD_dataNC);
  f->GetObject("dataFDCC", FD_dataCC); assert(FD_dataCC);
  f->GetObject("dataNDNC", ND_dataNC); assert(ND_dataNC);
  f->GetObject("dataNDCC", ND_dataCC); assert(ND_dataCC);

  RdataCC = (TH1D *) FD_dataCC->Clone();
  RdataCC->Divide(ND_dataCC);

  RdataNC = (TH1D *) FD_dataNC->Clone();
  RdataNC->Divide(ND_dataNC);

  // Copy data from ROOT histograms to simple arrays for more efficient access
  int n_bins_true_nd = NDCC_TrueNC->GetXaxis()->GetNbins();
  int n_bins_true_fd = FDCC_TrueNC->GetXaxis()->GetNbins();
  int n_bins_reco_cc = FD_dataCC->GetXaxis()->GetNbins();
  int n_bins_reco_nc = FD_dataNC->GetXaxis()->GetNbins();
  sm_nd_nc_truenc    = new double[n_bins_true_nd * n_bins_reco_nc];  assert(sm_nd_nc_truenc);
  sm_nd_nc_numu      = new double[n_bins_true_nd * n_bins_reco_nc];  assert(sm_nd_nc_numu);
  sm_nd_nc_beamnue   = new double[n_bins_true_nd * n_bins_reco_nc];  assert(sm_nd_nc_beamnue);
  sm_nd_nc_appnue    = new double[n_bins_true_nd * n_bins_reco_nc];  assert(sm_nd_nc_appnue);
  sm_nd_nc_appnutau  = new double[n_bins_true_nd * n_bins_reco_nc];  assert(sm_nd_nc_appnutau);
  sm_fd_nc_truenc    = new double[n_bins_true_fd * n_bins_reco_nc];  assert(sm_fd_nc_truenc);
  sm_fd_nc_numu      = new double[n_bins_true_fd * n_bins_reco_nc];  assert(sm_fd_nc_numu);
  sm_fd_nc_beamnue   = new double[n_bins_true_fd * n_bins_reco_nc];  assert(sm_fd_nc_beamnue);
  sm_fd_nc_appnue    = new double[n_bins_true_fd * n_bins_reco_nc];  assert(sm_fd_nc_appnue);
  sm_fd_nc_appnutau  = new double[n_bins_true_fd * n_bins_reco_nc];  assert(sm_fd_nc_appnutau);
  sm_nd_cc_truenc    = new double[n_bins_true_nd * n_bins_reco_cc];  assert(sm_nd_cc_truenc);
  sm_nd_cc_numu      = new double[n_bins_true_nd * n_bins_reco_cc];  assert(sm_nd_cc_numu);
  sm_nd_cc_beamnue   = new double[n_bins_true_nd * n_bins_reco_cc];  assert(sm_nd_cc_beamnue);
  sm_nd_cc_appnue    = new double[n_bins_true_nd * n_bins_reco_cc];  assert(sm_nd_cc_appnue);
  sm_nd_cc_appnutau  = new double[n_bins_true_nd * n_bins_reco_cc];  assert(sm_nd_cc_appnutau);
  sm_fd_cc_truenc    = new double[n_bins_true_fd * n_bins_reco_cc];  assert(sm_fd_cc_truenc);
  sm_fd_cc_numu      = new double[n_bins_true_fd * n_bins_reco_cc];  assert(sm_fd_cc_numu);
  sm_fd_cc_beamnue   = new double[n_bins_true_fd * n_bins_reco_cc];  assert(sm_fd_cc_beamnue);
  sm_fd_cc_appnue    = new double[n_bins_true_fd * n_bins_reco_cc];  assert(sm_fd_cc_appnue);
  sm_fd_cc_appnutau  = new double[n_bins_true_fd * n_bins_reco_cc];  assert(sm_fd_cc_appnutau);
  int k=0;
  for (int i=0; i < n_bins_reco_nc; i++)
    for (int j=0; j < n_bins_true_nd; j++)
    {
      sm_nd_nc_truenc[k]   = NDNC_TrueNC->GetBinContent(j+1, i+1);
      sm_nd_nc_numu[k]     = NDNC_NuMu->GetBinContent(j+1, i+1);
      sm_nd_nc_beamnue[k]  = NDNC_BeamNue->GetBinContent(j+1, i+1);
      sm_nd_nc_appnue[k]   = NDNC_AppNue->GetBinContent(j+1, i+1);
      sm_nd_nc_appnutau[k] = NDNC_AppNuTau->GetBinContent(j+1, i+1);
      k++;
    }
  k = 0;
  for (int i=0; i < n_bins_reco_nc; i++)
    for (int j=0; j < n_bins_true_fd; j++)
    {
      sm_fd_nc_truenc[k]   = FDNC_TrueNC->GetBinContent(j+1, i+1);
      sm_fd_nc_numu[k]     = FDNC_NuMu->GetBinContent(j+1, i+1);
      sm_fd_nc_beamnue[k]  = FDNC_BeamNue->GetBinContent(j+1, i+1);
      sm_fd_nc_appnue[k]   = FDNC_AppNue->GetBinContent(j+1, i+1);
      sm_fd_nc_appnutau[k] = FDNC_AppNuTau->GetBinContent(j+1, i+1);
      k++;
    }
  k = 0;
  for (int i=0; i < n_bins_reco_cc; i++)
    for (int j=0; j < n_bins_true_nd; j++)
    {
      sm_nd_cc_truenc[k]   = NDCC_TrueNC->GetBinContent(j+1, i+1);
      sm_nd_cc_numu[k]     = NDCC_NuMu->GetBinContent(j+1, i+1);
      sm_nd_cc_beamnue[k]  = NDCC_BeamNue->GetBinContent(j+1, i+1);
      sm_nd_cc_appnue[k]   = NDCC_AppNue->GetBinContent(j+1, i+1);
      sm_nd_cc_appnutau[k] = NDCC_AppNuTau->GetBinContent(j+1, i+1);
      k++;
    }
  k = 0;
  for (int i=0; i < n_bins_reco_cc; i++)
    for (int j=0; j < n_bins_true_fd; j++)
    {
      sm_fd_cc_truenc[k]   = FDCC_TrueNC->GetBinContent(j+1, i+1);
      sm_fd_cc_numu[k]     = FDCC_NuMu->GetBinContent(j+1, i+1);
      sm_fd_cc_beamnue[k]  = FDCC_BeamNue->GetBinContent(j+1, i+1);
      sm_fd_cc_appnue[k]   = FDCC_AppNue->GetBinContent(j+1, i+1);
      sm_fd_cc_appnutau[k] = FDCC_AppNuTau->GetBinContent(j+1, i+1);
      k++;
    }
  data_r_nc = RdataNC->GetArray() + 1;
  data_r_cc = RdataCC->GetArray() + 1;

  inv_cov_matrix_nc = CoVarNCinvert->GetMatrixArray();
  inv_cov_matrix_cc = CoVarCCinvert->GetMatrixArray();

  // Compute expected unoscillated rates at the ND for the ND normalization
  // penalty term
  double pred_nd_nc_unosc=0.0, pred_nd_cc_unosc=0.0;

  for (int i=0; i < n_bins_reco_nc; i++)
    for (int j=0; j < n_bins_true_nd; j++)
      pred_nd_nc_unosc += NDNC_TrueNC->GetBinContent(j+1,i+1)
                        + NDNC_NuMu->GetBinContent(j+1,i+1)
                        + NDNC_BeamNue->GetBinContent(j+1,i+1)
                        + NDNC_AppNue->GetBinContent(j+1,i+1)
                        + NDNC_AppNuTau->GetBinContent(j+1,i+1);
  for (int i=0; i < n_bins_reco_cc; i++)
    for (int j=0; j < n_bins_true_nd; j++)
      pred_nd_cc_unosc += NDCC_TrueNC->GetBinContent(j+1,i+1)
                        + NDCC_NuMu->GetBinContent(j+1,i+1)
                        + NDCC_BeamNue->GetBinContent(j+1,i+1)
                        + NDCC_AppNue->GetBinContent(j+1,i+1)
                        + NDCC_AppNuTau->GetBinContent(j+1,i+1);

  sigma_nd_norm_nc = 0.5 * pred_nd_nc_unosc;
  sigma_nd_norm_cc = 0.5 * pred_nd_cc_unosc;

  data_total_nd_nc = ND_dataNC->Integral(1, ND_dataNC->GetXaxis()->GetNbins());
  data_total_nd_cc = ND_dataCC->Integral(1, ND_dataCC->GetXaxis()->GetNbins());

  return 0;
}


// -------------------------------------------------------------------------
double MINOS_2016_prior(const glb_params params)
// -------------------------------------------------------------------------
// A new implementation of the MINOS fit, based on the code provided in the
// data release (http://www-numi.fnal.gov/PublicInfo/forscientists.html).
// Bypasses GLoBES except for the computation of the oscillation
// probabilities.
// -------------------------------------------------------------------------
{
  const double L_nd   = 1.04;     // km
  const double L_fd   = 735.0;    // km
  const double rho_nd = 0.0;      // grams/cm^3
  const double rho_fd = 2.8;      // grams/cm^3

  int n_bins_true_nd = NDCC_TrueNC->GetXaxis()->GetNbins();
  int n_bins_true_fd = FDCC_TrueNC->GetXaxis()->GetNbins();
  int n_bins_reco_cc = FD_dataCC->GetXaxis()->GetNbins();
  int n_bins_reco_nc = FD_dataNC->GetXaxis()->GetNbins();
  double P_nd_truenc[n_bins_true_nd], P_nd_numu[n_bins_true_nd], P_nd_beamnue[n_bins_true_nd];
  double P_nd_appnue[n_bins_true_nd], P_nd_appnutau[n_bins_true_nd];
  double P_fd_truenc[n_bins_true_fd], P_fd_numu[n_bins_true_fd], P_fd_beamnue[n_bins_true_fd];
  double P_fd_appnue[n_bins_true_fd], P_fd_appnutau[n_bins_true_fd];
  double pred_nd_nc[n_bins_reco_nc], pred_nd_cc[n_bins_reco_cc];
  double pred_fd_nc[n_bins_reco_nc], pred_fd_cc[n_bins_reco_cc];
  double pred_r_nc[n_bins_reco_nc], pred_r_cc[n_bins_reco_cc];

  // Compute oscillation probabilities  
  snu_set_oscillation_parameters(params, NULL);
  for (int j=0; j < n_bins_true_nd; j++)
  {
    double P_nd[3][3];
    snu_probability_matrix(P_nd, +1, L_nd / NDCC_TrueNC->GetXaxis()->GetBinCenter(j+1),
                           1, &L_nd, &rho_nd, 0.3, NULL);
    P_nd_truenc[j]   = P_nd[NU_MU][NU_E] + P_nd[NU_MU][NU_MU] + P_nd[NU_MU][NU_TAU];
    P_nd_numu[j]     = P_nd[NU_MU][NU_MU];
    P_nd_beamnue[j]  = P_nd[NU_E][NU_E];
    P_nd_appnue[j]   = P_nd[NU_MU][NU_E];
    P_nd_appnutau[j] = P_nd[NU_MU][NU_TAU];
      // FIXME More oscillation channels?
  }
  for (int j=0; j < n_bins_true_fd; j++)
  {
    double P_fd[3][3];
    snu_probability_matrix(P_fd, +1, L_fd / FDCC_TrueNC->GetXaxis()->GetBinCenter(j+1),
                           1, &L_fd, &rho_fd, 0.3, NULL);
    P_fd_truenc[j]   = P_fd[NU_MU][NU_E] + P_fd[NU_MU][NU_MU] + P_fd[NU_MU][NU_TAU];
    P_fd_numu[j]     = P_fd[NU_MU][NU_MU];
    P_fd_beamnue[j]  = P_fd[NU_E][NU_E];
    P_fd_appnue[j]   = P_fd[NU_MU][NU_E];
    P_fd_appnutau[j] = P_fd[NU_MU][NU_TAU];
      // FIXME More oscillation channels?
  }

  // Compute predicted event rates
  int k=0;
  for (int i=0; i < n_bins_reco_nc; i++)
  {
    pred_nd_nc[i] = 0.0;
    for (int j=0; j < n_bins_true_nd; j++)
    {
      pred_nd_nc[i] += sm_nd_nc_truenc[k] * P_nd_truenc[j]
                     + sm_nd_nc_numu[k] * P_nd_numu[j]
                     + sm_nd_nc_beamnue[k] * P_nd_beamnue[j]
                     + sm_nd_nc_appnue[k] * P_nd_appnue[j]
                     + sm_nd_nc_appnutau[k] * P_nd_appnutau[j];
      k++;
    }
  }
  k = 0;
  for (int i=0; i < n_bins_reco_nc; i++)
  {
    pred_fd_nc[i] = 0.0;
    for (int j=0; j < n_bins_true_fd; j++)
    {
      pred_fd_nc[i] += sm_fd_nc_truenc[k] * P_fd_truenc[j]
                     + sm_fd_nc_numu[k] * P_fd_numu[j]
                     + sm_fd_nc_beamnue[k] * P_fd_beamnue[j]
                     + sm_fd_nc_appnue[k] * P_fd_appnue[j]
                     + sm_fd_nc_appnutau[k] * P_fd_appnutau[j];
      k++;
    }
    pred_r_nc[i] = pred_fd_nc[i] / pred_nd_nc[i];
  }
  k = 0;
  for (int i=0; i < n_bins_reco_cc; i++)
  {
    pred_nd_cc[i] = 0.0;
    for (int j=0; j < n_bins_true_nd; j++)
    {
      pred_nd_cc[i] += sm_nd_cc_truenc[k] * P_nd_truenc[j]
                     + sm_nd_cc_numu[k] * P_nd_numu[j]
                     + sm_nd_cc_beamnue[k] * P_nd_beamnue[j]
                     + sm_nd_cc_appnue[k] * P_nd_appnue[j]
                     + sm_nd_cc_appnutau[k] * P_nd_appnutau[j];
      k++;
    }
  }
  k = 0;
  for (int i=0; i < n_bins_reco_cc; i++)
  {
    pred_fd_cc[i] = 0.0;
    for (int j=0; j < n_bins_true_fd; j++)
    {
      pred_fd_cc[i] += sm_fd_cc_truenc[k] * P_fd_truenc[j]
                     + sm_fd_cc_numu[k] * P_fd_numu[j]
                     + sm_fd_cc_beamnue[k] * P_fd_beamnue[j]
                     + sm_fd_cc_appnue[k] * P_fd_appnue[j]
                     + sm_fd_cc_appnutau[k] * P_fd_appnutau[j];
      k++;
    }
    pred_r_cc[i] = pred_fd_cc[i] / pred_nd_cc[i];
  }

  // Compute chi^2
  double diff_nc[n_bins_reco_nc], diff_cc[n_bins_reco_cc];
  for (int i=0; i < n_bins_reco_nc; i++)
    diff_nc[i] = data_r_nc[i] - pred_r_nc[i];
  for (int i=0; i < n_bins_reco_cc; i++)
    diff_cc[i] = data_r_cc[i] - pred_r_cc[i];

  double chi2_nc = 0.0;
  double chi2_cc = 0.0;
  int l=0;
  for (int i=0; i < n_bins_reco_nc; i++)
  {
    double t = 0.0;
    for (int k=0; k < n_bins_reco_nc; k++)
      t += inv_cov_matrix_nc[l++] * diff_nc[k];
    chi2_nc += diff_nc[i] * t;
  }
  l = 0;
  for (int i=0; i < n_bins_reco_cc; i++)
  {
    double t = 0.0;
    for (int k=0; k < n_bins_reco_cc; k++)
      t += inv_cov_matrix_cc[l++] * diff_cc[k];
    chi2_cc += diff_cc[i] * t;
  }
  double chi2 = chi2_nc + chi2_cc;

  // Add penalty if ND normalization is too far off
  double pred_nd_nc_total = 0.0;
  double pred_nd_cc_total = 0.0;
  for (int i=0; i < n_bins_reco_nc; i++)
    pred_nd_nc_total += pred_nd_nc[i];
  for (int i=0; i < n_bins_reco_cc; i++)
    pred_nd_cc_total += pred_nd_cc[i];
  chi2 += SQR(pred_nd_nc_total - data_total_nd_nc)/SQR(sigma_nd_norm_nc);
  chi2 += SQR(pred_nd_cc_total - data_total_nd_cc)/SQR(sigma_nd_norm_cc);

  return chi2;
}


// -------------------------------------------------------------------------
//   C O M B I N E D   M I N O S / M I N O S +   A N A L Y S I S   2 0 1 7
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
int MINOS_2017_init()
// -------------------------------------------------------------------------
// Initialization of the MINOS/MINOS+ 2017 fit, based on the code provided
// in the data release (https://arxiv.org/abs/1710.06488)
// Bypasses GLoBES except for the computation of the oscillation
// probabilities.
// -------------------------------------------------------------------------
{
#ifdef MINOS_2017_MINIMIZE
  printf("#   MINOS/MINOS+ minimizing over nuisance parameters\n");
#else
  printf("#   MINOS/MINOS+ no minimization over nuisance parameters\n");
#endif

#ifdef MINOS_2017_FREE_NORM
  printf("#   MINOS/MINOS+ using free normalization\n");
#else
  printf("#   MINOS/MINOS+ using fixed normalization\n");
#endif

#ifdef MINOS_2017_USE_GLOBES_LOWPASS
  printf("#   MINOS/MINOS+ antialiasing using GLoBES low-pass filter\n");
#else
  printf("#   MINOS/MINOS+ antialiasing using MINOS method\n");
#endif

  f = new TFile("glb/minos-plus/data/dataRelease.root");
  CoVarCC_relative = (TMatrixD*)f->Get("TotalCCCovar"); assert(CoVarCC_relative);
  CoVarNC_relative = (TMatrixD*)f->Get("TotalNCCovar"); assert(CoVarNC_relative);

  // Read MC histograms for MINOS
  f->GetObject("hRecoToTrueNDNCSelectedTrueNC_minos",   NDNC_TrueNC_minos);   assert(NDNC_TrueNC_minos);
  f->GetObject("hRecoToTrueNDNCSelectedNuMu_minos",     NDNC_NuMu_minos);     assert(NDNC_NuMu_minos);
  f->GetObject("hRecoToTrueNDNCSelectedBeamNue_minos",  NDNC_BeamNue_minos);  assert(NDNC_BeamNue_minos);
  f->GetObject("hRecoToTrueNDNCSelectedAppNue_minos",   NDNC_AppNue_minos);   assert(NDNC_AppNue_minos);
  f->GetObject("hRecoToTrueNDNCSelectedAppNuTau_minos", NDNC_AppNuTau_minos); assert(NDNC_AppNuTau_minos);

  f->GetObject("hRecoToTrueNDCCSelectedTrueNC_minos",   NDCC_TrueNC_minos);   assert(NDCC_TrueNC_minos);
  f->GetObject("hRecoToTrueNDCCSelectedNuMu_minos",     NDCC_NuMu_minos);     assert(NDCC_NuMu_minos);
  f->GetObject("hRecoToTrueNDCCSelectedBeamNue_minos",  NDCC_BeamNue_minos);  assert(NDCC_BeamNue_minos);
  f->GetObject("hRecoToTrueNDCCSelectedAppNue_minos",   NDCC_AppNue_minos);   assert(NDCC_AppNue_minos);
  f->GetObject("hRecoToTrueNDCCSelectedAppNuTau_minos", NDCC_AppNuTau_minos); assert(NDCC_AppNuTau_minos);

  f->GetObject("hRecoToTrueFDNCSelectedTrueNC_minos",   FDNC_TrueNC_minos);   assert(FDNC_TrueNC_minos);
  f->GetObject("hRecoToTrueFDNCSelectedNuMu_minos",     FDNC_NuMu_minos);     assert(FDNC_NuMu_minos);
  f->GetObject("hRecoToTrueFDNCSelectedBeamNue_minos",  FDNC_BeamNue_minos);  assert(FDNC_BeamNue_minos);
  f->GetObject("hRecoToTrueFDNCSelectedAppNue_minos",   FDNC_AppNue_minos);   assert(FDNC_AppNue_minos);
  f->GetObject("hRecoToTrueFDNCSelectedAppNuTau_minos", FDNC_AppNuTau_minos); assert(FDNC_AppNuTau_minos);

  f->GetObject("hRecoToTrueFDCCSelectedTrueNC_minos",   FDCC_TrueNC_minos);   assert(FDCC_TrueNC_minos);
  f->GetObject("hRecoToTrueFDCCSelectedNuMu_minos",     FDCC_NuMu_minos);     assert(FDCC_NuMu_minos);
  f->GetObject("hRecoToTrueFDCCSelectedBeamNue_minos",  FDCC_BeamNue_minos);  assert(FDCC_BeamNue_minos);
  f->GetObject("hRecoToTrueFDCCSelectedAppNue_minos",   FDCC_AppNue_minos);   assert(FDCC_AppNue_minos);
  f->GetObject("hRecoToTrueFDCCSelectedAppNuTau_minos", FDCC_AppNuTau_minos); assert(FDCC_AppNuTau_minos);

  // Read MC histograms for MINOS+
  f->GetObject("hRecoToTrueNDNCSelectedTrueNC_minosPlus",   NDNC_TrueNC_mplus);   assert(NDNC_TrueNC_mplus);
  f->GetObject("hRecoToTrueNDNCSelectedNuMu_minosPlus",     NDNC_NuMu_mplus);     assert(NDNC_NuMu_mplus);
  f->GetObject("hRecoToTrueNDNCSelectedBeamNue_minosPlus",  NDNC_BeamNue_mplus);  assert(NDNC_BeamNue_mplus);
  f->GetObject("hRecoToTrueNDNCSelectedAppNue_minosPlus",   NDNC_AppNue_mplus);   assert(NDNC_AppNue_mplus);
  f->GetObject("hRecoToTrueNDNCSelectedAppNuTau_minosPlus", NDNC_AppNuTau_mplus); assert(NDNC_AppNuTau_mplus);

  f->GetObject("hRecoToTrueNDCCSelectedTrueNC_minosPlus",   NDCC_TrueNC_mplus);   assert(NDCC_TrueNC_mplus);
  f->GetObject("hRecoToTrueNDCCSelectedNuMu_minosPlus",     NDCC_NuMu_mplus);     assert(NDCC_NuMu_mplus);
  f->GetObject("hRecoToTrueNDCCSelectedBeamNue_minosPlus",  NDCC_BeamNue_mplus);  assert(NDCC_BeamNue_mplus);
  f->GetObject("hRecoToTrueNDCCSelectedAppNue_minosPlus",   NDCC_AppNue_mplus);   assert(NDCC_AppNue_mplus);
  f->GetObject("hRecoToTrueNDCCSelectedAppNuTau_minosPlus", NDCC_AppNuTau_mplus); assert(NDCC_AppNuTau_mplus);

  f->GetObject("hRecoToTrueFDNCSelectedTrueNC_minosPlus",   FDNC_TrueNC_mplus);   assert(FDNC_TrueNC_mplus);
  f->GetObject("hRecoToTrueFDNCSelectedNuMu_minosPlus",     FDNC_NuMu_mplus);     assert(FDNC_NuMu_mplus);
  f->GetObject("hRecoToTrueFDNCSelectedBeamNue_minosPlus",  FDNC_BeamNue_mplus);  assert(FDNC_BeamNue_mplus);
  f->GetObject("hRecoToTrueFDNCSelectedAppNue_minosPlus",   FDNC_AppNue_mplus);   assert(FDNC_AppNue_mplus);
  f->GetObject("hRecoToTrueFDNCSelectedAppNuTau_minosPlus", FDNC_AppNuTau_mplus); assert(FDNC_AppNuTau_mplus);

  f->GetObject("hRecoToTrueFDCCSelectedTrueNC_minosPlus",   FDCC_TrueNC_mplus);   assert(FDCC_TrueNC_mplus);
  f->GetObject("hRecoToTrueFDCCSelectedNuMu_minosPlus",     FDCC_NuMu_mplus);     assert(FDCC_NuMu_mplus);
  f->GetObject("hRecoToTrueFDCCSelectedBeamNue_minosPlus",  FDCC_BeamNue_mplus);  assert(FDCC_BeamNue_mplus);
  f->GetObject("hRecoToTrueFDCCSelectedAppNue_minosPlus",   FDCC_AppNue_mplus);   assert(FDCC_AppNue_mplus);
  f->GetObject("hRecoToTrueFDCCSelectedAppNuTau_minosPlus", FDCC_AppNuTau_mplus); assert(FDCC_AppNuTau_mplus);

  // Read data histograms
  f->GetObject("dataFDNC_minos", FD_dataNC_minos); assert(FD_dataNC_minos);
  f->GetObject("dataFDCC_minos", FD_dataCC_minos); assert(FD_dataCC_minos);
  f->GetObject("dataNDNC_minos", ND_dataNC_minos); assert(ND_dataNC_minos);
  f->GetObject("dataNDCC_minos", ND_dataCC_minos); assert(ND_dataCC_minos);

  f->GetObject("dataFDNC_minosPlus", FD_dataNC_mplus); assert(FD_dataNC_mplus);
  f->GetObject("dataFDCC_minosPlus", FD_dataCC_mplus); assert(FD_dataCC_mplus);
  f->GetObject("dataNDNC_minosPlus", ND_dataNC_mplus); assert(ND_dataNC_mplus);
  f->GetObject("dataNDCC_minosPlus", ND_dataCC_mplus); assert(ND_dataCC_mplus);

  // Read histograms for beam miscalibration parameters
  CC_IMiscal_ND_minos = (TH1D*)f->Get("HornCurrentMiscalibration_weights_ND_CC_minos"); assert(CC_IMiscal_ND_minos);
  CC_IDist_ND_minos   = (TH1D*)f->Get("HornCurrentDistribution_weights_ND_CC_minos");   assert(CC_IDist_ND_minos);
  CC_IMiscal_FD_minos = (TH1D*)f->Get("HornCurrentMiscalibration_weights_FD_CC_minos"); assert(CC_IMiscal_FD_minos);
  CC_IDist_FD_minos   = (TH1D*)f->Get("HornCurrentDistribution_weights_FD_CC_minos");   assert(CC_IDist_FD_minos);

  NC_IMiscal_ND_minos = (TH1D*)f->Get("HornCurrentMiscalibration_weights_ND_NC_minos"); assert(NC_IMiscal_ND_minos);
  NC_IDist_ND_minos   = (TH1D*)f->Get("HornCurrentDistribution_weights_ND_NC_minos");   assert(NC_IDist_ND_minos);
  NC_IMiscal_FD_minos = (TH1D*)f->Get("HornCurrentMiscalibration_weights_FD_NC_minos"); assert(NC_IMiscal_FD_minos);
  NC_IDist_FD_minos   = (TH1D*)f->Get("HornCurrentDistribution_weights_FD_NC_minos");   assert(NC_IDist_FD_minos);

  CC_IMiscal_ND_mplus = (TH1D*)f->Get("HornCurrentMiscalibration_weights_ND_CC_minosPlus"); assert(CC_IMiscal_ND_mplus);
  CC_IDist_ND_mplus   = (TH1D*)f->Get("HornCurrentDistribution_weights_ND_CC_minosPlus");   assert(CC_IDist_ND_mplus);
  CC_IMiscal_FD_mplus = (TH1D*)f->Get("HornCurrentMiscalibration_weights_FD_CC_minosPlus"); assert(CC_IMiscal_FD_mplus);
  CC_IDist_FD_mplus   = (TH1D*)f->Get("HornCurrentDistribution_weights_FD_CC_minosPlus");   assert(CC_IDist_FD_mplus);

  NC_IMiscal_ND_mplus = (TH1D*)f->Get("HornCurrentMiscalibration_weights_ND_NC_minosPlus"); assert(NC_IMiscal_ND_mplus);
  NC_IDist_ND_mplus   = (TH1D*)f->Get("HornCurrentDistribution_weights_ND_NC_minosPlus");   assert(NC_IDist_ND_mplus);
  NC_IMiscal_FD_mplus = (TH1D*)f->Get("HornCurrentMiscalibration_weights_FD_NC_minosPlus"); assert(NC_IMiscal_FD_mplus);
  NC_IDist_FD_mplus   = (TH1D*)f->Get("HornCurrentDistribution_weights_FD_NC_minosPlus");   assert(NC_IDist_FD_mplus);


  // Copy data from ROOT histograms to simple arrays for more efficient access
  // -------------------------------------------------------------------------

  // Get binning
  int n_bins_true_nd = NDCC_TrueNC_minos->GetXaxis()->GetNbins(); // exploit that true binning is
  int n_bins_true_fd = FDCC_TrueNC_minos->GetXaxis()->GetNbins(); // the same for all FD and all ND samples
  int n_bins_reco_nd_cc = ND_dataCC_minos->GetXaxis()->GetNbins();
  int n_bins_reco_nd_nc = ND_dataNC_minos->GetXaxis()->GetNbins();
  int n_bins_reco_fd_cc = FD_dataCC_minos->GetXaxis()->GetNbins();
  int n_bins_reco_fd_nc = FD_dataNC_minos->GetXaxis()->GetNbins();

  // Vectors for predicted event rates
  pred_nd_nc_minos = new double[n_bins_reco_nd_nc];
  pred_nd_cc_minos = new double[n_bins_reco_nd_cc];
  pred_fd_nc_minos = new double[n_bins_reco_fd_nc];
  pred_fd_cc_minos = new double[n_bins_reco_fd_cc];
  pred_nd_nc_mplus = new double[n_bins_reco_nd_nc];
  pred_nd_cc_mplus = new double[n_bins_reco_nd_cc];
  pred_fd_nc_mplus = new double[n_bins_reco_fd_nc];
  pred_fd_cc_mplus = new double[n_bins_reco_fd_cc];
  pred_nc = gsl_vector_alloc(n_bins_reco_fd_nc + n_bins_reco_nd_nc);
  pred_cc = gsl_vector_alloc(n_bins_reco_fd_cc + n_bins_reco_nd_cc);
  temp_nc = gsl_vector_alloc(n_bins_reco_fd_nc + n_bins_reco_nd_nc);
  temp_cc = gsl_vector_alloc(n_bins_reco_fd_cc + n_bins_reco_nd_cc);

  // Vectors for smearing matrices
  sm_nd_nc_truenc_minos    = new double[n_bins_true_nd * n_bins_reco_nd_nc];  assert(sm_nd_nc_truenc_minos);
  sm_nd_nc_numu_minos      = new double[n_bins_true_nd * n_bins_reco_nd_nc];  assert(sm_nd_nc_numu_minos);
  sm_nd_nc_beamnue_minos   = new double[n_bins_true_nd * n_bins_reco_nd_nc];  assert(sm_nd_nc_beamnue_minos);
  sm_nd_nc_appnue_minos    = new double[n_bins_true_nd * n_bins_reco_nd_nc];  assert(sm_nd_nc_appnue_minos);
  sm_nd_nc_appnutau_minos  = new double[n_bins_true_nd * n_bins_reco_nd_nc];  assert(sm_nd_nc_appnutau_minos);
  sm_fd_nc_truenc_minos    = new double[n_bins_true_fd * n_bins_reco_fd_nc];  assert(sm_fd_nc_truenc_minos);
  sm_fd_nc_numu_minos      = new double[n_bins_true_fd * n_bins_reco_fd_nc];  assert(sm_fd_nc_numu_minos);
  sm_fd_nc_beamnue_minos   = new double[n_bins_true_fd * n_bins_reco_fd_nc];  assert(sm_fd_nc_beamnue_minos);
  sm_fd_nc_appnue_minos    = new double[n_bins_true_fd * n_bins_reco_fd_nc];  assert(sm_fd_nc_appnue_minos);
  sm_fd_nc_appnutau_minos  = new double[n_bins_true_fd * n_bins_reco_fd_nc];  assert(sm_fd_nc_appnutau_minos);
  sm_nd_cc_truenc_minos    = new double[n_bins_true_nd * n_bins_reco_nd_cc];  assert(sm_nd_cc_truenc_minos);
  sm_nd_cc_numu_minos      = new double[n_bins_true_nd * n_bins_reco_nd_cc];  assert(sm_nd_cc_numu_minos);
  sm_nd_cc_beamnue_minos   = new double[n_bins_true_nd * n_bins_reco_nd_cc];  assert(sm_nd_cc_beamnue_minos);
  sm_nd_cc_appnue_minos    = new double[n_bins_true_nd * n_bins_reco_nd_cc];  assert(sm_nd_cc_appnue_minos);
  sm_nd_cc_appnutau_minos  = new double[n_bins_true_nd * n_bins_reco_nd_cc];  assert(sm_nd_cc_appnutau_minos);
  sm_fd_cc_truenc_minos    = new double[n_bins_true_fd * n_bins_reco_fd_cc];  assert(sm_fd_cc_truenc_minos);
  sm_fd_cc_numu_minos      = new double[n_bins_true_fd * n_bins_reco_fd_cc];  assert(sm_fd_cc_numu_minos);
  sm_fd_cc_beamnue_minos   = new double[n_bins_true_fd * n_bins_reco_fd_cc];  assert(sm_fd_cc_beamnue_minos);
  sm_fd_cc_appnue_minos    = new double[n_bins_true_fd * n_bins_reco_fd_cc];  assert(sm_fd_cc_appnue_minos);
  sm_fd_cc_appnutau_minos  = new double[n_bins_true_fd * n_bins_reco_fd_cc];  assert(sm_fd_cc_appnutau_minos);

  sm_nd_nc_truenc_mplus    = new double[n_bins_true_nd * n_bins_reco_nd_nc];  assert(sm_nd_nc_truenc_mplus);
  sm_nd_nc_numu_mplus      = new double[n_bins_true_nd * n_bins_reco_nd_nc];  assert(sm_nd_nc_numu_mplus);
  sm_nd_nc_beamnue_mplus   = new double[n_bins_true_nd * n_bins_reco_nd_nc];  assert(sm_nd_nc_beamnue_mplus);
  sm_nd_nc_appnue_mplus    = new double[n_bins_true_nd * n_bins_reco_nd_nc];  assert(sm_nd_nc_appnue_mplus);
  sm_nd_nc_appnutau_mplus  = new double[n_bins_true_nd * n_bins_reco_nd_nc];  assert(sm_nd_nc_appnutau_mplus);
  sm_fd_nc_truenc_mplus    = new double[n_bins_true_fd * n_bins_reco_fd_nc];  assert(sm_fd_nc_truenc_mplus);
  sm_fd_nc_numu_mplus      = new double[n_bins_true_fd * n_bins_reco_fd_nc];  assert(sm_fd_nc_numu_mplus);
  sm_fd_nc_beamnue_mplus   = new double[n_bins_true_fd * n_bins_reco_fd_nc];  assert(sm_fd_nc_beamnue_mplus);
  sm_fd_nc_appnue_mplus    = new double[n_bins_true_fd * n_bins_reco_fd_nc];  assert(sm_fd_nc_appnue_mplus);
  sm_fd_nc_appnutau_mplus  = new double[n_bins_true_fd * n_bins_reco_fd_nc];  assert(sm_fd_nc_appnutau_mplus);
  sm_nd_cc_truenc_mplus    = new double[n_bins_true_nd * n_bins_reco_nd_cc];  assert(sm_nd_cc_truenc_mplus);
  sm_nd_cc_numu_mplus      = new double[n_bins_true_nd * n_bins_reco_nd_cc];  assert(sm_nd_cc_numu_mplus);
  sm_nd_cc_beamnue_mplus   = new double[n_bins_true_nd * n_bins_reco_nd_cc];  assert(sm_nd_cc_beamnue_mplus);
  sm_nd_cc_appnue_mplus    = new double[n_bins_true_nd * n_bins_reco_nd_cc];  assert(sm_nd_cc_appnue_mplus);
  sm_nd_cc_appnutau_mplus  = new double[n_bins_true_nd * n_bins_reco_nd_cc];  assert(sm_nd_cc_appnutau_mplus);
  sm_fd_cc_truenc_mplus    = new double[n_bins_true_fd * n_bins_reco_fd_cc];  assert(sm_fd_cc_truenc_mplus);
  sm_fd_cc_numu_mplus      = new double[n_bins_true_fd * n_bins_reco_fd_cc];  assert(sm_fd_cc_numu_mplus);
  sm_fd_cc_beamnue_mplus   = new double[n_bins_true_fd * n_bins_reco_fd_cc];  assert(sm_fd_cc_beamnue_mplus);
  sm_fd_cc_appnue_mplus    = new double[n_bins_true_fd * n_bins_reco_fd_cc];  assert(sm_fd_cc_appnue_mplus);
  sm_fd_cc_appnutau_mplus  = new double[n_bins_true_fd * n_bins_reco_fd_cc];  assert(sm_fd_cc_appnutau_mplus);

  int k=0;
  for (int i=0; i < n_bins_reco_nd_nc; i++)
    for (int j=0; j < n_bins_true_nd; j++)
    {
      sm_nd_nc_truenc_minos[k]   = NDNC_TrueNC_minos->GetBinContent(j+1, i+1);
      sm_nd_nc_numu_minos[k]     = NDNC_NuMu_minos->GetBinContent(j+1, i+1);
      sm_nd_nc_beamnue_minos[k]  = NDNC_BeamNue_minos->GetBinContent(j+1, i+1);
      sm_nd_nc_appnue_minos[k]   = NDNC_AppNue_minos->GetBinContent(j+1, i+1);
      sm_nd_nc_appnutau_minos[k] = NDNC_AppNuTau_minos->GetBinContent(j+1, i+1);

      sm_nd_nc_truenc_mplus[k]   = NDNC_TrueNC_mplus->GetBinContent(j+1, i+1);
      sm_nd_nc_numu_mplus[k]     = NDNC_NuMu_mplus->GetBinContent(j+1, i+1);
      sm_nd_nc_beamnue_mplus[k]  = NDNC_BeamNue_mplus->GetBinContent(j+1, i+1);
      sm_nd_nc_appnue_mplus[k]   = NDNC_AppNue_mplus->GetBinContent(j+1, i+1);
      sm_nd_nc_appnutau_mplus[k] = NDNC_AppNuTau_mplus->GetBinContent(j+1, i+1);
      k++;
    }
  k = 0;
  for (int i=0; i < n_bins_reco_fd_nc; i++)
    for (int j=0; j < n_bins_true_fd; j++)
    {
      sm_fd_nc_truenc_minos[k]   = FDNC_TrueNC_minos->GetBinContent(j+1, i+1);
      sm_fd_nc_numu_minos[k]     = FDNC_NuMu_minos->GetBinContent(j+1, i+1);
      sm_fd_nc_beamnue_minos[k]  = FDNC_BeamNue_minos->GetBinContent(j+1, i+1);
      sm_fd_nc_appnue_minos[k]   = FDNC_AppNue_minos->GetBinContent(j+1, i+1);
      sm_fd_nc_appnutau_minos[k] = FDNC_AppNuTau_minos->GetBinContent(j+1, i+1);

      sm_fd_nc_truenc_mplus[k]   = FDNC_TrueNC_mplus->GetBinContent(j+1, i+1);
      sm_fd_nc_numu_mplus[k]     = FDNC_NuMu_mplus->GetBinContent(j+1, i+1);
      sm_fd_nc_beamnue_mplus[k]  = FDNC_BeamNue_mplus->GetBinContent(j+1, i+1);
      sm_fd_nc_appnue_mplus[k]   = FDNC_AppNue_mplus->GetBinContent(j+1, i+1);
      sm_fd_nc_appnutau_mplus[k] = FDNC_AppNuTau_mplus->GetBinContent(j+1, i+1);
      k++;
    }
  k = 0;
  for (int i=0; i < n_bins_reco_nd_cc; i++)
    for (int j=0; j < n_bins_true_nd; j++)
    {
      sm_nd_cc_truenc_minos[k]   = NDCC_TrueNC_minos->GetBinContent(j+1, i+1);
      sm_nd_cc_numu_minos[k]     = NDCC_NuMu_minos->GetBinContent(j+1, i+1);
      sm_nd_cc_beamnue_minos[k]  = NDCC_BeamNue_minos->GetBinContent(j+1, i+1);
      sm_nd_cc_appnue_minos[k]   = NDCC_AppNue_minos->GetBinContent(j+1, i+1);
      sm_nd_cc_appnutau_minos[k] = NDCC_AppNuTau_minos->GetBinContent(j+1, i+1);

      sm_nd_cc_truenc_mplus[k]   = NDCC_TrueNC_mplus->GetBinContent(j+1, i+1);
      sm_nd_cc_numu_mplus[k]     = NDCC_NuMu_mplus->GetBinContent(j+1, i+1);
      sm_nd_cc_beamnue_mplus[k]  = NDCC_BeamNue_mplus->GetBinContent(j+1, i+1);
      sm_nd_cc_appnue_mplus[k]   = NDCC_AppNue_mplus->GetBinContent(j+1, i+1);
      sm_nd_cc_appnutau_mplus[k] = NDCC_AppNuTau_mplus->GetBinContent(j+1, i+1);
      k++;
    }
  k = 0;
  for (int i=0; i < n_bins_reco_fd_cc; i++)
    for (int j=0; j < n_bins_true_fd; j++)
    {
      sm_fd_cc_truenc_minos[k]   = FDCC_TrueNC_minos->GetBinContent(j+1, i+1);
      sm_fd_cc_numu_minos[k]     = FDCC_NuMu_minos->GetBinContent(j+1, i+1);
      sm_fd_cc_beamnue_minos[k]  = FDCC_BeamNue_minos->GetBinContent(j+1, i+1);
      sm_fd_cc_appnue_minos[k]   = FDCC_AppNue_minos->GetBinContent(j+1, i+1);
      sm_fd_cc_appnutau_minos[k] = FDCC_AppNuTau_minos->GetBinContent(j+1, i+1);

      sm_fd_cc_truenc_mplus[k]   = FDCC_TrueNC_mplus->GetBinContent(j+1, i+1);
      sm_fd_cc_numu_mplus[k]     = FDCC_NuMu_mplus->GetBinContent(j+1, i+1);
      sm_fd_cc_beamnue_mplus[k]  = FDCC_BeamNue_mplus->GetBinContent(j+1, i+1);
      sm_fd_cc_appnue_mplus[k]   = FDCC_AppNue_mplus->GetBinContent(j+1, i+1);
      sm_fd_cc_appnutau_mplus[k] = FDCC_AppNuTau_mplus->GetBinContent(j+1, i+1);
      k++;
    }

  // Covariance matrices
  cov_matrix_nc_rel = CoVarNC_relative->GetMatrixArray();
  cov_matrix_cc_rel = CoVarCC_relative->GetMatrixArray();

  cov_matrix_nc_abs = gsl_matrix_calloc(n_bins_reco_fd_nc+n_bins_reco_nd_nc,
                                        n_bins_reco_fd_nc+n_bins_reco_nd_nc);
  cov_matrix_cc_abs = gsl_matrix_calloc(n_bins_reco_fd_cc+n_bins_reco_nd_cc,
                                        n_bins_reco_fd_cc+n_bins_reco_nd_cc);
  cov_matrix_perm_nc = gsl_permutation_calloc(n_bins_reco_fd_nc+n_bins_reco_nd_nc);
  cov_matrix_perm_cc = gsl_permutation_calloc(n_bins_reco_fd_cc+n_bins_reco_nd_cc);

  // Data
  data_nd_nc = new double[n_bins_reco_nd_nc];  assert(data_nd_nc);
  data_fd_nc = new double[n_bins_reco_fd_nc];  assert(data_fd_nc);
  data_nd_cc = new double[n_bins_reco_nd_cc];  assert(data_nd_cc);
  data_fd_cc = new double[n_bins_reco_fd_cc];  assert(data_fd_cc);
  for (int i=0; i < n_bins_reco_nd_nc; i++)
    data_nd_nc[i] = ND_dataNC_minos->GetBinContent(i+1) + ND_dataNC_mplus->GetBinContent(i+1);
  for (int i=0; i < n_bins_reco_fd_nc; i++)
    data_fd_nc[i] = FD_dataNC_minos->GetBinContent(i+1) + FD_dataNC_mplus->GetBinContent(i+1);
  for (int i=0; i < n_bins_reco_nd_cc; i++)
    data_nd_cc[i] = ND_dataCC_minos->GetBinContent(i+1) + ND_dataCC_mplus->GetBinContent(i+1);
  for (int i=0; i < n_bins_reco_fd_cc; i++)
    data_fd_cc[i] = FD_dataCC_minos->GetBinContent(i+1) + FD_dataCC_mplus->GetBinContent(i+1);

  // Weights for systematic biases
  wgt_miscal_nd_cc_minos  = CC_IMiscal_ND_minos->GetArray() + 1;
  wgt_dist_nd_cc_minos    = CC_IDist_ND_minos->GetArray() + 1;
  wgt_miscal_fd_cc_minos  = CC_IMiscal_FD_minos->GetArray() + 1;
  wgt_dist_fd_cc_minos    = CC_IDist_FD_minos->GetArray() + 1;
  wgt_miscal_nd_nc_minos  = NC_IMiscal_ND_minos->GetArray() + 1;
  wgt_dist_nd_nc_minos    = NC_IDist_ND_minos->GetArray() + 1;
  wgt_miscal_fd_nc_minos  = NC_IMiscal_FD_minos->GetArray() + 1;
  wgt_dist_fd_nc_minos    = NC_IDist_FD_minos->GetArray() + 1;
  wgt_miscal_nd_cc_mplus  = CC_IMiscal_ND_mplus->GetArray() + 1;
  wgt_dist_nd_cc_mplus    = CC_IDist_ND_mplus->GetArray() + 1;
  wgt_miscal_fd_cc_mplus  = CC_IMiscal_FD_mplus->GetArray() + 1;
  wgt_dist_fd_cc_mplus    = CC_IDist_FD_mplus->GetArray() + 1;
  wgt_miscal_nd_nc_mplus  = NC_IMiscal_ND_mplus->GetArray() + 1;
  wgt_dist_nd_nc_mplus    = NC_IDist_ND_mplus->GetArray() + 1;
  wgt_miscal_fd_nc_mplus  = NC_IMiscal_FD_mplus->GetArray() + 1;
  wgt_dist_fd_nc_mplus    = NC_IDist_FD_mplus->GetArray() + 1;

  return 0;
}


// -------------------------------------------------------------------------
double MINOS_2017_prior_callback(const gsl_vector *x, void *params)
// -------------------------------------------------------------------------
// Callback function for GSL minimizer which we use to fit the nuisance
// parameters.  Returns the chi^2 value for the nuisance parameters given
// in x. The nuisance parameters are:
//   x[0]: horn current miscalibration MINOS
//   x[1]: horn current distribution MINOS
//   x[2]: horn current miscalibration MINOS+
//   x[3]: horn current distribution MINOS+
// -------------------------------------------------------------------------
{
  double x_miscal_minos = x->data[0];
  double x_dist_minos   = x->data[1];
  double x_miscal_mplus = x->data[2];
  double x_dist_mplus   = x->data[3];

  // Number of bins (we exploit here that the binning is the same for MINOS and MINOS+
  int n_bins_reco_nd_cc = ND_dataCC_minos->GetXaxis()->GetNbins();
  int n_bins_reco_nd_nc = ND_dataNC_minos->GetXaxis()->GetNbins();
  int n_bins_reco_fd_cc = FD_dataCC_minos->GetXaxis()->GetNbins();
  int n_bins_reco_fd_nc = FD_dataNC_minos->GetXaxis()->GetNbins();

  // Shift event rates by systematic biases, arrange FD+ND rates into single array
#ifdef MINOS_2017_FREE_NORM
  for (int i=0; i < n_bins_reco_fd_nc; i++) //FIXME remove
    pred_nc->data[i] = (1.+x->data[4]) * ( pred_fd_nc_minos[i] 
                         * (1.0 + x_miscal_minos * (wgt_miscal_fd_nc_minos[i] - 1.0))
                         * (1.0 + x_dist_minos * (wgt_dist_fd_nc_minos[i] - 1.0))
                     + pred_fd_nc_mplus[i]
                         * (1.0 + x_miscal_mplus * (wgt_miscal_fd_nc_mplus[i] - 1.0))
                         * (1.0 + x_dist_mplus * (wgt_dist_fd_nc_mplus[i] - 1.0)) );
  for (int i=0; i < n_bins_reco_nd_nc; i++)
    pred_nc->data[i+n_bins_reco_fd_nc] = (1.+x->data[4]) * ( pred_nd_nc_minos[i]
                         * (1.0 + x_miscal_minos * (wgt_miscal_nd_nc_minos[i] - 1.0))
                         * (1.0 + x_dist_minos * (wgt_dist_nd_nc_minos[i] - 1.0))
                     + pred_nd_nc_mplus[i]
                         * (1.0 + x_miscal_mplus * (wgt_miscal_nd_nc_mplus[i] - 1.0))
                         * (1.0 + x_dist_mplus * (wgt_dist_nd_nc_mplus[i] - 1.0)) );
  for (int i=0; i < n_bins_reco_fd_cc; i++)
    pred_cc->data[i] = (1.+x->data[4]) * ( pred_fd_cc_minos[i] 
                         * (1.0 + x_miscal_minos * (wgt_miscal_fd_cc_minos[i] - 1.0))
                         * (1.0 + x_dist_minos * (wgt_dist_fd_cc_minos[i] - 1.0))
                     + pred_fd_cc_mplus[i]
                         * (1.0 + x_miscal_mplus * (wgt_miscal_fd_cc_mplus[i] - 1.0))
                         * (1.0 + x_dist_mplus * (wgt_dist_fd_cc_mplus[i] - 1.0)) );
  for (int i=0; i < n_bins_reco_nd_cc; i++)
    pred_cc->data[i+n_bins_reco_fd_cc] = (1.+x->data[4]) * ( pred_nd_cc_minos[i]
                         * (1.0 + x_miscal_minos * (wgt_miscal_nd_cc_minos[i] - 1.0))
                         * (1.0 + x_dist_minos * (wgt_dist_nd_cc_minos[i] - 1.0))
                     + pred_nd_cc_mplus[i]
                         * (1.0 + x_miscal_mplus * (wgt_miscal_nd_cc_mplus[i] - 1.0))
                         * (1.0 + x_dist_mplus * (wgt_dist_nd_cc_mplus[i] - 1.0)) );
#else
  for (int i=0; i < n_bins_reco_fd_nc; i++)
    pred_nc->data[i] = pred_fd_nc_minos[i] 
                         * (1.0 + x_miscal_minos * (wgt_miscal_fd_nc_minos[i] - 1.0))
                         * (1.0 + x_dist_minos * (wgt_dist_fd_nc_minos[i] - 1.0))
                     + pred_fd_nc_mplus[i]
                         * (1.0 + x_miscal_mplus * (wgt_miscal_fd_nc_mplus[i] - 1.0))
                         * (1.0 + x_dist_mplus * (wgt_dist_fd_nc_mplus[i] - 1.0));
  for (int i=0; i < n_bins_reco_nd_nc; i++)
    pred_nc->data[i+n_bins_reco_fd_nc] = pred_nd_nc_minos[i]
                         * (1.0 + x_miscal_minos * (wgt_miscal_nd_nc_minos[i] - 1.0))
                         * (1.0 + x_dist_minos * (wgt_dist_nd_nc_minos[i] - 1.0))
                     + pred_nd_nc_mplus[i]
                         * (1.0 + x_miscal_mplus * (wgt_miscal_nd_nc_mplus[i] - 1.0))
                         * (1.0 + x_dist_mplus * (wgt_dist_nd_nc_mplus[i] - 1.0));
  for (int i=0; i < n_bins_reco_fd_cc; i++)
    pred_cc->data[i] = pred_fd_cc_minos[i] 
                         * (1.0 + x_miscal_minos * (wgt_miscal_fd_cc_minos[i] - 1.0))
                         * (1.0 + x_dist_minos * (wgt_dist_fd_cc_minos[i] - 1.0))
                     + pred_fd_cc_mplus[i]
                         * (1.0 + x_miscal_mplus * (wgt_miscal_fd_cc_mplus[i] - 1.0))
                         * (1.0 + x_dist_mplus * (wgt_dist_fd_cc_mplus[i] - 1.0));
  for (int i=0; i < n_bins_reco_nd_cc; i++)
    pred_cc->data[i+n_bins_reco_fd_cc] = pred_nd_cc_minos[i]
                         * (1.0 + x_miscal_minos * (wgt_miscal_nd_cc_minos[i] - 1.0))
                         * (1.0 + x_dist_minos * (wgt_dist_nd_cc_minos[i] - 1.0))
                     + pred_nd_cc_mplus[i]
                         * (1.0 + x_miscal_mplus * (wgt_miscal_nd_cc_mplus[i] - 1.0))
                         * (1.0 + x_dist_mplus * (wgt_dist_nd_cc_mplus[i] - 1.0));
#endif // MINOS_2017_FREE_NORM?

  // Rescale covariance matrices
  int k=0;
  for (int i=0; i < n_bins_reco_fd_nc+n_bins_reco_nd_nc; i++)
    for (int j=0; j < n_bins_reco_fd_nc+n_bins_reco_nd_nc; j++)
    {
      double stat = ( (i==j) ? pred_nc->data[i] : 0.0 );
      double syst = cov_matrix_nc_rel[k] * pred_nc->data[i] * pred_nc->data[j];
      cov_matrix_nc_abs->data[k] = stat + syst;
      k++;
    }
  k=0;
  for (int i=0; i < n_bins_reco_fd_cc+n_bins_reco_nd_cc; i++)
    for (int j=0; j < n_bins_reco_fd_cc+n_bins_reco_nd_cc; j++)
    {
      double stat = ( (i==j) ? pred_cc->data[i] : 0.0 );
      double syst = cov_matrix_cc_rel[k] * pred_cc->data[i] * pred_cc->data[j];
      cov_matrix_cc_abs->data[k] = stat + syst;
      k++;
    }

  // Form differences of data and prediction
  for (int i=0; i<n_bins_reco_fd_nc; i++) pred_nc->data[i] -= data_fd_nc[i];
  for (int i=0; i<n_bins_reco_nd_nc; i++) pred_nc->data[i+n_bins_reco_fd_nc] -= data_nd_nc[i];
  for (int i=0; i<n_bins_reco_fd_cc; i++) pred_cc->data[i] -= data_fd_cc[i];
  for (int i=0; i<n_bins_reco_nd_cc; i++) pred_cc->data[i+n_bins_reco_fd_cc] -= data_nd_cc[i];

  // Compute chi^2
  double chi2_nc=0.0, chi2_cc=0.0;
  int sgn_nc, sgn_cc;
  gsl_linalg_LU_decomp(cov_matrix_nc_abs, cov_matrix_perm_nc, &sgn_nc); // Solve M.y = x
  gsl_linalg_LU_solve(cov_matrix_nc_abs, cov_matrix_perm_nc, pred_nc, temp_nc);
  gsl_blas_ddot(pred_nc, temp_nc, &chi2_nc); // Multiply y.x = x.M^{-1}.x

  gsl_linalg_LU_decomp(cov_matrix_cc_abs, cov_matrix_perm_cc, &sgn_cc);
  gsl_linalg_LU_solve(cov_matrix_cc_abs, cov_matrix_perm_cc, pred_cc, temp_cc);
  gsl_blas_ddot(pred_cc, temp_cc, &chi2_cc);

  // Compute penalty (pull) terms for beam focusing parameters
  double chi2_pull = SQR(x_miscal_minos) + SQR(x_miscal_mplus)
                   + SQR(x_dist_minos) + SQR(x_dist_mplus);

  // Return total chi^2
  return chi2_nc + chi2_cc + chi2_pull;
}


// -------------------------------------------------------------------------
double MINOS_2017_prior(const glb_params params, unsigned flags)
// -------------------------------------------------------------------------
// chi^2 for the MINOS/MINOS+ 2017 fit, based on the code provided
// in the data release (https://arxiv.org/abs/1710.06488)
// Bypasses GLoBES except for the computation of the oscillation
// probabilities.
// If flags contains MINOS_2017_PRINT_RATES, the event rates are printed.
// -------------------------------------------------------------------------
{
  const double L_nd   = 1.04;     // km
  const double L_fd   = 735.0;    // km
  const double rho_nd = 0.0;      // grams/cm^3
  const double rho_fd = 2.8;      // grams/cm^3
  int k;

//  const double sig_miscal_minos = 5.04310835513857469e-01;
//  const double sig_dist_minos = 4.99672999073187185e-01;
//  const double sig_miscal_mplus = 1.43020989799832421e+00;
//  const double sig_dist_mplus = -4.71471682040894180e-01;

  // Number of bins (we exploit here that the binning is the same for MINOS and MINOS+,
  // and the true energy binning is moreover the same for CC and NC events)
  int n_bins_true_nd = NDCC_TrueNC_minos->GetXaxis()->GetNbins();
  int n_bins_true_fd = FDCC_TrueNC_minos->GetXaxis()->GetNbins();
  int n_bins_reco_nd_cc = ND_dataCC_minos->GetXaxis()->GetNbins();
  int n_bins_reco_nd_nc = ND_dataNC_minos->GetXaxis()->GetNbins();
  int n_bins_reco_fd_cc = FD_dataCC_minos->GetXaxis()->GetNbins();
  int n_bins_reco_fd_nc = FD_dataNC_minos->GetXaxis()->GetNbins();

  double P_nd_truenc[n_bins_true_nd], P_nd_numu[n_bins_true_nd], P_nd_beamnue[n_bins_true_nd];
  double P_nd_appnue[n_bins_true_nd], P_nd_appnutau[n_bins_true_nd];
  double P_fd_truenc[n_bins_true_fd], P_fd_numu[n_bins_true_fd], P_fd_beamnue[n_bins_true_fd];
  double P_fd_appnue[n_bins_true_fd], P_fd_appnutau[n_bins_true_fd];

  // Prepare oscillation engine
  snu_set_oscillation_parameters(params, NULL);

#if defined NU_USE_NUSQUIDS
  // Compute probabilities in nuSQuIDS
  // ---------------------------------
  //
  // Note:
  // - we need to build up the initial state vectors backwards as NuSQuIDS
  //   expects energies to be ascending
  // - To avoid having to run NuSQuIDS separately for each flux component,
  //   we compute approximate "oscillation+decay+regeneration" probabilities.
  //   In doing so we neglect for instance the disappearance of beam \nu_e
  // - Some initial state bins may be empty, but may contain events after neutrino
  //   decay. To avoid inifinite "oscillation probabilities", we set all initial
  //   state bins to a minimum value
  double P_fd[n_bins_true_fd][2][3];
  double bin_centers_fd[n_bins_true_fd];
  double ini_state_nu_fd[n_bins_true_fd][3];
  double ini_state_nubar_fd[n_bins_true_fd][3];
  double ini_state_total_fd = 0.0;
  memset(ini_state_nu_fd,    0, n_bins_true_fd * 3 * sizeof(ini_state_nu_fd[0][0]));
  memset(ini_state_nubar_fd, 0, n_bins_true_fd * 3 * sizeof(ini_state_nubar_fd[0][0]));
  k = 0;
  for (int i=0; i < n_bins_reco_fd_nc; i++)
    for (int j=0; j < n_bins_true_fd; j++)
    {
      int m = n_bins_true_fd - j - 1;
      double flux_nu_e  = sm_fd_cc_beamnue_minos[k] + sm_fd_cc_beamnue_mplus[k];
      double flux_nu_mu = sm_fd_cc_numu_minos[k] + sm_fd_cc_numu_mplus[k];
      ini_state_nu_fd[m][NU_E]  += flux_nu_e;
      ini_state_nu_fd[m][NU_MU] += flux_nu_mu;
      ini_state_total_fd        += flux_nu_e + flux_nu_mu;
      k++;
    }
  double eps_fd = 1e-10 * ini_state_total_fd;
  for (int j=0; j < n_bins_true_fd; j++)
  {
    for (int f=0; f < 3; f++)
      if (ini_state_nu_fd[j][f] < eps_fd)
        ini_state_nu_fd[j][f] = eps_fd;
    bin_centers_fd[j] = L_fd / FDCC_TrueNC_minos->GetXaxis()->GetBinCenter(n_bins_true_fd-j);
  }
  nu_hook_probability_matrix_nusquids(P_fd, n_bins_true_fd, bin_centers_fd,
                                      ini_state_nu_fd, ini_state_nubar_fd, 1, &L_fd, &rho_fd);
  for (int j=0; j < n_bins_true_fd; j++)
  {
    int m = n_bins_true_fd - j - 1;
    double ini_state_nu_e  = ini_state_nu_fd[m][NU_E]  + ini_state_nubar_fd[m][NU_E];
    double ini_state_nu_mu = ini_state_nu_fd[m][NU_MU] + ini_state_nubar_fd[m][NU_MU];
    if (ini_state_nu_mu > 1e-100)
    {
      P_fd_truenc[j]   = ( P_fd[m][0][NU_E] + P_fd[m][0][NU_MU] + P_fd[m][0][NU_TAU]
                         + P_fd[m][1][NU_E] + P_fd[m][1][NU_MU] + P_fd[m][1][NU_TAU] )
                         / (ini_state_nu_e + ini_state_nu_mu);
      P_fd_numu[j]     = ( P_fd[m][0][NU_MU] + P_fd[m][1][NU_MU] ) / ini_state_nu_mu;
      P_fd_appnue[j]   = ( P_fd[m][0][NU_E] + P_fd[m][1][NU_E] ) / ini_state_nu_mu;
      P_fd_appnutau[j] = ( P_fd[m][0][NU_TAU] + P_fd[m][1][NU_TAU] ) / ini_state_nu_mu;
    }
    else
      P_fd_truenc[j] = P_fd_numu[j] = P_fd_appnue[j] = P_fd_appnutau[j] = 0.0;

    P_fd_beamnue[j] = 1.0;
  }

  double P_nd[n_bins_true_nd][2][3];
  double bin_centers_nd[n_bins_true_nd];
  double ini_state_nu_nd[n_bins_true_nd][3];
  double ini_state_nubar_nd[n_bins_true_nd][3];
  double ini_state_total_nd = 0.0;
  memset(ini_state_nu_nd,    0, n_bins_true_nd * 3 * sizeof(ini_state_nu_nd[0][0]));
  memset(ini_state_nubar_nd, 0, n_bins_true_nd * 3 * sizeof(ini_state_nubar_nd[0][0]));
  k = 0;
  for (int i=0; i < n_bins_reco_nd_nc; i++)
    for (int j=0; j < n_bins_true_nd; j++)
    {
      int m = n_bins_true_nd - j - 1;
      double flux_nu_e  = sm_nd_cc_beamnue_minos[k] + sm_nd_cc_beamnue_mplus[k];
      double flux_nu_mu = sm_nd_cc_numu_minos[k] + sm_nd_cc_numu_mplus[k];
      ini_state_nu_nd[m][NU_E]  += flux_nu_e;
      ini_state_nu_nd[m][NU_MU] += flux_nu_mu;
      ini_state_total_nd        += flux_nu_e + flux_nu_mu;
      k++;
    }
  double eps_nd = 1e-10 * ini_state_total_nd;
  for (int j=0; j < n_bins_true_nd; j++)
  {
    for (int f=0; f < 3; f++)
      if (ini_state_nu_nd[j][f] < eps_nd)
        ini_state_nu_nd[j][f] = eps_nd;
    bin_centers_nd[j] = L_nd / NDCC_TrueNC_minos->GetXaxis()->GetBinCenter(n_bins_true_nd-j);
  }
  nu_hook_probability_matrix_nusquids(P_nd, n_bins_true_nd, bin_centers_nd,
                                      ini_state_nu_nd, ini_state_nubar_nd, 1, &L_nd, &rho_nd);
  for (int j=0; j < n_bins_true_nd; j++)
  {
    int m = n_bins_true_nd - j - 1;
    double ini_state_nu_e  = ini_state_nu_nd[m][NU_E]  + ini_state_nubar_nd[m][NU_E];
    double ini_state_nu_mu = ini_state_nu_nd[m][NU_MU] + ini_state_nubar_nd[m][NU_MU];
    if (ini_state_nu_mu > 1e-100)
    {
      P_nd_truenc[j]   = ( P_nd[m][0][NU_E] + P_nd[m][0][NU_MU] + P_nd[m][0][NU_TAU]
                         + P_nd[m][1][NU_E] + P_nd[m][1][NU_MU] + P_nd[m][1][NU_TAU] )
                         / (ini_state_nu_e + ini_state_nu_mu);
      P_nd_numu[j]     = ( P_nd[m][0][NU_MU] + P_nd[m][1][NU_MU] ) / ini_state_nu_mu;
      P_nd_appnue[j]   = ( P_nd[m][0][NU_E] + P_nd[m][1][NU_E] ) / ini_state_nu_mu;
      P_nd_appnutau[j] = ( P_nd[m][0][NU_TAU] + P_nd[m][1][NU_TAU] ) / ini_state_nu_mu;
    }
    else
      P_nd_truenc[j] = P_nd_numu[j] = P_nd_appnue[j] = P_nd_appnutau[j] = 0.0;

    P_nd_beamnue[j] = 1.0;
  }

#elif defined MINOS_2017_USE_GLOBES_LOWPASS

  // Compute probabilities in GLoBES using GLoBES lowpass filter
  // -----------------------------------------------------------
  for (int j=0; j < n_bins_true_nd; j++)
  {
    double P_nd[3][3];
    snu_probability_matrix(P_nd, +1, L_nd / NDCC_TrueNC_minos->GetXaxis()->GetBinCenter(j+1),
                           1, &L_nd, &rho_nd, 0.3, NULL);
    P_nd_truenc[j]   = P_nd[NU_MU][NU_E] + P_nd[NU_MU][NU_MU] + P_nd[NU_MU][NU_TAU];
    P_nd_numu[j]     = P_nd[NU_MU][NU_MU];
    P_nd_beamnue[j]  = P_nd[NU_E][NU_E];
    P_nd_appnue[j]   = P_nd[NU_MU][NU_E];
    P_nd_appnutau[j] = P_nd[NU_MU][NU_TAU];
  }
  for (int j=0; j < n_bins_true_fd; j++)
  {
    double P_fd[3][3];
    snu_probability_matrix(P_fd, +1, L_fd / FDCC_TrueNC_minos->GetXaxis()->GetBinCenter(j+1),
                           1, &L_fd, &rho_fd, 0.3, NULL);
    P_fd_truenc[j]   = P_fd[NU_MU][NU_E] + P_fd[NU_MU][NU_MU] + P_fd[NU_MU][NU_TAU];
    P_fd_numu[j]     = P_fd[NU_MU][NU_MU];
    P_fd_beamnue[j]  = P_fd[NU_E][NU_E];
    P_fd_appnue[j]   = P_fd[NU_MU][NU_E];
    P_fd_appnutau[j] = P_fd[NU_MU][NU_TAU];
  }

#else

  // Compute probabilities in GLoBES using MINOS' method to deal with aliasing
  // -------------------------------------------------------------------------
  const double dmsq43 = glbGetOscParamByName(params, "DM41")
                          - glbGetOscParamByName(params, "DM31");
  for (int j=0; j < n_bins_true_nd; j++)
  {
    // Averaging of oscillation probabilities follows MINOS code
    // FIXME this is valid only in 3+1 models
    double P_nd[2][3][3];
    double LoverE[2];
    double bc = NDCC_TrueNC_minos->GetXaxis()->GetBinCenter(j+1);
    double w = NDCC_TrueNC_minos->GetXaxis()->GetBinWidth(j+1);
    double arg = dmsq43 * w * 0.25 * KM/GEV;
    double sample = 0.5/M_SQRT_3 * w;
    if (arg != 0)
      sample = acos(sin(arg)/arg) / arg * 0.5*w;
    LoverE[0] = bc - sample;
    LoverE[1] = bc + sample;
    for (int k=0; k < 2; k++)
      snu_probability_matrix(P_nd[k], +1, L_nd / LoverE[k], 1, &L_nd, &rho_nd, -1.0, NULL);
    for (int k=0; k < 3; k++)
      for (int l=0; l < 3; l++)
        P_nd[0][k][l] = 0.5 * (P_nd[0][k][l] + P_nd[1][k][l]);

    P_nd_truenc[j]   = P_nd[0][NU_MU][NU_E] + P_nd[0][NU_MU][NU_MU] + P_nd[0][NU_MU][NU_TAU];
    P_nd_numu[j]     = P_nd[0][NU_MU][NU_MU];
    P_nd_beamnue[j]  = P_nd[0][NU_E][NU_E];
    P_nd_appnue[j]   = P_nd[0][NU_MU][NU_E];
    P_nd_appnutau[j] = P_nd[0][NU_MU][NU_TAU];
  }

  for (int j=0; j < n_bins_true_fd; j++)
  {
    // Averaging of oscillation probabilities follows MINOS code
    // FIXME this is valid only in 3+1 models
    double P_fd[2][3][3];
    double LoverE[2];
    double bc = FDCC_TrueNC_minos->GetXaxis()->GetBinCenter(j+1);
    double w = FDCC_TrueNC_minos->GetXaxis()->GetBinWidth(j+1);
    double arg = dmsq43 * w * 0.25 * KM/GEV;
    double sample = 0.5/M_SQRT_3 * w;
    if (arg != 0)
      sample = acos(sin(arg)/arg) / arg * 0.5*w;
    LoverE[0] = bc - sample;
    LoverE[1] = bc + sample;
    for (int k=0; k < 2; k++)
      snu_probability_matrix(P_fd[k], +1, L_fd / LoverE[k], 1, &L_fd, &rho_fd, -1.0, NULL);
    for (int k=0; k < 3; k++)
      for (int l=0; l < 3; l++)
        P_fd[0][k][l] = 0.5 * (P_fd[0][k][l] + P_fd[1][k][l]);

    P_fd_truenc[j]   = P_fd[0][NU_MU][NU_E] + P_fd[0][NU_MU][NU_MU] + P_fd[0][NU_MU][NU_TAU];
    P_fd_numu[j]     = P_fd[0][NU_MU][NU_MU];
    P_fd_beamnue[j]  = P_fd[0][NU_E][NU_E];
    P_fd_appnue[j]   = P_fd[0][NU_MU][NU_E];
    P_fd_appnutau[j] = P_fd[0][NU_MU][NU_TAU];
  }
#endif // MINOS_2017_USE_GLOBES_LOWPASS

  // Compute predicted event rates
  k = 0;
  for (int i=0; i < n_bins_reco_nd_nc; i++)
  {
    pred_nd_nc_minos[i] = 0.0;
    pred_nd_nc_mplus[i] = 0.0;
    for (int j=0; j < n_bins_true_nd; j++)
    {
//      if (i==0) //FIXME
//        printf("#### P-ND #### %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g\n", 
//            L_nd / NDNC_TrueNC_minos->GetXaxis()->GetBinCenter(n_bins_true_nd-j),
//            P_nd_numu[j], P_nd_truenc[j], P_nd_beamnue[j], P_nd_appnue[j],
//            P_nd_appnutau[j]);//FIXME
      pred_nd_nc_minos[i] += sm_nd_nc_truenc_minos[k] * P_nd_truenc[j]
                           + sm_nd_nc_numu_minos[k] * P_nd_numu[j]
                           + sm_nd_nc_beamnue_minos[k] * P_nd_beamnue[j]
                           + sm_nd_nc_appnue_minos[k] * P_nd_appnue[j]
                           + sm_nd_nc_appnutau_minos[k] * P_nd_appnutau[j];
      pred_nd_nc_mplus[i] += sm_nd_nc_truenc_mplus[k] * P_nd_truenc[j]
                           + sm_nd_nc_numu_mplus[k] * P_nd_numu[j]
                           + sm_nd_nc_beamnue_mplus[k] * P_nd_beamnue[j]
                           + sm_nd_nc_appnue_mplus[k] * P_nd_appnue[j]
                           + sm_nd_nc_appnutau_mplus[k] * P_nd_appnutau[j];
////      if (j==n_bins_true_nd-1) //FIXME
//      if (i==0) // FIXME
//        printf("#### pred ND NC #### %d   %10.5g %10.5g\n", j,
//            L_nd / NDNC_TrueNC_minos->GetXaxis()->GetBinCenter(n_bins_true_nd-j),
//            pred_nd_nc_minos[i], pred_nd_nc_mplus[i]);//FIXME
      k++;
    }
  }
  k = 0;
  for (int i=0; i < n_bins_reco_fd_nc; i++)
  {
    pred_fd_nc_minos[i] = 0.0;
    pred_fd_nc_mplus[i] = 0.0;
    for (int j=0; j < n_bins_true_fd; j++)
    {
      pred_fd_nc_minos[i] += sm_fd_nc_truenc_minos[k] * P_fd_truenc[j]
                           + sm_fd_nc_numu_minos[k] * P_fd_numu[j]
                           + sm_fd_nc_beamnue_minos[k] * P_fd_beamnue[j]
                           + sm_fd_nc_appnue_minos[k] * P_fd_appnue[j]
                           + sm_fd_nc_appnutau_minos[k] * P_fd_appnutau[j];
      pred_fd_nc_mplus[i] += sm_fd_nc_truenc_mplus[k] * P_fd_truenc[j]
                           + sm_fd_nc_numu_mplus[k] * P_fd_numu[j]
                           + sm_fd_nc_beamnue_mplus[k] * P_fd_beamnue[j]
                           + sm_fd_nc_appnue_mplus[k] * P_fd_appnue[j]
                           + sm_fd_nc_appnutau_mplus[k] * P_fd_appnutau[j];
      k++;
    }
  }
  k = 0;
  for (int i=0; i < n_bins_reco_nd_cc; i++)
  {
    pred_nd_cc_minos[i] = 0.0;
    pred_nd_cc_mplus[i] = 0.0;
    for (int j=0; j < n_bins_true_nd; j++)
    {
      pred_nd_cc_minos[i] += sm_nd_cc_truenc_minos[k] * P_nd_truenc[j]
                           + sm_nd_cc_numu_minos[k] * P_nd_numu[j]
                           + sm_nd_cc_beamnue_minos[k] * P_nd_beamnue[j]
                           + sm_nd_cc_appnue_minos[k] * P_nd_appnue[j]
                           + sm_nd_cc_appnutau_minos[k] * P_nd_appnutau[j];
      pred_nd_cc_mplus[i] += sm_nd_cc_truenc_mplus[k] * P_nd_truenc[j]
                           + sm_nd_cc_numu_mplus[k] * P_nd_numu[j]
                           + sm_nd_cc_beamnue_mplus[k] * P_nd_beamnue[j]
                           + sm_nd_cc_appnue_mplus[k] * P_nd_appnue[j]
                           + sm_nd_cc_appnutau_mplus[k] * P_nd_appnutau[j];
      k++;
    }
  }
  k = 0;
  for (int i=0; i < n_bins_reco_fd_cc; i++)
  {
    pred_fd_cc_minos[i] = 0.0;
    pred_fd_cc_mplus[i] = 0.0;
    for (int j=0; j < n_bins_true_fd; j++)
    {
//      if (i==0) //FIXME
//        printf("#### P-FD #### %10.5g %10.5g\n", 
//            L_fd / FDCC_TrueNC_minos->GetXaxis()->GetBinCenter(n_bins_true_fd-j),
//            P_fd_numu[j]);//FIXME
      pred_fd_cc_minos[i] += sm_fd_cc_truenc_minos[k] * P_fd_truenc[j]
                           + sm_fd_cc_numu_minos[k] * P_fd_numu[j]
                           + sm_fd_cc_beamnue_minos[k] * P_fd_beamnue[j]
                           + sm_fd_cc_appnue_minos[k] * P_fd_appnue[j]
                           + sm_fd_cc_appnutau_minos[k] * P_fd_appnutau[j];
      pred_fd_cc_mplus[i] += sm_fd_cc_truenc_mplus[k] * P_fd_truenc[j]
                           + sm_fd_cc_numu_mplus[k] * P_fd_numu[j]
                           + sm_fd_cc_beamnue_mplus[k] * P_fd_beamnue[j]
                           + sm_fd_cc_appnue_mplus[k] * P_fd_appnue[j]
                           + sm_fd_cc_appnutau_mplus[k] * P_fd_appnutau[j];
      k++;
    }
  }


  // Invoke GSL minimizer to fit beam focus parameters
#ifdef MINOS_2017_MINIMIZE
  gsl_multimin_fminimizer *m =
    gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, MINOS_2017_N_NUISANCE);
  assert(m);
  gsl_vector *x = gsl_vector_alloc(MINOS_2017_N_NUISANCE);  assert(x);
  gsl_vector *s = gsl_vector_alloc(MINOS_2017_N_NUISANCE);  assert(s);
  gsl_multimin_function f;
  f.n = MINOS_2017_N_NUISANCE;
  f.f = &MINOS_2017_prior_callback;
  f.params = NULL;
  gsl_vector_set_all(x, 0.0);   // Starting values
  gsl_vector_set_all(s, 0.1);   // Initial step size

  gsl_multimin_fminimizer_set(m, &f, x, s);
  size_t n_iter = 0;
  int status;
  double chi2 = NAN;
  do
  {
    n_iter++;
    if ((status = gsl_multimin_fminimizer_iterate(m)) != GSL_SUCCESS)
      break;
    else
      chi2 = m->fval;

    double size = gsl_multimin_fminimizer_size(m);
    status = gsl_multimin_test_size(size, 1e-2); // Check for convergence
  }
  while (status == GSL_CONTINUE && n_iter < 100);  

  if (x) { gsl_vector_free(x);  x = NULL; }
  if (s) { gsl_vector_free(s);  s = NULL; }
  if (m) { gsl_multimin_fminimizer_free(m);  m = NULL; }
#else
  #warning "Using MINOS/MINOS+ 2017 analysis without minimization over nuisance parameters."
  gsl_vector *x = gsl_vector_calloc(4);  assert(x);
  double chi2 = MINOS_2017_prior_callback(x, NULL);
  if (x)  { gsl_vector_free(x);  x = NULL; }
#endif


  // Print event rates if requested
  if (flags & MINOS_2017_PRINT_RATES)
  {
    int i;

    printf("(* MINOS/MINOS+ 2017: FD NC *)\n");
    printf("{ {\n");
    for (i=0; i < n_bins_reco_fd_nc-1; i++)
      printf("{%12.7g, %12.7g},\n", FD_dataNC_minos->GetXaxis()->GetBinCenter(i+1),
             pred_nc->data[i] + data_fd_nc[i]);
    printf("{%12.7g, %12.7g}", FD_dataNC_minos->GetXaxis()->GetBinCenter(i+1),
             pred_nc->data[i] + data_fd_nc[i]);
    printf("},\n");

    printf("(* MINOS/MINOS+ 2017: ND NC *)\n");
    printf("{\n");
    for (i=0; i < n_bins_reco_nd_nc-1; i++)
      printf("{%12.7g, %12.7g},\n", ND_dataNC_minos->GetXaxis()->GetBinCenter(i+1),
             pred_nc->data[i+n_bins_reco_fd_nc] + data_nd_nc[i]);
    printf("{%12.7g, %12.7g}", ND_dataNC_minos->GetXaxis()->GetBinCenter(i+1),
           pred_nc->data[i+n_bins_reco_fd_nc] + data_nd_nc[i]);
    printf("},\n");

    printf("(* MINOS/MINOS+ 2017: FD CC *)\n");
    printf("{\n");
    for (i=0; i < n_bins_reco_fd_cc-1; i++)
      printf("{%12.7g, %12.7g},\n", FD_dataCC_minos->GetXaxis()->GetBinCenter(i+1),
             pred_cc->data[i] + data_fd_cc[i]);
    printf("{%12.7g, %12.7g}", FD_dataCC_minos->GetXaxis()->GetBinCenter(i+1),
           pred_cc->data[i] + data_fd_cc[i]);
    printf("},\n");

    printf("(* MINOS/MINOS+ 2017: ND CC *)\n");
    printf("{\n");
    for (i=0; i < n_bins_reco_nd_cc-1; i++)
      printf("{%12.7g, %12.7g},\n", ND_dataCC_minos->GetXaxis()->GetBinCenter(i+1),
             pred_cc->data[i+n_bins_reco_fd_cc] + data_nd_cc[i]);
    printf("{%12.7g, %12.7g}", ND_dataCC_minos->GetXaxis()->GetBinCenter(i+1),
           pred_cc->data[i+n_bins_reco_fd_cc] + data_nd_cc[i]);
    printf("} }\n");

    printf("# chi^2 = %12.7g\n", chi2);
  }

  return chi2;
}


