/***************************************************************************
 * Functions for MiniBooNE fit                                             *
 ***************************************************************************
 * Author: Joachim Kopp                                                    *
 ***************************************************************************/
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <globes/globes.h>
#include "nu.h"

// Use MC events to compute signal rate (very slow)?
//#define USE_MC_EVENTS

#define NE    11            /* Number of E_nue bins      */
#define NMU    8            /* Number of E_numue bins    */
#define NCOV  (2*NE+NMU)    /* Size of covariance matrix */

static gsl_matrix *M3        = NULL;
static gsl_matrix *M2        = NULL;
static gsl_matrix *M2inv     = NULL;
static gsl_permutation *perm = NULL;
#ifdef USE_MC_EVENTS
  static long Nevents;
  static double mc_events[128000][4];
  static glb_params osc_params = NULL;
#endif


/***************************************************************************
 * Initialize GSL data structures required for MiniBooNE chi^2 calculation *
 ***************************************************************************/
int chiMB_init()
{
  M3    = gsl_matrix_alloc(NCOV, NCOV);
  M2    = gsl_matrix_alloc(NE+NMU, NE+NMU);
  M2inv = gsl_matrix_alloc(NE+NMU, NE+NMU);
  perm  = gsl_permutation_alloc(NE+NMU);
#ifdef USE_MC_EVENTS
  osc_params = glbAllocParams();
  FILE *f = fopen("glb/miniboone/data/miniboone_numubarnuebarfullosc_ntuple.txt", "r");
//  FILE *f = fopen("glb/miniboone/data/miniboone_nunubarfullosc_ntuple.txt", "r");
  Nevents = 0;
  while (!feof(f))
  {
    int status;
    status = fscanf(f, "%lg %lg %lg %lg", &mc_events[Nevents][0], &mc_events[Nevents][1],
                    &mc_events[Nevents][2], &mc_events[Nevents][3]);
    if (status == EOF)
      break;
    Nevents++;
  }
  if (f) fclose(f);
#endif
  return 0;
}


/***************************************************************************
 * Cleanup GSL data structures required for MiniBooNE chi^2 calculation    *
 ***************************************************************************/
int chiMB_clear()
{
#ifdef USE_MC_EVENTS
  if (osc_params) { glbFreeParams(osc_params); osc_params = NULL; }
#endif
  if (perm)  { gsl_permutation_free(perm); perm = NULL; }
  if (M2inv) { gsl_matrix_free(M2inv);     M2   = NULL; }
  if (M2)    { gsl_matrix_free(M2);        M2   = NULL; }
  if (M3)    { gsl_matrix_free(M3);        M3   = NULL; }
  return 0;
}


/***************************************************************************
 * Calculate chi^2 for the MiniBooNE anti-neutrino analysis described at   *
 * http://www-boone.fnal.gov/for_physicists/data_release/nuebar2010/       *
 ***************************************************************************/
double chiMBanti_nu2010(int exp, int rule, int n_params, double *x, double *errors,
                        void *user_data)
{
  const double data[NE+NMU] = { 49, 41, 29, 25, 39, 21, 13, 12, 12, 12, 24,     /* \bar\nu_e   */
                                4583, 6715, 5451, 3783, 2253, 1145, 462, 176 }; /* \bar\nu_\mu */
  const double pred_numu[NMU] = { 4952, 6637, 5537, 3715, 2115, 1042, 439.3, 167.5 };
#include "glb/miniboone/cov_matrix.dat"

  double (*_M3)[NCOV]      = (double (*)[NCOV])   gsl_matrix_ptr(M3,    0, 0);
  double (*_M2)[NE+NMU]    = (double (*)[NE+NMU]) gsl_matrix_ptr(M2,    0, 0);
  double (*_M2inv)[NE+NMU] = (double (*)[NE+NMU]) gsl_matrix_ptr(M2inv, 0, 0);
  double *sig  = glbGetSignalFitRatePtr(exp, rule);
  double *bg   = glbGetBGFitRatePtr(exp, rule);
  int i, j;

  /* Compute signal using Monte Carlo events */
#ifdef USE_MC_EVENTS
  double *binwidths = glbGetBinSizeListPtr(exp);
  glbGetOscillationParameters(osc_params);
  double dmsq  = glbGetOscParams(osc_params, GLB_DM_21);
  double s22th = SQR(sin(2.0 * glbGetOscParams(osc_params, GLB_THETA_12)));
  for (int i=0; i < NE; i++)
    sig[i] = 0.0;
  for (long k=0; k < Nevents; k++)
  {
    double Ereco, Etrue, L, w;
    Ereco = mc_events[k][0];
    Etrue = mc_events[k][1];
    L     = mc_events[k][2];
    w     = mc_events[k][3];
    double prob = sin(1.267*0.01*dmsq*L/Etrue);
    prob = s22th * SQR(prob);
    double Eup = 200.0 + 1000.*binwidths[0];
    i = 0;
    while (Ereco >= Eup && i < NE)
      Eup += 1000.*binwidths[++i];
    if (i < NE)
      sig[i] += prob * w;
  }
  for (int i=0; i < NE; i++)
    sig[i] /= Nevents;
#endif

  /* Construct covariance matrix in 3x3 block form */
  double P[NCOV];         /* Vector of predicted nu_e signal, nu_e bg, nu_mu signal */
  P[0] = Mfrac_nubar[0][0]; /* This is too avoid a compiler warning about an unsed variable */
  P[0] = Mfrac_nunubar[0][0];
  for (int i=0; i < NE; i++)
  {
    P[i]    = sig[i];
    P[i+NE] = bg[i];
  }
  for (int i=0; i < NMU; i++)
    P[i+2*NE] = pred_numu[i];
  for (int i=0; i < NCOV; i++)
    for (int j=0; j < NCOV; j++)
      _M3[i][j] = Mfrac_nubar[i][j] * P[i] * P[j];
  for (i=0; i < NE; i++)        /* Add statistical error of signal */
    _M3[i][i] += sig[i];

  /* Collapse covariance matrix to 2x2 block form */
  for (i=0; i < NE; i++)        /* Upper left block */
    for (j=0; j < NE; j++)
      _M2[i][j] = _M3[i][j] + _M3[i+NE][j] + _M3[i][j+NE] + _M3[i+NE][j+NE];
  for (i=0; i < NMU; i++)       /* Lower left block */
    for (j=0; j < NE; j++)
      _M2[i+NE][j] = _M3[i+2*NE][j] + _M3[i+2*NE][j+NE];
  for (i=0; i < NE; i++)       /* Upper right block */
    for (j=0; j < NMU; j++)
      _M2[i][j+NE] = _M3[i][j+2*NE] + _M3[i+NE][j+2*NE];
  for (i=0; i < NMU; i++)      /* Lower right block */
    for (j=0; j < NMU; j++)
      _M2[i+NE][j+NE] = _M3[i+2*NE][j+2*NE];

  /* Invert covariance matrix and compute log-likelihood */
  int signum;
  gsl_linalg_LU_decomp(M2, perm, &signum);
  gsl_linalg_LU_invert(M2, perm, M2inv);

  double chi2 = 0.0;
  double P2[NE+NMU];
  for (i=0; i < NE; i++)
    P2[i] = data[i] - (sig[i] + bg[i]);
  for (i=0; i < NMU; i++)
    P2[i+NE] = data[i+NE] - pred_numu[i];
  for (i=0; i < NE+NMU; i++)
    for (j=0; j < NE+NMU; j++)
      chi2 += P2[i] * _M2inv[i][j] * P2[j];
  chi2 += gsl_linalg_LU_lndet(M2);

  return chi2;
}


