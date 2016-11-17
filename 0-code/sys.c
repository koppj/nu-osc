/***************************************************************************
 * Chi^2 functions for reactor and superbeam experiments                   *
 ***************************************************************************
 * Author: Joachim Kopp, Toshihiko Ota                                     *
 ***************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <globes/globes.h>   /* GLoBES library */
#include "nu.h"
#include "const.h"

double sys_startval_beam_plus[MAX_SYS];
double sys_startval_beam_minus[MAX_SYS];
double sys_startval_reactor[MAX_SYS];

extern double true_ldm;

/* Square of real number */
inline double square(double x)
{
  return x*x;
} 

/* Gauss likelihood (this is sufficient for reactor experiments due to the large event
   numbers; for other setups, one should use Poisson statistics) */
inline double gauss_likelihood(double true_rate, double fit_rate, double sqr_sigma)
{
  if (sqr_sigma > 0)
    return square(true_rate - fit_rate) / sqr_sigma;
  else
    return 0.0;
}

/* Poisson likelihood */
inline double poisson_likelihood(double true_rate, double fit_rate)
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


// -------------------------------------------------------------------------
void glbShiftAbsoluteEnergy(double b, double *rates_in, double *rates_out,
                            int n_bins, double emin, double emax)
// -------------------------------------------------------------------------
// Shifts the energy scale of all bins by a fixed amount of b GeV.
// -------------------------------------------------------------------------
{
  int i, k;
  double x = -b * n_bins/(emax - emin);
  double delta;
  for (i=0; i < n_bins; i++)
  {
    delta = x + i;
    k     = (int) floor(delta);

    if (k < -1 || k > n_bins - 1)
      rates_out[i] = 0.0;
    else if (k == -1)         // Assume out-of-bounds bins to contain 0 events
      rates_out[i] = rates_in[k+1] * (delta - k);
    else if (k == n_bins - 1)
      rates_out[i] = -rates_in[k] * (delta - k) + rates_in[k];
    else
      rates_out[i] = (rates_in[k+1] - rates_in[k]) * (delta - k) + rates_in[k];
  }
}


///***************************************************************************
// * Calculate chi^2 for T2K, including the following systematical errors:   *
// *   x[ 0]: Correlated beam flux normalization error                       *
// *   x[ 1]: Correlated normalization of intrinsic nu_e background          *
// *   x[ 2]: Uncorr. signal norm error for nu_mu QE events in the ND        *
// *   x[ 3]: Uncorr. signal norm error for nu_mu QE events in the FD        *
// *   x[ 4]: Uncorr. background norm error for nu_mu QE events in the ND    *
// *   x[ 5]: Uncorr. background norm error for nu_mu QE events in the FD    *
// *   x[ 6]: Uncorr. signal norm error for nu_e QE events in the ND         *
// *   x[ 7]: Uncorr. signal norm error for nu_e QE events in the FD         *
// *   x[ 8]: Uncorr. background norm error for nu_e QE events in the ND     *
// *   x[ 9]: Uncorr. background norm error for nu_e QE events in the FD     *
// *   x[10]: Uncorr. background tilt error for nu_e QE events in the ND     *
// *   x[11]: Uncorr. background tilt error for nu_e QE events in the FD     *
// *   x[12]: Uncorr. signal norm error for nu_e CC events in the ND         *
// *   x[13]: Uncorr. signal norm error for nu_e CC events in the FD         *
// *   x[14]: Uncorr. background norm error for nu_e CC events in the ND     *
// *   x[15]: Uncorr. background norm error for nu_e CC events in the FD     *
// ***************************************************************************/
//double chiT2K(int exp, int rule, int n_params, double *x, double *errors,
//              void *user_data)
//{
//  double *true_rates_mu_N, *true_rates_mu_F;
//  double *true_rates_QE_N, *true_rates_QE_F;
//  double *true_rates_CC_N, *true_rates_CC_F;
//  double *signal_fit_rates_mu_N, *signal_fit_rates_mu_F, *bg_fit_rates_mu_N, *bg_fit_rates_mu_F;
//  double *signal_fit_rates_QE_N, *signal_fit_rates_QE_F, *bg_fit_rates_QE_N, *bg_fit_rates_QE_F;
//  double *signal_fit_rates_CC_N, *signal_fit_rates_CC_F, *bg_fit_rates_CC_N, *bg_fit_rates_CC_F;
//  if (rule == RULE_T2K_NUMU_QE)
//  {
//    true_rates_mu_N       = glbGetRuleRatePtr(EXP_BEAM_NEAR, RULE_T2K_NUMU_QE);
//    true_rates_mu_F       = glbGetRuleRatePtr(EXP_BEAM_FAR,  RULE_T2K_NUMU_QE);
//    true_rates_QE_N       = glbGetRuleRatePtr(EXP_BEAM_NEAR, RULE_T2K_NUE_QE);
//    true_rates_QE_F       = glbGetRuleRatePtr(EXP_BEAM_FAR,  RULE_T2K_NUE_QE);
//    true_rates_CC_N       = glbGetRuleRatePtr(EXP_BEAM_NEAR, RULE_T2K_NUE_CC);
//    true_rates_CC_F       = glbGetRuleRatePtr(EXP_BEAM_FAR,  RULE_T2K_NUE_CC);
//    
//    signal_fit_rates_mu_N = glbGetSignalFitRatePtr(EXP_BEAM_NEAR, RULE_T2K_NUMU_QE);
//    signal_fit_rates_mu_F = glbGetSignalFitRatePtr(EXP_BEAM_FAR,  RULE_T2K_NUMU_QE);
//    bg_fit_rates_mu_N     = glbGetBGFitRatePtr    (EXP_BEAM_NEAR, RULE_T2K_NUMU_QE);
//    bg_fit_rates_mu_F     = glbGetBGFitRatePtr    (EXP_BEAM_FAR,  RULE_T2K_NUMU_QE);
//    signal_fit_rates_QE_N = glbGetSignalFitRatePtr(EXP_BEAM_NEAR, RULE_T2K_NUE_QE);
//    signal_fit_rates_QE_F = glbGetSignalFitRatePtr(EXP_BEAM_FAR,  RULE_T2K_NUE_QE);
//    bg_fit_rates_QE_N     = glbGetBGFitRatePtr    (EXP_BEAM_NEAR, RULE_T2K_NUE_QE);
//    bg_fit_rates_QE_F     = glbGetBGFitRatePtr    (EXP_BEAM_FAR,  RULE_T2K_NUE_QE);
//    signal_fit_rates_CC_N = glbGetSignalFitRatePtr(EXP_BEAM_NEAR, RULE_T2K_NUE_CC);
//    signal_fit_rates_CC_F = glbGetSignalFitRatePtr(EXP_BEAM_FAR,  RULE_T2K_NUE_CC);
//    bg_fit_rates_CC_N     = glbGetBGFitRatePtr    (EXP_BEAM_NEAR, RULE_T2K_NUE_CC);
//    bg_fit_rates_CC_F     = glbGetBGFitRatePtr    (EXP_BEAM_FAR,  RULE_T2K_NUE_CC);
//  }
//  else
//  {
//    true_rates_mu_N       = glbGetRuleRatePtr(EXP_BEAM_NEAR, RULE_T2K_NUMU_BAR_QE);
//    true_rates_mu_F       = glbGetRuleRatePtr(EXP_BEAM_FAR,  RULE_T2K_NUMU_BAR_QE);
//    true_rates_QE_N       = glbGetRuleRatePtr(EXP_BEAM_NEAR, RULE_T2K_NUE_BAR_QE);
//    true_rates_QE_F       = glbGetRuleRatePtr(EXP_BEAM_FAR,  RULE_T2K_NUE_BAR_QE);
//    true_rates_CC_N       = glbGetRuleRatePtr(EXP_BEAM_NEAR, RULE_T2K_NUE_BAR_CC);
//    true_rates_CC_F       = glbGetRuleRatePtr(EXP_BEAM_FAR,  RULE_T2K_NUE_BAR_CC);
//    
//    signal_fit_rates_mu_N = glbGetSignalFitRatePtr(EXP_BEAM_NEAR, RULE_T2K_NUMU_BAR_QE);
//    signal_fit_rates_mu_F = glbGetSignalFitRatePtr(EXP_BEAM_FAR,  RULE_T2K_NUMU_BAR_QE);
//    bg_fit_rates_mu_N     = glbGetBGFitRatePtr    (EXP_BEAM_NEAR, RULE_T2K_NUMU_BAR_QE);
//    bg_fit_rates_mu_F     = glbGetBGFitRatePtr    (EXP_BEAM_FAR,  RULE_T2K_NUMU_BAR_QE);
//    signal_fit_rates_QE_N = glbGetSignalFitRatePtr(EXP_BEAM_NEAR, RULE_T2K_NUE_BAR_QE);
//    signal_fit_rates_QE_F = glbGetSignalFitRatePtr(EXP_BEAM_FAR,  RULE_T2K_NUE_BAR_QE);
//    bg_fit_rates_QE_N     = glbGetBGFitRatePtr    (EXP_BEAM_NEAR, RULE_T2K_NUE_BAR_QE);
//    bg_fit_rates_QE_F     = glbGetBGFitRatePtr    (EXP_BEAM_FAR,  RULE_T2K_NUE_BAR_QE);
//    signal_fit_rates_CC_N = glbGetSignalFitRatePtr(EXP_BEAM_NEAR, RULE_T2K_NUE_BAR_CC);
//    signal_fit_rates_CC_F = glbGetSignalFitRatePtr(EXP_BEAM_FAR,  RULE_T2K_NUE_BAR_CC);
//    bg_fit_rates_CC_N     = glbGetBGFitRatePtr    (EXP_BEAM_NEAR, RULE_T2K_NUE_BAR_CC);
//    bg_fit_rates_CC_F     = glbGetBGFitRatePtr    (EXP_BEAM_FAR,  RULE_T2K_NUE_BAR_CC);
//  }
//  
//  double signal_norm_mu_N, bg_norm_mu_N, signal_norm_mu_F, bg_norm_mu_F;
//  double signal_norm_QE_N, bg_norm_QE_N, signal_norm_QE_F, bg_norm_QE_F;
//  double signal_norm_CC_N, bg_norm_CC_N, signal_norm_CC_F, bg_norm_CC_F;
//  double bg_tilt_QE_N, bg_tilt_QE_F;
//  double emin, emax, ecenter;
//  double *bin_centers;
//  int ew_low, ew_high;
//  double fit_rate;
//  double total_fit_rate_CC_N, total_fit_rate_CC_F;
//  double total_true_rate_CC_N, total_true_rate_CC_F;
//  double chi2 = 0.0;
//  int i;
//
//
//  /* Request energy window */
//  glbGetEminEmax(exp, &emin, &emax);
//  glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
//  bin_centers = glbGetBinCentersListPtr(exp);
//  ecenter = 0.5 * (emax + emin);
//  bg_tilt_QE_N = x[10] / (emax - emin); 
//  bg_tilt_QE_F = x[11] / (emax - emin); 
//  
//  /* Calculate nuisance parameters */
//  signal_norm_mu_N = 1.0 + x[0] + x[2];
//  signal_norm_mu_F = 1.0 + x[0] + x[3];
//  bg_norm_mu_N     = 1.0 + x[0] + x[4];
//  bg_norm_mu_F     = 1.0 + x[0] + x[5];
//  signal_norm_QE_N = 1.0 + x[0] + x[6];
//  signal_norm_QE_F = 1.0 + x[0] + x[7];
//  bg_norm_QE_N     = 1.0 + x[0] + x[1] + x[8];
//  bg_norm_QE_F     = 1.0 + x[0] + x[1] + x[9];
//  signal_norm_CC_N = 1.0 + x[0] + x[12];
//  signal_norm_CC_F = 1.0 + x[0] + x[13];
//  bg_norm_CC_N     = 1.0 + x[0] + x[1] + x[14];
//  bg_norm_CC_F     = 1.0 + x[0] + x[1] + x[15];
//
//  /* Loop over all bins in energy window */
//  total_fit_rate_CC_N  = 0.0;
//  total_fit_rate_CC_F  = 0.0;
//  total_true_rate_CC_N = 0.0;
//  total_true_rate_CC_F = 0.0;
//  for (i=ew_low; i <= ew_high; i++)
//  {
//    /* Spectral analysis for QE nu_mu disappearance events */
//    fit_rate = signal_norm_mu_N * signal_fit_rates_mu_N[i] + bg_norm_mu_N * bg_fit_rates_mu_N[i];
//    chi2    += poisson_likelihood(true_rates_mu_N[i], fit_rate);
//
//    fit_rate = signal_norm_mu_F * signal_fit_rates_mu_F[i] + bg_norm_mu_F * bg_fit_rates_mu_F[i];
//    chi2    += poisson_likelihood(true_rates_mu_F[i], fit_rate);
//
//    /* Spectral analysis for QE nu_e appearance events */
//    fit_rate = signal_norm_QE_N * signal_fit_rates_QE_N[i]
//      + bg_norm_QE_N * bg_fit_rates_QE_N[i]
//      + bg_tilt_QE_N * (bin_centers[i]-ecenter) * bg_fit_rates_QE_N[i];
//    chi2 += poisson_likelihood(true_rates_QE_N[i], fit_rate);
//
//    fit_rate = signal_norm_QE_F * signal_fit_rates_QE_F[i]
//      + bg_norm_QE_F * bg_fit_rates_QE_F[i]
//      + bg_tilt_QE_F * (bin_centers[i]-ecenter) * bg_fit_rates_QE_F[i];
//    chi2 += poisson_likelihood(true_rates_QE_F[i], fit_rate);
//
//    /* Total rates analysis for CC nu_e appearance events */
//    total_fit_rate_CC_N  += signal_norm_CC_N * signal_fit_rates_CC_N[i]
//                                +  bg_norm_CC_N * bg_fit_rates_CC_N[i];
//    total_true_rate_CC_N += true_rates_CC_N[i];
//
//    total_fit_rate_CC_F  += signal_norm_CC_F * signal_fit_rates_CC_F[i]
//                                +  bg_norm_CC_F * bg_fit_rates_CC_F[i];
//    total_true_rate_CC_F += true_rates_CC_F[i];
//  }
//  chi2 += poisson_likelihood(total_true_rate_CC_N, total_fit_rate_CC_N);
//  chi2 += poisson_likelihood(total_true_rate_CC_F, total_fit_rate_CC_F);
//
//  /* Systematics priors */
//  for (i=0; i < n_params; i++)
//    chi2 += square(x[i] / errors[i]);
//
//  /* Save the systematics parameters as starting values for the next step */
//  if (rule == RULE_T2K_NUMU_QE)
//    for (i=0; i < n_params; i++)
//      sys_startval_beam_plus[i] = x[i];
//  else
//    for (i=0; i < n_params; i++)
//      sys_startval_beam_minus[i] = x[i];
//  
//  return chi2;
//}
//

/***************************************************************************
 * Calculate chi^2 for NOvA, including the following systematical errors:  *
 *   x[0]:         Correlated beam flux normalization error                *
 *   x[1]:         Correlated normalization of intrinsic nu_e background   *
 *   x[ 2], x[ 3]: Uncorr. signal norm and cal. errors for ND nu_e events  *
 *   x[ 4], x[ 5]: Uncorr. signal norm and cal. errors for FD nu_e events  *
 *   x[ 6], x[ 7]: Uncorr. BG norm and cal. errors for ND nu_e events      *
 *   x[ 8], x[ 9]: Uncorr. BG norm and cal. errors for FD nu_e events      *
 *   x[10], x[11]: Uncorr. signal norm and cal. errors for ND nu_mu events *
 *   x[12], x[13]: Uncorr. signal norm and cal. errors for FD nu_mu events *
 *   x[14], x[15]: Uncorr. BG norm and cal. errors for ND nu_mu events     *
 *   x[16], x[17]: Uncorr. BG norm and cal. errors for FD nu_mu events     *
 ***************************************************************************/
double chiNOvA(int exp, int rule, int n_params, double *x, double *errors,
              void *user_data)
{
  int n_bins = glbGetNumberOfBins(EXP_BEAM_FAR);
  double *true_rates_e_N, *true_rates_e_F, *true_rates_mu_N, *true_rates_mu_F;
  double signal_fit_rates_e_N[n_bins], signal_fit_rates_e_F[n_bins];
  double bg_fit_rates_e_N[n_bins], bg_fit_rates_e_F[n_bins];
  double signal_fit_rates_mu_N[n_bins], signal_fit_rates_mu_F[n_bins];
  double bg_fit_rates_mu_N[n_bins], bg_fit_rates_mu_F[n_bins];
  double signal_norm_e_N, bg_norm_e_N, signal_norm_e_F, bg_norm_e_F;
  double signal_norm_mu_N, bg_norm_mu_N, signal_norm_mu_F, bg_norm_mu_F;
  int ew_low, ew_high;
  double emin, emax;
  double fit_rate;
  double chi2 = 0.0;
  int i;

  /* Request simulated energy interval and analysis energy window */
  glbGetEminEmax(exp, &emin, &emax);
  glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
   
  if (rule == RULE_NOVA_NUE)
  {
    /* Get true event rates */
    true_rates_e_N        = glbGetRuleRatePtr(EXP_BEAM_NEAR, RULE_NOVA_NUE);
    true_rates_e_F        = glbGetRuleRatePtr(EXP_BEAM_FAR,  RULE_NOVA_NUE);
    true_rates_mu_N       = glbGetRuleRatePtr(EXP_BEAM_NEAR, RULE_NOVA_NUMU);
    true_rates_mu_F       = glbGetRuleRatePtr(EXP_BEAM_FAR,  RULE_NOVA_NUMU);

    /* Apply energy calibration errors */
    glbShiftEnergyScale(x[ 3], glbGetSignalFitRatePtr(EXP_BEAM_NEAR, RULE_NOVA_NUE),
                        signal_fit_rates_e_N, n_bins, emin, emax);
    glbShiftEnergyScale(x[ 5], glbGetSignalFitRatePtr(EXP_BEAM_FAR,  RULE_NOVA_NUE),
                        signal_fit_rates_e_F, n_bins, emin, emax);
    glbShiftEnergyScale(x[ 7], glbGetBGFitRatePtr    (EXP_BEAM_NEAR, RULE_NOVA_NUE),
                        bg_fit_rates_e_N, n_bins, emin, emax);
    glbShiftEnergyScale(x[ 9], glbGetBGFitRatePtr    (EXP_BEAM_FAR,  RULE_NOVA_NUE),
                        bg_fit_rates_e_F, n_bins, emin, emax);
    glbShiftEnergyScale(x[11], glbGetSignalFitRatePtr(EXP_BEAM_NEAR, RULE_NOVA_NUMU),
                        signal_fit_rates_mu_N, n_bins, emin, emax);
    glbShiftEnergyScale(x[13], glbGetSignalFitRatePtr(EXP_BEAM_FAR,  RULE_NOVA_NUMU),
                        signal_fit_rates_mu_F, n_bins, emin, emax);
    glbShiftEnergyScale(x[15], glbGetBGFitRatePtr    (EXP_BEAM_NEAR, RULE_NOVA_NUMU),
                        bg_fit_rates_mu_N, n_bins, emin, emax);
    glbShiftEnergyScale(x[17], glbGetBGFitRatePtr    (EXP_BEAM_FAR,  RULE_NOVA_NUMU),
                        bg_fit_rates_mu_F, n_bins, emin, emax);
  }
  else
  {
    /* Get true event rates */
    true_rates_e_N        = glbGetRuleRatePtr(EXP_BEAM_NEAR, RULE_NOVA_NUE_BAR);
    true_rates_e_F        = glbGetRuleRatePtr(EXP_BEAM_FAR,  RULE_NOVA_NUE_BAR);
    true_rates_mu_N       = glbGetRuleRatePtr(EXP_BEAM_NEAR, RULE_NOVA_NUMU_BAR);
    true_rates_mu_F       = glbGetRuleRatePtr(EXP_BEAM_FAR,  RULE_NOVA_NUMU_BAR);
    
    /* Apply energy calibration errors */
    glbShiftEnergyScale(x[ 3], glbGetSignalFitRatePtr(EXP_BEAM_NEAR, RULE_NOVA_NUE_BAR),
                        signal_fit_rates_e_N, n_bins, emin, emax);
    glbShiftEnergyScale(x[ 5], glbGetSignalFitRatePtr(EXP_BEAM_FAR,  RULE_NOVA_NUE_BAR),
                        signal_fit_rates_e_F, n_bins, emin, emax);
    glbShiftEnergyScale(x[ 7], glbGetBGFitRatePtr    (EXP_BEAM_NEAR, RULE_NOVA_NUE_BAR),
                        bg_fit_rates_e_N, n_bins, emin, emax);
    glbShiftEnergyScale(x[ 9], glbGetBGFitRatePtr    (EXP_BEAM_FAR,  RULE_NOVA_NUE_BAR),
                        bg_fit_rates_e_F, n_bins, emin, emax);
    glbShiftEnergyScale(x[11], glbGetSignalFitRatePtr(EXP_BEAM_NEAR, RULE_NOVA_NUMU_BAR),
                        signal_fit_rates_mu_N, n_bins, emin, emax);
    glbShiftEnergyScale(x[13], glbGetSignalFitRatePtr(EXP_BEAM_FAR,  RULE_NOVA_NUMU_BAR),
                        signal_fit_rates_mu_F, n_bins, emin, emax);
    glbShiftEnergyScale(x[15], glbGetBGFitRatePtr    (EXP_BEAM_NEAR, RULE_NOVA_NUMU_BAR),
                        bg_fit_rates_mu_N, n_bins, emin, emax);
    glbShiftEnergyScale(x[17], glbGetBGFitRatePtr    (EXP_BEAM_FAR,  RULE_NOVA_NUMU_BAR),
                        bg_fit_rates_mu_F, n_bins, emin, emax);
  }
  
   /* Calculate nuisance parameters */
  signal_norm_e_N  = 1.0 + x[0] + x[ 2];
  signal_norm_e_F  = 1.0 + x[0] + x[ 4];
  bg_norm_e_N      = 1.0 + x[0] + x[1] + x[ 6];
  bg_norm_e_F      = 1.0 + x[0] + x[1] + x[ 8];
  signal_norm_mu_N = 1.0 + x[0] + x[10];
  signal_norm_mu_F = 1.0 + x[0] + x[12];
  bg_norm_mu_N     = 1.0 + x[0] + x[14];
  bg_norm_mu_F     = 1.0 + x[0] + x[16];

  /* Loop over all bins in energy window */
  for (i=ew_low; i <= ew_high; i++)
  {
    /* Statistical contribution from the nu_e channels */
    fit_rate = signal_norm_e_N * signal_fit_rates_e_N[i] + bg_norm_e_N * bg_fit_rates_e_N[i];
    chi2    += poisson_likelihood(true_rates_e_N[i], fit_rate);

    fit_rate = signal_norm_e_F * signal_fit_rates_e_F[i] + bg_norm_e_F * bg_fit_rates_e_F[i];
    chi2    += poisson_likelihood(true_rates_e_F[i], fit_rate);

    /* Statistical contribution from the nu_mu channels */
    fit_rate = signal_norm_mu_N * signal_fit_rates_mu_N[i] + bg_norm_mu_N * bg_fit_rates_mu_N[i];
    chi2    += poisson_likelihood(true_rates_mu_N[i], fit_rate);

    fit_rate = signal_norm_mu_F * signal_fit_rates_mu_F[i] + bg_norm_mu_F * bg_fit_rates_mu_F[i];
    chi2    += poisson_likelihood(true_rates_mu_F[i], fit_rate);
  }

  /* Systematics priors */
  for (i=0; i < n_params; i++)
    chi2 += square(x[i] / errors[i]);

  /* Save the systematics parameters as starting values for the next step */
  if (rule == RULE_NOVA_NUE)
    for (i=0; i < n_params; i++)
      sys_startval_beam_plus[i] = x[i];
  else
    for (i=0; i < n_params; i++)
      sys_startval_beam_minus[i] = x[i];

  return chi2;
}


///***************************************************************************
// * Calculate chi^2 for WBB with WC detector, including the following       *
// * systematical errors:                                                    *
// *   x[ 0]:        Correlated beam flux normalization error                *
// *   x[ 1]:        Correlated normalization of background in nu_e sample   *
// *   x[ 2]:        Correlated tilt of background in nu_e sample            *
// *   x[ 3], x[ 4]: Uncorr. signal norm/tilt for ND nu_e events             *
// *   x[ 5], x[ 6]: Uncorr. signal norm/tilt for FD nu_e events             *
// *   x[ 7], x[ 8]: Uncorr. BG norm/tilt errors for ND nu_e events          *
// *   x[ 9], x[10]: Uncorr. BG norm/tilt errors for FD nu_e events          *
// *   x[11], x[12]: Uncorr. signal norm/tilt errors for ND nu_mu events     *
// *   x[13], x[14]: Uncorr. signal norm/tilt errors for FD nu_mu events     *
// *   x[15], x[16]: Uncorr. BG norm/tilt errors for ND nu_mu events         *
// *   x[17], x[18]: Uncorr. BG norm/tilt errors for FD nu_mu events         *
// * and a fixed 1.0% bin-to-bin error in the near detector to avoid overly  *
// * optimistic NSI sensitivities                                            *
// *                                                                         *
// * IMPORTANT: The function assumed the ND experiment to be loaded directly *
// *            before the FD experiment, so that exp points to the FD,      *
// *            and exp-1 to the ND                                          *
// ***************************************************************************/
//double chiWBB_WC(int exp, int rule, int n_params, double *x, double *errors,
//                 void *user_data)
//{
//  double *true_rates_e_N, *true_rates_e_F, *true_rates_mu_N, *true_rates_mu_F;
//  double *sig_fit_rates_e_N, *sig_fit_rates_e_F;
//  double *bg_fit_rates_e_N, *bg_fit_rates_e_F;
//  double *sig_fit_rates_mu_N, *sig_fit_rates_mu_F;
//  double *bg_fit_rates_mu_N, *bg_fit_rates_mu_F;
//  double sig_norm_e_N, bg_norm_e_N, sig_norm_e_F, bg_norm_e_F;
//  double sig_norm_mu_N, bg_norm_mu_N, sig_norm_mu_F, bg_norm_mu_F;
//  double sig_tilt_e_N, bg_tilt_e_N, sig_tilt_e_F, bg_tilt_e_F;
//  double sig_tilt_mu_N, bg_tilt_mu_N, sig_tilt_mu_F, bg_tilt_mu_F;
//
//  int ew_low, ew_high;
//  double emin, emax, inv_delta_e, e_center;
//  double *bin_centers;
//
//  const int EXP_BEAM_NEAR = exp - 1;
//  const int EXP_BEAM_FAR  = exp;
//  const double sigma_binbin = 0.01;
//  double fit_rate;
//  double chi2 = 0.0;
//  int i;
//
//  /* Request simulated energy interval and analysis energy window */
//  glbGetEminEmax(exp, &emin, &emax);
//  glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
//  bin_centers = glbGetBinCentersListPtr(exp);
//  inv_delta_e = 1.0 / (emax - emin);
//  e_center    = 0.5 * (emax + emin);
//   
//  /* Get true event rates */
//  true_rates_e_N     = glbGetRuleRatePtr(EXP_BEAM_NEAR, RULE_WBB_WC_NUE);
//  true_rates_e_F     = glbGetRuleRatePtr(EXP_BEAM_FAR,  RULE_WBB_WC_NUE);
//  true_rates_mu_N    = glbGetRuleRatePtr(EXP_BEAM_NEAR, RULE_WBB_WC_NUMU);
//  true_rates_mu_F    = glbGetRuleRatePtr(EXP_BEAM_FAR,  RULE_WBB_WC_NUMU);
//
//  /* Get simulated event rates */
//  sig_fit_rates_e_N  = glbGetSignalFitRatePtr(EXP_BEAM_NEAR, RULE_WBB_WC_NUE);
//  sig_fit_rates_e_F  = glbGetSignalFitRatePtr(EXP_BEAM_FAR,  RULE_WBB_WC_NUE);
//  bg_fit_rates_e_N   = glbGetBGFitRatePtr    (EXP_BEAM_NEAR, RULE_WBB_WC_NUE);
//  bg_fit_rates_e_F   = glbGetBGFitRatePtr    (EXP_BEAM_FAR,  RULE_WBB_WC_NUE);
//  sig_fit_rates_mu_N = glbGetSignalFitRatePtr(EXP_BEAM_NEAR, RULE_WBB_WC_NUMU);
//  sig_fit_rates_mu_F = glbGetSignalFitRatePtr(EXP_BEAM_FAR,  RULE_WBB_WC_NUMU);
//  bg_fit_rates_mu_N  = glbGetBGFitRatePtr    (EXP_BEAM_NEAR, RULE_WBB_WC_NUMU);
//  bg_fit_rates_mu_F  = glbGetBGFitRatePtr    (EXP_BEAM_FAR,  RULE_WBB_WC_NUMU);
//
//  /* Calculate nuisance parameters */
//  sig_norm_e_N  = 1.0 + x[0] + x[3];
//  sig_norm_e_F  = 1.0 + x[0] + x[5];
//  bg_norm_e_N   = 1.0 + x[0] + x[1] + x[7];
//  bg_norm_e_F   = 1.0 + x[0] + x[1] + x[9];
//  sig_norm_mu_N = 1.0 + x[0] + x[11];
//  sig_norm_mu_F = 1.0 + x[0] + x[13];
//  bg_norm_mu_N  = 1.0 + x[0] + x[15];
//  bg_norm_mu_F  = 1.0 + x[0] + x[17];
//
//  sig_tilt_e_N  = x[ 4] * inv_delta_e;
//  sig_tilt_e_F  = x[ 6] * inv_delta_e;
//  bg_tilt_e_N   = (x[ 2] + x[ 8]) * inv_delta_e;
//  bg_tilt_e_F   = (x[ 2] + x[10]) * inv_delta_e;
//  sig_tilt_mu_N = x[12] * inv_delta_e;
//  sig_tilt_mu_F = x[14] * inv_delta_e;
//  bg_tilt_mu_N  = x[16] * inv_delta_e;
//  bg_tilt_mu_F  = x[18] * inv_delta_e;
//
//  /* Loop over all bins in energy window */
//  for (i=ew_low; i <= ew_high; i++)
//  {
//    double x = bin_centers[i] - e_center;
//
//    /* Statistical contribution from the nu_e channels */
//    fit_rate = (sig_norm_e_N + x * sig_tilt_e_N) * sig_fit_rates_e_N[i]
//                + (bg_norm_e_N + x * bg_tilt_e_N) * bg_fit_rates_e_N[i];
//    chi2 += gauss_likelihood( true_rates_e_N[i], fit_rate,
//              true_rates_e_N[i] * (1.0 + true_rates_e_N[i]*square(sigma_binbin)) );
//
//    fit_rate = (sig_norm_e_F + x * sig_tilt_e_F) * sig_fit_rates_e_F[i]
//                + (bg_norm_e_F + x * bg_tilt_e_F) * bg_fit_rates_e_F[i];
//    chi2    += poisson_likelihood(true_rates_e_F[i], fit_rate);
//
//    /* Statistical contribution from the nu_mu channels */
//    fit_rate = (sig_norm_mu_N + x * sig_tilt_mu_N) * sig_fit_rates_mu_N[i]
//                + (bg_norm_mu_N + x * bg_tilt_mu_N) * bg_fit_rates_mu_N[i];
//    chi2 += gauss_likelihood( true_rates_mu_N[i], fit_rate,
//              true_rates_mu_N[i] * (1.0 + true_rates_mu_N[i]*square(sigma_binbin)) );
//
//    fit_rate = (sig_norm_mu_F + x * sig_tilt_mu_F) * sig_fit_rates_mu_F[i]
//                 + (bg_norm_mu_F + x * bg_tilt_mu_F) * bg_fit_rates_mu_F[i];
//    chi2    += poisson_likelihood(true_rates_mu_F[i], fit_rate);
//  }
//
//  /* Systematics priors */
//  for (i=0; i < n_params; i++)
//    chi2 += square(x[i] / errors[i]);
//
//  return chi2;
//}


/***************************************************************************
 * A faster version of chiWBB_WC, omitting the tilt errors                 *
 *   x[0]: Correlated beam flux normalization error                        *
 *   x[1]: Correlated normalization of background in nu_e sample           *
 *   x[2]: Uncorr. signal norm for ND nu_e events                          *
 *   x[3]: Uncorr. signal norm for FD nu_e events                          *
 *   x[4]: Uncorr. BG norm errors for ND nu_e events                       *
 *   x[5]: Uncorr. BG norm errors for FD nu_e events                       *
 *   x[6]: Uncorr. signal norm errors for ND nu_mu events                  *
 *   x[7]: Uncorr. signal norm errors for FD nu_mu events                  *
 *   x[8]: Uncorr. BG norm errors for ND nu_mu events                      *
 *   x[9]: Uncorr. BG norm errors for FD nu_mu events                      *
 * and a fixed 1.0% bin-to-bin error in the near detector to avoid overly  *
 * optimistic NSI sensitivities                                            *
 *                                                                         *
 * IMPORTANT: The function assumed the ND experiment to be loaded directly *
 *            before the FD experiment, so that exp points to the FD,      *
 *            and exp-1 to the ND                                          *
 ***************************************************************************/
double chiWBB_WCfast(int exp, int rule, int n_params, double *x, double *errors,
                     void *user_data)
{
  wbb_params_type *wbb_params = (wbb_params_type *) user_data;
  double *sig_true_rates_e_N, *sig_true_rates_e_F;
  double *bg_true_rates_e_N, *bg_true_rates_e_F;
  double *sig_true_rates_mu_N, *sig_true_rates_mu_F;
  double *bg_true_rates_mu_N, *bg_true_rates_mu_F;

  double *sig_fit_rates_e_N, *sig_fit_rates_e_F;
  double *bg_fit_rates_e_N, *bg_fit_rates_e_F;
  double *sig_fit_rates_mu_N, *sig_fit_rates_mu_F;
  double *bg_fit_rates_mu_N, *bg_fit_rates_mu_F;
  double sig_norm_e_N, bg_norm_e_N, sig_norm_e_F, bg_norm_e_F;
  double sig_norm_mu_N, bg_norm_mu_N, sig_norm_mu_F, bg_norm_mu_F;

  int ew_low, ew_high;
  double *bin_centers;

  const int EXP_BEAM_NEAR = exp - 1;
  const int EXP_BEAM_FAR  = exp;
  const double sigma_binbin = 0.01;

  
  /* Get true event rates */
  sig_true_rates_e_N  = glbGetSignalRatePtr(EXP_BEAM_NEAR, RULE_WBB_WC_NUE);
  sig_true_rates_e_F  = glbGetSignalRatePtr(EXP_BEAM_FAR,  RULE_WBB_WC_NUE);
  bg_true_rates_e_N   = glbGetBGRatePtr    (EXP_BEAM_NEAR, RULE_WBB_WC_NUE);
  bg_true_rates_e_F   = glbGetBGRatePtr    (EXP_BEAM_FAR,  RULE_WBB_WC_NUE);
  sig_true_rates_mu_N = glbGetSignalRatePtr(EXP_BEAM_NEAR, RULE_WBB_WC_NUMU);
  sig_true_rates_mu_F = glbGetSignalRatePtr(EXP_BEAM_FAR,  RULE_WBB_WC_NUMU);
  bg_true_rates_mu_N  = glbGetBGRatePtr    (EXP_BEAM_NEAR, RULE_WBB_WC_NUMU);
  bg_true_rates_mu_F  = glbGetBGRatePtr    (EXP_BEAM_FAR,  RULE_WBB_WC_NUMU);

  /* Get simulated event rates */
  sig_fit_rates_e_N   = glbGetSignalFitRatePtr(EXP_BEAM_NEAR, RULE_WBB_WC_NUE);
  sig_fit_rates_e_F   = glbGetSignalFitRatePtr(EXP_BEAM_FAR,  RULE_WBB_WC_NUE);
  bg_fit_rates_e_N    = glbGetBGFitRatePtr    (EXP_BEAM_NEAR, RULE_WBB_WC_NUE);
  bg_fit_rates_e_F    = glbGetBGFitRatePtr    (EXP_BEAM_FAR,  RULE_WBB_WC_NUE);
  sig_fit_rates_mu_N  = glbGetSignalFitRatePtr(EXP_BEAM_NEAR, RULE_WBB_WC_NUMU);
  sig_fit_rates_mu_F  = glbGetSignalFitRatePtr(EXP_BEAM_FAR,  RULE_WBB_WC_NUMU);
  bg_fit_rates_mu_N   = glbGetBGFitRatePtr    (EXP_BEAM_NEAR, RULE_WBB_WC_NUMU);
  bg_fit_rates_mu_F   = glbGetBGFitRatePtr    (EXP_BEAM_FAR,  RULE_WBB_WC_NUMU);

  /* Calculate nuisance parameters */
  sig_norm_e_N  = 1.0 + x[0] + x[2];
  sig_norm_e_F  = 1.0 + x[0] + x[3];
  bg_norm_e_N   = 1.0 + x[0] + x[1] + x[4];
  bg_norm_e_F   = 1.0 + x[0] + x[1] + x[5];
  sig_norm_mu_N = 1.0 + x[0] + x[6];
  sig_norm_mu_F = 1.0 + x[0] + x[7];
  bg_norm_mu_N  = 1.0 + x[0] + x[8];
  bg_norm_mu_F  = 1.0 + x[0] + x[9];

  /* Subroutine to compute chi2 using spectral information */
  double wbb_wc_sum_chi2(int ew_low, int ew_high, double eff_factor)
  {
    if (ew_low > ew_high)
      return 0.0;

    int i;
    double fit_rate, obs_rate;
    double eff_sig=eff_factor, eff_bg=1.0;
    double chi2 = 0.0;

    /* Loop over all bins in energy window */
    for (i=ew_low; i <= ew_high; i++)
    {
      /* Spectral analysis for nu_mu disappearance events */
      fit_rate = eff_sig*sig_norm_mu_N*sig_fit_rates_mu_N[i] + eff_bg*bg_norm_mu_N*bg_fit_rates_mu_N[i];
      obs_rate = eff_sig*sig_true_rates_mu_N[i] + eff_bg*bg_true_rates_mu_N[i];
      chi2 += gauss_likelihood( obs_rate, fit_rate,
                obs_rate * (1.0 + obs_rate*square(sigma_binbin)) );

      fit_rate = eff_sig*sig_norm_mu_F*sig_fit_rates_mu_F[i] + eff_bg*bg_norm_mu_F*bg_fit_rates_mu_F[i];
      obs_rate = eff_sig*sig_true_rates_mu_F[i] + eff_bg*bg_true_rates_mu_F[i];
      chi2 += poisson_likelihood(obs_rate, fit_rate);

      /* Spectral analysis for nu_e appearance events */
      fit_rate = eff_sig*sig_norm_e_N*sig_fit_rates_e_N[i] + eff_bg*bg_norm_e_N*bg_fit_rates_e_N[i];
      obs_rate = eff_sig*sig_true_rates_e_N[i] + eff_bg*bg_true_rates_e_N[i];
      chi2 += gauss_likelihood( obs_rate, fit_rate,
                obs_rate * (1.0 + obs_rate*square(sigma_binbin)) );

      fit_rate = eff_sig*sig_norm_e_F*sig_fit_rates_e_F[i] + eff_bg*bg_norm_e_F*bg_fit_rates_e_F[i];
      obs_rate = eff_sig*sig_true_rates_e_F[i] + eff_bg*bg_true_rates_e_F[i];
      chi2 += poisson_likelihood(obs_rate, fit_rate);
    }

    return chi2;
  }

  /* Subroutine to compute chi2 using total rates */
  double wbb_wc_chi2_total_rates(int ew_low, int ew_high, double eff_factor)
  {
    if (ew_low > ew_high)
      return 0.0;

    int i;
    double fit_rate, obs_rate;
    double eff_sig=eff_factor, eff_bg=1.0;
    double chi2 = 0.0;

    fit_rate = obs_rate = 0.0;
    for (i=ew_low; i <= ew_high; i++)
    {
      fit_rate += eff_sig*sig_norm_mu_N*sig_fit_rates_mu_N[i]
                    + eff_bg*bg_norm_mu_N*bg_fit_rates_mu_N[i];
      obs_rate += eff_sig*sig_true_rates_mu_N[i] + eff_bg*bg_true_rates_mu_N[i];
    }
    chi2 += gauss_likelihood( obs_rate, fit_rate,
                obs_rate * (1.0 + obs_rate*square(sigma_binbin)) );

    fit_rate = obs_rate = 0.0;
    for (i=ew_low; i <= ew_high; i++)
    {
      fit_rate += eff_sig*sig_norm_mu_F*sig_fit_rates_mu_F[i]
                    + eff_bg*bg_norm_mu_F*bg_fit_rates_mu_F[i];
      obs_rate += eff_sig*sig_true_rates_mu_F[i] + eff_bg*bg_true_rates_mu_F[i];
    }
    chi2 += poisson_likelihood(obs_rate, fit_rate);

    fit_rate = obs_rate = 0.0;
    for (i=ew_low; i <= ew_high; i++)
    {
      fit_rate += eff_sig*sig_norm_e_N*sig_fit_rates_e_N[i]
                    + eff_bg*bg_norm_e_N*bg_fit_rates_e_N[i];
      obs_rate += eff_sig*sig_true_rates_e_N[i] + eff_bg*bg_true_rates_e_N[i];
    }
    chi2 += gauss_likelihood( obs_rate, fit_rate,
              obs_rate * (1.0 + obs_rate*square(sigma_binbin)) );

    fit_rate = obs_rate = 0.0;
    for (i=ew_low; i <= ew_high; i++)
    {
      fit_rate += eff_sig*sig_norm_e_F*sig_fit_rates_e_F[i]
                    + eff_bg*bg_norm_e_F*bg_fit_rates_e_F[i];
      obs_rate += eff_sig*sig_true_rates_e_F[i] + eff_bg*bg_true_rates_e_F[i];
    }
    chi2 += poisson_likelihood(obs_rate, fit_rate);

    return chi2;
  }

  /* Request simulated energy interval and analysis energy window */
  glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
  bin_centers = glbGetBinCentersListPtr(exp);

  /* Compute position of 1st oscillation minimum */
  double L = glbGetBaselineInExperiment(EXP_BEAM_FAR);
  double E_1st_min = L*KM * fabs(true_ldm) / (4*M_PI) / GEV;
  int i_1st_min = 0;
  while (bin_centers[i_1st_min] <= E_1st_min)
    i_1st_min++;
  double chi2 = 0.0;
  int i;

  /* Find out if this chi^2 computation is for neutrino or anti-neutrino running   */
  /* We do this by checking which flux is used in the signal channel for this rule */
  double eff1, eff2;
  if (glbGetChannelFlux(exp, glbGetChannelInRule(exp, rule, 0, GLB_SIG)) == 0) // nu
  {
    eff1 = wbb_params->eff_1st_max_nu;
    eff2 = wbb_params->eff_2nd_max_nu;
  }
  else
  {
    eff1 = wbb_params->eff_1st_max_nubar;
    eff2 = wbb_params->eff_2nd_max_nubar;
  }

  if (!(wbb_params->flags & WBB_NO_2ND_MAX))
  {
    if (wbb_params->flags & WBB_2ND_MAX_TOTAL_RATES)
      chi2 += wbb_wc_chi2_total_rates(ew_low, MIN(ew_high, i_1st_min), eff2);
    else
      chi2 += wbb_wc_sum_chi2(ew_low, MIN(ew_high, i_1st_min), eff2);
  }

  if (!(wbb_params->flags & WBB_NO_1ST_MAX))
  {
    if (wbb_params->flags & WBB_1ST_MAX_TOTAL_RATES)
      chi2 += wbb_wc_chi2_total_rates(MAX(ew_low, i_1st_min+1), ew_high, eff1);
    else
      chi2 += wbb_wc_sum_chi2(MAX(ew_low, i_1st_min+1), ew_high, eff1);
  }

  /* Systematics priors */
  for (i=0; i < n_params; i++)
    chi2 += square(x[i] / errors[i]);

  return chi2;
}


/***************************************************************************
 * Calculate chi^2 for WBB with LAr detector, including the following      *
 * systematical errors:                                                    *
 *                                                                         *
 *   x[0]: Correlated beam flux normalization error                        *
 *   x[1]: Correlated normalization of background in nu_e sample           *
 *   x[2]: Uncorr. signal norm for ND nu_e events                          *
 *   x[3]: Uncorr. signal norm for FD nu_e events                          *
 *   x[4]: Uncorr. BG norm errors for ND nu_e events                       *
 *   x[5]: Uncorr. BG norm errors for FD nu_e events                       *
 *   x[6]: Uncorr. signal norm errors for ND nu_mu events                  *
 *   x[7]: Uncorr. signal norm errors for FD nu_mu events                  *
 *   x[8]: Uncorr. BG norm errors for ND nu_mu events                      *
 *   x[9]: Uncorr. BG norm errors for FD nu_mu events                      *
 * and a fixed 1.0% bin-to-bin error in the near detector to avoid overly  *
 * optimistic NSI sensitivities                                            *
 *                                                                         *
 * IMPORTANT: The function assumed the ND experiment to be loaded directly *
 *            before the FD experiment, so that exp points to the FD,      *
 *            and exp-1 to the ND                                          *
 ***************************************************************************/
double chiWBB_LAr(int exp, int rule, int n_params, double *x, double *errors,
                  void *user_data)
{
  wbb_params_type *wbb_params = (wbb_params_type *) user_data;
  double *sig_true_rates_mu_N, *sig_true_rates_mu_F;
  double *bg_true_rates_mu_N, *bg_true_rates_mu_F;
  double *sig_true_rates_e_N, *sig_true_rates_e_F;
  double *bg_true_rates_e_N, *bg_true_rates_e_F;
  double *sig_fit_rates_mu_N, *sig_fit_rates_mu_F, *bg_fit_rates_mu_N, *bg_fit_rates_mu_F;
  double *sig_fit_rates_e_N, *sig_fit_rates_e_F, *bg_fit_rates_e_N, *bg_fit_rates_e_F;

  const int EXP_NEAR = exp - 1;
  const int EXP_FAR  = exp;
  const double sigma_binbin = 0.01;

  if (rule == RULE_WBB_LAR_E_CC)
  {
    sig_true_rates_mu_N  = glbGetSignalRatePtr(EXP_NEAR, RULE_WBB_LAR_MU_CC);
    sig_true_rates_mu_F  = glbGetSignalRatePtr(EXP_FAR,  RULE_WBB_LAR_MU_CC);
    bg_true_rates_mu_N   = glbGetBGRatePtr    (EXP_NEAR, RULE_WBB_LAR_MU_CC);
    bg_true_rates_mu_F   = glbGetBGRatePtr    (EXP_FAR,  RULE_WBB_LAR_MU_CC);
    sig_true_rates_e_N   = glbGetSignalRatePtr(EXP_NEAR, RULE_WBB_LAR_E_CC);
    sig_true_rates_e_F   = glbGetSignalRatePtr(EXP_FAR,  RULE_WBB_LAR_E_CC);
    bg_true_rates_e_N    = glbGetBGRatePtr    (EXP_NEAR, RULE_WBB_LAR_E_CC);
    bg_true_rates_e_F    = glbGetBGRatePtr    (EXP_FAR,  RULE_WBB_LAR_E_CC);

    sig_fit_rates_mu_N   = glbGetSignalFitRatePtr(EXP_NEAR, RULE_WBB_LAR_MU_CC);
    sig_fit_rates_mu_F   = glbGetSignalFitRatePtr(EXP_FAR,  RULE_WBB_LAR_MU_CC);
    bg_fit_rates_mu_N    = glbGetBGFitRatePtr    (EXP_NEAR, RULE_WBB_LAR_MU_CC);
    bg_fit_rates_mu_F    = glbGetBGFitRatePtr    (EXP_FAR,  RULE_WBB_LAR_MU_CC);
    sig_fit_rates_e_N    = glbGetSignalFitRatePtr(EXP_NEAR, RULE_WBB_LAR_E_CC);
    sig_fit_rates_e_F    = glbGetSignalFitRatePtr(EXP_FAR,  RULE_WBB_LAR_E_CC);
    bg_fit_rates_e_N     = glbGetBGFitRatePtr    (EXP_NEAR, RULE_WBB_LAR_E_CC);
    bg_fit_rates_e_F     = glbGetBGFitRatePtr    (EXP_FAR,  RULE_WBB_LAR_E_CC);
  }
  else if (rule == RULE_WBB_LAR_EBAR_CC)
  {
    sig_true_rates_mu_N  = glbGetSignalRatePtr(EXP_NEAR, RULE_WBB_LAR_MUBAR_CC);
    sig_true_rates_mu_F  = glbGetSignalRatePtr(EXP_FAR,  RULE_WBB_LAR_MUBAR_CC);
    bg_true_rates_mu_N   = glbGetBGRatePtr    (EXP_NEAR, RULE_WBB_LAR_MUBAR_CC);
    bg_true_rates_mu_F   = glbGetBGRatePtr    (EXP_FAR,  RULE_WBB_LAR_MUBAR_CC);
    sig_true_rates_e_N   = glbGetSignalRatePtr(EXP_NEAR, RULE_WBB_LAR_EBAR_CC);
    sig_true_rates_e_F   = glbGetSignalRatePtr(EXP_FAR,  RULE_WBB_LAR_EBAR_CC);
    bg_true_rates_e_N    = glbGetBGRatePtr    (EXP_NEAR, RULE_WBB_LAR_EBAR_CC);
    bg_true_rates_e_F    = glbGetBGRatePtr    (EXP_FAR,  RULE_WBB_LAR_EBAR_CC);
    
    sig_fit_rates_mu_N   = glbGetSignalFitRatePtr(EXP_NEAR, RULE_WBB_LAR_MUBAR_CC);
    sig_fit_rates_mu_F   = glbGetSignalFitRatePtr(EXP_FAR,  RULE_WBB_LAR_MUBAR_CC);
    bg_fit_rates_mu_N    = glbGetBGFitRatePtr    (EXP_NEAR, RULE_WBB_LAR_MUBAR_CC);
    bg_fit_rates_mu_F    = glbGetBGFitRatePtr    (EXP_FAR,  RULE_WBB_LAR_MUBAR_CC);
    sig_fit_rates_e_N    = glbGetSignalFitRatePtr(EXP_NEAR, RULE_WBB_LAR_EBAR_CC);
    sig_fit_rates_e_F    = glbGetSignalFitRatePtr(EXP_FAR,  RULE_WBB_LAR_EBAR_CC);
    bg_fit_rates_e_N     = glbGetBGFitRatePtr    (EXP_NEAR, RULE_WBB_LAR_EBAR_CC);
    bg_fit_rates_e_F     = glbGetBGFitRatePtr    (EXP_FAR,  RULE_WBB_LAR_EBAR_CC);
  }
  else
    return NAN;
  
  double sig_norm_mu_N, bg_norm_mu_N, sig_norm_mu_F, bg_norm_mu_F;
  double sig_norm_e_N, bg_norm_e_N, sig_norm_e_F, bg_norm_e_F;
  double *bin_centers;
  int ew_low, ew_high;

 
  /* Calculate nuisance parameters */
  sig_norm_e_N  = 1.0 + x[0] + x[2];
  sig_norm_e_F  = 1.0 + x[0] + x[3];
  bg_norm_e_N   = 1.0 + x[0] + x[1] + x[4];
  bg_norm_e_F   = 1.0 + x[0] + x[1] + x[5];
  sig_norm_mu_N = 1.0 + x[0] + x[6];
  sig_norm_mu_F = 1.0 + x[0] + x[7];
  bg_norm_mu_N  = 1.0 + x[0] + x[8];
  bg_norm_mu_F  = 1.0 + x[0] + x[9];

  /* Subroutine to compute chi2 using spectral information */
  double wbb_lar_sum_chi2(int ew_low, int ew_high, double eff_factor)
  {
    if (ew_low > ew_high)
      return 0.0;

    int i;
    double fit_rate, obs_rate;
    double eff_sig=eff_factor, eff_bg=1.0;
    double chi2 = 0.0;

    /* Loop over all bins in energy window */
    for (i=ew_low; i <= ew_high; i++)
    {
      /* Spectral analysis for nu_mu disappearance events */
      fit_rate = eff_sig*sig_norm_mu_N*sig_fit_rates_mu_N[i] + eff_bg*bg_norm_mu_N*bg_fit_rates_mu_N[i];
      obs_rate = eff_sig*sig_true_rates_mu_N[i] + eff_bg*bg_true_rates_mu_N[i];
      chi2 += gauss_likelihood( obs_rate, fit_rate,
                obs_rate * (1.0 + obs_rate*square(sigma_binbin)) );

      fit_rate = eff_sig*sig_norm_mu_F*sig_fit_rates_mu_F[i] + eff_bg*bg_norm_mu_F*bg_fit_rates_mu_F[i];
      obs_rate = eff_sig*sig_true_rates_mu_F[i] + eff_bg*bg_true_rates_mu_F[i];
      chi2 += poisson_likelihood(obs_rate, fit_rate);

      /* Spectral analysis for nu_e appearance events */
      fit_rate = eff_sig*sig_norm_e_N*sig_fit_rates_e_N[i] + eff_bg*bg_norm_e_N*bg_fit_rates_e_N[i];
      obs_rate = eff_sig*sig_true_rates_e_N[i] + eff_bg*bg_true_rates_e_N[i];
      chi2 += gauss_likelihood( obs_rate, fit_rate,
                obs_rate * (1.0 + obs_rate*square(sigma_binbin)) );

      fit_rate = eff_sig*sig_norm_e_F*sig_fit_rates_e_F[i] + eff_bg*bg_norm_e_F*bg_fit_rates_e_F[i];
      obs_rate = eff_sig*sig_true_rates_e_F[i] + eff_bg*bg_true_rates_e_F[i];
      chi2 += poisson_likelihood(obs_rate, fit_rate);
    }

    return chi2;
  }

  /* Subroutine to compute chi2 using total rates */
  double wbb_lar_chi2_total_rates(int ew_low, int ew_high, double eff_factor)
  {
    if (ew_low > ew_high)
      return 0.0;

    int i;
    double fit_rate, obs_rate;
    double eff_sig=eff_factor, eff_bg=1.0;
    double chi2 = 0.0;

    fit_rate = obs_rate = 0.0;
    for (i=ew_low; i <= ew_high; i++)
    {
      fit_rate += eff_sig*sig_norm_mu_N*sig_fit_rates_mu_N[i]
                    + eff_bg*bg_norm_mu_N*bg_fit_rates_mu_N[i];
      obs_rate += eff_sig*sig_true_rates_mu_N[i] + eff_bg*bg_true_rates_mu_N[i];
    }
    chi2 += gauss_likelihood( obs_rate, fit_rate,
                obs_rate * (1.0 + obs_rate*square(sigma_binbin)) );

    fit_rate = obs_rate = 0.0;
    for (i=ew_low; i <= ew_high; i++)
    {
      fit_rate += eff_sig*sig_norm_mu_F*sig_fit_rates_mu_F[i]
                    + eff_bg*bg_norm_mu_F*bg_fit_rates_mu_F[i];
      obs_rate += eff_sig*sig_true_rates_mu_F[i] + eff_bg*bg_true_rates_mu_F[i];
    }
    chi2 += poisson_likelihood(obs_rate, fit_rate);

    fit_rate = obs_rate = 0.0;
    for (i=ew_low; i <= ew_high; i++)
    {
      fit_rate += eff_sig*sig_norm_e_N*sig_fit_rates_e_N[i]
                    + eff_bg*bg_norm_e_N*bg_fit_rates_e_N[i];
      obs_rate += eff_sig*sig_true_rates_e_N[i] + eff_bg*bg_true_rates_e_N[i];
    }
    chi2 += gauss_likelihood( obs_rate, fit_rate,
              obs_rate * (1.0 + obs_rate*square(sigma_binbin)) );

    fit_rate = obs_rate = 0.0;
    for (i=ew_low; i <= ew_high; i++)
    {
      fit_rate += eff_sig*sig_norm_e_F*sig_fit_rates_e_F[i]
                    + eff_bg*bg_norm_e_F*bg_fit_rates_e_F[i];
      obs_rate += eff_sig*sig_true_rates_e_F[i] + eff_bg*bg_true_rates_e_F[i];
    }
    chi2 += poisson_likelihood(obs_rate, fit_rate);

    return chi2;
  }

  /* Request energy window */
  glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
  bin_centers = glbGetBinCentersListPtr(exp);

  /* Compute position of 1st oscillation minimum */
  double L = glbGetBaselineInExperiment(EXP_FAR);
  double E_1st_min = L*KM * fabs(true_ldm) / (4*M_PI) / GEV;
  int i_1st_min = 0;
  while (bin_centers[i_1st_min] <= E_1st_min)
    i_1st_min++;
  double chi2 = 0.0;
  int i;

  /* Find out if this chi^2 computation is for neutrino or anti-neutrino running   */
  /* We do this by checking which flux is used in the signal channel for this rule */
  double eff1, eff2;
  if (glbGetChannelFlux(exp, glbGetChannelInRule(exp, rule, 0, GLB_SIG)) == 0) // nu
  {
    eff1 = wbb_params->eff_1st_max_nu;
    eff2 = wbb_params->eff_2nd_max_nu;
  }
  else
  {
    eff1 = wbb_params->eff_1st_max_nubar;
    eff2 = wbb_params->eff_2nd_max_nubar;
  }

  if (!(wbb_params->flags & WBB_NO_2ND_MAX))
  {
    if (wbb_params->flags & WBB_2ND_MAX_TOTAL_RATES)
      chi2 += wbb_lar_chi2_total_rates(ew_low, MIN(ew_high, i_1st_min), eff2);
    else
      chi2 += wbb_lar_sum_chi2(ew_low, MIN(ew_high, i_1st_min), eff2);
  }

  if (!(wbb_params->flags & WBB_NO_1ST_MAX))
  {
    if (wbb_params->flags & WBB_1ST_MAX_TOTAL_RATES)
      chi2 += wbb_lar_chi2_total_rates(MAX(ew_low, i_1st_min), ew_high, eff1);
    else
      chi2 += wbb_lar_sum_chi2(MAX(ew_low, i_1st_min), ew_high, eff1);
  }

  /* Systematics priors */
  for (i=0; i < n_params; i++)
    chi2 += square(x[i] / errors[i]);

  return chi2;
}


/***************************************************************************
 * Calculate chi^2 for Double Chooz, including the following systematical  *
 * errors:                                                                 *
 *   x[0]: Flux normalization of reactor                                   *
 *   x[1]: Fiducial mass error - far detector                              *
 *   x[2]: Fiducial mass error - near detector                             *
 *   x[3]: Energy calibration error - far detector                         *
 *   x[4]: Energy calibration error - near detector                        *
 * and a fixed 0.5% bin-to-bin error                                       *
 ***************************************************************************/
double chiDCNorm(int exp, int rule, int n_params, double *x, double *errors,
                 void *user_data)
{ 
  int n_bins = glbGetNumberOfBins(EXP_REACTOR_FAR);
  const double sigma_binbin = 0.005;
  double *true_rates_N = glbGetRuleRatePtr(EXP_REACTOR_NEAR, 0);
  double *true_rates_F = glbGetRuleRatePtr(EXP_REACTOR_FAR, 0);
  double signal_fit_rates_N[n_bins];
  double signal_fit_rates_F[n_bins];
  double signal_norm_N, signal_norm_F;
  int ew_low, ew_high;
  double emin, emax;
  double fit_rate;
  double chi2 = 0.0;
  int i;

  /* Request simulated energy interval and analysis energy window */
  glbGetEminEmax(exp, &emin, &emax);
  glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
  
  /* Apply energy calibration error */
  glbShiftEnergyScale(x[3], glbGetSignalFitRatePtr(EXP_REACTOR_FAR, 0),
                      signal_fit_rates_F, n_bins, emin, emax);
  glbShiftEnergyScale(x[4], glbGetSignalFitRatePtr(EXP_REACTOR_NEAR, 0),
                      signal_fit_rates_N, n_bins, emin, emax);
  
  /* Loop over all bins in energy window */
  signal_norm_F = 1.0 + x[0] + x[1];
  signal_norm_N = 1.0 + x[0] + x[2];
  for (i=ew_low; i <= ew_high; i++)
  { 
    /* Statistical part of chi^2 for far detector
     * Normalization is affected by flux error x[0] and fiducial mass error x[1] */
    fit_rate = signal_norm_F * signal_fit_rates_F[i];
    chi2 += gauss_likelihood( true_rates_F[i], fit_rate,
              true_rates_F[i] * (1.0 + true_rates_F[i]*square(sigma_binbin)) );

    /* Statistical part of chi^2 for near detector
     * Normalization is affected by flux error x[0] and fiducial mass error x[2] */
    fit_rate = signal_norm_N * signal_fit_rates_N[i];
    chi2 += gauss_likelihood( true_rates_N[i], fit_rate,
              true_rates_N[i] * (1.0 + true_rates_N[i]*square(sigma_binbin)) );
  }
  
  /* Systematical part of chi^2 (= priors) */
  for (i=0; i < n_params; i++)
    chi2 += square(x[i] / errors[i]);

  /* Save the systematics parameters as starting values for the next step */
  for (i=0; i < n_params; i++)
    sys_startval_reactor[i] = x[i];

  return chi2;
}


// -------------------------------------------------------------------------
double chiKamLAND(int exp, int rule, int n_params, double *x, double *errors,
                  void *user_data)
// -------------------------------------------------------------------------
// Pedro Machado's chi^2 function for KamLAND analysis. Includes the
// following nuisance parameters:
//   x[0]: Overall normalization uncertainty
// -------------------------------------------------------------------------
{
  // Observed KamLAND rates, backgrounds and error bars
  // from http://kamland.lbl.gov/FiguresPlots/
  const double kamland_data[17] = { 97,  190, 177, 159, 121, 148,  158, 130, 143, 106,
                                    63,   47,  44,  13,   8,   2,    4 };
  const double kamland_bg[17]   = { 84.83755, 102.88808, 58.66426, 40.61372, 21.66065,
                                     4.51264,   0.00000,  0.90253,  4.51264,  2.70758,
                                     0.00000,   6.31769, 15.34296,  1.80505,  0.90253,
                                     0.90253,   0.90253 };
  const double kamland_sq_sigma[17] = {
    136.96692248317362,  346.63475597465106, 310.87724372053236, 265.76545960934311,
    182.98828914312122,  241.30339698428230, 264.20853879774438, 202.19410917350456, 
    229.43130196016733,  154.28832699937561,  80.064779671683027, 56.652328405465063,
     52.045806067816770,  12.579097447898485,  8.2514638571098473, 1.3378547794750271,
      3.3726788462200283  };

  double *signal_fit_rates[KAMLAND_N_REACT];
  double signal_norm;
  int ew_low, ew_high;
  double fit_rate;
  double chi2 = 0.0;
  int i;

  glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
  for(int j=0; j < KAMLAND_N_REACT; j++)
    signal_fit_rates[j] = glbGetSignalFitRatePtr(exp+j, 0);

  chi2 = 0;
  signal_norm = 1.0 + x[0];
  for (i=ew_low; i <= ew_high; i++)
  {
    fit_rate = 0.0;
    for(int j=0; j < KAMLAND_N_REACT; j++)
      fit_rate += signal_norm*signal_fit_rates[j][i];
    fit_rate += kamland_bg[i];
    chi2 += gauss_likelihood(kamland_data[i], fit_rate, kamland_sq_sigma[i]);
  }

  for (i=0; i < n_params; i++)
    chi2 += square(x[i] / errors[i]);

  return chi2;
}


// -------------------------------------------------------------------------
double chiLSNDspectrum(int exp, int rule, int n_params, double *x, double *errors,
                       void *user_data)
// -------------------------------------------------------------------------
// chi^2 function for the LSND experiment. Nuisance parameters are:
//   x[0]: Signal normalization (no prior!)
//   x[1]: Background normalization
// Author: Patrick Huber
// -------------------------------------------------------------------------
{
  // FIXME: JK - I think there are several oddities in this function, e.g.
  // why is the error multiplied by 2?
  // Also, check the glb file!
  return NAN;

//  const double error_bars[] = { 1.36364, 0.909091, 1.81818, 2.72727, 2.72727, 3.18182,
//                                3.40909, 3.40909,  3.40909, 3.86364, 2.27273 };
//  double *true_rates_N = glbGetRuleRatePtr(exp, rule);
//  double *signal,*bg;
//  double signal_norm_N;
//  int ew_low, ew_high;
//  double emin, emax;
//  double fit_rate;
//  double chi2 = 0.0;
//  int i;
//
//  /* Request simulated energy interval and analysis energy window */
//  glbGetEminEmax(exp, &emin, &emax);
//  glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
//
//  /* Apply energy calibration error 
//     glbShiftEnergyScale(x[3], glbGetSignalFitRatePtr(EXP_FAR, 0),
//     signal_fit_rates_F, n_bins, emin, emax);
//     glbShiftEnergyScale(x[4], glbGetSignalFitRatePtr(EXP_NEAR, 0),
//     signal_fit_rates_N, n_bins, emin, emax);
//  */
//  /* Loop over all bins in energy window */
//  signal=glbGetSignalFitRatePtr(exp, rule);
//  bg=glbGetBGFitRatePtr(exp, rule);
//  signal_norm_N = 1.0 + x[0];
//  for (i=ew_low; i <= ew_high; i++)
//  {
//    /* Statistical part of chi^2 for near detector
//     * Normalization is affected by flux error x[0] and fiducial mass error x[2] */
//    fit_rate  = signal_norm_N * signal[i] + x[1]*bg[i];
//    chi2 += gauss_likelihood(true_rates_N[i], fit_rate, 2.0*square(error_bars[i]));
//  }
//
//  /* Systematical part of chi^2 (= priors) */
//  for (i=1; i < n_params; i++)
//    chi2 += square(x[i] / errors[i]);
//
//  return chi2;
}


///***************************************************************************
// * Calculate chi^2 for T2K, including the following systematical errors:   *
// *   x[ 0]: Correlated beam flux normalization error                       *
// *   x[ 1]: Energy scale error for mu-like events                          *
// *   x[ 2]: Energy scale error for e-like events                           *
// *   x[ 3]: Spectral tilt for mu-like events                               *
// *   x[ 4]: Spectral tilt for e-like events                                *
// *   x[ 5]: Normalization mu-like events                                   *
// *   x[ 6]: Normalization e-like events                                    *
// * The energy scale errors are included even in SYS_OFF simulations to     *
// * help resolving degeneracies between systematics params and osc. params  *
// ***************************************************************************/
double chiT2K(int exp, int rule, int n_params, double *x, double *errors,
              void *user_data)
{
  double exp_F = exp;
  double exp_N = exp - 1;
  if (exp_N < 0)
    return NAN;

  // T2K far detector data from Run 1-4
  // from Lorena Escudero's thesis, http://www.t2k.org/docs/thesis/070, fig 5.15 (left)
  static const double data_mu_F[] = {
    0, 0, 0, 0, 0,  3, 3, 8, 6, 4,
    4, 2, 3, 1, 4,  1, 1, 4, 4, 2,
    2, 1, 4, 5, 2,  2, 3, 1, 0, 5,
    3, 0, 3, 0, 1,  0, 1, 1, 0, 0,
    1, 1, 0, 2, 1,  0, 1, 2, 0, 1,
    0, 0, 0, 1, 0,  3, 0, 1, 1, 2, 
    3, 1, 3, 2,
    3, 4,  1, 1,
    1
  };
  static const double data_e_F[] = {
    0, 0, 1, 0, 0,  0, 0, 2, 2, 2,
    3, 3, 3, 0, 1,  4, 2, 2, 1, 1,
    0, 0, 1, 0, 0
  };

  int n_bins_F = glbGetNumberOfBins(exp_F);
  int n_bins_N = glbGetNumberOfBins(exp_N);
  double emin_F, emax_F, emin_N, emax_N;
  int ew_low, ew_high;
  double fit_rate;
  double chi2 = 0.0;

  glbGetEminEmax(exp_F, &emin_F, &emax_F);
  glbGetEminEmax(exp_N, &emin_N, &emax_N);

  // Systematics ON
  // --------------
  if (n_params > 0)
  {
    double *bin_centers_F = glbGetBinCentersListPtr(exp_F);
    double E_center_F = 0.5 * (emax_F + emin_F);
    double signal_mu_F[n_bins_F], signal_e_F[n_bins_F];
    double bg_mu_F[n_bins_F], bg_e_F[n_bins_F];
    double signal_mu_N[n_bins_N], bg_mu_N[n_bins_N];
    double signal_mu_nosc_N[n_bins_N], bg_mu_nosc_N[n_bins_N];
    glbShiftEnergyScale(x[1], glbGetSignalFitRatePtr(exp_F, RULE_T2K_NUMU),
                        signal_mu_F, n_bins_F, emin_F, emax_F);
    glbShiftEnergyScale(x[1], glbGetBGFitRatePtr(exp_F, RULE_T2K_NUMU),
                        bg_mu_F, n_bins_F, emin_F, emax_F);
    glbShiftEnergyScale(x[2], glbGetSignalFitRatePtr(exp_F, RULE_T2K_NUE),
                        signal_e_F, n_bins_F, emin_F, emax_F);
    glbShiftEnergyScale(x[2], glbGetBGFitRatePtr(exp_F, RULE_T2K_NUE),
                        bg_e_F, n_bins_F, emin_F, emax_F);

    glbShiftEnergyScale(0.0, glbGetSignalFitRatePtr(exp_N, RULE_T2K_NUMU),
                        signal_mu_N, n_bins_N, emin_N, emax_N);
    glbShiftEnergyScale(0.0, glbGetBGFitRatePtr(exp_N, RULE_T2K_NUMU),
                        bg_mu_N, n_bins_N, emin_N, emax_N);
    glbShiftEnergyScale(0.0, glbGetSignalFitRatePtr(exp_N, RULE_T2K_NUMU_NOSC),
                        signal_mu_nosc_N, n_bins_N, emin_N, emax_N);
    glbShiftEnergyScale(0.0, glbGetBGFitRatePtr(exp_N, RULE_T2K_NUMU_NOSC),
                        bg_mu_nosc_N, n_bins_N, emin_N, emax_N);

    // Alter far detector prediction according to near detector measurement
    glbGetEnergyWindowBins(exp_F, RULE_T2K_NUMU, &ew_low, &ew_high);
    for (int i=ew_low; i <= ew_high; i++)
    {
      signal_mu_F[i] *= (signal_mu_nosc_N[i]+bg_mu_nosc_N[i]) / (signal_mu_N[i]+bg_mu_N[i]);
      signal_e_F[i]  *= (signal_mu_nosc_N[i]+bg_mu_nosc_N[i]) / (signal_mu_N[i]+bg_mu_N[i]);
    }

    // Compare FD data to ND-corrected prediction
    glbGetEnergyWindowBins(exp_F, RULE_T2K_NUMU, &ew_low, &ew_high);
    for (int i=ew_low; i <= ew_high; i++)
    {
      fit_rate = (1 + x[0] + x[5]) * ( signal_mu_F[i] + bg_mu_F[i]
                 + x[3] * (bin_centers_F[i] - E_center_F) * (signal_mu_F[i] + bg_mu_F[i]) );
      chi2    += poisson_likelihood(data_mu_F[i], fit_rate);
    }

    glbGetEnergyWindowBins(exp_F, RULE_T2K_NUE, &ew_low, &ew_high);
    for (int i=ew_low; i <= ew_high; i++)
    {
      fit_rate = (1 + x[0] + x[6]) * (signal_e_F[i] + bg_e_F[i]
                 + x[4] * (bin_centers_F[i] - E_center_F) * (signal_e_F[i] + bg_e_F[i]) );
      chi2    += poisson_likelihood(data_e_F[i], fit_rate);
    }

    // Systematics priors
    for (int i=0; i < n_params; i++)
      chi2 += square(x[i] / errors[i]);
  }

  // Systematics OFF
  // ---------------
  else
  {
    double *signal_mu_F      = glbGetSignalFitRatePtr(exp_F, RULE_T2K_NUMU);
    double *bg_mu_F          = glbGetBGFitRatePtr(exp_F, RULE_T2K_NUMU);
    double *signal_e_F       = glbGetSignalFitRatePtr(exp_F, RULE_T2K_NUE);
    double *bg_e_F           = glbGetBGFitRatePtr(exp_F, RULE_T2K_NUE);
    double *signal_mu_N      = glbGetSignalFitRatePtr(exp_N, RULE_T2K_NUMU);
    double *bg_mu_N          = glbGetBGFitRatePtr(exp_N, RULE_T2K_NUMU);
    double *signal_mu_nosc_N = glbGetSignalFitRatePtr(exp_N, RULE_T2K_NUMU_NOSC);
    double *bg_mu_nosc_N     = glbGetBGFitRatePtr(exp_N, RULE_T2K_NUMU_NOSC);

    // Alter far detector prediction according to near detector measurement
    glbGetEnergyWindowBins(exp_F, RULE_T2K_NUMU, &ew_low, &ew_high);
    for (int i=ew_low; i <= ew_high; i++)
    {
      signal_mu_F[i] *= (signal_mu_nosc_N[i]+bg_mu_nosc_N[i]) / (signal_mu_N[i]+bg_mu_N[i]);
      signal_e_F[i]  *= (signal_mu_nosc_N[i]+bg_mu_nosc_N[i]) / (signal_mu_N[i]+bg_mu_N[i]);
    }

    // Compare FD data to ND-corrected prediction
    glbGetEnergyWindowBins(exp_F, RULE_T2K_NUMU, &ew_low, &ew_high);
    for (int i=ew_low; i <= ew_high; i++)
    {
      fit_rate = signal_mu_F[i] + bg_mu_F[i];
      chi2    += poisson_likelihood(data_mu_F[i], fit_rate);
    }

    glbGetEnergyWindowBins(exp_F, RULE_T2K_NUE, &ew_low, &ew_high);
    for (int i=ew_low; i <= ew_high; i++)
    {
      fit_rate = signal_e_F[i] + bg_e_F[i];
      chi2    += poisson_likelihood(data_e_F[i], fit_rate);
    }
  }

  return chi2;
}



