/***************************************************************************
 * Functions for the ICARUS 2014 analysis                                  *
 * Based on http://arxiv.org/abs/1209.0122 and Christian Farnese's talk at *
 * Nu2014, see https://indico.fnal.gov/conferenceOtherViews.py?confId=8022 *
 ***************************************************************************
 * Author: J. Kopp, based on work by Pedro Machado                         *
 * Data:   2016                                                            *
 ***************************************************************************/
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_min.h>
#include <globes/globes.h>
#include "nu.h"

//#define FELDMAN_COUSINS_TEST

static const int CH_MU_E   = 0;  /* Channel numbers */
static const int CH_E_E    = 1;
static const int CH_MU_TAU = 2;

static int ICARUS_EXP=99;
static double norm_sig    = 1E20;
static double norm_e_bg   = 1E20;
static double norm_tau_bg = 1E20;
static const double eff = 0.74;
static const double eff_e_bg = 0.65;
    /* Efficiencies for signal and for the intrinsic \nu_e BG, see slide 5
     * of Nu2014 slides */


inline static double glb_prior(double x, double center, double sigma)
{
      double tmp = (x - center)/sigma;
      return tmp*tmp;
}


inline static double glb_likelihood(double true_rate, double fit_rate)
{
      double res;
      res = fit_rate - true_rate;
      if (true_rate > 0)
      {
            if (fit_rate <= 0.0)
	      res = 1e100;
            else
	      res += true_rate * log( true_rate / fit_rate );
      }
      else
            res = fabs(res);
      return 2.0 * res;
}


/***************************************************************************
 * Function chi_icarus_2014                                                *
 ***************************************************************************
 * The chi^2 function for the ICARUS 2014 analysis                         *
 ***************************************************************************/
double chi_icarus_2014(int exp, int rule, int n_params, double *x, double *errors,
		       void *user_data)
{
  double *signal_fit = glbGetChannelFitRatePtr(ICARUS_EXP, CH_MU_E,   GLB_POST);
  double *bg_e_fit   = glbGetChannelFitRatePtr(ICARUS_EXP, CH_E_E,    GLB_POST);
  double *bg_tau_fit = glbGetChannelFitRatePtr(ICARUS_EXP, CH_MU_TAU, GLB_POST);
  double fit_rate;

  fit_rate = (norm_sig    * eff * signal_fit[0]
            + norm_e_bg   * eff_e_bg * bg_e_fit[0] 
            + norm_tau_bg * eff * bg_tau_fit[0]);

  /* Use the following to compare with official ICARUS plot, not taking into
   * account BG oscillation */
#warning "Using ICARUS code in debug mode"
  fit_rate = eff * signal_fit[0] * norm_sig
           + eff_e_bg * 7.0
           + eff * 1.6;
  fit_rate = eff * signal_fit[0] * norm_sig //FIXME
           + eff * 2.9
           + eff_e_bg * 7.0
           + eff * 1.6;

      /* Expected \nu_e background (from Nu2014 slides):
       * 7.0 (beam contamination) + 2.9 (\theta_{13}) + 1.6 (\nu_\tau) */
//  printf("  fit %g\n", fit_rate); //FIXME

#ifdef FELDMAN_COUSINS_TEST
  double *signal_rates = glbGetChannelRatePtr(ICARUS_EXP, CH_MU_E,   GLB_POST);
  double *bg_e_rates   = glbGetChannelRatePtr(ICARUS_EXP, CH_E_E,    GLB_POST);
  double *bg_tau_rates = glbGetChannelRatePtr(ICARUS_EXP, CH_MU_TAU, GLB_POST);

  double true_rate = (eff * signal_rates[0] * norm_sig
                    + eff_e_bg * 7.0
                    + eff * 1.6);
  double chi2 = glb_likelihood(true_rate, fit_rate);
#else
  double chi2 = glb_likelihood(6, fit_rate); /* 6 nue events observed */
#endif

  /* We have very low event number here -> do exact Poisson statistics and
   * convert to chi^2 value afterwards ("Gaussify") */
//  double p = gsl_ran_poisson_pdf(6, rate);
//  double chi2 = gsl_cdf_chisq_Pinv(1.-p, 1);

  /* the impact of the 5% normalization error is negligible */
  /* chi2 += glb_prior(x[0], 0.0, errors[0]); */
  return chi2;
}


#ifdef FELDMAN_COUSINS_TEST
/***************************************************************************
 * Function icarus_2014_fc_callback                                        *
 ***************************************************************************
 * Callback function for GSL minimizer in Feldman Cousins method           *
 ***************************************************************************/
double icarus_2014_fc_callback(double logs22thmue_test, void *params)
{
  glb_params test_values = (glb_params) params;
  glbSetOscParamByName(test_values, POW10(logs22thmue_test), "s22thmue");
  return glbChiSys(test_values, ICARUS_EXP, GLB_ALL);
}


/***************************************************************************
 * Function icarus_2014_feldman_cousins                                    *
 ***************************************************************************
 * Carry out Feldman-Cousins analysis of ICARUS 2014 data                  *
 ***************************************************************************/
int icarus_2014_feldman_cousins()
{
  const unsigned long n_pseudo_exp = 500;
  const double dm41 = 1.0; /* eV^2 */
  const double logs22thmue_min = -4.;
  const double logs22thmue_max = -0.5;
  const int logs22thmue_steps  = 35;
  int gsl_status;

  gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
  gsl_min_fminimizer *m = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);

  glb_params true_values = glbAllocParams();
  glb_params test_values = glbAllocParams();
  for (int i=0; i < glbGetNumOfOscParams(); i++)
    glbSetOscParams(true_values, 0.0, i);
  glbSetOscParamByName(true_values, asin(sqrt(0.0242)), "TH13");
  glbSetOscParamByName(true_values, M_PI/4,             "TH23");
  glbSetOscParamByName(true_values, 2.47E-3,            "DM31");
  glbSetOscParamByName(true_values, dm41,               "DM41");
  glbSetOscParamByName(true_values, 0.9,                "TH24");
  glbSetDensityParams(true_values,  0, GLB_ALL);

  double *signal_rates = glbGetChannelRatePtr(ICARUS_EXP, CH_MU_E,   GLB_POST);
  double *bg_e_rates   = glbGetChannelRatePtr(ICARUS_EXP, CH_E_E,    GLB_POST);
  double *bg_tau_rates = glbGetChannelRatePtr(ICARUS_EXP, CH_MU_TAU, GLB_POST);

  // Loop over "true" values of \sin^2 2\theta_{\mu e}
  for (int i=0; i <= logs22thmue_steps; i++)
  {
    double s22thmue_true = POW10(logs22thmue_min
        + i * (logs22thmue_max - logs22thmue_min)/logs22thmue_steps);
    glbSetOscParamByName(true_values, s22thmue_true, "s22thmue");
    glbSetOscillationParameters(true_values); 
    glbSetRates();
    double signal = eff * norm_sig
                      * glbGetChannelRatePtr(ICARUS_EXP, CH_MU_E,   GLB_POST)[0];
    double bg_e   = eff_e_bg * norm_e_bg
                      * glbGetChannelRatePtr(ICARUS_EXP, CH_E_E,    GLB_POST)[0];
    double bg_tau = eff * norm_tau_bg
                      * glbGetChannelRatePtr(ICARUS_EXP, CH_MU_TAU, GLB_POST)[0];

    if (debug_level > 1)
    {
      printf("FC: log(s22thmue_true) = %g\n", log10(s22thmue_true));
      printf("FC: Predicted rates:\n");
      printf("FC:   signal = %g\n", signal);
      printf("FC:   bg_e   = %g\n", bg_e);
      printf("FC:   bg_tau = %g\n", bg_tau);
    }

    // Loop over pseudo-experiments
    for (int n=0; n < n_pseudo_exp; n++)
    {
      // Apply statistical jitter
      signal_rates[0] = gsl_ran_poisson(rng, signal) / (eff * norm_sig);
      bg_e_rates[0]   = gsl_ran_poisson(rng, bg_e) / (eff * norm_e_bg);
      bg_tau_rates[0] = gsl_ran_poisson(rng, bg_tau) / (eff * norm_tau_bg);

      // Perform fit
      glbCopyParams(true_values, test_values);
      gsl_function F = { &icarus_2014_fc_callback, test_values };
      gsl_status = gsl_min_fminimizer_set(m, &F, log10(s22thmue_true),
                                          logs22thmue_min, logs22thmue_max);
      int iter = 0;
      do
      {
        iter++;
        gsl_status = gsl_min_fminimizer_iterate(m);
        if (gsl_status != GSL_SUCCESS)
        {
          fprintf(stderr, "icarus_2014_feldman_cousins: Minimization failed at "
                          "log(s22thmue_true) = %g, pseudo-exp #%d\n.",
                          log10(s22thmue_true), n);
          return -1;
        }
      } while (gsl_min_fminimizer_x_upper(m) - gsl_min_fminimizer_x_lower(m) > 0.001
            && iter < 1000);

      printf("FC_data  %2d %4d %10.7g %10.7g\n", i, n, log10(s22thmue_true),
             gsl_min_fminimizer_x_minimum(m));
    }
  }

  // Now do the fit to the real data
  signal_rates[0] = 6.0 / (eff * norm_sig);
  bg_e_rates[0]   = 0.0;
  bg_tau_rates[0] = 0.0;
  glbCopyParams(true_values, test_values);
  gsl_function F = { &icarus_2014_fc_callback, test_values };
  gsl_status = gsl_min_fminimizer_set(m, &F, -3.0, logs22thmue_min, logs22thmue_max);
  int iter = 0;
  do
  {
    iter++;
    gsl_status = gsl_min_fminimizer_iterate(m);
    if (gsl_status != GSL_SUCCESS)
    {
      fprintf(stderr, "icarus_2014_feldman_cousins: fit to data failed.\n");
      return -1;
    }
  } while (gsl_min_fminimizer_x_upper(m) - gsl_min_fminimizer_x_lower(m) > 0.001
        && iter < 1000);
  printf("FC_bf    %10.7g\n", gsl_min_fminimizer_x_minimum(m));

  gsl_rng_free(rng);          rng = NULL;
  gsl_min_fminimizer_free(m); m   = NULL;
  exit(0);
  return 0;
}
#endif


/***************************************************************************
 * Function init_icarus_2014                                               *
 ***************************************************************************
 * The chi^2 function for the ICARUS 2014 analysis                         *
 ***************************************************************************/
int init_icarus_2014()
{
  /* Defining chi square function */
  glbDefineChiFunction(&chi_icarus_2014,0,"chi_icarus_2014",NULL);
  ICARUS_EXP = glb_num_of_exps;
  glbInitExperiment("icarus-2014.glb", &glb_experiment_list[0], &glb_num_of_exps);

  /* Simulate \theta_{13}-induced \nu_e events to determine normalization */
  glb_params init_values = glbAllocParams();
  for (int i=0; i < glbGetNumOfOscParams(); i++)
    glbSetOscParams(init_values,0.0,i);
  glbSetOscParamByName(init_values, asin(sqrt(0.0242)), "TH13");
  glbSetOscParamByName(init_values, M_PI/4,             "TH23");
  glbSetOscParamByName(init_values, 2.47E-3,            "DM31");
  glbSetDensityParams(init_values,  1, GLB_ALL);
  glbSetOscillationParameters(init_values); 
  glbSetRates();
  double *signal_rates = glbGetChannelRatePtr(ICARUS_EXP, CH_MU_E,   GLB_POST);
  double *bg_e_rates   = glbGetChannelRatePtr(ICARUS_EXP, CH_E_E,    GLB_POST);
  double *bg_tau_rates = glbGetChannelRatePtr(ICARUS_EXP, CH_MU_TAU, GLB_POST);
  norm_sig    = 2.9 / signal_rates[0]; //FIXME
  norm_e_bg   = 7.0 / bg_e_rates[0];
  norm_tau_bg = 1.6 / bg_tau_rates[0];
      /* Expected background before applying efficiencies (from Nu2014 slides):
       * 7.0 (beam contamination) + 2.9 (\theta_{13}) + 1.6 (\nu_\tau) */

  //FIXME
//  double *mu_rates = glbGetChannelRatePtr(ICARUS_EXP, 3, GLB_POST);
//  printf(" ICARUS \nu_\mu rate: %g\n", mu_rates[0] * norm_sig);


  printf("# Initializing ICARUS code, Neutrino 2014 version ...\n");
  printf("#   Signal      normalization: %g\n", norm_sig);
  printf("#   BG \\nu_e    normalization: %g\n", norm_e_bg);
  printf("#   BG \\nu_\\tau normalization: %g\n", norm_tau_bg);

#ifdef FELDMAN_COUSINS_TEST
  icarus_2014_feldman_cousins();
  getchar();
#endif

  glbFreeParams(init_values);
  return 0;
}

