/***************************************************************************
 * Functions for the OPERA 2013 analysis                                   *
 * Based on http://arxiv.org/abs/1303.3953                                 *
 ***************************************************************************
 * Author: Pedro Machado                                                   *
 * Data:   2016                                                            *
 ***************************************************************************/
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_min.h>
#include <globes/globes.h>
#include "nu.h"

/* GLOBAL VAR */
/* Density correlation list */
static int density_corr[32];
static int OPERA_EXP=99;
static double NORMALIZATION=1E20;
static int OPERA_bins=3;


static double glb_prior(double x, double center, double sigma)
{
      double tmp = (x - center)/sigma;
      return tmp*tmp;
}


static double glb_likelihood(double true_rate, double fit_rate)
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


double chi_OPERA(int exp, int rule, int n_params, double *x, double *errors,
                 void *user_data)
{
      double *signal_fit = glbGetSignalFitRatePtr(exp, 0);
      double chi2 = 0.0;
      double rate;

      rate=0;
      /* All read from slide 29 of Neutrino 2016 Opera talk */
      double opera_data[3]  = {1., 6., 6.};
      double opera_bg[3]    = {0.58022, 3.40993, 5.15097}; /* bg WITHOUT SM oscillation */
      double opera_bg_sm[3] = {0.783446, 4.04865, 5.76063}; /* bg WITH SM oscillation */
      for(int i=0; i<OPERA_bins; i++){
        rate  = signal_fit[i]*NORMALIZATION + opera_bg_sm[i];
        chi2 += glb_likelihood(opera_data[i],(1+x[0])*rate);
      }
      /* 15% normalization error  */
      chi2 += glb_prior(x[0], 0.0, errors[0]);
      return chi2;
}


int init_OPERA()
{
  /* Defining chi square */
  glbDefineChiFunction(&chi_OPERA,1,"chi_OPERA",NULL);
  OPERA_EXP = glb_num_of_exps;
  glbInitExperiment("opera.glb", &glb_experiment_list[0], &glb_num_of_exps);
  NORMALIZATION = 2.7e3*17.97/5.25;
//  printf("# Initializing OPERA code, May 2016 version ...\n");
//  printf("#   Normalization: %g\n",NORMALIZATION);

  return 0;
}


