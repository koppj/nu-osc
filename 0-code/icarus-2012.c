/**********************        ICARUS        *********************
 *
 * ICARUS simulation: bounds on sterile neutrinos
 *
 * Author: PAN Machado
 * Date:   2012-Oct
 *
 * Comments:
 * Based on 1209.0122v3
 *
 **********************        ICARUS        *********************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>   /* GLoBES library */


/* Density correlation list */
static int ICARUS_EXP=99;
static double NORMALIZATION=1E20;


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

double chi_ICARUS_2012(int exp, int rule, int n_params, double *x, double *errors,
		  void *user_data)
{
      double *signal_fit = glbGetSignalFitRatePtr(exp, 0);
      double chi2 = 0.0;
      double rate;

      rate = /* (1+x[0])* */(0.74*signal_fit[0]*NORMALIZATION + 3.7-1.0);
        /* nue bg = 3.7, of which 1.0 come from th13 */
      chi2 += glb_likelihood(2.0,rate); /* 2 nue events observed */

      /* the impact of the 5% normalization error is negligible */
      /* chi2 += glb_prior(x[0], 0.0, errors[0]); */
      return chi2;
}

int init_icarus_2012()
{
  /* Defining chi square */
  glbDefineChiFunction(&chi_ICARUS_2012,0,"chi_ICARUS_2012",NULL);
  ICARUS_EXP = glb_num_of_exps;
  glbInitExperiment("icarus-2012.glb",  &glb_experiment_list[0], &glb_num_of_exps);

  glb_params init_values = glbAllocParams();

  /* Defining GLoBES parameters */
  for (int i=0; i < glbGetNumOfOscParams(); i++)
    glbSetOscParams(init_values,0.0,i);
  glbSetDensityParams(init_values,1,GLB_ALL);

  /* The simulated data are computed */
//  glbSetCentralValues(init_values);
  glbSetOscillationParameters(init_values); 
  glbSetRates();

  double *signal_rates  = glbGetRuleRatePtr(ICARUS_EXP, 0);
  NORMALIZATION = 627.0/signal_rates[0]; /* they expect 627 numu CC */

  printf("# Initializing ICARUS code, Oct 2012 version ...\n");
  printf("#   Normalization: %g\n",NORMALIZATION);

  glbSetOscParams(init_values,0.157, 1);
  glbSetOscParams(init_values,M_PI/4.0, 2);
  glbSetOscParams(init_values,2.3E-3,5);

  /* glbDefineParams(init_values, */
  /* 		  asin(sqrt(0.023))/2., 0, 0., */
  /* 		  0,0.3,2.3e-3); */
  /* glbSetDensityParams(init_values,1,GLB_ALL); */

  /* The simulated data are computed */
//  glbSetCentralValues(init_values);
  glbSetOscillationParameters(init_values);
  glbSetRates();
  signal_rates  = glbGetRuleRatePtr(ICARUS_EXP, 0);
  printf("# Induced theta 13 events in ICARUS: %f (%f after cuts)\n",
   signal_rates[0]*NORMALIZATION-627,0.74*(signal_rates[0]*NORMALIZATION-627));
  printf("#\n");

  glbFreeParams(init_values);
  return 0;
}

