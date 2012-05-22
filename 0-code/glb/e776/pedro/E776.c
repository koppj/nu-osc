/**********************        E776        *********************
 *
 * E776 simulation: bounds on sterile neutrinos
 *
 * Author: PAN Machado
 * Date:   2012-mar
 *
 * Comments:
 * See reference Phys Rev Lett 68 274 (1992)
 *
 **********************        E776        *********************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>   /* GLoBES library */
#include "myheader.h"             /* Non-standard interactions */
#include "glb_probability.h"

#define M_PI 3.14159265358979323846264338327950288419716939937510

/* Options for degfinder (flags) */
#define DEG_NO_NH      0x01   /* Omit normal hierarchy fits                    */
#define DEG_NO_IH      0x02   /* Omit inverted hierarchy fits                  */
#define DEG_NO_SYS     0x04   /* Switch off systematics                        */
#define DEG_NO_CORR    0x10   /* Switch off correlations                       */
#define DEG_NO_DEG     0x20   /* Switch off all degeneracies                   */
#define DEG_NO_NSI_DEG 0x40   /* Switch off degeneracies in the NSI parameters */
#define MAX_DEG    100    /* Maximum number of degeneracies to expect */
#define MAX_PARAMS  20

#define E776BINS (14) /* matches E776.glb */
/* GLOBAL VAR */
/* Density correlation list */
int density_corr[32];

int degfinder(const glb_params base_values, const int n_prescan_params,
      const int *prescan_params, const double *prescan_min,
      const double *prescan_max, const int *prescan_steps,
      const glb_projection prescan_proj, const glb_projection fit_proj,
      int *n_deg, glb_params *deg_pos, double *deg_chi2, const long flags);

double glb_prior(double x, double center, double sigma)
{
      double tmp = (x - center)/sigma;
      return tmp*tmp;
}

double glb_likelihood(double true_rate, double fit_rate)
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

double chi_E776(int exp, int rule, int n_params, double *x,
		double *errors, void *user_data)
{
      int n_bins = glbGetNumberOfBins(exp);
      double *data       = (double *) user_data; /* for both modes */
      /* this is pointing to E776data[4][E776BINS], so attention when reading */
      double *signal_fit_rates_pos = glbGetSignalFitRatePtr(exp, 0);
      double *signal_fit_rates_neg = glbGetSignalFitRatePtr(exp, 1);
      double *bg_fit_rates_pos     = glbGetBGFitRatePtr(exp, 0);
      double *bg_fit_rates_neg     = glbGetBGFitRatePtr(exp, 1);
      double pi0_norm_pos,pi0_norm_neg, nue_beam_norm_pos, nue_beam_norm_neg;
      double fit_rate_pos, fit_rate_neg;
      int ew_low, ew_high;
      double chi2 = 0.0;
      int i;

      glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
      pi0_norm_pos  = 1.0 + x[0]; /* pi0 bkgd error + */
      pi0_norm_neg  = 1.0 + x[1]; /* pi0 bkgd error - */
      nue_beam_norm_pos = 1.0 + x[2]; /* beam error + */
      nue_beam_norm_neg = 1.0 + x[3]; /* beam error - (not correlated?) */
      for (i=2; i <= ew_high; i++) /* i=2 to avoid first 2 bins (Enu > 1 GeV) */
	                           /* I cannot reproduce these very well */
      {
            fit_rate_pos = signal_fit_rates_pos[i]
	      + pi0_norm_pos*data[2*E776BINS+i] /* read pi0 + polarity */
      	      + nue_beam_norm_pos*bg_fit_rates_pos[i];
            fit_rate_neg = signal_fit_rates_neg[i]
	      + pi0_norm_neg*data[3*E776BINS+i] /* read pi0 - polarity */
      	      + nue_beam_norm_neg*bg_fit_rates_neg[i];
           chi2 += glb_likelihood(data[0*E776BINS+i], fit_rate_pos);
           chi2 += glb_likelihood(data[1*E776BINS+i], fit_rate_neg); /* antinu */
      }

      /* I have 4 errors: pi0 background and beam nue for each polarity */
      for(i=0; i<4; i++)chi2 += glb_prior(x[i], 0.0, errors[i]);


      {
          double *signal_fit_rates_pos = glbGetSignalFitRatePtr(0, 0);
          double *signal_fit_rates_neg = glbGetSignalFitRatePtr(0, 1);
          double *bg_fit_rates_pos     = glbGetBGFitRatePtr(0, 0);
          double *bg_fit_rates_neg     = glbGetBGFitRatePtr(0, 1);
          double events=0,events2=0;
          double bk=0,bk2=0;
          for(int i=0; i<14; i++)
            {
              events += signal_fit_rates_pos[i] + bg_fit_rates_pos[i]
                + data[2*E776BINS+i];
              events2 += signal_fit_rates_neg[i] + bg_fit_rates_neg[i]
                + data[3*E776BINS+i];
              bk += bg_fit_rates_pos[i];
              bk2 += bg_fit_rates_neg[i];
            }
          printf("chi=%e  evt+ = %e  evt- = %e\nbkgd+ = %e  bkgd- = %e\n",chi2,events,events2,bk,bk2);
      }









      return chi2;
}

/* I am not using this chi2 */
double chi_E776_rates(int exp, int rule, int n_params, double *x,
		double *errors, void *user_data)
{
      int n_bins = glbGetNumberOfBins(exp);
      double *data       = (double *) user_data; /* for both modes */
      double *signal_fit_rates_pos = glbGetSignalFitRatePtr(exp, 0);
      double *signal_fit_rates_neg = glbGetSignalFitRatePtr(exp, 1);
      double *bg_fit_rates_pos     = glbGetBGFitRatePtr(exp, 0);
      double *bg_fit_rates_neg     = glbGetBGFitRatePtr(exp, 1);
      double pi0_norm_pos,pi0_norm_neg, nue_beam_norm_pos, nue_beam_norm_neg;
      int ew_low, ew_high;
      double fit_rate_pos, fit_rate_neg;
      double chi2 = 0.0;
      int i;
      fit_rate_pos = 0;
      fit_rate_neg = 0;
      double sum_data_pos = 0;
      double sum_data_neg = 0;

      glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
      pi0_norm_pos  = 1.0 + x[0]; /* pi0 bkgd error + */
      pi0_norm_neg  = 1.0 + x[1]; /* pi0 bkgd error - */
      nue_beam_norm_pos = 1.0 + x[2]; /* beam error + */
      nue_beam_norm_neg = 1.0 + x[3]; /* beam error - (not correlated?) */
      for (i=0; i <= ew_high; i++)
      {
            fit_rate_pos += signal_fit_rates_pos[i]
	      + /* pi0_norm_pos* */data[2*E776BINS+i] /* read pi0 + polarity */
      	      + /* nue_beam_norm_pos* */bg_fit_rates_pos[i];
            fit_rate_neg += signal_fit_rates_neg[i]
	      + /* pi0_norm_neg* */data[3*E776BINS+i] /* read pi0 - polarity */
      	      + /* nue_beam_norm_neg* */bg_fit_rates_neg[i];
	    sum_data_pos +=data[0*E776BINS+i];
	    sum_data_neg +=data[1*E776BINS+i];
      }

      chi2 += glb_likelihood(sum_data_pos, fit_rate_pos);
      chi2 += glb_likelihood(sum_data_neg, fit_rate_neg); /* antinu */
      for(i=0; i<4; i++)chi2 += glb_prior(x[i], 0.0, errors[i]);
      return chi2;
}


/***************************************************************************
 *                                                                         *
 *                            M A I N   P R O G R A M                      *
 *                                                                         *
 *                                                                         *
 ***************************************************************************/
int main(int argc, char *argv[])
{
/* Default density correlation list: Every experiment independent */
  for (int i=0; i < 32; i++)
    density_corr[i] = i;

  /* Initialize libglobes */
  glbInit(argv[0]);
  glbSelectMinimizer(GLB_MIN_POWELL);

  /* E776 data from PRL 68 274 (1992):
     observed events for + (0) and - (1) polarities,
     pi0 background for + (2) and - (3) polarities  */
  double E776data[4][E776BINS] =
   {{ 35, 45, 23, 6, 6, 4, 1, 5, 5, 3, 1, 2, 0,0 },   /* + polarity, 136 total */
    {  9, 18,  6, 3, 4, 3, 1, 1, 0, 1, 0, 1, 0,0 },   /* - polarity, 47 total */
    {23.1, 43.4, 15.0, 6.28, 3.14, 0.92, 2.77, 0,0,0,0,0,0,0 }, /* + pol, pi0 94.61 */
    {8.01, 25.8, 4.69, 2.85, 0,0,0,0,0,0,0,0,0,0 } };           /* - pol, pi0 41.35 */
  //{9.01, 26.8, 4.69, 2.85, 0,0,0,0,0,0,0,0,0,0 } };           /* - pol, pi0 43.35 */

  /* Defining chi square */
  glbDefineChiFunction(&chi_E776,4,"chi_E776",&E776data);
  /* Initialize E776 */
  glbInitExperiment("E776.glb",&glb_experiment_list[0],&glb_num_of_exps);

  /* Define standard oscillation parameters */
  /* Based on hep-ph/1001.4524 */
  double theta12 = 0;
  double theta13 = 0;
  double bigtheta13 = 0;
  double theta23 = 0;
  double deltacp = 0.;
  double sdm = 0;
  double ldm = 0;

  /* Register non-standard probability engine. This has to be done
   * before any calls to glbAllocParams() or glbAllocProjections() */
  /* glbRegisterProbabilityEngine(6,     /\* Number of parameters *\/ */
  /*                              &twoflv_probability_matrix, */
  /*                              &twoflv_set_oscillation_parameters, */
  /*                              &twoflv_get_oscillation_parameters, */
  /*                              NULL); */

  /* Initialize parameter vector(s) */
  glb_params true_values = glbAllocParams();
  glb_params test_values = glbAllocParams();
  glb_params no_oscillation = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_projection standard_params_proj=glbAllocProjection();
  glb_projection prescan_proj=glbAllocProjection();

  /* Defining GLoBES parameters */
  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
  glbSetDensityParams(true_values,1,GLB_ALL);
  glbDefineParams(test_values,theta12,bigtheta13,theta23,deltacp,sdm,ldm);
  glbSetDensityParams(test_values,1,GLB_ALL);

  /* The simulated data are computed */
  glbSetCentralValues(true_values);

  /* Set starting values and input errors for all projections */
  glbDefineParams(input_errors, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );

  glbSetDensityParams(input_errors,0.05,GLB_ALL);
  glbSetCentralValues(true_values);
  glbSetInputErrors(input_errors);

  /* Setting projection: */
  glbDefineProjection(standard_params_proj,GLB_FIXED,GLB_FIXED,GLB_FIXED,
  		      GLB_FIXED,GLB_FIXED,GLB_FIXED);
  glbSetDensityProjectionFlag(standard_params_proj, GLB_FIXED, GLB_ALL);

  glbDefineProjection(prescan_proj,GLB_FIXED,GLB_FIXED,GLB_FIXED,
  		      GLB_FIXED,GLB_FIXED,GLB_FIXED);
  glbSetDensityProjectionFlag(prescan_proj, GLB_FIXED, GLB_ALL);

  glbSetOscillationParameters(true_values);
  glbSetRates();

  /* Iteration over all values to be computed */
  double degchimin,chimin; int globalmin;

  /***************************************************************
   *   CHI^2 CALCULATION - PROJECTING AND FINDING DEGENERACIES   *
   ***************************************************************/
  /* Data structures for degfinder -- prescan over th13, delta and arg(eps) */
  const double logs22th13min = -4.0; /* Limits (log10) in sin^2 2\theta_13 */
  const double logs22th13max = -1.0;
  const int logs22th13steps  = (logs22th13max - logs22th13min) * 4;
  /* # steps in s22th13 */
  /* degfinder inputs*/
  int n_prescan_params           = 0; /* # of params to be prescanned(PS) */
  int prescan_params[MAX_PARAMS] = { GLB_THETA_12     }; /* Index of params to be PS */
  double prescan_min[MAX_PARAMS] = { logs22th13min    }; /* Limits (min) */
  double prescan_max[MAX_PARAMS] = { logs22th13max    }; /* Limits (max) */
  int prescan_steps[MAX_PARAMS]  = { logs22th13steps  }; /* # steps */
  /* degfinder outputs */
  int n_deg; /* # found degeneracies */
  double deg_chi2[MAX_DEG]; /* chi2 for each degeneracy */
  glb_params deg_pos[MAX_DEG]; /* Parameters set for each degeneracy */
  for (int i=0; i < MAX_DEG; i++)
    deg_pos[i] = glbAllocParams();

  double ss2t;  /* sterile sin^2 2theta */
  double dm;    /* sterile dm^2 */
  double theta; /* sterile mixing angle */

  /******************************************************************
   * chi2 in the dm x ss2t plane
   ******************************************************************/
  FILE* chi2file=fopen("E776chi.dat","w");

//  glbSetOscParams(true_values,0.1,GLB_THETA_12);
//  glbSetOscParams(true_values,1.0,GLB_DM_21);
//  glbSetOscillationParameters(true_values);
//  printf("*** %g\n", glbConstantDensityProbability(1, 1, +1, 1.0, 1.0, 0.0));
  for(ss2t=0.001; ss2t <= 1.1; ss2t*=1.1)
    {
      theta=M_PI/4.;
      if(ss2t<0.999)theta = asin(sqrt(ss2t))/2.;
      glbSetOscParams(true_values,theta,GLB_THETA_12);
      for(dm=pow(10.0,-1.8)/*5E-2*/; dm <= 3E2; dm*=1.2)
  	{
  	  glbSetOscParams(true_values,dm,GLB_DM_21);
      
      /* Finding candidates for the global minimum */
  	  degfinder(true_values, n_prescan_params, prescan_params, prescan_min, prescan_max,
  		    prescan_steps, prescan_proj, standard_params_proj, &n_deg, deg_pos, deg_chi2,
                    DEG_NO_IH);
      
  	  /* Getting global minimum and respective parameters */
  	  for(int kk=0; kk<=n_deg-1; kk++)
  	    {
  	      if(kk==0){degchimin=deg_chi2[kk]; globalmin=kk;}
  	      if(deg_chi2[kk]<degchimin){degchimin=deg_chi2[kk]; globalmin=kk;}
  	    }
	  
  	  /* Writing chi square and parameters to file */
  	  fprintf(chi2file,"%10.7e  %10.7e  %10.7f\n",ss2t,dm,deg_chi2[globalmin]);
	  
	  double *signal_fit_rates_pos = glbGetSignalFitRatePtr(0, 0);
	  double *signal_fit_rates_neg = glbGetSignalFitRatePtr(0, 1);
	  double *bg_fit_rates_pos     = glbGetBGFitRatePtr(0, 0);
	  double *bg_fit_rates_neg     = glbGetBGFitRatePtr(0, 1);
	  double events=0,events2=0;
	  double bk=0,bk2=0;
	  for(int i=0; i<14; i++)
	    {
	      events += signal_fit_rates_pos[i] + bg_fit_rates_pos[i]
		+ E776data[2][i];
	      events2 += signal_fit_rates_neg[i] + bg_fit_rates_neg[i] 
		+ E776data[3][i];
	      bk += bg_fit_rates_pos[i];
	      bk2 += bg_fit_rates_neg[i];
	    }
	  printf("ss13=%e  dm=%e  chi=%e  evt+ = %e  evt- = %e\nbkgd+ = %e  bkgd- = %e\n",ss2t,dm,deg_chi2[globalmin],events,events2,bk,bk2);
          getchar();
  	}
    }
      
  fclose(chi2file);

  glbFreeParams(true_values);
  glbFreeParams(test_values);
  glbFreeParams(no_oscillation);
  glbFreeParams(input_errors);
  for (int i=0; i < MAX_DEG; i++)
    glbFreeParams(deg_pos[i]);
  glbFreeProjection(standard_params_proj);
  glbFreeProjection(prescan_proj);
  
  exit(0);
}
