/**********************        KARMEN        *********************
 *
 * KARMEN simulation: Carbon cross section
 *
 * Author: PAN Machado
 * Date:   Apr-2012
 *
 * Comments:
 * See J. Reichenbacher thesis for details, especially chapters 3 and 4
 *
 **********************        KARMEN        *********************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>   /* GLoBES library */
#include "myheader.h"
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

/******* Double Chooz simulation variables ********/
#define MAX_SYS 100
static double KM_sys_errors[MAX_SYS];       /* Uncertainties of
					       systematics params */
static double KM_sys_startval[MAX_SYS];     /* Starting values for
					       systematics
					       minimizer */
/*************************************************/

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

/* this is used to reproduce their analysis */
double JR(double true_rate, double fit_rate, double bkgd, double err)
{
      double tmp = (true_rate -bkgd - fit_rate);
      double tmp2 = true_rate + SQR(66.1) + SQR(0.7) /* SQR(err*fit_rate) */; 
      return tmp/sqrt(tmp2);
}

double KMchiJR(int exp, int rule, int n_params, 
		      double *x, double *errors, void *user_data)
{
      double KM_data_tot = 860;
      /* data from J. Reichenbacher thesis */
      double *signal_fit_rates = glbGetSignalFitRatePtr(exp, 0);
      double signal_rates[26];
      double signal_norm;
      double chi2 = 0.0;
      int  i,ew_low, ew_high;

      glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
      signal_norm = 1.0 + x[0];
      for(i=ew_low; i<=ew_high; i++)
	signal_rates[i] = signal_fit_rates[i];

      chi2 = JR(KM_data_tot,signal_rates[0], 13.9, 0.0825);
      return chi2;
}

double KMchi(int exp, int rule, int n_params, 
		      double *x, double *errors, void *user_data)
{
      /* double KM_data[26] = { 14, 20, 18, 33, 24, 34, 40, 36, 45, 34,  */
      /* 			     49, 60, 54, 56, 51, 39, 49, 53, 32, 36, */
      /* 			     27, 17, 21, 10, 5, 3}; */
      double KM_data_tot = 860;
      /* data from J. Reichenbacher thesis */
      double *signal_fit_rates = glbGetSignalFitRatePtr(exp, 0);
      double signal_rates[26],signal_tot;
      double signal_norm;
      double chi2 = 0.0;
      int  i,ew_low, ew_high;

      glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
      signal_norm = 1.0 + x[0];

      signal_tot = 0;
      for(i=ew_low; i<=ew_high; i++)
	signal_tot += signal_fit_rates[i];
      signal_tot *= signal_norm;
      signal_tot += 13.9; /* background */

      chi2 = glb_likelihood(KM_data_tot,signal_tot);
      chi2 += glb_prior(x[0], 0.0, errors[0]);
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

  /* Register 2 flavor probability engine. This has to be done 
   * before any calls to glbAllocParams() or glbAllocProjections() */
  /* glbRegisterProbabilityEngine(6,       /\* # of parameters - actually 2 *\/ */
  /*                              &twoflv_probability_matrix, */
  /*                              &twoflv_set_oscillation_parameters, */
  /*                              &twoflv_get_oscillation_parameters, */
  /*                              NULL); */

  /* Initialize KARMEN */
  glbDefineChiFunction(&KMchiJR,   1,    "KMchi", NULL);
  glbInitExperiment("KARMEN-nuecarbon.glb",&glb_experiment_list[0],&glb_num_of_exps);

  /* Define standard oscillation parameters */
  /* Based on hep-ph/1001.4524 */
  double th = asin(sqrt(0.13))/2.;
  double dm = 40.;
  
  /* Initialize parameter vector(s) */
  glb_params true_values = glbAllocParams();
  glb_params test_values = glbAllocParams();
  glb_params no_oscillation = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_projection standard_params_proj=glbAllocProjection();
  glb_projection prescan_proj=glbAllocProjection();

  /* Defining GLoBES parameters */
  glbDefineParams(true_values, th, 0, 0, 0, dm, 0 );
  glbSetDensityParams(true_values,1,GLB_ALL);
  glbDefineParams(test_values, th, 0, 0, 0, dm, 0 );
  glbSetDensityParams(test_values,1,GLB_ALL);

  /* The simulated data are computed */
  glbSetCentralValues(true_values);

  glbSetOscillationParameters(test_values);
  glbSetRates();
  double *signal_fit_rates = glbGetRuleRatePtr(0,0);
  double signal = signal_fit_rates[0];
  printf("num of events = %f  ==> %f oscillated\n",signal,806.2-signal);

  /* Set starting values and input errors for all projections */  
  glbDefineParams(input_errors, 0, 0, 0, 0, 0, 0 );

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

  /*************************************************************************
   * CHI^2 - th13 x chi2
   *************************************************************************/
  FILE* chi2file=fopen("KM-chi.dat","w");
  for(dm=1E-1; dm<105.; dm*=1.05)
    {
      glbSetOscParams(true_values,dm,GLB_DM_21);
      double lastchi = 0;
      double cutchi  = 10;
      for(double ss=0; ss < 1.01; ss+=0.01)
	{
	  th=M_PI/4.;
	  if(ss<0.999)th = asin(sqrt(ss))/2.;
	  glbSetOscParams(true_values,th,GLB_THETA_12);

	  if(lastchi>cutchi && ss>0.3 && dm>5E-1){
	    fprintf(chi2file,"%10.7e  %10.7e  %10.7f\n",dm,ss,lastchi);
            printf("lastchi2 = %f\n",lastchi);
	    continue;    }
	  /* Finding candidates for the global minimum */
	  degfinder(true_values, n_prescan_params, prescan_params, prescan_min, prescan_max,
		    prescan_steps, prescan_proj, standard_params_proj, &n_deg, deg_pos, deg_chi2,DEG_NO_IH);
	  
	  /* Getting global minimum and respective parameters */
	  for(int kk=0; kk<=n_deg-1; kk++)
	    {
	      if(kk==0){degchimin=deg_chi2[kk]; globalmin=kk;}
	      if(deg_chi2[kk]<degchimin){degchimin=deg_chi2[kk]; globalmin=kk;}
	    }
	  
	  /* Writing chi square and parameters to file */	
	  fprintf(chi2file,"%10.7e  %10.7e  %10.7f\n",dm,ss,deg_chi2[globalmin]);
	  printf("loop dm=%e ss=%e chi=%e\n",dm,ss,deg_chi2[globalmin]);
	  if( deg_chi2[globalmin]>cutchi )lastchi=deg_chi2[globalmin];
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
