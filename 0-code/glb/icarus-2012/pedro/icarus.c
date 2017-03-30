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
#include "myheader.h"             /* Non-standard interactions */
#include "snu.h"             /* Non-standard interactions */
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

/* GLOBAL VAR */
/* Density correlation list */
int density_corr[32];
int ICARUS_EXP=99;
double NORMALIZATION=1E20;


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

double chi_ICARUS(int exp, int rule, int n_params, double *x, double *errors,
		  void *user_data)
{
      double *signal_fit = glbGetSignalFitRatePtr(exp, 0);
      double chi2 = 0.0;
      double rate;
      int i;

      rate = /* (1+x[0])* */(0.74*signal_fit[0]*NORMALIZATION + 3.7); /* nue bg = 3.7 */
      chi2 += glb_likelihood(2.0,rate); /* 2 nue events observed */

      /* the impact of the 5% normalization error is negligible */
      /* chi2 += glb_prior(x[0], 0.0, errors[0]); */
      return chi2;
}

int init_icarus()
{
  /* Defining chi square */
  glbDefineChiFunction(&chi_ICARUS,0,"chi_ICARUS",NULL);
  ICARUS_EXP = glb_num_of_exps;
  glbInitExperiment("./files/icarus.glb",  &glb_experiment_list[0], &glb_num_of_exps);

  glb_params init_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_projection standard_params_proj=glbAllocProjection();
  glb_projection prescan_proj=glbAllocProjection();

  /* Defining GLoBES parameters */
  for (int i=0; i < glbGetNumOfOscParams(); i++)
    glbSetOscParams(init_values,0.0,i);
  glbDefineParams(init_values,0,0,0,0,0,0);
  glbSetDensityParams(init_values,1,GLB_ALL);

  /* The simulated data are computed */
  glbSetCentralValues(init_values);
  snu_set_oscillation_parameters(init_values, NULL);
  /* glbSetOscillationParameters(init_values); */
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
  glbSetCentralValues(init_values);
  /* glbSetOscillationParameters(init_values); */
  snu_set_oscillation_parameters(init_values,NULL);
  glbSetRates();
  signal_rates  = glbGetRuleRatePtr(ICARUS_EXP, 0);
  printf("Induced theta 13 events: %f (%f after cuts)\n",
   signal_rates[0]*NORMALIZATION-627,0.74*(signal_rates[0]*NORMALIZATION-627));
  /* exit(0); */

  glbFreeParams(init_values);
  return 0;
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

  /* TAKEN FROM JOACHIM'S NU.CC CODE */
  // Initialize and register non-standard probability engine. This has to be done
  // before any calls to glbAllocParams() or glbAllocProjections()
  int n_flavors=4;
  int default_rotation_order[][2] = { {3,4}, {2,4}, {1,4}, {2,3}, {1,3}, {1,2} };
  int default_phase_order[] = { -1,  1, -1, -1,  0,   2};
  snu_init_probability_engine(n_flavors, default_rotation_order, default_phase_order);
  glbRegisterProbabilityEngine(6*SQR(n_flavors)-n_flavors+3, &snu_probability_matrix,
			       &snu_set_oscillation_parameters, &snu_get_oscillation_parameters, NULL);

  for (int j=0; j < glbGetNumOfOscParams(); j++)
    glbSetParamName(snu_param_strings[j], j);
  /************************************/

  init_icarus();		/* initializing ICARUS experiment */

  /* exit(0); */

  /* Define standard oscillation parameters */
  /* Based on hep-ph/1001.4524 */
  double theta12 = 0;
  double theta13 = 0.157;
  double theta23 = 0;
  double deltacp = 0.;
  double sdm = 0;
  double ldm = 2.4e-3;

  /* Initialize parameter vector(s) */
  glb_params true_values = glbAllocParams();
  glb_params test_values = glbAllocParams();
  glb_params no_oscillation = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_projection standard_params_proj=glbAllocProjection();
  glb_projection prescan_proj=glbAllocProjection();

  /* Defining GLoBES parameters */
  for (int i=0; i < glbGetNumOfOscParams(); i++)
    glbSetOscParams(true_values,0.0,i);
  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
  glbSetDensityParams(true_values,1,GLB_ALL);

  for (int i=0; i < glbGetNumOfOscParams(); i++)
    glbSetOscParams(test_values,0.0,i);
  glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);
  glbSetDensityParams(test_values,1,GLB_ALL);

  /* The simulated data are computed */
  glbSetCentralValues(true_values);
  glbSetOscillationParameters(true_values);
  glbSetRates();

  /* Set starting values and input errors for all projections */
  for (int i=0; i < glbGetNumOfOscParams(); i++)
    glbSetOscParams(input_errors,0.0,i);
  glbDefineParams(input_errors, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
  glbSetDensityParams(input_errors,0.05,GLB_ALL);

  glbSetCentralValues(true_values);
  glbSetInputErrors(input_errors);

  /* Setting projection: */
  for (int i=0; i < glbGetNumOfOscParams(); i++)
    glbSetProjectionFlag(standard_params_proj,0.0,i);
  glbDefineProjection(standard_params_proj,GLB_FIXED,GLB_FIXED,GLB_FIXED,
  		      GLB_FIXED,GLB_FIXED,GLB_FIXED);
  glbSetDensityProjectionFlag(standard_params_proj, GLB_FIXED, GLB_ALL);

  for (int i=0; i < glbGetNumOfOscParams(); i++)
    glbSetProjectionFlag(prescan_proj,0.0,i);
  glbDefineProjection(prescan_proj,GLB_FIXED,GLB_FIXED,GLB_FIXED,
  		      GLB_FIXED,GLB_FIXED,GLB_FIXED);
  glbSetDensityProjectionFlag(prescan_proj, GLB_FIXED, GLB_ALL);

  /* printf("we have a total of %d parameters\n",glbGetNumOfOscParams()); */
  /* exit(0); */

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

  /* glb_params bf_find=glbAllocParams();   */
  /* glbDefineParams(bf_find,M_PI/4., 0, 0, 0, 4E-2, 0); */
  /* glbSetDensityParams(bf_find,1,GLB_ALL); */

  /* degfinder(bf_find, n_prescan_params, prescan_params, prescan_min, prescan_max, */
  /* 	    prescan_steps, prescan_proj, standard_params_proj, &n_deg, deg_pos, deg_chi2,DEG_NO_IH); */
  /* for(int kk=0; kk<=n_deg-1; kk++) */
  /*   { if(kk==0){degchimin=deg_chi2[kk]; globalmin=kk;} */
  /*     if(deg_chi2[kk]<degchimin){degchimin=deg_chi2[kk]; globalmin=kk;} } */

  /* FILE* bffile = fopen("MBbf.dat","w"); */
  /* fprintf(bffile,"%f %f %f\n", */
  /* 	  SQR(sin(2.*glbGetOscParams(deg_pos[globalmin],GLB_THETA_12))), */
  /* 	  glbGetOscParams(deg_pos[globalmin],GLB_DM_21),deg_chi2[globalmin]); */
  /* fprintf(bffile,"# best fit glb parameters:\n"); */
  /* glbPrintParams(bffile,deg_pos[globalmin]); */

  /* glbSetProjectionFlag(standard_params_proj, GLB_FIXED, GLB_THETA_12); */
  /* glbSetProjectionFlag(standard_params_proj, GLB_FIXED, GLB_DM_21); */
  /* glbDefineParams(bf_find, 0, 0, 0, 0, 0, 0); */
  /* glbSetDensityParams(bf_find,1,GLB_ALL); */

  /* degfinder(bf_find, n_prescan_params, prescan_params, prescan_min, prescan_max, */
  /* 	    prescan_steps, prescan_proj, standard_params_proj, &n_deg, deg_pos, deg_chi2,DEG_NO_IH); */
  /* for(int kk=0; kk<=n_deg-1; kk++) */
  /*   { if(kk==0){degchimin=deg_chi2[kk]; globalmin=kk;} */
  /*     if(deg_chi2[kk]<degchimin){degchimin=deg_chi2[kk]; globalmin=kk;} } */
  /* fprintf(bffile,"\n\n"); */
  /* fprintf(bffile,"%f %f %f\n", */
  /* 	  SQR(sin(2.*glbGetOscParams(deg_pos[globalmin],GLB_THETA_12))), */
  /* 	  glbGetOscParams(deg_pos[globalmin],GLB_DM_21),deg_chi2[globalmin]); */
  /* fprintf(bffile,"# no oscillation glb parameters:\n"); */
  /* glbPrintParams(bffile,deg_pos[globalmin]); */
  /* fclose(bffile); */

  /* exit(0); */

  double ss2t;  /* sterile sin^2 2theta */
  double dm;    /* sterile dm^2 */
  double theta; /* sterile mixing angle */

  /******************************************************************
   * chi2 in the dm x ss2t plane
   ******************************************************************/
  glbSetProjectionFlag(standard_params_proj, GLB_FIXED, GLB_THETA_12);
  glbSetProjectionFlag(standard_params_proj, GLB_FIXED, GLB_DM_21);

  FILE* chi2file=fopen("./output/icarus_chi.dat","w");
  for(ss2t=3E-4; ss2t <= 1.2; ss2t*=1.2)
    {
      theta=M_PI/4.;
      if(ss2t<0.999)theta = asin(sqrt(ss2t))/2.;
      /* glbSetOscParams(true_values,theta,GLB_THETA_12); */
      glbSetOscParams(true_values,theta,0);
      for(dm=1E-2; dm <= 120; dm*=1.2)
  	{
  	  /* glbSetOscParams(true_values,dm,GLB_DM_21); */
  	  glbSetOscParams(true_values,dm,4);
	  snu_set_oscillation_parameters(true_values,NULL);
	  glbSetRates();
      
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
  	  fprintf(chi2file,"%10.7e  %10.7e  %10.7f\n",ss2t,dm,deg_chi2[globalmin]);
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
