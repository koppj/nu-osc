/**********************        MINIBOONE        *********************
 *
 * MINIBOONE simulation: bounds on sterile neutrinos
 * NEUTRINO MODE: 6.46E20 POT
 *
 * Author: PAN Machado
 * Date:   2012-Jun
 *
 * Comments:
 * Based on 
 * http://www-boone.fnal.gov/for_physicists/data_release/lowe/
 *
 **********************        MINIBOONE        *********************/

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

/* GLOBAL VAR */
/* Density correlation list */
int density_corr[32];

int chiMB_init(int numode, int nubarmode);
int chiMB_clear();
double chiMB_nu(int exp, int rule, int n_params, double *x, double *errors,
		void *user_data);
double chiMB(int exp, int rule, int n_params, double *x, double *errors,
		   void *user_data);

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

  /* Initializing MB */
  chiMB_init(1,1);

  /* Define standard oscillation parameters */
  /* Based on hep-ph/1001.4524 */
  double theta12 = asin(sqrt(1E-2))/2.;
  double theta13 = 0;
  double bigtheta13 = 0;
  double theta23 = 0;
  double deltacp = 0.;
  double sdm = 1E0;
  double ldm = 0;

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
  glbSetOscillationParameters(true_values);
  glbSetRates();

  /* Set starting values and input errors for all projections */
  glbDefineParams(input_errors, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );

  glbSetDensityParams(input_errors,0.05,GLB_ALL);
  glbSetCentralValues(true_values);
  glbSetInputErrors(input_errors);

  /* Setting projection: */
  glbDefineProjection(standard_params_proj,GLB_FREE,GLB_FIXED,GLB_FIXED,
  		      GLB_FIXED,GLB_FREE,GLB_FIXED);
  glbSetDensityProjectionFlag(standard_params_proj, GLB_FIXED, GLB_ALL);

  glbDefineProjection(prescan_proj,GLB_FIXED,GLB_FIXED,GLB_FIXED,
  		      GLB_FIXED,GLB_FIXED,GLB_FIXED);
  glbSetDensityProjectionFlag(prescan_proj, GLB_FIXED, GLB_ALL);

  /* Iteration over all values to be computed */
  double degchimin,chimin; int globalmin;

  /* exit(0); */

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

  glb_params bf_find=glbAllocParams();  
  glbDefineParams(bf_find,0.015, 0, 0, 0, 3.0, 0);
  glbSetDensityParams(bf_find,1,GLB_ALL);

  degfinder(bf_find, n_prescan_params, prescan_params, prescan_min, prescan_max,
	    prescan_steps, prescan_proj, standard_params_proj, &n_deg, deg_pos, deg_chi2,DEG_NO_IH);
  for(int kk=0; kk<=n_deg-1; kk++)
    { if(kk==0){degchimin=deg_chi2[kk]; globalmin=kk;}
      if(deg_chi2[kk]<degchimin){degchimin=deg_chi2[kk]; globalmin=kk;} }

  FILE* bffile = fopen("MBbf.dat","w");
  fprintf(bffile,"%f %f %f\n",
	  SQR(sin(2.*glbGetOscParams(deg_pos[globalmin],GLB_THETA_12))),
  	  glbGetOscParams(deg_pos[globalmin],GLB_DM_21),deg_chi2[globalmin]);
  fprintf(bffile,"# best fit glb parameters:\n");
  glbPrintParams(bffile,deg_pos[globalmin]);

  glbPrintParams(stdout,deg_pos[globalmin]);

  glbSetProjectionFlag(standard_params_proj, GLB_FIXED, GLB_THETA_12);
  glbSetProjectionFlag(standard_params_proj, GLB_FIXED, GLB_DM_21);
  glbDefineParams(bf_find, 0, 0, 0, 0, 0, 0);
  glbSetDensityParams(bf_find,1,GLB_ALL);

  degfinder(bf_find, n_prescan_params, prescan_params, prescan_min, prescan_max,
	    prescan_steps, prescan_proj, standard_params_proj, &n_deg, deg_pos, deg_chi2,DEG_NO_IH);
  for(int kk=0; kk<=n_deg-1; kk++)
    { if(kk==0){degchimin=deg_chi2[kk]; globalmin=kk;}
      if(deg_chi2[kk]<degchimin){degchimin=deg_chi2[kk]; globalmin=kk;} }
  fprintf(bffile,"\n\n");
  fprintf(bffile,"%f %f %f\n",
	  SQR(sin(2.*glbGetOscParams(deg_pos[globalmin],GLB_THETA_12))),
  	  glbGetOscParams(deg_pos[globalmin],GLB_DM_21),deg_chi2[globalmin]);
  fprintf(bffile,"# no oscillation glb parameters:\n");
  glbPrintParams(bffile,deg_pos[globalmin]);
  fclose(bffile);

  /* exit(0); */

  double ss2t;  /* sterile sin^2 2theta */
  double dm;    /* sterile dm^2 */
  double theta; /* sterile mixing angle */

  /******************************************************************
   * chi2 in the dm x ss2t plane
   ******************************************************************/
  glbSetProjectionFlag(standard_params_proj, GLB_FIXED, GLB_THETA_12);
  glbSetProjectionFlag(standard_params_proj, GLB_FIXED, GLB_DM_21);

  FILE* chi2file=fopen("MBtest.dat","w");
  for(ss2t=1E-4; ss2t <= 1.1E0; ss2t*=1.2)
    {
      theta=M_PI/4.;
      if(ss2t<0.999)theta = asin(sqrt(ss2t))/2.;
      glbSetOscParams(true_values,theta,GLB_THETA_12);
      for(dm=1E-2; dm <= 110; dm*=1.2)
  	{
  	  glbSetOscParams(true_values,dm,GLB_DM_21);
      
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
	  printf("ss2t = %f    dm = %f   chi2 = %f\n",
		 ss2t,dm,deg_chi2[globalmin]);
  	}
      fprintf(chi2file,"\n");
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
