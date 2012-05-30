#include "def-reactors.h"

extern int old_new_main;

Fit fit;

double chi2reactor(params p)
{
  return fit.chisq(p);
}

void set_chisq_table(params &prm, double cff[NBIN_CHISQ][NPULLS+1])
{
  for(int i = 0; i < NBIN_CHISQ; i++)
    for(int j = 0; j <= NPULLS; j++)
      cff[i][j] = 0.;
#ifndef NO_CHOOZ
  set_table_chooz(prm, cff);
#endif  
  set_table_sbl(prm, cff);
  set_table_PV(prm, cff);
#ifndef NO_KAML
  set_table_kaml(prm, cff);
#endif
#ifndef BUGEY_TOTAL_RATE
  set_table_bugey(prm, cff);
#endif
  return;
}



/************/
/* Fit::Fit */
/************/

Fit::Fit(void)
{
   for(int i = 0; i < NPULLS; i++)
     pull_status[i] = ACTIVE;

   // overall scaling of fluxes in addition to individ. isotope uncert.
   pull_status[FLUX_NORM] = FIXED;

   // fix all the pulls for the isotope fractions of individ. react
   // (makes no difference in the result but speeds up things a lot)
#ifndef NO_KAML
   for(int i = PULL_KAML_0+1+KAML_NREACT; i < PULL_KAML_0+1+KAML_NREACT+KAML_NREACT*NISO; i++)
     pull_status[i] = FIXED;
#endif
  
   /* Initialize the correlation matrix for the pulls */

#ifdef BUGEY_TOTAL_RATE
   double pull_errors[PULL_KAML_0];
#else
   double pull_errors[PULL_KAML_0-1];
#endif

   // errors from tab I of Mention et al. 1101.xxx
   pull_errors[PULL_U235] = (old_new_main == OLD ? 0.019 : 0.0211);
   pull_errors[PULL_U238] = (old_new_main == OLD ? 0.1   : 0.0815);
   pull_errors[PULL_P239] = (old_new_main == OLD ? 0.024 : 0.0245);
   pull_errors[PULL_P241] = (old_new_main == OLD ? 0.021 : 0.0215);
  
   pull_errors[FLUX_NORM] = 0.02;   // invent just some number

   pull_errors[NORM_CHOOZ] = 0.01;  // adjusted to reproduce the old chooz bound

   /* set pull covariance matrix to 0 */
   for(int i = 0; i < NPULLS; i++)
     for(int j = 0; j < NPULLS; j++)
       S_pull[i][j] = 0.;

   // correlate U235, P239, P241 because of the common 1.8%  normalization
   // uncertainty from the Schreckenbach electron measurments
   S_pull[PULL_U235][PULL_P239] = S_pull[PULL_P239][PULL_U235] = norm(0.018);
   S_pull[PULL_U235][PULL_P241] = S_pull[PULL_P241][PULL_U235] = norm(0.018);
   S_pull[PULL_P239][PULL_P241] = S_pull[PULL_P241][PULL_P239] = norm(0.018);
  
   /* set pull errors apart from KamLAND and Bugey */
#ifdef BUGEY_TOTAL_RATE
   for(int i = 0; i < PULL_KAML_0; i++)
#else
   for(int i = 0; i < PULL_KAML_0-1; i++)
#endif
     S_pull[i][i] = norm(pull_errors[i]);

   /* set covariance matrix to 0 */
   for(int i = 0; i < NBIN_CHISQ; i++)
     for(int j = 0; j < NBIN_CHISQ; j++)
       S_data[i][j] = 0.;
}


/*****************/
/* Fit::invert_S */
/*****************/

// has to be called after S_data and S_pull have been set

void Fit::invert_S(void)
{
  if(invert(NBIN_CHISQ, S_data, S_data_inv))
    error("Fit.invert_S: cannot invert covariance matirx");

   if(invert(NPULLS, S_pull, S_pull_inv))
     error("Fit.invert_S: cannot invert pull covariance matrix");
  return;
}

/**************/
/* Fit::chisq */
/**************/

// if(use_xi_in) the fixed pull values from xi are used
// else fixed pulls are set to zero

double Fit::chisq(params &prm, double xi[NPULLS], bool use_xi_in)
{  
   double cff[NBIN_CHISQ][NPULLS+1];
   set_chisq_table(prm, cff);

   int n_active = 0;
   for(int a = 0; a < NPULLS; a++)
     if(pull_status[a] != FIXED) n_active++;
   
   // do nothing if all pulls fixed
   if(n_active > 0){

     double A[n_active][n_active+1];
     double pred[NBIN_CHISQ];

     for(int i = 0; i < NBIN_CHISQ; i++){
       pred[i] = cff[i][NPULLS];

       if(use_xi_in)
	 for(int a = 0; a < NPULLS; a++)
	   if(pull_status[a] == FIXED)
	     pred[i] += xi[a] * cff[i][a];
     }


     int aa = -1;
     for(int a = 0; a < NPULLS; a++){
       if(pull_status[a] != FIXED) {
	 aa++;

	 int bb = -1;
	 for(int b = 0; b < NPULLS; b++){
	   if(pull_status[b] != FIXED) {
	     bb++;

	     if(pull_status[a] == FREE || pull_status[b] == FREE)
	       A[aa][bb] = 0.;
	     else
	       A[aa][bb] = S_pull_inv[a][b];

	     for(int i = 0; i < NBIN_CHISQ; i++)
	       for(int j = 0; j < NBIN_CHISQ; j++)
		 A[aa][bb] += cff[i][a] * S_data_inv[i][j] * cff[j][b];
	   }
	 }
	 A[aa][n_active] = 0.;
	 for(int i = 0; i < NBIN_CHISQ; i++)
	   for(int j = 0; j < NBIN_CHISQ; j++)
	     A[aa][n_active] += cff[i][a] * S_data_inv[i][j] * (pred[j] - Data[j]);

       }
     }

     // solving the linear system
     double xxi[n_active];
     if(singsolve(n_active, A, xxi))
       error("Fit.chisq: cannot solve system");

     aa = -1;
     for(int a = 0; a < NPULLS; a++){

       if(pull_status[a] != FIXED){ 
	 aa++;
	 xi[a] = xxi[aa];
       }else 
	 if(!use_xi_in) xi[a] = 0.;
     }

   }else if(!use_xi_in) // all pulls fixed
     for(int a = 0; a < NPULLS; a++)
       xi[a] = 0.;

   /* build the chisq with the solution for the pulls */ 

   double t[NBIN_CHISQ];
   for(int i = 0; i < NBIN_CHISQ; i++){

     t[i] = cff[i][NPULLS];
     for(int a = 0; a < NPULLS; a++)
       t[i] += xi[a] * cff[i][a];
   }

   double chq = 0.;
   for(int i = 0; i < NBIN_CHISQ; i++)
     for(int j = 0; j < NBIN_CHISQ; j++)
       chq += (t[i] - Data[i]) * S_data_inv[i][j] * (t[j] - Data[j]);

   for(int a = 0; a < NPULLS; a ++)
     for(int b = 0; b < NPULLS; b++)
       if(pull_status[a] == ACTIVE && pull_status[b] == ACTIVE)
	 chq += xi[a] * S_pull_inv[a][b] * xi[b];

  
   /* print chisq terms 
   for(int i = 0; i < NBIN_CHISQ; i++){
     for(int j = 0; j < NBIN_CHISQ; j++)
       fprintf(stderr, "%+.3f  ", (t[i] - Data[i]) * S_data_inv[i][j] * (t[j] - Data[j]));
     fprintf(stderr, "\n");
   }
   double w = 0.;
   for(int a = 0; a < NPULLS; a ++)
     for(int b = 0; b < NPULLS; b++)
       if(pull_status[a] != FREE && pull_status[b] != FREE)
	 w += xi[a] * S_pull_inv[a][b] * xi[b];
   fprintf(stderr, "pull term: %f\n", w);  
  */
   return chq;
}
