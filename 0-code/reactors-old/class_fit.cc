#include "definitions.h"

namespace ns_reactor
{

extern int old_new_main;

void set_chisq_table(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1])
{
  for(int i = 0; i < NBIN_CHISQ; i++)
    for(int j = 0; j <= NPULLS; j++)
      cff[i][j] = 0.;

#ifdef USE_SBL
  set_table_sbl(prm, cff);
#endif  
#ifdef USE_CHOOZ
  set_table_chooz(prm, cff);
#endif
#ifdef USE_PV  
  set_table_PV(prm, cff);
#endif  
#ifdef USE_KAML
  set_table_kaml(prm, cff);
#endif
#ifdef USE_DC
  set_table_dc(prm, cff);
#endif
#ifdef USE_DB
  set_table_DB(prm, cff);
#endif
#ifdef USE_RENO
  set_table_RENO(prm, cff);
#endif
#ifdef USE_BUGEY_SP
  set_table_bugey(prm, cff);
#endif
  return;
}


/************/
/* Fit::Fit */
/************/

Fit::Fit(void)
{
   // setting the first bin for each experiment
   first_bin[SBL]   = 0;
   first_bin[CHOOZ] = NBIN_SBL;
   first_bin[PV]    = NBIN_SBL + NBIN_CHOOZ;
   first_bin[KAML]  = NBIN_SBL + NBIN_CHOOZ + NBIN_PV;
   first_bin[DC]    = NBIN_SBL + NBIN_CHOOZ + NBIN_PV + NBIN_KAML;
   first_bin[DB]    = NBIN_SBL + NBIN_CHOOZ + NBIN_PV + NBIN_KAML + NBIN_DC;
   first_bin[RENO]  = NBIN_SBL + NBIN_CHOOZ + NBIN_PV + NBIN_KAML + NBIN_DC + NBIN_DB;
   first_bin[BUG_SP]= NBIN_SBL + NBIN_CHOOZ + NBIN_PV + NBIN_KAML + NBIN_DC + NBIN_DB + NBIN_RENO;

   // setting the first pull for the private pulls for each experiment
   first_pull[SBL]   = PULL_GLOBAL;
   first_pull[CHOOZ] = PULL_GLOBAL + NPULL_SBL;
   first_pull[PV]    = PULL_GLOBAL + NPULL_SBL + NPULL_CHOOZ;
   first_pull[KAML]  = PULL_GLOBAL + NPULL_SBL + NPULL_CHOOZ + NPULL_PV;
   first_pull[DC]    = PULL_GLOBAL + NPULL_SBL + NPULL_CHOOZ + NPULL_PV + NPULL_KAML;
   first_pull[DB]    = PULL_GLOBAL + NPULL_SBL + NPULL_CHOOZ + NPULL_PV + NPULL_KAML + NPULL_DC;
   first_pull[RENO]  = PULL_GLOBAL + NPULL_SBL + NPULL_CHOOZ + NPULL_PV + NPULL_KAML + NPULL_DC + NPULL_DB;
   first_pull[BUG_SP]= PULL_GLOBAL + NPULL_SBL + NPULL_CHOOZ + NPULL_PV + NPULL_KAML + NPULL_DC + NPULL_DB + NPULL_RENO;

   for(int i = 0; i < NPULLS; i++)
     pull_status[i] = ACTIVE;

   // overall scaling of fluxes in addition to individ. isotope uncert.
   pull_status[FLUX_NORM] = FIXED;

   /*
   pull_status[PULL_U235_1] = FIXED;
   pull_status[PULL_U235_2] = FIXED;
   pull_status[PULL_P239_1] = FIXED;
   pull_status[PULL_P239_2] = FIXED;
   pull_status[PULL_P241_1] = FIXED;
   pull_status[PULL_P241_2] = FIXED;
   */
   // fix all the pulls for the isotope fractions of individ. react
   // (makes no difference in the result but speeds up things a lot)
#ifdef USE_KAML
   for(int i = first_pull[KAML]+1+KAML_NREACT; i < first_pull[KAML]+1+KAML_NREACT+KAML_NREACT*NISO; i++)
     pull_status[i] = FIXED;
#endif
  
   /* Initialize the correlation matrix for the pulls */

   /* set pull covariance matrix to 0 */
   for(int i = 0; i < NPULLS; i++)
     for(int j = 0; j < NPULLS; j++)
       S_pull[i][j] = 0.;

   /* set global pull errors 
    * private ones are set in the corresponding experiment initialization 
    */

   // pulls from uncorrelated errors of U235, P239, P241 
   // (2nd order polynomial fit to uncorr. errors from Patrick's table)
   const int n_e = 25;    // number of Enu bins
   double S_inv[N_CO_UNC][N_CO_UNC], S[N_CO_UNC][N_CO_UNC], delta[n_e];

   // U235
   FILE *fp = fopen(REACTOR_PATH"Dat/Patrick-U235-err_unc.dat", "r");
   for(int i = 0; i < n_e; i++)
     if(fscanf(fp, "%*f %lf", &delta[i]) != 1)
       error("Flux::Flux cannot read uncorr error file\n"); 
   fclose(fp);

   for(int i = 0; i < N_CO_UNC; i++){
     for(int j = 0; j < N_CO_UNC; j++){

       S_inv[i][j] = 0.;
       for(int k = 0; k < n_e; k++){
	 const double e = 2. + 0.25 * k;
	 S_inv[i][j] += pow(e, i) * pow(e, j) / norm(delta[k]);
       }
     }
   }
   if(invert(N_CO_UNC, S_inv, S))
     error("Fit::Fit: cannot invert covariance matrix");

   for(int i = 0; i < N_CO_UNC; i++)
     for(int j = 0; j < N_CO_UNC; j++)
       S_pull[PULL_U235_0 + i][PULL_U235_0 + j] = S[i][j];


   // P239
   fp = fopen(REACTOR_PATH"Dat/Patrick-Pu239-err_unc.dat", "r");
   for(int i = 0; i < n_e; i++)
     if(fscanf(fp, "%*f %lf", &delta[i]) != 1)
       error("Flux::Flux cannot read uncorr error file\n"); 
   fclose(fp);

   for(int i = 0; i < N_CO_UNC; i++){
     for(int j = 0; j < N_CO_UNC; j++){

       S_inv[i][j] = 0.;
       for(int k = 0; k < n_e; k++){
	 const double e = 2. + 0.25 * k;
	 S_inv[i][j] += pow(e, i) * pow(e, j) / norm(delta[k]);
       }
     }
   }
   if(invert(N_CO_UNC, S_inv, S))
     error("Fit::Fit: cannot invert covariance matrix");

   for(int i = 0; i < N_CO_UNC; i++)
     for(int j = 0; j < N_CO_UNC; j++)
       S_pull[PULL_P239_0 + i][PULL_P239_0 + j] = S[i][j];

   // P241
   fp = fopen(REACTOR_PATH"Dat/Patrick-Pu241-err_unc.dat", "r");
   for(int i = 0; i < n_e; i++)
     if(fscanf(fp, "%*f %lf", &delta[i]) != 1)
       error("Flux::Flux cannot read uncorr error file\n"); 
   fclose(fp);

   for(int i = 0; i < N_CO_UNC; i++){
     for(int j = 0; j < N_CO_UNC; j++){

       S_inv[i][j] = 0.;
       for(int k = 0; k < n_e; k++){
	 const double e = 2. + 0.25 * k;
	 S_inv[i][j] += pow(e, i) * pow(e, j) / norm(delta[k]);
       }
     }
   }
   if(invert(N_CO_UNC, S_inv, S))
     error("Fit::Fit: cannot invert covariance matrix");

   for(int i = 0; i < N_CO_UNC; i++)
     for(int j = 0; j < N_CO_UNC; j++)
       S_pull[PULL_P241_0 + i][PULL_P241_0 + j] = S[i][j];


   /* the remaining global pulls */
   S_pull[FLUX_COR] [FLUX_COR]  = 1.;            // error included in the pull couplings
   S_pull[PULL_U238][PULL_U238] = norm(0.0815);  // error from tab I of Mention et al. 1101.2755 
   S_pull[FLUX_NORM][FLUX_NORM] = norm(0.027);   // invent just some number (not used -> either fixed or free)


   /* set data covariance matrix to 0 */
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

double Fit::chisq(Param_5nu &prm, double xi[NPULLS], bool use_xi_in)
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

}
