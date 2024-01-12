#include "definitions.h"

namespace ns_reactor
{

extern int old_new_main;
bool fit_use_ex[N_EXP];

void set_chisq_table(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1])
{
  for(int i = 0; i < NBIN_CHISQ; i++)
    for(int j = 0; j <= NPULLS; j++)
      cff[i][j] = 0.;

#ifdef USE_SBL
  if(fit_use_ex[SBL])   set_table_sbl(prm, cff);
#endif
#ifdef USE_CHOOZ
  if(fit_use_ex[CHOOZ]) set_table_chooz(prm, cff);
#endif
#ifdef USE_PV
  if(fit_use_ex[PV])    set_table_PV(prm, cff);
#endif
#ifdef USE_KAML
  if(fit_use_ex[KAML])  set_table_kaml(prm, cff);
#endif
#ifdef USE_DC
  if(fit_use_ex[DC])    set_table_dc(prm, cff);
#endif
#ifdef USE_DB
  if(fit_use_ex[DB])    set_table_DB(prm, cff);
#endif
#ifdef USE_RENO
  if(fit_use_ex[RENO])  set_table_RENO(prm, cff);
#endif
#ifdef USE_BUGEY_SP
  if(fit_use_ex[BUG_SP])set_table_bugey(prm, cff);
#endif
#ifdef USE_DANSS
  if(fit_use_ex[DANSS]) set_table_danss(prm, cff);
#endif
#ifdef USE_GAL
  if(fit_use_ex[GAL])   set_table_gallium(prm, cff);
#endif
#ifdef USE_DB_FLUX
  if(fit_use_ex[DB_FLUX])set_table_db_flux(prm, cff);
#endif
  return;
}


/************/
/* Fit::Fit */
/************/

Fit::Fit(void)
{
   // setting number of bins for each experiment
   n_bin_ex[SBL]   = NBIN_SBL;
   n_bin_ex[CHOOZ] = NBIN_CHOOZ;
   n_bin_ex[PV]    = NBIN_PV;
   n_bin_ex[KAML]  = NBIN_KAML;
   n_bin_ex[DC]    = NBIN_DC;
   n_bin_ex[DB]    = NBIN_DB;
   n_bin_ex[RENO]  = NBIN_RENO;
   n_bin_ex[BUG_SP]= NBIN_BUG_SP;
   n_bin_ex[DANSS] = NBIN_DANSS;
   n_bin_ex[GAL]   = NBIN_GAL;
   n_bin_ex[DB_FLUX]=NBIN_DB_FLUX;
     
   // setting the first bin for each experiment
   first_bin[SBL]   = 0;
   for(int ex = 1; ex < N_EXP; ex++)
     first_bin[ex] = first_bin[ex-1] + n_bin_ex[ex-1];
    
   // setting the first pull for the private pulls for each experiment
   first_pull[SBL]   = PULL_GLOBAL;
   first_pull[CHOOZ] = first_pull[SBL]   + NPULL_SBL;
   first_pull[PV]    = first_pull[CHOOZ] + NPULL_CHOOZ;
   first_pull[KAML]  = first_pull[PV]    + NPULL_PV;
   first_pull[DC]    = first_pull[KAML]  + NPULL_KAML;
   first_pull[DB]    = first_pull[DC]    + NPULL_DC;
   first_pull[RENO]  = first_pull[DB]    + NPULL_DB;
   first_pull[BUG_SP]= first_pull[RENO]  + NPULL_RENO;
   first_pull[DANSS] = first_pull[BUG_SP]+ NPULL_BUG_SP;
   first_pull[GAL]   = first_pull[DANSS] + NPULL_DANSS;
   first_pull[DB_FLUX]=first_pull[GAL]   + NPULL_GAL;
   
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
   S_pull[FLUX_NORM][FLUX_NORM] = norm(0.030);   // invent just some number (not used -> either fixed or free)


#ifdef FLUX_FREE  
  S_pull[PULL_P241_0][PULL_P241_0] = norm(0.2);
  S_pull[PULL_U238][PULL_U238]   = norm(0.2);
  //pull_status[PULL_U238]   = pull_status[PULL_P241_0] = FREE;
  
  pull_status[PULL_U235_0] = pull_status[PULL_P239_0] = FREE;
#endif
  
   
   /* set data covariance matrix to 0 */
   for(int i = 0; i < NBIN_CHISQ; i++)
     for(int j = 0; j < NBIN_CHISQ; j++)
       S_data[i][j] = 0.;
    
   // first initialization of S_data_inv
   // (to be undone by set_experiments)
   n_bin_active = NBIN_CHISQ;
   S_data_inv = new double*[n_bin_active];
   for(int i = 0; i < n_bin_active; i++) 
     S_data_inv[i] = new double[n_bin_active];
}


/*****************/
/* Fit::invert_S */
/*****************/

// has to be called after S_data and S_pull have been set

void Fit::invert_S(void)
{
   if(invert(NPULLS, S_pull, S_pull_inv))
     error("Fit.invert_S: cannot invert pull covariance matrix");
  return;
}

   
/************************/
/* Fit::set_experiments */
/************************/

void Fit::set_experiments(bool use_ex[N_EXP])
{
  if(use_ex[DB] && use_ex[DB_FLUX])
    fprintf(stderr, "WARNING [set_experiments]: you are using DB and DB_FLUX together\n"); 
  
  // free memory of S_data_inv;
  for(int i = 0; i < n_bin_active; i++)
    delete[] S_data_inv[i];
  delete[] S_data_inv;
    
  n_bin_active = 0;
  for(int ex = 0; ex < N_EXP; ex++)
  {
    fit_use_ex[ex] = use_ex[ex];
    for(int b = 0; b < n_bin_ex[ex]; b++){
      
      if(use_ex[ex]) n_bin_active++;
      use_bin[first_bin[ex]+b] = use_ex[ex];
    }
  }
  // if Bugey spectrum is used, the Bugey3 rate has to be removed from the SBL data
  if(use_ex[SBL] && use_ex[BUG_SP]){
    use_bin[first_bin[SBL] + BUGEY3_1] = 
    use_bin[first_bin[SBL] + BUGEY3_2] = 
    use_bin[first_bin[SBL] + BUGEY3_3] = false;
    n_bin_active -= 3;
  }
  
  // checking
  // fprintf(stderr, "number of used data points: %d\n", n_bin_active);
  int k=0;
  for(int i = 0; i < NBIN_CHISQ; i++)
    if(use_bin[i]) k++;
  if(k != n_bin_active){
    fprintf(stderr, "fit.set_experiments: error: %d %d\n", n_bin_active, k);
    exit(0);
  }
  
    
  // allocate the memory for S_data_inv
  S_data_inv = new double*[n_bin_active];
  for(int i = 0; i < n_bin_active; i++) 
    S_data_inv[i] = new double[n_bin_active];
  
  // construct the covariance matrix for active bins
  double S[n_bin_active][n_bin_active], Sinv[n_bin_active][n_bin_active];
      
  int a = -1;
  for(int i = 0; i < NBIN_CHISQ; i++){
    if(use_bin[i]){ 
      a++;
      Data_active[a] = Data[i];
    
      int b = -1;
      for(int j = 0; j <NBIN_CHISQ; j++){
        if(use_bin[j]){ 
	  b++;
	  S[a][b] = S_data[i][j];
	}
      }
    }
  }
        
  if(invert(n_bin_active, S, Sinv))
    error("Fit.set_experiments: cannot invert covariance matrix");

  for(int i = 0; i < n_bin_active; i++)
    for(int j = 0; j < n_bin_active; j++)
      S_data_inv[i][j] = Sinv[i][j];
  
  return;
}

   
/**************/
/* Fit::chisq */
/**************/

// if(use_xi_in) the fixed pull values from xi are used
// else fixed pulls are set to zero

double Fit::chisq(Param_5nu &prm, double xi[NPULLS], bool use_xi_in)
{  
   double cff_all[NBIN_CHISQ][NPULLS+1];
   set_chisq_table(prm, cff_all);
  
   double cff[n_bin_active][NPULLS+1];
   int ii = -1;
   for(int i = 0; i < NBIN_CHISQ; i++){
     if(use_bin[i]){
       ii++;
       
       for(int a = 0; a < NPULLS+1; a++)
         cff[ii][a] = cff_all[i][a];
     }
   }

   int n_active = 0;
   for(int a = 0; a < NPULLS; a++)
     if(pull_status[a] != FIXED) n_active++;
   
   // do nothing if all pulls fixed
   if(n_active > 0){

     double A[n_active][n_active+1];
     double pred[n_bin_active];

     for(int i = 0; i < n_bin_active; i++){
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

         double t[n_bin_active];
         for(int i = 0; i < n_bin_active; i++){
             t[i] = 0.;
	     
             for(int j = 0; j < n_bin_active; j++){
                 t[i] += S_data_inv[i][j] * cff[j][a];
	     }
         }
	 	 
	 int bb = -1;
	 for(int b = 0; b < NPULLS; b++){
	   if(pull_status[b] != FIXED) {
	     bb++;

             // use that A has to be symmetric
             if(bb >= aa){
	     
	       if(pull_status[a] == FREE || pull_status[b] == FREE)
	         A[aa][bb] = 0.;
	       else
	         A[aa][bb] = S_pull_inv[a][b];

	       for(int i = 0; i < n_bin_active; i++){
		   A[aa][bb] += cff[i][b] * t[i];
	       }
	     }else
	       A[aa][bb] = A[bb][aa];
	   }
	 }
	 A[aa][n_active] = 0.;
	 for(int j = 0; j < n_bin_active; j++)
           A[aa][n_active] += t[j] * (pred[j] - Data_active[j]);
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

#define BGT_FRAC_175     (0.34/8.52)
   /* make sure that contribution of 175 keV excited state to 
    * cross section is always positive  
    */
   if(fit_use_ex[GAL]){
     const int p = first_pull[GAL];
     if(pull_status[p] != FIXED && xi[p] < -BGT_FRAC_175) xi[p] = -BGT_FRAC_175;
   }
#undef BGT_FRAC_175

   /* build the chisq with the solution for the pulls */ 

   double t[n_bin_active];
   for(int i = 0; i < n_bin_active; i++){
       
       t[i] = cff[i][NPULLS];
       for(int a = 0; a < NPULLS; a++)
         t[i] += xi[a] * cff[i][a];
   }

   double chq = 0., w[n_bin_active];
   for(int i = 0; i < n_bin_active; i++){
     w[i] = 0.;
     for(int j = 0; j < n_bin_active; j++)     
       w[i] += (t[j] - Data_active[j]) * S_data_inv[i][j];
   }
   for(int j = 0; j < n_bin_active; j++)
     chq += w[j] * (t[j] - Data_active[j]);

   double v[NPULLS];
   for(int a = 0; a < NPULLS; a ++){
     v[a] = 0.;
     for(int b = 0; b < NPULLS; b++)
       if(pull_status[b] == ACTIVE)
	 v[a] += xi[b] * S_pull_inv[a][b];
   }
  
   for(int a = 0; a < NPULLS; a ++)
     if(pull_status[a] == ACTIVE)
       chq += xi[a] * v[a];
  
   return chq;
}

}
