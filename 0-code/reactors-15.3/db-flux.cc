#include "definitions.h"

// Daya Bay 235 and 239 measurment 1704.01082

namespace ns_reactor
{

#ifdef USE_DB_FLUX

extern Fit fit;
extern int old_new_main;
extern Rate_pull_coupl rate_pull_coupl;  // defined in class_flux.cc
extern Rate_coef rate;                   // defined in class_flux.cc

#define NREACT_DB 6
#define NDET_DB 8
#define NDET_DB_FLUX 4 // using only detectors from EH1 and EH2 

double iso_DB_flux[NBIN_DB_FLUX][NISO];
double cov_sys_DBflux[NBIN_DB_FLUX][NBIN_DB_FLUX];  
  
// predicted effective cross sections, values from 1704.01082  
const double sigma_iso_pred[NISO] = {6.69, 10.1, 4.36, 6.05};

// baselines from table 2 of 1210.6327
const double basel_DB[NDET_DB][NREACT_DB] = {
  {0.362, 0.372, 0.903, 0.817, 1.354, 1.265},   // EH1
  {0.358, 0.368, 0.903, 0.817, 1.354, 1.266},   // EH1
  {1.332, 1.358, 0.468, 0.490, 0.558, 0.499},   // EH2
  {1.348, 1.348, 0.480, 0.480, 0.528, 0.528},   // EH2 inventing numbers!
  {1.920, 1.894, 1.533, 1.534, 1.551, 1.525},   // EH3
  {1.918, 1.892, 1.535, 1.535, 1.555, 1.528},   // EH3
  {1.925, 1.900, 1.539, 1.539, 1.556, 1.530},   // EH3
  {1.912, 1.912, 1.540, 1.540, 1.548, 1.548}    // EH3 inventing numbers!
};

const double DAQ_DB[NDET_DB] = {
  1117.178, 1117.178, 1114.337, 924.933, 1106.915, 1106.915, 1106.915, 917.417
};

// multiplying eps_mu * eps_m
const double eff_DB[NDET_DB] = {
  0.8255 * 0.9744,
  0.8221 * 0.9747,
  0.8573 * 0.9757,
  0.8571 * 0.9757,
  0.9824 * 0.9759,
  0.9823 * 0.9758,
  0.9821 * 0.9756,
  0.9826 * 0.9758
};

double weight_DB_flux[NDET_DB_FLUX][NREACT_DB];  
double weight_DB_flux_sum;
  
/**************** calc the table for class fit ***************/

void set_table_db_flux(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1])
{
  for(int k = 0; k < NBIN_DB_FLUX; k++){
    int b = fit.first_bin[DB_FLUX] + k;
  
    // probability weighted over reactors and detectors
    
    // as an approximation we take the F_i as f_i for
    // the calculation of the probability (see 1704.01082)

    // further approximation: we average over the near hall detectors
    // not clear how this is done in 1704.01082
    
    double pr = 0.;
    for(int i = 0; i < NDET_DB_FLUX; i++)
      for(int j = 0; j < NREACT_DB; j++)
	pr += weight_DB_flux[i][j] * rate.prob(prm, basel_DB[i][j] * 1000., iso_DB_flux[k]);
    pr /= weight_DB_flux_sum;

    // the prediction for the effective cross section   
    cff[b][NPULLS] = 0.;
    for(int i = 0; i < NISO; i++)
      cff[b][NPULLS] += iso_DB_flux[k][i] * sigma_iso_pred[i];   
    cff[b][NPULLS] *= pr;

    // the flux pulls
    cff[b][PULL_U235_0] = iso_DB_flux[k][U235] * sigma_iso_pred[U235] * pr;
    cff[b][PULL_U238]   = iso_DB_flux[k][U238] * sigma_iso_pred[U238] * pr;
    cff[b][PULL_P239_0] = iso_DB_flux[k][P239] * sigma_iso_pred[P239] * pr;
    cff[b][PULL_P241_0] = iso_DB_flux[k][P241] * sigma_iso_pred[P241] * pr;
    
    for(int i = 1; i < N_CO_UNC; i++){
      cff[b][PULL_U235_0 + i] = cff[b][NPULLS] * rate_pull_coupl.unc(U235, i, iso_DB_flux[k]);
      cff[b][PULL_P239_0 + i] = cff[b][NPULLS] * rate_pull_coupl.unc(P239, i, iso_DB_flux[k]); 
      cff[b][PULL_P241_0 + i] = cff[b][NPULLS] * rate_pull_coupl.unc(P241, i, iso_DB_flux[k]); 
    }

    cff[b][FLUX_COR]  = cff[b][NPULLS] * rate_pull_coupl.cor(iso_DB_flux[k]);
    cff[b][FLUX_NORM] = cff[b][NPULLS];

    // correlated detector efficiency
    cff[b][fit.first_pull[DB_FLUX]] = sqrt(cov_sys_DBflux[k][k]) * cff[b][NPULLS];
  }
  return;
}

/***************** init *************************************/

void init_db_flux(void)
{  
  if(old_new_main == OLD){
    fprintf(stderr, "cannot use DB_FLUX with old flux\n"); exit(0);
  }

  // read the data
  double data[NBIN_DB_FLUX];  
 
  // isotope fractions for the 8 bins from fig 1 of 1704.01082

  FILE *fp = fopen(REACTOR_PATH"Data_DB_flux/PRLsuppl/flux_data.dat","r");
  if(fp == NULL){
    fprintf(stderr, "ERROR init_db_flux: cannot open flux_data.dat\n"); exit(0);
  }
  for(int i = 0; i < NBIN_DB_FLUX; i++){
    if(fscanf(fp, "%*f %*f %lf %lf %lf %lf %lf\n",
	      &iso_DB_flux[i][P239], &iso_DB_flux[i][U235],
	      &iso_DB_flux[i][U238], &iso_DB_flux[i][P241], &data[i]) != 5){
      fprintf(stderr, "ERROR init_db_flux: cannot read flux_data.dat\n");exit(0);}
  }
  fclose(fp);

  // diagonal entries from Data_DB_flux/PRLsuppl/cov_stat.txt
  const double sigm_sq_stat[NBIN_DB_FLUX] = {
    1.22e-4, 1.37e-4, 1.34e-4, 1.37e-4, 1.14e-4, 1.11e-4, 1.10e-4, 1.97e-4
  };

  // systematic covariance matrix
  fp = fopen(REACTOR_PATH"Data_DB_flux/PRLsuppl/cov_syst.dat","r");
  if(fp == NULL){
    fprintf(stderr, "ERROR init_db_flux: cannot open cov_syst.dat\n"); exit(0);
  }
  for(int i = 0; i < NBIN_DB_FLUX; i++){
    for(int j = 0; j < NBIN_DB_FLUX; j++){
      if(fscanf(fp, "%lf", &cov_sys_DBflux[i][j]) != 1){
        fprintf(stderr, "ERROR init_db_flux: cannot read cov_syst.dat\n");exit(0);}
    }
  }
  fclose(fp);  
  // end reading data

  
 // calc the weights for averaged probability
  weight_DB_flux_sum = 0.;
  for(int i = 0; i < NDET_DB_FLUX; i++){
    for(int j = 0; j < NREACT_DB; j++){
      weight_DB_flux[i][j] = eff_DB[i] * DAQ_DB[i] / norm(basel_DB[i][j]);
      weight_DB_flux_sum += weight_DB_flux[i][j];
    }
  }
   
  // take into account theta13 effects
  Param_5nu p;
  p.theta[I12] = asin(sqrt(0.306));
  p.theta[I13] = asin(sqrt(0.0217));
  p.theta[I14] = p.theta[I15] = 0.;
  p.set_ang();
  p.dmq[0] = 0.;          // should always be zero!
  p.dmq[1] = 7.5e-5;
  p.dmq[2] = 2.51e-3;     
  p.dmq[3] = 1.;
  p.dmq[4] = 2.;

  // 3flavour probability weighted over reactors and detectors
  double pr[NBIN_DB_FLUX];
  for(int k = 0; k < NBIN_DB_FLUX; k++){
    pr[k] = 0.;
    for(int i = 0; i < NDET_DB_FLUX; i++)
      for(int j = 0; j < NREACT_DB; j++)
        pr[k] += weight_DB_flux[i][j] * rate.prob(p, basel_DB[i][j] * 1000., iso_DB_flux[k]);        
    pr[k] /= weight_DB_flux_sum;
  }

  const int b = fit.first_bin[DB_FLUX];
  
  // set the data points - rescaled by 3fl probability  
  for(int k = 0; k < NBIN_DB_FLUX; k++){
    fit.Data[b+k] = data[k] * pr[k];
    fit.S_data[b+k][b+k] = sigm_sq_stat[k] * norm(pr[k]);
  }

  
  // fully correlated part of the systematic covariance matrix is taken into account by pull
  // add the remaining part to the experimental error
  const double fudge = 0.5;
  // multiply by fudge factor 0.5 to match DB result of 7.9 for Delta chi^2 of global rescaling of fluxes
  
  for(int i = 0; i < NBIN_DB_FLUX; i++)
    for(int j = 0; j < NBIN_DB_FLUX; j++)
      fit.S_data[b+j][b+i] += data[j]*pr[j] * data[i]*pr[i] * fudge *
	(cov_sys_DBflux[i][j] - sqrt(cov_sys_DBflux[i][i] * cov_sys_DBflux[j][j]));
  /*
  fprintf(stderr, "\n");
  for(int i = 0; i < NBIN_DB_FLUX; i++){
    for(int j = 0; j < NBIN_DB_FLUX; j++)
      fprintf(stderr, "%e  ", cov_sys_DBflux[i][j] - sqrt(cov_sys_DBflux[i][i] * cov_sys_DBflux[j][j])); 
    fprintf(stderr, "\n");
  }
  exit(0);
  */
  
  // pull error
  const int pp = fit.first_pull[DB_FLUX]; 
  fit.S_pull[pp][pp] = 1.;
  return;  
}

#undef NREACT_DB 
#undef NDET_DB 
#undef NDET_DB_FLUX 

#endif // USE_DB_FLUX

}
