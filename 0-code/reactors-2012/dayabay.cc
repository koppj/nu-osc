#include "definitions.h"

namespace ns_reactor
{

#ifdef USE_DB

// source: 1st Daya Bay paper + talk at Neutrino2012

#define NREACT_DB 6

extern Fit fit;
extern Rate_pull_coupl rate_pull_coupl;  // defined in class_flux.cc
extern Rate_coef rate;                   // defined in class_flux.cc
extern int old_new_main;

// assume same values as in DC
const double isofract_DB[NISO] = {0.488, 0.087, 0.359, 0.076};

/************** data **************************/

// baselines from table I
const double basel_DB[NBIN_DB][NREACT_DB] = {
  {0.364, 0.364, 0.857, 0.857, 1.307, 1.307},   // EH1
  {0.364, 0.364, 0.857, 0.857, 1.307, 1.307},   // EH1
  {1.348, 1.348, 0.480, 0.480, 0.528, 0.528},   // EH2
  {1.912, 1.912, 1.540, 1.540, 1.548, 1.548},   // EH3
  {1.912, 1.912, 1.540, 1.540, 1.548, 1.548},   // EH3
  {1.912, 1.912, 1.540, 1.540, 1.548, 1.548}    // EH3
};

// data from slide 13 of Nu2012 talk
const double IDB_cand[NBIN_DB] = {
  69121., 69714., 66473., 9788., 9669., 9452.
};
// days
const double DAQ[NBIN_DB] = {
  127.547, 127.547, 127.3763, 127.3763, 126.2646, 126.2646
};
const double eff_DB[NBIN_DB] = {
  0.8015, 0.7986, 0.8364, 0.9544, 0.9552, 0.9547
};
 
#define NBG_DB 5
// bg rate per day
const double bg_DB_data[NBG_DB][NBIN_DB] = {
  {9.73, 9.61, 7.55, 3.05, 3.04, 2.93},  // accidentals
  {0.77, 0.77, 0.58, 0.05, 0.05, 0.05},  // fast neutron
  {2.9, 2.9, 2.0, 0.22, 0.22, 0.22},     // Li, He
  {0.2,0.2,0.2,0.2,0.2,0.2},             // Am-C
  {0.08, 0.07, 0.05, 0.04, 0.04, 0.04}   // CO
};
// bg error per day
const double bg_err_DB[NBG_DB][NBIN_DB] = {
  {0.1, 0.1, 0.08, 0.04, 0.04, 0.03},
  {0.24, 0.24, 0.33, 0.02, 0.02, 0.02},
  {1.5, 1.5, 1.1, 0.12, 0.12, 0.12},
  {0.2, 0.2, 0.2, 0.2, 0.2, 0.2},
  {0.04, 0.04, 0.03, 0.02, 0.02, 0.02}
}; 
// slide 16: measured/expected reactor neutrinos
const double ratios_DB[NBIN_DB] = {
  0.981935, 0.987097, 0.984516, 0.943226, 0.929677, 0.918710
};

// number of bg events and relative error
double bg_DB[NBIN_DB][2];
double no_osc_react_DB[NBIN_DB];
double weight_DB[NBIN_DB][NREACT_DB]; 

/**************** calc the table for class fit ***************/

void set_table_DB(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1])
{  
  for(int i = 0; i < NBIN_DB; i++){

    const int b = fit.first_bin[DB] + i;
  
    cff[b][NPULLS] = 0.;

    for(int j = 0; j < NREACT_DB; j++){

      // reactor norm
      const int pr = fit.first_pull[DB] + 1 + j;
      cff[b][pr] = no_osc_react_DB[i] * weight_DB[i][j] * rate.prob(prm, basel_DB[i][j] * 1.e3, isofract_DB);
      cff[b][NPULLS] += cff[b][pr];
    }

    // overall normalization of reactor signal
    cff[b][fit.first_pull[DB]] = cff[b][NPULLS];

    // the fluxes
    cff[b][FLUX_NORM] = cff[b][NPULLS];

    // the pulls
    for(int j = 0; j < N_CO_UNC; j++){
      cff[b][PULL_U235_0 + j] = rate_pull_coupl.unc(U235, j, isofract_DB);
      cff[b][PULL_P239_0 + j] = rate_pull_coupl.unc(P239, j, isofract_DB); 
      cff[b][PULL_P241_0 + j] = rate_pull_coupl.unc(P241, j, isofract_DB); 
    }

    cff[b][PULL_U238] = rate_pull_coupl.unc(U238, 0, isofract_DB);
    cff[b][FLUX_COR]  = rate_pull_coupl.cor(isofract_DB);

    // uncorr detector normalizations
    const int pd = fit.first_pull[DB] + 1 + NREACT_DB + i;
    cff[b][pd] = cff[b][NPULLS];

    // adding background to prediction
    cff[b][NPULLS] += bg_DB[i][0];
        
    // the backgrounds
    const int pb = fit.first_pull[DB] + 1 + NREACT_DB + NBIN_DB + i;
    cff[b][pb] = bg_DB[i][0];
  }
  return;
}

/***************** init *************************************/

void DB_init(void)
{
  // set the background per detector
  for(int i = 0; i < NBIN_DB; i++){
    bg_DB[i][0] = bg_DB[i][1] = 0.;

    for(int j = 0; j < NBG_DB; j++){
      bg_DB[i][0] += bg_DB_data[j][i];
      bg_DB[i][1] += norm(bg_err_DB[j][i]);
    }
    bg_DB[i][1] = sqrt(bg_DB[i][1]) / bg_DB[i][0]; // relat. error

    bg_DB[i][0] *= DAQ[i] * eff_DB[i];
  }

  for(int i = 0; i < NBIN_DB; i++)
    no_osc_react_DB[i] = (IDB_cand[i] - bg_DB[i][0]) / ratios_DB[i];

  // set weight_DB: relat contribution of each react to each detector
  for(int i = 0; i < NBIN_DB; i++){

    double w = 0;
    for(int j = 0; j < NREACT_DB; j++)
      w += 1./norm(basel_DB[i][j]);

    for(int j = 0; j < NREACT_DB; j++)
      weight_DB[i][j] = 1./norm(basel_DB[i][j]) / w;
  }
  
  // set data and stat. errors in class fit
  for(int i = 0; i < NBIN_DB; i++){
  
    const int b = fit.first_bin[DB] + i;  
    fit.Data[b] = fit.S_data[b][b] = IDB_cand[i];    
  }

  /*** set pull errors in class fit ***/
  int p = fit.first_pull[DB];

  // normalization (free)
  fit.S_pull[p][p] = 100.; // just some number
  fit.pull_status[p] = FREE;

  // reactor norm
  for(int i = 0; i < NREACT_DB; i++){
    p++;
    fit.S_pull[p][p] = norm(0.008);
  }
  // detector norm
  for(int i = 0; i < NBIN_DB; i++){
    p++;
    fit.S_pull[p][p] = norm(0.002);
  }
  // backgrounds
  for(int i = 0; i < NBIN_DB; i++){
    p++;
    fit.S_pull[p][p] = norm(bg_DB[i][1]);
  }

  return;  
}

#endif

} // namespace
