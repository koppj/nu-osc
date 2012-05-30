#include "definitions.h"

extern Fit fit;
extern int old_new_main;
extern Rate_pull_coupl rate_pull_coupl;  // defined in class_flux.cc
extern Rate_coef rate;                   // defined in class_flux.cc

void calc_coef_PV(void);

// baseline
#define BL_PV 820.

// estimated (in paper only 11% U238 is given)
const double isofract_PV[NISO] = {0.53, 0.11, 0.31, 0.05};

/**************** calc the table for class fit ***************/

void set_table_PV(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1])
{  
  const int b = fit.first_bin[PV];
  
  // the prediction
  cff[b][NPULLS] = rate.prob(prm, BL_PV, isofract_PV); 

  // the pulls
  for(int i = 0; i < N_CO_UNC; i++){
    cff[b][PULL_U235_0 + i] = rate_pull_coupl.unc(U235, i, isofract_PV);
    cff[b][PULL_P239_0 + i] = rate_pull_coupl.unc(P239, i, isofract_PV); 
    cff[b][PULL_P241_0 + i] = rate_pull_coupl.unc(P241, i, isofract_PV); 
  }

  cff[b][PULL_U238] = rate_pull_coupl.unc(U238, 0, isofract_PV);
  cff[b][FLUX_COR]  = rate_pull_coupl.cor(isofract_PV);
  cff[b][FLUX_NORM] = cff[b][NPULLS];
  return;
}

/***************** init *************************************/

void PV_init(void)
{
  const int b = fit.first_bin[PV];
  
  if(old_new_main == OLD){
    
    fit.Data[b] = 1.011;
    fit.S_data[b][b] = norm(0.053);
    
  }else{
    
    fit.Data[b] = 0.975;
    fit.S_data[b][b] = norm(0.056);
  }
  return;  
}

