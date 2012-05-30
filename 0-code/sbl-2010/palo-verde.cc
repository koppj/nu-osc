#include "def-reactors.h"

#define BIN_PV (NBIN_CHOOZ + N_SBL_R)

struct PV_coeff{
  double rate;
  double iso[NISO];
  double osc;
} PV_co;

extern Flux flux[NISO][2];  // old and new fluxes for each isotope
extern Fit fit;

extern int old_new_main;

void calc_coef_PV(void);

/********** global vars *******************/

int iso_gl_PV;

// estimated (in paper only 11% U238 is given)
const double isofract_PV[NISO] = {0.53, 0.11, 0.31, 0.05};


/**************** calc the table for class fit ***************/

#ifndef NO_PV
bool use_PV = true;
#else
bool use_PV = false;
#endif

void set_table_PV(params &prm, double cff[NBIN_CHISQ][NPULLS+1])
{  
  const double w = 1. - norm(prm.Ue3) - norm(prm.Ue[I4]) - norm(prm.Ue[I5]);
  const double prob = 1. - 2. * w * (norm(prm.Ue[I4]) + norm(prm.Ue[I5])) -
    4. * w * norm(prm.Ue3) * PV_co.osc / PV_co.rate -
    2. * norm(prm.Ue3 * prm.Ue[I4]) - 
    2. * norm(prm.Ue3 * prm.Ue[I5]) - 
    2. * norm(prm.Ue[I4] * prm.Ue[I5]);

  // the prediction
  cff[BIN_PV][NPULLS] = (use_PV ? prob : fit.Data[BIN_PV]);

  // the pulls
  cff[BIN_PV][PULL_U235] = (use_PV ? PV_co.iso[U235] / PV_co.rate : 0.);
  cff[BIN_PV][PULL_U238] = (use_PV ? PV_co.iso[U238] / PV_co.rate : 0.);
  cff[BIN_PV][PULL_P239] = (use_PV ? PV_co.iso[P239] / PV_co.rate : 0.);
  cff[BIN_PV][PULL_P241] = (use_PV ? PV_co.iso[P241] / PV_co.rate : 0.);

  cff[BIN_PV][FLUX_NORM] = (use_PV ? cff[BIN_PV][NPULLS] : 0.);
  return;
}

/***************** init *************************************/

void PV_init(void)
{
  // calc coefficients
  calc_coef_PV();
  
  if(old_new_main == OLD){
    
    fit.Data[BIN_PV] = 1.011;
    fit.S_data[BIN_PV][BIN_PV] = norm(0.053);
    
  }else{
    
    fit.Data[BIN_PV] = 0.975;
    fit.S_data[BIN_PV][BIN_PV] = norm(0.056);
  }
  return;  
}

/********************************************************/
/*    calculation of coefficients                       */
/********************************************************/

double PV_fold0(double eNu)
{
  const double EposKin=eNu-ME-DELTA;
  return isofract_PV[iso_gl_PV] * flux[iso_gl_PV][old_new_main].f(eNu) * crossSect(EposKin);
}

#define BL_PV 0.820

double PV_fold2(double eNu)
{
  const double EposKin = eNu - ME - DELTA;
  const double arg_atm = 1.27e3 * DMQ_31 * BL_PV / eNu;

  return isofract_PV[iso_gl_PV] * flux[iso_gl_PV][old_new_main].f(eNu) * 
    crossSect(EposKin) * norm(sin(arg_atm));
}


void calc_coef_PV(void)
{
  PV_co.rate = 0.;
  
  for(iso_gl_PV = 0; iso_gl_PV < NISO; iso_gl_PV++){
    
    PV_co.iso[iso_gl_PV] = qromb1(PV_fold0, EnuMIN, EnuMAX, 1.e-5);
    PV_co.rate += PV_co.iso[iso_gl_PV];
  }

  PV_co.osc = 0.;
    
  for(iso_gl_PV = 0; iso_gl_PV < NISO; iso_gl_PV++)
    PV_co.osc += qromb1(PV_fold2, EnuMIN, EnuMAX, 1.e-5);

  return;
}

