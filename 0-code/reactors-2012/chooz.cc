#include "definitions.h"

namespace ns_reactor
{

#define NBIN_CH_PROB 7

#define EPR_MIN 0.8
#define EPR_MAX 6.4 
#define DE ( (EPR_MAX - EPR_MIN)/NBIN_CH_PROB )
#define E_RESOL 0.05


#define ACCURACY 1.e-5

struct coeff{
  double n0;
  double c0[NBIN_CH_PROB];
  double c0I[NBIN_CH_PROB][NISO][N_CO_UNC];
  double c0_cor[NBIN_CH_PROB];
  double c1[NBIN_CH_PROB][DM2_NUM];
};

coeff co[2];

extern Fit fit;
 
double gauss_int(double eNu);
void calc_coef(void);

/********** global vars *******************/

double dm2_glob, L_gl;
int cur_bin, iso_gl, old_new_gl;

const double isofract_chooz[NISO]={0.568, 0.078, 0.297, 0.057};

/**************** calc the table for class fit ***************/

void set_table_chooz(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1])
{
  // for interpolation
  const double logDmq = log10(fabs(prm.dmq[2]));
  const int idmq = int( (logDmq - ATM_MIN) / D_ATM);
  const double x = (logDmq - idmq * D_ATM - ATM_MIN) / D_ATM;

  if(logDmq <  ATM_MIN) error("set_table_chooz::Dmq too small");
  if(logDmq >= ATM_MAX) error("set_table_chooz::Dmq too large");
  if(idmq < 0 || idmq >= DM2_NUM-1) error("set_table_chooz: index out of range");
  if(x < 0. || x > 1.) error("set_table_chooz: x out of range");

  // coeffs for the probability
  const double sq_eff = 4. * norm(prm.Ue[2]) * (norm(prm.Ue[0]) + norm(prm.Ue[1])); 
  double w = 0;
  for(int i = 0; i < 3; i++) w += norm(prm.Ue[i]);
  const double cons = 2. * (norm(prm.Ue[3]) + norm(prm.Ue[4])) * w + 2. * norm(prm.Ue[3]) * norm(prm.Ue[4]);

  // set the table
  for(int i = 0; i < NBIN_CH_PROB; i++){
    for(int r = 0; r < 2; r++){
      const int b = fit.first_bin[CHOOZ] + i + r*NBIN_CH_PROB;

      // the prediction
      double sinq = co[r].c1[i][idmq] + (co[r].c1[i][idmq+1] - co[r].c1[i][idmq]) * x;
      sinq /= co[r].c0[i];

      cff[b][NPULLS] = 1. - sq_eff * sinq - cons;

      // the pulls
      for(int k = 0; k  < N_CO_UNC; k++){
        cff[b][PULL_U235_0 + k] = co[r].c0I[i][U235][k] / co[r].c0[i];
        cff[b][PULL_P239_0 + k] = co[r].c0I[i][P239][k] / co[r].c0[i];
        cff[b][PULL_P241_0 + k] = co[r].c0I[i][P241][k] / co[r].c0[i];
      }
      cff[b][PULL_U238] = co[r].c0I[i][U238][0] / co[r].c0[i];
      cff[b][FLUX_COR] = co[r].c0_cor[i] / co[r].c0[i];
      cff[b][fit.first_pull[CHOOZ]] = cff[b][FLUX_NORM] = cff[b][NPULLS];
     }
   }  
  return;
}

/***************** init *************************************/



void chooz_init(const int old_new)
{
  // calc coefficients

  old_new_gl = OLD;
  calc_coef();

  double rescale_data[NBIN_CHOOZ];

  if(old_new == NEW){

    for(int r = 0; r < 2; r++)
      for(int i = 0; i < NBIN_CH_PROB; i++)     
	rescale_data[i + r*NBIN_CH_PROB] = co[r].c0[i];

    old_new_gl = NEW;
    calc_coef();

    for(int r = 0; r < 2; r++)
      for(int i = 0; i < NBIN_CH_PROB; i++)
	rescale_data[i + r*NBIN_CH_PROB] /= co[r].c0[i];

  }else
    for(int r = 0; r < 2; r++)
      for(int i = 0; i < NBIN_CH_PROB; i++)     
	rescale_data[i + r*NBIN_CH_PROB] = 1.;

  // set data and covariance matrix in class fit

  // table 8 of hep-ex/0301017
  double input_data[NBIN_CH_PROB][7] = {     
    {1.2,  0.151, 0.031, 0.176, 0.035, 0.172, -2.2},
    {2.0,  0.490, 0.039, 0.510, 0.047, 0.532, -1.5},
    {2.8,  0.656, 0.041, 0.610, 0.049, 0.632, -3.5},
    {3.6,  0.515, 0.036, 0.528, 0.044, 0.530, -3.3},
    {4.4,  0.412, 0.033, 0.408, 0.040, 0.379, -2.0},
    {5.2,  0.248, 0.030, 0.231, 0.034, 0.208, -0.7},
    {6.0,  0.102, 0.023, 0.085, 0.026, 0.101, -1.3},
  };
  
  
  for(int i = 0; i < NBIN_CH_PROB; i++){
    const int ii = fit.first_bin[CHOOZ] + i;
    
    fit.Data[ii]              = input_data[i][1] / input_data[i][5];
    fit.Data[ii+NBIN_CH_PROB] = input_data[i][3] / input_data[i][5];        

    fit.S_data[ii][ii] = norm(input_data[i][2] / input_data[i][5]);
    fit.S_data[ii+NBIN_CH_PROB][ii+NBIN_CH_PROB] = norm(input_data[i][4] / input_data[i][5]);

    fit.S_data[ii][ii+NBIN_CH_PROB] = fit.S_data[ii+NBIN_CH_PROB][ii] = 
      input_data[i][6] * 1.e-4 / norm(input_data[i][5]);
  }
  
  for(int i = 0; i < NBIN_CHOOZ; i++)
    fit.Data[fit.first_bin[CHOOZ] + i] *= rescale_data[i];
  
  
  // set the CHOOZ pull
  // value adjusted to reproduce the old chooz bound
  fit.S_pull[fit.first_pull[CHOOZ]][fit.first_pull[CHOOZ]] = norm(0.01); 
  return;
}

/********************************************************/
/*    calculation of coefficients                       */
/********************************************************/

  
double my_flux(const double Enu)
{
  if(iso_gl == -1){

    double w = 0.;
    for(int i = 0; i < NISO; i++)
      w += isofract_chooz[i] * global_flux(i, Enu, old_new_gl);
    return w;
  }
  return isofract_chooz[iso_gl] * global_flux(iso_gl, Enu, old_new_gl);
}


double fold0(double eNu){
  const double EposKin=eNu-ME-DELTA;
  return my_flux(eNu) * crossSect(EposKin)*gauss_int(eNu);
}
double fold1(double eNu){
  const double EposKin=eNu-ME-DELTA;
  return eNu * my_flux(eNu) * crossSect(EposKin)*gauss_int(eNu);
}
double fold2(double eNu){
  const double EposKin=eNu-ME-DELTA;
  return eNu * eNu * my_flux(eNu) * crossSect(EposKin)*gauss_int(eNu);
}

double fold_corr(double eNu){
  iso_gl = U235;
  double w = read(REACTOR_PATH"Dat/Patrick-U235-err_cor.dat", eNu) * my_flux(eNu);
  iso_gl = P239;
  w += read(REACTOR_PATH"Dat/Patrick-Pu239-err_cor.dat", eNu) * my_flux(eNu);
  iso_gl = P241;
  w += read(REACTOR_PATH"Dat/Patrick-Pu241-err_cor.dat", eNu) * my_flux(eNu);

  const double EposKin=eNu-ME-DELTA;
  return w * crossSect(EposKin)*gauss_int(eNu);
}

double fold_osc(double eNu)
{
  const double EposKin = eNu - ME - DELTA;
  const double arg_atm = 1.27e3 * dm2_glob * L_gl / eNu;
  return my_flux(eNu) * crossSect(EposKin) * gauss_int(eNu) * norm(sin(arg_atm));
}


#define NSIGM 4.

void calc_coef(void)
{
  const double bl[2] = {.9979, 1.1146};
  
  for(int r = 0; r < 2; r++){

    L_gl = bl[r];      
    double low,up,l,u;
  
    cur_bin = -1;
    iso_gl = -1;

    co[r].n0 = qromb1(fold0,EnuMIN,EnuMAX,ACCURACY);

    for(cur_bin = 0; cur_bin < NBIN_CH_PROB; cur_bin++){

      l = EPR_MIN + cur_bin*DE;
      u = l + DE;
      low = l + 0.8 - NSIGM*E_RESOL*sqrt(l);
      if(low<EnuMIN) low = EnuMIN;
      up = u + 0.8 + NSIGM*E_RESOL*sqrt(u);
      if(up>EnuMAX) up=EnuMAX;

      co[r].c0[cur_bin] = 0.;

      for(iso_gl = 0; iso_gl < NISO; iso_gl++){

	co[r].c0I[cur_bin][iso_gl][0] = qromb1(fold0, low, up, ACCURACY);
        co[r].c0[cur_bin] += co[r].c0I[cur_bin][iso_gl][0];

	if(iso_gl != U238 && old_new_gl == NEW){
 	  co[r].c0I[cur_bin][iso_gl][1] = qromb1(fold1, low, up, ACCURACY);
 	  co[r].c0I[cur_bin][iso_gl][2] = qromb1(fold2, low, up, ACCURACY);
	}
      }

      if(old_new_gl == NEW){

        co[r].c0_cor[cur_bin] = qromb1(fold_corr, low, up, ACCURACY);

        iso_gl = -1;
        Parameters p;
        for(int j = 0; j < DM2_NUM; j++){

	  p.idmq_31 = j;
          dm2_glob = p.dmq31();
          co[r].c1[cur_bin][j] = qromb1(fold_osc, low, up, ACCURACY);
	}
      }
    }
  }
  return;
}



/****************************************************/
/* SUBROUTINE gauss_int(T_min, T_max, T_e, integ_T) */
/****************************************************/
      
double gauss_int(double eNu)
{   
   double T_min, T_max;
   const double T_e = eNu + ME - DELTA;
   double sigma = E_RESOL*sqrt(T_e);

   if(cur_bin == -1){
     T_min = EPR_MIN;
     T_max = EPR_MIN + NBIN_CH_PROB*DE;
   }else{
     T_min = EPR_MIN + cur_bin*DE;
     T_max = T_min + DE;
   }

   double x_max = (T_max - T_e) / (sigma*M_SQRT2);
   double x_min = (T_min - T_e) / (sigma*M_SQRT2);

   double res = (erf(x_max) - erf(x_min)) / 2.;
   if(res < 0.) {
      fprintf(stderr, "BLOODY HELL!!\n");
      exit(1);
   }
   
   return res;
}

}
