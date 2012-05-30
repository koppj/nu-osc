#include "def-reactors.h"

#define NBIN_CH_PROB 7

#define SIGMA  0.039 // 0.01 

#define EPR_MIN 0.8
#define EPR_MAX 6.4 
#define DE ( (EPR_MAX - EPR_MIN)/NBIN_CH_PROB )
#define E_RESOL 0.05


#define ACCURACY 1.e-5

struct coeff{
  double n0;
  double c0[NBIN_CH_PROB];
  double c0I[NBIN_CH_PROB][NISO];
  double c1[NBIN_CH_PROB];
};

coeff co[2];

extern Flux flux[NISO][2];  // old and new fluxes for each isotope
extern Fit fit;

 
double gauss_int(double eNu);
void calc_coef(void);

/********** global vars *******************/

double L_gl;
int cur_bin, iso_gl, old_new_gl;

const double isofract_chooz[NISO]={0.568, 0.078, 0.297, 0.057};

/**************** calc the table for class fit ***************/

void set_table_chooz(params &prm, double cff[NBIN_CHISQ][NPULLS+1])
{
  const double w = 1. - norm(prm.Ue3) - norm(prm.Ue[I4]) - norm(prm.Ue[I5]);
  const double w_prob = 1. - 2. * w * (norm(prm.Ue[I4]) + norm(prm.Ue[I5])) -
    2. * norm(prm.Ue3 * prm.Ue[I4]) - 
    2. * norm(prm.Ue3 * prm.Ue[I5]) - 
    2. * norm(prm.Ue[I4] * prm.Ue[I5]);

  for(int i = 0; i < NBIN_CH_PROB; i++){
    for(int r = 0; r < 2; r++){
      const int b = i + r*NBIN_CH_PROB;

      const double prob = w_prob - 4. * w * norm(prm.Ue3) * co[r].c1[i] / co[r].c0[i];

      // the prediction
      cff[b][NPULLS] = prob;

      // the pulls
      cff[b][PULL_U235] = co[r].c0I[i][U235] / co[r].c0[i];
      cff[b][PULL_U238] = co[r].c0I[i][U238] / co[r].c0[i];
      cff[b][PULL_P239] = co[r].c0I[i][P239] / co[r].c0[i];
      cff[b][PULL_P241] = co[r].c0I[i][P241] / co[r].c0[i];

      cff[b][NORM_CHOOZ] = cff[b][FLUX_NORM] = prob;
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
    
    fit.Data[i]              = input_data[i][1] / input_data[i][5];
    fit.Data[i+NBIN_CH_PROB] = input_data[i][3] / input_data[i][5];        

    fit.S_data[i][i] = norm(input_data[i][2] / input_data[i][5]);
    fit.S_data[i+NBIN_CH_PROB][i+NBIN_CH_PROB] = norm(input_data[i][4] / input_data[i][5]);

    fit.S_data[i][i+NBIN_CH_PROB] = fit.S_data[i+NBIN_CH_PROB][i] = 
      input_data[i][6] * 1.e-4 / norm(input_data[i][5]);
  }
  
  for(int i = 0; i < NBIN_CHOOZ; i++)
    fit.Data[i] *= rescale_data[i];
}

/********************************************************/
/*    calculation of coefficients                       */
/********************************************************/

  
double my_flux(const double Enu)
{
  if(iso_gl == -1){

    double w = 0.;
    for(int i = 0; i < NISO; i++)
      w += isofract_chooz[i] * flux[i][old_new_gl].f(Enu);
    return w;
  }

  return isofract_chooz[iso_gl] * flux[iso_gl][old_new_gl].f(Enu);
}


double fold0(double eNu)
{
  const double EposKin=eNu-ME-DELTA;
  return my_flux(eNu) * crossSect(EposKin)*gauss_int(eNu);
}

double fold2(double eNu)
{
  const double EposKin = eNu - ME - DELTA;
  const double arg_atm = 1.27e3 * DMQ_31 * L_gl / eNu;

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

	co[r].c0I[cur_bin][iso_gl] = qromb1(fold0, low, up, ACCURACY);
        co[r].c0[cur_bin] += co[r].c0I[cur_bin][iso_gl];
      }

      iso_gl = -1;
      co[r].c1[cur_bin] = qromb1(fold2, low, up, ACCURACY);
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


