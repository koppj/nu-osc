#include "definitions.h"

/*****************************************************
 * using results presented at Cern in September 2016
 * https://indico.cern.ch/event/548805/
 *****************************************************/

namespace ns_reactor
{

#ifdef USE_DC

#define ACCURACY 1.e-5

/*** energy resolution from 1406.7763 ***/   

const double EresA = norm(0.077);
const double EresB = norm(0.018);
const double EresC = norm(0.017);

inline double e_res_dc(double e){
  return sqrt(EresA / e + EresB + EresC / norm(e)) * e;
}

extern Fit fit;
 
/********** global vars *******************/

double gl_DmqL_dc;
int bin_gl_dc, iso_gl_dc, old_new_gl_dc;

double sys_DC_e[NBIN_DC];   // energy scale uncert (blue histogram)
double sys_DC[NBIN_DC];     // reduced detector uncert (green histogram)
  
const double isofract_dc[NISO]={0.488, 0.087, 0.359, 0.076};

/* binning in positron energy */
static double bin[NBIN_DC/2][2];

const double bl_dc[2] = {.4, 1.05};  

enum DC_DET {DC_NEAR, DC_FAR};   
  
/**************** class probability ***************/

class DC_PROB{
 private:
  double log_DmqL_min, log_DmqL_max;
  double DmqL_min, DmqL_max;
  double Del;

  double sinq_av[NBIN_DC/2][N_RATE_COEF];
  double sinq(const double DmqL, const int bin); // interpolating on table

 public:
  void calc_integrals(void);
  double P(Param_5nu &prm, const int det, const int bin);

} DC_prob;


double DC_PROB::P(Param_5nu &p, const int det, const int bin)
{
  double P = 1.;
  for(int i = 0; i < N_NU-1; i++){
    for(int j = i+1; j < N_NU; j++){
      const double DmqL = fabs(p.Dmq(j,i)) * bl_dc[det] * 1.e3;
      if(DmqL > DmqL_min)
        P -= 4. * norm(p.Ue[i] * p.Ue[j]) * sinq(DmqL, bin);
    }
  }
  return P; 
} 

// interpolating on the table
double DC_PROB::sinq(const double DmqL, const int bin)
{
  if(DmqL <  DmqL_min) error("DC_PROB::sinq: DmqL too small");
  if(DmqL >= DmqL_max) error("DC_PROB::sinq: DmqL too large");

  const double logD = log10(DmqL);

  int i = int( (logD - log_DmqL_min) / Del);
  if(i < 0 || i >= N_RATE_COEF-1) error("DC_PROB::sinq: index out of range");

  const double x = (logD - i * Del - log_DmqL_min) / Del;
  if(x < 0. || x > 1.) error("DC_PROB::sinq: x out of range");

  return sinq_av[bin][i] + (sinq_av[bin][i+1] - sinq_av[bin][i]) * x;
}


  
/***************** init *************************************/

// For systematics, see slide 62 of Cern talk
void dc_init(const int old_new)
{
  // set bins
  for(int i = 0; i < NBIN_DC/2; i++){
    bin[i][0] = 1.0 + i * 0.25;
    bin[i][1] = bin[i][0] + 0.25;
  }

  // reading data based on slide 26 of Moriond 2016 talk  
  double data[NBIN_DC][2];       // data and stat error
  double ww, vv;

  // data points for FD-I data
  FILE *fp_data      = fopen(REACTOR_PATH"Data_DC-cern-2016/fd-i-data.csv","r");
  FILE *fp_data_err  = fopen(REACTOR_PATH"Data_DC-cern-2016/fd-i-lower.csv","r");
  FILE *fp_sys       = fopen(REACTOR_PATH"Data_DC-cern-2016/fd-i-yellowhist.csv","r");
  FILE *fp_sys_e     = fopen(REACTOR_PATH"Data_DC-cern-2016/fd-i-oldbluehist.csv","r"); 
  
  for(int i = 0; i < NBIN_DC/2; i++){
    fscanf(fp_data,     "%*f, %lf", &data[i][0]);
    fscanf(fp_data_err, "%*f, %lf", &data[i][1]);
    data[i][1] = (data[i][0] - data[i][1]);
    
    fscanf(fp_sys_e, "%*f, %lf", &ww);
    fscanf(fp_sys_e, "%*f, %lf", &vv);
    sys_DC_e[i] = 0.5*(ww - vv);
    
    fscanf(fp_sys, "%*f, %lf", &ww);
    fscanf(fp_sys, "%*f, %lf", &vv);
    sys_DC[i] = 0.5*(ww - vv);    
  }
  fclose(fp_data);
  fclose(fp_data_err);
  fclose(fp_sys_e);
  fclose(fp_sys);

  // data points for FD-II/ND data
  fp_data      = fopen(REACTOR_PATH"Data_DC-cern-2016/ratio-data.csv","r");
  fp_data_err  = fopen(REACTOR_PATH"Data_DC-cern-2016/ratio-lower.csv","r");
  fp_sys       = fopen(REACTOR_PATH"Data_DC-cern-2016/ratio-greenhist.csv","r");
  fp_sys_e     = fopen(REACTOR_PATH"Data_DC-cern-2016/ratio-oldbluehist.csv","r"); 

  for(int i = NBIN_DC/2; i < NBIN_DC; i++){
    
    fscanf(fp_data,     "%*f, %lf", &data[i][0]);
    fscanf(fp_data_err, "%*f, %lf", &data[i][1]);
    data[i][1] = (data[i][0] - data[i][1]);
    
    fscanf(fp_sys_e, "%*f, %lf", &ww);
    fscanf(fp_sys_e, "%*f, %lf", &vv);
    sys_DC_e[i] = 0.5*(ww - vv);
 
    fscanf(fp_sys, "%*f, %lf", &ww);
    fscanf(fp_sys, "%*f, %lf", &vv);
    sys_DC[i] = 0.5*(ww - vv);
  }
  fclose(fp_data);
  fclose(fp_data_err);
  fclose(fp_sys_e);
  fclose(fp_sys);
  
  // calc coefficients
  old_new_gl_dc = old_new;
  DC_prob.calc_integrals();
    
  // set data and covariance matrix in class fit    
  for(int i = 0; i < NBIN_DC; i++){
    const int ii = fit.first_bin[DC] + i;    
    fit.Data[ii] = data[i][0];
    fit.S_data[ii][ii] = norm(data[i][1]);
  }
    
  /*** set the pull errors ***/
  const int p0 = fit.first_pull[DC];
  for(int p = 0; p < NPULL_DC; p++){
    fit.S_pull[p0+p][p0+p] = 1.;
  }
  fit.S_pull[p0+2][p0+2] = fit.S_pull[p0+3][p0+3] = fit.S_pull[p0+4][p0+4]
    = fit.S_pull[p0+5][p0+5] = norm(0.012) + norm(0.008) + norm(0.007) + norm(0.004);
  
  return;
}

/**************** calc the table for class fit ***************/
   
void set_table_dc(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1])
{
  Param_5nu prm_3f(prm); // Parameter set for standard oscillation for undoing
                         // ND-based normalization
  for (int i=I14; i < N_NU; i++)
    prm_3f.theta[i] = 0.0;
  prm_3f.set_ang();

  for(int i = 0; i < NBIN_DC/2; i++){
    
    // prediction FD-I
    int b = fit.first_bin[DC] + i;
    const double w = DC_prob.P(prm, DC_FAR, i);
    cff[b][NPULLS] = w * DC_prob.P(prm_3f, DC_NEAR, i) / DC_prob.P(prm, DC_NEAR, i);

    // prediction FD-II/ND
    b += NBIN_DC/2;
    cff[b][NPULLS] = w / DC_prob.P(prm, DC_NEAR, i);
  }
  
  // pulls
  for(int i = 0; i < NBIN_DC/2; i++){
    const int b = fit.first_bin[DC] + i;
    const int sign = (i < 10 ? 1. : -1.); 
                // oposite sign of energy scale uncert below/above event maximum
                // this number is fudged
    cff[b]          [fit.first_pull[DC]]   = sign * sys_DC_e[i];           // E-scale FD-I
    cff[b+NBIN_DC/2][fit.first_pull[DC]+1] = sign * sys_DC_e[i+NBIN_DC/2]; // E-scale FD-II/ND
  }

  const int bin_low_high  = 16;  // separates low and high energy systematic 
  for(int i = 0; i < bin_low_high; i++){
    const int b = fit.first_bin[DC] + i;
    cff[b]          [fit.first_pull[DC]+2] = cff[b][NPULLS] + 0*sys_DC[i];           // norm FD-I
    cff[b+NBIN_DC/2][fit.first_pull[DC]+3] = cff[b][NPULLS] + 0*sys_DC[i+NBIN_DC/2]; // norm NF
  }
  for(int i = bin_low_high; i < NBIN_DC/2; i++){
    const int b = fit.first_bin[DC] + i;
    cff[b]          [fit.first_pull[DC]+4] = cff[b][NPULLS] + 0*sys_DC[i];           // norm FD-I
    cff[b+NBIN_DC/2][fit.first_pull[DC]+5] = cff[b][NPULLS] + 0*sys_DC[i+NBIN_DC/2]; // norm NF
  }

  return;
}

  
/********************************************************/
/*    calculation of coefficients                       */
/********************************************************/

double gauss_int_dc(double eNu);
  
double dc_flux(const double Enu)
{
  if(iso_gl_dc == -1){
    double w = 0.;
    for(int i = 0; i < NISO; i++)
      w += isofract_dc[i] * global_flux(i, Enu, old_new_gl_dc);
    return w;
  }
  return isofract_dc[iso_gl_dc] * global_flux(iso_gl_dc, Enu, old_new_gl_dc);
}


double fold_dc0(double eNu){
  const double EposKin = eNu - ME - DELTA;
  return dc_flux(eNu) * crossSect(EposKin) * gauss_int_dc(eNu);
}
double fold_sinq_dc(double eNu)
{
  double Cos;
  if(gl_DmqL_dc <= 0.) Cos = 1.;
  else
  {
    const double a1 = 2.54 * gl_DmqL_dc / eNu * 0.95;
    const double a2 = 2.54 * gl_DmqL_dc / eNu * 1.05;
    Cos = (sin(a2) - sin(a1))/(a2 - a1);
  }

  const double EposKin = eNu - ME - DELTA;
  return crossSect(EposKin) * 0.5*(1.-Cos) * dc_flux(eNu) * gauss_int_dc(eNu);
}



void DC_PROB::calc_integrals(void)
{
  // DmqL in units of eV^2 m 
  log_DmqL_min = ATM_MIN + log10(10.);      
  log_DmqL_max = STE_MAX + log10(2000.);   

  DmqL_min = exp10(log_DmqL_min);
  DmqL_max = exp10(log_DmqL_max);

  Del = (log_DmqL_max - log_DmqL_min) / (N_RATE_COEF - 1.);

  double low,up;
  iso_gl_dc = -1;
  
  for(int i = 0; i < NBIN_DC/2; i++){
    bin_gl_dc = i;

    low = bin[bin_gl_dc][0] + 0.8 - 3. * e_res_dc(bin[bin_gl_dc][0]);
    if(low<EnuMIN) low = EnuMIN;
    up = bin[bin_gl_dc][1] + 0.8 + 3. * e_res_dc(bin[bin_gl_dc][1]);
    if(up>EnuMAX) up=EnuMAX;
    
    const double n0 = qromb1(fold_dc0, low, up, 1.e-6);

    for(int j = 0; j  < N_RATE_COEF; j++){

      gl_DmqL_dc = exp10(log_DmqL_min + j * Del);
      sinq_av[i][j] = qromb1(fold_sinq_dc, low, up, 1.e-6) / n0;
    }
  }
  return;
}


/****************************************************/
/* SUBROUTINE gauss_int_dc(T_min, T_max, T_e, integ_T) */
/****************************************************/
      
double gauss_int_dc(double eNu)
{   
   double T_min, T_max;
   double T_e = eNu + ME - DELTA;
  
   //correction motivated by migration matrix from Guillaume
   //T_e += 0.1/3*(T_e - 5.);

   double sigma = e_res_dc(T_e);

   if(bin_gl_dc == -1){
     T_min = bin[0][0];
     T_max = bin[NBIN_DC/2-1][1];
   }else{
     T_min = bin[bin_gl_dc][0];
     T_max = bin[bin_gl_dc][1];
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


/********************************************
 *  plotting
 ********************************************/

void plot_dc_pred(Param_5nu &prm)
{    
  double cff[NBIN_CHISQ][NPULLS+1];
  set_table_dc(prm, cff);

  for(int i = 0; i < NBIN_DC/2; i++){
    int b = fit.first_bin[DC] + i;
    printf("%e %e %e\n", bin[i][0], cff[b][NPULLS], cff[b+NBIN_DC/2][NPULLS]);
    printf("%e %e %e\n", bin[i][1], cff[b][NPULLS], cff[b+NBIN_DC/2][NPULLS]);
  }
  return;
}

  
#undef E_RESOL
#undef ACCURACY

#endif // USE_DC
   
} // end namespace
   
