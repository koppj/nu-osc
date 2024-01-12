#include "definitions.h"

/*********************************************
 * using results presented at CERN Sept 2016
 *********************************************/

namespace ns_reactor
{

#ifdef USE_DC

#define ACCURACY 1.e-5

/*** energy resolution from slide 51 CERN talk (averaged) ***/   

const double EresA = norm(0.081);
const double EresB = norm(0.017);
const double EresC = norm(0.023);

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
double bin[NBIN_DC/2][2];

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

void dc_init(const int old_new)
{
  // set bins
  for(int i = 0; i < NBIN_DC/2; i++){
    bin[i][0] = 1. + i * 0.25;
    bin[i][1] = bin[i][0] + 0.25;
  }

  // reading data based on slide 59 of CERN 2016 talk  
  double data[NBIN_DC/2][6];       // event spectra and obs/no-osc-pred for ND, FD-I, FD-II

  // event spectra
  FILE *fp_data_ND   = fopen(REACTOR_PATH"Data_DC-CERN-2016/ND-spectrum.csv","r");
  FILE *fp_data_FDI  = fopen(REACTOR_PATH"Data_DC-CERN-2016/FD-I-spectrum.csv","r");
  FILE *fp_data_FDII = fopen(REACTOR_PATH"Data_DC-CERN-2016/FD-II-spectrum.csv","r");
  
  for(int i = 0; i < NBIN_DC/2; i++){
    
    if(fscanf(fp_data_ND,     "%*f, %lf", &data[i][0]) != 1){
      fprintf(stderr, "DC error in reading ND spectrum\n");
      exit(0);
    }
    if(fscanf(fp_data_FDI,    "%*f, %lf", &data[i][1]) != 1){
      fprintf(stderr, "DC error in reading FDI spectrum\n");
      exit(0);
    }
    if(fscanf(fp_data_FDII,   "%*f, %lf", &data[i][2]) != 1){
      fprintf(stderr, "DC error in reading FDII spectrum\n");
      exit(0);
    }
  }
  fclose(fp_data_ND);
  fclose(fp_data_FDI);
  fclose(fp_data_FDII);

  // ratio of observation to no-oscillation prediction
  fp_data_ND   = fopen(REACTOR_PATH"Data_DC-CERN-2016/ND-ratio.csv","r");
  fp_data_FDI  = fopen(REACTOR_PATH"Data_DC-CERN-2016/FD-I-ratio.csv","r");
  fp_data_FDII = fopen(REACTOR_PATH"Data_DC-CERN-2016/FD-II-ratio.csv","r");
  
  for(int i = 0; i < NBIN_DC/2; i++){
    
    if(fscanf(fp_data_ND,     "%*f, %lf", &data[i][3]) != 1){
      fprintf(stderr, "DC error in reading ND ratio\n");
      exit(0);
    }
    if(fscanf(fp_data_FDI,    "%*f, %lf", &data[i][4]) != 1){
      fprintf(stderr, "DC error in reading FDI ratio\n");
      exit(0);
    }
    if(fscanf(fp_data_FDII,   "%*f, %lf", &data[i][5]) != 1){
      fprintf(stderr, "DC error in reading FDII ratio\n");
      exit(0);
    }
  }
  fclose(fp_data_ND);
  fclose(fp_data_FDI);
  fclose(fp_data_FDII);
    
  // calc coefficients
  old_new_gl_dc = old_new;
  DC_prob.calc_integrals();
    
  // set data and covariance matrix in class fit    
  for(int i = 0; i < NBIN_DC/2; i++){
    int ii = fit.first_bin[DC] + i;    
    fit.Data[ii] = data[i][4]/data[i][3];  // FD-I / ND
    fit.S_data[ii][ii] = norm(fit.Data[ii]) * (1./data[i][0] + 1./data[i][1]);
    //fit.S_data[ii][ii] = norm(fit.Data[ii]) * 10.;

    ii += NBIN_DC/2;
    fit.Data[ii] = data[i][5]/data[i][3];  // FD-II / ND
    fit.S_data[ii][ii] = norm(fit.Data[ii]) * (1./data[i][0] + 1./data[i][2]);

    // correlation of statistical error due to using the same ND data for FD-I and FD-II
    const int jj = ii - NBIN_DC/2;
    fit.S_data[jj][ii] = fit.S_data[ii][jj] = fit.Data[ii] * fit.Data[jj] / data[i][0];    
  }
  
    
  /*** set the pull errors ***/
  // normalization
  int pp = fit.first_pull[DC];
  fit.S_pull[pp][pp] = norm(0.0075);  // slide 62 CERN talk
  fit.pull_status[pp] = FIXED;        // not used!!

  pp++;
  // energy scale --> tilt
  fit.S_pull[pp][pp] = 1.;
  fit.pull_status[pp] = ACTIVE;

  pp++;
  // FD-I / FD-II norm
  fit.S_pull[pp][pp] = norm(0.01); // guessing a number
  fit.pull_status[pp] = ACTIVE;

  return;
}

/**************** calc the table for class fit ***************/
   
void set_table_dc(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1])
{
  for(int i = 0; i < NBIN_DC/2; i++){
    
    // prediction FD-I/ND = FD-II/ND  
    const int b = fit.first_bin[DC] + i;
    cff[b][NPULLS] =  cff[b+NBIN_DC/2][NPULLS] = DC_prob.P(prm, DC_FAR, i) / DC_prob.P(prm, DC_NEAR, i);
  }
  
  // pulls
  int pp = fit.first_pull[DC];

  // normalization
  for(int i = 0; i < NBIN_DC; i++){
    const int b = fit.first_bin[DC] + i;
    cff[b][pp] = cff[b][NPULLS]; 
  }
  
  pp++;
  // "energy scale"
  for(int i = 0; i < NBIN_DC/2; i++){
    const int b = fit.first_bin[DC] + i;
    const double b_0 = 12.;
    cff[b][pp] = cff[b+NBIN_DC/2][pp] = 0.05 * (b_0 - i)/b_0; // 5% at first bin, minimum at b_0, for b>b_0 opposite sign
  }

  pp++;
  // FD-I / FD-II norm
  for(int i = 0; i < NBIN_DC/2; i++){
    int b = fit.first_bin[DC] + i;
    cff[b][pp] =  cff[b][NPULLS];
    b += NBIN_DC/2;
    cff[b][pp] = -cff[b][NPULLS];
  }

  
  
  // for(int i = 0; i < NBIN_DC/2; i++){
  //   const int b = fit.first_bin[DC] + i;
  //   const int sign = (i < 10 ? 1. : -1.); // oposite sign of energy scale uncert below/above event maximum
  //   cff[b]          [fit.first_pull[DC]]   = sign * sys_DC_e[i];            // E-scale FD-I
  //   cff[b+NBIN_DC/2][fit.first_pull[DC]+1] = sign * sys_DC_e[i+NBIN_DC/2];  // E-scale FD-II/ND

  //   cff[b][fit.first_pull[DC]+6] = cff[b][NPULLS];   // flux normalization for FD-I
  // }

  // const int bin_low_high  = 14;  // separates low and high energy systematic 
  
  // for(int i = 0; i < bin_low_high; i++){
  //   const int b = fit.first_bin[DC] + i;
  //   cff[b]          [fit.first_pull[DC]+2] = sys_DC[i];              // low-E sys FD-I
  //   cff[b+NBIN_DC/2][fit.first_pull[DC]+3] = sys_DC[i+NBIN_DC/2];    // low-E sys FD-II/ND
  // }
  // for(int i = bin_low_high; i < NBIN_DC/2; i++){
  //   const int b = fit.first_bin[DC] + i;  
  //   cff[b]          [fit.first_pull[DC]+4] = sys_DC[i];              // high-E sys FD-I
  //   cff[b+NBIN_DC/2][fit.first_pull[DC]+5] = sys_DC[i+NBIN_DC/2];    // high-E sys FD-II/ND
  // }
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
     T_max = bin[NBIN_DC-1][1];
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
   
