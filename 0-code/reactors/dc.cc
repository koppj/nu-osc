#include "definitions.h"

namespace ns_reactor
{

#ifdef USE_DC

#define ACCURACY 1.e-5

// 0.32 gives good agreement with DC no-osc MC 
#define E_RESOL 0.3

inline double e_res_dc(double e)
{
  //return 0.7;
  //return 0.3 + 0.16*sqrt(e);
  return E_RESOL*sqrt(e);
}


#define NBG_DC 3

struct coeff{
  double n0;
  double c0[NBIN_DC];
  double c0I[NBIN_DC][NISO][N_CO_UNC];
  double c0_cor[NBIN_DC];
  double c1[NBIN_DC][DM2_NUM];
} co_dc;

extern Fit fit;
extern double sbl_data[N_SBL_R][2];
 
double gauss_int_dc(double eNu);
void calc_coef_dc(void);

/********** global vars *******************/

double dm2_glob_dc, L_gl_dc;
int bin_gl_dc, iso_gl_dc, old_new_gl_dc;

const double isofract_dc[NISO]={0.488, 0.087, 0.359, 0.076};

double dc_no_osc_react[NBIN_DC];
double bg_dc[NBIN_DC][NBG_DC];

/**************** input data ***************/

/* binning in positron energy */
double bin[NBIN_DC][2] = {
  {0.7, 1.2},
  {1.2, 1.7},
  {1.7, 2.2},
  {2.2, 2.7},
  {2.7, 3.2},
  {3.2, 3.7},
  {3.7, 4.2},
  {4.2, 4.7},
  {4.7, 5.2},
  {5.2, 5.7},
  {5.7, 6.2},
  {6.2, 6.7},
  {6.7, 7.2},
  {7.2, 7.7},
  {7.7, 8.2},
  {8.2, 9.2},
  {9.2, 10.2},
  {10.2, 12.2}
};

// bg rate per day and error 
const double bg_dc_rate[NBG_DC][2] = {
  {0.261, 0.002}, // acc
  {1.2,   0.59},   // Li9
  {0.67,  0.2}    // fast n
};

/***************** init *************************************/

void dc_init(const int old_new)
{
  double pred_tot_dc[NBIN_DC];    // total pred for no osc including bg
  double data_dc[NBIN_DC];        // events per bin
  double FN_ACC_BG[NBIN_DC];      // sum of fast n and acc
  double bg_dc_tot[NBIN_DC];      // total BG
  double bf_spect[NBIN_DC];       // spectrum for sin^22th13 = 0.109, Dmq = 0.00232
  double ratio[NBIN_DC];          // ratio of data/prediction

  // read the data, extracted from DC spectrum plot
  FILE *fp_fn_ac = fopen(REACTOR_PATH"Data_DC/FN+acc.txt","r");
  FILE *fp_bf_sp = fopen(REACTOR_PATH"Data_DC/best-fit-spect.txt","r");
  FILE *fp_data  = fopen(REACTOR_PATH"Data_DC/data.txt","r");
  FILE *fp_no_os = fopen(REACTOR_PATH"Data_DC/no-osc-pred-tot.txt","r");
  FILE *fp_totbg = fopen(REACTOR_PATH"Data_DC/tot-bg.txt","r");
  FILE *fp_ratio = fopen(REACTOR_PATH"Data_DC/dc-ratio.dat","r");

  for(int i = 0; i < NBIN_DC; i++){
    fscanf(fp_totbg, "%lf", &bg_dc_tot[i]);
    fscanf(fp_fn_ac, "%lf", &FN_ACC_BG[i]);
    fscanf(fp_data,  "%lf", &data_dc[i]);
    fscanf(fp_no_os, "%lf", &pred_tot_dc[i]);
    fscanf(fp_bf_sp, "%lf", &bf_spect[i]);
    fscanf(fp_ratio, "%lf", &ratio[i]);
  }
  fclose(fp_totbg);
  fclose(fp_fn_ac);
  fclose(fp_data);
  fclose(fp_no_os);
  fclose(fp_bf_sp);
  fclose(fp_ratio);

  // calc coefficients
  old_new_gl_dc = old_new;
  calc_coef_dc();
  
  // total number of observed neutrinos
  const double tot_ev_obs = 8249.;

  // normalize the data read off from the plot 
  double n_nu_dc = 0.; 
  
  for(int i = 0; i < NBIN_DC; i++){
    data_dc[i] *= (bin[i][1] - bin[i][0])/0.5;
    n_nu_dc += data_dc[i];
  }  
  
  // correct last bins
  data_dc[NBIN_DC-1] -= 15.;
  data_dc[NBIN_DC-2] -= 7.;
  n_nu_dc -= 22.;
  
  //fprintf(stderr, "DC neutrinos = %f (should be %f)\n", n_nu_dc, tot_ev_obs);

  for(int i = 0; i < NBIN_DC; i++)
    data_dc[i] *= tot_ev_obs / n_nu_dc;
  
 // plotting the data
 /* 
  for(int i = 0; i < NBIN_DC; i++)
    fprintf(stderr, "%e %e %e\n", 0.5*(bin[i][0]+bin[i][1]),data_dc[i],sqrt(data_dc[i]));
  exit(0);
  */
  
  // set data and covariance matrix in class fit    
  for(int i = 0; i < NBIN_DC; i++){
    const int ii = fit.first_bin[DC] + i;    
    fit.Data[ii] = fit.S_data[ii][ii] = data_dc[i];
  }

  /*** set the backgrounds ***/
  double sum_bg[NBG_DC] = {0., 0., 0.};

  for(int i = 0; i < NBIN_DC; i++){
    // accidental 
    bg_dc[i][0] = (FN_ACC_BG[i] - FN_ACC_BG[NBIN_DC-1]) * (bin[i][1] - bin[i][0])/0.5; 
    sum_bg[0] += bg_dc[i][0];

    // Li9
    bg_dc[i][1] = (bg_dc_tot[i] - FN_ACC_BG[i])*(bin[i][1] - bin[i][0])/0.5; 
    sum_bg[1] += bg_dc[i][1];

    // fast n 
    bg_dc[i][2] = FN_ACC_BG[NBIN_DC-1] * (bin[i][1] - bin[i][0])/0.5;
    sum_bg[2] += bg_dc[i][2];
  }
  
  
  /*** set the pull errors ***/

  // DC normalization
#ifndef CHECK_DC_ERRORS 
  // only detector norm (increased from 0.01-->0.02 to improve fit)
  fit.S_pull[fit.first_pull[DC]][fit.first_pull[DC]] = norm(0.02);  
#else  
  fit.S_pull[fit.first_pull[DC]][fit.first_pull[DC]] = norm(0.02) + norm(0.017);  // detector + reactor  
  for(int i = 0; i < PULL_GLOBAL; i++)
    fit.pull_status[i] = FIXED;
#endif  

  // backgrounds
  for(int j = 0; j < NBG_DC; j++){
    const int p = fit.first_pull[DC] + 1 + j;
    fit.S_pull[p][p] = norm(bg_dc_rate[j][1] / bg_dc_rate[j][0]); 
  }

  // energy scale
  fit.S_pull[fit.first_pull[DC] + 4][fit.first_pull[DC] + 4] = norm(0.015);
    
  /**** the no oscillation prediction ****************
   * dc_no_osc_react contains the no-oscillation 
   * prediction scaled to the Bugey4 measurement
   ***************************************************/ 
    
  // subtracting background, scaling to events per bin
  for(int i = 0; i < NBIN_DC; i++){
    dc_no_osc_react[i] = pred_tot_dc[i] * (bin[i][1] - bin[i][0])/0.5;
    bf_spect[i] *= (bin[i][1] - bin[i][0])/0.5;
    
    for(int j = 0; j < NBG_DC; j++){
      dc_no_osc_react[i] -= bg_dc[i][j];
      bf_spect[i] -= bg_dc[i][j];
    }
    if(dc_no_osc_react[i] < 0.) dc_no_osc_react[i] = 0.;
    if(bf_spect[i] < 0.) bf_spect[i] = 0.;
  }

  // for last 5 bins use the number obtained from the ratio
  for(int i = NBIN_DC-6; i < NBIN_DC; i++){
    dc_no_osc_react[i] = data_dc[i] / ratio[i];
    for(int j = 0; j < NBG_DC; j++)
      dc_no_osc_react[i] -= bg_dc[i][j];
  }
  
  
  //  plotting the original no-osc and best fit spectra
  /*
  for(int i = 0; i < NBIN_DC; i++){
    printf("%e %e %e\n", bin[i][0], bf_spect[i], dc_no_osc_react[i]);
    printf("%e %e %e\n", bin[i][1], bf_spect[i], dc_no_osc_react[i]);
  }
  exit(0);
  */
#define USE_DC_PRED  
#ifndef USE_DC_PRED  
  // using my own no-osc prediction  
  double dc_expect_ev_tot = 0.; 
  for(int i = 0; i < NBIN_DC; i++)  
    dc_expect_ev_tot += dc_no_osc_react[i];
  
  double w = 0.;
  for(int i = 0; i < NBIN_DC; i++)
    w += co_dc.c0[i];

  for(int i = 0; i < NBIN_DC; i++)    
    dc_no_osc_react[i] = co_dc.c0[i] * dc_expect_ev_tot / w;
#endif  
  return;
}


/**************** calc the table for class fit ***************/

   
void set_table_dc(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1])
{
  // for interpolation
  const double logDmq = log10(fabs(prm.dmq[2]));
  const int idmq = int( (logDmq - ATM_MIN) / D_ATM);
  const double x = (logDmq - idmq * D_ATM - ATM_MIN) / D_ATM;

  if(logDmq <  ATM_MIN) error("set_table_dc::Dmq too small");
  if(logDmq >= ATM_MAX) error("set_table_dc::Dmq too large");
  if(idmq < 0 || idmq >= DM2_NUM-1) error("set_table_dc: index out of range");
  if(x < 0. || x > 1.) error("set_table_dc: x out of range");

  // coeffs for the probability
  const double sq_eff = 4. * norm(prm.Ue[2]) * (norm(prm.Ue[0]) + norm(prm.Ue[1])); 
  double w = 0;
  for(int i = 0; i < 3; i++) w += norm(prm.Ue[i]);
  const double cons = 2. * (norm(prm.Ue[3]) + norm(prm.Ue[4])) * w + 2. * norm(prm.Ue[3]) * norm(prm.Ue[4]);
      
  // set the table
  for(int i = 0; i < NBIN_DC; i++){
    const int b = fit.first_bin[DC] + i;

    // the prediction
    const double sinq = co_dc.c1[i][idmq] + (co_dc.c1[i][idmq+1] - co_dc.c1[i][idmq]) * x;
    cff[b][NPULLS] = 1. - sq_eff * sinq - cons;

    // factor out Bugey4 measurement
    cff[b][NPULLS] *= dc_no_osc_react[i] / sbl_data[BUGEY4][old_new_gl_dc];

    // DC normalization and flux normalization
    cff[b][fit.first_pull[DC]] = cff[b][FLUX_NORM] = cff[b][NPULLS];

    for(int j = 0; j < NBG_DC; j++){

      // add backgrounds to prediction
      cff[b][NPULLS] += bg_dc[i][j];

      // background pulls
      const int p = fit.first_pull[DC] + 1 + j;
      cff[b][p] = bg_dc[i][j];
    }

    // isotopes

    // the pulls
    for(int k = 0; k  < N_CO_UNC; k++){
        cff[b][PULL_U235_0 + k] = co_dc.c0I[i][U235][k] * dc_no_osc_react[i];
        cff[b][PULL_P239_0 + k] = co_dc.c0I[i][P239][k] * dc_no_osc_react[i];
        cff[b][PULL_P241_0 + k] = co_dc.c0I[i][P241][k] * dc_no_osc_react[i];
    }				
    cff[b][PULL_U238] = co_dc.c0I[i][U238][0] * dc_no_osc_react[i];
    cff[b][FLUX_COR] = co_dc.c0_cor[i] * dc_no_osc_react[i];

    // energy scale
    const double n_i = cff[b][FLUX_NORM];

    if(i < NBIN_DC-1){
      const double n_ip1 = cff[b+1][FLUX_NORM];		
      cff[b][fit.first_pull[DC] + 4] = (bin[i][1]*n_ip1 - bin[i][0]*n_i)/(bin[i][1]-bin[i][0]);
	
    }else
      cff[b][fit.first_pull[DC] + 4] = cff[b-1][fit.first_pull[DC] + 4];
  }
  return;
}

/********************************************************/
/*    calculation of coefficients                       */
/********************************************************/

  
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

#define DC_PE_MEV 214.  // PE/MEV

/* correct energy scale non-linearity */
double eNu_corr(double eNu)
{
  return eNu; // not used

  const double Evis = eNu - DELTA + ME;
  const double PE = Evis * DC_PE_MEV;    

  // non-linearity of E-scale
  if(PE <= 56.1478) error("PE too small");
  double PE_corr = PE * (0.0287 * log(PE - 56.1478) + 0.8423);
  //double PE_corr = PE * (0.05 * log(PE - 56.1478) + 0.8423);

  // additional correction due to z-dependence
  PE_corr *= 0.95;

  const double eNu_corr = PE_corr / DC_PE_MEV + DELTA - ME;
  return (eNu_corr < EnuMIN ? EnuMIN : eNu_corr); 
}

double fold_dc0(double eNu0){
  const double eNu = eNu_corr(eNu0);
  const double EposKin = eNu - ME - DELTA;
  return dc_flux(eNu) * crossSect(EposKin)*gauss_int_dc(eNu);
}
double fold_dc1(double eNu0){
  const double eNu = eNu_corr(eNu0);
  const double EposKin = eNu - ME - DELTA;
  return eNu * dc_flux(eNu) * crossSect(EposKin)*gauss_int_dc(eNu);
}
double fold_dc2(double eNu0){
  const double eNu = eNu_corr(eNu0);
  const double EposKin = eNu - ME - DELTA;
  return eNu * eNu * dc_flux(eNu) * crossSect(EposKin)*gauss_int_dc(eNu);
}

double fold_dc_corr(double eNu0)
{
  const double eNu = eNu_corr(eNu0);
  iso_gl_dc = U235;
  double w = read(REACTOR_PATH"Dat/Patrick-U235-err_cor.dat", eNu) * dc_flux(eNu);
  iso_gl_dc = P239;
  w += read(REACTOR_PATH"Dat/Patrick-Pu239-err_cor.dat", eNu) * dc_flux(eNu);
  iso_gl_dc = P241;
  w += read(REACTOR_PATH"Dat/Patrick-Pu241-err_cor.dat", eNu) * dc_flux(eNu);

  const double EposKin=eNu-ME-DELTA;
  return w * crossSect(EposKin)*gauss_int_dc(eNu);
}


double fold_dc_osc(double eNu0)
{
  const double eNu = eNu_corr(eNu0);
  const double EposKin = eNu - ME - DELTA;
  const double arg_atm = 1.27e3 * dm2_glob_dc * L_gl_dc / eNu;

  return dc_flux(eNu) * crossSect(EposKin) * gauss_int_dc(eNu) * norm(sin(arg_atm));
}


#define NSIGM 4.

void calc_coef_dc(void)
{
  const double bl[2] = {.9979, 1.1146};
  
  double low,up;
  
  bin_gl_dc = -1;
  iso_gl_dc = -1;

  co_dc.n0 = qromb1(fold_dc0,EnuMIN,EnuMAX,ACCURACY);

  for(bin_gl_dc = 0; bin_gl_dc < NBIN_DC; bin_gl_dc++){

    low = bin[bin_gl_dc][0] + 0.8 - NSIGM * e_res_dc(bin[bin_gl_dc][0]);
    if(low<EnuMIN) low = EnuMIN;
    up = bin[bin_gl_dc][1] + 0.8 + NSIGM * e_res_dc(bin[bin_gl_dc][1]);
    if(up>EnuMAX) up=EnuMAX;

    co_dc.c0[bin_gl_dc] = 0.;

    for(iso_gl_dc = 0; iso_gl_dc < NISO; iso_gl_dc++){

      co_dc.c0I[bin_gl_dc][iso_gl_dc][0] = qromb1(fold_dc0, low, up, ACCURACY);
      co_dc.c0[bin_gl_dc] += co_dc.c0I[bin_gl_dc][iso_gl_dc][0];

      if(iso_gl_dc != U238 && old_new_gl_dc == NEW){
 	  co_dc.c0I[bin_gl_dc][iso_gl_dc][1] = qromb1(fold_dc1, low, up, ACCURACY);
 	  co_dc.c0I[bin_gl_dc][iso_gl_dc][2] = qromb1(fold_dc2, low, up, ACCURACY);
      }
    }

    co_dc.c0_cor[bin_gl_dc] = qromb1(fold_dc_corr, low, up, ACCURACY);

    iso_gl_dc = -1;

    Parameters p;
    for(int j = 0; j < DM2_NUM; j++){

      p.idmq_31 = j;
      dm2_glob_dc = p.dmq31();
      co_dc.c1[bin_gl_dc][j] = 0.;

      // assume that both reactors contribute with same weight
      for(int r = 0; r < 2; r++){
        L_gl_dc = bl[r];      
        co_dc.c1[bin_gl_dc][j] += 0.5 * qromb1(fold_dc_osc, low, up, ACCURACY);
      }
    }
  }


  // convert to probability in each bin
  for(bin_gl_dc = 0; bin_gl_dc < NBIN_DC; bin_gl_dc++){

    co_dc.c0_cor[bin_gl_dc] /= co_dc.c0[bin_gl_dc];

    for(iso_gl_dc = 0; iso_gl_dc < NISO; iso_gl_dc++)
      for(int k = 0; k  < N_CO_UNC; k++)
        co_dc.c0I[bin_gl_dc][iso_gl_dc][k] /= co_dc.c0[bin_gl_dc];

    for(int j = 0; j < DM2_NUM; j++)
      co_dc.c1[bin_gl_dc][j] /= co_dc.c0[bin_gl_dc];
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

#undef E_RESOL
#undef ACCURACY

/****************************************
 * plotting functions
 ****************************************/

void plot_dc_pred(Parameters p, double f)
{
  fprintf(stderr, "dmq = %e\n", p.dmq31());
  
  for(int i = 0; i < NBIN_DC; i++){
    
    const double N_nos = dc_no_osc_react[i] * f  / sbl_data[BUGEY4][old_new_gl_dc];
    const double N_osc = N_nos * (1. - p.sq2_t13() * co_dc.c1[i][p.idmq_31]);

    printf("%f %f %f\n", bin[i][0], N_osc, N_nos);
    printf("%f %f %f\n", bin[i][1], N_osc, N_nos);

    /*	   
    printf("%e %e %e %e %e\n", bin[i][0], 
	   bg_dc[i][2], 
	   bg_dc[i][2]+bg_dc[i][0],
	   bg_dc[i][2]+bg_dc[i][0]+bg_dc[i][1], 
	   bg_dc[i][2]+bg_dc[i][0]+bg_dc[i][1] + w);
    printf("%e %e %e %e %e\n", bin[i][1], 
	   bg_dc[i][2], 
	   bg_dc[i][2]+bg_dc[i][0],
	   bg_dc[i][2]+bg_dc[i][0]+bg_dc[i][1], 
	   bg_dc[i][2]+bg_dc[i][0]+bg_dc[i][1] + w); 
    */ 
  }
  exit(0);  
}


/************************************************
 * a simple chisq with just overall norm uncert *
 ************************************************/
#define SIG_TEST sqrt(norm(0.017) + norm(0.013))

double dc_test_chisq(Parameters prm)
{
  double react[NBIN_DC];
  double data[NBIN_DC];

  for(int i = 0; i < NBIN_DC; i++){

    // the prediction
    react[i] = (1. - prm.sq2_t13() * co_dc.c1[i][prm.idmq_31]) * dc_no_osc_react[i];

    data[i] = fit.Data[fit.first_bin[DC]+i];

    for(int j = 0; j < NBG_DC; j++){
      data[i] -= bg_dc[i][j];
    }
  }

  double w, v;
  w = v = 1./norm(SIG_TEST);

  for(int i = 0; i < NBIN_DC; i++){
    w += norm(react[i]) / data[i];
    v += react[i];
  }

  double a = v/w;

  w = norm( (a - 1.) / SIG_TEST );
  for(int i = 0; i < NBIN_DC; i++)
    w += norm(a*react[i] - data[i]) / (data[i] + bg_dc[i][0] + bg_dc[i][1] + bg_dc[i][2]);

  return w;
}

#endif // USE_DC
   
} // end namespace
   
