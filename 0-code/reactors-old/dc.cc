#include "definitions.h"

namespace ns_reactor
{

#ifdef USE_DC

#define ACCURACY 1.e-5

// 0.35 gives good agreement with DC no-osc MC 
#define E_RESOL 0.35

inline double e_res_dc(double e)
{
  //return 0.4;
  //return 0.15 + 0.12*sqrt(e);
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

// total pred for no osc including bg, and background
const double pred_tot_dc[NBIN_DC][2] = {
  {122.377625, 18.4926186},    
  {308.935516, 18.4926186},
  {466.122772, 18.492618},
  {562.937073, 19.03651},
  {571.639465, 13.597513},
  {530.30304,  15.22921},
  {447.086243, 16.31701},
  {362.781677, 16.86091},
  {288.267303, 16.86091},
  {212.121216, 17.40481},
  {151.748245, 17.40481},
  {107.148407, 17.40481},
  {64.1802673, 16.31701},
  {39.1608391, 15.77311},
  {22.8438225, 14.14141},
  {26.1072254, 22.8438225},
  {17.4048176, 16.3170166},
  {17.4048176, 17.4048176}
};

// events per bin read off from plot
double data_dc[NBIN_DC] = {
  119.658119,
  304.040405,
  449.805756,
  513.986023,
  596.658875,
  478.088593,
  403.030304,
  335.586639,
  298.057495,
  207.770004,
  168.065262,
  82.6728821,
  54.9339561,
  23.9316235,
  25.0194244,
  28.2828274,
  13.0536127,
  15.2292156
};


// backgrounds read off from figure of DC 1st pub in ev / 0.5 MeV
const double FN_BG = 3.348;

// sum of fast n and acc
const double FN_ACC_BG[NBIN_DC] = {17.828, 12.127, 7.692, 6.878, 3.982, 
  3.529, 3.348,3.348,3.348,3.348,3.348,3.348,3.348,3.348,3.348,
  3.348,3.348,3.348};

// total BG
const double bg_dc_tot[NBIN_DC] = {24.434, 21.267, 16.199, 
  14.299, 13.213, 12.851, 13.937, 14.661,
  15.747, 14.751, 15.928, 15.385, 14.118, 12.851, 10.588, 
  7.873, 4.525, 4.525};


// bg rate per day and error [DC 1st pub]
const double bg_dc_rate[NBG_DC][2] = {
  {0.332, 0.03}, // acc
  {2.3,   1.2},   // Li9
  {0.83,  0.38}    // fast n
};


/*
double accidental(double e){
  return read("Data_DC/accidentals.dat", e);
}
*/

const double rate_nu = 42.6;  // neutrino rate per day


/**************** calc the table for class fit ***************/

void set_table_dc(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1])
{
  // for interpolation
  const double logDmq = log10(prm.dmq[2]);
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


    // scale to Bugey4 measurement
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

/***************** init *************************************/

void dc_init(const int old_new)
{
  // calc coefficients
  old_new_gl_dc = old_new;
  calc_coef_dc();

  // total number of neutrinos
  const double tot_ev_obs = 4121.;

  // normalize the data read off from the plot 
  double n_nu_dc = 0.; 

  // correct last bin
  data_dc[NBIN_DC-1] *= 5.5/4.;
  
  for(int i = 0; i < NBIN_DC; i++)
    n_nu_dc += data_dc[i];
  //fprintf(stderr, "DC neutrinos = %f\n", n_nu_dc);
  
  for(int i = 0; i < NBIN_DC; i++)
    data_dc[i] *= tot_ev_obs / n_nu_dc;

  // set data and covariance matrix in class fit    
  for(int i = 0; i < NBIN_DC; i++){
    const int ii = fit.first_bin[DC] + i;    
    fit.Data[ii] = fit.S_data[ii][ii] = data_dc[i];
  }

  /*** set the backgrounds ***/
  double sum_bg[NBG_DC] = {0., 0., 0.};

  for(int i = 0; i < NBIN_DC; i++){
    // accidental
    //bg_dc[i][0] = qromb1(accidental, bin[i][0], bin[i][1], ACCURACY);
    bg_dc[i][0] = (FN_ACC_BG[i] - FN_BG) * (bin[i][1] - bin[i][0])/0.5; 
    sum_bg[0] += bg_dc[i][0];

    // Li9
    bg_dc[i][1] = (bg_dc_tot[i] - FN_ACC_BG[i])*(bin[i][1] - bin[i][0])/0.5; 
    sum_bg[1] += bg_dc[i][1];

    // fast n: flat 
    bg_dc[i][2] = bin[i][1] - bin[i][0];
    sum_bg[2] += bg_dc[i][2];
  }
  
  // normalize the backgrounds relative to neutrino signal
  for(int i = 0; i < NBIN_DC; i++){
    for(int j = 0; j < NBG_DC; j++){
      bg_dc[i][j] *= bg_dc_rate[j][0] / sum_bg[j] / rate_nu * tot_ev_obs;
    }
    // print backgrounds
    //printf("%f %f %f %f\n", bg_dc[i][0], bg_dc[i][1], bg_dc[i][2], 
    //bg_dc[i][0]+bg_dc[i][1]+bg_dc[i][2]);    
  }

  double tot_bg = 0.;
  for(int j = 0; j < NBG_DC; j++)
    tot_bg += bg_dc_rate[j][0] / rate_nu * tot_ev_obs;
  //fprintf(stderr, "DC tot background = %f\n", tot_bg);
  

  
  /*** set the pull errors ***/

  // DC normalization
  fit.S_pull[fit.first_pull[DC]][fit.first_pull[DC]] = norm(0.021);  // only detector norm
  //fit.S_pull[fit.first_pull[DC]][fit.first_pull[DC]] = norm(0.021) + norm(0.018);  // detector + reactor

  // backgrounds
  for(int j = 0; j < NBG_DC; j++){
    const int p = fit.first_pull[DC] + 1 + j;
    fit.S_pull[p][p] = norm(bg_dc_rate[j][1] / bg_dc_rate[j][0]); 
  }

  // energy scale
  fit.S_pull[fit.first_pull[DC] + 4][fit.first_pull[DC] + 4] = norm(0.015);
  fit.pull_status[fit.first_pull[DC]+4] = FIXED; // energy scale error not used

    
  /**** the no oscillation prediction ****/
    
#ifdef USE_DC_PRED 
  // subtracting background, scaling to events per bin
  for(int i = 0; i < NBIN_DC; i++){
    dc_no_osc_react[i] = pred_tot_dc[i][0];
    
    for(int j = 0; j < NBG_DC; j++)
      dc_no_osc_react[i] -= bg_dc[i][j];
  }
  
#else // use my own prediction
  
  // expected # events incl bg normalized to Bugey4
  const double dc_expect_ev_tot = 4344.; // [DC 1st pub]
  double w = 0.;
  for(int i = 0; i < NBIN_DC; i++)
    w += co_dc.c0[i];

  for(int i = 0; i < NBIN_DC; i++)
    dc_no_osc_react[i] = co_dc.c0[i] * (dc_expect_ev_tot - tot_bg) / w;

#endif 
  
  /*
  for(int i = 0; i < NBIN_DC; i++){
    printf("%e %e %e\n", bin[i][0], co_dc.c0[i] * no_osc / w, dc_no_osc_react[i]);
    printf("%e %e %e\n", bin[i][1], co_dc.c0[i] * no_osc / w, dc_no_osc_react[i]);
  }
  exit(0);
  */
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

void plot_dc_data(void)
{
  for(int i = 0; i < NBIN_DC; i++)
    printf("%e %e %e\n", 0.5*(bin[i][0]+bin[i][1]), 
	   data_dc[i], sqrt(data_dc[i]));
  exit(0);
}

void plot_dc_pred(Parameters p, double f)
{
  fprintf(stderr, "dmq = %e\n", p.dmq31());
  
  for(int i = 0; i < NBIN_DC; i++){
    
    double w = (1. - p.sq2_t13() * co_dc.c1[i][p.idmq_31]);
    w *= dc_no_osc_react[i] * f  / sbl_data[BUGEY4][old_new_gl_dc];

    /*
    printf("%f %f %f %f %f\n",
	   bg_dc[i][2], 
	   bg_dc[i][0],
	   bg_dc[i][1],
	   w,
	   bg_dc[i][2]+bg_dc[i][0]+bg_dc[i][1] + w);
    /*/
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
