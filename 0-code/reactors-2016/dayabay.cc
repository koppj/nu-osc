#include "definitions.h"

namespace ns_reactor
{

#ifdef USE_DB_3F

// source: Daya Bay talk at Neutrino 2016 - 1230 days
// using data from spectrum on slide 11 
  
#define NREACT_DB 6
#define NDET_DB 8
#define NBG_DB 5

#define DB_BIN_WIDTH 0.2
inline double T_min(const int b){
  return (b == 0 ? 0.7 : 1.3 + DB_BIN_WIDTH * (b-1));
}
inline double T_max(const int b){
  return 1.3 + DB_BIN_WIDTH * b;
}

extern Fit fit;
extern Rate_pull_coupl rate_pull_coupl;  // defined in class_flux.cc
extern Rate_coef rate;                   // defined in class_flux.cc
extern int old_new_main;

double gl_DmqL_db;
int bin_gl_db;
double DB_data_ND[NBIN_DB], DB_accidentals[NBIN_DB], DB_other_backgr[NBIN_DB];
double DB_data[NBIN_DB], DB_noosc[NBIN_DB];
double DB_FD_weight_noosc, DB_ND_weight_noosc[NBIN_DB];

/************** data **************************/


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


// ***************** DB data ******************
// *** data from slide 12 of NuFact2016 talk ***
// *********************************************

#ifdef DB_CHECK_NORM 
const double IBD_cand_DB[NDET_DB] = {
  597618., 606351., 567196., 466013., 80479., 80742., 80067., 66862.
};
#endif
  
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
 
// bg rate per day
const double bg_DB_data[NBG_DB][NDET_DB] = {
  {8.46, 8.46, 6.29, 6.18, 1.27, 1.19, 1.20, 0.98},  // accidentals
  {0.79, 0.79, 0.57, 0.57, 0.05, 0.05, 0.05, 0.05},  // fast neutron
  {2.46, 2.46, 1.72, 1.72, 0.15, 0.15, 0.15, 0.15},  // Li, He
  {0.15, 0.16, 0.13, 0.15, 0.04, 0.03, 0.03, 0.05},  // Am-C - using 8-AD values
  {0.08, 0.07, 0.05, 0.07, 0.05, 0.05, 0.05, 0.05}   // CO
};
// bg error per day
const double bg_err_DB[NBG_DB][NDET_DB] = {
  {0.09, 0.09, 0.06, 0.06, 0.01, 0.01, 0.01, 0.01}, // acc
  {0.10, 0.10, 0.07, 0.07, 0.01, 0.01, 0.01, 0.01}, // FN
  {1.06, 1.06, 0.77, 0.77, 0.06, 0.06, 0.06, 0.06}, //Li/He
  {0.07, 0.07, 0.06, 0.07, 0.02, 0.02, 0.02, 0.02}, // Am-C
  {0.04, 0.04, 0.03, 0.04, 0.03, 0.03, 0.03, 0.03}  // CO
}; 



/************************************************************/

class DB_PROB{
 private:
  double log_DmqL_min, log_DmqL_max;
  double DmqL_min, DmqL_max;
  double Del;

  double sinq_av[NBIN_DB][N_RATE_COEF];
  double sinq(const double DmqL, const int bin); // interpolating on table

 public:
  void calc_integrals(void);
  double P(Param_5nu &prm, int det, int react, int bin);

} DB_prob;


double DB_PROB::P(Param_5nu &p, int det, int react, int bin)
{
  double P = 1.;
  for(int i = 0; i < N_NU-1; i++){
    for(int j = i+1; j < N_NU; j++){
      const double DmqL = fabs(p.Dmq(j,i)) * basel_DB[det][react] * 1.e3;
      if(DmqL > DmqL_min)
        P -= 4. * norm(p.Ue[i] * p.Ue[j]) * sinq(DmqL, bin);
    }
  }
  return P; 
} 

// re-scaling the no-osc FD spectrum with oscillation probability
// do not include background in prediction to match official results
// there seems to be an inconsistency between the data in the table on p12 and the figure on p11
//double pred_FD(Param_5nu &prm, int bin)
//{
//  // FD weights 
//  double w_osc = 0.;
//
//  for(int i = 4; i < NDET_DB; i++){
//    for(int j = 0; j < NREACT_DB; j++){
//      w_osc += eff_DB[i] * DAQ_DB[i] * DB_prob.P(prm, i, j, bin) / norm(basel_DB[i][j]); 
//    }
//  }
//  //const double b = DB_accidentals[bin] + DB_other_backgr[bin];
//  //return (DB_noosc[bin] - b) * w_osc / DB_FD_weight_noosc  + b;
//  return DB_noosc[bin] * w_osc / DB_FD_weight_noosc;
//}


// extrapolating the ND spectrum to the FD using oscillations in both  
// used for Neutrino14 data
double pred_FD(Param_5nu &prm, int bin)
{
  // ND and FD weights 
  double w_nd = 0., w_fd = 0.;

  for(int i = 0; i < 4; i++){
    for(int j = 0; j < NREACT_DB; j++){
      w_nd += eff_DB[i] * DAQ_DB[i] * DB_prob.P(prm, i, j, bin) / norm(basel_DB[i][j]);

      const int k = i+4;
      w_fd += eff_DB[k] * DAQ_DB[k] * DB_prob.P(prm, k, j, bin) / norm(basel_DB[k][j]);
    }
  }
//  double b = DB_accidentals[bin] + DB_other_backgr[bin];
  return DB_noosc[bin] * w_fd / DB_FD_weight_noosc * DB_ND_weight_noosc[bin]/w_nd;
}

/**************** calc the table for class fit ***************/

void set_table_DB(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1])
{
  const int fp = fit.first_pull[DB];
  
  for(int i = 0; i < NBIN_DB; i++){
    const int b = fit.first_bin[DB] + i;
  
    // spectrum includes background
    cff[b][NPULLS] = pred_FD(prm, i);

    // detector norm
    cff[b][fp] = cff[b][NPULLS];

    // accidentals
    cff[b][fp+1] = DB_accidentals[i];
  }

  // energy scale
  for(int i = 0; i < NBIN_DB; i++){
    const int b = fit.first_bin[DB] + i;

    const double n_i = cff[b][NPULLS];

    if(i < NBIN_DB-1){
      const double n_ip1 = cff[b+1][NPULLS];		
      cff[b][fp+2] = (T_max(i) * n_ip1 - T_min(i) * n_i) / (T_max(i) - T_min(i));
	
    }else
      cff[b][fp+2] = cff[b-1][fp+2]; 
  }
  return;
}



/***************** init *************************************/

// accidental shape of FD from fig. 2 of 1310.6732
double int_background_db(double e){
  return read(REACTOR_PATH"Data_DB-Neutrino16/DB_accidentals_shape.txt", e);
}  
// non-accidental background - read from yellow region in the insert-plot on p11 of Nu16 talk
double int_nonaccidental_db(double e){
  return read(REACTOR_PATH"Data_DB-Neutrino16/non-accidental-backgr.csv", e);
}
    
  
void DB_init(void)
{
  // reading data
  FILE *fp_data = fopen(REACTOR_PATH"Data_DB-Neutrino16/data.csv","r");
  FILE *fp_noos = fopen(REACTOR_PATH"Data_DB-Neutrino16/no-osc.csv","r");

  for(int i = 0; i < NBIN_DB; i++){
    if(fscanf(fp_data, "%*e %le", &DB_data[i]) != 1)
      error("cannot read DB data");
    if(fscanf(fp_noos, "%*e %le", &DB_noosc[i]) != 1)
      error("cannot read DB no oscillation data");
  }
  fclose(fp_data);
  fclose(fp_noos);

  // scale to number of events per bin (factor 0.4% to adjust event numbers)
  // this reproduces 308150 events according to table on p12
  for(int i = 0; i < NBIN_DB; i++){
    DB_data[i]  *= 1.004e5 * (i == 0 ? 0.6 : DB_BIN_WIDTH);   
    DB_noosc[i] *= 1.004e5 * (i == 0 ? 0.6 : DB_BIN_WIDTH);   
  }
  
  // FD weights for no oscillations (ND weights not needed for Neutrino16 data fit)
  double w_nd = 0.;
  DB_FD_weight_noosc = 0.;
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < NREACT_DB; j++){
      w_nd += eff_DB[i] * DAQ_DB[i] / norm(basel_DB[i][j]);

      const int k = i+4;
      DB_FD_weight_noosc += eff_DB[k] * DAQ_DB[k] / norm(basel_DB[k][j]);
    }
  }
  
  // setting the backgrounds

  // accidental background
  // normalizing the shape from fig. 2 of 1310.6732
  DB_accidentals[0] = qromb1(int_background_db, 0.75, T_max(0), 1.e-5);
  for(int i = 1; i < NBIN_DB; i++)
    DB_accidentals[i] = qromb1(int_background_db, T_min(i), T_max(i), 1.e-5);

  double sum_bg = 0.;
  for(int i = 0; i < NBIN_DB; i++)
    sum_bg += DB_accidentals[i];
  for(int i = 0; i < NBIN_DB; i++)
    DB_accidentals[i] /= sum_bg;

  // expected number of accidentals in the 4 FD
  sum_bg = 0.;
  double sum_uncert_sq = 0.;
  for(int i = 4; i < NDET_DB; i++){
    sum_bg += bg_DB_data[0][i] * DAQ_DB[i] * eff_DB[i];
    sum_uncert_sq += norm(bg_err_DB[0][i] * DAQ_DB[i] * eff_DB[i]);
  }
  for(int i = 0; i < NBIN_DB; i++)
    DB_accidentals[i] *= sum_bg;

  // non-accidental background - read from yellow region in the insert-plot on p11 of Nu16 talk

  double sum_non_acc = 0.;
  for(int i = 0; i < NBIN_DB; i++)
    sum_non_acc += DB_other_backgr[i] = qromb1(int_nonaccidental_db, T_min(i), T_max(i), 1.e-5);
  
  // fprintf(stderr, "sum_non_acc = %f\n", sum_non_acc); // gives only 215 events!

  double norm_from_table = 0.;
  for(int j = 1; j < NBG_DB; j++)
    for(int k = 4; k < NDET_DB; k++)
      norm_from_table += bg_DB_data[j][k] * DAQ_DB[k] * eff_DB[k];
  for(int i = 0; i < NBIN_DB; i++)
    DB_other_backgr[i] *= norm_from_table / sum_non_acc;
  

#ifdef DB_CHECK_NORM 
  // check data points
  double sum = 0.;
  for(int i = 0; i < NBIN_DB; i++)
    sum += DB_data[i];
  printf("\nsum spectrum = %f\n", sum); 

  double ev = 0.;
  sum = 0.;

  // check rate
  for(int i = 0; i < NDET_DB; i++){
    ev = IBD_cand_DB[i];
    if(i >= 4) sum += ev; // FD events
    for(int j = 0; j < NBG_DB; j++)
      ev -= bg_DB_data[j][i] * DAQ_DB[i] * eff_DB[i];
    
    
    printf("rate det %d = %f\n", i, ev / (DAQ_DB[i] * eff_DB[i]));
  }
  printf("sum FD candidates = %f\n", sum);

  // backgrounds
  double sum_oth = 0.;
  for(int i = 0; i < NBIN_DB; i++)
    sum_oth += DB_other_backgr[i];
  
  printf("accidentals = %f, other backgr = %f, B/S = %f\n", sum_bg, sum_oth, (sum_bg+sum_oth)/sum);
  exit(0);
#endif // end check norm

  
  // set data in far and near detectors
  for(int i = 0; i < NBIN_DB; i++){

    // using the no osc spectrum in FD to infer ND data per bin  (not used for Neutrino16 data fit!)   
    DB_data_ND[i] = (DB_noosc[i] - DB_accidentals[i] - DB_other_backgr[i]) * w_nd / DB_FD_weight_noosc;

    //DB_data_ND[i] *= 0.999; // introduced "by hand" to match th13 best fit point for Neutrino14 fit

    // set data and stat. errors in class fit
    const int b = fit.first_bin[DB] + i;  
    fit.Data[b] = fit.S_data[b][b] = DB_data[i]; 
  }

  // *** set pull errors in class fit ***
  int p = fit.first_pull[DB];

  // relative detector normalization: each detector 0.13%  - slide 9 of Nu16 talk
  // for 4 ND and 4 FD: 0.13% / sqrt(4)   
  // contribution of ND + FD: sqrt(2) ---> 0.0013 / 2 * sqrt(2)
  fit.S_pull[p][p] = norm(0.0027 / M_SQRT2);  // changed 0.0013 to 0.0027 by hand to get theta13 accuracy right
  p++;

  // accidental background
  fit.S_pull[p][p] = sum_uncert_sq / norm(sum_bg); 
  p++;

  // energy scale 
  fit.S_pull[p][p] = norm(0.0035); // relative energy scale uncert. 1310.6732 

  DB_prob.calc_integrals();

  // JK - compute FD weights. We do this only here as otherwise the data
  // JK - required for the oscillation probabilities would not be available yet
  Param_5nu db_bf;
  db_bf.theta[0] = NAN;
  db_bf.theta[I12] = asin(sqrt(0.31));
  db_bf.theta[I13] = 0.5*asin(sqrt(0.084));
  db_bf.theta[I14] = 0.0;
  db_bf.theta[I15] = 0.0;
  db_bf.set_ang();
  db_bf.dmq[0]     = 0.;
  db_bf.dmq[1]     = 7.59e-5;
  db_bf.dmq[2]     = 2.5e-3;
  db_bf.dmq[3]     = 1.;
  db_bf.dmq[4]     = 10.;
  for (int bin=0; bin < NBIN_DB; bin++){
    DB_ND_weight_noosc[bin] = 0.0;
    for(int i = 0; i < 4; i++){
      for(int j = 0; j < NREACT_DB; j++){
        DB_ND_weight_noosc[bin] += eff_DB[i] * DAQ_DB[i]
                                * DB_prob.P(db_bf, i, j, bin) / norm(basel_DB[i][j]);
      }
    }
  }
  return;  
}




// interpolating on the table
double DB_PROB::sinq(const double DmqL, const int bin)
{
  if(DmqL <  DmqL_min) error("DB_PROB::sinq: DmqL too small");
  if(DmqL >= DmqL_max) error("DB_PROB::sinq: DmqL too large");

  const double logD = log10(DmqL);

  int i = int( (logD - log_DmqL_min) / Del);
  if(i < 0 || i >= N_RATE_COEF-1) error("DB_PROB::sinq: index out of range");

  const double x = (logD - i * Del - log_DmqL_min) / Del;
  if(x < 0. || x > 1.) error("DB_PROB::sinq: x out of range");

  return sinq_av[bin][i] + (sinq_av[bin][i+1] - sinq_av[bin][i]) * x;
}


/****************************************************
 * plotting
 ****************************************************/

void DB_spectrum(Param_5nu &prm)
{
  double cff[NBIN_CHISQ][NPULLS+1];
  set_table_DB(prm, cff);
  
  for(int i = 0; i < NBIN_DB; i++){
    const double f = 1.e-5 / (T_max(i) - T_min(i));
    printf("%f %f %f %f %f\n", T_min(i), f*DB_noosc[i], f*DB_data[i],
	   f*cff[fit.first_bin[DB] + i][NPULLS], f*(DB_accidentals[i] + DB_other_backgr[i]));    
    printf("%f %f %f %f %f\n", T_max(i), f*DB_noosc[i], f*DB_data[i],
	   f*cff[fit.first_bin[DB] + i][NPULLS], f*(DB_accidentals[i] + DB_other_backgr[i]));    
  }
  return;
}
  
  

/****************************************************
 * integration routines
 ****************************************************/

// 1310.6732 gives sigma/E = 8% @ 1 MeV - assume sqrt(E) scaling
inline double e_res_db(const double T_e){ 
  return 0.08 * sqrt(T_e);
}
double gauss_int_db(double eNu)
{   
   double T_e = eNu + ME - DELTA;

   // take into account energy scale correction according to Fig. 1 of 1310.6732
   //T_e *= read(REACTOR_PATH"Data_DB-Neutrino14/energy-scale.txt", T_e);

   const double sigma = e_res_db(T_e);

   double x_max = (T_max(bin_gl_db) - T_e) / (sigma*M_SQRT2);
   double x_min = (T_min(bin_gl_db) - T_e) / (sigma*M_SQRT2);

   double res = (erf(x_max) - erf(x_min)) / 2.;
   if(res < 0.) {
      fprintf(stderr, "BLOODY HELL!!\n");
      exit(1);
   }   
   return res;
}

double db_flux(const double Enu)
{
  // slide 19 of Neutrino2014 talk
  const double isofract_DB[NISO] = {0.586, 0.076, 0.288, 0.050};

  double w = 0.;
  for(int i = 0; i < NISO; i++)
    w += isofract_DB[i] * global_flux(i, Enu, old_new_main);
  return w;
}

double int_n0_db(double e)
{
  const double EposKin = e - ME - DELTA;
  return crossSect(EposKin) * db_flux(e) * gauss_int_db(e);
}

double int_sinq_db(double e)
{
  double Cos;
  if(gl_DmqL_db <= 0.) Cos = 1.;
  else
  {
    const double a1 = 2.54 * gl_DmqL_db / e * 0.95;
    const double a2 = 2.54 * gl_DmqL_db / e * 1.05;
    Cos = (sin(a2) - sin(a1))/(a2 - a1);
  }

  const double EposKin = e - ME - DELTA;
  return crossSect(EposKin) * 0.5*(1.-Cos) * db_flux(e) * gauss_int_db(e);
}



void DB_PROB::calc_integrals(void)
{
  // DmqL in units of eV^2 m 
  log_DmqL_min = ATM_MIN + log10(10.);      
  log_DmqL_max = STE_MAX + log10(2000.);   

  DmqL_min = exp10(log_DmqL_min);
  DmqL_max = exp10(log_DmqL_max);

  Del = (log_DmqL_max - log_DmqL_min) / (N_RATE_COEF - 1.);

  for(int i = 0; i < NBIN_DB; i++){
    bin_gl_db = i;

    double emin = T_min(i) - 3. * e_res_db(T_min(i)) - ME + DELTA;
    if(emin < EnuMIN) emin = EnuMIN;
    double emax = T_max(i) + 3. * e_res_db(T_max(i)) - ME + DELTA;
    if(emax > EnuMAX) emax = EnuMAX;

    const double n0 = qromb1(int_n0_db, emin, emax, 1.e-6);

    for(int j = 0; j  < N_RATE_COEF; j++){

      gl_DmqL_db = exp10(log_DmqL_min + j * Del);
      sinq_av[i][j] = qromb1(int_sinq_db, emin, emax, 1.e-6) / n0;
    }
  }

  /*
  // printing the averaged sinq factor 
  for(double L = 100.; L < 2000.; L *= 1.01)
    printf("%e %e %e %e %e\n", L, sinq(7.5e-5 * L, 1), sinq(7.5e-5 * L, 35), 
	  sinq(2.3e-3 * L, 1), sinq(2.3e-3 * L, 35));
  exit(0);
  */
  return;
}


#endif

} // namespace
