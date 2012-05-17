#include "definitions.h"

/************************************
 *  ATTENTION L in km and E in GeV
 ************************************/

#define DIST (0.541) // target -> center of det.

// 50 m long decay tunnel, 10 m fiducial det. diameter
#define LU (DIST - 0.025)
#define LO (DIST + 0.005)

#define N_BINS_M 10
#define BIN_0 2     // lowest bin used in chisq

struct bin{
  double lw;
  double up;
} bins[N_BINS_M] = {
  {0.30,  0.375},
  {0.375, 0.475},
  {0.475, 0.55},
  {0.55,  0.675},
  {0.675, 0.8},
  {0.8,   0.95},
  {0.95,  1.1},
  {1.1,   1.3},     
  {1.3,   1.5},
  {1.5,   3.}
};

#define E_MIN 0.01
#define E_MAX 5.

inline double E_min_MB(const int i){
  return bins[i].lw;
}
inline double E_max_MB(const int i){
  return bins[i].up;
}
inline double DE_MB(const int i){
  return bins[i].up - bins[i].lw;
}


/* the data  1st col: data (excess) 
 *           2nd col: data+1sigma
 *           3rd col: nu_m background
 *           4th col: nu_e background
 *           5th col: data (events)
 */
const double MB_data[N_BINS_M][5] = {
  {0.710280,    0.943925, 1.544642, 0.400900,  2.651429},
  {0.429907,    0.593458, 0.896612, 0.365711,  1.704490},
  {0.163551,    0.313084, 0.531138, 0.424918,  1.112653},
  {0.070093,    0.182243, 0.259808, 0.389719,  0.722041},
  {0.009346,    0.093458, 0.189245, 0.330916,  0.508980},
  {-0.01   ,    0.074766, 0.130116, 0.272035,  0.390612},
  {-0.01   ,    0.074766, 0.071027, 0.260396,  0.331429},
  {0.004673,    0.070093, 0.047401, 0.177710,  0.236735},
  {0.004673,    0.065421, 0.047453, 0.142312,  0.177551},
  {0.000000,    0.084112, 0.035630, 0.023751,  0.047347}
};

// global variables
double cos_sin_MB[DM2ANZ][2][N_BINS_M];
double events_full[N_BINS_M];
double **s_inv_300, **s_inv_475;


/*********** flux and cross section *********************/


inline double nue_flux_mb(const double Enu)
{
  if(Enu <= 0.0) return 0.0;
  return pow(10., read(PATH"/Data-files/mbooneflux.dat", Enu));
}

inline double XS_mb(const double Enu){
  return read(PATH"/Data-files/XQE-nu_e.dat", log10(Enu)) * Enu;
}

/********************  the chisq *********************************/

void MB_events(params p, double no[N_BINS_M])
{
  for(int i = 0; i < N_BINS_M; i++){

    // interpolating on the table
    double Icos[3], Isin[3];

    for(int j = 0; j < 2; j++){  // j = 41, 51

      const double dmq = p.dmq[j];

#ifdef CHECK
      if(dmq <= 0.){
        Icos[j] = cos_sin_MB[0][COS][i];
        Isin[j] = 0.;
      }else
#endif
      {

        const double logD = log10(dmq);

	const int k = (logD < LDM2MIN ? 0 : int( (logD - LDM2MIN)/DM_STEP ) + 1);

	Icos[j] = cos_sin_MB[k][COS][i] + 
          (cos_sin_MB[k+1][COS][i] - cos_sin_MB[k][COS][i])/
	  (dm2k(k+1) - dm2k(k)) * (dmq - dm2k(k));

	Isin[j] = cos_sin_MB[k][SIN][i] + 
          (cos_sin_MB[k+1][SIN][i] - cos_sin_MB[k][SIN][i])/
	  (dm2k(k+1) - dm2k(k)) * (dmq - dm2k(k));
      }
    } // end for j 

    // 54
#ifndef Ip3pI
    const double dmq = fabs(p.dmq[I5] - p.dmq[I4]);
    const double logD = log10(dmq);
#else
    const double dmq = fabs(p.dmq[I5] + p.dmq[I4]);
    double logD = log10(dmq);
    if (logD > LDM2MAX-1.00001*DM_STEP) logD = LDM2MAX-1.00001*DM_STEP; // JK
#endif

    if(dmq == 0.){
        Icos[2] = cos_sin_MB[0][COS][i];
        Isin[2] = 0.;
    }else{

	const int k = (logD < LDM2MIN ? 0 : int( (logD - LDM2MIN)/DM_STEP ) + 1);

	Icos[2] = cos_sin_MB[k][COS][i] + 
          (cos_sin_MB[k+1][COS][i] - cos_sin_MB[k][COS][i])/
	  (dm2k(k+1) - dm2k(k)) * (dmq - dm2k(k));

	Isin[2] = cos_sin_MB[k][SIN][i] + 
          (cos_sin_MB[k+1][SIN][i] - cos_sin_MB[k][SIN][i])/
	  (dm2k(k+1) - dm2k(k)) * (dmq - dm2k(k));

#ifndef Ip3pI
	if(p.dmq[I5] - p.dmq[I4] < 0.)
	  Isin[2] *= -1.;
#endif
  
    } // end interpolation

    const double cosd = cos(p.delta);
    const double sind = sin(p.delta); 

    no[i] = cosd 
      + cosd * Icos[2]  + sind * Isin[2]   // 54 -term
      - cosd * Icos[I5] - sind * Isin[I5]  // 51 -term
#ifndef Ip3pI
      - cosd * Icos[I4] + sind * Isin[I4]; // 41 -term
#else
      - cosd * Icos[I4] - sind * Isin[I4]; // 41 -term
#endif

    no[i] *= 2. * p.Ue[I4] * p.Um[I4] * p.Ue[I5] * p.Um[I5];

    no[i] += 2. * sqr(p.Ue[I4] * p.Um[I4]) * (1. - Icos[I4]) +  
             2. * sqr(p.Ue[I5] * p.Um[I5]) * (1. - Icos[I5]);

    // events per MeV: 
    // multiplying probability with events for full osc per bin width
    no[i] *= events_full[i] / (E_max_MB(i) - E_min_MB(i)) * 2000.;

  } // end for i (bins)
  return;
}

/********************  the chisq *********************************/

double chisq_minib(const int bin_0, params p)
{
  double no[N_BINS_M];
  MB_events(p, no);

  double s = 0.;
  for(int i = bin_0; i < N_BINS_M; i++)
    for(int j = bin_0; j < N_BINS_M; j++)
      s += (no[i] - MB_data[i][0]) * (no[j] - MB_data[j][0]) 	
	* (bin_0 ? s_inv_475[i+1][j+1] : s_inv_300[i+1][j+1]);
	  
  return s;
}


/********************* alternative chisq **********************

#define SIG_NORM 0.3
double chisq_minib(const int bin_0, params p)
{
  double no[N_BINS_M], err[N_BINS_M];
  MB_events(p, no); 
  
  double s1, s2;
  s1 = s2 = 1./sqr(SIG_NORM);
  
  for(int i = bin_0; i < N_BINS_M; i++){
    err[i] = 0.9 * (MB_data[i][1] - MB_data[i][0]);
    
    s1 += sqr(no[i] / err[i]);
    s2 += no[i] * MB_data[i][0] / sqr(err[i]);
  }
  
  const double a = s2/s1;
  
  s1 = sqr( (a - 1.)/SIG_NORM );
  
  for(int i = bin_0; i < N_BINS_M; i++)
    s1 += sqr( (a * no[i] - MB_data[i][0]) / err[i] );
	  
  return s1;
}
*/

/******************** calc integrals ************************/

inline double eff_mb(const double e)
{
  return read(PATH"/Data-files/minib-efficiency.dat", e) * 1.4;
}

// energy resolution from MB talk Fermilab 11/04/07 slide 24
#define B 0.08
#define C 0.024

inline double sigma_mb(const double e){
  return sqrt(B * B * e + C * C);
}

// global vars
int bin_gl_MB, sin_cos_M;
double dm2_gl_MB;

double intFuncMB(double Enu)
{
  double i1 = erf((E_max_MB(bin_gl_MB) - Enu)/(M_SQRT2 * sigma_mb(Enu))); 
  i1 -= erf((E_min_MB(bin_gl_MB) - Enu)/(M_SQRT2 * sigma_mb(Enu))); 
  
  const double a = 2.53 * dm2_gl_MB / Enu;
  
  double sin_cos_mean;

  if(sin_cos_M == SIN)
  {
    if(a <= 0.)
      nrerror("[intFuncMB]: a should not be 0 in case of SIN");
    else
      sin_cos_mean = (cos(a*LU) - cos(a*LO))/(a * (LO-LU));
  }
  else
    sin_cos_mean = 
      (a <= 0. ? 1. : (sin(a*LO) - sin(a*LU))/(a * (LO-LU)) );
  
  return nue_flux_mb(Enu) * XS_mb(Enu) * i1 * (sin_cos_mean + 1.) * eff_mb(Enu);
}

#define N_SIGM_MB 4.

void initMiniboone(void)
{
  for(int k = 0; k < DM2ANZ; k++){
    
    dm2_gl_MB = dm2k(k);
    
    for(bin_gl_MB = 0; bin_gl_MB < N_BINS_M; bin_gl_MB++){
      
      double E_lo = E_min_MB(bin_gl_MB) - sigma_mb(E_min_MB(bin_gl_MB)) * N_SIGM_MB;
      if(E_lo < E_MIN)
	E_lo = E_MIN;
      double E_up = E_max_MB(bin_gl_MB) + sigma_mb(E_max_MB(bin_gl_MB)) * N_SIGM_MB;
      if(E_up > E_MAX)
	E_up = E_MAX;
      
      for(sin_cos_M = 0; sin_cos_M < 2; sin_cos_M++){

        if(k == 0 && sin_cos_M == SIN)
	  cos_sin_MB[k][SIN][bin_gl_MB] = 0.;
        else
          cos_sin_MB[k][sin_cos_M][bin_gl_MB] = qromb1(intFuncMB, E_lo, E_up, 1.0e-5);
      
        // correcting for the " + 1" 
        if(k == 0)
	  cos_sin_MB[k][sin_cos_M][bin_gl_MB] *= 0.5;
        else
	  cos_sin_MB[k][sin_cos_M][bin_gl_MB] -= cos_sin_MB[0][COS][bin_gl_MB];
      }
    }
  }    

  // normalization
  for(bin_gl_MB = 0; bin_gl_MB < N_BINS_M; bin_gl_MB++){
    for(int k = 1; k < DM2ANZ; k++){
      for(sin_cos_M = 0; sin_cos_M < 2; sin_cos_M++){	
	
	 cos_sin_MB[k][sin_cos_M][bin_gl_MB] /= cos_sin_MB[0][COS][bin_gl_MB];
      }
    }
    events_full[bin_gl_MB] = cos_sin_MB[0][COS][bin_gl_MB];
    cos_sin_MB[0][COS][bin_gl_MB] = 1.;
  }
  
  // init covariance matrix
  int i, j, *indx;
  double **s_300, **s_475, d, *col;

  s_300 = matrix(1, N_BINS_M, 1, N_BINS_M);
  s_475 = matrix(1, N_BINS_M, 1, N_BINS_M);
  s_inv_300 = matrix(1, N_BINS_M, 1, N_BINS_M);
  s_inv_475 = matrix(1, N_BINS_M, 1, N_BINS_M);
  col = vector(1, N_BINS_M);
  indx = ivector(1, N_BINS_M);
  
  double err[N_BINS_M+1];
  for(i = 0; i < N_BINS_M; i++)
    err[i+1] = MB_data[i][0] - MB_data[i][1];

  const double err_bg  = sqr(0.);
  const double err_cor = sqr(0.07);

  // covar matrices 
  for(i = 1; i <= N_BINS_M; i++){
    for(j = 1; j <= N_BINS_M; j++){
 
      s_475[i][j] = (i > 2 && j > 2 ? 
		 err_bg * (MB_data[i-1][2]*MB_data[j-1][2] + 
                           MB_data[i-1][3]*MB_data[j-1][3])+
                 err_cor * MB_data[i-1][4]*MB_data[j-1][4] : 0.);

      s_300[i][j] = 
		 err_bg * (MB_data[i-1][2]*MB_data[j-1][2] + 
                           MB_data[i-1][3]*MB_data[j-1][3])+
                 err_cor * MB_data[i-1][4]*MB_data[j-1][4];
    }    
    const double stretch_err = 1.1;
    if(s_300[i][i] > sqr(err[i] * stretch_err)){
      fprintf(stderr, "error %d too large: %f\n", i, s_300[i][i] / sqr(err[i] * stretch_err));
      //exit(0);
    }
    s_300[i][i] = sqr(err[i] * stretch_err);
    s_475[i][i] = sqr(err[i] * stretch_err);
  }
  
  // inverting the covar matrix for MB475
  ludcmp(s_475, N_BINS_M, indx, &d);
  
  for(j = 1; j <= N_BINS_M; j++){
    
    for(i = 1; i <= N_BINS_M; i++) col[i] = 0.0;
    
    col[j] = 1.0;
    lubksb(s_475, N_BINS_M, indx, col);
    
    for(i = 1; i <= N_BINS_M; i++) s_inv_475[i][j] = col[i];
  }
  // inverting the covar matrix for MB300
  ludcmp(s_300, N_BINS_M, indx, &d);
  
  for(j = 1; j <= N_BINS_M; j++){
    
    for(i = 1; i <= N_BINS_M; i++) col[i] = 0.0;
    
    col[j] = 1.0;
    lubksb(s_300, N_BINS_M, indx, col);
    
    for(i = 1; i <= N_BINS_M; i++) s_inv_300[i][j] = col[i];
  }
  
  return;
}


void plot_excess_mb(params p)
{   
  double no[N_BINS_M];
  MB_events(p, no);

  for(int i = 0; i < N_BINS_M; i++){
	
     printf("%e %e\n", E_min_MB(i), no[i]);
     printf("%e %e\n", E_max_MB(i), no[i]);
  }   
  exit(0);
}

void plot_excess_mb(char *name, params p)
{   
  double no[N_BINS_M];
  MB_events(p, no);
  
  FILE *fp = fopen(name, "w");

  for(int i = 0; i < N_BINS_M; i++){
	
     fprintf(fp, "%e %e\n", E_min_MB(i), no[i]);
     fprintf(fp, "%e %e\n", E_max_MB(i), no[i]);
  }
  fclose(fp);
  return;
}

/********************************************
 *   testing miniboone
 ********************************************/

#ifdef MINIBOONE

int main(void)
{
  initMiniboone();
  params p;

  p.Ue[0] = 0.;
  p.Um[0] = 0.;
  p.delta = 0.;
  p.dmq[0] = 0.;

  double min = 1.e4;
  params p_min;

  //for(int i = 0; i < 101; i++)
  int i = 20;
  {
    const double sq = pow(10., -4. + 4. * i/100.);
    p.Ue[1] = 0.01;
    p.Um[1] = sqrt(sq) / (p.Ue[1] * 2.);

    for(int j = 0; j < 101; j++){
      p.dmq[1] = pow(10., -2. + 4. * j/100.);

      printf("%e  %e  %e\n", sq, p.dmq[1], chi2mb475(p));
      /*
      double w = chi2mb300(p);
      if(w < min){
	min = w;
	p_min = p;
      }
       */ 
    }
  }
  exit(0);
  /*    
  double no[N_BINS_M];

  MB_events(p_min, no);
 
  for(int i = 0; i < N_BINS_M; i++){
    printf("%e  %e \n", E_min_MB(i), no[i]);
    printf("%e  %e \n", E_max_MB(i), no[i]);
  }

  return 0;
  */

  double no1[N_BINS_M], no2[N_BINS_M];

  /*
  double sq = 0.2;
  p.Ue[1] = 0.01;
  p.Um[1] = sqrt(sq) / (p.Ue[1] * 2.);
  p.dmq[1] = 0.1;  

  sq = 0.004;
  p.Um[1] = sqrt(sq) / (p.Ue[1] * 2.);
  p.dmq[1] = 1.;  
  */


  //read_param("out-5nu-fine/all_mb475.point", &p);
  read_param("out-5nu/app_mb475.point", &p);
  print_param_pretty(stderr, p);

  MB_events(p, no1);

  read_param("out-5nu/app_mb300.point", &p);
  MB_events(p, no2);
  print_param_pretty(stderr, p);
 

  for(int i = 0; i < N_BINS_M; i++){
    printf("%e  %e %e\n", E_min_MB(i), no1[i], no2[i]);
    printf("%e  %e %e\n", E_max_MB(i), no1[i], no2[i]);
  }
  
  return 0;
}
  
#endif

#undef E_MIN
#undef E_MAX
#undef DIST 
#undef LU 
#undef LO 
#undef N_BINS_M 
#undef BIN_0 
#undef B 
#undef C 
#undef N_SIGM_MB
