#include "definitions.h"

/************************************
 *  ATTENTION L in km and E in GeV
 ************************************/

#define DIST (0.541) // target -> center of det.

// 50 m long decay tunnel, 10 m fiducial det. diameter
#define LU (DIST - 0.025)
#define LO (DIST + 0.005)

#define N_BINS_MA 11
#define BIN_0 3     // lowest bin used in chisq

struct bin{
  double lw;
  double up;
} bins_mba[N_BINS_MA] = {
  {0.2,   0.3},
  {0.3,   0.375},
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

inline double E_mba_min(const int i){
  return bins_mba[i].lw;
}
inline double E_mba_max(const int i){
  return bins_mba[i].up;
}
inline double DE(const int i){
  return bins_mba[i].up - bins_mba[i].lw;
}



// global variables


/* the data  1st col: background 
 *           2nd col: data
 *           3rd col: 1sigma
 */
double MB_A_data[N_BINS_MA][3]; 

double cos_sin_MB_A[DM2ANZ][2][N_BINS_MA];
double **s_inv_200, **s_inv_475_mba;
double MB_anti_norm[N_BINS_MA];

#define REF_POINT3

#ifdef REF_POINT1
// prediction for sin^22theta = 0.004, Dmq = 1 eVq
const double MB_pred_anti[N_BINS_MA] = {
  0.020956,0.047426,0.059559,0.060662,0.052941,0.039706,
  0.026471,0.017647,0.009926, 0.005515,0.004412};
#endif

#ifdef REF_POINT2
// prediction for sin^22theta = 0.0061, Dmq = 4.42 eVq
const double MB_pred_anti[N_BINS_MA] = {
  0.030882,0.048529,0.044118,0.059559,0.082721,0.044118,
  0.017647,0.020956,0.028676,0.026471,0.001103};
#endif

#ifdef REF_POINT3
// prediction for sin^22theta = 0.2, Dmq = 0.1 eVq
const double MB_pred_anti[N_BINS_MA] = {
 0.091544,0.078309,0.063971,0.054044,0.038603,0.027574,
 0.017647,0.011029,0.005515,0.003309,0.0002};
#endif


/*********** flux and cross section *********************/


inline double nue_flux_anti(const double Enu)
{
  if(Enu <= 0.0) return 0.0;
  return pow(10., read(PATH"/Data-files/MB-anti-nu-flux.dat", Enu) - 5.);
}

inline double XS_anti(const double Enu){
  return read(PATH"/Data-files/XQE-anti-nu_e.dat", log10(Enu)) * Enu;
}

/********************  the chisq *********************************/

void MB_A_events(params p, double no[N_BINS_MA])
{
  for(int i = 0; i < N_BINS_MA; i++){

    // interpolating on the table
    double Icos[3], Isin[3];

    for(int j = 0; j < 2; j++){  // j = 41, 51

      const double dmq = p.dmq[j];

#ifdef CHECK
      if(dmq <= 0.){
        Icos[j] = cos_sin_MB_A[0][COS][i];
        Isin[j] = 0.;
      }else
#endif
      {

        const double logD = log10(dmq);

	const int k = (logD < LDM2MIN ? 0 : int( (logD - LDM2MIN)/DM_STEP ) + 1);

	Icos[j] = cos_sin_MB_A[k][COS][i] + 
          (cos_sin_MB_A[k+1][COS][i] - cos_sin_MB_A[k][COS][i])/
	  (dm2k(k+1) - dm2k(k)) * (dmq - dm2k(k));

	Isin[j] = cos_sin_MB_A[k][SIN][i] + 
          (cos_sin_MB_A[k+1][SIN][i] - cos_sin_MB_A[k][SIN][i])/
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
        Icos[2] = cos_sin_MB_A[0][COS][i];
        Isin[2] = 0.;
    }else{

	const int k = (logD < LDM2MIN ? 0 : int( (logD - LDM2MIN)/DM_STEP ) + 1);

	Icos[2] = cos_sin_MB_A[k][COS][i] + 
          (cos_sin_MB_A[k+1][COS][i] - cos_sin_MB_A[k][COS][i])/
	  (dm2k(k+1) - dm2k(k)) * (dmq - dm2k(k));

	Isin[2] = cos_sin_MB_A[k][SIN][i] + 
          (cos_sin_MB_A[k+1][SIN][i] - cos_sin_MB_A[k][SIN][i])/
	  (dm2k(k+1) - dm2k(k)) * (dmq - dm2k(k));

#ifndef Ip3pI
	if(p.dmq[I5] - p.dmq[I4] < 0.)
	  Isin[2] *= -1.;
#endif  
    } // end interpolation

    const double cosd = cos(p.delta);
    const double sind = -sin(p.delta); // minus because of ANTI-NU !!!

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

    // multiplying probability with scale factor to get events per MeV
    no[i] *= MB_anti_norm[i];

  } // end for i (bins)
  return;
}

/********************  the chisq *********************************/

double chisq_MB_anti(const int bin_0, params p)
{
  double no[N_BINS_MA], s = 0.;
  MB_A_events(p, no);

  for(int i = bin_0; i < N_BINS_MA; i++){
    double n = (no[i] + MB_A_data[i][0]) * (E_mba_max(i)-E_mba_min(i)) * 1000.;
    double d = MB_A_data[i][1] * (E_mba_max(i)-E_mba_min(i)) * 1000.;
    s += 2.*(n - d + d * log(d/n));
  }
  return s;

    

  for(int i = bin_0; i < N_BINS_MA; i++)
    no[i] += MB_A_data[i][0] - MB_A_data[i][1]; 
  
  for(int i = bin_0; i < N_BINS_MA; i++)
    for(int j = bin_0; j < N_BINS_MA; j++)
      s += no[i] * no[j] * (bin_0 ? s_inv_475_mba[i+1][j+1] : s_inv_200[i+1][j+1]);
	  
  return s;
}


/******************** calc integrals ************************/

// energy resolution from MB talk Fermilab 11/04/07 slide 24
#define B 0.08
#define C 0.024

inline double sigma_mba(const double e){
  return sqrt(B * B * e + C * C);
}

// global vars
int bin_gl_MB_A, sin_cos_MA;
double dm2_gl_MB_A;

double intFuncMB_A(double Enu)
{
  double i1 = erf((E_mba_max(bin_gl_MB_A) - Enu)/(M_SQRT2 * sigma_mba(Enu))); 
  i1 -= erf((E_mba_min(bin_gl_MB_A) - Enu)/(M_SQRT2 * sigma_mba(Enu))); 

  const double a = 2.53 * dm2_gl_MB_A / Enu;
  
  double sin_cos_mean = 0.;

  if(sin_cos_MA == SIN)
  {
    if(a <= 0.)
      nrerror("[intFuncMB_A]: a should not be 0 in case of SIN");
    else
      sin_cos_mean = (cos(a*LU) - cos(a*LO))/(a * (LO-LU));
  }
  else
    sin_cos_mean = 
      (a <= 0. ? 1. : (sin(a*LO) - sin(a*LU))/(a * (LO-LU)) );
  
  return nue_flux_anti(Enu) * XS_anti(Enu) * i1 * (sin_cos_mean + 1.);
}

#define N_SIGM_MBA 4.

/******************************************
 *   initialization
 ******************************************/ 

void initMBanti(void)
{
  // reading the data
  FILE *fp = fopen(PATH"/Data-files/MB-anti-data.dat","r");
  if(fp == NULL)
    nrerror("[initMinib-anti] cannot open PATH/Data-files/MB-anti-data.dat");
  
  for(int i = 0; i < N_BINS_MA; i++)
  {
    for(int j = 0; j < 3; j++)
      if(fscanf(fp, "%lf", &MB_A_data[i][j]) != 1)
        nrerror("[initMinib-anti] error in reading PATH/Data-files/MB-anti-data.dat");
    
    MB_A_data[i][2] -= MB_A_data[i][1];
  }
  fclose(fp);   

  
  for(int k = 0; k < DM2ANZ; k++){
    
    dm2_gl_MB_A = dm2k(k);
    
    for(bin_gl_MB_A = 0; bin_gl_MB_A < N_BINS_MA; bin_gl_MB_A++){
      
      double E_lo = E_mba_min(bin_gl_MB_A) - sigma_mba(E_mba_min(bin_gl_MB_A)) * N_SIGM_MBA;
      if(E_lo < E_MIN)
	E_lo = E_MIN;
      double E_up = E_mba_max(bin_gl_MB_A) + sigma_mba(E_mba_max(bin_gl_MB_A)) * N_SIGM_MBA;
      if(E_up > E_MAX)
	E_up = E_MAX;
      
      for(sin_cos_MA = 0; sin_cos_MA < 2; sin_cos_MA++){

        if(k == 0 && sin_cos_MA == SIN)
	  cos_sin_MB_A[k][SIN][bin_gl_MB_A] = 0.;
        else
          cos_sin_MB_A[k][sin_cos_MA][bin_gl_MB_A] = qromb1(intFuncMB_A, E_lo, E_up, 1.0e-5);
      
        // correcting for the " + 1" 
        if(k == 0)
	  cos_sin_MB_A[k][sin_cos_MA][bin_gl_MB_A] *= 0.5;
        else
	  cos_sin_MB_A[k][sin_cos_MA][bin_gl_MB_A] -= cos_sin_MB_A[0][COS][bin_gl_MB_A];
      }
    }
  }    

  // normalization
  for(bin_gl_MB_A = 0; bin_gl_MB_A < N_BINS_MA; bin_gl_MB_A++)
  {
    for(int k = 1; k < DM2ANZ; k++)
    {
      for(sin_cos_MA = 0; sin_cos_MA < 2; sin_cos_MA++)
      {		
	 cos_sin_MB_A[k][sin_cos_MA][bin_gl_MB_A] /= cos_sin_MB_A[0][COS][bin_gl_MB_A];
      }
    }
    cos_sin_MB_A[0][COS][bin_gl_MB_A] = 1.;
  }
  
  // norm factor to scale to MB prediction
  for(int i = 0; i < N_BINS_MA; i++)
    MB_anti_norm[i] = 1.;
  
  params p;
  p.Ue[I4] = p.Um[I4] = 0.;
  p.delta = 0.;
  p.dmq[I4] = 0.;
#ifdef REF_POINT1
  p.dmq[I5] = 1.;
  p.Ue[I5] = p.Um[I5] = .177827941; // sqrt( (sin 2theta)/2 ) for sin^22theta=0.004
#endif
#ifdef REF_POINT2   
  p.dmq[I5] = 4.42;
  p.Ue[I5] = p.Um[I5] = .197613887; // sqrt( (sin 2theta)/2 ) for sin^22theta=0.0061
#endif
#ifdef REF_POINT3   
  p.dmq[I5] = .1;
  p.Ue[I5] = p.Um[I5] = 0.4728708; // sqrt( (sin 2theta)/2 ) for sin^22theta=0.2
#endif
  
  double prob[N_BINS_MA];
  MB_A_events(p, prob);
  
  for(int i = 0; i < N_BINS_MA; i++)
  {
    if(prob[i] <= 0.)
      fprintf(stderr, "MB anti init: non-postive probability: %e\n", prob[i]);
    MB_anti_norm[i] = MB_pred_anti[i] / prob[i];
  }
  
  // init covariance matrix
  int i, j, *indx;
  double **s_200, **s_475, d, *col;

  s_200 = matrix(1, N_BINS_MA, 1, N_BINS_MA);
  s_475 = matrix(1, N_BINS_MA, 1, N_BINS_MA);
  s_inv_200 = matrix(1, N_BINS_MA, 1, N_BINS_MA);
  s_inv_475_mba = matrix(1, N_BINS_MA, 1, N_BINS_MA);
  col = vector(1, N_BINS_MA);
  indx = ivector(1, N_BINS_MA);
  
  const double err_bg  = sqr(0.);
  const double err_cor = sqr(0.);

  // covar matrices 
  for(i = 1; i <= N_BINS_MA; i++){
    for(j = 1; j <= N_BINS_MA; j++)
    {
      s_200[i][j] = err_bg  * MB_A_data[i-1][0]*MB_A_data[j-1][0] + 
                    err_cor * MB_A_data[i-1][1]*MB_A_data[j-1][1];
      
      s_475[i][j] = (i > 3 && j > 3 ? s_200[i][j] : 0.);

    }    
    s_200[i][i] += sqr(MB_A_data[i-1][2]);
    s_475[i][i] += sqr(MB_A_data[i-1][2]);
  }
  
  // inverting the covar matrix for MB_A475
  ludcmp(s_475, N_BINS_MA, indx, &d);
  
  for(j = 1; j <= N_BINS_MA; j++)
  {    
    for(i = 1; i <= N_BINS_MA; i++) col[i] = 0.0;
    
    col[j] = 1.0;
    lubksb(s_475, N_BINS_MA, indx, col);
    
    for(i = 1; i <= N_BINS_MA; i++) s_inv_475_mba[i][j] = col[i];
  }
  // inverting the covar matrix for MB_A200
  ludcmp(s_200, N_BINS_MA, indx, &d);
  
  for(j = 1; j <= N_BINS_MA; j++)
  {    
    for(i = 1; i <= N_BINS_MA; i++) col[i] = 0.0;
    
    col[j] = 1.0;
    lubksb(s_200, N_BINS_MA, indx, col);
    
    for(i = 1; i <= N_BINS_MA; i++) s_inv_200[i][j] = col[i];
  }
  
  return;
}


void plot_excess_mba(params p)
{   
  double no[N_BINS_MA];
  MB_A_events(p, no);

  for(int i = 0; i < N_BINS_MA; i++){
	
     printf("%e %e\n", E_mba_min(i), no[i]);
     printf("%e %e\n", E_mba_max(i), no[i]);
  }   
  exit(0);
}

void plot_excess_mba(const char *name, params p)
{   
  double no[N_BINS_MA];
  MB_A_events(p, no);
  
  FILE *fp = fopen(name, "w");

  for(int i = 0; i < N_BINS_MA; i++){
	
     fprintf(fp, "%e %e\n", E_mba_min(i), no[i]);
     fprintf(fp, "%e %e\n", E_mba_max(i), no[i]);
  }
  fclose(fp);
  return;
}


/********************************************
 *   testing miniboone
 ********************************************/

#ifdef MINIB_ANTI

void plot_data(void)
{   
  for(int i = 0; i < N_BINS_MA; i++)
     printf("%e %e %e %e\n", 0.5*(E_mba_min(i)+E_mba_max(i)), 
	    MB_A_data[i][1], MB_A_data[i][1] - MB_A_data[i][0], MB_A_data[i][2]);
  exit(0);
}

void plot_spect(params p)
{   
  double no[N_BINS_MA];
  MB_A_events(p, no);

  for(int i = 0; i < N_BINS_MA; i++){
	
     printf("%e %e\n", E_mba_min(i), no[i]);
     printf("%e %e\n", E_mba_max(i), no[i]);
  }   
  exit(0);
}


int main(void)
{
  initMBanti();
  //plot_data(); 
   
  params p;

      
  p.Ue[0] = 0.;
  p.Um[0] = 0.;
  p.delta = 0.;
  p.dmq[0] = 0.;
   
  //p.dmq[I5] = 4.42; p.Ue[I5] = p.Um[I5] = .197613887; // sqrt( (sin 2theta)/2 ) for sin^22theta=0.0061   
  //p.dmq[I5] = 1.;p.Ue[I5] = p.Um[I5] = .177827941; // sqrt( (sin 2theta)/2 ) for sin^22theta=0.004
  //plot_spect(p);   

  p.dmq[I5] = 1.;
  p.Ue[I5] = p.Um[I5] = 0.; 
  
  //double min = 1.e4;
  //params p_min;

  for(int i = 0; i < 101; i++){
    const double sq = pow(10., -4. + 4. * i/100.);
    p.Ue[1] = 0.01;
    p.Um[1] = sqrt(sq) / (p.Ue[1] * 2.);

    for(int j = 0; j < 101; j++){
      p.dmq[1] = pow(10., -2. + 4. * j/100.);

      printf("%e  %e  %e\n", sq, p.dmq[1], chi2_MBA_475(p));
      /*
      double w = chi2mb300(p);
      if(w < min){
	min = w;
	p_min = p;
      }
       */
    }
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
#undef N_SIGM_MBA
