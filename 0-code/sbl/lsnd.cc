#include "definitions.h"

#define LUNTEN 25.85
#define LOBEN 34.15

// mean L weighted by 1/L^2
#define L_MEAN ( 2. * LUNTEN * LOBEN / (LUNTEN + LOBEN) )

#define EUNTEN 20.0
#define EOBEN 60.0

// definitions for LoE bins
#define DATANZ 11
#define LoE_MIN 0.4
#define LoE_MAX 1.5
#define D_LoE ( (LoE_MAX - LoE_MIN)/DATANZ )

inline double LoE_min(int bin){
  return LoE_MIN + bin * D_LoE;
}
inline double LoE_max(int bin){
  return LoE_MIN + (bin + 1.) * D_LoE;
}

// converts positron kin en to nu energy
inline double eNuePos(double eKin){ 
  return eKin + ME + DELTA;
}
// converts nu energy to positron kin en
inline double ePoseNu(double eNu){
  return eNu - ME - DELTA;
}

double integral_lsnd[DM2ANZ][2][DATANZ];

// data read off from Fig. 24 of the LSND paper
const double data[DATANZ] = {
  3.538699,
  9.619963,
  6.762933,
  8.249668,
  7.813089,
  6.446993,
  3.470363,
  4.522533,
 -0.5 ,
  0.061916,
  1.299866
};
  
const double data_plus_sigm[DATANZ] = {
  5.773333,
  13.404586,
  10.236271,
  11.660280,
  11.224571,
  9.609074,
  6.135894,
  7.373440,
  1.300810,
  1.052685,
  2.661969
};
  
const double bg[DATANZ] = {
  0.558758,
  1.799822,
  2.605896,
  2.791126,
  2.666394,
  2.045618,
  1.673177,
  1.610749,
  0.990957,
  0.495338,
  0.061904
};

// global variables
double eUnt, eOb, dm2_glob, sqinv[DATANZ], norm_lsnd;


/********************  the chisq ******************************************/

void lsnd_events(params p, double no[DATANZ], double *prob)
{
  double rate = 0.;

  for(int i = 0; i < DATANZ; i++){

    // interpolating on the table
    double Icos[3], Isin[3];

    for(int j = 0; j < 2; j++){  // j = 41, 51

      const double dmq = p.dmq[j];

#ifdef CHECK
      if(dmq <= 0.){
        Icos[j] = integral_lsnd[0][COS][i];
        Isin[j] = 0.;
      }else
#endif
      {
        const double logD = log10(dmq);

	const int k = (logD < LDM2MIN ? 0 : int( (logD - LDM2MIN)/DM_STEP ) + 1);

	Icos[j] = integral_lsnd[k][COS][i] + 
          (integral_lsnd[k+1][COS][i] - integral_lsnd[k][COS][i])/
	  (dm2k(k+1) - dm2k(k)) * (dmq - dm2k(k));

	Isin[j] = integral_lsnd[k][SIN][i] + 
          (integral_lsnd[k+1][SIN][i] - integral_lsnd[k][SIN][i])/
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
        Icos[2] = integral_lsnd[0][COS][i];
        Isin[2] = 0.;
    }else{

        const int k = (logD < LDM2MIN ? 0 : int( (logD - LDM2MIN)/DM_STEP ) + 1);

	Icos[2] = integral_lsnd[k][COS][i] + 
          (integral_lsnd[k+1][COS][i] - integral_lsnd[k][COS][i])/
	  (dm2k(k+1) - dm2k(k)) * (dmq - dm2k(k));

	Isin[2] = integral_lsnd[k][SIN][i] + 
          (integral_lsnd[k+1][SIN][i] - integral_lsnd[k][SIN][i])/
	  (dm2k(k+1) - dm2k(k)) * (dmq - dm2k(k));

#ifndef Ip3pI
	if(p.dmq[I5] - p.dmq[I4] < 0.)
	  Isin[2] *= -1.;
#endif
    } // end interpolation

    const double n0 = integral_lsnd[0][COS][i];
    const double cosd = cos(p.delta);
    const double sind = -sin(p.delta); // minus because of ANTI-NU !!!

    no[i] = n0 * cosd 
      + cosd * Icos[2]  + sind * Isin[2]   // 54 -term
      - cosd * Icos[I5] - sind * Isin[I5]  // 51 -term
#ifndef Ip3pI
      - cosd * Icos[I4] + sind * Isin[I4]; // 41 -term
#else
      - cosd * Icos[I4] - sind * Isin[I4]; // 41 -term
#endif
    
    no[i] *= 2. * p.Ue[I4] * p.Um[I4] * p.Ue[I5] * p.Um[I5];

    no[i] += 2. * sqr(p.Ue[I4] * p.Um[I4]) * (n0 - Icos[I4]) +  
             2. * sqr(p.Ue[I5] * p.Um[I5]) * (n0 - Icos[I5]);
    
    rate += no[i];

  } // end for i (bins)

  *prob = rate/norm_lsnd;

  return;
}

#define PROB 0.264e-2
#define PROB_SIG_SQ ( (sqr(0.067e-2) + sqr(0.045e-2)) / 4. )

double chi2lsnd(params p)
{
  double no[DATANZ], prob;
  lsnd_events(p, no, &prob);

  // chisq rate only
  const double chisq_rate = sqr(prob - PROB) / PROB_SIG_SQ;
  printf("chi2-rate: %g\n", chisq_rate);//FIXME FIXME

  // chisq shape  
  double z = 0., n = 0.;
  for(int i = 0; i < DATANZ; i++){
      
      z += no[i] * sqinv[i] * (bg[i]-data[i]);
      n += no[i] * sqinv[i] * no[i];
  }

  if (n <= 0.)   // Temporary fix for th14=0.0 problem. JK, 2011-04-28
    return 1.e60;
  else
  {
    const double a = - z/n;
    
    double s = 0.;
    for(int i = 0; i < DATANZ; i++)      
      s += sqr(a * no[i]+bg[i]-data[i]) * sqinv[i]; 
    
    return s + chisq_rate;
  }

//  const double a = - z/n;
//  
//  double s = 0.;
//  for(int i = 0; i < DATANZ; i++)      
//    s += sqr(a * no[i]+bg[i]-data[i]) * sqinv[i]; 
//  
//  return s + chisq_rate;
}
  

/******************** calc integrals ************************/

// energy resolution of 7.7% at 52.8 MeV
// assume scaling with 1/sqrt(E)
 
inline double sigma_L(double Epos){
  const double rel = 0.077 * sqrt(52.8 / eNuePos(Epos));  
  return rel * eNuePos(Epos);
}

int sin_cos;

double intFuncLSND(double e)
{
  double i1 = erf( (eOb - e) / (M_SQRT2 * sigma_L(e)) );
  i1 -= erf( (eUnt - e) / (M_SQRT2 * sigma_L(e)) );
  
  const double a = 2.53 * dm2_glob / eNuePos(e);

  double sin_cos_mean;

  if(sin_cos == SIN)
  {
    if(a <= 0.)
      nrerror("[intFuncLSND]: a should not be 0 in case of SIN");
    else
      sin_cos_mean = (cos(a*LUNTEN) - cos(a*LOBEN))/(a * (LOBEN-LUNTEN));
  }
  else
    sin_cos_mean = 
      (a <= 0. ? 1. : (sin(a*LOBEN) - sin(a*LUNTEN))/(a * (LOBEN-LUNTEN)) );

  return flux_LK(eNuePos(e)) * crossSect(e) * i1 * (sin_cos_mean + 1.);
}

void intBerechnen(void)
{
  for(int k = 0; k < DM2ANZ; k++){
    
    dm2_glob = dm2k(k);
    
    for(int i = 0; i < DATANZ; i++){
      
      eUnt = ePoseNu(L_MEAN / LoE_max(i));
      eOb  = ePoseNu(L_MEAN / LoE_min(i));
            
      double eTr1 = eUnt - 3. * sigma_L(eUnt);
      double eTr2 = eOb  + 3. * sigma_L(eOb); 
      
      if(eTr1 < 0.0)  eTr1 = 0.1;
      if(eTr2 > 52.8) eTr2 = 52.8;
      
      for(sin_cos = 0; sin_cos < 2; sin_cos++){

        if(k == 0 && sin_cos == SIN)
	  integral_lsnd[k][SIN][i] = 0.;
        else
          integral_lsnd[k][sin_cos][i] = qromb1(intFuncLSND, eTr1, eTr2, 1.0e-5);
      
        // correcting for the " + 1" 
        if(k == 0)
	  integral_lsnd[k][sin_cos][i] *= 0.5;
        else
	  integral_lsnd[k][sin_cos][i] -= integral_lsnd[0][COS][i];
      }
    }
  }  
  
  // normalization
  norm_lsnd = 0.;
  for(int i = 0 ; i < DATANZ; i++)
    norm_lsnd += integral_lsnd[0][COS][i];

  return;
}

/******************* Initialisierung ********************/

void initLSND(void)
{
  intBerechnen();

  for(int i = 0; i < DATANZ; i++)
    sqinv[i] = 1./sqr(data_plus_sigm[i] - data[i]);

  return;
}


/********************  print the spectrum ********************************/

void spectrum_lsnd(const char *name, params p)
{
  double no[DATANZ], prob;
  lsnd_events(p, no, &prob);

  double z = 0., n = 0.;
  for(int i = 0; i < DATANZ; i++){
      
      z += no[i] * sqinv[i] * (bg[i]-data[i]);
      n += no[i] * sqinv[i] * no[i];
  }
  const double a = - z/n;
  
  FILE *fp = fopen(name, "w");
  
  for(int i = 0; i < DATANZ; i++){

    fprintf(fp, "%e %e\n", LoE_min(i), a* no[i] + bg[i]);
    fprintf(fp, "%e %e\n", LoE_max(i), a* no[i] + bg[i]);      
  }
  fclose(fp);
  fprintf(stderr, "lsnd prob: %e\n", prob);    
  return;
}


/********************************************
 *   testing LSND
 ********************************************/

#ifdef LSND

/********************  print the spectrum ********************************/

void spectrum_lsnd(params p)
{
  double no[DATANZ], prob;
  lsnd_events(p, no, &prob);

  double z = 0., n = 0.;
  for(int i = 0; i < DATANZ; i++){
      
      z += no[i] * sqinv[i] * (bg[i]-data[i]);
      n += no[i] * sqinv[i] * no[i];
  }
  const double a = - z/n;
  
  for(int i = 0; i < DATANZ; i++){

    printf("%e %e\n", LoE_min(i), a* no[i] + bg[i]);
    printf("%e %e\n", LoE_max(i), a* no[i] + bg[i]);      
  }
  return;
}


int main(void)
{
  initFlux();
  initLSND();
  
  params p;
  /*  
  p.Ue[0] = 0.;
  p.Um[0] = 0.;
  p.delta = 0.;
  p.dmq[0] = 0.;

  for(int i = 0; i < 101; i++){
    const double sq = pow(10., -3. + 3. * i/100.);
    p.Ue[1] = 0.01;
    p.Um[1] = sqrt(sq) / (p.Ue[1] * 2.);

    for(int j = 0; j < 101; j++){
      p.dmq[1] = pow(10., -2. + 4. * j/100.);

      printf("%e  %e  %e\n", sq, p.dmq[1], chi2lsnd(p));
    }
  }
  return 0;
  */
  read_param("out-5nu/app_mb475.point", &p);
  //p.Ue[I5] = 0.;

  print_param_pretty(stderr, p);

  spectrum_lsnd(p);

  double no[DATANZ], prob;
  lsnd_events(p, no, &prob);

  fprintf(stderr, "prob = %e, sig = %e, pull = %f\n", prob, PROB, (prob - PROB) / sqrt(PROB_SIG_SQ));
  return 0;
}
  
#endif
