#include "definitions.h"

#define L_MEAN 17.7
#define LU (L_MEAN - 2.)
#define LO (L_MEAN + 2.)

// numer of events expected for max mixing and large dmq 
// from Tab. IV of KARMEN paper
#define TOT_EVENTS 2913.

#define N_BINS_K 9

// prompt energy
#define E_MIN 16.0
#define E_MAX 52.0

#define D_E ( (E_MAX - E_MIN)/N_BINS_K )

inline double E_min(int bin){
  return E_MIN + bin * D_E;
}
inline double E_max(int bin){
  return E_MIN + (bin + 1.) * D_E;
}

// converts visible energy to nu energy
inline double Evis2Enu(double Evis){ 
  return Evis - ME + DELTA;
}
// converts nu energy to visible energy
inline double Enu2Evis(double Enu){
  return Enu + ME - DELTA;
}
// converts nu energy to kinetic energy
inline double Enu2Ekin(double Enu){
  return Enu2Evis(Enu) - 2.*ME;
}

// global variables
double integral[DM2ANZ][2][N_BINS_K];

// data read off from Fig. 11 of the KARMEN paper
const double data[N_BINS_K] = {3., 4., 1., 3., 3., 1., 0., 0., 0.};
    
const double bg[N_BINS_K] = {
  3.370370,
  3.681481,
  3.370370,
  2.333333,
  1.140741,
  0.777778,
  0.570370,
  0.362963,
  0.155556
};

/********************  the chisq *********************************/


void karmen_events(params p, double no[N_BINS_K])
{
  for(int i = 0; i < N_BINS_K; i++){

    // interpolating on the table
    double Icos[3], Isin[3];

    for(int j = 0; j < 2; j++){  // j = 41, 51

      const double dmq = p.dmq[j];

#ifdef CHECK
      if(dmq <= 0.){
        Icos[j] = integral[0][COS][i];
        Isin[j] = 0.;
      }else
#endif
      {

        const double logD = log10(dmq);

	const int k = (logD < LDM2MIN ? 0 : int( (logD - LDM2MIN)/DM_STEP ) + 1);

	Icos[j] = integral[k][COS][i] + 
          (integral[k+1][COS][i] - integral[k][COS][i])/
	  (dm2k(k+1) - dm2k(k)) * (dmq - dm2k(k));

	Isin[j] = integral[k][SIN][i] + 
          (integral[k+1][SIN][i] - integral[k][SIN][i])/
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
        Icos[2] = integral[0][COS][i];
        Isin[2] = 0.;
    }else{

	const int k = (logD < LDM2MIN ? 0 : int( (logD - LDM2MIN)/DM_STEP ) + 1);

	Icos[2] = integral[k][COS][i] + 
          (integral[k+1][COS][i] - integral[k][COS][i])/
	  (dm2k(k+1) - dm2k(k)) * (dmq - dm2k(k));

	Isin[2] = integral[k][SIN][i] + 
          (integral[k+1][SIN][i] - integral[k][SIN][i])/
	  (dm2k(k+1) - dm2k(k)) * (dmq - dm2k(k));

#ifndef Ip3pI
	if(p.dmq[I5] - p.dmq[I4] < 0.)
	  Isin[2] *= -1.;
#endif  
    } // end interpolation

    const double n0 = integral[0][COS][i];
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

  } // end for i (bins)
  return;
}


double chi2karmen(params p)
{
  double no[N_BINS_K];
  karmen_events(p, no);
  
  double s = 0.;

  for(int i = 0; i < N_BINS_K; i++)
    s += no[i] + bg[i] - data[i]; 
 
  // using only the bins where data[i] > 0
  for(int i = 0; i < 6; i++)
    s += data[i] * log(data[i] / (no[i] + bg[i]) );

  return 2. * s;
}
  

/******************** calc integrals ************************/

// energy resolution
inline double sigma(double e){
  return 0.0115 * sqrt(e);
}

// global vars
int bin_gl, sin_cos_K;
double dm2_gl;

double intFuncK(double Enu)
{
  double i1 = erf((E_max(bin_gl) - Enu2Evis(Enu))/(M_SQRT2 * sigma(Enu2Evis(Enu)))); 
  i1 -= erf((E_min(bin_gl) - Enu2Evis(Enu))/(M_SQRT2 * sigma(Enu2Evis(Enu)))); 
  
  const double a = 2.53 * dm2_gl / Enu;
  
  double sin_cos_mean;

  if(sin_cos_K == SIN)
  {
    if(a <= 0.)
      nrerror("[intFuncK]: a should not be 0 in case of SIN");
    else
      sin_cos_mean = (cos(a*LU) - cos(a*LO))/(a * (LO-LU));
  }
  else
    sin_cos_mean = 
      (a <= 0. ? 1. : (sin(a*LO) - sin(a*LU))/(a * (LO-LU)) );

  return flux_LK(Enu) * crossSect(Enu2Ekin(Enu)) * i1 * (sin_cos_mean + 1.);
}
  

#define N_SIG 4.

void calcKarmen(void)
{
  for(int k = 0; k < DM2ANZ; k++){
    
    dm2_gl = dm2k(k);
    
    for(bin_gl = 0; bin_gl < N_BINS_K; bin_gl++){
      
      const double eUnt = Evis2Enu(E_min(bin_gl) - N_SIG * sigma(E_min(bin_gl)));
      const double eOb  = Evis2Enu(E_max(bin_gl) + N_SIG * sigma(E_min(bin_gl)));
            
      
      for(sin_cos_K = 0; sin_cos_K < 2; sin_cos_K++){

        if(k == 0 && sin_cos_K == SIN)
	  integral[k][SIN][bin_gl] = 0.;
        else
          integral[k][sin_cos_K][bin_gl] = qromb1(intFuncK, eUnt, eOb, 1.0e-5);
      
        // correcting for the " + 1" 
        if(k == 0)
	  integral[k][sin_cos_K][bin_gl] *= 0.5;
        else
	  integral[k][sin_cos_K][bin_gl] -= integral[0][COS][bin_gl];
      }
    }
  }  
  
  // normalization
  double norm_karmen = 0.;
  for(int i = 0 ; i < N_BINS_K; i++)
    norm_karmen += 0.5 * integral[0][COS][i];
  
  for(int i = 0; i < N_BINS_K; i++)
    for(int k = 0; k < DM2ANZ; k++)
      for(int s = 0; s < 2; s++)
        integral[k][s][i] *= TOT_EVENTS / norm_karmen;
  
  return;
}

