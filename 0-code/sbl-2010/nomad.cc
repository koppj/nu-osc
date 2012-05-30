#include "definitions.h"


/************************************
 *  ATTENTION L in km and E in GeV
 ************************************/

#define L_MEAN 0.625

// 290 m long decay tunnel:
#define LU (L_MEAN - 0.145)
#define LO (L_MEAN + 0.145)

#define N_BINS_N 1

// numbers introduced "by hand" to fit the 
// original 2nu NOMAD contour (hep-ex/0306037)
#define ERR_EFF 8.e6
#define E_SCALE 2.5

// neutrino energy
#define E_MIN 0.1
#define E_MAX 170.

#define D_E ( (E_MAX - E_MIN)/N_BINS_N )

inline double E_min_nom(int bin){
  return E_MIN + bin * D_E;
}
inline double E_max_nom(int bin){
  return E_MIN + (bin + 1.) * D_E;
}

// global variables
double cos_sin_nomad[DM2ANZ][2][N_BINS_N];


/*********** flux and cross section *********************/

inline double nue_flux_nom(const double Enu){
  return pow(10., read(PATH"/Data-files/nomad-flux.dat", Enu));
}
inline double XS_nom(const double Enu){
  return read(PATH"/Data-files/XCC-nu_e.dat", log10(Enu)) * Enu;
}

/********************  the chisq *********************************/

void nomad_prob(params p, double no[N_BINS_N])
{
  for(int i = 0; i < N_BINS_N; i++){

    // interpolating on the table
    double Icos[3], Isin[3];

    for(int j = 0; j < 2; j++){  // j = 41, 51

      const double dmq = p.dmq[j];

#ifdef CHECK
      if(dmq <= 0.){
        Icos[j] = cos_sin_nomad[0][COS][i];
        Isin[j] = 0.;
      }else
#endif
      {

        const double logD = log10(dmq);

	const int k = (logD < LDM2MIN ? 0 : int( (logD - LDM2MIN)/DM_STEP ) + 1);

	Icos[j] = cos_sin_nomad[k][COS][i] + 
          (cos_sin_nomad[k+1][COS][i] - cos_sin_nomad[k][COS][i])/
	  (dm2k(k+1) - dm2k(k)) * (dmq - dm2k(k));

	Isin[j] = cos_sin_nomad[k][SIN][i] + 
          (cos_sin_nomad[k+1][SIN][i] - cos_sin_nomad[k][SIN][i])/
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
        Icos[2] = cos_sin_nomad[0][COS][i];
        Isin[2] = 0.;
    }else{

	const int k = (logD < LDM2MIN ? 0 : int( (logD - LDM2MIN)/DM_STEP ) + 1);

	Icos[2] = cos_sin_nomad[k][COS][i] + 
          (cos_sin_nomad[k+1][COS][i] - cos_sin_nomad[k][COS][i])/
	  (dm2k(k+1) - dm2k(k)) * (dmq - dm2k(k));

	Isin[2] = cos_sin_nomad[k][SIN][i] + 
          (cos_sin_nomad[k+1][SIN][i] - cos_sin_nomad[k][SIN][i])/
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

  } // end for i (bins)
  return;
}


double chi2nomad(params p)
{
  double no[N_BINS_N];
  nomad_prob(p, no);

  double s = 0.;
  for(int i = 0; i < N_BINS_N; i++)
    s += ERR_EFF * sqr( no[i] );
  
  return s;
}
  
/******************** calc integrals ************************/

// energy resolution
inline double sigma(double e){
  return 0.2 * sqrt(e);
}

// global vars
int bin_gl_nomad, sin_cos_N;
double dm2_gl_nomad;

double intFuncN(double Enu)
{
  double i1 = erf((E_max_nom(bin_gl_nomad) - Enu)/(M_SQRT2 * sigma(Enu))); 
  i1 -= erf((E_min_nom(bin_gl_nomad) - Enu)/(M_SQRT2 * sigma(Enu))); 
  
  const double a = 2.53 * dm2_gl_nomad / Enu;
  
  double sin_cos_mean;

  if(sin_cos_N == SIN)
  {
    if(a <= 0.)
      nrerror("[intFuncN]: a should not be 0 in case of SIN");
    else
      sin_cos_mean = (cos(a*LU) - cos(a*LO))/(a * (LO-LU));
  }
  else
    sin_cos_mean = 
      (a <= 0. ? 1. : (sin(a*LO) - sin(a*LU))/(a * (LO-LU)) );
  
  return nue_flux_nom(Enu*E_SCALE) * XS_nom(Enu*E_SCALE) * i1 * (sin_cos_mean + 1.);
}

#define N_SIGM 4.

void initNomad(void)
{
  for(int k = 0; k < DM2ANZ; k++){
    
    dm2_gl_nomad = dm2k(k);
    
    for(bin_gl_nomad = 0; bin_gl_nomad < N_BINS_N; bin_gl_nomad++){
      
      double E_lo = E_min_nom(bin_gl_nomad) - sigma(E_min_nom(bin_gl_nomad)) * N_SIGM;
      if(E_lo < E_MIN)
	E_lo = E_MIN;
      double E_up = E_max_nom(bin_gl_nomad) + sigma(E_max_nom(bin_gl_nomad)) * N_SIGM;
      if(E_up > E_MAX)
	E_up = E_MAX;
      
      for(sin_cos_N = 0; sin_cos_N < 2; sin_cos_N++){

        if(k == 0 && sin_cos_N == SIN)
	  cos_sin_nomad[k][SIN][bin_gl_nomad] = 0.;
        else
          cos_sin_nomad[k][sin_cos_N][bin_gl_nomad] = qromb1(intFuncN, E_lo, E_up, 1.0e-5);
      
        // correcting for the " + 1" 
        if(k == 0)
	  cos_sin_nomad[k][sin_cos_N][bin_gl_nomad] *= 0.5;
        else
	  cos_sin_nomad[k][sin_cos_N][bin_gl_nomad] -= cos_sin_nomad[0][COS][bin_gl_nomad];
      }
    }
  }    

  // normalization
  for(bin_gl_nomad = 0; bin_gl_nomad < N_BINS_N; bin_gl_nomad++){
    for(int k = 1; k < DM2ANZ; k++){
      for(sin_cos_N = 0; sin_cos_N < 2; sin_cos_N++){	
	
	 cos_sin_nomad[k][sin_cos_N][bin_gl_nomad] /= cos_sin_nomad[0][COS][bin_gl_nomad];
      }
    }
    cos_sin_nomad[0][COS][bin_gl_nomad] = 1.;
  }
  return;
}


#undef L_MEAN 
#undef LU
#undef LO
#undef N_BINS_N 
#undef ERR_EFF 
#undef E_SCALE 
#undef E_MIN 
#undef E_MAX 
#undef D_E 
