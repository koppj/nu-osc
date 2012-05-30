#include <cmath>
#include <cstdlib>
#include <string.h>
 
#include "mynrfunc.h"

#ifndef PATH
//# define PATH "/afs/cern.ch/user/s/schwetz/work-home/4nuFiles/5nu/code-04-miniboone"
//# define PATH "/home/lin/schwetz/work-home/4nuFiles/5nu/code-10-react_2011"
# define PATH "./sbl/"
#endif

// Uncomment the following to use a 1+3+1 model (otherwise, 3+2 is used)
//#define Ip3pI  // FIXME

// definitions of parameter range
#define LDM2MIN -2.
#define LDM2MAX 2. 

#define DM2ANZ 101
#define DM_STEP ( (LDM2MAX - LDM2MIN) / (DM2ANZ - 2.) )  // including dmq = 0

inline double dm2k(int k)
{
  if(k == 0)
    return 0.;

  if(k >= DM2ANZ)
  {
    char s[100];
    sprintf(s, "[dm2k]: k=%d out of range\n", k);
    nrerror(s);
  }

  return pow(10, LDM2MIN + (k-1) * DM_STEP);
}

#define I4 0
#define I5 1

#define SIN 1 // true
#define COS 0 // false

struct params{
  double Ue[2];
  double Um[2];
  double Ue3;
  double delta;
  double dmq[2];
};

/******************************************* 
 * DEFINITION OF DELTA:
 * delta = arg[ U_m4 U_e4^* U_m5^* U_e5 ]
 *******************************************/

struct fix_params{
  bool Ue[2];
  bool Um[2];
  bool Ue3;
  bool delta;
  bool dmq[2];
};

enum experiments{
  mb300,
  mb475,
  mba200,
  mba475,
  lsnd,
  karmen,
  nomad,
  reactor,
  cdhs,
  atm,
  NUM_EXP
};

#define MN 939.566
#define MP 938.272
#define DELTA 1.294
#define ME 0.511

double read(char *name, double x0);

// the minimization routine
double min_chisq(params *p, const fix_params f, const bool incl[NUM_EXP]);

// flux and cross section for LSND and KARMEN
void initFlux(void);
double flux_LK(double enu);
double crossSect(double epos);



/*****************************************************
 *              routines for the experiments
 *****************************************************/

void init_fluxes(void);          // defined in class_reactor-flux.cc
double chi2reactor(params p);
void calc_reactors(void);

void calcKarmen(void);
double chi2karmen(params p);

void initLSND(void);
double chi2lsnd(params p);
void spectrum_lsnd(const char *name, params p);

void initCDHS(void);
double chi2cdhs(params p);

void initNomad(void);
double chi2nomad(params p);

void initMiniboone(void);
double chisq_minib(const int bin_0, params p);
inline double chi2mb300(params p){
  return chisq_minib(0, p);
}
inline double chi2mb475(params p){
  return chisq_minib(2, p);
}
void plot_excess_mb(char *name, params p);

void initMBanti(void);
double chisq_MB_anti(const int bin_0, params p);
inline double chi2_MBA_200(params p){
  return chisq_MB_anti(0, p);
}
inline double chi2_MBA_475(params p){
  return chisq_MB_anti(3, p);
}
void plot_excess_mba(const char *name, params p);


// use 2-dim de-dm table for the atm chisq
// (compile the file atm-de-dm.cc)
//#define MINOS_NC 
// ifdefined MINOS_NC the MINOS chisq is added in this routine
// in an approximate way, see file atm-de-dm.cc for details
double chi2atm(params p);

// use 1-dim dm table for the atm chisq
// (do not compile the file atm-de-dm.cc)
/*
inline double chi2atm(params p){
  const double dmu = sqr(p.Um[I4]) + sqr(p.Um[I5]);
  return read(PATH"/Data-files/atm-data-dmu.dat", dmu);
}
*/

/*************************************************
 *    the chisq using the include flags          *
 *************************************************/

#if !defined(PC475) && !defined(PC300)

inline double chisq_main(params p, const bool incl[NUM_EXP])
{
  return
    (incl[mb300]  ? chi2mb300(p)   : 0.) +   
    (incl[mb475]  ? chi2mb475(p)   : 0.) +   
    (incl[mba200] ? chi2_MBA_200(p): 0.) +   
    (incl[mba475] ? chi2_MBA_475(p): 0.) +   
    (incl[lsnd]   ? chi2lsnd(p)    : 0.) +   
    (incl[karmen] ? chi2karmen(p)  : 0.) + 
    (incl[nomad]  ? chi2nomad(p)   : 0.) + 
    (incl[reactor]? chi2reactor(p) : 0.) +  
    (incl[cdhs]   ? chi2cdhs(p)    : 0.) +	
    (incl[atm]    ? chi2atm(p)     : 0.); 	
}

#else // chisq function to calculate the PC

# define MIN_DIS    60.780336  // to be corrected !!!
# define MIN_APP475 16.907701
# define MIN_APP300 18.526059

inline double chisq_main(params p, const bool incl[NUM_EXP])
{
  const double chisq_dis = 
    chi2reactor(p) + chi2cdhs(p) + chi2atm(p) - MIN_DIS;

  const double chisq_app = 
    chi2lsnd(p) + chi2karmen(p) + chi2nomad(p) +
# ifdef PC300
    chi2mb300(p) - MIN_APP300;
# else
    chi2mb475(p) - MIN_APP475;
# endif
  return 0.5 * (chisq_dis + chisq_app) + 200. * sqr(chisq_dis - chisq_app);
}

#endif

/***************************************************
 *      reading and printing parameters 
 ***************************************************/


inline void print_param_pretty(FILE *fp, params p)
{
  fprintf(fp, "Ue3 = %e\n", p.Ue3);
  fprintf(fp, "Ue4 = %e\n", p.Ue[I4]);
  fprintf(fp, "Um4 = %e\n", p.Um[I4]);
  fprintf(fp, "D41 = %e\n\n", p.dmq[I4]);
  fprintf(fp, "Ue5 = %e\n", p.Ue[I5]);
  fprintf(fp, "Um5 = %e\n", p.Um[I5]);
  fprintf(fp, "D51 = %e\n\n", p.dmq[I5]);
  fprintf(fp, "del/Pi = %e\n", p.delta/M_PI);
  return;
}
inline void print_param(char *name, params p)
{
  FILE *fp = fopen(name, "w");
  if(fp == NULL)
    nrerror("[print_param] cannot open file");

  fprintf(fp, "%e\n", p.Ue[I4]);
  fprintf(fp, "%e\n", p.Ue[I5]);
  fprintf(fp, "%e\n", p.Um[I4]);
  fprintf(fp, "%e\n", p.Um[I5]);
  fprintf(fp, "%e\n", p.Ue3);
  fprintf(fp, "%e\n", p.delta);
  fprintf(fp, "%e\n", p.dmq[I4]);
  fprintf(fp, "%e\n", p.dmq[I5]);
  fclose(fp);
  return;
}
inline void read_param(char *name, params *p)
{
  FILE *fp = fopen(name, "r");
  if(fp == NULL)
    nrerror("[read_param] cannot open file");

  if(fscanf(fp, "%le\n", &p->Ue[I4]) != 1 || 
     fscanf(fp, "%le\n", &p->Ue[I5])  != 1 || 
     fscanf(fp, "%le\n", &p->Um[I4]) != 1 || 
     fscanf(fp, "%le\n", &p->Um[I5]) != 1 || 
     fscanf(fp, "%le\n", &p->Ue3) != 1 || 
     fscanf(fp, "%le\n", &p->delta) != 1 || 
     fscanf(fp, "%le\n", &p->dmq[I4]) != 1 || 
     fscanf(fp, "%le\n", &p->dmq[I5]) != 1)
    nrerror("[read_param] cannot read parameters");

  fclose(fp);
  return;
}
