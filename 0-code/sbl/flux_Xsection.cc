#include "definitions.h"

double *fluxKoeff;

inline double pow1(double x, double y)
{
  return y == 0.0 ? 1. : pow(x,y);
}

double flux_LK(double enu)
{
  if(enu>52.8 || enu<0.0) return 0.0;
  
  double w = 0.;
  for (int i=1;i<=5;i++) 
    w += fluxKoeff[i] * pow1(enu, i-1);
  
  return w;
}


double crossSect(double eKin)
{
   const double f=1.0, g=1.26, f2=3.7;
   double e, p, psq, cs0, gamma;
   
   e = eKin + ME;
   psq = sqr(e) - sqr(ME);
   
   if(eKin <= 0. || psq <= 0)
     return 0.;
   
   p=sqrt(psq);
   cs0=(sqr(f)+3*sqr(g))*e*p;
   gamma=2.*(f+f2)*g*((2.*e+DELTA)-sqr(ME)/e);
   gamma+=(sqr(f)+sqr(g))*(DELTA+sqr(ME)/e);
   gamma+=(sqr(f)+3.*sqr(g))*e;
   gamma+= -2./3.*(sqr(f)-sqr(g))*(e+DELTA);
   return cs0-gamma/MP*e*p;
}



/******************* Initialisierung ********************/

#define ANZ 28

void initFlux(void)
{
  const double enuPkt[]={
    0,5,6,8,10,12,14,16,18,20,22,24,26,28,31,34,36,38,
    40,42,44,46,48,49,50,51,52,52.8},
  
  fluxPkt[]={
    0,2.2,3.1,5.5,9.2,13,17,22,27,32.5,38.5,44,49.5,
      56.5,65,73.5,79.5,84,88,92.5,96,99,101.2,
      102.2,103,103.5,103.8,104};
  
  int i,j,k,*indx;
  double **a,d;

  a=matrix(1,5,1,5);
  fluxKoeff=vector(1,5);
  indx=ivector(1,5);
  for(i=1;i<=5;i++){
    for(j=1;j<=5;j++){
      a[i][j]=0;
      for(k=0;k<ANZ;k++) a[i][j]+=pow1(enuPkt[k],i+j-2);
    }
    fluxKoeff[i]=0;
    for(k=0;k<ANZ;k++) fluxKoeff[i]+=fluxPkt[k]*pow1(enuPkt[k],i-1);
  }
  ludcmp(a,5,indx,&d);
  lubksb(a,5,indx,fluxKoeff);
  free_matrix(a,1,5,1,5);
  free_ivector(indx,1,5);
  return;
}
#undef ANZ

