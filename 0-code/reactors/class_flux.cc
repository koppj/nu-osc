#include "definitions.h"

namespace ns_reactor
{

extern Flux flux[NISO][2];  // old and new polynomaial fluxes for each isotope;
extern int old_new_main;

Rate_pull_coupl rate_pull_coupl;

double global_flux(const int iso, const double Enu, const int old_new)
{
  if(old_new == OLD)
    return flux[iso][old_new].f(Enu);

  // new fluxes

  if(Enu > 8.)
    return flux[iso][NEW].f(Enu);    // Patrick flux polynomial

  //fprintf(stderr, "here\n");

  switch(iso){

    case U235: return read(REACTOR_PATH"Dat/Patrick-U235-flux.dat", Enu);
    case P239: return read(REACTOR_PATH"Dat/Patrick-Pu239-flux.dat", Enu);
    case P241: return read(REACTOR_PATH"Dat/Patrick-Pu241-flux.dat", Enu);
    case U238: return flux[iso][NEW].f(Enu);  // Saclay flux (polynomial) for U238
  }
  return -1.;
}

void Flux::init(const char *filename)
{
   FILE *fp;
   fp = fopen(filename,"r");
   if(fscanf(fp,"%d", &Ncoef) != 1) error("cannot read isotope");
   if(Ncoef > MAX_COEF) error("Ncoef too big");

   for(int i = 0; i < Ncoef; i++)
     if(fscanf(fp, "%le", &coef[i]) != 1){
       fprintf(stderr, "%s %d ", filename, i); 
       error("cannot read isotope");
     }

   fclose(fp);
   return;
}

double Flux::f(const double Enu)
{
  double w = 0., pot = 1.;
  for(int i = 0; i < Ncoef; i++){
    w += coef[i] * pot;
    pot *= Enu;
  }
  return exp(w);
}


//Vogel and Beacom:
//http://arxiv.org/abs/hep-ph/9903554 
double crossSect(double eKin)
{
   const double f=1.0, g=1.26, f2=3.7;
   double e, p, psq, cs0, gamma;
   
   e = eKin + ME;
   psq = norm(e) - norm(ME);
   
   if(eKin <= 0. || psq <= 0.)
     return 0.;
   
   p=sqrt(psq);
   cs0=(norm(f)+3.*norm(g))*e*p;
   gamma=2.*(f+f2)*g*((2.*e+DELTA)-norm(ME)/e);
   gamma+=(norm(f)+norm(g))*(DELTA+norm(ME)/e);
   gamma+=(norm(f)+3.*norm(g))*e;
   gamma+= -2./3.*(norm(f)-norm(g))*(e+DELTA);
   return cs0-gamma/MP*e*p;
}


/*********************************************************************/

int iso_gl_init;
double gl_DmqL;

double int_func0(double e){  
  const double EposKin = e - ME - DELTA;
  return crossSect(EposKin) * global_flux(iso_gl_init, e, old_new_main);
}
double int_func1(double e){  
  const double EposKin = e - ME - DELTA;
  return e * crossSect(EposKin) * global_flux(iso_gl_init, e, old_new_main);
}
double int_func2(double e){  
  const double EposKin = e - ME - DELTA;
  return e * e * crossSect(EposKin) * global_flux(iso_gl_init, e, old_new_main);
}


double int_corr(double eNu){
  double w;
  switch(iso_gl_init){
  case U235:
    w = read(REACTOR_PATH"Dat/Patrick-U235-err_cor.dat", eNu);
    break;
  case P239:
    w = read(REACTOR_PATH"Dat/Patrick-Pu239-err_cor.dat", eNu);
    break;
  case P241:
    w = read(REACTOR_PATH"Dat/Patrick-Pu241-err_cor.dat", eNu);
    break;
  default:
    w = -1.;
  }

  const double EposKin=eNu-ME-DELTA;
  return w * crossSect(EposKin) * global_flux(iso_gl_init, eNu, old_new_main);
}


double int_sinq(double e)
{
  //const double EposKin = e - ME - DELTA;
  //return crossSect(EposKin) * norm(sin(1.27*gl_DmqL/e)) * global_flux(iso_gl_init, e, old_new_main);
  
  double Cos;
  if(gl_DmqL <= 0.) Cos = 1.;
  else
  {
    const double a1 = 2.54 * gl_DmqL / e * 0.95;
    const double a2 = 2.54 * gl_DmqL / e * 1.05;
    Cos = (sin(a2) - sin(a1))/(a2 - a1);
  }

  const double EposKin = e - ME - DELTA;
  return crossSect(EposKin) * 0.5*(1.-Cos) * global_flux(iso_gl_init, e, old_new_main);
}




void init_fluxes(void)
{
   const char *filename_old[NISO] = 
     {REACTOR_PATH"Dat/U235.dat", REACTOR_PATH"Dat/U238.dat", REACTOR_PATH"Dat/Pu239.dat", REACTOR_PATH"Dat/Pu241.dat"};
   const char *filename_new[NISO] = 
     {REACTOR_PATH"Dat/U235_patrick.dat", REACTOR_PATH"Dat/U238_new.dat", REACTOR_PATH"Dat/Pu239_patrick.dat",
      REACTOR_PATH"Dat/Pu241_patrick.dat"};

   for(int i = 0; i < NISO; i++){
     flux[i][OLD].init(filename_old[i]);
     flux[i][NEW].init(filename_new[i]);
   }

   // caluculated pull couplings for rate measurments

   for(iso_gl_init = 0; iso_gl_init < NISO; iso_gl_init++){

     rate_pull_coupl.unc1[iso_gl_init][0] = qromb1(int_func0, EnuMIN, EnuMAX, 1.e-5);

     if(iso_gl_init != U238){
       rate_pull_coupl.unc1[iso_gl_init][1] = qromb1(int_func1, EnuMIN, EnuMAX, 1.e-5);
       rate_pull_coupl.unc1[iso_gl_init][2] = qromb1(int_func2, EnuMIN, EnuMAX, 1.e-5);

       rate_pull_coupl.cor1[iso_gl_init] = qromb1(int_corr, EnuMIN, EnuMAX, 1.e-5);
     }
   }
   return;
}

/*********************************************************/
/*  probabilities for rate experiments                   */
/*********************************************************/

// used to interpolate on DMQ_STE and DMQ_ATM

void Rate_coef::init(void)
{
#ifndef LOG_ATM
  error("Rate_coef: works only for logarithmic Dmq_atm");
#endif
  log_DmqL_min = ATM_MIN + log10(9.);        // ILL has shortest BL of    9 m
  log_DmqL_max = STE_MAX + log10(2000.);     // DB  has longest  BL of 1912 m 

  DmqL_min = exp10(log_DmqL_min);
  DmqL_max = exp10(log_DmqL_max);

  Del = (log_DmqL_max - log_DmqL_min) / (N_RATE_COEF - 1.);

  // calculate the coefficients
  for(iso_gl_init = 0; iso_gl_init < NISO; iso_gl_init++){

    n0[iso_gl_init] = qromb1(int_func0, EnuMIN, EnuMAX, 1.e-6);

    for(int i = 0; i  < N_RATE_COEF; i++){
      gl_DmqL = exp10(log_DmqL_min + i * Del);
      Sinq[iso_gl_init][i] = qromb1(int_sinq, EnuMIN, EnuMAX, 1.e-6);
    }
  }
  return;
}

// interpolate on the table
double Rate_coef::sinq(const double DmqL, const double f[NISO])
{
  if(DmqL < DmqL_min) error("Rate_coef::sinq: DmqL too small");
  if(DmqL >= DmqL_max) error("Rate_coef::sinq: DmqL too large");

  const double logD = log10(DmqL);

  int i = int( (logD - log_DmqL_min) / Del);
  if(i < 0 || i >= N_RATE_COEF-1) error("Rate_coef:sinq: index out of range");

  const double x = (logD - i * Del - log_DmqL_min) / Del;
  if(x < 0. || x > 1.) error("Rate_coef:sinq: x out of range");

  double w = 0., v = 0.;
  for(int is = 0; is < NISO; is++){
    w += f[is] * (Sinq[is][i] + (Sinq[is][i+1] - Sinq[is][i]) * x);
    v += f[is] * n0[is];
  }
  return w / v;
} 


double Rate_coef::prob(Param_5nu &p, const double L, const double f[NISO])
{
  double P = 1.;
  for(int i = 0; i < N_NU-1; i++){
    for(int j = i+1; j < N_NU; j++){
      const double DmqL = fabs(p.Dmq(j,i)) * L;
      if(DmqL > DmqL_min)
        P -= 4. * norm(p.Ue[i] * p.Ue[j]) * sinq(DmqL, f);
    }
  }
  return P; 
} 


/*********************************************************/


inline int ind(const int i){
  return i == 0 ? U235 : i+1;
}

void calc_cs_per_fission(void)
{

  double scale = 6.39 / rate_pull_coupl.unc1[U235][0];

  for(int i = 0; i < NISO; i++){
    printf("%e\n", rate_pull_coupl.unc1[i][0] * scale);
  }
  printf("\n");

  double unc[3];
  unc[0] = unc[1] = unc[2] = 0.;

  for(int i = 0; i < 25; i++){
    const double e = 2. + 0.25 * i;
    const double EposKin=e-ME-DELTA;
    unc[0] +=  norm(0.25 * read(REACTOR_PATH"Dat/Patrick-U235-err_unc.dat", e)  * global_flux(U235, e, NEW) * crossSect(EposKin));
    unc[1] +=  norm(0.25 * read(REACTOR_PATH"Dat/Patrick-Pu239-err_unc.dat", e) * global_flux(P239, e, NEW) * crossSect(EposKin));
    unc[2] +=  norm(0.25 * read(REACTOR_PATH"Dat/Patrick-Pu241-err_unc.dat", e) * global_flux(P241, e, NEW) * crossSect(EposKin));
  }
  for(int i = 0; i < 3; i++){
    unc[i] = sqrt(unc[i]) / rate_pull_coupl.unc1[ind(i)][0];
    printf("%f\n", unc[i]);
  } 
  printf("\n");

  const char *name[3] = 
    {REACTOR_PATH"Dat/Patrick-U235-err_cor.dat",REACTOR_PATH"Dat/Patrick-Pu239-err_cor.dat",REACTOR_PATH"Dat/Patrick-Pu241-err_cor.dat"};

  double cor[3][3];
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      cor[i][j] = 0.;

      for(int a1 = 0; a1 < 25; a1++){
        const double e1 = 2. + 0.25 * a1;
        for(int a2 = 0; a2 < 25; a2++){
          const double e2 = 2. + 0.25 * a2;

	  cor[i][j] += 
            0.25 * read(name[i], e1) * global_flux(ind(i), e1, NEW) * crossSect(e1-ME-DELTA) *
            0.25 * read(name[j], e2) * global_flux(ind(j), e2, NEW) * crossSect(e2-ME-DELTA);
	}
      }
      cor[i][j] /= rate_pull_coupl.unc1[ind(i)][0] * rate_pull_coupl.unc1[ind(j)][0];
      printf("%e  ", cor[i][j]);
    }
    printf("\n");
  }
  printf("\n");

  
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++)
      printf("%e  ", cor[i][j]/sqrt(cor[i][i] * cor[j][j]));
    printf("\n");
  }
  printf("\n");
  
  
  for(int i = 0; i < 3; i++)
    printf("%f  %f\n", sqrt(cor[i][i]), sqrt(cor[i][i] + norm(unc[i])));

  exit(0);
}

}
