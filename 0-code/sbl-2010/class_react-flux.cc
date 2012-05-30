#include "def-reactors.h"

Flux flux[NISO][2];  // old and new reactor fluxes for each isotope

void init_fluxes(void)
{
   const char *filename_old[NISO] = 
     {PATH"/Data-files/U235.dat", PATH"/Data-files/U238.dat", PATH"/Data-files/Pu239.dat", PATH"/Data-files/Pu241.dat"};
   const char *filename_new[NISO] = 
     {PATH"/Data-files/U235_new.dat", PATH"/Data-files/U238_new.dat", PATH"/Data-files/Pu239_new.dat", PATH"/Data-files/Pu241_new.dat"};

   for(int i = 0; i < NISO; i++){
     flux[i][OLD].init(filename_old[i]);
     flux[i][NEW].init(filename_new[i]);
   }
   return;
}

void Flux::init(const char *filename)
{
   FILE *fp;
   fp = fopen(filename,"r");
   if(fscanf(fp,"%d", &Ncoef) != 1) error("cannot read isotope");
   if(Ncoef > MAX_COEF) error("Ncoef too big");

   for(int i = 0; i < Ncoef; i++)
     if(fscanf(fp, "%le", &coef[i]) != 1) error("cannot read isotope");

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

