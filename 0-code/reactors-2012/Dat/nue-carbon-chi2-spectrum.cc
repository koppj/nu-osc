#include <stdlib.h>
#include <stdio.h>
#include <math.h>


int main(void)
{

  // reading LSND/KARMEN table from Pedro
  double sq2t, dmq, chisq_L, chisq_K, chisq_LK;

  FILE *fp = fopen("nue-carbon-chi2-spectrum.dat","r");

  for(int i = 0; i < 14443; i++){
    if(fscanf(fp, "%lf %lf %lf %lf %lf", 
              &sq2t, &dmq, &chisq_L, &chisq_K, &chisq_LK) != 5){ 
      fprintf(stderr, "error in reading chisq data");
      exit(0);
    }

    printf("%le %le %le\n", sq2t, log10(dmq), chisq_LK);
  }
  return 0;
}
