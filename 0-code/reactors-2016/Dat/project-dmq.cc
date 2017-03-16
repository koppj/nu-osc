#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NM 86
#define NS 41

int main(void)
{
  double dmq[NM], chiq[NM], min = 1.e6;
  for(int i = 0; i < NM; i++) chiq[i] = 1.e6;

  FILE *fp = fopen("3p1-SBL+Gal+C12+Solar+Kaml.out", "r");

  for(int i = 0; i < NS; i++){
    for(int j = 0; j < NM; j++){

      double c;
      if(fscanf(fp, "%*f %lf %lf", &dmq[j], &c) != 2){
	fprintf(stderr, "error in reading\n");
	exit(0);
      }
      if(c < chiq[j]) chiq[j] = c;
      if(c < min) min = c;
    }
  }

  for(int i = 0; i < NM; i++) 
    printf("%e %e %e\n", dmq[i], chiq[i], chiq[i] - min);
  return 0;
}
