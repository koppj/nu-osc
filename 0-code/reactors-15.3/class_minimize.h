#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ACC 1.e-4
#define MAX_N_MIN 10000
#define MAX_N_LINE_MIN 10000


class Minimize
{
 public:
   int N_PAR_MIN;
   
   /* on initialization provide the function to be minimized
    * double f(double x[n])
    * and the number of parameters n */
   Minimize(double (* f)(double *x), int n){
      N_PAR_MIN = n;
      min_func = f;
   }
   
   /* minimize the function: give startpoint x[n]
    * the parameter x[i] is kept fixed if fixed[i] == true
    * the parameters at the minimum are returned in x[n] */
   double find_min(double *x, bool *fixed);
   /* all parameters free */
   double find_min(double *x);
	
   
 private:   
   double (* min_func)(double *x);
   double minimize(double *P, bool *fixed);
   double line_min(double *P, double *n, int *num);
   void error(const char *m){
      fprintf(stderr, "%s\n", m);
      exit(0);
   }
};

