#include <stdio.h>
#include <math.h>
#include <stdlib.h>

inline double sqr(double x){
   return x*x;
}

inline double min(double x, double y){
   return x < y ? x : y;
}


inline void nrerror(char error_text[]){
   fprintf(stderr, "%s\n", error_text);
   exit(1);
   return;
}


void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
double trapzd(double (*func)(double), double a, double b, int n);
double *vector(long nl, long nh);
void free_vector(double *v, long nl, long nh);
int *ivector(long nl, long nh);
void free_ivector(int *v, long nl, long nh);
double **matrix(long nrl, long nrh, long ncl, long nch);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
double qromb(double (*func)(double), double a, double b);
double rtbis(double (*func)(double), double x1, double x2, double xacc);
double qromb1(double (*func)(double), double a, double b, double eps);
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);
