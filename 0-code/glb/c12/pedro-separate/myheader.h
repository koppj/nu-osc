#include <globes/globes.h>
#include <gsl/gsl_math.h>

/* Effective electron numbers for calculation of matter profile */
#define Ne_MANTLE       0.497
#define Ne_CORE         0.468

#define REARTH           6371  /* km */
#define RCORE            3480  /* km */

#define GLB_M0 6
#define GLB_R  7

/* Macros */
#define SQR(x)      ((x)*(x))                              // x^2 
#define POW10(x)    (exp(M_LN10*(x)))                      // 10^x
#define MIN(X,Y)    ( ((X) < (Y)) ? (X) : (Y) )
#define MAX(X,Y)    ( ((X) > (Y)) ? (X) : (Y) )
#define SIGN(a,b)   ( (b) > 0.0 ? (fabs(a)) : (-fabs(a)) )


/* Function declarations */
/* --------------------- */

/* probability.c */               
int twoflv_set_oscillation_parameters(glb_params p, void *user_data);
int twoflv_get_oscillation_parameters(glb_params p, void *user_data);
int twoflv_probability_matrix(double P[3][3], int cp_sign, double E,
        int psteps, const double *length, const double *density,
        double filter_sigma, void *user_data);



