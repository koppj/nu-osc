#ifndef __NSI_H
#define __NSI_H

#include <globes/globes.h>

/* Names of non-standard parameters */
#define FIRST_NSI_INDEX  6
enum epsilon { ABS_EPSILON_S_EE=FIRST_NSI_INDEX, ARG_EPSILON_S_EE,
               ABS_EPSILON_S_MUE,    ARG_EPSILON_S_MUE,
               ABS_EPSILON_S_TAUE,   ARG_EPSILON_S_TAUE, /* 10 11 */
               ABS_EPSILON_S_EMU,    ARG_EPSILON_S_EMU, 
               ABS_EPSILON_S_MUMU,   ARG_EPSILON_S_MUMU,
               ABS_EPSILON_S_TAUMU,  ARG_EPSILON_S_TAUMU,
               ABS_EPSILON_S_ETAU,   ARG_EPSILON_S_ETAU,
               ABS_EPSILON_S_MUTAU,  ARG_EPSILON_S_MUTAU, /* 20 21 */
               ABS_EPSILON_S_TAUTAU, ARG_EPSILON_S_TAUTAU,

               EPSILON_M_EE, 
               ABS_EPSILON_M_EMU,    ARG_EPSILON_M_EMU, /* 25 26 */
               ABS_EPSILON_M_ETAU,   ARG_EPSILON_M_ETAU,
               EPSILON_M_MUMU,
               ABS_EPSILON_M_MUTAU,  ARG_EPSILON_M_MUTAU,
               EPSILON_M_TAUTAU,

               ABS_EPSILON_D_EE,     ARG_EPSILON_D_EE,  /* 33 34 */
               ABS_EPSILON_D_MUE,    ARG_EPSILON_D_MUE,
               ABS_EPSILON_D_TAUE,   ARG_EPSILON_D_TAUE,
               ABS_EPSILON_D_EMU,    ARG_EPSILON_D_EMU, 
               ABS_EPSILON_D_MUMU,   ARG_EPSILON_D_MUMU, /* 41 42 */
               ABS_EPSILON_D_TAUMU,  ARG_EPSILON_D_TAUMU,
               ABS_EPSILON_D_ETAU,   ARG_EPSILON_D_ETAU,
               ABS_EPSILON_D_MUTAU,  ARG_EPSILON_D_MUTAU,
               ABS_EPSILON_D_TAUTAU, ARG_EPSILON_D_TAUTAU }; /* 49 50 */
#define LAST_NSI_INDEX ARG_EPSILON_D_TAUTAU

/* Effective electron numbers for calculation of matter profile */
#define Ne_MANTLE       0.497
#define Ne_CORE         0.468

#define REARTH           6371  /* km */
#define RCORE            3480  /* km */

/* Macros */
#define SQR(x)      ((x)*(x))                              // x^2 
#define POW10(x)    (exp(M_LN10*(x)))                      // 10^x
#define MIN(X,Y)    ( ((X) < (Y)) ? (X) : (Y) )
#define MAX(X,Y)    ( ((X) > (Y)) ? (X) : (Y) )
#define SIGN(a,b)   ( (b) > 0.0 ? (fabs(a)) : (-fabs(a)) )


/* Function declarations */
/* --------------------- */

/* probability.c */               
int nsi_init_probability_engine();
int nsi_free_probability_engine();
int nsi_set_oscillation_parameters(glb_params p, void *user_data);
int nsi_get_oscillation_parameters(glb_params p, void *user_data);
int nsi_probability_matrix(double P[3][3], int cp_sign, double E,
        int psteps, const double *length, const double *density,
        double filter_sigma, void *user_data);

/* prem.c */
int LoadPREMProfile(const char *prem_file);
double GetPREMDensity(double t, double L);
double GetAvgPREMDensity(double L_tot, double L1, double L2);
int GetPREM3LayerApprox(double L, int *n_layers, double *lengths,
                        double *densities);

#endif

