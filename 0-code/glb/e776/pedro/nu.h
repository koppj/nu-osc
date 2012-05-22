#ifndef __NU_H
#define __NU_H

#define NU_PSEUDO_MPI

#ifdef NU_MPI
  #include <mpi.h>
#endif

#include <globes/globes.h>
#include "nsi.h"

/* MPI-related global variables and definition of parallelized for loop */
#if defined NU_MPI || defined NU_PSEUDO_MPI
  extern int mpi_rank;
  extern int mpi_size;
  #define MPIFOR(ii,min,max)                                               \
      for (int ii=min + (mpi_rank*(max+1-min))/mpi_size;               \
               ii < min + ((mpi_rank+1)*(max+1-min))/mpi_size; ii++)
#else
  #define MPIFOR(ii,min,max) for (int ii=min; ii <= max; ii++)
#endif

/* Mass hierarchies */
#define HIERARCHY_NORMAL     1
#define HIERARCHY_INVERTED  -1

/* Definition of command line or environment variable argument */
typedef struct
{
  char *name;
  int id;
} env_param;

/* The different scenarios we consider */
enum scenarios { SCENARIO_WBB_WC_60, SCENARIO_WBB_WC_120,
                 SCENARIO_WBB_LAR_60, SCENARIO_WBB_LAR_120,
                 SCENARIO_BBEAM, SCENARIO_DCHOOZ };

/* Different plot types */
enum nu_actions { NU_ACTION_TH13_DELTA, NU_ACTION_TH13_DISC_REACH,
                  NU_ACTION_TH13_SENS, NU_ACTION_TH13_PRECISION,
                  NU_ACTION_MH_SENS, NU_ACTION_CPV_SENS, NU_ACTION_OCT_SENS,
                  NU_ACTION_L_TH13, NU_ACTION_L_MH, NU_ACTION_L_CPV,
                  NU_ACTION_EXP_TH13, NU_ACTION_EXP_MH, NU_ACTION_EXP_CPV,
                  NU_ACTION_NSI_SENS, NU_ACTION_NSI_DISC_REACH, NU_ACTION_L_NSI,
                  NU_ACTION_TH13_FITS, NU_ACTION_TH13_FITS_RANDOM };

/* Experiment and rule numbers */
extern int EXP_BEAM_NEAR;
extern int EXP_BEAM_FAR;
extern int EXP_REACTOR_NEAR;
extern int EXP_REACTOR_FAR;

#define RULE_T2K_NUE_QE       0
#define RULE_T2K_NUE_BAR_QE   1
#define RULE_T2K_NUMU_QE      2
#define RULE_T2K_NUMU_BAR_QE  3
#define RULE_T2K_NUE_CC       4
#define RULE_T2K_NUE_BAR_CC   5

#define RULE_NOVA_NUE         0
#define RULE_NOVA_NUMU        1
#define RULE_NOVA_NUE_BAR     2
#define RULE_NOVA_NUMU_BAR    3

#define RULE_WBB_WC_NUE       0
#define RULE_WBB_WC_NUMU      1

#define RULE_WBB_LAR_E_CC     0
#define RULE_WBB_LAR_MU_CC    1
#define RULE_WBB_LAR_EBAR_CC  2
#define RULE_WBB_LAR_MUBAR_CC 3
#define RULE_WBB_LAR_E_QE     4
#define RULE_WBB_LAR_MU_QE    5
#define RULE_WBB_LAR_EBAR_QE  6
#define RULE_WBB_LAR_MUBAR_QE 7

/* Options for degfinder */
#define DEG_NO_NH      0x01   /* Omit normal hierarchy fits                    */
#define DEG_NO_IH      0x02   /* Omit inverted hierarchy fits                  */
#define DEG_NO_SYS     0x04   /* Switch off systematics                        */
#define DEG_NO_CORR    0x10   /* Switch off correlations                       */
#define DEG_NO_DEG     0x20   /* Switch off all degeneracies                   */
#define DEG_NO_NSI_DEG 0x40   /* Switch off degeneracies in the NSI parameters */


/* Starting values for systematics minimization */
#define MAX_SYS             200
extern double sys_startval_beam_plus[MAX_SYS];
extern double sys_startval_beam_minus[MAX_SYS];
extern double sys_startval_reactor[MAX_SYS];


/* Function declarations */
/* --------------------- */

/* degfinder.c */
int degfinder(const glb_params base_values, const int n_prescan_params,
      const int *prescan_params, const double *prescan_min,
      const double *prescan_max, const int *prescan_steps,
      const glb_projection prescan_proj, const glb_projection fit_proj,
      int *n_deg, glb_params *deg_pos, double *deg_chi2, const long flags);

/* sys.c */
double chiT2K(int exp, int rule, int n_params, double *x, double *errors,
              void *user_data);
double chiNOvA(int exp, int rule, int n_params, double *x, double *errors,
              void *user_data);
double chiWBB_WC(int exp, int rule, int n_params, double *x, double *errors,
                 void *user_data);
double chiWBB_WCfast(int exp, int rule, int n_params, double *x, double *errors,
                     void *user_data);
double chiWBB_LAr(int exp, int rule, int n_params, double *x, double *errors,
                  void *user_data);
double chiDCNorm(int exp, int rule, int n_params, double *x, double *errors,
              void *user_data);

/* sensitivities.c */
double ChiNPWrapper(glb_params base_values, double th12, double th13, double th23,
                    double delta, double dm21, double dm31, int hierarchy,
                    glb_params fit_values);
int baseline_opt(int (*sens_function)(const char *, double, double, int,
                                      double, double, int),
                 const char *key_string, double logs22th13min, double logs22th13max,
                 int logs22th13steps, double deltacp_min, double deltacp_max, int deltacp_steps);
int exposure_opt(int (*sens_function)(const char *, double, double, int,
                                      double, double, int),
                 const char *key_string, double logs22th13min, double logs22th13max,
                 int logs22th13steps, double deltacp_min, double deltacp_max, int deltacp_steps);

int best_so_fit(const char *key_string);
int s22th13_sens(const char *key_string, double logs22th13min, double logs22th13max,
       int logs22th13steps, double deltacp_min, double deltacp_max, int deltacp_steps);
int s22th13_delta_proj(const char *key_string, double logS22th13min, double logS22th13max,
       int logS22th13steps, double deltacp_min, double deltacp_max, int deltacp_steps);
int s22th13_disc_reach(const char *key_string, double logS22th13min, double logS22th13max,
       int logS22th13steps, double deltacp_min, double deltacp_max, int deltacp_steps);
int s22th13_precision(const char *key_string, double logS22th13min, double logS22th13max,
       int logS22th13steps, double deltacp_min, double deltacp_max, int deltacp_steps);
int mh_sens(const char *key_string, double logS22th13min, double logS22th13max, int logS22th13steps,
       double deltacp_min, double deltacp_max, int deltacp_steps);
int cpv_sens(const char *key_string, double logS22th13min, double logS22th13max,
       int logS22th13steps, double deltacp_min, double deltacp_max, int deltacp_steps);
int oct_sens(const char *key_string, double logS22th13min, double logS22th13max,
       int logS22th13steps, double deltacp_min, double deltacp_max, int deltacp_steps,
       double th23min, double th23max, int th23steps);

int th13_fits(const char *key_string);
int th13_fits_random(const char *key_string);
int nsi_sens(const char *key_string, double logepsmin, double logepsmax, int logepssteps,
             double argepsmin, double argepsmax, int argepssteps);
int nsi_disc_reach(const char *key_string, double logepsmin, double logepsmax, int logepssteps,
             double argepsmin, double argepsmax, int argepssteps);


#endif

