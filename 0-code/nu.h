#ifndef __NU_H
#define __NU_H

#define NU_PSEUDO_MPI

#ifdef NU_MPI
  #include <mpi.h>
#endif

#include <globes/globes.h>
#include "const.h"
#include "snu.h"

// Macros
#define SQR(x)      ((x)*(x))                        // x^2
#define SQR_ABS(x)  (SQR(creal(x)) + SQR(cimag(x)))  // |x|^2
#define POW10(x)    (exp(M_LN10*(x)))                // 10^x
#define MIN(X,Y)    ( ((X) < (Y)) ? (X) : (Y) )
#define MAX(X,Y)    ( ((X) > (Y)) ? (X) : (Y) )
#define SIGN(a,b)   ( (b) > 0.0 ? (fabs(a)) : (-fabs(a)) )
#define SGN(a)      ( ((a) >= 0.0) ? (1) : (-1) )
#define KRONECKER(i,j)  ( (i)==(j) ? 1 : 0 )
#define ROUND(x)    ( (int)((x) + 0.5) )


/* MPI-related global variables and definition of parallelized for loop */
#if defined NU_MPI || defined NU_PSEUDO_MPI
  extern int mpi_rank;
  extern int mpi_size;
  #define MPIFOR(ii,min,max)                                           \
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

/* Additional dynamic parameters for wide band beam */
typedef struct
{
  int flags;
  double eff_1st_max_nu;
  double eff_2nd_max_nu;
  double eff_1st_max_nubar;
  double eff_2nd_max_nubar;
  double E_1st_min;
} wbb_params_type;

/* Different plot types */
enum { NU_ACTION_SPECTRUM, NU_ACTION_PARAM_SCAN, NU_ACTION_EXPOSURE_SCAN,
       NU_ACTION_CHECK_BF };

/* External analysis routines */
enum { EXT_MB=0x001,  EXT_MBANTI=0x002, EXT_KARMEN=0x004, EXT_LSND=0x008,
       EXT_SBL=0x010, EXT_NOMAD=0x020,  EXT_CDHS=0x040,
       EXT_ATM_TABLE=0x080, EXT_ATM_COMP=0x100};

/* Experiment and rule numbers */
extern int EXP_BEAM_NEAR;
extern int EXP_BEAM_FAR;
extern int EXP_REACTOR_NEAR;
extern int EXP_REACTOR_FAR;


#define KAMLAND_N_REACT      16  // Number of reactors in KamLAND experiments

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
#define DEG_ONLY_STD_DEG 0x40   /* Switch off degeneracies in all but the std. osc. params */

#define DEG_LOGSCALE   0x01   /* Use logarithmic scale for a parameter         */
#define DEG_PLUS_MINUS 0x02   /* Consider positive and negative values for a   */
                              /* param (only in combination with DEG_LOGSCALE) */
#define DEG_S22        0x04   /* Scan over sin^2 2*param rather than param     */

/* Flags for WBB analysis */ 
#define WBB_NO_1ST_MAX          0x01 
#define WBB_NO_2ND_MAX          0x02 
#define WBB_1ST_MAX_TOTAL_RATES 0x04 
#define WBB_2ND_MAX_TOTAL_RATES 0x08


/* Starting values for systematics minimization */
#define MAX_SYS             200
extern double sys_startval_beam_plus[MAX_SYS];
extern double sys_startval_beam_minus[MAX_SYS];
extern double sys_startval_reactor[MAX_SYS];


/* Function declarations */
/* --------------------- */

#ifdef __cplusplus
extern "C" {
#endif

/* degfinder.c */
double ChiNPWrapper(glb_params base_values, double th12, double th13, double th23,
                    double delta, double dm21, double dm31, int hierarchy,
                    glb_params fit_values);
int degfinder(const glb_params base_values, const int n_prescan_params,
      const int *prescan_params, const double *prescan_min,
      const double *prescan_max, const int *prescan_steps,
      const glb_projection prescan_proj, const glb_projection fit_proj,
      int *n_deg, glb_params *deg_pos, double *deg_chi2, const unsigned long flags,
      const unsigned long *prescan_flags);

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
double chiMINOS(int exp, int rule, int n_params, double *x, double *errors,
              void *user_data);
double chiKamLAND(int exp, int rule, int n_params, double *x, double *errors,
              void *user_data);
double chiLSNDspectrum(int exp, int rule, int n_params, double *x, double *errors,
              void *user_data);

/* mb.c */
int chiMB_init();
int chiMB_clear();
double chiMBanti_nu2010(int exp, int rule, int n_params, double *x, double *errors,
              void *user_data);

/* sensitivities.c */
double sample(double min, double max, int steps, int i);
typedef int (*sens_func)(const char *, double, double, int, double, double, int);
int RestrictWBBEnergyWindow(wbb_params_type *wbb_params);
int print_rates();
int param_scan(const char *key_string, int n_p, char *params[], double p_min[], double p_max[],
       int p_steps[], unsigned long p_flags[], int n_min_params, char *min_params[],
       int prescan_n_p, char *prescan_params[], double prescan_p_min[],
       double prescan_p_max[], int prescan_p_steps[], unsigned long prescan_p_flags[]);

/* prem.c */
int LoadPREMProfile(const char *prem_file);
double GetPREMDensity(double t, double L);
double GetAvgPREMDensity(double L_tot, double L1, double L2);
int GetPREM3LayerApprox(double L, int *n_layers, double *lengths,
                        double *densities);

/* prior.cc */
int ext_init(int ext_flags);
double my_prior(const glb_params in, void* user_data);

/* iface.c */
int checkBF(int n_flavors);

/* nu.c */
int load_exps(const int n_exps, char **exps);

#ifdef __cplusplus
}
#endif

#endif

