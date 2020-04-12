/***************************************************************************
 * Non-standard interactions in reactor and superbeam experiments          *
 ***************************************************************************
 * Author: Joachim Kopp, Toshihiko Ota                                     *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <libgen.h>
#include <argp.h>
#include <gsl/gsl_complex_math.h>

#include <globes/globes.h>   // GLoBES library 
#include "glb_types.h"
BEGIN_C_DECLS
#include "glb_path.h"
END_C_DECLS
#include "glb_error.h"
#include "const.h"
#include "snu.h"
#include "nu.h"
#include "sbl/definitions.h"
#include "reactors/definitions.h"
#ifdef NU_USE_MONTECUBES
  #include <montecubes/montecubes.h>
#endif

#ifdef NU_MPI
  #include <mpi.h>
#endif

#if defined NU_MPI || defined NU_PSEUDO_MPI
  int mpi_rank = -1;
  int mpi_size = -1;
#endif

extern const env_param actions[];
extern const env_param scenarios[];
extern const env_param ext_scenarios[];

// Mode of operation 
int action = -1;
int use_nsi_constraints = 0; // Whether to use the constraint \eps^s = \eps^{d\dag} 
static int n_exps = 0;       // Number of experiments included in the current simulations
static char *exps[64];       // A list of experiments included in the current simulation
static char *true_param_def=NULL; // The definition of true osc. params (-t option)
static char *input_error_def=NULL; // The definition of external priors (-c option)

// Configuration of parameter scan
static int n_scan_params = 0;       // Number of parameters to scan over
static char *scan_params[32];       // Names of parameters to scan over
static double scan_p_min[32];       // Minimal parameter values
static double scan_p_max[32];       // Maximal parameter values
static int    scan_p_steps[32];     // Stepsizes
static unsigned long scan_p_flags[32]; // Extra flags (e.g. DEG_LOGSCALE)
static int n_prescan_params = 0;    // Number of parameters for degfinder prescanr
static char *prescan_params[32];    // Names of parameters for degfinder prescan
static double prescan_p_min[32];    // Minimal parameter values for degfinder pescan
static double prescan_p_max[32];    // Maximal parameter values for degfinder prescan
static int    prescan_p_steps[32];  // Stepsizes for degfinder precan
static unsigned long prescan_p_flags[32]; // Extra flags for prescan (e.g. DEG_LOGSCALE)
int n_min_params = 0;               // Names of parameters to marginalize over
static char *min_params[32];        // Names of parameters to marginalize over


// Experiment numbers 
int EXP_BEAM_NEAR    = -1;
int EXP_BEAM_FAR     = -1;
int EXP_REACTOR_NEAR = -1;
int EXP_REACTOR_FAR  = -1;

// True oscillation parameters and external priors 
double true_theta12;
double true_theta13;
double true_theta23;
double true_deltacp;
double true_sdm;
double true_ldm;
static int n_flavors = 5;

// Parameters for WBB analysis 
wbb_params_type wbb_params =
{
  0,         // flags
  1.0,       // eff_1st_max_nu
  1.0,       // eff_2nd_max_nu
  1.0,       // eff_1st_max_nubar
  1.0,       // eff_2nd_max_nubar
  0.0        // E_1st_min
};

// Parameters for Thomas' reactor analysis
namespace ns_reactor
{
  extern int old_new_main;
}

// Parameters for Michele's atmospheric neutrino analysis
int atm_decouple_e = 0; // Whether to decouple electron neutrinos by hand
                        // and instead include sterile neutrino matter effects

// Require mixing angles > 0?
int theta_positive = 0;

// Compute best fit point at end of parameter scan?
int compute_bf = 0;

/* Density correlations. Lists, for each experiment, the experiment with which
 * the corresponding matter density is correlated */
int density_corr[GLB_MAX_EXP];

// Defines which detectors should have their L varied in baseline_opt 
int L_opt[GLB_MAX_EXP];

// Parameter and projection vectors 
glb_params true_values;
glb_params test_values;
glb_params central_values;
glb_params input_errors;
glb_projection proj;
glb_projection prescan_proj;

// Misc 
int debug_level = 0;              // Verbosity level
long default_degfinder_flags = 0; // Default options for degeneracy finder 
char nu_flags[4096]   = "";       // Miscellaneous options                 
long ext_flags        = 0;        // Which external inputs should be used? 
static char output_file[FILENAME_MAX]=""; // Path and base name of output file in MCMC mode
char mb_tune[FILENAME_MAX] = "";  // MC tune for MiniBooNE backgrounds



// Definitions for argp
const char *argp_program_version = "nu 0.1";
const char *argp_program_bug_address = "<jkopp@fnal.gov>";
static char argp_doc[] = "nu neutrino oscillation simulation";
static char argp_option_doc[] = "[options]";

// The actual list of command line options we accept
#define OPT_NO_IH          1000
#define OPT_NO_NH          1001
#define OPT_ATM_DECOUPLE_E 1002
#define OPT_BF             1003
#define OPT_THETA_POSITIVE 1004
#define OPT_MB_TUNE        1005
static struct argp_option cmdline_options[] = {
  {"flavors",    'f',"NUMBER",0,"Number of flavors to use. Can be 3, 4, or 5)" },
  {"action",     'a',"ACTION",0,"What to do PARAM_SCAN, etc. (see const.c for more)" },
  {"experiments",'e',"EXPS",  0,"Which experiments to include in the fit" },
  {"parameter",  'p',"PARAMS",0,"Parameters to scan: <param_name>,<min>,<max>,<steps>,[<flags>]"},
  {"prescan",    's',"PARAMS",0,"Parameters for prescan in degfinder: "
                                "<param_name>,<min>,<max>,<steps>,[<flags>]"},
  {"minimize",   'm',"PARAMS",0,"Parameters to marginalize: <param_name>[,<param_name>[, ...]]"},
  {"true_params",'t',"PARAMS",0,"True oscillation parameters: \"NAME=VALUE, ...\""},
  {"cons",       'c',"PARAMS",0,"External 1\\sigma priors: \"NAME=VALUE, ...\""},
  {"no-ih",       OPT_NO_IH, NULL, 0,"Omit inverted hierarchy in fit"},
  {"no-nh",       OPT_NO_NH, NULL, 0,"Omit normal hierarchy in fit"},
  {"atm-decouple-e", OPT_ATM_DECOUPLE_E, NULL, 0,"Decouple electron neutrinos in ATM code"},
  {"best-fit",    OPT_BF,NULL,0,"Compute best fit point at the end of parameter scan"},
  {"theta-positive", OPT_THETA_POSITIVE, NULL, 0,"Require all mixing angles to be > 0"},
  {"mb-tune",     OPT_MB_TUNE,"TUNE",0,"MC tune to use for MB backgrounds"},
  {"verbose",    'v',NULL,    0,"Show debug output (use multiple times for more)"},
  {"outfile",    'o',"FILENAME",0,"Path/base name of output files in MCMC mode"},
  { 0 }
};



// -------------------------------------------------------------------------
//                     H E L P E R   F U N C T I O N S
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
int eval_param_def(char *p, char *param_name, double *param_value)
// -------------------------------------------------------------------------
// Parse parameter definition of the form "PARAM1=1.234"
// -------------------------------------------------------------------------
{
  char *sp;
  sp = strchr(p, '=');
  if (!sp)
  {
    fprintf(stderr, "eval_param_def: Invalid argument: %s.\n", p);
    return -1;
  }

  *sp++ = '\0';
  strcpy(param_name, p);
  *param_value = atof(sp);
  *(sp-1) = '=';

  return 0;
}


// -------------------------------------------------------------------------
int eval_osc_param_def(char *p, int *param_index, double *param_value)
// -------------------------------------------------------------------------
// Evaluate oscillation parameter definition of the form "NAME=VALUE"
// (similar to eval_param_def, but checks that the parameter is actually
// known to GLoBES)
// -------------------------------------------------------------------------
{
  char param_name[100];
  if (eval_param_def(p, param_name, param_value) < 0)
  {
    fprintf(stderr, "eval_osc_param_def: Invalid argument: %s.\n", p);
    return -1;
  }
  
  int j;
  if ((j=glbFindParamByName(param_name)) >= 0)
  {
    *param_index = j;
    return 0;
  }

  fprintf(stderr, "eval_osc_param_def: Unknown oscillation parameter: %s.\n", param_name);
  return -1;
}


// -------------------------------------------------------------------------
int eval_true_params(char arg[])
// -------------------------------------------------------------------------
// Parse true oscillation parameters (-t option). The expected format for p
// is "PARAM1=###,PARAM2=###,PARAM3=###" (allowed delimiters are ' ', ',',
// ';', and '\t')
// CAUTION: THIS FUNCTION MUST NOT BE CALLED BEFORE THE PARAMETER VECTORS
// HAVE BEEN ALLOCATED
// -------------------------------------------------------------------------
{
  const char delim[] = " ,;\t";
  char *this_param = strtok(arg, delim);
  while (this_param)
  {
    int param_index;
    double param_value;
    if (eval_osc_param_def(this_param, &param_index, &param_value) < 0)
      return -3;
    glbSetOscParams(true_values, param_value, param_index);
    this_param = strtok(NULL, delim);
  }
  return 0;
}


// -------------------------------------------------------------------------
int eval_input_errors(char arg[])
// -------------------------------------------------------------------------
// Parse external priors (-c option). The expected format for p
// is "PARAM1=###,PARAM2=###,PARAM3=###" (allowed delimiters are ' ', ',',
// ';', and '\t')
// CAUTION: THIS FUNCTION MUST NOT BE CALLED BEFORE THE PARAMETER VECTORS
// HAVE BEEN ALLOCATED
// -------------------------------------------------------------------------
{
  const char delim[] = " ,;\t";
  char *this_param = strtok(arg, delim);
  while (this_param)
  {
    int param_index;
    double param_value;
    if (eval_osc_param_def(this_param, &param_index, &param_value) < 0)
      return -3;
    glbSetOscParams(input_errors, param_value, param_index);
    this_param = strtok(NULL, delim);
  }
  return 0;
}


/*************************************************************************** 
 * Read AEDL variable definitions from environment variable AEDL_VARIABLES *
 * Format is AEDL_VARIABLES="PARAM1=###,PARAM2=###,PARAM3=###"             *
 * (allowed delimiters are ' ', ',', ';', '\t')                            *
 ***************************************************************************/
int eval_aedl_variables()
{
  if (!getenv("AEDL_VARIABLES"))
    return 0;
  if (debug_level > 1)
    printf("# Reading variable AEDL_VARIABLES = \"%s\"\n", getenv("AEDL_VARIABLES"));

  const char delim[] = " ,;\t";
  char env_string[300];
  strcpy(env_string, getenv("AEDL_VARIABLES"));
  char *this_param = strtok(env_string, delim);
  while (this_param)
  {
    char param_name[100];
    double param_value;
    if (eval_param_def(this_param, param_name, &param_value) < 0)
      return -6;
    glbDefineAEDLVariable(param_name, param_value);
    if (debug_level > 1)
      printf("# AEDL variable: %s = %10.7g\n", param_name, param_value);
    this_param = strtok(NULL, delim);
  }

  return 0;
}


/*************************************************************************** 
 * Print AEDL variables                                                    * 
 ***************************************************************************/
#include "glb_parser_type.h"
extern glb_symrec *pre_sym_table;
int print_aedl_variables()
{
  glb_symrec *ptr;
  printf("# AEDL variables:\n");
  for (ptr=pre_sym_table; ptr != NULL; ptr=(glb_symrec *)ptr->next)
    printf("#   %s = %g\n", ptr->name, ptr->value.var);

  return 0;
}


/*************************************************************************** 
 * Evaluate environment variables to select mode of operation etc.         * 
 ***************************************************************************/
int eval_env()
{
  // NU_FLAGS: Miscellaneous flags that affect the simulation in different ways 
  if (debug_level > 1)
    printf("# Reading variable NU_FLAGS = \"%s\"\n", getenv("NU_FLAGS"));
  if (getenv("NU_FLAGS"))
    strcpy(nu_flags, getenv("NU_FLAGS"));

  return 0;
}


// -------------------------------------------------------------------------
double parse_expr(const char *str)
// -------------------------------------------------------------------------
// Calls Python to interpret str as an expression and returns the resulting
// floating point value.
// -------------------------------------------------------------------------
{
  FILE *f = NULL;
  char cmd[strlen(str) + 100];
  double result;

  sprintf(cmd, "python -c 'from math import *; print %s'", str);
  f = popen(cmd, "r");
  fscanf(f, "%lf", &result);
  pclose(f);

  return result;
}


// -------------------------------------------------------------------------
error_t parse_opt(int key, char *arg, struct argp_state *state)
// -------------------------------------------------------------------------
// The main command line parser
// -------------------------------------------------------------------------
{
  const char delim[] = " ,;\t";
  env_param *p = NULL;
  switch (key)
  {
    // -------------------------------------------------
    case 'f':
      if (atoi(arg) < 3  ||  atoi(arg) > 5)
      {
        fprintf(stderr, "Invalid number of flavors: %s.\n", arg);
        return -1;
      }
      else
      {
        n_flavors = atoi(arg);
      }
      break;

    // -------------------------------------------------
    case 'a':
      for (p=(env_param *) actions; p->name != NULL; p++)
      {
        if (strcasecmp(arg, p->name) == 0)
        {
          action = p->id;
          break;
        }
      } 
      if (p->name == NULL)
      {
        fprintf(stderr, "Not a valid action: %s.\n", arg);
        return -1;
      }
      break;

    // -------------------------------------------------
    case 'e':
    {
      char *this_exp = strtok(arg, delim);
      while (this_exp)
      {
        exps[n_exps++] = strdup(this_exp);
        this_exp = strtok(NULL, delim);
      }
      break;
    }

    // -------------------------------------------------
    case 'p':
    {
      char *p = strtok(arg, delim);
      if (p)
      {
        scan_params[n_scan_params] = strdup(p);
        p = strtok(NULL, delim);
        if (p)
        {
          scan_p_min[n_scan_params] = parse_expr(p);
          p = strtok(NULL, delim);
          if (p)
          {
            scan_p_max[n_scan_params] = parse_expr(p);
            p = strtok(NULL, delim);
            if (p)
            {
              scan_p_steps[n_scan_params] = atoi(p);
              scan_p_flags[n_scan_params] = 0;
              while ((p=strtok(NULL, delim)) != NULL)
              {
                if (strstr(p, "LOGSCALE"))
                  scan_p_flags[n_scan_params] |= DEG_LOGSCALE;
              }
              n_scan_params++;
              break;
            }
          }
        }
      }
      fprintf(stderr, "Invalid parameter range: %s.\n", arg);
      return -2;
    }

    // -------------------------------------------------
    case 's':
    {
      char *p = strtok(arg, delim);
      if (p)
      {
        prescan_params[n_prescan_params] = strdup(p);
        p = strtok(NULL, delim);
        if (p)
        {
          prescan_p_min[n_prescan_params] = parse_expr(p);
          p = strtok(NULL, delim);
          if (p)
          {
            prescan_p_max[n_prescan_params] = parse_expr(p);
            p = strtok(NULL, delim);
            if (p)
            {
              prescan_p_steps[n_prescan_params] = atoi(p);
              prescan_p_flags[n_prescan_params] = 0;
              while ((p=strtok(NULL, delim)) != NULL)
              {
                if (strstr(p, "LOGSCALE"))
                  prescan_p_flags[n_prescan_params] |= DEG_LOGSCALE;
                if (strstr(p, "S22"))
                  prescan_p_flags[n_prescan_params] |= DEG_S22;
                if (strstr(p, "PM"))
                  prescan_p_flags[n_prescan_params] |= DEG_PLUS_MINUS;
              }
              n_prescan_params++;
              break;
            }
          }
        }
      }
      fprintf(stderr, "Invalid parameter range for prescan: %s.\n", arg);
      return -2;
    }

    // -------------------------------------------------
    case 'm':
    {
      char *p = strtok(arg, delim);
      while (p)
      {
        min_params[n_min_params++] = strdup(p);
        p = strtok(NULL, delim);
      }
      break;
    }

    // -------------------------------------------------
    case 't':
    {
      if (!true_param_def)
        true_param_def = strdup(arg);
      else
      {
        char s[strlen(true_param_def)+strlen(arg)+1];
        strcpy(s, true_param_def);
        strcat(s, ";");
        free(true_param_def);
        true_param_def = strdup(strcat(s, arg));
      }
      break;
    }

    // -------------------------------------------------
    case 'c':
    {
      if (!input_error_def)
        input_error_def = strdup(arg);
      else
      {
        char s[strlen(input_error_def)+strlen(arg)+1];
        strcpy(s, input_error_def);
        strcat(s, ";");
        free(input_error_def);
        input_error_def = strdup(strcat(s, arg));
      }
      break;
//      n_constrained_params = 0;
//      char *p = strtok(arg, delim);
//      while (p)
//      {
//        constrained_params[n_constrained_params++] = strdup(p);
//        p = strtok(NULL, delim);
//      }
//      break;
    }

    // -------------------------------------------------
    case 'v':
    {
      debug_level++;
      break;
    }

    case 'o':
    {
      if (strlen(arg) > FILENAME_MAX-1)
      {
        fprintf(stderr, "Name of output file too long.\n");
        return -3;
      }
      strncpy(output_file, arg, FILENAME_MAX);
      output_file[FILENAME_MAX-1] = '\0';
      break;
    }

    // -------------------------------------------------
    case OPT_NO_IH:
      default_degfinder_flags |= DEG_NO_IH;
      break;

    case OPT_NO_NH:
      default_degfinder_flags |= DEG_NO_NH;
      break;

    case OPT_ATM_DECOUPLE_E:
      atm_decouple_e = 1;
      break;

    case OPT_BF:
      compute_bf = 1;
      break;

    case OPT_THETA_POSITIVE:
      theta_positive = 1;
      break;

    case OPT_MB_TUNE:
    {
      if (strlen(arg) > FILENAME_MAX-1)
      {
        fprintf(stderr, "Name of MB tune too long.\n");
        return -3;
      }
      strncpy(mb_tune, arg, FILENAME_MAX);
      mb_tune[FILENAME_MAX-1] = '\0';
      break;
    }
 
    // -------------------------------------------------
    default:
      return ARGP_ERR_UNKNOWN;
  }

  return 0;
}


// -------------------------------------------------------------------------
int load_exps(const int n_exps, char **exps)
// -------------------------------------------------------------------------
// Load experiments
// -------------------------------------------------------------------------
{
  int karmen_c12_loaded = 0;
  int lsnd_c12_loaded   = 0;
  int c12_combi_loaded  = 0;
  int mb_loaded         = 0;

  // Default density correlations: Every experiment independent 
  for (int i=0; i < GLB_MAX_EXP; i++)
    density_corr[i] = i;

  // Default for baseline optimization - all exps at same baseline, and all included 
  for (int i=0; i < GLB_MAX_EXP; i++)
    L_opt[i] = 1;

  // Define AEDL variables specified via environment variables 
  if (eval_aedl_variables() < 0)
    return -1;

  // Load experiments 
  glbClearExperimentList();
  for (int i=0; i < n_exps; i++)
  {
    // Set some paths (for some experiments the path is set only when they are
    // actually loaded to avoid conflicts)
    setenv("GLB_PATH", "glb/nova:glb/t2k:glb/bbeam:glb/dchooz:glb/wbb_wc:"
                       "glb/wbb_lar:glb/nufact:glb/mb-2012:"
                       "glb/kamland:glb/lsnd:glb/e776:glb/c12:"
                       "glb/icarus-2012:glb/icarus-2014:glb/opera", 1);
    glb_setup_path();

    // --------------------------------------------------------------
    // Wide band beams
    // --------------------------------------------------------------

    // LBNE-like wide band beams (1 detector)
    if (strcasecmp(exps[i], "WBB_WC_60") == 0)
    {
      if (GLB_ISNAN(glbGetAEDLVariable("MASS")))
        glbDefineAEDLVariable("MASS", 200.0);
      if (GLB_ISNAN(glbGetAEDLVariable("BEAM_POWER")))
        glbDefineAEDLVariable("BEAM_POWER", 0.525);
      glbInitExperiment("wbb_wc_60_near.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbInitExperiment("wbb_wc_60_far.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbInitExperiment("wbb_wc_60_anti_near.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbInitExperiment("wbb_wc_60_anti_far.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbSetFilterStateInExperiment(1, GLB_ON);
      glbSetFilterStateInExperiment(3, GLB_ON);
      density_corr[3] = 1; // Matter density for far detector correlated in \nu/\bar{\nu} runs 
      L_opt[0] = L_opt[2] = 0; // Do not include near detectors in baseline optimization 
      EXP_BEAM_NEAR = 0;
      EXP_BEAM_FAR  = 1;
      struct glb_experiment *ND  = glb_experiment_list[0];
      struct glb_experiment *NDa = glb_experiment_list[2];
      for (int i=0; i < ND->numofchannels; i++)  // Set NOSC-flags for all ND channels -> faster
      {
        ND->listofchannels[2][i] = (ND->listofchannels[2][i] % 10) + 10;
        ND->listofchannels[3][i] = (ND->listofchannels[3][i] % 10) + 10;
      }
      for (int i=0; i < NDa->numofchannels; i++)
      {
        NDa->listofchannels[2][i] = (NDa->listofchannels[2][i] % 10) + 10;
        NDa->listofchannels[3][i] = (NDa->listofchannels[3][i] % 10) + 10;
      }
    }

    else if (strcasecmp(exps[i], "WBB_WC_120") == 0)
    {
      if (GLB_ISNAN(glbGetAEDLVariable("MASS")))
        glbDefineAEDLVariable("MASS", 200.0);
      if (GLB_ISNAN(glbGetAEDLVariable("BEAM_POWER")))
        glbDefineAEDLVariable("BEAM_POWER", 0.7);
      glbInitExperiment("wbb_wc_120_near.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbInitExperiment("wbb_wc_120_far.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbInitExperiment("wbb_wc_120_anti_near.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbInitExperiment("wbb_wc_120_anti_far.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbSetFilterStateInExperiment(1, GLB_ON);
      glbSetFilterStateInExperiment(3, GLB_ON);
      density_corr[3] = 1; // Matter density for far detector correlated in \nu/\bar{\nu} runs 
      L_opt[0] = L_opt[2] = 0;
      EXP_BEAM_NEAR = 0;
      EXP_BEAM_FAR  = 1;
      struct glb_experiment *ND  = glb_experiment_list[0];
      struct glb_experiment *NDa = glb_experiment_list[2];
      for (int i=0; i < ND->numofchannels; i++)
      {
        ND->listofchannels[2][i] = (ND->listofchannels[2][i] % 10) + 10;
        ND->listofchannels[3][i] = (ND->listofchannels[3][i] % 10) + 10;
      }
      for (int i=0; i < NDa->numofchannels; i++)
      {
        NDa->listofchannels[2][i] = (NDa->listofchannels[2][i] % 10) + 10;
        NDa->listofchannels[3][i] = (NDa->listofchannels[3][i] % 10) + 10;
      }
    }

    else if (strcasecmp(exps[i], "WBB_LAR_60") == 0)
    {
      if (GLB_ISNAN(glbGetAEDLVariable("MASS")))
        glbDefineAEDLVariable("MASS", 34.0);
      if (GLB_ISNAN(glbGetAEDLVariable("BEAM_POWER")))
        glbDefineAEDLVariable("BEAM_POWER", 0.525);
      glbInitExperiment("wbb_lar_60_near.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbInitExperiment("wbb_lar_60_far.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbSetFilterStateInExperiment(1, GLB_ON);
      L_opt[0] = 0;
      EXP_BEAM_NEAR = 0;
      EXP_BEAM_FAR  = 1;
      struct glb_experiment *ND  = glb_experiment_list[0];
      for (int i=0; i < ND->numofchannels; i++)
      {
        ND->listofchannels[2][i] = (ND->listofchannels[2][i] % 10) + 10;
        ND->listofchannels[3][i] = (ND->listofchannels[3][i] % 10) + 10;
      }
    }

    else if (strcasecmp(exps[i], "WBB_LAR_120") == 0)
    {
      if (GLB_ISNAN(glbGetAEDLVariable("MASS")))
        glbDefineAEDLVariable("MASS", 34.0);
      if (GLB_ISNAN(glbGetAEDLVariable("BEAM_POWER")))
        glbDefineAEDLVariable("BEAM_POWER", 0.7);
      glbInitExperiment("wbb_lar_120_near.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbInitExperiment("wbb_lar_120_far.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbSetFilterStateInExperiment(1, GLB_ON);
      L_opt[0] = 0;
      EXP_BEAM_NEAR = 0;
      EXP_BEAM_FAR  = 1;
      struct glb_experiment *ND  = glb_experiment_list[0];
      for (int i=0; i < ND->numofchannels; i++)
      {
        ND->listofchannels[2][i] = (ND->listofchannels[2][i] % 10) + 10;
        ND->listofchannels[3][i] = (ND->listofchannels[3][i] % 10) + 10;
      }
    }

    // LBNE-like wide band beams (2 detectors at different baselines)
    else if (strcasecmp(exps[i], "WBB_WC_120_2BL") == 0)
    {
      if (GLB_ISNAN(glbGetAEDLVariable("MASS")))
        glbDefineAEDLVariable("MASS", 100.0);
      if (GLB_ISNAN(glbGetAEDLVariable("BEAM_POWER")))
        glbDefineAEDLVariable("BEAM_POWER", 0.7);
      glbInitExperiment("wbb_wc_120.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbInitExperiment("wbb_wc_120_anti.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbInitExperiment("wbb_wc_120.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbInitExperiment("wbb_wc_120_anti.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbSetBaselineInExperiment(2, 3*glbGetBaselineInExperiment(0));
      glbSetBaselineInExperiment(3, 3*glbGetBaselineInExperiment(0));
      glbSetFilterStateInExperiment(1, GLB_ON);
      glbSetFilterStateInExperiment(3, GLB_ON);
      density_corr[3] = 1; // Matter density for far detector correlated in \nu/\bar{\nu} runs 
      L_opt[0] = L_opt[1] = 0; // Baseline optimization for 2nd detector 
      EXP_BEAM_NEAR = -1;
      EXP_BEAM_FAR  = 1;
    }


    // --------------------------------------------------------------
    // Reactor experiments
    // --------------------------------------------------------------

    // Double Chooz (near + far)
//    else if (strcasecmp(exps[i], "DCHOOZ") == 0)
//    {
//      glbInitExperiment("D-Chooz_near.glb",&glb_experiment_list[0],&glb_num_of_exps);
//      glbInitExperiment("D-Chooz_far.glb",&glb_experiment_list[0],&glb_num_of_exps);
//      L_opt[0] = 0;
//      EXP_REACTOR_NEAR = 0;
//      EXP_REACTOR_FAR  = 1;
//      prior_th23 = 0.10 * true_theta23;
//      prior_ldm  = 0.05 * true_ldm;
//    }

    // KamLAND
    else if (strcasecmp(exps[i], "KAMLAND_Pedro") == 0)
    {
      // Baseline and thermal power for the 16 most important KamLAND reactors
      const double distance[KAMLAND_N_REACT] =  { 160.0,  179.0,  191.0,  214.0,  139.0,
                                                   87.7,   145.0,  349.0,  345.0,  295.0,
                                                  401.0,  561.0,  755.0,  430.0,  783.0,
                                                  830.0 };
      const double power[KAMLAND_N_REACT] = { 24.317, 13.692, 10.200, 10.600, 4.5,
                                               1.6,    4.927,  14.2,   13.172, 3.293,
                                               3.8,    5.96,   10.146, 6.465,  3.3,
                                               5.32 };
      for(int i=0; i < KAMLAND_N_REACT; i++)
      {
        glbDefineAEDLVariable("setpower",  power[i]);
        glbDefineAEDLVariable("setlength", distance[i]);
        glbInitExperiment("kamland.glb",  &glb_experiment_list[0], &glb_num_of_exps);
        if (i != 0)
          glbSetChiFunction(glb_num_of_exps-1, GLB_ALL, GLB_ON, "chiZero", NULL);     
      }
    }


    // --------------------------------------------------------------
    // Neutrino factory
    // --------------------------------------------------------------

    // IDS-NF neutrino factory
    else if (strcasecmp(exps[i], "NUFACT") == 0)
    {
      glbDefineAEDLVariable("EMAX",       25.0);
      glbDefineAEDLVariable("BASELINE", 4000.0);
      glbDefineAEDLVariable("MASS",      100.0);
      glbInitExperiment("NF-Gold-jk.glb", &glb_experiment_list[0], &glb_num_of_exps);
      L_opt[0] = 1;

      glbClearAEDLVariables();
      glbDefineAEDLVariable("EMAX",       25.0);
      glbDefineAEDLVariable("BASELINE", 7000.0);
      glbDefineAEDLVariable("MASS",       50.0);
      glbInitExperiment("NF-Gold-jk.glb", &glb_experiment_list[0], &glb_num_of_exps);
      L_opt[1] = 0;
    }


    // --------------------------------------------------------------
    // MINOS
    // --------------------------------------------------------------

    // MINOS NC analysis (http://arxiv.org/abs/1607.01176)
//    else if (strcasecmp(exps[i], "MINOS_NC_2016") == 0)
//    {
//      MINOS_2016_init();
//      setenv("GLB_PATH", "glb/minos-2016", 1);
//      glb_setup_path();
//      glbInitExperiment("minos-nc.glb", &glb_experiment_list[0], &glb_num_of_exps);
//      glbInitExperiment("minos-cc.glb", &glb_experiment_list[0], &glb_num_of_exps);
//      L_opt[1] = L_opt[3] = 0;
//      
//      for (int i=glb_num_of_exps-3; i < glb_num_of_exps; i+=2)  // Loop over ND for CC, NC
//      {
//        glb_experiment *e = glb_experiment_list[i];
//        e->probability_user_data  = new int(MINOS_ND_PROBABILITY);
//      }
//    }
//
//    // MINOS CC \nu_\mu analysis (based on http://arxiv.org/abs/1607.01176)
//    else if (strcasecmp(exps[i], "MINOS_CC_2016") == 0)
//    {
//      setenv("GLB_PATH", "glb/minos-2016", 1);
//      glb_setup_path();
//      glbInitExperiment("minos-cc.glb", &glb_experiment_list[0], &glb_num_of_exps);
//      L_opt[1] = 0;
//    }


    // MINOS Neutral Current analysis (http://arxiv.org/abs/1001.0336, Nu2010, and 1103.0340)
    else if (strcasecmp(exps[i], "MINOS_NC_2011_noDecayPipe") == 0)
    {
      setenv("GLB_PATH", "glb/minos-2011", 1);
      glb_setup_path();
      glbInitExperiment("minos-nc.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbInitExperiment("minos-cc.glb", &glb_experiment_list[0], &glb_num_of_exps);
      L_opt[1] = L_opt[3] = 0;
    }

    // MINOS Neutral Current analysis (http://arxiv.org/abs/1001.0336, Nu2010, and 1103.0340)
    else if (strcasecmp(exps[i], "MINOS_NC_2011") == 0)
    {
      setenv("GLB_PATH", "glb/minos-2011", 1);
      glb_setup_path();
      glbInitExperiment("minos-nc.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbInitExperiment("minos-cc.glb", &glb_experiment_list[0], &glb_num_of_exps);
      L_opt[1] = L_opt[3] = 0;
      
      for (int i=glb_num_of_exps-3; i < glb_num_of_exps; i+=2)  // Loop over ND for CC and NC analysis
      {
        glb_experiment *e = glb_experiment_list[i];
        e->probability_user_data  = new int(MINOS_ND_PROBABILITY);
//        for (int j=0; j < e->num_of_sm; j++)
//          minos_smear(e->smear_data[j], e->smear[j], e->lowrange[j], e->uprange[j]);
//        glbPrintExp(i);
//        getchar();
      }
    }

    // MINOS Neutral Current analysis (http://arxiv.org/abs/1001.0336, Nu2010, and 1103.0340)
    // without oscillations in near detector
    else if (strcasecmp(exps[i], "MINOS_NC_2011_noNDosc") == 0)
    {
      setenv("GLB_PATH", "glb/minos-2011", 1);
      glb_setup_path();
      glbInitExperiment("minos-nc.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbInitExperiment("minos-cc.glb", &glb_experiment_list[0], &glb_num_of_exps);
      L_opt[1] = L_opt[3] = 0;

      for (int i=glb_num_of_exps-3; i < glb_num_of_exps; i+=2)  // Loop over ND for CC and NC analysis
      {
        glb_experiment *e = glb_experiment_list[i];
        for (int j=0; j < e->numofchannels; j++) // Loop over channels
        {
          if (e->listofchannels[2][j] < 10)
            e->listofchannels[2][j] += 10;
          if (e->listofchannels[3][j] < 10)
            e->listofchannels[3][j] += 10;
        }
      }
    }

    // MINOS CC \nu_\mu analysis (http://arxiv.org/abs/1103.0340)
    else if (strcasecmp(exps[i], "MINOS_CC_2011") == 0)
    {
      setenv("GLB_PATH", "glb/minos-2011", 1);
      glb_setup_path();
      glbInitExperiment("minos-cc.glb", &glb_experiment_list[0], &glb_num_of_exps);
      L_opt[1] = 0;
    }

    // 2010 MINOS NC+CC analysis (http://arxiv.org/abs/1001.0336 and Neutrino 2010)
    else if (strcasecmp(exps[i], "MINOS_2010") == 0)
    {
      setenv("GLB_PATH", "glb/minos-2011/1001.0336", 1);
      glb_setup_path();
      glbInitExperiment("MINOS-NC-far.glb",  &glb_experiment_list[0], &glb_num_of_exps);
      glbInitExperiment("MINOS-NC-near.glb", &glb_experiment_list[0], &glb_num_of_exps);
      L_opt[1] = 0;
    }

    // MINOS ND nu_e test
    else if (strcasecmp(exps[i], "MINERVA_TEST") == 0)
    {
      setenv("GLB_PATH", "glb/minerva-test", 1);
      glb_setup_path();
      glbInitExperiment("minerva-test.glb",  &glb_experiment_list[0], &glb_num_of_exps);
    }


    // --------------------------------------------------------------
    // Other beam experiments
    // --------------------------------------------------------------

    // T2K
    else if (strcasecmp(exps[i], "T2K") == 0)
    {
//      glbInitExperiment("t2k-nd280.glb", &glb_experiment_list[0], &glb_num_of_exps);
//      glbInitExperiment("t2k-nd-generic.glb", &glb_experiment_list[0], &glb_num_of_exps);
      glbInitExperiment("t2k-fd.glb", &glb_experiment_list[0], &glb_num_of_exps);
    }

    // E776 \nu_e appearance search
    else if (strcasecmp(exps[i], "E776") == 0)
      glbInitExperiment("E776.glb", &glb_experiment_list[0], &glb_num_of_exps);

    // KARMEN \nu_e--C-12 scattering data
    else if (strcasecmp(exps[i], "KARMEN-C12") == 0)
    {
      glbInitExperiment("KARMEN-nuecarbon.glb", &glb_experiment_list[0], &glb_num_of_exps);
      karmen_c12_loaded = 1;
    }
    else if (strcasecmp(exps[i], "KARMEN-C12-JR") == 0)
    {
      glbInitExperiment("KARMEN-nuecarbon-JR.glb", &glb_experiment_list[0], &glb_num_of_exps);
      karmen_c12_loaded = 1;
    }

    // LSND \nu_e--C-12 scattering data
    else if (strcasecmp(exps[i], "LSND-C12") == 0)
    {
      glbInitExperiment("LSND-nuecarbon.glb", &glb_experiment_list[0], &glb_num_of_exps);
      lsnd_c12_loaded = 1;
    }

    // Joint LSND + KARMEN C-12 scattering analysis
    else if (strcasecmp(exps[i], "C12") == 0)
    {
      init_nue_carbon(1);
      c12_combi_loaded = 1;
    }

    // ICARUS (2012 analysis)
    else if (strcasecmp(exps[i], "ICARUS-2012") == 0)
    {
      init_icarus_2012();
    }

    // ICARUS (2014 analysis)
    else if (strcasecmp(exps[i], "ICARUS-2014") == 0)
    {
      init_icarus_2014();
    }

    // Pedro's OPERA 2013 analysis
    else if (strcasecmp(exps[i], "OPERA") == 0)
    {
      init_OPERA();
    }

    // Pedro's MiniBooNE simulation
    #ifdef NU_USE_NUSQUIDS
    else if (strcasecmp(exps[i], "MBall200") == 0)
    {
      chiMB_init(0);
      mb_loaded++;
    }
    else if (strcasecmp(exps[i], "MBall475") == 0)
    {
      chiMB_init(1);
      mb_loaded++;
    }
    #endif

//    else if (strcasecmp(exps[i], "MBall") == 0  ||  strcasecmp(exps[i], "MBneutrino200.glb") == 0)
//    {
//      chiMB_init(1, 1); // 0=off, 1=full E range, 2=only E > 475 MeV, 1st number for nu, 2nd for nubar
//      mb_loaded++;
//    }
//    else if (strcasecmp(exps[i], "MBall475") == 0)
//    {
//      chiMB_init(2, 2);
//      mb_loaded++;
//    }
//    else if (strcasecmp(exps[i], "MB") == 0  ||  strcasecmp(exps[i], "MBnu") == 0)
//    {
//      chiMB_init(2, 0);
//      mb_loaded++;
//    }
//    else if (strcasecmp(exps[i], "MBanti") == 0  ||  strcasecmp(exps[i], "MBnubar") == 0)
//    {
//      chiMB_init(0, 2);
//      mb_loaded++;
//    }
//    else if (strcasecmp(exps[i], "MB300") == 0  ||  strcasecmp(exps[i], "MBnu300") == 0)
//    {
//      chiMB_init(1, 0);
//      mb_loaded++;
//    }
//    else if (strcasecmp(exps[i], "MBanti200") == 0  ||  strcasecmp(exps[i], "MBnubar200") == 0)
//    {
//      chiMB_init(0, 1);
//      mb_loaded++;
//    }

    // Dummy scenario that doesn't do anything
    else if (strcasecmp(exps[i], "DUMMY") == 0)
    {
    }

    // Generic glb file
    else if (strstr(exps[i], ".glb") != NULL)
    {
      glbInitExperiment(exps[i], &glb_experiment_list[0], &glb_num_of_exps);
    }

    // No glb file, no known experiment - see if we have an external code that can handle it
    else
    {
      env_param *p = NULL;
      for (p=(env_param *) ext_scenarios; p->name != NULL; p++)
      {
        if (strcasecmp(exps[i], p->name) == 0)
        {
          ext_flags |= p->id;
          break;
        }
      } 
      // Unknown experiment -> error
      if (!p->name)
      {
        fprintf(stderr, "Not a valid experiment: %s.\n", exps[i]);
        return -2;
      }

      // A few sanity checks
      if (ext_flags & EXT_MB  && ext_flags & EXT_MB_300)
      {
        fprintf(stderr, "Do not use MB and MB_300 analyses simultaneously.\n");
        return -3;
      }
      if (ext_flags & EXT_MBANTI  && ext_flags & EXT_MBANTI_200)
      {
        fprintf(stderr, "Do not use MBANTI and MBANTI200 analyses simultaneously.\n");
        return -4;
      }
    }
  } // for(i)

  // More sanity checks
  if (c12_combi_loaded  &&  (karmen_c12_loaded || lsnd_c12_loaded))
  {
    fprintf(stderr, "Do not use combined C-12 analysis together with individual analyses.\n");
    return -5;
  }
  if (mb_loaded > 1)
  {
    fprintf(stderr, "Do not load MiniBooNE multiple times - use MBall for combined nu+nubar fit.\n");
    return -6;
  }

  for (int i=0; i < glb_num_of_exps; i++)
    glbOptimizeSmearingMatrixInExperiment(i);

  // Restrict energy window for wide band beam if requested 
  if (wbb_params.flags & WBB_NO_1ST_MAX  &&  wbb_params.flags & WBB_NO_2ND_MAX)
  {
    fprintf(stderr, "Cannot exclude 1st and 2nd maximum - no data left.\n");
    return -7;
  }
  else if (wbb_params.flags & WBB_NO_1ST_MAX)
    printf("# Including only 2nd oscillation maximum.\n");
  else if (wbb_params.flags & WBB_NO_2ND_MAX)
    printf("# Including only 1st oscillation maximum.\n");

  if (wbb_params.flags & WBB_1ST_MAX_TOTAL_RATES)
    printf("# Using total rates analysis for 1st maximum (if included).\n");
  if (wbb_params.flags & WBB_2ND_MAX_TOTAL_RATES)
    printf("# Using total rates analysis for 2nd maximum (if included).\n");

  RestrictWBBEnergyWindow(&wbb_params);  

  return 0;
}


// -------------------------------------------------------------------------
int main(int argc, char *argv[])
// -------------------------------------------------------------------------
//                            M A I N   P R O G R A M
// -------------------------------------------------------------------------
{
  int status;
  true_theta12 = asin(sqrt(0.32));
  true_theta13 = asin(sqrt(0.092))/2.0;
  true_theta23 = M_PI/4;
  true_deltacp = 3.0*M_PI/2.0;
  true_sdm = 7.6e-5;
  true_ldm = 2.47e-3;

  // NuFit 3.2
  true_theta12 = 33.62 * M_PI/180.;
  true_theta13 =  8.54 * M_PI/180.;
  true_theta23 = 47.2  * M_PI/180.;
  true_deltacp = 3.0*M_PI/2.0;
  true_sdm     = 7.4e-5;
  true_ldm     = 2.494e-3;

  time_t start_time, end_time;
  time(&start_time);
  printf("# GLoBES neutrino oscillation simulation\n");
  printf("# --------------------------------------\n");
  printf("#\n");
  printf("# Run started on %s", asctime(localtime(&start_time)));
  printf("#\n");
  printf("# Command line:   ");
  for (int i=0; i < argc; i++)
    printf("%s ", argv[i]);
  printf("\n");
  printf("#\n");

  // Initialize MPI 
  #ifdef NU_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    printf("# Initializing thread %d of %d.\n", mpi_rank+1, mpi_size);
  #elif defined NU_PSEUDO_MPI
    /* Read environment variables NU_RANK and NU_SIZE that tell us which part
     * of the job this instance of the program should carry out */
    if (!getenv("NU_RANK") || !getenv("NU_SIZE"))
    {
      mpi_size = 1;
      mpi_rank = 0;
      printf("# NU_RANK and/or NU_SIZE not defined. Assuming no parallelization.\n");
    }
    else
    {
      mpi_size = atoi(getenv("NU_SIZE"));
      mpi_rank = atoi(getenv("NU_RANK"));
      if (mpi_rank >= mpi_size)
      {
        fprintf(stderr, "Error: NU_RANK must be smaller than NU_SIZE.\n");
        return -1;
      }
      printf("# Computing slice %d of %d in parallel environment.\n", mpi_rank+1, mpi_size);
    }
  #endif

  // Switch off GSL error handling
  gsl_set_error_handler_off();

  // Parse command line arguments
  struct argp argp = { cmdline_options, parse_opt, argp_option_doc, argp_doc };
  if ((status=argp_parse(&argp, argc, argv, 0, NULL, NULL) != 0))
    return status;

  // Check if all required command line arguments were given
  if (action < 0)
  {
    fprintf(stderr, "Please use -a option to specifiy the type of simulation.\n");
    return -1;
  }
  if (n_exps <= 0)
  {
    fprintf(stderr, "Please use -e option to specifiy the experiments to simulate.\n");
    return -2;
  }
  if (action != NU_ACTION_PARAM_SCAN  &&
      action != NU_ACTION_MCMC  &&
      action != NU_ACTION_PROB_TABLE  &&
      action != NU_ACTION_EXPOSURE_SCAN  &&  n_scan_params > 0)
  {
    fprintf(stderr, "Warning: -p is used only when -a PARAM_SCAN is given.\n");
    fprintf(stderr, "Will ignore all -p options\n");
  }
//  else if (action == NU_ACTION_PARAM_SCAN  &&  n_scan_params <= 0)
//  {
//    fprintf(stderr, "Error: -a PARAM_SCAN requires at least one -p specification.\n");
//    return -3;
//  }
  if (action!=NU_ACTION_PARAM_SCAN && action!=NU_ACTION_EXPOSURE_SCAN && n_min_params > 0)
  {
    fprintf(stderr, "Warning: -m is used only when -a PARAM_SCAN is given.\n");
    fprintf(stderr, "Will ignore all -m options\n");
  }

  // Initialize matter profile data 
  if (LoadPREMProfile("data/prem-profile.dat") != 0)
  {
    fprintf(stderr, "Failed to load PREM profile.\n");
    return -1;
  }

  // Initialize libglobes 
  glbInit(argv[0]);
  glbSelectMinimizer(GLB_MIN_POWELL); // Parts of code work ONLY with GLB_MIN_POWELL !!! 

  // Initialize and register non-standard probability engine. This has to be done
  // before any calls to glbAllocParams() or glbAllocProjections()
  int n_params = 0;
  if (n_flavors == 3)
  {
    int default_rotation_order[][2] = { {2,3}, {1,3}, {1,2} };
    int default_phase_order[] = { -1, 0, -1};
    n_params = snu_init_probability_engine(n_flavors,
                                           default_rotation_order, default_phase_order);
//    glbRegisterProbabilityEngine(6*SQR(n_flavors)-n_flavors+1, &snu_probability_matrix,
//      &snu_set_oscillation_parameters, &snu_get_oscillation_parameters, NULL);
  }
  else if (n_flavors == 4)
  {
    int default_rotation_order[][2] = { {3,4}, {2,4}, {1,4}, {2,3}, {1,3}, {1,2} };
    int default_phase_order[] = { -1,  1, -1, -1,  0,   2};
    n_params = snu_init_probability_engine(n_flavors,
                                           default_rotation_order, default_phase_order);
  }
  else if (n_flavors == 5)
  {
    int default_rotation_order[][2] = { {4, 5}, {3,5}, {2,5}, {1,5}, {3,4}, {2,4}, {1,4},
                                        {2,3}, {1,3}, {1,2} };
    int default_phase_order[] = { -1, -1, -1,  2, -1, -1,  1, -1,  0, -1};
    n_params = snu_init_probability_engine(n_flavors,
                                           default_rotation_order, default_phase_order);
  }
  glbRegisterProbabilityEngine(n_params, &snu_probability_matrix,
    &snu_set_oscillation_parameters, &snu_get_oscillation_parameters, NULL);

#ifdef NU_USE_NUSQUIDS
  glbRegisterNuSQuIDSEngine(nu_hook_probability_matrix_nusquids);
  printf("# Using NuSQuIDS probability engine. Majorana=%d\n", OSC_DECAY_MAJORANA);
#else
  printf("# Using SNU probability engine.\n");
#endif

  for (int j=0; j < glbGetNumOfOscParams(); j++)
  {
    glbSetParamName(snu_param_strings[j], j);
#ifdef NU_USE_MONTECUBES
    mcb_setVarName(j, snu_param_strings[j]);
#endif
  }
  
  // Initialize user-defined chi^2 functions (chiMB_init has to be called after the
  // number of oscillation parameters has been fixed
  int minos_nc = MINOS_NC;
  int minos_cc = MINOS_CC;
  glbDefineChiFunction(&chiNOvA,         18, "chiNOvA",          NULL);
  glbDefineChiFunction(&chiWBB_WCfast,   10, "chiWBB_WCfast",    &wbb_params);
  glbDefineChiFunction(&chiWBB_LAr,      10, "chiWBB_LAr",       &wbb_params);
  glbDefineChiFunction(&chiDCNorm,        5, "chiDCNorm",        NULL);
  glbDefineChiFunction(&chiKamLAND,       1, "chiKamLAND",       NULL);

  glbDefineChiFunction(&chiMINOS_2016,    5, "chiMINOS-NC-2016", &minos_nc);
  glbDefineChiFunction(&chiMINOS_2016,    5, "chiMINOS-CC-2016", &minos_cc);
  glbDefineChiFunction(&chiMINOS_2016,    0, "chiMINOS-nosys-NC-2016",&minos_nc);
  glbDefineChiFunction(&chiMINOS_2016,    0, "chiMINOS-nosys-CC-2016",&minos_cc);

  glbDefineChiFunction(&chiMINOS_2011,    5, "chiMINOS-NC-2011", &minos_nc);
  glbDefineChiFunction(&chiMINOS_2011,    5, "chiMINOS-CC-2011", &minos_cc);
  glbDefineChiFunction(&chiMINOS_2011,    0, "chiMINOS-nosys-NC-2011",&minos_nc);
  glbDefineChiFunction(&chiMINOS_2011,    0, "chiMINOS-nosys-CC-2011",&minos_cc);

  glbDefineChiFunction(&chiMINOS_2010,    5, "chiMINOS-NC-2010", &minos_nc);
  glbDefineChiFunction(&chiMINOS_2010,    5, "chiMINOS-CC-2010", &minos_cc);
  glbDefineChiFunction(&chiMINOS_2010,    0, "chiMINOS-nosys-NC-2010",&minos_nc);
  glbDefineChiFunction(&chiMINOS_2010,    0, "chiMINOS-nosys-CC-2010",&minos_cc);

  glbDefineChiFunction(&chiT2K_FDonly,    7, "chiT2K",       NULL);
  glbDefineChiFunction(&chiT2K_FDonly,    0, "chiT2K-nosys", NULL);

  glbDefineChiFunction(&chiLSNDspectrum,  2, "chiLSNDspectrum",  NULL);
//  chiMB_init(); // for 2010 version of MiniBooNE code
//  glbDefineChiFunction(&chiMBanti_nu2010, 0, "chiMBanti_nu2010", NULL);
  glbDefineChiFunction(&chi_E776,         4, "chi_E776",         NULL);
  glbDefineChiFunction(&chi_E776_rates,   4, "chi_E776_rates",   NULL);
  glbDefineChiFunction(&chi_karmen_c12,   1, "chi_karmen_c12",   NULL);
  glbDefineChiFunction(&chi_karmen_c12_JR,1, "chi_karmen_c12_JR",NULL);
  glbDefineChiFunction(&chi_lsnd_c12,     1, "chi_lsnd_c12",     NULL);

  // Evaluate first bunch of environment variables to determine mode of operation, etc. 
  if (eval_env() < 0)
    return -1;

  // Load experiments 
  if (load_exps(n_exps, exps) < 0)
    return -2;

  // Initialize parameter and projection vector(s) 
  true_values     = glbAllocParams();
  test_values     = glbAllocParams();
  central_values  = glbAllocParams();
  input_errors    = glbAllocParams();
  proj            = glbAllocProjection();
  prescan_proj    = glbAllocProjection();

  // Define oscillation parameters 
  for (int i=0; i < glbGetNumOfOscParams(); i++)
    glbSetOscParams(true_values,0.0,i);
  glbDefineParams(true_values,true_theta12,true_theta13,true_theta23,
                              true_deltacp,true_sdm,true_ldm);
  glbSetDensityParams(true_values, 1.0, GLB_ALL);

//  // Define input errors
//  //
//  // External priors are dealt with in the following way:
//  // - by default th12 and dm21 carry priors
//  // - if any -c option is given, only the parameters given there will have a prior
//  // - th12 and dm21 priors are omitted if solar neutrinos are included in the fit
//  // - Right now only priors on th12, dm21, th13, th23, dm31 are supported
////  prior_th12    = 0.05 * true_theta12;
////  prior_th13    = 0.0;
////  prior_th23    = 0.0;
////  prior_deltacp = 0.0;
////  prior_sdm     = 0.05 * true_sdm;
////  prior_ldm     = 0.0;
//  prior_th12 = prior_th13 = prior_th23 = prior_deltacp = prior_sdm = prior_ldm = 0.0;
//  for (int i=0; i < glbGetNumOfOscParams(); i++)
//    glbSetOscParams(input_errors, 0.0, i);
//  for (int i=0; i < n_constrained_params; i++)
//  {
////    if (strcmp(constrained_params[i], "TH12") == 0)
////      prior_th12    = 0.05 * true_theta12;
////    else if (strcmp(constrained_params[i], "DM21") == 0)
////      prior_sdm     = 0.05 * true_sdm;
////    else if (strcmp(constrained_params[i], "TH23") == 0)
////      prior_th23    = 0.09 * true_theta23;
////    else if (strcmp(constrained_params[i], "DM31") == 0)
////      prior_ldm     = 0.09e-3; // 0.05 * true_ldm;
//
//    // NuFit 3.2
//    if (strcmp(constrained_params[i], "TH12") == 0)
//      prior_th12    = 0.78 * M_PI/180.;
//    else if (strcmp(constrained_params[i], "DM21") == 0)
//      prior_sdm     = 0.21e-5;
//    else if (strcmp(constrained_params[i], "TH13") == 0)
//      prior_th13    = 0.15 * M_PI/180.;
//    else if (strcmp(constrained_params[i], "TH23") == 0)
//      prior_th23    = 3.9  * M_PI/180.;
//    else if (strcmp(constrained_params[i], "DM31") == 0)
//      prior_ldm     = 0.033e-3;
//
//    // rough constraints on sterile neutrinos
//    else if (strcmp(constrained_params[i], "Ue4") == 0)
//      glbSetOscParamByName(input_errors, 0.164, "Ue4");
//                                 // solar curve from 1803.10661, fig. 3, 1 dof
//    else if (strcmp(constrained_params[i], "Um4") == 0)
//      glbSetOscParamByName(input_errors, 0.0596, "Um4");
//                                 // based on noevid_noIC2 run, 1 dof
//  }
//  if (ext_flags & EXT_SOLAR) // Use prior on solar parameters only if solar data not included in fit
//  {
//    printf("# *** Ignoring **** priors on th12 and m21 since solar neutrinos "
//           "are included in fit.\n");
//    prior_th12 = prior_sdm = 0.0;
//  }
//  glbDefineParams(input_errors, prior_th12, prior_th13, prior_th23, prior_deltacp,
//                  prior_sdm, prior_ldm);
//  glbSetDensityParams(input_errors, 0.05, GLB_ALL);


  // Evaluate true parameters given on the command line (-t option)
  if (true_param_def  &&  eval_true_params(true_param_def) < 0)
    return -2;

  // Evaluate external priors given on the command line (-c option)
  for (int i=0; i < glbGetNumOfOscParams(); i++)
    glbSetOscParams(input_errors, 0.0, i);
  glbSetDensityParams(input_errors, 0.05, GLB_ALL);
  if (input_error_def  &&  eval_input_errors(input_error_def) < 0)
    return -3;

  // Copy true params to other parameter vectors AFTER NSI params have been
  // read from the environment variables
  glbCopyParams(true_values, test_values);
  glbCopyParams(true_values, central_values);

  // --------------------------------------------------------------------
  // The following code is used for debugging the implementation of the
  // 4- and 5-flavor parameterizations TODO: Remove
//  {
//    gsl_matrix_complex *U = snu_get_U();
//    glbSetOscParamByName(true_values, 0.9*M_PI/4, "TH23");
//    glbSetOscParamByName(true_values, 0.33, "TH13");
//    glbSetOscParamByName(true_values, 0.5, "TH14");
//    glbSetOscParamByName(true_values, 0.2, "TH24");
//    glbSetOscParamByName(true_values, 0.75, "TH15");
//    glbSetOscParamByName(true_values, 0.77, "TH25");
//    glbSetOscParamByName(true_values, M_PI/3, "DELTA_0");
//    glbSetOscParamByName(true_values, M_PI/4, "DELTA_1");
//    glbSetOscParamByName(true_values, M_PI/5, "DELTA_2");
//    glbSetOscillationParameters(true_values);
//    snu_print_gsl_matrix_complex(U);
//
//    glbSetOscParamByName(true_values, 0.0, "TH14");
//    glbSetOscParamByName(true_values, 0.0, "TH24");
//    glbSetOscParamByName(true_values, gsl_complex_abs(gsl_matrix_complex_get(U, 0, 3)), "Ue4");
////    glbSetOscParamByName(true_values, gsl_complex_abs(gsl_matrix_complex_get(U, 1, 3)), "Um4");
//    glbSetOscParamByName(true_values, gsl_complex_abs(gsl_matrix_complex_get(U, 1, 3))
//                            * gsl_complex_abs(gsl_matrix_complex_get(U, 0, 3)), "Ue4Um4");
//    glbSetOscParamByName(true_values, 0.0, "TH15");
//    glbSetOscParamByName(true_values, 0.0, "TH25");
//    glbSetOscParamByName(true_values, gsl_complex_abs(gsl_matrix_complex_get(U, 0, 4)), "Ue5");
////    glbSetOscParamByName(true_values, gsl_complex_abs(gsl_matrix_complex_get(U, 1, 4)), "Um5");
//    glbSetOscParamByName(true_values, gsl_complex_abs(gsl_matrix_complex_get(U, 1, 4))
//                            * gsl_complex_abs(gsl_matrix_complex_get(U, 0, 4)), "Ue5Um5");
//    printf("\n");
//    glbSetOscillationParameters(true_values);
//    snu_print_gsl_matrix_complex(U);
//
//    glbPrintParams(stdout, true_values);
//    glbGetOscillationParameters(true_values);
//    glbPrintParams(stdout, true_values);
//
//    exit(1);
//  }
  // --------------------------------------------------------------------

  // Print true parameter values and various kinds of meta information 
  printf("# Using %d flavor model\n", n_flavors);
  printf("# Simulating the following experiments:\n");
  for (int i=0; i < n_exps; i++)
    printf("#   %s\n", exps[i]);
  printf("# Including the following AEDL files:\n");
  for (int i=0; i < glb_num_of_exps; i++)
  {
    printf("#   %s (version %s)", glbGetFilenameOfExperiment(i),
           glbVersionOfExperiment(i));
    #ifdef NU_USE_NUSQUIDS
      if (strcasecmp(glbGetFilenameOfExperiment(i), "MBneutrino200.glb") == 0)
      {
//        extern double Eminnu, Eminbar;
//        printf("; nu mode E > %g, nu-bar mode E > %g", Eminnu, Eminbar);
        extern double Emin;
        printf("; E > %g MeV", Emin);
      }
    #endif
    printf("\n");
  }
  printf("#\n");
  print_aedl_variables();
  printf("#\n");
  printf("# External analysis routines included:\n");
  if (ext_flags & EXT_MB)
    printf("#   Thomas' MiniBooNE code (neutrino run, E > 475 MeV)\n");
  if (ext_flags & EXT_MB_300)
    printf("#   Thomas' MiniBooNE code (neutrino run, E > 300 MeV)\n");
  if (ext_flags & EXT_MBANTI)
    printf("#   Thomas' MiniBooNE code (anti-neutrino run, E > 475 MeV)\n");
  if (ext_flags & EXT_MBANTI_200)
    printf("#   Thomas' MiniBooNE code (anti-neutrino run, E > 200 MeV)\n");
  if (ext_flags & EXT_KARMEN)
    printf("#   KARMEN\n");
  if (ext_flags & EXT_LSND)
    printf("#   LSND (Thomas' code)\n");
  if (ext_flags & EXT_LSND_IVAN)
    printf("#   LSND (Ivan's code)\n");
  if (ext_flags & EXT_KARMEN_IVAN)
    printf("#   KARMEN (Ivan's code)\n");
  if (ext_flags & EXT_MB_JK)
    printf("#   MiniBooNE (Joachim's code)\n");

  if (ext_flags & EXT_REACTORS)
  {
    printf("#   \\nu_e disappearance searches: \n");
    #ifdef USE_SBL
      printf("#     Bugey, ROVNO, Goesgen, ILL, Kranoyarsk, SRP, Rovno rates\n");
    #endif
    #ifdef USE_CHOOZ
      printf("#     Chooz\n");
    #endif
    #ifdef USE_PV
      printf("#     Palo Verde\n");
    #endif
    #ifdef USE_KAML
      printf("#     KamLAND\n");
    #endif
    #ifdef USE_DC
      printf("#     Double Chooz\n");
    #endif
    #ifdef USE_DB
      printf("#     Daya Bay (sterile neutrino code)\n");
    #endif
    #ifdef USE_DB_3F
      printf("#     Daya Bay (3-flavor code)\n");
    #endif
    #ifdef USE_RENO
      printf("#     RENO\n");
    #endif
    #ifdef USE_BUGEY_SP
      printf("#     Bugey spectrum\n");
    #endif
    #ifdef USE_DANSS
      printf("#     DANSS\n");
    #endif
    #ifdef USE_NEOS
      printf("#     NEOS\n");
    #endif
    #ifdef USE_GAL
      printf("#     Gallium\n");
    #endif

    #ifdef USE_NEOS_DB_ALVARO
      printf("#     Combined NEOS+ Daya Bay analysis\n");
    #endif
    #ifdef USE_DB_ALVARO
      printf("#     Daya Bay (Alvaro's implementation)\n");
    #endif
  }

  if (ext_flags & EXT_NOMAD)
    printf("#   NOMAD\n");
  if (ext_flags & EXT_CDHS)
    printf("#   CDHS\n");
  if (ext_flags & EXT_ATM_TABLE)
    printf("#   SuperK Atmospheric neutrinos (tabulated chi^2)\n");
  if (ext_flags & EXT_ATM_COMP)
    printf("#   SueprK Atmospheric neutrinos (Michele's simulation)\n");
  if (ext_flags & EXT_DEEPCORE)
    printf("#   IceCube Deep Core (Michele's simulation)\n");
  if (ext_flags & EXT_SOLAR)
    printf("#   Solar neutrinos (Michele's simulation)\n");
  if (ext_flags & EXT_DECAY_KINEMATICS)
    printf("#   Beta decay kinematics\n");
  if (ext_flags & EXT_FREE_STREAMING)
    printf("#   Cosmology (free-streaming constraint)\n");
  printf("#\n");

  printf("# True oscillation parameters (only nonzero params shown):\n");
  for (int i=0; i < glbGetNumOfOscParams(); i++)
    if (glbGetOscParams(true_values,i) != 0.0)
      printf("#   %s = %g\n", snu_param_strings[i], glbGetOscParams(true_values,i));

  printf("# External priors:\n");
  for (int i=0; i < glbGetNumOfOscParams(); i++)
    if (glbGetOscParams(input_errors,i) != 0.0)
      printf("#   %s = %g +- %g (%g%%)\n", snu_param_strings[i],
             glbGetOscParams(central_values,i), glbGetOscParams(input_errors,i),
             100.0*glbGetOscParams(input_errors,i)/glbGetOscParams(central_values,i));

  printf("#\n");
  printf("# Miscellaneous options:\n");
  printf("#   Constraint \\eps^s = \\eps^{d\\dag}: %s\n", use_nsi_constraints ? "YES" : "NO");
  printf("#   New reactor fluxes (1101.2663):   %s\n",
         ns_reactor::old_new_main==NEW ? "YES" : "NO");
  if (n_flavors >= 5)
  {
#ifdef Ip3pI
    printf("#   5-flavor scheme:                  1+3+1\n");
#else
    printf("#   5-flavor scheme:                  3+2\n");
#endif
  }
  printf("#\n");

  printf("# Parameter scan:\n");
  if (n_scan_params > 0)
  {
    for (int i=0; i < n_scan_params; i++)
      printf("#   %10s, %10.7g, %10.7g, %5d  0x%lx\n", scan_params[i],
             scan_p_min[i], scan_p_max[i], scan_p_steps[i], scan_p_flags[i]);
  }
  else
    printf("#   Computing only best fit point\n");
  printf("# Marginalizing over:\n");
  printf("#   ");
  for (int i=0; i < n_min_params; i++)
    printf("%s ", min_params[i]);
  printf("\n");
  printf("# Prescan for degeneracy finding:\n");
  for (int i=0; i < n_prescan_params; i++)
    printf("#   %10s, %10.7g, %10.7g, %5d  0x%lx\n", prescan_params[i],
           prescan_p_min[i], prescan_p_max[i], prescan_p_steps[i], prescan_p_flags[i]);
  printf("# Degfinder flags: 0x%lx\n", default_degfinder_flags);
  printf("#\n");

  if (load_exps(n_exps, exps) < 0)  // Load experiments 
    return -5;
  ext_init(ext_flags);              // Initialize external codes 

  switch (action)
  {
    // Print event rates 
    case NU_ACTION_SPECTRUM:
      printf("# Printing event rates.\n");
      printf("#\n");
      print_rates(ext_flags);
      break;

    // General parameter scan
    case NU_ACTION_PARAM_SCAN:
    {
      printf("# General parameter scan\n");
      printf("#\n");

      // User-defined prior function to include external input
      prior_params pp;
      pp.ext_flags = ext_flags;
      pp.n_scan_params = n_scan_params;
      for (int i=0; i < n_scan_params; i++)
        pp.scan_params[i] = glbFindParamByName(scan_params[i]);
      memcpy(pp.scan_p_flags, scan_p_flags, n_scan_params*sizeof(*scan_p_flags));
      glbRegisterPriorFunction(&my_prior, NULL, NULL, &pp);

      // Do the scan
      param_scan("", n_scan_params, scan_params, scan_p_min, scan_p_max, scan_p_steps,
                 scan_p_flags, n_min_params, min_params, n_prescan_params, prescan_params,
                 prescan_p_min, prescan_p_max, prescan_p_steps, prescan_p_flags);

      // Reset prior function to default
      glbRegisterPriorFunction(NULL, NULL, NULL, NULL);
      break;
    }

    // Markov Chain Monte Carlo
    case NU_ACTION_MCMC:
    {
      printf("# Markov Chain Monte Carlo\n");
      printf("#\n");

      // Sample over parameters indicated as "scan" parameters (-p) and as
      // "minimization" parameters (-m) on the command line
      int n_mcmc_params = n_scan_params + n_min_params;
      char *mcmc_params[n_mcmc_params];
      double mcmc_p_min[n_mcmc_params], mcmc_p_max[n_mcmc_params];
      unsigned long mcmc_p_flags[n_mcmc_params];
      for (int i=0; i < n_scan_params; i++)
        mcmc_params[i] = scan_params[i];
      for (int i=0; i < n_min_params; i++)
      {
        mcmc_params[n_scan_params+i] = min_params[i];
        mcmc_p_min[n_scan_params+i] = 0.0;
        mcmc_p_max[n_scan_params+i] = 0.0;
        mcmc_p_flags[n_scan_params+i] = 0;
      }
      memcpy(mcmc_p_min, scan_p_min, sizeof(*mcmc_p_min) * n_scan_params);
      memcpy(mcmc_p_max, scan_p_max, sizeof(*mcmc_p_max) * n_scan_params);
      memcpy(mcmc_p_flags, scan_p_flags, sizeof(*mcmc_p_flags) * n_scan_params);

      // User-defined prior function to include external input
      prior_params pp;
      pp.ext_flags = ext_flags;
      pp.n_scan_params = n_mcmc_params;
      for (int i=0; i < n_mcmc_params; i++)
        pp.scan_params[i] = glbFindParamByName(mcmc_params[i]);
      memcpy(pp.scan_p_flags, mcmc_p_flags, n_mcmc_params*sizeof(*mcmc_p_flags));
      glbRegisterPriorFunction(&my_prior, NULL, NULL, &pp);

      // Run Markov chains
      mcmc_deg(output_file, n_mcmc_params, mcmc_params, mcmc_p_min, mcmc_p_max,
                 mcmc_p_flags, n_prescan_params, prescan_params,
                 prescan_p_min, prescan_p_max, prescan_p_steps, prescan_p_flags);

      // Reset prior function to default
      glbRegisterPriorFunction(NULL, NULL, NULL, NULL);
      break;
    }

    // Tabulate oscillation probabilities
    case NU_ACTION_PROB_TABLE:
    {
      printf("# Tabulating oscillation probabilities\n");
      printf("#\n");

      for (int i=0; i < glb_num_of_exps; i++)
      {
        // Prepare data structures
        char fname[1024];
        sprintf(fname, "%s-%s", output_file, basename(glb_experiment_list[i]->filename));
        struct snu_probability_table p;
        p.default_values = glbAllocParams();
        glbCopyParams(true_values, p.default_values);
        p.n_p = n_scan_params;
        for (int j=0; j < p.n_p; j++)
        {
          p.params[j]  = strdup(scan_params[j]);
          p.p_min[j]   = scan_p_min[j];
          p.p_max[j]   = scan_p_max[j];
          p.p_steps[j] = scan_p_steps[j];
          p.p_flags[j] = scan_p_flags[j];
        }

        // Compute probability table
        snu_compute_probability_table(i, &p, fname);


        // FIXME FIXME Read probability table again
        struct snu_probability_table *p2 = snu_alloc_probability_table();
        snu_load_probability_table(fname, p2);
      }
      break;
    }

    // Scan over parameters and over exposure
    case NU_ACTION_EXPOSURE_SCAN:
    {
      printf("# EXPOSURE SCAN NOT IMPLEMENTED YET!\n");
      break;
    }

//    // Check Thomas' best fit points
//    case NU_ACTION_CHECK_BF:
//      printf("# Checking Thomas' best fit point\n");
//      printf("#\n");
//      checkBF(n_flavors);
//      break;
  }

  time(&end_time);
  printf("# Run ended on %s", asctime(localtime(&end_time)));
  printf("# Running time (wallclock): %g sec\n\n", difftime(end_time, start_time));
 
  // Destroy parameter and projection vector(s) 
  glbFreeParams(true_values);
  glbFreeParams(test_values);
  glbFreeParams(central_values);
  glbFreeParams(input_errors);
  glbFreeProjection(prescan_proj);
  glbFreeProjection(proj);

  #ifdef NU_USE_NUSQUIDS
    chiMB_clear();
  #endif
  chiMB_jk_clear();
    
  // Cleanup MPI 
  #ifdef NU_MPI
    MPI_Finalize();
  #endif

  exit(0);
}

