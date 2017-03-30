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
#include <argp.h>

#include <globes/globes.h>   // GLoBES library 
#include "glb_types.h"
#include "glb_error.h"
#include "const.h"
#include "snu.h"
#include "nu.h"
#include "sbl/def-reactors.h"

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
static char *true_param_def=NULL; // The definition of true osc. params given via the -t option

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
double true_theta12, prior_th12;
double true_theta13, prior_th13;
double true_theta23, prior_th23;
double true_deltacp, prior_deltacp;
double true_sdm, prior_sdm;
double true_ldm, prior_ldm;
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
extern int old_new_main;

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


// Definitions for argp
const char *argp_program_version = "nu 0.1";
const char *argp_program_bug_address = "<jkopp@fnal.gov>";
static char argp_doc[] = "nu neutrino oscillation simulation";
static char argp_option_doc[] = "[options]";

// The actual list of command line options we accept
static struct argp_option cmdline_options[] = {
  {"flavors",    'f',"NUMBER",0,"Number of flavors to use. Can be 3, 4, or 5)" },
  {"action",     'a',"ACTION",0,"Which parameters to scan (allowed values defined in const.c)" },
  {"experiments",'e',"EXPS",  0,"Which experiments to include in the fit" },
  {"parameter",  'p',"PARAMS",0,"Parameters to scan: <param_name>,<min>,<max>,<steps>,[<flags>]"},
  {"prescan",    's',"PARAMS",0,"Parameters for prescan in degfinder: "
                                "<param_name>,<min>,<max>,<steps>,[<flags>]"},
  {"minimize",   'm',"PARAMS",0,"Parameters to marginalize: <param_name>[,<param_name[, ...]]"},
  {"true_params",'t',"PARAMS",0,"True oscillation parameters: \"NAME=VALUE, ...\""},
  {"verbose",    'v',NULL,    0,"Show debug output (use multiple times for more)"},
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
          scan_p_min[n_scan_params] = atof(p);
          p = strtok(NULL, delim);
          if (p)
          {
            scan_p_max[n_scan_params] = atof(p);
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
          prescan_p_min[n_prescan_params] = atof(p);
          p = strtok(NULL, delim);
          if (p)
          {
            prescan_p_max[n_prescan_params] = atof(p);
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
        char s[strlen(true_param_def)+strlen(arg)];
        strcpy(s, true_param_def);
        free(true_param_def);
        true_param_def = strdup(strcat(s, arg));
      }
      break;
    }

    // -------------------------------------------------
    case 'v':
    {
      debug_level++;
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

    // Double Chooz (near + far)
    else if (strcasecmp(exps[i], "DCHOOZ") == 0)
    {
      glbInitExperiment("D-Chooz_near.glb",&glb_experiment_list[0],&glb_num_of_exps);
      glbInitExperiment("D-Chooz_far.glb",&glb_experiment_list[0],&glb_num_of_exps);
      L_opt[0] = 0;
      EXP_REACTOR_NEAR = 0;
      EXP_REACTOR_FAR  = 1;
      prior_th23 = 0.10 * true_theta23;
      prior_ldm  = 0.05 * true_ldm;
    }

    // KamLAND
    else if (strcasecmp(exps[i], "KAMLAND") == 0)
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

    // MINOS Neutral Current analysis (http://arxiv.org/abs/1001.0336)
    else if (strcasecmp(exps[i], "MINOS_NC") == 0)
    {
      glbInitExperiment("MINOS-NC-near.glb",&glb_experiment_list[0],&glb_num_of_exps);
      glbInitExperiment("MINOS-NC-far.glb",&glb_experiment_list[0],&glb_num_of_exps);
      L_opt[0] = 0;
    }

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
    }
  } // for(i)

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
  true_theta13 = 0.0;
  true_theta23 = M_PI/4;
  true_deltacp = 3.0*M_PI/2.0;
  true_sdm = 7.6e-5;
  true_ldm = 2.4e-3;

  prior_th12    = 0.05 * true_theta12;
  prior_th13    = 0.0;
  prior_th23    = 0.0;
  prior_deltacp = 0.0;
  prior_sdm     = 0.05 * true_sdm;
  prior_ldm     = 0.0;

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
  if (action!=NU_ACTION_PARAM_SCAN && action!=NU_ACTION_EXPOSURE_SCAN  &&  n_scan_params > 0)
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

  time_t start_time;
  time(&start_time);
  printf("# GLoBES neutrino oscillation simulation\n");
  printf("# --------------------------------------\n");
  printf("#\n");
  printf("# Run started on %s", asctime(localtime(&start_time)));
  printf("#\n");

  // Initialize MPI 
  #ifdef NU_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    printf("# Initializing thread %d of %d\n.", mpi_rank+1, mpi_size);
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

  // Initialize matter profile data 
  if (LoadPREMProfile("data/prem-profile.dat") != 0)
  {
    fprintf(stderr, "Failed to load PREM profile.\n");
    return -1;
  }

  // Initialize libglobes 
  setenv("GLB_PATH", "glb/nova:glb/t2k:glb/bbeam:glb/dchooz:glb/wbb_wc:"
                     "glb/wbb_lar:glb/nufact:glb/minos-nc:glb/miniboone:"
                     "glb/kamland:glb/lsnd", 1);
  glbInit(argv[0]);
  glbSelectMinimizer(GLB_MIN_POWELL); // Parts of code work ONLY with GLB_MIN_POWELL !!! 

  // Initialize and register non-standard probability engine. This has to be done
  // before any calls to glbAllocParams() or glbAllocProjections()
  if (n_flavors == 3)
  {
    int default_rotation_order[][2] = { {2,3}, {1,3}, {1,2} };
    int default_phase_order[] = { -1, 0, -1};
    snu_init_probability_engine(n_flavors, default_rotation_order, default_phase_order);
    glbRegisterProbabilityEngine(6*SQR(n_flavors)-n_flavors+1, &snu_probability_matrix,
      &snu_set_oscillation_parameters, &snu_get_oscillation_parameters, NULL);
  }
  else if (n_flavors == 4)
  {
    int default_rotation_order[][2] = { {3,4}, {2,4}, {1,4}, {2,3}, {1,3}, {1,2} };
    int default_phase_order[] = { -1,  1, -1, -1,  0,   2};
    snu_init_probability_engine(n_flavors, default_rotation_order, default_phase_order);
    glbRegisterProbabilityEngine(6*SQR(n_flavors)-n_flavors+3, &snu_probability_matrix,
      &snu_set_oscillation_parameters, &snu_get_oscillation_parameters, NULL);
  }
  else if (n_flavors == 5)
  {
    int default_rotation_order[][2] = { {4, 5}, {3,5}, {2,5}, {1,5}, {3,4}, {2,4}, {1,4},
                                        {2,3}, {1,3}, {1,2} };
    int default_phase_order[] = { -1, -1, -1,  2, -1, -1,  1, -1,  0, -1};
    snu_init_probability_engine(n_flavors, default_rotation_order, default_phase_order);
    glbRegisterProbabilityEngine(6*SQR(n_flavors)-n_flavors+5, &snu_probability_matrix,
      &snu_set_oscillation_parameters, &snu_get_oscillation_parameters, NULL);
  }

  for (int j=0; j < glbGetNumOfOscParams(); j++)
    glbSetParamName(snu_param_strings[j], j);
  
  // Initialize user-defined chi^2 functions (chiMB_init has to be called after the
  // number of oscillation parameters has been fixed
  glbDefineChiFunction(&chiT2K,          16, "chiT2K",           NULL);
  glbDefineChiFunction(&chiNOvA,         18, "chiNOvA",          NULL);
  glbDefineChiFunction(&chiWBB_WCfast,   10, "chiWBB_WCfast",    &wbb_params);
  glbDefineChiFunction(&chiWBB_LAr,      10, "chiWBB_LAr",       &wbb_params);
  glbDefineChiFunction(&chiDCNorm,        5, "chiDCNorm",        NULL);
  glbDefineChiFunction(&chiKamLAND,       1, "chiKamLAND",       NULL);
  glbDefineChiFunction(&chiMINOS,         5, "chiMINOS",         NULL);
  glbDefineChiFunction(&chiMINOS,         0, "chiMINOS-nosys",   NULL);
  glbDefineChiFunction(&chiLSNDspectrum,  2, "chiLSNDspectrum",  NULL);
  chiMB_init();
  glbDefineChiFunction(&chiMBanti_nu2010, 0, "chiMBanti_nu2010", NULL);

  // Evaluate first bunch of environment variables to determine mode of operation, etc. 
  if (eval_env() < 0)
    return -1;

  // Load experiments 
  if (load_exps(n_exps, exps) < 0)
    return -2;

  // User-defined prior function to include external input 
  glbRegisterPriorFunction(&my_prior, NULL, NULL, &ext_flags);


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

  // Define input errors   
  for (int i=0; i < glbGetNumOfOscParams(); i++)
    glbSetOscParams(input_errors, 0.0, i);
  glbDefineParams(input_errors, prior_th12, prior_th13, prior_th23, prior_deltacp,
                  prior_sdm, prior_ldm);
  glbSetDensityParams(input_errors, 0.05, GLB_ALL);

  // Evaluate true parameters given on the command line
  if (true_param_def  &&  eval_true_params(true_param_def) < 0)
    return -2;

  // Copy true params to other parameter vectors AFTER NSI params have been
  // read from the environment variables
  glbCopyParams(true_values, test_values);
  glbCopyParams(true_values, central_values);

  // Print true parameter values and various kinds of meta information 
  printf("# Using %d flavor model\n", n_flavors);
  printf("# Simulating the following experiments:\n");
  for (int i=0; i < n_exps; i++)
    printf("#   %s\n", exps[i]);
  printf("# Including the following AEDL files:\n");
  for (int i=0; i < glb_num_of_exps; i++)
  {
    printf("#   %s (version %s)\n", glbGetFilenameOfExperiment(i),
           glbVersionOfExperiment(i));
  }
  printf("#\n");
  print_aedl_variables();
  printf("#\n");
  printf("# External analysis routines included:\n");
  if (ext_flags & EXT_MB)
    printf("#   MiniBooNE (neutrino run)\n");
  if (ext_flags & EXT_MBANTI)
    printf("#   MiniBooNE (anti-neutrino run)\n");
  if (ext_flags & EXT_KARMEN)
    printf("#   KARMEN\n");
  if (ext_flags & EXT_LSND)
    printf("#   LSND\n");
  if (ext_flags & EXT_SBL)
    printf("#   Bugey, Rovno, Krasnoyarsk, ILL, Goesgen, Chooz, Palo Verde\n");
  if (ext_flags & EXT_NOMAD)
    printf("#   NOMAD\n");
  if (ext_flags & EXT_CDHS)
    printf("#   CDHS\n");
  if (ext_flags & EXT_ATM_TABLE)
    printf("#   Atmospheric neutrinos (tabulated chi^2)\n");
  if (ext_flags & EXT_ATM_COMP)
    printf("#   Atmospheric neutrinos (Michele's simulation)\n");
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
  printf("#   New reactor fluxes (1101.2663):   %s\n", old_new_main==NEW ? "YES" : "NO");
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
      print_rates();
      break;

    // General parameter scan
    case NU_ACTION_PARAM_SCAN:
      printf("# General parameter scan\n");
      printf("#\n");
      param_scan("", n_scan_params, scan_params, scan_p_min, scan_p_max, scan_p_steps,
                 scan_p_flags, n_min_params, min_params, n_prescan_params, prescan_params,
                 prescan_p_min, prescan_p_max, prescan_p_steps, prescan_p_flags);
      break;

    // Scan over parameters and over exposure
    case NU_ACTION_EXPOSURE_SCAN:
    {
      printf("# EXPOSURE SCAN NOT IMPLEMENTED YET!\n");
      break;
    }

    // Check Thomas' best fit points
    case NU_ACTION_CHECK_BF:
      printf("# Checking Thomas' best fit point\n");
      printf("#\n");
      checkBF(n_flavors);
      break;
  }

 
  // Destroy parameter and projection vector(s) 
  glbFreeParams(true_values);
  glbFreeParams(test_values);
  glbFreeParams(central_values);
  glbFreeParams(input_errors);
  glbFreeProjection(prescan_proj);
  glbFreeProjection(proj);

  chiMB_clear();
    
  // Cleanup MPI 
  #ifdef NU_MPI
    MPI_Finalize();
  #endif

  exit(0);
}
