/***************************************************************************
 * Non-standard interactions in a neutrino factory with silver channel     *
 ***************************************************************************
 * Author: Joachim Kopp, Toshihiko Ota, Walter Winter                      *
 ***************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <globes/globes.h>
//#include <montecubes/montecubes.h>
#include "glb_types.h"
#include "glb_prior.h"
#include "const.h"
#include "nu.h"
#include "sbl/definitions.h"

// Global variables 
// ---------------- 

// Misc 
extern const double param_bounds[];
extern int scenario;
extern int mode;
extern int nsi_index;
extern const int debug_level;
extern int use_nsi_constraints;

extern int n_nsi_fit_params;
extern int nsi_fit_params[];

extern int default_degfinder_flags;
extern char key_string[100];

// Define standard oscillation parameters (cf. hep-ph/0405172v5) 
extern double true_theta12;
extern double true_theta13;
extern double true_theta23;
extern double true_deltacp;
extern double true_sdm;
extern double true_ldm;
extern int n_flavors;

extern double solar_dm21_min;
extern double solar_dm21_max;

extern int density_corr[];
extern wbb_params_type wbb_params;

// Parameter and projection vectors 
extern glb_params true_values;
extern glb_params test_values;
extern glb_params central_values;
extern glb_params input_errors;
extern glb_params Fit;
extern glb_params Fit_NH;
extern glb_params Fit_IH;
extern glb_params prescan_params1;
extern glb_params prescan_params2;
extern glb_projection proj;
extern glb_projection prescan_proj;

#define MAX_DEG   1000    // Maximum number of degeneracies to expect 
#define MAX_PARAMS  20


/***************************************************************************
 *                     H E L P E R   F U N C T I O N S                     *
 ***************************************************************************/

/*************************************************************************** 
 * Sampling of range of values                                             *
 ***************************************************************************/
double sample(double min, double max, int steps, int i)
{
  if (steps > 0)
    return min + i * (max - min)/steps;
  else
    return min;
}


/*************************************************************************** 
 * Setting density projection flag only for detectos with L > 100 km       *
 ***************************************************************************/
int DecideOnDensityProjection(glb_projection proj)
{
  for (int j=0; j < glb_num_of_exps; j++)
  {
    if(glbGetBaselineInExperiment(j) < 100 || density_corr[j] != j)
      glbSetDensityProjectionFlag(proj, GLB_FIXED, j);
    else
      glbSetDensityProjectionFlag(proj, GLB_FREE, j);
  }

  return 0;
}


/*************************************************************************** 
 * Restrict energy window to 1st max. only or 2nd max. only for LBNE       *
 ***************************************************************************/
int RestrictWBBEnergyWindow(wbb_params_type *wbb_params)
{
  for (int i=0; i < glb_num_of_exps; i++)
  {
    double L = glbGetBaselineInExperiment(i);
    if (L > 100.0) // First/Second maximum only visible in far detectors
    {
      wbb_params->E_1st_min = L*KM * fabs(true_ldm) / (4*M_PI) / GEV;
      break;
    }
  }

  if (wbb_params->flags & WBB_NO_1ST_MAX)
  {
    printf("# Including only 2nd oscillation maximum.\n");
    for (int i=0; i < glb_num_of_exps; i++)
    {
      for (int j=0; j < glbGetNumberOfRules(i); j++)
      {
        double lower, upper;
        glbGetEnergyWindow(i, j, &lower, &upper);
        glbSetEnergyWindow(i, j, lower, MIN(upper, wbb_params->E_1st_min));
      }
    }
  }
  else if (wbb_params->flags & WBB_NO_2ND_MAX)
  {
    printf("# Including only 1st oscillation maximum.\n");
    for (int i=0; i < glb_num_of_exps; i++)
    {
      for (int j=0; j < glbGetNumberOfRules(i); j++)
      {
        double lower, upper;
        glbGetEnergyWindow(i, j, &lower, &upper);
        glbSetEnergyWindow(i, j, MAX(lower, wbb_params->E_1st_min), upper);
      }
    }
  }

  return 0;  
}


/*************************************************************************** 
 * Implement the constraint \eps^s = \eps^{d\dag}                          *
 * If \eps^s is completely zero, \eps^{d\dag} will be copied into \eps^2,  *
 * otherwise vice-versa
 ***************************************************************************/
int implement_nsi_constraints(glb_params p)
{
  if (glbGetOscParamByName(p,"ABS_EPS_S_EE")     == 0.0 && 
      glbGetOscParamByName(p,"ABS_EPS_S_EMU")    == 0.0 && 
      glbGetOscParamByName(p,"ABS_EPS_S_ETAU")   == 0.0 && 
      glbGetOscParamByName(p,"ABS_EPS_S_MUE")    == 0.0 && 
      glbGetOscParamByName(p,"ABS_EPS_S_MUMU")   == 0.0 && 
      glbGetOscParamByName(p,"ABS_EPS_S_MUTAU")  == 0.0 && 
      glbGetOscParamByName(p,"ABS_EPS_S_TAUE")   == 0.0 && 
      glbGetOscParamByName(p,"ABS_EPS_S_TAUMU")  == 0.0 && 
      glbGetOscParamByName(p,"ABS_EPS_S_TAUTAU") == 0.0)
  {
    glbSetOscParamByName(p, glbGetOscParamByName(p,"ABS_EPS_D_EE"),    "ABS_EPS_S_EE");
    glbSetOscParamByName(p, glbGetOscParamByName(p,"ABS_EPS_D_MUE"),   "ABS_EPS_S_EMU");
    glbSetOscParamByName(p, glbGetOscParamByName(p,"ABS_EPS_D_TAUE"),  "ABS_EPS_S_ETAU");
    glbSetOscParamByName(p,-glbGetOscParamByName(p,"ARG_EPS_D_EE"),    "ARG_EPS_S_EE");
    glbSetOscParamByName(p,-glbGetOscParamByName(p,"ARG_EPS_D_MUE"),   "ARG_EPS_S_EMU");
    glbSetOscParamByName(p,-glbGetOscParamByName(p,"ARG_EPS_D_TAUE"),  "ARG_EPS_S_ETAU");
                                                              
    glbSetOscParamByName(p, glbGetOscParamByName(p,"ABS_EPS_D_EMU"),   "ABS_EPS_S_MUE");
    glbSetOscParamByName(p, glbGetOscParamByName(p,"ABS_EPS_D_MUMU"),  "ABS_EPS_S_MUMU");
    glbSetOscParamByName(p, glbGetOscParamByName(p,"ABS_EPS_D_TAUMU"), "ABS_EPS_S_MUTAU");
    glbSetOscParamByName(p,-glbGetOscParamByName(p,"ARG_EPS_D_EMU"),   "ARG_EPS_S_MUE");
    glbSetOscParamByName(p,-glbGetOscParamByName(p,"ARG_EPS_D_MUMU"),  "ARG_EPS_S_MUMU");
    glbSetOscParamByName(p,-glbGetOscParamByName(p,"ARG_EPS_D_TAUMU"), "ARG_EPS_S_MUTAU");
                                                              
    glbSetOscParamByName(p, glbGetOscParamByName(p,"ABS_EPS_D_ETAU"),  "ABS_EPS_S_TAUE");
    glbSetOscParamByName(p, glbGetOscParamByName(p,"ABS_EPS_D_MUTAU"), "ABS_EPS_S_TAUMU");
    glbSetOscParamByName(p, glbGetOscParamByName(p,"ABS_EPS_D_TAUTAU"),"ABS_EPS_S_TAUTAU");
    glbSetOscParamByName(p,-glbGetOscParamByName(p,"ARG_EPS_D_ETAU"),  "ARG_EPS_S_TAUE");
    glbSetOscParamByName(p,-glbGetOscParamByName(p,"ARG_EPS_D_MUTAU"), "ARG_EPS_S_TAUMU");
    glbSetOscParamByName(p,-glbGetOscParamByName(p,"ARG_EPS_D_TAUTAU"),"ARG_EPS_S_TAUTAU");
  }
  else
  {
    glbSetOscParamByName(p, glbGetOscParamByName(p,"ABS_EPS_S_EE"),    "ABS_EPS_D_EE");
    glbSetOscParamByName(p, glbGetOscParamByName(p,"ABS_EPS_S_EMU"),   "ABS_EPS_D_MUE");
    glbSetOscParamByName(p, glbGetOscParamByName(p,"ABS_EPS_S_ETAU"),  "ABS_EPS_D_TAUE");
    glbSetOscParamByName(p,-glbGetOscParamByName(p,"ARG_EPS_S_EE"),    "ARG_EPS_D_EE");
    glbSetOscParamByName(p,-glbGetOscParamByName(p,"ARG_EPS_S_EMU"),   "ARG_EPS_D_MUE");
    glbSetOscParamByName(p,-glbGetOscParamByName(p,"ARG_EPS_S_ETAU"),  "ARG_EPS_D_TAUE");

    glbSetOscParamByName(p, glbGetOscParamByName(p,"ABS_EPS_S_MUE"),   "ABS_EPS_D_EMU");
    glbSetOscParamByName(p, glbGetOscParamByName(p,"ABS_EPS_S_MUMU"),  "ABS_EPS_D_MUMU");
    glbSetOscParamByName(p, glbGetOscParamByName(p,"ABS_EPS_S_MUTAU"), "ABS_EPS_D_TAUMU");
    glbSetOscParamByName(p,-glbGetOscParamByName(p,"ARG_EPS_S_MUE"),   "ARG_EPS_D_EMU");
    glbSetOscParamByName(p,-glbGetOscParamByName(p,"ARG_EPS_S_MUMU"),  "ARG_EPS_D_MUMU");
    glbSetOscParamByName(p,-glbGetOscParamByName(p,"ARG_EPS_S_MUTAU"), "ARG_EPS_D_TAUMU");

    glbSetOscParamByName(p, glbGetOscParamByName(p,"ABS_EPS_S_TAUE"),  "ABS_EPS_D_ETAU");
    glbSetOscParamByName(p, glbGetOscParamByName(p,"ABS_EPS_S_TAUMU"), "ABS_EPS_D_MUTAU");
    glbSetOscParamByName(p, glbGetOscParamByName(p,"ABS_EPS_S_TAUTAU"),"ABS_EPS_D_TAUTAU");
    glbSetOscParamByName(p,-glbGetOscParamByName(p,"ARG_EPS_S_TAUE"),  "ARG_EPS_D_ETAU");
    glbSetOscParamByName(p,-glbGetOscParamByName(p,"ARG_EPS_S_TAUMU"), "ARG_EPS_D_MUTAU");
    glbSetOscParamByName(p,-glbGetOscParamByName(p,"ARG_EPS_S_TAUTAU"),"ARG_EPS_D_TAUTAU");
  }

  return 0;
}


/***************************************************************************
 *                  D E B U G G I N G   F U N C T I O N S                  *
 ***************************************************************************/

static void my_channel_printf_mathematica(FILE *stream, const double *energy,
                                          const double **res, size_t l,size_t c)
{
  size_t i,k;
  glbPrintDelimiter(stream,'l');
  for(k=0;k<c;k++)
  {
    glbPrintDelimiter(stream,'l');
    for(i=0;i<l;i++)
    {
      glbPrintDelimiter(stream,'l');
      fprintf(stream,"%f",energy[i]);
      glbPrintDelimiter(stream,'m');
      fprintf(stream,"%f",res[k][i]);
      glbPrintDelimiter(stream,'r');
      if(i<l-1)  glbPrintDelimiter(stream,'m');
      fprintf(stream,"\n");
    }
    glbPrintDelimiter(stream,'r');
    if(k<c-1)  glbPrintDelimiter(stream,'m');
    fprintf(stream,"\n");
  }
  glbPrintDelimiter(stream,'r');
}


/*************************************************************************** 
 * Print event rates                                                       *
 ***************************************************************************/
int print_rates()
{
  FILE *stream = stdout;
  int exp, rule;

  // Compute true rates 
  glbSetOscillationParameters(true_values);
  glbSetRates();

  // Setup output in Mathematica format 
  glbSetPrintDelimiters("{", ", ", "}");
  glbSetChannelPrintFunction((void *)my_channel_printf_mathematica);

  // Print lists 
  glbPrintDelimiter(stream,'l');
  fprintf(stream, " (* Begin list of experiments *)\n");
  for (exp=0; exp < glb_num_of_exps; exp++)
  {
    glbPrintDelimiter(stream,'l');
    fprintf(stream, " (* Begin experiment %d -- %s *)\n", exp, glbGetFilenameOfExperiment(exp));
    for (rule=0; rule < glbGetNumberOfRules(exp); rule++)
    {
      glbPrintDelimiter(stream,'l');
      fprintf(stream, " (* Rule %d -- signal *) \n", rule);
      glbShowRuleRates(stream, exp, rule, GLB_ALL, GLB_W_EFF, GLB_W_BG, GLB_W_COEFF, GLB_SIG);
      glbPrintDelimiter(stream,'m');
      fprintf(stream,"\n");
      fprintf(stream, " (* Rule %d -- background *) \n", rule);
      glbShowRuleRates(stream, exp, rule, GLB_ALL, GLB_W_EFF, GLB_W_BG, GLB_W_COEFF, GLB_BG);
      glbPrintDelimiter(stream,'r');
      if (rule < glbGetNumberOfRules(exp) - 1)
        glbPrintDelimiter(stream,'m');
      fprintf(stream,"\n");
    }
    glbPrintDelimiter(stream,'r');

    if (exp < glb_num_of_exps-1)
      glbPrintDelimiter(stream,'m');
    fprintf(stream,"\n");
  }
  glbPrintDelimiter(stream,'r');
  fprintf(stream, "\n");

  return 0;
}


// -------------------------------------------------------------------------
int my_print_params(glb_params p)
// -------------------------------------------------------------------------
// Print oscillation parameter vector
// -------------------------------------------------------------------------
{
  int k = 0;
  for (int i=0; i < glbGetNumOfOscParams(); i++)
  {
    double x = glbGetOscParams(p, i);
    if (x != 0.0)
    {
      printf("%s %g   ", glbGetParamName(i), x);
      k++;
      if (k % 4 == 0)
        printf("\n");
    }
  }
  printf("\n");

  return 0;
}


// -------------------------------------------------------------------------
//                    P H Y S I C S   F U N C T I O N S
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
int param_scan(const char *key_string, int n_p, char *params[], double p_min[], double p_max[],
       int p_steps[], unsigned long p_flags[], int n_min_params, char *min_params[],
       int prescan_n_p, char *prescan_params[], double prescan_p_min[],
       double prescan_p_max[], int prescan_p_steps[], unsigned long prescan_p_flags[])
// -------------------------------------------------------------------------
// Perform a scan over an arbitrary subset of the parameter scan,
// in addition marginalizing over arbitrary parameters
// -------------------------------------------------------------------------
// Parameters:
//   key_string: Prefix for output
//   n_p: Number of parameters to scan
//   params: Names of parameters to scan
//   p_min, p_max, p_steps: Min/max values and number of steps for the scan
//     over each parameter (# of scan point = # of steps + 1)
//   p_flags: Extra options for the scan over each parameter
//     (e.g. DEG_LOGSCALE)
//   n_min_params: Number of parameters to marginalize over
//   min_params: Names of parameters to marginalize over
//   prescan_n_p: Number of parameters in degfinder prescan
//   prescan_params: Names of parameters in degfinder prescan
//   prescan_p_min, prescan_p_max, prescan_p_steps: Scan ranges for prescan
//   prescan_p_flags: Extra flags for params in prescan (DEG_LOGSCALE,
//      DEG_S22, DEG_PM)
// -------------------------------------------------------------------------
{
  // Check for invalid oscillation parameter names
  for (int i=0; i < n_p; i++)
  {
    if (glbFindParamByName(params[i]) < 0)
    {
      fprintf(stderr, "Invalid oscillation parameter: %s.\n", params[i]);
      return -1;
    }
    else if (glbFindParamByName(params[i]) == GLB_DM_21) // Set range for solar code
    {
      solar_dm21_min = MAX(1.e-6, MIN(p_min[i], solar_dm21_min));
      solar_dm21_max = MAX(p_max[i], solar_dm21_max);
    }
  }

  for (int i=0; i < prescan_n_p; i++)
  {
    if (glbFindParamByName(prescan_params[i]) < 0)
    {
      fprintf(stderr, "Invalid oscillation parameter: %s.\n", prescan_params[i]);
      return -2;
    }
    else if (glbFindParamByName(prescan_params[i]) == GLB_DM_21) // Set range for solar code
    {
      solar_dm21_min = MAX(1.e-6, MIN(prescan_p_min[i], solar_dm21_min));
      solar_dm21_max = MAX(prescan_p_max[i], solar_dm21_max);
    }
  }

  for (int i=0; i < n_min_params; i++)
  {
    if (glbFindParamByName(min_params[i]) < 0)
    {
      fprintf(stderr, "Invalid oscillation parameter: %s.\n", min_params[i]);
      return -3;
    }
    else if (glbFindParamByName(min_params[i]) == GLB_DM_21) // Set range for solar code
    {
      solar_dm21_min = 1.e-6;
      solar_dm21_max = 1.e-3;
    }
  }

  if (debug_level > 0)
  {
    printf("# \\Delta m_{21}^2 range (if used): %g eV^2 -- %g eV^2\n", solar_dm21_min, solar_dm21_max);
    printf("#\n");
  }

  // Print header
  printf("#       ");
  for (int i=0; i < n_p; i++)
  {
    printf("%15s ", params[i]);
  }
  printf("      chi2_NH          chi2_IH\n");

  // Prepare for prescan in degfinder 
  int n_deg = MAX_DEG;
  double deg_chi2[MAX_DEG];
  glb_params deg_pos[MAX_DEG];
  for (int i=0; i < MAX_DEG; i++)
    deg_pos[i] = glbAllocParams();
  int prescan_param_indices[prescan_n_p];
  for (int i=0; i < prescan_n_p; i++)
    prescan_param_indices[i] = glbFindParamByName(prescan_params[i]);

  // Define projection for pre-scan 
  for (int i=0; i < glbGetNumOfOscParams(); i++)
    glbSetProjectionFlag(prescan_proj, GLB_FIXED, i);
  glbSetDensityProjectionFlag(prescan_proj, GLB_FIXED, GLB_ALL);

  // Define projection for final fit 
  for (int i=0; i < glbGetNumOfOscParams(); i++)
    glbSetProjectionFlag(proj, GLB_FIXED, i);
  for (int i=0; i < n_min_params; i++)
    glbSetProjectionFlagByName(proj, GLB_FREE, min_params[i]);
  glbSetDensityProjectionFlag(proj, GLB_FIXED, GLB_ALL);
//  DecideOnDensityProjection(proj);

  // Compute true rates 
  glbSetOscillationParameters(true_values);
  glbCopyParams(true_values, central_values);
  glbSetCentralValues(central_values);
  glbSetInputErrors(input_errors);
  glbSetRates();


  // Preparations for scan 
  unsigned long n_points = 1;
  for (int i=0; i < n_p; i++)
    n_points *= p_steps[i] + 1;

  // The main loop 
  MPIFOR(j, 0, (int) n_points-1)
  {
    // Compute parameter values at current grid point 
    glbCopyParams(true_values, test_values);
    double p_test_values[n_p];
    for (int i=0; i < n_p; i++)
    {
      // Convert 1d index to a multi-dimensional index for the i-th dimension
      int k = n_points;
      int m = j;
      for (int n=0; n <= i; n++)
      {
        m %= k;
        k /= p_steps[n] + 1;
      }
      m /= k;

      if (p_flags[i] & DEG_LOGSCALE)
      {
        if (p_steps[i] == 0)
          p_test_values[i] = POW10(p_min[i]);
        else
          p_test_values[i] = POW10(p_min[i] + m * (p_max[i]-p_min[i])/p_steps[i]);
      }
      else
      {
        if (p_steps[i] == 0)
          p_test_values[i] = p_min[i];
        else
          p_test_values[i] = p_min[i] + m * (p_max[i]-p_min[i])/p_steps[i];
      }
      glbSetOscParamByName(test_values, p_test_values[i], params[i]);
    }

    // For 5-neutrino scenarios, the distinction between 3+2 and 1+3+1 is hardcoded
    // for compatibility with Thomas' code FIXME
    if (n_flavors >= 5)
    {
#ifndef Ip3pI
      glbSetOscParamByName(test_values,  fabs(glbGetOscParamByName(test_values, "DM41")), "DM41");
      glbSetOscParamByName(test_values,  fabs(glbGetOscParamByName(test_values, "DM51")), "DM51");
#else
      glbSetOscParamByName(test_values, -fabs(glbGetOscParamByName(test_values, "DM41")), "DM41");
      glbSetOscParamByName(test_values,  fabs(glbGetOscParamByName(test_values, "DM51")), "DM51");
#endif
    }

    // Test if GLoBES will accept the chosen oscillation parameters (it may not,
    // for instance if the chosen values of s22thmue and Um4 are inconsistent),
    // and if so, run degfinder
    if (glbSetOscillationParameters(test_values) == 0)
    {
      n_deg = MAX_DEG;
      degfinder(test_values, prescan_n_p, prescan_param_indices, prescan_p_min, prescan_p_max,
                prescan_p_steps, prescan_proj, proj, &n_deg, deg_pos, deg_chi2,
                default_degfinder_flags, prescan_p_flags);
    }
    else
      n_deg = 0;

    // Determine best fit 
    double chi2_NH = 1.0e50;
    double chi2_IH = 1.0e50;
    glb_params Fit_NH = NULL;
    glb_params Fit_IH = NULL;
    for (int i=0; i < n_deg; i++)
    {
      if (glbGetOscParams(deg_pos[i], GLB_DM_31) >= 0.0  &&  deg_chi2[i] < chi2_NH)
      {
        chi2_NH = deg_chi2[i];
        Fit_NH  = deg_pos[i];
      }
      else if (glbGetOscParams(deg_pos[i], GLB_DM_31) < 0.0  &&  deg_chi2[i] < chi2_IH)
      {
        chi2_IH = deg_chi2[i];
        Fit_IH  = deg_pos[i];
      }
    }

    printf("%-7s ", key_string);
    for (int i=0; i < n_p; i++)
    {
      if (p_flags[i] & DEG_LOGSCALE)
        printf("%15.10g ", log10(p_test_values[i]));
      else
        printf("%15.10g ", p_test_values[i]);
    }
    printf("%15.10g %15.10g\n", chi2_NH, chi2_IH);

    // Print parameters at best fit points
    if (debug_level > 0)
    {
      if (Fit_NH)
      {
        printf("PARAMS_NH ");
        my_print_params(Fit_NH);
      }
      if (Fit_IH)
      {
        printf("PARAMS_IH ");
        my_print_params(Fit_IH);
      }
    }
  }

  for (int i=0; i < MAX_DEG; i++)
    if (deg_pos[i] != NULL)
      glbFreeParams(deg_pos[i]);

  return 0;
}


