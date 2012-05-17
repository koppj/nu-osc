// -------------------------------------------------------------------------
// Some routines to allow GLoBES to be called from Thomas' fitting code
// -------------------------------------------------------------------------
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

#include "sbl/definitions.h"
using namespace std;

extern int default_degfinder_flags;

extern glb_params true_values;
extern glb_params test_values;
extern glb_params central_values;
extern glb_params input_errors;
extern glb_projection proj;
extern glb_projection prescan_proj;

#define MAX_DEG    100       // Maximum number of degeneracies to expect
#define MAX_PARAMS  20



// -------------------------------------------------------------------------
int initGlb(int n_flavors)
// -------------------------------------------------------------------------
// Initialize GLoBES to allow for MINOS simulation
// -------------------------------------------------------------------------
{
  double true_theta12 = asin(sqrt(0.32));
  double true_theta13 = 0.0;
  double true_theta23 = M_PI/4;
  double true_deltacp = 3.0*M_PI/2.0;
  double true_sdm = 7.6e-5;
  double true_ldm = 2.4e-3;

  double prior_th12    = 0.05 * true_theta12;
  double prior_th13    = 0.0;
  double prior_th23    = 0.0;
  double prior_deltacp = 0.0;
  double prior_sdm     = 0.05 * true_sdm;
  double prior_ldm     = 0.0;

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


  // Initialize libglobes 
  setenv("GLB_PATH", "glb/nova:glb/t2k:glb/bbeam:glb/dchooz:glb/wbb_wc:"
                     "glb/wbb_lar:glb/nufact:glb/minos-nc:glb/miniboone", 1);
  glbInit("nu-code");
  glbSelectMinimizer(GLB_MIN_POWELL); // Parts of code work ONLY with GLB_MIN_POWELL !!! 

  // Initialize and register non-standard probability engine. This has to be done
  // before any calls to glbAllocParams() or glbAllocProjections()
  if (n_flavors == 4)
  {
    int default_rotation_order[][2] = { {3,4}, {2,4}, {1,4}, {2,3}, {1,3}, {1,2} };
    int default_phase_order[] = { -1,  1, -1, -1,  0,   2};
    snu_init_probability_engine(n_flavors, default_rotation_order, default_phase_order);
    glbRegisterProbabilityEngine(6*SQR(n_flavors)-n_flavors, &snu_probability_matrix,
      &snu_set_oscillation_parameters, &snu_get_oscillation_parameters, NULL);
  }
  else if (n_flavors == 5)
  {
    // FIXME - Where to put the extra phases?
    int default_rotation_order[][2] = { {4, 5}, {3,5}, {2,5}, {1,5}, {3,4}, {2,4}, {1,4},
                                        {2,3}, {1,3}, {1,2} };
    int default_phase_order[] = { -1, -1, -1, -1, -1,  1, -1, -1,  0,   2};
    snu_init_probability_engine(n_flavors, default_rotation_order, default_phase_order);
    glbRegisterProbabilityEngine(6*SQR(n_flavors)-n_flavors, &snu_probability_matrix,
      &snu_set_oscillation_parameters, &snu_get_oscillation_parameters, NULL);
  }
  else
  {
    printf("This code only supports 4- or 5-neutrino scenarios.\n");
    return -1;
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
  glbDefineChiFunction(&chiMINOS,         5, "chiMINOS",         NULL);
  chiMB_init();
  glbDefineChiFunction(&chiMBanti_nu2010, 0, "chiMBanti_nu2010", NULL);

  // Load experiments
  int n_exps = 1;
  char *exps[] = { "MINOS_NC" };
  if (load_exps(n_exps, exps) < 0)  // Load experiments 

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

  // Copy true params to other parameter vectors AFTER NSI params have been
  // read from the environment variables
  glbCopyParams(true_values, test_values);
  glbCopyParams(true_values, central_values);

  return 0;
}


// -------------------------------------------------------------------------
double chi2MINOS(params p, int n_flavors, int n_min_params, char **min_params,
                 glb_params bf_params=NULL)
// -------------------------------------------------------------------------
// Compute chi^2 using Thomas' parameterization of active--sterile neutrino
// mixing
// -------------------------------------------------------------------------
{
  // Parameters of prescan in degfinder
  int n_prescan_params           = 0;
  int prescan_params[MAX_PARAMS] = {};
  double prescan_min[MAX_PARAMS] = {};
  double prescan_max[MAX_PARAMS] = {};
  int prescan_steps[MAX_PARAMS]  = {};
  
//  int n_prescan_params           = 3; // For 4-flavor BF
//  int prescan_params[MAX_PARAMS] = { glbFindParamByName("TH24"),
//                                     glbFindParamByName("TH34"),
//                                     glbFindParamByName("TH13") };
//  double prescan_min[MAX_PARAMS] = { -0.7,   -0.7,  -2  };
//  double prescan_max[MAX_PARAMS] = {  0.7,    0.7,  -0.2  };
//  int prescan_steps[MAX_PARAMS]  = { 14,     14,    10    };

//  int n_prescan_params           = 3; // For 5-flavor BF
//  int prescan_params[MAX_PARAMS] = { glbFindParamByName("TH34"),
//                                     glbFindParamByName("TH13"),
//                                     glbFindParamByName("TH35") };
//  double prescan_min[MAX_PARAMS] = { -0.7,   -2,     -0.7};
//  double prescan_max[MAX_PARAMS] = {  0.7,   -0.2,    0.7};
//  int prescan_steps[MAX_PARAMS]  = {  7,     10,      7, };

  int n_deg = MAX_DEG;
  double deg_chi2[MAX_DEG];
  glb_params deg_pos[MAX_DEG];
  for (int i=0; i < MAX_DEG; i++)
    deg_pos[i] = glbAllocParams();

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

  glbSetOscillationParameters(true_values);
  glbCopyParams(true_values, central_values);
  glbSetCentralValues(central_values);
  glbSetInputErrors(input_errors);
  glbSetRates();

  // Convert Thomas' parametrization of U to ours
  glbCopyParams(true_values, test_values);
  if (n_flavors == 4)
  {
    double th14 = asin(p.Ue[I4]);
    double th24 = asin(p.Um[I4]/cos(th14));
    double th13 = asin(p.Ue3/(cos(th14)));
    glbSetOscParamByName(test_values, th14, "TH14");
    glbSetOscParamByName(test_values, th24, "TH24");
    glbSetOscParamByName(test_values, th13, "TH13");
    glbSetOscParamByName(test_values, p.dmq[I4], "DM41");
  }
  else if (n_flavors == 5)
  {
    double th15 = asin(p.Ue[I5]);
    double th25 = asin(p.Um[I5]/cos(th15));
    double th14 = asin(p.Ue[I4]/cos(th15));
    double th24 = asin((p.Um[I4] + sin(th14)*sin(th15)*sin(th25))/(cos(th14)*cos(th25)));
    double th13 = asin(p.Ue3/(cos(th14)*cos(th15)));
    glbSetOscParamByName(test_values, th15, "TH15");
    glbSetOscParamByName(test_values, th25, "TH25");
    glbSetOscParamByName(test_values, th14, "TH14");
    glbSetOscParamByName(test_values, th24, "TH24");
    glbSetOscParamByName(test_values, th13, "TH13");
    glbSetOscParamByName(test_values, p.dmq[I4], "DM41");
    glbSetOscParamByName(test_values, p.dmq[I5], "DM51");

//    printf("%g %g %g %g %g %g %g\n", th15, th25, th14, th24, th13, p.dmq[I4], p.dmq[I5]);
  }

  // Run degfinder
  n_deg = MAX_DEG;
  degfinder(test_values, n_prescan_params, prescan_params, prescan_min, prescan_max,
            prescan_steps, prescan_proj, proj, &n_deg, deg_pos, deg_chi2,
            default_degfinder_flags, NULL);

  // Determine best fit
  double chi2min = 1.0e10;
  for (int i=0; i < n_deg; i++)
  {
    if (deg_chi2[i] < chi2min)
    {
      chi2min = deg_chi2[i];
      if (bf_params)
        glbCopyParams(deg_pos[i], bf_params);
    }
  }

  for (int i=0; i < MAX_DEG; i++)
    if (deg_pos[i] != NULL)
      glbFreeParams(deg_pos[i]);

  return chi2min;
}


// -------------------------------------------------------------------------
int checkBF(int n_flavors)
// -------------------------------------------------------------------------
// Compute chi^2 at Thomas' best fit point
// -------------------------------------------------------------------------
{
  //FIXME Don't forget to switch on/off prescan if necessary!
  char *min_params3[]     = { "TH23", "DM31", "TH13" };
  char *min_params4[]     = { "TH23", "DM31", "TH34", "TH24", "TH13" };
  char *min_params4_all[] = { "TH23", "DM31", "TH34", "TH14", "TH24", "TH13", "DM41" };
  char *min_params5[]     = { "TH23", "DM31", "TH13", "TH34", "TH35" };
  char *min_params5_all[] = { "TH23", "DM31", "TH13", "TH14", "TH24", "TH15", "TH25",
                              "TH34", "TH35", "DM41", "DM51" };
  params p;
  double chi2=0.0, chi2all=0.0;
  glb_params bf_all    = glbAllocParams();
  glb_params bf_thomas = glbAllocParams();

  memset(&p, 0, sizeof(p));
  printf("# Best fit point for standard oscillations:\n");
  printf("#   chi^2 = %15.10g\n", chi2MINOS(p, n_flavors, 3, min_params3));
  printf("#\n");

  if (n_flavors == 4)
  {
    p.Ue[I4]  = 0.117;
    p.Um[I4]  = 0.0;
    p.Ue[I5]  = 0.0;
    p.Um[I5]  = 0.0;
    p.Ue3     = 0.0;
    p.delta   = 0.0;
    p.dmq[I4] = 0.45;
    p.dmq[I5] = 0.0;

    chi2 = chi2MINOS(p, n_flavors, 5, min_params4, bf_thomas);
//    memset(&p, 0, sizeof(p));
    chi2all = chi2MINOS(p, n_flavors, 7, min_params4_all, bf_all);
  }
//  else if (n_flavors == 5) // Table II, draft version, 10 Mar 2011
//  {
//    p.Ue[I4]  = 0.131;
//    p.Um[I4]  = 0.170;
//    p.Ue[I5]  = 0.115;
//    p.Um[I5]  = 0.142;
//    p.Ue3     = 0.0;
//    p.delta   = 1.62;
//    p.dmq[I4] = 0.47;
//    p.dmq[I5] = 0.93;
//
//    chi2 = chi2MINOS(p, n_flavors, 5, min_params5, bf_thomas);
////    memset(&p, 0, sizeof(p));
//    chi2all = chi2MINOS(p, n_flavors, 11, min_params5_all, bf_all);
//  }
  else if (n_flavors == 5) // Table II, final paper 1103.4570 (3+2 case)
  {
    p.Ue[I4]  = 0.128;
    p.Um[I4]  = 0.165;
    p.Ue[I5]  = 0.138;
    p.Um[I5]  = 0.148;
    p.Ue3     = 0.0;
    p.delta   = 1.64;
    p.dmq[I4] = 0.47;
    p.dmq[I5] = 0.87;

    chi2 = chi2MINOS(p, n_flavors, 5, min_params5, bf_thomas);
//    memset(&p, 0, sizeof(p));
    chi2all = chi2MINOS(p, n_flavors, 11, min_params5_all, bf_all);
  }
  else
  {
    glbFreeParams(bf_thomas);
    glbFreeParams(bf_all);
    return -1;
  }

  printf("# Thomas' best fit point (%d flavors):\n", n_flavors);
  printf("#   Ue4   = %g\n", p.Ue[I4]);
  printf("#   Um4   = %g\n", p.Um[I4]);
  printf("#   Ue5   = %g\n", p.Ue[I5]);
  printf("#   Um5   = %g\n", p.Um[I5]);
  printf("#   delta = %g\n", p.delta);
  printf("#   dm41  = %g\n", p.dmq[I4]);
  printf("#   dm51  = %g\n", p.dmq[I5]);
  printf("#   chi^2 = %15.10g\n", chi2);
  glbPrintParams(stdout, bf_thomas);
  printf("#\n");

  printf("# MINOS overall best fit:\n");
  printf("#   chi^2 = %15.10g\n", chi2all);
  printf("# Best fit parameters:\n");
  glbPrintParams(stdout, bf_all);
  printf("#\n");

  glbFreeParams(bf_thomas);
  glbFreeParams(bf_all);
  return 0;
}


