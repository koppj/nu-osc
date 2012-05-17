/***************************************************************************
 * Interface to other codes (e.g. Thomas Schwetz's)                        *
 ***************************************************************************
 * Author: Joachim Kopp                                                    *
 ***************************************************************************/
//#include <stdlib.h>
//#include <stdio.h>
//#include <string.h>
//#include <math.h>
using namespace std;

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <globes/globes.h>   /* GLoBES library */
#include "nu.h"
#include "const.h"
//#include "external/definitions.h"
#include "sbl/def-reactors.h"
#include "atm/LibWrap/out.interface.hh"

extern gsl_matrix_complex *U;
extern int n_flavors;

// Provide some global variables for Thomas' code
extern Fit fit;
int old_new_main = NEW; // Use OLD or NEW reactor neutrino fluxes? FIXME


/***************************************************************************
 * Initialize external experiments                                         *
 ***************************************************************************/
int ext_init(int ext_flags)
{
  if (ext_flags & EXT_MB)
    initMiniboone();                 // MiniBooNE (\nu)

  if (ext_flags & EXT_MBANTI)
    initMBanti();                    // MiniBooNE (\bar\nu)

  if (ext_flags & EXT_KARMEN  ||  ext_flags & EXT_LSND)
  {
    initFlux();                      // LSND/KARMEN
    initLSND();
    calcKarmen();
  }

  if (ext_flags & EXT_NOMAD)
    initNomad();                     // Nomad

  if (ext_flags & EXT_SBL)
  {
    init_fluxes();                   // Reactor flux
    chooz_init(old_new_main);        // CHOOZ
    init_sbl_reactors(old_new_main); // Bugey 4, Rovno, Krasnoyarsk, ILL, Goesgen
    PV_init();                       // Palo Verde
#ifndef BUGEY_TOTAL_RATE
    bugey_init(old_new_main);        // Bugey
#endif

    fit.invert_S();
    fit.pull_status[FLUX_NORM] = FIXED;
  }

  if (ext_flags & EXT_CDHS)
    initCDHS();

  if (ext_flags & EXT_ATM_COMP)
    atm_init(0x01);

  return 0;
}


/***************************************************************************
 * User-defined prior function---supposed to do the same as the default    *
 * one, but add contributions from external simulations                    *
 ***************************************************************************/
double my_prior(const glb_params in, void* user_data)
{
  glb_params params = glbAllocParams();
  glb_params central_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_projection p = glbAllocProjection();
  glbGetCentralValues(central_values);
  glbGetInputErrors(input_errors);
  glbGetProjection(p);
  glbGetOscillationParameters(params);
    // NOTE: We don't use the "in" parameters because they may contain inconsistent
    // values for Ue3 <-> th13, Ue4 <-> th14, etc. On the other, since GLoBES calls
    // the prior function immediately after a rate calculation, we can be sure that
    // the oscillation engine still has the correct parameter values stored

  int i;
  int ext_flags = *((int *) user_data);
  double pv = 0.0;
  double fitvalue,centralvalue,inputerror;

  // Force the active-sterile mixing angles to be < \pi/4 to prevent the minimizer
  // from jumping between the 3+2 and 1+3+1 cases FIXME
  if (n_flavors >= 4  &&
      (fabs(glbGetOscParamByName(params, "TH14")) > M_PI/4  ||
       fabs(glbGetOscParamByName(params, "TH24")) > M_PI/4  ||
       fabs(glbGetOscParamByName(params, "TH34")) > M_PI/4))
  {
    pv += 1e11;
    goto my_prior_end;
  }
  if (n_flavors >= 5  &&
      (fabs(glbGetOscParamByName(params, "TH15")) > M_PI/4  ||
       fabs(glbGetOscParamByName(params, "TH25")) > M_PI/4  ||
       fabs(glbGetOscParamByName(params, "TH35")) > M_PI/4))
  {
    pv += 2e11;
    goto my_prior_end;
  }

  // For 5-neutrino scenarios, the distinction between 3+2 and 1+3+1 is hardcoded
  // for compatibility with Thomas' code FIXME
  if (n_flavors >= 5)
  {
#ifndef Ip3pI
    if (glbGetOscParamByName(params, "DM41") < 0 ||
        glbGetOscParamByName(params, "DM51") < 0)
    {
      pv += 3e11;
      goto my_prior_end;
    }
#else
    if (glbGetOscParamByName(params, "DM41") > 0 ||
        glbGetOscParamByName(params, "DM51") < 0)
    {
      pv += 4e11;
      goto my_prior_end;
    }
#endif
  }

  // Add chi^2 from external codes
  if (ext_flags)
  {
    // Workaround for Thomas' code requiring 0 < dm41sq, dm51sq <~ 10^2
    // (actually slightly less than 2 due to a bug)
    if (n_flavors >= 4)
    {
      double dm41 = fabs(glbGetOscParamByName(params, "DM41"));
      if (dm41 < 0 || dm41 > 70)
      {
        pv += 5e11;
        goto my_prior_end;
      }
    }
    if (n_flavors >= 5)
    {
      double dm51 = fabs(glbGetOscParamByName(params, "DM51"));
      if (dm51 < 0 || dm51 > 70)
      {
        pv += 6e11;
        goto my_prior_end;
      }
    }

    struct params sbl_params;    // Parameter data structure for Thomas' code
    gsl_complex (*_U)[n_flavors] =
      (gsl_complex (*)[n_flavors]) gsl_matrix_complex_ptr(snu_get_U(), 0, 0);
    sbl_params.Ue3     = gsl_complex_abs(_U[0][2]);
    sbl_params.Ue[I4]  = 0.0;
    sbl_params.Um[I4]  = 0.0;
    sbl_params.dmq[I4] = 0.0;
    sbl_params.Ue[I5]  = 0.0;
    sbl_params.Um[I5]  = 0.0;
    sbl_params.dmq[I5] = 0.0;
    sbl_params.delta   = 0.0;
    if (n_flavors >= 4)
    {
      sbl_params.Ue[I4]  = gsl_complex_abs(_U[0][3]);
      sbl_params.Um[I4]  = gsl_complex_abs(_U[1][3]);
      sbl_params.Ue3     = gsl_complex_abs(_U[0][2]);
      sbl_params.dmq[I4] = fabs(glbGetOscParamByName(params, "DM41"));
    }
    if (n_flavors >= 5)
    {
      sbl_params.Ue[I5]  = gsl_complex_abs(_U[0][4]);
      sbl_params.Um[I5]  = gsl_complex_abs(_U[1][4]);
      sbl_params.dmq[I5] = fabs(glbGetOscParamByName(params, "DM51"));
      sbl_params.delta   = gsl_complex_arg(
          gsl_complex_mul(gsl_complex_mul(_U[1][3], gsl_complex_conjugate(_U[0][3])),
                          gsl_complex_mul(_U[0][4], gsl_complex_conjugate(_U[1][3]))) );
    }

    if (ext_flags & EXT_MB)
      pv += chi2mb475(sbl_params);
    if (ext_flags & EXT_MBANTI)
      pv += chi2_MBA_475(sbl_params);
    if (ext_flags & EXT_KARMEN)
      pv += chi2karmen(sbl_params);
    if (ext_flags & EXT_LSND)
      pv += chi2lsnd(sbl_params);
    if (ext_flags & EXT_SBL)
      pv += chi2reactor(sbl_params);
    if (ext_flags & EXT_NOMAD)
      pv += chi2nomad(sbl_params);
    if (ext_flags & EXT_CDHS)
      pv += chi2cdhs(sbl_params);
    if (ext_flags & EXT_ATM_TABLE)
      pv += chi2atm(sbl_params);
    if (ext_flags & EXT_ATM_COMP)
    {
      // Workaround for Michele's code requiring dm31sq < 1e-2
      if (fabs(glbGetOscParamByName(params, "DM31")) > 0.99e-2)
        pv += 1e15;
      else if (glbGetNumOfOscParams() == 51+1)  // 3 flavors
        pv += atm_chisq(glbGetOscParamByName(params, "TH23"),
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        glbGetOscParamByName(params, "DM31"));
      else if (glbGetNumOfOscParams() == 92+3) // 4 flavors
        pv += atm_chisq(glbGetOscParamByName(params, "TH23"),
                        glbGetOscParamByName(params, "TH24"),
                        0.0,
                        glbGetOscParamByName(params, "TH34"),
                        0.0,
                        0.0, //FIXME Include phases
                        0.0,
                        glbGetOscParamByName(params, "DM31"));
      else if (glbGetNumOfOscParams() == 145+5) // 5 flavors
        pv += atm_chisq(glbGetOscParamByName(params, "TH23"),
                        glbGetOscParamByName(params, "TH24"),
                        glbGetOscParamByName(params, "TH25"),
                        glbGetOscParamByName(params, "TH34"),
                        glbGetOscParamByName(params, "TH35"),
                        0.0,
                        0.0,
                        glbGetOscParamByName(params, "DM31"));
      else
        pv -= 1.e30;
    } // ext_flags & EXT_ATM_COMP
  } // if (ext_flags)

  // Add oscillation parameter priors
  for(i=0; i < glbGetNumOfOscParams(); i++)
    if(glbGetProjectionFlag(p,i)==GLB_FREE)
    {
      fitvalue     = glbGetOscParams(params,i);
      centralvalue = glbGetOscParams(central_values,i);
      inputerror   = glbGetOscParams(input_errors,i);
      if (inputerror > 1e-12)
        pv += SQR((centralvalue-fitvalue)/inputerror);
    }

  /* Add matter parameter priors */
  for(i=0; i < glb_num_of_exps; i++)
  if(glbGetDensityProjectionFlag(p,i) == GLB_FREE)
  {
    fitvalue     = glbGetDensityParams(params,i);
    centralvalue = 1.0;
    inputerror   = glbGetDensityParams(input_errors,i);
    if(inputerror > 1e-12)
      pv += SQR((centralvalue-fitvalue)/inputerror);
  }

my_prior_end:
  glbFreeParams(central_values);
  glbFreeParams(input_errors);
  glbFreeParams(params);
  glbFreeProjection(p);
  return pv;
}

