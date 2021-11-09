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

#include <complex>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <globes/globes.h>   /* GLoBES library */
#include "LSND.h"
#include "KARMEN.h"
#include "constraints.h"
#include "sbl/definitions.h"
#include "reactors/definitions.h"
#include "db-neos/db-neos.h"
#include "atm/LibWrap/out.interface.hh"
#include "solar/LibWrap/out.interface.hh"
#include "nu.h"
#include "const.h"


// Global variables
extern gsl_matrix_complex *U;
extern int n_flavors;
extern int action;
extern char mb_tune[];

// Minimum/maximum values of dm21^2 for the solar neutrino code
extern double true_sdm;
double solar_dm21_min = NAN;
double solar_dm21_max = NAN;

// Provide some global variables for Thomas' code
namespace ns_reactor
{
  extern Fit fit;
  extern Rate_coef rate;
}

// Parameters for Michele's atmospheric neutrino analysis
extern int atm_decouple_e; // Whether to decouple electron neutrinos by hand
                           // and instead include sterile neutrino matter effects

// Require mixing angles > 0?
extern int theta_positive;

// ATTENTION: Since the 2012 update of Thomas' reactor code, I'm not
// sure if this flag is still used consistently. Therefore, DO NOT
// change this!!! Fits with the old fluxes are no longer possible.
int old_new_main = NEW; // Use OLD or NEW reactor neutrino fluxes?


// Daya Bay and NEOS objects for Alvaro's code
DB_class DB_fitter;
Neos_DB_class Neos_fitter;

// LSND object for Ivan's code
#ifdef NU_USE_NUSQUIDS
ns_sbl::LSND        *LSND_ivan_fitter   = NULL;
ns_sbl::KARMEN      *KARMEN_ivan_fitter = NULL;
ns_sbl::Constraints *ext_constraints    = NULL;
#endif


/***************************************************************************
 * Two helper functions for solar_parameters, adapted from Michele's       *
 * test-rnfchisq.cc)                                                       *
 ***************************************************************************/
int rephase_R(const int f, const double phase, gsl_matrix_complex *U)
{
  if (!U)
    return -1;

  const int n_flavors = U->size1;
  if(f < 0 || f >= n_flavors)
  {
    fprintf(stderr, "rephase_R: Invalid flavor index");
    return -2;
  }

  complex<double> (*_U)[n_flavors]
    = (complex<double> (*)[n_flavors]) gsl_matrix_complex_ptr(U, 0, 0);
  complex<double> pol = polar(1., phase);
  for (int i=0; i < n_flavors; i++)
    _U[i][f] *= pol;

  return 0;
}


int rotate_R(const int f1, const int f2, const double angle, const double phase,
             gsl_matrix_complex *U)
{
  if (!U)
    return -1;

  const int n_flavors = U->size1;
  if(f1 < 0 || f1 >= n_flavors || f2 < 0 || f2 >= n_flavors || f1 == f2)
  {
    fprintf(stderr, "rotate_R: Invalid flavor index");
    return -2;
  }

  complex<double> (*_U)[n_flavors]
    = (complex<double> (*)[n_flavors]) gsl_matrix_complex_ptr(U, 0, 0);
  double cs = cos(angle);
  complex<double> sp = polar(sin(angle), phase);
  complex<double> sm = -conj(sp);
  for (int i=0; i < n_flavors; i++)
  {
    complex<double> t = sp * _U[i][f1] + cs * _U[i][f2];
    _U[i][f1] = cs * _U[i][f1] + sm * _U[i][f2];
    _U[i][f2] = t;
  }

  return 0;
}


/***************************************************************************
 * Convert mixing matrix to Michele's parameterization for solar neutrinos *
 * (function adapted from Michele's test-rnfchisq.cc)                      *
 * The function first removes any Majorana phases (which do not affect     *
 * neutrino oscillations, then separates off the 12-rotation and a phase   *
 * which is chosen such that \sum_j U_{j1} U_{j2}^* is real.               *
 * See notes attached to Michele's email from 2012-05-20 for details       *
 ***************************************************************************/
int solar_parameters(gsl_matrix_complex *U, Angles *a)
{
  if (!U)
    return -1;

  const int n_flavors = U->size1;
  gsl_matrix_complex *T = gsl_matrix_complex_alloc(n_flavors, n_flavors);
  if (!T)
    return -2;
  complex<double> (*_T)[n_flavors]
    = (complex<double> (*)[n_flavors]) gsl_matrix_complex_ptr(T, 0, 0);
  gsl_matrix_complex_memcpy(T, U);

  // Remove Majorana phases
  rephase_R(0, -arg(_T[0][0]), T);
  rephase_R(1, -arg(_T[0][1]), T);
  if (fabs(imag(_T[0][0])) + fabs(imag(_T[0][1])) > 1e-10)
  {
    fprintf(stderr, "solar_parameters: Removal of Majorana phases failed");
    return -3;
  }

  // Extract theta-12
  a->the12 = atan2(real(_T[0][1]), real(_T[0][0]));
  rotate_R(0, 1, -a->the12, 0., T);
  if (abs(_T[0][1]) > 1e-10)
  {
    fprintf(stderr, "solar_parameters: Extraction of th12 angle failed");
    return -4;
  }

  // 12-phase
  complex<double> t = 0.0;
  for (int i=3; i < n_flavors; i++)
    t += _T[i][0] * conj(_T[i][1]);
  a->dlt12 = arg(t);
  rephase_R(0, -a->dlt12, T);

  t = 0.0;
  for (int i=3; i < n_flavors; i++)
    t += _T[i][0] * conj(_T[i][1]);
  if (fabs(imag(t)) > 1e-10)
  {
    fprintf(stderr, "solar_parameters: Extraction of dlt12 phase failed");
    return -5;
  }

  // Remaining parameters
  a->eta_e = norm(_T[0][0]);

  a->ste_D = a->ste_N = 0.0;
  for (int i=3; i < n_flavors; i++)
  {
    a->ste_D += norm(_T[i][1]) - norm(_T[i][0]);
    a->ste_N += 2. * real(_T[i][0] * conj(_T[i][1]));
  }

  a->cst_e = 0.;
  a->cst_a = 1.;
  for (int j=0; j < n_flavors; j++)
  {
    double t = norm(_T[0][j]);
    a->cst_e += SQR(t);
    for (int i=3; i < n_flavors; i++)
      a->cst_a -= t * norm(_T[i][j]);
  }

  gsl_matrix_complex_free(T);

  return 0;
}


/***************************************************************************
 * Initialize external experiments                                         *
 ***************************************************************************/
int ext_init(int ext_flags, int use_feldman_cousins)
{
  int status;

  if (ext_flags & EXT_MB  ||  ext_flags & EXT_MB_300)
  {
    printf("# Initializing MiniBooNE code (neutrinos) ...\n");
    initMiniboone();                 // MiniBooNE (\nu)
  }

  if (ext_flags & EXT_MBANTI  ||  ext_flags & EXT_MBANTI_200)
  {
    printf("# Initializing MiniBooNE code (antineutrinos) ...\n");
    initMBanti();                    // MiniBooNE (\bar\nu)
  }

  if (ext_flags & EXT_KARMEN  ||  ext_flags & EXT_LSND)
  {
    printf("# Initializing KARMEN and LSND codes ...\n");
    initFlux();                      // LSND/KARMEN
    initLSND();
    calcKarmen();
  }

  if (ext_flags & EXT_NOMAD)
  {
    printf("# Initializing NOMAD code ...\n");
    initNomad();                     // Nomad
  }


  if (ext_flags & EXT_REACTORS)
  {
    bool sbl_exps[ns_reactor::N_EXP];  // Which experiments to use in Thomas' \nu_e codes
    for (int i=0; i < ns_reactor::N_EXP; i++)
      sbl_exps[i] = false;
    ns_reactor::init_fluxes();
    ns_reactor::rate.init();

    #ifdef USE_SBL
      printf("# Initializing SBL reactor code ...\n");
      ns_reactor::init_sbl_reactors(old_new_main); // Bugey 4, Rovno, Krasnoyarsk, ILL, Goesgen
      sbl_exps[ns_reactor::SBL] = true;
    #endif
    #ifdef USE_CHOOZ
      printf("# Initializing Chooz code ...\n");
      ns_reactor::chooz_init(old_new_main);        // Chooz
      sbl_exps[ns_reactor::CHOOZ] = true;
    #endif
    #ifdef USE_PV
      printf("# Initializing Palo Verde code ...\n");
      ns_reactor::PV_init();                       // Palo Verde
      sbl_exps[ns_reactor::PV] = true;
    #endif
    #ifdef USE_KAML
      printf("# Initializing KamLAND code ...\n");
      ns_reactor::init_kaml();                     // KamLAND
      sbl_exps[ns_reactor::KAML] = true;
    #endif
    #ifdef USE_DC
      printf("# Initializing Double Chooz code ...\n");
      ns_reactor::dc_init(old_new_main);           // Double Chooz
      sbl_exps[ns_reactor::DC] = true;
    #endif
//    #ifdef USE_DB
//      printf("# Initializing Daya Bay sterile neutrino code ...\n");
//      ns_reactor::DB_init();                       // Daya Bay
//      sbl_exps[ns_reactor::DB] = true;
//    #endif
//    #ifdef USE_DB_3F
//      printf("# Initializing Daya Bay 3-flavor code ...\n");
//      ns_reactor::DB_init();                       // Daya Bay
//      sbl_exps[ns_reactor::DB] = true;
//    #endif
    #ifdef USE_RENO
      printf("# Initializing RENO code ...\n");
      ns_reactor::RENO_init();                     // RENO
      sbl_exps[ns_reactor::RENO] = true;
    #endif
    #ifdef USE_BUGEY_SP
      printf("# Initializing Bugey code ...\n");
      ns_reactor::bugey_init(old_new_main);        // Bugey
      sbl_exps[ns_reactor::BUG_SP] = true;
    #endif
    #ifdef USE_DANSS
      printf("# Initializing DANSS code ...\n");
      ns_reactor::danss_init(old_new_main);        // DANSS
      sbl_exps[ns_reactor::DANSS] = true;
    #endif
//    #ifdef USE_NEOS
//      printf("# Initializing NEOS code ...\n");
//      ns_reactor::neos_init(old_new_main);         // NEOS
//      sbl_exps[ns_reactor::NEOS] = true;
//    #endif
    #ifdef USE_GAL
      printf("# Initializing Gallium code ...\n");
      ns_reactor::gallium_init();                  // Gallium radioactive source experiments
      sbl_exps[ns_reactor::GAL] = true;
    #endif

    #ifdef USE_DB_ALVARO
      printf("# Initializing Alvaro's Daya Bay code ...\n");
      DB_fitter.DB_init();                         // Daya Bay
    #endif
    #ifdef USE_NEOS_ALVARO
      printf("# Initializing Alvaro's NEOS code ...\n");
      Neos_fitter.Neos_DB_init();                  // NEOS
    #endif

    printf("# Setting up pull matrices for SBL codes ...\n");
    ns_reactor::fit.invert_S();
    ns_reactor::fit.set_experiments(sbl_exps);
    ns_reactor::fit.pull_status[ns_reactor::FLUX_NORM]
      = ns_reactor::FIXED; // See Thomas' email from 2012-05-23 for explanation
  }

  if (ext_flags & EXT_CDHS)          // CDHS
  {
    printf("# Initializing CDHS code ...\n");
    initCDHS();
  }

  if (ext_flags & EXT_ATM_COMP  ||  ext_flags & EXT_DEEPCORE) // Michele's atmospherics code
  {
#ifdef USE_ATM
    if (ext_flags & EXT_ATM_COMP  &&  ext_flags & EXT_DEEPCORE)
    {
      printf("# Initializing Michele's DeepCore + SuperK code ...\n");
      atm_init(0x03);
    }
    else if (ext_flags & EXT_ATM_COMP)
    {
      printf("# Initializing Michele's SuperK code ...\n");
      atm_init(0x01);
    }
    else if (ext_flags & EXT_DEEPCORE)
    {
      printf("# Initializing Michele's DeepCore code ...\n");
      atm_init(0x02);
    }
#else
    printf("ext_init: Compiled without atmospherics support.\n");
    exit(-1);
#endif
  }

  if (ext_flags & EXT_SOLAR)         // Michele's solar neutrino code
  {
#ifdef USE_SOLAR
    printf("# Initializing solar neutrino code ...\n");
    solar_dm21_min = true_sdm;
    solar_dm21_max = true_sdm;

    enum {
      EXP_Chlorine  = 1 <<  0,
      EXP_Gallium   = 1 <<  1,
      EXP_SK1_nadir = 1 <<  2,
      EXP_SK1_spect = 1 <<  3,
      EXP_SK1_whole = 1 <<  4,
      EXP_SK2_spect = 1 <<  5,
      EXP_SK3_spect = 1 <<  6,
      EXP_SK4_spect = 1 <<  7,
      EXP_SNO1_pure = 1 <<  8,
      EXP_SNO2_salt = 1 <<  9,
      EXP_SNO3_henc = 1 << 10,
      EXP_SNO_full  = 1 << 11,
      EXP_BX1_lower = 1 << 12,
      EXP_BX1_upper = 1 << 13,
      EXP_BX2_infra = 1 << 14
    };

    const uint solar_exp_mask = EXP_Chlorine | EXP_Gallium | EXP_SK1_nadir |
                                EXP_SK2_spect | EXP_SK3_spect | EXP_SK4_spect |
                                EXP_SNO1_pure | EXP_SNO2_salt | EXP_SNO3_henc |
                                EXP_BX1_lower | EXP_BX1_upper | EXP_BX2_infra;
    sun_init(solar_exp_mask);

    // This was for the fall 2012 version of the solar code
//    enum {
//      EXP_Chlorine = 1 << 0,
//      EXP_Gallium  = 1 << 1,
//      EXP_SuperK   = 1 << 2,
//      EXP_SNO_pure = 1 << 3,
//      EXP_SNO_salt = 1 << 4,
//      EXP_SNO_henc = 1 << 5,
//      EXP_SNO_full = 1 << 6,
//      EXP_BX_lower = 1 << 7,
//      EXP_BX_upper = 1 << 8
//    };
//    const uint solar_exp_mask = EXP_Chlorine | EXP_Gallium  | EXP_SuperK |
//                                EXP_BX_upper | EXP_BX_lower |
//                                EXP_SNO_pure | EXP_SNO_salt | EXP_SNO_henc;
//    sun_init(solar_exp_mask);

#else
    printf("ext_init: Compiled without solar support.\n");
    exit(-1);
#endif
  }

  if (ext_flags & EXT_MINOS2016)      // MINOS 2016
  {
    printf("# Initializing MINOS 2016 code ...\n");
    MINOS_2016_init();
  }

  if (ext_flags & EXT_MINOS2017)      // MINOS 2017
  {
    printf("# Initializing MINOS/MINOS+ 2017 code ...\n");
    MINOS_2017_init();
  }

#ifdef NU_USE_NUSQUIDS
  if (ext_flags & EXT_LSND_IVAN)
  {
    printf("# Initializing Ivan's LSND code ...\n");
    LSND_ivan_fitter = new ns_sbl::LSND();
    if (!LSND_ivan_fitter)
      return -1;
  }
  if (ext_flags & EXT_KARMEN_IVAN)
  {
    printf("# Initializing Ivan's KARMEN code ...\n");
    KARMEN_ivan_fitter = new ns_sbl::KARMEN();
    if (!KARMEN_ivan_fitter)
      return -1;
  }
  if (ext_flags & EXT_DECAY_KINEMATICS  ||  ext_flags & EXT_FREE_STREAMING)
  {
    printf("# Initializing external constraints code (beta decay/cosmo) ...\n");
    ext_constraints = new ns_sbl::Constraints();
    if (!ext_constraints)
      return -1;
  }
#endif

  if (ext_flags & EXT_MB_JK)
  {
    printf("# Initializing Joachim's MiniBooNE code ...\n");
    if ((status=chiMB_jk_init(mb_tune)) != 0)
      return status;
  }
  if (ext_flags & EXT_MUBOONE)
  {
    printf("# Initializing microBooNE code ...\n");
    if ((status=chiMuBooNE_init("", use_feldman_cousins)) != 0)
      return status;
  }
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
    // values for Ue3 <-> th13, Ue4 <-> th14, etc. On the other hand, since GLoBES calls
    // the prior function immediately after a rate calculation, we can be sure that
    // the oscillation engine still has the correct parameter values stored

  int i;
  prior_params *pp = (prior_params *) user_data;
  double pv = 0.0;
  double fitvalue,centralvalue,inputerror;

  // We do, however check that none of the parameters in "in" are NAN because
  // that is snu_set_oscillation_parameter's way of indicating inconsistencies
  for (int i=0; i < glbGetNumOfOscParams(); i++)
    if (isnan(glbGetOscParams(params, i)))
    {
      pv += 1.e25;
      goto my_prior_end;
    }

  // Force the active-sterile mixing angles to be < \pi/4 to prevent the minimizer
  // from jumping between the 3+2 and 1+3+1 cases TODO: Find a better solution here?
  // FIXME Check if this should be fabs or no fabs
  if (n_flavors >= 5  &&
      (fabs(glbGetOscParamByName(params, "TH14")) > M_PI/4 ||
       fabs(glbGetOscParamByName(params, "TH24")) > M_PI/4 ||
       fabs(glbGetOscParamByName(params, "TH34")) > M_PI/4 ||
       fabs(glbGetOscParamByName(params, "TH15")) > M_PI/4 ||
       fabs(glbGetOscParamByName(params, "TH25")) > M_PI/4 ||
       fabs(glbGetOscParamByName(params, "TH35")) > M_PI/4))
  {
    pv += 1e11;
    goto my_prior_end;
  }

  // For 5-neutrino scenarios, the distinction between 3+2 and 1+3+1 is hardcoded
  // for compatibility with Thomas' code TODO: Find a better solution here?
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

  // If the user wants only positive mixing angles, add large penalty
  // in case of negative ones
  if (theta_positive)
  {
    if (n_flavors >= 3 &&
       (glbGetOscParamByName(params, "TH12") < 0 ||
        glbGetOscParamByName(params, "TH13") < 0 ||
        glbGetOscParamByName(params, "TH23") < 0))
    {
      pv += 1e15;
      goto my_prior_end;
    }
    if (n_flavors >= 4 &&
       (glbGetOscParamByName(params, "TH14") < 0 ||
        glbGetOscParamByName(params, "TH24") < 0 ||
        glbGetOscParamByName(params, "TH34") < 0))
    {
      pv += 2e15;
      goto my_prior_end;
    }
    if (n_flavors >= 5 &&
       (glbGetOscParamByName(params, "TH15") < 0 ||
        glbGetOscParamByName(params, "TH25") < 0 ||
        glbGetOscParamByName(params, "TH35") < 0 ||
        glbGetOscParamByName(params, "TH45") < 0))
    {
      pv += 3e15;
      goto my_prior_end;
    }
  }

#ifdef NU_USE_NUSQUIDS
  if (glbGetOscParamByName(params, "MA_OVER_M4") >= 1. ||
      glbGetOscParamByName(params, "MA_OVER_M4") < 0)
  {
    pv += 4e15;
    goto my_prior_end;
  }
  if (glbGetOscParamByName(params, "M4_GAMMA") < 0.)
  {
    pv += 5e15;
    goto my_prior_end;
  }
#endif

  // Add chi^2 from external codes
  // -----------------------------
  if (pp->ext_flags)
  {
    if (pp->ext_flags & EXT_REACTORS)
    {
      // Workaround for Thomas' reactor code allowing only certain ranges for
      // dm21sq, dm31sq, dm41sq, dm51sq
      if (n_flavors >= 3)
      {
        #ifndef LOG_ATM
          #error "my_prior works only for logarithmic Dmq_atm in reactor code."
        #endif
        double dm31 = fabs(glbGetOscParamByName(params, "DM31"));
        if (dm31 < POW10(ATM_MIN)+1e-10 || dm31 > POW10(ATM_MAX)-1e-10)
        {
          pv += 5e11;
          goto my_prior_end;
        }

        double dm21 = fabs(glbGetOscParamByName(params, "DM21"));
      #ifdef LOG_SOL
        if (dm21 < POW10(SOL_MIN)+1e-10 || dm21 > POW10(SOL_MAX)-1e-10)
      #else
        if (dm21 < SOL_MIN+1e-10 || dm21 > SOL_MAX-1e-10)
      #endif
        {
          pv += 6e11;
          goto my_prior_end;
        }
      }
      if (n_flavors >= 4)
      {
        double dm41 = fabs(glbGetOscParamByName(params, "DM41"));
      #ifdef Ip3pI // See Thomas' email from 2013-03-08
        if (dm41 < POW10(ATM_MIN)+1e-10 || dm41 > 0.5*POW10(STE_MAX)-1e-10)
      #else
        if (dm41 < POW10(ATM_MIN)+1e-10 || dm41 > POW10(STE_MAX)-1e-10)
      #endif
        {
          pv += 7e11;
          goto my_prior_end;
        }
      }
      if (n_flavors >= 5)
      {
        double dm51 = fabs(glbGetOscParamByName(params, "DM51"));
      #ifdef Ip3pI
        if (dm51 < POW10(ATM_MIN)+1e-10 || dm51 > 0.5*POW10(STE_MAX)-1e-10)
      #else
        if (dm51 < POW10(ATM_MIN)+1e-10 || dm51 > POW10(STE_MAX)-1e-10)
      #endif
        {
          pv += 8e11;
          goto my_prior_end;
        }
      } // if (n_flavors >= 5)
    } // if (ext_flags & EXT_REACTORS)

    // Prepare parameter data structure for Thomas' 2011 code
    struct params sbl_params;
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
      // Don't allow zero Ue4 since Thomas' code cannot handle this
      sbl_params.Ue[I4]  = MAX(1.e-7, gsl_complex_abs(_U[0][3]));
      sbl_params.Um[I4]  = gsl_complex_abs(_U[1][3]);
      sbl_params.Ue3     = gsl_complex_abs(_U[0][2]);
      sbl_params.dmq[I4] = fabs(glbGetOscParamByName(params, "DM41"));
    }
    if (n_flavors >= 5)
    {
      sbl_params.Ue[I5]  = MAX(1.e-7, gsl_complex_abs(_U[0][4]));
      sbl_params.Um[I5]  = gsl_complex_abs(_U[1][4]);
      sbl_params.dmq[I5] = fabs(glbGetOscParamByName(params, "DM51"));
      sbl_params.delta   = gsl_complex_arg(
          gsl_complex_mul(gsl_complex_mul(_U[1][3], gsl_complex_conjugate(_U[0][3])),
                          gsl_complex_mul(_U[0][4], gsl_complex_conjugate(_U[1][4]))) );
    }


    // Compute chi^2 using external codes
    if (pp->ext_flags & EXT_MB)
      pv += chi2mb475(sbl_params);
    if (pp->ext_flags & EXT_MB_300)
      pv += chi2mb300(sbl_params);
    if (pp->ext_flags & EXT_MBANTI)
      pv += chi2_MBA_475(sbl_params);
    if (pp->ext_flags & EXT_MBANTI_200)
      pv += chi2_MBA_200(sbl_params);
    if (pp->ext_flags & EXT_KARMEN)
      pv += chi2karmen(sbl_params);
    if (pp->ext_flags & EXT_LSND)
      pv += chi2lsnd(sbl_params);
    if (pp->ext_flags & EXT_NOMAD)
      pv += chi2nomad(sbl_params);
    if (pp->ext_flags & EXT_CDHS)
      pv += chi2cdhs(sbl_params);
    if (pp->ext_flags & EXT_ATM_TABLE)
      pv += chi2atm(sbl_params);
    if (pp->ext_flags & EXT_MINOS2016)
      pv += MINOS_2016_prior(in);
    if (pp->ext_flags & EXT_MINOS2017)
      pv += MINOS_2017_prior(in, 0);

    // Ivan's codes
    // ------------
  #ifdef NU_USE_NUSQUIDS
    if (pp->ext_flags & EXT_LSND_IVAN)
    {
      if (LSND_ivan_fitter)
      {
        extern regeneration::Param p_oscdecay;
        pv += LSND_ivan_fitter->getChisq(p_oscdecay);
      }
      else
        pv += 2e33;
    }
    if (pp->ext_flags & EXT_KARMEN_IVAN)
    {
      if (KARMEN_ivan_fitter)
      {
        extern regeneration::Param p_oscdecay;
        pv += KARMEN_ivan_fitter->getChisq(p_oscdecay);
      }
      else
        pv += 3e33;
    }
    if (pp->ext_flags & EXT_DECAY_KINEMATICS)
    {
      if (ext_constraints)
      {
        extern regeneration::Param p_oscdecay;
        pv += ext_constraints->getChisq_beta(p_oscdecay);
      }
      else
        pv += 4e33;
    }
    if (pp->ext_flags & EXT_FREE_STREAMING)
    {
      if (ext_constraints)
      {
        extern regeneration::Param p_oscdecay;
        pv += ext_constraints->getChisq_cosmology(p_oscdecay);
      }
      else
        pv += 5e33;
    }
  #endif

    // Joachim's MiniBooNE code
    // ------------------------
    if (pp->ext_flags & EXT_MB_JK)
    {
      if (debug_level > 0)
        pv += chiMB_jk(2);
      else
        pv += chiMB_jk(0);
    }

    // MicroBooNE
    // ----------
    if (pp->ext_flags & EXT_MUBOONE)
    {
      if (debug_level > 0)
        pv += chiMuBooNE(2);
      else
        pv += chiMuBooNE(0);
    }

    // SBL reactor experiments + gallium experiments
    // ---------------------------------------------
    if (pp->ext_flags & EXT_REACTORS)
    {
      // Prepare parameter data structure for Thomas' 2012 reactor code
      struct ns_reactor::Param_5nu reactor_params =
        { { NAN, NAN, NAN, 0, 0 },  // params for 4th and 5th neutrino need to be initialized
          { NAN, NAN, NAN, 0, 0 },  // even if we do only a 3-flavor fit!
          { NAN, NAN, NAN, 0, 0 } };
      static struct ns_reactor::Param_5nu last_reactor_params =
        { { NAN, NAN, NAN, NAN, NAN },
          { NAN, NAN, NAN, NAN, NAN },
          { NAN, NAN, NAN, NAN, NAN } };
      static double last_chi2_reactor = NAN;

      reactor_params.dmq[0] = 0.0;
      reactor_params.dmq[1] = glbGetOscParamByName(params, "DM21");
      double dmsq31 = glbGetOscParamByName(params, "DM31");
      if (dmsq31 < 0) // Thomas' code works only for NH (and the hierarchy doesn't matter for it)
        dmsq31 = -dmsq31 + reactor_params.dmq[1];
      reactor_params.dmq[2] = dmsq31;
      reactor_params.theta[ns_reactor::I12] = glbGetOscParamByName(params, "TH12");
      reactor_params.theta[ns_reactor::I13] = glbGetOscParamByName(params, "TH13");
      if (n_flavors >= 4)
      {
        reactor_params.dmq[3] = glbGetOscParamByName(params, "DM41");
        reactor_params.theta[ns_reactor::I14] = glbGetOscParamByName(params, "TH14");
      }
      if (n_flavors >= 5)
      {
        reactor_params.dmq[4] = glbGetOscParamByName(params, "DM51");
        reactor_params.theta[ns_reactor::I15] = glbGetOscParamByName(params, "TH15");
      }
      reactor_params.set_ang();

      if (memcmp(&last_reactor_params, &reactor_params, sizeof(reactor_params)) != 0)
      {     // Recompute reactor chi^2 only if the relevant parameters have changed
        last_reactor_params = reactor_params;
//        ns_reactor::fit.pull_status[ns_reactor::FLUX_NORM] = ns_reactor::FIXED;
//        ns_reactor::fit.set_experiments(sbl_exps);
        last_chi2_reactor = ns_reactor::fit.chisq(reactor_params);
        #ifdef USE_NEOS_ALVARO
          last_chi2_reactor += Neos_fitter.chi2(&reactor_params);
        #endif
        #ifdef USE_DB_ALVARO
          last_chi2_reactor += DB_fitter.chi2_syst(&reactor_params);
        #endif
        pv += last_chi2_reactor;
      }
      else
        pv += last_chi2_reactor;
    }


//    // Atmospheric neutrinos -- 2011 code (v43)
//    // ----------------------------------------
//    if (pp->ext_flags & EXT_ATM_COMP)
//    {
//      // Workaround for Michele's code requiring dm31sq < 1e-2
//      if (fabs(glbGetOscParamByName(params, "DM31")) > 0.99e-2)
//        pv += 1e15;
//      else if (glbGetNumOfOscParams() == 51+1)  // 3 flavors
//                               // Remember to change when adding parameters!
//        pv += atm_chisq(glbGetOscParamByName(params, "TH23"),
//                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
//                        glbGetOscParamByName(params, "DM31"));
//      else if (glbGetNumOfOscParams() == 92+5) // 4 flavors
//        pv += atm_chisq(glbGetOscParamByName(params, "TH23"),
//                        glbGetOscParamByName(params, "TH24"),
//                        0.0,
//                        glbGetOscParamByName(params, "TH34"),
//                        0.0,
//                        0.0, //TODO: Include phases
//                        0.0,
//                        glbGetOscParamByName(params, "DM31"));
//      else if (glbGetNumOfOscParams() == 145+7) // 5 flavors
//        pv += atm_chisq(glbGetOscParamByName(params, "TH23"),
//                        glbGetOscParamByName(params, "TH24"),
//                        glbGetOscParamByName(params, "TH25"),
//                        glbGetOscParamByName(params, "TH34"),
//                        glbGetOscParamByName(params, "TH35"),
//                        0.0,
//                        0.0,
//                        glbGetOscParamByName(params, "DM31"));
//      else
//        pv -= 1.e30;
//    } // pp->ext_flags & EXT_ATM_COMP


    // Atmospheric neutrinos -- 2012 code (v54) and 2017 code (v70)
    // ------------------------------------------------------------
    if (pp->ext_flags & EXT_ATM_COMP  ||  pp->ext_flags & EXT_DEEPCORE)
    {
#ifdef USE_ATM
      complx (*_U)[n_flavors] =
        (complx (*)[n_flavors]) gsl_matrix_complex_ptr(snu_get_U(), 0, 0);
      complx Uatm[5][5];
      for (int i=0; i < n_flavors; i++)
      {
        for (int j=0; j < n_flavors; j++)
          Uatm[i][j] = _U[i][j];
        for (int j=n_flavors; j < 5; j++)
          Uatm[i][j] = 0.0;
      }
      for (int i=n_flavors; i < 5; i++)
      {
        for (int j=0; j < 5; j++)
          Uatm[i][j] = 0.0;
        Uatm[i][i] = 1.0;
      }

      // Workaround for Michele's code requiring dm31sq < 1e-2
      if (fabs(glbGetOscParamByName(params, "DM31")) > 0.99e-2)
        pv += 1e15;
      else if (atm_decouple_e)
        pv += atm_chisq(Uatm, glbGetOscParamByName(params, "DM31"), 1.051);
      else
        pv += atm_chisq(Uatm, glbGetOscParamByName(params, "DM31"), 0.0);
#else
    printf("my_prior: Compiled without atmospherics support.\n");
    exit(-2);
#endif
    } // pp->ext_flags & EXT_ATM_COMP  ||  pp->ext_flags & EXT_DEEPCORE


    // Interface to Michele's solar neutrino code
    // ------------------------------------------
    if (pp->ext_flags & EXT_SOLAR)
    {
#ifdef USE_SOLAR
      static Angles last_a = { NAN, NAN, NAN, NAN, NAN, NAN, NAN };
      static double last_dm21 = NAN;
      static double last_chi2_solar = NAN;
      Angles a;
      double chi2_solar = -1.e25;
      double dm21 = glbGetOscParamByName(params, "DM21");

      if (dm21 < solar_dm21_min || dm21 > solar_dm21_max)
      {
        pv += 1.e26;
        goto my_prior_end;
      }

      solar_parameters(snu_get_U(), &a); // Determine parameters of Michele's parameterization
      if (memcmp(&last_a, &a, sizeof(a)) != 0) // Recompute probabilities only if mixing angles
      {                                        // have changed
        // If this is a degfinder prescan (indicated by systematics OFF), compute solar
        // probabilities only once TODO: Find a better solution here
        if (glb_num_of_exps == 0  ||  glbGetSysOnOffState(0,0)==GLB_ON  ||  isnan(last_a.the12))
        {
          if (debug_level > 2)
          {
            printf("# Recomputing solar probabilities ...\n");
            printf("# Old params: %g %g %g %g %g %g %g\n", last_a.the12, last_a.dlt12,
                   last_a.eta_e, last_a.ste_D, last_a.ste_N, last_a.cst_e, last_a.cst_a);
            printf("# New params: %g %g %g %g %g %g %g\n", a.the12, a.dlt12,
                   a.eta_e, a.ste_D, a.ste_N, a.cst_e, a.cst_a);
          }
          sun_probs(a, solar_dm21_min, solar_dm21_max);
          memcpy(&last_a, &a, sizeof(a));
        }

//        sun_chisq(dm21, &chi2_solar); // This is for the slower, non-adiabatic solar code
        chi2_solar = sun_chisq(dm21); // This is for the new, adiabatic solar code
        last_dm21 = dm21;
        last_chi2_solar = chi2_solar;
      }
      else if (last_dm21 != dm21)     // Recompute chi^2 only if the relevant angles/phases
      {                               // or dm21 have changed
//        sun_chisq(dm21, &chi2_solar);
        chi2_solar = sun_chisq(dm21);
        last_dm21 = dm21;
        last_chi2_solar = chi2_solar;
      }
      else
        chi2_solar = last_chi2_solar;

      pv += chi2_solar;
#else
    printf("my_prior: Compiled without solar support.\n");
#endif
    }
  } // if (pp->ext_flags)


  // Add oscillation parameter priors
  for (i=0; i < glbGetNumOfOscParams(); i++)
//FIXME    if (glbGetProjectionFlag(p,i) == GLB_FREE)
    {
      fitvalue     = glbGetOscParams(params,i);
      centralvalue = glbGetOscParams(central_values,i);
      inputerror   = glbGetOscParams(input_errors,i);
      if (inputerror > 1e-12)
        pv += SQR((centralvalue-fitvalue)/inputerror);
    }

  // Add matter parameter priors
  for (i=0; i < glb_num_of_exps; i++)
    if (glbGetDensityProjectionFlag(p,i) == GLB_FREE)
    {
      fitvalue     = glbGetDensityParams(params,i);
      centralvalue = 1.0;
      inputerror   = glbGetDensityParams(input_errors,i);
      if(inputerror > 1e-12)
        pv += SQR((centralvalue-fitvalue)/inputerror);
    }

  // Extra priors for MCMC mode
  // --------------------------
  if (action == NU_ACTION_MCMC)
  {
    for (i=0; i < pp->n_scan_params; i++)
    {
      // Add non-constant priors only if this parameter doesn't have a user-defined prior yet
      // NOTE: these priors are not properly normalized (but this shouldn't hurt the MCMC)
      if (glbGetOscParams(input_errors, pp->scan_params[i]) <= 1e-12)
      {
        double x = glbGetOscParams(params,pp->scan_params[i]);
        if (pp->scan_p_flags[i] & DEG_PLUS_MINUS)
          x = fabs(x);
        else if (x < 0)
          pv += 1.e12;
        else
        {
          if ((pp->scan_p_flags[i] & DEG_S22) && (pp->scan_p_flags[i] & DEG_LOGSCALE))
          {
            double s2x  = sin(2.*x);
            if (s2x != 0.)
              pv += -2. * log( (1./SQR(s2x)) * 2. * s2x * sqrt(1. - SQR(s2x)) * 2. );
          }
          else if (pp->scan_p_flags[i] & DEG_LOGSCALE)
          {
            if (x != 0.)
              pv += -2. * log(1./x);
          }
          else if (pp->scan_p_flags[i] & DEG_S22)
          {
            double s2x  = sin(2.*x);
            pv += -2. * log(  2. * s2x * sqrt(1. - SQR(s2x)) * 2. );
          }
        } // if (DEG_PLUS_MINUS)
      } // if (input_errors <= 1e-12)
    } // for (i)
  } // if (NU_ACTION_MCMC)

my_prior_end:
  glbFreeParams(central_values);
  glbFreeParams(input_errors);
  glbFreeParams(params);
  glbFreeProjection(p);
  return pv;
}

