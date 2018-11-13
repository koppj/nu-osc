/***************************************************************************
 * Functions for KARMEN \nu_e--C-12 scattering fit                         *
 ***************************************************************************
 * Author: PAN Machado                                                     *
 ***************************************************************************/
#include <stdio.h>
#include <math.h>
#include <globes/globes.h>
#include <gsl/gsl_cdf.h>
#include "nu.h"


// Static global variables
static int NUC_first = -1; // Index of KARMEN experiment in GLoBES experiment list
static double NUC_error[] = { 0.075, 0.099, 0.12};


// Helper functions
inline static double glb_prior(double x, double center, double sigma)
{
      double tmp = (x - center)/sigma;
      return tmp*tmp;
}


inline static double glb_likelihood(double true_rate, double fit_rate)
{
      double res;
      res = fit_rate - true_rate;
      if (true_rate > 0)
      {
            if (fit_rate <= 0.0)
              res = 1e100;
            else
              res += true_rate * log( true_rate / fit_rate );
      }
      else
            res = fabs(res);
      return 2.0 * res;
}


/* this is used to reproduce JR's analysis */
static double JR(double true_rate, double fit_rate, double bkgd, double err)
{
      double tmp = (true_rate -bkgd - fit_rate);
      double tmp2 = true_rate + SQR(66.1) + SQR(0.7) /* SQR(err*fit_rate) */;
//      return tmp/sqrt(tmp2);  // Compute # of standard deviations
      return SQR(tmp)/tmp2;  // Compute chi^2
}


// -------------------------------------------------------------------------
double chi_karmen_c12_JR(int exp, int rule, int n_params,
                      double *x, double *errors, void *user_data)
// -------------------------------------------------------------------------
// chi^2 function for KARMEN's C-12 analysis, based on J. Reichenbacher's
// thesis
// -------------------------------------------------------------------------
{
      double KM_data_tot = 860;
      /* data from J. Reichenbacher thesis */
      double *signal_fit_rates = glbGetSignalFitRatePtr(exp, 0);
      double signal_rates[26];
//      double signal_norm; // FIXME Why is this not used?
      double chi2 = 0.0;
      int  i,ew_low, ew_high;

      glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
//      signal_norm = 1.0 + x[0];
      for(i=ew_low; i<=ew_high; i++)
        signal_rates[i] = signal_fit_rates[i];

      chi2 = JR(KM_data_tot,signal_rates[0], 13.9, 0.0825);
      return chi2;
}


// -------------------------------------------------------------------------
double chi_karmen_c12(int exp, int rule, int n_params,
                      double *x, double *errors, void *user_data)
// -------------------------------------------------------------------------
// chi^2 function for KARMEN's C-12 analysis - Pedro's likelihood fit
// -------------------------------------------------------------------------
{
      /* double KM_data[26] = { 14, 20, 18, 33, 24, 34, 40, 36, 45, 34,  */
      /*                             49, 60, 54, 56, 51, 39, 49, 53, 32, 36, */
      /*                             27, 17, 21, 10, 5, 3}; */
      double KM_data_tot = 860;
      /* data from J. Reichenbacher thesis */
      double *signal_fit_rates = glbGetSignalFitRatePtr(exp, 0);
      double signal_tot;
      double signal_norm;
      double chi2 = 0.0;
      int  i,ew_low, ew_high;

      glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
      signal_norm = 1.0 + x[0];

      signal_tot = 0;
      for(i=ew_low; i<=ew_high; i++)
        signal_tot += signal_fit_rates[i];
      signal_tot *= signal_norm;
      signal_tot += 13.9; /* background */

      chi2 = glb_likelihood(KM_data_tot,signal_tot);
      chi2 += glb_prior(x[0], 0.0, errors[0]);
      return chi2;
}


// -------------------------------------------------------------------------
double chi_lsnd_c12(int exp, int rule, int n_params,
                      double *x, double *errors, void *user_data)
// -------------------------------------------------------------------------
// chi^2 function for LSND's C-12 analysis - Pedro's likelihood fit
// -------------------------------------------------------------------------
{
      /* double LS_data[5] = { 11.25, 17.50, 21.25, 25.00, 35.00 }; */
      const double LS_data[12] = { 67, 75, 100, 106, 103, 98, 66, 66, 30, 23, 2, 0 };
//      const double LS_data_tot = 733;

      /* data from hep-ex/0105068 */
      double *signal_fit_rates = glbGetSignalFitRatePtr(exp, 0);
      double signal_rates[12];
      double signal_norm;
      double chi2 = 0.0;
      int  i,ew_low, ew_high;

      glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
      signal_norm = 1.0 + x[0];
      for(i=ew_low; i<=ew_high; i++)
        {
          signal_rates[i] = signal_norm * signal_fit_rates[i];
          /* chi2 += glb_likelihood(LS_data[i],signal_rates[i]); */
        }

 /* We have 12 bins. I will add some of them because the analysis with
    12 bins is rather poor (there should be some distortion which we
    cannot take into account. Set B = 2 for 6 bins, B = 3 for 4 bins,
    or B = 4, 6, 12. When we decide our analysis, this can be optimized. */
      int B=2;
      for(i=ew_low; i<=ew_high; i+=B)
        {
          double a=0,b=0;
          for(int j=i; j<i+B; j++)
            {
              a+=LS_data[j];
              b+=signal_rates[j];
            }
          chi2 += glb_likelihood(a,b);
        }

      chi2 += glb_prior(x[0], 0.0, errors[0]);
      return chi2;
}


// -------------------------------------------------------------------------
double chi_nue_carbon(int exp, int rule, int n_params,
                      double *x, double *errors, void *user_data)
// -------------------------------------------------------------------------
// chi^2 function for combined analysis of KARMEN's and LSND's C-12 data,
// using KARMEN rate and LSND spectral information.
// taken from Pedro's code
// -------------------------------------------------------------------------
{
  /* KARMEN: data from J. Reichenbacher thesis */
  double KM_data_tot = 860;
  double *KMsignal_fit_rates = glbGetSignalFitRatePtr(NUC_first, 0);
  double KMsignal_tot;
  double KMsignal_norm;

  /* LSND: data from hep-ex/0105068 */
  double LS_data[12] = { 67, 75, 100, 106, 103, 98, 66, 66, 30, 23, 2, 0 };
  double *LSsignal_fit_rates = glbGetSignalFitRatePtr(NUC_first+1, 0);
  double LSsignal_rates[12];
  double LSsignal_norm;

  /* common */
  double chi2 = 0.0;
  int  i,ew_low, ew_high;

  /* KARMEN */
  glbGetEnergyWindowBins(NUC_first, rule, &ew_low, &ew_high);
  KMsignal_norm = 1.0 + x[0] + x[2];
  KMsignal_tot = 0;
  for(i=ew_low; i<=ew_high; i++)
    KMsignal_tot += KMsignal_fit_rates[i];
  KMsignal_tot *= KMsignal_norm;
  KMsignal_tot += 13.9; /* KARMEN background */
  chi2 = glb_likelihood(KM_data_tot,KMsignal_tot);

  /* LSND */
  glbGetEnergyWindowBins(NUC_first+1, rule, &ew_low, &ew_high);
  LSsignal_norm = 1.0 + x[1] + x[2];
  for(i=ew_low; i<=ew_high; i++)
    LSsignal_rates[i] = LSsignal_norm * LSsignal_fit_rates[i];
  /* We have 12 bins. I will add some of them because the analysis with
     12 bins is rather poor (there should be some distortion which we
     cannot take into account. Set B = 2 for 6 bins, B = 3 for 4 bins,
     or B = 4, 6, 12. When we decide our analysis, this can be optimized. */
  int B=2;
  for(i=ew_low; i<=ew_high; i+=B)
    {
      double a=0,b=0;
      for(int j=i; j<i+B; j++)
        {
          a+=LS_data[j];
          b+=LSsignal_rates[j];
        }
      chi2 += glb_likelihood(a,b);
    }

  /* systematics priors */
  for(i=0; i<n_params; i++)
    chi2 += glb_prior(x[i], 0.0, errors[i]);
  return chi2;
}


// -------------------------------------------------------------------------
double chi_nue_carbon_spectrum(int exp, int rule, int n_params,
                               double *x, double *errors, void *user_data)
// -------------------------------------------------------------------------
// chi^2 function for combined analysis of KARMEN's and LSND's C-12 data,
// using spectral data for both KARMEN and LSND.
// taken from Pedro's code
// -------------------------------------------------------------------------
{
  /* KARMEN: data from J. Reichenbacher thesis */
  double KM_data[26] = { 14, 20, 18, 33, 24, 34, 40, 36, 45, 34,
                         49, 60, 54, 56, 51, 39, 49, 53, 32, 36,
                         27, 17, 21, 10, 5, 3};
  double KM_bkgd[26] = { 1.8193,   1.55312,  1.02658,  0.63072,  0.887069,
                         0.621092, 0.615869, 0.480769, 0.476157, 0.471734,
                         0.466314, 0.4611,   0.325999, 0.452668, 0.447651,
                         0.442833, 0.437215, 0.432594, 0.427974, 0.423156,
                         0.156577, 0.152163, 0.278427, 0.141924, 0.137906,
                         0.133088 };
  double *KMsignal_fit_rates = glbGetSignalFitRatePtr(NUC_first, 0);
  double KMsignal_rates[26];
  double KMsignal_norm;

  /* LSND: data from hep-ex/0105068 */
  double LS_data[12] = { 67, 75, 100, 106, 103, 98, 66, 66, 30, 23, 2, 0 };
  double *LSsignal_fit_rates = glbGetSignalFitRatePtr(NUC_first+1, 0);
  double LSsignal_rates[12];
  double LSsignal_norm;

  /* common */
  double chi2 = 0.0;
  int  i,ew_low, ew_high;

  /* KARMEN */
  glbGetEnergyWindowBins(NUC_first, rule, &ew_low, &ew_high);
  KMsignal_norm = 1.0 + x[0] + x[2];

  for(i=ew_low; i<=ew_high; i++)
    {
      KMsignal_rates[i] = KMsignal_norm * KMsignal_fit_rates[i] + KM_bkgd[i];
      chi2 += glb_likelihood(KM_data[i],KMsignal_rates[i]);
    }

  /* LSND */
  glbGetEnergyWindowBins(NUC_first+1, rule, &ew_low, &ew_high);
  LSsignal_norm = 1.0 + x[1] + x[2];
  for(i=ew_low; i<=ew_high; i++)
    LSsignal_rates[i] = LSsignal_norm * LSsignal_fit_rates[i];
  /* We have 12 bins. I will add some of them because the analysis with
     12 bins is rather poor (there should be some distortion which we
     cannot take into account. Set B = 2 for 6 bins, B = 3 for 4 bins,
     or B = 4, 6, 12. When we decide our analysis, this can be optimized. */
  int B=2;
  for(i=ew_low; i<=ew_high; i+=B)
    {
      double a=0,b=0;
      for(int j=i; j<i+B; j++)
        {
          a+=LS_data[j];
          b+=LSsignal_rates[j];
        }
      chi2 += glb_likelihood(a,b);
    }

  /* systematics priors */
  for(i=0; i<n_params; i++)
    chi2 += glb_prior(x[i], 0.0, errors[i]);
  return chi2;
}


// -------------------------------------------------------------------------
void init_nue_carbon(int KM_spectrum)
// -------------------------------------------------------------------------
// Initialize KARMEN/LSND joint C-12 analysis by defining chi^2 function
// and setting appropriate number of bins. Call this *before* loading the
// AEDL file. KM_spectrum = 1 ... use spectral info KARMEN and LSND,
// KM_spectrum = 0 ... use spectrum only for LSND.
// -------------------------------------------------------------------------
{
  NUC_first = glb_num_of_exps;

  /* KM_spectrum = 0 is for KARMEN rate only analysis */
  if(KM_spectrum==0)
  {
    glbDefineChiFunction(&chi_nue_carbon, 3, "chi_nue_carbon", NULL);
    glbDefineAEDLVariable("KMBINS",1);
  }
  else
  {
    glbDefineChiFunction(&chi_nue_carbon_spectrum, 3,"chi_nue_carbon", NULL);
    glbDefineAEDLVariable("KMBINS",26);
  }

  glbInitExperiment("c12-combi.glb",&glb_experiment_list[0],&glb_num_of_exps);

  glbSetChiFunction(NUC_first, 0, GLB_ON, "chi_nue_carbon", NUC_error);
  return;
}

