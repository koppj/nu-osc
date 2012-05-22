/***************************************************************************
 * Functions for KARMEN \nu_e--C-12 scattering fit                         *
 ***************************************************************************
 * Author: PAN Machado                                                     *
 ***************************************************************************/
#include <stdio.h>
#include <math.h>
#include <globes/globes.h>
#include "nu.h"


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


/* this is used to reproduce their analysis */
static double JR(double true_rate, double fit_rate, double bkgd, double err)
{
      double tmp = (true_rate -bkgd - fit_rate);
      double tmp2 = true_rate + SQR(66.1) + SQR(0.7) /* SQR(err*fit_rate) */;
      return tmp/sqrt(tmp2);
}


double chi_karmen_c12_JR(int exp, int rule, int n_params,
                      double *x, double *errors, void *user_data)
{
      double KM_data_tot = 860;
      /* data from J. Reichenbacher thesis */
      double *signal_fit_rates = glbGetSignalFitRatePtr(exp, 0);
      double signal_rates[26];
      double signal_norm;
      double chi2 = 0.0;
      int  i,ew_low, ew_high;

      glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
      signal_norm = 1.0 + x[0];
      for(i=ew_low; i<=ew_high; i++)
        signal_rates[i] = signal_fit_rates[i];

      chi2 = JR(KM_data_tot,signal_rates[0], 13.9, 0.0825);
      return chi2;
}

double chi_karmen_c12(int exp, int rule, int n_params,
                      double *x, double *errors, void *user_data)
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

