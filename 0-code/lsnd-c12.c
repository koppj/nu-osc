/***************************************************************************
 * Functions for LSND \nu_e--C-12 scattering fit                           *
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


double chi_lsnd_c12(int exp, int rule, int n_params,
                      double *x, double *errors, void *user_data)
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

