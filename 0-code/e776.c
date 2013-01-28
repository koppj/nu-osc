/***************************************************************************
 * Functions for E776 fit                                                  *
 ***************************************************************************
 * Author: PAN Machado                                                     *
 ***************************************************************************/
#include <stdio.h>
#include <math.h>
#include <globes/globes.h>
#include "nu.h"

#define E776BINS (14) /* matches E776.glb */


/* E776 data from PRL 68 274 (1992):
   observed events for + (0) and - (1) polarities,
   pi0 background for + (2) and - (3) polarities  */
double E776data[4][E776BINS] =
 {{ 35, 45, 23, 6, 6, 4, 1, 5, 5, 3, 1, 2, 0,0 },   /* + polarity, 136 total */
  {  9, 18,  6, 3, 4, 3, 1, 1, 0, 1, 0, 1, 0,0 },   /* - polarity, 47 total */
  {23.1, 43.4, 15.0, 6.28, 3.14, 0.92, 2.77, 0,0,0,0,0,0,0 }, /* + pol, pi0 94.61 */
  {9.01, 26.8, 4.69, 2.85, 0,0,0,0,0,0,0,0,0,0 } };           /* - pol, pi0 43.35 */

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

double chi_E776(int exp, int rule, int n_params, double *x,
                double *errors, void *user_data)
{
  double *data = (double *) E776data; /* for both modes */
  /* this is pointing to E776data[4][E776BINS], so attention when reading */
  double *signal_fit_rates_pos = glbGetSignalFitRatePtr(exp, 0);
  double *signal_fit_rates_neg = glbGetSignalFitRatePtr(exp, 1);
  double *bg_fit_rates_pos     = glbGetBGFitRatePtr(exp, 0);
  double *bg_fit_rates_neg     = glbGetBGFitRatePtr(exp, 1);
  double pi0_norm_pos,pi0_norm_neg, nue_beam_norm_pos, nue_beam_norm_neg;
  double fit_rate_pos=0, fit_rate_neg=0;
  double datapos=0,dataneg=0;
  int ew_low, ew_high;
  double chi2 = 0.0;
  int i;

  glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
  pi0_norm_pos  = 1.0 + x[0]; /* pi0 bkgd error + */
  pi0_norm_neg  = 1.0 + x[1]; /* pi0 bkgd error - */
  nue_beam_norm_pos = 1.0 + x[2]; /* beam error + */
  nue_beam_norm_neg = 1.0 + x[3]; /* beam error - (not correlated?) */

  /* I am summing the first bins. The minimum chi2 gets high
     when I fit them separetely. */
  for(i=1; i<3; i++) /* Avoid the first bin... */
  {
    fit_rate_pos += signal_fit_rates_pos[i]
        + pi0_norm_pos*data[2*E776BINS+i] /* read pi0 + polarity */
        + nue_beam_norm_pos*bg_fit_rates_pos[i];
      fit_rate_neg += signal_fit_rates_neg[i]
        + pi0_norm_neg*data[3*E776BINS+i] /* read pi0 - polarity */
        + nue_beam_norm_neg*bg_fit_rates_neg[i];
      datapos += data[0*E776BINS+i];
      dataneg += data[1*E776BINS+i];
  }
  chi2 += glb_likelihood(datapos, fit_rate_pos);
  chi2 += glb_likelihood(dataneg, fit_rate_neg);

  for (i=3; i <= ew_high; i++)
  {
      fit_rate_pos = signal_fit_rates_pos[i]
        + pi0_norm_pos*data[2*E776BINS+i] /* read pi0 + polarity */
        + nue_beam_norm_pos*bg_fit_rates_pos[i];
      fit_rate_neg = signal_fit_rates_neg[i]
        + pi0_norm_neg*data[3*E776BINS+i] /* read pi0 - polarity */
        + nue_beam_norm_neg*bg_fit_rates_neg[i];
      chi2 += glb_likelihood(data[0*E776BINS+i], fit_rate_pos);
      chi2 += glb_likelihood(data[1*E776BINS+i], fit_rate_neg); /* antinu */
  }

  /* I have 4 errors: pi0 background and beam nue for each polarity */
  for(i=0; i<4; i++)chi2 += glb_prior(x[i], 0.0, errors[i]);

  if (isnan(chi2))
    chi2=1.e102;
  return chi2;
}


/* I am not using this chi2 */
double chi_E776_rates(int exp, int rule, int n_params, double *x,
                double *errors, void *user_data)
{
      double *data = (double *) E776data; /* for both modes */
      /* this is pointing to E776data[4][E776BINS], so attention when reading */
      double *signal_fit_rates_pos = glbGetSignalFitRatePtr(exp, 0);
      double *signal_fit_rates_neg = glbGetSignalFitRatePtr(exp, 1);
      double *bg_fit_rates_pos     = glbGetBGFitRatePtr(exp, 0);
      double *bg_fit_rates_neg     = glbGetBGFitRatePtr(exp, 1);
      double pi0_norm_pos,pi0_norm_neg, nue_beam_norm_pos, nue_beam_norm_neg;
      int ew_low, ew_high;
      double fit_rate_pos, fit_rate_neg;
      double chi2 = 0.0;
      int i;
      fit_rate_pos = 0;
      fit_rate_neg = 0;
      double sum_data_pos = 0;
      double sum_data_neg = 0;

      glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
      pi0_norm_pos  = 1.0 + x[0]; /* pi0 bkgd error + */
      pi0_norm_neg  = 1.0 + x[1]; /* pi0 bkgd error - */
      nue_beam_norm_pos = 1.0 + x[2]; /* beam error + */
      nue_beam_norm_neg = 1.0 + x[3]; /* beam error - (not correlated?) */
      for (i=2; i <= ew_high; i++)
      {
            fit_rate_pos += signal_fit_rates_pos[i]
              + pi0_norm_pos*data[2*E776BINS+i] /* read pi0 + polarity */
              + nue_beam_norm_pos*bg_fit_rates_pos[i];
            fit_rate_neg += signal_fit_rates_neg[i]
              + pi0_norm_neg*data[3*E776BINS+i] /* read pi0 - polarity */
              + nue_beam_norm_neg*bg_fit_rates_neg[i];
            sum_data_pos +=data[0*E776BINS+i];
            sum_data_neg +=data[1*E776BINS+i];
      }

      chi2 += glb_likelihood(sum_data_pos, fit_rate_pos);
      chi2 += glb_likelihood(sum_data_neg, fit_rate_neg); /* antinu */
      for(i=0; i<4; i++)chi2 += glb_prior(x[i], 0.0, errors[i]);
      return chi2;
}

