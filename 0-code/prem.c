/* ---------------------------------------------------------------------------- */
/* Handling of matter profiles                                                  */
/* ---------------------------------------------------------------------------- */
/* Author: Joachim Kopp                                                         */
/* ---------------------------------------------------------------------------- */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#include <globes/globes.h>
#include "nu.h"

#define PREM_RES 10 /* km */

/* Global variables */
extern const int debug_level;

/* Local variables */
static double prem_data[(int)(REARTH/PREM_RES)+2][2];


/***************************************************************************
 *         H A N D L I N G   O F   M A T T E R   P R O F I L E S           *
 **************************************************************************/

/***************************************************************************
 * Read the PREM earth matter density profile from a file                  *
 ***************************************************************************/
int LoadPREMProfile(const char *prem_file)
{
  /* Load PREM profile from file */
  FILE *f = fopen(prem_file, "r");
  if (!f)
    return -1;

  int i = 0;
  char line[100];
  double raw_data[1000][2];
  while (fgets(line, 100, f))
  {
    if (line[0] != '#')       /* Ignore comments */
    {
      sscanf(line, "%lg %lg", &raw_data[i][0], &raw_data[i][1]);
      i++;
    }
  };
  fclose(f);

  /* Interpolate PREM profile */
  i = 0;
  int k;
  double r;
  for (k=0, r=0.0;  r <= REARTH;  k++, r += PREM_RES)
  {
    while (raw_data[i][0] <= r)
      i++;
    prem_data[k][0] = r;
    prem_data[k][1] = raw_data[i-1][1] + (raw_data[i][1] - raw_data[i-1][1])
                  * (r - raw_data[i-1][0])/(raw_data[i][0] - raw_data[i-1][0]);
  }
  prem_data[k][0] = r;                 /* We need one point a R > REarth          */
  prem_data[k][1] = prem_data[k-1][1]; /* This approximation should be sufficient */

  if (debug_level > 2)
    for (int i=0; i <= k; i++)
      printf("# PREM: %d %10.5g %10.5g\n", i, prem_data[i][0], prem_data[i][1]);

  return 0;
}


/***************************************************************************
 * Calculate the matter density at position t along a trajectory of        *
 * length L, assuming the PREM earth matter density profile                *
 ***************************************************************************/
double GetPREMDensity(double t, double L)
{
  if (L > 2*REARTH || t > L)
    return -1.0;
  
  double r = sqrt(SQR(t) + SQR(REARTH) - t*L);  /* Distance from earth center */
  int m = (int)(r/PREM_RES) + 1;
  double rho = prem_data[m-1][1] + (prem_data[m][1] - prem_data[m-1][1]) *
                 (r - prem_data[m-1][0])/(prem_data[m][0] - prem_data[m-1][0]);

  if (r > RCORE)   /* Use different effictive electron numbers for mantle and core */
    return rho;
  else
    return rho * Ne_CORE / Ne_MANTLE;
}


/***************************************************************************
 * Wrapper function for GetPREMDensity, required for numerical integration *
 ***************************************************************************/
static double GSLCallbackPREM(double t, void *L)
{
  return GetPREMDensity(t, *((double *)L));
}


/***************************************************************************
 * Return the averaged matter potential over part of the neutrino          *
 * trajectory                                                              *
 ***************************************************************************/
double GetAvgPREMDensity(double Ltot, double Lstart, double Lend)
{
  if (Ltot > 2*REARTH || Lstart > Lend || Lstart < 0.0)
    return -1.0;
  
  const int max_intervals = 1000;
  double eps_abs = 1000.0;
  double eps_rel = 1e-2;
  double quad_result, quad_error;
  gsl_integration_workspace *ws = gsl_integration_workspace_alloc(max_intervals);
  gsl_function f;
  f.function = &GSLCallbackPREM;
  f.params   = &Ltot;

  gsl_integration_qag(&f, Lstart, Lend, eps_abs, eps_rel, max_intervals,
      GSL_INTEG_GAUSS21, ws, &quad_result, &quad_error);

  return quad_result / (Lend - Lstart);
}

/***************************************************************************
 * Compute density data for a 3-layer approximation to the PREM profile    *
 ***************************************************************************/
int GetPREM3LayerApprox(double L, int *n_layers, double *lengths,
                        double *densities)
{
  /* If the trajectory does not pass the core, use only one layer */
  if (SQR(REARTH) - SQR(L)/4.0 > SQR(RCORE))
  {
    *n_layers    = 1;
    lengths[0]   = L;
    densities[0] = GetAvgPREMDensity(L, 0.0, L);
  }
  else
  {
    double LMantle = L/2.0 - sqrt(SQR(L)/4.0 - SQR(REARTH) + SQR(RCORE));
    double LCore   = L - 2.0*LMantle;
    *n_layers      = 3;
    lengths[0]     = LMantle;
    lengths[1]     = LCore;
    lengths[2]     = LMantle;
    densities[0]   = GetAvgPREMDensity(L, 0.0, LMantle);
    densities[1]   = GetAvgPREMDensity(L, LMantle, LMantle+LCore);
    densities[2]   = densities[0];
  }

  return 0;
}



