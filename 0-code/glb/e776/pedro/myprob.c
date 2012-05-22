/***************************************************************************
 * Probability engine for simulating large extra dimensions                *
 ***************************************************************************
 * Author: Pedro A N Machado                                               *
 * Date: 2010-09-21                                                        *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_eigen.h>
#include "globes/globes.h"
#include "myheader.h"

/* Constants */
#define GLB_V_FACTOR        7.5e-14   /* Conversion factor for matter potentials */
#define GLB_Ne_MANTLE       0.5        /* Effective electron numbers for calculation */
#define GLB_Ne_CORE         0.468      /*   of MSW potentials                        */
#define V_THRESHOLD 0.001*GLB_V_FACTOR*GLB_Ne_MANTLE  /* The minimum matter potential below */
                                                      /* which vacuum algorithms are used   */
#define M_SQRT3  1.73205080756887729352744634151      /* sqrt(3) */
#define matter_factor 18.9431827e-15 /* matter_factor * density (in g/cm3) is Vcc/2 in eV*/

/* Macros */
#define SQR(x)      ((x)*(x))                        /* x^2   */
#define SQR_ABS(x)  (SQR(creal(x)) + SQR(cimag(x)))  /* |x|^2 */

/* Fundamental oscillation parameters */
static double th12; // Mixing angles
static double ss2t;
static double c2;
static double s2;
static double dm;            // Squared masses


/***************************************************************************
 * Function twoflv_set_oscillation_parameters                                 *
 ***************************************************************************
 * Sets the fundamental oscillation parameters and precomputes the mixing  *
 * matrix and part of the Hamiltonian.                                     *
 ***************************************************************************/
int twoflv_set_oscillation_parameters(glb_params p, void *user_data)
{
  /* Copy standard oscillation parameters */
  th12  = glbGetOscParams(p, GLB_THETA_12);
  dm    = glbGetOscParams(p, GLB_DM_21);
  
  /* Compute vacuum mixing matrix */
  c2 = SQR(cos(th12));
  s2 = SQR(sin(th12));
  //ss2t = SQR(sin(2.0*th12));
  //printf("%f\n",dm);
  return 0;
}

/***************************************************************************
 * function led_get_oscillation_parameters                                 *
 ***************************************************************************
 * Returns the current set of oscillation parameters.                      *
 ***************************************************************************/
int twoflv_get_oscillation_parameters(glb_params p, void *user_data)
{
  glbDefineParams(p, th12, 0, 0, 0, dm, 0);
  return 0;
}

/***************************************************************************
 * Function twoflv_probability_matrix                                      *
 ***************************************************************************
 * Calculates the neutrino oscillation probability matrix.                 *
 ***************************************************************************
 * Parameters:                                                             *
 *   P:       Buffer for the storage of the matrix                         *
 *   cp_sign: does not matter here                                         *
 *   E:       Neutrino energy (in GeV)                                     *
 *   psteps:  Number of layers in the matter density profile               *
 *   length:  Lengths of the layers in the matter density profile in km    *
 *   density: not used                                                     *
 *   filter_sigma: Width of low-pass filter or <0 for no filter            *
 *   user_data: Unused here, should be NULL                                *
 ***************************************************************************/
int twoflv_probability_matrix(double P[3][3], int cp_sign, double energy,
    int psteps, const double *length, const double *density,
    double filter_sigma, void *user_data)
{
  int i,j;
  /* Convert energy to eV */
  double E = energy*1.0e9;
  double L = 0.0; for(int i=0; i < psteps; i++) L += length[i];
  L = GLB_KM_TO_EV(L);

  double prob=0, phase;
  phase = dm*L/2.0/E;
  prob = filter_sigma<0 ? 
    1.0 - 2.0*c2*s2*(1.0 - cos(phase)) :
    1.0 - 2.0*c2*s2*(1.0 - cos(phase)*exp(-SQR(phase/filter_sigma))) ;
    /* 2.0*c2*(c2 + s2*cos(phase)) : */
    /* 2.0*c2*(c2 + s2*cos(phase)*exp(-SQR(phase/filter_sigma))); */

  for(i=0; i<2; i++)
    for(j=0; j<2; j++)
      P[i][j] = i==j ? prob : 1-prob;
  if(P[1][0]<0)printf("P = %e\n",P[1][0]);

  return 0;
}
