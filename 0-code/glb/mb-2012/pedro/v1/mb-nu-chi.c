/***************************************************************************
 * Functions for MiniBooNE numu to nue fit                                 *
 ***************************************************************************
 * Author: Joachim Kopp                                                    *
 * Adapted/modified by PAN Machado - 2012-Jul                              *
 * See http://www-boone.fnal.gov/for_physicists/data_release/lowe/         *
 ***************************************************************************/
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <globes/globes.h>
#include "nu.h"
#include "glb_probability.h"

#define NU_FLAVOURS 3
// Use MC events to compute signal rate (very slow)?
/* I am not using this anymore (always use the MC events) */
#define USE_MC_EVENTS

/* Neutrino mode parameters */
int NEnu    = 0;             /* Number of E_nue bins      */
int NMUnu   = 8;             /* Number of E_numue bins    */
int NCOVnu  = 0;             /* Size of covariance matrix */
/* same for antineutrino mode */
int NEbar   = 0;             
int NMUbar  = 8;
int NCOVbar = 0;   

int MB_EXP  = -1;

#define TOTBINS (19)
#define NUEBINS (11)
#define MUBINS (8)

/************************ MINIBOONE DATA ************************/
/* Common to both modes */
const double numubin[MUBINS] = { 500, 700, 900, 1100, 1300, 1500, 1700, 1900 };

/* NU MODE: */
double ALLdatanu[TOTBINS]   = { 232,  156, 156, 79, 81, 71, 64, 66, 63, 35,  70,
				19189, 29943, 26748, 18692, 11123, 5830, 2674, 1273 };
double pred_numu[MUBINS] = {19123.3,30011.7,26840.2,19364.,11753.3,6288.31,2953.79,1362.23};
double ALLbgnu[NUEBINS]     = { 181.06, 108.41, 120.36, 64.22, 90.32, 
				67.69,  70.43,  57.48,  52.34, 39.01, 70.23 };
/* separate beam nue background from others: */
double OTHERbgnu[NUEBINS]  = {167.393,92.009,92.342,45.7692,60.4225,
			      37.7925,40.703,33.903,31.84,25.342,49.725};
double NUEbgnu[NUEBINS]  = {13.667,16.401,28.018,18.4508,29.8975,
			      29.8975,29.727,23.577,20.5,13.668,20.505};

/* NUBAR MODE: */
double ALLdatabar[TOTBINS]  = { 122, 70, 65, 43, 57, 39, 38, 23, 27, 26, 49,
  				9149.81, 13406.3, 10882.7, 7552.63, 4498.04, 2285.95, 922.367, 351.378 };
double pred_numubar[MUBINS] = { 4952, 6637, 5537, 3715, 
				2115, 1042, 439.3, 167.5 };
double ALLbgbar[NUEBINS]    = { 44.38, 26.6, 28.45, 16.46, 21.41,
				17.23, 17.73, 14.3, 14.11, 10.00, 21.03 };
double OTHERbgbar[NUEBINS]  = {37.2549,19.1002,16.575,7.46022,6.7228,
			      4.41736,3.48032,2.48763,2.36024,1.24954,11.6534};
double NUEbgbar[NUEBINS]  = {7.12509,7.49975,11.875,8.99978,14.6872,
			    12.8126,14.2497,11.8124,11.7498,8.75046,9.37657};
/* fraction of nue/nue+nuebar of beam bg */
/* FIXME - I NEED BETTER ESTIMATIONS FOR THIS!!! */
double nuefraction[NUEBINS] = {0.492247,0.503388,0.540815,0.548179,0.56201,
			       0.602036,0.619732,0.657237,0.730518,0.77557,1.01635};
/* fraction of numu/numu+numubar (of numu oscillates, evts*numufraction*prob should be added to the signal) */
double numufraction[NUEBINS] = {0.20, 0.20, 0.20, 0.20, 0.20, 
				0.20, 0.20, 0.20, 0.20, 0.20, 0.20 };
/****************************************************************/

double *datanu, *bgnu, *nosc_bgnu, *nuebeam_bgnu, *bincenternu;
double *databar, *bgbar, *nosc_bgbar, *nuebeam_bgbar, *bincenterbar,
  *nuefrac,*numufrac;

double binwidths[NUEBINS] = { 0.1, 0.075, 0.1, 0.075, 0.125, 0.125, 0.15, 0.15, 0.2, 0.2, 1.5 };
double bincenter[NUEBINS];
double *binwnu,*binwbar;

/* for binning the MC data */
#define LBINS (15) 
#define EBINS (70)
double Eminnu,Eminbar,Euptrue[EBINS],Ecenter[EBINS],Lup[LBINS],Lcenter[LBINS];

static gsl_matrix *M3nu        = NULL;
static gsl_matrix *M2nu        = NULL;
static gsl_matrix *M2invnu     = NULL;
static gsl_permutation *permnu = NULL;

static gsl_matrix *M3bar        = NULL;
static gsl_matrix *M2bar        = NULL;
static gsl_matrix *M2invbar     = NULL;
static gsl_permutation *permbar = NULL;
/* #ifdef USE_MC_EVENTS */
  static long Neventsnu;
  static double mc_eventsnu[128000][4];
  static long Neventsbar;
  static double mc_eventsbar[128000][4];
  static glb_params osc_params = NULL;
/* #endif */

int mb_mode;
double chiMB(int exper, int rule, int n_params, double *x, double *errors,
	     void *user_data);

/***************************************************************************
 * Initialize GSL data structures required for MiniBooNE chi^2 calculation *
 * numode and nubarmode = 0,1,2:                                           *
 * 0 - mode off                                                            *
 * 1 - Enu > 200 MeV                                                       *
 * 2 - Enu > 475 MeV                                                       *
 ***************************************************************************/
int chiMB_init(int numode, int nubarmode)
{
  int i;

  switch(numode)
    {
    case 0:			/* nu mode off */
      NEnu   = -1;
      Eminnu = 1E10;
      break;
    case 1:			/* nu mode E>200 MeV */
      NEnu   = 11;
      Eminnu = 200;
      break;
    case 2:			/* nu mode E>475 MeV */
      NEnu   = 8;
      Eminnu = 475;
      break;
    default:
      printf("chiMB_init::invalid numode: %d.\n",numode);
      exit(0);
    }

  switch(nubarmode)
    {
    case 0:			/* nubar mode off */
      NEbar   = -1;
      Eminbar = 1E10;
      break;
    case 1:			/* nubar mode E>200 MeV */
      NEbar   = 11;
      Eminbar = 200;
      break;
    case 2:			/* nubar mode E>475 MeV */
      NEbar   = 8;
      Eminbar = 475;
      break;
    default:
      printf("chiMB_init::invalid nubarmode: %d.\n",nubarmode);
      exit(0);
    }

  NCOVnu  = 2*NEnu + NMUnu;
  NCOVbar = 2*NEbar + NMUbar;   

  /* Pointers to signal, bkgd, binning */
  if(0<NEnu && NEnu<12)
    {
      datanu        = &ALLdatanu[11-NEnu];
      bgnu          = &ALLbgnu[11-NEnu];
      nosc_bgnu     = &OTHERbgnu[11-NEnu];
      nuebeam_bgnu  = &NUEbgnu[11-NEnu];
      binwnu        = &binwidths[11-NEnu];
      bincenternu   = &bincenter[11-NEnu];
    }
  
  if(0<NEbar && NEbar<12)
    {
      databar       = &ALLdatabar[11-NEbar];
      bgbar         = &ALLbgbar[11-NEbar];
      nosc_bgbar    = &OTHERbgbar[11-NEbar];
      nuebeam_bgbar = &NUEbgbar[11-NEbar];
      binwbar       = &binwidths[11-NEbar];
      bincenterbar  = &bincenter[11-NEbar];
      nuefrac       = &nuefraction[11-NEbar];
      numufrac      = &numufraction[11-NEbar];
    }

  /* Matrices to be used as covariance matrices */
  M3nu    = gsl_matrix_alloc(NCOVnu, NCOVnu);
  M2nu    = gsl_matrix_alloc(NEnu+NMUnu, NEnu+NMUnu);
  M2invnu = gsl_matrix_alloc(NEnu+NMUnu, NEnu+NMUnu);
  permnu  = gsl_permutation_alloc(NEnu+NMUnu);

  M3bar    = gsl_matrix_alloc(NCOVbar, NCOVbar);
  M2bar    = gsl_matrix_alloc(NEbar+NMUbar, NEbar+NMUbar);
  M2invbar = gsl_matrix_alloc(NEbar+NMUbar, NEbar+NMUbar);
  permbar  = gsl_permutation_alloc(NEbar+NMUbar);

  /* Nothing is taken from the glb files */
  glbDefineChiFunction(&chiMB,0,"chiMB",NULL);
  glbInitExperiment("MBneutrino200.glb",&glb_experiment_list[0],&glb_num_of_exps);
  MB_EXP = glb_num_of_exps;

/* #ifdef USE_MC_EVENTS */
  osc_params = glbAllocParams();

  FILE *f = fopen("miniboone_numunuefullosc_ntuple.txt", "r");
  Neventsnu = 0;
  while (!feof(f))
  {
    int status;
    status = fscanf(f, "%lg %lg %lg %lg", &mc_eventsnu[Neventsnu][0], &mc_eventsnu[Neventsnu][1],
                    &mc_eventsnu[Neventsnu][2], &mc_eventsnu[Neventsnu][3]);
    if (status == EOF)
      break;
    Neventsnu++;
  }
  if (f) fclose(f);

  f = fopen("miniboone_numubarnuebarfullosc_ntuple.txt", "r");
  Neventsbar = 0;
  while (!feof(f))
  {
    int status;
    status = fscanf(f, "%lg %lg %lg %lg", &mc_eventsbar[Neventsbar][0], &mc_eventsbar[Neventsbar][1],
                    &mc_eventsbar[Neventsbar][2], &mc_eventsbar[Neventsbar][3]);
    if (status == EOF)
      break;
    Neventsbar++;
  }
  if (f) fclose(f);
/* #endif */

  /* setting up L and E bins (or set up any other binning) */
  for(i=0; i<LBINS; i++)Lup[i]     = 48500 + i*(55000-48500)/LBINS;
  for(i=0; i<LBINS; i++)Lcenter[i] = 48500 + (i+0.5)*(55000-48500)/LBINS;
  for(i=0; i<EBINS; i++)Euptrue[i] = 200   + i*(3000-200)/EBINS;
  for(i=0; i<EBINS; i++)Ecenter[i] = 200   + (i+0.5)*(3000-200)/EBINS;

  /* fix units (cm, MeV to km, GeV) */
  for(i=0; i<LBINS; i++)Lup[i]     *= 1E-5;
  for(i=0; i<LBINS; i++)Lcenter[i] *= 1E-5;
  for(i=0; i<EBINS; i++)Euptrue[i] *= 1E-3;
  for(i=0; i<EBINS; i++)Ecenter[i] *= 1E-3;

  for(i=0; i<NUEBINS; i++)bincenter[i] = i==0 ? (0.2 + binwidths[i])/2. : 
			    (binwidths[i-1] + binwidths[i])/2.;

  return 0;
}
 
/***************************************************************************
 * Cleanup GSL data structures required for MiniBooNE chi^2 calculation    *
 ***************************************************************************/
int chiMB_clear()
{
#ifdef USE_MC_EVENTS
  if (osc_params) { glbFreeParams(osc_params); osc_params = NULL; }
#endif
  if (permnu)  { gsl_permutation_free(permnu); permnu = NULL; }
  if (M2invnu) { gsl_matrix_free(M2invnu);     M2nu   = NULL; }
  if (M2nu)    { gsl_matrix_free(M2nu);        M2nu   = NULL; }
  if (M3nu)    { gsl_matrix_free(M3nu);        M3nu   = NULL; }

  if (permbar)  { gsl_permutation_free(permbar); permbar = NULL; }
  if (M2invbar) { gsl_matrix_free(M2invbar);     M2bar   = NULL; }
  if (M2bar)    { gsl_matrix_free(M2bar);        M2bar   = NULL; }
  if (M3bar)    { gsl_matrix_free(M3bar);        M3bar   = NULL; }
  return 0;
}

/***************************************************************************
 * Chi^2 for the MiniBooNE NEUTRINO analysis (6.46E20 POT) described at    *
 * http://www-boone.fnal.gov/for_physicists/data_release/lowe/             *
 ***************************************************************************/
double chiMB(int exper, int rule, int n_params, double *x, double *errors,
	     void *user_data)
{
  /* The covariance matrices */
#include "cov_matrix_numode.dat"
#include "cov_matrix_nubarmode.dat"
  /* The numu and numubar spectra */
#include "barmode_numuspectrum.dat"
#include "barmode_numubarspectrum.dat"
  
  double sig[11],bg[11],chi2 = 0.0;
  double numu[MUBINS],numuosc[MUBINS];

  double (*_M3nu)[NCOVnu]        = (double (*)[NCOVnu])     gsl_matrix_ptr(M3nu,    0, 0);
  double (*_M2nu)[NEnu+NMUnu]    = (double (*)[NEnu+NMUnu]) gsl_matrix_ptr(M2nu,    0, 0);
  double (*_M2invnu)[NEnu+NMUnu] = (double (*)[NEnu+NMUnu]) gsl_matrix_ptr(M2invnu, 0, 0);
  double (*_M3bar)[NCOVbar]         = (double (*)[NCOVbar])     gsl_matrix_ptr(M3bar,    0, 0);
  double (*_M2bar)[NEbar+NMUbar]    = (double (*)[NEbar+NMUbar]) gsl_matrix_ptr(M2bar,    0, 0);
  double (*_M2invbar)[NEbar+NMUbar] = (double (*)[NEbar+NMUbar]) gsl_matrix_ptr(M2invbar, 0, 0);

  double prob_nu[EBINS][LBINS][NU_FLAVOURS][NU_FLAVOURS];
  double prob_bar[EBINS][LBINS][NU_FLAVOURS][NU_FLAVOURS];
  double prob_temp[NU_FLAVOURS][NU_FLAVOURS];
  double p = 0, p2=0, p3=0;
  double *length,*density,*pp;
  double dens[1] = {0}; density = &dens[0];
  double filter_sigma = 5.5E-2;
  int i, j, k;

  /* Generating NU probability table in L and E bins */
  /* This will be used in the NUBAR mode */
  for(i=0; i<EBINS; i++)
    for(j=0; j<LBINS; j++)
      {
	length = &Lcenter[j];
	pp = &prob_temp;
	glb_probability_matrix(pp, 1, Ecenter[i],
			       1, length, density, filter_sigma, NULL);
	
	for(int A=0; A<NU_FLAVOURS; A++)
	  for(int B=0; B<NU_FLAVOURS; B++)
	    prob_nu[i][j][A][B] = prob_temp[A][B];
      }

  for(i=0; i<MUBINS; i++)
    {numu[i] = 0; numuosc[i] = 0;}

  /************* NEUTRINO MODE *************/
  if(NEnu > 7)
    {
      /* Compute signal using Monte Carlo events */
      /* #ifdef USE_MC_EVENTS */
      glbGetOscillationParameters(osc_params);
      double dmsq  = glbGetOscParams(osc_params, GLB_DM_21);
      double s22th = SQR(sin(2.0 * glbGetOscParams(osc_params, GLB_THETA_12)));
      for (int i=0; i < NEnu; i++)
	sig[i] = 0.0;
      
      for (long k=0; k < Neventsnu; k++)
	{
	  double Ereco, Etrue, L, w;
	  Ereco = mc_eventsnu[k][0];
	  Etrue = 1E-3*mc_eventsnu[k][1];
	  L     = 1E-5*mc_eventsnu[k][2];
	  w     = mc_eventsnu[k][3];
	  
	  /* Localize event in Etrue -- L grid */
	  i = 0; j = 0;
	  while (Euptrue[i] <= Etrue && i <= EBINS-2) i++;
	  while (Lup[j] <= L && j <= LBINS-2) j++;
	  p = prob_nu[i][j][1][0];  /* mu -> e */
	  p2 = prob_nu[i][j][1][1]; /* mu -> mu */
	  double Eup = Eminnu + 1000.*binwnu[0];
	  i = 0;
	  while (Ereco >= Eup && i < NEnu)
	    Eup += 1000.*binwnu[++i];
	  if (i < NEnu && Ereco >= Eminnu)
	    sig[i] += p * w;

	  /* NUMU disappearance */
	  Eup = numubin[0];
	  i = 0;
	  while (Ereco >= Eup && i < MUBINS)
	    Eup = numubin[++i];
	  if (i < MUBINS && Ereco >= 0)
	    { numu[i] += w;
	      numuosc[i] += p2 * w; }
	}

      for(i=0; i < NEnu; i++)
	sig[i] = sig[i]/Neventsnu;
      for(i=0; i < MUBINS; i++)
	numuosc[i] = numuosc[i]/numu[i];

/* #endif */
    
      /* beam nue background oscillation */
      for(int k=0; k<NEnu; k++) 
	{
	  /* Localize event in Etrue -- L grid */
	  i = 0; j = 0;
	  while (Euptrue[i] <= bincenternu[k] && i <= EBINS-2) i++;
	  while (Lup[j] <= 0.520 && j <= LBINS-2) j++;
	  bg[k] = nosc_bgnu[k] + prob_nu[i][j][0][0]*nuebeam_bgnu[k]; 
	}

  /* Construct covariance matrix in 3x3 block form */
  double P[NCOVnu];         /* Vector of pred nu_e signal, nu_e bg, nu_mu signal */
  for (i=0; i < NEnu; i++)
  {
    P[i]      = sig[i];
    P[i+NEnu] = bg[i];
  }
  for (i=0; i < NMUnu; i++)
    P[i+2*NEnu] = pred_numu[i]*numuosc[i];

  for (i=0; i < NCOVnu; i++)
    for (j=0; j < NCOVnu; j++)
      {
	_M3nu[i][j] = NEnu==11 ? Mfrac_nu_200[i][j] * P[i] * P[j] :
	                       Mfrac_nu_475[i][j] * P[i] * P[j];
      }
  for (i=0; i < NEnu; i++)        /* Add statistical error of signal */
    _M3nu[i][i] += sig[i];

  /* Collapse covariance matrix to 2x2 block form */
  for (i=0; i < NEnu; i++)        /* Upper left block */
    for (j=0; j < NEnu; j++)
      _M2nu[i][j] = _M3nu[i][j] + _M3nu[i+NEnu][j] + _M3nu[i][j+NEnu] + _M3nu[i+NEnu][j+NEnu];
  for (i=0; i < NMUnu; i++)       /* Lower left block */
    for (j=0; j < NEnu; j++)
      _M2nu[i+NEnu][j] = _M3nu[i+2*NEnu][j] + _M3nu[i+2*NEnu][j+NEnu];
  for (i=0; i < NEnu; i++)       /* Upper right block */
    for (j=0; j < NMUnu; j++)
      _M2nu[i][j+NEnu] = _M3nu[i][j+2*NEnu] + _M3nu[i+NEnu][j+2*NEnu];
  for (i=0; i < NMUnu; i++)      /* Lower right block */
    for (j=0; j < NMUnu; j++)
      _M2nu[i+NEnu][j+NEnu] = _M3nu[i+2*NEnu][j+2*NEnu];

  /* Invert covariance matrix and compute log-likelihood */
  int signum;
  gsl_linalg_LU_decomp(M2nu, permnu, &signum);
  gsl_linalg_LU_invert(M2nu, permnu, M2invnu);

  double P2[NEnu+NMUnu];
  for (i=0; i < NEnu; i++)
    P2[i] = datanu[i] - (sig[i] + bg[i]);
  for (i=0; i < NMUnu; i++)
    P2[i+NEnu] = datanu[i+NEnu] - pred_numu[i]*numuosc[i];
  for (i=0; i < NEnu+NMUnu; i++)
    for (j=0; j < NEnu+NMUnu; j++)
      chi2 += P2[i] * _M2invnu[i][j] * P2[j];
  /* chi2 += gsl_linalg_LU_lndet(M2); */
    }  

  for(i=0; i<MUBINS; i++)
    {numu[i] = 0; numuosc[i] = 0;}

  double scale2012=11.3/5.66,scale2=1.015*11.3/5.66;
  /************* ANTINEUTRINO MODE *************/
  if(NEbar > 7)
    {
      /* Compute signal using Monte Carlo events */
      /* #ifdef USE_MC_EVENTS */
      glbGetOscillationParameters(osc_params);
      double dmsq  = glbGetOscParams(osc_params, GLB_DM_21);
      double s22th = SQR(sin(2.0 * glbGetOscParams(osc_params, GLB_THETA_12)));
      for (int i=0; i < NEbar; i++)
	sig[i] = 0.0;
      /* Generating probability table in L and E bins */
      for(i=0; i<EBINS; i++)
	for(j=0; j<LBINS; j++)
	  {
	    length = &Lcenter[j];
	    pp = &prob_temp;
	    glb_probability_matrix(pp, -1, Ecenter[i],
	    			   1, length, density, filter_sigma, NULL);

	    for(int A=0; A<NU_FLAVOURS; A++)
 	      for(int B=0; B<NU_FLAVOURS; B++)
	    	prob_bar[i][j][A][B] = prob_temp[A][B];
	  }
      
      for (long k=0; k < Neventsbar; k++)
	{
	  double Ereco, Etrue, L, w;
	  Ereco = mc_eventsbar[k][0];
	  Etrue = 1E-3*mc_eventsbar[k][1];
	  L     = 1E-5*mc_eventsbar[k][2];
	  w     = mc_eventsbar[k][3];
	  
	  /* Localize event in Etrue -- L grid */
	  i = 0; j = 0;
	  while (Euptrue[i] <= Etrue && i <= EBINS-2) i++;
	  while (Lup[j] <= L && j <= LBINS-2) j++;
	  p = prob_bar[i][j][1][0];  /* mub -> eb  */
	  p2 = prob_nu[i][j][1][0];  /* mu  -> e   */
	  p3 = prob_bar[i][j][1][1]; /* mub -> mub */
	  double Eup = Eminbar + 1000.*binwbar[0];
	  i = 0;
	  while (Ereco >= Eup && i < NEbar)
	    Eup += 1000.*binwbar[++i];
	  if (i < NEbar && Ereco >= Eminbar)
	    sig[i] += (p + p2*numufrac[i]) * w;
	  /* Here, the numu fraction in the nubar mode is assumed to oscillate! */

	  /* NUMUBAR disappearance - I am neglecting WS nu in numu to numu sample */
	  Eup = numubin[0];
	  i = 0;
	  while (Ereco >= Eup && i < MUBINS)
	    Eup = numubin[++i];
	  if (i < MUBINS && Ereco >= 0)
	    { 
	      numu[i] += w;
	      numuosc[i] += p3 * w; }

	}

      for (i=0; i < NEbar; i++)
	sig[i] = sig[i]*scale2012/Neventsbar;
      for(i=0; i < MUBINS; i++)
	numuosc[i] = numuosc[i]/numu[i];

/* #endif */
    
      /* beam nue background oscillation */
      for(int k=0; k<NEbar; k++) 
	{
	  /* Localize event in Etrue -- L grid */
	  i = 0; j = 0;
	  while (Euptrue[i] <= bincenterbar[k] && i <= EBINS-2) i++;
	  while (Lup[j] <= 0.520 && j <= LBINS-2) j++;
	  bg[k] = nosc_bgbar[k] + ((1.0 - nuefrac[k])*prob_bar[i][j][0][0]
	     + nuefrac[k]*prob_nu[i][j][0][0])*nuebeam_bgbar[k];
	}

  /* Construct covariance matrix in 3x3 block form */
  double P[NCOVbar];         /* Vector of pred bar_e signal, nu_e bg, nu_mu signal */
  for (i=0; i < NEbar; i++)
  {
    P[i]      = sig[i];
    P[i+NEbar] = bg[i]*scale2;
  }
  for (i=0; i < NMUbar; i++)
    P[i+2*NEbar] = pred_numubar[i]*scale2012;

  for (i=0; i < NCOVbar; i++)
    for (j=0; j < NCOVbar; j++)
      {
	_M3bar[i][j] = NEbar==11 ? Mfrac_bar_200[i][j] * P[i] * P[j] :
	                       Mfrac_bar_475[i][j] * P[i] * P[j];
      }
  for (i=0; i < NEbar; i++)        /* Add statistical error of signal */
    _M3bar[i][i] += sig[i];

  double new_bg_cov_diag[11] = { 0.0352733, 0.0370013, 0.0408353, 0.0501032,
  				 0.0524085, 0.0547091, 0.0551101, 0.0662563,
  				 0.0701908, 0.108929,  0.0776587  };
  for(i=NEbar; i < 2*NEbar; i++)      /* correcting the bg statistical error */
      _M3bar[i][i] = new_bg_cov_diag[i-NEbar] * P[i] * P[i];


  /* Collapse covariance matrix to 2x2 block form */
  for (i=0; i < NEbar; i++)        /* Upper left block */
    for (j=0; j < NEbar; j++)
      _M2bar[i][j] = _M3bar[i][j] + _M3bar[i+NEbar][j] + _M3bar[i][j+NEbar] + _M3bar[i+NEbar][j+NEbar];
  for (i=0; i < NMUbar; i++)       /* Lower left block */
    for (j=0; j < NEbar; j++)
      _M2bar[i+NEbar][j] = _M3bar[i+2*NEbar][j] + _M3bar[i+2*NEbar][j+NEbar];
  for (i=0; i < NEbar; i++)       /* Upper right block */
    for (j=0; j < NMUbar; j++)
      _M2bar[i][j+NEbar] = _M3bar[i][j+2*NEbar] + _M3bar[i+NEbar][j+2*NEbar];
  for (i=0; i < NMUbar; i++)      /* Lower right block */
    for (j=0; j < NMUbar; j++)
      _M2bar[i+NEbar][j+NEbar] = _M3bar[i+2*NEbar][j+2*NEbar];

  /* Invert covariance matrix and compute log-likelihood */
  int signum;
  gsl_linalg_LU_decomp(M2bar, permbar, &signum);
  gsl_linalg_LU_invert(M2bar, permbar, M2invbar);

  double P2[NEbar+NMUbar];
  for (i=0; i < NEbar; i++)
    P2[i] = databar[i] - (sig[i] + bg[i]*scale2);
  for (i=0; i < NMUbar; i++)
    P2[i+NEbar] = databar[i+NEbar] - pred_numubar[i]*scale2012*numuosc[i];
  for (i=0; i < NEbar+NMUbar; i++)
    for (j=0; j < NEbar+NMUbar; j++)
      chi2 += P2[i] * _M2invbar[i][j] * P2[j];
  /* chi2 += gsl_linalg_LU_lndet(M2bar); */
    }  



  return chi2;
}

