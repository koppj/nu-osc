/***************************************************************************
 * Functions for MiniBooNE numu to nue fit                                 *
 ***************************************************************************
 * Author: PAN Machado, 2012                                               *
 * See http://www-boone.fnal.gov/for_physicists/data_release/lowe/         *
 ***************************************************************************/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <globes/globes.h>
#include "nu.h"
#include "snu.h"
#include "glb_wrapper.h"
#include "glb_error.h"

// Size of probability arrays - we care only about the 3 active flavours here
#define NU_FLAVOURS 3

// Path for external data files. Note: #include's need to be adjusted separately
#define MB_DATA_PATH "./glb/mb-2012/mb-data/"

// TAG_DEFS - This signals ./compile-and-run where to insert new #define's
#define NU_APP               // TAG_DEF - this signals ./compile-and-run
#define NU_DISAPP            // TAG_DEF   which #define's to replace
#define NUBAR_APP            // TAG_DEF
#define NUBAR_DISAPP         // TAG_DEF
#define COMBINED_APP         // TAG_DEF
#define FULL_OSC             // TAG_DEF

#define DIS_RANK 42

/* External variables */
extern int n_flavors;

/* Neutrino mode parameters */
int NE    = 0;             /* Number of E_nue bins      */
int NMU   = 8;             /* Number of E_numue bins    */
int NCOV  = 0;             /* Size of covariance matrix */

int APP_FIRST = 0;
int APP_LAST  = 11;

int MB_EXP  = -1;

#define MAX_MC_EVT 1267007
/* #define NN 40 */

/************************ MINIBOONE DATA ************************/
/* Common to both modes */
//const double numubin[8] = { 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9 }; // JK - NOT USED
const double disbin[16] = { 0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  
                            1.2,  1.3,  1.4,  1.5,  1.6,  1.7,  1.8,  1.9 };

/* NU MODE: */
double ALLdatanu[19]   = { 232,  156, 156, 79, 81, 70, 63, 65, 62, 34,  70,
                           19189, 29943, 26748, 18692, 11123, 5830, 2674, 1273 };
double pred_numu[8]    = { 19402.247,   29854.561,   26700.085,   18701.467,
                           11572.734,    6250.383,   3025.1541,   1409.3744 };
double ALLbgnu[11]     = { 180.80171,   108.22448,   120.03353,   63.887782,   89.806966,
                           67.249431,   69.855878,   57.014477,   51.846417,   38.586738,   69.381391 };
/* separate beam nue background from others: */
double OTHERbgnu[11]  = { 163.326, 86.1045, 80.638, 34.8518, 42.7113, 
                          21.355,  22.551,  16.401, 13.668,  9.568,   20.505};
double NUEbgnuMU[11]  = { 13.3757, 15.9692, 27.7785, 18.2728, 30.012, 
                          29.6644, 28.8534, 23.188, 19.0444, 12.6187, 18.1264 };
double NUEbgnuK[11]    = { 4.1,   6.15075, 11.617,  10.7632, 17.0838, 
                           16.23, 18.4515, 17.4255, 19.134,  16.4, 30.75};
double nudisdata[16]  = {  8852, 16237, 22150, 24367, 23827, 
                          21895, 19142, 14996, 11673,  8963, 
                           6348,  4669,  3094,  2070,  1286,  885 };

/* NUBAR MODE: */
double ALLdatabar[19]  = { 122, 70, 65, 43, 57, 39, 37, 23, 26, 26, 43,
                           9481, 13581, 11308, 7667, 4682, 2371, 985, 380 };
double pred_numubar[8] = { 9998.9575,   13461.302,    11298.24,   7604.9604,
                           4331.8869,   2125.5371,   891.22261,   336.98711 };
double ALLbgbar[11]    = { 90.289907,   53.077595,   57.098801,   32.937945,   43.159072,
                           34.174322,   36.383542,   28.737807,   28.750297,   20.098079,   42.697791 };
double OTHERbgbar[11]  = { 75.73, 38.3175, 33.617, 14.913, 14.6988, 
                           8.4325, 7.2705, 4.8855, 4.57, 1.998, 14.19 };
double NUEbgbarMU[11]  = { 8.105, 7.71375, 11.844, 8.649, 14.805, 
                           12.0775, 12.6225, 10.752, 9.974, 6.856, 4.71};
double NUEbgbarK[11]  = { 6.233, 7.47975, 11.843, 9.58425, 14.805, 
                          14.0262, 15.897, 13.0935, 13.716, 10.596, 14.085 };
double nubardisdata[16] = { 1565, 2476, 3121, 3203, 3117,
                            2919, 2496, 2006, 1761, 1322,
                            997,   750,  505,  377,  257, 181 };
/* fraction of nue/nue+nuebar of beam bg */
double nuefraction[11];
/* fraction of numu/numu+numubar (of numu oscillates, evts*numufraction*prob should be added to the signal) */
double numufraction[11];
/****************************************************************/

/********************* MB-SB DIS VARIABLES **********************/
double bins_MBSB[21][2], binw_MBSB[21][2], bincenter_MBSB[21][2];
double MBSB_data[21][2]; /* 0-MB, 1-SB */
double MBSB_pred[21][2][2];     /* 0-MB, 1-SB, 0-RS, 2-WS */
double MfracMBSB[21][21][2][2][2][2]; /* MB-SB, MB-SB, RS-WS, RS-WS */
double totalerrmatrix[DIS_RANK][DIS_RANK];

static gsl_matrix *M_MBSB        = NULL;
static gsl_matrix *M_MBSB_inv     = NULL;
static gsl_permutation *perm_MBSB = NULL;
/****************************************************************/


/* Bin widths for appearance analysis, see
 * http://www-boone.fnal.gov/for_physicists/data_release/nue_nuebar_2012/ */
double binwidths[11] = { 0.1, 0.075, 0.1, 0.075, 0.125, 0.125, 0.15, 0.15, 0.2, 0.2, 1.5 };
double bincenter[11];

/* for binning the MC data */
//#define LBINS (10)
//#define EBINS (60)
#define LBINS (7)
#define EBINS (40)
double Emin,Euptrue[EBINS],Ecenter[EBINS],Lup[LBINS],Lcenter[LBINS],
       Lup_SB[LBINS],Lcenter_SB[LBINS];

/* spectrum variables - signal, other, nue from K, nue from mu */
/* call the chi2 and the spectrum will be available */
double app_nu_spectrum[11][4],app_bar_spectrum[11][4],
  dis_nu_spectrum[16],dis_bar_spectrum[16];

static gsl_matrix *M3nu        = NULL;
static gsl_matrix *M2nu        = NULL;
static gsl_matrix *M2invnu     = NULL;
static gsl_permutation *permnu = NULL;

static gsl_matrix *M3bar        = NULL;
static gsl_matrix *M2bar        = NULL;
static gsl_matrix *M2invbar     = NULL;
static gsl_permutation *permbar = NULL;

static gsl_matrix *Mdis         = NULL;
static gsl_matrix *Minvdis      = NULL;
static gsl_permutation *permdis = NULL;

static gsl_matrix *Mdisbar         = NULL;
static gsl_matrix *Minvdisbar      = NULL;
static gsl_permutation *permdisbar = NULL;

static gsl_matrix *M3        = NULL;
static gsl_matrix *M2        = NULL;
static gsl_matrix *M2inv     = NULL;
static gsl_permutation *perm = NULL;

static long Neventsnu;  /* app neutrino MC events */
static double mc_eventsnu[128000][4];
static long Neventsbar; /* app antineutrino MC events */
static double mc_eventsbar[128000][4];
static long long Neventsnumu;    /* disapp neutrino MC events */
static double mc_eventsnumu[1267007][4];
static long long Neventsnumubar; /* disapp antineutrino MC events */
static double mc_eventsnumubar[700000][4];
static long long Nevents_MB; /* disapp MB-SB analysis MC events */
static long long Nevents_SB; /* disapp MB-SB analysis MC events */
static double mc_events_MB[2000000][4];
static double mc_events_SB[700000][4];
static int mc_events_MB_type[2000000];
static int mc_events_SB_type[700000];

  /* Probability arrays */
  double prob_nu[EBINS][LBINS][NU_FLAVOURS][NU_FLAVOURS];
  double prob_bar[EBINS][LBINS][NU_FLAVOURS][NU_FLAVOURS];
  double prob_nu_SB2[EBINS][LBINS][NU_FLAVOURS][NU_FLAVOURS];
  double prob_bar_SB2[EBINS][LBINS][NU_FLAVOURS][NU_FLAVOURS];
  double prob_nu_SB[11][NU_FLAVOURS][NU_FLAVOURS];
  double prob_bar_SB[11][NU_FLAVOURS][NU_FLAVOURS];

double chiMB(int exper, int rule, int n_params, double *x, double *errors,
             void *user_data);

/***************************************************************************
 * Compute beam contaminations from MC data                                *
 ***************************************************************************/
int contaminations()
{
  int i;
  double numu[11],numubar[11],nue[11],nuebar[11];
  for(i=0; i<11; i++)
    {
      numu[i]    = 0;
      numubar[i] = 0;
      nue[i]     = 0;
      nuebar[i]  = 0;
    }

  /* MB antineutrino disappearance MC events: */
  for (long k=0; k < Nevents_MB; k++)
    {
      int type; double Ereco, Etrue, L, w;
      type  = mc_events_MB_type[k];
      Etrue = mc_events_MB[k][0]; /* GeV */
      Ereco = mc_events_MB[k][1]; /* GeV */
      L     = 1E-3*mc_events_MB[k][2];
      w     = mc_events_MB[k][3];

      double Eup = 0.2 + binwidths[0];
      i = 0;
      /* if(Ereco>1.9&&Etrue>1.9)printf("%f %f\n",Ereco,Etrue); */
      while (Ereco >= Eup && i < 11)
        Eup += binwidths[++i];

      if (i < 11 && Ereco >= 0.2)
        {
      if(type==14)  numu[i]    += w;  /* numu     */
      if(type==-14) numubar[i] += w;  /* numubar  */
      if(type==12)  nue[i]     += w;  /* nue      */
      if(type==-12) nuebar[i]  += w;  /* nuebar   */
        }
    }

  /* the fractions */
  for(i=0; i<11; i++)
    {
      numufraction[i] = numu[i]+numubar[i]==0 ? -1 : numu[i]/(numu[i]+numubar[i]);
      nuefraction[i]  = nue[i]+nuebar[i]==0   ? -1 : nue[i]/(nue[i]+nuebar[i]);
    }
  if(numufraction[0]<0) numufraction[0] = numufraction[1];
  if(nuefraction[0]<0)  nuefraction[0]  = nuefraction[1];
  /* for(i=0; i<11; i++) printf("fractions (numu,nue) %d: %f %f\n",i,numufraction[i],nuefraction[i]); */

  return 0;
}


/***************************************************************************
 * Retrieve event spectra                                                  *
 ***************************************************************************/
int getspectrum(const char *fname)
{
  FILE *file= fopen(fname, "w");

#if defined(NU_APP) || defined(COMBINED_APP)
  for(int i=0; i<11; i++) 
    fprintf(file,"%f  %f  %f  %f\n",app_nu_spectrum[i][0],
	    app_nu_spectrum[i][1],app_nu_spectrum[i][2],app_nu_spectrum[i][3]);
  fprintf(file,"\n");
#endif

#if defined(NUBAR_APP) || defined(COMBINED_APP)
  for(int i=0; i<11; i++) 
    fprintf(file,"%f  %f  %f  %f\n",app_bar_spectrum[i][0],app_bar_spectrum[i][1],
	    app_bar_spectrum[i][2],app_bar_spectrum[i][3]);
#endif

  return 0;
}


/***************************************************************************
 * Initialize GSL data structures required for MiniBooNE chi^2 calculation *
 * threshold = 0,1:                                                        *
 * 0 - Enu > 200 MeV                                                       *
 * 1 - Enu > 475 MeV                                                       *
 ***************************************************************************/
int chiMB_init(int threshold)
{
  int i;
  int status;

  printf("# Initializing MiniBooNE code, Sep 2012 version ...\n");
  printf("#   Energy threshold: %d MeV\n", threshold==0 ? 200 : 475);
  printf("#   MiniBooNE flags: ");
#ifdef NU_APP
  printf("NU_APP ");
#endif
#ifdef NU_DISAPP
  printf("NU_DISAPP ");
#endif
#ifdef COMBINED_APP
  printf("COMBINED_APP ");
#endif
#ifdef NUBAR_APP
  printf("NUBAR_APP ");
#endif
#ifdef NUBAR_DISAPP
  printf("NUBAR_DISAPP ");
#endif
#ifdef FULL_OSC
  printf("FULL_OSC ");
#endif
  printf("\n");

  switch(threshold)
    {
    case 0:                     /* 200 MeV */
      NE   = 11;
      Emin = 200;
      APP_FIRST = 0;
      APP_LAST = 11;
      break;
    case 1:                     /* 475 MeV */
      NE   = 8;
      Emin = 475;
      APP_FIRST = 3;
      APP_LAST = 11;
      break;
    default:
      printf("chiMB_init::invalid threshold mode: %d.\n",threshold);
      exit(0);
    }

  NCOV  = 2*NE + NMU;

  /* Matrices to be used as covariance matrices */
  M3nu    = gsl_matrix_alloc(NCOV, NCOV);
  M2nu    = gsl_matrix_alloc(NE+NMU, NE+NMU);
  M2invnu = gsl_matrix_alloc(NE+NMU, NE+NMU);
  permnu  = gsl_permutation_alloc(NE+NMU);

  M3bar    = gsl_matrix_alloc(NCOV, NCOV);
  M2bar    = gsl_matrix_alloc(NE+NMU, NE+NMU);
  M2invbar = gsl_matrix_alloc(NE+NMU, NE+NMU);
  permbar  = gsl_permutation_alloc(NE+NMU);

  Mdis     = gsl_matrix_alloc(16, 16);
  Minvdis  = gsl_matrix_alloc(16, 16);
  permdis  = gsl_permutation_alloc(16);

  Mdisbar     = gsl_matrix_alloc(16, 16);
  Minvdisbar  = gsl_matrix_alloc(16, 16);
  permdisbar  = gsl_permutation_alloc(16);

  /* matrices for combined the analysis */
  M3    = gsl_matrix_alloc(2*NCOV, 2*NCOV);
  M2    = gsl_matrix_alloc(2*(NE+NMU), 2*(NE+NMU));
  M2inv = gsl_matrix_alloc(2*(NE+NMU), 2*(NE+NMU));
  perm  = gsl_permutation_alloc(2*(NE+NMU));
  /* matrices for MB-SB numubar disapp analysis */
  M_MBSB     = gsl_matrix_alloc(DIS_RANK, DIS_RANK);
  M_MBSB_inv = gsl_matrix_alloc(DIS_RANK, DIS_RANK);
  perm_MBSB  = gsl_permutation_alloc(DIS_RANK);

  /* Nothing is taken from the glb files */
  glbDefineChiFunction(&chiMB,0,"chiMB",NULL);
  glbInitExperiment("MBneutrino200.glb",&glb_experiment_list[0],&glb_num_of_exps);
  MB_EXP = glb_num_of_exps;


  /* The following 2 data files are from the MiniBooNE combined \nu_e and \bar{nu}_e
   * appearance analysis, see
   * http://www-boone.fnal.gov/for_physicists/data_release/nue_nuebar_2012/combined.html */

  /* FILE *f = fopen(MB_DATA_PATH"/miniboone_numunuefullosc_ntuple.txt", "r"); /\* app neutrinos *\/ */
  FILE *f = fopen(MB_DATA_PATH"/miniboone_nufullosc_ntuple.txt", "r"); /* app neutrinos */
  Neventsnu = 0;
  while (!feof(f))
  {
    status = fscanf(f, "%lg %lg %lg %lg", &mc_eventsnu[Neventsnu][0], &mc_eventsnu[Neventsnu][1],
                    &mc_eventsnu[Neventsnu][2], &mc_eventsnu[Neventsnu][3]);
    if (status == EOF) //|| Neventsnu > MAX_MC_EVT)
      break;
    Neventsnu++;
  }
  if (f) fclose(f);

  f = fopen(MB_DATA_PATH"/miniboone_nubarfullosc_ntuple.txt", "r"); /* app antinu */
  Neventsbar = 0;
  while (!feof(f))
  {
    status = fscanf(f, "%lg %lg %lg %lg", &mc_eventsbar[Neventsbar][0], &mc_eventsbar[Neventsbar][1],
                    &mc_eventsbar[Neventsbar][2], &mc_eventsbar[Neventsbar][3]);
    if (status == EOF || Neventsbar > MAX_MC_EVT)
      break;
    Neventsbar++;
  }
  if (f) fclose(f);


  /* The next 2 data files are from the 2009 MiniBooNE \nu_\mu and \bar\nu_\mu
   * disappearance analysis, see
   * http://www-boone.fnal.gov/for_physicists/data_release/numu_numubar/ 
   * We use this only for \nu_\mu disappearance - in \bar\nu_\mu mode, we use
   * the combined MB/SB analysis instead */

  double dummy;
  f = fopen(MB_DATA_PATH"/numudisap_ntuple.txt", "r"); /* disapp neutrino */
  Neventsnumu = 0;
  while (!feof(f))
  {
    status = fscanf(f, "%lg %lg %lg %lg %lg", &dummy,&mc_eventsnumu[Neventsnumu][0], &mc_eventsnumu[Neventsnumu][1],
                    &mc_eventsnumu[Neventsnumu][2], &mc_eventsnumu[Neventsnumu][3]);
    if (status == EOF || Neventsnumu > MAX_MC_EVT )
      break;
    Neventsnumu++;
  }
  if (f) fclose(f);

  f = fopen(MB_DATA_PATH"/numubardisap_ntuple.txt", "r"); /* disapp antineutrino, not used */
  Neventsnumubar = 0;
  while (!feof(f))
  {
    status = fscanf(f, "%lg %lg %lg %lg %lg", &dummy,&mc_eventsnumubar[Neventsnumubar][0], &mc_eventsnumubar[Neventsnumubar][1],
                    &mc_eventsnumubar[Neventsnumubar][2], &mc_eventsnumubar[Neventsnumubar][3]);
    if (status == EOF || Neventsnumubar > MAX_MC_EVT )
      break;
    Neventsnumubar++;
  }
  if (f) fclose(f);


  /******* reading MB-SB files *******/
  /* BINS */
  /* for MB/SB combined disappearance analysis, see
   * http://www-sciboone.fnal.gov/data_release/joint_numubar_disap/ */
  f = fopen(MB_DATA_PATH"/bins.txt", "r"); i=0;
  /* while (!feof(f)) */
  for(i=0; i<21; i++)
  { status = fscanf(f, "%lg %lg", &bins_MBSB[i][0],&bins_MBSB[i][1]);
    bins_MBSB[i][0] = bins_MBSB[i][1];
    for(int j=0; j<2; j++)
      binw_MBSB[i][j] = i==0 ? bins_MBSB[i][j]-0.3 : bins_MBSB[i][j] -bins_MBSB[i-1][j];
    /* i++; */
    if(status == EOF || i==21) break; }
  if (f) fclose(f);

  /* for(i=0; i<21; i++) */
  /*   printf("bin: %d  %f  %f\n",i,bins_MBSB[i][0],bins_MBSB[i][1]); */

  /* DATA */
  /* for MB/SB combined disappearance analysis, see
   * http://www-sciboone.fnal.gov/data_release/joint_numubar_disap/ */
  for(int j=0; j<2; j++)
    {
      if(j==0) f = fopen(MB_DATA_PATH"/mb_data.txt", "r"); i=0;
      if(j==1) f = fopen(MB_DATA_PATH"/sb_data.txt", "r"); i=0;
      /* while (!feof(f)) */
      for(i=0; i<21; i++)
        { status = fscanf(f, "%lg", &MBSB_data[i][j]);
          /* printf("%d %d %g\n",i,j,MBSB_data[i][j]); */
          /* i++; */
          if(status == EOF) break; }
      if (f) fclose(f);
    }

  /* for(i=0; i<21; i++) */
  /*   printf("data: %d   MB: %f   SB: %f\n",i,MBSB_data[i][0],MBSB_data[i][1]); */

  /* PREDICTED EVENTS */
  /* for MB/SB combined disappearance analysis, see
   * http://www-sciboone.fnal.gov/data_release/joint_numubar_disap/ */
  for(int j=0; j<2; j++)
    for(int jj=0; jj<2; jj++)
    {
      if(j==0 && jj==0) f = fopen(MB_DATA_PATH"/mb_rs.txt", "r"); i=0;
      if(j==0 && jj==1) f = fopen(MB_DATA_PATH"/mb_ws.txt", "r"); i=0;
      if(j==1 && jj==0) f = fopen(MB_DATA_PATH"/sb_rs.txt", "r"); i=0;
      if(j==1 && jj==1) f = fopen(MB_DATA_PATH"/sb_ws.txt", "r"); i=0;
      /* while (!feof(f)) */
  for(i=0; i<21; i++)
        { status = fscanf(f, "%lg", &MBSB_pred[i][j][jj]);
          /* i++; */
          if(status == EOF || i==21) break; }
      if (f) fclose(f);
    }

  /* for(i=0; i<21; i++) */
  /*   printf("pred RS: %d   MB: %f   SB: %f\n",i,MBSB_pred[i][0][0],MBSB_pred[i][1][0]); */
  /* for(i=0; i<21; i++) */
  /*   printf("pred WS: %d   MB: %f   SB: %f\n",i,MBSB_pred[i][0][1],MBSB_pred[i][1][1]); */

  /* FRACTIONAL COVARIANCE MATRICES */
  /* for MB/SB combined disappearance analysis, see
   * http://www-sciboone.fnal.gov/data_release/joint_numubar_disap/ */
  for(int j=0; j<2; j++)
    for(int k=0; k<2; k++)
      for(int l=0; l<2; l++)
        for(int m=0; m<2; m++)
          {
            char filename[50];
            strcpy(filename,MB_DATA_PATH"/");
            if(j==0) strcat(filename,"mb");
            if(j==1) strcat(filename,"sb");
            if(k==0) strcat(filename,"mb");
            if(k==1) strcat(filename,"sb");
            if(l==0) strcat(filename,"_rs");
            if(l==1) strcat(filename,"_ws");
            if(m==0) strcat(filename,"rs_frac.txt");
            if(m==1) strcat(filename,"ws_frac.txt");
            f = fopen(filename, "r"); i=0;
            /* printf("%s\n",filename); */
            /* while (!feof(f)) */
            for(i=0; i<21; i++)
              { for(int ii=0; ii<21; ii++) 
                  status = fscanf(f, "%lg", &MfracMBSB[i][ii][j][k][l][m]);
                /* i++; */
                if(status == EOF || i==21) break; }
            if (f) fclose(f);
          }

  /* for(i=0; i<21; i++) */
  /*   { */
  /*     for(int ii=0; ii<21; ii++) */
  /*    printf("%f   ",MfracMBSB[i][ii][1][1][1][1]); */
  /*     printf("\n   "); */
  /*   } */

  /* /\* TOTAL COVARIANCE *\/ */
  /* for(int j=0; j<2; j++) */
  /*   for(int k=0; k<2; k++) */
  /*     for(int l=0; l<2; l++) */
  /*    for(int m=0; m<2; m++) */
  /*      { */
  /*        f = fopen(MB_DATA_PATH"/total_error_matrix.txt", "r"); i=0; */
  /*        while (!feof(f)) */
  /*          { for(int ii=0; ii<21; ii++)  */
  /*              status = fscanf(f, "%lg", &MfracMBSB[i][ii][j][k][l][m]); */
  /*            i++; */
  /*            if(status == EOF || i==21) break; } */
  /*        if (f) fclose(f); */
  /*      } */
  
  /* MC events */
  f = fopen(MB_DATA_PATH"/MB_mc_events.txt", "r"); /* disapp antineutrino */
  Nevents_MB = 0;
  while (!feof(f))
  {
    status = fscanf(f, "%d %lg %lg %lg %lg", &mc_events_MB_type[Nevents_MB], 
                    &mc_events_MB[Nevents_MB][0], &mc_events_MB[Nevents_MB][1], 
                    &mc_events_MB[Nevents_MB][2], &mc_events_MB[Nevents_MB][3]);
    if (status == EOF)
      break;
    Nevents_MB++;
  }
  if (f) fclose(f);

  f = fopen(MB_DATA_PATH"/SB_mc_events.txt", "r"); /* disapp antineutrino */
  Nevents_SB = 0;
  while (!feof(f))
  {
    status = fscanf(f, "%d %lg %lg %lg %lg", &mc_events_SB_type[Nevents_SB],
                    &mc_events_SB[Nevents_SB][0], &mc_events_SB[Nevents_SB][1],
                    &mc_events_SB[Nevents_SB][2], &mc_events_SB[Nevents_SB][3]);
    if (status == EOF)
      break;
    Nevents_SB++;
  }
  if (f) fclose(f);

  /* printf("%d %d\n",Nevents_SB,Nevents_MB);  exit(0); */

  /* setting up L and E bins (or set up any other binning) */
  /* for(i=0; i<LBINS; i++)Lup[i]     = 48500 + i*(55000-48500)/LBINS; */
  /* for(i=0; i<LBINS; i++)Lcenter[i] = 48500 + (i-0.5)*(55000-48500)/LBINS; */
  /* for(i=0; i<LBINS; i++)Lup_SB[i]     = 40 + i*(110-40)/LBINS; */
  /* for(i=0; i<LBINS; i++)Lcenter_SB[i] = 40 + (i-0.5)*(110-40)/LBINS; */
  /* for(i=0; i<EBINS; i++)Euptrue[i] = 200   + i*(3000-200)/EBINS; */
  /* for(i=0; i<EBINS; i++)Ecenter[i] = 200   + (i-0.5)*(3000-200)/EBINS; */
  for(i=0; i<LBINS; i++)Lup[i]     = 48550 + i*(54670-48550)/LBINS;
  for(i=0; i<LBINS; i++)Lcenter[i] = i==0 ? Lup[0] : (Lup[i]+Lup[i-1])/2.;
  for(i=0; i<LBINS; i++)Lup_SB[i]     = 44.7 + i*(101.6-44.7)/LBINS;
  for(i=0; i<LBINS; i++)Lcenter_SB[i] = i==0 ? Lup_SB[0] : (Lup_SB[i]+Lup_SB[i-1])/2.;
  /* non-uniform binning: */
  for(i=0; i<EBINS; i++)Euptrue[i] = 124   + i*i*(3000-124)/EBINS/EBINS;
  for(i=0; i<EBINS; i++)Ecenter[i] = i==0 ? Euptrue[0] : (Euptrue[i]+Euptrue[i-1])/2.;

  /* fix units (cm, MeV to km, GeV) */
  for(i=0; i<LBINS; i++)Lup[i]     *= 1E-5;
  for(i=0; i<LBINS; i++)Lcenter[i] *= 1E-5;
  for(i=0; i<LBINS; i++)Lup_SB[i]     *= 1E-3;
  for(i=0; i<LBINS; i++)Lcenter_SB[i] *= 1E-3;
  for(i=0; i<EBINS; i++)Euptrue[i] *= 1E-3;
  for(i=0; i<EBINS; i++)Ecenter[i] *= 1E-3;

  for(i=0; i<11; i++)bincenter[i] = i==0 ? (0.2 + binwidths[i])/2. : 
                            (binwidths[i-1] + binwidths[i])/2.;
  for(int j=0; j<2; j++)
    for(i=0; i<21; i++)bincenter_MBSB[i][j] = i==0 ? bincenter_MBSB[i][j]/2. : 
                         (bincenter_MBSB[i-1][j] + bincenter_MBSB[i][j])/2.;

  /* Seting the beam contaminations using MC events from nubar disapp data: */
  contaminations();

  return 0;
}


/***************************************************************************
 * Cleanup GSL data structures required for MiniBooNE chi^2 calculation    *
 ***************************************************************************/
int chiMB_clear()
{
  /* if (osc_params) { glbFreeParams(osc_params); osc_params = NULL; } */
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
 * Probability matrix routine that uses the filter only when fast          *
 * oscillations are expected at the considered baseline and energy         *
 ***************************************************************************/
inline int mb_probability_matrix(double max_dmsq, double _P[3][3], int cp_sign, double Egev,
    int psteps, const double *length, const double *density,
    double filter_sigma, void *user_data)
{
  /* My adaptation for the filter (PANM) */
  if ((max_dmsq*length[0]*(1./Egev - 1/(Egev+0.2)) > 0.1*M_PI  ||  max_dmsq/Egev*0.060>0.3*M_PI) )
    return snu_probability_matrix(_P, cp_sign, Egev, psteps, length, density,
                                  filter_sigma, user_data);
  else
    return snu_probability_matrix(_P, cp_sign, Egev, psteps, length, density,
                                  0.0, user_data);
}


/***************************************************************************
 * Chi^2 for the MiniBooNE NEUTRINO analysis (6.46E20 POT) described at    *
 * http://www-boone.fnal.gov/for_physicists/data_release/lowe/             *
 ***************************************************************************/
double chiMB(int exper, int rule, int n_params, double *x, double *errors,
             void *user_data)
{
  /* The covariance matrices */
#include "./glb/mb-2012/mb-data/miniboone_full_fractcovmatrix_combined.txt"
#include "./glb/mb-2012/mb-data/miniboone_full_fractcovmatrix_combined_lowe.txt"
#ifndef COMBINED_APP
# include "./glb/mb-2012/mb-data/miniboone_full_fractcovmatrix_nu.txt"
# include "./glb/mb-2012/mb-data/miniboone_full_fractcovmatrix_nu_lowe.txt"
# include "./glb/mb-2012/mb-data/miniboone_full_fractcovmatrix_nubar.txt"
# include "./glb/mb-2012/mb-data/miniboone_full_fractcovmatrix_nubar_lowe.txt"
#endif
#include "./glb/mb-2012/mb-data/cov_matrix_numode_disap.dat"
//#include "./glb/mb-2012/mb-data/cov_matrix_nubarmode_disap.dat"

  double sig[11],sigbar[11],bg[11],bgbar[11],chi2 = 0.0, dummy[11][2];
  double sig_dis[16];
//  double numu[8],numuosc[8];

#ifndef COMBINED_APP
  double (*_M3nu)[NCOV]       = (double (*)[NCOV])     gsl_matrix_ptr(M3nu,    0, 0);
  double (*_M2nu)[NE+NMU]     = (double (*)[NE+NMU]) gsl_matrix_ptr(M2nu,    0, 0);
  double (*_M2invnu)[NE+NMU]  = (double (*)[NE+NMU]) gsl_matrix_ptr(M2invnu, 0, 0);
  double (*_M3bar)[NCOV]      = (double (*)[NCOV])     gsl_matrix_ptr(M3bar,    0, 0);
  double (*_M2bar)[NE+NMU]    = (double (*)[NE+NMU]) gsl_matrix_ptr(M2bar,    0, 0);
  double (*_M2invbar)[NE+NMU] = (double (*)[NE+NMU]) gsl_matrix_ptr(M2invbar, 0, 0);
#endif
  double (*_Mdis)[16]       = (double (*)[16]) gsl_matrix_ptr(Mdis,    0, 0);
  double (*_Minvdis)[16]    = (double (*)[16]) gsl_matrix_ptr(Minvdis, 0, 0);
//  double (*_Mdisbar)[16]    = (double (*)[16]) gsl_matrix_ptr(Mdisbar, 0, 0);
//  double (*_Minvdisbar)[16] = (double (*)[16]) gsl_matrix_ptr(Minvdisbar, 0, 0);
  double (*_M3)[2*NCOV]        = (double (*)[2*NCOV])     gsl_matrix_ptr(M3,    0, 0);
  double (*_M2)[2*(NE+NMU)]    = (double (*)[2*(NE+NMU)]) gsl_matrix_ptr(M2,    0, 0);
  double (*_M2inv)[2*(NE+NMU)] = (double (*)[2*(NE+NMU)]) gsl_matrix_ptr(M2inv, 0, 0);
  double (*_M_MBSB)[DIS_RANK]        = (double (*)[DIS_RANK]) gsl_matrix_ptr(M_MBSB,     0, 0);
  double (*_M_MBSB_inv)[DIS_RANK]    = (double (*)[DIS_RANK]) gsl_matrix_ptr(M_MBSB_inv, 0, 0);

  double prob_temp[NU_FLAVOURS][NU_FLAVOURS];
  double p = 0, p2=0, p3=0,p4=0;
  double *length, *density,(*pp)[NU_FLAVOURS];
  double dens[1] = {0}; density = &dens[0];
  double filter_sigma  = .2E-2;
  double filter_sigma2 = 1E-2;
  int i, j, k;
  double Pdis[16],P2dis[16],tot;
  int signum;
  int FIRST,TOT,LAST;
  FIRST = APP_FIRST;
  LAST = APP_LAST;
  TOT  = 11-FIRST;

#ifndef COMBINED_APP
  double P[NCOV]; /* Vector of pred nu_e signal, nu_e bg, nu_mu signal */
  double P2[NE+NMU]; /* Collapsed P */
#endif
  
  /* Check if oscillation parameters have changed since last call and if not, reuse the
   * old chi^2 */
  static glb_params last_params  = NULL;
  static glb_params these_params = NULL;
  static double last_chi2_mb = NAN;
  if (!these_params)
  {
    these_params = glbAllocParams();
    if (!these_params)
      return 1e35;
  }
  glbGetOscillationParameters(these_params);

  if (!last_params)
  {
    last_params = glbAllocParams();
    if (!last_params)
      return 1e36;
  }
  else if (!GLB_ISNAN(last_chi2_mb))
  {
    if (memcmp(these_params->osc->osc_params,
          last_params->osc->osc_params,
          sizeof(these_params->osc->osc_params[0]) * these_params->osc->length) == 0 &&
        memcmp(these_params->density->density_params,
          last_params->density->density_params,
          sizeof(these_params->density->density_params[0]) * these_params->density->length) == 0)
    {
      return last_chi2_mb;
    }
  }

  glb_params tmp_params = these_params;
  these_params = last_params;
  last_params  = tmp_params;
 

  /* Generating probability arrays in L and E bins */
  double max_dmsq = 0.0;  /* Find largest dmsq to determine if low-pass filter is used */
  for (i=2; i <= n_flavors; i++)
  {
    char s[4];
    sprintf(s, "DM%d1", i);
    double dmsq = fabs(glbGetOscParamByName(these_params, s));
    if (dmsq > max_dmsq)
      max_dmsq = dmsq;
  }
  for(i=0; i<EBINS; i++)
    for(j=0; j<LBINS; j++)
      {
        /* MiniBoone */
        length = &Lcenter[j];
        pp = prob_temp;
        mb_probability_matrix(max_dmsq, pp, 1, Ecenter[i],
                              1, length, density, filter_sigma, NULL);
        for(int A=0; A<NU_FLAVOURS; A++)
          for(int B=0; B<NU_FLAVOURS; B++)
            prob_nu[i][j][A][B] = prob_temp[A][B];

        /* MiniBoone */
        mb_probability_matrix(max_dmsq, pp, -1, Ecenter[i],
                              1, length, density, filter_sigma, NULL);
        for(int A=0; A<NU_FLAVOURS; A++)
          for(int B=0; B<NU_FLAVOURS; B++)
            prob_bar[i][j][A][B] = prob_temp[A][B];

        /* SciBoone */
        length = &Lcenter_SB[j];
        pp = prob_temp;
        mb_probability_matrix(max_dmsq, pp, 1, Ecenter[i],
                              1, length, density, filter_sigma, NULL);
        for(int A=0; A<NU_FLAVOURS; A++)
          for(int B=0; B<NU_FLAVOURS; B++)
            prob_nu_SB2[i][j][A][B] = prob_temp[A][B];

        /* SciBoone */
        mb_probability_matrix(max_dmsq, pp, -1, Ecenter[i],
                              1, length, density, filter_sigma, NULL);
        for(int A=0; A<NU_FLAVOURS; A++)
          for(int B=0; B<NU_FLAVOURS; B++)
            prob_bar_SB2[i][j][A][B] = prob_temp[A][B];

      }

  /* Generating probability for SciBooNE baseline 
     (this should be used just for the appearance background) */
  double LSB[1] = {0.100};
  for(i=0; i<11; i++)
      {
        length = &LSB[0];
        pp = prob_temp;
        mb_probability_matrix(max_dmsq, pp, 1, bincenter[i],
                              1, length, density, filter_sigma2, NULL);
        for(int A=0; A<NU_FLAVOURS; A++)
          for(int B=0; B<NU_FLAVOURS; B++)
            prob_nu_SB[i][A][B] = prob_temp[A][B];

        mb_probability_matrix(max_dmsq, pp, -1, bincenter[i],
                              1, length, density, filter_sigma2, NULL);
        for(int A=0; A<NU_FLAVOURS; A++)
          for(int B=0; B<NU_FLAVOURS; B++)
            prob_bar_SB[i][A][B] = prob_temp[A][B];
      }

#if defined(COMBINED_APP) || defined(NU_APP)
  /************* NEUTRINO EVENTS *************/
  /* APPEARANCE */
  for (i=0; i < 11; i++) sig[i] = 0.0;
  
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
      p  = prob_nu[i][j][1][0];  /* mu -> e */
      p2 = prob_nu[i][j][1][1];  /* mu -> mu */
      double Eup = Emin + 1000.*binwidths[FIRST];
      i = 0;
      while (Ereco >= Eup && i < TOT)
        Eup += 1000.*binwidths[++i];
      if (i < TOT && Ereco >= Emin)
#ifdef FULL_OSC
        sig[i] += p/p2 * w; /* Pme/Pmm - see discussion with W Louis */
          /* Flux is normalized to \nu_\mu rate, i.e. if there is \nu_\mu
           * disappearance, the flux is underestimated by 1 / Pmm */
#else
        sig[i] += p * w; /* Pme/Pmm - see discussion with W Louis */
#endif
    }
  
  for(i=0; i < TOT; i++)
    sig[i] = sig[i]/Neventsnu;

  /* beam nue background oscillation */
  for(int k=0; k<TOT; k++)
    {
      /* Localize event in Etrue -- L grid */
      i = 0; j = 0;
      while (Euptrue[i] <= bincenter[FIRST+k] && i <= EBINS-2) i++;
      while (Lup[j] <= 0.520 && j <= LBINS-2) j++;
      /* prob rescaling - see discussion with W Louis */
#ifdef FULL_OSC
      dummy[k][0] = prob_nu[i][j][0][0]/prob_nu[i][j][1][1]*NUEbgnuMU[FIRST+k];
      dummy[k][1] = prob_nu[i][j][0][0]/prob_nu_SB[FIRST+k][0][0]*NUEbgnuK[FIRST+k];
      bg[k] = OTHERbgnu[k] + dummy[k][0] + dummy[k][1];
      /* bg[k] = OTHERbgnu[FIRST+k]  */
      /*   + prob_nu[i][j][0][0]/prob_nu[i][j][1][1]*NUEbgnuMU[FIRST+k]  */
      /*   + prob_nu[i][j][0][0]/prob_nu_SB[FIRST+k][0][0]*NUEbgnuK[FIRST+k]; */
      /* Rescale beam \nu_e backgrounds by Pee/Pmm (flux is normalized
       * to \nu_\mu flux, i.e. if the latter oscillates, a mistake is made);
       * Kaon backgrounds are normalized toe SB data, so use Pmm at the SciBooNE
       * baseline for them. */

  /* Exporting background */
  for(i=0; i<TOT; i++)
    {
      app_nu_spectrum[i][1] = OTHERbgnu[i];    
      app_nu_spectrum[i][2] = dummy[i][1]; /* nue from K */
      app_nu_spectrum[i][3] = dummy[i][0]; /* nue from mu */
    }
#else
      bg[k] = ALLbgnu[FIRST+k]; 
#endif
    }

  /* Exporting signal */
  for(i=0; i<TOT; i++) app_nu_spectrum[i][0] = sig[i];
  
#ifndef COMBINED_APP
  /* Construct covariance matrix in 3x3 block form */
  for (i=0; i < TOT; i++)
    { P[i]      = sig[i];
      P[i+NE] = bg[i];  }
  for (i=0; i < NMU; i++)
    P[i+2*NE] = pred_numu[i];
  for (i=0; i < NCOV; i++)
    for (j=0; j < NCOV; j++)
      _M3nu[i][j] = NE==11 ? Mfrac_nu_200[i][j] * P[i] * P[j]
                           : Mfrac_nu_475[i][j] * P[i] * P[j];
  for (i=0; i < NE; i++)        /* Add statistical error of signal */
    _M3nu[i][i] += sig[i];
      
  /* Collapse covariance matrix to 2x2 block form */
  for (i=0; i < NE; i++)        /* Upper left block */
    for (j=0; j < NE; j++)
      _M2nu[i][j] = _M3nu[i][j] + _M3nu[i+NE][j] 
                  + _M3nu[i][j+NE] + _M3nu[i+NE][j+NE];
  for (i=0; i < NMU; i++)       /* Lower left block */
    for (j=0; j < NE; j++)
      _M2nu[i+NE][j] = _M3nu[i+2*NE][j] + _M3nu[i+2*NE][j+NE];
  for (i=0; i < NE; i++)       /* Upper right block */
    for (j=0; j < NMU; j++)
      _M2nu[i][j+NE] = _M3nu[i][j+2*NE] + _M3nu[i+NE][j+2*NE];
  for (i=0; i < NMU; i++)      /* Lower right block */
    for (j=0; j < NMU; j++)
      _M2nu[i+NE][j+NE] = _M3nu[i+2*NE][j+2*NE];
      
  /* Invert covariance matrix and compute log-likelihood */
  gsl_linalg_LU_decomp(M2nu, permnu, &signum);
  gsl_linalg_LU_invert(M2nu, permnu, M2invnu);

  for (i=0; i < NE; i++)
    P2[i] = ALLdatanu[FIRST+i] - (sig[i] + bg[i]);
  for (i=0; i < NMU; i++)
    P2[i+NE] = ALLdatanu[i+11] - pred_numu[i];
  for (i=0; i < NE+NMU; i++)
    for (j=0; j < NE+NMU; j++)
      chi2 += P2[i] * _M2invnu[i][j] * P2[j];
#endif  /* not def COMBINED_APP */
#endif  /* any APP defined */

#ifdef NU_DISAPP
      /* DISAPPEARANCE */
      for (i=0; i<16; i++) sig_dis[i] = 0;
      for (long k=0; k < Neventsnumu; k++)
        {
          double Ereco, Etrue, L, w;
          Ereco = mc_eventsnumu[k][0];
          Etrue = mc_eventsnumu[k][1];
          L     = mc_eventsnumu[k][2];
          w     = mc_eventsnumu[k][3];
          
          /* Localize event in Etrue -- L grid */
          i = 0; j = 0;
          while (Euptrue[i] <= Etrue && i <= EBINS-2) i++;
          while (Lup[j] <= L && j <= LBINS-2) j++;
          p = prob_nu[i][j][1][1]; /* mu -> mu */
          double Eup = disbin[0];
          i = 0;
          while (Ereco >= Eup && i < 16)
            Eup = disbin[++i];
          if (i < 16 && Ereco >= 0)
            sig_dis[i] += p * w;
        }
      
      tot=0;
      for(i=0; i<16; i++)
        tot += sig_dis[i];
      
      /* Construct covariance matrix for disap */
      for (i=0; i < 16; i++)
        /* this could be a possibility for avoiding minimization: */
        Pdis[i] = 190454.0/tot*sig_dis[i];
      /* Pdis[i] = x[0]>-1 ?  (1.0+x[0])*190454.0/tot*sig_dis[i] */
      /*                   : -(1.0+x[0])*190454.0/tot*sig_dis[i]; */
      for (i=0; i < 16; i++)
        for (j=0; j < 16; j++)
          _Mdis[i][j] = M_dis[i][j] * Pdis[i] * Pdis[j];
      for (i=0; i < 16; i++)        /* Add statistical error of signal */
        _Mdis[i][i] += Pdis[i];
      
      /* Invert covariance matrix and compute log-likelihood */
      gsl_linalg_LU_decomp(Mdis, permdis, &signum);
      gsl_linalg_LU_invert(Mdis, permdis, Minvdis);
      
      for (i=0; i < 16; i++)
        P2dis[i] = nudisdata[i] - Pdis[i];
      for (i=0; i < 16; i++)
        for (j=0; j < 16; j++)
          chi2 += P2dis[i] * _Minvdis[i][j] * P2dis[j];
#endif

  /************* ANTINEUTRINO MODE *************/
#if defined(NUBAR_APP) || defined(COMBINED_APP)
  /* Compute signal using Monte Carlo events */
  for (i=0; i < 11; i++) sigbar[i] = 0.0;
  
  for (long k=0; k < Neventsbar; k++)
    {
      double Ereco, Etrue, L, w;
      Ereco = mc_eventsbar[k][0];
      Etrue = 1E-3*mc_eventsbar[k][1]; /* in GeV */
      L     = 1E-5*mc_eventsbar[k][2]; /* in km */
      w     = mc_eventsbar[k][3];
      
      /* Localize event in Etrue -- L grid */
      i = 0; j = 0;
      while (Euptrue[i] <= Etrue && i <= EBINS-2) i++;
      while (Lup[j] <= L && j <= LBINS-2) j++;
      p = prob_bar[i][j][1][0];  /* mub -> eb  */
      p2 = prob_nu[i][j][1][0];  /* mu  -> e   */
      p3 = prob_bar[i][j][1][1]; /* mub -> mub   */
      p4 = prob_nu[i][j][1][1];  /* mu  -> mu   */
      double Eup = Emin + 1000.*binwidths[FIRST];
      i = 0;
      while (Ereco >= Eup && i < TOT)
        Eup += 1000.*binwidths[++i];
      if (i < TOT && Ereco >= Emin) /* Pme/Pmm - see discussion with W Louis */
#ifdef FULL_OSC
        sigbar[i] += (p/p3 + p2/p4*numufraction[FIRST+i]) * w; 
#else
        sigbar[i] += p * w; 
#endif
      /* Here, the numu fraction in the nubar mode is assumed to oscillate! */
    }

  for (i=0; i < TOT; i++)
    sigbar[i] = sigbar[i]/Neventsbar;
    
  /* beam nue background oscillation */
  for(k=0; k<TOT; k++)
    {
      /* Localize event in Etrue -- L grid */
      i = 0; j = 0;
      while (Euptrue[i] <= bincenter[FIRST+k] && i <= EBINS-2) i++;
      while (Lup[j] <= 0.520 && j <= LBINS-2) j++;
      /* prob rescaling - see discussion with W Louis */
#ifdef FULL_OSC
      dummy[k][0] = ( (1.0 - nuefraction[FIRST+k])*prob_bar[i][j][0][0]/prob_bar[i][j][1][1]
		      + nuefraction[FIRST+k]*prob_nu[i][j][0][0]/prob_nu[i][j][1][1] )*NUEbgbarMU[FIRST+k];
      dummy[k][1] = ( (1.0 - nuefraction[FIRST+k])*prob_bar[i][j][0][0]/prob_bar_SB[FIRST+k][0][0]
	    + nuefraction[FIRST+k]*prob_nu[i][j][0][0]/prob_nu_SB[FIRST+k][0][0] )*NUEbgbarK[FIRST+k];
      bgbar[k] = OTHERbgbar[k] + dummy[k][0] + dummy[k][1];
      /* bgbar[k] = OTHERbgbar[FIRST+k] */
      /*   + ( (1.0 - nuefraction[FIRST+k])*prob_bar[i][j][0][0]/prob_bar[i][j][1][1] */
      /*       + nuefraction[FIRST+k]*prob_nu[i][j][0][0]/prob_nu[i][j][1][1] )*NUEbgbarMU[FIRST+k] */
      /*   + ( (1.0 - nuefraction[FIRST+k])*prob_bar[i][j][0][0]/prob_bar_SB[FIRST+k][0][0] */
      /*       + nuefraction[FIRST+k]*prob_nu[i][j][0][0]/prob_nu_SB[FIRST+k][0][0] )*NUEbgbarK[FIRST+k]; */

  /* Exporting background */
  for(i=0; i<TOT; i++)
    {
      app_bar_spectrum[i][1] = OTHERbgbar[i];    
      app_bar_spectrum[i][2] = dummy[i][1]; /* nue from K */
      app_bar_spectrum[i][3] = dummy[i][0]; /* nue from mu */
    }
#else
      bgbar[k] = ALLbgbar[FIRST+k];
#endif
    }

  /* Exporting signal */
  for(i=0; i<TOT; i++) app_bar_spectrum[i][0] = sigbar[i];
  
#ifndef COMBINED_APP
  /* Construct covariance matrix in 3x3 block form */
  for (i=0; i < TOT; i++)
    { P[i]       = sigbar[i];
      P[i+NE] = bgbar[i];  }
  for (i=0; i < NMU; i++)
    P[i+2*NE] = pred_numubar[i];
  for (i=0; i < NCOV; i++)
    for (j=0; j < NCOV; j++)
      {
        _M3bar[i][j] = NE==11 ? Mfrac_nubar_200[i][j] * P[i] * P[j]
                              : Mfrac_nubar_475[i][j] * P[i] * P[j];
      }
  for (i=0; i < NE; i++)        /* Add statistical error of signal */
    _M3bar[i][i] += sigbar[i];

  /* Collapse covariance matrix to 2x2 block form */
  for (i=0; i < NE; i++)        /* Upper left block */
    for (j=0; j < NE; j++)
      _M2bar[i][j] = _M3bar[i][j] + _M3bar[i+NE][j] + _M3bar[i][j+NE] + _M3bar[i+NE][j+NE];
  for (i=0; i < NMU; i++)       /* Lower left block */
    for (j=0; j < NE; j++)
      _M2bar[i+NE][j] = _M3bar[i+2*NE][j] + _M3bar[i+2*NE][j+NE];
  for (i=0; i < NE; i++)       /* Upper right block */
    for (j=0; j < NMU; j++)
      _M2bar[i][j+NE] = _M3bar[i][j+2*NE] + _M3bar[i+NE][j+2*NE];
  for (i=0; i < NMU; i++)      /* Lower right block */
    for (j=0; j < NMU; j++)
      _M2bar[i+NE][j+NE] = _M3bar[i+2*NE][j+2*NE];

  /* Invert covariance matrix and compute log-likelihood */
  gsl_linalg_LU_decomp(M2bar, permbar, &signum);
  gsl_linalg_LU_invert(M2bar, permbar, M2invbar);
  
  for (i=0; i < NE; i++)
    P2[i] = ALLdatabar[FIRST+i] - (sigbar[i] + bgbar[i]);
  for (i=0; i < NMU; i++)
    P2[i+NE] = ALLdatabar[i+LAST] - pred_numubar[i];
  for (i=0; i < NE+NMU; i++)
    for (j=0; j < NE+NMU; j++)
      chi2 += P2[i] * _M2invbar[i][j] * P2[j];
  /* chi2 += gsl_linalg_LU_lndet(M2bar); */
#endif
#endif

#ifdef NUBAR_DISAPP
  /* DISAPPEARANCE ANTINEUTRINO - MBSB */
  double sig_MBSB[21][2][2];    /* MB-SB, RS-WS */
  for (i=0; i<21; i++) 
    for (j=0; j<2; j++) 
      for (k=0; k<2; k++) sig_MBSB[i][j][k] = 0;

  /* MB events: */
  for (long k=0; k < Nevents_MB; k++)
    {
      int type; double Ereco, Etrue, L, w;
      type  = mc_events_MB_type[k];
      Etrue = mc_events_MB[k][0];
      Ereco = mc_events_MB[k][1];
      L     = 1E-3*mc_events_MB[k][2];
      w     = mc_events_MB[k][3];
      
      /* Localize event in Etrue -- L grid */
      i = 0; j = 0;
      while (Euptrue[i] <= Etrue && i <= EBINS-2) i++;
      while (Lup[j] <= L && j <= LBINS-2) j++;
#ifdef FULL_OSC
      if(type==14)  p = prob_nu[i][j][1][1];  /* mu    -> mu */
      if(type==-14) p = prob_bar[i][j][1][1]; /* mubar -> mubar */
      if(type==12)  p = prob_nu[i][j][0][1];  /* e     -> mu */
      if(type==-12) p = prob_bar[i][j][0][1]; /* ebar  -> mubar */
#else
      if(type==14)  p = 1.0;  /* mu    -> mu */
      if(type==-14) p = prob_bar[i][j][1][1]; /* mubar -> mubar */
      if(type==12)  p = 0;  /* e     -> mu */
      if(type==-12) p = 0; /* ebar  -> mubar */
#endif
      double Eup = bins_MBSB[0][0];
      i = 0;
      while (Ereco >= Eup && i < 21)
        Eup = bins_MBSB[++i][0];
      if (i < 21 && Ereco >= 0.3 && Ereco < 1.9)
        { if(type<0)sig_MBSB[i][0][0] += p * w; /* RS */
          if(type>0)sig_MBSB[i][0][1] += p * w; } /* WS */
    }

  /* CORRECT LAST BINS (WHY?) */
  /* sig_MBSB[19][0][1] *= 240.9/248.6; */
  /* sig_MBSB[20][0][1] *= 210.6/258.9; */
  /* sig_MBSB[19][0][0] *= 340.6/358.7; */
  /* sig_MBSB[20][0][0] *= 216.4/301.1; */

  /* SB events: */
  for (long k=0; k < Nevents_SB; k++)
    {
      int type; double Ereco, Etrue, L, w;
      type  = mc_events_SB_type[k];
      Etrue = mc_events_SB[k][0];
      Ereco = mc_events_SB[k][1];
      L     = 1E-3*mc_events_SB[k][2];
      w     = mc_events_SB[k][3];
      
      /* Localize event in Etrue -- L grid */
      i = 0; j = 0;
      while (Euptrue[i] <= Etrue && i <= EBINS-2) i++;
      while (Lup_SB[j] <= L && j <= LBINS-2) j++;
      /* if(L>0.105 || L<0.040)printf("%g!!!!!!!!!!\n",L); */
#ifdef FULL_OSC
      if(type==14)  p = prob_nu_SB2[i][j][1][1];  /* mu    -> mu */
      if(type==-14) p = prob_bar_SB2[i][j][1][1]; /* mubar -> mubar */
      if(type==12)  p = prob_nu_SB2[i][j][0][1];  /* e     -> mu */
      if(type==-12) p = prob_bar_SB2[i][j][0][1]; /* ebar  -> mubar */
#else
      if(type==14)  p = 1.0;  /* mu    -> mu */
      if(type==-14) p = prob_bar_SB2[i][j][1][1]; /* mubar -> mubar */
      if(type==12)  p = 0;  /* e     -> mu */
      if(type==-12) p = 0; /* ebar  -> mubar */
#endif
      double Eup = bins_MBSB[0][1];
      i = 0;
      while (Ereco >= Eup && i < 21)
        Eup = bins_MBSB[++i][1];
      if (i < 21 && Ereco >= 0.3 && Ereco < 1.9)
      { if(type<0)sig_MBSB[i][1][0] += p * w; /* RS */
        if(type>0)sig_MBSB[i][1][1] += p * w; } /* WS */
    }

  /* CORRECT FIRST BINS (WHY?) */
  /* sig_MBSB[0][1][1] *= 338.9/320.8; */
  /* sig_MBSB[1][1][1] *= 304.9/294.3; */
  /* sig_MBSB[0][1][0] *= 658.4/644.26; */
  /* sig_MBSB[1][1][0] *= 610.7/602.3; */

  /* for(i=0; i<21; i++) */
  /*   printf("%d   %f   %f\n",i,sig_MBSB[i][0][0],sig_MBSB[i][0][1]); */
  /* for(i=0; i<21; i++) */
  /*   printf("%d   %f   %f\n",i,sig_MBSB[i][1][0],sig_MBSB[i][1][1]); */
  /* exit(0); */

  for (i=0; i < 42; i++)
    for (j=0; j < 42; j++)
      _M_MBSB[i][j] = 0;
  for (i=0; i < 21; i++)
    for (j=0; j < 21; j++)
      for (int E1=0; E1 < 2; E1++)
        for (int E2=0; E2 < 2; E2++)
          for (int S1=0; S1 < 2; S1++)
            for (int S2=0; S2 < 2; S2++)
              _M_MBSB[i+21*E1][j+21*E2] += MfracMBSB[i][j][E1][E2][S1][S2] 
                * sig_MBSB[i][E1][S1] * sig_MBSB[j][E2][S2];
  for (i=0; i < 21; i++)        /* Add statistical error of signal */
    {
      _M_MBSB[i][i] += sig_MBSB[i][0][0]+sig_MBSB[i][0][1];
      _M_MBSB[i+21][i+21] += sig_MBSB[i][1][0]+sig_MBSB[i][1][1];
    }

  /* Invert covariance matrix and compute log-likelihood */
  gsl_linalg_LU_decomp(M_MBSB, perm_MBSB, &signum);
  gsl_linalg_LU_invert(M_MBSB, perm_MBSB, M_MBSB_inv);
  
  double V[DIS_RANK];
  for (i=0; i < 21; i++)
    {
      V[i]    = MBSB_data[i][0] - sig_MBSB[i][0][0] - sig_MBSB[i][0][1];
      V[i+21] = MBSB_data[i][1] - sig_MBSB[i][1][0] - sig_MBSB[i][1][1];
    }

  for (i=0; i < 42; i++)
    for (j=0; j < 42; j++)
      chi2 += V[i] * _M_MBSB_inv[i][j] * V[j];
#endif

  /* combined analysis for appearance */
#ifdef COMBINED_APP
  double Pcomb[2*NCOV];
  /* Construct covariance matrix in 6x6 block form */
  for (i=0; i < TOT; i++)
    {
      Pcomb[i]         = sig[i];
      Pcomb[i+NE]      = bg[i];  
      Pcomb[i+NCOV]    = sigbar[i];  
      Pcomb[i+NCOV+NE] = bgbar[i];  
     }
  for (i=0; i < NMU; i++)
    {
      Pcomb[i+2*NE]      = pred_numu[i];
      Pcomb[i+NCOV+2*NE] = pred_numubar[i];
    }
  for (i=0; i < 2*NCOV; i++)
    for (j=0; j < 2*NCOV; j++)
      {
        _M3[i][j] = NE==11 ? Mfrac_comb_200[i][j] * Pcomb[i] * Pcomb[j]
                           : Mfrac_comb_475[i][j] * Pcomb[i] * Pcomb[j];
      }
  for (i=0; i < NE; i++)        /* Add statistical error of signal */
    {
      _M3[i][i] += sig[i];
      _M3[NCOV+i][NCOV+i] += sigbar[i];
    }

  /* Collapse covariance matrix to 4x4 block form */
  int i1=-1,i2=-1,j1,j2;
  for(i=0; i<2*(NE+NMU); i++)
    for(j=0; j<2*(NE+NMU); j++)
      {
        if(i<NE)                /* nue */
          { i1=i;      i2=i+NE; }
        else if(i<NE+NMU)       /* numu */
          { i1=i+NE;   i2=1000; }
        else if(i<2*NE+NMU)     /* nuebar */
          { i1=i+NE;   i2=i+2*NE; }
        else if(i<2*NE+2*NMU)   /* numubar */
          { i1=i+2*NE; i2=1000; }

        if(j<NE)
          { j1=j;     j2=j+NE; }
        else if(j<NE+NMU)
          { j1=j+NE;  j2=1000; }
        else if(j<2*NE+NMU)
          { j1=j+NE;  j2=j+2*NE; }
        else if(j<2*NE+2*NMU)
          { j1=j+2*NE;j2=1000; }

        _M2[i][j] = _M3[i1][j1];
        if(i2<100) _M2[i][j] += _M3[i2][j1];
        if(j2<100) _M2[i][j] += _M3[i1][j2];
        if(i2<100 && j2<100) _M2[i][j] += _M3[i2][j2];

        /* printf("%d  %d  %d  %d  %d  %d  ===  %g  %g\n",i,j,i1,j1,i2,j2,_M2[i][j],_M3[i1][j1]); */
      }

  /* exit(0); */

  /* Invert covariance matrix and compute log-likelihood */
  gsl_linalg_LU_decomp(M2, perm, &signum);
  gsl_linalg_LU_invert(M2, perm, M2inv);
  
  double P2comb[2*(NE+NMU)];
  for (i=0; i < NE; i++)
    {
      P2comb[i]        = ALLdatanu[FIRST+i] - (sig[i] + bg[i]);
      P2comb[i+NE+NMU] = ALLdatabar[FIRST+i] - (sigbar[i] + bgbar[i]);
    }
  for (i=0; i < NMU; i++)
    {
      P2comb[i+NE] = ALLdatanu[i+LAST] - pred_numu[i];
      P2comb[i+2*NE+NMU] = ALLdatabar[i+LAST] - pred_numubar[i];
    }
  for (i=0; i < 2*(NE+NMU); i++)
    for (j=0; j < 2*(NE+NMU); j++)
      chi2 += P2comb[i] * _M2inv[i][j] * P2comb[j];
#endif

//  printf("MB chi^2: %g\n",chi2);
  last_chi2_mb = chi2;
  return chi2;
}


