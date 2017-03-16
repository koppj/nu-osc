#include "definitions.h"

namespace ns_reactor
{

// SBL data updated from sterile neutrino white paper Mar 2012

extern Fit  fit;
extern Rate_pull_coupl rate_pull_coupl;  // defined in class_flux.cc
extern Rate_coef rate;                   // defined in class_flux.cc
extern int old_new_main;

bool use_sbl_react[N_SBL_R];

/********* the data *********************/

// tab. XIX of white paper
// rate_old, rate_new
double sbl_data[N_SBL_R][2] = {
  {0.987, 0.926},         // BUGEY4
  {0.985, 0.924},	  // ROVNO   
  {0.988, 0.930}, 	  // BUGEY3_1
  {0.994, 0.936}, 	  // BUGEY3_2
  {0.915, 0.861},	  // BUGEY3_3
  {1.018, 0.949}, 	  // GOSGEN_1
  {1.045, 0.975}, 	  // GOSGEN_2
  {0.975, 0.909},	  // GOSGEN_3
  {0.832, 0.788},	  // ILL     
  {1.013, 0.920}, 	  // KRASN_1 
  {1.031, 0.937},	  // KRASN_2 
  {0.989, 0.931},	  // KRASN_3 
  {0.987, 0.936},	  // SRP_1 
  {1.055, 1.001},	  // SRP_2 
  {0.969, 0.901},	  // ROV_1I 
  {1.001, 0.932},	  // ROV_2I 
  {1.026, 0.955},	  // ROV_1S 
  {1.013, 0.943},	  // ROV_2S 
  {0.990, 0.922}	  // ROV_3S 
};


// uncorrelated error in percent
double sbl_err_uncor[N_SBL_R] = {
  1.09,   // BUGEY4
  2.10,	  // ROVNO   
  2.05,   // BUGEY3_1
  2.06,   // BUGEY3_2
  14.6,	  // BUGEY3_3
  2.38,	  // GOSGEN_1
  2.31,	  // GOSGEN_2
  4.81,	  // GOSGEN_3
  8.52,	  // ILL     
  3.55,	  // KRASN_1 
  19.8,	  // KRASN_2 
  2.67,	  // KRASN_3 
  1.95,	  // SRP_1 
  2.11,	  // SRP_2 
  4.24,	  // ROV_1I 
  4.24,	  // ROV_2I 
  4.95,	  // ROV_1S 
  4.95,	  // ROV_2S 
  4.53	  // ROV_3S 
};


// isotope fractions [1101.2755]
// enum Isotopes {U235, U238, P239, P241, NISO};
// Goesgen numbers from PRD 34, 2636
const double sbl_iso_fract[N_SBL_R][NISO] = {
  {0.538, 0.078, 0.328, 0.056},           // BUGEY4
  {0.614, 0.074, 0.275, 0.031},           // ROVNO   
  {0.538, 0.078, 0.328, 0.056}, 	  // BUGEY3_1
  {0.538, 0.078, 0.328, 0.056}, 	  // BUGEY3_2
  {0.538, 0.078, 0.328, 0.056}, 	  // BUGEY3_3
  {0.619, 0.067, 0.272, 0.042}, 	  // GOSGEN_1
  {0.584, 0.068, 0.298, 0.050}, 	  // GOSGEN_2
  {0.543, 0.070, 0.329, 0.058}, 	  // GOSGEN_3
  {0.930, 0.000, 0.070, 0.000},	          // ILL 
  {1., 0., 0., 0.}, 	                  // KRASN_1 
  {1., 0., 0., 0.},             	  // KRASN_2 
  {1., 0., 0., 0.},              	  // KRASN_3 
  {1., 0., 0., 0.},              	  // SRP_1 
  {1., 0., 0., 0.},              	  // SRP_2 
  {0.607, 0.074, 0.277, 0.042},	          // ROV_1I 
  {0.603, 0.076, 0.276, 0.045},	          // ROV_2I 
  {0.606, 0.074, 0.277, 0.043},	          // ROV_1S 
  {0.557, 0.076, 0.313, 0.054},	          // ROV_2S 
  {0.606, 0.074, 0.274, 0.046}	          // ROV_3S 
};


// the baselines
const double sbl_L[N_SBL_R] = {
  15.,    // BUGEY4
  18.,    // ROVNO   
  15., 	  // BUGEY3_1
  40., 	  // BUGEY3_2
  95.,	  // BUGEY3_3
  38.,	  // GOSGEN_1
  45.,	  // GOSGEN_2
  65.,	  // GOSGEN_3
  9.,     // ILL 
  33.,    // KRASN_1 
  92.,	  // KRASN_2 
  57.,	  // KRASN_3 
  18.,	  // SRP_1 
  24.,	  // SRP_2 
  18.,	  // ROV_1I 
  18.,	  // ROV_2I 
  18.,	  // ROV_1S 
  25.,	  // ROV_2S 
  18.	  // ROV_3S 
};

/**********************************************************/
/***************** set data and errors in Fit *************/
/**********************************************************/

int init_sbl_reactors(int old_new)
{
  // set data and errors in fit
  for(int i = 0; i < N_SBL_R; i++){
    int j = fit.first_bin[SBL] + i;

    fit.Data[j] = sbl_data[i][old_new];
    fit.S_data[j][j] = norm(0.01 * sbl_err_uncor[i] * fit.Data[j]);
  }

  // correlate syst. errors

  // construct correlation matrix (relat. errors)
  double S[N_SBL_R][N_SBL_R];
  for(int i = 0; i < N_SBL_R; i++)
    for(int j = 0; j < N_SBL_R; j++)
      S[i][j] = 0.;

  // BUGEY4 and ROVNO
  S[BUGEY4][BUGEY4] = norm(0.0084); 
  S[BUGEY4][ROVNO] = S[ROVNO][BUGEY4] = 0.0084 * 0.018; 
  S[ROVNO][ROVNO] = norm(0.018);

  // BUGEY3
  for(int i = BUGEY3_1; i <= BUGEY3_3; i++)
    for(int j = BUGEY3_1; j <= BUGEY3_3; j++)
      S[i][j] = norm(0.039);   

  // GOSGEN 
  for(int i = GOSGEN_1; i <= GOSGEN_3; i++)
    for(int j = GOSGEN_1; j <= GOSGEN_3; j++)
      S[i][j] = norm(0.048);

  // GOSGEN & ILL
  for(int i = GOSGEN_1; i <= GOSGEN_3; i++)
    S[i][ILL] = S[ILL][i] = 0.0436 * 0.08;

  S[ILL][ILL] = norm(0.08);

  // Krasnoyarsk
  S[KRASN_1][KRASN_1] = norm(0.0484);
  S[KRASN_1][KRASN_2] = S[KRASN_2][KRASN_1] = 0.0484 * 0.0476;
  S[KRASN_1][KRASN_3] = S[KRASN_3][KRASN_1] = 0.03 * 0.034;
  S[KRASN_2][KRASN_2] = norm(0.0476);
  S[KRASN_2][KRASN_3] = S[KRASN_3][KRASN_2] = 0.03 * 0.034;
  S[KRASN_3][KRASN_3] = norm(0.034);

  // SRP
  S[SRP_1][SRP_1] = S[SRP_1][SRP_2] = S[SRP_2][SRP_1] = S[SRP_2][SRP_2] = 
    norm(0.02);

  // ROVNO88
  for(int i = ROV_1I; i <= ROV_3S; i++)
    for(int j = ROV_1I; j <= ROV_3S; j++)
      S[i][j] = norm(0.022);

  for(int i = ROV_1I; i <= ROV_2I; i++)
    for(int j = ROV_1I; j <= ROV_2I; j++)
      S[i][j] += sbl_err_uncor[i] * sbl_err_uncor[j] * 1.e-4;

  for(int i = ROV_1S; i <= ROV_3S; i++)
    for(int j = ROV_1S; j <= ROV_3S; j++)
      S[i][j] += sbl_err_uncor[i] * sbl_err_uncor[j] * 1.e-4;


  // set the covar matrix in fit
  for(int i = 0; i < N_SBL_R; i++){
    int ii = fit.first_bin[SBL] + i;

    for(int j = 0; j < N_SBL_R; j++){
      int jj = fit.first_bin[SBL] + j;

      fit.S_data[ii][jj] += S[i][j] * fit.Data[ii] * fit.Data[jj];
    }
  }
   
  /* total correlation of all reactors
  for(int i = 0; i < N_SBL_R; i++){
    const int ii = fit.first_bin[SBL] + i;
    for(int j = 0; j < N_SBL_R; j++){
      const int jj = fit.first_bin[SBL] + j;

      fit.S_data[ii][jj] += norm(0.02) * fit.Data[ii] * fit.Data[jj];
    }
  }
  */

  for(int i = 0; i < N_SBL_R; i++)
  {
#ifdef USE_SBL
    use_sbl_react[i] = true;
#else
    use_sbl_react[i] = false;
#endif
  }

#ifdef USE_BUGEY_SP
  use_sbl_react[BUGEY3_1] = use_sbl_react[BUGEY3_2] = use_sbl_react[BUGEY3_3] = false;
#endif

  int nn = 0;
  for(int i = 0; i < N_SBL_R; i++)
    if(use_sbl_react[i]) nn++;
//  fprintf(stderr, "SBL (n_sbl_react = %d) ", nn);

  return nn;
}

/**********************************************************************/
/*************** routine called by fit.chisq **************************/
/**********************************************************************/

void set_table_sbl(Param_5nu &p, double cff[NBIN_CHISQ][NPULLS+1])
{
  for(int i = 0; i < N_SBL_R; i++){
    int ii = fit.first_bin[SBL] + i;

    if(use_sbl_react[i]){

      cff[ii][NPULLS] = cff[ii][FLUX_NORM] = rate.prob(p, sbl_L[i], sbl_iso_fract[i]);

      // the pulls
      for(int j = 0; j < N_CO_UNC; j++){
        cff[ii][PULL_U235_0 + j] = rate_pull_coupl.unc(U235, j, sbl_iso_fract[i]);
        cff[ii][PULL_P239_0 + j] = rate_pull_coupl.unc(P239, j, sbl_iso_fract[i]); 
        cff[ii][PULL_P241_0 + j] = rate_pull_coupl.unc(P241, j, sbl_iso_fract[i]); 
      }

      cff[ii][PULL_U238] = rate_pull_coupl.unc(U238, 0, sbl_iso_fract[i]);
      cff[ii][FLUX_COR]  = rate_pull_coupl.cor(sbl_iso_fract[i]);

    }else{

      for(int p = 0; p < PULL_GLOBAL; p++)
	cff[ii][p] = 0.;

      cff[ii][NPULLS] = fit.Data[ii];
    }
  }
  return;
}

/**********************************************************************/

void print_sbl_data(void)
{
  double L18 = 17.;
  for(int i = 0; i < N_SBL_R; i++){
    double l = sbl_L[i];
    if(sbl_L[i] == 18.){ L18 += 0.4; l = L18;}  
    const int ii = fit.first_bin[SBL] + i;
    printf("%d %e %e %e %e\n", i+1, l, sbl_data[i][NEW], 
           0.01*sbl_err_uncor[i] * sbl_data[i][NEW], sqrt(fit.S_data[ii][ii])); 
  }
  exit(0);
  return;
}

void print_sbl_pred(Param_5nu p, FILE *fp)
{
  for(double L = 5.; L < 200.; L *= 1.04){
    fprintf(fp, "%e %e\n", L, rate.prob(p, L, sbl_iso_fract[BUGEY4]));
  }
  exit(0);
  return;
}
  

}
