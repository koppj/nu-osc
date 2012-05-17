#include "def-reactors.h"

extern Fit  fit;
extern Flux flux[NISO][2];  // old and new fluxes for each isotope

bool use_sbl_react[N_SBL_R];
double sbl_flux_pulls[N_SBL_R][NISO];   // pull couplings
double sbl_cos_delta[N_SBL_R][DM2ANZ];  // averaged cosines
double prob_sbl[N_SBL_R];

int iso_gl_sbl, sbl_r_gl, idm2_sbl;
extern int old_new_main;

double int_func(double e), int_func_osc(double e);

/********* the data *********************/

// tab. II of 1101.2755
// rate_old, rate_new
const double sbl_data[N_SBL_R][2] = {
#ifdef BUGEY_TOTAL_RATE
  {0.988, 0.946},         // BUGEY3-15
  {0.994, 0.952},         // BUGEY3-40
  {0.915, 0.876},         // BUGEY3-95
#endif
  {0.987, 0.943},         // BUGEY4
  {0.985, 0.940},	  // ROVNO   
  {1.018, 0.971}, 	  // GOSGEN_1
  {1.045, 0.997}, 	  // GOSGEN_2
  {0.975, 0.930},	  // GOSGEN_3
  {0.832, 0.800},	  // ILL     
  {1.013, 0.944}, 	  // KRASN_1 
  {1.031, 0.960},	  // KRASN_2 
  {0.989, 0.954}	  // KRASN_3 
};


// columns err and corr from tab. II of 1101.2755
const double sbl_error[N_SBL_R][2] = {
#ifdef BUGEY_TOTAL_RATE
  {  4.8, 4.8  },         // BUGEY3-15
  {  4.9, 4.8  },         // BUGEY3-40
  { 14.1, 4.8  },         // BUGEY3-95
#endif
  {  3.0, 3.0},   // BUGEY4     
  {  3.9, 3.0},   // ROVNO     
  {  6.5, 6.0},   // GOSGEN_1  
  {  6.5, 6.0},   // GOSGEN_2  
  {  7.6, 6.0},	  // GOSGEN_3  
  {  9.5, 6.0},	  // ILL       
  {  5.1, 4.1},   // KRASN_1   
  { 20.3, 4.1},	  // KRASN_2   
  {  4.1, 4.1}	  // KRASN_3           
};
double sbl_err_stat[N_SBL_R], sbl_err_syst[N_SBL_R];


// isotope fractions [1101.2755]
// enum Isotopes {U235, U238, P239, P241, NISO};
// Goesgen numbers from PRD 34, 2636
const double sbl_iso_fract[N_SBL_R][NISO] = {
#ifdef BUGEY_TOTAL_RATE
  {0.538, 0.078, 0.328, 0.056},           // BUGEY3-15
  {0.538, 0.078, 0.328, 0.056},           // BUGEY3-40
  {0.538, 0.078, 0.328, 0.056},           // BUGEY3-95
#endif
  {0.538, 0.078, 0.328, 0.056},           // BUGEY4
  {0.614, 0.074, 0.275, 0.031},           // ROVNO   
  {0.619, 0.067, 0.272, 0.042}, 	  // GOSGEN_1
  {0.584, 0.068, 0.298, 0.050}, 	  // GOSGEN_2
  {0.543, 0.070, 0.329, 0.058}, 	  // GOSGEN_3
  {0.930, 0.000, 0.070, 0.000},	          // ILL 
  {1., 0., 0., 0.}, 	                  // KRASN_1 
  {1., 0., 0., 0.},             	  // KRASN_2 
  {1., 0., 0., 0.}              	  // KRASN_3 
};


// the baselines
const double sbl_L[N_SBL_R] = {
#ifdef BUGEY_TOTAL_RATE
  15.,    // BUGEY3-15
  40.,    // BUGEY3-40
  95.,    // BUGEY3-95
#endif
  14.,    // BUGEY4
  18.,    // ROVNO   
  38.,	  // GOSGEN_1
  45.,	  // GOSGEN_2
  65.,	  // GOSGEN_3
  9.,     // ILL 
  33.,    // KRASN_1 
  92.,	  // KRASN_2 
  57.	  // KRASN_3 
};


/**********************************************************************/
/*************** routine called by fit.chisq **************************/
/**********************************************************************/


// calc the probabilities and store in global var prob_sbl
void surv_sbl(params p)
{
 for(int i = 0; i < N_SBL_R; i++){
 
      // interpolating on the table
      double Icos[3];
      for(int j = 0; j < 3; j++){  // j = 41, 51, 54
#ifndef Ip3pI
        const double dmq = (j == 2 ? fabs(p.dmq[I5] - p.dmq[I4]) : p.dmq[j]);
#else
        const double dmq = (j == 2 ? fabs(p.dmq[I5] + p.dmq[I4]) : p.dmq[j]);
#endif

        if(dmq <= 0.)
          Icos[j] = sbl_cos_delta[i][0];

        else{

          double logD = log10(dmq);
#ifdef Ip3pI
          if (logD > LDM2MAX-1.00001*DM_STEP) logD = LDM2MAX-1.00001*DM_STEP; // JK
#endif
	  const int l = (logD < LDM2MIN ? 0 : int( (logD - LDM2MIN)/DM_STEP ) + 1);

	  Icos[j] = sbl_cos_delta[i][l] + 
	    (sbl_cos_delta[i][l+1] - sbl_cos_delta[i][l])/
	    (dm2k(l+1) - dm2k(l)) * (dmq - dm2k(l));
	}
      } // end for j (interpolation)
    
      const double q4 = sqr(p.Ue[I4]);
      const double q5 = sqr(p.Ue[I5]);
        
      prob_sbl[i] = 
        2.*(1. - q4 - q5) * (q4*Icos[I4] + q5*Icos[I5]);  // 41 and 51     
      prob_sbl[i] += 2. * q4 * q5 * Icos[2];              // 54
      prob_sbl[i] += (sqr(1. - q4 - q5) + sqr(q4) + sqr(q5) );    
  } 
  return;
}

void set_table_sbl(params &prm, double cff[NBIN_CHISQ][NPULLS+1])
{
  surv_sbl(prm);

  for(int i = 0; i < N_SBL_R; i++){
    int ii = NBIN_CHOOZ + i;

    if(use_sbl_react[i]){

      for(int p = 0; p < NISO; p++)
	cff[ii][p] = sbl_flux_pulls[i][p] * prob_sbl[i];

      cff[ii][NPULLS] = cff[ii][FLUX_NORM] = prob_sbl[i];

    }else{

      for(int p = 0; p < NISO; p++)
	cff[ii][p] = 0.;
      cff[ii][FLUX_NORM] = 0.;

      cff[ii][NPULLS] = fit.Data[ii];
    }
  }
  return;
}


/***************************************************/
/**************       initialize       *************/
/***************************************************/

void init_sbl_reactors(int old_new)
{
  for(int i = 0; i < N_SBL_R; i++){
    
    // absolute statistical error
    sbl_err_stat[i] = 0.01 * sbl_data[i][old_new] *
      sqrt(norm(sbl_error[i][0]) - norm(sbl_error[i][1]));
    if(sbl_err_stat[i] < 1.e-3) sbl_err_stat[i] = 1.e-3;
    
    // subtract the flux  error of 2.5% from the syst error
    sbl_err_syst[i] = 0.01 * sbl_data[i][old_new] *
      sqrt(norm(sbl_error[i][1]) - norm(2.5));
  }

  // set data and errors in fit
  for(int i = 0; i < N_SBL_R; i++){
    int j = NBIN_CHOOZ + i;

    fit.Data[j] = sbl_data[i][old_new];
    fit.S_data[j][j] = norm(sbl_err_stat[i]) + norm(sbl_err_syst[i]);
  }
  
#define SBL_COR
#ifdef SBL_COR // correlate syst. errors

  // BUGEY4 and ROVNO 
  
  fit.S_data[NBIN_CHOOZ+BUGEY4][NBIN_CHOOZ+ROVNO] = 
    fit.S_data[NBIN_CHOOZ+ROVNO][NBIN_CHOOZ+BUGEY4] =
    sbl_err_syst[BUGEY4] * sbl_err_syst[ROVNO];
  
  // GOSGEN 
  for(int i = GOSGEN_1; i <= GOSGEN_3; i++){
    int ii = NBIN_CHOOZ + i;

    for(int j = GOSGEN_1; j <= GOSGEN_3; j++){
      int jj = NBIN_CHOOZ + j;

      fit.S_data[ii][jj] = sbl_err_syst[i] * sbl_err_syst[j];
    }
    fit.S_data[ii][ii] += norm(sbl_err_stat[i]);
  }

  // Krasnoyarsk
  for(int i = KRASN_1; i <= KRASN_3; i++){
    int ii = NBIN_CHOOZ + i;

    for(int j = KRASN_1; j <= KRASN_3; j++){
      int jj = NBIN_CHOOZ + j;

      fit.S_data[ii][jj] = sbl_err_syst[i] * sbl_err_syst[j];
    }
    fit.S_data[ii][ii] += norm(sbl_err_stat[i]);
  }
#endif
 
  for(int i = 0; i < N_SBL_R; i++)
  {
#ifndef NO_SBL
    use_sbl_react[i] = true;
#else
    use_sbl_react[i] = false;
#endif
  }

  for(int i = 0; i < N_SBL_R; i++){  

    // set pull couplings for isotopes
    for(int j = 0; j < NISO; j++){
      
      iso_gl_sbl = j;      
      sbl_flux_pulls[i][j] = (sbl_iso_fract[i][j] > 1.e-20 ? 
        sbl_iso_fract[i][j] * qromb1(int_func, EnuMIN, EnuMAX, 1.e-5) : 0.);
    }        

    // averaging the cosine
    sbl_r_gl = i;
    for(idm2_sbl = 0; idm2_sbl < DM2ANZ; idm2_sbl++)
      sbl_cos_delta[i][idm2_sbl] = qromb1(int_func_osc, EnuMIN, EnuMAX, 1.e-5);

    
    // normalize
    double s = 0.;
    for(int j = 0; j < NISO; j++)
      s += sbl_flux_pulls[i][j];

    for(int j = 0; j < NISO; j++)
      sbl_flux_pulls[i][j] /= s;
    for(idm2_sbl = 0; idm2_sbl < DM2ANZ; idm2_sbl++)
      sbl_cos_delta[i][idm2_sbl] /= s;
  }    
  return;
}


/**********************************************************************/

double int_func(double e)
{  
  const double EposKin = e - ME - DELTA;
  return crossSect(EposKin) * flux[iso_gl_sbl][old_new_main].f(e);
}

#define DELTAL 3.0

double int_func_osc(double e)
{
  double ff = 0.;

  for(int i = 0; i < NISO; i++)  
    ff += sbl_iso_fract[sbl_r_gl][i] * flux[i][old_new_main].f(e);

  double Cos;
  if(idm2_sbl == 0)
    Cos = 1.;
  else
  {
    const double a = 2.53 * dm2k(idm2_sbl) / e;
    const double l1 = sbl_L[sbl_r_gl] - DELTAL/2.;
    const double l2 = sbl_L[sbl_r_gl] + DELTAL/2.;

    Cos = (sin(a * l2) - sin(a * l1))/(a * (l2 - l1));
  }

  const double EposKin = e - ME - DELTA;
  return crossSect(EposKin) * Cos * ff;
}

#undef DELTAL

/**********************************************************************/

void calc_reactors(void)
{
  fprintf(stderr, "calculating reactors only\n");
  
  bool incl[NUM_EXP];  
  for(int i = 0; i < NUM_EXP; i++)
    incl[i] = (i == reactor);
  
  params p;
  p.Ue3 = 0.;
  p.Um[I4] = p.Um[I5] = 0.;
  p.delta = 0.;

  fix_params f;  
  f.Ue[I4] = f.Ue[I5] = false; 
  f.Um[I4] = f.Um[I5] = true;
  f.dmq[I4] = f.dmq[I5] = true;
  f.delta = f.Ue3 = true;

  FILE *fp_proj = fopen("dat.reactors-proj1.out","w");
  if(fp_proj == NULL)
    nrerror("[calc_reactors] cannot open file");

  FILE *fp_grid = fopen("dat.reactors-grid1.out","w");
  if(fp_grid == NULL)
    nrerror("[calc_reactors] cannot open file");

#define N_DMQ 161  
 
  // the loop over the dmq's
  for(int i = 0; i < N_DMQ; i++)
  {
    p.dmq[I5] = pow(10., -1. + 2. * i / (N_DMQ-1.));
        
    double min = 1.e4;
    params p_min;

    for(int j = 0; j < N_DMQ; j++)
    {
      p.dmq[I4] = pow(10., -1. + 2. * j / (N_DMQ-1.));
	
      p.Ue[I4] = p.Ue[I5] = 0.1;      
      const double cc = min_chisq(&p, f, incl);

      fprintf(fp_grid, "%e  %e  %e  %e  %e\n", log10(p.dmq[I4]), log10(p.dmq[I5]), cc,
	      p.Ue[I4], p.Ue[I5]);

      if(cc < min){
	min = cc;
	p_min = p;
      }
    }
    
    fprintf(fp_proj, "%e %e\n", p.dmq[I5], min);
    
    fflush(fp_grid);
    fflush(fp_proj);
  }
  fclose(fp_grid);
  fclose(fp_proj);
  exit(0);
}

#undef N_DMQ

/**********************************************************************/

void print_sbl_data(void)
{
  for(int i = 0; i < N_SBL_R; i++){
    printf("%d %e %f %e %e %e\n", i+1, sbl_data[i][OLD], i+1.2, 
	   sbl_data[i][NEW], sbl_err_stat[i], 
	   sqrt(norm(sbl_err_stat[i]) + norm(sbl_err_syst[i])));
  }
  return;
}

void print_sbl_pred(double xi[NPULLS])
{
  for(int i = 0; i < N_SBL_R; i++){
    double pred = 1.;
    for(int j = 0; j < NISO; j++)
      pred += xi[j] * sbl_flux_pulls[i][j];
    
    printf("%f %f\n", i+0.5, pred);
    printf("%f %f\n", i+1.5, pred);
  }
  return;
}
  
