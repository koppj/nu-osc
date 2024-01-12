#include "definitions.h"

namespace ns_reactor
{
//#ifdef USE_GAL

extern Fit fit;

#define N_COEF_GAL 101
#define N_INT_GAL 500

struct point
{
  int n_lines;
  double E_line[4];
  double br_sig[4];    // branching times cross section
  double z_0;          // location of source
  double R;            // radius of zylinder
  double H;            // height of zylinder
  double Rexp;         // ratio of observed over expected  
  double sig;          // experimental error on Rexp
  double coef[N_COEF_GAL];
  double coef_no_osc;
} gal[NBIN_GAL];


const double log_dmq_E_min = log(0.05/0.9);
const double log_dmq_E_max = log(50./0.4);
const double d_log_dmq_E = (log_dmq_E_max - log_dmq_E_min)/(N_COEF_GAL-1.);

void calc_coef_gal(void);
void gallium_anomaly(void);


/* constants for the cross sections */

// reference cross sections 
#define SIG_BAHCALL_CR 58.1     // eq. 12 of Bahcall hep-ph/9710491
#define SIG_BAHCALL_AR 70.0     // eq. 43 of Bahcall hep-ph/9710491

// ground state cross sections
#define SIG_GS_CR (SIG_BAHCALL_CR*0.95)  // eq. 13 of Bahcall hep-ph/9710491
#define SIG_GS_AR 66.2                   // eq. 44 of Bahcall hep-ph/9710491

// phase space factors for excited states eq. 14 of Bahcall hep-ph/9710491
#define PS_175_CR 0.669 
#define PS_500_CR 0.220

// eq. 44 of Bahcall hep-ph/9710491
#define PS_175_AR (46.0/66.2)
#define PS_500_AR (17.4/66.2)

// measured values of the excited state contriutions from
// Frekers et al, PLB 706(2011)134, table 3
#define BGT_FRAC_175     (0.34/8.52)
#define BGT_FRAC_175_ERR (0.26/8.52)
#define BGT_FRAC_500     (1.76/8.52)
#define BGT_FRAC_500_ERR (0.14/8.52)

// corrected cross section
#define SIG_COR_CR (SIG_GS_CR*(1.+PS_175_CR*BGT_FRAC_175+PS_500_CR*BGT_FRAC_500))
#define SIG_COR_AR (SIG_GS_AR*(1.+PS_175_AR*BGT_FRAC_175+PS_500_AR*BGT_FRAC_500))

/******************************
 *     initialize
 ******************************/

void gallium_init(void)
{
  // data from tabs 1 and 2 of 0711.4222

  // the three Cr exps
  for(int i = 0; i < 3; i++){
    gal[i].n_lines = 4;

    gal[i].E_line[0] = 0.747;
    gal[i].E_line[1] = 0.752;
    gal[i].E_line[2] = 0.427;
    gal[i].E_line[3] = 0.432;

    gal[i].br_sig[0] = 0.8163 * 60.8;
    gal[i].br_sig[1] = 0.0849 * 61.5;
    gal[i].br_sig[2] = 0.0895 * 26.7;
    gal[i].br_sig[3] = 0.0093 * 27.1;
  }

  // the Ar exp
  gal[3].n_lines = 2;

  gal[3].E_line[0] = 0.811;
  gal[3].E_line[1] = 0.813;

  gal[3].br_sig[0] = 0.902 * 70.1;
  gal[3].br_sig[1] = 0.098 * 70.3;

  // Gallex
  gal[0].R = gal[1].R = 1.9;
  gal[0].H = gal[1].H = 5.0;
  gal[0].z_0 = 2.7;
  gal[1].z_0 = 2.38;

  // SAGE
  gal[2].R = gal[3].R = 0.7;
  gal[2].H = gal[3].H = 1.47;
  gal[2].z_0 = gal[3].z_0 = 0.72;

  // the observed ratio relative to Bahcall cross sect
  // eqs 1-4 from 1006.3244
  gal[0].Rexp = 0.953;
  gal[0].sig = 0.11;

  gal[1].Rexp = 0.812;
  gal[1].sig = 0.10;

  gal[2].Rexp = 0.95;
  gal[2].sig = 0.12;

  gal[3].Rexp = 0.791;
  gal[3].sig = 0.084;

  // correcting Cr experiments
  for(int i = 0; i < 3; i++){
    gal[i].Rexp *= SIG_BAHCALL_CR / SIG_COR_CR;
    gal[i].sig  *= SIG_BAHCALL_CR / SIG_COR_CR;
  }

  // correcting Ar values 
  gal[3].Rexp *= SIG_BAHCALL_AR / SIG_COR_AR;
  gal[3].sig  *= SIG_BAHCALL_AR / SIG_COR_AR;

  //fprintf(stderr, "gallium correction factors for cr = %f and ar = %f\n",
  //	  SIG_BAHCALL_CR / SIG_COR_CR, SIG_BAHCALL_AR / SIG_COR_AR);

  // set data and error in fit
  for(int i = 0; i < NBIN_GAL; i++){

    const int ii = fit.first_bin[GAL] + i;

    fit.Data[ii] = gal[i].Rexp;
    fit.S_data[ii][ii] = norm(gal[i].sig);
  }

  // errors on excited state contributions
  int p = fit.first_pull[GAL]; 
  fit.S_pull[p][p] = norm(BGT_FRAC_175_ERR);
  p++;
  fit.S_pull[p][p] = norm(BGT_FRAC_500_ERR);

#ifndef CHECK_GAL
# ifndef NO_GAL
  calc_coef_gal();
# endif
#else
  gallium_anomaly();
#endif
  return;
}


#ifndef CHECK_GAL

/**************************************
 *     the routine called by fit.chisq
 ***************************************/

void set_table_gallium(Param_5nu &p, double cff[NBIN_CHISQ][NPULLS+1])
{
  for(int k = 0; k < NBIN_GAL; k++){

    // calculate the prediction
    double P = 0.;

    for(int l = 0; l < gal[k].n_lines; l++){

      for(int i = 0; i < N_NU-1; i++){
	for(int j = i+1; j < N_NU; j++){
       
	  const double log_dmq_E = log(fabs(p.Dmq(j,i)) / gal[k].E_line[l]);
	  if(log_dmq_E > log_dmq_E_min){

	    if(log_dmq_E >= log_dmq_E_max) 
	      error("set_table_gal: dmq_E too big");

	    const int ind = int((log_dmq_E - log_dmq_E_min) / d_log_dmq_E);
	    if(ind < 0 || ind >= N_COEF_GAL-1)
	      error("set_table_gal: index out of range");

	    const double x = log_dmq_E - ind * d_log_dmq_E - log_dmq_E_min;
	    const double w = gal[k].coef[ind] + 
	      (gal[k].coef[ind+1] - gal[k].coef[ind]) * x / d_log_dmq_E;

	    P += 4. * norm(p.Ue[i] * p.Ue[j]) * w * gal[k].br_sig[l];
	  }
	}
      }
    }

    // set the table
    const int b = fit.first_bin[GAL] + k;

    // the prediction
    cff[b][NPULLS] = (gal[k].coef_no_osc - P) / gal[k].coef_no_osc;

    if(k < NBIN_GAL-1){
      // cross section error for Cr experiments
      cff[b][fit.first_pull[GAL]]   = cff[b][NPULLS] * PS_175_CR;
      cff[b][fit.first_pull[GAL]+1] = cff[b][NPULLS] * PS_500_CR;
    }else{
      // cross section error for Ar exp
      cff[b][fit.first_pull[GAL]]   = cff[b][NPULLS] * PS_175_AR;
      cff[b][fit.first_pull[GAL]+1] = cff[b][NPULLS] * PS_500_AR;
    }
  }		       
  return;
}

#else // CHECK_GAL

void set_table_gallium(Param_5nu &p, double cff[NBIN_CHISQ][NPULLS+1])
{
  for(int k = 0; k < NBIN_GAL; k++){
    // set the table
    const int b = fit.first_bin[GAL] + k;

    // the prediction
    cff[b][NPULLS] = p.Ue[4];

    if(k < NBIN_GAL-1){
      // cross section error for Cr experiments
      cff[b][fit.first_pull[GAL]]   = cff[b][NPULLS] * PS_175_CR;
      cff[b][fit.first_pull[GAL]+1] = cff[b][NPULLS] * PS_500_CR;
    }else{
      // cross section error for Ar exp
      cff[b][fit.first_pull[GAL]]   = cff[b][NPULLS] * PS_175_AR;
      cff[b][fit.first_pull[GAL]+1] = cff[b][NPULLS] * PS_500_AR;
    }
  }		       
  return;
}

void gallium_anomaly(void)
{
  fit.invert_S();
  //fit.pull_status[fit.first_pull[GAL]] = FIXED;
  //fit.pull_status[fit.first_pull[GAL]+1] = FIXED;

  Param_5nu p;
  double fmin = 0., min = 1000.;
  for(int i = 0; i < 101; i++){
    p.Ue[4] = 0.5 + 0.5 * i / 100.;
    double c = fit.chisq(p);
    printf("%f %f\n", p.Ue[4], c);
    if(c < min){
      min = c;
      fmin = p.Ue[4];
    }
  }
  p.Ue[4] = 1.;
  fprintf(stderr, "min = %f, fmin = %f, delta = %f\n",
	  min, fmin, fit.chisq(p) - min); 
  exit(0);
  return;
}
#endif // CHECK_GAL

/******************************
 *      calc the coefficients *
 ******************************/

void calc_coef_gal(void)
{
  for(int j = 0; j < NBIN_GAL; j++){

    for(int i = 0; i < N_COEF_GAL; i++) gal[j].coef[i] = 0.;
    gal[j].coef_no_osc = 0.;

    // integration over the volume
    for(int kz = 0; kz < N_INT_GAL; kz++){
      const double z = gal[j].H * kz / (N_INT_GAL-1.) - gal[j].z_0;

      for(int kr = 0; kr < N_INT_GAL; kr++){
	const double r = 0.01 + (gal[j].R-0.01) * kr / (N_INT_GAL-1.);
	const double L = sqrt(r*r + z*z);

	for(int i = 0; i < N_COEF_GAL; i++){

	  const double dmq_E = exp(log_dmq_E_min + i * d_log_dmq_E);
          gal[j].coef[i] += r/(L*L) * norm(sin(1.27 * dmq_E * L));
	}
	gal[j].coef_no_osc += r/(L*L);
      }
    }

    double w = 0.;
    for(int i = 0; i < gal[j].n_lines; i++)
      w += gal[j].br_sig[i];
    gal[j].coef_no_osc *= w;
  }
  return;
}

//#endif
}
