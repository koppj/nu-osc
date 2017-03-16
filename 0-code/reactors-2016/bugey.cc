#include "definitions.h"

namespace ns_reactor
{

#ifdef USE_BUGEY_SP


// definitions of parameter range for Dmq
#define LDM2MIN -2.0
#define LDM2MAX  2.0 

#define DM2ANZ 101
#define DM_STEP ( (LDM2MAX - LDM2MIN) / (DM2ANZ - 2.) )  // including dmq = 0

inline double dm2k(int k){
  if(k == 0)
    return 0.;

  if(k >= DM2ANZ)
    error("[dm2k]: k out of range\n");

  return pow(10, LDM2MIN + (k-1) * DM_STEP);
}



#define DELTAL 3.0
#define ERES_BUG 0.4

const double iso_fract_bug[NISO] = {0.538, 0.078, 0.328, 0.056};  

const int datanzBug[]={25,25,10};
const double entf[]={15.0, 40.0, 95.0}, deltaE[]={0.2,0.2,0.5};
const double daten115[]={1.0, 0.1, -1.0, -0.75, 0.0, -0.8, -1.5, -1.5, -0.9, -0.6,
			 0.25, -0.1, -2.0, 0.9, 0.0, -1.25, -0.6, -2.75, 0.3, -6,
			 0.9, -4.0, -7.8, 0.0, -5.6},
  sigDat115[]={2.0, 2.0, 2.0, 1.75, 1.75, 1.75, 1.75, 1.75, 1.8, 2.2, 2.5,
	       2.5, 2.25, 2.75, 2.9, 2.75, 3.5, 3.5, 4.0, 4.25, 5.5, 5.75,
	       6.0, 8.25, 8.75},
    daten140[]={-2.25, 1.9, -0.5, 0.5, 0.0, 1.6, -0.9, 1.8, -0.75, -3.1,
                -1.0, 0.7, 0.3, -5.2, -0.7, 1.6, 0.7, -6.1, 2.2, -0.4,
                3.3, -3.1, -4.8, -0.4, -4.7},
      sigDat140[]={7.5, 5.5, 4.2, 3.5, 3.7, 3.5, 3.5, 3.5, 3.5, 3.3, 3.7,
		   4.0, 4.2, 4.1, 5.0, 5.5, 6.0, 6.2, 7.0, 8.0, 8.8,
		   10.0, 12.5, 15.0, 20.0},
	daten195[]={-5.9, -1.3, 1.6, 2.2, -1.5, 1.1, 2.1, 2.2, 1.9, -3.0},
	  sigDat195[]={8.0, 4.0, 3.2, 2.7, 3.0, 4.0, 5.0, 4.2, 11.5, 13.0};


// global variables
double cosGem[DM2ANZ][3][25], mittlEnu[3][25], prob_bug[3][25], sigSqInv[3][25];

double spect_bug[NBIN_BUG_SP][NISO+1][2];
double spect_bug1[NBIN_BUG_SP][NISO];
double spect_bug2[NBIN_BUG_SP][NISO];
double spect_bug_cor[NBIN_BUG_SP];

double lInt, eInt, deInt, s2norm, dm2_bugey;
int old_new_bug, iso_bug;

extern Fit fit;

/*********************************************************
 *********************************************************/

// calc the probabilities and store in global var prob_bug
void survBug(Param_5nu &p)
{
 for(int i = 0; i < 3; i++){
    for(int k = 0; k < datanzBug[i]; k++){
 
      // interpolating on the table
      double Icos[3];
      for(int j = 0; j < 3; j++){  // j = 41, 51, 54
        double dmq;
	switch(j){
	case 0: dmq = fabs(p.dmq[3]); break; 
	case 1: dmq = fabs(p.dmq[4]); break; 
	case 2: dmq = fabs(p.dmq[3] - p.dmq[4]); break; 
	}

        if(dmq <= 0.)
          Icos[j] = cosGem[0][i][k];

        else{

          double logD = log10(dmq);
	  const int l = (logD < LDM2MIN ? 0 : int( (logD - LDM2MIN)/DM_STEP ) + 1);

	  Icos[j] = cosGem[l][i][k] + 
	    (cosGem[l+1][i][k] - cosGem[l][i][k])/
	    (dm2k(l+1) - dm2k(l)) * (dmq - dm2k(l));
	}
      } // end for j (interpolation)
    
      const double q4 = norm(p.Ue[3]);
      const double q5 = norm(p.Ue[4]);
        
      prob_bug[i][k] = 
        2.*(1. - q4 - q5) * (q4*Icos[0] + q5*Icos[1]);  // 41 and 51     
      prob_bug[i][k] += 2. * q4 * q5 * Icos[2];           // 54
      prob_bug[i][k] += (norm(1. - q4 - q5) + norm(q4) + norm(q5) );    
    }
  } 
  return;
}


/**************** calc the table for class fit ***************/

void set_table_bugey(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1])
{
  survBug(prm);
  
  const int on = old_new_bug;
  int b = 0;
  for(int i = 0; i < 3; i++){
    for(int k = 0; k < datanzBug[i]; k++){
      
      const int bb = fit.first_bin[BUG_SP] + b;

      // the prediction
      cff[bb][NPULLS] = cff[bb][FLUX_NORM] = prob_bug[i][k];

      // the pulls
      cff[bb][PULL_U235_0] = prob_bug[i][k] * spect_bug[b][U235][on] / spect_bug[b][NISO][on];
      cff[bb][PULL_P239_0] = prob_bug[i][k] * spect_bug[b][P239][on] / spect_bug[b][NISO][on];
      cff[bb][PULL_P241_0] = prob_bug[i][k] * spect_bug[b][P241][on] / spect_bug[b][NISO][on];

      cff[bb][PULL_U235_1] = prob_bug[i][k] * spect_bug1[b][U235] / spect_bug[b][NISO][on];
      cff[bb][PULL_P239_1] = prob_bug[i][k] * spect_bug1[b][P239] / spect_bug[b][NISO][on];
      cff[bb][PULL_P241_1] = prob_bug[i][k] * spect_bug1[b][P241] / spect_bug[b][NISO][on];

      cff[bb][PULL_U235_2] = prob_bug[i][k] * spect_bug2[b][U235] / spect_bug[b][NISO][on];
      cff[bb][PULL_P239_2] = prob_bug[i][k] * spect_bug2[b][P239] / spect_bug[b][NISO][on];
      cff[bb][PULL_P241_2] = prob_bug[i][k] * spect_bug2[b][P241] / spect_bug[b][NISO][on];

      cff[bb][FLUX_COR] = prob_bug[i][k] * spect_bug_cor[b] / spect_bug[b][NISO][on];

      cff[bb][PULL_U238] = prob_bug[i][k] * spect_bug[b][U238][on] / spect_bug[b][NISO][on];

      // energy scale 
      cff[bb][fit.first_pull[BUG_SP]] = (mittlEnu[i][k]-2.8) * prob_bug[i][k];
      
      b++;
     }
   }  
  return;
}


/*********************************************************
 *    initialization
 *********************************************************/

void mittelBugey(int dmq_i);
double mittelLBug(double);
double intFuncBug(double);
double fold_bug(double eNu);
double fold_bug1(double eNu);
double fold_bug2(double eNu);
double fold_bug_cor(double eNu);


void bugey_init(const int old_new)
{
  double datenBug[3][25];
  int i,j;

  // some initialization
  for (i=0;i<datanzBug[0];i++){
    datenBug[0][i] = daten115[i]*0.017964 + 1.0;
    sigSqInv[0][i] = 1./norm(sigDat115[i]*0.008982);
  }
  for (i=0;i<datanzBug[1];i++){
    datenBug[1][i] = daten140[i]*0.017778 + 1.0;
    sigSqInv[1][i] = 1./norm(sigDat140[i]*0.008889);
  }
  for (i=0;i<datanzBug[2];i++){
    datenBug[2][i] = daten195[i]*0.139304 + 1.0;
    sigSqInv[2][i] = 1./norm(sigDat195[i]*0.069652);
  }
  for (j=0;j<3;j++){
    for (i=0;i<datanzBug[j];i++){
      mittlEnu[j][i]= 1.8 + 1.0 + deltaE[j]/2 + i*deltaE[j];
    }
  }
  
  int b = fit.first_bin[BUG_SP];
  for(int i = 0; i < 3; i++){
    for (int j = 0; j < datanzBug[i]; j++){
      fit.Data[b] = datenBug[i][j];
      b++;
    }
  }
    
  
  // calc no-osc spectra for the isotopes
  // and coeffs for flux pulls
  b = 0;
  dm2_bugey = 0.;
    
  for (int i = 0; i < 3; i++){

    lInt = entf[i];
    deInt = deltaE[i];

    for (int j = 0; j < datanzBug[i]; j++){

      eInt = mittlEnu[i][j];
      const double e1 = eInt - deInt / 2. - 3. * ERES_BUG;
      const double e2 = eInt + deInt / 2. + 3. * ERES_BUG;
      
      spect_bug[b][NISO][OLD] = spect_bug[b][NISO][NEW] = 0.;
      for(iso_bug = 0; iso_bug < NISO; iso_bug++){

	old_new_bug = OLD;
        spect_bug[b][iso_bug][OLD] = qromb1(fold_bug, e1, e2, 1.0e-4);
	spect_bug[b][NISO][OLD] += spect_bug[b][iso_bug][OLD];

	old_new_bug = NEW;
        spect_bug[b][iso_bug][NEW] = qromb1(fold_bug, e1, e2, 1.0e-4);
	spect_bug[b][NISO][NEW] += spect_bug[b][iso_bug][NEW];

	if(iso_bug != U238)
        {
	  spect_bug1[b][iso_bug] = qromb1(fold_bug1, e1, e2, 1.e-4);
	  spect_bug2[b][iso_bug] = qromb1(fold_bug2, e1, e2, 1.e-4);
	}
      }
      spect_bug_cor[b] = qromb1(fold_bug_cor, e1, e2, 1.e-4);

      b++;
    }
  }
  old_new_bug = old_new;
  
  // rescale the data to new fluxes
  if(old_new == NEW) 
  {
    b = 0;
    for (int i = 0; i < 3; i++){
      for (int j = 0; j < datanzBug[i]; j++){
	fit.Data[fit.first_bin[BUG_SP] + b] *= spect_bug[b][NISO][OLD] / spect_bug[b][NISO][NEW];

	/* begin check rescaling  
        double s_new = 0., s_old = 0.;

        for(int is = 0; is < NISO; is++){
          s_new += iso_fract_bug[is] * flux[is][NEW].f(mittlEnu[i][j]);
          s_old += iso_fract_bug[is] * flux[is][OLD].f(mittlEnu[i][j]);
        }
        printf("%d %e %e \n", b, spect_bug[b][NISO][OLD] / spect_bug[b][NISO][NEW], s_old/s_new);
	 end check */

	b++;
      }
    }
  }  

  // systematic uncertainties
  const int b0 = fit.first_bin[BUG_SP];

  // fully correlated error (subtracting flux error)
  for(int i = b0; i < b0 + NBIN_BUG_SP; i++)
    for(int j = b0; j < b0 + NBIN_BUG_SP; j++)
      fit.S_data[i][j] = norm(0.039);
  
  // the individual detectors
  for(int i = b0; i < b0 + 25; i++)
    for(int j = b0; j < b0 + 25; j++)
      fit.S_data[i][j] += 4.e-4;

  for(int i = b0 + 25; i < b0 + 50; i++)
    for(int j = b0 + 25; j < b0 + 50; j++)
      fit.S_data[i][j] += 4.e-4;

  for(int i = b0 + 50; i < b0 + 60; i++)
    for(int j = b0 + 50; j < b0 + 60; j++)
      fit.S_data[i][j] += 4.e-4;

  // energy scale error
  fit.S_pull[fit.first_pull[BUG_SP]][fit.first_pull[BUG_SP]] = norm(0.02);
  
  // statistical error
  b = fit.first_bin[BUG_SP];
  for(int i = 0; i < 3; i++){
    for (int j = 0; j < datanzBug[i]; j++){
      fit.S_data[b][b] += 1./sigSqInv[i][j];
      b++;
    }
  }
  
  // calc the integrals
  for(i = 0; i < DM2ANZ; i++)
    mittelBugey(i);

  return;
}

/*************** end init ***************************************/

void mittelBugey(int dmq_i)
{
  extern double eInt, deInt, lInt, dm2_bugey;
  dm2_bugey = dm2k(dmq_i);

  for (int i = 0; i < 3; i++){

    lInt = entf[i];
    deInt = deltaE[i];

    for (int j = 0; j < datanzBug[i]; j++){

      eInt = mittlEnu[i][j];
      const double e1 = eInt - deInt / 2. - 3. * ERES_BUG;
      const double e2 = eInt + deInt / 2. + 3. * ERES_BUG;

      const double int1 = qromb1(intFuncBug, e1, e2, 1.0e-4);
      cosGem[dmq_i][i][j]= int1 / (2. * deltaE[i]);
    }
  }
  return;
}

double intFuncBug(double e)
{
  extern double eInt, deInt;
  double e1,e2,w;

  e1=eInt-deInt/2;
  e2=eInt+deInt/2;
  w = erff((e-e1)/(sqrt(2.)*ERES_BUG)) - erff((e-e2)/(sqrt(2.)*ERES_BUG));

  return mittelLBug(e) * w;
}

        
double mittelLBug(double e)
{
  extern double lInt, dm2_bugey;

  if(dm2_bugey == 0.)
    return 1.;

  const double a = 2.53 * dm2_bugey / e;
  const double l1 = lInt - DELTAL/2;
  const double l2 = lInt + DELTAL/2;

  return (sin(a * l2) - sin(a * l1))/(a * (l2 - l1));
}


/******************************************************
 * the spectrum for rescaling the data to new fluxes
 ******************************************************/ 


double fold_bug(double eNu)
{
  const double EposKin = eNu - ME - DELTA;
  
  return iso_fract_bug[iso_bug] * global_flux(iso_bug, eNu, old_new_bug) * 
    crossSect(EposKin) * intFuncBug(eNu);
}
double fold_bug1(double eNu)
{
  const double EposKin = eNu - ME - DELTA;
  
  return eNu * iso_fract_bug[iso_bug] * global_flux(iso_bug, eNu, old_new_bug) * 
    crossSect(EposKin) * intFuncBug(eNu);
}
double fold_bug2(double eNu)
{
  const double EposKin = eNu - ME - DELTA;
  
  return eNu * eNu * iso_fract_bug[iso_bug] * global_flux(iso_bug, eNu, old_new_bug) * 
    crossSect(EposKin) * intFuncBug(eNu);
}

double fold_bug_cor(double eNu)
{
  double w = read(REACTOR_PATH"Dat/Patrick-U235-err_cor.dat", eNu) * 
    iso_fract_bug[U235] * global_flux(U235, eNu, NEW);
  w += read(REACTOR_PATH"Dat/Patrick-Pu239-err_cor.dat", eNu) * 
    iso_fract_bug[P239] * global_flux(P239, eNu, NEW);
  w += read(REACTOR_PATH"Dat/Patrick-Pu241-err_cor.dat", eNu) * 
    iso_fract_bug[P241] * global_flux(P241, eNu, NEW);

  const double EposKin=eNu-ME-DELTA;
  return w * crossSect(EposKin) * intFuncBug(eNu);
}

/******************************************************
 * plotting Bugey data
 ******************************************************/ 

void plot_bugey_data(void)
{
  const char *file_name[3] = {REACTOR_PATH"dat.bugey.15", REACTOR_PATH"dat.bugey.40", REACTOR_PATH"dat.bugey.95"};

  int b = 0;
  
  for(int i = 0; i < 3; i++)
  { 
    FILE *fp = fopen(file_name[i], "w");
    
    for(int j = 0; j < datanzBug[i]; j++)
    {
      fprintf(fp, "%e %e %e\n", mittlEnu[i][j], fit.Data[fit.first_bin[BUG_SP] + b], sqrt(1./sigSqInv[i][j]));
  	  b++;
    }
    fclose(fp);
  }
  return;
}
    
#undef DELTAL
#undef ERES_BUG

#endif // USE_BUGEY_SP
}
