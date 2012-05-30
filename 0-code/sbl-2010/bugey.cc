#include "def-reactors.h"

#ifndef BUGEY_TOTAL_RATE

#define DELTAL 3.0
#define ERES_BUG 0.4


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
double cosGem[DM2ANZ][3][25], mittlEnu[3][25], prob_bug[3][25];
double spect_bug[NBIN_BUG][NISO+1][2];
double lInt, eInt, deInt, s2norm, dm2_bugey;
int old_new_bug, iso_bug;


extern Flux flux[NISO][2];  // old and new fluxes for each isotope
extern Fit fit;

/*********************************************************
 *********************************************************/

// calc the probabilities and store in global var prob_bug
void survBug(params p)
{
 for(int i = 0; i < 3; i++){
    for(int k = 0; k < datanzBug[i]; k++){
 
      // interpolating on the table
      double Icos[3];
      for(int j = 0; j < 3; j++){  // j = 41, 51, 54
#ifndef Ip3pI
        const double dmq = (j == 2 ? fabs(p.dmq[I5] - p.dmq[I4]) : p.dmq[j]);
#else
        const double dmq = (j == 2 ? fabs(p.dmq[I5] + p.dmq[I4]) : p.dmq[j]);
#endif

        if(dmq <= 0.)
          Icos[j] = cosGem[0][i][k];

        else{

          double logD = log10(dmq);
#ifdef Ip3pI
          if (logD > LDM2MAX-1.00001*DM_STEP) logD = LDM2MAX-1.00001*DM_STEP; // JK
#endif
	  const int l = (logD < LDM2MIN ? 0 : int( (logD - LDM2MIN)/DM_STEP ) + 1);

	  Icos[j] = cosGem[l][i][k] + 
	    (cosGem[l+1][i][k] - cosGem[l][i][k])/
	    (dm2k(l+1) - dm2k(l)) * (dmq - dm2k(l));
	}
      } // end for j (interpolation)
    
      const double q4 = sqr(p.Ue[I4]);
      const double q5 = sqr(p.Ue[I5]);
        
      prob_bug[i][k] = 
        2.*(1. - q4 - q5) * (q4*Icos[I4] + q5*Icos[I5]);  // 41 and 51     
      prob_bug[i][k] += 2. * q4 * q5 * Icos[2];           // 54
      prob_bug[i][k] += (sqr(1. - q4 - q5) + sqr(q4) + sqr(q5) );    
    }
  } 
  return;
}


/**************** calc the table for class fit ***************/

void set_table_bugey(params &prm, double cff[NBIN_CHISQ][NPULLS+1])
{
  survBug(prm);
  
  const int on = old_new_bug;
  int b = 0;
  for(int i = 0; i < 3; i++){
    for(int k = 0; k < datanzBug[i]; k++){
      
      const int bb = BIN_BUG_0 + b;

      // the prediction
      cff[bb][NPULLS] = cff[bb][FLUX_NORM] = prob_bug[i][k];

      // the pulls
      cff[bb][PULL_U235] = prob_bug[i][k] * spect_bug[b][U235][on] / spect_bug[b][NISO][on];
      cff[bb][PULL_U238] = prob_bug[i][k] * spect_bug[b][U238][on] / spect_bug[b][NISO][on];
      cff[bb][PULL_P239] = prob_bug[i][k] * spect_bug[b][P239][on] / spect_bug[b][NISO][on];
      cff[bb][PULL_P241] = prob_bug[i][k] * spect_bug[b][P241][on] / spect_bug[b][NISO][on];

      cff[bb][E_SCALE_BUG] = (mittlEnu[i][k]-2.8) * prob_bug[i][k];
      
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


void bugey_init(const int old_new)
{
  double datenBug[3][25], sigSqInv[3][25];
  int i,j;

  // some initialization
  for (i=0;i<datanzBug[0];i++){
    datenBug[0][i] = daten115[i]*0.017964 + 1.0;
    sigSqInv[0][i] = 1./sqr(sigDat115[i]*0.008982);
  }
  for (i=0;i<datanzBug[1];i++){
    datenBug[1][i] = daten140[i]*0.017778 + 1.0;
    sigSqInv[1][i] = 1./sqr(sigDat140[i]*0.008889);
  }
  for (i=0;i<datanzBug[2];i++){
    datenBug[2][i] = daten195[i]*0.139304 + 1.0;
    sigSqInv[2][i] = 1./sqr(sigDat195[i]*0.069652);
  }
  for (j=0;j<3;j++){
    for (i=0;i<datanzBug[j];i++){
      mittlEnu[j][i]= 1.8 + 1.0 + deltaE[j]/2 + i*deltaE[j];
    }
  }
  
  int b = BIN_BUG_0;
  for(int i = 0; i < 3; i++){
    for (int j = 0; j < datanzBug[i]; j++){
      fit.Data[b] = datenBug[i][j];
      b++;
    }
  }
    
  
  // calc no-osc spectra for the isotopes
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
      }
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
	fit.Data[BIN_BUG_0 + b] *= spect_bug[b][NISO][OLD] / spect_bug[b][NISO][NEW];
	b++;
      }
    }
  }  

  // systematic uncertainties

  // fully correlated error (subtracting flux error)
  for(int i = BIN_BUG_0; i < BIN_BUG_0 + NBIN_BUG; i++)
    for(int j = BIN_BUG_0; j < BIN_BUG_0 + NBIN_BUG; j++)
      fit.S_data[i][j] = norm(0.048) - norm(0.027);
  
  // the individual detectors
  for(int i = BIN_BUG_0; i < BIN_BUG_0 + 25; i++)
    for(int j = BIN_BUG_0; j < BIN_BUG_0 + 25; j++)
      fit.S_data[i][j] += 2.e-4;

  for(int i = BIN_BUG_0 + 25; i < BIN_BUG_0 + 50; i++)
    for(int j = BIN_BUG_0 + 25; j < BIN_BUG_0 + 50; j++)
      fit.S_data[i][j] += 2.e-4;

  for(int i = BIN_BUG_0 + 50; i < BIN_BUG_0 + 60; i++)
    for(int j = BIN_BUG_0 + 50; j < BIN_BUG_0 + 60; j++)
      fit.S_data[i][j] += 2.e-4;

  fit.S_pull[E_SCALE_BUG][E_SCALE_BUG] = norm(0.02);
  
  // statistical error
  b = BIN_BUG_0;
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
  const double iso_fract[NISO] = {0.538, 0.078, 0.328, 0.056};  
  const double EposKin = eNu - ME - DELTA;
  
  return iso_fract[iso_bug] * flux[iso_bug][old_new_bug].f(eNu) * 
    crossSect(EposKin)*intFuncBug(eNu);
}


#undef DELTAL
#undef ERES_BUG

#endif
