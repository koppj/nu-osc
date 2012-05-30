#include "definitions.h"

double dm2_gl_cdhs;

#define LBACK 0.885
#define DLBACK 0.124
#define LFRON 0.13
#define DLFRON 0.074
#define EMUINF 30.0
#define DATANZ 15
#define SIGMNORM2 0.000625

const double rCorr[]={0.985,1.006,0.968,1.148,1.0,1.137,1.155,0.887,
		     1.123,0.973,1.039,1.187,1.196,1.006,0.936};
const double sigmaRcorr[]={0.066,0.055,0.070,0.075,0.087,0.104,0.135,0.2,
			  0.078,0.075,0.085,0.103,0.12,0.132,0.11};

double sR2I[DATANZ]; //  1 / sigma_i^2 

double EmuBin[2][8];
const int binAnz[] = {7, 6};

// some integrals
double I1[2][9], IPback[2][9], IPfron[2][9];

// new integrals
double int_cdhs_back[DM2ANZ][DATANZ], int_cdhs_fron[DM2ANZ][DATANZ];

double norma(double s2);
double bdfP(double s2, int i);
double mittelLCD(double e);
double funcIP(double e);
double funcI1(double e);
double flux_cdhs(double e);

/*****************************************************
 * the chisq
 *****************************************************/

void back_fron_R(params p, double R[DATANZ])
{
  for(int i = 0; i < DATANZ; i++){

    // interpolating on the table
    double Icos_back[3], Icos_fron[3];

    for(int j = 0; j < 2; j++){  // j = 41, 51

      const double dmq = p.dmq[j];

#ifdef CHECK
      if(dmq <= 0.){
        Icos_back[j] = int_cdhs_back[0][i];
        Icos_fron[j] = int_cdhs_fron[0][i];
      }else
#endif
      {
        const double logD = log10(dmq);

	const int k = (logD < LDM2MIN ? 0 : int( (logD - LDM2MIN)/DM_STEP ) + 1);

	Icos_back[j] = int_cdhs_back[k][i] + 
          (int_cdhs_back[k+1][i] - int_cdhs_back[k][i])/
	  (dm2k(k+1) - dm2k(k)) * (dmq - dm2k(k));
	
	Icos_fron[j] = int_cdhs_fron[k][i] + 
          (int_cdhs_fron[k+1][i] - int_cdhs_fron[k][i])/
	  (dm2k(k+1) - dm2k(k)) * (dmq - dm2k(k));
      }
    } // end for j 

    // 54
#ifndef Ip3pI
    const double dmq = fabs(p.dmq[I5] - p.dmq[I4]);
    const double logD = log10(dmq);
#else
    const double dmq = fabs(p.dmq[I5] + p.dmq[I4]);
    double logD = log10(dmq);
    if (logD > LDM2MAX-1.00001*DM_STEP) logD = LDM2MAX-1.00001*DM_STEP; // JK
#endif

    if(dmq == 0.){
        Icos_back[2] = int_cdhs_back[0][i];
        Icos_fron[2] = int_cdhs_fron[0][i];
    }else{

	const int k = (logD < LDM2MIN ? 0 : int( (logD - LDM2MIN)/DM_STEP ) + 1);

	Icos_back[2] = int_cdhs_back[k][i] + 
          (int_cdhs_back[k+1][i] - int_cdhs_back[k][i])/
	  (dm2k(k+1) - dm2k(k)) * (dmq - dm2k(k));
	
	Icos_fron[2] = int_cdhs_fron[k][i] + 
          (int_cdhs_fron[k+1][i] - int_cdhs_fron[k][i])/
	  (dm2k(k+1) - dm2k(k)) * (dmq - dm2k(k));

    } // end interpolation

    const double n0_back = int_cdhs_back[0][i];
    const double n0_fron = int_cdhs_fron[0][i];
    
    const double q4 = sqr(p.Um[I4]);
    const double q5 = sqr(p.Um[I5]);
        
    double P_back = q4 * Icos_back[I4] + q5 * Icos_back[I5]; // 41 and 51
    P_back *= 2. * (1. - q4 - q5);      
    P_back += 2. * q4 * q5 * Icos_back[2];                   // 54
    P_back += (sqr(1. - q4 - q5) + sqr(q4) + sqr(q5) ) * n0_back;
    
    double P_fron = q4 * Icos_fron[I4] + q5 * Icos_fron[I5]; // 41 and 51
    P_fron *= 2. * (1. - q4 - q5);      
    P_fron += 2. * q4 * q5 * Icos_fron[2];                   // 54
    P_fron += (sqr(1. - q4 - q5) + sqr(q4) + sqr(q5) ) * n0_fron;

    R[i] = P_back / P_fron;

  } // end for i (bins)

  return;
}

double chi2cdhs(params p)
{
  double R[DATANZ];  
  back_fron_R(p, R);
  
  // normalization uncertainty
  double sum1 = 0., sum2 = 0.;
  for (int i = 0; i < DATANZ; i++){
    
    sum1 += rCorr[i] * sR2I[i] * R[i];
    sum2 += sR2I[i] * sqr(R[i]);
  }  
  const double a = (sum1 + 1./SIGMNORM2)/(sum2 + 1./SIGMNORM2);

  // the chisq
  double sum = 0.;
  for (int i=0; i < DATANZ; i++)    
    sum += sR2I[i] * sqr(rCorr[i] - a * R[i]);
  
  return sum + sqr(1.-a)/SIGMNORM2;
}



/*****************************************************
 * calc the integrals
 *****************************************************/

double l1, l2, ecom;

void mittelCdhs(int dm2_i)
{
  dm2_gl_cdhs = dm2k(dm2_i);

  for (int i = 0; i < 2; i++){
    for(int j = 0; j <= binAnz[i]; j++){
      
      if(dm2_i == 0)
	IPback[i][j] = IPfron[i][j] = I1[i][j];
      else{
      
        ecom = EmuBin[i][j];
      
        l1 = LBACK - DLBACK/2;
        l2 = LBACK + DLBACK/2;
      
        IPback[i][j] = qromb1(funcIP,ecom,EMUINF,1.0e-6);
	IPback[i][j] -= I1[i][j];
      
        l1 = LFRON - DLFRON/2;
        l2 = LFRON + DLFRON/2;
      
        IPfron[i][j] = qromb1(funcIP,ecom,EMUINF,1.0e-6);
	IPfron[i][j] -= I1[i][j];
      }	
    }
    I1[i][binAnz[i]+1] = IPback[i][binAnz[i]+1] = IPfron[i][binAnz[i]+1] = 0.0;
  }
  return;
}

void calcIntegralCdhs(void)
{
  for(int k = 0; k < DM2ANZ; k++){
    
    mittelCdhs(k);
    
    for(int jj = 0; jj < DATANZ; jj++){
      
      int i = 0, j = jj;      
      if(j >= 8){
        i = 1; 
        j -= 8;
      }      
      
      int_cdhs_back[k][jj] = IPback[i][j] - IPback[i][j+1];
      int_cdhs_fron[k][jj] = IPfron[i][j] - IPfron[i][j+1];
    }
  }
  return;
}
      

double mittelLCD(double e)
{
  extern double l1, l2, dm2_gl_cdhs;
  const double a = 2.53 * dm2_gl_cdhs / e;

  // integration over detector dimensions
  // calc mean cos(aL)
  if(a <= 0.)
    return 1.;
  
  return (sin(a * l2) - sin(a * l1))/(a * (l2 - l1));
}

#define ACS3 0.0464133
double funcIP(double e)
{
  return flux_cdhs(e)*e*(1.-ecom/e+ACS3*(1.-pow(ecom/e,3)))*(1. + mittelLCD(e));
}
double funcI1(double e)
{
  return flux_cdhs(e)*e*(1.-ecom/e+ACS3*(1.-pow(ecom/e,3)));
}
#undef ACS3

double flux_cdhs(double eNu)
{
  return exp(-eNu);
}

#define KAZ 7.9866
#define ME 0.511
#define MMU 105.7
#define I 0.000286         /* Anregungsenergie in MeV  */
#define X0 -0.0012
#define X1 3.1531
#define ADELTA 0.1468
#define MDELTA 2.9632
#define CDELTA -4.291
#define FDELTA 4.6052
#define AFE 55.845           /* Massenzahl A von Fe */
#define RHOFE 7.87             /* Dichte von Eisen in g/cm^3 */
#define COSI 1.1477527     /*  1/(cos(22ø)*cos(20ø))  */


double mdEdx(double bgSq)
{
  double betaSq, tmax, delta, x;

  betaSq=bgSq/(1.0+bgSq);
  tmax= 2.0*ME*bgSq/(1.0+2.0*sqrt(bgSq/betaSq)*ME/MMU+sqr(ME/MMU));
  x=0.5*log10(bgSq);
  delta=FDELTA*x+CDELTA;
  if(x > X0) 
    delta += ADELTA*pow(X1-x,MDELTA);
  else{
    delta = FDELTA*X0 + CDELTA + ADELTA*pow(X1-X0,MDELTA);
    delta *= pow(10,2*x-X0);
  }
  return KAZ/(AFE*betaSq)*(0.5*log(2*ME*bgSq*tmax/sqr(I))-betaSq-0.5*delta);
}


double intFunc(double Emu)
{
  double bgSq;

  bgSq=sqr(Emu/MMU)-1.0;
  return 1/mdEdx(bgSq);
}

double rcom;
double rangeBis(double emu)    /* emu:  Myonenergie in MeV  */
{
  return rcom - qromb1(intFunc,108,emu,1.0e-6);
}

double EmuR(double rcm)
{
  rcom = rcm * RHOFE;
  return 0.001 * rtbis(rangeBis, 600.0, 8000.0, 0.1);
}

    
/*******************************
 *  initialisation
 *******************************/
    
void initCDHS(void)
{
  const double prBin[2][8]={{37.5,50,75,100,137.5,187.5,250,337.5},
			   {40,50,70,100,125,170,225}};
  int i,j;

  for(i=0;i<=1;i++)
    for(j=0;j<=binAnz[i];j++) 
      EmuBin[i][j] = EmuR(prBin[i][j]*COSI);
 
  for(i=0;i<DATANZ;i++) 
    sR2I[i]=1/sqr(sigmaRcorr[i]);
  
  for(i=0; i<2; i++){
    for(j=0;j<=binAnz[i];j++){
      
      ecom = EmuBin[i][j];
      I1[i][j] = qromb1(funcI1,ecom,EMUINF,1.0e-6);
    }
  }
  
  calcIntegralCdhs();
  return;
}

#undef KAZ
#undef ME
#undef MMU
#undef I
#undef X0
#undef X1
#undef ADELTA
#undef MDELTA
#undef CDELTA
#undef FDELTA
#undef RHOFE
#undef AFE
#undef COSI
#undef LBACK
#undef DLBACK
#undef LFRON
#undef DLFRON
#undef DATANZ
#undef SIGMNORM2
#undef EMUINF
