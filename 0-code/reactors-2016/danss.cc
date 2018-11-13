#include <gsl/gsl_sf_expint.h>
#include "definitions.h"

/*****************************************************
 * DANSS results shown at Moriond 2017
 * https://indico.in2p3.fr/event/13763/
 *****************************************************/

namespace ns_reactor
{

#ifdef USE_DANSS

#define DANSS_CORE_HEIGHT  (3.5+0.5) // [meter] height of reactor core
                                     //           + approx. height of detector

extern Fit fit;

// Energy resolution based on fig. 5 in https://arxiv.org/abs/1412.0817
// e is given in MeV
inline double e_res_danss(double e){
  double sqrt_e = sqrt(e);
  return (0.320403 + 0.0221976*e - 0.1345*sqrt_e + 0.0367406/sqrt_e) * e;
}


/********** global vars *******************/
int this_bin_danss;     // current bin, needed during integration of probabilities
double this_DmqL_danss; // current \Delta m^2 * L value
double this_L_danss;    // current detector position
int this_old_new_danss; // old or new fluxes?

// binning in positron energy
static double bin[NBIN_DANSS][2];

// Baselines of the detector at the lower and upper positions [meters]
const double L_danss[] = { 10.7+0.2, 12.7+0.15  };
        // Baselines slightly increased to account for the average horizontal
        // distance between the neutrino production and detection points.
        // Possible, the distribution of neutrino emission from the core could
        // also be asymmetric. Numbers here were obtained by tuning to match
        // the predicted lower/upper ration from slide 10 of the Moriond 2017 talk
enum DANSS_BASELINES { DANSS_UPPER, DANSS_LOWER, DANSS_N_BL };

// Isotopic composition of DANSS reactor, from Moriond 2017 talk
const double isofract_danss[NISO] = { 0.58,  0.30,  0.07,   0.05 };
                                   // U-235, U-238, Pu-239, Pu-241

/**************** class probability ***************/

class DANSS_PROB{
 private:
  double log_DmqL_min, log_DmqL_max;
  double DmqL_min, DmqL_max;
  double Del;

  double sinq_av[NBIN_DANSS][N_RATE_COEF][DANSS_N_BL];
  double sinq(const double DmqL, const int bin, const int bl); // interpolating on table

 public:
  void calc_integrals(void);
  double P(Param_5nu &prm, const int det, const int bin);
} DANSS_prob;


double DANSS_PROB::P(Param_5nu &p, const int bl, const int _bin)
{
  const double L1 = L_danss[bl] - 0.5*DANSS_CORE_HEIGHT;
  const double L2 = L_danss[bl] + 0.5*DANSS_CORE_HEIGHT;
  double P = (1./DANSS_CORE_HEIGHT) * (1./L1 - 1./L2); // \int_{L1}^{L2} dL 1/L^2 / (L2-L1)
//  double P = 1.; //FIXME
  for(int i = 0; i < N_NU-1; i++){
    for(int j = i+1; j < N_NU; j++){
      const double DmqL = fabs(p.Dmq(j,i)) * L_danss[bl];
      if(DmqL > DmqL_min)
        P -= 4. * norm(p.Ue[i] * p.Ue[j]) * sinq(DmqL, bl, _bin);
    }
  }
  return P;
//  return P / (L_danss[bl]*L_danss[bl]); //FIXME

//  double P = 1.; //FIXME FIXME
//  for(int i = 0; i < N_NU-1; i++){
//    for(int j = i+1; j < N_NU; j++){
//      const double DmqL = fabs(p.Dmq(j,i)) * L_danss[bl];
//      if(DmqL > DmqL_min)
//      {
//        double E = 0.5 * (bin[_bin][0] + bin[_bin][1]);
//        P -= 4. * norm(p.Ue[i] * p.Ue[j]) * sin(DmqL/(4.*E)) * sin(DmqL/(4.*E));
//      }
//    }
//  }
//  return P / (L_danss[bl]*L_danss[bl]);
}


// interpolating on the table
double DANSS_PROB::sinq(const double DmqL, const int bl, const int bin)
{
  if(DmqL <  DmqL_min) error("DANSS_PROB::sinq: DmqL too small");
  if(DmqL >= DmqL_max) error("DANSS_PROB::sinq: DmqL too large");

  const double logD = log10(DmqL);

  int i = int( (logD - log_DmqL_min) / Del);
  if(i < 0 || i >= N_RATE_COEF-1) error("DANSS_PROB::sinq: index out of range");

  const double x = (logD - i * Del - log_DmqL_min) / Del;
  if(x < 0. || x > 1.) error("DANSS_PROB::sinq: x out of range");

  return sinq_av[bin][i][bl] + (sinq_av[bin][i+1][bl] - sinq_av[bin][i][bl]) * x;
}


/*************************** init *****************************/

void danss_init(const int old_new)
{
  this_old_new_danss = old_new;

  // set bin edges in MeV
  for(int i = 0; i < NBIN_DANSS; i++){
    bin[i][0] = 1.0 + i * 0.2;
    bin[i][1] = bin[i][0] + 0.2;
  }

  // read data based on slide 10 of Moriond 2017 talk
  double data[NBIN_DANSS][2];       // data and statistically dominated error
  double err_up, err_lo;

  // data file
  FILE *fp_data = fopen(REACTOR_PATH"Data_DANSS/danss-up-down-moriond2017.dat","r");
  if (!fp_data)  {
    error("danss_init: cannot open DANSS data file.");
    return;
  }

  for(int i = 0; i < NBIN_DANSS; i++){
    fscanf(fp_data, "%*f %lf %lf %lf", &data[i][0], &err_lo, &err_up);
    data[i][1] = 0.5 * (err_up - err_lo);
  }
  fclose(fp_data);

  // calc coefficients
  DANSS_prob.calc_integrals();

  // set data and covariance matrix in class fit
  for(int i = 0; i < NBIN_DANSS; i++){
    const int ii = fit.first_bin[DANSS] + i;
    fit.Data[ii] = data[i][0];
    fit.S_data[ii][ii] = norm(data[i][1]);
  }

  // set the pull errors
  const int p0 = fit.first_pull[DANSS];
  for(int p = 0; p < NPULL_DANSS; p++){
    fit.S_pull[p0+p][p0+p] = 1.;
  }
  fit.S_pull[p0][p0] = norm(0.02);

  return;
}

/**************** calc the table for class fit ***************/

void set_table_danss(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1])
{
  for(int i = 0; i < NBIN_DANSS; i++){
    int b = fit.first_bin[DANSS] + i;
    cff[b][NPULLS] = DANSS_prob.P(prm,DANSS_LOWER,i) / DANSS_prob.P(prm,DANSS_UPPER,i);
  }

  // pulls
  for(int i = 0; i < NBIN_DANSS; i++){
    int b = fit.first_bin[DANSS] + i;
    cff[b][fit.first_pull[DANSS]] = cff[b][NPULLS]; // Relative norm up/down
  }

  return;
}


/********************************************************/
/*    calculation of coefficients                       */
/********************************************************/


double danss_flux(const double Enu)
{
  double w = 0.;
  for(int i = 0; i < NISO; i++)
    w += isofract_danss[i] * global_flux(i, Enu, this_old_new_danss);
  return w;
}


// Compute contribution of neutrinos with energy eNu to the positron
// energy bin given by this_bin_danss
double gauss_int_danss(double eNu)
{
   double T_min, T_max;
   double T_e = eNu - DELTA - ME; // DANSS reports positron *kinetic* energy

   //correction motivated by migration matrix from Guillaume
   //T_e += 0.1/3*(T_e - 5.);

   double sigma = e_res_danss(T_e);

   if(this_bin_danss == -1){
     T_min = bin[0][0];
     T_max = bin[NBIN_DANSS-1][1];
   }else{
     T_min = bin[this_bin_danss][0];
     T_max = bin[this_bin_danss][1];
   }

   double x_max = (T_max - T_e) / (sigma*M_SQRT2);
   double x_min = (T_min - T_e) / (sigma*M_SQRT2);

   double res = (erf(x_max) - erf(x_min)) / 2.;
   if(res < 0.) {
      fprintf(stderr, "BLOODY HELL!!\n");
      exit(1);
   }
   return res;
}


// The integrand of the integral over true neutrino energy eNu
double fold_sinq_danss(double eNu)
{
  if (this_DmqL_danss <= 0.)
    return 0.0;
  else
  {
//    const double EposKin = eNu - ME - DELTA;
//    const double a = 0.5 * 2.538 * this_DmqL_danss / eNu;
//    return sin(a)*sin(a) * danss_flux(eNu) * crossSect(EposKin) * gauss_int_danss(eNu);

//    const double EposKin = eNu - ME - DELTA;
//    const double a1 = 2.538 * this_DmqL_danss / eNu * 0.95;
//    const double a2 = 2.538 * this_DmqL_danss / eNu * 1.05;
//    double Cos = (sin(a2) - sin(a1))/(a2 - a1);
//    return 0.5 * (1.-Cos) * danss_flux(eNu) * crossSect(EposKin) * gauss_int_danss(eNu);//FIXME

    // Analytically integrate over assumed size of neutrino production region
    const double EposKin = eNu - ME - DELTA;
    const double dmsq = this_DmqL_danss / this_L_danss;
    const double L1   = this_L_danss - 0.5*DANSS_CORE_HEIGHT;
    const double L2   = this_L_danss + 0.5*DANSS_CORE_HEIGHT;
    const double phi1 = 2.538 * dmsq * L1 / eNu; // 2.538 = 0.5*eV^2*meter/MeV
    const double phi2 = 2.538 * dmsq * L2 / eNu;
    double s1 = sin(0.5*phi1);
    double s2 = sin(0.5*phi2);
    return (1./DANSS_CORE_HEIGHT) * ( 1.269*dmsq/eNu * (gsl_sf_Si(phi2) - gsl_sf_Si(phi1))
                 + s1*s1 / L1 - s2*s2 / L2 )
             * danss_flux(eNu) * crossSect(EposKin) * gauss_int_danss(eNu);
                               // 1.269 = 0.25* eV^2*meter/MeV
  }
}


// The integrand of the normalization integral over true neutrino energy eNu
double fold_sinq_danss_0(double eNu)
{
  const double EposKin = eNu - ME - DELTA;
  return danss_flux(eNu) * crossSect(EposKin) * gauss_int_danss(eNu);
}


void DANSS_PROB::calc_integrals(void)
{
  // DmqL in units of eV^2 m
  log_DmqL_min = ATM_MIN + log10(10.);   // DANSS baselines are 10.7 m and 12.7 m
  log_DmqL_max = STE_MAX + log10(13.);

  DmqL_min = exp10(log_DmqL_min);
  DmqL_max = exp10(log_DmqL_max);

  Del = (log_DmqL_max - log_DmqL_min) / (N_RATE_COEF - 1.);

  double low,up;

  for(int i = 0; i < NBIN_DANSS; i++){   // Loop over all positron energy bins
    this_bin_danss = i;

    // Range of neutrino energies contributing to current positron energy bin
    low = bin[this_bin_danss][0] + ME + DELTA - 3. * e_res_danss(bin[this_bin_danss][0]);
    if(low<EnuMIN) low = EnuMIN;
    up = bin[this_bin_danss][1] + ME + DELTA + 3. * e_res_danss(bin[this_bin_danss][1]);
    if(up>EnuMAX) up=EnuMAX;

    const double n0 = qromb1(fold_sinq_danss_0, low, up, 1.e-6); // Normalization

    // Compute, for each dmsq*L sampling point, the contributions to
    // the i-th positron energy bin
    for(int j = 0; j < N_RATE_COEF; j++){  // loop over all sampling point in dmsq*L
      this_DmqL_danss = exp10(log_DmqL_min + j * Del);

      for (int k = 0; k < DANSS_N_BL; k++){
        this_L_danss = L_danss[k];
        sinq_av[i][j][k] = qromb1(fold_sinq_danss, low, up, 1.e-6) / n0; // \int dE_nu^true
      } // N_BL
    } // N_RATE_COEFF
  } // NBIN_DANSS
  return;
}


// Print DANSS rates
void plot_danss_pred(Param_5nu &prm)
{    
  double cff[NBIN_CHISQ][NPULLS+1];
  set_table_danss(prm, cff);

  for(int i = 0; i < NBIN_DANSS; i++){
    int b = fit.first_bin[DANSS] + i;
    printf("%e %e\n", bin[i][0], cff[b][NPULLS]);
    printf("%e %e\n", bin[i][1], cff[b][NPULLS]);
  }
  return;
}


#endif // USE_DANSS

} // end namespace

