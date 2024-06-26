#include "definitions.h"

namespace ns_reactor
{

#ifdef USE_RENO


// source: RENO talk in Venice, NeuTel 2013 
// updated from WIN2013 talk
// updated from Neutrino2014

#define NREACT_RENO 6

extern Fit fit;
extern Rate_pull_coupl rate_pull_coupl;  // defined in class_flux.cc
extern Rate_coef rate;                   // defined in class_flux.cc
extern int old_new_main;


/************** data **************************/

// assume same values as in DC
const double isofract_RENO[NISO] = {0.488, 0.087, 0.359, 0.076};

// talbe 1.2 of RENO TDR
const double basel_RENO[NBIN_RENO][NREACT_RENO] = {
  {.6679,.4518,.3048,.3361,.5139,.7391},
  {1.556, 1.456, 1.395, 1.381, 1.413, 1.490}
};

// data from slide 10 Neutrino2014
const double IDB_cand[NBIN_RENO] = {433196., 50750.};
const double DAQ_RENO[NBIN_RENO] = {761.11, 794.72};

#define NBG_RENO 3
// bg rate per day
const double bg_RENO_data[NBG_RENO][NBIN_RENO] = {
  { 1.82, 1.98},  // accidentals for ND, Cf for FD (neglect Cf for ND and acc for FD) 
  { 8.28, 1.85},  // Li, He    
  { 2.09, 0.44}   // fast neutron
};
// bg error per day
const double bg_err_RENO[NBG_RENO][NBIN_RENO] = {
  {0.11, 0.27},
  {0.66, 0.2}, 
  {0.06, 0.02}
}; 

// on slide 31 WIN2013 they give R = 0.929 for FD 
// R for ND is not available -  use 0.99 - "guessed" by Concha to reproduce their chisq curve
// R for FD not given in Neutrino2014 talk - assume same ratio as at WIN2013 !!
double ratios_RENO[NBIN_RENO] = {0.99, 0.929}; 

// number of bg events and relative error
double bg_RENO[NBIN_RENO][2];
double no_osc_react_RENO[NBIN_RENO];

// tab I of RENO paper (only in version 1)
double weight_RENO[NBIN_RENO][NREACT_RENO] = {
  {0.0678, 0.1493, 0.3419, 0.2701, 0.1150, 0.0558},
  {0.1373, 0.1574, 0.1809, 0.1856, 0.1780, 0.1608}
};

/**************** calc the table for class fit ***************/

void set_table_RENO(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1])
{  
  for(int i = 0; i < NBIN_RENO; i++){

    const int b = fit.first_bin[RENO] + i;
  
    cff[b][NPULLS] = 0.;

    for(int j = 0; j < NREACT_RENO; j++){

      // reactor norm
      const int pr = fit.first_pull[RENO] + 1 + j;
      cff[b][pr] = no_osc_react_RENO[i] * weight_RENO[i][j] *
	rate.prob(prm, basel_RENO[i][j] * 1.e3, isofract_RENO);

      cff[b][NPULLS] += cff[b][pr];
    }

    // overall normalization of reactor signal
    cff[b][fit.first_pull[RENO]] = cff[b][NPULLS];

    // the fluxes
    cff[b][FLUX_NORM] = cff[b][NPULLS];

    // the pulls - only flux normalizations
    cff[b][PULL_U235_0] = cff[b][NPULLS]*rate_pull_coupl.unc(U235, 0, isofract_RENO);
    cff[b][PULL_P239_0] = cff[b][NPULLS]*rate_pull_coupl.unc(P239, 0, isofract_RENO); 
    cff[b][PULL_P241_0] = cff[b][NPULLS]*rate_pull_coupl.unc(P241, 0, isofract_RENO); 
    cff[b][PULL_U238]   = cff[b][NPULLS]*rate_pull_coupl.unc(U238, 0, isofract_RENO);
    
    cff[b][FLUX_COR]  = rate_pull_coupl.cor(isofract_RENO); // ????

    // uncorr detector normalizations
    const int pd = fit.first_pull[RENO] + 1 + NREACT_RENO + i;
    cff[b][pd] = cff[b][NPULLS];

    // adding background to prediction
    cff[b][NPULLS] += bg_RENO[i][0];
        
    // the backgrounds
    const int pb = fit.first_pull[RENO] + 1 + NREACT_RENO + NBIN_RENO + i;
    cff[b][pb] = bg_RENO[i][0];
  }
  return;
}

/***************** init *************************************/

void RENO_init(void)
{
  // set the background per detector
  for(int i = 0; i < NBIN_RENO; i++){
    bg_RENO[i][0] = bg_RENO[i][1] = 0.;

    for(int j = 0; j < NBG_RENO; j++){
      bg_RENO[i][0] += bg_RENO_data[j][i];
      bg_RENO[i][1] += norm(bg_err_RENO[j][i]);
    }
    bg_RENO[i][1] = sqrt(bg_RENO[i][1]) / bg_RENO[i][0]; // relat. error

    bg_RENO[i][0] *= DAQ_RENO[i];
  }


  // set weight_RENO: relat contribution of each react to each detector
  // if MY_RENO_WEIGHT is not defined use weights from Tab. I of RENO paper
#ifdef MY_RENO_WEIGHT
  // react 1 and 2 off for one month
  const double time_weight[NREACT_RENO] = {1.-30./229., 1.-30./229., 1., 1., 1., 1.};

  for(int i = 0; i < NBIN_RENO; i++){

    double w = 0;
    for(int j = 0; j < NREACT_RENO; j++)
      w += time_weight[j]/norm(basel_RENO[i][j]);

    for(int j = 0; j < NREACT_RENO; j++)
      weight_RENO[i][j] = time_weight[j]/norm(basel_RENO[i][j]) / w;

  }
#endif

  /*
  // calculate the ratio from the observed event rate (bkg subtracted)
  const double rates[NBIN_RENO] = {569.16, 63.86}; //  slide 10 from Neutrino2014
  double w_n = 0., w_f = 0.;
  for(int i = 0; i < NREACT_RENO; i++){
    w_n += 1./norm(basel_RENO[0][i]);
    w_f += 1./norm(basel_RENO[1][i]);
  }

  ratios_RENO[0] = 1.;
  ratios_RENO[1] = rates[1] / rates[0] * w_n / w_f;
  fprintf(stderr, "%f\n", ratios_RENO[1]); exit(0);
  */

  for(int i = 0; i < NBIN_RENO; i++)
    no_osc_react_RENO[i] = (IDB_cand[i] - bg_RENO[i][0]) / ratios_RENO[i];

  /*
  for(int j = 0; j < NREACT_RENO; j++)
    fprintf(stderr, "%f %f\n", weight_RENO[0][j], weight_RENO[1][j]);
  */
  
  // set data and stat. errors in class fit
  for(int i = 0; i < NBIN_RENO; i++){
  
    const int b = fit.first_bin[RENO] + i;  
    fit.Data[b] = fit.S_data[b][b] = IDB_cand[i];    
  }

  /*** set pull errors in class fit ***/
  int p = fit.first_pull[RENO];

  // normalization (free)
  fit.S_pull[p][p] = norm(0.025); // used in RENO paper --> set pull ACTIVE
  fit.pull_status[p] = FREE;      // use free norm (different from RENO paper)

  // reactor norm
  for(int i = 0; i < NREACT_RENO; i++){
    p++;
    fit.S_pull[p][p] = norm(0.009);
  }
  // detector norm
  for(int i = 0; i < NBIN_RENO; i++){
    p++;
    fit.S_pull[p][p] = norm(0.002);
  }
  // backgrounds
  for(int i = 0; i < NBIN_RENO; i++){
    p++;
    fit.S_pull[p][p] = norm(bg_RENO[i][1]);
  }

  return;  
}

#endif // USE_RENO

}
