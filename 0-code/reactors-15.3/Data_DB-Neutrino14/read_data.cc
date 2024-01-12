#include <stdio.h>
#include <stdlib.h>
#include <math.h>

inline void error(const char *text) {
   fputs(text, stderr);
   fputc('\n', stderr);
   exit(1);
}

inline double norm(const double x){
  return x*x;
}

#define NBIN_DB 36
#define NDET_DB 8
#define NBG_DB 5
#define NREACT_DB 6

// baselines from table 2 of 1210.6327
const double basel_DB[NDET_DB][NREACT_DB] = {
  {0.362, 0.372, 0.903, 0.817, 1.354, 1.265},   // EH1
  {0.358, 0.368, 0.903, 0.817, 1.354, 1.266},   // EH1
  {1.332, 1.358, 0.468, 0.490, 0.558, 0.499},   // EH2
  {1.348, 1.348, 0.480, 0.480, 0.528, 0.528},   // EH2 inventing numbers!
  {1.920, 1.894, 1.533, 1.534, 1.551, 1.525},   // EH3
  {1.918, 1.892, 1.535, 1.535, 1.555, 1.528},   // EH3
  {1.925, 1.900, 1.539, 1.539, 1.556, 1.530},   // EH3
  {1.912, 1.912, 1.540, 1.540, 1.548, 1.548}    // EH3 inventing numbers!
};

double weight_DB[NDET_DB][NREACT_DB]; 


// ***************** DB8 data ******************
// *** data from slide 16 of NuFact2014 talk ***
// *********************************************

const double IBD_cand_DB8[NDET_DB] = {
  202461., 206217., 193356., 190046., 27067., 27389., 27032., 27419.
};
// days
const double DAQ_DB8[NDET_DB] = {
  374.447, 374.447, 378.407,  378.407, 372.685, 372.685, 372.685, 372.685 
};
// multiplying eps_mu * eps_m
const double eff_DB8[NDET_DB] = {
  0.8255 * 0.9746,
  0.8223 * 0.9749,
  0.8574 * 0.9759,
  0.8577 * 0.9756,
  0.9811 * 0.9762,
  0.9811 * 0.976,
  0.9808 * 0.9757,
  0.9811 * 0.9758
};
 
// bg rate per day
const double bg_DB8_data[NBG_DB][NDET_DB] = {
  {8.62, 8.76, 6.43, 6.86, 1.07, 0.94, 0.94, 1.26},  // accidentals
  {0.78, 0.78, 0.54, 0.54, 0.05, 0.05, 0.05, 0.05},  // fast neutron
  {2.80, 2.80, 1.70, 1.70, 0.27, 0.27, 0.27, 0.27},  // Li, He
  {0.20, 0.21, 0.18, 0.22, 0.06, 0.04, 0.04, 0.07},  // Am-C
  {0.08, 0.07, 0.05, 0.07, 0.05, 0.05, 0.05, 0.05}   // CO
};


// ***************** DB6 data ******************
// *** data from slide 7 of NuFact2013 talk  ***
// *********************************************

const double IBD_cand_DB6[NDET_DB] = {
  101290., 102519., 92912., 0., 13964., 13894., 13731., 0.
};
// days
const double DAQ_DB6[NDET_DB] = {
  191.001, 191.001, 189.645, 0., 189.779, 189.779, 189.779, 0.
};
const double eff_DB6[NDET_DB] = {
  0.7957, 0.7927, 0.8282, 1., 0.9577, 0.9568, 0.9566, 1.
};
 
// bg rate per day
const double bg_DB6_data[NBG_DB][NDET_DB] = {
  {9.54, 9.36, 7.44, 0., 2.96, 2.92, 2.87, 0.},  // accidentals
  {0.92, 0.92, 0.62, 0., 0.04, 0.04, 0.04, 0.},  // fast neutron
  {2.40, 2.40, 1.20, 0., 0.22, 0.22, 0.22, 0.},  // Li, He
  {0.26, 0.26, 0.26, 0., 0.26, 0.26, 0.26, 0.},  // Am-C
  {0.08, 0.07, 0.05, 0., 0.04, 0.04, 0.04, 0.}   // CO
};





int main(void)
{
  // set weight_DB: relat contribution of each react to each detector
  for(int i = 0; i < NDET_DB; i++){

    double w = 0;
    for(int j = 0; j < NREACT_DB; j++)
      w += 1./norm(basel_DB[i][j]);

    for(int j = 0; j < NREACT_DB; j++)
      weight_DB[i][j] = 1./norm(basel_DB[i][j]) / w;
  }

  // reading data
  FILE *fp_data = fopen("DB_observed.txt","r");
  FILE *fp_noos = fopen("DB_no-osc.txt","r");

  double DB_data[NBIN_DB], DB_noosc[NBIN_DB];

  for(int i = 0; i < NBIN_DB; i++){
    if(fscanf(fp_data, "%le", &DB_data[i]) != 1)
      error("cannot read DB data");
    if(fscanf(fp_noos, "%le", &DB_noosc[i]) != 1)
      error("cannot read DB no oscillation data");
  }
  fclose(fp_data);
  fclose(fp_noos);

  // check data points
  double sum = 0.;
  for(int i = 0; i < NBIN_DB; i++)
    sum += DB_data[i];

  printf("sum spectrum = %f\n", sum); 

  double ev6 = 0., ev8 = 0., rate = 0.;

  // sum over far detectors
  for(int i = 4; i < NDET_DB; i++){

    ev6 = IBD_cand_DB6[i];
    ev8 = IBD_cand_DB8[i];

    for(int j = 0; j < NBG_DB; j++){
      ev6 -= bg_DB6_data[j][i] * DAQ_DB6[i] * eff_DB6[i];
      ev8 -= bg_DB8_data[j][i] * DAQ_DB8[i] * eff_DB8[i];
    }
    rate += (ev6 + ev8) / (DAQ_DB6[4] + DAQ_DB8[4]);
  }
  printf("%f\n\n", rate  );


  // check no osc prediction
  sum = 0.;
  for(int i = 0; i < NBIN_DB; i++)
    sum += DB_noosc[i];

  printf("sum no osc = %f\n", sum); 

  ev6 = 0., ev8 = 0., rate = 0.;
  double w_nd = 0., w_fd = 0.;

  // sum over near detectors
  for(int i = 0; i < 4; i++){

    ev6 = IBD_cand_DB6[i];
    ev8 = IBD_cand_DB8[i];

    for(int j = 0; j < NBG_DB; j++){
      ev6 -= bg_DB6_data[j][i] * DAQ_DB6[i] * eff_DB6[i];
      ev8 -= bg_DB8_data[j][i] * DAQ_DB8[i] * eff_DB8[i];
    }
    rate += (ev6 + ev8) / (DAQ_DB6[0] + DAQ_DB8[0]);

    for(int j = 0; j < NREACT_DB; j++){
      w_nd += eff_DB6[i] * DAQ_DB6[i] / norm(basel_DB[i][j]);
      w_nd += eff_DB8[i] * DAQ_DB8[i] / norm(basel_DB[i][j]);

      const int k = i+4;
      w_fd += eff_DB6[k] * DAQ_DB6[k] / norm(basel_DB[k][j]);
      w_fd += eff_DB8[k] * DAQ_DB8[k] / norm(basel_DB[k][j]);
    }
  }

  printf("%f\n", rate * w_fd / w_nd );

  return 0;
}
