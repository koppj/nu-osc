#include "definitions.h"

namespace ns_reactor
{

// global classes
Flux flux[NISO][2];  // old and new fluxes for each isotope using polynomials
Fit  fit;
Rate_coef rate;

int old_new_main = NEW;

void plot_fluxes(void);
void calc_cs_per_fission(void);  // in class_flux.cc



/********************
 *    main          *
 ********************/

int reactor_main(void)
{
  fprintf(stderr, "For each experiment to be included define USE_ex with ex = SBL, CHOOZ, PV, KAML, DC, DB, RENO, BUGEY_SP\n");  

  /**** initialization *****/

  init_fluxes();
  rate.init();
  
  //calc_cs_per_fission(); return 0;
  //plot_fluxes(); return 0;

  int n_data_points = NBIN_CHISQ;
  
  fprintf(stderr, "using: ");
#ifdef USE_SBL
  const int n_sbl = init_sbl_reactors(old_new_main);
  n_data_points += n_sbl - N_SBL_R;
#endif
#ifdef USE_BUGEY_SP
  bugey_init();
  fprintf(stderr, "Bugey-spectr ");
#endif
#ifdef USE_CHOOZ
  chooz_init(old_new_main);
  fprintf(stderr, "Chooz ");
#endif  
#ifdef USE_PV
  PV_init();
  fprintf(stderr, "PaloV ");
#endif
#ifdef USE_KAML
  init_kaml();
  fprintf(stderr, "Kamland ");
#endif
#ifdef USE_DC
  dc_init(old_new_main);
  fprintf(stderr, "DC ");
#endif
#ifdef USE_DB
  DB_init();
  fprintf(stderr, "DayaB ");
#endif
#ifdef USE_RENO
  RENO_init();
  fprintf(stderr, "RENO ");
#endif
  fprintf(stderr, "\ntotal number of data points: %d\n", n_data_points);
  
  fit.invert_S();

  /**** end initialization *****/
  
  Param_5nu p;
  p.theta[I12] = asin(sqrt(0.31));
  p.theta[I13] = 0.5 * asin(sqrt(0.09));
  p.theta[I14] = p.theta[I15] = 0.;
  p.set_ang();
  p.dmq[0] = 0.;          // should always be zero!
  p.dmq[1] = 7.59e-5;
  p.dmq[2] = 2.32e-3;
  p.dmq[3] = 1.;
  p.dmq[4] = 2.;

  fit.pull_status[FLUX_NORM] = FREE;

  for(double sq2t = 0.01; sq2t < 1.; sq2t *= 1.05)
  {
    p.theta[I14] = 0.5 * asin(sqrt(sq2t));
    p.set_ang();

    for(p.dmq[3] = 0.1; p.dmq[3] < 12.; p.dmq[3] *= 1.05)
      printf("%f %f %f\n", sq2t, p.dmq[3], fit.chisq(p));
    fflush(stdout);
  }
  exit(0);

  /*
  // solar parameters
  for(int j = 0; j < 51; j++){
    double sq12 = 0.2 + 0.3 * j / 50.;
    p.theta[I12] = asin(sqrt(sq12));
    p.set_ang();

    for(int i = 0; i < 51; i++){
      p.dmq[1] = 6.8e-5 + (8.4e-5 - 6.8e-5) * i / 50.;
      printf("%e %e %e\n", sq12, p.dmq[1], fit.chisq(p));
    }
  }
  */
  return 0;  /*** end main ***/
} 

/*******************************
 *    function definitions     *
 *******************************/

void plot_fluxes(void)
{
  for(double e = 2.; e < 12.; e += 0.01){
    printf("%e  ", e);
    for(int i = 0; i < NISO; i++)
      printf("%e  ", global_flux(i, e, NEW));
    //printf("%e ", (flux[i][NEW].f(e) - flux[i][OLD].f(e))/flux[i][OLD].f(e) );
    printf("\n");
  }
  return;
} 

}
