#include "definitions.h"

using namespace ns_reactor;

namespace ns_reactor{

// global classes
Flux flux[NISO][2];  // old and new fluxes for each isotope using polynomials
Fit  fit;
Rate_coef rate;

int old_new_main = NEW;

void plot_fluxes(void);
void calc_cs_per_fission(void);  // in class_flux.cc


class Read_2dim{
 public:
  // at initialization provide number of x and y values and file name
  // set x_first to true if x changes and y is constant (not checked!!)
  Read_2dim(const int nx, const int ny, const char *file_name, bool x_first = false);
  ~Read_2dim(void);
  double interp(const double x, const double y);

 private:
  double **data;
  double xmin, xmax, ymin, ymax, dx, dy;
};


}

void reactor_anomaly(void);
void fit3p2(void);
void project_app(void);

/********************
 *    main          *
 ********************/

int reactors_main(void)
{
  fprintf(stderr, "For each experiment to be included define USE_ex with ex = SBL, CHOOZ, PV, KAML, DC, DB, RENO, BUGEY_SP, GAL\n");  

  /**** initialization *****/

  init_fluxes();
  rate.init();
  
  //calc_cs_per_fission(); return 0;
  //plot_fluxes(); return 0;

  int n_data_points = NBIN_CHISQ;
  
  fprintf(stderr, "using: ");
#ifdef USE_SBL
  fprintf(stderr, "SBL ");
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
  fprintf(stderr, "Kamland ");
  init_kaml();
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
#ifdef USE_GAL
  fprintf(stderr, "Gallium ");
  gallium_init();
#endif
  fprintf(stderr, "\ntotal number of data points: %d\n", n_data_points);
  
  fit.invert_S();

  // the chisq tables for the interpolation
  //Read_2dim chisq_C12(101, 143, "Dat/nue-carbon-chi2-spectrum_KL.dat");
  //Read_2dim chisq_sol(21, 31, "solar/Sterile-SU-27c-20/th13-th14_2.out", true);

  /**** end initialization *****/

  //print_sbl_data();
  //reactor_anomaly();
  //fit3p2();
  //project_app();

  
  Param_5nu p;
  p.theta[I12] = asin(sqrt(0.31));
  p.theta[I13] = 0.5 * asin(sqrt(0.09));
  p.theta[I14] = p.theta[I15] = 0.;
  p.set_ang();
  p.dmq[0] = 0.;          // should always be zero!
  p.dmq[1] = 7.59e-5;
  p.dmq[2] = 2.32e-3;
  p.dmq[3] = 0.9;
  p.dmq[4] = 2.;

  double min = 1.e6, s_m = 0., dm_m = 0.;

  for(int i = 0; i < 51; i++)
  {
    const double sq2t = exp10(-2.01+2.01*i/50.);
    //const double sq2t = 0.3*i/50.;
    p.theta[I14] = 0.5 * asin(sqrt(sq2t));
    p.set_ang();

    for(p.dmq[3] = 0.1; p.dmq[3] < 31.; p.dmq[3] *= 1.04)
    {
      /*
    for(int j = 0; j < 41; j++){

      const double sq2th13 = 0.3*j/40.;
      p.theta[I13] = 0.5 * asin(sqrt(sq2th13));
      p.set_ang();
      */

      /*
      double min_th13 = 1.e5;
      for(int j = 0; j < 21; j++){
	const double sq2th13 = 0.04 + (0.16-0.04)*j/20.;
	p.theta[I13] = 0.5 * asin(sqrt(sq2th13)); 
	p.set_ang();

	const double cc = fit.chisq(p) + chisq_sol.interp(sq2th13, sq2t);
	if(cc < min_th13) min_th13 = cc;
      }
      //const double c = min_th13 + chisq_C12.interp(sq2t, log10(p.dmq[3]));
      const double c = min_th13;
      */

      //const double c = fit.chisq(p) + chisq_sol.interp(sq2th13,sq2t); 
      const double c = fit.chisq(p);

      printf("%f %f %f\n", sq2t, p.dmq[3], c);
      //printf("%f %f %f\n", sq2th13, sq2t, c);

      if(c < min){
	min = c;
	s_m = sq2t;
	dm_m = p.dmq[3];
      }
    }
    fflush(stdout);
  }
  fprintf(stderr, "chisq_min = %f, sq2t = %f, dmq = %f  ", min, s_m, dm_m);

  // chisq for no oscillations
  p.dmq[3] = dm_m;
  p.theta[I14] = 0.;
  p.set_ang();
  const double c = fit.chisq(p);
  fprintf(stderr, "chisq_no-osc = %f, Dchisq_no-osc = %f\n", c, c-min);

  // plot the SBL probability for the bfp
  p.theta[I14] = 0.5 * asin(sqrt(s_m));
  p.set_ang();
  FILE *fp = fopen("tt_probL.dat", "w");
  print_sbl_pred(p, fp);
  fclose(fp);

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

void reactor_anomaly(void)
{
  Param_5nu p;
  p.theta[I12] = p.theta[I13] = 0.;
  p.theta[I14] = p.theta[I15] = 0.;
  p.set_ang();
  p.dmq[0] = 0.;          // should always be zero!
  p.dmq[1] = 7.59e-5;
  p.dmq[2] = 2.32e-3;
  p.dmq[3] = 10.;
  p.dmq[4] = 2.;

  fit.pull_status[FLUX_NORM] = FIXED;
  //fit.pull_status[FLUX_NORM] = ACTIVE;
  //printf("%f\n", fit.chisq(p));
  //exit(0);

  double xi[NPULLS];
  for(int i = 0; i < NPULLS; i++) xi[i] = 0.;

  double min = 1.e6, f = 0.; 
  const double c1 = fit.chisq(p);

  for(int i = 0; i < 51; i++){
    const double sq2t = 0.5*i/50.;
    p.theta[I14] = 0.5 * asin(sqrt(sq2t));
    p.set_ang();
    
    const double c = fit.chisq(p);
    printf("%e %f\n", 1.-0.5*sq2t, c);

    if(c < min) {
      min = c;
      f = 1.-0.5*sq2t;
    }
  }

  double fd = 1.;
  if(c1 - min > 1.){

    double del = 1. - f;

    while(del > 1.e-4){
      del *= 0.5;

      fd -= del;
      double sq2t = 2.*(1.-fd);
      p.theta[I14] = 0.5 * asin(sqrt(sq2t));
      p.set_ang();
      if(fit.chisq(p) - min < 1.) fd += del;
    }
  }

  fprintf(stderr, "chisq_min = %f, f = %f +/- %f, Dchisq = %f\n",
	  min, f, fd - f, c1 - min);

  exit(0);
  return;
}


/*******************************
 *    3+2                      *
 *******************************/

void fit3p2(void)
{
  Param_5nu p;
  p.theta[I12] = asin(sqrt(0.31));
  p.theta[I13] = 0.5 * asin(sqrt(0.09));
  p.theta[I14] = p.theta[I15] = 0.;
  p.set_ang();
  p.dmq[0] = 0.;          // should always be zero!
  p.dmq[1] = 7.59e-5;
  p.dmq[2] = 2.32e-3;
  p.dmq[3] = 1.694;
  p.dmq[4] = 2.;

#define NN 41

  FILE *fp = fopen("tt_dmq", "w");

  double min_abs = 1.e6;

  for(int i = 0; i < NN; i++)
  {
    p.dmq[3] = exp10(-1. + 2. * i/(NN-1.));

    double min_dmq = 1.e6;

    for(int j = 0; j < NN; j++)
    {
      p.dmq[4] = exp10(-1. + 2. * j/(NN-1.));

      double min = 1.e6;
      for(int k = 0; k < NN; k++){
	p.theta[I14] = M_PI/2. * k/(NN-1.);
        for(int l = 0; l < NN; l++){
	  p.theta[I15] = M_PI/2. * l/(NN-1.);
	  p.set_ang();

	  const double c = fit.chisq(p);
	  if(c < min) min = c;
	}
      }
      printf("%e %e %f\n", p.dmq[3], p.dmq[4], min);
      if(min < min_dmq) min_dmq = min;
    }		       
    fprintf(fp, "%e %f\n", p.dmq[3], min_dmq);
    fflush(fp);
    fflush(stdout);

    if(min_dmq < min_abs){ 
      min_abs = min_dmq;
      p.save("bfp_3p2");
    }
  }

  fclose(fp);
  exit(0);
  return;
}

#undef NN

/*******************************
 *    project                  *
 *******************************/

void project_app(void)
{
  Read_2dim chisq_LMB(101, 51, "Dat/3p1-LSND+MBa.out");

  Param_5nu p;
  p.theta[I12] = asin(sqrt(0.31));
  p.theta[I13] = 0.5 * asin(sqrt(0.09));
  p.theta[I14] = p.theta[I15] = 0.;
  p.set_ang();
  p.dmq[0] = 0.;          // should always be zero!
  p.dmq[1] = 7.59e-5;
  p.dmq[2] = 2.32e-3;
  p.dmq[3] = 1.694;
  p.dmq[4] = 2.;

#define NN 51

  double chisq_app[NN][NN], chisq_mu4[NN][NN];

  for(int i = 0; i < NN; i++)
    for(int j = 0; j < NN; j++)
      chisq_app[i][j] = chisq_mu4[i][j] = 1.e6;

  double min_app = 1.e6, min_dis = 1.e6, min_gl = 1.e6;

  for(int i = 0; i < NN; i++){
    const double sq_app = exp10(-3. + 3.*i/(NN-1.));

    for(int j = 0; j < NN; j++){
      const double Um4q = exp10(-2.5 + 2.5*j/(NN-1.));

      p.Ue[3] = sqrt(sq_app / (4. * Um4q));

      if(norm(p.Ue[3]) < 0.5 && norm(p.Ue[3]) + Um4q < 1.){

        for(int k = 0; k < NN; k++){
	  p.dmq[3] = exp10(-1. + 2.4 * k/(NN-1.));

	  const double c_dis = fit.chisq(p); 
	  const double c_app = chisq_LMB.interp(log10(p.dmq[3]), log10(sq_app));
	  const double c = c_dis + c_app;
	  //const double c = c_dis;

	  if(c < chisq_app[i][k]) chisq_app[i][k] = c;
	  if(c < chisq_mu4[j][k]) chisq_mu4[j][k] = c;

	  if(c_dis < min_dis) min_dis = c_dis;
	  if(c_app < min_app) min_app = c_app;
	  if(c < min_gl) min_gl = c;
	}
      }
    }
  }

  fprintf(stderr, "min_dis = %f, min_app = %f, min = %f, chisq_PG = %f\n", 
	  min_dis, min_app, min_gl, min_gl - min_dis - min_app);

  FILE *fp_app = fopen("tt_chisq_app", "w");
  FILE *fp_Um4 = fopen("tt_chisq_Um4", "w");

  for(int i = 0; i < NN; i++){
    const double sq_app = exp10(-3. + 3.*i/(NN-1.));
    const double Um4q   = exp10(-2.5 + 2.5*i/(NN-1.));

    for(int k = 0; k < NN; k++){
      p.dmq[3] = exp10(-1.+ 2.4 * k/(NN-1.));

      fprintf(fp_app, "%e %e %e\n", sq_app, p.dmq[3], chisq_app[i][k]);
      fprintf(fp_Um4, "%e %e %e\n", Um4q,   p.dmq[3], chisq_mu4[i][k]);
    }
  }
  exit(0);
  return;
}
