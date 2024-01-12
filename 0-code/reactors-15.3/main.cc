#include "definitions.h"
#include "class_minimize.h"

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


} // namespace

void reactor_anomaly(void);
void reactor_anomaly_th13(void);
void fit3p2(void);
void project_app(void);
void plot_sbl_prob(void);
void fit_3p1(void);
void fit_3p1_th13fix(void);
void calc_flux_table(void);

/*************************************
 * 3 flavour parameter and routines  *
 *************************************/

struct Param_3nu{
  double sq12;
  double sq13;
  double dmq21;
  double dmq31;  // sign of Dmq31 determines hierarchy
};

// flux norm free: makes only sense with using SBL data!
double chisq_3nu(Param_3nu p);

// treats flux norm as external param
double chisq_3nu(Param_3nu p, double f);

void convert_params(Param_3nu *p3nu, Param_5nu *p);

struct Parab{
  double x_min;
  double min;
  double sigm;
};

void fit_parabola(Param_3nu p, Parab *res, bool flux_free = false);
void fit_parabola(double xy[3][2], Parab *res);


/********************
 *    main          *
 ********************/

//int main(void)
//{
//  /**** initialization *****/
//#ifdef NO_GAL
//  fprintf(stderr, "ATTENTION: do not use GAL! If you want to use GAL remove #define NO_GAL in definitions.h\n");
//#endif
//
//  init_fluxes();
//  rate.init();
//  
//  //calc_cs_per_fission(); return 0;
//  //plot_fluxes(); return 0;
//
//  int n_data_points = NBIN_CHISQ;
//  
//  fprintf(stderr, "initializing: ");
//  fprintf(stderr, "DayaB-flux ");
//  init_db_flux();
//  fprintf(stderr, "SBL ");
//  const int n_sbl = init_sbl_reactors(old_new_main);
//  n_data_points += n_sbl - N_SBL_R;
//  fprintf(stderr, "Bugey-spectr ");
//  bugey_init();
//  fprintf(stderr, "Chooz ");
//  chooz_init(old_new_main);
//  fprintf(stderr, "PaloV ");
//  PV_init();
//  fprintf(stderr, "DANSS ");
//  danss_init(old_new_main);
//  fprintf(stderr, "DC ");
//  dc_init(old_new_main);
//  fprintf(stderr, "DayaB ");
//  DB_init();
//  fprintf(stderr, "Kamland ");
//  init_kaml();
//  fprintf(stderr, "RENO ");
//  RENO_init();
//  fprintf(stderr, "Gallium ");
//  gallium_init();
//  
//  fprintf(stderr, "\n");  
//  //fprintf(stderr, "\ntotal number of data points: %d\n", n_data_points);
//  
//  
//  //fit_3p1_th13fix(); exit(0); // to be called before fit.invert_S();  
//  
//  // fit.pull_status[PULL_U235_1] = FIXED;
//  // fit.pull_status[PULL_U235_2] = FIXED;
//  // fit.pull_status[PULL_P239_1] = FIXED;
//  // fit.pull_status[PULL_P239_2] = FIXED;
//  // fit.pull_status[PULL_P241_1] = FIXED;
//  // fit.pull_status[PULL_P241_2] = FIXED;
//
//  // fit.pull_status[FLUX_COR] = FIXED;
//
//  //fit.S_pull[PULL_P241_0][PULL_P241_0] = norm(0.10);  
//  //fit.S_pull[PULL_U238]  [PULL_U238]   = norm(0.10);
//
//  //fit.S_pull[FLUX_NORM][FLUX_NORM] = norm(0.10) / 2.;
//
//  fit.invert_S();
//
//  /*********** which experiments to use in chisq ***********************/
//  bool use_ex[N_EXP];
//  
//  // to include one of the following exps set the corresponding entry to true
//  // SBL, CHOOZ, PV, KAML, DC, DB, RENO, BUG_SP, GAL, DB_FLUX
//
//  for(int i = 0; i < N_EXP; i++) use_ex[i] = false;
//  
//  //use_ex[SBL] = true;
//  // use_ex[CHOOZ] = false;
//  // use_ex[PV] = false;
//  // use_ex[DB] = false;
//  // use_ex[DC] = false;
//  // use_ex[KAML] = false;
//  // use_ex[GAL] = false;
//  // use_ex[RENO] = false;
//  // use_ex[BUG_SP] = false;
//  use_ex[DB_FLUX] = true;
//
//  fit.set_experiments(use_ex);
//
//#ifdef FIT_FLUXES  
//  // fit U235 and P239 fluxes
//
//  Param_5nu p;
//  p.theta[I12] = asin(sqrt(0.306));
//  p.theta[I13] = asin(sqrt(0.0217));
//  p.theta[I14] = p.theta[I15] = 0.;
//  p.set_ang();
//  p.dmq[0] = 0.;          // should always be zero!
//  p.dmq[1] = 7.5e-5;
//  p.dmq[2] = 2.51e-3;     
//  p.dmq[3] = 1.;
//  p.dmq[4] = 2.;
//  
//  fit.pull_status[PULL_U235_0] = fit.pull_status[PULL_P239_0] = FIXED;
//  double xi[NPULLS];  
//  for(int k = 0; k < NPULLS; k++) xi[k] = 0.;
//
//  // Huber flux pred
//  const double sig235 = 6.69; 
//  const double sig239 = 4.36;
//
//  double min = 1000.;
//  
//  for(int i = 0; i < 41; i++){
//    const double f235 = 0.8 + 0.3 * i / 40.;
//    xi[PULL_U235_0] = f235 - 1.;
//
//    for(int j = 0; j < 41; j++){
//      const double f239 = 0.55 + 0.65 * j / 40.;
//      xi[PULL_P239_0] = f239 - 1.;
//      //double c = fit.chisq(p, xi, true) + norm( (1.-f235)/0.05) + norm( (1.-f239)/0.05);
//      double c = fit.chisq(p, xi, true);
//      if(c < min) min = c;
//      fprintf(stdout, "%f %f %f\n", sig235*f235, sig239*f239, c);
//    }
//  }
//  fprintf(stderr, "chi2_min from scan: %f\n", min);
//
//  // check result by flux minimization
//  fit.pull_status[PULL_U235_0] = fit.pull_status[PULL_P239_0] = FREE;
//  min = fit.chisq(p, xi);  
//  fprintf(stderr, "chi2_min from pull minimiziation: %f, sig_235 = %f, sig_239 = %f\n",
//	  min, (xi[PULL_U235_0] + 1.) * sig235,  (xi[PULL_P239_0] + 1.) * sig239);
//
//  // global rescaling of all fluxes
//  fit.pull_status[PULL_U235_0] = fit.pull_status[PULL_P239_0] = FIXED;
//  fit.pull_status[PULL_U238]   = fit.pull_status[PULL_P241_0] = FIXED;
//  for(int k = 0; k < NPULLS; k++) xi[k] = 0.;
//  fit.pull_status[FLUX_NORM] = ACTIVE;
//  double c = fit.chisq(p, xi);
//  fprintf(stderr, "chi2_min global normalization: %f (Delta chi^2 = %f)\n", c, c-min);
//
//  // chisq of H-M prediction
//  for(int i = 0; i < FLUX_NORM; i++)
//    fit.pull_status[i] = ACTIVE;
//  fit.pull_status[FLUX_NORM] = FIXED;
//  c = fit.chisq(p);
//  fprintf(stderr, "chi2_HuberMueller: %f (Delta chi^2 = %f)\n", c, c-min);
//  
//  fit.pull_status[FLUX_NORM] = FREE;
//  c = fit.chisq(p);
//  fprintf(stderr, "chi2_HuberMueller + free norm: %f (Delta chi^2 = %f)\n", c, c-min);
//
//  exit(0);
//#endif // FIT_FLUXES
//
//  
//  
//  //Reactor_anomaly(); exit(0);
//  //reactor_anomaly_th13(); exit(0);
//
//  Param_3nu pp;
//  pp.sq12 = 0.308;
//  pp.sq13 = 0.5*(1. - sqrt(1. - 0.0841));
//  pp.dmq21 = 7.59e-5;
//  pp.dmq31 = 2.5e-3 + pp.sq12 * pp.dmq21;
//  
//  Param_5nu p5;
//  convert_params(&pp, &p5);
//  //DB_spectrum(p5);
//
// 
////#define CHECK_ORDER  
//#ifdef CHECK_ORDER
//  {
//  const int NN = 51;
//  double min = 1.e4;
//  double chisq_dmq[NN][2], cc[NN][NN];
//  for(int i = 0; i < NN; i++)
//    chisq_dmq[i][1] = 1.e5;
//  
//  for(int j = 0; j < NN; j++){
//    pp.dmq31 = -(2.1e-3 + (3.e-3 - 2.1e-3) * j/(NN-1.));
//    
//    for(int i = 0; i < NN; i++){
//      double sq2t = 0.065 + (0.1-0.065) * i / (NN-1.);
//      pp.sq13 = 0.5 * (1. - sqrt(1. - sq2t));
//
//      cc[i][j] = chisq_3nu(pp, 1.);
//      if(cc[i][j] < chisq_dmq[j][1]){
//	chisq_dmq[j][1] = cc[i][j];
//      }
//    }
//    double dmqee = pp.dmq31 - pp.sq12 * pp.dmq21;
//    if(chisq_dmq[j][1] < min) min = chisq_dmq[j][1];
//    printf("%f %f %f\n", fabs(pp.dmq31) * 1.e3, fabs(dmqee) * 1.e3, chisq_dmq[j][1]);
//  }
//  fprintf(stderr, "min = %f\n", min);
//  exit(0);
//  }
//#endif // CHECK_ORDER
//
//  const int NN = 51;
//  int sq_min = 0, dm_min = 0; 
//  double chisq_sq2[NN][2], chisq_dmqee[NN][2], cc[NN][NN];
//  for(int i = 0; i < NN; i++)
//    chisq_sq2[i][1] = chisq_dmqee[i][1] = 1.e5;
//  
//  for(int i = 0; i < NN; i++){
//    double sq2t = 0.06 + (0.2-0.06) * i / (NN-1.);
//    pp.sq13 = 0.5 * (1. - sqrt(1. - sq2t));
//    chisq_sq2[i][0] = sq2t;
//      
//    for(int j = 0; j < NN; j++){
//      double dmqee = 2.1e-3 + (3.e-3 - 2.1e-3) * j/(NN-1.);
//      pp.dmq31 = dmqee + pp.sq12 * pp.dmq21;
//
//      cc[i][j] = chisq_3nu(pp, 1.) + norm((dmqee - 2.44e-3)/(0.09e-3));  // add prior on Dmq!!
//
//      printf("%f %f %f\n", sq2t, dmqee * 1.e3, cc[i][j]);
//      
//      if(cc[i][j] < chisq_sq2[i][1]){
//	chisq_sq2[i][1] = cc[i][j];
//      }
//    }
//  }
//  
//
//  double c_min = 1.e5;
//  for(int j = 0; j < NN; j++){
//    chisq_dmqee[j][0] = 2.1e-3 + (3.e-3 - 2.1e-3) * j/(NN-1.);
//    
//    for(int i = 0; i < NN; i++){
//      if(cc[i][j] < chisq_dmqee[j][1]) chisq_dmqee[j][1] = cc[i][j];
//      if(cc[i][j] < c_min){
//	c_min = cc[i][j];
//	dm_min = j;
//	sq_min = i;
//      }
//    }
//  }
//  fprintf(stderr, "min = %f\n", c_min);
//
//  FILE *fp1 = fopen("tt-out-sq2","w");
//  FILE *fp2 = fopen("tt-out-dmq","w");
//  for(int i = 0; i < NN; i++){
//    fprintf(fp1, "%f %f\n", chisq_sq2[i][0], chisq_sq2[i][1]-c_min);
//    fprintf(fp2, "%f %f\n", chisq_dmqee[i][0], chisq_dmqee[i][1]-c_min);
//  }
//
//  
//  Parab ch;
//  double xy[3][2];
//  for(int j = 0; j < 3; j++){
//    int i = sq_min + 5 * (j-1); 
//    xy[j][0] = chisq_sq2[i][0];
//    xy[j][1] = chisq_sq2[i][1];
//  }
//  fit_parabola(xy, &ch);
//  fprintf(stderr, "sq2th13 = %f +/- %f\n", ch.x_min, ch.sigm);
//
//  for(int j = 0; j < 3; j++){
//    int i = dm_min + 5 * (j-1); 
//    xy[j][0] = chisq_dmqee[i][0];
//    xy[j][1] = chisq_dmqee[i][1];
//  }
//  fit_parabola(xy, &ch);
//  fprintf(stderr, "dmqee = %f +/- %f * 1e-3\n", ch.x_min*1.e3, ch.sigm*1.e3);
//  exit(0);
//
//  
//  fit_parabola(pp, &ch, true);
//  double min_free = ch.min;
//  fit_parabola(pp, &ch, false);
//  double min_fix = ch.min;
//  fprintf(stderr, "flux free = %f, fix = %f, delta = %f\n", min_free, min_fix, min_fix-min_free); 
//  exit(0);  
//
//  for(pp.dmq31 = 2.3e-3; pp.dmq31 < 2.501e-3; pp.dmq31 += 0.1e-3){  
//    fit_parabola(pp, &ch);
//    printf("%e %e %e %e\n", pp.dmq31, ch.x_min, ch.x_min - ch.sigm, ch.x_min + ch.sigm);
//  }
//  exit(0);
//
//  fit_parabola(pp, &ch);
//  fprintf(stderr, "sq2t = %f +/- %f, chisq_min = %f\n", ch.x_min, ch.sigm, ch.min);
//  
//  for(int i = 0; i < 51; i++){
//    double sq2t = 0.08 + 0.06 * i / 50.;
//    pp.sq13 = 0.5 * (1. - sqrt(1. - sq2t));
//
//    const double w = norm(sq2t - ch.x_min) / norm(ch.sigm) + ch.min;
//    const double v = norm(sq2t - ch.x_min) / norm(0.018) + ch.min;
//    
//    printf("%f %f %f %f\n", sq2t, chisq_3nu(pp, 1.), w, v);
//  }
//  exit(0);
//  
//  // Param_5nu p5;
//  // convert_params(&pp, &p5);
//  // plot_dc_pred(p5);
//  // exit(0);
//
//  // 3nu theta13 fitting
//  /*
//  Param_3nu p;
//  p.sq12 = 0.3;
//  p.sq13 = 0.;
//  p.dmq21 = 7.59e-5;
//  p.dmq31 = 2.32e-3;
//
//  const int n_dm = 60;
//  const int n_sq = 60;
//
//  double chisq_dm[n_dm], chisq_sq[n_sq];
//  for(int i = 0; i < n_dm; i++) chisq_dm[i] = 1.e6;
//  for(int i = 0; i < n_sq; i++) chisq_sq[i] = 1.e6;
//
//  double dm=0., sm=0., m=1000.;
//
//  for(int i = 0; i < n_dm; i++){
//    double dmqee = 2.e-3 + 0.9e-3 * double(i)/(n_dm-1.);
//    p.dmq31 = dmqee + p.sq12 * p.dmq21; 
//
//    for(int j = 0; j < n_sq; j++){
//
//      double sq2 = 0.06 + 0.05 * double(j)/(n_sq-1.);
//      p.sq13 = norm(sin(.5*asin(sqrt(sq2))));
//      const double c = chisq_3nu(p, 1.);
//      printf("%f %f %f\n", sq2, dmqee*1.e3, c);
//
//      if(c < m){ 
//        m = c; 
//	dm = dmqee;
//	sm = sq2;
//      }
//      if(c < chisq_dm[i]) chisq_dm[i] = c; 
//      if(c < chisq_sq[j]) chisq_sq[j] = c; 
//    }  
//  }
//  fprintf(stderr, "min = %f, sq2t = %f, dmq = %f\n", m, sm, dm);
//
//  FILE *fp = fopen("out.chisq.dmq", "w");
//  for(int i = 0; i < n_dm; i++){
//    double dmqee = 2. + .9 * double(i)/(n_dm-1.);
//    fprintf(fp, "%e %e\n", dmqee, chisq_dm[i]-m);
//  }
//  fclose(fp);
//  
//  fp = fopen("out.chisq.sq", "w");
//  for(int i = 0; i < n_sq; i++){
//    double sq2 = 0.06 + 0.05 * double(i)/(n_sq-1.);
//    fprintf(fp, "%e %e\n", sq2, chisq_sq[i]-m);
//  }
//  fclose(fp);
//  */
//  return 0;  /*** end main ***/
//} 

/*******************************
 *    function definitions     *
 *******************************/

// converts 3nu parameter to 5nu parameter
void convert_params(Param_3nu *p3nu, Param_5nu *p)
{
  p->theta[I12] = asin(sqrt(p3nu->sq12));
  p->theta[I13] = asin(sqrt(p3nu->sq13));
  p->theta[I14] = 0.; 
  p->theta[I15] = 0.; 
  p->set_ang();

  p->dmq[0] = 0.;          // should always be zero!
  p->dmq[1] = p3nu->dmq21;
  p->dmq[2] = p3nu->dmq31; 
  p->dmq[3] = 1.;          // some number (irrelevant)
  p->dmq[4] = 2.;          // some number (irrelevant)
  return;
}

// flux norm free: makes only sense with using SBL data!
double chisq_3nu(Param_3nu p3nu)
{
  fit.pull_status[FLUX_NORM] = FREE;
  Param_5nu p;
  convert_params(&p3nu, &p);
  return fit.chisq(p);
}

// treats flux norm as external param
double chisq_3nu(Param_3nu p3nu, double f)
{
  Param_5nu p;
  convert_params(&p3nu, &p);

  fit.pull_status[FLUX_NORM] = FIXED;
  double xi[NPULLS];
  for(int p = 0; p < NPULLS; p++) xi[p] = 0.;
  xi[FLUX_NORM] = f - 1.;

  return fit.chisq(p, xi, true);
}

void fit_parabola(double xy[3][2], Parab *res)
{
  double A[3][4], result[3];

  for(int i = 0; i < 3; i++){
    A[i][0] = norm(xy[i][0]);
    A[i][1] = xy[i][0];
    A[i][2] = 1.;
    A[i][3] = -xy[i][1];
  }

  singsolve(3, &A, &result);

  if(result[0] < 0.){
    fprintf(stderr, "fit_parabola: x^2 coefficient negative\n"); exit(0);
  }
  
  res->sigm = 1./sqrt(result[0]);
  res->x_min = -0.5 * result[1] / result[0];
  res->min = result[2] - norm(result[1]) / (4. * result[0]);
  return;
}


void fit_parabola(Param_3nu p, Parab *res, bool flux_free)
{
  double A[3][4], result[3];
  const double sq2t[3] = {0.075, 0.09, 0.105};

  for(int i = 0; i < 3; i++){
    p.sq13 = 0.5*(1. - sqrt(1. - sq2t[i]));
    A[i][0] = norm(sq2t[i]);
    A[i][1] = sq2t[i];
    A[i][2] = 1.;
    A[i][3] = (flux_free ? -chisq_3nu(p) : -chisq_3nu(p, 1.));
  }

  singsolve(3, &A, &result);

  if(result[0] < 0.){
    fprintf(stderr, "fit_parabola: x^2 coefficient negative\n"); exit(0);
  }
  
  res->sigm = 1./sqrt(result[0]);
  res->x_min = -0.5 * result[1] / result[0];
  res->min = result[2] - norm(result[1]) / (4. * result[0]);
  return;
}


/********** caclulate flux th13 table  *******************/


const double sq13v[35] = {
 0.0000, 0.0050, 0.0100, 0.0120, 0.0140, 0.0160, 0.0170, 0.0180, 0.0190, 0.0195,
  0.0200, 0.0205, 0.0210, 0.0215, 0.0220, 0.0225, 0.0230, 0.0235, 0.0240, 0.0245,
  0.0250, 0.0255, 0.0260, 0.0265, 0.0270, 0.0280, 0.0290, 0.0300, 0.0320, 0.0350,
  0.0400, 0.0450, 0.0500, 0.0600, 0.0700
};

void calc_flux_table(void)
{ 
  Param_3nu p;
  p.sq12 = 0.304;
  p.sq13 = 0.;
  p.dmq21 = 7.50e-5;
  p.dmq31 = -2.448e-3 + p.dmq21;

  /*********** which experiments to use in chisq ***********************/
  bool use_all[N_EXP], use_1[N_EXP], use_2[N_EXP], use_3[N_EXP];

  // to include one of the following exps set the corresponding entry to true
  // SBL, CHOOZ, PV, KAML, DC, DB, RENO, BUG_SP, GAL
  
  for(int i = 0; i < N_EXP; i++){
    use_all[i] = true;
    use_1[i] = false;
    use_2[i] = false;
    use_3[i] = false;
  }

  use_1[SBL] = use_1[BUG_SP] = true;
  use_2[CHOOZ] = use_2[PV] = use_2[DC] = use_2[SBL] = use_2[BUG_SP] = true;
  use_3[DB] = use_3[RENO] = true;
  use_all[KAML] = use_all[GAL] = false;

  for(int i = 0; i < 35; i++){
    p.sq13 = sq13v[i];
    for(int j = 0; j < 41; j++){
      const double f = 0.85 + 0.2 * j/40.;

      fit.set_experiments(use_1);
      printf("%e %e %e ", sq13v[i], f, chisq_3nu(p, f)); 

      fit.set_experiments(use_2);
      printf(" %f ", chisq_3nu(p, f)); 

      fit.set_experiments(use_3);
      printf(" %f ", chisq_3nu(p, f)); 

      fit.set_experiments(use_all);
      printf(" %f\n", chisq_3nu(p, f)); 
    }
  }
  exit(0);
  return;
}

/********** end 3 flavour *******************/


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

void plot_sbl_prob(void)
{
  Param_5nu p;
  p.theta[I12] = 0.; 
  p.theta[I13] = 0.5*asin(sqrt(0.085));
  //p.theta[I13] = 0.;
  //p.theta[I14] = 0.5*asin(sqrt(0.057));
  p.theta[I14] = 0.;
  p.theta[I15] = 0.;
  p.set_ang();
  p.dmq[0] = 0.;          // should always be zero!
  p.dmq[1] = 7.59e-5;
  p.dmq[2] = 2.32e-3;
  p.dmq[3] = .9;
  p.dmq[4] = 20.;

  print_sbl_pred(p, stdout);
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

void reactor_anomaly_th13(void)
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

  for(double sq2t13 = 0.05; sq2t13 < 0.12; sq2t13 += 0.001){

    p.theta[I13] = 0.5 * asin(sqrt(sq2t13));
    p.theta[I14] = 0.;
    p.set_ang();

    double min = 1.e6, f = 0.; 
    const double c1 = fit.chisq(p);

    for(int i = 0; i < 51; i++){
      const double sq2t = 0.5*i/50.;
      p.theta[I14] = 0.5 * asin(sqrt(sq2t));
      p.set_ang();
    
      const double c = fit.chisq(p);
      if(c < min) {
	min = c;
	f = 1.-0.5*sq2t;
      }
    }
    printf("%f %f %f\n", sq2t13, c1 - min, f);
  }

  exit(0);
  return;
}


/*******************************
 *    3+1                      *
 *******************************/


// plot th14 - dmq41 plane with th13 fixed
//void fit_3p1_th13fix(void)
//{
//  // to be called here!
//  fit.invert_S();
//
//  /*********** which experiments to use in chisq ***********************/
//  bool use_ex[N_EXP];
//  
//  // to include one of the following exps set the corresponding entry to true
//  // SBL, CHOOZ, PV, KAML, DC, DB, RENO, BUG_SP, GAL, DB_FLUX
//
//  for(int i = 0; i < N_EXP; i++) use_ex[i] = false;
//  
//  // use_ex[CHOOZ] = false; 
//  // use_ex[DC] = true;
//  // use_ex[DB] = true;
//
//  
//  //use_ex[SBL] = true;
//  //use_ex[PV] = true;
//  //use_ex[RENO] = true;
//  //use_ex[BUG_SP] = true;
//  //use_ex[DANSS] = true;
//  use_ex[DB_FLUX] = true;
//  //use_ex[KAML] = true;
//  //use_ex[GAL] = true;
//  
//  fit.set_experiments(use_ex);
//
//  Param_5nu p;
//  p.theta[I12] = asin(sqrt(0.306));
//  p.theta[I13] = asin(sqrt(0.0217));
//  p.theta[I14] = p.theta[I15] = 0.;
//  p.set_ang();
//  p.dmq[0] = 0.;          // should always be zero!
//  p.dmq[1] = 7.5e-5;
//  p.dmq[2] = 2.51e-3;     
//  p.dmq[3] = 1.;
//  p.dmq[4] = 2.;
//
//  const double chisq_noosc = fit.chisq(p);
//
//  for(int i = 0; i < FLUX_NORM; i++)
//    fit.pull_status[i] = FIXED;
//  fit.pull_status[FLUX_NORM] = FREE;
//  fprintf(stdout, "chisq = %f chiq_norm_free = %f\n", chisq_noosc, fit.chisq(p));
//  exit(0);
//  
//  double min = 1.e6;
//  
//  for(int i = 0; i < 51; i++){
//    const double sq = exp10(-3.+3.*i/50.);
//    p.theta[I14] = asin(sqrt(sq));
//    p.set_ang();
//
//    for(int j = 0; j < 201; j++)
//    {
//      //p.dmq[3] = exp10(-1.4 + 2.4*j/100.);
//      p.dmq[3] = exp10(-2. + 3.*j/200.);
//
//      double c = fit.chisq(p);
//      fprintf(stdout, "%e %e %f\n", sq, p.dmq[3], c);
//      if(c < min) min = c;
//    }
//  }
//  fprintf(stderr, "chisq_min = %f, Delta chisq_noosc = %f\n", min, chisq_noosc - min);
//  exit(0);
//  return;
//}


double min_funct_3p1(double x[3])
{
  Param_5nu p;
  p.theta[I12] = asin(sqrt(0.31));
  p.theta[I13] = x[0];
  p.theta[I14] = x[1];
  p.theta[I15] = 0.;
  p.set_ang();
  p.dmq[0] = 0.;          // should always be zero!
  p.dmq[1] = 7.59e-5;
  p.dmq[2] = 2.32e-3;
  p.dmq[3] = x[2];
  p.dmq[4] = 10.;
  
  return fit.chisq(p);
}

// plot th14 - dmq41 plane using minimization in th13
//void fit_3p1(void)
//{
//  /*********** which experiments to use in chisq ***********************/
//  bool use_ex[N_EXP];
//
//  // to include one of the following exps set the corresponding entry to true
//  // SBL, CHOOZ, PV, KAML, DC, DB, RENO, BUG_SP, GAL
//  
//  for(int i = 0; i < N_EXP; i++)
//    use_ex[i] = true;
//  use_ex[KAML] = use_ex[GAL] = false;
//  fit.set_experiments(use_ex);
//
//  Minimize chisq_min(min_funct_3p1, 3);
//  bool fix[3] = {false, true, true}; // fix dmq's
//
//  const int nn = 81;
//  const double th13 = 0.5*asin(sqrt(0.085));
//
//  for(int i = 0; i < nn; i++)
//  {
//    double Ue4q = exp10(-3.5 + 2.8 * i/(nn-1.));
//
//    for(int j = 0; j < nn; j++)
//    {
//      double dmq41 = exp10(-3. + 3. * j/(nn-1.));
//
//      double x[3] = {th13, asin(sqrt(Ue4q)), dmq41};      
//
//      use_ex[GAL] = false;
//      fit.set_experiments(use_ex);
//      const double c1 = chisq_min.find_min(x, fix);
//
//      use_ex[GAL] = true;
//      fit.set_experiments(use_ex);
//      const double c2 = chisq_min.find_min(x, fix);
//
//      printf("%e %e %f  %f\n", Ue4q, dmq41, c1, c2);
//    }		       
//    fflush(stdout);
//  }
//  exit(0);
//  return;
//}
//
//
//void fit_3p1_th13_th14(void)
//{
//  /*********** which experiments to use in chisq ***********************/
//  bool use_all[N_EXP], use_sbl[N_EXP], use_lbl[N_EXP];
//
//  // to include one of the following exps set the corresponding entry to true
//  // SBL, CHOOZ, PV, KAML, DC, DB, RENO, BUG_SP, GAL
//  
//  for(int i = 0; i < N_EXP; i++){
//    use_all[i] = true;
//    use_sbl[i] = false;
//    use_lbl[i] = true;
//  }
//
//  use_all[KAML] = use_all[GAL] = false;
//  use_sbl[SBL] = use_sbl[BUG_SP] = true;
//  use_lbl[SBL] = use_lbl[KAML] = use_lbl[BUG_SP] = use_lbl[GAL] = false;
//
//
//  Param_5nu p;
//  p.theta[I12] = asin(sqrt(0.31));
//  p.theta[I13] = 0.;
//  p.theta[I14] = p.theta[I15] = 0.;
//  p.set_ang();
//  p.dmq[0] = 0.;          // should always be zero!
//  p.dmq[1] = 7.59e-5;
//  p.dmq[2] = 2.4e-3;
//  p.dmq[3] = 10.;
//  p.dmq[4] = 2.;
//
//  fit.pull_status[FLUX_NORM] = FIXED;
//  //fit.pull_status[FLUX_NORM] = ACTIVE;
//
//
//  double min = 1.e6; 
//
//  for(int i = 0; i < 51; i++){
//    const double sq2t13 = 0.05 + 0.05*i/50.;
//    p.theta[I13] = 0.5 * asin(sqrt(sq2t13));
//
//    for(int j = 0; j < 51; j++){
//      const double Ue4q = exp10(-3. + 2.* j / 50.);
//      p.theta[I14] = asin(sqrt(Ue4q));
//      p.set_ang();
//
//      fit.set_experiments(use_all);
//    
//      const double c = fit.chisq(p);
//      printf("%e %e %f  ", Ue4q, sq2t13, c);
//
//      if(c < min) {
//        min = c;
//      }
//
//      fit.set_experiments(use_sbl);
//      printf(" %f  ", fit.chisq(p));
//      fit.set_experiments(use_lbl);
//      printf(" %f\n", fit.chisq(p));
//    }
//    fflush(stdout);
//  }
//  return;
//}

/*******************************
 *    3+2                      *
 *******************************/

const double th12 = asin(sqrt(0.31));
const double th13 = 0.5 * asin(sqrt(0.09));
const double thst = 0.5 * asin(sqrt(0.1));

double min_funct(double x[4])
{
  double w = fmin(x[0], x[1]);
  if(w < 0.) return 1000. + norm(w);
  w = fmax(x[0], x[1]);
  if(w > M_PI/2.) return 1000. + norm(w);
  
  Param_5nu p;
  p.theta[I12] = th12;
  p.theta[I13] = th13;
  p.theta[I14] = x[0];
  p.theta[I15] = x[1];
  p.set_ang();
  p.dmq[0] = 0.;          // should always be zero!
  p.dmq[1] = 7.59e-5;
  p.dmq[2] = 2.32e-3;
  p.dmq[3] = x[2];
  p.dmq[4] = x[3];
  
  return fit.chisq(p);
}

// using minimization in angles
//void fit3p2(void)
//{
//  Minimize chisq_min(min_funct, 4);
//  bool fix[4] = {false, false, true, true}; // fix dmq's
//
//#define NN 201
//
//  FILE *fp = fopen("out.dmq_sbl_min1", "w");
//
//  double min_abs = 1.e6;
//
//  for(int i = 0; i < NN; i++)
//  {
//    double dmq41 = exp10(-1. + 2. * i/(NN-1.));
//
//    double min_dmq = 1.e6;
//
//    for(int j = 0; j < NN; j++)
//    {
//      double dmq51 = exp10(-1. + 2. * j/(NN-1.));
//
//      double x[4] = {thst, thst, dmq41, dmq51};
//      
//      const double c = chisq_min.find_min(x, fix);
//      
//      if(c < min_abs){ 
//	min_abs = c;
//	FILE *fp1 = fopen("out.bfp_3p2_sbl_min1", "w");
//	fprintf(fp1, "%e %e %f %f\n", x[0],x[1],x[2],x[3]);
//	fprintf(fp1, "chisq = %f\n", min_abs);
//        fclose(fp1);
//      }
//      printf("%e %e %f\n", dmq41, dmq51, c);
//      if(c < min_dmq) min_dmq = c;
//    }		       
//    fprintf(fp, "%e %f\n", dmq41, min_dmq);
//    fflush(fp);
//    fflush(stdout);
//  }
//
//  fclose(fp);
//  exit(0);
//  return;
//}

#undef NN

/*******************************
 *    project                  *
 *******************************/

//void project_app(void)
//{
//  Read_2dim chisq_LMB(101, 51, "Dat/3p1-LSND+MBa.out");
//
//  Param_5nu p;
//  p.theta[I12] = asin(sqrt(0.31));
//  p.theta[I13] = 0.5 * asin(sqrt(0.09));
//  p.theta[I14] = p.theta[I15] = 0.;
//  p.set_ang();
//  p.dmq[0] = 0.;          // should always be zero!
//  p.dmq[1] = 7.59e-5;
//  p.dmq[2] = 2.32e-3;
//  p.dmq[3] = 1.694;
//  p.dmq[4] = 2.;
//
//#define NN 51
//
//  double chisq_app[NN][NN], chisq_mu4[NN][NN];
//
//  for(int i = 0; i < NN; i++)
//    for(int j = 0; j < NN; j++)
//      chisq_app[i][j] = chisq_mu4[i][j] = 1.e6;
//
//  double min_app = 1.e6, min_dis = 1.e6, min_gl = 1.e6;
//
//  for(int i = 0; i < NN; i++){
//    const double sq_app = exp10(-3. + 3.*i/(NN-1.));
//
//    for(int j = 0; j < NN; j++){
//      const double Um4q = exp10(-2.5 + 2.5*j/(NN-1.));
//
//      p.Ue[3] = sqrt(sq_app / (4. * Um4q));
//
//      if(norm(p.Ue[3]) < 0.5 && norm(p.Ue[3]) + Um4q < 1.){
//
//        for(int k = 0; k < NN; k++){
//	  p.dmq[3] = exp10(-1. + 2.4 * k/(NN-1.));
//
//	  const double c_dis = fit.chisq(p); 
//	  const double c_app = chisq_LMB.interp(log10(p.dmq[3]), log10(sq_app));
//	  const double c = c_dis + c_app;
//	  //const double c = c_dis;
//
//	  if(c < chisq_app[i][k]) chisq_app[i][k] = c;
//	  if(c < chisq_mu4[j][k]) chisq_mu4[j][k] = c;
//
//	  if(c_dis < min_dis) min_dis = c_dis;
//	  if(c_app < min_app) min_app = c_app;
//	  if(c < min_gl) min_gl = c;
//	}
//      }
//    }
//  }
//
//  fprintf(stderr, "min_dis = %f, min_app = %f, min = %f, chisq_PG = %f\n", 
//	  min_dis, min_app, min_gl, min_gl - min_dis - min_app);
//
//  FILE *fp_app = fopen("tt_chisq_app", "w");
//  FILE *fp_Um4 = fopen("tt_chisq_Um4", "w");
//
//  for(int i = 0; i < NN; i++){
//    const double sq_app = exp10(-3. + 3.*i/(NN-1.));
//    const double Um4q   = exp10(-2.5 + 2.5*i/(NN-1.));
//
//    for(int k = 0; k < NN; k++){
//      p.dmq[3] = exp10(-1.+ 2.4 * k/(NN-1.));
//
//      fprintf(fp_app, "%e %e %e\n", sq_app, p.dmq[3], chisq_app[i][k]);
//      fprintf(fp_Um4, "%e %e %e\n", Um4q,   p.dmq[3], chisq_mu4[i][k]);
//    }
//  }
//  exit(0);
//  return;
//}
