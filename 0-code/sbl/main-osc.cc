#include "def-reactors.h"  // includes also "definitions.h"

#define LENGTH 200
#define N_DMQ 101

#define DMQ_MIN 0.0316
#define DMQ_MAX 31.622776  // exp10(1.5)

#define DEL_LOG ( (log10(DMQ_MAX) - log10(DMQ_MIN))/(N_DMQ-1.) )


/* for reactor experiments */

Flux flux[NISO][2];  // old and new reactor fluxes for each isotope
Fit  fit;            // chisq for the reactors

#ifdef OLD_FLUX
int old_new_main = OLD;
#else
int old_new_main = NEW;
#endif


/***************************************************
 *      define the data sets
 ***************************************************/

struct Read_Data_set{
  int nr;        // number of the set
  double min;    // minimum of the set
  double min_gl; // chisq of the set at the point of the global minimum
};

enum Sets{
  set_mb475,
  set_mb300,
  set_lsnd,
  set_l_k_n,
  set_nev_app,   // KARMEN + NOMAD
  set_nev,       // without MB
  set_nev_mb475, 
  set_app_mb300,      
  set_app_mb475,      
  set_dis,
  set_all_mb300,
  set_all_mb475,
  set_nev_mb300,
  set_K_N_mb300,
  set_K_N_mb475,
  set_evid,      // LSND + MBA
  NUM_SETS
};

struct Data_set{
  bool incl[NUM_EXP];
  char name[LENGTH];

} data_set[NUM_SETS];

#ifdef PC300
# define DATA_SET set_all_mb300
#endif

#ifndef DATA_SET
# define DATA_SET set_all_mb475
#endif

void define_data_sets(Data_set data_set[NUM_SETS]);
void print_param(char *name, params p);
void print_param_pretty(FILE *fp, params p);
void read_param(char *name, params *p);

void four_neutrinos(int num, char **args);
void minimize(char *param_file_name);
void calc_chisq(char *param_file_name);

// calc bounds on d_e and d_mu in 4nu
void disapp_bounds(void);


void check_points(void);

void calc_points_pilar(void);

/***************************************************
 *      main
 ***************************************************/

int main(int num, char **args)
{
  // define the data sets 
  define_data_sets(data_set);
    
#ifdef READ_DATA
  if(num != 5 && num != 6){
    fprintf(stderr, "[main] READ_DATA: wrong num or arguments\n");
    fprintf(stderr, "[main] give data directory and name of 2 or 3 sets, and global. Available sets are:\n");
    for(int j = 0; j < NUM_SETS; j++)
      fprintf(stderr, "%s\n", data_set[j].name);
    exit(1);    
  }
#elif defined(MINIMIZE) || defined(CHISQ)
  if(num != 2)
    nrerror("[main] MINIMIZE: wrong num of arguments, give start point");
#else
//  if(num != 2)
//    nrerror("[main] wrong num of arguments, give output directory"); 
#endif  


  /**************** initialize the experiments **************************/

  //fprintf(stderr, "initializing...\n");

  // LSND/KARMEN ------------------------------
  initFlux();    
  initLSND();
  calcKarmen();

  // reactor experiments ----------------------

  // some defines:
  // fprintf(stderr, "defines: NO_CHOOZ, NO_SBL, NO_PV, NO_KAML, OLD_FLUX\n");
  // fprintf(stderr, "defines: SBL_ONLY = NO_CHOOZ + NO_PV + NO_KAML\n");

  init_fluxes(); 

#ifdef OLD_FLUX
  fprintf(stderr, "using old reactor fluxes\n");
#else
  fprintf(stderr, "using 2011 reactor fluxes\n");
#endif

#ifdef MINOS_NC
  fprintf(stderr, "using MINOS_NC for 0.1 < dmq < 1\n");
#else
  fprintf(stderr, "to use MINOS NC data define MINOS_NC\n");
#endif

#ifndef NO_CHOOZ
  chooz_init(old_new_main);
#endif  
  init_sbl_reactors(old_new_main);
  PV_init();
#ifndef BUGEY_TOTAL_RATE
  bugey_init(old_new_main);
#endif
#ifndef NO_KAML
  init_kaml();
#endif

  fit.invert_S();
  fit.pull_status[FLUX_NORM] = FIXED;

  // init other experiments ------------------------------
  initCDHS();
  initNomad();
  initMiniboone();
  initMBanti();
  
  //fprintf(stderr, "... end init\n");

  /******* end initialization *****************************/
  
  //calc_reactors();
  calc_points_pilar();
  

#ifdef CHECK_POINTS  
  check_points();
#endif  
  
#ifdef FOUR_NU
  four_neutrinos(num, args);
  return 0;
#endif

#ifdef MINIMIZE
  minimize(args[1]);
  return 0;
#endif

#ifdef CHISQ
  calc_chisq(args[1]);
  return 0;
#endif

#ifdef DISA_BOUNDS
  disapp_bounds();
  return 0;
#endif
  
#ifdef PLOT_SPECTRA
  if(num != 2)
    nrerror("[main] PLOT_SPECTRA: give parameter point");

  params pp;
  read_param(args[1], &pp);
  plot_excess_mb("spect.mb.osc.dat", pp);
  plot_excess_mba("spect.mba.osc.dat", pp);
  spectrum_lsnd("spect.lsnd.osc.dat", pp);
  return 0;
#endif

 /******************************************
  *  read minima and calculate PG
  ******************************************/ 
#ifdef READ_DATA
  
  const int num_sets = num - 2; // n sets + global (n = 2,3)
  
  // find the number of the data sets
  Read_Data_set sets[num_sets];
  
  int i = 0;
  for(int j = 0; j < NUM_SETS; j++){
    for(int k = 0; k < num_sets; k++){
    
      if(strcmp(args[2+k], data_set[j].name) == 0){
        sets[k].nr = j;
        i++;
      }
    }
  }
  if(i != num_sets){
    fprintf(stderr, "[main read] ERROR: could only find %d sets\n", i);
    exit(1);
  }

  // read the parameter values at the global min
  params p_min;
  char file_name[LENGTH];
  snprintf(file_name, LENGTH, "%s/%s/%s.point", 
	   PATH, args[1], data_set[sets[num_sets-1].nr].name);
  read_param(file_name, &p_min);
  
  // read the minima
  for(int k = 0; k < num_sets; k++){
    int sett, exp_tot;
  
    snprintf(file_name, LENGTH, "%s/%s/%s.chisq-min",
	     PATH, args[1], args[2+k]);
    FILE *fp = fopen(file_name, "r");
    if(fp == NULL){
      fprintf(stderr, "[main read] cannot open file for set%d", k+1);
      exit(1);
    }
    
    if(fscanf(fp, "%d %d %lf", &sett, &exp_tot, &sets[k].min) != 3){
      fprintf(stderr, "[main read] error in reading min for set%d", k+1);
      exit(1);
    }
    fclose(fp);
  
    if(sett != sets[k].nr){
      fprintf(stderr, "[main read] wrong set number for set%d", k+1);
      exit(1);
    }
    if(exp_tot != NUM_EXP)
      fprintf(stderr, "[main read] WARNING: num of experiments differ (set%d)", k+1);
  }
  
  // calculate PG
  double PG = 0.;
  double check = 0.;
  for(int k = 0; k < num_sets-1; k++){
    sets[k].min_gl = chisq_main(p_min, data_set[sets[k].nr].incl);
  
    printf("Set %d: %s:\tmin = %f,\tDelta chisq = %f\n",
	   k+1, data_set[sets[k].nr].name, 
	   sets[k].min, sets[k].min_gl - sets[k].min);
    
    PG += sets[k].min_gl - sets[k].min;
    check += sets[k].min;
  }
  
      
  printf("Global: %s:\tmin = %f,\t chisq_PG = %f\t\t(check: PG+min1+min2 = %f)\n", 
	 data_set[sets[num_sets-1].nr].name, sets[num_sets-1].min, 
	 PG, PG + check);
  
  printf("Best fit point:\n");
  print_param_pretty(stdout, p_min);

  /**********************************
   *  end read
   **********************************/ 
     
#else     

# ifndef Ip3pI
  fprintf(stderr, "Define Ip3pI for 1+3+1 spectrum\n");
# endif
  
  params p;
  p.Ue3 = 0.;

  fix_params f;  
  f.Ue[I4] = f.Um[I4] = false;
  f.Ue[I5] = f.Um[I5] = false;
  f.dmq[I4] = f.dmq[I5] = true;
  f.Ue3 = true;

  bool appearance = false;
  // only appearance exps
  if(DATA_SET == set_mb475 ||
     DATA_SET == set_mb300 ||
     DATA_SET == set_lsnd ||
     DATA_SET == set_l_k_n ||
     DATA_SET == set_nev_app ||
     DATA_SET == set_app_mb300 ||      
     DATA_SET == set_app_mb475){
    appearance = true;      
    f.Ue[I4] = f.Ue[I5] = true;
  }

# ifdef CP_SCAN
#  define N_CP 50
  fprintf(stderr, "scanning CP!\n");
  f.delta = true;
#  else
  f.delta = false;
#  endif

  double min = 1.e4;
  params p_min;

  // open the output file
  char file_name[LENGTH];
  snprintf(file_name, LENGTH, "%s/%s/%s.chisq-grid", PATH, args[1], data_set[DATA_SET].name);
  
  FILE *fp = fopen(file_name,"w");
  if(fp == NULL)
    nrerror("[main] cannot open file");

  double chisq[N_DMQ][N_DMQ];
 
  // the loop over the dmq's
  for(int i = 0; i < N_DMQ; i++)
  {
    p.dmq[I5] = DMQ_MIN * pow(10., i * DEL_LOG);

    for(int j = 0; j < N_DMQ; j++)
    {
      p.dmq[I4] = DMQ_MIN * pow(10., j * DEL_LOG);

      double chq = 1.e6;
      
      if(j <= i){

	params pp_min = p;
	
# ifdef CP_SCAN
        for(int d = 0; d < N_CP; d++){
	  p.delta = double(d) * 2. * M_PI / N_CP;
# else
	for(double del = 0.23 * M_PI; del < 2. * M_PI; del += 0.4 * M_PI)
        {
	  p.delta = del;
# endif
          p.Ue[I4] = p.Ue[I5] = (appearance ? 1./sqrt(2.) * 0.98 : 0.1);
          p.Um[I4] = p.Um[I5] = 0.01;
	  //if(chq > min + 30.)
          {
            const double cc = min_chisq(&p, f, data_set[DATA_SET].incl);
	    if(cc < chq){
	      chq = cc;
	      pp_min = p;
	    }
	  }
	}
	p = pp_min;
      }
      fprintf(fp, "%f  %f  %.3e\n", log10(p.dmq[I4]), log10(p.dmq[I5]), chq);

      chisq[i][j] = chisq[j][i] = chq;

      if(chq < min){
	min = chq;
	p_min = p;
      }
    }
    fflush(fp);
  }
  fclose(fp);

  /**********************
   * output
   **********************/

  // write the minimum
  snprintf(file_name, LENGTH, "%s/%s/%s.chisq-min", PATH, args[1], data_set[DATA_SET].name);
  fp = fopen(file_name, "w");
  fprintf(fp, "%d  %d  %f\n", DATA_SET, NUM_EXP, min);
  fclose(fp);

  // write the point
  snprintf(file_name, LENGTH, "%s/%s/%s.point", PATH, args[1], data_set[DATA_SET].name);
  print_param(file_name, p_min);

  // write again the full dmq plane
  snprintf(file_name, LENGTH, "%s/%s/%s.chisq-grid-sym", PATH, args[1], data_set[DATA_SET].name);
  fp = fopen(file_name,"w");

  // file for the projection
  char file_name2[LENGTH];
  snprintf(file_name2, LENGTH, "%s/%s/%s.chisq-proj", PATH, args[1], data_set[DATA_SET].name);
  FILE *fp2 = fopen(file_name2,"w");

  for(int i = 0; i < N_DMQ; i++){
    const double dmq1 = DMQ_MIN * pow(10., i * DEL_LOG);
    min = 1.e8;
    double min2 = 1.e8;

    for(int j = 0; j < N_DMQ; j++){
      const double dmq2 = DMQ_MIN * pow(10., j * DEL_LOG);

      fprintf(fp, "%.4e  %.4e  %.3e\n", log10(dmq1), log10(dmq2), chisq[i][j]);

      if(chisq[i][j] < min)
	min = chisq[i][j];
      if(j <= i && chisq[i][j] < min2)
	min2 = chisq[i][j];
    }
    fprintf(fp2, "%e  %.3e  %.3e\n", dmq1, min, min2);
  }
  fclose(fp);
  fclose(fp2);
#endif
  return 0;      /** end main **/
}



/*****************************************************
 *          four neutrinos
 *****************************************************/

#define N_LDMU 101

/* the chisq for disapp for 4nu */

double chi2disapp_3p1(double dm2, double s2)
{
  params p;  
  p.Ue[I4] = p.Um[I4] = p.delta = p.dmq[I4] = p.Ue3 = 0.;

  p.dmq[I5] = dm2;

  double min = 1.e8, w = 1.e8;

  if(s2 >= 1.){
    p.Ue[I5] = p.Um[I5] = 1./M_SQRT2;
    return chisq_main(p, data_set[set_dis].incl);
  
  }else{

    const double ldmu_min = log10(s2 / 4.);

    for(int i = 0; i < N_LDMU; i++){

      const double ldm = ldmu_min * (1. - i / (N_LDMU-1.));
      const double dm = pow(10., ldm);
      const double de = s2 / (4. * dm);

      p.Um[I5] = sqrt(dm);
      p.Ue[I5] = sqrt(de);

      if(de + dm <= 1.)
	w = chisq_main(p, data_set[set_dis].incl);

      if(w < min) min = w;
    }
  }
  return min;
}

#ifdef FOUR_NU

# define N_Ame 51
# define Ame_MIN 1.e-4
# define Ame_MAX 1.
# define DEL_LOG_A ( (log10(Ame_MAX) - log10(Ame_MIN))/(N_Ame-1.) )

void four_neutrinos(int num, char **args)
{
  if(num > NUM_EXP+1)
    nrerror("too many experiments");
  
  if(num < 2){
    fprintf(stderr, "give at least one experiment (if 7,8,9 is given all dissapp exps are included)\n");     
    fprintf(stderr, " 0: mb300\n");
    fprintf(stderr, " 1: mb475\n");
    fprintf(stderr, " 2: mba200\n");
    fprintf(stderr, " 3: mba475\n");
    fprintf(stderr, " 4: lsnd\n");
    fprintf(stderr, " 5: karmen\n");
    fprintf(stderr, " 6: nomad\n");
    fprintf(stderr, " 7: reactors\n");
    fprintf(stderr, " 8: cdhs\n");
    fprintf(stderr, " 9: atm\n");
    exit(0);
  }
  
  int expers[num-1];  
  for(int i = 1; i < num; i++){
     
    if(sscanf(args[i], "%d", &expers[i-1]) != 1)
      nrerror("cannot read arguents");
  }
   
   
  bool use_dis = false;
  bool inc[NUM_EXP];
  for(int i = 0; i < NUM_EXP; i++) inc[i] = false;
  for(int i = 0; i < num-1; i++){
     
    inc[expers[i]] = true;     
    if(expers[i] >= 7){
       use_dis = true;
       inc[expers[i]] = false;
    }     
  }
   
  params p;
  p.Ue[I4] = p.Um[I4] = p.delta = p.dmq[I4] = p.Ue3 = 0.;

#ifdef FOUR_NU_PG
  p.dmq[I5] = pow(10, -1.230285);
  const double Ame = pow(10, -8.000000e-02);

  p.Ue[I5] = p.Um[I5] = pow(Ame / 4., 0.25);
      
  double chisq = chisq_main(p, inc);       
  if(use_dis)
    chisq += chi2disapp_3p1(p.dmq[I5], Ame);

  printf("chisq = %e\n", chisq);
  exit(0);
#endif  

  //FILE *fp = fopen("out.dmq-new1.dat","w");

  for(int i = 0; i < N_DMQ; i++){
    p.dmq[I5] = DMQ_MIN * pow(10., i * DEL_LOG);

    double min = 1.e6;
    
    for(int j = 0; j < N_Ame; j++){      
      const double Ame = Ame_MIN * pow(10., j * DEL_LOG_A);

      p.Ue[I5] = p.Um[I5] = pow(Ame / 4., 0.25);
      
      double chisq = chisq_main(p, inc);       
      if(use_dis)
	 chisq += chi2disapp_3p1(p.dmq[I5], Ame);
      
      if(chisq < min) min = chisq;
      
      //fprintf(stdout, "%e  %e %e\n", log10(p.dmq[I5]), log10(Ame), chisq);
    }
    
    fprintf(stdout, "%e %e\n", p.dmq[I5], min);
    //fflush(fp);
    
    fflush(stdout);
  }
  //fclose(fp);
  return;
} 
#endif /** end four neutrinos **/


/*************************
 * minimize all params
 *************************/

void minimize(char *param_file_name)
{
  // read the parameter values at the global min
  params p;
  //read_param(param_file_name, &p);

  p.Ue[I4] = 0.1;
  p.Um[I4] = 0.01;
  p.Ue[I5] = 0.1;
  p.Um[I5] = 0.01;
  p.Ue3 = 0.;
  p.delta = M_PI/2.234;
  p.dmq[I4] = .9;
  p.dmq[I5] = 6.3;
  
  
  fix_params f;  
  f.Ue3 = true;
  f.Ue[I4] = f.Um[I4] = false;
  f.Ue[I5] = f.Um[I5] = false;
  f.delta  = false;
  f.dmq[I4] = f.dmq[I5] = true;

  /*** check minimizer ************************
  params pp = p;
 
  for(double u = 0.; u < 0.2; u += 0.002){
    
    pp.Ue[I4] = u;
    f.Ue[I4] = true;      
    const double che4 = min_chisq(&pp, f, data_set[DATA_SET].incl);
    pp = p;
    f.Ue[I4] = false;

    pp.Ue[I5] = u;
    f.Ue[I5] = true;      
    const double che5 = min_chisq(&pp, f, data_set[DATA_SET].incl);
    pp = p;
    f.Ue[I5] = false;
    
    pp.Um[I4] = u;
    f.Um[I4] = true;      
    const double chm4 = min_chisq(&pp, f, data_set[DATA_SET].incl);
    pp = p;
    f.Um[I4] = false;
    
    pp.Um[I5] = u;
    f.Um[I5] = true;      
    const double chm5 = min_chisq(&pp, f, data_set[DATA_SET].incl);
    pp = p;
    f.Um[I5] = false;
    
    printf("%e  %f %f %f %f\n", u, che4, che5, chm4, chm5);
    //printf("%e  %f \n", u, che4);
    fflush(stdout);
  }    
  return; ** end check ****************/

  /*
  printf("start point:\n");
  print_param_pretty(stdout, p);
  printf("chisq at start point = %f\n", chisq_main(p, data_set[DATA_SET].incl));

  plot_excess_mb("spect.mb.osc.dat", p);
  plot_excess_mba("spect.mba.osc.dat", p);
  spectrum_lsnd("spect.lsnd.osc.dat", p);
  
  return;
  */

  
  double chq = 1.e6;
     
  params pp_min = p;
	
# ifdef CP_SCAN
  f.delta = true;
  for(int d = 0; d < N_CP; d++){
	  p.delta = double(d) * 2. * M_PI / N_CP;
# else
  f.delta = false;
  for(double del = 0.23 * M_PI; del < 2. * M_PI; del += 0.4 * M_PI){
	  p.delta = del;
# endif
          p.Ue[I4] = p.Ue[I5] = 0.1;
          p.Um[I4] = p.Um[I5] = 0.01;
	  //if(chq > min + 30.)
          {
            const double cc = min_chisq(&p, f, data_set[DATA_SET].incl);
	    if(cc < chq){
	      chq = cc;
	      pp_min = p;
	    }
	  }
  }
  
  const double chq_min = chq; //min_chisq(&p, f, data_set[DATA_SET].incl);
  p = pp_min;

  printf("new minimum (saved in point.min):\n");
  print_param("point.min", p);
  print_param_pretty(stdout, p);
  printf("chisq_min = %f\n", chq_min);
    
  return;
}

/*************************
 * calculate the chisq
 *************************/

void calc_chisq(char *param_file_name)
{
  fix_params f;  
  params pp;

  pp.Ue[I4] = 0.1;
  pp.Um[I4] = 0.01;
  pp.Ue[I5] = 0.1;
  pp.Um[I5] = 0.01;
  pp.Ue3 = 0.;
  pp.delta = M_PI/2.234;
  pp.dmq[I4] = .1;
  
  f.Ue[I4] = false;
  f.Ue[I5] = false; 
  f.Um[I4] = false;
  f.Um[I5] = false;
  f.dmq[I4] = f.dmq[I5] = true;
  f.delta = false;
  f.Ue3 = true;
  
  for(int i = 0; i < N_DMQ; i++){
    pp.dmq[I5] = DMQ_MIN * pow(10., i * DEL_LOG);

    double chq = 1.e6;
      
    params pp_min = pp;
	
# ifdef CP_SCAN
    for(int d = 0; d < N_CP; d++){
	  pp.delta = double(d) * 2. * M_PI / N_CP;
# else
    for(double del = 0.23 * M_PI; del < 2. * M_PI; del += 0.4 * M_PI)
        {
	  pp.delta = del;
# endif
          pp.Ue[I4] = pp.Ue[I5] = 0.1;
          pp.Um[I4] = pp.Um[I5] = 0.01;
	  //if(chq > min + 30.)
          {
            const double cc = min_chisq(&pp, f, data_set[DATA_SET].incl);
	    if(cc < chq){
	      chq = cc;
	      pp_min = pp;
	    }
	  }
	}
    printf("%e  %e\n", pp.dmq[I5], chq);
    fflush(stdout);
  }
  return;
  

  // read the parameter values
  params p;
  //read_param(param_file_name, &p);
    
  bool appearance = false;
  // only appearance exps
  if(DATA_SET == set_mb475 ||
     DATA_SET == set_mb300 ||
     DATA_SET == set_lsnd ||
     DATA_SET == set_l_k_n ||
     DATA_SET == set_nev_app ||
     DATA_SET == set_app_mb300 ||      
     DATA_SET == set_app_mb475){
    appearance = true;      
  }

  f.Ue[I4] = f.Ue[I5] = appearance;   // fix for appearance exps

  p.Ue[I4] = p.Ue[I5] = (appearance ? 1./sqrt(2.) * 0.98 : 0.01);  
  p.Um[I4] = p.Um[I5] = 0.01;

  f.Um[I4] = f.Um[I5] = false;
  f.delta  = true;
#ifdef LOOP_DMQ
  f.dmq[I4] = f.dmq[I5] = true;
#else
  f.dmq[I4] = f.dmq[I5] = false;

# define N_START_P 6
  const double dmq_start[N_START_P] = {0.01, 0.46, 0.9, 1.9, 6.5, 12.};
#endif

  fprintf(stderr, "Using data set '%s'\n", data_set[DATA_SET].name);
    
  for(int d = 0; d < 51; d++){

    p.delta = double(d)/50. * 2. * M_PI;
    //p.delta = 1.2 * M_PI - double(d)/100. * .3 * M_PI;

    params pmin = p;
    double min = 1.e6;

#ifdef LOOP_DMQ
    // loop over the dmq's
    for(int i = 10; i < N_DMQ-10; i++){
      p.dmq[I5] = DMQ_MIN * pow(10., i * DEL_LOG);

      for(int j = 10; j <= i; j++){  // require Dmq51 >= Dmq 41
        p.dmq[I4] = DMQ_MIN * pow(10., j * DEL_LOG);
#else
    for(int j = 0; j < N_START_P; j++){
      p.dmq[I5] = dmq_start[j];

      for(int k = 0; k <= j; k++){
	p.dmq[I4] = dmq_start[k];

        p.Ue[I4] = p.Ue[I5] = (appearance ? 1./sqrt(2.) * 0.98 : 0.1);  
        p.Um[I4] = 0.1;
        p.Um[I5] = 0.1;
#endif
        const double chq = min_chisq(&p, f, data_set[DATA_SET].incl);

	if(chq < min){ 
	  min = chq;
	  pmin = p;
	}
      }
    }
    printf("%f  %f  %e  %e  %e  %e  %e\n", p.delta/M_PI, min, pmin.dmq[I4], pmin.dmq[I5],
                                           pmin.dmq[I5] - pmin.dmq[I4], 
	   pmin.Ue[I4] * pmin.Um[I4], pmin.Ue[I5] * pmin.Um[I5]);
    fflush(stdout);
  }
  return;
}

    
/*************************
 * check some chisq points
 *************************/
    
void check_points(void)
{
  params p_min_dis, p_min_all, p_min_app;
    
  char file_name[LENGTH];
#ifdef OLD_FLUX  
  snprintf(file_name, LENGTH, "Out-old/disa.point"); 
  read_param(file_name, &p_min_dis);
  snprintf(file_name, LENGTH, "Out-old/all_mb475.point"); 
  read_param(file_name, &p_min_all);
  snprintf(file_name, LENGTH, "Out-old/app_mb475.point"); 
  read_param(file_name, &p_min_app);
#else
# ifdef MINOS_NC
  snprintf(file_name, LENGTH, "Out-minos_nc/disa.point");
  read_param(file_name, &p_min_dis);
  snprintf(file_name, LENGTH, "Out-minos_nc/all_mb475.point"); 
  read_param(file_name, &p_min_all);
  snprintf(file_name, LENGTH, "Out-minos_nc/app_mb475.point"); 
  read_param(file_name, &p_min_app);
# else  
  snprintf(file_name, LENGTH, "Out-new/disa.point");
  read_param(file_name, &p_min_dis);
  snprintf(file_name, LENGTH, "Out-new/all_mb475.point"); 
  read_param(file_name, &p_min_all);
  snprintf(file_name, LENGTH, "Out-new/app_mb475.point"); 
  read_param(file_name, &p_min_app);
# endif  
#endif  
/*  
    (incl[mb300]  ? chi2mb300(p)   : 0.) +   
    (incl[mb475]  ? chi2mb475(p)   : 0.) +   
    (incl[mba200] ? chi2_MBA_200(p): 0.) +   
    (incl[mba475] ? chi2_MBA_475(p): 0.) +   
    (incl[lsnd]   ? chi2lsnd(p)    : 0.) +   
    (incl[karmen] ? chi2karmen(p)  : 0.) + 
    (incl[nomad]  ? chi2nomad(p)   : 0.) + 
    (incl[reactor]? chi2reactor(p) : 0.) +  
    (incl[cdhs]   ? chi2cdhs(p)    : 0.) +	
    (incl[atm]    ? chi2atm(p)     : 0.); 	
*/  
  
  double s_dis, s_all;
  printf("react: dis = %e, all = %e, delta = %e\n", chi2reactor(p_min_dis), chi2reactor(p_min_all)  , -chi2reactor(p_min_dis)+ chi2reactor(p_min_all)  );
  printf("cdhs:  dis = %e, all = %e, delta = %e\n", chi2cdhs(p_min_dis), chi2cdhs(p_min_all)        , -chi2cdhs(p_min_dis)+ chi2cdhs(p_min_all)        );
  printf("atm:   dis = %e, all = %e, delta = %e\n", chi2atm(p_min_dis), chi2atm(p_min_all)          , -chi2atm(p_min_dis)+ chi2atm(p_min_all)          );
  printf("lsnd:  app = %e, all = %e, delta = %e\n", chi2lsnd(p_min_app), chi2lsnd(p_min_all)        , -chi2lsnd(p_min_app)+ chi2lsnd(p_min_all)        );
  printf("karmen:app = %e, all = %e, delta = %e\n", chi2karmen(p_min_app), chi2karmen(p_min_all)    , -chi2karmen(p_min_app)+ chi2karmen(p_min_all)    );
  printf("mb:    app = %e, all = %e, delta = %e\n", chi2mb475(p_min_app), chi2mb475(p_min_all)      , -chi2mb475(p_min_app)+ chi2mb475(p_min_all)      );
  printf("mba:   app = %e, all = %e, delta = %e\n", chi2_MBA_475(p_min_app), chi2_MBA_475(p_min_all), -chi2_MBA_475(p_min_app)+ chi2_MBA_475(p_min_all));
  
  s_dis = chi2reactor(p_min_dis) + chi2cdhs(p_min_dis) + chi2atm(p_min_dis);
  s_all = chi2reactor(p_min_all) + chi2cdhs(p_min_all) + chi2atm(p_min_all);
  printf("sumdis:dis = %e, all = %e, delta = %e\n", s_dis, s_all, s_all - s_dis);
  
  exit(0);
}
    
    
/*************************
 * define sets
 *************************/


void define_data_sets(Data_set data_set[NUM_SETS])
{
  fprintf(stderr, "WARNING define sets: check for MBA before trusting a set!!\n");
  
  // NEV_APP
  for(int i = 0; i < NUM_EXP; i++)
    data_set[set_nev_app].incl[i] = (i == nomad || i == karmen);
  snprintf(data_set[set_nev_app].name, LENGTH, "nev_app");
  // LSND
  for(int i = 0; i < NUM_EXP; i++)
    data_set[set_lsnd].incl[i] = (i == lsnd);
  snprintf(data_set[set_lsnd].name, LENGTH, "lsnd");

  // Miniboone 300
  for(int i = 0; i < NUM_EXP; i++)
    data_set[set_mb300].incl[i] = (i == mb300);
  snprintf(data_set[set_mb300].name, LENGTH, "minib300");

  // Miniboone 475
  for(int i = 0; i < NUM_EXP; i++)
    data_set[set_mb475].incl[i] = (i == mb475);
  snprintf(data_set[set_mb475].name, LENGTH, "minib475");

  // NEV (without MB)
  for(int i = 0; i < NUM_EXP; i++)
    data_set[set_nev].incl[i] = (i != lsnd && i != mb300 && i != mb475);
  snprintf(data_set[set_nev].name, LENGTH, "noev");
  
  // NEV incl MB 475
  for(int i = 0; i < NUM_EXP; i++)
    data_set[set_nev_mb475].incl[i] = (i != lsnd && i != mb300 && i != mba475 && i != mba200);
  snprintf(data_set[set_nev_mb475].name, LENGTH, "noev_mb475");

  // EVID
  for(int i = 0; i < NUM_EXP; i++)
    data_set[set_evid].incl[i] = (i == lsnd || i == mba475);
  snprintf(data_set[set_evid].name, LENGTH, "evid");
  
  // LSND + KARMEN + NOMAD
  for(int i = 0; i < NUM_EXP; i++)
    data_set[set_l_k_n].incl[i] = (i == lsnd || i == karmen || i == nomad);
  snprintf(data_set[set_l_k_n].name, LENGTH, "L_K_N");
  // APP MB 300
  for(int i = 0; i < NUM_EXP; i++)
    data_set[set_app_mb300].incl[i] = (i == lsnd || i == karmen || i == nomad || i == mb300 || i == mba200);
  snprintf(data_set[set_app_mb300].name, LENGTH, "app_mb300");
  // APP MB 475
  for(int i = 0; i < NUM_EXP; i++)
    data_set[set_app_mb475].incl[i] = (i == lsnd || i == karmen || i == nomad || i == mb475 || i == mba475);
  snprintf(data_set[set_app_mb475].name, LENGTH, "app_mb475");
  // DIS
  for(int i = 0; i < NUM_EXP; i++)
    data_set[set_dis].incl[i] =  (i == reactor || i == cdhs || i == atm);
  //data_set[set_dis].incl[cdhs] = false;
  //fprintf(stderr, "WARNING: no cdhs in dis\n");
  snprintf(data_set[set_dis].name, LENGTH, "disa");
  // ALL MB 300
  for(int i = 0; i < NUM_EXP; i++)
    data_set[set_all_mb300].incl[i] = (i != mb475);
  data_set[set_all_mb300].incl[mba475] = false;  
  snprintf(data_set[set_all_mb300].name, LENGTH, "all_mb300");
  // ALL MB 475
  for(int i = 0; i < NUM_EXP; i++)
    data_set[set_all_mb475].incl[i] = (i != mb300);
  data_set[set_all_mb475].incl[mba200] = false;    
  //data_set[set_all_mb475].incl[cdhs] = false;
  //fprintf(stderr, "WARNING: no cdhs in all_mb475\n");
  snprintf(data_set[set_all_mb475].name, LENGTH, "all_mb475");

  // no LSND
  for(int i = 0; i < NUM_EXP; i++){
    data_set[set_nev_mb300].incl[i] = (i != lsnd && i != mb475);
    data_set[set_K_N_mb300].incl[i] = (i == karmen || i == nomad || i == mb300);
    data_set[set_K_N_mb475].incl[i] = (i == karmen || i == nomad || i == mb475);
  }
  
  snprintf(data_set[set_nev_mb300].name, LENGTH, "noev_mb300");
  snprintf(data_set[set_K_N_mb300].name, LENGTH, "K_N_mb300");
  snprintf(data_set[set_K_N_mb475].name, LENGTH, "K_N_mb475");
  
  return;
}

#ifdef DISA_BOUNDS

void disapp_bounds(void)
{
  params p;  
  p.Ue[I4] = p.Um[I4] = p.delta = p.dmq[I4] = p.Ue3 = 0.;


  for(int i = 0; i < 81; i++){
    p.dmq[I5] = exp10(-2. + 4. * i / 80.);
    
    for(int j = 0; j < 81; j++){
      const double d = exp10(-3. + 3. * j / 80.);
      
      printf("%e  %e  ", p.dmq[I5], d);
      
      p.Ue[I5] = sqrt(d);
      p.Um[I5] = 0.;
      printf("%f  ", chisq_main(p, data_set[set_dis].incl));
      
      p.Um[I5] = sqrt(d);
      p.Ue[I5] = 0.;
      printf("%f\n", chisq_main(p, data_set[set_dis].incl));
    }
  }
  return;
}

#endif


// caclulates the chisq for points provided by Pilar H.
void calc_points_pilar(void)
{
  params p;  
  p.Ue3 = 0.;

  // testing bf from tab. II of our paper:
  /*
  p.dmq[I4] = 0.47;
  p.Ue[I4] = 0.128;
  p.Um[I4] = 0.165;
  p.dmq[I5] = 0.87;
  p.Ue[I5] = 0.138;
  p.Um[I5] = 0.148;
  p.delta = 1.64*M_PI;

  printf("%f\n", chisq_main(p, data_set[DATA_SET].incl));
  */

  /*  
  p.dmq[I4] = 0.47;
  p.Ue[I4] = 1.e-8;
  p.Um[I4] = 1.e-8;
  p.dmq[I5] = 0.9;
  p.Ue[I5] = 1.e-8;
  p.Um[I5] = 1.e-8;
  p.delta = 1.62*M_PI;

  printf("%f\n", chisq_main(p, data_set[DATA_SET].incl));
  exit(0);
  */
  FILE *fp = fopen("Data-Pilar/NH_2.dat", "r");
  if(fp == NULL){
    fprintf(stderr, "cannot open data file\n");
    exit(0);
  }

  int k;
  for(int i = 0; i < 10000; i++){

    if(fscanf(fp, "%d %le %le %le  %le %le %le %le", &k,
	      &p.dmq[I4], &p.dmq[I5], &p.Ue[I4], &p.Ue[I5], 
	      &p.Um[I4], &p.Um[I5], &p.delta) != 8 || k != i){
      fprintf(stderr, "error in reading data file\n");
      exit(0);
    }

    printf("%d  %f\n", i, chisq_main(p, data_set[DATA_SET].incl));
  }

  fclose(fp);
  exit(0);
}

