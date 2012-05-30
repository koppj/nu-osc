#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "qrom.h"

#define EnuMIN 1.8
#define EnuMAX 10.

#define USE_SBL
#define USE_BUGEY_SP

#ifdef USE_ALL
# define USE_SBL
# define USE_CHOOZ
# define USE_PV
# define USE_KAML
# define USE_DC
# define USE_DB
# define USE_RENO
# define USE_BUGEY_SP
#endif

enum EXPERIMENTS {SBL, CHOOZ, PV, KAML, DC, DB, RENO, BUG_SP, N_EXP};

enum SBL_Reactors {BUGEY4, ROVNO, BUGEY3_1, BUGEY3_2, BUGEY3_3, 
		   GOSGEN_1, GOSGEN_2, GOSGEN_3, ILL, 
		   KRASN_1, KRASN_2, KRASN_3, 
                   SRP_1, SRP_2, 
                   ROV_1I, ROV_2I, ROV_1S, ROV_2S, ROV_3S,
                   N_SBL_R};

/* defining the number of bins and private pulls for each experiment */
#ifdef USE_SBL
# define NBIN_SBL N_SBL_R
#else
# define NBIN_SBL 0
#endif
#define NPULL_SBL 0

#ifdef USE_CHOOZ
# define NBIN_CHOOZ 14
# define NPULL_CHOOZ 1
#else
# define NBIN_CHOOZ 0
# define NPULL_CHOOZ 0
#endif

#ifdef USE_PV
# define NBIN_PV 1
#else
# define NBIN_PV 0
#endif
#define NPULL_PV 0

#ifdef USE_KAML
# define NBIN_KAML 17
# define NPULL_KAML 109
#else
# define NBIN_KAML 0
# define NPULL_KAML 0
#endif
#define KAML_NREACT 21
#define PULL_E_SCALE 107 

/* KamLand pulls:                       */
/* 1              detector              */
/* NREACT         powers                */
/* NREACT*NISO    isotope compositions  */
/* 1              background            */
/* 1              energy scale          */
/* 1              geo nus               */

#ifdef USE_DC
# define NBIN_DC 18
# define NPULL_DC 5 
#else
# define NBIN_DC 0
# define NPULL_DC 0
#endif

#ifdef USE_DB
# define NBIN_DB 6
# define NPULL_DB 19 // 1(norm) + Nreact + Ndet + Ndet(bg) 
#else
# define NBIN_DB 0
# define NPULL_DB 0
#endif

#ifdef USE_RENO
# define NBIN_RENO 2
# define NPULL_RENO 11 // 1(norm) + Nreact + Ndet + Ndet(bg) 
#else
# define NBIN_RENO 0
# define NPULL_RENO 0
#endif

#ifdef USE_BUGEY_SP
# define NBIN_BUG_SP 60
# define NPULL_BUG_SP 1
#else
# define NBIN_BUG_SP 0
# define NPULL_BUG_SP 0
#endif

/* total bins in chisq */
#define NBIN_CHISQ (NBIN_SBL + NBIN_CHOOZ + NBIN_PV + NBIN_KAML + NBIN_DC + NBIN_DB + NBIN_RENO + NBIN_BUG_SP)


#define N_CO_UNC 3   // 2nd order polynomial for bin-to-bin uncorrelated flux error

/* global pulls common to all experiments */
enum Pulls {
  PULL_U235_0, PULL_U235_1, PULL_U235_2,     // polynomial for shape uncertainty
  PULL_P239_0, PULL_P239_1, PULL_P239_2,     // polynomial for shape uncertainty
  PULL_P241_0, PULL_P241_1, PULL_P241_2,     // polynomial for shape uncertainty
  FLUX_COR,                                  // flux error correlated btw U235, P239, P241
  PULL_U238,                                 // norm for U238
  FLUX_NORM,                                 // overall flux uncertainty
  PULL_GLOBAL
};


/* total pulls in chisq */
#define NPULLS (PULL_GLOBAL + NPULL_SBL + NPULL_CHOOZ + NPULL_PV + NPULL_KAML + NPULL_DC + NPULL_DB + NPULL_RENO + NPULL_BUG_SP)



/************************
 *  range of Dmq's      *
 ************************/

// logarithm of eV-scale dmqs
#define STE_MAX 1.5

#define LOG_ATM

#ifdef LOG_ATM
# define ATM_MIN (-3.5)
# define ATM_MAX (-1.5)
#else
# define ATM_MIN (0.)
# define ATM_MAX (4.e-3)
#endif

#define DM2_NUM 201
#define D_ATM ( (ATM_MAX - ATM_MIN)/(DM2_NUM-1.) )

#ifdef LOG_SOL
# define SOL_MIN  (-5.) 
# define SOL_MAX  (-4.) 
#else
# define SOL_MIN  (0.) 
# define SOL_MAX  (1.8e-4) 
#endif

#define DM2_NUM_21 301 
#define D_SOL ((SOL_MAX-SOL_MIN) / (DM2_NUM_21 - 1.))


/************************
 *  constants           *
 ************************/

#define MN 939.566
#define MP 938.272
#define DELTA (MN-MP)
#define ME 0.511

// fluxes
#define OLD 0
#define NEW 1

/********************/
/* Inline functions */
/********************/

inline double norm(double x) {
   return x*x;
}

inline void error(const char *text) {
   fputs(text, stderr);
   fputc('\n', stderr);
   exit(1);
}

inline double exp10(const double x){
  return pow(10, x);
}

/********************/
/* the parameters   */
/********************/

#define N_NU 5

enum angle_index {XX, I12, I13, I14, I15};

struct Param_5nu{
  double Ue[N_NU];
  double theta[N_NU];   // theta[0] not used, theta[i] = theta_{0i}

  void set_ang(void){
    for(int i = 0; i < N_NU; i++) 
      Ue[i] = 1.;

    for(int i = 1; i < N_NU; i++){
      for(int j = 0; j < i; j++)
	Ue[j] *= cos(theta[i]);
      Ue[i] *= sin(theta[i]);
    }
  }

  double dmq[N_NU];    // Dm^2_{j0} j = 0,...,N_NU-1
  double Dmq(const int j, const int i){
    return (j == 0 ? 0. : dmq[j]) - (i == 0 ? 0. : dmq[i]);
  }
};

struct Parameters{
  double sq12;
  double sq13;

  double sq2_t13(void){
    return 4.*sq13*(1.-sq13);
  }

  int idmq_21;
  int idmq_31;

  double dmq31(void){
#ifdef LOG_ATM
    return exp10(ATM_MIN + D_ATM * idmq_31);
#else
    return ATM_MIN + D_ATM * idmq_31;
#endif
  }
  double dmq21(void) {
#ifdef LOG_SOL
    return exp10(SOL_MIN + D_SOL * idmq_21);
#else
    return SOL_MIN + D_SOL * idmq_21;
#endif
  }
}; 


/************************/
/* function definitions */
/************************/

double read(const char *name, double x0);

// defined in class_flux.cc
double crossSect(double eKin);   
void init_fluxes(void);

bool invert(const int size, void *B_in, void *C_out);
bool singsolve(const int size, void *B_in, void *L_out);

/******** Bugey spectrum **************/
void bugey_init(int old_new = NEW);
void set_table_bugey(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1]);
/******** RENO **************/
void RENO_init(void);
void set_table_RENO(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1]);
/******** DB **************/
void DB_init(void);
void set_table_DB(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1]);
/******** DC **************/
void dc_init(const int old_new);
void set_table_dc(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1]);
void plot_dc_pred(Parameters p, double f);
void plot_dc_data(void);
double dc_test_chisq(Parameters prm);

/******** Chooz **************/
void chooz_init(const int old_new);
void set_table_chooz(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1]);

/******** Palo Verde **************/
void PV_init(void);
void set_table_PV(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1]);

/******** SBL reactors *******/
int init_sbl_reactors(int old_new);
void set_table_sbl(Param_5nu &p, double cff[NBIN_CHISQ][NPULLS+1]);
void print_sbl_data(void);
void print_sbl_pred(double xi[NPULLS]);

/********** KamLAND **********/
void init_kaml(void);  
void set_table_kaml(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1]);

/************************/
/* fluxes               */
/************************/

// flux function used by the experiments
double global_flux(const int iso, const double Enu, const int old_new);


// class for polynomial fit to fluxes
#define MAX_COEF 6

class Flux {
 private:
  int Ncoef;
  double coef[MAX_COEF];

 public:
  void init(const char *filename);
  double f(double Enu);
};

enum Isotopes {U235, U238, P239, P241, NISO};


/* struct four rate pull couplings */

struct Rate_pull_coupl{
  double unc1[NISO][N_CO_UNC]; 
  double cor1[NISO];

  double unc(int iso, int coef, const double f[NISO])
  {
    if(iso == U238 && coef != 0)
      error("rate_pull_coupl: for U238 only coef = 0 can be used!\n");
    double norm = 0;
    for(int i = 0; i < NISO; i++)
      norm += f[i] * unc1[i][0];
    return f[iso] * unc1[iso][coef] / norm;
  }
  double cor(const double f[NISO])
  {
    double norm = 0;
    for(int i = 0; i < NISO; i++)
      norm += f[i] * unc1[i][0];

    return (f[U235] * cor1[U235] + f[P239] * cor1[P239] + f[P241] * cor1[P241]) / norm;
  }
};

/* probabilities of rate experiments */
/* defined in class_flux.cc          */

#define N_RATE_COEF 2001

class Rate_coef{
 private:
  double DmqL_min, DmqL_max, log_DmqL_min, log_DmqL_max, Del;
  double Sinq[NISO][N_RATE_COEF];
  double n0[NISO];

 public:
  void init(void);
  double sinq(const double DmqL, const double f[NISO]);
  double prob(Param_5nu &p, const double L, const double f[NISO]);
};


/*************/
/* class Fit */
/*************/

enum Pull_status {ACTIVE, FIXED, FREE};

class Fit {
 public:
   Fit(void);
   void invert_S(void);
   int first_bin[N_EXP];  
   int first_pull[N_EXP];  

   // if(use_xi_in) the fixed pull values from xi are used
   // else fixed pulls are set to zero
   double chisq(Param_5nu &prm, double xi[NPULLS], bool use_xi_in = false);
   double chisq(Param_5nu &prm){
     double xi[NPULLS];
     return chisq(prm, xi);
   }

   double Data [NBIN_CHISQ];
   double S_data[NBIN_CHISQ][NBIN_CHISQ]; 
   double S_pull[NPULLS][NPULLS];
   double S_data_inv[NBIN_CHISQ][NBIN_CHISQ]; 
   double S_pull_inv[NPULLS][NPULLS];
   int pull_status[NPULLS];
};




