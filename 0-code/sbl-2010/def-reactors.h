#include "definitions.h"

#define EnuMIN 1.8
#define EnuMAX 12.

#define NO_KAML
//#define BUGEY_TOTAL_RATE

#ifdef SBL_ONLY
# define NO_CHOOZ
# define NO_KAML
# define NO_PV
#endif


#ifdef BUGEY_TOTAL_RATE
enum SBL_Reactors {BUGEY3_15, BUGEY3_40, BUGEY3_95, BUGEY4, ROVNO, 
		   GOSGEN_1, GOSGEN_2, GOSGEN_3, ILL, 
		   KRASN_1, KRASN_2, KRASN_3, N_SBL_R};
#else
enum SBL_Reactors {BUGEY4, ROVNO, 
		   GOSGEN_1, GOSGEN_2, GOSGEN_3, ILL, 
		   KRASN_1, KRASN_2, KRASN_3, N_SBL_R};
#endif

#ifndef NO_CHOOZ
# define NBIN_CHOOZ 14
#else
# define NBIN_CHOOZ 0
#endif

// Bugey
#ifdef BUGEY_TOTAL_RATE
  #define NBIN_BUG 0
#else
  #define NBIN_BUG 60
#endif


/* uncertainties on fluxes include just a normalization error for each isotope
 * (no shape uncertainties included)
 * FLUX_NORM includes in addition an overall flux uncertainty
 */
#ifdef BUGEY_TOTAL_RATE
enum Pulls {PULL_U235, PULL_U238, PULL_P239, PULL_P241, FLUX_NORM, 
            NORM_CHOOZ, PULL_KAML_0};
#else
enum Pulls {PULL_U235, PULL_U238, PULL_P239, PULL_P241, FLUX_NORM, 
            NORM_CHOOZ, E_SCALE_BUG, PULL_KAML_0};
#endif

#define KAML_NREACT 21

#ifndef NO_KAML
# define NBIN_KAML 17
# define NPULL_KAML 109
# define PULL_E_SCALE 107 
#else
# define NBIN_KAML 0
# define NPULL_KAML 0
# define PULL_E_SCALE 0 
#endif

/* KamLand pulls:                       */
/* 1              detector              */
/* NREACT         powers                */
/* NREACT*NISO    isotope compositions  */
/* 1              background            */
/* 1              energy scale          */
/* 1              geo nus               */



#define NPULLS (PULL_KAML_0 + NPULL_KAML)

// "1" for Palo Verde
#define NBIN_CHISQ (NBIN_CHOOZ + N_SBL_R + 1 + NBIN_BUG + NBIN_KAML)
#define BIN_KAML_0 (NBIN_CHOOZ + N_SBL_R + 1 + NBIN_BUG)
#define BIN_BUG_0  (NBIN_CHOOZ + N_SBL_R + 1)


#define DMQ_31 2.5e-3
#define DMQ_21 7.7e-5


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

/************************/
/* function definitions */
/************************/

double crossSect(double eKin);   // defined in main.cc

bool invert(const int size, void *B_in, void *C_out);
bool singsolve(const int size, void *B_in, void *L_out);

/******** Chooz **************/
void chooz_init(const int old_new);
void set_table_chooz(params &prm, double cff[NBIN_CHISQ][NPULLS+1]);

/******** Bugey **************/
void bugey_init(const int old_new);
void set_table_bugey(params &prm, double cff[NBIN_CHISQ][NPULLS+1]);

/******** Palo Verde **************/
void PV_init(void);
void set_table_PV(params &prm, double cff[NBIN_CHISQ][NPULLS+1]);

/******** SBL reactors *******/
void init_sbl_reactors(int old_new);
void set_table_sbl(params &prm, double cff[NBIN_CHISQ][NPULLS+1]);
void print_sbl_data(void);
void print_sbl_pred(double xi[NPULLS]);

/********** KamLAND **********/
void init_kaml(void);  
void set_table_kaml(const params &prm, double cff[NBIN_CHISQ][NPULLS+1]);

/************************/
/* fluxes               */
/************************/

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

#define OLD 0
#define NEW 1


/*************/
/* class Fit */
/*************/

enum Pull_status {ACTIVE, FIXED, FREE};

class Fit {
 public:
   Fit(void);
   void invert_S(void);

   // if(use_xi_in) the fixed pull values from xi are used
   // else fixed pulls are set to zero
   double chisq(params &prm, double xi[NPULLS], bool use_xi_in = false);
   double chisq(params &prm){
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




