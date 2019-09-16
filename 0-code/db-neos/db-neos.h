#define N_NU   5

#define DB_NEOS_PATH "db-neos/"

enum angle_index {XX, I12, I13, I14, I15};

// Parameter structure
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

  void save(const char *name){
    FILE *fp = fopen(name, "w");
    for(int i = 0; i < N_NU; i++)
      fprintf(fp, "%e %e %e\n", Ue[i], theta[i], dmq[i]);
    fclose(fp);
    return;
  }
  void load(const char *name){
    FILE *fp = fopen(name, "r");
    for(int i = 0; i < N_NU; i++)
      if(fscanf(fp, "%le %le %le\n", &Ue[i], &theta[i], &dmq[i]) != 3){
        fprintf(stderr, "cannot read parameter %s\n", name);
        exit(0);
      }
    fclose(fp);
    return;
  }
};


// DB only fit
#define NBINS_DB 35

class DB_class{
  public:
    
    void DB_init();
        double Emin[NBINS_DB],Emax[NBINS_DB],Ec[NBINS_DB];
        double O1[NBINS_DB],O2[NBINS_DB],O3[NBINS_DB],B1[NBINS_DB],B2[NBINS_DB],B3[NBINS_DB],Q[NBINS_DB];
      void read_DB_data(const char* data, double x[], double y[], double z[], double w[], double t[], double s[],
                        double a[], double b[], double c[], double d[], double e[], double f[]);
        double effmu[8],effm[8],DN[8],lifet[8],ef[8];
      void read_eff(const char* data, double x[], double y[], double z[], double t[]);
        double L_DB[8][6];
      void read_baselines_DB(const char* data, double x[][6]);
        double norm[3];
      void read_norm_DB(const char* data, double x[]);
        double II_DB_avout[NBINS_DB],I_DB[8][NBINS_DB],J_DB[NBINS_DB][2000],mL[2000];
	double IIe_DB_avout[NBINS_DB],Ie_DB[8][NBINS_DB],Je_DB[NBINS_DB][2000];
      void read_DB_integrals(const char* data1, const char* data2, const char* data3, double x1[][NBINS_DB], double x2[][2000], double x3[2000], double x4[NBINS_DB]);
        double cov[6][6];
      void systematics_DB(double x[][6]);
      
    double chi2(void *params);
    double chi2_syst(void *params);
      double int_DB_interpol(double a,int k,double x[],double y[][2000]);
};


// Daya Bay + NEOS combined fit
#define NBINS_NEOS_DB 61	// number of bins
#define N_COV_NEOS_DB 60	// number of used bins when computing the chi2
#define M_MAX_NEOS_DB 500

class Neos_DB_class{
  public:
	
  void Neos_DB_init();
	  double ON[NBINS_NEOS_DB],bgN[NBINS_NEOS_DB],RN[NBINS_NEOS_DB],Q[NBINS_NEOS_DB];
	void read_Neos_data(const char* data,double x[],double y[],double z[],double w[],double t[],double r[]);
	  double OR[NBINS_NEOS_DB],ersys[NBINS_NEOS_DB],erstat[NBINS_NEOS_DB];
	void read_Neos_ratio(const char* data, double x[], double y[], double z[], double t[]);
	  double effmu[8],effm[8],DN[8],lifet[8],ef[8];
	void read_eff(const char* data, double x[], double y[], double z[], double t[]);
	  double L_DB[8][6];
	void read_baselines_DB(const char* data, double x[][6]);
	  double norm_N[1],norm_DB[3]; 
	void read_norms_Neos_DB(const char* data_neos, const char* data_DB, double x[], double y[]);
	  double II_N[NBINS_NEOS_DB],J41_N[NBINS_NEOS_DB][M_MAX_NEOS_DB],mvec[M_MAX_NEOS_DB];
	void read_Neos_integrals(const char* data1, const char* data2, double x1[], double x2[][M_MAX_NEOS_DB], double x3[]);
      double II_DB_avout[NBINS_NEOS_DB],I_DB[8][NBINS_NEOS_DB],J_DB[NBINS_NEOS_DB][2000],mL[2000];
    void read_DB_integrals(const char* data1, const char* data2, const char* data3, double x1[][NBINS_NEOS_DB], double x2[][2000], double x3[2000], double x4[NBINS_NEOS_DB]);
      double cov[N_COV_NEOS_DB][N_COV_NEOS_DB];
    void inverse_cov(const char* data,double x[][N_COV_NEOS_DB],double y[]);
    
    
  double chi2(void *params);
    double int_Neos_intpol(double a,int k,double x[],double y[][M_MAX_NEOS_DB]);
    double int_DB_interpol(double a,int k,double x[],double y[][2000]);
    
};



