// Neos Chi^2 calculation

// Standard Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// GSL libraries
#include <gsl/gsl_integration.h>        // Library for doing integrations
#include <gsl/gsl_sf.h>                 // Library for special functions (e.g. error function)
#include <gsl/gsl_matrix.h>     //
#include <gsl/gsl_linalg.h>     //
#include <gsl/gsl_permutation.h>    // Libraries for soving systems of linears equations

// Reactor library
#include "db-neos.h"


///----------------- Initialization ------------------------------------
// Neos data
void Neos_DB_class::read_Neos_data(const char* data, double x[], double y[], double z[], double w[], double t[], double r[]){
  int i;
  FILE *fin;
  
  fin = fopen(data,"r");
  if(fin==NULL){
	fprintf(stdout,"Error: file '%s' not found in Neos_DB_init() \n",data);
  }
    for (i=0; i<NBINS_NEOS_DB; i++){
      fscanf (fin, "%lf %lf %lf %lf %lf %lf ", &x[i],&y[i],&z[i],&w[i],&t[i],&r[i]);
    }
  fclose(fin);
}


// Neos ratio  Neos_obs/Neos_pred
void Neos_DB_class::read_Neos_ratio(const char* data, double x[], double y[], double z[], double t[]){
  int i;
  FILE *fin;
  
  fin = fopen(data,"r");
  if(fin==NULL){
	fprintf(stdout,"Error: file '%s' not found in Neos_DB_init() \n",data);
  }
    for(i=0;i<NBINS_NEOS_DB;i++){
	  fscanf(fin,"%lf %lf %lf %lf",&x[i],&y[i],&z[i],&t[i]);
    }
  fclose(fin);
}


// Daya Bay efficiencies 
void Neos_DB_class::read_eff(const char* data, double x[], double y[], double z[], double t[]){
  int i;
  FILE *fin;
  
  fin = fopen(data,"r");
  if(fin==NULL){
	fprintf(stdout,"Error: file '%s' not found in Neos_DB_init() \n",data);
  }
    for(i=0;i<8;i++){
      fscanf (fin,"%lf %lf %lf %lf", &x[i], &y[i], &z[i], &t[i]);
    }
  fclose(fin); 
}


// Daya Bay baselines
void Neos_DB_class::read_baselines_DB(const char* data, double x[][6]){
  int i,j;
  FILE *fin;
  
  fin = fopen(data,"r");
  if(fin==NULL){
	fprintf(stdout,"Error: file '%s' not found in Neos_DB_init() \n",data);
  }
    for(i=0;i<8;i++){		// we read for each detector 'i' the relative baseline to the reactor 'j'
      for(j=0;j<6;j++){		//
	    fscanf (fin,"%lf", &x[i][j]);
      }
    }
  fclose(fin); 
}


// Normalization Neos and DB
void Neos_DB_class::read_norms_Neos_DB(const char* data_neos, const char* data_DB, double x[], double y[]){
	int i;
	FILE *fin1,*fin2;
  
  // we read the normalization of Neos
    fin1=fopen(data_neos,"r");
    if(fin1==NULL){
	fprintf(stdout,"Error: file '%s' not found in Neos_DB_init() \n",data_neos);
    }
      fscanf(fin1,"%le",&x[0]);
    fclose(fin1);
  
  // Normalization of DayaBay
    fin2=fopen(data_DB,"r");
    if(fin2==NULL){
	fprintf(stdout,"Error: file '%s' not found in Neos_DB_init() \n",data_DB);
    }
      for(i=0;i<3;i++){
        fscanf(fin2,"%le ",&y[i]);
      }
    fclose(fin2);
}


// Neos integrals
void Neos_DB_class::read_Neos_integrals(const char* data1, const char* data2,double x1[], double x2[][M_MAX_NEOS_DB], double x3[]){
  int i,m;
  FILE *fin1,*fin2;
  
  fin1=fopen(data1,"r");
  if(fin1==NULL){
	fprintf(stdout,"Error: file '%s' not found in Neos_DB_init() \n",data1);
  }
    for(i=0;i<NBINS_NEOS_DB;i++){
      fscanf(fin1,"%le \n",&x1[i]);
    }
  fclose(fin1);

  fin2=fopen(data2,"r");
  if(fin2==NULL){
	fprintf(stdout,"Error: file '%s' not found in Neos_DB_init() \n",data2);
  }
    for(m=0;m<M_MAX_NEOS_DB;m++){
      for(i=0;i<NBINS_NEOS_DB;i++){
        fscanf(fin2,"%lf %lf \n",&x3[m],&x2[i][m]);
      }
    }
  fclose(fin2);	
  		   
}


// BayaBay integrals
void Neos_DB_class::read_DB_integrals(const char* data1, const char* data2, const char* data3, double x1[][NBINS_NEOS_DB], double x2[][2000], double x3[2000], double x4[NBINS_NEOS_DB]){
  int i,j;
  FILE *fin;
  
  fin=fopen(data1,"r");
  if(fin==NULL){
	fprintf(stdout,"Error: file '%s' not found in Neos_DB_init() \n",data1);
  }
    for(j=0;j<8;j++){
      for(i=0;i<NBINS_NEOS_DB;i++){
        fscanf(fin,"%lf \n",&x1[j][i]);
      }
    }
  fclose(fin);
  fin=fopen(data2,"r");
  if(fin==NULL){
	fprintf(stdout,"Error: file '%s' not found in Neos_DB_init() \n",data2);
  }
    for(i=0;i<NBINS_NEOS_DB;i++){
      for(j=0;j<2000;j++){
        fscanf(fin,"%le %le \n",&x3[j],&x2[i][j]);
      }
    }
  fclose(fin);
  
  fin=fopen(data3,"r");
  if(fin==NULL){
	fprintf(stdout,"Error: file '%s' not found in Neos_DB_init() \n",data3);
  }
    for(i=0;i<NBINS_NEOS_DB;i++){
      fscanf(fin,"%lf \n",&x4[i]);
    }
  fclose(fin);
}


// inverse covariance matrix
void Neos_DB_class::inverse_cov(const char* data,double x[][N_COV_NEOS_DB],double y[]){
  int i,j;
  int signum;
  FILE *fin;
  ///-----------------
  gsl_matrix * matrix = gsl_matrix_alloc(N_COV_NEOS_DB,N_COV_NEOS_DB);
  gsl_matrix * inv = gsl_matrix_alloc(N_COV_NEOS_DB,N_COV_NEOS_DB);
  gsl_permutation * per = gsl_permutation_alloc(N_COV_NEOS_DB);
  ///-----------------
  // we read the covariance matrix
  fin=fopen(data,"r");
  if(fin==NULL){
	fprintf(stdout,"Error: file '%s' not found in Neos_DB_init() \n",data);
  }
    for(i=0;i<N_COV_NEOS_DB;i++){
      for(j=0;j<N_COV_NEOS_DB;j++){
        fscanf(fin,"%lf ",&x[i][j]);
        if(i==j) x[i][j] = x[i][j] + pow(y[i],2);
        gsl_matrix_set(matrix, i, j, x[i][j]);
      }
    }
  fclose(fin);
  ///----------------
  gsl_linalg_LU_decomp(matrix,per,&signum);
  gsl_linalg_LU_invert(matrix,per,inv);
  ///----------------
  for(i=0;i<N_COV_NEOS_DB;i++){
    for(j=0;j<N_COV_NEOS_DB;j++){
      x[i][j] = gsl_matrix_get(inv,i,j);
    }
  }
  ///-----------------
  gsl_permutation_free(per);
  gsl_matrix_free(matrix);
  gsl_matrix_free(inv);
  ///-----------------
}
///---------------------------------------------------------------------


// Initialization Neos_DB ----------------------------------------------
void Neos_DB_class::Neos_DB_init(){
  int i;
  
  /// we read data from Neos
  read_Neos_data("in_files_Neos/neos_data.dat",Q,Q,Q,ON,bgN,RN);
  
  /// we read data of the ratio Neos_obs/Neos_pred arXiv:1610.05134v3
  read_Neos_ratio("in_files_Neos/ratio_neos.dat",Q,OR,erstat,ersys);
  
  /// we rad Daya Bay efficiencies 
  read_eff("in_files_Neos/efficiencies.dat",effmu,effm,DN,lifet);
  for(i=0;i<8;i++){
    ef[i]=effmu[i]*effm[i]*DN[i]*lifet[i];
  }
  
  /// we read DayaBay baselines 
  read_baselines_DB(DB_NEOS_PATH"in_files_Neos/baselines_reactor_vs_AD.dat",L_DB);
  
  /// we read the normalizations
  read_norms_Neos_DB(DB_NEOS_PATH"in_files_Neos/norm_neos.dat",
                     DB_NEOS_PATH"in_files_Neos/norms_DB.dat",norm_N,norm_DB);
  
  /// we read the Neos integrals
  read_Neos_integrals(DB_NEOS_PATH"integral_results_Neos/II_neos_61.dat",
                      DB_NEOS_PATH"integral_results_Neos/J41_neos_61.dat",II_N,J41_N,mvec);
  
  /// we read the DayaBay integrals
  read_DB_integrals(DB_NEOS_PATH"integral_results_Neos/II_DB_61.dat",
                    DB_NEOS_PATH"integral_results_Neos/JJJ_DB_61.dat",
                    DB_NEOS_PATH"integral_results_Neos/III_DB_61.dat",I_DB,J_DB,mL,II_DB_avout);  
  
  /// we construct the inverse covariance matrix
  inverse_cov(DB_NEOS_PATH"in_files_Neos/db-cov-matrix-for-neos.dat",cov,erstat);

//  printf("Neos initialized \n");
}
//----------------------------------------------------------------------


///----------- Chi2 calculation ----------------------------------------
double Neos_DB_class::chi2(void *params){

  int i,j,k;
  Param_5nu p = *(Param_5nu *) params;
    
  double s12,c12,s13,c13,s14,c14;
    s12 = pow(sin(p.theta[I12]),2);
    c12 = pow(cos(p.theta[I12]),2);
    s13 = pow(sin(p.theta[I13]),2);
    c13 = pow(cos(p.theta[I13]),2);
    s14 = pow(sin(p.theta[I14]),2);
    c14 = pow(cos(p.theta[I14]),2);
  
  double NN[NBINS_NEOS_DB],NN_3[NBINS_NEOS_DB];    
  /// Neos --------------
    for(i=0;i<NBINS_NEOS_DB;i++){	// Neos
      NN[i] = II_N[i]-4.0*s14*c14*int_Neos_intpol(p.dmq[3],i,mvec,J41_N);
      NN_3[i] = II_N[i];
    }
  ///--------------------  
  double J21[8],J31[8],J32[8],J41[8],J42[8],J43[8];
  double mx;
  double N1[NBINS_NEOS_DB],N2[NBINS_NEOS_DB],N1_3[NBINS_NEOS_DB],N2_3[NBINS_NEOS_DB];
  /// DayaBay -----------
    for(i=0;i<NBINS_NEOS_DB;i++){
      for(j=0;j<4;j++){		// Only the near detectors
        J21[j]=0.0;J31[j]=0.0;J32[j]=0.0;
        J41[j]=0.0;J42[j]=0.0;J43[j]=0.0;
        for(k=0;k<6;k++){
          mx=p.dmq[1]*L_DB[j][k];
	      J21[j] = J21[j] + int_DB_interpol(mx,i,mL,J_DB)/pow(L_DB[j][k],2);
	      mx=p.dmq[2]*L_DB[j][k];
	      J31[j] = J31[j] + int_DB_interpol(mx,i,mL,J_DB)/pow(L_DB[j][k],2);
	      mx=(p.dmq[2]-p.dmq[1])*L_DB[j][k];
	      J32[j] = J32[j] + int_DB_interpol(mx,i,mL,J_DB)/pow(L_DB[j][k],2);
	      mx=p.dmq[3]*L_DB[j][k];
	      J41[j] = J41[j] + int_DB_interpol(mx,i,mL,J_DB)/pow(L_DB[j][k],2);
	      mx=fabs(p.dmq[3]-p.dmq[1])*L_DB[j][k];
	      J42[j] = J42[j] + int_DB_interpol(mx,i,mL,J_DB)/pow(L_DB[j][k],2);
	      mx=fabs(p.dmq[3]-p.dmq[2])*L_DB[j][k];
	      J43[j] = J43[j] + int_DB_interpol(mx,i,mL,J_DB)/pow(L_DB[j][k],2);
        }
      }
 
	  N1[i]=0.0;N1_3[i]=0.0;
	  for(j=0;j<2;j++){	// EH1
		/**N1[i] = N1[i] + norm_DB[0]*ef[j]*(I_DB[j][i] - 4.0*(c14*c14*c12*c13*s12*c13*J21[j] + c14*c14*c13*s13*J31[j]
											          +c14*s14*c13*J41[j] + c14*s13*s14*J43[j]));
	    N1_3[i] = N1_3[i] + norm_DB[0]*ef[j]*(I_DB[j][i] - 4.0*c13*c13*c12*s12*J21[j] - 4.0*c13*s13*J31[j]);	// 3 neutrino scenario*/
	    N1[i] = N1[i] + norm_DB[0]*ef[j]*(I_DB[j][i] - 4.0*(c14*c14*c12*c13*s12*c13*J21[j] + c14*c14*c12*c13*s13*J31[j]
											               +c14*s14*c12*c13*J41[j] + c14*c14*s12*c13*s13*J32[j]
											               +c14*s14*s12*c13*J42[j] + c14*s13*s14*J43[j]));
	    N1_3[i] = N1_3[i] + norm_DB[0]*ef[j]*(I_DB[j][i] - 4.0*c13*c13*c12*s12*J21[j] - 4.0*c13*s13*(c12*J31[j]+s12*J32[j]));
      }
	  N2[i]=0.0;N2_3[i]=0.0;
	  for(j=2;j<4;j++){	// EH2
		/**N2[i] = N2[i] + norm_DB[1]*ef[j]*(I_DB[j][i] - 4.0*(c14*c14*c12*c13*s12*c13*J21[j] + c14*c14*c13*s13*J31[j]
											          +c14*s14*c13*J41[j] + c14*s13*s14*J43[j]));
	    N2_3[i] = N2_3[i] + norm_DB[1]*ef[j]*(I_DB[j][i] - 4.0*c13*c13*c12*s12*J21[j] - 4.0*c13*s13*J31[j]); // 3 neutrino scenario*/
	    N2[i] = N2[i] + norm_DB[1]*ef[j]*(I_DB[j][i] - 4.0*(c14*c14*c12*c13*s12*c13*J21[j] + c14*c14*c12*c13*s13*J31[j]
											               +c14*s14*c12*c13*J41[j] + c14*c14*s12*c13*s13*J32[j]
											               +c14*s14*s12*c13*J42[j] + c14*s13*s14*J43[j]));
	    N2_3[i] = N2_3[i] + norm_DB[1]*ef[j]*(I_DB[j][i] - 4.0*c13*c13*c12*s12*J21[j] - 4.0*c13*s13*(c12*J31[j]+s12*J32[j]));
      }
	}
  ///--------------------
  
  /// Chi2 calculation ----------------------------------------------  
    double chi2;
    chi2=0.0;
      for(i=0;i<N_COV_NEOS_DB;i++){
        for(j=0;j<N_COV_NEOS_DB;j++){
         chi2 = chi2 + (OR[i]-NN[i]/NN_3[i]/((N1[i]+N2[i])/(N1_3[i]+N2_3[i])))*cov[i][j]*(OR[j]-NN[j]/NN_3[j]/((N1[j]+N2[j])/(N1_3[j]+N2_3[j])));
	    }
      }
  ///----------------------------------------------------------------
  
  return chi2;
}
///---------------------------------------------------------------------


// Neos integrals interpolation
double Neos_DB_class::int_Neos_intpol(double a,int k,double x[],double y[][M_MAX_NEOS_DB]){

  int i;
  double m,n;
  
  if(a<=0.0001023) return 0.0;			// 10^(-4.0+0.01)
  if(a>=10.0) return II_DB_avout[k]/2.0;	// 10^(2.0)
  
  i=int((log10(a)+4.0)/0.01)-1;	// ESTE -1 ES SEGUN HE HECHO LAS INTEGRALES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  m=(y[k][i]-y[k][i+1])/(x[i]-x[i+1]);
  n=y[k][i]-(y[k][i]-y[k][i+1])/(x[i]-x[i+1])*x[i];

  return a*m+n; 
}

// DayaBay integrals interpolation
double Neos_DB_class::int_DB_interpol(double a,int k,double x[],double y[][2000]){
  
  int i;
  double m,n;
  
  if(a<=0.002) return 0.0;			// 10^(-1.7)
  if(a>=19952.62) return II_DB_avout[k]/2.0;	// 10^(4.3)

  i=int((log10(a)+1.7)/0.003);
  m=(y[k][i]-y[k][i+1])/(x[i]-x[i+1]);
  n=y[k][i]-(y[k][i]-y[k][i+1])/(x[i]-x[i+1])*x[i];
  
  return a*m+n;
}


#undef NBINS_NEOS_DB
#undef N_COV_NEOS_DB
#undef M_MAX_NEOS_DB
