// DB Chi^2 calculation

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
// Daya Bay data in each Experiemntal Hall (EH)
void DB_class::read_DB_data(const char* data, double x[], double y[], double z[], double w[], double t[], double s[],
                  double a[], double b[], double c[], double d[], double e[], double f[]){
  int i;
  FILE *fin;
  
  fin = fopen(data,"r");
  if(fin==NULL){
    fprintf(stdout,"Error: file '%s' not found in DB_init() \n",data);
  }
    for (i=0; i<NBINS_DB; i++){
      fscanf (fin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
	      &x[i],&y[i],&z[i],&w[i],&t[i],&s[i],&a[i],&b[i],&c[i],&d[i],&e[i],&f[i]);
    }
  fclose(fin);
}


// Daya Bay efficiencies 
void DB_class::read_eff(const char* data, double x[], double y[], double z[], double t[]){
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
void DB_class::read_baselines_DB(const char* data, double x[][6]){
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

// Daya Bay normalization
void DB_class::read_norm_DB(const char* data, double x[]){
  int i;
  FILE *fin;
    fin=fopen(data,"r");
    if(fin==NULL){
	fprintf(stdout,"Error: file '%s' not found in Neos_DB_init() \n",data);
    }
      for(i=0;i<3;i++){
        fscanf(fin,"%le ",&x[i]);
      }
    fclose(fin);
}

// BayaBay integrals
void DB_class::read_DB_integrals(const char* data1, const char* data2, const char* data3, double x1[][NBINS_DB], double x2[][2000], double x3[2000], double x4[NBINS_DB]){
  int i,j;
  FILE *fin;
  
  fin=fopen(data1,"r");
  if(fin==NULL){
	fprintf(stdout,"Error: file '%s' not found in Neos_DB_init() \n",data1);
  }
    for(j=0;j<8;j++){
      for(i=0;i<NBINS_DB;i++){
        fscanf(fin,"%lf \n",&x1[j][i]);
      }
    }
  fclose(fin);
  fin=fopen(data2,"r");
  if(fin==NULL){
	fprintf(stdout,"Error: file '%s' not found in Neos_DB_init() \n",data2);
  }
    for(i=0;i<NBINS_DB;i++){
      for(j=0;j<2000;j++){
        fscanf(fin,"%le %le \n",&x3[j],&x2[i][j]);
      }
    }
  fclose(fin);
  
  fin=fopen(data3,"r");
  if(fin==NULL){
	fprintf(stdout,"Error: file '%s' not found in Neos_DB_init() \n",data3);
  }
    for(i=0;i<NBINS_DB;i++){
      fscanf(fin,"%lf \n",&x4[i]);
    }
  fclose(fin);
}


// correlation matrix for the energy scale and detector efficiencies
void DB_class::systematics_DB(double cov[][6]){
  int k,l;
  int signum,nx=6;
  
  for(k=0;k<nx;k++){
    for(l=0;l<nx;l++){
      cov[k][l]=0.0;
    }
  }
  cov[0][0]=pow(0.0013,2)/2.0;	    //
  cov[1][1]=pow(0.0013,2)/2.0;		//
  cov[2][2]=pow(0.0013,2)/4.0;		// Efficiency error
  cov[3][3]=pow(0.002,2)/2.0;  //
  cov[4][4]=pow(0.002,2)/2.0;  //
  cov[5][5]=pow(0.002,2)/4.0;  // Relative energy scale
  // correlations
    cov[0][3]=0.54*0.0013*0.002/2.0;
    cov[1][4]=0.54*0.0013*0.002/2.0;
    cov[2][5]=0.54*0.0013*0.002/4.0;
    cov[3][0]=0.54*0.0013*0.002/2.0;
    cov[4][1]=0.54*0.0013*0.002/2.0;
    cov[5][2]=0.54*0.0013*0.002/4.0;
    
  ///-----------------------------
  gsl_matrix * B = gsl_matrix_alloc(nx,nx);
  gsl_matrix * Binv = gsl_matrix_alloc(nx,nx);
  gsl_permutation * perB = gsl_permutation_alloc(nx);    
  for(k=0;k<nx;k++){
    for(l=0;l<nx;l++){
      gsl_matrix_set(B,k,l,cov[k][l]);
    }
  }

  gsl_linalg_LU_decomp(B, perB, &signum);
  gsl_linalg_LU_invert (B, perB, Binv);
  
  for(k=0;k<nx;k++){
    for(l=0;l<nx;l++){
      cov[k][l]=gsl_matrix_get(Binv,k,l);
    }
  }
  gsl_permutation_free(perB);
  gsl_matrix_free(B);
  gsl_matrix_free(Binv);  
  ///-----------------------------
}
///---------------------------------------------------------------------


// Initialization DB ---------------------------------------------------
void DB_class::DB_init(){
  int i;
  
  // we read the data of the EH's in DB
  read_DB_data(DB_NEOS_PATH"in_files_DB/EH1.dat",Emin,Emax,Ec,O1,Q,Q,B1,Q,Q,Q,Q,Q);
  read_DB_data(DB_NEOS_PATH"in_files_DB/EH2.dat",Q,Q,Q,O2,Q,Q,B2,Q,Q,Q,Q,Q);
  read_DB_data(DB_NEOS_PATH"in_files_DB/EH3.dat",Q,Q,Q,O3,Q,Q,B3,Q,Q,Q,Q,Q);
  
  /// we rad Daya Bay efficiencies 
  read_eff(DB_NEOS_PATH"in_files_DB/efficiencies.dat",effmu,effm,DN,lifet);
  for(i=0;i<8;i++){
    ef[i]=effmu[i]*effm[i]*DN[i]*lifet[i];
  }
  
  /// we read DayaBay baselines 
  read_baselines_DB(DB_NEOS_PATH"in_files_DB/baselines_reactor_vs_AD.dat",L_DB);
  
  /// we read normalization of DayaBay
  read_norm_DB(DB_NEOS_PATH"in_files_DB/norms_DB.dat",norm);
  
  /// we read the DayaBay integrals
  read_DB_integrals(DB_NEOS_PATH"integrals_DB/II.dat",
                    DB_NEOS_PATH"integrals_DB/JJJ.dat",
                    DB_NEOS_PATH"integrals_DB/III.dat",I_DB,J_DB,mL,II_DB_avout);  
  read_DB_integrals(DB_NEOS_PATH"integrals_DB/II_es.dat",
                    DB_NEOS_PATH"integrals_DB/JJJ_es.dat",
                    DB_NEOS_PATH"integrals_DB/III_es.dat",Ie_DB,Je_DB,mL,IIe_DB_avout);
 
  systematics_DB(cov);
  
//  printf("Daya Bay initialized \n");
}
//----------------------------------------------------------------------


///----------- Chi2 calculation ----------------------------------------
double DB_class::chi2(void *params){	// Chi^2 without systematics
  int i,j,k;
  Param_5nu p = *(Param_5nu *) params;
  
  double s12,c12,s13,c13,s14,c14;
    s12 = pow(sin(p.theta[I12]),2);
    c12 = pow(cos(p.theta[I12]),2);
    s13 = pow(sin(p.theta[I13]),2);
    c13 = pow(cos(p.theta[I13]),2);
    s14 = pow(sin(p.theta[I14]),2);
    c14 = pow(cos(p.theta[I14]),2);
    
  double J21[8],J31[8],J32[8],J41[8],J42[8],J43[8];
  double mx;
  double N1[NBINS_DB],N2[NBINS_DB],N3[NBINS_DB];
  /// DayaBay -----------
    for(i=0;i<NBINS_DB;i++){
      for(j=0;j<8;j++){
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
 
	  N1[i]=0.0;
	  for(j=0;j<2;j++){	// EH1
	    N1[i] = N1[i] + norm[0]*ef[j]*(I_DB[j][i] - 4.0*(c14*c14*c12*c13*s12*c13*J21[j] + c14*c14*c12*c13*s13*J31[j]
							                            +c14*s14*c12*c13*J41[j] + c14*c14*s12*c13*s13*J32[j]
							                            +c14*s14*s12*c13*J42[j] + c14*s13*s14*J43[j]));
      }
	  N2[i]=0.0;
	  for(j=2;j<4;j++){	// EH2
	    N2[i] = N2[i] + norm[1]*ef[j]*(I_DB[j][i] - 4.0*(c14*c14*c12*c13*s12*c13*J21[j] + c14*c14*c12*c13*s13*J31[j]
							                            +c14*s14*c12*c13*J41[j] + c14*c14*s12*c13*s13*J32[j]
							                            +c14*s14*s12*c13*J42[j] + c14*s13*s14*J43[j]));
      }
	  N3[i]=0.0;
	  for(j=4;j<8;j++){	// EH3
	    N3[i] = N3[i] + norm[2]*ef[j]*(I_DB[j][i] - 4.0*(c14*c14*c12*c13*s12*c13*J21[j] + c14*c14*c12*c13*s13*J31[j]
							                            +c14*s14*c12*c13*J41[j] + c14*c14*s12*c13*s13*J32[j]
							                            +c14*s14*s12*c13*J42[j] + c14*s13*s14*J43[j]));
	  }
	}

  // Chi2 calculation ----------------------------------------------  
    double chi2;
    chi2=0.0; 
      for(i=0;i<NBINS_DB;i++){
        chi2 = chi2 + pow((O3[i]-B3[i])/(O1[i]-B1[i])-N3[i]/N1[i],2)
                         /(O3[i]/pow(O1[i]-B1[i],2)+O1[i]*pow(O3[i]-B3[i],2)/pow(O1[i]-B1[i],4))
                         + pow((O2[i]-B2[i])/(O1[i]-B1[i])-N2[i]/N1[i],2)
                         /(O2[i]/pow(O1[i]-B1[i],2)+O1[i]*pow(O2[i]-B2[i],2)/pow(O1[i]-B1[i],4));
      }
  ///----------------------------------------------------------------
  
  return chi2;
}


double DB_class::chi2_syst(void *params){	// Chi^2 without systematics
  int i,j,k,l;
  int nx=6;
  int signum;
  Param_5nu p = *(Param_5nu *) params;
  
  double s12,c12,s13,c13,s14,c14;
    s12 = pow(sin(p.theta[I12]),2);
    c12 = pow(cos(p.theta[I12]),2);
    s13 = pow(sin(p.theta[I13]),2);
    c13 = pow(cos(p.theta[I13]),2);
    s14 = pow(sin(p.theta[I14]),2);
    c14 = pow(cos(p.theta[I14]),2);
    
  double J21[8],J31[8],J32[8],J41[8],J42[8],J43[8];
  double J21e[8],J31e[8],J32e[8],J41e[8],J42e[8],J43e[8];
  double mx;
  double N1[NBINS_DB],N2[NBINS_DB],N3[NBINS_DB];
  double N1e[NBINS_DB],N2e[NBINS_DB],N3e[NBINS_DB];
  double sc[NBINS_DB][6][2],sys,ind[NBINS_DB][2],indep,d3[NBINS_DB],d2[NBINS_DB];
  /// DayaBay -----------
    for(i=0;i<NBINS_DB;i++){
      for(j=0;j<8;j++){
        J21[j]=0.0;J31[j]=0.0;J32[j]=0.0;
        J41[j]=0.0;J42[j]=0.0;J43[j]=0.0;
	J21e[j]=0.0;J31e[j]=0.0;J32e[j]=0.0;
        J41e[j]=0.0;J42e[j]=0.0;J43e[j]=0.0;
        for(k=0;k<6;k++){
          mx=p.dmq[1]*L_DB[j][k];
	      J21[j] = J21[j] + int_DB_interpol(mx,i,mL,J_DB)/pow(L_DB[j][k],2);
	      J21e[j] = J21e[j] + int_DB_interpol(mx,i,mL,Je_DB)/pow(L_DB[j][k],2);
	      mx=p.dmq[2]*L_DB[j][k];
	      J31[j] = J31[j] + int_DB_interpol(mx,i,mL,J_DB)/pow(L_DB[j][k],2);
	      J31e[j] = J31e[j] + int_DB_interpol(mx,i,mL,Je_DB)/pow(L_DB[j][k],2);
	      mx=(p.dmq[2]-p.dmq[1])*L_DB[j][k];
	      J32[j] = J32[j] + int_DB_interpol(mx,i,mL,J_DB)/pow(L_DB[j][k],2);
	      J32e[j] = J32e[j] + int_DB_interpol(mx,i,mL,Je_DB)/pow(L_DB[j][k],2);
	      mx=p.dmq[3]*L_DB[j][k];
	      J41[j] = J41[j] + int_DB_interpol(mx,i,mL,J_DB)/pow(L_DB[j][k],2);
	      J41e[j] = J41e[j] + int_DB_interpol(mx,i,mL,Je_DB)/pow(L_DB[j][k],2);
	      mx=fabs(p.dmq[3]-p.dmq[1])*L_DB[j][k];
	      J42[j] = J42[j] + int_DB_interpol(mx,i,mL,J_DB)/pow(L_DB[j][k],2);
	      J42e[j] = J42e[j] + int_DB_interpol(mx,i,mL,Je_DB)/pow(L_DB[j][k],2);
	      mx=fabs(p.dmq[3]-p.dmq[2])*L_DB[j][k];
	      J43[j] = J43[j] + int_DB_interpol(mx,i,mL,J_DB)/pow(L_DB[j][k],2);
	      J43e[j] = J43e[j] + int_DB_interpol(mx,i,mL,Je_DB)/pow(L_DB[j][k],2);
        }
      }
 
	  N1[i]=0.0;
	  N1e[i]=0.0;
	  for(j=0;j<2;j++){	// EH1
	    N1[i] = N1[i] + norm[0]*ef[j]*(I_DB[j][i] - 4.0*(c14*c14*c12*c13*s12*c13*J21[j] + c14*c14*c12*c13*s13*J31[j]
							                            +c14*s14*c12*c13*J41[j] + c14*c14*s12*c13*s13*J32[j]
							                            +c14*s14*s12*c13*J42[j] + c14*s13*s14*J43[j]));
	    N1e[i] = N1e[i] + norm[0]*ef[j]*(Ie_DB[j][i] - 4.0*(c14*c14*c12*c13*s12*c13*J21e[j] + c14*c14*c12*c13*s13*J31e[j]
							                            +c14*s14*c12*c13*J41e[j] + c14*c14*s12*c13*s13*J32e[j]
							                            +c14*s14*s12*c13*J42e[j] + c14*s13*s14*J43e[j]));
      }
	  N2[i]=0.0;
	  N2e[i]=0.0;
	  for(j=2;j<4;j++){	// EH2
	    N2[i] = N2[i] + norm[1]*ef[j]*(I_DB[j][i] - 4.0*(c14*c14*c12*c13*s12*c13*J21[j] + c14*c14*c12*c13*s13*J31[j]
							                            +c14*s14*c12*c13*J41[j] + c14*c14*s12*c13*s13*J32[j]
							                            +c14*s14*s12*c13*J42[j] + c14*s13*s14*J43[j]));
	    N2e[i] = N2e[i] + norm[1]*ef[j]*(Ie_DB[j][i] - 4.0*(c14*c14*c12*c13*s12*c13*J21e[j] + c14*c14*c12*c13*s13*J31e[j]
							                            +c14*s14*c12*c13*J41e[j] + c14*c14*s12*c13*s13*J32e[j]
							                            +c14*s14*s12*c13*J42e[j] + c14*s13*s14*J43e[j]));
      }
	  N3[i]=0.0;
	  N3e[i]=0.0;
	  for(j=4;j<8;j++){	// EH3
	    N3[i] = N3[i] + norm[2]*ef[j]*(I_DB[j][i] - 4.0*(c14*c14*c12*c13*s12*c13*J21[j] + c14*c14*c12*c13*s13*J31[j]
							                            +c14*s14*c12*c13*J41[j] + c14*c14*s12*c13*s13*J32[j]
							                            +c14*s14*s12*c13*J42[j] + c14*s13*s14*J43[j]));
	    N3e[i] = N3e[i] + norm[2]*ef[j]*(Ie_DB[j][i] - 4.0*(c14*c14*c12*c13*s12*c13*J21e[j] + c14*c14*c12*c13*s13*J31e[j]
							                            +c14*s14*c12*c13*J41e[j] + c14*c14*s12*c13*s13*J32e[j]
							                            +c14*s14*s12*c13*J42e[j] + c14*s13*s14*J43e[j]));
	  }
	  
	  ///----------------- Systematics errors coefficients --------------
          d2[i] = O2[i]/pow(O1[i]-B1[i],2)+O1[i]*pow(O2[i]-B2[i],2)/pow(O1[i]-B1[i],4);
          d3[i] = O3[i]/pow(O1[i]-B1[i],2)+O1[i]*pow(O3[i]-B3[i],2)/pow(O1[i]-B1[i],4);
  
          sc[i][0][0] = N2[i]/N1[i];						//
          sc[i][0][1] = N3[i]/N1[i];						//
          sc[i][1][0] = -N2[i]/N1[i];					//
          sc[i][1][1] = 0.0;								//
          sc[i][2][0] = 0.0;								//
          sc[i][2][1] = -N3[i]/N1[i];					//
          sc[i][3][0] = N2[i]*N1e[i]/pow(N1[i],2);		//
          sc[i][3][1] = N3[i]*N1e[i]/pow(N1[i],2);		//
          sc[i][4][0] = -N2e[i]/N1[i];					//
          sc[i][4][1] = 0.0;								//
          sc[i][5][0] = 0.0;								//
          sc[i][5][1] = -N3e[i]/N1[i];					// energy scale and efficiency errors

          ind[i][0] = (O2[i]-B2[i])/(O1[i]-B1[i])-N2[i]/N1[i]-N2[i]*N1e[i]/pow(N1[i],2)+N2e[i]/N1[i];
          ind[i][1] = (O3[i]-B3[i])/(O1[i]-B1[i])-N3[i]/N1[i]-N3[i]*N1e[i]/pow(N1[i],2)+N3e[i]/N1[i];
      ///-----------------------------------------------------------------
        }
        ///---------- Solving system of linear equations -------- 
      ///------------- For linear system of equations ------------------------
      gsl_matrix * A = gsl_matrix_alloc(nx,nx);
      gsl_vector * BQ = gsl_vector_alloc(nx);		
      gsl_vector * x = gsl_vector_alloc(nx);		
      gsl_permutation * per = gsl_permutation_alloc(nx); 
      ///---------------------------------------------------------------------  
      for(k=0;k<nx;k++){
        for(l=0;l<nx;l++){
          sys=0.0;
		  for(i=0;i<NBINS_DB;i++){
            sys = sys + sc[i][k][0]*sc[i][l][0]/d2[i] + sc[i][k][1]*sc[i][l][1]/d3[i];
	      }
	      sys = sys + cov[k][l];
	      gsl_matrix_set(A,k,l,sys);
        }
        indep=0.0;
        for(i=0;i<NBINS_DB;i++){
          indep = indep + ind[i][0]*sc[i][k][0]/d2[i] + ind[i][1]*sc[i][k][1]/d3[i];
	    }
	    for(j=0;j<nx;j++){
	      indep = indep - cov[j][k];
	    }
        gsl_vector_set(BQ,k,-indep);
      }
     
      gsl_linalg_LU_decomp(A, per, &signum);
      gsl_linalg_LU_solve(A, per, BQ, x);
      ///------------------------------------------------------  
  // Chi2 calculation ----------------------------------------------  
    double chi2;
    chi2=0.0;   
      for(k=0;k<nx;k++){
        for(l=0;l<nx;l++){
          chi2 = chi2 + (gsl_vector_get(x,k)-1.0)*cov[k][l]*(gsl_vector_get(x,l)-1.0);
        }
      }
      for(i=0;i<NBINS_DB;i++){
        chi2 = chi2 + pow( ind[i][0]
                         + gsl_vector_get(x,0)*sc[i][0][0] + gsl_vector_get(x,1)*sc[i][1][0] 
                         + gsl_vector_get(x,3)*sc[i][3][0] + gsl_vector_get(x,4)*sc[i][4][0],2)/d2[i]
                         + pow( ind[i][1]
                         + gsl_vector_get(x,0)*sc[i][0][1] + gsl_vector_get(x,2)*sc[i][2][1] 
                         + gsl_vector_get(x,3)*sc[i][3][1] + gsl_vector_get(x,5)*sc[i][5][1],2)/d3[i];
      }
  ///----------------------------------------------------------------
  ///--------------------------------------------------------------------- 
      gsl_permutation_free(per);
      gsl_vector_free(x);
      gsl_vector_free(BQ);
      gsl_matrix_free(A);  
  ///--------------------------------------------------------------------- 
  
  return chi2;
}
///---------------------------------------------------------------------


// DayaBay integrals interpolation
double DB_class::int_DB_interpol(double a,int k,double x[],double y[][2000]){
  
  int i;
  double m,n;
  
  if(a<=0.002) return 0.0;			// 10^(-1.7)
  if(a>=19952.62) return II_DB_avout[k]/2.0;	// 10^(4.3)

  i=int((log10(a)+1.7)/0.003);
  m=(y[k][i]-y[k][i+1])/(x[i]-x[i+1]);
  n=y[k][i]-(y[k][i]-y[k][i+1])/(x[i]-x[i+1])*x[i];
  
  return a*m+n;
}

#undef NBINS_DB
