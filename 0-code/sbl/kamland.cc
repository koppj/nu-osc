#include "def-reactors.h"

extern Fit  fit;
extern Flux flux[NISO][2];  // old and new fluxes for each isotope
extern int old_new_main;

#define EnuMAX_KL 10.

#define E_RES_KL 0.065

#define BG_ERROR 0.11  // error on cosmogenic background
                       // estimated from tab. II of 1009.4771
#define EN_SCALE 0.011 
#define GEO_NU   10.0

#define D_E 0.425

inline double limlow(int bin){
   return 0.9 + bin * D_E;
}

inline double limup(int bin){
   return limlow(bin) + D_E;
}


/************* global variables: *************************************/

struct Coeff_kaml {
   double c0[NBIN_KAML][KAML_NREACT][NISO];
   double c1[NBIN_KAML][KAML_NREACT][NISO];
   double c2[NBIN_KAML][KAML_NREACT][NISO];
} gl_co;

int old_new_kaml;

// the background and normalization factor used by calc_noe.cc
double bg_kaml[NBIN_KAML], norm_kaml[NBIN_KAML], geo_nu[NBIN_KAML], bg_acc[NBIN_KAML];

/************* static variables: *************************************/

static int cur_bin, is, r_gl;

/* reactor data obtained from Michele (email on 15/10/2010 - fwd 29/10/2010)
 * Each line corresponds to a reactor complex, 
 * the first data is essentially the power divided by the square of the distance 
 * (which should represent the number of neutrinos arriving at Kamioka) 
 * while the second is the distance.
 */
const double reactor_data[KAML_NREACT][2] = 
{
   { 37.84,  87.9 },
   { 30.15, 138.4 },
   { 27.72, 145.8 },
   { 80.81, 159.3 },
   { 64.74, 179.0 },
   { 41.13, 191.5 },
   { 25.24, 214.6 },
   {  5.18, 295.1 },
   { 12.18, 345.3 },
   { 12.79, 349.2 },
   {  3.26, 400.8 },
   {  3.74, 430.4 },
   {  2.91, 561.0 },
   {  0.58, 642.2 },
   {  2.24, 709.2 },
   {  4.82, 710.5 },
   {  2.36, 733.4 },
   {  2.79, 754.5 },
   {  0.85, 784.5 },
   {  1.20, 830.9 },
   {  2.45, 987.3 }
};

static double power[KAML_NREACT], dist[KAML_NREACT];

/*************** functions: *******************************************/

double flux_kaml(double eNe);
double crossSect(double eKin);
double gauss_int_kaml(double eNu);

double n0(double eNu);
double nosc1(double eNu);
double nosc2(double eNu);

/******************************************************************/
/*                function init                                   */
/******************************************************************/

#define NSIGM 4.

void init_kaml(void)  
{
   // data file from fig. 1 of 0801.4589  (Data_KamL/data.out)
   // data file from fig. 1 of 1009.4771  (Data_KamL/data2010.out)
   FILE *fp = fopen("Data_KamL/data2010.out","r");
   if(!fp) error("cannot open Data_KamL/data2010.out");
   double no_osc[NBIN_KAML];
   
   double sum = 0., bgtot = 0., acc = 0.;
   double Data_KamL[NBIN_KAML];
  
   for(int i = 0; i < NBIN_KAML; i++){
      double geo;
      if(fscanf(fp, "%lf %lf %lf %lf %lf\n", &Data_KamL[i], &bg_kaml[i], &no_osc[i], &geo, &bg_acc[i]) != 5)
	error("error in reading Data_KamL/data.out");
      
      // abs event numbers
      //Data_KamL[i] *= D_E;
      //bg_kaml[i] *= D_E;
      //no_osc[i] *= D_E;
      geo_nu[i] = geo - bg_kaml[i];
      bg_kaml[i] -= bg_acc[i];
     
      sum += no_osc[i];
      bgtot += bg_kaml[i];
      acc += bg_acc[i];
     
      // set the data and stat. error
      const int ii = BIN_KAML_0+i;
      fit.Data[ii] = fit.S_data[ii][ii] = Data_KamL[i];
   }
   //fprintf(stderr, "no osc events: %e, background cosm: %e, accidentals: %e\n", sum, bgtot, acc);
   //exit(0);  

   // init powers
   for(int i = 0; i < KAML_NREACT; i++)
   {
     dist[i]  = reactor_data[i][1];
     power[i] = reactor_data[i][0] * dist[i] * dist[i];
   }
  
   // init pulls in class fit
   int p = PULL_KAML_0;
#ifndef NO_KAML     
   fit.S_pull[p][p] = norm(0.018);           // detector normalization
   p++; 
  
   for(int r = 0; r < KAML_NREACT; r++, p++)
     fit.S_pull[p][p] = norm(0.021);            // power of reactors

   for(int r = 0; r < KAML_NREACT; r++)
      for(int i = 0; i < NISO; i++, p++)
	 fit.S_pull[p][p] = norm(0.01);         // isotope fractions
  
   fit.S_pull[p][p] = norm(BG_ERROR);
   p++;
   fit.S_pull[p][p] = norm(EN_SCALE);   
   p++;
   fit.S_pull[p][p] = norm(GEO_NU);
#endif
   if(p+1 != NPULLS) error("wrong number of pulls!");

   // calc integrals
   old_new_kaml = old_new_main;
   double low, up;
  
   for(cur_bin = 0; cur_bin < NBIN_KAML; cur_bin++){

      low = limlow(cur_bin) + 0.8 - NSIGM * E_RES_KL * sqrt(limlow(cur_bin));
      if(low < EnuMIN) low = EnuMIN;
      up = limup(cur_bin) + 0.8 + NSIGM * E_RES_KL * sqrt(limup(cur_bin));
      if(up > EnuMAX_KL) up = EnuMAX_KL;

      for(r_gl = 0; r_gl < KAML_NREACT; r_gl++){
	for(is = 0; is < NISO; is++){

	  gl_co.c0[cur_bin][r_gl][is] = qromb1(n0, low, up, 1e-5);

	  gl_co.c1[cur_bin][r_gl][is] = qromb1(nosc1, low, up, 1e-4);
#ifdef NO_MAT // without matter effect
	  gl_co.c2[cur_bin][r_gl][is] = 0.;
#else		  
	  gl_co.c2[cur_bin][r_gl][is] = qromb1(nosc2, low, up, 1e-4);
#endif		  
	}
      }
   }

   if(old_new_main == NEW)
     old_new_kaml = OLD;
       
   // scale to no_osc prediction 
   for(int i = 0; i < NBIN_KAML; i++){
      
      double nn = 0.;
      for(int iso = 0; iso < NISO; iso++){
	    for(int r = 0; r < KAML_NREACT; r++){
	       nn += gl_co.c0[i][r][iso];
	    }
      }
      norm_kaml[i] = no_osc[i] / nn;
     
      // correction of scaling factor in next-to-last bin
      if(i==NBIN_KAML-2) norm_kaml[i] = norm_kaml[i-1];

          
      // scaling of the predition due to the new fluxes
      if(old_new_main == NEW){
	cur_bin = i;

        low = limlow(cur_bin) + 0.8 - NSIGM * E_RES_KL * sqrt(limlow(cur_bin));
        if(low < EnuMIN) low = EnuMIN;
        up = limup(cur_bin) + 0.8 + NSIGM * E_RES_KL * sqrt(limup(cur_bin));
        if(up > EnuMAX_KL) up = EnuMAX_KL;

	double n_old = 0.;
        for(r_gl = 0; r_gl < KAML_NREACT; r_gl++){
	   for(is = 0; is < NISO; is++){
	     n_old += qromb1(n0, low, up, 1e-5);
	   }
	}
	
	norm_kaml[i] *= nn / n_old;
      }
     
      //fprintf(stderr, "norm fact = %e\n", norm_kaml[i]);         
   }
   return;
}


/********************* structure definition ***************************/

struct Noe_Coeff {
   double nbin[NBIN_KAML];                 // N_i: noe in bin i
   double nr[NBIN_KAML][KAML_NREACT];           // Phi_r * d N_i / d Phi_r
   double nf[NBIN_KAML][KAML_NREACT][NISO];     // f_is,r  * d N_i / d f_is,r
   double nc[NBIN_KAML][NISO];             // neutrino flux of isotope
};

/**************************************************************************/
/*     calc number of events and derivatives                              */
/**************************************************************************/

// P = 1 - 2 sq13 cq13  - c^4_13 * sq2t * sin^2 D
//       - c^6_13 * sq2t * c2t [ 2*A sin D * (sin D - D cos D) ]

void calc_noe_coeff(struct Noe_Coeff *noe_c, double sqt, double sq13)
{
   const double sq2t = 4. * sqt * (1. - sqt);
   const double c2t = 1. - 2.*sqt;
   const double cq13 = 1. - sq13;

   for(int i = 0; i < NBIN_KAML; i++){

      double w[KAML_NREACT][NISO];
      for(int is = 0; is < NISO; is++){
	 for(int r = 0; r < KAML_NREACT; r++){
 	    w[r][is] = gl_co.c0[i][r][is] * (1. - 2. * sq13 * cq13); 
            w[r][is] -= norm(cq13) * sq2t * gl_co.c1[i][r][is];
	    w[r][is] -= pow(cq13, 3) * sq2t * c2t * gl_co.c2[i][r][is];
	    w[r][is] *= norm_kaml[i];
	    (*noe_c).nf[i][r][is] = w[r][is];
	 }
      }

      (*noe_c).nbin[i] = 0.;
      for(int is = 0; is < NISO; is++){
	(*noe_c).nc[i][is] = 0.;
	
        for(int r = 0; r < KAML_NREACT; r++){
	  (*noe_c).nbin[i] += w[r][is];    // summed over react and iso
	  (*noe_c).nc[i][is] += w[r][is];  // summed over react
	}
      }
     
      (*noe_c).nbin[i] += bg_kaml[i] + geo_nu[i] + bg_acc[i];

      for(int r = 0; r < KAML_NREACT; r++){
	 (*noe_c).nr[i][r] = 0.;
	 for(int is = 0; is < NISO; is++){
	    (*noe_c).nr[i][r] += w[r][is];
	 }
      }
   }   
   return;
}

/**************************************************************************/
/*     function kamland                                                   */
/**************************************************************************/

#define SQ12 0.3

void set_table_kaml(const params &prm, double cff[NBIN_CHISQ][NPULLS+1])
{
   struct Noe_Coeff noe_c;
   calc_noe_coeff(&noe_c, SQ12, sqr(prm.Ue3));

   for(int i = 0; i < NBIN_KAML; i++){
      const int ii = BIN_KAML_0 + i;

      int p = PULL_KAML_0; 
      cff[ii][NPULLS] = cff[ii][p] = cff[ii][FLUX_NORM] = noe_c.nbin[i];    // norm
      p++;
     
      for(int r = 0; r < KAML_NREACT; r++, p++)
	cff[ii][p] = noe_c.nr[i][r];                  // power of reactors

      for(int r = 0; r < KAML_NREACT; r++){
	 for(int is = 0; is < NISO; is++, p++){

	    cff[ii][p] = noe_c.nf[i][r][is];            // isotope fractions
	 }
      }

      for(int is = 0; is < NISO; is++)
	cff[ii][is] = noe_c.nc[i][is];                 // neutrino flux for each iso
     
      cff[ii][p] = bg_kaml[i];                               // background (only cosmog bg)
                                                       // accidentals not included (well known)

      p += 2;                                          // geo nus
      cff[ii][p] = geo_nu[i];
   }

   // energy scale
   const int p_e = PULL_KAML_0 + PULL_E_SCALE;  
   for(int i = 0; i < NBIN_KAML; i++){
      const int ii = BIN_KAML_0 + i;

      cff[ii][p_e] = - cff[ii][NPULLS] * limlow(i)/(limup(i) - limlow(i));
      if(i < NBIN_KAML-1)
        cff[ii][p_e] += cff[ii+1][NPULLS] * limlow(i+1)/(limup(i+1) - limlow(i+1));
   }
   return;
}


/**********************************************************************************/





#define TWO_V_MATT -2.268e-7   

// P = 1 - sq2t * sin^2 D - sq2t * c2t [ 2*A sin D * (sin D - D cos D) ]

double n0(double eNu)
{
   double EposKin=eNu-ME-DELTA;

   double w = power[r_gl] * flux_kaml(eNu) / norm(dist[r_gl]);
   return w * crossSect(EposKin) * gauss_int_kaml(eNu);
}

double nosc1(double eNu)
{
   double EposKin=eNu-ME-DELTA;

   double w = norm( sin(1.27e3 * DMQ_21 * dist[r_gl] / eNu) );
   w *= power[r_gl] * flux_kaml(eNu) / norm(dist[r_gl]);
   return w * crossSect(EposKin) * gauss_int_kaml(eNu);
}

double nosc2(double eNu)
{
   const double EposKin=eNu-ME-DELTA;
   const double A = TWO_V_MATT * eNu / DMQ_21;
   const double D = 1.27e3 * DMQ_21 * dist[r_gl] / eNu;

   double w = 2. * A * sin(D) * ( sin(D) - D * cos(D) );
   w *= power[r_gl] * flux_kaml(eNu) / norm(dist[r_gl]);
   return w * crossSect(EposKin) * gauss_int_kaml(eNu);
}


double flux_kaml(double eNu)
{
   const double isoenergy[NISO]={201.7, 210.0, 205.0, 212.4};
   const double isofract[NISO]={0.568, 0.078, 0.297, 0.057};
  
   return 1.e8 * isofract[is]/isoenergy[is] * flux[is][old_new_kaml].f(eNu);
}


/****************************************************/
/* SUBROUTINE gauss_int_kaml(T_min, T_max, T_e, integ_T) */
/****************************************************/

double gauss_int_kaml(double eNu)
{
   const double T_min = limlow(cur_bin);
   const double T_max = limup(cur_bin);
   const double T_e = eNu + ME - DELTA;

   double sigma = E_RES_KL * sqrt(T_e);

   double x_max = (T_max - T_e) / (sigma*M_SQRT2);
   double x_min = (T_min - T_e) / (sigma*M_SQRT2);

   double res = (erf(x_max) - erf(x_min)) / 2.;
   if(res < 0.) {
      fprintf(stderr, "BLOODY HELL!!\n");
      exit(1);
   }

   return res;
}

