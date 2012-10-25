#include "definitions.h"

namespace ns_reactor
{

extern Fit  fit;
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
   double c1[NBIN_KAML][DM2_NUM_21][KAML_NREACT][NISO];
   double c0_cor[NBIN_KAML];
   double c1_cor[NBIN_KAML][DM2_NUM_21];
} gl_co;

int old_new_kaml;

// the background and normalization factor used by calc_noe.cc
double bg_kaml[NBIN_KAML], norm_kaml[NBIN_KAML], geo_nu[NBIN_KAML], bg_acc[NBIN_KAML];

/************* static variables: *************************************/

static int cur_bin, is, r_gl;
static double dmq21_glob;

/* reactor data obtained from Michele (email on 12/04/2012)
 * Each line corresponds to a reactor complex, 
 * first col: eactor fraction, second col: distance (km). 
 * The reactor fraction tells which fraction of the total neutrino events comes from the
 * given reactor, and already incorporates the 1/dist^2 suppression term
 */
const double reactor_data[KAML_NREACT][2] = 
{  
    { 0.10810,  87.76 },
    { 0.08340, 138.47 },
    { 0.07551, 145.79 },
    { 0.21729, 159.93 },
    { 0.17332, 178.93 },
    { 0.11391, 191.48 },
    { 0.06945, 214.17 },
    { 0.01423, 295.37 },
    { 0.03354, 345.44 },
    { 0.03554, 349.46 },
    { 0.00892, 401.16 },
    { 0.01023, 430.55 },
    { 0.00793, 560.78 },
    { 0.00164, 635.87 },
    { 0.00675, 709.93 },
    { 0.01349, 711.31 },
    { 0.00683, 733.83 },
    { 0.00750, 754.47 },
    { 0.00223, 782.99 },
    { 0.00330, 830.33 },
    { 0.00691, 988.35 }
};

static double power[KAML_NREACT], dist[KAML_NREACT];

/*************** functions: *******************************************/

double flux_kaml(double eNe);
double crossSect(double eKin);
double gauss_int_kaml(double eNu);

double n0(double eNu);
double nosc1(double eNu);
double n0_cor(double eNu);
double nosc_cor(double eNu);

/******************************************************************/
/*                function init                                   */
/******************************************************************/

#define NSIGM 4.

void init_kaml(void)  
{
   // data file from fig. 1 of 0801.4589  (Data_KamL/data.out)
   // data file from fig. 1 of 1009.4771  (Data_KamL/data2010.out)
   FILE *fp = fopen(REACTOR_PATH"Data_KamL/data2010.out","r");
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
      const int ii = fit.first_bin[KAML] + i;
      fit.Data[ii] = fit.S_data[ii][ii] = Data_KamL[i];
   }
   //fprintf(stderr, "no osc events: %e, background cosm: %e, accidentals: %e\n", sum, bgtot, acc);
   //exit(0);  

   // init powers
   for(int i = 0; i < KAML_NREACT; i++)
   {
     dist[i]  = reactor_data[i][1] * 1.003; // correction factor from Michele email 12.4.2012
     power[i] = reactor_data[i][0] * dist[i] * dist[i];
   }
  
   // init pulls in class fit
   int p = fit.first_pull[KAML];

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


   // calc integrals
   old_new_kaml = old_new_main;
   double low, up;
  
   for(cur_bin = 0; cur_bin < NBIN_KAML; cur_bin++){

      low = limlow(cur_bin) + 0.8 - NSIGM * E_RES_KL * sqrt(limlow(cur_bin));
      if(low < EnuMIN) low = EnuMIN;
      up = limup(cur_bin) + 0.8 + NSIGM * E_RES_KL * sqrt(limup(cur_bin));
      if(up > EnuMAX_KL) up = EnuMAX_KL;

      gl_co.c0_cor[cur_bin] = qromb1(n0_cor, low, up, 1.e-5);

      Parameters p;
      for(int i = 0; i < DM2_NUM_21; i++)
      {		 
	p.idmq_21 = i; 
	dmq21_glob = p.dmq21();
	gl_co.c1_cor[cur_bin][i] = (dmq21_glob > 0. ? qromb1(nosc_cor, low, up, 1.e-5) : 0.);
      }

      for(r_gl = 0; r_gl < KAML_NREACT; r_gl++){
	for(is = 0; is < NISO; is++){

	  gl_co.c0[cur_bin][r_gl][is] = qromb1(n0, low, up, 1e-5);

	  for(int i = 0; i < DM2_NUM_21; i++){
		 
	    p.idmq_21 = i; 
	    dmq21_glob = p.dmq21();

	    if(dmq21_glob > 0.)
	      gl_co.c1[cur_bin][i][r_gl][is] = qromb1(nosc1, low, up, 1e-4);
	    else  
	      gl_co.c1[cur_bin][i][r_gl][is] = 0.;	       
	  }
	}
      }
   }

   if(old_new_main == NEW)
     old_new_kaml = OLD;
   double prev = 0.;
       
   // scale to no_osc prediction 
   for(int i = 0; i < NBIN_KAML; i++){
      
      double nn = 0.;
      for(int iso = 0; iso < NISO; iso++){
	    for(int r = 0; r < KAML_NREACT; r++){
	       nn += gl_co.c0[i][r][iso];
	    }
      }
      norm_kaml[i] = no_osc[i] / nn * 0.993; // correction factor from Michele email 12.04.2012
     
      // correction of scaling factor in next-to-last bin
      if(i==NBIN_KAML-2) norm_kaml[i] = prev;
      prev = norm_kaml[i];

      //printf("%e  ", norm_kaml[i]);
          
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
	//printf("%e\n", nn / n_old);
      }     
   }
   return;
}


/********************* structure definition ***************************/

struct Noe_Coeff {
   double nbin[NBIN_KAML];                      // N_i: noe in bin i
   double nr[NBIN_KAML][KAML_NREACT];           // Phi_r * d N_i / d Phi_r
   double nf[NBIN_KAML][KAML_NREACT][NISO];     // f_is,r  * d N_i / d f_is,r
   double nc[NBIN_KAML][NISO];                  // neutrino flux of isotope
   double ncor[NBIN_KAML];                     // correlated flux error
};

/**************************************************************************/
/*     calc number of events and derivatives                              */
/**************************************************************************/

void calc_noe_coeff(struct Noe_Coeff *noe_c, Param_5nu &prm)
{
  // for interpolation
#ifdef LOG_SOL
  const double x = log10(prm.dmq[1]);
#else
  const double x = prm.dmq[1];
#endif
  const int idmq = int( (x - SOL_MIN) / D_SOL);
  const double dx = (x - idmq * D_SOL - SOL_MIN) / D_SOL;

  if(x <  SOL_MIN) error("kamland: calc_noe_coeff::Dmq too small");
  if(x >= SOL_MAX) error("kamland: calc_noe_coeff::Dmq too large");
  if(idmq < 0 || idmq >= DM2_NUM_21-1) error("kamland: calc_noe_coeff: index out of range");
  if(dx < 0. || dx > 1.) error("kamland: calc_noe_coeff: x out of range");

  // coeffs for the probability
  const double sq_eff = 4. * norm(prm.Ue[0] * prm.Ue[1]);
  double cons = 1;
  for(int i = 2; i < N_NU; i++) 
    for(int j = 0; j < i; j++)
      cons -= 2. * norm(prm.Ue[i] * prm.Ue[j]);

  for(int i = 0; i < NBIN_KAML; i++){

      double w[KAML_NREACT][NISO];
      for(int is = 0; is < NISO; is++){
	 for(int r = 0; r < KAML_NREACT; r++){
 	    w[r][is] = gl_co.c0[i][r][is] * cons; 
            w[r][is] -= sq_eff * (gl_co.c1[i][idmq][r][is] + 
				  (gl_co.c1[i][idmq+1][r][is] - gl_co.c1[i][idmq][r][is]) * dx);
	    w[r][is] *= norm_kaml[i];
	    (*noe_c).nf[i][r][is] = w[r][is];
	 }
      }

      // correlated flux error 
      (*noe_c).ncor[i] =  gl_co.c0_cor[i] * cons;
      (*noe_c).ncor[i] -= sq_eff * (gl_co.c1_cor[i][idmq] +
				    (gl_co.c1_cor[i][idmq+1] - gl_co.c1_cor[i][idmq]) * dx);
      (*noe_c).ncor[i] *= norm_kaml[i];


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

void set_table_kaml(Param_5nu &prm, double cff[NBIN_CHISQ][NPULLS+1])
{
   struct Noe_Coeff noe_c;
   calc_noe_coeff(&noe_c, prm);

   for(int i = 0; i < NBIN_KAML; i++){
      const int ii = fit.first_bin[KAML] + i;

      int p = fit.first_pull[KAML]; 
      cff[ii][NPULLS] = cff[ii][p] = cff[ii][FLUX_NORM] = noe_c.nbin[i];    // norm
      p++;
     
      for(int r = 0; r < KAML_NREACT; r++, p++)
	cff[ii][p] = noe_c.nr[i][r];                  // power of reactors

      for(int r = 0; r < KAML_NREACT; r++)
	 for(int is = 0; is < NISO; is++, p++)
	    cff[ii][p] = noe_c.nf[i][r][is];            // isotope fractions

      cff[ii][p] = bg_kaml[i];                               // background (only cosmog bg)
                                                       // accidentals not included (well known)
      p += 2;                                          // geo nus
      cff[ii][p] = geo_nu[i];

      // flux uncertainties (only norm from uncorr uncertainties)
      cff[ii][PULL_U235_0] = noe_c.nc[i][U235];
      cff[ii][PULL_P239_0] = noe_c.nc[i][P239];
      cff[ii][PULL_P241_0] = noe_c.nc[i][P241];
      cff[ii][PULL_U238]   = noe_c.nc[i][U238];
      cff[ii][FLUX_COR]    = noe_c.ncor[i];
   }

   // energy scale
   const int p_e = fit.first_pull[KAML] + PULL_E_SCALE;  
   for(int i = 0; i < NBIN_KAML; i++){
      const int ii = fit.first_bin[KAML] + i;

      cff[ii][p_e] = - cff[ii][NPULLS] * limlow(i)/(limup(i) - limlow(i));
      if(i < NBIN_KAML-1)
        cff[ii][p_e] += cff[ii+1][NPULLS] * limlow(i+1)/(limup(i+1) - limlow(i+1));
   }
   return;
}


/**********************************************************************************/





// P = 1 - sq2t * sin^2 D 

double n0(double eNu)
{
   const double EposKin=eNu-ME-DELTA;

   double w = power[r_gl] * flux_kaml(eNu) / norm(dist[r_gl]);
   return w * crossSect(EposKin) * gauss_int_kaml(eNu);
}

double nosc1(double eNu)
{
   const double EposKin=eNu-ME-DELTA;

   double w = norm( sin(1.27e3 * dmq21_glob * dist[r_gl] / eNu) );
   w *= power[r_gl] * flux_kaml(eNu) / norm(dist[r_gl]);
   return w * crossSect(EposKin) * gauss_int_kaml(eNu);
}

double n0_cor(double eNu)
{
  double w = 0.;
  for(int r = 0; r < KAML_NREACT; r++)
    w += power[r] / norm(dist[r]);

  is = U235;
  double v = read(REACTOR_PATH"Dat/Patrick-U235-err_cor.dat", eNu) * flux_kaml(eNu);
  is = P239;
  v += read(REACTOR_PATH"Dat/Patrick-Pu239-err_cor.dat", eNu) * flux_kaml(eNu);
  is = P241;
  v += read(REACTOR_PATH"Dat/Patrick-Pu241-err_cor.dat", eNu) * flux_kaml(eNu);

  const double EposKin=eNu-ME-DELTA;
  return v * w * crossSect(EposKin)*gauss_int_kaml(eNu);
}

double nosc_cor(double eNu)
{
  double w = 0.;
  for(int r = 0; r < KAML_NREACT; r++)
    w += power[r] / norm(dist[r]) * norm( sin(1.27e3 * dmq21_glob * dist[r] / eNu) );

  is = U235;
  double v = read(REACTOR_PATH"Dat/Patrick-U235-err_cor.dat", eNu) * flux_kaml(eNu);
  is = P239;
  v += read(REACTOR_PATH"Dat/Patrick-Pu239-err_cor.dat", eNu) * flux_kaml(eNu);
  is = P241;
  v += read(REACTOR_PATH"Dat/Patrick-Pu241-err_cor.dat", eNu) * flux_kaml(eNu);

  const double EposKin=eNu-ME-DELTA;
  return v * w * crossSect(EposKin)*gauss_int_kaml(eNu);
}


double flux_kaml(double eNu)
{
   const double isoenergy[NISO]={201.7, 210.0, 205.0, 212.4};
   const double isofract[NISO]={0.568, 0.078, 0.297, 0.057};
  
   return 1.e8 * isofract[is]/isoenergy[is] * global_flux(is, eNu, old_new_kaml);
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

}
