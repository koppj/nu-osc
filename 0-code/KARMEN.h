#ifndef H_KARMEN_
#define H_KARMEN_
#include <gsl/gsl_spline.h>
#include <fstream>
#include <utility>
#include <algorithm>
//JK#include "prob_v15_scalar.h"
//#include "snu.h"    // JK

#ifdef NU_USE_NUSQUIDS // JK

using namespace regeneration;

namespace ns_sbl { // JK

class KARMEN{
public:
  /* Experimental setup */
  static constexpr int N_BINS = 9;

  const double Lmin = 15.7*units::meter;
  const double Lmax = 19.7*units::meter;
  const double bins[N_BINS][2] = {{16.0*units::MeV, 20.0*units::MeV},
				  {20.0*units::MeV, 24.0*units::MeV},
				  {24.0*units::MeV, 28.0*units::MeV},
				  {28.0*units::MeV, 32.0*units::MeV},
				  {32.0*units::MeV, 36.0*units::MeV},
				  {36.0*units::MeV, 40.0*units::MeV},
				  {40.0*units::MeV, 44.0*units::MeV},
				  {44.0*units::MeV, 48.0*units::MeV},
				  {48.0*units::MeV, 52.0*units::MeV}};

  const double bkg[N_BINS] = {3.370370, 3.681481, 3.370370, 2.333333, 1.140741, 0.777778, 0.570370, 0.362963, 0.155556};
  const double data[N_BINS] = {3.0, 4.0, 1.0, 3.0, 1.0, 0.0, 0.0, 0.0};

private:
  gsl_interp_accel *acc_finalFlux;
  gsl_spline *spl_finalFlux_L01, *spl_finalFlux_L02, *spl_finalFlux_L11, *spl_finalFlux_L12;
  
  static constexpr double m_mu = 105.66*units::MeV;
  static constexpr double MN = 939.566*units::MeV;
  static constexpr double MP = 938.272*units::MeV;
  static constexpr double DELTA = (MN-MP);
  static constexpr double ME = 0.511*units::MeV;
  static constexpr double karmen_lowpass_width = 5e5; // [eV] JK
  regProb prb; // Probability engine
  
  unsigned int N_PTS_INTERP = 50; // Number of points in which the final spectrum will be interpolated

  /* True energy limits where the energy resolution Gaussian falls by +/- 3 sigma */

  const double Elims[N_BINS][2] = {{15.862*units::MeV, 20.1543*units::MeV},
				   {19.8457*units::MeV, 24.169*units::MeV},
				   {23.831*units::MeV, 28.1826*units::MeV},
				   {27.8174*units::MeV, 32.1952*units::MeV},
				   {31.8048*units::MeV, 36.207*units::MeV},
				   {35.793*units::MeV, 40.2182*units::MeV},
				   {39.7818*units::MeV, 44.2288*units::MeV},
				   {43.7712*units::MeV, 48.239*units::MeV},
				   {47.761*units::MeV, 52.2488*units::MeV}};
  
  /* Flux functions [THE ARGUMENT IS IN eV] */
  
  inline double getFlux_e(double E){
    if (E > m_mu/2)
      return 0;
    return 96/SQR(SQR(m_mu)) * (m_mu*SQR(E) - 2*E*SQR(E));
  }

  double getFlux_ebar(double E){
    return 0;
  }

  double getFlux_mu(double E){
    return 0;
  }

  double getFlux_mubar(double E){
    if (E > m_mu/2)
      return 0;
    return 16/SQR(SQR(m_mu)) * (3*m_mu*SQR(E) - 4*E*SQR(E));
  }

  // converts visible energy to nu energy
  inline double Evis2Enu(double Evis){ 
    return Evis - ME + DELTA;
  }
  // converts nu energy to visible energy
  inline double Enu2Evis(double Enu){
    return Enu + ME - DELTA;
  }
  // converts nu energy to kinetic energy
  inline double Enu2Ekin(double Enu){
    return Enu2Evis(Enu) - 2.*ME;
  }

  // energy resolution
  inline double sigma(double Evis){
    return 0.0115 * sqrt(Evis/units::MeV);
  }

  double crossSect(double eKin){
    /* Returns the IBD cross section
       eKin is the KINETIC energy
     */
    const double f=1.0, g=1.26, f2=3.7;
    double e, p, psq, cs0, gamma;
   
    e = eKin + ME;
    psq = SQR(e) - SQR(ME);
   
    if(eKin <= 0. || psq <= 0)
      return 0.;
   
    p=sqrt(psq);
    cs0=(SQR(f)+3*SQR(g))*e*p;
    gamma=2.*(f+f2)*g*((2.*e+DELTA)-SQR(ME)/e);
    gamma+=(SQR(f)+SQR(g))*(DELTA+SQR(ME)/e);
    gamma+=(SQR(f)+3.*SQR(g))*e;
    gamma+= -2./3.*(SQR(f)-SQR(g))*(e+DELTA);
    return (cs0-gamma/MP*e*p)* 2913 / 3.73313e+15;; //Added normalisation factor so that for maximal mixing we predict 2913 events
  }
  
  /* Utilities for numerical integration */
  const double x32[16] = {0.0483076656877383162348126,0.1444719615827964934851864,0.2392873622521370745446032,0.3318686022821276497799168,0.4213512761306353453641194,0.5068999089322293900237475,0.5877157572407623290407455,0.6630442669302152009751152,0.7321821187402896803874267,0.7944837959679424069630973,0.8493676137325699701336930,0.8963211557660521239653072,0.9349060759377396891709191,0.9647622555875064307738119,0.9856115115452683354001750,0.9972638618494815635449811};
  const double w32[16] = {0.0965400885147278005667648,0.0956387200792748594190820,0.0938443990808045656391802,0.0911738786957638847128686,0.0876520930044038111427715,0.0833119242269467552221991,0.0781938957870703064717409,0.0723457941088485062253994,0.0658222227763618468376501,0.0586840934785355471452836,0.0509980592623761761961632,0.0428358980222266806568786,0.0342738629130214331026877,0.0253920653092620594557526,0.0162743947309056706051706,0.0070186100094700966004071};
  const double x4[2] = {0.3399810435848562648026658,0.8611363115940525752239465};
  const double w4[2] = {0.6521451548625461426269361,0.3478548451374538573730639};

  /* Returns the energy resolution factor */
  inline double eres(int n_bin, double Enu){
    return 0.5 * (erf((bins[n_bin][1] - Enu2Evis(Enu)) / (M_SQRT2 * sigma(Enu2Evis(Enu)))) -
		  erf((bins[n_bin][0] - Enu2Evis(Enu)) / (M_SQRT2 * sigma(Enu2Evis(Enu)))));
  }

public:
  /* Constructor */
  KARMEN(): prb(OSC_DECAY_MAJORANA){ //JK
    /* Load the interpolator */
    acc_finalFlux = gsl_interp_accel_alloc();

    spl_finalFlux_L01 = gsl_spline_alloc(gsl_interp_linear, N_PTS_INTERP),
    spl_finalFlux_L02 = gsl_spline_alloc(gsl_interp_linear, N_PTS_INTERP),
    spl_finalFlux_L11 = gsl_spline_alloc(gsl_interp_linear, N_PTS_INTERP),
    spl_finalFlux_L12 = gsl_spline_alloc(gsl_interp_linear, N_PTS_INTERP);    
    
    /* Load the initial state */
    static constexpr int N_PTS = 5000;
    marray<double,1> E_range = regeneration::linspace(1.5*units::MeV, 55*units::MeV, N_PTS); //JK
    marray<double, 3> inistate{N_PTS, 2, 4};

    for(uint i_E=0; i_E<N_PTS; ++i_E){
      inistate[i_E][0][0] = getFlux_e(E_range[i_E]);
      inistate[i_E][0][1] = getFlux_mu(E_range[i_E]);
      inistate[i_E][0][2] = inistate[i_E][0][3] = 0.;

      inistate[i_E][1][0] = getFlux_ebar(E_range[i_E]);
      inistate[i_E][1][1] = getFlux_mubar(E_range[i_E]);
      inistate[i_E][1][2] = inistate[i_E][0][3] = 0.;
    }

    prb.setIniState(E_range, inistate);
  }

  std::array<double, N_BINS> getSpectrum(const Param &prm, bool includeBkg = true){
    prb.setParam(prm);

    /* Precompute the spectrum for the 4 different values of L that are used in the integral over L*/
    double A_L = (Lmax-Lmin)/2.;
    double B_L = (Lmax+Lmin)/2.;
    double L_01 = A_L*x4[0] + B_L, L_02 = -A_L*x4[0] + B_L, L_11 = A_L*x4[1] + B_L, L_12 = -A_L*x4[1] + B_L;

    double energies[N_PTS_INTERP];
    double fluxes_01[N_PTS_INTERP], fluxes_02[N_PTS_INTERP], fluxes_11[N_PTS_INTERP], fluxes_12[N_PTS_INTERP]; //Electron antineutrino fluxes
    
    for(unsigned int i=0; i<N_PTS_INTERP; ++i){
      energies[i] = Evis2Enu( Elims[0][0] + (Elims[N_BINS-1][1] - Elims[0][0])*i/((double) N_PTS_INTERP - 1) );
      fluxes_01[i] = prb.getFinalFlux(energies[i], L_01, karmen_lowpass_width)[1][0];
      fluxes_02[i] = prb.getFinalFlux(energies[i], L_02, karmen_lowpass_width)[1][0];
      fluxes_11[i] = prb.getFinalFlux(energies[i], L_11, karmen_lowpass_width)[1][0];
      fluxes_12[i] = prb.getFinalFlux(energies[i], L_12, karmen_lowpass_width)[1][0];
    }
      
    gsl_spline_init(spl_finalFlux_L01, energies, fluxes_01, N_PTS_INTERP);
    gsl_spline_init(spl_finalFlux_L02, energies, fluxes_02, N_PTS_INTERP);
    gsl_spline_init(spl_finalFlux_L11, energies, fluxes_11, N_PTS_INTERP);
    gsl_spline_init(spl_finalFlux_L12, energies, fluxes_12, N_PTS_INTERP);

    std::array<double, N_BINS> result;
    
    for(size_t n=0; n<N_BINS; ++n){
      double theor;
      
      /* Do the integral: convolve final flux*xsec*energy resolution, averaging over L */
      double A_E = (Elims[n][1] - Elims[n][0])/2.;
      double B_E = (Elims[n][1] + Elims[n][0])/2.;
      double A_L = (Lmax-Lmin)/2.;
//JK      double B_L = (Lmax+Lmin)/2.;
      theor = 0;

      for(int i_E=0; i_E<16; ++i_E){
	double E_1 = A_E*x32[i_E] + B_E; double Enu_1 = Evis2Enu(E_1); double Ekin_1 = Enu2Ekin(Enu_1);
	double E_2 = -A_E*x32[i_E] + B_E; double Enu_2 = Evis2Enu(E_2); double Ekin_2 = Enu2Ekin(Enu_2);
	double flx_11_0 = gsl_spline_eval(spl_finalFlux_L01, Enu_1, acc_finalFlux),
	  flx_12_0 = gsl_spline_eval(spl_finalFlux_L02, Enu_1, acc_finalFlux),
	  flx_11_1 = gsl_spline_eval(spl_finalFlux_L11, Enu_1, acc_finalFlux),
	  flx_12_1 = gsl_spline_eval(spl_finalFlux_L12, Enu_1, acc_finalFlux),
	  flx_21_1 = gsl_spline_eval(spl_finalFlux_L11, Enu_2, acc_finalFlux),
	  flx_22_1 = gsl_spline_eval(spl_finalFlux_L12, Enu_2, acc_finalFlux),
	  flx_21_0 = gsl_spline_eval(spl_finalFlux_L01, Enu_2, acc_finalFlux),
	  flx_22_0 = gsl_spline_eval(spl_finalFlux_L02, Enu_2, acc_finalFlux);
	    
	theor += A_E * A_L * w4[0] * w32[i_E] *
	  ( flx_11_0 * crossSect(Ekin_1) * eres(n, Enu_1) / SQR(L_01) +
	    flx_12_0 * crossSect(Ekin_1) * eres(n, Enu_1) / SQR(L_02) +
	    flx_21_0 * crossSect(Ekin_2) * eres(n, Enu_2) / SQR(L_01) +
	    flx_22_0 * crossSect(Ekin_2) * eres(n, Enu_2) / SQR(L_02) ) +
	  A_E * A_L * w4[1] * w32[i_E] *
	  ( flx_11_1 * crossSect(Ekin_1) * eres(n, Enu_1) / SQR(L_11) +
	    flx_12_1 * crossSect(Ekin_1) * eres(n, Enu_1) / SQR(L_12) +
	    flx_21_1 * crossSect(Ekin_2) * eres(n, Enu_2) / SQR(L_11) +
	    flx_22_1 * crossSect(Ekin_2) * eres(n, Enu_2) / SQR(L_12) );
      }
      
      theor /= (1/Lmin - 1/Lmax);

      /* Add the background */
      if(includeBkg)
	theor += bkg[n];

      result[n] = theor;
    }

    return result;    
  }
  double getChisq(const Param &prm){
    std::array<double, N_BINS> pred = getSpectrum(prm, true);
    
    constexpr double ERR_SYS = 0.065;

    double the_tot = 0, exp_tot = 0;
    for(int n=0; n<N_BINS; ++n){
      the_tot += pred[n];
      exp_tot += data[n];
    }
    
    double cA = 1. / SQR(ERR_SYS);
    double cB = cA + the_tot;
    double cC = the_tot - exp_tot;

    double delta = SQR(cB) - 4.*cA*cC;
    double a_min = (sqrt(delta) - cB) / (2.*cA);
    
    double chq = SQR(a_min / ERR_SYS);

    for(int n = 0; n<N_BINS; ++n){
      double theor = (1 + a_min)*pred[n];
      if (data[n] > 1e-7)
	chq += theor - data[n] + data[n]*log(data[n] / theor);
      else
	chq += theor - data[n];
    }
    chq *= 2;
    
    return chq;
  }

//  int print_spectrum(const Param &prm){   // JK
//    std::array<double, 11> pred = getSpectrum(prm, true);
//    for (size_t n=0; n < pred.size(); n++){
//      printf("# KARMENSPECT  %ld  %g\n", n, pred[n]);
//    }
//    return 0;
//  }

  /* Destructor */
  ~KARMEN(){
    gsl_spline_free(spl_finalFlux_L01);
    gsl_spline_free(spl_finalFlux_L02);
    gsl_spline_free(spl_finalFlux_L11);
    gsl_spline_free(spl_finalFlux_L12);
    gsl_interp_accel_free(acc_finalFlux);
  }

  /* No copy (for simplicity) */
  KARMEN(const KARMEN&) = delete;
  KARMEN &operator=(const KARMEN&) = delete;
};

} // end namespace ns_sbl (JK)

#endif

#endif

