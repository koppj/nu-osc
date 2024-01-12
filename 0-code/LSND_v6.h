#ifndef H_LSND_
#define H_LSND_
#include <gsl/gsl_spline.h>
#include <fstream>
#include <utility>
// JK #include <array>
#include <algorithm>
#include "osc_decay/osc_decay.h" // JK
#include "snu.h"    // JK

#ifdef NU_USE_NUSQUIDS // JK

#define SQR(x) ((x)*(x))

using namespace regeneration;

namespace ns_sbl { // JK

class LSND{
public:
  /* Experimental setup */
  const double Lmin = 25.85*units::meter;
  const double Lmax = 34.15*units::meter;
  const double FCTR = units::meter/units::MeV;
  const double bins[11][2] = {{0.4*FCTR, 0.5*FCTR}, {0.5*FCTR, 0.6*FCTR},
			      {0.6*FCTR, 0.7*FCTR}, {0.7*FCTR, 0.8*FCTR},
			      {0.8*FCTR, 0.9*FCTR}, {0.9*FCTR, 1.0*FCTR},
			      {1.0*FCTR, 1.1*FCTR}, {1.1*FCTR, 1.2*FCTR},
			      {1.2*FCTR, 1.3*FCTR}, {1.3*FCTR, 1.4*FCTR},
			      {1.4*FCTR, 1.5*FCTR}
  };

  const double bkg[11] = {0.32723, 0.66329, 0.62780, 1.01812, 1.30893, 1.24529, 1.18773,
			  1.23634, 0.74152, 0.36418, 0.01557};
  const double data[11] = {3.51486, 9.62374, 6.68933, 8.23408, 7.77990, 6.39282, 3.41284,
			   4.51452, -0.53821, 0.00861, 1.24661};
  const double intrinsic[11] = {0.15756, 1.11457, 1.89502, 1.75170, 1.28377, 0.77955,
				0.52717, 0.34606, 0.17877, 0.07187, 0.03788};
  const double eff[11] = {28.71048, 17.10480, 14.68284, 13.90392, 13.80024, 12.41772,
			  12.69828, 12.54444, 9.617724, 5.6433, 4.803708};
  const double no_beam[11] = {1.48514, 4.37626, 5.31067, 3.76592, 4.22010, 3.60718,
			      3.58716, 3.48548, 3.53821, 0.99139, 0.75339};
  const double lsnd_lowpass_width = 5e5; // [eV]
  
private:
  gsl_interp_accel *acc_e, *acc_ebar, *acc_mu, *acc_mubar_DAR, *acc_mubar_DIF, *acc_finalFlux;
  gsl_spline *spl_e, *spl_ebar, *spl_mu, *spl_mubar_DAR, *spl_mubar_DIF,
    *spl_finalFlux_L01, *spl_finalFlux_L02, *spl_finalFlux_L11, *spl_finalFlux_L12;
  double E_e[249], E_ebar[249], E_mu[240], E_mubar_DIF[250], E_mubar_DAR[53];
  double f_e[249], f_ebar[249], f_mu[240], f_mubar_DIF[250], f_mubar_DAR[53];

  regProb prb; // Probability engine

  unsigned int N_PTS_INTERP = 50;

  /* Flux functions [THE ARGUMENT IS IN eV] */
  
  double getFlux_e(double E){
    if(E<E_e[0] || E>E_e[248])
      return 0;
    else
      return gsl_spline_eval(spl_e, E, acc_e);
  }

  double getFlux_ebar(double E){
    if(E<E_ebar[0] || E>E_ebar[248])
      return 0;
    else
      return gsl_spline_eval(spl_ebar, E, acc_ebar);
  }

  double getFlux_mu(double E){
    if(E<E_mu[0] || E>E_mu[239])
      return 0;
    if(E>29.5 && E<30.5) // Manually add DAR
      return 120.212 + gsl_spline_eval(spl_mu, E, acc_mu);
    else
      return gsl_spline_eval(spl_mu, E, acc_mu);
  }

  double getFlux_mubar(double E){
    if(E<E_mubar_DIF[0] || E>E_mubar_DIF[249])
      return 0;
    else if(E>E_mubar_DAR[52])
      return gsl_spline_eval(spl_mubar_DIF, E, acc_mubar_DIF);
    else
      return gsl_spline_eval(spl_mubar_DIF, E, acc_mubar_DIF) +
	gsl_spline_eval(spl_mubar_DAR, E, acc_mubar_DAR);
  }

  inline double getXsec(double E){
    constexpr double alpha1 = 5.85624;
    constexpr double beta1 = -13.462;
    constexpr double gamma1 = -42.137;
    constexpr double Delta = 1.293*units::MeV;
    constexpr double sigma0 = 0.01604/SQR(units::MeV)*1.0e-7;
    constexpr double M = 938.92*units::MeV;
    if(E<Delta) //Kinematic limit
      return 0;
    else
      return sigma0 * (alpha1 + beta1*Delta/M + gamma1*(E-Delta)/M) * SQR(E-Delta);
  }
  
  /* Utilities for numerical integration */
  const double x32[16] = {0.0483076656877383162348126,0.1444719615827964934851864,0.2392873622521370745446032,0.3318686022821276497799168,0.4213512761306353453641194,0.5068999089322293900237475,0.5877157572407623290407455,0.6630442669302152009751152,0.7321821187402896803874267,0.7944837959679424069630973,0.8493676137325699701336930,0.8963211557660521239653072,0.9349060759377396891709191,0.9647622555875064307738119,0.9856115115452683354001750,0.9972638618494815635449811};
  const double w32[16] = {0.0965400885147278005667648,0.0956387200792748594190820,0.0938443990808045656391802,0.0911738786957638847128686,0.0876520930044038111427715,0.0833119242269467552221991,0.0781938957870703064717409,0.0723457941088485062253994,0.0658222227763618468376501,0.0586840934785355471452836,0.0509980592623761761961632,0.0428358980222266806568786,0.0342738629130214331026877,0.0253920653092620594557526,0.0162743947309056706051706,0.0070186100094700966004071};
  const double x4[2] = {0.3399810435848562648026658,0.8611363115940525752239465};
  const double w4[2] = {0.6521451548625461426269361,0.3478548451374538573730639};

  /* Checks that the energies are in the [18.707 MeV, 58.707 MeV] range */
  void putInLimits(double &E){
    if(E < 18.707*units::MeV)
      E = 18.707*units::MeV;
    else if(E > 58.707*units::MeV)
      E = 58.707*units::MeV;
  }
  

  /* True energy limits where the energy resolution Gaussian falls by +/- 3 sigma */
  const double Elims[11][2] = {{4.2565619141404e+07, 7.0459185017045e+07},
			       {3.4821552113311e+07, 7.0459185017045e+07},
			       {2.9344336689788e+07, 6.8504707634606e+07},
			       {2.5273907756708e+07, 5.9593843191358e+07},
			       {2.2135171003842e+07, 5.2867407476971e+07},
			       {1.9644582664609e+07, 4.7604356632775e+07},
			       {1.7622607816678e+07, 4.3370323615519e+07},
			       {1.5950124398539e+07, 3.9887854972482e+07},
			       {1.4545036833160e+07, 3.6971328105424e+07},
			       {1.3552771987811e+07, 3.4491805056825e+07},
			       {1.3552771987811e+07, 3.2356893994149e+07}
  };

  /* Returns the energy resolution factor */
  inline double eres(int n_bin, double L, double Etrue){
    double Erecomin = L / bins[n_bin][1],
      Erecomax = L / bins[n_bin][0];

    putInLimits(Erecomin);
    putInLimits(Erecomax);

    /* Work in MeV */
    Erecomin /= units::MeV;
    Erecomax /= units::MeV;
    Etrue /= units::MeV;

    double Escale = 0.975;
    return 0.5*(erf(1.51515*(Etrue-Escale*Erecomin)/sqrt(Etrue))
		- erf(1.51515*(Etrue-Escale*Erecomax)/sqrt(Etrue)));
  }
  
public:
  /* Constructor */
  LSND(): prb(OSC_DECAY_MAJORANA){ // JK
// JK  LSND(bool pvector, bool majorana): prb(pvector, majorana){
    acc_e = gsl_interp_accel_alloc();
    acc_ebar = gsl_interp_accel_alloc();
    acc_mu = gsl_interp_accel_alloc();
    acc_mubar_DIF = gsl_interp_accel_alloc();
    acc_mubar_DAR = gsl_interp_accel_alloc();
    acc_finalFlux = gsl_interp_accel_alloc();

    /* Read the data files */
    std::ifstream ifileDAR, ifileDIF;
    double EDAR, yDAR, EDIF, yDIF;
      
    ifileDAR.open("glb/lsnd-ivan/DAR_e.dat"); // JK
    ifileDIF.open("glb/lsnd-ivan/DIF_e.dat"); // JK
    for(int i=0; i<249; ++i){
      ifileDAR>>EDAR>>yDAR;
      ifileDIF>>EDIF>>yDIF;

      E_e[i] = EDAR*units::MeV;
      f_e[i] = yDAR + yDIF*1.0e-3;
    }
    ifileDAR.close();
    ifileDIF.close();

    ifileDAR.open("glb/lsnd-ivan/DAR_ebar.dat"); // JK
    ifileDIF.open("glb/lsnd-ivan/DIF_ebar.dat"); // JK
    for(int i=0; i<249; ++i){
      ifileDAR>>EDAR>>yDAR;
      ifileDIF>>EDIF>>yDIF;

      E_ebar[i] = EDAR*units::MeV;
      f_ebar[i] = yDAR + yDIF*1.0e-3;
    }
    ifileDAR.close();
    ifileDIF.close();

    /* For the muon channel, just include DIF (DAR is monochromatic) */
    ifileDIF.open("glb/lsnd-ivan/DIF_mu.dat");
    for(int i=0; i<240; ++i){
      ifileDIF>>EDIF>>yDIF;

      E_mu[i] = EDIF*units::MeV;
      f_mu[i] = yDIF*1.0e-3;
    }
    ifileDIF.close();

    ifileDAR.open("glb/lsnd-ivan/DAR_mubar.dat");
    ifileDIF.open("glb/lsnd-ivan/DIF_mubar.dat");
    for(int i=0; i<250; ++i){
      ifileDIF>>EDIF>>yDIF;

      E_mubar_DIF[i] = EDIF*units::MeV;
      f_mubar_DIF[i] = yDIF*1.0e-3;
      if(i<53){
	ifileDAR>>EDAR>>yDAR;
	
	E_mubar_DAR[i] = EDAR*units::MeV;
	f_mubar_DAR[i] = yDAR;
      }
    }
    ifileDAR.close();
    ifileDIF.close();

    /* Load the interpolators */
    spl_e = gsl_spline_alloc(gsl_interp_cspline, 249);
    gsl_spline_init(spl_e, E_e, f_e, 249);

    spl_ebar = gsl_spline_alloc(gsl_interp_cspline, 249);
    gsl_spline_init(spl_ebar, E_ebar, f_ebar, 249);

    spl_mu = gsl_spline_alloc(gsl_interp_cspline, 240);
    gsl_spline_init(spl_mu, E_mu, f_mu, 240);

    spl_mubar_DIF = gsl_spline_alloc(gsl_interp_cspline, 250);
    gsl_spline_init(spl_mubar_DIF, E_mubar_DIF, f_mubar_DIF, 250);

    spl_mubar_DAR = gsl_spline_alloc(gsl_interp_cspline, 53);
    gsl_spline_init(spl_mubar_DAR, E_mubar_DAR, f_mubar_DAR, 53);

    spl_finalFlux_L01 = gsl_spline_alloc(gsl_interp_linear, N_PTS_INTERP),
    spl_finalFlux_L02 = gsl_spline_alloc(gsl_interp_linear, N_PTS_INTERP),
    spl_finalFlux_L11 = gsl_spline_alloc(gsl_interp_linear, N_PTS_INTERP),
    spl_finalFlux_L12 = gsl_spline_alloc(gsl_interp_linear, N_PTS_INTERP);    

    /* Load the initial state */
    static constexpr int N_PTS = 5000;
    marray<double,1> E_range = regeneration::linspace(1*units::MeV, 250*units::MeV, N_PTS);
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

  std::array<double, 11> getSpectrum(const Param &prm, bool includeBkg = true){
    prb.setParam(prm);

    /* Precompute the spectrum for the 4 different values of L that are used in the integral over L*/
    double A_L = (Lmax-Lmin)/2.;
    double B_L = (Lmax+Lmin)/2.;
    double L_01 = A_L*x4[0] + B_L, L_02 = -A_L*x4[0] + B_L, L_11 = A_L*x4[1] + B_L, L_12 = -A_L*x4[1] + B_L;

    double energies[N_PTS_INTERP];
    double fluxes_01[N_PTS_INTERP], fluxes_02[N_PTS_INTERP], fluxes_11[N_PTS_INTERP], fluxes_12[N_PTS_INTERP]; //Electron antineutrino fluxes
    
    for(unsigned int i=0; i<N_PTS_INTERP; ++i){
      energies[i] = Elims[10][0] + (Elims[0][1] - Elims[10][0])*i/((double) N_PTS_INTERP - 1);
      fluxes_01[i] = prb.getFinalFlux(energies[i], L_01, lsnd_lowpass_width)[1][0];
      fluxes_02[i] = prb.getFinalFlux(energies[i], L_02, lsnd_lowpass_width)[1][0];
      fluxes_11[i] = prb.getFinalFlux(energies[i], L_11, lsnd_lowpass_width)[1][0];
      fluxes_12[i] = prb.getFinalFlux(energies[i], L_12, lsnd_lowpass_width)[1][0];
    }
      
    gsl_spline_init(spl_finalFlux_L01, energies, fluxes_01, N_PTS_INTERP);
    gsl_spline_init(spl_finalFlux_L02, energies, fluxes_02, N_PTS_INTERP);
    gsl_spline_init(spl_finalFlux_L11, energies, fluxes_11, N_PTS_INTERP);
    gsl_spline_init(spl_finalFlux_L12, energies, fluxes_12, N_PTS_INTERP);

    std::array<double, 11> result;

    for(size_t n=0; n<11; ++n){
      double theor;
      
      /* Do the integral: convolve final flux*xsec*energy resolution, averaging over L */
      double A_E = (Elims[n][1] - Elims[n][0])/2.;
      double B_E = (Elims[n][1] + Elims[n][0])/2.;
      double A_L = (Lmax-Lmin)/2.;
      theor = 0;

      for(int i_E=0; i_E<16; ++i_E){
	double E_1 = A_E*x32[i_E] + B_E;
	double E_2 = -A_E*x32[i_E] + B_E;
	double flx_11_0 = gsl_spline_eval(spl_finalFlux_L01, E_1, acc_finalFlux),
	  flx_12_0 = gsl_spline_eval(spl_finalFlux_L02, E_1, acc_finalFlux),
	  flx_11_1 = gsl_spline_eval(spl_finalFlux_L11, E_1, acc_finalFlux),
	  flx_12_1 = gsl_spline_eval(spl_finalFlux_L12, E_1, acc_finalFlux),
	  flx_21_1 = gsl_spline_eval(spl_finalFlux_L11, E_2, acc_finalFlux),
	  flx_22_1 = gsl_spline_eval(spl_finalFlux_L12, E_2, acc_finalFlux),
	  flx_21_0 = gsl_spline_eval(spl_finalFlux_L01, E_2, acc_finalFlux),
	  flx_22_0 = gsl_spline_eval(spl_finalFlux_L02, E_2, acc_finalFlux);
	    
	theor += A_E * A_L * w4[0] * w32[i_E] *
	  ( flx_11_0 * getXsec(E_1) * eres(n, L_01, E_1) / SQR(L_01) +
	    flx_12_0 * getXsec(E_1) * eres(n, L_02, E_1) / SQR(L_02) +
	    flx_21_0 * getXsec(E_2) * eres(n, L_01, E_2) / SQR(L_01) +
	    flx_22_0 * getXsec(E_2) * eres(n, L_02, E_2) / SQR(L_02) ) +
	  A_E * A_L * w4[1] * w32[i_E] *
	  ( flx_11_1 * getXsec(E_1) * eres(n, L_11, E_1) / SQR(L_11) +
	    flx_12_1 * getXsec(E_1) * eres(n, L_12, E_1) / SQR(L_12) +
	    flx_21_1 * getXsec(E_2) * eres(n, L_11, E_2) / SQR(L_11) +
	    flx_22_1 * getXsec(E_2) * eres(n, L_12, E_2) / SQR(L_12) );
      }
      
      theor /= (1/Lmin - 1/Lmax);
      
      /* Multiply by the efficiency */
      theor *= eff[n];
      
      /* Add the background */
      if(includeBkg)
	theor += bkg[n];

      result[n] = theor;
    }

    return result;    
  }
  double getChisq(const Param &prm){
    std::array<double, 11> pred = getSpectrum(prm, false);

    double chqMin = 9.99e9;
    /* Marginalise over eta */

    for(double eta = 1-0.7; eta <= 1+0.7; eta += 0.05)
      for(double rho = 1-0.6; rho <= 1+0.1; rho += 0.02){
	double chq = 0;
	double tot = 0;
	for(size_t n = 0; n<pred.size(); ++n){
	  double theor = rho*pred[n];
	  chq += 2*(eta*theor + bkg[n] - data[n] + (data[n]+no_beam[n])*
		    log((data[n]+no_beam[n]) / (eta*theor + bkg[n] + no_beam[n])));
	  tot += theor- intrinsic[n];
	}

	chq += SQR((tot/16713.2 - 0.264e-2)/0.040e-2); // JK
	  
	chq += SQR((rho-1)/0.2); // JK

	if(chq < chqMin)
	  chqMin = chq;
      }
    
    return chqMin;
  }

  int print_spectrum(const Param &prm){   // JK
    std::array<double, 11> pred = getSpectrum(prm, true);
    for (size_t n=0; n < pred.size(); n++){
      printf("# LSNDSPECT  %ld  %g\n", n, pred[n]);
    }
    return 0;
  }

  /* Destructor */
  ~LSND(){
    gsl_spline_free(spl_e);
    gsl_interp_accel_free(acc_e);
    gsl_spline_free(spl_ebar);
    gsl_interp_accel_free(acc_ebar);
    gsl_spline_free(spl_mu);
    gsl_interp_accel_free(acc_mu);
    gsl_spline_free(spl_mubar_DIF);
    gsl_interp_accel_free(acc_mubar_DIF);
    gsl_spline_free(spl_mubar_DAR);
    gsl_interp_accel_free(acc_mubar_DAR);    
  }

  /* No copy (for simplicity) */
  LSND(const LSND&) = delete;
  LSND &operator=(const LSND&) = delete;
};

} // end namespace ns_sbl (JK)

#endif

#endif