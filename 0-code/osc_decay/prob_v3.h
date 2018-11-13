#ifndef regProb_H
#define regProb_H

#include <array>
#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <nusquids_decay/nusquids_decay.h>

using namespace nusquids;
namespace regeneration{
  using cmplx = std::complex<double>;
  using cmplxMat = std::array<std::array<cmplx, 4>, 4>;
  using uint = unsigned int;

  /*
   * The code works internally in eV. Here are some useful conversion factors.
   */
  struct units{
    static constexpr double meter = 5.06773093741e6;        // [eV^-1/m]
    static constexpr double km = 1.0e3*meter;               // [eV^-1/km]

    static constexpr double KeV = 1.0e3;                    // [eV/KeV]
    static constexpr double MeV = 1.0e6;                    // [eV/MeV]
    static constexpr double GeV = 1.0e9;                    // [eV/GeV]
  };

  /*
   * Class that encapsulates all the oscillation parameters
   */
  class Param{
  public:
    double snq12, snq13, snq23, snq14, snq24, snq34;
    double dlt24, dlt13, dlt12;
    double dmq41; /* eV^2 */
    std::vector<double> m_nu; /* Vector of neutrino masses (eV) */

    /* If we denote as g_s the bar(nu)_s nu_s coupling constant, then
         g = g_s m_4/m_A
       */
    double g; 
    cmplxMat U;
    
    void loadU(){
      cmplxMat R12 = getEulerMatrix(0, 1, asin(sqrt(snq12)), dlt12);
      cmplxMat R13 = getEulerMatrix(0, 2, asin(sqrt(snq13)), dlt13);
      cmplxMat R23 = getEulerMatrix(1, 2, asin(sqrt(snq23)));
      cmplxMat R14 = getEulerMatrix(0, 3, asin(sqrt(snq14)));
      cmplxMat R24 = getEulerMatrix(1, 3, asin(sqrt(snq24)), dlt24);
      cmplxMat R34 = getEulerMatrix(2, 3, asin(sqrt(snq34)));

      U = prod(R13, R12);
      U = prod(R23, U);
      U = prod(R14, U);
      U = prod(R24, U);
      U = prod(R34, U);
    }

    inline cmplx getProjector(uint i, uint j, int i_nu){
      /* If the decay product is lambda_k |k>, then the projector
	 is given by lambda_k lambda_l* |k><l|. */
      if(i==3 || j==3)
	return 0.;      

      double normalisation = 0;
      for(uint j=0; j<3; ++j)
	normalisation += std::norm(U[3][j]);

      if(i_nu == 0)
	return U[3][i]*std::conj(U[3][j])/normalisation;
      else
	return std::conj(U[3][i])*U[3][j]/normalisation;
    }

    double get_g(double m4Gamma, bool pscalar = true){
      loadU();
      double m4GammaOverG2 = 0;
      for(uint j=0; j<3; ++j){
	if(pscalar){
	  double g_4j = std::abs(U[3][j] * U[3][3]) * (m_nu[3]+m_nu[j])/m_nu[3];
	  m4GammaOverG2 += m_nu[3] * 1.0/(16.0*M_PI)*g_4j*g_4j*
	    (_g(m_nu[3],m_nu[j]) + _k(m_nu[3],m_nu[j]));
	}
	else{
	  double g_4j = std::abs(U[3][j] * U[3][3]) * (m_nu[3]-m_nu[j])/m_nu[3];
	  m4GammaOverG2 += m_nu[3] * 1.0/(16.0*M_PI)*g_4j*g_4j*
	    (_f(m_nu[3],m_nu[j]) + _k(m_nu[3],m_nu[j]));
	}
      }
      return sqrt(m4Gamma/m4GammaOverG2);
    }
    
  private:
    inline cmplxMat getEulerMatrix(int i, int j, double theta, double delta=0){
      cmplxMat U;
      for(uint i=0; i<4; ++i)
	U[i][i] = 1.;
    
      U[i][i] = U[j][j] = cos(theta);
      U[i][j] = std::polar(sin(theta), -delta);
      U[j][i] = std::polar(sin(theta), delta + M_PI);
      return U;
    }

    inline cmplxMat prod(cmplxMat A, cmplxMat B){
      cmplxMat C;
      for(int i=0; i<4; ++i)
	for(int j=0; j<4; ++j)
	  for(int k=0; k<4; ++k)
	    C[i][j] += A[i][k]*B[k][j];
    
      return C;
    }
    
    // Auxiliary functions (see eq. 4 in arXiv:1711.05921)
    inline double _f(double m_i, double m_j) const {
      double result = m_i/2.0 + 2.0*m_j + (2.0*m_j/m_i)*(m_j*log(m_i) - log(pow(m_j,m_j))) - (2.0*m_j*m_j*m_j/(m_i*m_i)) - (m_j*m_j*m_j*m_j/(2.0*m_i*m_i*m_i));
      return result;
    }
    inline double _g(double m_i, double m_j) const {
      double result = m_i/2.0 - 2.0*m_j + (2.0*m_j/m_i)*(m_j*log(m_i) - log(pow(m_j,m_j))) + (2.0*m_j*m_j*m_j/(m_i*m_i)) - (m_j*m_j*m_j*m_j/(2.0*m_i*m_i*m_i));
      return result;
    }
    inline double _k(double m_i, double m_j) const {
      double result = m_i/2.0 - (2.0*m_j/m_i)*(m_j*log(m_i) - log(pow(m_j,m_j))) - (m_j*m_j*m_j*m_j/(2.0*m_i*m_i*m_i));
      return result;
    }
  
  };

  /*
   * Main class for calculating oscillation probabilities in a 3+1 framework,
   * where the mass squared differences among the lightest states as well as the matter
   * effects are considered to be negligible
   */
  class regProb{
  private:
    // Whether neutrinos are Majorana (i.e.,the helicity flipping decays give antineutrinos)
    const bool majorana;
    // Pseudoscalar (true) or scalar (false) coupling
    const bool pscalar;

    // Oscillation parameters
    Param prm;

    /* Initial state.
     * inistate[i_E][i_nu][i_fla] gives the initial flux of:
     *  *neutrinos(i_nu=0) / antineutrinos(i_nu=1)
     *  *in the energy node i_E
     *  *with flavour i_fla:
     *      i_fla==0: [anti]nu_e
     *      i_fla==1: [anti]nu_mu
     *      i_fla==2: [anti]nu_tau
     *      i_fla==3: [anti]nu_s
     */ 
    marray<double, 3> inistate;

    /* Initial density matrix.
     * 
     * inirho[i_E][i_nu][i][j] is the (i,j) component, in the mass basis
     *   of the density matrix:
     *    *for neutrinos(i_nu=0) / antineutrinos(i_nu)
     *    *in the energy node i_E
     *
     * only j<=i is filled and used
     */
    marray<cmplx, 4> inirho;

    // 1-D marray containing the energy corresponding to each energy node
    marray<double,1> E_range;

    // nu_4->nu_i rest frame decay widths
    double gammaCPP[3], gammaCVP[3];

    // Linearly interpolates a density matrix to get it at an arbitrary energy E
    inline cmplx interpolateRho(const marray<cmplx, 4> &rho, double E, uint i_nu, uint i, uint j) const{
      //find bracketing state entries
      auto xit=std::lower_bound(E_range.begin(),E_range.end(),E);
      if(xit!=E_range.begin())
	xit--;
      size_t xid=std::distance(E_range.begin(),xit);
      //linearly interpolate between the two states
      double f2=((E-E_range[xid])/(E_range[xid+1]-E_range[xid]));
      double f1=1-f2;
      return f1*rho[xid][i_nu][i][j] + f2*rho[xid+1][i_nu][i][j];      
    }

    // Linearly interpolates a flavour state to get it at an arbitrary energy E
    inline double interpolateState(const marray<double, 3> & state, double E, uint i_nu, uint i_fla){
      //find bracketing state entries
      auto xit=std::lower_bound(E_range.begin(),E_range.end(),E);
      if(xit!=E_range.begin())
	xit--;
      size_t xid=std::distance(E_range.begin(),xit);
      //linearly interpolate between the two states
      double f2=((E-E_range[xid])/(E_range[xid+1]-E_range[xid]));
      double f1=1-f2;
      return f1*state[xid][i_nu][i_fla] + f2*state[xid+1][i_nu][i_fla];
    }

    // Set the initial state density matrix
    void setIniRho() {
      inirho = marray<cmplx, 4>({inistate.extent(0), 2, 4, 4});

      for(size_t i_E = 0; i_E<inirho.extent(0); ++i_E)
	for(size_t i_nu = 0; i_nu<inirho.extent(1); ++i_nu)
	  for(int i=0; i<4; ++i)
	    for(int j=0; j<=i; ++j){
	      inirho[i_E][i_nu][i][j] = 0;
	      for(int i_fla=0; i_fla<4; ++i_fla)
		inirho[i_E][i_nu][i][j] += (i_nu==0)? std::conj(prm.U[i_fla][i])*prm.U[i_fla][j] * inistate[i_E][i_nu][i_fla] :
		  prm.U[i_fla][i]*std::conj(prm.U[i_fla][j]) * inistate[i_E][i_nu][i_fla];
	    }
    }

    // Auxiliary functions (see eq. 4 in arXiv:1711.05921)
    inline double f(double m_i, double m_j) const {
      double result = m_i/2.0 + 2.0*m_j + (2.0*m_j/m_i)*(m_j*log(m_i) - log(pow(m_j,m_j))) - (2.0*m_j*m_j*m_j/(m_i*m_i)) - (m_j*m_j*m_j*m_j/(2.0*m_i*m_i*m_i));
      return result;
    }
    inline double g(double m_i, double m_j) const {
      double result = m_i/2.0 - 2.0*m_j + (2.0*m_j/m_i)*(m_j*log(m_i) - log(pow(m_j,m_j))) + (2.0*m_j*m_j*m_j/(m_i*m_i)) - (m_j*m_j*m_j*m_j/(2.0*m_i*m_i*m_i));
      return result;
    }
    inline double k(double m_i, double m_j) const {
      double result = m_i/2.0 - (2.0*m_j/m_i)*(m_j*log(m_i) - log(pow(m_j,m_j))) - (m_j*m_j*m_j*m_j/(2.0*m_i*m_i*m_i));
      return result;
    }

    // Computes the **rest frame** decay widths
    void compute_Gammas(){
      if (!pscalar){
	//CPP, scalar
	for(uint j=0; j<3; ++j) {
	  double g_4j = prm.g * std::abs(prm.U[3][j] * prm.U[3][3])
	    * (prm.m_nu[3]-prm.m_nu[j])/prm.m_nu[3];
	  gammaCPP[j] = 1.0/(16.0*M_PI)*g_4j*g_4j*f(prm.m_nu[3],prm.m_nu[j]);	  
	}
	//CVP, scalar
	for(uint j=0; j<3; ++j) {
	  double g_4j = prm.g * std::abs(prm.U[3][j] * prm.U[3][3])
	    * (prm.m_nu[3]-prm.m_nu[j])/prm.m_nu[3];
	  gammaCVP[j] = 1.0/(16.0*M_PI)*g_4j*g_4j*k(prm.m_nu[3],prm.m_nu[j]);	  
	}
      }
      if (pscalar){
	//CPP, pscalar
	for(uint j=0; j<3; ++j) {
	  double g_4j = prm.g * std::abs(prm.U[3][j] * prm.U[3][3])
	    * (prm.m_nu[3]+prm.m_nu[j])/prm.m_nu[3];
	  gammaCPP[j] = 1.0/(16.0*M_PI)*g_4j*g_4j*g(prm.m_nu[3],prm.m_nu[j]);
	}
	//CVP, pscalar
	for(uint j=0; j<3; ++j) {
	  double g_4j = prm.g * std::abs(prm.U[3][j] * prm.U[3][3])
	    * (prm.m_nu[3]+prm.m_nu[j])/prm.m_nu[3];
	  gammaCVP[j] = 1.0/(16.0*M_PI)*g_4j*g_4j*k(prm.m_nu[3],prm.m_nu[j]);	  
	}
      }
    }

    // Gives the total rest-frame nu_4 decay width
    inline double getGamma() const{
      double res = 0;
      for(uint j=0; j<3; ++j) //Sum incoherently over final states
	res += gammaCPP[j] + gammaCVP[j];
      return res;
    }

    //! Returns the spectral integrand without the projector
    /*!
      If majorana is true, there are regeneration contributions from both CVP and CPP.
      Otherwise, there is no contribution because right-handed neutrinos are sterile.
      Note that the contribution from CVP in the majorana case converts neutrinos to antineutrinos.
    */
    double spectrum(double eparent, double L, uint j, double edaughter, const marray<cmplx, 4> &inirho, uint i_nu) const{
      if(getGamma()<1.0e-14)
	return 0.;
      
      uint i=3;
      double gamma = eparent/prm.m_nu[i];

      double res = 0;
      // Interpolate
      double rho44 = std::real(interpolateRho(inirho, eparent, i_nu, 3, 3));
      uint i_nu_other = (i_nu==0)? 1 : 0;
      double rho44_other = std::real(interpolateRho(inirho, eparent, i_nu_other, 3, 3));
      if(fabs(prm.m_nu[j]-0.0)>1e-6){
	double xij = prm.m_nu[i]/prm.m_nu[j];
	if (!pscalar)
	  res += rho44 * eparent/(prm.m_nu[3]*getGamma()) *
	    (1.-exp(-prm.m_nu[3]*getGamma()*L/eparent))*
	    (xij*xij/(xij*xij-1))*
	    (1/(eparent*eparent*edaughter))*
	    (gammaCPP[j]/gamma)*
	    pow(eparent+xij*edaughter,2)/pow(xij+1,2);
	if (pscalar)
	  res += rho44 * eparent/(prm.m_nu[3]*getGamma()) *
	    (1.-exp(-prm.m_nu[3]*getGamma()*L/eparent))*
	    (xij*xij/(xij*xij-1))*
	    (1/(eparent*eparent*edaughter))*
	    (gammaCPP[j]/gamma)*
	    pow(eparent-xij*edaughter,2)/pow(xij-1,2);
	if(majorana){
	  //Include chirality-violating term. Majorana CVP sends nu to nubar.
	  //The procedure is the same, but the parent i_nu index is inverted.
	  if (!pscalar)
	    res += rho44_other * eparent/(prm.m_nu[3]*getGamma()) *
	      (1.-exp(-prm.m_nu[3]*getGamma()*L/eparent))*
	      (xij*xij/(xij*xij-1))*
	      ((eparent-edaughter)/(eparent*eparent*edaughter))*
	      ((gammaCVP[j]/gamma)*
	       (edaughter*pow(xij,2)-eparent)/pow(xij+1,2));
	  if (pscalar)
	    res += rho44_other * eparent/(prm.m_nu[3]*getGamma()) *
	      (1.-exp(-prm.m_nu[3]*getGamma()*L/eparent))*
	      (xij*xij/(xij*xij-1))*
	      ((eparent-edaughter)/(eparent*eparent*edaughter))*
	      ((gammaCVP[j]/gamma)*
	       (edaughter*pow(xij,2)-eparent)/pow(xij-1,2));
	}
      } else{
	//If prm.m_nu[j] is too close to zero, xij diverges, and we switch to an alternative
	//form for the differential decay rates, in terms of yij=1/xij.
	double yij = prm.m_nu[j]/prm.m_nu[i];
	if (!pscalar)
	  res += rho44 * eparent/(prm.m_nu[3]*getGamma()) *
	    (1.-exp(-prm.m_nu[3]*getGamma()*L/eparent))*
	    (1/(1-yij*yij))*
	    (1/(eparent*eparent*edaughter))*
	    (gammaCPP[j]/gamma)*
	    pow(eparent*yij+edaughter,2)/pow(yij+1,2);
	if (pscalar)
	  res += rho44 * eparent/(prm.m_nu[3]*getGamma()) *
	    (1.-exp(-prm.m_nu[3]*getGamma()*L/eparent))*
	    (1/(1-yij*yij))*
	    (1/(eparent*eparent*edaughter))*
	    (gammaCPP[j]/gamma)*
	    pow(eparent*yij-edaughter,2)/pow(1-yij,2);
	if(majorana){
	  //Include chirality-violating term. Majorana CVP sends nu to nubar.
	  //The procedure is the same, but the parent irho index is inverted.
	  if (!pscalar)
	    res += rho44_other * eparent/(prm.m_nu[3]*getGamma()) *
	      (1.-exp(-prm.m_nu[3]*getGamma()*L/eparent))*
	      (1/(1-yij*yij))*
	      ((eparent-edaughter)/(eparent*eparent*edaughter))*
	      ((gammaCVP[j]/gamma)*
	       (edaughter-pow(yij,2)*eparent)/pow(yij+1,2));
	  if (pscalar)
	    res += rho44_other * eparent/(prm.m_nu[3]*getGamma()) *
	      (1.-exp(-prm.m_nu[3]*getGamma()*L/eparent))*
	      (1/(1-yij*yij))*
	      ((eparent-edaughter)/(eparent*eparent*edaughter))*
	      ((gammaCVP[j]/gamma)*
	       (edaughter-pow(yij,2)*eparent)/pow(1-yij,2));
	}
      }
      return res;
    }

    /* Abscisae and weights for the 16-point rule integration */
    const double x16[8] = {0.0950125098376374401853193,0.2816035507792589132304605,0.4580167776572273863424194,0.6178762444026437484466718,0.7554044083550030338951012,0.8656312023878317438804679,0.9445750230732325760779884,0.9894009349916499325961542};
    const double w16[8] = {0.1894506104550684962853967,0.1826034150449235888667637,0.1691565193950025381893121,0.1495959888165767320815017,0.1246289712555338720524763,0.0951585116824927848099251,0.0622535239386478928628438,0.0271524594117540948517806};
    
    //! Returns the spectral integral without the projector
    double getIntegral(double edaughter, double L, const marray<cmplx, 4> &inirho, uint i_nu) const{
      // Sum over final states
      double result = 0;
      for(uint j= 0; j < 3; ++j){
	uint i = 3;
	double xij = prm.m_nu[i]/prm.m_nu[j];

	double eparent_low = edaughter;
	double eparent_high = xij*xij*edaughter;
	//If upper limit is larger than last element, cut the integral in the last element
	if(eparent_high > E_range[E_range.size()-1]) 
	  eparent_high = E_range[E_range.size()-1];

	double A = (eparent_high-eparent_low)/2.;
	double B = (eparent_high+eparent_low)/2.;

	for(int i_E = 0; i_E<8; ++i_E)
	  result += A * w16[i_E] * (spectrum(A*x16[i_E] + B, L, j, edaughter, inirho, i_nu) +
				    spectrum(-A*x16[i_E] + B, L, j, edaughter, inirho, i_nu));
	    
      }

      return result;
    }
  public:
    /* Constructor. The parameters are:
     *  * E_range_: 1-D marray containing the energy corresponding to each energy node
     *  * inistate_: 3-D marray. inistate[i_E][i_nu][i_fla] gives the initial flux of:
     *               *neutrinos(i_nu=0) / antineutrinos(i_nu=1)
     *               *in the energy node i_E
     *               *with flavour i_fla:
     *                     i_fla==0: [anti]nu_e
     *                     i_fla==1: [anti]nu_mu
     *                     i_fla==2: [anti]nu_tau
     *                     i_fla==3: [anti]nu_s
     *  * pscalar_: pseudoscalar (true) or scalar (false) decay
     *  * majorana_: whether neutrinos are Majorana or not.
     *               If true, there is a chirality-flipping regeneration contribution
     */
    regProb(marray<double,1> E_range_, marray<double, 3>& inistate_, bool pscalar_, bool majorana_):
      majorana(majorana_), pscalar(pscalar_),
      inistate(inistate_), E_range(E_range_){
      if(inistate.extent(1)!=2 || inistate.extent(2) != 4){
	std::cerr<<"Wrong dimensions for the initial state!"<<std::endl;
	exit(1);
      }
    }

    /* Constructor **WITHOUT** the initial state and energy range.
     * If this constructor is used, one **MUST** call setIniState before using
     * the object
     */
    regProb(bool pscalar_, bool majorana_):
      majorana(majorana_), pscalar(pscalar_){}

    // Modifies the initial state
    void setIniState(const marray<double, 1>& E_range_, const marray<double, 3>& inistate_){
      E_range = E_range_;
      inistate = inistate_;
      if(inistate.extent(1)!=2 || inistate.extent(2) != 4){
	std::cerr<<"Wrong dimensions for the initial state!"<<std::endl;
	exit(1);
      }
      setIniRho();
    }

    // Set the oscillation parameters
    void setParam(Param prm_){
      prm = prm_;
      prm.loadU();
      compute_Gammas();
      setIniRho();
    }

    /* Returns the final state flux for a given energy E.
     *  output[i_nu][i_fla] gives the final flux of:
     *   *neutrinos(i_nu=0) / antineutrinos(i_nu=1)
     *   *with flavour i_fla:
     *      i_fla==0: [anti]nu_e
     *      i_fla==1: [anti]nu_mu
     *      i_fla==2: [anti]nu_tau
     *      i_fla==3: [anti]nu_s
     */
    marray<double, 2> getInitialFlux(double E){
      marray<double, 2> output{2,4};
      
      for(uint i_nu = 0; i_nu<2; ++i_nu)
	for(uint i_fla = 0; i_fla<4; ++i_fla)
	  output[i_nu][i_fla] = interpolateState(inistate, E, i_nu, i_fla);
      
      return output;
    }
    
    /* Returns the final state flux for a given energy E and baseline L.
     *  output[i_nu][i_fla] gives the final flux of:
     *   *neutrinos(i_nu=0) / antineutrinos(i_nu=1)
     *   *with flavour i_fla:
     *      i_fla==0: [anti]nu_e
     *      i_fla==1: [anti]nu_mu
     *      i_fla==2: [anti]nu_tau
     *      i_fla==3: [anti]nu_s
     *  IMPORTANT: THERE IS NO 1/L^2 FLUX ATTENUATION 
     */
    marray<double, 2> getFinalFlux(double E, double L){
      marray<double, 2> finstate{2, 4};
      
      /* Evolve the initial state density matrix */
      marray<cmplx, 3> finrho{2, 4, 4};

      for(uint i_nu = 0; i_nu<2; ++i_nu){
	/* First, we start with rho_{44} */
	if(prm.g > 1.0e-7)
	  finrho[i_nu][3][3] = interpolateRho(inirho, E, i_nu, 3, 3) *
	    exp(-prm.m_nu[3] * getGamma() * L / E);
	else
	  finrho[i_nu][3][3] = interpolateRho(inirho, E, i_nu, 3, 3);
	
	/* Then, rho_{4k} */
	for(int j=0; j<3; ++j){
	  if(prm.g > 1.0e-7)
	    finrho[i_nu][3][j] = interpolateRho(inirho, E, i_nu, 3, j) *
	      std::polar(exp(-prm.m_nu[3] * getGamma() * L /(2*E) ),
			 -prm.dmq41 * L /(2.*E));
	  else
	    finrho[i_nu][3][j] = interpolateRho(inirho, E, i_nu, 3, j) *
	      std::polar(1., -prm.dmq41 * L /(2.*E));
	}
	  
	/* And, finally, the 3x3 subblock */
	for(int i=0; i<3; ++i)
	  for(int j=0; j<=i; ++j){
	    if(prm.g > 1.0e-7)
	      finrho[i_nu][i][j] = interpolateRho(inirho, E, i_nu, i, j)
		+ getIntegral(E, L, inirho, i_nu)*prm.getProjector(i, j, i_nu);
	    else
	      finrho[i_nu][i][j] = interpolateRho(inirho, E, i_nu, i, j);
	  }
	  
	/* We can now convert to the flavour basis */
	for(uint i_fla=0; i_fla<inistate.extent(2); ++i_fla){
	  finstate[i_nu][i_fla] = 0;
	  for(uint i=0; i<4; ++i){
	    for(uint j=0; j<i; ++j)
	      finstate[i_nu][i_fla] += (i_nu==0)? 2*std::real( prm.U[i_fla][i]*std::conj(prm.U[i_fla][j]) *finrho[i_nu][i][j] ):
		2*std::real( std::conj(prm.U[i_fla][i])*prm.U[i_fla][j] *finrho[i_nu][i][j] );
	    
	    finstate[i_nu][i_fla] += std::norm(prm.U[i_fla][i])*std::real(finrho[i_nu][i][i]);
	  }
	}
      }

      return finstate;
    }

    /*
     * Returns the probability matrix
     * 
     * P[alpha][beta] = Prob(nu_alpha->nu_beta)
     *
     * cp_sign is +1 for neutrinos and -1 for antineutrinos
     *
     * Careful!! This probability is defined as the number of final nu_beta with energy E
     * divided by the initial nu_alpha with energy E. It has therefore two drawbacks:
     * 
     * *Probabilities can be greater than 1 if there is energy migration
     * *It is not appropriate if there are neutrino->antineutrino transitions!
     *
     * ALSO, THERE IS NO 1/L^2 FLUX ATTENUATION
     */
    inline void getProbabilityMatrix(double P[3][3], int cp_sign, double E, double L){
      uint i_nu = (cp_sign==+1)? 0 : 1;

      marray<double, 2> finstate = getFinalFlux(E, L);
      
      for(uint i_fla=0; i_fla<3; ++i_fla)
	for(uint j_fla=0; j_fla<3; ++j_fla)
	  P[i_fla][j_fla] = finstate[i_nu][j_fla] / interpolateState(inistate, E, i_nu, i_fla);
    }
  
  };

  // Simple utility giving a 1-D array of equally spaced points
  marray<double,1> linspace(double Emin,double Emax,unsigned int div){
    if(div==0)
        throw std::length_error("number of samples requested from linspace must be nonzero");
    marray<double,1> linpoints{div};
    double step_lin = (Emax - Emin)/double(div-1);
    
    double EE = Emin;
    for(unsigned int i=0; i<div-1; i++, EE+=step_lin)
        linpoints[i] = EE;
    linpoints[div-1] = Emax;
	
    return linpoints;
  }
}

#endif
