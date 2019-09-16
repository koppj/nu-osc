#ifndef regProb_H_11
#define regProb_H_11

#include <nusquids_decay/nusquids_decay.h> //JK
#include <array>
#include <cmath>
#include <complex.h> 
#include <iostream>
#include <stdexcept>
#include <gsl/gsl_sf_expint.h>
#include "snu.h" //JK

#ifdef NU_USE_NUSQUIDS

#define SQR(x)  ((x)*(x))  // square of a number
#define _I csqrt(-1)
using namespace nusquids;  //JK
namespace regeneration{
  using cmplx    = double _Complex;
  using cmplxMat = double _Complex[4][4];
  using uint     = unsigned int;

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
    double m_nu[4]; /* Vector of neutrino masses (eV) */

    double m_A; /* Mass of the boson (eV) */
    double g;
    bool pscalar;
    
    cmplxMat U;

    Param(){                            // JK - constructor
      pscalar = false;
    }
    
    void loadU(){
      cmplxMat R12, R13, R23, R14, R24, R34;
      getEulerMatrix(0, 1, asin(sqrt(snq12)), dlt12, R12);
      getEulerMatrix(0, 2, asin(sqrt(snq13)), dlt13, R13);
      getEulerMatrix(1, 2, asin(sqrt(snq23)), 0.0,   R23);
      getEulerMatrix(0, 3, asin(sqrt(snq14)), 0.0,   R14);
      getEulerMatrix(1, 3, asin(sqrt(snq24)), dlt24, R24);
      getEulerMatrix(2, 3, asin(sqrt(snq34)), 0.0,   R34);

      prod(R13, R12, U);
      prod(R23, U, U);
      prod(R14, U, U);
      prod(R24, U, U);
      prod(R34, U, U);
    }

    inline cmplx getProjector(uint i, uint j, int i_nu) const{
      /* If the decay product is lambda_k |k>, then the projector
         is given by lambda_k lambda_l* |k><l|. */
      if(i==3 || j==3)
        return 0.;      

      double normalisation = 0;
      for(uint j=0; j<3; ++j)
        normalisation += SQR(cabs(U[3][j]));

      if(i_nu == 0)
        return U[3][i]*conj(U[3][j])/normalisation;
      else
        return conj(U[3][i])*U[3][j]/normalisation;
    }

    /* Some auxiliary functions */
    inline double g_4j(uint j) const{
      return g * cabs(U[3][j] * U[3][3]);
    }
    inline double xj4(uint j) const{
      return m_nu[j]/m_nu[3];
    }
    inline double xA4() const{
      return m_A / m_nu[3];
    }
    inline double xj4_sq(uint j) const{
      return SQR(m_nu[j]) / SQR(m_nu[3]);
    }
    inline double xA4_sq() const{
      return SQR(m_A) / SQR(m_nu[3]);
    }
    inline double ALPHA(uint j) const{
      if(xA4()+xj4(j) > 1){
        std::cerr << "Unphysical parameters: mA+m_j > m_4 at m_A=" << m_A
                  << ", m4=" << m_nu[3] << ", mj=" << m_nu[j] << std::endl;
        exit(1);
      }
      return 1 - xA4_sq() + xj4_sq(j);
    }

    /* Model-dependent decay functions */

    // Returns the value of g corresponding to a given m4*Gamma
    double get_g(double m4Gamma){
      loadU();
      double m4GammaOverG2 = 0;
      for(uint j=0; j<3; ++j){
        if (m_A+m_nu[j] >= m_nu[3]) // enforce kinematic constraint (JK)
            continue;
        double g_4j = cabs(U[3][j] * U[3][3]);
        double xj4 = this->xj4(j);
        double xj4_sq = this->xj4_sq(j);
        double ALPHA = this->ALPHA(j);
        double sqr = sqrt(SQR(ALPHA)-4*xj4_sq);

        if(pscalar)
          m4GammaOverG2 += m_nu[3] * g_4j*g_4j*m_nu[3]/(16.0*M_PI) * (ALPHA - 2*xj4)*sqr;
        else
          m4GammaOverG2 += m_nu[3] * g_4j*g_4j*m_nu[3]/(16.0*M_PI) * (ALPHA + 2*xj4)*sqr;
      }
      return sqrt(m4Gamma/m4GammaOverG2);
    }

    // Returns the nu_4 -> nu_j **rest frame** partial width
    double getGamma(int j) const{
      if (m_A+m_nu[j] >= m_nu[3]) // enforce kinematic constraint (JK)
          return 0.0;
      double g_4j = this->g_4j(j);
      double xj4 = this->xj4(j);
      double xj4_sq = this->xj4_sq(j);
      double ALPHA = this->ALPHA(j);
      double sqr = sqrt(SQR(ALPHA)-4*xj4_sq);

      if(pscalar)
        return g_4j*g_4j*m_nu[3]/(16.0*M_PI) * (ALPHA - 2*xj4)*sqr;
      else
        return g_4j*g_4j*m_nu[3]/(16.0*M_PI) * (ALPHA + 2*xj4)*sqr;
    }

    // Returns the total nu_4 **rest frame** decay width
    inline double getGamma() const{
      return getGamma(0) + getGamma(1) + getGamma(2);
    }

    // Returns the total **rest frame** decay width of the scalar phi
    // to active neutrinos, assuming the latter to be massless
    // see eq. (25) in Analytic_w_scalar_decay.pdf, summed over i and j and using the
    // unitarity of the mixing matrix
    inline double getGammaPhi() const{
      double Usq = cabs(U[3][3]*U[3][3]);
      return g*g * m_A / (8.*M_PI) * SQR(1. - Usq);
    }

    // Returns the LAB FRAME nu_4 -> nu_j ChP decay spectrum
    inline double dGdELL(uint j, double eparent, double edaughter) const{
      double g_4j  = this->g_4j(j);
      double xj4   = this->xj4(j);
      if(pscalar)
        return g_4j*g_4j*m_nu[3]*m_nu[3]/(16.0*M_PI*eparent*eparent) *
          SQR(edaughter/eparent - xj4)*eparent/edaughter;
      else
        return g_4j*g_4j*m_nu[3]*m_nu[3]/(16.0*M_PI*eparent*eparent) *
          SQR(edaughter/eparent + xj4)*eparent/edaughter;
    }

    // Returns the LAB FRAME nu_4 -> nu_j ChV decay spectrum
    inline double dGdELR(uint j, double eparent, double edaughter) const{
      double g_4j = this->g_4j(j);
      double xj4_sq = this->xj4_sq(j);
      double ALPHA = this->ALPHA(j);

      return g_4j*g_4j*m_nu[3]*m_nu[3]/(16.0*M_PI*eparent*eparent) *
        (ALPHA - edaughter/eparent - eparent/edaughter*xj4_sq);
    }    

  private:
    inline void getEulerMatrix(int i, int j, double theta, double delta, cmplxMat U){
      for(uint i=0; i<4; ++i)
        for(uint j=0; j<4; ++j)
          U[i][j] = (i==j)? 1 : 0;
    
      U[i][i] = U[j][j] = cos(theta);
      U[i][j] = sin(theta) * cexp(-_I * delta);
      U[j][i] = sin(theta) * cexp(_I * (delta + M_PI));
    }

    inline void prod(cmplxMat A, cmplxMat B, cmplxMat C){
      cmplxMat T;
      for(int i=0; i<4; ++i)
        for(int j=0; j<4; ++j){
          T[i][j] = 0;
          for(int k=0; k<4; ++k)
            T[i][j] += A[i][k]*B[k][j];
        }
      memcpy(C, T, 16*sizeof(T[0][0]));
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
    double (*inistate)[2][4];
    uint E_range_size; //Number of energy nodes

    // 1 array containing the energy corresponding to each energy node
    double *E_range;

    // 1d array containing the inverse step sizes between energy nodes
    double *inv_E_steps;

    // nu_4->nu_i rest frame decay widths
    double gamma[3], totalGamma;

    // \phi decay width
    double GammaPhi;

    // find index of energy bin into which a given energy falls
    // the function returns the index of the sampling point just *above* E
    inline int findEnergyBin(double E) const{
      //reject E too far outside the support of E_range to avoid undue extrapolation
      if (E < 0.9*E_range[0] || E > 1.1*E_range[E_range_size-1])
        return -1;

      uint lo=0, hi=E_range_size-1; // find appropriate energy bin using binary search
      while (hi - lo > 1)
      {
        uint k = (hi + lo) / 2;
        if (E < E_range[k])
          hi = k;
        else
          lo = k;
      }
      return hi;
    }

    // Linearly interpolates a flavour state to get it at an arbitrary energy E
    inline double interpolateState(double state[][2][4], double E, uint i_nu, uint i_fla) const{
      int k = findEnergyBin(E);
      if (k < 0)
        return 0.0;
      return interpolateState(state, E, k, i_nu, i_fla);
    }

    // Linearly interpolate a flavour state at energy E if is alreay known
    // that the energy E falls into the energy bin [k-1,k]
    inline double interpolateState(double state[][2][4], double E, int k,
                                   uint i_nu, uint i_fla) const{
      double t = (E - E_range[k-1]) * inv_E_steps[k-1];
      return (1.-t) * state[k-1][i_nu][i_fla] + t * state[k][i_nu][i_fla];
    }

    // Returns the initial density matrix, interpolated to an arbitrary energy E
    inline cmplx getIniRho(double E, uint i_nu, uint i, uint j) const{
      int k = findEnergyBin(E);
      if (k < 0)
        return 0.0;
      return getIniRho(E, k, i_nu, i, j);
    }

    // the same if the index of the energy bin into which E falls is already known
    inline cmplx getIniRho(double E, int k, uint i_nu, uint i, uint j) const{
      cmplx res = 0;
      double t = (E - E_range[k-1]) * inv_E_steps[k-1];
      if (i_nu == 0)
        for(int i_fla=0; i_fla<4; ++i_fla)
          res += conj(prm.U[i_fla][i]) * prm.U[i_fla][j]
                   * ((1.-t) * inistate[k-1][i_nu][i_fla] + t * inistate[k][i_nu][i_fla]);
      else
        for(int i_fla=0; i_fla<4; ++i_fla)
          res += prm.U[i_fla][i] * conj(prm.U[i_fla][j])
                   * ((1.-t) * inistate[k-1][i_nu][i_fla] + t * inistate[k][i_nu][i_fla]);
      return res;
    }
    
    // Loads the **rest frame** decay widths
    void load_Gammas(){
      totalGamma = 0;
      for(uint j=0; j<3; ++j){
        gamma[j] = prm.getGamma(j);
        totalGamma += gamma[j];
      }
      GammaPhi = prm.getGammaPhi();
    }

    //! Returns the spectral integrand without the projector
    /*!
      If majorana is true, there are regeneration contributions from both CVP and CPP.
      Otherwise, there is no contribution because right-handed neutrinos are sterile.
      Note that the contribution from CVP in the majorana case converts neutrinos to antineutrinos.
    */
    double spectrum(double eparent, double L, uint j, double edaughter, uint i_nu) const{
      if(totalGamma<1.0e-14)
        return 0.;
      
      double res = 0;
      // Interpolate
      int k = findEnergyBin(eparent);
      if (k < 0)
        return 0.0;
      uint i_nu_other    = (i_nu==0)? 1 : 0;
      double rho44       = creal(getIniRho(eparent, k, i_nu, 3, 3));
      double rho44_other = creal(getIniRho(eparent, k, i_nu_other, 3, 3));

      res += rho44 * eparent/(prm.m_nu[3]*totalGamma) *
        (1.-exp(-prm.m_nu[3]*totalGamma*L/eparent))*
        prm.dGdELL(j, eparent, edaughter);
      if(majorana){
        //Include chirality-violating term. Majorana CVP sends nu to nubar.
        //The procedure is the same, but the parent i_nu index is inverted.
        res += rho44_other * eparent/(prm.m_nu[3]*totalGamma) *
          (1.-exp(-prm.m_nu[3]*totalGamma*L/eparent))*
          prm.dGdELR(j, eparent, edaughter);
      }
      return res;
    }

    /* Abscisae and weights for the 16-point rule integration */
    const double x16[8] = {0.0950125098376374401853193,0.2816035507792589132304605,0.4580167776572273863424194,0.6178762444026437484466718,0.7554044083550030338951012,0.8656312023878317438804679,0.9445750230732325760779884,0.9894009349916499325961542};
    const double w16[8] = {0.1894506104550684962853967,0.1826034150449235888667637,0.1691565193950025381893121,0.1495959888165767320815017,0.1246289712555338720524763,0.0951585116824927848099251,0.0622535239386478928628438,0.0271524594117540948517806};
    
    //! Returns the function F(t) given by eq. (30) in Analytic_w_scalar_decay.pdf.
    // This function appears in the phase space integral for neutrinos from \phi decay.
    double getPhiDecayF(double Eq, double t, double L) const{
      double m_phi = prm.m_A;
      double a_phi = m_phi * GammaPhi * L;
      double a_4   = prm.m_nu[3] * totalGamma * L;

      double exp_a4Eq = exp(-a_4/Eq);
      double result = gsl_sf_expint_Ei(-a_phi/t)  +  log(t/Eq)
                    + exp_a4Eq * log(fabs(a_phi/t - a_4/Eq));
      if (a_4/Eq - a_phi/t < 500.) // exp_a4Eq * Ei is close to zero at large a$/Eq
        result -= exp_a4Eq * gsl_sf_expint_Ei(a_4/Eq - a_phi/t);
      return result;
    }

    //! Returns the integrand of the Eq integrals in eq. (29)
    double getPhiDecayIntegrand(double Eq, double E, double L) const{
      int k = findEnergyBin(Eq);
      if (k < 0)
        return 0.0;
      double x_phi4_sq   = prm.xA4_sq();
      double rho44_nu    = creal(getIniRho(Eq, k, 0, 3, 3));
      double rho44_nubar = creal(getIniRho(Eq, k, 1, 3, 3));
      return (rho44_nu + rho44_nubar) / ((1.-x_phi4_sq) * Eq)
                  * (getPhiDecayF(Eq, Eq, L) - getPhiDecayF(Eq, E,  L));
    }

    //! Returns the spectral integral without the projector
    double getIntegral(double edaughter, double L, uint i_nu) const{

      // 1. \nu_{1,2,3} from \nu_4 decay (first line in eq. (9) in Analytic_w_scalar_decay.pdf)
      // Sum over final states
      double result = 0;
      for(uint j= 0; j < 3; ++j){
        if (prm.m_A + prm.m_nu[j] > prm.m_nu[3]) // enforce kinematic constraint (JK)
          continue;
        double xj4_sq = prm.xj4_sq(j);
        double ALPHA = prm.ALPHA(j);
        
        double eparent_low = edaughter / (0.5 * (ALPHA + sqrt(SQR(ALPHA) - 4*xj4_sq)));
        double eparent_high = edaughter / (0.5 * (ALPHA - sqrt(SQR(ALPHA) - 4*xj4_sq)));
        
        //If upper limit is larger than last element, cut the integral in the last element
        if(eparent_high > E_range[E_range_size-1]) 
          eparent_high = E_range[E_range_size-1];
        if(eparent_low < E_range[0])
          eparent_low = E_range[0];
        if(eparent_low > eparent_high)
          continue;

        double A = (eparent_high-eparent_low)/2.;
        double B = (eparent_high+eparent_low)/2.;

        for(int i_E = 0; i_E<8; ++i_E)
          result += A * w16[i_E] * (spectrum(A*x16[i_E] + B, L, j, edaughter, i_nu) +
                                    spectrum(-A*x16[i_E] + B, L, j, edaughter, i_nu));
      }

      // 2. \nu_{1,2,3} from \phi decay (second line in eq. (9) in Analytic_w_scalar_decay.pdf,
      // using eq. (29) for simplification)
      double x_phi4_sq = prm.xA4_sq();
      double result_phi_decay = 0.0;
      double A, B;
      A = 0.5 * (edaughter/x_phi4_sq - edaughter);
      B = 0.5 * (edaughter/x_phi4_sq + edaughter);
      for(int i_E = 0; i_E<8; ++i_E){
        double Eq1 =  A * x16[i_E] + B;
        double Eq2 = -A * x16[i_E] + B;
        result_phi_decay += A*w16[i_E] * (getPhiDecayIntegrand(Eq1,edaughter,L)
                                        + getPhiDecayIntegrand(Eq2,edaughter,L));
      }
      A = 0.5; // we map \int_{E/x_phi4_sq}^\infty dE_q  onto the interval [0,1] and evaluate
      B = 0.5; // \int_0^1 dt f(E/x_phi4_sq + E*(1-t)/t) * E / t^2
      for(int i_E = 0; i_E<8; ++i_E){
        double t1 =  A * x16[i_E] + B;
        double t2 = -A * x16[i_E] + B;
        double Eq1 = edaughter/x_phi4_sq + edaughter*(1.-t1) / t1;
        double Eq2 = edaughter/x_phi4_sq + edaughter*(1.-t2) / t2;
        result_phi_decay += A*w16[i_E]
             * (getPhiDecayIntegrand(Eq1,x_phi4_sq*Eq1,L)*edaughter/(t1*t1)
              + getPhiDecayIntegrand(Eq2,x_phi4_sq*Eq2,L)*edaughter/(t2*t2));
      }

      return result + result_phi_decay;
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
     *  * majorana_: whether neutrinos are Majorana or not.
     *               If true, there is a chirality-flipping regeneration contribution
     */
    regProb(marray<double,1> E_range_, marray<double, 3>& inistate_, bool majorana_):
      majorana(majorana_){
      if(inistate_.extent(1)!=2 || inistate_.extent(2) != 4){
        std::cerr<<"Wrong dimensions for the initial state!"<<std::endl;
        exit(1);
      }
      E_range_size = E_range_.extent(0);

      /* Initialize the arrays */
      inistate    = new double [E_range_size][2][4];
      E_range     = new double[E_range_size];
      inv_E_steps = new double[E_range_size-1];

      /* And fill them */
      for(uint i=0; i<E_range_size; ++i){
        E_range[i] = E_range_[i];
        for(int j=0; j<2; ++j)
          for(int k=0; k<4; ++k)
            inistate[i][j][k] = inistate_[i][j][k];
      }
      for(uint i=0; i<E_range_size-1; ++i)
        inv_E_steps[i] = 1./(E_range[i+1] - E_range[i]);
    }

    /* Constructor **WITHOUT** the initial state and energy range.
     * If this constructor is used, one **MUST** call setIniState before using
     * the object
     */
    regProb(bool majorana_):
      majorana(majorana_){

      inistate    = NULL;
      E_range     = NULL;
      inv_E_steps = NULL;
    }

    // Modifies the initial state
    void setIniState(const marray<double, 1>& E_range_, const marray<double, 3>& inistate_){
      if(inistate_.extent(1)!=2 || inistate_.extent(2) != 4){
        std::cerr<<"Wrong dimensions for the initial state!"<<std::endl;
        exit(1);
      }

      E_range_size = E_range_.extent(0);

      /* Initialize the arrays */
      if (inistate)    delete[] inistate;
      if (E_range)     delete[] E_range;
      if (inv_E_steps) delete[] inv_E_steps;
      inistate    = new double [E_range_size][2][4];
      E_range     = new double[E_range_size];
      inv_E_steps = new double[E_range_size-1];

      /* And fill them */
      for(uint i=0; i<E_range_size; ++i){
        E_range[i] = E_range_[i];
        for(int j=0; j<2; ++j)
          for(int k=0; k<4; ++k)
            inistate[i][j][k] = inistate_[i][j][k];
      }
      for(uint i=0; i<E_range_size-1; ++i)
        inv_E_steps[i] = 1./(E_range[i+1] - E_range[i]);
    }

    // Set the oscillation parameters
    void setParam(Param prm_){
      prm = prm_;
      prm.loadU();
      load_Gammas();
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
     *   *filter: width of GLoBES-style low-pass filter, see GLoBES manual for details
     *  IMPORTANT: THERE IS NO 1/L^2 FLUX ATTENUATION 
     */
    inline marray<double, 2> getFinalFlux(double E, double L, const double filter_value=0.){
      marray<double, 2> finstate{2, 4};

      /* Evolve the initial state density matrix */
      marray<cmplx, 3> finrho{2, 4, 4};

      for(uint i_nu = 0; i_nu<2; ++i_nu){
        /* First, we start with rho_{44} */
        if(prm.g > 1.0e-7)
          finrho[i_nu][3][3] = getIniRho(E, i_nu, 3, 3) *
            exp(-prm.m_nu[3] * totalGamma * L / E);
        else
          finrho[i_nu][3][3] = getIniRho(E, i_nu, 3, 3);

        /* Then, rho_{4k} */
        for(int j=0; j<3; ++j){
          if(prm.g > 1.0e-7)
            finrho[i_nu][3][j] = getIniRho(E, i_nu, 3, j) *
              exp(-prm.m_nu[3] * totalGamma * L /(2*E) ) *
              cexp(_I*-prm.dmq41 * L /(2.*E));// JK include GLoBES-style low-pass filter
          else
            finrho[i_nu][3][j] = getIniRho(E, i_nu, 3, j) *
              exp(-0.5 * SQR(filter_value) / SQR(E) * SQR(L*0.5*prm.dmq41/E)) * 
              cexp(_I*-prm.dmq41 * L /(2.*E));
        }

        /* And, finally, the 3x3 subblock */
        for(int i=0; i<3; ++i)
          for(int j=0; j<=i; ++j){
            if(prm.g > 1.0e-7)
              finrho[i_nu][i][j] = getIniRho(E, i_nu, i, j)
                + getIntegral(E, L, i_nu)*prm.getProjector(i, j, i_nu);
            else
              finrho[i_nu][i][j] = getIniRho(E, i_nu, i, j);
          }
        for(int i=0; i<3; ++i){ // JK
          for(int j=0; j<=i; ++j){ // JK - implement oscillations and GLoBES-style low-pass filter
            double dmsq = SQR(prm.m_nu[i]) - SQR(prm.m_nu[j]);
            finrho[i_nu][i][j] *= exp(-0.5 * SQR(filter_value) / SQR(E)
                                      * SQR(L*0.5*dmsq/E)) * cexp(_I*-dmsq * L /(2.*E));
          }
        }
          
        /* We can now convert to the flavour basis */
        for(uint i_fla=0; i_fla<4; ++i_fla){
          finstate[i_nu][i_fla] = 0;
          for(uint i=0; i<4; ++i){
            for(uint j=0; j<i; ++j)
              finstate[i_nu][i_fla] += (i_nu==0)? 2*creal( prm.U[i_fla][i]*conj(prm.U[i_fla][j]) *finrho[i_nu][i][j] ):
                2*creal( conj(prm.U[i_fla][i])*prm.U[i_fla][j] *finrho[i_nu][i][j] );
            
            finstate[i_nu][i_fla] += SQR(cabs(prm.U[i_fla][i]))*creal(finrho[i_nu][i][i]);
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

    inline double getInitialTrace(double E){
      marray<double, 2> flx = getInitialFlux(E);
      double res = 0;
      for(uint i_nu = 0; i_nu<2; ++i_nu)
        for(uint i_fla = 0; i_fla<4; ++i_fla)
          res += flx[i_nu][i_fla];
      
      return res;
    }

    inline double getFinalTrace(double E, double L){
      marray<double, 2> flx = getFinalFlux(E, L);
      double res = 0;
      for(uint i_nu = 0; i_nu<2; ++i_nu)
        for(uint i_fla = 0; i_fla<4; ++i_fla)
          res += flx[i_nu][i_fla];
      
      return res;
    }

    /* Destructor */
    ~regProb(){
      if (E_range)     { delete[] E_range;      E_range  = NULL;    }
      if (inistate)    { delete[] inistate;     inistate = NULL;    }
      if (inv_E_steps) { delete[] inv_E_steps;  inv_E_steps = NULL; }
    }
  };

  // Simple utility giving a 1-D array of equally spaced points
  inline marray<double,1> linspace(double Emin,double Emax,unsigned int div){ //JK
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

#endif // ifdef NU_USE_NUSQUIDS

#endif
