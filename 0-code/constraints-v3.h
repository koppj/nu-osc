/*
 * *** COSMOLOGY [1704.00626] ***
 *    I assume that ms and GX are uncorrelated. In that case, the joint PDF factorises and
 *    integrating over ms is identical to taking ms=0.
 *    In addition, for a Gaussian posterior
 *        P(GX) = N * exp[-(GX/GF)^2 / 2sigma^2] ;  N = 2/sqrt(2pi sigma^2)
 *    and GX between 0 and infinity, the Bayesian 95% CL is at
 *        GX/GF = 1.95996 sigma
 *    which means
 *        sigma = 2.0e10/1.95996
 *
 * *** BETA DECAY ***
 *    I assume the chisq for each experiment to be Gaussian, centered at 0
 *     
 */

// JK: I have modified the last sampling point for Re-187 from 1000.0 to 1000.00001
//     to allow for meaningful results at m4 = 1 keV.

#ifndef H_CONSTRAINTS_
#define H_CONSTRAINTS_

#ifdef NU_USE_NUSQUIDS // JK

//JK #include "prob_v15_scalar.h"
#include "osc_decay/osc_decay.h" // JK
#include <gsl/gsl_spline.h>
#include <fstream>

using namespace regeneration;

namespace ns_sbl { // JK

class Constraints{
private:
  /* Utilities for interpolation */
  gsl_interp_accel *acc_Re187, *acc_1210, *acc_1307, *acc_1703, *acc_Hiddemann, *acc_Ni63;
  gsl_spline *spl_Re187, *spl_1210, *spl_1307, *spl_1703, *spl_Hiddemann, *spl_Ni63;
  
  static constexpr int N_Re187 = 17, N_1210 = 25, N_1307 = 63, N_1703 = 8, N_Hiddemann = 15, N_Ni63 = 300; // Number of data points
  double m_Re187[N_Re187], m_1210[N_1210], m_1307[N_1307], m_1703[N_1703], m_Hiddemann[N_Hiddemann], m_Ni63[N_Ni63]; // m4
  double ue4sq_Re187[N_Re187], ue4sq_1210[N_1210], ue4sq_1307[N_1307], ue4sq_1703[N_1703], ue4sq_Hiddemann[N_Hiddemann], ue4sq_Ni63[N_Ni63]; // Ue4sq
  static constexpr double dChq_Re187 = 3.84146 , dChq_1210 = 2.70554, dChq_1307 = 3.84146, dChq_1703 = 3.84146, dChq_Hiddemann = 3.84146, dChq_Ni63 = 3.84146; // dChq at the table value

  /* Fermi constant [eV^-2]*/
  const double GF = 1.166379e-23;
  
public:
  double getChisq_cosmology(Param &prm){
    prm.loadU(); // Make sure that the mixing matrix is properly loaded
    double Us1 = cabs(prm.U[3][0]);
    double GX = 0.125 * SQR(prm.g) * SQR(SQR(Us1)) / SQR(prm.m_A); // JK
    return SQR(GX/GF) / SQR(2e10/1.95996);
  }

  double getChisq_beta(Param &prm){
    prm.loadU(); // Make sure that the mixing matrix is properly loaded
    double m4 = prm.m_nu[3];
    double Ue4sq = SQR(cabs(prm.U[0][3]));

    double chisq_Re187 = ((m4 < m_Re187[0]) || (m4 > m_Re187[N_Re187-1])) ? 0. :
      dChq_Re187 * SQR(Ue4sq/gsl_spline_eval(spl_Re187, m4, acc_Re187));
    double chisq_1210 = ((m4 < m_1210[0]) || (m4 > m_1210[N_1210-1])) ? 0. :
      dChq_1210 * SQR(Ue4sq/gsl_spline_eval(spl_1210, m4, acc_1210));
    double chisq_1307 = ((m4 < m_1307[0]) || (m4 > m_1307[N_1307-1])) ? 0. :
      dChq_1307 * SQR(Ue4sq/gsl_spline_eval(spl_1307, m4, acc_1307));
    double chisq_1703 = ((m4 < m_1703[0]) || (m4 > m_1703[N_1703-1])) ? 0. :
      dChq_1703 * SQR(Ue4sq/gsl_spline_eval(spl_1703, m4, acc_1703));
    double chisq_Hiddemann = ((m4 < m_Hiddemann[0]) || (m4 > m_Hiddemann[N_Hiddemann-1])) ? 0. :
      dChq_Hiddemann * SQR(Ue4sq/gsl_spline_eval(spl_Hiddemann, m4, acc_Hiddemann));
    double chisq_Ni63 = ((m4 < m_Ni63[0]) || (m4 > m_Ni63[N_Ni63-1])) ? 0. :
      dChq_Ni63 * SQR(Ue4sq/gsl_spline_eval(spl_Ni63, m4, acc_Ni63));

    // JK enforce mass within the range for which we have the limits implemented
    if (chisq_Re187 > 0.  ||  chisq_1210 > 0.  ||  chisq_1307 > 0. || chisq_1703 > 0. || chisq_Hiddemann > 0. || chisq_Ni63 > 0.)
      return chisq_Re187 + chisq_1210 + chisq_1307 + chisq_1703 + chisq_Hiddemann + chisq_Ni63;
    else if ((m4 < m_1210[0]) && (m4 < m_1307[0])) // Return 0 for lower masses
      return 0;
    else
      return 1e17;
  }
  
  double getChisq(Param &prm){
    if (prm.g > 3.5449077) // IE Perturbativity
      return 1e17;
    return getChisq_cosmology(prm) + getChisq_beta(prm);
  }

  /* Constructor */
  Constraints(){
    /* Read the data files */
    std::ifstream ifile;
    double m4, ue4sq;
    
    ifile.open("glb/ext-constraints/Re-187_unraw.dat");
    for(int i=0; i<N_Re187; ++i){
      ifile>>m4>>ue4sq;      
      m_Re187[i] = m4;
      ue4sq_Re187[i] = ue4sq;
    }
    ifile.close();

    ifile.open("glb/ext-constraints/1210.4194_unraw.dat");
    for(int i=0; i<N_1210; ++i){
      ifile>>m4>>ue4sq;            
      m_1210[i] = m4;
      ue4sq_1210[i] = ue4sq;
    }
    ifile.close();

    ifile.open("glb/ext-constraints/1307.5687_unraw.dat");
    for(int i=0; i<N_1307; ++i){
      ifile>>m4>>ue4sq;            
      m_1307[i] = m4;
      ue4sq_1307[i] = ue4sq;
    }
    ifile.close();

    ifile.open("glb/ext-constraints/1703.10779_unraw.dat");
    for(int i=0; i<N_1703; ++i){
      ifile>>m4>>ue4sq;            
      m_1703[i] = m4;
      ue4sq_1703[i] = ue4sq;
    }
    ifile.close();

    ifile.open("glb/ext-constraints/Hiddemann.dat");
    for(int i=0; i<N_Hiddemann; ++i){
      ifile>>m4>>ue4sq;            
      m_Hiddemann[i] = m4;
      ue4sq_Hiddemann[i] = ue4sq;
    }
    ifile.close();

    ifile.open("glb/ext-constraints/Ni-63.dat");
    for(int i=0; i<N_Ni63; ++i){
      ifile>>m4>>ue4sq;            
      m_Ni63[i] = m4;
      ue4sq_Ni63[i] = ue4sq;
    }
    ifile.close();

    /* Load the interpolators */
    acc_Re187 = gsl_interp_accel_alloc();
    spl_Re187 = gsl_spline_alloc(gsl_interp_linear, N_Re187);
    gsl_spline_init(spl_Re187, m_Re187, ue4sq_Re187, N_Re187);

    acc_1210 = gsl_interp_accel_alloc();
    spl_1210 = gsl_spline_alloc(gsl_interp_linear, N_1210);
    gsl_spline_init(spl_1210, m_1210, ue4sq_1210, N_1210);
    
    acc_1307 = gsl_interp_accel_alloc();
    spl_1307 = gsl_spline_alloc(gsl_interp_linear, N_1307);
    gsl_spline_init(spl_1307, m_1307, ue4sq_1307, N_1307);

    acc_1703 = gsl_interp_accel_alloc();
    spl_1703 = gsl_spline_alloc(gsl_interp_linear, N_1703);
    gsl_spline_init(spl_1703, m_1703, ue4sq_1703, N_1703);

    acc_Hiddemann = gsl_interp_accel_alloc();
    spl_Hiddemann = gsl_spline_alloc(gsl_interp_linear, N_Hiddemann);
    gsl_spline_init(spl_Hiddemann, m_Hiddemann, ue4sq_Hiddemann, N_Hiddemann);

    acc_Ni63 = gsl_interp_accel_alloc();
    spl_Ni63 = gsl_spline_alloc(gsl_interp_linear, N_Ni63);
    gsl_spline_init(spl_Ni63, m_Ni63, ue4sq_Ni63, N_Ni63); 
  }

  /* Destructor */
  ~Constraints(){
    gsl_spline_free(spl_Re187);
    gsl_interp_accel_free(acc_Re187);
    gsl_spline_free(spl_1210);
    gsl_interp_accel_free(acc_1210);
    gsl_spline_free(spl_1307);
    gsl_interp_accel_free(acc_1307);
    gsl_spline_free(spl_1703);
    gsl_interp_accel_free(acc_1703);
    gsl_spline_free(spl_Hiddemann);
    gsl_interp_accel_free(acc_Hiddemann);
    gsl_spline_free(spl_Ni63);
    gsl_interp_accel_free(acc_Ni63);
  }

  /* No copy (for simplicity) */
  Constraints(const Constraints&) = delete;
  Constraints &operator=(const Constraints&) = delete;
};

} // end namespace ns_sbl (JK)

#endif

#endif
