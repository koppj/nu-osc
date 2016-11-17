#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <sys/stat.h>
#include <math.h>
#include <assert.h>
#include <complex>

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "THStack.h"
#include "TMatrixD.h"
#include "TFile.h"
#include "TDecompSVD.h"
#include "TMarker.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TMarker.h"
#include "TGraph.h"
#include "TGraphSmooth.h"
#include "TMath.h"


using namespace std;


const double k1267 = 1.26693276;
const double kKmUnits = 1000.;
double chisqCC;
double chisqNC;
double Penalty_NormCC;
double Penalty_NormNC;
double Penalty_dm232;
double totalChi2_CC;
double totalChi2_NC;
double totalChi2;
  
// Near MC
TH2D *NDNC_TrueNC, *NDNC_NuMu, *NDNC_BeamNue, *NDNC_AppNue, *NDNC_AppNuTau;
TH2D *NDCC_TrueNC, *NDCC_NuMu, *NDCC_BeamNue, *NDCC_AppNue, *NDCC_AppNuTau;

// Far MC
TH2D *FDNC_TrueNC, *FDNC_NuMu, *FDNC_BeamNue, *FDNC_AppNue, *FDNC_AppNuTau;
TH2D *FDCC_TrueNC, *FDCC_NuMu, *FDCC_BeamNue, *FDCC_AppNue, *FDCC_AppNuTau;

// Data (as well as fake) 
TH1D *FD_dataNC, *FD_dataCC, *ND_dataNC, *ND_dataCC;
TH1D *NDUnOscNC_MC, *NDUnOscCC_MC;

// MC (Oscillated)
TH1D *FDOscCC_MC, *NDOscCC_MC, *FDOscNC_MC, *NDOscNC_MC;

TMatrixD* CoVarCCinvert;
TMatrixD* CoVarNCinvert;
  
//---------------------------------------------------------------------------------
struct params
{
  double Dm232;
  double Dm221;
  double th23;
  double th12;
  double th13;
  double deltaCP;
  double Dm241;
  double th24;
  double th34;
  double th14;
  double delta24;
};
//---------------------------------------------------------------------------------
void Zombie(TFile* f);
Double_t PenaltyTermDm232(Double_t dm232);
Double_t NDNormPenaltyTerm(TH1D* data, TH1D* pred, TH1D* unosc);
Double_t ChiSqFunction(TH1D* rPred, TH1D* rData, TMatrixD* CoVarInvert);
Double_t ComparePredWithData(TH1D* RpredCC, 
			     TH1D* ND_predCC, 
			     TH1D* ND_UnOscCC, 
			     TH1D* RdataCC, 
			     TH1D* ND_dataCC, 
			     TMatrixD* CoVarCCinvert, 
			     TH1D* RpredNC, 
			     TH1D* ND_predNC, 
			     TH1D* ND_UnOscNC, 
			     TH1D* RdataNC, 
			     TH1D* ND_dataNC, 
			     TMatrixD* CoVarNCinvert,
			     Double_t Dm2);
TH1D* CreateTotalSpectrum(params my_pars,
			  TH2D* TrueNC,
			  TH2D* NuMu,
			  TH2D* BeamNue,
			  TH2D* AppNue,
			  TH2D* AppNuTau,
			  double baseline);
double FourFlavourNuMuToNuSProbability
(const double energy, 
 double dm232, const double theta23, double dm221, 
 double dm243, const double theta12, 
 const double theta13, const double theta14,
 const double theta24, const double theta34,
 double delta1, double delta2, double delta3,
 const double baseline);
double FourFlavourDisappearanceWeight
(const double energy, 
 double dm232, const double theta23, double dm221, 
 double dm243, const double theta12, 
 const double theta13, const double theta14,
 const double theta24, const double theta34,
 double delta1, double delta2, double delta3,
 const double baseline);
double FourFlavourNuESurvivalProbability
(const double energy, 
 double dm232, const double theta23, double dm221, 
 double dm243, const double theta12, 
 const double theta13, const double theta14,
 const double theta24, const double theta34,
 double delta1, double delta2, double delta3,
 const double baseline);
double FourFlavourNuMuToNuEProbability
(const double energy, 
 double dm232, const double theta23, double dm221, 
 double dm243, const double theta12, 
 const double theta13, const double theta14,
 const double theta24, const double theta34,
 double delta1, double delta2, double delta3,
 const double baseline);
double FourFlavourNuMuToNuTauProbability
(const double energy, 
 double dm232, const double theta23, double dm221, 
 double dm243, const double theta12, 
 const double theta13, const double theta14,
 const double theta24, const double theta34,
 double delta1, double delta2, double delta3,
 const double baseline);
TH1D* CreateSpectrumComponent(params my_pars, TString OscType, TH2D* oscDummy, Double_t baseline);
void PrintParStatus(params my_pars);
//---------------------------------------------------------------------------------
void Zombie(TFile* f){

  if(f->IsZombie() && (!f->IsOpen())){
    std::cout << "File " << f->GetName() << " failed to open." << std::endl;
    assert(false);
  }
  else{
    std::cout << "File " << f->GetName() << " opened successfully" << std::endl;
  }
}
//---------------------------------------------------------------------------------
void PrintParStatus(params my_pars)
{
  std::cout << "" << std::endl;
  std::cout << "/==========Parameter Status==========/" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "Dm232        =      " << my_pars.Dm232 << std::endl;
  std::cout << "Dm221        =      " << my_pars.Dm221 << std::endl;
  std::cout << "Theta23      =      " << my_pars.th23 << std::endl;
  std::cout << "Theta12      =      " << my_pars.th12 << std::endl;
  std::cout << "Theta13      =      " << my_pars.th13 << std::endl;
  std::cout << "DeltaCP      =      " << my_pars.deltaCP << std::endl;
  std::cout << "Dm241        =      " << my_pars.Dm241 << std::endl;
  std::cout << "Theta24      =      " << my_pars.th24 << std::endl;
  std::cout << "Theta34      =      " << my_pars.th34 << std::endl;
  std::cout << "Theta14      =      " << my_pars.th14 << std::endl;
  std::cout << "Delta24      =      " << my_pars.delta24 << std::endl;
  std::cout << "" << std::endl;
  std::cout << "/====================================/" << std::endl;
  std::cout << "" << std::endl;
}
//---------------------------------------------------------------------------------
Double_t PenaltyTermDm232(Double_t dm232)
{
  Double_t dm232_pen = 0.0;
  dm232_pen = TMath::Power( (TMath::Abs(dm232) - 0.0025) , 2); // numerator
  dm232_pen /= TMath::Power( 0.0005, 2);

  return dm232_pen;
}
//---------------------------------------------------------------------------------
Double_t NDNormPenaltyTerm(TH1D* data, TH1D* pred, TH1D* unosc)
{
  Int_t NumberOfBins = data->GetNbinsX();
  
  Double_t UnOsc = unosc->Integral(1, NumberOfBins );
  Double_t TotalNDData = data->Integral(1, NumberOfBins );
  Double_t TotalOscNDMC = pred->Integral(1, NumberOfBins );

  // This is a 50% uncertainty 
  Double_t sigma = 0.5*UnOsc;
  Double_t Norm = (TotalNDData - TotalOscNDMC) * (TotalNDData - TotalOscNDMC);

  Norm = Norm / (sigma * sigma);

  return Norm;
}
//---------------------------------------------------------------------------------
Double_t ChiSqFunction(TH1D* rPred, TH1D* rData, TMatrixD* CoVarInvert)
{
  if(!(rPred->GetNbinsX() == rData->GetNbinsX())){ 
    std::cout << "Incorrect Binning agreement. Asserting" << std::endl;
    assert(false);
  }

  Int_t NumberOfBins = rPred->GetNbinsX();

  TVectorD Difference(NumberOfBins);

  for(Int_t i=1; i<=NumberOfBins; ++i){
    Difference(i-1) = (rData->GetBinContent(i) - rPred->GetBinContent(i));
  }

  TVectorD temp = Difference;
  temp *= (*CoVarInvert);

  Double_t TotalChiSq = temp*Difference;

  return TotalChiSq;

}
//---------------------------------------------------------------------------------
//Double_t ComparePredWithData(const NuMMParameters &pars)
Double_t ComparePredWithData(TH1D* RpredCC, 
			     TH1D* ND_predCC, 
			     TH1D* ND_UnOscCC, 
			     TH1D* RdataCC, 
			     TH1D* ND_dataCC, 
			     TMatrixD* CoVarCCinvert, 
			     TH1D* RpredNC, 
			     TH1D* ND_predNC, 
			     TH1D* ND_UnOscNC, 
			     TH1D* RdataNC, 
			     TH1D* ND_dataNC, 
			     TMatrixD* CoVarNCinvert,
			     Double_t Dm2
			    )
{
  chisqCC    = ChiSqFunction(RpredCC, RdataCC, CoVarCCinvert);
  chisqNC    = ChiSqFunction(RpredNC, RdataNC, CoVarNCinvert);
  
  Penalty_NormCC  = NDNormPenaltyTerm(ND_dataCC,ND_predCC,ND_UnOscCC);
  Penalty_NormNC  = NDNormPenaltyTerm(ND_dataNC,ND_predNC,ND_UnOscNC);
  Penalty_dm232   = PenaltyTermDm232(Dm2);

  // CC and NC parts
  totalChi2_CC = chisqCC + Penalty_NormCC;
  totalChi2_NC = chisqNC + Penalty_NormNC;

  totalChi2 = totalChi2_CC + totalChi2_NC + Penalty_dm232;

  std::cout << "" << std::endl;
  std::cout << "/=======Chi2 Calculator Output=======/" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "NC Chi2         =     " << chisqNC << std::endl;
  std::cout << "NC ND Norm Pen  =     " << Penalty_NormNC << std::endl;
  std::cout << "CC Chi2         =     " << chisqCC << std::endl;
  std::cout << "CC ND Norm Pen  =     " << Penalty_NormCC << std::endl;
  std::cout << "Penalty Dm232   =     " << Penalty_dm232 << std::endl;
  std::cout << "Total Chi2 CC   =     " << totalChi2_CC << std::endl;
  std::cout << "Total Chi2 NC   =     " << totalChi2_NC << std::endl;
  std::cout << "Total Chi2      =     " << totalChi2 << std::endl;
  std::cout << "" << std::endl;
  std::cout << "/====================================/" << std::endl;
  std::cout << "" << std::endl;

  return totalChi2;
}
//---------------------------------------------------------------------------------
TH1D* CreateTotalSpectrum(params my_pars,
			  TH2D* TrueNC,
			  TH2D* NuMu,
			  TH2D* BeamNue,
			  TH2D* AppNue,
			  TH2D* AppNuTau,
			  double baseline
			 )
{
  TH1D* vtruenc   = (TH1D*)CreateSpectrumComponent(my_pars, "TrueNC",   TrueNC,   baseline);
  TH1D* vnumu     = (TH1D*)CreateSpectrumComponent(my_pars, "NuMu",     NuMu,     baseline);
  TH1D* vbeamnue  = (TH1D*)CreateSpectrumComponent(my_pars, "BeamNue",  BeamNue,  baseline);
  TH1D* vappnue   = (TH1D*)CreateSpectrumComponent(my_pars, "AppNue",   AppNue,   baseline);
  TH1D* vappnutau = (TH1D*)CreateSpectrumComponent(my_pars, "AppNuTau", AppNuTau, baseline);

  TH1D* hTotal = new TH1D(*vtruenc);
  hTotal->Add(vnumu);
  hTotal->Add(vbeamnue);
  hTotal->Add(vappnue);
  hTotal->Add(vappnutau);

  return hTotal;
}
//---------------------------------------------------------------------------------
double FourFlavourNuMuToNuSProbability
(const double energy, 
 double dm232, const double theta23, double dm221, 
 double dm243, const double theta12, 
 const double theta13, const double theta14,
 const double theta24, const double theta34,
 double delta1, double delta2, double delta3,
 const double baseline)
{ 

  // Calculate other mass splittings
  const double dm231 = dm221 + dm232;
  const double dm241 = dm231 + dm243;

  const double c12 = cos(theta12); const double s12 = sin(theta12);
  const double c13 = cos(theta13); const double s13 = sin(theta13);
  const double c14 = cos(theta14); const double s14 = sin(theta14);
  const double c23 = cos(theta23); const double s23 = sin(theta23);
  const double c24 = cos(theta24); const double s24 = sin(theta24);
  const double c34 = cos(theta34); const double s34 = sin(theta34);

  complex<double> expNegCP13 = complex<double>(cos(delta1), -sin(delta1));
  complex<double> expNegCP14 = complex<double>(cos(delta2), -sin(delta2));
  complex<double> expNegCP24 = complex<double>(cos(delta3), -sin(delta3));

  complex<double> Umu2  =  c12 * c23 * c24
                        -  c24 * s12 * s13 * s23 * conj(expNegCP13)
                        -  c13 * s12 * s14 * s24 * expNegCP24 * conj(expNegCP14);

  complex<double> Umu3  =  c13 * c24 * s23
                        -  s13 * s14 * s24 * expNegCP13 * expNegCP24 * conj(expNegCP14);
  
  complex<double> Umu4  =  c14 * s24 * expNegCP24;
  

  complex<double> Us2   =  -c13 * c24 * c34 * s12 * s14 * conj(expNegCP14)
                           -c12 * c23 * c34 * s24 * conj(expNegCP24)
                           +c34 * s12 * s13 * s23 * s24 * conj(expNegCP13 * expNegCP24)
                           +c23 * s12 * s13 * s34 * conj(expNegCP13)
                           +c12 * s23 * s34;
  
  complex<double> Us3   =  -c24 * c34 * s13 * s14 * expNegCP13 * conj(expNegCP14)
                           -c13 * c34 * s23 * s24 * conj(expNegCP24)
                           -c13 * c23 * s34;
  
  complex<double> Us4   =  c14 * c24 * c34;
  
  
  complex<double> i(0.0, 1.0);
  
  double DeltaM21 = (k1267 * dm221 * baseline / kKmUnits) / energy;
  double DeltaM31 = (k1267 * dm231 * baseline / kKmUnits) / energy;
  double DeltaM41 = (k1267 * dm241 * baseline / kKmUnits) / energy;
  
  complex<double> expDeltaM21 = complex<double>(cos(DeltaM21), -sin(DeltaM21));
  complex<double> expDeltaM31 = complex<double>(cos(DeltaM31), -sin(DeltaM31));
  complex<double> expDeltaM41 = complex<double>(cos(DeltaM41), -sin(DeltaM41));
  
  double oscProb = norm(-2.0 * i * conj(Umu2) * Us2 * sin(DeltaM21) * expDeltaM21
                        -2.0 * i * conj(Umu3) * Us3 * sin(DeltaM31) * expDeltaM31
			-2.0 * i * conj(Umu4) * Us4 * sin(DeltaM41) * expDeltaM41);

  return oscProb;
}
//---------------------------------------------------------------------------------
double FourFlavourDisappearanceWeight
(const double energy, 
 double dm232, const double theta23, double dm221, 
 double dm243, const double theta12, 
 const double theta13, const double theta14,
 const double theta24, const double theta34,
 double delta1, double delta2, double delta3,
 const double baseline)
{

  // Calculate other mass splittings
  const double dm231 = dm221 + dm232;
  const double dm241 = dm231 + dm243;

  const double c12 = cos(theta12); const double s12 = sin(theta12);
  const double c13 = cos(theta13); const double s13 = sin(theta13);
  const double c14 = cos(theta14); const double s14 = sin(theta14);
  const double c23 = cos(theta23); const double s23 = sin(theta23);
  const double c24 = cos(theta24); const double s24 = sin(theta24);

  complex<double> expNegCP13 = complex<double>(cos(delta1), -sin(delta1));
  complex<double> expNegCP14 = complex<double>(cos(delta2), -sin(delta2));
  complex<double> expNegCP24 = complex<double>(cos(delta3), -sin(delta3));

  complex<double> Umu2  =  c12 * c23 * c24
                        -  c24 * s12 * s13 * s23 * conj(expNegCP13)
                        -  c13 * s12 * s14 * s24 * expNegCP24 * conj(expNegCP14);

  complex<double> Umu3  =  c13 * c24 * s23
                        -  s13 * s14 * s24 * expNegCP13 * expNegCP24 * conj(expNegCP14);
  
  complex<double> Umu4  =  c14 * s24 * expNegCP24;
  
  complex<double> i(0.0, 1.0);
  
  double DeltaM21 = (k1267 * dm221 * baseline / kKmUnits) / energy;
  double DeltaM31 = (k1267 * dm231 * baseline / kKmUnits) / energy;
  double DeltaM41 = (k1267 * dm241 * baseline / kKmUnits) / energy;
  
  complex<double> expDeltaM21 = complex<double>(cos(DeltaM21), -sin(DeltaM21));
  complex<double> expDeltaM31 = complex<double>(cos(DeltaM31), -sin(DeltaM31));
  complex<double> expDeltaM41 = complex<double>(cos(DeltaM41), -sin(DeltaM41));
  
  double oscProb  = norm(1.0 
			 - 2.0 * i * conj(Umu2) * Umu2 * sin(DeltaM21) * expDeltaM21 
			 - 2.0 * i * conj(Umu3) * Umu3 * sin(DeltaM31) * expDeltaM31 
			 - 2.0 * i * conj(Umu4) * Umu4 * sin(DeltaM41) * expDeltaM41);
  
  return oscProb;
}
//---------------------------------------------------------------------------------
double FourFlavourNuESurvivalProbability
(const double energy, 
 double dm232, const double theta23, double dm221, 
 double dm243, const double theta12, 
 const double theta13, const double theta14,
 const double theta24, const double theta34,
 double delta1, double delta2, double delta3,
 const double baseline)
{

  // Calculate other mass splittings
  const double dm231 = dm221 + dm232;
  const double dm241 = dm231 + dm243;

  const double s12 = sin(theta12);
  const double c13 = cos(theta13); const double s13 = sin(theta13);
  const double c14 = cos(theta14); const double s14 = sin(theta14);

  complex<double> expNegCP13 = complex<double>(cos(delta1), -sin(delta1));
  complex<double> expNegCP14 = complex<double>(cos(delta2), -sin(delta2));

  complex<double> Ue2   =  c13 * c14 * s12;
  complex<double> Ue3   =  c14 * s13 * expNegCP13;
  complex<double> Ue4   =  s14 * expNegCP14;
  
  complex<double> i(0.0, 1.0);
  
  double DeltaM21 = (k1267 * dm221 * baseline / kKmUnits) / energy;
  double DeltaM31 = (k1267 * dm231 * baseline / kKmUnits) / energy;
  double DeltaM41 = (k1267 * dm241 * baseline / kKmUnits) / energy;
  
  complex<double> expDeltaM21 = complex<double>(cos(DeltaM21), -sin(DeltaM21));
  complex<double> expDeltaM31 = complex<double>(cos(DeltaM31), -sin(DeltaM31));
  complex<double> expDeltaM41 = complex<double>(cos(DeltaM41), -sin(DeltaM41));
  
  double oscProb = norm(1.0 
			- 2.0 * i * conj(Ue2) * Ue2 * sin(DeltaM21) * expDeltaM21 
			- 2.0 * i * conj(Ue3) * Ue3 * sin(DeltaM31) * expDeltaM31 
			- 2.0 * i * conj(Ue4) * Ue4 * sin(DeltaM41) * expDeltaM41);

  return oscProb;
}
//---------------------------------------------------------------------------------
double FourFlavourNuMuToNuEProbability
(const double energy, 
 double dm232, const double theta23, double dm221, 
 double dm243, const double theta12, 
 const double theta13, const double theta14,
 const double theta24, const double theta34,
 double delta1, double delta2, double delta3,
 const double baseline)
{

  // Calculate other mass splittings
  const double dm231 = dm221 + dm232;
  const double dm241 = dm231 + dm243;

  const double c12 = cos(theta12); const double s12 = sin(theta12);
  const double c13 = cos(theta13); const double s13 = sin(theta13);
  const double c14 = cos(theta14); const double s14 = sin(theta14);
  const double c23 = cos(theta23); const double s23 = sin(theta23);
  const double c24 = cos(theta24); const double s24 = sin(theta24);

  complex<double> expNegCP13 = complex<double>(cos(delta1), -sin(delta1));
  complex<double> expNegCP14 = complex<double>(cos(delta2), -sin(delta2));
  complex<double> expNegCP24 = complex<double>(cos(delta3), -sin(delta3));
  
  complex<double> Umu2  =  c12 * c23 * c24
                        -  c24 * s12 * s13 * s23 * conj(expNegCP13)
                        -  c13 * s12 * s14 * s24 * expNegCP24 * conj(expNegCP14);

  complex<double> Umu3  =  c13 * c24 * s23
                        -  s13 * s14 * s24 * expNegCP13 * expNegCP24 * conj(expNegCP14);
  
  complex<double> Umu4  =  c14 * s24 * expNegCP24;


  complex<double> Ue2   =  c13 * c14 * s12;
  complex<double> Ue3   =  c14 * s13 * expNegCP13;
  complex<double> Ue4   =  s14 * expNegCP14;
  
  complex<double> i(0.0, 1.0);
  
  double DeltaM21 = (k1267 * dm221 * baseline / kKmUnits) / energy;
  double DeltaM31 = (k1267 * dm231 * baseline / kKmUnits) / energy;
  double DeltaM41 = (k1267 * dm241 * baseline / kKmUnits) / energy;
  
  complex<double> expDeltaM21 = complex<double>(cos(DeltaM21), -sin(DeltaM21));
  complex<double> expDeltaM31 = complex<double>(cos(DeltaM31), -sin(DeltaM31));
  complex<double> expDeltaM41 = complex<double>(cos(DeltaM41), -sin(DeltaM41));
  
  double oscProb = norm(-2.0 * i * conj(Umu2) * Ue2 * sin(DeltaM21) * expDeltaM21       
			-2.0 * i * conj(Umu3) * Ue3 * sin(DeltaM31) * expDeltaM31
			-2.0 * i * conj(Umu4) * Ue4 * sin(DeltaM41) * expDeltaM41);

  return oscProb;
}
//---------------------------------------------------------------------------------
double FourFlavourNuMuToNuTauProbability
(const double energy, 
 double dm232, const double theta23, double dm221, 
 double dm243, const double theta12, 
 const double theta13, const double theta14,
 const double theta24, const double theta34,
 double delta1, double delta2, double delta3,
 const double baseline)
{

  // Calculate other mass splittings
  const double dm231 = dm221 + dm232;
  const double dm241 = dm231 + dm243;

  const double c12 = cos(theta12); const double s12 = sin(theta12);
  const double c13 = cos(theta13); const double s13 = sin(theta13);
  const double c14 = cos(theta14); const double s14 = sin(theta14);
  const double c23 = cos(theta23); const double s23 = sin(theta23);
  const double c24 = cos(theta24); const double s24 = sin(theta24);
  const double c34 = cos(theta34); const double s34 = sin(theta34);

  complex<double> expNegCP13 = complex<double>(cos(delta1), -sin(delta1));
  complex<double> expNegCP14 = complex<double>(cos(delta2), -sin(delta2));
  complex<double> expNegCP24 = complex<double>(cos(delta3), -sin(delta3));

  complex<double> Umu2  =  c12 * c23 * c24
                        -  c24 * s12 * s13 * s23 * conj(expNegCP13)
                        -  c13 * s12 * s14 * s24 * expNegCP24 * conj(expNegCP14);

  complex<double> Umu3  =  c13 * c24 * s23
                        -  s13 * s14 * s24 * expNegCP13 * expNegCP24 * conj(expNegCP14);
  
  complex<double> Umu4  =  c14 * s24 * expNegCP24;


  complex<double> Utau2 =  -c12 * c34 * s23
                           -c23 * c34 * s12 * s13 * conj(expNegCP13)
                           -c13 * c24 * s12 * s14 * s34 * conj(expNegCP14)
                           -c12 * c23 * s24 * s34 * conj(expNegCP24)
                           +s12 * s13 * s23 * s24 * s34 * conj(expNegCP13 * expNegCP24);

  complex<double> Utau3 =  c13 * c23 * c34
                        -  c24 * s13 * s14 * s34 * expNegCP13 * conj(expNegCP14)
                        -  c13 * s23 * s24 * s34 * conj(expNegCP24);

  complex<double> Utau4 =  c14 * c24 * s34;

  complex<double> i(0.0, 1.0);
  
  double DeltaM21 = (k1267 * dm221 * baseline / kKmUnits) / energy;
  double DeltaM31 = (k1267 * dm231 * baseline / kKmUnits) / energy;
  double DeltaM41 = (k1267 * dm241 * baseline / kKmUnits) / energy;
  
  complex<double> expDeltaM21 = complex<double>(cos(DeltaM21), -sin(DeltaM21));
  complex<double> expDeltaM31 = complex<double>(cos(DeltaM31), -sin(DeltaM31));
  complex<double> expDeltaM41 = complex<double>(cos(DeltaM41), -sin(DeltaM41));
  
  double oscProb =  norm(-2.0 * i * conj(Umu2) * Utau2 * sin(DeltaM21) * expDeltaM21     
			 -2.0 * i * conj(Umu3) * Utau3 * sin(DeltaM31) * expDeltaM31
			 -2.0 * i * conj(Umu4) * Utau4 * sin(DeltaM41) * expDeltaM41);

  return oscProb;
}
//---------------------------------------------------------------------------------
TH1D* CreateSpectrumComponent(params my_pars, TString OscType, TH2D* oscDummy, Double_t baseline)
{
  TH1D* bintemplate = oscDummy->ProjectionY();
  bintemplate->Reset();
 
  const double k1267 = 1.26693276;

  // Loop over every true energy bin in the reco vs. true matrices, then loop over every reco energy in that bin                                                 
  // to calculate an oscillation weight for that reco energy based on the true energy. 
  TAxis *Yaxis = oscDummy->GetYaxis();
  TAxis *Xaxis = oscDummy->GetXaxis();

  // Define Dm243 such that its actually Dm241 being altered.
  //41 = 43 + 32 + 21 
  //43 = 41 - 32 - 21
  Double_t dm243 = 0.0;

  dm243 = my_pars.Dm241 - my_pars.Dm232 - my_pars.Dm221;

  for(Int_t x = 1; x <= Xaxis->GetNbins(); x++){
    Double_t OscWeight = 0.0;
    
    if(baseline > 0){
      
      // Default iterations (1 at bin center)
      Int_t n_LoverE = 1;
      Double_t LoverE[5];
      LoverE[0] = Xaxis->GetBinCenter(x);
      
      // This is averaging oscialltions in true energy bins - see Technical Note http://minos-docdb.fnal.gov/cgi-bin/RetrieveFile?docid=10203&version=2
      const Double_t W = Xaxis->GetBinWidth(x);
      const Double_t arg = k1267*dm243*W; // half-period of oscillation
      Double_t sample = W/2/sqrt(3);

      if(arg!=0) sample = TMath::ACos(TMath::Sin(arg)/arg)/arg*W/2;

      n_LoverE = 2;
      Double_t bc = LoverE[0]; // bin center
      LoverE[0] = bc - sample;
      LoverE[1] = bc + sample;

      const Double_t E = 1.0;

      for(int i = 0; i < n_LoverE; i++){
	
	// each Osctype has a different probability function
	if(OscType == "TrueNC"){

	  OscWeight += FourFlavourNuMuToNuSProbability( E, 
  							my_pars.Dm232, 
							my_pars.th23, 
							my_pars.Dm221, 
							dm243, 
							my_pars.th12, 
							my_pars.th13, 
							my_pars.th14, 
							my_pars.th24, 
							my_pars.th34, 
							my_pars.deltaCP,  
							0, 
							my_pars.delta24,
							LoverE[i]*kKmUnits);
	}
	if(OscType == "NuMu"){

	  OscWeight += FourFlavourDisappearanceWeight( E, 
  							my_pars.Dm232, 
							my_pars.th23, 
							my_pars.Dm221, 
							dm243, 
							my_pars.th12, 
							my_pars.th13, 
							my_pars.th14, 
							my_pars.th24, 
							my_pars.th34, 
							my_pars.deltaCP,  
							0, 
							my_pars.delta24,
							LoverE[i]*kKmUnits);
	}
	if(OscType == "BeamNue"){
	  
	  OscWeight += FourFlavourNuESurvivalProbability( E, 
  							my_pars.Dm232, 
							my_pars.th23, 
							my_pars.Dm221, 
							dm243, 
							my_pars.th12, 
							my_pars.th13, 
							my_pars.th14, 
							my_pars.th24, 
							my_pars.th34, 
							my_pars.deltaCP,  
							0, 
							my_pars.delta24,
							LoverE[i]*kKmUnits);
	}
	if(OscType == "AppNue"){
	  
	  OscWeight += FourFlavourNuMuToNuEProbability( E, 
  							my_pars.Dm232, 
							my_pars.th23, 
							my_pars.Dm221, 
							dm243, 
							my_pars.th12, 
							my_pars.th13, 
							my_pars.th14, 
							my_pars.th24, 
							my_pars.th34, 
							my_pars.deltaCP,  
							0, 
							my_pars.delta24,
							LoverE[i]*kKmUnits);
	}
	if(OscType == "AppNuTau"){

	  OscWeight += FourFlavourNuMuToNuTauProbability( E, 
  							my_pars.Dm232, 
							my_pars.th23, 
							my_pars.Dm221, 
							dm243, 
							my_pars.th12, 
							my_pars.th13, 
							my_pars.th14, 
							my_pars.th24, 
							my_pars.th34, 
							my_pars.deltaCP,  
							0, 
							my_pars.delta24,
							LoverE[i]*kKmUnits);
	}
      }
      // Now average this
      OscWeight /= n_LoverE;
    }
    else // if baseline < 0
      {
	
        if(OscType == "TrueNC")   OscWeight = 0.0;
        if(OscType == "NuMu")     OscWeight = 1.0;
        if(OscType == "BeamNue")  OscWeight = 1.0;
        if(OscType == "AppNue")   OscWeight = 0.0;
        if(OscType == "AppNuTau") OscWeight = 0.0;
      }

    // using the oscillation weight, fill a 1d histogram for each type of event with the oscillated reco energy 
    for(Int_t y = 1; y <= Yaxis->GetNbins(); y++){
      
      Double_t sumWeights = 0;
      
      if(OscType == "TrueNC"){
	sumWeights += oscDummy->GetBinContent(x,y)*(1.0-OscWeight);
      }
      else{
	sumWeights += oscDummy->GetBinContent(x,y)*(OscWeight);
      }
      Double_t currBinContents = bintemplate->GetBinContent( y );
      bintemplate->SetBinContent( y, sumWeights + currBinContents);
    }
  }

  return bintemplate;
}
//---------------------------------------------------------------------------------
void dataRelease_chi2Calc_compile(string path = "dataRelease.root",
					double Dm232 = -0.00246058,
					double Dm221 = 0.0000754,
					double th23 = 0.908967,
					double th12 = 0.5540758073,
					double th13 = 0.149116,
					double deltaCP = 0.0,
					double Dm241 = 0.000102356,
					double th24 = 0.0120973,
					double th34 = 1.91265E-7,
					double th14 = 0.0,
					double delta24 = 0.0
				       )
{
  TString fileName = path;
  
  TFile *f  = new TFile(fileName);

  Zombie(f); 

  params my_pars;
  my_pars.Dm232   = Dm232;
  my_pars.Dm221   = Dm221;
  my_pars.th23    = th23;
  my_pars.th12    = th12;
  my_pars.th13    = th13;
  my_pars.deltaCP = deltaCP;
  my_pars.Dm241   = Dm241;
  my_pars.th24    = th24;
  my_pars.th34    = th34;
  my_pars.th14    = th14;
  my_pars.delta24 = delta24;

  PrintParStatus(my_pars);

  params my_pars0;
  my_pars0.Dm232 = 0.0;
  my_pars0.Dm221 = 0.0;
  my_pars0.th23  = 0.0;
  my_pars0.th13  = 0.0;
  my_pars0.th12  = 0.0;

////////////////////////Begin Calculation//////////////////////////////////////
  CoVarCCinvert = (TMatrixD*)f->Get("TotalInvertCC"); assert(CoVarCCinvert);
  CoVarNCinvert = (TMatrixD*)f->Get("TotalInvertNC"); assert(CoVarNCinvert);
  
  f->GetObject("hRecoToTrueNDNCSelectedTrueNC",   NDNC_TrueNC);   assert(NDNC_TrueNC);
  f->GetObject("hRecoToTrueNDNCSelectedNuMu",     NDNC_NuMu);     assert(NDNC_NuMu);
  f->GetObject("hRecoToTrueNDNCSelectedBeamNue",  NDNC_BeamNue);  assert(NDNC_BeamNue);
  f->GetObject("hRecoToTrueNDNCSelectedAppNue",   NDNC_AppNue);   assert(NDNC_AppNue);
  f->GetObject("hRecoToTrueNDNCSelectedAppNuTau", NDNC_AppNuTau); assert(NDNC_AppNuTau);

  f->GetObject("hRecoToTrueNDCCSelectedTrueNC",   NDCC_TrueNC);   assert(NDCC_TrueNC);
  f->GetObject("hRecoToTrueNDCCSelectedNuMu",     NDCC_NuMu);     assert(NDCC_NuMu);
  f->GetObject("hRecoToTrueNDCCSelectedBeamNue",  NDCC_BeamNue);  assert(NDCC_BeamNue);
  f->GetObject("hRecoToTrueNDCCSelectedAppNue",   NDCC_AppNue);   assert(NDCC_AppNue);
  f->GetObject("hRecoToTrueNDCCSelectedAppNuTau", NDCC_AppNuTau); assert(NDCC_AppNuTau);

  f->GetObject("hRecoToTrueFDNCSelectedTrueNC",   FDNC_TrueNC);   assert(FDNC_TrueNC);
  f->GetObject("hRecoToTrueFDNCSelectedNuMu",     FDNC_NuMu);     assert(FDNC_NuMu);
  f->GetObject("hRecoToTrueFDNCSelectedBeamNue",  FDNC_BeamNue);  assert(FDNC_BeamNue);
  f->GetObject("hRecoToTrueFDNCSelectedAppNue",   FDNC_AppNue);   assert(FDNC_AppNue);
  f->GetObject("hRecoToTrueFDNCSelectedAppNuTau", FDNC_AppNuTau); assert(FDNC_AppNuTau);

  f->GetObject("hRecoToTrueFDCCSelectedTrueNC",   FDCC_TrueNC);   assert(FDCC_TrueNC);
  f->GetObject("hRecoToTrueFDCCSelectedNuMu",     FDCC_NuMu);     assert(FDCC_NuMu);
  f->GetObject("hRecoToTrueFDCCSelectedBeamNue",  FDCC_BeamNue);  assert(FDCC_BeamNue);
  f->GetObject("hRecoToTrueFDCCSelectedAppNue",   FDCC_AppNue);   assert(FDCC_AppNue);
  f->GetObject("hRecoToTrueFDCCSelectedAppNuTau", FDCC_AppNuTau); assert(FDCC_AppNuTau);

  f->GetObject("dataFDNC", FD_dataNC); assert(FD_dataNC);
  f->GetObject("dataFDCC", FD_dataCC); assert(FD_dataCC);

  f->GetObject("dataNDNC", ND_dataNC); assert(ND_dataNC);
  f->GetObject("dataNDCC", ND_dataCC); assert(ND_dataCC);

  TH1D* RdataCC = (TH1D*)FD_dataCC->Clone();
  RdataCC->Divide(ND_dataCC);

  TH1D* RdataNC = (TH1D*)FD_dataNC->Clone();
  RdataNC->Divide(ND_dataNC);

  TH1D* NDUnOscCC_MC = (TH1D*)CreateTotalSpectrum(my_pars0, NDCC_TrueNC, NDCC_NuMu, NDCC_BeamNue, NDCC_AppNue, NDCC_AppNuTau, 1.04*kKmUnits);
  TH1D* NDUnOscNC_MC = (TH1D*)CreateTotalSpectrum(my_pars0, NDNC_TrueNC, NDNC_NuMu, NDNC_BeamNue, NDNC_AppNue, NDNC_AppNuTau, 1.04*kKmUnits);
  
  FDOscCC_MC = (TH1D*)CreateTotalSpectrum(my_pars, FDCC_TrueNC, FDCC_NuMu, FDCC_BeamNue, FDCC_AppNue, FDCC_AppNuTau, 735.0*kKmUnits);
  NDOscCC_MC = (TH1D*)CreateTotalSpectrum(my_pars, NDCC_TrueNC, NDCC_NuMu, NDCC_BeamNue, NDCC_AppNue, NDCC_AppNuTau, 1.04*kKmUnits);
  TH1D* RpredCC = (TH1D*)FDOscCC_MC->Clone();
  RpredCC->Divide(NDOscCC_MC);
  
  FDOscNC_MC = (TH1D*)CreateTotalSpectrum(my_pars, FDNC_TrueNC, FDNC_NuMu, FDNC_BeamNue, FDNC_AppNue, FDNC_AppNuTau, 735.0*kKmUnits);
  NDOscNC_MC = (TH1D*)CreateTotalSpectrum(my_pars, NDNC_TrueNC, NDNC_NuMu, NDNC_BeamNue, NDNC_AppNue, NDNC_AppNuTau, 1.04*kKmUnits);
  TH1D* RpredNC = (TH1D*)FDOscNC_MC->Clone();
  RpredNC->Divide(NDOscNC_MC);  

  Double_t chi2 = ComparePredWithData(RpredCC, 
	  			      NDOscCC_MC, 
				      NDUnOscCC_MC, 
				      RdataCC, 
				      ND_dataCC, 
				      CoVarCCinvert, 
				      RpredNC, 
				      NDOscNC_MC, 
				      NDUnOscNC_MC, 
				      RdataNC, 
				      ND_dataNC, 
				      CoVarNCinvert,
				      my_pars.Dm232);

  TString outputfilename = "./chi2Calc_output.root";

  TFile* save = new TFile(outputfilename,"RECREATE");
  save->cd();
  
  TTree* tree = new TTree("tree", "tree");
  tree->Branch("chi2", &chi2 );
  tree->Branch("output_chisqCC", &chisqCC );
  tree->Branch("output_chisqNC", &chisqNC );
  tree->Branch("output_NDnormPenCC", &Penalty_NormCC );
  tree->Branch("output_NDnormPenNC", &Penalty_NormNC );
  tree->Branch("output_PenDm232", &Penalty_dm232 );
  tree->Branch("output_totalChi2CC", &totalChi2_CC );
  tree->Branch("output_totalChi2NC", &totalChi2_NC );
  tree->Branch("param_dm232", &my_pars.Dm232 );
  tree->Branch("param_dm221", &my_pars.Dm221 );
  tree->Branch("param_theta23", &my_pars.th23 );
  tree->Branch("param_theta12", &my_pars.th12 );
  tree->Branch("param_theta13", &my_pars.th13 );
  tree->Branch("param_deltaCP", &my_pars.deltaCP );
  tree->Branch("param_dm241", &my_pars.Dm241 );
  tree->Branch("param_th24", &my_pars.th14 );
  tree->Branch("param_th24", &my_pars.th24 );
  tree->Branch("param_th34", &my_pars.th34 );
  tree->Branch("param_delta24", &my_pars.delta24 );
  tree->Fill();
  tree->Write();
 
  RpredCC->Write("FarOverNearRatio_prediction_CC"); 
  RpredNC->Write("FarOverNearRatio_prediction_NC"); 
  RdataCC->Write("FarOverNearRatio_data_CC"); 
  RdataNC->Write("FarOverNearRatio_data_NC");
  NDOscCC_MC->Write("NearDetector_prediction_CC"); 
  NDUnOscCC_MC->Write("NearDetector_unoscillated_CC"); 
  ND_dataCC->Write("NearDetector_data_CC"); 
  NDOscNC_MC->Write("NearDetector_prediction_NC"); 
  NDUnOscNC_MC->Write("NearDetector_unoscillated_NC"); 
  ND_dataNC->Write("NearDetector_data_NC");
  CoVarCCinvert->Write("Inverse_CovarianceMatrix_CC");
  CoVarNCinvert->Write("Inverse_CovarianceMatrix_NC"); 
  
  save->Close();
 
  std::cout << " " << std::endl;
}
