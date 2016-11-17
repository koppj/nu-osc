// ----------------------------------------------------------------------------
// minos.c: Extracting the data from the MINOS data release and converting to
//          GLoBES format
// ----------------------------------------------------------------------------
// Author: Joachim Kopp (PRISMA Cluster of Excellence, Mainz)
// Date:   2016
// ----------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include "root/TApplication.h"
#include "root/TSystem.h"
#include "root/TROOT.h"
#include "root/TFile.h"
#include "root/TTree.h"
#include "root/TCanvas.h"
#include "root/TH1.h"
#include "root/TH2.h"
using namespace std;

#include <PDG/PDGUtils.h>
#include <GHEP/GHepParticle.h>
#include <Ntuple/NtpMCTreeHeader.h>
#include <Ntuple/NtpMCEventRecord.h>


#define SQR(x)  ((x)*(x))

const double L_ND =   0.7; // km
const double L_FD = 735.;  // km

// Binning in true energy
const double E_min =  0.0;      // GeV
const double E_max =  5.0;      // GeV
const int E_bins   = 40;
const double E_binsize[E_bins] = {
  0.125, 0.125, 0.125, 0.125,   0.125, 0.125, 0.125, 0.125,
  0.125, 0.125, 0.125, 0.125,   0.125, 0.125, 0.125, 0.125,
  0.125, 0.125, 0.125, 0.125,   0.125, 0.125, 0.125, 0.125,
  0.125, 0.125, 0.125, 0.125,   0.125, 0.125, 0.125, 0.125,
  0.125, 0.125, 0.125, 0.125,   0.125, 0.125, 0.125, 0.125
};
double E_bin_centers[E_bins];
double E_bin_edges[E_bins+1];

TApplication *app = NULL;
TFile *f = NULL;


// -------------------------------------------------------------------------------------
void process_smearing_matrix(TH2D *sm, double L, const char *sm_file,
                             const char *spect_file)
// -------------------------------------------------------------------------------------
// Extract information from MINOS smearing matrix and output in GLoBES-compatible form
// -------------------------------------------------------------------------------------
{
  // Binning of smearing matrices
  int E_reco_bins = sm->GetYaxis()->GetNbins();
  TAxis *x = sm->GetXaxis();
  TAxis *y = sm->GetYaxis();
  printf("Bin widths in reconstructed energy (Emin=%g, Emax=%g):\n",
         y->GetBinLowEdge(1), y->GetBinUpEdge(E_reco_bins));
  for (int i=1; i <= E_reco_bins; i++)
  {
    printf("%7.3g, ", y->GetBinWidth(i));
    if ((i-1) % 10 == 9)
      printf("\n");
  }
  printf("\n");
  printf("\n");

  // Construct smearing matrix
  double smearing_matrix[E_reco_bins][E_bins];
  double norm[E_bins];
  for (int i=0; i < E_reco_bins; i++)
    for (int j=0; j < E_bins; j++)
    {
      double LE = L / E_bin_centers[j];
      if (x->GetBinUpEdge(x->GetNbins()) > LE && LE > x->GetBinLowEdge(1))
        smearing_matrix[i][j] = sm->Interpolate(L/E_bin_centers[j], y->GetBinCenter(i+1));
      else
        smearing_matrix[i][j] = 0.0;
      norm[j] += smearing_matrix[i][j];
    }
  for (int i=0; i < E_reco_bins; i++)
    for (int j=0; j < E_bins; j++)
      if (norm[j] > 0)
        smearing_matrix[i][j] /= norm[j];

  // Write GLoBES smearing matrices to file
  FILE *fOut = fopen(sm_file, "w");
  fprintf(fOut, "energy(#ERES_ND_CC)<\n");
  fprintf(fOut, "  @energy =\n");
  for (int i=0; i < E_reco_bins; i++)
  {
    fprintf(fOut, "{ %d, %d ", 0, E_bins-1);
    for (int j=0; j < E_bins; j++)
      fprintf(fOut, ", %g", smearing_matrix[i][j]);
    if (i != E_reco_bins-1)
      fprintf(fOut, " }:\n");
  }
  fprintf(fOut, " };\n");
  fprintf(fOut, ">\n");
  fclose(fOut);

  // Write GLoBES neutrino spectrum
//  fOut = fopen(spect_file, "w");
//  fprintf(fOut, "# MINOS neutrino spectrum based on October 2016 data release\n");
//  fprintf(fOut, "#\n");
//  fprintf(fOut, "# E [GeV]    \\nu_e       \\nu_\\mu     \\nu_\\tau    "
//                "\\bar\\nu_e   \\bar\\nu_\\mu  \\bar\\nu_\\tau\n");
//  for (int i=0; i < E_bins; i++)
//  {
//    fprintf(fOut, "   %5.3g        0     %10.5g        0             "
//                  "0            0           0\n",
//            E_bin_centers[i], hNuSpectrum->GetBinContent(i+1));
//  }
//  fclose(fOut);
}


// -------------------------------------------------------------------------------------
int main(int argc, char **argv)
// -------------------------------------------------------------------------------------
{
  // Initializa root
  int appargc = 1;
  char *appargv[] = { argv[0] };
  app = new TApplication(appargv[0], &appargc, appargv);

  // Open MINOS data file
  f = new TFile("dataRelease.root");

  // Read and print binning information
  TH1D *dataFDCC = (TH1D*) f->Get("dataFDCC");
  TH1D *dataFDNC = (TH1D*) f->Get("dataFDNC");
  TH1D *dataNDCC = (TH1D*) f->Get("dataNDCC");
  TH1D *dataNDNC = (TH1D*) f->Get("dataNDNC");

  // Create arrays of bin centers and bin edges
  E_bin_edges[0] = E_min;
  for (int i=0; i < E_bins; i++)
  {
    E_bin_centers[i] = E_bin_edges[i] + 0.5 * E_binsize[i];
    E_bin_edges[i+1] = E_bin_edges[i] + E_binsize[i];
  }

  // Read and print data
  printf("Data ND CC:\n");
  for (int i=1; i <= dataNDCC->GetXaxis()->GetNbins(); i++)
  {
    printf("%5.lf, ", dataNDCC->GetBinContent(i));
    if ((i-1) % 10 == 9)
      printf("\n");
  }
  printf("\n");
  printf("\n");

  printf("Data ND NC:\n");
  for (int i=1; i <= dataNDNC->GetXaxis()->GetNbins(); i++)
  {
    printf("%5.lf, ", dataNDNC->GetBinContent(i));
    if ((i-1) % 10 == 9)
      printf("\n");
  }
  printf("\n");
  printf("\n");

  printf("Data FD CC:\n");
  for (int i=1; i <= dataFDCC->GetXaxis()->GetNbins(); i++)
  {
    printf("%5.lf, ", dataFDCC->GetBinContent(i));
    if ((i-1) % 10 == 9)
      printf("\n");
  }
  printf("\n");
  printf("\n");

  printf("Data FD NC:\n");
  for (int i=1; i <= dataFDNC->GetXaxis()->GetNbins(); i++)
  {
    printf("%5.lf, ", dataFDNC->GetBinContent(i));
    if ((i-1) % 10 == 9)
      printf("\n");
  }
  printf("\n");
  printf("\n");

  // Smearing matrices
  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueNDCCSelectedTrueNC"), L_ND,
                           "smear-ndcc-truenc.dat", "");
  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueNDCCSelectedNuMu"), L_ND,
                           "smear-ndcc-truenumu.dat", "");
  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueNDCCSelectedBeamNue"), L_ND,
                           "smear-ndcc-beamnue.dat", "");
  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueNDCCSelectedAppNue"), L_ND,
                           "smear-ndcc-appnue.dat", "");
  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueNDCCSelectedAppNuTau"), L_ND,
                           "smear-ndcc-appnutau.dat", "");

  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueNDNCSelectedTrueNC"), L_ND,
                           "smear-ndnc-truenc.dat", "");
  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueNDNCSelectedNuMu"), L_ND,
                           "smear-ndnc-truenumu.dat", "");
  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueNDNCSelectedBeamNue"), L_ND,
                           "smear-ndnc-beamnue.dat", "");
  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueNDNCSelectedAppNue"), L_ND,
                           "smear-ndnc-appnue.dat", "");
  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueNDNCSelectedAppNuTau"), L_ND,
                           "smear-ndnc-appnutau.dat", "");

  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueFDCCSelectedTrueNC"), L_ND,
                           "smear-fdcc-truenc.dat", "");
  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueFDCCSelectedNuMu"), L_ND,
                           "smear-fdcc-truenumu.dat", "");
  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueFDCCSelectedBeamNue"), L_ND,
                           "smear-fdcc-beamnue.dat", "");
  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueFDCCSelectedAppNue"), L_ND,
                           "smear-fdcc-appnue.dat", "");
  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueFDCCSelectedAppNuTau"), L_ND,
                           "smear-fdcc-appnutau.dat", "");

  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueFDNCSelectedTrueNC"), L_ND,
                           "smear-fdnc-truenc.dat", "");
  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueFDNCSelectedNuMu"), L_ND,
                           "smear-fdnc-truenumu.dat", "");
  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueFDNCSelectedBeamNue"), L_ND,
                           "smear-fdnc-beamnue.dat", "");
  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueFDNCSelectedAppNue"), L_ND,
                           "smear-fdnc-appnue.dat", "");
  process_smearing_matrix((TH2D *) f->Get("hRecoToTrueFDNCSelectedAppNuTau"), L_ND,
                           "smear-fdnc-appnutau.dat", "");

  // Clean up
  f->Close();
  delete f;   f = NULL;
}


