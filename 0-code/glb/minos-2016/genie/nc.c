// ----------------------------------------------------------------------------
// nc.c: Create GLoBES smearing matrices from GENIE event simulation
// ----------------------------------------------------------------------------
// Author: Joachim Kopp (Fermilab)
// Date:   2011
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
using namespace std;

#include <PDG/PDGUtils.h>
#include <GHEP/GHepParticle.h>
#include <Ntuple/NtpMCTreeHeader.h>
#include <Ntuple/NtpMCEventRecord.h>

// Show histograms?
//#define GRAPHICS_ON

// Maximum number of events to analyze
#define MAX_EVENTS 10000

const double M_MU = 105.658367e6; // eV

const double E_min =   0.0;      // GeV
const double E_max = 120.0;      // GeV
const int E_bins   = 170;
const double E_binsize[E_bins] = {
  0.25, 0.25, 0.25, 0.25,  0.25, 0.25, 0.25, 0.25,  0.25, 0.25, 0.25, 0.25,
  0.25, 0.25, 0.25, 0.25,  0.25, 0.25, 0.25, 0.25,  0.25, 0.25, 0.25, 0.25,
  0.25, 0.25, 0.25, 0.25,  0.25, 0.25, 0.25, 0.25,  0.25, 0.25, 0.25, 0.25,
  0.25, 0.25, 0.25, 0.25,  0.5, 0.5,   0.5, 0.5,    0.5, 0.5,   0.5, 0.5,
  0.5, 0.5,   0.5, 0.5,    0.5, 0.5,   0.5, 0.5,    0.5, 0.5,   0.5, 0.5,
  0.5, 0.5,   0.5, 0.5,    0.5, 0.5,   0.5, 0.5,    0.5, 0.5,   0.5, 0.5,
  0.5, 0.5,   0.5, 0.5,    0.5, 0.5,   0.5, 0.5,    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,    1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};
double E_bin_centers[E_bins];
double E_bin_edges[E_bins+1];

const double E_reco_min =  0.0; // GeV
const double E_reco_max = 40.0; // GeV
//const int E_reco_bins   = 43;
//const double E_reco_binsize[E_reco_bins] = { // 1607.01176 CC selection
//  1,
//  0.2, 0.2, 0.2, 0.2, 0.2,  0.2, 0.2, 0.2, 0.2, 0.2,  0.2, 0.2, 0.2, 0.2, 0.2,
//  0.2, 0.2, 0.2, 0.2, 0.2,  1, 1, 1, 1, 1,  1, 1, 1, 1, 1,  1, 1, 1, 1, 1,
//  2, 2, 2, 2, 2,  5, 5
//};
const int E_reco_bins   = 29;
const double E_reco_binsize[E_reco_bins] = { // 1607.01176 NC selection
  1,
  0.2, 0.2, 0.2, 0.2, 0.2,  0.2, 0.2, 0.2, 0.2, 0.2,  0.2, 0.2, 0.2, 0.2, 0.2,
  0.2, 0.2, 0.2, 0.2, 0.2,  3, 3, 3, 3, 3,
  5, 5,  10
};
double E_reco_bin_centers[E_bins];
double E_reco_bin_edges[E_bins+1];

double smearing_matrix[E_reco_bins][E_bins];

const int debug_level = 0;

TApplication *app = NULL;
TCanvas *c1       = NULL;
TPad *pad1        = NULL;
TH1F *hEreco      = NULL;


// -------------------------------------------------------------------------------------
void generate_events(int flavor, double E)
// -------------------------------------------------------------------------------------
// Generate events for a monoenergetic neutrino beam.
// -------------------------------------------------------------------------------------
{
  printf("------------ GENERATING EVENTS: flavor = %d, E = %g -----------\n", flavor, E);

  int run_number = int(1000.0*E+0.5);

  // Generate events
  char command[1024];
  sprintf(command, "gevgen -e %g -p %d -r %d -n 50000 " // 50000 ev. at given E, flavor, run #
                   "--cross-sections gxspl/gxspl-fe56.xml "        // Cross section file
                   "--event-generator-list NC "                   // Generate only NC events
                   "--message-thresholds Messenger_laconic.xml "  // Be less verbose
                   "-t 1000260560",
          E, flavor, run_number);
  printf("%s\n", command);
  system(command);

  sprintf(command, "mv gntp.%d.ghep.root local/ev-%d-%5.3lf.root", run_number, flavor, E);
  system(command);
}


// -------------------------------------------------------------------------------------
void process_events(int flavor, double E)
// -------------------------------------------------------------------------------------
// Analyze events for a monoenergetic neutrino beam.
// -------------------------------------------------------------------------------------
{
  printf("------------ PROCESSING EVENTS: flavor %d, E = %g -----------\n", flavor, E);

  // Histogram of visible energy
  if (hEreco) { delete hEreco; hEreco = NULL; }
  char plot_label[100];
  sprintf(plot_label, "E_reco distribution for E_nu = %g", E);
  hEreco = new TH1F("hEreco", plot_label, E_reco_bins, E_reco_bin_edges);
//  hEreco = new TH1F("hEreco", plot_label, E_reco_bins, E_reco_min, E_reco_bins);

  // Open input file
  char infile[255];
  sprintf(infile, "local/ev-%d-%5.3lf.root", flavor, E);
  TFile *f = new TFile(infile);
  TTree *t = (TTree *) f->Get("gtree");
  genie::NtpMCEventRecord *mcrec = NULL;
  t->SetBranchAddress("gmcrec", &mcrec);

  // Loop over all events --- fill histograms and apply cuts
  for (int i=0; i < t->GetEntries(); i++)
  {
    if (i > MAX_EVENTS)
      break;
    t->GetEntry(i);

    genie::NtpMCRecHeader rec_header = mcrec->hdr;
    genie::EventRecord *ev           = mcrec->event;
    genie::GHepParticle *p           = NULL;

    if (debug_level > 0)
    {
      printf("Event %d\n", i);
      cout << *ev;
      printf("Final state: ");
      for (int j=0; j < ev->GetEntries(); j++)
      {
        p = ev->Particle(j);
        if (!p->HasDaughters())
          printf("%s ", p->Name().c_str());
      }
      printf("\n");
    }

    // Cuts
    if (!ev->Summary()->ProcInfo().IsWeakNC())
    {
      printf("WARNING: Non-NC event type (%s) found in input file.\n",
             ev->Summary()->ProcInfo().AsString().c_str());
      continue;
    }

    // Determine visible energy in this event
    double Ei    = 0.0;                  // Initial energy in the detector (without E_nu)
    double Ef    = 0.0;                  // Final energy in the detector (without E_nu)
    for (int j=0; j < ev->GetEntries(); j++)
    {
      p = ev->Particle(j);
      if (genie::pdg::IsNeutralLepton(p->Pdg()))       // Ignore neutrinos
        continue;
      if (p->FirstMother() < 0 && p->LastMother() < 0) // Initial state particle?
        Ei += p->Energy();
      if (!p->HasDaughters())                          // Final state particle?
        Ef += p->Energy();
    }
    if (Ef-Ei < 0)
    {
      fprintf(stderr, "WARNING: Negative energy deposit.\n");
      continue;
    }
    hEreco->Fill(Ef - Ei);
    if (debug_level > 0)
    {
      printf("Initial energy: %g GeV\n", Ei);
      printf("Final energy: %g GeV\n", Ef);
      printf("Visible energy: %g GeV\n", Ef-Ei);
      getchar();
    }
  }
  hEreco->Sumw2();

  // Normalize by number of events in input file
  hEreco->Scale(1./hEreco->GetEntries());

  double E_tmp = E_min; // Determine the energy bin corresponding to input energy E
  int j;
  for (j=0; j < E_bins; j++)
  {
    if (E_tmp <= E && E < E_tmp + E_binsize[j])
      break;
    E_tmp += E_binsize[j];
  }
//  int j = int(E_bins * (E - E_min) / (E_max - E_min));
  for (int i=0; i < E_reco_bins; i++)
  {
    if (j >= 0 && j < E_bins)
      smearing_matrix[i][j] = hEreco->GetBinContent(i+1);
  }

  #ifdef GRAPHICS_ON
    c1->cd();
    pad1->Draw();
    pad1->cd();
    hEreco->Draw();
    c1->Update();
    app->Run(kTRUE);
  #endif

  delete t;  t = NULL;
  f->Close();
  delete f;  f = NULL;
}


// -------------------------------------------------------------------------------------
int main(int argc, char **argv)
// -------------------------------------------------------------------------------------
{
  // Initializa root
  int appargc = 1;
  char *appargv[] = { argv[0] };
  app = new TApplication(appargv[0], &appargc, appargv);
  #ifdef GRAPHICS_ON
    c1     = new TCanvas("c1", "GENIE neutrino events",1000,500);
    pad1   = new TPad("pad1", "E_reco distribution", 0.03, 0.03, 0.97, 0.97);
  #endif

  // Create arrays of bin centers and bin edges
  E_bin_edges[0] = E_min;
  for (int i=0; i < E_bins; i++)
  {
    E_bin_centers[i] = E_bin_edges[i] + 0.5 * E_binsize[i];
    E_bin_edges[i+1] = E_bin_edges[i] + E_binsize[i];
  }
  E_reco_bin_edges[0] = E_reco_min;
  for (int i=0; i < E_reco_bins; i++)
  {
    E_reco_bin_centers[i] = E_reco_bin_edges[i] + 0.5 * E_reco_binsize[i];
    E_reco_bin_edges[i+1] = E_reco_bin_edges[i] + E_reco_binsize[i];
  }

  // Generate events
//  omp_set_num_threads(8);
//#pragma omp parallel for
//  for (int k=0; k < E_bins; k++)
//  {
//    generate_events(14, E_bin_centers[k]);  // FIXME We should do that separately for anti-nu!
//  }

  // Process events
  for (int k=0; k < E_bins; k++)
  {
    process_events(14, E_bin_centers[k]);
  }

  // Consistency check: does every column of the smearing matrix add up to 1?
  for (int j=0; j < E_bins; j++)
  {
    double t = 0.0;
    for (int i=0; i < E_reco_bins; i++)
      t += smearing_matrix[i][j];
    printf("CHECK: %g\n", t);
  }

  // Write output as GLoBES smearing matrix
  const char file_out[] = "smear-nc.dat";
  FILE *fOut = fopen(file_out, "w");
  fprintf(fOut, "energy(#ERES_NC)<\n");
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
}


