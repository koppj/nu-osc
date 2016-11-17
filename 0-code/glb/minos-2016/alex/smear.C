// -------------------------------------------------------------------
// Read histogram from root file and export as GLoBES smearing matrix
// -------------------------------------------------------------------
//#define FAR
#define CC
int smear()
{
#ifdef FAR
  TFile f("FD_NC_histos.root");
#else
  TFile f("ND_NC_histos.root");
#endif
  TH2F *hRecoTrueNC = (TH2F *) f.Get("RecoTrueNC");

#ifdef FAR
  #ifdef CC
    FILE *fOut = fopen("smear_nc_cc_analysis-far.dat", "w");
  #else
    FILE *fOut = fopen("smear_nc_nc_analysis-far.dat", "w");
  #endif
#else
  #ifdef CC
    FILE *fOut = fopen("smear_nc_cc_analysis-near.dat", "w");
  #else
    FILE *fOut = fopen("smear_nc_nc_analysis-near.dat", "w");
  #endif
#endif
  fprintf(fOut, "// Mapping of true neutrino energy to reconstructed energy\n");
  fprintf(fOut, "// in NC interactions in MINOS\n");
  fprintf(fOut, "// Extracted from histograms courtesy of Alex Sousa\n");
  fprintf(fOut, "energy(#ERES_NC_TEST)<\n");
  fprintf(fOut, "  @energy =\n");


  // Compute normalization constants
#ifdef CC
  const int nbins_reco =  43; // Binning for CC analysis from 1607.01176
  double reco_binwidth[] = {  // Analysis bins from glb file
    1,
    0.2, 0.2, 0.2, 0.2, 0.2,  0.2, 0.2, 0.2, 0.2, 0.2,  0.2, 0.2, 0.2, 0.2, 0.2,
    0.2, 0.2, 0.2, 0.2, 0.2,  1, 1, 1, 1, 1,  1, 1, 1, 1, 1,  1, 1, 1, 1, 1,
    2, 2, 2, 2, 2,  5, 5
  };
#else
  const int nbins_reco =  29; // Binning for NC analysis from 1607.01176
  double reco_binwidth[] = {  // Analysis bins from glb file
    1,
    0.2, 0.2, 0.2, 0.2, 0.2,  0.2, 0.2, 0.2, 0.2, 0.2,  0.2, 0.2, 0.2, 0.2, 0.2,
    0.2, 0.2, 0.2, 0.2, 0.2,  3, 3, 3, 3, 3,
    5, 5,  10
  };
#endif

  const int nbins_true = 170;
  double true_binwidth[] = { // Stepsizes for NC analysis (from glb file)
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

  // The smearing matrix - note that bin boundaries must coincide with
  // bin boundaries in the original histogram, though bins may be combined
  double smearing_matrix[nbins_reco][nbins_true];

//  fprintf(stderr, "%d %d\n", hRecoTrueNC->GetYaxis()->GetNbins(),
//                             hRecoTrueNC->GetXaxis()->GetNbins());
//  for (int i=0; i < hRecoTrueNC->GetYaxis()->GetNbins(); i++)
//    fprintf(stderr, "%g\n", hRecoTrueNC->GetYaxis()->GetBinWidth(i));
//  fprintf(stderr, "\n");
//  for (int i=0; i < hRecoTrueNC->GetXaxis()->GetNbins(); i++)
//    fprintf(stderr, "%g\n", hRecoTrueNC->GetXaxis()->GetBinWidth(i));

  // Loops over AEDL bins and sampling points
  double E_reco = 0.0; // Left edge of current reco bin in smearing_matrix
  for (int i=0; i < nbins_reco; i++)
  {
    double E_true = 0.0;  // Left edge of current true bin in smearing_matrix
    for (int j=0; j < nbins_true; j++) 
    {
      smearing_matrix[i][j] = 0.0;

      int k = hRecoTrueNC->GetYaxis()->FindFixBin(E_true + 1e-10);
      while (hRecoTrueNC->GetYaxis()->GetBinUpEdge(k) <= E_true + true_binwidth[j] + 1e-10)
      {
        int l = hRecoTrueNC->GetXaxis()->FindFixBin(E_reco + 1e-10);
        while (hRecoTrueNC->GetXaxis()->GetBinUpEdge(l) <= E_reco + reco_binwidth[i] + 1e-10)
        {
          smearing_matrix[i][j] += hRecoTrueNC->GetBinContent(l, k) / true_binwidth[j];
          l++;
        }
        k++;
      }

      E_true += true_binwidth[j];
    }
    E_reco += reco_binwidth[i];
  }


  for (int i=0; i < nbins_reco; i++)
  {
    fprintf(fOut, "{ %d, %d ", 0, nbins_true-1);
    for (int j=0; j < nbins_true; j++) 
      fprintf(fOut, ", %10.7g", smearing_matrix[i][j]);
    if (i != nbins_reco-1)
      fprintf(fOut, " }:\n");
  }
  fprintf(fOut, " };\n");
  fprintf(fOut, ">\n");
  fclose(fOut);
  f.Close();
  return 0;
}
