// -------------------------------------------------------------------
// Read histogram from root file and export as GLoBES smearing matrix
// -------------------------------------------------------------------
//#define FAR
int smear()
{
#ifdef FAR
  TFile f("FD_NC_histos.root");
#else
  TFile f("ND_NC_histos.root");
#endif
  TH2F *hRecoTrueNC = (TH2F *) f.Get("RecoTrueNC");

#ifdef FAR
  FILE *fOut = fopen("smear_nc_mc-far.dat", "w");
#else
  FILE *fOut = fopen("smear_nc_mc-near.dat", "w");
#endif
  fprintf(fOut, "// Mapping of true neutrino energy to reconstructed energy\n");
  fprintf(fOut, "// in NC interactions in MINOS\n");
  fprintf(fOut, "// Extracted from histograms courtesy of Alex Sousa\n");
  fprintf(fOut, "energy(#ERES_NC_TEST)<\n");
  fprintf(fOut, "  @energy =\n");


  // Compute normalization constants
//  const int nbins_reco =  39; // Binning for CC analysis from 1103.0340
//  double reco_binwidth[] = {  // Analysis bins from glb file
//    1,
//    0.25, 0.25, 0.25, 0.25,  0.25, 0.25, 0.25, 0.25,  0.25, 0.25, 0.25, 0.25,
//    0.25, 0.25, 0.25, 0.25,  0.5, 0.5,   0.5, 0.5,    0.5, 0.5,   0.5, 0.5,
//    0.5, 0.5,   1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  10,  20
//  };
  const int nbins_reco =  20; // Binning for NC analysis from 1001.0336
                              // and for comparing ND data to fig. 1 of 1103.0340
  double reco_binwidth[] = {  // Analysis bins from glb file
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,    1, 1, 1, 1, 1, 1, 1, 1, 1, 1
  };

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
