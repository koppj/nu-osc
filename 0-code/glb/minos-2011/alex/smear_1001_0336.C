// -------------------------------------------------------------------
// Read histogram from root file and export as GLoBES smearing matrix
// -------------------------------------------------------------------
//#define FAR
int smear_1001_0336()
{
#ifdef FAR
  TFile f("FD_NC_histos.root");
#else
  TFile f("ND_NC_histos.root");
#endif
  TH2F *hRecoTrueNC = (TH2F *) f.Get("RecoTrueNC");
  hRecoTrueNC->Rebin2D(4, 1);

#ifdef FAR
  FILE *fOut = fopen("smear_NC_MC-far.dat", "w");
#else
  FILE *fOut = fopen("smear_NC_MC-near.dat", "w");
#endif
  fprintf(fOut, "// Mapping of true neutrino energy to reconstructed energy\n");
  fprintf(fOut, "// in NC interactions in MINOS\n");
  fprintf(fOut, "// Extracted from histograms courtesy of Alex Sousa\n");
  fprintf(fOut, "//\n");
//  fprintf(fOut, "//      E_true [GeV]    E_reco [GeV]    weight\n");
  fprintf(fOut, "energy(#ERES_NC)<\n");
  fprintf(fOut, "  @energy =\n");


  // Compute normalization constants
  int nbins_reco =  20; //hRecoTrueNC->GetNbinsX();
  int nbins_true = hRecoTrueNC->GetNbinsY();
  int n_sampling_points = 170;
  double sampling_stepsize[] = { // Stepsizes for NC analysis (from glb file)
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
  double *norm = new double[nbins_true];
  for (int j=1; j <= nbins_true; j++)
  {
    norm[j-1] = 0.0;
    for (int i=1; i <= nbins_reco; i++)
      norm[j-1] += hRecoTrueNC->GetBinContent(i, j);
  }

  for (int i=1; i <= nbins_reco; i++)
  {
    fprintf(fOut, "{ %d, %d ", 0, n_sampling_points-1);
    
    int k = 0;
    int j = 1;
    while (j <= nbins_true)      // Loop over E_true bins
    {
      double x = 0.0;
      double w = 0.0;
      while ((w += hRecoTrueNC->GetYaxis()->GetBinWidth(j)) <= sampling_stepsize[k]+1e-10)
      {
        x += hRecoTrueNC->GetBinContent(i, j);
        j++;
      }
      fprintf(fOut, ", %10.7g", x / sampling_stepsize[k]);
      k++;
    }
//    for (int j=1; j <= nbins_true; j++) // Loop over E_true bins
//      fprintf(fOut, ", %10.7g", hRecoTrueNC->GetBinContent(i, j));
    if (i != nbins_reco)
      fprintf(fOut, " }:\n");
  }
  fprintf(fOut, " };\n");
  fprintf(fOut, ">\n");
  fclose(fOut);
  f.Close();
  if (norm) { delete[] norm; norm = NULL; }
  return 0;
}
