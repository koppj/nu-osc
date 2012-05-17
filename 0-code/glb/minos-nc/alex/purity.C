// -------------------------------------------------------------------
// Read histograms from root file and compute purity of sample
// -------------------------------------------------------------------
#define FAR
int purity()
{
#ifdef FAR
  TFile f("FD_NC_histos.root");
#else
  TFile f("ND_NC_histos.root");
#endif
  TH1F *hReco   = (TH1F *) f.Get("Reco");
  TH1F *hRecoNC = (TH1F *) f.Get("RecoNC");

#ifdef FAR
  FILE *fOut = fopen("NC-pur-far-nu2010.dat", "w");
  fprintf(fOut, "# MINOS far detector purity for NC selected events\n");
#else
  FILE *fOut = fopen("NC-pur-near-nu2010.dat", "w");
  fprintf(fOut, "# MINOS near detector purity for NC selected events\n");
#endif
  fprintf(fOut, "# Extracted from histograms courtesy of Alex Sousa\n");
  fprintf(fOut, "#\n");
  fprintf(fOut, "# E_reco [GeV]      purity\n");

  for (int j=1; j <= hReco->GetNbinsX(); j++)
    fprintf(fOut, "%10.7g    %15.10g\n", hReco->GetBinCenter(j),
            hRecoNC->GetBinContent(j) / hReco->GetBinContent(j));
  fclose(fOut);
  f.Close();
  return 0;
}
