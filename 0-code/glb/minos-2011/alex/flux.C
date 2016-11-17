// -------------------------------------------------------------------
// Extract initial neutrino flux from Alex Sousa's histograms
// -------------------------------------------------------------------
//#define FAR
int flux()
{
#ifdef FAR
  TFile f("FD_NC_histos.root");
#else
  TFile f("ND_NC_histos.root");
#endif
  TH2F *hRecoTrueNC = (TH2F *) f.Get("RecoTrueNC");
  hRecoTrueNC->Rebin2D(hRecoTrueNC->GetNbinsX(), 4);

#ifdef FAR
  FILE *fOut = fopen("MINOS-total-flux-far.dat", "w");
  fprintf(fOut, "# MINOS total neutrino flux (all flavors) as seen by the far detector\n");
#else
  FILE *fOut = fopen("MINOS-total-flux-near.dat", "w");
  fprintf(fOut, "# MINOS total neutrino flux (all flavors) as seen by the near detector\n");
#endif
  fprintf(fOut, "# Extracted from histograms courtesy of Alex Sousa\n");
  fprintf(fOut, "#\n");
  fprintf(fOut, "# E_true [GeV]      flux (a.u.)\n");

  for (int j=1; j <= hRecoTrueNC->GetNbinsY(); j++)
    fprintf(fOut, "%10.7g    %15.10g\n", hRecoTrueNC->GetYaxis()->GetBinCenter(j),
            hRecoTrueNC->GetBinContent(1,j));
  fclose(fOut);
  f.Close();
  return 0;
}
