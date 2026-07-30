// acceptance_e10_kpi.C
//
// Draw the S-2S acceptance for 9Be(K-,pi+)9_Sigma_He as three separate
// (theta,p) maps: generated, accepted, and the acceptance-probability
// ratio, each its own clean COLZ canvas (not subplots of one combined
// figure) -- following the visual style of geant4/e07 (g4-kurama)'s
// AcceptanceDraw.C.
//
// Outputs:
//   fig/acc_generated_e10_kpi.pdf
//   fig/acc_accepted_e10_kpi.pdf
//   fig/acc_ratio_e10_kpi.pdf
//
// Usage:  root -l -b -q 'ana/acceptance_e10_kpi.C("path/to/run.root")'

void acceptance_e10_kpi(const char* fname)
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(99);

  TFile* file = TFile::Open(fname);
  if(!file || file->IsZombie()){
    std::cerr << "cannot open " << fname << std::endl; return;
  }
  TH2D* hGen = (TH2D*)file->Get("PRMPThetaGen");
  TH2D* hAcc = (TH2D*)file->Get("PRMPThetaAcc");
  if(!hGen || !hAcc){
    std::cerr << "PRMPThetaGen/PRMPThetaAcc not found in "
              << fname << std::endl; return;
  }

  const char* axt = ";Emission angle #theta [deg];Momentum p [GeV/#it{c}]";

  // ---- acceptance-probability map  acc/gen  (per theta,p cell) ----
  TH2D* hRatio = (TH2D*)hAcc->Clone("hRatio");
  hRatio->SetTitle(Form("Acceptance probability%s", axt));
  hRatio->Divide(hAcc, hGen, 1., 1., "B");

  // Paint genuine zero-acceptance cells (generated but never accepted) at
  // the bottom of the palette: COLZ leaves content-zero bins white, which
  // looks like missing data, so nudge them to a tiny epsilon.  Cells with
  // no generated events stay white.
  const Double_t eps = 1e-12;
  for(Int_t ix=1; ix<=hRatio->GetNbinsX(); ++ix)
    for(Int_t iy=1; iy<=hRatio->GetNbinsY(); ++iy)
      if(hGen->GetBinContent(ix,iy)>0 && hRatio->GetBinContent(ix,iy)==0)
        hRatio->SetBinContent(ix, iy, eps);

  // ---- draw the three 2D maps ----
  auto draw2d = [&](TH2* h, const char* out, const char* opt="COLZ"){
    TCanvas* c = new TCanvas(Form("c_%s", h->GetName()), h->GetName(),
                             760, 620);
    c->SetRightMargin(0.15);
    h->Draw(opt);
    c->SaveAs(out);
  };
  hGen->SetTitle(Form("Generated #pi^{+}%s", axt));
  hAcc->SetTitle(Form("Accepted #pi^{+}%s", axt));
  draw2d(hGen,   "ana/fig/acc_generated_e10_kpi.pdf");
  draw2d(hAcc,   "ana/fig/acc_accepted_e10_kpi.pdf");
  hRatio->SetMinimum(0.); hRatio->SetMaximum(1.);
  draw2d(hRatio, "ana/fig/acc_ratio_e10_kpi.pdf");

  std::cout << "generated = " << (Long64_t)hGen->GetEntries()
            << ", accepted = " << (Long64_t)hAcc->GetEntries() << std::endl;
}
