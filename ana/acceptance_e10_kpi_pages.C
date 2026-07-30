// acceptance_e10_kpi_pages.C
//
// Same three quantities as acceptance.C (efficiency map, solid angle
// vs momentum, angular acceptance at a fixed momentum), but each on
// its own full page of a single multi-page PDF instead of three
// cramped subplots of one canvas.
//
// Output: fig/acc_pages_e10_kpi.pdf (page 1: efficiency map,
//         page 2: solid angle vs p, page 3: angular acceptance)
//
// Usage:  root -l -b -q 'ana/acceptance_e10_kpi_pages.C("run.root")'

void acceptance_e10_kpi_pages(const char* fname)
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(99);

  TFile* file = TFile::Open(fname);
  if(!file || file->IsZombie()){
    std::cerr << "cannot open " << fname << std::endl; return;
  }
  TH2D* h_gen = (TH2D*)file->Get("PRMPThetaGen");
  TH2D* h_acc = (TH2D*)file->Get("PRMPThetaAcc");
  if(!h_gen || !h_acc){
    std::cerr << "PRMPThetaGen/PRMPThetaAcc not found in "
              << fname << std::endl; return;
  }
  TH2D* h_eff = (TH2D*)h_acc->Clone("h_eff");
  h_eff->Divide(h_gen);
  h_eff->SetTitle(
    ";Emission angle #theta [deg];Momentum p [GeV/#it{c}]");
  h_eff->SetMinimum(0.); h_eff->SetMaximum(1.);

  const Double_t pmin = h_eff->GetYaxis()->GetXmin();
  const Double_t pmax = h_eff->GetYaxis()->GetXmax();
  const Double_t tmin = h_eff->GetXaxis()->GetXmin();
  const Double_t tmax = h_eff->GetXaxis()->GetXmax();
  const Double_t dp = (pmax - pmin)/h_eff->GetNbinsY();
  const Double_t dt = (tmax - tmin)/h_eff->GetNbinsX();

  auto h_msr = new TH1D("h_msr", ";Momentum p [GeV/#it{c}];Solid angle [msr]",
                        h_eff->GetNbinsY(), pmin, pmax);
  const Int_t ip_cent = h_eff->GetNbinsY()/2;
  const Double_t p_cent = h_eff->GetYaxis()->GetBinCenter(ip_cent+1);
  auto h_ang = new TH1D("h_ang",
    Form("Angular acceptance at p = %.2f GeV/#it{c};"
         "Emission angle #theta [deg];Acceptance probability",
         p_cent),
    h_eff->GetNbinsX(), tmin, tmax);

  for(Int_t ip=0, np=h_eff->GetNbinsY(); ip<np; ++ip){
    auto p = h_eff->GetYaxis()->GetBinCenter(ip+1);
    Double_t srsum = 0;
    for(Int_t it=0, nt=h_eff->GetNbinsX(); it<nt; ++it){
      auto t = h_eff->GetXaxis()->GetBinCenter(it+1);
      auto sr = 2*TMath::Pi()*(TMath::Cos((t - dt/2)*TMath::DegToRad()) -
                               TMath::Cos((t + dt/2)*TMath::DegToRad()));
      auto eff = h_eff->GetBinContent(it+1, ip+1);
      srsum += sr*eff;
      if(ip == ip_cent){
        h_ang->SetBinContent(it+1, eff);
      }
    }
    h_msr->SetBinContent(ip+1, srsum*1e3);
  }

  TString out = "ana/fig/acc_pages_e10_kpi.pdf";
  TCanvas c("c", "c", 900, 700);

  c.SetRightMargin(0.15);
  h_eff->SetTitle(Form("Acceptance probability%s", h_eff->GetTitle()));
  h_eff->Draw("COLZ");
  c.Print(out + "(");

  c.Clear();
  c.SetRightMargin(0.10); c.SetLeftMargin(0.13);
  h_msr->SetTitle("S-2S acceptance;Momentum p [GeV/#it{c}];Solid angle [msr]");
  h_msr->SetMinimum(0.);
  h_msr->SetMarkerStyle(20);
  h_msr->SetMarkerSize(0.8);
  h_msr->SetLineColor(kAzure+2);
  h_msr->SetMarkerColor(kAzure+2);
  h_msr->Draw("E1");
  c.Print(out);

  c.Clear();
  h_ang->SetMinimum(0.);
  h_ang->SetMaximum(1.05);
  h_ang->SetMarkerStyle(20);
  h_ang->SetMarkerSize(0.8);
  h_ang->SetLineColor(kAzure+2);
  h_ang->SetMarkerColor(kAzure+2);
  h_ang->Draw("E1");
  c.Print(out + ")");

  std::cout << "Saved: " << out << std::endl;
}
