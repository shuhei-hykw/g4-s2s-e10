// -*- C++ -*-

// namespace fs = std::filesystem;

void acceptance()
{
  if(!gFile || !gFile->IsOpen()){
    std::cout << "Usage: " << std::endl
              << "   $ root foo.root " << __func__ << ".C" << std::endl;
    gSystem->Exit(EXIT_FAILURE);
  }
  gStyle->SetOptStat(0);
  gStyle->SetStatX(0.900);
  gStyle->SetStatY(0.940);

  TString fig_dir = TString(gSystem->DirName(__FILE__))+"/fig";
  TString fig_path = fig_dir + "/" + gSystem->BaseName(gFile->GetName());
  fig_path.ReplaceAll(".root", "_acc.pdf");

  auto c1 = new TCanvas("c1", "c1", 1200, 400);

  auto h_gen = gFile->Get<TH2D>("PRMPThetaGen");
  auto h_acc = gFile->Get<TH2D>("PRMPThetaAcc");
  auto h_eff = dynamic_cast<TH2D*>(h_acc->Clone("h_eff"));

  const Double_t pmin = h_eff->GetYaxis()->GetXmin();
  const Double_t pmax = h_eff->GetYaxis()->GetXmax();
  const Double_t tmin = h_eff->GetXaxis()->GetXmin();
  const Double_t tmax = h_eff->GetXaxis()->GetXmax();
  const Double_t dp = (pmax - pmin)/h_eff->GetNbinsY();
  const Double_t dt = (tmax - tmin)/h_eff->GetNbinsX();
  auto h_msr = new TH1D("h_msr", "Solid Angle; [GeV/c]; [msr]",
                        h_eff->GetNbinsY(), pmin, pmax);
  auto h_ang = new TH1D("h_ang", "Anglar Acceptance at 0.9 GeV/c; [deg];",
                        h_eff->GetNbinsX(), tmin, tmax);
  h_eff->Divide(h_gen);
  h_eff->SetStats(0);
  h_eff->SetTitle("PRM P%Theta (Eff)");
  h_msr->SetStats(0);
  h_ang->SetStats(0);

  for(Int_t ip=0, np=h_eff->GetNbinsY(); ip<np; ++ip){
    auto p = h_eff->GetYaxis()->GetBinCenter(ip+1);
    Double_t srsum = 0;
    for(Int_t it=0, nt=h_eff->GetNbinsX(); it<nt; ++it){
      auto t = h_eff->GetXaxis()->GetBinCenter(it+1);
      auto sr = 2*TMath::Pi()*(TMath::Cos((t - dt/2)*TMath::DegToRad()) -
                               TMath::Cos((t + dt/2)*TMath::DegToRad()));
      auto eff = h_eff->GetBinContent(it+1, ip+1);
      srsum += sr*eff;
      if(ip == np/2){
        h_ang->Fill(t, eff);
      }
    }
    // std::cout << p << "\t" << srsum << std::endl;
    h_msr->Fill(p, srsum*1e3);
  }

  c1->Divide(3, 1);
  c1->cd(1);
  // h_gen->Draw("colz");
  // h_acc->Draw("colz");
  h_eff->Draw("colz");
  c1->cd(2);
  h_msr->Draw("hist");
  c1->cd(3);
  h_ang->Draw("hist");
  // c1->cd(4);
  c1->Print(fig_path);
}
