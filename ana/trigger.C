// -*- C++ -*-

#include "DetectorID.hh"

std::map<Int_t, TString> pdg_map;
std::map<TString, Int_t> pdg_inv;

Int_t code(const TParticle& p)
{
  auto pdg = p.GetPDG();
  if(pdg){
    TString n = pdg->GetName();
    if(pdg_inv.count(n) > 0){
      return pdg_inv[n];
    }
    else{
      std::cout << n << std::endl;
      return pdg_inv.at("Others");
    }
  }
  else{
    auto code = p.GetPdgCode();
    TString n = "Others";
    if(code == 1000010020) return pdg_inv.at("d");
    else if(code == 1000010030) return pdg_inv.at("t");
    else if(code == 1000020040) return pdg_inv.at("{}^{4}He");
    else if(code == 1000070140) return pdg_inv.at("{}^{14}N");
    else if(code == 1000080160) return pdg_inv.at("{}^{16}O");
    else{
      std::cout << n << std::endl;
      return pdg_inv.at("Others");
    }
  }
}

void trigger()
{
  if(!gFile || !gFile->IsOpen()){
    std::cout << "Usage: " << std::endl
              << "   $ root foo.root " << __func__ << ".C" << std::endl;
    gSystem->Exit(EXIT_FAILURE);
  }

  pdg_map[0] = "pi-"; pdg_inv["pi-"] = 0;
  pdg_map[1] = "pi+"; pdg_inv["pi+"] = 1;
  pdg_map[2] = "mu-"; pdg_inv["mu-"] = 2;
  pdg_map[3] = "mu+"; pdg_inv["mu+"] = 3;
  pdg_map[4] = "proton"; pdg_inv["proton"] = 4;
  pdg_map[5] = "e-"; pdg_inv["e-"] = 5;
  pdg_map[6] = "e+"; pdg_inv["e+"] = 6;
  pdg_map[7] = "K+"; pdg_inv["K+"] = 7;
  pdg_map[8] = "d"; pdg_inv["d"] = 8;
  pdg_map[9] = "t"; pdg_inv["t"] = 9;
  pdg_map[10] = "{}^{4}He"; pdg_inv["{}^{4}He"] = 10;
  pdg_map[11] = "{}^{14}N"; pdg_inv["{}^{14}N"] = 11;
  pdg_map[12] = "{}^{16}O"; pdg_inv["{}^{16}O"] = 12;
  pdg_map[13] = "Others"; pdg_inv["Others"] = 13;

  TString fig_dir = TString(gSystem->DirName(__FILE__))+"/fig";
  TString fig_path = fig_dir + "/" + gSystem->BaseName(gFile->GetName());
  fig_path.ReplaceAll(".root", TString("_")+__func__+".pdf");

  TTreeReader reader("g4s2s", gFile);
  TTreeReaderValue<Int_t> evnum(reader, "evnum");
  TTreeReaderValue<std::vector<TParticle>> PRM(reader, "PRM");
  TTreeReaderValue<std::vector<TParticle>> SDC(reader, "SDC");
  TTreeReaderValue<std::vector<TParticle>> TOF(reader, "TOF");
  TTreeReaderValue<std::vector<TParticle>> AC1(reader, "AC1");
  TTreeReaderValue<std::vector<TParticle>> WC(reader, "WC");
  TTreeReaderValue<std::vector<TParticle>> VP(reader, "VP");

  auto hTOFHitPat = new TH1D("hTOFHitPat", "TOF HitPat",
                             NumOfSegTOF, -0.5, NumOfSegTOF-0.5);
  auto hTOFMulti = new TH1D("hTOFMulti", "TOF Multiplicity",
                            NumOfSegTOF+1, -0.5, NumOfSegTOF+0.5);
  auto hAC1HitPat = new TH1D("hAC1HitPat", "AC1 HitPat",
			     1, -0.5, 0.5);
  auto hAC1Multi = new TH1D("hAC1Multi", "AC1 Multiplicity",
                            2, -0.5, 1.5);
  auto hWCHitPat = new TH1D("hWCHitPat", "WC HitPat",
			    NumOfSegWC, -0.5, NumOfSegWC-0.5);
  auto hWCMulti = new TH1D("hWCMulti", "WC Multiplicity",
			   NumOfSegWC+1, -0.5, NumOfSegWC+0.5);

  auto hTOFWCHitPat = new TH2D("hTOFWCHitPat", "WC % TOF HitPat",
			       NumOfSegTOF, -0.5, NumOfSegTOF-0.5,
			       NumOfSegWC, -0.5, NumOfSegWC-0.5);

  auto hTOFPdgPat = new TH1D("hTOFPdgPat", "TOF PdgPat",
			     pdg_map.size(), 0, pdg_map.size());
  auto hAC1PdgPat = new TH1D("hAC1PdgPat", "AC1 PdgPat",
			     pdg_map.size(), 0, pdg_map.size());
  auto hWCPdgPat = new TH1D("hWCPdgPat", "WC PdgPat",
			    pdg_map.size(), 0, pdg_map.size());
  for(auto& pair: pdg_map){
    hTOFPdgPat->GetXaxis()->SetBinLabel(pair.first+1, pair.second);
    hAC1PdgPat->GetXaxis()->SetBinLabel(pair.first+1, pair.second);
    hWCPdgPat->GetXaxis()->SetBinLabel(pair.first+1, pair.second);
  }

  Int_t n_pi_mu = 0;
  Int_t n_proton = 0;
  while(reader.Next()){
    Bool_t tof_is_pi_mu = false;
    Bool_t tof_is_proton = false;
    Int_t tof_multiplicity = 0;
    for(const auto& p : (*TOF)){
      hTOFHitPat->Fill(p.GetMother(1));
      Int_t c = code(p);
      hTOFPdgPat->Fill(c);
      tof_multiplicity++;
      if(pdg_map[c].Contains("pi") ||
	 pdg_map[c].Contains("mu")){
	tof_is_pi_mu = true;
      }
      if(pdg_map[c].Contains("proton")){
	tof_is_proton = true;
      }
    }
    Bool_t ac1_is_pi_mu = false;
    Bool_t ac1_is_proton = false;
    Int_t ac1_multiplicity = 0;
    for(const auto& p : (*WC)){
      hAC1HitPat->Fill(p.GetMother(1));
      Int_t c = code(p);
      hAC1PdgPat->Fill(c);
      ac1_multiplicity++;
      if(pdg_map[c].Contains("pi") ||
	 pdg_map[c].Contains("mu")){
	ac1_is_pi_mu = true;
      }
      if(pdg_map[c].Contains("proton")){
	ac1_is_proton = true;
      }
    }
    Bool_t wc_is_pi_mu = false;
    Bool_t wc_is_proton = false;
    Int_t wc_multiplicity = 0;
    for(const auto& p : (*WC)){
      hWCHitPat->Fill(p.GetMother(1));
      Int_t c = code(p);
      hWCPdgPat->Fill(c);
      wc_multiplicity++;
      if(pdg_map[c].Contains("pi") ||
	 pdg_map[c].Contains("mu")){
	wc_is_pi_mu = true;
      }
      if(pdg_map[c].Contains("proton")){
	wc_is_proton = true;
      }
    }
    hTOFMulti->Fill(tof_multiplicity);
    hAC1Multi->Fill(ac1_multiplicity);
    hWCMulti->Fill(wc_multiplicity);

    if(tof_is_pi_mu && ac1_is_pi_mu && wc_is_pi_mu)
      n_pi_mu++;
    if(tof_is_proton && ac1_is_proton && wc_is_proton)
      n_proton++;

    // TOFxWC
    for(const auto& pTOF : (*TOF)){
      for(const auto& pWC : (*WC)){
	hTOFWCHitPat->Fill(pTOF.GetMother(1), pWC.GetMother(1));
      }
    }
  }

  std::cout << "n = " << reader.GetEntries() << std::endl;
  std::cout << "n_pi_mu = " << n_pi_mu << std::endl;
  std::cout << "n_proton = " << n_proton << std::endl;

  auto c1 = new TCanvas("c1", "c1", 800, 600);
  c1->Divide(3, 3);
  Int_t i=0;
  c1->cd(++i); hTOFMulti->Draw();
  c1->cd(++i); hTOFHitPat->Draw();
  c1->cd(++i); hTOFPdgPat->Draw("colz");
  c1->cd(++i); hAC1Multi->Draw();
  c1->cd(++i); hAC1HitPat->Draw();
  c1->cd(++i); hAC1PdgPat->Draw("colz");
  c1->cd(++i); hWCMulti->Draw();
  c1->cd(++i); hWCHitPat->Draw();
  c1->cd(++i); hWCPdgPat->Draw("colz");
  c1->Print(fig_path);
}
