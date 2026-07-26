// k0_pip_thetap.C
//
// Dump the TRUTH (theta, p) distribution of the K0_S -> pi+ pi- decay
// pi+ for a coherent (K-,K0bar) channel, at GENERATION level (no
// spectrometer cut), so the same pi+ source can be folded through
// EITHER the S-2S or the SKS acceptance map in python (2026-07-24:
// the angular acceptance is just a uniform-pi+ map, independent of
// the reaction; no need to port the K0 generator into g4-sks-e05).
//
// Reads the primary K0bar 4-vector (PRM branch) from our Geant4 K0
// production, re-decays it as K0_S -> pi+ pi- isotropically in the
// K0_S rest frame (N_DECAY samples/event), boosts to the lab, and
// fills a fine 2D histogram h_pip(theta_pi [deg], p_pi [GeV]). The
// K0_S/K0_L split and BR(pi+pi-) are applied later as an overall
// factor in python (here every K0bar is decayed to pi+ pi-).
//
// Usage:
//   root -l -b -q 'ana/k0_pip_thetap.C("root/kpi/local_coh_gs_L2_p180.root",
//                                       "cohgsL2_p180")'

//______________________________________________________________________________
void k0_pip_thetap(const char* fname, const char* tag)
{
  gSystem->Load("bin/libG4S2S_lib.so");
  const int N_DECAY = 5;               // decays sampled per K0bar
  const double M_K0  = 497.611;        // MeV (K0_S ~ K0 mass)
  const double M_PIP = 139.57039;      // MeV

  TFile* file = TFile::Open(fname);
  if(!file || file->IsZombie()){
    std::cerr << "cannot open " << fname << std::endl; return;
  }
  TTree* tree = (TTree*)file->Get("g4s2s");
  std::vector<TParticle>* prm = nullptr;
  tree->SetBranchAddress("PRM", &prm);

  // theta 0-30 deg (covers both S-2S and SKS), p 0-2 GeV
  auto h = new TH2D("h_pip",
    "K0S decay pi+ (truth);#theta_{#pi} [deg];p_{#pi} [GeV/c]",
    300, 0., 30., 400, 0., 2.0);

  const double pstar = 0.5 * std::sqrt(M_K0*M_K0 - 4.*M_PIP*M_PIP);
  TRandom3 rng(12345);
  const Long64_t nEntries = tree->GetEntries();
  for(Long64_t ie = 0; ie < nEntries; ++ie){
    tree->GetEntry(ie);
    if(!prm || prm->empty()) continue;
    const auto& k0 = prm->at(0);
    TLorentzVector k0lv(k0.Px()*1000., k0.Py()*1000., k0.Pz()*1000.,
                        k0.Energy()*1000.); // GeV->MeV
    TVector3 boost = k0lv.BoostVector();
    for(int d = 0; d < N_DECAY; ++d){
      const double cth = 2.*rng.Rndm() - 1.;
      const double sth = std::sqrt(1. - cth*cth);
      const double phi = 2.*M_PI*rng.Rndm();
      TLorentzVector pip(pstar*sth*std::cos(phi),
                         pstar*sth*std::sin(phi),
                         pstar*cth,
                         std::sqrt(pstar*pstar + M_PIP*M_PIP));
      pip.Boost(boost);
      const double p = pip.Vect().Mag();
      const double th = pip.Vect().Theta()*180./M_PI;
      h->Fill(th, p/1000.);            // MeV->GeV
    }
  }
  h->Scale(1.0 / N_DECAY);             // per K0bar

  TString outName = TString::Format(
    "/gpfs/group/had/sks/Users/shuhei/work/e10/reaction-kinematics/"
    "geant4-scan/data/k0pip_%s.root", tag);
  auto out = new TFile(outName, "RECREATE");
  h->Write();
  TParameter<Long64_t>("n_generated", nEntries).Write();
  out->Close();
  std::cout << tag << ": " << nEntries << " K0bar -> h_pip saved "
            << outName << std::endl;
}
