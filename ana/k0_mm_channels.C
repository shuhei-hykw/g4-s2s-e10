// k0_mm_channels.C
//
// Per-channel extraction for the K0bar multipole study (2026-07-24):
// same acceptance criterion as k0_scan_missing_mass.C (SDC2-5
// tracking + TOF/AC1/WC + 100 ns coincidence, SDC2 momentum vector
// as the reconstructed pi+ candidate), but
//   - reads a TChain GLOB (segment files, no merge needed),
//   - also records the ACCEPTED pi+ momentum spectrum and the
//     (theta, p) scatter, for the "measured pi+ momentum
//     distribution" requested in the 2026-07-24 discussion.
// Missing mass is reconstructed under the (K-,pi+) hypothesis on a
// 9Be target at rest; -B_Lambda referencing (M(8He)+M_Lambda) is
// done by the plotting python, not here.
//
// Usage:
//   root -l -b -q 'ana/k0_mm_channels.C("root/kpi/seg_*.root", 1.8,
//                                        "cohgsL0")'
//
// Output: reaction-kinematics/geant4-scan/data/k0chan_<tag>.root
// with h_mm (2 MeV bins), h_ppi, h_thp (raw MC entries; absolute
// normalization applied by the calling python).

//______________________________________________________________________________
void k0_mm_channels(const char* glob, double p_beam_gevc,
                    const char* tag)
{
  gSystem->Load("bin/libG4S2S_lib.so");

  TChain* tree = new TChain("g4s2s");
  const int nfiles = tree->Add(glob);
  if(nfiles <= 0){
    std::cerr << "no files match " << glob << std::endl; return;
  }
  std::cout << nfiles << " files chained from " << glob << std::endl;

  std::vector<TParticle>* tofHits = nullptr;
  std::vector<TParticle>* ac1Hits = nullptr;
  std::vector<TParticle>* wcHits  = nullptr;
  std::vector<TParticle>* sdc2Hits = nullptr;
  std::vector<TParticle>* sdc3Hits = nullptr;
  std::vector<TParticle>* sdc4Hits = nullptr;
  std::vector<TParticle>* sdc5Hits = nullptr;
  tree->SetBranchAddress("TOF", &tofHits);
  tree->SetBranchAddress("AC1", &ac1Hits);
  tree->SetBranchAddress("WC",  &wcHits);
  tree->SetBranchAddress("SDC2", &sdc2Hits);
  tree->SetBranchAddress("SDC3", &sdc3Hits);
  tree->SetBranchAddress("SDC4", &sdc4Hits);
  tree->SetBranchAddress("SDC5", &sdc5Hits);

  const double M_KM   = 493.677;    // K- (PDG), MeV
  const double M_PIP  = 139.57039;  // pi+ (PDG), MeV
  const double M_9BE  = 8392.75;    // 9Be nuclear mass, MeV
  const double P_BEAM = p_beam_gevc*1000.; // MeV/c

  const double eBeam = std::sqrt(P_BEAM*P_BEAM + M_KM*M_KM);
  const double eInit = eBeam + M_9BE;

  // Missing-mass zoom: -B_Lambda in [-50, +150] MeV around the
  // 9_Lambda He threshold M(8He)+M_Lambda = 8598.224 MeV, 2 MeV/bin
  // (same axis convention as k0_scan_missing_mass.C).
  const double xmin = 8.548224, xmax = 8.748224;
  const int nbins = 100;
  auto h_mm = new TH1D("h_mm",
    TString::Format("Missing mass, p_{beam}=%.2f GeV/c;M_{X} "
                    "[GeV/c^{2}];raw MC entries", p_beam_gevc),
    nbins, xmin, xmax);
  auto h_ppi = new TH1D("h_ppi",
    "accepted #pi^{+};p_{#pi} [GeV/c];raw MC entries",
    200, 0.6, 1.6);
  auto h_thp = new TH2D("h_thp",
    "accepted #pi^{+};#theta_{#pi} [deg];p_{#pi} [GeV/c]",
    100, 0., 25., 100, 0.6, 1.6);

  int nAccepted = 0;
  const Long64_t nEntries = tree->GetEntries();
  for(Long64_t ie = 0; ie < nEntries; ++ie){
    tree->GetEntry(ie);
    if(!tofHits || !ac1Hits || !wcHits
       || !sdc2Hits || !sdc3Hits || !sdc4Hits || !sdc5Hits) continue;
    if(tofHits->empty() || ac1Hits->empty() || wcHits->empty()) continue;
    if(sdc2Hits->empty() || sdc3Hits->empty()
       || sdc4Hits->empty() || sdc5Hits->empty()) continue;

    std::map<int,double> tofT, ac1T, wcT;
    std::set<int> s3Ids, s4Ids, s5Ids;
    for(const auto& p : *tofHits) tofT[p.GetStatusCode()] = p.T();
    for(const auto& p : *ac1Hits) ac1T[p.GetStatusCode()] = p.T();
    for(const auto& p : *wcHits)  wcT[p.GetStatusCode()]  = p.T();
    for(const auto& p : *sdc3Hits) s3Ids.insert(p.GetStatusCode());
    for(const auto& p : *sdc4Hits) s4Ids.insert(p.GetStatusCode());
    for(const auto& p : *sdc5Hits) s5Ids.insert(p.GetStatusCode());

    for(const auto& p : *sdc2Hits){
      const int trackId = p.GetStatusCode();
      const auto itTof = tofT.find(trackId);
      const auto itAc1 = ac1T.find(trackId);
      const auto itWc  = wcT.find(trackId);
      if(itTof == tofT.end()) continue;
      if(itAc1 == ac1T.end()) continue;
      if(itWc  == wcT.end())  continue;
      if(!s3Ids.count(trackId)) continue;
      if(!s4Ids.count(trackId)) continue;
      if(!s5Ids.count(trackId)) continue;
      const double tmin = std::min({itTof->second, itAc1->second,
                                    itWc->second});
      const double tmax = std::max({itTof->second, itAc1->second,
                                    itWc->second});
      if(tmax - tmin > 100.) continue; // 100 ns coincidence

      const double px = p.Px(), py = p.Py(), pz = p.Pz();
      const double pAbs = std::sqrt(px*px + py*py + pz*pz);
      const double thDeg = std::acos(pz/pAbs)*180./M_PI;
      const double eReco = std::sqrt(pAbs*pAbs + M_PIP*M_PIP);
      const double eX = eInit - eReco;
      const double pxX = -px, pyX = -py, pzX = P_BEAM - pz;
      const double mm2 = eX*eX - (pxX*pxX + pyX*pyX + pzX*pzX);
      if(mm2 <= 0.) continue;
      h_mm->Fill(std::sqrt(mm2)/1000.); // MeV -> GeV
      h_ppi->Fill(pAbs/1000.);
      h_thp->Fill(thDeg, pAbs/1000.);
      ++nAccepted;
      break;
    }
  }
  std::cout << "tag=" << tag
            << "  Entries=" << nEntries
            << "  accepted=" << nAccepted
            << "  frac=" << 100.*nAccepted/nEntries << "%"
            << std::endl;

  TString outName = TString::Format(
    "/gpfs/group/had/sks/Users/shuhei/work/e10/reaction-kinematics/"
    "geant4-scan/data/k0chan_%s.root", tag);
  auto outFile = new TFile(outName, "RECREATE");
  h_mm->Write();
  h_ppi->Write();
  h_thp->Write();
  // record generated-entry count for absolute normalization
  TParameter<Long64_t>("n_generated", nEntries).Write();
  outFile->Close();
  std::cout << "Saved: " << outName << std::endl;
}
