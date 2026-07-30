// k0_scan_missing_mass.C
//
// Generalized version of k0_missing_mass_e10_kpi.C for the beam-
// momentum scan (see estimate.md): takes p_beam as an argument
// (instead of hardcoding 1.5 GeV/c) since both the beam energy and
// the missing-mass formula depend on it. Same acceptance criterion
// (SDC2-5 tracking + TOF/AC1(Cherenkov)/WC + 100 ns coincidence,
// using the SDC2 momentum vector as the reconstructed candidate) as
// the validated p_beam=1.5 GeV/c analysis.
//
// Usage:
//   root -l -b -q 'ana/k0_scan_missing_mass.C("run_tree.root", 1.5)'
//
// Returns/saves a histogram of raw (unnormalized) missing-mass
// entries; absolute normalization (beam flux, target, Bricman cross
// section) is applied by the calling Python script.

//______________________________________________________________________________
void k0_scan_missing_mass(const char* fname, double p_beam_gevc,
                           const char* tag = "")
{
  gSystem->Load("bin/libG4S2S_lib.so");
  gStyle->SetOptStat(0);

  TFile* file = TFile::Open(fname);
  if(!file || file->IsZombie()){
    std::cerr << "cannot open " << fname << std::endl; return;
  }
  TTree* tree = (TTree*)file->Get("g4s2s");
  if(!tree){
    std::cerr << "tree 'g4s2s' not found in " << fname << std::endl;
    return;
  }

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

  // Zoomed to the region requested: -B_Lambda = -50 to +150 MeV
  // (i.e. B_Lambda = +50 to -150 MeV) around the B_Lambda=0 threshold
  // M(8HE)+m_Lambda = 7482.541 + 1115.683 = 8598.224 MeV. The final
  // Lambda hypernucleus 9_Lambda He has an 8HE core (charge
  // conservation: 9Be(Z=4) + K-(-1) - pi+(+1) leaves Z=2), NOT the
  // 8Li core of the one-step Sigma- doorway -- earlier versions of
  // this macro wrongly used M(8Li)+m_Lambda = 8587.053 MeV, placing
  // the threshold 11.17 MeV too low (bug fixed 2026-07-14, see
  // analysis-note.md). -B_Lambda=-50 -> M_X=8548.224 MeV;
  // -B_Lambda=+150 -> M_X=8748.224 MeV; 2 MeV/bin.
  const double xmin = 8.548224, xmax = 8.748224;
  const int nbins = 100;
  TString hname = TString::Format("h_mm_p%.0f",
                                   std::round(p_beam_gevc*100.));
  auto h_mm = new TH1D(hname,
    TString::Format("Missing mass, p_{beam}=%.2f GeV/c;M_{X} "
                     "[GeV/c^{2}];raw MC entries", p_beam_gevc),
    nbins, xmin, xmax);

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
      const double eReco = std::sqrt(px*px + py*py + pz*pz
                                      + M_PIP*M_PIP);
      const double eX = eInit - eReco;
      const double pxX = -px, pyX = -py, pzX = P_BEAM - pz;
      const double mm2 = eX*eX - (pxX*pxX + pyX*pyX + pzX*pzX);
      if(mm2 <= 0.) continue;
      const double mm = std::sqrt(mm2);
      h_mm->Fill(mm/1000.); // MeV -> GeV
      ++nAccepted;
      break;
    }
  }
  std::cout << "p_beam=" << p_beam_gevc
            << "  Entries=" << nEntries
            << "  accepted=" << nAccepted
            << "  frac=" << 100.*nAccepted/nEntries << "%"
            << std::endl;

  TString tagSuffix = (TString(tag).Length() > 0)
    ? TString("_") + tag : TString("");
  TString outName = TString::Format(
    "/tmp/claude-14069/-gpfs-group-had-sks-Users-shuhei-work-e10"
    "/2d25f105-e83a-421a-9df7-4aa77f196d17/scratchpad/"
    "k0scan_mm_zoom_p%.0f%s.root", std::round(p_beam_gevc*100.),
    tagSuffix.Data());
  auto outFile = new TFile(outName, "RECREATE");
  h_mm->Write();
  outFile->Close();
  std::cout << "Saved: " << outName << std::endl;
}
