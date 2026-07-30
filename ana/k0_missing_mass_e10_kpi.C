// k0_missing_mass_e10_kpi.C
//
// Predicted missing-mass spectrum for 9Be(K-,pi+)9_Lambda_He,
// 500,000 K-/spill (1 spill = 4.2 s), 7 days of beam:
//  - K0 background (Generator 1001, e10_k0_background_kpi_tree.conf,
//    TREE:1): reconstruct the "fake" missing mass that results when
//    the actual accepted track (pi+/pi-/mu-/e- from the K0bar decay)
//    is mis-identified as the signal pi+ and run through the
//    standard 9Be(K-,pi+) missing-mass formula -- exactly what a
//    real analysis would do, since the trigger cannot distinguish
//    species. Absolute-normalized using the K-p->K0bar n cross
//    section (Conforto et al., Nucl. Phys. B105 (1976) 189-221) and
//    the 9Be target/beam parameters below.
//  - Signal (9_Lambda_He peak): reference Gaussian at
//    M(8He)+m_Lambda-B_Lambda (B_Lambda=0-10 MeV illustrative band,
//    see ana/kinematics_overlay_e10_kpi.C), with a missing-mass
//    resolution derived from docs/p10-update-dcx.pdf (P10-update DCX
//    proposal): SKS achieved 2.5 MeV FWHM excitation-energy
//    resolution with a 3.5 g/cm^2 9Be/6Li target, explicitly stated
//    to be DOMINATED by energy-loss straggling in the target
//    (scaling as sqrt(thickness in g/cm^2); other, thickness-
//    independent contributions are small at this thickness). Our
//    own target (DetSize_E10_9hesigma_kpi: 2.0 cm, Be9 rho=1.85
//    g/cm3) has an areal density of 3.7 g/cm^2 -- close enough to
//    3.5 that the sqrt-scaling correction is a ~3% effect:
//    FWHM = 2.5*sqrt(3.7/3.5) = 2.57 MeV -> sigma = 1.09 MeV. This
//    borrows a real measured value (different spectrometer/reaction,
//    but the same dominant physics -- target straggling) rather than
//    an arbitrary guess; see analysis-note.md 2026-07-13 entry. NOT
//    absolutely normalized (no elementary production cross section
//    is available in this repo -- see
//    hypernuclear-production/scripts/reaction_config.py sigma_elem,
//    which is a placeholder/None). Scaled here to a fixed reference
//    peak height for visual comparison only.
//
// REVISION (2026-07-13): the first version of this macro accepted
// any event with a track hitting TOF+AC1(Cherenkov)+WC and used that
// track's momentum AT WC directly in the missing-mass formula. Two
// problems, both caught by asking "is this really consistent with
// the actual S-2S momentum acceptance, and did you check the real
// Geant4 output for it?":
//  1. The momentum vector AT WC is expressed in the POST-BEND local
//     direction (S-2S bends the whole reference trajectory by a
//     large angle to separate it from the beam line) -- it is NOT
//     directly comparable to the 0-20 deg/1.0-1.6 GeV/c emission-
//     angle acceptance window used everywhere else in this repo, and
//     is NOT the right momentum vector to plug into the missing-mass
//     formula (which assumes a vector defined relative to the
//     ORIGINAL beam axis).
//  2. "TOF+AC1+WC hit" alone does not guarantee the track was
//     actually reconstructable: checking the SAME track id against
//     the SDC2-SDC5 tracking chambers shows only ~45% of
//     TOF+AC1+WC-matched tracks are consistently present in all four
//     (SDC1 has zero overlap for every matched track in this
//     geometry -- K0 decay-in-flight secondaries are simply born
//     downstream of it). The rest take some other, likely
//     non-analyzable path through the large-area trigger counters.
//
// Fix: require the SAME track id in SDC2, SDC3, SDC4, SDC5 (the
// tracking chambers) in addition to TOF+AC1+WC, and use the SDC2 hit
// (upstream of the main D1 bend, so its local direction is close to
// the true lab-frame emission direction/momentum -- verified to peak
// at theta~7.9 deg, matching the expected forward-peaked acceptance
// region, rather than the ~46 deg the raw WC-frame check gave) for
// the missing-mass reconstruction. This is a proxy for a full
// optics-matrix momentum reconstruction (not implemented here), but
// is far more defensible than using the post-bend WC vector.
//
// Track matching: VHitInfo stores the true G4 track id in each hit's
// TParticle::GetStatusCode() for Experiment==10 (see src/VHitInfo.cc).
//
// REVISION (2026-07-13, later): added a 100 ns trigger coincidence
// window (TOF/AC1/WC hit times must all fall within 100 ns of each
// other) -- the hit's global time-since-the-reaction is already
// stored as TParticle::T() (see VHitInfo::VHitInfo, "v" 4-vector),
// so this needed no new Geant4 output, just an additional cut here.
// Also added a zoomed view of the near-threshold ("bound") region.
//
// Usage:
//   root -l -b -q 'ana/k0_missing_mass_e10_kpi.C("run_tree.root")'
//
// See analysis-note.md 2026-07-12/13 entries ("predicted missing-mass
// spectrum").

//______________________________________________________________________________
void k0_missing_mass_e10_kpi(const char* fname)
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

  // Masses (MeV; TParticle momenta are stored in G4's MeV convention
  // -- see VHitInfo::VHitInfo).
  const double M_KM   = 493.677;    // K- (PDG)
  const double M_PIP  = 139.57039;  // pi+ (PDG)
  const double M_9BE  = 8392.75;    // 9Be nuclear mass (AME atomic
                                     // mass - Z*m_e, Z=4)
  // 8He nuclear mass (AME2020 atomic 7483.563 - 2*m_e): the final
  // Lambda hypernucleus 9_Lambda He has an 8HE core (Z=2 by charge
  // conservation), not the Sigma- doorway's 8Li -- threshold bug
  // fixed 2026-07-14 (old M(8Li)+M(Lambda) sat 11.17 MeV too low).
  const double M_8HE  = 7482.541;
  const double M_LAM  = 1115.683;   // Lambda (PDG)
  const double P_BEAM = 1500.;      // MeV/c, K- beam, along +z

  const double eBeam = std::sqrt(P_BEAM*P_BEAM + M_KM*M_KM);
  const double eInit = eBeam + M_9BE;

  const double xmin = 8.2, xmax = 9.6;
  const int nbins = 100;
  auto h_mm = new TH1D("h_mm_k0bg",
    "Predicted missing-mass spectrum, ^{9}Be(K^{-},#pi^{+});"
    "M_{X} [GeV/#it{c}^{2}];events / 7 days / "
    + TString::Format("%.1f MeV", 1000.*(xmax-xmin)/nbins),
    nbins, xmin, xmax);

  // Trigger coincidence window: a real DAQ trigger only fires if the
  // TOF/AC1/WC pulses arrive within a set gate. TParticle::T() already
  // stores each hit's time since the reaction (event start), so this
  // needs no extra Geant4 output -- just a cut on the spread of the
  // three trigger-counter times for the matched track.
  const double coincidenceWindow_ns = 100.;

  int nAccepted = 0;
  int nFailedCoincidence = 0;
  const Long64_t nEntries = tree->GetEntries();
  for(Long64_t ie = 0; ie < nEntries; ++ie){
    tree->GetEntry(ie);
    if(!tofHits || !ac1Hits || !wcHits
       || !sdc2Hits || !sdc3Hits || !sdc4Hits || !sdc5Hits) continue;
    if(tofHits->empty() || ac1Hits->empty() || wcHits->empty()) continue;
    if(sdc2Hits->empty() || sdc3Hits->empty()
       || sdc4Hits->empty() || sdc5Hits->empty()) continue;

    std::map<int,double> tofTime, ac1Time, wcTime;
    std::set<int> s3Ids, s4Ids, s5Ids;
    for(const auto& p : *tofHits) tofTime[p.GetStatusCode()] = p.T();
    for(const auto& p : *ac1Hits) ac1Time[p.GetStatusCode()] = p.T();
    for(const auto& p : *wcHits)  wcTime[p.GetStatusCode()]  = p.T();
    for(const auto& p : *sdc3Hits) s3Ids.insert(p.GetStatusCode());
    for(const auto& p : *sdc4Hits) s4Ids.insert(p.GetStatusCode());
    for(const auto& p : *sdc5Hits) s5Ids.insert(p.GetStatusCode());

    // Require the SAME track in TOF, AC1, WC and every tracking
    // chamber SDC2-SDC5 (SDC1 is excluded -- it has zero overlap for
    // any K0-decay secondary in this geometry, since the decay vertex
    // sits downstream of it). Reconstruct missing mass from the
    // SDC2 hit (upstream of the D1 bend -- see file-level comment).
    for(const auto& p : *sdc2Hits){
      const int trackId = p.GetStatusCode();
      const auto itTof = tofTime.find(trackId);
      const auto itAc1 = ac1Time.find(trackId);
      const auto itWc  = wcTime.find(trackId);
      if(itTof == tofTime.end()) continue;
      if(itAc1 == ac1Time.end()) continue;
      if(itWc  == wcTime.end())  continue;
      if(!s3Ids.count(trackId)) continue;
      if(!s4Ids.count(trackId)) continue;
      if(!s5Ids.count(trackId)) continue;

      const double tMin = std::min({itTof->second, itAc1->second,
                                     itWc->second});
      const double tMax = std::max({itTof->second, itAc1->second,
                                     itWc->second});
      if(tMax - tMin > coincidenceWindow_ns){
        ++nFailedCoincidence;
        continue;
      }

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
      break; // one reconstructed track per event
    }
  }
  std::cout << "Entries=" << nEntries
            << "  accepted(single-track SDC2-5+TOF+AC1+WC, "
            << coincidenceWindow_ns << " ns coincidence)=" << nAccepted
            << "  (rejected by coincidence window: "
            << nFailedCoincidence << ")"
            << std::endl;

  // --- absolute normalization: 500 k K-/spill, 1 spill=4.2 s, 7 days
  const double kaonsPerSpill = 5.0e5;
  const double spillLength_s = 4.2;
  const double runDays       = 7.0;
  const double nSpills = runDays*24.*3600./spillLength_s;
  const double nKaons  = kaonsPerSpill*nSpills;

  // 9Be target: 75x75x20 mm box (DetSize_E10_9hesigma_kpi), beam
  // along the 20 mm (z) dimension; rho=1.85 g/cm3 (MaterialList Be9).
  const double targetThick_cm = 2.0;
  const double rho_gcm3 = 1.85;
  const double molarMass = 9.012182; // g/mol
  const double avogadro = 6.02214076e23;
  const double nBeAreal = rho_gcm3*targetThick_cm*avogadro/molarMass; // /cm^2
  const double Z_Be = 4.;
  const double nProtonAreal = nBeAreal*Z_Be; // quasi-free approx.

  // K-p->K0bar n cross section (Conforto et al., highest measured
  // point 1.320-1.355 GeV/c, extrapolated to p_beam=1.5 GeV/c --
  // likely somewhat lower given the monotonic decrease with
  // momentum; treat as an upper-bound-ish estimate).
  const double sigma_mb = 1.8;
  const double sigma_cm2 = sigma_mb*1.e-27;

  const double nK0barProduced = nKaons*nProtonAreal*sigma_cm2;
  const double mcAcceptedFrac = double(nAccepted)/double(nEntries);
  const double nBgExpected = nK0barProduced*mcAcceptedFrac;

  std::cout << "Beam: " << nKaons << " K- over " << runDays
            << " days\n";
  std::cout << "Target areal protons (quasi-free, Z=4): "
            << nProtonAreal << " /cm^2\n";
  std::cout << "K0bar produced (sigma=" << sigma_mb << " mb): "
            << nK0barProduced << "\n";
  std::cout << "MC accepted fraction (this sample): "
            << mcAcceptedFrac*100. << " %\n";
  std::cout << "==> Expected background events / 7 days: "
            << nBgExpected << std::endl;

  if(h_mm->Integral() > 0.){
    h_mm->Scale(nBgExpected/h_mm->Integral());
  }

  // --- signal reference peak (NOT absolutely normalized -- see
  // file-level comment; scaled to a visual reference height).
  const double bLambdaLo = 0., bLambdaHi = 10.; // MeV, illustrative
  const double mSigCenter = (M_8HE + M_LAM - 0.5*(bLambdaLo+bLambdaHi))
                             /1000.; // GeV
  const double mSigLo = (M_8HE + M_LAM - bLambdaHi)/1000.;
  const double mSigHi = (M_8HE + M_LAM - bLambdaLo)/1000.;
  // Realistic resolution derived from docs/p10-update-dcx.pdf (see
  // file-level comment): 2.5 MeV FWHM at 3.5 g/cm^2, scaled by
  // sqrt(areal density) to our 3.7 g/cm^2 target.
  const double targetArealDensity_gcm2 = 3.7;
  const double referenceArealDensity_gcm2 = 3.5;
  const double referenceFWHM_GeV = 0.0025;
  const double fwhmAssumed_GeV = referenceFWHM_GeV *
    std::sqrt(targetArealDensity_gcm2/referenceArealDensity_gcm2);
  const double sigmaAssumed_GeV = fwhmAssumed_GeV/(2.*std::sqrt(2.*std::log(2.)));
  auto f_sig = new TF1("f_sig", "gaus", xmin, xmax);
  f_sig->SetNpx(2000); // sigma~1.1 MeV is far narrower than the
                        // 12 MeV/bin background histogram -- needs
                        // many points to render as a visible spike
  const double refHeight = h_mm->GetMaximum()>0. ? h_mm->GetMaximum()
                                                  : 1.;
  f_sig->SetParameters(refHeight, mSigCenter, sigmaAssumed_GeV);
  f_sig->SetLineColor(kBlue+1);
  f_sig->SetLineWidth(2);
  f_sig->SetLineStyle(1);

  auto c = new TCanvas("c_k0mm", "c_k0mm", 900, 650);
  h_mm->SetLineColor(kRed+1);
  h_mm->SetLineWidth(2);
  h_mm->SetMinimum(0.);
  h_mm->Draw("hist");
  f_sig->Draw("same");

  auto box = new TBox(mSigLo, 0., mSigHi, h_mm->GetMaximum()*1.15);
  box->SetFillColorAlpha(kBlue, 0.08);
  box->SetLineColor(kBlue);
  box->SetLineStyle(3);
  box->Draw("same");

  auto leg = new TLegend(0.55, 0.72, 0.89, 0.89);
  leg->SetBorderSize(0);
  leg->AddEntry(h_mm, "K0 background (absolute, 7 days)", "l");
  leg->AddEntry(f_sig,
    "^{9}_{#Lambda}He signal (position only, height NOT to scale)",
    "l");
  leg->Draw();

  // --- zoom on the bound region: the Lambda-separation threshold is
  // M(8He)+m_Lambda = 8.598 GeV/c^2 (B_Lambda=0); bound states sit at
  // or just below it (B_Lambda=0-10 MeV here), so zoom to bracket
  // that threshold and see the signal peak against the background's
  // rising edge in detail.
  const double zoomLo = 8.45, zoomHi = 8.75;
  auto c2 = new TCanvas("c_k0mm_zoom", "c_k0mm_zoom", 900, 650);
  h_mm->GetXaxis()->SetRangeUser(zoomLo, zoomHi);
  h_mm->SetTitle(
    "Predicted missing-mass spectrum (bound region), "
    "^{9}Be(K^{-},#pi^{+})");
  h_mm->Draw("hist");
  f_sig->Draw("same");
  box->Draw("same");
  leg->Draw();
  c2->Print("ana/fig/k0_missing_mass_e10_kpi_zoom.pdf");
  c2->Print("ana/fig/k0_missing_mass_e10_kpi_zoom.png");
  h_mm->GetXaxis()->SetRangeUser(xmin, xmax); // restore full range
  h_mm->SetTitle(
    "Predicted missing-mass spectrum, ^{9}Be(K^{-},#pi^{+})");

  c->Print("ana/fig/k0_missing_mass_e10_kpi.pdf");
  c->Print("ana/fig/k0_missing_mass_e10_kpi.png");

  auto outFile = new TFile("ana/fig/k0_missing_mass_e10_kpi_hist.root",
                            "RECREATE");
  h_mm->Write();
  f_sig->Write();
  outFile->Close();
  std::cout << "Saved: ana/fig/k0_missing_mass_e10_kpi.{pdf,png}"
            << " and _hist.root" << std::endl;
}
