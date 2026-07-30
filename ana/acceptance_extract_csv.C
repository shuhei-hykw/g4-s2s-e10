// acceptance_extract_csv.C
//
// Dump the S-2S (theta,p) acceptance-probability map (PRMPThetaAcc /
// PRMPThetaGen, same ratio as acceptance_e10_kpi.C's hRatio) to a CSV
// grid in the exact format consumed by
// reaction-kinematics/geant4-scan/beam_momentum_scan.py's lookup_eff
// (RegularGridInterpolator over theta_deg,p_gevc -> eff), i.e. the
// same format as data/heff_grid.csv. That file has no surviving
// extraction script (rescued from a volatile scratchpad, see
// analysis-note.md 2026-07-15 09:15) -- this macro is the fresh
// replacement, used first for the p_beam=1.6078 GeV/c rigidity-
// ceiling boundary point (k0-background.md section 7).
//
// Usage:
//   root -l -b -q 'ana/acceptance_extract_csv.C("path/to/run.root",
//                                               "out.csv")'

void acceptance_extract_csv(const char* fname, const char* out_csv)
{
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

  TH2D* hRatio = (TH2D*)hAcc->Clone("hRatio");
  hRatio->Divide(hAcc, hGen, 1., 1., "B");

  // Row order matters: beam_momentum_scan.py reshapes the flat CSV as
  // (n_theta, n_p) via np.unique-sorted axes, so theta must be the
  // OUTER (slow-varying) loop and p the INNER (fast-varying) one --
  // matching data/heff_grid.csv's existing row order (theta fixed
  // across a run of increasing p).
  std::ofstream ofs(out_csv);
  ofs << "theta_deg,p_gevc,eff\n";
  for(Int_t ix=1; ix<=hRatio->GetNbinsX(); ++ix){
    const Double_t theta = hRatio->GetXaxis()->GetBinCenter(ix);
    for(Int_t iy=1; iy<=hRatio->GetNbinsY(); ++iy){
      const Double_t p = hRatio->GetYaxis()->GetBinCenter(iy);
      const Double_t eff = hRatio->GetBinContent(ix, iy);
      ofs << theta << "," << p << "," << eff << "\n";
    }
  }
  ofs.close();

  std::cout << "generated = " << (Long64_t)hGen->GetEntries()
            << ", accepted = " << (Long64_t)hAcc->GetEntries()
            << ", saved " << out_csv << std::endl;
}
