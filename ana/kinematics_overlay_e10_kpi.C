// kinematics_overlay_e10_kpi.C
//
// Overlay the pi+ kinematics locus for 9Be(K-,pi+)9_Lambda_He
// (p_beam = 1.5 GeV/c) on the S-2S acceptance-probability (theta,p)
// map, so the physical pi+ production locus can be compared directly
// against the region where the spectrometer actually has acceptance.
//
// PHYSICS PICTURE (corrected twice after review -- see analysis-note.md
// 2026-07-12 entries; this is the final version):
//
// The Sigma- is an OFF-SHELL doorway/intermediate configuration, not
// a real asymptotic particle -- the actual observed final state of
// this reaction is pi+ + 9_Lambda_He directly (the Sigma- component
// only mixes virtually into the physical eigenstate via Lambda-Sigma
// coupling in the nuclear medium; it never appears as a free
// propagating particle). Energy-momentum conservation therefore
// applies between the INITIAL state (K- beam + 9Be target at rest)
// and the FINAL asymptotic state (pi+ + recoiling 9_Lambda_He, mass
// M(8HE) + m_Lambda - B_Lambda; 8He core by charge conservation,
// corrected 2026-07-14 from the earlier wrong 8Li core) -- m_Sigma
// does not enter the
// asymptotic kinematics at all. This is the PRIMARY curve below
// ("Lambda final state").
//
// Two earlier, less complete pictures are kept for reference/contrast
// only (both shown in grey, not the physics prediction):
//  - "Sigma- on-shell (hypothetical)": full 2-body kinematics but
//    swapping in m_Sigma instead of m_Lambda, as if Sigma- were a
//    real final-state particle. ~85-95 MeV/c LOWER than the Lambda
//    curve (m_Sigma - m_Lambda = 81.8 MeV shows up almost directly
//    as a momentum shift). Wrong for an off-shell doorway, but useful
//    to see how much the on-shell-vs-virtual distinction matters.
//  - "Quasi-free single-nucleon proxy": target = bound proton alone
//    (m_p-S_p), recoil = hyperon alone, core treated as a
//    non-recoiling spectator. This is what actually feeds the
//    Neff/form-factor calculation elsewhere in this repo
//    (momentum_transfer.q_lab, ReactionConfig) -- a model input for
//    the cross-section magnitude/shape, not a kinematics prediction.
//
// The Lambda-final-state band uses B_Lambda = 0 to 10 MeV (typical
// p-shell Lambda-hypernuclear binding scale; no specific measured
// value for 9_Lambda_He, so illustrative).
//
// Output: fig/kinematics_overlay_e10_kpi.pdf
//
// Usage:  root -l -b -q 'ana/kinematics_overlay_e10_kpi.C("run.root")'

void kinematics_overlay_e10_kpi(const char* fname)
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
  TH2D* hEff = (TH2D*)hAcc->Clone("hEff");
  hEff->Divide(hGen);
  hEff->SetMinimum(0.); hEff->SetMaximum(1.);
  hEff->SetTitle(
    "Acceptance probability with"
    " ^{9}Be(K^{-},#pi^{+})^{9}_{#Lambda}He kinematics"
    ";Emission angle #theta [deg];Momentum p [GeV/#it{c}]");

  const int nK = 41;
  double thetaDeg[nK] = {
    0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
    5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,
    10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5,
    15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0
  };

  // PRIMARY: final state pi+ + 9_Lambda_He directly (Sigma- virtual,
  // off-shell doorway only). Full 2-body kinematics, 9Be target at
  // rest, (8He+Lambda) recoiling as one body. B_Lambda = 0-10 MeV.
  double pLam_B0[nK] = {
    1.3655, 1.3655, 1.3655, 1.3654, 1.3654, 1.3653, 1.3652, 1.3651,
    1.3649, 1.3648, 1.3646, 1.3644, 1.3642, 1.3639, 1.3637, 1.3634,
    1.3631, 1.3628, 1.3625, 1.3622, 1.3618, 1.3614, 1.3611, 1.3606,
    1.3602, 1.3598, 1.3593, 1.3588, 1.3583, 1.3578, 1.3573, 1.3567,
    1.3562, 1.3556, 1.3550, 1.3543, 1.3537, 1.3531, 1.3524, 1.3517,
    1.3510
  };
  double pLam_B10[nK] = {
    1.3757, 1.3757, 1.3757, 1.3756, 1.3756, 1.3755, 1.3754, 1.3753,
    1.3751, 1.3750, 1.3748, 1.3746, 1.3744, 1.3741, 1.3739, 1.3736,
    1.3733, 1.3730, 1.3727, 1.3724, 1.3720, 1.3716, 1.3712, 1.3708,
    1.3704, 1.3699, 1.3695, 1.3690, 1.3685, 1.3680, 1.3674, 1.3669,
    1.3663, 1.3657, 1.3651, 1.3645, 1.3638, 1.3632, 1.3625, 1.3618,
    1.3611
  };

  // REFERENCE ONLY: hypothetical on-shell Sigma- (full 2-body,
  // 9Be target, 8Li+Sigma- recoil), B_Sigma=0.
  double pSigOnShell_B0[nK] = {
    1.2931, 1.2931, 1.2931, 1.2930, 1.2930, 1.2929, 1.2928, 1.2927,
    1.2925, 1.2924, 1.2922, 1.2920, 1.2918, 1.2916, 1.2914, 1.2911,
    1.2909, 1.2906, 1.2903, 1.2900, 1.2896, 1.2893, 1.2889, 1.2885,
    1.2881, 1.2877, 1.2872, 1.2868, 1.2863, 1.2858, 1.2853, 1.2848,
    1.2842, 1.2837, 1.2831, 1.2825, 1.2819, 1.2813, 1.2807, 1.2800,
    1.2794
  };

  // REFERENCE ONLY: quasi-free single-nucleon proxy, B_Sigma=0.
  double pQF_B0[nK] = {
    1.2744, 1.2743, 1.2741, 1.2737, 1.2732, 1.2726, 1.2717, 1.2708,
    1.2697, 1.2684, 1.2671, 1.2655, 1.2639, 1.2621, 1.2601, 1.2580,
    1.2558, 1.2534, 1.2510, 1.2483, 1.2456, 1.2427, 1.2397, 1.2366,
    1.2334, 1.2300, 1.2265, 1.2229, 1.2192, 1.2154, 1.2115, 1.2075,
    1.2034, 1.1992, 1.1949, 1.1905, 1.1860, 1.1814, 1.1767, 1.1720,
    1.1671
  };

  auto lamB0 = new TGraph(nK, thetaDeg, pLam_B0);
  lamB0->SetLineColor(kRed+1);
  lamB0->SetLineWidth(3);
  lamB0->SetLineStyle(2);

  auto lamB10 = new TGraph(nK, thetaDeg, pLam_B10);
  lamB10->SetLineColor(kRed+1);
  lamB10->SetLineWidth(3);
  lamB10->SetLineStyle(1);

  // Shaded band: the primary (Lambda final state) prediction.
  auto band = new TGraph(2*nK);
  for(int i = 0; i < nK; ++i)
    band->SetPoint(i, thetaDeg[i], pLam_B10[i]);
  for(int i = 0; i < nK; ++i)
    band->SetPoint(nK+i, thetaDeg[nK-1-i], pLam_B0[nK-1-i]);
  band->SetFillColorAlpha(kRed+1, 0.3);
  band->SetLineColor(kRed+1);

  auto sigRef = new TGraph(nK, thetaDeg, pSigOnShell_B0);
  sigRef->SetLineColor(kGray+2);
  sigRef->SetLineWidth(2);
  sigRef->SetLineStyle(2);

  auto qfRef = new TGraph(nK, thetaDeg, pQF_B0);
  qfRef->SetLineColor(kGray+2);
  qfRef->SetLineWidth(2);
  qfRef->SetLineStyle(3);

  auto c = new TCanvas("c", "c", 900, 700);
  c->SetRightMargin(0.15);
  hEff->Draw("COLZ");
  band->Draw("F SAME");
  lamB0->Draw("L SAME");
  lamB10->Draw("L SAME");
  sigRef->Draw("L SAME");
  qfRef->Draw("L SAME");

  auto leg = new TLegend(0.13, 0.13, 0.72, 0.33);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(lamB0,
    "#pi^{+}+^{9}_{#Lambda}He final state (#Sigma^{-} off-shell"
    " doorway), B_{#Lambda}=0", "l");
  leg->AddEntry(lamB10,
    "same, B_{#Lambda}=10 MeV (illustrative bound side)", "l");
  leg->AddEntry(sigRef,
    "ref. only: hypothetical on-shell #Sigma^{-}, B_{#Sigma}=0", "l");
  leg->AddEntry(qfRef,
    "ref. only: quasi-free single-nucleon proxy, B_{#Sigma}=0", "l");
  leg->Draw();

  c->Print("ana/fig/kinematics_overlay_e10_kpi.pdf");
  c->Print("ana/fig/kinematics_overlay_e10_kpi.png");
  std::cout << "Saved: ana/fig/kinematics_overlay_e10_kpi.{pdf,png}"
            << std::endl;
}
