// k0_background_e10_kpi.C
//
// Estimate the K0 background: K- + p (quasi-free, bound proton)
// -> K0bar n at p_beam=1.5 GeV/c, using the MEASURED c.m. angular
// distribution (Conforto et al., Nucl. Phys. B105 (1976) 189-221,
// Rutherford-Imperial College collaboration, K-p reactions
// 0.960-1.355 GeV/c -- see refs/1-s2.0-0550321376902625-main.pdf),
// followed by K0_S -> pi+ pi- decay (isotropic in the K0 rest frame).
//
// REVISION (2026-07-12): the first version of this macro sampled K0
// production only via the "forward c.m." lab-kinematics branch at a
// flat-in-theta_lab weighting -- i.e. it implicitly assumed the K0
// c.m. angular distribution is forward-peaked. Conforto's actual
// measured Legendre coefficients show K-p->K0bar n is BACKWARD-peaked
// in the K0bar c.m. angle throughout their whole momentum range
// (backward/forward dsigma/dcosTheta* ratio ~1.5-5.8; ~3.5 at their
// highest point, 1.355 GeV/c) -- i.e. the RECOILING NEUTRON is
// preferentially forward, not the K0bar. This macro now samples
// cosTheta*_K0 in the c.m. frame from the actual Legendre expansion
// (using the 1.355 GeV/c coefficients -- the closest measured point
// to our 1.5 GeV/c target; Conforto's paper does not extend high
// enough to reach 1.5 GeV/c exactly, hence "usually a bit short",
// per the collaboration's own characterization) and boosts properly
// from c.m. to lab, instead of assuming a lab-angle range a priori.
//
// This is a kinematics+acceptance overlap fraction, not an absolute
// rate: converting to a S/B ratio would additionally need the
// absolute hypernuclear-production cross section, which is out of
// scope for this repo (see hypernuclear-production/README).
//
// K0_S cw = 2.68 cm; at p_K0 ~1.2-1.6 GeV/c (gamma*beta ~2.4-3.2) the
// decay length is only ~6-9 cm, so K0_S decays essentially at the
// target -- the production kinematics (not a spread-out decay
// position) sets the decay-pi+ kinematics.
//
// Usage: root -l -b -q 'ana/k0_background_e10_kpi.C("run.root")'

//______________________________________________________________________________
double LegendreP(int l, double x)
{
  if(l == 0) return 1.;
  if(l == 1) return x;
  double p0 = 1., p1 = x, pl = 0.;
  for(int n = 2; n <= l; ++n){
    pl = ((2*n-1)*x*p1 - (n-1)*p0)/n;
    p0 = p1; p1 = pl;
  }
  return p1;
}

//______________________________________________________________________________
void k0_background_e10_kpi(const char* fname)
{
  gRandom->SetSeed(12345);

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
  TH2D* hEff = (TH2D*)hAcc->Clone("hEff_k0bg");
  hEff->Divide(hGen);

  const double M_PI_C = 139.57039;
  const double M_K0   = 497.611;
  const double M_KM    = 493.677;
  const double M_N    = 939.5654205;
  const double M_P_EFF = 921.3850010499991; // bound proton, quasi-free
  const double EstarPi = M_K0/2.0;
  const double PstarPi = TMath::Sqrt(EstarPi*EstarPi - M_PI_C*M_PI_C);
  const double P_BEAM  = 1500.0; // MeV/c

  // Conforto et al. Legendre A_l coefficients for K-p->K0bar n at
  // p_K-=1.355 GeV/c (their highest point; ~11% short of our 1.5
  // GeV/c beam, used as the closest available proxy for the c.m.
  // angular SHAPE -- only relative weights matter here, so the
  // overall normalization of A_l is irrelevant).
  // cos(theta*) convention: +1 = K0bar forward (beam direction).
  const int nA = 10;
  double Al[nA] = {
    0.156, -0.042, 0.198, -0.181, 0.100, -0.160, 0.210, -0.097,
    -0.019, -0.023
  };

  auto dSigmaDCosTheta = [&](double x){
    double s = 0.;
    for(int l = 0; l < nA; ++l) s += Al[l]*LegendreP(l, x);
    return std::max(s, 0.); // clamp small negative dips from data noise
  };
  // Find envelope max for rejection sampling.
  double fMax = 0.;
  for(int i = 0; i <= 2000; ++i){
    double x = -1. + 2.*i/2000.;
    fMax = std::max(fMax, dSigmaDCosTheta(x));
  }
  fMax *= 1.05; // safety margin

  // c.m. kinematics for K-(beam) + p(bound, at rest) -> K0bar + n.
  double eBeam = TMath::Sqrt(P_BEAM*P_BEAM + M_KM*M_KM);
  TLorentzVector beamLab(0., 0., P_BEAM, eBeam);
  TLorentzVector targetLab(0., 0., 0., M_P_EFF);
  TLorentzVector cmSystem = beamLab + targetLab;
  TVector3 cmBoost = cmSystem.BoostVector();
  double sqrtS = cmSystem.M();
  // c.m. breakup momentum for K0bar(497.611)+n(939.565).
  double eStarK0 = (sqrtS*sqrtS + M_K0*M_K0 - M_N*M_N)/(2*sqrtS);
  double pStarK0 = TMath::Sqrt(std::max(0., eStarK0*eStarK0 - M_K0*M_K0));
  std::cout << "sqrt(s) = " << sqrtS << " MeV, p*_K0bar = "
            << pStarK0 << " MeV/c" << std::endl;

  const int nGen = 2000000; // K0 production events (c.m. sampled)
  const int nDecayPerK0 = 20; // decay MC samples per K0

  auto h_gen = new TH1D("h_gen_k0", "K0 production;cos#theta*_{K0bar}",
                         50, -1., 1.);
  auto h_acc = new TH1D("h_acc_k0",
    "K0_{S}#rightarrow#pi^{+}#pi^{-} decay-#pi^{+} acceptance"
    " weight;cos#theta*_{K0bar} (c.m., + = K0bar fwd)",
    50, -1., 1.);
  auto h_thlab = new TH1D("h_thlab_k0",
    "K0bar lab angle (production);#theta_{K0bar}^{lab} [deg]",
    60, 0., 30.);

  double sumWeight = 0., sumAccWeight = 0.;
  int nSampled = 0;
  while(nSampled < nGen){
    double x = gRandom->Uniform(-1., 1.);
    double y = gRandom->Uniform(0., fMax);
    if(y > dSigmaDCosTheta(x)) continue;
    ++nSampled;
    h_gen->Fill(x);

    double thStarK0 = TMath::ACos(x);
    double phiStarK0 = gRandom->Uniform(0., 2*TMath::Pi());
    TLorentzVector k0Cm(
      pStarK0*TMath::Sin(thStarK0)*TMath::Cos(phiStarK0),
      pStarK0*TMath::Sin(thStarK0)*TMath::Sin(phiStarK0),
      pStarK0*TMath::Cos(thStarK0),
      eStarK0);
    TLorentzVector k0Lab(k0Cm);
    k0Lab.Boost(cmBoost);
    h_thlab->Fill(k0Lab.Theta()/TMath::DegToRad());
    TVector3 k0BoostVec = k0Lab.BoostVector();

    double accHere = 0.;
    for(int id = 0; id < nDecayPerK0; ++id){
      double cosThStar = gRandom->Uniform(-1., 1.);
      double thStar = TMath::ACos(cosThStar);
      double phiStar = gRandom->Uniform(0., 2*TMath::Pi());
      TLorentzVector piRest(
        PstarPi*TMath::Sin(thStar)*TMath::Cos(phiStar),
        PstarPi*TMath::Sin(thStar)*TMath::Sin(phiStar),
        PstarPi*TMath::Cos(thStar),
        EstarPi);
      TLorentzVector piLab(piRest);
      piLab.Boost(k0BoostVec);

      double pPi = piLab.P()/1000.0; // GeV/c
      double thPi = piLab.Theta()/TMath::DegToRad();

      double eff = 0.;
      if(thPi >= 0. && thPi <= 20. && pPi >= 1.0 && pPi <= 1.6){
        int ix = hEff->GetXaxis()->FindBin(thPi);
        int iy = hEff->GetYaxis()->FindBin(pPi);
        eff = hEff->GetBinContent(ix, iy);
      }
      accHere += eff;
    }
    accHere /= nDecayPerK0;
    h_acc->Fill(x, accHere);
    sumWeight += 1.;
    sumAccWeight += accHere;
  }

  double overallFrac = sumAccWeight/sumWeight;
  const double BR_K0S_pipi = 0.692;  // PDG K0_S -> pi+pi-
  const double frac_K0S = 0.50;      // K0bar = 50% K0_S + 50% K0_L

  std::cout << "\n=== K0 background estimate"
               " (Conforto c.m. angular distribution, 1.355 GeV/c"
               " Legendre coefficients, applied at p_beam=1.5 GeV/c) ===\n";
  std::cout << "Overall fraction of K0_S->pi+pi- decays with pi+"
            << " accepted by S-2S: " << overallFrac << std::endl;
  std::cout << "Including K0bar->K0_S branching (50%) and K0_S->pi+pi-"
            << " branching (69.2%):\n  fraction of ALL produced K0bar"
            << " yielding an accepted background pi+ via this channel: "
            << overallFrac*frac_K0S*BR_K0S_pipi << std::endl;

  auto c = new TCanvas("c", "c", 1200, 500);
  c->Divide(3, 1);
  c->cd(1);
  h_gen->SetLineColor(kBlack);
  h_gen->Draw("hist");
  c->cd(2);
  h_thlab->SetLineColor(kBlack);
  h_thlab->Draw("hist");
  c->cd(3);
  auto h_ratio = (TH1D*)h_acc->Clone("h_ratio_k0");
  h_ratio->Divide(h_gen);
  h_ratio->SetTitle("decay-#pi^{+} S-2S acceptance fraction vs"
                     " K0bar production angle;cos#theta*_{K0bar}"
                     " (c.m.);fraction accepted");
  h_ratio->SetMinimum(0.);
  h_ratio->SetLineColor(kBlue+2);
  h_ratio->SetLineWidth(2);
  h_ratio->Draw("hist");
  c->Print("ana/fig/k0_background_e10_kpi.pdf");
  c->Print("ana/fig/k0_background_e10_kpi.png");
  std::cout << "\nSaved: ana/fig/k0_background_e10_kpi.{pdf,png}"
            << std::endl;
}
