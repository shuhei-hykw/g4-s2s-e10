// -*- C++ -*-

#include "S2SPrimaryGeneratorAction.hh"

#include <algorithm>
#include <array>
#include <cmath>

#include <G4Event.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4LorentzVector.hh>
#include <G4ThreeVector.hh>
#include <Randomize.hh>
#include <G4RandomDirection.hh> // E63

#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGenPhaseSpace.h>

#include "BeamMan.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "FuncName.hh"
#include "S2SAnaManager.hh"
#include "DetSizeMan.hh"
#include "FdComplex.hh"
#include "FermiMotion.hh"
#include "CMSMomentum.hh"
#include "GeneratorParticleBranches.hh"


namespace
{
  using CLHEP::mm;
  using CLHEP::deg;
  using CLHEP::radian;
  using CLHEP::GeV;
  const auto& beamMan = BeamMan::GetInstance();
  const auto& confMan = ConfMan::GetInstance();
  const auto& geomMan = DCGeomMan::GetInstance();
  const auto& sizeMan = DetSizeMan::GetInstance();
  auto& anaMan = S2SAnaManager::GetInstance();
  const auto particleTable = G4ParticleTable::GetParticleTable();
  const auto& zK18Target = geomMan.LocalZ("K18Target");
  BeamInfo beam;
  namespace GenBranch = GeneratorParticleBranches;

  G4LorentzVector ToG4Lorentz(const TLorentzVector& lv)
  {
    return G4LorentzVector(lv.Px()*CLHEP::GeV,
                           lv.Py()*CLHEP::GeV,
                           lv.Pz()*CLHEP::GeV,
                           lv.E()*CLHEP::GeV);
  }

  G4LorentzVector ToG4Lorentz(const TVector3& vec, G4double energy_mev)
  {
    return G4LorentzVector(vec.X()*CLHEP::GeV,
                           vec.Y()*CLHEP::GeV,
                           vec.Z()*CLHEP::GeV,
                           energy_mev);
  }

  G4double TwoBodyMomentum(G4double parent_mass,
                           G4double child1_mass,
                           G4double child2_mass)
  {
    const G4double term1 = parent_mass*parent_mass
      - (child1_mass + child2_mass)*(child1_mass + child2_mass);
    const G4double term2 = parent_mass*parent_mass
      - (child1_mass - child2_mass)*(child1_mass - child2_mass);
    if(term1 <= 0. || term2 <= 0.)
      return 0.;
    return 0.5*std::sqrt(term1*term2)/parent_mass;
  }

  inline G4int ToGeneratorId(GenBranch::ParticleId id)
  {
    return static_cast<G4int>(id);
  }

  // GenPID (optional conf key): 1:K+ 2:K- 3:pi+ 4:pi- 5:p 6:e- 7:mu- 8:xi-
  // Returns "" for 0/unrecognized so callers fall back to their
  // historical per-experiment default species.
  G4String GenPidToName(G4int gen_pid)
  {
    switch(gen_pid){
    case 1: return "kaon+";
    case 2: return "kaon-";
    case 3: return "pi+";
    case 4: return "pi-";
    case 5: return "proton";
    case 6: return "e-";
    case 7: return "mu-";
    case 8: return "xi-";
    default: return "";
    }
  }

  // Legendre polynomial, used for the K0bar c.m.-frame angular
  // distribution in GenerateK0Production (Conforto et al., Nucl.
  // Phys. B105 (1976) 189-221, K-p->K0bar n Legendre A_l
  // coefficients at 1.355 GeV/c -- see analysis-note.md 2026-07-12
  // entries).
  G4double LegendreP_K0(G4int l, G4double x)
  {
    if(l == 0) return 1.;
    if(l == 1) return x;
    G4double p0 = 1., p1 = x, pl = 0.;
    for(G4int n = 2; n <= l; ++n){
      pl = ((2*n-1)*x*p1 - (n-1)*p0)/n;
      p0 = p1; p1 = pl;
    }
    return p1;
  }

  // Uniform-in-volume vertex sampling within a box-shaped target
  // (half_size = full extents/2, matching the
  // sizeMan.GetSize("Target")*mm/2 convention used throughout this
  // file). Added 2026-07-19: Generator 2/1001/1002 previously placed
  // every vertex at the single fixed point target_pos even though
  // the target volume is a real (now 50x50x20mm, sized to fit the
  // shared E90/HypTPC target holder's phi80mm bore -- see
  // param/DSIZE/DetSize_E10_9hesigma_kpi) G4Box: every generated
  // particle therefore traversed the same path length through target
  // material, missing the extra straggling/multiple-scattering
  // spread a real uniform-in-volume vertex distribution produces.
  // See k0-background.md section 7/8 and analysis-note.md 2026-07-19
  // entries.
  G4ThreeVector SampleUniformTargetVertex(const G4ThreeVector& target_pos,
                                          const G4ThreeVector& half_size)
  {
    const G4double x0 = G4RandFlat::shoot(-half_size.x(), half_size.x());
    const G4double y0 = G4RandFlat::shoot(-half_size.y(), half_size.y());
    const G4double z0 = G4RandFlat::shoot(-half_size.z(), half_size.z());
    return target_pos + G4ThreeVector(x0, y0, z0);
  }

  void RecordGeneratedParticle(const G4String& branch,
                               GenBranch::ParticleId mother_id,
                               G4int pdg,
                               const TLorentzVector& lv,
                               const G4LorentzVector& vertex)
  {
    anaMan.SetGeneratedParticle(branch,
                                ToGeneratorId(mother_id),
                                pdg,
                                ToG4Lorentz(lv),
                                vertex);
  }

  void RecordGeneratedParticle(const G4String& branch,
                               GenBranch::ParticleId mother_id,
                               G4int pdg,
                               const TVector3& vec,
                               G4double energy_mev,
                               const G4LorentzVector& vertex)
  {
    anaMan.SetGeneratedParticle(branch,
                                ToGeneratorId(mother_id),
                                pdg,
                                ToG4Lorentz(vec, energy_mev),
                                vertex);
  }
  void RecordGeneratedParticle(const G4String& branch,
                               GenBranch::ParticleId mother_id,
                               G4int pdg,
                               const G4LorentzVector& lv,
                               const G4LorentzVector& vertex)
  {
    anaMan.SetGeneratedParticle(branch,
                                ToGeneratorId(mother_id),
                                pdg,
                                lv,
                                vertex);
  }

}

//_____________________________________________________________________________
S2SPrimaryGeneratorAction::S2SPrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(),
    m_particleGun(nullptr),
    m_generator(confMan.Get<G4int>("Generator"))
{
  auto igene = confMan.Get<G4int>("Generator");
  if(igene==4 ||igene==7501 ){
    auto ifsK18name = confMan.Get<G4String>("K18ROOT");
    profileK18 = new TFile(ifsK18name);
    if(!profileK18){G4cout<< ifsK18name << " is not found." <<G4endl;}
    k18track = (TTree*)profileK18 ->Get("k18track");
    k18track ->SetBranchAddress("p_3rd",&p_3rd);
    k18track ->SetBranchAddress("xtgtK18",&xtgt);
    k18track ->SetBranchAddress("ytgtK18",&ytgt);
    k18track ->SetBranchAddress("utgtK18",&utgt);
    k18track ->SetBranchAddress("vtgtK18",&vtgt);
    k18track ->SetBranchAddress("trigflag",&trigflag);
    k18track ->SetBranchAddress("CBtof0",&CBtof0);
    k18track ->SetBranchAddress("ntK18",&ntK18);
    k18track ->SetBranchAddress("chisqrK18",&chisqrK18);
  }

  if(igene==7501 ){
    auto ifsTheoname = confMan.Get<G4String>("TheoCalc");
    gr = new TGraph(ifsTheoname, "%lg %*lg %*lg %*lg %lg");
    for(Int_t i=0;i<gr->GetN();++i){
      G4double x,y;
      gr-> GetPoint(i,x,y);
      Sum += y;
    }
  }
}

//_____________________________________________________________________________
S2SPrimaryGeneratorAction::~S2SPrimaryGeneratorAction()
{
  if(m_particleGun) delete m_particleGun;
}

//_____________________________________________________________________________
void
S2SPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  beam = beamMan.Get();

  if(m_particleGun) delete m_particleGun;
  switch(m_generator){
  case 0: GenerateDemo(anEvent); break;
  case 1: GenerateMonochromeBeam(anEvent); break;
  case 2: GenerateUniformSpherical(anEvent); break;
  case 3: GenerateBeam(anEvent); break;
  case 4: GenerateBeamThrough(anEvent); break;
  case 5: GenerateBeamGausProfile(anEvent); break;
  case 6: GenerateBeamFixSeed(anEvent); break;
  case 7: GenerateScatParticles(anEvent); break;
  case 7001: Generate12XiBeryllium(anEvent); break;
  case 7002: GenerateElementaryXiMinus(anEvent); break;
  case 7003: GenerateElementarySigmaMinus(anEvent); break;
  case 7004: GenerateElementarySigmaPlus(anEvent); break;
    //case 7005: Generate12XiBePeakStructure(anEvent); break;
  case 7501: GenerateKH7XiHSpectrum(anEvent); break;
  case 6301: GenerateE63_7LambdaLi(anEvent, 7); break;
  case 6302: GenerateE63_7LambdaLi(anEvent, 10); break;
  case 6303: GenerateE63_7LambdaLi(anEvent, 12); break;
  case 6374: GenerateE63_LiProfileMonoGamma(anEvent); break;
  case 6375: GenerateE63_4LHGamma(anEvent); break;
  case 9000: GenerateProton(anEvent); break;
  case 9001: GenerateSigmaNCusp(anEvent); break;
  case 9002: GenerateQFLambda(anEvent); break;
  case 9003: GenerateQFSigmaZ(anEvent); break;
  case 9004: GenerateQFSigmaP(anEvent); break;
  case 1001: GenerateK0Production(anEvent); break; // E10
  case 1002: GenerateK0CoherentProduction(anEvent); break; // E10
  default:
    G4cerr << " * Generator number error : " << m_generator << G4endl;
    break;
  }
}

//_____________________________________________________________________________
void // 0
S2SPrimaryGeneratorAction::GenerateDemo(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  // static const G4String name = "proton";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_size = sizeMan.GetSize("Target")*mm/2;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  const auto evnum = anEvent->GetEventID();
  G4double p0 = (experiment == 10) ? 0.9*GeV : 1.4*GeV;
  // p0 = (experiment == 90) ?
  p0 = (evnum%3 == 0) ? p0
    : (evnum%3 == 1) ? p0*1.075
    : p0*0.925;
  auto n = G4RandFlat::shootInt(5);
  G4double theta = (n%5 == 0) ? 0*deg
    : (n%5 == 1) ? 2*deg
    : (n%5 == 2) ? 4*deg
    : (n%5 == 3) ? -2*deg
    : -4*deg;
  G4double phi = 0*deg; // G4RandFlat::shoot(0., 360.)*deg;
  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setRThetaPhi(p0, theta, phi);
  // G4double x0 =  G4RandFlat::shoot(-target_size.x(), target_size.x());
  // G4double y0 =  G4RandFlat::shoot(-target_size.y(), target_size.y());
  // G4double z0 =  G4RandFlat::shoot(-target_size.z(), target_size.z());
  G4LorentzVector v(target_pos, 0);
  // G4LorentzVector v(target_pos + G4ThreeVector(0, 0, 300*CLHEP::mm), 0);
#if 0
  G4cout << FUNC_NAME << G4endl
         << " " << p0 << " " << theta/deg << " " << phi/deg << G4endl
         << " " << p << " " << p.theta()/deg << " " << v << G4endl;
#endif
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(p.v());
  m_particleGun->SetParticleEnergy(p.e() - m0);
  m_particleGun->SetParticlePosition(v.v());
  m_particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
}

//_____________________________________________________________________________
void // 1
S2SPrimaryGeneratorAction::GenerateMonochromeBeam(G4Event* anEvent)
{
  const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  // static const G4String name = "proton";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  const G4double m0 = particle->GetPDGMass();
  G4double p0 = (experiment == 10) ? 0.9*GeV : 1.4*GeV;
  const auto& target_pos = geomMan.GetGlobalPosition("Target");
  beam.pos.setX(target_pos.x());
  beam.pos.setY(target_pos.y());
  beam.pos.setZ(target_pos.z()-956.*mm);
  G4LorentzVector p(0, 0, p0, TMath::Sqrt(p0*p0 + m0*m0));
  G4LorentzVector v(beam.pos, 0);
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(p.v());
  m_particleGun->SetParticleEnergy(p.e() - m0);
  m_particleGun->SetParticlePosition(v.v());
  m_particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
  anaMan.SetPrimaryData(beam.pos.x(),beam.pos.y(),beam.pos.z(),
			0.,0.,0.,0.,p0,0.,9999);
}

//_____________________________________________________________________________
void // 2
S2SPrimaryGeneratorAction::GenerateUniformSpherical(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  // GenPID (optional conf key, see GenPidToName above); 0/unset keeps
  // the historical default of kaon+.
  static const G4int gen_pid = confMan.Get<G4int>("GenPID");
  static const G4String gen_pid_name = GenPidToName(gen_pid);
  static const G4String name = gen_pid_name.empty() ? "kaon+" : gen_pid_name;
  // static const G4String name = "proton";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_size = sizeMan.GetSize("Target")*mm/2;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  // GenMomCent/GenMomBite (GeV/c, optional conf keys): flat momentum
  // window [GenMomCent-GenMomBite/2, GenMomCent+GenMomBite/2]. Unset
  // (either <= 0) falls back to the historical per-experiment ranges.
  static const G4double gen_mom_cent = confMan.Get<G4double>("GenMomCent");
  static const G4double gen_mom_bite = confMan.Get<G4double>("GenMomBite");
  G4double p0;
  if(gen_mom_cent > 0. && gen_mom_bite > 0.){
    p0 = G4RandFlat::shoot(gen_mom_cent - gen_mom_bite/2.,
                           gen_mom_cent + gen_mom_bite/2.)*GeV;
  }else{
    p0 = (experiment == 10)
      ? G4RandFlat::shoot(0.6, 1.2)*GeV
      : G4RandFlat::shoot(1.37, 1.38)*GeV;
    if(experiment==90)
      p0 = G4RandFlat::shoot(0.9, 1.5)*GeV;
  }

  // GenTheta (deg, optional conf key): max polar angle. Unset (<= 0)
  // falls back to 20 deg.
  static const G4double gen_theta_max = confMan.Get<G4double>("GenTheta");
  static const G4double theta_max =
    (gen_theta_max > 0.) ? gen_theta_max*deg : 20.*deg;
  // GenThetaFlat (0/1, optional conf key, default 0): sampling in
  // theta. 0 = solid-angle-uniform (uniform in cos theta, the
  // historical default). 1 = FLAT in theta -- gives EQUAL statistics
  // per theta bin, so a (theta,p) acceptance MAP is not starved at
  // small theta (where dN/dtheta ~ sin theta -> ~0 for cos-uniform,
  // yet that is exactly the high-acceptance forward region). The
  // acceptance is a per-bin Acc/Gen ratio, so it is INVARIANT under
  // this sampling choice -- only the map's statistical noise changes.
  // See analysis-note.md 2026-07-23 acceptance-map entry.
  static const G4int gen_theta_flat = confMan.Get<G4int>("GenThetaFlat");
  // theta_max is already in CLHEP internal angle units (radian==1), so
  // a flat draw over [0, theta_max] is directly the polar angle.
  G4double theta =
    (gen_theta_flat != 0)
      ? G4RandFlat::shoot(0., theta_max)
      : std::acos(G4RandFlat::shoot(std::cos(0.*radian),
                                    std::cos(theta_max)))*radian;
  G4double phi = G4RandFlat::shoot(0., 360.)*deg;
  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setRThetaPhi(p0, theta, phi);
  G4double u0, v0;
  u0 = TMath::Tan(theta)*TMath::Cos(phi);
  v0 = TMath::Tan(theta)*TMath::Sin(phi);
  beam.VO(zK18Target);

  if(experiment!=90){
    const auto vtx = SampleUniformTargetVertex(target_pos, target_size);
    beam.pos.setX(vtx.x());
    beam.pos.setY(vtx.y());
    beam.pos.setZ(vtx.z());
  }
  if(experiment==90){
    double beam_x = G4RandGauss::shoot(target_pos.x(),23.);
    while(1){
      if(fabs(beam_x)<(54./2.))
	break;
      else
	beam_x = G4RandGauss::shoot(target_pos.x(),23.);
    }
    beam.pos.setX(beam_x);
    beam.pos.setY(target_pos.y());
    beam.pos.setZ(target_pos.z());
  }

  G4LorentzVector v(beam.pos, 0);
#if 0
  beam.Print();
  G4cout << FUNC_NAME << G4endl
         << " " << p0 << " " << theta/deg << " " << phi/deg << G4endl
         << " " << p << " " << p.theta()/deg << " " << v << G4endl;
#endif
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(p.v());
  m_particleGun->SetParticleEnergy(p.e() - m0);
  m_particleGun->SetParticlePosition(v.v());
  m_particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
  anaMan.SetPrimaryData(beam.pos.x(),beam.pos.y(),beam.pos.z(),u0,v0,phi,theta,p0,1.8*GeV,9999);
}

//_____________________________________________________________________________
void // 3
S2SPrimaryGeneratorAction::GenerateBeam(G4Event* anEvent)
{
  const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  G4String name;
  const G4double p0 = confMan.Get<G4double>("PK18")*CLHEP::GeV;
  if(experiment == 10){
    if(p0 < 0) name = "pi-";
    if(p0 > 0) name = "pi+";
  }else{
    if(p0 < 0) name = "kaon-";
    if(p0 > 0) name = "kaon+";
  }
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  const G4double m0 = particle->GetPDGMass();
  const auto& target_pos = geomMan.GetGlobalPosition("Target");
  beam.VO(zK18Target);
  beam.pos.setX(target_pos.x());
  beam.pos.setY(target_pos.y());
  beam.pos.setZ(target_pos.z());
  beam.mom.setMag(p0);
  G4LorentzVector p(beam.mom, TMath::Sqrt(p0*p0 + m0*m0));
  G4LorentzVector v(beam.pos, 0);
#if 0
  beam.Print();
  G4cout << FUNC_NAME << G4endl
         << " " << p << " " << v << G4endl;
#endif
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(p.v());
  m_particleGun->SetParticleEnergy(p.e() - m0);
  m_particleGun->SetParticlePosition(v.v());
  m_particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
}

void // 4
S2SPrimaryGeneratorAction::GenerateBeamThrough(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& bac1_pos = geomMan.GetGlobalPosition("BAC1")*mm;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  G4double theta;
  //  std::acos(G4RandFlat::shoot(std::cos(0*deg), std::cos(20*deg)))*radian;
  //
  G4double p0, x0, y0, z0;
  G4double u0, v0;
  G4int count =0;
  while(count<1){
    G4int i = int(G4UniformRand()*(k18track->GetEntries()));
    k18track->GetEntry(i);
    if(/*trigflag[21]>0 &&*/ ntK18==1 && abs(CBtof0[0])<0.2 && chisqrK18[0]<10.){
      z0 = bac1_pos.z() - target_pos.z() -556.*mm;
      // z0 = 556 mm upstrm from BAC1 (approximately dwnstr surface of BH2)
      x0 = xtgt[0] + utgt[0]*z0; // Horizontal direction
      y0 = ytgt[0] + vtgt[0]*z0; // vertical direction
      u0 = utgt[0];
      v0 = vtgt[0];
      p0 = p_3rd[0]*GeV;
      count++;
    }
  }
  beam.VO(zK18Target);
  beam.pos.setX(x0+target_pos.x());
  beam.pos.setY(y0+target_pos.y());
  beam.pos.setZ(z0+target_pos.z());
  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setX(p0*u0/TMath::Sqrt(1+u0*u0+v0*v0));
  p.setY(p0*v0/TMath::Sqrt(1+u0*u0+v0*v0));
  p.setZ(p0/TMath::Sqrt(1+u0*u0+v0*v0));
  G4LorentzVector v(beam.pos, 0);
#if 0
  beam.Print();
  G4cout << FUNC_NAME << G4endl
         << " " << p << " " << v << G4endl;
#endif
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(p.v());
  m_particleGun->SetParticleEnergy(p.e() - m0);
  m_particleGun->SetParticlePosition(v.v());
  m_particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
  anaMan.SetPrimaryData(beam.pos.x(),beam.pos.y(),beam.pos.z(),u0,v0,0.,0.,p0,p0,9999);
}

void // 5
S2SPrimaryGeneratorAction::GenerateBeamGausProfile(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");

  G4double BeamSizeX = 10.*mm;
  G4double BeamSizeY = 5.*mm;
  G4double BeamSizeU = 5.e-3;
  G4double BeamSizeV = 3.e-3;
  G4double MomCenter = 1.4*GeV ;
  G4double MomSize = 0.1*GeV ;

  G4double x0 = 0.0;
  G4double y0 = 0.0;
  G4double z0 = 0.0;
  G4double u0 = 0.0;
  G4double v0 = 0.0;
  G4double p0 = 0.0;

  x0 = G4RandGauss::shoot(target_pos.x(),BeamSizeX);
  y0 = G4RandGauss::shoot(target_pos.y(),BeamSizeY);
  z0 = target_pos.z();

  u0 = G4RandGauss::shoot(0.0,BeamSizeU);
  v0 = G4RandGauss::shoot(0.0,BeamSizeV);

  p0 = G4RandFlat::shoot(MomCenter-MomSize, MomCenter+MomSize);

  beam.VO(zK18Target);
  beam.pos.setX(x0);
  beam.pos.setY(y0);
  beam.pos.setZ(z0);
  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setX(p0*u0/TMath::Sqrt(1+u0*u0+v0*v0));
  p.setY(p0*v0/TMath::Sqrt(1+u0*u0+v0*v0));
  p.setZ(p0/TMath::Sqrt(1+u0*u0+v0*v0));
  G4LorentzVector v(beam.pos, 0);
#if 0
  beam.Print();
  G4cout << FUNC_NAME << G4endl
         << " " << p << " " << v << G4endl;
#endif
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(p.v());
  m_particleGun->SetParticleEnergy(p.e() - m0);
  m_particleGun->SetParticlePosition(v.v());
  m_particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
  anaMan.SetPrimaryData(beam.pos.x(),beam.pos.y(),beam.pos.z(),u0,v0,0.,0.,p0,p0,9999);
}

void // 6
S2SPrimaryGeneratorAction::GenerateBeamFixSeed(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");

  // G4double BeamSizeX = 10.*mm;
  // G4double BeamSizeY = 5.*mm;
  G4double BeamSizeX = 30.*mm;
  G4double BeamSizeY = 30.*mm;

  // G4double BeamSizeU = 5.e-3;
  // G4double BeamSizeV = 3.e-3;
  G4double BeamSizeU = 0.1;
  G4double BeamSizeV = 0.1;
  G4double ThetaSize =3.;
  G4double MomCenter = 1.35*GeV ;
  G4double MomSize = 0.*GeV ;

  G4double x0 = 0.0;
  G4double y0 = 0.0;
  G4double z0 = 0.0;
  G4double u0 = 0.0;
  G4double v0 = 0.0;
  G4double p0 = 0.0;

  ++inum;
  G4Random::setTheSeed(inum);
  //  G4cout << "   Seed = " << G4Random::getTheSeed() << G4endl;
  x0 = G4RandGauss::shoot(target_pos.x(),BeamSizeX);
  y0 = G4RandGauss::shoot(target_pos.y(),BeamSizeY);
  z0 = target_pos.z();

  // G4double theta =
  //   std::acos(G4RandFlat::shoot(std::cos(0*deg), std::cos(ThetaSize*deg)))*radian;
  // G4double phi = G4RandFlat::shoot(0., 360.)*deg;

  G4double theta =G4RandFlat::shoot(0., ThetaSize)*deg;
  G4double phi = G4RandFlat::shoot(0., 360.)*deg;

  //  std::cout<<"x:"<<x0<<", y:"<<y0<<", theta:"<<theta/deg<<", phi:"<<phi/deg<<std::endl;
  // u0 = G4RandGauss::shoot(0.0,BeamSizeU);
  // v0 = G4RandGauss::shoot(0.0,BeamSizeV);

  u0 = TMath::Tan(theta)*TMath::Cos(phi);
  v0 = TMath::Tan(theta)*TMath::Sin(phi);

  p0 = G4RandFlat::shoot(MomCenter-MomSize, MomCenter+MomSize);


  beam.VO(zK18Target);
  beam.pos.setX(x0);
  beam.pos.setY(y0);
  beam.pos.setZ(z0);
  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setX(p0*u0/TMath::Sqrt(1+u0*u0+v0*v0));
  p.setY(p0*v0/TMath::Sqrt(1+u0*u0+v0*v0));
  p.setZ(p0/TMath::Sqrt(1+u0*u0+v0*v0));
  //  G4LorentzVector v(beam.pos+target_pos, 0);
  G4LorentzVector v(beam.pos, 0);
  //  std::cout<<"target:"<<target_pos<<", beam:"<<beam.pos<<", v:"<<v.v()<<std::endl;
#if 0
  beam.Print();
  G4cout << FUNC_NAME << G4endl
         << " " << p << " " << v << G4endl;
#endif
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(p.v());
  m_particleGun->SetParticleEnergy(p.e() - m0);
  m_particleGun->SetParticlePosition(v.v());
  m_particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
  //anaMan.SetPrimaryData(x0,y0,z0,u0,v0,0.,0.,p0,p0,9999);
  anaMan.SetPrimaryData(x0,y0,z0,u0,v0,phi,theta,p0,p0,9999);
}

void // 7
S2SPrimaryGeneratorAction::GenerateScatParticles(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");

  G4double BeamSizeX = 10.*mm;
  G4double BeamSizeY = 5.*mm;
  // G4double BeamSizeU = 5.e-3;
  // G4double BeamSizeV = 3.e-3;
  G4double BeamSizeU = 0.1;
  G4double BeamSizeV = 0.1;
  G4double ThetaSize =15.;
  G4double MomCenter = 1.2*GeV ;
  G4double MomSize = 0.3*GeV ;

  G4double x0 = 0.0;
  G4double y0 = 0.0;
  G4double z0 = 0.0;
  G4double u0 = 0.0;
  G4double v0 = 0.0;
  G4double p0 = 0.0;

  x0 = G4RandGauss::shoot(target_pos.x(),BeamSizeX);
  y0 = G4RandGauss::shoot(target_pos.y(),BeamSizeY);
  z0 = target_pos.z();

  // G4double theta =
  //   std::acos(G4RandFlat::shoot(std::cos(0*deg), std::cos(ThetaSize*deg)))*radian;
  // G4double phi = G4RandFlat::shoot(0., 360.)*deg;

  G4double theta =G4RandFlat::shoot(0., ThetaSize)*deg;
  G4double phi = G4RandFlat::shoot(0., 360.)*deg;

  // u0 = G4RandGauss::shoot(0.0,BeamSizeU);
  // v0 = G4RandGauss::shoot(0.0,BeamSizeV);

  u0 = TMath::Tan(theta)*TMath::Cos(phi);
  v0 = TMath::Tan(theta)*TMath::Sin(phi);

  p0 = G4RandFlat::shoot(MomCenter-MomSize, MomCenter+MomSize);

  beam.VO(zK18Target);
  beam.pos.setX(x0);
  beam.pos.setY(y0);
  beam.pos.setZ(z0);
  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setX(p0*u0/TMath::Sqrt(1+u0*u0+v0*v0));
  p.setY(p0*v0/TMath::Sqrt(1+u0*u0+v0*v0));
  p.setZ(p0/TMath::Sqrt(1+u0*u0+v0*v0));
  //  G4LorentzVector v(beam.pos+target_pos, 0);
  G4LorentzVector v(beam.pos, 0);
  //  std::cout<<"target:"<<target_pos<<", beam:"<<beam.pos<<", v:"<<v.v()<<std::endl;
#if 0
  beam.Print();
  G4cout << FUNC_NAME << G4endl
         << " " << p << " " << v << G4endl;
#endif
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(p.v());
  m_particleGun->SetParticleEnergy(p.e() - m0);
  m_particleGun->SetParticlePosition(v.v());
  m_particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
  //anaMan.SetPrimaryData(x0,y0,z0,u0,v0,0.,0.,p0,p0,9999);
  anaMan.SetPrimaryData(x0,y0,z0,u0,v0,phi,theta,p0,p0,9999);
}

void // 7001 E70 12C(KK)12XiBe kinematics
S2SPrimaryGeneratorAction::Generate12XiBeryllium(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_size = sizeMan.GetSize("Target")*mm/2;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  G4double theta;
  //  std::acos(G4RandFlat::shoot(std::cos(0*deg), std::cos(20*deg)))*radian;
  G4double phi = G4RandFlat::shoot(0., 360.)*deg;
  G4double pB = 1.8*GeV;

  G4double mass_12C = 11.177929*GeV; // 7Li
  G4double mass_12XiBe = (10.2551+1.32171)*GeV; // 6He + Xi
  G4double m_tgt = mass_12C; // mass of target nucleus (7Li)
  G4double m_hyp = mass_12XiBe; // mass of target nucleus (7Li)
  G4double p0;
  G4double cost = 0.92;

  while(cost>1. || cost<0.93){
    p0 = G4RandFlat::shoot(1.33, 1.39)*GeV; // scatter momentum
    G4double Energy_B = TMath::Sqrt(pB*pB + m0*m0);
    G4double Energy_S = TMath::Sqrt(p0*p0 + m0*m0);
    G4double Energy_hyp = Energy_B + m_tgt - Energy_S;

    cost = (m_hyp*m_hyp - Energy_hyp*Energy_hyp + pB*pB + p0*p0)/(2.*pB*p0);
  }

  theta = std::acos(cost)*radian;

  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setRThetaPhi(p0, theta, phi);
  beam.VO(zK18Target);
  beam.pos.setX(target_pos.x());
  beam.pos.setY(target_pos.y());
  beam.pos.setZ(target_pos.z());
  G4LorentzVector v(beam.pos, 0);
#if 0
  beam.Print();
  G4cout << FUNC_NAME << G4endl
         << " " << p0 << " " << theta/deg << " " << phi/deg << G4endl
         << " " << p << " " << p.theta()/deg << " " << v << G4endl;
#endif
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(p.v());
  m_particleGun->SetParticleEnergy(p.e() - m0);
  m_particleGun->SetParticlePosition(v.v());
  m_particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
}

void // 7002 E70 p(KK)Xi kinematics
S2SPrimaryGeneratorAction::GenerateElementaryXiMinus(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_size = sizeMan.GetSize("Target")*mm/2;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  G4double theta;
  //  std::acos(G4RandFlat::shoot(std::cos(0*deg), std::cos(20*deg)))*radian;
  G4double phi = G4RandFlat::shoot(0., 360.)*deg;
  G4double pB = 1.8*GeV;

  G4double m_tgt = particleTable->FindParticle("proton") ->GetPDGMass(); // mass of target nucleus (p)
  G4double m_hyp = particleTable->FindParticle("xi-") ->GetPDGMass(); // mass of target nucleus (Xi)
  G4double p0;
  G4double cost = 0.92;

  while(cost>1. || cost<0.93){
    p0 = G4RandFlat::shoot(1.1, 1.32)*GeV; // scatter momentum
    G4double Energy_B = TMath::Sqrt(pB*pB + m0*m0);
    G4double Energy_S = TMath::Sqrt(p0*p0 + m0*m0);
    G4double Energy_hyp = Energy_B + m_tgt - Energy_S;

    cost = (m_hyp*m_hyp - Energy_hyp*Energy_hyp + pB*pB + p0*p0)/(2.*pB*p0);
  }

  theta = std::acos(cost)*radian;

  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setRThetaPhi(p0, theta, phi);
  beam.VO(zK18Target);
  beam.pos.setX(target_pos.x());
  beam.pos.setY(target_pos.y());
  beam.pos.setZ(target_pos.z());
  G4LorentzVector v(beam.pos, 0);
#if 0
  beam.Print();
  G4cout << FUNC_NAME << G4endl
         << " " << p0 << " " << theta/deg << " " << phi/deg << G4endl
         << " " << p << " " << p.theta()/deg << " " << v << G4endl;
#endif
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(p.v());
  m_particleGun->SetParticleEnergy(p.e() - m0);
  m_particleGun->SetParticlePosition(v.v());
  m_particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
}

void // 7003 E70 p(Kpi)SigmaMinus kinematics
S2SPrimaryGeneratorAction::GenerateElementarySigmaMinus(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "pi+";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_size = sizeMan.GetSize("Target")*mm/2;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  G4double theta;
  //  std::acos(G4RandFlat::shoot(std::cos(0*deg), std::cos(20*deg)))*radian;
  G4double phi = G4RandFlat::shoot(0., 360.)*deg;
  static const G4double pB = 1.8*GeV;
  static const G4double mB = particleTable->FindParticle("kaon-") ->GetPDGMass(); // mass of target nucleus (p)

  static const G4double m_tgt = particleTable->FindParticle("proton") ->GetPDGMass(); // mass of target nucleus (p)
  static const G4double m_hyp = particleTable->FindParticle("sigma-") ->GetPDGMass(); // mass of target nucleus (Sigma)
  G4double p0;
  G4double cost = 0.92;

  while(cost>1. || cost<0.93){
    p0 = G4RandFlat::shoot(1.4, 1.6)*GeV; // scatter momentum
    G4double Energy_B = TMath::Sqrt(pB*pB + mB*mB);
    G4double Energy_S = TMath::Sqrt(p0*p0 + m0*m0);
    G4double Energy_hyp = Energy_B + m_tgt - Energy_S;

    cost = (m_hyp*m_hyp - Energy_hyp*Energy_hyp + pB*pB + p0*p0)/(2.*pB*p0);
  }

  theta = std::acos(cost)*radian;

  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setRThetaPhi(p0, theta, phi);
  beam.VO(zK18Target);
  beam.pos.setX(target_pos.x());
  beam.pos.setY(target_pos.y());
  beam.pos.setZ(target_pos.z());
  G4LorentzVector v(beam.pos, 0);
#if 0
  beam.Print();
  G4cout << FUNC_NAME << G4endl
         << " " << p0 << " " << theta/deg << " " << phi/deg << G4endl
         << " " << p << " " << p.theta()/deg << " " << v << G4endl;
#endif
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(p.v());
  m_particleGun->SetParticleEnergy(p.e() - m0);
  m_particleGun->SetParticlePosition(v.v());
  m_particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
}

void // 7004 E70 p(pi,K)SigmaPlus kinematics
S2SPrimaryGeneratorAction::GenerateElementarySigmaPlus(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_size = sizeMan.GetSize("Target")*mm/2;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  G4double theta;
  //  std::acos(G4RandFlat::shoot(std::cos(0*deg), std::cos(20*deg)))*radian;
  G4double phi = G4RandFlat::shoot(0., 360.)*deg;
  G4double pB = 1.8*GeV;
  G4double mB = particleTable->FindParticle("pi+") ->GetPDGMass(); // mass of target nucleus (p)

  G4double m_tgt = particleTable->FindParticle("proton") ->GetPDGMass(); // mass of target nucleus (p)
  G4double m_hyp = particleTable->FindParticle("sigma+") ->GetPDGMass(); // mass of target nucleus (Sigma)
  G4double p0;
  G4double cost = 0.92;

  while(cost>1. || cost<0.93){
    p0 = G4RandFlat::shoot(1.2, 1.42)*GeV; // scatter momentum
    G4double Energy_B = TMath::Sqrt(pB*pB + mB*mB);
    G4double Energy_S = TMath::Sqrt(p0*p0 + m0*m0);
    G4double Energy_hyp = Energy_B + m_tgt - Energy_S;

    cost = (m_hyp*m_hyp - Energy_hyp*Energy_hyp + pB*pB + p0*p0)/(2.*pB*p0);
  }

  theta = std::acos(cost)*radian;

  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setRThetaPhi(p0, theta, phi);
  beam.VO(zK18Target);
  beam.pos.setX(target_pos.x());
  beam.pos.setY(target_pos.y());
  beam.pos.setZ(target_pos.z());
  G4LorentzVector v(beam.pos, 0);
#if 0
  beam.Print();
  G4cout << FUNC_NAME << G4endl
         << " " << p0 << " " << theta/deg << " " << phi/deg << G4endl
         << " " << p << " " << p.theta()/deg << " " << v << G4endl;
#endif
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(p.v());
  m_particleGun->SetParticleEnergy(p.e() - m0);
  m_particleGun->SetParticlePosition(v.v());
  m_particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
}


void // 7501 E75 phase-1 7XiH+6XiH spectrum kinematics
S2SPrimaryGeneratorAction::GenerateKH7XiHSpectrum(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const G4String name = "kaon+";
  static const auto particle = particleTable->FindParticle(name);
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_size = sizeMan.GetSize("Target")*mm/2;
  static const G4double m0 = particle->GetPDGMass();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  G4double theta;
  //  std::acos(G4RandFlat::shoot(std::cos(0*deg), std::cos(20*deg)))*radian;
  //
  G4double pB, x0, y0, z0;
  G4double u0, v0;
  G4int count =0;
  while(count<1){
    G4int i = int(G4UniformRand()*(k18track->GetEntries()));
    k18track->GetEntry(i);
    if(trigflag[21]>0 && ntK18>0 && abs(CBtof0[0])<0.2 && chisqrK18[0]<10.){
      //z0 = G4RandFlat::shoot(-target_size.z()/2, target_size.z()/2);
      z0 = G4RandFlat::shoot(-187.*mm/2, 187*mm/2);
      x0 = xtgt[0] + utgt[0]*z0; // Horizontal direction
      y0 = ytgt[0] + vtgt[0]*z0; // vertical direction
      pB = p_3rd[0]*GeV;
      count++;
    }
  }
  G4double phi = G4RandFlat::shoot(0., 360.)*deg;
  G4double p0=0.;
  G4double BE=0.;

  G4double r = Sum*G4UniformRand();
  G4double yint =0;
  for(Int_t i=0;i<gr->GetN();++i){
    Double_t x, y;
    gr->GetPoint(i,x,y);
    yint += y;
    if(yint>r){
      BE=x;
      break;
    }
  }

  G4double mass_7Li = 6.53468*GeV; // 7Li
  G4double mass_7XiH = (5.61153+1.32171+BE*1e-3)*GeV; // 6He + Xi
  G4double mT = mass_7Li; // mass of target nucleus (7Li)
  G4double mY = mass_7XiH; // mass of target nucleus (7Li)

  G4double cost = 0.92;
  G4int count_pt=0;

  while(cost>1. || cost<0.93){
    p0 = G4RandFlat::shoot(1.27, 1.47)*GeV; // scatter momentum
    G4double EB = TMath::Sqrt(pB*pB + m0*m0);
    G4double ES = TMath::Sqrt(p0*p0 + m0*m0);
    G4double EY = EB + mT - ES;
    cost = (mY*mY - EY*EY + pB*pB + p0*p0)/(2.*pB*p0);
    count_pt ++;
    if(count_pt >10){
      cost =0.5;
      p0 = 0.0*GeV;
      break;
    }
  }
  std::cout <<"BE: "<< BE <<", pB: "<< pB<<", pS: "<< p0 <<", Theta: "<<cost<< std::endl;
  theta = std::acos(cost)*radian;

  G4LorentzVector p(0, 0, 0, TMath::Sqrt(p0*p0 + m0*m0));
  p.setRThetaPhi(p0, theta, phi);
  beam.VO(zK18Target);

  beam.pos.setX(x0+target_pos.x());
  beam.pos.setY(y0+target_pos.y());
  beam.pos.setZ(z0+target_pos.z());
  G4LorentzVector v(beam.pos, 0);
  u0 = TMath::Tan(theta)*TMath::Cos(phi);
  v0 = TMath::Tan(theta)*TMath::Sin(phi);
#if 0
  beam.Print();
  G4cout << FUNC_NAME << G4endl
         << " " << p0 << " " << theta/deg << " " << phi/deg << G4endl
         << " " << p << " " << p.theta()/deg << " " << v << G4endl;
#endif
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(p.v());
  m_particleGun->SetParticleEnergy(p.e() - m0);
  m_particleGun->SetParticlePosition(v.v());
  m_particleGun->GeneratePrimaryVertex(anEvent);
  anaMan.SetPrimaryParticle(0, pdg, p, v);
  anaMan.SetPrimaryData(beam.pos.x(),beam.pos.y(),beam.pos.z(),u0,v0,phi,theta,p0,pB,9999);
}

//_____________________________________________________________________________
// E63
void // 6301~  [ E63 A(K-,pi-)lambda_hyper kinematics ]
S2SPrimaryGeneratorAction::GenerateE63_7LambdaLi(G4Event* anEvent, G4int MassNum)
{
  static const G4int n_particle = 1;  // should 1 even if you generate weak pion
  m_particleGun = new G4ParticleGun(n_particle);
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_size = sizeMan.GetSize("Target")*mm/2;
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  auto WeakParticle = confMan.Get<G4String>("WeakDecayParticle");
  int dummy_event_flag = 0;
  int dummy_event_flag_out_of_tgt = 0; // out of target
  int dummy_event_flag_out_of_cs = 0; // out of cross section

  // ***************
  // *** beam K- ***
  static const auto beam_particle = particleTable->FindParticle("kaon-");
  static const auto beam_pdg =      beam_particle->GetPDGEncoding();
  static const G4double m_beam =    beam_particle->GetPDGMass();
  // momentum
  G4double p_beam = confMan.Get<G4double>("PK18")*CLHEP::GeV;
  if(1){ // include BeamMomentuBite
    p_beam += G4RandGauss::shoot( 0.0,  p_beam*0.0134 );  // momentum bite [202501 data]
  }
  // momentum direction
  G4double beam_u_rms = 0.0164; // [202501 data]
  G4double beam_v_rms = 0.0045;
  G4double beam_u = G4RandGauss::shoot( 0.0, beam_u_rms );
  G4double beam_v = G4RandGauss::shoot( 0.0, beam_v_rms );
  G4double beam_z = 1./sqrt(1. + beam_u*beam_u + beam_v*beam_v);
  G4ThreeVector BeamMomDir = G4ThreeVector( beam_u*beam_z, beam_v*beam_z, beam_z );
  // profile at target center -> vertex at terget coodinate
  G4double beam_x_rms = 16.94; // [202501 data]
  G4double beam_y_rms =  7.01;
  G4double vertex_x = G4RandGauss::shoot( 0.0, beam_x_rms );
  G4double vertex_y = G4RandGauss::shoot( 0.0, beam_y_rms );
  G4double vertex_z = G4RandFlat::shoot( -target_size.z(), target_size.z() );
  G4ThreeVector VertexPos = G4ThreeVector( vertex_x, vertex_y, vertex_z );
  G4LorentzVector VertexLv(VertexPos, 0);

  if(1){  // InsideTarget cut
    if( fabs(vertex_x)>fabs(target_size.x()) || fabs(vertex_y)>fabs(target_size.y()) ){
      //G4cerr << "out of target" << G4endl;
      dummy_event_flag_out_of_tgt = 1;  // set dummy vertex if vertex is outof target_size
    }
  }

  // ****************
  // *** scat pi- ***
  static const auto scat_particle = particleTable->FindParticle("pi-");
  static const auto scat_pdg =      scat_particle->GetPDGEncoding();
  static const G4double m_scat =    scat_particle->GetPDGMass();
  G4double costLab;
  // limit of scat angle
  {
    //costLab = 1. - G4RandFlat::shoot(0., 0.003805); // (0- 5 deg.)
    costLab = 1. - G4RandFlat::shoot(0., 0.015192); // (0-10 deg.)
    //costLab = 1. - G4RandFlat::shoot(0., 0.034074); // (0-15 deg.)
    //costLab = 1. - G4RandFlat::shoot(0., 0.060307); // (0-20 deg.)
    //costLab = 1. - G4RandFlat::shoot(0., 0.093692); // (0-25 deg.)
    //costLab = 1. - G4RandFlat::shoot(0., 0.133974); // (0-30 deg.)
  }

  if(1){ // cross section shape cut
    G4int DeltaL = 0;  // need 0 or 1 or 2
    G4double p[6]; // f_cross = pol(6)
    G4double MaxCrossSection = 1100; // [a.u.]
    if(DeltaL==0){  // table for 0.9 GeV/c Li (k.pi)
      p[0]=1033; p[1]=70.55; p[2]=-59.28;
      p[3]=7.2325; p[4]=-0.37978; p[5]=0.00927192;
      p[6]=-8.56922e-5;
    }
    else if(DeltaL==1){
      p[0]=259.9; p[1]=-38.1581; p[2]=27.6301;
      p[3]=-3.38422; p[4]=0.149127; p[5]=-0.00233459;
      p[6]=3.37185e-6;
    }
    else if(DeltaL==2){
      p[0]=140.35; p[1]=31.2678; p[2]=-17.7024;
      p[3]=4.27016; p[4]=-0.3732; p[5]=0.0135775;
      p[6]=-0.000177677;
    }
    else{
      G4cerr << "DeltaL setting is wrong" << G4endl;
      exit(-1);
    }

    G4double thetaLab = acos(costLab)*(180./3.141592); // degree
    G4double RandValue = G4RandFlat::shoot(0., MaxCrossSection);
    G4double cross_section = p[0];
    for(int n=1; n<7; n++){
      cross_section += p[n]*pow(thetaLab,n);
    }
    if( RandValue>cross_section ){
      dummy_event_flag_out_of_cs = 1; // set dummy vertex if RandValue > cross section table
    }
  } // cross section cut

  // ****************
  // *** weak pi- ***
  static const auto weak_particle = particleTable->FindParticle("pi-");
  static const auto weak_pdg =      weak_particle->GetPDGEncoding();
  static const G4double m_weak =    weak_particle->GetPDGMass();
  G4double WeakT = 0.;
  if(WeakParticle == "3LH") WeakT = 40.88 *CLHEP::MeV ; // for 3LH
  if(WeakParticle == "4LH") WeakT = 53.25 *CLHEP::MeV ; // for 4LH
  if(WeakParticle == "6LH") WeakT = 37.20 *CLHEP::MeV ; // for 6LHe
  if(WeakT <= 0.){
    G4cerr << "Unknown WeakDecayParticle: " << WeakParticle << G4endl;
    exit(-1);
  }
  G4double WeakMom = sqrt( pow(WeakT+m_weak,2) -m_scat*m_scat );


  // *************************************
  // *** target nuclei and hypernuclei ***
  G4double AtomicMassUnit = 0.93149432;
  G4double LambdaMass = particleTable->FindParticle("lambda")->GetPDGMass();
  G4double m_tgt, m_hyp;
  {
    // 7Li -----------
    G4double mass_7Li = (7.0*AtomicMassUnit+0.014908)*GeV;
    //G4double mass_7LambdaLi = (6.0*AtomicMassUnit+0.014086-0.00522+0.000)*GeV + LambdaMass; // Ex=0 MeV
    G4double mass_7LambdaLi = (6.0*AtomicMassUnit+0.014086-0.00522+0.020)*GeV + LambdaMass; // Ex=20 MeV

    // 10B -----------
    G4double mass_10B = (10.0*AtomicMassUnit+0.0120508)*GeV;
    G4double mass_10LambdaB = (9.0*AtomicMassUnit+0.0113477-0.0081+0.000)*GeV + LambdaMass; // Ex=0 MeV

    // 12C -----------
    G4double mass_12C = (12.0*AtomicMassUnit+0.0)*GeV;
    G4double mass_12LambdaC = (11.0*AtomicMassUnit+0.010650-0.0108+0.000)*GeV + LambdaMass; // Ex=0 MeV

    if(MassNum==7){
      m_tgt = mass_7Li;
      m_hyp = mass_7LambdaLi;
    }
    else if(MassNum==10){
      m_tgt = mass_10B;
      m_hyp = mass_10LambdaB;
    }
    else if(MassNum==12){
      m_tgt = mass_12C;
      m_hyp = mass_12LambdaC;
    }
    else{
      G4cerr << "not regstered MassNumber  MassNum = " << MassNum << G4endl;
      exit(-1);
    }
  }


  // ************************
  // **** calculate scat ****
  // ************************

  //Kaon 1.5GeV/c
  G4LorentzVector BeamLv( p_beam*BeamMomDir,
			  sqrt( m_beam*m_beam+p_beam*p_beam ) );
  //Neutron 0.0GeV/c
  G4double NuclMom = 0.0;
  G4ThreeVector NuclMomDir( 0., 0., NuclMom );
  G4LorentzVector NuclLv( NuclMom*NuclMomDir,
			  sqrt( m_tgt*m_tgt+NuclMom*NuclMom )  );

  //Primary frame
  G4LorentzVector PrimaryLv =  BeamLv+NuclLv;
  G4double TotalEnergyCM = PrimaryLv.mag();
  G4ThreeVector beta( PrimaryLv.vect()/PrimaryLv.e() );

  //scat CM
  G4double ScatMomCM
    = 0.5*sqrt(( TotalEnergyCM*TotalEnergyCM
  		 -( m_scat+m_hyp )*( m_scat+m_hyp ))
  	       *( TotalEnergyCM*TotalEnergyCM
  		  -( m_scat-m_hyp )*( m_scat-m_hyp )))/TotalEnergyCM;

  G4double cottLab=costLab/sqrt(1.-costLab*costLab);
  G4double bt=beta.mag(), gamma=1./sqrt(1.-bt*bt);
  G4double gbep=gamma*bt*sqrt(ScatMomCM*ScatMomCM+m_scat*m_scat)/ScatMomCM;
  G4double a  = gamma*gamma+cottLab*cottLab;
  G4double bp = gamma*gbep;
  G4double c  = gbep*gbep-cottLab*cottLab;

  G4double dd=bp*bp-a*c;
  if( dd<0. ){
    G4cerr << "dd<0." << G4endl;
    exit(-1);
  }

  G4double costCM=(sqrt(dd)-bp)/a;
  if( costCM>1. || costCM<-1. ){
    G4cerr << "costCM>1. || costCM<-1." << G4endl;
    exit(-1);
  }

  G4double sintCM=sqrt(1.-costCM*costCM);
  G4double phiCM=G4RandFlat::shoot(0., 360.)*deg;
  G4ThreeVector ScatMomCM_vector( ScatMomCM*sintCM*cos(phiCM),
				  ScatMomCM*sintCM*sin(phiCM),
				  ScatMomCM*costCM );
  //ScatMomCM_vector.rotateY(KaonMomDir.theta());
  //ScatMomCM_vector.rotateZ(KaonMomDir.phi());
  ScatMomCM_vector.rotateUz(BeamMomDir);

  G4LorentzVector ScatLv( ScatMomCM_vector,
			  sqrt( ScatMomCM*ScatMomCM + m_scat*m_scat ));
  ScatLv.boost(beta);

  G4ThreeVector ScatMom_vector = ScatLv.vect();
  G4double ScatMom = ScatMom_vector.mag();
  G4ThreeVector ScatMomDir = ScatMom_vector/ScatMom;
  G4double ScatT = sqrt( ScatMom*ScatMom + m_scat*m_scat ) - m_scat;


  // ************************
  // **** calculate weak ****
  // ************************
  G4ThreeVector WeakMomDir = G4RandomDirection();
  G4ThreeVector WeakMom_vector = WeakMom * WeakMomDir;
  // NOTE: unused. The 4th argument must be the TOTAL energy (WeakT + m_weak),
  // not the kinetic energy WeakT — fix before using this Lorentz vector.
  G4LorentzVector WeakLv( WeakMom_vector, WeakT);


  // ***************************
  // **** particle generate ****
  // ***************************

  // ****************************
  // ** scat particle generate **
  if(dummy_event_flag_out_of_tgt || dummy_event_flag_out_of_cs) dummy_event_flag = 1;
  if(dummy_event_flag){ // set primary vertex at out of world
    ScatMomDir = G4ThreeVector(0, 0, -1);
    ScatMom = 0; ScatT = 0;
    VertexPos = G4ThreeVector(-200*CLHEP::m, -200*CLHEP::m, -200*CLHEP::m);
  }
  m_particleGun->SetParticleDefinition(scat_particle);
  m_particleGun->SetParticleMomentumDirection(ScatMomDir);
  m_particleGun->SetParticleEnergy(ScatT);
  m_particleGun->SetParticlePosition(VertexPos + target_pos);
  m_particleGun->GeneratePrimaryVertex(anEvent);
  // analyzer fill : scat particle
  {
    G4double x0 = VertexPos.x(); G4double y0 = VertexPos.y(); G4double z0 = VertexPos.z();
    G4double u0 = ScatMomDir.x()/ScatMomDir.z(); G4double v0 = ScatMomDir.y()/ScatMomDir.z();
    G4double phi = ScatMomDir.phi(); G4double theta = ScatMomDir.theta();
    G4double p0 = ScatMom; G4double pB = ScatT;
    anaMan.SetPrimaryData(x0,y0,z0,u0,v0,phi,theta,p0,pB,9999);
    // unit: x0,y0,z0=mm, phi,theta=rad, p0[momentum]=MeV/c, pB[KineticEnergy]=MeV
    anaMan.SetPrimaryParticle(0, scat_pdg, ScatLv, VertexLv);
  }

  // **********************************
  // ** weak decay particle generate **
  if(1){
    if(dummy_event_flag){ // set primary vertex at out of world
      WeakMomDir = G4ThreeVector(0, 0, -1);
      WeakMom = 0; WeakT = 0;
      if(dummy_event_flag_out_of_tgt){
	VertexPos = G4ThreeVector(-200*CLHEP::m, -200*CLHEP::m, -200*CLHEP::m);
	if(dummy_event_flag_out_of_cs)
	VertexPos = G4ThreeVector(-190*CLHEP::m, -190*CLHEP::m, -190*CLHEP::m);
      }
      if(dummy_event_flag_out_of_cs)
	VertexPos = G4ThreeVector(-180*CLHEP::m, -180*CLHEP::m, -180*CLHEP::m);
    }

    m_particleGun->SetParticleDefinition(weak_particle);
    m_particleGun->SetParticleMomentumDirection(WeakMomDir);
    m_particleGun->SetParticleEnergy(WeakT);
    m_particleGun->SetParticlePosition(VertexPos + target_pos);
    m_particleGun->GeneratePrimaryVertex(anEvent);
    // analyzer fill : weak decay particle
    {
      G4double x1 = VertexPos.x(); G4double y1 = VertexPos.y(); G4double z1 = VertexPos.z();
      G4double px1 = WeakMom_vector.x(); G4double py1 = WeakMom_vector.y(); G4double pz1 = WeakMom_vector.z();
      G4double p1 = WeakMom; G4double t1 = WeakT;
      anaMan.SetSecondaryData(x1,y1,z1,px1,py1,pz1,p1,t1);
      // unit: x1,y1,z1=mm, px1,py1,pz1,p1[momentum]=MeV/c, t1[KineticEnergy]=MeV
    }
  } // if weak decay particle

}

//_____________________________________________________________________________
// 6375: 4LambdaH gamma template with the Generator 6301 production model.
// Includes the K- beam bite/profile, target-inside rejection, 0-10 deg pi-
// production acceptance, and DeltaL=0 cross-section rejection before the
// 7LambdaLi* -> 4LambdaH* + 3He -> gamma decay chain.  The 7LambdaLi*
// excitation energy is sampled event-by-event for the 4LH window study.
void
S2SPrimaryGeneratorAction::GenerateE63_4LHGamma(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_size = sizeMan.GetSize("Target")*mm/2;
  int dummy_event_flag = 0;
  int dummy_event_flag_out_of_tgt = 0;
  int dummy_event_flag_out_of_cs = 0;

  static const auto beam_particle = particleTable->FindParticle("kaon-");
  static const G4double m_beam = beam_particle->GetPDGMass();
  G4double p_beam = confMan.Get<G4double>("PK18")*CLHEP::GeV;
  p_beam += G4RandGauss::shoot(0.0, p_beam*0.0134);

  const G4double beam_u = G4RandGauss::shoot(0.0, 0.0164);
  const G4double beam_v = G4RandGauss::shoot(0.0, 0.0045);
  const G4double beam_z = 1./std::sqrt(1. + beam_u*beam_u + beam_v*beam_v);
  const G4ThreeVector BeamMomDir(beam_u*beam_z, beam_v*beam_z, beam_z);

  G4double vertex_x = G4RandGauss::shoot(0.0, 16.94);
  G4double vertex_y = G4RandGauss::shoot(0.0, 7.01);
  G4double vertex_z = G4RandFlat::shoot(-target_size.z(), target_size.z());
  G4ThreeVector VertexPos(vertex_x, vertex_y, vertex_z);
  G4LorentzVector VertexLv(VertexPos, 0);

  if(std::abs(vertex_x) > std::abs(target_size.x())
     || std::abs(vertex_y) > std::abs(target_size.y()))
    dummy_event_flag_out_of_tgt = 1;

  static const auto scat_particle = particleTable->FindParticle("pi-");
  static const auto scat_pdg = scat_particle->GetPDGEncoding();
  static const G4double m_scat = scat_particle->GetPDGMass();
  static const auto gamma_particle = particleTable->FindParticle("gamma");
  static const auto gamma_pdg = gamma_particle->GetPDGEncoding();

  const G4double costLab = 1. - G4RandFlat::shoot(0., 0.015192);
  {
    const std::array<G4double, 7> p = {{
      1033., 70.55, -59.28, 7.2325, -0.37978, 0.00927192, -8.56922e-5
    }};
    const G4double thetaLab = std::acos(costLab)*(180./TMath::Pi());
    G4double cross_section = p[0];
    for(std::size_t n=1; n<p.size(); ++n)
      cross_section += p[n]*std::pow(thetaLab, static_cast<G4double>(n));
    if(G4RandFlat::shoot(0., 1100.) > cross_section)
      dummy_event_flag_out_of_cs = 1;
  }

  const G4double AtomicMassUnit = 0.93149432;
  const G4double LambdaMass = particleTable->FindParticle("lambda")->GetPDGMass();
  const G4double hyper_ex = G4RandFlat::shoot(20.0, 39.0)*CLHEP::MeV;
  const G4double gamma_rest_input = confMan.Get<G4double>("GammaRestE");
  const G4double gamma_rest_e = ((gamma_rest_input > 0.) ? gamma_rest_input : 1.090)
    * CLHEP::MeV;
  const G4double binding_4lh_input = confMan.Get<G4double>("FourLambdaHBinding");
  const G4double binding_4lh = ((binding_4lh_input > 0.) ? binding_4lh_input : 2.157)
    * CLHEP::MeV;

  const G4double mass_7Li = (7.0*AtomicMassUnit + 0.014908)*GeV;
  const G4double mass_7LambdaLi =
    (6.0*AtomicMassUnit + 0.014086 - 0.00522)*GeV + LambdaMass + hyper_ex;
  const G4double mass_3H = (3.0*AtomicMassUnit + 0.014949806)*GeV;
  const G4double mass_3He = (3.0*AtomicMassUnit + 0.014931214)*GeV;
  const G4double mass_4LambdaH_gs = mass_3H + LambdaMass - binding_4lh;
  const G4double mass_4LambdaH_star =
    gamma_rest_e + std::sqrt(mass_4LambdaH_gs*mass_4LambdaH_gs
                             + gamma_rest_e*gamma_rest_e);
  static constexpr G4int pdg_7LambdaLi = 1010030070;
  static constexpr G4int pdg_4LambdaH = 1010010040;
  static constexpr G4int pdg_3He = 1000020030;

  if(mass_7LambdaLi <= mass_4LambdaH_star + mass_3He){
    G4cerr << "GenerateE63_4LHGamma below threshold: HyperNucleusEx="
           << hyper_ex/CLHEP::MeV << " MeV, GammaRestE="
           << gamma_rest_e/CLHEP::MeV << " MeV, FourLambdaHBinding="
           << binding_4lh/CLHEP::MeV << " MeV" << G4endl;
    exit(-1);
  }

  const G4LorentzVector BeamLv(
    p_beam*BeamMomDir, std::sqrt(m_beam*m_beam + p_beam*p_beam));
  const G4LorentzVector NuclLv(
    G4ThreeVector(0., 0., 0.), mass_7Li);
  const G4LorentzVector PrimaryLv = BeamLv + NuclLv;
  const G4double TotalEnergyCM = PrimaryLv.mag();
  const G4ThreeVector beta(PrimaryLv.vect()/PrimaryLv.e());

  const G4double ScatMomCM =
    0.5*std::sqrt((TotalEnergyCM*TotalEnergyCM
                   - (m_scat + mass_7LambdaLi)*(m_scat + mass_7LambdaLi))
                  *(TotalEnergyCM*TotalEnergyCM
                    - (m_scat - mass_7LambdaLi)*(m_scat - mass_7LambdaLi)))
    / TotalEnergyCM;
  const G4double cottLab = costLab/std::sqrt(1. - costLab*costLab);
  const G4double bt = beta.mag();
  const G4double gamma = 1./std::sqrt(1. - bt*bt);
  const G4double gbep =
    gamma*bt*std::sqrt(ScatMomCM*ScatMomCM + m_scat*m_scat)/ScatMomCM;
  const G4double a = gamma*gamma + cottLab*cottLab;
  const G4double bp = gamma*gbep;
  const G4double c = gbep*gbep - cottLab*cottLab;
  const G4double dd = bp*bp - a*c;
  if(dd < 0.){
    G4cerr << "GenerateE63_4LHGamma dd<0." << G4endl;
    exit(-1);
  }
  const G4double costCM = (std::sqrt(dd) - bp)/a;
  if(costCM > 1. || costCM < -1.){
    G4cerr << "GenerateE63_4LHGamma costCM outside [-1,1]" << G4endl;
    exit(-1);
  }

  const G4double sintCM = std::sqrt(1. - costCM*costCM);
  const G4double phiCM = G4RandFlat::shoot(0., 360.)*deg;
  G4ThreeVector ScatMomCM_vector(ScatMomCM*sintCM*std::cos(phiCM),
                                 ScatMomCM*sintCM*std::sin(phiCM),
                                 ScatMomCM*costCM);
  ScatMomCM_vector.rotateUz(BeamMomDir);
  G4LorentzVector ScatLv(
    ScatMomCM_vector, std::sqrt(ScatMomCM*ScatMomCM + m_scat*m_scat));
  ScatLv.boost(beta);

  G4ThreeVector ScatMom_vector = ScatLv.vect();
  G4double ScatMom = ScatMom_vector.mag();
  G4ThreeVector ScatMomDir = ScatMom_vector/ScatMom;
  G4double ScatT = std::sqrt(ScatMom*ScatMom + m_scat*m_scat) - m_scat;
  const G4LorentzVector HyperLv = BeamLv + NuclLv - ScatLv;

  if(dummy_event_flag_out_of_tgt || dummy_event_flag_out_of_cs)
    dummy_event_flag = 1;
  if(dummy_event_flag){
    ScatMomDir = G4ThreeVector(0., 0., -1.);
    ScatMom = 0.;
    ScatT = 0.;
    VertexPos = G4ThreeVector(-200*CLHEP::m, -200*CLHEP::m, -200*CLHEP::m);
  }

  m_particleGun->SetParticleDefinition(scat_particle);
  m_particleGun->SetParticleMomentumDirection(ScatMomDir);
  m_particleGun->SetParticleEnergy(ScatT);
  m_particleGun->SetParticlePosition(VertexPos + target_pos);
  m_particleGun->GeneratePrimaryVertex(anEvent);

  {
    const G4double u0 = ScatMomDir.x()/ScatMomDir.z();
    const G4double v0 = ScatMomDir.y()/ScatMomDir.z();
    anaMan.SetPrimaryData(VertexPos.x(), VertexPos.y(), VertexPos.z(),
                          u0, v0, ScatMomDir.phi(), ScatMomDir.theta(),
                          ScatMom, ScatT, 9999);
    anaMan.SetPrimaryParticle(0, scat_pdg, ScatLv, VertexLv);
  }

  G4ThreeVector GammaMomDir(0., 0., -1.);
  G4double GammaE = 0.;
  if(!dummy_event_flag){
    RecordGeneratedParticle(GenBranch::kPrimPi,
                            GenBranch::ParticleId::None,
                            scat_pdg,
                            ScatLv,
                            VertexLv);
    RecordGeneratedParticle(GenBranch::kHypNucleus,
                            GenBranch::ParticleId::None,
                            pdg_7LambdaLi,
                            HyperLv,
                            VertexLv);

    const G4double breakup_mom =
      TwoBodyMomentum(mass_7LambdaLi, mass_4LambdaH_star, mass_3He);
    const G4ThreeVector frag_dir_hyp = G4RandomDirection();
    G4LorentzVector FragStarLvHyp(
      frag_dir_hyp*breakup_mom,
      std::sqrt(breakup_mom*breakup_mom
                + mass_4LambdaH_star*mass_4LambdaH_star));
    G4LorentzVector RecoilLvHyp(
      -frag_dir_hyp*breakup_mom,
      std::sqrt(breakup_mom*breakup_mom + mass_3He*mass_3He));
    G4LorentzVector FragStarLvLab = FragStarLvHyp;
    G4LorentzVector RecoilLvLab = RecoilLvHyp;
    const G4ThreeVector beta_hyp = HyperLv.vect()/HyperLv.e();
    FragStarLvLab.boost(beta_hyp);
    RecoilLvLab.boost(beta_hyp);

    const G4double gamma_mom =
      TwoBodyMomentum(mass_4LambdaH_star, mass_4LambdaH_gs, 0.);
    const G4ThreeVector gamma_dir_frag = G4RandomDirection();
    G4LorentzVector GammaLvFrag(gamma_mom*gamma_dir_frag, gamma_mom);
    G4LorentzVector GammaLvLab = GammaLvFrag;
    GammaLvLab.boost(FragStarLvLab.vect()/FragStarLvLab.e());

    GammaMomDir = GammaLvLab.vect().unit();
    GammaE = GammaLvLab.e();

    RecordGeneratedParticle(GenBranch::kHypFragment,
                            GenBranch::ParticleId::HypNucleus,
                            pdg_4LambdaH,
                            FragStarLvLab,
                            VertexLv);
    RecordGeneratedParticle(GenBranch::kRecoilIon,
                            GenBranch::ParticleId::HypNucleus,
                            pdg_3He,
                            RecoilLvLab,
                            VertexLv);
    RecordGeneratedParticle(GenBranch::kDecGamma,
                            GenBranch::ParticleId::HypFragment,
                            gamma_pdg,
                            GammaLvLab,
                            VertexLv);
  } else {
    if(dummy_event_flag_out_of_tgt && dummy_event_flag_out_of_cs)
      VertexPos = G4ThreeVector(-190*CLHEP::m, -190*CLHEP::m, -190*CLHEP::m);
    else if(dummy_event_flag_out_of_cs)
      VertexPos = G4ThreeVector(-180*CLHEP::m, -180*CLHEP::m, -180*CLHEP::m);
  }

  m_particleGun->SetParticleDefinition(gamma_particle);
  m_particleGun->SetParticleMomentumDirection(GammaMomDir);
  m_particleGun->SetParticleEnergy(GammaE);
  m_particleGun->SetParticlePosition(VertexPos + target_pos);
  m_particleGun->GeneratePrimaryVertex(anEvent);
}

//_____________________________________________________________________________
// 6374: mono-energetic gamma source from the Li target volume.
// x/y are sampled from truncated Gaussian beam profiles and z is sampled
// uniformly over the requested full target length.
void
S2SPrimaryGeneratorAction::GenerateE63_LiProfileMonoGamma(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);
  static const auto particle = particleTable->FindParticle("gamma");
  static const auto pdg = particle->GetPDGEncoding();
  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_half = sizeMan.GetSize("Target")*mm/2.;

  static G4double gamma_e = 0.;
  static G4double sigma_x = 0.;
  static G4double sigma_y = 0.;
  static G4double length_z = 0.;
  static G4double source_dx = 0.;
  static G4double source_dy = 0.;
  static G4double source_dz = 0.;
  static G4bool loaded = false;
  if(!loaded){
    gamma_e = confMan.Get<G4double>("LiGammaEnergy")*CLHEP::MeV;
    if(gamma_e <= 0.) gamma_e = 1.1*CLHEP::MeV;
    sigma_x = confMan.Get<G4double>("LiGammaSigmaX")*mm;
    sigma_y = confMan.Get<G4double>("LiGammaSigmaY")*mm;
    length_z = confMan.Get<G4double>("LiGammaLengthZ")*mm;
    source_dx = confMan.Get<G4double>("LiGammaSourceDX")*mm;
    source_dy = confMan.Get<G4double>("LiGammaSourceDY")*mm;
    source_dz = confMan.Get<G4double>("LiGammaSourceDZ")*mm;
    if(sigma_x < 0.) sigma_x = 0.;
    if(sigma_y < 0.) sigma_y = 0.;
    if(length_z <= 0.) length_z = 2.*target_half.z();
    G4cout << "GenerateE63_LiProfileMonoGamma: E=" << gamma_e/CLHEP::MeV
           << " MeV sigma=(" << sigma_x/mm << ", " << sigma_y/mm << ") mm"
           << " length_z=" << length_z/mm << " mm"
           << " source_offset=(" << source_dx/mm << ", "
           << source_dy/mm << ", " << source_dz/mm << ") mm"
           << " target_half=" << target_half/mm << " mm"
           << G4endl;
    loaded = true;
  }

  const G4double z_half = std::min(0.5*length_z, target_half.z());
  G4ThreeVector local(source_dx, source_dy, source_dz);
  for(G4int i=0; i<10000; ++i){
    const G4double x = source_dx + ((sigma_x > 0.)
      ? G4RandGauss::shoot(0., sigma_x) : 0.);
    const G4double y = source_dy + ((sigma_y > 0.)
      ? G4RandGauss::shoot(0., sigma_y) : 0.);
    const G4double z = source_dz + G4RandFlat::shoot(-z_half, z_half);
    if(std::abs(x) <= target_half.x()
       && std::abs(y) <= target_half.y()
       && std::abs(z) <= target_half.z()){
      local.set(x, y, z);
      break;
    }
  }
  local.setX(std::clamp(local.x(), -target_half.x(), target_half.x()));
  local.setY(std::clamp(local.y(), -target_half.y(), target_half.y()));
  local.setZ(std::clamp(local.z(), -target_half.z(), target_half.z()));

  const G4ThreeVector vertex = target_pos + local;
  const G4LorentzVector vertex_lv(vertex, 0.);
  const G4ThreeVector dir = G4RandomDirection();
  const G4LorentzVector p4(dir*gamma_e, gamma_e);

  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(dir);
  m_particleGun->SetParticleEnergy(gamma_e);
  m_particleGun->SetParticlePosition(vertex);
  m_particleGun->GeneratePrimaryVertex(anEvent);

  anaMan.SetPrimaryParticle(0, pdg, p4, vertex_lv);
  RecordGeneratedParticle(GenBranch::kDecGamma,
                          GenBranch::ParticleId::None,
                          pdg,
                          p4,
                          vertex_lv);
  anaMan.SetPrimaryData(vertex.x(), vertex.y(), vertex.z(),
                        0., 0., dir.phi(), dir.theta(),
                        gamma_e, 0., pdg);
}

//_____________________________________________________________________________
// 9000 proton beam
void S2SPrimaryGeneratorAction::GenerateProton(G4Event* anEvent)
{
  static const G4int n_particle = 1;
  m_particleGun = new G4ParticleGun(n_particle);

  const auto proton = particleTable->FindParticle("proton");
  const auto pdg = proton->GetPDGEncoding();
  const G4double proton_mass = proton->GetPDGMass();
  const auto& target_pos = geomMan.GetGlobalPosition("Target") * mm;
  const G4double momentum = confMan.Get<G4double>("Momentum") * CLHEP::GeV;

  G4LorentzVector p(0, 0, momentum,
                    std::sqrt(momentum * momentum + proton_mass * proton_mass));
  G4LorentzVector v(target_pos, 0);

  m_particleGun->SetParticleDefinition(proton);
  m_particleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
  m_particleGun->SetParticleEnergy(p.e() - proton_mass);
  m_particleGun->SetParticlePosition(target_pos);
  m_particleGun->GeneratePrimaryVertex(anEvent);

  anaMan.SetPrimaryParticle(0, pdg, p, v);
  anaMan.SetPrimaryData(target_pos.x(), target_pos.y(), target_pos.z(),
                        0., 0., 0., 0., momentum, momentum, 9999);
}

//_____________________________________________________________________________
// 9001 E90 SigmaN Cusp Scattering
void S2SPrimaryGeneratorAction::GenerateSigmaNCusp(G4Event* anEvent)
{
  m_particleGun = new G4ParticleGun();

  // Particle definitions
  const auto pTable = G4ParticleTable::GetParticleTable();
  const auto kaon_minus = pTable->FindParticle("kaon-");
  const auto pi_minus   = pTable->FindParticle("pi-");
  const auto proton     = pTable->FindParticle("proton");
  const auto deuteron   = pTable->FindParticle("deuteron");
  const auto lambda     = pTable->FindParticle("lambda");
  const auto sigma_plus = pTable->FindParticle("sigma+");
  const auto neutron    = pTable->FindParticle("neutron");

  const G4double M_Kaon     = kaon_minus->GetPDGMass(); // MeV
  const G4double M_PiM      = pi_minus->GetPDGMass();   // MeV
  const G4double M_Proton   = proton->GetPDGMass();     // MeV
  const G4double M_Deuteron = deuteron->GetPDGMass();   // MeV
  const G4double M_Lambda   = lambda->GetPDGMass();     // MeV
  const G4double M_SigmaP   = sigma_plus->GetPDGMass(); // MeV
  const G4double M_Neutron  = neutron->GetPDGMass();    // MeV
  const auto pdg_pi_minus = pi_minus->GetPDGEncoding();
  const auto pdg_lambda = lambda->GetPDGEncoding();
  const auto pdg_proton = proton->GetPDGEncoding();

  // ΣN cusp sampling window
  const G4double scat_a = 2.06;
  const G4double scat_b = 4.64;
  const G4double threshold_mass = M_SigmaP + M_Neutron; // MeV
  const G4double mass_min = threshold_mass - 100. * CLHEP::MeV; // MeV
  const G4double mass_max = threshold_mass + 100. * CLHEP::MeV; // MeV

  G4double fcusp_max = 0.;
  for(int i=0; i<200; ++i){
    double mass = mass_min + (mass_max-mass_min)*i/199.; // MeV
    double val = FdComplex::f_single(mass/CLHEP::MeV, scat_a, scat_b, 186.);
    if(val > fcusp_max) fcusp_max = val;
  }

  while(true){

    // Sample ΣN (X) mass
    G4double CuspM = 0.; // MeV
    while(true){
      const G4double MM = G4RandFlat::shoot(mass_min, mass_max); // MeV
      const G4double ds_MM = FdComplex::f_single(MM/CLHEP::MeV, scat_a, scat_b, 186.);
      const G4double Rand = G4RandFlat::shoot(0., fcusp_max * 1.3);
      if(Rand <= ds_MM){
        CuspM = MM; // MeV
        break;
      }
    }

    // K- + d -> π- + X
    const G4double beam_mom_mean = 1.4 * CLHEP::GeV;   // MeV/c
    const G4double beam_mom_sigma = 0.0 * CLHEP::GeV;; // MeV/c
    const G4double beam_mom = G4RandGauss::shoot(beam_mom_mean, beam_mom_sigma) / CLHEP::GeV; // GeV/c
    TLorentzVector beam_lv(0, 0, beam_mom, sqrt(beam_mom*beam_mom + (M_Kaon/CLHEP::GeV)*(M_Kaon/CLHEP::GeV))); // GeV
    TLorentzVector target_lv(0., 0., 0., M_Deuteron/CLHEP::GeV); // GeV
    TLorentzVector W = beam_lv + target_lv; // GeV

    Double_t masses[2] = { M_PiM/CLHEP::GeV, CuspM/CLHEP::GeV }; // GeV
    TGenPhaseSpace event;

    if (!event.SetDecay(W, 2, masses) || event.Generate() == 0) continue;

    TLorentzVector *pi_lv = event.GetDecay(0); // GeV
    TLorentzVector *X_lv  = event.GetDecay(1); // GeV

    if (pi_lv->Theta() > 15. * CLHEP::deg) continue;

    // Target vertex
    static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
    static const auto& target_size = sizeMan.GetSize("Target");

    const G4double target_radius = target_size.y() / 2.0 * mm;
    const G4double target_height_half = target_size.z() / 2.0 * mm;
    G4double r_vtx_x, r_vtx_y, r_vtx_z;
    const G4double beam_x_mean = -14.5 * CLHEP::mm;
    const G4double beam_y_mean = 1.3 * CLHEP::mm;
    const G4double beam_x_sigma = 24.3 * CLHEP::mm;
    const G4double beam_y_sigma = 4.2 * CLHEP::mm;
    while (true) {
        r_vtx_x = G4RandGauss::shoot(beam_x_mean, beam_x_sigma);
        r_vtx_y = G4RandGauss::shoot(beam_y_mean, beam_y_sigma);
        r_vtx_z = G4RandFlat::shoot(-target_radius, target_radius);
        if (r_vtx_x*r_vtx_x + r_vtx_z*r_vtx_z < target_radius*target_radius && r_vtx_y < target_height_half && r_vtx_y > -target_height_half) break;
    }
    const G4ThreeVector primary_vertex_pos = target_pos + G4ThreeVector(r_vtx_x, r_vtx_y, r_vtx_z);
    const auto vertex_lv = G4LorentzVector(primary_vertex_pos, 0);

    // Primary π-
    m_particleGun->SetParticleDefinition(pi_minus);
    m_particleGun->SetParticleMomentumDirection(G4ThreeVector(pi_lv->Px(), pi_lv->Py(), pi_lv->Pz()).unit());
    m_particleGun->SetParticleEnergy(pi_lv->E()*CLHEP::GeV - M_PiM); // MeV
    m_particleGun->SetParticlePosition(primary_vertex_pos); // mm
    m_particleGun->GeneratePrimaryVertex(anEvent);

    // X -> Λ + p
    Double_t masses_X[2] = { M_Lambda/CLHEP::GeV, M_Proton/CLHEP::GeV }; // GeV
    TGenPhaseSpace event_X;

    if (!event_X.SetDecay(*X_lv, 2, masses_X) || event_X.Generate() == 0) continue;

    TLorentzVector *lambda_lv = event_X.GetDecay(0); // GeV
    TLorentzVector *p_lv      = event_X.GetDecay(1); // GeV

    m_particleGun->SetParticleDefinition(lambda);
    m_particleGun->SetParticleMomentumDirection(G4ThreeVector(lambda_lv->Px(), lambda_lv->Py(), lambda_lv->Pz()).unit());
    m_particleGun->SetParticleEnergy(lambda_lv->E()*CLHEP::GeV - M_Lambda); // MeV
    m_particleGun->SetParticlePosition(primary_vertex_pos); // mm
    m_particleGun->GeneratePrimaryVertex(anEvent);

    m_particleGun->SetParticleDefinition(proton);
    m_particleGun->SetParticleMomentumDirection(G4ThreeVector(p_lv->Px(), p_lv->Py(), p_lv->Pz()).unit());
    m_particleGun->SetParticleEnergy(p_lv->E()*CLHEP::GeV - M_Proton); // MeV
    m_particleGun->SetParticlePosition(primary_vertex_pos); // mm
    m_particleGun->GeneratePrimaryVertex(anEvent);

    RecordGeneratedParticle(GenBranch::kPrimPi,
                            GenBranch::ParticleId::None,
                            pdg_pi_minus,
                            *pi_lv,
                            vertex_lv);
    RecordGeneratedParticle(GenBranch::kLambda,
                            GenBranch::ParticleId::None,
                            pdg_lambda,
                            *lambda_lv,
                            vertex_lv);
    RecordGeneratedParticle(GenBranch::kPrimProt,
                            GenBranch::ParticleId::None,
                            pdg_proton,
                            *p_lv,
                            vertex_lv);
    {
      Double_t decay_masses[2] = { M_Proton/CLHEP::GeV, M_PiM/CLHEP::GeV };
      TGenPhaseSpace lambda_decay;
      if(lambda_decay.SetDecay(*lambda_lv, 2, decay_masses) && lambda_decay.Generate()!=0){
        const auto decay_p = lambda_decay.GetDecay(0);
        const auto decay_pi = lambda_decay.GetDecay(1);
        RecordGeneratedParticle(GenBranch::kDecProt,
                                GenBranch::ParticleId::Lambda,
                                pdg_proton,
                                *decay_p,
                                vertex_lv);
        RecordGeneratedParticle(GenBranch::kDecPi,
                                GenBranch::ParticleId::Lambda,
                                pdg_pi_minus,
                                *decay_pi,
                                vertex_lv);
      }
    }

    anaMan.SetPrimaryParticle(0, pdg_pi_minus,
                              ToG4Lorentz(*pi_lv),
                              vertex_lv);
    const auto theta_pi = pi_lv->Theta();
    const auto phi_pi = pi_lv->Phi();
    const auto u0 = std::tan(theta_pi) * std::cos(phi_pi);
    const auto v0 = std::tan(theta_pi) * std::sin(phi_pi);
    const auto p0 = pi_lv->P() * CLHEP::GeV; // MeV/c
    const auto pB = beam_mom * CLHEP::GeV; // beam K- momentum (MeV/c)
    anaMan.SetPrimaryData(primary_vertex_pos.x(), primary_vertex_pos.y(), primary_vertex_pos.z(),
                          u0, v0, phi_pi, theta_pi, p0, pB, 9999); // x0,y0,z0,u0,v0,phi,theta,p0,pB,ParIdNb

    break;
  }
}

//_____________________________________________________________________________
// 9002 E90 QF Lambda reaction
void S2SPrimaryGeneratorAction::GenerateQFLambda(G4Event* anEvent)
{
  m_particleGun = new G4ParticleGun();

  // Particle definitions
  const auto pTable = G4ParticleTable::GetParticleTable();
  const auto kaon_minus = pTable->FindParticle("kaon-");
  const auto pi_minus   = pTable->FindParticle("pi-");
  const auto proton     = pTable->FindParticle("proton");
  const auto lambda     = pTable->FindParticle("lambda");
  const auto deuteron   = pTable->FindParticle("deuteron");

  const G4double M_Kaon     = kaon_minus->GetPDGMass(); // MeV
  const G4double M_PiM      = pi_minus->GetPDGMass();   // MeV
  const G4double M_Proton   = proton->GetPDGMass();     // MeV
  const G4double M_Lambda   = lambda->GetPDGMass();     // MeV
  const G4double M_Deuteron = deuteron->GetPDGMass();   // MeV
  const auto pdg_pi_minus = pi_minus->GetPDGEncoding();
  const auto pdg_lambda = lambda->GetPDGEncoding();
  const auto pdg_proton = proton->GetPDGEncoding();

  while(true){
    // Fermi motion for the struck neutron
    TVector3 p_fermi_vec = FermiMotion::GetMomentum(); // GeV/c
    G4double E_spectator_fermi = sqrt((M_Proton/CLHEP::GeV)*(M_Proton/CLHEP::GeV) + p_fermi_vec.Mag2()); // GeV
    G4double E_target = (M_Deuteron/CLHEP::GeV) - E_spectator_fermi; // GeV
    TLorentzVector target_lv(p_fermi_vec, E_target); // GeV

    // K- + "n" -> Λ + π-
    G4double beam_mom_mean = 1.4 * CLHEP::GeV;  // MeV/c
    G4double beam_mom_sigma = 0.0 * CLHEP::GeV; // MeV/c
    G4double beam_mom_gev = G4RandGauss::shoot(beam_mom_mean, beam_mom_sigma)/CLHEP::GeV; // GeV/c
    TLorentzVector beam_lv(0, 0, beam_mom_gev, sqrt(beam_mom_gev*beam_mom_gev + (M_Kaon/CLHEP::GeV)*(M_Kaon/CLHEP::GeV))); // GeV

    TLorentzVector W = beam_lv + target_lv;

    Double_t masses[2] = { M_Lambda/CLHEP::GeV, M_PiM/CLHEP::GeV }; // GeV
    TGenPhaseSpace event;

    if (!event.SetDecay(W, 2, masses) || event.Generate() == 0) continue;

    TLorentzVector *lambda_lv = event.GetDecay(0); // GeV
    TLorentzVector *pi_lv     = event.GetDecay(1); // GeV

    if (pi_lv->Theta() > 15. * CLHEP::deg) continue;

    // Target vertex
    static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm; // mm
    static const auto& target_size = sizeMan.GetSize("Target"); // mm
    G4double target_radius = target_size.y() / 2.0 * mm;        // mm
    G4double target_height_half = target_size.z() / 2.0 * mm;   // mm
    G4double r_vtx_x, r_vtx_y, r_vtx_z;       // mm
    G4double beam_x_mean = -14.5 * CLHEP::mm; // mm
    G4double beam_y_mean = 1.3 * CLHEP::mm;   // mm
    G4double beam_x_sigma = 24.3 * CLHEP::mm; // mm
    G4double beam_y_sigma = 4.2 * CLHEP::mm;  // mm
    while (true) {
        r_vtx_x = G4RandGauss::shoot(0, beam_x_sigma); // mm
        r_vtx_y = G4RandGauss::shoot(0, beam_y_sigma); // mm
        r_vtx_z = G4RandFlat::shoot(-target_radius, target_radius); // mm
        if (r_vtx_x*r_vtx_x + r_vtx_z*r_vtx_z < target_radius*target_radius && r_vtx_y < target_height_half && r_vtx_y > -target_height_half) break;
    }
    G4ThreeVector primary_vertex_pos = target_pos + G4ThreeVector(r_vtx_x, r_vtx_y, r_vtx_z); // mm
    const auto vertex_lv = G4LorentzVector(primary_vertex_pos, 0);

    // Primary π-
    m_particleGun->SetParticleDefinition(pi_minus);
    m_particleGun->SetParticleMomentumDirection(G4ThreeVector(pi_lv->Px(), pi_lv->Py(), pi_lv->Pz()).unit());
    m_particleGun->SetParticleEnergy(pi_lv->E()*CLHEP::GeV - M_PiM); // MeV
    m_particleGun->SetParticlePosition(primary_vertex_pos);          // mm
    m_particleGun->GeneratePrimaryVertex(anEvent);

    // Λ
    m_particleGun->SetParticleDefinition(lambda);
    m_particleGun->SetParticleMomentumDirection(G4ThreeVector(lambda_lv->Px(), lambda_lv->Py(), lambda_lv->Pz()).unit());
    m_particleGun->SetParticleEnergy(lambda_lv->E()*CLHEP::GeV - M_Lambda); // MeV
    m_particleGun->SetParticlePosition(primary_vertex_pos);                 // mm
    m_particleGun->GeneratePrimaryVertex(anEvent);

    // Spectator proton
    const TVector3 spectator_vec = -p_fermi_vec; // GeV/c
    const G4ThreeVector spectator_mom(spectator_vec.X(), spectator_vec.Y(), spectator_vec.Z()); // GeV/c

    G4double E_spectator = sqrt(M_Proton*M_Proton + spectator_mom.mag2()*(CLHEP::GeV*CLHEP::GeV)); // MeV
    m_particleGun->SetParticleDefinition(proton);
    m_particleGun->SetParticleMomentumDirection(spectator_mom.unit());
    m_particleGun->SetParticleEnergy(E_spectator - M_Proton); // MeV
    m_particleGun->SetParticlePosition(primary_vertex_pos);   // mm
    m_particleGun->GeneratePrimaryVertex(anEvent);

    RecordGeneratedParticle(GenBranch::kPrimPi,
                            GenBranch::ParticleId::None,
                            pdg_pi_minus,
                            *pi_lv,
                            vertex_lv);
    RecordGeneratedParticle(GenBranch::kLambda,
                            GenBranch::ParticleId::None,
                            pdg_lambda,
                            *lambda_lv,
                            vertex_lv);
    RecordGeneratedParticle(GenBranch::kPrimProt,
                            GenBranch::ParticleId::None,
                            pdg_proton,
                            spectator_vec,
                            E_spectator,
                            vertex_lv);

    {
      Double_t decay_masses[2] = { M_Proton/CLHEP::GeV, M_PiM/CLHEP::GeV };
      TGenPhaseSpace lambda_decay;
      if(lambda_decay.SetDecay(*lambda_lv, 2, decay_masses) && lambda_decay.Generate()!=0){
        const auto decay_p = lambda_decay.GetDecay(0);
        const auto decay_pi = lambda_decay.GetDecay(1);
        RecordGeneratedParticle(GenBranch::kDecProt,
                                GenBranch::ParticleId::Lambda,
                                pdg_proton,
                                *decay_p,
                                vertex_lv);
        RecordGeneratedParticle(GenBranch::kDecPi,
                                GenBranch::ParticleId::Lambda,
                                pdg_pi_minus,
                                *decay_pi,
                                vertex_lv);
      }
    }

    anaMan.SetPrimaryParticle(0, pdg_pi_minus,
                              ToG4Lorentz(*pi_lv),
                              vertex_lv);
    const auto theta_pi = pi_lv->Theta();
    const auto phi_pi = pi_lv->Phi();
    const auto u0 = std::tan(theta_pi) * std::cos(phi_pi);
    const auto v0 = std::tan(theta_pi) * std::sin(phi_pi);
    const auto p0 = pi_lv->P() * CLHEP::GeV; // MeV/c
    const auto pB = beam_mom_gev * CLHEP::GeV; // beam K- momentum (MeV/c)
    anaMan.SetPrimaryData(primary_vertex_pos.x(), primary_vertex_pos.y(), primary_vertex_pos.z(),
                          u0, v0, phi_pi, theta_pi, p0, pB, 9999); // x0,y0,z0,u0,v0,phi,theta,p0,pB,ParIdNb
    break;
  }
}

//_____________________________________________________________________________
// 9003 E90 QF Sigma0 reaction
void S2SPrimaryGeneratorAction::GenerateQFSigmaZ(G4Event* anEvent)
{
  m_particleGun = new G4ParticleGun();

  // Particle definitions
  const auto pTable = G4ParticleTable::GetParticleTable();
  const auto kaon_minus = pTable->FindParticle("kaon-");
  const auto pi_minus   = pTable->FindParticle("pi-");
  const auto proton     = pTable->FindParticle("proton");
  const auto sigma0     = pTable->FindParticle("sigma0");
  const auto lambda     = pTable->FindParticle("lambda");
  const auto deuteron   = pTable->FindParticle("deuteron");

  const G4double M_Kaon     = kaon_minus->GetPDGMass(); // MeV
  const G4double M_PiM      = pi_minus->GetPDGMass();   // MeV
  const G4double M_Proton   = proton->GetPDGMass();     // MeV
  const G4double M_Sigma0   = sigma0->GetPDGMass();     // MeV
  const G4double M_Lambda   = lambda->GetPDGMass();     // MeV
  const G4double M_Deuteron = deuteron->GetPDGMass();   // MeV
  const auto pdg_pi_minus = pi_minus->GetPDGEncoding();
  const auto pdg_sigma0 = sigma0->GetPDGEncoding();
  const auto pdg_lambda = lambda->GetPDGEncoding();
  const auto pdg_proton = proton->GetPDGEncoding();

  while(true){
    // Fermi motion for the struck neutron
    TVector3 p_fermi_vec = FermiMotion::GetMomentum(); // GeV/c
    G4double E_spectator_fermi = sqrt((M_Proton/CLHEP::GeV)*(M_Proton/CLHEP::GeV) + p_fermi_vec.Mag2()); // GeV
    G4double E_target = (M_Deuteron/CLHEP::GeV) - E_spectator_fermi; // GeV
    TLorentzVector target_lv(p_fermi_vec, E_target); // GeV

    // K- + "n" -> Σ0 + π-
    G4double beam_mom_mean = 1.4 * CLHEP::GeV;  // MeV/c
    G4double beam_mom_sigma = 0.0 * CLHEP::GeV; // MeV/c
    G4double beam_mom_gev = G4RandGauss::shoot(beam_mom_mean, beam_mom_sigma)/CLHEP::GeV; // GeV/c
    TLorentzVector beam_lv(0, 0, beam_mom_gev, sqrt(beam_mom_gev*beam_mom_gev + (M_Kaon/CLHEP::GeV)*(M_Kaon/CLHEP::GeV))); // GeV
    TLorentzVector W = beam_lv + target_lv;

    Double_t masses[2] = { M_Sigma0/CLHEP::GeV, M_PiM/CLHEP::GeV }; // GeV
    TGenPhaseSpace event;

    if (!event.SetDecay(W, 2, masses) || event.Generate() == 0) continue;

    TLorentzVector *sigma0_lv = event.GetDecay(0); // GeV
    TLorentzVector *pi_lv     = event.GetDecay(1); // GeV

    if (pi_lv->Theta() > 15. * CLHEP::deg) continue;

    // Target vertex
    static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm; // mm
    static const auto& target_size = sizeMan.GetSize("Target"); // mm
    G4double target_radius = target_size.y() / 2.0 * mm; // mm
    G4double target_height_half = target_size.z() / 2.0 * mm; // mm
    G4double r_vtx_x, r_vtx_y, r_vtx_z;       // mm
    G4double beam_x_mean = -14.5 * CLHEP::mm; // mm
    G4double beam_y_mean =   1.3 * CLHEP::mm; // mm
    G4double beam_x_sigma = 24.3 * CLHEP::mm; // mm
    G4double beam_y_sigma =  4.2 * CLHEP::mm; // mm
    while (true) {
        r_vtx_x = G4RandGauss::shoot(0, beam_x_sigma); // mm
        r_vtx_y = G4RandGauss::shoot(0, beam_y_sigma); // mm
        r_vtx_z = G4RandFlat::shoot(-target_radius, target_radius); // mm
        if (r_vtx_x*r_vtx_x + r_vtx_z*r_vtx_z < target_radius*target_radius && r_vtx_y < target_height_half && r_vtx_y > -target_height_half) break;
    }
    G4ThreeVector primary_vertex_pos = target_pos + G4ThreeVector(r_vtx_x, r_vtx_y, r_vtx_z);
    const auto vertex_lv = G4LorentzVector(primary_vertex_pos, 0);

    // Primary π-
    m_particleGun->SetParticleDefinition(pi_minus);
    m_particleGun->SetParticleMomentumDirection(G4ThreeVector(pi_lv->Px(), pi_lv->Py(), pi_lv->Pz()).unit());
    m_particleGun->SetParticleEnergy(pi_lv->E()*CLHEP::GeV - M_PiM); // MeV
    m_particleGun->SetParticlePosition(primary_vertex_pos); // mm
    m_particleGun->GeneratePrimaryVertex(anEvent);

    // Σ0
    m_particleGun->SetParticleDefinition(sigma0);
    m_particleGun->SetParticleMomentumDirection(G4ThreeVector(sigma0_lv->Px(), sigma0_lv->Py(), sigma0_lv->Pz()).unit());
    m_particleGun->SetParticleEnergy(sigma0_lv->E()*CLHEP::GeV - M_Sigma0); // MeV
    m_particleGun->SetParticlePosition(primary_vertex_pos); // mm
    m_particleGun->GeneratePrimaryVertex(anEvent);

    // Spectator proton
    const TVector3 spectator_vec = -p_fermi_vec; // GeV/c
    const G4ThreeVector spectator_mom(spectator_vec.X(), spectator_vec.Y(), spectator_vec.Z()); // GeV/c

    G4double E_spectator = sqrt(M_Proton*M_Proton + spectator_mom.mag2()*(CLHEP::GeV*CLHEP::GeV)); // MeV
    m_particleGun->SetParticleDefinition(proton);
    m_particleGun->SetParticleMomentumDirection(spectator_mom.unit());
    m_particleGun->SetParticleEnergy(E_spectator - M_Proton); // MeV
    m_particleGun->SetParticlePosition(primary_vertex_pos); // mm
    m_particleGun->GeneratePrimaryVertex(anEvent);

    TLorentzVector lambda_from_sigma = *sigma0_lv;
    {
      Double_t sigma_decay_masses[2] = { M_Lambda/CLHEP::GeV, 0.0 };
      TGenPhaseSpace sigma_decay;
      if(sigma_decay.SetDecay(*sigma0_lv, 2, sigma_decay_masses) && sigma_decay.Generate()!=0){
        lambda_from_sigma = *sigma_decay.GetDecay(0);
      }
    }

    RecordGeneratedParticle(GenBranch::kPrimPi,
                            GenBranch::ParticleId::None,
                            pdg_pi_minus,
                            *pi_lv,
                            vertex_lv);
    RecordGeneratedParticle(GenBranch::kSigma0,
                            GenBranch::ParticleId::None,
                            pdg_sigma0,
                            *sigma0_lv,
                            vertex_lv);
    RecordGeneratedParticle(GenBranch::kLambda,
                            GenBranch::ParticleId::Sigma0,
                            pdg_lambda,
                            lambda_from_sigma,
                            vertex_lv);
    RecordGeneratedParticle(GenBranch::kPrimProt,
                            GenBranch::ParticleId::None,
                            pdg_proton,
                            spectator_vec,
                            E_spectator,
                            vertex_lv);
    {
      Double_t lambda_decay_masses[2] = { M_Proton/CLHEP::GeV, M_PiM/CLHEP::GeV };
      TGenPhaseSpace lambda_decay;
      if(lambda_decay.SetDecay(lambda_from_sigma, 2, lambda_decay_masses) && lambda_decay.Generate()!=0){
        const auto decay_p = lambda_decay.GetDecay(0);
        const auto decay_pi = lambda_decay.GetDecay(1);
        RecordGeneratedParticle(GenBranch::kDecProt,
                                GenBranch::ParticleId::Lambda,
                                pdg_proton,
                                *decay_p,
                                vertex_lv);
        RecordGeneratedParticle(GenBranch::kDecPi,
                                GenBranch::ParticleId::Lambda,
                                pdg_pi_minus,
                                *decay_pi,
                                vertex_lv);
      }
    }

    anaMan.SetPrimaryParticle(0, pdg_pi_minus,
                              ToG4Lorentz(*pi_lv),
                              vertex_lv);
    const auto theta_pi = pi_lv->Theta();
    const auto phi_pi = pi_lv->Phi();
    const auto u0 = std::tan(theta_pi) * std::cos(phi_pi);
    const auto v0 = std::tan(theta_pi) * std::sin(phi_pi);
    const auto p0 = pi_lv->P() * CLHEP::GeV; // MeV/c
    const auto pB = beam_mom_gev * CLHEP::GeV; // beam K- momentum (MeV/c)
    anaMan.SetPrimaryData(primary_vertex_pos.x(), primary_vertex_pos.y(), primary_vertex_pos.z(),
                          u0, v0, phi_pi, theta_pi, p0, pB, 9999); // x0,y0,z0,u0,v0,phi,theta,p0,pB,ParIdNb
    break;
  }
}

//_____________________________________________________________________________
// 9004 E90 QF Sigma+ reaction
void S2SPrimaryGeneratorAction::GenerateQFSigmaP(G4Event* anEvent)
{
  m_particleGun = new G4ParticleGun();

  // Particle definitions
  const auto pTable = G4ParticleTable::GetParticleTable();
  const auto kaon_minus = pTable->FindParticle("kaon-");
  const auto pi_minus   = pTable->FindParticle("pi-");
  const auto pi_zero    = pTable->FindParticle("pi0");
  const auto sigma_plus = pTable->FindParticle("sigma+");
  const auto proton     = pTable->FindParticle("proton");
  const auto neutron    = pTable->FindParticle("neutron");
  const auto deuteron   = pTable->FindParticle("deuteron");

  const G4double M_Kaon     = kaon_minus->GetPDGMass(); // MeV
  const G4double M_PiM      = pi_minus->GetPDGMass(); // MeV
  const G4double M_SigmaP   = sigma_plus->GetPDGMass(); // MeV
  const G4double M_Pi0      = pi_zero->GetPDGMass();    // MeV
  const G4double M_Proton   = proton->GetPDGMass(); // MeV
  const G4double M_Neutron  = neutron->GetPDGMass(); // MeV
  const G4double M_Deuteron = deuteron->GetPDGMass(); // MeV
  const auto pdg_pi_minus = pi_minus->GetPDGEncoding();
  const auto pdg_sigma_plus = sigma_plus->GetPDGEncoding();
  const auto pdg_proton = proton->GetPDGEncoding();
  const auto pdg_neutron = neutron->GetPDGEncoding();
  const auto pdg_pi_zero = pi_zero->GetPDGEncoding();

  while(true){

    // Fermi motion for the struck proton
    TVector3 p_fermi_vec = FermiMotion::GetMomentum();
    G4double E_spectator_fermi = sqrt((M_Neutron/CLHEP::GeV)*(M_Neutron/CLHEP::GeV) + p_fermi_vec.Mag2()); // GeV
    G4double E_target = (M_Deuteron/CLHEP::GeV) - E_spectator_fermi; // GeV
    TLorentzVector target_lv(p_fermi_vec, E_target); //GeV

    // K- + "p" -> Σ+ + π-
    G4double beam_mom_mean = 1.4 * CLHEP::GeV; // MeV/c
    G4double beam_mom_sigma = 0.0 * CLHEP::GeV; // MeV/c
    G4double beam_mom_gev = G4RandGauss::shoot(beam_mom_mean, beam_mom_sigma)/CLHEP::GeV; // GeV/c
    TLorentzVector beam_lv(0, 0, beam_mom_gev, sqrt(beam_mom_gev*beam_mom_gev + (M_Kaon/CLHEP::GeV)*(M_Kaon/CLHEP::GeV))); // GeV
    TLorentzVector W = beam_lv + target_lv;

    Double_t masses[2] = { M_SigmaP/CLHEP::GeV, M_PiM/CLHEP::GeV }; // GeV
    TGenPhaseSpace event;

    if (!event.SetDecay(W, 2, masses) || event.Generate() == 0) continue;

    TLorentzVector *sigma_plus_lv = event.GetDecay(0); // GeV
    TLorentzVector *pi_lv     = event.GetDecay(1); // GeV

    if (pi_lv->Theta() > 15. * CLHEP::deg) continue;

    // Target vertex
    static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm; // mm
    static const auto& target_size = sizeMan.GetSize("Target"); // mm
    G4double target_radius = target_size.y() / 2.0 * mm; // mm
    G4double target_height_half = target_size.z() / 2.0 * mm; // mm
    G4double r_vtx_x, r_vtx_y, r_vtx_z;       // mm
    G4double beam_x_mean = -14.5 * CLHEP::mm; // mm
    G4double beam_y_mean = 1.3 * CLHEP::mm;   // mm
    G4double beam_x_sigma = 24.3 * CLHEP::mm; // mm
    G4double beam_y_sigma = 4.2 * CLHEP::mm;  // mm
    while (true) {
        r_vtx_x = G4RandGauss::shoot(0, beam_x_sigma); // mm
        r_vtx_y = G4RandGauss::shoot(0, beam_y_sigma); // mm
        r_vtx_z = G4RandFlat::shoot(-target_radius, target_radius); // mm
        if (r_vtx_x*r_vtx_x + r_vtx_z*r_vtx_z < target_radius*target_radius && r_vtx_y < target_height_half && r_vtx_y > -target_height_half) break;
    }
    G4ThreeVector primary_vertex_pos = target_pos + G4ThreeVector(r_vtx_x, r_vtx_y, r_vtx_z);
    const auto vertex_lv = G4LorentzVector(primary_vertex_pos, 0);

    // Primary π-
    m_particleGun->SetParticleDefinition(pi_minus);
    m_particleGun->SetParticleMomentumDirection(G4ThreeVector(pi_lv->Px(), pi_lv->Py(), pi_lv->Pz()).unit());
    m_particleGun->SetParticleEnergy(pi_lv->E()*CLHEP::GeV - M_PiM); // MeV
    m_particleGun->SetParticlePosition(primary_vertex_pos); // mm
    m_particleGun->GeneratePrimaryVertex(anEvent);

    // Σ+
    m_particleGun->SetParticleDefinition(sigma_plus);
    m_particleGun->SetParticleMomentumDirection(G4ThreeVector(sigma_plus_lv->Px(), sigma_plus_lv->Py(), sigma_plus_lv->Pz()).unit());
    m_particleGun->SetParticleEnergy(sigma_plus_lv->E()*CLHEP::GeV - M_SigmaP); // MeV
    m_particleGun->SetParticlePosition(primary_vertex_pos); // mm
    m_particleGun->GeneratePrimaryVertex(anEvent);

    // Spectator neutron
    const TVector3 spectator_vec = -p_fermi_vec; // GeV/c
    const G4ThreeVector spectator_mom(spectator_vec.X(), spectator_vec.Y(), spectator_vec.Z()); // GeV/c

    G4double E_spectator = sqrt(M_Neutron*M_Neutron + spectator_mom.mag2()*(CLHEP::GeV*CLHEP::GeV)); // MeV
    m_particleGun->SetParticleDefinition(neutron);
    m_particleGun->SetParticleMomentumDirection(spectator_mom.unit());
    m_particleGun->SetParticleEnergy(E_spectator - M_Neutron); // MeV
    m_particleGun->SetParticlePosition(primary_vertex_pos); // mm
    m_particleGun->GeneratePrimaryVertex(anEvent);

    // Σ+ -> pπ0
    RecordGeneratedParticle(GenBranch::kPrimPi,
                            GenBranch::ParticleId::None,
                            pdg_pi_minus,
                            *pi_lv,
                            vertex_lv);
    RecordGeneratedParticle(GenBranch::kSigmaP,
                            GenBranch::ParticleId::None,
                            pdg_sigma_plus,
                            *sigma_plus_lv,
                            vertex_lv);
    RecordGeneratedParticle(GenBranch::kPrimNeut,
                            GenBranch::ParticleId::None,
                            pdg_neutron,
                            spectator_vec,
                            E_spectator,
                            vertex_lv);
    {
      Double_t decay_masses[2] = { M_Proton/CLHEP::GeV, M_Pi0/CLHEP::GeV };
      TGenPhaseSpace sigma_plus_decay;
      if(sigma_plus_decay.SetDecay(*sigma_plus_lv, 2, decay_masses) && sigma_plus_decay.Generate()!=0){
        const auto decay_p = sigma_plus_decay.GetDecay(0);
        const auto decay_pi = sigma_plus_decay.GetDecay(1);
        RecordGeneratedParticle(GenBranch::kDecProt,
                                GenBranch::ParticleId::SigmaP,
                                pdg_proton,
                                *decay_p,
                                vertex_lv);
        RecordGeneratedParticle(GenBranch::kDecPi,
                                GenBranch::ParticleId::SigmaP,
                                pdg_pi_zero,
                                *decay_pi,
                                vertex_lv);
      }
    }

    anaMan.SetPrimaryParticle(0, pdg_pi_minus,
                              ToG4Lorentz(*pi_lv),
                              vertex_lv);
    const auto theta_pi = pi_lv->Theta();
    const auto phi_pi = pi_lv->Phi();
    const auto u0 = std::tan(theta_pi) * std::cos(phi_pi);
    const auto v0 = std::tan(theta_pi) * std::sin(phi_pi);
    const auto p0 = pi_lv->P() * CLHEP::GeV; // MeV/c
    const auto pB = beam_mom_gev * CLHEP::GeV; // beam K- momentum (MeV/c)
    anaMan.SetPrimaryData(primary_vertex_pos.x(), primary_vertex_pos.y(), primary_vertex_pos.z(),
                          u0, v0, phi_pi, theta_pi, p0, pB, 9999); // x0,y0,z0,u0,v0,phi,theta,p0,pB,ParIdNb
    break;
  }
}

//_____________________________________________________________________________
// Generator 1001 (E10 namespace: E63 uses 63xx, E70 uses 70xx/75xx,
// E90 uses 90xx -- keeping E10's own generators out of those blocks
// avoids ID collisions in this shared repo): K0 production for
// 9Be(K-,pi+) background studies.
//
// K- + p (bound proton, quasi-free m_p-S_p) -> K0bar + n, with the
// K0bar c.m.-frame angle sampled from the MEASURED angular
// distribution (Conforto et al., Nucl. Phys. B105 (1976) 189-221,
// Rutherford-Imperial College collaboration, Legendre A_l
// coefficients at their highest point, p_K-=1.355 GeV/c -- the
// closest available data to our p_beam=1.5 GeV/c; their range does
// not quite reach 1.5 GeV/c). K-p->K0bar n is BACKWARD-peaked in the
// K0bar c.m. angle throughout Conforto's whole 0.960-1.355 GeV/c
// range (backward/forward dsigma/dcosTheta* ratio ~1.5-5.8).
//
// REVISION (2026-07-12, later): rather than hand-computing the
// K0_S/K0_L species mixing, decay channel and flight length
// ourselves in C++ (as the two earlier revisions of this function
// did), the K0bar is generated as an actual Geant4 primary
// ("anti_kaon0") and everything downstream -- the K0_S/K0_L mixing,
// each species' full decay table (pi+pi-, pi0pi0, pi+pi-pi0,
// pi0pi0pi0, Ke3, Kmu3, ... with their real PDG branching ratios and
// matrix elements) and the flight length before decay (from each
// species' real proper lifetime) -- is handled by Geant4's own
// G4Decay physics (conf: DECAY:1) acting on the real lab-frame
// track, through the real geometry/field/materials. This is more
// complete (every channel, not just an importance-sampled pi+
// subset) and removes the need to hand-model the World-volume
// boundary for K0_L's long flight length.
//
// Because the resulting charged decay products are Geant4
// *secondaries* (parent track = the K0bar/K0_S/K0_L, not the
// event's primary vertex), they never satisfy VHitInfo::IsPrimary()
// (which requires track_id==1), so the historical
// trigger_flag[kTOF/kAC1/kWC] logic in S2SAnaManager (which gates on
// IsPrimary()) cannot see them. Acceptance for this generator (and,
// more generally, for all Experiment==10 studies) is therefore
// evaluated via the common e10 downstream-counter criterion in
// S2SAnaManager::EndOfEvent: a charged hit in TOF, a charged hit in
// AC1 exceeding its aerogel Cherenkov threshold (n=1.05), and a
// charged hit in WC, from ANY track -- matching how the real
// hardware trigger works (it does not know which track is "primary"
// or what species it is). See analysis-note.md for the derivation
// and results.
void
S2SPrimaryGeneratorAction::GenerateK0Production(G4Event* anEvent)
{
  m_particleGun = new G4ParticleGun();

  const auto pTable = G4ParticleTable::GetParticleTable();
  const auto anti_kaon0 = pTable->FindParticle("anti_kaon0");
  const auto kaon_minus = pTable->FindParticle("kaon-");
  const auto neutron    = pTable->FindParticle("neutron");

  const G4double M_K0 = anti_kaon0->GetPDGMass()/CLHEP::GeV;  // GeV
  const G4double M_KM = kaon_minus->GetPDGMass()/CLHEP::GeV;  // GeV
  const G4double M_N  = neutron->GetPDGMass()/CLHEP::GeV;     // GeV
  const auto pdg_anti_kaon0 = anti_kaon0->GetPDGEncoding();

  // Bound-proton effective mass (quasi-free): M_target - M_spectator
  // (9Be->8Li+p, or 12C->11B+p for the E05 SKS cross-check).
  // Selected by conf key "K0TargetA" (9 default / 12); missing key
  // reads 0 -> treated as 9.
  const G4int k0_target_a = confMan.Get<G4int>("K0TargetA");
  const G4bool k0_is_c12 = (k0_target_a == 12);
  const G4double M_P_EFF =
    k0_is_c12 ? (11.174862342 - 10.252547297) : 0.9213850010499991;
  // K- beam momentum (GeV/c): optional conf key "K0BeamMom" for a
  // beam-momentum scan (see analysis-note.md 2026-07-13 entries);
  // defaults to 1.5 GeV/c if unset (ConfMan::Get<G4double> silently
  // returns 0. for a missing key, so treat 0 as "not set").
  const G4double k0_beam_mom_conf = confMan.Get<G4double>("K0BeamMom");
  const G4double P_BEAM = (k0_beam_mom_conf > 0.) ? k0_beam_mom_conf : 1.5;

  // Optional Fermi motion of the struck (bound) proton: conf key
  // "K0FermiMotion" (0/1, default 0 = off, backward compatible).
  //
  // Momentum distribution: the exact 1p harmonic-oscillator momentum
  // density |phi~_1p(p)|^2 p^2 ~ p^4 exp(-b~^2 p^2), which is a
  // chi-distribution with 5 degrees of freedom in |p| (scale sigma =
  // hbarc/(b~ sqrt(2)) per Cartesian component), direction isotropic.
  // b~ = 1.716 fm is the CM-CORRECTED 9Be 1p oscillator parameter
  // (neff_calc.py B_9BE_N): the struck nucleon's momentum in the
  // nucleus rest frame IS the nucleon-core relative momentum (they
  // are back to back), and the relative-coordinate wave function
  // carries b~, not the bare mean-field b. <p^2> = (5/2)(hbarc/b~)^2
  // -> p_rms ~ 182 MeV/c, consistent with the |p|~Gauss(163,50)
  // MeV/c model in reaction-kinematics/module/channels.py (~170).
  //
  // Energy: spectator-recoil energy conservation (same form as
  // reaction-kinematics/module/generator.py `_initial_state`): the
  // nucleus at rest splits into the struck proton (+p) and an
  // on-shell 8Li spectator core (-p), so
  //   E_p(p) = M(9Be) - sqrt(p^2 + M(8Li)^2),
  // which reduces to E_p(0) = M(9Be)-M(8Li) = M_P_EFF exactly. An
  // earlier version (2026-07-13) instead used the WRONG
  // sqrt(p^2 + M_P_EFF^2) -- giving the proton its Fermi momentum
  // "for free" (energy INCREASES with p, violating the nuclear
  // energy budget) and a per-component Gaussian of sigma=hbarc/b_bare
  // (a 1s-shaped distribution with the wrong width). Both fixed
  // 2026-07-14; see analysis-note.md and
  // reaction-kinematics/CHANGELOG.md (2026-07-14 entries).
  const G4int k0_fermi_motion = confMan.Get<G4int>("K0FermiMotion");
  const G4double kHbarc_GeVfm = 0.197327;   // GeV*fm
  // CM-corrected 1p oscillator length b~ = b_bare*sqrt(A/(A-1)):
  // 9Be 1.617*sqrt(9/8)=1.716, 12C 1.669*sqrt(12/11)=1.743 fm.
  const G4double kB9BeTilde_fm = k0_is_c12 ? 1.743 : 1.716;
  const G4double kFermiSigma_GeV =
    kHbarc_GeVfm/(kB9BeTilde_fm*std::sqrt(2.)); // ~0.0813 GeV/c
  // Target/spectator nuclear masses (AME atomic - Z*m_e), consistent
  // with M_P_EFF above: M_target - M_spectator == M_P_EFF.
  const G4double M_9BE_GEV = k0_is_c12 ? 11.174862342 : 8.39275100376;
  const G4double M_8LI_GEV = k0_is_c12 ? 10.252547297 : 7.47136600282;
  // Cap |p| so the off-shell proton's invariant mass^2 stays
  // positive (P(chi5 tail beyond this) ~ 1e-13; physically the HO
  // density there is negligible anyway).
  const G4double kFermiPMax_GeV = 0.7;
  G4double targetPx = 0., targetPy = 0., targetPz = 0.;
  G4double eTargetFermi = M_P_EFF;
  if(k0_fermi_motion != 0){
    G4double pMag;
    do {
      G4double chi2 = 0.;
      for(G4int i = 0; i < 5; ++i){
        const G4double g = G4RandGauss::shoot(0., 1.);
        chi2 += g*g;
      }
      pMag = kFermiSigma_GeV*std::sqrt(chi2);
    } while(pMag > kFermiPMax_GeV);
    const G4double cosT = G4RandFlat::shoot(-1., 1.);
    const G4double sinT = std::sqrt(std::max(0., 1. - cosT*cosT));
    const G4double phi = G4RandFlat::shoot(0., 2.*M_PI);
    targetPx = pMag*sinT*std::cos(phi);
    targetPy = pMag*sinT*std::sin(phi);
    targetPz = pMag*cosT;
    eTargetFermi =
      M_9BE_GEV - std::sqrt(pMag*pMag + M_8LI_GEV*M_8LI_GEV);
  }

  // Measured K-p -> K0bar n angular distribution, Legendre A_l with
  // dsigma/dcosTheta* = sum_l A_l P_l(cosTheta*), cosTheta*=+1 =
  // K0bar forward. Nucl. Phys. B91 (1970) 12 at p_K- = 1.478, 1.519,
  // 1.800 GeV/c (data/legendre_kpk0n.yaml), SELECTED BY BEAM MOMENTUM
  // (linear interpolation, clamped to the measured endpoints).
  //
  // CRUCIAL FIX (2026-07-23, codex E10-2nd handoff sec.17): the shape
  // REVERSES with momentum. The old code used a FIXED Conforto 1.355
  // GeV/c shape (forward/backward ~0.12, strongly BACKWARD-peaked) at
  // every beam momentum -- but at 1.8 GeV/c the real reaction is
  // strongly FORWARD-peaked (fwd/bwd ~6.2). Since forward K0bar ->
  // forward K0_S -> forward pi+ -> into the signal window, the old
  // shape badly UNDER-produced the forward K0bar and thus
  // under-estimated the quasi-free window background at high p_beam.
  // Only the SHAPE matters here; the absolute rate is applied at
  // analysis time via sigma_k0_tot_mb (Bricman).
  static const G4int nA = 10;
  static const G4int nMom = 3;
  static const G4double momGrid[nMom] = {1.478, 1.519, 1.800};
  static const G4double AlGrid[nMom][nA] = {
    {0.168, 0.101, 0.113,-0.116,-0.003,-0.079, 0.110,-0.134,-0.043, 0.036},
    {0.154, 0.169, 0.159,-0.044, 0.072,-0.044, 0.081,-0.143, 0.022,-0.012},
    {0.169, 0.213, 0.282, 0.253, 0.185, 0.094,-0.013,-0.075, 0.031,-0.013}
  };
  G4double Al[nA];
  if(P_BEAM <= momGrid[0]){
    for(G4int l=0;l<nA;++l) Al[l]=AlGrid[0][l];
  } else if(P_BEAM >= momGrid[nMom-1]){
    for(G4int l=0;l<nA;++l) Al[l]=AlGrid[nMom-1][l];
  } else {
    G4int k=0; while(k<nMom-1 && P_BEAM>momGrid[k+1]) ++k;
    const G4double f=(P_BEAM-momGrid[k])/(momGrid[k+1]-momGrid[k]);
    for(G4int l=0;l<nA;++l) Al[l]=(1.-f)*AlGrid[k][l]+f*AlGrid[k+1][l];
  }
  auto dSigmaDCosTheta = [&Al](G4double x){
    G4double s = 0.;
    for(G4int l = 0; l < nA; ++l) s += Al[l]*LegendreP_K0(l, x);
    return std::max(s, 0.);
  };
  static G4double fMax = -1.;
  if(fMax < 0.){
    fMax = 0.;
    for(G4int i = 0; i <= 2000; ++i){
      G4double x = -1. + 2.*i/2000.;
      fMax = std::max(fMax, dSigmaDCosTheta(x));
    }
    fMax *= 1.05; // safety margin for the rejection-sampling envelope
  }

  // c.m. kinematics: K-(beam) + p(bound, at rest unless K0FermiMotion
  // is enabled above) -> K0bar + n.
  const G4double eBeam = std::sqrt(P_BEAM*P_BEAM + M_KM*M_KM);
  TLorentzVector beamLab(0., 0., P_BEAM, eBeam);
  TLorentzVector targetLab(targetPx, targetPy, targetPz, eTargetFermi);
  TLorentzVector cmSystem = beamLab + targetLab;
  const TVector3 cmBoost = cmSystem.BoostVector();
  // Approximation: the Conforto angular distribution below is still
  // sampled relative to the lab z-axis (nominal beam direction), not
  // the exact K-p c.m.-frame beam axis (which tilts slightly once
  // the target proton has transverse Fermi momentum). This is a
  // second-order effect on top of an already-approximate Fermi
  // motion model and is neglected here.
  const G4double sqrtS = cmSystem.M();
  const G4double eStarK0 = (sqrtS*sqrtS + M_K0*M_K0 - M_N*M_N)/(2*sqrtS);
  const G4double pStarK0 =
    std::sqrt(std::max(0., eStarK0*eStarK0 - M_K0*M_K0));

  G4double x, y;
  do {
    x = G4RandFlat::shoot(-1., 1.);
    y = G4RandFlat::shoot(0., fMax);
  } while(y > dSigmaDCosTheta(x));
  const G4double thStarK0 = std::acos(x);
  const G4double phiStarK0 = G4RandFlat::shoot(0., 2.*M_PI);
  TLorentzVector k0Cm(
    pStarK0*std::sin(thStarK0)*std::cos(phiStarK0),
    pStarK0*std::sin(thStarK0)*std::sin(phiStarK0),
    pStarK0*std::cos(thStarK0),
    eStarK0);
  TLorentzVector k0Lab(k0Cm);
  k0Lab.Boost(cmBoost);

  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_size = sizeMan.GetSize("Target")*mm/2;
  const auto vtx = SampleUniformTargetVertex(target_pos, target_size);
  const G4LorentzVector vertex_lv(vtx, 0.);

  m_particleGun->SetParticleDefinition(anti_kaon0);
  m_particleGun->SetParticleMomentumDirection(
    G4ThreeVector(k0Lab.Px(), k0Lab.Py(), k0Lab.Pz()).unit());
  m_particleGun->SetParticleEnergy(
    k0Lab.E()*CLHEP::GeV - anti_kaon0->GetPDGMass());
  m_particleGun->SetParticlePosition(vtx);
  m_particleGun->GeneratePrimaryVertex(anEvent);

  anaMan.SetPrimaryParticle(0, pdg_anti_kaon0, ToG4Lorentz(k0Lab), vertex_lv);
}

//_____________________________________________________________________________
// K- + 9Be(rest) -> K0bar + 9Li(g.s.), the COHERENT charge-exchange
// channel: unlike GenerateK0Production (quasi-free, K- scatters off
// ONE bound, Fermi-moving proton with an 8Li spectator core), here
// K0bar recoils against the WHOLE nucleus in a single bound final
// state. Deliberately kept as a SEPARATE generator/reaction number
// (1002, not folded into 1001's K0FermiMotion machinery): the two
// mechanisms have different targets (single proton vs. whole
// nucleus), different final states (continuum vs. one discrete
// state), different angular distributions (measured Conforto shape
// vs. a steep forward nuclear form factor) and, physically, do not
// interfere (quasi-free populates states above the n+8Li breakup
// threshold; coherent produces a genuine bound state below it) --
// so they are combined as an incoherent sum of independently
// generated, independently normalized samples, exactly as the
// signal and (K-,K0bar) background are already combined at the
// analysis level (beam_scan_absolute_signal.py). Running this as
// its own generator also means its (much rarer) yield gets its own
// dedicated event budget instead of being diluted inside 1001's
// statistics -- see k0-background.md section 5 and
// reaction-kinematics/geant4-scan/coherent_channel.py, which this
// function ports to Geant4 (that script's numbers are the reference
// this generator should reproduce once cross-checked).
//
// dsigma/dOmega(theta_K0) is NOT a measured shape (nothing has ever
// been measured for this exact reaction): it is modelled as a
// Gaussian nuclear form factor |F(q)|^2 = exp(-q_fm^2 R^2/6), R=2.52
// fm (9Be rms radius). The overall normalization (B(GT), spin-flip
// fraction, distortion) only sets the ABSOLUTE rate, which this
// generator does not need to know: like GenerateK0Production, it
// only shapes the KINEMATIC distribution; the absolute weight per
// event is applied at analysis time (matching sigma_k0_tot_mb's role
// for the quasi-free channel).
void
S2SPrimaryGeneratorAction::GenerateK0CoherentProduction(G4Event* anEvent)
{
  m_particleGun = new G4ParticleGun();

  const auto pTable = G4ParticleTable::GetParticleTable();
  const auto anti_kaon0 = pTable->FindParticle("anti_kaon0");
  const auto kaon_minus = pTable->FindParticle("kaon-");

  const G4double M_K0 = anti_kaon0->GetPDGMass()/CLHEP::GeV;
  const G4double M_KM = kaon_minus->GetPDGMass()/CLHEP::GeV;
  const auto pdg_anti_kaon0 = anti_kaon0->GetPDGEncoding();

  // Target and coherent bound recoil (ground state), AME2020 atomic
  // masses - Z*m_e. Selected by conf key "K0TargetA" (9 = 9Be->9Li
  // default, backward compatible; 12 = 12C->12B for the E05 SKS
  // cross-check). Missing key reads as 0 -> treated as 9.
  const G4int k0_target_a = confMan.Get<G4int>("K0TargetA");
  const G4bool k0_is_c12 = (k0_target_a == 12);
  const G4double M_9BE_GEV = k0_is_c12 ? 11.174862342 : 8.39274983976;
  const G4double M_9LI_GEV = k0_is_c12 ? 11.188742741 : 8.40686715782;

  // K- beam momentum: same conf key as the quasi-free channel
  // (K0BeamMom), so the two can be scanned together point by point.
  const G4double k0_beam_mom_conf = confMan.Get<G4double>("K0BeamMom");
  const G4double P_BEAM = (k0_beam_mom_conf > 0.) ? k0_beam_mom_conf : 1.5;

  // Multipole channel selection (2026-07-24, ff_multipole_study.py):
  //   K0CoherentEx [MeV]: daughter excitation energy (default 0 =
  //     ground state; e.g. 2.691 for 9Li* 1/2-).
  //   K0CoherentL: transferred orbital angular momentum, 0/1/2
  //     (default 0 = GT-like spin-flip, the pre-existing behaviour).
  // Missing keys read as 0 via ConfMan's map access, so old conf
  // files keep producing the g.s. GT channel unchanged.
  const G4double coh_ex_mev = confMan.Get<G4double>("K0CoherentEx");
  const G4int coh_ell = confMan.Get<G4int>("K0CoherentL");

  // 1p-shell harmonic-oscillator transition form factors F_L(q),
  // closed forms validated against numerical radial integrals in
  // reaction-kinematics/geant4-scan/ff_multipole_study.py:
  //   x = (q b / hbarc)^2, b = 1.617 fm for 9Be
  //   F_0 = (1 - x/6) exp(-x/4)              (1p->1p, monopole)
  //   F_1 = sqrt(2/45) sqrt(x) (5/2 - x/4) exp(-x/4)  (1p->1d)
  //   F_2 = (x/6) exp(-x/4)                  (1p->1p, quadrupole)
  // These replace the earlier Gaussian-matter exp(-q^2 R^2/6): the
  // multipole structure (not the overall size) decides WHICH bound
  // state dominates at each angle, which is the question here.
  // b = 1.617 fm (9Be) / 1.669 fm (12C), hbar_omega = 45A^-1/3 -
  // 25A^-2/3 (see ff_multipole_study.osc_length_fm).
  const G4double kOscB_fm = k0_is_c12 ? 1.669 : 1.617;
  const G4double kHbarc_GeVfm = 0.197327;

  // Lab-frame two-body momentum of the K0bar at angle theta_lab
  // (rad), for K-(p_beam,+z) + 9Be(rest) -> K0bar(theta_lab) + 9Li --
  // same quadratic-formula solution as
  // hypernuclear-production/scripts/momentum_transfer.two_body_p_out
  // and reaction-kinematics/geant4-scan/coherent_channel.
  // two_body_lab_p (already cross-checked there).
  const G4double mDaughter = M_9LI_GEV + coh_ex_mev*1e-3;
  auto k0MomentumAt = [&](G4double thetaLab) -> G4double {
    const G4double eIn = std::sqrt(P_BEAM*P_BEAM + M_KM*M_KM);
    const G4double eTot = eIn + M_9BE_GEV;
    const G4double s = M_KM*M_KM + M_9BE_GEV*M_9BE_GEV
      + 2.*M_9BE_GEV*eIn;
    const G4double f = s + M_K0*M_K0 - mDaughter*mDaughter;
    const G4double g = 2.*P_BEAM*std::cos(thetaLab);
    const G4double denom = 4.*eTot*eTot - g*g;
    const G4double disc = f*f - M_K0*M_K0*denom;
    return (f*g + 2.*eTot*std::sqrt(std::max(0., disc)))/denom;
  };
  auto formFactorL2 = [&](G4double thetaLab, G4double pK0) -> G4double {
    const G4double q2 = P_BEAM*P_BEAM + pK0*pK0
      - 2.*P_BEAM*pK0*std::cos(thetaLab);
    const G4double qb = std::sqrt(std::max(0., q2))
      / kHbarc_GeVfm * kOscB_fm;
    const G4double x = qb*qb;
    const G4double gauss = std::exp(-x/4.);
    G4double fl = 0.;
    switch(coh_ell){
    case 0: fl = (1. - x/6.)*gauss; break;
    case 1: fl = std::sqrt(2./45.)*qb*(2.5 - x/4.)*gauss; break;
    case 2: fl = (x/6.)*gauss; break;
    default:
      G4Exception("GenerateK0CoherentProduction", "K0COH", // NOLINT
                  FatalException, "K0CoherentL must be 0, 1 or 2");
    }
    return fl*fl;
  };

  // Angular weight for the coherent (GT-driven, spin-flip) K0bar:
  //   dsigma/dOmega ~ |t_spin-flip|^2 * |F(q)|^2,  with the meson
  //   (spin-0) SPIN-FLIP amplitude vanishing at theta*=0 (~ sin^2
  //   theta* near forward -- see k0-background.md 5.3 point 2). The
  //   old version used |F(q)|^2 alone (peaked at theta=0), which
  //   OVER-produced the deepest-wall K0bar exactly where it lands in
  //   the signal window. Multiplying by sin^2(theta) suppresses
  //   theta=0 and pushes the K0bar to a small non-zero angle, so the
  //   coherent background right at the wall is reduced (favourable).
  //   theta here is the lab K0bar angle; for coherent recoil off a
  //   heavy nucleus theta_lab ~ theta*_cm at these small angles, so
  //   sin^2(theta_lab) is used as the forward-suppression factor (the
  //   exact spin-flip angular shape needs the elementary PWA -- codex
  //   E10-2nd handoff sec.7.2/17, documented as an uncertainty).
  //   Added 2026-07-23. Full theta-sampling weight = |F|^2 * sin^2
  //   (spin-flip) * sin (solid-angle Jacobian) = |F|^2 * sin^3(theta).
  //   2026-07-24: the explicit sin^2 applies to the L=0 GT channel
  //   only (F_0(0)=1 needs it); F_1 ~ q and F_2 ~ q^2 already vanish
  //   at theta=0 by themselves, so L=1/2 use just |F_L|^2 * sin.
  auto angWeight = [&](G4double th) -> G4double {
    const G4double s = std::sin(th);
    const G4double sf = (coh_ell == 0) ? s*s : 1.;
    return formFactorL2(th, k0MomentumAt(th)) * sf * s;
  };
  static const G4double thetaMaxRad = 20.*CLHEP::deg/CLHEP::rad;
  static G4double envelope = -1.;
  if(envelope < 0.){
    envelope = 0.;
    for(G4int i = 0; i <= 400; ++i){
      const G4double th = thetaMaxRad*i/400.;
      envelope = std::max(envelope, angWeight(th));
    }
    envelope *= 1.05; // safety margin for the rejection-sampling envelope
  }
  G4double trialTheta, trialY;
  do {
    trialTheta = G4RandFlat::shoot(0., thetaMaxRad);
    trialY = G4RandFlat::shoot(0., envelope);
  } while(trialY > angWeight(trialTheta));
  const G4double thStarK0 = trialTheta;
  const G4double phiK0 = G4RandFlat::shoot(0., 2.*M_PI);
  const G4double pK0 = k0MomentumAt(thStarK0);
  const G4double eK0 = std::sqrt(pK0*pK0 + M_K0*M_K0);

  TLorentzVector k0Lab(
    pK0*std::sin(thStarK0)*std::cos(phiK0),
    pK0*std::sin(thStarK0)*std::sin(phiK0),
    pK0*std::cos(thStarK0),
    eK0);

  static const auto& target_pos = geomMan.GetGlobalPosition("Target")*mm;
  static const auto& target_size = sizeMan.GetSize("Target")*mm/2;
  const auto vtx = SampleUniformTargetVertex(target_pos, target_size);
  const G4LorentzVector vertex_lv(vtx, 0.);

  m_particleGun->SetParticleDefinition(anti_kaon0);
  m_particleGun->SetParticleMomentumDirection(
    G4ThreeVector(k0Lab.Px(), k0Lab.Py(), k0Lab.Pz()).unit());
  m_particleGun->SetParticleEnergy(
    k0Lab.E()*CLHEP::GeV - anti_kaon0->GetPDGMass());
  m_particleGun->SetParticlePosition(vtx);
  m_particleGun->GeneratePrimaryVertex(anEvent);

  anaMan.SetPrimaryParticle(0, pdg_anti_kaon0, ToG4Lorentz(k0Lab), vertex_lv);
}
