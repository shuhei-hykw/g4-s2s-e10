// -*- C++ -*-

#include "S2SAnaManager.hh"

#include <bitset>
#include <algorithm>

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Step.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
//#include "G4Poisson.h"

#include "Randomize.hh"

#include <TMath.h>
#include <TRandom3.h>
#include <TSystem.h>

#include "RootHelper.hh"
#include "GeneratorParticleBranches.hh"
#include "TPCMlFeature.hh"

#include <iomanip>
#include <cmath>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "DCHit.hh"
#include "TOFHit.hh"
#include "HTOFHit.hh"
#include "TPCHit.hh"
#include "ACHit.hh"
#include "WCHit.hh"
#include "VPHit.hh"
#include "HistMan.hh"

#include "GeHit.hh"
#include "RCHit.hh"
#include "PDHit.hh"


namespace
{
using CLHEP::mm;
using CLHEP::MeV;
using CLHEP::ns;
const auto& confMan = ConfMan::GetInstance();
const auto& histMan = HistMan::GetInstance();
const auto qnan = TMath::QuietNaN();
Event event;
std::map<TString, TH1*> hmap;
std::vector<G4int> n_acc(kTriggerFlagSize, 0);
std::size_t kMlTrackCount = 0;
const bool noHist = (confMan.Get<G4String>("BranchStyle") == "E90ML");

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

const std::vector<G4String>&
e63_branch_names()
{
  static const std::vector<G4String> names{
    "RC-X", "RC+X",
    "RC-X-PDY", "RC-X-PDZ",
    "RC+X-PDY", "RC+X-PDZ",
    "Ge", "BGO"
  };
  return names;
}

const std::vector<std::pair<G4String, ERCTriggerFlag>>&
e63_rc_collection_names()
{
  static const std::vector<std::pair<G4String, ERCTriggerFlag>> names{
    {"RC-X", kRCMinusX},
    {"RC+X", kRCPlusX}
  };
  return names;
}

const std::vector<std::pair<G4String, ERCTriggerFlag>>&
e63_pd_collection_names()
{
  static const std::vector<std::pair<G4String, ERCTriggerFlag>> names{
    {"RC-X-PDY", kRCMinusXPDY},
    {"RC-X-PDZ", kRCMinusXPDZ},
    {"RC+X-PDY", kRCPlusXPDY},
    {"RC+X-PDZ", kRCPlusXPDZ}
  };
  return names;
}

}

//_____________________________________________________________________________
S2SAnaManager::S2SAnaManager()
  : m_file_name("tmp.root"),
    fActive_(true),
    fTriggered(false),
    m_file(),
    m_tree()
{
}

//_____________________________________________________________________________
S2SAnaManager::~S2SAnaManager()
{
}

//_____________________________________________________________________________
void
S2SAnaManager::BeginOfRun( const G4Run* /* aRun */)
{
  fActive_=true;
  kMlTrackCount = static_cast<std::size_t>(confMan.Get<G4int>("TPCMt") + 1);
  m_file = new TFile(m_file_name, "recreate");
  //static auto obj = new TNamed("conf", confMan.ConfPath()+confMan.ConfBuf());
  static auto obj = new TNamed("conf", confMan.ConfBuf()); // [seong]
  obj->Write();
  static auto git = new TNamed
    ("git", ("\n"+gSystem->GetFromPipe("git log -1")).Data());
  git->Write();
  m_tree = new TTree("g4s2s", "S-2S simulation");
  event.hits.clear();
  event.evnum = -1;
  event.trig.assign(kTriggerFlagSize, false);
  event.rctrig.assign(kRCTriggerFlagSize, false); //[seong]
  ResetMlFeatures(event, confMan, qnan, kMlTrackCount);
  DefineTree();
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  for(const auto& sd_name : std::vector<G4String>{
      "PRM", "SDC1","SDC2","SDC3","SDC4","SDC5", "TOF", "AC1", "WC", "VP"}
        // S2SDetectorConstruction::GetSDList()
    ){
    G4cout << "   make branch : " << sd_name << G4endl;
    MakeBranch(sd_name);
    if (!noHist) MakeHistogram(sd_name);
  }
  if(experiment == 63){
    std::vector<G4String> e63_branches = e63_branch_names();
    for(const auto& sd_name : e63_branches){
      G4cout << "   make branch : " << sd_name << G4endl;
      MakeBranch(sd_name);
      if (!noHist && sd_name != "BGO") MakeHistogram(sd_name);
    }
  }
  // E10-2nd shares the E90/HypTPC setup, so record TPC/HTOF/SAC hits
  // when HypTPCSetup:1 (or E90) -- these are the K0S->pi+pi-
  // displaced-vertex / pi- rejection handles.
  static const G4bool hyptpc_setup = (confMan.Get<G4int>("HypTPCSetup") != 0);
  if(experiment == 90 || hyptpc_setup){
    for(const auto& sd_name : std::vector<G4String>{"TPC", "HTOF", "SAC"})
    {
      G4cout << "   make branch : " << sd_name << G4endl;
      MakeBranch(sd_name);
      if (!noHist) MakeHistogram(sd_name);
    }
  }
  {
    const auto generator = confMan.Get<G4int>("Generator");
    for(const auto& branch : GeneratorParticleBranches::BranchList(generator)){
      G4cout << "   make branch : " << branch << G4endl;
      MakeBranch(branch);
    }
  }
  for(auto& h: hmap){
    h.second->Reset();
  }
  InitializeEvent();
  n_acc.clear();
  n_acc.resize(kTriggerFlagSize);
}

//_____________________________________________________________________________
void
S2SAnaManager::EndOfRun(const G4Run* /* aRun */)
{
  m_file->cd();
  if(confMan.Get<G4bool>("TREE"))
    m_tree->Write();
  if (confMan.Get<G4String>("BranchStyle") != "E90ML") {
    for (auto& h : hmap) {
      h.second->Write();
    }
  }
  m_file->Close();
}

//_____________________________________________________________________________
void
S2SAnaManager::BeginOfPrimaryAction()
{
  event.x0In = qnan;
  event.y0In = qnan;
  event.u0In = qnan;
  event.v0In = qnan;
  event.uDeg = qnan;
  event.vDeg = qnan;
  event.p0   = qnan;
  event.pB   = qnan;
  event.t0   = qnan;
  event.phi0 = qnan;
  event.theta0 = qnan;
  // event.Id = qnan;
  //  G4cout<<"BeginOfPrimaryAction"<<G4endl;

#if 0
  // E63 for weak pi-
  event.x1 = qnan;
  event.y1 = qnan;
  event.z1 = qnan;
  event.px1 = qnan;
  event.py1 = qnan;
  event.pz1 = qnan;
  event.p1  = qnan;
  event.t1  = qnan;
#endif
}

//_____________________________________________________________________________
void
S2SAnaManager::MakeBranch(const G4String& sd_name)
{
  static const Int_t bufsize = 32000;
  // Always ensure the container exists for histogram filling.
  (void)event.hits[sd_name];
  if (confMan.Get<G4String>("BranchStyle") == "E90ML")
    return;
  m_tree->Branch(sd_name.data(),
                 "std::vector<TParticle>",
                 &event.hits[sd_name], bufsize, -1);
}

//_____________________________________________________________________________
void
S2SAnaManager::MakeHistogram(const G4String& sd_name)
{
  if (noHist) return;
  if(sd_name == "PRM"){
    const auto& params = histMan.Get("PTheta");
    TString key = sd_name + "PThetaGen";
    TString title = sd_name + " P%Theta (Generate); [deg]; [GeV/c]";
    hmap[key] = new TH2D(key, title,
                         params.at(0), params.at(1), params.at(2),
                         params.at(3), params.at(4), params.at(5));
    key = sd_name + "PThetaAcc";
    title = sd_name + " P%Theta (Accept); [deg]; [GeV/c]";
    hmap[key] = new TH2D(key, title,
                         params.at(0), params.at(1), params.at(2),
                         params.at(3), params.at(4), params.at(5));
    for(G4int i=0, n=TriggerFlag.size(); i<n; ++i){
      key = sd_name + "PThetaAcc" + TriggerFlag.at(i);
      title = sd_name + " P%Theta (Accept at " + TriggerFlag.at(i) + ");"
        + " [deg]; [GeV/c]";
      hmap[key] = new TH2D(key, title,
                           params.at(0), params.at(1), params.at(2),
                           params.at(3), params.at(4), params.at(5));
    }
  }else{
    for(const auto& suffix: std::vector<G4String>
          { "Nhits", "HitPat", "X", "Y", "Z", "U", "V",
            "Y%X", "V%U", "U%X", "V%Y" }){
      TString key = sd_name + suffix;
      const auto& params = histMan.Get(key);
      TString title = sd_name + " " + suffix;
      key.ReplaceAll("%", "");
      if(G4StrUtil::contains(suffix, "%")){
        hmap[key] = new TH2D(key, title,
                             params.at(0), params.at(1), params.at(2),
                             params.at(3), params.at(4), params.at(5));
      }else{
        hmap[key] = new TH1D(key, title,
                             params.at(0), params.at(1), params.at(2));
      }
    }
  }
}

//_____________________________________________________________________________
void
S2SAnaManager::SetPrimaryData(double x0, double y0, double z0,
                              double u0, double v0, double phi, double theta,
                              double p0,double pB, int /* ParIdNb */)
{
  event.x0In = x0; // generated position (x)
  event.y0In = y0; // generated position (y)
  event.z0In = z0; // generated position (z)
  //event.u0In = u0; // x' in rad
  event.u0In = -1.0 * u0; // x' in rad
  event.v0In = v0; // y' in rad
  //event.uDeg = atan(u0)*TMath::RadToDeg(); // x' in deg
  //event.vDeg = atan(v0)*TMath::RadToDeg(); // y' in deg
  //event.uDeg = -1.0 * u0 * TMath::RadToDeg(); // x' in deg
  //event.vDeg = v0*TMath::RadToDeg(); // y' in deg
  event.uDeg = u0;
  event.vDeg = v0;
  event.phi0 = phi; // Phi in rad
  event.theta0 = theta; // Theta in rad
  event.p0   = p0; // Momentum
  event.pB   = pB; // Kinetic energy
  //  event.Id = ParIdNb;
  //  G4cout<<"setPrimaryData"<<G4endl;
}

// E63 for weak pi
//_____________________________________________________________________________
void
S2SAnaManager::SetSecondaryData(double x1, double y1, double z1,
				double px1, double py1, double pz1,
				double p1,double t1)
{
  event.x1 = x1; // generated position (x)
  event.y1 = y1; // generated position (y)
  event.z1 = z1; // generated position (z)
  event.px1 = px1;
  event.py1 = py1;
  event.pz1 = pz1;
  event.p1   = p1; // Momentum
  event.t1   = t1; // Kinetic energy
  //  event.Id = ParIdNb;
  //  G4cout<<"setPrimaryData"<<G4endl;
}

//_____________________________________________________________________________
void S2SAnaManager::SetProcessData(G4int nP, G4int nN, G4int nL,
				 G4int nSm, G4int nSz, G4int nSp,
				 G4int nXm, G4int nXz, G4int nXsm,
				 G4int nXsz,G4int nPim,G4int nPiz,
				 G4int nPip,G4int nKm,G4int nKp)
{
  event.nP = nP;
  event.nN = nN;
  event.nL = nL;
  event.nSm = nSm;
  event.nSz = nSz;
  event.nSp = nSp;
  event.nXm = nXm;
  event.nXz = nXz;
  event.nXsm = nXsm;
  event.nXsz = nXsz;
  event.nPim = nPim;
  event.nPiz = nPiz;
  event.nPip = nPip;
  event.nKm = nKm;
  event.nKp = nKp;
}

void S2SAnaManager::BeginOfEvent(const G4Event *anEvent)
{
  event.evnum = anEvent->GetEventID();
}

void S2SAnaManager::EndOfEvent(const G4Event *anEvent)
{
  auto HCE = anEvent->GetHCofThisEvent();
  auto SDMan = G4SDManager::GetSDMpointer();
  std::bitset<kTriggerFlagSize> trigger_flag;
  // Common e10 acceptance definition: is there a SINGLE track that
  // (a) leaves a charged hit in TOF, a Cherenkov-threshold hit in AC1
  // (n=1.05, enforced by ACSD::ProcessHits itself) and a charged hit
  // in WC, all within a 100 ns trigger coincidence window, AND (b) is
  // also consistently present in the tracking chambers SDC2-SDC5 (a
  // real analysis only reconstructs a momentum for a track threading
  // all four planes)? This mirrors the real hardware trigger plus a
  // minimal track-quality requirement, regardless of species or
  // whether the track is a primary or a decay secondary -- unlike
  // the historical trigger_flag[kTOF/kAC1/kWC] above, which requires
  // hit->IsPrimary() and a GenPID species match.
  //
  // SDC1 is intentionally excluded: for a generator whose primary
  // decays in flight (e.g. K0 background), the decay vertex is
  // typically downstream of SDC1, so requiring it would reject every
  // such event by construction (verified empirically: zero overlap
  // for any K0-decay secondary in this geometry -- see
  // analysis-note.md 2026-07-13 entries).
  //
  // Scoped to Experiment==10 so other experiments' analyses
  // (63/90/...) are unaffected. Originally validated offline in
  // ana/k0_missing_mass_e10_kpi.C on a TREE:1 sample; moved here so
  // routine TREE:0 productions (e.g. a beam-momentum scan) get the
  // corrected PRMPThetaGen/Acc ratio directly, without needing large
  // per-event tree output.
  const G4double kE10CoincidenceWindow_ns = 100.;
  std::map<G4int, G4double> e10_tofTime, e10_ac1Time, e10_wcTime;
  std::set<G4int> e10_sdc2Tracks, e10_sdc3Tracks,
    e10_sdc4Tracks, e10_sdc5Tracks;
  std::bitset<kRCTriggerFlagSize> rc_trigger_flag;
  static const G4int experiment = confMan.Get<G4int>("Experiment");
  static const G4int requiredTPCMt = confMan.Get<G4int>("TPCMt");
  static const G4int minTPCPadHits = std::max<G4int>(1, confMan.Get<G4int>("MinHit"));
  static const G4double truncateRate = confMan.Get<G4double>("Trunc");
  std::vector<MlTrackFeature> mlTrackFeatures;
  // GenPID (optional conf key, see GenPidToName), if set, takes
  // precedence over the historical per-experiment species below.
  static const G4int gen_pid = confMan.Get<G4int>("GenPID");
  static const G4String gen_pid_name = GenPidToName(gen_pid);
  G4String particle_name = gen_pid_name;
  if(particle_name.empty()){
    particle_name = "kaon+";
    if(experiment == 63) particle_name = "pi-"; //for E63
    //particle_name = "kaon-"; //for E63
    if(experiment == 90) particle_name = "pi-";
  }

  {
    //G4String name = "SDC"+std::to_string(k);
    G4String name = "SDC1";
    static const auto id = SDMan->GetCollectionID(name);
    if(id >= 0){
      //auto HC = dynamic_cast<SDCHitsCollection*>(HCE->GetHC(id));
      auto HC = dynamic_cast<DCHitsCollection*>(HCE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        SetHitData((*HC)[i]);
      }
      SetNhits(name, HC->entries());
    }
  }
  {
    G4String name = "SDC2";
    static const auto id = SDMan->GetCollectionID(name);
    if(id >= 0){
      //auto HC = dynamic_cast<SDCHitsCollection*>(HCE->GetHC(id));
      auto HC = dynamic_cast<DCHitsCollection*>(HCE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        SetHitData((*HC)[i]);
        e10_sdc2Tracks.insert((*HC)[i]->GetTrackID());
      }
      SetNhits(name, HC->entries());
    }
  }
  {
    G4String name = "SDC3";
    static const auto id = SDMan->GetCollectionID(name);
    if(id >= 0){
      //auto HC = dynamic_cast<SDCHitsCollection*>(HCE->GetHC(id));
      auto HC = dynamic_cast<DCHitsCollection*>(HCE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        SetHitData((*HC)[i]);
        e10_sdc3Tracks.insert((*HC)[i]->GetTrackID());
      }
      SetNhits(name, HC->entries());
    }
  }
  {
    G4String name = "SDC4";
    static const auto id = SDMan->GetCollectionID(name);
    if(id >= 0){
      //auto HC = dynamic_cast<SDCHitsCollection*>(HCE->GetHC(id));
      auto HC = dynamic_cast<DCHitsCollection*>(HCE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        SetHitData((*HC)[i]);
        e10_sdc4Tracks.insert((*HC)[i]->GetTrackID());
      }
      SetNhits(name, HC->entries());
    }
  }
  {
    G4String name = "SDC5";
    static const auto id = SDMan->GetCollectionID(name);
    if(id >= 0){
      //auto HC = dynamic_cast<SDCHitsCollection*>(HCE->GetHC(id));
      auto HC = dynamic_cast<DCHitsCollection*>(HCE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        SetHitData((*HC)[i]);
        e10_sdc5Tracks.insert((*HC)[i]->GetTrackID());
      }
      SetNhits(name, HC->entries());
    }
  }
  /*static const auto id = SDMan->GetCollectionID("SDC");
  if(id >= 0){
    //auto HC = dynamic_cast<SDCHitsCollection*>(HCE->GetHC(id));
    auto HC = dynamic_cast<DCHitsCollection*>(HCE->GetHC(id));
    for(G4int i=0, n=HC->entries(); i<n; ++i){
      SetHitData((*HC)[i]);
    }
    SetNhits("SDC", HC->entries());
   }*/
  {
    static const auto id = SDMan->GetCollectionID("TOF");
    if(id >= 0){
      auto HC = dynamic_cast<TOFHitsCollection*>(HCE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        auto hit = (*HC)[i];
        if(hit->Is(particle_name) && hit->IsPrimary()) trigger_flag[kTOF] = true;
        if(experiment == 90) trigger_flag[kE90TOF] = true;
        if(hit->GetCharge() != 0.) e10_tofTime[hit->GetTrackID()] = hit->GetTime();
        SetHitData(hit);
      }
      SetNhits("TOF", HC->entries());
    }
  }
  {
    static const auto id = SDMan->GetCollectionID("AC1");
    if(id >= 0){
      auto HC = dynamic_cast<ACHitsCollection*>(HCE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        auto hit = (*HC)[i];
        if(hit->Is(particle_name) && hit->IsPrimary()) trigger_flag[kAC1] = true;
        if(experiment == 90 && hit->Is("pi-")) trigger_flag[kE90AC1] = true;
        // ACSD only inserts a hit once the track's velocity exceeds
        // the aerogel's Cherenkov threshold (see ACSD::ProcessHits),
        // so any hit here already satisfies that condition.
        if(hit->GetCharge() != 0.) e10_ac1Time[hit->GetTrackID()] = hit->GetTime();
        SetHitData(hit);
      }
      SetNhits("AC1", HC->entries());
    }
  }
  {
    static const auto id = SDMan->GetCollectionID("WC");
    if(id >= 0){
      auto HC = dynamic_cast<WCHitsCollection*>(HCE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        auto hit = (*HC)[i];
        if(hit->Is(particle_name) && hit->IsPrimary()) trigger_flag[kWC] = true;
        // WCHit declares its own (unused, never-Set) time_/GetTime()
        // left over from dead legacy code (see the #if 0 block in
        // WCSD::ProcessHits), which HIDES VHitInfo::GetTime() by
        // name lookup -- must qualify explicitly to get the real
        // global hit time (WCHit::GetTime() returns uninitialized
        // garbage). TOFHit/ACHit/DCHit have no such shadowing member.
        if(hit->GetCharge() != 0.){
          e10_wcTime[hit->GetTrackID()] = hit->VHitInfo::GetTime();
        }
        SetHitData(hit);
      }
      SetNhits("WC", HC->entries());
    }
  }
  {
    static const auto id = SDMan->GetCollectionID("VP");
    if(id >= 0){
      auto HC = dynamic_cast<VPHitsCollection*>(HCE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
        auto hit = (*HC)[i];
        if(hit->Is(particle_name) && hit->IsPrimary()){
          trigger_flag[kVP1-1+hit->GetCopyNumber()] = true;
        }
        SetHitData(hit);
      }
      SetNhits("VP", HC->entries());
    }
  }
  if(experiment == 63)
  {
    {
      G4String name = "Ge";
      static const auto id = SDMan->GetCollectionID(name);
      if(id >= 0){
	auto HC = dynamic_cast<GeHitsCollection*>(HCE->GetHC(id));
	for(G4int i=0, n=HC->entries(); i<n; ++i){
	  SetHitData((*HC)[i]);
	}
	SetNhits(name, HC->entries());
      }
    }
    {
      G4String name = "BGO";
      static const auto id = SDMan->GetCollectionID(name);
      if(id >= 0){
	auto HC = dynamic_cast<GeHitsCollection*>(HCE->GetHC(id));
	// BGO is branch-only for veto diagnostics; no histogram template is defined.
	if(HC){
	  for(G4int i=0, n=HC->entries(); i<n; ++i){
	    auto hit = (*HC)[i];
	    if(hit && hit->GetParticle())
	      event.hits.at(name).push_back(*hit->GetParticle());
	  }
	}
      }
    }
    for(const auto& spec : e63_rc_collection_names()){
      const auto& name = spec.first;
      const auto flag = spec.second;
      const auto id = SDMan->GetCollectionID(name);
      if(id >= 0){
	auto HC = dynamic_cast<RCHitsCollection*>(HCE->GetHC(id));
	for(G4int i=0, n=HC->entries(); i<n; ++i){
	  auto hit = (*HC)[i];
	  if(hit->Is(particle_name) && hit->IsWeakPi()) rc_trigger_flag[flag] = true;
	  SetHitData((*HC)[i]);
	}
	SetNhits(name, HC->entries());
      }
    }
    for(const auto& spec : e63_pd_collection_names()){
      const auto& name = spec.first;
      const auto flag = spec.second;
      const auto id = SDMan->GetCollectionID(name);
      if(id >= 0){
	auto HC = dynamic_cast<PDHitsCollection*>(HCE->GetHC(id));
	for(G4int i=0, n=HC->entries(); i<n; ++i){
	  auto hit = (*HC)[i];
	  if(hit->Is(particle_name) && hit->IsWeakPi()) rc_trigger_flag[flag] = true;
	  SetHitData((*HC)[i]);
	}
	SetNhits(name, HC->entries());
      }
    }
  }

  if (experiment == 90)
  {
    {
      static const auto id = SDMan->GetCollectionID("HTOF");
      if(id >= 0){
        auto HC = dynamic_cast<HTOFHitsCollection*>(HCE->GetHC(id));
        for(G4int i=0, n=HC->entries(); i<n; ++i){
          auto hit = (*HC)[i];
          SetHitData(hit);
        }
        SetNhits("HTOF", HC->entries());
      }
    }
    {
      static const auto id = SDMan->GetCollectionID("TPC");
      if(id >= 0){
        auto HC = dynamic_cast<TPCHitsCollection*>(HCE->GetHC(id));
        if(HC){
          std::unordered_map<G4int, G4int> trackHitCounts;
          std::unordered_map<G4int, G4int> trackPdg;
          std::unordered_set<G4int> primaryPiMinusTracks;
          std::unordered_map<G4int, std::vector<TpcMlHit>> tpcHitMap;
          for(G4int i=0, n=HC->entries(); i<n; ++i){
            auto hit = (*HC)[i];
            SetHitData(hit);
            const auto trackId = hit->GetTrackID();
            const auto pdg = hit->GetPDGEncoding();
            trackHitCounts[trackId]++;       // count TPC hits per track
            trackPdg.emplace(trackId, pdg);   // store PDG
            if(pdg == -211 && hit->IsPrimary()){
              primaryPiMinusTracks.insert(trackId); // reject primary beam pi-
            }
            TpcMlHit mlhit;
            mlhit.trackId = trackId;
            mlhit.pdg     = pdg;
            const auto* ptcl = hit->GetParticle();
            mlhit.parentId = ptcl ? ptcl->GetMother(0) : -1;
            mlhit.time    = hit->GetTime()/ns;
            mlhit.edep    = hit->GetEnergyDeposit();
            mlhit.pos     = hit->GetPosition();
            mlhit.mom     = hit->GetMomentum();
            tpcHitMap[trackId].push_back(std::move(mlhit));
          }

          std::unordered_set<G4int> multiplicityTracks;
          for(const auto& kv : trackHitCounts){
            const auto trackId = kv.first;
            const auto nhit = kv.second;
            const auto pdg = trackPdg.at(trackId);
            const bool isAcceptedParticle = (pdg == 2212 || pdg == 211 || pdg == -211);
            const bool isPrimaryPiMinus = primaryPiMinusTracks.count(trackId) > 0;
            if(isAcceptedParticle && !isPrimaryPiMinus && nhit >= minTPCPadHits){
              multiplicityTracks.insert(trackId);
            }
          }
          event.TPCMt = static_cast<G4int>(multiplicityTracks.size());

          for (auto& kv : tpcHitMap) {
            if (kv.second.size() < static_cast<std::size_t>(minTPCPadHits))
              continue;
            mlTrackFeatures.push_back(CalculateMlFeature(kv.second, truncateRate));
          }

          SetNhits("TPC", HC->entries());
        }
      }
    }
    {
      static const auto id = SDMan->GetCollectionID("SAC");
      if(id >= 0){
        auto HC = dynamic_cast<ACHitsCollection*>(HCE->GetHC(id));
        for(G4int i=0, n=HC->entries(); i<n; ++i){
          auto hit = (*HC)[i];
          if(hit->Is("pi-")) trigger_flag[kE90SAC] = true;
          SetHitData(hit);
        }
        SetNhits("SAC", HC->entries());
      }
    }
  }

  {
    // G4cout << "Acc eff." << G4endl;
    auto particle = event.hits.at("PRM").at(0);
    if (!noHist) {
      for(G4int i=0, n=TriggerFlag.size(); i<n; ++i){
        if(trigger_flag[i]){
          n_acc[i]++;
          hmap.at("PRMPThetaAcc"+TriggerFlag.at(i))->
            Fill(particle.Theta()/CLHEP::deg, particle.P()/CLHEP::GeV);
        }
      }

      G4bool e10_trackMatched = false;
      if(experiment == 10){
        for(const auto& tofPair : e10_tofTime){
          const G4int trackId = tofPair.first;
          const auto itAc1 = e10_ac1Time.find(trackId);
          const auto itWc  = e10_wcTime.find(trackId);
          if(itAc1 == e10_ac1Time.end()) continue;
          if(itWc  == e10_wcTime.end())  continue;
          if(!e10_sdc2Tracks.count(trackId)) continue;
          if(!e10_sdc3Tracks.count(trackId)) continue;
          if(!e10_sdc4Tracks.count(trackId)) continue;
          if(!e10_sdc5Tracks.count(trackId)) continue;
          const G4double tMin = std::min({tofPair.second, itAc1->second,
                                           itWc->second});
          const G4double tMax = std::max({tofPair.second, itAc1->second,
                                           itWc->second});
          if(tMax - tMin > kE10CoincidenceWindow_ns) continue;
          e10_trackMatched = true;
          break;
        }
      }
      const bool isAccepted = (experiment == 10)
        ? e10_trackMatched
        : (true
           && trigger_flag[kVP1]
           && trigger_flag[kVP2]
           && trigger_flag[kVP3]
           && trigger_flag[kVP4]
           && trigger_flag[kVP5]
           && trigger_flag[kVP6]
           && trigger_flag[kVP7]
           && trigger_flag[kVP8]
           && trigger_flag[kVP9]
           && trigger_flag[kVP10]
           && trigger_flag[kTOF]
           && trigger_flag[kWC]);
      if(isAccepted){
        hmap.at("PRMPThetaAcc")->
          Fill(particle.Theta()/CLHEP::deg, particle.P()/CLHEP::GeV);
      }
    }
  }

  {
    for(G4int i=0; i<kTriggerFlagSize; ++i){
      event.trig[i] = trigger_flag[i];
    }
    for(G4int i=0; i<kRCTriggerFlagSize; ++i){ //[seong]
      event.rctrig[i] = rc_trigger_flag[i];
    }
  }

  bool storeEvent = true;
  if (experiment == 90 && confMan.Get<G4String>("BranchStyle") == "E90ML") {
    const bool isPiTrigger = trigger_flag[kE90TOF] && trigger_flag[kE90SAC];
    const bool hasExpectedTracks = (mlTrackFeatures.size() == kMlTrackCount);
    if (isPiTrigger && hasExpectedTracks && event.TPCMt == requiredTPCMt) {
      std::sort(mlTrackFeatures.begin(), mlTrackFeatures.end(),
                [](const MlTrackFeature& a, const MlTrackFeature& b)
                { return a.dedx > b.dedx; });
      for (std::size_t i = 0; i < kMlTrackCount; ++i) {
        event.mlUx[i]    = static_cast<float>(mlTrackFeatures[i].ux);
        event.mlUy[i]    = static_cast<float>(mlTrackFeatures[i].uy);
        event.mlUz[i]    = static_cast<float>(mlTrackFeatures[i].uz);
        event.mlDedx[i]  = static_cast<float>(mlTrackFeatures[i].dedx);
        // event.mlPid[i]   = mlTrackFeatures[i].pidCode;
      }
      event.mlMM = static_cast<float>(
        CalculateMissingMass(event.pB, event.p0, event.theta0, event.phi0));
      storeEvent = true;
    } else {
      storeEvent = false;
    }
  }

  if(confMan.Get<G4bool>("TREE") && storeEvent)
    m_tree->Fill();

  InitializeEvent();
}

//_____________________________________________________________________________
void
S2SAnaManager::SetNhits(const G4String& sd_name, G4int nhits)
{
  if (noHist) return;
  hmap[sd_name + "Nhits"]->Fill(nhits);
}

//_____________________________________________________________________________
void
S2SAnaManager::SetHitData(const VHitInfo* hit)
{
  if(hit && hit->GetParticle()){
    const auto& name = hit->GetDetectorName();
    const auto& p = hit->GetParticle();
    event.hits.at(name).push_back(*p);
    if (!noHist) {
      hmap[name + "HitPat"]->Fill(p->GetMother(1));
      hmap[name + "X"]->Fill(p->Vx());
      hmap[name + "Y"]->Fill(p->Vy());
      hmap[name + "Z"]->Fill(p->Vz());
      hmap[name + "U"]->Fill(p->Px()/p->Pz());
      hmap[name + "V"]->Fill(p->Py()/p->Pz());
      hmap[name + "YX"]->Fill(p->Vx(), p->Vy());
      hmap[name + "VU"]->Fill(p->Px()/p->Pz(), p->Py()/p->Pz());
      hmap[name + "UX"]->Fill(p->Vx(), p->Px()/p->Pz());
      hmap[name + "VY"]->Fill(p->Vy(), p->Py()/p->Pz());
    }
  }
}

void S2SAnaManager::InitializeEvent()
{
  for(auto& pair: event.hits){
    pair.second.clear();
  }
  ResetMlFeatures(event, confMan, qnan, kMlTrackCount);
  event.TPCMt = 0;
}

void S2SAnaManager::DefineTree()
{
  if (confMan.Get<G4String>("BranchStyle") == "E90ML") {
    const auto mlTrackCount = kMlTrackCount; // Mt + (scat pi-)
    m_tree->Branch("label", &event.label, "label/I");
    m_tree->Branch("mm", &event.mlMM, "mm/F");
    for (std::size_t i = 0; i < mlTrackCount; ++i) {
      m_tree->Branch(Form("t%zu_ux", i),   &event.mlUx[i],   Form("t%zu_ux/F", i));
      m_tree->Branch(Form("t%zu_uy", i),   &event.mlUy[i],   Form("t%zu_uy/F", i));
      m_tree->Branch(Form("t%zu_uz", i),   &event.mlUz[i],   Form("t%zu_uz/F", i));
      m_tree->Branch(Form("t%zu_dedx", i), &event.mlDedx[i], Form("t%zu_dedx/F", i));
      // m_tree->Branch(Form("t%zu_pdg", i),  &event.mlPdg[i],  Form("t%zu_pdg/I", i));
    }
    return;
  }

  m_tree->Branch("evnum", &event.evnum, "evnum/I");
  m_tree->Branch("trig", &event.trig);
  m_tree->Branch("x0",&event.x0In, "x0/D");
  m_tree->Branch("y0",&event.y0In, "y0/D");
  m_tree->Branch("z0",&event.z0In, "z0/D");
  //   m_tree->Branch("u0In",&event.u0In, "u0In/D");
  //   m_tree->Branch("v0In",&event.v0In, "v0In/D");
  m_tree->Branch("xp0",&event.uDeg, "xp0/D");
  m_tree->Branch("yp0",&event.vDeg, "yp0/D");
  m_tree->Branch("phi0",&event.phi0, "phi0/D");
  m_tree->Branch("theta0",&event.theta0, "theta0/D");
  m_tree->Branch("p0",&event.p0,   "p0/D");
  m_tree->Branch("pB",&event.pB,   "pB/D");

  // E63 for weak pi
#if 0
  m_tree->Branch("rctrig", &event.rctrig);
  m_tree->Branch("x1",  &event.x1,  "x1/D"); // [mm]
  m_tree->Branch("y1",  &event.y1,  "y1/D");
  m_tree->Branch("z1",  &event.z1,  "z1/D");
  m_tree->Branch("px1", &event.px1, "px1/D"); // [MeV/c]
  m_tree->Branch("py1", &event.py1, "py1/D");
  m_tree->Branch("pz1", &event.pz1, "pz1/D");
  m_tree->Branch("p1",  &event.p1,  "p1/D");  // [MeV/c]
  m_tree->Branch("t1",  &event.t1,  "t1/D");  // [MeV]
#endif

  if(confMan.Get<G4int>("Experiment") == 90) m_tree->Branch("Mt", &event.TPCMt, "Mt/I");

  return;
  //  m_tree->Branch("t0",&event.t0,   "t0/D");
  // m_tree->Branch("Id",&event.Id, "Id/I");
  //   m_tree->Branch("nP",&event.nP,"nP/I");
  //   m_tree->Branch("nN",&event.nN,"nN/I");
  //   m_tree->Branch("nL",&event.nL,"nL/I");
  //   m_tree->Branch("nSm",&event.nSm,"nSm/I");
  //   m_tree->Branch("nSz",&event.nSz,"nSz/I");
  //   m_tree->Branch("nSp",&event.nSp,"nSp/I");
  //   m_tree->Branch("nXm",&event.nXm,"nXm/I");
  //   m_tree->Branch("nXz",&event.nXz,"nXz/I");
  //   m_tree->Branch("nXsm",&event.nXsm,"nXsm/I");
  //   m_tree->Branch("nXsz",&event.nXsz,"nXsz/I");
  //   m_tree->Branch("nPim",&event.nPim,"nPim/I");
  //   m_tree->Branch("nPiz",&event.nPiz,"nPiz/I");
  //   m_tree->Branch("nPip",&event.nPip,"nPip/I");
  //   m_tree->Branch("nKm",&event.nKm,"nKm/I");
  //   m_tree->Branch("nKp",&event.nKp,"nKp/I");

  m_tree->Branch("vdxp", event.SlituDeg, "vdxp[11]/D");
  m_tree->Branch("vdyp", event.SlitvDeg, "vdyp[11]/D");
  /*
    m_tree->Branch("SlitX", event.SlitX, "SlitX[8]/D");
    m_tree->Branch("SlitY", event.SlitY, "SlitY[8]/D");
    m_tree->Branch("Slitt", event.Slitt, "Slitt[8]/D");
    m_tree->Branch("SlitMom", event.SlitMom, "SlitMom[8]/D");
    m_tree->Branch("Slitp", event.Slitp, "Slitp[8]/D");
    m_tree->Branch("SlitNh", event.SlitNh, "SlitNh[8]/I");
    m_tree->Branch("SlitNP", event.SlitNP, "SlitNP[8]/I");
    m_tree->Branch("SlitNK", event.SlitNK, "SlitNK[8]/I");
    m_tree->Branch("SlitNPi", event.SlitNPi, "SlitNPi[8]/I");
    m_tree->Branch("SlitF", event.SlitF, "SlitF[8]/I");
  */

  // Order -->
  // vd8 | TOF | vd9 | AC | vd10 | WC | vd11 (vd[10]) | SDC1
  m_tree->Branch("vdx",   event.SlitX, Form("vdx[11][%d]/D", MaxHits));
  m_tree->Branch("vdy",   event.SlitY, "vdy[11]/D");
  m_tree->Branch("vdtime",event.Slitt, "vdtime[11]/D");
  m_tree->Branch("vdmom", event.SlitMom, "vdmom[11]/D");
  m_tree->Branch("vdpath",event.Slitp, "vdpath[11]/D");
  m_tree->Branch("vdNh",  event.SlitNh, "vdNh[11]/I");
  m_tree->Branch("vdNP",  event.SlitNP, "vdNP[11]/I");
  m_tree->Branch("vdNK",  event.SlitNK, "vdNK[11]/I");
  m_tree->Branch("vdNPi", event.SlitNPi, "vdNPi[11]/I");
  m_tree->Branch("vdPID", event.SlitF, "vdPID[11]/I");
  // m_tree->Branch("TOFAll", &event.TOFAll, "TOFAll/D");
  // m_tree->Branch("TOFt",event.TOFt,   "TOFt[17]/D");
  // m_tree->Branch("TOFtObs",event.TOFtObs,   "TOFtObs[17]/D");
  //m_tree->Branch("tofco",  &event.tofco,   "tofco[17]/I");

  m_tree->Branch("tofnhits",&event.TOFNhits, "tofnhits/I");
  m_tree->Branch("toftime", &event.toftime,  "toftime[18]/D");
  m_tree->Branch("tofdE",   &event.tofdE,    "tofdE[18]/D");

  m_tree->Branch("wcnhits",&event.WCNhits,  "wcnhits/I");
  //m_tree->Branch("wctime", &event.wctime,   "wctime[12]/D");
  //m_tree->Branch("wcdE",   &event.wcdE,     "wcdE[12]/D");
  //m_tree->Branch("wcnpe",  &event.wcnpe,    "wcnpe[12]/D");
  m_tree->Branch("wctime1", &event.wctime1,   "wctime1[6]/D");
  m_tree->Branch("wcdE1",   &event.wcdE1,     "wcdE1[6]/D");
  m_tree->Branch("wcnpe1",  &event.wcnpe1,    "wcnpe1[6]/D");
  m_tree->Branch("wctime2", &event.wctime2,   "wctime2[6]/D");
  m_tree->Branch("wcdE2",   &event.wcdE2,     "wcdE2[6]/D");
  m_tree->Branch("wcnpe2",  &event.wcnpe2,    "wcnpe2[6]/D");

  m_tree->Branch("Q1Trig", &event.Q1Trig,  "Q1Trig/B");
  m_tree->Branch("Q2Trig", &event.Q2Trig,  "Q2Trig/B");
  m_tree->Branch("VDTrig", &event.VDTrig,  "VDTrig/B");
  // m_tree->Branch("TOFTrig",&event.TOFTrig, "TOFTrig/B");

  m_tree->Branch("WCTrig", &event.WCTrig,  "WCTrig/B");
  // m_tree->Branch("TOFHit", &event.TOFHit, "TOFHit/I");
  // m_tree->Branch("ACX", &event.ACX, "ACX/D");
  // m_tree->Branch("ACY", &event.ACY, "ACY/D");
  // m_tree->Branch("ACt", &event.ACt, "ACt/D");
  // m_tree->Branch("ACp", &event.ACp, "ACp/D");
  // m_tree->Branch("ACNhits", &event.ACNhits, "ACNhits/I");
  // m_tree->Branch("ACHit", &event.ACHit, "ACHit/I");
  // m_tree->Branch("WCt", event.WCt,  "WCt[16]/D");
  // m_tree->Branch("WCp", event.WCp,  "WCp[16]/D");
  // m_tree->Branch("WCNhits", &event.WCNhits, "WCNhits/I");
  // m_tree->Branch("WCHit", &event.WCHit, "WCHit/I");
  //G4cout<<"tree and Branch is defined"<<G4endl;

}

//_____________________________________________________________________________
void
S2SAnaManager::SetPrimaryParticle(G4int id, G4int pdg,
                                  const G4LorentzVector& p,
                                  const G4LorentzVector& v,
                                  G4bool is_virtual_beam)
{
  G4int id1 = is_virtual_beam ? -1 : 1;
  G4int id2 = id;
  for (const auto& ptcl: event.hits.at("PRM")) {
    if (ptcl.GetMother(0) == id1 && ptcl.GetMother(1) == id2) {
      G4cerr << FUNC_NAME << " id1=" << id1 << ", id2=" << id2
             << " is already set" << G4endl;
    }
  }
  TParticle particle(pdg,
                     0, // fStatus
                     id1, // fMother[0]
                     id2, // fMother[1]
                     0, // fDaughter[0]
                     0, // fDaughter[1]
                     TLorentzVector(p.px(), p.py(), p.pz(), p.e()),
                     TLorentzVector(v.x(), v.y(), v.z(), v.t()));
  event.hits.at("PRM").push_back(particle);

  if(id == 0 && !noHist){
    hmap.at("PRMPThetaGen")->Fill(p.theta()/CLHEP::degree,
                                  p.v().mag()/CLHEP::GeV);
  }
}

//_____________________________________________________________________________
void
S2SAnaManager::SetGeneratedParticle(const G4String& branch_name,
                                    G4int mother_id, G4int pdg,
                                    const G4LorentzVector& p,
                                    const G4LorentzVector& v)
{
  const TString key(branch_name.c_str());
  auto it = event.hits.find(key);
  if(it == event.hits.end()){
    G4cerr << FUNC_NAME << " unknown branch : " << branch_name << G4endl;
    return;
  }
  const G4int id1 = (mother_id >= 0) ? 1 : -1;
  const G4int id2 = mother_id;
  TParticle particle(pdg,
                     0,            // fStatus
                     id1,          // fMother[0]
                     id2,          // fMother[1]
                     0,            // fDaughter[0]
                     0,            // fDaughter[1]
                     TLorentzVector(p.px(), p.py(), p.pz(), p.e()),
                     TLorentzVector(v.x(), v.y(), v.z(), v.t()));
  it->second.push_back(particle);
}

//_____________________________________________________________________________
void S2SAnaManager::PrintHitsInformation(const G4Event* /* anEvent */,
                                         std::ostream &ost) const
{
  //   G4cout<<"PrintHits is called"<<G4endl;

  //  int GeomFlag = confMan->GeomFlag();

  // G4HCofThisEvent *HCE = anEvent->GetHCofThisEvent();
  // G4SDManager *SDMan = G4SDManager::GetSDMpointer();

  std::ios::fmtflags oldFlags = ost.flags();
  std::size_t preSiz = ost.precision();
  ost.setf( std::ios::fixed );

  ost.precision(5);
  ost << event.x0In << std::setw(10)
      << event.y0In << std::setw(10)
      <<" "<< -4320.62-event.z0In << std::setw(10);
  ost.precision(5);

  ost << event.p0     << std::setw(10)
      << event.theta0 << std::setw(10)
      << event.phi0   << std::setw(10);
  ost.precision(7);

  ost << 1.8 <<std::setw(10);
  ost << 0   <<std::setw(10);
  ost << 0   <<std::setw(10);

  ost.precision(5);

  //G4cout<<"Printed"<<G4endl;
  //G4cout<<"osf="<<ost<<G4endl;
  //G4cout<<"x0="<<event.x0In<<" y0="<<event.y0In<<G4endl;

  // const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  // G4double nhTOF=0;
  // TOFHitsCollection    *TOFHC;
  // G4int colIdTOF = SDMan->GetCollectionID( "TOFSD"/*"DCCollection"*/ );
  // TOFHC =  dynamic_cast<TOFHitsCollection *>(  HCE->GetHC( colIdTOF ) );
  // if( TOFHC )     nhTOF     = TOFHC ->entries();

  // for( int i=0; i<nhTOF; ++i ){
  //   TOFHit *aHit = (*TOFHC)[i];
  //   double time = aHit->GetTime();
  //   double de = aHit->GetEdep();
  //   if(i==nhTOF-1){
  //     ost << std::setw(12) << 61
  //         << std::setw(12) << time;
  //     de=1.0;
  //     ost << std::setw(12) << 62 <<" "
  //         << std::setw(12) << de;
  //   }
  // }

  ost << std::setw(5) << -1 << std::endl;

  ost.flags( oldFlags );
  ost.precision( preSiz );

  ost << std::endl;

}
