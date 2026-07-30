// -*- C++ -*-

#include "VHitInfo.hh"

#include <G4ParticleTable.hh>
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4TouchableHistory.hh>
#include <G4TouchableHandle.hh>

#include <TLorentzVector.h>
#include <TParticle.h>

#include "ConfMan.hh"
#include "FuncName.hh"
#include "PrintHelper.hh"

namespace
{
const auto particleTable = G4ParticleTable::GetParticleTable();
}

//_____________________________________________________________________________
VHitInfo::VHitInfo(const G4String& name, G4Step* step,
                   G4bool use_center_point)
  : m_detector_name(name),
    m_particle_name(),
    m_position(),
    //m_lposition(),
    m_momentum(),
    m_time(),
    m_energy_deposit(),
    m_track_id(),
    m_pdg_encoding(),
    m_copy_number(),
    m_step_length(),
    m_mass(),
    m_charge(),
    m_parent_id(),
    m_track_length(),
    m_vertex_position(),
    m_vertex_momentum(),
    m_vertex_kinetic_energy(),
    m_particle()
{
  if (!step) return;
  const auto track = step->GetTrack();
  const auto pre = step->GetPreStepPoint();
  const auto post = step->GetPostStepPoint();
  const auto dp = track->GetDynamicParticle();
  G4TouchableHandle theTouchable = pre->GetTouchableHandle();
  m_particle_name = track->GetDefinition()->GetParticleName();
  if (use_center_point) {
    m_position = (pre->GetPosition() + post->GetPosition())/2.;
    m_lposition = theTouchable->GetHistory()->GetTopTransform().TransformPoint(m_position);
    m_momentum = (pre->GetMomentum() + post->GetMomentum())/2.;
  } else {
    m_position = pre->GetPosition();
    m_lposition = theTouchable->GetHistory()->GetTopTransform().TransformPoint(m_position);
    m_momentum = pre->GetMomentum();
  }
  m_time = pre->GetGlobalTime();
  m_energy_deposit = step->GetTotalEnergyDeposit();
  m_track_id = track->GetTrackID();
  m_pdg_encoding = track->GetDefinition()->GetPDGEncoding();
  m_copy_number = pre->GetPhysicalVolume()->GetCopyNo();
  m_step_length = step->GetStepLength();
  m_mass = track->GetDynamicParticle()->GetMass();
  m_charge = track->GetDynamicParticle()->GetCharge();
  m_parent_id = track->GetParentID();
  m_track_length = track->GetTrackLength();
  m_vertex_position = track->GetVertexPosition();
  m_vertex_momentum = track->GetVertexMomentumDirection();
  m_vertex_kinetic_energy = track->GetVertexKineticEnergy();
  TLorentzVector p(m_momentum.x(), m_momentum.y(), m_momentum.z(),
                   dp->GetTotalEnergy());
  //TLorentzVector v(m_position.x(), m_position.y(), m_position.z(),
  //                 pre->GetGlobalTime());
  TLorentzVector v(m_lposition.x(), m_lposition.y(), m_lposition.z(),
                   pre->GetGlobalTime());
  const auto experiment = ConfMan::GetInstance().Get<G4int>("Experiment");
  // E10 also needs the track id (fStatusCode) offline, to match a
  // single physical track across TOF/AC1/WC hits (e.g. for
  // reconstructing a missing-mass spectrum from a background
  // generator's decay secondaries, which are not primaries).
  const auto status_code =
    (experiment == 90 || experiment == 10) ? m_track_id : 0;
  m_particle = new TParticle(m_pdg_encoding,
                             status_code,      // fStatusCode / track id for E90/E10
                             m_parent_id,      // fMother[0]
                             m_copy_number,    // fMother[1]
                             0,                // fDaughter[0]
                             0,                // fDaughter[1]
                             p, v);
  m_particle->SetWeight(m_energy_deposit);
}

//_____________________________________________________________________________
VHitInfo::~VHitInfo()
{
  if (m_particle) delete m_particle;
}

//_____________________________________________________________________________
void
VHitInfo::SetEnergyDeposit(G4double edep)
{
  m_energy_deposit = edep;
  if (m_particle) m_particle->SetWeight(edep);
}

//_____________________________________________________________________________
G4bool
VHitInfo::Is(const G4String& particle_name) const
{
  return Is(particleTable->FindParticle(particle_name)->GetPDGEncoding());
}

//_____________________________________________________________________________
G4bool
VHitInfo::Is(G4int pdg_encoding) const
{
  return pdg_encoding == m_pdg_encoding;
}

//_____________________________________________________________________________
G4bool
VHitInfo::IsPrimary() const
{
  return (m_parent_id == 0) && (m_track_id == 1);
}

//_____________________________________________________________________________
G4bool
VHitInfo::IsWeakPi() const // E63 for weak pi-
{
  return (m_parent_id == 0) && (m_track_id == 2);
}

//_____________________________________________________________________________
void
VHitInfo::Print() const
{
  PrintHelper helper(3, std::ios::fixed, G4cout);
  G4cout << FUNC_NAME << " " << m_detector_name << G4endl
	 << "   ParticleName        = " << m_particle_name << G4endl
	 << "   Position            = " << m_position*(1/CLHEP::mm)
	 << " mm" << G4endl
	 << "   Momentum            = " << m_momentum*(1/CLHEP::GeV)
	 << " GeV/c" << G4endl
	 << "   Time                = " << m_time/CLHEP::ns << " ns" << G4endl
	 << "   EnergyDeposit       = " << m_energy_deposit/CLHEP::MeV
	 << " MeV" << G4endl
	 << "   TrackID             = " << m_track_id << G4endl
	 << "   PDGEncoding         = " << m_pdg_encoding << G4endl
	 << "   CopyNumber          = " << m_copy_number << G4endl
	 << "   StepLength          = " << m_step_length/CLHEP::mm
	 << " mm" << G4endl
	 << "   Mass                = " << m_mass/CLHEP::MeV
	 << " MeV" << G4endl
	 << "   Charge              = " << m_charge << G4endl
	 << "   ParentID            = " << m_parent_id << G4endl
	 << "   TrackLength         = " << m_track_length/CLHEP::mm
	 << " mm" << G4endl
	 << "   VertexPosition      = " << m_vertex_position*(1/CLHEP::mm)
	 << " mm" << G4endl
	 << "   VertexMomentum      = " << m_vertex_momentum
	 << " GeV/c" << G4endl
	 << "   VertexKineticEnergy = " << m_vertex_kinetic_energy/CLHEP::MeV
	 << " MeV" << G4endl;
}
