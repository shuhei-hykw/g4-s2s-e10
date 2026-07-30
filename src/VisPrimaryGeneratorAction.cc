// -*- C++ -*-

#include "VisPrimaryGeneratorAction.hh"

#include <G4Event.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4SystemOfUnits.hh>

//_____________________________________________________________________________
VisPrimaryGeneratorAction::VisPrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(),
    m_particleGun(new G4ParticleGun(1))
{
  auto particle =
    G4ParticleTable::GetParticleTable()->FindParticle("geantino");
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
  m_particleGun->SetParticleEnergy(1.*GeV);
  m_particleGun->SetParticlePosition(G4ThreeVector(0., 0., -10.*m));
}

//_____________________________________________________________________________
VisPrimaryGeneratorAction::~VisPrimaryGeneratorAction()
{
  delete m_particleGun;
}

//_____________________________________________________________________________
void
VisPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  m_particleGun->GeneratePrimaryVertex(anEvent);
}
