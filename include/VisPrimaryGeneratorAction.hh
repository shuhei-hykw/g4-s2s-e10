// -*- C++ -*-
//
// Trivial primary generator for the visualization-only binary
// (G4S2SVis). No analysis-manager interaction, unlike
// S2SPrimaryGeneratorAction -- this only exists to satisfy Geant4's
// mandatory user-action requirement while rendering static geometry.

#ifndef VIS_PRIMARY_GENERATOR_ACTION_HH
#define VIS_PRIMARY_GENERATOR_ACTION_HH

#include <G4VUserPrimaryGeneratorAction.hh>

class G4ParticleGun;
class G4Event;

//_____________________________________________________________________________
class VisPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  VisPrimaryGeneratorAction();
  ~VisPrimaryGeneratorAction() override;
  void GeneratePrimaries(G4Event* anEvent) override;

private:
  G4ParticleGun* m_particleGun;
};

#endif
