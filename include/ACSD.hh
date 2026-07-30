// -*- C++ -*-

#ifndef AC_SD_HH
#define AC_SD_HH 1

#include <G4VSensitiveDetector.hh>

#include "ACHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

//_____________________________________________________________________________
class ACSD : public G4VSensitiveDetector
{
public:
  ACSD(const G4String& name);
  ~ACSD();

private:
  ACHitsCollection *ACCollection;
  G4double m_refractive_index;

public:
  void EndOfEvent(G4HCofThisEvent* HCE);
  void Initialize(G4HCofThisEvent* HCE);
  G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist);
  void SetRefractiveIndex(G4double n){ m_refractive_index = n; }
};

#endif
