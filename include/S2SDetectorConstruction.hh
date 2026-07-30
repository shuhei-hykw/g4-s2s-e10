// -*- C++ -*-

#ifndef S2S_DETECTOR_CONSTRUCTION_HH
#define S2S_DETECTOR_CONSTRUCTION_HH 1

#include <map>

#include <G4VUserDetectorConstruction.hh>
#include <G4String.hh>
#include <globals.hh>

class G4Box;
class G4Element;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class G4VSensitiveDetector;

struct MaterialList;

//_____________________________________________________________________________
class S2SDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  static G4String ClassName();
  S2SDetectorConstruction();
  ~S2SDetectorConstruction();

public:
  G4VPhysicalVolume* Construct();

private:
  G4int                           m_experiment;
  G4bool                          m_check_overlaps;
  // E10-2nd uses the shared E90/HypTPC setup (HypTPC + HTOF + E90SAC,
  // solid target, SDC1 -- NOT the old (pi-,K+)-era SSD/SFT). Enabled
  // by conf key "HypTPCSetup:1"; 0/unset keeps the legacy E10 geometry.
  G4bool                          m_hyptpc_setup;
  G4LogicalVolume*                m_world_lv;
  static std::vector<G4String>    s_detector_list;

  void AddNewDetector(G4VSensitiveDetector* sd);
  void ConstructField();
  //     void ConstructCalorimeter();
  void ConstructBAC();
  void ConstructSAC();
  void ConstructTarget ();
  void ConstructTgtHeBag();
  void ConstructQ1();
  void ConstructQ2();
  void ConstructD1();
  void ConstructSDC1();
  void ConstructSDC2();
  void ConstructSDC3();
  void ConstructSDC4();
  void ConstructSDC5();
  void ConstructSFT();
  void ConstructSSD();
  void ConstructTOF();
  void ConstructAC1();
  void ConstructWC();
  void ConstructVP();
  void ConstructKLChamber(G4int i);
  void ConstructKLChamberMylar(G4int i);

  void ConstructHTOF();
  void ConstructHypTPC();
  void ConstructTargetE90();
  void ConstructTargetE10();

  void ConstructHBXXGe();
  void ConstructNaI();
  void ConstructRC(G4int i);
  void ConstructPD(G4int i, G4String YZ);
};

//_____________________________________________________________________________
inline G4String
S2SDetectorConstruction::ClassName()
{
  static G4String s_name("S2SDetectorConstruction");
  return s_name;
}

#endif
