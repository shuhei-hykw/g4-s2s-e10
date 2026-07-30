// -*- C++ -*-

#include "S2SDetectorConstruction.hh"

#include <string>

#include <G4FieldManager.hh>
#include <G4ChordFinder.hh>
#include <G4TransportationManager.hh>

#include <G4Polyhedra.hh>
#include <G4Material.hh>
#include <G4Tubs.hh>
#include <G4EllipticalTube.hh>
#include <G4Box.hh>
#include <G4Trd.hh>
#include <G4ExtrudedSolid.hh>
#include <G4Hype.hh>
#include <G4IntersectionSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4UnionSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4PVReplica.hh>
#include <G4UniformMagField.hh>
#include <G4Mag_UsualEqRhs.hh>
#include <G4SimpleRunge.hh>
#include <G4VisAttributes.hh>

#include <G4GeometryManager.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4SolidStore.hh>
#include <G4UserLimits.hh>
#include <G4SystemOfUnits.hh>

#include <G4Colour.hh>
#include <G4SDManager.hh>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "DetSizeMan.hh"
#include "MagnetConstant.hh"
#include "MaterialList.hh"
#include "S2SField.hh"

#include "DCSD.hh"
#include "VPSD.hh"
#include "TOFSD.hh"
#include "ACSD.hh"
#include "WCSD.hh"

namespace
{
const auto& confMan = ConfMan::GetInstance();
const auto& geomMan = DCGeomMan::GetInstance();
const auto& sizeMan = DetSizeMan::GetInstance();
const G4bool use_Tgthebag = false;
const G4bool use_Q1hebag = true;
const G4bool use_Q2hebag = true;
const G4bool use_D1hebag = true;
const G4bool use_SDCouthebag = true;
const MaterialList mlist;
G4PVPlacement* physWorld;
// G4Material* matlHeBag = mlist.at("HeGas");
inline G4Material* mat(const G4String& name) {
  return G4Material::GetMaterial(name, /*warning*/ true);
}

inline G4bool is_e63(G4int experiment)
{
  return experiment == 63;
}

inline G4int generator_id()
{
  return confMan.Get<G4int>("Generator");
}

inline G4bool should_construct_hbxx_ge(G4int experiment)
{
  return is_e63(experiment);
}

}

std::vector<G4String> S2SDetectorConstruction::s_detector_list;

//_____________________________________________________________________________
S2SDetectorConstruction::S2SDetectorConstruction()
  : G4VUserDetectorConstruction(),
    m_experiment(),
    m_world_lv(),
    m_check_overlaps(false)
{
}

//_____________________________________________________________________________
S2SDetectorConstruction::~S2SDetectorConstruction()
{
}

//_____________________________________________________________________________
G4VPhysicalVolume* S2SDetectorConstruction::Construct()
{
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  m_experiment = confMan.Get<G4int>("Experiment");
  m_check_overlaps = confMan.Get<G4bool>("CheckOverlap");
  // "HypTPCSetup:1" -> E10-2nd shared E90/HypTPC geometry. Missing key
  // returns 0 (legacy E10). E90 always uses this setup regardless.
  m_hyptpc_setup = (confMan.Get<G4int>("HypTPCSetup") != 0);

  ///// World
  const auto& half_size = sizeMan.GetSize("World")*mm/2.;
  auto solidWorld = new G4Box
    ("World", half_size.x(), half_size.y(), half_size.z());
  m_world_lv = new G4LogicalVolume(solidWorld, mlist.at("Air"), "World");
  m_world_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  physWorld = new G4PVPlacement
    (nullptr, G4ThreeVector(), m_world_lv, "World",
     nullptr, false, 0, m_check_overlaps);

#if 1
  ConstructBAC();
#endif

#if 1
  if(m_experiment == 63){
    ConstructPD(1, "Y");
    ConstructPD(1, "Z");
    ConstructRC(1);
    ConstructRC(2);
    ConstructPD(2, "Y");
    ConstructPD(2, "Z");
  }
#endif

#if 1
  // E10-2nd shares the E90/HypTPC setup (HypTPC + HTOF barrel + E90SAC),
  // differing only in the target -- built when HypTPCSetup:1 (or E90).
  // Requires the DCGEO/DSIZE to carry the HypTPC/HTOF/SAC/TargetHolder
  // keys (see param/DCGEO/DCGeomParam_e10_9hesigma_kpi_hyptpc and
  // param/DSIZE/DetSize_E10_9hesigma_kpi_hyptpc). See analysis-note.md
  // 2026-07-23 HypTPC-for-E10 entry.
  if(m_experiment == 90 || m_hyptpc_setup){
    ConstructHTOF();
    ConstructHypTPC();
    ConstructSAC();
  }
#endif

#if 1
  ConstructTarget();
  ConstructTgtHeBag();
#endif

#if 1
  if(should_construct_hbxx_ge(m_experiment))
    ConstructHBXXGe();
#endif

#if 1
  ConstructQ1();
  ConstructQ2();
  ConstructD1();
  ConstructField();
#endif

#if 1
  // Near-target tracking. Legacy E10 (pi-,K+ era) used SSD+SFT; the
  // E10-2nd HypTPC setup (and E90) instead tracks near the target with
  // HypTPC and uses SDC1, matching the shared setup -- SSD/SFT would
  // physically overlap the HypTPC/HTOF/SAC barrel and are dropped.
  if(m_experiment == 10 && !m_hyptpc_setup){
    ConstructSSD();
    ConstructSFT();
  }else{
    ConstructSDC1();
  }
  ConstructSDC2();
#endif

#if 1
  ConstructSDC3();
  ConstructSDC4();
  ConstructSDC5();
#endif

#if 1
  ConstructTOF();
  ConstructAC1();
  ConstructWC();
#endif

#if 1
  ConstructVP();
#endif

  return physWorld;
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::AddNewDetector(G4VSensitiveDetector* sd)
{
  G4SDManager::GetSDMpointer()->AddNewDetector(sd);
  s_detector_list.push_back(sd->GetName());
}

//_____________________________________________________________________________
void S2SDetectorConstruction::ConstructField()
{
  S2SField *field = new S2SField(confMan.Get<G4String>("FLDMAP"));
  auto fieldManager =
    G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldManager->SetDetectorField(field);
  fieldManager->CreateChordFinder(field);

  // G4Mag_UsualEqRhs* fEquation = new G4Mag_UsualEqRhs(field);
  // G4MagIntegratorStepper *pStepper = new G4SimpleRunge(fEquation);
  // G4ChordFinder *pChordFinder = new G4ChordFinder(field, 1.0E-3*mm, pStepper);
  // fieldManager->SetChordFinder(pChordFinder);

  fieldManager->GetChordFinder()->SetDeltaChord(1.e-3*mm);
  //fieldManager->SetDeltaIntersection(1.0E-6*mm);
  //fieldManager->SetDeltaOneStep(1.0E-4*mm);
  //fieldManager->SetMaximumEpsilonStep(1.0E-3);
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructBAC()
{
  if(use_Tgthebag) return;
  const auto& ra2 = geomMan.GetRotAngle2("BAC1") * deg;
  const auto& frame_size = sizeMan.GetSize("BacFrame") * 0.5 * mm;
  const auto& radiator_size = sizeMan.GetSize("Bac1Radiator") * 0.5 * mm;
  const G4bool e63_bac_pair = (m_experiment == 63);
  // Mother
  auto mother_solid = new G4Box("BacMotherSolid",
                                frame_size.x() + 5.*mm,
                                frame_size.y() + 5.*mm,
                                frame_size.z() + 5.*mm);
  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       mlist.at("Air"),
                                       "BacMotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  auto pos = 0.5*(geomMan.GetGlobalPosition("BAC1") + geomMan.GetGlobalPosition("BAC2"));
  G4ThreeVector offset(0., 0., 0.);
  if(!e63_bac_pair)
    offset.setZ(frame_size.z()-0.6*mm - radiator_size.z());
  offset.rotateY(ra2);
  new G4PVPlacement(rot, pos + offset,
                    "BacMotherPHYS", mother_lv, physWorld, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // Frame
  auto frame_solid = new G4Box("BacFrameSolid", frame_size.x(),
                               frame_size.y(), frame_size.z());
  auto frame_lv = new G4LogicalVolume(frame_solid,
                                      mlist.at("Air"),
                                      "BacFrameLV");
  frame_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  pos.setMag(0.);
  new G4PVPlacement(nullptr, pos, frame_lv,
                    "BacFramePHYS", mother_lv, false, 0, m_check_overlaps);
  const G4double reflector_thickness = 0.3*mm;
  const G4double black_sheet_thickness = 0.2*mm; // both tmp value
  const G4double wrap_clearance = 0.02*mm;

  auto radiator_solid = new G4Box("BacRadiatorSolid", radiator_size.x(),
                                  radiator_size.y(), radiator_size.z());
  auto radiator_material =
    (m_experiment == 63) ? mlist.AerogelBAC1E63 : mlist.Aerogel;
  auto radiator_lv = new G4LogicalVolume(radiator_solid,
                                         radiator_material, "BacRadiatorLV");
  radiator_lv->SetVisAttributes(G4Color::Magenta());
  auto place_rectangular_wrap =
    [&](const G4String& prefix, G4Material* material, const G4Colour& colour,
        const G4ThreeVector& center, const G4ThreeVector& inner_half,
        G4double thickness)
  {
    auto x_solid = new G4Box(prefix + "XSideSolid",
                             0.5*thickness,
                             inner_half.y() + thickness,
                             inner_half.z() + thickness);
    auto y_solid = new G4Box(prefix + "YSideSolid",
                             inner_half.x(),
                             0.5*thickness,
                             inner_half.z() + thickness);
    auto z_solid = new G4Box(prefix + "ZSideSolid",
                             inner_half.x(),
                             inner_half.y(),
                             0.5*thickness);
    auto x_lv = new G4LogicalVolume(x_solid, material, prefix + "XSideLV");
    auto y_lv = new G4LogicalVolume(y_solid, material, prefix + "YSideLV");
    auto z_lv = new G4LogicalVolume(z_solid, material, prefix + "ZSideLV");
    x_lv->SetVisAttributes(colour);
    y_lv->SetVisAttributes(colour);
    z_lv->SetVisAttributes(colour);

    for(G4int iside=0; iside<2; ++iside){
      const G4double sign = (iside == 0) ? -1. : 1.;
      new G4PVPlacement(
        nullptr,
        center + G4ThreeVector(sign*(inner_half.x() + 0.5*thickness),
                               0.0, 0.0),
        x_lv, prefix + "XSidePHYS", frame_lv, false, iside,
        m_check_overlaps);
      new G4PVPlacement(
        nullptr,
        center + G4ThreeVector(0.0,
                               sign*(inner_half.y() + 0.5*thickness),
                               0.0),
        y_lv, prefix + "YSidePHYS", frame_lv, false, iside,
        m_check_overlaps);
      new G4PVPlacement(
        nullptr,
        center + G4ThreeVector(0.0, 0.0,
                               sign*(inner_half.z() + 0.5*thickness)),
        z_lv, prefix + "ZSidePHYS", frame_lv, false, iside,
        m_check_overlaps);
    }
  };

  if(e63_bac_pair){
    // The 2mm center gap has a cardboard-like divider in the real detector;
    // modeled as air in G4.
    const G4double gap_half = 1.0*mm;
    const G4double end_offset = 21.0*mm;
    const G4double radiator_center_dz = gap_half + radiator_size.z();

    // BAC1 (upstream, copy 0) at -z, BAC2 (downstream, copy 1) at +z
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -radiator_center_dz),
                      radiator_lv, "Bac1RadiatorPHYS", frame_lv, false, 0,
                      m_check_overlaps);
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., +radiator_center_dz),
                      radiator_lv, "Bac2RadiatorPHYS", frame_lv, false, 1,
                      m_check_overlaps);

    // Teflon reflector: one wrap per radiator block (not shared)
    const G4ThreeVector teflon_inner(
      radiator_size.x() + wrap_clearance,
      radiator_size.y() + wrap_clearance,
      radiator_size.z() + wrap_clearance);

    place_rectangular_wrap(
      "Bac1Reflector", mlist.Teflon, G4Color::White(),
      G4ThreeVector(0., 0., +radiator_center_dz),
      teflon_inner,
      reflector_thickness - wrap_clearance);

    place_rectangular_wrap(
      "Bac2Reflector", mlist.Teflon, G4Color::White(),
      G4ThreeVector(0., 0., -radiator_center_dz),
      teflon_inner,
      reflector_thickness - wrap_clearance);

    // Black sheet: one envelope covering entire assembly incl. 21mm end offsets
    const G4double assembly_half_z =
      radiator_center_dz + radiator_size.z() + end_offset;
    const G4ThreeVector black_inner(
      radiator_size.x() + reflector_thickness + wrap_clearance,
      radiator_size.y() + reflector_thickness + wrap_clearance,
      assembly_half_z + wrap_clearance);

    place_rectangular_wrap(
      "BacBlackSheet", mlist.PVC, G4Colour(0.02, 0.02, 0.02),
      G4ThreeVector(), black_inner,
      black_sheet_thickness - wrap_clearance);

  } else {
    const G4ThreeVector radiator_pos(
      0., 0., -frame_size.z()+0.6*mm + radiator_size.z());
    new G4PVPlacement(nullptr, radiator_pos, radiator_lv,
                      "Bac1RadiatorPHYS", frame_lv, false, 0,
                      m_check_overlaps);
    place_rectangular_wrap(
      "Bac1Reflector", mlist.Teflon, G4Color::White(),
      radiator_pos,
      radiator_size + G4ThreeVector(wrap_clearance,
                                    wrap_clearance,
                                    wrap_clearance),
      reflector_thickness - wrap_clearance);

    place_rectangular_wrap(
      "Bac1BlackSheet", mlist.PVC, G4Colour(0.02, 0.02, 0.02),
      radiator_pos,
      radiator_size + G4ThreeVector(reflector_thickness + wrap_clearance,
                                    reflector_thickness + wrap_clearance,
                                    reflector_thickness + wrap_clearance),
      black_sheet_thickness - wrap_clearance);
  }

}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructSAC()
{
  // E90SAC is part of the shared E90/HypTPC setup -- build it for the
  // E10-2nd HypTPC setup too, not only E90. See Construct() gating.
  if(m_experiment != 90 && !m_hyptpc_setup) return;
  auto sac_sd = new ACSD("SAC");
  sac_sd->SetRefractiveIndex(1.05);
  G4SDManager::GetSDMpointer()->AddNewDetector(sac_sd);
  const auto& ra2 = geomMan.GetRotAngle2("SAC") * deg;
  const auto& frame_size = sizeMan.GetSize("SacFrame") * 0.5 * mm;
  const auto& radiator_size = sizeMan.GetSize("SacRadiator") * 0.5 * mm;
  // Mother
  auto mother_solid = new G4Box("SacMotherSolid",
                                frame_size.x() + 5.*mm,
                                frame_size.y() + 5.*mm,
                                frame_size.z() + 5.*mm);
  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       mlist.at("Air"),
                                       "SacMotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  auto pos = geomMan.GetGlobalPosition("SAC");
  G4ThreeVector offset(0., 0., frame_size.z()-0.6*mm - radiator_size.z());
  offset.rotateY(ra2);
  new G4PVPlacement(rot, pos + offset, "SacMotherPV",
                    mother_lv, physWorld, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // Frame
  auto frame_solid = new G4Box("SacFrameSolid", frame_size.x(),
                               frame_size.y(), frame_size.z());
  auto frame_lv = new G4LogicalVolume(frame_solid,
                                      mlist.at("Air"),
                                      "SacFrameLV");
  // frame_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  frame_lv->SetVisAttributes(G4Color::Gray());
  pos.setMag(0.);
  new G4PVPlacement(nullptr, pos, frame_lv,
                    "SacFramePV", mother_lv, false, 0, m_check_overlaps);
  // Radiator
  auto radiator_solid = new G4Box("SacRadiatorSolid", radiator_size.x(),
                                  radiator_size.y(), radiator_size.z());
  auto radiator_lv = new G4LogicalVolume(radiator_solid,
                                         mlist.Aerogel,
                                         "SacRadiatorLV");
  radiator_lv->SetSensitiveDetector(sac_sd);
  radiator_lv->SetVisAttributes(G4Color::Magenta());
  // pos.set(0., 0., -frame_size.z()+0.6*mm + radiator_size.z());
  new G4PVPlacement(nullptr, pos, radiator_lv,
                    "SacRadiatorPV", frame_lv, false, 0);
#if 0
  // Foward Black sheet
  const auto sheet_thickness = 0.2*mm/2;
  G4VSolid* solidBlackSheet1;
  solidBlackSheet1 = new G4Box("solidSacBlackSheet1",frame_size.x(),
                               frame_size.y(), sheet_thickness);
  solidBlackSheet1 = new G4IntersectionSolid("solidSacBlackSheet1", solidBlackSheet1, frame_solid);
  auto logicBlackSheet1 = new G4LogicalVolume(solidBlackSheet1, mlist.PVC, "logicBlackSheet1");
  G4RotationMatrix rot_sheet;
  logicBlackSheet1->SetVisAttributes(G4Color::Gray());
  new G4PVPlacement(G4Transform3D(rot_sheet, G4ThreeVector(0, 0, -frame_size.z()+sheet_thickness)),
                    logicBlackSheet1, "physSacBlackSheet1", frame_lv, false, 0, m_check_overlaps);

  // Foreward Reflector w/ teflon
  const G4double reflector_thickness = 0.3*mm/2;
  const G4ThreeVector reflector1_size(frame_size.x(),
					frame_size.y(),
                                  	reflector_thickness);
  auto reflector1_solid = new G4Box("SacReflector1Solid", reflector1_size.x(),
                                 reflector1_size.y(), reflector1_size.z());
  auto reflector1_lv = new G4LogicalVolume(reflector1_solid,
                                        mlist.Teflon,
					"SacReflector1LV");
  reflector1_lv->SetVisAttributes(G4Color::White());
  pos.set(0.,0., -frame_size.z()+sheet_thickness*2 + reflector_thickness);
  auto rot_reflector1 = new G4RotationMatrix;
  new G4PVPlacement(rot_reflector1, pos, reflector1_lv,
                      "SacReflector1PHYS", frame_lv, false, 0);

  // Backward Black sheet  w/ PVC
  const G4ThreeVector triangle_size(frame_size.x(), frame_size.y(), frame_size.z()-radiator_size.z()-0.3*mm);
  const G4double BlackSheet2_angle = std::atan2(triangle_size.z(),
                                           triangle_size.y());
  const G4ThreeVector BlackSheet2_size(triangle_size.x(),
				  std::hypot(triangle_size.y(),
                                   triangle_size.z()),
                                  sheet_thickness);
  auto solidBlackSheet2 = new G4Box("solidSacBlackSheet2", BlackSheet2_size.x(),
                                 BlackSheet2_size.y(), BlackSheet2_size.z());
  auto logicBlackSheet2 = new G4LogicalVolume(solidBlackSheet2,
                                        mlist.PVC,
					"logicSacBlackSheet2");
  logicBlackSheet2->SetVisAttributes(G4Color::Gray());
  pos.set(0.,0., -frame_size.z()+radiator_size.z()*2 + triangle_size.z());
  auto rot_BlackSheet2 = new G4RotationMatrix;
  rot_BlackSheet2->rotateX(BlackSheet2_angle);
  new G4PVPlacement(rot_BlackSheet2, pos, logicBlackSheet2,
                      "physSacBlackSheet1", frame_lv, false, 0);


  // Backward Reflector  w/ Teflon
  const G4ThreeVector reflector2_size(triangle_size.x(),
				  std::hypot(triangle_size.y(),
                                   triangle_size.z()),
                                  reflector_thickness);
  auto reflector2_solid = new G4Box("SacReflector2Solid", reflector2_size.x(),
                                 reflector2_size.y(), reflector2_size.z());
  auto reflector2_lv = new G4LogicalVolume(reflector2_solid,
                                        mlist.Teflon,
					"SacReflector2LV");
  reflector2_lv->SetVisAttributes(G4Color::White());
  pos.set(0.,0., -frame_size.z()+radiator_size.z()*2 + triangle_size.z()-0.5*mm);
  new G4PVPlacement(rot_BlackSheet2, pos, reflector2_lv,
                      "SacReflector2PHYS", frame_lv, false, 0);
#endif
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructTarget()
{
  if(use_Tgthebag) return;
  const auto& half_size = sizeMan.GetSize("Target")*mm/2.;
  G4Material *TargetMater = nullptr;
  auto Target = confMan.Get<G4String>("TargetMaterial");
  if(Target == "Be"){
    TargetMater = mlist.at("Be9");}
  else if(Target == "natLi"){
    TargetMater = mlist.at("natLi");}
  else if(Target == "CH2"){
    TargetMater = mlist.Polyethylene;}
  else if(Target == "HeGas"){
    TargetMater = mlist.at("HeGas");}
  else if(Target == "AFT"){
    TargetMater = mlist.at("Air");}
  else if(Target == "Air") { TargetMater = mat("Air"); }
  else if(Target == "LD2") { TargetMater = mat("LD2"); }
  else if(Target == "C") { TargetMater = mlist.at("C"); }
  else{
    G4cout << " Sorry, Target: " << Target
           << " is not defined. So Air will be used. " << G4endl;
    TargetMater = mlist.at("Air");
  }

  switch(m_experiment){
  case 90:
    ConstructTargetE90();
    break;
  case 10:
    // E10-2nd shares the E90/HypTPC setup and differs ONLY in the
    // target: a solid 9Be plate (sized to fit the HypTPC target
    // holder bore) at the SAME position as the E90 target, not the
    // E90 LD2 tube. See ConstructTargetE10() and analysis-note.md
    // 2026-07-23 HypTPC-for-E10 entry.
    ConstructTargetE10();
    break;
  default:
    G4double targetGapHalfZ = 0.0;
    try {
      targetGapHalfZ = confMan.Get<G4double>("TargetGapHalfZ") * mm;
    } catch(...) {
      targetGapHalfZ = 0.0;
    }
    G4RotationMatrix rotTarget;
    const auto& pos = geomMan.GetGlobalPosition("Target");
    if (targetGapHalfZ > 0.0) {
      // 2-block target with central gap (oldli physical setup).
      // Total target half-z = half_size.z(); each block half-z is reduced by
      // gapHalfZ. Block centers sit at +/- (gapHalfZ + block_half_z).
      G4double blockHalfZ = (half_size.z() - targetGapHalfZ) / 2.0;
      G4double blockCenterOffset = targetGapHalfZ + blockHalfZ;
      auto TargetBlock1 = new G4Box("TargetBlock1",
                                    half_size.x(), half_size.y(), blockHalfZ);
      auto TargetBlock2 = new G4Box("TargetBlock2",
                                    half_size.x(), half_size.y(), blockHalfZ);
      auto logTarget1 = new G4LogicalVolume(TargetBlock1, TargetMater, "logTargetBlock1");
      auto logTarget2 = new G4LogicalVolume(TargetBlock2, TargetMater, "logTargetBlock2");
      G4ThreeVector pos1(pos.x(), pos.y(), pos.z() - blockCenterOffset);
      G4ThreeVector pos2(pos.x(), pos.y(), pos.z() + blockCenterOffset);
      new G4PVPlacement(G4Transform3D(rotTarget, pos1), "physTargetBlock1",
                        logTarget1, physWorld, false, 0, m_check_overlaps);
      new G4PVPlacement(G4Transform3D(rotTarget, pos2), "physTargetBlock2",
                        logTarget2, physWorld, false, 1, m_check_overlaps);
      logTarget1->SetVisAttributes(G4Color::Gray());
      logTarget2->SetVisAttributes(G4Color::Gray());
      // Print summary so the run log records the geometry actually placed.
      G4cout << "[Target] 2-block geometry: gapHalfZ=" << targetGapHalfZ/mm
             << " mm, blockHalfZ=" << blockHalfZ/mm << " mm, centers at z="
             << (pos.z() - blockCenterOffset)/mm << " and "
             << (pos.z() + blockCenterOffset)/mm << " mm" << G4endl;
    }
    auto TargetBox = new G4Box
      ("TargetBox", half_size.x(), half_size.y(), half_size.z());
    auto logTarget = new G4LogicalVolume(TargetBox, TargetMater, "logTarget");
    if (targetGapHalfZ <= 0.0) {
      new G4PVPlacement(G4Transform3D(rotTarget, pos),
                        "physTarget",
                        logTarget,
                        physWorld,
                        false,
                        0,
                        m_check_overlaps);
      logTarget->SetVisAttributes(G4Color::Gray());
    }
    if(Target == "AFT"){
      // ===== Define fiber targets ================================
      G4Material *Fiber_Core_Material =mlist.Scin;   // Polystylene
      G4Material *Fiber_Clad_Material =mlist.Acrylic; // PMMA
      // ~~~~~~~ Create logical volumes of x layers ~~~~~~~~~~
      double fiber_phi_core = 2.91*mm;
      double fiber_phi_clad = 0.09*mm;
      double fiber_phi = fiber_phi_core + fiber_phi_clad;


      G4EllipticalTube* fiber_core_solid_x
        = new G4EllipticalTube("fiber_core_solid_x",
                               fiber_phi_core/2.0, fiber_phi_core/2.0, 50.0/2.0*mm );
      G4Tubs* fiber_clad_solid_x
        = new G4Tubs("fiber_clad_solid_x",
                     fiber_phi_core/2.0, fiber_phi/2.0, //Min and Max radii
                     50.*mm/2.0,//height
                     0.0*deg,
                     360.*deg);

      G4LogicalVolume* FiberXLV
        =  new G4LogicalVolume( fiber_core_solid_x, Fiber_Core_Material, "Fiber x LV");
      G4LogicalVolume* CladXLV
        =  new G4LogicalVolume( fiber_clad_solid_x, Fiber_Clad_Material, "Clad x LV");
      G4RotationMatrix* rotFiberx = new G4RotationMatrix();//Rotaion of fiber target
      rotFiberx->rotateX(90.0 * deg);//Rotation of fiber target

      // ~~~~~~~ Create logical volumes of y layers ~~~~~~~~~~
      G4EllipticalTube* fiber_core_solid_y
        = new G4EllipticalTube("fiber_core_solid_y",
                               fiber_phi_core/2.0, fiber_phi_core/2.0, 100./2.0*mm );
      G4Tubs* fiber_clad_solid_y
        = new G4Tubs("fiber_clad_solid_y",
                     fiber_phi_core/2.0, fiber_phi/2.0, //Min and Max radii
                     100.*mm/2.0,//height
                     0.0*deg,
                     360.*deg);

      G4LogicalVolume* FiberYLV
        =  new G4LogicalVolume( fiber_core_solid_y, Fiber_Core_Material, "Fiber y LV");
      G4LogicalVolume* CladYLV
        =  new G4LogicalVolume( fiber_clad_solid_y, Fiber_Clad_Material, "Clad y LV");
      G4RotationMatrix* rotFibery = new G4RotationMatrix();//Rotaion of fiber target
      rotFibery->rotateY(90.0 * deg);//Rotation of fiber target

      //double fiber_x = 150.0*mm;
      //double fiber_y =  50.0*mm;
      //double fiber_z = 100.0*mm;

      int fiber_tot_layers = 9; // 1-layer = xx'yy' // default
      double fiber_posz   = 0.0;
      double fiber_posz_p = 0.0;

      int fiber_nx = 32;
      double fiber_startx = -1.0 * fiber_phi*((double)fiber_nx-1.0) / 2.0;

      int fiber_ny = 16;
      double fiber_starty = -1.0 * fiber_phi*((double)fiber_ny-1.0) / 2.0;

      int count_fiber = 0;
      int countx=0;//, countxp=0;
      int county=0;//, countyp=0;
      double target_shift = 2.0 * (fiber_phi+fiber_phi/2.0*sqrt(3.0)) * fiber_tot_layers;
      target_shift = target_shift / 2.0 ;


      //if (TargetMaterial != "Vacuum"){
      for(int layer=0 ; layer<fiber_tot_layers ; layer++){
        // ~~~~ fiber z position ~~~~
        if(layer!=0) fiber_posz   = fiber_posz_p + fiber_phi; // xy zpos
        else fiber_posz=fiber_phi/2.0;
        fiber_posz_p = fiber_posz   + fiber_phi/2.0*sqrt(3.0);// x'y'-zpos
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for(int i=0 ; i<fiber_nx ; i++){
          countx++;
          new G4PVPlacement( rotFiberx
                             , G4ThreeVector(fiber_startx+fiber_phi*i, 0*cm, fiber_posz-target_shift)
                             ,FiberXLV
                             ,"FiberX"
                             , logTarget
                             , false
                             , countx-1 );
          //, 1000 );
          new G4PVPlacement( rotFiberx
                             //, G4ThreeVector(fiber_startx+fiber_phi*i, 0*cm, fiber_posz)
                             , G4ThreeVector(fiber_startx+fiber_phi*i, 0*cm, fiber_posz-target_shift)
                             , CladXLV
                             ,"CladX"
                             , logTarget
                             , false
                             , countx-1 );

        }
        for(int i=0 ; i<fiber_nx ; i++){
          // -- x' --
          countx++;
          new G4PVPlacement( rotFiberx
                             , G4ThreeVector(fiber_phi/2.0 + fiber_startx + fiber_phi*i,
                                             0*cm,
                                             fiber_posz_p-target_shift)
                             ,FiberXLV
                             ,"FiberX"
                             , logTarget
                             , false
                             , countx-1 );
          //, 1001 );
          new G4PVPlacement( rotFiberx
                             , G4ThreeVector(fiber_phi/2.0 + fiber_startx + fiber_phi*i,
                                             0*cm,
                                             fiber_posz_p-target_shift)
                             ,CladXLV
                             ,"CladX"
                             , logTarget
                             , false
                             , countx-1 );

        }

        // ~~~~ fiber z position ~~~~
        fiber_posz   = fiber_posz_p + fiber_phi;
        fiber_posz_p = fiber_posz + fiber_phi/2.0*sqrt(3.0);

        for(int i=0 ; i<fiber_ny ; i++){
          county++;
          // -- y --
          new G4PVPlacement( rotFibery
                             , G4ThreeVector( 0*cm,fiber_starty+fiber_phi*i, fiber_posz-target_shift)
                             ,FiberYLV
                             ,"FiberY"
                             , logTarget
                             , false
                             , county-1 );
          new G4PVPlacement( rotFibery
                             , G4ThreeVector( 0*cm,fiber_starty+fiber_phi*i, fiber_posz-target_shift)
                             ,CladYLV
                             ,"CladY"
                             , logTarget
                             , false
                             , county-1 );
        }
        // -- y' ------------------------------------
        for(int i=0 ; i<fiber_ny ; i++){
          //if(county!=0) county++;
          county++;
          new G4PVPlacement( rotFibery
                             , G4ThreeVector( 0*cm,fiber_phi/2.0+fiber_starty+fiber_phi*i,
                                              fiber_posz_p-target_shift)
                             ,FiberYLV
                             ,"FiberY"
                             , logTarget
                             , false
                             , county-1 );
          new G4PVPlacement( rotFibery
                             , G4ThreeVector( 0*cm,fiber_phi/2.0+fiber_starty+fiber_phi*i,
                                              fiber_posz_p-target_shift)
                             ,CladYLV
                             ,"CladY"
                             , logTarget
                             , false
                             , county-1 );
        }
      }
    }
    break;
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructTargetE10()
{
  // E10-2nd solid target: a plain box of the configured TargetMaterial
  // (9Be for the (K-,pi+) hypernuclear run) placed at the "Target"
  // position, which for E10 sits at the SAME location as the E90
  // target -- ~142 mm downstream of the HypTPC centre (the target
  // holder is NOT at the TPC centre; see param/DCGEO comments). The
  // box half-size and position come from the DSIZE/DCGEO "Target"
  // keys. Deliberately kept separate from the E90 LD2-tube builder
  // (ConstructTargetE90) so the two setups share everything EXCEPT the
  // target volume itself.
  const auto& half_size = sizeMan.GetSize("Target")*mm/2.;
  const auto& pos = geomMan.GetGlobalPosition("Target");
  G4Material* TargetMater = nullptr;
  const auto Target = confMan.Get<G4String>("TargetMaterial");
  if(Target == "Be"){        TargetMater = mlist.at("Be9"); }
  else if(Target == "natLi"){ TargetMater = mlist.at("natLi"); }
  else if(Target == "C"){     TargetMater = mlist.at("C"); }
  else if(Target == "CH2"){   TargetMater = mlist.Polyethylene; }
  else if(Target == "Air"){   TargetMater = mat("Air"); }
  else{
    G4cout << "[ConstructTargetE10] Target material '" << Target
           << "' not handled; using Air." << G4endl;
    TargetMater = mlist.at("Air");
  }

  G4RotationMatrix rotTarget;
  auto TargetBox = new G4Box("TargetBox",
                             half_size.x(), half_size.y(), half_size.z());
  auto logTarget = new G4LogicalVolume(TargetBox, TargetMater, "logTarget");
  new G4PVPlacement(G4Transform3D(rotTarget, pos), "physTarget",
                    logTarget, physWorld, false, 0, m_check_overlaps);
  logTarget->SetVisAttributes(G4Color::Gray());
  G4cout << "[ConstructTargetE10] " << Target << " box, half-size ("
         << half_size.x()/mm << ", " << half_size.y()/mm << ", "
         << half_size.z()/mm << ") mm at z=" << pos.z()/mm << " mm"
         << G4endl;
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructTgtHeBag()
{
  if(!use_Tgthebag) return;
  const auto& radius = 250./2.*mm;
  const auto& thickness = 900./2.*mm;
  G4Material *TgtHeBagMater = mlist.at("HeGas");


//  auto TgtHeBagBox = new G4Box
//    ("TgtHeBagBox", radius,radius, thickness);
  auto TgtHeBagBox = new G4Tubs("TgtHeBagBox",
                                0*mm, radius, thickness,
                                0*deg, 360*deg);
  auto logTgtHeBag = new G4LogicalVolume(TgtHeBagBox, TgtHeBagMater, "logTgtHeBag");
  G4RotationMatrix rotTgtHeBag;

  auto pos = geomMan.GetGlobalPosition("BAC1");
  new G4PVPlacement(G4Transform3D(rotTgtHeBag, pos),
                    "physTgtHeBag",
                    logTgtHeBag,
                    physWorld,
                    false,
                    0,
                    m_check_overlaps);

  const auto mylar_thickness = sizeMan.Get("HeBagMylarThickness")*mm;
  G4VSolid* solidTgtHeMylar;
  solidTgtHeMylar = new G4Tubs("solidTgtHeMylar",0.*mm,radius, mylar_thickness, 0*deg,360*deg);
  solidTgtHeMylar = new G4IntersectionSolid("solidTgtHeMylar", solidTgtHeMylar,  TgtHeBagBox);
  auto logicTgtHeMylar = new G4LogicalVolume(solidTgtHeMylar, mlist.at("Mylar"), "logicTgtHeMylar");
  G4RotationMatrix rot;
  new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(0, 0, -thickness+mylar_thickness/2)),
                    logicTgtHeMylar, "physTgtHeMylar", logTgtHeBag, false, 0, m_check_overlaps);
  new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(0, 0,  thickness-mylar_thickness/2)),
                    logicTgtHeMylar, "physTgtHeMylar", logTgtHeBag, false, 1, m_check_overlaps);
  logTgtHeBag->SetVisAttributes(G4Color::White());
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructQ1()
{
  const G4double a0 = 155*mm; // phi=310mm
  const G4ThreeVector pos(0, 0, -3885.35*mm);
  ///// Pole
  G4VSolid* solidQ1Gap;
  auto solidQ1Box = new G4Box("solidQ1Gap", 540*mm/2, 540*mm/2, 1240*mm/2);
  const G4double rCorner = (540*mm*std::sqrt(2)-2*a0)/2.;
  auto solidCorner = new G4Tubs("solidCorner",
                                0*mm, rCorner, 900*mm,
                                0*deg, 360*deg);
  const G4double rInner[] = { 0., 0. };
  G4double zPlaneCorner[] = { -630*mm, 630*mm };
  G4double rOuterCorner[] = { rCorner, rCorner };
  // auto solidCorner = new G4Polyhedra("solidCorner", 0*deg, 360*deg, 4, 2,
  //                                    zPlaneCorner, rInner, rOuterCorner);
  solidQ1Gap = new G4SubtractionSolid("solidQ1Gap", solidQ1Box,
                                      solidCorner, nullptr, G4ThreeVector(270*mm, 270*mm, 0));
  solidQ1Gap = new G4SubtractionSolid("solidQ1Gap", solidQ1Gap,
                                      solidCorner, nullptr, G4ThreeVector(270*mm, -270*mm, 0));
  solidQ1Gap = new G4SubtractionSolid("solidQ1Gap", solidQ1Gap,
                                      solidCorner, nullptr, G4ThreeVector(-270*mm, 270*mm, 0));
  solidQ1Gap = new G4SubtractionSolid("solidQ1Gap", solidQ1Gap,
                                      solidCorner, nullptr, G4ThreeVector(-270*mm, -270*mm, 0));
  G4RotationMatrix rotQ1Box;
  rotQ1Box.rotateZ(45*deg);
  G4double Q1Box3V = 900+1000/sqrt(2);
  G4VSolid* solidQ1Pole = nullptr;
  solidQ1Pole = new G4Box("solidQ1Pole", 2400*mm/2, 2400*mm/2, 880*mm/2);
  auto solidQ1Corner = new G4Box("tmp", 2000*mm/2, 2000*mm/2, 2000*mm/2);
  solidQ1Pole = new G4SubtractionSolid("tmp", solidQ1Pole, solidQ1Corner,
                                       G4Transform3D(rotQ1Box, G4ThreeVector(Q1Box3V,Q1Box3V,0)));
  solidQ1Pole = new G4SubtractionSolid("tmp", solidQ1Pole, solidQ1Corner,
                                       G4Transform3D(rotQ1Box, G4ThreeVector(-Q1Box3V,Q1Box3V,0)));
  solidQ1Pole = new G4SubtractionSolid("tmp", solidQ1Pole, solidQ1Corner,
                                       G4Transform3D(rotQ1Box, G4ThreeVector(Q1Box3V,-Q1Box3V,0)));
  solidQ1Pole = new G4SubtractionSolid("tmp", solidQ1Pole, solidQ1Corner,
                                       G4Transform3D(rotQ1Box, G4ThreeVector(-Q1Box3V,-Q1Box3V,0)));
  solidQ1Pole = new G4SubtractionSolid("solidQ1Pole", solidQ1Pole,
                                       solidQ1Gap, nullptr, G4ThreeVector());
  auto logicQ1Pole = new G4LogicalVolume
    (solidQ1Pole, mlist.at("Fe"), "logicQ1Pole");
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos),
                    "physQ1Pole", logicQ1Pole, physWorld, false, 0, m_check_overlaps);
  logicQ1Pole->SetVisAttributes(G4Colour::Cyan());
  ///// Coil
  G4VSolid* solidQ1Coil = nullptr;
  solidQ1Coil = new G4Box("solidQ1Coil", 1400*mm/2, 1400*mm/2, 180*mm/2);
  solidQ1Coil = new G4SubtractionSolid("solidQ1Coil", solidQ1Coil, solidQ1Box,
                                       nullptr, G4ThreeVector());
  auto logicQ1Coil = new G4LogicalVolume
    (solidQ1Coil, mlist.Cu, "logicQ1Coil");
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos - G4ThreeVector(0, 0, 530*mm)),
                    "physQ1CoilU", logicQ1Coil, physWorld, false, 0, m_check_overlaps);
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos + G4ThreeVector(0, 0, 530*mm)),
                    "physQ1CoilD", logicQ1Coil, physWorld, false, 1, m_check_overlaps);
  logicQ1Coil->SetVisAttributes(G4Color::Brown());

  ///// HeBag
  if(!use_Q1hebag) return;
  auto logicQ1HeBag = new G4LogicalVolume(solidQ1Gap, mlist.at("HeGas"), "logicQ1HeBag");
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos),
                    "physQ1HeBag", logicQ1HeBag, physWorld, false, 0, m_check_overlaps);
  const auto mylar_thickness = sizeMan.Get("HeBagMylarThickness")*mm;
  G4VSolid* solidQ1Mylar;
  solidQ1Mylar = new G4Box("solidQ1Mylar", 1*m/2, 1*m/2, mylar_thickness/2);
  solidQ1Mylar = new G4IntersectionSolid("solidQ1Mylar", solidQ1Mylar, solidQ1Gap);
  auto logicQ1Mylar = new G4LogicalVolume(solidQ1Mylar, mlist.at("Mylar"), "logicQ1Mylar");
  G4RotationMatrix rot;
  new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(0, 0, -620*mm+mylar_thickness/2)),
                    logicQ1Mylar, "physQ1Mylar", logicQ1HeBag, false, 0, m_check_overlaps);
  new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(0, 0,  620*mm-mylar_thickness/2)),
                    logicQ1Mylar, "physQ1Mylar", logicQ1HeBag, false, 1, m_check_overlaps);
  ///// Virtual Plane
  auto sdVP = G4SDManager::GetSDMpointer()->FindSensitiveDetector("VP", false);
  if(!sdVP){
    sdVP = new VPSD("VP");
    G4SDManager::GetSDMpointer()->AddNewDetector(sdVP);
  }
  {
    G4VSolid* solidVP;
    solidVP = new G4Box("solidVP1", 1*m/2, 1*m/2, 1*CLHEP::um/2);
    solidVP = new G4IntersectionSolid("solidVP1", solidVP, solidQ1Gap);
    auto logicVP = new G4LogicalVolume(solidVP, mlist.at("HeGas"), "logicVP1");
    logicVP->SetSensitiveDetector(sdVP);
    logicVP->SetVisAttributes(G4Color::Red());
    // logicVP->SetVisAttributes(G4VisAttributes::GetInvisible());
    auto rotVP1 = new G4RotationMatrix;
    rotVP1->rotateY(-geomMan.GetRotAngle2("VP1")*deg);
    // auto pos = geomMan.GetGlobalPosition("VP1");
    new G4PVPlacement(rotVP1, G4ThreeVector(), logicVP, "physVP1",
                      logicQ1HeBag, false, 1, m_check_overlaps);
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructQ2()
{
  const G4double a0 = 360*mm/2; // phi=360mm
  const G4ThreeVector pos(0, 0, -2776.5*mm);
  ///// Pole
  G4VSolid* solidQ2Pole = nullptr;
  solidQ2Pole = new G4Box("solidQ2Pole", 2100*mm/2, 1540*mm/2, 540*mm/2);
  const G4double zlength = 970*mm/2;

  // auto solidQ2Gap = new G4Box("solidQ2Gap", a0, a0, zlength);
  G4VSolid* solidQ2Gap;
  solidQ2Gap = new G4Box("solidQ2Gap", 1186*mm/2, 606*mm/2, zlength);
  auto solidCorner = new G4EllipticalTube("solidCorner",
                                          383*mm, 210*mm, zlength+1*mm);
  G4TwoVector cornerPos(450*mm, 245*mm);
  solidQ2Gap = new G4SubtractionSolid("solidQ1Gap", solidQ2Gap,
                                      solidCorner, nullptr,
                                      G4ThreeVector(cornerPos.x(), cornerPos.y(), 0));
  solidQ2Gap = new G4SubtractionSolid("solidQ1Gap", solidQ2Gap,
                                      solidCorner, nullptr,
                                      G4ThreeVector(-cornerPos.x(), cornerPos.y(), 0));
  solidQ2Gap = new G4SubtractionSolid("solidQ1Gap", solidQ2Gap,
                                      solidCorner, nullptr,
                                      G4ThreeVector(cornerPos.x(), -cornerPos.y(), 0));
  solidQ2Gap = new G4SubtractionSolid("solidQ1Gap", solidQ2Gap,
                                      solidCorner, nullptr,
                                      G4ThreeVector(-cornerPos.x(), -cornerPos.y(), 0));
  G4RotationMatrix rot;
  // rot.rotateZ(45*deg);
  solidQ2Pole = new G4SubtractionSolid("solidQ2Pole", solidQ2Pole, solidQ2Gap,
                                       G4Transform3D(rot, G4ThreeVector()));
  auto logicQ2Pole = new G4LogicalVolume(solidQ2Pole, mlist.at("Fe"), "logicQ2Pole");
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos),
                    "physQ2Pole", logicQ2Pole, physWorld, false, 0, m_check_overlaps);
  logicQ2Pole->SetVisAttributes(G4Color::Cyan());
  ///// Coil
  G4VSolid* solidQ2Coil;
  solidQ2Coil = new G4Box("solidQ2Coil", 1320*mm/2, 880*mm/2, 130*mm/2);
  auto solidQ2Box = new G4Box("solidQ2Gap", 1200*mm/2, 610*mm/2, 900*mm/2);
  solidQ2Coil = new G4SubtractionSolid("solidQ2Coil", solidQ2Coil, solidQ2Box,
                                       nullptr, G4ThreeVector());
  auto logicQ2Coil = new G4LogicalVolume
    (solidQ2Coil, mlist.Cu, "logicQ2Coil");
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos - G4ThreeVector(0, 0, 335*mm)),
                    "physQ2CoilU", logicQ2Coil, physWorld, false, 0, m_check_overlaps);
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos + G4ThreeVector(0, 0, 335*mm)),
                    "physQ2CoilD", logicQ2Coil, physWorld, false, 1, m_check_overlaps);
  logicQ2Coil->SetVisAttributes(G4Color::Brown());
  ///// HeBag
  if(!use_Q2hebag) return;
  auto logicQ2HeBag = new G4LogicalVolume(solidQ2Gap, mlist.at("HeGas"), "logicQ2HeBag");
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(rot), pos),
                    "physQ2HeBag", logicQ2HeBag, physWorld, false, 0, m_check_overlaps);
  const auto mylar_thickness = sizeMan.Get("HeBagMylarThickness")*mm;
  G4VSolid* solidQ2Mylar;
  solidQ2Mylar = new G4Box("solidQ2Mylar", 2*m/2, 2*m/2, mylar_thickness/2);
  solidQ2Mylar = new G4IntersectionSolid("solidQ2Mylar", solidQ2Mylar, solidQ2Gap);
  auto logicQ2Mylar = new G4LogicalVolume(solidQ2Mylar, mlist.at("Mylar"), "logicQ2Mylar");
  new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(0, 0, -zlength+mylar_thickness/2)),
                    logicQ2Mylar, "physQ2Mylar", logicQ2HeBag, false, 0, m_check_overlaps);
  new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(0, 0,  zlength-mylar_thickness/2)),
                    logicQ2Mylar, "physQ2Mylar", logicQ2HeBag, false, 1, m_check_overlaps);
  ///// Virtual Plane
  auto sdVP = G4SDManager::GetSDMpointer()->FindSensitiveDetector("VP", false);
  if(!sdVP){
    sdVP = new VPSD("VP");
    G4SDManager::GetSDMpointer()->AddNewDetector(sdVP);
  }
  {
    auto solidTmp = new G4Box("solidTmp", 2*m/2, 2*m/2, 1*CLHEP::um/2);
    auto solidVP = new G4IntersectionSolid("solidVP2", solidTmp, solidQ2Gap);
    auto logicVP = new G4LogicalVolume(solidVP, mlist.at("HeGas"), "logicVP2");
    logicVP->SetSensitiveDetector(sdVP);
    logicVP->SetVisAttributes(G4Color::Red());
    // logicVP->SetVisAttributes(G4VisAttributes::GetInvisible());
    auto rotVP2 = new G4RotationMatrix;
    rotVP2->rotateY(-geomMan.GetRotAngle2("VP2")*deg);
    new G4PVPlacement(rotVP2, G4ThreeVector(0, 0, -469.7), logicVP, "physVP2",
                      logicQ2HeBag, false, 2, m_check_overlaps);
    auto rotVP3 = new G4RotationMatrix;
    rotVP3->rotateY(-geomMan.GetRotAngle2("VP3")*deg);
    new G4PVPlacement(rotVP3, G4ThreeVector(), logicVP, "physVP3",
                      logicQ2HeBag, false, 3, m_check_overlaps);
    auto rotVP4 = new G4RotationMatrix;
    rotVP4->rotateY(-geomMan.GetRotAngle2("VP4")*deg);
    new G4PVPlacement(rotVP4, G4ThreeVector(0, 0, 472.95), logicVP, "physVP4",
                      logicQ2HeBag, false, 4, m_check_overlaps);
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructD1()
{
  G4Material *mD1Gap = use_D1hebag ? mlist.at("HeGas") : mlist.at("Air");
  // D magnet surface
  G4Tubs *solD1Tub = new G4Tubs("solD1Tub", 0, Dfr2, DfHalfGap,
                                0.*degree, bendAngleDf*degree);
  // D magnet
  G4Tubs *solD1Gap = new G4Tubs("solD1Gap", Dr1, Dr2, DHalfGap,
                                0.*degree, bendAngleD*degree);
  G4Box *tmpD1 = new G4Box("tmp", 1000, 1000, 1000);
  G4RotationMatrix tmpD1rot;
  tmpD1rot.rotateZ(35*deg);
  G4double tmpD1r = 1380-1000;
  G4ThreeVector tmpD13Y(tmpD1r*std::cos(35*deg),tmpD1r*std::sin(35*deg),0);
  ///// Yoke
  auto tmpsolD1Pole = new G4SubtractionSolid("tmp", solD1Tub, solD1Gap);
  auto solD1Pole = new G4SubtractionSolid
    ("solidDy", tmpsolD1Pole, tmpD1, G4Transform3D(tmpD1rot, tmpD13Y));
  auto logicD1Gap = new G4LogicalVolume(solD1Gap, mD1Gap, "logicD1Gap");
  auto logicD1Pole = new G4LogicalVolume(solD1Pole, mlist.at("Fe"), "logicD1Pole");
  G4RotationMatrix rotD1;
  rotD1.rotateX(90.*deg);
  rotD1.rotateY(-110.*deg);
  G4ThreeVector pos(3*m, 0, -3*m*std::tan(35*deg));
  new G4PVPlacement(G4Transform3D(rotD1, pos),
                    "physD1Pole", logicD1Pole, physWorld, false, 0, m_check_overlaps);
  new G4PVPlacement(G4Transform3D(rotD1, pos),
                    "physD1Gap", logicD1Gap, physWorld, false, 0, m_check_overlaps);
  ///// Endguard
  G4Box *solidEG1 = new G4Box("solidEG1", 1880*mm/2., 1600*mm/2., 76*mm/2.);
  G4Box *solidEG2 = new G4Box("solidEG2", 810*mm/2., 330*mm/2., 80*mm/2.);
  auto solidD1EG = new G4SubtractionSolid("solEnd", solidEG1, solidEG2);
  auto logicD1EG = new G4LogicalVolume(solidD1EG, mlist.at("Fe"), "logicD1EG");
  G4ThreeVector posEG(3*m*std::tan(35*deg) + 278.5*mm);
  posEG.rotateY(-20.*deg);
  G4RotationMatrix rotEG;
  rotEG.rotateY(70.*deg);
  new G4PVPlacement(G4Transform3D(rotEG, posEG),
        	    "physD1EG", logicD1EG, physWorld, false, 0, m_check_overlaps);
  logicD1Pole->SetVisAttributes(G4Color::Blue());
  // logicD1Gap->SetVisAttributes(G4Color::Blue());
  logicD1EG->SetVisAttributes(G4Color::Blue());
  logicD1Gap->SetUserLimits(new G4UserLimits(3.*mm));
  ///// Coil
  G4VSolid* solidD1Coil = nullptr;
  solidD1Coil = new G4Box("solidD1Coil", 1600*mm/2, 1600*mm/2, 165*mm/2);
  auto solidD1Box = new G4Box("solidD1Gap", 800*mm/2, 500*mm/2, 900*mm/2);
  solidD1Coil = new G4SubtractionSolid("solidD1Coil", solidD1Coil, solidD1Box,
                                       nullptr, G4ThreeVector());
  auto logicD1Coil = new G4LogicalVolume
    (solidD1Coil, mlist.Cu, "logicD1Coil");
  G4ThreeVector posCoil(0, 0, -3*m*std::tan(35*deg)-165*mm/2);
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), posCoil),
                    "physD1CoilU", logicD1Coil, physWorld, false, 0, m_check_overlaps);
  posCoil.rotateY(-110*deg);
  G4RotationMatrix rotCoil;
  rotCoil.rotateY(-110*deg);
  new G4PVPlacement(G4Transform3D(rotCoil, posCoil),
                    "physD1CoilD", logicD1Coil, physWorld, false, 1, m_check_overlaps);
  logicD1Coil->SetVisAttributes(G4Color::Brown());

  ///// HeBag UpStream
  if(!use_D1hebag){ return;}
  else{
    const G4double zlengthU = 180*mm;
    auto solidHeBagU = new G4Box("solidD1HeBag", (Dr2-Dr1)/2, DHalfGap, zlengthU/2);
    auto logicD1HeBagU = new G4LogicalVolume(solidHeBagU, mlist.at("HeGas"), "logicD1HeBag");
    G4ThreeVector posD1HeBag(0, 0, -3*m*std::tan(35*deg));
    posD1HeBag.setX(0);
    posD1HeBag += G4ThreeVector(0, 0, -zlengthU/2);
    G4RotationMatrix rot;
    new G4PVPlacement(G4Transform3D(rot, posD1HeBag),
                      "physD1HeBagU", logicD1HeBagU, physWorld, false, 0, m_check_overlaps);
    const auto mylar_thickness = sizeMan.Get("HeBagMylarThickness")*mm;
    auto solidD1MylarU = new G4Box("solidD1MylarU", (Dr2-Dr1)/2, DHalfGap, mylar_thickness/2);
    auto logicD1MylarU = new G4LogicalVolume(solidD1MylarU, mlist.at("Mylar"), "logicD1MylarU");
    new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(0, 0, -zlengthU/2+mylar_thickness/2)),
                      logicD1MylarU, "physD1Mylar", logicD1HeBagU, false, 0, m_check_overlaps);

    ///// HeBag Downstream
    const G4double zlengthD = 320*mm;
    auto solidHeBagD = new G4Box("solidD1HeBag", (Dr2-Dr1)/2, DHalfGap, zlengthD/2);
    auto logicD1HeBagD = new G4LogicalVolume(solidHeBagD, mlist.at("HeGas"), "logicD1HeBag");
    posD1HeBag.set(3*m*std::tan(35*deg) + zlengthD/2, 0, 0);
    posD1HeBag.rotateY(-20.*deg);
    rot.rotateY(70.*deg);
    new G4PVPlacement(G4Transform3D(rot, posD1HeBag),
                      "physD1HeBagD", logicD1HeBagD, physWorld, false, 0, m_check_overlaps);
    auto solidD1MylarD = new G4Box("solidD1MylarD", (Dr2-Dr1)/2, DHalfGap, mylar_thickness/2);
    auto logicD1MylarD = new G4LogicalVolume(solidD1MylarD, mlist.at("Mylar"), "logicD1MylarD");
    rot.rotateY(-70.*deg);
    new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(0, 0, zlengthU/2-mylar_thickness/2)),
                      logicD1MylarD, "physD1MylarD", logicD1HeBagD, false, 1, m_check_overlaps);

    ///// Virtual Plane
    auto sdVP = G4SDManager::GetSDMpointer()->FindSensitiveDetector("VP", false);
    if(!sdVP){
      sdVP = new VPSD("VP");
      G4SDManager::GetSDMpointer()->AddNewDetector(sdVP);
    }
    {
      auto solidVP = new G4Box("solidTmp", 800*mm/2, 320*mm/2, 1*CLHEP::um/2);
      auto logicVP = new G4LogicalVolume(solidVP, mlist.at("HeGas"), "logicVP5");
      logicVP->SetSensitiveDetector(sdVP);
      logicVP->SetVisAttributes(G4Color::Red());
      // logicVP->SetVisAttributes(G4VisAttributes::GetInvisible());
      auto rotVP5 = new G4RotationMatrix;
      rotVP5->rotateY(-geomMan.GetRotAngle2("VP5")*deg);
      new G4PVPlacement(rotVP5, G4ThreeVector(0, 0, zlengthU/2-1*CLHEP::um/2), logicVP,
                        "physVP5", logicD1HeBagU, false, 5, m_check_overlaps);
      for(G4int i=6; i<=8; ++i){
        auto rotVP = new G4RotationMatrix;
        rotVP->rotateX(90*deg);
        rotVP->rotateY(70*deg);
        const auto& ra2VP = geomMan.GetRotAngle2("VP"+std::to_string(i))*deg;
        rotVP->rotateY(-ra2VP);
        G4ThreeVector posVP(rhoD, 0, 0);
        posVP.rotateZ(70*deg-ra2VP);
        new G4PVPlacement(rotVP, posVP, logicVP, "physVP"+std::to_string(i),
                          logicD1Gap, false, i, m_check_overlaps);
      }
      new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -zlengthD/2+1*CLHEP::um/2),
                        logicVP, "physVP9", logicD1HeBagD, false, 9, m_check_overlaps);
      new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zlengthD/2-1*CLHEP::um/2),
                      logicVP, "physVP10", logicD1HeBagD, false, 10, m_check_overlaps);
    }
  }

  ///// HeBag between SDC34
  if(!use_SDCouthebag){ return;}
  else{
    const G4double thickness = 600*mm;
    const G4double width = 2000*mm;
    const G4double height = 1600*mm;

    auto solidSDCoutHeBag = new G4Box("solidSDCoutHeBag", width/2, height/2, thickness/2);
    auto logicSDCoutHeBag = new G4LogicalVolume(solidSDCoutHeBag, mlist.at("HeGas"), "logicSDCoutHeBag");
    const auto& posSDCoutHeBag = (geomMan.GetGlobalPosition("SDC3-Y2") +
                       geomMan.GetGlobalPosition("SDC4-Y1"))/2;
    G4RotationMatrix rot;
    rot.rotateY(70.*deg);
    new G4PVPlacement(G4Transform3D(rot, posSDCoutHeBag),
                      "physSDCoutHeBag", logicSDCoutHeBag, physWorld, false, 0, m_check_overlaps);

    const auto mylar_thickness = sizeMan.Get("HeBagMylarThickness")*mm;
    G4VSolid* solidSDCoutHeMylar;
    solidSDCoutHeMylar = new G4Box("solidSDCoutHeMylar", width/2, height/2, mylar_thickness/2);
    solidSDCoutHeMylar = new G4IntersectionSolid("solidSDCoutHeMylar", solidSDCoutHeMylar, solidSDCoutHeBag);
    auto logicSDCoutHeMylar = new G4LogicalVolume(solidSDCoutHeMylar, mlist.at("Mylar"), "logicSDCoutHeMylar");
    rot.rotateY(-70.*deg);
    new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(0, 0, -thickness/2+mylar_thickness/2)),
                      logicSDCoutHeMylar, "physSDCoutHeMylar", logicSDCoutHeBag, false, 0, m_check_overlaps);

    new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(0, 0, thickness/2-mylar_thickness/2)),
                      logicSDCoutHeMylar, "physSDCoutHeMylar", logicSDCoutHeBag, false, 1, m_check_overlaps);
  }

}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructSDC1()
{
  auto sdSDC1 = new DCSD("SDC1");
  G4SDManager::GetSDMpointer()->AddNewDetector(sdSDC1);
  const auto& sdc1_pos = (geomMan.GetGlobalPosition("SDC1-X1") +
                          geomMan.GetGlobalPosition("SDC1-X2"))/2;
  const auto& frame_size = sizeMan.GetSize("Sdc1Frame")*mm/2;
  const auto& drift_size = sizeMan.GetSize("Sdc1Drift")*mm/2;
  auto sdc1_solid = new G4Box("Sdc1Solid", frame_size.x(),
                              frame_size.y(), frame_size.z());
  auto sdc1_logic = new G4LogicalVolume(sdc1_solid, mlist.at("Ar80IsoButane20Gas"),
                                        "Sdc1LOGIC", 0, 0, 0);
  sdc1_logic->SetVisAttributes(G4Colour::Green());
  new G4PVPlacement(nullptr, sdc1_pos,
                    "Sdc1PHYS", sdc1_logic, physWorld, false, 0, m_check_overlaps);
  auto sdc1pl_solid = new G4Box("Sdc1PlSolid", drift_size.x(),
                                drift_size.y(), drift_size.z());
  G4String plane_name[] = { "SDC1-U1", "SDC1-U2",
			    "SDC1-X1", "SDC1-X2",
			    "SDC1-V1", "SDC1-V2" };
  for(G4int i=0; i<NumOfLayersSDC1; ++i){
    G4ThreeVector pos = geomMan.GetGlobalPosition(plane_name[i]);
    pos -= sdc1_pos;
    auto sdc1pl_logic = new G4LogicalVolume(sdc1pl_solid,
                                            mlist.at("Ar80IsoButane20Gas"),
                                            plane_name[i] + "LOGIC", 0, 0, 0);
    sdc1pl_logic->SetSensitiveDetector(sdSDC1);
    new G4PVPlacement(nullptr, pos, sdc1pl_logic, plane_name[i] + "PHYS",
                      sdc1_logic, false, 101+i, m_check_overlaps);
  }
  ///// Mylar
  const G4double almylar_thickness = sizeMan.Get("AlMylarThickness")*mm;
  const G4double alplate_thickness = sizeMan.Get("AlPlateThickness")*mm;
  const G4double mylar_thickness = sizeMan.Get("Sdc1MylarThickness")*mm;
  auto solidMylar = new G4Box("solidMylar", frame_size.x(),
                              frame_size.y(), mylar_thickness/2);
  auto logicMylar = new G4LogicalVolume
    (solidMylar, mlist.at("Mylar"), "logicMylar");
  G4ThreeVector pos = sdc1_pos;
  pos.setZ(sdc1_pos.z()-frame_size.z()-mylar_thickness/2);
  new G4PVPlacement(nullptr, pos, "physMylar", logicMylar,
                    physWorld, false, 0, m_check_overlaps);
  pos.setZ(sdc1_pos.z()+frame_size.z()+mylar_thickness/2);
  new G4PVPlacement(nullptr, pos, "physMylar", logicMylar,
                    physWorld, false, 1, m_check_overlaps);
  ///// Al mylar
  auto solidAlMylar = new G4Box("solidAlMylar", frame_size.x(),
                                frame_size.y(), almylar_thickness/2);
  auto logicAlMylar = new G4LogicalVolume
    (solidAlMylar, mlist.at("Mylar"), "logicAlMylar");
  pos = sdc1_pos;
  pos.setZ(sdc1_pos.z()-frame_size.z()-mylar_thickness-almylar_thickness/2);
  new G4PVPlacement(nullptr, pos, "physAlMylar", logicAlMylar,
                    physWorld, false, 0, m_check_overlaps);
  pos.setZ(sdc1_pos.z()+frame_size.z()+mylar_thickness+almylar_thickness/2);
  new G4PVPlacement(nullptr, pos, "physAlMylar", logicAlMylar,
                    physWorld, false, 1, m_check_overlaps);
  ///// Al plate
  auto solidAlPlate = new G4Box("solidAlMylar", frame_size.x(),
                                frame_size.y(), alplate_thickness/2);
  auto logicAlPlate = new G4LogicalVolume
    (solidAlPlate, mlist.at("Al"), "logicAlPlate");
  pos = sdc1_pos;
  pos.setZ(sdc1_pos.z()-frame_size.z()-mylar_thickness-almylar_thickness-alplate_thickness/2);
  new G4PVPlacement(nullptr, pos, "physAlPlate", logicAlPlate,
                    physWorld, false, 0, m_check_overlaps);
  pos.setZ(sdc1_pos.z()+frame_size.z()+mylar_thickness+almylar_thickness+alplate_thickness/2);
  new G4PVPlacement(nullptr, pos, "physAlPlate", logicAlPlate,
                    physWorld, false, 1, m_check_overlaps);
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructSDC2()
{
  auto sdSDC2 = new DCSD("SDC2");
  G4SDManager::GetSDMpointer()->AddNewDetector(sdSDC2);
  const auto& sdc2_pos = (geomMan.GetGlobalPosition("SDC2-V2") +
                          geomMan.GetGlobalPosition("SDC2-U1"))/2;
  const auto& frame_size = sizeMan.GetSize("Sdc2Frame")*mm/2;
  const auto& drift_size = sizeMan.GetSize("Sdc2Drift")*mm/2;
  auto sdc2_solid = new G4Box("Sdc2Solid", frame_size.x(),
                              frame_size.y(), frame_size.z());
  auto sdc2_logic = new G4LogicalVolume(sdc2_solid, mlist.at("Ar80IsoButane20Gas"),
                                        "Sdc2LOGIC", 0, 0, 0);
  sdc2_logic->SetVisAttributes(G4Colour::Green());
  new G4PVPlacement(nullptr, sdc2_pos,
                    "Sdc2PHYS", sdc2_logic, physWorld, false, 0, m_check_overlaps);
  auto sdc2pl_solid = new G4Box("Sdc2PlSolid", drift_size.x(),
                                drift_size.y(), drift_size.z());
  G4String plane_name[] = { "SDC2-V1", "SDC2-V2",
			    "SDC2-U1", "SDC2-U2" };
  for(G4int i=0; i<NumOfLayersSDC2; ++i){
    G4ThreeVector pos = geomMan.GetGlobalPosition(plane_name[i]);
    pos -= sdc2_pos;
    auto sdc2pl_logic = new G4LogicalVolume(sdc2pl_solid,
                                            mlist.at("Ar80IsoButane20Gas"),
                                            plane_name[i] + "LOGIC", 0, 0, 0);
    sdc2pl_logic->SetSensitiveDetector(sdSDC2);
    new G4PVPlacement(nullptr, pos, sdc2pl_logic, plane_name[i] + "PHYS",
                      sdc2_logic, false, 101+i, m_check_overlaps);
  }
  ///// Mylar
  const G4double mylar_thickness = sizeMan.Get("Sdc2MylarThickness")*mm;
  const G4double alplate_thickness = sizeMan.Get("AlPlateThickness")*mm;
  const G4double cplate_thickness = sizeMan.Get("CathodePlateThickness")*mm;
  auto solidMylar = new G4Box("solidMylar", frame_size.x(),
                              frame_size.y(), mylar_thickness/2);
  auto logicMylar = new G4LogicalVolume
    (solidMylar, mlist.at("Mylar"), "logicMylar");
  for(G4int i=0; i<2; ++i){
    G4ThreeVector pos;
    if(i==0) pos.setZ(-frame_size.z()+alplate_thickness+mylar_thickness/2);
    if(i==1) pos.setZ(+frame_size.z()-alplate_thickness-mylar_thickness/2);
    new G4PVPlacement(nullptr, pos, logicMylar, "physSdc2WindowMylar",
                      sdc2_logic, false, i, m_check_overlaps);
  }
  ///// Al plate
  auto solidAlPlate = new G4Box("solidAlPlate", frame_size.x(),
                                frame_size.y(), alplate_thickness/2);
  auto logicAlPlate = new G4LogicalVolume
    (solidAlPlate, mlist.at("Al"), "logicAlPlate");
  for(G4int i=0; i<2; ++i){
    G4ThreeVector pos;
    if(i==0) pos.setZ(-frame_size.z()+alplate_thickness/2);
    if(i==1) pos.setZ(+frame_size.z()-alplate_thickness/2);
    new G4PVPlacement(nullptr, pos, logicAlPlate, "physSdc2WindowAlPlate",
                      sdc2_logic, false, i, m_check_overlaps);
    if(i==0) pos.setZ(-frame_size.z()+alplate_thickness*3/2+mylar_thickness);
    if(i==1) pos.setZ(+frame_size.z()-alplate_thickness*3/2-mylar_thickness);
    new G4PVPlacement(nullptr, pos, logicAlPlate, "physSdc2WindowAlPlate",
                      sdc2_logic, false, i, m_check_overlaps);
  }
  ///// Cathode
  auto solidCathodePlate = new G4Box("solidCathodePlate", frame_size.x(),
                                     frame_size.y(), cplate_thickness/2);
  auto logicCarbonPlate = new G4LogicalVolume
    (solidCathodePlate, mlist.at("C"), "logicCarbonPlate");
  logicCarbonPlate->SetVisAttributes(G4Color::Gray());
  const G4int NumOfCathode = 6;
  G4double zoffset[NumOfCathode] = { -14.8*mm, -10*mm, -5.2*mm,
                                     5.2*mm, 10*mm, 14.8*mm };
  for(G4int i=0; i<NumOfCathode; ++i){
    G4ThreeVector pos(0, 0, zoffset[i]);
    new G4PVPlacement(nullptr, pos, logicMylar, "physSdc2CathodeMylar",
                      sdc2_logic, false, i, m_check_overlaps);
    new G4PVPlacement(nullptr, pos-G4ThreeVector(0, 0, mylar_thickness/2+cplate_thickness/2),
                      logicCarbonPlate, "physSdc2CarbonPlate",
                      sdc2_logic, false, 2*i, m_check_overlaps);
    new G4PVPlacement(nullptr, pos+G4ThreeVector(0, 0, mylar_thickness/2+cplate_thickness/2),
                      logicCarbonPlate, "physSdc2CarbonPlate",
                      sdc2_logic, false, 2*i+1, m_check_overlaps);
  }

}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructSDC3()
{
  ConstructKLChamber(3);
  ConstructKLChamberMylar(3);
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructSDC4()
{
  ConstructKLChamber(4);
  ConstructKLChamberMylar(4);
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructSDC5()
{
  ConstructKLChamber(5);
  ConstructKLChamberMylar(5);
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructKLChamber(G4int i)
{
  const G4String name = "Sdc"+std::to_string(i);
  const G4String uname = "SDC"+std::to_string(i);
  auto sd = new DCSD("SDC"+std::to_string(i));
  G4SDManager::GetSDMpointer()->AddNewDetector(sd);
  G4int NumOfLayers;
  std::vector<G4String> layer_name;
  switch(i){
  case 3:
    NumOfLayers = NumOfLayersSDC3;
    layer_name.push_back("X1");
    layer_name.push_back("X2");
    layer_name.push_back("Y1");
    layer_name.push_back("Y2");
    break;
  case 4:
    NumOfLayers = NumOfLayersSDC4;
    layer_name.push_back("Y1");
    layer_name.push_back("Y2");
    layer_name.push_back("X1");
    layer_name.push_back("X2");
    break;
  case 5:
    NumOfLayers = NumOfLayersSDC5;
    layer_name.push_back("Y1");
    layer_name.push_back("Y2");
    layer_name.push_back("X1");
    layer_name.push_back("X2");
    break;
  default:
    throw std::out_of_range("no such key : " + name);
    break;
  }
  const auto& frame_size = sizeMan.GetSize(name+"Frame")*mm/2;
  const auto& ra2 = geomMan.GetRotAngle2(uname+"-X1")*deg;
  auto solidGas = new G4Box("solid"+name+"Gas", frame_size.x(),
                            frame_size.y(), frame_size.z());
  auto logicGas = new G4LogicalVolume(solidGas, mlist.at("Ar50Ethane50Gas"),
                                      "logic"+name+"Gas");
  logicGas->SetVisAttributes(G4Colour::Green());
  { ///// Frame
    auto rot = new G4RotationMatrix;
    rot->rotateY(-ra2);
    const auto& pos = (geomMan.GetGlobalPosition(uname+"-X1") +
		       geomMan.GetGlobalPosition(uname+"-Y2"))/2;
    new G4PVPlacement(rot, pos, "phys"+name+"Gas",
		      logicGas, physWorld, false, 0, m_check_overlaps);
  }
  const auto& sdc_local_z = (geomMan.GetLocalZ(uname+"-X1") +
			     geomMan.GetLocalZ(uname+"-Y2"))/2;
  for(G4int l=0; l<NumOfLayers; ++l){
    const auto& layer_size = sizeMan.GetSize(name+layer_name[l]+"Layer")*mm/2;
    auto solidLayer = new G4Box("solid"+name+layer_name[l], layer_size.x(),
                                layer_size.y(), layer_size.z());
    const auto& local_z = geomMan.GetLocalZ(uname+"-"+layer_name[l])*mm;
    G4double zoffset = local_z-sdc_local_z;
    G4ThreeVector pos(0, 0, zoffset);
    auto logicLayer = new G4LogicalVolume(solidLayer,
                                          mlist.at("Ar50Ethane50Gas"),
                                          "logic"+name+layer_name[l]);
    logicLayer->SetSensitiveDetector(sd);
    new G4PVPlacement(nullptr, pos, logicLayer, "phys"+name+layer_name[l],
                      logicGas, false, 101+l, m_check_overlaps);
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructKLChamberMylar(G4int i)
{
  const G4String name = "Sdc"+std::to_string(i);
  const G4String uname = "SDC"+std::to_string(i);
  const auto& frame_size = sizeMan.GetSize(name+"Frame")*mm/2;
  const auto& centerL = ((geomMan.GetGlobalPosition(uname+"-X1") +
                          geomMan.GetGlobalPosition(uname+"-Y2"))/2).mag();
  const auto& ra2 = geomMan.GetRotAngle2(uname+"-X1")*deg;
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  ///// Mylar
  const G4double space = 0.01*mm;
  const G4double almylar_thickness = sizeMan.Get("AlMylarThickness")*mm;
  const G4double alplate_thickness = sizeMan.Get("AlPlateThickness")*mm;
  const G4double mylar_thickness = sizeMan.Get(name+"MylarThickness")*mm;
  auto solidMylar = new G4Box("solidMylar", frame_size.x(),
                              frame_size.y(), mylar_thickness/2);
  auto logicMylar = new G4LogicalVolume
    (solidMylar, mlist.at("Mylar"), "logicMylar");
  G4double d = frame_size.z() + space + mylar_thickness/2;
  {
    G4ThreeVector pos(0, 0, centerL-d);
    pos.rotateY(ra2);
    new G4PVPlacement(rot, pos, "physMylar", logicMylar,
                      physWorld, false, 100*i+0, m_check_overlaps);
  }
  {
    G4ThreeVector pos(0, 0, centerL+d);
    pos.rotateY(ra2);
    new G4PVPlacement(rot, pos, "physMylar", logicMylar,
                      physWorld, false, 100*i+1, m_check_overlaps);
  }
  ///// Al mylar
  auto solidAlMylar = new G4Box("solidAlMylar", frame_size.x(),
                                frame_size.y(), almylar_thickness/2);
  auto logicAlMylar = new G4LogicalVolume
    (solidAlMylar, mlist.at("Mylar"), "logicAlMylar");
  d += mylar_thickness/2 + almylar_thickness/2;
  {
    G4ThreeVector pos(0, 0, centerL-d);
    pos.rotateY(ra2);
    new G4PVPlacement(rot, pos, "physAlMylar", logicAlMylar,
                      physWorld, false, 100*i+0, m_check_overlaps);
  }
  {
    G4ThreeVector pos(0, 0, centerL+d);
    pos.rotateY(ra2);
    new G4PVPlacement(rot, pos, "physAlMylar", logicAlMylar,
                      physWorld, false, 100*i+1, m_check_overlaps);
  }
  ///// Al plate
  auto solidAlPlate = new G4Box("solidAlMylar", frame_size.x(),
                                frame_size.y(), alplate_thickness/2);
  auto logicAlPlate = new G4LogicalVolume
    (solidAlPlate, mlist.at("Al"), "logicAlPlate");
  d += almylar_thickness/2 + alplate_thickness/2;
  {
    G4ThreeVector pos(0, 0, centerL-d);
    pos.rotateY(ra2);
    new G4PVPlacement(rot, pos, "physAlPlate", logicAlPlate,
                      physWorld, false, 100*i+0, m_check_overlaps);
  }
  {
    G4ThreeVector pos(0, 0, centerL+d);
    pos.rotateY(ra2);
    new G4PVPlacement(rot, pos, "physAlPlate", logicAlPlate,
                      physWorld, false, 100*i+1, m_check_overlaps);
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructSSD()
{
  const auto& base_size = sizeMan.GetSize("SsdBase")*mm/2;
  const auto& area_size = sizeMan.GetSize("SsdArea")*mm/2;
  G4VSolid* solidBase = nullptr;
  solidBase = new G4Box("solidSsdBase", base_size.x(),
                        base_size.y(), base_size.z());
  auto solidWindow = new G4Box("solidSsdWindow", area_size.x(),
                               area_size.y(), base_size.z()+1*mm);
  solidBase = new G4SubtractionSolid("solidSsdBase", solidBase, solidWindow);
  auto solidArea = new G4Box("solidSsdArea", area_size.x(),
                             area_size.y(), area_size.z());

  auto logicBase = new G4LogicalVolume(solidBase, mlist.at("G10"), "logicSsdBase");
  auto logicArea = new G4LogicalVolume(solidArea, mlist.at("Si"), "logicSsdArea");
  logicBase->SetVisAttributes(G4Color::Green());
  logicArea->SetVisAttributes(G4Color::Yellow());
  // SSD1StripLV->SetSensitiveDetector( ssdSD );
  { ///// SSD1
    G4String name[NumOfLayersSSD1]
      = { "SSD1-Y1", "SSD1-X1", "SSD1-Y2", "SSD1-X2" };
    {
      auto pos = (geomMan.GetGlobalPosition(name[0])+
                  geomMan.GetGlobalPosition(name[3]))*0.5;
      new G4PVPlacement(nullptr, pos, "physSsd1Base", logicBase,
                        physWorld, false, 0, m_check_overlaps);
    }
    for(G4int i=0; i<NumOfLayersSSD1; ++i){
      const auto& pos = geomMan.GetGlobalPosition(name[i]);
      new G4PVPlacement(nullptr, pos, "physSsd1Area"+std::to_string(i),
                        logicArea, physWorld, false, i, m_check_overlaps);
    }
  }
  { ///// SSD2
    G4String name[NumOfLayersSSD2]
      = { "SSD2-X1", "SSD2-Y1", "SSD2-X2", "SSD2-Y2" };
    {
      auto pos = (geomMan.GetGlobalPosition(name[0])+
                  geomMan.GetGlobalPosition(name[3]))*0.5;
      new G4PVPlacement(nullptr, pos, "physSsd2Base", logicBase,
                        physWorld, false, 0, m_check_overlaps);
    }
    for(G4int i=0; i<NumOfLayersSSD2; ++i){
      const auto& pos = geomMan.GetGlobalPosition(name[i]);
      new G4PVPlacement(nullptr, pos, "physSsd2Area"+std::to_string(i),
                        logicArea, physWorld, false, i, m_check_overlaps);
    }
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructSFT()
{
  const auto& frame_size = sizeMan.GetSize("SftFrame")*mm/2;
  const auto& xarea_size = sizeMan.GetSize("SftXArea")*mm/2;
  const auto& uvarea_size = sizeMan.GetSize("SftUVArea")*mm/2;
  G4VSolid* solidFrame = nullptr;
  solidFrame = new G4Box("solidSftFrame", frame_size.x(),
                         frame_size.y(), frame_size.z());
  auto solidWindow = new G4Box("solidSftWindow", uvarea_size.x(),
                               uvarea_size.y(), frame_size.z()+1*mm);
  solidFrame = new G4SubtractionSolid("solidSftFrame", solidFrame, solidWindow);
  auto solidXArea = new G4Box("solidSftXArea", xarea_size.x(),
                              xarea_size.y(), xarea_size.z());
  auto solidUVArea = new G4Box("solidSftUVArea", uvarea_size.x(),
                               uvarea_size.y(), uvarea_size.z());
  auto logicFrame = new G4LogicalVolume(solidFrame, mlist.at("Al"), "logicSftFrame");
  auto logicXArea = new G4LogicalVolume(solidXArea, mlist.Scin, "logicSftXArea");
  auto logicUVArea = new G4LogicalVolume(solidUVArea, mlist.Scin, "logicSftUVArea");
  logicFrame->SetVisAttributes(G4Color::Green());
  logicXArea->SetVisAttributes(G4Color::Cyan());
  logicUVArea->SetVisAttributes(G4Color::Cyan());
  // SSD1StripLV->SetSensitiveDetector( ssdSD );
  G4ThreeVector pos = geomMan.GetGlobalPosition("SFT-Frame");
  { ///// SFT-Frame
    new G4PVPlacement(nullptr, pos, "physSftFrame", logicFrame,
                      physWorld, false, 0, m_check_overlaps);
  }
  { ///// SFT-X
    // const auto pos = geomMan.GetGlobalPosition("SFT-X");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, pos.z()+(-86.5+37+1.9)*mm), "physSftXArea",
                      logicXArea, physWorld, false, 0, m_check_overlaps);
  }
  { ///// SFT-U
    // const auto& pos = geomMan.GetGlobalPosition("SFT-U");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, pos.z()+(+86.5-53-16.1)*mm), "physSftUVArea",
                      logicUVArea, physWorld, false, 0, m_check_overlaps);
  }
  { ///// SFT-V
    // const auto& pos = geomMan.GetGlobalPosition("SFT-V");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, pos.z()+(+86.5-53+15.9)*mm), "physSftUVArea",
                      logicUVArea, physWorld, false, 1, m_check_overlaps);
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructTOF()
{
  const auto& ra2 = geomMan.GetRotAngle2("TOF") * deg;
  const auto& half_size = sizeMan.GetSize("TofSeg")*mm/2;
  const G4double pitch = geomMan.GetWirePitch("TOF")*mm;
  auto sdTOF = new TOFSD("TOF");
  G4SDManager::GetSDMpointer()->AddNewDetector(sdTOF);
  // Mother
  auto mother_solid = new G4Box("TofMotherSolid",
                                half_size.x()*NumOfSegTOF + 50.*mm,
                                half_size.y() + 50.*mm,
                                half_size.z()*2 + 50.*mm);
  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       mlist.at("Air"),
                                       "TofMotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  auto pos = geomMan.GetGlobalPosition("TOF");
  new G4PVPlacement(rot, pos, "TofMotherPHYS", mother_lv, physWorld, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // Segment
  auto segment_solid = new G4Box("TofSegmentSolid", half_size.x(),
                                 half_size.y(), half_size.z());
  auto segment_lv = new G4LogicalVolume(segment_solid,
                                        mlist.Scin,
                                        "TofSegmentLV");
  for(G4int i=0; i<NumOfSegTOF; ++i){
    segment_lv->SetVisAttributes(G4Colour::Cyan());
    segment_lv->SetSensitiveDetector(sdTOF);
    pos = G4ThreeVector((-NumOfSegTOF/2 + i)*pitch,
                        0.0,
                        2.*(- i%2 + 0.5)*half_size.z());
    new G4PVPlacement(nullptr, pos, segment_lv,
                      "TofSegmentPHYS", mother_lv, false, i, m_check_overlaps);
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructAC1()
{
  auto ac1_sd = new ACSD("AC1");
  ac1_sd->SetRefractiveIndex(1.05);
  G4SDManager::GetSDMpointer()->AddNewDetector(ac1_sd);
  const auto& ra2 = geomMan.GetRotAngle2("AC1") * deg;
  const auto& frame_size = sizeMan.GetSize("Ac1Frame") * 0.5 * mm;
  const auto& radiator_size = sizeMan.GetSize("Ac1Radiator") * 0.5 * mm;
  // Mother
  auto mother_solid = new G4Box("Ac1MotherSolid",
                                frame_size.x() + 5.*mm,
                                frame_size.y() + 5.*mm,
                                frame_size.z() + 5.*mm);
  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       mlist.at("Air"),
                                       "Ac1MotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  auto pos = geomMan.GetGlobalPosition("AC1");
  G4ThreeVector offset(0., 0., frame_size.z() - radiator_size.z());
  offset.rotateY(ra2);
  new G4PVPlacement(rot, pos + offset,
                    "Ac1MotherPHYS", mother_lv, physWorld, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // Frame
  auto frame_solid = new G4Box("Ac1FrameSolid", frame_size.x(),
                               frame_size.y(), frame_size.z());
  auto frame_lv = new G4LogicalVolume(frame_solid,
                                      mlist.at("Air"),
                                      "Ac1FrameLV");
  pos.setMag(0.);
  new G4PVPlacement(nullptr, pos, frame_lv,
                    "Ac1FramePHYS", mother_lv, false, 0, m_check_overlaps);
  // Radiator
  auto radiator_solid = new G4Box("Ac1RadiatorSolid", radiator_size.x(),
                                  radiator_size.y(), radiator_size.z());
  auto radiator_lv = new G4LogicalVolume(radiator_solid,
                                         mlist.Aerogel,
                                         "Ac1RadiatorLV");
  radiator_lv->SetSensitiveDetector(ac1_sd);
  radiator_lv->SetVisAttributes(G4Color::Magenta());
  pos.set(0., 0., -frame_size.z() + radiator_size.z());
  new G4PVPlacement(nullptr, pos, radiator_lv,
                    "Ac1RadiatorPHYS", frame_lv, false, 0);
  // Mirror
  const G4double mirror_thickness = 1.*mm/2.;
  const G4double mirror_space = 20.*mm;
  const G4ThreeVector triangle_size(1100*mm/2, frame_size.y(), 340.*mm);
  const G4double mirror_angle = std::atan2(triangle_size.z(),
                                           triangle_size.x());
  const G4ThreeVector mirror1_size((frame_size.x() - triangle_size.x())/2.,
                                   triangle_size.y(), mirror_thickness);
  const G4ThreeVector mirror2_size(std::hypot(triangle_size.x(),
                                              triangle_size.z())/2.,
                                   triangle_size.y(),
                                   mirror_thickness);
  auto mirror1_solid = new G4Box("Ac1Mirror1Solid", mirror1_size.x(),
                                 mirror1_size.y(), mirror1_size.z());
  auto mirror1_lv = new G4LogicalVolume(mirror1_solid,
                                        mlist.at("Al"),
					"Ac1Mirror1LV");
  auto mirror2_solid = new G4Box("Ac1Mirror2Solid", mirror2_size.x(),
                                 mirror2_size.y(), mirror2_size.z());
  auto mirror2_lv = new G4LogicalVolume(mirror2_solid,
                                        mlist.at("Al"),
					"Ac1Mirror2LV");
  for(G4int i=0; i<2; ++i){
    pos.set((triangle_size.x() + mirror1_size.x()) * (i*2 - 1),
            0., frame_size.z() - mirror_space);
    new G4PVPlacement(nullptr, pos, mirror1_lv,
                      "Ac1MirrorPHYS", frame_lv, false, 0);
    pos.set(triangle_size.x()/2 * (i*2 - 1),
            0., frame_size.z() - triangle_size.z()/2 - mirror_space);
    rot = new G4RotationMatrix;
    rot->rotateY(mirror_angle * (i*2 - 1));
    new G4PVPlacement(rot, pos, mirror2_lv,
                      "Ac1MirrorPHYS", frame_lv, false, 0);
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructWC()
{
  const auto& ra2 = geomMan.GetRotAngle2("WC") * deg;
  const auto& half_size_In = sizeMan.GetSize("WcSegIn") * 0.5 * mm;
  const auto& half_size_Out = sizeMan.GetSize("WcSegOut") * 0.5 * mm;
  const G4double pitch = geomMan.GetWirePitch("WC");
  auto wcSD = new WCSD("WC");
  // wcSD->SetRefractiveIndex(1.33);
  G4SDManager::GetSDMpointer()->AddNewDetector(wcSD);
  // Mother
  auto mother_solid = new G4Box("WcMotherSolid",
                                half_size_Out.x()*NumOfSegWC+10.*mm,
                                half_size_Out.y()+10.*mm,
                                half_size_Out.z()*2+10.*mm);

  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       mlist.at("Air"),
                                       "WcMotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  auto pos = geomMan.GetGlobalPosition("WC");
  new G4PVPlacement(rot, pos,
                    "WcMotherPHYS",  mother_lv, physWorld, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // Segment
  auto segment_solid = new G4Box("WcSegmentSolid", half_size_In.x(),
                                 half_size_In.y(), half_size_In.z());
  auto segment_lv = new G4LogicalVolume(segment_solid,
                                        mlist.Water,
                                        "WcSegmentLV");
  auto WCContainer     = new G4Box("WCContainer",
				   half_size_Out.x(),
				   half_size_Out.y(),
				   half_size_Out.z());
  auto WCContainer_gap = new G4Box("WCContainer_gap",
				   half_size_In.x(),
				   half_size_In.y(),
				   half_size_In.z());
  //G4RotationMatrix* rot_wccontainer_gap;
  auto rot_wccontainer_gap = new G4RotationMatrix;
  G4ThreeVector pos_wccontainer_gap(0.0, 0.0 ,0.0);
  auto solid_WCContainer
    = new G4SubtractionSolid("solid_WCContainer",
   			     WCContainer, WCContainer_gap,
   			     rot_wccontainer_gap,
			     pos_wccontainer_gap);
  auto logWCContainer = new G4LogicalVolume(solid_WCContainer,
					    mlist.Acrylic,
					    "logWCContainer");
  for(G4int i=0; i<NumOfSegWC; ++i){
    pos = G4ThreeVector((-NumOfSegWC/2 + i)*pitch,
                        0.0,
                        2.*(i%2 - 0.5)*half_size_Out.z());
    //for Vessel
    //    logWCContainer->SetVisAttributes(G4Colour::White());
    logWCContainer->SetVisAttributes(G4Colour::Cyan());
    new G4PVPlacement(nullptr, pos, logWCContainer,
                      "WcSegmentContainerPHYS", mother_lv, false, i);
    //for Water
    segment_lv->SetVisAttributes(G4Colour::Cyan());
    segment_lv->SetSensitiveDetector(wcSD);
    new G4PVPlacement(nullptr, pos, segment_lv,
                      "WcSegmentPHYS", mother_lv, false, i);

  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructVP()
{
  auto sdVP = G4SDManager::GetSDMpointer()->FindSensitiveDetector("VP", false);
  if(!sdVP){
    sdVP = new VPSD("VP");
    G4SDManager::GetSDMpointer()->AddNewDetector(sdVP);
  }
  auto solidVP = new G4Box("solidVP", 2*m/2, 2*m/2, 1*CLHEP::um/2);
  auto logicVP = new G4LogicalVolume(solidVP, mlist.at("Air"), "logicVP");
  logicVP->SetSensitiveDetector(sdVP);
  logicVP->SetVisAttributes(G4Color::Red());
  // logicVP->SetVisAttributes(G4VisAttributes::GetInvisible());
  for(G4int i=11; i<=11; ++i){
    auto pos = geomMan.GetGlobalPosition("VP"+std::to_string(i));
    auto rot = new G4RotationMatrix;
    rot->rotateY(-geomMan.GetRotAngle2("VP"+std::to_string(i))*deg);
    new G4PVPlacement(rot, pos, "physVP", logicVP, physWorld, false, i,
                      m_check_overlaps);
  }
}
