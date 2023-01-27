/* This file is the CONFIDENTIAL and PROPRIETARY information of
 * Varian Medical Systems, Inc., and is not to be distributed.
 * 
 * Copyright (c) 2017 Varian Medical Systems, Inc.
 * 
 * For information, contact Daren Sawkey  daren.sawkey@varian.com
 */

#include "TB02_TrueBeamDetectorConstruction.hh"

//#include "TB02_DetectorMessenger.hh"

//#include "TB02_PhspWorldConstruction.hh"

//#include "TB02_PhaseSpaceWriter.hh"
//#include "G4PSDirectionFlag.hh"

#include "G4PSDoseDeposit3D.hh"
#include "G4PSDoseDeposit.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4GenericPolycone.hh"
#include "G4Polycone.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Polyhedra.hh"
#include "G4ExtrudedSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4TessellatedSolid.hh"
//#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4GeometryManager.hh"
#include "G4RunManager.hh"
#include "G4TransportationManager.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4GDMLParser.hh"

#include "G4PVParameterised.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4ios.hh"
// CADMESH //
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! //
// !!! SET THIS BEFORE INCLUDING CADMESH.HH TO USE TETGEN !!! //
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! //
//#define USE_CADMESH_TETGEN
//#define CADMESH_LEXER_VERBOSE
#include "CADMesh.hh"

#include <fstream>
#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;

struct Leaf {
  G4PVPlacement *placement;
  G4int index;
  G4String bank;
} NDS120[120];

VolumeNameAndTraversalFlag TrueBeamTable[] = {
{"world", TraversedGeometry::NONE},
{"gantryMother", TraversedGeometry::NONE},
{"vacuum_LV", TraversedGeometry::NONE},
{"backscatter", TraversedGeometry::NONE},
{"BeWindow_LV", TraversedGeometry::NONE},
{"BeWindow6X_LV", TraversedGeometry::NONE},
{"VacuumChamber_LV", TraversedGeometry::NONE},
{"TargetButton_LV", TraversedGeometry::TARGET},
{"TargetNicoro_LV", TraversedGeometry::TARGET},
{"TargetBlockTop_LV", TraversedGeometry::TARGET},
{"XrayWin_LV", TraversedGeometry::NONE},
{"PrimColl_LV", TraversedGeometry::PRIMARY_COLLIMATOR},
{"carousel_bp_LV", TraversedGeometry::NONE},
{"ystageShield_LV", TraversedGeometry::YSTAGE_SHIELD},
{"FFCarousel_LV", TraversedGeometry::NONE},
{"FFNew_LV", TraversedGeometry::FLATTENING_FILTER},
{"FrictionRing_LV", TraversedGeometry::FLATTENING_FILTER},
{"openPort_LV", TraversedGeometry::FLATTENING_FILTER},
{"ICMother_LV", TraversedGeometry::IC},
{"ic_body_LV", TraversedGeometry::IC},
{"ICWindow_LV", TraversedGeometry::IC},
{"ICWinCU_LV", TraversedGeometry::IC},
{"IC_topring_LV", TraversedGeometry::IC},
{"IC_centerring_LV", TraversedGeometry::IC},
{"IC_HV_LV", TraversedGeometry::IC},
{"ICElCU_LV", TraversedGeometry::IC},
{"ICElCU_sig1_LV", TraversedGeometry::IC},
{"IC_windowring_LV", TraversedGeometry::IC},
{"ic_dose_lv", TraversedGeometry::IC},
{"coll_mother", TraversedGeometry::NONE},
{"upperJaw_LV", TraversedGeometry::UPPER_JAW},
{"lowerJaw_LV", TraversedGeometry::LOWER_JAW},
{"jaw_surround_LV", TraversedGeometry::NONE},
{"BasePlate_LogVol", TraversedGeometry::BASEPLATE},
{"BasePlate_01_LogVol", TraversedGeometry::BASEPLATE},
{"BasePlate_02_LogVol", TraversedGeometry::BASEPLATE02},
{"CornerClipper_LogVol", TraversedGeometry::BASEPLATE},
{"leaf_lv",TraversedGeometry::MLC},
{"mlc_hd120_half_target_B_leaf", TraversedGeometry::MLC},
{"mlc_hd120_half_target_A_leaf", TraversedGeometry::MLC},
{"mlc_hd120_half_iso_B_leaf", TraversedGeometry::MLC},
{"mlc_hd120_half_iso_A_leaf", TraversedGeometry::MLC},
{"mlc_hd120_quarter_target_A_leaf", TraversedGeometry::MLC},
{"mlc_hd120_quarter_target_B_leaf", TraversedGeometry::MLC},
{"mlc_hd120_quarter_iso_A_leaf", TraversedGeometry::MLC},
{"mlc_hd120_quarter_iso_B_leaf", TraversedGeometry::MLC},
{"mlc_hd120_outboard1_A_leaf", TraversedGeometry::MLC},
{"mlc_hd120_outboard60_A_leaf", TraversedGeometry::MLC},
{"mlc_hd120_outboard1_B_leaf", TraversedGeometry::MLC},
{"mlc_hd120_outboard60_B_leaf", TraversedGeometry::MLC},
{"shieldColl_LV", TraversedGeometry::SHIELD_COLLIMATOR},
{"WaterPhantom", TraversedGeometry::NONE},
{"World", TraversedGeometry::NONE},
{"MylarWindow_LV", TraversedGeometry::NONE},
};

//=======================================================================
//  TB02_TrueBeamDetectorConstruction
//=======================================================================

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TB02_TrueBeamDetectorConstruction::TB02_TrueBeamDetectorConstruction(G4String gdml_folder,G4String stl_folder, G4bool construct_phantom, G4ThreeVector phantom_position, G4ThreeVector phantom_box_size)
 //: TB02_BaseDetectorConstruction(),
 //  fLogicWorld(nullptr),
 //  fPhysiWorld(nullptr)

{
  m_gdml_folder = gdml_folder;
  m_stl_folder = stl_folder;
  fSAD = 100.*cm;

  fBeamType = "xray";

  fSimulateCollimators = true;
  fBuildBeWindow = true;
  fJawOffset = 1.5*mm;

  fFieldX1 = -20.*cm;
  fFieldX2 =  20.*cm;
  fFieldY1 = -20.*cm;
  fFieldY2 =  20.*cm;

  fJawX1 = nullptr;
  fJawX2 = nullptr;
  fJawY1 = nullptr;
  fJawY2 = nullptr;
  fUpperJawTrap_LV = nullptr;
  fLowerJawTrap_LV = nullptr;
  fBasePlate = nullptr;
  fBasePlate_LV = nullptr;
  fBasePlate1A = nullptr;
  fBasePlate1B = nullptr;
  fBasePlate2  = nullptr;
  fCornerClipper1 = nullptr;
  fCornerClipper2 = nullptr;
  fCornerClipper3 = nullptr;
  fCornerClipper4 = nullptr;
  fMLC1 = nullptr;
  fMLC2 = nullptr;
  fMLC_LV = nullptr;
  
  fMlc_type = "NDS120CAD"; // or "NDS120HD";
  fMlc_material = "W95";
  fMlc_material_NDS = "W92_5";
  
  fMlc_PV = nullptr;
  fMlc_x_size = 135.*mm;  // TODO based on HD120 drawing 
  fMlc_half_tar_A_LV = nullptr;
  fMlc_half_tar_B_LV = nullptr;
  fMlc_half_iso_A_LV = nullptr;
  fMlc_half_iso_B_LV = nullptr;
  fMlc_quarter_tar_A_LV = nullptr;
  fMlc_quarter_tar_B_LV = nullptr;
  fMlc_quarter_iso_A_LV = nullptr;
  fMlc_quarter_iso_B_LV = nullptr;
  fMlc_outboard1_A_LV   = nullptr;
  fMlc_outboard60_A_LV  = nullptr;
  fMlc_outboard1_B_LV   = nullptr;
  fMlc_outboard60_B_LV  = nullptr;


  fMylarWindow = nullptr;
  fMylarWindow_LV = nullptr;

  fSimulateVault = false;
  fVault_LV = nullptr;
  fVault = nullptr;

  fSimulateShielding = false;
  fSkullcap_pb_LV   = nullptr;
  fSkullcap_pb      = nullptr;
  fSkullcap_w_LV    = nullptr;
  fSkullcap_w       = nullptr;
  fYoke_LV          = nullptr;
  fYoke             = nullptr;
  fEars_LV          = nullptr;
  fEars             = nullptr;
  fMagnet_pole_LV   = nullptr;
  fMagnet_pole1     = nullptr;
  fMagnet_pole2     = nullptr;
  fOrbit_chamber_LV = nullptr;
  fOrbit_chamber    = nullptr;
  fGantrySteel_LV   = nullptr; // the "material" of the gantry i.e. steel
  fGantrySteel      = nullptr; // not the entire thing that rotates

  fFoil1Name      = "6E";
  fFoil1Material  = "G4_Cu";
  fFoil1Thickness = 1. *mm;
  fFoil1Radius    = 10.*mm;
  fFoil1Position  = 0. *mm;
  fFoil2Name      = "6E";
  fFoil2Material  = "G4_Cu";
  //fFoil2Thickness = 1. *mm;
  //fFoil2Radius    = 50.*mm;
  fFoil2Position  = 0. *mm;
  fApplicatorName = "None";
  fBuildCutOut    = false;
  fCutOutThickness =  13.9*mm;
  fCutOutBevelFactor = 1.0;
  fCutOutMaterial = "cerrotru";

  //target
  fTargetName     = "LowEnergy";
  //mother
  fTargetMother_LV = nullptr;
  fTargetMother    = nullptr;
  //custom
  fTargetThickness = 1. *mm;
  fTargetRadius    = 10.*mm;
  fTargetPosition  = 0. *mm;
  fTargetMaterial  = "W95";
  fTargetCustomTubs = nullptr;
  fTargetCustom_LV  = nullptr;
  fTargetCustom     = nullptr;
  fTargetCustomPC2  = nullptr;
  fTargetCustom2_LV = nullptr;
  fTargetCustom2    = nullptr;
  fTargetMaterial2  = "W95";
  //lowenergy
  fTargetLow_LV     = nullptr;
  fTargetLow        = nullptr;
  fBeWinLow_LV      = nullptr;
  fBeWinLow         = nullptr;
  fVacuumChamber_LV = nullptr;
  fVacuumChamber    = nullptr;
  fTargetNicoroLow_LV  = nullptr;
  fTargetNicoroLow     = nullptr;
  fTargetBlockTopLow_LV = nullptr; 
  fTargetBlockTopLow    = nullptr; 
  fXrayWin_LV       = nullptr;
  fXrayWin          = nullptr;
  //medium energy 
  fTargetBlockTopMed_LV = nullptr; 
  fTargetBlockTopMed    = nullptr; 
  //high energy
  fTargetHigh_LV            = nullptr;
  fTargetHigh               = nullptr;
  fTargetNicoroHigh_LV      = nullptr;
  fTargetNicoroHigh         = nullptr;
  fTargetBlockTopHigh_LV    = nullptr; 
  fTargetBlockTopHigh       = nullptr; 
  fTargetBlockBottomHigh_LV = nullptr; 
  fTargetBlockBottomHigh    = nullptr;
  fTargetWafer_LV           = nullptr;
  fTargetWafer              = nullptr;

  fShieldColl = nullptr;
  fPrimColl = nullptr;

  fGantry_LV = nullptr;  // the entire part that rotates
  fGantry    = nullptr;
  fGantryPos = 690.*mm; //fGantryPos = 865*mm;
  fGantryRot = 0.*deg;
  fGantryTheta = 0.*deg;

  fColl_LV = nullptr;
  fColl    = nullptr;
  fCollRot = 0.*deg;
  fCollPos = 430.*mm;//620*mm;


  fFlatteningFilterName = "6X";
  fFFOffset = G4ThreeVector(0.,0.,0.);

  fICMother_LV = nullptr;
  fICMother    = nullptr;
  fIC_dose_LV  = nullptr;
  fIC_dose     = nullptr;
  fUseIonChamber = true;
  fFoil1ThicknessFactor = 1.;
  fFoil2ThicknessFactor = 1.;
  fExitWindowThicknessFactor = 1.;
  //fKaptonThicknessFactor = 1.;
  fICCuHThick     = 0.0025*mm;
  fCopperThicknessFactor = 1.;

  fBuildBeWindow = true;

  fEfoil1_tubs = nullptr;
  fEfoil1_LV   = nullptr;
  fEfoil1      = nullptr;
  fEfoil2_LV   = nullptr;
  fEfoil2      = nullptr;

  fVis1 = true;
  fVis2 = true;
  fVis3 = true;
  fVis4 = true;

  fVerbosity = 2;
  fSDManager = G4SDManager::GetSDMpointer();
  fSDManager->SetVerboseLevel(fVerbosity);

  //fMessenger = new TB02_DetectorMessenger(this);

  fTargetRegion = nullptr;

  fOutputFilename = "dc_default";
  //DefineMaterials();

  fBuildWaterPhantom = construct_phantom;
  fPhantomPosition = phantom_position;
  fPhantomBoxSize = phantom_box_size;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TB02_TrueBeamDetectorConstruction::~TB02_TrueBeamDetectorConstruction()
{
  //delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* TB02_TrueBeamDetectorConstruction::Construct()
{
  //============================================================================
  //      Definitions of Solids, Logical Volumes, Physical Volumes 
  //============================================================================

  BuildWorld();
  BuildGantryMother();

  BuildBeWindow();

  // select electron or photon beam
  if (fBeamType=="electron") {
    BuildElectronGeometry();
  }
  else if (fBeamType == "xray") {
    BuildPhotonGeometry();
  }
  else {  // fBeamType neither "electron" nor "xray"
    G4ExceptionDescription ed;
    ed << "Detector Construction: Beam type " << fBeamType << "unknown!";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac100",
                FatalException, ed);
  }
  
  if (fUseIonChamber) { 
    BuildIonChamber(156.9835*mm);
  }

  BuildCollimatorMother();

  if (fSimulateCollimators) {
    BuildCollimators();
  }

  BuildShieldingCollimator();

  if (fSimulateShielding) {
    BuildShielding();
  } else {
    BuildBackscatterShield();
  }

  if (fSimulateVault) {
    BuildVault();
  }

  return fPhysiWorld;
}

void TB02_TrueBeamDetectorConstruction::GenerateTraversedTable() {
  m_id_to_traversed = generateVolumeToTraversed(TrueBeamTable, sizeof(TrueBeamTable) / sizeof(TrueBeamTable[0]), fLogicWorld);
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::BuildBeWindow() {
  //============================================================================
  // Vacuum before Be window
  //============================================================================

  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4String> physVols;
  physVols.push_back("vacuum");
  physVols.push_back("BeWindow_orbit");
  std::vector<G4String>::iterator it;
  for (it = physVols.begin(); it != physVols.end(); ++it) {
    G4VPhysicalVolume* myVol = store->GetVolume(*it, false);
    if (myVol) store->DeRegister(myVol);
  }

  G4LogicalVolumeStore* lv_store = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* lv = lv_store->GetVolume("vacuum_LV", false);
  if (lv) lv_store->DeRegister(lv);
  lv = lv_store->GetVolume("BeWindow_LV", false);
  if (lv) lv_store->DeRegister(lv);

  G4double BeThick = 0.0254*cm;
  G4double vacuumLength = 11.*mm;

  G4Tubs* vacuum = 
    new G4Tubs("vacuum", 0.*mm, 2.*cm, vacuumLength/2., 0.*deg, 360.*deg);

  fOrbitVacuum_LV =
    new G4LogicalVolume(vacuum,
                        G4Material::GetMaterial("G4_Galactic"),
                        "vacuum_LV", 0, 0, 0);

  G4VisAttributes* vacuum_VisAtt = 
    new G4VisAttributes(G4Colour(0, 0.2, 0.8, 0.2));
  vacuum_VisAtt->SetVisibility(fVis1);
  fOrbitVacuum_LV->SetVisAttributes(vacuum_VisAtt);

  fOrbitVacuum = new G4PVPlacement(0,
        G4ThreeVector(0., 0., fSAD + 2.7538*mm + BeThick/2. + vacuumLength/2. 
                      - fGantryPos),
        fOrbitVacuum_LV, "vacuum", fGantry_LV, false, 0);

  //============================================================================
  // Be Window
  //============================================================================

  if (!fBuildBeWindow) return;

  G4Tubs* BeWindow = 
    new G4Tubs("BeWindow", 0., 6.477*mm, BeThick/2., 0.*deg, 360.*deg);
  G4LogicalVolume* BeWin_LV = 
    new G4LogicalVolume(BeWindow,
                        G4Material::GetMaterial("G4_Be"), 
                        "BeWindow_LV", 0, 0, 0);

  G4VisAttributes* BeWindow_VisAtt = 
    new G4VisAttributes(G4Colour(0.5, 0.5, 1.0, 0.5));
  BeWindow_VisAtt->SetVisibility(fVis1);
  BeWindow_VisAtt->SetForceSolid(false);
  BeWin_LV->SetVisAttributes(BeWindow_VisAtt);

  //Be Window 2  (orbit chamber)
  new G4PVPlacement(0,
        G4ThreeVector(0., 0., fSAD + 2.7538*mm - fGantryPos),
        BeWin_LV, "BeWindow_orbit", fGantry_LV, false, 0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::BuildPhotonGeometry(G4bool build) {
  ////==========================================================================
  // Xray beams  (only components specific to them)
  ////==========================================================================

  BuildTarget(build);
  BuildPrimaryCollimator(build);
  BuildYStageShield(build);
  BuildFlatteningFilter(build);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::BuildElectronGeometry(G4bool build) {
  ////==========================================================================
  //// Electron beam components (only components specific to them)
  ////==========================================================================

  BuildFoil1(build);
  BuildFoil2(build);
  BuildApplicator();
  BuildCutOut();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::BuildCollimatorMother() {
  //============================================================================
  // Mother volume for shielding collimator, jaws, mlc
  //============================================================================
  G4Tubs* coll = 
    new G4Tubs("coll_mother", 0.*mm, 314.*mm, 390.*mm, 0.*deg, 360.*deg);

  fColl_LV = 
    new G4LogicalVolume(coll, G4Material::GetMaterial("G4_AIR"), "coll_mother",
        0, 0, 0);

  G4VisAttributes* VisAtt_CollM =
              new G4VisAttributes(G4Colour(0.0, 0.8, 0.8, 0.0001));
  fColl_LV->SetVisAttributes(VisAtt_CollM);
 
  G4RotationMatrix* collMotherRot = new G4RotationMatrix();
  collMotherRot->rotateZ(fCollRot);

  fColl =
    new G4PVPlacement(collMotherRot, 
      G4ThreeVector(0., 0., fCollPos - fGantryPos), fColl_LV, 
      "coll_mother", fGantry_LV, false, 0);
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::BuildShieldingCollimator() {
  //============================================================================
  // Shielding collimator
  //============================================================================
  // the steel backbone is ignored  (should be added above?, thickness of 0.597mm with open beamline width 152.4mm)
  // the coll plate is part of BuildShielding()

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  
  G4double extraThick   = 100.0*mm;

  G4double sb_temp = 225.2339*mm;
  G4double sb_offset = sb_temp+extraThick;
  
  // main part
  // the center of this is the center of the entire part
  G4Tubs* shieldcoll1 =
          new G4Tubs("shieldcoll1", 0., 141.*mm, 23.25*mm, 0.*deg, 360.*deg);

  //middle part -------------------------------

  G4Tubs*middleshieldcoll1 =
          new G4Tubs("middleshieldcoll1", 0., 203.9*mm, 6.75*mm, 0.*deg, 360.*deg);

  G4Box*middleshieldcoll2 = 
          new G4Box("middleshieldcoll2", 30.3*mm, 146.0*mm, 6.75*mm);

  G4SubtractionSolid*middleshieldcoll3 = new G4SubtractionSolid("middleshieldcoll3", 
    middleshieldcoll1, middleshieldcoll2, 0, G4ThreeVector(173.6*mm, 0., 0.));

  //-------------------------------
  
   G4UnionSolid* shieldcoll3 = 
    new G4UnionSolid("shieldcoll3", shieldcoll1, middleshieldcoll3, 0,
                     G4ThreeVector(0., 0., 30.*mm));

  // top 
  G4Tubs* shieldcoll2 =
          new G4Tubs("shieldcoll2", 0., 234.5*mm, 10.0*mm, 0.*deg, 360.*deg);
 
  G4UnionSolid* shieldcoll4 = 
    new G4UnionSolid("shieldcoll4", shieldcoll3, shieldcoll2, 0,
                     G4ThreeVector(0., 0., 40.*mm));

  
  // now try a Trd with cuts on the edges/corners of the shieldcoll

  G4Trd* scCenterTrd = new G4Trd("sc_tmp_2", 55.0*mm, 39.*mm, 55.0*mm,
                                             39.*mm, 42.0*mm);

  G4Box* scCenterBox = new G4Box("scCenter1", 80.*mm, 80.*mm, 20.*mm);

  // (39+55)/2 = 46.99  (this is half of the size of side in Trd plane at center)
  // + added 20mm (box size)  = 66.99 (now in the edge)
  // position : ctr is (39+55)/2 = 46.99
  // half width is 20, /sqrt(2)
  // 6.2 is the cut approximated from blueprints
  G4double p1 = 46.99*mm + 20./std::sqrt(2.)*mm - 6.2*mm;

  G4RotationMatrix* cR1 = new G4RotationMatrix();
  cR1->rotateZ(pi/4.);
  cR1->rotateX(pi/2. - std::atan(std::sqrt(2.)*16./80.));
 
  G4SubtractionSolid* scCenter1 = new G4SubtractionSolid("shieldcoll5", 
    scCenterTrd, scCenterBox, cR1, G4ThreeVector(p1, p1, 0.));

  G4RotationMatrix* cR2 = new G4RotationMatrix();
  cR2->rotateZ(3.*pi/4.);
  cR2->rotateX(pi/2. + std::atan(std::sqrt(2.)*16./80.));
 
  G4SubtractionSolid* scCenter2 = new G4SubtractionSolid("shieldcoll5", 
    scCenter1, scCenterBox, cR2, G4ThreeVector(-p1, p1, 0.));

  G4RotationMatrix* cR3 = new G4RotationMatrix();
  cR3->rotateZ(5.*pi/4.);
  cR3->rotateX(pi/2. - std::atan(std::sqrt(2.)*16./80.));
 
  G4SubtractionSolid* scCenter3 = new G4SubtractionSolid("shieldcoll5", 
    scCenter2, scCenterBox, cR3, G4ThreeVector(-p1, -p1, 0.));

  G4RotationMatrix* cR4 = new G4RotationMatrix();
  cR4->rotateZ(7.*pi/4.);
  cR4->rotateX(pi/2. + std::atan(std::sqrt(2.)*16./80.));
 
  G4SubtractionSolid* scCenter = new G4SubtractionSolid("shieldcoll5", 
    scCenter3, scCenterBox, cR4, G4ThreeVector(p1, -p1, 0.));


  G4SubtractionSolid* shieldcoll5 = new G4SubtractionSolid("shieldcoll5", 
              shieldcoll4,
              scCenter, 0, G4ThreeVector(0., 0., 16.75*mm));

  G4Box* shieldcoll6 = new G4Box("shieldcoll6",120.*mm, 100.*mm, 20.*mm);

  G4RotationMatrix* sc4Rot1 = new G4RotationMatrix();
  sc4Rot1->rotateX(-11.5*deg);
 
  G4SubtractionSolid* shieldcoll7 = 
    new G4SubtractionSolid("shieldcoll7", shieldcoll5, shieldcoll6, sc4Rot1, 
                     G4ThreeVector(0., 100.*mm, -40.5*mm));

  G4RotationMatrix* sc4Rot2 = new G4RotationMatrix();
  sc4Rot2->rotateX(11.5*deg);
 
  G4SubtractionSolid* shieldcoll8 = 
    new G4SubtractionSolid("shieldcoll8", shieldcoll7, shieldcoll6, sc4Rot2, 
                     G4ThreeVector(0., -100.*mm, -40.5*mm));


  G4LogicalVolume* shieldColl_LV = new G4LogicalVolume(shieldcoll8,
          G4Material::GetMaterial("Lead97Antimony"), "shieldColl_LV", 0, 0, 0);
  // subtraction solids don't always show up in vis.
  // Use this line to see solid shield coll:
  // the above statement is no longer true with newer versions of Geant4
  //G4LogicalVolume* shieldColl_LV =
  //new G4LogicalVolume(shieldcoll,Lead97Antimony,"shieldColl",0,0,0);

  G4VisAttributes* VisAtt_scSubt =
              new G4VisAttributes(G4Colour(0.0, 0.8, 0.8, 0.4));
  shieldColl_LV->SetVisAttributes(VisAtt_scSubt);

  fShieldColl = new G4PVPlacement(0,
      G4ThreeVector(0., 0., fSAD + extraThick + 0.2339*mm - sb_offset 
                    - 16.75*mm - fCollPos),
      shieldColl_LV, "shieldcoll", fColl_LV, false, 0);

  //G4Region* shieldColl_region = new G4Region("shieldcoll");
  //shieldColl_LV->SetRegion(shieldColl_region);
  //shieldColl_region->AddRootLogicalVolume(shieldColl_LV);

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::BuildCollimators(G4bool build) {
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  //TODO the MLC aren't in this yet.
  //TODO maybe should never delete these anyway?
  
  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4String> physVols {
    "UpperJaw1",            "UpperJaw2",          
    "LowerJaw1",            "LowerJaw2",
    "jaw_surround",
    "GDML_SimpleBasePlate", "GDML_BasePlate_01A", "GDML_BasePlate_01B",
    "GDML_BasePlate_02",    "GDML_MLC_corner_1",  "GDML_MLC_corner_2",
    "GDML_MLC_corner_3",    "GDML_MLC_corner_4",  
    "shieldingcorner",
    "MylarWindow", 
  };
  std::vector<G4String>::iterator it;
  for (it = physVols.begin(); it != physVols.end(); ++it) {
    G4VPhysicalVolume* myVol = store->GetVolume(*it, false);
    while (myVol) {  // need to delete all PV with same name
      store->DeRegister(myVol);
      myVol = store->GetVolume(*it, false);
    }
  }

  // This will be done with BuildMLC
  //RemoveMLC();

  if (!build) return;

  //============================================================================
  //  jaws
  //============================================================================

  //****************************
  // upper jaws
  //****************************
  G4Box* upperJawTrap = new G4Box("boxUp", 59.69*mm, 93.98*mm, 38.862*mm);
  fUpperJawTrap_LV = new G4LogicalVolume(upperJawTrap,
            G4Material::GetMaterial("W95"), "upperJaw_LV", 0, 0, 0);

  G4VisAttributes* VisAtt_upperJawTrap =
            new G4VisAttributes(G4Colour(0.2, 0.3, 1.0, 0.5));
  fUpperJawTrap_LV->SetVisAttributes(VisAtt_upperJawTrap);

  fJawY1 = new G4PVPlacement(0, G4ThreeVector(), fUpperJawTrap_LV,
    "UpperJaw1", fColl_LV, false, 0);
  
  SetJawPositionY(1, fFieldY1, true);

  fJawY2 = new G4PVPlacement(0, G4ThreeVector(), fUpperJawTrap_LV,
    "UpperJaw2", fColl_LV, false, 1);
  
  SetJawPositionY(2, fFieldY2, true);

  //****************************
  // lower jaws
  //****************************
  G4Trap* lowerJawTrap = new G4Trap("trapUp", 218.44*mm, 77.724*mm, 
                                    134.6*mm, 119.12*mm);
  fLowerJawTrap_LV = new G4LogicalVolume(lowerJawTrap,
            G4Material::GetMaterial("W95"), "lowerJaw_LV", 0, 0, 0);

  G4VisAttributes* VisAtt_lowerJawTrap = 
    new G4VisAttributes(G4Colour(0.2, 0.3, 1.0, 0.5));
  fLowerJawTrap_LV->SetVisAttributes(VisAtt_lowerJawTrap);

  fJawX1 = new G4PVPlacement(0, G4ThreeVector(), fLowerJawTrap_LV,
    "LowerJaw1", fColl_LV, false, 0);
  
  SetJawPositionX(1, fFieldX1, true);

  fJawX2 = new G4PVPlacement(0, G4ThreeVector(), fLowerJawTrap_LV,
    "LowerJaw2", fColl_LV, false, 1);
  
  SetJawPositionX(2, fFieldX2, true);

  //============================================================================
  // miscellaneous stuff around jaws
  //============================================================================
 
  G4Box* jaw_surround1 = 
    new G4Box("jaw_surround1", 215.*mm, 215.*mm, 103.*mm);

  G4Box* jaw_surround2 = 
    new G4Box("jaw_surround2", 195.*mm, 195.*mm, 104.*mm);
  
  G4SubtractionSolid* jaw_surround3 =
    new G4SubtractionSolid("jaw_surround3", jaw_surround1, jaw_surround2,
      0, G4ThreeVector());

  G4Box* jaw_surround4 = 
    new G4Box("jaw_surround4", 210.*mm, 120.*mm, 60.*mm);

  G4SubtractionSolid* jaw_surround5 =
    new G4SubtractionSolid("jaw_surround5", jaw_surround3, jaw_surround4,
      0, G4ThreeVector(0., 0., -52.*mm));

  G4LogicalVolume* jawSurround_LV = 
    new G4LogicalVolume(jaw_surround5,
                        G4Material::GetMaterial("misc"),
                        "jaw_surround_LV", 0, 0, 0);

  G4VisAttributes* VisAtt_jawsurround = 
    new G4VisAttributes(G4Colour(0.4, 0.4, 0.1, 0.1));
  jawSurround_LV->SetVisAttributes(VisAtt_jawsurround);

  new G4PVPlacement(0, G4ThreeVector(0., 0., 652.0*mm - fCollPos),
    jawSurround_LV, "jaw_surround", fColl_LV, false, 1);
  
  //============================================================================
  // Base Plate  (the baseplate itself)
  //============================================================================

  G4GDMLParser* parser = new G4GDMLParser();
  parser->Read(m_gdml_folder +  "/BasePlate.gdml", false);
  fBasePlate_LV = 
    parser->GetVolume("BasePlate_LogVol");
  delete parser;

  G4VisAttributes *VisAtt_BasePlate = 
    new G4VisAttributes(G4Colour(0.6, 0.0, 0.2, 0.4));
  fBasePlate_LV->SetVisAttributes(VisAtt_BasePlate);

  G4double BasePlate_yRef =   1.9493*mm;       // from Fastrad GDML file
  G4double BasePlate_zRef =  -4.6241*mm;       // from Fastrad GDML file
  G4double BasePlate_bottom = 533.02*mm;
  //  top:
  G4double BasePlate_z      = BasePlate_bottom + 15.24*mm;//100034445-2 + gdml
  // in z the central surface is 0.020" above center 
  // (this is not in data package!)
  G4double BasePlate_center = 0.508*mm;

  fBasePlate = new G4PVPlacement(0,
        G4ThreeVector(0.,
                      BasePlate_yRef,
                      BasePlate_z + BasePlate_zRef - fCollPos),
        fBasePlate_LV, "GDML_SimpleBasePlate", fColl_LV,
        false, 0);

  //============================================================================
  // Base Plate 01 "tungsten shielding corner"--fig12 of MCDataPackage
  //============================================================================

  G4double BasePlate1_zRef =   3.8100*mm;       // from Fastrad GDML file

  //coords of mounting holes in baseplate
  G4double bp_x = 134.62*mm; // 5.3", from MC data package
  G4double bp_y =  88.90*mm; // 3.5", from MC data package

  //coords of mounting holes in shield (long side), rel. to center
  G4double holepos_x = 22.86*mm; // 0.9" , from MC data package
  G4double holepos_y =  8.89*mm; // 0.35", from MC data package

  G4double bp_x1 = 45.72*mm;
  G4double bp_y1 = 16.51*mm;
  G4double bp_chamfer = 5.08*mm;  // the small chamfers
  G4double bp_corner = 30.861*mm; // the large diagonal cut
  G4TwoVector v11( bp_x1 - bp_corner,  bp_y1);
  G4TwoVector v21( bp_x1, +bp_y1 - bp_corner);
  G4TwoVector v22( bp_x1, -bp_y1);
  G4TwoVector v31(-bp_x1 + bp_chamfer, -bp_y1);  // these next 4 points are the
  G4TwoVector v32(-bp_x1, -bp_y1 + bp_chamfer);  // chamfered corners
  G4TwoVector v41(-bp_x1,  bp_y1 - bp_chamfer);
  G4TwoVector v42(-bp_x1 + bp_chamfer,  bp_y1);
 
  std::vector<G4TwoVector> shieldingcorner_poly 
                          {v11, v21, v22, v31, v32, v41, v42};
  G4TwoVector offset(0.0, 0.0);
  
  G4ExtrudedSolid* bp01_box = 
    new G4ExtrudedSolid("shieldingcorner", shieldingcorner_poly, 3.81*mm,
      offset, 1.0, offset, 1.0);


  G4LogicalVolume* baseplate01_LV = new G4LogicalVolume(bp01_box,
                    G4Material::GetMaterial("W95"),
                    "BasePlate_01_LogVol", 0, 0, 0);

  G4VisAttributes* VisAtt_BP1 = 
    new G4VisAttributes(G4Colour(0.0, 0.3, 0.3, 0.6));

  baseplate01_LV->SetVisAttributes(VisAtt_BP1);

  //  first one 
  G4RotationMatrix* bprot1 = new G4RotationMatrix();
  new G4PVPlacement(bprot1, 
        G4ThreeVector(-bp_x + holepos_x,
                      -bp_y + holepos_y,
                      BasePlate_z - BasePlate1_zRef + BasePlate_center 
                                  - fCollPos), 
        baseplate01_LV, "shieldingcorner", fColl_LV, false, 0);

  //  second one  
  G4RotationMatrix* bprot2 = new G4RotationMatrix();
  bprot2->rotateX(180.*deg);
  new G4PVPlacement(bprot2, 
        G4ThreeVector(-bp_x + holepos_x,
                      +bp_y - holepos_y,
                      BasePlate_z - BasePlate1_zRef + BasePlate_center 
                                  - fCollPos), 
        baseplate01_LV, "shieldingcorner", fColl_LV, false, 1);

  //============================================================================
  // Base Plate 02  "tungsten lower shielding component"--fig11 of MCdatapackage
  //============================================================================

  G4GDMLParser *ParserBasePlate_02 = new G4GDMLParser();
  ParserBasePlate_02->Read(m_gdml_folder + "/BasePlate_02.gdml", false);
  G4LogicalVolume* BasePlate2_LV = 
    ParserBasePlate_02->GetVolume("BasePlate_02_LogVol");
  delete ParserBasePlate_02;

  G4VisAttributes *VisAtt_BasePlate2 = 
    new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.6));
  BasePlate2_LV->SetVisAttributes(VisAtt_BasePlate2);

  // gdml is symmetric about 0
  G4double BasePlate2_xpos = 133.985*mm;//from MCdatapackage, and gdml symmetric

  G4double BasePlate2_zRef = 8.89*mm;  //.35"
  G4double BasePlate2_zpos = 7.62*mm;  //.3"

  G4RotationMatrix *BasePlateRot2 = new G4RotationMatrix();
  BasePlateRot2->rotateZ( 90.0*deg);
  BasePlateRot2->rotateY(180.0*deg);
  
  fBasePlate2 = new G4PVPlacement(BasePlateRot2,
      G4ThreeVector(BasePlate2_xpos,
                    0.,
                    BasePlate_z - BasePlate2_zpos + BasePlate2_zRef + 
                      BasePlate_center - fCollPos),               
      BasePlate2_LV, "GDML_BasePlate_02", fColl_LV, false, 0);

  //============================================================================
  // MLC_corner
  //============================================================================

  G4GDMLParser *ParserMLC_corner = new G4GDMLParser();
  ParserMLC_corner->Read(m_gdml_folder + "/CornerClipper.gdml", false);
  G4LogicalVolume* MLCcorner_LV = 
    ParserMLC_corner->GetVolume("CornerClipper_LogVol");
  delete ParserMLC_corner;

  G4VisAttributes *VisAtt_MLCcorner = 
    new G4VisAttributes(G4Colour(0.0, 1.0, 1.0, 0.9));
  MLCcorner_LV->SetVisAttributes(VisAtt_MLCcorner);

  //5.7" on diagonal, from MC data package
  G4double bp_hole_pos = 144.78*mm/std::sqrt(2.);
  // positions of mounting hole, relative to square corner
  G4double xpos = 22.631*mm;
  G4double ypos = 29.058*mm;// 7.112*mm;
 
  G4double thet;

  //gdml ref point is at the square corner (x,y); center for z
  G4double cornX = -39.7319*mm;       // from Fastrad GDML file
  G4double cornY = -24.3422*mm;       // from Fastrad GDML file
  G4double cornZ =   3.6830*mm;       // from Fastrad GDML file

  G4double cornZpos = BasePlate_bottom + cornZ + BasePlate_center;

  //-----------------------------------------------
  G4RotationMatrix *MLCcorner1_rot = new G4RotationMatrix();
  thet = 45.*deg;
  MLCcorner1_rot->rotateZ(thet);

  // this is in the -x-y quadrant
  fCornerClipper1 = new G4PVPlacement(MLCcorner1_rot,
      G4ThreeVector(
        // do translation in three parts
        // (1)put gdml ref point on beam axis; 
        // (2)move ref pt to the hole in baseplate
        // (3)move hole in corner
        // ref point on beam axis         // ref pt to hole   // align holes
        cornX*cos(thet)+cornY*sin(thet) - bp_hole_pos + 
          ( xpos*cos(thet)+ypos*sin(thet)),
       -cornX*sin(thet)+cornY*cos(thet) - bp_hole_pos + 
          (-xpos*sin(thet)+ypos*cos(thet)),
        cornZpos - fCollPos),
      MLCcorner_LV, "GDML_MLC_corner_1", fColl_LV, false, 0);

  //-----------------------------------------------
  G4RotationMatrix *MLCcorner2_rot = new G4RotationMatrix();
  thet = 45.*deg;
  MLCcorner2_rot->rotateX(180.*deg);
  MLCcorner2_rot->rotateZ(thet);

  fCornerClipper2 = new G4PVPlacement(MLCcorner2_rot,
      G4ThreeVector(
        cornX*cos(thet)+cornY*sin(thet) - bp_hole_pos + 
          ( xpos*cos(thet)+ypos*sin(thet)),
        // for 1st and 3rd terms, sign of y is opposite to above
        // (because of 180 deg rot about X)
       +cornX*sin(thet)-cornY*cos(thet) + bp_hole_pos - 
          (-xpos*sin(thet)+ypos*cos(thet)),
        cornZpos - fCollPos),
      MLCcorner_LV, "GDML_MLC_corner_2", fColl_LV, false, 1);

  //-----------------------------------------------
  G4RotationMatrix *MLCcorner3_rot = new G4RotationMatrix();
  thet = 225.*deg;
  MLCcorner3_rot->rotateZ(thet);

  fCornerClipper3 = new G4PVPlacement(MLCcorner3_rot,
      G4ThreeVector(
        cornX*cos(thet)+cornY*sin(thet) + bp_hole_pos + 
          ( xpos*cos(thet)+ypos*sin(thet)),
       -cornX*sin(thet)+cornY*cos(thet) + bp_hole_pos + 
          (-xpos*sin(thet)+ypos*cos(thet)),
        cornZpos - fCollPos),
      MLCcorner_LV, "GDML_MLC_corner_3", fColl_LV, false, 2);

  //-----------------------------------------------
  G4RotationMatrix *MLCcorner4_rot = new G4RotationMatrix();
  thet = 225.*deg;
  MLCcorner4_rot->rotateX(180*deg);
  MLCcorner4_rot->rotateZ(thet);

  fCornerClipper4 = new G4PVPlacement(MLCcorner4_rot,
      G4ThreeVector(
        cornX*cos(thet)+cornY*sin(thet) + bp_hole_pos + 
          ( xpos*cos(thet)+ypos*sin(thet)),
        // for 1st and 3rd terms, sign of y is opposite to above
        // (because of 180 deg rot about X)
       +cornX*sin(thet)-cornY*cos(thet) - bp_hole_pos - 
          (-xpos*sin(thet)+ypos*cos(thet)),
        cornZpos - fCollPos),
      MLCcorner_LV, "GDML_MLC_corner_4", fColl_LV, false, 3);

  //if      (fMlc_type == "NDS120CAD")   BuildCADMLC(); 
  //else if (fMlc_type == "NDS120HD") BuildHDMLC();

  //============================================================================
  //  Mylar window
  //============================================================================
  G4Box* MylarWindow =
    new G4Box("MylarWindow", 140.*mm, 140.*mm, 0.05*mm);
  fMylarWindow_LV =
    new G4LogicalVolume(MylarWindow, G4Material::GetMaterial("G4_MYLAR"),
      "MylarWindow_LV", 0, 0, 0);

  G4VisAttributes* MylarWindow_VisAtt =
    new G4VisAttributes(G4Colour(0.5, 0.5, 1.0, 0.5));
  fMylarWindow_LV->SetVisAttributes(MylarWindow_VisAtt);

  fMylarWindow = new G4PVPlacement(0, 
    G4ThreeVector(0., 0., fSAD - 55.7*cm - fCollPos),
        fMylarWindow_LV, "MylarWindow", fColl_LV, false, 0);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::BuildApplicator() { 
  //============================================================================
  //  electron applicator  (square only)
  //============================================================================
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  
  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4String> physVols;
  physVols.push_back("scraper1");
  physVols.push_back("scraper1_support");
  physVols.push_back("scraper2");
  physVols.push_back("scraper3");
  std::vector<G4String>::iterator it;
  for (it = physVols.begin(); it != physVols.end(); ++it) {
    G4VPhysicalVolume* myVol = store->GetVolume(*it, false);
    if (myVol) store->DeRegister(myVol);
  }

  G4LogicalVolumeStore* lv_store = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* lv = lv_store->GetVolume("scraper1", false);
  if (lv) lv_store->DeRegister(lv);
  lv = lv_store->GetVolume("scraper1support", false);
  if (lv) lv_store->DeRegister(lv);
  lv = lv_store->GetVolume("scraper2", false);
  if (lv) lv_store->DeRegister(lv);
  lv = lv_store->GetVolume("scraper3", false);
  if (lv) lv_store->DeRegister(lv);

  if (fApplicatorName == "None") return;

  //
  G4double delta = 1.*mm;
  // for scraper 1
  G4double pos_1 = 653.7*mm; 
  G4double thick_1 = 16.2*mm;
  // nb. values are from data package. convert to radius, then 
  // multiply by sqrt(2) because corners are at 45 deg
  G4double rtop_1 = -1.; // inside of top of scraper 
  G4double rbot_1 = -1.; // inside of bottom of scraper
  G4double r11 = -1., r12 = -1., r13 = -1., t11 = -1.;
  G4double al_width = -.1;  // width of al support for scraper 1
  // for scraper 2
  G4double pos_2 = 785.*mm;
  G4double thick_2 = 20.*mm;
  G4double sc2_width[4] = {12.*mm, 16.*mm, 32.*mm, 60.*mm};
  G4double rtop_2 = -1., rbot_2 = -1.;
  G4double r21 = -1., r22 = -1., r23 = -1., r24 = -1.; 
  G4double t21 = -1., t22 = -1., t23 = -1., t24 = -1.;
  // for scraper 3
  G4double pos_3 = 949.9*mm;
  G4double thick_3 = 13.9*mm;
  G4double rtop_3 = -1., rbot_3 = -1.;
  G4double sc3_width = 60.*mm;  // x/y dimension of scraper
  G4double r31 = -1., r32 = -1.;
 
  if (fApplicatorName == "6x6") {
    //scraper1
    rtop_1  = 70.2*mm/2.*std::sqrt(2.);
    r11     = rtop_1 +  5.*mm;
    r12     = rtop_1 + 30.*mm;
    r13     = rtop_1 + 80.*mm;
    t11     = thick_1/2. - 2.*mm;
    rbot_1  = 71.9*mm/2*std::sqrt(2.);
    al_width = 80.*mm;

    //scraper2
    rtop_2   = 65.1*mm/2.*std::sqrt(2.);
    r21      = rtop_2 + sc2_width[0];
    r22      = rtop_2 + sc2_width[1];
    r23      = rtop_2 + sc2_width[2];
    r24      = rtop_2 + sc2_width[3];
    t21      = thick_2/2. + 3.*mm;
    t22      = thick_2/2. + 1.*mm;
    t23      = thick_2/2. - 1.*mm;
    t24      = thick_2/2. - 3.*mm;
    rbot_2   = 66.7*mm/2.*std::sqrt(2.);

    //scraper3
    rtop_3   = 56.5*mm/2.*std::sqrt(2.);
    rbot_3   = 57.0*mm/2.*std::sqrt(2.);
    r31      = rtop_3 + sc3_width;
    r32      = rbot_3 + sc3_width;
  } 
  else if (fApplicatorName == "10x10") {
    //scraper1
    rtop_1  = 95.7*mm/2.*std::sqrt(2.);
    r11     = rtop_1 +  5.*mm;
    r12     = rtop_1 + 30.*mm;
    r13     = rtop_1 + 80.*mm;
    t11     = thick_1/2. - 2.*mm;
    rbot_1  = 98.1*mm/2.*std::sqrt(2.);
    al_width = 80.*mm;

    //scraper2
    rtop_2   = 95.6*mm/2.*std::sqrt(2.);
    r21      = rtop_2 + sc2_width[0];
    r22      = rtop_2 + sc2_width[1];
    r23      = rtop_2 + sc2_width[2];
    r24      = rtop_2 + sc2_width[3];
    t21      = thick_2/2. + 3.*mm;
    t22      = thick_2/2. + 1.*mm;
    t23      = thick_2/2. - 1.*mm;
    t24      = thick_2/2. - 3.*mm;
    rbot_2   = 98.1*mm/2.*std::sqrt(2.);

    //scraper3
    rtop_3   = 94.0*mm/2.*std::sqrt(2.);
    rbot_3   = 95.0*mm/2.*std::sqrt(2.);
    r31      = rtop_3 + sc3_width;
    r32      = rbot_3 + sc3_width;
  } 
  else if (fApplicatorName == "15x15") {
    //scraper1
    rtop_1  = 127.6*mm/2.*std::sqrt(2.);
    r11     = rtop_1 +  5.*mm;
    r12     = rtop_1 + 30.*mm;
    r13     = rtop_1 + 40.*mm;
    t11     = thick_1/2. - 2.*mm;
    rbot_1  = 130.8*mm/2.*std::sqrt(2.);
    al_width = 80.*mm;

    //scraper2
    rtop_2   = 133.9*mm/2.*std::sqrt(2.);
    r21      = rtop_2 + sc2_width[0];
    r22      = rtop_2 + sc2_width[1];
    r23      = rtop_2 + sc2_width[2];
    r24      = rtop_2 + sc2_width[3];
    t21      = thick_2/2. + 3.*mm;
    t22      = thick_2/2. + 1.*mm;
    t23      = thick_2/2. - 1.*mm;
    t24      = thick_2/2. - 3.*mm;
    rbot_2   = 137.4*mm/2.*std::sqrt(2.);

    //scraper3
    rtop_3   = 140.8*mm/2.*std::sqrt(2.);
    rbot_3   = 142.5*mm/2.*std::sqrt(2.);
    r31      = rtop_3 + sc3_width;
    r32      = rbot_3 + sc3_width;
  } 
  else if (fApplicatorName == "20x20") {
    //scraper1
    rtop_1  = 160.3*mm/2.*std::sqrt(2.);
    r11     = rtop_1 +  5.*mm;
    r12     = rtop_1 + 30.*mm;
    r13     = rtop_1 + 40.*mm;
    t11     = thick_1/2. - 2.*mm;
    rbot_1  = 163.5*mm/2.*std::sqrt(2.);
    al_width = 30.*mm;

    //scraper2
    rtop_2   = 172.1*mm/2.*std::sqrt(2.);
    r21      = rtop_2 + sc2_width[0];
    r22      = rtop_2 + sc2_width[0];
    r23      = rtop_2 + sc2_width[0];
    r24      = rtop_2 + sc2_width[0];
    t21      = thick_2/2. + 3.*mm;
    t22      = thick_2/2. + 1.*mm;
    t23      = thick_2/2. - 1.*mm;
    t24      = thick_2/2. - 3.*mm;
    rbot_2   = 176.6*mm/2.*std::sqrt(2.);

    //scraper3
    rtop_3   = 187.8*mm/2.*std::sqrt(2.);
    rbot_3   = 190.0*mm/2.*std::sqrt(2.);
    r31      = rtop_3 + sc3_width;
    r32      = rbot_3 + sc3_width;
  } 
  else if (fApplicatorName == "25x25") {
    //scraper1
    rtop_1  = 191.4*mm/2.*std::sqrt(2.);
    r11     = rtop_1 + 8.*mm;
    r12     = rtop_1 + 40.*mm;
    r13     = rtop_1 + 45.*mm;
    t11     = thick_1/2. - 2.*mm;
    rbot_1  = 196.2*mm/2.*std::sqrt(2.);
    al_width = 30.*mm;

    //scraper2
    rtop_2   = 210.4*mm/2.*std::sqrt(2.);
    r21      = rtop_2 + sc2_width[0];
    r22      = rtop_2 + sc2_width[1];
    r23      = rtop_2 + sc2_width[2];
    r24      = rtop_2 + sc2_width[3];
    t21      = thick_2/2. + 3.*mm;
    t22      = thick_2/2. + 1.*mm;
    t23      = thick_2/2. - 1.*mm;
    t24      = thick_2/2. - 3.*mm;
    rbot_2   = 215.9*mm/2.*std::sqrt(2.);

    //scraper3
    rtop_3   = 234.6*mm/2.*std::sqrt(2.);
    rbot_3   = 237.5*mm/2.*std::sqrt(2.);
    r31      = rtop_3 + sc3_width;
    r32      = rbot_3 + sc3_width;


  } else {
    //TODO how to do 10x6 applicator?
    //need r/z points to be function of theta
    //could use extruded polygon, new function BuildRectApplicator()
    G4ExceptionDescription ed;
    ed << "Applicator " << fApplicatorName << " not found!";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac101",
                FatalException, ed);
  }

  // scraper 1
  const G4int numRZ_1 = 8; 
  const G4double scr_r_1[numRZ_1] = 
    {rtop_1, r11, r12,   r13, r13,
     r11,    r11, rbot_1};
  const G4double scr_z_1[numRZ_1] = 
    {thick_1, thick_1, thick_1/2., thick_1/2., t11, 
     t11, 0.,  0.};

  G4Polyhedra* scraper1_poly = 
    new G4Polyhedra("scraper1", pi/4., 9.*pi/4., 4, numRZ_1, scr_r_1, scr_z_1);
    
  G4LogicalVolume* scraper1_LV =
      new G4LogicalVolume(scraper1_poly,
                          G4Material::GetMaterial("ZincZA8"),
                          "scraper1", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0., 0., fSAD - pos_1 - fCollPos),
      scraper1_LV, "scraper1", fColl_LV, false, 0);

  // Al support for scraper 1
  G4Box* scraper1_support1 = 
    new G4Box("scraper_support1", rbot_1/std::sqrt(2.) + al_width, 
              rbot_1/std::sqrt(2.) + al_width, t11/2.); 

  G4Box* scraper1_support2 = 
    new G4Box("scraper_support2", rbot_1/std::sqrt(2.) + 5.*mm, 
              rbot_1/std::sqrt(2.) + 5.*mm,
              t11/2. + delta);

  G4SubtractionSolid* scraper1_support3 =
    new G4SubtractionSolid("scraper_support3", scraper1_support1, 
      scraper1_support2, 0, G4ThreeVector());

  G4LogicalVolume* scraper1Support_LV =
    new G4LogicalVolume(scraper1_support3,
                        G4Material::GetMaterial("Aluminum6061"),
                        "scraper1_support", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0., 0., fSAD - pos_1 + t11/2. - fCollPos),
                scraper1Support_LV, "scraper1_support", fColl_LV, false, 0);

  //scraper 2
  const G4int numRZ_2 = 10;
  const G4double scr_r_2[numRZ_2] = 
    {rtop_2, r21, r22, r23,
     r24, r24, r23, r22,
     r21, rbot_2};
  const G4double scr_z_2[numRZ_2] = 
    {thick_2, thick_2, t21, t22,
     t22, t23, t23, t24,
     0., 0.};

  G4Polyhedra* scraper2poly = 
    new G4Polyhedra("scraper2", pi/4., 9.*pi/4., 4, numRZ_2, scr_r_2, scr_z_2);

  G4LogicalVolume* scraper2_LV =
      new G4LogicalVolume(scraper2poly,
                          G4Material::GetMaterial("ZincZA8"),
                          "scraper2", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0., 0., fSAD - pos_2 - fCollPos),
      scraper2_LV, "scraper2", fColl_LV, false, 0);

  //scraper 3
  const G4int numRZ_3 = 4;
  const G4double scr_r_3[numRZ_3] = {rtop_3, r31, r32, rbot_3};
  const G4double scr_z_3[numRZ_3] = {thick_3, thick_3, 0., 0.};

  G4Polyhedra* scraper3poly = 
    new G4Polyhedra("scraper3", pi/4., 9.*pi/4., 4, numRZ_3, scr_r_3, scr_z_3);
  
  G4LogicalVolume* scraper3_LV =
      new G4LogicalVolume(scraper3poly,
                          G4Material::GetMaterial("cerrotru"),
                          "scraper3", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0., 0., fSAD - pos_3 - fCollPos),
      scraper3_LV, "scraper3", fColl_LV, false, 0);

  G4VisAttributes* applicator1_VisAtt =
    new G4VisAttributes(G4Colour(0.5, 0.5, 1.0, 0.9));
  G4VisAttributes* applicator2_VisAtt =
    new G4VisAttributes(G4Colour(1.0, 0.5, 0.2, 0.9));
  scraper1_LV-> SetVisAttributes(applicator1_VisAtt);
  scraper1Support_LV->SetVisAttributes(applicator2_VisAtt);
  scraper2_LV-> SetVisAttributes(applicator1_VisAtt);
  scraper3_LV-> SetVisAttributes(applicator1_VisAtt);
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::BuildCutOut() {
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4String> physVols;
  physVols.push_back("cutout");
  std::vector<G4String>::iterator it;
  for (it = physVols.begin(); it != physVols.end(); ++it) {
    G4VPhysicalVolume* myVol = store->GetVolume(*it, false);
    if (myVol) store->DeRegister(myVol);
  }

  G4LogicalVolumeStore* lv_store = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* lv = lv_store->GetVolume("cutout", false);
  if (lv) lv_store->DeRegister(lv);

  // no cutout, if there is no applicator
  if (fApplicatorName == "None") return;
  if (fBuildCutOut    == false ) return;

  //get the applicator name, so we know how large to make the outer edge
  G4double cutout_pos   = 949.9*mm;  // bottom
  G4double delta        =   1.0*mm;  // extra space for subtraction of solids
  
  // the outside (dimension of lower surface)
  G4double x1 = 47.5*cm;
  G4double y1 = 47.5*cm;
  // ratio of upper/lower sizes of scraper 3 inside
  G4double ratio = 0.98947; 
  G4double small = 0.01*mm;  // to avoid overlaps between scraper and cutout
  //TODO 10x6
  if (fApplicatorName == "6x6")   {x1 =  28.5 *mm - small; ratio = 0.99123;}
  if (fApplicatorName == "10x10") {x1 =  47.5 *mm - small; ratio = 0.98947;}
  if (fApplicatorName == "15x15") {x1 =  71.25*mm - small; ratio = 0.98807;}
  if (fApplicatorName == "20x20") {x1 =  95.0 *mm - small; ratio = 0.98842;}
  if (fApplicatorName == "25x25") {x1 = 118.75*mm - small; ratio = 0.98779;}
  y1 = x1;
  G4TwoVector v11(x1, y1);
  G4TwoVector v21(x1, -y1);
  G4TwoVector v31(-x1, -y1);
  G4TwoVector v41(-x1, y1);
  std::vector<G4TwoVector> cutout_out_poly {v11, v21, v31, v41};
  G4TwoVector off(0.0, 0.0);
 
  G4ExtrudedSolid* cutout_out = 
    new G4ExtrudedSolid("cutout_out", cutout_out_poly, fCutOutThickness/2.,
      off, 1.0, off, 1.-((1.-ratio) * fCutOutThickness/(13.9*mm)));

  // the cutout
  // max cutout size is size of scraper
 
  G4ExtrudedSolid* cutout_in = 
    new G4ExtrudedSolid("cutout_in", fCutOutVertices, 
      fCutOutThickness/2. + delta,
      off, 1.0, off, fCutOutBevelFactor);
      // NB. because of rotations the proximal side of cutout is 2nd

  G4SubtractionSolid* cutout_subt =
    new G4SubtractionSolid("cutout_subt", cutout_out, cutout_in, 0,
      G4ThreeVector());

  fCutOut_LV = 
    new G4LogicalVolume(cutout_subt,
                        G4Material::GetMaterial(fCutOutMaterial),
                        "cutout", 0, 0, 0);

  fCutOut = new G4PVPlacement(0, 
      G4ThreeVector(0., 0., fSAD - cutout_pos + fCutOutThickness/2. - fCollPos),
      fCutOut_LV, "cutout", fColl_LV, false, 0);

  G4VisAttributes* cutout_VisAtt =
    new G4VisAttributes(G4Colour(0.0, 0.8, 0.1, 0.7));
  fCutOut_LV-> SetVisAttributes(cutout_VisAtt);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::BuildShielding(G4bool build) {

  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4String> physVols;
  physVols.push_back("skullcap_pb");
  physVols.push_back("skullcap_w");
  physVols.push_back("yoke");
  physVols.push_back("ears");
  physVols.push_back("orbit");
  physVols.push_back("gantry");
  physVols.push_back("main_coil1");
  physVols.push_back("main_coil2");
  physVols.push_back("magnet_pole1");
  physVols.push_back("magnet_pole2");
  physVols.push_back("magnet_shield1");
  physVols.push_back("magnet_shield2");
  physVols.push_back("collplate");
  physVols.push_back("collplate_out");
  physVols.push_back("cover");
  physVols.push_back("cover_in");

  std::vector<G4String>::iterator it;
  for (it = physVols.begin(); it != physVols.end(); ++it) {
    //G4cout << "Deleting volume: " << *it << G4endl;
    G4VPhysicalVolume* myVol = store->GetVolume(*it, false);
    if (myVol) store->DeRegister(myVol);
  }

  if (!fSimulateShielding || !build) {
    return;
  } else {

    //==========================================================================
    // covers
    //==========================================================================
    //TODO how to make this rotate with collimator?
    G4double cover_thick = 3.175*mm;

    G4Box* cover1 = 
      new G4Box("cover1", 320.*mm + cover_thick, 355.*mm + cover_thick, 
                280.*mm + cover_thick);

    G4Box* cover2 = 
      new G4Box("cover2", 320.*mm, 355.*mm, 280.*mm + 20.*mm);
    
    G4SubtractionSolid* cover3 =
      new G4SubtractionSolid("cover3", cover1, cover2, 0, 
        G4ThreeVector(0., 0., -20.*mm));
   
    fCover_LV = 
      new G4LogicalVolume(cover3, 
            G4Material::GetMaterial("G4_NYLON-6-6"), "cover_LV", 0, 0, 0);
  
    G4VisAttributes* cover_VisAtt = 
      new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.02));
    fCover_LV->SetVisAttributes(cover_VisAtt);

    fCover = 
      new G4PVPlacement(0, 
        G4ThreeVector(0, -40.*mm, fSAD + 30.*mm - fGantryPos), 
        fCover_LV, "cover", fGantry_LV, false, 0);

    //---------------------------------
    // part of cover inside collimator mother volume
    //---------------------------------
    G4double cover_in_h = 156.*mm;
    G4Tubs* cover_in_1 = 
      new G4Tubs("cover_in_1", 0.*mm, 310.*mm + cover_thick, 
                cover_in_h + cover_thick, 0., twopi);

    G4Tubs* cover_in_2 = 
      new G4Tubs("cover_in_2", 0.*mm, 310.*mm, cover_in_h + 10.*mm, 0., twopi);

    G4SubtractionSolid* cover_in_3 =
      new G4SubtractionSolid("cover_in_3", cover_in_1, cover_in_2, 0,
        G4ThreeVector(0., 0., 10.*mm));

    G4Box* cover_in_4 =
      new G4Box("cover_in_4", 140.*mm, 140.*mm, 10.*mm);

    G4SubtractionSolid* cover_in_5 =
      new G4SubtractionSolid("cover_in_5", cover_in_3, cover_in_4, 0, 
        G4ThreeVector(0., 0., -cover_in_h));
  
    G4LogicalVolume* cover_in_LV = 
      new G4LogicalVolume(cover_in_5,
            G4Material::GetMaterial("G4_NYLON-6-6"), "cover_in_LV", 0, 0, 0);

    cover_in_LV->SetVisAttributes(cover_VisAtt);
    
    new G4PVPlacement(0, 
        G4ThreeVector(0., 0., 390.*mm - 52.*mm - 10.*mm - cover_in_h), 
        cover_in_LV, "cover_in", fColl_LV, false, 0);

    //==========================================================================
    // skull cap
    //==========================================================================
    G4VisAttributes* skullcap1_VisAtt = 
      new G4VisAttributes(G4Colour(0.2,0.8,0.5,0.5));
    G4VisAttributes* skullcap2_VisAtt = 
      new G4VisAttributes(G4Colour(0.9,0.3,0.2,0.5));

    G4RotationMatrix* skullRot = new G4RotationMatrix();
    skullRot->rotateY(90.*deg);
   
    // cso of skull cap is cso of skullcap_pb
    G4Box* skullcap_pb = 
      new G4Box("skullcap_pb", 180.*mm, 180.*mm, 42.*mm);
    
    G4Tubs* skullcap_1 = 
      new G4Tubs("skullcap_1", 87.*mm, 160.*mm, 88.*mm, -57.*deg, 114.*deg);

    G4UnionSolid* skullcap_2 = 
      new G4UnionSolid("skullcap_2", skullcap_pb, skullcap_1, skullRot,
        G4ThreeVector(0., 0., -129.*mm));
    
    G4Box* skullcap_3 = 
      new G4Box("skullcap_3", 88.*mm, 28.*mm, 40.*mm);

   
    G4UnionSolid* skullcap_4 = 
      new G4UnionSolid("skullcap_4", skullcap_2, skullcap_3, 0,
        G4ThreeVector(0.*mm, 98.*mm, -45.*mm));

    G4UnionSolid* skullcap_5 = 
      new G4UnionSolid("skullcap_5", skullcap_4, skullcap_3, 0,
        G4ThreeVector(0.*mm, -98.*mm, -45.*mm));

    G4double orb_rad = 800.*mm;
    G4Orb* skullcap_6 =
      new G4Orb("skullcap_6", orb_rad);

    G4IntersectionSolid* skullcap_7 =
      new G4IntersectionSolid("skullcap_7", skullcap_5, skullcap_6, 0,
        G4ThreeVector(0., 0., -orb_rad + 18.*mm));

    G4Tubs* skullcap_8 = 
      new G4Tubs("skullcap_main_coil", 84.*mm, 133.*mm, 36.*mm, 0.*deg, 360.*deg);

    G4RotationMatrix* skullcap8_rot = new G4RotationMatrix();
    skullcap8_rot->rotateY(90.*deg);

    // shave some out for main coil
    G4SubtractionSolid* skullcap_9 = 
      new G4SubtractionSolid("skullcap9", skullcap_7, skullcap_8,
        skullcap8_rot,
        G4ThreeVector(125.5*mm + 1.*mm, 0.*mm, -156.*mm));

    G4SubtractionSolid* skullcap_10 = 
      new G4SubtractionSolid("skullcap10", skullcap_9, skullcap_8,
        skullcap8_rot,
        G4ThreeVector(-(125.5*mm + 1.*mm), 0.*mm, -156.*mm));

    fSkullcap_pb_LV = 
      new G4LogicalVolume(skullcap_10,
                          G4Material::GetMaterial("Lead97Antimony"),
                          "skullcap_pb_LV", 0, 0, 0);
    fSkullcap_pb_LV->SetVisAttributes(skullcap1_VisAtt);

    fSkullcap_pb = new G4PVPlacement(0,// TODO y pos is same as yoke
          //G4ThreeVector(0, -50*mm, fSAD + 244*mm + 42*mm), 
          G4ThreeVector(0., -86.*mm, fSAD + 244.*mm + 42.*mm - fGantryPos), 
          fSkullcap_pb_LV, "skullcap_pb", fGantry_LV, false, 0);

    G4Tubs* skullcap_w = 
      //new G4Tubs("skullcap_w", 92*mm, 142*mm, 82*mm, -36*deg, 72*deg);
      new G4Tubs("skullcap_w", 92.*mm, 142.*mm, 82.*mm, -36.*deg, 80.*deg);
    
    fSkullcap_w_LV =
      new G4LogicalVolume(skullcap_w,
                          G4Material::GetMaterial("W95"),
                          "skullcap_w_LV", 0, 0, 0);
    fSkullcap_w_LV->SetVisAttributes(skullcap2_VisAtt);

    fSkullcap_w = new G4PVPlacement(skullRot,
          G4ThreeVector(0., 0., -129.*mm),
          fSkullcap_w_LV, "skullcap_w", fSkullcap_pb_LV, false, 0);

    //==========================================================================
    // bend magnet yoke
    //==========================================================================

    // cso for entire part is cso of back part
    G4double yoke_delta = 1.*mm;
    // the back part
    G4Box* yoke1 =  
      new G4Box("yoke1", 295.*mm, 79.*mm, 122.*mm);
    // the sides
    G4Box* yoke2 = 
      new G4Box("yoke2", 66.5*mm, 214.*mm - yoke_delta, 122.*mm);
    G4UnionSolid* yoke3 = 
      new G4UnionSolid("yoke3", yoke1, yoke2, 0, 
        G4ThreeVector(295.*mm - 133.*mm/2., (428.*mm - 158.*mm)/2., 0.));
    G4UnionSolid* yoke4 = 
      new G4UnionSolid("yoke4", yoke3, yoke2, 0, 
        G4ThreeVector(-(295.*mm - 133.*mm/2.), (428.*mm - 158.*mm)/2., 0.));

    // shave the corners
    G4Orb* yoke_orb = 
      new G4Orb("yoke_orb", 380.*mm); 
    G4IntersectionSolid* yoke5 =
      new G4IntersectionSolid("yoke5", yoke4, yoke_orb, 0, 
        G4ThreeVector(0., (270.*mm - 79.*mm - 10.*mm), -198.*mm));

    // the round parts
    G4Tubs* yoke_tub =
      new G4Tubs("yoke_tube", 0., 83.5*mm, 36.5*mm + yoke_delta, 0.*deg, 360.*deg);
    G4RotationMatrix* yoke_tub_rot = new G4RotationMatrix();
    yoke_tub_rot->rotateY(90.*deg);

    G4UnionSolid* yoke6 =
      new G4UnionSolid("yoke6", yoke5, yoke_tub, yoke_tub_rot, 
        G4ThreeVector(125.5*mm + yoke_delta, 213.*mm, 8.*mm));

    G4UnionSolid* yoke7 =
      new G4UnionSolid("yoke7", yoke6, yoke_tub, yoke_tub_rot, 
        G4ThreeVector(-(125.5*mm + yoke_delta), 213.*mm, 8.*mm));

    // the lower part, going into gantry weldment
    G4Box* yoke8 = 
      new G4Box("yoke8", 230.*mm, 181.5*mm, 24.*mm);

    G4Box* yoke9 = 
      new G4Box("yoke8", 180.*mm, 181.5*mm, 25.*mm);
   
    G4SubtractionSolid* yoke10 =
      new G4SubtractionSolid("yoke10", yoke8, yoke9, 0, 
              G4ThreeVector(0., 75.*mm, 0.));

    G4UnionSolid* yoke11 = 
      new G4UnionSolid("yoke11", yoke7, yoke10, 0, 
        G4ThreeVector(0., 165.*mm + 26.5*mm - 25.*mm, -146.*mm));

    fYoke_LV = 
      new G4LogicalVolume(yoke11,
                          G4Material::GetMaterial("SSlowcarbon"),
                          "yoke_LV", 0, 0, 0);
    G4VisAttributes* yoke_VisAtt = 
      new G4VisAttributes(G4Colour(0.1, 0.4, 0.8, 0.5));
    fYoke_LV->SetVisAttributes(yoke_VisAtt);

    fYoke = new G4PVPlacement(0, // TODO get position in Y
        G4ThreeVector(0., -213.*mm - 86.*mm, fSAD + 122.*mm - fGantryPos),  
        fYoke_LV, "yoke", fGantry_LV, false, 0);

    //==========================================================================
    // main magnet coil (approx. as solid Cu)
    //==========================================================================

    G4Tubs* coil = 
      new G4Tubs("main_coil", 85.5*mm, 131.5*mm, 35.*mm, 0.*deg, 360.*deg);

    fMainCoil_LV = 
      new G4LogicalVolume(coil,
                          G4Material::GetMaterial("G4_Cu"),
                          "coil_LV", 0, 0, 0);

    G4VisAttributes* coil_VisAtt = 
      new G4VisAttributes(G4Colour(0.1, 0.9, 0.3, 0.3));
    fMainCoil_LV->SetVisAttributes(coil_VisAtt);

    G4RotationMatrix* coil_rot = new G4RotationMatrix();
    coil_rot->rotateY(90.*deg);

    fMainCoil1 = new G4PVPlacement(coil_rot,
        G4ThreeVector(125.5*mm + yoke_delta, -86.*mm, 
                      fSAD + 130.*mm - fGantryPos),
        fMainCoil_LV, "main_coil1", fGantry_LV, false, 0);

    fMainCoil2 = new G4PVPlacement(coil_rot,
        G4ThreeVector(-(125.5*mm + yoke_delta), -86.*mm, 
          fSAD + 130.*mm - fGantryPos),
        fMainCoil_LV, "main_coil2", fGantry_LV, false, 1);

    //==========================================================================
    // bend magnet pole
    //==========================================================================

    G4Tubs* pole = 
      new G4Tubs("bend_magnet_pole", 0.*mm, 85.*mm, 41.75*mm, 0.*deg, 360.*deg);

    fMagnet_pole_LV = 
      new G4LogicalVolume(pole,
                          G4Material::GetMaterial("SSlowcarbon"),
                          "magnetpole_LV", 0, 0, 0);

    G4VisAttributes* pole_VisAtt = 
      new G4VisAttributes(G4Colour(0.9, 0.1, 0.9, 0.6));
    fMagnet_pole_LV->SetVisAttributes(pole_VisAtt);

    G4RotationMatrix* pole_rot = new G4RotationMatrix();
    pole_rot->rotateY(90.*deg);

    fMagnet_pole1 = new G4PVPlacement(pole_rot,
        G4ThreeVector(45.5*mm + yoke_delta, -86.*mm, 
          fSAD + 130.*mm - fGantryPos),
        fMagnet_pole_LV, "magnet_pole1", fGantry_LV, false, 0);

    fMagnet_pole2 = new G4PVPlacement(pole_rot,
        G4ThreeVector(-(45.5*mm + yoke_delta), -86.*mm, 
          fSAD + 130.*mm - fGantryPos),
        fMagnet_pole_LV, "magnet_pole2", fGantry_LV, false, 1);

    //==========================================================================
    // bend magnet shield
    //==========================================================================
    
    G4Box* magnet_shield =
      new G4Box("magnet_shield", 40.*mm, 85.*mm, 9.*mm);

    G4LogicalVolume* magnetShield_LV =
      new G4LogicalVolume(magnet_shield,
                          G4Material::GetMaterial("Lead97Antimony"),
                          "magnet_shield_LV", 0, 0, 0);

    G4VisAttributes* ms_VisAtt =
      new G4VisAttributes(G4Colour(0.0, 0.4, 0.6, 0.5));
    magnetShield_LV->SetVisAttributes(ms_VisAtt);

    //G4PVPlacement* magnetShield1 = 
      new G4PVPlacement(0, 
        G4ThreeVector(45.5*mm + yoke_delta, -86.*mm, 
          fSAD + 34.*mm - fGantryPos),
        magnetShield_LV, "magnet_shield1", fGantry_LV, false, 0);
    
    //G4PVPlacement* magnetShield2 = 
      new G4PVPlacement(0, 
        G4ThreeVector(-45.5*mm - yoke_delta, -86.*mm, 
          fSAD + 34.*mm - fGantryPos),
        magnetShield_LV, "magnet_shield2", fGantry_LV, false, 1);
    //==========================================================================
    // shielding ears
    //==========================================================================


    G4Box* ears1 =
      new G4Box("shielding_ears1", 300.*mm, 50.*mm, 122.*mm);

    // add some material to seal up around the gantry
    G4Box* ears2 =
      new G4Box("ears2", 300.*mm, 70.5*mm, 50.*mm);

    G4UnionSolid* ears3 =
      new G4UnionSolid("ears3", ears1, ears2, 0, 
        G4ThreeVector(0, 120.*mm, -72.*mm));
   
    // a part extending below target
    // extends to z=973 mm
    G4Box* ears4 =
      new G4Box("ears4", 230.*mm, 50.*mm, 15.*mm);

    G4UnionSolid* ears5 =
      new G4UnionSolid("ears5", ears3, ears4, 0, 
          G4ThreeVector(0., 0., -134.*mm));

    fEars_LV =
      new G4LogicalVolume(ears5,
                          G4Material::GetMaterial("Lead97Antimony"),
                          "ears_LV", 0, 0, 0);
    G4VisAttributes* ears_VisAtt = 
      new G4VisAttributes(G4Colour(0.4, 0.4, 0.1, 0.5));
    fEars_LV->SetVisAttributes(ears_VisAtt);

    fEars = new G4PVPlacement(0,
        //G4ThreeVector(0, 136*mm, fSAD + 122*mm - 2.5*mm), 
        //G4ThreeVector(0, 99*mm, fSAD + 122*mm - 2.5*mm), 
        G4ThreeVector(0., 99.*mm, fSAD + 122.*mm - fGantryPos), 
            // the -5 mm compensates for adjustment of gantry position
        fEars_LV, "ears", fGantry_LV, false, 0);

   
    //==========================================================================
    // orbit chamber
    //==========================================================================

    G4Box* orbit_chamber1 =
      new G4Box("orbit_chamber1", 4.*mm, 60.*mm, 73.*mm);
  
    G4Box* orbit_chamber2 =
      new G4Box("orbit_chamber", 16.*mm, 60.*mm, 5.*mm);

    G4UnionSolid* orbit_chamber3 = 
      new G4UnionSolid("orbit_chamber3", orbit_chamber1, orbit_chamber2,
                       0, G4ThreeVector(0., 0., -69.*mm));
                       
    fOrbit_chamber_LV = 
      new G4LogicalVolume(orbit_chamber3,
                          G4Material::GetMaterial("G4_Cu"),
                          "orbit_chamber_LV", 0, 0, 0);
  
    G4VisAttributes* orbit_VisAtt = 
      new G4VisAttributes(G4Colour(0.5, 0.4, 0.8, 0.4));
    fOrbit_chamber_LV->SetVisAttributes(orbit_VisAtt);
    
    fOrbit_chamber = new G4PVPlacement(0, 
        G4ThreeVector(0., -40.*mm, fSAD + 88.*mm - fGantryPos),
        fOrbit_chamber_LV, "orbit", fGantry_LV, false, 0);
   
    // ***********************************************************************
    // collimator mounting plate
    G4Box* collplate1 = 
      new G4Box("collplate1", 313.99*mm, 348.99*mm, 25.4*mm);

    G4Tubs* collplate2 =
      new G4Tubs("collplate2", 0., 235.*mm, 26.*mm, 0.*deg, 360.*deg);

    G4SubtractionSolid* collplate3 =
      new G4SubtractionSolid("collplate3", collplate1, collplate2, 0, 
        G4ThreeVector(0., 39.*mm, 0.));

    G4Tubs* collplate4 =
      new G4Tubs("collplate4", 314.*mm, 600.*mm, 26.1*mm, 0.*deg, 360.*deg);

    //this subtracted part should not be part of fColl_LV
    //this is an approximation; this entire LV does not rotate but this way the
    //coll mother can be a G4Tubs
    //This piece added below
    G4SubtractionSolid* collplate5 =
      new G4SubtractionSolid("collplate5", collplate3, collplate4, 0, 
        G4ThreeVector(0., 39.*mm, 0.));

    G4LogicalVolume* collplate_LV = 
      new G4LogicalVolume(collplate5,
                          G4Material::GetMaterial("SS_A36"),
                          "collplate", 0, 0, 0);

    G4VisAttributes* VisAtt_collplate =
                new G4VisAttributes(G4Colour(1.0, 0.8, 0.8, 0.3));
    collplate_LV->SetVisAttributes(VisAtt_collplate);

    //TODO:  here, 100 and 325.2339 are local variables for the shielding 
    //  collimator
    new G4PVPlacement(0, 
      G4ThreeVector(0, -39*mm, 
        fSAD + 100.*mm + 0.2339*mm - 325.2339*mm - 12.5*mm + 25.4*mm - fCollPos),
      collplate_LV, "collplate", fColl_LV, false, 0);

    // ************************
    // part outside of fColl_LV
    // ************************
    G4Tubs* collplate_out_2 =
      new G4Tubs("collplate_out_2", 0.*mm, 314.01*mm, 26.*mm, 0.*deg, 360.*deg);

    G4SubtractionSolid* collplate_out_3 =
      new G4SubtractionSolid("collplate_out_3", collplate1, collplate_out_2,
        0, G4ThreeVector(0., 39.*mm, 0.));

    //TODO this part is inside the collimator mother, but doesn't rotate with
    //collimator!
    // approximate it as on the outside of coll mother
    G4Tubs* collplate_out_4 = 
      // reasonably correct geometry:
      //new G4Tubs("collplate_out_4", 250.*mm, 305.*mm, 25.4*mm, 220.*deg, 
      // approximated:
      new G4Tubs("collplate_out_4", 315.*mm, 370.*mm, 25.4*mm, 220.*deg, 
        100.*deg);

    G4UnionSolid* collplate_out_5 =
      new G4UnionSolid("collplate_out_5", collplate_out_3, collplate_out_4, 0,
        G4ThreeVector(0, 30.*mm, -50.8*mm));

    G4LogicalVolume* collplate_out_LV = 
      //new G4LogicalVolume(collplate_out_3,
      new G4LogicalVolume(collplate_out_5,
                          G4Material::GetMaterial("SS_A36"),
                          "collplate_LV", 0, 0, 0);

    G4VisAttributes* VisAtt_collplate_out =
                new G4VisAttributes(G4Colour(0.2, 0.8, 0.8, 0.7));
    collplate_out_LV->SetVisAttributes(VisAtt_collplate_out);

    new G4PVPlacement(0, 
      G4ThreeVector(0, -39.*mm, 
        fSAD + 100.*mm + 0.2339*mm - 325.2339*mm- 12.5*mm + 25.4*mm 
          - fGantryPos),
      collplate_out_LV, "collplate_out", fGantry_LV, false, 0);

    //==========================================================================
    // gantry
    //==========================================================================
    G4double gantry_thick  = 60.*mm;
    G4double gantry_height = 89.*mm;

    G4Box* gantry1 = 
      new G4Box("gantry1", 300.*mm, 336.*mm, gantry_height);
    G4Box* gantry2 = 
      new G4Box("gantry2", 300*mm - gantry_thick, 
                336.*mm - gantry_thick, 
                gantry_height + 5.*mm);

    G4SubtractionSolid* gantry = 
      new G4SubtractionSolid("gantry", gantry1, gantry2, 0, 
                             G4ThreeVector(0., 0., 0.));

    fGantrySteel_LV = 
      new G4LogicalVolume(gantry,
                          G4Material::GetMaterial("SS_A36"),
                          "gantrysteel_LV", 0, 0, 0);
    
    G4VisAttributes* gantry_VisAtt = 
      new G4VisAttributes(G4Colour(0.2, 0.4, 0.0, 0.05));
    fGantrySteel_LV->SetVisAttributes(gantry_VisAtt);

    fGantrySteel = new G4PVPlacement(0,
        // set the y position such that inner edge is at -315mm
        // to line up with yoke
        G4ThreeVector(0., -39.*mm, fSAD - gantry_height - fGantryPos),
        fGantrySteel_LV, "gantry", fGantry_LV, false, 0);

  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::BuildVault(G4bool build) {
  //============================================================================
  // Concrete vault
  //============================================================================

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  
  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4String> physVols;
  physVols.push_back("vault_walls");

  std::vector<G4String>::iterator it;
  for (it = physVols.begin(); it != physVols.end(); ++it) {
    //G4cout << "Deleting volume: " << *it << G4endl;
    G4VPhysicalVolume* myVol = store->GetVolume(*it, false);
    if (myVol) store->DeRegister(myVol);
  }

  if (!fSimulateVault || !build) {
    return;
  }
  G4double vault_size = 4.*m;
  G4double wall_thick = 1.*m;
  G4double vault_outer_size = vault_size + wall_thick;
  G4Box* vault_outer = 
    new G4Box("vault_outer", vault_outer_size/2., vault_outer_size/2., 
              vault_outer_size/2.);
  G4Box* vault_inner = 
    new G4Box("vault_inner", vault_size/2., vault_size/2., vault_size/2.);

  G4SubtractionSolid* walls = new G4SubtractionSolid("vault_walls", 
                      vault_outer, vault_inner, 0, G4ThreeVector());

  fVault_LV =
    new G4LogicalVolume(walls,
                        G4Material::GetMaterial("Concrete"),
                        "vault_LV", 0, 0, 0);

  G4VisAttributes* vault_VisAtt = 
    new G4VisAttributes(G4Colour(0.1, 0.5, 0.2, 0.2));
  vault_VisAtt->SetVisibility(fVis1);
  fVault_LV->SetVisAttributes(vault_VisAtt);

  fVault = new G4PVPlacement(0,
        G4ThreeVector(),
        fVault_LV, "vault_walls", fLogicWorld, false, 0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::BuildBackscatterShield(G4bool build) {
  //============================================================================
  // backscatter killer
  //============================================================================

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  
  if (!build) { 
    G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
    G4VPhysicalVolume* myVol = store->GetVolume("backscatter", false);
    if (myVol) store->DeRegister(myVol);
    return;
  }
  
  G4VisAttributes* backscatt_VisAtt = 
    new G4VisAttributes(G4Colour(0.2, 0.0, 0.8, 0.9));
  backscatt_VisAtt->SetVisibility(fVis2);

  G4Sphere* backScatter = 
    new G4Sphere("backscatter", 8.*mm, 9.*mm, 0.*deg, 180.*deg, 
                                0.*deg, 360.*deg);
  fBackScatter_LV = new G4LogicalVolume(backScatter,
                  G4Material::GetMaterial("G4_Galactic"),
                  "backscatter", 0, 0, 0);

  G4RotationMatrix* backscatterRot = new G4RotationMatrix();
  backscatterRot->rotateX(-90.0*deg);
  
  fBackScatter = new G4PVPlacement(backscatterRot,
        G4ThreeVector(0., 0., -5.*mm),
        fBackScatter_LV, "backscatter", fOrbitVacuum_LV, false, 0);

  fBackScatter_LV->SetVisAttributes(backscatt_VisAtt);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::BuildTarget(G4bool build) {
  //============================================================================
  // Target
  //============================================================================

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  
  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4String> physVols;
  physVols.push_back("customTarget");
  physVols.push_back("customTarget2");
  physVols.push_back("BeWindowLow");
  physVols.push_back("VacuumChamber");
  physVols.push_back("TargetButtonLow");
  physVols.push_back("NicoroLow");
  physVols.push_back("TargetBlockTopLow");
  physVols.push_back("XrayWindow");
  physVols.push_back("TargetBlockTopMed");
  physVols.push_back("TargetButtonHigh");
  physVols.push_back("TargetBlockTopHigh");
  physVols.push_back("NicoroHigh");
  physVols.push_back("brazeWafer");
  physVols.push_back("TargetBlockBottomHigh");
  physVols.push_back("TargetBlockTopImaging");
  //physVols.push_back("targetMotherPV");

  std::vector<G4String>::iterator it;
  for (it = physVols.begin(); it != physVols.end(); ++it) {
    //G4cout << "Deleting volume: " << *it << G4endl;
    G4VPhysicalVolume* myVol = store->GetVolume(*it, false);
    if (myVol) store->DeRegister(myVol);
  }

  if (fTargetRegion) {
    G4RegionStore::GetInstance()->DeRegister(fTargetRegion);
  }
  //G4Region* reg = G4RegionStore::GetInstance()->GetRegion("target");
  //G4cout << "pointer to target region: " << reg << G4endl;
  //if (reg) {
  //  G4cout << "name: " << reg->GetName() << G4endl;
  //}
  fTargetRegion = new G4Region("target");

  if (!build) return;
  if (fBeamType != "xray") return;

  G4VisAttributes* target_VisAtt = 
    new G4VisAttributes(G4Colour(0.2, 0.2, 0.8, 0.9));
  target_VisAtt->SetVisibility(fVis2);

  G4VisAttributes* target2_VisAtt = 
    new G4VisAttributes(G4Colour(0.5, 0.8, 0.2, 0.9));
  target2_VisAtt->SetVisibility(fVis2);
  
  G4VisAttributes* target_block_VisAtt = 
    new G4VisAttributes(G4Colour(0.5, 0.0, 0.5, 0.5));
  target_block_VisAtt->SetVisibility(fVis2);


  //==========================================================================
  if (fTargetName == "Custom") {
  //==========================================================================
    fTargetCustomTubs = 
      new G4Tubs("target", 0.*mm, fTargetRadius, fTargetThickness/2.,
                 0.*deg, 360.*deg);
    fTargetCustom_LV = new G4LogicalVolume(fTargetCustomTubs,
                    G4Material::GetMaterial(fTargetMaterial),
                    "customTarget_LV", 0, 0, 0);

    fTargetCustom = new G4PVPlacement(0,
          G4ThreeVector(0., 0., 
                        fSAD + fTargetPosition - fTargetThickness/2. 
                        - fGantryPos),
          fTargetCustom_LV, "customTarget", fGantry_LV, false, 0);

    fTargetCustom_LV->SetVisAttributes(target_VisAtt);

    // rotate so coordinates are "logical"
    G4RotationMatrix* targetRot2 = new G4RotationMatrix();
    targetRot2->rotateX(pi);
    if (fTargetCustom2_z.size() > 2) {
      fTargetCustomPC2 =
        new G4GenericPolycone("target2", 0., twopi, static_cast<int>(fTargetCustom2_z.size()),
                              &fTargetCustom2_r[0], &fTargetCustom2_z[0]);

      fTargetCustom2_LV = new G4LogicalVolume(fTargetCustomPC2,
                      G4Material::GetMaterial(fTargetMaterial2),
                      "customTarget2_LV", 0, 0, 0);

      fTargetCustom2 = new G4PVPlacement(targetRot2,
            G4ThreeVector(0., 0.,
                          fSAD + fTargetPosition
                          - fTargetThickness - fGantryPos),
            fTargetCustom2_LV, "customTarget2", fGantry_LV, false, 0);

      fTargetCustom2_LV->SetVisAttributes(target2_VisAtt);
    }

    fTargetCustom_LV->SetRegion(fTargetRegion);
    fTargetRegion->AddRootLogicalVolume(fTargetCustom_LV);

    if (fTargetCustom2_LV) {
      fTargetCustom2_LV->SetRegion(fTargetRegion);
      fTargetRegion->AddRootLogicalVolume(fTargetCustom2_LV);
    }
  }

  //==========================================================================
  else if (fTargetName=="LowEnergy") {
  //==========================================================================
    //------------------------------------------------------------------------
    // Be window
    //------------------------------------------------------------------------
    G4double BeThick6X = 0.0254*cm;  // same as BeThick for orbit chamber
    G4Tubs* BeWindow6X = 
      new G4Tubs("BeWindow6X", 0.*mm, 6.477*mm, BeThick6X/2., 0.*deg, 360.*deg);
    fBeWinLow_LV = 
      new G4LogicalVolume(BeWindow6X,
                          G4Material::GetMaterial("G4_Be"), 
                          "BeWindow6X_LV", 0, 0, 0);

    fBeWinLow = new G4PVPlacement(0,
          G4ThreeVector(0., 0., fSAD + 0.751*mm - fGantryPos),
          fBeWinLow_LV, "BeWindowLow", fGantry_LV, false, 0);

    G4VisAttributes* VisAtt_Be6X = 
      new G4VisAttributes(G4Colour(0.5, 0.5, 1.0, 0.5));
    VisAtt_Be6X->SetVisibility(fVis2);
    fBeWinLow_LV->SetVisAttributes(VisAtt_Be6X);

    //========================================================================
    // Vacuum Chamber between Be window and target button
    //========================================================================

    const G4int vc_points = 6;
    G4double r_vc[vc_points] = 
      { 0.   *mm,  3.785*mm,  3.785*mm,  4.9  *mm,  4.9  *mm,  0.   *mm};
    G4double z_vc[vc_points] = 
      { 0.261*mm,  0.261*mm,  0.324*mm,  0.324*mm,  0.624*mm,  0.624*mm};

    G4Polycone* vacuum_chamber = new G4Polycone("VacuumChamber_PhysVol",
              0.*deg, 360.*deg, vc_points, r_vc, z_vc);

    fVacuumChamber_LV =
      new G4LogicalVolume(vacuum_chamber, 
                          G4Material::GetMaterial("G4_Galactic"),
                          "VacuumChamber_LV",0,0,0);

    fVacuumChamber = new G4PVPlacement(0,
          G4ThreeVector(0., 0., fSAD - fGantryPos), fVacuumChamber_LV,
          "VacuumChamber", fGantry_LV, false, 0);

    G4VisAttributes* VisAtt_VacuumChamber = 
      new G4VisAttributes(G4Colour(1., 0., 0., 0.6));
    VisAtt_VacuumChamber->SetVisibility(fVis2);
    fVacuumChamber_LV->SetVisAttributes(VisAtt_VacuumChamber);

    //------------------------------------------------------------------------
    // 6X Target Button
    //------------------------------------------------------------------------

    G4double nicoroHThick = 0.0508*mm; // half the nicoro thickness

    const G4int tb6pts = 5;
    G4double r_tb[tb6pts] = 
      {0.0*mm, 3.7592*mm, 3.7592*mm,  3.5052*mm,  0.0  *mm};
    G4double z_tb[tb6pts] = 
      {0.0*mm,    0.0*mm, -0.381*mm, -0.635 *mm, -0.635*mm};

    G4double Target_z = -0.159*mm; 

    G4Polycone* target_button = new G4Polycone("TargetButton",0.*deg, 360.*deg,
                                        tb6pts, r_tb, z_tb);
    fTargetLow_LV = 
      new G4LogicalVolume(target_button,
                          G4Material::GetMaterial("G4_W"),
                          "TargetButton_LV", 0, 0, 0);

    fTargetLow = new G4PVPlacement(0,
        G4ThreeVector(0., 0., fSAD - Target_z + 2.*nicoroHThick - fGantryPos),
          fTargetLow_LV, "TargetButtonLow", fGantry_LV, 
          false, 0);

    fTargetLow_LV->SetVisAttributes(target_VisAtt);

    //------------------------------------------------------------------------
    // nicoro brazing sheet
    //------------------------------------------------------------------------
    // not sure what the material should be: take it to be nicoro (BAu-3) 
    // (other choice is nicoro 80)

    G4Tubs* TargetNicoro = 
      new G4Tubs("nicoro", 0.0*mm, 3.0734*mm, nicoroHThick, 0.*deg, 360.*deg);

    fTargetNicoroLow_LV = 
      new G4LogicalVolume(TargetNicoro,
                          G4Material::GetMaterial("Nicoro"), 
                          "TargetNicoro_LV",0,0,0);

    fTargetNicoroLow = new G4PVPlacement(0, 
      G4ThreeVector(0., 0., fSAD - Target_z - 0.635*mm + nicoroHThick 
                    - fGantryPos),
      fTargetNicoroLow_LV, "NicoroLow", fGantry_LV, false, 0);

    G4VisAttributes* VisAtt_Nicoro = 
      new G4VisAttributes(G4Colour(0.0, 0.5, 0.5));
    VisAtt_Nicoro->SetVisibility(fVis2);
    fTargetNicoroLow_LV->SetVisAttributes(VisAtt_Nicoro);

    //------------------------------------------------------------------------
    //  Target block top
    //------------------------------------------------------------------------

    const G4int tbt_pts = 31;
    G4double r_pc[tbt_pts] = {
      0.0   *mm, 3.785 *mm, 3.785 *mm, 4.9   *mm,
      4.9   *mm, 6.6   *mm, 6.6   *mm, 8.2   *mm,
      8.2   *mm, 7.1   *mm, 7.1   *mm, 6.6   *mm,
      6.6   *mm, 4.502 *mm,
      3.2137*mm, 3.1928*mm, 3.1660*mm, 3.1338*mm,
      3.0963*mm, 3.0538*mm, 3.0067*mm, 2.9552*mm,
      2.8999*mm, 2.8410*mm, 2.7791*mm, 2.7146*mm,
      2.6481*mm, 2.5799*mm, 2.5107*mm, 2.4410*mm, 0.0*mm};

    G4double z_pc[tbt_pts] = {
      -0.476 *mm, -0.476 *mm,  0.324 *mm,  0.324*mm,
       0.624 *mm,  0.624 *mm,  1.124 *mm,  1.124*mm,
      -5.876 *mm, -5.876 *mm, -8.876 *mm, -8.876*mm,
      -8.376 *mm, -8.376 *mm,
      -3.5689*mm, -3.5024*mm, -3.4379*mm, -3.3760*mm,
      -3.3171*mm, -3.2618*mm, -3.2103*mm, -3.1632*mm,
      -3.1207*mm, -3.0832*mm, -3.0510*mm, -3.0242*mm,
      -3.0033*mm, -2.9882*mm, -2.9790*mm, -2.9760*mm, -2.9760*mm};

    G4GenericPolycone* target_block_top = 
      new G4GenericPolycone("TargetBlockTop_PhysVol",
                            0.*deg, 360.*deg, tbt_pts, r_pc, z_pc);

    fTargetBlockTopLow_LV = 
      new G4LogicalVolume(target_block_top,
                          G4Material::GetMaterial("COPPER_GLIDCOP"),
                          "TargetBlockTop_LV", 0, 0, 0);

    fTargetBlockTopLow_LV->SetVisAttributes(target_block_VisAtt);

    fTargetBlockTopLow = 
      new G4PVPlacement(0, G4ThreeVector(0., 0., fSAD - fGantryPos),
                fTargetBlockTopLow_LV, "TargetBlockTopLow", 
                fGantry_LV, false, 0);

    //------------------------------------------------------------------------
    // X-ray window
    //------------------------------------------------------------------------

    G4Tubs* XrayWin = 
      new G4Tubs("XrayWin", 0.0*mm, 6.4135*mm, 0.0508*mm, 0.*deg, 360.*deg);
    fXrayWin_LV = 
      new G4LogicalVolume(XrayWin,
                          G4Material::GetMaterial("SS304"),
                          "XrayWin_LV",0,0,0);

    G4VisAttributes *VisAtt_XrayWindow = 
      new G4VisAttributes(G4Colour(1.0,0.0,1.0));
    VisAtt_XrayWindow->SetVisibility(fVis2);
    fXrayWin_LV->SetVisAttributes(VisAtt_XrayWindow);

    fXrayWin = new G4PVPlacement(0,
          G4ThreeVector(0.*m, 0.*m, fSAD - 8.503*mm - fGantryPos),
          fXrayWin_LV, "XrayWindow", fGantry_LV, false, 0);

    fBeWinLow_LV->SetRegion(fTargetRegion);
    fTargetLow_LV->SetRegion(fTargetRegion);
    fTargetNicoroLow_LV->SetRegion(fTargetRegion);
    fTargetBlockTopLow_LV->SetRegion(fTargetRegion);
    fTargetRegion->AddRootLogicalVolume(fTargetLow_LV);
    fTargetRegion->AddRootLogicalVolume(fTargetNicoroLow_LV);
    fTargetRegion->AddRootLogicalVolume(fTargetBlockTopLow_LV);
    fTargetRegion->AddRootLogicalVolume(fBeWinLow_LV);


  }

  //============================================================================
  else if (fTargetName=="MediumEnergy") {
  //============================================================================
  // NB The 10XFFF mode uses the 15X (HighEnergy) target!
    //--------------------------------------------------------------------------
    // Target Block Top 
    //--------------------------------------------------------------------------

    const G4int tbt_10x_points = 9;
    G4double r_pc[tbt_10x_points] = {
      0.0*mm, 3.1*mm, 3.1   *mm, 8.0   *mm,
      8.0*mm, 5.5*mm, 4.0677*mm, 3.7875*mm, 0.0*mm};

    G4double z_pc[tbt_10x_points] = {
        .924*mm,   .924*mm,  1.124*mm,  1.124*mm,
      -8.876*mm, -8.876*mm, -5.369*mm, -4.776*mm, -4.776*mm};

    G4GenericPolycone* target_block_top= 
      new G4GenericPolycone("TargetBlockTop",
                     0.*deg, 360.*deg, tbt_10x_points, r_pc, z_pc);

    fTargetBlockTopMed_LV = 
      new G4LogicalVolume(target_block_top,
                          G4Material::GetMaterial("COPPER_GLIDCOP"),
                          "TargetBlockTop_LV", 0, 0, 0);

    fTargetBlockTopMed_LV->SetVisAttributes(target_block_VisAtt);

    fTargetBlockTopMed = new G4PVPlacement(0,
          G4ThreeVector(0., 0., fSAD - fGantryPos), fTargetBlockTopMed_LV, 
          "TargetBlockTopMed", fGantry_LV, false, 0);

    fTargetBlockTopMed_LV->SetVisAttributes(target_block_VisAtt);

    fTargetBlockTopMed_LV->SetRegion(fTargetRegion);
    fTargetRegion->AddRootLogicalVolume(fTargetBlockTopMed_LV);

  } 

  //============================================================================
  else if (fTargetName=="HighEnergy") { //15X and higher, and 10FFF
  //============================================================================
    //--------------------------------------------------------------------------
    //  Target Button
    //--------------------------------------------------------------------------

    G4double nicoroHThick = 0.0508*mm; // half the brazing sheet thickness

    const G4int he_pts = 5;
    G4double r_tb[he_pts] = 
      {0.0*mm, 3.0734*mm,  3.0734*mm,  2.8194*mm,  0.0  *mm};
    G4double z_tb[he_pts] = 
      {0.0*mm, 0.0   *mm, -0.381 *mm, -0.635 *mm, -0.635*mm};

    G4double Target_z = -0.859*mm;

    G4Polycone* target_button = 
      new G4Polycone("TargetButton",0.*deg, 360.*deg, he_pts, r_tb, z_tb);

    fTargetHigh_LV = 
      new G4LogicalVolume(target_button,
                          G4Material::GetMaterial("G4_W"),
                          "TargetButton_LV", 0, 0, 0);

    fTargetHigh = new G4PVPlacement(0,
          G4ThreeVector(0., 0., fSAD - Target_z + 2.*nicoroHThick - fGantryPos),
            fTargetHigh_LV, "TargetButtonHigh", 
            fGantry_LV, false, 0);

    fTargetHigh_LV->SetVisAttributes(target_VisAtt);

    //------------------------------------------------------------------------
    // nicoro brazing sheet
    //------------------------------------------------------------------------
    // not sure what the material should be: take it to be nicoro (BAu-3) 
    // (other choice is nicoro 80)

    G4Tubs* TargetNicoro =
      new G4Tubs("nicoro", 0.*mm, 3.0734*mm, nicoroHThick, 0.*deg, 360.*deg);

    fTargetNicoroHigh_LV = 
      new G4LogicalVolume(TargetNicoro,
                          G4Material::GetMaterial("Nicoro"),
                          "TargetNicoro_LV", 0, 0, 0);

    //new G4PVPlacement(0,G4ThreeVector(0, 0, fSAD - Target_z - 0.635*mm 
    fTargetNicoroHigh = 
      new G4PVPlacement(0, 
        G4ThreeVector(0., 0., fSAD - Target_z - 0.635*mm + nicoroHThick 
          - fGantryPos),
        fTargetNicoroHigh_LV, "NicoroHigh", fGantry_LV, 
        false, 0);

    G4VisAttributes* VisAtt_Nicoro = 
      new G4VisAttributes(G4Colour(0.0, 0.5, 0.5));
    fTargetNicoroHigh_LV->SetVisAttributes(VisAtt_Nicoro);

    //------------------------------------------------------------------------
    // Target Block Top - also use the 15 MV target for the 10XFFF beam.
    //------------------------------------------------------------------------

    //split the target holder as in the ProE drawing, 
    //into target block top and target block bottom

    const G4int tbt_he_pts = 6;
    G4double r_pc[tbt_he_pts] = 
      { 0.0  *mm, 3.1  *mm, 3.1  *mm, 7.0  *mm,  7.0  *mm,  0.0  *mm};
    G4double z_pc[tbt_he_pts] = 
      { 0.224*mm, 0.224*mm, 1.124*mm, 1.124*mm, -5.876*mm, -5.876*mm}; 

    G4GenericPolycone* target_block_top =
        new G4GenericPolycone("TargetBlockTopInit", 0.*deg, 360.*deg, 
                       tbt_he_pts, r_pc, z_pc);

    fTargetBlockTopHigh_LV = 
      new G4LogicalVolume(target_block_top,
                          G4Material::GetMaterial("COPPER_GLIDCOP"),
                          "TargetBlockTop_LV", 0, 0, 0);

    fTargetBlockTopHigh_LV->SetVisAttributes(target_block_VisAtt);

    fTargetBlockTopHigh = 
      new G4PVPlacement(0, G4ThreeVector(0., 0., fSAD - fGantryPos),
                     fTargetBlockTopHigh_LV,
                     "TargetBlockTopHigh", fGantry_LV, false, 0);


    //------------------------------------------------------------------------
    // braze wafer between target block top and target block bottom
    //------------------------------------------------------------------------
    // 35% Au & 65% Cu, 0.076 mm thick

    G4double waferThick = 0.076*mm;
    G4Tubs* target_wafer =
        new G4Tubs("TargetTopWafer", 0.*mm, 7.*mm, waferThick/2., 
                   0.*deg, 360.*deg);

    G4double waferPos = 5.876*mm; //using the target block shift
    fTargetWafer_LV = 
      new G4LogicalVolume(target_wafer,
                          G4Material::GetMaterial("CuAuAlloy"),
                          "TargetWafer_LV", 0, 0, 0);

    fTargetWafer = 
      new G4PVPlacement(0, 
          G4ThreeVector(0., 0., fSAD - waferPos - waferThick/2. - fGantryPos),
          fTargetWafer_LV, "brazeWafer", fGantry_LV, false, 0);

    G4VisAttributes* VisAtt_wafer = 
      new G4VisAttributes(G4Colour(0.2, 0.5, 0.5, 0.3));
    VisAtt_wafer->SetVisibility(fVis2);
    fTargetWafer_LV->SetVisAttributes(VisAtt_wafer);


    //------------------------------------------------------------------------
    //target block bottom
    //------------------------------------------------------------------------

    G4double bottom_thick = 2.924*mm;
    G4Tubs* target_bottom =
      new G4Tubs("TargetBlockBottom", 0.*mm, 7.0*mm,
                            bottom_thick/2, 0.*deg, 360.*deg);

    G4double target_bottom_pos = waferPos + waferThick;

    fTargetBlockBottomHigh_LV = 
      new G4LogicalVolume(target_bottom,
                          G4Material::GetMaterial("COPPER_GLIDCOP"),
                          "TargetBottom_LV", 0, 0, 0);

    //new G4PVPlacement(0,G4ThreeVector(0, 0, fSAD - target_bottom_pos 
    fTargetBlockBottomHigh = 
      new G4PVPlacement(0, 
          G4ThreeVector(0., 0., fSAD -target_bottom_pos - bottom_thick/2. 
            - fGantryPos),
          fTargetBlockBottomHigh_LV, "TargetBlockBottomHigh",
          fGantry_LV, false, 0);

    fTargetBlockBottomHigh_LV->SetVisAttributes(target_block_VisAtt);

    fTargetHigh_LV->SetRegion(fTargetRegion);
    fTargetNicoroHigh_LV->SetRegion(fTargetRegion);
    fTargetBlockTopHigh_LV->SetRegion(fTargetRegion);
    fTargetWafer_LV->SetRegion(fTargetRegion);
    fTargetBlockBottomHigh_LV->SetRegion(fTargetRegion);
    fTargetRegion->AddRootLogicalVolume(fTargetHigh_LV);
    fTargetRegion->AddRootLogicalVolume(fTargetNicoroHigh_LV);
    fTargetRegion->AddRootLogicalVolume(fTargetBlockTopHigh_LV);
    fTargetRegion->AddRootLogicalVolume(fTargetWafer_LV);
    fTargetRegion->AddRootLogicalVolume(fTargetBlockBottomHigh_LV);

  }

  //==========================================================================
  else if (fTargetName=="Imaging") {
  //==========================================================================
    //------------------------------------------------------------------------
    // Target Block Top for imaging target
    //------------------------------------------------------------------------
    // 2 mm Cu

    const G4int tbt_imaging_points = 6;
    G4double r_pc[tbt_imaging_points] = {
      0.0  *mm, 8.0  *mm, 8.0*mm,
      5.679*mm, 3.535*mm, 0.0*mm
    };

    G4double z_pc[tbt_imaging_points] = {
        0.*mm,  0.*mm, -10.*mm,
      -10.*mm, -2.*mm,  -2.*mm
    };

    G4GenericPolycone* target_block_top = 
      new G4GenericPolycone("TargetBlockTop",
                     0*deg, 360*deg, tbt_imaging_points, r_pc, z_pc);

    fTargetBlockTopImage_LV = 
      new G4LogicalVolume(target_block_top,
                          G4Material::GetMaterial("COPPER_GLIDCOP"),
                          "TargetBlockTop_LV", 0, 0, 0);

    fTargetBlockTopImage_LV->SetVisAttributes(target_block_VisAtt);

    fTargetBlockTopImage = new G4PVPlacement(0,
        G4ThreeVector(0., 0., fSAD - fGantryPos), fTargetBlockTopImage_LV,
        "TargetBlockTopImaging", fGantry_LV, false, 0);

    fTargetBlockTopImage_LV->SetRegion(fTargetRegion);
    fTargetRegion->AddRootLogicalVolume(fTargetBlockTopImage_LV);

  } else {
    G4ExceptionDescription ed;
    ed << "Target " << fTargetName << " undefined.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac102",
                FatalException, ed);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::BuildPrimaryCollimator(G4bool build) {
  //============================================================================
  //  Primary Collimator
  //============================================================================
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
  G4VPhysicalVolume* myVol = store->GetVolume("PrimColl", false);
  if (myVol) store->DeRegister(myVol);
  // TODO this should be part of BuildShielding
  myVol = store->GetVolume("carousel_bp", false);
  if (myVol) store->DeRegister(myVol);

  if (!build) return;

  G4double pc_pos        = 53.25*mm;

  G4RotationMatrix* PrimCollRot = new G4RotationMatrix();
  PrimCollRot->rotateX(180.0*deg);
  PrimCollRot->rotateZ(90.0*deg);

  const G4int pc_points = 8;
  G4double r_pcoll[pc_points] =  {
      8.5 *mm, 112.55 *mm,112.55*mm, 100.0   *mm,
    100.0 *mm,   26.0 *mm, 26.0 *mm,  21.4924*mm
  };
  G4double z_pcoll[pc_points] = {
    -26.25   *mm, -26.25   *mm,  -0.25*mm,  -0.25*mm,
     12.74999*mm,  12.74999*mm,  26.24*mm,  26.24*mm
  };
  // value above is 12.74999 rather than 12.75 to solve 
  // "GeomTest problem: solid problem" when running overlap check
  
  G4GenericPolycone* PCCone = 
    new G4GenericPolycone("PrimColl", 0.*deg, 360.*deg, pc_points, 
                          r_pcoll, z_pcoll);
  G4LogicalVolume* primColl_LV = 
    new G4LogicalVolume(PCCone,
                        G4Material::GetMaterial("W95"),
                        "PrimColl_LV", 0, 0, 0);

  G4VisAttributes* VisAtt_PC = 
    new G4VisAttributes(G4Colour(0.8, 0.1, 0.8, 0.5));
  VisAtt_PC->SetVisibility(fVis2);
  primColl_LV->SetVisAttributes(VisAtt_PC);
  
  fPrimColl = new G4PVPlacement(PrimCollRot, 
    G4ThreeVector(0., 0., fSAD - pc_pos - fGantryPos),
    primColl_LV, "PrimColl", fGantry_LV, false, 0);

  //G4Region* primColl_region = new G4Region("primcoll");
  //primColl_LV->SetRegion(primColl_region);
  //primColl_region->AddRootLogicalVolume(primColl_LV);

  // carousel baseplate
  G4Tubs* carousel_bp = 
    new G4Tubs("carousel_bp", 0.*mm, 175.*mm, 12.5*mm, 0., twopi);

  G4Box* carousel_bp_box = 
    new G4Box("carousel_bp_1", 175.*mm, 120.*mm, 12.5*mm);

  G4UnionSolid* carousel_bp_union =
    new G4UnionSolid("carousel_union", carousel_bp, carousel_bp_box, 0,
        G4ThreeVector(0., -110.*mm, 0.));

  G4Tubs* carousel_bp_2 =
    new G4Tubs("carousel_bp_2", 0.*mm, 112.55*mm, 13.5*mm, 0, twopi);

  G4SubtractionSolid* carousel_bp_subt = 
    new G4SubtractionSolid("carousel_subt", carousel_bp_union, carousel_bp_2,
      0, G4ThreeVector(0., 0., 0.));

  G4LogicalVolume* carousel_bp_LV =
    new G4LogicalVolume(carousel_bp_subt,
                        G4Material::GetMaterial("SS_A36"),
                        "carousel_bp_LV", 0, 0, 0);

  G4VisAttributes* VisAtt_carousel_bp = 
    new G4VisAttributes(G4Colour(0.3, 0.7, 0.1, 0.3));
  
  carousel_bp_LV->SetVisAttributes(VisAtt_carousel_bp);

  new G4PVPlacement(0, 
    G4ThreeVector(0., 0., 
      fSAD - pc_pos + 26.25*mm - 12.5*mm - fGantryPos),
    carousel_bp_LV, "carousel_bp", fGantry_LV, false, 0);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::BuildYStageShield(G4bool build) {
  //============================================================================
  // Y stage shield
  //============================================================================
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();

  G4VPhysicalVolume* myVol = store->GetVolume("ystageAll", false);
  if (myVol) store->DeRegister(myVol);

  if (!build) return;

  G4Cons* ystage_1 = 
    new G4Cons("ystage_1", 23.5*mm, 145.*mm, 36.4325*mm, 145.*mm, 
                               16.5*mm, 0.*deg, 360.*deg);
  //not supposed to have overlapping surfaces because they create
  // "false" surfaces
  G4double trap_delta = 0.1*mm;

  G4Tubs* ystage_4 =
    new G4Tubs("ystage_4", 26.5*mm, 145.*mm, 6.5*mm, 0.*deg, 360.*deg);

  G4UnionSolid* ystage_5 =
    new G4UnionSolid("ystage_5", ystage_1, ystage_4,
                      0, G4ThreeVector(0., 0., -23.*mm + trap_delta));

  G4Box* ystage_6 =
    new G4Box("ystage_6", 145.*mm, 140.*mm, 23.*mm);

  G4UnionSolid* ystage_7 =
    new G4UnionSolid("ystage_7", ystage_5, ystage_6, 0,
      G4ThreeVector(0., 170.*mm, - 6.5*mm ));

  // cutout groove as an extruded polygon
  G4double point1 = 2.5*mm;
  G4double point2 = 19.*mm;
  G4double point3 = 11.*mm;  // half length in beam direction

  G4TwoVector v1( point1,  point3);
  G4TwoVector v2( point2, -point3);
  G4TwoVector v3(-point2, -point3);
  G4TwoVector v4(-point1,  point3);
  std::vector<G4TwoVector> polygon {v1, v2, v3, v4};

  G4double l1 = 30.*mm;
  G4TwoVector vv1(0.*mm, 0.*mm);
  G4ExtrudedSolid::ZSection zs1(-l1, vv1, 1.);
  G4ExtrudedSolid::ZSection zs2( l1, vv1, 1.);
  std::vector<G4ExtrudedSolid::ZSection> zsections {zs1, zs2};
  G4ExtrudedSolid* ystage_2 = 
          new G4ExtrudedSolid("scCenter", polygon, zsections);

  //subtract several of the extruded solid, around an arc
  std::vector<G4RotationMatrix*> ystage_rot;
  std::vector<G4SubtractionSolid*>  subtsolids;
  G4double ext_rad = 100.*mm;  // radius of groove
  G4double ext_angle;
  G4double ext_dx, ext_dy;

  const G4int num_grooves = 7;  // can't visualize more than 7!
  for (G4int i = 0; i <= num_grooves; i += 1) {
    ext_angle = -90.*deg + 180.*deg * G4double(i)/G4double(num_grooves-1);
    ystage_rot.push_back(new G4RotationMatrix());
    ystage_rot[i]->rotateX(90.*deg);
    ystage_rot[i]->rotateY(90.*deg);
    ystage_rot[i]->rotateY(ext_angle);

    ext_dx = ext_rad * std::sin(ext_angle);
    ext_dy = ext_rad * (1. - std::cos(ext_angle));
    if (i == 0) {
      subtsolids.push_back(
        new G4SubtractionSolid("ystage_3a", ystage_7, ystage_2,
               ystage_rot[i],
               G4ThreeVector(ext_dx, ext_dy, 5.5*mm + trap_delta)));
    } else {
      subtsolids.push_back(
        new G4SubtractionSolid("ystage_3a", subtsolids[i-1], ystage_2,
               ystage_rot[i],
               G4ThreeVector(ext_dx, ext_dy, 5.5*mm + trap_delta)));
    }
  }
  
  G4LogicalVolume* ystageAll_LV = 
    new G4LogicalVolume(subtsolids[num_grooves-1],
                        G4Material::GetMaterial("Lead97Antimony"),
                        "ystageShield_LV", 0, 0, 0);

  G4VisAttributes* VisAtt_ystageAll = 
    new G4VisAttributes(G4Colour(0.2, 0.2, 0.6, 0.8));
  VisAtt_ystageAll->SetVisibility(fVis2);
  VisAtt_ystageAll->SetForceSolid(false);
  ystageAll_LV->SetVisAttributes(VisAtt_ystageAll);

  G4RotationMatrix* ystage2_Rot2 = new G4RotationMatrix();
  ystage2_Rot2->rotateX(180.*deg);

  new G4PVPlacement(ystage2_Rot2,
    G4ThreeVector(0., 0., fSAD - 96.5*mm - fGantryPos),
                    ystageAll_LV, "ystageAll", fGantry_LV, false, 0); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::BuildFlatteningFilter(G4bool build) {
  //============================================================================
  // Flattening filters / open port
  //============================================================================
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4String> physVols;
  physVols.push_back("carousel");
  physVols.push_back("FF_4X");
  physVols.push_back("FF_6X");
  physVols.push_back("FF_8X");
  physVols.push_back("FF_10X");
  physVols.push_back("FF_15X");
  physVols.push_back("FF_15X_Ta");
  physVols.push_back("FF_15X_plate");
  physVols.push_back("FF_15X_bottom");
  physVols.push_back("FF_18X");
  physVols.push_back("FF_18X_Ta");
  physVols.push_back("FF_18X_plate");
  physVols.push_back("FF_18X_bottom");
  physVols.push_back("FF_20X");
  physVols.push_back("FF_20X_Ta");
  physVols.push_back("FF_20X_plate");
  physVols.push_back("FF_20X_bottom");
  physVols.push_back("openPort");
  //physVols.push_back("FRMother");
  physVols.push_back("FrictionRing");
  std::vector<G4String>::iterator it;
  for (it = physVols.begin(); it != physVols.end(); ++it) {
    G4VPhysicalVolume* myVol = store->GetVolume(*it, false);
    if (myVol) store->DeRegister(myVol);
  }

  if (!build) return;

  //============================================================================
  //  Aluminum carousel (holds flattening filters)
  //============================================================================

  G4double z_FFNew = 125.7*mm; 

  const G4int npts_carousel = 6;
  G4double rInner_carousel[npts_carousel] = 
    {34.*mm, 34.*mm, 38.1*mm, 38.1*mm, 44.*mm, 44.*mm};
  G4double rOuter_carousel[npts_carousel] = 
    {200.*mm, 200.*mm, 200.*mm, 200.*mm, 200.*mm, 200.*mm};
  G4double zPlane_carousel[npts_carousel] = 
    { 8.9*mm, 0.*mm, 0.*mm, -1.5*mm, -1.5*mm, -7.*mm}; // TODO signs

  G4Polycone* carousel = 
    new G4Polycone("carousel", 0.*deg, 360.*deg, npts_carousel,
      zPlane_carousel, rInner_carousel, rOuter_carousel);

  G4LogicalVolume* carousel_LV = 
    new G4LogicalVolume(carousel,
                        G4Material::GetMaterial("Aluminum6061"),
                        "FFCarousel_LV", 0, 0, 0);

  new G4PVPlacement(0, 
                    G4ThreeVector(fFFOffset.x(), fFFOffset.y(),
                            fSAD - z_FFNew + fFFOffset.z() - fGantryPos),
                    carousel_LV, "carousel", fGantry_LV, false, 0);

  G4VisAttributes* VisAtt_carousel = 
    new G4VisAttributes(G4Colour(0.5, 0.4, 0.0, 0.2));
  VisAtt_carousel->SetVisibility(fVis2);
  carousel_LV->SetVisAttributes(VisAtt_carousel);

  //============================================================================
  //  Flattening filters
  //============================================================================

  G4GenericPolycone* FF_polycone;
  G4LogicalVolume* FF_LV;
  
  G4RotationMatrix* FF_rot = new G4RotationMatrix();
  FF_rot->rotateX(180.0*deg);

  G4VisAttributes* VisAtt_FF = 
    new G4VisAttributes(G4Colour(1.0,1.0,0.0,0.5));
  VisAtt_FF->SetVisibility(fVis2);

  //==========================================================================
  if (fFlatteningFilterName == "4X") {
  //==========================================================================
    const G4int cone1size = 24;

    G4double r_FF[cone1size] = {
       0.0 *mm,   0.635*mm,  1.27 *mm,  1.905*mm,
       2.54 *mm,  3.81 *mm,  5.08 *mm,  6.35 *mm,
       7.62 *mm,  8.89 *mm, 10.16 *mm, 12.7  *mm,
      15.24 *mm, 17.78 *mm, 20.32 *mm, 22.86 *mm,
      25.4  *mm, 27.94 *mm, 30.607*mm, 33.02 *mm,
      33.655*mm, 38.1  *mm, 38.1  *mm,  0.0  *mm
    };

    G4double z_FF[cone1size] = {
     -12.014*mm,-11.836*mm,-11.633*mm,-11.43 *mm,
     -11.201*mm,-10.668*mm, -9.982*mm, -9.246*mm,
      -8.56 *mm, -7.798*mm, -7.112*mm, -5.537*mm,
      -4.14 *mm, -2.616*mm, -1.473*mm, -0.33 *mm,
       0.813*mm,  1.956*mm,  2.032*mm,  2.032*mm,
       0.   *mm,  0.   *mm,  3.175*mm,  3.175*mm
    };

    FF_polycone = 
      new G4GenericPolycone("FFNew_PhysVol", 0.*deg, 360.*deg, cone1size, 
                            r_FF, z_FF);

    FF_LV = 
      new G4LogicalVolume(FF_polycone,
                          G4Material::GetMaterial("G4_Cu"),
                          "FFNew_LV",0,0,0);

    FF_LV->SetVisAttributes(VisAtt_FF);

    new G4PVPlacement(FF_rot,
                      G4ThreeVector(fFFOffset.x(),fFFOffset.y(),
                              fSAD - z_FFNew + fFFOffset.z() - fGantryPos),
                      FF_LV, "FF_4X", fGantry_LV, false, 0);
  }
  //==========================================================================
  else if (fFlatteningFilterName == "6X") {
  //==========================================================================
    const G4int cone1size = 24;

    G4double r_FF[cone1size] = {
       0.0  *mm,  0.635*mm,  1.27 *mm,  1.905*mm,
       2.54 *mm,  3.81 *mm,  5.08 *mm,  6.35 *mm,
       7.62 *mm,  8.89 *mm, 10.16 *mm, 12.7  *mm,
      15.24 *mm, 17.78 *mm, 20.32 *mm, 22.86 *mm,
      25.4  *mm, 27.94 *mm, 30.607*mm, 33.02 *mm,
      33.655*mm, 38.1  *mm, 38.1  *mm,  0.   *mm
    };

    G4double z_FF[cone1size] = {
     -18.999*mm,-18.72 *mm,-18.44 *mm,-18.059*mm,
     -17.653*mm,-16.916*mm,-15.57 *mm,-14.478*mm,
     -13.386*mm,-12.268*mm,-11.227*mm, -9.144*mm,
      -7.239*mm, -5.385*mm, -3.708*mm, -2.159*mm,
      -0.737*mm,  0.559*mm,  1.016*mm,  1.016*mm,
       0.   *mm,  0.   *mm,  3.175*mm,  3.175*mm
    };

    FF_polycone = 
      new G4GenericPolycone("FFNew_PhysVol", 0.*deg, 360.*deg, cone1size, 
                            r_FF, z_FF);

    FF_LV = 
      new G4LogicalVolume(FF_polycone,
                          G4Material::GetMaterial("G4_Cu"),
                          "FFNew_LV",0,0,0);

    FF_LV->SetVisAttributes(VisAtt_FF);

    new G4PVPlacement(FF_rot,
                      G4ThreeVector(fFFOffset.x(),fFFOffset.y(),
                            fSAD - z_FFNew + fFFOffset.z() - fGantryPos),
                      FF_LV, "FF_6X", fGantry_LV, false, 0);
  }
  //==========================================================================
  else if (fFlatteningFilterName=="8X") {
  //==========================================================================

    const G4int cone1size = 25;

    G4double r_FF[cone1size] = {
       0.   *mm,  0.610*mm,  1.219*mm,  1.829*mm,
       2.438*mm,  3.658*mm,  4.851*mm,  7.290*mm,
       9.728*mm, 12.167*mm, 14.580*mm, 17.018*mm,
      19.456*mm, 21.869*mm, 24.308*mm, 25.527*mm,
      26.746*mm, 27.965*mm, 29.185*mm, 30.378*mm,
      31.598*mm, 34.290*mm, 38.1  *mm, 38.1  *mm,
       0.   *mm
    };

    G4double z_FF[cone1size] = {
     -27.559*mm,-27.28 *mm,-26.772*mm,-26.162*mm,
     -25.375*mm,-23.571*mm,-21.666*mm,-18.009*mm,
     -14.605*mm,-11.557*mm, -8.788*mm, -6.248*mm,
      -3.912*mm, -1.778*mm,  0.0  *mm,  0.762*mm,
       1.524*mm,  2.083*mm,  2.464*mm,  2.565*mm,
       2.667*mm,  0.0  *mm,  0.0  *mm,  3.175*mm,
       3.175*mm
    };

    FF_polycone = 
      new G4GenericPolycone("FFNew_PhysVol", 0.*deg, 360.*deg, cone1size, 
                            r_FF, z_FF);

    FF_LV = new G4LogicalVolume(FF_polycone,
                                G4Material::GetMaterial("G4_Cu"),
                                "FFNew_LV", 0, 0, 0);

    FF_LV->SetVisAttributes(VisAtt_FF);

    new G4PVPlacement(FF_rot,  // TODO define sign of ffOffset.z
                      G4ThreeVector(fFFOffset.x(),fFFOffset.y(),
                            fSAD - z_FFNew + fFFOffset.z() - fGantryPos),
                      FF_LV,  "FF_8X", fGantry_LV, false, 0);
  }
  //==========================================================================
  else if (fFlatteningFilterName=="10X") {
  //==========================================================================

    const G4int cone1size = 34;

    G4double r_FF[cone1size] = {
       0.0 *mm,  0.635*mm,  1.6 *mm,  2.54 *mm,
       3.81*mm,  5.08 *mm,  7.62*mm, 10.16 *mm,
      12.7 *mm, 15.24 *mm, 17.78*mm, 20.32 *mm,
      22.86*mm, 25.4  *mm, 27.94*mm, 30.48 *mm,
      //38.1 *mm, 38.1  *mm, 31.75*mm, 29.464*mm,
      38.1 *mm, 38.1  *mm, 31.00*mm, 29.464*mm,
      27.94*mm, 25.4  *mm, 22.86*mm, 20.32 *mm,
      17.78*mm, 15.24 *mm, 12.7 *mm, 10.16 *mm,
       7.62*mm,  5.08 *mm,  3.81*mm,  2.54 *mm,
       1.27*mm,  0.0  *mm
    };

    G4double z_FF[cone1size] = {
      -32.131*mm, -32.131*mm,-30.912*mm, -29.515*mm,
      -27.381*mm, -25.349*mm,-21.742*mm, -18.618*mm,
      -15.748*mm, -12.954*mm,-10.084*mm,  -7.341*mm,
       -5.004*mm,  -3.226*mm, -1.575*mm,   0.0  *mm,
        0.0  *mm,   2.921*mm,  2.921*mm,   3.81 *mm,
        3.81 *mm,   3.81 *mm,  3.81 *mm,   3.988*mm,
        4.191*mm,   4.699*mm,  5.588*mm,   6.68 *mm,
        7.823*mm,   9.093*mm,  9.652*mm,  10.16 *mm,
       10.541*mm,  10.77 *mm
    };

    FF_polycone =
      new G4GenericPolycone("FF10X_PhysVol", 0.*deg, 360.*deg, cone1size, 
                            r_FF, z_FF);

    FF_LV = 
      new G4LogicalVolume(FF_polycone,
                          G4Material::GetMaterial("G4_Cu"),
                          "FF10X_LV", 0, 0, 0);

    FF_LV->SetVisAttributes(VisAtt_FF);

    new G4PVPlacement(FF_rot,
                      G4ThreeVector(fFFOffset.x(),fFFOffset.y(),
                            fSAD - z_FFNew + fFFOffset.z() - fGantryPos),
                      FF_LV, "FF_10X", fGantry_LV, false, 0);
  }

  //==========================================================================
  else if (fFlatteningFilterName=="15X") {
  //==========================================================================
    //-----------------------------------------------------
    // outer steel field flattener, top, 15x -checked Magda
    //-----------------------------------------------------
    const G4int cone1size = 36;

    G4double r_FF[cone1size] = {
       0.0  *mm,  1.27 *mm, 1.905  *mm,  2.54  *mm, //1-4
       3.175*mm,  3.81 *mm, 4.445  *mm,  5.08  *mm, //5-8
       6.35 *mm,  7.62 *mm, 9.144  *mm,  9.398 *mm, //9-12
      10.16 *mm, 10.541*mm, 11.43  *mm,  12.7  *mm, //13-16
      13.97 *mm, 15.24 *mm,                         //17,18
      16.51 *mm, 17.78 *mm,                         //19,20
      19.05 *mm, 20.32 *mm, 21.59  *mm, 22.86  *mm, //21-24
      24.13 *mm, 25.4  *mm, 26.67  *mm, 27.94  *mm, //25-28
      38.094*mm, 38.097*mm, 35.5345*mm, 35.5345*mm, //29,30,31,32
      13.335*mm,  0.36 *mm,  0.36  *mm,  0.0   *mm  //33,34,35,36
    };

    G4double z_FF[cone1size] = {
     -22.403*mm,-22.352 *mm,-22.276 *mm,-21.336*mm, //1-4
     -20.65 *mm,-20.244 *mm,-19.634 *mm,-19.101*mm, //5-8
     -18.39 *mm,-17.831 *mm,-17.018 *mm,-16.942*mm, //9-12
     -17.018*mm,-17.018 *mm,-15.748 *mm,-14.046*mm, //13-16
     -12.319*mm,-11.151 *mm,                        //17,18
      -9.957*mm, -8.89  *mm,                        //19,20
      -8.204*mm, -7.468 *mm, -5.893 *mm, -5.258*mm, //21-24
      -3.861*mm, -2.438 *mm, -0.711 *mm, -0.0  *mm, //25-28
       0.0  *mm,  1.0922*mm,  1.0922*mm,  2.235*mm, //29,30,31,32
       2.235*mm,-10.72  *mm,-11.1   *mm, -11.1 *mm  //33,34,35,36
    };

    FF_polycone = 
      new G4GenericPolycone("FF15X_PhysVol", 0.*deg, 360.*deg, cone1size, 
                            r_FF, z_FF);

    FF_LV = 
      new G4LogicalVolume(FF_polycone,
                          G4Material::GetMaterial("SS12L14"),
                          "FF15X_LV", 0, 0, 0);

    FF_LV->SetVisAttributes(VisAtt_FF);

    new G4PVPlacement(FF_rot,
                      G4ThreeVector(fFFOffset.x(),fFFOffset.y(),
                              fSAD - z_FFNew + fFFOffset.z() - fGantryPos),
                      FF_LV, "FF_15X", fGantry_LV, false, 0);

    //------------------------------------------------------------------------
    // tantalum conical insert - proE shows that it has a smaller diameter
    //------------------------------------------------------------------------

    const G4int cone2size = 3;
    // changed the table value of 13.335 to ProE measured value of 13.1616 mm
    G4double rout2_15XFF[cone2size] = {0.0*mm, 13.1616*mm, 0.0*mm};
    // conical insert has a smaller height than 11.1 mm
    G4double z2_15XFF[cone2size] = {-11.0579*mm, 2.235*mm, 2.235*mm};

    G4Polycone* FF15X_insert = 
      new G4Polycone("FF15X2_PhysVol", 0.*deg, 360.*deg, cone2size, 
                     rout2_15XFF, z2_15XFF);

    G4LogicalVolume* FF15X_insert_LV = 
      new G4LogicalVolume(FF15X_insert,
                          G4Material::GetMaterial("G4_Ta"),
                          "FF15X_2_LV", 0, 0, 0);

    G4VisAttributes* VisAtt_FF15X_insert = 
      new G4VisAttributes(G4Colour(0.5,0.0,0.0,0.5));
    VisAtt_FF15X_insert->SetVisibility(fVis2);
    FF15X_insert_LV->SetVisAttributes(VisAtt_FF15X_insert);

    new G4PVPlacement(FF_rot,
                      G4ThreeVector(fFFOffset.x(),fFFOffset.y(),
                            fSAD - z_FFNew + fFFOffset.z() - fGantryPos),
                      FF15X_insert_LV, "FF_15X_Ta", fGantry_LV, 
                      false, 0);

    //------------------------------------------------------------------------
    // Steel plate - this plate has a deformation in the middle,
    // I don't know if it's relevant, not yet modeled
    //------------------------------------------------------------------------

    G4double FF15X_steel_hthick = 0.1905*mm;
    G4Tubs* FF15X_steel = 
      new G4Tubs("FF15Xsteel", 0.*mm, 34.9247*mm, FF15X_steel_hthick,
                 0.*deg, 360.*deg);

    G4LogicalVolume* FF15X_steel_LV = 
      new G4LogicalVolume(FF15X_steel,
                          G4Material::GetMaterial("SS12L14"),
                          "FF15Xsteel_LV", 0, 0, 0);

    G4VisAttributes* VisAtt_FF15X_steel = 
      new G4VisAttributes(G4Colour(0.2, 0.2, 0.7, 1.0));
    VisAtt_FF15X_steel->SetVisibility(fVis2);
    FF15X_steel_LV->SetVisAttributes(VisAtt_FF15X_steel);

    new G4PVPlacement(FF_rot,
                G4ThreeVector(fFFOffset.x(),fFFOffset.y(),
                    fSAD - z_FFNew - 2.997*mm + 
                    FF15X_steel_hthick + fFFOffset.z() - fGantryPos),
                FF15X_steel_LV, "FF_15X_plate", fGantry_LV, false, 0);

    //------------------------------------------------------------------------
    // Steel bottom - FIELD FLATTENER, BOTTOM 18X
    // nb we are still in the 15X flattener
    //------------------------------------------------------------------------

    const G4int cone3size = 35;
    G4double rout3_15XFF[cone3size] = {
      0.0   *mm, 35.5345*mm, 35.5345*mm, 38.094*mm,   //xx,yy,zz,ww
      38.094*mm, 34.798 *mm, 34.29  *mm, 27.9  *mm,   //vv,33
      25.15 *mm, 23.927 *mm, 22.708 *mm, 21.488*mm,   //34-37
      20.295*mm, 19.075 *mm, 17.882 *mm,              //38-40
      16.662*mm, 15.469 *mm, 14.732 *mm,              //41-43
      14.275*mm, 13.97  *mm, 13.716 *mm,              //44-46
      13.462*mm, 13.081 *mm, 11.887 *mm,              //47-49
      10.693*mm,  9.5   *mm, 8.306  *mm,              //50-52
      7.112 *mm,  5.918 *mm, 4.75   *mm,              //53-55
      3.556 *mm,  2.362 *mm, 1.194  *mm,              //56-59
      0.    *mm  
    };

    G4double z3_15XFF[cone3size] = {
       2.997 *mm,  2.997 *mm,  1.0922*mm,  1.0922*mm,
       5.029 *mm,  5.029 *mm,  4.536 *mm,  4.536 *mm,
       4.79  *mm,  5.199 *mm,  5.913 *mm,  6.853 *mm,
       7.894 *mm,  8.987 *mm,  9.924 *mm,
       10.635*mm,  11.14 *mm,  11.394*mm,
       11.547*mm,  11.25 *mm,  11.044*mm,
       10.914*mm,  10.762*mm,  10.556*mm,
       10.404*mm,  10.251*mm,  10.251*mm,
       10.328*mm,  10.505*mm,  10.683*mm,
       10.861*mm,  11.14 *mm,  11.42 *mm,
       11.75*mm 
    };

    G4GenericPolycone* FF15X_bottom = 
      new G4GenericPolycone("FF15X_bottom",
                     0.*deg, 360.*deg, cone3size, rout3_15XFF, z3_15XFF);

    G4LogicalVolume* FF15X_bottom_LV = 
      new G4LogicalVolume(FF15X_bottom,
                          G4Material::GetMaterial("SS12L14"),
                          "FF15X3_LV", 0, 0, 0);

    G4VisAttributes* VisAtt_FF15X_bottom = 
      new G4VisAttributes(G4Colour(0.0,0.0,1.0, 0.5));
    VisAtt_FF15X_bottom->SetVisibility(fVis2);
    FF15X_bottom_LV->SetVisAttributes(VisAtt_FF15X_bottom);

    new G4PVPlacement(FF_rot,
                      G4ThreeVector(fFFOffset.x(),fFFOffset.y(),
                              fSAD - z_FFNew + fFFOffset.z() - fGantryPos),
                      FF15X_bottom_LV, "FF_15X_bottom", fGantry_LV, 
                      false, 0);
  }

  //==========================================================================
  else if (fFlatteningFilterName=="18X") {
  //==========================================================================

    G4double FF18XTaPositionShift = -0.148*mm;
    G4double FF18XbottomShift     = -0.021*mm;
    //-----------------------------------------------------
    // outer steel field flattener, top, 18x
    //-----------------------------------------------------
    const G4int cone1size = 52;

    G4double r_FF[cone1size] = {
       0.0    *mm, 0.305  *mm, 0.635 *mm, 0.838 *mm,   //1-4
       1.02   *mm, 1.27   *mm, 1.473 *mm, 1.727 *mm,
       1.905  *mm, 2.286  *mm, 2.5   *mm, 3.175 *mm,
       3.81   *mm, 4.445  *mm, 5.1   *mm, 6.35  *mm,
       7.6    *mm, 8.89   *mm, 9.144 *mm, 9.398 *mm,
       9.652  *mm, 9.906  *mm, 10.2  *mm, 10.41 *mm,
      10.67   *mm, 10.92  *mm, 11.43 *mm, 12.065*mm,
      13.97   *mm, 15.2   *mm, 16.51 *mm, 17.8  *mm,
      19.05   *mm, 20.3   *mm, 21.59 *mm, 22.9  *mm,
      24.13   *mm, 25     *mm, 26.67 *mm, 27.9  *mm,
      38.094  *mm, 38.094 *mm,
      35.5345 *mm, 35.5345*mm,
      13.335  *mm, 0.38   *mm, 0.38  *mm, 0.    *mm 
    };

    G4double z_FF[cone1size] = {
     -32.304 *mm,-32.253*mm,-32.117*mm,-32.052*mm,
     -31.902 *mm,-31.674*mm,-31.473*mm,-31.168*mm,
     -30.919 *mm,-29.967*mm,-29.309*mm,-27.981*mm,
     -26.982 *mm,-26.142*mm,-25.428*mm,-24.056*mm,
     -22.885 *mm,-21.946*mm,-21.765*mm,-21.664*mm,
     -21.562 *mm,-21.46 *mm,-21.186*mm,-20.907*mm,
     -20.627 *mm,-20.155*mm,-19.215*mm,-18.194*mm,
     -15.197 *mm,-13.619*mm,-12.22 *mm,-10.95 *mm,
      -9.754 *mm,-8.435 *mm,-7.038 *mm,-5.591 *mm,
      -4.092 *mm,-2.593 *mm,-1.095 *mm, 0.    *mm,
       0.    *mm, 1.0922*mm-FF18XbottomShift,
       1.0922*mm -FF18XbottomShift,  2.383*mm,
       2.383 *mm,-10.752*mm,-10.952*mm,-10.952*mm 
    };

    FF_polycone = 
      new G4GenericPolycone("FF18X", 0*deg, 360*deg, cone1size, r_FF, z_FF);

    FF_LV = 
      new G4LogicalVolume(FF_polycone,
                          G4Material::GetMaterial("SS12L14"),
                          "FF18X_LV", 0, 0, 0);

    FF_LV->SetVisAttributes(VisAtt_FF);

    new G4PVPlacement(FF_rot,
                      G4ThreeVector(fFFOffset.x(),fFFOffset.y(),
                              fSAD - z_FFNew + fFFOffset.z() - fGantryPos),
                      FF_LV, "FF_18X", fGantry_LV, false, 0);

    //------------------------------------------------------------------------
    // tantalum conical insert - proE shows that it has a smaller diameter
    // this is the same as for 15X, except for the positioning
    //------------------------------------------------------------------------

    const G4int cone2size = 3;
    // changed the table value of 13.335 to ProE measured value of 13.1616 mm
    G4double rout2_15XFF[cone2size] = {0.0*mm, 13.1616*mm, 0.0*mm};
    // conical insert has a smaller height than 11.1 mm
    G4double z2_15XFF[cone2size] = {-11.0579*mm,  2.235*mm,  2.235*mm};

    G4GenericPolycone* FF15X_insert = 
      new G4GenericPolycone("FF15X2_PhysVol",
                            0.*deg, 360.*deg, cone2size, rout2_15XFF, 
                            z2_15XFF);

    G4LogicalVolume* FF15X_insert_LV = 
      new G4LogicalVolume(FF15X_insert,
                          G4Material::GetMaterial("G4_Ta"),
                          "FF15X2_LV", 0, 0, 0);

    G4VisAttributes* VisAtt_FF15X_insert = 
      new G4VisAttributes(G4Colour(0.5, 0.0, 0.0, 0.5));
    VisAtt_FF15X_insert->SetVisibility(fVis2);
    FF15X_insert_LV->SetVisAttributes(VisAtt_FF15X_insert);

    new G4PVPlacement(FF_rot,
                      G4ThreeVector(fFFOffset.x(),fFFOffset.y(),
                                fSAD - z_FFNew + fFFOffset.z()
                                  + FF18XTaPositionShift - fGantryPos),
                                // positive sign on FF18.. because of rotation
                      FF15X_insert_LV, "FF_18X_Ta", fGantry_LV, 
                      false, 0);

    //------------------------------------------------------------------------
    // Steel plate - this plate has a deformation in the middle,
    // I don't know if it's relevant, not yet modeled
    // Same as for 18X as for 15X, but the position is different
    //------------------------------------------------------------------------

    G4double FF15X_steel_hthick = 0.1905*mm;
    G4Tubs* FF15X_steel = 
      new G4Tubs("FF15Xsteel", 0.*mm, 34.9247*mm, FF15X_steel_hthick, 
                 0.*deg, 360.*deg);

    G4LogicalVolume* FF15X_steel_LV = 
      new G4LogicalVolume(FF15X_steel,
                          G4Material::GetMaterial("SS12L14"),
                          "FF15Xsteel_LV",0,0,0);

    G4VisAttributes* VisAtt_FF15X_steel = 
      new G4VisAttributes(G4Colour(0.2,0.2,0.7,1.0));
    VisAtt_FF15X_steel->SetVisibility(fVis2);
    FF15X_steel_LV->SetVisAttributes(VisAtt_FF15X_steel);

    new G4PVPlacement(FF_rot,
                      G4ThreeVector(fFFOffset.x(),fFFOffset.y(),
                            fSAD - z_FFNew - 2.997*mm 
                              + FF15X_steel_hthick 
                          + fFFOffset.z() + FF18XbottomShift - fGantryPos),
                  FF15X_steel_LV, "FF_18X_plate", fGantry_LV, false, 0);

    //------------------------------------------------------------------------
    // Steel bottom - FIELD FLATTENER, BOTTOM 18X
    // (same for 18X as for 15X, but the position is different)
    //------------------------------------------------------------------------

    const G4int cone3size = 35;
    G4double rout3_15XFF[cone3size] = {
       0.0  *mm, 35.5345*mm, 35.5345*mm, 38.094*mm,   //xx,yy,zz,ww
      38.094*mm, 34.798 *mm, 34.29  *mm, 27.9  *mm,   //vv,33
      25.15 *mm, 23.927 *mm, 22.708 *mm, 21.488*mm,   //34-37
      20.295*mm, 19.075 *mm, 17.882 *mm,              //38-40
      16.662*mm, 15.469 *mm, 14.732 *mm,              //41-43
      14.275*mm, 13.97  *mm, 13.716 *mm,              //44-46
      13.462*mm, 13.081 *mm, 11.887 *mm,              //47-49
      10.693*mm,  9.5   *mm, 8.306  *mm,              //50-52
       7.112*mm,  5.918 *mm, 4.75   *mm,              //53-55
       3.556*mm,  2.362 *mm, 1.194  *mm, 0.    *mm    //56-59
    };

    G4double z3_15XFF[cone3size] = {
        2.997*mm,  2.997 *mm,  1.0922*mm,  1.0922*mm,
        5.029*mm,  5.029 *mm,  4.536 *mm,  4.536 *mm,
        4.79 *mm,  5.199 *mm,  5.913 *mm,  6.853 *mm,
        7.894*mm,  8.987 *mm,  9.924 *mm,
       10.635*mm,  11.14 *mm,  11.394*mm,
       11.547*mm,  11.25 *mm,  11.044*mm,
       10.914*mm,  10.762*mm,  10.556*mm,
       10.404*mm,  10.251*mm,  10.251*mm,
       10.328*mm,  10.505*mm,  10.683*mm,
       10.861*mm,  11.14 *mm,  11.42 *mm,  11.75 *mm 
    };

    G4GenericPolycone* FF15X_bottom = 
      new G4GenericPolycone("FF15X_bottom",
                            0*deg, 360*deg, cone3size, rout3_15XFF, 
                            z3_15XFF);

    G4LogicalVolume* FF15X_bottom_LV = 
      new G4LogicalVolume(FF15X_bottom,
                          G4Material::GetMaterial("SS12L14"),
                          "FF15X3_LV", 0, 0, 0);

    G4VisAttributes* VisAtt_FF15X_bottom = 
      new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.5));
    VisAtt_FF15X_bottom->SetVisibility(fVis2);
    FF15X_bottom_LV->SetVisAttributes(VisAtt_FF15X_bottom);

    new G4PVPlacement(FF_rot,
                      G4ThreeVector(fFFOffset.x(),fFFOffset.y(),
                              fSAD - z_FFNew + fFFOffset.z() 
                                + FF18XbottomShift - fGantryPos),
                      FF15X_bottom_LV, "FF_18X_bottom", fGantry_LV, 
                      false, 0);
  }

  //==========================================================================
  else if (fFlatteningFilterName=="20X") {
  //==========================================================================

    G4double FF20XbottomShift   = -0.015*mm;
    //------------------------------------------------------------------------
    // outer steel field flattener, top, 20x
    //------------------------------------------------------------------------
    const G4int cone1size = 51;

    G4double r_FF[cone1size] = {
       0.0  *mm, 0.305 *mm, 0.635 *mm, 0.838 *mm,   //1
       1.016*mm, 1.27  *mm, 1.473 *mm, 1.727 *mm,
       1.905*mm, 2.286 *mm, 2.54  *mm, 3.175 *mm,
       3.81 *mm, 4.445 *mm, 5.08  *mm, 6.35  *mm,
       7.62 *mm, 8.89  *mm, 9.144 *mm, 9.398 *mm,
       9.652*mm, 9.906 *mm, 10.16 *mm, 10.414*mm,
      10.668*mm, 10.922*mm, 11.43 *mm, 12.065*mm,
      12.7  *mm, 13.97 *mm, 15.24 *mm, 15.748*mm,
      16.51 *mm, 17.78 *mm, 19.05 *mm, 20.32 *mm,
      21.59 *mm, 22.86 *mm, 24.13 *mm, 25.4  *mm,   //10
      30.48 *mm, 30.48 *mm, 38.1  *mm, 38.1  *mm,
      35.535*mm, 35.535*mm, 13.335*mm,  0.381*mm,
       0.381*mm,  0.   *mm,  0.   *mm 
    };

    // these include the nickel plating
    G4double z_FF[cone1size] = {
      -31.361*mm, -31.158*mm, -30.808*mm,-30.508*mm,
      -30.229*mm, -29.725*mm, -29.319*mm,-28.740*mm,
      -28.379*mm, -27.178*mm, -26.558*mm,-25.182*mm,
      -24.059*mm, -22.868*mm, -21.773*mm,-19.941*mm,
      -18.841*mm, -17.976*mm, -17.750*mm,-17.546*mm,
      -17.490*mm, -17.490*mm, -17.414*mm,-17.168*mm,
      -16.962*mm, -16.807*mm, -16.055*mm,-15.189*mm,
      -14.051*mm, -11.440*mm, -10.005*mm, -9.520*mm,
       -8.585*mm,  -7.010*mm,  -5.586*mm, -4.239*mm,
       -2.918*mm,  -1.547*mm,  -0.075*mm,  1.575*mm,
        1.585*mm,   0.0  *mm,   0.0  *mm,  1.092*mm,
        1.092*mm,   2.362*mm,   2.362*mm,-11.333*mm,
      -11.735*mm, -11.735*mm, -31.361*mm 
    };

    FF_polycone = 
      new G4GenericPolycone("FF20X_PhysVol", 0.*deg, 360.*deg, cone1size, 
                            r_FF, z_FF);

    FF_LV = 
      new G4LogicalVolume(FF_polycone,
                          G4Material::GetMaterial("SS12L14"),
                          "FF20X_LV", 0, 0, 0);

    FF_LV->SetVisAttributes(VisAtt_FF);

    new G4PVPlacement(FF_rot,
                G4ThreeVector(fFFOffset.x(),fFFOffset.y(),
                              fSAD - z_FFNew + fFFOffset.z() - fGantryPos),
                FF_LV, "FF_20X", fGantry_LV, false, 0);

    //--------------------------------------------------------------------------
    // tantalum conical insert
    // height is different from 15 and 18X
    //--------------------------------------------------------------------------

    const G4int cone2size = 3;
    G4double rout2_15XFF[cone2size] = { 0.0   *mm, 13.334*mm,  0.   *mm};
    G4double z2_15XFF   [cone2size] = {-11.735*mm,  2.362*mm,  2.362*mm};

    G4Polycone* FF15X_insert = 
      new G4Polycone("FF15X_insert",
              0.*deg, 360.*deg, cone2size, rout2_15XFF, z2_15XFF);

    G4LogicalVolume* FF15X_insert_LV = 
      new G4LogicalVolume(FF15X_insert,
                          G4Material::GetMaterial("G4_Ta"),
                          "FF15X2_LV", 0, 0, 0);

    G4VisAttributes* VisAtt2_FF15X = 
      new G4VisAttributes(G4Colour(0.5,0.0,0.0,0.5));
    VisAtt2_FF15X->SetVisibility(fVis2);
    FF15X_insert_LV->SetVisAttributes(VisAtt2_FF15X);

    new G4PVPlacement(FF_rot,
                G4ThreeVector(fFFOffset.x(),fFFOffset.y(),
                              fSAD - z_FFNew + fFFOffset.z() - fGantryPos),
                        // negative sign on FF20.. because of rotation
                FF15X_insert_LV, "FF_20X_Ta", fGantry_LV, 
                false, 0);

    //------------------------------------------------------------------------
    // Steel plate - this plate has a deformation in the middle,
    // I don't know if it's relevant, not yet modeled
    // Same as for 18X and 15X, but the position is different
    //------------------------------------------------------------------------

    G4double FF15X_steel_hthick = 0.1905*mm;
    G4Tubs* FF15X_steel =
      new G4Tubs("FF15Xsteel", 0.*mm, 34.9247*mm, FF15X_steel_hthick,
                 0.*deg,360.*deg);

    G4LogicalVolume* FF15X_steel_LV = 
      new G4LogicalVolume(FF15X_steel,
                          G4Material::GetMaterial("SS12L14"),
                          "FF15Xsteel_LV", 0, 0, 0);

    G4VisAttributes* VisAtt_FF15X_steel = 
      new G4VisAttributes(G4Colour(0.2, 0.2, 0.7, 1.0));
    VisAtt_FF15X_steel->SetVisibility(fVis2);
    FF15X_steel_LV->SetVisAttributes(VisAtt_FF15X_steel);

    new G4PVPlacement(FF_rot,
                      G4ThreeVector(fFFOffset.x(),fFFOffset.y(),
                                    fSAD - z_FFNew - 2.997*mm 
                                      + FF15X_steel_hthick + fFFOffset.z()
                                      + FF20XbottomShift - fGantryPos),
                  FF15X_steel_LV, "FF_20X_plate", fGantry_LV, false, 0);

    //------------------------------------------------------------------------
    // Steel bottom - FIELD FLATTENER, BOTTOM 20X
    //------------------------------------------------------------------------

    const G4int cone3size = 35;
    G4double rout3_20XFF[cone3size] = {
       0.0  *mm, 35.5345*mm, 35.5345*mm, 38.094*mm,   //xx,yy,zz,ww
      38.094*mm, 34.798 *mm, 34.29  *mm, 27.9  *mm,   //vv,33
      25.15 *mm, 23.927 *mm, 22.708 *mm, 21.488*mm,   //34-37
      20.295*mm, 19.075 *mm, 17.882 *mm,              //38-40
      16.662*mm, 15.469 *mm, 14.732 *mm,              //41-43
      14.275*mm, 13.97  *mm, 13.716 *mm,              //44-46
      13.462*mm, 13.081 *mm, 11.887 *mm,              //47-49
      10.693*mm, 9.5    *mm, 8.306  *mm,              //50-52
       7.112*mm, 5.918  *mm, 4.75   *mm,              //53-55
       3.556*mm, 2.362  *mm, 1.194  *mm, 0.    *mm    //56-59
    };

    G4double z3_20XFF[cone3size] = {
        2.997*mm,  2.997 *mm,  1.0922*mm,  1.0922*mm,
        5.029*mm,  5.029 *mm,  4.536 *mm,  4.536 *mm,
        4.79 *mm,  5.199 *mm,  5.913 *mm,  6.853 *mm,
        7.894*mm,  8.987 *mm,  9.924 *mm,
       10.635*mm,  11.14 *mm,  11.394*mm,
       11.547*mm,  11.25 *mm,  11.044*mm,
       10.914*mm,  10.762*mm,  10.556*mm,
       10.404*mm,  10.251*mm,  10.251*mm,
       10.328*mm,  10.505*mm,  10.683*mm,
       10.861*mm,  11.14 *mm,  11.42 *mm,  11.75 *mm 
    };

    G4GenericPolycone* FF20X_bottom = 
      new G4GenericPolycone("FF20X_bottom",
                            0*deg, 360*deg, cone3size, rout3_20XFF, 
                            z3_20XFF);

    G4LogicalVolume* FF20X_bottom_LV = 
      new G4LogicalVolume(FF20X_bottom,
                          G4Material::GetMaterial("SS12L14"),
                          "FF20X3_LV", 0, 0, 0);

    G4VisAttributes* VisAtt_FF20X_bottom = 
      new G4VisAttributes(G4Colour(0.0,0.0,1.0, 0.5));
    VisAtt_FF20X_bottom->SetVisibility(fVis2);
    FF20X_bottom_LV->SetVisAttributes(VisAtt_FF20X_bottom);

    new G4PVPlacement(FF_rot,
                      G4ThreeVector(fFFOffset.x(),fFFOffset.y(),
                                    fSAD - z_FFNew + fFFOffset.z()
                                      + FF20XbottomShift),
                      FF20X_bottom_LV, "FF_20X_bottom", fLogicWorld, 
                      false, 0);
  }


  //==========================================================================
  else if (fFlatteningFilterName == "open" ||
           fFlatteningFilterName == "open port") { // open port (FFF)
  //==========================================================================
    
    G4double openPortHThick = 0.405*mm;
    G4Tubs* openPort = 
      new G4Tubs("openPort", 0.*mm, 42.2*mm, openPortHThick, 0.*deg, 360.*deg);

    G4LogicalVolume* openPort_LV = 
      new G4LogicalVolume(openPort,
                          G4Material::GetMaterial("Brass"),
                          "openPort_LV",0,0,0);

    G4VisAttributes* openPort_VisAtt = 
      new G4VisAttributes(G4Colour(0.1,1.0,0.1));
    openPort_VisAtt->SetVisibility(fVis2);
    openPort_VisAtt->SetForceSolid(false);
    openPort_LV->SetVisAttributes(openPort_VisAtt);
    
    new G4PVPlacement(0,
                  G4ThreeVector(fFFOffset.x(),fFFOffset.y(),
                          fSAD - 128.01*mm + fFFOffset.z() - fGantryPos),
                  openPort_LV,"openPort", fGantry_LV, false, 0);
  }

  else if (fFlatteningFilterName.compare("none") == 0 ||
           fFlatteningFilterName.compare("None") == 0) {
   ; // do nothing!
  }
  else {

    G4ExceptionDescription ed;
    ed << "Flattening filter " << fFlatteningFilterName << " undefined!";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac103",
                FatalException, ed);
  }

  //==========================================================================
  // Friction Ring
  //==========================================================================
  G4double fr_zhthick = 2.375*mm;
  //G4double FR_zRef = + fr_zhthick;     
  //from Fastrad GDML file (must reverse the sign)
  //ATTENTION: the old value 128.75mm makes the ring collapse w/ 15x filter
  G4double FR_z = 0.*mm;
  // the 'zero' of the flattening filter is 125.7 mm; 
  // FR_z = 125.7 + the thickness of the edge of the FF
  // ie FR_z is the top (close to target) side of friction ring
  G4bool useFrictionRing = true;
  if      (fFlatteningFilterName == "4X" ) FR_z = 128.875*mm;
  else if (fFlatteningFilterName == "6X" ) FR_z = 128.875*mm;
  else if (fFlatteningFilterName == "8X" ) FR_z = 128.875*mm;
  else if (fFlatteningFilterName == "10X") FR_z = 128.621*mm;
  else if (fFlatteningFilterName == "15X") FR_z = 130.729*mm;
  else if (fFlatteningFilterName == "18X") FR_z = 130.75 *mm;
  else if (fFlatteningFilterName == "20X") FR_z = 130.744*mm;
  else if (fFlatteningFilterName == "open" || 
           fFlatteningFilterName == "open port")
          {useFrictionRing = false;}
  else if (fFlatteningFilterName.compare("none") == 0 ||
           fFlatteningFilterName.compare("None") == 0) {
      useFrictionRing = false;
  }
  else { 

    G4ExceptionDescription ed;
    ed << "Error for friction ring: undefined flattening filter name "
           << fFlatteningFilterName;
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac104",
                FatalException, ed);
  }

  if (useFrictionRing) {
    
    G4Tubs* FrictionRing1 = 
      new G4Tubs("Friction ring_1", 31.0*mm, 43.5*mm, fr_zhthick, 
                 0.*deg, 360.*deg);

    G4Box* FrictionRing2 = 
      new G4Box("Friction_ring_notch", 3.*mm, 4.*mm, fr_zhthick + 1.*mm);

    G4SubtractionSolid* FrictionRing = 
      new G4SubtractionSolid("Friction_ring", FrictionRing1, FrictionRing2,
                             0, G4ThreeVector(0., -33.5*mm, 0.));

    G4LogicalVolume* FrictionRing_LV = 
      new G4LogicalVolume(FrictionRing,
                          G4Material::GetMaterial("SS_A36"),
                          "FrictionRing_LV", 0, 0, 0);

    new G4PVPlacement(0,
                    G4ThreeVector(fFFOffset.x(), fFFOffset.y(), 
                    fSAD - fr_zhthick - FR_z + fFFOffset.z() - fGantryPos),
                FrictionRing_LV, "FrictionRing", fGantry_LV,
                false, 0);

    G4VisAttributes *VisAtt_FrictionRing = 
      new G4VisAttributes(G4Colour(0.7, 0.0, 0.0, 0.5));
    VisAtt_FrictionRing->SetVisibility(fVis2);
    VisAtt_FrictionRing->SetForceSolid(false);
    FrictionRing_LV->SetVisAttributes(VisAtt_FrictionRing);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::BuildFoil1(G4bool build) {
  ////========================================================================
  //// foil 1 
  ////========================================================================

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  
  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
  G4VPhysicalVolume* myVol = store->GetVolume("electronFoil1", false);
  if (myVol) store->DeRegister(myVol);
  myVol = store->GetVolume("surround", false);
  if (myVol) store->DeRegister(myVol);
  
  G4LogicalVolumeStore* lv_store = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* lv = lv_store->GetVolume("efoil1_LV", false);
  if (lv) lv_store->DeRegister(lv);

  if (!build) return;

  G4double foil1Thick = -1.;
  G4Material* material = nullptr;
  G4double foil1PositionBottom = 62.7*mm;

  G4VisAttributes* efoil1_VisAtt;
  efoil1_VisAtt = new G4VisAttributes(G4Colour(0.0, 0.6, 0.5, 0.8));

  if (fFoil1Name != "None") {
    G4double foil1Radius = 6.*mm;
    if (fFoil1Name == "6E" || fFoil1Name == "9E") {
      foil1Thick = 0.23*mm * fFoil1ThicknessFactor;
      material = G4Material::GetMaterial("Brass");
    }
    else if (fFoil1Name == "12E") {
      // thickness given as 0.0813 +0.0051 -0 mm
      // this is the only primary foil that has different + and -
      foil1Thick = 0.0838*mm * fFoil1ThicknessFactor;
      material = G4Material::GetMaterial("G4_Ta");
    }
    else if (fFoil1Name == "15E") {
      foil1Thick = 0.136*mm * fFoil1ThicknessFactor;
      material = G4Material::GetMaterial("G4_Ta");
    }
    else if (fFoil1Name == "16E") {
      foil1Thick = 0.1565*mm * fFoil1ThicknessFactor;
      material = G4Material::GetMaterial("G4_Ta");
    }
    else if (fFoil1Name == "18E") {
      foil1Thick = 0.189*mm * fFoil1ThicknessFactor;
      material = G4Material::GetMaterial("G4_Ta");
    }
    else if (fFoil1Name == "20E") {
      foil1Thick = 0.221*mm * fFoil1ThicknessFactor;
      material = G4Material::GetMaterial("G4_Ta");
    }
    else if (fFoil1Name == "22E") {
      foil1Thick = 0.292*mm * fFoil1ThicknessFactor;
      material = G4Material::GetMaterial("G4_Ta");
    }
    else if (fFoil1Name == "Custom") {
      foil1Thick = fFoil1Thickness;
      foil1Radius = fFoil1Radius;
      material = G4Material::GetMaterial(fFoil1Material);
    }
    else {
      G4ExceptionDescription ed;
      ed << "Foil 1 undefined!";
      G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac105",
                  FatalException, ed);
    }

    fEfoil1_tubs = new G4Tubs("efoil1", 0.0*mm, foil1Radius, foil1Thick/2.,
                        0.*deg, 360.*deg);
    fEfoil1_LV = new G4LogicalVolume(fEfoil1_tubs, material,
                                     "efoil1_LV", 0, 0, 0);

    G4double foil1_pos = fSAD - foil1PositionBottom + foil1Thick/2.;
    if (fFoil1Name == "Custom") foil1_pos += fFoil1Position;

    fEfoil1 = new G4PVPlacement(0,
          G4ThreeVector(0., 0., foil1_pos - fGantryPos),
          fEfoil1_LV, "electronFoil1", fGantry_LV, false, 0);

    fEfoil1_LV->SetVisAttributes(efoil1_VisAtt);
  }
  
  // the surrounding holder
 
  if (fFoil1Name != "None" && fFoil1Name != "Custom") {
    G4VisAttributes* efoil1_VisAtt2;
    efoil1_VisAtt2 = new G4VisAttributes(G4Colour(0.4, 0.1, 0.1, 0.8));
    efoil1_VisAtt2->SetVisibility(fVis1);
    efoil1_VisAtt2->SetForceSolid(false);
    
    G4double surround_thick = 4.826*mm;
    // approximate gently curved piece as straight
    G4Box* surround_out = new G4Box("foi1_surround_out", 30.*mm, 8.5*mm,
                          surround_thick/2.);
   
    G4Tubs* surround_in = new G4Tubs("foil1_surround_in", 0., 5.5*mm,
                          surround_thick/2., 0.*deg, 360.*deg);

    G4SubtractionSolid* surround = new G4SubtractionSolid("foil1_surround", 
                        surround_out, surround_in, 0, G4ThreeVector());

    fEfoil1_holder_LV = new G4LogicalVolume(surround,
              G4Material::GetMaterial("Aluminum5052"),"surround_LV", 0, 0, 0);
    fEfoil1_holder_LV->SetVisAttributes(efoil1_VisAtt2);

    new G4PVPlacement(0,
          G4ThreeVector(0., 0., 
            fSAD - foil1PositionBottom - surround_thick/2. - fGantryPos),
          fEfoil1_holder_LV, "surround", fGantry_LV, false, 0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::BuildFoil2(G4bool build) {
  //============================================================================
  // foil 2 
  //============================================================================
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  
  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
  //G4VPhysicalVolume* myVol = store->GetVolume("electronFoil2", false);
  //if (myVol) store->DeRegister(myVol);
  //myVol = store->GetVolume("electronFoil2_custom", false);
  //if (myVol) store->DeRegister(myVol);

  std::vector<G4String> physVols;
  physVols.push_back("electronFoil2");
  physVols.push_back("electronFoil2_custom");
  physVols.push_back("foil2");
  std::vector<G4String>::iterator it;
  for (it = physVols.begin(); it != physVols.end(); ++it) {
    G4VPhysicalVolume* myVol = store->GetVolume(*it, false);
    if (myVol) store->DeRegister(myVol);
  }
 
  G4LogicalVolumeStore* lv_store = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* lv = lv_store->GetVolume("efoil2_LV", false);
  if (lv) lv_store->DeRegister(lv);
  lv = lv_store->GetVolume("efoil2_LV_custom", false);
  if (lv) lv_store->DeRegister(lv);
 
  if (!build) return; 

  if (fFoil2Name != "None") {
    G4GenericPolycone* efoil2 = nullptr;
    G4VisAttributes* efoil2_VisAtt;
    efoil2_VisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.8));
    efoil2_VisAtt->SetVisibility(fVis1);
    efoil2_VisAtt->SetForceSolid(false);
    G4double foil2PositionBottom = 135.526*mm;

    G4double foil2Position = fSAD - foil2PositionBottom;
    G4Material* foil2material = G4Material::GetMaterial("Aluminum6061");
    G4RotationMatrix* foil2Rot = new G4RotationMatrix();

    //if (fFoil2Name != "Custom") {
    if (fFoil2Name == "6E") {
      const G4int foil2pts = 14;
      G4double r_foil2[foil2pts] = {
           0.0 *mm,
           4.06*mm,  4.06*mm, 8.15*mm,  8.15*mm,
          10.17*mm, 10.17*mm,
          19.5 *mm, 19.5 *mm, 25. *mm, 25.  *mm,
          40.  *mm, 40.  *mm,
           0.  *mm};
      G4double z_foil2[foil2pts] = {
          1.524*mm,
          1.524*mm, 1.27 *mm, 1.27 *mm, 1.143*mm,
          1.143*mm, 0.152*mm,
          0.152*mm, 0.254*mm, 0.254*mm, 3.   *mm,
          3.   *mm, 0.   *mm,
          0.   *mm};
      efoil2 = new G4GenericPolycone("foil2",
              0.*deg, 360.*deg, foil2pts, r_foil2, z_foil2);
    }
    else if (fFoil2Name == "9E") {
      const G4int foil2pts = 12;
      G4double r_foil2[foil2pts] = {
           0.0 *mm,
           8.15*mm,  8.15*mm, 10.17*mm, 10.17*mm,
          19.5 *mm, 19.5 *mm, 25.  *mm, 25.  *mm,
          80.  *mm, 80.  *mm,
           0.  *mm};
      G4double z_foil2[foil2pts] = {
          1.422*mm,
          1.422*mm, 1.168*mm, 1.168*mm, 0.152*mm,
          0.152*mm, 0.254*mm, 0.254*mm, 3.   *mm,
          3.   *mm, 0.   *mm,
          0.   *mm};
      efoil2 = new G4GenericPolycone("foil2",
              0.*deg, 360.*deg, foil2pts, r_foil2, z_foil2);
    }
    else if (fFoil2Name == "12E") {
      const G4int foil2pts = 12;
      G4double r_foil2[foil2pts] = {
           0.0 *mm,
           8.15*mm,  8.15*mm, 10.17*mm, 10.17*mm,
          19.5 *mm, 19.5 *mm, 25.  *mm, 25.  *mm,
          80.  *mm, 80.  *mm,
           0.  *mm};
      G4double z_foil2[foil2pts] = {
          1.803 *mm,
          1.803 *mm, 1.168*mm, 1.168*mm, 0.1655*mm,
          0.1655*mm, 0.305*mm, 0.305*mm, 3.    *mm,
          3.    *mm, 0.   *mm,
          0.    *mm};   //.1655,.305 are min values
      efoil2 = new G4GenericPolycone("foil2",
                0.*deg, 360.*deg, foil2pts, r_foil2, z_foil2);
    }
    else if (fFoil2Name == "15E" || fFoil2Name == "16E") {
      const G4int foil2pts = 16;
      G4double r_foil2[foil2pts] = {
           0.0 *mm,
           4.06*mm,  4.06*mm,  6.11*mm,  6.11*mm,
           8.15*mm,  8.15*mm, 10.17*mm, 10.17*mm,
          19.5 *mm, 19.5 *mm, 25.  *mm, 25.  *mm,
          80.  *mm, 80.  *mm,
           0.  *mm};
      G4double z_foil2[foil2pts] = {
          2.108 *mm,
          2.108 *mm, 1.600*mm, 1.600*mm, 1.092 *mm,
          1.092 *mm, 0.660*mm, 0.660*mm, 0.1655*mm,
          0.1655*mm, 0.305*mm, 0.305*mm, 3.    *mm,
          3.    *mm, 0.   *mm,
          0.    *mm};
      efoil2 = new G4GenericPolycone("foil2",
                0.*deg, 360.*deg, foil2pts, r_foil2, z_foil2);
    }
    else if (fFoil2Name == "18E") {
      const G4int foil2pts = 16;
      G4double r_foil2[foil2pts] = {
           0.0 *mm,
           4.06*mm,  4.06*mm,  6.11*mm,  6.11*mm,
           8.15*mm,  8.15*mm, 10.17*mm, 10.17*mm,
          19.5 *mm, 19.5 *mm, 25.  *mm, 25.  *mm,
          80.  *mm, 80.  *mm,
           0.  *mm};
      G4double z_foil2[foil2pts] = {
          2.311 *mm,
          2.311 *mm, 1.803*mm, 1.803*mm, 1.168 *mm,
          1.168 *mm, 0.660*mm, 0.660*mm, 0.1655*mm,
          0.1655*mm, 0.305*mm, 0.305*mm, 3.    *mm,
          3.    *mm, 0.   *mm,
          0.    *mm};
      efoil2 = new G4GenericPolycone("foil2",
                0.*deg, 360.*deg, foil2pts, r_foil2, z_foil2);
    }
    else if (fFoil2Name == "20E" || fFoil2Name == "22E") {
      const G4int foil2pts = 14;
      G4double r_foil2[foil2pts] = {
           0.0 *mm,
           6.11*mm,  6.11*mm,
           8.15*mm,  8.15*mm, 10.17*mm, 10.17*mm,
          19.5 *mm, 19.5 *mm, 25.  *mm, 25.  *mm,
          80.  *mm, 80.  *mm,
           0.  *mm};
      G4double z_foil2[foil2pts] = {
          2.184 *mm,
          2.184 *mm, 1.168*mm,
          1.168 *mm, 0.660*mm, 0.660*mm, 0.1655*mm,
          0.1655*mm, 0.305*mm, 0.305*mm, 3.    *mm,
          3.    *mm, 0.   *mm,
          0.    *mm};
      efoil2 = new G4GenericPolycone("foil2",
                0.*deg, 360.*deg, foil2pts, r_foil2, z_foil2);
    }
    else if (fFoil2Name == "Custom") {
      foil2Rot->rotateX(pi);
      if (fFoil2Custom_z.size() > 2) {
        efoil2 =
          new G4GenericPolycone("foil2", 0., twopi, static_cast<int>(fFoil2Custom_z.size()),
                                &fFoil2Custom_r[0], &fFoil2Custom_z[0]);
      }
      else {
        G4ExceptionDescription ed;
        ed << "custom foil 2 has insufficient r/z points!";
        G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac131",
                    FatalException, ed);
      }
    }
    else {
      G4ExceptionDescription ed;
      ed << "Foil 2 undefined!";
      G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac130",
                  FatalException, ed);
    }
    fEfoil2_LV = new G4LogicalVolume(efoil2,
      foil2material, "efoil2_LV", 0, 0, 0);

    fEfoil2 = new G4PVPlacement(foil2Rot,
          G4ThreeVector(0., 0., foil2Position - fFoil2Position - fGantryPos),
          fEfoil2_LV, "electronFoil2", fGantry_LV, false, 0);

  fEfoil2_LV->SetVisAttributes(efoil2_VisAtt);
  } // is not "None"
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// custom foil2, layer 2 polycone
void TB02_TrueBeamDetectorConstruction::AddCustomFoil2Vertex(G4double r, G4double z) {
  // TODO verify inputs
  G4cout << "Adding custom foil 2 vertex r: " << r/mm << " z: " << z/mm
         << G4endl;
  fFoil2Custom_r.push_back(r);
  fFoil2Custom_z.push_back(z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::ConstructSDandField() {
  //// ============= phantom ========================
  //   in parallel world
  // ============= monitor chamber dose ========================
  G4String doseSDname = "mon_chamber_SD";
  G4MultiFunctionalDetector*& monchamber = fMonitorChamberCache.Get();

  if (!monchamber) {
    monchamber = new G4MultiFunctionalDetector(doseSDname);
    G4SDManager::GetSDMpointer()->AddNewDetector(monchamber);
    SetSensitiveDetector("ic_dose_lv", monchamber);
  }
  if (monchamber->GetNumberOfPrimitives() == 0) {
    fMonChamberDose = new G4PSDoseDeposit("monchamber");
    monchamber->RegisterPrimitive(fMonChamberDose);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetBuildBeWindow(G4bool build) {
  if (build == fBuildBeWindow) return;
  fBuildBeWindow = build;
  G4cout << "Build exit window set to: " << fBuildBeWindow << "." << G4endl;
  BuildBeWindow();
}


void TB02_TrueBeamDetectorConstruction::SetMLCNDS120Position(G4int index, G4double pos) {
  // pos == position is positive for non-overtravel, for both leaf banks
  if (index < 0 || index > 119) {
    G4ExceptionDescription ed;
    ed << "MLC index " << index << " out of range.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac060",
                FatalException, ed);
  }
  if (pos < -20.*cm || pos > 20.*cm) {
    G4ExceptionDescription ed;
    ed << "MLC position " << pos/cm << " cm out of range.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac061",
                FatalException, ed);
  }
  const G4int num_pos = 41;
  G4double mlc_nominal[num_pos] = 
    { 20.0*cm,  19.0*cm,  18.0*cm,  17.0*cm,  16.0*cm,  15.0*cm,  14.0*cm,
      13.0*cm,  12.0*cm,  11.0*cm,  10.0*cm,   9.0*cm,   8.0*cm,   7.0*cm,
       6.0*cm,   5.0*cm,   4.0*cm,   3.0*cm,   2.0*cm,   1.0*cm,   0.0*cm,
      -1.0*cm,  -2.0*cm,  -3.0*cm,  -4.0*cm,  -5.0*cm,  -6.0*cm,  -7.0*cm,
      -8.0*cm,  -9.0*cm, -10.0*cm, -11.0*cm, -12.0*cm, -13.0*cm, -14.0*cm,
     -15.0*cm, -16.0*cm, -17.0*cm, -18.0*cm, -19.0*cm, -20.0*cm};
  G4double mlc_actual[num_pos] =
    { 20.3206*cm,  19.2892*cm,  18.2593*cm,  17.2312*cm,  16.2046*cm,
      15.1797*cm,  14.1564*cm,  13.1347*cm,  12.1147*cm,  11.0963*cm,
      10.0795*cm,   9.0643*cm,   8.0508*cm,   7.0388*cm,   6.0285*cm,
       5.0198*cm,   4.0126*cm,   3.0071*cm,   2.0031*cm,   1.0008*cm,
       0.0000*cm,  -0.9992*cm,  -1.9969*cm,  -2.9930*cm,  -3.9875*cm,
      -4.9806*cm,  -5.9721*cm,  -6.9620*cm,  -7.9505*cm,  -8.9375*cm,
      -9.9230*cm, -10.9070*cm, -11.8895*cm, -12.8706*cm, -13.8503*cm,
     -14.8285*cm, -15.8053*cm, -16.7807*cm, -17.7547*cm, -18.7273*cm,
     -19.6986*cm};
  G4double actual = -999.*cm;
  G4bool found = false;
  for (G4int i = 0; i < num_pos; ++i) {
    if (pos == mlc_nominal[i]) {
      actual = mlc_actual[i];
      found = true;
      break;
    } else if (pos > mlc_nominal[i]) {
      //Millenium leaf tip radius 80mm
      actual = pos + 80*mm *(std::sqrt(1 + std::pow(pos/1000*mm,2))-1) * (1000/509.952);
      // Interpolation option:
      //  actual = (pos - mlc_nominal[i]) /
      //            (mlc_nominal[i-1] - mlc_nominal[i]) *
      //            (mlc_actual[i-1] - mlc_actual[i])
      //            + mlc_actual[i];
      //else if (fMlc_type ==  "NDS120HD") {
      //  actual = (pos - mlc_nominal[i]) /
      //           (mlc_nominal[i-1] - mlc_nominal[i]) *
      //            (hd120_actual[i-1] - hd120_actual[i])
      //            + hd120_actual[i];
      //}

      found = true;
      break;

    }
  }
  if (!found) {
    G4ExceptionDescription ed;
    ed << "Error setting MLC position.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac062",
                FatalException, ed);
  }
  //rescale to actual leaf position
  actual *= 509.952/1000.;
  const G4ThreeVector orig = NDS120[index].placement->GetObjectTranslation();

  G4GeometryManager::GetInstance()->OpenGeometry(fColl);

  G4int opp_index = index < 60 ? index + 60 : index - 60;
  // calculate the position of the leaf in G4 coordinates
  G4double newpos;
  if (index < 60) newpos = -actual; 
  else          newpos =  actual; 
  if (newpos == orig.x())  return;  // leaf already in requested position

  NDS120[index].placement->SetTranslation(G4ThreeVector(newpos, orig.y(), orig.z()));

  //check to see if there is a collision
  //G4cout << "Setting MLC position: index " << index << " position: " << (newpos)/cm << G4endl;
  if (fVerbosity > 3) {
    G4cout << "Set MLC leaf position for leaf " << index << " to " << pos/mm 
           << " mm." << " newpos: " << newpos/mm << G4endl;
  }
  // closing the geometry is slow due to reoptimization that takes place
  // don't want to do on every leaf, so reoptimize after last leaf is positioned
  if (index == 120 - 1) {
      G4GeometryManager::GetInstance()->CloseGeometry(fColl);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::
                                SetMLCHDLeafPosition(G4int index, G4double pos) {
  // pos == position is positive for non-overtravel, for both leaf banks
  if (index < 0 || index > 119) {
    G4ExceptionDescription ed;
    ed << "MLC index " << index << " out of range.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac060",
                FatalException, ed);
  }
  if (pos < -20.*cm || pos > 20.*cm) {
    G4ExceptionDescription ed;
    ed << "MLC position " << pos/cm << " cm out of range.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac061",
                FatalException, ed);
  }
  // this seems to be necessary for vis:
  //G4RunManager::GetRunManager()->GeometryHasBeenModified();
  //
  // table of nominal vs actual  leaf positions, relative to isocenter
  // here, negative is overtravel. In the file, positive is overtravel
  // from ~/montecarlo/mlc/mlc_positions.pdf
  const G4int num_pos = 41;
  // TODO these should be class variables
  G4double mlc_nominal[num_pos] = 
    { 20.0*cm,  19.0*cm,  18.0*cm,  17.0*cm,  16.0*cm,  15.0*cm,  14.0*cm,
      13.0*cm,  12.0*cm,  11.0*cm,  10.0*cm,   9.0*cm,   8.0*cm,   7.0*cm,
       6.0*cm,   5.0*cm,   4.0*cm,   3.0*cm,   2.0*cm,   1.0*cm,   0.0*cm,
      -1.0*cm,  -2.0*cm,  -3.0*cm,  -4.0*cm,  -5.0*cm,  -6.0*cm,  -7.0*cm,
      -8.0*cm,  -9.0*cm, -10.0*cm, -11.0*cm, -12.0*cm, -13.0*cm, -14.0*cm,
     -15.0*cm, -16.0*cm, -17.0*cm, -18.0*cm, -19.0*cm, -20.0*cm};
  G4double hd120_actual[num_pos] =
    { 20.6627*cm,  19.5967*cm,  18.5343*cm,  17.4754*cm,  16.4201*cm,
      15.3683*cm,  14.3200*cm,  13.2752*cm,  12.2339*cm,  11.1960*cm,
      10.1616*cm,   9.1305*cm,   8.1028*cm,   7.0785*cm,   6.0575*cm,
       5.0398*cm,   4.0254*cm,   3.0142*cm,   2.0063*cm,   1.0016*cm,
       0.0000*cm,  -0.9984*cm,  -1.9938*cm,  -2.9860*cm,  -3.9752*cm,
      -4.9614*cm,  -5.9446*cm,  -6.9249*cm,  -7.9022*cm,  -8.8766*cm,
      -9.8482*cm, -10.8170*cm, -11.7830*cm, -12.7462*cm, -13.7067*cm,
     -14.6645*cm, -15.6196*cm, -16.5721*cm, -17.5220*cm, -18.4694*cm,
     -19.4142*cm};



  // find if requested position is in the nominal array, return the actual 
  // position. Otherwise interpolate between two values

  G4double actual = -999.*cm;
  G4bool found = false;
  //G4double tmppos = -pos;  
  // signs of mlc_nominal/actual are backwards to what we want

  for (G4int i = 0; i < num_pos; ++i) {
    if (pos == mlc_nominal[i]) {
      if (fMlc_type == "NDS120HD") actual = hd120_actual[i];
      else { 
        G4ExceptionDescription ed;
        ed << "MLC model " << fMlc_type << " not found.";
        G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac131",
                    FatalException, ed);
      } 
      found = true;
      break;
    } else if (pos > mlc_nominal[i]) {
      //interpolate between i-1, i
      if (fMlc_type ==  "NDS120HD") {
        actual = (pos - mlc_nominal[i]) /
                  (mlc_nominal[i-1] - mlc_nominal[i]) *
                  (hd120_actual[i-1] - hd120_actual[i])
                  + hd120_actual[i];
      }

      found = true;
      break;

    }
  }
  // return in coordinates such that non-overtravel position is positive
  //actual *= -1.;
  if (!found) {
    G4ExceptionDescription ed;
    ed << "Error setting MLC position.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac062",
                FatalException, ed);
  }

  //rescale to actual leaf position
  actual *= 510./1000.;

  G4PVPlacement* leaf = fMlc_PV[index];

  G4GeometryManager::GetInstance()->OpenGeometry(fColl);

  G4int opp_index = index < 60 ? index + 60 : index - 60;
  const G4ThreeVector orig = leaf->GetObjectTranslation();
  G4PVPlacement* opp_leaf = fMlc_PV[opp_index];

  // calculate the position of the leaf in G4 coordinates
  G4double newpos;
  if (index < 60) newpos = -actual - fMlc_x_size/2.;
  else newpos =  actual + fMlc_x_size/2.;

  if (newpos == orig.x())  return;  // leaf already in requested position

  leaf->SetTranslation(G4ThreeVector(newpos, orig.y(), orig.z()));

  //check to see if there is a collision
  //G4cout << "Setting MLC position: index " << index << " position: " << newpos/cm << G4endl;;
  if (opp_leaf) {  // may not have been build yet
    const G4ThreeVector opp_orig = opp_leaf->GetObjectTranslation();
    //G4cout << " newpos: " << newpos/cm << " actual: " << actual/cm << " opposite pos: " << opp_orig.x()/cm << G4endl;
    //G4cout << " newpos-opppos: " << (newpos - opp_orig.x())/cm << G4endl;
    //G4cout << "check:    " << GetMLCLeafPosition(index)/cm << " cm " << G4endl;
    //G4cout << "Opp.pos.: " << GetMLCLeafPosition(opp_index)/cm << " cm." << G4endl;
    G4double opp_orig_pos = opp_orig.x() + fMlc_x_size/2.;
    if (index < 60) opp_orig_pos *= -1.;

    G4double leaf_delta = 0.01*mm;
    if (index < 60) {
      if (opp_orig.x() - newpos < fMlc_x_size) {
      //  G4cout << " newpos: " << newpos/cm << " opposite pos: " << opp_orig.x()/cm << G4endl;
      //  G4cout << " newpos-opppos: " << (newpos - opp_orig.x())/cm << G4endl;
      //  G4cout << "Setting opposite leaf " << opp_index << " to " << -pos + leaf_delta << G4endl;
        SetMLCHDLeafPosition(opp_index, -pos - leaf_delta);
      }
    }
    else if (index >= 60) {
      if (newpos - opp_orig.x() < fMlc_x_size) {
      //  G4cout << " newpos: " << newpos/cm << " opposite pos: " << opp_orig.x()/cm << G4endl;
      //  G4cout << " newpos-opppos: " << (newpos - opp_orig.x())/cm << G4endl;
      //  G4cout << "Opp.pos.: " << GetMLCLeafPosition(opp_index)/cm << " cm." << G4endl;
      //  G4cout << "Setting opposite leaf " << opp_index << " to " << -pos + leaf_delta << G4endl;
        SetMLCHDLeafPosition(opp_index, -pos - leaf_delta);
      }
    }
  }
  if (fVerbosity > 3) {
    G4cout << "Set MLC leaf position for leaf " << index << " to " << pos/mm 
           << " mm." << " newpos: " << newpos/mm << G4endl;
  }
  
  G4GeometryManager::GetInstance()->CloseGeometry(fColl);
//#ifdef G4VIS_USE
//    G4RunManager::GetRunManager()->GeometryHasBeenModified();
//#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetJawPositionY(G4int jawint,
                                                G4double field, 
                                                G4bool force) {
  G4VPhysicalVolume* jaw = nullptr;
  // without overtravel, jaw1 is -ve, jaw2 +ve.
  // variable field1 is +ve for no overtravel for both jaws.
  // because of jaw focussing back to +-1.5mm, need to keep 
  //    jaw separation >~0.5mm.
  G4double field1 = field * (2.*(G4double)jawint - 3.);
  if ((jawint == 1 && field > fFieldY2 - 0.5*mm) || 
      (jawint == 2 && field < fFieldY1 + 0.5*mm)) {
    G4cout << "Jaw collision while setting jawY" << jawint << " to " 
           << field/cm << " cm." << G4endl
           << "Moving other jaw too." << G4endl;
    if      (jawint == 1) SetJawPositionY(2, field + 0.5*mm);
    else if (jawint == 2) SetJawPositionY(1, field - 0.5*mm);

  }

  if (jawint == 1) {
    if (field == fFieldY1 && !force) return;
    fFieldY1 = field;
    if (fSimulateCollimators) jaw = fJawY1;
  } else if (jawint == 2) {
    if (field == fFieldY2 && !force) return;
    fFieldY2 = field;
    if (fSimulateCollimators) jaw = fJawY2;
  }

  if (jaw && fSimulateCollimators) {
    G4GeometryManager::GetInstance()->OpenGeometry(fColl);
    // (center of face of) jaws moves in an arc with radius of curvature
    const G4double radius = 318.377*mm;
    const G4double ysize  = 59.69*mm; //half y-size of jaw
    // coords of center of jaw
    G4double ycenter;
    G4double zcenter;
    // coords of centre of face of jaw
    G4double yedge;
    G4double zedge;
    G4double theta;
    G4double theta_0;  // angle, with no fJawOffset
    //  fieldY1 is negative for open field
    //if (abs(field) <= fJawOffset) { //protection for closed fields
    //  theta = 0.;
    //  ycenter = ysize + fJawOffset;
    //  zcenter = fSAD - radius;
    //} else {
    
    theta   = std::atan((field1-fJawOffset)/fSAD);
    theta_0 = std::atan((field1)/fSAD);
    //yedge = radius * std::sin(theta); // + fJawOffset;
    yedge = radius * std::sin(theta_0); // + fJawOffset;
    zedge = fSAD - radius * std::cos(theta_0);
    ycenter = yedge + ysize * std::cos(theta_0);
    zcenter = zedge + ysize * std::sin(theta_0);
    
    //}

    G4RotationMatrix* rotate = new G4RotationMatrix();
    // +90 for jawint=1; -90 for jawint=2
    rotate->rotateZ(90.*deg - (G4double)(jawint-1) * 180.*deg);
    rotate->rotateY(theta);

    jaw->SetRotation(rotate);
    if (jawint == 1) ycenter *= -1.;
    jaw->SetTranslation(G4ThreeVector(0., ycenter, zcenter - fCollPos));

    G4GeometryManager::GetInstance()->CloseGeometry(fColl);
//#ifdef G4VIS_USE
//    G4RunManager::GetRunManager()->GeometryHasBeenModified();
//#endif

    G4cout << "Jaw Y" << jawint << " set to : " << field/cm << " cm." 
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetJawPositionX(G4int jawint,
                                                G4double field,
                                                G4bool force) {
  // allow overtravel. without overtravel, jaw1 is -ve, jaw2 +ve
  // convert to positive
  G4double field1 = field * (G4double)(2*jawint - 3);
  G4VPhysicalVolume* jaw = nullptr;
  if ((jawint == 1 && field > fFieldX2 - 0.5*mm) || 
      (jawint == 2 && field < fFieldX1 + 0.5*mm)){
    G4cout << "Jaw collision while setting jawX" << jawint << " to " 
           << field/cm << " cm." << G4endl
           << "Moving other jaw too." << G4endl;
    if      (jawint == 1) SetJawPositionX(2, field + 0.5*mm);
    else if (jawint == 2) SetJawPositionX(1, field - 0.5*mm);
  }

  if (jawint == 1) {
    if (field == fFieldX1 && !force) return;
    fFieldX1 = field;
    if (fSimulateCollimators) jaw = fJawX1;
  } else if (jawint == 2) {
    if (field == fFieldX2 && !force) return;
    fFieldX2 = field;
    if (fSimulateCollimators) jaw = fJawX2;
  }

  if (jaw && fSimulateCollimators) { 
    G4GeometryManager::GetInstance()->OpenGeometry(fColl);
    // jaw rotates about trunion. trunion z is constant
    // the z position of center of jaw, when closed (from drawing 10034445-2)
    // differs from MC package by 0.5 mm
    const G4double jaw_center_z = fSAD-405.415*mm;  
    const G4double xsize        = 63.43*mm;  // half-size of jaw 
    //const G4double zsize        = 77.724*mm/2;//not used
    G4double theta; //angle of jaw jace
    G4double theta_0; // position, without fJawOffset
    // position of trunion, relative to center (note sign; trunion is closer
    // to beam axis, and target, then center
    const G4double xtrun_rel    = 33.21*mm; 
    const G4double ztrun_rel    = 4.826*mm; 
    //coords of center of jaw
    G4double       xcenter;
    G4double       zcenter;
    //coords of trunion
    G4double       xtrunion;
    const G4double ztrunion = jaw_center_z + ztrun_rel;
    G4double x1;  //x coord of intersection of hypotenuse, line of constant z
                  // through trunion

    //if (abs(field) <= fJawOffset) {//protection for div by 0 for closed fields
    //  theta = 0.;
    //  xcenter = xsize + fJawOffset;
    //  zcenter = jaw_center_z;
    //} else {
    
    theta   = std::atan((field1 - fJawOffset)/fSAD);
    theta_0 = std::atan((field1)/fSAD);
    //x1 = (field1 - fJawOffset) * (fSAD - ztrunion)/fSAD + fJawOffset; 
    x1 = (field1) * (fSAD - ztrunion)/fSAD; // + fJawOffset; 
    xtrunion = x1 + (xsize - xtrun_rel)/std::cos(theta_0);
    xcenter = xtrunion + ( std::cos(theta_0)*xtrun_rel
                       + std::sin(theta_0)*ztrun_rel);
    zcenter = ztrunion - (-std::sin(theta_0)*xtrun_rel
                       + std::cos(theta_0)*ztrun_rel);
    
    //} 

    G4RotationMatrix* rotate = new G4RotationMatrix();
    rotate->rotateY(180.*deg - (G4double)(jawint-1)*180.*deg);
    rotate->rotateX(-90.*deg + (G4double)(jawint-1)*180.*deg);
    rotate->rotateZ(theta);
    jaw->SetRotation(rotate);
    if (jawint == 1) xcenter *= -1.;
    jaw->SetTranslation(G4ThreeVector(xcenter, 0., zcenter - fCollPos));
    
    G4GeometryManager::GetInstance()->CloseGeometry(fColl);
//#ifdef G4VIS_USE
//    G4RunManager::GetRunManager()->GeometryHasBeenModified();
//#endif

    G4cout << "Jaw X" << jawint << " set to : " << field/cm << " cm." 
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetSimulateCollimators(G4bool v) {
  if (v == fSimulateCollimators) return;
  fSimulateCollimators = v;
  BuildCollimators(fSimulateCollimators);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4cout << "Simulate collimators: " << fSimulateCollimators << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetSimulateShielding(G4bool v) {
  if (v == fSimulateShielding) return;
  fSimulateShielding = v;
  // put the logic whether to build/destroy into BuildShielding
  BuildShielding();
  BuildBackscatterShield(!fSimulateShielding);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4cout << "Simulate shielding: " << fSimulateShielding << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetJawOffset(G4double v) {
  if (v == fJawOffset) return;

  fJawOffset = v;
  SetJawPositionX(1, fFieldX1, true);
  SetJawPositionX(2, fFieldX2, true);
  SetJawPositionY(1, fFieldY1, true);
  SetJawPositionY(2, fFieldY2, true);
  
  //(void)v;
  //G4ExceptionDescription ed;
  //ed << "Jaw offset command not enabled.";
  //G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac111",
  //            FatalException, ed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetUseIonChamber(G4bool b) {
  (void)b;
  //TODO make this an option?
  G4ExceptionDescription ed;
  ed << "Use ion chamber command not enabled.";
  G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac110",
              FatalException, ed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetExitWindowThicknessFactor(
                                                                  G4double v) {
  (void)v;
  G4ExceptionDescription ed;
  ed << "Exit window thickness factor not enabled.";
  G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac109",
              FatalException, ed);
  //fExitWindowThicknessFactor=v;
  //  //TODO don't reinitialize
  //G4RunManager::GetRunManager()->ReinitializeGeometry();
  //G4cout << "Exit window thickness factor: "<<v<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetFoil2ThicknessFactor(G4double v) {
  G4cout << "Asked for: " << v << G4endl; // to avoid compiler warning

  G4ExceptionDescription ed;
  ed << "Foil2 thickness factor not in use.";
  G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac039",
              FatalException, ed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetCopperThicknessFactor(G4double v) {
  if (v == fCopperThicknessFactor) return;
  G4double originalFactor = fCopperThicknessFactor;
  fCopperThicknessFactor = v;
  G4cout << "Setting ion chamber metal thickness factor to: "
         << fCopperThicknessFactor << G4endl;

  // Window
  G4double thick = fICWinCu_tubs->GetZHalfLength();
  fICWinCu_tubs->SetZHalfLength(thick*fCopperThicknessFactor);

  G4ThreeVector pos = fICWinCu0->GetTranslation();
  G4double new_z = pos.z() + fICCuHThick *
                            (fCopperThicknessFactor - originalFactor);
  pos.setZ(new_z);
  fICWinCu0->SetTranslation(pos);

  pos = fICWinCu1->GetTranslation();
  new_z = pos.z() + fICCuHThick *
                            (fCopperThicknessFactor - originalFactor);
  pos.setZ(new_z);
  fICWinCu1->SetTranslation(pos);

  pos = fICWinCu2->GetTranslation();
  new_z = pos.z() + fICCuHThick *
                            (fCopperThicknessFactor - originalFactor);
  pos.setZ(new_z);
  fICWinCu2->SetTranslation(pos);

  // hv electrode
  thick = fICElCu_tubs->GetZHalfLength();
  fICElCu_tubs->SetZHalfLength(thick*fCopperThicknessFactor);

  pos = fICElCu0->GetTranslation();
  new_z = pos.z() + fICCuHThick *
                            (fCopperThicknessFactor - originalFactor);
  pos.setZ(new_z);
  fICElCu0->SetTranslation(pos);

  pos = fICElCu1->GetTranslation();
  new_z = pos.z() + fICCuHThick *
                            (fCopperThicknessFactor - originalFactor);
  pos.setZ(new_z);
  fICElCu1->SetTranslation(pos);

  // signal electrode
  thick = fICElCusig_tubs->GetZHalfLength();
  fICElCusig_tubs->SetZHalfLength(thick*fCopperThicknessFactor);

  pos = fICElCusig0->GetTranslation();
  new_z = pos.z() + fICCuHThick *
                            (fCopperThicknessFactor - originalFactor);
  pos.setZ(new_z);
  fICElCusig0->SetTranslation(pos);

  pos = fICElCusig1->GetTranslation();
  new_z = pos.z() + fICCuHThick *
                            (fCopperThicknessFactor - originalFactor);
  pos.setZ(new_z);
  fICElCusig1->SetTranslation(pos);

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//void TB02_TrueBeamDetectorConstruction::SetKaptonThicknessFactor(G4double v) {
//  G4cerr << "Kapton thickness factor deprecated." << G4endl;
//  exit(EXIT_FAILURE);
//  //if (v == fKaptonThicknessFactor) return;
//  //fKaptonThicknessFactor = v;
//  //BuildIonChamber();
//  //G4cout << "I.C. kapton thickness factor: " << v << G4endl;
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetFoil1Name(G4String v) {
  if (v == fFoil1Name) return;
  if      (v == "6E" ) fFoil1Name = "6E" ;
  else if (v == "9E" ) fFoil1Name = "9E" ;
  else if (v == "12E") fFoil1Name = "12E";
  else if (v == "15E") fFoil1Name = "15E";
  else if (v == "16E") fFoil1Name = "16E";
  else if (v == "18E") fFoil1Name = "18E";
  else if (v == "20E") fFoil1Name = "20E";
  else if (v == "22E") fFoil1Name = "22E";
  else if (v == "Custom") fFoil1Name = "Custom";
  else if (v == "None")   fFoil1Name = "None";
  else {
    G4ExceptionDescription ed;
    ed << "Foil 1 definition " << fFoil1Name << " not found.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac108",
                FatalException, ed);
  }
  if (fBeamType == "electron") BuildFoil1();
  G4cout << "Foil1: " << fFoil1Name<<G4endl;;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetFoil1ThicknessFactor(G4double v) {
  // TODO limits
  if (v == fFoil1ThicknessFactor) return;
  fFoil1ThicknessFactor = v;
  if (fBeamType == "electron") BuildFoil1();
  G4cout << "Foil1 thickness factor: " << v << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetFoil1Material(G4String v) {
  if (v == fFoil1Material) return;
  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(v);
  
  if (pttoMaterial) {
    fFoil1Material = v;
    if (fEfoil1_LV) fEfoil1_LV->SetMaterial(pttoMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4ExceptionDescription ed;
    ed << "Foil 1 custom material not allowed.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac032",
                FatalException, ed);
  }
  if (fVerbosity>0) {
    G4cout << "Custom foil 1 material set to: " << fFoil1Material << G4endl;
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetFoil1Thickness(G4double v) {
  if (v == fFoil1Thickness) return;
  if (v > 0. && v < 2.*cm) {
    fFoil1Thickness = v;
    if (fBeamType == "electron" && fFoil1Name == "Custom") BuildFoil1();
    G4cout << "Foil 1 custom thickness [mm]: " << fFoil1Thickness/mm << G4endl;
  } else {
    G4ExceptionDescription ed;
    ed << "Foil 1 custom thickness " << fFoil1Thickness/mm 
       << " mm out of range.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac106",
                FatalException, ed);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetFoil1Radius(G4double v) {
  if (v == fFoil1Radius) return;
  if (v > 0. && v < 5.*cm) {
    fFoil1Radius = v;
    if (fBeamType == "electron" && fFoil1Name == "Custom") BuildFoil1();
    G4cout << "Foil 1 custom radius [mm]: " << fFoil1Radius/mm << G4endl;
  } else {
    G4ExceptionDescription ed;
    ed << "Foil 1 custom radius " << fFoil1Radius/mm 
       << " mm out of range.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac107",
                FatalException, ed);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetFoil1Position(G4double v) {
  if (v == fFoil1Position) return;
  if (v > -5.*cm && v < 5.*cm) {
    fFoil1Position = v;
    if (fBeamType == "electron" && fFoil1Name == "Custom") BuildFoil1();
    G4cout << "Foil 1 custom position [mm]: " << fFoil1Position/mm << G4endl;
  } else {
    G4ExceptionDescription ed;
    ed << "Foil 1 custom position " << fFoil1Position/mm 
       << " mm out of range.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac108",
                FatalException, ed);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetFoil2Name(G4String v) {
  if (v == fFoil2Name) return;
  if      (v == "6E" ) fFoil2Name = "6E" ;
  else if (v == "9E" ) fFoil2Name = "9E" ;
  else if (v == "12E") fFoil2Name = "12E";
  else if (v == "15E") fFoil2Name = "15E";
  else if (v == "16E") fFoil2Name = "16E";
  else if (v == "18E") fFoil2Name = "18E";
  else if (v == "20E") fFoil2Name = "20E";
  else if (v == "22E") fFoil2Name = "22E";
  else if (v == "Custom") fFoil2Name = "Custom";
  else if (v == "None")   fFoil2Name = "None";
  else {
    G4ExceptionDescription ed;
    ed << "Foil 2 definition " << fFoil2Name << " not found.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac112",
                FatalException, ed);
  }
  if (fBeamType == "electron") BuildFoil2();
  G4cout << "Foil2: " << fFoil2Name << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetFoil2Material(G4String v) {
  if (v == fFoil1Material) return;
  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(v);
  
  if (pttoMaterial) {
    fFoil2Material = v;
    if (fEfoil2_LV) fEfoil2_LV->SetMaterial(pttoMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4ExceptionDescription ed;
    ed << "Foil 2 custom material not allowed.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac032",
                FatalException, ed);
  }

  if (fVerbosity>0) {
    G4cout << "Custom foil 2 material set to: " << fFoil2Material << G4endl;
  }
}

//TODO for all the SetFoilX, would be better to set the value of the existing
//LV and PV rather than rebuilding
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetFoil2Position(G4double v) {
  if (v == fFoil2Position) return;
  G4double old = fFoil2Position;
  if (v >- 3.*mm && v < 5.*cm) {
    fFoil2Position = v;
    if (fBeamType == "electron" && fFoil2Name == "Custom") {
      G4ThreeVector oldpos = fEfoil2->GetTranslation();
      fEfoil2->SetTranslation(
        G4ThreeVector(oldpos.x(), oldpos.y(),
                      oldpos.z() - old + fFoil2Position));
      G4RunManager::GetRunManager()->GeometryHasBeenModified();
      //BuildFoil2();
      G4cout << "Foil 2 custom position [mm]: " << fFoil2Position/mm << G4endl;
    }
  } else {
    G4ExceptionDescription ed;
    ed << "Foil 2 custom position " << v/mm << " mm out of range.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac033",
                FatalException, ed);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetTargetName(G4String v) {
  if (v == fTargetName) return;
  if (v == "LowEnergy"    ||
      v == "MediumEnergy" ||
      v == "HighEnergy"   ||
      v == "Imaging"      ||
      v == "Custom") {
    fTargetName = v;
    if (fBeamType == "xray") BuildTarget();
    if (fVerbosity > 0) {
      G4cout << "Target set to: " << fTargetName << G4endl;
    }
  } else {
    G4ExceptionDescription ed;
    ed << "Target " << v << " undefined.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac034",
                FatalException, ed);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetTargetPosition(G4double v) {
  if (v == fTargetPosition) return;
  //TODO limits should depend on thickness!
  //NB. this may collide with primary collimator. It is up to user to avoid.
  if (v >= -150.*mm && v <= 2.*mm) {
    fTargetPosition = v;
    if (fTargetName == "Custom" && fBeamType == "xray") {
      fTargetCustom->SetTranslation(G4ThreeVector(0., 0.,
          fSAD + fTargetPosition - fTargetThickness/2. - fGantryPos));
      G4RunManager::GetRunManager()->GeometryHasBeenModified();
    }
    if (fVerbosity > 0) {
      G4cout << "Custom target position set to: " << fTargetPosition/mm 
             << " mm. " << G4endl;
    }
  } else {
    G4ExceptionDescription ed;
    ed << "Custom target position " << v/mm << " mm out of range.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac035",
                FatalException, ed);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetTargetRadius(G4double v) {
  if (v == fTargetRadius) return;
  if (v >= 0. && v < 10.*cm) {
    fTargetRadius = v;
    if (fVerbosity > 0) {
      G4cout << "Custom target radius set to: " << fTargetRadius/mm 
             << " mm." << G4endl;
    }
    if (fTargetName == "Custom") {
      fTargetCustomTubs->SetOuterRadius(fTargetRadius);
      G4RunManager::GetRunManager()->GeometryHasBeenModified();
    }
  } else {
    G4ExceptionDescription ed;
    ed << "Custom target radius " << v/mm << " mm out of range.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac036",
                FatalException, ed);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetTargetThickness(G4double v) {
  if (v == fTargetThickness) return;
  // TODO large thicknesses may cause overlap with exit window
  if (v >= 0. && v < 5.*cm) {
    fTargetThickness = v;
    if (fVerbosity>0) {
      G4cout << "Custom target thickness set to: " << fTargetThickness/mm
             << " mm." << G4endl;
    }
    if (fTargetName == "Custom") {
      fTargetCustomTubs->SetZHalfLength(fTargetThickness/2.);
      fTargetCustom->SetTranslation(G4ThreeVector(0., 0.,
        fSAD + fTargetPosition - fTargetThickness/2. - fGantryPos));
      if (fTargetCustom2) {  // layer 2 abuts layer 1
        fTargetCustom2->SetTranslation(G4ThreeVector(0., 0.,
          fSAD + fTargetPosition - fTargetThickness - fGantryPos));
      }
      G4RunManager::GetRunManager()->GeometryHasBeenModified();
    }
  } else {
    G4ExceptionDescription ed;
    ed << "Custom target thickness " << v/mm << " mm out of range.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac037",
                FatalException, ed);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetTargetMaterial(G4String v) {
  if (v == fTargetMaterial) return;

  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(v);
  
  if (pttoMaterial) {
    fTargetMaterial = v;
    if (fTargetCustom_LV) fTargetCustom_LV->SetMaterial(pttoMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4ExceptionDescription ed;
    ed << "Target 1 custom material not allowed.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac031",
                FatalException, ed);
  }
  if (fVerbosity>0) {
    G4cout << "Custom target material set to: " << fTargetMaterial << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetTargetMaterial2(G4String v) {
  if (v == fTargetMaterial2) return;
  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(v);
  
  if (pttoMaterial) {
    fTargetMaterial2 = v;
    if (fTargetCustom2_LV) fTargetCustom2_LV->SetMaterial(pttoMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4ExceptionDescription ed;
    ed << "Target 2 custom material not allowed.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac030",
                FatalException, ed);
  }
  if (fVerbosity>0) {
    G4cout << "Custom target material 2 set to: " << fTargetMaterial2 << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// custom target, layer 2 polycone
void TB02_TrueBeamDetectorConstruction::AddCustomTargetVertex(
                                                      G4double r, G4double z) {
  // TODO verify inputs
  G4cout << "Adding custom target vertex r: " << r/mm << " z: " << z/mm
         << G4endl;
  fTargetCustom2_r.push_back(r);
  fTargetCustom2_z.push_back(z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetFFName(G4String v) {
  if (v == fFlatteningFilterName) return;
  fFlatteningFilterName = v;
  if (fBeamType == "xray") BuildFlatteningFilter();
  if (fVerbosity>0) {
    G4cout << "Flattening filter set to: " << fFlatteningFilterName << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetFFOffset(G4ThreeVector v) {
  if (v == fFFOffset) return;
  if (v.x() >= -1.*cm && v.x() <= 1.*cm &&
      v.y() >= -1.*cm && v.y() <= 1.*cm &&
      v.z() >= -1.*cm && v.z() <= 1.*cm) {
    fFFOffset = v;
    BuildFlatteningFilter();
    if (fVerbosity>0) {
      G4cout << "Flattening filter offset set to: " 
             << fFFOffset.x()/mm << " mm [X], "
             << fFFOffset.y()/mm << " mm [Y], "
             << fFFOffset.z()/mm << " mm [Z]. " << G4endl;
    }
  } else {
    G4ExceptionDescription ed;
    ed << "Flattening filter offset out of allowed range.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac038",
                FatalException, ed);
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double TB02_TrueBeamDetectorConstruction::GetJawPositionX(G4int jaw) {
  if (jaw==1) return fFieldX1;
  if (jaw==2) return fFieldX2;
  else return -9999.*cm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double TB02_TrueBeamDetectorConstruction::GetJawPositionY(G4int jaw) {
  if (jaw==1) return fFieldY1;
  if (jaw==2) return fFieldY2;
  else return -9999.*cm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetBeamType(G4String v) {
  if (v == fBeamType) return;
  fBeamType = v;
  if (v == "electron") {
    BuildPhotonGeometry(false);
    BuildElectronGeometry(true);
  }
  else if (v == "xray") {
    BuildElectronGeometry(false);
    BuildPhotonGeometry(true);
  }
  else {
    G4ExceptionDescription ed;
    ed << "Beam type " << v << " undefined";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac121",
                FatalException, ed);
  }
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetApplicatorName(G4String v) {
  if (v == fApplicatorName) return;
  if (v == "6x6" || v == "10x10" || v == "15x15" || v == "20x20" || 
      v == "25x25" || v == "None") {
    fApplicatorName = v;
    BuildApplicator();
    BuildCutOut();
    //TODO don't reinitialize
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    if (fVerbosity > 0) {
      G4cout << "Applicator set to: " << fApplicatorName << G4endl;
    }
  } else {
    G4ExceptionDescription ed;
    ed << "Applicator " << v << " undefined";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac120",
                FatalException, ed);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetBuildCutOut(G4bool v) {
  // call this once, after all the vertices have been defined
  if (v == fBuildCutOut) return;
  fBuildCutOut = v;
  BuildCutOut();
  if (fVerbosity > 0) {
    G4cout << "Build cut out set to: " << fBuildCutOut << G4endl;
  }
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::AddCutOutVertex(
                                                      G4double x, G4double y) {
  // expected to specify which applicator first
  G4double maxsize = 0.;
  if (fApplicatorName == "6x6")   maxsize =  28.25*mm;
  if (fApplicatorName == "10x10") maxsize =  47.0 *mm;
  if (fApplicatorName == "15x15") maxsize =  70.4 *mm;
  if (fApplicatorName == "20x20") maxsize =  93.9 *mm;
  if (fApplicatorName == "25x25") maxsize = 117.3 *mm;

  if (x>maxsize || y>maxsize) {
    G4ExceptionDescription ed;
    ed << "Cutout size too large for applicator.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac040",
                FatalException, ed);
  }

  fCutOutVertices.push_back(G4TwoVector(x,y));
  G4cout << "Added cutout vertex point at x: " << x/cm << " y: " << y/cm 
         << " cm." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetCutOutThickness(G4double v) {
  if (v == fCutOutThickness) return;
  if (v > 0. && v <= 30.*mm) {
    fCutOutThickness = v;
    // it does not look as if G4ExtrudedSolids can be modified after building!
    // thus need to call BuildCutOut again
    BuildCutOut();
  } else {
    G4ExceptionDescription ed;
    ed << "Cutout thickness out of range 0 < t <= 30*mm.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac051",
                FatalException, ed);
  }
  if (fVerbosity > 0) {
    G4cout << "Cutout thickness set to " << GetCutOutThickness()/mm 
           << " mm." << G4endl;
  }
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetCutOutBevelFactor(G4double v) {
  if (v == fCutOutBevelFactor) return;
  if (v > 0.5 && v <= 2.0) {
    fCutOutBevelFactor = v;
    BuildCutOut();
  } else {
    G4ExceptionDescription ed;
    ed << "Cutout bevel factor out of range 0.5 < f <= 2.0.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac052",
                FatalException, ed);
  }
  if (fVerbosity > 0) {
    G4cout << "Cutout bevel factor set to " << GetCutOutBevelFactor() 
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetCutOutMaterial(G4String v) {
  if (v == fCutOutMaterial) return;

  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(v);
 
  if (pttoMaterial) {
    fCutOutMaterial = v;
    BuildCutOut();
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4ExceptionDescription ed;
    ed << "Material " << v << " requested for cutout not allowed.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac053",
                FatalException, ed);
  }
  if (fVerbosity>0) {
    G4cout << "Cutout material set to: " << fCutOutMaterial << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetMLC(G4String type) {
  if (type == "NDS120CAD") {
    //RemoveMLC();
    BuildCADMLC();
  } else if (type == "NDS120HD") {
    //RemoveMLC();
    BuildHDMLC();
  } else {
    G4ExceptionDescription ed;
    ed << "Setting MLC to value not allowed (allowed are NDS120CAD and NDS120HD).";
    G4Exception("VirtuaLinac::SetMLC",
                "VirtuaLinac009", FatalException, ed);
  }
  fMlc_type = type;
  G4cout << "MLC model set to: " << fMlc_type << G4endl; 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::SetMLCDensity(G4double v) {
  std::ostringstream strs;
  strs << v;
  G4String mat = G4String("W95_") + strs.str();
  G4cout << "Will create material: " << mat << G4endl;

  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(mat);
 
  if (!pttoMaterial) {
    AddMaterial(mat, v*(g/cm3), "W95");
    pttoMaterial = 
      G4NistManager::Instance()->FindOrBuildMaterial(mat);
  }

  if (pttoMaterial) {
    fMlc_material = mat;
    RemoveMLC();
    if      (fMlc_type == "NDS120CAD")   BuildCADMLC(); 
    else if (fMlc_type == "NDS120HD") BuildHDMLC();
    else {
      G4ExceptionDescription ed;
      ed << "Error building MLC while setting density.";
      G4Exception("VirtuaLinac::SetMLC",
                  "VirtuaLinac309", FatalException, ed);
    }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4ExceptionDescription ed;
    ed << "Mlc custom material not allowed.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac131",
                FatalException, ed);
  }
  if (fVerbosity>0) {
    G4cout << "Custom MLC material set to: " << fMlc_material << G4endl;
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::RemoveMLC() {

  G4cout << "Removing MLC" << G4endl;
  G4PhysicalVolumeStore* pv_store = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4String> physVols;
  //mlc120
  physVols.push_back("NDS_leaf");
  //hd120
  physVols.push_back("hd120_outboard1_A");
  physVols.push_back("hd120_outboard60_A");
  physVols.push_back("hd120_half_target_A");
  physVols.push_back("hd120_half_iso_A");
  physVols.push_back("hd120_quarter_target_A");
  physVols.push_back("hd120_quarter_iso_A");
  physVols.push_back("hd120_outboard1_B");
  physVols.push_back("hd120_outboard60_B");
  physVols.push_back("hd120_half_target_B");
  physVols.push_back("hd120_half_iso_B");
  physVols.push_back("hd120_quarter_target_B");
  physVols.push_back("hd120_quarter_iso_B");

  std::vector<G4String>::iterator it;
  for (it = physVols.begin(); it != physVols.end(); ++it) {
    G4VPhysicalVolume* myVol = pv_store->GetVolume(*it, false);
    while (myVol) {
      pv_store->DeRegister(myVol);
      myVol = pv_store->GetVolume(*it, false);
    }
  }

  G4LogicalVolumeStore* lv_store = G4LogicalVolumeStore::GetInstance();
  std::vector<G4String> logVols;
  logVols.push_back("leaf_lv");
  //hd120
  logVols.push_back("mlc_hd120_half_target_A");
  logVols.push_back("mlc_hd120_half_target_B");
  logVols.push_back("mlc_hd120_half_iso_A");
  logVols.push_back("mlc_hd120_half_iso_B");
  logVols.push_back("mlc_hd120_quarter_target_A");
  logVols.push_back("mlc_hd120_quarter_target_B");
  logVols.push_back("mlc_hd120_quarter_iso_A");
  logVols.push_back("mlc_hd120_quarter_iso_B");
  logVols.push_back("mlc_hd120_outboard1_A");
  logVols.push_back("mlc_hd120_outboard1_B");
  logVols.push_back("mlc_hd120_outboard60_A");
  logVols.push_back("mlc_hd120_outboard60_B");
  

  for (it = logVols.begin(); it != logVols.end(); ++it) {
    G4LogicalVolume* myVol = lv_store->GetVolume(*it, false);
    while (myVol) {
      lv_store->DeRegister(myVol);
      myVol = lv_store->GetVolume(*it, false);
    }
  }

  for (G4int i = 0; i < 120; ++i) {
    fMlc_PV[i] = nullptr;
    //delete fMlc_PV[i];
  }

  if (fMlc_PV) delete fMlc_PV;
 
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void TB02_TrueBeamDetectorConstruction::BuildCADMLC() {
  // Leaves are imported from CAD-model as STL-files. By using CADMesh able to get solids of the leaves so that the centerpoint is in the target (0,0,0)
  // +-9mm shift
  // This means that there is no need for setting the positions in y,z directions but only in leaf movement direction
  // BankA, BankB

  //TETGEN: TetrahedralMesh (load meshes into Geant4, and fill them with tetrahedra. This requires adding `tetgen` to your project as a dependency.)
  //auto leaf_mesh = CADMesh::TetrahedralMesh::FromSTL(path_string);
  //leaf_mesh->SetMaterial(G4Material::GetMaterial(fMlc_material_NDS));
  //auto leaf_assembly = leaf_mesh->GetAssembly();
  //auto position = G4ThreeVector(0., 0., fSAD - fCollPos);
  //auto rotation = new G4RotationMatrix();
  //leaf_assembly->MakeImprint(fColl_LV, position, rotation);
  
  std::string directory = m_stl_folder + "/MLC120-LEAF-STL";
  std::vector<G4String> file_names;
  int dir_length;
  
  G4VisAttributes *visAtt;
  std::string path_string;
  
  // clear old geometries if they exist
  std::vector<G4String> volume_names;
  volume_names.push_back("NDS_leaf");
  volume_names.push_back("leaf_lv");
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  const G4int num_leaves = 120;
  std::string mlc_names[num_leaves] = {"100026959-02-ao60-60.stl", "100026959-04-bo60-60.stl", "100026960-01-ao1-01.stl", "100026960-02-bo1-01.stl",
  "1105333-15-at1-11.stl", "1105333-15-at1-19.stl", "1105333-15-at1-27.stl",
  "1105333-15-at1-35.stl", "1105333-15-at1-43.stl" ,"1105333-16-at2-15.stl",
  "1105333-16-at2-23.stl", "1105333-16-at2-31.stl", "1105333-16-at2-39.stl",
  "1105333-16-at2-47.stl", "1105333-17-at3-13.stl", "1105333-17-at3-21.stl",
  "1105333-17-at3-29.stl", "1105333-17-at3-37.stl", "1105333-17-at3-45.stl",
  "1105333-18-at4-17.stl", "1105333-18-at4-25.stl", "1105333-18-at4-33.stl",
  "1105333-18-at4-41.stl", "1105333-18-at4-49.stl", "1105333-19-bt1-11.stl",
  "1105333-19-bt1-19.stl", "1105333-19-bt1-27.stl", "1105333-19-bt1-35.stl",
  "1105333-19-bt1-43.stl", "1105333-20-bt2-15.stl", "1105333-20-bt2-23.stl",
  "1105333-20-bt2-31.stl", "1105333-20-bt2-39.stl", "1105333-20-bt2-47.stl",
  "1105333-21-bt3-13.stl", "1105333-21-bt3-21.stl", "1105333-21-bt3-29.stl",
  "1105333-21-bt3-37.stl", "1105333-21-bt3-45.stl", "1105333-22-bt4-17.stl",
  "1105333-22-bt4-25.stl", "1105333-22-bt4-33.stl", "1105333-22-bt4-41.stl",
  "1105333-22-bt4-49.stl", "1105334-14-ai1-12.stl", "1105334-14-ai1-16.stl",
  "1105334-14-ai1-20.stl", "1105334-14-ai1-24.stl", "1105334-14-ai1-28.stl",
  "1105334-14-ai1-32.stl", "1105334-14-ai1-36.stl", "1105334-14-ai1-40.stl",
  "1105334-14-ai1-44.stl", "1105334-14-ai1-48.stl", "1105334-15-ai2-26.stl",
  "1105334-15-ai2-34.stl", "1105334-15-ai2-42.stl", "1105334-15-ai2-50.stl",
  "1105334-15-at2-18.stl", "1105334-16-a13-14.stl", "1105334-16-ai3-22.stl",
  "1105334-16-ai3-30.stl", "1105334-16-ai3-38.stl", "1105334-16-ai3-46.stl",
  "1105334-17-bi1-12.stl", "1105334-17-bi1-16.stl", "1105334-17-bi1-20.stl",
  "1105334-17-bi1-24.stl", "1105334-17-bi1-28.stl", "1105334-17-bi1-32.stl",
  "1105334-17-bi1-36.stl", "1105334-17-bi1-40.stl", "1105334-17-bi1-44.stl",
  "1105334-17-bi1-48.stl", "1105334-18-bi2-18.stl", "1105334-18-bi2-26.stl",
  "1105334-18-bi2-34.stl", "1105334-18-bi2-42.stl", "1105334-18-bi2-50.stl",
  "1105334-19-bi1-14.stl", "1105334-19-bi3-22.stl", "1105334-19-bi3-30.stl",
  "1105334-19-bi3-38.stl", "1105334-19-bi3-46.stl", "1105335-14-af1-03.stl",
  "1105335-14-af1-06.stl", "1105335-14-af1-09.stl", "1105335-14-af1-52.stl",
  "1105335-14-af1-55.stl", "1105335-14-af1-58.stl", "1105335-15-af2-04.stl",
  "1105335-15-af2-07.stl", "1105335-15-af2-10.stl", "1105335-15-af2-53.stl",
  "1105335-15-af2-56.stl", "1105335-15-af2-59.stl", "1105335-16-af3-02.stl",
  "1105335-16-af3-05.stl", "1105335-16-af3-08.stl", "1105335-16-af3-51.stl",
  "1105335-16-af3-54.stl", "1105335-16-af3-57.stl", "1105335-17-bf1-03.stl",
  "1105335-17-bf1-06.stl", "1105335-17-bf1-09.stl", "1105335-17-bf1-52.stl",
  "1105335-17-bf1-55.stl", "1105335-17-bf1-58.stl", "1105335-18-bf2-04.stl",
  "1105335-18-bf2-07.stl", "1105335-18-bf2-10.stl", "1105335-18-bf2-53.stl",
  "1105335-18-bf2-56.stl", "1105335-18-bf2-59.stl", "1105335-19-bf3-02.stl",
  "1105335-19-bf3-05.stl", "1105335-19-bf3-08.stl", "1105335-19-bf3-51.stl",
  "1105335-19-bf3-54.stl", "1105335-19-bf3-57.stl"};

  for (G4int i = 0; i < num_leaves; ++i) {    
    bool found = false;
    dir_length = directory.length();
    for (const auto &entry : fs::directory_iterator(directory)) {
      path_string = entry.path().u8string();
      if ("\\" + mlc_names[i] == path_string.substr(dir_length)) {
        if (path_string[dir_length + 1] == '1') {
            G4int index_value;
            if (path_string[dir_length + 12] == 'a') {
                //create mesh from stl Tessellated
                auto leaf_mesh = CADMesh::TessellatedMesh::FromSTL(path_string);
                //+9mm comes from somewhere CAD geometry centerpoint
                leaf_mesh->SetOffset(9.0, 0.0, 0.0);
                G4VSolid *leaf_solid = leaf_mesh->GetSolid();
                G4LogicalVolume *leaf_lv = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial(fMlc_material_NDS), "leaf_lv", 0, 0, 0);
                visAtt = new G4VisAttributes(G4Colour(0.4, 0.8, 0.3, 0.9));
                leaf_lv->SetVisAttributes(visAtt);
                index_value = std::stoi(path_string.substr(path_string.length() - 6,2));
                index_value += num_leaves/2 - 1;
                NDS120[index_value].index = index_value;
                NDS120[index_value].bank = "BANKA";
                NDS120[index_value].placement = new G4PVPlacement(0, G4ThreeVector(0., 0., fSAD - fCollPos), leaf_lv, "NDS_leaf", fColl_LV, false, 0);
                found = true;
                break;
            } else if (path_string[dir_length + 12] == 'b') {
                auto leaf_mesh = CADMesh::TessellatedMesh::FromSTL(path_string);
                leaf_mesh->SetOffset(191.0, 0.0, 0.0);
                G4VSolid *leaf_solid = leaf_mesh->GetSolid();
                G4LogicalVolume *leaf_lv = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial(fMlc_material_NDS), "leaf_lv", 0, 0, 0);
                visAtt = new G4VisAttributes(G4Colour(0.2, 0.1, 1.0, 1.0));
                leaf_lv->SetVisAttributes(visAtt);
                index_value = std::stoi(path_string.substr(path_string.length() - 6,2));
                index_value -= 1;
                NDS120[index_value].index = index_value;
                NDS120[index_value].bank = "BANKB";
                NDS120[index_value].placement = new G4PVPlacement(0, G4ThreeVector(0., 0., fSAD - fCollPos), leaf_lv, "NDS_leaf", fColl_LV, false, 0);
                found = true;
                break;
            } else if (path_string[dir_length + 14] == 'a') {
                auto leaf_mesh = CADMesh::TessellatedMesh::FromSTL(path_string);
                leaf_mesh->SetOffset(9.0, 0.0, 0.0);
                G4VSolid *leaf_solid = leaf_mesh->GetSolid();
                G4LogicalVolume *leaf_lv = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial(fMlc_material_NDS), "leaf_lv", 0, 0, 0);
                visAtt = new G4VisAttributes(G4Colour(0.4, 0.8, 0.3, 0.9));
                leaf_lv->SetVisAttributes(visAtt);
                index_value = std::stoi(path_string.substr(path_string.length() - 6,2));
                index_value += num_leaves/2 - 1;
                NDS120[index_value].index = index_value;
                NDS120[index_value].bank = "BANKA";
                NDS120[index_value].placement = new G4PVPlacement(0, G4ThreeVector(0., 0., fSAD - fCollPos), leaf_lv, "NDS_leaf", fColl_LV, false, 0);
                found = true;
                break;
            } else if (path_string[dir_length + 14] == 'b') {
                auto leaf_mesh = CADMesh::TessellatedMesh::FromSTL(path_string);
                leaf_mesh->SetOffset(191.0, 0.0, 0.0);
                G4VSolid *leaf_solid = leaf_mesh->GetSolid();
                G4LogicalVolume *leaf_lv = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial(fMlc_material_NDS), "leaf_lv", 0, 0, 0);
                visAtt = new G4VisAttributes(G4Colour(0.2, 0.1, 1.0, 1.0));
                leaf_lv->SetVisAttributes(visAtt);
                index_value = std::stoi(path_string.substr(path_string.length() - 6,2));
                index_value -= 1;
                NDS120[index_value].index = index_value;
                NDS120[index_value].bank = "BANKB";
                NDS120[index_value].placement = new G4PVPlacement(0, G4ThreeVector(0., 0., fSAD - fCollPos), leaf_lv, "NDS_leaf", fColl_LV, false, 0);
                found = true;
                break;
            }  
         }
      }
    } if (!found) {
          G4ExceptionDescription ed;
          ed << "NDS120 MLC STL-file "<< mlc_names[i] <<" not found";
          G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac040", FatalException, ed);
    }
  }
    
}

   
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_TrueBeamDetectorConstruction::BuildHDMLC() {
  //============================================================================
  // HD 120 MLC
  //============================================================================
  //
  // refer to fig 18 of MC data package
  // x is direction of motion of leaves

  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  const G4int num_leaves = 120;
  const G4double delta = 0.01*mm;  // spacing for boolean solids
  const G4double leaf_trap_angle = 0.071*deg;

  // there are 6 leaf types:
  // half target, half iso, quarter target, quarter iso, outboard1, outboard60
  // different depending on which bank

  // for both bank A and B, leaf 1 is at negative y

  // TODO rounded corners for the edges
  // TODO hangers
  // TODO bolt holes
  // TODO length along X
  // TODO data package says leaf gap is 0.047 mm (how is that defined?)
  // TODO measure position of outboard leaf (should be 11?)
  // TODO tongue/groove is at 510 mm from target (approx)
  //   iso quarter: bottom of leaf 544.9 from target

  //---------------------------------------------------------------------------
  // the tip (starts common to all, but the corners aren't) TODO split
  //  by 'tip' I mean the rounded part at the field edge
  //---------------------------------------------------------------------------
  // the main rounding (16 cm radius) (used as intersection)
  // here, create a positive to be subtracted later
  // create a G4Box and subtract the G4Tubs
  // center of box: center of tip
  G4RotationMatrix* tipRot = new G4RotationMatrix();
  tipRot->rotateX(90.*deg);

  // use this to rotate the tip cutout for bank B
  G4RotationMatrix* tipRot_bankB = new G4RotationMatrix();
  tipRot_bankB->rotateZ(pi);

  G4double hd120_tip_center = 50.*mm;
  G4double tip_radius = 160.*mm;

  // center of final object:
  // it is offset by along x by hd120_tip_center
  G4Box* mlc_tip_1 =
    new G4Box("hd120_tip_1", hd120_tip_center, 20.*mm, 50.*mm);
  G4Tubs* mlc_tip_2 =
    new G4Tubs("hd120_tip_2", 0., tip_radius, 10.*mm, 0., twopi);
  G4SubtractionSolid* mlc_tip_3 =
    new G4SubtractionSolid("hd120_tip_3", mlc_tip_1, mlc_tip_2,
      tipRot,
      G4ThreeVector(- (tip_radius - hd120_tip_center), 0., 0.));

  // corners (description in 2d)
  // start with a square
  // quarter circle with radius equal to side of square
  // subtract quarter circle from square to form a negative quarter circle
  // then add the neg quarter circle to the main tip
  //
  // vertical size of leaf is different for quarter/half
  // thus the tip cutout needs to be different

  G4Box* mlc_tip_4 =
    new G4Box("hd120_tip_4", 5.*mm, 5.*mm, 4.3*mm); // tune Z size by vis
  G4Tubs* mlc_tip_5 =
    new G4Tubs("hd120_tip_5", 0., 5.*mm, 10.*mm, 0, 2.*pi);
  G4SubtractionSolid* mlc_tip_6 =
    new G4SubtractionSolid("hd120_tip_6", mlc_tip_4, mlc_tip_5,
      tipRot,
      G4ThreeVector(-5.*mm, 0., -5.*mm));

  G4UnionSolid* mlc_tip =
    new G4UnionSolid("hd120_tip_7", mlc_tip_3, mlc_tip_6,
      0,
      G4ThreeVector(hd120_tip_center - 2.62*mm, 0., 33.5*mm)); //for half leafs
      // tune x pos by vis; 2.6 too small, 2.65 too much
      //G4ThreeVector(round_box_size - 3.*cm + std::tan(sq_angle)*leaf_ztip/2.- std::sin(45*deg)*0.5*mm - 0.01*mm,0, -leaf_ztip/2. -  0.01*mm));

  //---------------------------------------------------------------------------
  // half target leaf
  //---------------------------------------------------------------------------

  G4double leaf_wl   =  2.49*mm;
  G4double leaf_wt   =  0.4 *mm;
  G4double leaf_ztip = 67.0 *mm;
  G4double leaf_zt   = 34.7 *mm;
  G4double leaf_zg   = 34.8 *mm;
  G4double leaf_wg   =  0.4 *mm;
  G4double leaf_width = leaf_wl + leaf_wt;
  G4double leaf_wtip1 = 0.64*mm;   // the amount of material removed
                                   // [wl-wt-wtip in book]

  G4double leaf_zbs  =  3.6 *mm;   // bottom of support rail
  G4double leaf_wbs  =  1.38*mm;   // support rail width
  G4double leaf_zts  =  1.0 *mm;   // top of support rail

  G4double leaf_zl = 1.2*mm;  //ztip-zl in data book

  // the main leaf object
  //G4Box* mlc_htar1 =
  //    new G4Box("hd120_halftar1", fMlc_x_size/2., leaf_width/2.,
  //              leaf_ztip/2.);

  G4Trd* mlc_htar1 =
    new G4Trd("hd120_halftar1", fMlc_x_size/2., fMlc_x_size/2.,
              3.051/2.*mm, 2.720/2.*mm, 66.96/2.*mm);

  //cutouts for (negative) tongue and groove
  // for the half target, zt and zg contain material
  // amount of material to remove is ztip - zt etc
  // half iso, zt and zg material is removed
  G4double zsize1 = leaf_ztip - leaf_zt;
  G4Box* mlc_htar2 =
    new G4Box("hd120_halftar2", fMlc_x_size/2. + delta, leaf_wt,
              zsize1/2.);

  G4RotationMatrix* ht_tongue_rot = new G4RotationMatrix();
  ht_tongue_rot->rotateX(-2.*leaf_trap_angle);

  G4SubtractionSolid* mlc_htar3 = // TODO sign on y displacement
    new G4SubtractionSolid("hd120_halftar3", mlc_htar1, mlc_htar2,
      ht_tongue_rot,
      G4ThreeVector(0., 1.486*mm, -leaf_ztip/2. + zsize1/2.));

  // groove
  zsize1 = leaf_ztip - leaf_zg;
  //G4Box* mlc_htar4 =
  //  new G4Box("hd120_halftar4", fMlc_x_size/2. + delta, leaf_wg,
  //            zsize1/2.);

  G4RotationMatrix* ht_groove_rot = new G4RotationMatrix();
  ht_groove_rot->rotateX(2.*leaf_trap_angle);
  
  G4SubtractionSolid* mlc_htar5 = // TODO sign on y displacement
    new G4SubtractionSolid("hd120_halftar5", mlc_htar3, mlc_htar2,
      ht_groove_rot,
      G4ThreeVector(0., -leaf_width/2., -leaf_ztip/2. + zsize1/2.));

  // support rail  (same side as groove)
  // TODO does this protrude into tip? if not, how far back?
  G4Box* mlc_htar6 =
    new G4Box("hd120_half_tar6", fMlc_x_size/2. + delta, leaf_wbs,
              (leaf_zbs-leaf_zts)/2.);

  G4SubtractionSolid* mlc_htar7 =
    new G4SubtractionSolid("hd120_half_tar7", mlc_htar5, mlc_htar6,
      0, G4ThreeVector(0., -leaf_width/2.,
                        leaf_ztip/2. - (leaf_zbs+leaf_zts)/2.));

  // edge closest to isocenter
  G4Box* mlc_htar8 = 
    new G4Box("hd120_half_tar8", fMlc_x_size/2. + delta, leaf_wtip1, leaf_zl);

  G4SubtractionSolid* mlc_htar9 =
    new G4SubtractionSolid("hd120_half_tar8", mlc_htar7, mlc_htar8, 0,
      G4ThreeVector(0., (leaf_wl-leaf_wg)/2, -leaf_ztip/2.));

  // tip
  G4SubtractionSolid* mlc_htar_10 =
    new G4SubtractionSolid("hd120_halftar_10", mlc_htar9, mlc_tip,
  //G4UnionSolid* mlc_htar_6 =
  //  new G4UnionSolid("hd120_halftar_6", mlc_htar5, mlc_tip,
                     0,
                     //G4ThreeVector(-(160.*mm - fMlc_x_size/2.), 0., 0.));
                     G4ThreeVector(fMlc_x_size/2. - hd120_tip_center, 0., 0.));

  G4SubtractionSolid* mlc_htar_11 =
    new G4SubtractionSolid("hd120_halftar_11", mlc_htar9, mlc_tip,
                     tipRot_bankB,
                     //G4ThreeVector(+(160.*mm - fMlc_x_size/2.), 0., 0.));
                     G4ThreeVector(-fMlc_x_size/2. + hd120_tip_center, 0., 0.));

  // TODO these are class variables; what to do about them?
  fMlc_half_tar_A_LV =
    new G4LogicalVolume(mlc_htar_10, G4Material::GetMaterial(fMlc_material),
    "mlc_hd120_half_target_A_leaf", 0, 0, 0);

  fMlc_half_tar_B_LV =
    new G4LogicalVolume(mlc_htar_11, G4Material::GetMaterial(fMlc_material),
    "mlc_hd120_half_target_B_leaf", 0, 0, 0);

  //---------------------------------------------------------------------------
  // half isocenter leaf
  //---------------------------------------------------------------------------

  leaf_wl =  2.50*mm;
  leaf_wt =  0.4 *mm;
  leaf_wg =  0.4 *mm;
  leaf_zt = 32.5 *mm;
  leaf_zg = 32.9 *mm;
  leaf_width = leaf_wl + leaf_wt;

  leaf_ztip = 67.0 *mm;
  leaf_width = leaf_wl + leaf_wt;
  leaf_wtip1 = 1.35*mm;   // the amount of material removed [wl-wt-wtip in book]

  leaf_zbs  = 66. *mm;   // bottom of support rail
  leaf_wbs  =  1.41*mm;   // support rail width
  leaf_zts  = 63.4 *mm;   // top of support rail

  leaf_zl = 1.2*mm;  //ztip-zl in data book

  // the main leaf object
  G4Trd* mlc_hiso1 =
    new G4Trd("hd120_halftar1", fMlc_x_size/2., fMlc_x_size/2.,
              3.064/2.*mm, 2.733/2.*mm, 66.95/2.*mm);
  
  //cutouts for (negative) tongue and groove
  // half iso, zt and zg material is removed
  zsize1 = leaf_zt;
  G4Box* mlc_hiso2 =
    new G4Box("hd120_halfiso2", fMlc_x_size/2. + delta, leaf_wt,
              (zsize1)/2.);

  G4RotationMatrix* hi_tongue_rot = new G4RotationMatrix();
  hi_tongue_rot->rotateX(-2.*leaf_trap_angle);
  
  G4SubtractionSolid* mlc_hiso3 = // TODO sign on y displacement
    new G4SubtractionSolid("hd120_halfiso3", mlc_hiso1, mlc_hiso2, hi_tongue_rot,
      G4ThreeVector(0., 1.407*mm, +leaf_ztip/2. - zsize1/2.));

  // groove
  zsize1 = leaf_zg;
  G4Box* mlc_hiso4 =
    new G4Box("hd120_halfiso4", fMlc_x_size/2. + delta, leaf_wg,
              zsize1/2.);

  G4RotationMatrix* hi_groove_rot = new G4RotationMatrix();
  hi_groove_rot->rotateX(-2.*leaf_trap_angle);
  
  G4SubtractionSolid* mlc_hiso5 =
    new G4SubtractionSolid("hd120_halfiso5", mlc_hiso3, mlc_hiso4, hi_groove_rot,
      G4ThreeVector(0., -1.407*mm, +leaf_ztip/2. - zsize1/2.));

  // support rail  (same side as groove)
  // TODO does this protrude into tip? if not, how far back?
  G4Box* mlc_hiso6 =
    new G4Box("hd120_half_iso6", fMlc_x_size/2. + delta, leaf_wbs,
              (leaf_zbs-leaf_zts)/2.);

  G4SubtractionSolid* mlc_hiso7 =
    new G4SubtractionSolid("hd120_half_iso7", mlc_hiso5, mlc_hiso6,
      0, G4ThreeVector(0., leaf_width/2.,
                        leaf_ztip/2. - (leaf_zbs+leaf_zts)/2.));

  // edge closest to isocenter
  G4Box* mlc_hiso8 = 
    new G4Box("hd120_half_iso8", fMlc_x_size/2. + delta, leaf_wtip1, leaf_zl);

  G4SubtractionSolid* mlc_hiso9 =
    new G4SubtractionSolid("hd120_half_iso9", mlc_hiso7, mlc_hiso8, 0,
      G4ThreeVector(0., -(leaf_wl-leaf_wg)/2, leaf_ztip/2.));

  // tip
  G4SubtractionSolid* mlc_hiso10 =
    new G4SubtractionSolid("hd120_halfiso_10", mlc_hiso9, mlc_tip,
                     0,
                     G4ThreeVector(fMlc_x_size/2. - hd120_tip_center, 0., 0.));

  G4SubtractionSolid* mlc_hiso11 =
    new G4SubtractionSolid("hd120_halfiso_11", mlc_hiso9, mlc_tip,
                   tipRot_bankB,
                   G4ThreeVector(-(fMlc_x_size/2. - hd120_tip_center),0., 0.));


  // TODO these are class variables; what to do about them?
  fMlc_half_iso_A_LV =
    new G4LogicalVolume(mlc_hiso10, G4Material::GetMaterial(fMlc_material),
    "mlc_hd120_half_iso_A_leaf", 0, 0, 0);

  fMlc_half_iso_B_LV =
    new G4LogicalVolume(mlc_hiso11, G4Material::GetMaterial(fMlc_material),
    "mlc_hd120_half_iso_B_leaf", 0, 0, 0);

  //---------------------------------------------------------------------------
  // quarter target leaf
  //---------------------------------------------------------------------------

  leaf_wl   =  1.22*mm;
  leaf_wt   =  0.4 *mm;
  leaf_wg   =  0.4 *mm;
  leaf_ztip = 67.5 *mm;
  leaf_zt   = 34.88*mm;
  leaf_zg   = 34.88*mm;
  leaf_zbs  =  3.3 *mm;   // bottom of support rail
  leaf_wbs  =  0.89*mm;   // support rail width   // wbs and wts nearly equal
  leaf_zts  =  1.0 *mm;   // top of support rail
  leaf_width = leaf_wl + leaf_wt;
  leaf_wtip1 = 0.32*mm;   // the amount of material removed [wl-wt-wtip in book]
  leaf_zl    = 1.2*mm;   // this is ztip - zl in the data book

  G4Trd* mlc_qtar1 =
    new G4Trd("hd120_quarttar1", fMlc_x_size/2., fMlc_x_size/2.,
              1.702/2.*mm, 1.538/2.*mm, 67.46/2.*mm);

  //cutouts for (negative) tongue and groove
  // for the quarter target, zt and zg contain material
  // amount of material to remove is ztip - zt etc
  zsize1 = leaf_ztip - leaf_zt;
  G4Box* mlc_qtar2 =
    new G4Box("hd120_quarter_tar2", fMlc_x_size/2. + delta, leaf_wt,
              zsize1/2. + delta);

  G4RotationMatrix* qt_tongue_rot = new G4RotationMatrix();
  qt_tongue_rot->rotateX(-leaf_trap_angle);

  G4SubtractionSolid* mlc_qtar3 =
    new G4SubtractionSolid("hd120_quarter_tar3", mlc_qtar1, mlc_qtar2,
      qt_tongue_rot, G4ThreeVector(0., .831*mm, -leaf_ztip/2. + zsize1/2.));

  // groove
  zsize1 = leaf_ztip - leaf_zg;

  G4RotationMatrix* groove_rot = new G4RotationMatrix();
  groove_rot->rotateX(leaf_trap_angle);

  // see note at quarter iso leaf for offset position calc
  G4SubtractionSolid* mlc_qtar5 = // TODO sign on y displacement
    new G4SubtractionSolid("hd120_quarter_tar5", mlc_qtar3, mlc_qtar2,
      groove_rot, G4ThreeVector(0., -.831*mm, -leaf_ztip/2. + zsize1/2.));

  // support rail  (same side as groove)
  // TODO does this protrude into tip? if not, how far back?
  G4Box* mlc_qtar6 =
    new G4Box("hd120_quarter_tar6", fMlc_x_size/2. + delta, leaf_wbs,
              (leaf_zbs-leaf_zts)/2.);

  G4SubtractionSolid* mlc_qtar7 =
    new G4SubtractionSolid("hd120_quarter_tar7", mlc_qtar5, mlc_qtar6,
      0, G4ThreeVector(0., -leaf_width/2.,
          leaf_ztip/2. - (leaf_zbs+leaf_zts)/2.));

  // edge closest to isocenter
  G4Box* mlc_qtar8 = 
    new G4Box("hd120_quarter_tar8", fMlc_x_size/2. + delta, leaf_wtip1,
              leaf_zl);

  G4SubtractionSolid* mlc_qtar9 =
    new G4SubtractionSolid("hd120_tar8", mlc_qtar7, mlc_qtar8, 0,
      G4ThreeVector(0., (leaf_wl-leaf_wg)/2, -leaf_ztip/2.));

  // tip
  G4SubtractionSolid* mlc_qtar10 =
    new G4SubtractionSolid("hd120_quarter_tar10", mlc_qtar9, mlc_tip,
             0, G4ThreeVector((fMlc_x_size/2. - hd120_tip_center), 0., 0.));

  G4SubtractionSolid* mlc_qtar11 =
    new G4SubtractionSolid("hd120_quarter_tar11", mlc_qtar9, mlc_tip,
                 tipRot_bankB,
                 G4ThreeVector(-(fMlc_x_size/2. - hd120_tip_center), 0., 0.));

  // TODO these are class variables; what to do about them?
  fMlc_quarter_tar_A_LV =
    new G4LogicalVolume(mlc_qtar10, G4Material::GetMaterial(fMlc_material),
    "mlc_hd120_quarter_target_A_leaf", 0, 0, 0);

  fMlc_quarter_tar_B_LV =
    new G4LogicalVolume(mlc_qtar11, G4Material::GetMaterial(fMlc_material),
    "mlc_hd120_quarter_target_B_leaf", 0, 0, 0);



  //---------------------------------------------------------------------------
  // quarter isocenter leaf
  //---------------------------------------------------------------------------
  leaf_wl   =  1.22*mm;
  leaf_wt   =  0.4 *mm;
  leaf_wg   =  0.4 *mm;
  leaf_ztip = 67.45*mm;
  //leaf_ztip = 67.5*mm;
  leaf_zt   = 32.95*mm;
  leaf_zg   = 32.83*mm;
  leaf_width = leaf_wl + leaf_wt;
  // keep leaf_zbs and leaf_zts from quarter target leaf

  // the main leaf object
  G4Trd* mlc_qiso1 =
    new G4Trd("hd120_quartiso1", fMlc_x_size/2., fMlc_x_size/2.,
              1.708/2.*mm, 1.542/2.*mm, 67.45/2.*mm);

  //cutouts for (negative) tongue and groove
  // quarter iso, zt and zg material is removed
  zsize1 = leaf_zt;
  G4Box* mlc_qiso2 =
    new G4Box("hd120_quarter_iso2", fMlc_x_size/2. + delta, leaf_wt,
              zsize1/2. + delta);

  G4double qi_trap_angle = 0.071*deg;

  G4RotationMatrix* qi_tongue_rot = new G4RotationMatrix();
  qi_tongue_rot->rotateX(-qi_trap_angle);
 
  // y value of offset calculated based on trapezoid
  // it is edge of trapezoid at midpoint of tongue
  G4SubtractionSolid* mlc_qiso3 =
    new G4SubtractionSolid("hd120_quarter_iso3", mlc_qiso1, mlc_qiso2,
      qi_tongue_rot, G4ThreeVector(0., .791*mm, +leaf_ztip/2. - zsize1/2.));

  // groove
  zsize1 = leaf_zg;
  G4Box* mlc_qiso4 =
    new G4Box("hd120_quarter_iso4", fMlc_x_size/2. + delta, leaf_wg,
              zsize1/2. + delta);

  G4RotationMatrix* qi_groove_rot = new G4RotationMatrix();
  qi_groove_rot->rotateX(qi_trap_angle);
  
  G4SubtractionSolid* mlc_qiso5 = // TODO sign on y displacement
    new G4SubtractionSolid("hd120_quarter_iso5", mlc_qiso3, mlc_qiso4,
      qi_groove_rot, G4ThreeVector(0., -.791*mm, +leaf_ztip/2. - zsize1/2.));

  // support rail  (opposite side as groove)
  // TODO does this protrude into tip? if not, how far back?
  G4Box* mlc_qiso6 =
    new G4Box("hd120_quarter_iso6", fMlc_x_size/2. + delta, leaf_wbs,
              (leaf_zbs-leaf_zts)/2.);

  G4SubtractionSolid* mlc_qiso7 =
    new G4SubtractionSolid("hd120_quarter_iso7", mlc_qiso5, mlc_qiso6,
      0, G4ThreeVector(0., leaf_width/2.,
          -leaf_ztip/2. + (leaf_zbs+leaf_zts)/2.));

  // edge closest to isocenter
  G4Box* mlc_qiso8 = 
    new G4Box("hd120_quarter_iso8", fMlc_x_size/2. + delta,
              leaf_wtip1, leaf_zl);

  G4SubtractionSolid* mlc_qiso9 =
    new G4SubtractionSolid("hd120_iso8", mlc_qiso7, mlc_qiso8, 0,
      G4ThreeVector(0., -(leaf_wl-leaf_wg)/2., leaf_ztip/2.));

  // tip
  G4SubtractionSolid* mlc_qiso10 =
    new G4SubtractionSolid("hd120_quarter_iso8", mlc_qiso9, mlc_tip,
                   0,
                   G4ThreeVector((fMlc_x_size/2. - hd120_tip_center), 0., 0.));

  G4SubtractionSolid* mlc_qiso11 =
    new G4SubtractionSolid("hd120_quarter_iso11", mlc_qiso9, mlc_tip,
                 tipRot_bankB,
                 G4ThreeVector(-(fMlc_x_size/2. - hd120_tip_center), 0., 0.));

  // TODO these are class variables; what to do about them?
  fMlc_quarter_iso_A_LV =
    new G4LogicalVolume(mlc_qiso10, G4Material::GetMaterial(fMlc_material),
    "mlc_hd120_quarter_iso_A_leaf", 0, 0, 0);

  fMlc_quarter_iso_B_LV =
    new G4LogicalVolume(mlc_qiso11, G4Material::GetMaterial(fMlc_material),
    "mlc_hd120_quarter_iso_B_leaf", 0, 0, 0);

  //---------------------------------------------------------------------------
  // outboard leaf #1
  //---------------------------------------------------------------------------
  // TODO top and bottom sides of outboard leaves have more material removed
  leaf_wl   =  2.72*mm;
  leaf_wt   =  1.58*mm;
  leaf_wg   =  0.4 *mm;
  leaf_ztip = 67.0 *mm;
  leaf_zt   = 63.7 *mm;
  leaf_zg   = 34.8 *mm;
  leaf_width = leaf_wl + leaf_wt;

  leaf_wbs =  1.37*mm;
  leaf_zbs =  1.0 *mm;
  leaf_zts =  3.6 *mm;

  leaf_wtip1 = 3.55*mm;  // the amount of material removed [wl-wt-wtip in book]

  // the main leaf object
  G4Box* mlc_outboard1_1 =
      new G4Box("hd120_outboard1_1", fMlc_x_size/2., leaf_width/2.,
                leaf_ztip/2.);

  //cutouts for (negative) tongue and groove
  // for the outboard1, zt and zg contain material(like a target leaf)
  // amount of material to remove is ztip - zt etc
  zsize1 = leaf_ztip - leaf_zt;
  G4Box* mlc_outboard1_2 =
    new G4Box("hd120_outboard1_2", fMlc_x_size/2. + delta, leaf_wt,
              zsize1/2.);

  G4SubtractionSolid* mlc_outboard1_3 =
    new G4SubtractionSolid("hd120_outboard1_3", mlc_outboard1_1, 
      mlc_outboard1_2, 0,
      G4ThreeVector(0., -leaf_width/2., -leaf_ztip/2. + zsize1/2.));

  // groove
  zsize1 = leaf_ztip - leaf_zg;
  G4Box* mlc_outboard1_4 =
    new G4Box("hd120_outboard1_4", fMlc_x_size/2. + delta, leaf_wg,
              zsize1/2.);

  G4SubtractionSolid* mlc_outboard1_5 =
    new G4SubtractionSolid("hd120_outboard1_5", mlc_outboard1_3, 
      mlc_outboard1_4, 0,
      G4ThreeVector(0., leaf_width/2., -leaf_ztip/2. + zsize1/2.));

  // support rail  (opposite side as groove)
  // TODO does this protrude into tip? if not, how far back?
  G4Box* mlc_outboard1_6 =
    new G4Box("hd120_outboard1_6", fMlc_x_size/2. + delta, leaf_wbs,
              (leaf_zts-leaf_zbs)/2.);

  G4SubtractionSolid* mlc_outboard1_7 =
    new G4SubtractionSolid("hd120_outboard1_7", mlc_outboard1_5, 
      mlc_outboard1_6,
      0, G4ThreeVector(0., -leaf_width/2.,
                       leaf_ztip/2. - (leaf_zbs+leaf_zts)/2.));

  // edge closest to isocenter
  G4Box* mlc_outboard1_8 = 
    new G4Box("hd120_outboard1_8", fMlc_x_size/2. + delta, leaf_wtip1, leaf_zl);

  G4SubtractionSolid* mlc_outboard1_9 =
    new G4SubtractionSolid("hd120_outboard1_9", mlc_outboard1_7,
      mlc_outboard1_8, 0,
      G4ThreeVector(0., -(leaf_wl-leaf_wg)/2., -leaf_ztip/2.));

  // tip
  G4SubtractionSolid* mlc_outboard1_10 =
    new G4SubtractionSolid("hd120_outboard1_10", mlc_outboard1_9, mlc_tip,
                   0,
                   G4ThreeVector((fMlc_x_size/2. - hd120_tip_center), 0., 0.));

  G4SubtractionSolid* mlc_outboard1_11 =
    new G4SubtractionSolid("hd120_outboard1_11", mlc_outboard1_9, mlc_tip,
                   tipRot_bankB,
                   G4ThreeVector(-(fMlc_x_size/2. - hd120_tip_center), 0., 0.));

  // TODO these are class variables; what to do about them?
  fMlc_outboard1_A_LV =
    new G4LogicalVolume(mlc_outboard1_10,
      G4Material::GetMaterial(fMlc_material), "mlc_hd120_outboard1_A_leaf",
        0, 0, 0);

  fMlc_outboard1_B_LV =
    new G4LogicalVolume(mlc_outboard1_11,
      G4Material::GetMaterial(fMlc_material), "mlc_hd120_outboard1_B_leaf",
        0, 0, 0);

  //---------------------------------------------------------------------------
  // outboard leaf #60
  //---------------------------------------------------------------------------
  leaf_wl   =  3.90*mm;
  leaf_wt   =  0.4 *mm;
  leaf_wg   =  1.52*mm;
  leaf_ztip = 66.9 *mm;
  leaf_zt   = 32.4 *mm;
  leaf_zg   = 60.9 *mm;
  leaf_width = leaf_wl + leaf_wt;

  // the main leaf object
  G4Box* mlc_outboard60_1 =
      new G4Box("hd120_outboard1", fMlc_x_size/2., (leaf_wl + leaf_wt)/2.,
                leaf_ztip/2.);

  //cutouts for (negative) tongue and groove
  // for outboard60, zg contains material but zt is coutout 
  // (this is unique for this leaf)
  zsize1 = leaf_zt;
  G4Box* mlc_outboard60_2 =
    new G4Box("hd120_outboard60_2", fMlc_x_size/2. + delta, leaf_wt,
              zsize1/2.);

  G4SubtractionSolid* mlc_outboard60_3 =
    new G4SubtractionSolid("hd120_outboard60_3", mlc_outboard60_1,
      mlc_outboard60_2, 0,
      G4ThreeVector(0., -leaf_width/2., +leaf_ztip/2. - zsize1/2.));

  // groove
  zsize1 = leaf_ztip - leaf_zg;
  G4Box* mlc_outboard60_4 =
    new G4Box("hd120_outboard60_4", fMlc_x_size/2. + delta, leaf_wg,
              (zsize1)/2.);

  G4SubtractionSolid* mlc_outboard60_5 =
    new G4SubtractionSolid("hd120_outboard60_5", mlc_outboard60_3,
      mlc_outboard60_4, 0,
      G4ThreeVector(0., -leaf_width/2., -leaf_ztip/2. + zsize1/2.));

  // support rail  (opposite side as groove)
  // TODO does this protrude into tip? if not, how far back?
  G4Box* mlc_outboard60_6 =
    new G4Box("hd120_outboard60_6", fMlc_x_size/2. + delta, leaf_wbs,
              (leaf_zts-leaf_zbs)/2.);

  G4SubtractionSolid* mlc_outboard60_7 =
    new G4SubtractionSolid("hd120_outboard60_7", mlc_outboard60_5, 
      mlc_outboard60_6,
      0, G4ThreeVector(0., leaf_width/2.,
                       -leaf_ztip/2. + (leaf_zbs+leaf_zts)/2.));

  // edge closest to isocenter
  G4Box* mlc_outboard60_8 = 
    new G4Box("hd120_outboard60_8", fMlc_x_size/2. + delta, leaf_wtip1,
              leaf_zl);

  G4SubtractionSolid* mlc_outboard60_9 =
    new G4SubtractionSolid("hd120_outboard60_9", mlc_outboard60_7,
      mlc_outboard60_8, 0,
      G4ThreeVector(0., +(leaf_wl-leaf_wg)/2., leaf_ztip/2.));

  // tip
  G4SubtractionSolid* mlc_outboard60_10 =
    new G4SubtractionSolid("hd120_outboard60_10", mlc_outboard60_9, mlc_tip,
                   0,
                   G4ThreeVector((fMlc_x_size/2. - hd120_tip_center), 0., 0.));

  G4SubtractionSolid* mlc_outboard60_11 =
    new G4SubtractionSolid("hd120_outboard60_11", mlc_outboard60_9, mlc_tip,
                   tipRot_bankB,
                   G4ThreeVector(-(fMlc_x_size/2. - hd120_tip_center), 0., 0.));

  // TODO these are class variables; what to do about them?
  fMlc_outboard60_A_LV = new G4LogicalVolume(mlc_outboard60_10,
      G4Material::GetMaterial(fMlc_material),
      "mlc_hd120_outboard60_A_leaf", 0, 0, 0);

  fMlc_outboard60_B_LV = new G4LogicalVolume(mlc_outboard60_11,
      G4Material::GetMaterial(fMlc_material),
      "mlc_hd120_outboard60_B_leaf", 0, 0, 0);

  //---------------------------------------------------------------------------
  // define the positions
  //   Aug 13, 2018: change to coordinates defined in 100029091-1.pdf
  //                 zfocus is 510 mm
  //   target are even numbers, iso odd
  //---------------------------------------------------------------------------

  // these are the values in the doc
  std::vector<G4double> mlc_pos = {
    -107.200*mm, -102.194*mm, -97.194*mm, -92.196*mm, -87.197*mm, -82.206*mm,
     -77.210*mm,  -72.224*mm, -67.230*mm, -62.249*mm, -57.257*mm, -52.281*mm,
     -47.290*mm,  -42.317*mm, -38.575*mm, -36.090*mm, -33.595*mm, -31.111*mm,
     -28.616*mm,  -26.132*mm, -23.638*mm, -21.154*mm, -18.660*mm, -16.177*mm,
     -13.683*mm,  -11.200*mm,  -8.706*mm,  -6.223*mm,  -3.730*mm,  -1.247*mm,
       1.247*mm,    3.730*mm,   6.223*mm,   8.706*mm,  11.200*mm,  13.683*mm,
      16.177*mm,   18.660*mm,  21.154*mm,  23.638*mm,  26.132*mm,  28.616*mm,
      31.111*mm,   33.595*mm,  36.090*mm,  38.575*mm,  42.317*mm,  47.290*mm,
      52.281*mm,   57.257*mm,  62.249*mm,  67.230*mm,  72.224*mm,  77.210*mm,
      82.206*mm,   87.197*mm,  92.196*mm,  97.194*mm, 102.194*mm, 107.200*mm
  };

  // corrections to iso leaves, set to make leaves centered w.r.t target leaves
  // target leaves look to be in pretty good alignment with lines drawn from 
  //   target to isocenter
  // values determined by visual inspection
  // leaves all need to go out a bit
  // quarter leaves
  mlc_pos[31] += 0.021*mm;
  mlc_pos[33] += 0.051*mm;
  mlc_pos[35] += 0.072*mm;
  mlc_pos[37] += 0.082*mm;
  mlc_pos[39] += 0.123*mm;
  mlc_pos[41] += 0.144*mm;
  mlc_pos[43] += 0.175*mm;
  mlc_pos[45] += 0.190*mm;

  // 29 is unchanged
  mlc_pos[27] -= 0.021*mm;
  mlc_pos[25] -= 0.046*mm;
  mlc_pos[23] -= 0.077*mm;
  mlc_pos[21] -= 0.101*mm;
  mlc_pos[19] -= 0.118*mm;
  mlc_pos[17] -= 0.144*mm;
  mlc_pos[15] -= 0.166*mm;
 
  // half leaves
  mlc_pos[47] += 0.219*mm + 0.06*mm;
  mlc_pos[49] += 0.256*mm + 0.12*mm;
  mlc_pos[51] += 0.292*mm + 0.13*mm;
  mlc_pos[53] += 0.347*mm + 0.15*mm;
  mlc_pos[55] += 0.401*mm + 0.15*mm;
  mlc_pos[57] += 0.420*mm + 0.16*mm;

  mlc_pos[13] -= 0.180*mm;
  mlc_pos[11] -= 0.220*mm;
  mlc_pos[9]  -= 0.250*mm;
  mlc_pos[7]  -= 0.300*mm;
  mlc_pos[5]  -= 0.360*mm;
  mlc_pos[3]  -= 0.420*mm;
  mlc_pos[1]  -= 0.460*mm;

  // outboard leaves
  mlc_pos[0]  -= 1.2*mm;
  mlc_pos[59] += 2.2*mm;
 
  // convert from isocenter to elevation
  for (int i = 0; i < num_leaves/2; ++i) {
    mlc_pos[i] *= 0.51;
  }

  //---------------------------------------------------------------------------
  // PV placements
  //---------------------------------------------------------------------------

  fMlc_PV = new G4PVPlacement*[num_leaves];

  for (size_t i = 0; i < num_leaves; ++i) {
    fMlc_PV[i] = nullptr;
  }
  G4int half_A_target_index    = 0;
  G4int half_A_iso_index       = 0;
  G4int half_B_target_index    = 0;
  G4int half_B_iso_index       = 0;
  G4int quarter_A_target_index = 0;
  G4int quarter_A_iso_index    = 0;
  G4int quarter_B_target_index = 0;
  G4int quarter_B_iso_index    = 0;
  G4PVPlacement* mlc_leaf = nullptr;
  G4RotationMatrix* mlc_rot[num_leaves];

  // as per 100029091-1. this aligns the step of the tongue/groove of 
  // target leaves with 510 mm, and sets the spacing target/iso by eye
  G4double ht_leaf_offset = 1.24*mm;
  G4double hi_leaf_offset = 1.4 *mm;
  G4double qt_leaf_offset = 1.14*mm;
  G4double qi_leaf_offset = 1.3 *mm;
  G4double outbd0_offset  = 2.3 *mm;
  G4double outbd60_offset = 2.3 *mm;

  for (size_t i = 0; i < num_leaves/2; ++i) {
    G4int ii = static_cast<int>(i < 60 ? i : i - 60); // ii is index common to both banks
    mlc_rot[i] = new G4RotationMatrix();
    mlc_rot[i]->rotateX(-std::atan(mlc_pos[ii]/510.*mm));
  }
  // adjustments by visual inspection
  mlc_rot[47]->rotateX(std::atan(0.065*deg));
  mlc_rot[49]->rotateX(std::atan(0.011*deg));
  mlc_rot[51]->rotateX(std::atan(0.013*deg));
  mlc_rot[53]->rotateX(std::atan(0.015*deg));
  mlc_rot[55]->rotateX(std::atan(0.020*deg));
  mlc_rot[57]->rotateX(std::atan(0.021*deg));

  mlc_rot[13]->rotateX( std::atan(0.012*deg));
  //mlc_rot[11]->rotateX(-std::atan(0.00*deg));
  //mlc_rot[9] ->rotateX(-std::atan(0.00*deg));
  //mlc_rot[7] ->rotateX(-std::atan(0.00*deg));
  //mlc_rot[5] ->rotateX(-std::atan(0.00*deg));
  mlc_rot[3] ->rotateX(-std::atan(0.01*deg));
  mlc_rot[1] ->rotateX(-std::atan(0.02*deg));

  mlc_rot[0] ->rotateX(-std::atan(0.18*deg));
  mlc_rot[59]->rotateX( std::atan(0.18*deg));

  for (size_t i = 0; i < num_leaves; ++i) {
    G4int ii = static_cast<int>(i < 60 ? i : i - 60); // ii is index common to both banks
    // BANK A
    if (i < 60) {
      if (i == 0) {
        //outboard1 
        mlc_leaf = new G4PVPlacement(
          mlc_rot[ii],
          G4ThreeVector(0., mlc_pos[ii],
                        fSAD - fCollPos - 510.*mm + outbd0_offset),
          fMlc_outboard1_A_LV, "hd120_outboard1_A_leaf", fColl_LV, false, 0);
      } else if (i == 59) {
        //outboard60
        mlc_leaf = new G4PVPlacement(
          mlc_rot[ii],
          G4ThreeVector(0., mlc_pos[ii],
                        fSAD - fCollPos - 510.*mm - outbd60_offset),
          fMlc_outboard60_A_LV, "hd120_outboard60_A_leaf", fColl_LV, false, 0);
      } else if ( (i < 14 ) || (i > 45) ) {
        //half leaves
        if (!(i%2)) {  // half target
          mlc_leaf = new G4PVPlacement(
            mlc_rot[ii],
            G4ThreeVector(0., mlc_pos[ii],
                          fSAD - fCollPos - 510.*mm + ht_leaf_offset),
            fMlc_half_tar_A_LV, "hd120_half_target_A_leaf", fColl_LV, false,
            half_A_target_index++);
        } else {  // half iso
          mlc_leaf = new G4PVPlacement(
            mlc_rot[ii],
            G4ThreeVector(0., mlc_pos[ii],
                          fSAD - fCollPos - 510.*mm - hi_leaf_offset),
            fMlc_half_iso_A_LV, "hd120_half_iso_A_leaf", fColl_LV, false,
            half_A_iso_index++);
        }
      } else {
        // quarter leaves
        if (!(i%2)) {  //quarter target
          mlc_leaf = new G4PVPlacement(
            mlc_rot[ii],
            G4ThreeVector(0., mlc_pos[ii],
                          fSAD - fCollPos - 510.*mm + qt_leaf_offset),
            fMlc_quarter_tar_A_LV, "hd120_quarter_target_A_leaf", fColl_LV,
            false, quarter_A_target_index++);
        } else { // quarter iso
          mlc_leaf = new G4PVPlacement(
            mlc_rot[ii],
            G4ThreeVector(0., mlc_pos[ii],
                          fSAD - fCollPos - 510.*mm - qi_leaf_offset),
            fMlc_quarter_iso_A_LV, "hd120_quarter_iso_A_leaf", fColl_LV, false,
            quarter_A_iso_index++);
        }
      }
    } else { // end of BANK A
      // BANK B
      if (ii == 0) {
        //outboard
        mlc_leaf = new G4PVPlacement(
          mlc_rot[ii],
          G4ThreeVector(0., mlc_pos[ii],
                        fSAD - fCollPos - 510.*mm + outbd0_offset),
          fMlc_outboard1_B_LV, "hd120_outboard1_B_leaf", fColl_LV, false, 0);
      } else if (ii == 59) {
        //outboard
        mlc_leaf = new G4PVPlacement(
          mlc_rot[ii],
          G4ThreeVector(0., mlc_pos[ii],
                        fSAD - fCollPos - 510.*mm - outbd60_offset),
          fMlc_outboard60_B_LV, "hd120_outboard60_B_leaf", fColl_LV, false, 0);
      } else if ( (ii < 14 ) || (ii > 45) ) {
        //half leaves
        if (!(i%2)) { //half target
          mlc_leaf = new G4PVPlacement(
            mlc_rot[ii],
            G4ThreeVector(0., mlc_pos[ii],
                          fSAD - fCollPos - 510.*mm + ht_leaf_offset),
            fMlc_half_tar_B_LV, "hd120_half_target_B_leaf", fColl_LV, false,
            half_B_target_index++);
        } else { //half iso
          mlc_leaf = new G4PVPlacement(
            mlc_rot[ii],
            G4ThreeVector(0., mlc_pos[ii],
                          fSAD - fCollPos - 510.*mm - hi_leaf_offset),
            fMlc_half_iso_B_LV, "hd120_half_iso_B_leaf", fColl_LV, false,
            half_B_iso_index++);
        }
      } else {
        // quarter leaves
        if (!(i%2)) {  //quarter target
          mlc_leaf = new G4PVPlacement(
            mlc_rot[ii],
            G4ThreeVector(0., mlc_pos[ii],
                          fSAD - fCollPos - 510.*mm + qt_leaf_offset),
            fMlc_quarter_tar_B_LV, "hd120_quarter_target_B_leaf", fColl_LV,
            false, quarter_B_target_index++);
        } else { // quarter iso
          mlc_leaf = new G4PVPlacement(
            mlc_rot[ii],
            G4ThreeVector(0., mlc_pos[ii],
                          fSAD - fCollPos - 510.*mm - qi_leaf_offset),
            fMlc_quarter_iso_B_LV, "hd120_quarter_iso_B_leaf", fColl_LV, false,
            quarter_B_iso_index++);
        }
      }
    } // end of BANK B

    fMlc_PV[i] = mlc_leaf;
    SetMLCHDLeafPosition(static_cast<int>(i), 200.*mm);
  }

  G4VisAttributes* VisAtt_half_tar_A_mlc =
    new G4VisAttributes(G4Colour(0.4, 0.8, 0.3, 0.9));
  G4VisAttributes* VisAtt_half_iso_A_mlc =
    new G4VisAttributes(G4Colour(0.3, 0.9, 0.3, 0.9));
  G4VisAttributes* VisAtt_half_tar_B_mlc =
    new G4VisAttributes(G4Colour(0.4, 0.8, 0.3, 0.9));
  G4VisAttributes* VisAtt_half_iso_B_mlc =
    new G4VisAttributes(G4Colour(0.3, 0.9, 0.3, 0.9));
  fMlc_half_tar_A_LV->SetVisAttributes(VisAtt_half_tar_A_mlc);
  fMlc_half_tar_B_LV->SetVisAttributes(VisAtt_half_tar_B_mlc);
  fMlc_half_iso_A_LV->SetVisAttributes(VisAtt_half_iso_A_mlc);
  fMlc_half_iso_B_LV->SetVisAttributes(VisAtt_half_iso_B_mlc);

  G4VisAttributes* VisAtt_quarter_tar_A_mlc =
    new G4VisAttributes(G4Colour(0.6, 0.5, 0.3, 0.9));
  G4VisAttributes* VisAtt_quarter_iso_A_mlc =
    new G4VisAttributes(G4Colour(0.5, 0.6, 0.3, 0.9));
  G4VisAttributes* VisAtt_quarter_tar_B_mlc =
    new G4VisAttributes(G4Colour(0.6, 0.5, 0.3, 0.9));
  G4VisAttributes* VisAtt_quarter_iso_B_mlc =
    new G4VisAttributes(G4Colour(0.5, 0.6, 0.3, 0.9));
  fMlc_quarter_tar_A_LV->SetVisAttributes(VisAtt_quarter_tar_A_mlc);
  fMlc_quarter_tar_B_LV->SetVisAttributes(VisAtt_quarter_tar_B_mlc);
  fMlc_quarter_iso_A_LV->SetVisAttributes(VisAtt_quarter_iso_A_mlc);
  fMlc_quarter_iso_B_LV->SetVisAttributes(VisAtt_quarter_iso_B_mlc);

  G4VisAttributes* VisAtt_outboard_A_mlc =
    new G4VisAttributes(G4Colour(0.4, 0.4, 0.7, 0.9));
  G4VisAttributes* VisAtt_outboard_B_mlc =
    new G4VisAttributes(G4Colour(0.4, 0.4, 0.7, 0.9));
  fMlc_outboard1_A_LV->SetVisAttributes(VisAtt_outboard_A_mlc);
  fMlc_outboard60_A_LV->SetVisAttributes(VisAtt_outboard_A_mlc);
  fMlc_outboard1_B_LV->SetVisAttributes(VisAtt_outboard_B_mlc);
  fMlc_outboard60_B_LV->SetVisAttributes(VisAtt_outboard_B_mlc);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double TB02_TrueBeamDetectorConstruction::GetMLCLeafPosition(G4int i) {
  // return the actual position of the leaf (ignoring the difference between
  // nominal and actual due to curvature
  G4PVPlacement* leaf = fMlc_PV[i];
  G4double newpos = leaf->GetTranslation().x();
 
  G4double actual;
  if (i < 60) actual =  -newpos - fMlc_x_size/2.;
  else        actual =   newpos - fMlc_x_size/2.;
 
  actual *= 100./51.;
  return actual;
}
