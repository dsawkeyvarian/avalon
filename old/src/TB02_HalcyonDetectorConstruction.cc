/* This file is the CONFIDENTIAL and PROPRIETARY information of
 * Varian Medical Systems, Inc., and is not to be distributed.
 * 
 * Copyright (c) 2017 Varian Medical Systems, Inc.
 * 
 * For information, contact Daren Sawkey  daren.sawkey@varian.com
 */

#include "TB02_BaseDetectorConstruction.hh"
#include "TB02_HalcyonDetectorConstruction.hh"

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
//#include "G4Region.hh"
#include "G4RegionStore.hh"

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
//CHECK THESE!
//alwindow "al_torus"
//"AlPlate"
// plate "Buildup"

VolumeNameAndTraversalFlag HalcyonTable[] = {
{"world", TraversedGeometry::NONE},
{"gantryMother", TraversedGeometry::NONE},
{"vacuum_LV", TraversedGeometry::NONE},
{"backscatter", TraversedGeometry::NONE},
{"BackscatterPlate", TraversedGeometry::NONE},
{"collimatorMother", TraversedGeometry::NONE},
{"TargetButton_LV", TraversedGeometry::TARGET},
{"TargetBlockTop_LV", TraversedGeometry::TARGET},
{"PrimColl", TraversedGeometry::PRIMARY_COLLIMATOR},
{"carousel_bp_LV", TraversedGeometry::NONE},
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
{"mlc_proximalLeaf_a", TraversedGeometry::MLC},
{"mlc_proximalLeaf_b", TraversedGeometry::MLC},
{"mlc_distalLeafA", TraversedGeometry::MLC},
{"mlc_distalLeafB", TraversedGeometry::MLC},
{"mlc_distalFixedLeaf", TraversedGeometry::MLC},
{"mlc_proximalShield", TraversedGeometry::MLC},
{"petgbox", TraversedGeometry::NONE},
{"SecondaryColl", TraversedGeometry::SHIELD_COLLIMATOR},
{"WaterPhantom", TraversedGeometry::NONE},
{"World", TraversedGeometry::NONE},
{"al_torus", TraversedGeometry::NONE},
{"AlPlate", TraversedGeometry::NONE},
{"Buildup", TraversedGeometry::NONE},
};

//=======================================================================
//  TB02_HalcyonDetectorConstruction
//=======================================================================

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TB02_HalcyonDetectorConstruction::TB02_HalcyonDetectorConstruction(G4String gdml_folder,G4bool construct_phantom,G4ThreeVector phantom_position, G4ThreeVector phantom_box_size)
// : G4VUserDetectorConstruction(),
//   fLogicWorld(nullptr),
//   fPhysiWorld(nullptr)

{
  m_gdml_folder = gdml_folder;
  fSAD = 100.*cm;

  fBeamType = "xray";

  fSimulateCollimators = true;
  fJawOffset = 1.5*mm;

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
  
  fMlc_type = "NDS120"; // or "NDS120HD";

  fMlc_PV = nullptr;
  fMlc_x_size = 132.*mm;  // TODO based on HD120 drawing
  fMlc_height    = 80.*mm;  //
  fMlc_tipradius = 234.*mm;

  //H
  fProximalZPos    = 349.008*mm;
  fLeafHLength   =  73.*mm;
  fDistalZPos    = 439.261*mm;

  fNProx = 58;
  fNDist = 56;
  // ---------------------------------------
  // vector of positions
  for (G4int i = 0; i < fNProx; ++i) {
    fProximalPos.push_back(140.*mm);  // default
  }
  // ---------------------------------------
  // vector of positions
  for (G4int i = 0; i < fNDist; ++i) {
    fDistalPos.push_back(140.*mm);  // default
  }



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

  //fFoil1Name      = "6E";
  //fFoil1Material  = "G4_Cu";
  //fFoil1Thickness = 1. *mm;
  //fFoil1Radius    = 10.*mm;
  //fFoil1Position  = 0. *mm;
  //fFoil2Name      = "6E";
  //fFoil2Material  = "G4_Cu";
  //fFoil2Thickness = 1. *mm;
  //fFoil2Radius    = 50.*mm;
  //fFoil2Position  = 0. *mm;
  //fApplicatorName = "None";
  //fBuildCutOut    = false;
  //fCutOutThickness =  13.9*mm;
  //fCutOutBevelFactor = 1.0;
  //fCutOutMaterial = "cerrotru";

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
  fCollPos = -8.*mm;


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

  fVis1 = true;
  fVis2 = true;
  fVis3 = true;
  fVis4 = true;

  fVerbosity = 0;
  fSDManager = G4SDManager::GetSDMpointer();
  fSDManager->SetVerboseLevel(fVerbosity);

  //fMessenger = new TB02_DetectorMessenger(this);

  fTargetRegion = nullptr;

  fOutputFilename = "dc_default";
  fBuildWaterPhantom = construct_phantom;
  fPhantomPosition = phantom_position;
  fPhantomBoxSize = phantom_box_size;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TB02_HalcyonDetectorConstruction::~TB02_HalcyonDetectorConstruction()
{
  //delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* TB02_HalcyonDetectorConstruction::Construct()
{
  //============================================================================
  //      Definitions of Solids, Logical Volumes, Physical Volumes 
  //============================================================================

  BuildWorld();

  BuildGantryMother();
  BuildCollimatorMother();

  BuildVacuum();

  // select electron or photon beam
  if (fBeamType=="electron") {
    //BuildElectronGeometry();
  }
  else if (fBeamType == "xray") {  // This is always true for Halcyon
    BuildPhotonGeometry();
  }
  else {  // fBeamType neither "electron" nor "xray"
    G4ExceptionDescription ed;
    ed << "Detector Construction: Beam type " << fBeamType << "unknown!";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac100",
                FatalException, ed);
  }
  
  if (fUseIonChamber) { 
    BuildIonChamber(129.*mm);
  }

  //if (fSimulateCollimators) {
    //BuildCollimators();
  //}

  //BuildShieldingCollimator();

  if (fSimulateShielding) {
    //BuildShielding();
  } else {
    BuildBackscatterKiller();
  }

  if (fSimulateVault) {
    //BuildVault();
  }
  return fPhysiWorld;
}

void TB02_HalcyonDetectorConstruction::GenerateTraversedTable() {
  m_id_to_traversed = generateVolumeToTraversed(HalcyonTable, sizeof(HalcyonTable) / sizeof(HalcyonTable[0]), fLogicWorld);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::BuildCollimatorMother() {
  //============================================================================
  // mother volume for collimators
  //============================================================================

  G4Tubs* collTub = new G4Tubs("collimator", 0., 30.*cm, 17.*cm, 0, twopi);
  fColl_LV =
    new G4LogicalVolume(collTub,
                        G4Material::GetMaterial("G4_AIR"),
                        "collimatorMother", 0, 0, 0);

  G4VisAttributes* coll_VisAtt = 
    new G4VisAttributes(G4Colour(0., 0.4, 0.4, 0.0));
  fColl_LV->SetVisAttributes(coll_VisAtt);

  fColl =
    new G4PVPlacement(0, G4ThreeVector(0., 0., fCollPos),
      fColl_LV, "CollimatorMother", fGantry_LV, false, 0);

  SetCollimatorRotation(fCollRot);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::BuildVacuum() {
  //============================================================================
  // Vacuum before target
  //============================================================================

  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4String> physVols;
  physVols.push_back("vacuum");
  for (auto pv : physVols) {
    G4VPhysicalVolume* myVol = store->GetVolume(pv, false);
    if (myVol) store->DeRegister(myVol);
  }

  G4LogicalVolumeStore* lv_store = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* lv = lv_store->GetVolume("vacuum_LV", false);
  if (lv) lv_store->DeRegister(lv);

  G4double vacuumLength = 11.*mm;

  G4Tubs* vacuum = 
    new G4Tubs("vacuum", 0.*mm, 1.*cm, vacuumLength/2., 0., twopi);

  fOrbitVacuum_LV =
    new G4LogicalVolume(vacuum,
                        G4Material::GetMaterial("G4_Galactic"),
                        "vacuum_LV", 0, 0, 0);

  G4VisAttributes* vacuum_VisAtt = 
    new G4VisAttributes(G4Colour(0, 0.2, 0.8, 0.2));
  vacuum_VisAtt->SetVisibility(fVis1);
  fOrbitVacuum_LV->SetVisAttributes(vacuum_VisAtt);

  fOrbitVacuum = new G4PVPlacement(nullptr,
        G4ThreeVector(0., 0., fSAD + vacuumLength/2.  - fGantryPos),
        fOrbitVacuum_LV, "vacuum", fGantry_LV, false, 0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::BuildPhotonGeometry(G4bool build) {
  ////==========================================================================
  // Xray beams  (only components specific to them)
  ////==========================================================================

  (void)build;

  BuildTarget();
  BuildPrimaryCollimator();
  BuildBuildupPlate();
  BuildMCBackscatterPlate();
  BuildSecondaryCollimator();
  BuildAlPlate();
  ProximalLeafBank();
  DistalLeafBank();
  BuildDistalFixed();
  BuildProximalShield();
  BuildPlasticMLCBox();
  BuildAlWindow();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::BuildBackscatterKiller(G4bool build) {
  //============================================================================
  // backscatter killer
  // down to top of primary collimator, at z=15 mm
  //============================================================================

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  
  if (!build) { 
    G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
    G4VPhysicalVolume* myVol = store->GetVolume("backscatter", false);
    if (myVol) store->DeRegister(myVol);
    return;
  }
  
  G4VisAttributes* backscatt_VisAtt = 
    new G4VisAttributes(G4Colour(0.2, 0.0, 0.8, 0.1));
  backscatt_VisAtt->SetVisibility(fVis2);

  G4Sphere* backScatter = 
    new G4Sphere("backscatter", 29.*mm, 30.*mm, 0.*deg, 180.*deg, 
                                0.*deg, 360.*deg);
  fBackScatter_LV = new G4LogicalVolume(backScatter,
                  G4Material::GetMaterial("G4_Galactic"),
                  "backscatter", 0, 0, 0);

  G4RotationMatrix* backscatterRot = new G4RotationMatrix();
  backscatterRot->rotateX(-90.0*deg);
  
  fBackScatter = new G4PVPlacement(backscatterRot,
        G4ThreeVector(0., 0., fSAD - fGantryPos - 15.*mm),
        fBackScatter_LV, "backscatter", fGantry_LV, false, 0);

  fBackScatter_LV->SetVisAttributes(backscatt_VisAtt);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::BuildTarget(G4bool build) {
  //============================================================================
  // Target
  //============================================================================

  (void)build;

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  
  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4String> physVols;
  physVols.push_back("TargetButton");
  physVols.push_back("TargetBlock");
  
  for (auto pv : physVols) {
    G4VPhysicalVolume* myVol = store->GetVolume(pv, false);
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

  G4VisAttributes* target_VisAtt = 
    new G4VisAttributes(G4Colour(0.2, 0.2, 0.8, 0.9));
  target_VisAtt->SetVisibility(fVis2);

  G4VisAttributes* target2_VisAtt = 
    new G4VisAttributes(G4Colour(0.5, 0.8, 0.2, 0.9));
  target2_VisAtt->SetVisibility(fVis2);
  
  //------------------------------------------------------------------------
  // Target Button
  //------------------------------------------------------------------------

  G4Tubs* target_button =
    new G4Tubs("TargetButton", 0., 5.*mm, 0.3175*mm, 0., twopi);

  fTargetLow_LV = 
    new G4LogicalVolume(target_button,
                        G4Material::GetMaterial("G4_W"),
                        "TargetButton_LV", 0, 0, 0);

  fTargetLow = new G4PVPlacement(0,
      G4ThreeVector(0., 0., fSAD - fGantryPos - 0.3175*mm),
        fTargetLow_LV, "TargetButton", fGantry_LV, 
        false, 0);

  fTargetLow_LV->SetVisAttributes(target_VisAtt);

  //------------------------------------------------------------------------
  //  Target block top
  //------------------------------------------------------------------------

  G4double bt_thick = 2.692*mm;

  G4Tubs* target_block_top = 
    new G4Tubs("TargetBlockTop", 0., 10.*mm, bt_thick/2., 0., twopi);

  fTargetBlockTopLow_LV = 
    new G4LogicalVolume(target_block_top,
                        G4Material::GetMaterial("COPPER_GLIDCOP"),
                        "TargetBlockTop_LV", 0, 0, 0);

  fTargetBlockTopLow_LV->SetVisAttributes(target2_VisAtt);

  fTargetBlockTopLow = 
    new G4PVPlacement(nullptr, 
          G4ThreeVector(0., 0., fSAD - fGantryPos - (bt_thick/2. + 0.635*mm)),
          fTargetBlockTopLow_LV, "TargetBlockTop", 
          fGantry_LV, false, 0);

  fTargetLow_LV->SetRegion(fTargetRegion);
  fTargetBlockTopLow_LV->SetRegion(fTargetRegion);
  fTargetRegion->AddRootLogicalVolume(fTargetLow_LV);
  fTargetRegion->AddRootLogicalVolume(fTargetBlockTopLow_LV);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::BuildPrimaryCollimator(G4bool build) {
  //============================================================================
  //  Primary Collimator
  //   this has a trapezoid aperture implemented as boolean solid
  //============================================================================

  (void)build;

  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
  G4VPhysicalVolume* myVol = store->GetVolume("PrimColl", false);
  if (myVol) store->DeRegister(myVol);

  G4double pc_length = 94.*mm;

  G4Tubs* pc_tubs = 
    new G4Tubs("primcoll_tub", 0., 100.*mm, pc_length/2., 0., twopi);

  G4double pc_pos        = 15.*mm + pc_length/2.;
  G4double pc_delta = 1.*mm;  // for booleans

  G4double pc_aperture = 141.022 *mm;
  G4double pc_isodist  = 1012.507*mm;
  G4double deltah      = 12.507  *mm;  // pc_isodist - 1000


  G4double pc_trd_x1 = (15. *mm + deltah - pc_delta)/pc_isodist * pc_aperture;
  G4double pc_trd_x2 = (109.*mm + deltah - pc_delta)/pc_isodist * pc_aperture;

  G4Trd* pc_trd = 
    new G4Trd("primcoll_trd", pc_trd_x2, pc_trd_x1, pc_trd_x2, pc_trd_x1,
              (pc_length+pc_delta)/2.);

  G4SubtractionSolid* pc_subt = 
    new G4SubtractionSolid("PrimaryColl", pc_tubs, pc_trd, nullptr, 
                           G4ThreeVector());

  G4LogicalVolume* primColl_LV = 
    new G4LogicalVolume(pc_subt,
                        G4Material::GetMaterial("W95"),
                        "PrimColl", 0, 0, 0);

  G4VisAttributes* VisAtt_PC = 
    new G4VisAttributes(G4Colour(0.8, 0.1, 0.8, 0.5));
  VisAtt_PC->SetVisibility(fVis2);
  primColl_LV->SetVisAttributes(VisAtt_PC);
  
  fPrimColl = new G4PVPlacement(nullptr, 
    G4ThreeVector(0., 0., fSAD - fGantryPos - pc_pos) ,
    primColl_LV, "PrimColl", fGantry_LV, false, 0);

  G4Region* primColl_region = new G4Region("primcoll");
  primColl_LV->SetRegion(primColl_region);
  primColl_region->AddRootLogicalVolume(primColl_LV);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::BuildSecondaryCollimator() {
  //============================================================================
  //  Secondary Collimator
  //   this has a trapezoid aperture implemented as boolean solid
  //   outer dimension larger than actual part--approx other shielding
  //============================================================================
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
  G4VPhysicalVolume* myVol = store->GetVolume("SecondaryColl", false);
  if (myVol) store->DeRegister(myVol);

  G4double sc_length = 140.*mm;

  G4Tubs* sc_tubs = 
    new G4Tubs("seccoll_tub", 0., 100.*mm, sc_length/2., 0., twopi);

  G4double sc_pos = 149.*mm + sc_length/2.;
  G4double delta  = 1.*mm;  // for booleans

  G4double sc_trd_x1 = 45.044/2.*mm;
  G4double sc_trd_x2 = 83.585/2.*mm;

  G4Trd* sc_trd = 
    new G4Trd("seccoll_trd", sc_trd_x2, sc_trd_x1, sc_trd_x2, sc_trd_x1,
              (sc_length + delta)/2.);

  G4SubtractionSolid* sc_subt = 
    new G4SubtractionSolid("SecondaryColl", sc_tubs, sc_trd, nullptr, 
                           G4ThreeVector());

  G4LogicalVolume* secColl_LV = 
    new G4LogicalVolume(sc_subt,
                        G4Material::GetMaterial("W95"),
                        "SecondaryColl", 0, 0, 0);

  G4VisAttributes* VisAtt_SC = 
    new G4VisAttributes(G4Colour(0.7, 0.2, 0.9, 0.5));
  VisAtt_SC->SetVisibility(fVis2);
  secColl_LV->SetVisAttributes(VisAtt_SC);
  
  new G4PVPlacement(nullptr, 
    G4ThreeVector(0., 0., fSAD - fGantryPos - fCollPos - sc_pos) ,
    secColl_LV, "SecondaryColl", fColl_LV, false, 0);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::BuildBuildupPlate() {
  //============================================================================
  //   Buildup plate upstream of monitor chamber
  //============================================================================
  
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4Tubs* buildup_tubs =
    new G4Tubs("buildup_plate", 0., 100.*mm, 0.4*mm, 0., twopi);

  G4LogicalVolume* buildup_LV =
    new G4LogicalVolume(buildup_tubs, G4Material::GetMaterial("Brass"), 
                        "Buildup", 0, 0, 0);

  new G4PVPlacement(nullptr, G4ThreeVector(0., 0., fSAD - fGantryPos - 112.2*mm),
                    buildup_LV, "Buildup", fGantry_LV, false, 0);


  G4VisAttributes* VisAtt_Buildup = 
    new G4VisAttributes(G4Colour(0.1, 0.7, 0.3, 0.5));
  VisAtt_Buildup->SetVisibility(fVis2);
  buildup_LV->SetVisAttributes(VisAtt_Buildup);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::BuildMCBackscatterPlate() {
  //============================================================================
  //   Backscatter plate downstream of monitor chamber
  //============================================================================
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4Tubs* backscat_tubs =
    new G4Tubs("backscatter_plate", 0., 100.*mm, 0.3*mm, 0., twopi);

  G4LogicalVolume* backscatter_LV =
    new G4LogicalVolume(backscat_tubs, G4Material::GetMaterial("SS304"), 
                        "BackscatterPlate", 0, 0, 0);

  new G4PVPlacement(nullptr,
                    G4ThreeVector(0., 0., fSAD - fGantryPos - 145.8*mm),
                    backscatter_LV, "BackscatterPlate", fGantry_LV, false, 0);

  G4VisAttributes* VisAtt_Backscatter = 
    new G4VisAttributes(G4Colour(0.1, 0.7, 0.3, 0.5));
  VisAtt_Backscatter->SetVisibility(fVis2);
  backscatter_LV->SetVisAttributes(VisAtt_Backscatter);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::BuildAlPlate() {
  //============================================================================
  //   Al plate upstream of MLC
  //============================================================================
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4Tubs* alplate_tubs =
    new G4Tubs("al_plate", 0., 100.*mm, 0.4*mm, 0., twopi);

  G4LogicalVolume* alplate_LV =
    new G4LogicalVolume(alplate_tubs, G4Material::GetMaterial("Aluminum5052"), 
                        "AlPlate", 0, 0, 0);

  new G4PVPlacement(nullptr,
                    G4ThreeVector(0., 0., 
                                 fSAD - fGantryPos - fCollPos - 304.2*mm),
                    alplate_LV, "AlPlate", fColl_LV, false, 0);

  G4VisAttributes* VisAtt_Alplate = 
    new G4VisAttributes(G4Colour(0.1, 0.7, 0.3, 0.5));
  VisAtt_Alplate->SetVisibility(fVis2);
  alplate_LV->SetVisAttributes(VisAtt_Alplate);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::ProximalLeafBank() {
  //============================================================================
  //   proximal leaf bank
  //============================================================================

  // the outer envelope of the leaf is actually (according to drawing) a 
  // symmetric trapezoid. Here it is approximated as a trapezoid with one two
  // right angles. According to the drawing, the tongue/groove cutouts are
  // vertical.
  // Leaves fixed according to drawing March 2021 by LL
  // TODO : Add screwholes?

  // where is the Z position defined for each leaf?: center of trapezoid

  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4PhysicalVolumeStore* pv_store = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4String> physVols;
  physVols.push_back("proximalLeafA");
  physVols.push_back("proximalLeafB");

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
  logVols.push_back("mlc_proximalLeaf_a");
  logVols.push_back("mlc_proximalLeaf_b");

  for (it = logVols.begin(); it != logVols.end(); ++it) {
    G4LogicalVolume* myVol = lv_store->GetVolume(*it, false);
    while (myVol) {
      lv_store->DeRegister(myVol);
      myVol = lv_store->GetVolume(*it, false);
    }
  }

  fProximalLeaves.clear();

  //define the center of the leaf (for positioning)
  G4double delta = 0.2*mm;  // for Booleans

  G4double tg = 0.4*mm;  // tongue/groove depth
  G4double rot = 0.569*deg; // angle of sides (half the total)
  //const G4int nleaves = 58;

  G4double leafHeight = 76.7*mm;
  //G4double lh_fact = 1.-fMlc_height/(80.*mm);
  const G4double sq_angle = 7.973*deg;

  // extruded solid is faceted, doesn't seem to work with boolean solids
  // thus need to build leaf as booleans of boxes

  // this is the envelope of leaf. Everything else subtracted.
  //  the center of G4Trap is such that the thickness is the same on both sides of the center. 
  //  Thus the flat surface is at -(x1+x2)/4
  G4Trap* leaf1 = new G4Trap("leaf1", fLeafHLength*2., leafHeight,  
                            4.218*mm, //0.37 = (4.218-3.478)/2
                           3.478*mm);

  ///  *********************** method so far
  //  overall cross-section is a trapezoid
  //  G4Trap has trap section in X, although we want it in Y
  //  thus for the rest of boolean components, need to use that coord. system
  //  component z is actual x; y is actual z; x is actual y
  //      (trap specified as z, y, x)
  //  i.e. actual x -> component z; y -> x; z -> y;

  // subtract out tongue, groove, rails etc.
  // top left of sec b-b
  G4Box* leaf2 = new G4Box("leaf2", tg, leafHeight/2., fLeafHLength+delta);
  G4SubtractionSolid* leaf3 =
    new G4SubtractionSolid("leaf3", leaf1, leaf2, nullptr, 
                           G4ThreeVector(-((3.478*mm+tg)/2. - std::tan(0.37*mm/leafHeight)*(leafHeight-35.987*mm)), 35.987*mm, 0.));

  //G4SubtractionSolid* leaf3 =
  //  new G4SubtractionSolid("leaf3", leaf1, leaf2, nullptr, 
  //                         G4ThreeVector(-1.924*mm, 39.287*mm, 0.));

  G4RotationMatrix* tg_rot = new G4RotationMatrix();
  tg_rot->rotateZ(-rot);

  G4SubtractionSolid* leaf5 = 
    new G4SubtractionSolid("leaf5", leaf3, leaf2, tg_rot,
                           G4ThreeVector((4.218*mm+tg)/2. - std::tan(0.37*mm/leafHeight)*(leafHeight-39.021*mm), -39.021*mm, 0.));

  //Rails:
  G4Box* leaf6 = new G4Box("leaf6", 2.718*mm/2., 2.2*mm/2., fLeafHLength+delta);
  G4SubtractionSolid* leaf7 = 
    new G4SubtractionSolid("leaf7", leaf5, leaf6, tg_rot,
                           G4ThreeVector(1.1*mm/2., -leafHeight/2.+ 1.1*mm + 2.2*mm/2., 0.));
  G4Box* leaf8 = new G4Box("leaf8", 0.778*mm/2., 2.2*mm/2., fLeafHLength+delta);
  G4SubtractionSolid* leaf9 = 
    new G4SubtractionSolid("leaf9", leaf7, leaf8, tg_rot,
                           G4ThreeVector(3.478*mm/2. - 0.778*mm/2., leafHeight/2. - 2.2*mm/2., 0.));
  G4Box* leaf10 = new G4Box("leaf10", 1.018*mm/2., 1.1*mm/2., fLeafHLength+delta);
  G4SubtractionSolid* leaf11 = 
    new G4SubtractionSolid("leaf11", leaf9, leaf10, tg_rot,
                           G4ThreeVector(3.818*mm/2. -1.018*mm/2., -leafHeight/2. + 1.1*mm/2., 0.));                    

  //Screwhole height: -leafHeight/2. + 35.987*mm/2.??
  G4Box* leaf12 = new G4Box("leaf12", 4.5*mm, 1.5*mm, 113*mm/2.);
  G4SubtractionSolid* leaf13 = 
    new G4SubtractionSolid("leaf13", leaf11, leaf12, tg_rot,
                           G4ThreeVector(0., -leafHeight/2. + 35.987*mm/2., 1.5*mm)); 


  // ---------------------------------------
  // tip. start with a box, trim parts out. Then subtract from extruded solid.
  // TODO this is only a circular tip. there are flat parts at the edges.
  // TODO milled groove (sec m)
  G4Box*  tip1 = new G4Box("tip1", 3.*mm, std::sin(sq_angle)*fMlc_tipradius, 50.*mm);
  G4Tubs* tip2 = new G4Tubs("tip2", 0., fMlc_tipradius, 3.1*mm, 0., twopi);
  G4RotationMatrix* tiprot1 = new G4RotationMatrix();
  tiprot1->rotateY(pi/2.);

  G4SubtractionSolid* tip4 =
    new G4SubtractionSolid("tip4", tip1, tip2, tiprot1,
                           G4ThreeVector(0., 0., fMlc_tipradius));
  // flat parts of tip (tipradius - 4.789*mm and z = +- 31.35*mm)
  G4Box*  tip3 = new G4Box("tip3", 3*mm, 10.*mm/2., 5.5*mm/2.); 

  G4RotationMatrix* box_rot1 = new G4RotationMatrix();
  G4RotationMatrix* box_rot2 = new G4RotationMatrix();
  box_rot1->rotateX(sq_angle);
  box_rot2->rotateX(-sq_angle);
  G4double box_x = std::tan(sq_angle)*fMlc_tipradius*std::sin(sq_angle) - std::cos(sq_angle)*(5.5*mm/2.);
  G4double box_z = fMlc_tipradius*std::sin(sq_angle) + 5.5*mm/std::sqrt(2.)*std::sin(pi/4.-sq_angle);
  
  G4UnionSolid* tip5 =
    new G4UnionSolid("tip5", tip4, tip3, box_rot1, G4ThreeVector(0., -box_z, box_x));
  G4UnionSolid* tip6 =
    new G4UnionSolid("tip6", tip5, tip3, box_rot2, G4ThreeVector(0., box_z, box_x));

  //---- corner rounding  (radius 0.5 mm)
  G4double round = 0.5*mm;
  G4Box* tip7 = 
    new G4Box("tip7", round + .01*mm, 5*mm, 0.35*mm);
  G4Tubs* tip8 = 
    new G4Tubs("tip8", 0.*cm, round, 5*mm, 0., twopi);
  G4RotationMatrix* mlc_sub_rot = new G4RotationMatrix();
  mlc_sub_rot->rotateX(90.*deg);
  G4SubtractionSolid* tip9 =
    new G4SubtractionSolid("tip9", tip7, tip8,
      mlc_sub_rot, G4ThreeVector(-round, 0., round));

  G4SubtractionSolid* tip12 =
    new G4SubtractionSolid("tip12", tip7, tip8,
      mlc_sub_rot, G4ThreeVector(-round, 0., -round));

  G4RotationMatrix* mlc_sub_rot_corner = new G4RotationMatrix();
  mlc_sub_rot_corner->rotateZ(90.*deg);
  mlc_sub_rot_corner->rotateY(-90.*deg);
  G4UnionSolid* tip10 =
    new G4UnionSolid("tip10", tip6, tip9, mlc_sub_rot_corner,
      G4ThreeVector(0, leafHeight/2. +  0.01*mm, std::tan(sq_angle)*leafHeight/2.- std::sin(45*deg)*0.5*mm- 0.01*mm)); 
  
  G4UnionSolid* tip11 =
    new G4UnionSolid("tip11", tip10, tip12, mlc_sub_rot_corner,
      G4ThreeVector(0, -leafHeight/2. -  0.01*mm, std::tan(sq_angle)*leafHeight/2. - std::sin(45*deg)*0.5*mm - 0.01*mm)); //std::sin(45*deg)*0.5*mm
  //----


  G4SubtractionSolid* tip6_a =
    new G4SubtractionSolid("tip6_a", leaf5, tip11, nullptr,
                           G4ThreeVector(0., 0., -fLeafHLength));

  G4RotationMatrix* tiprot2 = new G4RotationMatrix();
  tiprot2->rotateY(pi);

  G4SubtractionSolid* tip6_b =
    new G4SubtractionSolid("tip6_b", leaf5, tip11, tiprot2,
                           G4ThreeVector(0., 0., fLeafHLength));

  G4LogicalVolume* leaf1_a_LV =
    new G4LogicalVolume(tip6_a, G4Material::GetMaterial("W95"),
                        "mlc_proximalLeaf_a", 0, 0, 0);

  G4LogicalVolume* leaf1_b_LV =
    new G4LogicalVolume(tip6_b, G4Material::GetMaterial("W95"),
                        "mlc_proximalLeaf_b", 0, 0, 0);

  // ---------------------------------------
  // vectors of leaf rotations
  // Rotations are of form i * angle; but positions aren't

  std::vector<G4double> ypos { 
    -48.850*mm, -45.341*mm, -41.836*mm, -38.335*mm, -34.838*mm,
    -31.345*mm, -27.854*mm, -24.367*mm, -20.881*mm, -17.398*mm,
    -13.916*mm, -10.436*mm,  -6.957*mm,  -3.478*mm,   0.000*mm,
      3.478*mm,   6.957*mm,  10.436*mm,  13.917*mm,  17.399*mm,
     20.883*mm,  24.368*mm,  27.857*mm,  31.348*mm,  34.842*mm,
     38.340*mm,  41.842*mm,  45.347*mm,  48.858*mm };

  // center the leaf
  for (auto& v : ypos) { v += 0.2*mm; }
  
  std::vector<G4RotationMatrix*> proxLeafRot;
 
  for (G4int i = 0; i < fNProx; ++i) {
    G4int ii = (i < fNProx/2 ? i : i-fNProx/2);
    G4RotationMatrix* r = new G4RotationMatrix(); 
    r->rotateX(-pi/2.);
    r->rotateY(-pi/2.);
    r->rotateZ(-(G4float(ii-14)-0.5) * 0.569066*deg);
    proxLeafRot.push_back(r);
  }

  // create the PV placements
  for (G4int i = 0; i < fNProx/2; ++i) {
    G4double xpos = GetActualProximalPosition(i, fProximalPos[i]);
    fProximalLeaves.push_back(
      new G4PVPlacement(proxLeafRot[i], 
            G4ThreeVector(xpos, ypos[i]-0.2*mm,
                          fSAD - fGantryPos - fCollPos - fProximalZPos),
            leaf1_a_LV, "proximalLeafA", fColl_LV, false, i));
  }
  for (G4int i = 0; i < fNProx/2; ++i) {
    G4int ii = fNProx/2 + i;
    G4double xpos = GetActualProximalPosition(ii, fProximalPos[ii]);
    fProximalLeaves.push_back(
      new G4PVPlacement(proxLeafRot[ii], 
            G4ThreeVector(xpos, ypos[i], 
                          fSAD - fGantryPos - fCollPos - fProximalZPos),
            leaf1_b_LV, "proximalLeafB", fColl_LV, false, i));
  }
  
  G4VisAttributes* VisAtt_ProxLeaf = 
    new G4VisAttributes(G4Colour(0.5, 0.5, 0.1, 1.0));
  VisAtt_ProxLeaf->SetVisibility(fVis2);

  leaf1_a_LV->SetVisAttributes(VisAtt_ProxLeaf);
  leaf1_b_LV->SetVisAttributes(VisAtt_ProxLeaf);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double TB02_HalcyonDetectorConstruction::GetActualProximalPosition(G4int index,
                                                              G4double nom) {
  //============================================================================
  //  get the actual leaf position, given the nominal (isocenter)
  //  takes curvature of leaf into account
  //  also converts sign
  //  arguments are leaf index [index], and nominal position (i.e. position at
  //   isocenter) [nom]
  //    TODO  not quite exact yet (also distal)
  //============================================================================
  
  // TODO range check on index
 
  // implementation of S.Scollay's equation
  const G4double zp = 349.*mm;      // source-centre of leaf distance
  const G4double ziso = 1000.*mm;   // source-isocenter distance
  const G4double tiprad = 234.*mm;  // NB. leave this as 234; using 
    // fMlc_tipradius gives incorrect results
  
  // need the y position of the leaf
  //   y = ((index-28)-14) * 3.478*mm 
  G4int ii = (index <= 28 ? index : index-29);   // independent of bank
  // this is an approximation:
  G4double y = G4float(ii-14) * 3.478*mm;

  G4double trig = std::cos(std::atan(y/ziso));
  // note the sign difference from Stuart's calc. He must define x differently
  //                         v  this one
  G4double f = nom * zp/trig 
        + tiprad * (std::sqrt(nom*nom + y*y + ziso*ziso) - ziso/trig);

  f /= ziso/trig;

  f += fLeafHLength;
  if (index > 28) f *= -1.;
 
  //G4cout << "MLC index: " << index << " nominal: " << nom/mm 
  //       << " actual: ";
  //if (index <= 28) G4cout << f-fLeafHLength << G4endl; 
  //else             G4cout << f+fLeafHLength << G4endl; 

  return f;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::SetAllProximalPositions(G4double pos) {
  //============================================================================
  //   set positions of all proximal leaves to same value
  //============================================================================

  for (G4int i = 0; i < 58; ++i) {
    SetProximalPosition(i, pos);
  }
  return;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::SetProximalPosition(
                                                    G4int index, G4double pos) {
  //============================================================================
  //   set positions of proximal leaf
  //============================================================================

  if (pos == fProximalPos[index]) return;
  if (index < 0 || index > 57) {
    G4ExceptionDescription ed;
    ed << "MLC index " << index << " out of range.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac060",
                FatalException, ed);
  }
  fProximalPos[index] = pos;
  G4PVPlacement* leaf = fProximalLeaves[index];
  const G4ThreeVector orig = leaf->GetObjectTranslation();
  
  G4GeometryManager::GetInstance()->OpenGeometry(fColl);
  leaf->SetTranslation(G4ThreeVector(GetActualProximalPosition(index, pos),
                                     orig.y(), orig.z()));
  G4GeometryManager::GetInstance()->CloseGeometry(fColl);

  if (fVerbosity > 1) {
    G4cout << "Proximal MLC leaf: " << index << " set to: " << pos/cm << " cm."
           << G4endl;
  }

//#ifdef G4VIS_USE
//  G4RunManager::GetRunManager()->GeometryHasBeenModified();
//#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::DistalLeafBank() {
  //============================================================================
  //   distal leaf bank
  //    on base of Proximal 
  //============================================================================

  // Leaves fixed according to drawing March 2021 by LL

  // TODO is there interdigitation? if so need to define different leaves for 
  //   each bank

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  
  G4PhysicalVolumeStore* pv_store = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4String> physVols;
  physVols.push_back("distalLeafA");
  physVols.push_back("distalLeafB");

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
  logVols.push_back("mlc_distalLeafA");
  logVols.push_back("mlc_distalLeafB");

  for (it = logVols.begin(); it != logVols.end(); ++it) {
    G4LogicalVolume* myVol = lv_store->GetVolume(*it, false);
    while (myVol) {
      lv_store->DeRegister(myVol);
      myVol = lv_store->GetVolume(*it, false);
    }
  }

  fDistalLeaves.clear();

  //define the center of the leaf (for positioning)
  G4double delta = 0.2*mm;  // for Booleans

  G4double tg = 0.4*mm;  // tongue/groove depth
  G4double rot = 0.569*deg; // angle of sides (half the total)
  const G4double sq_angle = 8.107*deg;
  G4double leafHeight = 76.7*mm;
  G4double LeafHLength = 73.*mm; //SHOULD BE 155mm/2. according to drawings?
  //G4double lh_fact = 1.-fMlc_height/(80.*mm);

  //  as with prox,  actual x -> component z; y -> x; z -> y;

  // this is the envelope of leaf. Everything else subtracted.
  G4Trap* leaf1 =
    new G4Trap("leaf1", LeafHLength*2., leafHeight,
               5.119*mm,
               4.378*mm);

  // subtract out tongue, groove, rails etc.
  // top left sec c-c
  G4Box* leaf2 = new G4Box("leaf2", tg, leafHeight/2., LeafHLength+delta);

  G4SubtractionSolid* leaf3 =
    new G4SubtractionSolid("leaf3", leaf1, leaf2, nullptr,
                           G4ThreeVector(-((4.378*mm+tg)/2. - std::tan(0.3705*mm/leafHeight)*(leafHeight-31.25*mm)), -31.25*mm, 0.));

  G4RotationMatrix* tg_rot = new G4RotationMatrix();
  tg_rot->rotateZ(-rot);

  //the -.12*mm is offset due to angle of edge
  G4SubtractionSolid* leaf5 = 
    new G4SubtractionSolid("leaf3", leaf3, leaf2, tg_rot,
                           //G4ThreeVector(2.536*mm, -25.*mm, 0.));
                           G4ThreeVector((4.378*mm+tg)/2. - std::tan(0.3705*mm/leafHeight)*(leafHeight-34.227*mm), 34.227*mm, 0.));

  //Rails:
  G4Box* leaf6 = new G4Box("leaf6", 2.878*mm/2., 2.2*mm/2., LeafHLength+delta);
  G4SubtractionSolid* leaf7 = 
    new G4SubtractionSolid("leaf7", leaf5, leaf6, tg_rot,
                           G4ThreeVector(-1.1*mm/2., -leafHeight/2.+ 1.1*mm + 2.2*mm/2., 0.));
  G4Box* leaf8 = new G4Box("leaf8", 1.619*mm/2., 2.2*mm/2., LeafHLength+delta);
  G4SubtractionSolid* leaf9 = 
    new G4SubtractionSolid("leaf9", leaf7, leaf8, tg_rot,
                           G4ThreeVector(-5.119*mm/2. + 1.619*mm/2., leafHeight/2. - 2.2*mm/2., 0.));
  G4Box* leaf10 = new G4Box("leaf10", 1.178*mm/2., 1.1*mm/2., LeafHLength+delta);
  G4SubtractionSolid* leaf11 = 
    new G4SubtractionSolid("leaf11", leaf9, leaf10, tg_rot,
                           G4ThreeVector(-4.378*mm/2. + 1.178*mm/2., -leafHeight/2. + 1.1*mm/2., 0.));

  //Screwhole height: -leafHeight/2. + 35.987*mm/2.??
  G4Box* leaf12 = new G4Box("leaf12", 4.5*mm, 2*mm, 137.4*mm/2.);
  G4SubtractionSolid* leaf13 = 
    new G4SubtractionSolid("leaf13", leaf11, leaf12, tg_rot,
                           G4ThreeVector(0., -leafHeight/2. + 31.25*mm/2., 0.8*mm));
  // tail relief cut  // actual y = 0.2, z= 30.
  // doesn't help spacing issue. it is on distal side
  // FIXME positioning
  //G4Box* leaf6 = new G4Box("leaf6", 0.2*mm, 20.*mm, LeafHLength+delta);
  //G4SubtractionSolid* leaf7 = 
  //  new G4SubtractionSolid("leaf7", leaf5, leaf6, tg_rot,
  //    G4ThreeVector(2.136*mm, -40.*mm, 0.));

  // ---------------------------------------
  // tip. start with a box, trim parts out. Then subtract from extruded solid.
  // TODO this is only a circular tip. there are flat parts at the edges.
  G4Box*  tip1 = new G4Box("tip1", 3.*mm, std::sin(sq_angle)*fMlc_tipradius, 50.*mm);
  G4Tubs* tip2 = new G4Tubs("tip2", 0., fMlc_tipradius, 3.1*mm, 0., twopi);
  G4RotationMatrix* tiprot1 = new G4RotationMatrix();
  tiprot1->rotateY(pi/2.);

  G4SubtractionSolid* tip4 =
    new G4SubtractionSolid("tip4", tip1, tip2, tiprot1,
                           G4ThreeVector(0., 0., fMlc_tipradius));

  // flat parts of tip (tipradius - 3.7*mm and z = +- 31.35*mm) & 8.107*deg
  G4Box*  tip3 = new G4Box("tip3", 3*mm, 3.7*mm, 5.5*mm/2.);
  
  G4RotationMatrix* box_rot1 = new G4RotationMatrix();
  G4RotationMatrix* box_rot2 = new G4RotationMatrix();
  box_rot1->rotateX(sq_angle);
  box_rot2->rotateX(-sq_angle);
  
  G4double box_x = std::tan(sq_angle)*fMlc_tipradius*std::sin(sq_angle) - std::cos(sq_angle)*(5.5*mm/2.);
  G4double box_z = fMlc_tipradius*std::sin(sq_angle) + 5.5*mm/std::sqrt(2.)*std::sin(pi/4.-sq_angle);


  G4UnionSolid* tip5 =
    new G4UnionSolid("tip5", tip4, tip3, box_rot1, G4ThreeVector(0., -box_z, box_x));
  G4UnionSolid* tip6 =
    new G4UnionSolid("tip6", tip5, tip3, box_rot2, G4ThreeVector(0., box_z, box_x));

  //---- corner rounding  (radius 0.5 mm)
  G4double round = 0.5*mm;
  G4Box* tip7 = 
    new G4Box("tip7", round + .01*mm, 5*mm, 0.35*mm);
  G4Tubs* tip8 = 
    new G4Tubs("tip8", 0.*cm, round, 5*mm, 0., twopi);
  G4RotationMatrix* mlc_sub_rot = new G4RotationMatrix();
  mlc_sub_rot->rotateX(90.*deg);
  G4SubtractionSolid* tip9 =
    new G4SubtractionSolid("tip9", tip7, tip8,
      mlc_sub_rot, G4ThreeVector(-round, 0., round));

  G4SubtractionSolid* tip12 =
    new G4SubtractionSolid("tip12", tip7, tip8,
      mlc_sub_rot, G4ThreeVector(-round, 0., -round));

  G4RotationMatrix* mlc_sub_rot_corner = new G4RotationMatrix();
  mlc_sub_rot_corner->rotateZ(90.*deg);
  mlc_sub_rot_corner->rotateY(-90.*deg);
  G4UnionSolid* tip10 =
    new G4UnionSolid("tip10", tip6, tip9, mlc_sub_rot_corner,
      G4ThreeVector(0,leafHeight/2. +  0.01*mm , std::tan(sq_angle)*leafHeight/2.- std::sin(45*deg)*0.5*mm- 0.01*mm)); 
  
  G4UnionSolid* tip11 =
    new G4UnionSolid("tip11", tip10, tip12, mlc_sub_rot_corner,
      G4ThreeVector(0, -leafHeight/2. -  0.01*mm, std::tan(sq_angle)*leafHeight/2.- std::sin(45*deg)*0.5*mm- 0.01*mm)); //std::sin(45*deg)*0.5*mm

  //----

  G4SubtractionSolid* tip6_a =
    new G4SubtractionSolid("tip6_a", leaf5, tip11, nullptr,
                           G4ThreeVector(0., 0., -LeafHLength));

  G4RotationMatrix* tiprot2 = new G4RotationMatrix();
  tiprot2->rotateY(pi);

  G4SubtractionSolid* tip6_b =
    new G4SubtractionSolid("tip6_b", leaf5, tip11, tiprot2,
                           G4ThreeVector(0., 0., LeafHLength));
  
  G4LogicalVolume* leaf1_a_LV =
    new G4LogicalVolume(tip6_a, G4Material::GetMaterial("W95"),
                        "mlc_distalLeafA", 0, 0, 0);
  
  G4LogicalVolume* leaf1_b_LV =
    new G4LogicalVolume(tip6_b, G4Material::GetMaterial("W95"),
                        "mlc_distalLeafB", 0, 0, 0);

  // ---------------------------------------
  // vectors of leaf rotations

  std::vector<G4double> ypos {
   -59.286*mm, -54.871*mm, -50.461*mm, -46.056*mm, -41.656*mm,
   -37.260*mm, -32.868*mm, -28.479*mm, -24.093*mm, -19.709*mm,
   -15.327*mm, -10.947*mm,  -6.567*mm,  -2.189*mm,   2.189*mm,
     6.567*mm,  10.947*mm,  15.327*mm,  19.709*mm,  24.093*mm,
    28.479*mm,  32.868*mm,  37.260*mm,  41.656*mm,  46.056*mm,
    50.461*mm,  54.871*mm,  59.286*mm
  };

  std::vector<G4RotationMatrix*> distalLeafRot;
  for (G4int i = 0; i < fNDist; ++i) {
    G4int ii = (i < fNDist/2 ? i : i-fNDist/2);
    G4RotationMatrix* r = new G4RotationMatrix();
    r->rotateX(-pi/2.);
    r->rotateY(-pi/2.);
    //r->rotateZ((13.5-G4float(ii)) * 0.569351*deg);
    r->rotateZ((13.5-G4float(ii)) * 0.569351*deg);
    distalLeafRot.push_back(r);
  }

  // create the PV placements
  for (G4int i = 0; i < fNDist/2; ++i) {
    G4double xpos = GetActualDistalPosition(i, fDistalPos[i]);
    fDistalLeaves.push_back(
      new G4PVPlacement(distalLeafRot[i], 
            G4ThreeVector(xpos, ypos[i],
                          fSAD - fGantryPos - fCollPos - fDistalZPos),
            leaf1_a_LV, "distalLeafA", fColl_LV, false, i));
  }
  for (G4int i = 0; i < fNDist/2; ++i) {
    G4int ii = fNDist/2 + i;
    G4double xpos = GetActualDistalPosition(ii, fDistalPos[ii]);
    fDistalLeaves.push_back(
      new G4PVPlacement(distalLeafRot[ii], 
            G4ThreeVector(xpos, ypos[i],
                          fSAD - fGantryPos - fCollPos - fDistalZPos),
            leaf1_b_LV, "distalLeafB", fColl_LV, false, i));
  }
 
  G4VisAttributes* VisAtt_DistalLeaf = 
    new G4VisAttributes(G4Colour(0.0, 0.05, 0.95, 0.9));
  VisAtt_DistalLeaf->SetVisibility(fVis2);

  leaf1_a_LV->SetVisAttributes(VisAtt_DistalLeaf);
  leaf1_b_LV->SetVisAttributes(VisAtt_DistalLeaf);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::BuildDistalFixed(void) {
  //============================================================================
  // fixed leaves on field edges of distal
  //  take the distal leaf, make bigger
  //============================================================================
  G4double delta = 0.2*mm;  // for Booleans

  G4double tg = 0.4*mm;  // tongue/groove depth
  G4double rot = 0.569*deg; // angle of sides (half the total)

  // this is the envelope of leaf. Everything else subtracted.
  G4Trap* leaf1 =
    new G4Trap("leaf1", fLeafHLength*4., 80.*mm, 5.119*mm, 3.976*mm);

  // top left
  G4Box* leaf2 = new G4Box("leaf2", tg, 28.*mm, fLeafHLength*4.+delta);

  G4SubtractionSolid* leaf3 =
    new G4SubtractionSolid("leaf3", leaf1, leaf2, nullptr,
                           G4ThreeVector(-2.136*mm, 25.*mm, 0.));

  G4RotationMatrix* tg_rot = new G4RotationMatrix();
  tg_rot->rotateZ(-rot);

  G4SubtractionSolid* leaf5 = 
    new G4SubtractionSolid("leaf3", leaf3, leaf2, tg_rot,
                           G4ThreeVector(2.536*mm, -25.*mm, 0.));

  G4LogicalVolume* leaf1_LV =
    new G4LogicalVolume(leaf5, G4Material::GetMaterial("W95"),
                        "mlc_distalFixedLeaf", 0, 0, 0);

  G4double ypos = -63.707*mm;

  G4RotationMatrix* fixedLeafRot1 = new G4RotationMatrix();
  fixedLeafRot1->rotateX(-pi/2.);
  fixedLeafRot1->rotateY(-pi/2.);
  fixedLeafRot1->rotateZ(8.256*deg);
  G4RotationMatrix* fixedLeafRot2 = new G4RotationMatrix();
  fixedLeafRot2->rotateX(-pi/2.);
  fixedLeafRot2->rotateY(-pi/2.);
  fixedLeafRot2->rotateZ(-8.256*deg);

  new G4PVPlacement(fixedLeafRot1, 
        G4ThreeVector(0., ypos,
                      fSAD - fGantryPos - fCollPos - fDistalZPos),
        leaf1_LV, "fixedDistalLeaf", fColl_LV, false, 0);
 
  new G4PVPlacement(fixedLeafRot2, 
        G4ThreeVector(0., -ypos,
                      fSAD - fGantryPos - fCollPos - fDistalZPos),
        leaf1_LV, "fixedDistalLeaf", fColl_LV, false, 1);


  G4VisAttributes* VisAtt_FixedDistalLeaf = 
    new G4VisAttributes(G4Colour(0.5, 0.3, 0.3, 0.5));
  VisAtt_FixedDistalLeaf->SetVisibility(fVis2);

  leaf1_LV->SetVisAttributes(VisAtt_FixedDistalLeaf);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double TB02_HalcyonDetectorConstruction::GetActualDistalPosition(G4int index,
                                                              G4double nom) {
  //============================================================================
  //  get the actual leaf position, given the nominal (isocenter)
  //  takes curvature of leaf into account
  //  also converts sign
  //   TODO not exat yet
  //============================================================================
  
  // TODO range check on index

  // implementation of S.Scollay's equation
  const G4double zd = 439.*mm;      // source-centre of leaf distance
  const G4double ziso = 1000.*mm;   // source-isocenter distance
  const G4double tiprad = 234.*mm;  // NB. leave this as 234; using 
    // fMlc_tipradius gives incorrect results
  
  // need the y position of the leaf
  //   y = ((index-28)-14) * 3.478*mm 
  G4int ii = (index <= 27 ? index : index-28);   // independent of bank
  // this is an approximation:
  G4double y = (G4float(ii)-13.5) * 4.378*mm;

  G4double trig = std::cos(std::atan(y/ziso));
  // note the sign difference from Stuart's calc. He must define x differently
  //                         v  this one
  G4double f = nom * zd/trig 
          + tiprad * (std::sqrt(nom*nom + y*y + ziso*ziso) - ziso/trig);

  f /= ziso/trig;

  f += fLeafHLength;
  if (index > 27) f *= -1.;

  //G4cout << "MLC index: " << index << " nominal: " << nom/mm 
  //       << " actual: " ;
  //if (index <= 27) G4cout << f-fLeafHLength << G4endl; 
  //else             G4cout << f+fLeafHLength << G4endl; 

  return f;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::SetDistalPosition(G4int index, G4double pos) {
  //============================================================================
  //   set positions of distal leaf
  //============================================================================

  if (pos == fDistalPos[index]) return;
  if (index < 0 || index > 55) {
    G4ExceptionDescription ed;
    ed << "MLC index " << index << " out of range.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac061",
                FatalException, ed);
  }
  fDistalPos[index] = pos;
  G4PVPlacement* leaf = fDistalLeaves[index];
  const G4ThreeVector orig = leaf->GetObjectTranslation();
  
  G4GeometryManager::GetInstance()->OpenGeometry(fColl);
  leaf->SetTranslation(G4ThreeVector(GetActualDistalPosition(index, pos),
                                     orig.y(), orig.z()));
  G4GeometryManager::GetInstance()->CloseGeometry(fColl);

  if (fVerbosity > 1) {
    G4cout << "Distal MLC leaf: " << index << " set to: " << pos/cm << " cm."
           << G4endl;
  }

//#ifdef G4VIS_USE
//  G4RunManager::GetRunManager()->GeometryHasBeenModified();
//#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::SetAllDistalPositions(G4double pos) {
  //============================================================================
  //   set positions of all distal leaves to same value
  //============================================================================

  for (G4int i = 0; i < 56; ++i) {
    SetDistalPosition(i, pos);
  }
  return;
}

// FIXME add the fixed distal leaves
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::BuildProximalShield() {
  //============================================================================
  //  shielding to define field at edges of proximal leaves 
  //============================================================================

  G4Box* sbox = new G4Box("sbox", 100.*mm, 5.*mm, 42.*mm);

  G4LogicalVolume* sbox_LV =
    new G4LogicalVolume(sbox, G4Material::GetMaterial("W95"),
                        "mlc_proximalShield", 0, 0, 0);

  G4RotationMatrix* sbox_rot_1 = new G4RotationMatrix();
  sbox_rot_1->rotateX(-std::atan(.145));
  G4RotationMatrix* sbox_rot_2 = new G4RotationMatrix();
  sbox_rot_2->rotateX(std::atan(.145));

  new G4PVPlacement(sbox_rot_1,
    G4ThreeVector(0., 57.0*mm,
                  fSAD - fGantryPos - fCollPos - fProximalZPos),
    sbox_LV, "proximalLeafShield", fColl_LV, false, 0);

  new G4PVPlacement(sbox_rot_2,
    G4ThreeVector(0., -57.0*mm,
                  fSAD - fGantryPos - fCollPos - fProximalZPos),
    sbox_LV, "proximalLeafShield", fColl_LV, false, 1);

  G4VisAttributes* VisAtt_ProximalShield = 
    new G4VisAttributes(G4Colour(0.8, 0.0, 0.3, 0.2));
  VisAtt_ProximalShield->SetVisibility(fVis2);
  sbox_LV->SetVisAttributes(VisAtt_ProximalShield);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::BuildPlasticMLCBox() {
  //============================================================================
  //  PETG above/below MLC
  //============================================================================

  G4Box* petgbox = new G4Box("petg_box", 100.*mm, 100.*mm, 0.75*mm);

  G4LogicalVolume* petgbox_LV =
    new G4LogicalVolume(petgbox, G4Material::GetMaterial("PETG"), "petgbox",
                        0, 0, 0);

  new G4PVPlacement(nullptr, G4ThreeVector(0., 0.,
                              fSAD - fGantryPos - fCollPos - 475.55*mm),
                    petgbox_LV, "petgbox1", fColl_LV, false, 0);

  new G4PVPlacement(nullptr, G4ThreeVector(0., 0.,
                              fSAD - fGantryPos - fCollPos - 295.5*mm),
                    petgbox_LV, "petgbox1", fColl_LV, false, 1);

  G4VisAttributes* VisAtt_petg = 
    new G4VisAttributes(G4Colour(0.1, 0.7, 0.3, 0.2));
  VisAtt_petg->SetVisibility(fVis2);
  petgbox_LV->SetVisAttributes(VisAtt_petg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TB02_HalcyonDetectorConstruction::BuildAlWindow() {
  //============================================================================
  //  Al torus
  //============================================================================

  G4Box* albox = new G4Box("al_box", 200.*mm, 200.*mm, 1.4*mm);

  // FIXME correct material
  G4LogicalVolume* albox_LV =
    new G4LogicalVolume(albox, G4Material::GetMaterial("Aluminum6061"),
                        "al_torus", 0, 0, 0);

  new G4PVPlacement(nullptr, G4ThreeVector(0., 0., fSAD - fGantryPos - 499.*mm),
                  albox_LV, "al_torus", fGantry_LV, false, 0);

  G4VisAttributes* VisAtt_AlTorus = 
    new G4VisAttributes(G4Colour(0.8, 0.0, 0.3, 0.2));
  VisAtt_AlTorus->SetVisibility(fVis2);
  albox_LV->SetVisAttributes(VisAtt_AlTorus);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_HalcyonDetectorConstruction::ConstructSDandField() {
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
