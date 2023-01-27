/* This file is the CONFIDENTIAL and PROPRIETARY information of
 * Varian Medical Systems, Inc., and is not to be distributed.
 * 
 * Copyright (c) 2017 Varian Medical Systems, Inc.
 * 
 * For information, contact Daren Sawkey  daren.sawkey@varian.com
 */

#include "TB02_BaseDetectorConstruction.hh"

//#include "TB02_DetectorMessenger.hh"
//#include "TB02_PhspWorldConstruction.hh"

//#include "TB02_PhaseSpaceWriter.hh"
//#include "G4PSDirectionFlag.hh"

//#include "G4PSDoseDeposit3D.hh"
#include "G4PSDoseDeposit.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
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

//#include "G4GDMLParser.hh"

//#include "G4PVParameterised.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4ios.hh"

void walkVolumeHierarchy(G4LogicalVolume* volume, const std::unordered_map<std::string, TraversedGeometry>& name_to_traversed,
  std::unordered_map<int, TraversedGeometry>& id_to_traversed, std::vector<std::string>& unlisted_names) {
  const int id = volume->GetInstanceID();
  std::string name = volume->GetName();
  const auto result = name_to_traversed.find(volume->GetName());
  const bool contains = result != name_to_traversed.end();
  if (contains) {
    const TraversedGeometry traversed = result->second;
    id_to_traversed[id] = traversed;
  }
  else {
    unlisted_names.push_back(name);
  }
  const size_t num = volume->GetNoDaughters();
  for (size_t i = 0; i < num; ++i) {
    const auto daughter = volume->GetDaughter((int)i)->GetLogicalVolume();
    //std::cout << daughter->GetName() << " " << daughter->GetInstanceID() << std::endl;
    walkVolumeHierarchy(daughter, name_to_traversed, id_to_traversed, unlisted_names);
  }
}

std::unordered_map<int, TraversedGeometry> generateVolumeToTraversed(const VolumeNameAndTraversalFlag* _array, const size_t _array_size, G4LogicalVolume* _root_volume) {
  std::unordered_map<int, TraversedGeometry> id_to_traversed;
  std::unordered_map<std::string, TraversedGeometry> name_to_traversed;
  std::vector<std::string> unlisted_names;
  for (size_t i = 0; i < _array_size; ++i) {
    name_to_traversed[_array[i].name] = _array[i].traversed;
  }
  walkVolumeHierarchy(_root_volume, name_to_traversed, id_to_traversed, unlisted_names);
  if (unlisted_names.size() > 0) {
    G4cerr << "ParticleTracking Initialization Error" << G4endl;
    G4cerr << "Unable to proceed: Following geometry names are not listed in the traversal geometry table" << G4endl;
    for (size_t i = 0; i < unlisted_names.size(); ++i) {
      G4cerr << "\t" << unlisted_names[i] << G4endl;
    }
    G4cerr << "Stopping now" << G4endl;
    throw std::runtime_error("No all geometry names are in traversal table of the detector");
  }
  return id_to_traversed;
}


//=======================================================================
//  TB02_BaseDetectorConstruction
//=======================================================================

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TB02_BaseDetectorConstruction::TB02_BaseDetectorConstruction()
 : G4VUserDetectorConstruction(),
   fLogicWorld(nullptr),
   fPhysiWorld(nullptr)

{
  fSAD = 100.*cm;

  fBeamType = "xray";

  fSimulateCollimators = true;
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

  fMlc_PV = nullptr;
  fMlc_x_size = 132.*mm;  // TODO based on HD120 drawing
  fMlc_half_tar_A_LV = nullptr;
  fMlc_half_tar_B_LV = nullptr;
  fMlc_half_iso_A_LV = nullptr;
  fMlc_half_iso_B_LV = nullptr;
  fMlc_full_A_LV = nullptr;
  fMlc_full_B_LV = nullptr;
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

  fEfoil1_tubs = nullptr;
  fEfoil1_LV   = nullptr;
  fEfoil1      = nullptr;
  fEfoil2_LV   = nullptr;
  fEfoil2      = nullptr;

  // srs cone
  fCone_size  = 4.*mm;   // nominal diameter
  fCone_Zpos  = 37.5*cm;   // relative to isocenter, of proximal surface
  fCone_thick = 110.*mm;
  fCone_OuterRadius = 32.*mm;
  fBuildCone = false;

  // TSI scatterer
  fTSIPlateUse = false;
  fTSIPlateThickness = 1.*mm;
  fTSIPlateSideLength = 200.*mm;
  fTSIPlatePosition = 40.*cm;
  fTSIPlateMaterial = "G4_Al";
  

  fVis1 = true;
  fVis2 = true;
  fVis3 = true;
  fVis4 = true;

  fVerbosity = 2;
  fSDManager = G4SDManager::GetSDMpointer();
  //fSDManager->SetVerboseLevel(fVerbosity);

  //fMessenger = new TB02_DetectorMessenger(this);

  fTargetRegion = nullptr;

  fOutputFilename = "dc_default";
  DefineMaterials();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TB02_BaseDetectorConstruction::~TB02_BaseDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::DefineMaterials()
{
  //=====================
  // Material Definitions
  //=====================
  //  
  //-------- NIST Materials ----------------------------------------------------
  //  Material Information imported from NIST database.
  //
  G4NistManager* NISTman = G4NistManager::Instance();

  NISTman->FindOrBuildMaterial("G4_Be");
  NISTman->FindOrBuildMaterial("G4_C");
  NISTman->FindOrBuildMaterial("G4_Cu");
  NISTman->FindOrBuildMaterial("G4_W");
  NISTman->FindOrBuildMaterial("G4_Ta");
  NISTman->FindOrBuildMaterial("G4_Al");
  NISTman->FindOrBuildMaterial("G4_Pb");
  NISTman->FindOrBuildMaterial("G4_Ti");
  G4Element* elAl = NISTman->FindOrBuildElement("Al");
  G4Element* elAu = NISTman->FindOrBuildElement("Au");
  G4Element* elBi = NISTman->FindOrBuildElement("Bi");
  G4Element* elC  = NISTman->FindOrBuildElement("C" );
  G4Element* elCd = NISTman->FindOrBuildElement("Cd");
  G4Element* elCl = NISTman->FindOrBuildElement("Cl");
  G4Element* elCr = NISTman->FindOrBuildElement("Cr");
  G4Element* elCu = NISTman->FindOrBuildElement("Cu");
  G4Element* elFe = NISTman->FindOrBuildElement("Fe");
  G4Element* elH  = NISTman->FindOrBuildElement("H" );
  G4Element* elMg = NISTman->FindOrBuildElement("Mg");
  G4Element* elMn = NISTman->FindOrBuildElement("Mn");
  G4Element* elN  = NISTman->FindOrBuildElement("N" );
  G4Element* elNa = NISTman->FindOrBuildElement("Na");
  G4Element* elNi = NISTman->FindOrBuildElement("Ni");
  G4Element* elO  = NISTman->FindOrBuildElement("O" );
  G4Element* elP  = NISTman->FindOrBuildElement("P" );
  G4Element* elPb = NISTman->FindOrBuildElement("Pb");
  G4Element* elS  = NISTman->FindOrBuildElement("S" );
  G4Element* elSb = NISTman->FindOrBuildElement("Sb");
  G4Element* elSi = NISTman->FindOrBuildElement("Si");
  G4Element* elSn = NISTman->FindOrBuildElement("Sn");
  G4Element* elTi = NISTman->FindOrBuildElement("Ti");
  G4Element* elV  = NISTman->FindOrBuildElement("V" );
  G4Element* elW  = NISTman->FindOrBuildElement("W" );
  G4Element* elZn = NISTman->FindOrBuildElement("Zn");


  G4double density;
  G4int ncomponents;
  G4String name;

  NISTman->FindOrBuildMaterial("G4_AIR");
  NISTman->FindOrBuildMaterial("G4_WATER");
  NISTman->FindOrBuildMaterial("G4_MYLAR");
  NISTman->FindOrBuildMaterial("G4_NYLON-6-6");
  NISTman->FindOrBuildMaterial("G4_Galactic");

  // especially for phantom
  NISTman->FindOrBuildMaterial("G4_LUNG_ICRP");
  
  NISTman->FindOrBuildMaterial("G4_B-100_BONE");
  NISTman->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  NISTman->FindOrBuildMaterial("G4_BONE_CORTICAL_ICRP");
  
  NISTman->FindOrBuildMaterial("G4_BRAIN_ICRP");
  NISTman->FindOrBuildMaterial("G4_BLOOD_ICRP");
  NISTman->FindOrBuildMaterial("G4_EYE_LENS_ICRP");
  NISTman->FindOrBuildMaterial("G4_TESTIS_ICRP");
  NISTman->FindOrBuildMaterial("G4_MUSCLE_SKELETAL_ICRP");
  NISTman->FindOrBuildMaterial("G4_MUSCLE_STRIATED_ICRU");
  NISTman->FindOrBuildMaterial("G4_SKIN_ICRP");
  NISTman->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRP");
  NISTman->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRU-4");
  NISTman->FindOrBuildMaterial("G4_ADIPOSE_TISSUE_ICRP");
  NISTman->FindOrBuildMaterial("G4_A-150_TISSUE");
  NISTman->FindOrBuildMaterial("G4_MS20_TISSUE");
  NISTman->FindOrBuildMaterial("G4_MUSCLE_WITH_SUCROSE");
  NISTman->FindOrBuildMaterial("G4_MUSCLE_WITHOUT_SUCROSE");
  NISTman->FindOrBuildMaterial("G4_PLEXIGLASS");


  // mixture of stuff for misc. electronics etc.
  density = 2.*g/cm3;
  G4Material* MiscMaterial = 
    new G4Material(name="misc", density, ncomponents=5);
  MiscMaterial->AddElement(elAl, 20.*perCent);
  MiscMaterial->AddElement(elFe, 20.*perCent);
  MiscMaterial->AddElement(elCu, 20.*perCent);
  MiscMaterial->AddElement(elC,  20.*perCent);
  MiscMaterial->AddElement(elO,  20.*perCent);

  //aluminum 6061
  density = 2.70*g/cm3;
  G4Material* Aluminum6061 =
    new G4Material(name="Aluminum6061", density, ncomponents=5);
  //original definitions from xray phsp v2 and earlier
  //Aluminum6061->AddElement(elAl,97.7*perCent);
  //Aluminum6061->AddElement(elSi, 0.6*perCent);
  //Aluminum6061->AddElement(elMg, 1.0*perCent);
  //Aluminum6061->AddElement(elCu, 0.3*perCent);
  //Aluminum6061->AddElement(elCr, 0.2*perCent);
  //Aluminum6061->AddElement(elFe, 0.2*perCent);

  // newer, designed to reduce amount of highZ materials
  Aluminum6061->AddElement(elAl,98.01*perCent);
  Aluminum6061->AddElement(elSi, 0.6 *perCent);
  Aluminum6061->AddElement(elMg, 1.2 *perCent);
  Aluminum6061->AddElement(elCu, 0.15*perCent);
  Aluminum6061->AddElement(elCr, 0.04*perCent);

  //kapton
  density = 1.42*g/cm3;
  G4Material* kapton = new G4Material(name="kapton", density, ncomponents=4);
  kapton->AddElement(elH, 02.6362 *perCent);
  kapton->AddElement(elC, 69.1133 *perCent);
  kapton->AddElement(elN, 07.3270 *perCent);
  kapton->AddElement(elO, 20.9235 *perCent);

  //brass
  density = 8.53*g/cm3;
  G4Material* Brass = new G4Material(name="Brass", density, ncomponents=4);
  Brass->AddElement(elZn,29.88*perCent);
  Brass->AddElement(elPb, 0.07*perCent);
  Brass->AddElement(elFe, 0.05*perCent);
  Brass->AddElement(elCu,70.  *perCent);

  //  SS304  
  density = 8.0*g/cm3;
  G4Material* SS304 = new G4Material(name="SS304", density, ncomponents=9);
  SS304->AddElement(elS,  0.03*perCent);
  SS304->AddElement(elC,  0.08*perCent);
  SS304->AddElement(elSi, 0.75*perCent);
  SS304->AddElement(elP,  0.045*perCent);
  SS304->AddElement(elN,  0.10*perCent);
  SS304->AddElement(elCr,19.00*perCent);
  SS304->AddElement(elMn, 2.00*perCent);
  SS304->AddElement(elNi, 9.25*perCent);
  SS304->AddElement(elFe,68.745*perCent);

  // SS12L14
  density = 7.87*g/cm3;
  G4Material* SS12L14 = new G4Material(name="SS12L14", density, ncomponents=5);
  SS12L14->AddElement(elFe, 98.4*perCent);
  SS12L14->AddElement(elPb, 0.25*perCent);
  SS12L14->AddElement(elMn, 1.0 *perCent);
  SS12L14->AddElement(elP,  0.05*perCent);
  SS12L14->AddElement(elS,  0.3 *perCent);

  // SS G10080  low carbon steel
  density = 7.87*g/cm3;
  G4Material* SSlowcarbon = 
    new G4Material(name="SSlowcarbon", density, ncomponents=5);
  SSlowcarbon->AddElement(elFe, 99.4*perCent);
  SSlowcarbon->AddElement(elMn, 0.4 *perCent);
  SSlowcarbon->AddElement(elC,  0.1 *perCent);
  SSlowcarbon->AddElement(elS,  0.05*perCent);
  SSlowcarbon->AddElement(elP,  0.05*perCent);

  //aluminum 5052,  // from asm.matweb.com
  density = 2.68*g/cm3;
  G4Material* Aluminum5052 = 
    new G4Material(name="Aluminum5052", density, ncomponents=6);
  Aluminum5052->AddElement(elAl, 96.9 *perCent);
  Aluminum5052->AddElement(elMg,  2.5 *perCent);
  Aluminum5052->AddElement(elCr,  0.25*perCent);
  Aluminum5052->AddElement(elCu,  0.05*perCent);
  Aluminum5052->AddElement(elSi,  0.1 *perCent);
  Aluminum5052->AddElement(elFe,  0.2 *perCent);

  //aluminum 2024   //from wikipedia
  density = 2.78*g/cm3;
  G4Material* Aluminum2024 = 
    new G4Material(name="Aluminum2024", density, ncomponents=4);
  Aluminum2024->AddElement(elAl,93.7*perCent);
  Aluminum2024->AddElement(elMg, 1.4*perCent);
  Aluminum2024->AddElement(elMn, 0.5*perCent);
  Aluminum2024->AddElement(elCu, 4.4*perCent);

  // 95% density tungsten
  density = 18.0*g/cm3;
  G4Material* W95 = new G4Material(name="W95", density, ncomponents=3);
  W95->AddElement(elW , 95.0*perCent);
  W95->AddElement(elNi,  3.5*perCent);
  W95->AddElement(elFe,  1.5*perCent);
  //W95->AddElement(elCu,  1.0*perCent);

  // 92.5% density tungsten
  density = 17.6*g/cm3;
  G4Material* W92_5 = new G4Material(name="W92_5", density, ncomponents=3);
  W92_5->AddElement(elW , 92.5*perCent);
  W92_5->AddElement(elNi,  5.25*perCent);
  W92_5->AddElement(elFe,  2.25*perCent);
  //W95->AddElement(elCu,  1.0*perCent);

  // Lead-Antimony
  density = 11.0*g/cm3;
  G4Material* Lead97Antimony = 
    new G4Material(name="Lead97Antimony", density, ncomponents=2);
  Lead97Antimony->AddElement(elPb, 97.*perCent);
  Lead97Antimony->AddElement(elSb,  3.*perCent);

  // nicoro
  density = 10.9*g/cm3;
  G4Material* Nicoro = new G4Material(name="Nicoro", density, ncomponents=3);
  Nicoro->AddElement(elAu,35.*perCent);
  Nicoro->AddElement(elCu,62.*perCent);
  Nicoro->AddElement(elNi, 3.*perCent);

  // copper glidcop
  density = 8.88528*g/cm3;
  G4Material* copper_glidcop = 
    new G4Material(name="COPPER_GLIDCOP", density, ncomponents=3);
  copper_glidcop->AddElement(elCu, 99.600 *perCent);
  copper_glidcop->AddElement(elO,  0.1883 *perCent);
  copper_glidcop->AddElement(elAl, 0.2117 *perCent);

  // SS_A36
  density = 7.86112*g/cm3;
  G4Material* SS_A36 = new G4Material(name="SS_A36", density, ncomponents=5);
  SS_A36->AddElement(elC,   0.25*perCent);
  SS_A36->AddElement(elFe, 99.26*perCent);
  SS_A36->AddElement(elSi,  0.4 *perCent);
  SS_A36->AddElement(elP,   0.04*perCent);
  SS_A36->AddElement(elS,   0.05*perCent);

  // wafer alloy
  density = 11.0*g/cm3;
  G4Material* CuAuAlloy = 
    new G4Material(name="CuAuAlloy", density, ncomponents=2);
  CuAuAlloy->AddElement(elAu, 35.*perCent);
  CuAuAlloy->AddElement(elCu, 65.*perCent);

  // concrete
  //density = 2.03*g/cm3;
  //G4Material* Concrete = new G4Material("Concrete", density, 10);
  //Concrete->AddElement(elH ,  0.01);
  //Concrete->AddElement(elO ,  0.529);
  //Concrete->AddElement(elNa,  0.016);
  //Concrete->AddElement(elHg,  0.002);
  //Concrete->AddElement(elAl,  0.034);
  //Concrete->AddElement(elSi,  0.337);
  //Concrete->AddElement(elK ,  0.013);
  //Concrete->AddElement(elCa,  0.044);
  //Concrete->AddElement(elFe,  0.014);
  //Concrete->AddElement(elC ,  0.001);

  // Zinc ZA8
  density = 6.3*g/cm3;
  G4Material* ZincZA8 = new G4Material("ZincZA8", density, 8);
  ZincZA8->AddElement(elZn, 89.78*perCent);
  ZincZA8->AddElement(elAl,  8.8*perCent);
  ZincZA8->AddElement(elCu,  1.3*perCent);
  ZincZA8->AddElement(elFe,  0.075*perCent);
  ZincZA8->AddElement(elMn,  0.03*perCent);
  ZincZA8->AddElement(elPb,  0.006*perCent);
  ZincZA8->AddElement(elCd,  0.006*perCent);
  ZincZA8->AddElement(elSn,  0.003*perCent);

  // cerrotru   // from applicator drawing
  density = 6.3*g/cm3;  //  csalloys.com  Tru 281
  G4Material* cerrotru = new G4Material(name="cerrotru", density, 2);
  cerrotru->AddElement(elBi,58.*perCent);
  cerrotru->AddElement(elSn,42.*perCent);

  // PETG   PET is C_10 H_8 O_4; wikipedia says replace ethylene glycol with
  // cyclohexane dimethanal  (CH2OH)2 -> C6H10(CH2OH)2  i.e. add C6H10
  G4Material* PETG = new G4Material("PETG", 1.38*g/cm3, 3);
  PETG->AddElement(elC, 16);
  PETG->AddElement(elH, 18);
  PETG->AddElement(elO, 4);

  // sintered diamond, from Dave P.  // cf. my paper
  density = 3.6*g/cm3;  // 3.5
  G4Material* sinteredDiamond = new G4Material("sinteredDiamond", density, 3);
  sinteredDiamond->AddElement(elC,  88.5*perCent);  // 88.4
  sinteredDiamond->AddElement(elSi, 10.5*perCent);  // 9
  sinteredDiamond->AddElement(elTi,  1.0*perCent);  // 2.6

  // from Eclipse manual
  // Cartilage
  density = 1.1*g/cm3;
  G4Material* Cartilage = new G4Material("cartilage", density, 8);
  Cartilage->AddElement(elH ,  0.096);
  Cartilage->AddElement(elC ,  0.099);
  Cartilage->AddElement(elN ,  0.022);
  Cartilage->AddElement(elO ,  0.744);
  Cartilage->AddElement(elNa,  0.005);
  Cartilage->AddElement(elP,   0.022);
  Cartilage->AddElement(elS,   0.009);
  Cartilage->AddElement(elCl,  0.003);

  // Ti alloy, from Eclipse manual
  density = 4.42*g/cm3;
  G4Material* TiAlloy = new G4Material("TiAlloy", density, 3);
  TiAlloy->AddElement(elTi, 0.90);
  TiAlloy->AddElement(elAl, 0.06);
  TiAlloy->AddElement(elV,  0.04);

  // Stainless Steel, from Eclipse manual
  density = 8.00*g/cm3;
  G4Material* SSteel = new G4Material("StainlessSteel", density, 7);
  SSteel->AddElement(elC,  0.00080);
  SSteel->AddElement(elSi, 0.01);
  SSteel->AddElement(elP,  0.00045);
  SSteel->AddElement(elCr, 0.19);
  SSteel->AddElement(elMn, 0.02);
  SSteel->AddElement(elFe, 0.68375);
  SSteel->AddElement(elNi, 0.095);

  // Wood/cork
  density = 0.70*g/cm3;
  G4Material* Wood = new G4Material("wood", density, 3);
  Wood->AddElement(elH, 0.06216);
  Wood->AddElement(elC, 0.44445);
  Wood->AddElement(elO, 0.49339);
  
  // polystyrene
  density = 1.05*g/cm3;
  G4Material* Polystyrene = new G4Material("polystyrene", density, 2);
  Polystyrene->AddElement(elH, 0.07742);
  Polystyrene->AddElement(elC, 0.92258);

  // other plastics
  density = 1.04*g/cm3;  // density varies depending on material
  G4Material* Plastics = new G4Material("plastics", density, 4);
  Plastics->AddElement(elH, 0.05011);
  Plastics->AddElement(elC, 0.73282);
  Plastics->AddElement(elO, 0.14461);
  Plastics->AddElement(elS, 0.07246);

  // H2O
  density = 1.0*g/cm3;
  G4int nAtoms;
  G4Material* H2O = new G4Material("H2O", density, 2);
  H2O->AddElement(elH, nAtoms=2);
  H2O->AddElement(elO, nAtoms=1);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::AddMaterial(
    G4String name, G4double density, G4String basemat) {
  new G4Material(name, density, G4Material::GetMaterial(basemat));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildWorld() {
  //============================================================================
  // World Volume 
  //============================================================================

  G4ThreeVector worldSize = G4ThreeVector(1000.*cm, 1000.*cm, 1000.*cm);
  
  G4Box* solidWorld =
    new G4Box("world", worldSize.x()/2., worldSize.y()/2., worldSize.z()/2.);
  fLogicWorld =
    new G4LogicalVolume(solidWorld, G4Material::GetMaterial("G4_AIR"), 
      "World", 0, 0, 0);

  G4VisAttributes* world_VisAtt = 
    new G4VisAttributes(G4Colour(1., 1., 1., .0));
  world_VisAtt->SetForceWireframe(true);
  fLogicWorld->SetVisAttributes(world_VisAtt);  
  //fLogicWorld->SetVisAttributes(G4VisAttributes::Invisible);  

  fPhysiWorld =
    new G4PVPlacement(0, G4ThreeVector(), fLogicWorld, "World",         
                        0, false, 0);       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildWaterPhantom(G4ThreeVector phantom_position, G4ThreeVector phantom_box_size) {
    //============================================================================
    // Water Phantom Box
    //============================================================================

    const G4String phantom_name = "WaterPhantom";
    const double SAD = fSAD;
    G4Box* box = new G4Box(phantom_name, phantom_box_size.x() / 2.0, phantom_box_size.y() / 2.0, phantom_box_size.z() / 2.0);
    //G4LogicalVolume* logical = new G4LogicalVolume(box, G4Material::GetMaterial("G4_Pb"), phantom_name, 0, 0, 0);
    G4LogicalVolume* logical = new G4LogicalVolume(box, G4Material::GetMaterial("G4_WATER"),phantom_name, 0, 0, 0);
    G4VisAttributes* vis_attrib = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 1.0));
    vis_attrib->SetVisibility(true);
    logical->SetVisAttributes(vis_attrib);
    new G4PVPlacement(nullptr, phantom_position, logical, phantom_name, fLogicWorld, false, 0);
    auto region = new G4Region(phantom_name);
    logical->SetRegion(region);
    region->AddRootLogicalVolume(logical);
    G4cout << "WaterPhantom position " << phantom_position / mm << G4endl; 
    G4cout << "WaterPhantom size " << phantom_box_size / mm << G4endl;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildGantryMother() {
  //============================================================================
  // mother volume for everything that rotates with gantry 
  //============================================================================

  G4Box* gantryBox = new G4Box("gantry", 45.*cm, 45.*cm, 65.*cm);
  fGantry_LV = 
    new G4LogicalVolume(gantryBox, 
                        G4Material::GetMaterial("G4_AIR"), 
                        "gantryMother", 0, 0, 0);

  G4VisAttributes* gantry1_VisAtt = 
    new G4VisAttributes(G4Colour(0., 0.2, 0.8, 0.01));
  fGantry_LV->SetVisAttributes(gantry1_VisAtt);

  fGantry =
    new G4PVPlacement(0, G4ThreeVector(0., 0., fGantryPos),
      fGantry_LV, "GantryMother", fLogicWorld, false, 0);

  SetGantryRotation(fGantryRot); 
}
  
  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildPhotonGeometry(G4bool build) {
  ////==========================================================================
  // Xray beams  (only components specific to them)
  ////==========================================================================

  BuildTarget(build);
  BuildPrimaryCollimator(build);
  BuildYStageShield(build);
  BuildFlatteningFilter(build);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildElectronGeometry(G4bool build) {
  ////==========================================================================
  //// Electron beam components (only components specific to them)
  ////==========================================================================

  BuildFoil1(build);
  BuildFoil2(build);
  BuildApplicator();
  BuildCutOut();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildCollimatorMother() { }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildShieldingCollimator() { }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildCollimators(G4bool build) {
  (void)build;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildApplicator() { } 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildCutOut() { }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildIonChamber(G4double position) {
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
  G4VPhysicalVolume* myVol = store->GetVolume("ICMother", false);
  if (myVol) store->DeRegister(myVol);

  //TODO possibility of not using?
  //
  //  note that the top/bottom ring spacers are such that the amount of
  //  material is approx. correct. There are small extra bits of material
  //  that will change the spacing of the kapton windows, though, and
  //  these are not modelled. Therefore the kapton electrodes are in the
  //  wrong positions

  // centre of entire IC
  const G4double ic_pos             = fSAD - position; //156.9835*mm;
  // inner half thickness of ic body
  const G4double icbody_inthick     = 8.31*mm;
  // half-thickness of .005" kapton
  const G4double ICWin_5thou        = 0.0635*mm;
  // halfthickness of .002" kapton
  const G4double IC_electrodeHthick = 0.0254*mm;
  // half thickness of center holder
  const G4double IC_centralSupportHthick = 0.635*mm;
  // half thickness of center ring
  const G4double IC_centerringHthick = 0.508*mm;
  // half thickness of bulk of spacer
  const G4double IC_topringHthick    = 1.10*mm;
  const G4double IC_SE_pos = IC_centralSupportHthick + 2.*IC_topringHthick +
                             IC_electrodeHthick;
  const G4double IC_HV_pos = IC_centralSupportHthick + 2.*IC_electrodeHthick +
                             2.*IC_topringHthick + 2.*IC_centerringHthick +
                             IC_electrodeHthick;

  G4VisAttributes *VisAtt_IonChamber1 =
            new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.5));
  G4VisAttributes *VisAtt_IonChamber2 =
            new G4VisAttributes(G4Colour(0.0, 0.2, 1.0, 0.8));
  G4VisAttributes *VisAtt_IonChamber3 =
            new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.6));
  G4VisAttributes *VisAtt_IonChamber4 =
            new G4VisAttributes(G4Colour(0.2, 1.0, 0.0, 0.4));
  //G4VisAttributes *VisAtt_IonChamber5 =
  //          new G4VisAttributes(G4Colour(0.0, 0.5, 0.5, 0.4));

  // first create a mother volume to put the bits in
  // the LV needs to be a class variable because sensitive detectors
  // created in diff. function
  // TODO: air in ion chamber is at 1.2 bar
  G4Tubs* ICMotherTubs = new G4Tubs("ICMother", 0., 80.*mm, 16.0*mm,
            0.*deg,360.*deg);
  fICMother_LV = new G4LogicalVolume(ICMotherTubs,
            G4Material::GetMaterial("G4_AIR"), "ICMother_LV", 0, 0, 0);

  G4VisAttributes* ICMotherLV_VisAtt =
            new G4VisAttributes(G4Colour(0.2, 0.4, 0.7, 0.8));
  ICMotherLV_VisAtt->SetVisibility(false);
  fICMother_LV->SetVisAttributes(ICMotherLV_VisAtt);

  G4RotationMatrix* ICRot = new G4RotationMatrix();
  ICRot->rotateX(180.0*deg);
  
  fICMother = new G4PVPlacement(ICRot, 
                    G4ThreeVector(0.*m, 0.*m, ic_pos - fGantryPos), 
                    fICMother_LV,
                    "ICMother", fGantry_LV, false, 0);

  //============================== Ion Chamber body =========================
  const G4int icbody_pts = 6;
  G4double rInner_icbody[icbody_pts] =
            {50.16*mm, 49.02*mm, 66.68*mm, 66.68*mm, 49.02*mm, 50.16*mm};
  G4double rOuter_icbody[icbody_pts] =
            {76.2 *mm, 76.2 *mm, 76.2 *mm, 76.2 *mm, 76.2 *mm, 76.2 *mm};
  G4double zPlane_icbody[icbody_pts] =
            {-15.09*mm, -icbody_inthick, -icbody_inthick,
             icbody_inthick,  icbody_inthick, 15.09*mm};

  G4Polycone* IC_body_PV = new G4Polycone("IC_body", 0.*deg,
            360.*deg, icbody_pts, zPlane_icbody, rInner_icbody, rOuter_icbody);

  G4LogicalVolume* IC_body_LV =
            new G4LogicalVolume(IC_body_PV,G4Material::GetMaterial("SS304"),
            "ic_body_LV", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(),
            IC_body_LV, "icbody", fICMother_LV, false, 0);

  IC_body_LV->SetVisAttributes(VisAtt_IonChamber1);

  //============================== Ion Chamber Window - 3 instances===========
  //  Cu plating
  //  put it on downstream side for all, for simplicity
  fICWinCu_tubs = new G4Tubs("ICWinCu", 0.*mm, 48.*mm,
            fICCuHThick*fCopperThicknessFactor, 0.*deg, 360.*deg);
  fICWinCu_LV = new G4LogicalVolume(fICWinCu_tubs,
            G4Material::GetMaterial("G4_Cu"),"ICWinCU_LV", 0, 0, 0);
  fICWinCu_LV->SetVisAttributes(VisAtt_IonChamber3);

  G4Tubs* ICWin = new G4Tubs("ICWin", 0.*mm, 54.61*mm,
            ICWin_5thou, 0*deg, 360*deg);
  G4LogicalVolume* ICWin_LV = new G4LogicalVolume(ICWin,
            G4Material::GetMaterial("kapton"), "ICWindow_LV", 0, 0, 0);
  ICWin_LV->SetVisAttributes(VisAtt_IonChamber2);

  //first window
  new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, -icbody_inthick + ICWin_5thou),
            ICWin_LV,  "IonChamberWindow-1", fICMother_LV, false, 0);

  fICWinCu0 = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, -icbody_inthick +
            2.*ICWin_5thou + fICCuHThick * fCopperThicknessFactor),
            fICWinCu_LV, "ICWinCu", fICMother_LV, false, 0);

  //second window
  new G4PVPlacement(0, G4ThreeVector(),
            ICWin_LV, "IonChamberWindow-2", fICMother_LV, false, 1);

  fICWinCu1 = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, ICWin_5thou +
            fICCuHThick * fCopperThicknessFactor),
            fICWinCu_LV, "ICWinCu", fICMother_LV, false, 1);

  //third window
  new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, icbody_inthick - ICWin_5thou),
            ICWin_LV, "IonChamberWindow-3", fICMother_LV, false, 2);

  fICWinCu2 = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, icbody_inthick
            - 2.*ICWin_5thou - fICCuHThick * fCopperThicknessFactor),
            fICWinCu_LV, "ICWinCu", fICMother_LV, false, 2);


  //============================================= electrodes =================
  //============================================= spacer rings ===============

  // top and bottom rings (Al 2024)
  G4Tubs* IC_topring = new G4Tubs("IC_topring", 50.42*mm, 58.42*mm,
            IC_topringHthick, 0.*deg, 360.*deg);
  G4LogicalVolume* IC_topring_LV = new G4LogicalVolume(IC_topring,
            G4Material::GetMaterial("Aluminum2024"),"IC_topring_LV",0,0,0);

  IC_topring_LV->SetVisAttributes(VisAtt_IonChamber4);

  new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m,-IC_topringHthick -
            IC_centralSupportHthick),
            IC_topring_LV, "topOrBottomSpacer", fICMother_LV, false, 0);

  new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, +IC_topringHthick +
            IC_centralSupportHthick),
            IC_topring_LV, "topOrBottomSpacer", fICMother_LV, false, 1);

  new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, -IC_HV_pos -
            IC_electrodeHthick-IC_topringHthick),
            IC_topring_LV, "topOrBottomSpacer", fICMother_LV, false, 2);

  new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, +IC_HV_pos +
            IC_electrodeHthick+IC_topringHthick),
            IC_topring_LV, "topOrBottomSpacer", fICMother_LV, false, 3);

  // center ring
  G4Tubs* IC_centerring = new G4Tubs("IC_topring", 50.42*mm, 58.42*mm,
            IC_centerringHthick, 0.*deg, 360.*deg);
  G4LogicalVolume* IC_centerring_LV = new G4LogicalVolume(IC_centerring,
            G4Material::GetMaterial("G4_Cu"),"IC_centerring_LV",0,0,0);

  IC_centerring_LV->SetVisAttributes(VisAtt_IonChamber4);

  new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, +(IC_centralSupportHthick +
            2.*IC_electrodeHthick + 2.*IC_topringHthick +
            IC_centerringHthick)),
            IC_centerring_LV, "centerSpacer", fICMother_LV, false, 0);

  new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, -(IC_centralSupportHthick +
            2.*IC_electrodeHthick + 2.*IC_topringHthick +
            IC_centerringHthick)),
            IC_centerring_LV, "centerSpacer", fICMother_LV, false, 1);

  //============================================= HV Electrode ===============

  // kapton
  G4Tubs* IC_HV = new G4Tubs("IC_HV", 0.*mm, 54.61*mm,
            IC_electrodeHthick, 0.*deg, 360.*deg);
  G4LogicalVolume* IC_HV_LV = new G4LogicalVolume(IC_HV,
            G4Material::GetMaterial("kapton"), "IC_HV_LV", 0, 0, 0);

  IC_HV_LV->SetVisAttributes(VisAtt_IonChamber2);

  //  Cu plating
  fICElCu_tubs = new G4Tubs("ICElCu", 0.*mm, 48.0*mm,
            fICCuHThick*fCopperThicknessFactor, 0.*deg, 360.*deg);
  fICElCu_LV = new G4LogicalVolume(fICElCu_tubs,
            G4Material::GetMaterial("G4_Cu"), "ICElCU_LV", 0, 0, 0);
  fICElCu_LV->SetVisAttributes(VisAtt_IonChamber3);

  new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, -IC_HV_pos),
            IC_HV_LV, "IC_ElectrodeHV-1", fICMother_LV, false, 0);

  fICElCu0 = new G4PVPlacement(0,
            G4ThreeVector(0.*m, 0.*m, -IC_HV_pos+IC_electrodeHthick +
            fICCuHThick*fCopperThicknessFactor),
            fICElCu_LV, "ICElCu", fICMother_LV, false, 0);

  new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, IC_HV_pos),
            IC_HV_LV, "IC_ElectrodeHV-2", fICMother_LV, false, 1);

  fICElCu1 = new G4PVPlacement(0,
            G4ThreeVector(0.*m, 0.*m, IC_HV_pos + IC_electrodeHthick +
            fICCuHThick*fCopperThicknessFactor),
            fICElCu_LV, "ICElCu",fICMother_LV, false, 1);
  
  //============================================= Signal Electrode ===========

  //  Cu plating, ignoring the patterning
  fICElCusig_tubs = new G4Tubs("ICElCu", 0.*mm, 48.*mm,
            fICCuHThick*fCopperThicknessFactor, 0.*deg, 360.*deg);
  fICElCusig_LV = new G4LogicalVolume(fICElCusig_tubs,
            G4Material::GetMaterial("G4_Cu"), "ICElCU_sig1_LV", 0, 0, 0);
  fICElCusig_LV->SetVisAttributes(VisAtt_IonChamber3);

  // reuse the HV electrode LV for the kapton
  new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, -IC_SE_pos),
            IC_HV_LV, "IC_SignalElectrode-1", fICMother_LV, false, 0);

  fICElCusig0 = new G4PVPlacement(0,G4ThreeVector(0.*m, 0.*m, -IC_SE_pos +
            IC_electrodeHthick + fICCuHThick*fCopperThicknessFactor),
            fICElCusig_LV, "ICElCu_sig1", fICMother_LV, false, 0);

  new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, IC_SE_pos),
            IC_HV_LV, "IC_SignalElectrode-2", fICMother_LV, false, 1);

  fICElCusig1 = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, IC_SE_pos +
            IC_electrodeHthick + fICCuHThick*fCopperThicknessFactor),
            fICElCusig_LV, "ICElCu_sig1", fICMother_LV, false, 1);

  //=============================== ring for central kapton ==================

  const G4int IC_WinRing_pts = 4;
  G4double rInner_icr[IC_WinRing_pts] =
                              {55.37 *mm, 55.37 *mm, 48.89 *mm, 48.89 *mm};
  //rOuter lower than nominal; extra is in ic_body
  G4double rOuter_icr[IC_WinRing_pts] =
                              {66.67 *mm, 66.67 *mm, 66.67 *mm, 66.67 *mm};
  G4double zPlane_icr[IC_WinRing_pts] =
                              { -.635*mm,   .076*mm,   .076*mm,   .635*mm};

  G4Polycone* IC_WinRing_PV = new G4Polycone("IC_windowring", 0.*deg,
            360.*deg, IC_WinRing_pts, zPlane_icr, rInner_icr, rOuter_icr);

  G4LogicalVolume* IC_WinRing_LV =
          new G4LogicalVolume(IC_WinRing_PV,G4Material::GetMaterial("G4_Cu"),
          "IC_windowring_LV",0,0,0);

  IC_WinRing_LV->SetVisAttributes(VisAtt_IonChamber4);

  new G4PVPlacement(0, G4ThreeVector(),
            IC_WinRing_LV, "ic_winring", fICMother_LV, false, 0);

  //=============================== volume for sensitive detector ==============
 
  G4Tubs* dose_tubs = 
    new G4Tubs("IC_dose", 0.*mm, 48.*mm, 0.4*mm, 0.*deg, 360.*deg);

  fIC_dose_LV = 
    new G4LogicalVolume(dose_tubs, G4Material::GetMaterial("G4_AIR"),
    "ic_dose_lv", 0, 0, 0);

  fIC_dose_LV->SetVisAttributes(VisAtt_IonChamber2);

  fIC_dose = 
    new G4PVPlacement(0, 
      G4ThreeVector(0., 0., IC_SE_pos + (IC_HV_pos - IC_SE_pos)/2), 
      fIC_dose_LV, "ic_dose", fICMother_LV, false, 0);
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildShielding(G4bool build) {
  (void)build;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildVault(G4bool build) {
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
void TB02_BaseDetectorConstruction::BuildBackscatterShield(G4bool build) {
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
void TB02_BaseDetectorConstruction::BuildTarget(G4bool build) {
  (void)build;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildPrimaryCollimator(G4bool build) {
  (void)build;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildYStageShield(G4bool build) {
  (void)build;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildFlatteningFilter(G4bool build) {
  (void)build;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildFoil1(G4bool build) {
  (void)build;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildFoil2(G4bool build) {
  (void)build;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// custom foil2, layer 2 polycone
void TB02_BaseDetectorConstruction::AddCustomFoil2Vertex(
                                                      G4double r, G4double z) {
  (void)r;
  (void)z;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::ConstructSDandField() {
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
void TB02_BaseDetectorConstruction::BuildBeWindow() { }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetBuildBeWindow(G4bool build) {
  if (build == fBuildBeWindow) return;
  fBuildBeWindow = build;
  G4cout << "Build exit window set to: " << fBuildBeWindow << "." << G4endl;
  BuildBeWindow();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetProximalPosition(G4int, G4double) {}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetDistalPosition(G4int, G4double) {}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetMLCDensity(G4double) {}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetMLCNDS120Position(G4int index, G4double pos) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetMLCHDLeafPosition(G4int index, 
                                                       G4double pos) {
  // pos == position is positive for non-overtravel, for both leaf banks
  if (index < 0 || index > 119) {
    // TODO make this range specific to actual number of leaves
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
  else            newpos =  actual + fMlc_x_size/2.;

  if (newpos == orig.x())  return;  // leaf already in requested position

  leaf->SetTranslation(G4ThreeVector(newpos, orig.y(), orig.z()));

  //check to see if there is a collision
  //G4cout << "Setting MLC position: index " << index << " position: " 
  //       << pos/cm << G4endl;;
  if (opp_leaf) {  // may not have been build yet
    const G4ThreeVector opp_orig = opp_leaf->GetObjectTranslation();
    //G4cout << " newpos: " << newpos/cm
    //       << " actual: " << actual/cm
    //       << " opposite pos: " << opp_orig.x()/cm << G4endl
    //       << " newpos-opppos: " << (newpos - opp_orig.x())/cm << G4endl;
    //G4cout << "check:    " << GetMLCLeafPosition(index)/cm << " cm " ;
    //G4cout << "Opp.pos.: " << GetMLCLeafPosition(opp_index)/cm << " cm." 
    //       << G4endl;
    G4double opp_orig_pos = opp_orig.x() + fMlc_x_size/2.;
    if (index < 60) opp_orig_pos *= -1.;

    G4double leaf_delta = 0.01*mm;
    if (index < 60) {
      if (opp_orig.x() - newpos < fMlc_x_size) {
        //       << " newpos: " << newpos/cm
        //       << " opposite pos: " << opp_orig.x()/cm << G4endl
        //       << " newpos-opppos: " << (newpos - opp_orig.x())/cm << G4endl;
        //G4cout << "Setting opposite leaf " << opp_index << " to " 
        //       << -pos + leaf_delta << G4endl;
        SetMLCHDLeafPosition(opp_index, -pos + leaf_delta);
      }
    }
    else if (index >= 60) {
      if (newpos - opp_orig.x() < fMlc_x_size) {
        //       << " newpos: " << newpos/cm
        //       << " opposite pos: " << opp_orig.x()/cm << G4endl
        //       << " newpos-opppos: " << (newpos - opp_orig.x())/cm << G4endl;
        //G4cout << "Opp.pos.: " << GetMLCLeafPosition(opp_index)/cm << " cm." 
        //       << G4endl;
        //G4cout << "Setting opposite leaf " << opp_index << " to " 
        //       << -pos + leaf_delta << G4endl;
        SetMLCHDLeafPosition(opp_index, -pos + leaf_delta);
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
void TB02_BaseDetectorConstruction::SetJawPositionY(G4int jawint,
                                                G4double field, 
                                                G4bool force) {
  (void)jawint;
  (void)field;
  (void)force;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetJawPositionX(G4int jawint,
                                                G4double field,
                                                G4bool force) {
  (void)jawint;
  (void)field;
  (void)force;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetCollimatorRotation(G4double v) {
  if (v == fCollRot) return;
  fCollRot = v;

  G4RotationMatrix* collMotherRot = new G4RotationMatrix();
  collMotherRot->rotateZ(fCollRot);
  fColl->SetRotation(collMotherRot);
  
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4cout << "Collimator rotation set to: " << fCollRot/deg << " deg." 
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetGantryRotation(G4double v,
                                                  G4double theta) {
  // theta is the longitude (not physically possible with TrueBeam
  if (v == fGantryRot && theta == fGantryTheta) return;
  fGantryRot = v;
  fGantryTheta = theta;

  G4RotationMatrix* gantryRot = new G4RotationMatrix();
  gantryRot->rotateZ(-fGantryTheta);
  gantryRot->rotateY(-fGantryRot); // note minus sign! Due to volume/frame
  fGantry->SetRotation(gantryRot);
  G4ThreeVector defaultpos(0.,0.,fGantryPos);
  G4RotationMatrix* invRot = new G4RotationMatrix();
  invRot->rotateY(fGantryRot);  // this is the opposite of the rotation matrix
  invRot->rotateZ(fGantryTheta);
  G4ThreeVector newpos = *invRot * defaultpos;
  fGantry->SetTranslation(newpos);
  
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  //G4RunManager::GetRunManager()->ReOptimizeMotherOf(fGantry);
  G4cout << "Gantry rotation set to: " << fGantryRot/deg << " deg." 
         << G4endl;
  if (fGantryTheta != 0.) {
    G4cout << "Gantry theta rotation set to: " << fGantryTheta/deg << " deg." 
           << G4endl;
  }
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetSimulateCollimators(G4bool v) {
  if (v == fSimulateCollimators) return;
  fSimulateCollimators = v;
  BuildCollimators(fSimulateCollimators);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4cout << "Simulate collimators: " << fSimulateCollimators << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetSimulateVault(G4bool v) {
  if (v == fSimulateVault) return;
  fSimulateVault = v;
  // put the logic whether to build/destroy into BuildVault
  BuildVault(fSimulateVault);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4cout << "Simulate vault: " << fSimulateVault << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetSimulateShielding(G4bool v) {
  if (v == fSimulateShielding) return;
  fSimulateShielding = v;
  // put the logic whether to build/destroy into BuildShielding
  BuildShielding();
  BuildBackscatterShield(!fSimulateShielding);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4cout << "Simulate shielding: " << fSimulateShielding << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetJawOffset(G4double v) {
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
void TB02_BaseDetectorConstruction::SetUseIonChamber(G4bool b) {
  (void)b;
  //TODO make this an option?
  G4ExceptionDescription ed;
  ed << "Use ion chamber command not enabled.";
  G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac110",
              FatalException, ed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetExitWindowThicknessFactor(G4double v) {
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
void TB02_BaseDetectorConstruction::SetFoil2ThicknessFactor(G4double v) {
  G4cout << "Asked for: " << v << G4endl; // to avoid compiler warning

  G4ExceptionDescription ed;
  ed << "Foil2 thickness factor not in use.";
  G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac039",
              FatalException, ed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetCopperThicknessFactor(G4double v) {
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
//void TB02_BaseDetectorConstruction::SetKaptonThicknessFactor(G4double v) {
//  G4cerr << "Kapton thickness factor deprecated." << G4endl;
//  exit(EXIT_FAILURE);
//  //if (v == fKaptonThicknessFactor) return;
//  //fKaptonThicknessFactor = v;
//  //BuildIonChamber();
//  //G4cout << "I.C. kapton thickness factor: " << v << G4endl;
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetFoil1Name(G4String v) {
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
void TB02_BaseDetectorConstruction::SetFoil1ThicknessFactor(G4double v) {
  // TODO limits
  if (v == fFoil1ThicknessFactor) return;
  fFoil1ThicknessFactor = v;
  if (fBeamType == "electron") BuildFoil1();
  G4cout << "Foil1 thickness factor: " << v << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetFoil1Material(G4String v) {
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
void TB02_BaseDetectorConstruction::SetFoil1Thickness(G4double v) {
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
void TB02_BaseDetectorConstruction::SetFoil1Radius(G4double v) {
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
void TB02_BaseDetectorConstruction::SetFoil1Position(G4double v) {
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
void TB02_BaseDetectorConstruction::SetFoil2Name(G4String v) {
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
void TB02_BaseDetectorConstruction::SetFoil2Material(G4String v) {
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
void TB02_BaseDetectorConstruction::SetFoil2Position(G4double v) {
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
void TB02_BaseDetectorConstruction::SetTargetName(G4String v) {
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
void TB02_BaseDetectorConstruction::SetTargetPosition(G4double v) {
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
void TB02_BaseDetectorConstruction::SetTargetRadius(G4double v) {
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
void TB02_BaseDetectorConstruction::SetTargetThickness(G4double v) {
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
void TB02_BaseDetectorConstruction::SetTargetMaterial(G4String v) {
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
void TB02_BaseDetectorConstruction::SetTargetMaterial2(G4String v) {
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
void TB02_BaseDetectorConstruction::AddCustomTargetVertex(G4double r, 
                                                          G4double z) {
  // TODO verify inputs
  G4cout << "Adding custom target vertex r: " << r/mm << " z: " << z/mm
         << G4endl;
  fTargetCustom2_r.push_back(r);
  fTargetCustom2_z.push_back(z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetFFName(G4String v) {
  if (v == fFlatteningFilterName) return;
  fFlatteningFilterName = v;
  if (fBeamType == "xray") BuildFlatteningFilter();
  if (fVerbosity>0) {
    G4cout << "Flattening filter set to: " << fFlatteningFilterName << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetFFOffset(G4ThreeVector v) {
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
G4double TB02_BaseDetectorConstruction::GetJawPositionX(G4int jaw) {
  if (jaw==1) return fFieldX1;
  if (jaw==2) return fFieldX2;
  else return -9999.*cm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double TB02_BaseDetectorConstruction::GetJawPositionY(G4int jaw) {
  if (jaw==1) return fFieldY1;
  if (jaw==2) return fFieldY2;
  else return -9999.*cm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TB02_BaseDetectorConstruction::SetOutputFilename(G4String v) {
  if (v == fOutputFilename) return;
  fOutputFilename = v;
  //TB02_PhspWorldConstruction* phspWorld =
  //      (TB02_PhspWorldConstruction*)GetParallelWorld(1);
  //phspWorld->SetOutputFilename(fOutputFilename);

  if (fVerbosity > 0) {
    G4cout << "Output base file name: " << fOutputFilename << G4endl;
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetBeamType(G4String v) {
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
void TB02_BaseDetectorConstruction::SetApplicatorName(G4String v) {
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
void TB02_BaseDetectorConstruction::SetBuildCutOut(G4bool v) {
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
void TB02_BaseDetectorConstruction::AddCutOutVertex(G4double x, G4double y) {
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
void TB02_BaseDetectorConstruction::SetCutOutThickness(G4double v) {
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
void TB02_BaseDetectorConstruction::SetCutOutBevelFactor(G4double v) {
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
void TB02_BaseDetectorConstruction::SetCutOutMaterial(G4String v) {
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
void TB02_BaseDetectorConstruction::SetMLC(G4String) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::RemoveMLC() {

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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double TB02_BaseDetectorConstruction::GetMLCLeafPosition(G4int i) {
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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildCone() {
  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
  G4VPhysicalVolume* pv = store->GetVolume("cone", false);
  if (pv) store->DeRegister(pv);
  G4LogicalVolumeStore* lv_store = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* lv = lv_store->GetVolume("cone_LV", false);
  if (lv) lv_store->DeRegister(lv);

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4double delta = 0.001*mm;  // to avoid overlapping surfaces

  G4Tubs* cone1 = new G4Tubs("cone1", 0., fCone_OuterRadius, fCone_thick/2.,
                             0., twopi);

  // specify "upside down"
  G4Cons* cone2 = new G4Cons("cone2",
             0.*mm, ((fSAD - fCone_Zpos + fCone_thick)/fSAD) * fCone_size/2.,
             0.*mm, ((fSAD - fCone_Zpos)/fSAD) * fCone_size/2.,
             (fCone_thick + delta)/2.,
             0., twopi);

  G4SubtractionSolid* cone3 =
    new G4SubtractionSolid("cone3", cone1, cone2, 0, G4ThreeVector());

  G4LogicalVolume* cone_LV = new G4LogicalVolume(cone3,
          G4Material::GetMaterial("G4_W"), "cone_LV", 0, 0, 0);

  G4VisAttributes* VisAtt_cone =
              new G4VisAttributes(G4Colour(0.0, 0.5, 0.9, 0.6));
  cone_LV->SetVisAttributes(VisAtt_cone);

  //G4PVPlacement* cone = new G4PVPlacement(0,
  new G4PVPlacement(0,
      G4ThreeVector(0., 0., fCone_Zpos - fCollPos - fCone_thick/2.),
      cone_LV, "cone", fColl_LV, false, 0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetConeSize(G4double v) {
  if (v == fCone_size) return;
  if (v <= 0. || v > 50.*mm) {
    G4ExceptionDescription ed;
    ed << "Cone size " << v/mm << " mm out of range.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac500",
                FatalException, ed);
  }
  fCone_size = v;
  G4cout << "Cone size set to " << fCone_size/mm << " mm." << G4endl;
  BuildCone();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetConeThickness(G4double v) {
  if (v == fCone_thick) return;
  if (v <= 0. || v > 200.*mm) {
    G4ExceptionDescription ed;
    ed << "Cone thickness " << v/mm << " mm out of range.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac501",
                FatalException, ed);
  }
  fCone_thick = v;
  G4cout << "Cone thickness set to " << fCone_thick/mm << " mm." << G4endl;
  BuildCone();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetConeElevation(G4double v) {
  if (v == fCone_Zpos) return;
  if (v <= 20.*cm || v > 50.*cm) {
    G4ExceptionDescription ed;
    ed << "Cone elevation " << v/mm << " mm out of range.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac502",
                FatalException, ed);
  }
  fCone_Zpos = v;
  G4cout << "Cone elevation set to " << fCone_Zpos/mm << " mm." << G4endl;
  BuildCone();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetConeOuterRadius(G4double v) {
  if (v == fCone_OuterRadius) return;
  if (v <= 0. || v > 100.*mm) {
    G4ExceptionDescription ed;
    ed << "Cone outer radius " << v/mm << " mm out of range.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac503",
                FatalException, ed);
  }
  fCone_OuterRadius = v;
  G4cout << "Cone outer radius set to " << fCone_OuterRadius/mm << " mm."
         << G4endl;
  BuildCone();
}

// scatterer for TSI
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetTSIPlateUse(G4bool v) {
  if (v == fTSIPlateUse) return;
  fTSIPlateUse = v;
  if (v) G4cout << "Will build TSI scatterer." << G4endl;
  else   G4cout << "Will not build TSI scatterer." << G4endl;
  // this next line will remove if fTSIPlateUse = false
  BuildTSIPlate();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetTSIPlatePosition(G4double v) {
  if (v == fTSIPlatePosition) return;
  if (v < 0.*cm || v >= 50.*cm) {
    G4ExceptionDescription ed;
    ed << "TSI scatterer position " << v/cm << " cm out of range (0--50 cm).";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac504",
                FatalException, ed);
  }
  fTSIPlatePosition = v;
  G4cout << "TSI scatterer position set to " << fTSIPlatePosition/cm
         << " cm. " << G4endl;
  if (fTSIPlateUse) BuildTSIPlate();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetTSIPlateSideLength(G4double v) {
  if (v == fTSIPlateSideLength) return;
  if (v <= 0.*cm || v >= 50.*cm) {
    G4ExceptionDescription ed;
    ed << "TSI scatterer side length " << v/cm 
       << " cm out of range (0--50 cm).";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac505",
                FatalException, ed);
  }
  fTSIPlateSideLength = v;
  G4cout << "TSI scatterer side length set to " << fTSIPlateSideLength/cm
         << " cm. " << G4endl;
  if (fTSIPlateUse) BuildTSIPlate();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetTSIPlateThickness(G4double v) {
  if (v == fTSIPlateThickness) return;
  if (v <= 0.*cm || v >= 20.*mm) {
    G4ExceptionDescription ed;
    ed << "TSI scatterer thickness " << v/mm 
       << " mm out of range (0--20 mm).";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac506",
                FatalException, ed);
  }
  fTSIPlateThickness = v;
  G4cout << "TSI scatterer thickness set to " << fTSIPlateThickness/mm
         << " mm. " << G4endl;
  if (fTSIPlateUse) BuildTSIPlate();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::SetTSIPlateMaterial(G4String v) {
  if (v == fTSIPlateMaterial) return;
  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(v);
  
  if (pttoMaterial) {
    fTSIPlateMaterial = v;
  } else {
    G4ExceptionDescription ed;
    ed << "TSI scatterer material " << v << " not allowed.";
    G4Exception("VirtuaLinac:DetectorConstruction", "VirtuaLinac507",
                FatalException, ed);
  }
  G4cout << "TSI scatterer material set to " << fTSIPlateMaterial 
         << G4endl;
  if (fTSIPlateUse) BuildTSIPlate();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TB02_BaseDetectorConstruction::BuildTSIPlate() {
  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
  G4VPhysicalVolume* pv = store->GetVolume("tsiPlate", false);
  if (pv) store->DeRegister(pv);
  G4LogicalVolumeStore* lv_store = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* lv = lv_store->GetVolume("tsiPlate_LV", false);
  if (lv) lv_store->DeRegister(lv);

  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  if (!fTSIPlateUse) return; 
  G4Box* tsi_box1 = new G4Box("box1", fTSIPlateSideLength/2., fTSIPlateSideLength/2.,
                             fTSIPlateThickness/2.);

  G4LogicalVolume* tsi_LV = new G4LogicalVolume(tsi_box1,
          G4Material::GetMaterial(fTSIPlateMaterial), "tsiPlate_LV", 0, 0, 0);

  G4VisAttributes* VisAtt_tsi =
              new G4VisAttributes(G4Colour(0.3, 0.5, 0.9, 0.6));
  tsi_LV->SetVisAttributes(VisAtt_tsi);

  new G4PVPlacement(0,
      G4ThreeVector(0., 0., fTSIPlatePosition - fCollPos - fTSIPlateThickness/2.),
      tsi_LV, "tsiPlate", fColl_LV, false, 0);
}
