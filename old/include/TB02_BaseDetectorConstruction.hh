/* This file is the CONFIDENTIAL and PROPRIETARY information of
 * Varian Medical Systems, Inc., and is not to be distributed.
 *
 * Copyright (c) 2017 Varian Medical Systems, Inc.
 *
 * For information, contact Daren Sawkey  daren.sawkey@varian.com
 */

#ifndef TB02_BaseDetectorConstruction_h
#define TB02_BaseDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4Cache.hh"

#include "G4ExtrudedSolid.hh"  // provides G4TwoVector

#include <unordered_map>
#include "../linac.h"

class G4Box;
class G4Sphere;
class G4Tubs;
class G4GenericPolycone;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4PVPlacement;
class G4SDManager;
class G4VIStore;
class G4PVReplica;
class G4Region;

std::unordered_map<int, TraversedGeometry> generateVolumeToTraversed(const VolumeNameAndTraversalFlag* _array, const size_t _array_size, G4LogicalVolume* _root_volume);

//class TB02_PhaseSpaceWriter;
class G4PSDoseDeposit;

class TB02_BaseDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  // constructor and destructor.
  TB02_BaseDetectorConstruction();
  virtual ~TB02_BaseDetectorConstruction();

  // virtual method from G4VUserDetectorConstruction.
  virtual G4VPhysicalVolume* Construct() = 0;
  virtual void GenerateTraversedTable() = 0;
  virtual void DefineMaterials();
  virtual void ConstructSDandField();
  virtual void BuildWorld();
  virtual void BuildWaterPhantom(G4ThreeVector phantom_position, G4ThreeVector phantom_box_size);
  virtual void BuildGantryMother();  // mother volume for gantry,tx head etc.

  virtual void BuildBeWindow();

  virtual void BuildElectronGeometry(G4bool build = true);
  virtual void BuildPhotonGeometry(G4bool build = true);

  virtual void BuildBackscatterShield(G4bool build = true);
  virtual void BuildTarget(G4bool build = true);
  virtual void BuildPrimaryCollimator(G4bool build = true);
  virtual void BuildYStageShield(G4bool build = true);
  virtual void BuildFlatteningFilter(G4bool build = true);

  virtual void BuildFoil1(G4bool build = true);
  virtual void BuildFoil2(G4bool build = true);

  virtual void BuildShielding(G4bool build = true);
  virtual void BuildIonChamber(G4double position);
  virtual void BuildShieldingCollimator();
  virtual void BuildCollimatorMother(); //mother volume for shielding coll 
                                // and collim (jaw+mlc etc)
  virtual void BuildCollimators(G4bool build = true);  // builds jaws and tx head components after:
                            // incl base plate, MLC, mylar window
  virtual void BuildApplicator();
  virtual void BuildCutOut();
  //void BuildPhaseSpace();
  //void RemovePhaseSpace(G4int index);
  virtual void BuildVault(G4bool build = true);


  virtual void     SetBeamType(G4String v);
  G4String GetBeamType(void) const {return fBeamType;}

  virtual void   SetBuildBeWindow(G4bool build = true);
  G4bool GetBuildBeWindow(void) const {return fBuildBeWindow;}

  virtual void   SetUseIonChamber(G4bool v);
  G4bool GetUseIonChamber(void) const {return fUseIonChamber;}

  virtual void     SetExitWindowThicknessFactor(G4double v);
  G4double GetExitWindowThicknessFactor(void) const 
    {return fExitWindowThicknessFactor;}

  virtual void     SetFoil1Name(G4String s);
  G4String GetFoil1Name(void) const {return fFoil1Name;}

  virtual void     SetFoil1Material(G4String v);
  G4String GetFoil1Material(void) const {return fFoil1Material;}

  virtual void     SetFoil1Thickness(G4double v);
  G4double GetFoil1Thickness(void) const {return fFoil1Thickness;}

  virtual void     SetFoil1Radius(G4double v);
  G4double GetFoil1Radius(void) const {return fFoil1Radius;}

  virtual void     SetFoil1Position(G4double v);
  G4double GetFoil1Position(void) const {return fFoil1Position;}

  virtual void     SetFoil1ThicknessFactor(G4double v);
  G4double GetFoil1ThicknessFactor(void) {return fFoil1ThicknessFactor;}
 
  virtual void     SetFoil2Name(G4String s);
  G4String GetFoil2Name(void) {return fFoil2Name;}

  virtual void     SetFoil2Material(G4String v);
  G4String GetFoil2Material(void) const {return fFoil2Material;}

  virtual void     SetFoil2Position(G4double v);
  G4double GetFoil2Position(void) const {return fFoil2Position;}

  virtual void     SetFoil2ThicknessFactor(G4double v);
  G4double GetFoil2ThicknessFactor(void) const {return fFoil2ThicknessFactor;}

  virtual void     AddCustomFoil2Vertex(G4double r, G4double z);

  virtual void     SetTargetName(G4String v);
  G4String GetTargetName(void) const {return fTargetName;}

  virtual void     SetTargetPosition(G4double v);
  G4double GetTargetPosition(void) const {return fTargetPosition;}

  virtual void     SetTargetRadius(G4double v);
  G4double GetTargetRadius(void) const {return fTargetRadius;}

  virtual void     SetTargetThickness(G4double v);
  G4double GetTargetThickness(void) const {return fTargetThickness;}

  virtual void     SetTargetMaterial(G4String v);
  G4String GetTargetMaterial(void) const {return fTargetMaterial;}

  virtual void     SetTargetMaterial2(G4String v);
  G4String GetTargetMaterial2(void) const {return fTargetMaterial2;}

  virtual void AddCustomTargetVertex(G4double r, G4double z);

  virtual void     SetFFName(G4String v);
  G4String GetFFName(void) const {return fFlatteningFilterName;}

  virtual void          SetFFOffset(G4ThreeVector v);
  G4ThreeVector GetFFOffset(void) const {return fFFOffset;}

  virtual void     SetCopperThicknessFactor(G4double v);
  G4double GetCopperThicknessFactor(void) const {return fCopperThicknessFactor;}
  
  virtual void   SetSimulateCollimators(G4bool v);
  G4bool GetSimulateCollimators(void) const {return fSimulateCollimators;}

  virtual void   SetSimulateShielding(G4bool v);
  G4bool GetSimulateShielding(void) const {return fSimulateShielding;}

  virtual void   SetSimulateVault(G4bool v);
  G4bool GetSimulateVault(void) const {return fSimulateVault;}

  virtual void     SetCollimatorRotation(G4double v);
  G4double GetCollimatorRotation(void) const {return fCollRot;}

  virtual void     SetJawPositionX(G4int jawint, G4double v, G4bool force=false);
  G4double GetJawPositionX(G4int jawint);
  virtual void     SetJawPositionY(G4int jawint, G4double v, G4bool force=false);
  G4double GetJawPositionY(G4int jawint);

  virtual void     SetJawOffset(G4double v);
  G4double GetJawOffset(void) {return fJawOffset;}

  virtual void     SetMLCHDLeafPosition(G4int index, G4double pos);
  virtual void     SetMLCNDS120Position(G4int index, G4double pos);
  G4double GetMLCLeafPosition(G4int index);

  virtual void SetProximalPosition(G4int, G4double);
  virtual void SetDistalPosition(G4int, G4double);

  virtual void SetMLCDensity(G4double);

  virtual void BuildCone();

  virtual void SetConeSize(G4double);
  G4double     GetConeSize(void) {return fCone_size;}

  virtual void SetConeElevation(G4double);
  G4double     GetConeElevation(void) {return fCone_Zpos;}

  virtual void SetConeThickness(G4double);
  G4double     GetConeThickness(void) {return fCone_thick;}

  virtual void SetConeOuterRadius(G4double);
  G4double     GetConeOuterRadius(void) {return fCone_OuterRadius;}


  // scatterer for TSI
  void     SetTSIPlateUse(G4bool v);
  G4bool   GetTSIPlateUse(void) const {return fTSIPlateUse;}

  void     SetTSIPlatePosition(G4double v);
  G4double GetTSIPlatePosition(void) const {return fTSIPlatePosition;}

  void     SetTSIPlateSideLength(G4double v);
  G4double GetTSIPlateSideLength(void) const 
    {return fTSIPlateSideLength;}

  void     SetTSIPlateThickness(G4double v);
  G4double GetTSIPlateThickness(void) const 
    {return fTSIPlateThickness;}

  void     SetTSIPlateMaterial(G4String v);
  G4String GetTSIPlateMaterial(void) const {return fTSIPlateMaterial;}

  void BuildTSIPlate(void); 


  virtual void     SetGantryRotation(G4double v, G4double theta=0.);
  G4double GetGantryRotation(void) const {return fGantryRot;}
  G4double GetGantryTheta(void) const {return fGantryTheta;}
  
  virtual void     SetApplicatorName(G4String v);
  G4String GetApplicatorName(void) const {return fApplicatorName;}

  virtual void     SetBuildCutOut(G4bool v);
  G4bool   GetBuildButOut(void) const {return fBuildCutOut;}

  virtual void  AddCutOutVertex(G4double x, G4double y);
  
  virtual void     SetCutOutThickness(G4double v);
  G4double GetCutOutThickness(void) const {return fCutOutThickness;}

  virtual void     SetCutOutBevelFactor(G4double v);
  G4double GetCutOutBevelFactor(void) const {return fCutOutBevelFactor;}

  virtual void     SetCutOutMaterial(G4String v);
  G4String GetCutOutMaterial(void) const {return fCutOutMaterial;}

  virtual void     SetOutputFilename(G4String v);
  G4String GetOutputFilename(void) const {return fOutputFilename;}

  virtual void RemoveMLC(void);
  virtual void SetMLC(G4String type);

  virtual void AddMaterial(G4String name, G4double density, G4String basemat);

  const G4VPhysicalVolume* GetBackscatterKiller(void) {return fBackScatter;}
  const G4VPhysicalVolume* GetTargetBTLow(void) {return fTargetBlockTopLow;}
  const G4VPhysicalVolume* GetShieldColl(void) {return fShieldColl;}
  const G4VPhysicalVolume* GetPrimColl(void) {return fPrimColl;}

  const G4VPhysicalVolume* GetJawY1(void) {return fJawY1;}
  const G4VPhysicalVolume* GetJawY2(void) {return fJawY2;}
  const G4VPhysicalVolume* GetJawX1(void) {return fJawX1;}
  const G4VPhysicalVolume* GetJawX2(void) {return fJawX2;}

  void SetVerbosity(G4int v) {fVerbosity = v;}

  std::unordered_map<int, TraversedGeometry> m_id_to_traversed;

  G4LogicalVolume* fLogicWorld;
protected:
  //static G4ThreadLocal G4bool fConstructedSDandField;

  // Data members
  //TB02_DetectorMessenger* fMessenger; 

  G4bool        fSimulateCollimators;
  G4bool        fSimulateShielding;  //simulate shielding outside of beam
                                     // esp. for neutrons
  G4bool        fSimulateVault;      // simulate vault etc.
 
  G4String fBeamType;   // electron or xray

  //TODO group by element, rather than log/phys
  // Logical and physical Volumes 
  //G4LogicalVolume* fLogicWorld;
  G4LogicalVolume* fLogicplane;
  G4PVPlacement*   fPlane;

  // mother volume for linac (gantry rotation rotates this)
  G4LogicalVolume* fGantry_LV;
  G4PVPlacement*   fGantry;
  G4double         fGantryPos;  //distance from isocenter
  G4double         fGantryRot;  // angular position
  G4double         fGantryTheta;  // angular position, about the other axis
                                  // ie longitude

  G4Tubs*          fEfoil1_tubs;
  G4PVPlacement*   fEfoil1;
  G4LogicalVolume* fEfoil1_LV;
  G4LogicalVolume* fEfoil1_holder_LV;

  // backscatter killer
  G4LogicalVolume*   fBackScatter_LV;
  G4VPhysicalVolume* fBackScatter;


  // vacuum of orbit chamber
  G4LogicalVolume* fOrbitVacuum_LV;
  G4PVPlacement*   fOrbitVacuum;

  // Be Window
  G4bool    fBuildBeWindow;

  // the target
  G4String fTargetName;
  G4double fTargetThickness;
  G4double fTargetRadius;
  G4double fTargetPosition;
  G4String fTargetMaterial;
  //G4double fTargetpos;   // position of mother in world; 
                  // and daughters relatvie to mother


  G4LogicalVolume* fTargetMother_LV;
  G4PVPlacement*   fTargetMother;
  
  //low energy target 
  G4LogicalVolume* fBeWinLow_LV;
  G4PVPlacement*   fBeWinLow;
  G4LogicalVolume* fVacuumChamber_LV;
  G4PVPlacement*   fVacuumChamber;
  G4LogicalVolume* fTargetLow_LV;
  G4PVPlacement*   fTargetLow;
  G4LogicalVolume* fTargetNicoroLow_LV;
  G4PVPlacement*   fTargetNicoroLow;
  G4LogicalVolume* fTargetBlockTopLow_LV;
  //G4PVPlacement*   fTargetBlockTopLow;
  G4VPhysicalVolume*   fTargetBlockTopLow;
  G4LogicalVolume* fXrayWin_LV;
  G4PVPlacement*   fXrayWin;

  //medium energy target
  G4LogicalVolume* fTargetBlockTopMed_LV;
  G4PVPlacement*   fTargetBlockTopMed;

  //high energy target
  G4LogicalVolume* fTargetHigh_LV;
  G4PVPlacement*   fTargetHigh;
  G4LogicalVolume* fTargetNicoroHigh_LV;
  G4PVPlacement*   fTargetNicoroHigh;
  G4LogicalVolume* fTargetBlockTopHigh_LV;
  G4PVPlacement*   fTargetBlockTopHigh;
  G4LogicalVolume* fTargetWafer_LV;
  G4PVPlacement*   fTargetWafer;
  G4LogicalVolume* fTargetBlockBottomHigh_LV;
  G4PVPlacement*   fTargetBlockBottomHigh;
  
  //imaging target
  G4LogicalVolume* fTargetBlockTopImage_LV;
  G4PVPlacement*   fTargetBlockTopImage;

  //custom target
  G4Tubs*          fTargetCustomTubs;
  G4LogicalVolume* fTargetCustom_LV;
  G4PVPlacement*   fTargetCustom;
  //second layer of custom target
  std::vector<G4double> fTargetCustom2_r;
  std::vector<G4double> fTargetCustom2_z;
  G4GenericPolycone*    fTargetCustomPC2;
  G4LogicalVolume*      fTargetCustom2_LV;
  G4PVPlacement*        fTargetCustom2;
  G4String              fTargetMaterial2;

  //mother for collimators
  G4double          fCollRot;  // rotation angle
  G4double          fCollPos;
  G4LogicalVolume*  fColl_LV;
  G4PVPlacement*    fColl;

  //shield coll
  G4VPhysicalVolume* fShieldColl;
  //prim coll
  G4VPhysicalVolume* fPrimColl;

  // jaws, MLC
  G4LogicalVolume* fUpperJawTrap_LV;
  G4LogicalVolume* fLowerJawTrap_LV;
  G4LogicalVolume* fBasePlate_LV;
  G4LogicalVolume* fMLC_LV; // the fake one (solid)
  G4LogicalVolume* fMylarWindow_LV;

  G4PVPlacement**  fMlc_PV;
  G4double         fMlc_x_size;
  G4LogicalVolume* fMlc_half_tar_A_LV;
  G4LogicalVolume* fMlc_half_tar_B_LV;
  G4LogicalVolume* fMlc_half_iso_A_LV;
  G4LogicalVolume* fMlc_half_iso_B_LV;
  G4LogicalVolume* fMlc_full_A_LV;
  G4LogicalVolume* fMlc_full_B_LV;
  // for HD120
  G4LogicalVolume* fMlc_quarter_tar_A_LV;
  G4LogicalVolume* fMlc_quarter_tar_B_LV;
  G4LogicalVolume* fMlc_quarter_iso_A_LV;
  G4LogicalVolume* fMlc_quarter_iso_B_LV;
  G4LogicalVolume* fMlc_outboard1_A_LV;
  G4LogicalVolume* fMlc_outboard60_A_LV;
  G4LogicalVolume* fMlc_outboard1_B_LV;
  G4LogicalVolume* fMlc_outboard60_B_LV;
 
  G4String fMlc_type;  // allowed values: NDS120, NDS120HD
  G4String fMlc_material; 
  G4String fMlc_material_NDS; 

  G4double fMlc_tipradius;
  G4double fMlc_height;

  G4LogicalVolume* fVault_LV;

  G4LogicalVolume* fSkullcap_pb_LV;
  G4PVPlacement*   fSkullcap_pb;
  G4LogicalVolume* fSkullcap_w_LV;
  G4PVPlacement*   fSkullcap_w;
  G4LogicalVolume* fYoke_LV;
  G4PVPlacement*   fYoke;
  G4LogicalVolume* fMainCoil_LV;
  G4PVPlacement*   fMainCoil1;
  G4PVPlacement*   fMainCoil2;
  G4LogicalVolume* fEars_LV;
  G4PVPlacement*   fEars;
  G4LogicalVolume* fMagnet_pole_LV;
  G4PVPlacement*   fMagnet_pole1;
  G4PVPlacement*   fMagnet_pole2;
  G4LogicalVolume* fOrbit_chamber_LV;
  G4PVPlacement*   fOrbit_chamber;
  G4LogicalVolume* fGantrySteel_LV;
  G4PVPlacement*   fGantrySteel;
  G4LogicalVolume* fCover_LV;
  G4PVPlacement*   fCover;

  G4double fSAD;  // distance from virtual source to isocenter

  G4PVPlacement* fPhysiWorld;
  G4VPhysicalVolume* fJawX1;
  G4VPhysicalVolume* fJawX2;
  G4VPhysicalVolume* fJawY1;
  G4VPhysicalVolume* fJawY2;
  G4PVPlacement* fBasePlate;
  G4PVPlacement* fBasePlate1A;
  G4PVPlacement* fBasePlate1B;
  G4PVPlacement* fBasePlate2;
  G4PVPlacement* fCornerClipper1;
  G4PVPlacement* fCornerClipper2;
  G4PVPlacement* fCornerClipper3;
  G4PVPlacement* fCornerClipper4;
  G4PVPlacement* fMLC1;
  G4PVPlacement* fMLC2;
  G4PVPlacement* fMylarWindow;
  G4PVPlacement* fVault;

  G4double fProx_extra_rot;
  G4double fDist_extra_rot;
  G4double fProx_extra_leaf_rot;
  G4double fDist_extra_leaf_rot;
  
  G4double fProx_spacing;
  G4double fDist_spacing;
  //chamfer: _ratio used to set the angle of the cut. It is the ratio
  //of the length of cut in leaf face to length of cut in leaf side.
  // _frac is the fraction of the face (per chamfer) removed, at the midplane
  G4double fDist_chamfer_ratio;
  G4double fDist_chamfer_frac;

  // srs cone
  G4bool fBuildCone;  
  G4double fCone_size;  // nominal diameter
  G4double fCone_Zpos;  // elevation of proximal surface relative to iso
  G4double fCone_thick; // length of cone along z
  G4double fCone_OuterRadius; // outer radius of cone

  // TSI scatterer
  G4bool   fTSIPlateUse;
  G4double fTSIPlateThickness;
  G4double fTSIPlateSideLength;
  G4double fTSIPlatePosition;
  G4String fTSIPlateMaterial;


  G4PVPlacement* targetA_phys;
  G4PVPlacement* targetB_phys;
  G4PVPlacement* UpperCollimator_phys;
  G4PVPlacement* FFL1_1PV;
  G4PVPlacement* FFL2_1PV;

  G4String fFoil1Name;
  G4String fFoil1Material;
  G4double fFoil1Thickness;
  G4double fFoil1Radius;
  G4double fFoil1Position;
  G4String fFoil2Name;
  G4String fFoil2Material;
  //G4double fFoil2Thickness;
  //G4double fFoil2Radius;
  G4double fFoil2Position;
  std::vector<G4double> fFoil2Custom_r;
  std::vector<G4double> fFoil2Custom_z;
  G4LogicalVolume* fEfoil2_LV;
  G4PVPlacement*   fEfoil2;
  //G4GenericPolycone*    fFoilCustomPC;



  G4String fFlatteningFilterName;
  G4ThreeVector fFFOffset;

  G4double fExitWindowThicknessFactor;
  G4double fFoil1ThicknessFactor;
  G4double fFoil2ThicknessFactor;

  // ion chamber
  G4bool fUseIonChamber;
  G4double fCopperThicknessFactor;
  //G4double fKaptonThicknessFactor;
  G4LogicalVolume* fICMother_LV;
  G4PVPlacement*   fICMother;

  ////    cu window
  G4double         fICCuHThick;  // nominal thickness
  G4Tubs*          fICWinCu_tubs;
  G4LogicalVolume* fICWinCu_LV;
  G4PVPlacement*   fICWinCu0;
  G4PVPlacement*   fICWinCu1;
  G4PVPlacement*   fICWinCu2;

  ////    hv electrode
  G4Tubs*          fICElCu_tubs;
  G4LogicalVolume* fICElCu_LV;
  G4PVPlacement*   fICElCu0;
  G4PVPlacement*   fICElCu1;

  ////    signal electrode
  G4Tubs*          fICElCusig_tubs;
  G4LogicalVolume* fICElCusig_LV;
  G4PVPlacement*   fICElCusig0;
  G4PVPlacement*   fICElCusig1;

  G4LogicalVolume* fIC_dose_LV;
  G4PVPlacement*   fIC_dose;
  G4PSDoseDeposit* fMonChamberDose;
 
  // applicator and cutout
  G4String fApplicatorName;
  G4bool   fBuildCutOut;
  //G4double fCutOutX, fCutOutY;
  std::vector<G4TwoVector> fCutOutVertices;
  G4double fCutOutThickness;
  G4String fCutOutMaterial;
  G4double fCutOutBevelFactor;  // the ratio of size of upper part of cutout,
                                // to lower
  G4LogicalVolume* fCutOut_LV;
  G4PVPlacement*   fCutOut;


  G4String fOutputFilename;

  G4double fFieldX1, fFieldX2, fFieldY1, fFieldY2;   // field sizes (jaws)

  G4double fJawOffset;

  //necessary for multithreading
  const G4Cache<G4MultiFunctionalDetector*> fMonitorChamberCache;

  G4Region* fTargetRegion;

 
  G4SDManager* fSDManager;
  
  G4bool fVis1, fVis2, fVis3, fVis4, fVisOff;

  G4int fVerbosity;

};
#endif
