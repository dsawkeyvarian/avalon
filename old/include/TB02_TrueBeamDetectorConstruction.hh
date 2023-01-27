/* This file is the CONFIDENTIAL and PROPRIETARY information of
 * Varian Medical Systems, Inc., and is not to be distributed.
 * 
 * Copyright (c) 2017 Varian Medical Systems, Inc.
 * 
 * For information, contact Daren Sawkey  daren.sawkey@varian.com
 */

#ifndef TB02_TrueBeamDetectorConstruction_h
#define TB02_TrueBeamDetectorConstruction_h 1

#include "TB02_BaseDetectorConstruction.hh"
//#include "TB02_DetectorMessenger.hh"

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4Cache.hh"

#include "G4ExtrudedSolid.hh"  // provides G4TwoVector

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

//class TB02_PhaseSpaceWriter;
class G4PSDoseDeposit;


class TB02_TrueBeamDetectorConstruction : public TB02_BaseDetectorConstruction
{
public:
    G4bool fBuildWaterPhantom;
    G4ThreeVector fPhantomPosition;
    G4ThreeVector fPhantomBoxSize;
  // constructor and destructor.
  TB02_TrueBeamDetectorConstruction(G4String gdml_folder,G4String stl_folder, G4bool construct_phantom, G4ThreeVector phantom_position, G4ThreeVector phantom_box_size);
  virtual ~TB02_TrueBeamDetectorConstruction();

  // virtual method from G4VUserDetectorConstruction.
  virtual G4VPhysicalVolume* Construct() override;
  virtual void GenerateTraversedTable() override;
  void ConstructSDandField() override;

  void BuildBeWindow();

  void BuildElectronGeometry(G4bool build = true) override;
  void BuildPhotonGeometry(G4bool build = true) override;

  void BuildBackscatterShield(G4bool build = true) override;
  void BuildTarget(G4bool build = true) override;
  void BuildPrimaryCollimator(G4bool build = true) override;
  void BuildYStageShield(G4bool build = true) override;
  void BuildFlatteningFilter(G4bool build = true) override;

  void BuildFoil1(G4bool build = true) override;
  void BuildFoil2(G4bool build = true) override;

  void BuildShielding(G4bool build = true) override;
  void BuildShieldingCollimator() override;
  void BuildCollimatorMother() override; //mother volume for shielding coll 
                                // and collim (jaw+mlc etc)
  void BuildCollimators(G4bool build = true) override;
                            // builds jaws and tx head components after:
                            // incl base plate, MLC, mylar window
  void BuildApplicator() override;
  void BuildCutOut() override;
  //void BuildPhaseSpace();
  //void RemovePhaseSpace(G4int index);
  void BuildVault(G4bool build = true) override;


  void     SetBeamType(G4String v) override;
  G4String GetBeamType(void) const {return fBeamType;}

  void   SetBuildBeWindow(G4bool build = true) override;
  G4bool GetBuildBeWindow(void) const {return fBuildBeWindow;}

  void   SetUseIonChamber(G4bool v) override;
  G4bool GetUseIonChamber(void) const {return fUseIonChamber;}

  void     SetExitWindowThicknessFactor(G4double v) override;
  G4double GetExitWindowThicknessFactor(void) const 
    {return fExitWindowThicknessFactor;}

  void     SetFoil1Name(G4String s) override;
  G4String GetFoil1Name(void) const {return fFoil1Name;}

  void     SetFoil1Material(G4String v) override;
  G4String GetFoil1Material(void) const {return fFoil1Material;}

  void     SetFoil1Thickness(G4double v) override;
  G4double GetFoil1Thickness(void) const {return fFoil1Thickness;}

  void     SetFoil1Radius(G4double v) override;
  G4double GetFoil1Radius(void) const {return fFoil1Radius;}

  void     SetFoil1Position(G4double v) override;
  G4double GetFoil1Position(void) const {return fFoil1Position;}

  void     SetFoil1ThicknessFactor(G4double v) override;
  G4double GetFoil1ThicknessFactor(void) {return fFoil1ThicknessFactor;}
 
  void     SetFoil2Name(G4String s) override;
  G4String GetFoil2Name(void) {return fFoil2Name;}

  void     SetFoil2Material(G4String v) override;
  G4String GetFoil2Material(void) const {return fFoil2Material;}

  //void     SetFoil2Thickness(G4double v);
  //G4double GetFoil2Thickness(void) const {return fFoil2Thickness;}

  //void     SetFoil2Radius(G4double v);
  //G4double GetFoil2Radius(void) const {return fFoil2Radius;}

  void     SetFoil2Position(G4double v) override;
  G4double GetFoil2Position(void) const {return fFoil2Position;}

  void     SetFoil2ThicknessFactor(G4double v) override;
  G4double GetFoil2ThicknessFactor(void) const {return fFoil2ThicknessFactor;}

  void     AddCustomFoil2Vertex(G4double r, G4double z) override;

  void     SetTargetName(G4String v) override;
  G4String GetTargetName(void) const {return fTargetName;}

  void     SetTargetPosition(G4double v) override;
  G4double GetTargetPosition(void) const {return fTargetPosition;}

  void     SetTargetRadius(G4double v) override;
  G4double GetTargetRadius(void) const {return fTargetRadius;}

  void     SetTargetThickness(G4double v) override;
  G4double GetTargetThickness(void) const {return fTargetThickness;}

  void     SetTargetMaterial(G4String v) override;
  G4String GetTargetMaterial(void) const {return fTargetMaterial;}

  void     SetTargetMaterial2(G4String v) override;
  G4String GetTargetMaterial2(void) const {return fTargetMaterial2;}

  void AddCustomTargetVertex(G4double r, G4double z) override;

  void     SetFFName(G4String v) override;
  G4String GetFFName(void) const {return fFlatteningFilterName;}

  void          SetFFOffset(G4ThreeVector v) override;
  G4ThreeVector GetFFOffset(void) const {return fFFOffset;}

  void     SetCopperThicknessFactor(G4double v) override;
  G4double GetCopperThicknessFactor(void) const {return fCopperThicknessFactor;}
  
  void   SetSimulateCollimators(G4bool v) override;
  G4bool GetSimulateCollimators(void) const {return fSimulateCollimators;}

  void   SetSimulateShielding(G4bool v) override;
  G4bool GetSimulateShielding(void) const {return fSimulateShielding;}

  void     SetJawPositionX(
              G4int jawint, G4double v, G4bool force=false) override;
  G4double GetJawPositionX(G4int jawint);
  void     SetJawPositionY(
              G4int jawint, G4double v, G4bool force=false) override;
  G4double GetJawPositionY(G4int jawint);

  void     SetJawOffset(G4double v) override;
  G4double GetJawOffset(void) {return fJawOffset;}

  void     SetMLCHDLeafPosition(G4int index, G4double pos) override;
  void     SetMLCNDS120Position(G4int index, G4double pos) override;
  G4double GetMLCLeafPosition(G4int index);

  void     SetMLCDensity(G4double v) override;

  void     SetApplicatorName(G4String v) override;
  G4String GetApplicatorName(void) const {return fApplicatorName;}

  void     SetBuildCutOut(G4bool v) override;
  G4bool   GetBuildButOut(void) const {return fBuildCutOut;}

  void  AddCutOutVertex(G4double x, G4double y) override;
  
  void     SetCutOutThickness(G4double v) override;
  G4double GetCutOutThickness(void) const {return fCutOutThickness;}

  void     SetCutOutBevelFactor(G4double v) override;
  G4double GetCutOutBevelFactor(void) const {return fCutOutBevelFactor;}

  void     SetCutOutMaterial(G4String v) override;
  G4String GetCutOutMaterial(void) const {return fCutOutMaterial;}

  void BuildCADMLC(void);
  void BuildHDMLC(void);
  void RemoveMLC(void) override;
  void SetMLC(G4String type) override;

  const G4VPhysicalVolume* GetBackscatterKiller(void) {return fBackScatter;}
  const G4VPhysicalVolume* GetTargetBTLow(void) {return fTargetBlockTopLow;}
  const G4VPhysicalVolume* GetShieldColl(void) {return fShieldColl;}
  const G4VPhysicalVolume* GetPrimColl(void) {return fPrimColl;}

  const G4VPhysicalVolume* GetJawY1(void) {return fJawY1;}
  const G4VPhysicalVolume* GetJawY2(void) {return fJawY2;}
  const G4VPhysicalVolume* GetJawX1(void) {return fJawX1;}
  const G4VPhysicalVolume* GetJawX2(void) {return fJawX2;}


private:
  G4String m_gdml_folder;
  G4String m_stl_folder;

};
#endif
