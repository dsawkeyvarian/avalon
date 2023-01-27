/* This file is the CONFIDENTIAL and PROPRIETARY information of
 * Varian Medical Systems, Inc., and is not to be distributed.
 * 
 * Copyright (c) 2017 Varian Medical Systems, Inc.
 * 
 * For information, contact Daren Sawkey  daren.sawkey@varian.com
 */

#ifndef TB02_HalcyonDetectorConstruction_h
#define TB02_HalcyonDetectorConstruction_h 1

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

class TB02_HalcyonDetectorConstruction : public TB02_BaseDetectorConstruction
{
public:
  G4bool fBuildWaterPhantom;
  G4ThreeVector fPhantomPosition;
  G4ThreeVector fPhantomBoxSize;
  // constructor and destructor.
  TB02_HalcyonDetectorConstruction(G4String gdml_folder,G4bool construct_phantom,G4ThreeVector phantom_position, G4ThreeVector phantom_box_size);
  virtual ~TB02_HalcyonDetectorConstruction();

  // virtual method from G4VUserDetectorConstruction.
  virtual G4VPhysicalVolume* Construct() override;
  virtual void GenerateTraversedTable() override;
  void ConstructSDandField() override;

  void BuildVacuum();

  virtual void BuildPhotonGeometry(G4bool build = true);

  virtual void BuildBackscatterKiller(G4bool build = true);
  void BuildTarget(G4bool build = true) override;
  void BuildPrimaryCollimator(G4bool build = true) override;
  void BuildBuildupPlate();
  void BuildMCBackscatterPlate();
  void BuildAlPlate();

  void ProximalLeafBank();
  void DistalLeafBank();
  void BuildDistalFixed();
  void BuildProximalShield();
  void BuildAlWindow();

  void BuildSecondaryCollimator();
  virtual void BuildCollimatorMother(); //mother volume for shielding coll 
                                  // and collim (jaw+mlc etc)
  void BuildPlasticMLCBox();

  G4double GetActualProximalPosition(G4int index, G4double nom);
  virtual void SetProximalPosition(G4int index, G4double nom) override;
  G4double GetProximalPosition(G4int index) {return fProximalPos[index];}
  void     SetAllProximalPositions(G4double nom);

  G4double GetActualDistalPosition(G4int index, G4double nom);
  virtual void SetDistalPosition(G4int index, G4double nom) override;
  G4double GetDistalPosition(G4int index) {return fDistalPos[index];}
  void     SetAllDistalPositions(G4double nom);

  //void RemoveMLC(void) override;
  //void SetMLC(G4String type) override;

  const G4VPhysicalVolume* GetBackscatterKiller(void) {return fBackScatter;}
  const G4VPhysicalVolume* GetShieldColl(void) {return fShieldColl;}
  const G4VPhysicalVolume* GetPrimColl(void) {return fPrimColl;}


private:
  G4String m_gdml_folder;
  G4double         fMlc_height;
  G4double         fMlc_tipradius;

  //// for H, a vector of leaves PV
  std::vector<G4PVPlacement*>  fProximalLeaves;
  std::vector<G4double>        fProximalPos;  // their positions
      // positive is open, negative is overtravel
  G4double fProximalZPos;  // position of center in Z
  G4double fLeafHLength;  // half length, along direction of travel 
                          // for both leaves
  G4int fNProx;  // number of proximal leaves
  G4int fNDist;  // number of distal leaves 

  std::vector<G4PVPlacement*>  fDistalLeaves;
  std::vector<G4double>        fDistalPos;  // their positions
      // positive is open, negative is overtravel
  G4double fDistalZPos;  // position of center in Z
  G4double fDistalHLength;  // half length, along direction of travel


};
#endif
