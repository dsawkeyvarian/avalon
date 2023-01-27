#include "treatment_heads.h"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"

namespace nohead {
G4VPhysicalVolume* TreatmentHeadDetector::Construct() {
	//We actually want to store physical volumes for rotations
	G4NistManager* NISTman = G4NistManager::Instance();
	auto world_material = NISTman->FindOrBuildMaterial(m_background_material);

	G4LogicalVolume* world_logical = nullptr;
	G4VPhysicalVolume* world_physical = nullptr;
	//World volume
	{
		G4ThreeVector world_size = G4ThreeVector(1000. * cm, 1000. * cm, 1000. * cm);
		G4Box* world_box = new G4Box("world", world_size.x() / 2., world_size.y() / 2., world_size.z() / 2.);
		world_logical = new G4LogicalVolume(world_box, world_material, "world", 0, 0, 0);
		world_physical = new G4PVPlacement(0, G4ThreeVector(), world_logical, "world", 0, false, 0);
	}
	VolumeNameAndTraversalFlag NoHeadTable[] = {{"world", TraversedGeometry::NONE}};
	m_id_to_traversed = generateVolumeToTraversed(NoHeadTable, sizeof(NoHeadTable) / sizeof(NoHeadTable[0]), world_logical);
	return world_physical;
}

void TreatmentHeadDetector::ConstructSDandField() {}
}