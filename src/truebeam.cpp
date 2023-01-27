#include "G4PhysListFactory.hh"
#include "G4SystemOfUnits.hh"

#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"

#include "G4SDManager.hh"
#include "G4PSDoseDeposit3D.hh"
#include "G4PVReplica.hh"
#include "G4NistManager.hh"
#include "G4UnitsTable.hh"
#include "G4Tubs.hh"

#include <iostream>
#include <fstream>
#include <chrono>
#include <unordered_map>
#include <atomic>
#include <filesystem>
namespace fs = std::filesystem;

#include "G4PhysicalConstants.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4GeometryManager.hh"
#include "G4GDMLParser.hh"
#include "G4GenericPolycone.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "stl.h"

#include "treatment_heads.h"
#include "plan.h"

static const double g_SAD = 100.0 * cm;

//These are general, not only truebeam
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

namespace truebeam {

enum class Target
{
	LOW_ENERGY,
	HIGH_ENERGY,
	NUM_TARGETS
};

enum class FlatteningFilter
{
    OPEN,
	ENERGY_6X,
    NUM_FLATTENING_FILTERS
};

//HD MLC 
MLCLeafVolumes buildHDMLC(G4LogicalVolume* parent_logical, const double parent_world_z, const std::filesystem::path& stl_folder)
{
	//Visualization
	const auto bank_X2_vis = new G4VisAttributes(G4Colour(0.3, 0.8, 0.3, 0.9));
	bank_X2_vis->SetForceSolid(true);
	const auto bank_X1_vis = new G4VisAttributes(G4Colour(0.3, 0.3, 0.8, 0.9));
	bank_X1_vis->SetForceSolid(true);
	//Order of files is meaningful so that leaf position indices are correct
	static const char* filenames_bank_X2[] =
	{	"060_100025506-13.stl",
		"059_100022597-12.stl",
		"058_100022596-14.stl",
		"057_100022597-19.stl",
		"056_100022596-11.stl",
		"055_100022597-14.stl",
		"054_100022596-12A.stl",
		"053_100022597-12.stl",
		"052_100022596-14.stl",
		"051_100022597-19.stl",
		"050_100022596-11.stl",
		"049_100022597-14.stl",
		"048_100022596-13.stl",
		"047_100022597-12.stl",
		"046_p1000882004.stl",
		"045_p1000883001.stl",
		"044_p1000882002.stl",
		"043_p1000883004.stl",
		"042_p1000882003.stl",
		"041_p1000883002.stl",
		"040_p1000882001.stl",
		"039_p1000883003.stl",
		"038_p1000882004.stl",
		"037_p1000883001.stl",
		"036_p1000882002.stl",
		"035_p1000883004.stl",
		"034_p1000882003.stl",
		"033_p1000883002.stl",
		"032_p1000882001.stl",
		"031_p1000883003.stl",
		"030_p1000882004.stl",
		"029_p1000883001.stl",
		"028_p1000882002.stl",
		"027_p1000883004.stl",
		"026_p1000882003.stl",
		"025_p1000883002.stl",
		"024_p1000882001.stl",
		"023_p1000883003.stl",
		"022_p1000882004.stl",
		"021_p1000883001.stl",
		"020_p1000882002.stl",
		"019_p1000883004.stl",
		"018_p1000882003.stl",
		"017_p1000883002.stl",
		"016_p1000882001.stl",
		"015_p1000883003.stl",
		"014_100022596-14.stl",
		"013_100022597-19.stl",
		"012_100022596-12.stl",
		"011_100022597-12.stl",
		"010_100022596-15.stl",
		"009_100022597-19.stl",
		"008_100022596-11.stl",
		"007_100022597-13.stl",
		"006_100022596-12.stl",
		"005_100022597-12.stl",
		"004_100022596-14.stl",
		"003_100022597-19.stl",
		"002_100022596-11.stl",
		"001_100025507-13.stl"
		
	};

	static const char* filenames_bank_X1[] = 
	{	"160_100025506-14.stl",
		"159_100022597-16.stl",
		"158_100022596-19.stl",
		"157_100022597-20.stl",
		"156_100022596-16.stl",
		"155_100022597-18.stl",
		"154_100022596-17.stl",
		"153_100022597-16.stl",
		"152_100022596-19.stl",
		"151_100022597-20.stl",
		"150_100022596-16.stl",
		"149_100022597-18.stl",
		"148_100022596-18.stl",
		"147_100022597-16.stl",
		"146_100022025-26.stl",
		"145_100022026-23.stl",
		"144_100022025-24.stl",
		"143_100022026-26.stl",
		"142_100022025-25.stl",
		"141_100022026-24.stl",
		"140_100022025-23.stl",
		"139_100022026-25.stl",
		"138_100022025-26.stl",
		"137_100022026-23.stl",
		"136_100022025-24.stl",
		"135_100022026-26.stl",
		"134_100022025-25.stl",
		"133_100022026-24.stl",
		"132_100022025-23.stl",
		"131_100022026-25.stl",
		"130_100022025-26.stl",
		"129_100022026-23.stl",
		"128_100022025-24.stl",
		"127_100022026-26.stl",
		"126_100022025-25.stl",
		"125_100022026-24.stl",
		"124_100022025-23.stl",
		"123_100022026-25.stl",
		"122_100022025-26.stl",
		"121_100022026-23.stl",
		"120_100022025-24.stl",
		"119_100022026-26.stl",
		"118_100022025-25.stl",
		"117_100022026-24.stl",
		"116_100022025-23.stl",
		"115_100022026-25.stl",
		"114_100022597-19.stl",
		"113_100022597-20.stl",
		"112_100022596-17.stl",
		"111_100022597-16.stl",
		"110_100022596-20.stl",
		"109_100022597-20.stl",
		"108_100022596-16.stl",
		"107_100022597-17.stl",
		"106_100022596-17.stl",
		"105_100022597-16.stl",
		"104_100022596-19.stl",
		"103_100022597-20.stl",
		"102-100022596-16.stl",
		"101_100025507-14.stl"
	};



	auto mlc_volumes = MLCLeafVolumes();
	//Bank X2 / A
	{
		size_t num_leaves = sizeof(filenames_bank_X2) / sizeof(*filenames_bank_X2);
		for (size_t i = 0; i < num_leaves; ++i) {
			std::filesystem::path filepath = stl_folder / "A" / filenames_bank_X2[i]; 
			const G4ThreeVector offset(-115.0, 0.0, 0.0);
            G4VSolid* leaf_solid = getTessalatedSolidFromSTL(filenames_bank_X2[i], filepath, offset);
			G4LogicalVolume* leaf_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("W95"), "mlc_leaf", 0, 0, 0);
			leaf_logical->SetVisAttributes(bank_X2_vis);
			auto volume = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_world_z),
				leaf_logical, "hd_leaf", parent_logical, false, 0);
			mlc_volumes.bank_X2.push_back(volume);
		}
	}
	//Bank X1 / B
	{
		size_t num_leaves = sizeof(filenames_bank_X1) / sizeof(*filenames_bank_X1);
		for (size_t i = 0; i < num_leaves; ++i) {
			std::filesystem::path filepath = stl_folder / "B" / filenames_bank_X1[i]; 
			const G4ThreeVector offset(110.497, 0.0, 0.0);
            G4VSolid* leaf_solid = getTessalatedSolidFromSTL(filenames_bank_X1[i], filepath, offset);
			G4LogicalVolume* leaf_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("W95"), "mlc_leaf", 0, 0, 0);
			leaf_logical->SetVisAttributes(bank_X1_vis);
			auto volume = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_world_z),
				leaf_logical, "hd_leaf", parent_logical, false, 0);
			mlc_volumes.bank_X1.push_back(volume);
		}
	}
	return mlc_volumes;
}

void setHDMLCLeafPositions(const MLCLeafVolumes& volumes, const LeafPositions& positions)
{
	// position is positive for non-overtravel, for both leaf banks
	const size_t num_leaves_in_bank = 60;
	const bool valid = (volumes.bank_X2.size() == num_leaves_in_bank) || (volumes.bank_X1.size() == num_leaves_in_bank) 
		|| (positions.bank_X2.size() == num_leaves_in_bank) || (positions.bank_X1.size() == num_leaves_in_bank);
	if (!valid)
		throw std::runtime_error("Number of leaf volumes and leaf positions is inconsistent");
	
	const std::vector<double>* bank_positions[] = { &positions.bank_X2, &positions.bank_X1 };
	const std::vector<G4VPhysicalVolume*>* bank_volumes[] = { &volumes.bank_X2, &volumes.bank_X1 };
	const G4int num_pos = 41;
	double physical_position = 0.0;
	/* Tables from machine
	G4double mlc_nominal[num_pos] = 
		{ 20.0*cm,  19.0*cm,  18.0*cm,  17.0*cm,  16.0*cm,  15.0*cm,  14.0*cm,
		13.0*cm,  12.0*cm,  11.0*cm,  10.0*cm,   9.0*cm,   8.0*cm,   7.0*cm,
		6.0*cm,   5.0*cm,   4.0*cm,   3.0*cm,   2.0*cm,   1.0*cm,   0.0*cm,
		-1.0*cm,  -2.0*cm,  -3.0*cm,  -4.0*cm,  -5.0*cm,  -6.0*cm,  -7.0*cm,
		-8.0*cm,  -9.0*cm, -10.0*cm, -11.0*cm, -12.0*cm, -13.0*cm, -14.0*cm,
		-15.0*cm, -16.0*cm, -17.0*cm, -18.0*cm, -19.0*cm, -20.0*cm};
	G4double mlc_actual[num_pos] =
		{ 20.6627*cm,  19.5967*cm,  18.5343*cm,  17.4754*cm,  16.4201*cm,
		15.3683*cm,  14.3200*cm,  13.2752*cm,  12.2339*cm,  11.1960*cm,
		10.1616*cm,   9.1305*cm,   8.1028*cm,   7.0785*cm,   6.0575*cm,
		5.0398*cm,   4.0254*cm,   3.0142*cm,   2.0063*cm,   1.0016*cm,
		0.0000*cm,  -0.9984*cm,  -1.9938*cm,  -2.9860*cm,  -3.9752*cm,
		-4.9614*cm,  -5.9446*cm,  -6.9249*cm,  -7.9022*cm,  -8.8766*cm,
		-9.8482*cm, -10.8170*cm, -11.7830*cm, -12.7462*cm, -13.7067*cm,
		-14.6645*cm, -15.6196*cm, -16.5721*cm, -17.5220*cm, -18.4694*cm,
		-19.4142*cm};
	*/
	for (size_t j = 0; j < 2; ++j) {
		for (size_t i = 0; i < num_leaves_in_bank; ++i) {
			double nominal_position = bank_positions[j]->at(i);
			if (j == 1) nominal_position *= (-1);
			physical_position = nominal_position + 160 * mm * (std::sqrt(1 + std::pow(nominal_position / g_SAD * mm, 2)) - 1)* (1000/509.952);
			if (j == 1) physical_position *= (-1);
			physical_position *= 509.952/1000.; 
			auto volume = bank_volumes[j]->at(i);
			auto translation = volume->GetTranslation();
			translation.setX(physical_position);
			volume->SetTranslation(translation);
		}
	}
}

//MLC leaf z-coordinate origin is at isocenter in local_to_world
//parent_to_world transform is used to compute local_to_parent = parent_to_world.inverse()*local_to_world
//where local_to_parent transform is used to place leafs inside parent volume
MLCLeafVolumes buildMilleniumMLC(G4LogicalVolume* parent_logical, const double parent_world_z, const std::filesystem::path& stl_folder)
{
	//Visualization
	const auto bank_X2_vis = new G4VisAttributes(G4Colour(0.3, 0.8, 0.3, 0.9));
	bank_X2_vis->SetForceSolid(true);
	const auto bank_X1_vis = new G4VisAttributes(G4Colour(0.3, 0.3, 0.8, 0.9));
	bank_X1_vis->SetForceSolid(true);
	//Millenium MLC material
	G4NistManager* NISTman = G4NistManager::Instance();
	G4Element* elW = NISTman->FindOrBuildElement("W");
	G4Element* elNi = NISTman->FindOrBuildElement("Ni");
	G4Element* elFe = NISTman->FindOrBuildElement("Fe");
	// 92.5% density tungsten
	const double density = 17.6 * g / cm3;
	G4Material* W92_5 = new G4Material("W92_5", density, 3);
	W92_5->AddElement(elW, 92.5 * perCent);
	W92_5->AddElement(elNi, 5.25 * perCent);
	W92_5->AddElement(elFe, 2.25 * perCent);
	//Order of files is meaningful so that leaf position indices are correct
  static const char* filenames_bank_X2[] =
  {
		"100026960-01-ao1-01.stl",
		"1105335-16-af3-02.stl",
		"1105335-14-af1-03.stl",
		"1105335-15-af2-04.stl",
		"1105335-16-af3-05.stl",
		"1105335-14-af1-06.stl",
		"1105335-15-af2-07.stl",
		"1105335-16-af3-08.stl",
		"1105335-14-af1-09.stl",
		"1105335-15-af2-10.stl",
		"1105333-15-at1-11.stl",
		"1105334-14-ai1-12.stl",
		"1105333-17-at3-13.stl",
		"1105334-16-a13-14.stl",
		"1105333-16-at2-15.stl",
		"1105334-14-ai1-16.stl",
		"1105333-18-at4-17.stl",
		"1105334-15-at2-18.stl",
		"1105333-15-at1-19.stl",
		"1105334-14-ai1-20.stl",
		"1105333-17-at3-21.stl",
		"1105334-16-ai3-22.stl",
		"1105333-16-at2-23.stl",
		"1105334-14-ai1-24.stl",
		"1105333-18-at4-25.stl",
		"1105334-15-ai2-26.stl",
		"1105333-15-at1-27.stl",
		"1105334-14-ai1-28.stl",
		"1105333-17-at3-29.stl",
		"1105334-16-ai3-30.stl",
		"1105333-16-at2-31.stl",
		"1105334-14-ai1-32.stl",
		"1105333-18-at4-33.stl",
		"1105334-15-ai2-34.stl",
		"1105333-15-at1-35.stl",
		"1105334-14-ai1-36.stl",
		"1105333-17-at3-37.stl",
		"1105334-16-ai3-38.stl",
		"1105333-16-at2-39.stl",
		"1105334-14-ai1-40.stl",
		"1105333-18-at4-41.stl",
		"1105334-15-ai2-42.stl",
		"1105333-15-at1-43.stl",
		"1105334-14-ai1-44.stl",
		"1105333-17-at3-45.stl",
		"1105334-16-ai3-46.stl",
		"1105333-16-at2-47.stl",
		"1105334-14-ai1-48.stl",
		"1105333-18-at4-49.stl",
		"1105334-15-ai2-50.stl",
		"1105335-16-af3-51.stl",
		"1105335-14-af1-52.stl",
		"1105335-15-af2-53.stl",
		"1105335-16-af3-54.stl",
		"1105335-14-af1-55.stl",
		"1105335-15-af2-56.stl",
		"1105335-16-af3-57.stl",
		"1105335-14-af1-58.stl",
		"1105335-15-af2-59.stl",
		"100026959-02-ao60-60.stl",
  };

	static const char* filenames_bank_X1[] = 
	{
		"100026960-02-bo1-01.stl",
		"1105335-19-bf3-02.stl",
		"1105335-17-bf1-03.stl",
		"1105335-18-bf2-04.stl",
		"1105335-19-bf3-05.stl",
		"1105335-17-bf1-06.stl",
		"1105335-18-bf2-07.stl",
		"1105335-19-bf3-08.stl",
		"1105335-17-bf1-09.stl",
		"1105335-18-bf2-10.stl",
		"1105333-19-bt1-11.stl",
		"1105334-17-bi1-12.stl",
		"1105333-21-bt3-13.stl",
		"1105334-19-bi1-14.stl",
		"1105333-20-bt2-15.stl",
		"1105334-17-bi1-16.stl",
		"1105333-22-bt4-17.stl",
		"1105334-18-bi2-18.stl",
		"1105333-19-bt1-19.stl",
		"1105334-17-bi1-20.stl",
		"1105333-21-bt3-21.stl",
		"1105334-19-bi3-22.stl",
		"1105333-20-bt2-23.stl",
		"1105334-17-bi1-24.stl",
		"1105333-22-bt4-25.stl",
		"1105334-18-bi2-26.stl",
		"1105333-19-bt1-27.stl",
		"1105334-17-bi1-28.stl",
		"1105333-21-bt3-29.stl",
		"1105334-19-bi3-30.stl",
		"1105333-20-bt2-31.stl",
		"1105334-17-bi1-32.stl",
		"1105333-22-bt4-33.stl",
		"1105334-18-bi2-34.stl",
		"1105333-19-bt1-35.stl",
		"1105334-17-bi1-36.stl",
		"1105333-21-bt3-37.stl",
		"1105334-19-bi3-38.stl",
		"1105333-20-bt2-39.stl",
		"1105334-17-bi1-40.stl",
		"1105333-22-bt4-41.stl",
		"1105334-18-bi2-42.stl",
		"1105333-19-bt1-43.stl",
		"1105334-17-bi1-44.stl",
		"1105333-21-bt3-45.stl",
		"1105334-19-bi3-46.stl",
		"1105333-20-bt2-47.stl",
		"1105334-17-bi1-48.stl",
		"1105333-22-bt4-49.stl",
		"1105334-18-bi2-50.stl",
		"1105335-19-bf3-51.stl",
		"1105335-17-bf1-52.stl",
		"1105335-18-bf2-53.stl",
		"1105335-19-bf3-54.stl",
		"1105335-17-bf1-55.stl",
		"1105335-18-bf2-56.stl",
		"1105335-19-bf3-57.stl",
		"1105335-17-bf1-58.stl",
		"1105335-18-bf2-59.stl",
		"100026959-04-bo60-60.stl",
	};



	auto mlc_volumes = MLCLeafVolumes();
	//Bank X2 / A
	{
		size_t num_leaves = sizeof(filenames_bank_X2) / sizeof(*filenames_bank_X2);
		for (size_t i = 0; i < num_leaves; ++i) {
			std::filesystem::path filepath = stl_folder / filenames_bank_X2[i];
			const G4ThreeVector offset(9.0, 0.0, 0.0);
            G4VSolid* leaf_solid = getTessalatedSolidFromSTL(filenames_bank_X2[i], filepath, offset);
			G4LogicalVolume* leaf_logical = new G4LogicalVolume(leaf_solid, W92_5, "mlc_leaf", 0, 0, 0);
			leaf_logical->SetVisAttributes(bank_X2_vis);
			auto volume = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_world_z),
				leaf_logical, "millennium_leaf", parent_logical, false, 0);
			mlc_volumes.bank_X2.push_back(volume);
		}
	}
	//Bank X1 / B
	{
		size_t num_leaves = sizeof(filenames_bank_X1) / sizeof(*filenames_bank_X1);
		for (size_t i = 0; i < num_leaves; ++i) {
			std::filesystem::path filepath = stl_folder / filenames_bank_X1[i];
			const G4ThreeVector offset(191.0, 0.0, 0.0);
            G4VSolid* leaf_solid = getTessalatedSolidFromSTL(filenames_bank_X1[i], filepath, offset);
			G4LogicalVolume* leaf_logical = new G4LogicalVolume(leaf_solid, W92_5, "mlc_leaf", 0, 0, 0);
			leaf_logical->SetVisAttributes(bank_X1_vis);
			auto volume = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_world_z),
				leaf_logical, "millennium_leaf", parent_logical, false, 0);
			mlc_volumes.bank_X1.push_back(volume);
		}
	}
	return mlc_volumes;
}
void setMilleniumMLCLeafPositions(const MLCLeafVolumes& volumes, const LeafPositions& positions)
{
	// position is positive for non-overtravel, for both leaf banks
	const size_t num_leaves_in_bank = 60;
	const bool valid = (volumes.bank_X2.size() == num_leaves_in_bank) || (volumes.bank_X1.size() == num_leaves_in_bank) 
		|| (positions.bank_X2.size() == num_leaves_in_bank) || (positions.bank_X1.size() == num_leaves_in_bank);
	if (!valid)
		throw std::runtime_error("Number of leaf volumes and leaf positions is inconsistent");
	const std::vector<double>* bank_positions[] = { &positions.bank_X2, &positions.bank_X1 };
	const std::vector<G4VPhysicalVolume*>* bank_volumes[] = { &volumes.bank_X2, &volumes.bank_X1 };
	const G4int num_pos = 41;
	double physical_position = 0.0;
	/*
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
	*/	
	for (size_t j = 0; j < 2; ++j) {
		for (size_t i = 0; i < num_leaves_in_bank; ++i) {
			double nominal_position = bank_positions[j]->at(i);
			if (j == 1) nominal_position *= (-1);
			//if (std::abs(nominal_position) - 10.0*cm < 0)
			physical_position = nominal_position + 80 * mm * (std::sqrt(1 + std::pow(nominal_position / g_SAD * mm, 2)) - 1)* (1000/509.952);
			//else {
			//	for (size_t k = 0; k < num_pos; ++k) {
			//		if (std::abs(nominal_position - mlc_nominal[k]) < 1e-7) {
			//			physical_position = mlc_actual[k];
			//		break;
			//		}else physical_position = (nominal_position - mlc_nominal[k]) / (mlc_nominal[k-1] - mlc_nominal[k]) * (mlc_actual[k-1] - mlc_actual[k]) + mlc_actual[k];
			//	}
			//}
			if (j == 1) physical_position *= (-1);
			physical_position *= 509.952/1000.;
			auto volume = bank_volumes[j]->at(i);
			auto translation = volume->GetTranslation();
			translation.setX(physical_position);
			volume->SetTranslation(translation);
		}
	}
}


void buildTarget(G4LogicalVolume* parent_logical, double parent_z_world, Target target_type)
{

	const double target_z_world = g_SAD;
	auto target_region = new G4Region("target");
	auto target_vis_attributes = new G4VisAttributes(G4Colour(0.2, 0.2, 0.8, 0.9));
	G4VisAttributes* target_block_vis_attributes = new G4VisAttributes(G4Colour(0.5, 0.0, 0.5, 0.5));


	//Materials for target
	{
		G4NistManager* NISTman = G4NistManager::Instance();
		G4Element* elAl = NISTman->FindOrBuildElement("Al");
		G4Element* elAu = NISTman->FindOrBuildElement("Au");
		G4Element* elCu = NISTman->FindOrBuildElement("Cu");
		G4Element* elNi = NISTman->FindOrBuildElement("Ni");
		G4Element* elO = NISTman->FindOrBuildElement("O");

		{
			// nicoro
			const double density = 10.9 * g / cm3;
			G4Material* nicoro = new G4Material("Nicoro", density, 3);
			nicoro->AddElement(elAu, 35. * perCent);
			nicoro->AddElement(elCu, 62. * perCent);
			nicoro->AddElement(elNi, 3. * perCent);
		}
		{
			// copper glidcop
			const double density = 8.88528 * g / cm3;
			G4Material* copper_glidcop = new G4Material("COPPER_GLIDCOP", density, 3);
			copper_glidcop->AddElement(elCu, 99.600 * perCent);
			copper_glidcop->AddElement(elO, 0.1883 * perCent);
			copper_glidcop->AddElement(elAl, 0.2117 * perCent);
		}
	}


	if (target_type == Target::LOW_ENERGY) {
		// Be window
		G4double BeThick6X = 0.0254 * cm;  // same as BeThick for orbit chamber
		G4Tubs* BeWindow6X = new G4Tubs("BeWindow6X", 0. * mm, 6.477 * mm, BeThick6X / 2., 0. * deg, 360. * deg);
		auto be_window_logical = new G4LogicalVolume(BeWindow6X, 
			G4Material::GetMaterial("G4_Be"), "BeWindow6X_LV", 0, 0, 0);
		be_window_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.5, 0.5, 1.0, 0.5)));

		new G4PVPlacement(0, G4ThreeVector(0., 0., target_z_world + 0.751 * mm - parent_z_world),
			be_window_logical, "BeWindowLow", parent_logical, false, 0);

		// Vacuum Chamber between Be window and target button
		const G4int vc_points = 6;
		G4double r_vc[vc_points] =
		{ 0. * mm,  3.785 * mm,  3.785 * mm,  4.9 * mm,  4.9 * mm,  0. * mm };
		G4double z_vc[vc_points] =
		{ 0.261 * mm,  0.261 * mm,  0.324 * mm,  0.324 * mm,  0.624 * mm,  0.624 * mm };

		G4Polycone* vacuum_chamber = new G4Polycone("vacuum_chamber_polycone",
			0. * deg, 360. * deg, vc_points, r_vc, z_vc);
		auto vacuum_chamber_logical = new G4LogicalVolume(vacuum_chamber,
				G4Material::GetMaterial("G4_Galactic"), "vacuum_chamber", 0, 0, 0);
		vacuum_chamber_logical->SetVisAttributes(new G4VisAttributes(G4Colour(1., 0., 0., 0.6)));

		new G4PVPlacement(0,
			G4ThreeVector(0., 0., target_z_world - parent_z_world), vacuum_chamber_logical,
			"vacuum_chamber", parent_logical, false, 0);

		// 6X Target Button
		G4double nicoroHThick = 0.0508 * mm; // half the nicoro thickness

		const G4int tb6pts = 5;
		G4double r_tb[tb6pts] =
		{ 0.0 * mm, 3.7592 * mm, 3.7592 * mm,  3.5052 * mm,  0.0 * mm };
		G4double z_tb[tb6pts] =
		{ 0.0 * mm,    0.0 * mm, -0.381 * mm, -0.635 * mm, -0.635 * mm };

		G4double Target_z = -0.159 * mm;

		G4Polycone* target_button = new G4Polycone("target_button", 0. * deg, 360. * deg, tb6pts, r_tb, z_tb);

		auto target_logical =
			new G4LogicalVolume(target_button, G4Material::GetMaterial("G4_W"),
				"target_button", 0, 0, 0);
		target_logical->SetVisAttributes(target_vis_attributes);

		new G4PVPlacement(0,
			G4ThreeVector(0., 0., target_z_world - Target_z + 2. * nicoroHThick - parent_z_world),
			target_logical, "target_button", parent_logical,
			false, 0);

		// nicoro brazing sheet
		// not sure what the material should be: take it to be nicoro (BAu-3) 
		// (other choice is nicoro 80)
		G4Tubs* target_nicoro =
			new G4Tubs("nicoro", 0.0 * mm, 3.0734 * mm, nicoroHThick, 0. * deg, 360. * deg);

		auto target_nicoro_logical =
			new G4LogicalVolume(target_nicoro, G4Material::GetMaterial("Nicoro"),
				"target_nicoro", 0, 0, 0);
		target_nicoro_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.0, 0.5, 0.5)));


		new G4PVPlacement(0,
			G4ThreeVector(0., 0., target_z_world - Target_z - 0.635 * mm + nicoroHThick - parent_z_world),
			target_nicoro_logical, "target_nicoro", parent_logical, false, 0);

		//  Target block top
		const G4int tbt_pts = 31;
		G4double r_pc[tbt_pts] = {
			0.0 * mm, 3.785 * mm, 3.785 * mm, 4.9 * mm,
			4.9 * mm, 6.6 * mm, 6.6 * mm, 8.2 * mm,
			8.2 * mm, 7.1 * mm, 7.1 * mm, 6.6 * mm,
			6.6 * mm, 4.502 * mm,
			3.2137 * mm, 3.1928 * mm, 3.1660 * mm, 3.1338 * mm,
			3.0963 * mm, 3.0538 * mm, 3.0067 * mm, 2.9552 * mm,
			2.8999 * mm, 2.8410 * mm, 2.7791 * mm, 2.7146 * mm,
			2.6481 * mm, 2.5799 * mm, 2.5107 * mm, 2.4410 * mm, 0.0 * mm };

		G4double z_pc[tbt_pts] = {
			-0.476 * mm, -0.476 * mm,  0.324 * mm,  0.324 * mm,
			 0.624 * mm,  0.624 * mm,  1.124 * mm,  1.124 * mm,
			-5.876 * mm, -5.876 * mm, -8.876 * mm, -8.876 * mm,
			-8.376 * mm, -8.376 * mm,
			-3.5689 * mm, -3.5024 * mm, -3.4379 * mm, -3.3760 * mm,
			-3.3171 * mm, -3.2618 * mm, -3.2103 * mm, -3.1632 * mm,
			-3.1207 * mm, -3.0832 * mm, -3.0510 * mm, -3.0242 * mm,
			-3.0033 * mm, -2.9882 * mm, -2.9790 * mm, -2.9760 * mm, -2.9760 * mm };

		G4GenericPolycone* target_block_top =
			new G4GenericPolycone("target_block_top", 0. * deg, 360. * deg, tbt_pts, r_pc, z_pc);

		auto target_block_top_logical =
			new G4LogicalVolume(target_block_top, G4Material::GetMaterial("COPPER_GLIDCOP"),
				"target_block_top", 0, 0, 0);
		target_block_top_logical->SetVisAttributes(target_block_vis_attributes);

		new G4PVPlacement(0, G4ThreeVector(0., 0., target_z_world - parent_z_world),
				target_block_top_logical, "target_block_top", parent_logical, false, 0);

		// X-ray window
		G4Tubs* xray_window_tubs = new G4Tubs("xray_window", 0.0 * mm, 6.4135 * mm, 0.0508 * mm, 0. * deg, 360. * deg);
		auto xray_window_logical =
			new G4LogicalVolume(xray_window_tubs, G4Material::GetMaterial("SS304"), "xray_window", 0, 0, 0);
		xray_window_logical->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 0.0, 1.0)));

		new G4PVPlacement(0,
			G4ThreeVector(0. * m, 0. * m, target_z_world - 8.503 * mm - parent_z_world),
			xray_window_logical, "xray_window", parent_logical, false, 0);
		be_window_logical->SetRegion(target_region);
		target_logical->SetRegion(target_region);
		target_nicoro_logical->SetRegion(target_region);
		target_block_top_logical->SetRegion(target_region);
		target_region->AddRootLogicalVolume(target_logical);
		target_region->AddRootLogicalVolume(target_nicoro_logical);
		target_region->AddRootLogicalVolume(target_block_top_logical);
		target_region->AddRootLogicalVolume(be_window_logical);
	}
	else {
		throw std::runtime_error("Target type not supported in buildTarget");
	}
}

//Returns the logical volume where the dose should be recorded
G4LogicalVolume* buildIonisationChamber(G4LogicalVolume* parent_logical, double parent_z_world)
{
	const double ic_position_z_world = g_SAD - 156.9835 * mm;
	// ic_position_z_world is the centre of entire IC
	// inner half thickness of ic body
	const G4double icbody_inthick = 8.31 * mm;
	// half-thickness of .005" kapton
	const G4double ICWin_5thou = 0.0635 * mm;
	// halfthickness of .002" kapton
	const G4double IC_electrodeHthick = 0.0254 * mm;
	// half thickness of center holder
	const G4double IC_centralSupportHthick = 0.635 * mm;
	// half thickness of center ring
	const G4double IC_centerringHthick = 0.508 * mm;
	// half thickness of bulk of spacer
	const G4double IC_topringHthick = 1.10 * mm;
	const G4double IC_SE_pos = IC_centralSupportHthick + 2. * IC_topringHthick +
		IC_electrodeHthick;
	const G4double IC_HV_pos = IC_centralSupportHthick + 2. * IC_electrodeHthick +
		2. * IC_topringHthick + 2. * IC_centerringHthick +
		IC_electrodeHthick;

	G4VisAttributes* VisAtt_IonChamber1 =
		new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.5));
	G4VisAttributes* VisAtt_IonChamber2 =
		new G4VisAttributes(G4Colour(0.0, 0.2, 1.0, 0.8));
	G4VisAttributes* VisAtt_IonChamber3 =
		new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.6));
	G4VisAttributes* VisAtt_IonChamber4 =
		new G4VisAttributes(G4Colour(0.2, 1.0, 0.0, 0.4));

	const double cu_thickness_factor = 1.0;
	const double ic_cu_thick = 0.0025 * mm;

	// first create a mother volume to put the bits in
	// mother volume is used for the sensitive detector
	// TODO: air in ion chamber is at 1.2 bar
	G4Tubs* ic_mother_tubs = new G4Tubs("ic_mother", 0., 80. * mm, 16.0 * mm,
		0. * deg, 360. * deg);
	auto ic_mother_logical = new G4LogicalVolume(ic_mother_tubs,
		G4Material::GetMaterial("G4_AIR"), "ic_mother", 0, 0, 0);
	ic_mother_logical->SetVisAttributes(new G4VisAttributes(false, G4Colour(0.2, 0.4, 0.7, 0.8)));
	
	{
		G4RotationMatrix* rot = new G4RotationMatrix();
		rot->rotateX(180.0 * deg);
		new G4PVPlacement(rot,
			G4ThreeVector(0. * m, 0. * m, ic_position_z_world - parent_z_world),
			ic_mother_logical, "ic_mother", parent_logical, false, 0);
	}
	{
		// Ion Chamber body 
		const G4int icbody_pts = 6;
		G4double rInner_icbody[icbody_pts] =
		{ 50.16 * mm, 49.02 * mm, 66.68 * mm, 66.68 * mm, 49.02 * mm, 50.16 * mm };
		G4double rOuter_icbody[icbody_pts] =
		{ 76.2 * mm, 76.2 * mm, 76.2 * mm, 76.2 * mm, 76.2 * mm, 76.2 * mm };
		G4double zPlane_icbody[icbody_pts] =
		{ -15.09 * mm, -icbody_inthick, -icbody_inthick,
		 icbody_inthick,  icbody_inthick, 15.09 * mm };

		G4Polycone* IC_body_PV = new G4Polycone("ic_body", 0. * deg,
			360. * deg, icbody_pts, zPlane_icbody, rInner_icbody, rOuter_icbody);

		G4LogicalVolume* ic_body_logical =
			new G4LogicalVolume(IC_body_PV, G4Material::GetMaterial("SS304"), "ic_body", 0, 0, 0);
		
		new G4PVPlacement(0, G4ThreeVector(), ic_body_logical, "ic_body", ic_mother_logical, false, 0);
		ic_body_logical->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.5)));
	}
	{
		// Ion Chamber Window - 3 instances
		//  Cu plating
		//  put it on downstream side for all, for simplicity

		auto fICWinCu_tubs = new G4Tubs("ic_window_cu", 0. * mm, 48. * mm,
			ic_cu_thick * cu_thickness_factor, 0. * deg, 360. * deg);
		auto ic_window_cu_logical = new G4LogicalVolume(fICWinCu_tubs,
			G4Material::GetMaterial("G4_Cu"), "ic_window_cu", 0, 0, 0);
		ic_window_cu_logical->SetVisAttributes(VisAtt_IonChamber3);
		

		G4Tubs* ICWin = new G4Tubs("ic_window", 0. * mm, 54.61 * mm,
			ICWin_5thou, 0 * deg, 360 * deg);
		G4LogicalVolume* ic_window = new G4LogicalVolume(ICWin,
			G4Material::GetMaterial("kapton"), "ic_window", 0, 0, 0);
		ic_window->SetVisAttributes(VisAtt_IonChamber2);

		//first window
		new G4PVPlacement(0, G4ThreeVector(0. * m, 0. * m, -icbody_inthick + ICWin_5thou),
			ic_window, "ion_chamber_window", ic_mother_logical, false, 0);

		new G4PVPlacement(0, G4ThreeVector(0. * m, 0. * m, -icbody_inthick +
			2. * ICWin_5thou + ic_cu_thick * cu_thickness_factor),
			ic_window_cu_logical, "ic_window_cu", ic_mother_logical, false, 0);

		//second window
		new G4PVPlacement(0, G4ThreeVector(),
			ic_window, "ion_chamber_window", ic_mother_logical, false, 1);

		new G4PVPlacement(0, G4ThreeVector(0. * m, 0. * m, ICWin_5thou +
			ic_cu_thick * cu_thickness_factor),
			ic_window_cu_logical, "ic_window_cu", ic_mother_logical, false, 1);

		//third window
		new G4PVPlacement(0, G4ThreeVector(0. * m, 0. * m, icbody_inthick - ICWin_5thou),
			ic_window, "ion_chamber_window", ic_mother_logical, false, 2);

		new G4PVPlacement(0, G4ThreeVector(0. * m, 0. * m, icbody_inthick
			- 2. * ICWin_5thou - ic_cu_thick * cu_thickness_factor),
			ic_window_cu_logical, "ic_window_cu", ic_mother_logical, false, 2);
	}
	{
		// electrodes 
		// spacer rings 

		// top and bottom rings (Al 2024)
		G4Tubs* IC_topring = new G4Tubs("ic_topring", 50.42 * mm, 58.42 * mm,
			IC_topringHthick, 0. * deg, 360. * deg);
		G4LogicalVolume* ic_topring_logical = new G4LogicalVolume(IC_topring,
			G4Material::GetMaterial("Aluminum2024"), "ic_topring", 0, 0, 0);

		ic_topring_logical->SetVisAttributes(VisAtt_IonChamber4);

		new G4PVPlacement(0, G4ThreeVector(0. * m, 0. * m, -IC_topringHthick -
			IC_centralSupportHthick),
			ic_topring_logical, "top_or_bottom_spacer", ic_mother_logical, false, 0);

		new G4PVPlacement(0, G4ThreeVector(0. * m, 0. * m, +IC_topringHthick +
			IC_centralSupportHthick),
			ic_topring_logical, "top_or_bottom_spacer", ic_mother_logical, false, 1);

		new G4PVPlacement(0, G4ThreeVector(0. * m, 0. * m, -IC_HV_pos -
			IC_electrodeHthick - IC_topringHthick),
			ic_topring_logical, "top_or_bottom_spacer", ic_mother_logical, false, 2);

		new G4PVPlacement(0, G4ThreeVector(0. * m, 0. * m, +IC_HV_pos +
			IC_electrodeHthick + IC_topringHthick),
			ic_topring_logical, "top_or_bottom_spacer", ic_mother_logical, false, 3);

		// center ring
		G4Tubs* ic_center_ring = new G4Tubs("ic_center_ring", 50.42 * mm, 58.42 * mm,
			IC_centerringHthick, 0. * deg, 360. * deg);
		G4LogicalVolume* ic_center_ring_logical = new G4LogicalVolume(ic_center_ring,
			G4Material::GetMaterial("G4_Cu"), "IC_center_ring", 0, 0, 0);

		ic_center_ring_logical->SetVisAttributes(VisAtt_IonChamber4);
		

		new G4PVPlacement(0, G4ThreeVector(0. * m, 0. * m, 
			+(IC_centralSupportHthick + 2. * IC_electrodeHthick + 2. * IC_topringHthick + IC_centerringHthick)),
			ic_center_ring_logical, "center_spacer", ic_mother_logical, false, 0);

		new G4PVPlacement(0, G4ThreeVector(0. * m, 0. * m, 
			-(IC_centralSupportHthick + 2. * IC_electrodeHthick + 2. * IC_topringHthick + IC_centerringHthick)),
			ic_center_ring_logical, "center_spacer", ic_mother_logical, false, 1);
	}
	{
		// HV Electrode 
		// kapton
		G4Tubs* ic_hv = new G4Tubs("ic_electrode_hv", 0. * mm, 54.61 * mm,
			IC_electrodeHthick, 0. * deg, 360. * deg);
		G4LogicalVolume* ic_hv_logical = new G4LogicalVolume(ic_hv,
			G4Material::GetMaterial("kapton"), "ic_electrode_hv", 0, 0, 0);

		ic_hv_logical->SetVisAttributes(VisAtt_IonChamber2);
		
		{
			//  Cu plating
			auto ic_cu_plating_tubs = new G4Tubs("ic_cu_plating", 0. * mm, 48.0 * mm,
				ic_cu_thick * cu_thickness_factor, 0. * deg, 360. * deg);
			auto ic_cu_plating_logical = new G4LogicalVolume(ic_cu_plating_tubs,
				G4Material::GetMaterial("G4_Cu"), "ic_cu_plating", 0, 0, 0);
			ic_cu_plating_logical->SetVisAttributes(VisAtt_IonChamber3);
			

			new G4PVPlacement(0,
				G4ThreeVector(0. * m, 0. * m, -IC_HV_pos),
				ic_hv_logical, "ic_electrode_hv", ic_mother_logical, false, 0);

			new G4PVPlacement(0,
				G4ThreeVector(0. * m, 0. * m, -IC_HV_pos + IC_electrodeHthick + ic_cu_thick * cu_thickness_factor),
				ic_cu_plating_logical, "ic_cu_plating", ic_mother_logical, false, 0);

			new G4PVPlacement(0, G4ThreeVector(0. * m, 0. * m, IC_HV_pos),
				ic_hv_logical, "ic_electrode_hv", ic_mother_logical, false, 1);

			new G4PVPlacement(0,
				G4ThreeVector(0. * m, 0. * m, IC_HV_pos + IC_electrodeHthick + ic_cu_thick * cu_thickness_factor),
				ic_cu_plating_logical, "ic_cu_plating", ic_mother_logical, false, 1);
		}
		{
			// Signal Electrode 

			//  Cu plating, ignoring the patterning
			auto ic_signal_electrode_tubs = new G4Tubs("ic_cu_sig", 0. * mm, 48. * mm,
				ic_cu_thick * cu_thickness_factor, 0. * deg, 360. * deg);
			auto ic_signal_electrode_logical = new G4LogicalVolume(ic_signal_electrode_tubs,
				G4Material::GetMaterial("G4_Cu"), "ic_cu_sig", 0, 0, 0);
			ic_signal_electrode_logical->SetVisAttributes(VisAtt_IonChamber3);

			// reuse the HV electrode LV for the kapton
			new G4PVPlacement(0, G4ThreeVector(0. * m, 0. * m, -IC_SE_pos),
				ic_hv_logical, "ic_signal_electrode", ic_mother_logical, false, 0);

			new G4PVPlacement(0, G4ThreeVector(0. * m, 0. * m, -IC_SE_pos +
				IC_electrodeHthick + ic_cu_thick * cu_thickness_factor),
				ic_signal_electrode_logical, "ic_cu_sig", ic_mother_logical, false, 0);

			new G4PVPlacement(0, G4ThreeVector(0. * m, 0. * m, IC_SE_pos),
				ic_hv_logical, "ic_signal_electrode", ic_mother_logical, false, 1);

			new G4PVPlacement(0, G4ThreeVector(0. * m, 0. * m, IC_SE_pos +
				IC_electrodeHthick + ic_cu_thick * cu_thickness_factor),
				ic_signal_electrode_logical, "ic_cu_sig", ic_mother_logical, false, 1);
		}
	}
	// ring for central kapton 
	{
		const G4int IC_WinRing_pts = 4;
		G4double rInner_icr[IC_WinRing_pts] =
		{ 55.37 * mm, 55.37 * mm, 48.89 * mm, 48.89 * mm };
		//rOuter lower than nominal; extra is in ic_body
		G4double rOuter_icr[IC_WinRing_pts] =
		{ 66.67 * mm, 66.67 * mm, 66.67 * mm, 66.67 * mm };
		G4double zPlane_icr[IC_WinRing_pts] =
		{ -.635 * mm,   .076 * mm,   .076 * mm,   .635 * mm };

		G4Polycone* IC_WinRing_PV = new G4Polycone("IC_windowring", 0. * deg,
			360. * deg, IC_WinRing_pts, zPlane_icr, rInner_icr, rOuter_icr);

		G4LogicalVolume* IC_WinRing_LV =
			new G4LogicalVolume(IC_WinRing_PV, G4Material::GetMaterial("G4_Cu"),
				"IC_windowring_LV", 0, 0, 0);

		IC_WinRing_LV->SetVisAttributes(VisAtt_IonChamber4);
		
		new G4PVPlacement(0, G4ThreeVector(),
			IC_WinRing_LV, "ic_winring", ic_mother_logical, false, 0);

		// volume for sensitive detector 

		G4Tubs* dose_tubs =
			new G4Tubs("ic_dose", 0. * mm, 48. * mm, 0.4 * mm, 0. * deg, 360. * deg);

		auto ic_dose_logical =
			new G4LogicalVolume(dose_tubs, G4Material::GetMaterial("G4_AIR"),
				"ic_dose", 0, 0, 0);
		ic_dose_logical->SetVisAttributes(VisAtt_IonChamber2);
		

			new G4PVPlacement(0,
				G4ThreeVector(0., 0., IC_SE_pos + (IC_HV_pos - IC_SE_pos) / 2),
				ic_dose_logical, "ic_dose", ic_mother_logical, false, 0);

			return ic_dose_logical;
	}
}

void buildShieldingCollimator(G4LogicalVolume* parent, double parent_z_world, const std::filesystem::path& stl_folder)
{
	// Shielding collimator
	//Visualization
	const auto shield_vis = new G4VisAttributes(G4Colour(0.0, 0.8, 0.8, 0.4));
	
	static const char* filename_shieldcoll = "100022238-68-10.stl";
	std::filesystem::path filepath = stl_folder / filename_shieldcoll;
	const G4ThreeVector offset(0.0, 0.0, 0.0);
    G4VSolid* shieldcoll_solid = getTessalatedSolidFromSTL(filename_shieldcoll, filepath, offset);
	G4LogicalVolume *shieldColl_LV = new G4LogicalVolume(shieldcoll_solid, G4Material::GetMaterial("Lead97Antimony"), "shieldColl_LV", 0, 0, 0);
	shieldColl_LV->SetVisAttributes(shield_vis);
	auto ShieldColl = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_z_world), shieldColl_LV, "shieldColl", parent, false, 0);
}


void buildMylarWindow(G4LogicalVolume* parent, double parent_z_world)
{
	const double mylar_window_z_world = g_SAD - 557 * mm;
	G4NistManager* NISTman = G4NistManager::Instance();
	G4Box* mylar_window = new G4Box("mylar_window", 140. * mm, 140. * mm, 0.05 * mm);
	auto mylar_window_logical =
		new G4LogicalVolume(mylar_window, NISTman->FindOrBuildMaterial("G4_MYLAR"), "mylar_window", 0, 0, 0);
	mylar_window_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.5, 0.5, 1.0, 0.5)));

	new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, mylar_window_z_world - parent_z_world), 
		mylar_window_logical, "mylar_window", parent, false, 0);
}

//Build common materials used in detector
void buildMaterials() {
	G4NistManager* NISTman = G4NistManager::Instance();

	NISTman->FindOrBuildMaterial("G4_Galactic");
	NISTman->FindOrBuildMaterial("G4_W");
	NISTman->FindOrBuildMaterial("G4_Be");
	NISTman->FindOrBuildMaterial("G4_Cu");

	G4Element* elH = NISTman->FindOrBuildElement("H");
	G4Element* elO = NISTman->FindOrBuildElement("O");
	G4Element* elAl = NISTman->FindOrBuildElement("Al");
	G4Element* elW = NISTman->FindOrBuildElement("W");
	G4Element* elNi = NISTman->FindOrBuildElement("Ni");
	G4Element* elFe = NISTman->FindOrBuildElement("Fe");
	G4Element* elS = NISTman->FindOrBuildElement("S");
	G4Element* elSi = NISTman->FindOrBuildElement("Si");
	G4Element* elC = NISTman->FindOrBuildElement("C");
	G4Element* elP = NISTman->FindOrBuildElement("P");
	G4Element* elN = NISTman->FindOrBuildElement("N");
	G4Element* elCr = NISTman->FindOrBuildElement("Cr");
	G4Element* elMn = NISTman->FindOrBuildElement("Mn");
	G4Element* elMg = NISTman->FindOrBuildElement("Mg");
	G4Element* elCu = NISTman->FindOrBuildElement("Cu");
	G4Element* elPb = NISTman->FindOrBuildElement("Pb");
	G4Element* elSb = NISTman->FindOrBuildElement("Sb");


	{
		const double density = 18.0 * g / cm3;
		G4Material* W95 = new G4Material("W95", density, 3);
		W95->AddElement(elW, 95.0 * perCent);
		W95->AddElement(elNi, 3.5 * perCent);
		W95->AddElement(elFe, 1.5 * perCent);
	}
	{
		const double density = 8.0 * g / cm3;
		G4Material* SS304 = new G4Material("SS304", density, 9);
		SS304->AddElement(elS, 0.03 * perCent);
		SS304->AddElement(elC, 0.08 * perCent);
		SS304->AddElement(elSi, 0.75 * perCent);
		SS304->AddElement(elP, 0.045 * perCent);
		SS304->AddElement(elN, 0.10 * perCent);
		SS304->AddElement(elCr, 19.00 * perCent);
		SS304->AddElement(elMn, 2.00 * perCent);
		SS304->AddElement(elNi, 9.25 * perCent);
		SS304->AddElement(elFe, 68.745 * perCent);
	}
	{
		const double density = 7.86112 * g / cm3;
		G4Material* SS_A36 = new G4Material("SS_A36", density, 5);
		SS_A36->AddElement(elC, 0.25 * perCent);
		SS_A36->AddElement(elFe, 99.26 * perCent);
		SS_A36->AddElement(elSi, 0.4 * perCent);
		SS_A36->AddElement(elP, 0.04 * perCent);
		SS_A36->AddElement(elS, 0.05 * perCent);
	}
	{
		const double density = 1.42 * g / cm3;
		G4Material* kapton = new G4Material("kapton", density, 4);
		kapton->AddElement(elH, 02.6362 * perCent);
		kapton->AddElement(elC, 69.1133 * perCent);
		kapton->AddElement(elN, 07.3270 * perCent);
		kapton->AddElement(elO, 20.9235 * perCent);
	}
	{
		//aluminum 2024   //from wikipedia
		const double density = 2.78 * g / cm3;
		G4Material* Aluminum2024 = new G4Material("Aluminum2024", density, 4);
		Aluminum2024->AddElement(elAl, 93.7 * perCent);
		Aluminum2024->AddElement(elMg, 1.4 * perCent);
		Aluminum2024->AddElement(elMn, 0.5 * perCent);
		Aluminum2024->AddElement(elCu, 4.4 * perCent);
	}
	{
		const double density = 11.0 * g / cm3;
		G4Material* Lead97Antimony = new G4Material("Lead97Antimony", density, 2);
		Lead97Antimony->AddElement(elPb, 97. * perCent);
		Lead97Antimony->AddElement(elSb, 3. * perCent);
	}
	{
		// mixture of stuff for misc. electronics etc.
		const double density = 2. * g / cm3;
		G4Material* MiscMaterial = new G4Material("misc", density, 5);
		MiscMaterial->AddElement(elAl, 20. * perCent);
		MiscMaterial->AddElement(elFe, 20. * perCent);
		MiscMaterial->AddElement(elCu, 20. * perCent);
		MiscMaterial->AddElement(elC, 20. * perCent);
		MiscMaterial->AddElement(elO, 20. * perCent);
	}
}

//flattening_filter_z_world = g_SAD - 125.7 * mm
void buildFlatteningFilter(G4LogicalVolume* parent_logical, double parent_z_world, 
		FlatteningFilter flattening_filter)
{
	//  Aluminum carousel (holds flattening filters)
	const double flattening_filter_z_world = g_SAD - 125.7 * mm;
	{
		G4NistManager* NISTman = G4NistManager::Instance();
		G4Element* elAl = NISTman->FindOrBuildElement("Al");
		G4Element* elSi = NISTman->FindOrBuildElement("Si");
		G4Element* elMg = NISTman->FindOrBuildElement("Mg");
		G4Element* elCu = NISTman->FindOrBuildElement("Cu");
		G4Element* elCr = NISTman->FindOrBuildElement("Cr");

		const double density = 2.70 * g / cm3;
		G4Material* Aluminum6061 =
			new G4Material("Aluminum6061", density, 5);
		Aluminum6061->AddElement(elAl, 98.01 * perCent);
		Aluminum6061->AddElement(elSi, 0.6 * perCent);
		Aluminum6061->AddElement(elMg, 1.2 * perCent);
		Aluminum6061->AddElement(elCu, 0.15 * perCent);
		Aluminum6061->AddElement(elCr, 0.04 * perCent);
	}

	const G4int npts_carousel = 6;
	G4double rInner_carousel[npts_carousel] =
	{ 34. * mm, 34. * mm, 38.1 * mm, 38.1 * mm, 44. * mm, 44. * mm };
	G4double rOuter_carousel[npts_carousel] =
	{ 200. * mm, 200. * mm, 200. * mm, 200. * mm, 200. * mm, 200. * mm };
	G4double zPlane_carousel[npts_carousel] =
	{ 8.9 * mm, 0. * mm, 0. * mm, -1.5 * mm, -1.5 * mm, -7. * mm };

	G4Polycone* carousel =
		new G4Polycone("ff_carousel", 0. * deg, 360. * deg, npts_carousel,
			zPlane_carousel, rInner_carousel, rOuter_carousel);

	G4LogicalVolume* carousel_logical =
		new G4LogicalVolume(carousel,
			G4Material::GetMaterial("Aluminum6061"),
			"ff_carousel", 0, 0, 0);

	//In previous version it was possible to tune flattening filter offset from default position
	auto ff_offset = G4ThreeVector(0.0, 0.0, 0.0);
	new G4PVPlacement(0,
		G4ThreeVector(ff_offset.x(), ff_offset.y(), flattening_filter_z_world + ff_offset.z() - parent_z_world),
		carousel_logical, "ff_carousel", parent_logical, false, 0);
	carousel_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.5, 0.4, 0.0, 0.2)));

	//  Flattening filters
	//G4GenericPolycone* FF_polycone;
	//G4LogicalVolume* FF_LV;

	//G4VisAttributes* VisAtt_FF = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.5));

	if (flattening_filter == FlatteningFilter::OPEN) {
		{
			G4NistManager* NISTman = G4NistManager::Instance();
			G4Element* elZn = NISTman->FindOrBuildElement("Zn");
			G4Element* elPb = NISTman->FindOrBuildElement("Pb");
			G4Element* elFe = NISTman->FindOrBuildElement("Fe");
			G4Element* elCu = NISTman->FindOrBuildElement("Cu");

			const double density = 8.53 * g / cm3;
			G4Material* Brass = new G4Material("Brass", density, 4);
			Brass->AddElement(elZn, 29.88 * perCent);
			Brass->AddElement(elPb, 0.07 * perCent);
			Brass->AddElement(elFe, 0.05 * perCent);
			Brass->AddElement(elCu, 70. * perCent);
		}

		const double offset_from_carousel_z = -2.31 * mm;
		G4double openPortHThick = 0.405 * mm;
		G4Tubs* openPort =
			new G4Tubs("open_port", 0. * mm, 42.2 * mm, openPortHThick, 0. * deg, 360. * deg);

		G4LogicalVolume* open_port_logical =
			new G4LogicalVolume(openPort, G4Material::GetMaterial("Brass"),
				"open_port", 0, 0, 0);
		open_port_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.1, 1.0, 0.1)));

		new G4PVPlacement(0,
			G4ThreeVector(ff_offset.x(), ff_offset.y(),
				flattening_filter_z_world + offset_from_carousel_z + ff_offset.z() - parent_z_world),
			open_port_logical, "openPort", parent_logical, false, 0);
	} 
	else if (flattening_filter == FlatteningFilter::ENERGY_6X) {
			
		const G4int cone1size = 24;
		G4RotationMatrix* FF_rot = new G4RotationMatrix();
		FF_rot->rotateX(180.0*deg);
		//Values from drawings
		G4double r_FF[cone1size] = {
		0.0*mm,  0.635*mm,  1.27 *mm,  1.905*mm,
		2.54 *mm,  3.81 *mm,  5.08 *mm,  6.35 *mm,
		7.62 *mm,  8.89 *mm, 10.16 *mm, 12.7  *mm,
		15.24 *mm, 17.78 *mm, 20.32 *mm, 22.86 *mm,
		25.4  *mm, 27.94 *mm, 30.607*mm, 33.02 *mm,
		33.655*mm, 38.1  *mm, 38.1  *mm,  0.*mm
		};
		G4double z_FF[cone1size] = {
		-18.999*mm,-18.72 *mm,-18.44 *mm,-18.059*mm,
		-17.653*mm,-16.916*mm,-15.57 *mm,-14.478*mm,
		-13.386*mm,-12.268*mm,-11.227*mm, -9.144*mm,
		-7.239*mm, -5.385*mm, -3.708*mm, -2.159*mm,
		-0.737*mm,  0.559*mm,  1.016*mm,  1.016*mm,
		0.*mm,  0.*mm,  3.175*mm,  3.175*mm
		};
		G4GenericPolycone* ff_6X = new G4GenericPolycone("ff_6X", 0.*deg, 360.*deg, cone1size, r_FF, z_FF);
		G4LogicalVolume* FF_6X_logical = new G4LogicalVolume(ff_6X, G4Material::GetMaterial("G4_Cu"), "ff_6X",0,0,0);
		FF_6X_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.1, 1.0, 0.1)));
		new G4PVPlacement(FF_rot, G4ThreeVector(ff_offset.x(), ff_offset.y(), flattening_filter_z_world + ff_offset.z() - parent_z_world),
							FF_6X_logical, "ff_6X", parent_logical, false, 0);
  }
	else {
		throw std::runtime_error("Given flattening filter not supported in buildFlatteningFilter");
	}
}

//base_plate_bottom_z_world is location of base plate bottom in world coordinates
//parent_z_world is the location of parent volume center in world coordinates
void buildBasePlate(G4LogicalVolume* parent_logical, double parent_z_world, fs::path gdml_folder)
{
	//Three parts in the BasePlate construction
	//	* Base Plate (from GDML)
	//	* Base Plates 01
	//  * Base Plate 02 (from GDML)
	//	* MLC Corner Clippers (from GDML)

	const double base_plate_bottom_z_world = 533.02 * mm;

	G4double BasePlate_yRef = 1.9493 * mm;       // from GDML file
	G4double BasePlate_zRef = -4.6241 * mm;       // from GDML file
	//  top:
	G4double BasePlate_z = base_plate_bottom_z_world + 15.24 * mm; // 100034445-2 + gdml
	// in z the central surface is 0.020" above center (this is not in data package!)
	G4double BasePlate_center = 0.508 * mm;

	{
		G4GDMLParser* parser = new G4GDMLParser();
		auto base_plate_file = gdml_folder / "BasePlate.gdml";
		parser->Read(base_plate_file.generic_u8string(), false);
		auto base_plate_logical = parser->GetVolume("BasePlate_LogVol");
		base_plate_logical->SetName("base_plate_gdml");
		base_plate_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.6, 0.0, 0.2, 0.4)));
		delete parser;

		new G4PVPlacement(0, G4ThreeVector(0., BasePlate_yRef, BasePlate_z + BasePlate_zRef - parent_z_world),
			base_plate_logical, "base_plate_gdml", parent_logical, false, 0);
	}

	{
		// Base Plate 01 "tungsten shielding corner" -- fig12 of MCDataPackage
		G4double BasePlate1_zRef = 3.8100*mm;       // from Fastrad GDML file
		
		//coords of mounting holes in baseplate
		G4double bp_x = 134.62*mm; // 5.3", from MC data package
		G4double bp_y = 88.90*mm; // 3.5", from MC data package

		//coords of mounting holes in shield (long side), rel. to center
		G4double holepos_x = 22.86*mm; // 0.9" , from MC data package
		G4double holepos_y = 8.89*mm; // 0.35", from MC data package
		
		G4double bp_x1 = 45.72*mm;
		G4double bp_y1 = 16.51*mm;
		G4double bp_chamfer = 5.08*mm;  // the small chamfers
		G4double bp_corner = 30.861*mm; // the large diagonal cut
		G4TwoVector v11(bp_x1 - bp_corner, bp_y1);
		G4TwoVector v21(bp_x1, +bp_y1 - bp_corner);
		G4TwoVector v22(bp_x1, -bp_y1);
		G4TwoVector v31(-bp_x1 + bp_chamfer, -bp_y1);  // these next 4 points are the
		G4TwoVector v32(-bp_x1, -bp_y1 + bp_chamfer);  // chamfered corners
		G4TwoVector v41(-bp_x1, bp_y1 - bp_chamfer);
		G4TwoVector v42(-bp_x1 + bp_chamfer, bp_y1);

		std::vector<G4TwoVector> shieldingcorner_poly
		{ v11, v21, v22, v31, v32, v41, v42 };
		G4TwoVector offset(0.0, 0.0);

		G4ExtrudedSolid* bp01_box =
			new G4ExtrudedSolid("base_plate_shielding_corner", shieldingcorner_poly, 3.81 * mm, offset, 1.0, offset, 1.0);


		G4LogicalVolume* base_plate_01_logical = new G4LogicalVolume(bp01_box, G4Material::GetMaterial("W95"),
			"base_plate_01", 0, 0, 0);
		base_plate_01_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.0, 0.3, 0.3, 0.6)));

		{
			//  first one 
			G4RotationMatrix* rot = new G4RotationMatrix();
			new G4PVPlacement(rot,
				G4ThreeVector(-bp_x + holepos_x, -bp_y + holepos_y,
					BasePlate_z - BasePlate1_zRef + BasePlate_center - parent_z_world),
				base_plate_01_logical, "base_plate_shielding_corner", parent_logical, false, 0);
		}
		{
			//  second one  
			G4RotationMatrix* rot = new G4RotationMatrix();
			rot->rotateX(180. * deg);
			new G4PVPlacement(rot,
				G4ThreeVector(-bp_x + holepos_x, +bp_y - holepos_y,
					BasePlate_z - BasePlate1_zRef + BasePlate_center - parent_z_world),
				base_plate_01_logical, "base_plate_shielding_corner", parent_logical, false, 1);
		}
	}
	{
		// Base Plate 02  "tungsten lower shielding component" -- fig11 of MCdatapackage
		G4GDMLParser* parser = new G4GDMLParser();
		auto base_plate_02_file = gdml_folder / "BasePlate_02.gdml";
		parser->Read(base_plate_02_file.generic_u8string(), false);
		G4LogicalVolume* base_plate_02_logical = parser->GetVolume("BasePlate_02_LogVol");
		base_plate_02_logical->SetName("base_plate_02_gdml");
		base_plate_02_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.6)));
		delete parser;

		// gdml is symmetric about 0
		G4double BasePlate2_xpos = 133.985 * mm;//from MCdatapackage, and gdml symmetric

		G4double BasePlate2_zRef = 8.89 * mm;  //.35"
		G4double BasePlate2_zpos = 7.62 * mm;  //.3"

		G4RotationMatrix* rot = new G4RotationMatrix();
		rot->rotateZ(90.0 * deg);
		rot->rotateY(180.0 * deg);

		new G4PVPlacement(rot,
			G4ThreeVector(BasePlate2_xpos, 0.,
				BasePlate_z - BasePlate2_zpos + BasePlate2_zRef + BasePlate_center - parent_z_world),
			base_plate_02_logical, "base_plate_02_gdml", parent_logical, false, 0);
	}
	{
		G4GDMLParser* parser = new G4GDMLParser();
		auto mlc_corner_file = gdml_folder / "CornerClipper.gdml";
		parser->Read(mlc_corner_file.generic_u8string(), false);
		G4LogicalVolume* mlc_corner_logical = parser->GetVolume("CornerClipper_LogVol");
		mlc_corner_logical->SetName("mlc_corner_clipper_gdml");
		mlc_corner_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.0, 1.0, 1.0, 0.9)));
		delete parser;
		//5.7" on diagonal, from MC data package
		G4double bp_hole_pos = 144.78 * mm / std::sqrt(2.);
		// positions of mounting hole, relative to square corner
		G4double xpos = 22.631 * mm;
		G4double ypos = 29.058 * mm;// 7.112*mm;

		G4double thet;

		//gdml ref point is at the square corner (x,y); center for z
		G4double cornX = -39.7319 * mm;       // from Fastrad GDML file
		G4double cornY = -24.3422 * mm;       // from Fastrad GDML file
		G4double cornZ = 3.6830 * mm;       // from Fastrad GDML file

		G4double cornZpos = base_plate_bottom_z_world + cornZ + BasePlate_center;

		{
			G4RotationMatrix* rot = new G4RotationMatrix();
			thet = 45. * deg;
			rot->rotateZ(thet);

			// this is in the -x-y quadrant
			new G4PVPlacement(rot,
				G4ThreeVector(
					// do translation in three parts
					// (1)put gdml ref point on beam axis; 
					// (2)move ref pt to the hole in baseplate
					// (3)move hole in corner
					// ref point on beam axis         // ref pt to hole   // align holes
					cornX * cos(thet) + cornY * sin(thet) - bp_hole_pos +
					(xpos * cos(thet) + ypos * sin(thet)),
					-cornX * sin(thet) + cornY * cos(thet) - bp_hole_pos +
					(-xpos * sin(thet) + ypos * cos(thet)),
					cornZpos - parent_z_world),
				mlc_corner_logical, "mlc_corner_clipper_gdml_1", parent_logical, false, 0);
		}
		{
			G4RotationMatrix* rot = new G4RotationMatrix();
			thet = 45. * deg;
			rot->rotateX(180. * deg);
			rot->rotateZ(thet);

			new G4PVPlacement(rot,
				G4ThreeVector(
					cornX * cos(thet) + cornY * sin(thet) - bp_hole_pos +
					(xpos * cos(thet) + ypos * sin(thet)),
					// for 1st and 3rd terms, sign of y is opposite to above
					// (because of 180 deg rot about X)
					+cornX * sin(thet) - cornY * cos(thet) + bp_hole_pos -
					(-xpos * sin(thet) + ypos * cos(thet)),
					cornZpos - parent_z_world),
				mlc_corner_logical, "mlc_corner_clipper_gdml_2", parent_logical, false, 1);
		}
		{
			G4RotationMatrix* rot = new G4RotationMatrix();
			thet = 225. * deg;
			rot->rotateZ(thet);

			new G4PVPlacement(rot,
				G4ThreeVector(
					cornX * cos(thet) + cornY * sin(thet) + bp_hole_pos +
					(xpos * cos(thet) + ypos * sin(thet)),
					-cornX * sin(thet) + cornY * cos(thet) + bp_hole_pos +
					(-xpos * sin(thet) + ypos * cos(thet)),
					cornZpos - parent_z_world),
				mlc_corner_logical, "mlc_corner_clipper_gdml_3", parent_logical, false, 2);
		}
		{
			G4RotationMatrix* rot = new G4RotationMatrix();
			thet = 225. * deg;
			rot->rotateX(180 * deg);
			rot->rotateZ(thet);

			new G4PVPlacement(rot,
				G4ThreeVector(
					cornX * cos(thet) + cornY * sin(thet) + bp_hole_pos +
					(xpos * cos(thet) + ypos * sin(thet)),
					// for 1st and 3rd terms, sign of y is opposite to above
					// (because of 180 deg rot about X)
					+cornX * sin(thet) - cornY * cos(thet) - bp_hole_pos -
					(-xpos * sin(thet) + ypos * cos(thet)),
					cornZpos - parent_z_world),
				mlc_corner_logical, "mlc_corner_clipper_gdml_4", parent_logical, false, 3);
		}
	}
}

void buildBackscatterShield(G4LogicalVolume* vacuum_logical)
{
	// backscatter killer
	G4Sphere* backScatter =
		new G4Sphere("backscatter_shield", 8. * mm, 9. * mm, 0. * deg, 180. * deg,
			0. * deg, 360. * deg);
	auto backscatter_shield_logical = new G4LogicalVolume(backScatter,
		G4Material::GetMaterial("G4_Galactic"),
		"backscatter", 0, 0, 0);

	G4RotationMatrix* backscatter_rotation = new G4RotationMatrix();
	backscatter_rotation->rotateX(-90.0 * deg);

	new G4PVPlacement(backscatter_rotation,
		G4ThreeVector(0., 0., -5. * mm),
		backscatter_shield_logical, "backscatter", vacuum_logical, false, 0);

	backscatter_shield_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.2, 0.0, 0.8, 0.9)));
	
}

void buildBeWindow(G4LogicalVolume* parent_logical, double parent_z_world, bool build_backscatter_shield)
{
	G4double BeThick = 0.0254 * cm;
	G4double vacuumLength = 11. * mm;
	const double offset = 2.7538 * mm;
	const double be_window_position_world_z = g_SAD + offset + BeThick / 2. + vacuumLength / 2;


	//Build vacuum
	G4Tubs* vacuum =
		new G4Tubs("vacuum", 0. * mm, 2. * cm, vacuumLength / 2., 0. * deg, 360. * deg);

	auto orbit_vacuum_logical =
		new G4LogicalVolume(vacuum, G4Material::GetMaterial("G4_Galactic"),
			"vacuum", 0, 0, 0);

	orbit_vacuum_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0, 0.2, 0.8, 0.2)));

	new G4PVPlacement(0, G4ThreeVector(0., 0., be_window_position_world_z - parent_z_world),
		orbit_vacuum_logical, "vacuum", parent_logical, false, 0);

	//Build Be window itself
	G4Tubs* BeWindow =
		new G4Tubs("BeWindow", 0., 6.477 * mm, BeThick / 2., 0. * deg, 360. * deg);
	G4LogicalVolume* be_window_logical =new G4LogicalVolume(BeWindow,
			G4Material::GetMaterial("G4_Be"), "BeWindow_LV", 0, 0, 0);

	be_window_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.5, 0.5, 1.0, 0.5)));
	
	//Be Window 2  (orbit chamber)
	new G4PVPlacement(0, G4ThreeVector(0., 0.,g_SAD + offset - parent_z_world),
		be_window_logical, "BeWindow_orbit", parent_logical, false, 0);
	if (build_backscatter_shield) {
		buildBackscatterShield(orbit_vacuum_logical);
	}
}

void buildPrimaryCollimator(G4LogicalVolume* parent_logical, double parent_z_world) 
{
  G4double pc_pos = 53.25 * mm;
	const double primary_collimator_world_z = g_SAD - pc_pos;
	const double carousel_world_z = g_SAD - pc_pos + 26.25 * mm - 12.5 * mm;

  G4RotationMatrix* PrimCollRot = new G4RotationMatrix();
  PrimCollRot->rotateX(180.0 * deg);
  PrimCollRot->rotateZ(90.0 * deg);

  const G4int pc_points = 8;
  G4double r_pcoll[pc_points] = {
      8.5 * mm, 112.55 * mm,112.55 * mm, 100.0 * mm,
    100.0 * mm,   26.0 * mm, 26.0 * mm,  21.4924 * mm
  };
  G4double z_pcoll[pc_points] = {
    -26.25 * mm, -26.25 * mm,  -0.25 * mm,  -0.25 * mm,
     12.74999 * mm,  12.74999 * mm,  26.24 * mm,  26.24 * mm
  };
  // value above is 12.74999 rather than 12.75 to solve 
  // "GeomTest problem: solid problem" when running overlap check

  G4GenericPolycone* PCCone =
    new G4GenericPolycone("primary_collimator", 0. * deg, 360. * deg, pc_points, r_pcoll, z_pcoll);
  G4LogicalVolume* primary_collimator_logical = new G4LogicalVolume(PCCone, G4Material::GetMaterial("W95"), "primary_collimator", 0, 0, 0);

  primary_collimator_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.8, 0.1, 0.8, 0.5)));

  new G4PVPlacement(PrimCollRot,
    G4ThreeVector(0., 0., primary_collimator_world_z - parent_z_world),
    primary_collimator_logical, "primary_collimator", parent_logical, false, 0);

  // carousel baseplate
  G4Tubs* carousel_bp =
    new G4Tubs("carousel_bp", 0. * mm, 175. * mm, 12.5 * mm, 0., twopi);

  G4Box* carousel_bp_box =
    new G4Box("carousel_bp_1", 175. * mm, 120. * mm, 12.5 * mm);

  G4UnionSolid* carousel_bp_union =
    new G4UnionSolid("carousel_union", carousel_bp, carousel_bp_box, 0,
      G4ThreeVector(0., -110. * mm, 0.));

  G4Tubs* carousel_bp_2 =
    new G4Tubs("carousel_bp_2", 0. * mm, 112.55 * mm, 13.5 * mm, 0, twopi);

  G4SubtractionSolid* carousel_bp_subt =
    new G4SubtractionSolid("carousel_subt", carousel_bp_union, carousel_bp_2,
      0, G4ThreeVector(0., 0., 0.));

  G4LogicalVolume* carousel_bp_logical =
    new G4LogicalVolume(carousel_bp_subt, G4Material::GetMaterial("SS_A36"),
      "carousel_bp_LV", 0, 0, 0);
  carousel_bp_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.3, 0.7, 0.1, 0.3)));

  new G4PVPlacement(0, G4ThreeVector(0., 0., carousel_world_z - parent_z_world),
    carousel_bp_logical, "carousel_bp", parent_logical, false, 0);
	
}

void buildCADPrimaryCollimator(G4LogicalVolume* parent_logical, double parent_z_world, const std::filesystem::path& stl_folder) 
{
	//Primary collimator CAD highest resolution
	//Visualization
	const auto PriColl_vis = new G4VisAttributes(G4Colour(0.8, 0.1, 0.8, 0.5));
	const auto carousel_vis = new G4VisAttributes(G4Colour(0.3, 0.7, 0.1, 0.3));
	
	static const char* filename_PriColl[] = {"100016533-03.stl","100022239-02-mod.stl"};
	std::filesystem::path filepath_pricoll = stl_folder / filename_PriColl[0];
    G4VSolid* Pricoll_solid = getTessalatedSolidFromSTL(filename_PriColl[0], filepath_pricoll, G4ThreeVector(0.0, 0.0, 0.0));
	
	G4LogicalVolume *primary_collimator_logical = new G4LogicalVolume(Pricoll_solid, G4Material::GetMaterial("W95"), "primary_collimator", 0, 0, 0);
	primary_collimator_logical->SetVisAttributes(PriColl_vis);
	new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_z_world), primary_collimator_logical, "primary_collimator", parent_logical, false, 0);
	
	std::filesystem::path filepath_carousel = stl_folder / filename_PriColl[1];
    G4VSolid* carousel_solid = getTessalatedSolidFromSTL(filename_PriColl[1], filepath_carousel, G4ThreeVector(0.0, 0.0, 0.0));
	G4LogicalVolume *carousel_bp_logical = new G4LogicalVolume(carousel_solid, G4Material::GetMaterial("SS_A36"), "carousel_bp_LV", 0, 0, 0);
	carousel_bp_logical->SetVisAttributes(carousel_vis);
	new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_z_world), carousel_bp_logical, "carousel_bp", parent_logical, false, 0);	
}

void buildCADYStageShield(G4LogicalVolume* parent_logical, double parent_z_world, const std::filesystem::path& stl_folder)
{
	//Visualization
	const auto Ystage_vis = new G4VisAttributes(G4Colour(0.2, 0.2, 0.6, 0.8));
	
	static const char* filename_ystageshield = "100022240-01.stl";
	std::filesystem::path filepath = stl_folder / filename_ystageshield;
	const G4ThreeVector offset(0.0, 0.0, 0.0);
    G4VSolid* ystageshield_solid = getTessalatedSolidFromSTL(filename_ystageshield, filepath, offset);
	G4LogicalVolume *ystage_shield_logical = new G4LogicalVolume(ystageshield_solid, G4Material::GetMaterial("Lead97Antimony"), "ystage_shield", 0, 0, 0);
	ystage_shield_logical->SetVisAttributes(Ystage_vis);
	new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_z_world), ystage_shield_logical, "ystage_shield", parent_logical, false, 0);

}

void buildYStageShield(G4LogicalVolume* parent_logical, double parent_z_world)
{
	const double ystage_shield_world_z = g_SAD - 96.5 * mm;

	G4Cons* ystage_1 =
		new G4Cons("ystage_1", 23.5 * mm, 145. * mm, 36.4325 * mm, 145. * mm,
			16.5 * mm, 0. * deg, 360. * deg);
	//not supposed to have overlapping surfaces because they create
	// "false" surfaces
	G4double trap_delta = 0.1 * mm;

	G4Tubs* ystage_4 =
		new G4Tubs("ystage_4", 26.5 * mm, 145. * mm, 6.5 * mm, 0. * deg, 360. * deg);

	G4UnionSolid* ystage_5 =
		new G4UnionSolid("ystage_5", ystage_1, ystage_4,
			0, G4ThreeVector(0., 0., -23. * mm + trap_delta));

	G4Box* ystage_6 =
		new G4Box("ystage_6", 145. * mm, 140. * mm, 23. * mm);

	G4UnionSolid* ystage_7 =
		new G4UnionSolid("ystage_7", ystage_5, ystage_6, 0,
			G4ThreeVector(0., 176.4325 * mm, -6.5 * mm));

	// cutout groove as an extruded polygon
	G4double point1 = 2.5 * mm;
	G4double point2 = 19. * mm;
	G4double point3 = 11. * mm;  // half length in beam direction

	G4TwoVector v1(point1, point3);
	G4TwoVector v2(point2, -point3);
	G4TwoVector v3(-point2, -point3);
	G4TwoVector v4(-point1, point3);
	std::vector<G4TwoVector> polygon{ v1, v2, v3, v4 };

	G4double l1 = 30. * mm;
	G4TwoVector vv1(0. * mm, 0. * mm);
	G4ExtrudedSolid::ZSection zs1(-l1, vv1, 1.);
	G4ExtrudedSolid::ZSection zs2(l1, vv1, 1.);
	std::vector<G4ExtrudedSolid::ZSection> zsections{ zs1, zs2 };
	G4ExtrudedSolid* ystage_2 =
		new G4ExtrudedSolid("scCenter", polygon, zsections);

	//subtract several of the extruded solid, around an arc
	std::vector<G4RotationMatrix*> ystage_rot;
	std::vector<G4SubtractionSolid*>  subtsolids;
	G4double ext_rad = 100. * mm;  // radius of groove
	G4double ext_angle;
	G4double ext_dx, ext_dy;

	const G4int num_grooves = 7;  // can't visualize more than 7!
	for (G4int i = 0; i <= num_grooves; i += 1) {
		ext_angle = -90. * deg + 180. * deg * G4double(i) / G4double(num_grooves - 1);
		ystage_rot.push_back(new G4RotationMatrix());
		ystage_rot[i]->rotateX(90. * deg);
		ystage_rot[i]->rotateY(90. * deg);
		ystage_rot[i]->rotateY(ext_angle);

		ext_dx = ext_rad * std::sin(ext_angle);
		ext_dy = ext_rad * (1. - std::cos(ext_angle));
		if (i == 0) {
			subtsolids.push_back(
				new G4SubtractionSolid("ystage_3a", ystage_7, ystage_2,
					ystage_rot[i],
					G4ThreeVector(ext_dx, ext_dy, 5.5 * mm + trap_delta)));
		}
		else {
			subtsolids.push_back(
				new G4SubtractionSolid("ystage_3a", subtsolids[i - 1], ystage_2,
					ystage_rot[i],
					G4ThreeVector(ext_dx, ext_dy, 5.5 * mm + trap_delta)));
		}
	}

	G4LogicalVolume* ystage_shield_logical = new G4LogicalVolume(subtsolids[num_grooves - 1],
			G4Material::GetMaterial("Lead97Antimony"), "ystage_shield", 0, 0, 0);
	ystage_shield_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.2, 0.2, 0.6, 0.8)));

	G4RotationMatrix* ystage2_Rot2 = new G4RotationMatrix();
	ystage2_Rot2->rotateX(180. * deg);

	new G4PVPlacement(ystage2_Rot2,
		G4ThreeVector(0., 0., ystage_shield_world_z - parent_z_world),
		ystage_shield_logical, "ystage_shield", parent_logical, false, 0);
}

JawVolumes buildJaws(G4LogicalVolume* parent_logical, double parent_z_world, const std::filesystem::path& stl_folder) {
	JawVolumes volumes{};
	volumes.parent_z_world = parent_z_world;
	static const char* filenames_jaws[] = {"00876386-02r-a-10.stl","00876386-02r-b-10.stl","00885474-01r-b-10_fixed.stl","00885474-01r-a-10_fixed.stl"};
	//****************************
	// upper jaws
	//****************************
	//Visualization
	const auto upperjaws_vis = new G4VisAttributes(G4Colour(0.2, 0.3, 1.0, 0.5));
	std::filesystem::path filepath = stl_folder / filenames_jaws[0];
    G4VSolid* jawupperY1_solid = getTessalatedSolidFromSTL(filenames_jaws[0], filepath, G4ThreeVector(0.0, 0.0, 0.0));
	G4LogicalVolume *fUpperJawY1_LV = new G4LogicalVolume(jawupperY1_solid,  G4Material::GetMaterial("W95"), "upperJawY1_LV", 0, 0, 0);
	fUpperJawY1_LV->SetVisAttributes(upperjaws_vis);
	volumes.y1 = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_z_world), fUpperJawY1_LV, "UpperJaw1", parent_logical, false, 0);

	filepath = stl_folder / filenames_jaws[1];
    G4VSolid* jawupperY2_solid = getTessalatedSolidFromSTL(filenames_jaws[1], filepath, G4ThreeVector(0.0, 0.0, 0.0));
	G4LogicalVolume *fUpperJawY2_LV = new G4LogicalVolume(jawupperY2_solid,  G4Material::GetMaterial("W95"), "upperJawY2_LV", 0, 0, 0);
	volumes.y2 = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_z_world), fUpperJawY2_LV, "UpperJaw2", parent_logical, false, 0);

	//****************************
	// lower jaws
	//****************************
	//Visualization
	const auto lowerjaws_vis = new G4VisAttributes(G4Colour(0.2, 0.3, 1.0, 0.5));
	filepath = stl_folder / filenames_jaws[2];
    G4VSolid* jawlowerX1_solid = getTessalatedSolidFromSTL(filenames_jaws[2], filepath, G4ThreeVector(0.0, 0.0, 0.0));
	G4LogicalVolume *fLowerJawX1_LV = new G4LogicalVolume(jawlowerX1_solid,  G4Material::GetMaterial("W95"), "lowerJawX1_LV", 0, 0, 0);
	fLowerJawX1_LV->SetVisAttributes(lowerjaws_vis);
	volumes.x1 = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_z_world), fLowerJawX1_LV, "LowerJaw1", parent_logical, false, 0);

	filepath = stl_folder / filenames_jaws[3];
    G4VSolid* jawlowerX2_solid = getTessalatedSolidFromSTL(filenames_jaws[3], filepath, G4ThreeVector(0.0, 0.0, 0.0));
	G4LogicalVolume *fLowerJawX2_LV = new G4LogicalVolume(jawlowerX2_solid,  G4Material::GetMaterial("W95"), "lowerJawX2_LV", 0, 0, 0);
	fLowerJawX2_LV->SetVisAttributes(G4Colour(0.3, 0.7, 0.1, 0.3));
	volumes.x2 = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_z_world), fLowerJawX2_LV, "LowerJaw2", parent_logical, false, 0);

	// miscellaneous stuff around jaws

	G4Box* jaw_surround1 =
		new G4Box("jaw_surround1", 215. * mm, 215. * mm, 103. * mm);

	G4Box* jaw_surround2 =
		new G4Box("jaw_surround2", 195. * mm, 195. * mm, 104. * mm);

	G4SubtractionSolid* jaw_surround3 =
		new G4SubtractionSolid("jaw_surround3", jaw_surround1, jaw_surround2,
			0, G4ThreeVector());

	G4Box* jaw_surround4 =
		new G4Box("jaw_surround4", 210. * mm, 120. * mm, 60. * mm);

	G4SubtractionSolid* jaw_surround5 =
		new G4SubtractionSolid("jaw_surround5", jaw_surround3, jaw_surround4,
			0, G4ThreeVector(0., 0., -52. * mm));

	G4LogicalVolume* jawSurround_LV =
		new G4LogicalVolume(jaw_surround5,
			G4Material::GetMaterial("misc"),
			"jaw_surround_LV", 0, 0, 0);

	G4VisAttributes* VisAtt_jawsurround =
		new G4VisAttributes(G4Colour(0.4, 0.4, 0.1, 0.1));
	jawSurround_LV->SetVisAttributes(VisAtt_jawsurround);

	new G4PVPlacement(0, G4ThreeVector(0., 0., 652.0 * mm - parent_z_world),
		jawSurround_LV, "jaw_surround", parent_logical, false, 1);
	return volumes;
}

void setJawPositions(const JawVolumes& volumes, const JawPositions& pos, double parent_z_world) {
	//y-jaws (upper) moves in an arc with radius of curvature
    const G4double R_y = 27.9515*cm;
    const G4double ysize  = 5.969*cm; //half y-size of jaw
    const G4double xsize = 63.43*mm;  // half-size of jaw 
    const G4double zsize = 38.862*mm;  // half-size of jaw
    const G4double y_CAD_offset = 2.78*mm; 
	G4double theta_orig_Y1, theta_0_Y1, theta_0_Y2;
    theta_0_Y1 = std::atan((-pos.y1)/g_SAD);
    theta_0_Y2 = std::atan((-pos.y2)/g_SAD);
    theta_orig_Y1 = std::atan((y_CAD_offset)/R_y);

	G4RotationMatrix* rotateY1 = new G4RotationMatrix();
    G4RotationMatrix* rotateY2 = new G4RotationMatrix();
    rotateY1->rotateX(theta_0_Y1 );
    rotateY2->rotateX(theta_0_Y2+ theta_orig_Y1);
    volumes.y1->SetRotation(rotateY1);
    volumes.y2->SetRotation(rotateY2);

	//Move so that z-position tr stays constant while moving in x-direction
	//and rotating s.t. jaw y-z plane creates angle theta with beam cone
	//jaw rotates about trunion. trunion z is constant
	G4double theta_X1, theta_X2;
    G4double R_x = 367*mm; //633.477*mm; // Lower X jaws face moves in an arc with radius of curvature
    G4double radius = 7.7*mm; //radius of the screwhole, edges to screwhole measured from STL file
    G4double S_x = 21.5*mm + radius; //measured from STL file
    G4double S_z_up = 12.70*mm + 13.66*mm + radius; //measured from STL file
    G4double S_z_down = 12.70*mm + 23.13*mm + radius; //measured from STL file
    G4double screw_height = g_SAD - R_x - S_z_up; //height of the center of the screw
    theta_X1  = std::atan((pos.x1)/g_SAD);
    theta_X2 = std::atan((pos.x2)/g_SAD);
	// jaw rotates about screw. z is constant
    G4double zcenterX1 =  (S_z_up + R_x) - (S_z_up + R_x) *std::cos(theta_X1);
    G4double zcenterX2 = (S_z_up + R_x) - (S_z_up + R_x) *std::cos(theta_X2);
	G4double xcenterX1 = ((S_z_up + R_x) - (S_z_up + R_x) *std::cos(theta_X1)) * std::tan((pos.x1)/g_SAD);
	G4double xcenterX2 = ((S_z_up + R_x) - (S_z_up + R_x)*std::cos(theta_X2)) * std::tan((pos.x2)/g_SAD);

	G4RotationMatrix* rotateX1 = new G4RotationMatrix();
    G4RotationMatrix* rotateX2 = new G4RotationMatrix();
    rotateX1->rotateY(theta_X1);
    rotateX2->rotateY(theta_X2);
    volumes.x1->SetRotation(rotateX1);
    volumes.x2->SetRotation(rotateX2);
    volumes.x1->SetTranslation(G4ThreeVector(xcenterX1, 0.0, g_SAD - parent_z_world - zcenterX1)); 
    volumes.x2->SetTranslation(G4ThreeVector(xcenterX2, 0.0, g_SAD - parent_z_world - zcenterX2)); 
}
void buildMCBackscatterPlate(G4LogicalVolume* parent_logical, double parent_z_world)
{
    const double backscatter_plate_world_z = 824.124*mm;
    G4Tubs* backscatter_plate_tubs = new G4Tubs("backscatter_plate", 0., 100.*mm, 0.3*mm, 0., 360 * deg);
    G4LogicalVolume* backscatter_plate_logical = new G4LogicalVolume(backscatter_plate_tubs, G4Material::GetMaterial("SS304"), "backscatter_plate", 0, 0, 0);
    backscatter_plate_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.1, 0.7, 0.3, 0.5)));
    new G4PVPlacement(0, G4ThreeVector(0., 0., backscatter_plate_world_z - parent_z_world), backscatter_plate_logical, "backscatter_plate", parent_logical, false, 0);

}

VolumeNameAndTraversalFlag TrueBeamTable[] = {
	{"world", TraversedGeometry::NONE},
	{"gantry", TraversedGeometry::NONE},
	{"collimator", TraversedGeometry::NONE},
	{"mylar_window", TraversedGeometry::NONE},
	{"upperJawY1_LV", TraversedGeometry::UPPER_JAW},
	{"upperJawY2_LV", TraversedGeometry::UPPER_JAW},
	{"lowerJawX1_LV", TraversedGeometry::LOWER_JAW},
	{"lowerJawX2_LV", TraversedGeometry::LOWER_JAW},
	{"jaw_surround_LV", TraversedGeometry::NONE},
	{"BeWindow6X_LV", TraversedGeometry::NONE},
	{"vacuum_chamber", TraversedGeometry::NONE},
	{"target_button", TraversedGeometry::TARGET},
	{"target_nicoro", TraversedGeometry::TARGET},
	{"target_block_top", TraversedGeometry::TARGET},
	{"xray_window", TraversedGeometry::NONE},
	{"ff_carousel", TraversedGeometry::FLATTENING_FILTER},
	{"open_port", TraversedGeometry::FLATTENING_FILTER},
	{"ff_6X", TraversedGeometry::FLATTENING_FILTER},
	{"primary_collimator", TraversedGeometry::PRIMARY_COLLIMATOR},
	{"carousel_bp_LV", TraversedGeometry::NONE},
	{"ystage_shield", TraversedGeometry::YSTAGE_SHIELD},
	{"shieldColl_LV", TraversedGeometry::SHIELD_COLLIMATOR},
	{"shield_collimator", TraversedGeometry::SHIELD_COLLIMATOR},
	{"vacuum", TraversedGeometry::NONE},
	{"backscatter", TraversedGeometry::NONE},
	{"BeWindow_LV", TraversedGeometry::NONE},
	{"ic_mother", TraversedGeometry::IC},
	{"ic_body", TraversedGeometry::IC},
	{"ic_window", TraversedGeometry::IC},
	{"ic_window_cu", TraversedGeometry::IC},
	{"ic_topring", TraversedGeometry::IC},
	{"IC_center_ring", TraversedGeometry::IC},
	{"ic_electrode_hv", TraversedGeometry::IC},
	{"ic_cu_plating", TraversedGeometry::IC},
	{"ic_cu_sig", TraversedGeometry::IC},
	{"ic_electrode_hv", TraversedGeometry::IC},
	{"ic_cu_sig", TraversedGeometry::IC},
	{"IC_windowring_LV", TraversedGeometry::IC},
	{"ic_dose", TraversedGeometry::IC},
	{"mlc_leaf", TraversedGeometry::MLC},
	{"base_plate_gdml", TraversedGeometry::BASEPLATE},
	{"base_plate_01", TraversedGeometry::BASEPLATE},
	{"base_plate_02_gdml", TraversedGeometry::BASEPLATE02},
	{"mlc_corner_clipper_gdml", TraversedGeometry::BASEPLATE},
	{"backscatter_plate", TraversedGeometry::NONE}
};

TreatmentHeadDetector::TreatmentHeadDetector(MLCType mlc_type, EnergyMode energy_mode, 
	const std::string& sd_monitor_chamber_name, const fs::path& gdml_path, const fs::path& stl_path) 
		: m_mlc_type(mlc_type), m_energy_mode(energy_mode),
		m_gdml_path(gdml_path), m_stl_path(stl_path), m_sd_monitor_chamber_name(sd_monitor_chamber_name),
		m_monitor_chamber(nullptr)
{}

G4VPhysicalVolume* TreatmentHeadDetector::Construct() {
	//We actually want to store physical volumes for rotations
	buildMaterials();
	G4NistManager* NISTman = G4NistManager::Instance();
	auto world_material = NISTman->FindOrBuildMaterial("G4_AIR");
	const double collimator_position_z = 430.0 * mm;
	const double gantry_position_z = 690.0 * mm;

	G4LogicalVolume* world_logical = nullptr;
	G4VPhysicalVolume* world_physical = nullptr;
	//World volume
	{
		G4ThreeVector world_size = G4ThreeVector(1000. * cm, 1000. * cm, 1000. * cm);
		G4Box* world_box = new G4Box("world", world_size.x() / 2., world_size.y() / 2., world_size.z() / 2.);
		world_logical = new G4LogicalVolume(world_box, world_material, "world", 0, 0, 0);
		world_physical = new G4PVPlacement(0, G4ThreeVector(), world_logical, "world", 0, false, 0);
	}
	//Gantry volume
	G4LogicalVolume* gantry_logical = nullptr;
	{
		G4Box* gantry_box = new G4Box("gantry", 45. * cm, 45. * cm, 65. * cm);
		gantry_logical = new G4LogicalVolume(gantry_box, world_material, "gantry", 0, 0, 0);
		m_gantry = new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, gantry_position_z), 
			gantry_logical, "gantry", world_logical, false, 0);
	}
	//Collimator volume
	G4LogicalVolume* collimator_logical = nullptr;
	{
		
		G4Tubs* collimator_tubs = new G4Tubs("collimator", 0. * mm, 314. * mm, 390. * mm, 0. * deg, 360. * deg);
		collimator_logical = new G4LogicalVolume(collimator_tubs, world_material, "collimator", 0, 0, 0);
		m_collimator = new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, collimator_position_z - gantry_position_z), 
			collimator_logical, "collimator", gantry_logical, false, 0);
	}
	const bool build_photon = true;
	if (build_photon) {
		Target target;
		if((m_energy_mode == EnergyMode::X06) || (m_energy_mode == EnergyMode::X06_FFF))
			target = Target::LOW_ENERGY;
		else
			throw std::runtime_error("Unknown energy mode for Truebeam target");
		buildTarget(gantry_logical, gantry_position_z, target);

		FlatteningFilter flattening_filter;
		//With lot of energies, this can be cleaned up 
		if(m_energy_mode == EnergyMode::X06_FFF)
			flattening_filter = FlatteningFilter::OPEN;
		else if(m_energy_mode == EnergyMode::X06)
			flattening_filter = FlatteningFilter::ENERGY_6X;
		else
			throw std::runtime_error("Unknown energy mode for Truebeam flattening filter");
		buildFlatteningFilter(gantry_logical, gantry_position_z, flattening_filter);

		buildPrimaryCollimator(gantry_logical, gantry_position_z);
		buildYStageShield(gantry_logical, gantry_position_z);
		fs::path pricoll_path = m_stl_path / "TB_UPPERCOLLIMATOR" / "binary";
		//buildCADPrimaryCollimator(gantry_logical, gantry_position_z, pricoll_path);
		//buildCADYStageShield(gantry_logical, gantry_position_z, pricoll_path);
	}
	fs::path shieldcoll_path = m_stl_path / "TB_SHIELDCOLLIMATOR" / "binary";
	buildShieldingCollimator(collimator_logical, collimator_position_z, shieldcoll_path);
	const bool build_backscatter_shield = true;
	buildBeWindow(gantry_logical, gantry_position_z, build_backscatter_shield);
	m_monitor_chamber = buildIonisationChamber(gantry_logical, gantry_position_z);
	buildMCBackscatterPlate(gantry_logical, gantry_position_z);

	if (m_mlc_type == MLCType::MILLENNIUM120) {
		fs::path mlc_path = m_stl_path / "MLC120-LEAF-STL" / "binary";
		m_mlc_leaf_volumes = buildMilleniumMLC(collimator_logical, collimator_position_z, mlc_path);
	}
	else if (m_mlc_type == MLCType::HD){
		fs::path mlc_path = m_stl_path / "HD-MLC120-LEAF-STL" / "binary";
		m_mlc_leaf_volumes = buildHDMLC(collimator_logical, collimator_position_z, mlc_path);
	}
	else if (m_mlc_type == MLCType::NONE) {
		//Do nothing
	}
	else { 
		static_assert(truebeam::mlcs.size() == 3);
	}

	buildMylarWindow(collimator_logical, collimator_position_z);
	fs::path jaw_path = m_stl_path / "TB_JAWS" / "binary";
	m_jaw_volumes = buildJaws(collimator_logical, collimator_position_z, jaw_path);
	setJawPositions(m_jaw_volumes, JawPositions{ -20.0*cm, 20.0*cm, -20.0*cm, 20.0*cm },collimator_position_z); 
	buildBasePlate(collimator_logical, collimator_position_z, m_gdml_path);
	m_id_to_traversed = generateVolumeToTraversed(TrueBeamTable, sizeof(TrueBeamTable) / sizeof(TrueBeamTable[0]), world_logical);
	return world_physical;
}

void TreatmentHeadDetector::setState(const Plan& plan, size_t pt_idx)
{
	G4GeometryManager::GetInstance()->OpenGeometry(m_gantry);
	//auto gantry_rotation = G4RotateY3D(state.th_rotations.gantry_angle);
	const auto& state = plan.getJawAndMLC(TreatmentHeadType::TRUEBEAM, pt_idx);
	const double gantry_position_z = 690.0 * mm;
	const double collimator_position_z = 430.0 * mm;
	auto gantry_transform =  G4RotateY3D(state.th_rotations.gantry_angle)*G4TranslateZ3D(gantry_position_z);
	m_gantry->SetTranslation(gantry_transform.getTranslation());
	m_gantry->SetRotation(new G4RotationMatrix(gantry_transform.getRotation().inverse()));
	auto collimator_transform = G4RotateZ3D(state.th_rotations.collimator_angle);
	m_collimator->SetRotation(new G4RotationMatrix(collimator_transform.getRotation().inverse())); //Could just use -angle
	setJawPositions(m_jaw_volumes, state.jaw_positions, collimator_position_z);
	if (m_mlc_type == MLCType::MILLENNIUM120) {
		setMilleniumMLCLeafPositions(m_mlc_leaf_volumes, state.leaf_positions);
	}
	else if (m_mlc_type == MLCType::HD) {
		setHDMLCLeafPositions(m_mlc_leaf_volumes, state.leaf_positions);
	}
	else if(m_mlc_type == MLCType::NONE) {
		//Do nothing
	}
	else {
		static_assert(truebeam::mlcs.size() == 3);
	}
	
	G4GeometryManager::GetInstance()->CloseGeometry(m_gantry);
}

void TreatmentHeadDetector::ConstructSDandField() {
	auto sd_monitor_chamber = new G4MultiFunctionalDetector(m_sd_monitor_chamber_name);
	G4SDManager::GetSDMpointer()->AddNewDetector(sd_monitor_chamber);
	SetSensitiveDetector(m_monitor_chamber, sd_monitor_chamber);
	auto monitor_chamber_dose = new G4PSDoseDeposit("monitor_chamber_dose");
	sd_monitor_chamber->RegisterPrimitive(monitor_chamber_dose);
}

}
