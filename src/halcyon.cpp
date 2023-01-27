//HALCYON MACHINE

#include "G4PhysListFactory.hh"
#include "G4SystemOfUnits.hh"


#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"

#include "G4SDManager.hh"
#include "G4PSDoseDeposit3D.hh"
#include "G4PVReplica.hh"
#include "G4Material.hh"
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

namespace halcyon {
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
}


MLCLeafVolumesProximal buildProximalLeafBank(G4LogicalVolume* parent_logical, const double parent_world_z, double parent_gantry_z, const std::filesystem::path& stl_folder) 
{
    // Rotation problem with stl files, must convert (Y -> Z)
    auto rotation_proximal = new G4RotationMatrix();
    rotation_proximal->rotateX(-90.*deg);
    //Visualization
	const auto bank_X1_vis = new G4VisAttributes(G4Colour(0.3, 0.3, 0.8, 0.9));
	const auto bank_X2_vis = new G4VisAttributes(G4Colour(0.3, 0.3, 0.8, 0.9));
    const auto outerleaves_vis = new G4VisAttributes(G4Colour(0.3, 0.8, 0.3, 0.9));

    //Order of files is meaningful so that leaf position indices are correct
    static const char* mlc_proximal_bank_X1[] = 
    {"P1018965015_AP_001.STL",
    "P1018965009_AP_002.STL",
    "P1018965013_AP_003.STL",
    "P1018965011_AP_004.STL",
    "P1018965015_AP_005.STL",
    "P1018965009_AP_006.STL",
    "P1018965013_AP_007.STL",
    "P1018965011_AP_008.STL",
    "P1018965015_AP_009.STL",
    "P1018965009_AP_010.STL",
    "P1018965013_AP_011.STL",
    "P1018965011_AP_012.STL",
    "P1018965015_AP_013.STL",
    "P1018965009_AP_014.STL",
    "P1018965013_AP_015.STL",
    "P1018965011_AP_016.STL",
    "P1018965015_AP_017.STL",
    "P1018965009_AP_018.STL",
    "P1018965013_AP_019.STL",
    "P1018965011_AP_020.STL",
    "P1018965015_AP_021.STL",
    "P1018965009_AP_022.STL",
    "P1018965013_AP_023.STL",
    "P1018965011_AP_024.STL",
    "P1018965015_AP_025.STL",
    "P1018965009_AP_026.STL",
    "P1018965013_AP_027.STL", 
    "P1018965011_AP_028.STL",
    "P1018965015_AP_029.STL"};
    static const char* mlc_proximal_bank_X2[] = 
    {"P1018965016_BP_001.STL",
    "P1018965010_BP_002.STL",
    "P1018965014_BP_003.STL",
    "P1018965012_BP_004.STL",
    "P1018965016_BP_005.STL",
    "P1018965010_BP_006.STL",
    "P1018965014_BP_007.STL", 
    "P1018965012_BP_008.STL",
    "P1018965016_BP_009.STL", 
    "P1018965010_BP_010.STL",
    "P1018965014_BP_011.STL",
    "P1018965012_BP_012.STL",
    "P1018965016_BP_013.STL", 
    "P1018965010_BP_014.STL",
    "P1018965014_BP_015.STL",
    "P1018965012_BP_016.STL",
    "P1018965016_BP_017.STL", 
    "P1018965010_BP_018.STL",
    "P1018965014_BP_019.STL",
    "P1018965012_BP_020.STL",
    "P1018965016_BP_021.STL",
    "P1018965010_BP_022.STL",
    "P1018965014_BP_023.STL", 
    "P1018965012_BP_024.STL",
    "P1018965016_BP_025.STL",
    "P1018965010_BP_026.STL", 
    "P1018965014_BP_027.STL", 
    "P1018965012_BP_028.STL",
    "P1018965016_BP_029.STL"};

    static const char* mlc_proximal_outerleaves[] = 
    {"P1018769001_2_1.stl",
    "P1018769001_3_1.stl"};

    
    auto mlc_volumes_proximal = MLCLeafVolumesProximal();
    // BANK A
    {
        size_t num_bank_X1 = sizeof(mlc_proximal_bank_X1) / sizeof(*mlc_proximal_bank_X1);
        for (size_t i = 0; i < num_bank_X1; ++i) {
		    std::filesystem::path filepath = stl_folder / mlc_proximal_bank_X1[i];
            const G4ThreeVector offset(-3.2, 0.0, 0.0); 
            G4VSolid* leaf_solid = getTessalatedSolidFromSTL(mlc_proximal_bank_X1[i], filepath, offset);
            G4LogicalVolume *leaf_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("W95"), "mlc_proximal_leaf", 0, 0, 0);
            leaf_logical->SetVisAttributes(bank_X1_vis);
            auto volume = new G4PVPlacement(rotation_proximal, G4ThreeVector(0., 0., g_SAD - parent_gantry_z - parent_world_z), leaf_logical, "mlc_proximal_leaf", parent_logical, false, 0);
            mlc_volumes_proximal.bank_X1.push_back(volume);
        }
    }
    //BANK B
    {
        size_t num_bank_X2 = sizeof(mlc_proximal_bank_X2) / sizeof(*mlc_proximal_bank_X2);
        for (size_t i = 0; i < num_bank_X2; ++i){  
            std::filesystem::path filepath = stl_folder / mlc_proximal_bank_X2[i];
            const G4ThreeVector offset(3.2, 0.0, 0.0);
            G4VSolid* leaf_solid = getTessalatedSolidFromSTL(mlc_proximal_bank_X2[i], filepath, offset);
            G4LogicalVolume *leaf_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("W95"), "mlc_proximal_leaf", 0, 0, 0);
            leaf_logical->SetVisAttributes(bank_X2_vis);
            auto volume = new G4PVPlacement(rotation_proximal, G4ThreeVector(0., 0., g_SAD - parent_gantry_z - parent_world_z), leaf_logical, "mlc_proximal_leaf", parent_logical, false, 0);
            mlc_volumes_proximal.bank_X2.push_back(volume);
        }
    }
    {
        size_t num_outerleaves = sizeof(mlc_proximal_outerleaves) / sizeof(*mlc_proximal_outerleaves);
        for (size_t i = 0; i < num_outerleaves; ++i) {
		    std::filesystem::path filepath = stl_folder / mlc_proximal_outerleaves[i];
            const G4ThreeVector offset(0.0, 0.0, 0.0);
            G4VSolid *leaf_solid = getTessalatedSolidFromSTL(mlc_proximal_outerleaves[i], filepath, offset);
            G4LogicalVolume *leaf_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("W95"), "mlc_proximal_outer_leaf", 0, 0, 0);
            leaf_logical->SetVisAttributes(outerleaves_vis);
            auto volume = new G4PVPlacement(rotation_proximal, G4ThreeVector(0., 0.,  g_SAD - parent_gantry_z - parent_world_z), leaf_logical, "mlc_proximal_outer_leaf", parent_logical, false, 0);
        }
    }  
    return mlc_volumes_proximal;
}

void SetProximalPositions(const MLCLeafVolumesProximal& volumes, const LeafPositions& positions) 
{
    // position is positive for non-overtravel, for both leaf banks
    const size_t num_leaves_in_bank = 29;
    const bool valid = (volumes.bank_X2.size() == num_leaves_in_bank) || (volumes.bank_X1.size() == num_leaves_in_bank) 
    || (positions.bank_X2.size() == num_leaves_in_bank) || (positions.bank_X1.size() == num_leaves_in_bank);
    if (!valid)
        throw std::runtime_error("Number of Proximal leaf volumes and leaf positions is inconsistent");
    
    const std::vector<double>* bank_positions[] = { &positions.bank_X1, &positions.bank_X2 };
    const std::vector<G4VPhysicalVolume*>* bank_volumes[] = { &volumes.bank_X1, &volumes.bank_X2};
    
    const G4double zp = 349.*mm;         // source-centre of leaf distance
    const G4double ziso = 1000.*mm;         // source-isocenter distance
    const G4double R = 234.*mm;
    double physical_position = 0.0;
    
    for (size_t j = 0; j < 2; ++j) {
        for (size_t i = 0; i < num_leaves_in_bank; ++i) {
            double nominal_position = bank_positions[j]->at(i);
            if (j == 0) {
                physical_position = nominal_position*zp/ziso - R* (std::sqrt(std::pow(nominal_position / ziso ,2) + 1) - 1);
            } else if (j == 1) {
                physical_position = nominal_position*zp/ziso + R* (std::sqrt(std::pow(nominal_position / ziso ,2) + 1) - 1);
            }
            physical_position *= -1;
            auto volume = bank_volumes[j]->at(i);
            auto translation = volume->GetTranslation();
            translation.setX(physical_position);
            volume->SetTranslation(translation);
        }
    }
}

MLCLeafVolumesDistal buildDistalLeafBank(G4LogicalVolume* parent_logical, const double parent_world_z, double parent_gantry_z, const std::filesystem::path& stl_folder) 
{
    // Rotation problem with stl files, must convert (Y -> Z)
    auto rotation_distal = new G4RotationMatrix();
    rotation_distal->rotateX(-90.*deg);
    //Visualization
	const auto bank_X1_vis = new G4VisAttributes(G4Colour(0.0, 0.8, 0.8, 0.4));
	const auto bank_X2_vis = new G4VisAttributes(G4Colour(0.0, 0.8, 0.8, 0.4));
    const auto outerleaves_vis = new G4VisAttributes(G4Colour(0.3, 0.8, 0.3, 0.9));

    //Order of files is meaningful so that leaf position indices are correct
    static const char* mlc_distal_bank_X1[] = 
    {"P1018964011_AD_001.STL",
    "P1018964009_AD_002.STL",
    "P1018964007_AD_003.STL",
    "P1018964011_AD_004.STL",
    "P1018964009_AD_005.STL",
    "P1018964007_AD_006.STL",
    "P1018964011_AD_007.STL",
    "P1018964009_AD_008.STL",
    "P1018964007_AD_009.STL",
    "P1018964011_AD_010.STL",
    "P1018964009_AD_011.STL",
    "P1018964007_AD_012.STL",
    "P1018964011_AD_013.STL",
    "P1018964009_AD_014.STL",
    "P1018964007_AD_015.STL",
    "P1018964011_AD_016.STL",
    "P1018964009_AD_017.STL",
    "P1018964007_AD_018.STL",
    "P1018964011_AD_019.STL",
    "P1018964009_AD_020.STL",
    "P1018964007_AD_021.STL",
    "P1018964011_AD_022.STL",
    "P1018964009_AD_023.STL",
    "P1018964007_AD_024.STL",
    "P1018964011_AD_025.STL",
    "P1018964009_AD_026.STL",
    "P1018964007_AD_027.STL",
    "P1018964011_AD_028.STL"};
    static const char* mlc_distal_bank_X2[] = 
    {"P1018964012_BD_001.STL",
    "P1018964010_BD_002.STL",
    "P1018964008_BD_003.STL",
    "P1018964012_BD_004.STL", 
    "P1018964010_BD_005.STL",
    "P1018964008_BD_006.STL",
    "P1018964012_BD_007.STL",
    "P1018964010_BD_008.STL",
    "P1018964008_BD_009.STL",
    "P1018964012_BD_010.STL",
    "P1018964010_BD_011.STL",
    "P1018964008_BD_012.STL",
    "P1018964012_BD_013.STL",
    "P1018964010_BD_014.STL",
    "P1018964008_BD_015.STL",
    "P1018964012_BD_016.STL",
    "P1018964010_BD_017.STL",
    "P1018964008_BD_018.STL",
    "P1018964012_BD_019.STL",
    "P1018964010_BD_020.STL",
    "P1018964008_BD_021.STL",
    "P1018964012_BD_022.STL",
    "P1018964010_BD_023.STL",
    "P1018964008_BD_024.STL",
    "P1018964012_BD_025.STL", 
    "P1018964010_BD_026.STL",
    "P1018964008_BD_027.STL",
    "P1018964012_BD_028.STL"};
    static const char* mlc_distal_outerleaves[] = 
    {"P1018981003_AO_000.STL",
    "P1018981003_AO_029.STL",
    "P1018981004_BO_000.STL",
    "P1018981004_BO_029.STL"};
    
    auto mlc_volumes_distal = MLCLeafVolumesDistal();
    // BANK A
    {
        size_t num_bank_X1 = sizeof(mlc_distal_bank_X1) / sizeof(*mlc_distal_bank_X1);
        for (size_t i = 0; i < num_bank_X1; ++i) {
            std::filesystem::path filepath = stl_folder / mlc_distal_bank_X1[i];
            const G4ThreeVector offset(-15.75, 0.0, 0.0);
            G4VSolid* leaf_solid = getTessalatedSolidFromSTL(mlc_distal_bank_X1[i], filepath, offset);
            G4LogicalVolume *leaf_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("W95"), "mlc_distal_leaf", 0, 0, 0);
            leaf_logical->SetVisAttributes(bank_X1_vis);
            auto volume = new G4PVPlacement(rotation_distal, G4ThreeVector(0., 0., g_SAD - parent_gantry_z - parent_world_z), leaf_logical, "mlc_distal_leaf", parent_logical, false, 0);
            mlc_volumes_distal.bank_X1.push_back(volume);
        }
    }
    // BANK B
    {
        size_t num_bank_X2 = sizeof(mlc_distal_bank_X2) / sizeof(*mlc_distal_bank_X2);
        for (size_t i = 0; i < num_bank_X2; ++i) {
		    std::filesystem::path filepath = stl_folder / mlc_distal_bank_X2[i];
            const G4ThreeVector offset(15.75, 0.0, 0.0);
            G4VSolid* leaf_solid = getTessalatedSolidFromSTL(mlc_distal_bank_X2[i], filepath, offset);
            G4LogicalVolume *leaf_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("W95"), "mlc_distal_leaf", 0, 0, 0);
            leaf_logical->SetVisAttributes(bank_X2_vis);
            auto volume = new G4PVPlacement(rotation_distal, G4ThreeVector(0., 0.,  g_SAD - parent_gantry_z - parent_world_z), leaf_logical, "mlc_distal_leaf", parent_logical, false, 0);
            mlc_volumes_distal.bank_X2.push_back(volume);
        }
    }
    {
        size_t num_outerleaves = sizeof(mlc_distal_outerleaves) / sizeof(*mlc_distal_outerleaves);
        for (size_t i = 0; i < num_outerleaves; ++i) {
		    std::filesystem::path filepath = stl_folder / mlc_distal_outerleaves[i];
            const G4ThreeVector offset(0.0, 0.0, 0.0);
            G4VSolid* leaf_solid = getTessalatedSolidFromSTL(mlc_distal_bank_X2[i], filepath, offset);
            G4LogicalVolume *leaf_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("W95"), "mlc_distal_outer_leaf", 0, 0, 0);
            leaf_logical->SetVisAttributes(outerleaves_vis);
            auto volume = new G4PVPlacement(rotation_distal, G4ThreeVector(0., 0.,  g_SAD - parent_gantry_z - parent_world_z), leaf_logical, "mlc_distal_outer_leaf", parent_logical, false, 0);
        }
    }  

    return mlc_volumes_distal;
}

void SetDistalPositions(const MLCLeafVolumesDistal& volumes, const LeafPositions& positions) 
{
    // position is positive for non-overtravel, for both leaf bank
    const size_t num_leaves_in_bank = 28;
    const bool valid = (volumes.bank_X1.size() == num_leaves_in_bank) || (volumes.bank_X2.size() == num_leaves_in_bank) 
        || (positions.bank_X2.size() == num_leaves_in_bank) || (positions.bank_X1.size() == num_leaves_in_bank);
    if (!valid)
        throw std::runtime_error("Number of Distal leaf volumes and leaf positions is inconsistent");
    
    const std::vector<double>* bank_positions[] = { &positions.bank_X1, &positions.bank_X2 };
    const std::vector<G4VPhysicalVolume*>* bank_volumes[] = { &volumes.bank_X1, &volumes.bank_X2};
    
    const G4double zd = 438.399*mm;             // source-centre of leaf distance
    const G4double ziso = 1000.*mm;             // source-isocenter distance
    const G4double R = 234.*mm;
    double physical_position = 0.0;
    for (size_t j = 0; j < 2; ++j) {
        for (size_t i = 0; i < num_leaves_in_bank; ++i) {
            double nominal_position = bank_positions[j]->at(i); 
            if (j == 0) {
                physical_position = nominal_position*zd/ziso - R* (std::sqrt(std::pow(nominal_position / ziso ,2) + 1) - 1);
            } else if (j == 1) {
                physical_position = nominal_position*zd/ziso + R* (std::sqrt(std::pow(nominal_position / ziso ,2) + 1) - 1);
            }
            physical_position *= -1;
            auto volume = bank_volumes[j]->at(i);
            auto translation = volume->GetTranslation();
            translation.setX(physical_position);
            volume->SetTranslation(translation);
        }
    }
}

void buildTarget(G4LogicalVolume* parent_logical, double parent_z_world)
{

	const double target_z_world = g_SAD;
	auto target_region = new G4Region("target");
	auto target_vis_attributes = new G4VisAttributes(G4Colour(0.2, 0.2, 0.8, 0.9));
	auto target_block_vis_attributes = new G4VisAttributes(G4Colour(0.5, 0.0, 0.5, 0.5));
    auto vacuum_vis_attributes = new G4VisAttributes(G4Colour(0, 0.2, 0.8, 0.2));


	//Materials for target
    G4NistManager* NISTman = G4NistManager::Instance();
    G4Element* elAl = NISTman->FindOrBuildElement("Al");
    G4Element* elO = NISTman->FindOrBuildElement("O");
    G4Element* elCu = NISTman->FindOrBuildElement("Cu");
    NISTman->FindOrBuildMaterial("G4_W");
    // copper glidcop
    const double density = 8.88528 * g / cm3;
    G4Material* copper_glidcop = new G4Material("COPPER_GLIDCOP", density, 3);
    copper_glidcop->AddElement(elCu, 99.600 * perCent);
    copper_glidcop->AddElement(elO, 0.1883 * perCent);
    copper_glidcop->AddElement(elAl, 0.2117 * perCent);
	
    G4double vacuum_length = 11.*mm;
    G4Tubs* vacuum = new G4Tubs("vacuum", 0.*mm, 10.*mm, vacuum_length/2., 0., 360. * deg);
    G4LogicalVolume* vacuum_logical = new G4LogicalVolume(vacuum, G4Material::GetMaterial("G4_Galactic"), "vacuum", 0, 0, 0);
    vacuum_logical->SetVisAttributes(vacuum_vis_attributes);
    new G4PVPlacement(0, G4ThreeVector(0., 0., target_z_world - parent_z_world + vacuum_length/2.), vacuum_logical, "vacuum", parent_logical, false, 0);

    G4Tubs* target_button = new G4Tubs("target_button", 0., 5.*mm, 0.3175*mm, 0., 360. * deg);
    G4LogicalVolume* target_logical = new G4LogicalVolume(target_button, G4Material::GetMaterial("G4_W"), "target_button", 0, 0, 0);
    target_logical->SetVisAttributes(target_vis_attributes);
    new G4PVPlacement(0, G4ThreeVector(0., 0., target_z_world - parent_z_world - 0.3175*mm), target_logical, "target_button", parent_logical, false, 0);


    G4double bt_thick = 2.692*mm;
    G4Tubs* target_block_top = new G4Tubs("target_block_top", 0., 10.*mm, bt_thick/2., 0., 360. * deg);
    G4LogicalVolume* target_block_top_logical = new G4LogicalVolume(target_block_top, copper_glidcop, "target_block_top", 0, 0, 0);
    target_block_top_logical->SetVisAttributes(target_block_vis_attributes);
    new G4PVPlacement(0, G4ThreeVector(0., 0., target_z_world - parent_z_world - (bt_thick/2. + 0.635*mm)), target_block_top_logical, "target_block_top", parent_logical, false, 0); 
    
    target_logical->SetRegion(target_region);
	target_block_top_logical->SetRegion(target_region);
	target_region->AddRootLogicalVolume(target_logical);
	target_region->AddRootLogicalVolume(target_block_top_logical);
}

void buildPrimaryCollimator(G4LogicalVolume* parent_logical, double parent_z_world) 
{   
    //checked from CAD model
    G4double r_1 = 172.8*mm/2.;
    G4double r_2 = 141.7*mm/2.;
    G4double length_1 = 21.5*mm;
    G4double length_2 = 71.*mm;
    const double primary_collimator_world_z = g_SAD - 15.*mm - length_1 - length_2/2.;
    G4double trd_x1 = 7.582*mm/2.; 
    G4double trd_x2 = 33.766*mm/2.;
    G4double pc_delta = 1.*mm;

    G4Tubs* pc_tubs_1 = new G4Tubs("primarycoll_tubs_1", 0., r_1, length_1/2., 0., 360.0 * deg);
    G4Tubs* pc_tubs_2 = new G4Tubs("primarycoll_tubs_2", 0., r_2, length_2/2., 0., 360.0 * deg);
    G4UnionSolid* pc_tubs = new G4UnionSolid("primarycoll_tubs", pc_tubs_2, pc_tubs_1, 0, G4ThreeVector(0.,0., length_1/2. + length_2/2.));

    G4Trd* pc_trd = new G4Trd("primarycoll_trd", trd_x2, trd_x1, trd_x2, trd_x1, (length_1 + length_2 + pc_delta)/2.);

    G4SubtractionSolid* pc_subt = new G4SubtractionSolid("primary_collimator", pc_tubs, pc_trd, 0, G4ThreeVector(0., 0., length_1/2.));
    G4LogicalVolume* primary_collimator_logical = new G4LogicalVolume(pc_subt, G4Material::GetMaterial("W95"), "primary_collimator", 0, 0, 0);
    primary_collimator_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.8, 0.1, 0.8, 0.5)));

    new G4PVPlacement(0, G4ThreeVector(0., 0., primary_collimator_world_z - parent_z_world), primary_collimator_logical, "primary_collimator", parent_logical, false, 0);
}

void buildUpPlate(G4LogicalVolume* parent_logical, double parent_z_world)
{
    //UpPlate material
    G4NistManager* NISTman = G4NistManager::Instance();
	G4Element* elZn = NISTman->FindOrBuildElement("Zn");
	G4Element* elPb = NISTman->FindOrBuildElement("Pb");
	G4Element* elFe = NISTman->FindOrBuildElement("Fe");
    G4Element* elCu = NISTman->FindOrBuildElement("Cu");
    //Brass
    const double density = 8.47*g/cm3;
    G4Material* Brass = new G4Material("Brass", density, 4);
    Brass->AddElement(elZn,37.77*perCent);
    Brass->AddElement(elPb, 0.08*perCent);
    Brass->AddElement(elFe, 0.15*perCent);
    Brass->AddElement(elCu,62.  *perCent);
    const double buildupplate_world_z = g_SAD - 112.2*mm;
    G4Tubs* buildup_tubs = new G4Tubs("builup_plate", 0., 100.*mm, 0.4*mm, 0., 360 * deg);
    G4LogicalVolume* buildup_logical = new  G4LogicalVolume(buildup_tubs, Brass, "buildup", 0, 0, 0);
    buildup_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.1, 0.7, 0.3, 0.5)));
    new G4PVPlacement(0, G4ThreeVector(0., 0., buildupplate_world_z - parent_z_world), buildup_logical, "buildup", parent_logical, false, 0);
}

void buildMCBackscatterPlate(G4LogicalVolume* parent_logical, double parent_z_world)
{
    const double backscatter_plate_world_z = g_SAD - 145.8*mm;
    G4Tubs* backscatter_plate_tubs = new G4Tubs("backscatter_plate", 0., 100.*mm, 0.3*mm, 0., 360 * deg);
    G4LogicalVolume* backscatter_plate_logical = new G4LogicalVolume(backscatter_plate_tubs, G4Material::GetMaterial("SS304"), "backscatter_plate", 0, 0, 0);
    backscatter_plate_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.1, 0.7, 0.3, 0.5)));
    new G4PVPlacement(0, G4ThreeVector(0., 0., backscatter_plate_world_z - parent_z_world), backscatter_plate_logical, "backscatter_plate", parent_logical, false, 0);

}

void buildSecondaryCollimator(G4LogicalVolume* parent_logical, double parent_z_world, double parent_gantry_z)
{   //Values checked from HALCYON CAD model
    G4double sc_length = 140.*mm;
    G4double sc_pos = 149.*mm + sc_length/2.;
    const double secondary_collimator_world_z = g_SAD - sc_pos;
    G4double sc_delta = 1.*mm;
    G4double sc_trd_x1 = 44.647/2.*mm; 
    G4double sc_trd_x2 = 83.189/2.*mm;
    G4Tubs* sc_tubs = new G4Tubs("secondarycoll_tubs", 0., 128.*mm, sc_length/2., 0., 360.0 * deg);
    G4Trd* sc_trd = new G4Trd("secondarycoll_trd", sc_trd_x2, sc_trd_x1, sc_trd_x2, sc_trd_x1, (sc_length + sc_delta)/2.);

    G4SubtractionSolid* sc_subt = new G4SubtractionSolid("secondary_collimator", sc_tubs, sc_trd, 0, G4ThreeVector());
    G4LogicalVolume* secondary_collimator_logical = new G4LogicalVolume(sc_subt, G4Material::GetMaterial("W95"), "secondary_collimator", 0, 0, 0);
    secondary_collimator_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.7, 0.2, 0.9, 0.5)));

    new G4PVPlacement(0, G4ThreeVector(0., 0., secondary_collimator_world_z - parent_gantry_z - parent_z_world), secondary_collimator_logical, "secodary_collimator", parent_logical, false, 0);
}

void buildAlPlate(G4LogicalVolume* parent_logical, double parent_z_world, double parent_gantry_z)
{   
    G4NistManager* NISTman = G4NistManager::Instance();
    G4Element* elAl = NISTman->FindOrBuildElement("Al");
    G4Element* elMg = NISTman->FindOrBuildElement("Mg");
	G4Element* elCr = NISTman->FindOrBuildElement("Cr");
    G4Element* elCu = NISTman->FindOrBuildElement("Cu");
	G4Element* elSi = NISTman->FindOrBuildElement("Si");
	G4Element* elFe = NISTman->FindOrBuildElement("Fe");
    //Aluminum5052
    const double density = 2.68*g/cm3;
    G4Material* Aluminum5052 = new G4Material("Aluminum5052", density, 6);
    Aluminum5052->AddElement(elAl, 96.9 *perCent);
    Aluminum5052->AddElement(elMg,  2.5 *perCent);
    Aluminum5052->AddElement(elCr,  0.25*perCent);
    Aluminum5052->AddElement(elCu,  0.05*perCent);
    Aluminum5052->AddElement(elSi,  0.1 *perCent);
    Aluminum5052->AddElement(elFe,  0.2 *perCent);

    const double al_plate_world_z = g_SAD - 301.25*mm;
    G4Tubs* al_plate_tubs = new G4Tubs("backscatter_plate", 0., 100.*mm, 0.25*mm, 0., 360 * deg);
    G4LogicalVolume* al_plate_logical = new G4LogicalVolume(al_plate_tubs, Aluminum5052, "al_plate", 0, 0, 0);
    al_plate_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.1, 0.7, 0.3, 0.5)));
    new G4PVPlacement(0, G4ThreeVector(0., 0., al_plate_world_z - parent_gantry_z - parent_z_world), al_plate_logical, "al_plate", parent_logical, false, 0);

}


void buildPlasticMLCBox(G4LogicalVolume* parent_logical, double parent_z_world, double parent_gantry_z)
{
    G4NistManager* NISTman = G4NistManager::Instance();
    G4Element* elC = NISTman->FindOrBuildElement("C");
    G4Element* elH = NISTman->FindOrBuildElement("H");
	G4Element* elO = NISTman->FindOrBuildElement("O");
    // PETG   PET is C_10 H_8 O_4; wikipedia says replace ethylene glycol with cyclohexane dimethanal  (CH2OH)2 -> C6H10(CH2OH)2  i.e. add C6H10
    G4Material* PETG = new G4Material("PETG", 1.38*g/cm3, 3);
    PETG->AddElement(elC, 16);
    PETG->AddElement(elH, 18);
    PETG->AddElement(elO, 4);
    G4double plastic_box_1_world_z = g_SAD - 481.3*mm;
    G4double plastic_box_2_world_z = g_SAD - 295.5*mm;
    G4Box* plastic_box = new G4Box("plastic_box", 100.*mm, 100.*mm, 0.75*mm);
    G4LogicalVolume* plastic_box_logical = new G4LogicalVolume(plastic_box, PETG, "plastic_box", 0, 0, 0);
    plastic_box_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.1, 0.7, 0.3, 0.2)));
    new G4PVPlacement(0, G4ThreeVector(0., 0., plastic_box_1_world_z - parent_gantry_z - parent_z_world), plastic_box_logical, "plastic_box", parent_logical, false, 0);
    //new G4PVPlacement(0, G4ThreeVector(0., 0., plastic_box_2_world_z - parent_gantry_z - parent_z_world), plastic_box_logical, "plastic_box", parent_logical, false, 1);

}

void buildWindow(G4LogicalVolume* parent_logical, double parent_z_world){

    //TWARON 1000 TWILL 180-190 GSM/AGMP 3612 PREPREG, C14H14N2O4
    G4NistManager* NISTman = G4NistManager::Instance();
    auto kevlar = NISTman->FindOrBuildMaterial("G4_KEVLAR");
    G4double window_world_z = g_SAD - 498.602*mm;
    G4Box* window_box = new G4Box("window_box", 200.*mm, 200.*mm, 0.7*mm);
    G4LogicalVolume* window_logical = new G4LogicalVolume(window_box, kevlar , "window", 0, 0, 0);
    window_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.8, 0.0, 0.3, 0.2)));
    new G4PVPlacement(0, G4ThreeVector(0., 0., window_world_z - parent_z_world), window_logical, "window", parent_logical, false, 0);
}



//Returns the logical volume where the dose should be recorded
G4LogicalVolume* buildIonisationChamber(G4LogicalVolume* parent_logical, double parent_z_world)
{
	const double ic_position_z_world = g_SAD - 129. * mm;
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

	//G4VisAttributes* VisAtt_IonChamber1 =
	//	new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.5));
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

VolumeNameAndTraversalFlag HalcyonTable[] = {
    {"world", TraversedGeometry::NONE},
    {"gantry", TraversedGeometry::NONE},
    {"collimator", TraversedGeometry::NONE},
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
    {"vacuum", TraversedGeometry::NONE},
    {"backscatter_plate", TraversedGeometry::NONE},
    {"target_button", TraversedGeometry::TARGET},
    {"target_block_top", TraversedGeometry::TARGET},
    {"primary_collimator", TraversedGeometry::PRIMARY_COLLIMATOR},
    {"mlc_proximal_leaf", TraversedGeometry::MLC},
    {"mlc_distal_leaf", TraversedGeometry::MLC},
    {"mlc_distal_outer_leaf", TraversedGeometry::MLC},
    {"mlc_proximal_outer_leaf", TraversedGeometry::MLC},
    {"proximal_shield", TraversedGeometry::MLC},
    {"plastic_box", TraversedGeometry::NONE},
    {"secondary_collimator", TraversedGeometry::SHIELD_COLLIMATOR},
    {"al_plate", TraversedGeometry::NONE},
    {"window", TraversedGeometry::NONE},
    {"buildup", TraversedGeometry::NONE},
};


TreatmentHeadDetector::TreatmentHeadDetector(const std::string& sd_monitor_chamber_name, 
    const fs::path& gdml_path, const fs::path& stl_path) 
    :   m_monitor_chamber(nullptr), m_gdml_path(gdml_path), 
        m_stl_path(stl_path), m_sd_monitor_chamber_name(sd_monitor_chamber_name)
{
	//auto test = getTessalatedSolidFromSTL("test", m_stl_path / "HALCYON-MLC-DISTAL" / "P1018964007_AD_003.STL", 
    //G4ThreeVector(0.0, 0.0, 0.0));
}

G4VPhysicalVolume* TreatmentHeadDetector::Construct() {
	//We actually want to store physical volumes for rotations
	buildMaterials();
	G4NistManager* NISTman = G4NistManager::Instance();
	auto world_material = NISTman->FindOrBuildMaterial("G4_AIR");
	const double collimator_position_z = -8.*mm; 
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
	G4Transform3D gantry_transform = G4Translate3D(0., 0., gantry_position_z);
	{
		G4Box* gantry_box = new G4Box("gantry", 45. * cm, 45. * cm, 65. * cm);
		gantry_logical = new G4LogicalVolume(gantry_box, world_material, "gantry", 0, 0, 0);
		m_gantry = new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, gantry_position_z), 
			gantry_logical, "gantry", world_logical, false, 0);
	}
	//Collimator volume
	G4LogicalVolume* collimator_logical = nullptr;
	{
		G4Tubs* collimator_tubs = new G4Tubs("collimator", 0. * mm, 300. * mm, 170. * mm, 0. * deg, 360. * deg);
		collimator_logical = new G4LogicalVolume(collimator_tubs, world_material, "collimator", 0, 0, 0);
		m_collimator = new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, collimator_position_z), 
			collimator_logical, "collimator", gantry_logical, false, 0);
	}
	const bool build_photon = true;
	if (build_photon) {
		buildTarget(gantry_logical, gantry_position_z);
		buildPrimaryCollimator(gantry_logical, gantry_position_z);
        buildUpPlate(gantry_logical, gantry_position_z);
        buildMCBackscatterPlate(gantry_logical, gantry_position_z);
        buildSecondaryCollimator(collimator_logical, collimator_position_z, gantry_position_z);
        buildAlPlate(collimator_logical, collimator_position_z, gantry_position_z);
        buildPlasticMLCBox(collimator_logical, collimator_position_z, gantry_position_z);
        buildWindow(gantry_logical, gantry_position_z);
	}
	
	m_monitor_chamber = buildIonisationChamber(gantry_logical, gantry_position_z);
    fs::path mlc_proximal_path = m_stl_path / "HALCYON-MLC-PROXIMAL" / "binary";
    m_mlc_proximal_volumes = buildProximalLeafBank(collimator_logical, collimator_position_z, gantry_position_z, mlc_proximal_path);
	fs::path mlc_distal_path = m_stl_path / "HALCYON-MLC-DISTAL" / "binary";
    m_mlc_distal_volumes = buildDistalLeafBank(collimator_logical, collimator_position_z, gantry_position_z, mlc_distal_path);

	m_id_to_traversed = generateVolumeToTraversed(HalcyonTable, sizeof(HalcyonTable) / sizeof(HalcyonTable[0]), world_logical);
	return world_physical;
}

void TreatmentHeadDetector::setState(const Plan& plan, size_t pt_idx)
{
	G4GeometryManager::GetInstance()->OpenGeometry(m_gantry);
	const double gantry_position_z = 690.0 * mm;
    const auto& state = plan.getDualLayerMLC(TreatmentHeadType::HALCYON, pt_idx);
	auto gantry_transform =  G4RotateY3D(state.th_rotations.gantry_angle)*G4TranslateZ3D(gantry_position_z);
	m_gantry->SetTranslation(gantry_transform.getTranslation());
	m_gantry->SetRotation(new G4RotationMatrix(gantry_transform.getRotation().inverse()));
	auto collimator_transform = G4RotateZ3D(state.th_rotations.collimator_angle);
	m_collimator->SetRotation(new G4RotationMatrix(collimator_transform.getRotation().inverse())); //Could just use -angle
    SetProximalPositions(m_mlc_proximal_volumes, state.proximal_leaf_positions);
    SetDistalPositions(m_mlc_distal_volumes, state.distal_leaf_positions);
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