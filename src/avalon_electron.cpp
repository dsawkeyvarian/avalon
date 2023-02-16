//AVALON MACHINE

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
#include "external/CADMesh.hh"
#include "G4GDMLParser.hh"
#include "G4GenericPolycone.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "treatment_heads.h"
#include "plan.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include "math_utils.h"

static const double g_SAD = 100.0 * cm;

namespace avalon_electron {

enum class Target
{
	LOW_ENERGY,
  MEDIUM_ENERGY,
	HIGH_ENERGY,
	NUM_TARGETS
};
enum class FlatteningFilter
{
    OPEN,
  	ENERGY_6X,
    ENERGY_10X,
    NUM_FLATTENING_FILTERS
};

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
    G4Element* elTi = NISTman->FindOrBuildElement("Ti");
    G4Element* elZn = NISTman->FindOrBuildElement("Zn");

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
    //Al6061-T651
    {
        const double density = 2.7 * g / cm3;
		G4Material* Al6061 = new G4Material("Al6061-T651", density, 9);
		Al6061->AddElement(elAl, 97.2 * perCent); // 95.8-98.6
		Al6061->AddElement(elMg, 1.0 * perCent); // 0.8-1.2
		Al6061->AddElement(elMn, 0.1 * perCent); // max 0.15
		Al6061->AddElement(elCu, 0.3 * perCent); //0.15-0.4
        Al6061->AddElement(elCr, 0.2 * perCent); //0.04-0.35
        Al6061->AddElement(elSi, 0.6 * perCent); //0.4-0.8
        Al6061->AddElement(elFe, 0.3 * perCent); //max 0.7
        Al6061->AddElement(elTi, 0.1 * perCent); //max 0.15
        Al6061->AddElement(elZn, 0.2 * perCent); //max 0.25
    }
    //1010-1020HRS
    {
		const double density = 7.86112 * g / cm3;
		G4Material* SS_A36 = new G4Material("1010-1020HRS", density, 5);
		SS_A36->AddElement(elC, 0.2 * perCent); //0.17 - 0.230 %
		SS_A36->AddElement(elFe, 99.305 * perCent); //99.08 - 99.53 %
		SS_A36->AddElement(elMn, 0.405 * perCent); //0.30 - 0.60 %
		SS_A36->AddElement(elP, 0.04 * perCent); //≤ 0.040 %
		SS_A36->AddElement(elS, 0.05 * perCent); //≤ 0.050 %
	}	

}


MLCLeafVolumesProximal buildProximalLeafBank(G4LogicalVolume* parent_logical, const double parent_world_z, double parent_gantry_z, const std::filesystem::path& stl_folder) 
{
    //Visualization
	const auto bank_X1_vis = new G4VisAttributes(G4Colour(0.3, 0.3, 0.8, 0.9));
	const auto bank_X2_vis = new G4VisAttributes(G4Colour(0.3, 0.3, 0.8, 0.9));
    const auto outerleaves_vis = new G4VisAttributes(G4Colour(0.3, 0.8, 0.3, 0.9));

    //Order of files is meaningful so that leaf position indices are correct
    static const char* mlc_proximal_bank_X2[] = 
    {
        "p1040242008_p01a_ruby.stl",
        "p1040242006_p02a_ruby.stl",
        "p1040242007_p03a_ruby.stl",
        "p1040242005_p04a_ruby.stl",
        "p1040242008_p05a_ruby.stl",
        "p1040242006_p06a_ruby.stl",
        "p1040242007_p07a_ruby.stl",
        "p1040242005_p08a_ruby.stl",
        "p1040242008_p09a_ruby.stl",
        "p1040242006_p10a_ruby.stl",
        "p1040242007_p11a_ruby.stl",
        "p1040239002_p12a_ruby.stl",
        "p1040142009_p13a_ruby.stl",
        "p1040142007_p14a_ruby.stl",
        "p1040142010_p15a_ruby.stl",
        "p1040142008_p16a_ruby.stl",
        "p1040142006_p17a_ruby.stl",
        "p1040142009_p18a_ruby.stl",
        "p1040142007_p19a_ruby.stl",
        "p1040142010_p20a_ruby.stl",
        "p1040142008_p21a_ruby.stl",
        "p1040142006_p22a_ruby.stl",
        "p1040142009_p23a_ruby.stl",
        "p1040142007_p24a_ruby.stl",
        "p1040142010_p25a_ruby.stl",
        "p1040142008_p26a_ruby.stl",
        "p1040142006_p27a_ruby.stl",
        "p1040142009_p28a_ruby.stl",
        "p1040142007_p29a_ruby.stl",
        "p1040142010_p30a_ruby.stl",
        "p1040142008_p31a_ruby.stl",
        "p1040142006_p32a_ruby.stl",
        "p1040142009_p33a_ruby.stl",
        "p1040142007_p34a_ruby.stl",
        "p1040142010_p35a_ruby.stl",
        "p1040239002_p36a_ruby.stl",
        "p1040242007_p37a_ruby.stl",
        "p1040242006_p38a_ruby.stl",
        "p1040242008_p39a_ruby.stl",
        "p1040242005_p40a_ruby.stl",
        "p1040242007_p41a_ruby.stl",
        "p1040242006_p42a_ruby.stl",
        "p1040242008_p43a_ruby.stl",
        "p1040242005_p44a_ruby.stl",
        "p1040242007_p45a_ruby.stl",
        "p1040242006_p46a_ruby.stl",
        "p1040242007_p47a_ruby.stl"
    };
    static const char* mlc_proximal_bank_X1[] = 
    {
        "p1040242008_p01b_ruby.stl",
        "p1040242006_p02b_ruby.stl",
        "p1040242007_p03b_ruby.stl",
        "p1040242005_p04b_ruby.stl",
        "p1040242008_p05b_ruby.stl",
        "p1040242006_p06b_ruby.stl",
        "p1040242007_p07b_ruby.stl",
        "p1040242005_p08b_ruby.stl",
        "p1040242008_p09b_ruby.stl",
        "p1040242006_p10b_ruby.stl",
        "p1040242007_p11b_ruby.stl",
        "p1040239002_p12b_ruby.stl",
        "p1040142010_p13b_ruby.stl",
        "p1040142007_p14b_ruby.stl",
        "p1040142009_p15b_ruby.stl",
        "p1040142006_p16b_ruby.stl",
        "p1040142008_p17b_ruby.stl",
        "p1040142010_p18b_ruby.stl",
        "p1040142007_p19b_ruby.stl",
        "p1040142009_p20b_ruby.stl",
        "p1040142006_p21b_ruby.stl",
        "p1040142008_p22b_ruby.stl",
        "p1040142010_p23b_ruby.stl",
        "p1040142007_p24b_ruby.stl",
        "p1040142009_p25b_ruby.stl",
        "p1040142006_p26b_ruby.stl",
        "p1040142008_p27b_ruby.stl",
        "p1040142010_p28b_ruby.stl",
        "p1040142007_p29b_ruby.stl",
        "p1040142009_p30b_ruby.stl",
        "p1040142006_p31b_ruby.stl",
        "p1040142008_p32b_ruby.stl",
        "p1040142010_p33b_ruby.stl",
        "p1040142007_p34b_ruby.stl",
        "p1040142009_p35b_ruby.stl",
        "p1040239002_p36b_ruby.stl",
        "p1040242007_p37b_ruby.stl",
        "p1040242006_p38b_ruby.stl",
        "p1040242008_p39b_ruby.stl",
        "p1040242005_p40b_ruby.stl",
        "p1040242007_p41b_ruby.stl",
        "p1040242006_p42b_ruby.stl",
        "p1040242008_p43b_ruby.stl",
        "p1040242005_p44b_ruby.stl",
        "p1040242007_p45b_ruby.stl",
        "p1040242006_p46b_ruby.stl",
        "p1040242008_p47b_ruby.stl"
    };

    static const char* mlc_proximal_outerleaves[] = 
    {
        "p1047532002_bdy-leaf_a_ava.stl",
        "p1047532002_bdy-leaf_b_ava.stl"
    };

    auto mlc_volumes_proximal = MLCLeafVolumesProximal();
    // BANK A
    {
        size_t num_bank_X1 = sizeof(mlc_proximal_bank_X1) / sizeof(*mlc_proximal_bank_X1);
        for (size_t i = 0; i < num_bank_X1; ++i) {
		    std::filesystem::path filepath = stl_folder / "Helsinki-RAD_COLL_RUBY-MLC"/ mlc_proximal_bank_X1[i];
		    auto leaf_mesh = CADMesh::TessellatedMesh::FromSTL(filepath.u8string());
            leaf_mesh->SetOffset(4.8146, 0.0, 0.0);
            G4VSolid *leaf_solid = leaf_mesh->GetSolid();
            G4LogicalVolume *leaf_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("W95"), "mlc_proximal_leaf", 0, 0, 0);
            leaf_logical->SetVisAttributes(bank_X1_vis);
            auto volume = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_gantry_z - parent_world_z), leaf_logical, "mlc_proximal_leaf", parent_logical, false, 0);
            mlc_volumes_proximal.bank_X1.push_back(volume);
        }
    }
    //BANK B
    {
        size_t num_bank_X2 = sizeof(mlc_proximal_bank_X2) / sizeof(*mlc_proximal_bank_X2);
        for (size_t i = 0; i < num_bank_X2; ++i){  
            std::filesystem::path filepath = stl_folder/ "Helsinki-RAD_COLL_RUBY-MLC" / mlc_proximal_bank_X2[i];
            auto leaf_mesh = CADMesh::TessellatedMesh::FromSTL(filepath.u8string());
            leaf_mesh->SetOffset(-88.2507, 0.0, 0.0);
            G4VSolid *leaf_solid = leaf_mesh->GetSolid();
            G4LogicalVolume *leaf_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("W95"), "mlc_proximal_leaf", 0, 0, 0);
            leaf_logical->SetVisAttributes(bank_X2_vis);
            auto volume = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_gantry_z - parent_world_z), leaf_logical, "mlc_proximal_leaf", parent_logical, false, 0);
            mlc_volumes_proximal.bank_X2.push_back(volume);
        }
    }

    {
        size_t num_outerleaves = sizeof(mlc_proximal_outerleaves) / sizeof(*mlc_proximal_outerleaves);
        for (size_t i = 0; i < num_outerleaves; ++i) {
		    std::filesystem::path filepath = stl_folder / "MLC" / mlc_proximal_outerleaves[i];
            auto leaf_mesh = CADMesh::TessellatedMesh::FromSTL(filepath.u8string());
            G4VSolid *leaf_solid = leaf_mesh->GetSolid();
            G4LogicalVolume *leaf_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("W95"), "mlc_proximal_outer_leaf", 0, 0, 0);
            leaf_logical->SetVisAttributes(outerleaves_vis);
            auto volume = new G4PVPlacement(0, G4ThreeVector(0., 0.,  g_SAD - parent_gantry_z - parent_world_z), leaf_logical, "mlc_proximal_outer_leaf", parent_logical, false, 0);
        }
    }
    return mlc_volumes_proximal;
}

void SetProximalPositions(const MLCLeafVolumesProximal& volumes, const LeafPositions& positions) 
{
    const size_t num_leaves_in_bank = 47;
    const bool valid = (volumes.bank_X2.size() == num_leaves_in_bank) || (volumes.bank_X1.size() == num_leaves_in_bank) 
    || (positions.bank_X2.size() == num_leaves_in_bank) || (positions.bank_X1.size() == num_leaves_in_bank);
    if (!valid)
        throw std::runtime_error("Number of Proximal leaf volumes and leaf positions is inconsistent");
    
    const std::vector<double>* bank_positions[] = { &positions.bank_X1, &positions.bank_X2 };
    const std::vector<G4VPhysicalVolume*>* bank_volumes[] = { &volumes.bank_X1, &volumes.bank_X2};
    G4double physical_position = 0.0*mm;
    const int num_proximal_pos = 401;
    G4double mlc_proximal_nominal[num_proximal_pos] = 
		{  -200.0*mm, -199.0*mm, -198.0*mm, -197.0*mm, -196.0*mm, -195.0*mm, -194.0*mm, -193.0*mm, -192.0*mm, -191.0*mm, -190.0*mm,
        -189.0*mm, -188.0*mm, -187.0*mm, -186.0*mm, -185.0*mm, -184.0*mm, -183.0*mm, -182.0*mm, -181.0*mm, -180.0*mm, -179.0*mm, 
        -178.0*mm, -177.0*mm, -176.0*mm, -175.0*mm, -174.0*mm, -173.0*mm, -172.0*mm, -171.0*mm, -170.0*mm, -169.0*mm, -168.0*mm,
        -167.0*mm, -166.0*mm, -165.0*mm, -164.0*mm, -163.0*mm, -162.0*mm, -161.0*mm, -160.0*mm, -159.0*mm, -158.0*mm, -157.0*mm,
        -156.0*mm, -155.0*mm, -154.0*mm, -153.0*mm, -152.0*mm, -151.0*mm, -150.0*mm, -149.0*mm, -148.0*mm, -147.0*mm, -146.0*mm, 
        -145.0*mm, -144.0*mm, -143.0*mm, -142.0*mm, -141.0*mm, -140.0*mm, -139.0*mm, -138.0*mm, -137.0*mm, -136.0*mm, -135.0*mm, 
        -134.0*mm, -133.0*mm, -132.0*mm, -131.0*mm, -130.0*mm, -129.0*mm, -128.0*mm, -127.0*mm, -126.0*mm, -125.0*mm, -124.0*mm, 
        -123.0*mm, -122.0*mm, -121.0*mm, -120.0*mm, -119.0*mm, -118.0*mm, -117.0*mm, -116.0*mm, -115.0*mm, -114.0*mm, -113.0*mm, 
        -112.0*mm, -111.0*mm, -110.0*mm, -109.0*mm, -108.0*mm, -107.0*mm, -106.0*mm, -105.0*mm, -104.0*mm, -103.0*mm, -102.0*mm, 
        -101.0*mm, -100.0*mm, -99.0*mm, -98.0*mm, -97.0*mm, -96.0*mm, -95.0*mm, -94.0*mm, -93.0*mm, -92.0*mm, -91.0*mm, -90.0*mm, 
        -89.0*mm, -88.0*mm, -87.0*mm, -86.0*mm, -85.0*mm, -84.0*mm, -83.0*mm, -82.0*mm, -81.0*mm, -80.0*mm, -79.0*mm, -78.0*mm, 
        -77.0*mm, -76.0*mm, -75.0*mm, -74.0*mm, -73.0*mm, -72.0*mm, -71.0*mm, -70.0*mm, -69.0*mm, -68.0*mm, -67.0*mm, -66.0*mm, 
        -65.0*mm, -64.0*mm, -63.0*mm, -62.0*mm, -61.0*mm, -60.0*mm, -59.0*mm, -58.0*mm, -57.0*mm, -56.0*mm, -55.0*mm, -54.0*mm, 
        -53.0*mm, -52.0*mm, -51.0*mm, -50.0*mm, -49.0*mm, -48.0*mm, -47.0*mm, -46.0*mm, -45.0*mm, -44.0*mm, -43.0*mm, -42.0*mm, 
        -41.0*mm, -40.0*mm, -39.0*mm, -38.0*mm, -37.0*mm, -36.0*mm, -35.0*mm, -34.0*mm, -33.0*mm, -32.0*mm, -31.0*mm, -30.0*mm, 
        -29.0*mm, -28.0*mm, -27.0*mm, -26.0*mm, -25.0*mm, -24.0*mm, -23.0*mm, -22.0*mm, -21.0*mm, -20.0*mm, -19.0*mm, -18.0*mm, 
        -17.0*mm, -16.0*mm, -15.0*mm, -14.0*mm, -13.0*mm, -12.0*mm, -11.0*mm, -10.0*mm, -9.0*mm, -8.0*mm, -7.0*mm, -6.0*mm, 
        -5.0*mm, -4.0*mm, -3.0*mm, -2.0*mm, -1.0*mm, 0.0*mm, 1.0*mm, 2.0*mm, 3.0*mm, 4.0*mm, 5.0*mm, 6.0*mm, 7.0*mm, 8.0*mm, 
        9.0*mm, 10.0*mm, 11.0*mm, 12.0*mm, 13.0*mm, 14.0*mm, 15.0*mm, 16.0*mm, 17.0*mm, 18.0*mm, 19.0*mm, 20.0*mm, 21.0*mm, 
        22.0*mm, 23.0*mm, 24.0*mm, 25.0*mm, 26.0*mm, 27.0*mm, 28.0*mm, 29.0*mm, 30.0*mm, 31.0*mm, 32.0*mm, 33.0*mm, 34.0*mm, 
        35.0*mm, 36.0*mm, 37.0*mm, 38.0*mm, 39.0*mm, 40.0*mm, 41.0*mm, 42.0*mm, 43.0*mm, 44.0*mm, 45.0*mm, 46.0*mm, 47.0*mm, 
        48.0*mm, 49.0*mm, 50.0*mm, 51.0*mm, 52.0*mm, 53.0*mm, 54.0*mm, 55.0*mm, 56.0*mm, 57.0*mm, 58.0*mm, 59.0*mm, 60.0*mm, 
        61.0*mm, 62.0*mm, 63.0*mm, 64.0*mm, 65.0*mm, 66.0*mm, 67.0*mm, 68.0*mm, 69.0*mm, 70.0*mm, 71.0*mm, 72.0*mm, 73.0*mm, 
        74.0*mm, 75.0*mm, 76.0*mm, 77.0*mm, 78.0*mm, 79.0*mm, 80.0*mm, 81.0*mm, 82.0*mm, 83.0*mm, 84.0*mm, 85.0*mm, 86.0*mm, 
        87.0*mm, 88.0*mm, 89.0*mm, 90.0*mm, 91.0*mm, 92.0*mm, 93.0*mm, 94.0*mm, 95.0*mm, 96.0*mm, 97.0*mm, 98.0*mm, 99.0*mm, 
        100.0*mm, 101.0*mm, 102.0*mm, 103.0*mm, 104.0*mm, 105.0*mm, 106.0*mm, 107.0*mm, 108.0*mm, 109.0*mm, 110.0*mm, 111.0*mm, 
        112.0*mm, 113.0*mm, 114.0*mm, 115.0*mm, 116.0*mm, 117.0*mm, 118.0*mm, 119.0*mm, 120.0*mm, 121.0*mm, 122.0*mm, 123.0*mm, 
        124.0*mm, 125.0*mm, 126.0*mm, 127.0*mm, 128.0*mm, 129.0*mm, 130.0*mm, 131.0*mm, 132.0*mm, 133.0*mm, 134.0*mm, 135.0*mm, 
        136.0*mm, 137.0*mm, 138.0*mm, 139.0*mm, 140.0*mm, 141.0*mm, 142.0*mm, 143.0*mm, 144.0*mm, 145.0*mm, 146.0*mm, 147.0*mm, 
        148.0*mm, 149.0*mm, 150.0*mm, 151.0*mm, 152.0*mm, 153.0*mm, 154.0*mm, 155.0*mm, 156.0*mm, 157.0*mm, 158.0*mm, 159.0*mm, 
        160.0*mm, 161.0*mm, 162.0*mm, 163.0*mm, 164.0*mm, 165.0*mm, 166.0*mm, 167.0*mm, 168.0*mm, 169.0*mm, 170.0*mm, 171.0*mm, 
        172.0*mm, 173.0*mm, 174.0*mm, 175.0*mm, 176.0*mm, 177.0*mm, 178.0*mm, 179.0*mm, 180.0*mm, 181.0*mm, 182.0*mm, 183.0*mm, 
        184.0*mm, 185.0*mm, 186.0*mm, 187.0*mm, 188.0*mm, 189.0*mm, 190.0*mm, 191.0*mm, 192.0*mm, 193.0*mm, 194.0*mm, 195.0*mm, 
        196.0*mm, 197.0*mm, 198.0*mm, 199.0*mm, 200.0*mm};
	G4double mlc_proximal_actual[num_proximal_pos] =
		{-72.142958*mm,
        -71.770646*mm,
        -71.398334*mm,
        -71.026022*mm,
        -70.65371*mm,
        -70.281399*mm,
        -69.909087*mm,
        -69.536775*mm,
        -69.164463*mm,
        -68.792151*mm,
        -68.419839*mm,
        -68.047527*mm,
        -67.675215*mm,
        -67.302903*mm,
        -66.930592*mm,
        -66.55828*mm,
        -66.185968*mm,
        -65.813656*mm,
        -65.441344*mm,
        -65.069627*mm,
        -64.70162*mm,
        -64.333613*mm,
        -63.965606*mm,
        -63.597598*mm,
        -63.229591*mm,
        -62.861584*mm,
        -62.493577*mm,
        -62.12557*mm,
        -61.757563*mm,
        -61.389556*mm,
        -61.021549*mm,
        -60.653542*mm,
        -60.285534*mm,
        -59.917527*mm,
        -59.54952*mm,
        -59.181513*mm,
        -58.813506*mm,
        -58.445499*mm,
        -58.077492*mm,
        -57.709485*mm,
        -57.341478*mm,
        -56.97347*mm,
        -56.605463*mm,
        -56.237456*mm,
        -55.869449*mm,
        -55.501442*mm,
        -55.133435*mm,
        -54.765428*mm,
        -54.397421*mm,
        -54.029414*mm,
        -53.661407*mm,
        -53.293399*mm,
        -52.925498*mm,
        -52.560415*mm,
        -52.195332*mm,
        -51.830249*mm,
        -51.465166*mm,
        -51.100083*mm,
        -50.735*mm,
        -50.369917*mm,
        -50.004835*mm,
        -49.639752*mm,
        -49.274669*mm,
        -48.909586*mm,
        -48.544503*mm,
        -48.17942*mm,
        -47.814337*mm,
        -47.449254*mm,
        -47.084171*mm,
        -46.719089*mm,
        -46.354006*mm,
        -45.988923*mm,
        -45.62384*mm,
        -45.259282*mm,
        -44.895567*mm,
        -44.531851*mm,
        -44.168136*mm,
        -43.80442*mm,
        -43.440704*mm,
        -43.076989*mm,
        -42.713273*mm,
        -42.349558*mm,
        -41.986155*mm,
        -41.62404*mm,
        -41.261925*mm,
        -40.89981*mm,
        -40.537694*mm,
        -40.175579*mm,
        -39.813464*mm,
        -39.451348*mm,
        -39.089233*mm,
        -38.727118*mm,
        -38.365002*mm,
        -38.002887*mm,
        -37.640772*mm,
        -37.278657*mm,
        -36.916541*mm,
        -36.554426*mm,
        -36.192311*mm,
        -35.830195*mm,
        -35.46808*mm,
        -35.105965*mm,
        -34.743849*mm,
        -34.383406*mm,
        -34.024302*mm,
        -33.665198*mm,
        -33.306094*mm,
        -32.94699*mm,
        -32.587885*mm,
        -32.228781*mm,
        -31.869677*mm,
        -31.510573*mm,
        -31.151469*mm,
        -30.792365*mm,
        -30.43326*mm,
        -30.074156*mm,
        -29.715052*mm,
        -29.355948*mm,
        -28.996844*mm,
        -28.637739*mm,
        -28.278635*mm,
        -27.919531*mm,
        -27.560427*mm,
        -27.201323*mm,
        -26.842218*mm,
        -26.483114*mm,
        -26.12401*mm,
        -25.764906*mm,
        -25.406304*mm,
        -25.050662*mm,
        -24.695021*mm,
        -24.33938*mm,
        -23.983738*mm,
        -23.628097*mm,
        -23.272456*mm,
        -22.916814*mm,
        -22.561173*mm,
        -22.205532*mm,
        -21.849891*mm,
        -21.494249*mm,
        -21.138608*mm,
        -20.782967*mm,
        -20.427325*mm,
        -20.071684*mm,
        -19.716043*mm,
        -19.360401*mm,
        -19.00476*mm,
        -18.649119*mm,
        -18.293478*mm,
        -17.937836*mm,
        -17.582195*mm,
        -17.226554*mm,
        -16.870912*mm,
        -16.515271*mm,
        -16.15963*mm,
        -15.803988*mm,
        -15.448347*mm,
        -15.092706*mm,
        -14.738221*mm,
        -14.386121*mm,
        -14.034022*mm,
        -13.681922*mm,
        -13.329822*mm,
        -12.977722*mm,
        -12.625622*mm,
        -12.273523*mm,
        -11.921423*mm,
        -11.569323*mm,
        -11.217223*mm,
        -10.865123*mm,
        -10.513024*mm,
        -10.160924*mm,
        -9.808824*mm,
        -9.456724*mm,
        -9.104625*mm,
        -8.752525*mm,
        -8.400425*mm,
        -8.048325*mm,
        -7.696225*mm,
        -7.344126*mm,
        -6.992026*mm,
        -6.639926*mm,
        -6.287826*mm,
        -5.935726*mm,
        -5.583627*mm,
        -5.231527*mm,
        -4.879427*mm,
        -4.530437*mm,
        -4.181916*mm,
        -3.833394*mm,
        -3.484873*mm,
        -3.136352*mm,
        -2.78783*mm,
        -2.439309*mm,
        -2.090788*mm,
        -1.742267*mm,
        -1.393814*mm,
        -1.04536*mm,
        -0.696907*mm,
        -0.348453*mm,
        0.0*mm,
        0.348453*mm,
        0.696907*mm,
        1.04536*mm,
        1.393814*mm,
        1.742267*mm,
        2.09072*mm,
        2.439174*mm,
        2.787627*mm,
        3.136081*mm,
        3.484534*mm,
        3.832987*mm,
        4.179157*mm,
        4.52391*mm,
        4.868662*mm,
        5.213414*mm,
        5.558167*mm,
        5.902919*mm,
        6.247671*mm,
        6.592424*mm,
        6.937176*mm,
        7.281928*mm,
        7.62668*mm,
        7.971433*mm,
        8.316185*mm,
        8.660937*mm,
        9.00569*mm,
        9.350442*mm,
        9.695194*mm,
        10.039946*mm,
        10.384699*mm,
        10.729451*mm,
        11.074203*mm,
        11.418956*mm,
        11.763708*mm,
        12.10846*mm,
        12.453212*mm,
        12.797965*mm,
        13.142717*mm,
        13.487401*mm,
        13.828284*mm,
        14.169168*mm,
        14.510052*mm,
        14.850935*mm,
        15.191819*mm,
        15.532702*mm,
        15.873586*mm,
        16.21447*mm,
        16.555353*mm,
        16.896237*mm,
        17.23712*mm,
        17.578004*mm,
        17.918888*mm,
        18.259771*mm,
        18.600655*mm,
        18.941539*mm,
        19.282422*mm,
        19.623306*mm,
        19.964189*mm,
        20.305073*mm,
        20.645957*mm,
        20.98684*mm,
        21.327724*mm,
        21.668607*mm,
        22.009491*mm,
        22.347277*mm,
        22.684163*mm,
        23.021049*mm,
        23.357934*mm,
        23.69482*mm,
        24.031706*mm,
        24.368592*mm,
        24.705478*mm,
        25.042363*mm,
        25.379249*mm,
        25.716135*mm,
        26.053021*mm,
        26.389907*mm,
        26.726792*mm,
        27.063678*mm,
        27.400564*mm,
        27.73745*mm,
        28.074336*mm,
        28.411221*mm,
        28.748107*mm,
        29.084993*mm,
        29.421879*mm,
        29.758765*mm,
        30.09565*mm,
        30.430469*mm,
        30.764062*mm,
        31.097655*mm,
        31.431248*mm,
        31.76484*mm,
        32.098433*mm,
        32.432026*mm,
        32.765619*mm,
        33.099212*mm,
        33.432804*mm,
        33.766397*mm,
        34.09999*mm,
        34.433583*mm,
        34.767176*mm,
        35.100768*mm,
        35.434361*mm,
        35.767954*mm,
        36.101547*mm,
        36.43514*mm,
        36.768732*mm,
        37.102325*mm,
        37.435918*mm,
        37.769503*mm,
        38.099673*mm,
        38.429842*mm,
        38.760011*mm,
        39.090181*mm,
        39.42035*mm,
        39.750519*mm,
        40.080689*mm,
        40.410858*mm,
        40.741027*mm,
        41.071197*mm,
        41.401366*mm,
        41.731535*mm,
        42.061704*mm,
        42.391874*mm,
        42.722043*mm,
        43.052212*mm,
        43.382382*mm,
        43.712551*mm,
        44.04272*mm,
        44.37289*mm,
        44.702769*mm,
        45.029385*mm,
        45.356*mm,
        45.682616*mm,
        46.009232*mm,
        46.335847*mm,
        46.662463*mm,
        46.989078*mm,
        47.315694*mm,
        47.64231*mm,
        47.968925*mm,
        48.295541*mm,
        48.622156*mm,
        48.948772*mm,
        49.275388*mm,
        49.602003*mm,
        49.928619*mm,
        50.255234*mm,
        50.58185*mm,
        50.90845*mm,
        51.231764*mm,
        51.554772*mm,
        51.877705*mm,
        52.200639*mm,
        52.523572*mm,
        52.846505*mm,
        53.169439*mm,
        53.492372*mm,
        53.815306*mm,
        54.138239*mm,
        54.461172*mm,
        54.784106*mm,
        55.107039*mm,
        55.429973*mm,
        55.752906*mm,
        56.075839*mm,
        56.398773*mm,
        56.719246*mm,
        57.038973*mm,
        57.3587*mm,
        57.678427*mm,
        57.998154*mm,
        58.317881*mm,
        58.637607*mm,
        58.957334*mm,
        59.277061*mm,
        59.596282*mm,
        59.9154*mm,
        60.234518*mm,
        60.553636*mm,
        60.872753*mm,
        61.191871*mm,
        61.510989*mm,
        61.82983*mm,
        62.145895*mm,
        62.461961*mm,
        62.778026*mm,
        63.094091*mm,
        63.410156*mm,
        63.726221*mm,
        64.042286*mm,
        64.358352*mm,
        64.674417*mm,
        64.990482*mm,
        65.306547*mm,
        65.622088*mm,
        65.937251*mm,
        66.252415*mm,
        66.567578*mm};
    

    for (size_t j = 0; j < 2; ++j) {
        for (size_t i = 0; i < num_leaves_in_bank; ++i) {
            double nominal_position = bank_positions[j]->at(i);
            if (j == 1) {nominal_position *= (-1);}
            for (size_t k = 0; k < num_proximal_pos; ++k) {
                if (std::abs(nominal_position - mlc_proximal_nominal[k]) <= 0.5) {
					physical_position = (nominal_position - mlc_proximal_nominal[k]) / (mlc_proximal_nominal[k+1] - mlc_proximal_nominal[k]) * (mlc_proximal_actual[k+1] - mlc_proximal_actual[k]) + mlc_proximal_actual[k];
                    break;
                }
            }
            if (j == 1) {physical_position *= (-1);}
            auto volume = bank_volumes[j]->at(i);
            auto translation = volume->GetTranslation();
            translation.setX(physical_position);
            volume->SetTranslation(translation);
        }
    }
}

MLCLeafVolumesDistal buildDistalLeafBank(G4LogicalVolume* parent_logical, const double parent_world_z, double parent_gantry_z, const std::filesystem::path& stl_folder) 
{
    //Visualization
	const auto bank_X1_vis = new G4VisAttributes(G4Colour(0.0, 0.8, 0.8, 0.4));
	const auto bank_X2_vis = new G4VisAttributes(G4Colour(0.0, 0.8, 0.8, 0.4));
    const auto outerleaves_vis = new G4VisAttributes(G4Colour(0.3, 0.8, 0.3, 0.9));

    //Order of files is meaningful so that leaf position indices are correct
    static const char* mlc_distal_bank_X2[] = 
    {"p1040250005_d01a_ruby.stl",
    "p1040250008_d02a_ruby.stl",
    "p1040250006_d03a_ruby.stl",
    "p1040250007_d04a_ruby.stl",
    "p1040250005_d05a_ruby.stl",
    "p1040250008_d06a_ruby.stl",
    "p1040250006_d07a_ruby.stl",
    "p1040250007_d08a_ruby.stl",
    "p1040250005_d09a_ruby.stl",
    "p1040250008_d10a_ruby.stl",
    "p1040250006_d11a_ruby.stl",
    "p1040247007_d12a_ruby.stl",
    "p1040247005_d13a_ruby.stl",
    "p1040247008_d14a_ruby.stl",
    "p1040247006_d15a_ruby.stl",
    "p1040247007_d16a_ruby.stl",
    "p1040247005_d17a_ruby.stl",
    "p1040247008_d18a_ruby.stl",
    "p1040247006_d19a_ruby.stl",
    "p1040247007_d20a_ruby.stl",
    "p1040247005_d21a_ruby.stl",
    "p1040247008_d22a_ruby.stl",
    "p1040247006_d23a_ruby.stl",
    "p1040247007_d24a_ruby.stl",
    "p1040247005_d25a_ruby.stl",
    "p1040247008_d26a_ruby.stl",
    "p1040247006_d27a_ruby.stl",
    "p1040247007_d28a_ruby.stl",
    "p1040247005_d29a_ruby.stl",
    "p1040247008_d30a_ruby.stl",
    "p1040247006_d31a_ruby.stl",
    "p1040247007_d32a_ruby.stl",
    "p1040247005_d33a_ruby.stl",
    "p1040247008_d34a_ruby.stl",
    "p1040247006_d35a_ruby.stl",
    "p1040250007_d36a_ruby.stl",
    "p1040250005_d37a_ruby.stl",
    "p1040250008_d38a_ruby.stl",
    "p1040250006_d39a_ruby.stl",
    "p1040250007_d40a_ruby.stl",
    "p1040250005_d41a_ruby.stl",
    "p1040250008_d42a_ruby.stl",
    "p1040250006_d43a_ruby.stl",
    "p1040250007_d44a_ruby.stl",
    "p1040250005_d45a_ruby.stl",
    "p1040250008_d46a_ruby.stl"};
    static const char* mlc_distal_bank_X1[] = 
    {"p1040250008_d01b_ruby.stl",
    "p1040250005_d02b_ruby.stl",
    "p1040250007_d03b_ruby.stl",
    "p1040250006_d04b_ruby.stl",
    "p1040250008_d05b_ruby.stl",
    "p1040250005_d06b_ruby.stl",
    "p1040250007_d07b_ruby.stl",
    "p1040250006_d08b_ruby.stl",
    "p1040250008_d09b_ruby.stl",
    "p1040250005_d10b_ruby.stl",
    "p1040250007_d11b_ruby.stl",
    "p1040247006_d12b_ruby.stl",
    "p1040247008_d13b_ruby.stl",
    "p1040247005_d14b_ruby.stl",
    "p1040247007_d15b_ruby.stl",
    "p1040247006_d16b_ruby.stl",
    "p1040247008_d17b_ruby.stl",
    "p1040247005_d18b_ruby.stl",
    "p1040247007_d19b_ruby.stl",
    "p1040247006_d20b_ruby.stl",
    "p1040247008_d21b_ruby.stl",
    "p1040247005_d22b_ruby.stl",
    "p1040247007_d23b_ruby.stl",
    "p1040247006_d24b_ruby.stl",
    "p1040247008_d25b_ruby.stl",
    "p1040247005_d26b_ruby.stl",
    "p1040247007_d27b_ruby.stl",
    "p1040247006_d28b_ruby.stl",
    "p1040247008_d29b_ruby.stl",
    "p1040247005_d30b_ruby.stl",
    "p1040247007_d31b_ruby.stl",
    "p1040247006_d32b_ruby.stl",
    "p1040247008_d33b_ruby.stl",
    "p1040247005_d34b_ruby.stl",
    "p1040247007_d35b_ruby.stl",
    "p1040250006_d36b_ruby.stl",
    "p1040250008_d37b_ruby.stl",
    "p1040250005_d38b_ruby.stl",
    "p1040250007_d39b_ruby.stl",
    "p1040250006_d40b_ruby.stl",
    "p1040250008_d41b_ruby.stl",
    "p1040250005_d42b_ruby.stl",
    "p1040250007_d43b_ruby.stl",
    "p1040250006_d44b_ruby.stl",
    "p1040250008_d45b_ruby.stl",
    "p1040250005_d46b_ruby.stl"
    };
    
    static const char* mlc_distal_outerleaves[] = 
        {"p1040451002_bdy-leaf_a_ava.stl",
        "p1040451002_bdy-leaf_b_ava.stl"
        };

    auto mlc_volumes_distal = MLCLeafVolumesDistal();
    // BANK A
    {
        size_t num_bank_X1 = sizeof(mlc_distal_bank_X1) / sizeof(*mlc_distal_bank_X1);
        for (size_t i = 0; i < num_bank_X1; ++i) {
            std::filesystem::path filepath = stl_folder / "Helsinki-RAD_COLL_RUBY-MLC" / mlc_distal_bank_X1[i];
		    auto leaf_mesh = CADMesh::TessellatedMesh::FromSTL(filepath.u8string());
            leaf_mesh->SetOffset(12.394, 0.0, 0.0);
            G4VSolid *leaf_solid = leaf_mesh->GetSolid();
            G4LogicalVolume *leaf_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("W95"), "mlc_distal_leaf", 0, 0, 0);
            leaf_logical->SetVisAttributes(bank_X1_vis);
            auto volume = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_gantry_z - parent_world_z), leaf_logical, "mlc_distal_leaf", parent_logical, false, 0);
            mlc_volumes_distal.bank_X1.push_back(volume);
        }
    }
    // BANK B
    {
        size_t num_bank_X2 = sizeof(mlc_distal_bank_X2) / sizeof(*mlc_distal_bank_X2);
        for (size_t i = 0; i < num_bank_X2; ++i) {
		    std::filesystem::path filepath = stl_folder/ "Helsinki-RAD_COLL_RUBY-MLC" / mlc_distal_bank_X2[i];
            auto leaf_mesh = CADMesh::TessellatedMesh::FromSTL(filepath.u8string());
            leaf_mesh->SetOffset(-95.8306, 0.0, 0.0);
            G4VSolid *leaf_solid = leaf_mesh->GetSolid();
            G4LogicalVolume *leaf_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("W95"), "mlc_distal_leaf", 0, 0, 0);
            leaf_logical->SetVisAttributes(bank_X2_vis);
            auto volume = new G4PVPlacement(0, G4ThreeVector(0., 0.,  g_SAD - parent_gantry_z - parent_world_z), leaf_logical, "mlc_distal_leaf", parent_logical, false, 0);
            mlc_volumes_distal.bank_X2.push_back(volume);
        }
    }
    {
        size_t num_outerleaves = sizeof(mlc_distal_outerleaves) / sizeof(*mlc_distal_outerleaves);
        for (size_t i = 0; i < num_outerleaves; ++i) {
		    std::filesystem::path filepath = stl_folder / "MLC" /  mlc_distal_outerleaves[i];
            auto leaf_mesh = CADMesh::TessellatedMesh::FromSTL(filepath.u8string());
            G4VSolid *leaf_solid = leaf_mesh->GetSolid();
            G4LogicalVolume *leaf_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("W95"), "mlc_distal_outer_leaf", 0, 0, 0);
            leaf_logical->SetVisAttributes(outerleaves_vis);
            auto volume = new G4PVPlacement(0, G4ThreeVector(0., 0.,  g_SAD - parent_gantry_z - parent_world_z), leaf_logical, "mlc_distal_outer_leaf", parent_logical, false, 0);
        }
    }  

    return mlc_volumes_distal;
}

void SetDistalPositions(const MLCLeafVolumesDistal& volumes, const LeafPositions& positions) 
{
    const size_t num_leaves_in_bank = 46;
    const bool valid = (volumes.bank_X1.size() == num_leaves_in_bank) || (volumes.bank_X2.size() == num_leaves_in_bank) 
        || (positions.bank_X2.size() == num_leaves_in_bank) || (positions.bank_X1.size() == num_leaves_in_bank);
    if (!valid)
        throw std::runtime_error("Number of Distal leaf volumes and leaf positions is inconsistent");
    
    const std::vector<double>* bank_positions[] = { &positions.bank_X1, &positions.bank_X2 };
    const std::vector<G4VPhysicalVolume*>* bank_volumes[] = { &volumes.bank_X1, &volumes.bank_X2};
    G4double physical_position = 0.0*mm;

    const int num_distal_pos = 401;
    G4double mlc_distal_nominal[num_distal_pos] = 
		{-200.0*mm, -199.0*mm, -198.0*mm, -197.0*mm, -196.0*mm, -195.0*mm, -194.0*mm, -193.0*mm, -192.0*mm, -191.0*mm, -190.0*mm,
        -189.0*mm, -188.0*mm, -187.0*mm, -186.0*mm, -185.0*mm, -184.0*mm, -183.0*mm, -182.0*mm, -181.0*mm, -180.0*mm, -179.0*mm, 
        -178.0*mm, -177.0*mm, -176.0*mm, -175.0*mm, -174.0*mm, -173.0*mm, -172.0*mm, -171.0*mm, -170.0*mm, -169.0*mm, -168.0*mm,
        -167.0*mm, -166.0*mm, -165.0*mm, -164.0*mm, -163.0*mm, -162.0*mm, -161.0*mm, -160.0*mm, -159.0*mm, -158.0*mm, -157.0*mm,
        -156.0*mm, -155.0*mm, -154.0*mm, -153.0*mm, -152.0*mm, -151.0*mm, -150.0*mm, -149.0*mm, -148.0*mm, -147.0*mm, -146.0*mm, 
        -145.0*mm, -144.0*mm, -143.0*mm, -142.0*mm, -141.0*mm, -140.0*mm, -139.0*mm, -138.0*mm, -137.0*mm, -136.0*mm, -135.0*mm, 
        -134.0*mm, -133.0*mm, -132.0*mm, -131.0*mm, -130.0*mm, -129.0*mm, -128.0*mm, -127.0*mm, -126.0*mm, -125.0*mm, -124.0*mm, 
        -123.0*mm, -122.0*mm, -121.0*mm, -120.0*mm, -119.0*mm, -118.0*mm, -117.0*mm, -116.0*mm, -115.0*mm, -114.0*mm, -113.0*mm, 
        -112.0*mm, -111.0*mm, -110.0*mm, -109.0*mm, -108.0*mm, -107.0*mm, -106.0*mm, -105.0*mm, -104.0*mm, -103.0*mm, -102.0*mm, 
        -101.0*mm, -100.0*mm, -99.0*mm, -98.0*mm, -97.0*mm, -96.0*mm, -95.0*mm, -94.0*mm, -93.0*mm, -92.0*mm, -91.0*mm, -90.0*mm, 
        -89.0*mm, -88.0*mm, -87.0*mm, -86.0*mm, -85.0*mm, -84.0*mm, -83.0*mm, -82.0*mm, -81.0*mm, -80.0*mm, -79.0*mm, -78.0*mm, 
        -77.0*mm, -76.0*mm, -75.0*mm, -74.0*mm, -73.0*mm, -72.0*mm, -71.0*mm, -70.0*mm, -69.0*mm, -68.0*mm, -67.0*mm, -66.0*mm, 
        -65.0*mm, -64.0*mm, -63.0*mm, -62.0*mm, -61.0*mm, -60.0*mm, -59.0*mm, -58.0*mm, -57.0*mm, -56.0*mm, -55.0*mm, -54.0*mm, 
        -53.0*mm, -52.0*mm, -51.0*mm, -50.0*mm, -49.0*mm, -48.0*mm, -47.0*mm, -46.0*mm, -45.0*mm, -44.0*mm, -43.0*mm, -42.0*mm, 
        -41.0*mm, -40.0*mm, -39.0*mm, -38.0*mm, -37.0*mm, -36.0*mm, -35.0*mm, -34.0*mm, -33.0*mm, -32.0*mm, -31.0*mm, -30.0*mm, 
        -29.0*mm, -28.0*mm, -27.0*mm, -26.0*mm, -25.0*mm, -24.0*mm, -23.0*mm, -22.0*mm, -21.0*mm, -20.0*mm, -19.0*mm, -18.0*mm, 
        -17.0*mm, -16.0*mm, -15.0*mm, -14.0*mm, -13.0*mm, -12.0*mm, -11.0*mm, -10.0*mm, -9.0*mm, -8.0*mm, -7.0*mm, -6.0*mm, 
        -5.0*mm, -4.0*mm, -3.0*mm, -2.0*mm, -1.0*mm, 0.0*mm, 1.0*mm, 2.0*mm, 3.0*mm, 4.0*mm, 5.0*mm, 6.0*mm, 7.0*mm, 8.0*mm, 
        9.0*mm, 10.0*mm, 11.0*mm, 12.0*mm, 13.0*mm, 14.0*mm, 15.0*mm, 16.0*mm, 17.0*mm, 18.0*mm, 19.0*mm, 20.0*mm, 21.0*mm, 
        22.0*mm, 23.0*mm, 24.0*mm, 25.0*mm, 26.0*mm, 27.0*mm, 28.0*mm, 29.0*mm, 30.0*mm, 31.0*mm, 32.0*mm, 33.0*mm, 34.0*mm, 
        35.0*mm, 36.0*mm, 37.0*mm, 38.0*mm, 39.0*mm, 40.0*mm, 41.0*mm, 42.0*mm, 43.0*mm, 44.0*mm, 45.0*mm, 46.0*mm, 47.0*mm, 
        48.0*mm, 49.0*mm, 50.0*mm, 51.0*mm, 52.0*mm, 53.0*mm, 54.0*mm, 55.0*mm, 56.0*mm, 57.0*mm, 58.0*mm, 59.0*mm, 60.0*mm, 
        61.0*mm, 62.0*mm, 63.0*mm, 64.0*mm, 65.0*mm, 66.0*mm, 67.0*mm, 68.0*mm, 69.0*mm, 70.0*mm, 71.0*mm, 72.0*mm, 73.0*mm, 
        74.0*mm, 75.0*mm, 76.0*mm, 77.0*mm, 78.0*mm, 79.0*mm, 80.0*mm, 81.0*mm, 82.0*mm, 83.0*mm, 84.0*mm, 85.0*mm, 86.0*mm, 
        87.0*mm, 88.0*mm, 89.0*mm, 90.0*mm, 91.0*mm, 92.0*mm, 93.0*mm, 94.0*mm, 95.0*mm, 96.0*mm, 97.0*mm, 98.0*mm, 99.0*mm, 
        100.0*mm, 101.0*mm, 102.0*mm, 103.0*mm, 104.0*mm, 105.0*mm, 106.0*mm, 107.0*mm, 108.0*mm, 109.0*mm, 110.0*mm, 111.0*mm, 
        112.0*mm, 113.0*mm, 114.0*mm, 115.0*mm, 116.0*mm, 117.0*mm, 118.0*mm, 119.0*mm, 120.0*mm, 121.0*mm, 122.0*mm, 123.0*mm, 
        124.0*mm, 125.0*mm, 126.0*mm, 127.0*mm, 128.0*mm, 129.0*mm, 130.0*mm, 131.0*mm, 132.0*mm, 133.0*mm, 134.0*mm, 135.0*mm, 
        136.0*mm, 137.0*mm, 138.0*mm, 139.0*mm, 140.0*mm, 141.0*mm, 142.0*mm, 143.0*mm, 144.0*mm, 145.0*mm, 146.0*mm, 147.0*mm, 
        148.0*mm, 149.0*mm, 150.0*mm, 151.0*mm, 152.0*mm, 153.0*mm, 154.0*mm, 155.0*mm, 156.0*mm, 157.0*mm, 158.0*mm, 159.0*mm, 
        160.0*mm, 161.0*mm, 162.0*mm, 163.0*mm, 164.0*mm, 165.0*mm, 166.0*mm, 167.0*mm, 168.0*mm, 169.0*mm, 170.0*mm, 171.0*mm, 
        172.0*mm, 173.0*mm, 174.0*mm, 175.0*mm, 176.0*mm, 177.0*mm, 178.0*mm, 179.0*mm, 180.0*mm, 181.0*mm, 182.0*mm, 183.0*mm, 
        184.0*mm, 185.0*mm, 186.0*mm, 187.0*mm, 188.0*mm, 189.0*mm, 190.0*mm, 191.0*mm, 192.0*mm, 193.0*mm, 194.0*mm, 195.0*mm, 
        196.0*mm, 197.0*mm, 198.0*mm, 199.0*mm, 200.0*mm};
	G4double mlc_distal_actual[num_distal_pos] =
		{ -89.28967*mm,
        -88.830952*mm,
        -88.372866*mm,
        -87.915144*mm,
        -87.457423*mm,
        -86.999702*mm,
        -86.541981*mm,
        -86.08426*mm,
        -85.626538*mm,
        -85.168817*mm,
        -84.711096*mm,
        -84.253951*mm,
        -83.797275*mm,
        -83.340599*mm,
        -82.883923*mm,
        -82.427247*mm,
        -81.970571*mm,
        -81.513894*mm,
        -81.057218*mm,
        -80.600542*mm,
        -80.143866*mm,
        -79.688172*mm,
        -79.232592*mm,
        -78.777012*mm,
        -78.321432*mm,
        -77.865852*mm,
        -77.410273*mm,
        -76.954693*mm,
        -76.499113*mm,
        -76.043533*mm,
        -75.588068*mm,
        -75.133638*mm,
        -74.679208*mm,
        -74.224778*mm,
        -73.770348*mm,
        -73.315918*mm,
        -72.861488*mm,
        -72.407058*mm,
        -71.952628*mm,
        -71.498198*mm,
        -71.043768*mm,
        -70.590315*mm,
        -70.136958*mm,
        -69.683711*mm,
        -69.230488*mm,
        -68.777266*mm,
        -68.324043*mm,
        -67.87082*mm,
        -67.417597*mm,
        -66.964374*mm,
        -66.511151*mm,
        -66.058165*mm,
        -65.606131*mm,
        -65.154175*mm,
        -64.70222*mm,
        -64.250264*mm,
        -63.798309*mm,
        -63.346354*mm,
        -62.894398*mm,
        -62.442443*mm,
        -61.990487*mm,
        -61.538532*mm,
        -61.086577*mm,
        -60.635826*mm,
        -60.185202*mm,
        -59.734578*mm,
        -59.283954*mm,
        -58.833331*mm,
        -58.382707*mm,
        -57.932083*mm,
        -57.48146*mm,
        -57.030836*mm,
        -56.580212*mm,
        -56.129589*mm,
        -55.679742*mm,
        -55.230517*mm,
        -54.781293*mm,
        -54.332069*mm,
        -53.882845*mm,
        -53.43362*mm,
        -52.984396*mm,
        -52.535172*mm,
        -52.085947*mm,
        -51.636723*mm,
        -51.187499*mm,
        -50.738274*mm,
        -50.289771*mm,
        -49.842019*mm,
        -49.394266*mm,
        -48.946514*mm,
        -48.498762*mm,
        -48.05101*mm,
        -47.603257*mm,
        -47.155505*mm,
        -46.707753*mm,
        -46.26*mm,
        -45.812248*mm,
        -45.364496*mm,
        -44.916744*mm,
        -44.46956*mm,
        -44.023359*mm,
        -43.577158*mm,
        -43.130956*mm,
        -42.684755*mm,
        -42.238554*mm,
        -41.792353*mm,
        -41.346151*mm,
        -40.89995*mm,
        -40.453749*mm,
        -40.007736*mm,
        -39.561879*mm,
        -39.117084*mm,
        -38.67251*mm,
        -38.227937*mm,
        -37.783364*mm,
        -37.338791*mm,
        -36.894217*mm,
        -36.449644*mm,
        -36.005071*mm,
        -35.560497*mm,
        -35.115924*mm,
        -34.671704*mm,
        -34.228289*mm,
        -33.784873*mm,
        -33.341458*mm,
        -32.898042*mm,
        -32.454627*mm,
        -32.011211*mm,
        -31.567796*mm,
        -31.12438*mm,
        -30.680965*mm,
        -30.2381*mm,
        -29.795856*mm,
        -29.353611*mm,
        -28.911367*mm,
        -28.469123*mm,
        -28.026879*mm,
        -27.584635*mm,
        -27.14239*mm,
        -26.700146*mm,
        -26.257902*mm,
        -25.816654*mm,
        -25.375595*mm,
        -24.934535*mm,
        -24.493476*mm,
        -24.052416*mm,
        -23.611357*mm,
        -23.170298*mm,
        -22.729238*mm,
        -22.288179*mm,
        -21.847659*mm,
        -21.407799*mm,
        -20.967938*mm,
        -20.528077*mm,
        -20.088216*mm,
        -19.648355*mm,
        -19.208494*mm,
        -18.768633*mm,
        -18.328772*mm,
        -17.889236*mm,
        -17.450588*mm,
        -17.011941*mm,
        -16.573293*mm,
        -16.134645*mm,
        -15.695998*mm,
        -15.25735*mm,
        -14.818702*mm,
        -14.380054*mm,
        -13.941963*mm,
        -13.504548*mm,
        -13.067133*mm,
        -12.629718*mm,
        -12.192303*mm,
        -11.754888*mm,
        -11.317473*mm,
        -10.880058*mm,
        -10.442642*mm,
        -10.006115*mm,
        -9.569946*mm,
        -9.133777*mm,
        -8.697608*mm,
        -8.261439*mm,
        -7.82527*mm,
        -7.389101*mm,
        -6.952933*mm,
        -6.516764*mm,
        -6.081661*mm,
        -5.646745*mm,
        -5.211828*mm,
        -4.776911*mm,
        -4.341995*mm,
        -3.907079*mm,
        -3.472162*mm,
        -3.037246*mm,
        -2.602329*mm,
        -2.168327*mm,
        -1.734661*mm,
        -1.300996*mm,
        -0.867331*mm,
        -0.433665*mm,
        0.0*mm,
        0.43365*mm,
        0.8673*mm,
        1.30095*mm,
        1.733764*mm,
        2.166167*mm,
        2.59857*mm,
        3.030974*mm,
        3.463377*mm,
        3.895781*mm,
        4.328184*mm,
        4.760588*mm,
        5.192991*mm,
        5.624376*mm,
        6.055532*mm,
        6.486689*mm,
        6.917845*mm,
        7.349002*mm,
        7.780159*mm,
        8.211315*mm,
        8.642472*mm,
        9.073628*mm,
        9.503782*mm,
        9.933687*mm,
        10.363592*mm,
        10.793497*mm,
        11.223402*mm,
        11.653307*mm,
        12.083213*mm,
        12.513118*mm,
        12.943023*mm,
        13.371973*mm,
        13.800617*mm,
        14.229262*mm,
        14.657906*mm,
        15.08655*mm,
        15.515194*mm,
        15.943839*mm,
        16.372483*mm,
        16.801127*mm,
        17.228774*mm,
        17.656145*mm,
        18.083516*mm,
        18.510886*mm,
        18.938257*mm,
        19.365627*mm,
        19.792998*mm,
        20.220369*mm,
        20.647702*mm,
        21.073788*mm,
        21.499875*mm,
        21.925962*mm,
        22.352048*mm,
        22.778135*mm,
        23.204222*mm,
        23.630308*mm,
        24.056395*mm,
        24.481976*mm,
        24.906769*mm,
        25.331563*mm,
        25.756357*mm,
        26.18115*mm,
        26.605944*mm,
        27.030737*mm,
        27.455531*mm,
        27.880325*mm,
        28.303952*mm,
        28.727443*mm,
        29.150934*mm,
        29.574425*mm,
        29.997916*mm,
        30.421407*mm,
        30.844898*mm,
        31.26839*mm,
        31.691216*mm,
        32.113396*mm,
        32.535575*mm,
        32.957755*mm,
        33.379934*mm,
        33.802113*mm,
        34.224293*mm,
        34.646112*mm,
        35.067126*mm,
        35.48814*mm,
        35.909154*mm,
        36.330168*mm,
        36.751183*mm,
        37.172197*mm,
        37.593211*mm,
        38.014225*mm,
        38.434388*mm,
        38.854222*mm,
        39.274056*mm,
        39.693891*mm,
        40.113725*mm,
        40.533559*mm,
        40.953393*mm,
        41.373227*mm,
        41.792661*mm,
        42.2113*mm,
        42.62994*mm,
        43.04858*mm,
        43.467219*mm,
        43.885858*mm,
        44.304498*mm,
        44.723137*mm,
        45.141582*mm,
        45.559012*mm,
        45.976442*mm,
        46.393872*mm,
        46.811302*mm,
        47.228732*mm,
        47.646162*mm,
        48.063592*mm,
        48.480803*mm,
        48.897009*mm,
        49.313215*mm,
        49.729421*mm,
        50.145627*mm,
        50.561833*mm,
        50.978039*mm,
        51.394245*mm,
        51.809956*mm,
        52.224923*mm,
        52.63989*mm,
        53.054857*mm,
        53.469825*mm,
        53.884792*mm,
        54.299759*mm,
        54.714454*mm,
        55.128118*mm,
        55.541676*mm,
        55.955234*mm,
        56.368793*mm,
        56.782351*mm,
        57.195909*mm,
        57.609467*mm,
        58.022541*mm,
        58.434992*mm,
        58.847443*mm,
        59.259893*mm,
        59.672344*mm,
        60.084795*mm,
        60.49703*mm,
        60.90841*mm,
        61.319582*mm,
        61.730754*mm,
        62.141925*mm,
        62.553097*mm,
        62.964268*mm,
        63.37544*mm,
        63.785698*mm,
        64.195704*mm,
        64.60571*mm,
        65.015716*mm,
        65.425722*mm,
        65.835728*mm,
        66.245734*mm,
        66.654862*mm,
        67.063687*mm,
        67.472512*mm,
        67.881337*mm,
        68.290162*mm,
        68.698988*mm,
        69.107813*mm,
        69.516029*mm,
        69.923657*mm,
        70.331284*mm,
        70.738912*mm,
        71.14654*mm,
        71.554167*mm,
        71.961391*mm,
        72.368598*mm,
        72.775456*mm,
        73.181869*mm,
        73.588282*mm,
        73.994621*mm,
        74.400685*mm,
        74.806749*mm,
        75.212814*mm,
        75.61794*mm,
        76.022843*mm,
        76.427746*mm,
        76.832649*mm,
        77.237552*mm,
        77.642366*mm,
        78.046088*mm,
        78.449811*mm,
        78.853533*mm,
        79.257255*mm,
        79.660978*mm,
        80.0647*mm,
        80.46736*mm,
        80.869881*mm,
        81.272403*mm,
        81.674924*mm,
        82.077446*mm,
        82.479967*mm,
        82.882122*mm,
        83.283421*mm,
        83.68472*mm};


    for (size_t j = 0; j < 2; ++j) {
        for (size_t i = 0; i < num_leaves_in_bank; ++i) {
            double nominal_position = bank_positions[j]->at(i);
            if (j == 1) {nominal_position *= (-1);}
            for (size_t k = 0; k < num_distal_pos; ++k) {
                if (std::abs(nominal_position - mlc_distal_nominal[k]) <= 0.5) {
				    physical_position = (nominal_position - mlc_distal_nominal[k]) / (mlc_distal_nominal[k+1] - mlc_distal_nominal[k]) * (mlc_distal_actual[k+1] - mlc_distal_actual[k]) + mlc_distal_actual[k];
                    break;
                }
            }
            if (j == 1) {physical_position *= (-1);}
            auto volume = bank_volumes[j]->at(i);
            auto translation = volume->GetTranslation();
            translation.setX(physical_position);
            volume->SetTranslation(translation);
        }
    }
}
void buildMCBackscatterPlate(G4LogicalVolume* parent_logical, double parent_z_world)
{
    const double backscatter_plate_world_z = g_SAD - 145.8*mm;
    G4Tubs* backscatter_plate_tubs = new G4Tubs("backscatter_plate", 0., 100.*mm, 0.3*mm, 0., 360 * deg);
    G4LogicalVolume* backscatter_plate_logical = new G4LogicalVolume(backscatter_plate_tubs, G4Material::GetMaterial("SS304"), "backscatter_plate", 0, 0, 0);
    backscatter_plate_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.1, 0.7, 0.3, 0.5)));
    new G4PVPlacement(0, G4ThreeVector(0., 0., backscatter_plate_world_z - parent_z_world), backscatter_plate_logical, "backscatter_plate", parent_logical, false, 0);

}
void buildSecondaryCollimator(G4LogicalVolume* parent_logical, const double parent_world_z, double parent_gantry_z, const std::filesystem::path& stl_folder)
{   
     static const char* secondary_collimator[] = 
    {
        "p1040226067_ava.stl",
        "p1040226068_ava.stl", 
        "p1040226069_ava.stl"// 1010-1020 HRS
    };
    {
        std::filesystem::path filepath = stl_folder / secondary_collimator[0];
        auto leaf_mesh = CADMesh::TessellatedMesh::FromSTL(filepath.u8string());
        G4VSolid *leaf_solid = leaf_mesh->GetSolid();
        G4LogicalVolume *secondary_collimator_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("W95"), "secondary_collimator", 0, 0, 0);
        secondary_collimator_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.7, 0.2, 0.9, 0.5)));
        auto volume = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_gantry_z - parent_world_z), secondary_collimator_logical, "secondary_collimator", parent_logical, false, 0);
    }
    {
        std::filesystem::path filepath = stl_folder / secondary_collimator[1];
        auto leaf_mesh = CADMesh::TessellatedMesh::FromSTL(filepath.u8string());
        G4VSolid *leaf_solid = leaf_mesh->GetSolid();
        G4LogicalVolume *secondary_collimator_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("Lead97Antimony"), "secondary_collimator_shield", 0, 0, 0);
        secondary_collimator_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.7, 0.2, 0.9, 0.5)));
        auto volume = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_gantry_z - parent_world_z), secondary_collimator_logical, "secondary_collimator_shield", parent_logical, false, 0);
    }
    {
        std::filesystem::path filepath = stl_folder / secondary_collimator[2];
        auto leaf_mesh = CADMesh::TessellatedMesh::FromSTL(filepath.u8string());
        G4VSolid *leaf_solid = leaf_mesh->GetSolid();
        G4LogicalVolume *secondary_collimator_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("1010-1020HRS"), "secondary_collimator_shield", 0, 0, 0);
        secondary_collimator_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.7, 0.2, 0.9, 0.5)));
        auto volume = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_gantry_z - parent_world_z), secondary_collimator_logical, "secondary_collimator_shield", parent_logical, false, 0);
    
    }
}

void buildFoils(G4LogicalVolume* parent_logical, double parent_z_world, EnergyMode energyMode)
{
	// TODO be window?

  G4NistManager* NISTman = G4NistManager::Instance();
  NISTman->FindOrBuildMaterial("G4_Ta");
  G4Element* elAl = NISTman->FindOrBuildElement("Al");
  G4Element* elSi = NISTman->FindOrBuildElement("Si");
  G4Element* elMg = NISTman->FindOrBuildElement("Mg");
  G4Element* elCu = NISTman->FindOrBuildElement("Cu");
  G4Element* elCr = NISTman->FindOrBuildElement("Cr");
  G4Element* elZn = NISTman->FindOrBuildElement("Zn");
  G4Element* elPb = NISTman->FindOrBuildElement("Pb");
  G4Element* elFe = NISTman->FindOrBuildElement("Fe");

    //aluminum 6061
  double density = 2.70*g/cm3;
  G4Material* Aluminum6061 =
    new G4Material("Aluminum6061", density, 5);

  Aluminum6061->AddElement(elAl,98.01*perCent);
  Aluminum6061->AddElement(elSi, 0.6 *perCent);
  Aluminum6061->AddElement(elMg, 1.2 *perCent);
  Aluminum6061->AddElement(elCu, 0.15*perCent);
  Aluminum6061->AddElement(elCr, 0.04*perCent);

  //brass
  density = 8.53*g/cm3;
  G4Material* Brass = new G4Material("Brass", density, 4);
  Brass->AddElement(elZn,29.88*perCent);
  Brass->AddElement(elPb, 0.07*perCent);
  Brass->AddElement(elFe, 0.05*perCent);
  Brass->AddElement(elCu,70.  *perCent);


  //G4String beamName = "9E";

  // FOIL 1
  G4double foil1Radius = 6.*mm;
  G4double foil1PositionBottom = 62.7*mm;
  G4double foil1Thick = 0.23*mm;
  auto material = G4Material::GetMaterial("Brass");
  
  if (energyMode == EnergyMode::E09) {
    foil1Thick = 0.23*mm;
  }
  else if (energyMode == EnergyMode::E15) {
    foil1Thick = 0.136*mm;
    material = G4Material::GetMaterial("G4_Ta");
  }
  auto fEfoil1_tubs = new G4Tubs("efoil1", 0.0*mm, foil1Radius, foil1Thick/2.,
                      0.*deg, 360.*deg);
  G4LogicalVolume* fEfoil1_LV = new G4LogicalVolume(fEfoil1_tubs, material,
                                   "efoil1_LV", 0, 0, 0);

  auto foil1_pos = g_SAD - foil1PositionBottom + foil1Thick/2.;

  new G4PVPlacement(0,
        G4ThreeVector(0., 0., foil1_pos - parent_z_world),
        fEfoil1_LV, "electronFoil1", parent_logical, false, 0);


  // FOIL 2
  auto foil2material = G4Material::GetMaterial("Aluminum6061");
  double foil2Position = g_SAD - 135.526*mm;
  G4GenericPolycone* efoil2 = nullptr;
  if (energyMode == EnergyMode::E09) {
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
  else if (energyMode == EnergyMode::E15) {
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

  auto fEfoil2_LV = new G4LogicalVolume(efoil2,
    foil2material, "efoil2_LV", 0, 0, 0);

  auto fEfoil2 = new G4PVPlacement(0,
        G4ThreeVector(0., 0., foil2Position - parent_z_world),
        fEfoil2_LV, "electronFoil2", parent_logical, false, 0);

}

void buildApplicator(G4LogicalVolume* parent_logical, double parent_z_world)
{
	G4NistManager* NISTman = G4NistManager::Instance();
	G4Element* elZn = NISTman->FindOrBuildElement("Zn");
	G4Element* elAl = NISTman->FindOrBuildElement("Al");
	G4Element* elCu = NISTman->FindOrBuildElement("Cu");
	G4Element* elFe = NISTman->FindOrBuildElement("Fe");
	G4Element* elMn = NISTman->FindOrBuildElement("Mn");
	G4Element* elPb = NISTman->FindOrBuildElement("Pb");
	G4Element* elCd = NISTman->FindOrBuildElement("Cd");
	G4Element* elSn = NISTman->FindOrBuildElement("Sn");
  // Zinc ZA8
  double density = 6.3*g/cm3;
  G4Material* ZincZA8 = new G4Material("ZincZA8", density, 8);
  ZincZA8->AddElement(elZn, 89.78*perCent);
  ZincZA8->AddElement(elAl,  8.8*perCent);
  ZincZA8->AddElement(elCu,  1.3*perCent);
  ZincZA8->AddElement(elFe,  0.075*perCent);
  ZincZA8->AddElement(elMn,  0.03*perCent);
  ZincZA8->AddElement(elPb,  0.006*perCent);
  ZincZA8->AddElement(elCd,  0.006*perCent);
  ZincZA8->AddElement(elSn,  0.003*perCent);

	G4Element* elSi = NISTman->FindOrBuildElement("Si");
	G4Element* elMg = NISTman->FindOrBuildElement("Mg");
	G4Element* elCr = NISTman->FindOrBuildElement("Cr");

	//aluminum 6061
	//density = 2.70*g/cm3;
	//G4Material* Aluminum6061 =
	//	new G4Material("Aluminum6061", density, 5);

	//Aluminum6061->AddElement(elAl,98.01*perCent);
	//Aluminum6061->AddElement(elSi, 0.6 *perCent);
	//Aluminum6061->AddElement(elMg, 1.2 *perCent);
	//Aluminum6061->AddElement(elCu, 0.15*perCent);
	//Aluminum6061->AddElement(elCr, 0.04*perCent);

	G4Element* elBi = NISTman->FindOrBuildElement("Bi");
  
	// cerrotru   // from applicator drawing
  density = 6.3*g/cm3;  //  csalloys.com  Tru 281
  G4Material* cerrotru = new G4Material("cerrotru", density, 2);
  cerrotru->AddElement(elBi,58.*perCent);
  cerrotru->AddElement(elSn,42.*perCent);


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
  G4double sc3_width = 42.*mm;  // x/y dimension of scraper
  G4double r31 = -1., r32 = -1.;

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
  r22      = rtop_2 + sc2_width[1];
  r23      = rtop_2 + sc2_width[2];
  r24      = rtop_2 + sc2_width[3];
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

	//G4double collPos = 430.*mm;

  // scraper 1
  const G4int numRZ_1 = 8;
  const G4double scr_r_1[numRZ_1] =
    {rtop_1, r11, r12,   r13, r13,
     r11,    r11, rbot_1};
  const G4double scr_z_1[numRZ_1] =
    {thick_1, thick_1, thick_1/2., thick_1/2., t11,
     t11, 0.,  0.};

  G4Polyhedra* scraper1_poly =
    new G4Polyhedra("scraper1", pi/4., 9.*pi/4., 4., numRZ_1, scr_r_1, scr_z_1);

  G4LogicalVolume* scraper1_LV =
      new G4LogicalVolume(scraper1_poly,
                          G4Material::GetMaterial("ZincZA8"),
                          "scraper1", 0, 0, 0);

 // 0 is 682 mm from isocenter minus 245 equals 437
 // bottom of scraper is 346 mm from iso
 //new G4PVPlacement(0, G4ThreeVector(0., 0., -336.*mm),
 new G4PVPlacement(0, G4ThreeVector(0., 0., -101.*mm),
      scraper1_LV, "scraper1", parent_logical, false, 0);

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

  new G4PVPlacement(0, G4ThreeVector(0., 0., + t11/2. - 101.*mm),
                scraper1Support_LV, "scraper1_support", parent_logical, false, 0);

  //scraper 2
  const G4int numRZ_2 = 10;
  const G4double scr_r_2[numRZ_2] =
    //{rtop_2, rtop_2 + 20.*mm, rtop_2 + 50.*mm, rtop_2 + 20.*mm, rbot_2};
    {rtop_2, r21, r22, r23, r24, r24, r23, r22, r21, rbot_2};
  const G4double scr_z_2[numRZ_2] =
    {thick_2, thick_2, t21, t22, t22, t23, t23, t24, 0., 0.};
    //{thick_2, thick_2, thick_2/2., 0., 0.};

  G4Polyhedra* scraper2poly =
    new G4Polyhedra("scraper2", pi/4., 9.*pi/4., 4., numRZ_2, scr_r_2, scr_z_2);

  G4LogicalVolume* scraper2_LV =
      new G4LogicalVolume(scraper2poly,
                          G4Material::GetMaterial("ZincZA8"),
                          "scraper2", 0, 0, 0);

  // bottom is 215 mm from isocentre
  new G4PVPlacement(0, G4ThreeVector(0., 0., -222.*mm),
      scraper2_LV, "scraper2", parent_logical, false, 0);

  ////scraper 3
  const G4int numRZ_3 = 4;
  const G4double scr_r_3[numRZ_3] = {rtop_3, r31, r32, rbot_3};
  const G4double scr_z_3[numRZ_3] = {thick_3, thick_3, 0., 0.};

  G4Polyhedra* scraper3poly =
    new G4Polyhedra("scraper3", pi/4., 9.*pi/4., 4., numRZ_3, scr_r_3, scr_z_3);

  G4LogicalVolume* scraper3_LV =
      new G4LogicalVolume(scraper3poly,
                          G4Material::GetMaterial("ZincZA8"),
                          "scraper3", 0, 0, 0);

  // bottom is 50 mm from isocenter
  new G4PVPlacement(0, G4ThreeVector(0., 0., -387.*mm),
      scraper3_LV, "scraper3", parent_logical, false, 0);

  G4VisAttributes* scraper1_VisAtt =
    new G4VisAttributes(G4Colour(0.5, 0.5, 1.0, 0.9));
  G4VisAttributes* scraper1Support_VisAtt =
    new G4VisAttributes(G4Colour(0.5, 1.0, 0.5, 0.9));
  G4VisAttributes* scraper2_VisAtt =
    new G4VisAttributes(G4Colour(1.0, 0.5, 0.2, 0.9));
  G4VisAttributes* scraper3_VisAtt =
    new G4VisAttributes(G4Colour(1.0, 1.0, 0.2, 0.9));
  scraper1_LV->SetVisAttributes(scraper1_VisAtt);
  scraper1Support_LV->SetVisAttributes(scraper1Support_VisAtt);
  scraper2_LV->SetVisAttributes(scraper2_VisAtt);
  scraper3_LV->SetVisAttributes(scraper3_VisAtt);



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
        {
            //wafer alloy
            const double density = 11.0 * g / cm3;
            G4Material* CuAuAlloy = new G4Material("CuAuAlloy", density, 2);
            CuAuAlloy->AddElement(elAu, 35.*perCent);
            CuAuAlloy->AddElement(elCu, 65.*perCent);
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

		// Target Button
		G4double nicoroHThick = 0.0508 * mm; // half the nicoro brazing
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
		G4Tubs* xray_window_tubs = new G4Tubs("xray_window", 0.0 * mm, 6.4135 * mm, 0.5 * mm / 2., 0. * deg, 360. * deg);
		auto xray_window_logical =
			new G4LogicalVolume(xray_window_tubs, G4Material::GetMaterial("SS304"), "xray_window", 0, 0, 0);
		xray_window_logical->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 0.0, 1.0)));

		new G4PVPlacement(0,
			G4ThreeVector(0. * m, 0. * m, target_z_world - 8.326 * mm - 0.5 * mm / 2. - parent_z_world),
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
    else if (target_type == Target::MEDIUM_ENERGY)
    {
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

		// nicoro 
		G4Tubs* target_nicoro =
			new G4Tubs("nicoro", 0.0 * mm, 3.0734 * mm, nicoroHThick, 0. * deg, 360. * deg);

		auto target_nicoro_logical =
			new G4LogicalVolume(target_nicoro, G4Material::GetMaterial("Nicoro"),
				"target_nicoro", 0, 0, 0);
		target_nicoro_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.0, 0.5, 0.5)));


		new G4PVPlacement(0,
			G4ThreeVector(0., 0., target_z_world - Target_z - 0.635 * mm + nicoroHThick - parent_z_world),
			target_nicoro_logical, "target_nicoro", parent_logical, false, 0);


        const G4int tbt_10x_points = 9;
        G4double r_pc[tbt_10x_points] = {0.0*mm, 3.1*mm, 3.1   *mm, 8.0   *mm, 8.0*mm, 5.5*mm, 4.0677*mm, 3.7875*mm, 0.0*mm};
        G4double z_pc[tbt_10x_points] = {.924*mm,   .924*mm,  1.124*mm,  1.124*mm, -8.876*mm, -8.876*mm, -5.369*mm, -4.776*mm, -4.776*mm};

        G4GenericPolycone* target_block_top =
			new G4GenericPolycone("target_block_top", 0. * deg, 360. * deg, tbt_10x_points, r_pc, z_pc);

		auto target_block_top_logical =
			new G4LogicalVolume(target_block_top, G4Material::GetMaterial("COPPER_GLIDCOP"),
				"target_block_top", 0, 0, 0);
		target_block_top_logical->SetVisAttributes(target_block_vis_attributes);

		new G4PVPlacement(0, G4ThreeVector(0., 0., target_z_world - parent_z_world),
				target_block_top_logical, "target_block_top", parent_logical, false, 0);

        target_block_top_logical->SetRegion(target_region);
		target_region->AddRootLogicalVolume(target_block_top_logical);

        // X-ray window
		G4Tubs* xray_window_tubs = new G4Tubs("xray_window", 0.0 * mm, 6.4135 * mm, 0.5 * mm / 2., 0. * deg, 360. * deg);
		auto xray_window_logical =
			new G4LogicalVolume(xray_window_tubs, G4Material::GetMaterial("SS304"), "xray_window", 0, 0, 0);
		xray_window_logical->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 0.0, 1.0)));
        new G4PVPlacement(0,
			G4ThreeVector(0. * m, 0. * m, target_z_world - 8.326 * mm - 0.5 * mm / 2. - parent_z_world),
			xray_window_logical, "xray_window", parent_logical, false, 0);

    }
    else if (target_type == Target::HIGH_ENERGY){
        
        const G4int he_pts = 5;
        G4double r_tb[he_pts] = {0.0*mm, 3.0734*mm,  3.0734*mm,  2.8194*mm,  0.0  *mm};
        G4double z_tb[he_pts] = {0.0*mm, 0.0   *mm, -0.381 *mm, -0.635 *mm, -0.635*mm};

        G4double Target_z = -0.859*mm;

        G4Polycone* target_button = new G4Polycone("TargetButton",0.*deg, 360.*deg, he_pts, r_tb, z_tb);
        auto target_logical = new G4LogicalVolume(target_button, G4Material::GetMaterial("G4_W"), "target_button", 0, 0, 0);
		target_logical->SetVisAttributes(target_vis_attributes);

		new G4PVPlacement(0, G4ThreeVector(0., 0., target_z_world - Target_z - parent_z_world),
			target_logical, "target_button", parent_logical,
			false, 0);

        const G4int tbt_he_pts = 6;
        G4double r_pc[tbt_he_pts] = { 0.0  *mm, 3.1  *mm, 3.1  *mm, 7.0  *mm,  7.0  *mm,  0.0  *mm};
        G4double z_pc[tbt_he_pts] = { 0.224*mm, 0.224*mm, 1.124*mm, 1.124*mm, -5.876*mm, -5.876*mm}; 
        G4GenericPolycone* target_block_top =
			new G4GenericPolycone("target_block_top", 0. * deg, 360. * deg, tbt_he_pts, r_pc, z_pc);

		auto target_block_top_logical =
			new G4LogicalVolume(target_block_top, G4Material::GetMaterial("COPPER_GLIDCOP"),
				"target_block_top", 0, 0, 0);
		target_block_top_logical->SetVisAttributes(target_block_vis_attributes);

		new G4PVPlacement(0, G4ThreeVector(0., 0., target_z_world - parent_z_world),
				target_block_top_logical, "target_block_top", parent_logical, false, 0);

        
        //------------------------------------------------------------------------
        // braze wafer between target block top and target block bottom
        //------------------------------------------------------------------------
        // 35% Au & 65% Cu, 0.076 mm thick
        G4double waferThick = 0.076*mm;
        G4Tubs* target_top_wafer = new G4Tubs("target_top_wafer", 0.*mm, 7.*mm, waferThick/2., 0.*deg, 360.*deg);
        G4double waferPos = 5.876*mm; //using the target block shift
        auto target_top_wafer_logical = new G4LogicalVolume(target_top_wafer, G4Material::GetMaterial("CuAuAlloy"), "target_top_wafer", 0, 0, 0);
        target_top_wafer_logical->SetVisAttributes(target_block_vis_attributes);
        new G4PVPlacement(0, G4ThreeVector(0., 0., target_z_world - waferPos - waferThick/2. - parent_z_world),
            target_top_wafer_logical, "target_top_wafer", parent_logical, false, 0);

        //------------------------------------------------------------------------
        //target block bottom
        //------------------------------------------------------------------------

        G4double bottom_thick = 2.924*mm;
        G4Tubs* target_bottom = new G4Tubs("target_bottom", 0.*mm, 7.0*mm, bottom_thick/2, 0.*deg, 360.*deg);
        G4double target_bottom_pos = waferPos + waferThick;
        auto target_bottom_logical = new G4LogicalVolume(target_bottom, G4Material::GetMaterial("COPPER_GLIDCOP"), "target_bottom", 0, 0, 0);
        new G4PVPlacement(0, G4ThreeVector(0., 0., target_z_world -target_bottom_pos - bottom_thick/2. - parent_z_world),
            target_bottom_logical, "target_bottom",
            parent_logical, false, 0);
        target_bottom_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.0, 0.5, 0.5)));
       
        target_logical->SetRegion(target_region);
		target_block_top_logical->SetRegion(target_region);
        target_bottom_logical->SetRegion(target_region);
        target_top_wafer_logical->SetRegion(target_region);
		target_region->AddRootLogicalVolume(target_logical);
		target_region->AddRootLogicalVolume(target_block_top_logical);
		target_region->AddRootLogicalVolume(target_bottom_logical);
        target_region->AddRootLogicalVolume(target_top_wafer_logical);

    }
	else {
		throw std::runtime_error("Target type not supported in buildTarget");
	}
}


void buildPrimaryCollimator(G4LogicalVolume* parent_logical, const double parent_world_z, const std::filesystem::path& stl_folder) 
{   
    static const char* primary_collimator[] = 
    {
        "p1051568069_carr-ava.stl", //AL 6061-T651
        "p1051927001_carr-ava.stl", // ASTM A36 STEEL
        "p1055251001_carr-ava.stl", // TUNGSTEN AMS-T-21014, CLASS 3
        "p1058166069_carr-ava.stl", //LEAD 3% +/- 1% ANTIMONY
        "p1058169069_carr-ava.stl" //LEAD 3% +/- 1% ANTIMONY
    };
    {
        std::filesystem::path filepath = stl_folder / primary_collimator[0];
        auto leaf_mesh = CADMesh::TessellatedMesh::FromSTL(filepath.u8string());
        G4VSolid *leaf_solid = leaf_mesh->GetSolid();
        G4LogicalVolume *primary_collimator_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("Al6061-T651"), "primary_collimator", 0, 0, 0);
        primary_collimator_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.8, 0.1, 0.8, 0.5)));
        auto volume = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_world_z), primary_collimator_logical, "primary_collimator", parent_logical, false, 0);
    
    }
    {
        std::filesystem::path filepath = stl_folder / primary_collimator[1];
        auto leaf_mesh = CADMesh::TessellatedMesh::FromSTL(filepath.u8string());
        G4VSolid *leaf_solid = leaf_mesh->GetSolid();
        G4LogicalVolume *primary_collimator_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("SS_A36"), "primary_collimator", 0, 0, 0);
        primary_collimator_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.8, 0.1, 0.8, 0.5)));
        auto volume = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_world_z), primary_collimator_logical, "primary_collimator", parent_logical, false, 0);
    
    }
    {
        std::filesystem::path filepath = stl_folder / primary_collimator[2];
        auto leaf_mesh = CADMesh::TessellatedMesh::FromSTL(filepath.u8string());
        G4VSolid *leaf_solid = leaf_mesh->GetSolid();
        G4LogicalVolume *primary_collimator_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("W95"), "primary_collimator", 0, 0, 0);
        primary_collimator_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.8, 0.1, 0.8, 0.5)));
        auto volume = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_world_z), primary_collimator_logical, "primary_collimator", parent_logical, false, 0);
    
    }
    {
        std::filesystem::path filepath = stl_folder / primary_collimator[3];
        auto leaf_mesh = CADMesh::TessellatedMesh::FromSTL(filepath.u8string());
        G4VSolid *leaf_solid = leaf_mesh->GetSolid();
        G4LogicalVolume *primary_collimator_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("Lead97Antimony"), "primary_collimator_shield", 0, 0, 0);
        primary_collimator_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.8, 0.1, 0.8, 0.5)));
        auto volume = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_world_z), primary_collimator_logical, "primary_collimator_shield", parent_logical, false, 0);
    
    }
    {
        std::filesystem::path filepath = stl_folder / primary_collimator[4];
        auto leaf_mesh = CADMesh::TessellatedMesh::FromSTL(filepath.u8string());
        G4VSolid *leaf_solid = leaf_mesh->GetSolid();
        G4LogicalVolume *primary_collimator_logical = new G4LogicalVolume(leaf_solid, G4Material::GetMaterial("Lead97Antimony"), "primary_collimator_shield", 0, 0, 0);
        primary_collimator_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.8, 0.1, 0.8, 0.5)));
        auto volume = new G4PVPlacement(0, G4ThreeVector(0., 0., g_SAD - parent_world_z), primary_collimator_logical, "primary_collimator_shield", parent_logical, false, 0);
    
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
    else if (flattening_filter == FlatteningFilter::ENERGY_10X) {
        G4RotationMatrix* FF_rot = new G4RotationMatrix();
		FF_rot->rotateX(180.0*deg);
        const G4int cone1size = 34;
        G4double r_FF[cone1size] = {
        0.0 *mm,  0.635*mm,  1.6 *mm,  2.54 *mm,
        3.81*mm,  5.08 *mm,  7.62*mm, 10.16 *mm,
        12.7 *mm, 15.24 *mm, 17.78*mm, 20.32 *mm,
        22.86*mm, 25.4  *mm, 27.94*mm, 30.48 *mm,
        38.1 *mm, 38.1  *mm, 31.00*mm, 29.464*mm,
        27.94*mm, 25.4  *mm, 22.86*mm, 20.32 *mm,
        17.78*mm, 15.24 *mm, 12.7 *mm, 10.16 *mm,
        7.62*mm,  5.08 *mm,  3.81*mm,  2.54 *mm,
        1.27*mm,  0.0  *mm};

        G4double z_FF[cone1size] = {
        -32.131*mm, -32.131*mm,-30.912*mm, -29.515*mm,
        -27.381*mm, -25.349*mm,-21.742*mm, -18.618*mm,
        -15.748*mm, -12.954*mm,-10.084*mm,  -7.341*mm,
        -5.004*mm,  -3.226*mm, -1.575*mm,   0.0  *mm,
            0.0  *mm,   2.921*mm,  2.921*mm,   3.81 *mm,
            3.81 *mm,   3.81 *mm,  3.81 *mm,   3.988*mm,
            4.191*mm,   4.699*mm,  5.588*mm,   6.68 *mm,
            7.823*mm,   9.093*mm,  9.652*mm,  10.16 *mm,
        10.541*mm,  10.77 *mm};

    G4GenericPolycone* ff_10X = new G4GenericPolycone("ff_10X", 0.*deg, 360.*deg, cone1size, r_FF, z_FF);
	G4LogicalVolume* FF_10X_logical = new G4LogicalVolume(ff_10X, G4Material::GetMaterial("G4_Cu"), "ff_10X",0,0,0);
	FF_10X_logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.1, 1.0, 0.1)));
	new G4PVPlacement(FF_rot, G4ThreeVector(ff_offset.x(), ff_offset.y(), flattening_filter_z_world + ff_offset.z() - parent_z_world),
							FF_10X_logical, "ff_10X", parent_logical, false, 0);
    }
	else {
		throw std::runtime_error("Given flattening filter not supported in buildFlatteningFilter");
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

VolumeNameAndTraversalFlag AvalonTable[] = {
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
    {"BeWindow6X_LV", TraversedGeometry::NONE},
    {"xray_window", TraversedGeometry::NONE},
    {"target_button", TraversedGeometry::TARGET},
    {"target_block_top", TraversedGeometry::TARGET},
    {"target_nicoro", TraversedGeometry::TARGET},
    {"target_bottom", TraversedGeometry::TARGET},
    {"target_top_wafer", TraversedGeometry::TARGET},
    {"vacuum_chamber", TraversedGeometry::TARGET},
    {"primary_collimator", TraversedGeometry::PRIMARY_COLLIMATOR},
    {"primary_collimator_shield", TraversedGeometry::PRIMARY_COLLIMATOR},
    {"mlc_proximal_leaf", TraversedGeometry::MLC},
    {"mlc_distal_leaf", TraversedGeometry::MLC},
    {"mlc_proximal_outer_leaf", TraversedGeometry::MLC},
    {"mlc_distal_outer_leaf", TraversedGeometry::MLC},
    {"secondary_collimator", TraversedGeometry::SHIELD_COLLIMATOR},
    {"secondary_collimator_shield", TraversedGeometry::SHIELD_COLLIMATOR},
    {"ff_carousel", TraversedGeometry::FLATTENING_FILTER},
    {"open_port", TraversedGeometry::FLATTENING_FILTER},
    {"ff_6X", TraversedGeometry::FLATTENING_FILTER},
    {"ff_10X", TraversedGeometry::FLATTENING_FILTER},
    {"efoil1_LV", TraversedGeometry::EFOIL1},
    {"efoil2_LV", TraversedGeometry::EFOIL2},
    {"scraper1", TraversedGeometry::SCRAPER1},
    {"scraper1_support", TraversedGeometry::SCRAPER1_SUPPORT},
    {"scraper2", TraversedGeometry::SCRAPER2},
    {"scraper3", TraversedGeometry::SCRAPER3},


};


TreatmentHeadDetector::TreatmentHeadDetector(EnergyMode energy_mode, const std::string& sd_monitor_chamber_name, 
    const fs::path& gdml_path, const fs::path& stl_path) 
    :   m_energy_mode(energy_mode),m_monitor_chamber(nullptr), m_gdml_path(gdml_path), 
        m_stl_path(stl_path), m_sd_monitor_chamber_name(sd_monitor_chamber_name)
{}

G4VPhysicalVolume* TreatmentHeadDetector::Construct() {
	//We actually want to store physical volumes for rotations
	buildMaterials();
	G4NistManager* NISTman = G4NistManager::Instance();
	auto world_material = NISTman->FindOrBuildMaterial("G4_AIR");

  const G4double c1 = 245.*mm;//345.*mm;
  const double collimator_position_z = -8.*mm - c1;
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
		//G4Tubs* collimator_tubs = new G4Tubs("collimator", 0. * mm, 300. * mm, 170. * mm + 462.*mm, 0. * deg, 360. * deg);
		G4Tubs* collimator_tubs = new G4Tubs("collimator", 0. * mm, 300. * mm, 170. * mm + 462.*mm - c1, 0. * deg, 360. * deg);
		collimator_logical = new G4LogicalVolume(collimator_tubs, world_material, "collimator", 0, 0, 0);
		m_collimator = new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, collimator_position_z), 
			collimator_logical, "collimator", gantry_logical, false, 0);
	}
	const bool build_photon = false;
	const bool build_electron = true;

    fs::path cad_path = m_stl_path / "AVALON";
	if (build_photon) {
        Target target;
        if((m_energy_mode == EnergyMode::X06_FFF) | (m_energy_mode == EnergyMode::X06))
			target = Target::LOW_ENERGY;
		else if ((m_energy_mode == EnergyMode::X10) | (m_energy_mode == EnergyMode::X08) | (m_energy_mode == EnergyMode::X08_FFF))
			target = Target::MEDIUM_ENERGY;
        else if ((m_energy_mode == EnergyMode::X10_FFF) |(m_energy_mode == EnergyMode::X15))
			target = Target::HIGH_ENERGY;
		else
			throw std::runtime_error("Unknown energy mode for Avalon target");
		buildTarget(gantry_logical, gantry_position_z, target);

        FlatteningFilter flattening_filter;
		if ((m_energy_mode == EnergyMode::X10_FFF) || (m_energy_mode == EnergyMode::X06_FFF))
			flattening_filter = FlatteningFilter::OPEN;
        else if(m_energy_mode == EnergyMode::X06)
			flattening_filter = FlatteningFilter::ENERGY_6X;
		else if(m_energy_mode == EnergyMode::X10)
			flattening_filter = FlatteningFilter::ENERGY_10X;
		else
			throw std::runtime_error("Unknown energy mode for Truebeam flattening filter");
		buildFlatteningFilter(gantry_logical, gantry_position_z, flattening_filter);

		buildPrimaryCollimator(gantry_logical, gantry_position_z, cad_path);
        buildMCBackscatterPlate(gantry_logical, gantry_position_z);
        buildSecondaryCollimator(collimator_logical, collimator_position_z, gantry_position_z, cad_path);
	}
  else if (build_electron) {
    buildSecondaryCollimator(collimator_logical, collimator_position_z, gantry_position_z, cad_path);
    buildFoils(gantry_logical, gantry_position_z, m_energy_mode);
		buildApplicator(collimator_logical, collimator_position_z);
  }
	
	m_monitor_chamber = buildIonisationChamber(gantry_logical, gantry_position_z);
    m_mlc_proximal_volumes = buildProximalLeafBank(collimator_logical, collimator_position_z, gantry_position_z, cad_path);
    m_mlc_distal_volumes = buildDistalLeafBank(collimator_logical, collimator_position_z, gantry_position_z, cad_path);

	m_id_to_traversed = generateVolumeToTraversed(AvalonTable, sizeof(AvalonTable) / sizeof(AvalonTable[0]), world_logical);
	return world_physical;
}

void TreatmentHeadDetector::setState(const Plan& plan, size_t pt_idx)
{
	G4GeometryManager::GetInstance()->OpenGeometry(m_gantry);
	const double gantry_position_z = 690.0 * mm;
    const auto& state = plan.getDualLayerMLC(TreatmentHeadType::AVALON_ELECTRON, pt_idx);
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
