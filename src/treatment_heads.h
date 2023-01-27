#pragma once
#include "linac.h"
#include <filesystem>
#include <unordered_map>
namespace fs = std::filesystem;

#include "G4VUserDetectorConstruction.hh"

class Plan;
class G4LogicalVolume;
class LinacDetector : public G4VUserDetectorConstruction {
public:
  virtual std::optional<std::string> getMonitorChamberName() const = 0;
  virtual void setState(const Plan& plan, size_t pt_idx) = 0;
  virtual const std::unordered_map<int, TraversedGeometry>& getTraversedGeometryMapping() const = 0;
  virtual ~LinacDetector() = default;
};

std::unordered_map<int, TraversedGeometry> generateVolumeToTraversed(
  const VolumeNameAndTraversalFlag* _array, 
  const size_t _array_size, 
  G4LogicalVolume* _root_volume);

namespace truebeam {

  struct JawVolumes {
    G4VPhysicalVolume* y1;
    G4VPhysicalVolume* y2;
    G4VPhysicalVolume* x1;
    G4VPhysicalVolume* x2;
    double parent_z_world;
  };

  struct MLCLeafVolumes {
    std::vector<G4VPhysicalVolume*> bank_X1;
    std::vector<G4VPhysicalVolume*> bank_X2;
  };

  class TreatmentHeadDetector : public LinacDetector {
  public:
    TreatmentHeadDetector(MLCType mlc_type, EnergyMode energy_mode,  
      const std::string& sd_monitor_chamber_name, 
      const fs::path& gdml_path, const fs::path& stl_path);
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    void setState(const Plan& plan, size_t pt_idx) override;
    const std::unordered_map<int, TraversedGeometry>& getTraversedGeometryMapping() const override {
      return m_id_to_traversed;
    }
    std::optional<std::string> getMonitorChamberName() const override {
      return m_sd_monitor_chamber_name;
    }
    //virtual const G4RotationMatrix* getGantryRotation() { return m_gantry->GetRotation(); }
    //virtual const G4RotationMatrix* getCollimatorRotation() {  return m_collimator->GetRotation(); };

    MLCType m_mlc_type;
    EnergyMode m_energy_mode;
    fs::path m_gdml_path;
    fs::path m_stl_path;
    std::string m_sd_monitor_chamber_name;
    std::unordered_map<int, TraversedGeometry> m_id_to_traversed;
    G4VPhysicalVolume* m_gantry;
    G4VPhysicalVolume* m_collimator;
    G4LogicalVolume* m_monitor_chamber;
    JawVolumes m_jaw_volumes;
    MLCLeafVolumes m_mlc_leaf_volumes;
  };

}

namespace halcyon {
  struct MLCLeafVolumesProximal {
    std::vector<G4VPhysicalVolume*> bank_X1;
    std::vector<G4VPhysicalVolume*> bank_X2;
  };
  struct MLCLeafVolumesDistal {
    std::vector<G4VPhysicalVolume*> bank_X1;
    std::vector<G4VPhysicalVolume*> bank_X2;
  };
  class TreatmentHeadDetector : public LinacDetector {
  public:
    TreatmentHeadDetector(const std::string& sd_monitor_chamber_name, 
      const fs::path& gdml_path, const fs::path& stl_path);
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    void setState(const Plan& plan, size_t pt_idx) override;
    const std::unordered_map<int, TraversedGeometry>& getTraversedGeometryMapping() const override {
      return m_id_to_traversed;
    }
    std::optional<std::string> getMonitorChamberName() const override {
      return m_sd_monitor_chamber_name;
    }
    //virtual const G4RotationMatrix* getGantryRotation() { return m_gantry->GetRotation(); }
    //virtual const G4RotationMatrix* getCollimatorRotation() {  return m_collimator->GetRotation(); };

    fs::path m_gdml_path;
    fs::path m_stl_path;
    std::string m_sd_monitor_chamber_name;
    std::unordered_map<int, TraversedGeometry> m_id_to_traversed;
    G4VPhysicalVolume* m_gantry;
    G4VPhysicalVolume* m_collimator;
    G4LogicalVolume* m_monitor_chamber;
    MLCLeafVolumesProximal m_mlc_proximal_volumes;
    MLCLeafVolumesDistal m_mlc_distal_volumes;
  };
}
namespace avalon {
  struct MLCLeafVolumesProximal {
    std::vector<G4VPhysicalVolume*> bank_X1;
    std::vector<G4VPhysicalVolume*> bank_X2;
  };
  struct MLCLeafVolumesDistal {
    std::vector<G4VPhysicalVolume*> bank_X1;
    std::vector<G4VPhysicalVolume*> bank_X2;
  };
  class TreatmentHeadDetector : public LinacDetector {
  public:
    TreatmentHeadDetector(EnergyMode energy_mode,const std::string& sd_monitor_chamber_name, 
      const fs::path& gdml_path, const fs::path& stl_path);
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    void setState(const Plan& plan, size_t pt_idx) override;
    const std::unordered_map<int, TraversedGeometry>& getTraversedGeometryMapping() const override {
      return m_id_to_traversed;
    }
    std::optional<std::string> getMonitorChamberName() const override {
      return m_sd_monitor_chamber_name;
    }
    //virtual const G4RotationMatrix* getGantryRotation() { return m_gantry->GetRotation(); }
    //virtual const G4RotationMatrix* getCollimatorRotation() {  return m_collimator->GetRotation(); };
    
    EnergyMode m_energy_mode;
    fs::path m_gdml_path;
    fs::path m_stl_path;
    std::string m_sd_monitor_chamber_name;
    std::unordered_map<int, TraversedGeometry> m_id_to_traversed;
    G4VPhysicalVolume* m_gantry;
    G4VPhysicalVolume* m_collimator;
    G4LogicalVolume* m_monitor_chamber;
    MLCLeafVolumesProximal m_mlc_proximal_volumes;
    MLCLeafVolumesDistal m_mlc_distal_volumes;
  };
}

namespace nohead {
  class TreatmentHeadDetector : public LinacDetector {
  public:
    TreatmentHeadDetector(const std::string& background_material)
      : m_background_material(background_material) {}
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;
    void setState(const Plan&, size_t) override {}
    const std::unordered_map<int, TraversedGeometry>& getTraversedGeometryMapping() const override
    { return m_id_to_traversed; }
    std::optional<std::string> getMonitorChamberName() const override {
      return std::nullopt;
    }

    std::string m_background_material;
    std::unordered_map<int, TraversedGeometry> m_id_to_traversed;
  };
}