#include "linac.h"
#include <stdlib.h>
#include <filesystem>
namespace fs = std::filesystem;
#include <iostream>
#include <fstream>
#include <chrono>
#include <variant>
#include <thread>
#include <random>
#include "logger.h"
#include "source.h"
#include "plan.h"
#include "treatment_heads.h"

#include <fmt/core.h>
#include "external/popl.hpp"
#define JSON_DIAGNOSTICS 1
#include "external/json.hpp"
using json = nlohmann::json;
#include "parse.h"

#include "G4SystemOfUnits.hh"

void from_json(const json& j, TreatmentHeadType& type) {
  const size_t num = static_cast<size_t>(TreatmentHeadType::NUM_TREATMENT_HEAD_TYPES);
  const auto& str = j.get<std::string>();
  for(size_t i = 0; i < num; ++i) {
    if(str.compare(th_type_names[i]) == 0) {
      type = static_cast<TreatmentHeadType>(i);
      return;
    }
  }
  throw nlohmann::detail::parse_error::create(302, 0, fmt::format("string '{}' is not valid TreatmentHeadType", str), j);
}

void from_json(const json& j, EnergyMode& type) {
  const size_t num = static_cast<size_t>(EnergyMode::NUM_ENERGY_MODES);
  const auto& str = j.get<std::string>();
  for(size_t i = 0; i < num; ++i) {
    if(str.compare(energy_mode_names[i]) == 0) {
      type = static_cast<EnergyMode>(i);
      return;
    }
  }
  throw nlohmann::detail::parse_error::create(302, 0, fmt::format("string '{}' is not valid EnergyMode", str), j);
}

void from_json(const json& j, MLCType& type) {
  const size_t num = static_cast<size_t>(MLCType::NUM_MLC);
  const auto& str = j.get<std::string>();
  for(size_t i = 0; i < num; ++i) {
    if(str.compare(mlc_type_names[i]) == 0) {
      type = static_cast<MLCType>(i);
      return;
    }
  }
  throw nlohmann::detail::parse_error::create(302, 0, fmt::format("string '{}' is not valid MLCType",str), j);
}

void from_json(const json& j, SourceType& type) {
  const size_t num = static_cast<size_t>(SourceType::NUM_SOURCE_TYPES);
  const auto& str = j.get<std::string>();
  for(size_t i = 0; i < num; ++i) {
    if(str.compare(source_type_names[i]) == 0) {
      type = static_cast<SourceType>(i);
      return;
    }
  }
  throw nlohmann::detail::parse_error::create(302, 0, fmt::format("string '{}' is not valid SourceType", str), j);
}

struct TreatmentHeadAndSourceInfo
{
  TreatmentHeadConfiguration th_config;
  SourceType source_type;
  using CustomSource = std::variant<GaussianRandomSource, PhaseSpaceSourceInfo, UniformPointSourceInfo>;
  CustomSource custom_source;
  std::optional<std::string> background_material;
};

GaussianRandomSource parseGaussianRandomSource(const json& obj);
PhaseSpaceSourceInfo parsePhaseSpaceSource(const json& obj);
UniformPointSourceInfo parseUniformPointSource(const json& obj);

TreatmentHeadAndSourceInfo parseTreatmentHeadAndSource(const json& obj) {
  TreatmentHeadAndSourceInfo th;
  th.source_type = SourceType::DEFAULT;
  {
    const auto& j = obj.at("treatment_head");
    th.th_config.th_type = j.at("type").get<TreatmentHeadType>();
    if(th.th_config.th_type == TreatmentHeadType::HALCYON) {
      th.th_config.mlc_type = MLCType::RDS;
      th.th_config.energy_mode = EnergyMode::X06_FFF;
    }
    if(th.th_config.th_type == TreatmentHeadType::AVALON) {
      th.th_config.mlc_type = MLCType::AVALON;
    }
    if(th.th_config.th_type == TreatmentHeadType::AVALON_ELECTRON) {
      th.th_config.mlc_type = MLCType::AVALON;
    }
    if(th.th_config.th_type == TreatmentHeadType::NOHEAD) {
      th.th_config.mlc_type = MLCType::NONE;
      th.th_config.energy_mode = EnergyMode::NONE;
      th.source_type = SourceType::PHASESPACE;
      th.background_material = "G4_Galactic";
    }
    if(j.contains("mlc")) {
      th.th_config.mlc_type = j["mlc"].get<MLCType>();
    }
    if(j.contains("energy_mode")) {
      th.th_config.energy_mode = j["energy_mode"].get<EnergyMode>();
    }
    if(j.contains("source")) {
      th.source_type = j["source"].get<SourceType>();
    }
    if(j.contains("background_material")) {
      th.background_material = j["background_material"].get<std::string>();
    }
  }
  if(th.source_type == SourceType::GAUSSIAN) {
    if(!obj.contains("gaussian_source"))
      throw std::runtime_error("No parameters provided for gaussian_source");
    th.custom_source = parseGaussianRandomSource(obj["gaussian_source"]);
  }
  if(th.source_type == SourceType::PHASESPACE) {
    if(!obj.contains("phasespace_source"))
      throw std::runtime_error("No parameters provided for phasespace_source");
    th.custom_source = parsePhaseSpaceSource(obj.at("phasespace_source"));
  }
  if(th.source_type == SourceType::UNIFORM_POINT) {
    if(!obj.contains("uniform_point_source"))
      throw std::runtime_error("No parameters provided for uniform_point_source");
    th.custom_source = parseUniformPointSource(obj.at("uniform_point_source"));
  }
  //Validate treatment head, energy mode and source combinations
  if(!isValidCombination(th.th_config))
    throw std::runtime_error("Treatment head configuration is not valid");
  return th;
}

GaussianRandomSource parseGaussianRandomSource(const json& obj) {
  GaussianRandomSource rs{};
  rs.energy = parseArrayWithUnit<1, EnergyUnit>(obj.at("energy"))[0];
  rs.energy_sigma = parseArrayWithUnit<1, EnergyUnit>(obj.at("energy_sigma"))[0];
  rs.spot_position = parseArrayWithUnit<2, LengthUnit>(obj.at("spot_position"));
  rs.spot_size = parseArrayWithUnit<2, LengthUnit>(obj.at("spot_size"));
  rs.angle = parseArrayWithUnit<2, AngleUnit>(obj.at("angle"));
  rs.angle_divergence = parseArrayWithUnit<2, AngleUnit>(obj.at("angle_divergence"));
  return rs;
}

PhaseSpaceSourceInfo parsePhaseSpaceSource(const json& obj) {
  PhaseSpaceSourceInfo pss{};
  pss.phsp_filename = obj.at("filename");
  return pss;
}

void from_json(const json& j, Particle::Type& t) {
  auto const& name = j.get<std::string>();
  t = Particle::UNKNOWN;
  for(size_t i = 0; i < Particle::NUM_TYPES; ++i) {
    auto cand = static_cast<Particle::Type>(i);
    if(name.compare(Particle::getParticleName(cand)) == 0)
      t = cand;
  }
  if(t == Particle::UNKNOWN)
    throw nlohmann::detail::parse_error::create(302, 0, fmt::format("String {} is not valid Particle::Type", name), j);
}

UniformPointSourceInfo parseUniformPointSource(const json& obj) {
  UniformPointSourceInfo ups{};
  ups.particle = obj.at("particle").get<Particle::Type>();
  ups.z = parseArrayWithUnit<1, LengthUnit>(obj.at("z"))[0];
  {
    auto const& fs = obj.at("field_size");
    if(fs.is_array()) {
      double symmetric_fs = parseArrayWithUnit<1, LengthUnit>(fs)[0];
       ups.x1 = -symmetric_fs/2.0;
       ups.x2 = symmetric_fs/2.0;
       ups.y1 = -symmetric_fs/2.0;
       ups.y2 = symmetric_fs/2.0;
    }
    else {
      auto lu = fs.at("unit").get<LengthUnit>();
      ups.x1 = lu.scaler * fs.at("x1").get<double>();
      ups.x2 = lu.scaler * fs.at("x2").get<double>();
      ups.y1 = lu.scaler * fs.at("y1").get<double>();
      ups.y2 = lu.scaler * fs.at("y2").get<double>();
    }
  }
  {
    auto const& spectrum = obj.at("spectrum");
    EnergyUnit eu = spectrum.at("unit").get<EnergyUnit>();
    ups.spectrum_energies = spectrum.at("energies").get<std::vector<double>>();
    ups.spectrum_intensities = spectrum.at("intensities").get<std::vector<double>>();
    for(size_t i = 0; i < ups.spectrum_energies.size(); ++i) {
      ups.spectrum_energies[i] *= eu.scaler;
    }
  }
  if(ups.spectrum_energies.size() != ups.spectrum_intensities.size())
    throw nlohmann::detail::parse_error::create(302, 0, "spectrum_energies and spectrum_intensities have different lengths", obj);
  return ups;
}

DoseGrid parseDoseGrid(const json& obj) {
  DoseGrid dg{};
  dg.filename = obj.at("filename").get<std::string>();
  dg.center_position = parseArrayWithUnit<3, LengthUnit>(obj.at("center_position"));
  dg.size = parseArrayWithUnit<3, LengthUnit>(obj.at("size"));
  dg.num_voxels = parseArray<3, int>(obj.at("num_voxels"));
  dg.binary_output = true;
  if(obj.contains("binary_output"))
    dg.binary_output = obj.at("binary_output").get<bool>();
  if (obj.contains("rotation")) {
    const auto& j = obj.at("rotation");
    Rotation rot;
    rot.axis = parseArray<3, double>(j.at("vector"));
    rot.angle = parseArrayWithUnit<1, AngleUnit>(j.at("angle"))[0];
    dg.rotation = rot;
  }
  return dg;
}

void from_json(const json& j, PhantomType& type) {
  const size_t num = static_cast<size_t>(PhantomType::NUM_PHANTOM_TYPES);
  const auto& str = j.get<std::string>();
  for(size_t i = 0; i < num; ++i) {
    if(str.compare(phantom_type_names[i]) == 0) {
      type = static_cast<PhantomType>(i);
      return;
    }
  }
  throw nlohmann::detail::parse_error::create(302, 0, fmt::format("string '{}' is not valid PhantomType", str), j);
}

void from_json(const json& j, Axis& type) {
  const size_t num = static_cast<size_t>(Axis::NUM_AXES);
  const auto& str = j.get<std::string>();
  for(size_t i = 0; i < num; ++i) {
    if(str.compare(axis_names[i]) == 0) {
      type = static_cast<Axis>(i);
      return;
    }
  }
  throw nlohmann::detail::parse_error::create(302, 0, fmt::format("string '{}' is not valid Axis", str), j);
}

void from_json(const json& j, BiasingMode& mode) {
  const size_t num = static_cast<size_t>(BiasingMode::NUM_BIASING_MODES);
  const auto& str = j.get<std::string>();
  for(size_t i = 0; i < num; ++i) {
    if(str.compare(biasing_mode_names[i]) == 0) {
      mode = static_cast<BiasingMode>(i);
      return;
    }
  }
  throw nlohmann::detail::parse_error::create(302, 0, fmt::format("string '{}' is not valid BiasingMode", str), j);
}

WaterBoxPhantom parseWaterBoxPhantom(const json& obj) {
  WaterBoxPhantom p{};
  p.center_position = parseArrayWithUnit<3, LengthUnit>(obj.at("center_position"));
  p.size = parseArrayWithUnit<3, LengthUnit>(obj.at("size"));
  return p;
}

WaterCylinderPhantom parseWaterCylinderPhantom(const json& obj) {
  WaterCylinderPhantom p;
  p.center_position = parseArrayWithUnit<3, LengthUnit>(obj.at("center_position"));
  p.height_axis = obj.at("height_axis").get<Axis>();
  p.height = parseArrayWithUnit<1, LengthUnit>(obj.at("height"))[0];
  p.radius = parseArrayWithUnit<1, LengthUnit>(obj.at("radius"))[0];
  return p;
}

VoxelBoxPhantom parseVoxelBoxPhantom(const json& obj) {
  VoxelBoxPhantom p;
  p.voxel_sizes = parseArrayWithUnit<3, LengthUnit>(obj.at("sizes"));
  p.voxel_counts = parseArray<3, size_t>(obj.at("counts"));
  p.center_position = parseArrayWithUnit<3, LengthUnit>(obj.at("center_position"));
  //params.geant4_commands = j.at("geant_commands").get<std::vector<std::string>>();
  p.material_idxs = obj.at("materials").get<std::vector<size_t>>();
  if(obj.contains("densities")) {
    p.densities = obj.at("densities").get<std::vector<double>>();
    if(p.densities.size() != p.material_idxs.size())
      throw nlohmann::detail::parse_error::create(302, 0, "material indices and density counts should match", obj);
  }
  else {
    p.densities = std::vector<double>(p.material_idxs.size(), 0.0);
  }
  if(p.voxel_counts[0]*p.voxel_counts[1]*p.voxel_counts[2] != p.material_idxs.size())
    throw nlohmann::detail::parse_error::create(302, 0, "material indices and voxel counts should match", obj);
  {
    const auto& material_dict_json = obj.at("material_dict");
    for(const auto& mat : material_dict_json.items()) {
      size_t idx = std::stoi(mat.key());
      p.idx_to_material[idx] = VoxelBoxPhantom::UserMaterial{};
      VoxelBoxPhantom::Material& material = p.idx_to_material[idx];
      if(mat.value().is_string()) {
        material = mat.value().get<std::string>();
      }
      else if(mat.value().is_object()) {
        const auto& j_it = mat.value();
        VoxelBoxPhantom::UserMaterial user_mat;
        if(!j_it.size())
          throw nlohmann::detail::parse_error::create(302, 0, "material object should have at least one component", j_it);
        user_mat.default_density = j_it.at("default_density").get<double>() * g / cm3;
        double total_mass_fraction = 0.0;
        for(const auto& component : j_it.items()) {
          if(component.key().compare("default_density") == 0)
            continue;
          user_mat.names.emplace_back(component.key());
          const double mass_fraction = component.value().get<double>();
          user_mat.mass_fractions.emplace_back(mass_fraction);
          total_mass_fraction += mass_fraction;
        }
        if(abs(total_mass_fraction - 1.0) > 1e-15)
          throw nlohmann::detail::parse_error::create(302, 0, "material mass fractions don't add up to one", j_it);
        material = user_mat;
      }
      else {
        throw nlohmann::detail::parse_error::create(302, 0, "material is not string or material object", mat.value());
      }
    }
  }
  //Validate as much as possible
  for(size_t i = 0; i < p.densities.size(); ++i) {
    if(p.densities[i] < 0.0)
      throw std::runtime_error(fmt::format("Density provided in VoxelBoxPhantom is negative at position {}", i));
  }
  return p;
}

Phantom parsePhantom(const json& obj) {
  Phantom p;
  p.type = obj.at("type").get<PhantomType>();
  if(p.type == PhantomType::WATER_BOX)
    p.def = parseWaterBoxPhantom(obj);
  else if(p.type == PhantomType::WATER_CYLINDER)
    p.def = parseWaterCylinderPhantom(obj);
  else if(p.type == PhantomType::VOXEL_BOX)
    p.def = parseVoxelBoxPhantom(obj);
  else {
    static_assert(static_cast<size_t>(PhantomType::NUM_PHANTOM_TYPES) == 3);
  }
  return p;
}

PhaseSpacePlane parsePhaseSpacePlane(const json& obj) {
  PhaseSpacePlane pl{};
  pl.phsp_filename = obj.at("phsp_filename").get<std::string>();
  pl.center_position = parseArrayWithUnit<3, LengthUnit>(obj.at("center_position"));
  pl.size = parseArrayWithUnit<3, LengthUnit>(obj.at("size"));
  pl.kill_particles = obj.at("kill_particles").get<bool>();
  return pl;
}

std::unique_ptr<PhaseSpaceSourceData> readPhaseSpaceFile(fs::path filename) {
  std::ifstream file(filename, std::ios::binary);
  const size_t file_size = fs::file_size(filename);
  file.seekg(0, std::ios::beg);
  if (file_size < sizeof(PhaseSpaceDataHeader))
    throw std::runtime_error("PhaseSpace file does not contain header");
  PhaseSpaceDataHeader header{};
  file.read(reinterpret_cast<char*>(&header), sizeof(header));
  if (header.magic_number != g_PHASESPACE_HEADER_MAGIC)
    throw std::runtime_error("PhaseSpace file header is not valid");
  if (header.version != 1)
    throw std::runtime_error("PhaseSpace file version is not supported");
  size_t expected_stride = sizeof(PhaseSpaceParticle);
  if (header.has_tracking_information)
    expected_stride += sizeof(PhaseSpaceTrack);
  if (header.data_stride != expected_stride)
    throw std::runtime_error("PhaseSpace file data stride is not consistent with structs");
  if (header.data_size % header.data_stride != 0)
    throw std::runtime_error("Phase Space file particle data is not consistent");
  const size_t num_particles = header.data_size / header.data_stride;
  auto ptr = std::make_unique<PhaseSpaceSourceData>();
  auto& data = *ptr;
  data.data_buffer.resize(num_particles*header.data_stride);
  file.seekg(header.data_offset, std::ios::beg);
  file.read(reinterpret_cast<char*>(data.data_buffer.data()), header.data_stride * num_particles);
  data.num_particles = num_particles;
  data.num_histories = header.num_particle_histories;
  data.particle_data_offset = static_cast<uint32_t>(header.particle_offset);
  data.track_data_offset = header.track_offset;
  data.has_tracking_information = static_cast<bool>(header.has_tracking_information);
  data.stride = header.data_stride;
  return ptr;
}

SimulationInputParams parseSimulationParamsWithDefaults(const json& obj)
{
  SimulationInputParams params;
  std::random_device r;
  //defaults
  params.num_threads = std::thread::hardware_concurrency();
  params.biasing_mode = BiasingMode::NONE;
  params.seed = r();
  params.track_particles = false;
  params.physics_model = "QGSP_BIC_EMZ";

  const auto& j = obj;
  if(j.contains("seed"))
    params.seed = j.at("seed").get<int>();
  if(j.contains("physics_model"))
    params.physics_model = j.at("physics_model").get<std::string>();
  if(j.contains("num_threads"))
    params.num_threads = j.at("num_threads").get<int>();
  if(j.contains("track_particles"))
    params.track_particles = j.at("track_particles").get<bool>();
  if(j.contains("biasing_mode"))
    params.biasing_mode = j.at("biasing_mode").get<BiasingMode>();
  if (j.contains("geant_commands"))
      params.geant4_commands = j.at("geant_commands").get<std::vector<std::string>>();
  if(j.contains("phantom")) {
    params.phantom = parsePhantom(j.at("phantom"));
  }
  return params;
}

SimulationOutputParams parseSimulationOutputParams(const json& obj)
{
  SimulationOutputParams params;
  if (obj.contains("dose_grids"))
  {
    const auto& grids = obj.at("dose_grids");
    for (size_t i = 0; i < grids.size(); ++i) {
      params.dose_grids.emplace_back(parseDoseGrid(grids.at(i)));
    }
  }
  if (obj.contains("phase_space_recorders"))
  {
    const auto phsps = obj.at("phase_space_recorders");
    for (size_t i = 0; i < phsps.size(); ++i) {
      const auto& phsp = phsps.at(i);
      const auto& str = phsp.at("type").get<std::string>();
      if(str.compare("plane") != 0)
          throw nlohmann::detail::parse_error::create(302, 0, std::string("Unknown PhaseSpace recorder type: ") + str, phsps.at("type"));
     params.phase_space_plane.emplace_back(parsePhaseSpacePlane(phsp));
    }
  }
  return params;
}


int doLinacSimulation(
  const CommandLineAndEnvironmentArguments& args,
  RunSimulation run_simulation,
  ILogger* logger) 
{
  auto& log = *logger;
  try {
  Timing timing{};
  const auto start_total = std::chrono::high_resolution_clock::now();
  if(!fs::exists(args.virtuallinac_root))
    throw std::runtime_error(fmt::format("Virtuallinac root location {} does not exist", args.virtuallinac_root.string()));
  json simulation_config = json::parse(args.input_json, nullptr, true, true);
  if(args.additional_input_json.has_value()) {
    json additional_config = json::parse(args.additional_input_json.value(), nullptr, true, true);
    simulation_config.merge_patch(additional_config);
  }
  //Set input parameters
  uint64_t num_particles = 0;
  uint64_t num_simulation_points = 1000;
  if(simulation_config.contains("num_particles"))
    num_particles = simulation_config.at("num_particles").get<uint64_t>();
  num_particles = args.num_particles.value_or(num_particles);
  uint64_t num_original_histories = num_particles;
  if(simulation_config.contains("num_simulation_points"))
    num_simulation_points = simulation_config.at("num_simulation_points").get<uint64_t>();

  SimulationInputParams input_params = parseSimulationParamsWithDefaults(simulation_config);
  SimulationOutputParams output_params = parseSimulationOutputParams(simulation_config);
  TreatmentHeadAndSourceInfo th_source_info = parseTreatmentHeadAndSource(simulation_config);
  //Incorporate command line information
  //commandline overrides all other parameters by default
  input_params.num_threads = args.num_threads.value_or(input_params.num_threads);
  input_params.seed = args.seed.value_or(input_params.seed);
  if (args.phsp_source_file.has_value()) {
    th_source_info.source_type = SourceType::PHASESPACE;
    th_source_info.custom_source = PhaseSpaceSourceInfo{args.phsp_source_file.value()};
  }
  //Handle output redirecting here if necessary!
  if (args.output_directory.has_value()) {
    fs::path out_dir = args.output_directory.value();
    fs::create_directories(out_dir);
    for (auto& dosegrid : output_params.dose_grids) {
      dosegrid.filename = (out_dir / dosegrid.filename).generic_u8string();
    }
    for (auto& phsp_plane : output_params.phase_space_plane) {
      phsp_plane.phsp_filename = (out_dir / phsp_plane.phsp_filename).generic_u8string();
    }
  }
  auto gdml_path = args.virtuallinac_root / "gdml";
  auto stl_path = args.virtuallinac_root / "stls";
  //Setup source and treatment head
  //Combinations should be validated in parse so this should result in valid treatment head and source
  std::unique_ptr<LinacSource> linac_source;
  std::unique_ptr<LinacDetector> treatment_head;
  if(th_source_info.th_config.th_type == TreatmentHeadType::TRUEBEAM) {
    if(th_source_info.source_type == SourceType::DEFAULT) {
      if(th_source_info.th_config.energy_mode == EnergyMode::X06_FFF || 
          th_source_info.th_config.energy_mode == EnergyMode::X06)
        linac_source = getTrueBeam6MeVSource();
    }
    treatment_head = std::make_unique<truebeam::TreatmentHeadDetector>(
      th_source_info.th_config.mlc_type, th_source_info.th_config.energy_mode,
      "sd_monitor_chamber", gdml_path, stl_path
    );
  }
  else if(th_source_info.th_config.th_type == TreatmentHeadType::HALCYON) {
    if(th_source_info.source_type == SourceType::DEFAULT)
      linac_source = getHalcyonSource();
    treatment_head = std::make_unique<halcyon::TreatmentHeadDetector>(
      "sd_monitor_chamber", gdml_path, stl_path
    );
  }
  else if(th_source_info.th_config.th_type == TreatmentHeadType::AVALON) {
    if(th_source_info.th_config.energy_mode == EnergyMode::X06_FFF || 
          th_source_info.th_config.energy_mode == EnergyMode::X06)
        linac_source = getTrueBeam6MeVSource();
    else if(th_source_info.th_config.energy_mode == EnergyMode::X10_FFF)
        linac_source = getAvalon10XFFFSource();
    else if (th_source_info.th_config.energy_mode == EnergyMode::X10)
        linac_source = getAvalon10XSource();
    treatment_head = std::make_unique<avalon::TreatmentHeadDetector>(th_source_info.th_config.energy_mode,
      "sd_monitor_chamber", gdml_path, stl_path
    );
  }

  else if(th_source_info.th_config.th_type == TreatmentHeadType::AVALON_ELECTRON) {
    if(th_source_info.th_config.energy_mode == EnergyMode::E09) {
        linac_source = getAvalon09ESource();
    }
    else if(th_source_info.th_config.energy_mode == EnergyMode::E15) {
        linac_source = getAvalon15ESource();
    }
    treatment_head = std::make_unique<avalon_electron::TreatmentHeadDetector>(th_source_info.th_config.energy_mode,
    //treatment_head = std::make_unique<avalon::TreatmentHeadDetector>(th_source_info.th_config.energy_mode,
      "sd_monitor_chamber", gdml_path, stl_path
    );
  }

  else if(th_source_info.th_config.th_type == TreatmentHeadType::NOHEAD) {
    treatment_head = std::make_unique<nohead::TreatmentHeadDetector>(
      th_source_info.background_material.value()
    );
  }
  else {
    static_assert(truebeam::mlcs.size() == 3);
    throw std::runtime_error("Unknown treatment head type");
  }
  if(th_source_info.source_type == SourceType::GAUSSIAN) {
    linac_source = std::make_unique<GaussianElectronSource>(
      std::get<GaussianRandomSource>(th_source_info.custom_source));
  }
  else if (th_source_info.source_type == SourceType::PHASESPACE) {
    const auto& phps_info = std::get<PhaseSpaceSourceInfo>(th_source_info.custom_source);
    auto phsp_data = readPhaseSpaceFile(phps_info.phsp_filename);
    num_particles = phsp_data->num_particles;
    num_original_histories = phsp_data->num_histories;
    linac_source = std::make_unique<PhaseSpaceSource>(std::move(phsp_data));
  }
  else if(th_source_info.source_type == SourceType::UNIFORM_POINT) {
    const auto& info = std::get<UniformPointSourceInfo>(th_source_info.custom_source);
    linac_source = std::make_unique<LinacUniformPointSource>(info);
  }
  if(linac_source == nullptr)
    throw std::logic_error("Linac source not initialized");
  //Parse plan. Num particles should be set here. After this plan contains correct
  //information on number of particles
  std::unique_ptr<Plan> plan = nullptr;
  {
    if(num_particles == 0)
      throw std::runtime_error("No particle number is provided");
    bool has_control_points = simulation_config.contains("control_points");
    bool has_simulation_points = simulation_config.contains("simulation_points");
    if(!(has_control_points || has_simulation_points))
      throw std::runtime_error("Simulation should have at least one control/simulation point");
    if(has_control_points && has_simulation_points)
      throw std::runtime_error("Simulation should have either control points or simulation points but not both");
    if(has_control_points) {
      plan = Plan::parseControlPointPlan(
        th_source_info.th_config.th_type,
        th_source_info.th_config.mlc_type,
        num_particles,
        simulation_config.at("control_points"),
        log,
        static_cast<uint64_t>(input_params.seed),
        num_simulation_points);
    }
    else if(has_simulation_points) {
      plan = Plan::parseSimulationPointPlan(
        th_source_info.th_config.th_type,
        th_source_info.th_config.mlc_type,
        num_particles,
        simulation_config.at("simulation_points"),
        log
      );
    }
    else {
      throw std::runtime_error("No simulation points or control_points found");
    }
  }
  //Correct the number of histories (planning might have changed total particle numbers)
  num_original_histories = static_cast<uint64_t>((static_cast<double>(plan->getTotalNumParticles())/num_particles) * num_original_histories);
  
  log.info(fmt::format("Number of particles: {}", num_particles));
  log.info(fmt::format("Number of histories: {}", num_original_histories));
  log.info(fmt::format("Number of threads: {}", input_params.num_threads));
  log.info(fmt::format("Number of simulation points: {}", plan->getNumSimulationPoints()));

  //Prepare Geant4 for simulation
  const auto config_duration = std::chrono::high_resolution_clock::now() - start_total;
  timing.configuration = std::chrono::duration_cast<std::chrono::milliseconds>(config_duration).count();
  //Run Geant4 simulation
  bool start_simulation = !args.dry_run;
  if(start_simulation) {
    run_simulation(
      input_params, output_params, treatment_head.get(),
      linac_source.get(), *plan, num_original_histories,
      args.visualize, timing);
    const auto total_duration = std::chrono::high_resolution_clock::now() - start_total;
    timing.total = std::chrono::duration_cast<std::chrono::milliseconds>(total_duration).count();
    log.log("Simulation finished\n")
      .log("Timings (milliseconds)\n")
      .log(fmt::format("\tConfiguration:  {}\n", timing.configuration))
      .log(fmt::format("\tInitialization: {}\n", timing.initialization))
      .log(fmt::format("\tTH movement:    {}\n", timing.state_transitions))
      .log(fmt::format("\tBeam on:        {}\n", timing.beam_on))
      .log(fmt::format("\tTotal:          {}\n", timing.total));
  }
  else {
    log.info("Dry run complete. No simulation is run");
  } 
  {
    auto detector = treatment_head.release();
    //delete crashes for unknown reason, probably related to G4
    //delete detector; 
  }
  }
  catch (const std::exception& e) {
    log.error(e.what()).info("Aborting");
    return -1;
  }
  return 0;
}
