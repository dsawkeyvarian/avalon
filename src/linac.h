#pragma once
#include "machine_information.h"
#include <string>
#include <optional>
#include <variant>
#include <array>
#include <unordered_map>
#include <functional>
#include <filesystem>
namespace fs = std::filesystem;

struct CommandLineAndEnvironmentArguments
{
  std::string input_json;
  std::optional<std::string> additional_input_json;
  fs::path virtuallinac_root;
  std::optional<int> seed;
  std::optional<int> num_threads;
  std::optional<size_t> num_particles;
  std::optional<std::string> phsp_source_file;
  std::optional<fs::path> output_directory;
  bool dry_run = false;
  bool visualize = false;
};

struct SimulationInputParams;
struct SimulationOutputParams;
class LinacDetector;
class LinacSource;
class Plan;
struct Timing;
class ILogger;

using RunSimulation = std::function<void(
  const SimulationInputParams&,
  const SimulationOutputParams&,
  LinacDetector*,
  LinacSource*,
  const Plan&,
  uint64_t,
  bool,
  Timing&)>;
int doLinacSimulation(
  const CommandLineAndEnvironmentArguments& args,
  RunSimulation run_simulation,
  ILogger* log);

struct Vec3
{
  double x;
  double y;
  double z;

  Vec3() : x(0), y(0), z(0) {}

  Vec3(const std::array<double, 3> arr)
    : x(arr[0]), y(arr[1]), z(arr[2]) {}

  Vec3(double _x, double _y, double _z)
    : x(_x), y(_y), z(_z) {}
};

struct Vec2
{
  double x;
  double y;

  Vec2() : x(0), y(0) {}

  Vec2(const std::array<double, 2> arr)
    : x(arr[0]), y(arr[1]) {}

  Vec2(double _x, double _y)
    : x(_x), y(_y) {}
};

struct Particle {
  enum Type {
    UNKNOWN = 0,
    PHOTON = 1,
    ELECTRON,
    POSITRON,
    PROTON,
    NEUTRON,
    NUM_TYPES
  };
  static Type getParticleType(int particlePDGEncoding) {
    switch (particlePDGEncoding) {
    case 22: return PHOTON;
    case 11: return ELECTRON;
    case -11: return POSITRON;
    case 2212: return PROTON;
    case 2112: return NEUTRON;
    default:
      throw std::runtime_error("Unknown PDGE particle encoding: " + std::to_string(particlePDGEncoding));
    }
  }
  static int getPDGEncoding(Type _type) {
    switch (_type) {
    case PHOTON:
      return 22;
    case ELECTRON:
      return 11;
    case POSITRON:
      return -11;
    case PROTON:
      return 2212;
    case NEUTRON:
      return 2112;
    default:
      throw std::runtime_error("Unknown particle type. Unable to find PDF encoding");
    }
  }
  static const char* getParticleName(Type _type) {
    switch (_type) {
    case PHOTON:
      return "PHOTON";
    case ELECTRON:
      return "ELECTRON";
    case POSITRON:
      return "POSITRON";
    case PROTON:
      return "PROTON";
    case NEUTRON:
      return "NEUTRON";
    case UNKNOWN:
      return "UNKNOWN";
    default:
      throw std::runtime_error("Unknown particle type: " + std::to_string(_type));
    }
  }
};

struct Rotation
{
  Vec3 axis;
  double angle;
};

struct DoseGrid
{
  std::string filename;
  bool binary_output;
  Vec3 center_position;
  Vec3 size;
  std::array<int, 3> num_voxels;
  std::optional<Rotation> rotation;
};

struct PhaseSpacePlane
{
  std::string phsp_filename;
  Vec3 center_position;
  Vec3 size;
  bool kill_particles;
};

enum class PhantomType {
  WATER_BOX,
  WATER_CYLINDER,
  VOXEL_BOX,
  NUM_PHANTOM_TYPES
};

static const char* phantom_type_names[] = {
  "water_box",
  "water_cylinder",
  "voxel_box"
};

static_assert(static_cast<size_t>(PhantomType::NUM_PHANTOM_TYPES) == SIZEOFARRAY(phantom_type_names));

enum class Axis {
  X, Y, Z,
  NUM_AXES
};

static const char* axis_names[] = {
  "x", "y", "z"
};

static_assert(static_cast<size_t>(Axis::NUM_AXES) == SIZEOFARRAY(axis_names));


struct WaterBoxPhantom {
  Vec3 center_position;
  Vec3 size;
};

struct WaterCylinderPhantom {
  Axis height_axis;
  double radius;
  double height;
  Vec3 center_position;
};

struct VoxelBoxPhantom {
  struct UserMaterial {
    std::vector<std::string> names;
    std::vector<double> mass_fractions;
    double default_density;
  };
  using Material = std::variant<std::string, UserMaterial>;
  Vec3 voxel_sizes;
  std::array<size_t, 3> voxel_counts;
  Vec3 center_position;
  std::vector<size_t> material_idxs;
  std::vector<double> densities;
  std::unordered_map<size_t, Material> idx_to_material;
};

struct Phantom {
  PhantomType type;
  using PhantomDefinition = std::variant<WaterBoxPhantom, WaterCylinderPhantom, VoxelBoxPhantom>;
  PhantomDefinition def;
};

struct GaussianRandomSource
{
  double energy;
  double energy_sigma;
  Vec2 spot_position;
  Vec2 spot_size;
  Vec2 angle;
  Vec2 angle_divergence;
};

struct UniformPointSourceInfo
{
  Particle::Type particle;
  double z;
  double x1;
  double x2;
  double y1;
  double y2;
  std::vector<double> spectrum_energies;
  std::vector<double> spectrum_intensities;
};

struct PhaseSpaceSourceInfo
{
  std::string phsp_filename;
};

enum class BiasingMode {
  NONE,
  TARGET,
  LIGHTFIELD,
  DIRECTIONAL_SPLITTING,
  NUM_BIASING_MODES
};

static const char* biasing_mode_names[] = {
  "none",
  "target",
  "lightfield",
  "directional_splitting"
};
static_assert(static_cast<size_t>(BiasingMode::NUM_BIASING_MODES) == SIZEOFARRAY(biasing_mode_names));

enum SourceType {
  DEFAULT,
  GAUSSIAN,
  PHASESPACE,
  UNIFORM_POINT,
  NUM_SOURCE_TYPES
};

static const char* source_type_names[] = {
  "default",
  "gaussian",
  "phasespace",
  "uniform_point"
};

static_assert(static_cast<size_t>(SourceType::NUM_SOURCE_TYPES) == SIZEOFARRAY(source_type_names));

struct SimulationInputParams
{
  int seed;
  int num_threads;
  std::string physics_model;
  BiasingMode biasing_mode;
  std::vector<std::string> geant4_commands;
  bool track_particles;
  std::optional<Phantom> phantom;
};

struct SimulationOutputParams
{
  std::vector<DoseGrid> dose_grids;
  std::vector<PhaseSpacePlane> phase_space_plane;
};

struct PhaseSpaceSourceData {
  uint64_t num_particles;
  uint64_t num_histories;

  std::vector<uint8_t> data_buffer;
  uint32_t stride;
  uint32_t particle_data_offset;
  uint32_t track_data_offset;
  bool has_tracking_information;
};

struct Timing {
  uint64_t total;
  uint64_t configuration;
  uint64_t initialization;
  uint64_t state_transitions;
  uint64_t beam_on;
};

//These values are used for dual-purpose:
//1. For geometry traversal tracking
//2. For indicating creation geometry of a particles
//OTHER and SOURCE are used for 2. but *not* used for 1.
//NONE is not used for 2
enum class TraversedGeometry : uint32_t {
  TARGET = 0,
  PRIMARY_COLLIMATOR = 1,
  SHIELD_COLLIMATOR = 2,
  YSTAGE_SHIELD = 3,
  FLATTENING_FILTER = 4,
  IC = 5,
  UPPER_JAW = 6,
  LOWER_JAW = 7,
  BASEPLATE = 8,
  BASEPLATE02 = 9,
  EFOIL1 = 10,
  EFOIL2 = 11,
  COLL_PLATE = 12,
  MLC = 13,
  SOURCE = 14,
  SCRAPER1 = 15,
  SCRAPER1_SUPPORT = 16,
  SCRAPER2 = 17,
  SCRAPER3 = 18,
  OTHER = 19,
  NUM_TRAVERSED,
  NONE
};

static const char* g_traversed_geometry_names[] = {
  "TARGET",
  "PRIMARY_COLLIMATOR",
  "SHIELD_COLLIMATOR",
  "YSTAGE_SHIELD",
  "FLATTENING_FILTER",
  "IC",
  "UPPER_JAW",
  "LOWER_JAW",
  "BASEPLATE",
  "BASEPLATE02",
  "EFOIL1",
  "EFOIL2",
  "COLL_PLATE",
  "MLC",
  "SOURCE",
  "SCRAPER1",
  "SCRAPER1_SUPPORT",
  "SCRAPER2",
  "SCRAPER3",
  "OTHER"
};

static_assert(static_cast<size_t>(TraversedGeometry::NUM_TRAVERSED) == SIZEOFARRAY(g_traversed_geometry_names));

struct VolumeNameAndTraversalFlag {
  char name[50];
  TraversedGeometry traversed;
};

enum class PhysicsProcess : uint16_t {
  UNKNOWN = 0,
  PRIMARY,
  PHOTON_PHOT,
  PHOTON_COMPT,
  PHOTON_CONV,
  PHOTON_RAYL,
  ELECTRON_MSC,
  ELECTRON_IONI,
  ELECTRON_EBREM,
  ELECTRON_PAIRPROD,
  ELECTRON_COULOMBSCAT,
  NUM_PROCESSES
};

static const char* g_physics_process_names[] = {
  "DUMMY0",
  "DUMMY1",
  "phot",
  "compt",
  "conv",
  "Rayl",
  "msc",
  "eIoni",
  "eBrem",
  "ePairProd",
  "CoulombScat"
};

static_assert(static_cast<size_t>(PhysicsProcess::NUM_PROCESSES) == SIZEOFARRAY(g_physics_process_names));

//Binary format of the output files
//Member order is important, do not touch

const uint32_t g_PHASESPACE_HEADER_MAGIC = 0xAFAFA0A0u;


//PhaseSpace file has a header. After header (data offset) it contains array of PhaseSpaceParticle+PhaseSpaceTrack
//where PhaseSpaceTrack is optional (has_tracking_information).
struct PhaseSpaceParticle {
  uint32_t particle_type;
  float energy;
  float x, y, z;
  float mom_x, mom_y;
  float weight;
};

struct PhaseSpaceTrack {
  uint32_t traversed_geometry; //Add stuff already here
  uint16_t creation_process;
  uint16_t creation_geometry;
};

struct PhaseSpaceDataHeader {
  uint32_t magic_number; // magic stamp to check that file is "valid"
  uint32_t version;
  uint64_t data_offset;
  uint64_t data_size;
  uint32_t has_tracking_information;
  //particle_offset and track_offset relative to data_offset
  uint32_t particle_offset;
  uint32_t track_offset;
  uint32_t data_stride;
  uint64_t num_particle_histories;
  char hash[64];
  char geant_version[64];
};

//Dose file has a header. After the header (data_offset) dose data is a double array.
struct DoseBinaryHeader {
  uint32_t version;
  uint32_t data_offset;
  uint32_t num_bins[3];
  uint32_t pad0; //Align to 8-bytes for double
  double size_in_mm[3];
  double center_position_in_mm[3];
  uint32_t is_rotated;
  uint32_t pad1; //Align to 8-bytes for double
  double rotation_axis[3];
  double rotation_angle_in_rad;
  double monitor_chamber_dose;
  char hash[64];
  char geant_version[64];
};
