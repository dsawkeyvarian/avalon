#include "logger.h"
#include "linac.h"
#include "plan.h"
#include "test.h"

#include "external/span.hpp"
#include <fmt/core.h>
#include <fmt/format.h>
#include <optional>
#include <fstream>
#define JSON_DIAGNOSTICS 1
#include "external/json.hpp"
using json = nlohmann::json;

#define DOCTEST_CONFIG_IMPLEMENT
#include "external/doctest.h"

static std::optional<fs::path> g_virtuallinac_root = {};
static std::optional<fs::path> g_output_folder = {};

fs::path getVirtualLinacRoot() {
  return g_virtuallinac_root.value();
}

fs::path getOutputFolder() {
  if(g_output_folder.has_value())
    return g_output_folder.value();
  throw std::runtime_error("Define output folder for this test");
}

int main(int argc, char** argv) {
  doctest::Context context;
  context.applyCommandLine(argc, argv);
  int custom_argument_idx = -1;
  for(int i = 0; i < argc; ++i) {
    if(strcmp(argv[i], "--") == 0)
      custom_argument_idx = i;
  }
  if(custom_argument_idx > 0) {
    if(custom_argument_idx == argc - 1) {
      std::cout << "Custom arguments for tests:" << std::endl;
      std::cout << "  --root <virtuallinac_root_path>: Path where to find stl/gdml/examples" << std::endl;
      std::cout << "  --out <output_path>: Directory for test outputs (doses/phasespaces)" << std::endl;
      return 0;
    }
    for(int i = custom_argument_idx+1; i < argc; ++i) {
      if(strcmp(argv[i], "--root") == 0) {
        if(i + 1 == argc) {
          std::cerr << "No value for --out" << std::endl;
          return -1;
        }
        g_virtuallinac_root = fs::path(argv[i+1]);
        i += 1;
      }
      else if(strcmp(argv[i], "--out") == 0) {
        if(i + 1 == argc) {
          std::cerr << "No value for --out" << std::endl;
          return -1;
        }
        g_output_folder = fs::path(argv[i+1]);
        i += 1;
      }
    }
  }
  const char* virtuallinac_root_env = getenv("VIRTUALLINAC_ROOT");
  if (!virtuallinac_root_env && !g_virtuallinac_root.has_value()) {
    std::cerr << "Define VIRTUALLINAC_ROOT and point it to VirtualLinac repository or set --root argument. VirtualLinac needs the root path to find .gdml and .stl files." << std::endl;
    std::cerr << "Aborting" << std::endl;
    return -1;
  }
  if(!g_virtuallinac_root.has_value() && (virtuallinac_root_env != nullptr))
    g_virtuallinac_root = virtuallinac_root_env;
  int res = context.run();
  return res;
}

TEST_CASE("Examples"
  * doctest::description("Dry run all files in examples folder"))
{
  auto root = getVirtualLinacRoot();
  auto examples_folder = root / "examples";

  CommandLineAndEnvironmentArguments args;
  args.virtuallinac_root = root;
  args.dry_run = true;

  for(const auto& entry : fs::recursive_directory_iterator(examples_folder)) {
    if(entry.is_directory())
      continue;
    auto filepath = entry.path();
    auto extension = filepath.extension();
    bool validInputFile = (extension.compare(".json") == 0) || (extension.compare(".jsonc") == 0);
    if(!validInputFile)
      continue;
    std::ifstream input(filepath);
    CHECK_MESSAGE(input.good(), fmt::format("Unable to read file {}", filepath.string()));
    std::stringstream buffer;
    buffer << input.rdbuf();
    args.input_json = buffer.str();
    SinkLogger log;
    int return_value = doLinacSimulation(args, nullptr, &log);
    CHECK_MESSAGE(return_value == 0, fmt::format("Input JSON file {} failed", filepath.string()));
    INFO(fmt::format("Input JSON file {} passed", filepath.string()));
  }
}

TEST_CASE("CommandLineOverride")
{
  auto root = getVirtualLinacRoot();
  auto examples_folder = root / "examples";
  auto filename = examples_folder / "truebeam.jsonc";
  std::ifstream input(filename);
  if(!input.good())
    throw std::runtime_error(fmt::format("Unable to read file {}", filename.string()));
  std::stringstream buffer;
  buffer << input.rdbuf();

  CommandLineAndEnvironmentArguments args;
  args.virtuallinac_root = root;
  args.num_particles = 20;
  args.num_threads = 2;
  args.seed = 643576;
  args.dry_run = false;
  args.input_json = buffer.str();

  auto runTestSimulation = [=](const SimulationInputParams& in_params, const SimulationOutputParams&,
    LinacDetector*, LinacSource*, const Plan& plan, uint64_t num_histories, bool, Timing&) {
      CHECK(plan.getTotalNumParticles() == args.num_particles);
      CHECK(num_histories == args.num_particles);
      CHECK(in_params.seed == args.seed);
      CHECK(in_params.num_threads == args.num_threads);
      CHECK(in_params.track_particles == true);
      CHECK(in_params.phantom.value().type == PhantomType::WATER_BOX);
      CHECK(plan.getNumSimulationPoints() == 1);
      CHECK(plan.getNumParticles(0) == args.num_particles);
    };
  SinkLogger log;
  CHECK(doLinacSimulation(args, runTestSimulation, &log) == 0);
}

TEST_CASE("ControlPointInput") {
  auto root = getVirtualLinacRoot();
  auto examples_folder = root / "examples";
  auto filename = examples_folder / "sweeping_gap.jsonc";
  std::ifstream input(filename);
  if(!input.good())
    throw std::runtime_error(fmt::format("Unable to read file {}", filename.string()));
  std::stringstream buffer;
  buffer << input.rdbuf();

  CommandLineAndEnvironmentArguments args;
  args.virtuallinac_root = root;
  args.dry_run = false;
  args.input_json = buffer.str();
  SinkLogger log;
  json simulation_config = json::parse(args.input_json, nullptr, true, true);

  SUBCASE("CheckPlan") {
    uint64_t expected_num_sp = simulation_config.at("num_simulation_points").get<uint64_t>();
    uint64_t expected_num_particles = simulation_config.at("num_particles").get<uint64_t>();
    auto runTestSimulation = [=](const SimulationInputParams& in_params, const SimulationOutputParams&,
    LinacDetector*, LinacSource*, const Plan& plan, uint64_t num_histories, bool, Timing&) {
      CHECK(plan.getTotalNumParticles() == expected_num_particles);
      CHECK(num_histories == expected_num_particles);
      CHECK(plan.getNumSimulationPoints() == expected_num_sp);
      CHECK(in_params.biasing_mode == BiasingMode::DIRECTIONAL_SPLITTING);
      CHECK(in_params.track_particles == false);
    };
    CHECK(doLinacSimulation(args, runTestSimulation, &log) == 0);
  }
  SUBCASE("OverrideValues") {
    args.num_particles = 10000;
    uint64_t expected_num_sp = 200;
    int expected_seed = 236476;
    std::string additional_input = fmt::format("{{\"num_simulation_points\": {}, \"seed\": {} }}", expected_num_sp, expected_seed);
    args.additional_input_json = additional_input;
    auto runTestSimulation = [=](const SimulationInputParams& in_params, const SimulationOutputParams&,
    LinacDetector*, LinacSource*, const Plan& plan, uint64_t num_histories, bool, Timing&) {
      CHECK(plan.getTotalNumParticles() == args.num_particles);
      CHECK(num_histories == args.num_particles);
      CHECK(plan.getNumSimulationPoints() == expected_num_sp);
      CHECK(in_params.seed == expected_seed);
    };
    CHECK(doLinacSimulation(args, runTestSimulation, &log) == 0);
  }

}