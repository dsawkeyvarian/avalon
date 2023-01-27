#include "logger.h"
#include "external/popl.hpp"
#include <fmt/core.h>
#include "linac.h"
#include "simulation.h"


void printHelp(const popl::OptionParser& op) {
  std::cout << "CLI for VirtualLinac" << std::endl;
  std::cout << "Usage: linac.exe <simulation_config.json> <optional.json> <options>" << std::endl;
  std::cout << "It is possible to provide only <simulation_config.json>. Additional options can be given in <optional.json>." << std::endl;
  std::cout << "Additional options will be added to the configuration and duplicate values are replaced with optional ones." << std::endl;
  std::cout << "Use command line options to override/define default values or values in .json configuration files." << std::endl;
  std::cout << std::endl;
  std::cout << op << std::endl;
  std::cout << printValidCombinations() << std::endl;
}

int main(int argc, char** argv)
{
  try {
  popl::OptionParser op("Allowed options");
  auto help_option = op.add<popl::Switch>("h", "help", "Produce help message");
  auto visualize_option = op.add<popl::Switch>("", "visualize", "Visualize treatment head (no simulation is run)");
  auto starting_seed_opt = op.add<popl::Value<int>>("", "seed", "Starting seed for the simulations");
  auto num_threads_opt = op.add<popl::Value<int>>("", "num_threads", "Number of threads for the simulation");
  auto num_particles_opt = op.add<popl::Value<int>>("", "num_particles", "Number of particles for the simulation. Only if random_source is used");
  auto phsp_source_opt = op.add<popl::Value<std::string>>("", "phsp_source", "Enables phsp source and uses given file as source");
  auto output_directory_opt = op.add<popl::Value<std::string>>("", "output_dir", "Create directory and direct output files there");
  auto dry_run_option = op.add<popl::Switch>("", "dry_run", "Do dry run (no simulation is run)");
  auto root_option = op.add<popl::Value<std::string>>("", "root", "Set data root where gdml and stl fields are searched");
  op.parse(argc, argv);

  if (help_option->count() > 0) {
    printHelp(op);
    return 0;
  }
  if (op.non_option_args().size() == 0 || op.non_option_args().size() > 2) {
    std::cout << "Please provide correct input files. See help below." << std::endl << std::endl;
    printHelp(op);
    std::cerr << "Aborting" << std::endl;
    return -1;
  }
  COutAndFileLogger log;
  auto input_args = CommandLineAndEnvironmentArguments{};
  {
    const auto& filename = op.non_option_args()[0];
    std::ifstream input(filename);
    if(!input.good())
      throw std::runtime_error(fmt::format("Unable to read file {}", filename));
    std::stringstream buffer;
    buffer << input.rdbuf();
    input_args.input_json = buffer.str();
  }
  input_args.additional_input_json = std::nullopt;
  if(op.non_option_args().size() > 1) {
    const auto& filename = op.non_option_args()[1];
    std::ifstream input(filename);
    if(!input.good())
      throw std::runtime_error(fmt::format("Unable to read file {}", filename));
    std::stringstream buffer;
    buffer << input.rdbuf();
    input_args.additional_input_json = buffer.str();
  }

  input_args.seed = starting_seed_opt->is_set() ? std::make_optional(starting_seed_opt->value()) : std::nullopt;
  input_args.num_threads = num_threads_opt->is_set() ? std::make_optional(num_threads_opt->value()) : std::nullopt;
  input_args.num_particles = num_particles_opt->is_set() ? std::make_optional(num_particles_opt->value()) : std::nullopt;
  input_args.phsp_source_file = phsp_source_opt->is_set() ? std::make_optional(phsp_source_opt->value()) : std::nullopt;
  input_args.output_directory = output_directory_opt->is_set() ? std::make_optional(output_directory_opt->value()) : std::nullopt;
  input_args.dry_run = dry_run_option->is_set() ? true : false;
  input_args.visualize = visualize_option->is_set() ? true : false;

  const char* virtuallinac_root_env = getenv("VIRTUALLINAC_ROOT");
  if (!virtuallinac_root_env && !root_option->is_set()) {
    std::cerr << "Define VIRTUALLINAC_ROOT and point it to VirtualLinac repository or set --root argument. VirtualLinac needs the root path to find .gdml and .stl files." << std::endl;
    std::cerr << "Aborting" << std::endl;
    return -1;
  }
  if(root_option->is_set())
    input_args.virtuallinac_root = root_option->value();
  else
    input_args.virtuallinac_root = virtuallinac_root_env;
  
  return doLinacSimulation(input_args, &runG4Simulation, &log);
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Aborting" << std::endl;
    return -1;
  }
}