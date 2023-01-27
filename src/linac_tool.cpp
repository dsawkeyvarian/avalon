#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <array>
#include <bitset>
#include "external/popl.hpp"

#include "linac.h"

void printHelp() {
  std::cout << "CLI helper tool for VirtualLinac data" << std::endl;
  std::cout << "Usage: LinacTool.exe <command> <command related options>" << std::endl;
  std::cout << "To get more information about specific command, use -h" << std::endl;
  std::cout << "Available commands: " << std::endl;
  std::cout << "\tinspect" << std::endl;
}

void inspectPhsp(const std::string& filename);

int main(int argc, char** argv) {
  if (argc == 1) {
    printHelp();
    return -1;
  }
  if (std::string(argv[1]).compare("inspect") == 0) {
    popl::OptionParser op("Inspect options");
    auto help_option = op.add<popl::Switch>("h", "help", "Produce help message");
    auto phsp_filename = op.add<popl::Value<std::string>>("p", "phsp", "Inspect PhaseSpace file");
    op.parse(argc, argv);
    if (help_option->count() > 0) {
      std::cout << "VirtualLinac file inspector" << std::endl;
      std::cout << op << std::endl;
      return 0;
    }
    if (phsp_filename->is_set()) {
      try {
        inspectPhsp(phsp_filename->value());
      }
      catch (const std::runtime_error& err) {
        std::cerr << "Error encountered while reading phsp file" << std::endl;
        std::cerr << err.what() << std::endl;
      }
    }

  }
}


void inspectPhsp(const std::string& filename) {
  std::cout << "Reading PhaseSpace file " << filename << std::endl;
  const auto path = std::filesystem::path(filename);
  std::ifstream f_in(path, std::ios::binary);
  const size_t file_size = std::filesystem::file_size(path);
  PhaseSpaceDataHeader header;
  f_in.read(reinterpret_cast<char*>(&header), sizeof(header));
  if (header.magic_number != g_PHASESPACE_HEADER_MAGIC)
    throw std::runtime_error("Header magic number does not match");
  if (header.version != 1)
    throw std::runtime_error("Unable to process given version");
  std::cout << "Has tracking information: " << (header.has_tracking_information ? "true" : "false") << std::endl;
  std::cout << "Number of original histories: " << header.num_particle_histories << std::endl;

  std::cout << file_size << std::endl;
  if (header.data_size % header.data_stride != 0)
    throw std::runtime_error("Data size and stride are inconsistent");
  const size_t num_particles = header.data_size / header.data_stride;
  std::cout << "Number of particles: " << num_particles << std::endl;
  std::array<size_t, Particle::NUM_TYPES> num_particle_types = {};
  uint8_t* data = new uint8_t[file_size];
  f_in.seekg(header.data_offset);
  //f_in.seekg(0, f_in.beg);
  f_in.read((char*)data, header.data_size);
  double total_weight = 0.0;
  for (size_t i = 0; i < num_particles; ++i) {
    const PhaseSpaceParticle* phsp_particle = (PhaseSpaceParticle*)(&data[i * header.data_stride + header.particle_offset]);
    num_particle_types[phsp_particle->particle_type] += 1;
    total_weight += phsp_particle->weight;
  }
  std::cout << "Counts by particle type:" << std::endl;
  for (size_t i = 0; i < Particle::NUM_TYPES; ++i) {
    std::cout << Particle::getParticleName(static_cast<Particle::Type>(i)) << ": " << num_particle_types[i] << std::endl;
  }
  std::cout << std::endl;
  if (header.has_tracking_information) {
    const size_t num_traversed = 8;
    std::array<TraversedGeometry, num_traversed> traversed{
        TraversedGeometry::TARGET,
        TraversedGeometry::FLATTENING_FILTER,
        TraversedGeometry::IC,
        TraversedGeometry::PRIMARY_COLLIMATOR,
        TraversedGeometry::UPPER_JAW,
        TraversedGeometry::LOWER_JAW,
        TraversedGeometry::BASEPLATE,
        TraversedGeometry::MLC
    };
    std::array<size_t, num_traversed> traversed_counts = {};
    std::array<double, num_traversed> traversed_weights = {};
    for (size_t i = 0; i < num_particles; ++i) {
      for (size_t j = 0; j < num_traversed; ++j) {
        const uint32_t mask = 1 << static_cast<uint32_t>(traversed[j]);
        const auto phsp_particle = reinterpret_cast<PhaseSpaceParticle*>(&data[i * header.data_stride + header.particle_offset]);
        const auto track = reinterpret_cast<PhaseSpaceTrack*>(&data[i * header.data_stride + header.track_offset]);
        if (mask & track->traversed_geometry) {
          traversed_counts[j] += 1;
          traversed_weights[j] += phsp_particle->weight;
        }
      }
    }
    std::cout << "Particle tracking statistics" << std::endl;
    std::ios old_state(nullptr);
    old_state.copyfmt(std::cout);
    std::cout << std::fixed << std::setprecision(2);
    for (size_t i = 0; i < num_traversed; ++i) {
      size_t name_idx = static_cast<size_t>(traversed[i]);
      //const double percentage = 100.0 * (static_cast<double>(traversed_counts[i]) / static_cast<double>(num_particles));
      const double percentage = 100.0 * (static_cast<double>(traversed_weights[i]) / static_cast<double>(total_weight));
      std::cout << std::setw(25) << g_traversed_geometry_names[name_idx] << ": " << std::setw(6) << percentage << " %" << std::endl;
    }
    std::cout.copyfmt(old_state);
  }
  delete[] data;
  f_in.close();
}
