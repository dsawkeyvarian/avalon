#pragma once
#include <filesystem>
namespace fs = std::filesystem;
#include <array>
#include "external/span.hpp"
#include "linac.h"


class Dose {
public:
    Dose(const fs::path& filepath);
    size_t getTotalVoxelCount() const {return m_num_voxels;}
    std::array<size_t, 3> getVoxelCounts() const;
    double getPointDose(size_t x, size_t y, size_t z) {
        size_t idx = x*m_header->num_bins[1]*m_header->num_bins[2] + y*m_header->num_bins[2] + z;
        return m_doses[idx];
    }
    //axis == X, indices are (y,z)
    //axis == Y, indices are (x,z)
    //axis == Z, indices are (x,y)
    std::vector<double> getSlice1D(Axis axis, size_t i, size_t j);

    const DoseBinaryHeader* m_header;
    tcb::span<const double> m_doses;
private:
    size_t m_num_voxels;
    std::vector<uint8_t> m_data;
};

class PhaseSpace {
public:
    PhaseSpace(const fs::path& filepath);

    bool hasTrack() const {return m_header->has_tracking_information;}
    size_t getNumParticles() const {return m_num_particles;}
    const PhaseSpaceParticle* getParticle(size_t idx) const;
    const PhaseSpaceTrack* getParticleTrack(size_t idx) const;
private:
    const PhaseSpaceDataHeader* m_header;
    size_t m_num_particles;
    std::vector<uint8_t> m_data;
};