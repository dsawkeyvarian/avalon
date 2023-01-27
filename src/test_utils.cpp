#include "test_utils.h"
#include <fstream>
#include <fmt/core.h>

Dose::Dose(const fs::path& filepath)
{
    std::ifstream file(filepath, std::ios::binary);
    size_t file_size = fs::file_size(filepath);
    m_data.resize(file_size);
    file.read(reinterpret_cast<char*>(m_data.data()), file_size);
    if(!file.good())
        throw std::runtime_error(fmt::format("Unable to read file {} for Dose", filepath.string()));
    m_header = reinterpret_cast<DoseBinaryHeader*>(&m_data[0]);
    m_num_voxels = m_header->num_bins[0]*m_header->num_bins[1]*m_header->num_bins[2];
    m_doses = tcb::span(reinterpret_cast<double*>(&m_data[m_header->data_offset]), m_num_voxels);
}

std::array<size_t, 3> Dose::getVoxelCounts() const {
    return {m_header->num_bins[0], m_header->num_bins[1], m_header->num_bins[2]};
}

std::vector<double> Dose::getSlice1D(Axis axis, size_t i, size_t j) {
    size_t stride = 0;
    size_t offset = 0;
    auto counts = getVoxelCounts();
    if(axis == Axis::X) {
        stride = counts[1]*counts[2];
        offset = i*counts[0] + j;
    }
    else if(axis == Axis::Y) {
        stride = counts[2];
        offset = i*counts[1]*counts[2] + j;
    }
    else {
        stride = 1;
        offset = i*counts[1]*counts[2] + j*counts[0];
    }
    size_t num_points = counts[static_cast<size_t>(axis)];
    std::vector<double> slice(num_points);
    for(size_t k = 0; k < num_points; ++k) {
        slice[k] = m_doses[offset + k*stride];
    }
    return slice;
}

PhaseSpace::PhaseSpace(const fs::path& filepath)
{
    std::ifstream file(filepath, std::ios::binary);
    size_t file_size = fs::file_size(filepath);
    m_data.resize(file_size);
    file.read(reinterpret_cast<char*>(m_data.data()), file_size);
    if(!file.good())
        throw std::runtime_error(fmt::format("Unable to read file {} for PhaseSpace", filepath.string()));
    m_header = reinterpret_cast<PhaseSpaceDataHeader*>(&m_data[0]);
    if(m_header->magic_number != g_PHASESPACE_HEADER_MAGIC)
        throw std::runtime_error(fmt::format("File {} does not look like a Linac PhaseSpace file", filepath.string()));
    m_num_particles = m_header->data_size / m_header->data_stride;
}

const PhaseSpaceParticle* PhaseSpace::getParticle(size_t idx) const
{
    size_t data_idx = m_header->data_offset + idx*m_header->data_stride + m_header->particle_offset;
    return reinterpret_cast<const PhaseSpaceParticle*>(&m_data.at(data_idx));
}
const PhaseSpaceTrack* PhaseSpace::getParticleTrack(size_t idx) const
{
    if(!m_header->has_tracking_information)
        return nullptr;
    size_t data_idx = m_header->data_offset + idx*m_header->data_stride + m_header->track_offset;
    return reinterpret_cast<const PhaseSpaceTrack*>(&m_data.at(data_idx));
}