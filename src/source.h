#pragma once
#include "linac.h"
#include "external/span.hpp"
#include <atomic>

namespace CLHEP{
    class HepRotation;
}
class G4PrimaryVertex;
class G4ParticleDefinition;

class LinacSource {
public:
    LinacSource(bool rotate_with_collimator) : m_rotate_with_collimator(rotate_with_collimator), 
        m_rotation_matrix(nullptr), m_gantry_angle(0), m_collimator_angle(0) 
    {
        setAngles(m_gantry_angle, m_collimator_angle);
    }
    virtual ~LinacSource() {}
    virtual void setAngles(double gantry_angle, double collimator_angle);
    G4PrimaryVertex* generate();
protected:
    virtual G4PrimaryVertex* _generate() = 0;
    bool m_rotate_with_collimator;
    CLHEP::HepRotation* m_rotation_matrix;
    double m_gantry_angle;
    double m_collimator_angle;
};

class PhaseSpaceSource : public LinacSource
{
    std::unique_ptr<PhaseSpaceSourceData> m_phsp_data;
    std::atomic<uint64_t> m_source_data_index = 0;
public:
    PhaseSpaceSource(std::unique_ptr<PhaseSpaceSourceData> phsp_data)
        : LinacSource(false), m_phsp_data(std::move(phsp_data)) {};
protected:
    G4PrimaryVertex* _generate() override;
};

class GaussianElectronSource : public LinacSource
{
    GaussianRandomSource m_config;
public:
    GaussianElectronSource(const GaussianRandomSource& config);
protected:
    G4PrimaryVertex* _generate() override;
};

//Spectrum is given as cumulative probability distribution function
class GaussianElectronSourceWithSpectrum : public LinacSource
{
    GaussianRandomSource m_config;
    std::vector<double> m_energy_g;
    std::vector<double> m_cdf_g;
public:
    GaussianElectronSourceWithSpectrum(
        const GaussianRandomSource& config,
        tcb::span<const double> energy_g,  //energy grid
        tcb::span<const double> spectrum_cdf_g //cumulative distribution function of spectrum
    );
protected:
    G4PrimaryVertex* _generate() override;
};

class LinacUniformPointSource : public LinacSource
{
    double m_z;
    std::vector<double> m_energy_g;
    std::vector<double> m_cdf_g;
    Particle::Type m_particle_type;
    double m_x1;
    double m_x2;
    double m_y1;
    double m_y2;
public:
    LinacUniformPointSource(const UniformPointSourceInfo& config);
protected:
    G4PrimaryVertex* _generate() override;
};

std::unique_ptr<LinacSource> getTrueBeam6MeVSource();
std::unique_ptr<LinacSource> getHalcyonSource();
std::unique_ptr<LinacSource> getAvalon10XFFFSource();
std::unique_ptr<LinacSource> getAvalon10XSource();
std::unique_ptr<LinacSource> getAvalon15ESource();
