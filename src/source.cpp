#include "source.h"
#include "linac.h"
#include "simulation.h"
//Impl
#include "G4SystemOfUnits.hh"
#include "G4PrimaryVertex.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"
#include "G4RotationMatrix.hh"
#define _USE_MATH_DEFINES
#include <math.h>
#include "math_utils.h"

void LinacSource::setAngles(double gantry_angle, double collimator_angle)
{
    m_gantry_angle = gantry_angle;
    m_collimator_angle = collimator_angle;
    if(m_rotation_matrix) {
        delete m_rotation_matrix;
    }
    m_rotation_matrix = new G4RotationMatrix();
    if(this->m_rotate_with_collimator)
        m_rotation_matrix->rotateZ(m_collimator_angle);
    m_rotation_matrix->rotateY(m_gantry_angle);
}
G4PrimaryVertex* LinacSource::generate()
{
    auto vertex = this->_generate();
    {
        auto p = (*m_rotation_matrix)*vertex->GetPosition();
        vertex->SetPosition(p.x(), p.y(), p.z());
    }
    {
        auto primary = vertex->GetPrimary();
        auto m = (*m_rotation_matrix)*primary->GetMomentum();
        primary->SetMomentum(m.x(), m.y(), m.z());
    }
    return vertex;
}

G4PrimaryVertex* PhaseSpaceSource::_generate() {
    const size_t stride = m_phsp_data->stride;
    const uint8_t* buffer = m_phsp_data->data_buffer.data();
    uint64_t particle_idx = m_source_data_index.fetch_add(1);
    particle_idx = particle_idx % m_phsp_data->num_particles;
    const PhaseSpaceParticle& phsp_particle = *reinterpret_cast<const PhaseSpaceParticle*>(&buffer[m_phsp_data->particle_data_offset + particle_idx * stride]);
    G4ParticleTable* particle_table = G4ParticleTable::GetParticleTable();
    const int pdg_encoding = Particle::getPDGEncoding((Particle::Type)phsp_particle.particle_type);
    G4ParticleDefinition* particle_definition = particle_table->FindParticle(pdg_encoding);
    const double mass = particle_definition->GetPDGMass();
    const double energy = phsp_particle.energy * MeV + mass;
    const double momentum = sqrt(energy * energy - mass * mass);
    const double mom_sign = phsp_particle.weight > 0.0 ? 1.0 : -1.0;
    const double weight = abs(phsp_particle.weight);
    const double mom_z = mom_sign * sqrt(1.0 - phsp_particle.mom_x * phsp_particle.mom_x - phsp_particle.mom_y * phsp_particle.mom_y);
    const G4ThreeVector momentum_direction(momentum * phsp_particle.mom_x, momentum * phsp_particle.mom_y, momentum * mom_z);
    G4PrimaryVertex* vertex = new G4PrimaryVertex(G4ThreeVector(phsp_particle.x, phsp_particle.y, phsp_particle.z), 0.0);
    auto particle = new G4PrimaryParticle(particle_definition, momentum_direction.x(), momentum_direction.y(), momentum_direction.z());
    particle->SetWeight(weight);
    if (m_phsp_data->has_tracking_information) {
        const PhaseSpaceTrack& track_data = *reinterpret_cast<const PhaseSpaceTrack*>(&buffer[m_phsp_data->track_data_offset + particle_idx * stride]);
        particle->SetUserInformation(new LinacPrimaryParticleInformation(track_data));
    }
    vertex->SetPrimary(particle);
    return vertex;
}

inline G4PrimaryVertex* generateGaussianElectronWithEnergy(double energy, double z_pos, const GaussianRandomSource& config)
{
    const double position_z = z_pos;
    G4ThreeVector position(
		G4RandGauss::shoot(0., config.spot_size.x / mm) * mm + config.spot_position.x,
		G4RandGauss::shoot(0., config.spot_size.y / mm) * mm + config.spot_position.y,
		position_z);
    const double dir_x = G4RandGauss::shoot(0, config.angle_divergence.x) + config.angle.x;
    const double dir_y = G4RandGauss::shoot(0, config.angle_divergence.y) + config.angle.y;
    const double dir_z = -sqrt(1.0 - dir_x * dir_x - dir_y * dir_y);
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    auto particle_definition = particleTable->FindParticle("e-");
    const double mass = particle_definition->GetPDGMass();
	const double total_energy = energy * MeV + mass;
	const double momentum_length = sqrt(total_energy * total_energy - mass * mass);
    const auto momentum = momentum_length*G4ThreeVector(dir_x, dir_y, dir_z);
    G4PrimaryVertex* vertex = new G4PrimaryVertex(position, 0.0);
	auto particle = new G4PrimaryParticle(particle_definition, momentum.x(), momentum.y(), momentum.z());
    vertex->SetPrimary(particle);
    return vertex;
}

GaussianElectronSource::GaussianElectronSource(const GaussianRandomSource& config)
    : LinacSource(false), m_config(config)
{ } 

G4PrimaryVertex* GaussianElectronSource::_generate() {
	//TODO(jouko): why not use abs here?
	double energy = -1.0;
	while (energy < 0.0) {
		energy = G4RandGauss::shoot(m_config.energy / MeV, m_config.energy_sigma / MeV) * MeV;
	}
	return generateGaussianElectronWithEnergy(energy, 202.0 / 2.0 * cm, m_config);
}

GaussianElectronSourceWithSpectrum::GaussianElectronSourceWithSpectrum(
    const GaussianRandomSource& config,
    tcb::span<const double> energy_g,
    tcb::span<const double> spectrum_cdf_g)
    : LinacSource(false), m_config(config)
{
    m_energy_g.insert(m_energy_g.begin(), energy_g.begin(), energy_g.end());
    m_cdf_g.insert(m_cdf_g.begin(), spectrum_cdf_g.begin(), spectrum_cdf_g.end());
}

G4PrimaryVertex* GaussianElectronSourceWithSpectrum::_generate()
{
    //Inverse sample from energy CDF
    double value = G4RandFlat::shoot();
    double energy = interpolate(tcb::make_span(m_cdf_g), tcb::make_span(m_energy_g), value);
    return generateGaussianElectronWithEnergy(energy, 1002.0, m_config);
}


LinacUniformPointSource::LinacUniformPointSource(const UniformPointSourceInfo& config)
    : LinacSource(true)
{
    m_x1 = config.x1;
    m_x2 = config.x2;
    m_y1 = config.y1;
    m_y2 = config.y2;
    m_particle_type = config.particle;
    m_z = config.z;
    
    m_energy_g = config.spectrum_energies;
    m_cdf_g = calculateCDF(
        tcb::make_span(m_energy_g),
        tcb::make_span(config.spectrum_intensities));
}

G4PrimaryVertex* LinacUniformPointSource::_generate()
{
    //Rejection sample from lower sphere surface for a point inside the field
    const double sad = 1000.0;
    constexpr double sqrt2 = M_SQRT2;
    double max_extent = std::max(std::max(abs(m_x1), abs(m_x2)), std::max(abs(m_y1), abs(m_y2)));
    double max_radius = sqrt2*max_extent;
    const double max_z = -sad/sqrt(sad*sad + max_radius*max_radius);
    //Coordinates on unit sphere
    double x = 0;
    double y = 0;
    double z = -1.0;
    bool isInsideField = false;
    while(!isInsideField) {
        z = G4RandFlat::shoot(-1.0, max_z);
        double phi = G4RandFlat::shoot(-M_PI, M_PI);
        x = cos(phi)*sqrt(1 - z*z);
        y = sin(phi)*sqrt(1 - z*z);
        double r_iso = sad/abs(z);
        double x_iso = r_iso*x;
        double y_iso = r_iso*y;
        isInsideField = (m_x1 <= x_iso) && (x_iso <= m_x2) && (m_y1 <= y_iso) && (y_iso <= m_y2);
    }
    //Transfer particles to required height m_z
    double r_pos = (sad - m_z)/abs(z);
    double x_pos = r_pos*x;
    double y_pos = r_pos*y;
    double z_pos = m_z;
    
    //Inverse sample from energy CDF
    double value = G4RandFlat::shoot();
    double energy = interpolate(tcb::make_span(m_cdf_g), tcb::make_span(m_energy_g), value);
    //These can be only looked after G4 initialization
    G4ParticleTable* particle_table = G4ParticleTable::GetParticleTable();
    const int pdg_encoding = Particle::getPDGEncoding(m_particle_type);
    auto particle_definition = particle_table->FindParticle(pdg_encoding);
    const double mass = particle_definition->GetPDGMass();

	const double total_energy = energy * MeV + mass;
	const double momentum_length = sqrt(total_energy * total_energy - mass * mass);
    G4PrimaryVertex* vertex = new G4PrimaryVertex(x_pos, y_pos, z_pos, 0.0);
	auto particle = new G4PrimaryParticle(particle_definition, momentum_length*x, momentum_length*y, momentum_length*z);
    vertex->SetPrimary(particle);
    return vertex;
}

std::unique_ptr<LinacSource> getTrueBeam6MeVSource()
{
    const GaussianRandomSource truebeam_6MeV_src = {
        6.18, 0.053, //energy
        Vec2{0.0, 0.0}, //spot_position
        Vec2{0.6866, 0.7615}, //spot_size
        Vec2{0.0, 0.0}, //angle (radians)
        Vec2{0.001, 0.001} //angle_divergence (radians)
    };
    return std::make_unique<GaussianElectronSource>(truebeam_6MeV_src);
}
std::unique_ptr<LinacSource> getAvalon10XFFFSource()
{
    const GaussianRandomSource truebeam_10XFFF_src = {
        10.8, 0.0909, //energy
        Vec2{0.0, 0.0}, //spot_position
        Vec2{0.6778, 0.6778}, //spot_size
        Vec2{0.0, 0.0}, //angle (radians)
        Vec2{0.001, 0.001} //angle_divergence (radians)
    };
    return std::make_unique<GaussianElectronSource>(truebeam_10XFFF_src);
}
std::unique_ptr<LinacSource> getAvalon10XSource()
{
    const GaussianRandomSource truebeam_10X_src = {
        11.5, 0.0909, //energy
        Vec2{0.0, 0.0}, //spot_position
        Vec2{0.8345, 0.8710}, //spot_size
        Vec2{0.0, 0.0}, //angle (radians)
        Vec2{0.001, 0.001} //angle_divergence (radians)
    };
    return std::make_unique<GaussianElectronSource>(truebeam_10X_src);
}

std::unique_ptr<LinacSource> getHalcyonSource()
{
    //Halcyon electron histogram (left_edge, right_edge, intensity)
    std::vector<double> data = {
        0.00, 0.05, 3.015993e-02,
        0.05, 0.10, 2.407503e-02,
        0.10, 0.15, 2.169398e-02,
        0.15, 0.20, 3.121817e-02,
        0.20, 0.25, 3.756763e-02,
        0.25, 0.30, 3.545114e-02,
        0.30, 0.35, 3.148273e-02,
        0.35, 0.40, 4.418165e-02,
        0.40, 0.45, 2.989537e-02,
        0.45, 0.50, 3.439290e-02,
        0.50, 0.55, 2.857256e-02,
        0.55, 0.60, 2.724976e-02,
        0.60, 0.65, 2.698520e-02,
        0.65, 0.70, 2.592695e-02,
        0.70, 0.75, 1.957750e-02,
        0.75, 0.80, 2.222310e-02,
        0.80, 0.85, 2.777888e-02,
        0.85, 0.90, 2.169398e-02,
        0.90, 0.95, 1.957750e-02,
        0.95, 1.00, 2.063574e-02,
        1.00, 1.05, 2.645608e-02,
        1.05, 1.10, 2.407503e-02,
        1.10, 1.15, 2.777888e-02,
        1.15, 1.20, 2.486871e-02,
        1.20, 1.25, 3.148273e-02,
        1.25, 1.30, 2.486871e-02,
        1.30, 1.35, 2.989537e-02,
        1.35, 1.40, 2.539783e-02,
        1.40, 1.45, 2.857256e-02,
        1.45, 1.50, 2.619152e-02,
        1.50, 1.55, 2.698520e-02,
        1.55, 1.60, 2.883712e-02,
        1.60, 1.65, 3.201185e-02,
        1.65, 1.70, 3.862587e-02,
        1.70, 1.75, 3.227641e-02,
        1.75, 1.80, 3.545114e-02,
        1.80, 1.85, 3.756763e-02,
        1.85, 1.90, 3.095361e-02,
        1.90, 1.95, 3.227641e-02,
        1.95, 2.00, 4.074236e-02,
        2.00, 2.05, 4.365253e-02,
        2.05, 2.10, 5.053111e-02,
        2.10, 2.15, 4.709182e-02,
        2.15, 2.20, 5.211847e-02,
        2.20, 2.25, 5.185391e-02,
        2.25, 2.30, 5.132479e-02,
        2.30, 2.35, 4.920830e-02,
        2.35, 2.40, 4.788550e-02,
        2.40, 2.45, 5.502864e-02,
        2.45, 2.50, 5.079567e-02,
        2.50, 2.55, 4.285884e-02,
        2.55, 2.60, 4.973742e-02,
        2.60, 2.65, 4.285884e-02,
        2.65, 2.70, 4.788550e-02,
        2.70, 2.75, 4.603357e-02,
        2.75, 2.80, 5.158935e-02,
        2.80, 2.85, 4.603357e-02,
        2.85, 2.90, 3.650939e-02,
        2.90, 2.95, 3.677395e-02,
        2.95, 3.00, 4.391709e-02,
        3.00, 3.05, 4.338797e-02,
        3.05, 3.10, 3.809675e-02,
        3.10, 3.15, 4.021324e-02,
        3.15, 3.20, 4.418165e-02,
        3.20, 3.25, 4.232972e-02,
        3.25, 3.30, 4.550445e-02,
        3.30, 3.35, 4.232972e-02,
        3.35, 3.40, 5.132479e-02,
        3.40, 3.45, 4.920830e-02,
        3.45, 3.50, 5.688056e-02,
        3.50, 3.55, 6.058441e-02,
        3.55, 3.60, 6.508195e-02,
        3.60, 3.65, 6.534651e-02,
        3.65, 3.70, 7.513526e-02,
        3.70, 3.75, 8.413032e-02,
        3.75, 3.80, 6.640475e-02,
        3.80, 3.85, 7.566438e-02,
        3.85, 3.90, 7.037316e-02,
        3.90, 3.95, 7.830999e-02,
        3.95, 4.00, 7.830999e-02,
        4.00, 4.05, 9.074434e-02,
        4.05, 4.10, 8.254296e-02,
        4.10, 4.15, 7.487070e-02,
        4.15, 4.20, 8.360120e-02,
        4.20, 4.25, 7.619350e-02,
        4.25, 4.30, 8.386576e-02,
        4.30, 4.35, 7.857455e-02,
        4.35, 4.40, 8.651137e-02,
        4.40, 4.45, 9.206715e-02,
        4.45, 4.50, 9.259627e-02,
        4.50, 4.55, 8.862786e-02,
        4.55, 4.60, 7.539982e-02,
        4.60, 4.65, 8.756961e-02,
        4.65, 4.70, 1.103218e-01,
        4.70, 4.75, 1.246081e-01,
        4.75, 4.80, 1.296348e-01,
        4.80, 4.85, 1.388944e-01,
        4.85, 4.90, 1.325449e-01,
        4.90, 4.95, 1.293702e-01,
        4.95, 5.00, 1.298993e-01,
        5.00, 5.05, 1.476249e-01,
        5.05, 5.10, 1.582073e-01,
        5.10, 5.15, 1.608529e-01,
        5.15, 5.20, 1.568845e-01,
        5.20, 5.25, 1.658796e-01,
        5.25, 5.30, 1.772557e-01,
        5.30, 5.35, 1.851925e-01,
        5.35, 5.40, 2.296387e-01,
        5.40, 5.45, 2.394275e-01,
        5.45, 5.50, 2.587404e-01,
        5.50, 5.55, 2.658836e-01,
        5.55, 5.60, 3.172084e-01,
        5.60, 5.65, 3.214413e-01,
        5.65, 5.70, 2.910168e-01,
        5.70, 5.75, 3.047740e-01,
        5.75, 5.80, 3.283199e-01,
        5.80, 5.85, 3.777928e-01,
        5.85, 5.90, 4.370544e-01,
        5.90, 5.95, 5.000198e-01,
        5.95, 6.00, 7.299231e-01,
        6.00, 6.05, 2.369142e+00,
        6.05, 6.10, 1.420691e+00,
        6.10, 6.15, 1.330741e+00,
        6.15, 6.20, 1.625726e+00,
        6.20, 6.25, 1.013532e+00,
        6.25, 6.30, 1.216450e+00 };
      
    const GaussianRandomSource halcyon_src = {
        0.0, 0.0, //Gaussian energy is not used here
        Vec2{0.0, 0.0}, //spot_position
        Vec2{0.68517, 0.77031}, //spot_size
        Vec2{0.0, 0.0}, //angle (radians)
        Vec2{0.00153, 0.00199} //angle_divergence (radians)
    };
    //Compute CDF from histogram
    size_t num_bins = data.size() / 3;
    std::vector<double> energy_g(num_bins + 1, 0.0);
    std::vector<double> cdf_g(num_bins + 1, 0.0);
    energy_g[0] = 0.0;
    cdf_g[0] = 0.0;
    for(size_t i = 0; i < num_bins; ++i) {
        energy_g[i+1] = data[3*i + 1];
        double width = data[3*i + 1] - data[3*i + 0];
        cdf_g[i+1] = cdf_g[i] + data[3*i + 2]*width;
    }
    double I = cdf_g.back();
    for(size_t i = 0; i < cdf_g.size(); ++i)
        cdf_g[i] /= I;
    return std::make_unique<GaussianElectronSourceWithSpectrum>(
        halcyon_src, tcb::make_span(energy_g), tcb::make_span(cdf_g));
}