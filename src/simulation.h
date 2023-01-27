#pragma once
#include "linac.h"

#include "G4VUserTrackInformation.hh"
#include "G4VUserPrimaryParticleInformation.hh"

void runG4Simulation(
	const SimulationInputParams& input_params,
	const SimulationOutputParams& output_params,
	LinacDetector* detector,
	LinacSource* source,
	const Plan& plan,
	uint64_t num_original_histories,
	bool visualize, 
	Timing& timing); 

class LinacTrackInformation : public G4VUserTrackInformation
{
public:
	uint32_t traversed_geometry;
	PhysicsProcess creation_process;
	uint32_t creation_volume;
};

class LinacPrimaryParticleInformation : public G4VUserPrimaryParticleInformation
{
public:
	LinacPrimaryParticleInformation(const PhaseSpaceTrack& track_data) {
		m_track_info = LinacTrackInformation{};
		m_track_info.traversed_geometry = track_data.traversed_geometry;
		m_track_info.creation_process = static_cast<PhysicsProcess>(track_data.creation_process);
		m_track_info.creation_volume = static_cast<uint32_t>(track_data.creation_geometry);
	}
	void Print() const {
		G4cout << m_track_info.traversed_geometry << G4endl;
	}
	LinacTrackInformation m_track_info;
};