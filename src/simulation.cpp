#include "treatment_heads.h"
#include "source.h"
#include "plan.h"
#include "simulation.h"

#include "G4PhysListFactory.hh"
#include "G4SystemOfUnits.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
//#include "G4UIManager.hh"
#include "G4UImanager.hh"
#include "G4ScoringManager.hh"
#include "G4NistManager.hh"

//Visualization
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4SDManager.hh"
#include "G4Run.hh"
#include "G4UserRunAction.hh"
#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VUserPrimaryParticleInformation.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserParallelWorld.hh"
#include "G4VNestedParameterisation.hh"
#include "G4PVParameterised.hh"

#include "G4ParallelWorldPhysics.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4EllipticalTube.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSDoseDeposit3D.hh"

#include "G4UImessenger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"

#include <fmt/core.h>
#include <chrono>
#include <atomic>


G4ThreeVector Vec3ToG4(const Vec3& v) {
	return G4ThreeVector(v.x, v.y, v.z);
}

class LinacPhaseSpaceSphereScorer : public G4VPrimitiveScorer
{
public:
	LinacPhaseSpaceSphereScorer(G4String name) : G4VPrimitiveScorer(name, 0) {}

protected:
	G4bool ProcessHits(G4Step* step, G4TouchableHistory*)
	{
		G4Track* track = step->GetTrack();
		track->SetTrackStatus(G4TrackStatus::fStopAndKill);
		return true;
	}
};

//static G4Mutex g_phasespace_file_mutex = G4MUTEX_INITIALIZER; //Protects the phase space file from simultanous write access
static const size_t g_phsp_particle_write_buffer_size = 10000;

struct LinacPhaseSpaceWriterInfo {
	std::string phsp_filename;
	G4Mutex* file_mutex;
	bool m_track_particles;
	bool m_kill_particles;
};

class LinacPhaseSpaceWriter : public G4VPrimitiveScorer
{
public:

	LinacPhaseSpaceWriter(G4String name, LinacPhaseSpaceWriterInfo& info) : G4VPrimitiveScorer(name, 0)
		, m_buffer_capacity(g_phsp_particle_write_buffer_size)
		, m_info(info), m_buffer_idx(0), m_track_offset(0)
	{
		GetStrideAndOffset(m_info.m_track_particles, m_stride, m_particle_offset, m_track_offset);
		const size_t buffer_data_size = m_buffer_capacity * m_stride;
		m_buffer.resize(buffer_data_size);
	}

	//uint64_t m_num_particles[Particle::NUM_TYPES];
protected:
	const size_t m_buffer_capacity;
	LinacPhaseSpaceWriterInfo& m_info;
	size_t m_buffer_idx;
	size_t m_stride;
	size_t m_particle_offset;
	size_t m_track_offset;
	std::vector<uint8_t> m_buffer;


	G4bool ProcessHits(G4Step* step, G4TouchableHistory*)
	{
		//One could only use z momentum to check where the particle is going: edge case: The step happens *inside* of the phasespace volume
		auto post_step = step->GetPostStepPoint();
		//Should work if the capturing box is thin enough and it is in air (or relatively non-interacting material)
		bool is_exiting = post_step->GetStepStatus() == G4StepStatus::fGeomBoundary;
		if (!is_exiting) return false;
		G4StepPoint* step_point = step->GetPreStepPoint(); //Choose the point outside phase space volume
		const G4ThreeVector momentum_direction = step_point->GetMomentumDirection();
		G4Track* track = step->GetTrack();
		WriteToBuffer(step_point, track);
		if (m_info.m_track_particles)
			WriteTrackInfoToBuffer((LinacTrackInformation*)track->GetUserInformation());
		m_buffer_idx++;
		if (m_buffer_idx == m_buffer_capacity) {
			WriteBufferToFile();
			m_buffer_idx = 0;
		}
		if (m_info.m_kill_particles)
			track->SetTrackStatus(G4TrackStatus::fStopAndKill);

		return true;
	}


	void WriteToBuffer(G4StepPoint* step_point, G4Track* track)
	{
		const G4ThreeVector position = step_point->GetPosition();
		const G4ThreeVector momentum_direction = step_point->GetMomentumDirection();
		const float energy = static_cast<float>(step_point->GetKineticEnergy() / MeV);
		const float weight = (momentum_direction.z() > 0.0 ? 1.0f : -1.0f) * static_cast<float>(step_point->GetWeight());
		const uint32_t particle_type = static_cast<uint32_t>(Particle::getParticleType(track->GetDefinition()->GetPDGEncoding()));
		PhaseSpaceParticle& data = *reinterpret_cast<PhaseSpaceParticle*>(&m_buffer[m_buffer_idx*m_stride + m_particle_offset]);
		data.particle_type = particle_type;
		data.energy = energy;
		data.x = static_cast<float>(position.x());
		data.y = static_cast<float>(position.y());
		data.z = static_cast<float>(position.z());
		data.mom_x = static_cast<float>(momentum_direction.x());
		data.mom_y = static_cast<float>(momentum_direction.y());
		data.weight = weight;

    //if (fabs(position.x()) > 100.*mm) { // && fabs(position.y()) > 150.*mm) {
    //  //G4cout << "Keeping!" << G4endl;
    //  G4EventManager::GetEventManager()->KeepTheCurrentEvent();
    //}

	}

	void WriteTrackInfoToBuffer(LinacTrackInformation* info)
	{
		const uint32_t traversed = info->traversed_geometry;
		PhaseSpaceTrack* buffer = reinterpret_cast<PhaseSpaceTrack*>(&m_buffer[m_buffer_idx * m_stride + m_track_offset]);
		buffer->traversed_geometry = traversed;
		buffer->creation_geometry = static_cast<uint16_t>(info->creation_volume);
		buffer->creation_process = static_cast<uint16_t>(info->creation_process);
	}

public:
	void WriteBufferToFile()
	{
		G4AutoLock file_lock(m_info.file_mutex);
		std::ofstream f_out;
		f_out.open(m_info.phsp_filename, std::ios::app | std::ios::binary);
		if (!f_out) {
			G4cerr << "PhaseSpaceWriter: can't open file " << m_info.phsp_filename << G4endl;
			throw std::runtime_error("Error: can't open file " + m_info.phsp_filename);
		}
		f_out.write(reinterpret_cast<char*>(m_buffer.data()), m_buffer_idx * m_stride);
		f_out.close();
	}

	static size_t getDataOffset() {
		const size_t data_offset = 256;
		static_assert(sizeof(PhaseSpaceDataHeader) < data_offset);
		return data_offset;
	}

	static void GetStrideAndOffset(bool track_particles, size_t& stride, size_t& particle_offset, size_t& track_offset) {
		stride = sizeof(PhaseSpaceParticle);
		particle_offset = 0;
		track_offset = 0;
		if (track_particles) {
			track_offset = stride;
			stride += sizeof(PhaseSpaceTrack);
		}
	}

	static void InitFile(std::string filename) {
		std::ofstream f_out;
		f_out.open(filename,  std::ios::trunc | std::ios::binary);
		if(!f_out.good())
			throw std::runtime_error(fmt::format("Unable to open file {} for phasespace writing", filename));
		const size_t offset = getDataOffset();
		std::vector<uint8_t> data(offset, 0);
		f_out.write(reinterpret_cast<char*>(data.data()), offset);
		f_out.close();
	}

	static void WriteBinaryHeaderToFile(std::string filename, bool track_particles, uint64_t num_original_histories) 
	{
		std::fstream f_out(filename, std::ios::in | std::ios::out | std::ios::binary);
		const auto path = fs::path(filename);
		size_t file_size = fs::file_size(path);
		const size_t data_offset = getDataOffset();
		const size_t data_size = file_size - getDataOffset();
		PhaseSpaceDataHeader header{};
		header.num_particle_histories = num_original_histories;
		header.has_tracking_information = track_particles;
		header.magic_number = g_PHASESPACE_HEADER_MAGIC;
		header.version = 1;
		header.data_size = data_size;
		header.data_offset = data_offset;
		size_t stride, particle_offset, track_offset;
		GetStrideAndOffset(track_particles, stride, particle_offset, track_offset);
		header.data_stride = static_cast<uint32_t>(stride);
		header.particle_offset = static_cast<uint32_t>(particle_offset);
		header.track_offset = static_cast<uint32_t>(track_offset);
		f_out.seekp(0, std::ios::beg);
		f_out.write(reinterpret_cast<char*>(&header), sizeof(PhaseSpaceDataHeader));
		f_out.close();
	}
};

class LinacPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
	LinacSource* m_source;
public:
	LinacPrimaryGeneratorAction(LinacSource* source)
		: m_source(source) {}
		
	void GeneratePrimaries(G4Event* event) override {
		auto vertex = m_source->generate();
		event->AddPrimaryVertex(vertex);
	}
};

struct SensitiveDetectorNames {
	std::optional<std::string> monitor_chamber;
	std::vector<std::string> dose_grids;
	std::vector<std::string> phaseplanes;
};

class LinacRun : public G4Run
{
public:
	std::optional<size_t> m_monitor_chamber_idx;
	std::vector<size_t> m_dose_grid_idxs;
	std::vector<size_t> m_phase_plane_idxs;
	std::vector<std::string> m_sd_names;
	std::vector<G4THitsMap<G4double>*> m_hits;
	std::vector<G4int> m_collection_ids;

	static std::tuple<G4int, G4THitsMap<G4double>*> retrieveCollectionIDAndHitMap(const G4String& detector_name, bool debug = false)
	{
		G4SDManager* sd_manager = G4SDManager::GetSDMpointer();
		G4MultiFunctionalDetector* detector = dynamic_cast<G4MultiFunctionalDetector*>(sd_manager->FindSensitiveDetector(detector_name));
		if (detector == nullptr)
			throw std::runtime_error("Could not find detector with name " + detector_name);
		const G4int num_primitives = detector->GetNumberOfPrimitives();
		if (num_primitives != 1)
			throw std::runtime_error("Detector " + detector_name + " has wrong number of primitive scorers: " + std::to_string(num_primitives) + ". Expected: 1");
		G4VPrimitiveScorer* scorer = detector->GetPrimitive(0);
		const G4String collection_name = scorer->GetName();
		const G4String full_collection_name = detector_name + "/" + collection_name;
		const G4int collection_id = sd_manager->GetCollectionID(full_collection_name);
		if (debug) {
			G4cout << "Creating a hit map: " << G4endl << "\t" << full_collection_name << " " << num_primitives << G4endl;
			G4cout << "\t" << scorer << G4endl;
			G4cout << "\tCollection id: " << collection_id << G4endl;
		}
		return std::make_tuple(collection_id, new G4THitsMap<G4double>(detector_name, collection_name));
	}

	LinacRun(const SensitiveDetectorNames& names)
	{
		const bool debug = false;
		m_monitor_chamber_idx.reset();
		if(names.monitor_chamber.has_value()) {
			m_monitor_chamber_idx = m_sd_names.size();
			m_sd_names.push_back(names.monitor_chamber.value());
		}
		
		for (const auto& name : names.dose_grids) {
			m_dose_grid_idxs.push_back( m_sd_names.size());
			m_sd_names.push_back(name);
		}
		for (const auto& name : names.phaseplanes) {
			m_phase_plane_idxs.push_back(m_sd_names.size());
			m_sd_names.push_back(name);
		}
		try {
			for (const std::string& name : m_sd_names) {
				G4int collection_id = 0;
				G4THitsMap<G4double>* hits = NULL;
				std::tie(collection_id, hits) = retrieveCollectionIDAndHitMap(name, debug);
				m_hits.push_back(hits);
				m_collection_ids.push_back(collection_id);
			}
		}
		catch (const std::exception& e) {
			G4cerr << "Error while trying to retrieving collections and hits map in linac run" << G4endl;
			throw e;
		}
	}

	virtual void RecordEvent(const G4Event* event) {
		G4HCofThisEvent* pHCE = event->GetHCofThisEvent();
		if (!pHCE) return;
		const size_t num_detectors = m_collection_ids.size();
		for (size_t i = 0; i < num_detectors; ++i) {
			auto hit_map = reinterpret_cast<G4THitsMap<G4double>*>(pHCE->GetHC(m_collection_ids[i]));
			if(hit_map) //This test is required because it might be null pointer for phasespace planes
				*(m_hits[i]) += *hit_map;
		}
		G4Run::RecordEvent(event);
	}

	void Merge(const G4Run* run) {
		const LinacRun* local_run = static_cast<const LinacRun*>(run);
		const size_t num_detectors = m_collection_ids.size();
		for (size_t i = 0; i < num_detectors; ++i) {
			*(m_hits[i]) += *(local_run->m_hits[i]);
		}
		G4Run::Merge(run);
	}
};

void writeDoseGridBinary(const DoseGrid& config, uint64_t num_original_histories,
	std::vector<G4double>& dosegrid, double monitor_chamber_dose)
{
	const G4ThreeVector size = Vec3ToG4(config.size) / mm;
	const G4ThreeVector center_position = Vec3ToG4(config.center_position) / mm;
	const double original_histories = static_cast<double>(num_original_histories);
	DoseBinaryHeader header = DoseBinaryHeader{};
	const uint32_t version = 2;
	const size_t offset = 256;
	static_assert(offset >= sizeof(DoseBinaryHeader));
	header.version = version;
	header.data_offset = offset;
	header.num_bins[0] = config.num_voxels[0];
	header.num_bins[1] = config.num_voxels[1];
	header.num_bins[2] = config.num_voxels[2];
	header.size_in_mm[0] = size.x();
	header.size_in_mm[1] = size.y();
	header.size_in_mm[2] = size.z();
	header.center_position_in_mm[0] = center_position.x();
	header.center_position_in_mm[1] = center_position.y();
	header.center_position_in_mm[2] = center_position.z();
	if (config.rotation.has_value()) {
		const auto& rotation = config.rotation.value();
		header.is_rotated = true;
		header.rotation_axis[0] = rotation.axis.x;
		header.rotation_axis[1] = rotation.axis.y;
		header.rotation_axis[2] = rotation.axis.z;
		header.rotation_angle_in_rad = rotation.angle;
	}
	header.monitor_chamber_dose = monitor_chamber_dose / gray / original_histories;
	const size_t nx = (size_t)config.num_voxels[0];
	const size_t ny = (size_t)config.num_voxels[1];
	const size_t nz = (size_t)config.num_voxels[2];
	std::ofstream f_out;
	f_out.open(config.filename, std::ios::binary);
	f_out.write(reinterpret_cast<char*>(&header), sizeof(header));
	const uint8_t zeroes[offset] = {};
	f_out.write(reinterpret_cast<const char*>(zeroes), offset - sizeof(header)); //pad to offset
	for (size_t i = 0; i < nx * ny * nz; ++i) {
		const double dose = dosegrid[i];
		const double dose_per_original_history = dose / gray / original_histories;
		f_out.write(reinterpret_cast<const char*>(&dose_per_original_history), sizeof(dose_per_original_history));
	}
	f_out.close();
}

//Note that order of dose grids in SensitiveDetectorNames
//must match order dose_grid_configs for correct output
//Same for phsp and phsp_filenames
struct OutputInfo {
	std::vector<std::string> phsp_filenames;
	double monitor_chamber_dose;
	std::vector<std::vector<double>> dose_grids;
};

class LinacRunAction : public G4UserRunAction
{
	const SensitiveDetectorNames& m_sd_names;
	OutputInfo& m_output_info;

public:
	LinacRunAction(const SensitiveDetectorNames& sd_names, OutputInfo& output_info)
		: m_sd_names(sd_names), m_output_info(output_info)
	{}

	G4Run* GenerateRun()
	{
		return new LinacRun(m_sd_names);
	}

	void BeginOfRunAction(const G4Run*)
	{}

	void EndOfRunAction(const G4Run* run)
	{
		
		G4SDManager* sd_manager = G4SDManager::GetSDMpointer();
		const LinacRun* linac_run = dynamic_cast<const LinacRun*>(run);
		for (size_t i = 0; i < linac_run->m_phase_plane_idxs.size(); ++i) {
			const size_t phsp_idx = linac_run->m_phase_plane_idxs[i];
			auto phsp_detector = (G4MultiFunctionalDetector*)sd_manager->FindSensitiveDetector(linac_run->m_sd_names[phsp_idx]);
			auto phsp_writer = (LinacPhaseSpaceWriter*)phsp_detector->GetPrimitive(0);
			
			phsp_writer->WriteBufferToFile();
		}
		if (!IsMaster()) return;
		{	
			double monitor_dose = 0.0;
			if(linac_run->m_monitor_chamber_idx.has_value()) {
				auto dose_ptr = (*linac_run->m_hits[linac_run->m_monitor_chamber_idx.value()])[0];
				monitor_dose = dose_ptr ? *dose_ptr : 0.0;
			}
			m_output_info.monitor_chamber_dose += monitor_dose;
		}
		for (size_t i = 0; i < linac_run->m_dose_grid_idxs.size(); ++i) {
			const size_t dose_grid_idx = linac_run->m_dose_grid_idxs[i];
			G4THitsMap<G4double>& dosegrid_hit_map = *linac_run->m_hits[dose_grid_idx];
			auto& container = m_output_info.dose_grids.at(i);
			for (size_t j = 0; j < container.size(); ++j) {
				const double* dose_ptr = dosegrid_hit_map[(int)j];
				const double dose = dose_ptr ? *dose_ptr : 0.0;
				container[j] += dose;	
 			}
		}
	}
};

class LinacTrackingAction : public G4UserTrackingAction
{
public:
	LinacTrackingAction(const std::unordered_map<std::string, PhysicsProcess> & name_to_process, const std::unordered_map<int, TraversedGeometry>& id_to_traversed)
		: m_name_to_process(name_to_process), m_id_to_traversed(id_to_traversed) {}

	const std::unordered_map<int, TraversedGeometry>& m_id_to_traversed;
	const std::unordered_map<std::string, PhysicsProcess>& m_name_to_process;

	virtual void PreUserTrackingAction(const G4Track* track)
	{
		if (track->GetTrackID() != 1) {
			return;
		}
		const auto user_info = track->GetDynamicParticle()->GetPrimaryParticle()->GetUserInformation();
		LinacTrackInformation* track_info = nullptr;
		if (user_info) {
			const LinacPrimaryParticleInformation* info = reinterpret_cast<const LinacPrimaryParticleInformation*>(user_info);
			track_info = new LinacTrackInformation(info->m_track_info);
		}
		else {
			track_info = new LinacTrackInformation;
			*track_info = {};
			track_info->creation_process = PhysicsProcess::PRIMARY;
			track_info->creation_volume = static_cast<uint32_t>(TraversedGeometry::SOURCE);
		}
		track->SetUserInformation(track_info);
	}

	virtual void PostUserTrackingAction(const G4Track* track)
	{
		G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
		if (!secondaries) return;
		auto track_info = (LinacTrackInformation*)track->GetUserInformation();
		const size_t num = secondaries->size();
		for (size_t i = 0; i < num; ++i) {
			auto new_track_info = new LinacTrackInformation;
			*new_track_info = *track_info;
			const auto& secondary = (*secondaries)[i];
			secondary->SetUserInformation(new_track_info);
			const auto& process_name = secondary->GetCreatorProcess()->GetProcessName();
			const auto& it = m_name_to_process.find(process_name);
			if (it == m_name_to_process.end()) {
				new_track_info->creation_process = PhysicsProcess::UNKNOWN;
			}
			else {
				new_track_info->creation_process = it->second;
			}
			G4LogicalVolume* volume = secondary->GetTouchable()->GetVolume()->GetLogicalVolume();
			const TraversedGeometry traversed = m_id_to_traversed.find(volume->GetInstanceID())->second;
			uint32_t volume_id = traversed == TraversedGeometry::NONE ? static_cast<uint32_t>(TraversedGeometry::OTHER) : static_cast<uint32_t>(traversed);
			new_track_info->creation_volume = volume_id;
		}
	}
};

class LinacSteppingAction : public G4UserSteppingAction
{
public:
	LinacSteppingAction(const std::unordered_map<int, TraversedGeometry>& id_to_traversed)
		: m_id_to_traversed(id_to_traversed) {}
	virtual ~LinacSteppingAction() {}

	virtual void UserSteppingAction(const G4Step* step) {
		G4StepPoint* pre = step->GetPreStepPoint();
		G4LogicalVolume* volume = pre->GetTouchable()->GetVolume()->GetLogicalVolume();
		const TraversedGeometry traversed = m_id_to_traversed.find(volume->GetInstanceID())->second;
		auto track_info = (LinacTrackInformation*)step->GetTrack()->GetUserInformation();
		if (traversed != TraversedGeometry::NONE)
			track_info->traversed_geometry |= 1 << (uint32_t)traversed;
	}

	const std::unordered_map<int, TraversedGeometry>& m_id_to_traversed;
};

struct LinacInitializationInfo {
	const SensitiveDetectorNames& sd_names;
	LinacSource* source;
	//const std::variant<GaussianRandomSource, PhaseSpaceSourceData> source;
	OutputInfo* output_info;
	bool particle_tracking;
	const std::unordered_map<int, TraversedGeometry>* id_to_traversed;
	const std::unordered_map<std::string, PhysicsProcess>* name_to_physics_process;
};

class LinacActionInitialization : public G4VUserActionInitialization
{
	const LinacInitializationInfo& m_info;

public:
	LinacActionInitialization(const LinacInitializationInfo& info)
		: m_info(info) {}
	virtual ~LinacActionInitialization() {}

	virtual void Build() const {
		SetUserAction(new LinacPrimaryGeneratorAction(m_info.source));
		SetUserAction(new LinacRunAction(m_info.sd_names, *m_info.output_info));
		if (m_info.particle_tracking) {
			SetUserAction(new LinacTrackingAction(*m_info.name_to_physics_process, *m_info.id_to_traversed));
			SetUserAction(new LinacSteppingAction(*m_info.id_to_traversed));
		}
	}

	virtual void BuildForMaster() const {
		SetUserAction(new LinacRunAction(m_info.sd_names, *m_info.output_info));
	}
};

class LinacDoseGrid : public G4VUserParallelWorld
{
public:
	const DoseGrid& m_config;
	G4VPhysicalVolume* m_scoring_grid_physical;

	LinacDoseGrid(const G4String& parallelWorldName, const DoseGrid& dose_config)
		: G4VUserParallelWorld(parallelWorldName), m_config(dose_config), m_scoring_grid_physical(nullptr)
	{}

	void Construct()
	{

		G4VPhysicalVolume* world = GetWorld();
		G4LogicalVolume* logical_world = world->GetLogicalVolume();

		const int num_x = m_config.num_voxels[0];
		const int num_y = m_config.num_voxels[1];
		const int num_z = m_config.num_voxels[2];
		const G4ThreeVector box_size = Vec3ToG4(m_config.size);
		const G4ThreeVector pos = Vec3ToG4(m_config.center_position);
		//Scoring grid is divided into uniform cells
		auto cell_size = G4ThreeVector(box_size.x() / (double)num_x, box_size.y() / (double)num_y, box_size.z() / (double)num_z);
		//Create volume for replicas (replicas must exactly fit into the mother volume which is why we dont use mother_volume directly)
		G4cout << "Constructing DoseGrid with name " << fWorldName << G4endl;
		const G4String replica_mother_name = "DoseGridReplica";
		G4cout << num_x << " " << num_y << " " << num_z << G4endl,
			G4cout << "box_size: " << box_size << G4endl;
		G4cout << "box_pos: " << pos << G4endl;
		auto replica_box = new G4Box(replica_mother_name, box_size.x() / 2.0, box_size.y() / 2.0, box_size.z() / 2.0);
		auto logical_replica_box = new G4LogicalVolume(replica_box, nullptr, replica_mother_name);
		if (m_config.rotation.has_value()) {
			const auto& config_rot = m_config.rotation.value();
			G4RotationMatrix rotate = G4Rotate3D(config_rot.angle, Vec3ToG4(config_rot.axis)).getRotation();
			std::cout << "Rotation test: " << rotate * G4ThreeVector(1.0, 0.0, 1.0) << std::endl;
			new G4PVPlacement(new G4RotationMatrix(rotate), pos, logical_replica_box, replica_mother_name, logical_world, false, 0);
		}
		else {
			new G4PVPlacement(nullptr, pos, logical_replica_box, replica_mother_name, logical_world, false, 0);
		}
		const G4String replica_name_x = "ReplicaX";
		auto box_replica_x = new G4Box(replica_name_x, cell_size.x() / 2.0, box_size.y() / 2.0, box_size.z() / 2.0);
		auto logical_replica_x = new G4LogicalVolume(box_replica_x, nullptr, replica_name_x);
		new G4PVReplica(replica_name_x, logical_replica_x, logical_replica_box, kXAxis, num_x, cell_size.x());

		const G4String replica_name_y = "ReplicaY";
		auto box_replica_y = new G4Box(replica_name_y, cell_size.x() / 2.0, cell_size.y() / 2.0, box_size.z() / 2.0);
		auto logical_replica_y = new G4LogicalVolume(box_replica_y, nullptr, replica_name_y);
		new G4PVReplica(replica_name_y, logical_replica_y, logical_replica_x, kYAxis, num_y, cell_size.y());

		const G4String replica_name_z = "ReplicaZ";
		auto box_replica_z = new G4Box(replica_name_z, cell_size.x() / 2.0, cell_size.y() / 2.0, cell_size.z() / 2.0);
		auto logical_scoring = new G4LogicalVolume(box_replica_z, nullptr, replica_name_z);
		auto physical_scoring = new G4PVReplica(replica_name_z, logical_scoring, logical_replica_y, kZAxis, num_z, cell_size.z());

		G4cout << "Dose scoring grid size " << box_size / mm << " mm" << G4endl;
		G4cout << "Segmentation (" << num_x << ", " << num_y << ", " << num_z << ")" << G4endl;
		G4cout << "Cell size " << cell_size / mm << " mm" << G4endl;
		m_scoring_grid_physical = physical_scoring;
	}


	void ConstructSD()
	{
		//G4MultiFunctionalDetector* sd = new G4MultiFunctionalDetector("SD_" + fWorldName);
		G4MultiFunctionalDetector* sd = new G4MultiFunctionalDetector(fWorldName); //NOTE: This line has changed, check if not working!
		SetSensitiveDetector(m_scoring_grid_physical->GetLogicalVolume(), sd);
		G4VPrimitiveScorer* total_dose = new G4PSDoseDeposit3D(fWorldName, m_config.num_voxels[0], m_config.num_voxels[1], m_config.num_voxels[2]);
		sd->RegisterPrimitive(total_dose);
	}
};

class LinacWaterPhantom : public G4VUserParallelWorld
{
public:
	const WaterBoxPhantom m_config;

	LinacWaterPhantom(const G4String& parallel_world_name, const WaterBoxPhantom& phantom_config)
		: G4VUserParallelWorld(parallel_world_name), m_config(phantom_config)
	{}

	void Construct()
	{	
		//Visualization
		const auto waterphantom_vis = new G4VisAttributes(G4Colour(0.3, 0.8, 0.3, 0.9));
		waterphantom_vis->SetForceSolid(true);
		waterphantom_vis->SetVisibility(true);
		G4VPhysicalVolume* world = GetWorld();
		G4LogicalVolume* logical_world = world->GetLogicalVolume();
		const auto center_position = Vec3ToG4(m_config.center_position);
		const auto box_size = Vec3ToG4(m_config.size);
		G4cout << "Constructing WaterPhantom" << G4endl;
		G4cout << "\tPosition: " << center_position / mm << G4endl;
		G4cout << "\tSize: " << box_size / mm << G4endl;
		const G4String name = "WaterPhantom";
		G4NistManager* NISTman = G4NistManager::Instance();
		auto box = new G4Box(name, box_size.x() / 2.0, box_size.y() / 2.0, box_size.z() / 2.0);
		auto water = NISTman->FindOrBuildMaterial("G4_WATER");
		//auto logical = new G4LogicalVolume(box, G4Material::GetMaterial("G4_Pb"), name, 0, 0, 0);
		auto logical = new G4LogicalVolume(box, water, name, 0, 0, 0);
		logical->SetVisAttributes(waterphantom_vis);
		new G4PVPlacement(nullptr, center_position, logical, name, logical_world, false, 0);
	}
};

class LinacWaterCylinderPhantom : public G4VUserParallelWorld
{
public:
	const WaterCylinderPhantom m_config;

	LinacWaterCylinderPhantom(const G4String& parallel_world_name, const WaterCylinderPhantom& phantom_config)
		: G4VUserParallelWorld(parallel_world_name), m_config(phantom_config)
	{}

	void Construct()
	{	
		//Visualization
		const auto waterphantom_vis = new G4VisAttributes(G4Colour(0.3, 0.8, 0.3, 0.9));
		waterphantom_vis->SetForceSolid(true);
		waterphantom_vis->SetVisibility(true);
		G4VPhysicalVolume* world = GetWorld();
		G4LogicalVolume* logical_world = world->GetLogicalVolume();
		const auto center_position = Vec3ToG4(m_config.center_position);
		//const auto box_size = Vec3ToG4(m_config.size);
		G4cout << "Constructing WaterPhantom" << G4endl;
		G4cout << "\tPosition: " << center_position / mm << G4endl;
		G4cout << "\tRadius: " << m_config.radius << G4endl;
		G4cout << "\tHeight" << m_config.height / mm << G4endl;
		const G4String name = "WaterPhantom";
		G4NistManager* NISTman = G4NistManager::Instance();
		auto cylinder = new G4EllipticalTube(name, m_config.radius, m_config.radius, m_config.height/2.0);
		auto rotation = new G4RotationMatrix(G4RotationMatrix::IDENTITY);
		const double PI = 3.14159265358979323846;
		if(m_config.height_axis == Axis::Y)
			rotation->rotateX(PI/2.0);
		else if(m_config.height_axis == Axis::X)
			rotation->rotateY(PI/2.0);
		auto water = NISTman->FindOrBuildMaterial("G4_WATER");
		auto logical = new G4LogicalVolume(cylinder, water, name, 0, 0, 0);
		logical->SetVisAttributes(waterphantom_vis);
		new G4PVPlacement(rotation, center_position, logical, name, logical_world, false, 0);
	}
};

class VoxelPhantomNestedParametrisation : public G4VNestedParameterisation
{
	size_t m_num_cells_x;
	size_t m_num_cells_y;
	size_t m_num_cells_z;
	G4ThreeVector m_voxel_size;
	std::vector<G4Material*> m_material_in_cells;
	std::vector<G4Material*> m_unique_materials;
	std::vector<double> m_z_pos;
public:
	VoxelPhantomNestedParametrisation(
		const G4ThreeVector& voxel_size,
		const std::vector<G4Material*>& unique_materials,
		size_t num_cells_x,
		size_t num_cells_y, 
		size_t num_cells_z,
		const std::vector<G4Material*>& materials) 
		: G4VNestedParameterisation(), m_voxel_size(voxel_size), 
		m_material_in_cells(materials), m_num_cells_x(num_cells_x),
		m_num_cells_y(num_cells_y), m_num_cells_z(num_cells_z),
		m_unique_materials(unique_materials) 
	{	
		std::cout << "Z positions: " << m_num_cells_z << std::endl;
		std::cout << "Voxel size z: " << m_voxel_size.z() << std::endl;
		for(size_t i = 0; i < m_num_cells_z; ++i) {
			double z = (-static_cast<int>(m_num_cells_z) + 1 + 2*static_cast<int>(i))*(0.5*voxel_size.z());
			std::cout << z << std::endl;
			m_z_pos.emplace_back(z);
		}
	}
		

	G4Material* ComputeMaterial(G4VPhysicalVolume* currentVol,
        const G4int repNo, const G4VTouchable* parentTouch = nullptr) override 
	{
		//if(parentTouch == 0) return m_materials[0];
		G4int ix = parentTouch->GetReplicaNumber(1);
  		G4int iy = parentTouch->GetReplicaNumber(0);
  		G4int iz = repNo;
		//size_t idx = ix*m_num_cells_y*m_num_cells_z + iy*m_num_cells_z + iz;
		//auto material = m_material_in_cells.at(idx);
		//std::cout << fmt::format("({}, {}, {}): {}, {}\n", ix, iy, iz, idx, (void*)material);
		return m_material_in_cells.at(ix*m_num_cells_y*m_num_cells_z + iy*m_num_cells_z + iz);
	}
	G4int GetNumberOfMaterials() const override
	{
		return static_cast<int>(m_unique_materials.size());
	}
	G4Material* GetMaterial(G4int idx) const override
	{
		return m_unique_materials[idx];
	}
	virtual void ComputeTransformation(const G4int no, G4VPhysicalVolume* currentPV) const override
	{
		//std::cout << "Compute transformation: " << (uint64_t)currentPV << std::endl;
		G4ThreeVector position(0.0, 0.0, m_z_pos[no]);
		currentPV->SetTranslation(position);
		//std::cout << "Done" << std::endl;
	}
	void ComputeDimensions(G4Box& box, const G4int, const G4VPhysicalVolume* ) const override
	{
		box.SetXHalfLength(m_voxel_size.x()/2.0);
		box.SetYHalfLength(m_voxel_size.y()/2.0);
		box.SetZHalfLength(m_voxel_size.z()/2.0);
	}
};

class LinacVoxelPhantom : public G4VUserParallelWorld
{
public:
	const VoxelBoxPhantom m_config;
	std::unordered_map<size_t, G4Material*> m_materials;

	LinacVoxelPhantom(const G4String& parallel_world_name, const VoxelBoxPhantom& phantom_config)
		: G4VUserParallelWorld(parallel_world_name), m_config(phantom_config)
	{}

	static std::unordered_map<size_t, G4Material*> buildUniqueMaterialCompositions(
		//const std::vector<VoxelBoxPhantom::Material>& materials,
		const std::unordered_map<size_t, VoxelBoxPhantom::Material>& materials) 
	{
		std::unordered_map<size_t, G4Material*> material_compositions;
		G4NistManager* NISTman = G4NistManager::Instance();
		const size_t num_mat = materials.size();
		for(auto& it : materials) {
		//for(size_t i = 0; i < num_mat; ++i) {
			size_t i = it.first;
			if(std::holds_alternative<std::string>(it.second)) {
				const auto& mat_name = std::get<std::string>(it.second);
				auto material = NISTman->FindOrBuildMaterial(mat_name);
				if(!material)
					throw std::runtime_error(fmt::format("Unable to find Geant4 material '{}'", mat_name));
				material_compositions[i] = material;
			}
			else if(std::holds_alternative<VoxelBoxPhantom::UserMaterial>(it.second)) {
				std::cout << "Creating user material!" << std::endl;
				const auto& mat = std::get<VoxelBoxPhantom::UserMaterial>(it.second);
				const size_t num_components = mat.names.size();
				G4Material* g4_mat = new G4Material(fmt::format("material_{}", i), 
					mat.default_density, static_cast<int>(num_components));
				std::cout << fmt::format("density, components: {}, {}\n", mat.default_density, num_components);
				for(size_t j = 0; j < num_components; ++j) {
					auto element = NISTman->FindOrBuildElement(mat.names[j]);
					std::cout << fmt::format("Adding element {}\n", mat.names[j]);
					if(!element)
						throw std::runtime_error(fmt::format("Unable to find/build element {} in material index {}", mat.names[j], j));
					g4_mat->AddElement(element, mat.mass_fractions[j]);
				}
				material_compositions[i] = g4_mat;
			}
			else {
				throw std::runtime_error("Unknown variant in VoxelBoxPhantom materials");
			}
		}
		return material_compositions;
	}

	void Construct()
	{	
		//Visualization
		const auto waterphantom_vis = new G4VisAttributes(G4Colour(0.3, 0.8, 0.3, 0.9));
		waterphantom_vis->SetForceSolid(true);
		waterphantom_vis->SetVisibility(true);
		G4VPhysicalVolume* world = GetWorld();
		G4LogicalVolume* logical_world = world->GetLogicalVolume();
		const auto center_position = Vec3ToG4(m_config.center_position);
		const auto box_size = G4ThreeVector(
			m_config.voxel_sizes.x*m_config.voxel_counts[0],
			m_config.voxel_sizes.y*m_config.voxel_counts[1],
			m_config.voxel_sizes.z*m_config.voxel_counts[2]);
		G4cout << "Constructing VoxelPhantom" << G4endl;
		G4cout << "\tPosition: " << center_position / mm << G4endl;
		G4cout << "\tSize: " << box_size / mm << G4endl;
		const G4String name = "VoxelPhantom";
		G4NistManager* NISTman = G4NistManager::Instance();
		m_materials = buildUniqueMaterialCompositions(m_config.idx_to_material);
		auto box = new G4Box(name, box_size.x() / 2.0, box_size.y() / 2.0, box_size.z() / 2.0);
		//auto logical = new G4LogicalVolume(box, G4Material::GetMaterial("G4_Pb"), name, 0, 0, 0);
		auto water = NISTman->FindOrBuildMaterial("G4_WATER");
		auto logical = new G4LogicalVolume(box, water, name, 0, 0, 0);
		logical->SetVisAttributes(waterphantom_vis);
		new G4PVPlacement(nullptr, center_position, logical, name, logical_world, false, 0);
		size_t nxCells = m_config.voxel_counts[0];
  		size_t nyCells = m_config.voxel_counts[1];
  		size_t nzCells = m_config.voxel_counts[2];
		double voxel_size_x = box_size.x()/static_cast<double>(nxCells);
		double voxel_size_y = box_size.y()/static_cast<double>(nyCells);
		double voxel_size_z = box_size.z()/static_cast<double>(nzCells);
	
		G4String xRepName("RepX");
		G4VSolid* solXRep =
    		new G4Box(xRepName, voxel_size_x/2.0, box_size.y()/2., box_size.z()/2.);
		G4LogicalVolume* logXRep = new G4LogicalVolume(solXRep, water, xRepName);
		new G4PVReplica(xRepName, logXRep, logical, kXAxis, static_cast<int>(nxCells), voxel_size_x);
		G4String yRepName("RepY");
		G4VSolid* solYRep = new G4Box(yRepName, voxel_size_x/2., voxel_size_y/2., box_size.z()/2.);
		G4LogicalVolume* logYRep = new G4LogicalVolume(solYRep, water, yRepName);
		new G4PVReplica(yRepName, logYRep, logXRep, kYAxis, static_cast<int>(nyCells), voxel_size_y);
		G4String zVoxName("VoxZ");
		G4VSolid* solVoxel = new G4Box(zVoxName, voxel_size_x/2., voxel_size_y/2., voxel_size_z/2.);
		auto fLVPhantomSens = new G4LogicalVolume(solVoxel, water, zVoxName);
		std::vector<G4Material*> phantom_materials(nzCells*nyCells*nxCells, water);
		for(size_t i = 0; i < phantom_materials.size(); ++i) {
			//std::cout << fmt::format("{}: {}\n", i, m_config.material_idxs[i]);
			phantom_materials[i] = m_materials[m_config.material_idxs[i]];
		}
		std::vector<G4Material*> unique_materials;
		for(auto it : m_materials) {
			unique_materials.emplace_back(it.second);
		}
		//G4Material* lead = NISTman->FindOrBuildMaterial("G4_Pb");
		//G4Material* air = NISTman->FindOrBuildMaterial("G4_AIR");
		//phantom_materials[1] = air;
		std::cout << "Start constructing VoxelPhantomNestedParametrisation" << std::endl;
		auto param_phantom = new VoxelPhantomNestedParametrisation(
			G4ThreeVector(voxel_size_x, voxel_size_y, voxel_size_z),
			unique_materials,
			nxCells, nyCells, nzCells, 
			phantom_materials);
		std::cout << "Start constructing G4PVParameterised" << std::endl;
		new G4PVParameterised("PhantomSens",     // their name
            fLVPhantomSens,    // their logical volume
            logYRep,           // Mother logical volume
            kUndefined,        // Are placed along this axis 
            static_cast<int>(nzCells),           // Number of cells
            param_phantom);     // Parameterisation.
		std::cout << "Finished VoxelPhantom Construction" << std::endl;
	}

	void ConstructSDandField() {}
};

class LinacPhaseSpacePlane : public G4VUserParallelWorld
{
public:
	const PhaseSpacePlane& m_config;
	G4LogicalVolume* m_phase_space_logical;
	LinacPhaseSpaceWriterInfo& m_info;

	LinacPhaseSpacePlane(const G4String& parallel_world_name, const PhaseSpacePlane& config, LinacPhaseSpaceWriterInfo& info)
		: G4VUserParallelWorld(parallel_world_name), m_config(config), m_info(info) {}

	void Construct() {
		G4VPhysicalVolume* world = GetWorld();
		G4LogicalVolume* logical_world = world->GetLogicalVolume();

		G4LogicalVolume* gantry = nullptr;
		G4LogicalVolume* collimator = nullptr;
		//findGantryAndCollimator(logical_world, &gantry, &collimator);
		//std::cout << "PhaseSpacePlane: " << std::endl << "Gantry found: " << (gantry ? "true" : "false") << std::endl << "Collimator found: " << (collimator ? "true" : "false") << std::endl;

		const G4ThreeVector size = Vec3ToG4(m_config.size);
		const G4ThreeVector pos = Vec3ToG4(m_config.center_position);
		auto phase_space_box = new G4Box("PhaseSpaceBox", size.x() / 2.0, size.y() / 2.0, size.z() / 2.0);
		m_phase_space_logical = new G4LogicalVolume(phase_space_box, nullptr, "PhaseSpacePlaneLogical");
		new G4PVPlacement(nullptr, G4ThreeVector(pos.x(), pos.y(), pos.z()),
			m_phase_space_logical, "PhaseSpacePlane", logical_world, false, 0);
	}

	void ConstructSD() {
		G4MultiFunctionalDetector* sd = new G4MultiFunctionalDetector(fWorldName);
		SetSensitiveDetector(m_phase_space_logical, sd);
		LinacPhaseSpaceWriter* phsp_writer = new LinacPhaseSpaceWriter(fWorldName, m_info);
		sd->RegisterPrimitive(phsp_writer);
	}
};

/*
class LinacPhaseSpaceSphere : public G4VUserParallelWorld
{
public:
	G4LogicalVolume* m_phase_space_sphere_logical;

	LinacPhaseSpaceSphere(const G4String& parallel_world_name) : G4VUserParallelWorld(parallel_world_name) {}

	void Construct() {
		G4VPhysicalVolume* world = GetWorld();
		G4LogicalVolume* logical_world = world->GetLogicalVolume();

		G4LogicalVolume* gantry = nullptr;
		G4LogicalVolume* collimator = nullptr;
		findGantryAndCollimator(logical_world, &gantry, &collimator);
		std::cout << "PhaseSpaceSphere: " << std::endl << "Gantry found: " << (gantry ? "true" : "false") << std::endl << "Collimator found: " << (collimator ? "true" : "false") << std::endl;

		G4Sphere* phase_space_sphere = new G4Sphere("phasespace_sphere", 9.*cm, 10.*cm, 0. * deg, 180. * deg, 0. * deg, 180. * deg);
		m_phase_space_sphere_logical = new G4LogicalVolume(phase_space_sphere, nullptr, "PhaseSpaceSphereLogical");
		auto placement = new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 95.*cm), m_phase_space_sphere_logical, "PhaseSpaceSphere", logical_world, false, 0);
		placement->SetRotation(new G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), -90. *deg));
	}

	void ConstructSD() {
		G4MultiFunctionalDetector* sd = new G4MultiFunctionalDetector(fWorldName);
		SetSensitiveDetector(m_phase_space_sphere_logical, sd);
		LinacPhaseSpaceSphereScorer* phsp_sphere = new LinacPhaseSpaceSphereScorer(fWorldName);
		sd->RegisterPrimitive(phsp_sphere);
	}
};
*/

class G4VisualizationMessenger : public G4UImessenger {
public:
	G4VisualizationMessenger(G4UImanager* ui, const Plan& plan, LinacDetector* detector, LinacSource* source)
		: m_ui(ui),  m_plan(plan), m_detector(detector), m_source(source)
	{
		new G4UIdirectory("/simulation_points/");
		m_simulation_point_list_cmd = new G4UIcmdWithoutParameter("/simulation_points/list", this);
		m_simulation_point_list_cmd->SetGuidance("List all simulation points");
		m_simulation_point_set_cmd = new G4UIcmdWithAnInteger("/simulation_points/set", this);
		m_simulation_point_set_cmd->SetGuidance("Set TreatmenHead state to given simulation point index");
		m_simulation_point_set_cmd->SetParameterName("simulation_point_index", false);
		m_simulation_point_next_cmd = new G4UIcmdWithoutParameter("/simulation_points/next", this);
		m_simulation_point_next_cmd->SetGuidance("Set TreatmenHead state to next simulation point index");
		m_simulation_point_current_cmd = new G4UIcmdWithAnInteger("/simulation_points/current", this);
		m_simulation_point_current_cmd->SetGuidance("Print information about current simulation point");
		m_simulation_point_current_cmd->SetGuidance("Level 0: only index");
		m_simulation_point_current_cmd->SetGuidance("Level 1: + rotations + jaws/proximal mlc");
		m_simulation_point_current_cmd->SetGuidance("Level 2: + mlc/distal mlc");
		m_simulation_point_current_cmd->SetParameterName("print_level", true);
		m_simulation_point_current_cmd->SetDefaultValue(0);
		m_simulation_point_idx = 0;
		m_detector->setState(m_plan, m_simulation_point_idx);
		auto angles = m_plan.getRotations(m_simulation_point_idx);
		m_source->setAngles(angles.gantry_angle, angles.collimator_angle);
		new G4UIdirectory("/linac_debug/");
		m_linac_debug_trajectory_cmd = new G4UIcmdWithoutParameter("/linac_debug/rich_trajectories", this);
		m_linac_debug_trajectory_cmd->SetGuidance("Sets rich trajectories for debugging");
	}

	void setSimulationPoint(size_t point_idx) {
		if (!(point_idx < m_plan.getNumSimulationPoints())) {
			G4cerr << "Simulation point outside the range: " << std::to_string(point_idx) << G4endl;
			return;
		}
		m_detector->setState(m_plan, point_idx);
		auto angles = m_plan.getRotations(point_idx);
		m_source->setAngles(angles.gantry_angle, angles.collimator_angle);
		m_simulation_point_idx = point_idx;
		m_ui->ApplyCommand("/vis/drawVolume"); //refresh view
	}

	void printSimulationPoint(size_t point_idx, int level) {
		auto sp_string = m_plan.getSimulationPointDescription(level, point_idx);
		G4cout << sp_string << G4endl;
		/*
		G4cout << "Simulation point index: " << point_idx << G4endl;
		std::ios old_state(nullptr);
		old_state.copyfmt(std::cout);
		if (level < 1) return;
		if (std::holds_alternative<truebeam::TreatmentHeadDetector*>(m_treatment_head_detector)) {
			const auto& state = std::get<THDynamicStateJawAndMLC>(m_simulation_points.at(point_idx).th_state);
			G4cout << "Rotations (degrees):" << G4endl;
			G4cout << "\tCollimator: " << state.th_rotations.collimator_angle / deg << G4endl;
			G4cout << "\tGantry: " << state.th_rotations.gantry_angle / deg << G4endl;
			G4cout << G4endl << "Jaws (mm):" << G4endl;
			G4cout << "\tX1: " << state.jaw_positions.x1 / mm << G4endl;
			G4cout << "\tX2: " << state.jaw_positions.x2 / mm << G4endl;
			G4cout << "\tY1: " << state.jaw_positions.y1 / mm << G4endl;
			G4cout << "\tY2: " << state.jaw_positions.y2 / mm << G4endl;
			if (level > 1) {
				G4cout << "MLC leaf positions (mm): " << G4endl;
				G4cout << "Index BankA BankB" << G4endl;
				const auto& bank_X2 = state.leaf_positions.bank_X2;
				const auto& bank_X1 = state.leaf_positions.bank_X1;
				const size_t num_leaves = state.leaf_positions.bank_X2.size();
				std::cout << std::fixed << std::setprecision(1);
				for (size_t i = 0; i < num_leaves; ++i) {
					G4cout << '\t' << std::setfill('0') << std::setw(2) << i << " "  << std::setfill(' ') << std::setw(5) << bank_X2.at(i) << "  " << std::setw(5) << bank_X1.at(i) << std::endl;
				}
			}
			std::cout.copyfmt(old_state);
		}
		else if (std::holds_alternative<halcyon::TreatmentHeadDetector*>(m_treatment_head_detector)) {
			const auto& state = std::get<THDynamicStateDualLayerMLC>(m_simulation_points.at(point_idx).th_state);
			G4cout << "Rotations (degrees):" << G4endl;
			G4cout << "\tCollimator: " << state.th_rotations.collimator_angle / deg << G4endl;
			G4cout << "\tGantry: " << state.th_rotations.gantry_angle / deg << G4endl;
			if (level > 1) {
				G4cout << "Proximal MLC leaf positions (mm): " << G4endl;
				G4cout << "Index BankA BankB" << G4endl;
				const auto& proximal_bank_X2 = state.proximal_leaf_positions.bank_X2;
				const auto& proximal_bank_X1 = state.proximal_leaf_positions.bank_X1;
				const size_t num_leaves_proximal = state.proximal_leaf_positions.bank_X2.size();
				std::cout << std::fixed << std::setprecision(1);
				for (size_t i = 0; i < num_leaves_proximal; ++i) {
					G4cout << '\t' << std::setfill('0') << std::setw(2) << i << " "  << std::setfill(' ') << std::setw(5) << proximal_bank_X2.at(i) << "  " << std::setw(5) << proximal_bank_X1.at(i) << std::endl;
				}
				G4cout << "Distal MLC leaf positions (mm): " << G4endl;
				G4cout << "Index BankA BankB" << G4endl;
				const auto& distal_bank_X2 = state.distal_leaf_positions.bank_X2;
				const auto& distal_bank_X1 = state.distal_leaf_positions.bank_X1;
				const size_t num_leaves_distal = state.distal_leaf_positions.bank_X2.size();
				std::cout << std::fixed << std::setprecision(1);
				for (size_t i = 0; i < num_leaves_distal; ++i) {
					G4cout << '\t' << std::setfill('0') << std::setw(2) << i << " "  << std::setfill(' ') << std::setw(5) << distal_bank_X2.at(i) << "  " << std::setw(5) << distal_bank_X1.at(i) << std::endl;
				}
			}
			std::cout.copyfmt(old_state);
		}
		else {
			G4cout << "Treatment head not currently supported" << G4endl;
		}
		*/
	}

	void SetNewValue(G4UIcommand* command, G4String newValues) {
		const size_t num_pts = m_plan.getNumSimulationPoints();
		if (command == m_simulation_point_list_cmd) {
			G4cout << "Number of simulation points: " << num_pts << G4endl;
		}
		else if (command == m_simulation_point_set_cmd) {
			size_t idx = static_cast<size_t>(m_simulation_point_set_cmd->ConvertToInt(newValues));
			if (!(idx < num_pts)) {
				G4cerr << "Simulation point outside the range: " << std::to_string(idx) << G4endl;
				return;
			}
			m_detector->setState(m_plan, idx);
			m_simulation_point_idx = idx;
			m_ui->ApplyCommand("/vis/drawVolume"); //refresh view
		}
		else if (command == m_simulation_point_next_cmd) {
			setSimulationPoint(m_simulation_point_idx + 1);
		}
		else if (command == m_simulation_point_current_cmd) {
			const int level = m_simulation_point_set_cmd->ConvertToInt(newValues);
			printSimulationPoint(m_simulation_point_idx, level);
		}
		else if (command == m_linac_debug_trajectory_cmd) {
			m_ui->ApplyCommand("/vis/scene/add/trajectories");
			m_ui->ApplyCommand("/vis/scene/add/hits");
			m_ui->ApplyCommand("/vis/scene/endOfEventAction accumulate");
			m_ui->ApplyCommand("/vis/scene/add/trajectories rich");
			m_ui->ApplyCommand("/vis/modeling/trajectories/create/drawByCharge");
			m_ui->ApplyCommand("/vis/modeling/trajectories/drawByCharge-0/default/setDrawStepPts true");
			m_ui->ApplyCommand("/vis/modeling/trajectories/drawByCharge-0/default/setStepPtsSize 2");
		}
	}
	G4String GetCurrentValue(G4UIcommand* command) {
		if (command == m_simulation_point_set_cmd) {
			return std::to_string(m_simulation_point_idx);
		}
		return "";
	}

	G4UImanager* m_ui;
	const Plan& m_plan;
	LinacDetector* m_detector;
	LinacSource* m_source;
	size_t m_simulation_point_idx;
	G4UIcmdWithoutParameter* m_simulation_point_list_cmd;
	G4UIcmdWithAnInteger* m_simulation_point_set_cmd;
	G4UIcmdWithoutParameter* m_simulation_point_next_cmd;
	G4UIcmdWithAnInteger* m_simulation_point_current_cmd;
	G4UIcmdWithoutParameter* m_linac_debug_trajectory_cmd;
};

std::unordered_map<std::string, PhysicsProcess> createPhysicsProcessMapping() {
	std::unordered_map<std::string, PhysicsProcess> map;
	size_t num_processes = static_cast<size_t>(PhysicsProcess::NUM_PROCESSES);
	for (size_t i = 0; i < num_processes; ++i) {
		map[g_physics_process_names[i]] = static_cast<PhysicsProcess>(i);
	}
	return map;
}

std::string convertMillisecondsToReadable(uint64_t ms) {
	uint64_t seconds = ms/1000;
	uint64_t minutes = seconds / 60;
	uint64_t hours = minutes / 60;
	seconds -= 60*minutes;
	minutes -= hours*60;
	//seconds = seconds - static_cast<uint64_t>(minutes)*60.0;
	//double hours = minutes / 60.0;
	//minutes = minutes - static_cast<uint64_t>(hours)*60.0;
	return fmt::format("{} h : {} m : {} s", hours, minutes, seconds);
}

void runG4Simulation(
	const SimulationInputParams& input_params,
	const SimulationOutputParams& output_params,
	LinacDetector* detector,
	LinacSource* source,
	const Plan& plan,
	uint64_t num_original_histories,
	bool visualize, 
	Timing& timing) 
{
	auto const start_init = std::chrono::high_resolution_clock::now();
	//std::stringstream buffer;
	//G4cout.set_rdbuf(buffer.rdbuf());
	G4UIExecutive* ui = nullptr;
	G4Random::setTheSeed(input_params.seed);
	if(visualize)
	{
		char* str = "config.exe";
		ui = new G4UIExecutive(1, &str);
	}
#ifdef G4MULTITHREADED
	G4MTRunManager* runManager = new G4MTRunManager;
	runManager->SetNumberOfThreads(input_params.num_threads);
#else
	G4RunManager* runManager = new G4RunManager;
#endif
	G4UImanager* pUI = G4UImanager::GetUIpointer();
	const std::string geant_version = runManager->GetVersionString();
	G4String physics_name = input_params.physics_model;
	G4PhysListFactory factory;
	G4VModularPhysicsList* physics_list = factory.GetReferencePhysList(physics_name);
	//pUI->ApplyCommand("/process/em/verbose 0");
	//pUI->ApplyCommand("/process/had/verbose 0");
	
	const std::unordered_map<int, TraversedGeometry>* id_to_traversed = &detector->getTraversedGeometryMapping();

	if (input_params.phantom.has_value()) {
		const Phantom& phantom_config = input_params.phantom.value(); 
		G4VUserParallelWorld* parallel_world = nullptr;
		G4String name;
		if(phantom_config.type == PhantomType::WATER_BOX) {
			name = "WaterBoxPhantom";
			parallel_world = new LinacWaterPhantom(name, 
				std::get<WaterBoxPhantom>(phantom_config.def));
		}
		else if(phantom_config.type == PhantomType::WATER_CYLINDER) {
			name = "WaterCylinderPhantom";
			parallel_world = new LinacWaterCylinderPhantom(name,
				std::get<WaterCylinderPhantom>(phantom_config.def));
		}
		else if(phantom_config.type == PhantomType::VOXEL_BOX) {
			name = "VoxelBoxPhantom";
			parallel_world = new LinacVoxelPhantom(
				name, std::get<VoxelBoxPhantom>(phantom_config.def));
		}
		else {
			static_assert(static_cast<size_t>(PhantomType::NUM_PHANTOM_TYPES) == 3);
		}
		detector->RegisterParallelWorld(parallel_world);
		physics_list->RegisterPhysics(new G4ParallelWorldPhysics(name, true));
	}
	

	OutputInfo out_info{};
	out_info.monitor_chamber_dose = 0.0;
	for (const auto& dose_grid : output_params.dose_grids) {
		out_info.dose_grids.push_back({});
		auto& grid = out_info.dose_grids.back();
		const size_t num_voxels = static_cast<size_t>(dose_grid.num_voxels[0] * dose_grid.num_voxels[1] * dose_grid.num_voxels[2]);
		grid.resize(num_voxels, 0.0);
	}
	//We construct sensitive detectors
	SensitiveDetectorNames sd_names{};
	sd_names.monitor_chamber = detector->getMonitorChamberName();
	//DoseGrids
	for (size_t i = 0; i < output_params.dose_grids.size(); ++i) {
		const auto& config = output_params.dose_grids[i];
		G4String name = "DoseGrid_" + std::to_string(i);
		auto parallel_world = new LinacDoseGrid(name, config);
		detector->RegisterParallelWorld(parallel_world);
		physics_list->RegisterPhysics(new G4ParallelWorldPhysics(name, false));
		sd_names.dose_grids.push_back(name);
	}
	//PhaseSpacePlanes
	std::vector<LinacPhaseSpaceWriterInfo> phsp_writer_infos;
	std::vector<G4Mutex> phsp_writer_mutexes(output_params.phase_space_plane.size());
	for (size_t i = 0; i < output_params.phase_space_plane.size(); ++i) {
		const auto& config = output_params.phase_space_plane[i];
		const G4String name = "PhaseSpacePlane_" + std::to_string(i);
		phsp_writer_infos.emplace_back(LinacPhaseSpaceWriterInfo{
			config.phsp_filename, &phsp_writer_mutexes[i], input_params.track_particles, config.kill_particles
		});
		auto parallel_world = new LinacPhaseSpacePlane(name, config, phsp_writer_infos[i]);
		detector->RegisterParallelWorld(parallel_world);
		physics_list->RegisterPhysics(new G4ParallelWorldPhysics(name, false));
		sd_names.phaseplanes.push_back(name);
		out_info.phsp_filenames.push_back(config.phsp_filename);
	}
	//Phasespace killing sphere
	//const G4String sphere_name = "PhaseSpaceSphere";
	//auto parallel_world = new LinacPhaseSpaceSphere(sphere_name);
	//detector->RegisterParallelWorld(parallel_world);
	//physics_list->RegisterPhysics(new G4ParallelWorldPhysics(sphere_name, false));
	//sd_names.sphere_scorer = sphere_name;

	std::cout << "Initializing" << std::endl;
	auto name_to_physics_process = createPhysicsProcessMapping();
	//Gather linac info to structs
	const LinacInitializationInfo info{
		sd_names,
		source,
		&out_info,
		input_params.track_particles,
		id_to_traversed,
		&name_to_physics_process
	};

	runManager->SetUserInitialization(detector);
	runManager->SetUserInitialization(physics_list);
	runManager->SetUserInitialization(new LinacActionInitialization(info));
	runManager->Initialize();

	const auto bm = input_params.biasing_mode;
	if (bm == BiasingMode::NONE) {
		pUI->ApplyCommand("/run/setCut 100 um");
    pUI->ApplyCommand("/tracking/verbose 1");
		std::cout << "Biasing mode: NONE" << std::endl;
	}
	else if (bm == BiasingMode::TARGET) {
		pUI->ApplyCommand("/run/setCut 100 um");
		pUI->ApplyCommand("/process/em/setSecBiasing eBrem target 100 100 MeV");
		pUI->ApplyCommand("/process/em/setDirectionalSplitting true");
		pUI->ApplyCommand("/process/em/setDirectionalSplittingTarget 0 0 0 mm");
		pUI->ApplyCommand("/process/em/setDirectionalSplittingRadius 40 cm");
		std::cout << "Biasing mode: eBrem in TARGET" << std::endl;
	}
	else if (bm == BiasingMode::LIGHTFIELD){
		pUI->ApplyCommand("/run/setCut 1 m");
		pUI->ApplyCommand("/run/setCutForRegion target 100 um");
		std::cout << "Biasing mode: LIGHTFIELD" << std::endl;
	}
	else if (bm == BiasingMode::DIRECTIONAL_SPLITTING){
		pUI->ApplyCommand("/run/setCut 100 um");
		pUI->ApplyCommand("/process/em/setDirectionalSplitting true");
		pUI->ApplyCommand("/process/em/setDirectionalSplittingTarget 0 0 0 mm");
		pUI->ApplyCommand("/process/em/setDirectionalSplittingRadius 40 cm");
		pUI->ApplyCommand("/process/em/setSecBiasing eBrem world 100 100 MeV");
		pUI->ApplyCommand("/process/em/setSecBiasing Rayl world 100 100 MeV");
		pUI->ApplyCommand("/process/em/setSecBiasing phot world 100 100 MeV");
		pUI->ApplyCommand("/process/em/setSecBiasing compt world 100 100 MeV");
		pUI->ApplyCommand("/process/em/setSecBiasing annihil world 100 100 MeV");
		
		std::cout << "Biasing mode: Directional splitting" << std::endl;
	}
	else {
		throw std::runtime_error("Unable to determine biasing mode");
	}
	for (const auto& phsp_plane : output_params.phase_space_plane) {
		LinacPhaseSpaceWriter::InitFile(phsp_plane.phsp_filename);
	}

	const auto init_duration = std::chrono::high_resolution_clock::now() - start_init;
	timing.initialization = std::chrono::duration_cast<std::chrono::milliseconds>(init_duration).count();
	if (visualize) {
		std::cout << "Starting visualization" << std::endl;
		new G4VisualizationMessenger(pUI, plan, detector, source);
		G4VisManager* visManager = new G4VisExecutive;
		visManager->Initialize();
		pUI->ApplyCommand("/vis/open OGLSQt 1024x768-0+0");
		pUI->ApplyCommand("/vis/drawVolume");
    pUI->ApplyCommand("/vis/viewer/set/viewpointVector 1 0 0");
    pUI->ApplyCommand("/vis/scene/add/trajectories smooth");
    pUI->ApplyCommand("/vis/scene/endOfEventAction accumulate");
		ui->SessionStart();
		return;
	}
	for (const auto& cmd : input_params.geant4_commands) {
		pUI->ApplyCommand(cmd);
	}
	size_t num_simulation_points = plan.getNumSimulationPoints();
	uint64_t total_num_particles = plan.getTotalNumParticles();
	uint64_t done_num_particles = 0;
	double estimated_time_per_particle = 0.0;
	for (size_t pt_idx = 0; pt_idx < num_simulation_points; ++pt_idx) {
		auto const start_state = std::chrono::high_resolution_clock::now();
		detector->setState(plan, pt_idx);
		const auto& rotations = plan.getRotations(pt_idx);
		source->setAngles(rotations.gantry_angle, rotations.collimator_angle);
		const auto state_duration = std::chrono::high_resolution_clock::now() - start_state;
		timing.state_transitions += std::chrono::duration_cast<std::chrono::milliseconds>(state_duration).count();
		auto num_particles = plan.getNumParticles(pt_idx);
		std::string estimated_time_remaining_msg = "???";
		if(done_num_particles > 0) {
			double estimated_time_in_ms = static_cast<double>(total_num_particles - done_num_particles) * estimated_time_per_particle;
			estimated_time_remaining_msg = convertMillisecondsToReadable(static_cast<uint64_t>(estimated_time_in_ms));
		}
			
		std::cout << fmt::format("Simulation point {}, Number of particles: {}, Estimated time remaining {}\n",
			pt_idx, num_particles, estimated_time_remaining_msg);
		auto const start_beam = std::chrono::high_resolution_clock::now();
		runManager->BeamOn(static_cast<G4int>(num_particles));
		done_num_particles += num_particles;
		const auto beam_duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start_beam).count();
		double v = static_cast<double>(beam_duration) / static_cast<double>(num_particles);
		estimated_time_per_particle += (v - estimated_time_per_particle)/(pt_idx + 1);
		timing.beam_on += beam_duration;
	}
	std::cout << fmt::format("Estimated particle/time: {:.1f} particles/s\n", 1000.0/estimated_time_per_particle);
	for (size_t i = 0; i < output_params.dose_grids.size(); ++i) {
		writeDoseGridBinary(
			output_params.dose_grids[i], 
			num_original_histories, 
			out_info.dose_grids[i], 
			out_info.monitor_chamber_dose);
	}
	for (const auto& phsp_plane : output_params.phase_space_plane) {
		LinacPhaseSpaceWriter::WriteBinaryHeaderToFile(
			phsp_plane.phsp_filename, 
			input_params.track_particles, 
			num_original_histories);
	}
	delete runManager;
}
