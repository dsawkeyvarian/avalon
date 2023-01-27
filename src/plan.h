#pragma once
#include <vector>
#include <variant>
#include "linac.h"
#include "external/json_fwd.hpp"
using json = nlohmann::json;

class Logger;

struct LeafPositions {
  std::vector<double> bank_X1;
  std::vector<double> bank_X2;
};

struct JawPositions {
  double x1;
  double x2;
  double y1;
  double y2;
};

struct TreatmentHeadRotations {
  double gantry_angle;
  double collimator_angle;
};

struct THDynamicStateNoHead {};

struct THDynamicStateJawAndMLC {
  TreatmentHeadRotations th_rotations;
  JawPositions jaw_positions;
  LeafPositions leaf_positions;
};

struct THDynamicStateDualLayerMLC {
  TreatmentHeadRotations th_rotations;
  LeafPositions proximal_leaf_positions;
  LeafPositions distal_leaf_positions;
};

class Plan
{
public:
    
    size_t size() const {return m_num_particles_sp.size();}
    size_t getNumSimulationPoints() const {return size();}
    uint64_t getNumParticles(size_t idx) const {return m_num_particles_sp.at(idx);}
    uint64_t getTotalNumParticles() const;
    std::string getSimulationPointDescription(int, size_t) const { return ""; }
    const TreatmentHeadRotations& getRotations(size_t idx) const;
    const THDynamicStateJawAndMLC& getJawAndMLC(TreatmentHeadType type, size_t idx) const;
    const THDynamicStateDualLayerMLC& getDualLayerMLC(TreatmentHeadType type, size_t idx) const;
    //SimulationPointNone getSPNone(TreatmentHeadType type, size_t idx) const;

    static std::unique_ptr<Plan> parseSimulationPointPlan(
        TreatmentHeadType th_type, 
        MLCType mlc_type,
        uint64_t num_particles, 
        const json& sp_plan,
        ILogger& log);
    static std::unique_ptr<Plan> parseControlPointPlan(
        TreatmentHeadType th_type, 
        MLCType nlc_type,
        uint64_t num_particles, 
        const json& cp_plan, 
        ILogger& log,
        uint64_t seed,
        uint64_t num_simulation_points);

private:
    TreatmentHeadType m_type;
    MLCType m_mlc_type;
    std::vector<uint64_t> m_num_particles_sp;
    using THDynamicState = std::variant<THDynamicStateJawAndMLC, THDynamicStateDualLayerMLC, TreatmentHeadRotations>;
    std::vector<THDynamicState> m_th_states;

    Plan(
        TreatmentHeadType type,
        MLCType mlc_type,
        std::vector<uint64_t> num_particles_sp,
        const std::vector<THDynamicState>& th_states
    )
    : m_type(type), m_mlc_type(mlc_type), m_num_particles_sp(num_particles_sp),
    m_th_states(th_states)  {}

    static std::vector<THDynamicState> parseStates(
        const json &plan, TreatmentHeadType& th_type);
    static bool validateStates( 
        const std::vector<THDynamicState>& states, //Could be span
         TreatmentHeadType th_type,
         MLCType mlc_type,
         ILogger& log);

    static void interpolateControlPoint(
        TreatmentHeadType th_type, 
        double t, 
        Plan::THDynamicState& sp,
        const THDynamicState& cp0, 
        const THDynamicState& cp1);
};


