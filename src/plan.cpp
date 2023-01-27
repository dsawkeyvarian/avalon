#include "plan.h"
#include "machine_information.h"
#include "logger.h"
#include <random>
#include <algorithm>

#include <fmt/core.h>
#define JSON_DIAGNOSTICS 1
#include "external/json.hpp"
using json = nlohmann::json;
#include "parse.h"

namespace
{
  TreatmentHeadRotations parseRotations(const json &j)
  {
    TreatmentHeadRotations r = {};
    const auto unit = j["unit"].get<AngleUnit>().scaler;
    // TODO(jouko): wrap angles
    r.collimator_angle = unit * j["collimator"].get<double>();
    r.gantry_angle = unit * j["gantry"].get<double>();
    return r;
  }
  JawPositions parseJawPositions(const json &j) {
    JawPositions jaw = {};
    const auto unit = j["unit"].get<LengthUnit>().scaler;
    jaw.x1 = unit * j["x1"].get<double>();
    jaw.x2 = unit * j["x2"].get<double>();
    jaw.y1 = unit * j["y1"].get<double>();
    jaw.y2 = unit * j["y2"].get<double>();
    return jaw;
  }
  LeafPositions parseLeafPositions(const json &j) {
    LeafPositions leafs = {};
    const auto unit = j["unit"].get<LengthUnit>().scaler;
    leafs.bank_X2 = j["bank_X2"].get<std::vector<double>>();
    leafs.bank_X1 = j["bank_X1"].get<std::vector<double>>();
    for (size_t i = 0; i < leafs.bank_X2.size(); ++i)
      leafs.bank_X2[i] *= unit;
    for (size_t i = 0; i < leafs.bank_X1.size(); ++i)
      leafs.bank_X1[i] *= unit;
    return leafs;
  }
}

uint64_t Plan::getTotalNumParticles() const {
  uint64_t total = 0;
  for(auto n : m_num_particles_sp )
    total += n;
  return total;
}

const TreatmentHeadRotations& Plan::getRotations(size_t idx) const {
  const auto& pt = m_th_states.at(idx);
  if(std::holds_alternative<THDynamicStateJawAndMLC>(pt)) {
    return std::get<THDynamicStateJawAndMLC>(pt).th_rotations;
  }
  else if(std::holds_alternative<THDynamicStateDualLayerMLC>(pt)) {
    return std::get<THDynamicStateDualLayerMLC>(pt).th_rotations;
  }
  else if(std::holds_alternative<TreatmentHeadRotations>(pt)) {
    return std::get<TreatmentHeadRotations>(pt);
  }
  else {
    //Can we get variant count for static_assert ?
    throw std::logic_error("Trying to get rotations for state whcih has no rotations");
  }
}

const THDynamicStateJawAndMLC& Plan::getJawAndMLC(TreatmentHeadType type, size_t idx) const
{
  if(type != m_type)
    throw std::logic_error("Trying to get incompatible state with TreatmentHeadType from plan");
  return std::get<THDynamicStateJawAndMLC>(m_th_states.at(idx));
}
const THDynamicStateDualLayerMLC&  Plan::getDualLayerMLC(TreatmentHeadType type, size_t idx) const
{
  if(type != m_type)
    throw std::logic_error("Trying to get incompatible state with TreatmentHeadType from plan");
  return std::get<THDynamicStateDualLayerMLC>(m_th_states.at(idx));
}

std::unique_ptr<Plan> Plan::parseSimulationPointPlan(
  TreatmentHeadType th_type, 
  MLCType mlc_type,
  uint64_t num_particles, 
  const json& sp_plan, ILogger& log)
{
  const size_t num_sp = sp_plan.size();
  std::vector<double> weights_sp(num_sp);
  double total_weight = 0.0;
  for (size_t i = 0; i < num_sp; ++i) {
    double weight = sp_plan.at(i)["weight"].get<double>();
    total_weight += weight;
    weights_sp[i] = weight;
  }
  std::vector<uint64_t> num_particles_sp(num_sp);
  uint64_t total_num_particles = 0;
  for(size_t i = 0; i < num_sp; ++i) {
    double relative = weights_sp[i]/total_weight;
    uint64_t absolute = static_cast<uint64_t>(relative * num_particles);
    num_particles_sp[i] = absolute;
    total_num_particles += absolute;
  }
  //PhaseSpace contains only certain number of particles so we want always to
  //underestimate rather than overestimate
  if(total_num_particles > num_particles)
    throw std::logic_error("Unable to divide num_particles into simulation points");
  const auto states = parseStates(sp_plan, th_type);
  bool ok = validateStates(states, th_type, mlc_type, log);
  if(!ok)
    throw std::runtime_error("Given simulation points are not a valid plan");
  return std::unique_ptr<Plan>( 
    new Plan( th_type, mlc_type, num_particles_sp, 
    parseStates(sp_plan, th_type))
  );
}

std::unique_ptr<Plan> Plan::parseControlPointPlan(
  TreatmentHeadType th_type, 
  MLCType mlc_type,
  uint64_t num_particles, 
  const json& cp_plan, ILogger& log,
  uint64_t seed,
  uint64_t num_simulation_points)
{
  const size_t num_cp = cp_plan.size();
  if(num_cp < 2)
    throw std::runtime_error("Control point plan must have at least two points");
  const size_t num_sp = num_simulation_points;
  uint64_t particles_per_sample = num_particles/num_sp;
  //const size_t num_sp = static_cast<size_t>(num_particles/particles_per_sample);
  const auto control_points = parseStates(cp_plan, th_type);
  bool ok = validateStates(control_points, th_type, mlc_type, log);
  if(!ok)
    throw std::runtime_error("Given control points are not valid plan");
  std::vector<double> cmw_cp(num_cp);
  cmw_cp[0] = cp_plan.at(0).at("cumulative_meterset_weight").get<double>();
  for(size_t i = 1; i < num_cp; ++i) {
    double cmw = cp_plan.at(i).at("cumulative_meterset_weight").get<double>();
    cmw_cp[i] = cmw;
    if(cmw < cmw_cp[i-1])
      throw std::runtime_error("Control points are not properly ordered. Cumulative meterset weights are not increasing");
  }
  if(cmw_cp.front() != 0.0)
    throw std::runtime_error("First control point does not have cumulative meterset weight of 0.0");
  if(cmw_cp.back() != 1.0)
    throw std::runtime_error("Last control point does not have cumulate meterset weight of 1.0");

  //Sample simulation points from control points wrt cumulative meter weight
  std::vector<double> samples(num_sp);
  std::mt19937_64 re(seed);
  std::uniform_real_distribution<double> distrib(0.0, 1.0);
  for(size_t i = 0; i < num_sp; ++i) {
    samples[i] = distrib(re);
  }
  std::sort(samples.begin(), samples.end());
  std::vector<THDynamicState> simulation_points(num_sp);
  for(size_t i = 0; i < num_sp; ++i) {
    double s = samples[i];
    for(size_t j = 1; j < num_cp; ++j) {
      if((s >= cmw_cp[j-1]) && (s < cmw_cp[j])) {
        double width = cmw_cp[j] - cmw_cp[j-1];
        double t = (s - cmw_cp[j-1])/width;
        interpolateControlPoint(th_type, t, simulation_points[i], control_points[j-1], control_points[j]);
      }
    }
  }
  std::vector<uint64_t> num_particles_sp(num_sp, particles_per_sample);
  return std::unique_ptr<Plan>( 
    new Plan( th_type, mlc_type, num_particles_sp, 
    simulation_points)
  );
}

std::vector<Plan::THDynamicState> Plan::parseStates(const json& plan, TreatmentHeadType& th_type)
{
  std::vector<Plan::THDynamicState> states(plan.size());
  for (size_t i = 0; i < plan.size(); ++i) {
    const auto& sp = plan.at(i);
    if(th_type == TreatmentHeadType::TRUEBEAM) {
      states[i] = THDynamicStateJawAndMLC{
        parseRotations(sp.at("angles")),
        parseJawPositions(sp.at("jaw_positions")),
        parseLeafPositions(sp.at("mlc_leaf_positions"))
      };
    }
    else if(th_type == TreatmentHeadType::HALCYON)
    {
      states[i] = THDynamicStateDualLayerMLC{
        parseRotations(sp.at("angles")),
        parseLeafPositions(sp.at("mlc_proximal_leaf_positions")),
        parseLeafPositions(sp.at("mlc_distal_leaf_positions"))
      };
    }
    else if(th_type == TreatmentHeadType::AVALON)
    {
      states[i] = THDynamicStateDualLayerMLC{
        parseRotations(sp.at("angles")),
        parseLeafPositions(sp.at("mlc_proximal_leaf_positions")),
        parseLeafPositions(sp.at("mlc_distal_leaf_positions"))
      };
    }
    else if(th_type == TreatmentHeadType::AVALON_ELECTRON)
    {
      states[i] = THDynamicStateDualLayerMLC{
        parseRotations(sp.at("angles")),
        parseLeafPositions(sp.at("mlc_proximal_leaf_positions")),
        parseLeafPositions(sp.at("mlc_distal_leaf_positions"))
      };
    }
    else if(th_type == TreatmentHeadType::NOHEAD) {
      auto rotations = TreatmentHeadRotations{0.0, 0.0};
      if(sp.contains("angles"))
        rotations = parseRotations(sp.at("angles"));
      states[i] = rotations;
    }
    else {
      static_assert(static_cast<int>(TreatmentHeadType::NUM_TREATMENT_HEAD_TYPES) == 5);
    }
  }
  return states;
}

bool Plan::validateStates( 
  const std::vector<THDynamicState>& states, //Could be span
  TreatmentHeadType th_type, MLCType mlc_type, ILogger& log)
{
  bool ok = true;
  size_t num_pts = states.size();
  if(th_type == TreatmentHeadType::TRUEBEAM) {
    for(size_t i = 0; i < num_pts; ++i) {
      const auto& pt = std::get<THDynamicStateJawAndMLC>(states[i]);
      bool leaf_number_ok = true;
      if(mlc_type == MLCType::HD) {
        leaf_number_ok &= pt.leaf_positions.bank_X1.size() == truebeam::hd::num_leaves_in_bank;
        leaf_number_ok &= pt.leaf_positions.bank_X2.size() == truebeam::hd::num_leaves_in_bank;
      }
      else if(mlc_type == MLCType::MILLENNIUM120) {
        leaf_number_ok &= pt.leaf_positions.bank_X1.size() == truebeam::millennium120::num_leaves_in_bank;
        leaf_number_ok &= pt.leaf_positions.bank_X2.size() == truebeam::millennium120::num_leaves_in_bank;
      }
      else if(mlc_type == MLCType::NONE) {
        leaf_number_ok &= pt.leaf_positions.bank_X1.size() == 0;
        leaf_number_ok &= pt.leaf_positions.bank_X2.size() == 0;
      }
      else {
        static_assert(truebeam::mlcs.size() == 3);
      }
      if(!leaf_number_ok) {
          log.error(fmt::format("Wrong number of leaves at simulation point {}", i));
          ok = false;
      }
    }
  }
  else if(th_type == TreatmentHeadType::HALCYON) {
    for(size_t i = 0; i < num_pts; ++i) {
      const auto& pt = std::get<THDynamicStateDualLayerMLC>(states[i]);
      bool leaf_number_ok = true;
      leaf_number_ok &= pt.distal_leaf_positions.bank_X1.size() == halcyon::num_distal_leaves_in_bank;
      leaf_number_ok &= pt.distal_leaf_positions.bank_X2.size() == halcyon::num_distal_leaves_in_bank;
      leaf_number_ok &= pt.proximal_leaf_positions.bank_X1.size() == halcyon::num_proximal_leaves_in_bank;
      leaf_number_ok &= pt.proximal_leaf_positions.bank_X2.size() == halcyon::num_proximal_leaves_in_bank;
      if(!leaf_number_ok) {
          log.error(fmt::format("Wrong number of leaves at simulation point {}", i));
          ok = false;
      }
    }
  }
  else if(th_type == TreatmentHeadType::AVALON) {
    for(size_t i = 0; i < num_pts; ++i) {
      const auto& pt = std::get<THDynamicStateDualLayerMLC>(states[i]);
      bool leaf_number_ok = true;
      leaf_number_ok &= pt.distal_leaf_positions.bank_X1.size() == avalon::num_distal_leaves_in_bank;
      leaf_number_ok &= pt.distal_leaf_positions.bank_X2.size() == avalon::num_distal_leaves_in_bank;
      leaf_number_ok &= pt.proximal_leaf_positions.bank_X1.size() == avalon::num_proximal_leaves_in_bank;
      leaf_number_ok &= pt.proximal_leaf_positions.bank_X2.size() == avalon::num_proximal_leaves_in_bank;
      if(!leaf_number_ok) {
          log.error(fmt::format("Wrong number of leaves at simulation point {}", i));
          ok = false;
      }
    }
  }
  else if(th_type == TreatmentHeadType::AVALON_ELECTRON) {
    for(size_t i = 0; i < num_pts; ++i) {
      const auto& pt = std::get<THDynamicStateDualLayerMLC>(states[i]);
      bool leaf_number_ok = true;
      leaf_number_ok &= pt.distal_leaf_positions.bank_X1.size() == avalon::num_distal_leaves_in_bank;
      leaf_number_ok &= pt.distal_leaf_positions.bank_X2.size() == avalon::num_distal_leaves_in_bank;
      leaf_number_ok &= pt.proximal_leaf_positions.bank_X1.size() == avalon::num_proximal_leaves_in_bank;
      leaf_number_ok &= pt.proximal_leaf_positions.bank_X2.size() == avalon::num_proximal_leaves_in_bank;
      if(!leaf_number_ok) {
          log.error(fmt::format("Wrong number of leaves at simulation point {}", i));
          ok = false;
      }
    }
  }


  else if(th_type == TreatmentHeadType::NOHEAD) {
    for(size_t i = 0; i < num_pts; ++i) {
      if(!std::holds_alternative<TreatmentHeadRotations>(states[i]))
        ok = false;
    }
  }
  else {
    static_assert(static_cast<int>(TreatmentHeadType::NUM_TREATMENT_HEAD_TYPES) == 5);
  }
  return ok;
}

void interpolate(double t, TreatmentHeadRotations& r,
  const TreatmentHeadRotations& p0, const TreatmentHeadRotations& p1)
{
  r.collimator_angle = (1-t) * p0.collimator_angle + t*p1.collimator_angle;
  r.gantry_angle = (1-t) * p0.gantry_angle + t*p1.gantry_angle;
}

void interpolate(double t, JawPositions& r,
  const JawPositions& p0, const JawPositions& p1)
{
  r.x1 = (1-t) * p0.x1 + t*p1.x1; 
  r.x2 = (1-t) * p0.x2 + t*p1.x2;
  r.y1 = (1-t) * p0.y1 + t*p1.y1;
  r.y2 = (1-t) * p0.y2 + t*p1.y2;
}

void interpolate(double t, LeafPositions& r, 
  const LeafPositions& p0, const LeafPositions& p1)
{
  r.bank_X1.resize(p0.bank_X1.size());
  r.bank_X2.resize(p0.bank_X2.size());
  for(size_t i = 0; i < r.bank_X1.size(); ++i)
    r.bank_X1[i] = (1-t)*p0.bank_X1[i] + t*p1.bank_X1[i];
  for(size_t i = 0; i < r.bank_X2.size(); ++i)
    r.bank_X2[i] = (1-t)*p0.bank_X2[i] + t*p1.bank_X2[i];
}

void Plan::interpolateControlPoint(
  TreatmentHeadType th_type, 
  double t, 
  Plan::THDynamicState& _sp,
  const THDynamicState& cp0, 
  const THDynamicState& cp1)
{
  if(th_type == TreatmentHeadType::TRUEBEAM) {
    const auto& p0 = std::get<THDynamicStateJawAndMLC>(cp0);
    const auto& p1 = std::get<THDynamicStateJawAndMLC>(cp1);
    _sp = THDynamicStateJawAndMLC{};
    auto& sp = std::get<THDynamicStateJawAndMLC>(_sp);
    interpolate(t, sp.th_rotations, p0.th_rotations, p1.th_rotations);
    interpolate(t, sp.jaw_positions, p0.jaw_positions, p1.jaw_positions);
    interpolate(t, sp.leaf_positions, p0.leaf_positions, p1.leaf_positions);
  }
  else if(th_type == TreatmentHeadType::HALCYON)
  {
    const auto& p0 = std::get<THDynamicStateDualLayerMLC>(cp0);
    const auto& p1 = std::get<THDynamicStateDualLayerMLC>(cp1);
    _sp = THDynamicStateDualLayerMLC{};
    auto& sp = std::get<THDynamicStateDualLayerMLC>(_sp);
    interpolate(t, sp.th_rotations, p0.th_rotations, p1.th_rotations);
    interpolate(t, sp.proximal_leaf_positions, p0.proximal_leaf_positions, p1.proximal_leaf_positions);
    interpolate(t, sp.distal_leaf_positions, p0.distal_leaf_positions, p1.distal_leaf_positions);
  }
  else if(th_type == TreatmentHeadType::AVALON)
  {
    const auto& p0 = std::get<THDynamicStateDualLayerMLC>(cp0);
    const auto& p1 = std::get<THDynamicStateDualLayerMLC>(cp1);
    _sp = THDynamicStateDualLayerMLC{};
    auto& sp = std::get<THDynamicStateDualLayerMLC>(_sp);
    interpolate(t, sp.th_rotations, p0.th_rotations, p1.th_rotations);
    interpolate(t, sp.proximal_leaf_positions, p0.proximal_leaf_positions, p1.proximal_leaf_positions);
    interpolate(t, sp.distal_leaf_positions, p0.distal_leaf_positions, p1.distal_leaf_positions);
  }
  else if(th_type == TreatmentHeadType::AVALON_ELECTRON)
  {
    const auto& p0 = std::get<THDynamicStateDualLayerMLC>(cp0);
    const auto& p1 = std::get<THDynamicStateDualLayerMLC>(cp1);
    _sp = THDynamicStateDualLayerMLC{};
    auto& sp = std::get<THDynamicStateDualLayerMLC>(_sp);
    interpolate(t, sp.th_rotations, p0.th_rotations, p1.th_rotations);
    interpolate(t, sp.proximal_leaf_positions, p0.proximal_leaf_positions, p1.proximal_leaf_positions);
    interpolate(t, sp.distal_leaf_positions, p0.distal_leaf_positions, p1.distal_leaf_positions);
  }

  else if(th_type == TreatmentHeadType::NOHEAD) {
    const auto& p0 = std::get<TreatmentHeadRotations>(cp0);
    const auto& p1 = std::get<TreatmentHeadRotations>(cp1);
    _sp = TreatmentHeadRotations{};
    auto& sp = std::get<TreatmentHeadRotations>(_sp);
    interpolate(t, sp, p0, p1);
  }
  else {
    throw std::logic_error("Unknown treatment head type while interpolating control points");
  }
}
