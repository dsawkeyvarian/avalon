#pragma once
#include <array>
#include <vector>
#include <string>

#define SIZEOFARRAY(arr) (sizeof(arr)/sizeof(*arr))

enum class TreatmentHeadType
{
  NOHEAD,
  TRUEBEAM,
  HALCYON,
  AVALON,
  NUM_TREATMENT_HEAD_TYPES
};
static const char* th_type_names[] = {
  "no_head",
  "truebeam",
  "halcyon",
  "avalon"
};
static_assert(static_cast<size_t>(TreatmentHeadType::NUM_TREATMENT_HEAD_TYPES) == SIZEOFARRAY(th_type_names));

enum class MLCType {
  NONE,
  HD,
  MILLENNIUM120,
  RDS,
  AVALON,
  NUM_MLC
};
static const char* mlc_type_names[] = {
  "none",
  "hd",
  "millennium120",
  "rds",
  "avalon"
};
static_assert(static_cast<size_t>(MLCType::NUM_MLC) == SIZEOFARRAY(mlc_type_names));

enum class EnergyMode {
  NONE,
  X06_FFF,
  X06,
  X08,
  X08_FFF,
  X10,
  X10_FFF,
  X15,
  NUM_ENERGY_MODES
};
static const char* energy_mode_names[] = {
  "none",
  "6X-FFF",
  "6X",
  "08X",
  "08X-FFF",
  "10X",
  "10X-FFF",
  "15X"
};
static_assert(static_cast<size_t>(EnergyMode::NUM_ENERGY_MODES) == SIZEOFARRAY(energy_mode_names));

struct TreatmentHeadConfiguration {
  TreatmentHeadType th_type;
  MLCType mlc_type;
  EnergyMode energy_mode;
};

struct ValidCombination {
  TreatmentHeadType th_type;
  std::vector<MLCType> mlc_types;
  std::vector<EnergyMode> energy_modes;
};


static const std::array< ValidCombination, 4> validCombinations = {
   ValidCombination{TreatmentHeadType::TRUEBEAM, 
    {MLCType::NONE, MLCType::MILLENNIUM120, MLCType::HD}, 
    {EnergyMode::X06, EnergyMode::X06_FFF}},

   ValidCombination{TreatmentHeadType::HALCYON,
    {MLCType::RDS}, {EnergyMode::X06_FFF}},

   ValidCombination{TreatmentHeadType::AVALON,
    {MLCType::AVALON}, {EnergyMode::X06, EnergyMode::X06_FFF,EnergyMode::X10_FFF,EnergyMode::X10}},

   ValidCombination{TreatmentHeadType::NOHEAD,
    {MLCType::NONE}, {EnergyMode::NONE}}
};

void runtimeValidationCombinationCheck();
bool isValidCombination(TreatmentHeadConfiguration config);
std::string printValidCombinations();

namespace truebeam
{
  //These are separate lists due to C++ constraints, should match lists above
  static const std::array<MLCType, 3> mlcs = {MLCType::NONE, MLCType::MILLENNIUM120, MLCType::HD};
  static const std::array<EnergyMode, 2> energy_modes = {EnergyMode::X06, EnergyMode::X06_FFF};
  namespace millennium120
  {
    static const size_t num_leaves_in_bank = 60;
  }
  namespace hd
  {
    static const size_t num_leaves_in_bank = 60;
  }
  static const double max_jaw_overtravel = 20.0;
}

namespace halcyon
{
  const size_t num_proximal_leaves_in_bank = 29;
  const size_t num_distal_leaves_in_bank = 28;
}
namespace avalon
{
  static const std::array<EnergyMode, 4> energy_modes = {EnergyMode::X06, EnergyMode::X06_FFF,EnergyMode::X10, EnergyMode::X10_FFF};
  const size_t num_proximal_leaves_in_bank = 47;
  const size_t num_distal_leaves_in_bank = 46;
}