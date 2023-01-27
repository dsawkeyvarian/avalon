#include "machine_information.h"

#include <fmt/core.h>
#include <fmt/format.h>

void runtimeValidationCombinationCheck() {
  //These trigger tests if something goes awry with matching valid combinations 
  //and other compilation time constants used in other parts of the code
  //Would like to have as comptilation time checks but
  //no patience for figuring C++ template magic
  for(const auto& combination : validCombinations) {
    if(combination.th_type == TreatmentHeadType::TRUEBEAM) {
      bool ok = combination.mlc_types.size() == truebeam::mlcs.size();
      ok |= combination.energy_modes.size() == truebeam::energy_modes.size();
      if(!ok)
        throw std::logic_error("Runtime checks for validation combinations do not pass");
    }
  }
}

bool isValidCombination(TreatmentHeadConfiguration config) {
  bool valid = false;
  for(const auto& combination : validCombinations) {
    if(combination.th_type != config.th_type)
      continue;
    for(auto mlc_type : combination.mlc_types) {
      if(mlc_type != config.mlc_type)
        continue;
      for(auto energy_mode : combination.energy_modes) {
        if(energy_mode == config.energy_mode) {
          valid = true;
          break;
        }
      }
      if(valid) break;
    }
    if(valid) break;
  }
  return valid;
}

std::string printValidCombinations() {
  auto out = fmt::memory_buffer();
  format_to(std::back_inserter(out), "Valid treatment head combinations in (TreatmentHeadType, MLCTypes, EnergyModes) format");
  for(const auto& combination : validCombinations) {
    format_to(std::back_inserter(out), "\n  ( {},", th_type_names[static_cast<int>(combination.th_type)]);
    std::string mlc_choices;
    for(auto mlc_type : combination.mlc_types)
      mlc_choices += fmt::format(" {} |", mlc_type_names[static_cast<int>(mlc_type)]);
    mlc_choices.pop_back();
    mlc_choices.pop_back();
    std::string energy_mode_choices;
    for(auto energy_mode : combination.energy_modes)
      energy_mode_choices += fmt::format(" {} |", energy_mode_names[static_cast<int>(energy_mode)]);
    energy_mode_choices.pop_back();
    energy_mode_choices.pop_back();
    format_to(std::back_inserter(out), "{},{} )", mlc_choices, energy_mode_choices);
  }
  return to_string(out);
}