#include <string>
#include "G4SystemOfUnits.hh"

#define JSON_DIAGNOSTICS 1
#include "external/json.hpp"
using json = nlohmann::json;
#include "parse.h"
#include <fmt/core.h>

void from_json(const json& j, LengthUnit& unit) {
  const auto& str = j.get<std::string>();
  if (str.compare("mm") == 0)
    unit.scaler = mm;
  else if (str.compare("cm") == 0)
    unit.scaler = cm;
  else {
    throw nlohmann::detail::parse_error::create(302, 0, fmt::format("String {} is not LengthUnit", str), j);
  }
}

void from_json(const json& j, AngleUnit& unit) {
  const auto& str = j.get<std::string>();
  if (str.compare("deg") == 0)
    unit.scaler = degree;
  else if (str.compare("rad") == 0)
    unit.scaler = radian;
  else {
    throw nlohmann::detail::parse_error::create(302, 0, fmt::format("String {} is not AngleUnit", str), j);
  }
}

void from_json(const json& j, EnergyUnit& unit) {
  const auto& str = j.get<std::string>();
  if (str.compare("MeV") == 0)
    unit.scaler = MeV;
  else {
    throw nlohmann::detail::parse_error::create(302, 0, fmt::format("String {} is not EnergyUnit", str), j);
  }
}