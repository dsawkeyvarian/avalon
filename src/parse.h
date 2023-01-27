#pragma once
#include <array>
#include "external/json_fwd.hpp"

struct LengthUnit {
  double scaler;
};

void from_json(const json& j, LengthUnit& unit);

struct AngleUnit {
  double scaler;
};

void from_json(const json& j, AngleUnit& unit);

struct EnergyUnit {
  double scaler;
};

void from_json(const json& j, EnergyUnit& unit);

template <size_t N, typename Unit>
std::array<double, N> parseArrayWithUnit(const json& object) {
  const auto unit_scaler = object.at(N).get<Unit>().scaler;
  std::array<double, N> arr;
  for (size_t i = 0; i < N; ++i) {
    arr[i] = unit_scaler*object.at(i).get<double>();
  }
  return arr;
}

template <size_t N, typename T>
std::array<T, N> parseArray(const json& object) {
  std::array<T, N> arr;
  for (size_t i = 0; i < N; ++i) {
    arr[i] = object.at(i).get<T>();
  }
  return arr;
}