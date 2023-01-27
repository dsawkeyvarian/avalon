#pragma once
#include "external/span.hpp"
#include <fmt/core.h>
#include <vector>
#include <algorithm>

//Assumes that xs is sorted in ascending order
double interpolate(
  tcb::span<const double> xs,
  tcb::span<const double> ys,
  double x);
//Does not do anything to elements where out_xs[i] is out-of-bounds w.r.t. in_xs
//Assumes that curves in_xs and out_xs are sorted in ascending order
void interpolateCurve(
  tcb::span<const double> in_xs,
  tcb::span<const double> in_ys,
  tcb::span<const double> out_xs,
  tcb::span<double> out_ys);

//Assumes that curve is sorted in ascending order w.r.t pdf_xs
//CDF - Cumulative Distribution Function
//PDF - Probability Distribution Function
void calculateCDF(
  tcb::span<const double> pdf_xs,
  tcb::span<const double> pdf_ys,
  tcb::span<double> cdf);

std::vector<double> calculateCDF(
  tcb::span<const double> pdf_xs,
  tcb::span<const double> pdf_ys);



