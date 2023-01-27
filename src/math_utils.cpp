#include "math_utils.h"

double interpolate(
  tcb::span<const double> xs,
  tcb::span<const double> ys,
  double x) 
{
  if((x < xs[0]) || x > xs.back())
    throw std::runtime_error(fmt::format("Interpolation error: Value {} out-of-bounds [{}, {}]", x, xs[0], xs.back()));
  auto it = std::lower_bound(xs.begin(), xs.end(), x);
  double x0 = *(it - 1);
  double x1 = *it;
  size_t idx = it - xs.begin();
  double y0 = ys[idx - 1];
  double y1 = ys[idx];
  double t = (x - x0)/(x1 - x0);
  return (1-t)*y0 + t*y1;
}

void interpolateCurve(
  tcb::span<const double> in_xs,
  tcb::span<const double> in_ys,
  tcb::span<const double> out_xs,
  tcb::span<double> out_ys)
{
  size_t in_idx = 1;
  for(size_t i = 0; i < out_xs.size(); ++i) {
    double x = out_xs[i];
    for(; in_idx < in_xs.size(); ++in_idx) {
      if((in_xs[in_idx-1] <= x) && (x < in_xs[in_idx]))
        break;
    }
    if(in_idx == in_xs.size()) break;
    double x0 = in_xs[in_idx - 1];
    double x1 = in_xs[in_idx];
    double t = (x - x0)/(x1 - x0);
    out_ys[i] =  (1- t)*in_ys[in_idx-1] + t*in_ys[in_idx];
  }
}

void calculateCDF(
  tcb::span<const double> pdf_xs,
  tcb::span<const double> pdf_ys,
  tcb::span<double> cdf)
{
  //Use trapezoidal integration
  cdf[0] = 0.0;
  double integral = 0.0;
  const size_t num = pdf_ys.size();
  for(size_t i = 1; i < num; ++i) {
    double weight = 0.5*(pdf_xs[i] - pdf_xs[i-1]);
    integral += weight*(pdf_ys[i] + pdf_ys[i-1]);
    cdf[i] = integral;
  }
  const double I = cdf.back();
  for(size_t i = 0; i < num; ++i)
    cdf[i] = cdf[i] / I;
}

std::vector<double> calculateCDF(
  tcb::span<const double> pdf_xs,
  tcb::span<const double> pdf_ys)
{
  std::vector<double> cdf(pdf_xs.size());
  calculateCDF(pdf_xs, pdf_ys, tcb::span(cdf.data(), cdf.size()));
  return cdf;
}