#pragma once
// Physical constants, values as in WRF share/module_model_constants.F.
namespace cppwrf::cst {
inline constexpr double g = 9.81;
inline constexpr double r_d = 287.0;
inline constexpr double cp = 7.0 * r_d / 2.0;
inline constexpr double cv = cp - r_d;
inline constexpr double cvpm = -cv / cp;
inline constexpr double rcv = r_d / cv;
inline constexpr double rcp = r_d / cp;
inline constexpr double cpovcv = cp / (cp - r_d);
inline constexpr double p1000mb = 100000.0;
inline constexpr double t0 = 300.0;  // reference theta: t_* = theta - t0
}  // namespace cppwrf::cst
