#pragma once
// share/module_model_constants.F evaluated in the working precision, as WRF evaluates its
// PARAMETERs in the default REAL kind (real*8 with -r8): cp = 7.*r_d/2., rcp = r_d/cp, ...
#include "core/precision.hpp"

namespace cppwrf::phys::wrfc {

inline constexpr real g = real(9.81);
inline constexpr real r_d = real(287.);
inline constexpr real cp = real(7.) * r_d / real(2.);
inline constexpr real r_v = real(461.6);
inline constexpr real rcp = r_d / cp;
inline constexpr real rovg = r_d / g;
inline constexpr real p1000mb = real(100000.);
inline constexpr real t0 = real(300.);
inline constexpr real xlv = real(2.5e6);
inline constexpr real xlf = real(3.50e5);
inline constexpr real rhowater = real(1000.);
inline constexpr real svp1 = real(0.6112), svp2 = real(17.67), svp3 = real(29.65), svpt0 = real(273.15);
inline constexpr real ep_1 = r_v / r_d - real(1.);
inline constexpr real ep_2 = r_d / r_v;
inline constexpr real karman = real(0.4);
inline constexpr real stbolt = real(5.67051e-8);

}  // namespace cppwrf::phys::wrfc
