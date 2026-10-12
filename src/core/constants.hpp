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
inline constexpr double r_v = 461.6;           // water vapour gas constant
inline constexpr double rvovrd = r_v / r_d;    // theta_m = theta (1 + rvovrd qv)  (TN 2.2)
inline constexpr double cpv = 4.0 * r_v, cvv = cpv - r_v, cliq = 4190.0;
inline constexpr double xlv = 2.5e6;           // latent heat of vaporisation at 0 C
inline constexpr double ep_2 = r_d / r_v;
inline constexpr double svp1 = 0.6112, svp2 = 17.67, svp3 = 29.65, svpt0 = 273.15;  // Bolton-type e_s
}  // namespace cppwrf::cst
