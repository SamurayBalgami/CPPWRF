#pragma once
// Physical constants from WRF share/module_model_constants.F.
// Do not modify src/core/constants.hpp; physics routines use this header.

namespace cppwrf::phys::cst {

inline constexpr double g = 9.81;                 // gravity (m/s^2)
inline constexpr double r_d = 287.0;              // dry air gas constant (J/kg/K)
inline constexpr double r_v = 461.6;              // water vapor gas constant (J/kg/K)
inline constexpr double cp = 7.0 * r_d / 2.0;     // 1004.5 (J/kg/K)
inline constexpr double cv = cp - r_d;            // 717.5 (J/kg/K)
inline constexpr double cvpm = -cv / cp;
inline constexpr double rcv = r_d / cv;
inline constexpr double rcp = r_d / cp;
inline constexpr double cpovcv = cp / (cp - r_d);
inline constexpr double cvovcp = 1.0 / cpovcv;
inline constexpr double rvovrd = r_v / r_d;
inline constexpr double ep_1 = r_v / r_d - 1.0;
inline constexpr double ep_2 = r_d / r_v;
inline constexpr double p1000mb = 100000.0;       // standard reference pressure (Pa)
inline constexpr double t0 = 300.0;               // reference potential temperature (K)

inline constexpr double xlv = 2.5e6;              // latent heat of vaporization (J/kg)
inline constexpr double xls = 2.85e6;             // latent heat of sublimation (J/kg)
inline constexpr double xlf = 3.50e5;             // latent heat of fusion (J/kg)
inline constexpr double rhowater = 1000.0;        // density of liquid water (kg/m^3)
inline constexpr double rhosnow = 100.0;          // density of snow (kg/m^3)
inline constexpr double rhoair0 = 1.28;           // air density at 0 C, 1000 mb (kg/m^3)

inline constexpr double svp1 = 0.6112;            // constant for saturation vapor pressure
inline constexpr double svp2 = 17.67;             // constant for saturation vapor pressure
inline constexpr double svp3 = 29.65;             // constant for saturation vapor pressure (K)
inline constexpr double svpt0 = 273.15;           // constant for saturation vapor pressure (K)

inline constexpr double karman = 0.4;             // von Karman constant
inline constexpr double stbolt = 5.67051e-8;      // Stefan-Boltzmann constant (W/m^2/K^4)
inline constexpr double cliq = 4190.0;            // specific heat of liquid water (J/kg/K)
inline constexpr double cice = 2106.0;            // specific heat of ice (J/kg/K)

}  // namespace cppwrf::phys::cst
