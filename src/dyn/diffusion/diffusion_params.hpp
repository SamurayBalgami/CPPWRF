#pragma once

#include "core/precision.hpp"

namespace cppwrf::dyn::diffusion {

struct DiffusionParams {
  int diff_opt = 0;              // 0 = none, 1 = simple along eta, 2 = full stress tensor in physical space
  int km_opt = 1;                // 1 = constant/isotropic, 2 = TKE, 3 = 3D Smagorinsky, 4 = 2D Smagorinsky
  real khdif = 0;                // horizontal diffusion constant (m^2/s)
  real kvdif = 0;                // vertical diffusion constant (m^2/s)
  real c_s = real(0.25);         // Smagorinsky constant
  real prandtl = real(1.0 / 3.0);// turbulent Prandtl number (WRF default = 1/3)
  real pr_inv = real(3.0);       // 1 / prandtl

  // 6th order diffusion
  int diff_6th_opt = 0;          // 0 = none, 1 = linear, 2 = monotonic flux-limited
  real diff_6th_factor = 0;      // 6th order diffusion non-dimensional factor
  int diff_6th_slopeopt = 0;     // 0 = no slope damping

  // Vertical velocity damping
  int w_damping = 0;             // 0 = none, 1 = damp excessive w
  real w_crit_cfl = real(2.0);   // critical vertical CFL for damping
  real w_alpha = real(0.3);      // damping strength (m/s^2)
  real w_beta = real(1.0);       // activation CFL number

  // Boundary conditions
  bool periodic_x = false;
  bool periodic_y = false;
  bool specified = false;        // real-data lateral boundary conditions
  bool mix_full_fields = false;  // mix perturbations vs full fields (WRF default = false)
  bool open_xs = false;
  bool open_xe = false;
  bool open_ys = false;
  bool open_ye = false;
  bool polar = false;
  real fft_filter_lat = 0;
  int zadvect_implicit = 0;
};

}  // namespace cppwrf::dyn::diffusion
