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

  // TKE (km_opt = 2) and damping (damp_opt = 1, 2)
  real c_k = real(0.15);                  // TKE diffusion coefficient constant (WRF default = 0.15)
  real mix_upper_bound = real(0.1);       // non-dimensional diffusivity upper limit (WRF default = 0.1)
  int bl_pbl_physics = 0;                 // PBL scheme option (0 = none)
  int damp_opt = 0;                       // 0 = none, 1 = upper diffusive layer, 2 = Rayleigh
  real zdamp = real(5000.0);              // damping depth from model top (m)
  real dampcoef = real(0.0);              // damping coefficient
  int isotropic = 0;                      // 0 = anisotropic mixing length, 1 = isotropic
  int isfflx = 1;                         // surface flux option (WRF default = 1)
  int use_theta_m = 0;                    // moist theta flag for tke_buoyancy
  real tke_drag_coefficient = real(0.0);  // surface drag coefficient when isfflx = 0
  real tke_heat_flux = real(0.0);         // surface kinematic heat flux when isfflx = 0 or 2
  real tke_upper_bound = real(1000.0);    // upper bound on TKE (m^2/s^2)
  bool tke_mix2_off = false;              // disable 2nd-order diffusion of TKE
  bool tke_mix6_off = false;              // disable 6th-order diffusion of TKE
  bool moist_mix2_off = false;            // disable 2nd-order diffusion of moist species
  bool nested = false;                    // nested domain flag
  int spec_bdy_width = 5;                 // specified boundary zone width (WRF default = 5)
};

}  // namespace cppwrf::dyn::diffusion
