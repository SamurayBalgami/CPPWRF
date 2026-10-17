#pragma once
// WRF diffusion in the time step (diff_opt = 1 or 2, km_opt = 1..4, damp_opt = 1 or 2, diff_6th_opt),
// built from the kernels in src/dyn/diffusion/ (differential-tested against the WRF routines).
// Call points (all in RK stage 1, into the held *_tendf, after update_phy_ten):
//   prepare  = first_rk_step_part2: compute_diff_metrics, cal_deform_and_div, calculate_km_kh
//   dynamics = first_rk_step_part2 (diff_opt 2: tke_rhs, vertical_diffusion_2, horizontal_diffusion_2)
//              and rk_tendency (rk_step 1: diff_opt 1, sixth_order_diffusion, rk_rayleigh_damp)
//   scalar   = rk_scalar_tend (rk_step 1): diff_opt 1 horizontal/vertical diffusion, sixth_order_diffusion
//   bound_tke = solve_em after rk_update_scalar on tke_2 (every RK stage)
#include <memory>

#include "core/field.hpp"
#include "dyn/diffusion/diffusion_params.hpp"
#include "dyn/state.hpp"

namespace cppwrf::phys::driver {
class PhysicsSuite;
struct PhyPrepFields;
}  // namespace cppwrf::phys::driver

namespace cppwrf::dyn {

struct DiffusionConfig {
  int diff_opt = 0;           // 0 off, 1 along eta, 2 physical space (full metric)
  int km_opt = 1;             // 1 constant khdif/kvdif, 2 1.5-order TKE, 3 3D Smag, 4 2D Smag
  int damp_opt = 0;           // 0 off, 1 cal_dampkm, 2 rk_rayleigh_damp, 3 implicit w (in small_step)
  real zdamp = real(5000.0);
  real dampcoef = real(0.2);
  real khdif = 0, kvdif = 0;  // m2/s (km_opt 1; kvdif also the vertical K without PBL)
  real c_s = real(0.25);      // Smagorinsky constant (km_opt 2, 3, 4)
  real c_k = real(0.15);      // TKE mixing coefficient (km_opt 2)
  int mix_isotropic = 0;      // 0 anisotropic (dx, dy, dz), 1 isotropic
  bool mix_full_fields = false;
  real mix_upper_bound = real(0.1);
  real prandtl = real(1.0 / 3.0);
  real tke_drag_coefficient = 0;
  real tke_heat_flux = 0;
  real tke_upper_bound = real(1000.0);
  bool moist_mix2_off = false, chem_mix2_off = false, tracer_mix2_off = false, scalar_mix2_off = false, tke_mix2_off = false;
  bool moist_mix6_off = false, chem_mix6_off = false, tracer_mix6_off = false, scalar_mix6_off = false, tke_mix6_off = false;
  int diff_6th_opt = 0;       // 0 off, 1, 2 (flux limited)
  real diff_6th_factor = real(0.12);
  bool pbl = false;           // bl_pbl_physics != 0: no vertical diffusion in dynamics
  int sf_sfclay_physics = 0;
  int isfflx = 1;
  int use_theta_m = 1;
  int fft_filter_lat = 90;
  bool on() const { return diff_opt != 0 || diff_6th_opt != 0 || damp_opt == 2; }
};

class DiffusionDriver {
 public:
  DiffusionDriver(const Grid& g, const DiffusionConfig& c);
  ~DiffusionDriver();
  const DiffusionConfig& config() const { return c_; }
  const diffusion::DiffusionParams& params() const { return p_; }
  void set_physics_suite(phys::driver::PhysicsSuite* suite) { suite_ = suite; }
  void prepare(const State& s, const BaseState& b);
  void dynamics(const State& s, const BaseState& b, Tendencies& t) const;
  void scalar(const State& s, const BaseState& b, int n, Field3& tend, real dt_rk) const;
  void bound_tke(Field3& tke) const;

 private:
  Grid g_, gp_;  // gp_: one more level (fields on w levels at u / v points)
  DiffusionConfig c_;
  diffusion::DiffusionParams p_;
  Field3 z_, rdz_, rdzw_, zx_, zy_;
  Field3 d11_, d22_, d33_, d12_, d13_, d23_, div_;
  Field3 xkmh_, xkmv_, xkhh_, xkhv_, bn2_;
  mutable Field3 titau11_, titau22_, titau33_, titau12_, titau13_, titau23_;
  phys::driver::PhysicsSuite* suite_ = nullptr;
  std::unique_ptr<phys::driver::PhyPrepFields> prep_;
  Field2 ustm_, hfx_, qfx_;
};

}  // namespace cppwrf::dyn
