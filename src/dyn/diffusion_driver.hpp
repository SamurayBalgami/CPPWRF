#pragma once
// WRF diffusion in the time step (diff_opt = 1, km_opt = 1 or 4, diff_6th_opt), built from
// kernels in src/dyn/diffusion/ (differential-tested against the WRF routines).
// Call points (all in RK stage 1, into the held *_tendf, after update_phy_ten):
//   prepare  = first_rk_step_part2: compute_diff_metrics, cal_deform_and_div, calculate_km_kh
//   dynamics = rk_tendency (rk_step 1): horizontal_diffusion u, v, w, horizontal_diffusion_3dmp
//              theta (xkhh), vertical_diffusion_u/_v/'w'/_3dmp (kvdif, only without PBL),
//              sixth_order_diffusion u, v, w, theta (dt)
//   scalar   = rk_scalar_tend (rk_step 1): horizontal_diffusion 'm' (xkhh), vertical_diffusion_mp
//              (qv) / 'm' (3 kvdif, without PBL), sixth_order_diffusion 'm' (dt_rk = dt / 3)
#include "core/field.hpp"
#include "dyn/diffusion/diffusion_params.hpp"
#include "dyn/state.hpp"

namespace cppwrf::dyn {

struct DiffusionConfig {
  int diff_opt = 0;           // 0 off, 1 (diff_opt = 2 not wired yet)
  int km_opt = 1;             // 1 constant khdif/kvdif, 4 2D Smagorinsky
  real khdif = 0, kvdif = 0;  // m2/s (km_opt 1; kvdif also the vertical K without PBL)
  real c_s = real(0.25);      // Smagorinsky constant (km_opt 4)
  int diff_6th_opt = 0;       // 0 off, 1, 2 (flux limited)
  real diff_6th_factor = real(0.12);
  bool pbl = false;           // bl_pbl_physics != 0: no constant-K vertical diffusion
  bool on() const { return diff_opt != 0 || diff_6th_opt != 0; }
};

class DiffusionDriver {
 public:
  DiffusionDriver(const Grid& g, const DiffusionConfig& c);
  const DiffusionConfig& config() const { return c_; }
  void prepare(const State& s, const BaseState& b);
  void dynamics(const State& s, const BaseState& b, Tendencies& t) const;
  void scalar(const State& s, const BaseState& b, int n, Field3& tend, real dt_rk) const;

 private:
  Grid g_, gp_;  // gp_: one more level (fields on w levels at u / v points)
  DiffusionConfig c_;
  diffusion::DiffusionParams p_;
  Field3 z_, rdz_, rdzw_, zx_, zy_;
  Field3 d11_, d22_, d33_, d12_, d13_, d23_, div_;
  Field3 xkmh_, xkmv_, xkhh_, xkhv_;
};

}  // namespace cppwrf::dyn
