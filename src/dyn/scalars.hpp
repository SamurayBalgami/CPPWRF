#pragma once
// Scalar transport (moist species and passive tracers), WRF solve_em scalar section:
// rk_scalar_tend (advect_scalar / advect_scalar_pd with the acoustic-step-averaged mass fluxes
// ru_m, rv_m, ww_m; stage-1 diffusion into *_tendf) followed by rk_update_scalar
// (and rk_update_scalar_pd for moist_adv_opt = 1). TN section 3.1.3 / 4.1.
#include "dyn/large_step.hpp"
#include "dyn/state.hpp"

namespace cppwrf::dyn {

struct ScalarParams {
  int moist_adv_opt = 0;     // 0 = original (advect_scalar), 1 = positive definite (advect_scalar_pd)
  int scalar_adv_opt = 0;    // same for passive tracers (WRF scalar_adv_opt / tracer_adv_opt)
  real khdif = 0, kvdif = 0; // momentum K; scalars use K / prandtl (WRF rk_scalar_tend)
  real pr_inv = 3;           // 1 / prandtl (WRF module_model_constants: prandtl = 1/3)
  int h_adv_order = 5, v_adv_order = 3;
};

class LateralBC;

class DiffusionDriver;

class ScalarTransport {
 public:
  ScalarTransport(const Grid& g, const ScalarParams& p);
  // One RK stage, after small_step_finish. Uses s.ru_m, s.rv_m, s.ww_m, s.muts (end of stage),
  // s.mu_1, s.mu_2, s.alt, s.mut; updates s.q_2 from s.q_1. rk = 1..3, dt_rk = stage length.
  void run_stage(State& s, Tendencies& t, const BaseState& b, int rk, real dt_rk);
  const ScalarParams& params() const { return p_; }
  void set_lateral_bc(const LateralBC* bc) { bc_ = bc; }
  void set_diffusion(const DiffusionDriver* d) { dd_ = d; }

 private:
  Grid g_;
  ScalarParams p_;
  const LateralBC* bc_ = nullptr;
  const DiffusionDriver* dd_ = nullptr;
  Field3 q_new_, scale_;
};

}  // namespace cppwrf::dyn
