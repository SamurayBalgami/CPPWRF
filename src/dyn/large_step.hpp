#pragma once
// Large-step (RK stage) part of the ARW dynamical core (TN ch. 3.1-3.2, 4.1; WRF solve_em:
// rk_step_prep, rk_tendency, rk_addtend_dry), dry, flat or terrain-following eta.
#include "core/backend.hpp"
#include "core/field.hpp"
#include "dyn/state.hpp"

namespace cppwrf::dyn {

struct LargeStepParams {
  real khdif = 0, kvdif = 0;  // constant eddy viscosity (m2/s), diff_opt=1-like, km_opt=1
  real pr_inv = 3;            // theta/scalar K = K / prandtl; WRF rk_tendency: khdq = 3*khdif
  bool top_lid = false;       // WRF top_lid (default .false.); true: rigid lid, p' extrapolated to the top
  int h_adv_order = 5;        // 5 (default) or 6 (centred); WRF h_mom/h_sca_adv_order
  int v_adv_order = 3;        // 3 (default) or 4 (centred); WRF v_mom/v_sca_adv_order
  int w_damping = 0;          // WRF w_damping: 1 damps w where the vertical Courant number > 1
  real w_crit_cfl = 1;        // WRF w_crit_cfl (subtracted from the Courant number)
};

// 1D vertical coefficients packed for kernels.
struct VCoef {
  View1<const real> dnw, rdnw, dn, rdn, fnm, fnp, c1h, c2h, c1f, c2f;
  real cf1, cf2, cf3, cfn, cfn1;
  static VCoef from(const BaseState& b);
};

namespace kern { struct MoistCoef; }
class LateralBC;
// calc_cq coefficients for the stage state (needs LargeStep::prep to have filled s.qtot)
kern::MoistCoef moist_coef(const State& s);

class LargeStep {
 public:
  LargeStep(const Grid& g, const LargeStepParams& p) : g_(g), p_(p) {}
  // rk_step_prep: mut/muu/muv from mu_2; ru, rv, rw (couple_momentum); ww (calc_ww_cp);
  // alt = al + alb; php. Fills halos of the stage prognostics and of all derived fields.
  void prep(State& s, const BaseState& b) const;
  // rk_tendency + rk_addtend_dry. Diffusion is evaluated only when rk_step == 1 (*_tendf).
  void tendencies(const State& s, Tendencies& t, const BaseState& b, int rk_step) const;
  const LargeStepParams& params() const { return p_; }
  void set_lateral_bc(const LateralBC* bc) { bc_ = bc; }

 private:
  Grid g_;
  LargeStepParams p_;
  const LateralBC* bc_ = nullptr;
};

}  // namespace cppwrf::dyn
