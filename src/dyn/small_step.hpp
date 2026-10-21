#pragma once
// Acoustic small step of the ARW split-explicit scheme (TN 3.2; WRF module_small_step_em.F
// called from solve_em).
//
// Within one RK stage the acoustic variables are perturbations from the stage-start state
// (u'' = mu_t u_t - mu* u*, etc.; WRF small_step_prep). Each substep:
//   K_uv  : advance_uv     horizontal momentum (forward), external-mode filter
//   K_mut : advance_mu_t   column mass, omega, theta (forward-backward) + sumflux
//   K_w   : advance_w      vertically implicit w-phi, then calc_p_rho (+ smdiv) per column
// After the loop, small_step_finish re-couples. calc_p_rho_phi (p, al from the new state) is a
// separate call because WRF runs it after the scalar update (solve_em), with the new moisture.
#include "core/field.hpp"
#include "core/grid.hpp"
#include "dyn/state.hpp"

namespace cppwrf::dyn {

struct SmallStepParams {
  real epssm = 0.1, smdiv = 0.1, emdiv = 0.01;
  bool top_lid = false;  // WRF top_lid; true: rigid lid, w = 0 and phi'' fixed at the model top
  int damp_opt = 0;     // 3: w Rayleigh damping (advance_w), zdamp m deep, dampcoef 1/s
  int hypsometric_opt = 1;  // calc_p_rho_phi: 1 linear, 2 log-linear in p (real data)
  real zdamp = 5000, dampcoef = 0.2;
};

class LateralBC;

class SmallStep {
 public:
  SmallStep(const Grid& g, const SmallStepParams& p);
  // One RK stage (rk_step = 1..3): n_small substeps of length dts.
  void run_stage(State& s, const Tendencies& tend, const BaseState& base, real dt_stage, real dts, int n_small,
                 int rk_step);
  // calc_p_rho_phi (non-hydrostatic, hypsometric_opt = 1): al, p from mu_2, ph_2, t_2 (+ moisture).
  void calc_p_rho_phi(State& s, const BaseState& base);
  void calc_p_rho_phi_interior(State& s, const BaseState& b, int spec_zone);
  void set_lateral_bc(const LateralBC* bc) { bc_ = bc; }

 private:
  Grid g_;
  SmallStepParams p_;
  const LateralBC* bc_ = nullptr;
  int stage_ = 0;
  // stage-start state (WRF *_save), ww1, c2a, pm1, t_2save, and the w-solver coefficients
  Field3 u_save, v_save, w_save, t_save, ph_save, ww1, c2a, pm1, t_2save, a, alpha, gamma;
  Field2 mu_save, mudf, muave;
  Field3 w_rhs, t2a;  // advance_w phases: explicit phi'' (rhs) and t_2ave
  Field3 dvdxi;       // advance_mu_t phase: mass divergence
};

}  // namespace cppwrf::dyn
