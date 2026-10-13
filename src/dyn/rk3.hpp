#pragma once
// Wicker-Skamarock RK3 with split-explicit acoustic substeps (TN 3.1-3.2; WRF solve_em).
//   stage 1: dt/3, 1 acoustic step of dt/3
//   stage 2: dt/2, n/2 acoustic steps of dts = dt/n
//   stage 3: dt,   n   acoustic steps of dts
// Per stage: rk_step_prep, rk_tendency, small steps, small_step_finish, scalar transport
// (rk_scalar_tend + rk_update_scalar), calc_p_rho_phi.
#include <memory>

#include "dyn/large_step.hpp"
#include "dyn/physics_coupling.hpp"
#include "dyn/scalars.hpp"
#include "dyn/state.hpp"

namespace cppwrf::dyn {

struct AcousticParams {
  real epssm = 0.1, smdiv = 0.1, emdiv = 0.01;
  bool top_lid = false;                 // WRF top_lid (default .false.: constant-pressure top, w free)
  int damp_opt = 0;                     // 3: implicit Rayleigh damping of w in the top layer (WRF)
  int hypsometric_opt = 1;              // calc_p_rho_phi: 1 linear (ideal input), 2 log-linear (real input)
  real zdamp = 5000, dampcoef = 0.2;    // depth of the damping layer (m), coefficient (1/s)
};

// WRF set_w_surface: w(surface) = u . grad(ht) from u, v (fill: also w(k) = w(surface) znw(k)^2,
// as start_em does for an initial state without w).
void set_w_surface(Field3& w, Field3& u, Field3& v, const BaseState& b, const Grid& g, bool fill);

class SmallStepHolder;
class LateralBC;  // wraps SmallStep (src/dyn/small_step.hpp) when present

class RK3 {
 public:
  RK3(const Grid& g, const LargeStepParams& lp, const AcousticParams& ap, int time_step_sound,
      const ScalarParams& sp = {});
  ~RK3();
  // WRF solve_em auto rule when time_step_sound == 0 (n even, >= 4).
  static int auto_sound_steps(double dt, double dx, double dy, bool two_d);
  int n_sound() const { return n_sound_; }
  static bool small_step_available();
  // Column physics after the RK3 loop (WRF solve_em order). mp_physics = 0 leaves it off.
  void enable_physics(const PhysicsParams& pp) { phys_ = std::make_unique<PhysicsCoupling>(g_, pp); }
  PhysicsCoupling* physics() const { return phys_.get(); }
  void step(State& s, Tendencies& t, const BaseState& b);
  // Specified lateral boundaries (real data); the RK3 driver does not own the data.
  void set_lateral_bc(LateralBC* bc);
  // Column physics of first_rk_step_part1/2 (radiation, surface, PBL, cumulus); not owned.
  void set_column_physics(ColumnPhysics* cp) { cp_ = cp; large_.set_column_physics(cp); }
  // WRF diffusion (diff_opt 1, km_opt 1/4, sixth order) in RK stage 1; not owned.
  void set_diffusion(DiffusionDriver* d) { dd_ = d; large_.set_diffusion(d); scalars_.set_diffusion(d); }
  LateralBC* lateral_bc() const { return bc_; }

 private:
  Grid g_;
  LargeStep large_;
  ScalarTransport scalars_;
  std::unique_ptr<PhysicsCoupling> phys_;
  int n_sound_;
  std::unique_ptr<SmallStepHolder> small_;
  LateralBC* bc_ = nullptr;
  ColumnPhysics* cp_ = nullptr;
  DiffusionDriver* dd_ = nullptr;
};

}  // namespace cppwrf::dyn
