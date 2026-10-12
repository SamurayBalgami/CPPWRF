#pragma once
// Coupling of column physics to the dycore, as WRF solve_em does after the RK3 loop:
//   moist_physics_prep_em  -> microphysics_driver -> moist_physics_finish_em -> calc_p_rho_phi
// The microphysics latent heating is kept in State::h_diabatic. During the next step it
// enters the theta tendency in every RK stage (rk_addtend_dry), and it is removed again in
// the last small_step_finish, so the dynamics feel the heating but theta is not changed twice.
// Interface for physics schemes (src/phys/): docs/PHYSICS_INTERFACE.md.
#include "core/field.hpp"
#include "dyn/state.hpp"

namespace cppwrf::dyn {

struct PhysicsParams {
  int mp_physics = 0;        // 0 = none, 1 = Kessler (src/phys/mp_kessler)
  real mp_tend_lim = 10;     // K/s, WRF namelist mp_tend_lim (limit on the mp heating rate)
};

class PhysicsCoupling {
 public:
  PhysicsCoupling(const Grid& g, const PhysicsParams& p);
  // One call per time step, after the RK3 loop (before the final calc_p_rho_phi).
  void step(State& s, const BaseState& b, real dt);
  const PhysicsParams& params() const { return p_; }
  Field2 rainnc, rainncv;  // accumulated / last-step grid-scale precipitation [mm]

 private:
  Grid g_;
  PhysicsParams p_;
  // moist_physics_prep_em outputs (WRF names)
  Field3 th_phy, rho, pii, z, dz8w, qv_diabatic;
};

}  // namespace cppwrf::dyn
