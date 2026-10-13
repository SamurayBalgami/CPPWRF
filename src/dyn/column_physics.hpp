#pragma once
// Interface between the dycore time step and the column physics suite of a real-data run
// (radiation, surface layer, land surface, PBL, cumulus). The suite lives in src/phys/driver/
//; the dycore only calls it at WRF's points of the time step (solve_em):
//   RK stage 1, after rk_step_prep and init_zero_tendency (LargeStep::zero_tendf):
//     part1  = first_rk_step_part1: phy_prep, radiation_driver, surface_driver, pbl_driver,
//              cumulus_driver (updates the surface state, fills the uncoupled R*TEN fields)
//   RK stage 1, after the stage-1 diffusion into *_tendf, before rk_tendency:
//     part2  = first_rk_step_part2: calculate_phy_tend, update_phy_ten (adds into t.ru_tendf,
//              t.rv_tendf, t.t_tendf, t.q_tendf; the scalar diffusion is added after it)
//   after the RK loop, before microphysics:
//     after_rk = phy_prep_part2 (decouples the persistent R*TEN fields with the new mut)
#include "dyn/state.hpp"

namespace cppwrf::dyn {

class ColumnPhysics {
 public:
  virtual ~ColumnPhysics() = default;
  virtual void part1(const State& s, const BaseState& b, Tendencies& t) = 0;
  virtual void part2(const State& s, const BaseState& b, Tendencies& t) = 0;
  virtual void after_rk(const State& s, const BaseState& b) = 0;
};

}  // namespace cppwrf::dyn
