#pragma once
// WRF dyn_em/module_big_step_utilities_em.F phy_prep: diagnostic inputs of the column physics
// from the dynamics state (statement port, column kernel; C++ k = WRF k - 1).
#include "core/backend.hpp"
#include "core/field.hpp"
#include "core/grid.hpp"
#include "dyn/state.hpp"

namespace cppwrf::phys::driver {

struct PhyPrepFields {
  explicit PhyPrepFields(const Grid& g);
  // mass levels
  Field3 th_phy;       // theta (dry: theta_m / (1 + rv/rd qv) with use_theta_m = 1)
  Field3 th_phy_m_t0;  // th_phy - t0
  Field3 p_phy;        // p + pb
  Field3 pi_phy;       // (p_phy / p1000mb)**rcp
  Field3 t_phy;        // th_phy * pi_phy (WRF grid%t_phy)
  Field3 rho;          // 1/alt * (1 + qv) (WRF grid%rho)
  Field3 u_phy, v_phy; // winds at mass points (WRF grid%u_phy, grid%v_phy)
  Field3 z;            // height of mass levels (WRF grid%z)
  Field3 p_hyd;        // hydrostatic pressure (WRF grid%p_hyd)
  // w levels (nz + 1)
  Field3 z_at_w;       // (phb + ph) / g (WRF grid%z_at_w)
  Field3 dz8w;         // z_at_w(k+1) - z_at_w(k), 0 at the top level
  Field3 p8w, t8w;     // interface pressure / temperature
  Field3 p_hyd_w;      // hydrostatic pressure at w levels (WRF grid%p_hyd_w)
};

struct PhyPrepParams {
  int use_theta_m = 1;
};

void run_phy_prep(const dyn::State& s, const dyn::BaseState& b, PhyPrepFields& p, const PhyPrepParams& params = {});

}  // namespace cppwrf::phys::driver
