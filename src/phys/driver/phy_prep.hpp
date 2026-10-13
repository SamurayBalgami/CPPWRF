#pragma once

#include "core/backend.hpp"
#include "core/field.hpp"
#include "core/grid.hpp"
#include "dyn/state.hpp"

namespace cppwrf::phys::driver {

// Diagnostic fields produced by phy_prep (WRF dyn_em/module_big_step_utilities_em.F)
struct PhyPrepFields {
  explicit PhyPrepFields(const Grid& g);

  // Mass-point diagnostics (Stagger::M, nz levels)
  Field3 th_phy;       // dry potential temperature (or theta_m / (1 + rvovrd * qv)) [K]
  Field3 th_phy_m_t0;  // perturbation potential temperature th_phy - t0 [K]
  Field3 p_phy;        // full pressure [Pa]
  Field3 pi_phy;       // Exner function (p_phy / p1000mb)^rcp
  Field3 t_phy;        // full temperature [K]
  Field3 rho;          // air density [kg/m^3]
  Field3 u_phy;        // zonal wind at mass points [m/s]
  Field3 v_phy;        // meridional wind at mass points [m/s]
  Field3 z;            // height at mass levels [m]
  Field3 p_hyd;        // hydrostatic pressure at mass levels [Pa]

  // W-point diagnostics (Stagger::W, nz + 1 levels)
  Field3 z_at_w;       // height at interface/w levels [m]
  Field3 dz8w;         // layer thickness [m] (level nz is 0, levels 0..nz-1 are dz)
  Field3 p8w;          // full pressure at interface/w levels [Pa]
  Field3 t8w;          // temperature at interface/w levels [K]
  Field3 p_hyd_w;      // hydrostatic pressure at interface/w levels [Pa]
};

struct PhyPrepParams {
  int use_theta_m = 1;
};

// Main driver entry point for phy_prep using State and BaseState.
void run_phy_prep(const dyn::State& s, const dyn::BaseState& b,
                  PhyPrepFields& p, const PhyPrepParams& params = {});

}  // namespace cppwrf::phys::driver
