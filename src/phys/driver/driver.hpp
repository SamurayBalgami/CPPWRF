#pragma once

// CPPWRF Physics Driver Glue
// Ported statement by statement from WRF:
//   - dyn_em/module_big_step_utilities_em.F (phy_prep)
//   - phys/module_physics_addtendc.F (update_phy_ten, add_a2a, add_a2c_u, add_a2c_v)
//   - dyn_em/module_em.F (calculate_phy_tend, rk_addtend_dry map-factor divisions)
//
// Interfaces:
//   PhyPrepFields, run_phy_prep
//   PhysicsTendencies, calculate_phy_tend, update_phy_ten, apply_tendf_map_factors, couple_physics_tendencies

#include "phys/driver/landuse_tables.hpp"
#include "phys/driver/phy_prep.hpp"
#include "phys/driver/physics_init.hpp"
#include "phys/driver/tendency_coupling.hpp"
