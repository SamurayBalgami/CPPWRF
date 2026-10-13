#pragma once

#include "core/backend.hpp"
#include "core/field.hpp"
#include "core/grid.hpp"
#include "dyn/state.hpp"

namespace cppwrf::phys::driver {

// Container for uncoupled physics tendencies (A-grid, mass points)
struct PhysicsTendencies {
  explicit PhysicsTendencies(const Grid& g);

  // Radiation tendency [K/s]
  Field3 rthraten;

  // Boundary layer (PBL, e.g. YSU) tendencies
  Field3 rthblten;  // potential temperature [K/s]
  Field3 rublten;   // u-wind [m/s^2] (A-grid)
  Field3 rvblten;   // v-wind [m/s^2] (A-grid)
  Field3 rqvblten;  // vapor [kg/kg/s]
  Field3 rqcblten;  // cloud water [kg/kg/s]
  Field3 rqiblten;  // cloud ice [kg/kg/s]

  // Cumulus (CU, e.g. KF) tendencies
  Field3 rthcuten;  // potential temperature [K/s]
  Field3 rucuten;   // u-wind [m/s^2] (A-grid)
  Field3 rvcuten;   // v-wind [m/s^2] (A-grid)
  Field3 rqvcuten;  // vapor [kg/kg/s]
  Field3 rqccuten;  // cloud water [kg/kg/s]
  Field3 rqrcuten;  // rain [kg/kg/s]
  Field3 rqicuten;  // cloud ice [kg/kg/s]
  Field3 rqscuten;  // snow [kg/kg/s]

  // Shallow cumulus (SHCU) tendencies
  Field3 rthshten;  // potential temperature [K/s]
  Field3 rushten;   // u-wind [m/s^2] (A-grid)
  Field3 rvshten;   // v-wind [m/s^2] (A-grid)
  Field3 rqvshten;  // vapor [kg/kg/s]
  Field3 rqcshten;  // cloud water [kg/kg/s]
  Field3 rqrshten;  // rain [kg/kg/s]
  Field3 rqishten;  // cloud ice [kg/kg/s]
  Field3 rqsshten;  // snow [kg/kg/s]
  Field3 rqgshten;  // graupel [kg/kg/s]
};

// 1. A-to-A tendency accumulation (WRF module_physics_addtendc.F: add_a2a)
void add_a2a(Field3& lvar, const Field3& rvar, const Grid& g);
void add_a2a(View3<real> lvar, View3<const real> rvar, Range3D range);

// 2. A-to-C-grid U tendency accumulation (WRF module_physics_addtendc.F: add_a2c_u)
//    lvar(i,j,k) += 0.5 * (rvar(i,j,k) + rvar(i-1,j,k))
void add_a2c_u(Field3& lvar, const Field3& rvar, const Grid& g);
void add_a2c_u(View3<real> lvar, View3<const real> rvar, Range3D range_u);

// 3. A-to-C-grid V tendency accumulation (WRF module_physics_addtendc.F: add_a2c_v)
//    lvar(i,j,k) += 0.5 * (rvar(i,j,k) + rvar(i,j-1,k))
void add_a2c_v(Field3& lvar, const Field3& rvar, const Grid& g);
void add_a2c_v(View3<real> lvar, View3<const real> rvar, Range3D range_v);

// 4. Mu-coupling of physics tendencies (WRF dyn_em/module_em.F: calculate_phy_tend)
//    ten(i,j,k) = (c1h(k) * mut(i,j) + c2h(k)) * ten(i,j,k)
void calculate_phy_tend(PhysicsTendencies& pt, const dyn::State& s, const dyn::BaseState& b, const Grid& g);

// 5. Update model tendencies from mu-coupled physics tendencies (WRF module_physics_addtendc.F: update_phy_ten)
//    Adds mu-coupled tendencies to t_tendf, ru_tendf, rv_tendf, q_tendf using add_a2a, add_a2c_u, add_a2c_v
void update_phy_ten(dyn::Tendencies& t, const PhysicsTendencies& pt, const dyn::State& s, const Grid& g);

// 6. Map-factor coupling into large-step tendencies (WRF dyn_em/module_em.F: rk_addtend_dry)
//    ru_tend += ru_tendf / msfuy
//    rv_tend += rv_tendf * msfvx_inv
//    t_tend  += t_tendf  / msfty
void apply_tendf_map_factors(dyn::Tendencies& t, const dyn::BaseState& b, const Grid& g);

// 7. Combined helper for physics tendency coupling:
//    calculate_phy_tend -> update_phy_ten -> optionally apply_tendf_map_factors
void couple_physics_tendencies(dyn::Tendencies& t, PhysicsTendencies& pt,
                               const dyn::State& s, const dyn::BaseState& b, const Grid& g,
                               bool apply_map_factors = false);

}  // namespace cppwrf::phys::driver
