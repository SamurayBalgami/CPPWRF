#pragma once
// Coupling of the column-physics tendencies to the dynamics (statement ports, WRF v4.6.1):
//   calculate_phy_tend      dyn_em/module_em.F
//   update_phy_ten          phys/module_physics_addtendc.F (phy_ra_ten, phy_bl_ten; add_a2a,
//                           add_a2c_u, add_a2c_v with the specified-boundary ranges)
//   conv_t_tendf_to_moist   dyn_em/module_big_step_utilities_em.F (use_theta_m = 1)
//   phy_prep_part2          dyn_em/module_big_step_utilities_em.F (decoupling with muts)
// Only the radiation and PBL tendencies of the P1 suite are held here (cumulus and shallow
// cumulus coupling are added with those schemes).
#include "core/backend.hpp"
#include "core/field.hpp"
#include "core/grid.hpp"
#include "dyn/state.hpp"

namespace cppwrf::phys::driver {

struct PhysicsTendencies {
  explicit PhysicsTendencies(const Grid& g);
  Field3 rthraten;                                                 // radiation [K/s]
  Field3 rublten, rvblten, rthblten, rqvblten, rqcblten, rqiblten;  // PBL, mass points
};

struct CouplingFlags {
  bool radiation = false;  // ra_lw_physics > 0 .or. ra_sw_physics > 0
  bool pbl = false;        // bl_pbl_physics > 0 (YSU: phy_bl_ten adds th, u, v, qv, qc, qi)
};

void calculate_phy_tend(PhysicsTendencies& pt, const dyn::State& s, const dyn::BaseState& b, const Grid& g,
                        const CouplingFlags& f);
void update_phy_ten(dyn::Tendencies& t, PhysicsTendencies& pt, const dyn::State& s, const Grid& g,
                    const CouplingFlags& f);
void conv_t_tendf_to_moist(dyn::Tendencies& t, const dyn::State& s, const Grid& g);
void phy_prep_part2(PhysicsTendencies& pt, const dyn::State& s, const dyn::BaseState& b, const Grid& g,
                    const CouplingFlags& f);

}  // namespace cppwrf::phys::driver
