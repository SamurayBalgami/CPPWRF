#pragma once
// Coupling of column physics to the dycore, as WRF solve_em does after the RK3 loop:
//   moist_physics_prep_em  -> microphysics_driver -> moist_physics_finish_em -> calc_p_rho_phi
// The microphysics latent heating is kept in State::h_diabatic. During the next step it
// enters the theta tendency in every RK stage (rk_addtend_dry), and it is removed again in
// the last small_step_finish, so the dynamics feel the heating but theta is not changed twice.
// Interface for physics schemes (src/phys/): docs/PHYSICS_INTERFACE.md.
#include <memory>
#include <string>

#include "core/field.hpp"
#include "dyn/state.hpp"
#include "phys/mp_radar/radar.hpp"
#include "phys/mp_thompson/thompson_tables.hpp"
#include "phys/mp_wsm6/wsm6.hpp"

namespace cppwrf::dyn {

struct PhysicsParams {
  int mp_physics = 0;        // 0 = none, 1 = Kessler, 6 = WSM6, 8 = Thompson, 28 = aerosol-aware Thompson
  int hail_opt = 0;          // WSM6 graupel (0) or hail (1) parameters, WRF namelist hail_opt
  real mp_tend_lim = 10;     // K/s, WRF namelist mp_tend_lim (limit on the mp heating rate)
  std::string table_dir;     // Thompson: WRF's qr_acr_qg_V4.dat, qr_acr_qsV2.dat, freezeH2O.dat if there
  bool force_read_thompson = false;  // WRF namelist force_read_thompson: those files must be there
  bool microp_re = false;    // WSM6 / Thompson effective radii for RRTMG (WRF has_reqc = has_reqi = has_reqs = 1:
                             // use_mp_re = 1 with RRTMG LW and SW)
  int do_radar_ref = 0;      // WRF namelist do_radar_ref: REFL_10CM from WSM6 / Thompson at history steps
};

class PhysicsCoupling {
 public:
  PhysicsCoupling(const Grid& g, const PhysicsParams& p);
  // One call per time step, after the RK3 loop (before the final calc_p_rho_phi).
  void step(State& s, const BaseState& b, real dt);
  // mp_physics = 28: thompson_init's aerosol part (start_domain, also at a restart): without initial CCN
  // (MAXVAL(nwfa) < eps) the boundary-layer following nwfa profile and its emission qnwfa2d, likewise nifa
  void init_aerosols(State& s, const BaseState& b);
  const PhysicsParams& params() const { return p_; }
  bool thompson() const { return p_.mp_physics == 8 || p_.mp_physics == 28; }
  // WRF grid fields moist_physics_prep/finish_em leave (rho, qv_diabatic, qc_diabatic), by WRF name
  std::vector<const Field3*> fields() const;
  Field2 rainnc, rainncv;  // accumulated / last-step grid-scale precipitation [mm]
  Field2 snownc, snowncv, graupelnc, graupelncv, sr;  // WSM6/Thompson: snow, graupel [mm]; sr = frozen fraction
  Field3 re_cloud, re_ice, re_snow;  // microp_re: effective radii [m] (RE_QC_BG ... before the first call)
  Field2 qnwfa2d, qnifa2d;  // mp_physics = 28: surface emission of water- / ice-friendly aerosols [#/kg/s]
  // do_radar_ref = 1 with WSM6 / Thompson: 10-cm radar reflectivity [dBZ], 0 until the first history step and
  // outside the microphysics columns (Kessler and no microphysics: 0 throughout, as in WRF)
  Field3 refl_10cm;
  // solve_em's diag_flag: the history alarm rings at the end of this step (set before each step)
  bool diagflag = false;

 private:
  Grid g_;
  PhysicsParams p_;
  // moist_physics_prep_em outputs (WRF names)
  Field3 th_phy, rho, pii, z, dz8w, qv_diabatic, qc_diabatic, p_phy;
  phys::wsm6::Wsm6Params<real> wp_{};  // WSM6 constants (wsm6init), computed once
  // Thompson lookup tables (thompson_init), device-resident: one set shared by all domains, as WRF's
  std::shared_ptr<const phys::mp_thompson::ThompsonDeviceTables<real>> tp_;
  std::shared_ptr<const phys::radar::RadarDeviceTables<real>> radar_;  // radar_init of WSM6 / Thompson
};

}  // namespace cppwrf::dyn
