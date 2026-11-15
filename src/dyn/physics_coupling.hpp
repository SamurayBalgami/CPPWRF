// SPDX-License-Identifier: GPL-3.0-or-later
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
#include "phys/mp_morrison/morrison.hpp"
#include "phys/mp_p3/p3.hpp"
#include "phys/mp_radar/radar.hpp"
#include "phys/mp_thompson/thompson_tables.hpp"
#include "phys/mp_wsm5/wsm5.hpp"
#include "phys/mp_wsm6/wsm6.hpp"

namespace cppwrf::dyn {

struct PhysicsParams {
  int mp_physics = 0;        // 0 = none, 1 = Kessler, 3 = WSM3, 4 = WSM5, 6 = WSM6, 8 = Thompson, 10 = Morrison,
                             // 28 = aerosol-aware Thompson,
                             // 50 / 51 = P3 one ice category (specified / predicted droplet number)
  int hail_opt = 0;          // WSM6 graupel (0) or hail (1) parameters, WRF namelist hail_opt
  int morr_rimed_ice = 1;    // Morrison: dense precipitating ice is hail (1, WRF's default) or graupel (0)
  real mp_tend_lim = 10;     // K/s, WRF namelist mp_tend_lim (limit on the mp heating rate)
  std::string table_dir;     // Thompson: WRF's qr_acr_qg_V4.dat, qr_acr_qsV2.dat, freezeH2O.dat if there
  bool force_read_thompson = false;  // WRF namelist force_read_thompson: those files must be there
  bool microp_re = false;    // WSM3 / WSM5 / WSM6 / Thompson effective radii for RRTMG (WRF has_reqc = has_reqi = has_reqs = 1:
                             // use_mp_re = 1 with RRTMG LW and SW); P3: has_reqc = has_reqi = 1, has_reqs = 0
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
  bool p3() const { return p_.mp_physics == 50 || p_.mp_physics == 51; }
  // the schemes with snow and graupel accumulations (SNOWNC, GRAUPELNC) and SR: WSM3, WSM5, WSM6, Thompson, Morrison,
  // P3 (WSM3, WSM5, P3: no graupel, GRAUPELNC stays 0 as in WRF)
  bool ice_precip() const {
    return p_.mp_physics == 3 || p_.mp_physics == 4 || p_.mp_physics == 6 || p_.mp_physics == 10 || thompson() || p3();
  }
  // P3's grid%itimestep (incremented per step; the restart's step after a restart)
  long itimestep = 0;
  // Morrison: the cumulus scheme's decoupled tendencies RQRCUTEN, RQSCUTEN, RQICUTEN (phy_prep_part2 has run when
  // the microphysics runs); without cumulus WRF passes its zero arrays
  void set_cumulus_tendencies(const Field3* qr, const Field3* qs, const Field3* qi) {
    cu_qr_ = qr; cu_qs_ = qs; cu_qi_ = qi;
  }
  // WRF grid fields moist_physics_prep/finish_em leave (rho, qv_diabatic, qc_diabatic), by WRF name
  std::vector<const Field3*> fields() const;
  Field2 rainnc, rainncv;  // accumulated / last-step grid-scale precipitation [mm]
  Field2 snownc, snowncv, graupelnc, graupelncv, sr;  // ice_precip(): snow, graupel [mm]; sr = frozen fraction
  Field3 re_cloud, re_ice, re_snow;  // microp_re: effective radii [m] (RE_QC_BG ... before the first call)
  Field2 qnwfa2d, qnifa2d;  // mp_physics = 28: surface emission of water- / ice-friendly aerosols [#/kg/s]
  // do_radar_ref = 1 with WSM6 / Thompson: 10-cm radar reflectivity [dBZ], 0 until the first history step and
  // outside the microphysics columns (Kessler and no microphysics: 0 throughout, as in WRF)
  Field3 refl_10cm;
  // P3: TH_OLD / QV_OLD (the theta and qv after the previous call; restart fields), the diagnostics V_ICE, D_ICE,
  // RHO_ICE (history); REFL_10CM is P3's diag_ze, RE_CLOUD / RE_ICE its effective radii (every step)
  Field3 th_old, qv_old, vmi3d, di3d, rhopo3d;
  // solve_em's diag_flag: the history alarm rings at the end of this step (set before each step)
  bool diagflag = false;

 private:
  Grid g_;
  PhysicsParams p_;
  // moist_physics_prep_em outputs (WRF names)
  Field3 th_phy, rho, pii, z, dz8w, qv_diabatic, qc_diabatic, p_phy;
  phys::wsm6::Wsm6Params<real> wp_{};  // WSM6 constants (wsm6init), computed once
  phys::wsm5::Wsm5Params<real> w5_{};  // WSM5 / WSM3 constants (wsm5init = wsm3init), computed once
#if defined(CPPWRF_DOUBLE)
  Buffer wsm6_work_;  // FP64 builds: the values WSM6's column stages carry between them (Wsm6Work)
#endif
  // Thompson lookup tables (thompson_init), device-resident: one set shared by all domains, as WRF's
  std::shared_ptr<const phys::mp_thompson::ThompsonDeviceTables<real>> tp_;
  std::shared_ptr<const phys::radar::RadarDeviceTables<real>> radar_;  // radar_init of WSM6 / Thompson
  phys::mp_morrison::MorrisonParams<real> morr_{};  // morr_two_moment_init's constants
  phys::mp_p3::P3Params<real> p3c_{};  // p3_init's constants
  std::shared_ptr<const std::pair<Buffer, phys::mp_p3::P3Tables<real>>> p3t_;  // its tables on the device
  const Field3 *cu_qr_ = nullptr, *cu_qs_ = nullptr, *cu_qi_ = nullptr;  // Morrison: RQRCUTEN, RQSCUTEN, RQICUTEN
};

}  // namespace cppwrf::dyn
