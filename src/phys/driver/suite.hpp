#pragma once
// P1 physics suite: the column physics of WRF v4.6.1 first_rk_step_part1/part2 and
// phy_prep_part2 for
//   sf_sfclay_physics = 1 (sfclayrev), sf_surface_physics = 2 (Noah, 4 soil layers),
//   bl_pbl_physics = 1 (YSU), mp/cu/ra = 0, isfflx = 1, fractional_seaice = 0, use_theta_m = 1,
//   real data (specified boundaries), sf_urban_physics = 0, no lake model, no SST update.
// Statement-by-statement port of the WRF driver path (call map: docs/PHYSICS_SUITE.md); the
// column schemes are the verified ports in src/phys/{sf_sfclayrev,lsm_noah,bl_ysu}.
//
// Device resident: every field below is allocated once in the constructor (WRF grid% names,
// C++ k = WRF k - 1), the Noah and sfclayrev tables are read/built and uploaded once. Fields
// WRF keeps as locals of solve_em / surface_driver (th_phy, pi_phy, p_phy, dz8w, p8w, t8w,
// psim, psih, gz1oz0, chklowq, QGH, CHS, CHS2, CQS2, CPM, u_phytmp) are scratch members.
// The initial surface state (what wrfinput + phy_init provide) is filled by the caller.
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "core/grid.hpp"
#include "dyn/column_physics.hpp"
#include "dyn/state.hpp"
#include "phys/driver/phy_prep.hpp"
#include "phys/driver/tendency_coupling.hpp"
#include "phys/lsm_noah/noah_tables.hpp"
#include "phys/sf_sfclayrev/sfclayrev.hpp"

namespace cppwrf::phys::driver {

// Integer 2D field (ivgtyp, isltyp, kpbl), mass points, same layout as Field2.
class IntField2 {
 public:
  IntField2() = default;
  IntField2(const Grid& g, std::string name);
  ~IntField2();
  IntField2(IntField2&& o) noexcept;
  IntField2& operator=(IntField2&& o) noexcept;
  IntField2(const IntField2&) = delete;
  IntField2& operator=(const IntField2&) = delete;
  const std::string& name() const { return name_; }
  const Layout& layout() const { return L_; }
  int* data() const { return buf_ + L_.offset; }
  View2<int> view() const { return {data(), int(L_.sj)}; }
  View2<const int> cview() const { return {data(), int(L_.sj)}; }
  void upload(const std::vector<int>& h);   // whole buffer (layout total)
  void download(std::vector<int>& h) const;

 private:
  std::string name_;
  Layout L_;
  int* buf_ = nullptr;
};

// Namelist options of the P1 suite (WRF Registry defaults unless noted).
struct SuiteOptions {
  int sf_sfclay_physics = 1;     // SFCLAYREVSCHEME
  int sf_surface_physics = 2;    // LSMSCHEME (Noah)
  int bl_pbl_physics = 1;        // YSUSCHEME
  int num_soil_layers = 4;
  double bldt = 0;               // minutes; 0 = every step
  int isfflx = 1;
  int ifsnow = 1;
  int fractional_seaice = 0;
  int use_theta_m = 1;
  int ysu_topdown_pblmix = 1;
  int isftcflx = 0, iz0tlnd = 0, shalwater_z0 = 0, scm_force_flux = 0;
  bool rdlai2d = false, usemonalb = false;
  int seaice_albedo_opt = 0;
  double seaice_albedo_default = 0.65;
  int seaice_thickness_opt = 0;
  double seaice_thickness_default = 3.0;
  int seaice_snowdepth_opt = 0;
  double seaice_snowdepth_max = 1.e10, seaice_snowdepth_min = 0.001;
  // land-use / soil classification (global attributes of the input)
  std::string mminlu = "MODIFIED_IGBP_MODIS_NOAH";
  std::string mminsl = "STAS";
  int iswater = 17, isice = 15, isurban = 13;
  std::string table_dir = ".";   // VEGPARM.TBL, SOILPARM.TBL, GENPARM.TBL
};

class PhysicsSuite : public dyn::ColumnPhysics {
 public:
  PhysicsSuite(const Grid& g, const SuiteOptions& o);
  ~PhysicsSuite() override;

  void part1(const dyn::State& s, const dyn::BaseState& b, dyn::Tendencies& t) override;
  void part2(const dyn::State& s, const dyn::BaseState& b, dyn::Tendencies& t) override;
  void after_rk(const dyn::State& s, const dyn::BaseState& b) override;

  long itimestep = 0;  // WRF grid%itimestep (incremented in part1, as solve_em does per step)
  int stepbl = 1;      // WRF grid%stepbl (set_physics_rconfigs: nint(bldt*60/dt), >= 1)

  // ---- persistent 2D surface state (WRF grid% names) ----
  Field2 tsk, xland, xice, lakemask, water_depth, mavail, tmn, vegfra, shdmin, shdmax, snoalb, albbck, albedo,
      emiss, embck, znt, z0, ust, ustm, mol, rmol, pblh, hfx, qfx, lh, grdflx, qsfc, flhc, flqc, wspd, br, fm, fh,
      zol, regime, u10, v10, t2, th2, q2, psfc, ch, snow, snowh, snowc, canwat, lai, snotime, qz0, rainbl, raincv,
      rainncv, rainshv, sr, glw, swdown, gsw, smstav, smstot, sfcrunoff, udrunoff, acsnom, acsnow, snopcx, potevp,
      noahres, sfcevp, sfcexc, achfx, aclhf, acgrdflx, dx2d, uoce, voce, ctopo, ctopo2, ck, cka, cd, cda, wstar,
      delta, hol;
  IntField2 ivgtyp, isltyp, kpbl;
  Field3 tslb, smois, sh2o, smcrel;  // (i, j, soil layer)
  std::vector<double> dzs;            // soil layer thicknesses (m)
  Field3 exch_h, exch_m;              // M (WRF kme-sized; levels 0..nz-1 used)
  PhysicsTendencies pt;               // R*BLTEN (and RTHRATEN) at mass points
  PhyPrepFields prep;                 // phy_prep diagnostics (incl. t_phy, rho, u_phy, ...)

  // ---- locals of solve_em / surface_driver (scratch) ----
  Field2 psim, psih, gz1oz0, chklowq, qgh, chs, chs2, cqs2, cpm;
  Field3 u_phytmp, v_phytmp;

  const SuiteOptions& options() const { return o_; }

 private:
  Grid g_, gs_;
  SuiteOptions o_;
  noah::NoahTables noah_tables_;
  sfclayrev::SFClayRevTables<real>* d_sfclay_tables_ = nullptr;  // device copy, uploaded once
  // Noah driver scratch (inputs/outputs of SFLX per point)
  Field2 n_lwdn, n_soldn, n_solnet, n_sfcprs, n_prcp, n_sfctmp, n_th2, n_q2, n_q2sat, n_dqsdt2, n_ffrozp,
      n_shdfac, n_shmin, n_shmax, n_albbrd, n_snoalb1, n_tbot, n_z0brd, n_emissi, n_embrd, n_cmc, n_t1,
      n_snowhk, n_sneqv, n_chk, n_cm, n_albedok, n_z0k, n_xlai, n_snotime1, n_ribb, n_eta, n_sheat, n_eta_kin,
      n_fdown, n_ec, n_edir, n_ett, n_esnow, n_drip, n_dew, n_beta, n_etp, n_ssoil, n_flx1, n_flx2, n_flx3,
      n_snomlt, n_sncovr, n_runoff1, n_runoff2, n_runoff3, n_rc, n_pc, n_rsmin, n_soilw, n_soilm, n_q1;
  Field3 n_stc, n_smc, n_swc, n_et, n_smav;
  IntField2 n_vegtyp, n_soiltyp, n_slopetyp, n_ice, n_active;
  Field3 zero3_;  // qc/qi stand-ins for YSU (flag_qc = flag_qi = false)

  void surface_driver(const dyn::State& s);
  void pbl_driver(const dyn::State& s);
};

}  // namespace cppwrf::phys::driver
