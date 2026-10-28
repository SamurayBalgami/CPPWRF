#pragma once
// P1 physics suite: the column physics of WRF v4.8.0 first_rk_step_part1/part2 and
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
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "core/grid.hpp"
#include "dyn/column_physics.hpp"
#include "dyn/state.hpp"
#include "phys/bl_mynn/mynnedmf_common.hpp"
#include "phys/cu_kf/kf.hpp"
#include "phys/cu_kf/kf_tables.hpp"
#include "phys/driver/phy_prep.hpp"
#include "phys/driver/radiation.hpp"
#include "phys/driver/tendency_coupling.hpp"
#include "phys/lsm_noah/noah_tables.hpp"
#include "phys/lsm_noahmp/noahmp.hpp"
#include "phys/lsm_noahmp/noahmp_tables.hpp"
#include "phys/sf_myjsfc/myjsfc.hpp"
#include "phys/sf_mynn/mynnsfc_common.hpp"
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

// Namelist options of the P1/P2 suite (WRF Registry defaults unless noted).
struct SuiteOptions {
  int sf_sfclay_physics = 1;     // 1 = SFCLAYREVSCHEME, 2 = MYJSFCSCHEME, 5 = MYNNSFCSCHEME
  int sf_surface_physics = 2;    // 2 = LSMSCHEME (Noah), 4 = NOAHMPSCHEME (Noah-MP v5.1.2)
  int bl_pbl_physics = 1;        // 1 = YSUSCHEME, 2 = MYJPBLSCHEME, 5 = MYNNPBLSCHEME
  int icloud_bl = 1;             // MYNN subgrid clouds in the radiation (CLDFRA_BL, QC_BL, QI_BL)
  mynnedmf::Options mynn;        // bl_mynn_* options of the MYNN-EDMF PBL
  int bl_mynn_mixscalars = 0;    // MYNN mixing of WRF's generic scalar (not connected in WRF: no effect)
  int bl_mynn_output = 0;        // 1: EDMF_A, EDMF_W, EDMF_THL, EDMF_QT, EDMF_ENT, EDMF_QC (and zero SUB_/DET_*3D)
  int tke_budget = 0;            // 1: the MYNN TKE budget (QSHEAR, QBUOY, QDISS, QWT, DQKE)
  int cu_physics = 0;            // 0 = none, 1 = KFETASCHEME, 3 = GFLSCHEME, 6 = TIEDTKESCHEME, 16 = NTIEDTKESCHEME
  int gwd_opt = 0;                // 1 = KIM orographic gravity wave drag (bl_gwdo)
  double gwd_dx_factor = 2.0;     // effective grid factor of the gwd stress (WRF gwd_dx_factor)
  bool gwd_if_nonhyd = true;      // nonhydrostatic effect (WRF gwd_if_nonhyd)
  int ishallow = 0;              // Grell-Freitas-Li shallow convection (WRF ishallow)
  int gf_ichoice = 0;            // Grell-Freitas-Li closure (WRF clos_choice)
  int num_soil_layers = 4;
  double bldt = 0;               // minutes; 0 = every step
  double cudt = 0;               // minutes; 0 = every step
  double convtrans_avglen_m = 30.0;  // minutes (module_convtrans_prep.F averaging window)
  int kfeta_trigger = 1;         // 1 = default KF trigger, 3 = RH perturbation
  bool warm_rain = false;        // mp_init: true for mp_physics == 1 (Kessler)
  RadiationOptions rad;          // rad.ra_lw_physics/ra_sw_physics = 0: no radiation
  noahmp::NoahmpOptions noahmp;  // Noah-MP options (sf_surface_physics == 4)
  int isfflx = 1;
  int fractional_seaice = 0;
  bool tice2tsk_if2cold = false;  // fractional_seaice: ice temperature = min(TSK, 273.15)
  int sst_update = 0;             // 1: SST, SEAICE, VEGFRA, ALBBCK, LAI from wrflowinp (read_lowinp)
  int use_theta_m = 1;
  int ysu_topdown_pblmix = 1;
  int isftcflx = 0, iz0tlnd = 0, shalwater_z0 = 0, scm_force_flux = 0;
  bool kim_tofd = false;         // sfclayrev turbulent orographic form drag (WRF 4.8)
  double tofd_factor = 0.003;
  bool rdlai2d = false, usemonalb = false;
  int opt_thcnd = 1;  // Noah soil thermal conductivity (1 or 2)
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
  std::string table_dir = ".";   // VEGPARM.TBL, SOILPARM.TBL, GENPARM.TBL, NoahmpTable.TBL
};

class PhysicsSuite : public dyn::ColumnPhysics {
 public:
  PhysicsSuite(const Grid& g, const SuiteOptions& o);
  ~PhysicsSuite() override;

  void part1(dyn::State& s, const dyn::BaseState& b, dyn::Tendencies& t) override;
  void part2(const dyn::State& s, const dyn::BaseState& b, dyn::Tendencies& t) override;
  void after_rk(const dyn::State& s, const dyn::BaseState& b) override;
  Field3* th_forcing() override { return o_.cu_physics == 16 || o_.cu_physics == 3 ? &pt.rthften : nullptr; }
  // rk_scalar_tend sets RQVFTEN for Kain-Fritsch too (module_em.F), though KF does not use it
  Field3* qv_forcing() override {
    return o_.cu_physics == 16 || o_.cu_physics == 3 || o_.cu_physics == 1 || o_.cu_physics == 6 ? &pt.rqvften : nullptr;
  }

  long itimestep = 0;  // WRF grid%itimestep (incremented in part1, as solve_em does per step)
  // WRF solve_em's restart_flag: the restart alarm rings at the end of this step (Is_alarm_tstep); the MYNN
  // surface layer skips its z/L iteration at itimestep 1 when it is set (module_sf_mynnsfc_land.F90 etc.)
  bool restart_flag = false;
  int stepbl = 1;      // WRF grid%stepbl (set_physics_rconfigs: nint(bldt*60/dt), >= 1)
  int stepcu = 1;      // WRF grid%stepcu (cu_init: max(nint(cudt*60/dt), 1))
  // true (model): keep WRF's grid%t_1 for conv_t_tendf_to_moist, i.e. theta_m at the start of the
  // previous step (the dycore's s.t_1 is the current step start); false: s.t_1 already holds it
  // (the differential test passes WRF's t_1).
  bool wrf_t1_lag = true;

  // ---- persistent 2D surface and cumulus state (WRF grid% names) ----
  Field2 tsk, xland, xice, lakemask, water_depth, mavail, tmn, vegfra, shdmin, shdmax, snoalb, albbck, albedo,
      emiss, embck, znt, z0, ust, ustm, mol, rmol, pblh, hfx, qfx, lh, grdflx, qsfc, flhc, flqc, wspd, br, fm, fh,
      zol, regime, u10, v10, t2, th2, q2, psfc, ch, snow, snowh, snowc, canwat, lai, snotime, qz0, rainbl, raincv,
      raincv_a, raincv_b, rainncv, rainshv, rainc, rainsh, pratec, pratesh, nca, htop, hbot, cutop, cubot, cuppt, sr,
      glw, swdown, gsw,
      smstav, smstot, sfcrunoff, udrunoff, acsnom, acsnow, snopcx, potevp, noahres, sfcevp, sfcexc, achfx, aclhf,
      acgrdflx, dx2d, var2d, uoce, voce, ctopo, ctopo2, ck, cka, cd, cda, wstar, delta, hol,
      thz0, uz0, vz0, akhs, akms, ct, tshltr, th10, qshltr, q10, pshltr, u10e, v10e, mixht;
  IntField2 ivgtyp, isltyp, kpbl, lowlyr;
  IntField2 td_icum, td_row_icum;  // Tiedtke: the column's and the row's icum after the first ascent
  Field3 tslb, smois, sh2o, smcrel;  // (i, j, soil layer)
  // ---- Noah-MP persistent state and diagnostics (allocated when sf_surface_physics == 4) ----
  Field2 xlat_sfc, xlong_sfc, tvxy, tgxy, canicexy, canliqxy, eahxy, tahxy, cmxy, chxy, fwetxy, sneqvoxy, alboldxy,
      qsnowxy, qrainxy, wslakexy, zwtxy, waxy, wtxy, lfmassxy, rtmassxy, stmassxy, woodxy, stblcpxy, fastcpxy,
      xsaixy, taussxy, smcwtdxy, deeprechxy, rechxy, grainxy, gddxy, qtdrain, t2mvxy, t2mbxy, q2mvxy, q2mbxy,
      tradxy, neexy, gppxy, nppxy, fvegxy, runsfxy, runsbxy, ecanxy, edirxy, etranxy, fsaxy, firaxy, aparxy, psnxy,
      savxy, sagxy, rssunxy, rsshaxy, bgapxy, wgapxy, tgvxy, tgbxy, chvxy, chbxy, shgxy, shcxy, shbxy, evgxy, evbxy,
      ghvxy, ghbxy, irgxy, ircxy, irbxy, trxy, evcxy, chleafxy, chucxy, chv2xy, chb2xy, chstarxy, acc_ssoilxy,
      acc_qinsurxy, acc_qsevaxy, eflxbxy, acc_dwaterxy, acc_prcpxy, acc_ecanxy, acc_etranxy, acc_edirxy, acc_glaflwxy;
  IntField2 isnowxy, cropcat, pgsxy;
  IntField2 noahmp_ranland;  // not in the restart file (NoahmpIO is initialised again)
  Field3 tsnoxy, snicexy, snliqxy, zsnsoxy, albsoildirxy, albsoildifxy, smoiseq, acc_etranixy;
  std::vector<double> dzs;            // soil layer thicknesses (m)
  Field3 exch_h, exch_m;              // M (WRF kme-sized; levels 0..nz-1 used)
  Field3 tke_pbl, el_pbl;             // MYJ PBL 3D TKE and master length scale (el_pbl: MYNN too)
  // MYNN-EDMF state (bl_pbl_physics = 5): 2*TKE, variances, stability functions, subgrid clouds
  Field3 qke, tsq, qsq, cov, sh3d, sm3d, qc_bl, qi_bl, cldfra_bl;
  Field3 edmf_a, edmf_w, edmf_thl, edmf_qt, edmf_ent, edmf_qc;  // bl_mynn_output = 1
  Field3 qshear, qbuoy, qdiss, qwt, dqke;                        // tke_budget = 1
  Field2 cqs, maxwidth, maxmf, ztop_plume, excess_h, excess_q;
  // gwd_opt = 1: subgrid orography statistics (wrfinput CON, OA1-4, OL1-4, ELVMAX) and the drag diagnostics
  Field2 oc12d, oa1, oa2, oa3, oa4, ol1, ol2, ol3, ol4, elvmax, dusfcg, dvsfcg;
  Field3 dtaux3d, dtauy3d;
  // fractional_seaice = 1 or sst_update = 1: SST (the surface layer wrappers raise it to 271.4 under ice)
  // and TSK_SAVE; sst_update = 1: XICEM, the sea ice of the previous surface_driver call
  Field2 sst, tsk_save, xicem;
  // the open-water halves of the sea-ice points (module_surface_driver.F *_SEA and *_HOLD locals)
  struct SeaIceFields {
    explicit SeaIceFields(const Grid& g);
    Field2 tsk_local, tsk_sea, xland_sea, mavail_sea, znt_sea, z0base_sea, br_sea, chs_sea, chs2_sea, cqs2_sea,
        cpm_sea, flhc_sea, flqc_sea, gz1oz0_sea, hfx_sea, qfx_sea, lh_sea, mol_sea, psim_sea, psih_sea, fm_sea,
        fh_sea, qgh_sea, qsfc_sea, regime_sea, rmol_sea, ust_sea, ustm_sea, wspd_sea, zol_sea, u10_sea, v10_sea,
        th2_sea, t2_sea, q2_sea, ck_sea, cka_sea, cd_sea, cda_sea, qz0_sea, thz0_sea, uz0_sea, vz0_sea, pblh_sea,
        akhs_sea, akms_sea, ct_sea, tshltr_sea, th10_sea, qshltr_sea, q10_sea, pshltr_sea, u10e_sea, v10e_sea,
        regime_hold;
  };
  std::unique_ptr<SeaIceFields> sea_;
  Field3 w0avg, cldfra_dp, cldfra_sh, qc_cu, qi_cu, gd_cloud, gd_cloud2, gd_cloud_a, gd_cloud2_a, gd_cldfr;  // cumulus 3D state and diagnostics
  // every persistent field above, the physics tendencies and the phy_prep diagnostics (by WRF name)
  std::vector<const Field2*> fields2() const;
  std::vector<const Field3*> fields3() const;
  PhysicsTendencies pt;               // R*BLTEN, R*CUTEN, RTHRATEN at mass points
  PhyPrepFields prep;                 // phy_prep diagnostics (incl. t_phy, rho, u_phy, ...)

  // ---- locals of solve_em / surface_driver (scratch) ----
  Field2 psim, psih, gz1oz0, chklowq, qgh, chs, chs2, cqs2, cpm;
  Field3 u_phytmp, v_phytmp;
  Field3 t1_lag_;  // WRF grid%t_1 as seen by first_rk_step_part2 (wrf_t1_lag)
  Field2 coszen;   // radiation: cosine of the solar zenith angle (middle of the radiation interval)
  // radiation diagnostics (radiation_driver): LW/SW heating rates, fluxes and their accumulations (every step)
  Field3 rthratenlw, rthratensw, rthratenlwc, rthratenswc, o3rad;
  Field2 lwcf, swcf, hrang;
  // surface_driver diagnostics: uratx, vratx, tratx (obs-nudging ratios), ra (Noah), rs (Noah canopy
  // resistance), shall (KF shallow flag)
  Field2 uratx, vratx, tratx, ra, rs, shall;
  Field2 lwupt, lwuptc, lwdnt, lwdntc, lwupb, lwupbc, lwdnb, lwdnbc, olr, swupt, swuptc, swdnt, swdntc, swupb, swupbc, swdnb, swdnbc;
  // RRTMG_SWRAD's surface components (restart fields): visible / near-IR direct and diffuse, direct, direct normal,
  // diffuse, clear-sky direct and direct normal irradiance; the diffuse fraction (radiation_driver)
  Field2 swvisdir, swvisdif, swnirdir, swnirdif, swddir, swddni, swddif, swddirc, swddnic, diffuse_frac;
  Field2 aclwupt, aclwuptc, aclwdnt, aclwdntc, aclwupb, aclwupbc, aclwdnb, aclwdnbc, acswupt, acswuptc, acswdnt, acswdntc, acswupb, acswupbc, acswdnb, acswdnbc;
  std::unique_ptr<RadiationDriver> rad_;
  Field3 *re_cloud_ = nullptr, *re_ice_ = nullptr, *re_snow_ = nullptr;
  // radiation needs latitude/longitude (wrfinput XLAT/XLONG, i + nx*j); call before the first step
  void init_radiation(const std::vector<double>& xlat, const std::vector<double>& xlong);
  // Noah-MP initialization (NoahmpWRFinit -> NoahmpInitMain + NoahmpSnowInitMain)
  void init_noahmp(bool fndsnowh);
  // WSM6 effective radii (PhysicsCoupling, has_reqc = has_reqi = has_reqs = 1) seen by the radiation
  void set_effective_radii(Field3* re_cloud, Field3* re_ice, Field3* re_snow) {
    re_cloud_ = re_cloud; re_ice_ = re_ice; re_snow_ = re_snow;
  }
  // grid-scale precipitation of the last step (PhysicsCoupling): WRF's grid%rainncv and grid%sr, read by
  // surface_driver (RAINBL) and the Noah LSM (FFROZP = SR); sr = nullptr for schemes without it (WRF: 0)
  void set_microphysics_precip(const Field2* rainncv_mp, const Field2* sr_mp,
                               const Field2* snowncv_mp = nullptr, const Field2* graupelncv_mp = nullptr) {
    rainncv_mp_ = rainncv_mp; sr_mp_ = sr_mp;
    snowncv_mp_ = snowncv_mp; graupelncv_mp_ = graupelncv_mp;
  }
  const Field2 *rainncv_mp_ = nullptr, *sr_mp_ = nullptr, *snowncv_mp_ = nullptr, *graupelncv_mp_ = nullptr;
  Field3 cldfra;  // radiation cloud fraction (WRF CLDFRA)
  const RadiationDriver* radiation() const { return rad_.get(); }
  bool radiation_on() const { return o_.rad.ra_lw_physics != 0 || o_.rad.ra_sw_physics != 0; }
  bool t1_lag_set_ = false;

  const SuiteOptions& options() const { return o_; }

  // State WRF carries between steps (written to and read from restart files): the persistent
  // grid% fields above, R*BLTEN/RTHRATEN, EXCH_H/M, the soil state and the lagged t_1.
  std::vector<std::pair<std::string, Field2*>> carried2();
  std::vector<std::pair<std::string, Field3*>> carried3();
  std::vector<std::pair<std::string, IntField2*>> carried_int();
  void restarted(long step) {
    itimestep = step;
    t1_lag_set_ = true;
  }

 private:
  Grid g_, gs_;
  SuiteOptions o_;
  noah::NoahTables noah_tables_;
  noahmp::NoahmpTables noahmp_tables_;
  noahmp::NoahmpTables* d_noahmp_tables_ = nullptr;              // device copy, uploaded once
  sfclayrev::SFClayRevTables<real>* d_sfclay_tables_ = nullptr;  // device copy, uploaded once
  sf_myjsfc::MYJSFCTableManager myjsfc_tables_;                  // host/device copy, uploaded once
  cu_kf::KFTables<real>* d_kf_tables_ = nullptr;                 // device copy, uploaded once
  mynnsfc::PsiTables* d_mynn_psi_ = nullptr;                     // MYNN psi_init tables, device copy
  // Noah driver scratch (inputs/outputs of SFLX per point)
  Field2 n_lwdn, n_soldn, n_solnet, n_sfcprs, n_prcp, n_sfctmp, n_th2, n_q2, n_q2sat, n_dqsdt2, n_ffrozp,
      n_shdfac, n_shmin, n_shmax, n_albbrd, n_snoalb1, n_tbot, n_z0brd, n_emissi, n_embrd, n_cmc, n_t1,
      n_snowhk, n_sneqv, n_chk, n_cm, n_albedok, n_z0k, n_xlai, n_snotime1, n_ribb, n_eta, n_sheat, n_eta_kin,
      n_fdown, n_ec, n_edir, n_ett, n_esnow, n_drip, n_dew, n_beta, n_etp, n_ssoil, n_flx1, n_flx2, n_flx3,
      n_snomlt, n_sncovr, n_runoff1, n_runoff2, n_runoff3, n_rc, n_pc, n_rsmin, n_soilw, n_soilm, n_q1;
  Field3 n_stc, n_smc, n_swc, n_et, n_smav;
  IntField2 n_vegtyp, n_soiltyp, n_slopetyp, n_ice, n_active;
  Field3 zero3_;  // qc/qi stand-ins for YSU/MYJ (flag_qc = flag_qi = false)

  void surface_driver(const dyn::State& s, const dyn::BaseState& b);
  void pbl_driver(dyn::State& s, const dyn::BaseState& b);
  void cumulus_driver(const dyn::State& s, const dyn::BaseState& b);
  void convtrans_prep(const dyn::State& s);
  void advance_ppt(const dyn::State& s);
};

}  // namespace cppwrf::phys::driver
