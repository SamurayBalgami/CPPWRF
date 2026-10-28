#include "phys/driver/suite.hpp"

#include "core/vm_dump.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

#include "phys/bl_gwdo/gwdo.hpp"
#include "phys/bl_myj/myjpbl.hpp"
#include "phys/bl_mynn/mynnedmf.hpp"
#include "phys/bl_ysu/ysu.hpp"
#include "phys/cu_gfl/gfl_driver.hpp"
#include "phys/cu_ntiedtke/ntiedtke.hpp"
#include "phys/cu_tiedtke/tiedtke.hpp"
#include "phys/driver/noah_ice.hpp"
#include "phys/driver/suite_internal.hpp"
#include "phys/driver/wrf_constants.hpp"
#include "phys/lsm_noah/noah.hpp"

namespace cppwrf::phys::driver {

// ------------------------------------------------------------------ IntField2
IntField2::IntField2(const Grid& g, std::string name)
    : name_(std::move(name)), L_(Layout::make(g, Stagger::M, 2)),
      buf_(static_cast<int*>(dev_alloc(L_.total * sizeof(int)))) {}
IntField2::~IntField2() { dev_free(buf_); }
IntField2::IntField2(IntField2&& o) noexcept : name_(std::move(o.name_)), L_(o.L_), buf_(o.buf_) { o.buf_ = nullptr; }
IntField2& IntField2::operator=(IntField2&& o) noexcept {
  std::swap(name_, o.name_);
  std::swap(L_, o.L_);
  std::swap(buf_, o.buf_);
  return *this;
}
void IntField2::upload(const std::vector<int>& h) { copy_to_dev(buf_, h.data(), L_.total * sizeof(int)); }
void IntField2::download(std::vector<int>& h) const {
  h.resize(L_.total);
  copy_to_host(h.data(), buf_, L_.total * sizeof(int));
}

using namespace suite_detail;

namespace {

struct KAccum2 {  // ac = ac + x*dt
  V2 ac;
  C2 x;
  real dt;
  CPPWRF_HD void operator()(int i, int j) const { ac(i, j) = ac(i, j) + x(i, j) * dt; }
};

}  // namespace

// ------------------------------------------------------------------ PhysicsSuite
PhysicsSuite::PhysicsSuite(const Grid& g, const SuiteOptions& o)
    : tsk(g, Stagger::M, "tsk"), xland(g, Stagger::M, "xland"), xice(g, Stagger::M, "xice"),
      lakemask(g, Stagger::M, "lakemask"), water_depth(g, Stagger::M, "water_depth"), mavail(g, Stagger::M, "mavail"),
      tmn(g, Stagger::M, "tmn"), vegfra(g, Stagger::M, "vegfra"), shdmin(g, Stagger::M, "shdmin"),
      shdmax(g, Stagger::M, "shdmax"), snoalb(g, Stagger::M, "snoalb"), albbck(g, Stagger::M, "albbck"),
      albedo(g, Stagger::M, "albedo"), emiss(g, Stagger::M, "emiss"), embck(g, Stagger::M, "embck"),
      znt(g, Stagger::M, "znt"), z0(g, Stagger::M, "z0"), ust(g, Stagger::M, "ust"), ustm(g, Stagger::M, "ustm"),
      mol(g, Stagger::M, "mol"), rmol(g, Stagger::M, "rmol"), pblh(g, Stagger::M, "pblh"), hfx(g, Stagger::M, "hfx"),
      qfx(g, Stagger::M, "qfx"), lh(g, Stagger::M, "lh"), grdflx(g, Stagger::M, "grdflx"), qsfc(g, Stagger::M, "qsfc"),
      flhc(g, Stagger::M, "flhc"), flqc(g, Stagger::M, "flqc"), wspd(g, Stagger::M, "wspd"), br(g, Stagger::M, "br"),
      fm(g, Stagger::M, "fm"), fh(g, Stagger::M, "fh"), zol(g, Stagger::M, "zol"), regime(g, Stagger::M, "regime"),
      u10(g, Stagger::M, "u10"), v10(g, Stagger::M, "v10"), t2(g, Stagger::M, "t2"), th2(g, Stagger::M, "th2"),
      q2(g, Stagger::M, "q2"), psfc(g, Stagger::M, "psfc"), ch(g, Stagger::M, "ch"), snow(g, Stagger::M, "snow"),
      snowh(g, Stagger::M, "snowh"), snowc(g, Stagger::M, "snowc"), canwat(g, Stagger::M, "canwat"),
      lai(g, Stagger::M, "lai"), snotime(g, Stagger::M, "snotime"), qz0(g, Stagger::M, "qz0"),
      rainbl(g, Stagger::M, "rainbl"), raincv(g, Stagger::M, "raincv"), raincv_a(g, Stagger::M, "raincv_a"),
      raincv_b(g, Stagger::M, "raincv_b"), rainncv(g, Stagger::M, "rainncv"),
      rainshv(g, Stagger::M, "rainshv"), rainc(g, Stagger::M, "rainc"), rainsh(g, Stagger::M, "rainsh"),
      pratec(g, Stagger::M, "pratec"), pratesh(g, Stagger::M, "pratesh"), nca(g, Stagger::M, "nca"),
      htop(g, Stagger::M, "htop"), hbot(g, Stagger::M, "hbot"), cutop(g, Stagger::M, "cutop"),
      cubot(g, Stagger::M, "cubot"), cuppt(g, Stagger::M, "cuppt"), sr(g, Stagger::M, "sr"), glw(g, Stagger::M, "glw"),
      swdown(g, Stagger::M, "swdown"), gsw(g, Stagger::M, "gsw"), smstav(g, Stagger::M, "smstav"),
      smstot(g, Stagger::M, "smstot"), sfcrunoff(g, Stagger::M, "sfcrunoff"), udrunoff(g, Stagger::M, "udrunoff"),
      acsnom(g, Stagger::M, "acsnom"), acsnow(g, Stagger::M, "acsnow"), snopcx(g, Stagger::M, "snopcx"),
      potevp(g, Stagger::M, "potevp"), noahres(g, Stagger::M, "noahres"), sfcevp(g, Stagger::M, "sfcevp"),
      sfcexc(g, Stagger::M, "sfcexc"), achfx(g, Stagger::M, "achfx"), aclhf(g, Stagger::M, "aclhf"),
      acgrdflx(g, Stagger::M, "acgrdflx"), dx2d(g, Stagger::M, "dx2d"), var2d(g, Stagger::M, "var2d"), uoce(g, Stagger::M, "uoce"),
      voce(g, Stagger::M, "voce"), ctopo(g, Stagger::M, "ctopo"), ctopo2(g, Stagger::M, "ctopo2"),
      ck(g, Stagger::M, "ck"), cka(g, Stagger::M, "cka"), cd(g, Stagger::M, "cd"), cda(g, Stagger::M, "cda"),
      wstar(g, Stagger::M, "wstar"), delta(g, Stagger::M, "delta"), hol(g, Stagger::M, "hol"),
      thz0(g, Stagger::M, "thz0"), uz0(g, Stagger::M, "uz0"), vz0(g, Stagger::M, "vz0"),
      akhs(g, Stagger::M, "akhs"), akms(g, Stagger::M, "akms"), ct(g, Stagger::M, "ct"),
      tshltr(g, Stagger::M, "tshltr"), th10(g, Stagger::M, "th10"), qshltr(g, Stagger::M, "qshltr"),
      q10(g, Stagger::M, "q10"), pshltr(g, Stagger::M, "pshltr"), u10e(g, Stagger::M, "u10e"),
      v10e(g, Stagger::M, "v10e"), mixht(g, Stagger::M, "mixht"),
      ivgtyp(g, "ivgtyp"), isltyp(g, "isltyp"), kpbl(g, "kpbl"), lowlyr(g, "lowlyr"), td_icum(g, "td_icum"), td_row_icum(g, "td_row_icum"), oc12d(g, Stagger::M, "oc12d"), oa1(g, Stagger::M, "oa1"), oa2(g, Stagger::M, "oa2"), oa3(g, Stagger::M, "oa3"), oa4(g, Stagger::M, "oa4"), ol1(g, Stagger::M, "ol1"), ol2(g, Stagger::M, "ol2"), ol3(g, Stagger::M, "ol3"), ol4(g, Stagger::M, "ol4"), elvmax(g, Stagger::M, "elvmax"), dusfcg(g, Stagger::M, "dusfcg"), dvsfcg(g, Stagger::M, "dvsfcg"), dtaux3d(g, Stagger::M, "dtaux3d"), dtauy3d(g, Stagger::M, "dtauy3d"), pt(g), prep(g),
      psim(g, Stagger::M, "psim"), psih(g, Stagger::M, "psih"), gz1oz0(g, Stagger::M, "gz1oz0"),
      chklowq(g, Stagger::M, "chklowq"), qgh(g, Stagger::M, "qgh"), chs(g, Stagger::M, "chs"),
      chs2(g, Stagger::M, "chs2"), cqs2(g, Stagger::M, "cqs2"), cpm(g, Stagger::M, "cpm"),
      u_phytmp(g, Stagger::M, "u_phytmp"), v_phytmp(g, Stagger::M, "v_phytmp"), t1_lag_(g, Stagger::M, "t1_lag"), coszen(g, Stagger::M, "coszen"), rthratenlw(g, Stagger::M, "rthratenlw"), rthratensw(g, Stagger::M, "rthratensw"), rthratenlwc(g, Stagger::M, "rthratenlwc"), rthratenswc(g, Stagger::M, "rthratenswc"), o3rad(g, Stagger::M, "o3rad"), lwcf(g, Stagger::M, "lwcf"), swcf(g, Stagger::M, "swcf"), hrang(g, Stagger::M, "hrang"), uratx(g, Stagger::M, "uratx"), vratx(g, Stagger::M, "vratx"), tratx(g, Stagger::M, "tratx"), ra(g, Stagger::M, "ra"), rs(g, Stagger::M, "rs"), shall(g, Stagger::M, "shall"), lwupt(g, Stagger::M, "lwupt"), lwuptc(g, Stagger::M, "lwuptc"), lwdnt(g, Stagger::M, "lwdnt"), lwdntc(g, Stagger::M, "lwdntc"), lwupb(g, Stagger::M, "lwupb"), lwupbc(g, Stagger::M, "lwupbc"), lwdnb(g, Stagger::M, "lwdnb"), lwdnbc(g, Stagger::M, "lwdnbc"), olr(g, Stagger::M, "olr"), swupt(g, Stagger::M, "swupt"), swuptc(g, Stagger::M, "swuptc"), swdnt(g, Stagger::M, "swdnt"), swdntc(g, Stagger::M, "swdntc"), swupb(g, Stagger::M, "swupb"), swupbc(g, Stagger::M, "swupbc"), swdnb(g, Stagger::M, "swdnb"), swdnbc(g, Stagger::M, "swdnbc"), swvisdir(g, Stagger::M, "swvisdir"), swvisdif(g, Stagger::M, "swvisdif"), swnirdir(g, Stagger::M, "swnirdir"), swnirdif(g, Stagger::M, "swnirdif"), swddir(g, Stagger::M, "swddir"), swddni(g, Stagger::M, "swddni"), swddif(g, Stagger::M, "swddif"), swddirc(g, Stagger::M, "swddirc"), swddnic(g, Stagger::M, "swddnic"), diffuse_frac(g, Stagger::M, "diffuse_frac"), aclwupt(g, Stagger::M, "aclwupt"), aclwuptc(g, Stagger::M, "aclwuptc"), aclwdnt(g, Stagger::M, "aclwdnt"), aclwdntc(g, Stagger::M, "aclwdntc"), aclwupb(g, Stagger::M, "aclwupb"), aclwupbc(g, Stagger::M, "aclwupbc"), aclwdnb(g, Stagger::M, "aclwdnb"), aclwdnbc(g, Stagger::M, "aclwdnbc"), acswupt(g, Stagger::M, "acswupt"), acswuptc(g, Stagger::M, "acswuptc"), acswdnt(g, Stagger::M, "acswdnt"), acswdntc(g, Stagger::M, "acswdntc"), acswupb(g, Stagger::M, "acswupb"), acswupbc(g, Stagger::M, "acswupbc"), acswdnb(g, Stagger::M, "acswdnb"), acswdnbc(g, Stagger::M, "acswdnbc"), g_(g), gs_(g), o_(o) {
  if ((o.sf_sfclay_physics != 1 && o.sf_sfclay_physics != 2 && o.sf_sfclay_physics != 5) ||
      (o.sf_surface_physics != 2 && o.sf_surface_physics != 4) ||
      (o.bl_pbl_physics != 0 && o.bl_pbl_physics != 1 && o.bl_pbl_physics != 2 && o.bl_pbl_physics != 5))
    throw std::runtime_error("PhysicsSuite: only sf_sfclay_physics=1|2|5, sf_surface_physics=2|4, bl_pbl_physics=0|1|2|5");
  // module_physics_init: YSUSCHEME is fatal without sfclay (isfc = 1)
  if (o.bl_pbl_physics == 1 && o.sf_sfclay_physics != 1)
    throw std::runtime_error("PhysicsSuite: bl_pbl_physics=1 (YSU) needs sf_sfclay_physics=1, as WRF's module_physics_init");
  // module_physics_init: MYNNPBLSCHEME needs mynnsfc, sfclay or myjsfc; only mynnsfc is ported with it
  if (o.bl_pbl_physics == 5 && o.sf_sfclay_physics != 5)
    throw std::runtime_error("PhysicsSuite: bl_pbl_physics=5 is ported with sf_sfclay_physics=5 only");
  // module_physics_init.F pbl_select: YSU stops unless isfc = 1 (sfclay), MYJ unless isfc = 2 (myjsfc)
  if (o.bl_pbl_physics == 1 && o.sf_sfclay_physics != 1)
    throw std::runtime_error("PhysicsSuite: bl_pbl_physics=1 (YSU) needs sf_sfclay_physics=1, as WRF");
  if (o.bl_pbl_physics == 2 && o.sf_sfclay_physics != 2)
    throw std::runtime_error("PhysicsSuite: bl_pbl_physics=2 (MYJ) needs sf_sfclay_physics=2, as WRF");
  if (o.sf_sfclay_physics == 5 && o.fractional_seaice == 1)
    throw std::runtime_error("PhysicsSuite: sf_sfclay_physics=5 with fractional_seaice=1 (MYNN_SEAICE_WRAPPER) is not ported");
  if (o.cu_physics != 0 && o.cu_physics != 1 && o.cu_physics != 3 && o.cu_physics != 6 && o.cu_physics != 16)
    throw std::runtime_error("PhysicsSuite: only cu_physics=0, 1 (Kain-Fritsch), 3 (Grell-Freitas-Li), 6 (Tiedtke) or 16 (New Tiedtke)");
  if (o.cu_physics == 3 && g.nz > cu_gfl::kMaxLevels)
    throw std::runtime_error("PhysicsSuite: nz exceeds the Grell-Freitas-Li column limit");
  // ishallow: WRF's GF runs the shallow plume for ishallow == 1 only, any other value is off.
  // clos_choice above 16 reads past xf_ens (undefined in WRF too); a negative one is the full ensemble.
  if (o.cu_physics == 3 && o.gf_ichoice > 16)
    throw std::runtime_error("PhysicsSuite: clos_choice must be at most 16");
  if (o.cu_physics == 6 && g.nz > cu_tiedtke::kMaxLevels)
    throw std::runtime_error("PhysicsSuite: nz exceeds the Tiedtke column limit");
  if (o.cu_physics == 16 && g.nz > cu_ntiedtke::kMaxLevels)
    throw std::runtime_error("PhysicsSuite: nz exceeds the New Tiedtke column limit");
  if (o.opt_thcnd != 1 && o.opt_thcnd != 2)
    throw std::runtime_error("PhysicsSuite: opt_thcnd must be 1 or 2");
  if (o.cu_physics == 1 && o.kfeta_trigger != 1 && o.kfeta_trigger != 3)
    throw std::runtime_error("PhysicsSuite: kfeta_trigger=2 is not supported");
  if (o.fractional_seaice != 0 && o.fractional_seaice != 1)
    throw std::runtime_error("PhysicsSuite: fractional_seaice must be 0 or 1");
  if (o.sst_update != 0 && o.sst_update != 1) throw std::runtime_error("PhysicsSuite: sst_update must be 0 or 1");
  if (o.fractional_seaice == 1 || o.sst_update == 1) {
    sst = Field2(g, Stagger::M, "sst");
    tsk_save = Field2(g, Stagger::M, "tsk_save");
  }
  if (o.fractional_seaice == 1) sea_ = std::make_unique<SeaIceFields>(g);
  if (o.sst_update == 1) xicem = Field2(g, Stagger::M, "xicem");
  if (o.seaice_thickness_opt != 0 || o.seaice_snowdepth_opt != 0)
    throw std::runtime_error("PhysicsSuite: seaice_thickness_opt/seaice_snowdepth_opt 1 need ICEDEPTH/SNOWSI input");
  if (o.seaice_albedo_opt != 0)
    throw std::runtime_error("PhysicsSuite: seaice_albedo_opt 1 and 2 are not supported (2 needs ALBSI input)");
  if (o.num_soil_layers > noah::kMaxSoilLevels) throw std::runtime_error("PhysicsSuite: too many soil layers");
  if (o.sf_surface_physics == 4) {
    if (o.num_soil_layers != 4) throw std::runtime_error("PhysicsSuite: Noah-MP requires num_soil_layers=4");
    if (o.mminlu == "NLCD40") throw std::runtime_error("PhysicsSuite: Noah-MP does not support NLCD40 landuse");
    noahmp::validate_noahmp_options(o.noahmp);
  }
  if (o.bl_pbl_physics == 1 && g.nz > ysu::kMaxLevels)
    throw std::runtime_error("PhysicsSuite: nz exceeds the YSU column limit");
  if (o.bl_pbl_physics == 5 && g.nz > mynnedmf::kMaxLevels)
    throw std::runtime_error("PhysicsSuite: nz exceeds the MYNN-EDMF column limit");
  if (o.bl_pbl_physics == 2 && g.nz > myjpbl::kMaxLevels)
    throw std::runtime_error("PhysicsSuite: nz exceeds the MYJ column limit");
  if (o.sf_sfclay_physics == 2 && g.nz > sf_myjsfc::kMaxLevels)
    throw std::runtime_error("PhysicsSuite: nz exceeds the MYJSFC column limit");
  if (o.gwd_opt != 0 && o.gwd_opt != 1)
    throw std::runtime_error("PhysicsSuite: only gwd_opt=0 or 1 (GWDO)");
  if (o.gwd_opt == 1 && g.nz > bl_gwdo::kMaxLevels)
    throw std::runtime_error("PhysicsSuite: nz exceeds the GWDO column limit");
  if (o.cu_physics == 1 && g.nz > cu_kf::kMaxLevels)
    throw std::runtime_error("PhysicsSuite: nz exceeds the KF column limit");
  gs_.nz = o.num_soil_layers;
  auto mk3 = [&](Field3& f, const Grid& gg, const char* n) { f = Field3(gg, Stagger::M, n); };
  mk3(tslb, gs_, "tslb"); mk3(smois, gs_, "smois"); mk3(sh2o, gs_, "sh2o"); mk3(smcrel, gs_, "smcrel");
  mk3(exch_h, g, "exch_h"); mk3(exch_m, g, "exch_m");
  mk3(tke_pbl, g, "tke_pbl"); mk3(el_pbl, g, "el_pbl");
  mk3(w0avg, g, "w0avg"); mk3(cldfra_dp, g, "cldfra_dp"); mk3(cldfra_sh, g, "cldfra_sh");
  mk3(qc_cu, g, "qc_cu"); mk3(qi_cu, g, "qi_cu"); mk3(cldfra, g, "cldfra");
  mk3(gd_cloud, g, "gd_cloud"); mk3(gd_cloud2, g, "gd_cloud2");
  mk3(gd_cloud_a, g, "gd_cloud_a"); mk3(gd_cloud2_a, g, "gd_cloud2_a");
  mk3(gd_cldfr, g, "gd_cldfr");
  mk3(zero3_, g, "zero3");
  dzs.assign(o.num_soil_layers, 0.0);
  if (o.sf_surface_physics == 2) {
    mk3(n_stc, gs_, "n_stc"); mk3(n_smc, gs_, "n_smc"); mk3(n_swc, gs_, "n_swc"); mk3(n_et, gs_, "n_et");
    mk3(n_smav, gs_, "n_smav");
    for (Field2* f : {&n_lwdn, &n_soldn, &n_solnet, &n_sfcprs, &n_prcp, &n_sfctmp, &n_th2, &n_q2, &n_q2sat, &n_dqsdt2,
                      &n_ffrozp, &n_shdfac, &n_shmin, &n_shmax, &n_albbrd, &n_snoalb1, &n_tbot, &n_z0brd, &n_emissi,
                      &n_embrd, &n_cmc, &n_t1, &n_snowhk, &n_sneqv, &n_chk, &n_cm, &n_albedok, &n_z0k, &n_xlai,
                      &n_snotime1, &n_ribb, &n_eta, &n_sheat, &n_eta_kin, &n_fdown, &n_ec, &n_edir, &n_ett, &n_esnow,
                      &n_drip, &n_dew, &n_beta, &n_etp, &n_ssoil, &n_flx1, &n_flx2, &n_flx3, &n_snomlt, &n_sncovr,
                      &n_runoff1, &n_runoff2, &n_runoff3, &n_rc, &n_pc, &n_rsmin, &n_soilw, &n_soilm, &n_q1})
      *f = Field2(g, Stagger::M, "noah_scratch");
    n_vegtyp = IntField2(g, "n_vegtyp"); n_soiltyp = IntField2(g, "n_soiltyp"); n_slopetyp = IntField2(g, "n_slopetyp");
    n_ice = IntField2(g, "n_ice"); n_active = IntField2(g, "n_active");
    if (!noah::read_noah_tables(o.table_dir, o.mminlu, o.mminsl, noah_tables_))
      throw std::runtime_error("PhysicsSuite: cannot read the Noah tables from " + o.table_dir);
  } else if (o.sf_surface_physics == 4) {
    auto mk2 = [&](Field2& f, const char* n) { f = Field2(g, Stagger::M, n); };
    mk2(xlat_sfc, "xlat_sfc"); mk2(xlong_sfc, "xlong_sfc");
    mk2(tvxy, "tvxy"); mk2(tgxy, "tgxy"); mk2(canicexy, "canicexy"); mk2(canliqxy, "canliqxy");
    mk2(eahxy, "eahxy"); mk2(tahxy, "tahxy"); mk2(cmxy, "cmxy"); mk2(chxy, "chxy");
    mk2(fwetxy, "fwetxy"); mk2(sneqvoxy, "sneqvoxy"); mk2(alboldxy, "alboldxy");
    mk2(qsnowxy, "qsnowxy"); mk2(qrainxy, "qrainxy"); mk2(wslakexy, "wslakexy");
    mk2(zwtxy, "zwtxy"); mk2(waxy, "waxy"); mk2(wtxy, "wtxy");
    mk2(lfmassxy, "lfmassxy"); mk2(rtmassxy, "rtmassxy"); mk2(stmassxy, "stmassxy");
    mk2(woodxy, "woodxy"); mk2(stblcpxy, "stblcpxy"); mk2(fastcpxy, "fastcpxy");
    mk2(xsaixy, "xsaixy"); mk2(taussxy, "taussxy"); mk2(smcwtdxy, "smcwtdxy");
    mk2(deeprechxy, "deeprechxy"); mk2(rechxy, "rechxy"); mk2(grainxy, "grainxy");
    mk2(gddxy, "gddxy"); mk2(qtdrain, "qtdrain");
    mk2(t2mvxy, "t2mvxy"); mk2(t2mbxy, "t2mbxy"); mk2(q2mvxy, "q2mvxy"); mk2(q2mbxy, "q2mbxy");
    mk2(tradxy, "tradxy"); mk2(neexy, "neexy"); mk2(gppxy, "gppxy"); mk2(nppxy, "nppxy");
    mk2(fvegxy, "fvegxy"); mk2(runsfxy, "runsfxy"); mk2(runsbxy, "runsbxy");
    mk2(ecanxy, "ecanxy"); mk2(edirxy, "edirxy"); mk2(etranxy, "etranxy");
    mk2(fsaxy, "fsaxy"); mk2(firaxy, "firaxy"); mk2(aparxy, "aparxy"); mk2(psnxy, "psnxy");
    mk2(savxy, "savxy"); mk2(sagxy, "sagxy"); mk2(rssunxy, "rssunxy"); mk2(rsshaxy, "rsshaxy");
    mk2(bgapxy, "bgapxy"); mk2(wgapxy, "wgapxy"); mk2(tgvxy, "tgvxy"); mk2(tgbxy, "tgbxy");
    mk2(chvxy, "chvxy"); mk2(chbxy, "chbxy"); mk2(shgxy, "shgxy"); mk2(shcxy, "shcxy");
    mk2(shbxy, "shbxy"); mk2(evgxy, "evgxy"); mk2(evbxy, "evbxy"); mk2(ghvxy, "ghvxy");
    mk2(ghbxy, "ghbxy"); mk2(irgxy, "irgxy"); mk2(ircxy, "ircxy"); mk2(irbxy, "irbxy");
    mk2(trxy, "trxy"); mk2(evcxy, "evcxy"); mk2(chleafxy, "chleafxy"); mk2(chucxy, "chucxy");
    mk2(chv2xy, "chv2xy"); mk2(chb2xy, "chb2xy"); mk2(chstarxy, "chstarxy");
    mk2(acc_ssoilxy, "acc_ssoil"); mk2(acc_qinsurxy, "acc_qinsur"); mk2(acc_qsevaxy, "acc_qseva");
    mk2(eflxbxy, "eflxbxy");
    mk2(acc_dwaterxy, "acc_dwaterxy"); mk2(acc_prcpxy, "acc_prcpxy"); mk2(acc_ecanxy, "acc_ecanxy");
    mk2(acc_etranxy, "acc_etranxy"); mk2(acc_edirxy, "acc_edirxy"); mk2(acc_glaflwxy, "acc_glaflwxy");
    isnowxy = IntField2(g, "isnowxy"); cropcat = IntField2(g, "cropcat"); pgsxy = IntField2(g, "pgsxy");
    isnowxy.upload(std::vector<int>(isnowxy.layout().total, 0));
    cropcat.upload(std::vector<int>(cropcat.layout().total, 0));
    pgsxy.upload(std::vector<int>(pgsxy.layout().total, 0));
    noahmp_ranland = IntField2(g, "noahmp_ranland");
    noahmp_ranland.upload(std::vector<int>(noahmp_ranland.layout().total, 0));
    Grid gsno = g; gsno.nz = noahmp::kMaxSnowLayers;
    Grid gsnso = g; gsnso.nz = noahmp::kMaxSnowLayers + o.num_soil_layers;
    Grid gband = g; gband.nz = noahmp::kNumRadBands;
    mk3(tsnoxy, gsno, "tsnoxy"); mk3(snicexy, gsno, "snicexy"); mk3(snliqxy, gsno, "snliqxy");
    mk3(zsnsoxy, gsnso, "zsnsoxy");
    mk3(albsoildirxy, gband, "albsoildirxy"); mk3(albsoildifxy, gband, "albsoildifxy");
    mk3(smoiseq, gs_, "smoiseq"); mk3(acc_etranixy, gs_, "acc_etrani");
    if (!noahmp::read_noahmp_tables(o.table_dir, o.mminlu, noahmp_tables_))
      throw std::runtime_error("PhysicsSuite: cannot read NoahmpTable.TBL from " + o.table_dir);
    d_noahmp_tables_ = static_cast<noahmp::NoahmpTables*>(dev_alloc(sizeof(noahmp_tables_)));
    copy_to_dev(d_noahmp_tables_, &noahmp_tables_, sizeof(noahmp_tables_));
  }
  if (o.sf_sfclay_physics == 1) {
    sfclayrev::SFClayRevTables<real> h{};
    sfclayrev::init_sfclayrev_tables(h);
    d_sfclay_tables_ = static_cast<sfclayrev::SFClayRevTables<real>*>(dev_alloc(sizeof(h)));
    copy_to_dev(d_sfclay_tables_, &h, sizeof(h));
  }
  if (o.sf_sfclay_physics == 2) {
    lowlyr.upload(std::vector<int>(lowlyr.layout().total, 1));  // myjsfcinit: LOWLYR = 1
    ust.fill_host([](int, int, int) { return 0.1; });           // myjsfcinit: USTAR = 0.1
  }
  if (o.sf_sfclay_physics == 5) {  // mynnsfc_init -> psi_init
    mynnsfc::PsiTables* h = new mynnsfc::PsiTables;
    mynnsfc::psi_init(*h);
    d_mynn_psi_ = static_cast<mynnsfc::PsiTables*>(dev_alloc(sizeof(*h)));
    copy_to_dev(d_mynn_psi_, h, sizeof(*h));
    delete h;
    cqs = Field2(g, Stagger::M, "cqs");
  }
  if (o.bl_pbl_physics == 5) {
    mk3(qke, g, "qke"); mk3(tsq, g, "tsq"); mk3(qsq, g, "qsq"); mk3(cov, g, "cov"); mk3(sh3d, g, "sh3d");
    mk3(sm3d, g, "sm3d"); mk3(qc_bl, g, "qc_bl"); mk3(qi_bl, g, "qi_bl"); mk3(cldfra_bl, g, "cldfra_bl");
    maxwidth = Field2(g, Stagger::M, "maxwidth"); maxmf = Field2(g, Stagger::M, "maxmf");
    ztop_plume = Field2(g, Stagger::M, "ztop_plume"); excess_h = Field2(g, Stagger::M, "excess_h");
    excess_q = Field2(g, Stagger::M, "excess_q");
  }
  if (o.bl_pbl_physics == 2) {
    lowlyr.upload(std::vector<int>(lowlyr.layout().total, 1));
    tke_pbl.fill_host([](int, int, int) { return 0.2; });       // myjpblinit: TKE_MYJ = EPSQ2 = 0.2
  }
  if (o.cu_physics == 1) {
    using namespace wrfc;
    cu_kf::KFTables<real> h_kf{};
    cu_kf::init_kf_tables(h_kf, svp1, svp2, svp3, svpt0);
    d_kf_tables_ = static_cast<cu_kf::KFTables<real>*>(dev_alloc(sizeof(h_kf)));
    copy_to_dev(d_kf_tables_, &h_kf, sizeof(h_kf));
    nca.fill_host([](int, int, int) { return -100.0; });  // kf_eta_init default for non-restart
  }
  // bldt = 0: every step; otherwise WRF stepbl = nint(bldt*60/dt), at least 1
  stepbl = o.bldt == 0 ? 1 : int(std::lround(o.bldt * 60. / g.dt));
  if (stepbl < 1) stepbl = 1;
  // cu_init (module_physics_init.F:4025-4026): STEPCU = max(nint(CUDT*60./DT), 1)
  stepcu = int(std::lround(o.cudt * 60. / g.dt));
  if (stepcu < 1) stepcu = 1;
}

void PhysicsSuite::init_radiation(const std::vector<double>& xlat, const std::vector<double>& xlong) {
  if (!radiation_on()) return;
  rad_ = std::make_unique<RadiationDriver>(g_, o_.rad, xlat, xlong);
}

void PhysicsSuite::init_noahmp(bool fndsnowh) {
  if (o_.sf_surface_physics != 4) return;
  noahmp::NoahmpInitKernel k{};
  k.ivgtyp = ivgtyp.cview();
  k.isltyp = isltyp.cview();
  k.tsk = tsk.cview();
  k.xice = xice.cview();
  k.tmn = tmn.view();
  k.snow = snow.view();
  k.snowh = snowh.view();
  k.canwat = canwat.view();
  k.tslb = tslb.view();
  k.smois = smois.view();
  k.sh2o = sh2o.view();
  k.isnowxy = isnowxy.view();
  k.tvxy = tvxy.view();
  k.tgxy = tgxy.view();
  k.canicexy = canicexy.view();
  k.canliqxy = canliqxy.view();
  k.eahxy = eahxy.view();
  k.tahxy = tahxy.view();
  k.cmxy = cmxy.view();
  k.chxy = chxy.view();
  k.fwetxy = fwetxy.view();
  k.sneqvoxy = sneqvoxy.view();
  k.alboldxy = alboldxy.view();
  k.qsnowxy = qsnowxy.view();
  k.qrainxy = qrainxy.view();
  k.wslakexy = wslakexy.view();
  k.zwtxy = zwtxy.view();
  k.waxy = waxy.view();
  k.wtxy = wtxy.view();
  k.lfmassxy = lfmassxy.view();
  k.rtmassxy = rtmassxy.view();
  k.stmassxy = stmassxy.view();
  k.woodxy = woodxy.view();
  k.stblcpxy = stblcpxy.view();
  k.fastcpxy = fastcpxy.view();
  k.xsaixy = xsaixy.view();
  k.lai = lai.view();
  k.grainxy = grainxy.view();
  k.gddxy = gddxy.view();
  k.cropcat = cropcat.view();
  k.t2mvxy = t2mvxy.view();
  k.t2mbxy = t2mbxy.view();
  k.chstarxy = chstarxy.view();
  k.qtdrain = qtdrain.view();
  k.tsnoxy = tsnoxy.view();
  k.snicexy = snicexy.view();
  k.snliqxy = snliqxy.view();
  k.zsnsoxy = zsnsoxy.view();
  k.albsoildirxy = albsoildirxy.view();
  k.albsoildifxy = albsoildifxy.view();
  k.tables = d_noahmp_tables_;
  k.nsoil = o_.num_soil_layers;
  k.nsnow = noahmp::kMaxSnowLayers;
  for (int ns = 0; ns < o_.num_soil_layers; ++ns) k.dzs[ns] = real(dzs[ns]);
  k.fndsnowh = fndsnowh;
  k.restart = false;
  k.iopt_runsub = o_.noahmp.opt_runsub;
  k.iopt_crop = o_.noahmp.opt_crop;
  k.sf_urban_physics = o_.noahmp.sf_urban_physics;
  noahmp::run_noahmp_init(k, g_.columns());
}

PhysicsSuite::~PhysicsSuite() {
  if (d_noahmp_tables_) dev_free(d_noahmp_tables_);
  if (d_sfclay_tables_) dev_free(d_sfclay_tables_);
  if (d_kf_tables_) dev_free(d_kf_tables_);
  if (d_mynn_psi_) dev_free(d_mynn_psi_);
}

PhysicsSuite::SeaIceFields::SeaIceFields(const Grid& g) {
  for (Field2* f : {&tsk_local, &tsk_sea, &xland_sea, &mavail_sea, &znt_sea, &z0base_sea, &br_sea, &chs_sea, &chs2_sea,
                    &cqs2_sea, &cpm_sea, &flhc_sea, &flqc_sea, &gz1oz0_sea, &hfx_sea, &qfx_sea, &lh_sea, &mol_sea,
                    &psim_sea, &psih_sea, &fm_sea, &fh_sea, &qgh_sea, &qsfc_sea, &regime_sea, &rmol_sea, &ust_sea,
                    &ustm_sea, &wspd_sea, &zol_sea, &u10_sea, &v10_sea, &th2_sea, &t2_sea, &q2_sea, &ck_sea, &cka_sea,
                    &cd_sea, &cda_sea, &qz0_sea, &thz0_sea, &uz0_sea, &vz0_sea, &pblh_sea, &akhs_sea, &akms_sea,
                    &ct_sea, &tshltr_sea, &th10_sea, &qshltr_sea, &q10_sea, &pshltr_sea, &u10e_sea, &v10e_sea,
                    &regime_hold})
    *f = Field2(g, Stagger::M, "seaice_scratch");
}

std::vector<std::pair<std::string, Field2*>> PhysicsSuite::carried2() {
  std::vector<std::pair<std::string, Field2*>> v;
  for (Field2* f : {&tsk, &xland, &xice, &lakemask, &water_depth, &mavail, &tmn, &vegfra, &shdmin, &shdmax, &snoalb,
                    &albbck, &albedo, &emiss, &embck, &znt, &z0, &ust, &ustm, &mol, &rmol, &pblh, &hfx, &qfx, &lh,
                    &grdflx, &qsfc, &flhc, &flqc, &wspd, &br, &fm, &fh, &zol, &regime, &u10, &v10, &t2, &th2, &q2,
                    &psfc, &ch, &snow, &snowh, &snowc, &canwat, &lai, &snotime, &qz0, &rainbl, &raincv, &raincv_a,
                    &raincv_b, &rainncv,
                    &rainshv, &sr, &glw, &swdown, &gsw, &smstav, &smstot, &sfcrunoff, &udrunoff, &acsnom, &acsnow,
                    &snopcx, &potevp, &noahres, &sfcevp, &sfcexc, &achfx, &aclhf, &acgrdflx, &dx2d, &uoce, &voce,
                    &ctopo, &ctopo2, &ck, &cka, &cd, &cda, &wstar, &delta, &hol, &rainc, &rainsh, &pratec, &pratesh, &nca, &htop, &hbot,
                    &cutop, &cubot, &cuppt, &coszen})
    v.emplace_back("sfc_" + f->name(), f);
  if (o_.sf_sfclay_physics == 2 || o_.bl_pbl_physics == 2) {
    for (Field2* f : {&thz0, &uz0, &vz0, &akhs, &akms, &ct, &tshltr, &th10, &qshltr, &q10, &pshltr, &u10e, &v10e, &mixht})
      v.emplace_back("sfc_" + f->name(), f);
  }
  if (o_.sf_sfclay_physics == 5) v.emplace_back("sfc_cqs", &cqs);
  if (o_.bl_pbl_physics == 5)
    for (Field2* f : {&maxwidth, &maxmf, &ztop_plume, &excess_h, &excess_q}) v.emplace_back("sfc_" + f->name(), f);
  if (o_.fractional_seaice == 1 || o_.sst_update == 1) {
    v.emplace_back("sfc_sst", &sst);
    v.emplace_back("sfc_tsk_save", &tsk_save);
  }
  if (o_.sst_update == 1) v.emplace_back("sfc_xicem", &xicem);
  if (rad_)  // radiation diagnostics and their accumulations
    for (Field2* f : {&lwcf, &swcf, &lwupt, &lwuptc, &lwdnt, &lwdntc, &lwupb, &lwupbc, &lwdnb, &lwdnbc, &olr, &swupt, &swuptc, &swdnt, &swdntc, &swupb, &swupbc, &swdnb, &swdnbc, &swvisdir, &swvisdif, &swnirdir, &swnirdif, &swddir, &swddni, &swddif, &swddirc, &swddnic, &diffuse_frac, &aclwupt, &aclwuptc, &aclwdnt, &aclwdntc, &aclwupb, &aclwupbc, &aclwdnb, &aclwdnbc, &acswupt, &acswuptc, &acswdnt, &acswdntc, &acswupb, &acswupbc, &acswdnb, &acswdnbc})
      v.emplace_back("sfc_" + f->name(), f);
  if (o_.sf_surface_physics == 4) {
    for (Field2* f : {&xlat_sfc, &xlong_sfc, &tvxy, &tgxy, &canicexy, &canliqxy, &eahxy, &tahxy, &cmxy, &chxy,
                      &fwetxy, &sneqvoxy, &alboldxy, &qsnowxy, &qrainxy, &wslakexy, &zwtxy, &waxy, &wtxy, &lfmassxy,
                      &rtmassxy, &stmassxy, &woodxy, &stblcpxy, &fastcpxy, &xsaixy, &taussxy, &smcwtdxy, &deeprechxy,
                      &rechxy, &grainxy, &gddxy, &qtdrain, &t2mvxy, &t2mbxy, &q2mvxy, &q2mbxy, &tradxy, &neexy,
                      &gppxy, &nppxy, &fvegxy, &runsfxy, &runsbxy, &ecanxy, &edirxy, &etranxy, &fsaxy, &firaxy,
                      &aparxy, &psnxy, &savxy, &sagxy, &rssunxy, &rsshaxy, &bgapxy, &wgapxy, &tgvxy, &tgbxy, &chvxy,
                      &chbxy, &shgxy, &shcxy, &shbxy, &evgxy, &evbxy, &ghvxy, &ghbxy, &irgxy, &ircxy, &irbxy, &trxy,
                      &evcxy, &chleafxy, &chucxy, &chv2xy, &chb2xy, &chstarxy, &rs, &acc_ssoilxy, &acc_qinsurxy,
                      &acc_qsevaxy, &eflxbxy, &acc_dwaterxy, &acc_prcpxy, &acc_ecanxy, &acc_etranxy, &acc_edirxy, &acc_glaflwxy})
      v.emplace_back("sfc_" + f->name(), f);
  }
  return v;
}

std::vector<std::pair<std::string, Field3*>> PhysicsSuite::carried3() {
  std::vector<std::pair<std::string, Field3*>> v;
  for (Field3* f : {&tslb, &smois, &sh2o, &smcrel, &exch_h, &exch_m, &pt.rthraten, &pt.rublten, &pt.rvblten,
                    &pt.rthblten, &pt.rqvblten, &pt.rqcblten, &pt.rqiblten, &t1_lag_, &w0avg, &cldfra_dp, &cldfra_sh,
                    &qc_cu, &qi_cu, &gd_cloud, &gd_cloud2, &gd_cloud_a, &gd_cloud2_a, &gd_cldfr,
                    &pt.rucuten, &pt.rvcuten, &pt.rthcuten, &pt.rqvcuten, &pt.rqccuten, &pt.rqrcuten,
                    &pt.rqicuten, &pt.rqscuten})
    v.emplace_back("sfc_" + f->name(), f);
  if (o_.sf_sfclay_physics == 2 || o_.bl_pbl_physics == 2) {
    v.emplace_back("sfc_tke_pbl", &tke_pbl);
    v.emplace_back("sfc_el_pbl", &el_pbl);
  }
  if (o_.bl_pbl_physics == 5) {
    for (Field3* f : {&qke, &tsq, &qsq, &cov, &sh3d, &sm3d, &qc_bl, &qi_bl, &cldfra_bl})
      v.emplace_back("sfc_" + f->name(), f);
    v.emplace_back("sfc_el_pbl", &el_pbl);
  }
  if (o_.sf_surface_physics == 4) {
    for (Field3* f : {&tsnoxy, &snicexy, &snliqxy, &zsnsoxy, &albsoildirxy, &albsoildifxy, &smoiseq, &acc_etranixy})
      v.emplace_back("sfc_" + f->name(), f);
  }
  if (rad_) {
    v.emplace_back("sfc_cldfra", &cldfra);
    v.emplace_back("sfc_rthratenlw", &rthratenlw);
    v.emplace_back("sfc_rthratensw", &rthratensw);
    v.emplace_back("sfc_rthratenlwc", &rthratenlwc);
    v.emplace_back("sfc_rthratenswc", &rthratenswc);
  }
  if (o_.cu_physics == 16 || o_.cu_physics == 3) {
    v.emplace_back("sfc_rthften", &pt.rthften);
    v.emplace_back("sfc_rqvften", &pt.rqvften);
  }
  if (o_.cu_physics == 6) v.emplace_back("sfc_rqvften", &pt.rqvften);  // read by the next step's cu_tiedtke
  for (Field3* f : {re_cloud_, re_ice_, re_snow_})
    if (f) v.emplace_back("sfc_" + f->name(), f);
  return v;
}

std::vector<std::pair<std::string, IntField2*>> PhysicsSuite::carried_int() {
  std::vector<std::pair<std::string, IntField2*>> v = {{"sfc_ivgtyp", &ivgtyp}, {"sfc_isltyp", &isltyp}, {"sfc_kpbl", &kpbl}};
  if (o_.sf_sfclay_physics == 2 || o_.bl_pbl_physics == 2) v.emplace_back("sfc_lowlyr", &lowlyr);
  if (o_.sf_surface_physics == 4) {
    v.emplace_back("sfc_isnowxy", &isnowxy);
    v.emplace_back("sfc_cropcat", &cropcat);
    v.emplace_back("sfc_pgsxy", &pgsxy);
  }
  return v;
}

void PhysicsSuite::part1(dyn::State& s, const dyn::BaseState& b, dyn::Tendencies&) {
  ++itimestep;  // solve_em: grid%itimestep = grid%itimestep + 1
  if (rainncv_mp_) rainncv.copy_from(*rainncv_mp_);  // the microphysics' grid%rainncv and grid%sr
  if (sr_mp_) sr.copy_from(*sr_mp_);
  run_phy_prep(s, b, prep, PhyPrepParams{o_.use_theta_m});
  if (rad_) {  // radiation_driver: xtime (minutes) and julian of the step start
    const double xtime = double(itimestep - 1) * g_.dt / 60.;
    const double julian = wrf_julian(o_.rad.julyr, o_.rad.julday, o_.rad.gmt, double(itimestep - 1) * g_.dt);
    auto q = [&](int n) -> const Field3* { return n >= 0 ? &s.q_2[n] : nullptr; };
    RadiationInputs ri{&prep.t_phy, &prep.p_hyd, &prep.p_hyd_w, &prep.pi_phy, &prep.rho, &prep.dz8w,
                       &prep.t8w, &s.q_2[s.P_QV], q(s.P_QC), q(s.P_QI), q(s.P_QS), re_cloud_, re_ice_,
                       re_snow_, &tsk, &albedo, &emiss, &xland, &xice, &snow};
    // check_a_mundo: icloud_bl = 1 needs the MYNN PBL (reset to 0 otherwise)
    if (o_.icloud_bl == 1 && o_.bl_pbl_physics == 5) {
      ri.cldfra_bl = &cldfra_bl;
      ri.qc_bl = &qc_bl;
      ri.qi_bl = &qi_bl;
    }
    RadiationOutputs ro{&pt.rthraten, &glw, &gsw, &swdown, &coszen, &cldfra, &rthratenlw, &rthratensw,
                        &rthratenlwc, &rthratenswc, &o3rad, &lwcf, &swcf, &hrang,
                        {&lwupt, &lwuptc, &lwdnt, &lwdntc, &lwupb, &lwupbc, &lwdnb, &lwdnbc, &olr}, {&swupt, &swuptc, &swdnt, &swdntc, &swupb, &swupbc, &swdnb, &swdnbc},
                        {&swvisdir, &swvisdif, &swnirdir, &swnirdif, &swddir, &swddni, &swddif, &swddirc, &swddnic, &diffuse_frac}};
    rad_->run(itimestep, xtime, julian, ri, ro);
  }
  if (rad_) {  // radiation_driver: accumulated fluxes, every step (DTaccum = DT)
    const real dt = real(g_.dt);
    const Range2D cols{0, g_.nx, 0, g_.ny};
    for (auto [ac, x] : std::initializer_list<std::pair<Field2*, const Field2*>>{{&aclwupt, &lwupt}, {&aclwuptc, &lwuptc}, {&aclwdnt, &lwdnt}, {&aclwdntc, &lwdntc}, {&aclwupb, &lwupb}, {&aclwupbc, &lwupbc}, {&aclwdnb, &lwdnb}, {&aclwdnbc, &lwdnbc}, {&acswupt, &swupt}, {&acswuptc, &swuptc}, {&acswdnt, &swdnt}, {&acswdntc, &swdntc}, {&acswupb, &swupb}, {&acswupbc, &swupbc}, {&acswdnb, &swdnb}, {&acswdnbc, &swdnbc}})
      column_for(cols, KAccum2{ac->view(), x->cview(), dt});
  }
  vmdump::dump("vm_phy_prep", 1,
               {&prep.th_phy, &prep.pi_phy, &prep.t_phy, &prep.p_hyd, &prep.p_hyd_w, &prep.p8w, &prep.t8w, &prep.rho,
                &prep.dz8w, &prep.z, &prep.u_phy, &prep.v_phy, &prep.p_phy});
  surface_driver(s, b);
  vmdump::dump("vm_surface", 1, {},
               {&hfx, &qfx, &ust, &znt, &mol, &qsfc, &tsk, &flhc, &flqc, &br, &psim, &psih, &wspd, &u10, &v10, &t2,
                &th2, &q2, &lh, &grdflx, &pblh, &regime, &zol, &fm, &fh, &rmol, &canwat});
  pbl_driver(s, b);
  vmdump::dump("vm_pbl", 1, {&pt.rublten, &pt.rvblten, &pt.rthblten, &pt.rqvblten, &exch_h}, {&pblh});
  cumulus_driver(s, b);
  if (o_.cu_physics == 3) convtrans_prep(s);
}

void PhysicsSuite::part2(const dyn::State& s, const dyn::BaseState& b, dyn::Tendencies& t) {
  const CouplingFlags f{radiation_on(), o_.bl_pbl_physics > 0, o_.cu_physics, o_.use_theta_m == 1};
  calculate_phy_tend(pt, s, b, g_, f);
  update_phy_ten(t, pt, s, g_, f);
  if (o_.use_theta_m == 1) {
    vmdump::dump("vm_conv_in", 0, {&t.t_tendf, &t.q_tendf[s.P_QV], wrf_t1_lag && t1_lag_set_ ? &t1_lag_ : &s.t_1,
                                   &s.q_1[s.P_QV]});
    if (!wrf_t1_lag) {
      conv_t_tendf_to_moist(t, s.t_1, s, g_);
    } else {  // grid%t_1 lags one step behind the dycore's t_1 here (see conv_t_tendf_to_moist)
      if (!t1_lag_set_) t1_lag_.copy_from(s.t_1);
      conv_t_tendf_to_moist(t, t1_lag_, s, g_);
      t1_lag_.copy_from(s.t_1);
      t1_lag_set_ = true;
    }
    vmdump::dump("vm_conv_out", 0, {&t.t_tendf});
  }
}

void PhysicsSuite::after_rk(const dyn::State& s, const dyn::BaseState& b) {
  advance_ppt(s);
  phy_prep_part2(pt, s, b, g_, CouplingFlags{radiation_on(), o_.bl_pbl_physics > 0, o_.cu_physics, o_.use_theta_m == 1});
}

std::vector<const Field2*> PhysicsSuite::fields2() const {
  std::vector<const Field2*> v = {&lwcf, &swcf, &hrang, &uratx, &vratx, &tratx, &ra, &rs, &shall, &lwupt, &lwuptc, &lwdnt, &lwdntc, &lwupb, &lwupbc, &lwdnb, &lwdnbc, &swupt, &swuptc, &swdnt, &swdntc, &swupb, &swupbc, &swdnb, &swdnbc, &swvisdir, &swvisdif, &swnirdir, &swnirdif, &swddir, &swddni, &swddif, &swddirc, &swddnic, &diffuse_frac, &aclwupt, &aclwuptc, &aclwdnt, &aclwdntc, &aclwupb, &aclwupbc, &aclwdnb, &aclwdnbc, &acswupt, &acswuptc, &acswdnt, &acswdntc, &acswupb, &acswupbc, &acswdnb, &acswdnbc, &olr, &tsk, &xland, &xice, &lakemask, &water_depth, &mavail, &tmn, &vegfra, &shdmin, &shdmax, &snoalb, &albbck, &albedo, &emiss, &embck, &znt, &z0, &ust, &ustm, &mol, &rmol, &pblh, &hfx, &qfx, &lh, &grdflx, &qsfc, &flhc, &flqc, &wspd, &br, &fm, &fh, &zol, &regime, &u10, &v10, &t2, &th2, &q2, &psfc, &ch, &snow, &snowh, &snowc, &canwat, &lai, &snotime, &qz0, &rainbl, &raincv, &raincv_a, &raincv_b, &rainncv, &rainshv, &rainc, &rainsh, &pratec, &pratesh, &nca, &htop, &hbot, &cutop, &cubot, &cuppt, &sr, &glw, &swdown, &gsw, &smstav, &smstot, &sfcrunoff, &udrunoff, &acsnom, &acsnow, &snopcx, &potevp, &noahres, &sfcevp, &sfcexc, &achfx, &aclhf, &acgrdflx, &dx2d, &var2d, &uoce, &voce, &ctopo, &ctopo2, &ck, &cka, &cd, &cda, &wstar, &delta, &hol, &thz0, &uz0, &vz0, &akhs, &akms, &ct, &tshltr, &th10, &qshltr, &q10, &pshltr, &u10e, &v10e, &mixht};
  if (o_.sf_surface_physics == 4) {
    v.insert(v.end(), {&xlat_sfc, &xlong_sfc, &tvxy, &tgxy, &canicexy, &canliqxy, &eahxy, &tahxy, &cmxy, &chxy,
                       &fwetxy, &sneqvoxy, &alboldxy, &qsnowxy, &qrainxy, &wslakexy, &zwtxy, &waxy, &wtxy, &lfmassxy,
                       &rtmassxy, &stmassxy, &woodxy, &stblcpxy, &fastcpxy, &xsaixy, &taussxy, &smcwtdxy, &deeprechxy,
                       &rechxy, &grainxy, &gddxy, &qtdrain, &t2mvxy, &t2mbxy, &q2mvxy, &q2mbxy, &tradxy, &neexy,
                       &gppxy, &nppxy, &fvegxy, &runsfxy, &runsbxy, &ecanxy, &edirxy, &etranxy, &fsaxy, &firaxy,
                       &aparxy, &psnxy, &savxy, &sagxy, &rssunxy, &rsshaxy, &bgapxy, &wgapxy, &tgvxy, &tgbxy, &chvxy,
                       &chbxy, &shgxy, &shcxy, &shbxy, &evgxy, &evbxy, &ghvxy, &ghbxy, &irgxy, &ircxy, &irbxy, &trxy,
                       &evcxy, &chleafxy, &chucxy, &chv2xy, &chb2xy, &chstarxy, &acc_ssoilxy, &acc_qinsurxy,
                       &acc_qsevaxy, &eflxbxy, &acc_dwaterxy, &acc_prcpxy, &acc_ecanxy, &acc_etranxy, &acc_edirxy, &acc_glaflwxy});
  }
  if (o_.gwd_opt == 1) v.insert(v.end(), {&dusfcg, &dvsfcg});
  if (o_.sf_sfclay_physics == 5) v.push_back(&cqs);
  if (o_.bl_pbl_physics == 5) v.insert(v.end(), {&maxwidth, &maxmf, &ztop_plume, &excess_h, &excess_q});
  if (o_.fractional_seaice == 1 || o_.sst_update == 1) v.insert(v.end(), {&sst, &tsk_save});
  if (o_.sst_update == 1) v.push_back(&xicem);
  return v;
}

std::vector<const Field3*> PhysicsSuite::fields3() const {
  std::vector<const Field3*> v = {&rthratenlw, &rthratensw, &rthratenlwc, &rthratenswc, &o3rad, &tslb, &smois, &sh2o, &smcrel, &exch_h, &exch_m, &tke_pbl, &el_pbl, &w0avg, &cldfra_dp, &cldfra_sh, &qc_cu, &qi_cu, &gd_cloud, &gd_cloud2, &gd_cloud_a, &gd_cloud2_a, &gd_cldfr, &pt.rthraten, &pt.rublten, &pt.rvblten, &pt.rthblten, &pt.rqvblten, &pt.rqcblten, &pt.rqiblten, &pt.rucuten, &pt.rvcuten, &pt.rthcuten, &pt.rqvcuten, &pt.rqccuten, &pt.rqrcuten, &pt.rqicuten, &pt.rqscuten, &pt.rthften, &pt.rqvften, &prep.th_phy, &prep.th_phy_m_t0, &prep.p_phy, &prep.pi_phy, &prep.t_phy, &prep.rho, &prep.u_phy, &prep.v_phy, &prep.z, &prep.p_hyd, &prep.z_at_w, &prep.dz8w, &prep.p8w, &prep.t8w, &prep.p_hyd_w};
  if (o_.sf_surface_physics == 4) {
    v.insert(v.end(), {&tsnoxy, &snicexy, &snliqxy, &zsnsoxy, &albsoildirxy, &albsoildifxy, &smoiseq, &acc_etranixy});
  }
  if (o_.gwd_opt == 1) v.insert(v.end(), {&dtaux3d, &dtauy3d});
  if (o_.bl_pbl_physics == 5) v.insert(v.end(), {&qke, &tsq, &qsq, &cov, &sh3d, &sm3d, &qc_bl, &qi_bl, &cldfra_bl});
  return v;
}

}  // namespace cppwrf::phys::driver
