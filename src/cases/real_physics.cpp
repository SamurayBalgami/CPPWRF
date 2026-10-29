// Column physics of a real-data run (P1 suite: sfclayrev + Noah + YSU): the suite is created
// and its initial surface state filled as WRF's start_domain does for a non-restart run:
//   - the wrfinput fields (TSK, TMN, TSLB, SMOIS, SEAICE, LU_INDEX, ISLTYP, VEGFRA, ...),
//   - phy_init (module_physics_init.F): UST = 1e-4, MOL = PBLH = HFX = QFX = RAINBL = ACSNOW = 0,
//     IVGTYP = nint(LU_INDEX), compute_2d_dx_area (DX2D = dx), landuse_init (ALBBCK, ALBEDO,
//     THC, Z0, ZNT, EMBCK, EMISS, MAVAIL, XLAND from LANDUSE.TBL for the season of JULDAY),
//   - bl_init: ACHFX = ACLHF = ACGRDFLX = SFCEVP = 0, LSMINIT (module_sf_noahdrv.F): SH2O from
//     SMOIS/TSLB (FRH2O for frozen layers), SNOWH = 0.005*SNOW without input SNOWH, CANWAT = 0,
//     SNOALB = MAXALB/100 when rdmaxalb = false,
//   - start_em: CTOPO = CTOPO2 = 1 (topo_wind = 0).
// Every other surface field starts at zero (WRF allocates grid fields zeroed).
#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cases/model.hpp"
#include "phys/driver/landuse_tables.hpp"
#include "phys/driver/suite.hpp"
#include "phys/lsm_noah/noah.hpp"
#include "wrf_io/namelist.hpp"
#ifdef CPPWRF_HAVE_NETCDF
#include <netcdf.h>
#endif

namespace cppwrf::cases {

bool wants_column_physics(const io::Config& cfg) {
  return cfg.get_int("phys.sf_sfclay_physics", 0) != 0 || cfg.get_int("phys.sf_surface_physics", 0) != 0 ||
         cfg.get_int("phys.bl_pbl_physics", 0) != 0 || cfg.get_int("phys.cu_physics", 0) != 0;
}

LanduseInit landuse_init(const std::string& table_dir, const std::string& mminlu, int julday, double cen_lat,
                         int iswater, int isice, const std::vector<double>& lu, const std::vector<double>& xice,
                         const std::vector<double>& snowc, int fractional_seaice, bool usemonalb,
                         const std::vector<double>& albbck_in, const std::vector<double>& snoalb_in) {
  phys::driver::LanduseTables lt;
  if (!phys::driver::read_landuse_table(table_dir + "/LANDUSE.TBL", mminlu, lt))
    throw std::runtime_error("real physics: LANDUSE.TBL lacks " + mminlu);
  int isn = 1;
  if (julday < 105 || julday > 288) isn = 2;
  if (cen_lat < 0.0) isn = 3 - isn;
  if (lt.n_seasons == 1) isn = 1;
  const double xice_threshold = fractional_seaice == 0 ? 0.5 : 0.02;
  const std::size_t n2 = lu.size();
  LanduseInit r;
  std::vector<double>&albbck = r.albbck, &albedo = r.albedo, &z0 = r.z0, &embck = r.embck, &emiss = r.emiss,
                     &mavail = r.mavail, &xland = r.xland;
  for (auto* v : {&albbck, &albedo, &z0, &embck, &emiss, &mavail, &xland, &r.thc}) v->assign(n2, 0.0);
  r.isn = isn;
  r.lucats = lt.n_cats;
  r.luseas = lt.n_seasons;
  {  // landuse_init's saved table (max_cats = 100, max_seas = 12)
    r.lu_state.assign(7501, 0.0);
    std::size_t curs = 0;
    for (int cats = 1; cats <= 100; ++cats) {
      const bool in = cats < phys::driver::kMaxLanduseCats && cats <= lt.n_cats;
      // SCFX(LC) is read again for every season: the last season's value stays
      r.lu_state[curs++] = in ? double(lt.categories[std::max(0, lt.n_seasons - 1)][cats].scfx) : 0.0;
      for (int seas = 1; seas <= 12; ++seas) {
        const bool ok = in && seas <= lt.n_seasons;
        const phys::driver::LanduseCategory* c = ok ? &lt.categories[seas - 1][cats] : nullptr;
        for (real x : {c ? c->albd : real(0), c ? c->slmo : real(0), c ? c->sfem : real(0), c ? c->sfz0 : real(0),
                       c ? c->sfhc : real(0), c ? c->therin : real(0)})
          r.lu_state[curs++] = double(x);
      }
    }
  }
  for (std::size_t p = 0; p < n2; ++p) {
    int is = int(std::lround(lu[p]));
    if (is < 0 || is > lt.n_cats) throw std::runtime_error("real physics: LANDUSE OUTSIDE RANGE");
    if (is == 0) is = iswater;
    const auto& c = lt.categories[isn - 1][is];
    albbck[p] = usemonalb ? albbck_in[p] : c.albd / 100.;
    albedo[p] = albbck[p];
    // SCFX(IS) comes from lu_state: the last season read (see above), not season ISN
    if (snowc[p] > 0.5)
      albedo[p] = usemonalb ? snoalb_in[p] : albbck[p] * (1. + lt.categories[std::max(0, lt.n_seasons - 1)][is].scfx);
    z0[p] = c.sfz0 / 100.;
    embck[p] = c.sfem;
    emiss[p] = embck[p];
    mavail[p] = c.slmo;
    r.thc[p] = double(c.therin / real(100.));
    xland[p] = is != iswater ? 1.0 : 2.0;
    if (xice[p] >= xice_threshold) {
      const auto& ci = lt.categories[isn - 1][isice];
      xland[p] = 1.0;
      albbck[p] = ci.albd / 100.;
      embck[p] = ci.sfem;
      if (fractional_seaice == 1) {  // open water: albedo 0.08, emissivity 0.98
        albedo[p] = (xice[p] * albbck[p]) + ((1.0 - xice[p]) * 0.08);
        emiss[p] = (xice[p] * embck[p]) + ((1.0 - xice[p]) * 0.98);
      } else {
        albedo[p] = albbck[p];
        emiss[p] = embck[p];
      }
      z0[p] = ci.sfz0 / 100.;
      mavail[p] = ci.slmo;
      r.thc[p] = double(ci.therin / real(100.));
    }
  }
  return r;
}

#ifndef CPPWRF_HAVE_NETCDF
void init_real_physics(const io::Config&, Model&) { throw std::runtime_error("column physics needs NetCDF"); }
void read_wrflowinp(const std::string&, const std::string&, phys::driver::PhysicsSuite&) {
  throw std::runtime_error("sst_update needs NetCDF");
}
#else
namespace {
class NcIn {
 public:
  explicit NcIn(const std::string& path) : path_(path) {
    if (nc_open(path.c_str(), NC_NOWRITE, &id_) != NC_NOERR) throw std::runtime_error("real physics: cannot open " + path);
  }
  ~NcIn() { nc_close(id_); }
  NcIn(const NcIn&) = delete;
  NcIn& operator=(const NcIn&) = delete;
  bool has(const char* name) const {
    int v = 0;
    return nc_inq_varid(id_, name, &v) == NC_NOERR;
  }
  std::vector<double> var(const char* name, std::size_t n) const {
    int v = 0;
    std::vector<double> out(n);
    if (nc_inq_varid(id_, name, &v) != NC_NOERR || nc_get_var_double(id_, v, out.data()) != NC_NOERR)
      throw std::runtime_error(std::string("real physics: ") + path_ + " lacks " + name);
    return out;
  }
  std::string att_str(const char* name, const std::string& def) const {
    std::size_t len = 0;
    if (nc_inq_attlen(id_, NC_GLOBAL, name, &len) != NC_NOERR) return def;
    std::string s(len, '\0');
    if (nc_get_att_text(id_, NC_GLOBAL, name, s.data()) != NC_NOERR) return def;
    while (!s.empty() && (s.back() == '\0' || s.back() == ' ')) s.pop_back();
    return s;
  }
  double att_num(const char* name, double def) const {
    double x = 0;
    return nc_get_att_double(id_, NC_GLOBAL, name, &x) == NC_NOERR ? x : def;
  }

 private:
  std::string path_;
  int id_ = -1;
};
}  // namespace

// med_auxinput_in for auxinput4 (wrflowinp_d0N): the record at `date` overwrites the fields of the stream
// that the file holds (SST, SEAICE, VEGFRA, ALBBCK, LAI, UOCE, VOCE)
void read_wrflowinp(const std::string& path, const std::string& date, phys::driver::PhysicsSuite& su) {
  int id = -1;
  if (nc_open(path.c_str(), NC_NOWRITE, &id) != NC_NOERR) throw std::runtime_error("sst_update: cannot open " + path);
  int vt = -1, dims[2];
  std::size_t nt = 0, sl = 0;
  if (nc_inq_varid(id, "Times", &vt) != NC_NOERR) throw std::runtime_error("sst_update: " + path + " has no Times");
  nc_inq_vardimid(id, vt, dims);
  nc_inq_dimlen(id, dims[0], &nt);
  nc_inq_dimlen(id, dims[1], &sl);
  long rec = -1;
  for (std::size_t t = 0; t < nt && rec < 0; ++t) {
    std::string s(sl, '\0');
    const std::size_t st[2] = {t, 0}, ct[2] = {1, sl};
    nc_get_vara_text(id, vt, st, ct, s.data());
    if (s.compare(0, date.size(), date) == 0) rec = long(t);
  }
  if (rec < 0) {
    nc_close(id);
    throw std::runtime_error("sst_update: " + path + " has no record for " + date);
  }
  const Layout& L = su.sst.layout();
  const int nx = L.nx, ny = L.ny;
  const std::pair<const char*, Field2*> fields[] = {
      {"SST", &su.sst}, {"SEAICE", &su.xice}, {"VEGFRA", &su.vegfra}, {"ALBBCK", &su.albbck},
      {"LAI", &su.lai}, {"UOCE", &su.uoce}, {"VOCE", &su.voce}};
  std::vector<double> buf(std::size_t(nx) * ny);
  for (const auto& [name, f] : fields) {
    int v = -1;
    if (nc_inq_varid(id, name, &v) != NC_NOERR) continue;
    const std::size_t st[3] = {std::size_t(rec), 0, 0}, ct[3] = {1, std::size_t(ny), std::size_t(nx)};
    if (nc_get_vara_double(id, v, st, ct, buf.data()) != NC_NOERR)
      throw std::runtime_error(std::string("sst_update: cannot read ") + name + " of " + path);
    f->fill_host([&](int i, int j, int) { return buf[i + std::size_t(nx) * j]; });
  }
  nc_close(id);
}

void init_real_physics(const io::Config& cfg, Model& m) {
  using phys::driver::PhysicsSuite;
  using phys::driver::SuiteOptions;
  const Grid& g = m.g;
  const int nx = g.nx, ny = g.ny;
  const std::string dir = cfg.get_string("real.dir");
  const NcIn in(dir + "/" + cfg.get_string("real.wrfinput", "wrfinput_d01"));

  SuiteOptions o;
  o.sf_sfclay_physics = int(cfg.get_int("phys.sf_sfclay_physics", 0));
  o.sf_surface_physics = int(cfg.get_int("phys.sf_surface_physics", 0));
  o.bl_pbl_physics = int(cfg.get_int("phys.bl_pbl_physics", 0));
  o.icloud_bl = int(cfg.get_int("phys.icloud_bl", 1));
  o.mynn.cloudpdf = int(cfg.get_int("phys.bl_mynn_cloudpdf", 2));  // MYNN-EDMF options (WRF's Registry names)
  o.mynn.mixlength = int(cfg.get_int("phys.bl_mynn_mixlength", 1));
  o.mynn.edmf = int(cfg.get_int("phys.bl_mynn_edmf", 1));
  o.mynn.edmf_mom = int(cfg.get_int("phys.bl_mynn_edmf_mom", 1));
  o.mynn.edmf_tke = int(cfg.get_int("phys.bl_mynn_edmf_tke", 0));
  o.mynn.cloudmix = int(cfg.get_int("phys.bl_mynn_cloudmix", 1));
  o.mynn.mixqt = int(cfg.get_int("phys.bl_mynn_mixqt", 0));
  o.mynn.ess = int(cfg.get_int("phys.bl_mynn_ess", 1));
  o.mynn.tkeadvect = cfg.get_bool("phys.bl_mynn_tkeadvect", false);
  o.mynn.closure = real(cfg.get_double("phys.bl_mynn_closure", 2.6));
  o.bl_mynn_mixscalars = int(cfg.get_int("phys.bl_mynn_mixscalars", 0));
  o.bl_mynn_output = int(cfg.get_int("phys.bl_mynn_output", 0));
  o.tke_budget = int(cfg.get_int("phys.tke_budget", 0));
  o.grav_settling = int(cfg.get_int("phys.grav_settling", 0));
  o.cu_physics = int(cfg.get_int("phys.cu_physics", 0));
  o.num_soil_layers = int(cfg.get_int("phys.num_soil_layers", 4));
  o.bldt = cfg.get_double("phys.bldt", 0.0);
  o.cudt = cfg.get_double("phys.cudt", 0.0);
  // check_a_mundo (module_check_a_mundo.F:2156-2174): the adaptive time step runs BMJ, the SAS schemes and Tiedtke
  // every step (cudt = 0)
  if (cfg.get_bool("time.adaptive", false) && o.cu_physics == 6) o.cudt = 0.0;
  {
    const std::string nl_path = cfg.get_string("output.wrf_namelist", dir + "/namelist.input");
    const double nl_avg = std::filesystem::exists(nl_path)
                              ? wrf_io::Namelist::read(nl_path).get_real("physics", "convtrans_avglen_m", 30.0)
                              : 30.0;
    o.convtrans_avglen_m = cfg.get_double("phys.convtrans_avglen_m", nl_avg);
  }
  o.kfeta_trigger = int(cfg.get_int("phys.kfeta_trigger", 1));
  o.ishallow = int(cfg.get_int("phys.ishallow", 0));          // Grell-Freitas shallow convection
  o.gf_ichoice = int(cfg.get_int("phys.clos_choice", 0));     // Grell-Freitas closure
  o.warm_rain = cfg.get_int("phys.mp_physics", 0) == 1;  // mp_init: true only for KESSLERSCHEME
  o.isfflx = int(cfg.get_int("phys.isfflx", 1));
  o.opt_thcnd = int(cfg.get_int("phys.opt_thcnd", 1));
  o.rdlai2d = cfg.get_bool("phys.rdlai2d", false);      // Noah: LAI from the input instead of the table
  o.usemonalb = cfg.get_bool("phys.usemonalb", false);  // landuse_init and Noah: ALBBCK from the input
  o.ysu_topdown_pblmix = int(cfg.get_int("phys.ysu_topdown_pblmix", 1));
  o.kim_tofd = cfg.get_bool("phys.kim_tofd", false);
  o.isftcflx = int(cfg.get_int("phys.isftcflx", 0));  // sfclayrev: z0t / z0q and roughness over water
  o.iz0tlnd = int(cfg.get_int("phys.iz0tlnd", 0));    // sfclayrev: z0t over land (czil)
  o.sf_mynn_sfcflux_land = int(cfg.get_int("phys.sf_mynn_sfcflux_land", 0));    // MYNN surface layer
  o.sf_mynn_sfcflux_water = int(cfg.get_int("phys.sf_mynn_sfcflux_water", 1));
  o.fractional_seaice = int(cfg.get_int("phys.fractional_seaice", 0));
  o.tice2tsk_if2cold = cfg.get_bool("phys.tice2tsk_if2cold", false);
  o.sst_update = int(cfg.get_int("phys.sst_update", 0));
  o.gwd_opt = int(cfg.get_int("dyn.gwd_opt", 0));
  o.gwd_dx_factor = cfg.get_double("dyn.gwd_dx_factor", 2.0);
  o.gwd_if_nonhyd = cfg.get_bool("dyn.gwd_if_nonhyd", true);
  o.tofd_factor = cfg.get_double("phys.tofd_factor", 0.003);
  o.mminlu = in.att_str("MMINLU", o.mminlu);
  o.iswater = int(in.att_num("ISWATER", o.iswater));
  o.isice = int(in.att_num("ISICE", o.isice));
  o.isurban = int(in.att_num("ISURBAN", o.isurban));
  o.table_dir = cfg.get_string("phys.table_dir", dir);
  const bool rdmaxalb = cfg.get_bool("phys.rdmaxalb", true);
  {  // radiation (radiation_driver; RRTMG tables in phys.rrtmg_tables or data/rrtmg, ozone and CAMtr from phys.rrtmg_tables)
    auto& r = o.rad;
    r.ra_lw_physics = int(cfg.get_int("phys.ra_lw_physics", 0));
    r.ra_sw_physics = int(cfg.get_int("phys.ra_sw_physics", 0));
    r.radt = cfg.get_double("phys.radt", 0.0);
    r.icloud = int(cfg.get_int("phys.icloud", 1));
    r.cldovrlp = int(cfg.get_int("phys.cldovrlp", 2));
    r.idcor = int(cfg.get_int("phys.idcor", 0));
    r.o3input = int(cfg.get_int("phys.o3input", 2));
    r.oznini = !cfg.get_bool("real.nested", false);
    r.ghg_input = int(cfg.get_int("phys.ghg_input", 1));
    r.aer_opt = int(cfg.get_int("phys.aer_opt", 0));
    r.julday = int(in.att_num("JULDAY", 1));
    r.julyr = int(in.att_num("JULYR", 0));
    r.gmt = in.att_num("GMT", 0.0);
    r.p_top = g.p_top;
    r.table_dir = cfg.get_string("phys.rrtmg_tables", o.table_dir);
  }
  if (o.sf_surface_physics == 4) {
    auto& nm = o.noahmp;
    nm.dveg = int(cfg.get_int("phys.dveg", nm.dveg));
    nm.opt_crs = int(cfg.get_int("phys.opt_crs", nm.opt_crs));
    nm.opt_btr = int(cfg.get_int("phys.opt_btr", nm.opt_btr));
    nm.opt_runsrf = int(cfg.get_int("phys.opt_runsrf", nm.opt_runsrf));
    nm.opt_runsub = int(cfg.get_int("phys.opt_runsub", nm.opt_runsub));
    nm.opt_sfc = int(cfg.get_int("phys.opt_sfc", nm.opt_sfc));
    nm.opt_frz = int(cfg.get_int("phys.opt_frz", nm.opt_frz));
    nm.opt_inf = int(cfg.get_int("phys.opt_inf", nm.opt_inf));
    nm.opt_rad = int(cfg.get_int("phys.opt_rad", nm.opt_rad));
    nm.opt_alb = int(cfg.get_int("phys.opt_alb", nm.opt_alb));
    nm.opt_snf = int(cfg.get_int("phys.opt_snf", nm.opt_snf));
    nm.opt_tbot = int(cfg.get_int("phys.opt_tbot", nm.opt_tbot));
    nm.opt_stc = int(cfg.get_int("phys.opt_stc", nm.opt_stc));
    nm.opt_gla = int(cfg.get_int("phys.opt_gla", nm.opt_gla));
    nm.opt_rsf = int(cfg.get_int("phys.opt_rsf", nm.opt_rsf));
    nm.opt_soil = int(cfg.get_int("phys.opt_soil", nm.opt_soil));
    nm.opt_pedo = int(cfg.get_int("phys.opt_pedo", nm.opt_pedo));
    nm.opt_crop = int(cfg.get_int("phys.opt_crop", nm.opt_crop));
    nm.opt_irr = int(cfg.get_int("phys.opt_irr", nm.opt_irr));
    nm.opt_irrm = int(cfg.get_int("phys.opt_irrm", nm.opt_irrm));
    nm.opt_infdv = int(cfg.get_int("phys.opt_infdv", nm.opt_infdv));
    nm.opt_tdrn = int(cfg.get_int("phys.opt_tdrn", nm.opt_tdrn));
    nm.opt_tksno = int(cfg.get_int("phys.opt_tksno", nm.opt_tksno));
    nm.opt_compact = int(cfg.get_int("phys.opt_compact", nm.opt_compact));
    nm.opt_scf = int(cfg.get_int("phys.opt_scf", nm.opt_scf));
    nm.opt_wetland = int(cfg.get_int("phys.opt_wetland", nm.opt_wetland));
    nm.sf_urban_physics = int(cfg.get_int("phys.sf_urban_physics", 0));
    nm.iz0tlnd = o.iz0tlnd;
    nm.soiltstep = cfg.get_double("phys.soiltstep", 0.0);
  }
  m.suite = std::make_unique<PhysicsSuite>(g, o);
  PhysicsSuite& su = *m.suite;

  const std::size_t n2 = std::size_t(nx) * ny;
  const int ns = o.num_soil_layers;
  auto put2 = [&](Field2& f, const std::vector<double>& a) { f.fill_host([&](int i, int j, int) { return a[i + std::size_t(nx) * j]; }); };
  auto put_const = [&](Field2& f, double c) { f.fill_host([=](int, int, int) { return c; }); };
  auto rd2 = [&](const char* name) { return in.var(name, n2); };

  // ---- wrfinput fields
  const std::vector<double> lu = rd2("LU_INDEX"), xice = rd2("SEAICE"), isl = rd2("ISLTYP");
  const std::vector<double> tslb = in.var("TSLB", n2 * ns), smois = in.var("SMOIS", n2 * ns);
  std::vector<double> snow = rd2("SNOW"), snowc = rd2("SNOWC"), snoalb = rd2("SNOALB");
  put2(su.tsk, rd2("TSK")); put2(su.tmn, rd2("TMN")); put2(su.xice, xice);
  put2(su.vegfra, rd2("VEGFRA")); put2(su.shdmax, rd2("SHDMAX")); put2(su.shdmin, rd2("SHDMIN"));
  put2(su.snow, snow); put2(su.snowc, snowc); put2(su.lai, rd2("LAI"));
  put2(su.lakemask, in.has("LAKEMASK") ? rd2("LAKEMASK") : std::vector<double>(n2, 0.0));
  if (in.has("VAR")) put2(su.var2d, rd2("VAR"));  // subgrid orography std. dev. (kim_tofd; history VAR)
  if (o.gwd_opt == 1) {  // package gwd_used_1
    const std::pair<const char*, Field2*> gw[] = {{"CON", &su.oc12d}, {"OA1", &su.oa1}, {"OA2", &su.oa2}, {"OA3", &su.oa3},
                                                  {"OA4", &su.oa4}, {"OL1", &su.ol1}, {"OL2", &su.ol2}, {"OL3", &su.ol3},
                                                  {"OL4", &su.ol4}, {"ELVMAX", &su.elvmax}};
    for (const auto& [n, f] : gw) {
      if (!in.has(n)) throw std::runtime_error(std::string("gwd_opt = 1: wrfinput lacks ") + n);
      put2(*f, rd2(n));
    }
  }
  if (o.fractional_seaice == 1 || o.sst_update == 1) put2(su.sst, rd2("SST"));
  if (o.fractional_seaice == 1) put2(su.tsk_save, rd2("TSK"));  // phy_init: TSK_SAVE = TSK
  if (o.sst_update == 1) put2(su.xicem, xice);                    // landuse_init: XICEM = XICE
  put2(su.water_depth, in.has("WATER_DEPTH") ? rd2("WATER_DEPTH") : std::vector<double>(n2, 0.0));
  const std::pair<const char*, Field2*> diag[] = {{"T2", &su.t2}, {"TH2", &su.th2}, {"Q2", &su.q2},
                                                   {"U10", &su.u10}, {"V10", &su.v10}, {"PSFC", &su.psfc}};
  for (const auto& [n, f] : diag)
    if (in.has(n)) put2(*f, rd2(n));
  if (o.sf_sfclay_physics == 2) {  // phy_init (module_physics_init.F): u10e = u10, v10e = v10 for MYJSFCSCHEME
    if (in.has("U10")) put2(su.u10e, rd2("U10"));
    if (in.has("V10")) put2(su.v10e, rd2("V10"));
  }
  const std::vector<double> dzs = in.var("DZS", std::size_t(ns));
  for (int k = 0; k < ns; ++k) su.dzs[k] = dzs[k];

  // ---- phy_init / myjsfcinit
  put_const(su.ust, o.sf_sfclay_physics == 2 ? 0.1 : 0.0001);
  put_const(su.dx2d, g.dx);
  put_const(su.ctopo, 1.0); put_const(su.ctopo2, 1.0);
  {
    std::vector<int> h(su.ivgtyp.layout().total, 0), hs(su.isltyp.layout().total, 0);
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        h[su.ivgtyp.layout().index(i, j, 0)] = int(std::lround(lu[i + std::size_t(nx) * j]));
        hs[su.isltyp.layout().index(i, j, 0)] = int(std::lround(isl[i + std::size_t(nx) * j]));
      }
    su.ivgtyp.upload(h);
    su.isltyp.upload(hs);
  }
  // landuse_init
  const LanduseInit li = landuse_init(o.table_dir, o.mminlu, int(in.att_num("JULDAY", 0)), in.att_num("CEN_LAT", 0.0),
                                      o.iswater, o.isice, lu, xice, snowc, o.fractional_seaice, o.usemonalb,
                                      o.usemonalb ? rd2("ALBBCK") : std::vector<double>(),
                                      o.usemonalb ? snoalb : std::vector<double>());
  const std::vector<double>&albbck = li.albbck, &albedo = li.albedo, &z0 = li.z0, &embck = li.embck,
                            &emiss = li.emiss, &mavail = li.mavail, &xland = li.xland;
  put2(su.albbck, albbck); put2(su.albedo, albedo); put2(su.z0, z0); put2(su.znt, z0);
  put2(su.embck, embck); put2(su.emiss, emiss); put2(su.mavail, mavail); put2(su.xland, xland);

  auto put_soil = [&](Field3& f, const std::vector<double>& a) {
    f.fill_host([&](int i, int j, int k) { return a[i + std::size_t(nx) * j + n2 * k]; });
  };

  if (o.sf_surface_physics == 2) {
    // ---- bl_init / LSMINIT (Noah)
    phys::noah::NoahTables nt;
    if (!phys::noah::read_noah_tables(o.table_dir, o.mminlu, o.mminsl, nt))
      throw std::runtime_error("real physics: cannot read the Noah tables from " + o.table_dir);
    if (!rdmaxalb)
      for (std::size_t p = 0; p < n2; ++p) snoalb[p] = nt.maxalb[int(std::lround(lu[p])) - 1] * 0.01;
    put2(su.snoalb, snoalb);
    std::vector<double> sh2o(n2 * ns);
    for (std::size_t p = 0; p < n2; ++p) {
      const int st = int(std::lround(isl[p]));
      if (st < 1) throw std::runtime_error("real physics: lsminit: out of range value of ISLTYP");
      real bx = nt.bb[st - 1];
      const real smcmax = nt.maxsmc[st - 1], psisat = nt.satpsi[st - 1];
      for (int k = 0; k < ns; ++k) {
        const std::size_t q = p + n2 * k;
        if (bx > 0 && smcmax > 0 && psisat > 0) {
          if (tslb[q] < 273.149) {
            bx = nt.bb[st - 1];
            if (bx > real(5.5)) bx = real(5.5);
            const real t = real(tslb[q]);
            real fk = std::pow((real(3.335e5) / (real(9.81) * (-psisat))) * ((t - real(273.15)) / t), real(-1) / bx) * smcmax;
            if (fk < real(0.02)) fk = real(0.02);
            real s = std::min(fk, real(smois[q]));
            real free_w = 0;
            phys::noah::NoahKernel::frh2o(free_w, t, real(smois[q]), s, smcmax, bx, psisat);
            sh2o[q] = free_w;
          } else {
            sh2o[q] = smois[q];
          }
        } else {
          sh2o[q] = smois[q];
        }
      }
    }
    put_soil(su.tslb, tslb); put_soil(su.smois, smois); put_soil(su.sh2o, sh2o);
    if (in.has("SMCREL")) put_soil(su.smcrel, in.var("SMCREL", n2 * ns));
    // phy_init passes FNDSNOWH = .true. to LSMINIT: SNOWH is never derived from SNOW (0 when not in the input)
    if (in.has("SNOWH")) put2(su.snowh, rd2("SNOWH"));
  } else if (o.sf_surface_physics == 4) {
    // ---- bl_init / NoahmpWRFinit (Noah-MP)
    put2(su.snoalb, snoalb);
    for (std::size_t p = 0; p < n2; ++p)  // NoahmpInitMain: stop on any ISLTYP < 1
      if (std::lround(isl[p]) < 1) throw std::runtime_error("real physics: lsminit: out of range ISLTYP");
    put2(su.xlat_sfc, m.xlat);
    put2(su.xlong_sfc, m.xlong);
    put_soil(su.tslb, tslb);
    put_soil(su.smois, smois);
    if (in.has("SH2O")) put_soil(su.sh2o, in.var("SH2O", n2 * ns));
    if (in.has("SNOWH")) put2(su.snowh, rd2("SNOWH"));
    su.init_noahmp(in.has("SNOWH"));
  }
  // CANWAT = 0 (lsminit / NoahmpInitMain), the accumulations and fluxes stay zero
  su.init_radiation(m.xlat, m.xlong);
  if (auto* pc = m.rk->physics(); pc && pc->params().microp_re)
    su.set_effective_radii(&pc->re_cloud, &pc->re_ice, &pc->re_snow);
  if (auto* pc = m.rk->physics(); pc && pc->params().mp_physics != 0) {
    const bool has_ice_precip = pc->params().mp_physics == 6 || pc->params().mp_physics == 8;
    su.set_microphysics_precip(&pc->rainncv, has_ice_precip ? &pc->sr : nullptr,
                               has_ice_precip ? &pc->snowncv : nullptr,
                               has_ice_precip ? &pc->graupelncv : nullptr);
  }
  m.rk->set_column_physics(m.suite.get());
}
#endif

}  // namespace cppwrf::cases
