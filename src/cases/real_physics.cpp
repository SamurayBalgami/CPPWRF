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
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cases/model.hpp"
#include "phys/driver/landuse_tables.hpp"
#include "phys/driver/suite.hpp"
#include "phys/lsm_noah/noah.hpp"
#ifdef CPPWRF_HAVE_NETCDF
#include <netcdf.h>
#endif

namespace cppwrf::cases {

bool wants_column_physics(const io::Config& cfg) {
  return cfg.get_int("phys.sf_sfclay_physics", 0) != 0 || cfg.get_int("phys.sf_surface_physics", 0) != 0 ||
         cfg.get_int("phys.bl_pbl_physics", 0) != 0 || cfg.get_int("phys.cu_physics", 0) != 0;
}

#ifndef CPPWRF_HAVE_NETCDF
void init_real_physics(const io::Config&, Model&) { throw std::runtime_error("column physics needs NetCDF"); }
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
  o.cu_physics = int(cfg.get_int("phys.cu_physics", 0));
  o.num_soil_layers = int(cfg.get_int("phys.num_soil_layers", 4));
  o.bldt = cfg.get_double("phys.bldt", 0.0);
  o.cudt = cfg.get_double("phys.cudt", 0.0);
  o.kfeta_trigger = int(cfg.get_int("phys.kfeta_trigger", 1));
  o.warm_rain = cfg.get_int("phys.mp_physics", 0) == 1;  // mp_init: true only for KESSLERSCHEME
  o.isfflx = int(cfg.get_int("phys.isfflx", 1));
  o.ysu_topdown_pblmix = int(cfg.get_int("phys.ysu_topdown_pblmix", 1));
  o.mminlu = in.att_str("MMINLU", o.mminlu);
  o.iswater = int(in.att_num("ISWATER", o.iswater));
  o.isice = int(in.att_num("ISICE", o.isice));
  o.isurban = int(in.att_num("ISURBAN", o.isurban));
  o.table_dir = cfg.get_string("phys.table_dir", dir);
  const bool rdmaxalb = cfg.get_bool("phys.rdmaxalb", true);
  {  // radiation (radiation_driver; RRTMG tables dumped from WRF's rrtmg_lw/swinit in phys.rrtmg_tables)
    auto& r = o.rad;
    r.ra_lw_physics = int(cfg.get_int("phys.ra_lw_physics", 0));
    r.ra_sw_physics = int(cfg.get_int("phys.ra_sw_physics", 0));
    r.radt = cfg.get_double("phys.radt", 0.0);
    r.icloud = int(cfg.get_int("phys.icloud", 1));
    r.cldovrlp = int(cfg.get_int("phys.cldovrlp", 2));
    r.idcor = int(cfg.get_int("phys.idcor", 0));
    r.o3input = int(cfg.get_int("phys.o3input", 2));
    r.ghg_input = int(cfg.get_int("phys.ghg_input", 1));
    r.aer_opt = int(cfg.get_int("phys.aer_opt", 0));
    r.julday = int(in.att_num("JULDAY", 1));
    r.julyr = int(in.att_num("JULYR", 0));
    r.gmt = in.att_num("GMT", 0.0);
    r.p_top = g.p_top;
    r.table_dir = cfg.get_string("phys.rrtmg_tables", o.table_dir);
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
  put2(su.water_depth, in.has("WATER_DEPTH") ? rd2("WATER_DEPTH") : std::vector<double>(n2, 0.0));
  const std::pair<const char*, Field2*> diag[] = {{"T2", &su.t2}, {"TH2", &su.th2}, {"Q2", &su.q2},
                                                   {"U10", &su.u10}, {"V10", &su.v10}, {"PSFC", &su.psfc}};
  for (const auto& [n, f] : diag)
    if (in.has(n)) put2(*f, rd2(n));
  const std::vector<double> dzs = in.var("DZS", std::size_t(ns));
  for (int k = 0; k < ns; ++k) su.dzs[k] = dzs[k];

  // ---- phy_init
  put_const(su.ust, 0.0001);
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
  phys::driver::LanduseTables lt;
  if (!phys::driver::read_landuse_table(o.table_dir + "/LANDUSE.TBL", o.mminlu, lt))
    throw std::runtime_error("real physics: LANDUSE.TBL lacks " + o.mminlu);
  const int julday = int(in.att_num("JULDAY", 0));
  const double cen_lat = in.att_num("CEN_LAT", 0.0);
  int isn = 1;
  if (julday < 105 || julday > 288) isn = 2;
  if (cen_lat < 0.0) isn = 3 - isn;
  if (lt.n_seasons == 1) isn = 1;
  const double xice_threshold = 0.5;  // fractional_seaice = 0
  std::vector<double> albbck(n2), albedo(n2), z0(n2), embck(n2), emiss(n2), mavail(n2), xland(n2);
  for (std::size_t p = 0; p < n2; ++p) {
    int is = int(std::lround(lu[p]));
    if (is < 0 || is > lt.n_cats) throw std::runtime_error("real physics: LANDUSE OUTSIDE RANGE");
    if (is == 0) is = o.iswater;
    const auto& c = lt.categories[isn - 1][is];
    albbck[p] = c.albd / 100.;
    albedo[p] = albbck[p];
    if (snowc[p] > 0.5) albedo[p] = albbck[p] * (1. + c.scfx);
    z0[p] = c.sfz0 / 100.;
    embck[p] = c.sfem;
    emiss[p] = embck[p];
    mavail[p] = c.slmo;
    xland[p] = is != o.iswater ? 1.0 : 2.0;
    if (xice[p] >= xice_threshold) {
      const auto& ci = lt.categories[isn - 1][o.isice];
      xland[p] = 1.0;
      albbck[p] = ci.albd / 100.;
      embck[p] = ci.sfem;
      albedo[p] = albbck[p];
      emiss[p] = embck[p];
      z0[p] = ci.sfz0 / 100.;
      mavail[p] = ci.slmo;
    }
  }
  put2(su.albbck, albbck); put2(su.albedo, albedo); put2(su.z0, z0); put2(su.znt, z0);
  put2(su.embck, embck); put2(su.emiss, emiss); put2(su.mavail, mavail); put2(su.xland, xland);

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
  auto put_soil = [&](Field3& f, const std::vector<double>& a) {
    f.fill_host([&](int i, int j, int k) { return a[i + std::size_t(nx) * j + n2 * k]; });
  };
  put_soil(su.tslb, tslb); put_soil(su.smois, smois); put_soil(su.sh2o, sh2o);
  if (in.has("SMCREL")) put_soil(su.smcrel, in.var("SMCREL", n2 * ns));
  if (in.has("SNOWH")) put2(su.snowh, rd2("SNOWH"));
  else {
    for (double& x : snow) x *= 0.005;
    put2(su.snowh, snow);
  }
  // CANWAT = 0 (lsminit), the accumulations and fluxes stay zero
  su.init_radiation(m.xlat, m.xlong);
  m.rk->set_column_physics(m.suite.get());
}
#endif

}  // namespace cppwrf::cases
