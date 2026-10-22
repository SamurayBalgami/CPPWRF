#ifdef CPPWRF_HAVE_NETCDF
#include "cases/wrf_history.hpp"
#include "core/strings.hpp"

#include <netcdf.h>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <set>
#include <stdexcept>

#include "core/constants.hpp"
#include "dyn/physics_coupling.hpp"
#include "dyn/rk3.hpp"
#include "phys/driver/suite.hpp"
#include "wrf_io/global_atts.hpp"
#include "wrf_io/namelist.hpp"
#include "wrf_io/stream_file.hpp"
#include "wrf_io/wrf_schema.hpp"

#ifndef CPPWRF_DATA_DIR
#define CPPWRF_DATA_DIR "data"
#endif

namespace cppwrf::cases {

namespace {

std::string str(double v) {
  char b[64];
  std::snprintf(b, sizeof b, "%.17g", v);
  return b;
}

// field (i, j, k) -> NetCDF order (k, j, i) for dims {nk, nj, ni} / {nj, ni} / {nk}
template <class F>
bool copy_field(const F& f, const std::vector<std::size_t>& dims, std::vector<double>& out) {
  const Layout& L = f.layout();
  std::size_t nk = 1, nj = 1, ni = 1;
  if (dims.size() == 3) { nk = dims[0]; nj = dims[1]; ni = dims[2]; }
  else if (dims.size() == 2) { nj = dims[0]; ni = dims[1]; }
  else return false;
  // a mass-level field written on WRF's kms:kme levels (MYJ TKE_PBL, EL_PBL, EXCH_H, ...): the level above
  // kte stays 0 in WRF
  const bool pad = dims.size() == 3 && L.nz > 1 && int(nk) == L.nz + 1;
  // defor12/13/23: WRF's arrays are mass-sized, ours carry the extra staggered row (same indices)
  const bool sub = f.name() == "defor12" || f.name() == "defor13" || f.name() == "defor23";
  if (sub ? (int(ni) > L.nx || int(nj) > L.ny || int(nk) != L.nz)
          : (int(ni) != L.nx || int(nj) != L.ny || (int(nk) != std::max(1, L.nz) && !pad)))
    return false;
  std::vector<real> h;
  f.download(h);
  out.assign(nk * nj * ni, 0.0);
  for (std::size_t k = 0; k < (pad ? nk - 1 : nk); ++k)
    for (std::size_t j = 0; j < nj; ++j)
      for (std::size_t i = 0; i < ni; ++i)
        out[(k * nj + j) * ni + i] = double(h[L.index(int(i), int(j), int(k))]);
  return true;
}

// NetCDF order (k, j, i) -> field (i, j, k), the inverse of copy_field
template <class F>
bool store_field(F& f, const std::vector<std::size_t>& dims, const std::vector<double>& in) {
  const Layout& L = f.layout();
  std::size_t nk = 1, nj = 1, ni = 1;
  if (dims.size() == 3) { nk = dims[0]; nj = dims[1]; ni = dims[2]; }
  else if (dims.size() == 2) { nj = dims[0]; ni = dims[1]; }
  else return false;
  const bool pad = dims.size() == 3 && L.nz > 1 && int(nk) == L.nz + 1;
  if (int(ni) != L.nx || int(nj) != L.ny || (int(nk) != std::max(1, L.nz) && !pad)) return false;
  std::vector<real> h;
  f.download(h);
  for (std::size_t k = 0; k < (pad ? nk - 1 : nk); ++k)
    for (std::size_t j = 0; j < nj; ++j)
      for (std::size_t i = 0; i < ni; ++i) h[L.index(int(i), int(j), int(k))] = real(in[(k * nj + j) * ni + i]);
  copy_to_dev(f.data() - L.offset, h.data(), L.total * sizeof(real));
  return true;
}

// base-state field (horizontally uniform storage broadcast), dims as copy_field
template <class F>
bool copy_base(const F& f, const std::vector<std::size_t>& dims, int nx, int ny, int nzl, std::vector<double>& out) {
  const Layout& L = f.layout();
  std::size_t nk = 1, nj = 1, ni = 1;
  if (dims.size() == 3) { nk = dims[0]; nj = dims[1]; ni = dims[2]; }
  else if (dims.size() == 2) { nj = dims[0]; ni = dims[1]; }
  else return false;
  if (int(ni) != nx || int(nj) != ny || int(nk) != nzl) return false;
  std::vector<real> h;
  f.download(h);
  const bool u = f.uniform();
  out.resize(nk * nj * ni);
  for (std::size_t k = 0; k < nk; ++k)
    for (std::size_t j = 0; j < nj; ++j)
      for (std::size_t i = 0; i < ni; ++i)
        out[(k * nj + j) * ni + i] = double(h[L.index(u ? 0 : int(i), u ? 0 : int(j), int(k))]);
  return true;
}

}  // namespace

WrfHistory::WrfHistory(const io::Config& cfg, const Model& m, std::string out_dir, std::string domain_tag,
                       std::string start_date, std::vector<const Field2*> fields2, std::vector<const Field3*> fields3)
    : m_(m), out_dir_(std::move(out_dir)), domain_tag_(std::move(domain_tag)), start_date_(std::move(start_date)),
      f2_(std::move(fields2)), f3_(std::move(fields3)) {
  const std::string real_dir = cfg.get_string("real.dir", "");
  std::string nl_path = cfg.get_string("output.wrf_namelist", real_dir.empty() ? "" : real_dir + "/namelist.input");
  nl_ = std::make_unique<wrf_io::Namelist>(!nl_path.empty() && std::filesystem::exists(nl_path)
                                               ? wrf_io::Namelist::read(nl_path)
                                               : wrf_io::Namelist::parse(""));
  schema_ = std::make_unique<wrf_io::WrfSchema>(wrf_io::WrfSchema::load(
      cfg.get_string("output.wrf_schema", std::string(CPPWRF_DATA_DIR) + "/wrf_schema/wrf-4.8.0.tsv")));
  // the nest feedback (u= fields) does not depend on the I/O streams: the Registry masks before iofields_filename
  fb_schema_ = std::make_unique<wrf_io::WrfSchema>(*schema_);
  {  // iofields_filename: runtime additions to / removals from the stream masks (relative to the run dir)
    wrf_io::WrfConfig c0(*nl_, schema_->defaults(), domain_tag_.size() == 3 ? std::atoi(domain_tag_.c_str() + 1) : 1);
    std::string io = c0.get("iofields_filename");
    if (!io.empty() && io.front() == '"') io = io.substr(1, io.size() - 2);
    if (!io.empty() && io != "NONE_SPECIFIED") {
      const std::string nl_dir = std::filesystem::path(nl_path).parent_path().string();
      if (io.front() != '/' && !nl_dir.empty()) io = nl_dir + "/" + io;
      if (!schema_->apply_iofields(io) && !c0.logical("ignore_iofields_warning"))
        throw std::runtime_error("wrfout: cannot read iofields_filename " + io);
    }
  }
  const int dom = domain_tag_.size() == 3 ? std::atoi(domain_tag_.c_str() + 1) : 1;
  wc_ = std::make_unique<wrf_io::WrfConfig>(*nl_, schema_->defaults(), dom);
  // the model's options (INI configuration) where namelist.input does not set them
  struct ModelOpts {
    wrf_io::WrfConfig& c;
    const wrf_io::Namelist& nl;
    const std::map<std::string, std::pair<std::string, std::string>>& def;
    int dom;
    void set(const std::string& k, const std::string& v) {
      auto it = def.find(k);
      if (it == def.end() || !nl.raw(it->second.first, k, dom)) c.set(k, v);
    }
  } c{*wc_, *nl_, schema_->defaults(), dom};
  const Grid& g = m_.g;
  if (real_dir.empty() && g.two_d()) jrep_ = 2;
  c.set("e_we", std::to_string(g.nx + 1)); c.set("e_sn", std::to_string(g.ny * jrep_ + 1));
  c.set("e_vert", std::to_string(g.nz + 1)); c.set("dx", str(g.dx)); c.set("dy", str(g.dy));
  c.set("diff_opt", std::to_string(m_.diff_opt)); c.set("km_opt", std::to_string(m_.km_opt));
  c.set("hybrid_opt", std::to_string(m_.hybrid_opt));
  c.set("damp_opt", std::to_string(m_.ap.damp_opt)); c.set("dampcoef", str(double(m_.ap.dampcoef)));
  c.set("zdamp", str(double(m_.ap.zdamp))); c.set("hypsometric_opt", std::to_string(m_.ap.hypsometric_opt));
  c.set("khdif", str(double(m_.lp.khdif))); c.set("kvdif", str(double(m_.lp.kvdif)));
  c.set("w_damping", std::to_string(m_.lp.w_damping));
  c.set("use_theta_m", m_.st->P_QV >= 0 ? "1" : "0");
  c.set("mp_physics", std::to_string(m_.rk && m_.rk->physics() ? m_.rk->physics()->params().mp_physics : 0));
  // check_a_mundo: num_soil_layers = 5 without a land surface model, whatever the namelist says (the history
  // soil fields are then 5 layers of zeros); Noah / Noah-MP set it below
  if (!m_.suite) wc_->set("num_soil_layers", "5");
  if (m_.suite) {
    const auto& o = m_.suite->options();
    c.set("sf_sfclay_physics", std::to_string(o.sf_sfclay_physics));
    c.set("sf_surface_physics", std::to_string(o.sf_surface_physics));
    c.set("bl_pbl_physics", std::to_string(o.bl_pbl_physics));
    c.set("cu_physics", std::to_string(o.cu_physics));
    c.set("ra_lw_physics", std::to_string(o.rad.ra_lw_physics));
    c.set("ra_sw_physics", std::to_string(o.rad.ra_sw_physics));
    c.set("num_soil_layers", std::to_string(o.num_soil_layers));
    c.set("ishallow", std::to_string(o.ishallow)); c.set("isfflx", std::to_string(o.isfflx));
    c.set("convtrans_avglen_m", str(o.convtrans_avglen_m));
    c.set("isftcflx", std::to_string(o.isftcflx)); c.set("iz0tlnd", std::to_string(o.iz0tlnd));
    c.set("ysu_topdown_pblmix", std::to_string(o.ysu_topdown_pblmix));
    c.set("kim_tofd", o.kim_tofd ? ".true." : ".false.");
    if (o.sf_surface_physics == 4) {
      c.set("dveg", std::to_string(o.noahmp.dveg));
      c.set("opt_crs", std::to_string(o.noahmp.opt_crs));
      c.set("opt_btr", std::to_string(o.noahmp.opt_btr));
      c.set("opt_run", std::to_string(o.noahmp.opt_runsub));
      c.set("opt_sfc", std::to_string(o.noahmp.opt_sfc));
      c.set("opt_frz", std::to_string(o.noahmp.opt_frz));
      c.set("opt_inf", std::to_string(o.noahmp.opt_inf));
      c.set("opt_rad", std::to_string(o.noahmp.opt_rad));
      c.set("opt_alb", std::to_string(o.noahmp.opt_alb));
      c.set("opt_snf", std::to_string(o.noahmp.opt_snf));
      c.set("opt_tbot", std::to_string(o.noahmp.opt_tbot));
      c.set("opt_stc", std::to_string(o.noahmp.opt_stc));
      c.set("opt_gla", std::to_string(o.noahmp.opt_gla));
      c.set("opt_rsf", std::to_string(o.noahmp.opt_rsf));
      c.set("opt_soil", std::to_string(o.noahmp.opt_soil));
      c.set("opt_pedo", std::to_string(o.noahmp.opt_pedo));
      c.set("opt_crop", std::to_string(o.noahmp.opt_crop));
      c.set("opt_irr", std::to_string(o.noahmp.opt_irr));
      c.set("opt_irrm", std::to_string(o.noahmp.opt_irrm));
      c.set("opt_infdv", std::to_string(o.noahmp.opt_infdv));
      c.set("opt_tdrn", std::to_string(o.noahmp.opt_tdrn));
    }
  } else {
    for (const char* p : {"sf_sfclay_physics", "sf_surface_physics", "bl_pbl_physics", "cu_physics", "ra_lw_physics",
                          "ra_sw_physics"})
      c.set(p, "0");
  }
  wrf_io::derive_config(*wc_, std::max(1, wc_->integer("max_dom")), false);
  frames_per_outfile_ = std::max(1, wc_->integer("frames_per_outfile"));
  if (cfg.get_int("phys.sst_update", 0) == 1)
    auxinput4_interval_s_ = cfg.get_double("phys.auxinput4_interval", 360.0) * 60.0;
  outname_ = "wrfout_" + domain_tag_ + "_";
  if (!real_dir.empty()) {
    wrfinput_ = real_dir + "/" + cfg.get_string("real.wrfinput", "wrfinput_" + domain_tag_);
    if (domain_tag_ != "d01") wrfinput_ = real_dir + "/wrfinput_" + domain_tag_;
    if (!std::filesystem::exists(wrfinput_)) wrfinput_.clear();
  }
  if (!wrfinput_.empty() && nc_open(wrfinput_.c_str(), NC_NOWRITE, &wrfinput_id_) != NC_NOERR) wrfinput_id_ = -1;
  if (wrfinput_id_ >= 0) {  // MAPFAC_M/U/V are the x components of the projection (the model keeps those)
    auto same = [&](const char* a, const char* bname, std::size_t n) {
      std::vector<double> x, y;
      return from_wrfinput(a, n, x) && from_wrfinput(bname, n, y) && x == y;
    };
    same_m_ = same("MAPFAC_M", "MAPFAC_MX", std::size_t(g.nx) * g.ny);
    same_u_ = same("MAPFAC_U", "MAPFAC_UX", std::size_t(g.nx + 1) * g.ny);
    same_v_ = same("MAPFAC_V", "MAPFAC_VX", std::size_t(g.nx) * (g.ny + 1));
  }
  if (wrfinput_id_ >= 0) {  // phy_init: landuse_init (all fields without column physics), UST
    const std::size_t n2 = std::size_t(g.nx) * g.ny;
    std::vector<double> lu, xice, snowc;
    if (from_wrfinput("LU_INDEX", n2, lu) && from_wrfinput("SEAICE", n2, xice)) {
      if (!from_wrfinput("SNOWC", n2, snowc)) snowc.assign(n2, 0.0);
      auto att_i = [&](const char* a, int d) { int v = d; nc_get_att_int(wrfinput_id_, NC_GLOBAL, a, &v); return v; };
      double cen_lat = 0;
      nc_get_att_double(wrfinput_id_, NC_GLOBAL, "CEN_LAT", &cen_lat);
      std::string mminlu;
      std::size_t len = 0;
      if (nc_inq_attlen(wrfinput_id_, NC_GLOBAL, "MMINLU", &len) == NC_NOERR) {
        mminlu.assign(len, ' ');
        nc_get_att_text(wrfinput_id_, NC_GLOBAL, "MMINLU", mminlu.data());
        while (!mminlu.empty() && (mminlu.back() == ' ' || mminlu.back() == '\0')) mminlu.pop_back();
      }
      const std::string tables = cfg.get_string("phys.table_dir", real_dir);
      if (std::filesystem::exists(tables + "/LANDUSE.TBL")) {
        const LanduseInit li = landuse_init(tables, mminlu, att_i("JULDAY", 0), cen_lat, att_i("ISWATER", 17),
                                            att_i("ISICE", 15), lu, xice, snowc,
                                            int(cfg.get_int("phys.fractional_seaice", 0)));
        auto rl = [](const std::vector<double>& v) {  // WRF REAL
          std::vector<double> r(v.size());
          for (std::size_t p = 0; p < v.size(); ++p) r[p] = double(real(v[p]));
          return r;
        };
        if (!m_.suite) {
          phy_init_["ALBBCK"] = rl(li.albbck); phy_init_["ALBEDO"] = rl(li.albedo); phy_init_["ZNT"] = rl(li.z0);
          phy_init_["EMISS"] = rl(li.emiss); phy_init_["Z0"] = rl(li.z0); phy_init_["EMBCK"] = rl(li.embck);
          phy_init_["MAVAIL"] = rl(li.mavail);
        }
        phy_init_["THC"] = rl(li.thc);
        phy_init_["LU_STATE"] = rl(li.lu_state);
        phy_init_["LANDUSE_ISICE"] = {double(att_i("ISICE", 15))};
        phy_init_["LANDUSE_LUCATS"] = {double(li.lucats)};
        phy_init_["LANDUSE_LUSEAS"] = {double(li.luseas)};
        phy_init_["LANDUSE_ISN"] = {double(li.isn)};
      }
      if (!m_.suite) phy_init_["UST"] = std::vector<double>(n2, double(real(0.0001)));
    }
  }
  if (wrfinput_id_ < 0 && !m_.suite) {  // ideal: phy_init's landuse_init for MMINLU '    ', land (XLAND = 1)
    const std::size_t n2 = std::size_t(g.nx) * g.ny;
    auto c = [&](double v) { return std::vector<double>(n2, double(real(v))); };
    phy_init_["XLAND"] = c(1.0); phy_init_["ALBBCK"] = c(0.2); phy_init_["ALBEDO"] = c(0.2);
    phy_init_["EMBCK"] = c(0.85); phy_init_["EMISS"] = c(0.85); phy_init_["THC"] = c(0.04);
    phy_init_["ZNT"] = c(0.1); phy_init_["Z0"] = c(0.1); phy_init_["MAVAIL"] = c(0.3); phy_init_["UST"] = c(0.0001);
    for (const char* f : {"TSK", "TMN"}) {
      auto it = m_.init_out.find(f);
      if (it != m_.init_out.end()) phy_init_[f] = it->second;
    }
  }
  if (!m_.suite) {  // start_em: ctopo = ctopo2 = 1 (topo_wind = 0)
    phy_init_["CTOPO"] = std::vector<double>(std::size_t(g.nx) * g.ny, 1.0);
    phy_init_["CTOPO2"] = phy_init_["CTOPO"];
  }
  start_em_slopes();
  // integer state of the restart file: physics call intervals (radt = bldt = cudt = 0 -> every step),
  // adaptive time step (off), the serial patch size
  auto stepn = [&](const char* opt) {
    const double m = wc_->num(opt);
    return double(std::max(1L, std::lround(m * 60. / g.dt)));
  };
  phy_init_["STEPRA"] = {stepn("radt")};
  phy_init_["STEPBL"] = {m_.suite ? double(m_.suite->stepbl) : stepn("bldt")};
  phy_init_["STEPCU"] = {m_.suite ? double(m_.suite->stepcu) : stepn("cudt")};
  phy_init_["STEPFG"] = {1.0};
  phy_init_["LAST_STEP_UPDATED"] = {-1.0};
  phy_init_["MIN_PTCHSZ"] = {double(std::min(g.nx, g.ny) + 1)};
  {  // lbc_fcx_gcx (spec_exp = 0)
    const int w = std::max(1, wc_->integer("spec_bdy_width")), sz = wc_->integer("spec_zone"),
              rz = wc_->integer("relax_zone");
    std::vector<double> fcx(std::size_t(w), 0.0), gcx(std::size_t(w), 0.0);
    for (int loop = sz + 1; loop <= sz + rz - 1 && loop <= w; ++loop) {  // module_bc_em's order
      fcx[std::size_t(loop - 1)] = double(real(0.1) / real(g.dt) * real(sz + rz - loop) / real(rz - 1));
      gcx[std::size_t(loop - 1)] = double(real(1.0) / real(g.dt) / real(50.) * real(sz + rz - loop) / real(rz - 1));
    }
    phy_init_["FCX"] = fcx;
    phy_init_["GCX"] = gcx;
  }
}

void WrfHistory::start_em_slopes() {  // start_em: slope-radiation arrays and the topography laplacian
  const Grid& g = m_.g;
  const dyn::BaseState& b = *m_.base;
  const int nx = g.nx, ny = g.ny;
  std::vector<double> ht, mx, my, sa, ca;
  if (!copy_base(b.ht, {std::size_t(ny), std::size_t(nx)}, nx, ny, 1, ht)) return;
  const bool mapped = b.mapped;
  if (mapped) {
    copy_base(b.msftx, {std::size_t(ny), std::size_t(nx)}, nx, ny, 1, mx);
    copy_base(b.msfty, {std::size_t(ny), std::size_t(nx)}, nx, ny, 1, my);
    copy_base(b.sina, {std::size_t(ny), std::size_t(nx)}, nx, ny, 1, sa);
    copy_base(b.cosa, {std::size_t(ny), std::size_t(nx)}, nx, ny, 1, ca);
  }
  const std::size_t n2 = std::size_t(nx) * ny;
  std::vector<double> sx(n2), sy(n2), lap(n2), slope(n2), azi(n2);
  const bool px = wc_->logical("periodic_x"), py = wc_->logical("periodic_y");
  const real rdx = real(1.) / real(g.dx), rdy = real(1.) / real(g.dy);
  auto H = [&](int i, int j) { return real(ht[std::size_t(j) * nx + i]); };
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const int im1 = px ? (i - 1 + nx) % nx : std::max(i - 1, 0), ip1 = px ? (i + 1) % nx : std::min(i + 1, nx - 1);
      const int jm1 = py ? (j - 1 + ny) % ny : std::max(j - 1, 0), jp1 = py ? (j + 1) % ny : std::min(j + 1, ny - 1);
      const int di = px ? 2 : ip1 - im1, dj = py ? 2 : jp1 - jm1;
      const std::size_t q = std::size_t(j) * nx + i;
      const real msx = mapped ? real(mx[q]) : real(1), msy = mapped ? real(my[q]) : real(1);
      const real hx = (H(ip1, j) - H(im1, j)) * msx * rdx / real(di);
      const real hy = (H(i, jp1) - H(i, jm1)) * msy * rdy / real(dj);
      sx[q] = double(hx);
      sy[q] = double(hy);
      lap[q] = double((H(ip1, j) + H(im1, j) + H(i, jp1) + H(i, jm1) - H(i, j) * real(4.)) / real(4.));
      const real pi = real(4.) * std::atan(real(1.));
      real sl = std::atan(std::pow(hx * hx + hy * hy, real(.5)));
      real az = 0;
      if (sl < real(1.e-4)) {
        sl = 0;
      } else {
        az = std::atan2(hx, hy) + pi;
        const real sina = mapped ? real(sa[q]) : real(0), cosa = mapped ? real(ca[q]) : real(1);
        if (cosa >= 0) az = az - std::asin(sina);
        else az = az - (pi - std::asin(sina));
      }
      slope[q] = double(sl);
      azi[q] = double(az);
    }
  phy_init_["TOPOSLPX"] = sx;
  phy_init_["TOPOSLPY"] = sy;
  phy_init_["LAP_HGT"] = lap;
  phy_init_["SLOPE"] = slope;
  phy_init_["SLP_AZI"] = azi;
}

void WrfHistory::after_solve(long n) {
  std::vector<double> t;
  const long saved = step_;
  step_ = -1;  // no snapshot lookup: the current state
  const bool ok = value("T", "th_phy_m_t0", {std::size_t(m_.g.nz), std::size_t(m_.g.ny), std::size_t(m_.g.nx)}, t);
  step_ = saved;
  if (ok) t_snap_[n] = t;
}

void WrfHistory::before_step(long n) {
  if (m_.suite) return;
  prep_step_ = n;
}

void WrfHistory::first_rk_prep() {
  if (m_.suite || prep_step_ < 0) return;
  const long n = prep_step_;
  prep_step_ = -1;
  if (!prep_) prep_ = std::make_unique<phys::driver::PhyPrepFields>(m_.g);
  phys::driver::PhyPrepParams pp;
  pp.use_theta_m = m_.st->P_QV >= 0 ? 1 : 0;
  pp.mut_from_mu = true;
  phys::driver::run_phy_prep(*m_.st, *m_.base, *prep_, pp);
  std::vector<double> ph;
  if (copy_field(prep_->p_hyd, {std::size_t(m_.g.nz), std::size_t(m_.g.ny), std::size_t(m_.g.nx)}, ph))
    p_hyd_snap_[n] = std::move(ph);
}

WrfHistory::~WrfHistory() {
  file_.reset();
  if (wrfinput_id_ >= 0) nc_close(wrfinput_id_);
}

bool WrfHistory::from_wrfinput(const std::string& name, std::size_t n, std::vector<double>& out) {
  if (wrfinput_id_ < 0) return false;
  int id;
  if (nc_inq_varid(wrfinput_id_, name.c_str(), &id) != NC_NOERR) return false;
  int nd;
  nc_inq_varndims(wrfinput_id_, id, &nd);
  std::vector<int> dimids(static_cast<std::size_t>(nd));
  nc_inq_vardimid(wrfinput_id_, id, dimids.data());
  std::size_t total = 1;
  std::vector<std::size_t> st(static_cast<std::size_t>(nd), 0), ct(static_cast<std::size_t>(nd), 1);
  for (int d = 1; d < nd; ++d) {
    nc_inq_dimlen(wrfinput_id_, dimids[std::size_t(d)], &ct[std::size_t(d)]);
    total *= ct[std::size_t(d)];
  }
  if (total != n) return false;
  nc_type t;
  nc_inq_vartype(wrfinput_id_, id, &t);
  if (t == NC_CHAR) return false;
  out.resize(n);
  return nc_get_vara_double(wrfinput_id_, id, st.data(), ct.data(), out.data()) == NC_NOERR;
}

bool WrfHistory::out_value(const wrf_io::OutVar& v, std::vector<double>& out) {
  const std::string& vn = v.var ? v.var->varname : v.name;
  const std::size_t nd = v.dimnames.size();
  if (jrep_ == 1 || nd < 2 || v.dimnames[nd - 2].rfind("south_north", 0) != 0) return value(v.name, vn, v.dimlens, out);
  std::vector<std::size_t> d = v.dimlens;
  const std::size_t nj = d[nd - 2], ni = d[nd - 1];
  d[nd - 2] = std::size_t(m_.g.ny);  // also for south_north_stag: an x-z model has one v row
  std::vector<double> in;
  if (!value(v.name, vn, d, in)) return false;
  const std::size_t sj = d[nd - 2], outer = in.size() / (sj * ni);
  out.resize(outer * nj * ni);
  for (std::size_t o = 0; o < outer; ++o)
    for (std::size_t j = 0; j < nj; ++j)
      std::copy_n(in.begin() + std::ptrdiff_t((o * sj + j % sj) * ni), ni, out.begin() + std::ptrdiff_t((o * nj + j) * ni));
  return true;
}

bool WrfHistory::value(const std::string& name, const std::string& varname, const std::vector<std::size_t>& dims,
                       std::vector<double>& out) {
  const dyn::State& s = *m_.st;
  const dyn::BaseState& b = *m_.base;
  std::size_t n = 1;
  for (std::size_t l : dims) n *= l;
  auto mass3 = [&](auto f) {  // (k, j, i) over dims
    out.resize(n);
    const std::size_t nk = dims[0], nj = dims[1], ni = dims[2];
    for (std::size_t k = 0; k < nk; ++k)
      for (std::size_t j = 0; j < nj; ++j)
        for (std::size_t i = 0; i < ni; ++i) out[(k * nj + j) * ni + i] = f(int(i), int(j), int(k));
    return true;
  };
  if (name.size() > 2 && (name.compare(name.size() - 2, 2, "_2") == 0 || name.compare(name.size() - 2, 2, "_1") == 0) &&
      varname.size() > 2 && varname.compare(varname.size() - 2, 2, name.substr(name.size() - 2)) == 0) {
    const std::string base = name.substr(0, name.size() - 2);
    if (name.back() == '2') return value(base, varname, dims, out);  // the current time level
    if (base == "U") return copy_field(s.u_1, dims, out);
    if (base == "V") return copy_field(s.v_1, dims, out);
    if (base == "W") return copy_field(s.w_1, dims, out);
    if (base == "PH") return copy_field(s.ph_1, dims, out);
    if (base == "THM") return copy_field(s.t_1, dims, out);
    if (base == "MU") return copy_field(s.mu_1, dims, out);
    for (int q = 0; q < s.n_scalars(); ++q)
      if (s.q_names[std::size_t(q)] == "tke" && base == "TKE") return copy_field(s.q_1[std::size_t(q)], dims, out);
  }
  // landuse_init and the surface driver set XICEM = XICE (the previous step's sea ice)
  if (name == "XICEM" && !(m_.suite && m_.suite->options().sst_update == 1) && value("SEAICE", "SEAICE", dims, out))
    return true;
  if (name == "U") return copy_field(s.u_2, dims, out);
  if (name == "V") return copy_field(s.v_2, dims, out);
  if (name == "W") return copy_field(s.w_2, dims, out);
  if (name == "PH") return copy_field(s.ph_2, dims, out);
  if (name == "THM") return copy_field(s.t_2, dims, out);
  if (name == "AL") return copy_field(s.al, dims, out);
  if (name == "ALT") return copy_field(s.alt, dims, out);
  if (name == "PHP") return copy_field(s.php, dims, out);
  if (name == "WW") return copy_field(s.ww, dims, out);
  if (name == "WW_M") return copy_field(s.ww_m, dims, out);
  if (name == "H_DIABATIC") return copy_field(s.h_diabatic, dims, out);
  if (name == "P") return copy_field(s.p, dims, out);
  if (name == "MU") return copy_field(s.mu_2, dims, out);
  if (name == "T") {  // th_phy_m_t0 as moist_physics_finish_em left it at the end of the step
    auto it = t_snap_.find(step_);
    if (it != t_snap_.end() && it->second.size() == n) { out = it->second; return true; }
    if (step_ == 0 && from_wrfinput("T", n, out)) return true;  // initial: the input's th_phy_m_t0
    if (step_ == 0 && wrfinput_id_ < 0) {  // ideal: the dry theta of the initialisation
      auto t0 = m_.init_out.find("T");
      if (t0 != m_.init_out.end() && t0->second.size() == n) { out = t0->second; return true; }
    }
  }
  if (name == "T") {  // th_phy_m_t0: dry potential temperature - t0
    std::vector<real> ht, hq;
    s.t_2.download(ht);
    if (s.P_QV >= 0) s.q_2[s.P_QV].download(hq);
    const Layout& L = s.t_2.layout();
    return mass3([&](int i, int j, int k) {
      const double tm = double(ht[L.index(i, j, k)]) + b.t0;
      const double th = hq.empty() || m_.st->P_QV < 0 ? tm : tm / (1.0 + cst::rvovrd * double(hq[L.index(i, j, k)]));
      return th - b.t0;
    });
  }
  if (name == "ITIMESTEP") { out.assign(1, double(step_)); return true; }
  if (m_.suite && dims.size() == 2) {
    const phys::driver::IntField2* fi = nullptr;
    if (name == "KPBL") fi = &m_.suite->kpbl;
    else if (name == "IVGTYP") fi = &m_.suite->ivgtyp;
    else if (name == "ISLTYP") fi = &m_.suite->isltyp;
    else if (m_.suite->options().sf_surface_physics == 4) {
      if (name == "ISNOW" || name == "ISNOWXY") fi = &m_.suite->isnowxy;
      else if (name == "CROPCAT") fi = &m_.suite->cropcat;
      else if (name == "PGS" || name == "PGSXY") fi = &m_.suite->pgsxy;
    }
    if (fi) {
      std::vector<int> h;
      fi->download(h);
      const Layout& L = fi->layout();
      out.resize(n);
      for (std::size_t j = 0; j < dims[0]; ++j)
        for (std::size_t i = 0; i < dims[1]; ++i) out[j * dims[1] + i] = double(h[L.index(int(i), int(j), 0)]);
      return true;
    }
  }
  if (wrfinput_id_ < 0) {  // ideal runs: the vertical coordinate and constants of the model (real data: wrfinput)
    using BS = dyn::BaseState;
    static const std::map<std::string, Array1 BS::*> a1 = {
        {"ZNU", &BS::znu}, {"ZNW", &BS::znw}, {"DNW", &BS::dnw}, {"DN", &BS::dn}, {"RDNW", &BS::rdnw}, {"RDN", &BS::rdn},
        {"FNM", &BS::fnm}, {"FNP", &BS::fnp}, {"C1H", &BS::c1h}, {"C2H", &BS::c2h}, {"C3H", &BS::c3h}, {"C4H", &BS::c4h},
        {"C1F", &BS::c1f}, {"C2F", &BS::c2f}, {"C3F", &BS::c3f}, {"C4F", &BS::c4f}};
    if (auto it = a1.find(name); it != a1.end()) {
      const Array1& v = b.*(it->second);
      if (std::size_t(v.size()) < n) return false;
      out.resize(n);
      for (std::size_t k = 0; k < n; ++k) out[k] = v[int(k)];
      return true;
    }
    static const std::map<std::string, double BS::*> a0 = {
        {"P_TOP", &BS::p_top}, {"CF1", &BS::cf1}, {"CF2", &BS::cf2}, {"CF3", &BS::cf3}, {"CFN", &BS::cfn},
        {"CFN1", &BS::cfn1}, {"RDX", &BS::rdx}, {"RDY", &BS::rdy}};
    if (auto it = a0.find(name); it != a0.end() && n == 1) { out.assign(1, b.*(it->second)); return true; }
    if (name == "THIS_IS_AN_IDEAL_RUN" && n == 1) { out.assign(1, 1.0); return true; }
  }
  {  // base state (the model's: blend_terrain and nest feedback change it)
    const Grid& g = m_.g;
    const int nx = g.nx, ny = g.ny, nz = g.nz;
    if (name == "PHB") return copy_base(b.phb, dims, nx, ny, nz + 1, out);
    if (name == "PB") return copy_base(b.pb, dims, nx, ny, nz, out);
    if (name == "ALB") return copy_base(b.alb, dims, nx, ny, nz, out);
    if (name == "T_INIT") return copy_base(b.t_init, dims, nx, ny, nz, out);
    if (name == "MUB") return copy_base(b.mub, dims, nx, ny, 1, out);
    if (name == "HGT") return copy_base(b.ht, dims, nx, ny, 1, out);
    if (b.mapped || wrfinput_id_ < 0) {  // ideal: m = 1, f = e = 0 (module_initialize_ideal ideal_landmap)
      if (wrfinput_id_ < 0) {
        if (name == "COSALPHA") return copy_base(b.cosa, dims, nx, ny, 1, out);
        if (name == "SINALPHA") return copy_base(b.sina, dims, nx, ny, 1, out);
        if (name == "F") return copy_base(b.f, dims, nx, ny, 1, out);
        if (name == "E") return copy_base(b.e, dims, nx, ny, 1, out);
      }
      if (name == "MAPFAC_MX" || (name == "MAPFAC_M" && same_m_)) return copy_base(b.msftx, dims, nx, ny, 1, out);
      if (name == "MAPFAC_MY") return copy_base(b.msfty, dims, nx, ny, 1, out);
      if (name == "MAPFAC_UX" || (name == "MAPFAC_U" && same_u_)) return copy_base(b.msfux, dims, nx + 1, ny, 1, out);
      if (name == "MAPFAC_UY") return copy_base(b.msfuy, dims, nx + 1, ny, 1, out);
      const int nyv = g.nys(Stagger::V);  // ny + 1 (ny for an x-z run)
      if (name == "MAPFAC_VX" || (name == "MAPFAC_V" && same_v_)) return copy_base(b.msfvx, dims, nx, nyv, 1, out);
      if (name == "MAPFAC_VY") return copy_base(b.msfvy, dims, nx, nyv, 1, out);
      if (name == "MF_VX_INV") return copy_base(b.msfvx_inv, dims, nx, nyv, 1, out);
    }
  }
  if (name == "P_HYD" && step_ > 0) {
    auto it = p_hyd_snap_.find(step_);
    if (it != p_hyd_snap_.end() && it->second.size() == n) { out = it->second; return true; }
  }
  if (!model_only_) {  // phy_init / start_em values (nest feedback replaces them on the parent)
    auto it = phy_init_.find(name);
    if (it != phy_init_.end() && it->second.size() == n) {
      auto fb = fed_.find(name);
      out = fb != fed_.end() && fb->second.size() == n ? fb->second : it->second;
      return true;
    }
  }
  if (name == "P_HYD" && !(m_.suite && step_ > 0)) {  // start_em: p_hyd = p + pb until phy_prep sets it
    if (p_hyd0_.empty()) {
      std::vector<real> hp, hpb;
      s.p.download(hp);
      b.pb.download(hpb);
      const Layout &L = s.p.layout(), &Lb = b.pb.layout();
      mass3([&](int i, int j, int k) { return double(hp[L.index(i, j, k)] + hpb[Lb.index(i, j, k)]); });
      p_hyd0_ = out;
    }
    out = p_hyd0_;
    return out.size() == n;
  }
  if (name == "SAVE_TOPO_FROM_REAL") { out.assign(1, 0.0); return true; }
  if (name == "DTBC") {  // time since the last boundary update (spec_bdyupdate's dtbc)
    // nests: med_nest_force sets dtbc = 0 at every parent step
    const double iv = wc_->domain() > 1 ? m_.g.dt * std::max(1, wc_->integer("parent_time_step_ratio"))
                                        : wc_->num("interval_seconds");
    out.assign(1, iv > 0 ? std::fmod(time_sec_, iv) == 0 && time_sec_ > 0 ? iv : std::fmod(time_sec_, iv) : time_sec_);
    return true;
  }  // output_wrf: 0 in model output
  if (name == "AREA2D" || name == "DX2D") {  // compute_2d_dx_area: dx*dy and dx (REAL)
    const real dx = real(m_.g.dx), dy = real(m_.g.dy);
    out.assign(n, name == "AREA2D" ? double(dx * dy) : double(dx));
    return true;
  }
  if (name == "XTIME") { out.assign(1, time_sec_ / 60.0); return true; }
  for (int q = 0; q < s.n_scalars(); ++q) {  // moist species and scalars under WRF's names
    std::string qn = s.q_names[std::size_t(q)];
    static const std::map<std::string, std::string> wn = {{"qv", "QVAPOR"}, {"qc", "QCLOUD"}, {"qr", "QRAIN"},
                                                          {"qi", "QICE"},   {"qs", "QSNOW"},  {"qg", "QGRAUP"},
                                                          {"qni", "QNICE"}, {"qnr", "QNRAIN"}, {"tke", "TKE"}};
    auto it = wn.find(qn);
    if (it != wn.end() && it->second == name) return copy_field(s.q_2[std::size_t(q)], dims, out);
  }
  const std::string vn = to_lower(varname);
  if (vn == "rho" && dims.size() == 3) {
    // grid%rho: moist_physics_prep_em sets it on the tile without the spec zone; the spec-zone ring keeps
    // phy_prep's value of the step start
    const Field3* ring = m_.suite ? &m_.suite->prep.rho : prep_ ? &prep_->rho : nullptr;
    const bool spec = m_.g.bc_x == BC::Specified || m_.g.bc_y == BC::Specified;
    const Field3* inner = nullptr;
    for (const Field3* f : f3_)
      if (f->name() == vn) { inner = f; break; }
    std::vector<double> r;
    if (inner && ring && spec && copy_field(*inner, dims, out) && copy_field(*ring, dims, r)) {
      const std::size_t nk = dims[0], nj = dims[1], ni = dims[2];
      const bool px = m_.g.bc_x != BC::Periodic;
      for (std::size_t k = 0; k < nk; ++k)
        for (std::size_t j = 0; j < nj; ++j)
          for (std::size_t i = 0; i < ni; ++i)
            if (j == 0 || j + 1 == nj || (px && (i == 0 || i + 1 == ni)))
              out[(k * nj + j) * ni + i] = r[(k * nj + j) * ni + i];
      return true;
    }
  }
  for (const Field2* f : f2_)
    if (f->name() == vn && copy_field(*f, dims, out)) return true;
  for (const Field3* f : f3_)
    if (f->name() == vn && copy_field(*f, dims, out)) return true;
  if (prep_)  // phy_prep of a run without column physics (first_rk_step_part1 runs it every step)
    for (const Field3* f : {&prep_->t_phy, &prep_->u_phy, &prep_->v_phy, &prep_->rho, &prep_->p_hyd_w, &prep_->z,
                            &prep_->z_at_w, &prep_->p8w, &prep_->t8w})
      if (f->name() == vn && copy_field(*f, dims, out)) return true;
  if ((name == "DECLIN" || name == "SOLCON") && m_.suite && m_.suite->radiation()) {  // radiation_driver's radconst
    out.assign(1, name == "DECLIN" ? m_.suite->radiation()->declin() : m_.suite->radiation()->solcon());
    return true;
  }
  if (name == "HAS_REQC" || name == "HAS_REQI" || name == "HAS_REQS") {  // microphysics radii to the radiation
    out.assign(1, m_.rk && m_.rk->physics() && m_.rk->physics()->params().microp_re ? 1.0 : 0.0);
    return true;
  }
  if (name == "AKPBL" && m_.suite && m_.suite->options().cu_physics != 0 && dims.size() == 2) {  // cumulus_driver
    std::vector<int> h;
    m_.suite->kpbl.download(h);
    const Layout& L = m_.suite->kpbl.layout();
    out.resize(n);
    for (std::size_t j = 0; j < dims[0]; ++j)
      for (std::size_t i = 0; i < dims[1]; ++i) out[j * dims[1] + i] = double(h[L.index(int(i), int(j), 0)]);
    return true;
  }
  if (model_only_) return false;
  auto fb = fed_.find(name);  // nest feedback of a field the model does not carry
  if (fb != fed_.end() && fb->second.size() == n) { out = fb->second; return true; }
  return from_wrfinput(name, n, out);
}

namespace {

// share/interp_fcn.F feedback functions on NetCDF-ordered (k, j, i) arrays: coarse c (nk, cny, cnx),
// fine f (nk, fny, fnx); ipos, jpos = i/j_parent_start; staggered dims include the extra point
void feedback(const std::string& fn, std::vector<double>& c, int cnx, int cny, const std::vector<double>& f,
              int fnx, int fny, int nk, int ipos, int jpos, int nri, int nrj, bool xstag, bool ystag,
              int spec_zone, int fide, int fjde) {
  const int istag = xstag ? 0 : 1, jstag = ystag ? 0 : 1;
  auto C = [&](int ci, int k, int cj) -> double& {
    return c[(std::size_t(k) * cny + std::size_t(cj - 1)) * cnx + std::size_t(ci - 1)];
  };
  auto F = [&](int ni, int k, int nj) {
    return f[(std::size_t(k) * fny + std::size_t(nj - 1)) * fnx + std::size_t(ni - 1)];
  };
  const bool mark = fn == "mark_domain";
  const int sz = mark ? 1 : spec_zone;
  const int cj0 = jpos + sz, cj1 = std::min(jpos + (fjde - 1) / nrj - jstag - sz, cny);
  const int ci0 = ipos + sz, ci1 = std::min(ipos + (fide - 1) / nri - istag - sz, cnx);
  for (int cj = std::max(cj0, 1); cj <= cj1; ++cj)
    for (int k = 0; k < nk; ++k)
      for (int ci = std::max(ci0, 1); ci <= ci1; ++ci) {
        if (mark) {
          C(ci, k, cj) = 9021000.;
        } else if (fn == "copy_fcnm" || fn == "copy_fcni") {
          if (nrj % 2 != 0) C(ci, k, cj) = F((ci - ipos) * nri + istag + 1, k, (cj - jpos) * nrj + jstag + 1);
          else C(ci, k, cj) = F((ci - ipos) * nri + 1 + nri / 2 - 1, k, (cj - jpos) * nrj + 1 + nrj / 2 - 1);
        } else {  // copy_fcn: average of the fine points of the coarse cell / face
          double sum = 0.0;
          if (!xstag && !ystag) {
            const double inv = double(real(1.0) / real(nri * nrj));
            for (int dj = 0; dj < nrj; ++dj)
              for (int di = 0; di < nri; ++di)
                sum += double(real(inv) * real(F((ci - ipos) * nri + 1 + di, k, (cj - jpos) * nrj + 1 + dj)));
          } else if (xstag) {
            const double inv = double(real(1.0) / real(nrj));
            for (int dj = 0; dj < nrj; ++dj) sum += double(real(inv) * real(F((ci - ipos) * nri + 1, k, (cj - jpos) * nrj + 1 + dj)));
          } else {
            const double inv = double(real(1.0) / real(nri));
            for (int di = 0; di < nri; ++di) sum += double(real(inv) * real(F((ci - ipos) * nri + 1 + di, k, (cj - jpos) * nrj + 1)));
          }
          C(ci, k, cj) = sum;
        }
      }
}

}  // namespace

// Also without column physics: the microphysics accumulations and h_diabatic are fed back as well
void WrfHistory::feedback_physics(WrfHistory& child, int ipos, int jpos, int ratio) {
  const Grid& g = m_.g;
  const Grid& cg = child.m_.g;
  std::vector<wrf_io::OutVar> vars = fb_schema_->history_vars(*wc_, g.nx + 1, g.ny + 1, g.nz + 1, int(sizeof(real)));
  {
    std::set<std::string> have;
    for (const wrf_io::OutVar& v : vars) have.insert(v.name);
    for (wrf_io::OutVar& v : fb_schema_->restart_vars(*wc_, g.nx + 1, g.ny + 1, g.nz + 1, int(sizeof(real))))
      if (!have.count(v.name)) vars.push_back(std::move(v));
  }
  const int spec_zone = std::max(1, wc_->integer("spec_zone"));
  std::set<std::string> done;
  for (const wrf_io::OutVar& v : vars) {
    if (!v.var || v.var->nest_up.empty() || v.dimlens.size() < 2 || !done.insert(v.var->varname).second) continue;
    const std::string fn = v.var->nest_up.substr(0, v.var->nest_up.find(':'));
    if (fn == "mark_domain") continue;
    const std::string vn = to_lower(v.var->varname);
    // the suite's fields, the microphysics accumulations and h_diabatic (state of the dycore)
    auto find2 = [&](const Model& m, const std::vector<const Field2*>& l) -> Field2* {
      for (const Field2* f : l)
        if (f->name() == vn) return const_cast<Field2*>(f);
      if (m.rk && m.rk->physics()) {
        dyn::PhysicsCoupling& pc = *m.rk->physics();
        for (Field2* f : {&pc.rainnc, &pc.snownc, &pc.graupelnc})
          if (f->name() == vn && f->layout().total > 0) return f;
      }
      return nullptr;
    };
    auto find3 = [&](const Model& m, const std::vector<const Field3*>& l) -> Field3* {
      for (const Field3* f : l)
        if (f->name() == vn) return const_cast<Field3*>(f);
      if (vn == "h_diabatic" && m.st->n_moist > 0) return &m.st->h_diabatic;
      return nullptr;
    };
    Field2* p2 = find2(m_, f2_);
    Field3* p3 = p2 ? nullptr : find3(m_, f3_);
    if (!p2 && !p3) {  // integer fields (copy_fcni / copy_fcnm): land-use and soil categories, Noah-MP's
      auto ifield = [&](const Model& m) -> phys::driver::IntField2* {
        if (!m.suite) return nullptr;
        phys::driver::PhysicsSuite& su = *m.suite;
        if (v.name == "IVGTYP") return &su.ivgtyp;
        if (v.name == "ISLTYP") return &su.isltyp;
        if (su.options().sf_surface_physics == 4) {
          if (v.name == "ISNOWXY") return &su.isnowxy;
          if (v.name == "CROPCAT") return &su.cropcat;
          if (v.name == "PGSXY") return &su.pgsxy;
        }
        return nullptr;
      };
      phys::driver::IntField2* pi = ifield(m_);
      const phys::driver::IntField2* ci = ifield(child.m_);
      if (!pi || !ci || v.dimlens.size() != 2) continue;
      auto get = [](const phys::driver::IntField2& fi, int ni, int nj) {
        std::vector<int> h;
        fi.download(h);
        std::vector<double> o(std::size_t(ni) * nj);
        for (int j = 0; j < nj; ++j)
          for (int i = 0; i < ni; ++i) o[std::size_t(j) * ni + i] = double(h[fi.layout().index(i, j, 0)]);
        return o;
      };
      std::vector<double> c = get(*pi, g.nx, g.ny), f = get(*ci, cg.nx, cg.ny);
      feedback(fn, c, g.nx, g.ny, f, cg.nx, cg.ny, 1, ipos, jpos, ratio, ratio, false, false, spec_zone, cg.nx + 1,
               cg.ny + 1);
      std::vector<int> h;
      pi->download(h);
      for (int j = 0; j < g.ny; ++j)
        for (int i = 0; i < g.nx; ++i) h[pi->layout().index(i, j, 0)] = int(std::lround(c[std::size_t(j) * g.nx + i]));
      pi->upload(h);
      continue;
    }
    std::vector<double> c, f;
    if (!(p2 ? copy_field(*p2, v.dimlens, c) : copy_field(*p3, v.dimlens, c))) continue;
    const bool xs = v.stagger == "X", ys = v.stagger == "Y";
    std::vector<std::size_t> cd = v.dimlens;
    cd[cd.size() - 1] = std::size_t(cg.nx + (xs ? 1 : 0));
    cd[cd.size() - 2] = std::size_t(cg.ny + (ys ? 1 : 0));
    const Field2* c2 = p2 ? find2(child.m_, child.f2_) : nullptr;
    const Field3* c3 = p3 ? find3(child.m_, child.f3_) : nullptr;
    if (!(c2 ? copy_field(*c2, cd, f) : c3 && copy_field(*c3, cd, f))) continue;
    const int nk = v.dimlens.size() == 3 ? int(v.dimlens[0]) : 1;
    feedback(fn, c, int(v.dimlens.back()), int(v.dimlens[v.dimlens.size() - 2]), f, int(cd.back()),
             int(cd[cd.size() - 2]), nk, ipos, jpos, ratio, ratio, xs, ys, spec_zone, cg.nx + 1, cg.ny + 1);
    if (p2) store_field(*p2, v.dimlens, c); else store_field(*p3, v.dimlens, c);
  }
}

void WrfHistory::feedback_from(WrfHistory& child, int ipos, int jpos, int ratio) {
  const Grid& g = m_.g;
  const Grid& cg = child.m_.g;
  // history and restart fields (the restart stream carries more of them)
  std::vector<wrf_io::OutVar> vars = fb_schema_->history_vars(*wc_, g.nx + 1, g.ny + 1, g.nz + 1, int(sizeof(real)));
  {
    std::set<std::string> have;
    for (const wrf_io::OutVar& v : vars) have.insert(v.name);
    for (wrf_io::OutVar& v : fb_schema_->restart_vars(*wc_, g.nx + 1, g.ny + 1, g.nz + 1, int(sizeof(real))))
      if (!have.count(v.name)) vars.push_back(std::move(v));
  }
  const int spec_zone = std::max(1, wc_->integer("spec_zone"));
  for (const wrf_io::OutVar& v : vars) {
    if (v.name == "HT_COARSE") {  // med_feedback_domain: parent%ht_coarse = parent%ht
      std::vector<double> h;
      if (value("HGT", "ht", v.dimlens, h)) fed_[v.name] = std::move(h);
      continue;
    }
    if (!v.var || v.var->nest_up.empty() || v.dimlens.size() < 2) continue;
    const std::string fn = v.var->nest_up.substr(0, v.var->nest_up.find(':'));
    std::vector<double> tmp;
    model_only_ = true;
    const bool model = value(v.name, v.var->varname, v.dimlens, tmp);
    model_only_ = false;
    if (model) continue;  // the model feeds back its own fields
    std::vector<double> c;
    if (!value(v.name, v.var->varname, v.dimlens, c)) c.assign([&] { std::size_t n = 1; for (auto l : v.dimlens) n *= l; return n; }(), 0.0);
    const bool xs = v.stagger == "X", ys = v.stagger == "Y";
    std::vector<std::size_t> cd = v.dimlens;
    cd[cd.size() - 1] = std::size_t(cg.nx + (xs ? 1 : 0));
    cd[cd.size() - 2] = std::size_t(cg.ny + (ys ? 1 : 0));
    std::vector<double> f;
    if (!child.value(v.name, v.var->varname, cd, f)) {
      if (fn != "mark_domain") continue;
      f.assign(cd.size() == 3 ? cd[0] * cd[1] * cd[2] : cd[0] * cd[1], 0.0);  // mark_domain reads no fine values
    }
    const int nk = v.dimlens.size() == 3 ? int(v.dimlens[0]) : 1;
    feedback(fn, c, int(v.dimlens.back()), int(v.dimlens[v.dimlens.size() - 2]), f, int(cd.back()),
             int(cd[cd.size() - 2]), nk, ipos, jpos, ratio, ratio, xs, ys, spec_zone, cg.nx + 1, cg.ny + 1);
    fed_[v.name] = std::move(c);
  }
}

wrf_io::DomainAtts WrfHistory::domain_atts() {
  const Grid& g = m_.g;
  wrf_io::DomainAtts d;
  d.title = " OUTPUT FROM CPPWRF V4.8.0 MODEL";
  d.start_date = restart_sec_ > 0 ? wrf_io::add_seconds(start_date_, std::lround(restart_sec_)) : start_date_;
  d.gmt_date = start_date_;
  d.simulation_start_date = wc_->get("simulation_start_date").empty() ? start_date_ : start_date_;
  d.ide = g.nx + 1; d.jde = g.ny * jrep_ + 1; d.kde = g.nz + 1;
  d.grid_id = wc_->domain();
  d.parent_id = wc_->integer("parent_id"); d.i_parent_start = wc_->integer("i_parent_start");
  d.j_parent_start = wc_->integer("j_parent_start"); d.parent_grid_ratio = wc_->integer("parent_grid_ratio");
  if (d.grid_id == 1) { d.parent_id = 0; d.parent_grid_ratio = 1; }  // i/j_parent_start: the namelist value (WRF writes it)
  d.dt = g.dt;
  d.real_bytes = int(sizeof(real));
  d.ideal = wrfinput_id_ < 0;
  if (wrfinput_id_ >= 0) {  // grid metadata of the input file
    auto gatt_d = [&](const char* a, double def) {
      double v = def;
      nc_get_att_double(wrfinput_id_, NC_GLOBAL, a, &v);
      return v;
    };
    auto gatt_i = [&](const char* a, int def) {
      int v = def;
      nc_get_att_int(wrfinput_id_, NC_GLOBAL, a, &v);
      return v;
    };
    for (const char* a : {"CEN_LAT", "CEN_LON", "TRUELAT1", "TRUELAT2", "MOAD_CEN_LAT", "STAND_LON", "POLE_LAT", "POLE_LON"})
      d.proj[a] = gatt_d(a, 0);
    d.map_proj = gatt_i("MAP_PROJ", 0);
    std::size_t len = 0;
    if (nc_inq_attlen(wrfinput_id_, NC_GLOBAL, "MMINLU", &len) == NC_NOERR) {
      std::string s(len, ' ');
      nc_get_att_text(wrfinput_id_, NC_GLOBAL, "MMINLU", s.data());
      d.mminlu = s;
    }
    d.num_land_cat = gatt_i("NUM_LAND_CAT", 21); d.iswater = gatt_i("ISWATER", 17); d.islake = gatt_i("ISLAKE", -1);
    d.isice = gatt_i("ISICE", 15); d.isurban = gatt_i("ISURBAN", 13); d.isoilwater = gatt_i("ISOILWATER", 14);
    wc_->set("num_land_cat", std::to_string(d.num_land_cat));
  } else {  // ideal: no land-use attributes (the namelist defaults)
    d.mminlu = "";
    d.iswater = d.islake = d.isice = d.isurban = d.isoilwater = 0;
  }
  return d;
}

void WrfHistory::write_restart(long step, double time_sec, double history_interval_s, double restart_interval_s) {
  step_ = step;
  time_sec_ = time_sec;
  const std::string date = wrf_io::add_seconds(start_date_, std::lround(time_sec));
  std::filesystem::create_directories(out_dir_);
  const std::string path = out_dir_ + "/wrfrst_" + domain_tag_ + "_" + date;
  const Grid& g = m_.g;
  wrf_io::DomainAtts d = domain_atts();
  // output_wrf restart_only: the alarms WRF created (all but 26 and 54) at this time; the history alarm has
  // just rung (or rings every history_interval), the restart alarm is ringing, 55 counts from the start
  const long el = std::lround(time_sec);
  for (int a = 1; a <= 55; ++a) {
    if (a == 26 || a == 54 || (a == 52 && d.grid_id > 1)) continue;  // nests have no boundary alarm
    char n1[40], n2[40];
    std::snprintf(n1, sizeof n1, "WRF_ALARM_ISRINGING_%02d", a);
    std::snprintf(n2, sizeof n2, "WRF_ALARM_SECS_TIL_NEXT_RING_%02d", a);
    int ring = 0, secs = 1700000000;
    if (a == 1 && history_interval_s > 0) {
      const long iv = std::lround(history_interval_s);
      secs = int(iv - el % iv);
      // the last step: med_restart_out runs before med_last_solve_io writes that frame
      // nests: med_restart_out of the parent runs before the nest's next med_before_solve_io turns it off
      if ((step == m_.nsteps || d.grid_id > 1) && el % iv == 0) ring = 1;
    } else if (a == 30 && auxinput4_interval_s_ > 0) {  // read (and turned off) before med_restart_out
      const long iv = std::lround(auxinput4_interval_s_);
      secs = int(iv - el % iv);
      if (step == m_.nsteps && el % iv == 0) ring = 1;
    } else if (a == 51 && restart_interval_s > 0) {
      ring = 1;
      secs = int(std::lround(restart_interval_s));
    } else if (a == 55) {
      secs = int(-el);
    }
    d.restart_alarms.push_back({n1, ring});
    d.restart_alarms.push_back({n2, secs});
  }
  std::vector<wrf_io::OutVar> vars = schema_->restart_vars(*wc_, g.nx + 1, g.ny * jrep_ + 1, g.nz + 1, int(sizeof(real)));
  for (wrf_io::OutVar& v : vars)
    if (v.name == "XTIME") {
      std::string sd = start_date_.substr(0, 19);
      sd[10] = ' ';
      v.desc = v.units = "minutes since " + sd;
    }
  wrf_io::StreamFile f(path, std::move(vars), *wc_, d, wrf_io::Stream::restart, !wc_->logical("use_netcdf_classic"));
  f.write(date, [&](const wrf_io::OutVar& v, std::vector<double>& out) {
    return out_value(v, out);
  });
  last_restart_ = path;
}

void WrfHistory::write(long step, double time_sec) {
  step_ = step;
  time_sec_ = time_sec;
  const std::string date = wrf_io::add_seconds(start_date_, std::lround(time_sec));
  if (!file_ || file_->records() >= std::size_t(frames_per_outfile_)) {
    file_.reset();
    std::filesystem::create_directories(out_dir_);
    last_file_ = out_dir_ + "/" + outname_ + date;
    const Grid& g = m_.g;
    const wrf_io::DomainAtts d = domain_atts();
    std::vector<wrf_io::OutVar> vars = schema_->history_vars(*wc_, g.nx + 1, g.ny * jrep_ + 1, g.nz + 1, int(sizeof(real)));
    for (wrf_io::OutVar& v : vars)
      if (v.name == "XTIME") {  // output_wrf: minutes since the simulation start
        std::string sd = start_date_.substr(0, 19);
        sd[10] = ' ';
        v.desc = v.units = "minutes since " + sd;
      }
    file_ = std::make_unique<wrf_io::StreamFile>(last_file_, std::move(vars), *wc_, d, wrf_io::Stream::history,
                                                 !wc_->logical("use_netcdf_classic"));
  }
  file_->write(date, [&](const wrf_io::OutVar& v, std::vector<double>& out) {
    return out_value(v, out);
  });
  t_snap_.erase(t_snap_.begin(), t_snap_.lower_bound(step));
  p_hyd_snap_.erase(p_hyd_snap_.begin(), p_hyd_snap_.lower_bound(step));
}

}  // namespace cppwrf::cases
#endif
