#ifdef CPPWRF_HAVE_NETCDF
#include "cases/wrf_history.hpp"

#include <netcdf.h>

#include <cmath>
#include <cstdio>
#include <filesystem>
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

std::string lower(std::string s) {
  for (char& c : s) c = char(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

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
  if (int(ni) != L.nx || int(nj) != L.ny || (int(nk) != std::max(1, L.nz) && !pad)) return false;
  std::vector<real> h;
  f.download(h);
  out.assign(nk * nj * ni, 0.0);
  for (std::size_t k = 0; k < (pad ? nk - 1 : nk); ++k)
    for (std::size_t j = 0; j < nj; ++j)
      for (std::size_t i = 0; i < ni; ++i)
        out[(k * nj + j) * ni + i] = double(h[L.index(int(i), int(j), int(k))]);
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
  c.set("e_we", std::to_string(g.nx + 1)); c.set("e_sn", std::to_string(g.ny + 1));
  c.set("e_vert", std::to_string(g.nz + 1)); c.set("dx", str(g.dx)); c.set("dy", str(g.dy));
  c.set("diff_opt", std::to_string(m_.diff_opt)); c.set("km_opt", std::to_string(m_.km_opt));
  c.set("hybrid_opt", std::to_string(m_.hybrid_opt));
  c.set("damp_opt", std::to_string(m_.ap.damp_opt)); c.set("dampcoef", str(double(m_.ap.dampcoef)));
  c.set("zdamp", str(double(m_.ap.zdamp))); c.set("hypsometric_opt", std::to_string(m_.ap.hypsometric_opt));
  c.set("khdif", str(double(m_.lp.khdif))); c.set("kvdif", str(double(m_.lp.kvdif)));
  c.set("w_damping", std::to_string(m_.lp.w_damping));
  c.set("use_theta_m", m_.st->P_QV >= 0 ? "1" : "0");
  c.set("mp_physics", std::to_string(m_.rk && m_.rk->physics() ? m_.rk->physics()->params().mp_physics : 0));
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
    c.set("isftcflx", std::to_string(o.isftcflx)); c.set("iz0tlnd", std::to_string(o.iz0tlnd));
    c.set("ysu_topdown_pblmix", std::to_string(o.ysu_topdown_pblmix));
    c.set("kim_tofd", o.kim_tofd ? ".true." : ".false.");
  } else {
    for (const char* p : {"sf_sfclay_physics", "sf_surface_physics", "bl_pbl_physics", "cu_physics", "ra_lw_physics",
                          "ra_sw_physics"})
      c.set(p, "0");
  }
  wrf_io::derive_config(*wc_, std::max(1, wc_->integer("max_dom")), false);
  frames_per_outfile_ = std::max(1, wc_->integer("frames_per_outfile"));
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
  if (!m_.suite && wrfinput_id_ >= 0) {  // phy_init of a run without column physics: landuse_init, UST
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
                                            att_i("ISICE", 15), lu, xice, snowc);
        auto rl = [](const std::vector<double>& v) {  // WRF REAL
          std::vector<double> r(v.size());
          for (std::size_t p = 0; p < v.size(); ++p) r[p] = double(real(v[p]));
          return r;
        };
        phy_init_["ALBBCK"] = rl(li.albbck); phy_init_["ALBEDO"] = rl(li.albedo); phy_init_["ZNT"] = rl(li.z0);
        phy_init_["EMISS"] = rl(li.emiss);
      }
      phy_init_["UST"] = std::vector<double>(n2, double(real(0.0001)));
    }
  }
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
  if (name == "U") return copy_field(s.u_2, dims, out);
  if (name == "V") return copy_field(s.v_2, dims, out);
  if (name == "W") return copy_field(s.w_2, dims, out);
  if (name == "PH") return copy_field(s.ph_2, dims, out);
  if (name == "THM") return copy_field(s.t_2, dims, out);
  if (name == "H_DIABATIC") return copy_field(s.h_diabatic, dims, out);
  if (name == "P") return copy_field(s.p, dims, out);
  if (name == "MU") return copy_field(s.mu_2, dims, out);
  if (name == "T") {  // th_phy_m_t0 as moist_physics_finish_em left it at the end of the step
    auto it = t_snap_.find(step_);
    if (it != t_snap_.end() && it->second.size() == n) { out = it->second; return true; }
    if (step_ == 0 && from_wrfinput("T", n, out)) return true;  // initial: the input's th_phy_m_t0
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
  if (name == "KPBL" && m_.suite && dims.size() == 2) {  // integer PBL top level
    std::vector<int> h;
    m_.suite->kpbl.download(h);
    const Layout& L = m_.suite->kpbl.layout();
    out.resize(n);
    for (std::size_t j = 0; j < dims[0]; ++j)
      for (std::size_t i = 0; i < dims[1]; ++i) out[j * dims[1] + i] = double(h[L.index(int(i), int(j), 0)]);
    return true;
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
    if (b.mapped) {
      if (name == "MAPFAC_MX" || (name == "MAPFAC_M" && same_m_)) return copy_base(b.msftx, dims, nx, ny, 1, out);
      if (name == "MAPFAC_MY") return copy_base(b.msfty, dims, nx, ny, 1, out);
      if (name == "MAPFAC_UX" || (name == "MAPFAC_U" && same_u_)) return copy_base(b.msfux, dims, nx + 1, ny, 1, out);
      if (name == "MAPFAC_UY") return copy_base(b.msfuy, dims, nx + 1, ny, 1, out);
      if (name == "MAPFAC_VX" || (name == "MAPFAC_V" && same_v_)) return copy_base(b.msfvx, dims, nx, ny + 1, 1, out);
      if (name == "MAPFAC_VY") return copy_base(b.msfvy, dims, nx, ny + 1, 1, out);
      if (name == "MF_VX_INV") return copy_base(b.msfvx_inv, dims, nx, ny + 1, 1, out);
    }
  }
  if (name == "P_HYD" && step_ > 0) {
    auto it = p_hyd_snap_.find(step_);
    if (it != p_hyd_snap_.end() && it->second.size() == n) { out = it->second; return true; }
  }
  {
    auto it = phy_init_.find(name);
    if (it != phy_init_.end() && it->second.size() == n) { out = it->second; return true; }
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
  if (name == "SAVE_TOPO_FROM_REAL") { out.assign(1, 0.0); return true; }  // output_wrf: 0 in model output
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
  const std::string vn = lower(varname);
  for (const Field2* f : f2_)
    if (f->name() == vn && copy_field(*f, dims, out)) return true;
  for (const Field3* f : f3_)
    if (f->name() == vn && copy_field(*f, dims, out)) return true;
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

void WrfHistory::feedback_from(WrfHistory& child, int ipos, int jpos, int ratio) {
  const Grid& g = m_.g;
  const Grid& cg = child.m_.g;
  std::vector<wrf_io::OutVar> vars = schema_->history_vars(*wc_, g.nx + 1, g.ny + 1, g.nz + 1, int(sizeof(real)));
  const int spec_zone = std::max(1, wc_->integer("spec_zone"));
  for (const wrf_io::OutVar& v : vars) {
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

void WrfHistory::write(long step, double time_sec) {
  step_ = step;
  time_sec_ = time_sec;
  const std::string date = wrf_io::add_seconds(start_date_, std::lround(time_sec));
  if (!file_ || file_->records() >= std::size_t(frames_per_outfile_)) {
    file_.reset();
    std::filesystem::create_directories(out_dir_);
    last_file_ = out_dir_ + "/" + outname_ + date;
    const Grid& g = m_.g;
    wrf_io::DomainAtts d;
    d.title = " OUTPUT FROM CPPWRF V4.8.0 MODEL";
    d.start_date = start_date_;
    d.simulation_start_date = wc_->get("simulation_start_date").empty() ? start_date_ : start_date_;
    d.ide = g.nx + 1; d.jde = g.ny + 1; d.kde = g.nz + 1;
    d.grid_id = wc_->domain();
    d.parent_id = wc_->integer("parent_id"); d.i_parent_start = wc_->integer("i_parent_start");
    d.j_parent_start = wc_->integer("j_parent_start"); d.parent_grid_ratio = wc_->integer("parent_grid_ratio");
    if (d.grid_id == 1) { d.parent_id = 0; d.i_parent_start = 1; d.j_parent_start = 1; d.parent_grid_ratio = 1; }
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
    } else {
      d.mminlu = "";
    }
    std::vector<wrf_io::OutVar> vars = schema_->history_vars(*wc_, g.nx + 1, g.ny + 1, g.nz + 1, int(sizeof(real)));
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
    return value(v.name, v.var ? v.var->varname : v.name, v.dimlens, out);
  });
  t_snap_.erase(t_snap_.begin(), t_snap_.upper_bound(step));
  p_hyd_snap_.erase(p_hyd_snap_.begin(), p_hyd_snap_.upper_bound(step));
}

}  // namespace cppwrf::cases
#endif
