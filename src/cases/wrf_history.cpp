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
  if (int(ni) != L.nx || int(nj) != L.ny || int(nk) != std::max(1, L.nz)) return false;
  std::vector<real> h;
  f.download(h);
  out.resize(nk * nj * ni);
  for (std::size_t k = 0; k < nk; ++k)
    for (std::size_t j = 0; j < nj; ++j)
      for (std::size_t i = 0; i < ni; ++i)
        out[(k * nj + j) * ni + i] = double(h[L.index(int(i), int(j), int(k))]);
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
      if (io.front() != '/' && !real_dir.empty()) io = real_dir + "/" + io;
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
  if (name == "P") return copy_field(s.p, dims, out);
  if (name == "MU") return copy_field(s.mu_2, dims, out);
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
  return from_wrfinput(name, n, out);
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
}

}  // namespace cppwrf::cases
#endif
