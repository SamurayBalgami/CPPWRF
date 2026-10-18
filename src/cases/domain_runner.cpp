#include "cases/domain_runner.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

#ifdef CPPWRF_HAVE_NETCDF
#include <netcdf.h>
#endif

#include <algorithm>
#include <cctype>
#include <map>
#include <type_traits>

#include "core/constants.hpp"
#include "core/halo.hpp"
#include "core/precision.hpp"
#include "dyn/diagnostics.hpp"

namespace cppwrf::cases {

namespace {

const char* bc_name(BC b) {
  return b == BC::Periodic ? "periodic" : (b == BC::Symmetric ? "symmetric" : "specified");
}

#ifdef CPPWRF_HAVE_NETCDF
// A WRF restart file (wrfrst_d0N): FLAG_RESTART and WRF's Times variable (CPPWRF's own restart files
// have FLAG_RESTART without Times)
bool is_wrfrst(const std::string& path) {
  int id;
  if (nc_open(path.c_str(), NC_NOWRITE, &id) != NC_NOERR) return false;
  int flag = 0, vid;
  const bool r = nc_get_att_int(id, NC_GLOBAL, "FLAG_RESTART", &flag) == NC_NOERR &&
                 nc_inq_varid(id, "Times", &vid) == NC_NOERR;
  nc_close(id);
  return r;
}

// The WRF variable of a restart field (names of save_restart_file)
std::vector<std::string> wrf_names(const std::string& n) {
  auto up = [](std::string x) {
    std::transform(x.begin(), x.end(), x.begin(), [](unsigned char c) { return char(std::toupper(c)); });
    return x;
  };
  static const std::map<std::string, std::string> moist = {
      {"qv", "QVAPOR"}, {"qc", "QCLOUD"}, {"qr", "QRAIN"}, {"qi", "QICE"}, {"qs", "QSNOW"},
      {"qg", "QGRAUP"}, {"qni", "QNICE"}, {"qnr", "QNRAIN"}};
  if (n.rfind("q1_", 0) == 0 || n.rfind("q2_", 0) == 0) {  // moist arrays have one time level (q_1 = q_2
    const std::string q = n.substr(3);                      // at the step start); tke_1 / tke_2
    auto it = moist.find(q);
    if (it != moist.end()) return {it->second};
    return {up(q) + (n[1] == '1' ? "_1" : "_2"), up(q)};
  }
  static const std::map<std::string, std::string> renamed = {
      {"t_1", "THM_1"}, {"t_2", "THM_2"}, {"sfc_t1_lag", "THM_1"},  // grid%t_1 of first_rk_step_part2
      {"sfc_xice", "SEAICE"}, {"sfc_sfcrunoff", "SFROFF"}, {"sfc_udrunoff", "UDROFF"}};
  auto r = renamed.find(n);
  if (r != renamed.end()) return {r->second};
  if (n.rfind("sfc_", 0) == 0) return {up(n.substr(4))};
  return {up(n)};
}

// Restart fields from a WRF wrfrst file (Time = 0); WRF's kms:kme arrays give their first nz levels
void read_wrfrst(const std::string& path, io_netcdf::RestartMeta& meta, std::vector<io_netcdf::VarView>& vars) {
  int id;
  if (nc_open(path.c_str(), NC_NOWRITE, &id) != NC_NOERR) throw std::runtime_error("cannot open " + path);
  auto dimlen = [&](const char* d) {
    int did;
    std::size_t l = 0;
    if (nc_inq_dimid(id, d, &did) == NC_NOERR) nc_inq_dimlen(id, did, &l);
    return int(l);
  };
  meta.nx = dimlen("west_east"); meta.ny = dimlen("south_north"); meta.nz = dimlen("bottom_top");
  float dt = 0;
  nc_get_att_float(id, NC_GLOBAL, "DT", &dt);
  meta.dt = dt;
  int vid, it = 0;
  if (nc_inq_varid(id, "ITIMESTEP", &vid) == NC_NOERR) nc_get_var_int(id, vid, &it);
  meta.step = it;
  double xt = 0;
  if (nc_inq_varid(id, "XTIME", &vid) == NC_NOERR) nc_get_var_double(id, vid, &xt);
  meta.time_sec = xt * 60.;
  std::vector<std::string> missing;
  std::vector<double> buf;
  for (io_netcdf::VarView& v : vars) {
    int var = -1;
    for (const std::string& wn : wrf_names(v.name))
      if (nc_inq_varid(id, wn.c_str(), &var) == NC_NOERR) break; else var = -1;
    if (var < 0) {  // not in WRF's restart stream: WRF's initial value 0, or what phy_init sets on
      missing.push_back(v.name);  // restarts too (myjsfcinit: LOWLYR = 1), which the model's init has
      if (v.name == "sfc_lowlyr") continue;
      real* d = static_cast<real*>(const_cast<void*>(v.data));
      for (int k = 0; k < std::max(1, v.nz); ++k)
        for (int j = 0; j < v.ny; ++j)
          for (int i = 0; i < v.nx; ++i) d[k * v.sk + j * v.sj + i] = real(0);
      continue;
    }
    int nd;
    nc_inq_varndims(id, var, &nd);
    std::vector<int> dids(static_cast<std::size_t>(nd));
    nc_inq_vardimid(id, var, dids.data());
    std::vector<std::size_t> len(static_cast<std::size_t>(nd));
    for (int d = 0; d < nd; ++d) nc_inq_dimlen(id, dids[std::size_t(d)], &len[std::size_t(d)]);
    // (Time, [k,] j, i)
    const bool three = nd == 4;
    if (nd < 3 || int(len[std::size_t(nd - 1)]) != v.nx || int(len[std::size_t(nd - 2)]) != v.ny ||
        (three && int(len[1]) < v.nz) || (!three && v.nz > 1))
      throw std::runtime_error("wrfrst: " + std::string(path) + ": shape of " + v.name + " differs");
    const int nk = three ? v.nz : 1;
    std::vector<std::size_t> st(std::size_t(nd), 0), ct(len);
    ct[0] = 1;
    if (three) ct[1] = std::size_t(nk);
    buf.resize(std::size_t(nk) * std::size_t(v.ny) * std::size_t(v.nx));
    if (nc_get_vara_double(id, var, st.data(), ct.data(), buf.data()) != NC_NOERR)
      throw std::runtime_error("wrfrst: cannot read " + v.name);
    real* d = static_cast<real*>(const_cast<void*>(v.data));
    for (int k = 0; k < nk; ++k)
      for (int j = 0; j < v.ny; ++j)
        for (int i = 0; i < v.nx; ++i)
          d[k * v.sk + j * v.sj + i] = real(buf[(std::size_t(k) * std::size_t(v.ny) + std::size_t(j)) * std::size_t(v.nx) + std::size_t(i)]);
  }
  nc_close(id);
  if (!missing.empty()) {
    std::string m;
    for (const std::string& x : missing) m += " " + x;
    std::fprintf(stderr, "cppwrf: %s: not in the restart stream, start from 0:%s\n", path.c_str(), m.c_str());
  }
}
#endif

}  // namespace

DomainRunner::DomainRunner(const io::Config& cfg, cases::Model& m, const std::string& out_dir,
                           std::string domain_tag)
    : cfg_(cfg), m_(m), out_dir_(out_dir), domain_tag_(std::move(domain_tag)) {
  cname_ = cfg_.get_string("run.case");
  out_fmt_ = cfg_.get_string("output.format", "native");
  write_native_ = (out_fmt_ == "native" || out_fmt_ == "both");
  write_wrfout_ = (out_fmt_ == "wrfout" || out_fmt_ == "both");
  if (out_fmt_ != "native" && out_fmt_ != "wrfout" && out_fmt_ != "both") {
    throw std::runtime_error("unknown output.format: '" + out_fmt_ + "' (native|wrfout|both)");
  }
}

DomainRunner::~DomainRunner() = default;

#ifdef CPPWRF_HAVE_NETCDF
void DomainRunner::save_restart_file(const std::string& path, long step, double time_sec) {
  const Grid& g = m_.g;
  dyn::State& st = *m_.st;
  dyn::RK3* rk = m_.rk.get();

  io_netcdf::RestartMeta meta;
  meta.step = step;
  meta.time_sec = time_sec;
  meta.nx = g.nx; meta.ny = g.ny; meta.nz = g.nz;
  meta.dx = g.dx; meta.dy = g.dy; meta.dt = g.dt;
  meta.case_name = cname_;
  meta.date_str = io_netcdf::format_wrf_date(start_date_, time_sec);

  const io_netcdf::DType dt = sizeof(real) == 8 ? io_netcdf::DType::F64 : io_netcdf::DType::F32;
  std::vector<std::vector<real>> bufs;
  std::vector<io_netcdf::VarView> vars;

  auto add_field3 = [&](const std::string& name, const Field3& f) {
    bufs.emplace_back();
    f.download(bufs.back());
    const Layout& L = f.layout();
    char stag = f.stagger() == Stagger::U ? 'U' : (f.stagger() == Stagger::V ? 'V' : (f.stagger() == Stagger::W ? 'W' : 'M'));
    vars.emplace_back(name, stag, L.nx, L.ny, L.nz, bufs.back().data() + L.offset, dt, L.sj, L.sk);
  };

  auto add_field2 = [&](const std::string& name, const Field2& f) {
    bufs.emplace_back();
    f.download(bufs.back());
    const Layout& L = f.layout();
    vars.emplace_back(name, 'M', L.nx, L.ny, 1, bufs.back().data() + L.offset, dt, L.sj, 0);
  };

  add_field3("u_1", st.u_1); add_field3("u_2", st.u_2);
  add_field3("v_1", st.v_1); add_field3("v_2", st.v_2);
  add_field3("w_1", st.w_1); add_field3("w_2", st.w_2);
  add_field3("t_1", st.t_1); add_field3("t_2", st.t_2);
  add_field3("ph_1", st.ph_1); add_field3("ph_2", st.ph_2);
  add_field2("mu_1", st.mu_1); add_field2("mu_2", st.mu_2);
  add_field3("p", st.p);
  add_field3("al", st.al);

  for (int n = 0; n < st.n_scalars(); ++n) {
    add_field3("q1_" + st.q_names[n], st.q_1[n]);
    add_field3("q2_" + st.q_names[n], st.q_2[n]);
  }
  if (st.n_moist > 0) {
    add_field3("h_diabatic", st.h_diabatic);
  }
  if (rk && rk->physics()) {
    add_field2("rainnc", rk->physics()->rainnc);
    add_field2("rainncv", rk->physics()->rainncv);
    if (rk->physics()->params().mp_physics == 6 || rk->physics()->params().mp_physics == 8) {
      add_field2("snownc", rk->physics()->snownc);
      add_field2("snowncv", rk->physics()->snowncv);
      add_field2("graupelnc", rk->physics()->graupelnc);
      add_field2("graupelncv", rk->physics()->graupelncv);
      add_field2("sr", rk->physics()->sr);
    }
  }

  if (m_.suite) {  // column physics state (soil, surface, PBL, lagged t_1)
    for (auto& [name, f] : m_.suite->carried2()) add_field2(name, *f);
    for (auto& [name, f] : m_.suite->carried3()) add_field3(name, *f);
    for (auto& [name, f] : m_.suite->carried_int()) {  // integers stored as reals (exact)
      std::vector<int> hi;
      f->download(hi);
      bufs.emplace_back(hi.begin(), hi.end());
      const Layout& L = f->layout();
      vars.emplace_back(name, 'M', L.nx, L.ny, 1, bufs.back().data() + L.offset, dt, L.sj, 0);
    }
  }

  io_netcdf::write_restart(path, meta, vars);
}

void DomainRunner::load_restart_file(const std::string& path) {
  const Grid& g = m_.g;
  dyn::State& st = *m_.st;
  dyn::RK3* rk = m_.rk.get();

  io_netcdf::RestartMeta meta;
  const io_netcdf::DType dt = sizeof(real) == 8 ? io_netcdf::DType::F64 : io_netcdf::DType::F32;
  std::vector<std::vector<real>> bufs;
  std::vector<io_netcdf::VarView> vars;

  struct Target3 { std::string name; Field3* f; size_t buf_idx; };
  struct Target2 { std::string name; Field2* f; size_t buf_idx; };
  std::vector<Target3> targets3;
  std::vector<Target2> targets2;

  auto reg_field3 = [&](const std::string& name, Field3& f) {
    size_t idx = bufs.size();
    bufs.emplace_back();
    f.download(bufs.back());
    const Layout& L = f.layout();
    char stag = f.stagger() == Stagger::U ? 'U' : (f.stagger() == Stagger::V ? 'V' : (f.stagger() == Stagger::W ? 'W' : 'M'));
    vars.emplace_back(name, stag, L.nx, L.ny, L.nz, bufs.back().data() + L.offset, dt, L.sj, L.sk);
    targets3.push_back({name, &f, idx});
  };

  auto reg_field2 = [&](const std::string& name, Field2& f) {
    size_t idx = bufs.size();
    bufs.emplace_back();
    f.download(bufs.back());
    const Layout& L = f.layout();
    vars.emplace_back(name, 'M', L.nx, L.ny, 1, bufs.back().data() + L.offset, dt, L.sj, 0);
    targets2.push_back({name, &f, idx});
  };

  reg_field3("u_1", st.u_1); reg_field3("u_2", st.u_2);
  reg_field3("v_1", st.v_1); reg_field3("v_2", st.v_2);
  reg_field3("w_1", st.w_1); reg_field3("w_2", st.w_2);
  reg_field3("t_1", st.t_1); reg_field3("t_2", st.t_2);
  reg_field3("ph_1", st.ph_1); reg_field3("ph_2", st.ph_2);
  reg_field2("mu_1", st.mu_1); reg_field2("mu_2", st.mu_2);
  reg_field3("p", st.p);
  reg_field3("al", st.al);
  const bool wrf_file = is_wrfrst(path);
  // base-state fields of a WRF restart file (nest feedback changes the parent's)
  struct TargetB { BaseFieldT<3>* f3; BaseFieldT<2>* f2; size_t buf_idx; };
  std::vector<TargetB> targets_b;
  if (wrf_file) {  // WRF's restart stream also carries these diagnostics of the last step
    reg_field3("ww", st.ww);
    reg_field3("ww_m", st.ww_m);
    reg_field3("php", st.php);
    reg_field3("alt", st.alt);
    dyn::BaseState& b = *m_.base;
    auto reg_base = [&](const std::string& name, auto& f, char stag) {
      if (f.uniform()) return;
      const Layout& L = f.layout();
      size_t idx = bufs.size();
      bufs.emplace_back();
      f.download(bufs.back());
      vars.emplace_back(name, stag, L.nx, L.ny, L.nz, bufs.back().data() + L.offset, dt, L.sj, L.sk);
      if constexpr (std::is_same_v<std::decay_t<decltype(f)>, BaseFieldT<3>>) targets_b.push_back({&f, nullptr, idx});
      else targets_b.push_back({nullptr, &f, idx});
    };
    reg_base("PHB", b.phb, 'W'); reg_base("PB", b.pb, 'M'); reg_base("ALB", b.alb, 'M'); reg_base("T_INIT", b.t_init, 'M');
    reg_base("MUB", b.mub, 'M'); reg_base("HGT", b.ht, 'M');
    if (b.mapped) {
      reg_base("MAPFAC_MX", b.msftx, 'M'); reg_base("MAPFAC_MY", b.msfty, 'M'); reg_base("MAPFAC_UX", b.msfux, 'U');
      reg_base("MAPFAC_UY", b.msfuy, 'U'); reg_base("MAPFAC_VX", b.msfvx, 'V'); reg_base("MAPFAC_VY", b.msfvy, 'V');
      reg_base("MF_VX_INV", b.msfvx_inv, 'V'); reg_base("F", b.f, 'M'); reg_base("E", b.e, 'M');
      reg_base("SINALPHA", b.sina, 'M'); reg_base("COSALPHA", b.cosa, 'M');
    }
  }

  for (int n = 0; n < st.n_scalars(); ++n) {
    reg_field3("q1_" + st.q_names[n], st.q_1[n]);
    reg_field3("q2_" + st.q_names[n], st.q_2[n]);
  }
  if (st.n_moist > 0) {
    reg_field3("h_diabatic", st.h_diabatic);
  }
  if (rk && rk->physics()) {
    reg_field2("rainnc", rk->physics()->rainnc);
    reg_field2("rainncv", rk->physics()->rainncv);
    if (rk->physics()->params().mp_physics == 6 || rk->physics()->params().mp_physics == 8) {
      reg_field2("snownc", rk->physics()->snownc);
      reg_field2("snowncv", rk->physics()->snowncv);
      reg_field2("graupelnc", rk->physics()->graupelnc);
      reg_field2("graupelncv", rk->physics()->graupelncv);
      reg_field2("sr", rk->physics()->sr);
    }
  }

  struct TargetI { phys::driver::IntField2* f; size_t buf_idx; };
  std::vector<TargetI> targets_i;
  if (m_.suite) {
    for (auto& [name, f] : m_.suite->carried2()) reg_field2(name, *f);
    for (auto& [name, f] : m_.suite->carried3()) reg_field3(name, *f);
    for (auto& [name, f] : m_.suite->carried_int()) {
      std::vector<int> hi;
      f->download(hi);
      size_t idx = bufs.size();
      bufs.emplace_back(hi.begin(), hi.end());
      const Layout& L = f->layout();
      vars.emplace_back(name, 'M', L.nx, L.ny, 1, bufs.back().data() + L.offset, dt, L.sj, 0);
      targets_i.push_back({f, idx});
    }
  }

  if (wrf_file) {
    read_wrfrst(path, meta, vars);
    meta.dx = g.dx; meta.dy = g.dy;
  } else {
    io_netcdf::read_restart(path, meta, vars);
  }

  if (meta.nx != g.nx || meta.ny != g.ny || meta.nz != g.nz) {
    throw std::runtime_error("restart file grid dimension mismatch");
  }
  if (std::abs(meta.dt - g.dt) > 1e-6 * g.dt) {
    throw std::runtime_error("restart file dt (" + std::to_string(meta.dt) + ") does not match grid dt (" + std::to_string(g.dt) + ")");
  }

  start_step = meta.step;
  if (m_.nsteps < start_step) {
    throw std::runtime_error("run.run_time (" + std::to_string(m_.run_time) + "s / step " + std::to_string(m_.nsteps) +
                             ") is less than restart file step " + std::to_string(start_step));
  }

  for (const auto& t : targets3) {
    copy_to_dev(t.f->data() - t.f->layout().offset, bufs[t.buf_idx].data(), t.f->bytes());
    fill_halo(*t.f, g);
  }
  for (const auto& t : targets2) {
    copy_to_dev(t.f->data() - t.f->layout().offset, bufs[t.buf_idx].data(), t.f->bytes());
    fill_halo(*t.f, g);
  }
  for (const auto& t : targets_b) {
    if (t.f3) copy_to_dev(t.f3->data() - t.f3->layout().offset, bufs[t.buf_idx].data(), t.f3->layout().total * sizeof(real));
    else copy_to_dev(t.f2->data() - t.f2->layout().offset, bufs[t.buf_idx].data(), t.f2->layout().total * sizeof(real));
  }
  for (const auto& t : targets_i) {
    std::vector<int> hi(bufs[t.buf_idx].size());
    for (size_t k = 0; k < hi.size(); ++k) hi[k] = int(std::lround(bufs[t.buf_idx][k]));
    t.f->upload(hi);
  }
  if (m_.suite) m_.suite->restarted(start_step);
}
#endif

// Output field lists and the wrfout writer (also called by the nest driver before start(), so that
// the nest initialisation feedback reaches the parent's first history frame)
void DomainRunner::prepare_output() {
  if (output_prepared_) return;
  output_prepared_ = true;
  if (start_date_.empty())
    start_date_ = cfg_.get_string("run.start_date", !m_.start_date.empty() ? m_.start_date : "2000-01-01_00:00:00");
  if (m_.rk && m_.rk->physics()) {
    extra_2d_ = {&m_.rk->physics()->rainnc, &m_.rk->physics()->rainncv};
    if (m_.rk->physics()->params().mp_physics == 6 || m_.rk->physics()->params().mp_physics == 8) {
      extra_2d_.push_back(&m_.rk->physics()->snownc);
      extra_2d_.push_back(&m_.rk->physics()->snowncv);
      extra_2d_.push_back(&m_.rk->physics()->graupelnc);
      extra_2d_.push_back(&m_.rk->physics()->graupelncv);
      extra_2d_.push_back(&m_.rk->physics()->sr);
    }
  }
  if (m_.suite) {  // column physics surface and cumulus fields (WRF names in lower case)
    const phys::driver::PhysicsSuite& su = *m_.suite;
    for (const Field2* f : {&su.hfx, &su.qfx, &su.lh, &su.ust, &su.znt, &su.pblh, &su.tsk, &su.grdflx, &su.qsfc,
                            &su.flhc, &su.flqc, &su.t2, &su.th2, &su.q2, &su.u10, &su.v10, &su.psfc, &su.albedo,
                            &su.albbck, &su.emiss, &su.canwat, &su.snow, &su.snowh, &su.snowc, &su.noahres,
                            &su.acgrdflx, &su.achfx, &su.aclhf, &su.sfcevp, &su.potevp, &su.mol, &su.rmol, &su.br,
                            &su.wspd, &su.fm, &su.fh, &su.zol, &su.regime, &su.sfcexc, &su.lai, &su.xland,
                            &su.rainc, &su.raincv, &su.nca, &su.pratec, &su.cutop, &su.cubot, &su.htop, &su.hbot})
      extra_2d_.push_back(f);
    extra_3d_.insert(extra_3d_.end(), {&su.tslb, &su.smois, &su.sh2o, &su.pt.rublten, &su.pt.rvblten, &su.pt.rthblten,
                                       &su.pt.rqvblten, &su.exch_h});
    if (su.options().sf_sfclay_physics == 2 || su.options().bl_pbl_physics == 2) {
      for (const Field2* f : {&su.thz0, &su.uz0, &su.vz0, &su.akhs, &su.akms, &su.ct, &su.tshltr, &su.th10, &su.qshltr,
                              &su.q10, &su.pshltr, &su.u10e, &su.v10e, &su.mixht})
        extra_2d_.push_back(f);
      extra_3d_.insert(extra_3d_.end(), {&su.tke_pbl, &su.el_pbl});
    }
    if (su.radiation_on()) {
      for (const Field2* f : {&su.glw, &su.gsw, &su.swdown, &su.coszen}) extra_2d_.push_back(f);
      extra_3d_.push_back(&su.pt.rthraten);
      extra_3d_.push_back(&su.cldfra);
      if (auto* pc = m_.rk->physics(); pc && pc->params().microp_re)
        extra_3d_.insert(extra_3d_.end(), {&pc->re_cloud, &pc->re_ice, &pc->re_snow});
    }
    if (su.options().cu_physics != 0)
    {
      extra_3d_.insert(extra_3d_.end(), {&su.w0avg, &su.pt.rthcuten, &su.pt.rqvcuten, &su.pt.rqccuten});
      if (m_.st->P_QR >= 0) extra_3d_.push_back(&su.pt.rqrcuten);
      if (m_.st->P_QI >= 0) extra_3d_.push_back(&su.pt.rqicuten);
      if (m_.st->P_QS >= 0) extra_3d_.push_back(&su.pt.rqscuten);
      if (su.options().cu_physics == 16 || su.options().cu_physics == 3)
        extra_3d_.insert(extra_3d_.end(), {&su.pt.rucuten, &su.pt.rvcuten, &su.pt.rthften, &su.pt.rqvften});
    }
    if (m_.st->P_QC >= 0) extra_3d_.push_back(&su.pt.rqcblten);
    if (m_.st->P_QI >= 0) extra_3d_.push_back(&su.pt.rqiblten);
  }

#ifdef CPPWRF_HAVE_NETCDF
  if (write_wrfout_ && !wrfout_)
  {
    std::vector<const Field2*> f2 = extra_2d_;
    std::vector<const Field3*> f3;  // first match wins: the microphysics' grid%rho is the last writer
    if (m_.rk && m_.rk->physics())
      for (const Field3* f : m_.rk->physics()->fields()) f3.push_back(f);
    f3.insert(f3.end(), extra_3d_.begin(), extra_3d_.end());
    if (m_.suite) {
      for (const Field2* f : m_.suite->fields2()) f2.push_back(f);
      for (const Field3* f : m_.suite->fields3()) f3.push_back(f);
    }
    if (m_.diff)
      for (const Field3* f : m_.diff->fields()) f3.push_back(f);
    wrfout_ = std::make_unique<WrfHistory>(cfg_, m_, out_dir_, domain_tag_, start_date_, f2, f3);
    if (m_.rk) m_.rk->on_first_rk_prep = [this] { wrfout_->first_rk_prep(); };
  }
#endif
}

bool DomainRunner::start() {
  const Grid& g = m_.g;
  dyn::State& st = *m_.st;
  dyn::BaseState& base = *m_.base;
  const dyn::LargeStepParams& lp = m_.lp;
  const dyn::AcousticParams& ap = m_.ap;
  const cases::CaseDef& cd = m_.cd;

  diag_iv_ = cfg_.get_int("output.diag_interval", 10);
  out_iv_ = cfg_.get_int("output.out_interval", std::max(1L, m_.nsteps));
  if (diag_iv_ < 1) throw std::runtime_error("output.diag_interval must be >= 1 (got " + std::to_string(diag_iv_) + ")");
  if (out_iv_ < 1) throw std::runtime_error("output.out_interval must be >= 1 (got " + std::to_string(out_iv_) + ")");

  if (g.dt <= 0.0) throw std::runtime_error("time.dt must be > 0 (got " + std::to_string(g.dt) + ")");
  if (g.nx <= 0) throw std::runtime_error("grid.nx must be > 0 (got " + std::to_string(g.nx) + ")");
  if (g.ny <= 0) throw std::runtime_error("grid.ny must be > 0 (got " + std::to_string(g.ny) + ")");
  if (g.nz <= 0) throw std::runtime_error("grid.nz must be > 0 (got " + std::to_string(g.nz) + ")");
  if (m_.run_time < 0.0) throw std::runtime_error("run.run_time must be >= 0 (got " + std::to_string(m_.run_time) + ")");
  if (m_.run_time > 0.0) {
    double steps_exact = m_.run_time / g.dt;
    if (std::abs(steps_exact - std::round(steps_exact)) > 1e-4) {
      throw std::runtime_error("run.run_time (" + std::to_string(m_.run_time) + ") must be an integer multiple of time.dt (" + std::to_string(g.dt) + ")");
    }
  }

  std::string restart_file = cfg_.get_string("run.restart_file", "");
  double restart_interval = cfg_.get_double("restart.interval", cfg_.get_double("output.restart_interval", 0.0));
  if (restart_interval < 0.0) throw std::runtime_error("restart.interval must be >= 0 (got " + std::to_string(restart_interval) + ")");
  rst_iv_ = restart_interval > 0 ? std::max(1L, std::lround(restart_interval / g.dt)) : 0;
  start_date_ = cfg_.get_string("run.start_date", !m_.start_date.empty() ? m_.start_date : "2000-01-01_00:00:00");

#ifndef CPPWRF_HAVE_NETCDF
  if (write_wrfout_) {
    throw std::runtime_error("output.format=" + out_fmt_ + " requires NetCDF support (CPPWRF_HAVE_NETCDF)");
  }
  if (rst_iv_ > 0 || !restart_file.empty()) {
    throw std::runtime_error("restarts require NetCDF support (CPPWRF_HAVE_NETCDF)");
  }
#endif

  prepare_output();  // reads the output.* keys
  for (const auto& k : cfg_.unaccessed_keys()) {
    std::fprintf(stderr, "cppwrf: warning: unaccessed config key '%s'\n", k.c_str());
  }

  start_step = 0;
#ifdef CPPWRF_HAVE_NETCDF
  if (!restart_file.empty()) {
    load_restart_file(restart_file);
    if (wrfout_) wrfout_->set_restart_time(double(start_step) * g.dt);
  }
#endif

  std::filesystem::create_directories(out_dir_);

  // run.json
  {
    io::JsonWriter j;
    j.start_object();
    j.key("grid"); j.start_object();
    j.key_value("nx", (long long)g.nx); j.key_value("ny", (long long)g.ny); j.key_value("nz", (long long)g.nz);
    j.key_value("dx", g.dx); j.key_value("dy", g.dy); j.key_value("dt", g.dt);
    j.key_value("ztop", g.ztop); j.key_value("p_top", g.p_top);
    j.key_value("bc_x", std::string_view(bc_name(g.bc_x))); j.key_value("bc_y", std::string_view(bc_name(g.bc_y)));
    j.key("eta_levels"); j.start_array(); for (double e : m_.eta) j.write_double(e); j.end_array();
    j.end_object();
    j.key_value("case", std::string_view(cname_));
    j.key("case_params"); j.start_object();
    for (auto& [k, v] : cd.params) j.key_value(k, v);
    j.end_object();
    j.key("params"); j.start_object();
    j.key_value("epssm", double(ap.epssm)); j.key_value("smdiv", double(ap.smdiv)); j.key_value("emdiv", double(ap.emdiv));
    j.key_value("khdif", double(lp.khdif)); j.key_value("kvdif", double(lp.kvdif));
    j.key_value("top_lid", lp.top_lid); j.key_value("time_step_sound", (long long)m_.n_sound);
    j.key_value("run_time", m_.run_time); j.key_value("nsteps", (long long)m_.nsteps);
    j.key_value("h_adv_order", (long long)lp.h_adv_order); j.key_value("v_adv_order", (long long)lp.v_adv_order);
    j.key_value("hybrid_opt", (long long)m_.hybrid_opt);
    j.key_value("diff_opt", (long long)m_.diff_opt);
    j.key_value("km_opt", (long long)m_.km_opt);
    j.key_value("pr_inv", double(lp.pr_inv));
    j.key_value("moist_adv_opt", (long long)m_.sp.moist_adv_opt); j.key_value("scalar_adv_opt", (long long)m_.sp.scalar_adv_opt);
    j.key("scalars"); j.start_array(); for (auto& n : st.q_names) j.write_string(n); j.end_array();
    j.key_value("n_moist", (long long)st.n_moist);
    j.key_value("mp_physics", (long long)(m_.rk && m_.rk->physics() ? m_.rk->physics()->params().mp_physics : 0));
    j.end_object();
    j.key_value("precision", std::string_view(precision_name));
    j.key_value("backend", std::string_view(backend_name()));
    j.key_value("small_step", dyn::RK3::small_step_available());
    j.key_value("config_raw", "__CPPWRF_CONFIG_JSON_RAW__");
    j.end_object();
    std::string js = j.str();
    std::string cj = cfg_.to_json();
    std::string token = "\"config_raw\": \"__CPPWRF_CONFIG_JSON_RAW__\"";
    auto pos = js.find(token);
    if (pos != std::string::npos) {
      js.replace(pos, token.size(), "\"config\": " + cj);
    }
    io::write_text_file_atomic(out_dir_ + "/run.json", js + "\n");
  }

  // Preserve pre-restart rows in diag.csv verbatim (M7)
  std::vector<std::string> pre_diag_lines;
  if (start_step > 0 && std::filesystem::exists(out_dir_ + "/diag.csv")) {
    std::ifstream f(out_dir_ + "/diag.csv");
    std::string line;
    if (std::getline(f, line)) {
      while (std::getline(f, line)) {
        if (line.empty()) continue;
        auto comma = line.find(',');
        if (comma != std::string::npos) {
          try {
            long step_val = std::stol(line.substr(0, comma));
            if (step_val < start_step) {
              pre_diag_lines.push_back(line);
            }
          } catch (...) {}
        }
      }
    }
  }

  diag_ = std::make_unique<io::DiagWriter>(out_dir_ + "/diag.csv",
                                           std::vector<std::string>{"step", "time", "max_abs_u", "max_abs_v", "max_abs_w", "min_th_pert", "max_th_pert",
                                                                    "dry_mass", "theta_mass", "total_energy", "wall_s"});
  for (const auto& line : pre_diag_lines) {
    diag_->raw_line(line);
  }

  snap_ = std::make_unique<io::SnapshotWriter>(out_dir_);
  dg_ = std::make_unique<dyn::Diagnostics>(g);

  prepare_output();

  t0_ = std::chrono::steady_clock::now();
  step_time = 0.0;

  std::printf("cppwrf: case=%s grid=%dx%dx%d dx=%g dt=%g n_sound=%d steps=%ld backend=%s precision=%s\n",
              cname_.c_str(), g.nx, g.ny, g.nz, g.dx, g.dt, m_.n_sound, m_.nsteps, backend_name(), precision_name);

  if (!do_diag(start_step)) return false;
  if (start_step == 0) do_snap(start_step);
  return true;
}

void DomainRunner::after_solve(long n) {
#ifdef CPPWRF_HAVE_NETCDF
  if (write_wrfout_ && wrfout_ && (n % out_iv_ == 0 || n == m_.nsteps)) wrfout_->after_solve(n);
#else
  (void)n;
#endif
}

void DomainRunner::before_step(long n) {
#ifdef CPPWRF_HAVE_NETCDF
  if (write_wrfout_ && wrfout_ && (n % out_iv_ == 0 || n == m_.nsteps || (rst_iv_ > 0 && n % rst_iv_ == 0)))
    wrfout_->before_step(n);
#endif
  if (pending_ >= 0) {  // deferred frame of the previous step (set_defer_output)
    const long p = pending_;
    pending_ = -1;
    if (p % diag_iv_ == 0) do_diag(p);
    if (p % out_iv_ == 0) do_snap(p);
  }
  step_start_ = std::chrono::steady_clock::now();
  if (m_.before_step) m_.before_step(n);
}

bool DomainRunner::after_step(long n) {
  step_time += std::chrono::duration<double>(std::chrono::steady_clock::now() - step_start_).count();
  bool last = (n == m_.nsteps);
  if (defer_output_ && !last && (n % diag_iv_ == 0 || n % out_iv_ == 0)) {
    pending_ = n;
  } else if (n % diag_iv_ == 0 || last) {
    if (!do_diag(n)) return false;
  }
  if ((n % out_iv_ == 0 || last) && pending_ != n) {
    do_snap(n);
  }
#ifdef CPPWRF_HAVE_NETCDF
  if (rst_iv_ > 0 && (n % rst_iv_ == 0 || last)) {
    do_restart(n);
  }
#endif
  return true;
}

void DomainRunner::finish() {
  if (m_.nsteps > start_step) {
    double cells = double(m_.g.nx) * m_.g.ny * m_.g.nz;
    long taken = m_.nsteps - start_step;
    std::printf("cppwrf: done, %ld steps, stepping %.3f s, %.2f ns/cell/step\n", taken, step_time,
                step_time / (cells * taken) * 1e9);
  } else {
    std::printf("cppwrf: initialised, no time stepping (run.run_time = 0)\n");
  }
}

bool DomainRunner::do_diag(long n) {
  dyn::DiagRow r = dg_->compute(*m_.st, *m_.base);
  double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0_).count();
  if (!r.finite()) {
    std::fprintf(stderr, "ERROR: non-finite values at step %ld (t=%g s): max|u|=%g max|w|=%g th'=[%g,%g]\n", n,
                 n * m_.g.dt, r.max_abs_u, r.max_abs_w, r.min_th_pert, r.max_th_pert);
    return false;
  }
  diag_->row({double(n), n * m_.g.dt, r.max_abs_u, r.max_abs_v, r.max_abs_w, r.min_th_pert, r.max_th_pert,
              r.dry_mass, r.theta_mass, r.total_energy, wall});
  return true;
}

void DomainRunner::do_snap(long n) {
  const Grid& g = m_.g;
  dyn::State& st = *m_.st;
  dyn::BaseState& base = *m_.base;

  if (write_native_) {
    auto vars = dyn::snapshot_vars(st, base, g, store_, extra_2d_, extra_3d_);
    snap_->write(n, n * g.dt, vars);
  }

#ifdef CPPWRF_HAVE_NETCDF
  if (write_wrfout_ && wrfout_) wrfout_->write(n, n * g.dt);
#endif
}

#ifdef CPPWRF_HAVE_NETCDF
void DomainRunner::do_restart(long n) {
  char rst_name[64];
  std::snprintf(rst_name, sizeof(rst_name), "restart_%06ld.nc", n);
  std::string rst_path = out_dir_ + "/" + rst_name;
  save_restart_file(rst_path, n, n * m_.g.dt);
  if (wrfout_) wrfout_->write_restart(n, n * m_.g.dt, double(out_iv_) * m_.g.dt, double(rst_iv_) * m_.g.dt);
}
#endif

}  // namespace cppwrf::cases
