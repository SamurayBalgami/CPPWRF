#include "cases/domain_runner.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "core/constants.hpp"
#include "core/halo.hpp"
#include "core/precision.hpp"
#include "dyn/diagnostics.hpp"

namespace cppwrf::cases {

namespace {

const char* bc_name(BC b) {
  return b == BC::Periodic ? "periodic" : (b == BC::Symmetric ? "symmetric" : "specified");
}

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
    if (rk->physics()->params().mp_physics == 6) {
      add_field2("snownc", rk->physics()->snownc);
      add_field2("snowncv", rk->physics()->snowncv);
      add_field2("graupelnc", rk->physics()->graupelnc);
      add_field2("graupelncv", rk->physics()->graupelncv);
      add_field2("sr", rk->physics()->sr);
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
    if (rk->physics()->params().mp_physics == 6) {
      reg_field2("snownc", rk->physics()->snownc);
      reg_field2("snowncv", rk->physics()->snowncv);
      reg_field2("graupelnc", rk->physics()->graupelnc);
      reg_field2("graupelncv", rk->physics()->graupelncv);
      reg_field2("sr", rk->physics()->sr);
    }
  }

  io_netcdf::read_restart(path, meta, vars);

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
}
#endif

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

  for (const auto& k : cfg_.unaccessed_keys()) {
    std::fprintf(stderr, "cppwrf: warning: unaccessed config key '%s'\n", k.c_str());
  }

  start_step = 0;
#ifdef CPPWRF_HAVE_NETCDF
  if (!restart_file.empty()) {
    load_restart_file(restart_file);
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

  if (m_.rk && m_.rk->physics()) {
    extra_2d_ = {&m_.rk->physics()->rainnc, &m_.rk->physics()->rainncv};
    if (m_.rk->physics()->params().mp_physics == 6) {
      extra_2d_.push_back(&m_.rk->physics()->snownc);
      extra_2d_.push_back(&m_.rk->physics()->snowncv);
      extra_2d_.push_back(&m_.rk->physics()->graupelnc);
      extra_2d_.push_back(&m_.rk->physics()->graupelncv);
      extra_2d_.push_back(&m_.rk->physics()->sr);
    }
  }

#ifdef CPPWRF_HAVE_NETCDF
  if (write_wrfout_) {
    io_netcdf::WrfoutMeta wm;
    wm.nx = g.nx; wm.ny = g.ny; wm.nz = g.nz;
    wm.dx = g.dx; wm.dy = g.dy; wm.dt = g.dt;
    wm.ztop = g.ztop; wm.p_top = g.p_top;
    wm.start_date = start_date_;
    wm.simulation_start_date = start_date_;
    wm.diff_opt = m_.diff_opt;
    wm.km_opt = m_.km_opt;
    wm.khdif = double(lp.khdif); wm.kvdif = double(lp.kvdif);
    wm.mp_physics = int(m_.rk && m_.rk->physics() ? m_.rk->physics()->params().mp_physics : 0);
    wm.hybrid_opt = m_.hybrid_opt;
    wm.map_proj = m_.map_proj;
    wm.cen_lat = m_.cen_lat; wm.cen_lon = m_.cen_lon;
    wm.truelat1 = m_.truelat1; wm.truelat2 = m_.truelat2;
    wm.moad_cen_lat = m_.moad_cen_lat; wm.stand_lon = m_.stand_lon;
    wm.is_restart = (start_step > 0);
    wrfout_ = std::make_unique<io_netcdf::WrfoutWriter>(out_dir_, wm);
  }
#endif

  t0_ = std::chrono::steady_clock::now();
  step_time = 0.0;

  std::printf("cppwrf: case=%s grid=%dx%dx%d dx=%g dt=%g n_sound=%d steps=%ld backend=%s precision=%s\n",
              cname_.c_str(), g.nx, g.ny, g.nz, g.dx, g.dt, m_.n_sound, m_.nsteps, backend_name(), precision_name);

  if (!do_diag(start_step)) return false;
  if (start_step == 0) do_snap(start_step);
  return true;
}

void DomainRunner::before_step(long n) {
  step_start_ = std::chrono::steady_clock::now();
  if (m_.before_step) m_.before_step(n);
}

bool DomainRunner::after_step(long n) {
  step_time += std::chrono::duration<double>(std::chrono::steady_clock::now() - step_start_).count();
  bool last = (n == m_.nsteps);
  if (n % diag_iv_ == 0 || last) {
    if (!do_diag(n)) return false;
  }
  if (n % out_iv_ == 0 || last) {
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
    auto vars = dyn::snapshot_vars(st, base, g, store_, extra_2d_);
    snap_->write(n, n * g.dt, vars);
  }

#ifdef CPPWRF_HAVE_NETCDF
  if (write_wrfout_ && wrfout_) {
    const int nx = g.nx, ny = g.ny, nz = g.nz;
    const int nxs_u = g.nxs(Stagger::U), nys_v = g.nys(Stagger::V);

    std::vector<real> hu, hv, hw, ht, hp, hph, hmu;
    st.u_2.download(hu); st.v_2.download(hv); st.w_2.download(hw); st.t_2.download(ht);
    st.p.download(hp); st.ph_2.download(hph); st.mu_2.download(hmu);

    std::vector<real> bpb, bphb, bmub, bht;
    base.pb.download(bpb); base.phb.download(bphb); base.mub.download(bmub); base.ht.download(bht);

    std::vector<real> hqv;
    if (st.P_QV >= 0) st.q_2[st.P_QV].download(hqv);

    const Layout &Lu = st.u_2.layout(), &Lv = st.v_2.layout(), &Lw = st.w_2.layout(), &Lm = st.t_2.layout(),
                 &Lph = st.ph_2.layout(), &Lmu = st.mu_2.layout(), &Lpb = base.pb.layout(), &Lphb = base.phb.layout(),
                 &Lmub = base.mub.layout(), &Lht = base.ht.layout();

    auto at = [](const std::vector<real>& h, const Layout& L, int i, int j, int k) {
      return double(h[L.index(i, j, k)]);
    };

    std::vector<float> u_buf(size_t(nxs_u) * ny * nz);
    std::vector<float> v_buf(size_t(nx) * nys_v * nz);
    std::vector<float> w_buf(size_t(nx) * ny * (nz + 1));
    std::vector<float> ph_buf(size_t(nx) * ny * (nz + 1));
    std::vector<float> phb_buf(size_t(nx) * ny * (nz + 1));
    std::vector<float> p_buf(size_t(nx) * ny * nz);
    std::vector<float> pb_buf(size_t(nx) * ny * nz);
    std::vector<float> t_buf(size_t(nx) * ny * nz);
    std::vector<float> mu_buf(size_t(nx) * ny);
    std::vector<float> mub_buf(size_t(nx) * ny);
    std::vector<float> hgt_buf(size_t(nx) * ny);

    for (int k = 0; k < nz; ++k)
      for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nxs_u; ++i)
          u_buf[(size_t(k) * ny + j) * nxs_u + i] = float(at(hu, Lu, i, j, k));

    for (int k = 0; k < nz; ++k)
      for (int j = 0; j < nys_v; ++j)
        for (int i = 0; i < nx; ++i)
          v_buf[(size_t(k) * nys_v + j) * nx + i] = float(at(hv, Lv, i, j, k));

    for (int k = 0; k <= nz; ++k)
      for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
          w_buf[(size_t(k) * ny + j) * nx + i] = float(at(hw, Lw, i, j, k));
          ph_buf[(size_t(k) * ny + j) * nx + i] = float(at(hph, Lph, i, j, k));
          phb_buf[(size_t(k) * ny + j) * nx + i] = float(at(bphb, Lphb, i, j, k));
        }

    for (int k = 0; k < nz; ++k)
      for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
          p_buf[(size_t(k) * ny + j) * nx + i] = float(at(hp, Lm, i, j, k));
          pb_buf[(size_t(k) * ny + j) * nx + i] = float(at(bpb, Lpb, i, j, k));
          double tm = at(ht, Lm, i, j, k) + cst::t0;
          double th = hqv.empty() ? tm : tm / (1.0 + cst::rvovrd * at(hqv, Lm, i, j, k));
          t_buf[(size_t(k) * ny + j) * nx + i] = float(th - 300.0);
        }

    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        mu_buf[size_t(j) * nx + i] = float(at(hmu, Lmu, i, j, 0));
        mub_buf[size_t(j) * nx + i] = float(at(bmub, Lmub, i, j, 0));
        hgt_buf[size_t(j) * nx + i] = float(at(bht, Lht, i, j, 0));
      }

    std::vector<io_netcdf::VarView> nc_vars;
    const io_netcdf::DType f32 = io_netcdf::DType::F32;
    nc_vars.emplace_back("U", 'U', nxs_u, ny, nz, u_buf.data(), f32, nxs_u, (std::ptrdiff_t)nxs_u * ny);
    nc_vars.emplace_back("V", 'V', nx, nys_v, nz, v_buf.data(), f32, nx, (std::ptrdiff_t)nx * nys_v);
    nc_vars.emplace_back("W", 'W', nx, ny, nz + 1, w_buf.data(), f32, nx, (std::ptrdiff_t)nx * ny);
    nc_vars.emplace_back("PH", 'W', nx, ny, nz + 1, ph_buf.data(), f32, nx, (std::ptrdiff_t)nx * ny);
    nc_vars.emplace_back("PHB", 'W', nx, ny, nz + 1, phb_buf.data(), f32, nx, (std::ptrdiff_t)nx * ny);
    nc_vars.emplace_back("P", 'M', nx, ny, nz, p_buf.data(), f32, nx, (std::ptrdiff_t)nx * ny);
    nc_vars.emplace_back("PB", 'M', nx, ny, nz, pb_buf.data(), f32, nx, (std::ptrdiff_t)nx * ny);
    nc_vars.emplace_back("T", 'M', nx, ny, nz, t_buf.data(), f32, nx, (std::ptrdiff_t)nx * ny);
    nc_vars.emplace_back("MU", 'M', nx, ny, 1, mu_buf.data(), f32, nx, 0);
    nc_vars.emplace_back("MUB", 'M', nx, ny, 1, mub_buf.data(), f32, nx, 0);
    nc_vars.emplace_back("HGT", 'M', nx, ny, 1, hgt_buf.data(), f32, nx, 0);

    std::vector<float> xlat_buf(size_t(nx) * ny, 0.0f);
    std::vector<float> xlong_buf(size_t(nx) * ny, 0.0f);
    if (!m_.xlat.empty() && m_.xlat.size() == size_t(nx) * ny) {
      for (size_t idx = 0; idx < m_.xlat.size(); ++idx) xlat_buf[idx] = float(m_.xlat[idx]);
    }
    if (!m_.xlong.empty() && m_.xlong.size() == size_t(nx) * ny) {
      for (size_t idx = 0; idx < m_.xlong.size(); ++idx) xlong_buf[idx] = float(m_.xlong[idx]);
    }
    nc_vars.emplace_back("XLAT", 'M', nx, ny, 1, xlat_buf.data(), f32, nx, 0);
    nc_vars.emplace_back("XLONG", 'M', nx, ny, 1, xlong_buf.data(), f32, nx, 0);

    std::vector<float> znu_buf(nz), znw_buf(nz + 1);
    std::vector<float> c1h_buf(nz), c2h_buf(nz), c3h_buf(nz), c4h_buf(nz);
    std::vector<float> c1f_buf(nz + 1), c2f_buf(nz + 1), c3f_buf(nz + 1), c4f_buf(nz + 1);
    for (int k = 0; k < nz; ++k) {
      znu_buf[k] = float(base.znu[k]);
      c1h_buf[k] = float(base.c1h[k]);
      c2h_buf[k] = float(base.c2h[k]);
      c3h_buf[k] = float(base.c3h[k]);
      c4h_buf[k] = float(base.c4h[k]);
    }
    for (int k = 0; k <= nz; ++k) {
      znw_buf[k] = float(base.znw[k]);
      c1f_buf[k] = float(base.c1f[k]);
      c2f_buf[k] = float(base.c2f[k]);
      c3f_buf[k] = float(base.c3f[k]);
      c4f_buf[k] = float(base.c4f[k]);
    }
    nc_vars.emplace_back("ZNU", 'M', 1, 1, nz, znu_buf.data(), f32, 0, 1);
    nc_vars.emplace_back("ZNW", 'W', 1, 1, nz + 1, znw_buf.data(), f32, 0, 1);
    nc_vars.emplace_back("C1H", 'M', 1, 1, nz, c1h_buf.data(), f32, 0, 1);
    nc_vars.emplace_back("C2H", 'M', 1, 1, nz, c2h_buf.data(), f32, 0, 1);
    nc_vars.emplace_back("C3H", 'M', 1, 1, nz, c3h_buf.data(), f32, 0, 1);
    nc_vars.emplace_back("C4H", 'M', 1, 1, nz, c4h_buf.data(), f32, 0, 1);
    nc_vars.emplace_back("C1F", 'W', 1, 1, nz + 1, c1f_buf.data(), f32, 0, 1);
    nc_vars.emplace_back("C2F", 'W', 1, 1, nz + 1, c2f_buf.data(), f32, 0, 1);
    nc_vars.emplace_back("C3F", 'W', 1, 1, nz + 1, c3f_buf.data(), f32, 0, 1);
    nc_vars.emplace_back("C4F", 'W', 1, 1, nz + 1, c4f_buf.data(), f32, 0, 1);

    std::vector<std::vector<float>> scalar_bufs(st.n_scalars());
    for (int s = 0; s < st.n_scalars(); ++s) {
      std::vector<real> hq;
      st.q_2[s].download(hq);
      const Layout& Lq = st.q_2[s].layout();
      scalar_bufs[s].resize(size_t(nx) * ny * nz);
      for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
          for (int i = 0; i < nx; ++i)
            scalar_bufs[s][(size_t(k) * ny + j) * nx + i] = float(at(hq, Lq, i, j, k));
      std::string qname = st.q_names[s];
      if (qname == "qv") qname = "QVAPOR";
      else if (qname == "qc") qname = "QCLOUD";
      else if (qname == "qr") qname = "QRAIN";
      else if (qname == "qi") qname = "QICE";
      else if (qname == "qs") qname = "QSNOW";
      else if (qname == "qg") qname = "QGRAUP";
      nc_vars.emplace_back(qname, 'M', nx, ny, nz, scalar_bufs[s].data(), f32, nx, (std::ptrdiff_t)nx * ny);
    }

    std::vector<std::vector<float>> extra_bufs(extra_2d_.size());
    for (size_t e = 0; e < extra_2d_.size(); ++e) {
      const Field2* f = extra_2d_[e];
      std::vector<real> he;
      f->download(he);
      const Layout& Le = f->layout();
      extra_bufs[e].resize(size_t(nx) * ny);
      for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
          extra_bufs[e][size_t(j) * nx + i] = float(at(he, Le, i, j, 0));
      std::string ename = f->name();
      for (char& c : ename) c = std::toupper(c);
      nc_vars.emplace_back(ename, 'M', nx, ny, 1, extra_bufs[e].data(), f32, nx, 0);
    }

    std::string target_file = out_dir_ + "/wrfout_" + domain_tag_ + "_" + start_date_;
    wrfout_->write_frame(n, n * g.dt, nc_vars, target_file);
  }
#endif
}

#ifdef CPPWRF_HAVE_NETCDF
void DomainRunner::do_restart(long n) {
  char rst_name[64];
  std::snprintf(rst_name, sizeof(rst_name), "restart_%06ld.nc", n);
  std::string rst_path = out_dir_ + "/" + rst_name;
  save_restart_file(rst_path, n, n * m_.g.dt);
}
#endif

}  // namespace cppwrf::cases
