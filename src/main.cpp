// cppwrf: ./cppwrf <config.ini> [--set section.key=value ...]
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "cases/model.hpp"
#include "core/backend.hpp"
#include "core/halo.hpp"
#include "core/precision.hpp"
#include "dyn/diagnostics.hpp"
#include "dyn/rk3.hpp"
#include "io/config.hpp"
#include "io/output.hpp"

#ifdef CPPWRF_HAVE_NETCDF
#include "io_netcdf/restart.hpp"
#include "io_netcdf/wrfout_writer.hpp"
#endif

using namespace cppwrf;

static const char* bc_name(BC b) {
  return b == BC::Periodic ? "periodic" : (b == BC::Symmetric ? "symmetric" : "specified");
}

#ifdef CPPWRF_HAVE_NETCDF
static void save_restart(const std::string& path, long step, double time_sec,
                         const std::string& case_name, const std::string& date_str,
                         const Grid& g, dyn::State& st, dyn::RK3* rk) {
  io_netcdf::RestartMeta meta;
  meta.step = step;
  meta.time_sec = time_sec;
  meta.nx = g.nx; meta.ny = g.ny; meta.nz = g.nz;
  meta.dx = g.dx; meta.dy = g.dy; meta.dt = g.dt;
  meta.case_name = case_name;
  meta.date_str = date_str;

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

static void load_restart(const std::string& path, const Grid& g, dyn::State& st,
                         dyn::RK3* rk, io_netcdf::RestartMeta& meta) {
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

static int run(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <config.ini> [--set section.key=value ...]\n", argv[0]);
    return 2;
  }
  std::vector<std::string> ov;
  for (int a = 2; a < argc; ++a) {
    std::string s = argv[a];
    if (s == "--set" && a + 1 < argc) ov.push_back(argv[++a]);
    else if (s.rfind("--set=", 0) == 0) ov.push_back(s.substr(6));
    else throw std::runtime_error("unknown argument: " + s);
  }
  io::Config cfg = io::Config::load(argv[1], ov);
  const std::string cname = cfg.get_string("run.case");

  cases::Model M = cases::Model::build(cfg);
  const Grid& g = M.g;
  dyn::BaseState& base = *M.base;
  dyn::State& st = *M.st;
  dyn::Tendencies& tend = *M.tend;
  const dyn::LargeStepParams& lp = M.lp;
  const dyn::AcousticParams& ap = M.ap;
  const double run_time = M.run_time;
  const long nsteps = M.nsteps;
  const int n_sound = M.n_sound;
  const std::vector<double>& eta = M.eta;
  const cases::CaseDef& cd = M.cd;
  dyn::RK3* rk = M.rk.get();
  const long diag_iv = cfg.get_int("output.diag_interval", 10);
  const long out_iv = cfg.get_int("output.out_interval", std::max(1L, nsteps));
  const std::string out_dir = cfg.get_string("output.dir", "out/" + cname);

  const std::string out_fmt = cfg.get_string("output.format", "native");
  const bool write_native = (out_fmt == "native" || out_fmt == "both");
  const bool write_wrfout = (out_fmt == "wrfout" || out_fmt == "both");
  if (out_fmt != "native" && out_fmt != "wrfout" && out_fmt != "both") {
    throw std::runtime_error("unknown output.format: '" + out_fmt + "' (native|wrfout|both)");
  }
  const std::string restart_file = cfg.get_string("run.restart_file", "");
  const double restart_interval = cfg.get_double("restart.interval", cfg.get_double("output.restart_interval", 0.0));
  const long rst_iv = restart_interval > 0 ? std::max(1L, std::lround(restart_interval / g.dt)) : 0;
  const std::string start_date = cfg.get_string("run.start_date", "2000-01-01_00:00:00");

#ifndef CPPWRF_HAVE_NETCDF
  if (write_wrfout) {
    throw std::runtime_error("output.format=" + out_fmt + " requires NetCDF support (CPPWRF_HAVE_NETCDF)");
  }
  if (rst_iv > 0 || !restart_file.empty()) {
    throw std::runtime_error("restarts require NetCDF support (CPPWRF_HAVE_NETCDF)");
  }
#endif

  long start_step = 0;
#ifdef CPPWRF_HAVE_NETCDF
  if (!restart_file.empty()) {
    io_netcdf::RestartMeta rmeta;
    load_restart(restart_file, g, st, rk, rmeta);
    start_step = rmeta.step;
    if (nsteps < start_step) {
      throw std::runtime_error("run.run_time (" + std::to_string(run_time) + "s / step " + std::to_string(nsteps) +
                               ") is less than restart file step " + std::to_string(start_step));
    }
  }
#endif

  // run.json
  std::filesystem::create_directories(out_dir);
  {
    io::JsonWriter j;
    j.start_object();
    j.key("grid"); j.start_object();
    j.key_value("nx", (long long)g.nx); j.key_value("ny", (long long)g.ny); j.key_value("nz", (long long)g.nz);
    j.key_value("dx", g.dx); j.key_value("dy", g.dy); j.key_value("dt", g.dt);
    j.key_value("ztop", g.ztop); j.key_value("p_top", g.p_top);
    j.key_value("bc_x", std::string_view(bc_name(g.bc_x))); j.key_value("bc_y", std::string_view(bc_name(g.bc_y)));
    j.key("eta_levels"); j.start_array(); for (double e : eta) j.write_double(e); j.end_array();
    j.end_object();
    j.key_value("case", std::string_view(cname));
    j.key("case_params"); j.start_object();
    for (auto& [k, v] : cd.params) j.key_value(k, v);
    j.end_object();
    j.key("params"); j.start_object();
    j.key_value("epssm", double(ap.epssm)); j.key_value("smdiv", double(ap.smdiv)); j.key_value("emdiv", double(ap.emdiv));
    j.key_value("khdif", double(lp.khdif)); j.key_value("kvdif", double(lp.kvdif));
    j.key_value("top_lid", lp.top_lid); j.key_value("time_step_sound", (long long)n_sound);
    j.key_value("run_time", run_time); j.key_value("nsteps", (long long)nsteps);
    j.key_value("h_adv_order", (long long)lp.h_adv_order); j.key_value("v_adv_order", (long long)lp.v_adv_order); j.key_value("hybrid_opt", 0LL);
    j.key_value("pr_inv", double(lp.pr_inv));
    j.key_value("moist_adv_opt", (long long)M.sp.moist_adv_opt); j.key_value("scalar_adv_opt", (long long)M.sp.scalar_adv_opt);
    j.key("scalars"); j.start_array(); for (auto& n : st.q_names) j.write_string(n); j.end_array();
    j.key_value("n_moist", (long long)st.n_moist);
    j.key_value("mp_physics", (long long)(M.rk && M.rk->physics() ? M.rk->physics()->params().mp_physics : 0));
    j.end_object();
    j.key_value("precision", std::string_view(precision_name));
    j.key_value("backend", std::string_view(backend_name()));
    j.key_value("small_step", dyn::RK3::small_step_available());
    j.key("config"); j.start_object(); j.end_object();
    j.end_object();
    std::string js = j.str();
    // embed the full resolved config (already JSON) in place of the empty object
    std::string cj = cfg.to_json();
    auto pos = js.rfind("{}");
    if (pos != std::string::npos) js.replace(pos, 2, cj);
    io::write_text_file_atomic(out_dir + "/run.json", js + "\n");
  }

  // Preserve pre-restart rows in diag.csv if continuing in the same directory
  std::vector<std::vector<double>> pre_diag_rows;
  if (start_step > 0 && std::filesystem::exists(out_dir + "/diag.csv")) {
    std::ifstream f(out_dir + "/diag.csv");
    std::string line;
    if (std::getline(f, line)) {
      while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string val;
        std::vector<double> r;
        while (std::getline(ss, val, ',')) {
          r.push_back(std::stod(val));
        }
        if (!r.empty() && r[0] < start_step) {
          pre_diag_rows.push_back(r);
        }
      }
    }
  }

  io::DiagWriter diag(out_dir + "/diag.csv",
                      {"step", "time", "max_abs_u", "max_abs_v", "max_abs_w", "min_th_pert", "max_th_pert",
                       "dry_mass", "theta_mass", "total_energy", "wall_s"});
  for (const auto& r : pre_diag_rows) {
    diag.row(r);
  }

  io::SnapshotWriter snap(out_dir);
  dyn::Diagnostics dg(g);
  std::vector<std::vector<real>> store;
  auto t0 = std::chrono::steady_clock::now();
  double step_time = 0.0;
  auto do_diag = [&](long n) {
    dyn::DiagRow r = dg.compute(st, base);
    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (!r.finite()) {
      std::fprintf(stderr, "ERROR: non-finite values at step %ld (t=%g s): max|u|=%g max|w|=%g th'=[%g,%g]\n", n,
                   n * g.dt, r.max_abs_u, r.max_abs_w, r.min_th_pert, r.max_th_pert);
      return false;
    }
    diag.row({double(n), n * g.dt, r.max_abs_u, r.max_abs_v, r.max_abs_w, r.min_th_pert, r.max_th_pert, r.dry_mass,
              r.theta_mass, r.total_energy, wall});
    return true;
  };
  std::vector<const Field2*> extra_2d;
  if (rk && rk->physics()) {
    extra_2d = {&rk->physics()->rainnc, &rk->physics()->rainncv};
    if (rk->physics()->params().mp_physics == 6) {
      extra_2d.push_back(&rk->physics()->snownc);
      extra_2d.push_back(&rk->physics()->snowncv);
      extra_2d.push_back(&rk->physics()->graupelnc);
      extra_2d.push_back(&rk->physics()->graupelncv);
      extra_2d.push_back(&rk->physics()->sr);
    }
  }

#ifdef CPPWRF_HAVE_NETCDF
  std::unique_ptr<io_netcdf::WrfoutWriter> wrfout;
  std::vector<real> bpb, bphb, bmub;
  if (write_wrfout) {
    io_netcdf::WrfoutMeta wm;
    wm.nx = g.nx; wm.ny = g.ny; wm.nz = g.nz;
    wm.dx = g.dx; wm.dy = g.dy; wm.dt = g.dt;
    wm.ztop = g.ztop; wm.p_top = g.p_top;
    wm.start_date = start_date;
    wm.simulation_start_date = start_date;
    wm.diff_opt = int(lp.khdif > 0 ? 1 : 0);
    wm.km_opt = 1;
    wm.khdif = double(lp.khdif); wm.kvdif = double(lp.kvdif);
    wm.mp_physics = int(rk && rk->physics() ? rk->physics()->params().mp_physics : 0);
    wrfout = std::make_unique<io_netcdf::WrfoutWriter>(out_dir, wm);

    base.pb.download(bpb);
    base.phb.download(bphb);
    base.mub.download(bmub);
  }
#endif

  auto do_snap = [&](long n) {
    auto vars = dyn::snapshot_vars(st, base, g, store, extra_2d);
    if (write_native) {
      snap.write(n, n * g.dt, vars);
    }
#ifdef CPPWRF_HAVE_NETCDF
    if (write_wrfout && wrfout) {
      std::vector<io_netcdf::VarView> nc_vars;
      for (const auto& v : vars) {
        nc_vars.emplace_back(v.name, v.stagger, v.nx, v.ny, v.nz, v.data,
                             v.dtype == io::DType::F64 ? io_netcdf::DType::F64 : io_netcdf::DType::F32,
                             v.sj, v.sk);
      }
      const io_netcdf::DType nc_dt = sizeof(real) == 8 ? io_netcdf::DType::F64 : io_netcdf::DType::F32;
      nc_vars.emplace_back("PB", 'M', 1, 1, g.nz, bpb.data(), nc_dt, 0, 1);
      nc_vars.emplace_back("PHB", 'W', 1, 1, g.nz + 1, bphb.data(), nc_dt, 0, 1);
      nc_vars.emplace_back("MUB", 'M', 1, 1, 1, bmub.data(), nc_dt, 0, 0);
      wrfout->write_frame(n, n * g.dt, nc_vars);
    }
#endif
  };

#ifdef CPPWRF_HAVE_NETCDF
  auto do_restart = [&](long n) {
    char rst_name[64];
    std::snprintf(rst_name, sizeof(rst_name), "restart_%06ld.nc", n);
    std::string rst_path = out_dir + "/" + rst_name;
    std::string cur_date = io_netcdf::format_wrf_date(start_date, n * g.dt);
    save_restart(rst_path, n, n * g.dt, cname, cur_date, g, st, rk);
  };
#endif

  std::printf("cppwrf: case=%s grid=%dx%dx%d dx=%g dt=%g n_sound=%d steps=%ld backend=%s precision=%s\n",
              cname.c_str(), g.nx, g.ny, g.nz, g.dx, g.dt, n_sound, nsteps, backend_name(), precision_name);
  if (!do_diag(start_step)) return 3;
  do_snap(start_step);
  for (long n = start_step + 1; n <= nsteps; ++n) {
    auto a = std::chrono::steady_clock::now();
    rk->step(st, tend, base);
    dev_sync();
    step_time += std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
    bool last = n == nsteps;
    if (n % diag_iv == 0 || last)
      if (!do_diag(n)) return 3;
    if (n % out_iv == 0 || last) do_snap(n);
#ifdef CPPWRF_HAVE_NETCDF
    if (rst_iv > 0 && (n % rst_iv == 0 || last)) do_restart(n);
#endif
  }
  if (nsteps > start_step) {
    double cells = double(g.nx) * g.ny * g.nz;
    long taken = nsteps - start_step;
    std::printf("cppwrf: done, %ld steps, stepping %.3f s, %.2f ns/cell/step\n", taken, step_time,
                step_time / (cells * taken) * 1e9);
  } else {
    std::printf("cppwrf: initialised, no time stepping (run.run_time = 0)\n");
  }
  return 0;
}

int main(int argc, char** argv) {
  try {
    return run(argc, argv);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "cppwrf: error: %s\n", e.what());
    return 2;
  }
}
