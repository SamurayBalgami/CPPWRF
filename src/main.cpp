// cppwrf: ./cppwrf <config.ini> [--set section.key=value ...]
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "cases/cases.hpp"
#include "core/backend.hpp"
#include "core/precision.hpp"
#include "dyn/diagnostics.hpp"
#include "dyn/rk3.hpp"
#include "io/config.hpp"
#include "io/output.hpp"

using namespace cppwrf;

static BC parse_bc(const std::string& s) {
  if (s == "periodic") return BC::Periodic;
  if (s == "symmetric") return BC::Symmetric;
  throw std::runtime_error("unknown boundary condition '" + s + "' (periodic|symmetric)");
}
static const char* bc_name(BC b) { return b == BC::Periodic ? "periodic" : "symmetric"; }

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

  Grid g;
  g.nx = int(cfg.get_int("domain.nx"));
  g.ny = int(cfg.get_int("domain.ny", 1));
  g.nz = int(cfg.get_int("domain.nz"));
  g.dx = cfg.get_double("domain.dx");
  g.dy = cfg.get_double("domain.dy", g.dx);
  g.ztop = cfg.get_double("domain.ztop");
  g.bc_x = parse_bc(cfg.get_string("domain.bc_x", "periodic"));
  g.bc_y = parse_bc(cfg.get_string("domain.bc_y", "periodic"));
  g.dt = cfg.get_double("time.dt");
  const double run_time = cfg.get_double("run.run_time", 0.0);
  const long nsteps = std::lround(run_time / g.dt);
  if (std::fabs(nsteps * g.dt - run_time) > 1e-9 * std::max(1.0, run_time))
    throw std::runtime_error("run.run_time must be a multiple of time.dt");
  const long diag_iv = cfg.get_int("output.diag_interval", 10);
  const long out_iv = cfg.get_int("output.out_interval", std::max(1L, nsteps));
  const std::string out_dir = cfg.get_string("output.dir", "out/" + cname);

  // case, sounding, vertical levels, base state, initial state
  cases::CaseDef cd = cases::make_case(cname, cfg, g);
  cases::Sounding snd(cd.theta_env, cd.p_surf, g.ztop * 1.2 + 1000.0);
  g.p_top = snd.p(g.ztop);
  std::vector<double> eta = cases::make_eta_levels(cfg, g, snd);
  dyn::BaseState base(g, true);
  dyn::State st(g);
  dyn::Tendencies tend(g);
  cases::initialize(cd, snd, eta, g, base, st);

  dyn::LargeStepParams lp;
  lp.khdif = real(cfg.get_double("dyn.khdif", 0.0));
  lp.kvdif = real(cfg.get_double("dyn.kvdif", 0.0));
  lp.top_lid = cfg.get_bool("dyn.top_lid", true);
  lp.h_adv_order = int(cfg.get_int("dyn.h_adv_order", 5));
  lp.v_adv_order = int(cfg.get_int("dyn.v_adv_order", 3));
  if ((lp.h_adv_order != 5 && lp.h_adv_order != 6) || (lp.v_adv_order != 3 && lp.v_adv_order != 4))
    throw std::runtime_error("dyn.h_adv_order must be 5 or 6, dyn.v_adv_order 3 or 4");
  dyn::AcousticParams ap;
  ap.epssm = real(cfg.get_double("dyn.epssm", 0.1));
  ap.smdiv = real(cfg.get_double("dyn.smdiv", 0.1));
  ap.emdiv = real(cfg.get_double("dyn.emdiv", 0.01));
  const int tss = int(cfg.get_int("time.time_step_sound", 0));
  std::unique_ptr<dyn::RK3> rk;
  int n_sound = tss > 0 ? tss : dyn::RK3::auto_sound_steps(g.dt, g.dx, g.dy, g.two_d());
  if (nsteps > 0) {
    rk = std::make_unique<dyn::RK3>(g, lp, ap, tss);
    n_sound = rk->n_sound();
  }

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

  io::DiagWriter diag(out_dir + "/diag.csv",
                      {"step", "time", "max_abs_u", "max_abs_v", "max_abs_w", "min_th_pert", "max_th_pert",
                       "dry_mass", "theta_mass", "total_energy", "wall_s"});
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
  auto do_snap = [&](long n) { snap.write(n, n * g.dt, dyn::snapshot_vars(st, base, g, store)); };

  std::printf("cppwrf: case=%s grid=%dx%dx%d dx=%g dt=%g n_sound=%d steps=%ld backend=%s precision=%s\n",
              cname.c_str(), g.nx, g.ny, g.nz, g.dx, g.dt, n_sound, nsteps, backend_name(), precision_name);
  if (!do_diag(0)) return 3;
  do_snap(0);
  for (long n = 1; n <= nsteps; ++n) {
    auto a = std::chrono::steady_clock::now();
    rk->step(st, tend, base);
    dev_sync();
    step_time += std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
    bool last = n == nsteps;
    if (n % diag_iv == 0 || last)
      if (!do_diag(n)) return 3;
    if (n % out_iv == 0 || last) do_snap(n);
  }
  if (nsteps > 0) {
    double cells = double(g.nx) * g.ny * g.nz;
    std::printf("cppwrf: done, %ld steps, stepping %.3f s, %.2f ns/cell/step\n", nsteps, step_time,
                step_time / (cells * nsteps) * 1e9);
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
