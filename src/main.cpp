// cppwrf: ./cppwrf <config.ini> [--set section.key=value ...]
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "cases/model.hpp"
#include "core/backend.hpp"
#include "core/precision.hpp"
#include "dyn/diagnostics.hpp"
#include "dyn/rk3.hpp"
#include "io/config.hpp"
#include "io/output.hpp"

using namespace cppwrf;

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
  std::vector<const Field2*> extra_2d;
  if (rk && rk->physics()) extra_2d = {&rk->physics()->rainnc, &rk->physics()->rainncv};
  auto do_snap = [&](long n) { snap.write(n, n * g.dt, dyn::snapshot_vars(st, base, g, store, extra_2d)); };

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
