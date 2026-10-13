#include "cases/model.hpp"
#include "cases/hill2d.hpp"
#include "cases/quarter_ss.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

#include "core/constants.hpp"
#include "core/halo.hpp"

namespace cppwrf::cases {

static BC parse_bc(const std::string& s) {
  if (s == "periodic") return BC::Periodic;
  if (s == "symmetric") return BC::Symmetric;
  if (s == "specified") return BC::Specified;
  throw std::runtime_error("unknown boundary condition '" + s + "' (periodic|symmetric|specified)");
}

std::vector<std::string> split_list(const std::string& s) {
  std::vector<std::string> out;
  std::string cur;
  auto flush = [&] {
    size_t a = cur.find_first_not_of(" \t"), b = cur.find_last_not_of(" \t");
    if (a != std::string::npos) out.push_back(cur.substr(a, b - a + 1));
    cur.clear();
  };
  for (char ch : s) {
    if (ch == ',') flush();
    else cur += ch;
  }
  flush();
  return out;
}

// Passive tracer initial condition: tracers.<name>_init = const (value <name>_value, default 1)
// or blob (cos^2 in r/radius, amplitude <name>_amp, centre <name>_xc/_zc, radius <name>_radius).
static void init_tracer(const io::Config& cfg, const std::string& n, const Grid& G, const dyn::BaseState& b,
                        Field3& q1, Field3& q2) {
  const std::string pre = "tracers." + n + "_";
  const std::string kind = cfg.get_string(pre + "init", "const");
  std::vector<real> phb;
  b.phb.download(phb);
  const Layout& L = b.phb.layout();
  auto zm = [&](int i, int j, int k) { return 0.5 * (double(phb[L.index(i, j, k)]) + double(phb[L.index(i, j, k + 1)])) / cst::g; };
  std::function<double(int, int, int)> f;
  if (kind == "const") {
    const double v = cfg.get_double(pre + "value", 1.0);
    f = [v](int, int, int) { return v; };
  } else if (kind == "blob") {
    const double xc = cfg.get_double(pre + "xc", 0.5 * G.nx * G.dx), zc = cfg.get_double(pre + "zc", 2000.0);
    const double r = cfg.get_double(pre + "radius", 2000.0), amp = cfg.get_double(pre + "amp", 1.0);
    f = [=](int i, int j, int k) {
      const double x = (i + 0.5) * G.dx, z = zm(i, j, k);
      const double d = std::sqrt(std::pow((x - xc) / r, 2) + std::pow((z - zc) / r, 2));
      return d <= 1.0 ? amp * std::pow(std::cos(0.5 * std::numbers::pi * d), 2) : 0.0;
    };
  } else {
    throw std::runtime_error(pre + "init must be const or blob");
  }
  q2.fill_host(f);
  fill_halo(q2, G);
  q1.copy_from(q2);
}

static void init_ideal(const io::Config& cfg, Model& m, const std::string& cname) {
  Grid& g = m.g;
  g.nx = int(cfg.get_int("domain.nx"));
  g.ny = int(cfg.get_int("domain.ny", 1));
  g.nz = int(cfg.get_int("domain.nz"));
  g.dx = cfg.get_double("domain.dx");
  g.dy = cfg.get_double("domain.dy", g.dx);
  g.ztop = cfg.get_double("domain.ztop");
  g.bc_x = parse_bc(cfg.get_string("domain.bc_x", "periodic"));
  g.bc_y = parse_bc(cfg.get_string("domain.bc_y", "periodic"));
  g.dt = cfg.get_double("time.dt");
  m.run_time = cfg.get_double("run.run_time", 0.0);
  m.nsteps = std::lround(m.run_time / g.dt);
  if (std::fabs(m.nsteps * g.dt - m.run_time) > 1e-9 * std::max(1.0, m.run_time))
    throw std::runtime_error("run.run_time must be a multiple of time.dt");

  // Test options: a constant map factor (physical spacing dx / m; the case is laid out in physical
  // distance) and an f-plane (domain.f0, domain.e0 in 1/s).
  const double mapf = cfg.get_double("domain.map_factor", 1.0);
  if (!(mapf > 0)) throw std::runtime_error("domain.map_factor must be > 0");
  Grid gp = g;
  gp.dx /= mapf;
  gp.dy /= mapf;
  m.cd = make_case(cname, cfg, gp);
  Sounding snd(m.cd.theta_env, m.cd.p_surf, g.ztop * 1.2 + 1000.0);
  // domain.p_top overrides the exact sounding pressure at ztop (to match a WRF run, whose ideal.exe
  // interpolates p_top linearly between input-sounding levels)
  g.p_top = cfg.get_double("domain.p_top", snd.p(g.ztop));
  m.eta = make_eta_levels(cfg, g, snd);
  const bool flat = (cname != "hill2d");
  m.base = std::make_unique<dyn::BaseState>(g, flat);
  if (mapf != 1.0) {
    for (BaseField2* f : {&m.base->msftx, &m.base->msfty, &m.base->msfux, &m.base->msfuy, &m.base->msfvx,
                          &m.base->msfvy})
      f->fill_host_halo([=](int, int, int) { return mapf; });
    m.base->msfvx_inv.fill_host_halo([=](int, int, int) { return 1.0 / mapf; });
  }
  const double f0 = cfg.get_double("domain.f0", 0.0), e0 = cfg.get_double("domain.e0", 0.0);
  if (f0 != 0.0 || e0 != 0.0) {
    m.base->f.fill_host_halo([=](int, int, int) { return f0; });
    m.base->e.fill_host_halo([=](int, int, int) { return e0; });
    m.base->coriolis = true;
  }
  m.st = std::make_unique<dyn::State>(g);
  m.tend = std::make_unique<dyn::Tendencies>(g);
  // scalars: moist species (moist.species), then passive tracers (tracers.names)
  for (const std::string& n : split_list(cfg.get_string("moist.species", ""))) {
    m.st->add_scalar(g, n, true);
    m.tend->add_scalar(g, n);
  }
  std::vector<std::string> tracers = split_list(cfg.get_string("tracers.names", ""));
  for (const std::string& n : tracers) {
    m.st->add_scalar(g, n, false);
    m.tend->add_scalar(g, n);
  }
  if (m.cd.qv_env && m.st->P_QV < 0) throw std::runtime_error("case.qv0 needs qv in moist.species");
  std::unique_ptr<Sounding> moist;
  if (m.st->P_QV >= 0) {  // moist sounding for the full state (qv = 0 if the case has no vapour)
    auto qv = m.cd.qv_env ? m.cd.qv_env : std::function<double(double)>([](double) { return 0.0; });
    moist = std::make_unique<Sounding>(m.cd.theta_env, m.cd.p_surf, g.ztop * 1.2 + 1000.0, 0.5, qv);
  }
  gp = g;  // now with p_top
  gp.dx /= mapf;
  gp.dy /= mapf;
  if (cname == "hill2d") {
    init_hill2d(cfg, m.eta, gp, *m.base, *m.st);
  } else if (cname == "quarter_ss") {
    init_quarter_ss(cfg, m.eta, gp, *m.base, *m.st);
  } else {
    initialize(m.cd, snd, m.eta, gp, *m.base, *m.st, moist.get());
  }
  if (!m.base->ht.uniform()) {  // start_em: kinematic surface w, filled up with znw^2 (no input w)
    dyn::set_w_surface(m.st->w_1, m.st->u_1, m.st->v_1, *m.base, g, true);
    dyn::set_w_surface(m.st->w_2, m.st->u_2, m.st->v_2, *m.base, g, true);
  }
  for (const std::string& n : tracers) {
    int idx = 0;
    while (m.st->q_names[idx] != n) ++idx;
    init_tracer(cfg, n, g, *m.base, m.st->q_1[idx], m.st->q_2[idx]);
  }

}

Model Model::build(const io::Config& cfg) {
  Model m;
  Grid& g = m.g;
  const std::string cname = cfg.get_string("run.case");
  if (cname == "real") init_real_case(cfg, m);
  else init_ideal(cfg, m, cname);
  dyn::LargeStepParams& lp = m.lp;
  lp.khdif = real(cfg.get_double("dyn.khdif", 0.0));
  lp.kvdif = real(cfg.get_double("dyn.kvdif", 0.0));
  lp.pr_inv = real(cfg.get_double("dyn.pr_inv", 3.0));
  lp.top_lid = cfg.get_bool("dyn.top_lid", false);  // WRF default
  lp.h_adv_order = int(cfg.get_int("dyn.h_adv_order", 5));
  lp.v_adv_order = int(cfg.get_int("dyn.v_adv_order", 3));
  lp.w_damping = int(cfg.get_int("dyn.w_damping", 0));
  lp.w_crit_cfl = real(cfg.get_double("dyn.w_crit_cfl", 1.0));
  if ((lp.h_adv_order != 5 && lp.h_adv_order != 6) || (lp.v_adv_order != 3 && lp.v_adv_order != 4))
    throw std::runtime_error("dyn.h_adv_order must be 5 or 6, dyn.v_adv_order 3 or 4");
  m.ap.epssm = real(cfg.get_double("dyn.epssm", 0.1));
  m.ap.smdiv = real(cfg.get_double("dyn.smdiv", 0.1));
  m.ap.emdiv = real(cfg.get_double("dyn.emdiv", 0.01));
  m.ap.top_lid = lp.top_lid;
  // WRF takes hypsometric_opt from the input file: ideal.exe writes 1, real.exe 2 (namelist default)
  if (cname != "real") m.ap.hypsometric_opt = int(cfg.get_int("dyn.hypsometric_opt", 1));  // real: from the input
  if (m.ap.hypsometric_opt != 1 && m.ap.hypsometric_opt != 2) throw std::runtime_error("dyn.hypsometric_opt: 1 or 2");
  m.ap.damp_opt = int(cfg.get_int("dyn.damp_opt", 0));
  if (m.ap.damp_opt != 0 && m.ap.damp_opt != 3) throw std::runtime_error("dyn.damp_opt must be 0 or 3");
  m.ap.zdamp = real(cfg.get_double("dyn.zdamp", 5000.0));
  m.ap.dampcoef = real(cfg.get_double("dyn.dampcoef", 0.2));
  dyn::ScalarParams& sp = m.sp;
  sp.moist_adv_opt = int(cfg.get_int("dyn.moist_adv_opt", 1));
  sp.scalar_adv_opt = int(cfg.get_int("dyn.scalar_adv_opt", 1));
  if (sp.moist_adv_opt < 0 || sp.moist_adv_opt > 1 || sp.scalar_adv_opt < 0 || sp.scalar_adv_opt > 1)
    throw std::runtime_error("dyn.moist_adv_opt and dyn.scalar_adv_opt must be 0 or 1");
  if (m.st->n_scalars() > 0 && (sp.moist_adv_opt == 1 || sp.scalar_adv_opt == 1) && lp.h_adv_order != 5)
    throw std::runtime_error("positive-definite advection is implemented for h_adv_order = 5 only");
  sp.khdif = lp.khdif; sp.kvdif = lp.kvdif; sp.pr_inv = lp.pr_inv;
  sp.h_adv_order = lp.h_adv_order; sp.v_adv_order = lp.v_adv_order;
  const int tss = int(cfg.get_int("time.time_step_sound", 0));
  m.n_sound = tss > 0 ? tss : dyn::RK3::auto_sound_steps(g.dt, g.dx, g.dy, g.two_d());
  if (m.nsteps > 0) {
    m.rk = std::make_unique<dyn::RK3>(g, lp, m.ap, tss, sp);
    m.n_sound = m.rk->n_sound();
    dyn::PhysicsParams pp;
    pp.mp_physics = int(cfg.get_int("phys.mp_physics", 0));
    pp.mp_tend_lim = real(cfg.get_double("phys.mp_tend_lim", 10.0));
    pp.hail_opt = int(cfg.get_int("phys.hail_opt", 0));
    if (pp.mp_physics != 0) m.rk->enable_physics(pp);
    if (m.bc) m.rk->set_lateral_bc(m.bc.get());
  }
  return m;
}

}  // namespace cppwrf::cases
