#include "cases/model.hpp"
#include "cases/wrf_ideal.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

#include "core/constants.hpp"
#include "core/halo.hpp"

namespace cppwrf::cases {

BC parse_bc(const std::string& s) {
  if (s == "periodic") return BC::Periodic;
  if (s == "symmetric") return BC::Symmetric;
  if (s == "specified") return BC::Specified;
  if (s == "open") return BC::Open;
  throw std::runtime_error("unknown boundary condition '" + s + "' (periodic|symmetric|specified|open)");
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
  m.base = std::make_unique<dyn::BaseState>(g, true);
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
  if (cfg.get_int("dyn.km_opt", 1) == 2 && m.st->P_TKE < 0) {
    m.st->add_scalar(g, "tke", false);
    m.tend->add_scalar(g, "tke");
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
  initialize(m.cd, snd, m.eta, gp, *m.base, *m.st, moist.get());
  if (!m.base->ht.uniform()) {  // start_em: kinematic surface w, filled up with znw^2 (no input w)
    dyn::set_w_surface(m.st->w_1, m.st->u_1, m.st->v_1, *m.base, g, true);
    dyn::set_w_surface(m.st->w_2, m.st->u_2, m.st->v_2, *m.base, g, true);
  }
  for (const std::string& n : tracers) {
    int idx = 0;
    while (m.st->q_names[idx] != n) ++idx;
    init_tracer(cfg, n, g, *m.base, m.st->q_1[idx], m.st->q_2[idx]);
  }
  // module_initialize_ideal: u_base(k) = u_1(1,k,1), v_base(k) = v_1(1,k,1), t_base, qv_base, z_base (all ideal cases)
  std::vector<real> h;
  m.st->u_2.download(h);
  for (int k = 0; k < g.nz; ++k) m.base->u_base[k] = h[m.st->u_2.layout().index(0, 0, k)];
  m.st->v_2.download(h);
  for (int k = 0; k < g.nz; ++k) m.base->v_base[k] = h[m.st->v_2.layout().index(0, 0, k)];
  m.base->t_init.download(h);
  for (int k = 0; k < g.nz; ++k) m.base->t_base[k] = h[m.base->t_init.layout().index(0, 0, k)];
  m.base->qv_base.download(h);
  for (int k = 0; k < g.nz; ++k) m.base->qv_base_1d[k] = h[m.base->qv_base.layout().index(0, 0, k)];
  m.base->phb.download(h);
  for (int k = 0; k < g.nz; ++k)
    m.base->z_base[k] = real(0.5) * (h[m.base->phb.layout().index(0, 0, k)] + h[m.base->phb.layout().index(0, 0, k + 1)]) / cst::g;
  m.base->upload_levels();
}

// use_adaptive_time_step: start_em.F:932-1000 before phy_init (the physics take STEPRA, STEPBL, STEPCU from this dt)
static void init_adaptive(const io::Config& cfg, Model& m, const std::string& cname) {
  if (cname != "real") throw std::runtime_error("time.adaptive: the adaptive time step is ported for real-data runs");
  Grid& g = m.g;
  AdaptOptions o;
  o.target_cfl = real(cfg.get_double("time.target_cfl", 1.2));
  o.target_hcfl = real(cfg.get_double("time.target_hcfl", 0.84));
  o.max_step_increase_pct = int(cfg.get_int("time.max_step_increase_pct", 5));
  o.starting_time_step = int(cfg.get_int("time.starting_time_step", -1));
  o.starting_time_step_den = int(cfg.get_int("time.starting_time_step_den", 0));
  o.max_time_step = int(cfg.get_int("time.max_time_step", -1));
  o.max_time_step_den = int(cfg.get_int("time.max_time_step_den", 0));
  o.min_time_step = int(cfg.get_int("time.min_time_step", -1));
  o.min_time_step_den = int(cfg.get_int("time.min_time_step_den", 0));
  o.step_to_output_time = cfg.get_bool("time.step_to_output_time", true);
  // d01's history interval in seconds (history_interval_d/h/m/s, m = history_interval when 0)
  o.history_interval = 0;
  o.history_interval_m = 0;
  o.history_interval_s = int(cfg.get_int("time.history_interval_s", 3600));
  if (o.step_to_output_time && o.history_interval_s <= 0)
    throw std::runtime_error("time.adaptive with step_to_output_time needs time.history_interval_s > 0");
  o.interval_seconds = int(cfg.get_int("time.interval_seconds", 21600));
  if (o.interval_seconds <= 0) throw std::runtime_error("time.interval_seconds must be > 0");
  m.adapt = std::make_unique<AdaptTimeStep>(o);
  // MAXVAL(grid%msftx), MAXVAL(grid%msfty): the input's map factors (the memory halo of a serial WRF is 0)
  auto maxval = [](const BaseField2& f) {
    std::vector<real> h;
    f.download(h);
    const Layout& L = f.layout();
    if (f.uniform()) return h[L.index(0, 0, 0)];
    real mx = h[L.index(0, 0, 0)];
    for (int j = 0; j < L.ny; ++j)
      for (int i = 0; i < L.nx; ++i) mx = std::max(mx, h[L.index(i, j, 0)]);
    return mx;
  };
  const bool restart = !cfg.get_string("run.restart_file", "").empty();
  m.adapt->start_em(0, g.dx, g.dy, maxval(m.base->msftx), maxval(m.base->msfty), restart);
  g.dt = double(m.adapt->dt());
  m.run_length = wrf_time::make(std::llround(m.run_time));
  m.nsteps = m.run_time > 0 ? std::numeric_limits<long>::max() / 4 : 0;  // set when the last step is known
  if (m.bc) m.bc->set_dt(g.dt);  // start_em's lbc_fcx_gcx
}

void Model::set_dt(double dt) {
  g.dt = dt;
  n_sound = dyn::RK3::adaptive_sound_steps(real(dt), real(g.dx), real(g.dy), std::max(adapt->max_msftx(), adapt->max_msfty()));
  if (rk) rk->set_dt(dt, n_sound);
  if (bc) bc->set_dt(dt);
  if (diff) diff->set_dt(dt);
  if (suite) suite->set_dt(dt);
}

void Model::adapt_step(long n) {
  const wrf_time::Interval elapsed = sim_offset + now;
  if (before_step_time) before_step_time(elapsed);  // med_before_solve_io: med_latbound_in
  adapt->adapt(n - 1, now, real(bc ? bc->dtbc : 0.0));  // grid%itimestep before solve_em's increment
  set_dt(double(adapt->dt()));
  adapt->reset_cfl();
  step_end = now + adapt->time_step();
  if (suite) {
    suite->adapt_step_flag = true;
    suite->elapsed = elapsed;
  }
}

Model Model::build(const io::Config& cfg) {
  Model m;
  Grid& g = m.g;
  const std::string cname = cfg.get_string("run.case");
  const bool wrf_ideal = is_wrf_ideal_case(cname);
  if (cname == "real") init_real_case(cfg, m);
  else if (wrf_ideal) init_wrf_ideal(cfg, m, cname);
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
  // real and WRF ideal cases: from the input (init_*); ideal.exe writes the namelist value, default 2
  if (cname != "real" && !wrf_ideal) m.ap.hypsometric_opt = int(cfg.get_int("dyn.hypsometric_opt", 1));
  if (m.ap.hypsometric_opt != 1 && m.ap.hypsometric_opt != 2) throw std::runtime_error("dyn.hypsometric_opt: 1 or 2");
  m.ap.damp_opt = int(cfg.get_int("dyn.damp_opt", 0));
  if (m.ap.damp_opt < 0 || m.ap.damp_opt > 3) throw std::runtime_error("dyn.damp_opt must be 0, 1, 2 or 3");
  m.ap.zdamp = real(cfg.get_double("dyn.zdamp", 5000.0));
  m.ap.dampcoef = real(cfg.get_double("dyn.dampcoef", 0.2));
  dyn::ScalarParams& sp = m.sp;
  sp.moist_adv_opt = int(cfg.get_int("dyn.moist_adv_opt", 1));
  sp.scalar_adv_opt = int(cfg.get_int("dyn.scalar_adv_opt", 1));
  sp.tke_adv_opt = int(cfg.get_int("dyn.tke_adv_opt", 1));
  if (sp.moist_adv_opt < 0 || sp.moist_adv_opt > 1 || sp.scalar_adv_opt < 0 || sp.scalar_adv_opt > 1 ||
      sp.tke_adv_opt < 0 || sp.tke_adv_opt > 1)
    throw std::runtime_error("dyn.moist_adv_opt, dyn.scalar_adv_opt and dyn.tke_adv_opt must be 0 or 1");
  // positive-definite advection with h_adv_order 6: WRF's order-6 branch of advect_scalar_pd is its
  // order-5 branch with flux6 / flux4 for flux5 / flux3, i.e. the PD fluxes with hup = 0 (scalars.cpp)
  // Diffusion: WRF diff_opt / km_opt / diff_6th_opt. diff_opt defaults to 1 when khdif or kvdif
  // is set (older configs); diff_opt = 0 is no diffusion whatever khdif / kvdif (WRF).
  dyn::DiffusionConfig dc;
  dc.diff_opt = int(cfg.get_int("dyn.diff_opt", (lp.khdif > 0 || lp.kvdif > 0) ? 1 : 0));
  dc.km_opt = int(cfg.get_int("dyn.km_opt", 1));
  // restart files carry BN2, which WRF evaluates for every km_opt
  dc.n2_diagnostic = cfg.get_double("restart.interval", cfg.get_double("output.restart_interval", 0.0)) > 0;
  m.diff_opt = dc.diff_opt;
  m.km_opt = dc.km_opt;
  // real: from the input file; input_wrf: "hybrid_opt values must be consistent" (fatal)
  if (cname == "real" && cfg.has("dyn.hybrid_opt") && int(cfg.get_int("dyn.hybrid_opt")) != m.hybrid_opt)
    throw std::runtime_error("hybrid_opt values must be consistent: input file hybrid_opt = " +
                             std::to_string(m.hybrid_opt) + ", configuration dyn.hybrid_opt = " +
                             std::to_string(cfg.get_int("dyn.hybrid_opt")));
  m.hybrid_opt = int(cfg.get_int("dyn.hybrid_opt", m.hybrid_opt));
  dc.damp_opt = m.ap.damp_opt;
  dc.zdamp = m.ap.zdamp;
  dc.dampcoef = m.ap.dampcoef;
  dc.khdif = lp.khdif; dc.kvdif = lp.kvdif;
  dc.c_s = real(cfg.get_double("dyn.c_s", 0.25));
  dc.c_k = real(cfg.get_double("dyn.c_k", 0.15));
  dc.mix_isotropic = int(cfg.get_int("dyn.mix_isotropic", 0));
  dc.mix_full_fields = cfg.get_bool("dyn.mix_full_fields", false);
  dc.mix_upper_bound = real(cfg.get_double("dyn.mix_upper_bound", 0.1));
  dc.prandtl = real(1.0) / lp.pr_inv;
  dc.tke_drag_coefficient = real(cfg.get_double("dyn.tke_drag_coefficient", 0.0));
  dc.tke_heat_flux = real(cfg.get_double("dyn.tke_heat_flux", 0.0));
  dc.tke_upper_bound = real(cfg.get_double("dyn.tke_upper_bound", 1000.0));
  dc.moist_mix2_off = cfg.get_bool("dyn.moist_mix2_off", false);
  dc.scalar_mix2_off = cfg.get_bool("dyn.scalar_mix2_off", false);
  dc.tke_mix2_off = cfg.get_bool("dyn.tke_mix2_off", false);
  dc.moist_mix6_off = cfg.get_bool("dyn.moist_mix6_off", false);
  dc.scalar_mix6_off = cfg.get_bool("dyn.scalar_mix6_off", false);
  dc.tke_mix6_off = cfg.get_bool("dyn.tke_mix6_off", false);
  dc.diff_6th_opt = int(cfg.get_int("dyn.diff_6th_opt", 0));
  dc.diff_6th_factor = real(cfg.get_double("dyn.diff_6th_factor", 0.12));
  dc.pbl = cfg.get_int("phys.bl_pbl_physics", 0) != 0;
  dc.isfflx = int(cfg.get_int("phys.isfflx", 1));
  if (dc.isfflx != 0 && dc.isfflx != 1) throw std::runtime_error("phys.isfflx must be 0 or 1");
  dc.use_theta_m = int(cfg.get_int("dyn.use_theta_m", 1));
  if (dc.use_theta_m != 1) throw std::runtime_error("dyn.use_theta_m must be 1 (use_theta_m = 0 is not verified against WRF)");
  if (cfg.has("dyn.legacy_diffusion"))
    throw std::runtime_error("dyn.legacy_diffusion was removed; use dyn.diff_opt = 1 (WRF horizontal/vertical_diffusion)");
  sp.h_adv_order = lp.h_adv_order; sp.v_adv_order = lp.v_adv_order;
  if (cfg.get_bool("time.adaptive", false)) init_adaptive(cfg, m, cname);
  // the adaptive time step sets time_step_sound = 0 (start_em.F:979)
  const int tss = m.adapt ? 0 : int(cfg.get_int("time.time_step_sound", 0));
  m.n_sound = tss > 0 ? tss : dyn::RK3::auto_sound_steps(g.dt, g.dx, g.dy);
  if (m.nsteps > 0) {
    m.rk = std::make_unique<dyn::RK3>(g, lp, m.ap, tss, sp);
    m.n_sound = m.rk->n_sound();
    if (m.adapt) {
      m.n_sound = dyn::RK3::adaptive_sound_steps(real(g.dt), real(g.dx), real(g.dy),
                                                 std::max(m.adapt->max_msftx(), m.adapt->max_msfty()));
      m.rk->set_dt(g.dt, m.n_sound);
      AdaptTimeStep* a = m.adapt.get();
      m.rk->on_cfl = [a](real vert, real horiz) { a->add_cfl(vert, horiz); };
    }
    dyn::PhysicsParams pp;
    pp.mp_physics = int(cfg.get_int("phys.mp_physics", 0));
    pp.mp_tend_lim = real(cfg.get_double("phys.mp_tend_lim", 10.0));
    pp.hail_opt = int(cfg.get_int("phys.hail_opt", 0));
    pp.table_dir = cfg.get_string("phys.table_dir", cfg.get_string("real.dir", ""));
    pp.force_read_thompson = cfg.get_bool("phys.force_read_thompson", false);
    // phy_init: has_reqc = has_reqi = has_reqs = 1 with use_mp_re = 1, RRTMG LW + SW and WSM6 or Thompson
    pp.microp_re = cfg.get_int("phys.use_mp_re", 1) == 1 && cfg.get_int("phys.ra_lw_physics", 0) == 4 &&
                   cfg.get_int("phys.ra_sw_physics", 0) == 4 && (pp.mp_physics == 6 || pp.mp_physics == 8);
    if (pp.mp_physics != 0) m.rk->enable_physics(pp);
    if (wants_column_physics(cfg)) {
      if (cname != "real") throw std::runtime_error("phys.sf_*/bl_pbl_physics: real-data runs only");
      init_real_physics(cfg, m);
    } else if (cfg.get_int("phys.ra_lw_physics", 0) != 0 || cfg.get_int("phys.ra_sw_physics", 0) != 0) {
      // radiation is part of the column physics suite: without it ra_* would be silently ignored
      throw std::runtime_error("phys.ra_lw/ra_sw_physics without surface-layer, land-surface, PBL or cumulus physics "
                               "is not supported");
    }
    if (m.bc) m.rk->set_lateral_bc(m.bc.get());
    if (dc.on()) {
      m.diff = std::make_unique<dyn::DiffusionDriver>(g, dc);
      if (m.suite) m.diff->set_physics_suite(m.suite.get());
      m.rk->set_diffusion(m.diff.get());
    }
  }
  return m;
}

}  // namespace cppwrf::cases
