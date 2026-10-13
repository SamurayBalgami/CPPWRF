// Real-data start (WRF dyn_em/start_em.F for a single domain, no restart, use_theta_m = 1):
// grid, base state, map factors, Coriolis, terrain, initial state and lateral boundaries from
// real.exe's wrfinput_d01 / wrfbdy_d01 (read by io_netcdf readers).
//   - 1D vertical arrays (znw, znu, dnw, dn, rdnw, rdn, fnm, fnp, c1..c4, cf1..cfn) as stored.
//   - pb = c3h*mub + c4h + p_top; alb from t_init (start_em, max_dom = 1); phb, mub as stored.
//   - prognostic theta: THM = (T + t0)(1 + Rv/Rd qv) - t0 (Registry dataname "thm" for t).
//   - al from ph (hypsometric_opt from the input, 2 for real.exe), p from al and theta_m.
//   - w = 0 in the input -> set_w_surface with the znw^2 fill (start_em w_needs_to_be_set).
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>

#include "cases/model.hpp"
#include "core/constants.hpp"
#include "core/halo.hpp"
#include "dyn/lateral_bc.hpp"
#ifdef CPPWRF_HAVE_NETCDF
#include "io_netcdf/wrf_readers.hpp"
#endif

namespace cppwrf::cases {

#ifdef CPPWRF_HAVE_NETCDF
namespace {
// "YYYY-MM-DD_HH:MM:SS" -> seconds between two times (proleptic Gregorian, UTC)
double seconds_between(const std::string& a, const std::string& b) {
  auto secs = [](const std::string& s) {
    int Y, M, D, h, mi, se;
    if (std::sscanf(s.c_str(), "%d-%d-%d_%d:%d:%d", &Y, &M, &D, &h, &mi, &se) != 6) return 0.0;
    const int y = Y - (M <= 2), era = (y >= 0 ? y : y - 399) / 400, yoe = y - era * 400;
    const int doy = (153 * (M + (M > 2 ? -3 : 9)) + 2) / 5 + D - 1, doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    const double days = double(era) * 146097 + doe - 719468;
    return days * 86400.0 + h * 3600.0 + mi * 60.0 + se;
  };
  return secs(b) - secs(a);
}
}  // namespace
#endif

#ifndef CPPWRF_HAVE_NETCDF
void init_real_case(const io::Config&, Model&) { throw std::runtime_error("case = real needs NetCDF"); }
#else
void init_real_case(const io::Config& cfg, Model& m) {
  const std::string dir = cfg.get_string("real.dir");
  const io_netcdf::WrfInputData d = io_netcdf::read_wrfinput(dir + "/" + cfg.get_string("real.wrfinput", "wrfinput_d01"));
  Grid& g = m.g;
  g.nx = d.nx; g.ny = d.ny; g.nz = d.nz;
  g.dx = d.dx; g.dy = d.dy;
  g.dt = cfg.get_double("time.dt", d.dt);
  g.p_top = d.p_top;
  g.bc_x = g.bc_y = BC::Specified;
  m.run_time = cfg.get_double("run.run_time", 0.0);
  m.nsteps = std::lround(m.run_time / g.dt);
  const int nx = g.nx, ny = g.ny, nz = g.nz;

  // ---- base state and coefficients
  m.base = std::make_unique<dyn::BaseState>(g, false, true);
  dyn::BaseState& b = *m.base;
  for (int k = 0; k <= nz; ++k) { b.znw[k] = d.znw[k]; b.c1f[k] = d.c1f[k]; b.c2f[k] = d.c2f[k]; b.c3f[k] = d.c3f[k]; b.c4f[k] = d.c4f[k]; }
  for (int k = 0; k < nz; ++k) {
    b.znu[k] = d.znu[k]; b.dnw[k] = d.dnw[k]; b.rdnw[k] = d.rdnw[k];
    b.c1h[k] = d.c1h[k]; b.c2h[k] = d.c2h[k]; b.c3h[k] = d.c3h[k]; b.c4h[k] = d.c4h[k];
    // w-level arrays: file index = Fortran - 1 = C++ w level (dn(1) = 0 in WRF)
    b.dn[k] = d.dn[k]; b.rdn[k] = d.rdn[k]; b.fnm[k] = d.fnm[k]; b.fnp[k] = d.fnp[k];
  }
  b.dn[nz] = 0; b.rdn[nz] = 0; b.fnm[nz] = 0; b.fnp[nz] = 0;
  b.cf1 = d.cf1; b.cf2 = d.cf2; b.cf3 = d.cf3; b.cfn = d.cfn; b.cfn1 = d.cfn1;
  b.upload_levels();
  b.p_top = d.p_top;

  auto cl = [](int a, int lo, int hi) { return a < lo ? lo : (a > hi ? hi : a); };
  auto m2 = [&](const std::vector<double>& v) {  // mass-point 2D, flat halo
    return [&v, nx, ny, cl](int i, int j, int) { return v[std::size_t(cl(j, 0, ny - 1)) * nx + cl(i, 0, nx - 1)]; };
  };
  b.mub.fill_host_halo(m2(d.mub));
  b.ht.fill_host_halo(m2(d.hgt));
  b.msftx.fill_host_halo(m2(d.mapfac_mx));
  b.msfty.fill_host_halo(m2(d.mapfac_my));
  b.f.fill_host_halo(m2(d.f));
  b.e.fill_host_halo(m2(d.e));
  b.sina.fill_host_halo(m2(d.sinalpha));
  b.cosa.fill_host_halo(m2(d.cosalpha));
  b.coriolis = true;
  const int nxs = d.nx_stag, nys = d.ny_stag;
  auto u2 = [&](const std::vector<double>& v) {
    return [&v, nxs, ny, cl](int i, int j, int) { return v[std::size_t(cl(j, 0, ny - 1)) * nxs + cl(i, 0, nxs - 1)]; };
  };
  auto v2 = [&](const std::vector<double>& v) {
    return [&v, nx, nys, cl](int i, int j, int) { return v[std::size_t(cl(j, 0, nys - 1)) * nx + cl(i, 0, nx - 1)]; };
  };
  b.msfux.fill_host_halo(u2(d.mapfac_ux));
  b.msfuy.fill_host_halo(u2(d.mapfac_uy));
  b.msfvx.fill_host_halo(v2(d.mapfac_vx));
  b.msfvy.fill_host_halo(v2(d.mapfac_vy));
  b.msfvx_inv.fill_host_halo([&](int i, int j, int) { return 1.0 / d.mapfac_vx[std::size_t(cl(j, 0, nys - 1)) * nx + cl(i, 0, nx - 1)]; });

  auto m3 = [&](const std::vector<double>& v, int nk) {
    return [&v, nx, ny, nk, cl](int i, int j, int k) {
      return v[(std::size_t(cl(k, 0, nk - 1)) * ny + cl(j, 0, ny - 1)) * nx + cl(i, 0, nx - 1)];
    };
  };
  // start_em (max_dom = 1): pb from mub, alb from t_init
  std::vector<double> pb(std::size_t(nx) * ny * nz), alb(pb.size());
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        const std::size_t n = (std::size_t(k) * ny + j) * nx + i;
        pb[n] = d.c3h[k] * d.mub[std::size_t(j) * nx + i] + d.c4h[k] + d.p_top;
        alb[n] = (cst::r_d / cst::p1000mb) * (d.t_init[n] + cst::t0) * std::pow(pb[n] / cst::p1000mb, cst::cvpm);
      }
  b.pb.fill_host_halo(m3(pb, nz));
  b.alb.fill_host_halo(m3(alb, nz));
  b.t_init.fill_host_halo(m3(d.t_init, nz));
  b.phb.fill_host_halo(m3(d.phb, nz + 1));
  b.qv_base.fill_host_halo([](int, int, int) { return 0.0; });

  // ---- state
  m.st = std::make_unique<dyn::State>(g);
  m.tend = std::make_unique<dyn::Tendencies>(g);
  dyn::State& s = *m.st;
  for (const std::string& n : split_list(cfg.get_string("moist.species", "qv"))) {
    s.add_scalar(g, n, true);
    m.tend->add_scalar(g, n);
  }
  if (s.P_QV < 0) throw std::runtime_error("case = real needs qv in moist.species");
  const std::vector<double>& qv = d.qvapor;
  std::vector<double> thm(qv.size());
  for (std::size_t n = 0; n < thm.size(); ++n) thm[n] = (d.t[n] + cst::t0) * (1.0 + cst::rvovrd * qv[n]) - cst::t0;
  auto idx3 = [&](int i, int j, int k, int ni, int nj) { return (std::size_t(k) * nj + j) * ni + i; };
  s.u_2.fill_host([&](int i, int j, int k) { return d.u[idx3(i, j, k, nxs, ny)]; });
  s.v_2.fill_host([&](int i, int j, int k) { return d.v[idx3(i, j, k, nx, nys)]; });
  s.w_2.fill_host([](int, int, int) { return 0.0; });
  s.t_2.fill_host([&](int i, int j, int k) { return thm[idx3(i, j, k, nx, ny)]; });
  s.ph_2.fill_host([&](int i, int j, int k) { return d.ph[idx3(i, j, k, nx, ny)]; });
  s.mu_2.fill_host([&](int i, int j, int) { return d.mu[std::size_t(j) * nx + i]; });
  s.q_2[s.P_QV].fill_host([&](int i, int j, int k) { return qv[idx3(i, j, k, nx, ny)]; });
  // al (hypsometric_opt) and p from start_em, theta_m with qvf = 1
  const int hypso = int(cfg.get_int("dyn.hypsometric_opt", d.attrs_int.count("HYPSOMETRIC_OPT")
                                                                ? d.attrs_int.at("HYPSOMETRIC_OPT") : 2));
  std::vector<double> al(pb.size()), p(pb.size());
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        const std::size_t n = idx3(i, j, k, nx, ny), n1 = idx3(i, j, k + 1, nx, ny);
        const double mu = d.mu[std::size_t(j) * nx + i], mt = d.mub[std::size_t(j) * nx + i] + mu;
        double a;
        if (hypso == 2) {
          const double pfu = d.c3f[k + 1] * mt + d.c4f[k + 1] + d.p_top, pfd = d.c3f[k] * mt + d.c4f[k] + d.p_top;
          const double phm = d.c3h[k] * mt + d.c4h[k] + d.p_top;
          a = (d.ph[n1] - d.ph[n] + d.phb[n1] - d.phb[n]) / phm / std::log(pfd / pfu) - alb[n];
        } else {
          a = -1.0 / (d.c1h[k] * mt + d.c2h[k]) * (alb[n] * (d.c1h[k] * mu) + d.rdnw[k] * (d.ph[n1] - d.ph[n]));
        }
        al[n] = a;
        p[n] = cst::p1000mb * std::pow((cst::r_d * (cst::t0 + thm[n])) / (cst::p1000mb * (a + alb[n])), cst::cpovcv) - pb[n];
      }
  s.al.fill_host([&](int i, int j, int k) { return al[idx3(i, j, k, nx, ny)]; });
  s.p.fill_host([&](int i, int j, int k) { return p[idx3(i, j, k, nx, ny)]; });
  fill_halos({&s.u_2, &s.v_2, &s.w_2, &s.t_2, &s.ph_2, &s.al, &s.p, &s.q_2[s.P_QV]}, g);
  fill_halo(s.mu_2, g);
  s.u_1.copy_from(s.u_2); s.v_1.copy_from(s.v_2); s.t_1.copy_from(s.t_2); s.ph_1.copy_from(s.ph_2);
  s.mu_1.copy_from(s.mu_2); s.q_1[s.P_QV].copy_from(s.q_2[s.P_QV]);
  dyn::set_w_surface(s.w_1, s.u_1, s.v_1, b, g, true);
  dyn::set_w_surface(s.w_2, s.u_2, s.v_2, b, g, true);
  m.ap.hypsometric_opt = hypso;

  // ---- lateral boundaries: wrfbdy record r covers [t_r, t_r + interval); the record needed for
  // a step is loaded before it, and dtbc is the time since the start of that interval.
  auto bd = std::make_shared<io_netcdf::WrfBdyData>(
      io_netcdf::read_wrfbdy(dir + "/" + cfg.get_string("real.wrfbdy", "wrfbdy_d01")));
  if (bd->bdy_width != 5) throw std::runtime_error("real: spec_bdy_width must be 5");
  m.bc = std::make_unique<dyn::LateralBC>(g, g.dt, bd->bdy_width);
  double interval = cfg.get_double("real.bdy_interval", 0.0);
  if (interval <= 0 && !bd->this_bdy_times.empty() && !bd->next_bdy_times.empty())
    interval = seconds_between(bd->this_bdy_times[0], bd->next_bdy_times[0]);
  if (interval <= 0) throw std::runtime_error("real: cannot determine the boundary interval (set real.bdy_interval)");
  dyn::LateralBC* bc = m.bc.get();
  auto load = [bd, bc, nx, ny, nz, nxs, nys](int rec) {
    using LB = dyn::LateralBC;
    if (rec >= bd->n_times) throw std::runtime_error("real: run longer than the boundary file");
    auto put = [&](LB::Var v, int npx, int npy, int nk, const std::vector<double>* const (&a)[8]) {
      const int npos[4] = {npx, npx, npy, npy};
      for (int sd = 0; sd < 4; ++sd) {
        const std::size_t n = std::size_t(npos[sd]) * nk * bd->bdy_width;
        auto cut = [&](const std::vector<double>& all) {
          return std::vector<double>(all.begin() + std::ptrdiff_t(rec * n), all.begin() + std::ptrdiff_t((rec + 1) * n));
        };
        bc->set(v, LB::Side(sd), cut(*a[sd]), cut(*a[4 + sd]));
      }
    };
    const io_netcdf::WrfBdyData& d = *bd;
    { const std::vector<double>* a[8] = {&d.u_bxs, &d.u_bxe, &d.u_bys, &d.u_bye, &d.u_btxs, &d.u_btxe, &d.u_btys, &d.u_btye};
      put(LB::U, ny, nxs, nz, a); }
    { const std::vector<double>* a[8] = {&d.v_bxs, &d.v_bxe, &d.v_bys, &d.v_bye, &d.v_btxs, &d.v_btxe, &d.v_btys, &d.v_btye};
      put(LB::V, nys, nx, nz, a); }
    { const std::vector<double>* a[8] = {&d.t_bxs, &d.t_bxe, &d.t_bys, &d.t_bye, &d.t_btxs, &d.t_btxe, &d.t_btys, &d.t_btye};
      put(LB::T, ny, nx, nz, a); }
    { const std::vector<double>* a[8] = {&d.ph_bxs, &d.ph_bxe, &d.ph_bys, &d.ph_bye, &d.ph_btxs, &d.ph_btxe, &d.ph_btys, &d.ph_btye};
      put(LB::PH, ny, nx, nz + 1, a); }
    { const std::vector<double>* a[8] = {&d.mu_bxs, &d.mu_bxe, &d.mu_bys, &d.mu_bye, &d.mu_btxs, &d.mu_btxe, &d.mu_btys, &d.mu_btye};
      put(LB::MU, ny, nx, 1, a); }
    { const std::vector<double>* a[8] = {&d.qvapor_bxs, &d.qvapor_bxe, &d.qvapor_bys, &d.qvapor_bye,
                                         &d.qvapor_btxs, &d.qvapor_btxe, &d.qvapor_btys, &d.qvapor_btye};
      put(LB::QV, ny, nx, nz, a); }
  };
  load(0);
  auto rec_now = std::make_shared<int>(0);
  const double dt = g.dt;
  m.before_step = [load, bc, rec_now, interval, dt](long n) {
    const double t0 = double(n - 1) * dt;  // model time at the start of step n
    const int need = int(std::floor(t0 / interval + 1e-9));
    if (need != *rec_now) { load(need); *rec_now = need; }
    bc->dtbc = t0 - need * interval;  // RK3::step adds dt (solve_em: dtbc = dtbc + dt)
  };
}
#endif

}  // namespace cppwrf::cases
