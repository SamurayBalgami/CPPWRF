#include "cases/cases.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

#include "core/constants.hpp"
#include "core/halo.hpp"

namespace cppwrf::cases {

using namespace cppwrf::cst;

Sounding::Sounding(std::function<double(double)> theta, double p_surf, double zmax, double dz)
    : th_(std::move(theta)), dz_(dz) {
  int n = int(std::ceil(zmax / dz)) + 2;
  pi_.resize(n);
  pi_[0] = std::pow(p_surf / p1000mb, rcp);
  // d(pi)/dz = -g / (cp theta), integrated with Simpson's rule on each interval.
  for (int m = 1; m < n; ++m) {
    double z0 = (m - 1) * dz, z1 = m * dz;
    double f = (1.0 / th_(z0) + 4.0 / th_(0.5 * (z0 + z1)) + 1.0 / th_(z1)) / 6.0;
    pi_[m] = pi_[m - 1] - g / cp * f * dz;
  }
}

double Sounding::p(double z) const {
  double s = z / dz_;
  int m = int(std::floor(s));
  if (m < 0 || m + 1 >= int(pi_.size())) throw std::runtime_error("Sounding::p out of range");
  // integrate exactly from the tabulated point (Simpson on the partial interval)
  double z0 = m * dz_, h = z - z0;
  double f = (1.0 / th_(z0) + 4.0 / th_(z0 + 0.5 * h) + 1.0 / th_(z)) / 6.0;
  double pi = pi_[m] - g / cp * f * h;
  return p1000mb * std::pow(pi, 1.0 / rcp);
}

double Sounding::z_of_p(double pp) const {
  double lo = 0, hi = (pi_.size() - 2) * dz_;
  if (pp > p(lo) || pp < p(hi)) throw std::runtime_error("Sounding::z_of_p out of range");
  for (int it = 0; it < 200 && hi - lo > 1e-9; ++it) {
    double mid = 0.5 * (lo + hi);
    (p(mid) > pp ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}

static double xcenter(const Grid& G) { return 0.5 * G.nx * G.dx; }

CaseDef make_case(const std::string& name, const io::Config& cfg, const Grid& G) {
  CaseDef c;
  c.name = name;
  c.u_env = [](double) { return 0.0; };
  c.theta_pert = [](double, double, double, double) { return 0.0; };
  const double pi = std::numbers::pi;
  if (name == "rest_state") {
    // Stably stratified atmosphere at rest: constant N.
    double th0 = cfg.get_double("case.theta0", 300.0), N = cfg.get_double("case.n_bv", 0.01);
    c.theta_env = [=](double z) { return th0 * std::exp(N * N * z / g); };
    c.params = {{"theta0", th0}, {"n_bv", N}};
  } else if (name == "straka") {
    // Straka et al. (1993): neutral 300 K, cold bubble dT = -15 K (temperature perturbation).
    double xc = cfg.get_double("case.xc", xcenter(G)), zc = cfg.get_double("case.zc", 3000.0);
    double xr = cfg.get_double("case.xr", 4000.0), zr = cfg.get_double("case.zr", 2000.0);
    double dT = cfg.get_double("case.dt_bubble", -15.0);
    c.theta_env = [](double) { return 300.0; };
    c.theta_pert = [=](double x, double, double z, double exner) {
      double L = std::sqrt(std::pow((x - xc) / xr, 2) + std::pow((z - zc) / zr, 2));
      return L <= 1.0 ? dT * 0.5 * (std::cos(pi * L) + 1.0) / exner : 0.0;
    };
    c.params = {{"xc", xc}, {"zc", zc}, {"xr", xr}, {"zr", zr}, {"dt_bubble", dT}};
  } else if (name == "warm_bubble") {
    // Bryan & Fritsch (2002) dry thermal: neutral 300 K, theta' = 2 cos^2(pi L / 2).
    double xc = cfg.get_double("case.xc", xcenter(G)), zc = cfg.get_double("case.zc", 2000.0);
    double r = cfg.get_double("case.radius", 2000.0), dth = cfg.get_double("case.dtheta", 2.0);
    c.theta_env = [](double) { return 300.0; };
    c.theta_pert = [=](double x, double, double z, double) {
      double L = std::sqrt(std::pow((x - xc) / r, 2) + std::pow((z - zc) / r, 2));
      return L <= 1.0 ? dth * std::pow(std::cos(0.5 * pi * L), 2) : 0.0;
    };
    c.params = {{"xc", xc}, {"zc", zc}, {"radius", r}, {"dtheta", dth}};
  } else if (name == "igw") {
    // Skamarock & Klemp (1994) inertia-gravity wave (here without rotation, f = 0).
    double th0 = 300.0, N = cfg.get_double("case.n_bv", 0.01), U = cfg.get_double("case.u_mean", 20.0);
    double dth = cfg.get_double("case.dtheta", 0.01), a = cfg.get_double("case.a", 5000.0);
    double xc = cfg.get_double("case.xc", 100000.0), H = G.ztop;
    c.theta_env = [=](double z) { return th0 * std::exp(N * N * z / g); };
    c.u_env = [=](double) { return U; };
    c.theta_pert = [=](double x, double, double z, double) {
      return dth * std::sin(pi * z / H) / (1.0 + std::pow((x - xc) / a, 2));
    };
    c.params = {{"theta0", th0}, {"n_bv", N}, {"u_mean", U}, {"dtheta", dth}, {"a", a}, {"xc", xc}};
  } else {
    throw std::runtime_error("unknown case: " + name);
  }
  c.p_surf = cfg.get_double("case.p_surf", 1.0e5);
  return c;
}

std::vector<double> make_eta_levels(const io::Config& cfg, const Grid& G, const Sounding& snd) {
  std::string mode = cfg.get_string("domain.eta", "uniform_z");
  std::vector<double> e(G.nz + 1);
  if (mode == "uniform_z") {
    double ps = snd.p(0.0), pt = snd.p(G.ztop);
    for (int k = 0; k <= G.nz; ++k) e[k] = (snd.p(G.ztop * k / G.nz) - pt) / (ps - pt);
    e[0] = 1.0;
    e[G.nz] = 0.0;
  } else if (mode == "uniform_eta") {
    for (int k = 0; k <= G.nz; ++k) e[k] = 1.0 - double(k) / G.nz;
  } else if (mode == "list") {
    e = cfg.get_double_list("domain.eta_levels");
  } else {
    throw std::runtime_error("domain.eta must be uniform_z, uniform_eta or list");
  }
  return e;
}

void initialize(const CaseDef& c, const Sounding& snd, const std::vector<double>& eta, const Grid& G,
                dyn::BaseState& b, dyn::State& s) {
  const int nz = G.nz;
  b.set_levels(eta);
  b.upload_levels();
  b.p_top = G.p_top;
  const double mub = c.p_surf - G.p_top;
  // 1D base state (flat terrain): WRF module_initialize_ideal, dry.
  std::vector<double> pb(nz), tb(nz), alb(nz), phb(nz + 1), zm(nz), exn(nz);
  for (int k = 0; k < nz; ++k) {
    pb[k] = b.c3h[k] * mub + b.c4h[k] + G.p_top;
    tb[k] = snd.theta(snd.z_of_p(pb[k]));
    alb[k] = (r_d / p1000mb) * tb[k] * std::pow(pb[k] / p1000mb, cvpm);
    exn[k] = std::pow(pb[k] / p1000mb, rcp);
  }
  phb[0] = 0.0;  // flat terrain, ht = 0
  for (int k = 0; k < nz; ++k) phb[k + 1] = phb[k] - b.dnw[k] * (b.c1h[k] * mub + b.c2h[k]) * alb[k];
  for (int k = 0; k < nz; ++k) zm[k] = 0.5 * (phb[k] + phb[k + 1]) / g;
  b.mub.fill_host([&](int, int, int) { return mub; });
  b.pb.fill_host([&](int, int, int k) { return pb[k]; });
  b.alb.fill_host([&](int, int, int k) { return alb[k]; });
  b.t_init.fill_host([&](int, int, int k) { return tb[k] - t0; });
  b.phb.fill_host([&](int, int, int k) { return phb[k]; });

  // Perturbation state (dry): mu' = 0, p' = 0, theta = theta_b + theta'.
  const int nxh = G.nx, nyh = G.ny;
  std::vector<double> th((size_t)nxh * nyh * nz), al(th.size()), alt(th.size()), php((size_t)nxh * nyh * (nz + 1));
  auto I = [&](int i, int j, int k) { return (size_t(k) * nyh + j) * nxh + i; };
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < nyh; ++j)
      for (int i = 0; i < nxh; ++i) {
        double x = (i + 0.5) * G.dx, y = (j + 0.5) * G.dy;
        double t = tb[k] + c.theta_pert(x, y, zm[k], exn[k]);
        th[I(i, j, k)] = t;
        alt[I(i, j, k)] = (r_d / p1000mb) * t * std::pow(pb[k] / p1000mb, cvpm);
        al[I(i, j, k)] = alt[I(i, j, k)] - alb[k];
      }
  for (int j = 0; j < nyh; ++j)
    for (int i = 0; i < nxh; ++i) {
      double ph = 0.0;
      php[(size_t(0) * nyh + j) * nxh + i] = 0.0;
      for (int k = 0; k < nz; ++k) {
        ph -= b.dnw[k] * ((b.c1h[k] * mub + b.c2h[k]) * al[I(i, j, k)]);  // mu' = 0
        php[(size_t(k + 1) * nyh + j) * nxh + i] = ph;
      }
    }
  auto P = [&](int i, int j, int k) { return php[(size_t(k) * nyh + j) * nxh + i]; };
  for (dyn::State* ps : {&s}) {
    ps->t_2.fill_host([&](int i, int j, int k) { return th[I(i, j, k)] - t0; });
    ps->al.fill_host([&](int i, int j, int k) { return al[I(i, j, k)]; });
    ps->alt.fill_host([&](int i, int j, int k) { return alt[I(i, j, k)]; });
    ps->ph_2.fill_host(P);
    ps->u_2.fill_host([&](int i, int, int k) {
      if (G.bc_x == BC::Symmetric && (i == 0 || i == G.nx)) return 0.0;
      return c.u_env(zm[k]);
    });
    // p, mu_2, v_2, w_2, ww are zero (allocation zero-fills).
  }
  for (Field3* f : {&s.t_2, &s.al, &s.alt, &s.ph_2, &s.u_2, &s.v_2, &s.w_2, &s.p}) fill_halo(*f, G);
  fill_halo(s.mu_2, G);
  s.u_1.copy_from(s.u_2); s.v_1.copy_from(s.v_2); s.w_1.copy_from(s.w_2);
  s.t_1.copy_from(s.t_2); s.ph_1.copy_from(s.ph_2); s.mu_1.copy_from(s.mu_2);
}

}  // namespace cppwrf::cases
