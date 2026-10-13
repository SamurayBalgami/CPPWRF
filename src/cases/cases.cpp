#include "cases/cases.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <numbers>
#include <stdexcept>

#include "core/constants.hpp"
#include "core/halo.hpp"

namespace cppwrf::cases {

using namespace cppwrf::cst;

Sounding::Sounding(std::function<double(double)> theta, double p_surf, double zmax, double dz,
                   std::function<double(double)> qv)
    : th_(std::move(theta)), qv_(std::move(qv)), dz_(dz) {
  int n = int(std::ceil(zmax / dz)) + 2;
  pi_.resize(n);
  pi_[0] = std::pow(p_surf / p1000mb, rcp);
  // d(pi)/dz = -g / (cp theta_rho), integrated with Simpson's rule on each interval.
  for (int m = 1; m < n; ++m) {
    double z0 = (m - 1) * dz, z1 = m * dz;
    double f = (1.0 / theta_rho(z0) + 4.0 / theta_rho(0.5 * (z0 + z1)) + 1.0 / theta_rho(z1)) / 6.0;
    pi_[m] = pi_[m - 1] - g / cp * f * dz;
  }
  if (qv_) {  // vapour load, trapezoid from the top: rho_d = p / (R_d theta_m pi)
    auto rq = [&](int m) {
      const double z = m * dz, pi = pi_[m], pp = p1000mb * std::pow(pi, 1.0 / rcp);
      return pp / (r_d * th_(z) * (1.0 + rvovrd * qv_(z)) * pi) * qv_(z);
    };
    w_.assign(n, 0.0);
    for (int m = n - 2; m >= 0; --m) w_[m] = w_[m + 1] + 0.5 * g * dz * (rq(m) + rq(m + 1));
  }
}

double Sounding::theta_rho(double z) const {
  if (!qv_) return th_(z);
  const double q = qv_(z);
  return th_(z) * (1.0 + rvovrd * q) / (1.0 + q);
}

double Sounding::p(double z) const {
  double s = z / dz_;
  int m = int(std::floor(s));
  if (m < 0 || m + 1 >= int(pi_.size())) throw std::runtime_error("Sounding::p out of range");
  // integrate exactly from the tabulated point (Simpson on the partial interval)
  double z0 = m * dz_, h = z - z0;
  double f = (1.0 / theta_rho(z0) + 4.0 / theta_rho(z0 + 0.5 * h) + 1.0 / theta_rho(z)) / 6.0;
  double pi = pi_[m] - g / cp * f * h;
  return p1000mb * std::pow(pi, 1.0 / rcp);
}

double Sounding::vapour_load(double z) const {
  if (!qv_) return 0.0;
  double s = z / dz_;
  int m = int(std::floor(s));
  if (m < 0 || m + 1 >= int(w_.size())) throw std::runtime_error("Sounding::vapour_load out of range");
  double a = s - m;
  return (1.0 - a) * w_[m] + a * w_[m + 1];
}

double Sounding::pd(double z) const { return p(z) - vapour_load(z); }

double Sounding::z_of_p(double pp) const {
  double lo = 0, hi = (pi_.size() - 2) * dz_;
  if (pp > p(lo) || pp < p(hi)) throw std::runtime_error("Sounding::z_of_p out of range");
  for (int it = 0; it < 200 && hi - lo > 1e-9; ++it) {
    double mid = 0.5 * (lo + hi);
    (p(mid) > pp ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}

double Sounding::z_of_pd(double pp) const {
  if (!qv_) return z_of_p(pp);
  double lo = 0, hi = (pi_.size() - 2) * dz_;
  if (pp > pd(lo) || pp < pd(hi)) throw std::runtime_error("Sounding::z_of_pd out of range");
  for (int it = 0; it < 200 && hi - lo > 1e-9; ++it) {
    double mid = 0.5 * (lo + hi);
    (pd(mid) > pp ? lo : hi) = mid;
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
    const double u0 = cfg.get_double("case.u0", 0.0);  // uniform wind (inertial-oscillation test)
    c.u_env = [=](double) { return u0; };
    c.params = {{"theta0", th0}, {"n_bv", N}, {"u0", u0}};
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
  } else if (name == "squall2d") {
    // Weisman & Klemp (1982) analytic sounding (WRF em_quarter_ss / em_squall2d_x family) with
    // a warm bubble. theta: theta0 + (theta_tr - theta0) (z/z_tr)^(5/4) below z_tr, isothermal
    // T_tr above; RH = 1 - 0.75 (z/z_tr)^(5/4) (0.25 above z_tr); qv capped at qv_max.
    // u = u_shear tanh(z / 3000 m) - u_shift (default: no wind, so the storm stays symmetric).
    const double th0 = cfg.get_double("case.theta0", 300.0), thtr = cfg.get_double("case.theta_tr", 343.0);
    const double ztr = cfg.get_double("case.z_tr", 12000.0), ttr = cfg.get_double("case.t_tr", 213.0);
    const double qvmax = cfg.get_double("case.qv_max", 0.014);
    const double us = cfg.get_double("case.u_shear", 0.0), ushift = cfg.get_double("case.u_shift", 0.0);
    const double xc = cfg.get_double("case.xc", xcenter(G)), zc = cfg.get_double("case.zc", 1400.0);
    const double xr = cfg.get_double("case.xr", 10000.0), zr = cfg.get_double("case.zr", 1400.0);
    const double dth = cfg.get_double("case.dtheta", 2.0), psurf = cfg.get_double("case.p_surf", 1.0e5);
    auto theta = [=](double z) {
      return z <= ztr ? th0 + (thtr - th0) * std::pow(z / ztr, 1.25) : thtr * std::exp(g / (cp * ttr) * (z - ztr));
    };
    auto rh = [=](double z) { return z <= ztr ? 1.0 - 0.75 * std::pow(z / ztr, 1.25) : 0.25; };
    // qv depends on p(z), which depends on qv: fixed-point iteration on a 50 m table.
    const double zmax = G.ztop * 1.2 + 1000.0, dzt = 50.0;
    const int nt = int(zmax / dzt) + 2;
    auto qtab = std::make_shared<std::vector<double>>(nt, 0.0);
    auto qv_of = [qtab, dzt, nt](double z) {
      double f = std::clamp(z / dzt, 0.0, double(nt - 1) - 1e-9);
      int n = int(f);
      return (*qtab)[n] + (f - n) * ((*qtab)[n + 1] - (*qtab)[n]);
    };
    for (int it = 0; it < 8; ++it) {
      Sounding snd(theta, psurf, zmax, 0.5, it == 0 ? std::function<double(double)>{} : std::function<double(double)>(qv_of));
      std::vector<double> q(nt);
      for (int n = 0; n < nt; ++n) {
        const double z = std::min(n * dzt, zmax), pz = snd.p(z), T = theta(z) * std::pow(pz / p1000mb, rcp);
        const double es = 611.2 * std::exp(17.67 * (T - 273.15) / (T - 29.65));  // Pa (Kessler/WRF svp)
        q[n] = std::min(qvmax, rh(z) * ep_2 * es / (pz - es));
      }
      *qtab = q;
    }
    c.theta_env = theta;
    c.qv_env = qv_of;
    c.u_env = [=](double z) { return us * std::tanh(z / 3000.0) - ushift; };
    c.theta_pert = [=](double x, double, double z, double) {
      double L = std::sqrt(std::pow((x - xc) / xr, 2) + std::pow((z - zc) / zr, 2));
      return L <= 1.0 ? dth * std::pow(std::cos(0.5 * pi * L), 2) : 0.0;
    };
    c.params = {{"theta0", th0}, {"theta_tr", thtr}, {"z_tr", ztr}, {"t_tr", ttr}, {"qv_max", qvmax},
                {"u_shear", us}, {"u_shift", ushift}, {"xc", xc}, {"zc", zc}, {"xr", xr}, {"zr", zr},
                {"dtheta", dth}};
  } else {
    throw std::runtime_error("unknown case: " + name);
  }
  c.p_surf = cfg.get_double("case.p_surf", 1.0e5);
  // optional unsaturated vapour profile qv(z) = qv0 exp(-z / qv_hscale) for any case
  if (cfg.has("case.qv0")) {
    const double qv0 = cfg.get_double("case.qv0"), hq = cfg.get_double("case.qv_hscale", 2500.0);
    c.qv_env = [=](double z) { return qv0 * std::exp(-z / hq); };
    c.params["qv0"] = qv0;
    c.params["qv_hscale"] = hq;
  }
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

static void initialize_moist(const CaseDef& c, const Sounding& snd, const Grid& G, dyn::BaseState& b, dyn::State& s,
                             const std::vector<double>& pb, const std::vector<double>& alb,
                             const std::vector<double>& zm, const std::vector<double>& exn, double mub);

void initialize(const CaseDef& c, const Sounding& snd, const std::vector<double>& eta, const Grid& G,
                dyn::BaseState& b, dyn::State& s, const Sounding* moist) {
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
  if (moist && s.P_QV >= 0) initialize_moist(c, *moist, G, b, s, pb, alb, zm, exn, mub);
}

// Full moist state, WRF module_initialize_ideal (lines after "getting moist sounding for full
// state") with use_theta_m = 1. Fortran k -> C++ k-1; w level kk+1 -> C++ k+1.
static void initialize_moist(const CaseDef& c, const Sounding& snd, const Grid& G, dyn::BaseState& b, dyn::State& s,
                             const std::vector<double>& pb, const std::vector<double>& alb,
                             const std::vector<double>& zm, const std::vector<double>& exn, double mub) {
  const int nz = G.nz, nx = G.nx, ny = G.ny;
  const double pd_surf = c.p_surf - snd.vapour_load(0.0);
  const double mu1 = pd_surf - G.p_top - mub;  // dry-mass perturbation from the vapour load
  std::vector<double> qv(nz), thd(nz);
  for (int k = 0; k < nz; ++k) {
    const double p_level = b.c3h[k] * (pd_surf - G.p_top) + b.c4h[k] + G.p_top;
    const double z = snd.z_of_pd(p_level);
    qv[k] = snd.qv(z);
    thd[k] = snd.theta(z);
  }
  const size_t ncol = size_t(nx) * ny;
  auto I = [&](int i, int j, int k) { return (size_t(k) * ny + j) * nx + i; };
  std::vector<double> tm(ncol * nz), al(ncol * nz), alt(ncol * nz), pp(ncol * nz), php(ncol * (nz + 1));
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const double x = (i + 0.5) * G.dx, y = (j + 0.5) * G.dy;
      std::vector<double> th(nz);
      for (int k = 0; k < nz; ++k) th[k] = thd[k] + c.theta_pert(x, y, zm[k], exn[k]);
      auto eos = [&](int k, double p) {  // alt = R_d/p0 theta (1 + rvovrd qv) ((p + pb)/p0)^cvpm
        const double qvf = 1.0 + rvovrd * qv[k];
        alt[I(i, j, k)] = (r_d / p1000mb) * th[k] * qvf * std::pow((p + pb[k]) / p1000mb, cvpm);
        al[I(i, j, k)] = alt[I(i, j, k)] - alb[k];
        pp[I(i, j, k)] = p;
      };
      // hydrostatic p' from the top (pg_buoy_w balance with moisture loading)
      int k = nz - 1;
      double qvf1 = 0.5 * (qv[k] + qv[k]);
      double qvf2 = 1.0 / (1.0 + qvf1);
      qvf1 = qvf1 * qvf2;
      double p = -0.5 * ((b.c1f[nz] * mu1) + qvf1 * (b.c1f[nz] * mub + b.c2f[nz])) / b.rdnw[k] / qvf2;
      eos(k, p);
      for (k = nz - 2; k >= 0; --k) {
        qvf1 = 0.5 * (qv[k] + qv[k + 1]);
        qvf2 = 1.0 / (1.0 + qvf1);
        qvf1 = qvf1 * qvf2;
        p = p - ((b.c1f[k + 1] * mu1) + qvf1 * (b.c1f[k + 1] * mub + b.c2f[k + 1])) / qvf2 / b.rdn[k + 1];
        eos(k, p);
      }
      double ph = 0.0;
      php[I(i, j, 0)] = 0.0;
      for (k = 0; k < nz; ++k) {
        ph = ph - b.dnw[k] * (((b.c1h[k] * mub + b.c2h[k]) + b.c1h[k] * mu1) * al[I(i, j, k)] + b.c1h[k] * mu1 * alb[k]);
        php[(size_t(k + 1) * ny + j) * nx + i] = ph;
      }
      for (k = 0; k < nz; ++k) tm[I(i, j, k)] = th[k] * (1.0 + rvovrd * qv[k]) - t0;  // theta_m - t0
    }
  s.t_2.fill_host([&](int i, int j, int k) { return tm[I(i, j, k)]; });
  s.al.fill_host([&](int i, int j, int k) { return al[I(i, j, k)]; });
  s.alt.fill_host([&](int i, int j, int k) { return alt[I(i, j, k)]; });
  s.p.fill_host([&](int i, int j, int k) { return pp[I(i, j, k)]; });
  s.ph_2.fill_host([&](int i, int j, int k) { return php[(size_t(k) * ny + j) * nx + i]; });
  s.mu_2.fill_host([&](int, int, int) { return mu1; });
  s.q_2[s.P_QV].fill_host([&](int, int, int k) { return qv[k]; });
  b.qv_base.fill_host([&](int, int, int k) { return qv[k]; });
  for (Field3* f : {&s.t_2, &s.al, &s.alt, &s.ph_2, &s.p, &s.q_2[s.P_QV]}) fill_halo(*f, G);
  fill_halo(s.mu_2, G);
  s.t_1.copy_from(s.t_2); s.ph_1.copy_from(s.ph_2); s.mu_1.copy_from(s.mu_2);
  s.q_1[s.P_QV].copy_from(s.q_2[s.P_QV]);
}

}  // namespace cppwrf::cases
