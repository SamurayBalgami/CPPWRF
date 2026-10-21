#include "cases/wrf_ideal.hpp"
#include "cases/model.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "core/constants.hpp"
#include "core/halo.hpp"
#include "dyn/rk3.hpp"

namespace cppwrf::cases {

namespace {

using cst::cvpm;
using cst::g;
using cst::p1000mb;
using cst::r_d;
using cst::rcp;
using cst::rvovrd;
using cst::t0;

enum class Ideal { hill2d_x, quarter_ss, squall2d_x, squall2d_y, grav2d_x };

Ideal ideal_of(const std::string& name) {
  if (name == "em_hill2d_x") return Ideal::hill2d_x;
  if (name == "em_quarter_ss") return Ideal::quarter_ss;
  if (name == "em_squall2d_x") return Ideal::squall2d_x;
  if (name == "em_squall2d_y") return Ideal::squall2d_y;
  if (name == "em_grav2d_x") return Ideal::grav2d_x;
  throw std::runtime_error("unknown WRF ideal case " + name);
}

// 1-based Fortran-style vector: v(1..n)
struct F1 {
  std::vector<double> d;
  explicit F1(int n = 0) : d(std::size_t(n) + 1, 0.0) {}
  double& operator()(int i) { return d[std::size_t(i)]; }
  double operator()(int i) const { return d[std::size_t(i)]; }
};

// module_init_utilities interp_0
double interp_0(const F1& v_in, const F1& z_in, double z_out, int nz_in) {
  double height = z_out, w1, w2, r = 0.0;
  if (z_in(nz_in) > z_in(1)) {
    if (height > z_in(nz_in)) {
      w2 = (z_in(nz_in) - height) / (z_in(nz_in) - z_in(nz_in - 1));
      w1 = 1. - w2;
      r = w1 * v_in(nz_in) + w2 * v_in(nz_in - 1);
    } else if (height < z_in(1)) {
      w2 = (z_in(2) - height) / (z_in(2) - z_in(1));
      w1 = 1. - w2;
      r = w1 * v_in(2) + w2 * v_in(1);
    } else {
      bool interp = false;
      int kp = nz_in;
      while (!interp && kp >= 2) {
        if (z_in(kp) >= height && z_in(kp - 1) <= height) {
          w2 = (height - z_in(kp)) / (z_in(kp - 1) - z_in(kp));
          w1 = 1. - w2;
          r = w1 * v_in(kp) + w2 * v_in(kp - 1);
          interp = true;
        }
        kp = kp - 1;
      }
    }
  } else {
    if (height < z_in(nz_in)) {
      w2 = (z_in(nz_in) - height) / (z_in(nz_in) - z_in(nz_in - 1));
      w1 = 1. - w2;
      r = w1 * v_in(nz_in) + w2 * v_in(nz_in - 1);
    } else if (height > z_in(1)) {
      w2 = (z_in(2) - height) / (z_in(2) - z_in(1));
      w1 = 1. - w2;
      r = w1 * v_in(2) + w2 * v_in(1);
    } else {
      bool interp = false;
      int kp = nz_in;
      while (!interp && kp >= 2) {
        if (z_in(kp) <= height && z_in(kp - 1) >= height) {
          w2 = (height - z_in(kp)) / (z_in(kp - 1) - z_in(kp));
          w1 = 1. - w2;
          r = w1 * v_in(kp) + w2 * v_in(kp - 1);
          interp = true;
        }
        kp = kp - 1;
      }
    }
  }
  return r;
}

struct WrfSounding {
  int nl = 0;
  double th_surf = 0.0;
  F1 zk, p, p_dry, theta, rho, u, v, qv;
};

// read_sounding + get_sounding (module_initialize_ideal.F)
WrfSounding get_sounding(const std::string& path, bool dry) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("cannot open input_sounding " + path);
  std::vector<double> vals;
  std::string line;
  double p_surf = 0, th_surf = 0, qv_surf = 0;
  if (!std::getline(f, line)) throw std::runtime_error("empty input_sounding " + path);
  {
    std::istringstream s(line);
    if (!(s >> p_surf >> th_surf >> qv_surf)) throw std::runtime_error("bad input_sounding header " + path);
  }
  std::vector<double> h, th, qvi, ui, vi;
  while (std::getline(f, line)) {
    std::istringstream s(line);
    double a, b, c, d, e;
    if (!(s >> a >> b >> c >> d >> e)) break;
    h.push_back(a); th.push_back(b); qvi.push_back(c); ui.push_back(d); vi.push_back(e);
  }
  const int nl = int(h.size());
  if (nl < 2) throw std::runtime_error("input_sounding needs >= 2 levels: " + path);
  F1 h_input(nl), th_input(nl), qv_input(nl), u_input(nl), v_input(nl);
  for (int k = 1; k <= nl; ++k) {
    h_input(k) = h[k - 1]; th_input(k) = th[k - 1]; qv_input(k) = qvi[k - 1];
    u_input(k) = ui[k - 1]; v_input(k) = vi[k - 1];
  }
  if (dry)
    for (int k = 1; k <= nl; ++k) qv_input(k) = 0.;
  for (int k = 1; k <= nl; ++k) qv_input(k) = 0.001 * qv_input(k);
  p_surf = 100. * p_surf;
  const double r = r_d;
  double qvf = 1. + rvovrd * qv_input(1);
  const double rho_surf = 1. / ((r / p1000mb) * th_surf * qvf * std::pow(p_surf / p1000mb, cvpm));
  F1 rho_input(nl), pm_input(nl), p_input(nl);
  qvf = 1. + rvovrd * qv_input(1);
  double qvf1 = 1. + qv_input(1);
  rho_input(1) = rho_surf;
  double dz = h_input(1);
  for (int it = 1; it <= 10; ++it) {
    pm_input(1) = p_surf - 0.5 * dz * (rho_surf + rho_input(1)) * g * qvf1;
    rho_input(1) = 1. / ((r / p1000mb) * th_input(1) * qvf * std::pow(pm_input(1) / p1000mb, cvpm));
  }
  for (int k = 2; k <= nl; ++k) {
    rho_input(k) = rho_input(k - 1);
    dz = h_input(k) - h_input(k - 1);
    qvf1 = 0.5 * (2. + (qv_input(k - 1) + qv_input(k)));
    qvf = 1. + rvovrd * qv_input(k);
    for (int it = 1; it <= 10; ++it) {
      pm_input(k) = pm_input(k - 1) - 0.5 * dz * (rho_input(k) + rho_input(k - 1)) * g * qvf1;
      if (pm_input(k) <= 0.) throw std::runtime_error("Integrated pressure has gone negative - too cold for chosen height");
      rho_input(k) = 1. / ((r / p1000mb) * th_input(k) * qvf * std::pow(pm_input(k) / p1000mb, cvpm));
    }
  }
  p_input(nl) = pm_input(nl);
  for (int k = nl - 1; k >= 1; --k) {
    dz = h_input(k + 1) - h_input(k);
    p_input(k) = p_input(k + 1) + 0.5 * dz * (rho_input(k) + rho_input(k + 1)) * g;
  }
  WrfSounding s;
  s.nl = nl;
  s.th_surf = th_surf;
  s.zk = h_input; s.p = pm_input; s.p_dry = p_input; s.theta = th_input; s.rho = rho_input;
  s.u = u_input; s.v = v_input; s.qv = qv_input;
  return s;
}

}  // namespace

bool is_wrf_ideal_case(const std::string& name) {
  return name == "em_hill2d_x" || name == "em_quarter_ss" || name == "em_squall2d_x" || name == "em_squall2d_y" ||
         name == "em_grav2d_x";
}

void init_wrf_ideal(const io::Config& cfg, Model& m, const std::string& name) {
  const Ideal ic = ideal_of(name);
  Grid& gr = m.g;
  gr.nx = int(cfg.get_int("domain.nx"));
  gr.ny = int(cfg.get_int("domain.ny", 1));
  gr.nz = int(cfg.get_int("domain.nz"));
  gr.dx = cfg.get_double("domain.dx");
  gr.dy = cfg.get_double("domain.dy", gr.dx);
  gr.ztop = cfg.get_double("domain.ztop");
  gr.bc_x = parse_bc(cfg.get_string("domain.bc_x", "periodic"));
  gr.bc_y = parse_bc(cfg.get_string("domain.bc_y", "periodic"));
  gr.dt = cfg.get_double("time.dt");
  m.run_time = cfg.get_double("run.run_time", 0.0);
  m.nsteps = std::lround(m.run_time / gr.dt);
  if (std::fabs(m.nsteps * gr.dt - m.run_time) > 1e-9 * std::max(1.0, m.run_time))
    throw std::runtime_error("run.run_time must be a multiple of time.dt");
  const int hybrid_opt = int(cfg.get_int("dyn.hybrid_opt", 2));
  const double etac = cfg.get_double("dyn.etac", 0.2);
  m.hybrid_opt = hybrid_opt;
  // share/input_wrf.F: "Ideal cases always use hypsometric_opt=1, regardless of namelist setting"
  const int hypsometric_opt = 1;
  const std::string sounding = cfg.get_string("case.sounding_file", "input_sounding");

  // WRF domain indices (Fortran, 1-based): ide = e_we, jde = e_sn, kde = e_vert; 2D x-z: e_sn = 3
  const int nx = gr.nx, ny = gr.ny, nz = gr.nz;
  const int ids = 1, ide = nx + 1, jds = 1, jde = (ny == 1 ? 3 : ny + 1), kde = nz + 1, kte = kde;
  const double dx = gr.dx, dy = gr.dy;

  // ---- ideal_constants
  double delt = 0., z_scale = 8000. / gr.ztop, hm = 0., xa = 0.;
  int icm = 0, nxc = (ide - ids) / 2, nyc = (jde - jds) / 2;
  const bool stretch_grid = true;
  const double pi = 2. * std::asin(1.0);
  switch (ic) {
    case Ideal::hill2d_x: hm = 100.; xa = 5.0; icm = ide / 2; delt = 0.; break;
    case Ideal::quarter_ss: case Ideal::squall2d_x: case Ideal::squall2d_y: delt = 3.; break;
    case Ideal::grav2d_x: break;
  }
  hm = cfg.get_double("case.hm", hm);  // CPPWRF options (WRF hard-codes them)
  xa = cfg.get_double("case.xa", xa);

  // ---- ideal_levels (stretched: the k loop is vectorised in WRF, libmvec exp)
  F1 znw(kde), dnw(kde), rdnw(kde), znu(kde), dn(kde), rdn(kde), fnp(kde), fnm(kde);
  for (int k = 1; k <= kde; ++k) {
    if (stretch_grid)
      znw(k) = (vm::exp(-double(k - 1) / double(kde - 1) / z_scale, k - 1, kde) - std::exp(-1. / z_scale)) /
               (1. - std::exp(-1. / z_scale));
    else
      znw(k) = 1. - double(k - 1) / double(kde - 1);
  }
  for (int k = 1; k <= kde - 1; ++k) {
    dnw(k) = znw(k + 1) - znw(k);
    rdnw(k) = 1. / dnw(k);
    znu(k) = 0.5 * (znw(k + 1) + znw(k));
  }
  for (int k = 2; k <= kde - 1; ++k) {
    dn(k) = 0.5 * (dnw(k) + dnw(k - 1));
    rdn(k) = 1. / dn(k);
    fnp(k) = .5 * dnw(k) / dn(k);
    fnm(k) = .5 * dnw(k - 1) / dn(k);
  }
  const double cof1 = (2. * dn(2) + dn(3)) / (dn(2) + dn(3)) * dnw(1) / dn(2);
  const double cof2 = dn(2) / (dn(2) + dn(3)) * dnw(1) / dn(3);
  const double cf1 = fnp(2) + cof1, cf2 = fnm(2) - cof1 - cof2, cf3 = cof2;
  const double cfn = (.5 * dnw(kde - 1) + dn(kde - 1)) / dn(kde - 1);
  const double cfn1 = -.5 * dnw(kde - 1) / dn(kde - 1);
  const double rdx = 1. / dx, rdy = 1. / dy;

  // ---- dry sounding for the base state, p_top
  WrfSounding snd = get_sounding(sounding, true);
  int nl_in = snd.nl;
  const double p_top = interp_0(snd.p, snd.zk, gr.ztop, nl_in);
  gr.p_top = p_top;

  // ---- hybrid coordinate coefficients
  F1 c1f(kde), c2f(kde), c3f(kde), c4f(kde), c1h(kde), c2h(kde), c3h(kde), c4h(kde);
  for (int k = 1; k <= kde; ++k) {
    if (hybrid_opt == 0 || hybrid_opt == 1) {
      c3f(k) = znw(k);
    } else if (hybrid_opt == 2) {
      const double B1 = 2. * powi(etac, 2) * (1. - etac);
      const double B2 = -etac * (4. - 3. * etac - powi(etac, 3));
      const double B3 = 2. * (1. - powi(etac, 3));
      const double B4 = -(1. - powi(etac, 2));
      const double B5 = powi(1. - etac, 4);
      c3f(k) = (B1 + B2 * znw(k) + B3 * powi(znw(k), 2) + B4 * powi(znw(k), 3)) / B5;
      if (znw(k) < etac) c3f(k) = 0.;
      if (k == 1) c3f(k) = 1.;
      else if (k == kde) c3f(k) = 0.;
    } else {
      throw std::runtime_error("dyn.hybrid_opt must be 0, 1 or 2 for the WRF ideal cases");
    }
  }
  for (int k = 1; k <= kde; ++k) c4f(k) = (znw(k) - c3f(k)) * (p1000mb - p_top);
  for (int k = 1; k <= kde - 1; ++k) {
    c3h(k) = (c3f(k + 1) + c3f(k)) * 0.5;
    c4h(k) = (znu(k) - c3h(k)) * (p1000mb - p_top);
  }
  for (int k = 2; k <= kde - 1; ++k) c1f(k) = (c3h(k) - c3h(k - 1)) / (znu(k) - znu(k - 1));
  c1f(1) = 1.;
  c1f(kde) = (hybrid_opt == 0 || hybrid_opt == 1) ? 1. : 0.;
  for (int k = 1; k <= kde; ++k) c2f(k) = (1. - c1f(k)) * (p1000mb - p_top);
  for (int k = 1; k <= kde - 1; ++k) {
    c1h(k) = (c3f(k + 1) - c3f(k)) / (znw(k + 1) - znw(k));
    c2h(k) = (1. - c1h(k)) * (p1000mb - p_top);
  }

  // ---- 3D work arrays on the WRF mass/w columns: (i, k, j) Fortran, i = 1..ide-1, j = 1..jm
  const int jm = (ny == 1 ? 1 : ny);  // x-z: WRF's two identical rows, one kept
  const int im = nx;
  auto c2 = [&](int i, int j) { return std::size_t(j - 1) * im + (i - 1); };
  auto c3 = [&](int i, int k, int j) { return (std::size_t(k - 1) * jm + (j - 1)) * im + (i - 1); };
  std::vector<double> ht(std::size_t(im) * jm), mub(ht.size()), mu_1(ht.size());
  const std::size_t n3 = std::size_t(im) * jm * kde;
  std::vector<double> phb(n3), pb(n3), t_init(n3), alb(n3), ph_1(n3), t_1(n3), qv(n3), p(n3), alt(n3), al(n3);

  // ---- ideal_terrain
  for (int j = 1; j <= jm; ++j)
    for (int i = 1; i <= im; ++i) {
      if (ic == Ideal::hill2d_x) {
        const double r = double(i - icm) / xa;
        ht[c2(i, j)] = hm / (1. + r * r);
      } else {
        ht[c2(i, j)] = 0.;
      }
      phb[c3(i, 1, j)] = g * ht[c2(i, j)];
    }

  // ---- base state
  for (int j = 1; j <= jm; ++j)
    for (int i = 1; i <= im; ++i) {
      const double p_surf = interp_0(snd.p, snd.zk, phb[c3(i, 1, j)] / g, nl_in);
      mub[c2(i, j)] = p_surf - p_top;
      for (int k = 1; k <= kte - 1; ++k) {
        const double p_level = c3h(k) * (p_surf - p_top) + c4h(k) + p_top;
        pb[c3(i, k, j)] = p_level;
        t_init[c3(i, k, j)] = interp_0(snd.theta, snd.p, p_level, nl_in) - t0;
        alb[c3(i, k, j)] = (r_d / p1000mb) * (t_init[c3(i, k, j)] + t0) * std::pow(pb[c3(i, k, j)] / p1000mb, cvpm);
      }
      for (int kk = 2; kk <= kte; ++kk) {
        const int k = kk - 1;
        phb[c3(i, kk, j)] = phb[c3(i, kk - 1, j)] -
                            dnw(kk - 1) * (c1h(k) * mub[c2(i, j)] + c2h(k)) * alb[c3(i, kk - 1, j)];
      }
    }

  // ---- moist sounding for the full state
  snd = get_sounding(sounding, false);
  nl_in = snd.nl;
  for (int j = 1; j <= jm; ++j)
    for (int i = 1; i <= im; ++i) {
      const double pd_surf = interp_0(snd.p_dry, snd.zk, phb[c3(i, 1, j)] / g, nl_in);
      mu_1[c2(i, j)] = pd_surf - p_top - mub[c2(i, j)];
      for (int k = 1; k <= kde - 1; ++k) {
        const double p_level = c3h(k) * (pd_surf - p_top) + c4h(k) + p_top;
        qv[c3(i, k, j)] = interp_0(snd.qv, snd.p_dry, p_level, nl_in);
        t_1[c3(i, k, j)] = interp_0(snd.theta, snd.p_dry, p_level, nl_in) - t0;
      }
      const double mu1 = mu_1[c2(i, j)], mb = mub[c2(i, j)];
      int kk = kte - 1;
      int k = kk + 1;
      double qvf1 = 0.5 * (qv[c3(i, kk, j)] + qv[c3(i, kk, j)]);
      double qvf2 = 1. / (1. + qvf1);
      qvf1 = qvf1 * qvf2;
      p[c3(i, kk, j)] = -0.5 * ((c1f(k) * mu1) + qvf1 * (c1f(k) * mb + c2f(k))) / rdnw(kk) / qvf2;
      double qvf = 1. + rvovrd * qv[c3(i, kk, j)];
      alt[c3(i, kk, j)] = (r_d / p1000mb) * (t_1[c3(i, kk, j)] + t0) * qvf *
                          std::pow((p[c3(i, kk, j)] + pb[c3(i, kk, j)]) / p1000mb, cvpm);
      al[c3(i, kk, j)] = alt[c3(i, kk, j)] - alb[c3(i, kk, j)];
      for (kk = kte - 2; kk >= 1; --kk) {
        k = kk + 1;
        qvf1 = 0.5 * (qv[c3(i, kk, j)] + qv[c3(i, kk + 1, j)]);
        qvf2 = 1. / (1. + qvf1);
        qvf1 = qvf1 * qvf2;
        p[c3(i, kk, j)] = p[c3(i, kk + 1, j)] - ((c1f(k) * mu1) + qvf1 * (c1f(k) * mb + c2f(k))) / qvf2 / rdn(kk + 1);
        qvf = 1. + rvovrd * qv[c3(i, kk, j)];
        alt[c3(i, kk, j)] = (r_d / p1000mb) * (t_1[c3(i, kk, j)] + t0) * qvf *
                            std::pow((p[c3(i, kk, j)] + pb[c3(i, kk, j)]) / p1000mb, cvpm);
        al[c3(i, kk, j)] = alt[c3(i, kk, j)] - alb[c3(i, kk, j)];
      }
      ph_1[c3(i, 1, j)] = 0.;
      for (kk = 2; kk <= kte; ++kk) {
        k = kk - 1;
        ph_1[c3(i, kk, j)] = ph_1[c3(i, kk - 1, j)] -
                             (dnw(kk - 1)) * (((c1h(k) * mb + c2h(k)) + (c1h(k) * mu1)) * al[c3(i, kk - 1, j)] +
                                              (c1h(k) * mu1) * alb[c3(i, kk - 1, j)]);
      }
    }

  // ---- ideal_pert
  auto rebalance_column = [&](int i, int j) {
    const double mu1 = mu_1[c2(i, j)], mb = mub[c2(i, j)];
    for (int k = 2; k <= kte; ++k)
      ph_1[c3(i, k, j)] = ph_1[c3(i, k - 1, j)] -
                          (dnw(k - 1)) * (((c1h(k - 1) * mb + c2h(k - 1)) + (c1h(k - 1) * mu1)) * al[c3(i, k - 1, j)] +
                                          (c1h(k - 1) * mu1) * alb[c3(i, k - 1, j)]);
  };
  auto set_alt = [&](int i, int k, int j) {
    const double qvf = 1. + rvovrd * qv[c3(i, k, j)];
    alt[c3(i, k, j)] = (r_d / p1000mb) * (t_1[c3(i, k, j)] + t0) * qvf *
                       std::pow((p[c3(i, k, j)] + pb[c3(i, k, j)]) / p1000mb, cvpm);
    al[c3(i, k, j)] = alt[c3(i, k, j)] - alb[c3(i, k, j)];
  };
  auto zrad_of = [&](int i, int k, int j) {
    return 0.5 * (ph_1[c3(i, k, j)] + ph_1[c3(i, k + 1, j)] + phb[c3(i, k, j)] + phb[c3(i, k + 1, j)]) / g;
  };
  if (ic == Ideal::quarter_ss || ic == Ideal::squall2d_x || ic == Ideal::squall2d_y) {
    for (int j = 1; j <= jm; ++j) {
      const double yrad = ic == Ideal::quarter_ss   ? dy * double(j - nyc) / 10000.
                          : ic == Ideal::squall2d_y ? dy * double(j - nyc) / 4000.
                                                    : 0.;
      for (int i = 1; i <= im; ++i) {
        const double xrad = ic == Ideal::quarter_ss   ? dx * double(i - nxc) / 10000.
                            : ic == Ideal::squall2d_x ? dx * double(i - nxc) / 4000.
                                                      : 0.;
        for (int k = 1; k <= kte - 1; ++k) {
          double zrad = zrad_of(i, k, j);
          zrad = (zrad - 1500.) / 1500.;
          const double rad = std::sqrt(xrad * xrad + yrad * yrad + zrad * zrad);
          if (rad <= 1.) {
            const double c = std::cos(.5 * pi * rad);
            t_1[c3(i, k, j)] = t_1[c3(i, k, j)] + delt * (c * c);
            set_alt(i, k, j);
          }
        }
        rebalance_column(i, j);
      }
    }
  } else if (ic == Ideal::grav2d_x) {
    const double u_mean = 00.;
    const double xpos = dx * nxc - u_mean * 900.;
    const double xposml = xpos - dx * (ide - 1);
    const double xpospl = xpos + dx * (ide - 1);
    for (int j = 1; j <= jm; ++j)
      for (int i = 1; i <= im; ++i) {
        const double xrad = std::min({std::fabs(dx * double(i) - xpos), std::fabs(dx * double(i) - xposml),
                                      std::fabs(dx * double(i) - xpospl)}) / 4000.;
        for (int k = 1; k <= kte - 1; ++k) {
          double zrad = zrad_of(i, k, j);
          zrad = (zrad - 3000.) / 2000.;
          const double rad = std::sqrt(xrad * xrad + zrad * zrad);
          if (rad <= 1.) {
            const double dlt = -15.0 / std::pow((p[c3(i, k, j)] + pb[c3(i, k, j)]) / p1000mb, rcp);
            t_1[c3(i, k, j)] = t_1[c3(i, k, j)] + dlt * (std::cos(pi * rad) + 1.0) / 2.;
            set_alt(i, k, j);
          }
        }
        rebalance_column(i, j);
      }
  }

  // ---- winds: v at (i, j) v points from z_at_v, u at u points from z_at_u (sounding u, v)
  const int nxs = nx + 1, jvm = (ny == 1 ? 1 : ny + 1);
  std::vector<double> u1(std::size_t(nxs) * jm * nz), v1(std::size_t(nx) * jvm * nz);
  for (int j = 1; j <= jvm; ++j)
    for (int i = 1; i <= im; ++i) {
      double z_at_v;
      if (ny == 1) z_at_v = phb[c3(i, 1, 1)] / g;  // x-z: both rows identical
      else if (j == jds) z_at_v = phb[c3(i, 1, j)] / g;
      else if (j == ny + 1) z_at_v = phb[c3(i, 1, j - 1)] / g;
      else z_at_v = 0.5 * (phb[c3(i, 1, j)] + phb[c3(i, 1, j - 1)]) / g;
      const double p_surf = interp_0(snd.p, snd.zk, z_at_v, nl_in);
      for (int k = 1; k <= nz; ++k) {
        const double p_level = c3h(k) * (p_surf - p_top) + c4h(k) + p_top;
        v1[(std::size_t(k - 1) * jvm + (j - 1)) * nx + (i - 1)] = interp_0(snd.v, snd.p, p_level, nl_in);
      }
    }
  // x-z: the dycore does not advance v (WRF does, on its two periodic rows), so only v = 0 is WRF's run
  if (ny == 1)
    for (double v : v1)
      if (v != 0.) throw std::runtime_error("x-z ideal case: the input_sounding v column must be 0 (v is not advanced in 2-D)");
  for (int j = 1; j <= jm; ++j)
    for (int i = 1; i <= nxs; ++i) {
      double z_at_u;
      if (i == ids) z_at_u = phb[c3(i, 1, j)] / g;
      else if (i == ide) z_at_u = phb[c3(i - 1, 1, j)] / g;
      else z_at_u = 0.5 * (phb[c3(i, 1, j)] + phb[c3(i - 1, 1, j)]) / g;
      const double p_surf = interp_0(snd.p, snd.zk, z_at_u, nl_in);
      for (int k = 1; k <= nz; ++k) {
        const double p_level = c3h(k) * (p_surf - p_top) + c4h(k) + p_top;
        u1[(std::size_t(k - 1) * jm + (j - 1)) * nxs + (i - 1)] = interp_0(snd.u, snd.p, p_level, nl_in);
      }
    }

  // ---- 1D base profiles from column (1, 1)
  m.base = std::make_unique<dyn::BaseState>(gr, false);
  dyn::BaseState& b = *m.base;
  for (int k = 1; k <= kte - 1; ++k) {
    b.t_base[k - 1] = t_1[c3(1, k, 1)];
    b.qv_base_1d[k - 1] = qv[c3(1, k, 1)];
    b.u_base[k - 1] = u1[(std::size_t(k - 1) * jm) * nxs];
    b.v_base[k - 1] = v1[(std::size_t(k - 1) * jvm) * nx];
    b.z_base[k - 1] = 0.5 * (phb[c3(1, k, 1)] + phb[c3(1, k + 1, 1)] + ph_1[c3(1, k, 1)] + ph_1[c3(1, k + 1, 1)]) / g;
  }

  // ---- ideal_surfacet (module_initialize_ideal.F:1593): TSK from the lowest three levels (dry theta, p of the
  // initialisation), TMN = TSK - 0.5; and the dry theta of the first history frame
  {
    std::vector<double> tsk(std::size_t(im) * jm), tmn(tsk.size()), th(std::size_t(nz) * jm * im);
    for (int j = 1; j <= jm; ++j)
      for (int i = 1; i <= im; ++i) {
        double temp[3];
        for (int k = 1; k <= 3; ++k)
          temp[k - 1] = (t_1[c3(i, k, j)] + t0) * std::pow((p[c3(i, k, j)] + pb[c3(i, k, j)]) / p1000mb, rcp);
        const std::size_t ij = std::size_t(j - 1) * im + (i - 1);
        tsk[ij] = cf1 * temp[0] + cf2 * temp[1] + cf3 * temp[2];
        tmn[ij] = tsk[ij] - 0.5;
        for (int k = 1; k <= nz; ++k) th[(std::size_t(k - 1) * jm + (j - 1)) * im + (i - 1)] = t_1[c3(i, k, j)];
      }
    m.init_out["TSK"] = tsk; m.init_out["TMN"] = tmn; m.init_out["T"] = th;
  }

  // ---- use_theta_m = 1: t_2 = theta_m - t0 (t_1 keeps theta until start_em sets t_1 = t_2)
  std::vector<double> t_2(t_1);
  if (cfg.get_int("dyn.use_theta_m", 1) == 1)
    for (int j = 1; j <= jm; ++j)
      for (int k = 1; k <= kte - 1; ++k)
        for (int i = 1; i <= im; ++i)
          t_2[c3(i, k, j)] = (t_2[c3(i, k, j)] + t0) * (1. + (cst::r_v / r_d) * qv[c3(i, k, j)]) - t0;

  // ---- start_em, first trip of an ideal run (max_dom = 1): pb, alb from MUB and t_init; al from
  // ph_1 (hypsometric_opt), p from theta_m (qvf = 1)
  for (int j = 1; j <= jm; ++j)
    for (int k = 1; k <= kte - 1; ++k)
      for (int i = 1; i <= im; ++i) {
        pb[c3(i, k, j)] = c3h(k) * mub[c2(i, j)] + c4h(k) + p_top;
        alb[c3(i, k, j)] = (r_d / p1000mb) * (t_init[c3(i, k, j)] + t0) * std::pow(pb[c3(i, k, j)] / p1000mb, cvpm);
      }
  for (int j = 1; j <= jm; ++j)
    for (int k = 1; k <= kte - 1; ++k)
      for (int i = 1; i <= im; ++i) {
        const double mu = mu_1[c2(i, j)], mb = mub[c2(i, j)];
        if (hypsometric_opt == 1) {
          al[c3(i, k, j)] = -1. / ((c1h(k) * mb + c2h(k)) + (c1h(k) * mu)) *
                            (alb[c3(i, k, j)] * (c1h(k) * mu) + rdnw(k) * (ph_1[c3(i, k + 1, j)] - ph_1[c3(i, k, j)]));
        } else {
          const double pfu = c3f(k + 1) * (mb + mu) + c4f(k + 1) + p_top;
          const double pfd = c3f(k) * (mb + mu) + c4f(k) + p_top;
          const double phm = c3h(k) * (mb + mu) + c4h(k) + p_top;
          al[c3(i, k, j)] = (ph_1[c3(i, k + 1, j)] - ph_1[c3(i, k, j)] + phb[c3(i, k + 1, j)] - phb[c3(i, k, j)]) / phm /
                                std::log(pfd / pfu) -
                            alb[c3(i, k, j)];
        }
        const double qvf = 1.;
        p[c3(i, k, j)] = p1000mb * std::pow((r_d * (t0 + t_2[c3(i, k, j)]) * qvf) /
                                                (p1000mb * (al[c3(i, k, j)] + alb[c3(i, k, j)])),
                                            cst::cpovcv) -
                         pb[c3(i, k, j)];
      }

  // ---- CPPWRF base state and state
  for (int k = 0; k <= nz; ++k) {
    b.znw[k] = znw(k + 1); b.c1f[k] = c1f(k + 1); b.c2f[k] = c2f(k + 1); b.c3f[k] = c3f(k + 1); b.c4f[k] = c4f(k + 1);
  }
  for (int k = 0; k < nz; ++k) {
    b.znu[k] = znu(k + 1); b.dnw[k] = dnw(k + 1); b.rdnw[k] = rdnw(k + 1);
    b.c1h[k] = c1h(k + 1); b.c2h[k] = c2h(k + 1); b.c3h[k] = c3h(k + 1); b.c4h[k] = c4h(k + 1);
    b.dn[k] = dn(k + 1); b.rdn[k] = rdn(k + 1); b.fnm[k] = fnm(k + 1); b.fnp[k] = fnp(k + 1);  // w levels, dn(1) = 0
  }
  b.dn[nz] = 0; b.rdn[nz] = 0; b.fnm[nz] = 0; b.fnp[nz] = 0;
  b.cf1 = cf1; b.cf2 = cf2; b.cf3 = cf3; b.cfn = cfn; b.cfn1 = cfn1;
  b.p_top = p_top;
  b.rdx = rdx; b.rdy = rdy;
  b.upload_levels();

  auto wrap = [](int a, int n, BC bc) {
    if (bc == BC::Periodic) return ((a % n) + n) % n;
    return std::clamp(a, 0, n - 1);
  };
  auto m2 = [&](const std::vector<double>& v) {
    return [&v, &gr, nx, jm, wrap](int i, int j, int) {
      return v[std::size_t(jm == 1 ? 0 : wrap(j, jm, gr.bc_y)) * nx + wrap(i, nx, gr.bc_x)];
    };
  };
  auto m3 = [&](const std::vector<double>& v) {
    return [&v, &gr, nx, jm, wrap](int i, int j, int k) {
      return v[(std::size_t(k) * jm + (jm == 1 ? 0 : wrap(j, jm, gr.bc_y))) * nx + wrap(i, nx, gr.bc_x)];
    };
  };
  b.ht.fill_host_halo(m2(ht));
  b.mub.fill_host_halo(m2(mub));
  b.pb.fill_host_halo(m3(pb));
  b.alb.fill_host_halo(m3(alb));
  b.t_init.fill_host_halo(m3(t_init));
  b.phb.fill_host_halo(m3(phb));
  // WRF's 1D qv_base (vertical_diffusion_mp diffuses qv - qv_base)
  b.qv_base.fill_host_halo([&](int, int, int k) { return double(b.qv_base_1d[std::clamp(k, 0, nz - 1)]); });
  for (BaseField2* f : {&b.msftx, &b.msfty, &b.msfux, &b.msfuy, &b.msfvx, &b.msfvy, &b.msfvx_inv, &b.cosa})
    f->fill_host_halo([](int, int, int) { return 1.0; });
  for (BaseField2* f : {&b.f, &b.e, &b.sina}) f->fill_host_halo([](int, int, int) { return 0.0; });

  m.st = std::make_unique<dyn::State>(gr);
  m.tend = std::make_unique<dyn::Tendencies>(gr);
  dyn::State& s = *m.st;
  for (const std::string& n : split_list(cfg.get_string("moist.species", "qv"))) {
    s.add_scalar(gr, n, true);
    m.tend->add_scalar(gr, n);
  }
  for (const std::string& n : split_list(cfg.get_string("tracers.names", ""))) {
    s.add_scalar(gr, n, false);
    m.tend->add_scalar(gr, n);
  }
  if (cfg.get_int("dyn.km_opt", 1) == 2 && s.P_TKE < 0) {
    s.add_scalar(gr, "tke", false);
    m.tend->add_scalar(gr, "tke");
  }
  if (s.P_QV < 0) throw std::runtime_error(name + " needs qv in moist.species (WRF carries qv for every mp_physics)");
  auto mi = [&](int i, int j, int k) { return (std::size_t(k) * jm + (jm == 1 ? 0 : j)) * nx + i; };
  s.u_2.fill_host([&](int i, int j, int k) { return u1[(std::size_t(k) * jm + (jm == 1 ? 0 : j)) * nxs + i]; });
  s.v_2.fill_host([&](int i, int j, int k) { return v1[(std::size_t(k) * jvm + (jvm == 1 ? 0 : j)) * nx + i]; });
  s.w_2.fill_host([](int, int, int) { return 0.0; });
  s.t_2.fill_host([&](int i, int j, int k) { return t_2[mi(i, j, k)]; });
  s.ph_2.fill_host([&](int i, int j, int k) { return ph_1[mi(i, j, k)]; });
  s.mu_2.fill_host([&](int i, int j, int) { return mu_1[std::size_t(jm == 1 ? 0 : j) * nx + i]; });
  s.al.fill_host([&](int i, int j, int k) { return al[mi(i, j, k)]; });
  s.p.fill_host([&](int i, int j, int k) { return p[mi(i, j, k)]; });
  for (int q = 0; q < s.n_scalars(); ++q) {
    if (q == s.P_QV) s.q_2[q].fill_host([&](int i, int j, int k) { return qv[mi(i, j, k)]; });
    else s.q_2[q].fill_host([](int, int, int) { return 0.0; });
    fill_halo(s.q_2[q], gr);
    s.q_1[q].copy_from(s.q_2[q]);
  }
  fill_halos({&s.u_2, &s.v_2, &s.w_2, &s.t_2, &s.ph_2, &s.al, &s.p}, gr);
  fill_halo(s.mu_2, gr);
  s.u_1.copy_from(s.u_2); s.v_1.copy_from(s.v_2); s.t_1.copy_from(s.t_2); s.ph_1.copy_from(s.ph_2);
  s.mu_1.copy_from(s.mu_2);
  dyn::set_w_surface(s.w_1, s.u_1, s.v_1, b, gr, true);
  dyn::set_w_surface(s.w_2, s.u_2, s.v_2, b, gr, true);
  m.ap.hypsometric_opt = hypsometric_opt;
  m.start_date = "0001-01-01_00:00:00";  // the start date of WRF's test/em_* namelists (run.start_date overrides)
}

}  // namespace cppwrf::cases
