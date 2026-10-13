#include "cases/quarter_ss.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/constants.hpp"
#include "core/halo.hpp"

namespace cppwrf::cases {

namespace {

struct SoundingPoint {
  double h = 0.0;
  double th = 0.0;
  double qv = 0.0;
  double u = 0.0;
  double v = 0.0;
};

struct RawSounding {
  double ps = 1000.0;
  double ts = 300.0;
  double qvs = 14.0;
  std::vector<SoundingPoint> pts;
};

RawSounding read_sounding_file(const std::string& path) {
  std::ifstream f(path);
  if (!f.is_open()) {
    throw std::runtime_error("Cannot open input_sounding: " + path);
  }
  std::string line;
  RawSounding snd;
  bool header_read = false;
  while (std::getline(f, line)) {
    size_t first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos || line[first] == '#' || line[first] == '!') continue;
    std::istringstream iss(line);
    if (!header_read) {
      if (iss >> snd.ps >> snd.ts >> snd.qvs) {
        header_read = true;
      }
    } else {
      SoundingPoint pt;
      if (iss >> pt.h >> pt.th >> pt.qv >> pt.u >> pt.v) {
        snd.pts.push_back(pt);
      }
    }
  }
  if (!header_read || snd.pts.empty()) {
    throw std::runtime_error("Failed to parse sounding data from: " + path);
  }
  return snd;
}

RawSounding load_sounding(const io::Config& cfg) {
  std::vector<std::string> candidates;
  if (cfg.has("case.sounding_file")) {
    candidates.push_back(cfg.get_string("case.sounding_file"));
  }
  candidates.push_back("configs/input_sounding_qss");
  candidates.push_back("/home/ostriquetrum/wrf_ref/run_qss/input_sounding");
  candidates.push_back("/home/ostriquetrum/wrf_ref/WRF_qss/test/em_quarter_ss/input_sounding");
  candidates.push_back("input_sounding");

  for (const auto& c : candidates) {
    std::ifstream test(c);
    if (test.good()) {
      return read_sounding_file(c);
    }
  }
  throw std::runtime_error("Could not locate input_sounding file for quarter_ss");
}

// 1D linear interpolation matching WRF interp_0 in module_init_utilities.F
double interp_0(const std::vector<double>& v_in, const std::vector<double>& z_in, double z_out) {
  const int n = int(z_in.size());
  if (n < 2) throw std::runtime_error("interp_0: need >= 2 points");
  if (z_in.back() > z_in.front()) {
    // Increasing z coordinate
    if (z_out >= z_in.back()) {
      double w2 = (z_in.back() - z_out) / (z_in.back() - z_in[n - 2]);
      return (1.0 - w2) * v_in.back() + w2 * v_in[n - 2];
    } else if (z_out <= z_in.front()) {
      double w2 = (z_in[1] - z_out) / (z_in[1] - z_in[0]);
      return (1.0 - w2) * v_in[1] + w2 * v_in[0];
    } else {
      for (int kp = n - 1; kp >= 1; --kp) {
        if (z_in[kp] >= z_out && z_in[kp - 1] <= z_out) {
          double w2 = (z_out - z_in[kp]) / (z_in[kp - 1] - z_in[kp]);
          return (1.0 - w2) * v_in[kp] + w2 * v_in[kp - 1];
        }
      }
    }
  } else {
    // Decreasing z coordinate (such as pressure)
    if (z_out <= z_in.back()) {
      double w2 = (z_in.back() - z_out) / (z_in.back() - z_in[n - 2]);
      return (1.0 - w2) * v_in.back() + w2 * v_in[n - 2];
    } else if (z_out >= z_in.front()) {
      double w2 = (z_in[1] - z_out) / (z_in[1] - z_in[0]);
      return (1.0 - w2) * v_in[1] + w2 * v_in[0];
    } else {
      for (int kp = n - 1; kp >= 1; --kp) {
        if (z_in[kp] <= z_out && z_in[kp - 1] >= z_out) {
          double w2 = (z_out - z_in[kp]) / (z_in[kp - 1] - z_in[kp]);
          return (1.0 - w2) * v_in[kp] + w2 * v_in[kp - 1];
        }
      }
    }
  }
  return v_in.back();
}

}  // namespace

void init_quarter_ss(const io::Config& cfg, const std::vector<double>& eta, const Grid& g,
                     dyn::BaseState& b, dyn::State& s) {
  RawSounding raw = load_sounding(cfg);
  const int nl = int(raw.pts.size());
  std::vector<double> h(nl), th(nl), qv(nl), u(nl), v(nl);
  for (int k = 0; k < nl; ++k) {
    h[k] = raw.pts[k].h;
    th[k] = raw.pts[k].th;
    qv[k] = 0.001 * raw.pts[k].qv;
    u[k] = raw.pts[k].u;
    v[k] = raw.pts[k].v;
  }

  const double p_surf = 100.0 * raw.ps;
  const double th_surf = raw.ts;

  // 1. Dry hydrostatic sounding integration (for dry base state)
  const double rho_surf_dry = 1.0 / ((cst::r_d / cst::p1000mb) * th_surf * std::pow(p_surf / cst::p1000mb, cst::cvpm));
  std::vector<double> p_dry(nl, 0.0), rho_dry(nl, 0.0);
  rho_dry[0] = rho_surf_dry;
  double dz = h[0];
  for (int it = 0; it < 10; ++it) {
    p_dry[0] = p_surf - 0.5 * dz * (rho_surf_dry + rho_dry[0]) * cst::g;
    rho_dry[0] = 1.0 / ((cst::r_d / cst::p1000mb) * th[0] * std::pow(p_dry[0] / cst::p1000mb, cst::cvpm));
  }
  for (int k = 1; k < nl; ++k) {
    rho_dry[k] = rho_dry[k - 1];
    dz = h[k] - h[k - 1];
    for (int it = 0; it < 10; ++it) {
      p_dry[k] = p_dry[k - 1] - 0.5 * dz * (rho_dry[k] + rho_dry[k - 1]) * cst::g;
      rho_dry[k] = 1.0 / ((cst::r_d / cst::p1000mb) * th[k] * std::pow(p_dry[k] / cst::p1000mb, cst::cvpm));
    }
  }

  // 2. Moist hydrostatic sounding integration (for moist full state)
  const double qvf_sfc = 1.0 + cst::rvovrd * qv[0];
  const double rho_surf_moist = 1.0 / ((cst::r_d / cst::p1000mb) * th_surf * qvf_sfc * std::pow(p_surf / cst::p1000mb, cst::cvpm));
  std::vector<double> pm(nl, 0.0), rho_m(nl, 0.0);
  rho_m[0] = rho_surf_moist;
  dz = h[0];
  for (int it = 0; it < 10; ++it) {
    pm[0] = p_surf - 0.5 * dz * (rho_surf_moist + rho_m[0]) * cst::g * (1.0 + qv[0]);
    rho_m[0] = 1.0 / ((cst::r_d / cst::p1000mb) * th[0] * qvf_sfc * std::pow(pm[0] / cst::p1000mb, cst::cvpm));
  }
  for (int k = 1; k < nl; ++k) {
    rho_m[k] = rho_m[k - 1];
    dz = h[k] - h[k - 1];
    const double qvf1 = 0.5 * (2.0 + (qv[k - 1] + qv[k]));
    const double qvf_k = 1.0 + cst::rvovrd * qv[k];
    for (int it = 0; it < 10; ++it) {
      pm[k] = pm[k - 1] - 0.5 * dz * (rho_m[k] + rho_m[k - 1]) * cst::g * qvf1;
      rho_m[k] = 1.0 / ((cst::r_d / cst::p1000mb) * th[k] * qvf_k * std::pow(pm[k] / cst::p1000mb, cst::cvpm));
    }
  }

  // 3. Dry pressure profile pd integrated downward from top of moist sounding
  std::vector<double> pd(nl, 0.0);
  pd[nl - 1] = pm[nl - 1];
  for (int k = nl - 2; k >= 0; --k) {
    dz = h[k + 1] - h[k];
    pd[k] = pd[k + 1] + 0.5 * dz * (rho_m[k] + rho_m[k + 1]) * cst::g;
  }

  // Model top pressure
  const double p_top_calc = interp_0(p_dry, h, g.ztop);
  const double p_top = cfg.get_double("domain.p_top", p_top_calc);
  b.p_top = p_top;

  b.set_levels(eta);
  b.upload_levels();

  const int nx = g.nx, ny = g.ny, nz = g.nz;

  // Base state columns (flat terrain, horizontally uniform)
  const double p_surf_dry = interp_0(p_dry, h, 0.0);
  const double mub = p_surf_dry - p_top;

  std::vector<double> pb(nz), t_init(nz), alb(nz), phb(nz + 1);
  phb[0] = 0.0;
  for (int k = 0; k < nz; ++k) {
    pb[k] = b.znu[k] * mub + p_top;
    t_init[k] = interp_0(th, p_dry, pb[k]) - cst::t0;
    alb[k] = (cst::r_d / cst::p1000mb) * (t_init[k] + cst::t0) * std::pow(pb[k] / cst::p1000mb, cst::cvpm);
    phb[k + 1] = phb[k] - b.dnw[k] * mub * alb[k];
  }

  // Moist unperturbed 1D state
  const double pd_surf = interp_0(pd, h, 0.0);
  const double mu1 = pd_surf - p_top - mub;

  std::vector<double> qv_1d(nz), t_1d(nz);
  for (int k = 0; k < nz; ++k) {
    const double p_lvl = b.znu[k] * (pd_surf - p_top) + p_top;
    qv_1d[k] = interp_0(qv, pd, p_lvl);
    t_1d[k] = interp_0(th, pd, p_lvl) - cst::t0;
  }

  std::vector<double> p_1d(nz), alt_1d(nz), al_1d(nz), ph1_1d(nz + 1);
  // Integrate hydrostatic p' down from model top
  int kk = nz - 1;
  double qvf1 = qv_1d[kk];
  double qvf2 = 1.0 / (1.0 + qvf1);
  qvf1 = qvf1 * qvf2;
  p_1d[kk] = -0.5 * (mu1 + qvf1 * mub) / b.rdnw[kk] / qvf2;
  double qvf = 1.0 + cst::rvovrd * qv_1d[kk];
  alt_1d[kk] = (cst::r_d / cst::p1000mb) * (t_1d[kk] + cst::t0) * qvf * std::pow((p_1d[kk] + pb[kk]) / cst::p1000mb, cst::cvpm);
  al_1d[kk] = alt_1d[kk] - alb[kk];

  for (kk = nz - 2; kk >= 0; --kk) {
    qvf1 = 0.5 * (qv_1d[kk] + qv_1d[kk + 1]);
    qvf2 = 1.0 / (1.0 + qvf1);
    qvf1 = qvf1 * qvf2;
    p_1d[kk] = p_1d[kk + 1] - (mu1 + qvf1 * mub) / qvf2 / b.rdn[kk + 1];
    qvf = 1.0 + cst::rvovrd * qv_1d[kk];
    alt_1d[kk] = (cst::r_d / cst::p1000mb) * (t_1d[kk] + cst::t0) * qvf * std::pow((p_1d[kk] + pb[kk]) / cst::p1000mb, cst::cvpm);
    al_1d[kk] = alt_1d[kk] - alb[kk];
  }

  ph1_1d[0] = 0.0;
  for (int k = 0; k < nz; ++k) {
    ph1_1d[k + 1] = ph1_1d[k] - b.dnw[k] * ((mub + mu1) * al_1d[k] + mu1 * alb[k]);
  }

  // 4. 3D Thermal bubble perturbation (case quarter_ss)
  const double delt = cfg.get_double("case.delt", 3.0);
  const int nxc = int(cfg.get_int("case.nxc", nx / 2));  // Fortran nxc = ide/2 = 20 -> 1-based index 20 (0-based 19)
  const int nyc = int(cfg.get_int("case.nyc", ny / 2));
  const double pi = std::numbers::pi;

  const size_t ncells = size_t(nx) * ny * nz;
  const size_t ncells_w = size_t(nx) * ny * (nz + 1);
  auto I = [&](int i, int j, int k) { return (size_t(k) * ny + j) * nx + i; };

  std::vector<double> t_3d(ncells), alt_3d(ncells), al_3d(ncells), ph1_3d(ncells_w);

  for (int j = 0; j < ny; ++j) {
    const double yrad = g.dy * double(j + 1 - nyc) / 10000.0;
    for (int i = 0; i < nx; ++i) {
      const double xrad = g.dx * double(i + 1 - nxc) / 10000.0;
      for (int k = 0; k < nz; ++k) {
        t_3d[I(i, j, k)] = t_1d[k];
        alt_3d[I(i, j, k)] = alt_1d[k];
        al_3d[I(i, j, k)] = al_1d[k];

        const double zm = 0.5 * (ph1_1d[k] + ph1_1d[k + 1] + phb[k] + phb[k + 1]) / cst::g;
        const double zrad = (zm - 1500.0) / 1500.0;
        const double rad = std::sqrt(xrad * xrad + yrad * yrad + zrad * zrad);
        if (rad <= 1.0) {
          const double dt = delt * std::pow(std::cos(0.5 * pi * rad), 2);
          t_3d[I(i, j, k)] += dt;
          const double qvf_local = 1.0 + cst::rvovrd * qv_1d[k];
          alt_3d[I(i, j, k)] = (cst::r_d / cst::p1000mb) * (t_3d[I(i, j, k)] + cst::t0) * qvf_local *
                               std::pow((p_1d[k] + pb[k]) / cst::p1000mb, cst::cvpm);
          al_3d[I(i, j, k)] = alt_3d[I(i, j, k)] - alb[k];
        }
      }

      // Re-integrate ph1 hydrostatically for column (i, j)
      ph1_3d[I(i, j, 0)] = 0.0;
      for (int k = 0; k < nz; ++k) {
        ph1_3d[I(i, j, k + 1)] = ph1_3d[I(i, j, k)] - b.dnw[k] * ((mub + mu1) * al_3d[I(i, j, k)] + mu1 * alb[k]);
      }
    }
  }

  // 5. Winds from sounding
  const double p_surf_moist = interp_0(pm, h, 0.0);
  std::vector<double> u_lvl(nz), v_lvl(nz);
  for (int k = 0; k < nz; ++k) {
    const double p_lvl = b.znu[k] * (p_surf_moist - p_top) + p_top;
    u_lvl[k] = interp_0(u, pm, p_lvl);
    v_lvl[k] = interp_0(v, pm, p_lvl);
  }

  // 6. Fill BaseState (flat = true, uniform arrays)
  b.mub.fill_host([=](int, int, int) { return mub; });
  b.pb.fill_host([&](int, int, int k) { return pb[k]; });
  b.alb.fill_host([&](int, int, int k) { return alb[k]; });
  b.t_init.fill_host([&](int, int, int k) { return t_init[k]; });
  b.phb.fill_host([&](int, int, int k) { return phb[k]; });
  b.ht.fill_host([](int, int, int) { return 0.0; });
  b.qv_base.fill_host([&](int, int, int k) { return qv_1d[k]; });

  // 7. Fill State s (t_2 holds theta_m - t0 in CPPWRF moist formulation)
  s.t_2.fill_host([&](int i, int j, int k) {
    const double qvf_k = 1.0 + cst::rvovrd * qv_1d[k];
    return (t_3d[I(i, j, k)] + cst::t0) * qvf_k - cst::t0;
  });
  s.al.fill_host([&](int i, int j, int k) { return al_3d[I(i, j, k)]; });
  s.alt.fill_host([&](int i, int j, int k) { return alt_3d[I(i, j, k)]; });
  s.p.fill_host([&](int, int, int k) { return p_1d[k]; });
  s.ph_2.fill_host([&](int i, int j, int k) { return ph1_3d[I(i, j, k)]; });
  s.mu_2.fill_host([=](int, int, int) { return mu1; });
  s.u_2.fill_host([&](int, int, int k) { return u_lvl[k]; });
  s.v_2.fill_host([&](int, int, int k) { return v_lvl[k]; });
  s.w_2.fill_host([](int, int, int) { return 0.0; });

  if (s.P_QV >= 0) s.q_2[s.P_QV].fill_host([&](int, int, int k) { return qv_1d[k]; });
  if (s.P_QC >= 0) s.q_2[s.P_QC].fill_host([](int, int, int) { return 0.0; });
  if (s.P_QR >= 0) s.q_2[s.P_QR].fill_host([](int, int, int) { return 0.0; });

  // Fill halos
  for (Field3* f : {&s.t_2, &s.al, &s.alt, &s.ph_2, &s.p, &s.u_2, &s.v_2, &s.w_2}) fill_halo(*f, g);
  fill_halo(s.mu_2, g);
  if (s.P_QV >= 0) fill_halo(s.q_2[s.P_QV], g);
  if (s.P_QC >= 0) fill_halo(s.q_2[s.P_QC], g);
  if (s.P_QR >= 0) fill_halo(s.q_2[s.P_QR], g);

  // Copy stage values to time t (_1)
  s.u_1.copy_from(s.u_2);
  s.v_1.copy_from(s.v_2);
  s.w_1.copy_from(s.w_2);
  s.t_1.copy_from(s.t_2);
  s.ph_1.copy_from(s.ph_2);
  s.mu_1.copy_from(s.mu_2);
  if (s.P_QV >= 0) s.q_1[s.P_QV].copy_from(s.q_2[s.P_QV]);
  if (s.P_QC >= 0) s.q_1[s.P_QC].copy_from(s.q_2[s.P_QC]);
  if (s.P_QR >= 0) s.q_1[s.P_QR].copy_from(s.q_2[s.P_QR]);
}

}  // namespace cppwrf::cases
