#include "cases/hill2d.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
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
  double ts = 288.0;
  double qvs = 0.0;
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
  if (cfg.has("case.sounding_file")) {
    std::string path = cfg.get_string("case.sounding_file");
    std::ifstream test(path);
    if (!test.good()) {
      throw std::runtime_error("Cannot open specified sounding_file: " + path);
    }
    return read_sounding_file(path);
  }
  const std::vector<std::string> candidates = {"configs/input_sounding_hill2d", "input_sounding"};
  for (const auto& c : candidates) {
    std::ifstream test(c);
    if (test.good()) {
      return read_sounding_file(c);
    }
  }
  throw std::runtime_error("Could not locate input_sounding file for hill2d (tried configs/input_sounding_hill2d, input_sounding)");
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

void init_hill2d(const io::Config& cfg, const std::vector<double>& eta, const Grid& g,
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

  // Hydrostatic sounding integration as in WRF get_sounding
  const double p_surf = 100.0 * raw.ps;
  const double th_surf = raw.ts;
  const double qvf = 1.0 + cst::rvovrd * qv[0];
  const double rho_surf = 1.0 / ((cst::r_d / cst::p1000mb) * th_surf * qvf * std::pow(p_surf / cst::p1000mb, cst::cvpm));

  std::vector<double> rho_in(nl, 0.0);
  std::vector<double> pm_in(nl, 0.0);
  rho_in[0] = rho_surf;
  pm_in[0] = p_surf;

  double dz = h[0];
  for (int it = 0; it < 10; ++it) {
    pm_in[0] = p_surf - 0.5 * dz * (rho_surf + rho_in[0]) * cst::g * (1.0 + qv[0]);
    rho_in[0] = 1.0 / ((cst::r_d / cst::p1000mb) * th[0] * qvf * std::pow(pm_in[0] / cst::p1000mb, cst::cvpm));
  }

  for (int k = 1; k < nl; ++k) {
    rho_in[k] = rho_in[k - 1];
    dz = h[k] - h[k - 1];
    const double qvf1 = 0.5 * (2.0 + (qv[k - 1] + qv[k]));
    const double qvf_k = 1.0 + cst::rvovrd * qv[k];
    for (int it = 0; it < 10; ++it) {
      pm_in[k] = pm_in[k - 1] - 0.5 * dz * (rho_in[k] + rho_in[k - 1]) * cst::g * qvf1;
      rho_in[k] = 1.0 / ((cst::r_d / cst::p1000mb) * th[k] * qvf_k * std::pow(pm_in[k] / cst::p1000mb, cst::cvpm));
    }
  }

  // Model top pressure
  const double p_top_calc = interp_0(pm_in, h, g.ztop);
  const double p_top = cfg.get_double("domain.p_top", p_top_calc);
  b.p_top = p_top;

  b.set_levels(eta);
  b.upload_levels();

  // Terrain height formula: h(i) = hm / (1 + ((i - icm)/xa)^2)
  const double hm = cfg.get_double("case.hm", 100.0);
  const double xa = cfg.get_double("case.xa", 5.0);
  const double icm = cfg.get_double("case.icm", double(g.nx / 2));

  auto ht_func = [&](int i) {
    double x = double(i);
    if (g.bc_x == BC::Periodic) {
      x = double((i % g.nx + g.nx) % g.nx);
    }
    return hm / (1.0 + std::pow((x - icm) / xa, 2));
  };

  b.ht.fill_host_halo([&](int i, int, int) { return ht_func(i); });

  // Base state columns
  const int nx = g.nx, nz = g.nz;
  std::vector<double> mub_col(nx);
  std::vector<double> pb_col(nz * nx);
  std::vector<double> alb_col(nz * nx);
  std::vector<double> t_init_col(nz * nx);
  std::vector<double> phb_col((nz + 1) * nx);

  for (int i = 0; i < nx; ++i) {
    const double zsfc = ht_func(i);
    const double psurf = interp_0(pm_in, h, zsfc);
    mub_col[i] = psurf - p_top;
    phb_col[0 * nx + i] = cst::g * zsfc;

    for (int k = 0; k < nz; ++k) {
      const double p_lvl = b.znu[k] * mub_col[i] + p_top;
      pb_col[k * nx + i] = p_lvl;
      const double th_lvl = interp_0(th, pm_in, p_lvl);
      t_init_col[k * nx + i] = th_lvl - cst::t0;
      alb_col[k * nx + i] = (cst::r_d / cst::p1000mb) * th_lvl * std::pow(p_lvl / cst::p1000mb, cst::cvpm);
    }

    for (int k = 0; k < nz; ++k) {
      phb_col[(k + 1) * nx + i] = phb_col[k * nx + i] - b.dnw[k] * mub_col[i] * alb_col[k * nx + i];
    }
  }

  // Fill BaseState arrays with boundary-aware halo wrapping
  auto col_idx = [&](int i) {
    if (g.bc_x == BC::Periodic) return (i % nx + nx) % nx;
    return std::clamp(i, 0, nx - 1);
  };
  b.mub.fill_host_halo([&](int i, int, int) { return mub_col[col_idx(i)]; });
  b.pb.fill_host_halo([&](int i, int, int k) { return pb_col[k * nx + col_idx(i)]; });
  b.alb.fill_host_halo([&](int i, int, int k) { return alb_col[k * nx + col_idx(i)]; });
  b.t_init.fill_host_halo([&](int i, int, int k) { return t_init_col[k * nx + col_idx(i)]; });
  b.phb.fill_host_halo([&](int i, int, int k) { return phb_col[k * nx + col_idx(i)]; });
  b.qv_base.fill_host_halo([](int, int, int) { return 0.0; });

  // Sounding wind interpolation (WRF interp_0 from sounding u, or override by case.u0)
  const bool has_u0 = cfg.has("case.u0");
  const double u0_val = cfg.get_double("case.u0", 10.0);
  std::vector<double> u_stag((nx + 1) * nz);
  for (int i = 0; i <= nx; ++i) {
    const double zsfc_u = 0.5 * (ht_func(i) + ht_func(i > 0 ? i - 1 : (g.bc_x == BC::Periodic ? nx - 1 : 0)));
    const double psurf_u = interp_0(pm_in, h, zsfc_u);
    for (int k = 0; k < nz; ++k) {
      if (has_u0) {
        u_stag[k * (nx + 1) + i] = u0_val;
      } else {
        const double p_lvl = b.znu[k] * (psurf_u - p_top) + p_top;
        u_stag[k * (nx + 1) + i] = interp_0(u, pm_in, p_lvl);
      }
    }
  }

  // Fill State s
  s.t_2.fill_host([&](int i, int, int k) { return t_init_col[k * nx + i]; });
  s.al.fill_host([](int, int, int) { return 0.0; });
  s.alt.fill_host([&](int i, int, int k) { return alb_col[k * nx + i]; });
  s.ph_2.fill_host([](int, int, int) { return 0.0; });
  s.p.fill_host([](int, int, int) { return 0.0; });
  s.mu_2.fill_host([](int, int, int) { return 0.0; });
  s.u_2.fill_host([&](int i, int, int k) { return u_stag[k * (nx + 1) + i]; });
  s.v_2.fill_host([](int, int, int) { return 0.0; });
  s.w_2.fill_host([](int, int, int) { return 0.0; });

  for (Field3* f : {&s.t_2, &s.al, &s.alt, &s.ph_2, &s.u_2, &s.v_2, &s.w_2, &s.p}) fill_halo(*f, g);
  fill_halo(s.mu_2, g);

  s.u_1.copy_from(s.u_2);
  s.v_1.copy_from(s.v_2);
  s.w_1.copy_from(s.w_2);
  s.t_1.copy_from(s.t_2);
  s.ph_1.copy_from(s.ph_2);
  s.mu_1.copy_from(s.mu_2);
}

}  // namespace cppwrf::cases
