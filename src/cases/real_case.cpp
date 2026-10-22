// Real-data start (WRF dyn_em/start_em.F for a single domain, no restart, use_theta_m = 1):
// grid, base state, map factors, Coriolis, terrain, initial state and lateral boundaries from
// real.exe's wrfinput_d01 / wrfbdy_d01 (read by the io_netcdf readers).
//   - 1D vertical arrays (znw, znu, dnw, dn, rdnw, rdn, fnm, fnp, c1..c4, cf1..cfn) as stored.
//   - pb = c3h*mub + c4h + p_top; alb from t_init (start_em, max_dom = 1); phb, mub as stored.
//   - prognostic theta: THM = (T + t0)(1 + Rv/Rd qv) - t0 (Registry dataname "thm" for t).
//   - al from ph (hypsometric_opt of the namelist, default 2), p from al and theta_m.
//   - w = 0 in the input -> set_w_surface with the znw^2 fill (start_em w_needs_to_be_set).
#include <algorithm>
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
#include "wrf_io/global_atts.hpp"

#include <netcdf.h>
#endif

namespace cppwrf::cases {

#ifdef CPPWRF_HAVE_NETCDF
namespace {
// Base-state parameters (scalar variables T00, P00, TLP, TISO, TLP_STRAT, P_STRAT of wrfinput)
struct BaseParams { double t00, p00, a, tiso, a_strat, p_strat; };
BaseParams read_base_params(const std::string& path) {
  int nc = 0;
  if (nc_open(path.c_str(), NC_NOWRITE, &nc) != NC_NOERR) throw std::runtime_error("real: cannot open " + path);
  auto get = [&](const char* name) {
    int v = 0;
    double x = 0;
    if (nc_inq_varid(nc, name, &v) != NC_NOERR || nc_get_var_double(nc, v, &x) != NC_NOERR) {
      nc_close(nc);
      throw std::runtime_error(std::string("real: wrfinput lacks base-state parameter ") + name);
    }
    return x;
  };
  BaseParams b{get("T00"), get("P00"), get("TLP"), get("TISO"), get("TLP_STRAT"), get("P_STRAT")};
  nc_close(nc);
  return b;
}

// A whole (Time = 1) variable of wrfinput as doubles; empty if absent.
std::vector<double> read_var_if_present(const std::string& path, const char* name, std::size_t n) {
  int nc = 0, v = 0;
  std::vector<double> out;
  if (nc_open(path.c_str(), NC_NOWRITE, &nc) != NC_NOERR) throw std::runtime_error("real: cannot open " + path);
  if (nc_inq_varid(nc, name, &v) == NC_NOERR) {
    out.resize(n);
    if (nc_get_var_double(nc, v, out.data()) != NC_NOERR) out.clear();
  }
  nc_close(nc);
  return out;
}

// "YYYY-MM-DD_HH:MM:SS" -> seconds between two times (proleptic Gregorian, UTC)
double seconds_between(const std::string& a, const std::string& b) {
  auto secs = [](const std::string& s) {
    int Y, M, D, h, mi, se;
    if (std::sscanf(s.c_str(), "%d-%d-%d_%d:%d:%d", &Y, &M, &D, &h, &mi, &se) != 6) {
      throw std::runtime_error("seconds_between: invalid date string '" + s + "' (expected YYYY-MM-DD_HH:MM:SS)");
    }
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
  m.xlat = d.xlat;
  m.xlong = d.xlong;
  m.hybrid_opt = d.hybrid_opt;
  // map_proj 6 (lat-lon): WRF's curvature terms take the tan(lat) form (module_big_step_utilities_em.F curvature),
  // not ported
  if (d.map_proj == 6) throw std::runtime_error("MAP_PROJ = 6 (lat-lon) is not supported");
  m.map_proj = d.map_proj;
  m.cen_lat = d.cen_lat;
  m.cen_lon = d.cen_lon;
  m.truelat1 = d.truelat1;
  m.truelat2 = d.truelat2;
  m.moad_cen_lat = d.moad_cen_lat;
  m.stand_lon = d.stand_lon;
  if (!d.start_date.empty()) m.start_date = d.start_date;
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
  b.p_top = d.p_top;
  {  // WRF uses RDX, RDY and 1D *_BASE profiles of the input (registry 'i')
    const std::string wi = dir + "/" + cfg.get_string("real.wrfinput", "wrfinput_d01");
    const std::vector<double> rdx = read_var_if_present(wi, "RDX", 1), rdy = read_var_if_present(wi, "RDY", 1);
    if (!rdx.empty()) b.rdx = rdx[0];
    if (!rdy.empty()) b.rdy = rdy[0];
    auto load_base1d = [&](const char* name, Array1& dst) {
      const std::vector<double> v = read_var_if_present(wi, name, std::size_t(nz));
      for (int k = 0; k < nz && std::size_t(k) < v.size(); ++k) dst[k] = real(v[k]);
    };
    load_base1d("U_BASE", b.u_base);
    load_base1d("V_BASE", b.v_base);
    load_base1d("T_BASE", b.t_base);
    load_base1d("QV_BASE", b.qv_base_1d);
    load_base1d("Z_BASE", b.z_base);
  }
  b.upload_levels();

  auto cl = [](int a, int lo, int hi) { return a < lo ? lo : (a > hi ? hi : a); };
  auto m2 = [&](const std::vector<double>& v) {  // mass-point 2D, flat halo
    return [&v, nx, ny, cl](int i, int j, int) { return v[std::size_t(cl(j, 0, ny - 1)) * nx + cl(i, 0, nx - 1)]; };
  };
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
  {  // WRF reads MF_VX_INV from the input (registry 'i'); real.exe of WRF v4 writes it
    const std::string wi = dir + "/" + cfg.get_string("real.wrfinput", "wrfinput_d01");
    const std::vector<double> inv = read_var_if_present(wi, "MF_VX_INV", std::size_t(nx) * nys);
    if (inv.empty()) throw std::runtime_error("real: " + wi + " lacks MF_VX_INV (needs real.exe of WRF v4)");
    b.msfvx_inv.fill_host_halo(v2(inv));
  }

  auto m3 = [&](const std::vector<double>& v, int nk) {
    return [&v, nx, ny, nk, cl](int i, int j, int k) {
      return v[(std::size_t(cl(k, 0, nk - 1)) * ny + cl(j, 0, ny - 1)) * nx + cl(i, 0, nx - 1)];
    };
  };
  // start_em (input_from_file, not restart): the base state is recomputed from the terrain and
  // the base-state parameters (start_em.F 590-620), not taken from wrfinput: p_surf, pb, t_init,
  // alb, mub, phb. Then (max_dom = 1) pb, alb again from mub. hypsometric_opt is the namelist's (Registry
  // default 2); input_wrf.F resets it to 1 only when the input has no valid HYPSOMETRIC_OPT attribute.
  const auto ha = d.attrs_int.find("HYPSOMETRIC_OPT");
  const bool hypso_attr = ha != d.attrs_int.end() && (ha->second == 1 || ha->second == 2);
  const int hypso = hypso_attr ? int(cfg.get_int("dyn.hypsometric_opt", 2)) : 1;
  const BaseParams bp = read_base_params(dir + "/" + cfg.get_string("real.wrfinput", "wrfinput_d01"));
  // start_em / input_wrf: the namelist values must agree with the input file (wrf.exe stops otherwise)
  if (cfg.has("real.p_top_requested") && cfg.get_double("real.p_top_requested") != d.p_top)
    throw std::runtime_error("start_em: p_top from the namelist does not match p_top from the input file.");
  if (bp.t00 < 100. || bp.p00 < 10000.)
    throw std::runtime_error("start_em: did not find base state parameters in wrfinput.");
  if (cfg.has("real.iso_temp") && cfg.get_double("real.iso_temp") != bp.tiso)
    throw std::runtime_error("start_em: namelist iso_temp is not equal to iso_temp in wrfinput.");
  if (d.attrs_int.count("USE_THETA_M") && d.attrs_int.at("USE_THETA_M") != 1)
    throw std::runtime_error("use_theta_m values must be consistent: input file use_theta_m = " +
                             std::to_string(d.attrs_int.at("USE_THETA_M")) + ", CPPWRF use_theta_m = 1");
  std::vector<double> pb(std::size_t(nx) * ny * nz), alb(pb.size()), t_init(pb.size());
  std::vector<double> mub(std::size_t(nx) * ny), phb(std::size_t(nx) * ny * (nz + 1));
  {
    const double t00 = bp.t00, p00 = bp.p00, a = bp.a, tiso = bp.tiso, a_strat = bp.a_strat, p_strat = bp.p_strat;
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        const std::size_t ij = std::size_t(j) * nx + i;
        const double ht = d.hgt[ij];
        const double p_surf = p00 * std::exp(-t00 / a + std::pow((t00 / a) * (t00 / a) - 2. * cst::g * ht / a / cst::r_d, 0.5));  // WRF **0.5 is pow, not sqrt
        for (int k = 0; k < nz; ++k) {
          const std::size_t n = (std::size_t(k) * ny + j) * nx + i;
          pb[n] = d.c3h[k] * (p_surf - d.p_top) + d.c4h[k] + d.p_top;
          double temp = std::max(tiso, t00 + a * std::log(pb[n] / p00));
          if (pb[n] < p_strat) temp = tiso + a_strat * std::log(pb[n] / p_strat);
          t_init[n] = temp * std::pow(p00 / pb[n], cst::r_d / cst::cp) - cst::t0;
          alb[n] = (cst::r_d / cst::p1000mb) * (t_init[n] + cst::t0) * std::pow(pb[n] / cst::p1000mb, cst::cvpm);
        }
        mub[ij] = p_surf - d.p_top;
        phb[ij] = ht * cst::g;
        for (int k = 1; k <= nz; ++k) {
          const std::size_t n = (std::size_t(k) * ny + j) * nx + i, nm = (std::size_t(k - 1) * ny + j) * nx + i;
          if (hypso == 1) {
            phb[n] = phb[nm] - d.dnw[k - 1] * (d.c1h[k - 1] * mub[ij] + d.c2h[k - 1]) * alb[nm];
          } else {
            const double pfu = d.c3f[k] * mub[ij] + d.c4f[k] + d.p_top;
            const double pfd = d.c3f[k - 1] * mub[ij] + d.c4f[k - 1] + d.p_top;
            const double phm = d.c3h[k - 1] * mub[ij] + d.c4h[k - 1] + d.p_top;
            phb[n] = phb[nm] + alb[nm] * phm * std::log(pfd / pfu);
          }
        }
        for (int k = 0; k < nz; ++k) {  // max_dom = 1: pb, alb from mub and t_init
          const std::size_t n = (std::size_t(k) * ny + j) * nx + i;
          pb[n] = d.c3h[k] * mub[ij] + d.c4h[k] + d.p_top;
          alb[n] = (cst::r_d / cst::p1000mb) * (t_init[n] + cst::t0) * std::pow(pb[n] / cst::p1000mb, cst::cvpm);
        }
      }
  }
  b.mub.fill_host_halo(m2(mub));
  b.pb.fill_host_halo(m3(pb, nz));
  b.alb.fill_host_halo(m3(alb, nz));
  b.t_init.fill_host_halo(m3(t_init, nz));
  b.phb.fill_host_halo(m3(phb, nz + 1));
  b.qv_base.fill_host_halo([](int, int, int) { return 0.0; });

  // ---- state
  m.st = std::make_unique<dyn::State>(g);
  m.tend = std::make_unique<dyn::Tendencies>(g);
  dyn::State& s = *m.st;
  for (const std::string& n : split_list(cfg.get_string("moist.species", "qv"))) {
    s.add_scalar(g, n, true);
    m.tend->add_scalar(g, n);
  }
  for (const std::string& n : split_list(cfg.get_string("tracers.names", ""))) {
    s.add_scalar(g, n, false);
    m.tend->add_scalar(g, n);
  }
  if (cfg.get_int("dyn.km_opt", 1) == 2 && s.P_TKE < 0) {
    s.add_scalar(g, "tke", false);
    m.tend->add_scalar(g, "tke");
  }
  if (s.P_QV < 0) throw std::runtime_error("case = real needs qv in moist.species");
  const std::vector<double>& qv = d.qvapor;
  const std::string wi_path = dir + "/" + cfg.get_string("real.wrfinput", "wrfinput_d01");
  // t_2 = theta_m - t0: WRF reads it from wrfinput (Registry dataname THM)
  const std::vector<double> thm = read_var_if_present(wi_path, "THM", qv.size());
  if (thm.empty()) throw std::runtime_error("real: " + wi_path + " lacks THM (needs real.exe of WRF v4)");
  auto idx3 = [&](int i, int j, int k, int ni, int nj) { return (std::size_t(k) * nj + j) * ni + i; };
  s.u_2.fill_host([&](int i, int j, int k) { return d.u[idx3(i, j, k, nxs, ny)]; });
  s.v_2.fill_host([&](int i, int j, int k) { return d.v[idx3(i, j, k, nx, nys)]; });
  s.w_2.fill_host([](int, int, int) { return 0.0; });
  s.t_2.fill_host([&](int i, int j, int k) { return thm[idx3(i, j, k, nx, ny)]; });
  s.ph_2.fill_host([&](int i, int j, int k) { return d.ph[idx3(i, j, k, nx, ny)]; });
  s.mu_2.fill_host([&](int i, int j, int) { return d.mu[std::size_t(j) * nx + i]; });
  s.q_2[s.P_QV].fill_host([&](int i, int j, int k) { return qv[idx3(i, j, k, nx, ny)]; });
  for (int q = 0; q < s.n_scalars(); ++q) {  // other moist/scalar species from wrfinput if present, else 0
    if (q == s.P_QV) continue;
    const std::string& nm = s.q_names[q];
    std::vector<double> extra_q;
    if (nm == "qni") extra_q = read_var_if_present(wi_path, "QNICE", qv.size());
    else if (nm == "qnr") extra_q = read_var_if_present(wi_path, "QNRAIN", qv.size());
    else if (nm == "tke") extra_q = read_var_if_present(wi_path, "TKE", qv.size());
    const std::vector<double>* src = nm == "qc" ? &d.qcloud : nm == "qr" ? &d.qrain : nm == "qi" ? &d.qice
                                     : nm == "qs" ? &d.qsnow : nm == "qg" ? &d.qgraup
                                     : (!extra_q.empty() ? &extra_q : nullptr);
    const bool have = src && src->size() == qv.size();
    s.q_2[q].fill_host([&](int i, int j, int k) { return have ? (*src)[idx3(i, j, k, nx, ny)] : 0.0; });
    fill_halo(s.q_2[q], g);
    s.q_1[q].copy_from(s.q_2[q]);
  }
  // al (hypsometric_opt) and p from start_em, theta_m with qvf = 1
  std::vector<double> al(pb.size()), p(pb.size());
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        const std::size_t n = idx3(i, j, k, nx, ny), n1 = idx3(i, j, k + 1, nx, ny);
        const double mu = d.mu[std::size_t(j) * nx + i], mt = mub[std::size_t(j) * nx + i] + mu;
        double a;
        if (hypso == 2) {
          const double pfu = d.c3f[k + 1] * mt + d.c4f[k + 1] + d.p_top, pfd = d.c3f[k] * mt + d.c4f[k] + d.p_top;
          const double phm = d.c3h[k] * mt + d.c4h[k] + d.p_top;
          a = (d.ph[n1] - d.ph[n] + phb[n1] - phb[n]) / phm / std::log(pfd / pfu) - alb[n];
        } else {
          // start_em.F:812: (c1h mub + c2h) + c1h mu_1, not c1h (mub + mu_1) + c2h
          const double mb = mub[std::size_t(j) * nx + i];
          a = -1.0 / ((d.c1h[k] * mb + d.c2h[k]) + d.c1h[k] * mu) * (alb[n] * (d.c1h[k] * mu) + d.rdnw[k] * (d.ph[n1] - d.ph[n]));
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

  // input_wrf: "DX and DY do not match comparing namelist to the input file" (fatal; a nest's wrfinput too,
  // against check_a_mundo's dx(1) / accumulated parent_grid_ratio, which src/cases/nest_driver.cpp passes)
  for (const char* k : {"real.dx", "real.dy"})
    if (cfg.has(k)) {
      const double file = k[5] == 'x' ? d.dx : d.dy;
      if (std::fabs(file - cfg.get_double(k)) > 1.e-5 * file)
        throw std::runtime_error("DX and DY do not match comparing namelist to the input file");
    }

  // ---- nested domain (WRF config_flags%nested): no wrfbdy; the boundary slabs are written by the
  // nest forcing from the parent each parent step (src/cases/nest*), dtbc reset there.
  if (cfg.get_bool("real.nested", false)) {
    m.bc = std::make_unique<dyn::LateralBC>(g, g.dt, cfg.get_int("real.spec_bdy_width", 5), 1, 4, true);
    return;
  }

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
  // WRF's dtbc is a running sum of dt, reset only when the model time equals a boundary time exactly
  // (mediation_integrate.F); dtbc = t0 - need * interval below is that sum only when dt is an integer
  // (exact sums) that divides the interval
  if (g.dt != std::floor(g.dt) || std::fmod(interval, g.dt) != 0.0)
    throw std::runtime_error("real: time.dt must be an integer number of seconds that divides the boundary interval");
  // record r is used for [start + r * interval, start + (r + 1) * interval): WRF selects records by their times
  // (mediation_integrate.F), so the file must start at the simulation start with a uniform interval
  if (!d.start_date.empty())
    for (std::size_t r = 0; r < bd->this_bdy_times.size(); ++r)
      if (bd->this_bdy_times[r].substr(0, 19) != wrf_io::add_seconds(d.start_date, std::lround(double(r) * interval)).substr(0, 19))
        throw std::runtime_error("real: wrfbdy record " + std::to_string(r) + " starts at " + bd->this_bdy_times[r] +
                                 ", not at the simulation start " + d.start_date + " + " + std::to_string(r) +
                                 " boundary intervals");
  dyn::LateralBC* bc = m.bc.get();
  const int p_qv = m.st->P_QV;
  auto load = [bd, bc, nx, ny, nz, nxs, nys, p_qv](int rec) {
    using LB = dyn::LateralBC;
    if (rec >= bd->n_times) throw std::runtime_error("real: run longer than the boundary file");
    auto put = [&](LB::Var v, int npx, int npy, int nk, const std::vector<double>* const (&a)[8]) {
      const int npos[4] = {npx, npx, npy, npy};
      for (int sd = 0; sd < 4; ++sd) {
        const std::size_t n = std::size_t(npos[sd]) * nk * bd->bdy_width;
        auto cut = [&](const std::vector<double>& all) {
          if (std::size_t((rec + 1) * n) > all.size()) {
            throw std::runtime_error("real: boundary record " + std::to_string(rec) + " exceeds buffer size (" +
                                     std::to_string(all.size()) + ")");
          }
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
      for (int sd = 0; sd < 4; ++sd) {  // QVAPOR -> scalar P_QV
        const std::size_t n = std::size_t(sd < 2 ? ny : nx) * nz * bd->bdy_width;
        auto cut = [&](const std::vector<double>& all) {
          if (std::size_t((rec + 1) * n) > all.size()) throw std::runtime_error("real: QVAPOR boundary record too short");
          return std::vector<double>(all.begin() + std::ptrdiff_t(rec * n), all.begin() + std::ptrdiff_t((rec + 1) * n));
        };
        bc->set_scalar(p_qv, LB::Side(sd), cut(*a[sd]), cut(*a[4 + sd]));
      } }
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
