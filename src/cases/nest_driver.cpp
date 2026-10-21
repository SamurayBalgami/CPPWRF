#include "cases/nest_driver.hpp"

#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <sstream>
#include <stdexcept>

#include "core/backend.hpp"
#include "core/constants.hpp"
#include "core/halo.hpp"
#include "nest/nest_bdy.hpp"
#include "nest/nest_couple.hpp"
#include "nest/nest_feedback.hpp"
#include "nest/nest_init.hpp"
#include "nest/nest_interp.hpp"

#ifdef CPPWRF_HAVE_NETCDF
#include "io_netcdf/wrf_readers.hpp"
#include <netcdf.h>
#endif

namespace cppwrf::cases {

namespace {

std::string make_dom_tag(int dom_id) {
  char buf[16];
  std::snprintf(buf, sizeof(buf), "d%02d", dom_id);
  return std::string(buf);
}

#ifdef CPPWRF_HAVE_NETCDF
struct BaseParams {
  double t00, p00, a, tiso, a_strat, p_strat;
};

BaseParams read_base_params_nc(const std::string& path) {
  int nc = 0;
  if (nc_open(path.c_str(), NC_NOWRITE, &nc) != NC_NOERR)
    throw std::runtime_error("nest_driver: cannot open " + path);
  auto get = [&](const char* name) {
    int v = 0;
    double x = 0;
    if (nc_inq_varid(nc, name, &v) != NC_NOERR || nc_get_var_double(nc, v, &x) != NC_NOERR) {
      nc_close(nc);
      throw std::runtime_error(std::string("nest_driver: wrfinput lacks base-state parameter ") + name);
    }
    return x;
  };
  BaseParams b{get("T00"), get("P00"), get("TLP"), get("TISO"), get("TLP_STRAT"), get("P_STRAT")};
  nc_close(nc);
  return b;
}

// Re-run WRF start_em base-state, hydrostatic al/p, optional press_adj on mu_2, halo fills,
// and surface w initialization after terrain blending (child, press_adj=true) or feedback (parent, press_adj=false).
void recompute_start_domain(const io::Config& cfg, Model& m, const BaseField2* ht_fine_for_press_adj) {
  const std::string dir = cfg.get_string("real.dir");
  const std::string wi = dir + "/" + cfg.get_string("real.wrfinput", "wrfinput_d01");
  const BaseParams bp = read_base_params_nc(wi);
  const Grid& g = m.g;
  const int nx = g.nx, ny = g.ny, nz = g.nz;
  const int hypso = m.ap.hypsometric_opt;
  dyn::BaseState& b = *m.base;
  dyn::State& s = *m.st;

  std::vector<real> h_ht;
  b.ht.download(h_ht);
  const Layout& Lht = b.ht.layout();

  std::vector<double> ht_2d(std::size_t(nx) * ny);
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      ht_2d[std::size_t(j) * nx + i] = static_cast<double>(h_ht[Lht.index(i, j, 0)]);

  std::vector<double> pb(std::size_t(nx) * ny * nz), alb(pb.size()), t_init(pb.size());
  std::vector<double> mub(std::size_t(nx) * ny), phb(std::size_t(nx) * ny * (nz + 1));
  {
    const double t00 = bp.t00, p00 = bp.p00, a = bp.a, tiso = bp.tiso, a_strat = bp.a_strat, p_strat = bp.p_strat;
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        const std::size_t ij = std::size_t(j) * nx + i;
        const double ht = ht_2d[ij];
        const double p_surf =
            p00 * std::exp(-t00 / a + std::pow((t00 / a) * (t00 / a) - 2. * cst::g * ht / a / cst::r_d, 0.5));  // **0.5
        for (int k = 0; k < nz; ++k) {
          const std::size_t n = (std::size_t(k) * ny + j) * nx + i;
          pb[n] = double(b.c3h[k]) * (p_surf - b.p_top) + double(b.c4h[k]) + b.p_top;
          double temp = std::max(tiso, t00 + a * std::log(pb[n] / p00));
          if (pb[n] < p_strat) temp = tiso + a_strat * std::log(pb[n] / p_strat);
          t_init[n] = temp * std::pow(p00 / pb[n], cst::r_d / cst::cp) - cst::t0;
          alb[n] = (cst::r_d / cst::p1000mb) * (t_init[n] + cst::t0) * std::pow(pb[n] / cst::p1000mb, cst::cvpm);
        }
        mub[ij] = p_surf - b.p_top;
        phb[ij] = ht * cst::g;
        for (int k = 1; k <= nz; ++k) {
          const std::size_t n = (std::size_t(k) * ny + j) * nx + i;
          const std::size_t nm = (std::size_t(k - 1) * ny + j) * nx + i;
          if (hypso == 1) {
            phb[n] = phb[nm] - double(b.dnw[k - 1]) * (double(b.c1h[k - 1]) * mub[ij] + double(b.c2h[k - 1])) * alb[nm];
          } else {
            const double pfu = double(b.c3f[k]) * mub[ij] + double(b.c4f[k]) + b.p_top;
            const double pfd = double(b.c3f[k - 1]) * mub[ij] + double(b.c4f[k - 1]) + b.p_top;
            const double phm = double(b.c3h[k - 1]) * mub[ij] + double(b.c4h[k - 1]) + b.p_top;
            phb[n] = phb[nm] + alb[nm] * phm * std::log(pfd / pfu);
          }
        }
        for (int k = 0; k < nz; ++k) {
          const std::size_t n = (std::size_t(k) * ny + j) * nx + i;
          pb[n] = double(b.c3h[k]) * mub[ij] + double(b.c4h[k]) + b.p_top;
          alb[n] = (cst::r_d / cst::p1000mb) * (t_init[n] + cst::t0) * std::pow(pb[n] / cst::p1000mb, cst::cvpm);
        }
      }
  }

  auto cl = [](int a, int lo, int hi) { return a < lo ? lo : (a > hi ? hi : a); };
  auto m2 = [&](const std::vector<double>& v) {
    return [&v, nx, ny, cl](int i, int j, int) { return v[std::size_t(cl(j, 0, ny - 1)) * nx + cl(i, 0, nx - 1)]; };
  };
  auto m3 = [&](const std::vector<double>& v, int nk) {
    return [&v, nx, ny, nk, cl](int i, int j, int k) {
      return v[(std::size_t(cl(k, 0, nk - 1)) * ny + cl(j, 0, ny - 1)) * nx + cl(i, 0, nx - 1)];
    };
  };

  b.ht.fill_host_halo(m2(ht_2d));
  b.mub.fill_host_halo(m2(mub));
  b.pb.fill_host_halo(m3(pb, nz));
  b.alb.fill_host_halo(m3(alb, nz));
  b.t_init.fill_host_halo(m3(t_init, nz));
  b.phb.fill_host_halo(m3(phb, nz + 1));

  std::vector<real> h_ph, h_mu, h_thm;
  s.ph_2.download(h_ph);
  s.mu_2.download(h_mu);
  s.t_2.download(h_thm);
  const Layout& Lph = s.ph_2.layout();
  const Layout& Lmu = s.mu_2.layout();
  const Layout& Lm = s.t_2.layout();

  auto idx3 = [&](int i, int j, int k) { return (std::size_t(k) * ny + j) * nx + i; };
  std::vector<double> al(pb.size()), p(pb.size());
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        const std::size_t n = idx3(i, j, k), n1 = idx3(i, j, k + 1);
        const double mu = static_cast<double>(h_mu[Lmu.index(i, j, 0)]);
        const double mt = mub[std::size_t(j) * nx + i] + mu;
        const double ph0 = static_cast<double>(h_ph[Lph.index(i, j, k)]);
        const double ph1 = static_cast<double>(h_ph[Lph.index(i, j, k + 1)]);
        const double thm = static_cast<double>(h_thm[Lm.index(i, j, k)]);
        double a;
        if (hypso == 2) {
          const double pfu = double(b.c3f[k + 1]) * mt + double(b.c4f[k + 1]) + b.p_top;
          const double pfd = double(b.c3f[k]) * mt + double(b.c4f[k]) + b.p_top;
          const double phm = double(b.c3h[k]) * mt + double(b.c4h[k]) + b.p_top;
          a = (ph1 - ph0 + phb[n1] - phb[n]) / phm / std::log(pfd / pfu) - alb[n];
        } else {
          a = -1.0 / (double(b.c1h[k]) * mt + double(b.c2h[k])) *
              (alb[n] * (double(b.c1h[k]) * mu) + double(b.rdnw[k]) * (ph1 - ph0));
        }
        al[n] = a;
        p[n] = cst::p1000mb * std::pow((cst::r_d * (cst::t0 + thm)) / (cst::p1000mb * (a + alb[n])), cst::cpovcv) -
               pb[n];
      }

  s.al.fill_host([&](int i, int j, int k) { return al[idx3(i, j, k)]; });
  s.p.fill_host([&](int i, int j, int k) { return p[idx3(i, j, k)]; });

  // start_em.F lines 988-996: IF (press_adj) mu_2(i,j) += al(i,1,j) / (alt(i,1,j) * alb(i,1,j)) * g * (ht - ht_fine)
  if (ht_fine_for_press_adj != nullptr) {
    std::vector<real> h_ht_fine;
    ht_fine_for_press_adj->download(h_ht_fine);
    const Layout& Lhf = ht_fine_for_press_adj->layout();
    std::vector<double> mu_adj(std::size_t(nx) * ny);
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        const std::size_t ij = std::size_t(j) * nx + i;
        const std::size_t n0 = idx3(i, j, 0);
        const double mu = static_cast<double>(h_mu[Lmu.index(i, j, 0)]);
        const double htf = static_cast<double>(h_ht_fine[Lhf.index(i, j, 0)]);
        const double alt0 = al[n0] + alb[n0];
        mu_adj[ij] = mu + al[n0] / (alt0 * alb[n0]) * cst::g * (ht_2d[ij] - htf);
      }
    s.mu_2.fill_host([&](int i, int j, int) { return mu_adj[std::size_t(j) * nx + i]; });
  }

  fill_halos({&s.u_2, &s.v_2, &s.w_2, &s.t_2, &s.ph_2, &s.al, &s.p}, g);
  fill_halo(s.mu_2, g);
  for (int q = 0; q < s.n_scalars(); ++q) {
    fill_halo(s.q_2[q], g);
    s.q_1[q].copy_from(s.q_2[q]);
  }
  s.u_1.copy_from(s.u_2);
  s.v_1.copy_from(s.v_2);
  s.t_1.copy_from(s.t_2);
  s.ph_1.copy_from(s.ph_2);
  s.mu_1.copy_from(s.mu_2);
  dyn::set_w_surface(s.w_1, s.u_1, s.v_1, b, g, true);
  dyn::set_w_surface(s.w_2, s.u_2, s.v_2, b, g, true);
}
#endif

}  // namespace

bool has_nests(const io::Config& root_cfg) {
  return root_cfg.get_int("domains.max_dom", 1) > 1;
}

std::vector<NestDomainSpec> parse_nest_specs(const io::Config& root_cfg) {
  const int max_dom = static_cast<int>(root_cfg.get_int("domains.max_dom", 1));
  if (max_dom < 1) throw std::runtime_error("domains.max_dom must be >= 1");

  // Mark all [domains] and [d0N] keys accessed on root_cfg so DomainRunner::start() doesn't warn on them
  for (const std::string& k : root_cfg.keys()) {
    if (k.rfind("domains.", 0) == 0 || (k.size() >= 4 && k[0] == 'd' && std::isdigit(k[1]) && std::isdigit(k[2]) && k[3] == '.')) {
      root_cfg.mark_accessed(k);
    }
  }

  const int def_feedback = static_cast<int>(root_cfg.get_int("domains.feedback", 1));
  const int def_smooth = static_cast<int>(root_cfg.get_int("domains.smooth_option", 2));  // WRF Registry default
  const int def_save_topo = static_cast<int>(root_cfg.get_int("domains.save_topo_from_real", 1));
  const int def_spec_bdy = static_cast<int>(root_cfg.get_int("domains.spec_bdy_width", 5));
  const int def_blend = static_cast<int>(root_cfg.get_int("domains.blend_width", 5));

  std::vector<NestDomainSpec> specs(max_dom);

  for (int id = 1; id <= max_dom; ++id) {
    NestDomainSpec& sp = specs[id - 1];
    sp.dom_id = id;
    sp.tag = make_dom_tag(id);
    const std::string pfx = sp.tag + ".";
    const std::string dpfx = "domains." + sp.tag + "_";

    auto get_int_dom = [&](const std::string& name, int def_val) -> int {
      if (root_cfg.has(pfx + name)) return static_cast<int>(root_cfg.get_int(pfx + name));
      if (root_cfg.has(dpfx + name)) return static_cast<int>(root_cfg.get_int(dpfx + name));
      if (id == 2 && root_cfg.has("domains." + name)) return static_cast<int>(root_cfg.get_int("domains." + name));
      return def_val;
    };

    if (id == 1) {
      sp.parent_id = 0;
      sp.i_parent_start = 1;
      sp.j_parent_start = 1;
      sp.parent_grid_ratio = 1;
      sp.parent_time_step_ratio = 1;
      sp.feedback = def_feedback;
      sp.smooth_option = def_smooth;
    } else {
      sp.parent_id = get_int_dom("parent_id", id - 1);
      if (sp.parent_id < 1 || sp.parent_id >= id)
        throw std::runtime_error("invalid parent_id for domain " + sp.tag);
      sp.i_parent_start = get_int_dom("i_parent_start", 1);
      sp.j_parent_start = get_int_dom("j_parent_start", 1);
      sp.parent_grid_ratio = get_int_dom("parent_grid_ratio", 3);
      sp.parent_time_step_ratio = get_int_dom("parent_time_step_ratio", 1);
      sp.feedback = get_int_dom("feedback", def_feedback);
      sp.smooth_option = get_int_dom("smooth_option", def_smooth);
      sp.save_topo_from_real = get_int_dom("save_topo_from_real", def_save_topo);
      sp.spec_bdy_width = get_int_dom("spec_bdy_width", def_spec_bdy);
      sp.blend_width = get_int_dom("blend_width", def_blend);
    }

    // Build per-domain io::Config
    std::string ext_cfg_path;
    if (root_cfg.has(pfx + "config")) ext_cfg_path = root_cfg.get_string(pfx + "config");
    else if (root_cfg.has(dpfx + "config")) ext_cfg_path = root_cfg.get_string(dpfx + "config");

    std::vector<std::string> ov;
    if (id >= 2) {
      const NestDomainSpec& parent_sp = specs[sp.parent_id - 1];
      const double parent_dt = parent_sp.cfg.get_double("time.dt", 60.0);
      const double child_dt = parent_dt / sp.parent_time_step_ratio;
      const long long parent_out_iv = parent_sp.cfg.get_int("output.out_interval", 1);
      const long long parent_diag_iv = parent_sp.cfg.get_int("output.diag_interval", 10);

      ov.push_back("real.wrfinput=wrfinput_" + sp.tag);
      ov.push_back("real.nested=true");
      ov.push_back("real.spec_bdy_width=" + std::to_string(sp.spec_bdy_width));
      {
        std::ostringstream oss;
        oss.precision(17);
        oss << "time.dt=" << child_dt;
        ov.push_back(oss.str());
      }
      ov.push_back("output.out_interval=" + std::to_string(parent_out_iv * sp.parent_time_step_ratio));
      ov.push_back("output.diag_interval=" + std::to_string(parent_diag_iv * sp.parent_time_step_ratio));
      // a resumed run: this domain's restart file next to the root's (wrfrst_d01_<date> -> wrfrst_d0N_<date>)
      const std::string rf = root_cfg.get_string("run.restart_file", "");
      const std::string key = "wrfrst_d01_";
      const std::size_t at = rf.rfind(key);
      if (!rf.empty()) {
        if (at == std::string::npos)
          throw std::runtime_error("nested runs resume from WRF restart files (run.restart_file = wrfrst_d01_<date>)");
        ov.push_back("run.restart_file=" + rf.substr(0, at) + "wrfrst_" + sp.tag + "_" + rf.substr(at + key.size()));
      }
    }

    // Map shorthand keys in [d0N] or domains.d0N_* to standard Model/DomainRunner sections
    auto map_key = [&](const std::string& sub, const std::string& val) {
      if (sub == "config" || sub == "parent_id" || sub == "i_parent_start" || sub == "j_parent_start" ||
          sub == "parent_grid_ratio" || sub == "parent_time_step_ratio" || sub == "feedback" ||
          sub == "smooth_option" || sub == "save_topo_from_real" || sub == "blend_width") {
        return;
      }
      if (sub == "spec_bdy_width") ov.push_back("real.spec_bdy_width=" + val);
      else if (sub == "wrfinput" || sub == "wrfbdy" || sub == "dir" || sub == "bdy_interval" || sub == "nested")
        ov.push_back("real." + sub + "=" + val);
      else if (sub == "dt" || sub == "time_step_sound")
        ov.push_back("time." + sub + "=" + val);
      else if (sub == "out_interval" || sub == "diag_interval" || sub == "format" || sub == "out_dir" || sub == "restart_interval")
        ov.push_back("output." + (sub == "out_dir" ? std::string("dir") : sub) + "=" + val);
      else if (sub == "nx" || sub == "ny" || sub == "nz" || sub == "dx" || sub == "dy")
        ov.push_back("domain." + sub + "=" + val);
      else if (sub.rfind("phys_", 0) == 0)  // per-domain physics (WRF max_domains options: cu_physics, radt, ...)
        ov.push_back("phys." + sub.substr(5) + "=" + val);
      else if (sub.find('.') != std::string::npos)
        ov.push_back(sub + "=" + val);
      else
        ov.push_back("dyn." + sub + "=" + val);
    };

    for (const std::string& k : root_cfg.keys()) {
      if (k.rfind(pfx, 0) == 0) {
        map_key(k.substr(pfx.size()), root_cfg.get_string(k));
      } else if (k.rfind(dpfx, 0) == 0) {
        map_key(k.substr(dpfx.size()), root_cfg.get_string(k));
      }
    }

    if (!ext_cfg_path.empty()) {
      sp.cfg = io::Config::load(ext_cfg_path, ov);
    } else {
      // Serialize non-domain-specific keys from root_cfg as base INI
      std::ostringstream base_ini;
      for (const std::string& k : root_cfg.keys()) {
        if (k.rfind("domains.", 0) == 0) continue;
        if (k.size() >= 4 && k[0] == 'd' && std::isdigit(k[1]) && std::isdigit(k[2]) && k[3] == '.') continue;
        if (id >= 2 && (k == "real.wrfbdy" || k == "real.bdy_interval")) continue;
        auto dot = k.find('.');
        if (dot != std::string::npos) {
          base_ini << "[" << k.substr(0, dot) << "]\n" << k.substr(dot + 1) << " = " << root_cfg.get_string(k) << "\n";
        }
      }
      sp.cfg = io::Config::from_string(base_ini.str(), ov);
    }
    if (sp.cfg.has("output.dir")) sp.cfg.mark_accessed("output.dir");
  }

  return specs;
}

void med_nest_initial(const io::Config& parent_cfg, Model& parent,
                      const NestDomainSpec& child_spec, Model& child) {
#ifndef CPPWRF_HAVE_NETCDF
  (void)parent_cfg; (void)parent; (void)child_spec; (void)child;
  throw std::runtime_error("med_nest_initial requires NetCDF support");
#else
  const Grid& cg = child.g;
  const int nx = cg.nx, ny = cg.ny, nz = cg.nz;
  const int ipos = child_spec.i_parent_start;
  const int jpos = child_spec.j_parent_start;
  const int nri = child_spec.parent_grid_ratio;
  const int nrj = child_spec.parent_grid_ratio;

  // 1. med_interp_domain(parent, child): SINT-interpolate parent ht, mub, phb onto child
  BaseField2 ht_int(cg, Stagger::M, "ht_int", false), mub_int(cg, Stagger::M, "mub_int", false);
  BaseField3 phb_int(cg, Stagger::W, "phb_int", false);
  {
    auto c_ht = nest::const_view2(parent.base->ht);
    auto n_ht = nest::mut_view2(ht_int);
    parallel_for(Range3D{0, nx, 0, ny, 0, 1},
                 nest::KSintView2D<decltype(c_ht), decltype(n_ht)>{c_ht, n_ht, ipos, jpos, nri, nrj, false, false});

    auto c_mub = nest::const_view2(parent.base->mub);
    auto n_mub = nest::mut_view2(mub_int);
    parallel_for(Range3D{0, nx, 0, ny, 0, 1},
                 nest::KSintView2D<decltype(c_mub), decltype(n_mub)>{c_mub, n_mub, ipos, jpos, nri, nrj, false, false});

    auto c_phb = nest::const_view3(parent.base->phb);
    auto n_phb = nest::mut_view3(phb_int);
    parallel_for(Range3D{0, nx, 0, ny, 0, nz + 1},
                 nest::KSintView3D<decltype(c_phb), decltype(n_phb)>{c_phb, n_phb, ipos, jpos, nri, nrj, false, false});
  }

  // 2. Load raw wrfinput_d0N fields (before start_domain(child) recomputes base state / w)
  const std::string dir = child_spec.cfg.get_string("real.dir");
  const std::string wi = dir + "/" + child_spec.cfg.get_string("real.wrfinput", "wrfinput_" + child_spec.tag);
  const io_netcdf::WrfInputData d = io_netcdf::read_wrfinput(wi);

  auto cl = [](int a, int lo, int hi) { return a < lo ? lo : (a > hi ? hi : a); };
  auto m2 = [&](const std::vector<double>& v) {
    return [&v, nx, ny, cl](int i, int j, int) { return v[std::size_t(cl(j, 0, ny - 1)) * nx + cl(i, 0, nx - 1)]; };
  };
  auto m3 = [&](const std::vector<double>& v, int nk) {
    return [&v, nx, ny, nk, cl](int i, int j, int k) {
      return v[(std::size_t(cl(k, 0, nk - 1)) * ny + cl(j, 0, ny - 1)) * nx + cl(i, 0, nx - 1)];
    };
  };

  BaseField2 ht_fine(cg, Stagger::M, "ht_fine", false), save_mub(cg, Stagger::M, "save_mub", false);
  ht_fine.fill_host_halo(m2(d.hgt));
  save_mub.fill_host_halo(m2(d.mub));

  child.base->ht.fill_host_halo(m2(d.hgt));
  child.base->mub.fill_host_halo(m2(d.mub));
  child.base->phb.fill_host_halo(m3(d.phb, nz + 1));
  child.base->pb.fill_host_halo(m3(d.pb, nz));
  child.base->alb.fill_host_halo(m3(d.alb, nz));

  // p: kept from the child's own start (WRF: med_initialdata_input runs start_domain(nest) after reading
  // wrfinput, which recomputes p from theta_m and al before adjust_tempqv)
  child.st->w_2.fill_host([](int, int, int) { return 0.0; });
  child.st->w_1.fill_host([](int, int, int) { return 0.0; });

  // 3. blend_terrain(ht, mub, phb) when save_topo_from_real == 1
  if (child_spec.save_topo_from_real == 1) {
    nest::blend_terrain_field2(ht_int, child.base->ht, nx, ny, child_spec.spec_bdy_width, child_spec.blend_width);
    nest::blend_terrain_field2(mub_int, child.base->mub, nx, ny, child_spec.spec_bdy_width, child_spec.blend_width);
    nest::blend_terrain_field3(phb_int, child.base->phb, nx, ny, nz + 1, child_spec.spec_bdy_width,
                               child_spec.blend_width);
  }

  // 4. adjust_tempqv on child t_2 and qv
  const int use_theta_m = child.st->P_QV >= 0 ? 1 : 0;
  if (child.st->P_QV >= 0) {
    nest::adjust_tempqv_fields(child.base->mub, save_mub, child.base->c3h, child.base->c4h, child.base->p_top,
                               child.st->t_2, child.st->p, child.st->q_2[child.st->P_QV], use_theta_m, nx, ny, nz);
  }

  // 5. If feedback != 0, med_nest_feedback(parent, child) runs BEFORE start_domain(child) and start_domain(parent)
  if (child_spec.feedback != 0) {
    med_nest_feedback(parent, child, child_spec);
  }

  // 6. start_domain(child) with press_adj = true
  recompute_start_domain(child_spec.cfg, child, &ht_fine);

  // 7. start_domain(parent) with press_adj = false
  recompute_start_domain(parent_cfg, parent, nullptr);
  dev_sync();
#endif
}

void med_nest_force(Model& parent, Model& child, const NestDomainSpec& child_spec) {
  const Grid& pg = parent.g;
  const Grid& cg = child.g;
  const int ipos = child_spec.i_parent_start;
  const int jpos = child_spec.j_parent_start;
  const int nri = child_spec.parent_grid_ratio;
  const int nrj = child_spec.parent_grid_ratio;
  const double cdt = pg.dt;

  // 1. Couple parent in-place
  nest::couple_or_uncouple_state(pg, *parent.st, *parent.base, true);
  fill_halos({&parent.st->u_2, &parent.st->v_2, &parent.st->w_2, &parent.st->ph_2, &parent.st->t_2}, pg);
  for (int s = 0; s < parent.st->n_scalars(); ++s) {
    fill_halo(parent.st->q_2[s], pg);
  }
  // mu_2's halo is left as the dycore's halo schedule had it (WRF's uncouple reads it at ids/jds)

  // 2. Couple child in-place
  nest::couple_or_uncouple_state(cg, *child.st, *child.base, true);

  // 3. Interpolate parent -> child LateralBC slabs (bdy_interp1 / KBdySlabView3D / KBdySlabView2D)
  dyn::LateralBC& lbc = *child.bc;
  nest::force_nest_bdy_field3(lbc, dyn::LateralBC::U, parent.st->u_2, child.st->u_2, cg.nx, cg.ny, ipos, jpos, nri, nrj,
                              cdt);
  nest::force_nest_bdy_field3(lbc, dyn::LateralBC::V, parent.st->v_2, child.st->v_2, cg.nx, cg.ny, ipos, jpos, nri, nrj,
                              cdt);
  nest::force_nest_bdy_field3(lbc, dyn::LateralBC::W, parent.st->w_2, child.st->w_2, cg.nx, cg.ny, ipos, jpos, nri, nrj,
                              cdt);
  nest::force_nest_bdy_field3(lbc, dyn::LateralBC::PH, parent.st->ph_2, child.st->ph_2, cg.nx, cg.ny, ipos, jpos, nri,
                              nrj, cdt);
  nest::force_nest_bdy_field3(lbc, dyn::LateralBC::T, parent.st->t_2, child.st->t_2, cg.nx, cg.ny, ipos, jpos, nri, nrj,
                              cdt);
  nest::force_nest_bdy_field2(lbc, dyn::LateralBC::MU, parent.st->mu_2, child.st->mu_2, cg.nx, cg.ny, ipos, jpos, nri,
                              nrj, cdt);

  const int n_scalars = std::min(parent.st->n_scalars(), child.st->n_scalars());
  for (int s = 0; s < n_scalars; ++s) {
    // Registry: tke has no f (force) flag - the nest keeps flow-dependent boundaries for it (solve_em)
    if (s == child.st->P_TKE) continue;
    nest::force_nest_bdy_scalar_field3(lbc, s, parent.st->q_2[s], child.st->q_2[s], cg.nx, cg.ny, cg.nz, ipos, jpos,
                                       nri, nrj, cdt);
  }
  lbc.dtbc = 0.0;

  // o3rad (Registry rdf=(p2c)): the child's whole field from the parent's, which only domain 1 computes
  if (parent.suite && child.suite && parent.suite->radiation() && child.suite->radiation()) {
    fill_halo(parent.suite->o3rad, pg);
    auto c_o3 = nest::const_view3(parent.suite->o3rad);
    auto n_o3 = nest::mut_view3(child.suite->o3rad);
    parallel_for(Range3D{0, cg.nx, 0, cg.ny, 0, cg.nz},
                 nest::KSintView3D<decltype(c_o3), decltype(n_o3)>{c_o3, n_o3, ipos, jpos, nri, nrj, false, false});
  }

  // 4. Uncouple child and parent in-place (preserving exact reciprocal-multiply roundoff)
  nest::couple_or_uncouple_state(cg, *child.st, *child.base, false);
  nest::couple_or_uncouple_state(pg, *parent.st, *parent.base, false);

  fill_halos({&child.st->u_2, &child.st->v_2, &child.st->w_2, &child.st->ph_2, &child.st->t_2}, cg);
  for (int s = 0; s < child.st->n_scalars(); ++s) {
    fill_halo(child.st->q_2[s], cg);
  }
  fill_halos({&parent.st->u_2, &parent.st->v_2, &parent.st->w_2, &parent.st->ph_2, &parent.st->t_2}, pg);
  for (int s = 0; s < parent.st->n_scalars(); ++s) {
    fill_halo(parent.st->q_2[s], pg);
  }
}

void med_nest_feedback(Model& parent, const Model& child, const NestDomainSpec& child_spec) {
  const Grid& pg = parent.g;
  const Grid& cg = child.g;
  const int n_nx = cg.nx, n_ny = cg.ny, n_nz = cg.nz;
  const int ipos = child_spec.i_parent_start;
  const int jpos = child_spec.j_parent_start;
  const int nri = child_spec.parent_grid_ratio;
  const int nrj = child_spec.parent_grid_ratio;
  const int sm = child_spec.smooth_option;

  // Note: WRF med_feedback_domain operates directly on UNCOUPLED fields (no couple_or_uncouple_em call)
  nest::feedback_field3(pg, parent.st->u_2, child.st->u_2, n_nx, n_ny, n_nz, ipos, jpos, nri, nrj, true, false, sm);
  nest::feedback_field3(pg, parent.st->v_2, child.st->v_2, n_nx, n_ny, n_nz, ipos, jpos, nri, nrj, false, true, sm);
  nest::feedback_field3(pg, parent.st->w_2, child.st->w_2, n_nx, n_ny, n_nz + 1, ipos, jpos, nri, nrj, false, false, sm);
  nest::feedback_field3(pg, parent.st->ph_2, child.st->ph_2, n_nx, n_ny, n_nz + 1, ipos, jpos, nri, nrj, false, false, sm);
  nest::feedback_field3(pg, parent.st->t_2, child.st->t_2, n_nx, n_ny, n_nz, ipos, jpos, nri, nrj, false, false, sm);
  nest::feedback_field2(pg, parent.st->mu_2, child.st->mu_2, n_nx, n_ny, ipos, jpos, nri, nrj, sm);

  const int n_scalars = std::min(parent.st->n_scalars(), child.st->n_scalars());
  for (int s = 0; s < n_scalars; ++s) {
    if (s == parent.st->P_TKE) continue;  // tke is not fed back (Registry: no u flag)
    nest::feedback_field3(pg, parent.st->q_2[s], child.st->q_2[s], n_nx, n_ny, n_nz, ipos, jpos, nri, nrj, false, false, sm);
  }

  // Base-state fields fed back in inc/nest_feedbackup_interp.inc; inc/nest_feedbackup_smooth.inc smooths them
  // and ht as well (map factors, f, e, sina, cosa are fed back without smoothing)
  nest::feedback_field3(pg, parent.base->phb, child.base->phb, n_nx, n_ny, n_nz + 1, ipos, jpos, nri, nrj, false, false, sm);
  nest::feedback_field2(pg, parent.base->mub, child.base->mub, n_nx, n_ny, ipos, jpos, nri, nrj, sm);
  nest::feedback_field3(pg, parent.base->alb, child.base->alb, n_nx, n_ny, n_nz, ipos, jpos, nri, nrj, false, false, sm);
  nest::feedback_field3(pg, parent.base->pb, child.base->pb, n_nx, n_ny, n_nz, ipos, jpos, nri, nrj, false, false, sm);
  nest::feedback_field2(pg, parent.base->ht, child.base->ht, n_nx, n_ny, ipos, jpos, nri, nrj, sm);

  if (parent.base->mapped && child.base->mapped) {
    nest::feedback_field2_fcnm(parent.base->msftx, child.base->msftx, n_nx, n_ny, ipos, jpos, nri, nrj, false, false);
    nest::feedback_field2_fcnm(parent.base->msfty, child.base->msfty, n_nx, n_ny, ipos, jpos, nri, nrj, false, false);
    nest::feedback_field2_fcnm(parent.base->msfux, child.base->msfux, n_nx, n_ny, ipos, jpos, nri, nrj, true, false);
    nest::feedback_field2_fcnm(parent.base->msfuy, child.base->msfuy, n_nx, n_ny, ipos, jpos, nri, nrj, true, false);
    nest::feedback_field2_fcnm(parent.base->msfvx, child.base->msfvx, n_nx, n_ny, ipos, jpos, nri, nrj, false, true);
    nest::feedback_field2_fcnm(parent.base->msfvx_inv, child.base->msfvx_inv, n_nx, n_ny, ipos, jpos, nri, nrj, false, true);
    nest::feedback_field2_fcnm(parent.base->msfvy, child.base->msfvy, n_nx, n_ny, ipos, jpos, nri, nrj, false, true);
    nest::feedback_field2_fcnm(parent.base->f, child.base->f, n_nx, n_ny, ipos, jpos, nri, nrj, false, false);
    nest::feedback_field2_fcnm(parent.base->e, child.base->e, n_nx, n_ny, ipos, jpos, nri, nrj, false, false);
    nest::feedback_field2_fcnm(parent.base->sina, child.base->sina, n_nx, n_ny, ipos, jpos, nri, nrj, false, false);
    nest::feedback_field2_fcnm(parent.base->cosa, child.base->cosa, n_nx, n_ny, ipos, jpos, nri, nrj, false, false);
  }

  // update_after_feedback_em: only a DM halo exchange in WRF (no-op on one task): the parent's halos
  // keep their values (the dycore refreshes them on WRF's schedule)
}

int run_nested(const io::Config& root_cfg, const std::string& out_dir) {
  std::vector<NestDomainSpec> specs = parse_nest_specs(root_cfg);
  const int max_dom = static_cast<int>(specs.size());

  struct DomInstance {
    NestDomainSpec spec;
    Model model;
    std::unique_ptr<DomainRunner> runner;
    std::vector<int> children;  // 0-based indices into doms
    long step = 0;
    bool started = false;
  };

  std::vector<DomInstance> doms(max_dom);
  for (int i = 0; i < max_dom; ++i) {
    doms[i].spec = std::move(specs[i]);
    if (doms[i].spec.parent_id > 0) {
      doms[doms[i].spec.parent_id - 1].children.push_back(i);
    }
  }

  // 1. Build all domains and run med_nest_initial for each parent -> child link in tree order
  for (int i = 0; i < max_dom; ++i) {
    doms[i].model = Model::build(doms[i].spec.cfg);
  }
  for (int i = 0; i < max_dom; ++i) {
    for (int c_idx : doms[i].children) {
      med_nest_initial(doms[i].spec.cfg, doms[i].model, doms[c_idx].spec, doms[c_idx].model);
    }
  }

  // 2. Create DomainRunner per domain (out_dir/d01, out_dir/d02, ...)
  for (int i = 0; i < max_dom; ++i) {
    std::string dom_out = out_dir + "/" + doms[i].spec.tag;
    doms[i].runner = std::make_unique<DomainRunner>(doms[i].spec.cfg, doms[i].model, dom_out, doms[i].spec.tag);
    if (doms[i].spec.parent_id > 0) doms[i].runner->set_defer_output(true);
  }

#ifdef CPPWRF_HAVE_NETCDF
  // med_nest_initial's feedback (feedback != 0) of the history fields the model does not carry
  for (int i = 0; i < max_dom; ++i)
    for (int c_idx : doms[i].children)
      if (doms[c_idx].spec.feedback != 0) {
        doms[i].runner->set_need_history(true);
        doms[c_idx].runner->set_need_history(true);
      }
  for (int i = 0; i < max_dom; ++i) doms[i].runner->prepare_output();
  for (int i = 0; i < max_dom; ++i)
    for (int c_idx : doms[i].children)
      if (doms[c_idx].spec.feedback != 0 && doms[i].runner->history() && doms[c_idx].runner->history()) {
        const NestDomainSpec& cs = doms[c_idx].spec;
        doms[i].runner->history()->feedback_physics(*doms[c_idx].runner->history(), cs.i_parent_start,
                                                    cs.j_parent_start, cs.parent_grid_ratio);
        doms[i].runner->history()->feedback_from(*doms[c_idx].runner->history(), cs.i_parent_start,
                                                 cs.j_parent_start, cs.parent_grid_ratio);
      }
  // start_domain(parent) after the feedback (start_em, topo_wind = 0): CTOPO = CTOPO2 = 1 again (the fed-back
  // average of ones is 1 + 1 ulp)
  for (int i = 0; i < max_dom; ++i)
    if (doms[i].model.suite)
      for (int c_idx : doms[i].children)
        if (doms[c_idx].spec.feedback != 0) {
          doms[i].model.suite->ctopo.fill_host([](int, int, int) { return 1.0; });
          doms[i].model.suite->ctopo2.fill_host([](int, int, int) { return 1.0; });
        }
#endif

  // 3. Root domain (d01) writes its t=0 initial frame before its first step
  if (!doms[0].runner->start()) return 3;
  doms[0].started = true;
  doms[0].step = doms[0].runner->start_step;
  // a resumed run: the nests read their restart files before the first med_nest_force, which takes the
  // nest's boundary tendencies from its current state (no initial frames are written)
  if (doms[0].step > 0)
    for (int i = 1; i < max_dom; ++i) {
      if (!doms[i].runner->start()) return 3;
      doms[i].started = true;
      doms[i].step = doms[i].runner->start_step;
    }

  // Recursive WRF integrate(dom_idx):
  //   - One step of dom_idx
  //   - For each child c_idx:
  //       med_nest_force(dom_idx, c_idx)
  //       if (!c_idx.started) c_idx.runner->start() (writes child t=0 frame AFTER first med_nest_force, as in WRF!)
  //       for sub = 1 .. parent_time_step_ratio:
  //         advance_domain(c_idx)
  //       if (c_idx.feedback != 0) med_nest_feedback(dom_idx, c_idx)
  //   - dom_idx.runner->after_step(step) (writes parent history frame AFTER med_nest_feedback, as in WRF!)
  std::function<bool(int)> advance_domain = [&](int idx) -> bool {
    DomInstance& d = doms[idx];
    const long n = ++d.step;
    d.runner->before_step(n);
    d.model.rk->step(*d.model.st, *d.model.tend, *d.model.base);
    d.runner->after_solve(n);
    dev_sync();

    for (int c_idx : d.children) {
      DomInstance& c = doms[c_idx];
      med_nest_force(d.model, c.model, c.spec);
      dev_sync();

      if (!c.started) {
        if (!c.runner->start()) return false;
        c.started = true;
        c.step = c.runner->start_step;
      }

      const int ratio = c.spec.parent_time_step_ratio;
      for (int sub = 0; sub < ratio; ++sub) {
        if (!advance_domain(c_idx)) return false;
      }

      if (c.spec.feedback != 0) {
        med_nest_feedback(d.model, c.model, c.spec);
        dev_sync();
#ifdef CPPWRF_HAVE_NETCDF
        if (d.runner->history() && c.runner->history()) {  // the physics state, then history-only fields (Registry u=)
          d.runner->history()->feedback_physics(*c.runner->history(), c.spec.i_parent_start, c.spec.j_parent_start,
                                                c.spec.parent_grid_ratio);
          d.runner->history()->feedback_from(*c.runner->history(), c.spec.i_parent_start, c.spec.j_parent_start,
                                             c.spec.parent_grid_ratio);
        }
#endif
      }
    }

    if (!d.runner->after_step(n)) return false;
    return true;
  };

  // If run_time == 0 (0 steps), still start child runners after an initial force so t=0 snapshots exist
  if (doms[0].model.nsteps == doms[0].runner->start_step) {
    for (int i = 0; i < max_dom; ++i) {
      for (int c_idx : doms[i].children) {
        DomInstance& c = doms[c_idx];
        med_nest_force(doms[i].model, c.model, c.spec);
        dev_sync();
        if (!c.started) {
          if (!c.runner->start()) return 3;
          c.started = true;
        }
      }
    }
  } else {
    const long total_root_steps = doms[0].model.nsteps;
    while (doms[0].step < total_root_steps) {
      if (!advance_domain(0)) return 3;
    }
  }

  for (int i = 0; i < max_dom; ++i) {
    if (doms[i].started) doms[i].runner->finish();
  }

  return 0;
}

}  // namespace cppwrf::cases
