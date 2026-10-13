#include "dyn/diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "core/constants.hpp"
#include "dyn/large_step.hpp"

namespace cppwrf::dyn {
namespace {
using C3 = View3<const real>;
using C2 = View2<const real>;
using B3 = BView3<const real>;
using B2 = BView2<const real>;
constexpr double kInf = std::numeric_limits<double>::infinity();

CPPWRF_HD inline double absnan(double v) {  // NaN/Inf -> +inf so that max() exposes it
  return (v != v || v - v != 0.0) ? kInf : (v < 0 ? -v : v);
}
struct FMaxAbs { C3 f; CPPWRF_HD double operator()(int i, int j, int k) const { return absnan(double(f(i, j, k))); } };
struct FThp {  // dry theta' = theta - t_init, theta = theta_m / (1 + rvovrd qv) when moist
  C3 t; B3 ti; double sign; C3 qv; bool moist;
  CPPWRF_HD double operator()(int i, int j, int k) const {
    double th = double(t(i, j, k));
    if (moist) th = (th + cst::t0) / (1.0 + cst::rvovrd * double(qv(i, j, k))) - cst::t0;
    double v = th - double(ti(i, j, k));
    return v != v ? kInf : sign * v;
  }
};
struct FMass { C2 mu; B2 mub; CPPWRF_HD double operator()(int i, int j, int) const { return double(mu(i, j)) + double(mub(i, j)); } };
struct FThetaMass {
  C3 t; C2 mu; B2 mub; VCoef c;
  CPPWRF_HD double operator()(int i, int j, int k) const {
    double m = double(c.c1h(k)) * (double(mu(i, j)) + double(mub(i, j))) + double(c.c2h(k));
    return -double(c.dnw(k)) * m * (double(t(i, j, k)) + cst::t0);
  }
};
struct FEnergy {
  C3 u, v, w, t, p, ph; C2 mu; B2 mub; B3 pb, phb; VCoef c;
  CPPWRF_HD double operator()(int i, int j, int k) const {
    double m = double(c.c1h(k)) * (double(mu(i, j)) + double(mub(i, j))) + double(c.c2h(k));
    double uc = 0.5 * (double(u(i, j, k)) + double(u(i + 1, j, k)));
    double vc = 0.5 * (double(v(i, j, k)) + double(v(i, j + 1, k)));
    double wc = 0.5 * (double(w(i, j, k)) + double(w(i, j, k + 1)));
    double T = (double(t(i, j, k)) + cst::t0) * pow((double(p(i, j, k)) + double(pb(i, j, k))) / cst::p1000mb, cst::rcp);
    double phm = 0.5 * (double(ph(i, j, k)) + double(ph(i, j, k + 1)) + double(phb(i, j, k)) + double(phb(i, j, k + 1)));
    return -double(c.dnw(k)) * m * (0.5 * (uc * uc + vc * vc + wc * wc) + cst::cv * T + phm);
  }
};
}  // namespace

bool DiagRow::finite() const {
  for (double v : {max_abs_u, max_abs_v, max_abs_w, min_th_pert, max_th_pert, dry_mass, theta_mass, total_energy,
                   max_abs_q})
    if (!std::isfinite(v)) return false;
  return true;
}

DiagRow Diagnostics::compute(const State& s, const BaseState& b) {
  const Grid& g = g_;
  VCoef c = VCoef::from(b);
  DiagRow r{};
  r.max_abs_u = red_.max(g.cells(Stagger::U), FMaxAbs{s.u_2.cview()});
  r.max_abs_v = red_.max(g.cells(Stagger::V), FMaxAbs{s.v_2.cview()});
  r.max_abs_w = red_.max(g.cells(Stagger::W), FMaxAbs{s.w_2.cview()});
  for (int n = 0; n < s.n_scalars(); ++n)
    r.max_abs_q = std::max(r.max_abs_q, red_.max(g.mass_cells(), FMaxAbs{s.q_2[n].cview()}));
  const bool moist = s.P_QV >= 0;
  const C3 qv = moist ? s.q_2[s.P_QV].cview() : C3{};
  r.max_th_pert = red_.max(g.mass_cells(), FThp{s.t_2.cview(), b.t_init.view(), 1.0, qv, moist});
  r.min_th_pert = -red_.max(g.mass_cells(), FThp{s.t_2.cview(), b.t_init.view(), -1.0, qv, moist});
  const double area = g.dx * (g.two_d() ? 1.0 : g.dy) / cst::g;  // 2D: per metre in y
  r.dry_mass = area * red_.sum(Range3D{0, g.nx, 0, g.ny, 0, 1}, FMass{s.mu_2.cview(), b.mub.view()});
  r.theta_mass = area * red_.sum(g.mass_cells(), FThetaMass{s.t_2.cview(), s.mu_2.cview(), b.mub.view(), c});
  r.total_energy = area * red_.sum(g.mass_cells(), FEnergy{s.u_2.cview(), s.v_2.cview(), s.w_2.cview(), s.t_2.cview(),
                                                           s.p.cview(), s.ph_2.cview(), s.mu_2.cview(), b.mub.view(),
                                                           b.pb.view(), b.phb.view(), c});
  return r;
}

std::vector<io::VarView> snapshot_vars(const State& s, const BaseState& b, const Grid& g,
                                       std::vector<std::vector<real>>& store,
                                       const std::vector<const Field2*>& extra_2d) {
  std::vector<real> hu, hv, hw, ht, hp, hph, hmu, bpb, bphb, bti, bmub;
  s.u_2.download(hu); s.v_2.download(hv); s.w_2.download(hw); s.t_2.download(ht);
  s.p.download(hp); s.ph_2.download(hph); s.mu_2.download(hmu);
  b.pb.download(bpb); b.phb.download(bphb); b.t_init.download(bti); b.mub.download(bmub);
  const io::DType dt = sizeof(real) == 8 ? io::DType::F64 : io::DType::F32;
  std::vector<io::VarView> vars;
  store.clear();
  store.reserve(16);
  auto add = [&](const char* name, char st, int nx, int ny, int nz, auto f) {
    store.emplace_back(size_t(nx) * ny * nz);
    auto& a = store.back();
    for (int k = 0; k < nz; ++k)
      for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) a[(size_t(k) * ny + j) * nx + i] = real(f(i, j, k));
    vars.push_back(io::VarView{name, st, nx, ny, nz, a.data(), dt, nx, std::ptrdiff_t(nx) * ny});
  };
  auto at = [](const std::vector<real>& h, const Layout& L, int i, int j, int k) { return double(h[L.index(i, j, k)]); };
  const Layout &Lu = s.u_2.layout(), &Lv = s.v_2.layout(), &Lw = s.w_2.layout(), &Lm = s.t_2.layout(),
               &Lph = s.ph_2.layout(), &Lmu = s.mu_2.layout(), &Lpb = b.pb.layout(), &Lphb = b.phb.layout(),
               &Lti = b.t_init.layout(), &Lmub = b.mub.layout();
  const int nx = g.nx, ny = g.ny, nz = g.nz;
  add("u", 'U', g.nxs(Stagger::U), ny, nz, [&](int i, int j, int k) { return at(hu, Lu, i, j, k); });
  add("v", 'V', nx, g.nys(Stagger::V), nz, [&](int i, int j, int k) { return at(hv, Lv, i, j, k); });
  add("w", 'W', nx, ny, nz + 1, [&](int i, int j, int k) { return at(hw, Lw, i, j, k); });
  std::vector<real> hqv;
  if (s.P_QV >= 0) s.q_2[s.P_QV].download(hqv);
  auto theta = [&](int i, int j, int k) {  // dry theta (t_2 holds theta_m - t0 when moist)
    const double tm = at(ht, Lm, i, j, k) + cst::t0;
    return hqv.empty() ? tm : tm / (1.0 + cst::rvovrd * at(hqv, Lm, i, j, k));
  };
  add("th_pert", 'M', nx, ny, nz, [&](int i, int j, int k) { return theta(i, j, k) - cst::t0 - at(bti, Lti, i, j, k); });
  add("theta", 'M', nx, ny, nz, theta);
  if (!hqv.empty()) add("theta_m", 'M', nx, ny, nz, [&](int i, int j, int k) { return at(ht, Lm, i, j, k) + cst::t0; });
  add("p", 'M', nx, ny, nz, [&](int i, int j, int k) { return at(hp, Lm, i, j, k) + at(bpb, Lpb, i, j, k); });
  add("p_pert", 'M', nx, ny, nz, [&](int i, int j, int k) { return at(hp, Lm, i, j, k); });
  add("ph", 'W', nx, ny, nz + 1, [&](int i, int j, int k) { return at(hph, Lph, i, j, k) + at(bphb, Lphb, i, j, k); });
  add("mu", 'M', nx, ny, 1, [&](int i, int j, int) { return at(hmu, Lmu, i, j, 0) + at(bmub, Lmub, i, j, 0); });
  add("mu_pert", 'M', nx, ny, 1, [&](int i, int j, int) { return at(hmu, Lmu, i, j, 0); });
  add("z_m", 'M', nx, ny, nz, [&](int i, int j, int k) {
    return 0.5 * (at(hph, Lph, i, j, k) + at(hph, Lph, i, j, k + 1) + at(bphb, Lphb, i, j, k) + at(bphb, Lphb, i, j, k + 1)) / cst::g;
  });
  std::vector<real> hq;
  for (int n = 0; n < s.n_scalars(); ++n) {  // moist species and tracers, under their names
    s.q_2[n].download(hq);
    add(s.q_names[n].c_str(), 'M', nx, ny, nz, [&](int i, int j, int k) { return at(hq, Lm, i, j, k); });
  }
  for (const Field2* f : extra_2d) {
    f->download(hq);
    add(f->name().c_str(), 'M', nx, ny, 1, [&](int i, int j, int) { return at(hq, f->layout(), i, j, 0); });
  }
  return vars;
}

}  // namespace cppwrf::dyn
