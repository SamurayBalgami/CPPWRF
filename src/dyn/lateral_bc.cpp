#include "dyn/lateral_bc.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "core/halo.hpp"

namespace cppwrf::dyn {
namespace {

using V3 = View3<real>;
using C3 = View3<const real>;
using V2 = View2<real>;
using C2 = View2<const real>;
using B2 = BView2<const real>;
using A1 = View1<const real>;

struct Slab {  // (pos, k, b), pos fastest
  const real* p = nullptr;
  int npos = 0, nk = 0;
  CPPWRF_HD real operator()(int pos, int k, int b) const {
    return p[pos + std::size_t(k) * npos + std::size_t(b) * npos * nk];
  }
};
struct SlabSet {
  Slab v[4], t[4];  // values, tendencies; sides xs, xe, ys, ye
  CPPWRF_HD real val(int side, int pos, int k, int b, real dtbc) const {
    return v[side](pos, k, b) + dtbc * t[side](pos, k, b);
  }
};

// Boundary ring (loop structure): sides in the order south, north, west, east;
// south/north rows cover b..ibe-b, west/east columns b+1..jbe-b-1 (no overlap).
// Side index 0..3 here = S, N, W, E; slabs are stored xs, xe, ys, ye.
CPPWRF_HD inline int slab_of_side(int side) { return side < 2 ? side + 2 : side - 2; }  // S,N,W,E -> ys,ye,xs,xe
struct Ring {
  int ibe, jbe;  // inclusive last index in x / y for this staggering
  int b0;        // first b of this launch
  CPPWRF_HD bool map(int pos, int side, int b, int& i, int& j) const {
    switch (side) {
      case 0: j = b; i = pos; return pos >= b && pos <= ibe - b;
      case 1: j = jbe - b; i = pos; return pos >= b && pos <= ibe - b;
      case 2: i = b; j = pos; return pos >= b + 1 && pos <= jbe - b - 1;
      default: i = ibe - b; j = pos; return pos >= b + 1 && pos <= jbe - b - 1;
    }
  }
};

Ring ring_for(const Grid& g, char var, int b0) {
  return Ring{var == 'u' ? g.nx : g.nx - 1, var == 'v' ? g.ny : g.ny - 1, b0};
}
Range3D ring_range(const Ring& r, int nk, int nb) {
  return Range3D{0, std::max(r.ibe, r.jbe) + 1, 0, nk, 0, 4 * nb};
}

struct KSpecTend {  // spec_bdytend: tend = boundary tendency in the spec zone
  V3 t3; V2 t2; bool two; SlabSet s; Ring r;
  CPPWRF_HD void operator()(int pos, int k, int sb) const {
    const int side = sb % 4, b = r.b0 + sb / 4;
    int i, j;
    if (!r.map(pos, side, b, i, j)) return;
    const real v = s.t[slab_of_side(side)](pos, k, b);
    if (two) t2(i, j) = v; else t3(i, j, k) = v;
  }
};

// relax_bdytend_core: Newtonian relaxation and smoothing of (boundary - model) in rows
// spec_zone .. spec_zone + relax_zone - 2. The model field is mass-weighted if c1 is given.
struct KRelax {
  V3 o3; V2 o2; C3 f3; C2 f2; bool two;
  C2 mut; A1 c1, c2; bool massw; B2 msf; bool use_msf;
  SlabSet s; Ring r; A1 fcx, gcx; real dtbc;
  CPPWRF_HD real F(int i, int j, int k) const {
    const real x = two ? f2(i, j) : f3(i, j, k);
    return massw ? x * (c1(k) * mut(i, j) + c2(k)) : x;
  }
  CPPWRF_HD void operator()(int pos, int k, int sb) const {
    const int side = sb % 4, b = r.b0 + sb / 4;
    int i, j;
    if (!r.map(pos, side, b, i, j)) return;
    const int sl = slab_of_side(side);
    auto B = [&](int p, int bb) { return s.val(sl, p, k, bb, dtbc); };
    real fls0, fls1, fls2, fls3, fls4;
    fls0 = B(pos, b) - F(i, j, k);
    if (side <= 1) {  // south / north row
      const int im1 = i - 1 > 0 ? i - 1 : 0, ip1 = i + 1 < r.ibe ? i + 1 : r.ibe;
      const int jin = side == 0 ? j + 1 : j - 1, jout = side == 0 ? j - 1 : j + 1;
      fls1 = B(im1, b) - F(im1, j, k);
      fls2 = B(ip1, b) - F(ip1, j, k);
      fls3 = B(i, b - 1) - F(i, jout, k);
      fls4 = B(i, b + 1) - F(i, jin, k);
    } else {  // west / east column
      const int iin = side == 2 ? i + 1 : i - 1, iout = side == 2 ? i - 1 : i + 1;
      fls1 = B(j - 1, b) - F(i, j - 1, k);
      fls2 = B(j + 1, b) - F(i, j + 1, k);
      fls3 = B(j, b - 1) - F(iout, j, k);
      fls4 = B(j, b + 1) - F(iin, j, k);
    }
    const real val = fcx(b) * fls0 - gcx(b) * (fls1 + fls2 + fls3 + fls4 - real(4) * fls0);
    if (two) o2(i, j) += val;
    else o3(i, j, k) += use_msf ? val * msf(i, j) : val;  // rk_addtend_dry: *_tendf += *_save * msf
  }
};

struct KUpdate {  // spec_bdyupdate: field += dts * tend in the spec zone
  V3 f3; C3 t3; V2 f2; C2 t2; bool two; Ring r; real dts;
  CPPWRF_HD void operator()(int pos, int k, int sb) const {
    const int side = sb % 4, b = r.b0 + sb / 4;
    int i, j;
    if (!r.map(pos, side, b, i, j)) return;
    if (two) f2(i, j) += dts * t2(i, j); else f3(i, j, k) += dts * t3(i, j, k);
  }
};

struct KUpdatePh {  // spec_bdyupdate_ph
  V3 ph; C3 ph_save, ph_tend; C2 mu_tend, muts; A1 c1f, c2f; Ring r; real dts;
  CPPWRF_HD void operator()(int pos, int k, int sb) const {
    const int side = sb % 4, b = r.b0 + sb / 4;
    int i, j;
    if (!r.map(pos, side, b, i, j)) return;
    const real mu_old = muts(i, j) - dts * mu_tend(i, j);
    const real m_old = c1f(k) * mu_old + c2f(k), m_new = c1f(k) * muts(i, j) + c2f(k);
    const real ratio = m_old / m_new;
    ph(i, j, k) = ph(i, j, k) * ratio + dts * ph_tend(i, j, k) / m_new + ph_save(i, j, k) * (ratio - real(1));
  }
};

struct KZeroGrad {  // zero_grad_bdy on w: copy from the first interior row / column
  V3 w; Ring r;
  CPPWRF_HD void operator()(int pos, int k, int sb) const {
    const int side = sb % 4, b = r.b0 + sb / 4;
    int i, j;
    if (!r.map(pos, side, b, i, j)) return;
    const int ic = i < 1 ? 1 : (i > r.ibe - 1 ? r.ibe - 1 : i);
    const int jc = j < 1 ? 1 : (j > r.jbe - 1 ? r.jbe - 1 : j);
    switch (side) {
      case 0: w(i, j, k) = w(ic, 1, k); break;
      case 1: w(i, j, k) = w(ic, r.jbe - 1, k); break;
      case 2: w(i, j, k) = w(1, jc, k); break;
      default: w(i, j, k) = w(r.ibe - 1, jc, k); break;
    }
  }
};

struct KFlowDep {  // flow_dep_bdy
  V3 f; C3 u, v; Ring r;
  CPPWRF_HD void operator()(int pos, int k, int sb) const {
    const int side = sb % 4, b = r.b0 + sb / 4;
    int i, j;
    if (!r.map(pos, side, b, i, j)) return;
    const int ic = i < 1 ? 1 : (i > r.ibe - 1 ? r.ibe - 1 : i);
    const int jc = j < 1 ? 1 : (j > r.jbe - 1 ? r.jbe - 1 : j);
    switch (side) {
      case 0: f(i, j, k) = v(i, j, k) < real(0) ? f(ic, 1, k) : real(0); break;
      case 1: f(i, j, k) = v(i, j + 1, k) > real(0) ? f(ic, r.jbe - 1, k) : real(0); break;
      case 2: f(i, j, k) = u(i, j, k) < real(0) ? f(1, jc, k) : real(0); break;
      default: f(i, j, k) = u(i + 1, j, k) > real(0) ? f(r.ibe - 1, jc, k) : real(0); break;
    }
  }
};

struct KFinal {  // spec_bdy_final: exact boundary values at the end of the step
  V3 f3; V2 f2; bool two; C2 mu; A1 c1, c2; bool mucouple; B2 msf; bool msfcouple;
  SlabSet s; Ring r; real dtbc;
  CPPWRF_HD void operator()(int pos, int k, int sb) const {
    const int side = sb % 4, b = r.b0 + sb / 4;
    int i, j;
    if (!r.map(pos, side, b, i, j)) return;
    const real bval = s.val(slab_of_side(side), pos, k, b, dtbc);
    const real xmsf = msfcouple ? msf(i, j) : real(1);
    const real xmu = mucouple ? (c1(k) * mu(i, j) + c2(k)) : real(1);
    const real v = xmsf * bval / xmu;
    if (two) f2(i, j) = v; else f3(i, j, k) = v;
  }
};

}  // namespace

LateralBC::LateralBC(const Grid& g, double dt, int width, int spec_zone, int relax_zone, bool nested)
    : g_(g), width_(width), spec_zone_(spec_zone), relax_zone_(relax_zone), nested_(nested),
      coef_(2 * std::size_t(width)) {
  if (g.two_d()) throw std::runtime_error("LateralBC: specified boundaries need a 3D grid");
  // lbc_fcx_gcx (module_bc_em.F): the same ramp for specified (spongeweight = 1 with spec_exp = 0)
  // and nested domains
  fcx_h_.assign(width, 0.0);
  gcx_h_.assign(width, 0.0);
  for (int loop = spec_zone + 1; loop <= spec_zone + relax_zone; ++loop) {  // Fortran index
    const double fac = double(spec_zone + relax_zone - loop) / double(relax_zone - 1);
    fcx_h_[loop - 1] = 0.1 / dt * fac;
    gcx_h_[loop - 1] = 1.0 / dt / 50.0 * fac;
  }
  std::vector<real> h(2 * std::size_t(width));
  for (int b = 0; b < width; ++b) { h[b] = real(fcx_h_[b]); h[width + b] = real(gcx_h_[b]); }
  copy_to_dev(coef_.get(), h.data(), h.size() * sizeof(real));
  const int nz = g.nz;
  nk_ = {nz, nz, nz, nz + 1, 1, nz + 1};  // U V T PH MU W
  for (int v = 0; v < NVAR; ++v) {
    if (v == W && !nested) continue;
    std::size_t o = 0;
    for (int half = 0; half < 2; ++half)
      for (int s = 0; s < 4; ++s) {
        off_[v][half * 4 + s] = o;
        o += std::size_t(npos(Var(v), Side(s))) * nk_[v] * width;
      }
    data_[v] = Buffer(o);
  }
}

int LateralBC::npos(Var v, Side s) const {
  if (s == XS || s == XE) return v == V ? g_.ny + 1 : g_.ny;
  return v == U ? g_.nx + 1 : g_.nx;
}

void LateralBC::set(Var v, Side s, const std::vector<double>& values, const std::vector<double>& tends) {
  const std::size_t n = std::size_t(npos(v, s)) * nk_[v] * width_;
  if (values.size() != n || tends.size() != n) throw std::runtime_error("LateralBC::set: wrong slab size");
  std::vector<real> hv(values.begin(), values.end()), ht(tends.begin(), tends.end());
  copy_to_dev(data_[v].get() + off_[v][s], hv.data(), n * sizeof(real));
  copy_to_dev(data_[v].get() + off_[v][4 + s], ht.data(), n * sizeof(real));
  have_[v] = true;
}

real* LateralBC::scalar_slab(int n, Side s, bool tend) {
  if (n < 0) throw std::runtime_error("LateralBC: bad scalar index");
  if (n >= int(sdata_.size())) sdata_.resize(n + 1);
  ScalarData& d = sdata_[n];
  if (d.off.empty()) {
    const std::size_t per_x = std::size_t(g_.ny) * g_.nz * width_, per_y = std::size_t(g_.nx) * g_.nz * width_;
    std::size_t o = 0;
    for (int half = 0; half < 2; ++half)
      for (int sd = 0; sd < 4; ++sd) {
        d.off.push_back(o);
        o += sd < 2 ? per_x : per_y;
      }
    d.data = Buffer(o);
  }
  return d.data.get() + d.off[(tend ? 4 : 0) + s];
}

void LateralBC::set_scalar(int n, Side s, const std::vector<double>& values, const std::vector<double>& tends) {
  scalar_slab(n, s, false);  // allocate
  ScalarData& d = sdata_[n];
  const std::size_t per_x = std::size_t(g_.ny) * g_.nz * width_, per_y = std::size_t(g_.nx) * g_.nz * width_;
  const std::size_t nn = s < 2 ? per_x : per_y;
  if (values.size() != nn || tends.size() != nn) throw std::runtime_error("LateralBC::set_scalar: wrong slab size");
  std::vector<real> hv(values.begin(), values.end()), ht(tends.begin(), tends.end());
  copy_to_dev(d.data.get() + d.off[s], hv.data(), nn * sizeof(real));
  copy_to_dev(d.data.get() + d.off[4 + s], ht.data(), nn * sizeof(real));
}

namespace {
SlabSet slabs(const Buffer& d, const std::array<std::size_t, 8>& off, int nk, int nposx, int nposy) {
  SlabSet s;
  for (int half = 0; half < 2; ++half)
    for (int side = 0; side < 4; ++side) {
      Slab sl{d.get() + off[half * 4 + side], side < 2 ? nposx : nposy, nk};
      (half == 0 ? s.v : s.t)[side] = sl;
    }
  return s;
}
SlabSet scalar_slabs(const Buffer& d, const std::vector<std::size_t>& off, int nk, int nposx, int nposy) {
  std::array<std::size_t, 8> a{};
  for (int i = 0; i < 8; ++i) a[i] = off[i];
  return slabs(d, a, nk, nposx, nposy);
}
}  // namespace

void LateralBC::relax_dry(const State& s, Tendencies& t, const BaseState& b) const {
  const A1 fcx{coef_.get()}, gcx{coef_.get() + width_};
  const int nb = relax_zone_ - 1;  // rows spec_zone .. spec_zone + relax_zone - 2
  auto sl = [&](Var v) { return slabs(data_[v], off_[v], nk_[v], npos(v, XS), npos(v, YS)); };
  const real dtbc_r = real(dtbc);
  {  // ru -> ru_tendf (u_save * msfuy)
    const Ring r = ring_for(g_, 'u', spec_zone_);
    parallel_for(ring_range(r, g_.nz, nb), KRelax{t.ru_tendf.view(), {}, s.ru.cview(), {}, false, {}, {}, {}, false,
                                                  b.msfuy.view(), true, sl(U), r, fcx, gcx, dtbc_r});
  }
  {
    const Ring r = ring_for(g_, 'v', spec_zone_);
    parallel_for(ring_range(r, g_.nz, nb), KRelax{t.rv_tendf.view(), {}, s.rv.cview(), {}, false, {}, {}, {}, false,
                                                  b.msfvx.view(), true, sl(V), r, fcx, gcx, dtbc_r});
  }
  {  // ph, mass-weighted with c1f/c2f and mut (nz+1 levels)
    const Ring r = ring_for(g_, 'h', spec_zone_);
    parallel_for(ring_range(r, g_.nz + 1, nb),
                 KRelax{t.ph_tendf.view(), {}, s.ph_2.cview(), {}, false, s.mut.cview(), b.c1f.view(), b.c2f.view(),
                        true, {}, false, sl(PH), r, fcx, gcx, dtbc_r});
  }
  {  // t, mass-weighted with c1h/c2h
    const Ring r = ring_for(g_, 't', spec_zone_);
    parallel_for(ring_range(r, g_.nz, nb),
                 KRelax{t.t_tendf.view(), {}, s.t_2.cview(), {}, false, s.mut.cview(), b.c1h.view(), b.c2h.view(), true,
                        {}, false, sl(T), r, fcx, gcx, dtbc_r});
  }
  {  // mu -> mu_tend (2D)
    const Ring r = ring_for(g_, 'm', spec_zone_);
    parallel_for(ring_range(r, 1, nb), KRelax{{}, t.mu_tend.view(), {}, s.mu_2.cview(), true, {}, {}, {}, false, {},
                                              false, sl(MU), r, fcx, gcx, dtbc_r});
  }
  if (nested_ && have_[W]) {  // nests: w mass-weighted (c1f, c2f, mut) -> w_save -> rw_tendf (* msfty)
    const Ring r = ring_for(g_, 'h', spec_zone_);
    parallel_for(ring_range(r, g_.nz + 1, nb),
                 KRelax{t.rw_tendf.view(), {}, s.w_2.cview(), {}, false, s.mut.cview(), b.c1f.view(), b.c2f.view(),
                        true, b.msfty.view(), true, sl(W), r, fcx, gcx, dtbc_r});
  }
}

void LateralBC::spec_dry(Tendencies& t) const {
  auto sl = [&](Var v) { return slabs(data_[v], off_[v], nk_[v], npos(v, XS), npos(v, YS)); };
  const int nb = spec_zone_;
  auto r = [&](char c) { return ring_for(g_, c, 0); };
  parallel_for(ring_range(r('u'), g_.nz, nb), KSpecTend{t.ru_tend.view(), {}, false, sl(U), r('u')});
  parallel_for(ring_range(r('v'), g_.nz, nb), KSpecTend{t.rv_tend.view(), {}, false, sl(V), r('v')});
  parallel_for(ring_range(r('h'), g_.nz + 1, nb), KSpecTend{t.ph_tend.view(), {}, false, sl(PH), r('h')});
  parallel_for(ring_range(r('t'), g_.nz, nb), KSpecTend{t.t_tend.view(), {}, false, sl(T), r('t')});
  parallel_for(ring_range(r('m'), 1, nb), KSpecTend{{}, t.mu_tend.view(), true, sl(MU), r('m')});
  if (nested_ && have_[W])  // spec_bdytend rw_tend (nests)
    parallel_for(ring_range(r('h'), g_.nz + 1, nb), KSpecTend{t.rw_tend.view(), {}, false, sl(W), r('h')});
}

void LateralBC::update_uv(State& s, const Tendencies& t, real dts) const {
  const int nb = spec_zone_;
  const Ring ru = ring_for(g_, 'u', 0), rv = ring_for(g_, 'v', 0);
  parallel_for(ring_range(ru, g_.nz, nb), KUpdate{s.u_2.view(), t.ru_tend.cview(), {}, {}, false, ru, dts});
  parallel_for(ring_range(rv, g_.nz, nb), KUpdate{s.v_2.view(), t.rv_tend.cview(), {}, {}, false, rv, dts});
}

void LateralBC::update_mu_t(State& s, const Tendencies& t, real dts) const {
  const int nb = spec_zone_;
  const Ring rm = ring_for(g_, 't', 0);
  parallel_for(ring_range(rm, g_.nz, nb), KUpdate{s.t_2.view(), t.t_tend.cview(), {}, {}, false, rm, dts});
  parallel_for(ring_range(rm, 1, nb), KUpdate{{}, {}, s.mu_2.view(), t.mu_tend.cview(), true, rm, dts});
  parallel_for(ring_range(rm, 1, nb), KUpdate{{}, {}, s.muts.view(), t.mu_tend.cview(), true, rm, dts});
}

void LateralBC::update_ph_w(State& s, const Field3& ph_save, const Tendencies& t, const BaseState& b, real dts) const {
  const int nb = spec_zone_;
  const Ring rh = ring_for(g_, 'h', 0);
  parallel_for(ring_range(rh, g_.nz + 1, nb), KUpdatePh{s.ph_2.view(), ph_save.cview(), t.ph_tend.cview(),
                                                         t.mu_tend.cview(), s.muts.cview(), b.c1f.view(),
                                                         b.c2f.view(), rh, dts});
  if (nested_ && have_[W])  // nests: spec_bdyupdate(w, rw_tend) (solve_em 1603)
    parallel_for(ring_range(rh, g_.nz + 1, nb), KUpdate{s.w_2.view(), t.rw_tend.cview(), {}, {}, false, rh, dts});
  else
    parallel_for(ring_range(rh, g_.nz + 1, nb), KZeroGrad{s.w_2.view(), rh});
}

void LateralBC::relax_scalar(int n, Field3& q_tend, const Field3& q, const State& s, const BaseState& b) const {
  if (!has_scalar(n)) return;
  const A1 fcx{coef_.get()}, gcx{coef_.get() + width_};
  const Ring r = ring_for(g_, 'q', spec_zone_);
  const ScalarData& d = sdata_[n];
  parallel_for(ring_range(r, g_.nz, relax_zone_ - 1),
               KRelax{q_tend.view(), {}, q.cview(), {}, false, s.mut.cview(), b.c1h.view(), b.c2h.view(), true, {},
                      false, scalar_slabs(d.data, d.off, g_.nz, g_.ny, g_.nx), r, fcx, gcx, real(dtbc)});
}

void LateralBC::spec_scalar(int n, Field3& q_tend) const {
  if (!has_scalar(n)) return;
  const Ring r = ring_for(g_, 'q', 0);
  const ScalarData& d = sdata_[n];
  parallel_for(ring_range(r, g_.nz, spec_zone_),
               KSpecTend{q_tend.view(), {}, false, scalar_slabs(d.data, d.off, g_.nz, g_.ny, g_.nx), r});
}

void LateralBC::flow_dep(Field3& q, const State& s) const {
  const Ring r = ring_for(g_, 'q', 0);
  parallel_for(ring_range(r, g_.nz, spec_zone_), KFlowDep{q.view(), s.ru_m.cview(), s.rv_m.cview(), r});
}

void LateralBC::final(State& s, const BaseState& b) const {
  auto sl = [&](Var v) { return slabs(data_[v], off_[v], nk_[v], npos(v, XS), npos(v, YS)); };
  const int nb = spec_zone_;
  const real d = real(dtbc);
  auto r = [&](char c) { return ring_for(g_, c, 0); };
  parallel_for(ring_range(r('u'), g_.nz, nb), KFinal{s.u_2.view(), {}, false, s.muus.cview(), b.c1h.view(),
                                                     b.c2h.view(), true, b.msfuy.view(), true, sl(U), r('u'), d});
  parallel_for(ring_range(r('v'), g_.nz, nb), KFinal{s.v_2.view(), {}, false, s.muvs.cview(), b.c1h.view(),
                                                     b.c2h.view(), true, b.msfvx.view(), true, sl(V), r('v'), d});
  parallel_for(ring_range(r('t'), g_.nz, nb), KFinal{s.t_2.view(), {}, false, s.muts.cview(), b.c1h.view(),
                                                     b.c2h.view(), true, b.msfty.view(), false, sl(T), r('t'), d});
  parallel_for(ring_range(r('h'), g_.nz + 1, nb), KFinal{s.ph_2.view(), {}, false, s.muts.cview(), b.c1f.view(),
                                                         b.c2f.view(), true, b.msfty.view(), false, sl(PH), r('h'), d});
  parallel_for(ring_range(r('m'), 1, nb), KFinal{{}, s.mu_2.view(), true, s.muts.cview(), b.c1h.view(), b.c2h.view(),
                                                 false, b.msfty.view(), false, sl(MU), r('m'), d});
  for (int n = 0; n < s.n_scalars() && n < int(sdata_.size()); ++n) {  // moist species with data ('t' type)
    if (!has_scalar(n)) continue;
    const ScalarData& sd = sdata_[n];
    parallel_for(ring_range(r('t'), g_.nz, nb),
                 KFinal{s.q_2[n].view(), {}, false, s.muts.cview(), b.c1h.view(), b.c2h.view(), true, b.msfty.view(),
                        false, scalar_slabs(sd.data, sd.off, g_.nz, g_.ny, g_.nx), r('t'), d});
  }
  if (nested_ && have_[W])  // spec_bdy_final 'w': mu (c1f, c2f) and msfty coupled
    parallel_for(ring_range(r('h'), g_.nz + 1, nb), KFinal{s.w_2.view(), {}, false, s.muts.cview(), b.c1f.view(),
                                                           b.c2f.view(), true, b.msfty.view(), true, sl(W), r('h'), d});
}

}  // namespace cppwrf::dyn
