#include "dyn/scalars.hpp"

#include "core/vm_dump.hpp"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <vector>

#include "core/backend.hpp"
#include "core/halo.hpp"
#include "dyn/diffusion_driver.hpp"
#include "dyn/kernels_common.hpp"
#include "dyn/lateral_bc.hpp"

// WRF dyn_em/module_advect_em.F (advect_scalar, advect_scalar_pd) and module_em.F
// (rk_scalar_tend, rk_update_scalar, rk_update_scalar_pd). C++ k = Fortran k - 1.
// With map factors (msftx, msfty as in advect_scalar*) and periodic, symmetric or specified
// lateral boundaries (degraded flux orders and the spec zone as WRF for specified).

namespace cppwrf::dyn {
namespace {

using V3 = View3<real>;
using C3 = View3<const real>;
using C2 = View2<const real>;
using B2 = BView2<const real>;
constexpr real H = real(0.5);

// High-order (advect_scalar) and first-order upwind (advect_scalar_pd) face fluxes.
struct Fluxes {
  C3 q, qo, ru, rv, ww; C2 mut; VCoef c;
  real rdx, rdy, dt, hup, vup; int nz; bool do_y;
  B2 msftx, msfty; kern::Spec sp;
  bool pd = false;  // advect_scalar_pd's flux statement functions
  CPPWRF_HD real f5(real a, real b, real cc, real d, real e, real f, real v) const {
    return pd ? kern::flux5_pd(a, b, cc, d, e, f, v, hup) : kern::flux5(a, b, cc, d, e, f, v, hup);
  }
  CPPWRF_HD real f3(real a, real b, real cc, real d, real v, real up) const {
    return pd ? kern::flux3_pd(a, b, cc, d, v, up) : kern::flux3(a, b, cc, d, v, up);
  }
  CPPWRF_HD real fx(int i, int j, int k) const {  // west face of cell i
    const real vel = ru(i, j, k);
    const int o = kern::ord_face(i, sp.nx, sp.x);
    if (o == 2) return H * vel * (q(i, j, k) + q(i - 1, j, k));
    if (o == 3) return vel * f3(q(i - 2, j, k), q(i - 1, j, k), q(i, j, k), q(i + 1, j, k), vel, hup);
    return vel * f5(q(i - 3, j, k), q(i - 2, j, k), q(i - 1, j, k), q(i, j, k), q(i + 1, j, k), q(i + 2, j, k), vel);
  }
  CPPWRF_HD real fy(int i, int j, int k) const {
    const real vel = rv(i, j, k);
    const int o = kern::ord_face(j, sp.ny, sp.y);
    if (o == 2) return H * vel * (q(i, j, k) + q(i, j - 1, k));
    if (o == 3) return vel * f3(q(i, j - 2, k), q(i, j - 1, k), q(i, j, k), q(i, j + 1, k), vel, hup);
    return vel * f5(q(i, j - 3, k), q(i, j - 2, k), q(i, j - 1, k), q(i, j, k), q(i, j + 1, k), q(i, j + 2, k), vel);
  }
  CPPWRF_HD real fz(int i, int j, int k) const {  // w level k
    if (k == 0 || k == nz) return 0;
    const real vel = ww(i, j, k);
    if (k == 1 || k == nz - 1) return vel * (c.fnm(k) * q(i, j, k) + c.fnp(k) * q(i, j, k - 1));
    return vel * f3(q(i, j, k - 2), q(i, j, k - 1), q(i, j, k), q(i, j, k + 1), -vel, vup);
  }
  CPPWRF_HD real mh(int i, int j, int k) const { return c.c1h(k) * mut(i, j) + c.c2h(k); }
  CPPWRF_HD real fxl(int i, int j, int k) const {
    const real dx = real(2) / (msfty(i, j) + msfty(i - 1, j)) / rdx, mu = H * (mh(i, j, k) + mh(i - 1, j, k)),
               vel = ru(i, j, k);
    // advect_scalar_pd, degrade_xs, i == ids+1 only: vel = ru/mu, cr = vel*dt/dx
    const real cr = (sp.x && i == 1) ? (vel / mu) * dt / dx : vel * dt / dx / mu;
    return mu * (dx / dt) * kern::flux_upwind(qo(i - 1, j, k), qo(i, j, k), cr);
  }
  CPPWRF_HD real fyl(int i, int j, int k) const {
    const real dy = real(2) / (msftx(i, j) + msftx(i, j - 1)) / rdy, mu = H * (mh(i, j, k) + mh(i, j - 1, k)),
               vel = rv(i, j, k);
    const real cr = vel * dt / dy / mu;
    return mu * (dy / dt) * kern::flux_upwind(qo(i, j - 1, k), qo(i, j, k), cr);
  }
  CPPWRF_HD real fzl(int i, int j, int k) const {
    if (k == 0 || k == nz) return 0;
    const real dz = real(2) / (c.rdnw(k) + c.rdnw(k - 1)), mu = mh(i, j, k), vel = ww(i, j, k);
    const real cr = vel * dt / dz / mu;
    return mu * (dz / dt) * kern::flux_upwind(qo(i, j, k - 1), qo(i, j, k), cr);
  }
  // spec zone of a specified boundary: WRF advect_scalar(_pd) does not compute it (tendency 0)
  CPPWRF_HD bool spec_cell(int i, int j) const {
    return (sp.sx() && (i == 0 || i == sp.nx - 1)) || (sp.sy() && (j == 0 || j == sp.ny - 1));
  }
  // outermost row of a specified or open boundary (advect_scalar_pd does not limit it)
  CPPWRF_HD bool edge_cell(int i, int j) const {
    return (sp.x && (i == 0 || i == sp.nx - 1)) || (sp.y && (j == 0 || j == sp.ny - 1));
  }
  CPPWRF_HD bool open_col(int i) const { return sp.ox && (i == 0 || i == sp.nx - 1); }
  CPPWRF_HD bool open_row(int j) const { return sp.oy && (j == 0 || j == sp.ny - 1); }
  // open boundaries: WRF adds the open terms (x, then y) to the advection tendency; field_old is
  // scalar_old in advect_scalar_pd, scalar itself in advect_scalar
  CPPWRF_HD real open_terms(real r, int i, int j, int k, C3 fo) const {
    if (open_col(i)) r = r - kern::open_x_scalar(q, fo, ru, rdx, i, j, k, sp.nx);
    if (open_row(j)) r = r - kern::open_y_scalar(q, fo, rv, rdy, i, j, k, sp.ny);
    return r;
  }
  // advect_scalar's tendency in a spec-zone cell (set_tend only): the x flux divergence is skipped
  // in the outer columns, the y one in the outer rows, the vertical term is kept
  CPPWRF_HD real adv_spec(int i, int j, int k) const {
    const bool ocol = sp.sx() && (i == 0 || i == sp.nx - 1), orow = sp.sy() && (j == 0 || j == sp.ny - 1);
    real r = real(0);
    if (!ocol) r = -(msftx(i, j) * rdx) * (fx(i + 1, j, k) - fx(i, j, k));
    if (do_y && !orow) r -= (msftx(i, j) * rdy) * (fy(i, j + 1, k) - fy(i, j, k));
    return r - c.rdnw(k) * (fz(i, j, k + 1) - fz(i, j, k));
  }
  CPPWRF_HD real adv(int i, int j, int k) const {  // advect_scalar: horizontal terms times msftx
    if (spec_cell(i, j)) return 0;
    real r;
    if (open_col(i) || open_row(j)) {  // WRF order: y flux, x flux, open x, open y, vertical
      r = real(0);
      if (do_y && !open_row(j)) r = r - (msftx(i, j) * rdy) * (fy(i, j + 1, k) - fy(i, j, k));
      if (!open_col(i)) r = r - (msftx(i, j) * rdx) * (fx(i + 1, j, k) - fx(i, j, k));
      r = open_terms(r, i, j, k, q);
    } else {
      r = -(msftx(i, j) * rdx) * (fx(i + 1, j, k) - fx(i, j, k));
      if (do_y) r -= (msftx(i, j) * rdy) * (fy(i, j + 1, k) - fy(i, j, k));
    }
    return r - c.rdnw(k) * (fz(i, j, k + 1) - fz(i, j, k));
  }
};

CPPWRF_HD inline real m_full(const VCoef& c, int k, real mub, real mu) { return c.c1h(k) * (mub + mu) + c.c2h(k); }

// rk_update_scalar (non-PD stages / options): q_new = (m_old q_1 + dt (adv + tendf)) / m_new
struct KUpdate {
  Fluxes f; C3 q1, tendf; C2 mu1, mu2; B2 mub; V3 qn; real dt_rk; bool add_tendf;
  V3 ften; bool want_ften;  // set_tend: advect_tend * msfty
  CPPWRF_HD void operator()(int i, int j, int k) const {
    real tend = f.adv(i, j, k) * f.msfty(i, j);  // rk_update_scalar: advect_tend * msfty + sc_tend
    if (want_ften) ften(i, j, k) = f.spec_cell(i, j) ? f.adv_spec(i, j, k) * f.msfty(i, j) : tend;
    if (add_tendf) tend += tendf(i, j, k);
    qn(i, j, k) = (m_full(f.c, k, mub(i, j), mu1(i, j)) * q1(i, j, k) + dt_rk * tend) /
                  m_full(f.c, k, mub(i, j), mu2(i, j));
  }
};

// rk_update_scalar_pd (last stage, before advection): q_1 += dt tendf / m_old
struct KPreUpdatePD {
  V3 q1; C3 tendf; C2 mu1; B2 mub; VCoef c; real dt;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    const real m = m_full(c, k, mub(i, j), mu1(i, j));
    q1(i, j, k) = (m * q1(i, j, k) + dt * tendf(i, j, k)) / m;
  }
};

// advect_scalar_pd, step 1: every face flux once. dh = high-order minus first-order upwind
// flux, lo = upwind flux, on x (U points), y (V points) and eta (W levels) faces. The limiter and
// the update read them (the same expressions as evaluating Fluxes per cell: bitwise identical).
struct KFacePDX {
  Fluxes f; V3 dh, lo;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    const real l = f.fxl(i, j, k);
    dh(i, j, k) = f.fx(i, j, k) - l;
    lo(i, j, k) = l;
  }
};
struct KFacePDY {
  Fluxes f; V3 dh, lo;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    const real l = f.fyl(i, j, k);
    dh(i, j, k) = f.fy(i, j, k) - l;
    lo(i, j, k) = l;
  }
};
struct KFacePDZ {
  Fluxes f; V3 dh, lo;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    const real l = f.fzl(i, j, k);
    dh(i, j, k) = f.fz(i, j, k) - l;
    lo(i, j, k) = l;
  }
};

struct FacesPD {
  C3 dhx, lox, dhy, loy, dhz, loz;
};

// advect_scalar_pd limiter: scale factor of the outgoing high-minus-low fluxes of each cell
struct KScalePD {
  Fluxes f; FacesPD a; C2 mu_old; B2 mub; V3 scale;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    if (f.edge_cell(i, j)) { scale(i, j, k) = 1; return; }  // WRF limits interior cells only
    const VCoef& c = f.c;
    const real mxy = f.msftx(i, j) * f.msfty(i, j);
    real hlow = f.rdx * (a.lox(i + 1, j, k) - a.lox(i, j, k));
    if (f.do_y) hlow += f.rdy * (a.loy(i, j + 1, k) - a.loy(i, j, k));
    const real ph_low = ((c.c1h(k) * mub(i, j) + c.c2h(k)) + c.c1h(k) * mu_old(i, j)) * f.qo(i, j, k) -
                        f.dt * (mxy * hlow + f.msfty(i, j) * c.rdnw(k) * (a.loz(i, j, k + 1) - a.loz(i, j, k)));
    const real fxp = a.dhx(i + 1, j, k), fxm = a.dhx(i, j, k);
    real hout = f.rdx * ((fxp > 0 ? fxp : real(0)) - (fxm < 0 ? fxm : real(0)));
    if (f.do_y) {
      const real fyp = a.dhy(i, j + 1, k), fym = a.dhy(i, j, k);
      hout += f.rdy * ((fyp > 0 ? fyp : real(0)) - (fym < 0 ? fym : real(0)));
    }
    const real fzp = a.dhz(i, j, k + 1), fzm = a.dhz(i, j, k);
    const real out = mxy * hout + f.msfty(i, j) * c.rdnw(k) * ((fzp < 0 ? fzp : real(0)) - (fzm > 0 ? fzm : real(0)));
    const real flux_out = f.dt * out;
    real s = 1;
    if (flux_out > ph_low) {
      s = ph_low / (flux_out + real(1e-20));
      if (s < 0) s = 0;
    }
    scale(i, j, k) = s;
  }
};

// PD-limited flux divergence + rk_update_scalar (tendf already folded into q_1)
struct KUpdatePD {
  Fluxes f; FacesPD a; C3 scale, q1; C2 mu1, mu2; B2 mub; V3 qn; real dt_rk;
  // limited high-minus-low flux on the west face of (i,j,k): scaled by the upwind (outflow) cell
  CPPWRF_HD real dx(int i, int j, int k) const {
    const real d = a.dhx(i, j, k);
    return d * (d > 0 ? scale(i - 1, j, k) : (d < 0 ? scale(i, j, k) : real(1)));
  }
  CPPWRF_HD real dy(int i, int j, int k) const {
    const real d = a.dhy(i, j, k);
    return d * (d > 0 ? scale(i, j - 1, k) : (d < 0 ? scale(i, j, k) : real(1)));
  }
  CPPWRF_HD real dz(int i, int j, int k) const {  // eta decreases upward: fz > 0 leaves cell k
    if (k == 0 || k == f.nz) return 0;
    const real d = a.dhz(i, j, k);
    return d * (d > 0 ? scale(i, j, k) : (d < 0 ? scale(i, j, k - 1) : real(1)));
  }
  V3 ften; bool want_ften;  // set_tend: advect_tend * msfty
  CPPWRF_HD void operator()(int i, int j, int k) const {  // WRF advect_scalar_pd final tendency
    if (f.spec_cell(i, j)) {
      if (want_ften) {  // WRF still forms the tendency here: x term off in outer columns, y off in outer rows
        const bool ocol = f.sp.sx() && (i == 0 || i == f.sp.nx - 1), orow = f.sp.sy() && (j == 0 || j == f.sp.ny - 1);
        real tend = -f.c.rdnw(k) * (dz(i, j, k + 1) - dz(i, j, k) + a.loz(i, j, k + 1) - a.loz(i, j, k));
        if (!ocol) tend -= f.msftx(i, j) * (f.rdx * (dx(i + 1, j, k) - dx(i, j, k) + a.lox(i + 1, j, k) - a.lox(i, j, k)));
        if (f.do_y && !orow)
          tend -= f.msftx(i, j) * (f.rdy * (dy(i, j + 1, k) - dy(i, j, k) + a.loy(i, j + 1, k) - a.loy(i, j, k)));
        ften(i, j, k) = tend * f.msfty(i, j);
      }
      qn(i, j, k) = m_full(f.c, k, mub(i, j), mu1(i, j)) * q1(i, j, k) / m_full(f.c, k, mub(i, j), mu2(i, j));
      return;
    }
    real tend;
    if (f.open_col(i) || f.open_row(j)) {  // open boundary: open terms, vertical, x (not outer columns), y
      tend = f.open_terms(real(0), i, j, k, f.qo);
      tend = tend - f.c.rdnw(k) * (dz(i, j, k + 1) - dz(i, j, k) + a.loz(i, j, k + 1) - a.loz(i, j, k));
      if (!f.open_col(i))
        tend -= f.msftx(i, j) * (f.rdx * (dx(i + 1, j, k) - dx(i, j, k) + a.lox(i + 1, j, k) - a.lox(i, j, k)));
      if (f.do_y && !f.open_row(j))
        tend -= f.msftx(i, j) * (f.rdy * (dy(i, j + 1, k) - dy(i, j, k) + a.loy(i, j + 1, k) - a.loy(i, j, k)));
    } else {
      tend = -f.c.rdnw(k) * (dz(i, j, k + 1) - dz(i, j, k) + a.loz(i, j, k + 1) - a.loz(i, j, k));
      tend -= f.msftx(i, j) * (f.rdx * (dx(i + 1, j, k) - dx(i, j, k) + a.lox(i + 1, j, k) - a.lox(i, j, k)));
      if (f.do_y) tend -= f.msftx(i, j) * (f.rdy * (dy(i, j + 1, k) - dy(i, j, k) + a.loy(i, j + 1, k) - a.loy(i, j, k)));
    }
    if (want_ften) ften(i, j, k) = tend * f.msfty(i, j);
    qn(i, j, k) = (m_full(f.c, k, mub(i, j), mu1(i, j)) * q1(i, j, k) + dt_rk * (tend * f.msfty(i, j))) /
                  m_full(f.c, k, mub(i, j), mu2(i, j));
  }
};

}  // namespace

ScalarTransport::Work::Work(const Grid& g)
    : q_new(g, Stagger::M, "q_new"), scale(g, Stagger::M, "pd_scale"), dhx(g, Stagger::U, "pd_dhx"),
      lox(g, Stagger::U, "pd_lox"), dhy(g, Stagger::V, "pd_dhy"), loy(g, Stagger::V, "pd_loy"),
      dhz(g, Stagger::W, "pd_dhz"), loz(g, Stagger::W, "pd_loz") {}

ScalarTransport::ScalarTransport(const Grid& g, const ScalarParams& p) : g_(g), p_(p) {}

void ScalarTransport::reserve(int n_scalars) {
  const int n = n_scalars < kBatch ? n_scalars : kBatch;
  while (int(work_.size()) < n) work_.push_back(std::make_unique<Work>(g_));
}

// rk_scalar_tend + rk_update_scalar (pd = false) or rk_update_scalar_pd (pd = true) of the species ns, kBatch at a
// time: each kernel runs for all species of a batch in one launch. The species are independent (each reads its
// own q_1, q_2, tendency and work fields and the common mass fluxes / masses, and writes only its own), so the
// values are those of the species one after the other (WRF's loop over the species).
void ScalarTransport::update(State& s, Tendencies& t, const BaseState& b, const std::vector<int>& ns, bool pd, int rk,
                             real dt_rk) {
  const Grid& g = g_;
  const VCoef c = VCoef::from(b);
  const Range3D mr = g.mass_cells();
  const bool do_y = !g.two_d();
  const real hup = p_.h_adv_order == 5 ? real(1) : real(0), vup = p_.v_adv_order == 3 ? real(1) : real(0);
  // debug (CPPWRF_VM_DUMP): the qv advection tendency * msfty
  if (vmdump::on() && !qvf_ && !dbg_qvf_) dbg_qvf_ = std::make_unique<Field3>(g, Stagger::M, "dbg_qvf");
  Field3* fsink = qvf_ ? qvf_ : (vmdump::on() ? dbg_qvf_.get() : nullptr);
  reserve(int(ns.size()));
  for (std::size_t n0 = 0; n0 < ns.size(); n0 += kBatch) {
    const int nb = int(std::min<std::size_t>(kBatch, ns.size() - n0));
    auto species = [&](int m) { return ns[n0 + m]; };
    auto fluxes = [&](int n) {
      Fluxes fl{s.q_2[n].cview(), s.q_1[n].cview(), s.ru_m.cview(), s.rv_m.cview(), s.ww_m.cview(), s.muts.cview(), c,
                real(b.rdx), real(b.rdy), dt_rk, hup, vup, g.nz, do_y, b.msftx.view(), b.msfty.view(),
                kern::Spec::of(g)};
      fl.pd = pd;
      return fl;
    };
    auto ften = [&](int n) { return fsink != nullptr && n == s.P_QV ? fsink->view() : V3{}; };
    if (!pd) {
      LaunchBatch<KUpdate, kBatch> up;
      for (int m = 0; m < nb; ++m) {
        const int n = species(m);
        up.add(mr, KUpdate{fluxes(n), s.q_1[n].cview(), t.q_tendf[n].cview(), s.mu_1.cview(), s.mu_2.cview(),
                           b.mub.view(), work_[m]->q_new.view(), dt_rk, true, ften(n), ften(n).p != nullptr});
      }
      up.run();
    } else {
      LaunchBatch<KPreUpdatePD, kBatch> pre;
      std::vector<Field3*> q1s;
      for (int m = 0; m < nb; ++m) {
        const int n = species(m);
        pre.add(mr, KPreUpdatePD{s.q_1[n].view(), t.q_tendf[n].cview(), s.mu_1.cview(), b.mub.view(), c, dt_rk});
        q1s.push_back(&s.q_1[n]);
      }
      pre.run();
      fill_halos(q1s, g);
      {
        LaunchBatch<KFacePDX, kBatch> fx;
        LaunchBatch<KFacePDY, kBatch> fy;
        LaunchBatch<KFacePDZ, kBatch> fz;
        for (int m = 0; m < nb; ++m) {
          Work& w = *work_[m];
          const Fluxes fl = fluxes(species(m));
          fx.add(Range3D{0, g.nx + 1, 0, g.ny, 0, g.nz}, KFacePDX{fl, w.dhx.view(), w.lox.view()});
          if (do_y) fy.add(Range3D{0, g.nx, 0, g.ny + 1, 0, g.nz}, KFacePDY{fl, w.dhy.view(), w.loy.view()});
          fz.add(Range3D{0, g.nx, 0, g.ny, 0, g.nz + 1}, KFacePDZ{fl, w.dhz.view(), w.loz.view()});
        }
        fx.run();
        fy.run();
        fz.run();
      }
      auto faces = [&](const Work& w) {
        return FacesPD{w.dhx.cview(), w.lox.cview(), w.dhy.cview(), w.loy.cview(), w.dhz.cview(), w.loz.cview()};
      };
      LaunchBatch<KScalePD, kBatch> sc;
      std::vector<Field3*> scales;
      for (int m = 0; m < nb; ++m) {
        Work& w = *work_[m];
        sc.add(mr, KScalePD{fluxes(species(m)), faces(w), s.mu_1.cview(), b.mub.view(), w.scale.view()});
        scales.push_back(&w.scale);
      }
      sc.run();
      fill_halos(scales, g);
      LaunchBatch<KUpdatePD, kBatch> up;
      for (int m = 0; m < nb; ++m) {
        const int n = species(m);
        const Work& w = *work_[m];
        up.add(mr, KUpdatePD{fluxes(n), faces(w), w.scale.cview(), s.q_1[n].cview(), s.mu_1.cview(), s.mu_2.cview(),
                             b.mub.view(), w.q_new.view(), dt_rk, ften(n), ften(n).p != nullptr});
      }
      up.run();
    }
    std::vector<int> fd;  // species of this batch for flow_dep_bdy
    for (int m = 0; m < nb; ++m) {
      const int n = species(m);
      if (n == s.P_QV && fsink != nullptr && vmdump::on()) vmdump::dump("vm_qv_rk", rk, {&s.q_1[n], fsink});
      swap(s.q_2[n], work_[m]->q_new);
      if (n == s.P_TKE && dd_) dd_->bound_tke(s.q_2[n]);
      fd.push_back(n);
    }
    if (bc_) bc_->flow_dep(fd, s);  // the species without boundary data (have_bcs_moist/scalar = .false.)
  }
}

void ScalarTransport::run_stage(State& s, Tendencies& t, const BaseState& b, int rk, real dt_rk) {
  const Grid& g = g_;
  const int ns = s.n_scalars();
  if (ns == 0) return;
  if (int(t.q_tendf.size()) != ns) throw std::runtime_error("Tendencies: q_tendf not registered for all scalars");
  // WRF's species loop (rk_scalar_tend, then rk_update_scalar(_pd) of each species) as the same steps per kernel for
  // all species: every step of a species reads and writes only that species' fields (and the common masses and
  // fluxes), so each species gets the same values in the same order of its own operations
  std::vector<int> all(ns), pd_n, plain_n;
  for (int n = 0; n < ns; ++n) {
    all[n] = n;
    const bool moist = n < s.n_moist;
    const bool is_tke = n == s.P_TKE;
    const bool pd = rk == 3 && (is_tke ? p_.tke_adv_opt : (moist ? p_.moist_adv_opt : p_.scalar_adv_opt)) == 1;
    (pd ? pd_n : plain_n).push_back(n);
  }
  if (rk == 1) {  // rk_scalar_tend, rk_step == 1: diffusion into the held tendency (diff_opt = 1)
    // q_tendf was zeroed by zero_tendf and holds the physics tendencies (update_phy_ten)
    if (dd_)
      for (int n = 0; n < ns; ++n) dd_->scalar(s, b, n, t.q_tendf[n], dt_rk);  // diff_opt 1 / sixth order
    if (bc_) bc_->relax_scalars(all, t, s, b);  // relax_bdy_scalar (rk 1): scalars with boundary data
  }
  if (bc_) bc_->spec_scalars(all, t);  // spec_bdy_scalar (every stage)
  if (!plain_n.empty()) update(s, t, b, plain_n, false, rk, dt_rk);
  if (!pd_n.empty()) update(s, t, b, pd_n, true, rk, dt_rk);
  std::vector<Field3*> qs;  // halos of the new scalars, batched (one launch per 8 fields)
  for (int n = 0; n < ns; ++n) qs.push_back(&s.q_2[n]);
  fill_halos(qs, g);
}

}  // namespace cppwrf::dyn
