#include "dyn/small_step.hpp"

#include <cmath>

#include "core/backend.hpp"
#include "core/constants.hpp"
#include "core/halo.hpp"
#include "dyn/kernels_common.hpp"
#include "dyn/large_step.hpp"
#include "dyn/lateral_bc.hpp"

// Discrete forms follow WRF dyn_em/module_small_step_em.F (small_step_prep, calc_p_rho,
// calc_coef_w, advance_uv, advance_mu_t, advance_w, sumflux, small_step_finish) and
// calc_p_rho_phi (module_big_step_utilities_em.F, hypsometric_opt = 1 as ideal.exe writes it,
// or 2 for real-data input).
// C++ k = Fortran k - 1. Map factors = 1, flat terrain (ht = 0). Moisture enters through
// theta_m (t_2) and the calc_cq coefficients cqu, cqv and cqw (= cq1 after pg_buoy_w).
// Variables named *pp below are the acoustic perturbations (WRF keeps them in *_2, p, al).

namespace cppwrf::dyn {
namespace {

using V3 = View3<real>;
using C3 = View3<const real>;
using V2 = View2<real>;
using C2 = View2<const real>;
using B3 = BView3<const real>;
using B2 = BView2<const real>;
constexpr real Gr = real(cst::g);
constexpr real H = real(0.5);
constexpr real T0 = real(cst::t0);

// Specified lateral boundaries: WRF advance_uv / advance_mu_t / advance_w leave the outermost
// row (spec zone) to spec_bdyupdate*, zero_grad_bdy. True for mass point (i,j) in that row.
struct Edge {
  bool x = false, y = false; int nx = 0, ny = 0;
  CPPWRF_HD bool mass(int i, int j) const { return (x && (i == 0 || i == nx - 1)) || (y && (j == 0 || j == ny - 1)); }
  CPPWRF_HD bool uface(int i, int j) const { return (x && (i <= 0 || i >= nx)) || (y && (j == 0 || j == ny - 1)); }
  CPPWRF_HD bool vface(int i, int j) const { return (x && (i == 0 || i == nx - 1)) || (y && (j <= 0 || j >= ny)); }
};

// ------------------------------------------------------------------ small_step_prep
// Saves the stage state and forms the coupled perturbations from it. Full masses:
// mut/muu/muv = stage (from rk_step_prep), muts = mub + mu_1 (time t).
struct KPrep {
  V3 u, v, w, t, ph, us, vs, ws, ts, phs, ww1, c2a;
  C3 u1, v1, w1, t1, ph1, ww, p, alt;
  V2 mu, mus, muts, muus, muvs;
  C2 mu1, mut, muu, muv;
  B2 mub; B3 pb; VCoef c;
  int nx, ny, nyv, nz;
  B2 msfuy, msfvx_inv, msfty;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    const real mts = mub(i, j) + mu1(i, j);
    const real mus_u = H * (mts + mub(i - 1, j) + mu1(i - 1, j));
    const real mus_v = H * (mts + mub(i, j - 1) + mu1(i, j - 1));
    if (k == 0) {
      if (i < nx && j < ny) { muts(i, j) = mts; mus(i, j) = mu(i, j); mu(i, j) = mu1(i, j) - mu(i, j); }
      if (j < ny) muus(i, j) = mus_u;
      if (i < nx && j < nyv) muvs(i, j) = mus_v;
    }
    if (j < ny && k < nz) {
      const real ui = u(i, j, k);
      us(i, j, k) = ui;
      u(i, j, k) = ((c.c1h(k) * mus_u + c.c2h(k)) * u1(i, j, k) - (c.c1h(k) * muu(i, j) + c.c2h(k)) * ui) / msfuy(i, j);
    }
    if (i < nx && j < nyv && k < nz) {
      const real vi = v(i, j, k);
      vs(i, j, k) = vi;
      v(i, j, k) = ((c.c1h(k) * mus_v + c.c2h(k)) * v1(i, j, k) - (c.c1h(k) * muv(i, j) + c.c2h(k)) * vi) * msfvx_inv(i, j);
    }
    if (i < nx && j < ny) {
      if (k < nz) {
        const real ti = t(i, j, k);
        ts(i, j, k) = ti;
        t(i, j, k) = (c.c1h(k) * mts + c.c2h(k)) * t1(i, j, k) - (c.c1h(k) * mut(i, j) + c.c2h(k)) * ti;
        c2a(i, j, k) = real(cst::cpovcv) * (pb(i, j, k) + p(i, j, k)) / alt(i, j, k);
      }
      const real wi = w(i, j, k), phi = ph(i, j, k);
      ws(i, j, k) = wi;
      w(i, j, k) = ((c.c1f(k) * mts + c.c2f(k)) * w1(i, j, k) - (c.c1f(k) * mut(i, j) + c.c2f(k)) * wi) / msfty(i, j);
      phs(i, j, k) = phi;
      ph(i, j, k) = ph1(i, j, k) - phi;
      ww1(i, j, k) = (k == 0 || k == nz) ? real(0) : ww(i, j, k);
    }
  }
};

// ------------------------------------------------------------------ calc_p_rho
// Acoustic p'' and al'' from mu'', phi'', Theta''. step > 0 adds the 3D divergence damping.
CPPWRF_HD inline void p_rho_point(int i, int j, int k, V3 al, V3 p, V3 pm1, C3 ph, C3 alt, C3 t, C3 ts, C3 c2a, C2 mu,
                                  C2 muts, const VCoef& c, real smdiv, bool damp) {
  const real m = c.c1h(k) * muts(i, j) + c.c2h(k);
  const real a = -(alt(i, j, k) * (c.c1h(k) * mu(i, j)) + c.rdnw(k) * (ph(i, j, k + 1) - ph(i, j, k))) / m;
  al(i, j, k) = a;
  real pp = c2a(i, j, k) * (alt(i, j, k) * (t(i, j, k) - c.c1h(k) * mu(i, j) * ts(i, j, k)) / (m * (T0 + ts(i, j, k))) - a);
  if (damp) {
    const real ptmp = pp;
    pp = pp + smdiv * (pp - pm1(i, j, k));
    pm1(i, j, k) = ptmp;
  } else {
    pm1(i, j, k) = pp;
  }
  p(i, j, k) = pp;
}

struct KPRho {
  V3 al, p, pm1; C3 ph, alt, t, ts, c2a; C2 mu, muts; VCoef c;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    p_rho_point(i, j, k, al, p, pm1, ph, alt, t, ts, c2a, mu, muts, c, real(0), false);
  }
};

// ------------------------------------------------------------------ calc_coef_w
struct KCoefW {
  V3 a, alpha, gamma; C3 c2a; C2 mut; VCoef c; kern::MoistCoef mc; real cof; real lid; int nz;
  CPPWRF_HD real mh(int k, real m) const { return c.c1h(k) * m + c.c2h(k); }
  CPPWRF_HD real mf(int k, real m) const { return c.c1f(k) * m + c.c2f(k); }
  CPPWRF_HD void operator()(int i, int j) const {
    const real m = mut(i, j);
    a(i, j, 0) = 0; alpha(i, j, 0) = 0; gamma(i, j, 0) = 0;
    a(i, j, 1) = 0;
    a(i, j, nz) = real(-2) * cof * (c.rdnw(nz - 1) * c.rdnw(nz - 1)) * c2a(i, j, nz - 1) * lid /
                  (mh(nz - 1, m) * mf(nz - 1, m));
    for (int k = 2; k < nz; ++k)  // WRF quirk kept: c1f at the level below
      a(i, j, k) = -mc.cq1(i, j, k) * cof * c.rdn(k) * c.rdnw(k - 1) * c2a(i, j, k - 1) / (mh(k - 1, m) * mf(k - 1, m));
    for (int k = 1; k < nz; ++k) {
      const real cq1 = mc.cq1(i, j, k);
      const real b = real(1) + cq1 * cof * c.rdn(k) *
                                   (c.rdnw(k) * c2a(i, j, k) / (mh(k, m) * mf(k, m)) +
                                    c.rdnw(k - 1) * c2a(i, j, k - 1) / (mh(k - 1, m) * mf(k, m)));
      const real cc = -cq1 * cof * c.rdn(k) * c.rdnw(k) * c2a(i, j, k) / (mh(k, m) * mf(k + 1, m));
      alpha(i, j, k) = real(1) / (b - a(i, j, k) * gamma(i, j, k - 1));
      gamma(i, j, k) = cc * alpha(i, j, k);
    }
    const real b = real(1) + real(2) * cof * (c.rdnw(nz - 1) * c.rdnw(nz - 1)) * c2a(i, j, nz - 1) /
                                 (mh(nz - 1, m) * mf(nz, m));
    alpha(i, j, nz) = real(1) / (b - a(i, j, nz) * gamma(i, j, nz - 1));
    gamma(i, j, nz) = 0;
  }
};

// ------------------------------------------------------------------ advance_uv
struct KUV {
  V3 u, v; C3 ru_tend, rv_tend, p, al, alt, ph, php; C2 mu, muu, muv, mudf; B3 pb; VCoef c; kern::MoistCoef mc;
  real dts, rdx, rdy, emdx, emdy; int nz, ui0, ui1, vj0, vj1; bool do_y, top_lid;
  B2 msfux, msfuy, msfvx, msfvy, msfvx_inv; Edge ed;
  template <int DI, int DJ> CPPWRF_HD real ps(int i, int j, int k) const { return p(i, j, k) + p(i - DI, j - DJ, k); }
  template <int DI, int DJ> CPPWRF_HD real dpn(int i, int j, int k) const {
    if (k == 0) return H * (c.cf1 * ps<DI, DJ>(i, j, 0) + c.cf2 * ps<DI, DJ>(i, j, 1) + c.cf3 * ps<DI, DJ>(i, j, 2));
    if (k == nz)  // WRF advance_uv uses cf1..cf3 at the lid too
      return top_lid ? H * (c.cf1 * ps<DI, DJ>(i, j, nz - 1) + c.cf2 * ps<DI, DJ>(i, j, nz - 2) +
                            c.cf3 * ps<DI, DJ>(i, j, nz - 3))
                     : real(0);
    return H * (c.fnm(k) * ps<DI, DJ>(i, j, k) + c.fnp(k) * ps<DI, DJ>(i, j, k - 1));
  }
  template <int DI, int DJ> CPPWRF_HD real dp(int i, int j, int k, real mface, real rd, real rm) const {
    real d = rm * H * rd * (c.c1h(k) * mface + c.c2h(k)) *
             ((ph(i, j, k + 1) - ph(i - DI, j - DJ, k + 1) + ph(i, j, k) - ph(i - DI, j - DJ, k)) +
              (alt(i, j, k) + alt(i - DI, j - DJ, k)) * (p(i, j, k) - p(i - DI, j - DJ, k)) +
              (al(i, j, k) + al(i - DI, j - DJ, k)) * (pb(i, j, k) - pb(i - DI, j - DJ, k)));
    d += rm * rd * (php(i, j, k) - php(i - DI, j - DJ, k)) *
         (c.rdnw(k) * (dpn<DI, DJ>(i, j, k + 1) - dpn<DI, DJ>(i, j, k)) -
          H * (c.c1h(k) * mu(i - DI, j - DJ) + c.c1h(k) * mu(i, j)));
    return d;
  }
  CPPWRF_HD void operator()(int i, int j, int k) const {
    if (i >= ui0 && i < ui1 && !ed.uface(i, j)) {
      const real mudf_x = -emdx * (mudf(i, j) - mudf(i - 1, j)) / msfuy(i, j);
      u(i, j, k) = u(i, j, k) + dts * ru_tend(i, j, k) -
                   dts * mc.cqu(i, j, k) * dp<1, 0>(i, j, k, muu(i, j), rdx, msfux(i, j) / msfuy(i, j)) + c.c1h(k) * mudf_x;
    }
    if (do_y && j >= vj0 && j < vj1 && !ed.vface(i, j)) {
      const real mudf_y = -emdy * (mudf(i, j) - mudf(i, j - 1)) * msfvx_inv(i, j);
      v(i, j, k) = v(i, j, k) + dts * rv_tend(i, j, k) -
                   dts * mc.cqv(i, j, k) * dp<0, 1>(i, j, k, muv(i, j), rdy, msfvy(i, j) / msfvx(i, j)) + c.c1h(k) * mudf_y;
    }
  }
};

// ------------------------------------------------------------------ advance_mu_t + sumflux
// Split into point-parallel phases (same expressions and summation order as one column loop).
struct KMuTDiv {  // dvdxi at every mass point
  C3 u, v, us, vs; C2 muu, muv; VCoef c; real rdx, rdy; bool do_y;
  B2 msftx, msfty, msfuy, msfvx_inv; V3 dv;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    const real fx = (u(i + 1, j, k) + (c.c1h(k) * muu(i + 1, j) + c.c2h(k)) * us(i + 1, j, k) / msfuy(i + 1, j)) -
                    (u(i, j, k) + (c.c1h(k) * muu(i, j) + c.c2h(k)) * us(i, j, k) / msfuy(i, j));
    real d = rdx * fx;
    if (do_y)
      d += rdy * ((v(i, j + 1, k) + (c.c1h(k) * muv(i, j + 1) + c.c2h(k)) * vs(i, j + 1, k) * msfvx_inv(i, j + 1)) -
                  (v(i, j, k) + (c.c1h(k) * muv(i, j) + c.c2h(k)) * vs(i, j, k) * msfvx_inv(i, j)));
    dv(i, j, k) = msftx(i, j) * msfty(i, j) * d;
  }
};

struct KMuTCol {  // column: mu, mudf, muts, muave; omega'' (ww) and its time sum ww_m
  C3 dv, ww1; V3 ww, ww_m; V2 mu, mudf, muts, muave; C2 mut, mu_tend; B2 msfty; VCoef c;
  real dts, epssm; int nz, iter, nsteps; Edge ed;
  CPPWRF_HD void operator()(int i, int j) const {
    const bool first = iter == 1, last = iter == nsteps;
    const real rn = real(nsteps);  // WRF sumflux divides by number_of_small_timesteps
    if (ed.mass(i, j)) {  // spec zone: no advance_mu_t; sumflux still accumulates ww
      for (int k = 0; k <= nz; ++k) {
        const real sw = (first ? real(0) : ww_m(i, j, k)) + ww(i, j, k);
        ww_m(i, j, k) = last ? sw / rn + ww1(i, j, k) : sw;
      }
      return;
    }
    real dmdt = 0;
    for (int k = 0; k < nz; ++k) dmdt += c.dnw(k) * dv(i, j, k);
    const real mtend = mu_tend(i, j);
    const real mu_old = mu(i, j);
    const real mu_new = mu_old + dts * (dmdt + mtend);
    mu(i, j) = mu_new;
    mudf(i, j) = dmdt + mtend;
    muts(i, j) = mut(i, j) + mu_new;
    muave(i, j) = H * ((real(1) + epssm) * mu_new + (real(1) - epssm) * mu_old);
    auto acc_ww = [&](int k, real wpp) {
      const real sw = (first ? real(0) : ww_m(i, j, k)) + wpp;
      ww_m(i, j, k) = last ? sw / rn + ww1(i, j, k) : sw;
    };
    real ww_full = 0;  // running integral (full flux); ww'' at the bottom is 0
    ww(i, j, 0) = 0;
    acc_ww(0, real(0));
    for (int k = 0; k < nz; ++k) {
      real wpp_hi = 0;
      if (k + 1 < nz) {
        ww_full = ww_full - c.dnw(k) * (c.c1h(k) * dmdt + dv(i, j, k) + c.c1h(k) * mtend) / msfty(i, j);
        wpp_hi = ww_full - ww1(i, j, k + 1);
        ww(i, j, k + 1) = wpp_hi;
      } else {
        ww(i, j, nz) = 0;
      }
      acc_ww(k + 1, wpp_hi);
    }
  }
};

struct KMuTPoint {  // sumflux of ru_m, rv_m and the theta'' update, at every mass point
  C3 u, v, us, vs, ts, t_tend, ww; V3 t, t_ave, ru_m, rv_m;
  C2 muu, muv; VCoef c; real dts, rdx, rdy; int nz, iter, nsteps; bool do_y;
  B2 msftx, msfty, msfuy, msfvx_inv; Edge ed;
  CPPWRF_HD real wdtn(int i, int j, int k) const {  // ww'' (ts at w level), 0 at the bottom and top
    if (k == 0 || k == nz) return 0;
    return ww(i, j, k) * (c.fnm(k) * ts(i, j, k) + c.fnp(k) * ts(i, j, k - 1));
  }
  CPPWRF_HD void operator()(int i, int j, int k) const {
    const bool first = iter == 1, last = iter == nsteps;
    const real rn = real(nsteps);
    const real su = (first ? real(0) : ru_m(i, j, k)) + u(i, j, k);
    ru_m(i, j, k) = last ? su / rn + (c.c1h(k) * muu(i, j) + c.c2h(k)) * us(i, j, k) / msfuy(i, j) : su;
    if (do_y) {
      const real sv = (first ? real(0) : rv_m(i, j, k)) + v(i, j, k);
      rv_m(i, j, k) = last ? sv / rn + (c.c1h(k) * muv(i, j) + c.c2h(k)) * vs(i, j, k) * msfvx_inv(i, j) : sv;
    }
    const real tk = t(i, j, k);
    t_ave(i, j, k) = tk;
    if (ed.mass(i, j)) return;  // spec zone: theta'' by spec_bdyupdate
    real hflux = H * rdx * (u(i + 1, j, k) * (ts(i + 1, j, k) + ts(i, j, k)) - u(i, j, k) * (ts(i, j, k) + ts(i - 1, j, k)));
    if (do_y)
      hflux += H * rdy * (v(i, j + 1, k) * (ts(i, j + 1, k) + ts(i, j, k)) - v(i, j, k) * (ts(i, j, k) + ts(i, j - 1, k)));
    const real tn = tk + msfty(i, j) * dts * t_tend(i, j, k);
    t(i, j, k) = tn - dts * msfty(i, j) * (msftx(i, j) * hflux + c.rdnw(k) * (wdtn(i, j, k + 1) - wdtn(i, j, k)));
  }
};

// ------------------------------------------------------------------ advance_w + calc_p_rho
// Split into point-parallel phases (GPU: nz times more threads than one thread per column);
// every expression is the one of the original single column loop, so results are identical.
struct KWCommon {
  C3 w, ph, rw_tend, ph_tend, ww, phs, t, t_ave, ts, c2a, alt;
  C2 mut, muts, muave; B3 phb; VCoef c; kern::MoistCoef mc;
  real dts, epssm; int nz; bool top_lid;
  B2 msfty;
  CPPWRF_HD real mh(int k, real m) const { return c.c1h(k) * m + c.c2h(k); }
  CPPWRF_HD real mf(int k, real m) const { return c.c1f(k) * m + c.c2f(k); }
  // t_2ave of WRF advance_w (time-averaged, normalised theta'')
  CPPWRF_HD real t2ave(int i, int j, int k) const {
    const real tav = H * ((real(1) + epssm) * t(i, j, k) + (real(1) - epssm) * t_ave(i, j, k));
    return (tav + c.c1h(k) * muave(i, j) * T0) / (mh(k, muts(i, j)) * (T0 + ts(i, j, k)));
  }
  // vertical advection of the base+stage geopotential by omega'' (phi_adv_z = 1), level m = 1..nz
  CPPWRF_HD real wdwn(int i, int j, int m) const {
    return H * (ww(i, j, m) + ww(i, j, m - 1)) * c.rdnw(m - 1) *
           (phs(i, j, m) - phs(i, j, m - 1) + phb(i, j, m) - phb(i, j, m - 1));
  }
  // explicit part of phi'' (WRF rhs), with the old w'' and phi''
  CPPWRF_HD real rhs(int i, int j, int k) const {
    if (k == 0) return 0;
    if (top_lid && k == nz) return 0;
    real r = dts * (ph_tend(i, j, k) + H * Gr * (real(1) - epssm) * w(i, j, k));
    if (k < nz) r -= dts * (c.fnm(k) * wdwn(i, j, k + 1) + c.fnp(k) * wdwn(i, j, k));
    return ph(i, j, k) + msfty(i, j) * r / mf(k, mut(i, j));
  }
};

// phase 1a (w levels 0..nz): rhs -> R, t_2ave -> T2
struct KWRhs {
  KWCommon q; V3 R, T2; Edge ed;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    if (ed.mass(i, j)) return;
    R(i, j, k) = q.rhs(i, j, k);
    if (k < q.nz) T2(i, j, k) = q.t2ave(i, j, k);
  }
};

// phase 1b (w levels 0..nz): explicit w'' (in place; reads only old ph'' and R)
struct KWExpl {
  KWCommon q; C3 R, T2;
  C3 us2, vs2; B2 ht, msftx; real rdx, rdy; bool terrain, do_y;  // lower boundary w = u.grad(ht)
  V3 wout; Edge ed;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    if (ed.mass(i, j)) return;  // spec zone: zero_grad_bdy
    const VCoef& c = q.c;
    if (k == 0) {
      if (terrain) {  // WRF advance_w: w''(surface) from the small-step u'', v'' and the terrain slope
        real ws = msftx(i, j) * H * rdx *
                  ((ht(i + 1, j) - ht(i, j)) * (c.cf1 * us2(i + 1, j, 0) + c.cf2 * us2(i + 1, j, 1) + c.cf3 * us2(i + 1, j, 2)) +
                   (ht(i, j) - ht(i - 1, j)) * (c.cf1 * us2(i, j, 0) + c.cf2 * us2(i, j, 1) + c.cf3 * us2(i, j, 2)));
        if (do_y)
          ws = q.msfty(i, j) * H * rdy *
                   ((ht(i, j + 1) - ht(i, j)) * (c.cf1 * vs2(i, j + 1, 0) + c.cf2 * vs2(i, j + 1, 1) + c.cf3 * vs2(i, j + 1, 2)) +
                    (ht(i, j) - ht(i, j - 1)) * (c.cf1 * vs2(i, j, 0) + c.cf2 * vs2(i, j, 1) + c.cf3 * vs2(i, j, 2))) +
               ws;
        wout(i, j, 0) = ws;
      } else {
        wout(i, j, 0) = 0;
      }
      return;
    }
    const real m = q.mut(i, j), mav = q.muave(i, j), dts = q.dts;
    const real ope = real(1) + q.epssm, ome = real(1) - q.epssm;
    const real msft_inv = real(1) / q.msfty(i, j);
    const real r_k = R(i, j, k), r_km1 = R(i, j, k - 1), ph_k = q.ph(i, j, k), ph_km1 = q.ph(i, j, k - 1);
    if (k < q.nz) {
      const real r_kp1 = R(i, j, k + 1), ph_kp1 = q.ph(i, j, k + 1);
      wout(i, j, k) = q.w(i, j, k) + dts * q.rw_tend(i, j, k) +
                      msft_inv * q.mc.cq1(i, j, k) *
                          (H * dts * Gr * c.rdn(k) *
                           (q.c2a(i, j, k) * c.rdnw(k) / q.mh(k, m) * (ope * (r_kp1 - r_k) + ome * (ph_kp1 - ph_k)) -
                            q.c2a(i, j, k - 1) * c.rdnw(k - 1) / q.mh(k - 1, m) * (ope * (r_k - r_km1) + ome * (ph_k - ph_km1)))) +
                      dts * Gr * msft_inv *
                          (c.rdn(k) * (q.c2a(i, j, k) * q.alt(i, j, k) * T2(i, j, k) -
                                       q.c2a(i, j, k - 1) * q.alt(i, j, k - 1) * T2(i, j, k - 1)) -
                           c.c1f(k) * mav);
    } else {
      const int km = q.nz - 1;
      real wn = q.w(i, j, k) + dts * q.rw_tend(i, j, k) +
                msft_inv * (-H * dts * Gr / q.mh(km, m) * (c.rdnw(km) * c.rdnw(km)) * real(2) * q.c2a(i, j, km) *
                                (ope * (r_k - r_km1) + ome * (ph_k - ph_km1)) -
                            dts * Gr * (real(2) * c.rdnw(km) * q.c2a(i, j, km) * q.alt(i, j, km) * T2(i, j, km) + c.c1f(k) * mav));
      if (q.top_lid) wn = 0;
      wout(i, j, k) = wn;
    }
  }
};

// phase 2 (columns): tridiagonal solve (calc_coef_w coefficients), damp_opt = 3, phi'' = R + implicit w part
struct KWSolve {
  V3 w, ph; C3 R, a, alpha, gamma, phs, ws; C2 mut, muts; B3 phb; B2 msfty; VCoef c;
  real dts, epssm, dampmag, hdepth; int nz; Edge ed;
  CPPWRF_HD real mf(int k, real m) const { return c.c1f(k) * m + c.c2f(k); }
  CPPWRF_HD void operator()(int i, int j) const {
    if (ed.mass(i, j)) return;  // spec zone: spec_bdyupdate_ph
    const real m = mut(i, j), mts = muts(i, j), ope = real(1) + epssm;
    for (int k = 1; k <= nz; ++k) w(i, j, k) = (w(i, j, k) - a(i, j, k) * w(i, j, k - 1)) * alpha(i, j, k);
    for (int k = nz - 1; k >= 1; --k) w(i, j, k) = w(i, j, k) - gamma(i, j, k) * w(i, j, k + 1);
    if (dampmag != real(0)) {  // damp_opt = 3: implicit Rayleigh damping toward w = 0 near the top
      const real pi = real(4) * atan(real(1));
      const real htop = (phs(i, j, nz) + phb(i, j, nz)) / Gr, hbot = htop - hdepth;
      for (int k = nz; k >= 1; --k) {
        const real hk = (phs(i, j, k) + phb(i, j, k)) / Gr;
        real dampwt = 0;
        if (hk >= hbot) {
          const real sn = sin(H * pi * (hk - hbot) / hdepth);
          dampwt = dampmag * sn * sn;
        }
        w(i, j, k) = (w(i, j, k) - dampwt * mf(k, m) * ws(i, j, k)) / (real(1) + dampwt);
      }
    }
    for (int k = nz; k >= 1; --k) ph(i, j, k) = R(i, j, k) + msfty(i, j) * H * dts * Gr * ope * w(i, j, k) / mf(k, mts);
  }
};

// phase 3 (mass points): calc_p_rho with divergence damping
struct KWPRho {
  V3 al, p, pm1; C3 ph, alt, t, ts, c2a; C2 mu, muts; VCoef c; real smdiv;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    p_rho_point(i, j, k, al, p, pm1, ph, alt, t, ts, c2a, mu, muts, c, smdiv, true);
  }
};

// ------------------------------------------------------------------ small_step_finish
struct KFinish {
  V3 u, v, w, t, ph, ww; C3 us, vs, ws, ts, phs, ww1; C2 muu, muv, mut, muts; VCoef c;
  int nx, ny, nyv, nz; bool do_y;
  C3 hd; real hd_dt;  // rk_step == 3: dts*number_of_small_timesteps (0 otherwise or when dry)
  B2 msfuy, msfvx, msfty;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    if (j < ny && k < nz) {
      const real mus_u = H * (muts(i, j) + muts(i - 1, j));
      u(i, j, k) = (msfuy(i, j) * u(i, j, k) + us(i, j, k) * (c.c1h(k) * muu(i, j) + c.c2h(k))) / (c.c1h(k) * mus_u + c.c2h(k));
    }
    if (i < nx && j < nyv && k < nz) {  // also in 2D, as WRF (v'' is then constant)
      const real mus_v = H * (muts(i, j) + muts(i, j - 1));
      v(i, j, k) = (msfvx(i, j) * v(i, j, k) + vs(i, j, k) * (c.c1h(k) * muv(i, j) + c.c2h(k))) / (c.c1h(k) * mus_v + c.c2h(k));
    }
    if (i < nx && j < ny) {
      const real m = mut(i, j), mts = muts(i, j);
      w(i, j, k) = (msfty(i, j) * w(i, j, k) + ws(i, j, k) * (c.c1f(k) * m + c.c2f(k))) / (c.c1f(k) * mts + c.c2f(k));
      ph(i, j, k) = ph(i, j, k) + phs(i, j, k);
      ww(i, j, k) = ww(i, j, k) + ww1(i, j, k);
      if (k < nz)
        t(i, j, k) = (t(i, j, k) - (hd_dt != real(0) ? hd_dt * (c.c1h(k) * m + c.c2h(k)) * hd(i, j, k) : real(0)) +
                      ts(i, j, k) * (c.c1h(k) * m + c.c2h(k))) /
                     (c.c1h(k) * mts + c.c2h(k));
    }
  }
};

struct KFinishMu {
  V2 mu, muus, muvs; C2 mus, muts; int nx, ny, nyv;
  CPPWRF_HD void operator()(int i, int j, int) const {
    if (i < nx && j < ny) mu(i, j) = mu(i, j) + mus(i, j);
    if (j < ny) muus(i, j) = H * (muts(i, j) + muts(i - 1, j));
    if (i < nx && j < nyv) muvs(i, j) = H * (muts(i, j) + muts(i, j - 1));
  }
};

// calc_p_rho_phi, non-hydrostatic, dry, hypsometric_opt = 1
struct KPRhoPhi {
  V3 al, p; C3 ph, t; C2 mu, muts; B3 alb, pb; VCoef c;
  int hypso; B3 phb; View1<const real> c3h, c4h, c3f, c4f; real ptop;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    real a;
    if (hypso == 2) {  // log-linear in pressure
      const real m = muts(i, j);
      const real pfu = c3f(k + 1) * m + c4f(k + 1) + ptop, pfd = c3f(k) * m + c4f(k) + ptop;
      const real phm = c3h(k) * m + c4h(k) + ptop;
      a = (ph(i, j, k + 1) - ph(i, j, k) + phb(i, j, k + 1) - phb(i, j, k)) / phm / rlog(pfd / pfu) - alb(i, j, k);
    } else {
      a = real(-1) / (c.c1h(k) * muts(i, j) + c.c2h(k)) *
          (alb(i, j, k) * (c.c1h(k) * mu(i, j)) + c.rdnw(k) * (ph(i, j, k + 1) - ph(i, j, k)));
    }
    al(i, j, k) = a;
    // in the working precision, as WRF (REAL): FP64 pow on consumer GPUs runs at 1/32 rate
    const real base = real(cst::r_d) * (T0 + t(i, j, k)) / (real(cst::p1000mb) * (a + alb(i, j, k)));
    p(i, j, k) = real(cst::p1000mb) * rpow(base, real(cst::cpovcv)) - pb(i, j, k);
  }
};

}  // namespace

SmallStep::SmallStep(const Grid& g, const SmallStepParams& p)
    : g_(g), p_(p),
      u_save(g, Stagger::U, "u_save"), v_save(g, Stagger::V, "v_save"), w_save(g, Stagger::W, "w_save"),
      t_save(g, Stagger::M, "t_save"), ph_save(g, Stagger::W, "ph_save"), ww1(g, Stagger::W, "ww1"),
      c2a(g, Stagger::M, "c2a"), pm1(g, Stagger::M, "pm1"), t_2save(g, Stagger::M, "t_2save"),
      a(g, Stagger::W, "a"), alpha(g, Stagger::W, "alpha"), gamma(g, Stagger::W, "gamma"),
      mu_save(g, Stagger::M, "mu_save"), mudf(g, Stagger::M, "mudf"), muave(g, Stagger::M, "muave"),
      w_rhs(g, Stagger::W, "w_rhs"), t2a(g, Stagger::M, "t_2ave"), dvdxi(g, Stagger::M, "dvdxi") {}

void SmallStep::run_stage(State& s, const Tendencies& tend, const BaseState& b, real dt_stage, real dts, int n_small,
                          int rk_step) {
  (void)dt_stage;  // used by WRF only for h_diabatic (physics)
  const Grid& g = g_;
  const int nx = g.nx, ny = g.ny, nz = g.nz, nyv = g.nys(Stagger::V);
  const bool do_y = !g.two_d();
  const VCoef c = VCoef::from(b);
  const kern::MoistCoef mc = moist_coef(s);  // calc_cq from the stage state (rk_step_prep)
  const Edge ed{g.bc_x == BC::Specified, g.bc_y == BC::Specified && do_y, nx, ny};
  stage_ = rk_step > 0 ? rk_step : stage_ % 3 + 1;

  if (stage_ == 1) {  // WRF small_step_prep (rk_step == 1): reset the external-mode filter
    V2 md = mudf.view();
    parallel_for(Range3D{-Grid::halo, nx + Grid::halo, -g.hy(), ny + g.hy(), 0, 1},
                 [=] CPPWRF_HD(int i, int j, int) { md(i, j) = 0; });
  }

  parallel_for(Range3D{0, nx + 1, 0, nyv, 0, nz + 1},
               KPrep{s.u_2.view(), s.v_2.view(), s.w_2.view(), s.t_2.view(), s.ph_2.view(), u_save.view(),
                     v_save.view(), w_save.view(), t_save.view(), ph_save.view(), ww1.view(), c2a.view(),
                     s.u_1.cview(), s.v_1.cview(), s.w_1.cview(), s.t_1.cview(), s.ph_1.cview(), s.ww.cview(),
                     s.p.cview(), s.alt.cview(), s.mu_2.view(), mu_save.view(), s.muts.view(), s.muus.view(),
                     s.muvs.view(), s.mu_1.cview(), s.mut.cview(), s.muu.cview(), s.muv.cview(), b.mub.view(),
                     b.pb.view(), c, nx, ny, nyv, nz, b.msfuy.view(), b.msfvx_inv.view(), b.msfty.view()});
  fill_halos({&s.u_2, &s.v_2, &s.ph_2, &t_save, &u_save, &v_save}, g);
  fill_halos({&s.mu_2, &s.muts}, g);

  parallel_for(g.mass_cells(), KPRho{s.al.view(), s.p.view(), pm1.view(), s.ph_2.cview(), s.alt.cview(),
                                     s.t_2.cview(), t_save.cview(), c2a.cview(), s.mu_2.cview(), s.muts.cview(), c});
  fill_halos({&s.p, &s.al}, g);

  const real cof = real(std::pow(0.5 * double(dts) * cst::g * (1.0 + double(p_.epssm)), 2));
  column_for(g.columns(), KCoefW{a.view(), alpha.view(), gamma.view(), c2a.cview(), s.mut.cview(), c, mc, cof,
                                 p_.top_lid ? real(0) : real(1), nz});

  const Range3D uv_range{0, nx, 0, ny, 0, nz};
  for (int it = 1; it <= n_small; ++it) {
    parallel_for(uv_range, KUV{s.u_2.view(), s.v_2.view(), tend.ru_tend.cview(), tend.rv_tend.cview(), s.p.cview(),
                               s.al.cview(), s.alt.cview(), s.ph_2.cview(), s.php.cview(), s.mu_2.cview(),
                               s.muu.cview(), s.muv.cview(), mudf.cview(), b.pb.view(), c, mc, dts, real(b.rdx),
                               real(b.rdy), p_.emdiv * real(g.dx), do_y ? p_.emdiv * real(g.dy) : real(0), nz,
                               g.u_i0(), g.u_i1(), g.v_j0(), g.v_j1(), do_y, p_.top_lid, b.msfux.view(),
                               b.msfuy.view(), b.msfvx.view(), b.msfvy.view(), b.msfvx_inv.view(), ed});
    if (bc_) bc_->update_uv(s, tend, dts);  // spec_bdyupdate u, v
    if (do_y) fill_halos({&s.u_2, &s.v_2}, g);
    else fill_halo(s.u_2, g);

    parallel_for(g.mass_cells(), KMuTDiv{s.u_2.cview(), s.v_2.cview(), u_save.cview(), v_save.cview(), s.muu.cview(),
                                         s.muv.cview(), c, real(b.rdx), real(b.rdy), do_y, b.msftx.view(),
                                         b.msfty.view(), b.msfuy.view(), b.msfvx_inv.view(), dvdxi.view()});
    column_for(g.columns(), KMuTCol{dvdxi.cview(), ww1.cview(), s.ww.view(), s.ww_m.view(), s.mu_2.view(), mudf.view(),
                                    s.muts.view(), muave.view(), s.mut.cview(), tend.mu_tend.cview(), b.msfty.view(), c,
                                    dts, p_.epssm, nz, it, n_small, ed});
    parallel_for(g.mass_cells(), KMuTPoint{s.u_2.cview(), s.v_2.cview(), u_save.cview(), v_save.cview(),
                                           t_save.cview(), tend.t_tend.cview(), s.ww.cview(), s.t_2.view(),
                                           t_2save.view(), s.ru_m.view(), s.rv_m.view(), s.muu.cview(), s.muv.cview(),
                                           c, dts, real(b.rdx), real(b.rdy), nz, it, n_small, do_y, b.msftx.view(),
                                           b.msfty.view(), b.msfuy.view(), b.msfvx_inv.view(), ed});
    if (bc_) bc_->update_mu_t(s, tend, dts);  // spec_bdyupdate t, mu, muts
    fill_halos({&s.mu_2, &mudf, &s.muts}, g);

    {
      const KWCommon q{s.w_2.cview(), s.ph_2.cview(), tend.rw_tend.cview(), tend.ph_tend.cview(), s.ww.cview(),
                       ph_save.cview(), s.t_2.cview(), t_2save.cview(), t_save.cview(), c2a.cview(), s.alt.cview(),
                       s.mut.cview(), s.muts.cview(), muave.cview(), b.phb.view(), c, mc, dts, p_.epssm, nz,
                       p_.top_lid, b.msfty.view()};
      const Range3D wr{0, nx, 0, ny, 0, nz + 1};
      parallel_for(wr, KWRhs{q, w_rhs.view(), t2a.view(), ed});
      parallel_for(wr, KWExpl{q, w_rhs.cview(), t2a.cview(), s.u_2.cview(), s.v_2.cview(), b.ht.view(),
                              b.msftx.view(), real(b.rdx), real(b.rdy), !b.ht.uniform(), do_y, s.w_2.view(), ed});
      column_for(g.columns(), KWSolve{s.w_2.view(), s.ph_2.view(), w_rhs.cview(), a.cview(), alpha.cview(),
                                      gamma.cview(), ph_save.cview(), w_save.cview(), s.mut.cview(), s.muts.cview(),
                                      b.phb.view(), b.msfty.view(), c, dts, p_.epssm,
                                      p_.damp_opt == 3 ? dts * p_.dampcoef : real(0), p_.zdamp, nz, ed});
      if (bc_) bc_->update_ph_w(s, ph_save, tend, b, dts);  // spec_bdyupdate_ph, zero_grad_bdy(w)
      parallel_for(g.mass_cells(), KWPRho{s.al.view(), s.p.view(), pm1.view(), s.ph_2.cview(), s.alt.cview(),
                                          s.t_2.cview(), t_save.cview(), c2a.cview(), s.mu_2.cview(), s.muts.cview(),
                                          c, p_.smdiv});
    }
    fill_halos({&s.ph_2, &s.p, &s.al}, g);
  }

  fill_halo(s.muts, g);
  fill_halos({&s.ru_m, &s.rv_m, &s.ww_m}, g);
  parallel_for(Range3D{0, nx + 1, 0, nyv, 0, nz + 1},
               KFinish{s.u_2.view(), s.v_2.view(), s.w_2.view(), s.t_2.view(), s.ph_2.view(), s.ww.view(),
                       u_save.cview(), v_save.cview(), w_save.cview(), t_save.cview(), ph_save.cview(), ww1.cview(),
                       s.muu.cview(), s.muv.cview(), s.mut.cview(), s.muts.cview(), c, nx, ny, nyv, nz, do_y,
                       s.h_diabatic.cview(), (stage_ == 3 && s.n_moist > 0) ? dts * real(n_small) : real(0),
                       b.msfuy.view(), b.msfvx.view(), b.msfty.view()});
  parallel_for(Range3D{0, nx + 1, 0, nyv, 0, 1}, KFinishMu{s.mu_2.view(), s.muus.view(), s.muvs.view(),
                                                            mu_save.cview(), s.muts.cview(), nx, ny, nyv});
  fill_halos({&s.u_2, &s.v_2, &s.w_2, &s.t_2, &s.ph_2, &s.ww}, g);
  fill_halo(s.mu_2, g);
}

void SmallStep::calc_p_rho_phi(State& s, const BaseState& b) {
  const Grid& g = g_;
  const VCoef c = VCoef::from(b);
  parallel_for(g.mass_cells(), KPRhoPhi{s.al.view(), s.p.view(), s.ph_2.cview(), s.t_2.cview(), s.mu_2.cview(),
                                        s.muts.cview(), b.alb.view(), b.pb.view(), c, p_.hypsometric_opt,
                                        b.phb.view(), b.c3h.view(), b.c4h.view(), b.c3f.view(), b.c4f.view(),
                                        real(b.p_top)});
  fill_halo(s.p, g);
  fill_halo(s.al, g);
}

}  // namespace cppwrf::dyn
