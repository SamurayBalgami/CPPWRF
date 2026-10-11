#include "dyn/small_step.hpp"

#include <cmath>

#include "core/backend.hpp"
#include "core/constants.hpp"
#include "core/halo.hpp"
#include "dyn/large_step.hpp"

// Discrete forms follow WRF dyn_em/module_small_step_em.F (small_step_prep, calc_p_rho,
// calc_coef_w, advance_uv, advance_mu_t, advance_w, sumflux, small_step_finish) and
// calc_p_rho_phi (module_big_step_utilities_em.F, hypsometric_opt = 1, which matches the
// linear hydrostatic integration used by the base state and the case initialisation).
// C++ k = Fortran k - 1. Map factors = 1, dry (cqu = cqv = cqw = 1), flat terrain (ht = 0).
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
      u(i, j, k) = (c.c1h(k) * mus_u + c.c2h(k)) * u1(i, j, k) - (c.c1h(k) * muu(i, j) + c.c2h(k)) * ui;
    }
    if (i < nx && j < nyv && k < nz) {
      const real vi = v(i, j, k);
      vs(i, j, k) = vi;
      v(i, j, k) = (c.c1h(k) * mus_v + c.c2h(k)) * v1(i, j, k) - (c.c1h(k) * muv(i, j) + c.c2h(k)) * vi;
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
      w(i, j, k) = (c.c1f(k) * mts + c.c2f(k)) * w1(i, j, k) - (c.c1f(k) * mut(i, j) + c.c2f(k)) * wi;
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
  V3 a, alpha, gamma; C3 c2a; C2 mut; VCoef c; real cof; real lid; int nz;
  CPPWRF_HD real mh(int k, real m) const { return c.c1h(k) * m + c.c2h(k); }
  CPPWRF_HD real mf(int k, real m) const { return c.c1f(k) * m + c.c2f(k); }
  CPPWRF_HD void operator()(int i, int j) const {
    const real m = mut(i, j);
    a(i, j, 0) = 0; alpha(i, j, 0) = 0; gamma(i, j, 0) = 0;
    a(i, j, 1) = 0;
    a(i, j, nz) = real(-2) * cof * c.rdnw(nz - 1) * c.rdnw(nz - 1) * c2a(i, j, nz - 1) * lid /
                  (mh(nz - 1, m) * mf(nz - 1, m));
    for (int k = 2; k < nz; ++k)  // WRF quirk kept: c1f at the level below
      a(i, j, k) = -cof * c.rdn(k) * c.rdnw(k - 1) * c2a(i, j, k - 1) / (mh(k - 1, m) * mf(k - 1, m));
    for (int k = 1; k < nz; ++k) {
      const real b = real(1) + cof * c.rdn(k) *
                                   (c.rdnw(k) * c2a(i, j, k) / (mh(k, m) * mf(k, m)) +
                                    c.rdnw(k - 1) * c2a(i, j, k - 1) / (mh(k - 1, m) * mf(k, m)));
      const real cc = -cof * c.rdn(k) * c.rdnw(k) * c2a(i, j, k) / (mh(k, m) * mf(k + 1, m));
      alpha(i, j, k) = real(1) / (b - a(i, j, k) * gamma(i, j, k - 1));
      gamma(i, j, k) = cc * alpha(i, j, k);
    }
    const real b = real(1) + real(2) * cof * c.rdnw(nz - 1) * c.rdnw(nz - 1) * c2a(i, j, nz - 1) /
                                 (mh(nz - 1, m) * mf(nz, m));
    alpha(i, j, nz) = real(1) / (b - a(i, j, nz) * gamma(i, j, nz - 1));
    gamma(i, j, nz) = 0;
  }
};

// ------------------------------------------------------------------ advance_uv
struct KUV {
  V3 u, v; C3 ru_tend, rv_tend, p, al, alt, ph, php; C2 mu, muu, muv, mudf; B3 pb; VCoef c;
  real dts, rdx, rdy, emdx, emdy; int nz, ui0, ui1, vj0, vj1; bool do_y, top_lid;
  template <int DI, int DJ> CPPWRF_HD real ps(int i, int j, int k) const { return p(i, j, k) + p(i - DI, j - DJ, k); }
  template <int DI, int DJ> CPPWRF_HD real dpn(int i, int j, int k) const {
    if (k == 0) return H * (c.cf1 * ps<DI, DJ>(i, j, 0) + c.cf2 * ps<DI, DJ>(i, j, 1) + c.cf3 * ps<DI, DJ>(i, j, 2));
    if (k == nz)  // WRF advance_uv uses cf1..cf3 at the lid too
      return top_lid ? H * (c.cf1 * ps<DI, DJ>(i, j, nz - 1) + c.cf2 * ps<DI, DJ>(i, j, nz - 2) +
                            c.cf3 * ps<DI, DJ>(i, j, nz - 3))
                     : real(0);
    return H * (c.fnm(k) * ps<DI, DJ>(i, j, k) + c.fnp(k) * ps<DI, DJ>(i, j, k - 1));
  }
  template <int DI, int DJ> CPPWRF_HD real dp(int i, int j, int k, real mface, real rd) const {
    real d = H * rd * (c.c1h(k) * mface + c.c2h(k)) *
             ((ph(i, j, k + 1) - ph(i - DI, j - DJ, k + 1) + ph(i, j, k) - ph(i - DI, j - DJ, k)) +
              (alt(i, j, k) + alt(i - DI, j - DJ, k)) * (p(i, j, k) - p(i - DI, j - DJ, k)) +
              (al(i, j, k) + al(i - DI, j - DJ, k)) * (pb(i, j, k) - pb(i - DI, j - DJ, k)));
    d += rd * (php(i, j, k) - php(i - DI, j - DJ, k)) *
         (c.rdnw(k) * (dpn<DI, DJ>(i, j, k + 1) - dpn<DI, DJ>(i, j, k)) -
          H * (c.c1h(k) * mu(i - DI, j - DJ) + c.c1h(k) * mu(i, j)));
    return d;
  }
  CPPWRF_HD void operator()(int i, int j, int k) const {
    if (i >= ui0 && i < ui1) {
      const real mudf_x = -emdx * (mudf(i, j) - mudf(i - 1, j));
      u(i, j, k) = u(i, j, k) + dts * ru_tend(i, j, k) - dts * dp<1, 0>(i, j, k, muu(i, j), rdx) + c.c1h(k) * mudf_x;
    }
    if (do_y && j >= vj0 && j < vj1) {
      const real mudf_y = -emdy * (mudf(i, j) - mudf(i, j - 1));
      v(i, j, k) = v(i, j, k) + dts * rv_tend(i, j, k) - dts * dp<0, 1>(i, j, k, muv(i, j), rdy) + c.c1h(k) * mudf_y;
    }
  }
};

// ------------------------------------------------------------------ advance_mu_t + sumflux
struct KMuT {
  C3 u, v, us, vs, ts, t_tend, ww1; V3 t, t_ave, ww, ru_m, rv_m, ww_m;
  V2 mu, mudf, muts, muave; C2 mut, muu, muv, mu_tend; VCoef c;
  real dts, rdx, rdy, epssm; int nz, iter, nsteps; bool do_y;
  CPPWRF_HD real dvdxi(int i, int j, int k) const {
    const real fx = (u(i + 1, j, k) + (c.c1h(k) * muu(i + 1, j) + c.c2h(k)) * us(i + 1, j, k)) -
                    (u(i, j, k) + (c.c1h(k) * muu(i, j) + c.c2h(k)) * us(i, j, k));
    real d = rdx * fx;
    if (do_y)
      d += rdy * ((v(i, j + 1, k) + (c.c1h(k) * muv(i, j + 1) + c.c2h(k)) * vs(i, j + 1, k)) -
                  (v(i, j, k) + (c.c1h(k) * muv(i, j) + c.c2h(k)) * vs(i, j, k)));
    return d;
  }
  CPPWRF_HD void operator()(int i, int j) const {
    real dmdt = 0;
    for (int k = 0; k < nz; ++k) dmdt += c.dnw(k) * dvdxi(i, j, k);
    const real mtend = mu_tend(i, j);
    const real mu_old = mu(i, j);
    const real mu_new = mu_old + dts * (dmdt + mtend);
    mu(i, j) = mu_new;
    mudf(i, j) = dmdt + mtend;
    muts(i, j) = mut(i, j) + mu_new;
    muave(i, j) = H * ((real(1) + epssm) * mu_new + (real(1) - epssm) * mu_old);

    const bool first = iter == 1, last = iter == nsteps;
    const real rn = real(1) / real(nsteps);
    // sumflux for the west/south faces of this column (face nx / ny comes from the halo fill)
    for (int k = 0; k < nz; ++k) {
      const real su = (first ? real(0) : ru_m(i, j, k)) + u(i, j, k);
      ru_m(i, j, k) = last ? su * rn + (c.c1h(k) * muu(i, j) + c.c2h(k)) * us(i, j, k) : su;
      if (do_y) {
        const real sv = (first ? real(0) : rv_m(i, j, k)) + v(i, j, k);
        rv_m(i, j, k) = last ? sv * rn + (c.c1h(k) * muv(i, j) + c.c2h(k)) * vs(i, j, k) : sv;
      }
    }
    // omega'' and theta'' marching upward; ww_full is the running integral (full flux)
    real ww_full = 0, wpp_lo = 0;  // ww'' at the bottom is 0
    real wdtn_lo = 0;
    auto acc_ww = [&](int k, real wpp) {
      const real sw = (first ? real(0) : ww_m(i, j, k)) + wpp;
      ww_m(i, j, k) = last ? sw * rn + ww1(i, j, k) : sw;
    };
    ww(i, j, 0) = 0;
    acc_ww(0, real(0));
    for (int k = 0; k < nz; ++k) {
      real wpp_hi = 0, wdtn_hi = 0;
      if (k + 1 < nz) {
        ww_full = ww_full - c.dnw(k) * (c.c1h(k) * dmdt + dvdxi(i, j, k) + c.c1h(k) * mtend);
        wpp_hi = ww_full - ww1(i, j, k + 1);
        ww(i, j, k + 1) = wpp_hi;
        wdtn_hi = wpp_hi * (c.fnm(k + 1) * ts(i, j, k + 1) + c.fnp(k + 1) * ts(i, j, k));
      } else {
        ww(i, j, nz) = 0;
      }
      acc_ww(k + 1, wpp_hi);
      const real tk = t(i, j, k);
      t_ave(i, j, k) = tk;
      real hflux = H * rdx *
                   (u(i + 1, j, k) * (ts(i + 1, j, k) + ts(i, j, k)) - u(i, j, k) * (ts(i, j, k) + ts(i - 1, j, k)));
      if (do_y)
        hflux += H * rdy *
                 (v(i, j + 1, k) * (ts(i, j + 1, k) + ts(i, j, k)) - v(i, j, k) * (ts(i, j, k) + ts(i, j - 1, k)));
      t(i, j, k) = tk + dts * t_tend(i, j, k) - dts * (hflux + c.rdnw(k) * (wdtn_hi - wdtn_lo));
      wdtn_lo = wdtn_hi;
      wpp_lo = wpp_hi;
    }
    (void)wpp_lo;
  }
};

// ------------------------------------------------------------------ advance_w + calc_p_rho
struct KW {
  V3 w, ph, al, p, pm1;
  C3 rw_tend, ph_tend, ww, phs, t, t_ave, ts, c2a, alt, a, alpha, gamma;
  C2 mu, mut, muts, muave; B3 phb; VCoef c;
  real dts, epssm, smdiv; int nz; bool top_lid;
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
  CPPWRF_HD real rhs(int i, int j, int k, real w_old, real ph_old) const {
    if (k == 0) return 0;
    if (top_lid && k == nz) return 0;
    real r = dts * (ph_tend(i, j, k) + H * Gr * (real(1) - epssm) * w_old);
    if (k < nz) r -= dts * (c.fnm(k) * wdwn(i, j, k + 1) + c.fnp(k) * wdwn(i, j, k));
    return ph_old + r / mf(k, mut(i, j));
  }
  CPPWRF_HD void operator()(int i, int j) const {
    const real m = mut(i, j), mts = muts(i, j), mav = muave(i, j);
    const real ope = real(1) + epssm, ome = real(1) - epssm;
    // phase 1: explicit w'' (stored in w) and rhs (stored in ph); rolling old values
    real ph_km1 = ph(i, j, 0), ph_k = ph(i, j, 1);
    real r_km1 = 0, r_k = rhs(i, j, 1, w(i, j, 1), ph_k);
    w(i, j, 0) = 0;  // flat terrain
    for (int k = 1; k <= nz; ++k) {
      if (k < nz) {
        const real ph_kp1 = ph(i, j, k + 1);
        const real r_kp1 = rhs(i, j, k + 1, w(i, j, k + 1), ph_kp1);
        real wn = w(i, j, k) + dts * rw_tend(i, j, k) +
                  H * dts * Gr * c.rdn(k) *
                      (c2a(i, j, k) * c.rdnw(k) / mh(k, m) * (ope * (r_kp1 - r_k) + ome * (ph_kp1 - ph_k)) -
                       c2a(i, j, k - 1) * c.rdnw(k - 1) / mh(k - 1, m) * (ope * (r_k - r_km1) + ome * (ph_k - ph_km1))) +
                  dts * Gr *
                      (c.rdn(k) * (c2a(i, j, k) * alt(i, j, k) * t2ave(i, j, k) -
                                   c2a(i, j, k - 1) * alt(i, j, k - 1) * t2ave(i, j, k - 1)) -
                       c.c1f(k) * mav);
        w(i, j, k) = wn;
        ph(i, j, k) = r_k;
        ph_km1 = ph_k; ph_k = ph_kp1;
        r_km1 = r_k; r_k = r_kp1;
      } else {
        const int km = nz - 1;
        real wn = w(i, j, nz) + dts * rw_tend(i, j, nz) -
                  H * dts * Gr / mh(km, m) * c.rdnw(km) * c.rdnw(km) * real(2) * c2a(i, j, km) *
                      (ope * (r_k - r_km1) + ome * (ph_k - ph_km1)) -
                  dts * Gr * (real(2) * c.rdnw(km) * c2a(i, j, km) * alt(i, j, km) * t2ave(i, j, km) + c.c1f(nz) * mav);
        if (top_lid) wn = 0;
        w(i, j, nz) = wn;
        ph(i, j, nz) = r_k;
      }
    }
    // phase 2: tridiagonal solve (calc_coef_w coefficients), then phi'' = rhs + implicit w part
    for (int k = 1; k <= nz; ++k) w(i, j, k) = (w(i, j, k) - a(i, j, k) * w(i, j, k - 1)) * alpha(i, j, k);
    for (int k = nz - 1; k >= 1; --k) w(i, j, k) = w(i, j, k) - gamma(i, j, k) * w(i, j, k + 1);
    for (int k = nz; k >= 1; --k) ph(i, j, k) = ph(i, j, k) + H * dts * Gr * ope * w(i, j, k) / mf(k, mts);
    // calc_p_rho with divergence damping, this column
    for (int k = 0; k < nz; ++k)
      p_rho_point(i, j, k, al, p, pm1, C3{ph.p, ph.sj, ph.sk}, alt, t, ts, c2a, mu, muts, c, smdiv, true);
  }
};

// ------------------------------------------------------------------ small_step_finish
struct KFinish {
  V3 u, v, w, t, ph, ww; C3 us, vs, ws, ts, phs, ww1; C2 muu, muv, mut, muts; VCoef c;
  int nx, ny, nyv, nz; bool do_y;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    if (j < ny && k < nz) {
      const real mus_u = H * (muts(i, j) + muts(i - 1, j));
      u(i, j, k) = (u(i, j, k) + us(i, j, k) * (c.c1h(k) * muu(i, j) + c.c2h(k))) / (c.c1h(k) * mus_u + c.c2h(k));
    }
    if (i < nx && j < nyv && k < nz) {  // also in 2D, as WRF (v'' is then constant)
      const real mus_v = H * (muts(i, j) + muts(i, j - 1));
      v(i, j, k) = (v(i, j, k) + vs(i, j, k) * (c.c1h(k) * muv(i, j) + c.c2h(k))) / (c.c1h(k) * mus_v + c.c2h(k));
    }
    if (i < nx && j < ny) {
      const real m = mut(i, j), mts = muts(i, j);
      w(i, j, k) = (w(i, j, k) + ws(i, j, k) * (c.c1f(k) * m + c.c2f(k))) / (c.c1f(k) * mts + c.c2f(k));
      ph(i, j, k) = ph(i, j, k) + phs(i, j, k);
      ww(i, j, k) = ww(i, j, k) + ww1(i, j, k);
      if (k < nz)
        t(i, j, k) = (t(i, j, k) + ts(i, j, k) * (c.c1h(k) * m + c.c2h(k))) / (c.c1h(k) * mts + c.c2h(k));
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
  CPPWRF_HD void operator()(int i, int j, int k) const {
    const real a = -(alb(i, j, k) * (c.c1h(k) * mu(i, j)) + c.rdnw(k) * (ph(i, j, k + 1) - ph(i, j, k))) /
                   (c.c1h(k) * muts(i, j) + c.c2h(k));
    al(i, j, k) = a;
    const double base = double(cst::r_d) * (double(T0) + double(t(i, j, k))) /
                        (double(cst::p1000mb) * (double(a) + double(alb(i, j, k))));
    p(i, j, k) = real(double(cst::p1000mb) * pow(base, double(cst::cpovcv)) - double(pb(i, j, k)));
  }
};

}  // namespace

SmallStep::SmallStep(const Grid& g, const SmallStepParams& p)
    : g_(g), p_(p),
      u_save(g, Stagger::U, "u_save"), v_save(g, Stagger::V, "v_save"), w_save(g, Stagger::W, "w_save"),
      t_save(g, Stagger::M, "t_save"), ph_save(g, Stagger::W, "ph_save"), ww1(g, Stagger::W, "ww1"),
      c2a(g, Stagger::M, "c2a"), pm1(g, Stagger::M, "pm1"), t_2save(g, Stagger::M, "t_2save"),
      a(g, Stagger::W, "a"), alpha(g, Stagger::W, "alpha"), gamma(g, Stagger::W, "gamma"),
      mu_save(g, Stagger::M, "mu_save"), mudf(g, Stagger::M, "mudf"), muave(g, Stagger::M, "muave") {}

void SmallStep::run_stage(State& s, const Tendencies& tend, const BaseState& b, real dt_stage, real dts, int n_small,
                          int rk_step) {
  (void)dt_stage;  // used by WRF only for h_diabatic (physics)
  const Grid& g = g_;
  const int nx = g.nx, ny = g.ny, nz = g.nz, nyv = g.nys(Stagger::V);
  const bool do_y = !g.two_d();
  const VCoef c = VCoef::from(b);
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
                     b.pb.view(), c, nx, ny, nyv, nz});
  for (Field3* f : {&s.u_2, &s.v_2, &s.ph_2, &t_save, &u_save, &v_save}) fill_halo(*f, g);
  fill_halo(s.mu_2, g);
  fill_halo(s.muts, g);

  parallel_for(g.mass_cells(), KPRho{s.al.view(), s.p.view(), pm1.view(), s.ph_2.cview(), s.alt.cview(),
                                     s.t_2.cview(), t_save.cview(), c2a.cview(), s.mu_2.cview(), s.muts.cview(), c});
  fill_halo(s.p, g);
  fill_halo(s.al, g);

  const real cof = real(std::pow(0.5 * double(dts) * cst::g * (1.0 + double(p_.epssm)), 2));
  column_for(g.columns(), KCoefW{a.view(), alpha.view(), gamma.view(), c2a.cview(), s.mut.cview(), c, cof,
                                 p_.top_lid ? real(0) : real(1), nz});

  const Range3D uv_range{0, nx, 0, ny, 0, nz};
  for (int it = 1; it <= n_small; ++it) {
    parallel_for(uv_range, KUV{s.u_2.view(), s.v_2.view(), tend.ru_tend.cview(), tend.rv_tend.cview(), s.p.cview(),
                               s.al.cview(), s.alt.cview(), s.ph_2.cview(), s.php.cview(), s.mu_2.cview(),
                               s.muu.cview(), s.muv.cview(), mudf.cview(), b.pb.view(), c, dts, real(b.rdx),
                               real(b.rdy), p_.emdiv * real(g.dx), do_y ? p_.emdiv * real(g.dy) : real(0), nz,
                               g.u_i0(), g.u_i1(), g.v_j0(), g.v_j1(), do_y, p_.top_lid});
    fill_halo(s.u_2, g);
    if (do_y) fill_halo(s.v_2, g);

    column_for(g.columns(),
               KMuT{s.u_2.cview(), s.v_2.cview(), u_save.cview(), v_save.cview(), t_save.cview(), tend.t_tend.cview(),
                    ww1.cview(), s.t_2.view(), t_2save.view(), s.ww.view(), s.ru_m.view(), s.rv_m.view(),
                    s.ww_m.view(), s.mu_2.view(), mudf.view(), s.muts.view(), muave.view(), s.mut.cview(),
                    s.muu.cview(), s.muv.cview(), tend.mu_tend.cview(), c, dts, real(b.rdx), real(b.rdy), p_.epssm,
                    nz, it, n_small, do_y});
    fill_halo(s.mu_2, g);
    fill_halo(mudf, g);

    column_for(g.columns(),
               KW{s.w_2.view(), s.ph_2.view(), s.al.view(), s.p.view(), pm1.view(), tend.rw_tend.cview(),
                  tend.ph_tend.cview(), s.ww.cview(), ph_save.cview(), s.t_2.cview(), t_2save.cview(), t_save.cview(),
                  c2a.cview(), s.alt.cview(), a.cview(), alpha.cview(), gamma.cview(), s.mu_2.cview(), s.mut.cview(),
                  s.muts.cview(), muave.cview(), b.phb.view(), c, dts, p_.epssm, p_.smdiv, nz, p_.top_lid});
    for (Field3* f : {&s.ph_2, &s.p, &s.al}) fill_halo(*f, g);
  }

  fill_halo(s.muts, g);
  for (Field3* f : {&s.ru_m, &s.rv_m, &s.ww_m}) fill_halo(*f, g);
  parallel_for(Range3D{0, nx + 1, 0, nyv, 0, nz + 1},
               KFinish{s.u_2.view(), s.v_2.view(), s.w_2.view(), s.t_2.view(), s.ph_2.view(), s.ww.view(),
                       u_save.cview(), v_save.cview(), w_save.cview(), t_save.cview(), ph_save.cview(), ww1.cview(),
                       s.muu.cview(), s.muv.cview(), s.mut.cview(), s.muts.cview(), c, nx, ny, nyv, nz, do_y});
  parallel_for(Range3D{0, nx + 1, 0, nyv, 0, 1}, KFinishMu{s.mu_2.view(), s.muus.view(), s.muvs.view(),
                                                            mu_save.cview(), s.muts.cview(), nx, ny, nyv});
  for (Field3* f : {&s.u_2, &s.v_2, &s.w_2, &s.t_2, &s.ph_2, &s.ww}) fill_halo(*f, g);
  fill_halo(s.mu_2, g);
  parallel_for(g.mass_cells(), KPRhoPhi{s.al.view(), s.p.view(), s.ph_2.cview(), s.t_2.cview(), s.mu_2.cview(),
                                        s.muts.cview(), b.alb.view(), b.pb.view(), c});
  fill_halo(s.p, g);
  fill_halo(s.al, g);
}

}  // namespace cppwrf::dyn
