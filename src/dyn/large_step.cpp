#include "dyn/large_step.hpp"

#include "core/constants.hpp"
#include "core/halo.hpp"
#include "dyn/kernels_common.hpp"

// Formulas follow TN ch. 3 (flux-form, perturbation equations 3.3-3.8 evaluated at t*) in
// the discrete form of WRF dyn_em (module_big_step_utilities_em.F, module_advect_em.F).
// C++ k = Fortran k - 1. Map factors = 1 (idealized). Moisture: dry (cqu = cqv = cqw = 1).

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

using kern::flux3;
using kern::flux5;

// ---------------------------------------------------------------- rk_step_prep
struct KMu {  // calc_mu_uv: full column masses at stage
  C2 mu; B2 mub; V2 mut, muu, muv;
  int nx, ny, nyv;
  CPPWRF_HD void operator()(int i, int j, int) const {
    if (i < nx && j < ny) mut(i, j) = mu(i, j) + mub(i, j);
    if (j < ny) muu(i, j) = H * (mu(i, j) + mub(i, j) + mu(i - 1, j) + mub(i - 1, j));
    if (i < nx && j < nyv) muv(i, j) = H * (mu(i, j) + mub(i, j) + mu(i, j - 1) + mub(i, j - 1));
  }
};

struct KPrep3 {  // couple_momentum, calc_alt, calc_php
  C3 u, v, w, al, ph; C2 muu, muv, mut; B3 alb, phb; VCoef c;
  V3 ru, rv, rw, alt, php;
  int nx, ny, nyv, nz;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    if (j < ny && k < nz) ru(i, j, k) = (c.c1h(k) * muu(i, j) + c.c2h(k)) * u(i, j, k);
    if (i < nx && j < nyv && k < nz) rv(i, j, k) = (c.c1h(k) * muv(i, j) + c.c2h(k)) * v(i, j, k);
    if (i < nx && j < ny) rw(i, j, k) = (c.c1f(k) * mut(i, j) + c.c2f(k)) * w(i, j, k);
    if (i < nx && j < ny && k < nz) {
      alt(i, j, k) = al(i, j, k) + alb(i, j, k);
      php(i, j, k) = H * (phb(i, j, k) + phb(i, j, k + 1) + ph(i, j, k) + ph(i, j, k + 1));
    }
  }
};

struct KWW {  // calc_ww_cp: mass-coupled eta velocity from the column-integrated divergence
  C3 ru, rv; V3 ww; VCoef c; real rdx, rdy; int nz;
  CPPWRF_HD real divv(int i, int j, int k) const {
    return c.dnw(k) * (rdx * (ru(i + 1, j, k) - ru(i, j, k)) + rdy * (rv(i, j + 1, k) - rv(i, j, k)));
  }
  CPPWRF_HD void operator()(int i, int j) const {
    real dmdt = 0;
    for (int k = 0; k < nz; ++k) dmdt += divv(i, j, k);
    real w = 0;
    ww(i, j, 0) = 0;
    for (int k = 1; k < nz; ++k) {
      w = w - c.dnw(k - 1) * c.c1h(k - 1) * dmdt - divv(i, j, k - 1);
      ww(i, j, k) = w;
    }
    ww(i, j, nz) = 0;
  }
};

// ---------------------------------------------------------------- rk_tendency
// U: advect_u (5th horizontal, 3rd vertical) + horizontal_pressure_gradient (+ ru_tendf)
struct KTendU {
  C3 u, ru, rv, ww, p, ph, alt, al, php, tendf; C2 mu, muu; B3 pb; VCoef c; kern::MoistCoef mc;
  V3 out; real rdx, rdy; int nz; bool do_y, top_lid, use_tendf; real hup, vup;
  CPPWRF_HD real fx(int i, int j, int k) const {
    real vel = H * (ru(i, j, k) + ru(i - 1, j, k));
    return vel * flux5(u(i - 3, j, k), u(i - 2, j, k), u(i - 1, j, k), u(i, j, k), u(i + 1, j, k), u(i + 2, j, k), vel, hup);
  }
  CPPWRF_HD real fy(int i, int j, int k) const {
    real vel = H * (rv(i, j, k) + rv(i - 1, j, k));
    return vel * flux5(u(i, j - 3, k), u(i, j - 2, k), u(i, j - 1, k), u(i, j, k), u(i, j + 1, k), u(i, j + 2, k), vel, hup);
  }
  CPPWRF_HD real fz(int i, int j, int k) const {
    if (k == 0 || k == nz) return 0;
    real vel = H * (ww(i - 1, j, k) + ww(i, j, k));
    if (k == 1 || k == nz - 1) return vel * (c.fnm(k) * u(i, j, k) + c.fnp(k) * u(i, j, k - 1));
    return vel * flux3(u(i, j, k - 2), u(i, j, k - 1), u(i, j, k), u(i, j, k + 1), -vel, vup);
  }
  CPPWRF_HD real pu(int i, int j, int k) const { return p(i - 1, j, k) + p(i, j, k); }
  CPPWRF_HD real dpn(int i, int j, int k) const {  // p' at w level k, u column (x2 averaged)
    if (k == 0) return H * (c.cf1 * pu(i, j, 0) + c.cf2 * pu(i, j, 1) + c.cf3 * pu(i, j, 2));
    if (k == nz) return top_lid ? H * (c.cfn * pu(i, j, nz - 1) + c.cfn1 * pu(i, j, nz - 2)) : real(0);
    return H * (c.fnm(k) * pu(i, j, k) + c.fnp(k) * pu(i, j, k - 1));
  }
  CPPWRF_HD void operator()(int i, int j, int k) const {
    real t = -rdx * (fx(i + 1, j, k) - fx(i, j, k));
    if (do_y) t -= rdy * (fy(i, j + 1, k) - fy(i, j, k));
    t -= c.rdnw(k) * (fz(i, j, k + 1) - fz(i, j, k));
    // TN eq. 2.? / WRF horizontal_pressure_gradient (non-hydrostatic)
    real dpx = H * rdx * (c.c1h(k) * muu(i, j) + c.c2h(k)) *
               ((ph(i, j, k + 1) - ph(i - 1, j, k + 1) + ph(i, j, k) - ph(i - 1, j, k)) +
                (alt(i, j, k) + alt(i - 1, j, k)) * (p(i, j, k) - p(i - 1, j, k)) +
                (al(i, j, k) + al(i - 1, j, k)) * (pb(i, j, k) - pb(i - 1, j, k)));
    dpx += rdx * (php(i, j, k) - php(i - 1, j, k)) *
           (c.rdnw(k) * (dpn(i, j, k + 1) - dpn(i, j, k)) - H * (c.c1h(k) * mu(i - 1, j) + c.c1h(k) * mu(i, j)));
    t -= mc.cqu(i, j, k) * dpx;
    if (use_tendf) t += tendf(i, j, k);
    out(i, j, k) = t;
  }
};

// V: advect_v + horizontal_pressure_gradient (y part); 3D runs only
struct KTendV {
  C3 v, ru, rv, ww, p, ph, alt, al, php, tendf; C2 mu, muv; B3 pb; VCoef c; kern::MoistCoef mc;
  V3 out; real rdx, rdy; int nz; bool top_lid, use_tendf; real hup, vup;
  CPPWRF_HD real fx(int i, int j, int k) const {
    real vel = H * (ru(i, j, k) + ru(i, j - 1, k));
    return vel * flux5(v(i - 3, j, k), v(i - 2, j, k), v(i - 1, j, k), v(i, j, k), v(i + 1, j, k), v(i + 2, j, k), vel, hup);
  }
  CPPWRF_HD real fy(int i, int j, int k) const {
    real vel = H * (rv(i, j, k) + rv(i, j - 1, k));
    return vel * flux5(v(i, j - 3, k), v(i, j - 2, k), v(i, j - 1, k), v(i, j, k), v(i, j + 1, k), v(i, j + 2, k), vel, hup);
  }
  CPPWRF_HD real fz(int i, int j, int k) const {
    if (k == 0 || k == nz) return 0;
    real vel = H * (ww(i, j - 1, k) + ww(i, j, k));
    if (k == 1 || k == nz - 1) return vel * (c.fnm(k) * v(i, j, k) + c.fnp(k) * v(i, j, k - 1));
    return vel * flux3(v(i, j, k - 2), v(i, j, k - 1), v(i, j, k), v(i, j, k + 1), -vel, vup);
  }
  CPPWRF_HD real pv(int i, int j, int k) const { return p(i, j - 1, k) + p(i, j, k); }
  CPPWRF_HD real dpn(int i, int j, int k) const {
    if (k == 0) return H * (c.cf1 * pv(i, j, 0) + c.cf2 * pv(i, j, 1) + c.cf3 * pv(i, j, 2));
    if (k == nz) return top_lid ? H * (c.cfn * pv(i, j, nz - 1) + c.cfn1 * pv(i, j, nz - 2)) : real(0);
    return H * (c.fnm(k) * pv(i, j, k) + c.fnp(k) * pv(i, j, k - 1));
  }
  CPPWRF_HD void operator()(int i, int j, int k) const {
    real t = -rdx * (fx(i + 1, j, k) - fx(i, j, k)) - rdy * (fy(i, j + 1, k) - fy(i, j, k));
    t -= c.rdnw(k) * (fz(i, j, k + 1) - fz(i, j, k));
    real dpy = H * rdy * (c.c1h(k) * muv(i, j) + c.c2h(k)) *
               ((ph(i, j, k + 1) - ph(i, j - 1, k + 1) + ph(i, j, k) - ph(i, j - 1, k)) +
                (alt(i, j, k) + alt(i, j - 1, k)) * (p(i, j, k) - p(i, j - 1, k)) +
                (al(i, j, k) + al(i, j - 1, k)) * (pb(i, j, k) - pb(i, j - 1, k)));
    dpy += rdy * (php(i, j, k) - php(i, j - 1, k)) *
           (c.rdnw(k) * (dpn(i, j, k + 1) - dpn(i, j, k)) - H * (c.c1h(k) * mu(i, j - 1) + c.c1h(k) * mu(i, j)));
    t -= mc.cqv(i, j, k) * dpy;
    if (use_tendf) t += tendf(i, j, k);
    out(i, j, k) = t;
  }
};

// W (levels 1..nz): advect_w + pg_buoy_w (+ rw_tendf)
struct KTendW {
  C3 w, ru, rv, ww, p, tendf; C2 mu; B2 mub; VCoef c; kern::MoistCoef mc;
  V3 out; real rdx, rdy; int nz; bool do_y, use_tendf; real hup, vup;
  CPPWRF_HD real rux(int i, int j, int k) const {  // ru interpolated to w level k
    if (k < nz) return c.fnm(k) * ru(i, j, k) + c.fnp(k) * ru(i, j, k - 1);
    return (real(2) - c.fnm(nz - 1)) * ru(i, j, nz - 1) - c.fnp(nz - 1) * ru(i, j, nz - 2);
  }
  CPPWRF_HD real rvy(int i, int j, int k) const {
    if (k < nz) return c.fnm(k) * rv(i, j, k) + c.fnp(k) * rv(i, j, k - 1);
    return (real(2) - c.fnm(nz - 1)) * rv(i, j, nz - 1) - c.fnp(nz - 1) * rv(i, j, nz - 2);
  }
  CPPWRF_HD real fx(int i, int j, int k) const {
    real vel = rux(i, j, k);
    return vel * flux5(w(i - 3, j, k), w(i - 2, j, k), w(i - 1, j, k), w(i, j, k), w(i + 1, j, k), w(i + 2, j, k), vel, hup);
  }
  CPPWRF_HD real fy(int i, int j, int k) const {
    real vel = rvy(i, j, k);
    return vel * flux5(w(i, j - 3, k), w(i, j - 2, k), w(i, j - 1, k), w(i, j, k), w(i, j + 1, k), w(i, j + 2, k), vel, hup);
  }
  CPPWRF_HD real fz(int i, int j, int m) const {  // vertical flux at mass level m
    if (m == 0 || m == nz - 1) return real(0.25) * (ww(i, j, m + 1) + ww(i, j, m)) * (w(i, j, m + 1) + w(i, j, m));
    real vel = H * (ww(i, j, m + 1) + ww(i, j, m));
    return vel * flux3(w(i, j, m - 1), w(i, j, m), w(i, j, m + 1), w(i, j, m + 2), -vel, vup);
  }
  CPPWRF_HD void operator()(int i, int j, int k) const {
    real t = -rdx * (fx(i + 1, j, k) - fx(i, j, k));
    if (do_y) t -= rdy * (fy(i, j + 1, k) - fy(i, j, k));
    if (k < nz) {
      t -= c.rdn(k) * (fz(i, j, k) - fz(i, j, k - 1));
      // pg_buoy_w: cq1 = 1/(1+cqw), cq2 = cqw cq1 (moisture loading of the base state)
      const real cqw = mc.cqw(i, j, k), cq1 = real(1) / (real(1) + cqw), cq2 = cqw * cq1;
      t += Gr * (cq1 * c.rdn(k) * (p(i, j, k) - p(i, j, k - 1)) - c.c1f(k) * mu(i, j) -
                 cq2 * (c.c1f(k) * mub(i, j) + c.c2f(k)));
    } else {
      t += real(2) * c.rdn(nz - 1) * fz(i, j, nz - 1);
      const real cqw = mc.cqw(i, j, nz - 1), cq1 = real(1) / (real(1) + cqw), cq2 = cqw * cq1;
      t += Gr * (cq1 * real(2) * c.rdnw(nz - 1) * (-p(i, j, nz - 1)) - c.c1f(nz) * mu(i, j) -
                 cq2 * (c.c1f(nz) * mub(i, j) + c.c2f(nz)));
    }
    if (use_tendf) t += tendf(i, j, k);
    out(i, j, k) = t;
  }
};

// Theta: advect_scalar (+ t_tendf)
struct KTendT {
  C3 t, ru, rv, ww, tendf; VCoef c;
  V3 out; real rdx, rdy; int nz; bool do_y, use_tendf; real hup, vup;
  C3 hd; C2 mut; bool use_hd;  // rk_addtend_dry: + (c1h*mut + c2h) * h_diabatic
  CPPWRF_HD real fx(int i, int j, int k) const {
    real vel = ru(i, j, k);
    return vel * flux5(t(i - 3, j, k), t(i - 2, j, k), t(i - 1, j, k), t(i, j, k), t(i + 1, j, k), t(i + 2, j, k), vel, hup);
  }
  CPPWRF_HD real fy(int i, int j, int k) const {
    real vel = rv(i, j, k);
    return vel * flux5(t(i, j - 3, k), t(i, j - 2, k), t(i, j - 1, k), t(i, j, k), t(i, j + 1, k), t(i, j + 2, k), vel, hup);
  }
  CPPWRF_HD real fz(int i, int j, int k) const {
    if (k == 0 || k == nz) return 0;
    real vel = ww(i, j, k);
    if (k == 1 || k == nz - 1) return vel * (c.fnm(k) * t(i, j, k) + c.fnp(k) * t(i, j, k - 1));
    return vel * flux3(t(i, j, k - 2), t(i, j, k - 1), t(i, j, k), t(i, j, k + 1), -vel, vup);
  }
  CPPWRF_HD void operator()(int i, int j, int k) const {
    real r = -rdx * (fx(i + 1, j, k) - fx(i, j, k));
    if (do_y) r -= rdy * (fy(i, j + 1, k) - fy(i, j, k));
    r -= c.rdnw(k) * (fz(i, j, k + 1) - fz(i, j, k));
    if (use_tendf) r += tendf(i, j, k);
    if (use_hd) r += (c.c1h(k) * mut(i, j) + c.c2h(k)) * hd(i, j, k);
    out(i, j, k) = r;
  }
};

// Geopotential (levels 1..nz): rhs_ph with phi_adv_z = 1, 6th-order centred horizontal
// gradient (h_sca_adv_order = 5 branch), plus the g*w term (non-hydrostatic).
struct KTendPh {
  C3 u, v, w, ww, ph; C2 mut, muu, muv; B3 phb; VCoef c;
  V3 out; real rdx, rdy; int nz; bool do_y;
  CPPWRF_HD real phf(int i, int j, int k) const { return ph(i, j, k) + phb(i, j, k); }
  CPPWRF_HD real wdwn(int i, int j, int m) const {
    return H * (ww(i, j, m + 1) + ww(i, j, m)) * c.rdnw(m) * (phf(i, j, m + 1) - phf(i, j, m));
  }
  CPPWRF_HD void operator()(int i, int j, int k) const {
    real t = 0;
    if (k < nz) t -= c.fnm(k) * wdwn(i, j, k) + c.fnp(k) * wdwn(i, j, k - 1);
    t += (c.c1f(k) * mut(i, j) + c.c2f(k)) * Gr * w(i, j, k);
    real dx6 = (real(45) * (phf(i + 1, j, k) - phf(i - 1, j, k)) - real(9) * (phf(i + 2, j, k) - phf(i - 2, j, k)) +
                (phf(i + 3, j, k) - phf(i - 3, j, k))) / real(60);
    real mup = c.c1f(k) * muu(i + 1, j) + c.c2f(k), mum = c.c1f(k) * muu(i, j) + c.c2f(k);
    if (k < nz)
      t -= real(0.25) * rdx * (mup * (u(i + 1, j, k) + u(i + 1, j, k - 1)) + mum * (u(i, j, k) + u(i, j, k - 1))) * dx6;
    else
      t -= H * rdx *
           (mup * (c.cfn * u(i + 1, j, nz - 1) + c.cfn1 * u(i + 1, j, nz - 2)) +
            mum * (c.cfn * u(i, j, nz - 1) + c.cfn1 * u(i, j, nz - 2))) * dx6;
    if (do_y) {
      real dy6 = (real(45) * (phf(i, j + 1, k) - phf(i, j - 1, k)) - real(9) * (phf(i, j + 2, k) - phf(i, j - 2, k)) +
                  (phf(i, j + 3, k) - phf(i, j - 3, k))) / real(60);
      real mvp = c.c1f(k) * muv(i, j + 1) + c.c2f(k), mvm = c.c1f(k) * muv(i, j) + c.c2f(k);
      if (k < nz)
        t -= real(0.25) * rdy * (mvp * (v(i, j + 1, k) + v(i, j + 1, k - 1)) + mvm * (v(i, j, k) + v(i, j, k - 1))) * dy6;
      else
        t -= H * rdy *
             (mvp * (c.cfn * v(i, j + 1, nz - 1) + c.cfn1 * v(i, j + 1, nz - 2)) +
              mvm * (c.cfn * v(i, j, nz - 1) + c.cfn1 * v(i, j, nz - 2))) * dy6;
    }
    out(i, j, k) = t;
  }
};

// ---------------------------------------------------------------- constant-K diffusion
// Horizontal: along eta surfaces, mass-weighted (WRF horizontal_diffusion, diff_opt=1).
// Vertical: physical-space flux via dz = -alpha mu deta / g (WRF vertical_diffusion*), with
// zero-flux (free-slip / insulating) bottom and top for u, v, theta and w = 0 at both ends.
struct KDiffU {
  C3 u, alt; C2 mut, muu; VCoef c; V3 out; real kh, kv, rdx, rdy; int nz; bool do_y;
  CPPWRF_HD real m(int i, int j, int k) const { return c.c1h(k) * mut(i, j) + c.c2h(k); }
  CPPWRF_HD real vf(int i, int j, int k) const {  // vertical flux at w level k
    if (k == nz) return 0;
    if (k == 0) k = 1;  // WRF vertical_diffusion_u/v: vflux(0) = vflux(1)
    real a = real(0.25) * (alt(i, j, k) + alt(i - 1, j, k) + alt(i, j, k - 1) + alt(i - 1, j, k - 1));
    return kv * c.rdn(k) / a * (u(i, j, k) - u(i, j, k - 1));
  }
  CPPWRF_HD void operator()(int i, int j, int k) const {
    real t = rdx * rdx * kh * (m(i, j, k) * (u(i + 1, j, k) - u(i, j, k)) - m(i - 1, j, k) * (u(i, j, k) - u(i - 1, j, k)));
    if (do_y) {
      real mp = real(0.25) * (m(i, j, k) + m(i, j + 1, k) + m(i - 1, j + 1, k) + m(i - 1, j, k));
      real mm = real(0.25) * (m(i, j, k) + m(i, j - 1, k) + m(i - 1, j - 1, k) + m(i - 1, j, k));
      t += rdy * rdy * kh * (mp * (u(i, j + 1, k) - u(i, j, k)) - mm * (u(i, j, k) - u(i, j - 1, k)));
    }
    real mu_u = c.c1h(k) * muu(i, j) + c.c2h(k);
    if (k < nz - 1)  // WRF: the top level gets no vertical diffusion (tendency loop to ktf-1)
      t += Gr * Gr * c.rdnw(k) / (mu_u * H * (alt(i - 1, j, k) + alt(i, j, k))) * (vf(i, j, k + 1) - vf(i, j, k));
    out(i, j, k) = t;
  }
};

struct KDiffV {
  C3 v, alt; C2 mut, muv; VCoef c; V3 out; real kh, kv, rdx, rdy; int nz;
  CPPWRF_HD real m(int i, int j, int k) const { return c.c1h(k) * mut(i, j) + c.c2h(k); }
  CPPWRF_HD real vf(int i, int j, int k) const {
    if (k == nz) return 0;
    if (k == 0) k = 1;  // WRF vertical_diffusion_u/v: vflux(0) = vflux(1)
    real a = real(0.25) * (alt(i, j, k) + alt(i, j - 1, k) + alt(i, j, k - 1) + alt(i, j - 1, k - 1));
    return kv * c.rdn(k) / a * (v(i, j, k) - v(i, j, k - 1));
  }
  CPPWRF_HD void operator()(int i, int j, int k) const {
    real t = rdy * rdy * kh * (m(i, j, k) * (v(i, j + 1, k) - v(i, j, k)) - m(i, j - 1, k) * (v(i, j, k) - v(i, j - 1, k)));
    real mp = real(0.25) * (m(i, j, k) + m(i + 1, j, k) + m(i + 1, j - 1, k) + m(i, j - 1, k));
    real mm = real(0.25) * (m(i, j, k) + m(i - 1, j, k) + m(i - 1, j - 1, k) + m(i, j - 1, k));
    t += rdx * rdx * kh * (mp * (v(i + 1, j, k) - v(i, j, k)) - mm * (v(i, j, k) - v(i - 1, j, k)));
    real mu_v = c.c1h(k) * muv(i, j) + c.c2h(k);
    if (k < nz - 1)  // WRF: the top level gets no vertical diffusion (tendency loop to ktf-1)
      t += Gr * Gr * c.rdnw(k) / (mu_v * H * (alt(i, j - 1, k) + alt(i, j, k))) * (vf(i, j, k + 1) - vf(i, j, k));
    out(i, j, k) = t;
  }
};

struct KDiffW {  // levels 1..nz-1 (w = 0 at the rigid bottom and lid)
  C3 w, alt; C2 mut; VCoef c; V3 out; real kh, kv, rdx, rdy; int nz; bool do_y;
  CPPWRF_HD real m(int i, int j, int k) const { return c.c1f(k) * mut(i, j) + c.c2f(k); }
  CPPWRF_HD real vf(int i, int j, int mm) const {  // flux at mass level mm
    if (mm == nz - 1) return 0;  // WRF vertical_diffusion('w'): vflux(ktf) = 0
    return kv / alt(i, j, mm) * c.rdnw(mm) * (w(i, j, mm + 1) - w(i, j, mm));
  }
  CPPWRF_HD void operator()(int i, int j, int k) const {
    real t = rdx * rdx * kh *
             (H * (m(i + 1, j, k) + m(i, j, k)) * (w(i + 1, j, k) - w(i, j, k)) -
              H * (m(i, j, k) + m(i - 1, j, k)) * (w(i, j, k) - w(i - 1, j, k)));
    if (do_y)
      t += rdy * rdy * kh *
           (H * (m(i, j + 1, k) + m(i, j, k)) * (w(i, j + 1, k) - w(i, j, k)) -
            H * (m(i, j, k) + m(i, j - 1, k)) * (w(i, j, k) - w(i, j - 1, k)));
    t += c.rdn(k) * Gr * Gr / m(i, j, k) / (H * (alt(i, j, k) + alt(i, j, k - 1))) * (vf(i, j, k) - vf(i, j, k - 1));
    out(i, j, k) = t;
  }
};


}  // namespace

VCoef VCoef::from(const BaseState& b) {
  return VCoef{b.dnw.view(), b.rdnw.view(), b.dn.view(), b.rdn.view(), b.fnm.view(), b.fnp.view(),
               b.c1h.view(), b.c2h.view(), b.c1f.view(), b.c2f.view(),
               real(b.cf1), real(b.cf2), real(b.cf3), real(b.cfn), real(b.cfn1)};
}

void LargeStep::prep(State& s, const BaseState& b) const {
  const Grid& g = g_;
  const int nyv = g.nys(Stagger::V);
  for (Field3* f : {&s.u_2, &s.v_2, &s.w_2, &s.t_2, &s.ph_2, &s.p, &s.al}) fill_halo(*f, g);
  fill_halo(s.mu_2, g);
  VCoef c = VCoef::from(b);
  parallel_for(Range3D{0, g.nx + 1, 0, nyv, 0, 1},
               KMu{s.mu_2.cview(), b.mub.view(), s.mut.view(), s.muu.view(), s.muv.view(), g.nx, g.ny, nyv});
  for (Field2* f : {&s.mut, &s.muu, &s.muv}) fill_halo(*f, g);
  parallel_for(Range3D{0, g.nx + 1, 0, nyv, 0, g.nz + 1},
               KPrep3{s.u_2.cview(), s.v_2.cview(), s.w_2.cview(), s.al.cview(), s.ph_2.cview(), s.muu.cview(),
                      s.muv.cview(), s.mut.cview(), b.alb.view(), b.phb.view(), c, s.ru.view(), s.rv.view(),
                      s.rw.view(), s.alt.view(), s.php.view(), g.nx, g.ny, nyv, g.nz});
  for (Field3* f : {&s.ru, &s.rv, &s.rw, &s.alt, &s.php}) fill_halo(*f, g);
  column_for(g.columns(), KWW{s.ru.cview(), s.rv.cview(), s.ww.view(), c, real(b.rdx), real(b.rdy), g.nz});
  fill_halo(s.ww, g);
  if (s.n_moist > 0) {  // calc_cq input: qtot = sum of moist species (stage state)
    V3 qt = s.qtot.view();
    parallel_for(g.mass_cells(), [=] CPPWRF_HD(int i, int j, int k) { qt(i, j, k) = 0; });
    for (int n = 0; n < s.n_moist; ++n) {
      C3 q = s.q_2[n].cview();
      parallel_for(g.mass_cells(), [=] CPPWRF_HD(int i, int j, int k) { qt(i, j, k) += q(i, j, k); });
    }
    fill_halo(s.qtot, g);
  }
}

kern::MoistCoef moist_coef(const State& s) {
  return s.n_moist > 0 ? kern::MoistCoef{s.qtot.cview(), true} : kern::MoistCoef{};
}

void LargeStep::tendencies(const State& s, Tendencies& t, const BaseState& b, int rk_step) const {
  const Grid& g = g_;
  const VCoef c = VCoef::from(b);
  const real rdx = real(b.rdx), rdy = real(b.rdy);
  const bool do_y = !g.two_d();
  const bool diff = p_.khdif > 0 || p_.kvdif > 0;
  const int nz = g.nz;
  const Range3D ur{g.u_i0(), g.u_i1(), 0, g.ny, 0, nz};
  const Range3D vr{0, g.nx, g.v_j0(), g.v_j1(), 0, nz};
  const Range3D wr{0, g.nx, 0, g.ny, 1, nz + 1};
  const Range3D wri{0, g.nx, 0, g.ny, 1, nz};
  const Range3D mr = g.mass_cells();
  const real kh = p_.khdif, kv = p_.kvdif;
  // even orders (6 / 4) = the same fluxes without the upwind dissipation term (WRF flux6/flux4)
  const real hup = p_.h_adv_order == 5 ? real(1) : real(0), vup = p_.v_adv_order == 3 ? real(1) : real(0);

  if (diff && rk_step == 1) {  // WRF: diffusion only in the first RK stage, held in *_tendf
    parallel_for(ur, KDiffU{s.u_2.cview(), s.alt.cview(), s.mut.cview(), s.muu.cview(), c, t.ru_tendf.view(), kh, kv,
                            rdx, rdy, nz, do_y});
    if (do_y)
      parallel_for(vr, KDiffV{s.v_2.cview(), s.alt.cview(), s.mut.cview(), s.muv.cview(), c, t.rv_tendf.view(), kh,
                              kv, rdx, rdy, nz});
    parallel_for(wri, KDiffW{s.w_2.cview(), s.alt.cview(), s.mut.cview(), c, t.rw_tendf.view(), kh, kv, rdx, rdy, nz,
                             do_y});
    // theta: horizontal/vertical_diffusion_3dmp on theta - t_init with khdq = kvdq = K / prandtl
    parallel_for(mr, kern::KDiffM{s.t_2.cview(), s.alt.cview(), b.t_init.view(), s.mut.cview(), c, t.t_tendf.view(),
                                  kh * p_.pr_inv, kv * p_.pr_inv, rdx, rdy, nz, do_y, true, true});
  }
  parallel_for(ur, KTendU{s.u_2.cview(), s.ru.cview(), s.rv.cview(), s.ww.cview(), s.p.cview(), s.ph_2.cview(),
                          s.alt.cview(), s.al.cview(), s.php.cview(), t.ru_tendf.cview(), s.mu_2.cview(),
                          s.muu.cview(), b.pb.view(), c, moist_coef(s), t.ru_tend.view(), rdx, rdy, nz, do_y, p_.top_lid, diff, hup, vup});
  if (do_y)
    parallel_for(vr, KTendV{s.v_2.cview(), s.ru.cview(), s.rv.cview(), s.ww.cview(), s.p.cview(), s.ph_2.cview(),
                            s.alt.cview(), s.al.cview(), s.php.cview(), t.rv_tendf.cview(), s.mu_2.cview(),
                            s.muv.cview(), b.pb.view(), c, moist_coef(s), t.rv_tend.view(), rdx, rdy, nz, p_.top_lid, diff, hup, vup});
  parallel_for(wr, KTendW{s.w_2.cview(), s.ru.cview(), s.rv.cview(), s.ww.cview(), s.p.cview(), t.rw_tendf.cview(),
                          s.mu_2.cview(), b.mub.view(), c, moist_coef(s), t.rw_tend.view(), rdx, rdy, nz, do_y, diff, hup, vup});
  parallel_for(mr, KTendT{s.t_2.cview(), s.ru.cview(), s.rv.cview(), s.ww.cview(), t.t_tendf.cview(), c,
                          t.t_tend.view(), rdx, rdy, nz, do_y, diff, hup, vup, s.h_diabatic.cview(), s.mut.cview(),
                          s.n_moist > 0});
  parallel_for(wr, KTendPh{s.u_2.cview(), s.v_2.cview(), s.w_2.cview(), s.ww.cview(), s.ph_2.cview(), s.mut.cview(),
                           s.muu.cview(), s.muv.cview(), b.phb.view(), c, t.ph_tend.view(), rdx, rdy, nz, do_y});
  // mu_tend stays 0 (dry, no physics); allocated zero.
}

}  // namespace cppwrf::dyn
