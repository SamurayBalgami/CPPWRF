#pragma once
// Flux operators and the mass-point diffusion kernel shared by the large step (theta) and
// scalar transport (WRF module_advect_em.F statement functions, module_diffusion_em.F).
#include "core/constants.hpp"
#include "core/field.hpp"
#include "dyn/large_step.hpp"

namespace cppwrf::dyn::kern {

CPPWRF_HD inline real sgn1(real a) { return a >= real(0) ? real(1) : real(-1); }  // Fortran sign(1.,a)
// 5th-order upwind flux (WRF flux5 = flux6 - sign*dissipation), TN 3.3 (eq. 3.31/3.32 family)
CPPWRF_HD inline real flux5(real qm3, real qm2, real qm1, real q0, real qp1, real qp2, real ua, real up) {
  return (real(37) * (q0 + qm1) - real(8) * (qp1 + qm2) + (qp2 + qm3)) / real(60) -
         up * sgn1(ua) * ((qp2 - qm3) - real(5) * (qp1 - qm2) + real(10) * (q0 - qm1)) / real(60);
}
// 3rd-order upwind flux (WRF flux3 = flux4 + sign*dissipation)
CPPWRF_HD inline real flux3(real qm2, real qm1, real q0, real qp1, real ua, real up) {
  return (real(7) * (q0 + qm1) - (qp1 + qm2)) / real(12) +
         up * sgn1(ua) * ((qp1 - qm2) - real(3) * (q0 - qm1)) / real(12);
}
// advect_scalar_pd's own statement functions (module_advect_em.F): the coefficients are folded
// constants, (37./60.)*(..) - (2./15.)*(..) + (1./60.)*(..), so the roundoff differs from flux5/flux3
CPPWRF_HD inline real flux5_pd(real qm3, real qm2, real qm1, real q0, real qp1, real qp2, real ua, real up) {
  const real f6 = (real(37.) / real(60.)) * (q0 + qm1) - (real(2.) / real(15.)) * (qp1 + qm2) +
                  (real(1.) / real(60.)) * (qp2 + qm3);
  return f6 - up * sgn1(ua) * (real(1.) / real(60.)) * ((qp2 - qm3) - real(5.) * (qp1 - qm2) + real(10.) * (q0 - qm1));
}
CPPWRF_HD inline real flux3_pd(real qm2, real qm1, real q0, real qp1, real ua, real up) {
  const real f4 = (real(7.) / real(12.)) * (q0 + qm1) - (real(1.) / real(12.)) * (qp1 + qm2);
  return f4 + up * sgn1(ua) * (real(1.) / real(12.)) * ((qp1 - qm2) - real(3.) * (q0 - qm1));
}
// Specified (real-data) and open lateral boundaries: WRF degrades the horizontal advection near the
// edges (module_advect_em.F, degrade_xs/xe/ys/ye): 2nd order next to the boundary, 3rd one further in.
struct Spec {
  bool x = false, y = false;  // degraded in x / y (specified or open)
  int nx = 0, ny = 0;         // mass points
  bool ox = false, oy = false;  // open in x / y (WRF open_xs/xe, open_ys/ye): boundary cells are advected
  // specified only: the outermost row is set from the boundary data (no advection tendency)
  CPPWRF_HD bool sx() const { return x && !ox; }
  CPPWRF_HD bool sy() const { return y && !oy; }
  static Spec of(const Grid& g) {
    return Spec{g.flat_x(), g.flat_y() && !g.two_d(), g.nx, g.ny, g.open_x(), g.open_y()};
  }
};
// Order of a flux between mass points f-1 and f (field not staggered in this direction).
CPPWRF_HD inline int ord_face(int f, int n, bool on) {
  if (!on) return 5;
  if (f <= 1 || f >= n - 1) return 2;
  if (f == 2 || f == n - 2) return 3;
  return 5;
}
// Order of a flux at mass point f-1, between faces f-1 and f (field staggered in this direction).
CPPWRF_HD inline int ord_mass(int f, int n, bool on) {
  if (!on) return 5;
  if (f <= 1 || f >= n) return 2;
  if (f == 2 || f == n - 1) return 3;
  return 5;
}

// Open boundaries, advect_scalar / advect_scalar_pd (module_advect_em.F): the outermost column
// (row) takes, instead of the flux divergence, rdx*(ub*(q_old(inner)-q_old(outer)) + q*(ru(b)-ru(a)))
// with the inflow (ub < 0 west, > 0 east) mean normal mass flux of the cell. Returned without sign.
CPPWRF_HD inline real open_x_scalar(View3<const real> q, View3<const real> qo, View3<const real> ru, real rdx,
                                    int i, int j, int k, int nx) {
  if (i == 0) {
    real ub = real(0.5) * (ru(0, j, k) + ru(1, j, k));
    ub = ub < real(0) ? ub : real(0);
    return rdx * (ub * (qo(1, j, k) - qo(0, j, k)) + q(0, j, k) * (ru(1, j, k) - ru(0, j, k)));
  }
  real ub = real(0.5) * (ru(nx - 1, j, k) + ru(nx, j, k));
  ub = ub > real(0) ? ub : real(0);
  return rdx * (ub * (qo(nx - 1, j, k) - qo(nx - 2, j, k)) + q(nx - 1, j, k) * (ru(nx, j, k) - ru(nx - 1, j, k)));
}
CPPWRF_HD inline real open_y_scalar(View3<const real> q, View3<const real> qo, View3<const real> rv, real rdy,
                                    int i, int j, int k, int ny) {
  if (j == 0) {
    real vb = real(0.5) * (rv(i, 0, k) + rv(i, 1, k));
    vb = vb < real(0) ? vb : real(0);
    return rdy * (vb * (qo(i, 1, k) - qo(i, 0, k)) + q(i, 0, k) * (rv(i, 1, k) - rv(i, 0, k)));
  }
  real vb = real(0.5) * (rv(i, ny - 1, k) + rv(i, ny, k));
  vb = vb > real(0) ? vb : real(0);
  return rdy * (vb * (qo(i, ny - 1, k) - qo(i, ny - 2, k)) + q(i, ny - 1, k) * (rv(i, ny, k) - rv(i, ny - 1, k)));
}

// first-order upwind flux with Courant number cr (advect_scalar_pd flux_upwind)
CPPWRF_HD inline real flux_upwind(real qm1, real q0, real cr) {
  const real a = cr < real(0) ? -cr : cr;
  const real lo = cr + a, hi = cr - a;
  return real(0.5) * (lo < real(1) ? lo : real(1)) * qm1 + real(0.5) * (hi > real(-1) ? hi : real(-1)) * q0;
}

// Moisture coefficients of WRF calc_cq (cqu, cqv at U/V points, cqw at w levels 1..nz-1),
// precomputed per stage by LargeStep::prep in WRF's summation order. on = false (dry) gives
// exactly cqu = cqv = cq1 = 1, cqw = 0. cq1 = 1/(1+cqw) as pg_buoy_w stores it.
struct MoistCoef {
  View3<const real> u, v, w; bool on = false;
  CPPWRF_HD real cqu(int i, int j, int k) const { return on ? u(i, j, k) : real(1); }
  CPPWRF_HD real cqv(int i, int j, int k) const { return on ? v(i, j, k) : real(1); }
  CPPWRF_HD real cqw(int i, int j, int k) const { return on ? w(i, j, k) : real(0); }
  CPPWRF_HD real cq1(int i, int j, int k) const { return real(1) / (real(1) + cqw(i, j, k)); }
};

// Constant-K diffusion of a mass-point scalar f - base (WRF horizontal_diffusion[_3dmp] 'm' and
// vertical_diffusion[_3dmp|_mp], diff_opt = 1), zero-flux bottom and top. Mass-coupled output.
struct KDiffM {
  View3<const real> f, alt; BView3<const real> base; View2<const real> mut; VCoef c;
  View3<real> out; real kh, kv, rdx, rdy; int nz; bool do_y, hbase, vbase;  // base removed in h / v terms
  CPPWRF_HD real d(int i, int j, int k) const { return hbase ? f(i, j, k) - base(i, j, k) : f(i, j, k); }
  CPPWRF_HD real dv(int i, int j, int k) const { return vbase ? f(i, j, k) - base(i, j, k) : f(i, j, k); }
  CPPWRF_HD real m(int i, int j, int k) const { return c.c1h(k) * mut(i, j) + c.c2h(k); }
  CPPWRF_HD real vf(int i, int j, int k) const {
    // WRF vertical_diffusion(_3dmp): vflux(0) = vflux(1) (no vertical diffusion tendency in the
    // lowest level), vflux(top) = 0
    if (k == nz) return 0;
    if (k == 0) k = 1;
    return kv * c.rdn(k) / (real(0.5) * (alt(i, j, k) + alt(i, j, k - 1))) * (dv(i, j, k) - dv(i, j, k - 1));
  }
  CPPWRF_HD void operator()(int i, int j, int k) const {
    constexpr real H = real(0.5), Gr = real(cst::g);
    real r = rdx * rdx * kh *
             (H * (m(i + 1, j, k) + m(i, j, k)) * (d(i + 1, j, k) - d(i, j, k)) -
              H * (m(i, j, k) + m(i - 1, j, k)) * (d(i, j, k) - d(i - 1, j, k)));
    if (do_y)
      r += rdy * rdy * kh *
           (H * (m(i, j + 1, k) + m(i, j, k)) * (d(i, j + 1, k) - d(i, j, k)) -
            H * (m(i, j, k) + m(i, j - 1, k)) * (d(i, j, k) - d(i, j - 1, k)));
    r += Gr * Gr / (m(i, j, k) * alt(i, j, k)) * c.rdnw(k) * (vf(i, j, k + 1) - vf(i, j, k));
    out(i, j, k) += r;  // accumulates into the held tendency (zeroed at stage 1)
  }
};

}  // namespace cppwrf::dyn::kern
