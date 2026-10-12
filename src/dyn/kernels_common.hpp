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
// first-order upwind flux with Courant number cr (advect_scalar_pd flux_upwind)
CPPWRF_HD inline real flux_upwind(real qm1, real q0, real cr) {
  const real a = cr < real(0) ? -cr : cr;
  const real lo = cr + a, hi = cr - a;
  return real(0.5) * (lo < real(1) ? lo : real(1)) * qm1 + real(0.5) * (hi > real(-1) ? hi : real(-1)) * q0;
}

// Moisture coefficients of WRF calc_cq / pg_buoy_w from qtot (sum of moist species, mass
// points). on = false (dry) gives exactly cqu = cqv = cq1 = 1, cqw = cq2 = 0.
struct MoistCoef {
  View3<const real> qt; bool on = false;
  CPPWRF_HD real cqu(int i, int j, int k) const {
    return on ? real(1) / (real(1) + real(0.5) * (qt(i, j, k) + qt(i - 1, j, k))) : real(1);
  }
  CPPWRF_HD real cqv(int i, int j, int k) const {
    return on ? real(1) / (real(1) + real(0.5) * (qt(i, j, k) + qt(i, j - 1, k))) : real(1);
  }
  CPPWRF_HD real cqw(int i, int j, int k) const {  // w levels 1..nz-1 (calc_cq value)
    return on ? real(0.5) * (qt(i, j, k) + qt(i, j, k - 1)) : real(0);
  }
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
    if (k == 0 || k == nz) return 0;
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
    out(i, j, k) = r;
  }
};

}  // namespace cppwrf::dyn::kern
