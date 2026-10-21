#pragma once
// MYNN-EDMF PBL (WRF v4.8.0 phys/MYNN-EDMF, bl_pbl_physics = 5): module_bl_mynnedmf_common's constants and
// the helper functions of module_bl_mynnedmf (esat_blend, qsat_blend, xl_blend, tridiag2, get_pblh,
// scale_aware). REAL(kind_phys) is REAL*8 in WRF's -r8 build; its DOUBLE PRECISION locals are REAL*8 only in
// the -r8 -fdefault-double-8 build (the reference of this port, as for Thompson). x**2 is x*x, x**real pow.
#include <cmath>

#include "core/backend.hpp"
#include "core/precision.hpp"
#include "phys/driver/wrf_constants.hpp"

namespace cppwrf::phys::mynnedmf {

namespace wc = cppwrf::phys::wrfc;

// module_model_constants
inline constexpr real cp = wc::cp, r_d = wc::r_d, r_v = wc::r_v, rcp = wc::rcp, xlv = wc::xlv, xlf = wc::xlf,
                      xls = wc::xls, p1000mb = wc::p1000mb, karman = wc::karman, ep_2 = wc::ep_2;
inline constexpr real rvovrd = wc::r_v / wc::r_d;
inline constexpr real p608 = rvovrd - real(1.);
inline constexpr real cpv = real(4.) * wc::r_v;
inline constexpr real cliq = real(4190.), cice = real(2106.);

// module_bl_mynnedmf_common
inline constexpr real grav = wc::g;
inline constexpr real t0c = wc::svpt0;
inline constexpr real tref = real(300.0);
inline constexpr real tice = real(238.0);
inline constexpr real tliq = real(269.0);
inline constexpr real ep_3 = real(1.) - ep_2;
inline constexpr real gtr = grav / tref;
inline constexpr real tv0 = p608 * tref;
inline constexpr real xlscp = (xlv + xlf) / cp;
inline constexpr real xlvcp = xlv / cp;
inline constexpr real cphm_st = real(5.0), cphm_unst = real(16.0);
inline constexpr real cphh_st = real(5.0), cphh_unst = real(16.0);
inline constexpr real qmin = real(0.0);
inline constexpr real zmax = real(1.0);
inline constexpr real Sqfac = real(3.0);
inline constexpr real qkemin = real(1.e-5);
inline constexpr real pr = real(0.74);
inline constexpr real g1 = real(0.235);
inline constexpr real b1 = real(24.0);
inline constexpr real b2 = real(15.0);
inline constexpr real c2 = real(0.729);
inline constexpr real c3 = real(0.340);
inline constexpr real c5 = real(0.2);
inline constexpr real a1 = b1 * (real(1.0) - real(3.0) * g1) / real(6.0);
inline constexpr real c1 = g1 - real(1.0) / (real(3.0) * a1 * real(2.88449914061481660));
inline constexpr real a2 = a1 * (g1 - c1) / (g1 * pr);
inline constexpr real g2 = b2 / b1 * (real(1.0) - c3) + real(2.0) * a1 / b1 * (real(3.0) - real(2.0) * c2);
inline constexpr real cc2 = real(1.0) - c2;
inline constexpr real cc3 = real(1.0) - c3;
inline constexpr real e1c = real(3.0) * a2 * b2 * cc3;
inline constexpr real e2c = real(9.0) * a1 * a2 * cc2;
inline constexpr real e3c = real(9.0) * a2 * a2 * cc2 * (real(1.0) - c5);
inline constexpr real e4c = real(12.0) * a1 * a2 * cc2;
inline constexpr real e5c = real(6.0) * a1 * a1;
inline constexpr real zero = 0, one = 1, two = 2, three = 3, four = 4, five = 5, seven = 7, nine = 9, ten = 10,
                      twenty = 20, thirty = 30, forty = 40, fifty = 50, hundred = 100;
inline constexpr real p01 = real(0.01), p1 = real(0.1), p2 = real(0.2), p25 = real(0.25), p3 = real(0.3),
                      p333 = real(1.0) / real(3.0), p4 = real(0.4), p5 = real(0.5), p666 = real(2.0) / real(3.0),
                      p95 = real(0.95);

// esat_blend / qsat_blend: Flatau polynomials, blended between liquid and ice
inline constexpr real J0 = real(.611583699E03), J1 = real(.444606896E02), J2 = real(.143177157E01),
                      J3 = real(.264224321E-1), J4 = real(.299291081E-3), J5 = real(.203154182E-5),
                      J6 = real(.702620698E-8), J7 = real(.379534310E-11), J8 = real(-.321582393E-13);
inline constexpr real K0 = real(.609868993E03), K1 = real(.499320233E02), K2 = real(.184672631E01),
                      K3 = real(.402737184E-1), K4 = real(.565392987E-3), K5 = real(.521693933E-5),
                      K6 = real(.307839583E-7), K7 = real(.105785160E-9), K8 = real(.161444444E-12);

CPPWRF_HD inline real esl_poly(real XC) {
  return J0 + XC * (J1 + XC * (J2 + XC * (J3 + XC * (J4 + XC * (J5 + XC * (J6 + XC * (J7 + XC * J8)))))));
}
CPPWRF_HD inline real esi_poly(real XC) {
  return K0 + XC * (K1 + XC * (K2 + XC * (K3 + XC * (K4 + XC * (K5 + XC * (K6 + XC * (K7 + XC * K8)))))));
}

CPPWRF_HD inline real esat_blend(real t) {
  const real XC = c_max(real(-80.), t - t0c);
  if (t >= (t0c - real(6.))) return esl_poly(XC);
  if (t <= tice) return esi_poly(XC);
  const real ESL = esl_poly(XC), ESI = esi_poly(XC);
  const real chi = ((t0c - real(6.)) - t) / ((t0c - real(6.)) - tice);
  return (one - chi) * ESL + chi * ESI;
}

CPPWRF_HD inline real qsat_blend(real t, real P) {
  const real XC = c_max(real(-80.), t - t0c);
  if (t >= (t0c - real(6.))) {
    real ESL = esl_poly(XC);
    ESL = c_min(ESL, P * real(0.15));
    return real(0.622) * ESL / c_max(P - ESL, real(1e-5));
  }
  if (t <= tice) {
    real ESI = esi_poly(XC);
    ESI = c_min(ESI, P * real(0.15));
    return real(0.622) * ESI / c_max(P - ESI, real(1e-5));
  }
  real ESL = esl_poly(XC);
  ESL = c_min(ESL, P * real(0.15));
  real ESI = esi_poly(XC);
  ESI = c_min(ESI, P * real(0.15));
  const real RSLF = real(0.622) * ESL / c_max(P - ESL, real(1e-5));
  const real RSIF = real(0.622) * ESI / c_max(P - ESI, real(1e-5));
  const real chi = ((t0c - real(6.)) - t) / ((t0c - real(6.)) - tice);
  return (one - chi) * RSLF + chi * RSIF;
}

CPPWRF_HD inline real xl_blend(real t) {
  if (t >= t0c) return xlv + (cpv - cliq) * (t - t0c);
  if (t <= tice) return xls + (cpv - cice) * (t - t0c);
  const real xlvt = xlv + (cpv - cliq) * (t - t0c);
  const real xlst = xls + (cpv - cice) * (t - t0c);
  const real chi = (t0c - t) / (t0c - tice);
  return (one - chi) * xlvt + chi * xlst;
}

// tridiag2 (1-based a, b, c, d, x of kte elements; a(1) is not used)
template <int K>
CPPWRF_HD inline void tridiag2(int kte, const real* a, const real* b, const real* c, const real* d, real* x) {
  real cp1[K + 2], dp1[K + 2];
  cp1[1] = c[1] / b[1];
  dp1[1] = d[1] / b[1];
  for (int k = 2; k <= kte; ++k) {
    const real m = b[k] - cp1[k - 1] * a[k];
    cp1[k] = c[k] / m;
    dp1[k] = (d[k] - dp1[k - 1] * a[k]) / m;
  }
  x[kte] = dp1[kte];
  for (int k = kte - 1; k >= 1; --k) x[k] = dp1[k] - cp1[k] * x[k + 1];
}

// GET_PBLH (stable_method = 1)
CPPWRF_HD inline void get_pblh(int kte, real& pblh, const real* thv1, const real* zw1, const real* dz1, real landsea,
                               real ust, int& kpbl) {
  constexpr int kts = 1;
  constexpr real sbl_lim = real(200.), sbl_damp = real(170.);
  kpbl = 2;
  int k = kts + 1;
  int kthv = 1;
  real minthv = real(9.E9);
  while (zw1[k] <= real(200.)) {
    if (minthv > thv1[k]) {
      minthv = thv1[k];
      kthv = k;
    }
    k = k + 1;
  }
  (void)kthv;
  const real delt_thv = ((landsea - real(1.5)) >= zero) ? real(0.75) : real(1.25);
  pblh = zero;
  for (k = kts + 1; k <= kte - 1; ++k) {
    if (thv1[k] >= (minthv + delt_thv))
      pblh = zw1[k] - dz1[k - 1] * c_min((thv1[k] - (minthv + delt_thv)) / c_max(thv1[k] - thv1[k - 1], real(1E-6)), one);
    if (k == kte - 1) pblh = zw1[kts + 1];  // exit safeguard
    if (pblh > zero) break;
  }
  real wt = p5 * std::tanh((pblh - hundred) / hundred) + p5;
  const real cpblh = real(400.) * (one - wt) + real(700.) * wt;
  real pblh_stable = c_max(ust * cpblh, ten);
  pblh_stable = c_max(pblh_stable, pblh - sbl_lim);
  pblh_stable = c_min(pblh_stable, pblh + sbl_lim);
  wt = p5 * std::tanh((pblh - sbl_lim) / sbl_damp) + p5;
  pblh = pblh_stable * (one - wt) + pblh * wt;
  for (k = kts + 1; k <= kte - 1; ++k) {
    if (zw1[k] >= pblh) {
      kpbl = k - 1;
      break;
    }
  }
}

// SCALE_AWARE
CPPWRF_HD inline void scale_aware(real dx, real pblh, real& Psig_bl, real& Psig_shcu) {
  real dxdh = c_max(real(2.5) * dx, ten) / c_min(pblh, real(3000.));
  Psig_bl = ((dxdh * dxdh) + real(0.106) * rpow(dxdh, real(0.667))) /
            ((dxdh * dxdh) + real(0.066) * rpow(dxdh, real(0.667)) + real(0.071));
  dxdh = c_max(real(2.5) * dx, ten) / c_min(pblh + real(500.), real(3500.));
  Psig_shcu = ((dxdh * dxdh) + real(0.145) * rpow(dxdh, real(0.667))) /
              ((dxdh * dxdh) + real(0.172) * rpow(dxdh, real(0.667)) + real(0.170));
  if (Psig_bl > one) Psig_bl = one;
  if (Psig_bl < zero) Psig_bl = zero;
  if (Psig_shcu > one) Psig_shcu = one;
  if (Psig_shcu < zero) Psig_shcu = zero;
}

}  // namespace cppwrf::phys::mynnedmf
