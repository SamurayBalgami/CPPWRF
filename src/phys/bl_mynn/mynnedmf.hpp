#pragma once
// MYNN-EDMF PBL scheme (WRF v4.8.0 phys/MYNN-EDMF, bl_pbl_physics = 5): mynnedmf, mynn_tendencies,
// moisture_check, retrieve_exchange_coeffs and the column part of mynnedmf_driver (pre/post run) for the options
// in Options (bl_mynn_edmf 0 / 1, edmf_mom, edmf_tke, cloudpdf, mixlength, closure 2.5 / 2.6 / 3.0, cloudmix,
// mixqt, ess, tkeadvect) and tke_budget; bl_mynn_edmf_dd = 0, bl_mynn_mixnumcon = 0, spp_pbl = 0, no chemistry
// (bl_mynn_mixscalars mixes WRF's unconnected zero scalar only); a microphysics with qc, qi and qs (flag_qs: WRF
// passes the snow content uninitialised otherwise). Column state in specific contents, 1-based (k = 1 .. kte).
#include "phys/bl_mynn/mynnedmf_mf.hpp"
#include "phys/bl_mynn/mynnedmf_turb.hpp"

namespace cppwrf::phys::mynnedmf {

inline constexpr int kMaxLevels = 128;  // Column<kMaxLevels> in the model
inline constexpr real ckmod = real(1.0);
inline constexpr real taue = real(80.);

// moisture_check (1-based columns of kte levels)
CPPWRF_HD inline void moisture_check(int kte, real delt, const real* dp, const real* exner, real* qv, real* qc,
                                     real* qi, real* qs, real* th, real* dqv, real* dqc, real* dqi, real* dqs,
                                     real* dth) {
  constexpr real qvmin = real(1e-20), qcmin = zero, qimin = zero;
  real dqv2 = zero;
  for (int k = kte; k >= 1; --k) {
    const real dqc2 = c_max(zero, qcmin - qc[k]);
    const real dqi2 = c_max(zero, qimin - qi[k]);
    const real dqs2 = c_max(zero, qimin - qs[k]);
    dqc[k] = dqc[k] + dqc2 / delt;
    dqi[k] = dqi[k] + dqi2 / delt;
    dqs[k] = dqs[k] + dqs2 / delt;
    dqv[k] = dqv[k] - (dqc2 + dqi2 + dqs2) / delt;
    dth[k] = dth[k] + xlvcp / exner[k] * (dqc2 / delt) + xlscp / exner[k] * ((dqi2 + dqs2) / delt);
    qc[k] = qc[k] + dqc2;
    qi[k] = qi[k] + dqi2;
    qs[k] = qs[k] + dqs2;
    qv[k] = qv[k] - dqc2 - dqi2 - dqs2;
    th[k] = th[k] + xlvcp / exner[k] * dqc2 + xlscp / exner[k] * (dqi2 + dqs2);
    dqv2 = c_max(zero, qvmin - qv[k]);
    dqv[k] = dqv[k] + dqv2 / delt;
    qv[k] = qv[k] + dqv2;
    if (k != 1) {
      qv[k - 1] = qv[k - 1] - dqv2 * dp[k] / dp[k - 1];
      dqv[k - 1] = dqv[k - 1] - dqv2 * dp[k] / dp[k - 1] / delt;
    }
    qv[k] = c_max(qv[k], qvmin);
    qc[k] = c_max(qc[k], qcmin);
    qi[k] = c_max(qi[k], qimin);
    qs[k] = c_max(qs[k], qimin);
  }
  if (dqv2 > real(1.e-20)) {
    real sum = zero;
    for (int k = 1; k <= kte; ++k)
      if (qv[k] > two * qvmin) sum = sum + qv[k] * dp[k];
    const real aa = dqv2 * dp[1] / c_max(real(1.e-20), sum);
    if (aa < p5) {
      for (int k = 1; k <= kte; ++k)
        if (qv[k] > two * qvmin) {
          const real dum = aa * qv[k];
          qv[k] = qv[k] - dum;
          dqv[k] = dqv[k] - dum / delt;
        }
    }
  }
}

// mynn_tendencies, implicit mass-flux form (bl_mynn_edmf 0 or 1), bl_mynn_mixqt 0 / 1, bl_mynn_cloudmix 0 / 1,
// flag_qc = flag_qi = .true., no number concentrations / aerosols / ozone; the snow arguments are zero
// (mynnedmf passes kzero1); no downdraft fluxes (sd_* = 0)
template <int K>
CPPWRF_HD inline void mynn_tendencies(int kte, real delt, const real* dz, const real* rho, const real* u, const real* v,
                                      const real* th, const real* tk, const real* qv, real psfc, const real* p,
                                      const real* exner, real* thl, const real* sqv, const real* sqc, const real* sqi,
                                      const real* sqw, real ust, real flt, real flq, real flqv, real flqc, real wspd,
                                      real uoce, real voce, const real* tcd, const real* qcd, const real* dfm,
                                      const real* dfh, real* du, real* dv, real* dth, real* dqv, real* dqc, real* dqi,
                                      real* dqs, const real* diss_heat, const real* s_aw1, const real* s_awthl1,
                                      const real* s_awqt1, const real* s_awqv1, const real* s_awqc1,
                                      const real* s_awu1, const real* s_awv1, int bl_mynn_cloudmix, int bl_mynn_mixqt,
                                      int bl_mynn_edmf, int bl_mynn_edmf_mom) {
  constexpr int kts = 1;
  real dtz[K + 2], delp[K + 2], rhoinv[K + 2], a[K + 2], b[K + 2], c[K + 2], d[K + 2], x[K + 2], sqv2[K + 2],
      sqc2[K + 2], sqi2[K + 2], sqs2[K + 2], sqw2[K + 2], rhoz[K + 3], khdz[K + 3], kmdz[K + 3];
  // zero downdraft fluxes and the zero environmental subsidence / detrainment tendencies (env_subs = .false.)
  constexpr real sd = zero, sub = zero, det = zero;
  const real onoff = (bl_mynn_edmf_mom == 0) ? zero : one;
  const real ustovwsp = (wspd <= real(1e-6)) ? zero : ust / wspd;
  const real rhosfc = psfc / (r_d * (tk[kts] + p608 * qv[kts]));
  dtz[kts] = delt / dz[kts];
  rhoz[kts] = rho[kts];
  rhoinv[kts] = one / rho[kts];
  khdz[kts] = rhoz[kts] * dfh[kts];
  kmdz[kts] = rhoz[kts] * dfm[kts];
  for (int k = kts + 1; k <= kte; ++k) {
    dtz[k] = delt / dz[k];
    rhoz[k] = (rho[k] * dz[k - 1] + rho[k - 1] * dz[k]) / (dz[k - 1] + dz[k]);
    rhoz[k] = c_max(rhoz[k], real(1E-4));
    rhoinv[k] = one / c_max(rho[k], real(1E-4));
    khdz[k] = rhoz[k] * dfh[k];
    kmdz[k] = rhoz[k] * dfm[k];
  }
  rhoz[kte + 1] = rhoz[kte];
  khdz[kte + 1] = rhoz[kte + 1] * dfh[kte];
  kmdz[kte + 1] = rhoz[kte + 1] * dfm[kte];
  for (int k = kts; k <= kte; ++k) delp[k] = rho[k] * grav * dz[k];
  if (delp[kts] < p333 * delp[kts + 1]) delp[kts] = p333 * delp[kts + 1];
  if (bl_mynn_edmf == 1) {
    for (int k = kts + 1; k <= kte - 1; ++k) {
      khdz[k] = c_max(khdz[k], p5 * (s_aw1[k] + sd));
      khdz[k] = c_max(khdz[k], -p5 * (s_aw1[k] - s_aw1[k + 1]) - p5 * (sd - sd));
      kmdz[k] = c_max(kmdz[k], p5 * (s_aw1[k] + sd));
      kmdz[k] = c_max(kmdz[k], -p5 * (s_aw1[k] - s_aw1[k + 1]) - p5 * (sd - sd));
    }
  }
  for (int k = kts; k <= kte; ++k) dth[k] = zero;

  // u
  {
    int k = kts;
    a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
    b[k] = one + dtz[k] * (kmdz[k + 1] + rhosfc * ust * ustovwsp) * rhoinv[k] -
           p5 * dtz[k] * rhoinv[k] * s_aw1[k + 1] * onoff - p5 * dtz[k] * rhoinv[k] * sd * onoff;
    c[k] = -dtz[k] * kmdz[k + 1] * rhoinv[k] - p5 * dtz[k] * rhoinv[k] * s_aw1[k + 1] * onoff -
           p5 * dtz[k] * rhoinv[k] * sd * onoff;
    d[k] = u[k] + dtz[k] * uoce * (ust * ust) / wspd - dtz[k] * rhoinv[k] * s_awu1[k + 1] * onoff +
           dtz[k] * rhoinv[k] * sd * onoff + sub * delt + det * delt;
    for (k = kts + 1; k <= kte - 1; ++k) {
      a[k] = -dtz[k] * kmdz[k] * rhoinv[k] + p5 * dtz[k] * rhoinv[k] * s_aw1[k] * onoff +
             p5 * dtz[k] * rhoinv[k] * sd * onoff;
      b[k] = one + dtz[k] * (kmdz[k] + kmdz[k + 1]) * rhoinv[k] +
             p5 * dtz[k] * rhoinv[k] * (s_aw1[k] - s_aw1[k + 1]) * onoff + p5 * dtz[k] * rhoinv[k] * (sd - sd) * onoff;
      c[k] = -dtz[k] * kmdz[k + 1] * rhoinv[k] - p5 * dtz[k] * rhoinv[k] * s_aw1[k + 1] * onoff -
             p5 * dtz[k] * rhoinv[k] * sd * onoff;
      d[k] = u[k] + dtz[k] * rhoinv[k] * (s_awu1[k] - s_awu1[k + 1]) * onoff -
             dtz[k] * rhoinv[k] * (sd - sd) * onoff + sub * delt + det * delt;
    }
    a[kte] = zero;
    b[kte] = one;
    c[kte] = zero;
    d[kte] = u[kte];
    tridiag2<K>(kte, a, b, c, d, x);
    for (k = kts; k <= kte; ++k) du[k] = (x[k] - u[k]) / delt;
  }
  // v
  {
    int k = kts;
    a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
    b[k] = one + dtz[k] * (kmdz[k + 1] + rhosfc * (ust * ust) / wspd) * rhoinv[k] -
           p5 * dtz[k] * rhoinv[k] * s_aw1[k + 1] * onoff - p5 * dtz[k] * rhoinv[k] * sd * onoff;
    c[k] = -dtz[k] * kmdz[k + 1] * rhoinv[k] - p5 * dtz[k] * rhoinv[k] * s_aw1[k + 1] * onoff -
           p5 * dtz[k] * rhoinv[k] * sd * onoff;
    d[k] = v[k] + dtz[k] * voce * (ust * ust) / wspd - dtz[k] * rhoinv[k] * s_awv1[k + 1] * onoff +
           dtz[k] * rhoinv[k] * sd * onoff + sub * delt + det * delt;
    for (k = kts + 1; k <= kte - 1; ++k) {
      a[k] = -dtz[k] * kmdz[k] * rhoinv[k] + p5 * dtz[k] * rhoinv[k] * s_aw1[k] * onoff +
             p5 * dtz[k] * rhoinv[k] * sd * onoff;
      b[k] = one + dtz[k] * (kmdz[k] + kmdz[k + 1]) * rhoinv[k] +
             p5 * dtz[k] * rhoinv[k] * (s_aw1[k] - s_aw1[k + 1]) * onoff + p5 * dtz[k] * rhoinv[k] * (sd - sd) * onoff;
      c[k] = -dtz[k] * kmdz[k + 1] * rhoinv[k] - p5 * dtz[k] * rhoinv[k] * s_aw1[k + 1] * onoff -
             p5 * dtz[k] * rhoinv[k] * sd * onoff;
      d[k] = v[k] + dtz[k] * rhoinv[k] * (s_awv1[k] - s_awv1[k + 1]) * onoff -
             dtz[k] * rhoinv[k] * (sd - sd) * onoff + sub * delt + det * delt;
    }
    a[kte] = zero;
    b[kte] = one;
    c[kte] = zero;
    d[kte] = v[kte];
    tridiag2<K>(kte, a, b, c, d, x);
    for (k = kts; k <= kte; ++k) dv[k] = (x[k] - v[k]) / delt;
  }
  // thl
  {
    int k = kts;
    a[k] = -dtz[k] * khdz[k] * rhoinv[k];
    b[k] = real(1.) + dtz[k] * (khdz[k + 1] + khdz[k]) * rhoinv[k] - p5 * dtz[k] * rhoinv[k] * s_aw1[k + 1] -
           p5 * dtz[k] * rhoinv[k] * sd;
    c[k] = -dtz[k] * khdz[k + 1] * rhoinv[k] - p5 * dtz[k] * rhoinv[k] * s_aw1[k + 1] - p5 * dtz[k] * rhoinv[k] * sd;
    d[k] = thl[k] + dtz[k] * rhosfc * flt * rhoinv[k] + tcd[k] * delt - dtz[k] * rhoinv[k] * s_awthl1[k + 1] +
           dtz[k] * rhoinv[k] * sd + diss_heat[k] * delt + sub * delt + det * delt;
    for (k = kts + 1; k <= kte - 1; ++k) {
      a[k] = -dtz[k] * khdz[k] * rhoinv[k] + p5 * dtz[k] * rhoinv[k] * s_aw1[k] + p5 * dtz[k] * rhoinv[k] * sd;
      b[k] = real(1.) + dtz[k] * (khdz[k] + khdz[k + 1]) * rhoinv[k] + p5 * dtz[k] * rhoinv[k] * (s_aw1[k] - s_aw1[k + 1]) +
             p5 * dtz[k] * rhoinv[k] * (sd - sd);
      c[k] = -dtz[k] * khdz[k + 1] * rhoinv[k] - p5 * dtz[k] * rhoinv[k] * s_aw1[k + 1] - p5 * dtz[k] * rhoinv[k] * sd;
      d[k] = thl[k] + tcd[k] * delt + dtz[k] * rhoinv[k] * (s_awthl1[k] - s_awthl1[k + 1]) -
             dtz[k] * rhoinv[k] * (sd - sd) + diss_heat[k] * delt + sub * delt + det * delt;
    }
    a[kte] = zero;
    b[kte] = one;
    c[kte] = zero;
    d[kte] = thl[kte];
    tridiag2<K>(kte, a, b, c, d, x);
    for (k = kts; k <= kte; ++k) thl[k] = x[k];
  }
  // sqw (bl_mynn_mixqt > 0)
  if (bl_mynn_mixqt > 0) {
    int k = kts;
    a[k] = -dtz[k] * khdz[k] * rhoinv[k];
    b[k] = one + dtz[k] * (khdz[k + 1] + khdz[k]) * rhoinv[k] - p5 * dtz[k] * rhoinv[k] * s_aw1[k + 1] -
           p5 * dtz[k] * rhoinv[k] * sd;
    c[k] = -dtz[k] * khdz[k + 1] * rhoinv[k] - p5 * dtz[k] * rhoinv[k] * s_aw1[k + 1] - p5 * dtz[k] * rhoinv[k] * sd;
    d[k] = sqw[k] + dtz[k] * rhosfc * flq * rhoinv[k] + qcd[k] * delt - dtz[k] * rhoinv[k] * s_awqt1[k + 1] +
           dtz[k] * rhoinv[k] * sd;
    for (k = kts + 1; k <= kte - 1; ++k) {
      a[k] = -dtz[k] * khdz[k] * rhoinv[k] + p5 * dtz[k] * rhoinv[k] * s_aw1[k] + p5 * dtz[k] * rhoinv[k] * sd;
      b[k] = one + dtz[k] * (khdz[k] + khdz[k + 1]) * rhoinv[k] + p5 * dtz[k] * rhoinv[k] * (s_aw1[k] - s_aw1[k + 1]) +
             p5 * dtz[k] * rhoinv[k] * (sd - sd);
      c[k] = -dtz[k] * khdz[k + 1] * rhoinv[k] - p5 * dtz[k] * rhoinv[k] * s_aw1[k + 1] - p5 * dtz[k] * rhoinv[k] * sd;
      d[k] = sqw[k] + qcd[k] * delt + dtz[k] * rhoinv[k] * (s_awqt1[k] - s_awqt1[k + 1]) -
             dtz[k] * rhoinv[k] * (sd - sd);
    }
    a[kte] = zero;
    b[kte] = one;
    c[kte] = zero;
    d[kte] = sqw[kte];
    tridiag2<K>(kte, a, b, c, d, sqw2);
  } else {
    for (int k = kts; k <= kte; ++k) sqw2[k] = sqw[k];
  }
  // sqc (bl_mynn_mixqt = 0; mixed when bl_mynn_cloudmix > 0, flag_qc)
  if (bl_mynn_mixqt == 0 && bl_mynn_cloudmix == 0) {
    for (int k = kts; k <= kte; ++k) sqc2[k] = sqc[k];
  } else if (bl_mynn_mixqt == 0) {
    int k = kts;
    a[k] = -dtz[k] * khdz[k] * rhoinv[k];
    b[k] = one + dtz[k] * (khdz[k + 1] + khdz[k]) * rhoinv[k] - p5 * dtz[k] * rhoinv[k] * s_aw1[k + 1] -
           p5 * dtz[k] * rhoinv[k] * sd;
    c[k] = -dtz[k] * khdz[k + 1] * rhoinv[k] - p5 * dtz[k] * rhoinv[k] * s_aw1[k + 1] - p5 * dtz[k] * rhoinv[k] * sd;
    d[k] = sqc[k] + dtz[k] * rhosfc * flqc * rhoinv[k] + qcd[k] * delt - dtz[k] * rhoinv[k] * s_awqc1[k + 1] +
           dtz[k] * rhoinv[k] * sd + det * delt;
    for (k = kts + 1; k <= kte - 1; ++k) {
      a[k] = -dtz[k] * khdz[k] * rhoinv[k] + p5 * dtz[k] * rhoinv[k] * s_aw1[k] + p5 * dtz[k] * rhoinv[k] * sd;
      b[k] = one + dtz[k] * (khdz[k] + khdz[k + 1]) * rhoinv[k] + p5 * dtz[k] * rhoinv[k] * (s_aw1[k] - s_aw1[k + 1]) +
             p5 * dtz[k] * rhoinv[k] * (sd - sd);
      c[k] = -dtz[k] * khdz[k + 1] * rhoinv[k] - p5 * dtz[k] * rhoinv[k] * s_aw1[k + 1] - p5 * dtz[k] * rhoinv[k] * sd;
      d[k] = sqc[k] + qcd[k] * delt + dtz[k] * rhoinv[k] * (s_awqc1[k] - s_awqc1[k + 1]) -
             dtz[k] * rhoinv[k] * (sd - sd) + det * delt;
    }
    a[kte] = zero;
    b[kte] = one;
    c[kte] = zero;
    d[kte] = sqc[kte];
    tridiag2<K>(kte, a, b, c, d, sqc2);
  }
  // sqv (bl_mynn_mixqt = 0)
  if (bl_mynn_mixqt > 0) {
    for (int k = kts; k <= kte; ++k) sqv2[k] = sqv[k];
  } else {
    int k = kts;
    real qvflux = flqv;
    if (qvflux < zero) qvflux = c_max(qvflux, (c_min(real(0.9) * sqv[kts] - real(1e-8), zero) / dtz[kts]));
    a[k] = -dtz[k] * khdz[k] * rhoinv[k];
    b[k] = one + dtz[k] * (khdz[k + 1] + khdz[k]) * rhoinv[k] - p5 * dtz[k] * rhoinv[k] * s_aw1[k + 1] -
           p5 * dtz[k] * rhoinv[k] * sd;
    c[k] = -dtz[k] * khdz[k + 1] * rhoinv[k] - p5 * dtz[k] * rhoinv[k] * s_aw1[k + 1] - p5 * dtz[k] * rhoinv[k] * sd;
    d[k] = sqv[k] + dtz[k] * rhosfc * qvflux * rhoinv[k] + qcd[k] * delt - dtz[k] * rhoinv[k] * s_awqv1[k + 1] +
           dtz[k] * rhoinv[k] * sd + sub * delt + det * delt;
    for (k = kts + 1; k <= kte - 1; ++k) {
      a[k] = -dtz[k] * khdz[k] * rhoinv[k] + p5 * dtz[k] * rhoinv[k] * s_aw1[k] + p5 * dtz[k] * rhoinv[k] * sd;
      b[k] = one + dtz[k] * (khdz[k] + khdz[k + 1]) * rhoinv[k] + p5 * dtz[k] * rhoinv[k] * (s_aw1[k] - s_aw1[k + 1]) +
             p5 * dtz[k] * rhoinv[k] * (sd - sd);
      c[k] = -dtz[k] * khdz[k + 1] * rhoinv[k] - p5 * dtz[k] * rhoinv[k] * s_aw1[k + 1] - p5 * dtz[k] * rhoinv[k] * sd;
      d[k] = sqv[k] + qcd[k] * delt + dtz[k] * rhoinv[k] * (s_awqv1[k] - s_awqv1[k + 1]) -
             dtz[k] * rhoinv[k] * (sd - sd) + sub * delt + det * delt;
    }
    a[kte] = zero;
    b[kte] = one;
    c[kte] = zero;
    d[kte] = sqv[kte];
    tridiag2<K>(kte, a, b, c, d, sqv2);
  }
  // sqi (bl_mynn_cloudmix > 0, flag_qi): only downdraft fluxes
  if (bl_mynn_cloudmix == 0) {
    for (int k = kts; k <= kte; ++k) sqi2[k] = sqi[k];
  } else {
    int k = kts;
    a[k] = -dtz[k] * khdz[k] * rhoinv[k];
    b[k] = one + dtz[k] * (khdz[k + 1] + khdz[k]) * rhoinv[k] - p5 * dtz[k] * rhoinv[k] * sd;
    c[k] = -dtz[k] * khdz[k + 1] * rhoinv[k] - p5 * dtz[k] * rhoinv[k] * sd;
    d[k] = sqi[k] + dtz[k] * rhoinv[k] * sd;
    for (k = kts + 1; k <= kte - 1; ++k) {
      a[k] = -dtz[k] * khdz[k] * rhoinv[k] + p5 * dtz[k] * rhoinv[k] * sd;
      b[k] = one + dtz[k] * (khdz[k] + khdz[k + 1]) * rhoinv[k] + p5 * dtz[k] * rhoinv[k] * (sd - sd);
      c[k] = -dtz[k] * khdz[k + 1] * rhoinv[k] - p5 * dtz[k] * rhoinv[k] * sd;
      d[k] = sqi[k] - dtz[k] * rhoinv[k] * (sd - sd);
    }
    a[kte] = zero;
    b[kte] = one;
    c[kte] = zero;
    d[kte] = sqi[kte];
    tridiag2<K>(kte, a, b, c, d, sqi2);
  }
  for (int k = kts; k <= kte; ++k) sqs2[k] = zero;  // snow not mixed: sqs2 = sqs = kzero1
  if (bl_mynn_mixqt > 0) {  // partition the mixed total water
    for (int k = kts; k <= kte; ++k) {
      const real th_new = thl[k] + xlvcp / exner[k] * sqc[k] + xlscp / exner[k] * sqi[k];
      const real t = th_new * exner[k];
      const real qsat = qsat_blend(t, p[k]);
      if (sqc[k] > zero || sqi[k] > zero) {
        sqv2[k] = c_min(sqw2[k], qsat);
        const real portion_qc = sqc[k] / (sqc[k] + sqi[k]);
        const real portion_qi = sqi[k] / (sqc[k] + sqi[k]);
        const real condensate = c_max(sqw2[k] - qsat, zero);
        sqc2[k] = condensate * portion_qc;
        sqi2[k] = condensate * portion_qi;
      } else {
        sqv2[k] = sqw2[k];
        sqi2[k] = zero;
        sqc2[k] = zero;
      }
    }
  }
  for (int k = kts; k <= kte; ++k) dqv[k] = (sqv2[k] - sqv[k]) / delt;
  if (bl_mynn_cloudmix > 0) {
    for (int k = kts; k <= kte; ++k) dqc[k] = (sqc2[k] - sqc[k]) / delt;
    for (int k = kts; k <= kte; ++k) dqi[k] = (sqi2[k] - sqi[k]) / delt;
    for (int k = kts; k <= kte; ++k) dqs[k] = real(0.);
  } else {
    for (int k = kts; k <= kte; ++k) {
      dqc[k] = real(0.);
      dqi[k] = real(0.);
      dqs[k] = real(0.);
    }
  }
  moisture_check(kte, delt, delp, exner, sqv2, sqc2, sqi2, sqs2, thl, dqv, dqc, dqi, dqs, dth);
  for (int k = kts; k <= kte; ++k)
    dth[k] = (thl[k] + xlvcp / exner[k] * sqc[k] + xlscp / exner[k] * (sqi[k]) - th[k]) / delt;
}

// MYNN-EDMF state of one column: the interface arrays of mynnedmf_driver, 1-based, specific contents for
// sqv/sqc/sqi/sqs (mynnedmf_pre_run), tendencies in specific contents (before mynnedmf_post_run)
template <int K>
struct Column {
  int kte;
  // in
  real dz[K + 2], u[K + 2], v[K + 2], w[K + 2], th[K + 2], pres[K + 2], ex[K + 2], rho[K + 2], tk[K + 2];
  real sqv[K + 2], sqc[K + 2], sqi[K + 2], sqs[K + 2];
  // inout
  real qke[K + 2], qke_adv[K + 2], el[K + 2], sh[K + 2], sm[K + 2], kh[K + 2], km[K + 2], tsq[K + 2], qsq[K + 2],
      cov[K + 2], qc_bl[K + 2], qi_bl[K + 2], cldfra_bl[K + 2];
  // out
  real du[K + 2], dv[K + 2], dth[K + 2], dqv[K + 2], dqc[K + 2], dqi[K + 2];
  real edmf_a[K + 2], edmf_w[K + 2], edmf_qt[K + 2], edmf_thl[K + 2], edmf_ent[K + 2], edmf_qc[K + 2];
  real dqke[K + 2], qwt[K + 2], qshear[K + 2], qbuoy[K + 2], qdiss[K + 2];  // tke_budget = 1
  real pblh, maxwidth, maxmf, ztop_plume, excess_h, excess_q, maxmf_dd, maxwidth_dd, maxtkeprod, cldtop_cooling,
      ent_eff;
  int kpbl;
};

// mynnedmf (initflag: first step; restart = cycling = .false.); tke_budget = 1 fills c.dqke, qwt, qshear, qbuoy, qdiss
template <int K>
CPPWRF_HD inline void mynnedmf(Column<K>& c, const Options& op, int tke_budget, int initflag, real delt, real dx,
                               real xland, real ts, real ps, real ust, real hfx, real qfx, real wspd, real uoce,
                               real voce) {
  constexpr int kts = 1;
  const int kte = c.kte;
  const int bl_mynn_edmf = op.edmf, bl_mynn_edmf_mom = op.edmf_mom, bl_mynn_edmf_tke = op.edmf_tke,
            bl_mynn_ess = op.ess;
  const real closure = op.closure;
  real* qwt = tke_budget == 1 ? c.qwt : nullptr;
  real* qdiss = tke_budget == 1 ? c.qdiss : nullptr;
  real* qshear = tke_budget == 1 ? c.qshear : nullptr;
  real* qbuoy = tke_budget == 1 ? c.qbuoy : nullptr;
  real thl1[K + 2], thv1[K + 2], thlv1[K + 2], qv1[K + 2], qc1[K + 2], qi1[K + 2], sqw1[K + 2], thl_tot1[K + 2],
      qc_tot1[K + 2], qi_tot1[K + 2], dfm1[K + 2], dfh1[K + 2], dfq1[K + 2], tcd1[K + 2], qcd1[K + 2], pdk1[K + 2],
      pdt1[K + 2], pdq1[K + 2], pdc1[K + 2], vt1[K + 2], vq1[K + 2], sgm1[K + 2], diss_heat1[K + 2],
      tkeprod_up[K + 2], dqs1[K + 2];
  real s_aw1[K + 3], s_awthl1[K + 3], s_awqt1[K + 3], s_awqv1[K + 3], s_awqc1[K + 3], s_awu1[K + 3], s_awv1[K + 3],
      s_awqke1[K + 3], s_awqsq1[K + 3], zw1[K + 3];
  real psig_bl, psig_shcu, rmol, rmolh;
  for (int k = kts; k <= kte; ++k) {
    c.edmf_a[k] = zero;
    c.edmf_w[k] = zero;
    c.edmf_qt[k] = zero;
    c.edmf_thl[k] = zero;
    c.edmf_ent[k] = zero;
    c.edmf_qc[k] = zero;
  }
  c.ztop_plume = zero;
  c.maxwidth = zero;
  c.maxmf = zero;
  c.excess_h = zero;
  c.excess_q = zero;
  if (initflag > 0) {  // .not. restart
    const bool initialize_qke = true;
    for (int k = kts; k <= kte; ++k) {
      c.sh[k] = zero;
      c.sm[k] = zero;
      c.el[k] = zero;
      c.tsq[k] = zero;
      c.qsq[k] = zero;
      c.cov[k] = zero;
      c.cldfra_bl[k] = zero;
      c.qc_bl[k] = zero;
      c.qi_bl[k] = zero;
      c.qke[k] = zero;
      c.qke_adv[k] = zero;
      vt1[k] = zero;
      vq1[k] = zero;
      sgm1[k] = zero;
      c.km[k] = zero;
      c.kh[k] = zero;
    }
    if (tke_budget == 1)
      for (int k = kts; k <= kte; ++k) {
        c.qwt[k] = zero;
        c.qshear[k] = zero;
        c.qbuoy[k] = zero;
        c.qdiss[k] = zero;
        c.dqke[k] = zero;
      }
    zw1[kts] = zero;
    for (int k = kts; k <= kte; ++k) {
      sqw1[k] = c.sqv[k] + c.sqc[k] + c.sqi[k];
      thl1[k] = c.th[k] - xlvcp / c.ex[k] * c.sqc[k] - xlscp / c.ex[k] * (c.sqi[k]);
      thlv1[k] = thl1[k] * (one + p608 * c.sqv[k]);
      thv1[k] = c.th[k] * (one + p608 * c.sqv[k] - (c.sqc[k] + c.sqi[k]));
      zw1[k + 1] = zw1[k] + c.dz[k];
    }
    for (int k = kts; k <= kte; ++k)
      c.qke[k] = five * ust * c_max((ust * real(700.) - zw1[k]) / (c_max(ust, p01) * real(700.)), p01);
    get_pblh(kte, c.pblh, thv1, zw1, c.dz, xland, ust, c.kpbl);
    scale_aware(dx, c.pblh, psig_bl, psig_shcu);
    rmol = zero;
    rmolh = zero;
    mym_initialize<K>(kte, ckmod, taue, xland, c.pblh, c.dz, dx, zw1, c.pres, c.ex, c.u, c.v, thl1, c.sqv, c.th, thv1,
                      thlv1, c.sh, c.sm, ust, rmol, c.el, c.qke, c.tsq, c.qsq, c.cov, psig_bl, c.cldfra_bl,
                      op.mixlength, bl_mynn_ess, c.edmf_w, c.edmf_a, initialize_qke);
    if (op.tkeadvect)
      for (int k = kts; k <= kte; ++k) c.qke_adv[k] = c.qke[k];
  }
  if (op.tkeadvect)
    for (int k = kts; k <= kte; ++k) c.qke[k] = c.qke_adv[k];
  if (tke_budget == 1)
    for (int k = kts; k <= kte; ++k) c.dqke[k] = c.qke[k];
  real cldfra_bl1_old[K + 2], qc_bl1_old[K + 2];
  for (int k = kts; k <= kte; ++k) {
    cldfra_bl1_old[k] = c.cldfra_bl[k];
    qc_bl1_old[k] = c.qc_bl[k];
    c.cldfra_bl[k] = zero;
    c.qc_bl[k] = zero;
    c.qi_bl[k] = zero;
    dqs1[k] = zero;
  }
  (void)cldfra_bl1_old; (void)qc_bl1_old;
  for (int k = kts; k <= kte; ++k) {
    c.edmf_a[k] = zero;
    c.edmf_w[k] = zero;
    c.edmf_qc[k] = zero;
    c.edmf_thl[k] = zero;
  }
  for (int k = kts; k <= kte + 1; ++k) {
    s_aw1[k] = zero;
    s_awthl1[k] = zero;
    s_awqt1[k] = zero;
    s_awqv1[k] = zero;
    s_awqc1[k] = zero;
    s_awu1[k] = zero;
    s_awv1[k] = zero;
    s_awqke1[k] = zero;
    s_awqsq1[k] = zero;
  }
  zw1[kts] = zero;
  for (int k = kts; k <= kte; ++k) {
    zw1[k + 1] = zw1[k] + c.dz[k];
    qv1[k] = c.sqv[k] / (one - c.sqv[k]);
    qc1[k] = c.sqc[k] / (one - c.sqv[k]);
    qi1[k] = c.sqi[k] / (one - c.sqv[k]);
    sqw1[k] = c.sqv[k] + c.sqc[k] + c.sqi[k];
    thl1[k] = c.th[k] - xlvcp / c.ex[k] * c.sqc[k] - xlscp / c.ex[k] * (c.sqi[k]);
    thv1[k] = c.th[k] * (one + p608 * c.sqv[k] - (c.sqc[k] + c.sqi[k]));
  }
  (void)qc1; (void)qi1;
  get_pblh(kte, c.pblh, thv1, zw1, c.dz, xland, ust, c.kpbl);
  scale_aware(dx, c.pblh, psig_bl, psig_shcu);
  const real cpm = cp * (one + real(0.84) * c_max(c.sqv[kts], real(1e-8)));
  const real flqv = qfx / c.rho[kts];
  const real flqc = zero;
  const real th_sfc = ts / c.ex[kts];
  const real flq = flqv + flqc;
  const real flt = hfx / (c.rho[kts] * cpm) - xlvcp * flqc / c.ex[kts];
  const real fltv = flt + flqv * p608 * th_sfc;
  rmol = -karman * gtr * fltv / c_max(ust * ust * ust, real(1.0e-6));
  rmolh = -karman * gtr * flt / c_max(ust * ust * ust, real(1.0e-6));
  (void)rmolh;
  real zet = p5 * c.dz[kts] * rmol;
  zet = c_max(zet, real(-10.));
  zet = c_min(zet, real(10.));
  real pmz, phh;
  if ((xland - real(1.5)) >= zero) {
    if (zet >= zero) {
      pmz = one + real(3.0) * zet;
      phh = one + cphh_st * zet;
    } else {
      pmz = one / rpow(one - cphm_unst * zet, real(0.25)) - zet;
      phh = one / std::sqrt(one - cphh_unst * zet);
    }
  } else {
    if (zet >= zero) {
      pmz = real(0.95) + real(3.5) * zet;
      phh = one + cphh_st * zet;
    } else {
      pmz = real(0.95) / rpow(one - cphm_unst * zet, real(0.25)) - zet;
      phh = one / std::sqrt(one - cphh_unst * zet);
    }
  }
  mym_condensation<K>(kte, c.dz, zw1, xland, thl1, sqw1, c.sqc, c.sqi, c.sqs, c.pres, c.ex, c.tsq, c.qsq, c.cov,
                      c.sh, c.el, op.cloudpdf, c.qc_bl, c.qi_bl, c.cldfra_bl, c.pblh, vt1, vq1, c.th, sgm1, closure);
  for (int k = kts; k <= kte; ++k) tkeprod_up[k] = zero;
  MfOut mo{};
  if (bl_mynn_edmf > 0) {
    dmp_mf<K>(kte, delt, zw1, c.dz, c.pres, c.rho, bl_mynn_edmf_mom, bl_mynn_edmf_tke, closure, c.u, c.v, c.w, c.th,
              thl1, thv1, c.tk, sqw1, c.sqv, c.sqc, c.qke, c.qsq, c.ex, vt1, vq1, ust, flt, fltv, flq, c.pblh, c.kpbl,
              dx, xland, th_sfc, c.edmf_a, c.edmf_w, c.edmf_qt, c.edmf_thl, c.edmf_ent, c.edmf_qc, s_aw1, s_awthl1,
              s_awqt1, s_awqv1, s_awqc1, s_awu1, s_awv1, s_awqke1, s_awqsq1, c.qc_bl, c.cldfra_bl, psig_shcu, mo,
              tkeprod_up, c.el);
    c.maxwidth = mo.maxwidth;
    c.maxmf = mo.maxmf;
    c.ztop_plume = mo.ztop;
    c.excess_h = mo.excess_h;
    c.excess_q = mo.excess_q;
  }
  c.maxtkeprod = zero;
  c.cldtop_cooling = zero;
  c.ent_eff = zero;
  c.maxmf_dd = zero;
  c.maxwidth_dd = zero;
  const real delt2 = delt;
  for (int k = kts; k <= kte; ++k) {
    qc_tot1[k] = c_max(c.qc_bl[k], c.sqc[k]);
    qi_tot1[k] = c_max(c.qi_bl[k], c.sqi[k] + c.sqs[k]);
    thlv1[k] = (c.th[k] - xlvcp / c.ex[k] * qc_tot1[k] - xlscp / c.ex[k] * qi_tot1[k]) * (one + p608 * c.sqv[k]);
    thl_tot1[k] = c.th[k] - xlvcp / c.ex[k] * qc_tot1[k] - xlscp / c.ex[k] * qi_tot1[k];
    thv1[k] = c.th[k] * (one + p608 * c.sqv[k] - (qc_tot1[k] + qi_tot1[k]));
  }
  (void)thl_tot1;
  mym_turbulence<K>(kte, xland, closure, ckmod, taue, c.dz, dx, zw1, c.pres, c.ex, c.u, c.v, thl1, thv1, thlv1, sqw1,
                    c.qke, c.tsq, c.qsq, c.cov, vt1, vq1, rmol, flt, fltv, flq, c.pblh, c.th, c.sh, c.sm, c.el, dfm1,
                    dfh1, dfq1, tcd1, qcd1, pdk1, pdt1, pdq1, pdc1, qshear, qbuoy, psig_bl, c.cldfra_bl,
                    op.mixlength, bl_mynn_ess, c.edmf_w, c.edmf_a, tkeprod_up);
  mym_predict<K>(kte, closure, delt2, c.dz, ust, pmz, phh, c.el, dfq1, c.rho, pdk1, pdt1, pdq1, pdc1, c.qke, c.tsq,
                 c.qsq, c.cov, s_aw1, s_awqke1, s_awqsq1, bl_mynn_edmf, bl_mynn_edmf_tke, qwt, qdiss);
  // dissipative heating: gfortran vectorises this loop (glibc libmvec pow and exp, kte-1 trips)
  for (int k = kts; k <= kte - 1; ++k) {
    diss_heat1[k] = c_min(c_max(vm::pow(c.qke[k], real(1.5), k - kts, kte - 1) /
                                    (b1 * c_max(p5 * (c.el[k] + c.el[k + 1]), one)) / cp,
                                zero),
                          real(0.002));
    diss_heat1[k] = diss_heat1[k] * vm::exp(real(-10000.) / c_max(c.pres[k], one), k - kts, kte - 1);
  }
  diss_heat1[kte] = zero;
  mynn_tendencies<K>(kte, delt, c.dz, c.rho, c.u, c.v, c.th, c.tk, qv1, ps, c.pres, c.ex, thl1, c.sqv, c.sqc, c.sqi,
                     sqw1, ust, flt, flq, flqv, flqc, wspd, uoce, voce, tcd1, qcd1, dfm1, dfh1, c.du, c.dv, c.dth,
                     c.dqv, c.dqc, c.dqi, dqs1, diss_heat1, s_aw1, s_awthl1, s_awqt1, s_awqv1, s_awqc1, s_awu1, s_awv1,
                     op.cloudmix, op.mixqt, bl_mynn_edmf, bl_mynn_edmf_mom);
  c.km[kts] = zero;
  c.kh[kts] = zero;
  for (int k = kts + 1; k <= kte; ++k) {
    const real dzk = p5 * (c.dz[k] + c.dz[k - 1]);
    c.km[k] = dfm1[k] * dzk;
    c.kh[k] = dfh1[k] * dzk;
  }
  if (tke_budget == 1) {
    const real phi_m = phim(zet);
    const real ust3 = ust * ust * ust;
    int k = kts;
    c.qshear[k] = real(4.) * (ust3 * phi_m / (karman * c.dz[k])) - c.qshear[k + 1];
    c.qbuoy[k] = real(4.) * (-(ust3 * zet / (karman * c.dz[k]))) - c.qbuoy[k + 1];
    real dummy1[K + 2], dummy2[K + 2];
    for (k = kts; k <= kte - 1; ++k) {
      dummy1[k] = p5 * (c.qshear[k] + c.qshear[k + 1]);
      dummy2[k] = p5 * (c.qbuoy[k] + c.qbuoy[k + 1]);
      c.dqke[k] = p5 * (c.qke[k] - c.dqke[k]) / delt;
    }
    for (k = kts; k <= kte - 1; ++k) {
      c.qshear[k] = dummy1[k];
      c.qbuoy[k] = dummy2[k];
    }
    k = kte;
    c.qshear[k] = zero;
    c.qbuoy[k] = zero;
    c.qwt[k] = zero;
    c.qdiss[k] = zero;
    c.dqke[k] = zero;
  }
  for (int k = kts; k <= kte; ++k) c.qke_adv[k] = c.qke[k];
}

}  // namespace cppwrf::phys::mynnedmf
