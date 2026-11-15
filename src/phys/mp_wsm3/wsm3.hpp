// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// WSM3 (WRF mp_physics = 3): phys/module_mp_wsm3.F of WRF v4.8.0 (wsm3, wsm32D, slope_wsm3, nislfv_rain_plm,
// wsm3init, effectRad_wsm3). Three classes: qv, qci (cloud water above 0 C, cloud ice below; WRF's QCLOUD) and qrs
// (rain above 0 C, snow below; QRAIN). One thread per column, WRF's statements in WRF's order; gfortran vectorises two
// of wsm32D's i loops with libmvec (the ice number loop: log, exp; the saturation loop before the condensation: log,
// exp, exp), vm:: lanes over its..ite = i0..i0+ni-1. wsm3init computes the same constants as wsm5init (Wsm5Params; WSM3
// leaves pacrc unused). The effective radii (has_reqc = has_reqi = has_reqs = 1) are computed at the end of the column
// from WRF's split of qci / qrs by the temperature t of wsm32D.

#include <cmath>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "core/libm.hpp"
#include "phys/mp_wsm5/wsm5.hpp"

namespace cppwrf::phys::wsm3 {

constexpr int kMaxLevels = 128;

template <typename Real, int NL = kMaxLevels>
struct Wsm3Kernel {
  wsm5::Wsm5Params<Real> c;
  Real delt;
  int nz;
  View3<Real> th, q, qci, qrs;
  View3<const Real> w, den, pii, p, delz;
  View2<Real> rain, rainncv, snow, snowncv, sr;
  View3<Real> re_cloud, re_ice, re_snow;  // has_reqc = has_reqi = has_reqs = 1 when re_cloud.p is set
  int i0 = 0, ni = 1;                     // WRF's its..ite (the lanes of the vectorised i loops)

  CPPWRF_HD void operator()(int i, int j) const {
    using std::abs;
    using std::sqrt;
    const int li = i - i0;
    const int kte = nz;
    Real t[NL], qv[NL], qc1[NL], qr1[NL], dn[NL], pr[NL], dz[NL], qsat[NL], rh[NL], cpm[NL], xl[NL], denfac[NL],
        xni[NL], rslope[NL], rslopeb[NL], rslope2[NL], rslope3[NL], work1[NL], work2[NL], fall[NL], denqrs[NL],
        denqci[NL], work1c[NL];
    Real pgen[NL], pisd[NL], paut[NL], pacr[NL], pres[NL], pcon[NL];
    for (int k = 0; k < kte; ++k) {
      t[k] = th(i, j, k) * pii(i, j, k);
      qv[k] = q(i, j, k);
      qc1[k] = qci(i, j, k);
      qr1[k] = qrs(i, j, k);
      dn[k] = den(i, j, k);
      pr[k] = p(i, j, k);
      dz[k] = delz(i, j, k);
    }
    for (int k = 0; k < kte; ++k) {
      qc1[k] = c_max(qc1[k], Real(0.0));
      qr1[k] = c_max(qr1[k], Real(0.0));
    }
    for (int k = 0; k < kte; ++k) {  // cpmcal, xlcal
      cpm[k] = c.cpd * (Real(1.) - c_max(qv[k], c.qmin)) + c_max(qv[k], c.qmin) * c.cpv;
      xl[k] = c.xlv0 - c.xlv1 * (t[k] - c.t0c);
    }
    Real rainncv_v = Real(0.), snowncv_v = Real(0.), sr_v = Real(0.), tstepsnow = Real(0.);
    Real rain_v = rain(i, j), snow_v = snow(i, j);
    const int loops = c_max(int(std::round(delt / Real(c.dtcldcr))), 1);
    Real dtcld = delt / Real(loops);
    if (delt <= Real(c.dtcldcr)) dtcld = delt;
    const Real pi = c.pi, rv = c.rv;
    // statement functions of wsm32D: diffus, viscos, xka, diffac, venfac, conden
    auto diffus = [](Real x, Real y) { return Real(8.794e-5) * libm::exp(libm::log(x) * (Real(1.81))) / y; };
    auto viscos = [](Real x, Real y) { return Real(1.496e-6) * (x * sqrt(x)) / (x + Real(120.)) / y; };
    auto xka = [&](Real x, Real y) { return Real(1.414e3) * viscos(x, y) * y; };
    auto diffac = [&](Real a, Real b, Real cc, Real d, Real e) {
      return d * a * a / (xka(cc, d) * rv * cc * cc) + Real(1.) / (e * diffus(cc, b));
    };
    auto venfac = [&](Real a, Real b, Real cc) {
      return libm::exp(libm::log((viscos(b, cc) / diffus(b, a))) * ((Real(.3333333)))) / sqrt(viscos(b, cc)) *
             sqrt(sqrt(c.den0 / cc));
    };
    auto conden = [&](Real a, Real b, Real cc, Real d, Real e) {
      return (c_max(b, c.qmin) - cc) / (Real(1.) + d * d / (rv * e) * cc / (a * a));
    };

    for (int loop = 1; loop <= loops; ++loop) {
      for (int k = 0; k < kte; ++k) denfac[k] = sqrt(Real(1.0) / dn[k] * c.den0);  // vrec, vsqrt
      const Real cvap = c.cpv, hvap = c.xlv0, hsub = c.xls;
      const Real ttp = c.t0c + Real(0.01);
      const Real dldt = cvap - c.cliq;
      const Real xa = -dldt / rv;
      const Real xb = xa + hvap / (rv * ttp);
      const Real dldti = cvap - c.cice;
      const Real xai = -dldti / rv;
      const Real xbi = xai + hsub / (rv * ttp);
      for (int k = 0; k < kte; ++k) {
        const Real tr = ttp / t[k];
        if (t[k] < ttp)
          qsat[k] = c.psat * (libm::exp(libm::log(tr) * (xai))) * libm::exp(xbi * (Real(1.) - tr));
        else
          qsat[k] = c.psat * (libm::exp(libm::log(tr) * (xa))) * libm::exp(xb * (Real(1.) - tr));
        qsat[k] = c_min(qsat[k], Real(0.99) * pr[k]);
        qsat[k] = c.ep2 * qsat[k] / (pr[k] - qsat[k]);
        qsat[k] = c_max(qsat[k], c.qmin);
        rh[k] = c_max(qv[k] / qsat[k], c.qmin);
      }
      for (int k = 0; k < kte; ++k) {
        pres[k] = paut[k] = pacr[k] = pgen[k] = pisd[k] = pcon[k] = fall[k] = Real(0.);
        xni[k] = Real(1.e3);
      }
      for (int k = 0; k < kte; ++k)  // vectorised (libmvec log, exp)
        xni[k] = c_min(c_max(Real(5.38e7) * vm::exp(vm::log((dn[k] * c_max(qc1[k], c.qmin)), li, ni) * (Real(0.75)), li, ni),
                             Real(1.e3)),
                       Real(1.e6));
      // sedimentation of rain / snow
      slope(rslope, rslopeb, rslope2, rslope3, work1, qr1, dn, denfac, t, kte);
      for (int k = kte - 1; k >= 0; --k) denqrs[k] = dn[k] * qr1[k];
      Real delqrs, delqi;
      wsm5::wsm5_nislfv<Real, NL>(c, kte, dn, denfac, t, dz, work1, denqrs, delqrs, dtcld, 3, 1);
      for (int k = 0; k < kte; ++k) {
        qr1[k] = c_max(denqrs[k] / dn[k], Real(0.));
        fall[k] = denqrs[k] * work1[k] / dz[k];
      }
      fall[0] = delqrs / dz[0] / dtcld;
      // sedimentation of cloud ice
      for (int k = kte - 1; k >= 0; --k) {
        if (t[k] < c.t0c && qc1[k] > Real(0.)) {
          const Real xmi = dn[k] * qc1[k] / xni[k];
          const Real diameter = c_max(Real(c.dicon) * sqrt(xmi), Real(1.e-25));
          work1c[k] = Real(1.49e4) * libm::exp(libm::log(diameter) * (Real(1.31)));
        } else {
          work1c[k] = Real(0.);
        }
      }
      for (int k = kte - 1; k >= 0; --k) denqci[k] = dn[k] * qc1[k];
      wsm5::wsm5_nislfv<Real, NL>(c, kte, dn, denfac, t, dz, work1c, denqci, delqi, dtcld, 1, 0);
      for (int k = 0; k < kte; ++k) qc1[k] = c_max(denqci[k] / dn[k], Real(0.));
      const Real fallc1 = delqi / dz[0] / dtcld;
      // freezing / melting at the freezing level (Fortran levels: mstep, kwork1, kwork2)
      {
        int mstep = 0;
        for (int k = 1; k <= kte; ++k)
          if (t[k - 1] >= c.t0c) mstep = k;
        const int kwork2 = mstep;
        int kwork1 = mstep;
        if (mstep != 0)
          if (w(i, j, mstep - 1) > Real(0.)) kwork1 = mstep + 1;
        int k = kwork1;
        const int kk = kwork2;
        // k = kte + 1 (rising air at a warm model top) reads WRF's unused top level: not reproduced, kept at kte
        if (k > kte) k = kte;
        if (k * kk >= 1) {
          const Real qrsci = qr1[k - 1] + qc1[k - 1];
          if (qrsci > Real(0.) || fall[kk - 1] > Real(0.)) {
            const Real frzmlt = c_min(c_max(-w(i, j, k - 1) * qrsci / dz[k - 1], -qrsci / dtcld), qrsci / dtcld);
            const Real snomlt =
                c_min(c_max(fall[kk - 1] / dn[kk - 1], -qr1[k - 1] / dtcld), qr1[k - 1] / dtcld);
            if (k == kk) {
              t[k - 1] = t[k - 1] - c.xlf0 / cpm[k - 1] * (frzmlt + snomlt) * dtcld;
            } else {
              t[k - 1] = t[k - 1] - c.xlf0 / cpm[k - 1] * frzmlt * dtcld;
              t[kk - 1] = t[kk - 1] - c.xlf0 / cpm[kk - 1] * snomlt * dtcld;
            }
          }
        }
      }
      // the precipitation at the surface
      {
        Real fallsum = fall[0];
        Real fallsum_qsi = Real(0.);
        if ((c.t0c - t[0]) > Real(0)) {
          fallsum = fallsum + fallc1;
          fallsum_qsi = fall[0] + fallc1;
        }
        if (fallsum > Real(0.)) {
          rainncv_v = fallsum * dz[0] / c.denr * dtcld * Real(1000.) + rainncv_v;
          rain_v = fallsum * dz[0] / c.denr * dtcld * Real(1000.) + rain_v;
        }
        if (fallsum_qsi > Real(0.)) {
          tstepsnow = fallsum_qsi * dz[0] / c.denr * dtcld * Real(1000.) + tstepsnow;
          snowncv_v = fallsum_qsi * dz[0] / c.denr * dtcld * Real(1000.) + snowncv_v;
          snow_v = fallsum_qsi * dz[0] / c.denr * dtcld * Real(1000.) + snow_v;
        }
        if (fallsum > Real(0.)) sr_v = snowncv_v / (rainncv_v + Real(1.e-12));
      }
      slope(rslope, rslopeb, rslope2, rslope3, work1, qr1, dn, denfac, t, kte);
      for (int k = 0; k < kte; ++k) {
        if (t[k] >= c.t0c)
          work1[k] = diffac(xl[k], pr[k], t[k], dn[k], qsat[k]);
        else
          work1[k] = diffac(c.xls, pr[k], t[k], dn[k], qsat[k]);
        work2[k] = venfac(pr[k], t[k], dn[k]);
      }
      for (int k = 0; k < kte; ++k) {
        const Real supsat = c_max(qv[k], c.qmin) - qsat[k];
        const Real satdt = supsat / dtcld;
        if (t[k] >= c.t0c) {
          if (qc1[k] > c.qc0) {  // paut: auto conversion of cloud water
            paut[k] = c.qck1 * libm::exp(libm::log(qc1[k]) * ((Real(7. / 3.))));
            paut[k] = c_min(paut[k], qc1[k] / dtcld);
          }
          if (qr1[k] > Real(c.qcrmin) && qc1[k] > c.qmin)  // pacr: accretion of cloud water by rain
            pacr[k] = c_min(c.pacrr * rslope3[k] * rslopeb[k] * qc1[k] * denfac[k], qc1[k] / dtcld);
          if (qr1[k] > Real(0.)) {  // pres: evaporation / condensation of rain
            const Real coeres = rslope2[k] * sqrt(rslope[k] * rslopeb[k]);
            pres[k] = (rh[k] - Real(1.)) * (c.precr1 * rslope2[k] + c.precr2 * work2[k] * coeres) / work1[k];
            if (pres[k] < Real(0.)) {
              pres[k] = c_max(pres[k], -qr1[k] / dtcld);
              pres[k] = c_max(pres[k], satdt / Real(2));
            } else {
              pres[k] = c_min(pres[k], satdt / Real(2));
            }
          }
        } else {
          const Real supcol = c.t0c - t[k];
          const Real n0sfac = c_max(c_min(libm::exp(Real(c.alpha) * supcol), Real(c.n0smax / c.n0s)), Real(1.));
          int ifsat = 0;
          xni[k] = c_min(c_max(Real(5.38e7) * libm::exp(libm::log((dn[k] * c_max(qc1[k], c.qmin))) * (Real(0.75))),
                               Real(1.e3)),
                         Real(1.e6));
          const Real eacrs = libm::exp(Real(0.07) * (-supcol));
          if (qr1[k] > Real(c.qcrmin) && qc1[k] > c.qmin) {  // pacr: accretion of cloud ice by snow
            const Real xmi = dn[k] * qc1[k] / xni[k];
            const Real diameter = c_min(Real(c.dicon) * sqrt(xmi), Real(c.dimax));
            const Real vt2i = Real(1.49e4) * libm::pow(diameter, Real(1.31));
            const Real vt2s = c.pvts * rslopeb[k] * denfac[k];
            const Real acrfac = Real(2.) * rslope3[k] + Real(2.) * diameter * rslope2[k] + diameter * diameter * rslope[k];
            pacr[k] = c_min(pi * qc1[k] * eacrs * Real(c.n0s) * n0sfac * abs(vt2s - vt2i) * acrfac / Real(4.),
                            qc1[k] / dtcld);
          }
          if (qc1[k] > Real(0.)) {  // pisd: deposition / sublimation of cloud ice
            const Real xmi = dn[k] * qc1[k] / xni[k];
            const Real diameter = Real(c.dicon) * sqrt(xmi);
            pisd[k] = Real(4.) * diameter * xni[k] * (rh[k] - Real(1.)) / work1[k];
            if (pisd[k] < Real(0.)) {
              pisd[k] = c_max(pisd[k], satdt / Real(2));
              pisd[k] = c_max(pisd[k], -qc1[k] / dtcld);
            } else {
              pisd[k] = c_min(pisd[k], satdt / Real(2));
            }
            if (abs(pisd[k]) >= abs(satdt)) ifsat = 1;
          }
          if (qr1[k] > Real(0.) && ifsat != 1) {  // pres: deposition / sublimation of snow
            const Real coeres = rslope2[k] * sqrt(rslope[k] * rslopeb[k]);
            pres[k] = (rh[k] - Real(1.)) * n0sfac * (c.precs1 * rslope2[k] + c.precs2 * work2[k] * coeres) / work1[k];
            const Real supice = satdt - pisd[k];
            if (pres[k] < Real(0.)) {
              pres[k] = c_max(pres[k], -qr1[k] / dtcld);
              pres[k] = c_max(c_max(pres[k], satdt / Real(2)), supice);
            } else {
              pres[k] = c_min(c_min(pres[k], satdt / Real(2)), supice);
            }
            if (abs(pisd[k] + pres[k]) >= abs(satdt)) ifsat = 1;
          }
          if (supsat > Real(0) && ifsat != 1) {  // pgen: generation (nucleation) of ice
            const Real supice = satdt - pisd[k] - pres[k];
            const Real xni0 = Real(1.e3) * libm::exp(Real(0.1) * supcol);
            const Real roqi0 = Real(4.92e-11) * libm::exp(libm::log(xni0) * (Real(1.33)));
            pgen[k] = c_max(Real(0.), (roqi0 / dn[k] - c_max(qc1[k], Real(0.))) / dtcld);
            pgen[k] = c_min(c_min(pgen[k], satdt), supice);
          }
          if (qc1[k] > Real(0.)) {  // paut: conversion (aggregation) of ice to snow
            const Real qimax = c.roqimax / dn[k];
            paut[k] = c_max(Real(0.), (qc1[k] - qimax) / dtcld);
          }
        }
      }
      // check mass conservation of the generation terms and feedback to the large scale
      for (int k = 0; k < kte; ++k) {
        const Real qciik = c_max(c.qmin, qc1[k]);
        const Real delqci = (paut[k] + pacr[k] - pgen[k] - pisd[k]) * dtcld;
        if (delqci >= qciik) {
          const Real facqci = qciik / delqci;
          paut[k] = paut[k] * facqci;
          pacr[k] = pacr[k] * facqci;
          pgen[k] = pgen[k] * facqci;
          pisd[k] = pisd[k] * facqci;
        }
        const Real qik = c_max(c.qmin, qv[k]);
        const Real delq = (pres[k] + pgen[k] + pisd[k]) * dtcld;
        if (delq >= qik) {
          const Real facq = qik / delq;
          pres[k] = pres[k] * facq;
          pgen[k] = pgen[k] * facq;
          pisd[k] = pisd[k] * facq;
        }
        work2[k] = -pres[k] - pgen[k] - pisd[k];
        qv[k] = qv[k] + work2[k] * dtcld;
        qc1[k] = c_max(qc1[k] - (paut[k] + pacr[k] - pgen[k] - pisd[k]) * dtcld, Real(0.));
        qr1[k] = c_max(qr1[k] + (paut[k] + pacr[k] + pres[k]) * dtcld, Real(0.));
        if (t[k] < c.t0c)
          t[k] = t[k] - c.xls * work2[k] / cpm[k] * dtcld;
        else
          t[k] = t[k] - xl[k] * work2[k] / cpm[k] * dtcld;
      }
      for (int k = 0; k < kte; ++k) {  // vectorised (libmvec log, exp, exp); WRF's denfac here is not used again
        const Real tr = ttp / t[k];
        qsat[k] = c.psat * (vm::exp(vm::log(tr, li, ni) * (xa), li, ni)) * vm::exp(xb * (Real(1.) - tr), li, ni);
        qsat[k] = c_min(qsat[k], Real(0.99) * pr[k]);
        qsat[k] = c.ep2 * qsat[k] / (pr[k] - qsat[k]);
        qsat[k] = c_max(qsat[k], c.qmin);
      }
      // pcon: condensational / evaporational rate of cloud water
      for (int k = 0; k < kte; ++k) {
        work1[k] = conden(t[k], qv[k], qsat[k], xl[k], cpm[k]);
        work2[k] = qc1[k] + work1[k];
        pcon[k] = c_min(c_max(work1[k], Real(0.)), c_max(qv[k], Real(0.))) / dtcld;
        if (qc1[k] > Real(0.) && work1[k] < Real(0) && t[k] > c.t0c) pcon[k] = c_max(work1[k], -qc1[k]) / dtcld;
        qv[k] = qv[k] - pcon[k] * dtcld;
        qc1[k] = c_max(qc1[k] + pcon[k] * dtcld, Real(0.));
        t[k] = t[k] + pcon[k] * xl[k] / cpm[k] * dtcld;
      }
      for (int k = 0; k < kte; ++k) {
        if (qc1[k] <= c.qmin) qc1[k] = Real(0.0);
        if (qr1[k] <= Real(c.qcrmin)) qr1[k] = Real(0.0);
      }
    }
    for (int k = 0; k < kte; ++k) {
      th(i, j, k) = t[k] / pii(i, j, k);
      q(i, j, k) = qv[k];
      qci(i, j, k) = qc1[k];
      qrs(i, j, k) = qr1[k];
    }
    rain(i, j) = rain_v;
    rainncv(i, j) = rainncv_v;
    snow(i, j) = snow_v;
    snowncv(i, j) = snowncv_v;
    sr(i, j) = sr_v;
    if (re_cloud.p) {  // effectRad_wsm3 (= effectRad_wsm5 / WSM6's) on the split of qci / qrs by t
      constexpr Real R1 = Real(1.E-12), R2 = Real(1.E-6), obmr = Real(1.0) / Real(3.0), nc0 = Real(3.E8);
      for (int k = 0; k < kte; ++k) {
        const Real t1d = th(i, j, k) * pii(i, j, k), rho = dn[k];
        const bool warm = t[k] >= c.t0c;
        const Real qc1d = warm ? qc1[k] : Real(0.0), qi1d = warm ? Real(0.0) : qc1[k],
                   qs1d = warm ? Real(0.0) : qr1[k];
        Real rc = Real(2.49E-6), ri = Real(4.99E-6), rs = Real(9.99E-6);  // RE_QC_BG, RE_QI_BG, RE_QS_BG
        const Real rqc = c_max(R1, qc1d * rho);
        const Real rqi = c_max(R1, qi1d * rho);
        Real temp = (rho * c_max(qi1d, c.qmin));
        temp = sqrt(sqrt(temp * temp * temp));
        const Real nik = c_min(c_max(Real(5.38e7) * temp, Real(1.e3)), Real(1.e6));
        const Real rni = c_max(R2, nik * rho);
        const Real rqs = c_max(R1, qs1d * rho);
        if (rqc > R1) {
          const Real lamdac = libm::pow(c.pidnc * nc0 / rqc, obmr);
          rc = c_max(Real(2.51E-6), c_min(Real(1.5) * (Real(1.0) / lamdac), Real(50.E-6)));
        }
        if (rqi > R1 && rni > R2) {
          const Real diai = Real(11.9) * sqrt(rqi / nik);
          ri = c_max(Real(10.01E-6), c_min(Real(0.75) * Real(0.163) * diai, Real(125.E-6)));
        }
        if (rqs > R1) {
          const Real supcol = c.t0c - t1d;
          const Real n0sfac = c_max(c_min(libm::exp(Real(c.alpha) * supcol), Real(c.n0smax / c.n0s)), Real(1.));
          const Real lamdas = sqrt(sqrt(c.pidn0s * n0sfac / rqs));
          rs = c_max(Real(25.E-6), c_min(Real(0.5) * (Real(1.) / lamdas), Real(999.E-6)));
        }
        re_cloud(i, j, k) = c_max(Real(2.49E-6), c_min(rc, Real(50.E-6)));
        re_ice(i, j, k) = c_max(Real(4.99E-6), c_min(ri, Real(125.E-6)));
        re_snow(i, j, k) = c_max(Real(9.99E-6), c_min(rs, Real(999.E-6)));
      }
    }
  }

  // slope_wsm3: rain (t >= t0c) or snow slope parameters and fall speed vt
  CPPWRF_HD void slope(Real* rs, Real* rb, Real* r2, Real* r3, Real* vt, const Real* qr, const Real* dn,
                       const Real* denfac, const Real* t, int km) const {
    for (int k = 0; k < km; ++k) {
      Real pvt;
      if (t[k] >= Real(273.15)) {
        pvt = c.pvtr;
        if (qr[k] <= Real(c.qcrmin)) {
          rs[k] = c.rslopermax;
          rb[k] = c.rsloperbmax;
          r2[k] = c.rsloper2max;
          r3[k] = c.rsloper3max;
        } else {
          rs[k] = Real(1.) / std::sqrt(std::sqrt(c.pidn0r / (qr[k] * dn[k])));
          rb[k] = libm::exp(libm::log(rs[k]) * (Real(c.bvtr)));
          r2[k] = rs[k] * rs[k];
          r3[k] = r2[k] * rs[k];
        }
      } else {
        const Real supcol = Real(273.15) - t[k];
        const Real n0sfac = c_max(c_min(libm::exp(Real(c.alpha) * supcol), Real(c.n0smax / c.n0s)), Real(1.));
        pvt = c.pvts;
        if (qr[k] <= Real(c.qcrmin)) {
          rs[k] = c.rslopesmax;
          rb[k] = c.rslopesbmax;
          r2[k] = c.rslopes2max;
          r3[k] = c.rslopes3max;
        } else {
          rs[k] = Real(1.) / std::sqrt(std::sqrt(c.pidn0s * n0sfac / (qr[k] * dn[k])));
          rb[k] = libm::exp(libm::log(rs[k]) * (Real(c.bvts)));
          r2[k] = rs[k] * rs[k];
          r3[k] = r2[k] * rs[k];
        }
      }
      vt[k] = pvt * rb[k] * denfac[k];
      if (qr[k] <= Real(0.0)) vt[k] = Real(0.0);
    }
  }
};

}  // namespace cppwrf::phys::wsm3
