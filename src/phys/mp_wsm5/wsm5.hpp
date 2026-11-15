// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// WSM5 (WRF mp_physics = 4): phys/module_mp_wsm5.F of WRF v4.8.0 (wsm5, wsm52D, slope_wsm5, slope_rain,
// slope_snow, nislfv_rain_plm, wsm5init). WRF's wsm52D works on an i-k slab (i fastest); every column is independent,
// so a thread runs one column with the statements in WRF's order. gfortran (-O2 -ftree-vectorize) vectorises three of
// the i loops with libmvec (vm:: lanes over its..ite = i0..i0+ni-1): the saturation loop at the start of each
// microphysics loop (log, exp), the diffusion / ventilation loop (log(t), exp, log, exp) and the saturation loop
// before the condensation (log, exp); every other exp / log / pow is the scalar libm. effectRad_wsm5 is WSM6's
// (phys/mp_wsm6 Wsm6EffRadKernel: the same statements and constants). refl10cm_wsm5 (DOUBLE PRECISION, REAL*16 in
// the -r8 build) is not ported: do_radar_ref = 1 with WSM5 is an error.

#include <cmath>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "core/libm.hpp"
#include "phys/mp_wsm6/wsm6.hpp"

namespace cppwrf::phys::wsm5 {

constexpr int kMaxLevels = 128;

template <typename Real>
struct Wsm5Params {
  // the microphysics_driver arguments (module_model_constants)
  Real cpd{}, cpv{}, rv{}, t0c{}, ep2{}, qmin{}, xls{}, xlv0{}, xlf0{}, den0{}, denr{}, cliq{}, cice{}, psat{};
  // module_mp_wsm5.F parameters
  static constexpr double dtcldcr = 120., n0r = 8.e6, avtr = 841.9, bvtr = 0.8, r0 = .8e-5, peaut = .55,
                          xncr = 3.e8, xmyu = 1.718e-5, avts = 11.72, bvts = .41, n0smax = 1.e11, lamdarmax = 8.e4,
                          lamdasmax = 1.e5, dicon = 11.9, dimax = 500.e-6, n0s = 2.e6, alpha = .12, pfrz1 = 100.,
                          pfrz2 = 0.66, qcrmin = 1.e-9, eacrc = 1.0;
  // wsm5init (module SAVE variables)
  Real qc0{}, qck1{}, pidnc{}, bvtr1{}, bvtr2{}, bvtr3{}, bvtr4{}, g1pbr{}, g3pbr{}, g4pbr{}, g5pbro2{}, pvtr{},
      eacrr{}, pacrr{}, precr1{}, precr2{}, roqimax{}, bvts1{}, bvts2{}, bvts3{}, bvts4{}, g1pbs{}, g3pbs{}, g4pbs{},
      g5pbso2{}, pvts{}, pacrs{}, precs1{}, precs2{}, pidn0r{}, pidn0s{}, xlv1{}, pacrc{}, pi{}, rslopermax{},
      rslopesmax{}, rsloperbmax{}, rslopesbmax{}, rsloper2max{}, rslopes2max{}, rsloper3max{}, rslopes3max{};

  // wsm5init(den0, denr, dens, cl, cpv) on the host. gfortran folds the powers of parameters (r0**3, dimax**8,
  // avtr**.5, avts**.5) and those of the just-assigned rslopermax / rslopesmax (correctly rounded: libm::pow of the
  // constants, folded by gcc as well; checked against WRF's module variables by tests/phys/mp_wsm5)
  void init(Real p_den0, Real p_denr, Real dens, Real cl, Real p_cpv) {
    den0 = p_den0;
    denr = p_denr;
    pi = Real(4. * std::atan(1.));
    xlv1 = cl - p_cpv;
    constexpr double r03 = r0 * r0 * r0;  // r0**3, folded by gfortran: the correctly rounded cube, = r0*r0*r0 here
    qc0 = Real(4. / 3.) * pi * denr * Real(r03) * Real(xncr) / den0;
    qck1 = Real(.104 * 9.8 * peaut) / Real(std::pow(double(Real(xncr) * denr), 1. / 3.)) / Real(xmyu) *
           Real(std::pow(double(den0), 4. / 3.));
    pidnc = pi * denr / Real(6.);
    bvtr1 = Real(1. + bvtr);
    bvtr2 = Real(2.5 + .5 * bvtr);
    bvtr3 = Real(3. + bvtr);
    bvtr4 = Real(4. + bvtr);
    g1pbr = wsm6::rgmma_calc(bvtr1);
    g3pbr = wsm6::rgmma_calc(bvtr3);
    g4pbr = wsm6::rgmma_calc(bvtr4);
    g5pbro2 = wsm6::rgmma_calc(bvtr2);
    pvtr = Real(avtr) * g4pbr / Real(6.);
    eacrr = Real(1.0);
    pacrr = pi * Real(n0r) * Real(avtr) * g3pbr * Real(.25) * eacrr;
    precr1 = Real(2.) * pi * Real(n0r) * Real(.78);
    precr2 = Real(2.) * pi * Real(n0r) * Real(.31) * Real(std::sqrt(avtr)) * g5pbro2;
    roqimax = Real(2.08e22 * std::pow(dimax, 8.));
    bvts1 = Real(1. + bvts);
    bvts2 = Real(2.5 + .5 * bvts);
    bvts3 = Real(3. + bvts);
    bvts4 = Real(4. + bvts);
    g1pbs = wsm6::rgmma_calc(bvts1);
    g3pbs = wsm6::rgmma_calc(bvts3);
    g4pbs = wsm6::rgmma_calc(bvts4);
    g5pbso2 = wsm6::rgmma_calc(bvts2);
    pvts = Real(avts) * g4pbs / Real(6.);
    pacrs = pi * Real(n0s) * Real(avts) * g3pbs * Real(.25);
    precs1 = Real(4. * n0s * .65);
    precs2 = Real(4. * n0s * .44) * Real(std::sqrt(avts)) * g5pbso2;
    pidn0r = pi * denr * Real(n0r);
    pidn0s = pi * dens * Real(n0s);
    pacrc = pi * Real(n0s) * Real(avts) * g3pbs * Real(.25) * Real(eacrc);
    rslopermax = Real(1. / lamdarmax);
    rslopesmax = Real(1. / lamdasmax);
    rsloperbmax = Real(std::pow(1. / lamdarmax, bvtr));
    rslopesbmax = Real(std::pow(1. / lamdasmax, bvts));
    rsloper2max = rslopermax * rslopermax;
    rslopes2max = rslopesmax * rslopesmax;
    rsloper3max = rsloper2max * rslopermax;
    rslopes3max = rslopes2max * rslopesmax;
  }
};

// the microphysics_driver arguments (module_model_constants: CPV = 4 R_V, QMIN = epsilon = 1e-15, PSAT = 610.78)
// and wsm5init(rhoair0, rhowater, rhosnow, cliq, cpv) of phy_init
template <typename Real>
inline Wsm5Params<Real> wsm5_params() {
  namespace pc = phys::cst;
  Wsm5Params<Real> c;
  c.cpd = Real(pc::cp);
  c.cpv = Real(4.0 * pc::r_v);
  c.rv = Real(pc::r_v);
  c.t0c = Real(pc::svpt0);
  c.ep2 = Real(pc::ep_2);
  c.qmin = Real(1.e-15);
  c.xls = Real(pc::xls);
  c.xlv0 = Real(pc::xlv);
  c.xlf0 = Real(pc::xlf);
  c.cliq = Real(pc::cliq);
  c.cice = Real(pc::cice);
  c.psat = Real(610.78);
  c.init(Real(pc::rhoair0), Real(pc::rhowater), Real(pc::rhosnow), Real(pc::cliq), Real(4.0 * pc::r_v));
  return c;
}

// nislfv_rain_plm for one column (WRF's i_loop body; km levels): semi-Lagrangian sedimentation of rql = den*q with the
// fall speed ww; id 1: rain (slope_rain), 2: snow (slope_snow), 3: WSM3's rain or snow (slope_wsm3) in the iter
// fall-speed iterations; precip: the amount that leaves the column bottom
template <typename Real, int NL>
CPPWRF_HD inline void wsm5_nislfv(const Wsm5Params<Real>& c, int km, const Real* den, const Real* denfac,
                                  const Real* tk, const Real* dz, const Real* wwl, Real* rql, Real& precip, Real dt,
                                  int id, int iter) {
  precip = Real(0.0);
  Real allold = Real(0.0);
  for (int k = 0; k < km; ++k) allold = allold + rql[k];
  if (allold <= Real(0.0)) return;
  Real ww[NL], qq[NL], wd[NL], wa[NL], was[NL], qn[NL], qr[NL];
  Real wi[NL + 1], zi[NL + 1], za[NL + 1], dza[NL + 1], qa[NL + 1], qmi[NL + 1], qpi[NL + 1];
  for (int k = 0; k < km; ++k) {
    qq[k] = rql[k];
    ww[k] = wwl[k];
  }
  zi[0] = Real(0.0);
  for (int k = 0; k < km; ++k) zi[k + 1] = zi[k] + dz[k];
  for (int k = 0; k < km; ++k) wd[k] = ww[k];
  int n = 1;
  for (;;) {  // label 100
    // (WRF's first wi loop is overwritten by the next statements)
    const Real fa1 = Real(9.) / Real(16.), fa2 = Real(1.) / Real(16.);
    wi[0] = ww[0];
    wi[1] = Real(0.5) * (ww[1] + ww[0]);
    for (int k = 2; k < km - 1; ++k) wi[k] = fa1 * (ww[k] + ww[k - 1]) - fa2 * (ww[k + 1] + ww[k - 2]);
    wi[km - 1] = Real(0.5) * (ww[km - 1] + ww[km - 2]);
    wi[km] = ww[km - 1];
    for (int k = 1; k < km; ++k)
      if (ww[k] == Real(0.0)) wi[k] = ww[k - 1];
    const Real con1 = Real(0.05);
    for (int k = km - 1; k >= 0; --k) {
      const Real decfl = (wi[k + 1] - wi[k]) * dt / dz[k];
      if (decfl > con1) wi[k] = wi[k + 1] - con1 * dz[k] / dt;
    }
    for (int k = 0; k <= km; ++k) za[k] = zi[k] - wi[k] * dt;
    for (int k = 0; k < km; ++k) dza[k] = za[k + 1] - za[k];
    dza[km] = zi[km] - za[km];
    for (int k = 0; k < km; ++k) {
      qa[k] = qq[k] * dz[k] / dza[k];
      qr[k] = qa[k] / den[k];
    }
    qa[km] = Real(0.0);
    if (n <= iter) {
      for (int k = 0; k < km; ++k) {
        Real rslopeb;
        if (id == 1) {  // slope_rain
          if (qr[k] <= Real(c.qcrmin)) {
            rslopeb = c.rsloperbmax;
          } else {
            const Real rslope = Real(1.) / std::sqrt(std::sqrt(c.pidn0r / (qr[k] * den[k])));
            rslopeb = libm::pow(rslope, Real(c.bvtr));
          }
          wa[k] = c.pvtr * rslopeb * denfac[k];
        } else if (id == 3) {  // slope_wsm3 (module_mp_wsm3.F): rain at t >= t0c, snow below
          if (tk[k] >= Real(273.15)) {
            if (qr[k] <= Real(c.qcrmin)) {
              rslopeb = c.rsloperbmax;
            } else {
              const Real rslope = Real(1.) / std::sqrt(std::sqrt(c.pidn0r / (qr[k] * den[k])));
              rslopeb = libm::exp(libm::log(rslope) * Real(c.bvtr));
            }
            wa[k] = c.pvtr * rslopeb * denfac[k];
          } else {
            const Real supcol = Real(273.15) - tk[k];
            const Real n0sfac = c_max(c_min(libm::exp(Real(c.alpha) * supcol), Real(c.n0smax / c.n0s)), Real(1.));
            if (qr[k] <= Real(c.qcrmin)) {
              rslopeb = c.rslopesbmax;
            } else {
              const Real rslope = Real(1.) / std::sqrt(std::sqrt(c.pidn0s * n0sfac / (qr[k] * den[k])));
              rslopeb = libm::exp(libm::log(rslope) * Real(c.bvts));
            }
            wa[k] = c.pvts * rslopeb * denfac[k];
          }
        } else {  // slope_snow
          const Real supcol = Real(273.15) - tk[k];
          const Real n0sfac = c_max(c_min(libm::exp(Real(c.alpha) * supcol), Real(c.n0smax / c.n0s)), Real(1.));
          if (qr[k] <= Real(c.qcrmin)) {
            rslopeb = c.rslopesbmax;
          } else {
            const Real rslope = Real(1.) / std::sqrt(std::sqrt(c.pidn0s * n0sfac / (qr[k] * den[k])));
            rslopeb = libm::pow(rslope, Real(c.bvts));
          }
          wa[k] = c.pvts * rslopeb * denfac[k];
        }
        if (qr[k] <= Real(0.0)) wa[k] = Real(0.0);
      }
      if (n >= 2)
        for (int k = 0; k < km; ++k) wa[k] = Real(0.5) * (wa[k] + was[k]);
      for (int k = 0; k < km; ++k) ww[k] = Real(0.5) * (wd[k] + wa[k]);
      for (int k = 0; k < km; ++k) was[k] = wa[k];
      n = n + 1;
      continue;
    }
    break;
  }
  // piecewise linear profile in each layer, monotone
  for (int k = 1; k < km; ++k) {
    const Real dip = (qa[k + 1] - qa[k]) / (dza[k + 1] + dza[k]);
    const Real dim = (qa[k] - qa[k - 1]) / (dza[k - 1] + dza[k]);
    if (dip * dim <= Real(0.0)) {
      qmi[k] = qa[k];
      qpi[k] = qa[k];
    } else {
      qpi[k] = qa[k] + Real(0.5) * (dip + dim) * dza[k];
      qmi[k] = Real(2.0) * qa[k] - qpi[k];
      if (qpi[k] < Real(0.0) || qmi[k] < Real(0.0)) {
        qpi[k] = qa[k];
        qmi[k] = qa[k];
      }
    }
  }
  qpi[0] = qa[0];
  qmi[0] = qa[0];
  qmi[km] = qa[km];
  qpi[km] = qa[km];
  // interpolation to the regular grid (Fortran indices kb, kt, kk: 1-based)
  for (int k = 0; k < km; ++k) qn[k] = Real(0.0);
  int kb = 1, kt = 1;
  for (int k = 1; k <= km; ++k) {
    kb = kb - 1 > 1 ? kb - 1 : 1;
    kt = kt - 1 > 1 ? kt - 1 : 1;
    if (zi[k - 1] >= za[km]) break;
    for (int kk = kb; kk <= km; ++kk)
      if (zi[k - 1] <= za[kk]) {
        kb = kk;
        break;
      }
    for (int kk = kt; kk <= km; ++kk)
      if (zi[k] <= za[kk - 1]) {
        kt = kk;
        break;
      }
    kt = kt - 1;
    if (kt == kb) {
      const Real tl = (zi[k - 1] - za[kb - 1]) / dza[kb - 1];
      const Real th = (zi[k] - za[kb - 1]) / dza[kb - 1];
      const Real tl2 = tl * tl;
      const Real th2 = th * th;
      const Real qqd = Real(0.5) * (qpi[kb - 1] - qmi[kb - 1]);
      const Real qqh = qqd * th2 + qmi[kb - 1] * th;
      const Real qql = qqd * tl2 + qmi[kb - 1] * tl;
      qn[k - 1] = (qqh - qql) / (th - tl);
    } else if (kt > kb) {
      const Real tl = (zi[k - 1] - za[kb - 1]) / dza[kb - 1];
      const Real tl2 = tl * tl;
      Real qqd = Real(0.5) * (qpi[kb - 1] - qmi[kb - 1]);
      const Real qql = qqd * tl2 + qmi[kb - 1] * tl;
      const Real dql = qa[kb - 1] - qql;
      Real zsum = (Real(1.) - tl) * dza[kb - 1];
      Real qsum = dql * dza[kb - 1];
      if (kt - kb > 1)
        for (int m = kb + 1; m <= kt - 1; ++m) {
          zsum = zsum + dza[m - 1];
          qsum = qsum + qa[m - 1] * dza[m - 1];
        }
      const Real th = (zi[k] - za[kt - 1]) / dza[kt - 1];
      const Real th2 = th * th;
      qqd = Real(0.5) * (qpi[kt - 1] - qmi[kt - 1]);
      const Real dqh = qqd * th2 + qmi[kt - 1] * th;
      zsum = zsum + th * dza[kt - 1];
      qsum = qsum + dqh * dza[kt - 1];
      qn[k - 1] = qsum / zsum;
    }
  }
  // the precipitation that leaves the bottom
  for (int k = 0; k < km; ++k) {
    if (za[k] < Real(0.0) && za[k + 1] < Real(0.0)) {
      precip = precip + qa[k] * dza[k];
      continue;
    } else if (za[k] < Real(0.0) && za[k + 1] >= Real(0.0)) {
      precip = precip + qa[k] * (Real(0.0) - za[k]);
      break;
    }
    break;
  }
  for (int k = 0; k < km; ++k) rql[k] = qn[k];
}

// wsm5 + wsm52D for the column (i, j); level k = Fortran k - 1. qc/qi = WRF's qci(:,:,1/2), qr/qs = qrs(:,:,1/2)
template <typename Real, int NL = kMaxLevels>
struct Wsm5Kernel {
  Wsm5Params<Real> c;
  Real delt;
  int nz;
  View3<Real> th, q, qc, qr, qi, qs;
  View3<const Real> den, pii, p, delz;
  View2<Real> rain, rainncv, snow, snowncv, sr;
  int i0 = 0, ni = 1;  // WRF's its..ite (the lanes of the vectorised i loops)

  CPPWRF_HD void operator()(int i, int j) const {
    using std::abs;
    using std::sqrt;
    const int li = i - i0;
    const int kte = nz;
    Real t[NL], qv[NL], qci1[NL], qci2[NL], qrs1[NL], qrs2[NL], dn[NL], pr[NL], dz[NL];
    Real qs1[NL], qs2[NL], rh1[NL], rh2[NL], cpm[NL], xl[NL], denfac[NL], xni[NL], n0sfac[NL];
    Real rslope1[NL], rslope2[NL], rslopeb1[NL], rslopeb2[NL], rslope21[NL], rslope22[NL], rslope31[NL],
        rslope32[NL], work11[NL], work12[NL], work2[NL], workr[NL], works[NL], denqrs1[NL], denqrs2[NL],
        denqci[NL], work1c[NL];
    Real prevp[NL], psdep[NL], praut[NL], psaut[NL], pracw[NL], psaci[NL], psacw[NL], pigen[NL], pidep[NL],
        pcond[NL], psmlt[NL], psevp[NL];
    for (int k = 0; k < kte; ++k) {
      t[k] = th(i, j, k) * pii(i, j, k);
      qv[k] = q(i, j, k);
      qci1[k] = qc(i, j, k);
      qci2[k] = qi(i, j, k);
      qrs1[k] = qr(i, j, k);
      qrs2[k] = qs(i, j, k);
      dn[k] = den(i, j, k);
      pr[k] = p(i, j, k);
      dz[k] = delz(i, j, k);
    }
    for (int k = 0; k < kte; ++k) {
      qci1[k] = c_max(qci1[k], Real(0.0));
      qrs1[k] = c_max(qrs1[k], Real(0.0));
      qci2[k] = c_max(qci2[k], Real(0.0));
      qrs2[k] = c_max(qrs2[k], Real(0.0));
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
    const Real pi = c.pi;

    for (int loop = 1; loop <= loops; ++loop) {
      for (int k = 0; k < kte; ++k) denfac[k] = sqrt(Real(1.0) / dn[k] * c.den0);  // vrec, vsqrt
      const Real hsub = c.xls, hvap = c.xlv0, cvap = c.cpv;
      const Real ttp = c.t0c + Real(0.01);
      const Real dldt = cvap - c.cliq;
      const Real xa = -dldt / c.rv;
      const Real xb = xa + hvap / (c.rv * ttp);
      const Real dldti = cvap - c.cice;
      const Real xai = -dldti / c.rv;
      const Real xbi = xai + hsub / (c.rv * ttp);
      for (int k = 0; k < kte; ++k) {  // vectorised (libmvec log, exp)
        const Real xal = t[k] < ttp ? xai : xa;
        const Real xbl = t[k] < ttp ? xbi : xb;
        const Real tr = ttp / t[k];
        const Real logtr = vm::log(tr, li, ni);
        qs1[k] = c.psat * vm::exp(logtr * (xa) + xb * (Real(1.) - tr), li, ni);
        qs1[k] = c_min(qs1[k], Real(0.99) * pr[k]);
        qs1[k] = c.ep2 * qs1[k] / (pr[k] - qs1[k]);
        qs1[k] = c_max(qs1[k], c.qmin);
        rh1[k] = c_max(qv[k] / qs1[k], c.qmin);
        qs2[k] = c.psat * vm::exp(logtr * (xal) + xbl * (Real(1.) - tr), li, ni);
        qs2[k] = c_min(qs2[k], Real(0.99) * pr[k]);
        qs2[k] = c.ep2 * qs2[k] / (pr[k] - qs2[k]);
        qs2[k] = c_max(qs2[k], c.qmin);
        rh2[k] = c_max(qv[k] / qs2[k], c.qmin);
      }
      for (int k = 0; k < kte; ++k) {
        prevp[k] = psdep[k] = praut[k] = psaut[k] = pracw[k] = psaci[k] = psacw[k] = pigen[k] = pidep[k] =
            pcond[k] = psmlt[k] = psevp[k] = Real(0.);
        xni[k] = Real(1.e3);
      }
      for (int k = 0; k < kte; ++k) {
        Real temp = (dn[k] * c_max(qci2[k], c.qmin));
        temp = sqrt(sqrt(temp * temp * temp));
        xni[k] = c_min(c_max(Real(5.38e7) * temp, Real(1.e3)), Real(1.e6));
      }
      slope(rslope1, rslope2, rslopeb1, rslopeb2, rslope21, rslope22, rslope31, rslope32, work11, work12, qrs1, qrs2,
            dn, denfac, t, kte);
      for (int k = kte - 1; k >= 0; --k) {
        workr[k] = work11[k];
        works[k] = work12[k];
        denqrs1[k] = dn[k] * qrs1[k];
        denqrs2[k] = dn[k] * qrs2[k];
        if (qrs1[k] <= Real(0.0)) workr[k] = Real(0.0);
        if (qrs2[k] <= Real(0.0)) works[k] = Real(0.0);
      }
      Real delqrs1, delqrs2, delqi;
      wsm5_nislfv<Real, NL>(c, kte, dn, denfac, t, dz, workr, denqrs1, delqrs1, dtcld, 1, 1);
      wsm5_nislfv<Real, NL>(c, kte, dn, denfac, t, dz, works, denqrs2, delqrs2, dtcld, 2, 1);
      for (int k = 0; k < kte; ++k) {
        qrs1[k] = c_max(denqrs1[k] / dn[k], Real(0.));
        qrs2[k] = c_max(denqrs2[k] / dn[k], Real(0.));
      }
      const Real fall1 = delqrs1 / dz[0] / dtcld;
      const Real fall2 = delqrs2 / dz[0] / dtcld;
      slope(rslope1, rslope2, rslopeb1, rslopeb2, rslope21, rslope22, rslope31, rslope32, work11, work12, qrs1, qrs2,
            dn, denfac, t, kte);
      // melting of snow
      for (int k = kte - 1; k >= 0; --k) {
        const Real supcol = c.t0c - t[k];
        n0sfac[k] = c_max(c_min(libm::exp(Real(c.alpha) * supcol), Real(c.n0smax / c.n0s)), Real(1.));
        if (t[k] > c.t0c && qrs2[k] > Real(0.)) {
          const Real xlf = c.xlf0;
          const Real tk = t[k], tsq = tk * sqrt(tk);
          const Real visc = Real(1.496e-6) * tsq / (tk + Real(120.)) / dn[k];
          work2[k] = libm::exp(libm::log(visc / (Real(8.794e-5) * libm::exp(libm::log(tk) * Real(1.81)) / pr[k])) *
                               Real(.3333333)) /
                     sqrt(visc) * sqrt(sqrt(c.den0 / dn[k]));
          const Real coeres = rslope22[k] * sqrt(rslope2[k] * rslopeb2[k]);
          psmlt[k] = Real(1.414e3) * visc * dn[k] / xlf * (c.t0c - tk) * pi / Real(2.) * n0sfac[k] *
                     (c.precs1 * rslope22[k] + c.precs2 * work2[k] * coeres) / dn[k];
          psmlt[k] = c_min(c_max(psmlt[k] * dtcld / Real(1), -qrs2[k] / Real(1)), Real(0.));  // mstep = 1
          qrs2[k] = qrs2[k] + psmlt[k];
          qrs1[k] = qrs1[k] - psmlt[k];
          t[k] = t[k] + xlf / cpm[k] * psmlt[k];
        }
      }
      // sedimentation of cloud ice
      for (int k = kte - 1; k >= 0; --k) {
        if (qci2[k] <= Real(0.)) {
          work1c[k] = Real(0.);
        } else {
          const Real xmi = dn[k] * qci2[k] / xni[k];
          const Real diameter = c_max(c_min(Real(c.dicon) * sqrt(xmi), Real(c.dimax)), Real(1.e-25));
          work1c[k] = Real(1.49e4) * libm::exp(libm::log(diameter) * Real(1.31));
        }
      }
      for (int k = kte - 1; k >= 0; --k) denqci[k] = dn[k] * qci2[k];
      wsm5_nislfv<Real, NL>(c, kte, dn, denfac, t, dz, work1c, denqci, delqi, dtcld, 1, 0);
      for (int k = 0; k < kte; ++k) qci2[k] = c_max(denqci[k] / dn[k], Real(0.));
      const Real fallc1 = delqi / dz[0] / dtcld;
      // the precipitation at the surface
      {
        const Real fallsum = fall1 + fall2 + fallc1;
        const Real fallsum_qsi = fall2 + fallc1;
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
      // instantaneous melting / freezing of cloud ice, homogeneous and heterogeneous freezing
      for (int k = 0; k < kte; ++k) {
        const Real supcol = c.t0c - t[k];
        Real xlf = c.xls - xl[k];
        if (supcol < Real(0.)) xlf = c.xlf0;
        if (supcol < Real(0.) && qci2[k] > Real(0.)) {
          qci1[k] = qci1[k] + qci2[k];
          t[k] = t[k] - xlf / cpm[k] * qci2[k];
          qci2[k] = Real(0.);
        }
        if (supcol > Real(40.) && qci1[k] > Real(0.)) {
          qci2[k] = qci2[k] + qci1[k];
          t[k] = t[k] + xlf / cpm[k] * qci1[k];
          qci1[k] = Real(0.);
        }
        if (supcol > Real(0.) && qci1[k] > Real(0.)) {
          const Real supcolt = c_min(supcol, Real(50.));
          const Real pfrzdtc = c_min(Real(c.pfrz1) * (libm::exp(Real(c.pfrz2) * supcolt) - Real(1.)) * dn[k] / c.denr /
                                         Real(c.xncr) * qci1[k] * qci1[k] * dtcld,
                                     qci1[k]);
          qci2[k] = qci2[k] + pfrzdtc;
          t[k] = t[k] + xlf / cpm[k] * pfrzdtc;
          qci1[k] = qci1[k] - pfrzdtc;
        }
        if (supcol > Real(0.) && qrs1[k] > Real(0.)) {
          const Real supcolt = c_min(supcol, Real(50.));
          Real temp = rslope1[k];
          temp = temp * temp * temp * temp * temp * temp * temp;
          const Real pfrzdtr = c_min(Real(20.) * (pi * pi) * Real(c.pfrz1) * Real(c.n0r) * c.denr / dn[k] *
                                         (libm::exp(Real(c.pfrz2) * supcolt) - Real(1.)) * temp * dtcld,
                                     qrs1[k]);
          qrs2[k] = qrs2[k] + pfrzdtr;
          t[k] = t[k] + xlf / cpm[k] * pfrzdtr;
          qrs1[k] = qrs1[k] - pfrzdtr;
        }
      }
      slope(rslope1, rslope2, rslopeb1, rslopeb2, rslope21, rslope22, rslope31, rslope32, work11, work12, qrs1, qrs2,
            dn, denfac, t, kte);
      for (int k = 0; k < kte; ++k) {  // vectorised (libmvec log, exp): diffac, venfac
        const Real tk = t[k], tsq = tk * sqrt(tk);
        const Real dfs = Real(8.794e-5) * vm::exp(vm::log(tk, li, ni) * Real(1.81), li, ni);
        work11[k] = dn[k] * xl[k] * xl[k] * (tk + Real(120.)) * dn[k] /
                        (Real(1.414e3) * (Real(1.496e-6) * tsq) * dn[k] * (c.rv * tk * tk)) +
                    pr[k] / (qs1[k] * dfs);
        work12[k] = dn[k] * c.xls * c.xls * (tk + Real(120.)) * dn[k] /
                        (Real(1.414e3) * (Real(1.496e-6) * tsq) * dn[k] * (c.rv * tk * tk)) +
                    pr[k] / (qs2[k] * dfs);
        work2[k] = vm::exp(Real(.3333333) * vm::log(Real(1.496e-6) * tsq * pr[k] / ((tk + Real(120.)) * dn[k] * dfs), li, ni),
                           li, ni) *
                   sqrt(sqrt(c.den0 / dn[k])) / sqrt(Real(1.496e-6) * tsq / ((tk + Real(120.)) * dn[k]));
      }
      // warm rain processes
      for (int k = 0; k < kte; ++k) {
        const Real supsat = c_max(qv[k], c.qmin) - qs1[k];
        const Real satdt = supsat / dtcld;
        if (qci1[k] > c.qc0) {  // praut: auto conversion
          praut[k] = c.qck1 * libm::exp(libm::log(qci1[k]) * Real(7. / 3.));
          praut[k] = c_min(praut[k], qci1[k] / dtcld);
        }
        if (qrs1[k] > Real(c.qcrmin) && qci1[k] > c.qmin)  // pracw: accretion of cloud water by rain
          pracw[k] = c_min(c.pacrr * rslope31[k] * rslopeb1[k] * qci1[k] * denfac[k], qci1[k] / dtcld);
        if (qrs1[k] > Real(0.)) {  // prevp: evaporation / condensation of rain
          const Real coeres = rslope21[k] * sqrt(rslope1[k] * rslopeb1[k]);
          prevp[k] = (rh1[k] - Real(1.)) * (c.precr1 * rslope21[k] + c.precr2 * work2[k] * coeres) / work11[k];
          if (prevp[k] < Real(0.)) {
            prevp[k] = c_max(prevp[k], -qrs1[k] / dtcld);
            prevp[k] = c_max(prevp[k], satdt / Real(2));
          } else {
            prevp[k] = c_min(prevp[k], satdt / Real(2));
          }
        }
      }
      // cold rain processes
      const Real rdtcld = Real(1.) / dtcld;
      for (int k = 0; k < kte; ++k) {
        const Real supcol = c.t0c - t[k];
        n0sfac[k] = c_max(c_min(libm::exp(Real(c.alpha) * supcol), Real(c.n0smax / c.n0s)), Real(1.));
        const Real supsat = c_max(qv[k], c.qmin) - qs2[k];
        const Real satdt = supsat / dtcld;
        int ifsat = 0;
        Real temp = (dn[k] * c_max(qci2[k], c.qmin));
        temp = sqrt(sqrt(temp * temp * temp));
        xni[k] = c_min(c_max(Real(5.38e7) * temp, Real(1.e3)), Real(1.e6));
        const Real eacrs = libm::exp(Real(0.07) * (-supcol));
        if (supcol > Real(0)) {  // psaci: accretion of cloud ice by snow
          if (qrs2[k] > Real(c.qcrmin) && qci2[k] > c.qmin) {
            const Real xmi = dn[k] * qci2[k] / xni[k];
            const Real diameter = c_min(Real(c.dicon) * sqrt(xmi), Real(c.dimax));
            const Real vt2i = Real(1.49e4) * libm::pow(diameter, Real(1.31));
            const Real vt2s = c.pvts * rslopeb2[k] * denfac[k];
            const Real acrfac = Real(2.) * rslope32[k] + Real(2.) * diameter * rslope22[k] + diameter * diameter * rslope2[k];
            psaci[k] = pi * qci2[k] * eacrs * Real(c.n0s) * n0sfac[k] * abs(vt2s - vt2i) * acrfac / Real(4.);
          }
        }
        if (qrs2[k] > Real(c.qcrmin) && qci1[k] > c.qmin)  // psacw: accretion of cloud water by snow
          psacw[k] = c_min(c.pacrc * n0sfac[k] * rslope32[k] * rslopeb2[k] * qci1[k] * denfac[k], qci1[k] * rdtcld);
        if (supcol > Real(0)) {
          if (qci2[k] > Real(0) && ifsat != 1) {  // pidep: deposition / sublimation of cloud ice
            const Real xmi = dn[k] * qci2[k] / xni[k];
            const Real diameter = Real(c.dicon) * sqrt(xmi);
            pidep[k] = Real(4.) * diameter * xni[k] * (rh2[k] - Real(1.)) / work12[k];
            const Real supice = satdt - prevp[k];
            if (pidep[k] < Real(0.)) {
              pidep[k] = c_max(c_max(pidep[k], satdt * Real(.5)), supice);
              pidep[k] = c_max(pidep[k], -qci2[k] * rdtcld);
            } else {
              pidep[k] = c_min(c_min(pidep[k], satdt * Real(.5)), supice);
            }
            if (abs(prevp[k] + pidep[k]) >= abs(satdt)) ifsat = 1;
          }
          if (qrs2[k] > Real(0.) && ifsat != 1) {  // psdep: deposition / sublimation of snow
            const Real coeres = rslope22[k] * sqrt(rslope2[k] * rslopeb2[k]);
            psdep[k] = (rh2[k] - Real(1.)) * n0sfac[k] * (c.precs1 * rslope22[k] + c.precs2 * work2[k] * coeres) /
                       work12[k];
            const Real supice = satdt - prevp[k] - pidep[k];
            if (psdep[k] < Real(0.)) {
              psdep[k] = c_max(psdep[k], -qrs2[k] * rdtcld);
              psdep[k] = c_max(c_max(psdep[k], satdt * Real(.5)), supice);
            } else {
              psdep[k] = c_min(c_min(psdep[k], satdt * Real(.5)), supice);
            }
            if (abs(prevp[k] + pidep[k] + psdep[k]) >= abs(satdt)) ifsat = 1;
          }
          if (supsat > Real(0) && ifsat != 1) {  // pigen: generation (nucleation) of ice
            const Real supice = satdt - prevp[k] - pidep[k] - psdep[k];
            const Real xni0 = Real(1.e3) * libm::exp(Real(0.1) * supcol);
            const Real roqi0 = Real(4.92e-11) * libm::exp(libm::log(xni0) * Real(1.33));
            pigen[k] = c_max(Real(0.), (roqi0 / dn[k] - c_max(qci2[k], Real(0.))) * rdtcld);
            pigen[k] = c_min(c_min(pigen[k], satdt), supice);
          }
          if (qci2[k] > Real(0.)) {  // psaut: conversion (aggregation) of ice to snow
            const Real qimax = c.roqimax / dn[k];
            psaut[k] = c_max(Real(0.), (qci2[k] - qimax) * rdtcld);
          }
        }
        if (supcol < Real(0.)) {  // psevp: evaporation of melting snow
          if (qrs2[k] > Real(0.) && rh1[k] < Real(1.)) psevp[k] = psdep[k] * work12[k] / work11[k];
          psevp[k] = c_min(c_max(psevp[k], -qrs2[k] * rdtcld), Real(0.));
        }
      }
      // check mass conservation of the generation terms and feedback to the large scale
      for (int k = 0; k < kte; ++k) {
        if (t[k] <= c.t0c) {
          Real value = c_max(c.qmin, qci1[k]);
          Real source = (praut[k] + pracw[k] + psacw[k]) * dtcld;
          if (source > value) {
            const Real factor = value / source;
            praut[k] = praut[k] * factor;
            pracw[k] = pracw[k] * factor;
            psacw[k] = psacw[k] * factor;
          }
          value = c_max(c.qmin, qci2[k]);
          source = (psaut[k] + psaci[k] - pigen[k] - pidep[k]) * dtcld;
          if (source > value) {
            const Real factor = value / source;
            psaut[k] = psaut[k] * factor;
            psaci[k] = psaci[k] * factor;
            pigen[k] = pigen[k] * factor;
            pidep[k] = pidep[k] * factor;
          }
          value = c_max(c.qmin, qrs1[k]);
          source = (-praut[k] - pracw[k] - prevp[k]) * dtcld;
          if (source > value) {
            const Real factor = value / source;
            praut[k] = praut[k] * factor;
            pracw[k] = pracw[k] * factor;
            prevp[k] = prevp[k] * factor;
          }
          value = c_max(c.qmin, qrs2[k]);
          source = (-psdep[k] - psaut[k] - psaci[k] - psacw[k]) * dtcld;
          if (source > value) {
            const Real factor = value / source;
            psdep[k] = psdep[k] * factor;
            psaut[k] = psaut[k] * factor;
            psaci[k] = psaci[k] * factor;
            psacw[k] = psacw[k] * factor;
          }
          work2[k] = -(prevp[k] + psdep[k] + pigen[k] + pidep[k]);
          qv[k] = qv[k] + work2[k] * dtcld;
          qci1[k] = c_max(qci1[k] - (praut[k] + pracw[k] + psacw[k]) * dtcld, Real(0.));
          qrs1[k] = c_max(qrs1[k] + (praut[k] + pracw[k] + prevp[k]) * dtcld, Real(0.));
          qci2[k] = c_max(qci2[k] - (psaut[k] + psaci[k] - pigen[k] - pidep[k]) * dtcld, Real(0.));
          qrs2[k] = c_max(qrs2[k] + (psdep[k] + psaut[k] + psaci[k] + psacw[k]) * dtcld, Real(0.));
          const Real xlf = c.xls - xl[k];
          const Real xlwork2 = -c.xls * (psdep[k] + pidep[k] + pigen[k]) - xl[k] * prevp[k] - xlf * psacw[k];
          t[k] = t[k] - xlwork2 / cpm[k] * dtcld;
        } else {
          Real value = c_max(c.qmin, qci1[k]);
          Real source = (praut[k] + pracw[k] + psacw[k]) * dtcld;
          if (source > value) {
            const Real factor = value / source;
            praut[k] = praut[k] * factor;
            pracw[k] = pracw[k] * factor;
            psacw[k] = psacw[k] * factor;
          }
          value = c_max(c.qmin, qrs1[k]);
          source = (-praut[k] - pracw[k] - prevp[k] - psacw[k]) * dtcld;
          if (source > value) {
            const Real factor = value / source;
            praut[k] = praut[k] * factor;
            pracw[k] = pracw[k] * factor;
            prevp[k] = prevp[k] * factor;
            psacw[k] = psacw[k] * factor;
          }
          value = c_max(Real(c.qcrmin), qrs2[k]);
          source = (-psevp[k]) * dtcld;
          if (source > value) {
            const Real factor = value / source;
            psevp[k] = psevp[k] * factor;
          }
          work2[k] = -(prevp[k] + psevp[k]);
          qv[k] = qv[k] + work2[k] * dtcld;
          qci1[k] = c_max(qci1[k] - (praut[k] + pracw[k] + psacw[k]) * dtcld, Real(0.));
          qrs1[k] = c_max(qrs1[k] + (praut[k] + pracw[k] + prevp[k] + psacw[k]) * dtcld, Real(0.));
          qrs2[k] = c_max(qrs2[k] + psevp[k] * dtcld, Real(0.));
          const Real xlwork2 = -xl[k] * (prevp[k] + psevp[k]);
          t[k] = t[k] - xlwork2 / cpm[k] * dtcld;
        }
      }
      for (int k = 0; k < kte; ++k) {  // vectorised (libmvec log, exp)
        const Real tr = ttp / t[k];
        const Real logtr = vm::log(tr, li, ni);
        qs1[k] = c.psat * vm::exp(logtr * (xa) + xb * (Real(1.) - tr), li, ni);
        qs1[k] = c_min(qs1[k], Real(0.99) * pr[k]);
        qs1[k] = c.ep2 * qs1[k] / (pr[k] - qs1[k]);
        qs1[k] = c_max(qs1[k], c.qmin);
      }
      // pcond: condensational / evaporational rate of cloud water
      for (int k = 0; k < kte; ++k) {
        work11[k] = ((c_max(qv[k], c.qmin) - (qs1[k])) /
                     (Real(1.) + (xl[k]) * (xl[k]) / (c.rv * (cpm[k])) * (qs1[k]) / ((t[k]) * (t[k]))));
        work2[k] = qci1[k] + work11[k];
        pcond[k] = c_min(c_max(work11[k] / dtcld, Real(0.)), c_max(qv[k], Real(0.)) / dtcld);
        if (qci1[k] > Real(0.) && work11[k] < Real(0.)) pcond[k] = c_max(work11[k], -qci1[k]) / dtcld;
        qv[k] = qv[k] - pcond[k] * dtcld;
        qci1[k] = c_max(qci1[k] + pcond[k] * dtcld, Real(0.));
        t[k] = t[k] + pcond[k] * xl[k] / cpm[k] * dtcld;
      }
      for (int k = 0; k < kte; ++k) {
        if (qci1[k] <= c.qmin) qci1[k] = Real(0.0);
        if (qci2[k] <= c.qmin) qci2[k] = Real(0.0);
      }
    }
    for (int k = 0; k < kte; ++k) {
      th(i, j, k) = t[k] / pii(i, j, k);
      q(i, j, k) = qv[k];
      qc(i, j, k) = qci1[k];
      qi(i, j, k) = qci2[k];
      qr(i, j, k) = qrs1[k];
      qs(i, j, k) = qrs2[k];
    }
    rain(i, j) = rain_v;
    rainncv(i, j) = rainncv_v;
    snow(i, j) = snow_v;
    snowncv(i, j) = snowncv_v;
    sr(i, j) = sr_v;
  }

  // slope_wsm5: intercept / slope parameters and fall speeds (work1 = vt) of rain and snow
  CPPWRF_HD void slope(Real* rs1, Real* rs2, Real* rb1, Real* rb2, Real* r21, Real* r22, Real* r31, Real* r32,
                       Real* vt1, Real* vt2, const Real* q1, const Real* q2, const Real* dn, const Real* denfac,
                       const Real* t, int km) const {
    for (int k = 0; k < km; ++k) {
      const Real supcol = Real(273.15) - t[k];
      const Real n0sfac = c_max(c_min(libm::exp(Real(c.alpha) * supcol), Real(c.n0smax / c.n0s)), Real(1.));
      if (q1[k] <= Real(c.qcrmin)) {
        rs1[k] = c.rslopermax;
        rb1[k] = c.rsloperbmax;
        r21[k] = c.rsloper2max;
        r31[k] = c.rsloper3max;
      } else {
        rs1[k] = Real(1.) / std::sqrt(std::sqrt(c.pidn0r / (q1[k] * dn[k])));
        rb1[k] = libm::exp(libm::log(rs1[k]) * Real(c.bvtr));
        r21[k] = rs1[k] * rs1[k];
        r31[k] = r21[k] * rs1[k];
      }
      if (q2[k] <= Real(c.qcrmin)) {
        rs2[k] = c.rslopesmax;
        rb2[k] = c.rslopesbmax;
        r22[k] = c.rslopes2max;
        r32[k] = c.rslopes3max;
      } else {
        rs2[k] = Real(1.) / std::sqrt(std::sqrt(c.pidn0s * n0sfac / (q2[k] * dn[k])));
        rb2[k] = libm::exp(libm::log(rs2[k]) * Real(c.bvts));
        r22[k] = rs2[k] * rs2[k];
        r32[k] = r22[k] * rs2[k];
      }
      vt1[k] = c.pvtr * rb1[k] * denfac[k];
      vt2[k] = c.pvts * rb2[k] * denfac[k];
      if (q1[k] <= Real(0.0)) vt1[k] = Real(0.0);
      if (q2[k] <= Real(0.0)) vt2[k] = Real(0.0);
    }
  }
};

}  // namespace cppwrf::phys::wsm5
