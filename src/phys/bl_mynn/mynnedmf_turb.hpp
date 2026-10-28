#pragma once
// MYNN-EDMF turbulence closure (WRF v4.8.0 phys/MYNN-EDMF/module_bl_mynnedmf.F90): mym_level2 (bl_mynn_ess
// 0 and 1), boulac_length, mym_length (bl_mynn_mixlength 0, 1, 2), mym_initialize, mym_turbulence and
// mym_predict (bl_mynn_edmf <= 1, closure 2.5 / 2.6 / 3.0, tke_budget) and mym_condensation (bl_mynn_cloudpdf
// 0, 1, 2 and -1, -2). Columns are 1-based (k = 1 .. kte, interface arrays 1 .. kte+1); K bounds kte. The
// downdraft (edmf_dd) and SPP inputs are zero (bl_mynn_edmf_dd = 0, spp_pbl = 0).
#include "phys/bl_mynn/mynnedmf_common.hpp"

namespace cppwrf::phys::mynnedmf {

// mym_level2: dtl, dqw, dtv, gm, gh, sm, sh at the walls (k = 2 .. kte)
template <int K>
CPPWRF_HD inline void mym_level2(int kte, real ckmod, int bl_mynn_ess, const real* zw, const real* dz, real xland,
                                 real pblh, const real* u, const real* v, const real* thl, const real* thv,
                                 const real* thlv, const real* th, const real* p, const real* exner, const real* qke,
                                 const real* cldfra, const real* edmf_a, const real* edmf_w, const real* qw,
                                 const real* vt, const real* vq, real* dtl, real* dqw, real* dtv, real* gm, real* gh,
                                 real* sm, real* sh) {
  constexpr int kts = 1;
  constexpr real thvp = real(0.0);
  real ri[K + 2];
  dtl[kts] = zero;
  dqw[kts] = zero;
  dtv[kts] = zero;
  gm[kts] = zero;
  gh[kts] = zero;
  sm[kts] = zero;
  sh[kts] = zero;
  ri[kts] = zero;
  if (bl_mynn_ess == 0) {
    for (int k = kts + 1; k <= kte; ++k) {
      const real dzk = p5 * (dz[k] + dz[k - 1]);
      const real afk = dz[k] / (dz[k] + dz[k - 1]);
      const real abk = one - afk;
      real duz = (u[k] - u[k - 1]) * (u[k] - u[k - 1]) + (v[k] - v[k - 1]) * (v[k] - v[k - 1]);
      duz = duz / (dzk * dzk);
      const real dtz = (thl[k] - thl[k - 1]) / (dzk);
      const real dqz = (qw[k] - qw[k - 1]) / (dzk);
      const real vtt = one + vt[k] * abk + vt[k - 1] * afk;
      const real vqq = tv0 + vq[k] * abk + vq[k - 1] * afk;
      const real dtq = vtt * dtz + vqq * dqz;
      dtl[k] = dtz;
      dqw[k] = dqz;
      dtv[k] = dtq;
      gm[k] = duz;
      gh[k] = -dtq * gtr;
      ri[k] = -gh[k] / c_max(duz, real(1.0e-10));
    }
  } else {
    const real ugrid = std::sqrt(u[kts] * u[kts] + v[kts] * v[kts]);
    const real uonset = real(15.);
    const real taper = one - p95 * c_min(one, c_max(zero, ugrid - uonset) / twenty);
    const real clam0 = ((xland - real(1.5)) >= zero) ? real(1.8) * taper : real(0.08);
    for (int k = kts + 1; k <= kte; ++k) {
      const real dzk = p5 * (dz[k] + dz[k - 1]);
      real duz = (u[k] - u[k - 1]) * (u[k] - u[k - 1]) + (v[k] - v[k - 1]) * (v[k] - v[k - 1]);
      duz = duz / (dzk * dzk);
      const real dtz = (thl[k] - thl[k - 1]) / dzk;
      const real dqz = (qw[k] - qw[k - 1]) / dzk;
      const real wt = c_min(one, c_max(zero, zw[k] - (pblh + hundred)) / c_max(real(400.), p3 * pblh));
      real clam = clam0 * (one - wt) + p2 * wt;
      clam = clam * c_min(one, c_max(zero, zw[k] - twenty) / real(200.0));
      const real thv1 = thvp * thv[k] + (one - thvp) * thlv[k];
      const real thv0 = thvp * thv[k - 1] + (one - thvp) * thlv[k - 1];
      real tabs = th[k] * exner[k];
      real xl = xl_blend(tabs);
      real qsat = qsat_blend(tabs, p[k]);
      const real eth1 = th[k] * rexp(xl * qsat / (cp * tabs));
      tabs = th[k - 1] * exner[k - 1];
      xl = xl_blend(tabs);
      qsat = qsat_blend(tabs, p[k - 1]);
      const real eth0 = th[k - 1] * rexp(xl * qsat / (cp * tabs));
      const real mfi = c_max(zero, c_min(p2, (p5 * (edmf_a[k] + edmf_a[k - 1])) * (p5 * (edmf_w[k] + edmf_w[k - 1]))));
      real qkei = c_max(qkemin, p5 * (qke[k] + qke[k - 1]));
      qkei = c_min(p5, std::sqrt(p25 * qkei)) * p2;
      real cldfrai = p5 * (cldfra[k] + cldfra[k - 1]);
      const real cfac = c_max(real(0.01), c_min(real(2.7) * ((cldfrai - one) * (cldfrai - one)), one));
      cldfrai = c_min(p1, cldfrai * cfac) * ten;
      const real lambda = clam * c_max(mfi, qkei) * cldfrai;
      const real dtq = (thv1 - thv0) / dzk + lambda * c_min(zero, (eth1 - eth0) / dzk);
      dtl[k] = dtz;
      dqw[k] = dqz;
      dtv[k] = dtq;
      gm[k] = duz;
      gh[k] = -dtq * gtr;
      ri[k] = -gh[k] / c_max(duz, real(1.0e-10));
    }
  }
  for (int k = kts + 1; k <= kte; ++k) {
    const real a2fac = one / (one + ckmod * c_max(ri[k], zero));
    const real rfc = g1 / (g1 + g2);
    const real f1 = b1 * (g1 - c1) + three * a2 * a2fac * (one - c2) * (one - c5) + two * a1 * (three - two * c2);
    const real f2 = b1 * (g1 + g2) - three * a1 * (one - c2);
    const real rf1 = b1 * (g1 - c1) / f1;
    const real rf2 = b1 * g1 / f2;
    const real smc = a1 / (a2 * a2fac) * f1 / f2;
    const real shc = three * (a2 * a2fac) * (g1 + g2);
    const real ri1 = p5 / smc;
    const real ri2 = rf1 * smc;
    const real ri3 = four * rf2 * smc - two * ri2;
    const real ri4 = ri2 * ri2;
    const real rf = c_min(ri1 * (ri[k] + ri2 - std::sqrt(ri[k] * ri[k] - ri3 * ri[k] + ri4)), rfc);
    sh[k] = shc * (rfc - rf) / (one - rf);
    sm[k] = smc * (rf1 - rf) / (rf2 - rf) * sh[k];
  }
}

// boulac_length
template <int K>
CPPWRF_HD inline void boulac_length(int kte, const real* zw, const real* dz, const real* qtke, const real* theta,
                                    real* lb1, real* lb2) {
  constexpr int kts = 1;
  constexpr real Lmax = real(1500.);
  real dlu[K + 2], dld[K + 2];
  for (int iz = kts; iz <= kte; ++iz) {
    real zup = zero;
    dlu[iz] = zw[kte + 1] - zw[iz] - dz[iz] * p5;
    real zzz = zero;
    real zup_inf = zero;
    const real beta = gtr;
    real dzt, tl, bbb;
    if (iz < kte) {
      int found = 0;
      int izz = iz;
      while (found == 0) {
        if (izz < kte) {
          dzt = dz[izz];
          zup = zup - beta * theta[iz] * dzt;
          zup = zup + beta * (theta[izz + 1] + theta[izz]) * dzt * p5;
          zzz = zzz + dzt;
          if (qtke[iz] < zup && qtke[iz] >= zup_inf) {
            bbb = (theta[izz + 1] - theta[izz]) / dzt;
            if (bbb != real(0.)) {
              tl = (-beta * (theta[izz] - theta[iz]) +
                    std::sqrt(c_max(zero, (beta * (theta[izz] - theta[iz])) * (beta * (theta[izz] - theta[iz])) +
                                              two * bbb * beta * (qtke[iz] - zup_inf)))) /
                   bbb / beta;
            } else {
              if (theta[izz] != theta[iz]) tl = (qtke[iz] - zup_inf) / (beta * (theta[izz] - theta[iz]));
              else tl = zero;
            }
            dlu[iz] = zzz - dzt + tl;
            found = 1;
          }
          zup_inf = zup;
          izz = izz + 1;
        } else {
          found = 1;
        }
      }
    }
    real zdo = zero;
    real zdo_sup = zero;
    dld[iz] = zw[iz];
    zzz = zero;
    if (iz > kts) {
      int found = 0;
      int izz = iz;
      while (found == 0) {
        if (izz > kts) {
          dzt = dz[izz - 1];
          zdo = zdo + beta * theta[iz] * dzt;
          zdo = zdo - beta * (theta[izz - 1] + theta[izz]) * dzt * p5;
          zzz = zzz + dzt;
          if (qtke[iz] < zdo && qtke[iz] >= zdo_sup) {
            bbb = (theta[izz] - theta[izz - 1]) / dzt;
            if (bbb != real(0.)) {
              tl = (beta * (theta[izz] - theta[iz]) +
                    std::sqrt(c_max(zero, (beta * (theta[izz] - theta[iz])) * (beta * (theta[izz] - theta[iz])) +
                                              two * bbb * beta * (qtke[iz] - zdo_sup)))) /
                   bbb / beta;
            } else {
              if (theta[izz] != theta[iz]) tl = (qtke[iz] - zdo_sup) / (beta * (theta[izz] - theta[iz]));
              else tl = zero;
            }
            dld[iz] = zzz - dzt + tl;
            found = 1;
          }
          zdo_sup = zdo;
          izz = izz - 1;
        } else {
          found = 1;
        }
      }
    }
    dld[iz] = c_min(dld[iz], zw[iz + 1]);
    dlu[iz] = c_max(p1, dlu[iz] / (one + (dlu[iz] / Lmax)));
    dld[iz] = c_max(p1, dld[iz] / (one + (dld[iz] / Lmax)));
    lb1[iz] = c_min(dlu[iz], dld[iz]);
    lb2[iz] = std::sqrt(dlu[iz] * dld[iz]);
  }
  lb1[kte] = lb1[kte - 1];
  lb2[kte] = lb2[kte - 1];
}

// mym_length, bl_mynn_mixlength = 1 (nonlocal BouLac form)
template <int K>
CPPWRF_HD inline void mym_length1(int kte, real xland, const real* dz, real dx, const real* zw, real rmol, real fltv,
                                  const real* u1, const real* v1, const real* qke, const real* qpe, const real* dtv,
                                  real* el, real pblh, const real* theta, real* qkw, real* qtw, const real* edmf_w1,
                                  const real* edmf_a1) {
  constexpr int kts = 1;
  constexpr real qkw_elb_min = real(0.05);
  real qtke[K + 2], elBLmin[K + 2], elBLavg[K + 2], thetaw[K + 2];
  const real ugrid = std::sqrt(u1[kts] * u1[kts] + v1[kts] * v1[kts]);
  const real uonset = real(15.);
  const real wt_u1 = one - p2 * c_min(one, c_max(zero, ugrid - uonset) / forty);
  const real wt_u2 = one - p5 * c_min(one, c_max(zero, ugrid - uonset) / forty);
  const real wt_dx = one - c_min(one, (c_max(zero, dx - real(4000.)) / real(12000.)));
  const real cns = real(3.5);
  const real alp1 = real(0.23) * wt_dx + (one - wt_dx) * real(0.24);
  const real alp2 = p3 * wt_dx + (one - wt_dx) * real(0.40);
  const real alp3 = five * wt_u2;
  const real alp4 = ((xland - real(1.5)) >= zero) ? twenty * wt_u2 : real(15.0);
  const real alp5 = p3;
  const real alp6 = fifty;
  const real pblh2 = c_max(pblh, real(200.));
  real h1 = c_max(p3 * pblh2, real(300.));
  h1 = c_min(h1, real(600.));
  const real h2 = h1 * p5;
  qtke[kts] = c_max(p5 * qke[kts], p5 * qkemin);
  thetaw[kts] = theta[kts];
  qkw[kts] = std::sqrt(c_max(qke[kts], qkemin));
  qtw[kts] = std::sqrt(c_max(qke[kts] + qpe[kts], qkemin));
  for (int k = kts + 1; k <= kte; ++k) {
    const real afk = dz[k] / (dz[k] + dz[k - 1]);
    const real abk = one - afk;
    qkw[k] = std::sqrt(c_max(qke[k] * abk + qke[k - 1] * afk, qkemin));
    qtw[k] = std::sqrt(c_max((qke[k] + qpe[k]) * abk + (qke[k - 1] + qpe[k - 1]) * afk, qkemin));
    qtke[k] = c_max(p5 * (qkw[k] * qkw[k]), real(0.005));
    thetaw[k] = theta[k] * abk + theta[k - 1] * afk;
  }
  real elt = real(1.0e-5);
  real vsc = real(1.0e-5);
  int k = kts + 1;
  real zwk = zw[k];
  while (zwk <= pblh2 + c_max(real(200.), c_min(p3 * pblh2, real(600.)))) {
    const real dzk = p5 * (dz[k] + dz[k - 1]);
    const real qdz = c_min(c_max(qkw[k], real(0.01)), real(30.0)) * dzk;
    elt = elt + qdz * zwk;
    vsc = vsc + qdz;
    k = k + 1;
    zwk = zw[k];
  }
  real elt_max;
  if ((xland - real(1.5)) >= zero) elt_max = real(350.) + real(100.) * c_min(one, c_max(zero, ugrid - twenty) / thirty);
  else elt_max = real(400.);
  elt = c_min(c_max(alp1 * elt / vsc, real(8.)), elt_max);
  const real vflx = fltv - real(0.005);
  vsc = rpow(gtr * elt * c_max(vflx, zero), p333);
  el[kts] = zero;
  boulac_length<K>(kte, zw, dz, qtke, thetaw, elBLmin, elBLavg);
  for (k = kts + 1; k <= kte; ++k) {
    zwk = zw[k];
    // edmf_a_dd1 * edmf_w_dd1 = 0 (no downdrafts)
    real qkw_mf = c_max((p5 * ((edmf_a1[k] + edmf_a1[k - 1]))) * (p5 * ((edmf_w1[k] + edmf_w1[k - 1]))),
                        std::fabs(zero * zero));
    const real cpblh = c_min((zwk + p25 * pblh2) / pblh2, one);
    qkw_mf = cpblh * qkw_mf;
    const real alp3z = cpblh * alp3;
    real elb, elf;
    if (dtv[k] > real(0.0)) {
      const real bv = c_max(std::sqrt(gtr * dtv[k]), real(0.001));
      real elb_mf = alp6 * qkw_mf / bv;
      elb_mf = elb_mf / (one + (elb_mf / real(100.)));
      elb = alp2 * (c_max(qkw[k], qkw_elb_min)) / bv * (one + alp3z * std::sqrt(vsc / (bv * elt)));
      elb = c_max(elb, elb_mf);
      elb = c_min(elb, zwk);
      elf = one * c_max(qkw[k], qkw_elb_min) / bv;
      elf = c_max(elf, elb_mf);
      elBLavg[k] = c_max(alp5 * elBLavg[k], elb_mf);
    } else {
      elb = real(1.0e10);
      elf = elb;
      elBLavg[k] = p4 * elBLavg[k];
    }
    const real z_m = c_max(karman, zwk - four);
    real els, els1;
    if (rmol > real(0.0)) {
      els = karman * zwk / (one + cns * c_min(zwk * rmol, zmax));
      els1 = karman * z_m / (one + cns * c_min(zwk * rmol, zmax));
    } else {
      els = karman * zwk * rpow(one - alp4 * zwk * rmol, p2);
      els1 = karman * z_m * rpow(one - alp4 * zwk * rmol, p2);
    }
    const real wt = p5 * std::tanh((zwk - (pblh2 + h1)) / h2) + p5;
    const real el_unstab = els / (one + (els1 / elt));
    el[k] = c_min(el_unstab, elb);
    el[k] = c_min(el[k], elf);
    if ((xland - real(1.5)) >= zero) el[k] = el[k] * wt_u1;
    el[k] = el[k] * (one - wt) + elBLavg[k] * wt;
  }
}

// mym_length, bl_mynn_mixlength = 0 (original MYNN mixing length)
template <int K>
CPPWRF_HD inline void mym_length0(int kte, const real* dz, const real* zw, real rmol, real flt, real flq,
                                  const real* vt, const real* vq, const real* qke, const real* dtv, real* el, real pblh,
                                  real* qkw, real* qtw) {
  constexpr int kts = 1;
  constexpr real mindz = real(250.), maxdz = real(750.);
  const real cns = real(2.7), alp1 = real(0.23), alp2 = one, alp3 = five, alp4 = hundred;
  const real pblh2 = pblh + dz[kts];
  real h1 = c_max(p3 * pblh2, mindz);
  h1 = c_min(h1, maxdz);
  qkw[kts] = std::sqrt(c_max(qke[kts], qkemin));
  for (int k = kts + 1; k <= kte; ++k) {
    const real afk = dz[k] / (dz[k] + dz[k - 1]);
    const real abk = one - afk;
    qkw[k] = std::sqrt(c_max(qke[k] * abk + qke[k - 1] * afk, qkemin));
  }
  for (int k = kts; k <= kte; ++k) qtw[k] = qkw[k];
  real elt = real(1.0e-5);
  real vsc = real(1.0e-5);
  int k = kts + 1;
  real zwk = zw[k];
  while (zwk <= pblh2 + h1) {
    const real dzk = p5 * (dz[k] + dz[k - 1]);
    const real qdz = c_max(qkw[k] - qmin, real(0.03)) * dzk;
    elt = elt + qdz * zwk;
    vsc = vsc + qdz;
    k = k + 1;
    zwk = zw[k];
  }
  elt = alp1 * elt / vsc;
  const real vflx = (vt[kts] + one) * flt + (vq[kts] + tv0) * flq;
  vsc = rpow(gtr * elt * c_max(vflx, zero), p333);
  el[kts] = zero;
  for (k = kts + 1; k <= kte; ++k) {
    zwk = zw[k];
    real elb, elf;
    if (dtv[k] > zero) {
      const real bv = std::sqrt(gtr * dtv[k]);
      elb = alp2 * qkw[k] / bv * (one + alp3 / alp2 * std::sqrt(vsc / (bv * elt)));
      elf = alp2 * qkw[k] / bv;
    } else {
      elb = real(1.0e10);
      elf = elb;
    }
    real els;
    if (rmol > real(0.0)) els = karman * zwk / (one + cns * c_min(zwk * rmol, zmax));
    else els = karman * zwk * rpow(one - alp4 * zwk * rmol, p2);
    el[k] = c_min(elb / (elb / elt + elb / els + one), elf);
  }
}

// mym_length, bl_mynn_mixlength = 2 (local mixing length; the downdraft term edmf_a_dd1 * edmf_w_dd1 is 0)
template <int K>
CPPWRF_HD inline void mym_length2(int kte, const real* dz, const real* zw, real rmol, real fltv, const real* qke,
                                  const real* qpe, const real* dtv, real* el, real pblh, real* qkw, real* qtw,
                                  const real* edmf_w1, const real* edmf_a1) {
  constexpr int kts = 1;
  constexpr real ctau = real(1000.);
  const real cns = real(3.5), alp1 = real(0.22), alp2 = p3, alp3 = real(2.5), alp4 = five, alp5 = alp2, alp6 = fifty;
  real qtke[K + 2];
  const real pblh2 = c_max(pblh, real(200.));
  real h1 = c_max(p3 * pblh2, real(300.));
  h1 = c_min(h1, real(600.));
  const real h2 = h1 * p5;
  qtke[kts] = c_max(p5 * qke[kts], p5 * qkemin);
  qkw[kts] = std::sqrt(c_max(qke[kts], qkemin));
  qtw[kts] = std::sqrt(c_max(qke[kts] + qpe[kts], qkemin));
  for (int k = kts + 1; k <= kte; ++k) {
    const real afk = dz[k] / (dz[k] + dz[k - 1]);
    const real abk = one - afk;
    qkw[k] = std::sqrt(c_max(qke[k] * abk + qke[k - 1] * afk, qkemin));
    qtw[k] = std::sqrt(c_max((qke[k] + qpe[k]) * abk + (qke[k - 1] + qpe[k - 1]) * afk, qkemin));
    qtke[k] = p5 * (qkw[k] * qkw[k]);
  }
  real elt = real(1.0e-5);
  real vsc = real(1.0e-5);
  const real pblh_plus_ent = c_max(pblh + h1, real(100.));
  int k = kts + 1;
  real zwk = zw[k];
  while (zwk <= pblh_plus_ent) {
    const real dzk = p5 * (dz[k] + dz[k - 1]);
    const real qdz = c_min(c_max(qkw[k] - qmin, real(0.03)), real(30.0)) * dzk;
    elt = elt + qdz * zwk;
    vsc = vsc + qdz;
    k = k + 1;
    zwk = zw[k];
  }
  elt = c_min(c_max(alp1 * elt / vsc, ten), real(400.));
  const real vflx = fltv;
  vsc = rpow(gtr * elt * c_max(vflx, zero), p333);
  el[kts] = zero;
  for (k = kts + 1; k <= kte; ++k) {
    zwk = zw[k];
    const real dzk = p5 * (dz[k] + dz[k - 1]);
    const real qkw_mf = c_max((p5 * (edmf_a1[k] + edmf_a1[k - 1])) * (p5 * (edmf_w1[k] + edmf_w1[k - 1])),
                              std::fabs(zero * zero));
    real elb, elf, elb_mf;
    if (dtv[k] > real(0.0)) {
      const real bv = c_max(std::sqrt(gtr * dtv[k]), real(0.001));
      elb_mf = c_max(alp2 * qkw[k], alp6 * qkw_mf) / bv * (one + alp3 * std::sqrt(vsc / (bv * elt)));
      elb = c_min(c_max(alp5 * qkw[k], alp6 * qkw_mf) / bv, zwk);
      const real wstar = real(1.25) * rpow(gtr * pblh * c_max(vflx, real(1.0e-4)), p333);
      real tau_cloud = c_min(c_max(ctau * wstar / grav, real(30.)), real(150.));
      const real wt = p5 * std::tanh((zwk - (pblh2 + h1)) / h2) + p5;
      tau_cloud = tau_cloud * (real(1.) - wt) + fifty * wt;
      elf = c_min(c_max(tau_cloud * std::sqrt(c_min(qtke[k], forty)), alp6 * qkw_mf / bv), zwk);
    } else {
      const real wstar = real(1.25) * rpow(gtr * pblh * c_max(vflx, real(1.0e-4)), p333);
      real tau_cloud = c_min(c_max(ctau * wstar / grav, fifty), real(200.));
      const real wt = p5 * std::tanh((zwk - (pblh2 + h1)) / h2) + p5;
      tau_cloud = tau_cloud * (one - wt) + c_max(hundred, dzk * p25) * wt;
      elb = c_min(tau_cloud * std::sqrt(c_min(qtke[k], real(40.))), zwk);
      elf = elb;
      elb_mf = elb;
    }
    elf = elf / (one + (elf / real(800.)));
    elb_mf = c_max(elb_mf, real(0.01));
    real els;
    if (rmol > real(0.0)) els = karman * zwk / (one + cns * c_min(zwk * rmol, zmax));
    else els = karman * zwk * rpow(one - alp4 * zwk * rmol, p2);
    const real wt = p5 * std::tanh((zwk - (pblh2 + h1)) / h2) + p5;
    el[k] = std::sqrt((els * els) / (one + ((els * els) / (elt * elt)) + ((els * els) / (elb_mf * elb_mf))));
    el[k] = el[k] * (one - wt) + elf * wt;
  }
}

// mym_length (WRF's SELECT CASE (bl_mynn_mixlength), then the blend with the LES length scale); qkw, qtw inout.
// vt, vq, flt, flq are used by mixlength 0 only
template <int K>
CPPWRF_HD inline void mym_length(int kte, real xland, const real* dz, real dx, const real* zw, real rmol, real flt,
                                 real fltv, real flq, const real* vt, const real* vq, const real* u1, const real* v1,
                                 const real* qke, const real* qpe, const real* dtv, real* el, real pblh,
                                 const real* theta, real* qkw, real* qtw, real Psig_bl, int mixlength,
                                 const real* edmf_w1, const real* edmf_a1) {
  constexpr int kts = 1;
  if (mixlength == 0) {
    mym_length0<K>(kte, dz, zw, rmol, flt, flq, vt, vq, qke, dtv, el, pblh, qkw, qtw);
  } else if (mixlength == 2) {
    mym_length2<K>(kte, dz, zw, rmol, fltv, qke, qpe, dtv, el, pblh, qkw, qtw, edmf_w1, edmf_a1);
  } else if (mixlength == 1) {
    mym_length1<K>(kte, xland, dz, dx, zw, rmol, fltv, u1, v1, qke, qpe, dtv, el, pblh, theta, qkw, qtw, edmf_w1,
                   edmf_a1);
  }
  for (int k = kts + 1; k <= kte; ++k) {
    const real el_les = p25 * p5 * (dz[k] + dz[k - 1]);
    el[k] = el[k] * Psig_bl + (one - Psig_bl) * c_min(el_les, el[k]);
  }
}

// mym_initialize (pmz = phh = 1, flt = fltv = flq = 0 as its local defaults; rmol = 0)
template <int K>
CPPWRF_HD inline void mym_initialize(int kte, real ckmod, real taue, real xland, real pblh, const real* dz, real dx,
                                     const real* zw, const real* p, const real* exner, const real* u, const real* v,
                                     const real* thl, const real* qw, const real* theta, const real* thv,
                                     const real* thlv, real* sh, real* sm, real ust, real rmol, real* el, real* qke,
                                     real* tsq, real* qsq, real* cov, real Psig_bl, const real* cldfra_bl1,
                                     int bl_mynn_mixlength, int bl_mynn_ess, const real* edmf_w1, const real* edmf_a1,
                                     bool initialize_qke) {
  constexpr int kts = 1;
  constexpr real pmz = 1., phh = 1., flt = 0., fltv = 0., flq = 0.;
  real ql[K + 2], pdk[K + 2], pdt[K + 2], pdq[K + 2], pdc[K + 2], dtl[K + 2], dqw[K + 2], dtv[K + 2], gm[K + 2],
      gh[K + 2], qkw[K + 2], qtw[K + 2], vt[K + 2], vq[K + 2], qpe[K + 2];
  for (int k = kts; k <= kte; ++k) {
    ql[k] = zero;
    vt[k] = zero;
    vq[k] = zero;
    qkw[k] = zero;
    qtw[k] = zero;
  }
  (void)ql;
  const real ustm = c_max(p01, ust);
  mym_level2<K>(kte, ckmod, bl_mynn_ess, zw, dz, xland, pblh, u, v, thl, thv, thlv, theta, p, exner, qke, cldfra_bl1,
                edmf_a1, edmf_w1, qw, vt, vq, dtl, dqw, dtv, gm, gh, sm, sh);
  el[kts] = zero;
  if (initialize_qke) {
    qke[kts] = real(1.5) * (ustm * ustm) * rpow(b1 * pmz, p666);
    for (int k = kts + 1; k <= kte; ++k)
      qke[k] = qke[kts] * c_max((ust * real(700.) - zw[k]) / (ustm * real(700.)), p01);
  }
  real phm = phh * b2 / rpow(b1 * pmz, p333);
  tsq[kts] = c_min(three, phm * ((c_max(zero, flt) / ustm) * (c_max(zero, flt) / ustm)));
  qsq[kts] = c_min(real(5e-6), phm * ((c_max(zero, flq) / ustm) * (c_max(zero, flq) / ustm)));
  cov[kts] = phm * (flt / ustm) * (flq / ustm);
  qpe[kts] = c_min(three, c_max(zero, two * tsq[kts] * (gtr * gtr) * (taue * taue)));
  for (int k = kts + 1; k <= kte; ++k) {
    const real vkz = karman * zw[k];
    el[k] = vkz / (one + vkz / hundred);
    tsq[k] = zero;
    qsq[k] = zero;
    cov[k] = zero;
    qpe[k] = zero;
  }
  constexpr int lmax = 5;
  for (int l = 1; l <= lmax; ++l) {
    mym_length<K>(kte, xland, dz, dx, zw, rmol, flt, fltv, flq, vt, vq, u, v, qke, qpe, dtv, el, pblh, theta, qkw,
                  qtw, Psig_bl, bl_mynn_mixlength, edmf_w1, edmf_a1);
    for (int k = kts + 1; k <= kte; ++k) {
      const real elq = el[k] * qkw[k];
      pdk[k] = elq * (sm[k] * gm[k] + sh[k] * gh[k]);
      pdt[k] = elq * sh[k] * (dtl[k] * dtl[k]);
      pdq[k] = elq * sh[k] * (dqw[k] * dqw[k]);
      pdc[k] = elq * sh[k] * dtl[k] * dqw[k];
    }
    const real vkz = karman * p5 * dz[kts];
    const real elv = p5 * (el[kts + 1] + el[kts]) / vkz;
    if (initialize_qke) qke[kts] = one * (ustm * ustm) * rpow(b1 * pmz * elv, p666);
    phm = phh * b2 / rpow(b1 * pmz / (elv * elv), p333);
    tsq[kts] = c_min(three, phm * ((c_max(zero, flt) / ustm) * (c_max(zero, flt) / ustm)));
    qsq[kts] = c_min(real(5e-6), phm * ((c_max(zero, flq) / ustm) * (c_max(zero, flq) / ustm)));
    cov[kts] = phm * (flt / ustm) * (flq / ustm);
    for (int k = kts + 1; k <= kte - 1; ++k) {
      const real b1l = b1 * real(0.25) * (el[k + 1] + el[k]);
      const real tmpq = c_min(c_max(b1l * (pdk[k + 1] + pdk[k]), qkemin), real(125.));
      if (initialize_qke) qke[k] = rpow(tmpq, p666);
      real b2l;
      if (qke[k] <= zero) b2l = zero;
      else b2l = b2 * (b1l / b1) / std::sqrt(qke[k]);
      tsq[k] = c_min(three, c_max(zero, b2l * (pdt[k + 1] + pdt[k])));
      qsq[k] = c_min(real(5e-6), c_max(zero, b2l * (pdq[k + 1] + pdq[k])));
      cov[k] = b2l * (pdc[k + 1] + pdc[k]);
    }
  }
  if (initialize_qke) {
    qke[kts] = p5 * (qke[kts] + qke[kts + 1]);
    qke[kte] = qke[kte - 1];
  }
  tsq[kte] = tsq[kte - 1];
  qsq[kte] = qsq[kte - 1];
  cov[kte] = cov[kte - 1];
}

// mym_turbulence (tkeprod_dn = 0: no downdrafts); sh, sm are overwritten with the final stability functions.
// qshear / qbuoy (tke_budget = 1) may be null
template <int K>
CPPWRF_HD inline void mym_turbulence(int kte, real xland, real closure, real ckmod, real taue, const real* dz, real dx,
                                     const real* zw, const real* p, const real* exner, const real* u, const real* v,
                                     const real* thl, const real* thv, const real* thlv, const real* qw,
                                     const real* qke, const real* tsq, const real* qsq, const real* cov,
                                     const real* vt, const real* vq, real rmol, real flt, real fltv, real flq,
                                     real pblh, const real* theta, real* sh, real* sm, real* el, real* dfm, real* dfh,
                                     real* dfq, real* tcd, real* qcd, real* pdk, real* pdt, real* pdq, real* pdc,
                                     real* qshear, real* qbuoy, real Psig_bl, const real* cldfra_bl1,
                                     int bl_mynn_mixlength, int bl_mynn_ess, const real* edmf_w1,
                                     const real* edmf_a1, const real* TKEprod_up) {
  constexpr int kts = 1;
  constexpr real Prlimit_fre = real(6.0), Prlimit_sfc = real(4.0);
  real qkw[K + 2], qtw[K + 2], dtl[K + 2], dqw[K + 2], dtv[K + 2], gm[K + 2], gh[K + 2], qpe[K + 2];
  for (int k = kts; k <= kte; ++k) {
    qkw[k] = zero;
    qtw[k] = zero;
  }
  const real pblh2 = c_max(pblh, real(200.));
  if (closure < real(2.7)) {
    for (int k = kts; k <= kte; ++k) qpe[k] = zero;
  } else {
    for (int k = kts; k <= kte; ++k) {
      const real cpblh = c_min((zw[k] + p1 * pblh2) / (real(0.8) * pblh2), one);
      qpe[k] = c_min(three, c_max(zero, cpblh * two * tsq[k] * (gtr * gtr) * (taue * taue)));
    }
  }
  mym_level2<K>(kte, ckmod, bl_mynn_ess, zw, dz, xland, pblh, u, v, thl, thv, thlv, theta, p, exner, qke, cldfra_bl1,
                edmf_a1, edmf_w1, qw, vt, vq, dtl, dqw, dtv, gm, gh, sm, sh);
  mym_length<K>(kte, xland, dz, dx, zw, rmol, flt, fltv, flq, vt, vq, u, v, qke, qpe, dtv, el, pblh, theta, qkw, qtw,
                Psig_bl, bl_mynn_mixlength, edmf_w1, edmf_a1);
  for (int k = kts + 1; k <= kte; ++k) {
    const real dzk = p5 * (dz[k] + dz[k - 1]);
    const real afk = dz[k] / (dz[k] + dz[k - 1]);
    const real abk = one - afk;
    const real elsq = el[k] * el[k];
    real q3sq = qkw[k] * qkw[k];
    const real q2sq = b1 * elsq * (sm[k] * gm[k] + sh[k] * gh[k]);
    sh[k] = c_max(sh[k], real(1e-5));
    real duz = (u[k] - u[k - 1]) * (u[k] - u[k - 1]) + (v[k] - v[k - 1]) * (v[k] - v[k - 1]);
    duz = duz / (dzk * dzk);
    const real ri = -gh[k] / c_max(duz, real(1.0e-10));
    const real a2fac = one / (one + ckmod * c_max(ri, zero));
    real Prlim;
    if (ri >= one) Prlim = c_min(seven * ri, Prlimit_fre);
    else if (ri >= real(0.01) && ri <= one) Prlim = c_min(real(6.873) * ri + one / (real(6.873) * ri), Prlimit_fre);
    else Prlim = Prlimit_fre;
    const real wt = c_min(one, zw[k] / c_max(pblh, real(200.)));
    Prlim = wt * Prlim + (one - wt) * c_min(Prlimit_sfc, Prlim);
    const real gmel = gm[k] * elsq;
    const real ghel = gh[k] * elsq;
    real dlsq = elsq;
    if (q3sq / dlsq < -gh[k]) q3sq = -dlsq * gh[k];
    real qdiv, e1, e2, e3, e4, eden;
    if (q3sq < q2sq) {
      qdiv = std::sqrt(q3sq / q2sq);
      sh[k] = sh[k] * qdiv;
      sm[k] = sm[k] * qdiv;
      e1 = q3sq - e1c * ghel * a2fac * (qdiv * qdiv);
      e2 = q3sq - e2c * ghel * a2fac * (qdiv * qdiv);
      e3 = e1 + e3c * ghel * (a2fac * a2fac) * (qdiv * qdiv);
      e4 = e1 - e4c * ghel * a2fac * (qdiv * qdiv);
      eden = e2 * e4 + e3 * e5c * gmel * (qdiv * qdiv);
      eden = c_max(eden, real(1.0e-20));
    } else {
      e1 = q3sq - e1c * ghel * a2fac;
      e2 = q3sq - e2c * ghel * a2fac;
      e3 = e1 + e3c * ghel * (a2fac * a2fac);
      e4 = e1 - e4c * ghel * a2fac;
      eden = e2 * e4 + e3 * e5c * gmel;
      eden = c_max(eden, real(1.0e-20));
      qdiv = one;
      sm[k] = q3sq * a1 * (e3 - three * c1 * e4) / eden;
      sh[k] = q3sq * (a2 * a2fac) * (e2 + three * c1 * e5c * gmel) / eden;
    }
    constexpr real sh25max = four, sh25min = zero;
    if (sh[k] > sh25max) sh[k] = sh25max;
    if (sh[k] < sh25min) sh[k] = sh25min;
    const real shb = c_max(sh[k], real(0.004));
    sm[k] = c_min(sm[k], Prlim * shb);
    real gamt, gamq, gamv;
    if (closure >= real(3.0)) {  // level 3: the countergradient terms (Nakanishi and Niino 2009)
      real t2sq = qdiv * b2 * elsq * sh[k] * (dtl[k] * dtl[k]);
      real r2sq = qdiv * b2 * elsq * sh[k] * (dqw[k] * dqw[k]);
      real c2sq = qdiv * b2 * elsq * sh[k] * dtl[k] * dqw[k];
      real t3sq = c_max(tsq[k] * abk + tsq[k - 1] * afk, zero);
      real r3sq = c_max(qsq[k] * abk + qsq[k - 1] * afk, zero);
      real c3sq = cov[k] * abk + cov[k - 1] * afk;
      c3sq = std::copysign(c_min(std::fabs(c3sq), std::sqrt(t3sq * r3sq)), c3sq);
      const real vtt = one + vt[k] * abk + vt[k - 1] * afk;
      const real vqq = tv0 + vq[k] * abk + vq[k - 1] * afk;
      t2sq = vtt * t2sq + vqq * c2sq;
      r2sq = vtt * c2sq + vqq * r2sq;
      c2sq = c_max(vtt * t2sq + vqq * r2sq, real(0.0));
      t3sq = vtt * t3sq + vqq * c3sq;
      r3sq = vtt * c3sq + vqq * r3sq;
      c3sq = c_max(vtt * t3sq + vqq * r3sq, real(0.0));
      const real cw25 = e1 * (e2 + three * c1 * e5c * gmel * (qdiv * qdiv)) / (three * eden);
      dlsq = elsq;
      if (q3sq / dlsq < -gh[k]) q3sq = -dlsq * gh[k];
      // auh, aum, adh, adm, aeh, aem, Req, Rsl: no effect (Rsl is overwritten and unused)
      e2 = q3sq - e2c * ghel * a2fac * (qdiv * qdiv);
      e3 = q3sq + e3c * ghel * (a2fac * a2fac) * (qdiv * qdiv);
      e4 = q3sq - e4c * ghel * a2fac * (qdiv * qdiv);
      eden = e2 * e4 + e3 * e5c * gmel * (qdiv * qdiv);
      const real wden = cc3 * (gtr * gtr) * (dlsq * dlsq) / elsq * (qdiv * qdiv) *
                        (e2 * e4c * a2fac - e3c * e5c * gmel * (a2fac * a2fac) * (qdiv * qdiv));
      if (wden != zero) {
        const real clow = q3sq * (real(0.12) - cw25) * eden / wden;
        const real cupp = q3sq * (real(0.76) - cw25) * eden / wden;
        if (wden > zero) c3sq = c_min(c_max(c3sq, c2sq + clow), c2sq + cupp);
        else c3sq = c_max(c_min(c3sq, c2sq + clow), c2sq + cupp);
      }
      e1 = e2 + e5c * gmel * (qdiv * qdiv);
      eden = c_max(eden, real(1.0e-20));
      const real e6c = real(3.0) * (a2 * a2fac) * cc3 * gtr * dlsq / elsq;
      real enum_;
      if (t2sq >= zero) enum_ = c_max(qdiv * e6c * (t3sq - t2sq), real(0.0));
      else enum_ = c_min(qdiv * e6c * (t3sq - t2sq), real(0.0));
      gamt = -e1 * enum_ / eden;
      if (r2sq >= zero) enum_ = c_max(qdiv * e6c * (r3sq - r2sq), real(0.0));
      else enum_ = c_min(qdiv * e6c * (r3sq - r2sq), real(0.0));
      gamq = -e1 * enum_ / eden;
      enum_ = c_max(qdiv * e6c * (c3sq - c2sq), real(0.0));
      const real smd =
          dlsq * enum_ * gtr / eden * (qdiv * qdiv) * (e3c * (a2fac * a2fac) + e4c * a2fac) * a1 / (a2 * a2fac);
      gamv = e1 * enum_ * gtr / eden;
      sm[k] = sm[k] + smd;
      qdiv = one;
    } else {
      gamt = zero;
      gamq = zero;
      gamv = zero;
    }
    const real cldavg = p5 * (cldfra_bl1[k - 1] + cldfra_bl1[k]);
    const real mfmax = c_max(p5 * (edmf_a1[k - 1] + edmf_a1[k]) * p5 * (edmf_w1[k - 1] + edmf_w1[k]),
                             std::fabs(zero * zero));
    sm[k] = c_max(sm[k], real(0.04) * c_min(real(10.) * mfmax, one));
    sh[k] = c_max(sh[k], real(0.04) * c_min(real(10.) * mfmax, one));
    sm[k] = c_max(sm[k], real(0.04) * c_min(cldavg, p5));
    sh[k] = c_max(sh[k], real(0.04) * c_min(cldavg, p5));
    if (closure >= real(2.7) && closure < real(3.0)) sm[k] = c_max(sm[k], c_min(p2, c_max(zero, three * qpe[k])));
    const real elq = el[k] * qkw[k];
    const real elqt = el[k] * qtw[k];
    const real elh = elq * qdiv;
    pdk[k] = elq * (sm[k] * gm[k] + sh[k] * gh[k] + gamv) + p5 * zero + p5 * TKEprod_up[k];
    pdt[k] = elh * (sh[k] * dtl[k] + gamt) * dtl[k];
    pdq[k] = elh * (sh[k] * dqw[k] + gamq) * dqw[k];
    pdc[k] = elh * (sh[k] * dtl[k] + gamt) * dqw[k] * p5 + elh * (sh[k] * dqw[k] + gamq) * dtl[k] * p5;
    tcd[k] = elq * gamt;
    qcd[k] = elq * gamq;
    dfm[k] = elqt * sm[k] / dzk;
    dfh[k] = elq * sh[k] / dzk;
    dfq[k] = dfm[k];
    if (qshear) {  // tke_budget = 1
      qshear[k] = elq * sm[k] * gm[k];
      qbuoy[k] = elq * (sh[k] * gh[k] + gamv) + p5 * zero + p5 * TKEprod_up[k];
    }
  }
  dfm[kts] = zero;
  dfh[kts] = zero;
  dfq[kts] = zero;
  pdk[kts] = zero;
  pdt[kts] = zero;
  pdq[kts] = zero;
  pdc[kts] = zero;
  tcd[kts] = zero;
  qcd[kts] = zero;
  tcd[kte] = zero;
  qcd[kte] = zero;
  for (int k = kts; k <= kte - 1; ++k) {
    const real dzk = dz[k];
    tcd[k] = (tcd[k + 1] - tcd[k]) / (dzk);
    qcd[k] = (qcd[k + 1] - qcd[k]) / (dzk);
  }
}

// mym_predict (bl_mynn_edmf <= 1: implicit form); qwt / qdiss (tke_budget = 1) may be null
template <int K>
CPPWRF_HD inline void mym_predict(int kte, real closure, real delt, const real* dz, real ust, real pmz, real phh,
                                  const real* el, const real* dfq, const real* rho, real* pdk, real* pdt, real* pdq,
                                  real* pdc, real* qke, real* tsq, real* qsq, real* cov, const real* s_aw1,
                                  const real* s_awqke1, const real* s_awqsq1, int bl_mynn_edmf, int bl_mynn_edmf_tke,
                                  real* qwt, real* qdiss) {
  constexpr int kts = 1;
  real qkw[K + 2], bp[K + 2], rp[K + 2], df3q[K + 2], dtz[K + 2], a[K + 2], b[K + 2], c[K + 2], d[K + 2], x[K + 2],
      rhoinv[K + 2], rhoz[K + 3], kqdz[K + 3], kmdz[K + 3];
  const real onoff = (bl_mynn_edmf_tke == 0) ? zero : one;
  const real vkz = karman * p5 * dz[kts];
  for (int k = kts; k <= kte; ++k) {
    qkw[k] = std::sqrt(c_max(qke[k], qkemin));
    df3q[k] = Sqfac * dfq[k];
    dtz[k] = delt / dz[k];
  }
  rhoz[kts] = rho[kts];
  rhoinv[kts] = one / rho[kts];
  kqdz[kts] = rhoz[kts] * df3q[kts];
  kmdz[kts] = rhoz[kts] * dfq[kts];
  for (int k = kts + 1; k <= kte; ++k) {
    rhoz[k] = (rho[k] * dz[k - 1] + rho[k - 1] * dz[k]) / (dz[k - 1] + dz[k]);
    rhoz[k] = c_max(rhoz[k], real(1E-4));
    rhoinv[k] = one / c_max(rho[k], real(1E-4));
    kqdz[k] = rhoz[k] * df3q[k];
    kmdz[k] = rhoz[k] * dfq[k];
  }
  rhoz[kte + 1] = rhoz[kte];
  kqdz[kte + 1] = rhoz[kte + 1] * df3q[kte];
  kmdz[kte + 1] = rhoz[kte + 1] * dfq[kte];
  if (bl_mynn_edmf == 1) {
    for (int k = kts + 1; k <= kte - 1; ++k) {
      kqdz[k] = c_max(kqdz[k], p5 * s_aw1[k]);
      kqdz[k] = c_max(kqdz[k], -p5 * (s_aw1[k] - s_aw1[k + 1]));
      kmdz[k] = c_max(kmdz[k], p5 * s_aw1[k]);
      kmdz[k] = c_max(kmdz[k], -p5 * (s_aw1[k] - s_aw1[k + 1]));
    }
  }
  const real pdk1 = two * (ust * ust * ust) * pmz / (vkz);
  pdk[kts] = pdk1 - pdk[kts + 1];
  pdt[kts] = pdt[kts + 1];
  pdq[kts] = pdq[kts + 1];
  pdc[kts] = pdc[kts + 1];
  (void)phh;
  for (int k = kts; k <= kte - 1; ++k) {
    const real b1l = b1 * p5 * (el[k + 1] + el[k]);
    bp[k] = two * qkw[k] / b1l;
    rp[k] = pdk[k + 1] + pdk[k];
  }
  bp[kte] = zero;
  for (int k = kts; k <= kte - 1; ++k) {
    a[k] = -dtz[k] * kqdz[k] * rhoinv[k] + real(0.5) * dtz[k] * rhoinv[k] * s_aw1[k] * onoff;
    b[k] = real(1.) + dtz[k] * (kqdz[k] + kqdz[k + 1]) * rhoinv[k] +
           real(0.5) * dtz[k] * rhoinv[k] * (s_aw1[k] - s_aw1[k + 1]) * onoff + bp[k] * delt;
    c[k] = -dtz[k] * kqdz[k + 1] * rhoinv[k] - real(0.5) * dtz[k] * rhoinv[k] * s_aw1[k + 1] * onoff;
    d[k] = rp[k] * delt + qke[k] + dtz[k] * rhoinv[k] * (s_awqke1[k] - s_awqke1[k + 1]) * onoff;
  }
  a[kte] = real(0.);
  b[kte] = real(1.);
  c[kte] = real(0.);
  d[kte] = qke[kte];
  tridiag2<K>(kte, a, b, c, d, x);
  for (int k = kts; k <= kte; ++k) {
    qke[k] = c_max(x[k], qkemin);
    qke[k] = c_min(qke[k], real(150.));
  }
  if (qwt) {  // tke_budget = 1
    real tke_up[K + 2], dzinv[K + 2];
    for (int k = kts; k <= kte; ++k) {
      tke_up[k] = p5 * qke[k];
      dzinv[k] = one / dz[k];
    }
    int k = kts;
    qwt[k] = dzinv[k] * ((kqdz[k + 1] * (tke_up[k + 1] - tke_up[k]) - kqdz[k] * tke_up[k]) +
                         p5 * rhoinv[k] *
                             (s_aw1[k + 1] * tke_up[k + 1] + (s_aw1[k + 1] - s_aw1[k]) * tke_up[k] +
                              (s_awqke1[k] - s_awqke1[k + 1])) *
                             onoff);
    for (k = kts + 1; k <= kte - 1; ++k)
      qwt[k] = dzinv[k] * ((kqdz[k + 1] * (tke_up[k + 1] - tke_up[k]) - kqdz[k] * (tke_up[k] - tke_up[k - 1])) +
                           p5 * rhoinv[k] *
                               (s_aw1[k + 1] * tke_up[k + 1] + (s_aw1[k + 1] - s_aw1[k]) * tke_up[k] -
                                s_aw1[k] * tke_up[k - 1] + (s_awqke1[k] - s_awqke1[k + 1])) *
                               onoff);
    k = kte;
    qwt[k] = dzinv[k] * (-kqdz[k] * (tke_up[k] - tke_up[k - 1]) +
                         p5 * rhoinv[k] * (-s_aw1[k] * tke_up[k] - s_aw1[k] * tke_up[k - 1] + s_awqke1[k]) * onoff);
    for (k = kts; k <= kte; ++k) qdiss[k] = bp[k] * tke_up[k];
  }
  if (closure > real(2.5)) {
    for (int k = kts; k <= kte - 1; ++k) {
      const real b2l = b2 * p5 * (el[k + 1] + el[k]);
      bp[k] = two * qkw[k] / b2l;
      rp[k] = pdq[k + 1] + pdq[k];
    }
    bp[kte] = zero;
    for (int k = kts; k <= kte - 1; ++k) {
      a[k] = -dtz[k] * kmdz[k] * rhoinv[k] + real(0.5) * dtz[k] * rhoinv[k] * s_aw1[k] * onoff;
      b[k] = real(1.) + dtz[k] * (kmdz[k] + kmdz[k + 1]) * rhoinv[k] +
             real(0.5) * dtz[k] * rhoinv[k] * (s_aw1[k] - s_aw1[k + 1]) * onoff + bp[k] * delt;
      c[k] = -dtz[k] * kmdz[k + 1] * rhoinv[k] - real(0.5) * dtz[k] * rhoinv[k] * s_aw1[k + 1] * onoff;
      d[k] = rp[k] * delt + qsq[k] + dtz[k] * rhoinv[k] * (s_awqsq1[k] - s_awqsq1[k + 1]) * onoff;
    }
    a[kte] = -one;
    b[kte] = one;
    c[kte] = zero;
    d[kte] = zero;
    tridiag2<K>(kte, a, b, c, d, x);
    for (int k = kts; k <= kte; ++k) {
      qsq[k] = c_max(x[k], real(1e-17));
      qsq[k] = c_min(qsq[k], real(6e-6));
    }
  } else {
    for (int k = kts; k <= kte - 1; ++k) {
      real b2l;
      if (qkw[k] <= zero) b2l = zero;
      else b2l = b2 * real(0.25) * (el[k + 1] + el[k]) / qkw[k];
      qsq[k] = c_min(real(6e-6), c_max(real(1e-17), b2l * (pdq[k + 1] + pdq[k])));
    }
    qsq[kte] = qsq[kte - 1];
  }
  if (closure >= real(2.7)) {  // prognostic temperature variance
    for (int k = kts; k <= kte - 1; ++k) {
      const real b2l = b2 * p5 * (el[k + 1] + el[k]);
      bp[k] = two * qkw[k] / b2l;
      rp[k] = pdt[k + 1] + pdt[k];
    }
    for (int k = kts; k <= kte - 1; ++k) {
      a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
      b[k] = one + dtz[k] * (kmdz[k] + kmdz[k + 1]) * rhoinv[k] + bp[k] * delt;
      c[k] = -dtz[k] * kmdz[k + 1] * rhoinv[k];
      d[k] = rp[k] * delt + c_max(zero, tsq[k]);
    }
    a[kte] = real(-1.);
    b[kte] = one;
    c[kte] = zero;
    d[kte] = zero;
    tridiag2<K>(kte, a, b, c, d, x);
    for (int k = kts; k <= kte; ++k) tsq[k] = c_min(three, c_max(zero, x[k]));
  } else {
    for (int k = kts; k <= kte - 1; ++k) {
      real b2l;
      if (qkw[k] <= zero) b2l = zero;
      else b2l = b2 * real(0.25) * (el[k + 1] + el[k]) / qkw[k];
      tsq[k] = c_min(three, c_max(zero, b2l * (pdt[k + 1] + pdt[k])));
    }
    tsq[kte] = tsq[kte - 1];
  }
  if (closure >= real(3.0)) {  // prognostic covariance
    for (int k = kts; k <= kte - 1; ++k) {
      const real b2l = b2 * p5 * (el[k + 1] + el[k]);
      bp[k] = two * qkw[k] / b2l;
      rp[k] = pdc[k + 1] + pdc[k];
    }
    for (int k = kts; k <= kte - 1; ++k) {
      a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
      b[k] = one + dtz[k] * (kmdz[k] + kmdz[k + 1]) * rhoinv[k] + bp[k] * delt;
      c[k] = -dtz[k] * kmdz[k + 1] * rhoinv[k];
      d[k] = rp[k] * delt + cov[k];
    }
    a[kte] = real(-1.);
    b[kte] = one;
    c[kte] = zero;
    d[kte] = zero;
    tridiag2<K>(kte, a, b, c, d, x);
    for (int k = kts; k <= kte; ++k) cov[k] = x[k];
  } else {
    for (int k = kts; k <= kte - 1; ++k) {
      real b2l;
      if (qkw[k] <= zero) b2l = zero;
      else b2l = b2 * real(0.25) * (el[k + 1] + el[k]) / qkw[k];
      cov[k] = b2l * (pdc[k + 1] + pdc[k]);
    }
    cov[kte] = cov[kte - 1];
  }
}

// mym_condensation, bl_mynn_cloudpdf 0 (original MYNN), 1 / -1 (Kuwano), 2 / -2 (Chaboureau and Bechtold;
// negative values: no subgrid clouds); qs is the snow content, sh and el the previous step's (cloudpdf 1)
template <int K>
CPPWRF_HD inline void mym_condensation(int kte, const real* dz, const real* zw, real xland, const real* thl,
                                       const real* qw, const real* qc, const real* qi, const real* qs, const real* p,
                                       const real* exner, const real* tsq, const real* qsq, const real* cov,
                                       const real* sh, const real* el, int cloudpdf, real* qc_bl1, real* qi_bl1,
                                       real* cldfra_bl1, real pblh1, real* vt, real* vq, const real* th, real* sgm,
                                       real closure) {
  constexpr int kts = 1;
  constexpr real qlim_sfc = real(0.008), qlim_pbl = real(0.020), qlim_trp = real(0.025);
  constexpr real clim_sfc = real(0.010), clim_pbl = real(0.020), clim_trp = real(0.030);
  constexpr real rhcrit = real(0.83), rhmax = real(1.10);
  real alp[K + 2], bet[K + 2], a[K + 2], b[K + 2], q1[K + 2], rh[K + 2];
  int k;
  for (k = kte - 3; k >= kts; --k) {
    const real theta1 = th[k];
    const real theta2 = th[k + 2];
    const real ht1 = real(44307.692) * (one - rpow(p[k] / real(101325.), real(0.190)));
    const real ht2 = real(44307.692) * (one - rpow(p[k + 2] / real(101325.), real(0.190)));
    if ((((theta2 - theta1) / (ht2 - ht1)) < real(10.) / real(1500.)) && (ht1 < real(19000.)) && (ht1 > real(4000.)))
      break;
  }
  const int k_tropo = (kts + 2 > k + 2) ? kts + 2 : k + 2;
  if (cloudpdf == 0 || cloudpdf == 1 || cloudpdf == -1) {
    for (k = kts; k <= kte - 1; ++k) {
      const real t = th[k] * exner[k];
      const real esat = esat_blend(t);
      const real qsl = ep_2 * esat / c_max(real(1.e-4), (p[k] - ep_3 * esat));
      const real dqsl = qsl * ep_2 * xlv / (r_d * (t * t));
      alp[k] = one / (one + dqsl * xlvcp);
      bet[k] = dqsl * exner[k];
      if (cloudpdf == 0) {
        const real t3sq = c_max(tsq[k], zero);
        real r3sq = c_max(qsq[k], zero);
        real c3sq = cov[k];
        c3sq = std::copysign(c_min(std::fabs(c3sq), std::sqrt(t3sq * r3sq)), c3sq);
        r3sq = r3sq + (bet[k] * bet[k]) * t3sq - real(2.0) * bet[k] * c3sq;
        sgm[k] = std::sqrt(c_max(r3sq, real(1.0e-10)));
      } else {
        const real dzk = (k == kts) ? real(0.5) * dz[k] : dz[k];
        const int km1 = (k - 1 > kts) ? k - 1 : kts;
        const real dth = real(0.5) * (thl[k + 1] + thl[k]) - real(0.5) * (thl[k] + thl[km1]);
        const real dqw = real(0.5) * (qw[k + 1] + qw[k]) - real(0.5) * (qw[k] + qw[km1]);
        const real grad = dqw / dzk - bet[k] * (dth / dzk);
        sgm[k] = std::sqrt(c_max((alp[k] * alp[k] * c_max(el[k] * el[k], p1) * b2 * c_max(sh[k], real(0.03))) / four *
                                     (grad * grad),
                                 real(1.0e-10)));
      }
      const real qmq = qw[k] - qsl;
      q1[k] = qmq / sgm[k];
      cldfra_bl1[k] = real(0.5) * (one + std::erf(q1[k] * rr2));
      const real q1k = q1[k];
      const real eq1 = rrp * std::exp(real(-0.5) * q1k * q1k);
      const real qll = c_max(cldfra_bl1[k] * q1k + eq1, zero);
      const real ql = alp[k] * sgm[k] * qll;
      const real liq_frac = c_min(one, c_max(zero, (t - real(240.0)) / real(29.0)));
      qc_bl1[k] = liq_frac * ql;
      qi_bl1[k] = (one - liq_frac) * ql;
      const real q2p = xlvcp / exner[k];
      const real pt = thl[k] + q2p * ql;
      const real qt = one + p608 * qw[k] - (real(1.) + p608) * (qc_bl1[k] + qi_bl1[k]) * cldfra_bl1[k];
      const real rac = alp[k] * (cldfra_bl1[k] - qll * eq1) * (q2p * qt - (real(1.) + p608) * pt);
      vt[k] = qt - one - rac * bet[k];
      vq[k] = p608 * pt - tv0 + rac;
    }
  } else if (cloudpdf == 2 || cloudpdf == -2) {
  const real pblh2 = c_max(ten, pblh1);
  for (k = kts; k <= kte - 1; ++k) {
    const real zagl = zw[k] + p5 * dz[k];
    const real t = th[k] * exner[k];
    const real xl = xl_blend(t);
    const real qsat_tk = qsat_blend(t, p[k]);
    rh[k] = c_max(c_min(rhmax, qw[k] / c_max(real(1E-10), qsat_tk)), real(0.001));
    const real dqsl = qsat_tk * ep_2 * xlv / (r_d * (t * t));
    alp[k] = one / (one + dqsl * xlvcp);
    bet[k] = dqsl * exner[k];
    const real rsl = xl * qsat_tk / (r_v * (t * t));
    const real cpm = cp + qw[k] * cpv;
    a[k] = one / (one + xl * rsl / cpm);
    b[k] = a[k] * rsl;
    const real qw_pert = qw[k] + qw[k] * p5 * zero * real(0);  // spp_pbl = 0
    const real qmq = qw_pert - qsat_tk;
    real r3sq = c_max(qsq[k], zero);
    const real t3sq = c_max(tsq[k], zero);
    real c3sq = cov[k];
    c3sq = std::copysign(c_min(std::fabs(c3sq), std::sqrt(t3sq * r3sq)), c3sq);
    if (closure < real(2.7)) r3sq = r3sq;
    else if (closure < real(3.0)) r3sq = r3sq + (bet[k] * bet[k]) * t3sq;
    else r3sq = r3sq + (bet[k] * bet[k]) * t3sq - two * bet[k] * c3sq;
    sgm[k] = c_max(real(1e-13), std::sqrt(r3sq));
    sgm[k] = c_min(sgm[k], qw[k] * p5);
    real wt = c_min(one, c_max(zero, dz[k] - real(100.)) / real(500.));
    sgm[k] = sgm[k] + sgm[k] * p2 * wt;
    real sgmq = sgm[k];
    real sgmc = sgm[k];
    wt = c_min(one, c_max(zero, (zagl - (pblh2 + real(10.)))) / real(250.));
    real clim = clim_pbl * (one - wt) + clim_trp * wt;
    real zsl = c_min(real(150.), c_max(real(50.), p1 * pblh2));
    wt = c_min(one, c_max(zero, zagl - zsl) / real(150.));
    clim = clim_sfc * (one - wt) + clim * wt;
    sgmc = c_max(sgmc, qw[k] * clim);
    sgmc = c_max(real(1e-13), sgmc);
    if (qmq >= zero) sgmc = c_max(real(0.02) * qw[k], sgmc);
    q1[k] = qmq / sgmc;
    real wt2 = c_min(one, c_max(zero, zagl - pblh2) / real(300.));
    if ((qi[k] + qs[k]) > real(1.e-10) && (zagl > pblh2)) {
      real rh_adj = c_min(rhmax, rhcrit + wt2 * real(0.037) * (c_max(zero, ten + std::log10(qi[k] + qs[k]))));
      rh_adj = c_max(rh[k], rh_adj);
      const real q1_rh = real(-3.) + three * (rh_adj - rhcrit) / (one - rhcrit);
      q1[k] = c_max(q1_rh, q1[k]);
    }
    if (qc[k] > real(1.e-5) && (zagl > pblh2)) {
      real rh_adj = c_min(rhmax, rhcrit + wt2 * real(0.07) * (c_max(zero, five + std::log10(qc[k]))));
      rh_adj = c_max(rh[k], rh_adj);
      const real q1_rh = real(-3.) + three * (rh_adj - rhcrit) / (one - rhcrit);
      q1[k] = c_max(q1_rh, q1[k]);
    }
    real q1k = q1[k];
    wt2 = c_min(one, c_max(zero, (zagl - (pblh1 - real(100.))) / real(200.)));
    const real cldfra_qsq0 = c_max(zero, c_min(one, p5 + real(0.35) * std::atan(real(3.6) * (q1k + real(0.05)))));
    const real cldfra_qsq1 = c_max(zero, c_min(one, p5 + real(0.39) * std::atan(real(1.6) * (q1k + real(0.55)))));
    const real cldfra_qsq = cldfra_qsq0 * (one - wt2) + cldfra_qsq1 * wt2;
    const real cldfra_rh = zero;
    cldfra_bl1[k] = c_max(cldfra_qsq, cldfra_rh);
    wt = c_min(one, c_max(zero, (zagl - (pblh2 + ten))) / real(300.));
    real qlim = qlim_pbl * (one - wt) + qlim_trp * wt;
    zsl = c_min(real(150.), c_max(real(50.), p1 * pblh2));
    wt = c_min(one, c_max(zero, zagl - zsl) / real(200.));
    qlim = qlim_sfc * (one - wt) + qlim * wt;
    sgmq = c_max(sgmq, qw[k] * qlim);
    real ql_water = c_min(sgmq, real(0.025) * qw[k]) * cldfra_bl1[k];
    real ql_ice = c_min(sgmq, real(0.025) * qw[k]) * cldfra_bl1[k];
    if (cldfra_bl1[k] < real(0.001)) {
      ql_ice = zero;
      ql_water = zero;
      cldfra_bl1[k] = zero;
    }
    const real liq_frac = c_min(one, c_max(zero, (t - tice) / (tliq - tice)));
    qc_bl1[k] = liq_frac * ql_water;
    qi_bl1[k] = (one - liq_frac) * ql_ice;
    if (k >= k_tropo) {
      cldfra_bl1[k] = zero;
      qc_bl1[k] = zero;
      qi_bl1[k] = zero;
    }
    if ((xland - real(1.5)) >= zero) q1k = c_max(q1[k], real(-2.5));
    else q1k = c_max(q1[k], real(-2.0));
    real Fng;
    if (q1k >= one) Fng = one;
    else if (q1k >= real(-1.7) && q1k < one) Fng = rexp(real(-0.4) * (q1k - one));
    else if (q1k >= real(-2.5) && q1k < real(-1.7)) Fng = three + rexp(real(-3.8) * (q1k + real(1.7)));
    else Fng = c_min(real(23.9) + rexp(real(-1.6) * (q1k + real(2.5))), real(60.));
    real cfmax = c_min(cldfra_bl1[k], real(0.6));
    zsl = c_min(c_max(real(25.), p1 * pblh2), hundred);
    wt = c_min(zagl / zsl, one);
    cfmax = cfmax * wt;
    const real bb = b[k] * t / th[k];
    const real qww = one + p608 * qw[k];
    const real alpha = p608 * th[k];
    const real beta = (th[k] / t) * (xl / cp) - real(1.61) * th[k];
    vt[k] = qww - cfmax * beta * bb * Fng - one;
    vq[k] = alpha + cfmax * beta * a[k] * Fng - tv0;
  }
  }
  if (cloudpdf < 0) {
    for (k = kts; k <= kte - 1; ++k) {
      cldfra_bl1[k] = zero;
      qc_bl1[k] = zero;
      qi_bl1[k] = zero;
    }
  }
  vt[kte] = vt[kte - 1];
  vq[kte] = vq[kte - 1];
  qc_bl1[kte] = zero;
  qi_bl1[kte] = zero;
  cldfra_bl1[kte] = zero;
}

}  // namespace cppwrf::phys::mynnedmf
