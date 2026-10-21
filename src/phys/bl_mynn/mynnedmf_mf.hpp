#pragma once
// MYNN-EDMF dynamic multi-plume mass flux (WRF v4.8.0 module_bl_mynnedmf.F90 DMP_mf, condensation_edmf) for
// WRF's defaults: env_subs = .false. (no subsidence/detrainment tendencies), bl_mynn_edmf_tke = 0,
// bl_mynn_mixscalars = 0, bl_mynn_mixnumcon = 0, no aerosols / chemistry, spp_pbl = 0. The plume-size loops of
// NUP = 8 trips: the area loop uses glibc's libmvec pow in WRF's build (vm::pow with n = 8), the sum for cn
// (a reduction) the scalar pow.
#include "phys/bl_mynn/mynnedmf_common.hpp"

namespace cppwrf::phys::mynnedmf {

CPPWRF_HD inline void condensation_edmf(real QT, real THL, real P, real zagl, real& THV, real& QC) {
  constexpr int niter = 50;
  const real diff = real(1.e-6);
  const real EXN = rpow(P / p1000mb, rcp);
  QC = zero;
  for (int ni = 1; ni <= niter; ++ni) {
    const real T = EXN * THL + xlvcp / EXN * QC;
    const real QS = qsat_blend(T, P);
    const real QCOLD = QC;
    QC = p5 * QC + p5 * c_max((QT - QS), zero);
    if (std::fabs(QC - QCOLD) < diff) break;
  }
  const real T = EXN * THL + xlvcp / EXN * QC;
  const real QS = qsat_blend(T, P);
  QC = c_max(QT - QS, zero);
  if (zagl < real(100.)) QC = zero;
  THV = (THL + xlvcp / EXN * QC) * (one + QT * (rvovrd - one) - rvovrd * QC);
}

struct MfOut {
  real maxwidth, maxmf, ztop, excess_h, excess_q;
  int ktop;
};

// DMP_mf: qt1 = sqw, qv1 = sqv, qc1 = sqc (specific contents); ts = theta at the surface (th_sfc)
template <int K>
CPPWRF_HD inline void dmp_mf(int kte, real dt, const real* zw1, const real* dz1, const real* pres1, const real* rho1,
                             int momentum_opt, const real* u1, const real* v1, const real* w1, const real* th1,
                             const real* thl1, const real* thv1, const real* tk1, const real* qt1, const real* qv1,
                             const real* qc1, const real* qke1, const real* qsq1, const real* ex1, real* vt1, real* vq1,
                             real ust, real flt, real fltv, real flq, real pblh, int kpbl, real dx, real landsea,
                             real ts, real* edmf_a1, real* edmf_w1, real* edmf_qt1, real* edmf_thl1, real* edmf_ent1,
                             real* edmf_qc1, real* s_aw1, real* s_awthl1, real* s_awqt1, real* s_awqv1, real* s_awqc1,
                             real* s_awu1, real* s_awv1, real* qc_bl1, real* cldfra_bl1, real Psig_shcu, MfOut& o,
                             real* tkeprod_up, const real* el1) {
  constexpr int kts = 1;
  constexpr int nup = 8;
  constexpr real Atot = real(0.10), lmax = real(1000.), lmin = real(300.), dcut = real(1.2);
  constexpr real fluxportion = real(0.75);
  constexpr real pgfac = real(0.00);
  (void)dt; (void)ust; (void)th1;
  real UPW[K + 3][nup + 1], UPTHL[K + 3][nup + 1], UPQT[K + 3][nup + 1], UPQC[K + 3][nup + 1], UPA[K + 3][nup + 1],
      UPU[K + 3][nup + 1], UPV[K + 3][nup + 1], ENT[K + 2][nup + 1];
  real exneri[K + 2], dzi[K + 2], rhoi[K + 2], qti[K + 2], edmf_th1[K + 2];
  int ktop_plume[nup + 1];
  int ktop = 0;
  o.ztop = zero;
  o.maxmf = zero;
  o.maxwidth = zero;
  o.excess_h = zero;
  o.excess_q = zero;
  for (int k = kts; k <= kte + 1; ++k)
    for (int ip = 1; ip <= nup; ++ip) {
      UPW[k][ip] = zero;
      UPTHL[k][ip] = zero;
      UPQT[k][ip] = zero;
      UPA[k][ip] = zero;
      UPU[k][ip] = zero;
      UPV[k][ip] = zero;
      UPQC[k][ip] = zero;
    }
  for (int k = kts; k <= kte; ++k)
    for (int ip = 1; ip <= nup; ++ip) ENT[k][ip] = real(0.001);
  for (int k = kts; k <= kte; ++k) {
    edmf_a1[k] = zero;
    edmf_w1[k] = zero;
    edmf_qt1[k] = zero;
    edmf_thl1[k] = zero;
    edmf_ent1[k] = zero;
    edmf_qc1[k] = zero;
  }
  for (int k = kts; k <= kte + 1; ++k) {
    s_aw1[k] = zero;
    s_awthl1[k] = zero;
    s_awqt1[k] = zero;
    s_awqv1[k] = zero;
    s_awqc1[k] = zero;
    s_awu1[k] = zero;
    s_awv1[k] = zero;
  }
  real nup2 = nup;
  for (int k = kts; k <= kte; ++k) tkeprod_up[k] = zero;

  real maxw = zero;
  real cloud_base = real(9000.0);
  real qkebl = zero;
  int k50 = 0;  // WRF leaves it undefined when the first mass level is above 50 m
  for (int k = 1; k <= kte - 1; ++k) {
    const real zagl = zw1[k] + p5 * dz1[k];
    if (zagl > (pblh + real(500.))) break;
    real wpbl = w1[k];
    if (w1[k] < zero) wpbl = two * w1[k];
    maxw = c_max(maxw, std::fabs(wpbl));
    if (zagl <= real(50.)) k50 = k;
    if (k <= kpbl) qkebl = qkebl + qke1[k];
    const real qc_sgs = c_max(qc1[k], qc_bl1[k]);
    if ((qc_sgs > real(1E-5)) && (cldfra_bl1[k] >= p5) && cloud_base == real(9000.0)) cloud_base = zw1[k];
  }
  qkebl = qkebl / real(kpbl);
  maxw = c_max(zero, maxw - one);
  real Psig_w = c_max(zero, one - maxw);
  Psig_w = c_min(Psig_w, Psig_shcu);
  real fltv2 = fltv;
  if (Psig_w < real(1e-2) && fltv > zero) fltv2 = real(-1.) * fltv;

  bool superadiabatic = false;
  real hux = ((landsea - real(1.5)) >= zero) ? real(-0.0024) : real(-0.0030);
  const real tvs = ts * (one + p608 * qv1[kts]);
  if ((thv1[kts] - tvs) / (p5 * dz1[kts]) < hux) {
    superadiabatic = true;
    const int kend = (2 > k50) ? 2 : k50;
    for (int k = 2; k <= kend; ++k) {
      hux = real(-0.0018);
      const real dthvdz = (thv1[k] - thv1[k - 1]) / (p5 * (dz1[k] + dz1[k - 1]));
      if (dthvdz < hux) {
        superadiabatic = true;
      } else {
        superadiabatic = false;
        break;
      }
    }
  }
  const real maxwidth_dx = c_min(dx * dcut, lmax);
  const real maxwidth_pbl = c_min(real(1.1) * pblh, lmax);
  real maxwidth_cld;
  if ((landsea - real(1.5)) < zero) maxwidth_cld = c_min(lmax, c_max(real(0.5) * cloud_base, real(300.)));
  else maxwidth_cld = c_min(lmax, c_max(real(0.8) * cloud_base, real(300.)));
  const real wspd_pbl = std::sqrt(c_max(u1[kts] * u1[kts] + v1[kts] * v1[kts], real(0.01)));
  real maxwidth_flx;
  if ((landsea - real(1.5)) < zero)
    maxwidth_flx = c_max(c_min(real(1000.) * (real(0.6) * std::tanh((fltv - real(0.040)) / real(0.04)) + p5), real(1000.)), zero);
  else
    maxwidth_flx = c_max(c_min(real(1000.) * (real(0.6) * std::tanh((fltv - real(0.007)) / real(0.02)) + p5), real(1000.)), zero);
  real maxwidth = c_min(maxwidth_dx, maxwidth_pbl);
  maxwidth = c_min(maxwidth, maxwidth_cld);
  maxwidth = c_min(maxwidth, maxwidth_flx);
  const real minwidth = lmin;
  if (maxwidth <= minwidth) {
    nup2 = 0;
    maxwidth = zero;
  }
  real exc_heat = zero, exc_moist = zero;
  if (flt > zero && fltv2 > real(0.008) && (maxwidth > minwidth) && superadiabatic) {
    real cn = zero;
    const real d = real(-1.9);
    const real dl = (maxwidth - minwidth) / real(nup - 1);
    for (int ip = 1; ip <= nup; ++ip) {
      const real l = minwidth + dl * real(ip - 1);
      cn = cn + rpow(l, d) * (l * l) / (dx * dx) * dl;  // a reduction: not vectorised by gfortran
    }
    const real C = Atot / cn;
    real acfac;
    if ((landsea - real(1.5)) < zero) acfac = p5 * std::tanh((fltv2 - real(0.02)) / real(0.05)) + p5;
    else acfac = p5 * std::tanh((fltv2 - real(0.012)) / real(0.03)) + p5;
    const real ac_wsp = one - c_min((c_max(wspd_pbl - real(13.0), zero)) / real(12.), one);
    acfac = c_min(acfac, ac_wsp);
    for (int ip = 1; ip <= nup; ++ip) {
      const real l = minwidth + dl * real(ip - 1);
      const real N = C * vm::pow(l, d, ip - 1, nup);
      UPA[1][ip] = N * l * l / (dx * dx) * dl;
      UPA[1][ip] = UPA[1][ip] * acfac;
    }
    const real z0 = real(50.);
    const real pwmin = real(0.1), pwmax = real(0.4);
    const real wstar = c_max(real(1.E-2), rpow(gtr * fltv2 * pblh, p333));
    const real qstar = c_max(flq, real(1e-5)) / wstar;
    const real thstar = flt / wstar;
    const real csigma = real(1.34);
    real exc_fac = ((landsea - real(1.5)) >= zero) ? real(9.28) : real(0.58);
    exc_fac = exc_fac * ac_wsp;
    const real sigmaW = csigma * wstar * rpow(z0 / pblh, p333) * (one - real(0.8) * z0 / pblh);
    const real sigmaQT = csigma * qstar * rpow(z0 / pblh, p333);
    const real sigmaTH = csigma * thstar * rpow(z0 / pblh, p333);
    const real wmin0 = c_min(sigmaW * pwmin, real(0.1));
    const real wmax = c_min(sigmaW * pwmax, real(0.5));
    const real dzs = dz1[kts] + dz1[kts + 1];
    for (int ip = 1; ip <= nup; ++ip) {
      UPA[2][ip] = UPA[1][ip];
      UPW[2][ip] = wmin0 + real(ip) / real(nup) * (wmax - wmin0);
      UPU[2][ip] = (u1[kts] * dz1[kts + 1] + u1[kts + 1] * dz1[kts]) / dzs;
      UPV[2][ip] = (v1[kts] * dz1[kts + 1] + v1[kts + 1] * dz1[kts]) / dzs;
      UPQC[2][ip] = zero;
      exc_heat = exc_fac * UPW[2][ip] * sigmaTH / sigmaW;
      UPTHL[2][ip] = (thl1[kts] * dz1[kts + 1] + thl1[kts + 1] * dz1[kts]) / dzs + exc_heat;
      exc_moist = exc_fac * UPW[2][ip] * sigmaQT / sigmaW;
      UPQT[2][ip] = (qt1[kts] * dz1[kts + 1] + qt1[kts + 1] * dz1[kts]) / dzs + exc_moist;
    }
    rhoi[kts] = rho1[kts];
    dzi[kts] = p5 * dz1[kts];
    exneri[kts] = ex1[kts];
    for (int k = kts + 1; k <= kte; ++k) {
      rhoi[k] = (rho1[k - 1] * dz1[k] + rho1[k] * dz1[k - 1]) / (dz1[k - 1] + dz1[k]);
      dzi[k] = (dz1[k - 1] * dz1[k] + dz1[k] * dz1[k - 1]) / (dz1[k - 1] + dz1[k]);
      exneri[k] = (ex1[k - 1] * dz1[k] + ex1[k] * dz1[k - 1]) / (dz1[k - 1] + dz1[k]);
    }
    const real dxsa = one - c_min(c_max((real(12000.0) - dx) / (real(12000.0) - real(3000.0)), zero), one);
    for (int ip = 1; ip <= nup; ++ip) ktop_plume[ip] = 0;
    for (int ip = 1; ip <= nup; ++ip) {
      real QCn = zero;
      const real l = minwidth + dl * real(ip - 1);
      int k;
      for (k = kts + 1; k <= kte - 2; ++k) {
        const real zagl = zw1[k] + p5 * dz1[k];
        const real wmin = p3 + l * real(0.0005);
        real entfac = real(0.21) * c_min(real(1.57), c_max(real(1.30), std::sqrt(qkebl)));
        const real wt2 = c_min(one, c_max(zero, zagl - pblh) / real(500.));
        entfac = entfac * (one - wt2) + wt2 * real(0.33);
        ENT[k][ip] = entfac / (c_min(c_max(UPW[k - 1][ip], wmin), one) * l);
        ENT[k][ip] = c_max(ENT[k][ip], real(0.0003));
        if (zw1[k] >= c_min(pblh + real(1500.), real(4000.)))
          ENT[k][ip] = ENT[k][ip] + (zw1[k] - c_min(pblh + real(1500.), real(4000.))) * real(5.0E-6);
        ENT[k][ip] = ENT[k][ip] * (one - zero);  // spp_pbl = 0
        ENT[k][ip] = c_min(ENT[k][ip], real(0.9) / dz1[k]);
        const real Uk = (u1[k] * dz1[k + 1] + u1[k + 1] * dz1[k]) / (dz1[k + 1] + dz1[k]);
        const real Ukm1 = (u1[k - 1] * dz1[k] + u1[k] * dz1[k - 1]) / (dz1[k - 1] + dz1[k]);
        const real Vk = (v1[k] * dz1[k + 1] + v1[k + 1] * dz1[k]) / (dz1[k + 1] + dz1[k]);
        const real Vkm1 = (v1[k - 1] * dz1[k] + v1[k] * dz1[k - 1]) / (dz1[k - 1] + dz1[k]);
        const real EntExp = ENT[k][ip] * dz1[k];
        const real EntExm = EntExp * p333;
        const real QTn = UPQT[k][ip] * (one - EntExp) + qt1[k] * EntExp;
        const real THLn = UPTHL[k][ip] * (one - EntExp) + thl1[k] * EntExp;
        const real Un = UPU[k][ip] * (one - EntExm) + u1[k] * EntExm + dxsa * pgfac * (Uk - Ukm1);
        const real Vn = UPV[k][ip] * (one - EntExm) + v1[k] * EntExm + dxsa * pgfac * (Vk - Vkm1);
        const real Pk = (pres1[k] * dz1[k + 1] + pres1[k + 1] * dz1[k]) / (dz1[k] + dz1[k + 1]);
        real THVn;
        condensation_edmf(QTn, THLn, Pk, zw1[k + 1], THVn, QCn);
        const real THVk = (thv1[k] * dz1[k + 1] + thv1[k + 1] * dz1[k]) / (dz1[k + 1] + dz1[k]);
        const real B = grav * (THVn / THVk - one);
        const real BCOEFF = (B > real(0.)) ? real(0.15) : real(0.2);
        real Wn;
        if (UPW[k - 1][ip] < real(0.2))
          Wn = UPW[k][ip] + (-(real(2.) * ENT[k][ip] * UPW[k][ip]) + BCOEFF * B / c_max(UPW[k][ip], real(0.2))) *
                                c_min(dz1[k], real(250.));
        else
          Wn = UPW[k][ip] + (-(real(2.) * ENT[k][ip] * UPW[k][ip]) + BCOEFF * B / UPW[k][ip]) * c_min(dz1[k], real(250.));
        if (Wn > UPW[k][ip] + c_min(real(1.25) * dz1[k] / real(200.), two))
          Wn = UPW[k][ip] + c_min(real(1.25) * dz1[k] / real(200.), two);
        if (Wn < UPW[k][ip] - c_min(real(1.25) * dz1[k] / real(200.), two))
          Wn = UPW[k][ip] - c_min(real(1.25) * dz1[k] / real(200.), two);
        Wn = c_min(c_max(Wn, zero), real(3.0));
        if (k == kts + 1 && Wn < real(1e-8)) {
          nup2 = 0;
          break;
        }
        if (Wn > real(0.)) {
          UPW[k + 1][ip] = Wn;
          UPTHL[k + 1][ip] = THLn;
          UPQT[k + 1][ip] = QTn;
          UPQC[k + 1][ip] = QCn;
          UPU[k + 1][ip] = Un;
          UPV[k + 1][ip] = Vn;
          UPA[k + 1][ip] = UPA[k][ip];
          ktop = (ktop > k) ? ktop : k;
        } else {
          break;
        }
      }
      ktop_plume[ip] = k;
      o.excess_h = c_max(o.excess_h, exc_heat);
      o.excess_q = c_max(o.excess_q, exc_moist);
    }
  } else {
    nup2 = real(0.);
  }
  ktop = (ktop < kte - 1) ? ktop : kte - 1;
  if (ktop == 0) o.ztop = zero;
  else o.ztop = zw1[ktop];

  if (nup2 > 0) {
    for (int ip = 1; ip <= nup; ++ip)
      for (int k = kts; k <= kte - 1; ++k) {
        s_aw1[k] = s_aw1[k] + rhoi[k] * UPA[k][ip] * UPW[k][ip] * Psig_w;
        s_awthl1[k] = s_awthl1[k] + rhoi[k] * UPA[k][ip] * UPW[k][ip] * UPTHL[k][ip] * Psig_w;
        s_awqt1[k] = s_awqt1[k] + rhoi[k] * UPA[k][ip] * UPW[k][ip] * UPQT[k][ip] * Psig_w;
        const real qc_plume = UPQC[k][ip];
        s_awqc1[k] = s_awqc1[k] + rhoi[k] * UPA[k][ip] * UPW[k][ip] * qc_plume * Psig_w;
        s_awqv1[k] = s_awqt1[k] - s_awqc1[k];
      }
    if (momentum_opt > 0)
      for (int ip = 1; ip <= nup; ++ip)
        for (int k = kts; k <= kte - 1; ++k) {
          s_awu1[k] = s_awu1[k] + rhoi[k] * UPA[k][ip] * UPW[k][ip] * UPU[k][ip] * Psig_w;
          s_awv1[k] = s_awv1[k] + rhoi[k] * UPA[k][ip] * UPW[k][ip] * UPV[k][ip] * Psig_w;
        }
    real flx1;
    if (std::fabs(s_aw1[kts + 1]) > real(1e-8))
      flx1 = c_max(s_aw1[kts + 1] * (thv1[kts] - thv1[kts + 1]) / dzi[kts + 1], real(1.0e-6));
    else
      flx1 = zero;
    const real flt2 = c_max(fltv, zero);
    if (flx1 > fluxportion * flt2 / dz1[kts] && flx1 > zero) {
      const real adjustment = c_max(real(0.01), fluxportion * flt2 / dz1[kts] / flx1);
      for (int k = kts; k <= kte + 1; ++k) {
        s_aw1[k] = s_aw1[k] * adjustment;
        s_awthl1[k] = s_awthl1[k] * adjustment;
        s_awqt1[k] = s_awqt1[k] * adjustment;
        s_awqc1[k] = s_awqc1[k] * adjustment;
        s_awqv1[k] = s_awqv1[k] * adjustment;
      }
      if (momentum_opt > 0)
        for (int k = kts; k <= kte + 1; ++k) {
          s_awu1[k] = s_awu1[k] * adjustment;
          s_awv1[k] = s_awv1[k] * adjustment;
        }
      for (int k = kts; k <= kte + 1; ++k)
        for (int ip = 1; ip <= nup; ++ip) UPA[k][ip] = UPA[k][ip] * adjustment;
    }
    {
      const int k = kts;
      for (int ip = 1; ip <= nup; ++ip) {
        edmf_a1[k] = edmf_a1[k] + UPA[k + 1][ip];
        edmf_w1[k] = edmf_w1[k] + UPA[k + 1][ip] * p5 * UPW[k + 1][ip];
        edmf_qt1[k] = edmf_qt1[k] + UPA[k + 1][ip] * UPQT[k + 1][ip];
        edmf_thl1[k] = edmf_thl1[k] + UPA[k + 1][ip] * UPTHL[k + 1][ip];
        edmf_ent1[k] = edmf_ent1[k] + UPA[k + 1][ip] * ENT[k + 1][ip];
        edmf_qc1[k] = edmf_qc1[k] + UPA[k + 1][ip] * UPQC[k + 1][ip];
      }
    }
    for (int ip = 1; ip <= nup; ++ip)
      for (int k = kts + 1; k <= ktop_plume[ip] - 1; ++k) {
        const real dd = dzi[k + 1] + dzi[k];
        const real upak = (UPA[k + 1][ip] * dzi[k] + UPA[k][ip] * dzi[k + 1]) / dd;
        edmf_a1[k] = edmf_a1[k] + upak;
        edmf_w1[k] = edmf_w1[k] + upak * (UPW[k + 1][ip] * dzi[k] + UPW[k][ip] * dzi[k + 1]) / dd;
        edmf_qt1[k] = edmf_qt1[k] + upak * (UPQT[k + 1][ip] * dzi[k] + UPQT[k][ip] * dzi[k + 1]) / dd;
        edmf_thl1[k] = edmf_thl1[k] + upak * (UPTHL[k + 1][ip] * dzi[k] + UPTHL[k][ip] * dzi[k + 1]) / dd;
        edmf_ent1[k] = edmf_ent1[k] + upak * (ENT[k + 1][ip] * dzi[k] + ENT[k][ip] * dzi[k + 1]) / dd;
        if (std::fabs(UPQC[k][ip]) < real(1e-10) && UPQC[k + 1][ip] > zero)
          edmf_qc1[k] = edmf_qc1[k] + upak * UPQC[k + 1][ip];
        else
          edmf_qc1[k] = edmf_qc1[k] + upak * (UPQC[k + 1][ip] * dzi[k] + UPQC[k][ip] * dzi[k + 1]) / dd;
      }
    for (int ip = 1; ip <= nup; ++ip) {
      const int k = ktop_plume[ip];
      edmf_a1[k] = edmf_a1[k] + UPA[k][ip];
      edmf_w1[k] = edmf_w1[k] + UPA[k][ip] * UPW[k][ip];
      edmf_qt1[k] = edmf_qt1[k] + UPA[k][ip] * UPQT[k][ip];
      edmf_thl1[k] = edmf_thl1[k] + UPA[k][ip] * UPTHL[k][ip];
      edmf_ent1[k] = edmf_ent1[k] + UPA[k][ip] * ENT[k][ip];
      edmf_qc1[k] = edmf_qc1[k] + UPA[k][ip] * UPQC[k][ip];
    }
    for (int k = kts; k <= kte - 1; ++k) {
      if (edmf_a1[k] > real(0.)) {
        edmf_w1[k] = edmf_w1[k] / edmf_a1[k];
        edmf_qt1[k] = edmf_qt1[k] / edmf_a1[k];
        edmf_thl1[k] = edmf_thl1[k] / edmf_a1[k];
        edmf_ent1[k] = edmf_ent1[k] / edmf_a1[k];
        edmf_qc1[k] = edmf_qc1[k] / edmf_a1[k];
        edmf_a1[k] = edmf_a1[k] * Psig_w;
      }
    }
    for (int k = kts; k <= kte - 1; ++k) {
      if (edmf_a1[k] * edmf_w1[k] > o.maxmf) o.maxmf = edmf_a1[k] * edmf_w1[k];
      const real aw = std::fabs(edmf_w1[k]);
      tkeprod_up[k] = (aw * aw * aw) * edmf_a1[k] / (b1 * c_max(p5 * (el1[k] + el1[k + 1]), real(0.1)));
    }
    for (int k = kts; k <= kte - 1; ++k) {
      edmf_th1[k] = edmf_thl1[k] + xlvcp / ex1[k] * edmf_qc1[k];
      qti[k + 1] = (qt1[k] * dz1[k + 1] + qt1[k + 1] * dz1[k]) / (dz1[k + 1] + dz1[k]);
    }
    qti[kts] = qt1[kts];
    int cb_check = 0;
    real mf_at_cb = zero;
    for (int k = kts + 1; k <= kte - 2; ++k) {
      if (k > ktop) break;
      if (edmf_qc1[k] > zero) {
        if (cb_check == 0) {
          mf_at_cb = edmf_a1[k] * edmf_w1[k];
          cb_check = 1;
        }
        const real Aup = edmf_a1[k];
        const real QTp = edmf_qt1[k];
        const real QCp = edmf_qc1[k];
        const real dqwdz = (qti[k] - qti[k + 1]) / dz1[k];
        const real xl = xl_blend(tk1[k]);
        const real qsat_tk = qsat_blend(tk1[k], pres1[k]);
        const real rsl = xl * qsat_tk / (r_v * (tk1[k] * tk1[k]));
        const real cpm = cp + qt1[k] * cpv;
        const real a = one / (one + xl * rsl / cpm);
        const real b9 = a * rsl;
        const real q2p = xlvcp / ex1[k];
        const real pt = thl1[k] + q2p * QCp * Aup;
        const real bb = b9 * tk1[k] / pt;
        const real qww = one + p608 * qt1[k];
        const real alpha = p608 * pt;
        const real beta = pt * xl / (tk1[k] * cp) - real(1.61) * pt;
        const real tauc = real(1800.);
        real sigq = five * Aup * edmf_w1[k] * c_max(real(4e-3), tauc * c_max(zero, dqwdz));
        sigq = c_max(sigq, three * Aup * (QTp - qt1[k]));
        const real qmq = qt1[k] - qsat_tk;
        if (qmq > zero) sigq = c_min(qsat_tk * real(0.015), sigq);
        real Q1 = qmq / sigq;
        real cf_mf;
        if ((landsea - real(1.5)) >= zero) {
          const real cffac = c_min(c_max(zero, mf_at_cb - real(0.05)) / real(0.05), one);
          cf_mf = c_min(c_max(p5 + real(0.36) * std::atan(real(1.8) * (Q1 + p2 + cffac)), real(0.01)), one);
          cf_mf = c_max(cf_mf, real(1.2) * Aup);
        } else {
          cf_mf = c_min(c_max(p5 + real(0.36) * std::atan(real(1.8) * (Q1 + p2)), real(0.01)), one);
          cf_mf = c_max(cf_mf, real(1.8) * Aup);
        }
        const real qcfac = real(1.20) + real(1.0) * c_min(c_max(zero, mf_at_cb - real(0.04)) / real(0.05), one);
        real qc_mf;
        if ((landsea - real(1.5)) >= zero) {
          qc_mf = qcfac * (QCp * Aup);
        } else {
          if (QCp * Aup > real(5e-5)) qc_mf = real(1.86) * (QCp * Aup) - real(2.2e-5);
          else qc_mf = real(1.20) * (QCp * Aup);
        }
        const real wt2 = one - c_min(one, c_max(zero, (cf_mf - thirty)) / (hundred - thirty));
        const real qc_mf_min = wt2 * qsat_tk * real(0.01) * cf_mf + (one - wt2) * qsat_tk * real(0.025) * cf_mf;
        qc_mf = c_max(qc_mf, qc_mf_min);
        const real cf_strat = cldfra_bl1[k];
        const real qc_strat = qc_bl1[k];
        cldfra_bl1[k] = c_max(cf_mf, cf_strat);
        const real pct_mf = cf_mf / cldfra_bl1[k];
        qc_bl1[k] = qc_mf * pct_mf + (one - pct_mf) * qc_strat;
        Q1 = c_max(Q1, real(-2.25));
        real Fng;
        if (Q1 >= one) Fng = one;
        else if (Q1 >= real(-1.7) && Q1 < one) Fng = rexp(-p4 * (Q1 - one));
        else if (Q1 >= real(-2.5) && Q1 < real(-1.7)) Fng = three + rexp(real(-3.8) * (Q1 + real(1.7)));
        else Fng = c_min(real(23.9) + rexp(real(-1.6) * (Q1 + real(2.5))), real(60.));
        vt1[k] = qww - (real(1.5) * Aup) * beta * bb * Fng - one;
        vq1[k] = alpha + (real(1.5) * Aup) * beta * a * Fng - tv0;
      }
    }
  }
  if (ktop > 0) {
    real maxqc = edmf_qc1[1];
    for (int k = 2; k <= ktop; ++k) maxqc = (edmf_qc1[k] > maxqc) ? edmf_qc1[k] : maxqc;
    if (maxqc < real(1.E-8)) o.maxmf = real(-1.) * o.maxmf;
  }
  o.maxwidth = maxwidth;
  o.ktop = ktop;
}

}  // namespace cppwrf::phys::mynnedmf
