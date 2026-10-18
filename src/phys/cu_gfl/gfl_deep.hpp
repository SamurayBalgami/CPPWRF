#pragma once

// Grell-Freitas-Li (GFL) deep convection kernel for one column.
// Ported statement-by-statement from WRF v4.8.0 phys/GFL/module_cu_gfl_deep.F90
// (cu_gfl_deep_run, get_partition_liq_ice, get_melting_profile, cup_dd_moisture,
// cup_dd_edt, cup_forcing_ens_3d, cup_output_ens_3d, neg_check).

#include "phys/cu_gfl/gfl_common.hpp"

namespace cppwrf::phys::cu_gfl {

#define R(x) static_cast<real>(x)

// ---------------------------------------------------------------------------
// get_partition_liq_ice (lines 5816-5913 of module_cu_gfl_deep.F90)
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void get_partition_liq_ice(
    int ierr, const real* tn, const real* po_cup, real* p_liq_ice,
    real* melting_layer, int ktf, int kte, bool is_deep) {
  using C = GflConstants<real>;
  constexpr real t1 = R(276.16);

  for (int k = 1; k <= kte; ++k) {
    p_liq_ice[k] = R(1.0);
    melting_layer[k] = R(0.0);
  }

  if (C::melt_glac && is_deep) {
    if (ierr == 0) {
      for (int k = 1; k <= ktf; ++k) {
        if (tn[k] <= C::t_ice) {
          p_liq_ice[k] = R(0.0);
        } else if (tn[k] > C::t_ice && tn[k] < C::t_0) {
          real r = (tn[k] - C::t_ice) / (C::t_0 - C::t_ice);
          p_liq_ice[k] = r * r;
        } else {
          p_liq_ice[k] = R(1.0);
        }
      }

      for (int k = 1; k <= ktf; ++k) {
        if (tn[k] <= C::t_0 + R(1.0)) {
          melting_layer[k] = R(0.0);
        } else if (tn[k] > C::t_0 + R(1.0) && tn[k] < t1) {
          real r = (tn[k] - C::t_0 + R(1.0)) / (t1 - C::t_0 + R(1.0));
          melting_layer[k] = r * r;
        } else {
          melting_layer[k] = R(1.0);
        }
        melting_layer[k] = melting_layer[k] * (R(1.0) - melting_layer[k]);
      }

      real norm = R(0.0);
      for (int k = 1; k <= ktf - 1; ++k) {
        real dp = R(100.0) * (po_cup[k] - po_cup[k + 1]);
        norm = norm + melting_layer[k] * dp / C::g;
      }

      real scale = (R(100.0) * (po_cup[1] - po_cup[ktf])) / C::g;
      for (int k = 1; k <= kte; ++k) {
        melting_layer[k] = melting_layer[k] / (norm + R(1.e-6)) * scale;
      }
    }
  } else {
    for (int k = 1; k <= kte; ++k) {
      p_liq_ice[k] = R(1.0);
      melting_layer[k] = R(0.0);
    }
  }
}

// ---------------------------------------------------------------------------
// get_melting_profile (lines 5917-5991 of module_cu_gfl_deep.F90)
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void get_melting_profile(
    int ierr, const real* tn_cup, const real* po_cup, const real* p_liq_ice,
    const real* melting_layer, const real* qrco, const real* pwo, real edto,
    const real* pwdo, real* melting, int ktf, int kte, bool is_deep) {
  using C = GflConstants<real>;

  if (C::melt_glac && is_deep) {
    if (ierr > 0) {
      for (int k = 1; k <= kte; ++k) melting[k] = R(0.0);
    }
    real total_pwo_solid_phase = R(0.0);
    if (ierr == 0) {
      for (int k = 1; k <= ktf - 1; ++k) {
        real dp = R(100.0) * (po_cup[k] - po_cup[k + 1]);
        real pwo_eff =
            R(0.5) * (pwo[k] + pwo[k + 1] + edto * (pwdo[k] + pwdo[k + 1]));
        real pwo_solid_phase = (R(1.0) - p_liq_ice[k]) * pwo_eff;
        total_pwo_solid_phase =
            total_pwo_solid_phase + pwo_solid_phase * dp / C::g;
      }
      real denom = (R(100.0) * (po_cup[1] - po_cup[ktf])) / C::g;
      for (int k = 1; k <= ktf; ++k) {
        melting[k] = melting_layer[k] * (total_pwo_solid_phase / denom);
      }
    }
  } else {
    for (int k = 1; k <= kte; ++k) melting[k] = R(0.0);
  }
  (void)tn_cup;
  (void)qrco;
}

// ---------------------------------------------------------------------------
// cup_dd_moisture (lines 3035-3185 of module_cu_gfl_deep.F90)
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void cup_dd_moisture(
    const real* zd, const real* hcd, const real* hes_cup, real* qcd,
    const real* qes_cup, real* pwd, const real* q_cup, const real* z_cup,
    const real* dd_massentr, const real* dd_massdetr, int jmin, int& ierr,
    const real* gamma_cup, real& pwev, real& bu, real* qrcd, const real* p_cup,
    const real* q, const real* he, int iloop, int ktf) {
  using C = GflConstants<real>;

  bu = R(0.0);
  pwev = R(0.0);
  for (int k = 1; k <= ktf; ++k) {
    qcd[k] = R(0.0);
    qrcd[k] = R(0.0);
    pwd[k] = R(0.0);
  }
  if (ierr != 0) return;

  int k = jmin;
  real dz = z_cup[k + 1] - z_cup[k];
  real dp = -R(100.0) * (p_cup[k + 1] - p_cup[k]);
  qcd[k] = q_cup[k];
  real dh = hcd[k] - hes_cup[k];
  if (dh < R(0.0)) {
    qrcd[k] = qes_cup[k] + (R(1.0) / C::xlv) *
                               (gamma_cup[k] / (R(1.0) + gamma_cup[k])) * dh;
  } else {
    qrcd[k] = qes_cup[k];
  }
  pwd[jmin] = zd[jmin] * min(R(0.0), qcd[k] - qrcd[k]);
  qcd[k] = qrcd[k];
  pwev = pwev + pwd[jmin] * C::g / dp;
  bu = dz * dh;

  for (int ki = jmin - 1; ki >= 1; --ki) {
    dz = z_cup[ki + 1] - z_cup[ki];
    dp = -R(100.0) * (p_cup[ki + 1] - p_cup[ki]);
    real denom = zd[ki + 1] - R(0.5) * dd_massdetr[ki] + dd_massentr[ki];
    if (denom < R(1.e-16)) {
      ierr = 51;
      break;
    }
    qcd[ki] = (qcd[ki + 1] * zd[ki + 1] -
               R(0.5) * dd_massdetr[ki] * qcd[ki + 1] + dd_massentr[ki] * q[ki]) /
              (zd[ki + 1] - R(0.5) * dd_massdetr[ki] + dd_massentr[ki]);
    dh = hcd[ki] - hes_cup[ki];
    bu = bu + dz * dh;
    qrcd[ki] = qes_cup[ki] + (R(1.0) / C::xlv) *
                                 (gamma_cup[ki] / (R(1.0) + gamma_cup[ki])) *
                                 dh;
    real dqeva = qcd[ki] - qrcd[ki];
    if (dqeva > R(0.0)) {
      dqeva = R(0.0);
      qrcd[ki] = qcd[ki];
    }
    pwd[ki] = zd[ki] * dqeva;
    qcd[ki] = qrcd[ki];
    pwev = pwev + pwd[ki] * C::g / dp;
  }
  if (pwev == R(0.0) && iloop == 1) ierr = 7;
  if (bu >= R(0.0) && iloop == 1) ierr = 7;
  (void)he;
}

// ---------------------------------------------------------------------------
// cup_dd_edt (lines 2890-3032 of module_cu_gfl_deep.F90)
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void cup_dd_edt(
    int ierr, const real* us, const real* vs, const real* z, int ktop,
    int kbcon, real& edt, const real* p, real pwav, const real* pw, real ccn,
    real ccnclean, real pwev, real edtmax, real edtmin, real& edtc, real psum2,
    real psumh, const real* rho, int aeroevap, real& pefc, int xland1,
    int ktf) {
  pefc = R(0.0);
  edt = R(0.0);
  edtc = R(0.0);
  real vws = R(0.0), sdp = R(0.0), vshear = R(0.0);

  for (int kk = 1; kk <= ktf - 1; ++kk) {
    if (ierr != 0) continue;
    if (kk <= min(ktop, ktf) && kk >= kbcon) {
      vws = vws + (abs((us[kk + 1] - us[kk]) / (z[kk + 1] - z[kk])) +
                   abs((vs[kk + 1] - vs[kk]) / (z[kk + 1] - z[kk]))) *
                      (p[kk] - p[kk + 1]);
      sdp = sdp + p[kk] - p[kk + 1];
    }
    if (kk == ktf - 1) vshear = R(1.e3) * vws / sdp;
  }

  if (ierr == 0) {
    real pef = R(1.591) - R(0.639) * vshear + R(0.0953) * (vshear * vshear) -
               R(0.00496) * (vshear * vshear * vshear);
    if (pef > R(0.9)) pef = R(0.9);
    if (pef < R(0.1)) pef = R(0.1);

    real zkbc = z[kbcon] * R(3.281e-3);
    real prezk = R(0.02);
    if (zkbc > R(3.0)) {
      prezk = R(0.96729352) +
              zkbc * (R(-0.70034167) +
                      zkbc * (R(0.162179896) +
                              zkbc * (R(-1.2569798e-2) +
                                      zkbc * (R(4.2772e-4) - zkbc * R(5.44e-6)))));
    }
    if (zkbc > R(25.0)) prezk = R(2.4);
    real pefb = R(1.0) / (R(1.0) + prezk);
    if (pefb > R(0.9)) pefb = R(0.9);
    if (pefb < R(0.1)) pefb = R(0.1);
    pefb = pef;

    edt = R(1.0) - R(0.5) * (pefb + pef);
    if (aeroevap > 1) {
      constexpr real alpha3 = R(0.75);
      constexpr real beta3 = R(-0.15);
      pefb = (xland1 == 1) ? R(0.3) : R(0.5);
      if (psumh > R(0.0) && psum2 > R(0.0)) {
        real aeroadd = pow(ccnclean, beta3) * pow(psumh, alpha3 - R(1.0));
        real prop_c = pefb / aeroadd;
        aeroadd = pow(ccn, beta3) * pow(psum2, alpha3 - R(1.0));
        aeroadd = prop_c * aeroadd;
        pefc = aeroadd;
        if (pefc > R(0.9)) pefc = R(0.9);
        if (pefc < R(0.1)) pefc = R(0.1);
        edt = R(1.0) - pefc;
      }
    }
    edtc = edt;
  }

  if (ierr == 0) {
    edtc = -edtc * psum2 / pwev;
    if (edtc > edtmax) edtc = edtmax;
    if (edtc < edtmin) edtc = edtmin;
  }
  (void)pwav;
  (void)pw;
  (void)rho;
}

// ---------------------------------------------------------------------------
// cup_forcing_ens_3d (lines 3433-3839 of module_cu_gfl_deep.F90)
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void cup_forcing_ens_3d(
    real& closure_n, int xland, real aa0, real aa1, real xaa0, real mbdt,
    real dtime, int ierr, int ierr2, int ierr3, real* xf_ens, real axx,
    real mconv, const real* p_cup, int ktop, const real* omeg, const real* zd,
    const real* zu, const real* pr_ens, real edt, int kbcon, int ichoice,
    int dicycle, real tau_ecmwf, real aa1_bl, real& xf_dicycle) {
  using C = GflConstants<real>;
  constexpr int maxens3 = C::maxens3;

  real ens_adj = R(1.0);
  real xff_dicycle = R(0.0);
  real xff_ens3[maxens3 + 1] = {};

  if (ierr == 0) {
    real a_ave = axx;
    a_ave = max(R(0.0), a_ave);
    a_ave = min(a_ave, aa1);
    a_ave = max(R(0.0), a_ave);
    (void)a_ave;

    real xff0 = (aa1 - aa0) / dtime;
    xff_ens3[1] = max(R(0.0), (aa1 - aa0) / dtime);
    xff_ens3[2] = max(R(0.0), (aa1 - aa0) / dtime);
    xff_ens3[3] = max(R(0.0), (aa1 - aa0) / dtime);
    xff_ens3[16] = max(R(0.0), (aa1 - aa0) / dtime);

    real xomg = R(0.0);
    int kk = 0;
    xff_ens3[4] = R(0.0);
    xff_ens3[5] = R(0.0);
    xff_ens3[6] = R(0.0);
    for (int k = kbcon - 1; k <= kbcon + 1; ++k) {
      if (zu[k] > R(0.0)) {
        xomg = xomg - omeg[k] / R(9.81);
        kk = kk + 1;
      }
    }
    if (kk > 0) xff_ens3[4] = xomg / real(kk);
    xff_ens3[4] = C::betajb * xff_ens3[4];
    xff_ens3[5] = xff_ens3[4];
    xff_ens3[6] = xff_ens3[4];
    if (xff_ens3[4] < R(0.0)) xff_ens3[4] = R(0.0);
    if (xff_ens3[5] < R(0.0)) xff_ens3[5] = R(0.0);
    if (xff_ens3[6] < R(0.0)) xff_ens3[6] = R(0.0);
    xff_ens3[14] = xff_ens3[4];

    xff_ens3[7] = mconv;
    xff_ens3[8] = mconv;
    xff_ens3[9] = mconv;
    xff_ens3[15] = mconv;

    xff_ens3[10] = aa1 / tau_ecmwf;
    xff_ens3[11] = aa1 / tau_ecmwf;
    xff_ens3[12] = aa1 / tau_ecmwf;
    xff_ens3[13] = aa1 / tau_ecmwf;

    if (ichoice == 0) {
      if (xff0 < R(0.0)) {
        for (int n = 1; n <= maxens3; ++n) xff_ens3[n] = R(0.0);
        xff_ens3[1] = R(0.0);
        xff_ens3[2] = R(0.0);
        xff_ens3[3] = R(0.0);
        xff_ens3[10] = R(0.0);
        xff_ens3[11] = R(0.0);
        xff_ens3[12] = R(0.0);
        xff_ens3[13] = R(0.0);
        xff_ens3[16] = R(0.0);
      }
    }

    real xk = (xaa0 - aa1) / mbdt;
    if (xk < R(0.0) && xk > R(-0.01) * mbdt) xk = R(-0.01) * mbdt;
    if (xk >= R(0.0) && xk < R(1.e-2)) xk = R(1.e-2);

    if (xland < 0.1) {
      if (ierr2 > 0 || ierr3 > 0) {
        for (int n = 1; n <= 16; ++n) xff_ens3[n] = ens_adj * xff_ens3[n];
      }
    }

    if (xk < R(0.0)) {
      if (xff_ens3[1] > R(0.0)) xf_ens[1] = max(R(0.0), -xff_ens3[1] / xk);
      if (xff_ens3[2] > R(0.0)) xf_ens[2] = max(R(0.0), -xff_ens3[2] / xk);
      if (xff_ens3[3] > R(0.0)) xf_ens[3] = max(R(0.0), -xff_ens3[3] / xk);
      if (xff_ens3[16] > R(0.0)) xf_ens[16] = max(R(0.0), -xff_ens3[16] / xk);
    } else {
      xff_ens3[1] = R(0.0);
      xff_ens3[2] = R(0.0);
      xff_ens3[3] = R(0.0);
      xff_ens3[16] = R(0.0);
    }

    xf_ens[4] = max(R(0.0), xff_ens3[4]);
    xf_ens[5] = max(R(0.0), xff_ens3[5]);
    xf_ens[6] = max(R(0.0), xff_ens3[6]);
    xf_ens[14] = max(R(0.0), xff_ens3[14]);
    real a1 = max(R(1.e-3), pr_ens[7]);
    xf_ens[7] = max(R(0.0), xff_ens3[7] / a1);
    a1 = max(R(1.e-3), pr_ens[8]);
    xf_ens[8] = max(R(0.0), xff_ens3[8] / a1);
    a1 = max(R(1.e-3), pr_ens[9]);
    xf_ens[9] = max(R(0.0), xff_ens3[9] / a1);
    a1 = max(R(1.e-3), pr_ens[15]);
    xf_ens[15] = max(R(0.0), xff_ens3[15] / a1);

    if (xk < R(0.0)) {
      xf_ens[10] = max(R(0.0), -xff_ens3[10] / xk);
      xf_ens[11] = max(R(0.0), -xff_ens3[11] / xk);
      xf_ens[12] = max(R(0.0), -xff_ens3[12] / xk);
      xf_ens[13] = max(R(0.0), -xff_ens3[13] / xk);
    } else {
      xf_ens[10] = R(0.0);
      xf_ens[11] = R(0.0);
      xf_ens[12] = R(0.0);
      xf_ens[13] = R(0.0);
    }

    if (ichoice >= 1) {
      real v = xf_ens[ichoice];
      for (int n = 1; n <= 16; ++n) xf_ens[n] = v;
    }
  } else if (ierr != 20 && ierr != 0) {
    for (int n = 1; n <= maxens3; ++n) xf_ens[n] = R(0.0);
  }

  if (dicycle == 1) {
    xf_dicycle = R(0.0);
    if (ierr == 0) {
      real xk = (xaa0 - aa1) / mbdt;
      if (xk < R(0.0) && xk > R(-0.01) * mbdt) xk = R(-0.01) * mbdt;
      if (xk >= R(0.0) && xk < R(1.e-2)) xk = R(1.e-2);
      xff_dicycle = (aa1 - aa1_bl) / tau_ecmwf;
      if (xk < R(0.0)) xf_dicycle = max(R(0.0), -xff_dicycle / xk);
      xf_dicycle = xf_ens[10] - xf_dicycle;
    }
  } else {
    xf_dicycle = R(0.0);
  }
  (void)p_cup;
  (void)ktop;
  (void)zd;
  (void)edt;
}

// ---------------------------------------------------------------------------
// cup_output_ens_3d (lines 4343-4582 of module_cu_gfl_deep.F90)
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void cup_output_ens_3d(
    const real* xff_mid, real* xf_ens, int& ierr, const real* dellat,
    const real* dellaq, const real* dellaqc, real* outtem, real* outq,
    real* outqc, real* subten_h, real* subten_q, real* subten_t, real& pre,
    const real* pw, real& xmb, int ktop, const real* pr_ens, real sig,
    real closure_n, real xmbs_in, int ichoice, int imid, int ktf, int kte,
    int dicycle, real xf_dicycle) {
  using C = GflConstants<real>;
  constexpr int maxens3 = C::maxens3;

  for (int k = 1; k <= kte; ++k) {
    outtem[k] = R(0.0);
    outq[k] = R(0.0);
    outqc[k] = R(0.0);
  }
  pre = R(0.0);
  xmb = R(0.0);

  if (ierr == 0) {
    for (int n = 1; n <= maxens3; ++n) {
      if (pr_ens[n] <= R(0.0)) xf_ens[n] = R(0.0);
    }
  }

  if (imid == 0) {
    if (ierr == 0) {
      int k = 0;
      real xmb_ave = R(0.0);
      for (int n = 1; n <= maxens3; ++n) {
        k = k + 1;
        xmb_ave = xmb_ave + xf_ens[n];
      }
      xmb_ave = xmb_ave / real(k);
      if (dicycle == 2) {
        xmb_ave = xmb_ave - max(R(0.0), xmbs_in);
        xmb_ave = max(R(0.0), xmb_ave);
      } else if (dicycle == 1) {
        xmb_ave = xmb_ave - xf_dicycle;
        xmb_ave = max(R(0.0), xmb_ave);
      }
      real clos_wei = R(16.0) / max(R(1.0), closure_n);
      xmb_ave = min(xmb_ave, R(100.0));
      xmb = clos_wei * sig * xmb_ave;
      if (xmb < R(1.e-16)) ierr = 19;
    }
  } else {
    real xmb_ave = R(0.0);
    if (ierr == 0) {
      if (ichoice == 1 || ichoice == 2) {
        xmb_ave = sig * xff_mid[ichoice];
      } else if (ichoice > 2) {
        int k = 0;
        for (int n = 1; n <= maxens3; ++n) {
          k = k + 1;
          xmb_ave = xmb_ave + xf_ens[n];
        }
        xmb_ave = xmb_ave / real(k);
      } else if (ichoice == 0) {
        xmb_ave = R(0.5) * sig * (xff_mid[1] + xff_mid[2]);
      }
      if (dicycle == 2) {
        xmb = max(R(0.0), xmb_ave - xmbs_in);
      } else if (dicycle == 1) {
        xmb = xmb_ave - xf_dicycle;
        xmb = max(R(0.0), xmb_ave);
      } else if (dicycle == 0) {
        xmb = max(R(0.0), xmb_ave);
      }
    }
  }

  if (ierr == 0) {
    real dtpw = R(0.0);
    for (int k = 1; k <= ktop; ++k) {
      dtpw = dtpw + pw[k];
      outtem[k] = xmb * (dellat[k] - subten_t[k]);
      outq[k] = xmb * (dellaq[k] - subten_q[k]);
      subten_t[k] = subten_t[k] * xmb;
      subten_q[k] = subten_q[k] * xmb;
      subten_h[k] = subten_h[k] * xmb;
      outqc[k] = xmb * dellaqc[k];
    }
    pre = pre + xmb * dtpw;
  } else {
    for (int k = 1; k <= ktf; ++k) {
      outtem[k] = R(0.0);
      outq[k] = R(0.0);
      subten_t[k] = R(0.0);
      subten_q[k] = R(0.0);
      subten_h[k] = R(0.0);
      outqc[k] = R(0.0);
    }
    pre = R(0.0);
  }
}

// ---------------------------------------------------------------------------
// neg_check (lines 4202-4339 of module_cu_gfl_deep.F90)
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void neg_check(bool is_shallow_or_mid, real dt,
                                         const real* q, real* outq, real* sq,
                                         real* outt, real* st, real* outu,
                                         real* su, real* outv, real* sv,
                                         real* outqc, real& pret, int ktop) {
  if (ktop <= 2) return;

  real thresh = R(800.01);
  real names = R(1.0);
  if (is_shallow_or_mid) {
    thresh = R(148.01);
    names = R(1.0);
  }

  real qmemf = R(1.0);
  for (int k = 1; k <= ktop; ++k) {
    real qmem = (outt[k] + st[k]) * R(86400.0);
    if (qmem > thresh) {
      real qmem2 = thresh / qmem;
      qmemf = min(qmemf, qmem2);
    }
    if (qmem < R(-0.5) * thresh * names) {
      real qmem2 = R(-0.5) * names * thresh / qmem;
      qmemf = min(qmemf, qmem2);
    }
  }
  for (int k = 1; k <= ktop; ++k) {
    outq[k] = outq[k] * qmemf;
    sq[k] = sq[k] * qmemf;
    outt[k] = outt[k] * qmemf;
    st[k] = st[k] * qmemf;
    outu[k] = outu[k] * qmemf;
    su[k] = su[k] * qmemf;
    outv[k] = outv[k] * qmemf;
    sv[k] = sv[k] * qmemf;
    outqc[k] = outqc[k] * qmemf;
  }
  pret = pret * qmemf;

  thresh = R(1.e-32);
  qmemf = R(1.0);
  for (int k = 1; k <= ktop; ++k) {
    real qmem = outq[k] + sq[k];
    if (abs(qmem) > R(0.0) && q[k] > R(1.e-6)) {
      real qtest = q[k] + (outq[k] + sq[k]) * dt;
      if (qtest < thresh) {
        real qmem1 = abs(outq[k] + sq[k]);
        real qmem2 = abs((thresh - q[k]) / dt);
        qmemf = min(qmemf, qmem2 / qmem1);
        qmemf = max(R(0.0), qmemf);
      }
    }
  }
  for (int k = 1; k <= ktop; ++k) {
    outq[k] = outq[k] * qmemf;
    sq[k] = sq[k] * qmemf;
    outt[k] = outt[k] * qmemf;
    st[k] = st[k] * qmemf;
    outu[k] = outu[k] * qmemf;
    su[k] = su[k] * qmemf;
    outv[k] = outv[k] * qmemf;
    sv[k] = sv[k] * qmemf;
    outqc[k] = outqc[k] * qmemf;
  }
  pret = pret * qmemf;
}

template <typename real>
struct GflColumnIn {
  int nz;
  const real *t, *q, *tn, *qo, *po, *us, *vs, *rho, *omeg, *zo, *dhdt;
  real z1, psur, hfx, qfx, xland, dx, dtime;
  int kpbl, ichoice;
  real xmbs_in;
};

template <typename real>
struct GflDeepOut {
  int ierr, kbcon, ktop, k22, jmin;
  real frh_out;
};

// ---------------------------------------------------------------------------
// cu_gfl_deep_run (lines 176-2610 of module_cu_gfl_deep.F90)
// Called from cu_grell_freitas_li with imid = 0, sub3d = 0, dicycle = 0,
// nranflag = 0, do_capsuppress = 0, do_smoke_transport = .false.
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline GflDeepOut<real> cu_gfl_deep_run(
    const GflColumnIn<real>& in, real* outu, real* outv, real* outt, real* outq,
    real* outqc, real* subten_h, real* subten_q, real* subten_t, real* subten_u,
    real* subten_v, real* cupclw, real* cnvwt, real& pre, real& xmb_out) {
  using C = GflConstants<real>;
  const int kte = in.nz, ktf = in.nz;
  const real* t = in.t;
  const real* q = in.q;
  const real* tn = in.tn;
  const real* qo = in.qo;
  const real* po = in.po;
  const real* us = in.us;
  const real* vs = in.vs;
  const real* rho = in.rho;
  const real* omeg = in.omeg;
  const real z1 = in.z1, psur = in.psur, dx = in.dx, dtime = in.dtime;
  const int kpbl = in.kpbl;

  constexpr int imid = 0;
  constexpr int sub3d = 0;
  constexpr int dicycle = 0;
  int csum = 0;
  int ierr = 0, kbcon = 0, ktop = 0, k22 = 0, jmin = 0, kdet = 0, kbmax = 1;
  int kstabi = 0, kzdown = 0, ktopdby = 0, pmin_lev = 0;
  int ierr2 = 0, ierr3 = 0;

  real zo[KM] = {}, z[KM] = {}, xz[KM] = {};
  real zuo[KM] = {}, zdo[KM] = {}, zu[KM] = {}, xzu[KM] = {};
  real qes[KM] = {}, he[KM] = {}, hes[KM] = {};
  real qeso[KM] = {}, heo[KM] = {}, heso[KM] = {};
  real qes_cup[KM] = {}, q_cup[KM] = {}, he_cup[KM] = {}, hes_cup[KM] = {};
  real z_cup[KM] = {}, p_cup[KM] = {}, gamma_cup[KM] = {}, t_cup[KM] = {};
  real qeso_cup[KM] = {}, qo_cup[KM] = {}, heo_cup[KM] = {}, heso_cup[KM] = {};
  real zo_cup[KM] = {}, po_cup[KM] = {}, gammao_cup[KM] = {}, tn_cup[KM] = {};
  real u_cup[KM] = {}, v_cup[KM] = {}, uc[KM] = {}, vc[KM] = {};
  real ucd[KM] = {}, vcd[KM] = {};
  real hc[KM] = {}, dby[KM] = {}, hco[KM] = {}, dbyo[KM] = {}, dbyt[KM] = {};
  real cd[KM] = {}, cdd[KM] = {}, c1d[KM] = {}, c0t3d[KM] = {};
  real entr_rate_2d[KM] = {}, mentrd_rate_2d[KM] = {};
  real up_massentr[KM] = {}, up_massdetr[KM] = {};
  real up_massentro[KM] = {}, up_massdetro[KM] = {};
  real up_massentru[KM] = {}, up_massdetru[KM] = {};
  real dd_massentro[KM] = {}, dd_massdetro[KM] = {};
  real dd_massentru[KM] = {}, dd_massdetru[KM] = {};
  real hcdo[KM] = {}, qcdo[KM] = {}, qrcdo[KM] = {}, pwdo[KM] = {},
       dbydo[KM] = {};
  real qco[KM] = {}, qrco[KM] = {}, pwo[KM] = {}, clw_all[KM] = {};
  real dellu[KM] = {}, dellv[KM] = {}, dellah[KM] = {}, dellaq[KM] = {};
  real dellat[KM] = {}, dellaqc[KM] = {};
  real dellat_ens[KM] = {}, dellaq_ens[KM] = {}, dellaqc_ens[KM] = {},
       pwo_ens[KM] = {};
  real xhe[KM] = {}, xq[KM] = {}, xt[KM] = {}, xqes[KM] = {}, xhes[KM] = {};
  real xqes_cup[KM] = {}, xq_cup[KM] = {}, xhe_cup[KM] = {}, xhes_cup[KM] = {};
  real xz_cup[KM] = {}, xt_cup[KM] = {};
  real xhc[KM] = {}, xdby[KM] = {};
  real p_liq_ice[KM] = {}, melting_layer[KM] = {}, melting[KM] = {};
  real xf_ens[C::maxens3 + 1] = {}, pr_ens[C::maxens3 + 1] = {};
  real xff_mid[3] = {};

  for (int k = 1; k <= kte; ++k) zo[k] = in.zo[k];

  const real flux_tun = C::fluxtune;
  const real pmin = R(150.0);
  const real elocp = C::xlv / C::cp;
  const real el2orc = C::xlv * C::xlv / (C::r_v * C::cp);
  const real evfact = R(0.25);
  const real evfactl = R(0.25);
  const real pgcon = R(0.0);
  real lambau = R(2.0);

  real c0 = R(0.004);
  int xland1 = int(in.xland + R(0.0001));
  if (in.xland > R(1.5) || in.xland < R(0.5)) {
    xland1 = 0;
  }
  if (xland1 == 1) c0 = R(0.002);

  real ztexec = R(0.0), zqexec = R(0.0), zws = R(0.0);
  {
    real buo_flux =
        (in.hfx / C::cp + R(0.608) * t[1] * in.qfx / C::xlv) / rho[1];
    real pgeoh = zo[2] * C::g;
    zws = max(R(0.0), flux_tun * R(0.41) * buo_flux * zo[2] * C::g / t[1]);
    if (zws > std::numeric_limits<real>::min()) {
      (void)pgeoh;
      zws = R(1.2) * pow(zws, R(0.3333));
      ztexec = max(flux_tun * in.hfx / (rho[1] * zws * C::cp), R(0.0));
      zqexec = max(flux_tun * in.qfx / C::xlv / (rho[1] * zws), R(0.0));
    }
    zws = max(R(0.0),
              R(0.001) - flux_tun * R(0.41) * buo_flux * zo[kpbl] * C::g / t[kpbl]);
    zws = R(1.2) * pow(zws, R(0.3333));
    zws = zws * rho[kpbl];
  }

  const real cap_maxs = R(75.0);
  real edto = R(0.0);
  real closure_n = R(16.0);
  xmb_out = R(0.0);
  real cap_max = cap_maxs;
  real cap_max_increment = R(20.0);
  if (xland1 != 0) {
    if (ztexec > R(0.0)) cap_max = cap_max + R(25.0);
    if (ztexec < R(0.0)) cap_max = cap_max - R(25.0);
  }
  if (C::use_excess == 0) {
    ztexec = R(0.0);
    zqexec = R(0.0);
  }

  int start_level = kte;
  real entr_rate = C::entr_rate_deep;
  if (dx < C::dx_thresh) entr_rate = C::entr_rate_deep_hr;
  real radius = R(0.2) / entr_rate;
  real frh = min(R(1.0), R(3.14) * radius * radius / dx / dx);
  if (frh > C::frh_thresh) {
    if (sub3d == 0) {
      frh = C::frh_thresh;
      radius = sqrt(frh * dx * dx / R(3.14));
      entr_rate = R(0.2) / radius;
    }
  }
  real sig = (R(1.0) - frh) * (R(1.0) - frh);
  real frh_out = frh;
  real sig_init = sig;
  const real sig_thresh = (R(1.0) - C::frh_thresh) * (R(1.0) - C::frh_thresh);

  for (int k = 1; k <= ktf; ++k) {
    cnvwt[k] = R(0.0);
    zuo[k] = R(0.0);
    zdo[k] = R(0.0);
    z[k] = zo[k];
    xz[k] = zo[k];
    cupclw[k] = R(0.0);
    cd[k] = R(0.1) * entr_rate;
    cdd[k] = R(1.e-9);
    hcdo[k] = R(0.0);
    qrcdo[k] = R(0.0);
    dellaqc[k] = R(0.0);
  }

  const real edtmax = R(1.0);
  const real edtmin = R(0.1);
  const real depth_min = R(1000.0);
  real aa0 = R(0.0), aa1 = R(0.0), edt = R(0.0);
  const int kstabm = ktf - 1;
  const real zkbmax = R(4000.0), zcutdown = R(4000.0), z_detr = R(500.0);

  cup_env(ktf, kte, ierr, -1, psur, z1, po, t, z, he, hes, qes, q);
  cup_env(ktf, kte, ierr, -1, psur, z1, po, tn, zo, heo, heso, qeso, qo);
  cup_env_clev(ktf, kte, ierr, psur, z1, t, qes, q, he, hes, z, po, qes_cup,
               q_cup, he_cup, hes_cup, z_cup, p_cup, gamma_cup, t_cup);
  cup_env_clev(ktf, kte, ierr, psur, z1, tn, qeso, qo, heo, heso, zo, po,
               qeso_cup, qo_cup, heo_cup, heso_cup, zo_cup, po_cup, gammao_cup,
               tn_cup);

  get_partition_liq_ice(ierr, tn, po_cup, p_liq_ice, melting_layer, ktf, kte,
                        true);

  if (ierr == 0) {
    u_cup[1] = us[1];
    v_cup[1] = vs[1];
    for (int k = 2; k <= ktf; ++k) {
      u_cup[k] = R(0.5) * (us[k - 1] + us[k]);
      v_cup[k] = R(0.5) * (vs[k - 1] + vs[k]);
    }
    for (int k = 1; k <= ktf; ++k) {
      if (zo_cup[k] > zkbmax + z1) {
        kbmax = k;
        break;
      }
    }
    for (int k = 1; k <= ktf; ++k) {
      if (zo_cup[k] > z_detr + z1) {
        kdet = k;
        break;
      }
    }
  }

  const int start_k22 = 2;
  if (ierr == 0) {
    k22 = maxloc_abs(heo_cup, start_k22, kbmax + 2);
    if (k22 >= kbmax) {
      ierr = 2;
      ktop = 0;
      k22 = 0;
      kbcon = 0;
    }
  }

  real hkb = R(0.0), hkbo = R(0.0);
  if (ierr == 0) {
    real x_add = C::xlv * zqexec + C::cp * ztexec;
    get_cloud_bc(kte, he_cup, hkb, k22, x_add);
    get_cloud_bc(kte, heo_cup, hkbo, k22, x_add);
  }

  cup_kbcon(cap_max_increment, 1, k22, kbcon, heo_cup, heso_cup, hkbo, ierr,
            kbmax, po_cup, cap_max, ztexec, zqexec, kte, z_cup, entr_rate, heo,
            imid);
  cup_minimi(heso_cup, kbcon, kstabm, ierr, kstabi);

  if (ierr == 0) {
    frh = min(qo_cup[kbcon] / qeso_cup[kbcon], R(1.0));
    if (!(frh >= C::rh_thresh && sig <= sig_thresh)) {
      real x_add0 = R(0.0);
      for (int k = kbcon + 1; k <= ktf; ++k) {
        if (po[kbcon] - po[k] > pmin + x_add0) {
          pmin_lev = k;
          break;
        }
      }
      start_level = k22;
      real x_add = C::xlv * zqexec + C::cp * ztexec;
      get_cloud_bc(kte, he_cup, hkb, k22, x_add);
    }
  }

  if (kstabi < kbcon) {
    kbcon = 1;
    ierr = 1142;
  }
  for (int k = 1; k <= ktf; ++k) entr_rate_2d[k] = entr_rate;
  if (ierr == 0) {
    kbcon = max(2, kbcon);
    for (int k = 2; k <= ktf; ++k) {
      frh = min(qo_cup[k] / qeso_cup[k], R(1.0));
      entr_rate_2d[k] = entr_rate * (R(1.3) - frh);
    }
  }

  rates_up_pdf(1, ktf, kte, ierr, po_cup, ktop, k22, zo_cup, entr_rate_2d, heo,
               heso_cup, hkbo, xland1, kbcon, zuo, ktopdby, csum, pmin_lev);

  if (ierr == 0) {
    if (k22 > 1) {
      for (int k = 1; k <= k22 - 1; ++k) {
        zuo[k] = R(0.0);
        zu[k] = R(0.0);
        xzu[k] = R(0.0);
      }
    }
    real tot_clw = R(0.0);
    for (int k = k22; k <= ktop; ++k) {
      xzu[k] = zuo[k];
      zu[k] = zuo[k];
    }
    if (tot_clw <= R(0.0) && dx < C::dx_thresh) sig = R(1.0);
    sig_init = sig;
    for (int k = ktop + 1; k <= kte; ++k) {
      zuo[k] = R(0.0);
      zu[k] = R(0.0);
      xzu[k] = R(0.0);
    }
  }

  get_lateral_massflux(ktf, kte, ierr, ktop, zo_cup, zuo, cd, entr_rate_2d,
                       up_massentro, up_massdetro, up_massentr, up_massdetr, 1,
                       kbcon, k22, up_massentru, up_massdetru, lambau);

  for (int k = 1; k <= ktf; ++k) {
    uc[k] = R(0.0);
    vc[k] = R(0.0);
    hc[k] = R(0.0);
    dby[k] = R(0.0);
    hco[k] = R(0.0);
    dbyo[k] = R(0.0);
  }
  if (ierr == 0) {
    for (int k = 1; k <= start_level; ++k) {
      uc[k] = u_cup[k];
      vc[k] = v_cup[k];
    }
    for (int k = 1; k <= start_level - 1; ++k) {
      hc[k] = he_cup[k];
      hco[k] = heo_cup[k];
    }
    hc[start_level] = hkb;
    hco[start_level] = hkbo;
  }

  // 1st guess for moist static energy and dbyo (not including ice phase)
  int ktopkeep = 0;
  for (int k = 1; k <= kte; ++k) dbyt[k] = R(0.0);
  if (ierr == 0) {
    ktopkeep = ktop;
    for (int k = start_level + 1; k <= ktop; ++k) {
      real denom = zuo[k - 1] - R(0.5) * up_massdetro[k - 1] + up_massentro[k - 1];
      if (denom < R(1.e-8)) {
        ierr = 51;
        break;
      }
      hco[k] = (hco[k - 1] * zuo[k - 1] -
                R(0.5) * up_massdetro[k - 1] * hco[k - 1] +
                up_massentro[k - 1] * heo[k - 1]) /
               (zuo[k - 1] - R(0.5) * up_massdetro[k - 1] + up_massentro[k - 1]);
      dbyo[k] = hco[k] - heso_cup[k];
    }
    for (int k = ktop - 1; k >= kbcon; --k) {
      if (dbyo[k] > R(0.0)) {
        ktopkeep = k + 1;
        break;
      }
    }
  }

  kzdown = 0;
  if (ierr == 0) {
    real zktop = (zo_cup[ktop] - z1) * R(0.6);
    zktop = min(zktop + z1, zcutdown + z1);
    for (int k = 1; k <= ktf; ++k) {
      if (zo_cup[k] > zktop) {
        kzdown = k;
        kzdown = min(kzdown, kstabi - 1);
        break;
      }
    }
  }

  cup_minimi(heso_cup, k22, kzdown, ierr, jmin);
  if (ierr == 0) {
    int jmini = jmin;
    bool keep_going = true;
    while (keep_going) {
      keep_going = false;
      if (jmini - 1 < kdet) kdet = jmini - 1;
      if (jmini >= ktop - 1) jmini = ktop - 2;
      int ki = jmini;
      hcdo[ki] = heso_cup[ki];
      real dh = R(0.0);
      for (int k = ki - 1; k >= 1; --k) {
        hcdo[k] = heso_cup[jmini];
        real dz = zo_cup[k + 1] - zo_cup[k];
        dh = dh + dz * (hcdo[k] - heso_cup[k]);
        if (dh > R(0.0)) {
          jmini = jmini - 1;
          if (jmini > 5) {
            keep_going = true;
          } else {
            ierr = 9;
            break;
          }
        }
      }
    }
    jmin = jmini;
    if (jmini <= 5) ierr = 4;
  }

  if (ierr == 0) {
    for (int k = ktop + 1; k <= ktf; ++k) {
      hco[k] = heso_cup[k];
      dbyo[k] = R(0.0);
    }
  }

  real pwavo = R(0.0), psum = R(0.0), psumh = R(0.0);
  real ccn = R(1500.0);
  cup_up_moisture(false, ierr, zo_cup, qco, qrco, pwo, pwavo, p_cup, kbcon,
                  ktop, dbyo, clw_all, xland1, qo, gammao_cup, zuo, qeso_cup,
                  k22, qo_cup, c0, c0t3d, zqexec, ccn, C::ccnclean, rho, c1d,
                  tn_cup, C::autoconv, up_massentr, up_massdetr, psum, psumh, 1,
                  ktf, kte);

  ktopkeep = 0;
  for (int k = 1; k <= kte; ++k) dbyt[k] = R(0.0);
  if (ierr == 0) {
    ktopkeep = ktop;
    for (int k = start_level + 1; k <= ktop; ++k) {
      real denom = zuo[k - 1] - R(0.5) * up_massdetro[k - 1] + up_massentro[k - 1];
      if (denom < R(1.e-8)) {
        ierr = 51;
        break;
      }
      hc[k] = (hc[k - 1] * zu[k - 1] -
               R(0.5) * up_massdetr[k - 1] * hc[k - 1] +
               up_massentr[k - 1] * he[k - 1]) /
              (zu[k - 1] - R(0.5) * up_massdetr[k - 1] + up_massentr[k - 1]);
      uc[k] = (uc[k - 1] * zu[k - 1] -
               R(0.5) * up_massdetru[k - 1] * uc[k - 1] +
               up_massentru[k - 1] * us[k - 1] -
               pgcon * R(0.5) * (zu[k] + zu[k - 1]) * (u_cup[k] - u_cup[k - 1])) /
              (zu[k - 1] - R(0.5) * up_massdetru[k - 1] + up_massentru[k - 1]);
      vc[k] = (vc[k - 1] * zu[k - 1] -
               R(0.5) * up_massdetru[k - 1] * vc[k - 1] +
               up_massentru[k - 1] * vs[k - 1] -
               pgcon * R(0.5) * (zu[k] + zu[k - 1]) * (v_cup[k] - v_cup[k - 1])) /
              (zu[k - 1] - R(0.5) * up_massdetru[k - 1] + up_massentru[k - 1]);
      dby[k] = hc[k] - hes_cup[k];
      hco[k] = (hco[k - 1] * zuo[k - 1] -
                R(0.5) * up_massdetro[k - 1] * hco[k - 1] +
                up_massentro[k - 1] * heo[k - 1]) /
               (zuo[k - 1] - R(0.5) * up_massdetro[k - 1] + up_massentro[k - 1]);
      hc[k] = hc[k] + (R(1.0) - p_liq_ice[k]) * qrco[k] * C::xlf;
      hco[k] = hco[k] + (R(1.0) - p_liq_ice[k]) * qrco[k] * C::xlf;
      dby[k] = hc[k] - hes_cup[k];
      dbyo[k] = hco[k] - heso_cup[k];
      real dz = zo_cup[k + 1] - zo_cup[k];
      dbyt[k] = dbyt[k - 1] + dbyo[k] * dz;
    }
    for (int k = ktop - 1; k >= kbcon; --k) {
      if (dbyo[k] > R(0.0)) {
        ktopkeep = k + 1;
        break;
      }
    }
    (void)ktopkeep;
  }

  if (ierr == 0) {
    for (int k = ktop + 1; k <= ktf; ++k) {
      hc[k] = hes_cup[k];
      uc[k] = u_cup[k];
      vc[k] = v_cup[k];
      hco[k] = heso_cup[k];
      dby[k] = R(0.0);
      dbyo[k] = R(0.0);
      zu[k] = R(0.0);
      zuo[k] = R(0.0);
      cd[k] = R(0.0);
      entr_rate_2d[k] = R(0.0);
      up_massentr[k] = R(0.0);
      up_massdetr[k] = R(0.0);
      up_massentro[k] = R(0.0);
      up_massdetro[k] = R(0.0);
    }
  }
  if (ierr == 0) {
    if (ktop < kbcon + 2) {
      ierr = 5;
      ktop = 0;
    }
  }
  if (ierr == 0) {
    if (jmin - 1 < kdet) kdet = jmin - 1;
    if (-zo_cup[kbcon] + zo_cup[ktop] < depth_min) ierr = 6;
  }

  for (int k = 1; k <= ktf; ++k) {
    zdo[k] = R(0.0);
    cdd[k] = R(0.0);
    dd_massentro[k] = R(0.0);
    dd_massdetro[k] = R(0.0);
    dd_massentru[k] = R(0.0);
    dd_massdetru[k] = R(0.0);
    hcdo[k] = heso_cup[k];
    ucd[k] = u_cup[k];
    vcd[k] = v_cup[k];
    dbydo[k] = R(0.0);
    mentrd_rate_2d[k] = entr_rate;
  }

  if (ierr == 0) {
    real beta = max(R(0.025), R(0.055) - real(csum) * R(0.0015));
    real bud = R(0.0);
    for (int k = 1; k <= jmin; ++k) cdd[k] = R(0.1) * entr_rate;
    cdd[jmin] = R(0.0);
    real zd_local[KM] = {};
    real zuh2[40] = {};
    get_zu_zd_pdf_fim(0, po_cup, R(0.0), R(0.0), 0, xland1, zuh2, 4, ierr, kdet,
                      jmin + 1, zd_local, 1, kte, ktf, beta, kpbl, csum,
                      pmin_lev);
    for (int k = 1; k <= jmin; ++k) zdo[k] = zd_local[k];
    bool skip_dd = false;
    if (zdo[jmin] < R(1.e-8)) {
      zdo[jmin] = R(0.0);
      jmin = jmin - 1;
      for (int k = jmin; k <= ktf; ++k) cdd[k] = R(0.0);
      for (int k = jmin + 1; k <= ktf; ++k) zdo[k] = R(0.0);
      if (zdo[jmin] < R(1.e-8)) {
        ierr = 876;
        skip_dd = true;
      }
    }

    if (!skip_dd) {
      int itemp = maxloc_abs(zdo, 1, kte);
      for (int ki = jmin; ki >= itemp; --ki) {
        real dzo = zo_cup[ki + 1] - zo_cup[ki];
        dd_massdetro[ki] = cdd[ki] * dzo * zdo[ki + 1];
        dd_massentro[ki] = zdo[ki] - zdo[ki + 1] + dd_massdetro[ki];
        if (dd_massentro[ki] < R(0.0)) {
          dd_massentro[ki] = R(0.0);
          dd_massdetro[ki] = zdo[ki + 1] - zdo[ki];
          if (zdo[ki + 1] > R(0.0)) {
            cdd[ki] = dd_massdetro[ki] / (dzo * zdo[ki + 1]);
          }
        }
        if (zdo[ki + 1] > R(0.0)) {
          mentrd_rate_2d[ki] = dd_massentro[ki] / (dzo * zdo[ki + 1]);
        }
      }
      mentrd_rate_2d[1] = R(0.0);
      for (int ki = itemp - 1; ki >= 1; --ki) {
        real dzo = zo_cup[ki + 1] - zo_cup[ki];
        dd_massentro[ki] = mentrd_rate_2d[ki] * dzo * zdo[ki + 1];
        dd_massdetro[ki] = zdo[ki + 1] + dd_massentro[ki] - zdo[ki];
        if (dd_massdetro[ki] < R(0.0)) {
          dd_massdetro[ki] = R(0.0);
          dd_massentro[ki] = zdo[ki] - zdo[ki + 1];
          if (zdo[ki + 1] > R(0.0)) {
            mentrd_rate_2d[ki] = dd_massentro[ki] / (dzo * zdo[ki + 1]);
          }
        }
        if (zdo[ki + 1] > R(0.0)) {
          cdd[ki] = dd_massdetro[ki] / (dzo * zdo[ki + 1]);
        }
      }

      for (int k = 2; k <= jmin + 1; ++k) {
        dd_massentru[k - 1] =
            dd_massentro[k - 1] + lambau * dd_massdetro[k - 1];
        dd_massdetru[k - 1] =
            dd_massdetro[k - 1] + lambau * dd_massdetro[k - 1];
      }
      dbydo[jmin] = hcdo[jmin] - heso_cup[jmin];
      bud = dbydo[jmin] * (zo_cup[jmin + 1] - zo_cup[jmin]);
      ucd[jmin + 1] = R(0.5) * (uc[jmin + 1] + u_cup[jmin + 1]);
      for (int ki = jmin; ki >= 1; --ki) {
        real dzo = zo_cup[ki + 1] - zo_cup[ki];
        real h_entr = R(0.5) * (heo[ki] + R(0.5) * (hco[ki] + hco[ki + 1]));
        ucd[ki] = (ucd[ki + 1] * zdo[ki + 1] -
                   R(0.5) * dd_massdetru[ki] * ucd[ki + 1] +
                   dd_massentru[ki] * us[ki] -
                   pgcon * zdo[ki + 1] * (us[ki + 1] - us[ki])) /
                  (zdo[ki + 1] - R(0.5) * dd_massdetru[ki] + dd_massentru[ki]);
        vcd[ki] = (vcd[ki + 1] * zdo[ki + 1] -
                   R(0.5) * dd_massdetru[ki] * vcd[ki + 1] +
                   dd_massentru[ki] * vs[ki] -
                   pgcon * zdo[ki + 1] * (vs[ki + 1] - vs[ki])) /
                  (zdo[ki + 1] - R(0.5) * dd_massdetru[ki] + dd_massentru[ki]);
        hcdo[ki] = (hcdo[ki + 1] * zdo[ki + 1] -
                    R(0.5) * dd_massdetro[ki] * hcdo[ki + 1] +
                    dd_massentro[ki] * h_entr) /
                   (zdo[ki + 1] - R(0.5) * dd_massdetro[ki] + dd_massentro[ki]);
        dbydo[ki] = hcdo[ki] - heso_cup[ki];
        bud = bud + dbydo[ki] * dzo;
      }
      if (bud > R(0.0)) ierr = 7;
    }
  }

  real pwevo = R(0.0), bu = R(0.0);
  cup_dd_moisture(zdo, hcdo, heso_cup, qcdo, qeso_cup, pwdo, qo_cup, zo_cup,
                  dd_massentro, dd_massdetro, jmin, ierr, gammao_cup, pwevo, bu,
                  qrcdo, po_cup, qo, heo, 1, ktf);

  if (ierr == 0) {
    for (int k = 2; k <= ktop; ++k) {
      real dp = R(100.0) * (po_cup[1] - po_cup[2]);
      cupclw[k] = qrco[k];
      cnvwt[k] = zuo[k] * cupclw[k] * C::g / dp;
    }
  }

  cup_up_aa0(ktf, kte, aa0, z, zu, dby, gamma_cup, t_cup, kbcon, ktop, ierr);
  cup_up_aa0(ktf, kte, aa1, zo, zuo, dbyo, gammao_cup, tn_cup, kbcon, ktop,
             ierr);

  if (ierr == 0) {
    if (aa1 == R(0.0)) ierr = 17;
    if (aa1 <= aa0 && dx < C::dx_thresh) ierr = 18;
    if (dx < C::dx_thresh) {
      real trash = (aa1 - aa0) / dtime / frh_out;
      if (trash < R(0.0)) ierr = 2242;
    }
  }

  real aa1_bl = R(0.0), xf_dicycle = R(0.0), tau_ecmwf = R(0.0);
  if (ierr == 0) {
    real wmean = R(3.0);
    tau_ecmwf = (zo_cup[ktop] - zo_cup[kbcon]) / wmean;
    tau_ecmwf = max(tau_ecmwf, R(720.0));
    tau_ecmwf = tau_ecmwf * (R(1.0061) + R(1.23e-2) * (dx / R(1000.0)));
  }

  real axx = aa1;
  real edtc = R(0.0), pefc = R(0.0);
  cup_dd_edt(ierr, us, vs, zo, ktop, kbcon, edt, po, pwavo, pwo, ccn,
             C::ccnclean, pwevo, edtmax, edtmin, edtc, psum, psumh, rho,
             C::aeroevap, pefc, xland1, ktf);
  if (ierr == 0) edto = edtc;

  get_melting_profile(ierr, tn_cup, po_cup, p_liq_ice, melting_layer, qrco, pwo,
                      edto, pwdo, melting, ktf, kte, true);

  for (int k = 1; k <= kte; ++k) {
    dellu[k] = R(0.0);
    dellv[k] = R(0.0);
    dellah[k] = R(0.0);
    dellat[k] = R(0.0);
    dellaq[k] = R(0.0);
    dellaqc[k] = R(0.0);
    subten_h[k] = R(0.0);
    subten_q[k] = R(0.0);
    subten_t[k] = R(0.0);
    subten_u[k] = R(0.0);
    subten_v[k] = R(0.0);
  }

  if (sub3d == 0) {
    if (ierr == 0) {
      real dp = R(100.0) * (po_cup[1] - po_cup[2]);
      dellu[1] =
          C::pgcd * (edto * zdo[2] * ucd[2] - edto * zdo[2] * u_cup[2]) * C::g / dp -
          zuo[2] * (uc[2] - u_cup[2]) * C::g / dp;
      dellv[1] =
          C::pgcd * (edto * zdo[2] * vcd[2] - edto * zdo[2] * v_cup[2]) * C::g / dp -
          zuo[2] * (vc[2] - v_cup[2]) * C::g / dp;
      for (int k = 2; k <= ktop; ++k) {
        dp = R(100.0) * (po_cup[k] - po_cup[k + 1]);
        dellu[k] = -((zuo[k + 1] * (uc[k + 1] - u_cup[k + 1]) -
                      zuo[k] * (uc[k] - u_cup[k])) *
                     C::g / dp) +
                   (zdo[k + 1] * (ucd[k + 1] - u_cup[k + 1]) -
                    zdo[k] * (ucd[k] - u_cup[k])) *
                       C::g / dp * edto * C::pgcd;
        dellv[k] = -((zuo[k + 1] * (vc[k + 1] - v_cup[k + 1]) -
                      zuo[k] * (vc[k] - v_cup[k])) *
                     C::g / dp) +
                   (zdo[k + 1] * (vcd[k + 1] - v_cup[k + 1]) -
                    zdo[k] * (vcd[k] - v_cup[k])) *
                       C::g / dp * edto * C::pgcd;
      }

      dp = R(100.0) * (po_cup[1] - po_cup[2]);
      dellah[1] = (edto * zdo[2] * hcdo[2] - edto * zdo[2] * heo_cup[2]) * C::g / dp -
                  zuo[2] * (hco[2] - heo_cup[2]) * C::g / dp;
      dellaq[1] = (edto * zdo[2] * qcdo[2] - edto * zdo[2] * qo_cup[2]) * C::g / dp -
                  zuo[2] * (qco[2] - qo_cup[2]) * C::g / dp;
      real g_rain = R(0.5) * (pwo[1] + pwo[2]) * C::g / dp;
      real e_dn = -R(0.5) * (pwdo[1] + pwdo[2]) * C::g / dp * edto;
      dellaq[1] = dellaq[1] + e_dn - g_rain;

      for (int k = 2; k <= ktop; ++k) {
        dp = R(100.0) * (po_cup[k] - po_cup[k + 1]);
        dellah[k] = -((zuo[k + 1] * (hco[k + 1] - heo_cup[k + 1]) -
                       zuo[k] * (hco[k] - heo_cup[k])) *
                      C::g / dp) +
                    (zdo[k + 1] * (hcdo[k + 1] - heo_cup[k + 1]) -
                     zdo[k] * (hcdo[k] - heo_cup[k])) *
                        C::g / dp * edto;
        dellah[k] = dellah[k] - C::xlf * melting[k] * C::g / dp;

        real detup = up_massdetro[k];
        real dz = zo_cup[k] - zo_cup[k - 1];
        if (k < ktop) {
          dellaqc[k] = zuo[k] * c1d[k] * qrco[k] * dz / dp * C::g;
        } else {
          dellaqc[k] = detup * R(0.5) * (qrco[k + 1] + qrco[k]) * C::g / dp;
        }

        g_rain = R(0.5) * (pwo[k] + pwo[k + 1]) * C::g / dp;
        e_dn = -R(0.5) * (pwdo[k] + pwdo[k + 1]) * C::g / dp * edto;
        real c_up = dellaqc[k] +
                    (zuo[k + 1] * qrco[k + 1] - zuo[k] * qrco[k]) * C::g / dp +
                    g_rain;
        dellaq[k] = -((zuo[k + 1] * (qco[k + 1] - qo_cup[k + 1]) -
                       zuo[k] * (qco[k] - qo_cup[k])) *
                      C::g / dp) +
                    (zdo[k + 1] * (qcdo[k + 1] - qo_cup[k + 1]) -
                     zdo[k] * (qcdo[k] - qo_cup[k])) *
                        C::g / dp * edto -
                    c_up + e_dn;
      }
    }
  }

  const real mbdt = R(0.1);
  real xaa0_ens = R(0.0);
  if (ierr == 0) {
    for (int k = 1; k <= ktf; ++k) {
      xhe[k] = dellah[k] * mbdt + heo[k];
      xq[k] = max(R(1.e-16), dellaq[k] * mbdt + qo[k]);
      dellat[k] = (R(1.0) / C::cp) * (dellah[k] - C::xlv * dellaq[k]);
      subten_t[k] = (R(1.0) / C::cp) * (subten_h[k] - C::xlv * subten_q[k]);
      xt[k] = dellat[k] * mbdt + tn[k];
      xt[k] = max(R(190.0), xt[k]);
    }
    for (int k = 2; k <= ktf - 1; ++k) {
      xt[k] = tn[k] +
              R(0.25) * (dellat[k - 1] + R(2.0) * dellat[k] + dellat[k + 1]) *
                  mbdt;
      xt[k] = max(R(190.0), xt[k]);
      xq[k] = max(
          R(1.e-16),
          qo[k] + R(0.25) * (dellaq[k - 1] + R(2.0) * dellaq[k] + dellaq[k + 1]) *
                      mbdt);
      xhe[k] = heo[k] +
               R(0.25) * (dellah[k - 1] + R(2.0) * dellah[k] + dellah[k + 1]) *
                   mbdt;
    }
    xhe[ktf] = heo[ktf];
    xq[ktf] = qo[ktf];
    xt[ktf] = tn[ktf];
  }

  cup_env(ktf, kte, ierr, -1, psur, z1, po, xt, xz, xhe, xhes, xqes, xq);
  cup_env_clev(ktf, kte, ierr, psur, z1, xt, xqes, xq, xhe, xhes, xz, po,
               xqes_cup, xq_cup, xhe_cup, xhes_cup, xz_cup, po_cup, gamma_cup,
               xt_cup);

  for (int k = 1; k <= ktf; ++k) {
    xhc[k] = R(0.0);
    xdby[k] = R(0.0);
  }
  if (ierr == 0) {
    real x_add = C::xlv * zqexec + C::cp * ztexec;
    real xhkb = R(0.0);
    get_cloud_bc(kte, xhe_cup, xhkb, k22, x_add);
    for (int k = 1; k <= start_level - 1; ++k) xhc[k] = xhe_cup[k];
    xhc[start_level] = xhkb;
    for (int k = start_level + 1; k <= ktop; ++k) {
      xhc[k] = (xhc[k - 1] * xzu[k - 1] -
                R(0.5) * up_massdetro[k - 1] * xhc[k - 1] +
                up_massentro[k - 1] * xhe[k - 1]) /
               (xzu[k - 1] - R(0.5) * up_massdetro[k - 1] + up_massentro[k - 1]);
      xhc[k] = xhc[k] + C::xlf * (R(1.0) - p_liq_ice[k]) * qrco[k];
      xdby[k] = xhc[k] - xhes_cup[k];
    }
    for (int k = ktop + 1; k <= ktf; ++k) {
      xhc[k] = xhes_cup[k];
      xdby[k] = R(0.0);
    }
  }

  real xaa0 = R(0.0);
  cup_up_aa0(ktf, kte, xaa0, xz, xzu, xdby, gamma_cup, xt_cup, kbcon, ktop,
             ierr);
  if (ierr == 0) {
    xaa0_ens = xaa0;
    for (int k = 1; k <= ktop; ++k) {
      for (int n = 1; n <= C::maxens3; ++n) {
        pr_ens[n] = pr_ens[n] + pwo[k] + edto * pwdo[k];
      }
    }
    if (pr_ens[7] < R(1.e-6)) {
      ierr = 18;
      for (int n = 1; n <= C::maxens3; ++n) pr_ens[n] = R(0.0);
    }
    for (int n = 1; n <= C::maxens3; ++n) {
      if (pr_ens[n] < R(1.e-5)) pr_ens[n] = R(0.0);
    }
  }

  ierr2 = ierr;
  ierr3 = ierr;
  int k22x = k22, kbconx = 1;
  cup_maximi(heo_cup, 2, kbmax, ierr, k22x);
  {
    real hkbo2 = hkbo;
    cup_kbcon(cap_max_increment, 2, k22x, kbconx, heo_cup, heso_cup, hkbo2,
              ierr2, kbmax, po_cup, cap_max, ztexec, zqexec, kte, z_cup,
              entr_rate, heo, imid);
    cup_kbcon(cap_max_increment, 3, k22x, kbconx, heo_cup, heso_cup, hkbo2,
              ierr3, kbmax, po_cup, cap_max, ztexec, zqexec, kte, z_cup,
              entr_rate, heo, imid);
  }

  real mconv = R(0.0);
  if (ierr == 0) {
    for (int k = 1; k <= ktop; ++k) {
      real dq = qo_cup[k + 1] - qo_cup[k];
      mconv = mconv + omeg[k] * dq / C::g;
    }
    if (mconv < R(0.0) && dx < C::dx_thresh) ierr = 2242;
  }

  cup_forcing_ens_3d(closure_n, xland1, aa0, aa1, xaa0_ens, mbdt, dtime, ierr,
                     ierr2, ierr3, xf_ens, axx, mconv, po_cup, ktop, omeg, zdo,
                     zuo, pr_ens, edto, kbcon, in.ichoice, dicycle, tau_ecmwf,
                     aa1_bl, xf_dicycle);

  for (int k = 1; k <= ktf; ++k) {
    if (ierr == 0) {
      dellat_ens[k] = dellat[k];
      dellaq_ens[k] = dellaq[k];
      dellaqc_ens[k] = dellaqc[k];
      pwo_ens[k] = pwo[k] + edto * pwdo[k];
    } else {
      dellat_ens[k] = R(0.0);
      dellaq_ens[k] = R(0.0);
      dellaqc_ens[k] = R(0.0);
      pwo_ens[k] = R(0.0);
    }
  }

  real xmb = R(0.0);
  cup_output_ens_3d(xff_mid, xf_ens, ierr, dellat_ens, dellaq_ens, dellaqc_ens,
                    outt, outq, outqc, subten_h, subten_q, subten_t, pre,
                    pwo_ens, xmb, ktop, pr_ens, sig, closure_n, in.xmbs_in,
                    in.ichoice, imid, ktf, kte, dicycle, xf_dicycle);

  if (ierr == 0 && pre > R(0.0)) {
    pre = max(pre, R(0.0));
    xmb_out = xmb;
    outu[1] = dellu[1] * xmb;
    outv[1] = dellv[1] * xmb;
    for (int k = 2; k <= ktop; ++k) {
      outu[k] = R(0.25) *
                ((dellu[k - 1] - subten_u[k - 1]) +
                 R(2.0) * (dellu[k] - subten_u[k]) +
                 (dellu[k + 1] - subten_u[k + 1])) *
                xmb;
      outv[k] = R(0.25) *
                ((dellv[k - 1] - subten_v[k - 1]) +
                 R(2.0) * (dellv[k] - subten_v[k]) +
                 (dellv[k + 1] - subten_v[k + 1])) *
                xmb;
    }
    for (int k = 1; k <= ktop; ++k) {
      subten_u[k] = subten_u[k] * xmb;
      subten_v[k] = subten_v[k] * xmb;
    }
  } else if (ierr != 0 || pre == R(0.0)) {
    ktop = 0;
    pre = R(0.0);
    xmb_out = R(0.0);
    kbcon = 0;
    k22 = 0;
    for (int k = 1; k <= kte; ++k) {
      outt[k] = R(0.0);
      outq[k] = R(0.0);
      outqc[k] = R(0.0);
      outu[k] = R(0.0);
      outv[k] = R(0.0);
      cupclw[k] = R(0.0);
      cnvwt[k] = R(0.0);
      subten_h[k] = R(0.0);
      subten_q[k] = R(0.0);
      subten_t[k] = R(0.0);
      subten_u[k] = R(0.0);
      subten_v[k] = R(0.0);
    }
  }

  if (C::irainevap == 1) {
    real rntot = R(0.0), delqev = R(0.0), delq2 = R(0.0), rn = R(0.0);
    if (ierr == 0) {
      for (int k = ktop; k >= 1; --k) {
        real rain = pwo[k] + edto * pwdo[k];
        rntot = rntot + rain * xmb * R(0.001) * dtime;
      }
      real qevap = R(0.0);
      bool flg = true;
      real evef = edt * evfact * (sig_init * sig_init);
      if (in.xland > R(0.5) && in.xland < R(1.5)) {
        evef = edt * evfactl * (sig_init * sig_init);
      }
      for (int k = ktop; k >= 1; --k) {
        real rain = pwo[k] + edto * pwdo[k];
        rn = rn + rain * xmb * R(0.001) * dtime;
        if (flg) {
          real q1 = qo[k] + outq[k] * dtime;
          real t1 = tn[k] + outt[k] * dtime;
          real qcond = evef * (q1 - qeso[k]) /
                       (R(1.0) + el2orc * qeso[k] / (t1 * t1));
          real dp = -R(100.0) * (p_cup[k + 1] - p_cup[k]);
          if (rn > R(0.0) && qcond < R(0.0)) {
            qevap = -qcond * (R(1.0) - exp(R(-0.32) * sqrt(dtime * rn)));
            qevap = min(qevap, rn * R(1000.0) * C::g / dp);
            delq2 = delqev + R(0.001) * qevap * dp / C::g;
          }
          if (rn > R(0.0) && qcond < R(0.0) && delq2 > rntot) {
            qevap = R(1000.0) * C::g * (rntot - delqev) / dp;
            flg = false;
          }
          if (rn > R(0.0) && qevap > R(0.0)) {
            outq[k] = outq[k] + qevap / dtime;
            outt[k] = outt[k] - elocp * qevap / dtime;
            rn = max(R(0.0), rn - R(0.001) * qevap * dp / C::g);
            pre = pre - qevap * dp / C::g / dtime;
            pre = max(pre, R(0.0));
            delqev = delqev + R(0.001) * dp * qevap / C::g;
          }
        }
      }
    }
  }

  if (ierr == 0) {
    real dts = R(0.0), fpi = R(0.0);
    for (int k = 1; k <= ktop; ++k) {
      real dp = (po_cup[k] - po_cup[k + 1]) * R(100.0);
      dts = dts - (outu[k] * us[k] + outv[k] * vs[k]) * dp / C::g;
      fpi = fpi + sqrt(outu[k] * outu[k] + outv[k] * outv[k]) * dp;
    }
    if (fpi > R(0.0)) {
      for (int k = 1; k <= ktop; ++k) {
        real fp = sqrt(outu[k] * outu[k] + outv[k] * outv[k]) / fpi;
        outt[k] = outt[k] + fp * dts * C::g / C::cp;
      }
    }
  }

  return GflDeepOut<real>{ierr, kbcon, ktop, k22, jmin, frh_out};
}

#undef R

}  // namespace cppwrf::phys::cu_gfl
