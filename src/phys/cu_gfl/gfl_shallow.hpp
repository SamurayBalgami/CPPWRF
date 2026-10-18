#pragma once
// Grell-Freitas-Li shallow convection: WRF 4.8.0 phys/GFL/module_cu_gfl_sh.F90
// cu_gfl_sh_run for one column (ichoice = 3 from module_cu_gfl_wrf.F90,
// make_calc_for_xk = .true.).
//
// In WRF 4.8.0 GFL, get_inversion_layers initializes sec_deriv(:) = 0. at entry
// and resets sec_deriv(:) = 1.e9 inside the column loop, with ilev = max(4, kstart+1),
// so cu_gfl_sh_run is 100% column-local with no row-carried state.

#include <limits>

#include "phys/cu_gfl/gfl_common.hpp"

namespace cppwrf::phys::cu_gfl {

#define R(x) static_cast<real>(x)

template <typename real>
struct GflShallowIn {
  int nz;
  const real *us, *vs, *zo, *t, *q, *tn, *qo, *po, *dhdt, *rho;
  real z1, psur, hfx, qfx, xland, dx, tcrit, dtime;
  int kpbl, ichoice;
};

struct GflShallowOut {
  int ierr, kbcon, ktop, k22;
};

template <typename real>
CPPWRF_HOST_DEVICE inline GflShallowOut cu_gfl_sh_run(
    const GflShallowIn<real>& in, real& frhs, real* zuo, real& xmb_out,
    real* outt, real* outq, real* outqc, real* outu, real* outv, real* cnvwt,
    real& pre, real* cupclw) {
  using C = GflConstants<real>;
  const int kte = in.nz, ktf = in.nz;
  const real* us = in.us;
  const real* vs = in.vs;
  const real* t = in.t;
  const real* q = in.q;
  const real* tn = in.tn;
  const real* qo = in.qo;
  const real* po = in.po;
  const real* dhdt = in.dhdt;
  const real* rho = in.rho;
  const real z1 = in.z1, psur = in.psur, dtime = in.dtime, tcrit = in.tcrit;
  (void)tcrit;
  const int kpbl = in.kpbl;
  const int ichoice = in.ichoice;

  int ierr = 0, kbcon = 0, ktop = 0, k22 = 0, kbmax = 1, kstabi = 0, ktopx = 0;
  int pmin_lev = 0;

  real zo[KM] = {}, z[KM] = {}, xz[KM] = {}, zu[KM] = {}, xzu[KM] = {};
  real qes[KM] = {}, he[KM] = {}, hes[KM] = {};
  real qeso[KM] = {}, heo[KM] = {}, heso[KM] = {};
  real qes_cup[KM] = {}, q_cup[KM] = {}, he_cup[KM] = {}, hes_cup[KM] = {};
  real z_cup[KM] = {}, p_cup[KM] = {}, gamma_cup[KM] = {}, t_cup[KM] = {};
  real qeso_cup[KM] = {}, qo_cup[KM] = {}, heo_cup[KM] = {}, heso_cup[KM] = {};
  real zo_cup[KM] = {}, po_cup[KM] = {}, gammao_cup[KM] = {}, tn_cup[KM] = {};
  real u_cup[KM] = {}, v_cup[KM] = {}, uc[KM] = {}, vc[KM] = {};
  real hc[KM] = {}, dby[KM] = {}, hco[KM] = {}, dbyo[KM] = {}, dbyt[KM] = {};
  real qco[KM] = {}, qrco[KM] = {}, pwo[KM] = {}, c1d[KM] = {};
  real cd[KM] = {}, entr_rate_2d[KM] = {};
  real up_massentr[KM] = {}, up_massdetr[KM] = {};
  real up_massentro[KM] = {}, up_massdetro[KM] = {};
  real up_massentru[KM] = {}, up_massdetru[KM] = {};
  real dellah[KM] = {}, dellaq[KM] = {}, dellat[KM] = {}, dellaqc[KM] = {};
  real dellu[KM] = {}, dellv[KM] = {};
  real xhe[KM] = {}, xq[KM] = {}, xt[KM] = {}, xqes[KM] = {}, xhes[KM] = {};
  real xqes_cup[KM] = {}, xq_cup[KM] = {}, xhe_cup[KM] = {}, xhes_cup[KM] = {};
  real xz_cup[KM] = {}, xt_cup[KM] = {};
  real xhc[KM] = {}, xdby[KM] = {};

  for (int k = 1; k <= kte; ++k) zo[k] = in.zo[k];

  int start_level = 0;
  const real flux_tun = C::fluxtune;
  real lambau = R(2.0);
  int xland1 = static_cast<int>(in.xland + R(0.001));
  if (in.xland > R(1.5) || in.xland < R(0.5)) xland1 = 0;

  pre = R(0.0);
  xmb_out = R(0.0);
  const real cap_max_increment = R(25.0);
  const real entr_rate = R(1.e-3);
  const real radius = R(0.2) / entr_rate;
  frhs = min(R(1.0), R(3.14) * radius * radius / in.dx / in.dx);

  for (int k = 1; k <= ktf; ++k) {
    up_massentro[k] = R(0.0);
    up_massdetro[k] = R(0.0);
    up_massentru[k] = R(0.0);
    up_massdetru[k] = R(0.0);
    z[k] = zo[k];
    xz[k] = zo[k];
    qrco[k] = R(0.0);
    pwo[k] = R(0.0);
    cd[k] = R(0.75) * entr_rate;
    dellaqc[k] = R(0.0);
    cupclw[k] = R(0.0);
    c1d[k] = R(0.0);
  }

  const real cap_maxs = R(175.0);
  real aa0 = R(0.0), aa1 = R(0.0);
  real cap_max = cap_maxs;
  real ztexec = R(0.0), zqexec = R(0.0), zws = R(0.0);
  {
    const real buo_flux =
        (in.hfx / C::cp + R(0.608) * t[1] * in.qfx / C::xlv) / rho[1];
    const real pgeoh = zo[2] * C::g;
    (void)pgeoh;
    zws = max(R(0.0), flux_tun * R(0.41) * buo_flux * zo[2] * C::g / t[1]);
    if (zws > std::numeric_limits<real>::min()) {
      zws = R(1.2) * pow(zws, R(0.3333));
      ztexec = max(flux_tun * in.hfx / (rho[1] * zws * C::cp), R(0.0));
      zqexec = max(flux_tun * in.qfx / C::xlv / (rho[1] * zws), R(0.0));
    }
    zws =
        max(R(0.0), flux_tun * R(0.41) * buo_flux * zo[kpbl] * C::g / t[kpbl]);
    zws = R(1.2) * pow(zws, R(0.3333));
    zws = zws * rho[kpbl];
  }

  const real zkbmax = R(3000.0);
  cup_env(ktf, kte, ierr, -1, psur, z1, po, t, z, he, hes, qes, q);
  cup_env(ktf, kte, ierr, -1, psur, z1, po, tn, zo, heo, heso, qeso, qo);
  cup_env_clev(ktf, kte, ierr, psur, z1, t, qes, q, he, hes, z, po, qes_cup,
               q_cup, he_cup, hes_cup, z_cup, p_cup, gamma_cup, t_cup);
  cup_env_clev(ktf, kte, ierr, psur, z1, tn, qeso, qo, heo, heso, zo, po,
               qeso_cup, qo_cup, heo_cup, heso_cup, zo_cup, po_cup, gammao_cup,
               tn_cup);

  if (ierr == 0) {
    u_cup[1] = us[1];
    v_cup[1] = vs[1];
    for (int k = 2; k <= ktf; ++k) {
      u_cup[k] = R(0.5) * (us[k - 1] + us[k]);
      v_cup[k] = R(0.5) * (vs[k - 1] + vs[k]);
    }
  }
  if (ierr == 0) {
    for (int k = 1; k <= ktf; ++k) {
      if (zo_cup[k] > zkbmax + z1) {
        kbmax = k;
        break;
      }
    }
    kbmax = min(kbmax, ktf / 2);
  }

  if (kpbl > 3) cap_max = po_cup[kpbl];
  if (ierr == 0) {
    k22 = maxloc_abs(heo_cup, 2, kbmax) - 1;
    k22 = max(2, k22);
    if (k22 > kbmax) {
      ierr = 2;
      ktop = 0;
      k22 = 0;
      kbcon = 0;
    }
  }

  real hkb = R(0.0), hkbo = R(0.0);
  if (ierr == 0) {
    const real x_add = C::xlv * zqexec + C::cp * ztexec;
    get_cloud_bc(kte, he_cup, hkb, k22, x_add);
    get_cloud_bc(kte, heo_cup, hkbo, k22, x_add);
  }

  for (int k = 1; k <= ktf; ++k) dbyo[k] = R(0.0);

  cup_kbcon(cap_max_increment, 5, k22, kbcon, heo_cup, heso_cup, hkbo, ierr,
            kbmax, po_cup, cap_max, ztexec, zqexec, kte, z_cup, entr_rate, heo,
            0);
  cup_minimi(heso_cup, kbcon, kbmax, ierr, kstabi);

  int kinv1 = 0, kinv2 = 0;
  get_inversion_layers(ierr, p_cup, t_cup, z_cup, kbcon, kstabi, kte, kinv1,
                       kinv2);

  for (int k = 1; k <= kte; ++k) entr_rate_2d[k] = entr_rate;
  if (ierr == 0) {
    start_level = k22;
    const real x_add = C::xlv * zqexec + C::cp * ztexec;
    get_cloud_bc(kte, he_cup, hkb, k22, x_add);
    if (kbcon > ktf - 4) ierr = 231;
    for (int k = 1; k <= ktf; ++k) {
      entr_rate_2d[k] = entr_rate;
      cd[k] = R(0.75) * entr_rate_2d[k];
    }
    ktop = 1;
    int kstart = kpbl;
    if (kpbl < 5) kstart = kbcon;
    if (kinv1 > 0 && (po_cup[kstart] - po_cup[kinv1]) < R(200.0)) {
      ktop = kinv1;
    } else {
      for (int k = kbcon + 1; k <= ktf; ++k) {
        if ((po_cup[kstart] - po_cup[k]) > R(200.0)) {
          ktop = k;
          break;
        }
      }
    }
  }

  rates_up_pdf(2, ktf, kte, ierr, po_cup, ktop, k22, zo_cup, entr_rate_2d, heo,
               heso_cup, hkbo, xland1, kbcon, zuo, ktopx, kbcon, pmin_lev);

  if (ierr == 0) {
    if (k22 > 1) {
      for (int k = 1; k <= k22 - 1; ++k) {
        zuo[k] = R(0.0);
        zu[k] = R(0.0);
        xzu[k] = R(0.0);
      }
    }
    for (int k = maxloc_abs(zuo, 1, kte); k <= ktop; ++k) {
      if (zuo[k] < R(1.e-6)) {
        ktop = k - 1;
        break;
      }
    }
    for (int k = k22; k <= ktop; ++k) {
      xzu[k] = zuo[k];
      zu[k] = zuo[k];
    }
    for (int k = ktop + 1; k <= ktf; ++k) {
      zuo[k] = R(0.0);
      zu[k] = R(0.0);
      xzu[k] = R(0.0);
    }
    k22 = max(2, k22);
  }

  get_lateral_massflux(ktf, kte, ierr, ktop, zo_cup, zuo, cd, entr_rate_2d,
                       up_massentro, up_massdetro, up_massentr, up_massdetr, 2,
                       kbcon, k22, up_massentru, up_massdetru, lambau);

  for (int k = 1; k <= ktf; ++k) {
    hc[k] = R(0.0);
    qco[k] = R(0.0);
    qrco[k] = R(0.0);
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

  for (int k = 1; k <= kte; ++k) dbyt[k] = R(0.0);
  if (ierr == 0) {
    for (int k = start_level + 1; k <= ktop; ++k) {
      hc[k] = (hc[k - 1] * zu[k - 1] - R(0.5) * up_massdetr[k - 1] * hc[k - 1] +
               up_massentr[k - 1] * he[k - 1]) /
              (zu[k - 1] - R(0.5) * up_massdetr[k - 1] + up_massentr[k - 1]);
      uc[k] = (uc[k - 1] * zu[k - 1] - R(0.5) * up_massdetr[k - 1] * uc[k - 1] +
               up_massentr[k - 1] * us[k - 1]) /
              (zu[k - 1] - R(0.5) * up_massdetr[k - 1] + up_massentr[k - 1]);
      vc[k] = (vc[k - 1] * zu[k - 1] - R(0.5) * up_massdetr[k - 1] * vc[k - 1] +
               up_massentr[k - 1] * vs[k - 1]) /
              (zu[k - 1] - R(0.5) * up_massdetr[k - 1] + up_massentr[k - 1]);
      dby[k] = max(R(0.0), hc[k] - hes_cup[k]);
      hco[k] =
          (hco[k - 1] * zuo[k - 1] - R(0.5) * up_massdetro[k - 1] * hco[k - 1] +
           up_massentro[k - 1] * heo[k - 1]) /
          (zuo[k - 1] - R(0.5) * up_massdetro[k - 1] + up_massentro[k - 1]);
      dbyo[k] = hco[k] - heso_cup[k];
      const real dz = zo_cup[k + 1] - zo_cup[k];
      if (k >= kbcon) dbyt[k] = dbyt[k - 1] + dbyo[k] * dz;
    }
    const int ki = maxloc_abs(dbyt, 1, kte);
    if (ktop > ki + 1) {
      ktop = ki + 1;
      for (int k = ktop + 1; k <= ktf; ++k) {
        zuo[k] = R(0.0);
        zu[k] = R(0.0);
        cd[k] = R(0.0);
      }
      up_massdetro[ktop] = zuo[ktop];
      for (int k = ktop; k <= ktf; ++k) up_massentro[k] = R(0.0);
      for (int k = ktop + 1; k <= ktf; ++k) {
        up_massdetro[k] = R(0.0);
        entr_rate_2d[k] = R(0.0);
      }
    }
    bool go42 = false;
    if (ktop < kbcon + 1) {
      ierr = 5;
      go42 = true;
    }
    if (!go42 && ktop > ktf - 2) {
      ierr = 5;
      go42 = true;
    }
    if (!go42) {
      real qaver = R(0.0);
      get_cloud_bc(kte, qo_cup, qaver, k22, R(0.0));
      qaver = qaver + zqexec;
      for (int k = 1; k <= start_level - 1; ++k) qco[k] = qo_cup[k];
      qco[start_level] = qaver;
      for (int k = start_level + 1; k <= ktop; ++k) {
        const real trash = qeso_cup[k] + (R(1.0) / C::xlv) *
                                             (gammao_cup[k] /
                                              (R(1.0) + gammao_cup[k])) *
                                             dbyo[k];
        const real trash2 = qco[k - 1];
        qco[k] = (trash2 * (zuo[k - 1] - R(0.5) * up_massdetr[k - 1]) +
                  up_massentr[k - 1] * qo[k - 1]) /
                 (zuo[k - 1] - R(0.5) * up_massdetr[k - 1] +
                  up_massentr[k - 1]);
        if (qco[k] >= trash) {
          const real dz = z_cup[k] - z_cup[k - 1];
          c1d[k] = R(0.02) * up_massdetr[k - 1];
          qrco[k] = (qco[k] - trash) / (R(1.0) + (C::c0_shal + c1d[k]) * dz);
          if (qrco[k] < R(0.0)) {
            qrco[k] = R(0.0);
            c1d[k] = R(0.0);
          }
          pwo[k] = C::c0_shal * dz * qrco[k] * zuo[k];
          qco[k] = trash + qrco[k];
        } else {
          qrco[k] = R(0.0);
        }
        cupclw[k] = qrco[k];
      }
      for (int k = k22 + 1; k <= ktop; ++k) {
        const real dp = R(100.0) * (po_cup[k] - po_cup[k + 1]);
        cnvwt[k] = zuo[k] * cupclw[k] * C::g / dp;
        qco[k] = qco[k] - qrco[k];
      }
      for (int k = ktop + 1; k <= ktf - 1; ++k) {
        hc[k] = hes_cup[k];
        hco[k] = heso_cup[k];
        qco[k] = qeso_cup[k];
        uc[k] = u_cup[k];
        vc[k] = v_cup[k];
        qrco[k] = R(0.0);
        dby[k] = R(0.0);
        dbyo[k] = R(0.0);
        zu[k] = R(0.0);
        xzu[k] = R(0.0);
        zuo[k] = R(0.0);
      }
    }
  }

  cup_up_aa0(ktf, kte, aa0, z, zu, dby, gamma_cup, t_cup, kbcon, ktop, ierr);
  cup_up_aa0(ktf, kte, aa1, zo, zuo, dbyo, gammao_cup, tn_cup, kbcon, ktop,
             ierr);
  if (ierr == 0 && aa1 <= R(0.0)) ierr = 17;

  for (int k = 1; k <= kte; ++k) {
    dellah[k] = R(0.0);
    dellaq[k] = R(0.0);
    dellaqc[k] = R(0.0);
    dellu[k] = R(0.0);
    dellv[k] = R(0.0);
  }
  if (ierr == 0) {
    const real dp1 = R(100.0) * (po_cup[1] - po_cup[2]);
    dellu[1] = -zuo[2] * (uc[2] - u_cup[2]) * C::g / dp1;
    dellv[1] = -zuo[2] * (vc[2] - v_cup[2]) * C::g / dp1;
    dellah[1] = -zuo[2] * (hco[2] - heo_cup[2]) * C::g / dp1;
    dellaq[1] = -zuo[2] * (qco[2] - qo_cup[2]) * C::g / dp1;

    for (int k = k22; k <= ktop; ++k) {
      const real detup = up_massdetro[k];
      const real dp = R(100.0) * (po_cup[k] - po_cup[k + 1]);
      dellah[k] = -((zuo[k + 1] * (hco[k + 1] - heo_cup[k + 1]) -
                     zuo[k] * (hco[k] - heo_cup[k])) *
                    C::g / dp);
      const real dz = zo_cup[k + 1] - zo_cup[k];
      if (k < ktop && c1d[k] > R(0.0)) {
        dellaqc[k] = zuo[k] * c1d[k] * qrco[k] * dz / dp * C::g;
      } else {
        dellaqc[k] = detup * R(0.5) * (qrco[k + 1] + qrco[k]) * C::g / dp;
      }
      const real c_up =
          dellaqc[k] + (zuo[k + 1] * qrco[k + 1] - zuo[k] * qrco[k]) * C::g / dp;
      dellaq[k] = -((zuo[k + 1] * (qco[k + 1] - qo_cup[k + 1]) -
                     zuo[k] * (qco[k] - qo_cup[k])) *
                    C::g / dp) -
                  c_up - R(0.5) * (pwo[k] + pwo[k + 1]) * C::g / dp;
      dellu[k] = -((zuo[k + 1] * (uc[k + 1] - u_cup[k + 1]) -
                    zuo[k] * (uc[k] - u_cup[k])) *
                   C::g / dp);
      dellv[k] = -((zuo[k + 1] * (vc[k + 1] - v_cup[k + 1]) -
                    zuo[k] * (vc[k] - v_cup[k])) *
                   C::g / dp);
    }
  }

  const real mbdt = R(0.5);
  for (int k = 1; k <= ktf; ++k) {
    dellat[k] = R(0.0);
    if (ierr != 0) continue;
    xhe[k] = dellah[k] * mbdt + heo[k];
    xq[k] = max(R(1.e-16), (dellaq[k] + dellaqc[k]) * mbdt + qo[k]);
    dellat[k] = (R(1.0) / C::cp) * (dellah[k] - C::xlv * dellaq[k]);
    xt[k] = (-dellaqc[k] * C::xlv / C::cp + dellat[k]) * mbdt + tn[k];
    xt[k] = max(R(190.0), xt[k]);
  }
  if (ierr == 0) {
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
    const real x_add = C::xlv * zqexec + C::cp * ztexec;
    real xhkb = R(0.0);
    get_cloud_bc(kte, xhe_cup, xhkb, k22, x_add);
    for (int k = 1; k <= start_level - 1; ++k) xhc[k] = xhe_cup[k];
    xhc[start_level] = xhkb;
  }
  if (ierr == 0) {
    for (int k = 1; k <= ktf; ++k) xzu[k] = zuo[k];
    for (int k = start_level + 1; k <= ktop; ++k) {
      xhc[k] =
          (xhc[k - 1] * xzu[k - 1] - R(0.5) * up_massdetro[k - 1] * xhc[k - 1] +
           up_massentro[k - 1] * xhe[k - 1]) /
          (xzu[k - 1] - R(0.5) * up_massdetro[k - 1] + up_massentro[k - 1]);
      xdby[k] = xhc[k] - xhes_cup[k];
    }
    for (int k = ktop + 1; k <= ktf; ++k) {
      xhc[k] = xhes_cup[k];
      xdby[k] = R(0.0);
      xzu[k] = R(0.0);
    }
  }

  real xaa0 = R(0.0);
  cup_up_aa0(ktf, kte, xaa0, xz, xzu, xdby, gamma_cup, xt_cup, kbcon, ktop,
             ierr);

  real xmb = R(0.0);
  real xff_shal[4] = {};
  if (ierr == 0) {
    const real xmbmax = R(1.0);
    real xkshal = (xaa0 - aa1) / mbdt;
    if (xkshal <= R(0.0) && xkshal > R(-0.01) * mbdt) xkshal = R(-0.01) * mbdt;
    if (xkshal > R(0.0) && xkshal < R(1.e-2)) xkshal = R(1.e-2);
    xff_shal[1] = max(R(0.0), -(aa1 - aa0) / (xkshal * dtime));
    xff_shal[2] = R(0.03) * zws;
    real blqe = R(0.0);
    for (int k = 1; k <= kbcon; ++k) {
      blqe = blqe + R(100.0) * dhdt[k] * (po_cup[k] - po_cup[k + 1]) / C::g;
    }
    const real trash = max(hc[kbcon] - he_cup[kbcon], R(1.e1));
    xff_shal[3] = max(R(0.0), blqe / trash);
    xff_shal[3] = min(xmbmax, xff_shal[3]);
    xmb = (xff_shal[1] + xff_shal[2] + xff_shal[3]) / R(3.0);
    xmb = min(xmbmax, xmb);
    if (ichoice > 0) xmb = min(xmbmax, xff_shal[ichoice]);
    if (xmb <= R(0.0)) ierr = 21;
  }

  if (ierr != 0) {
    k22 = 0;
    kbcon = 0;
    ktop = 0;
    xmb = R(0.0);
    for (int k = 1; k <= kte; ++k) {
      outt[k] = R(0.0);
      outu[k] = R(0.0);
      outv[k] = R(0.0);
      outq[k] = R(0.0);
      outqc[k] = R(0.0);
    }
  } else {
    xmb_out = xmb;
    pre = R(0.0);
    for (int k = 2; k <= ktop; ++k) {
      outt[k] = dellat[k] * xmb;
      outq[k] = dellaq[k] * xmb;
      outqc[k] = dellaqc[k] * xmb;
      pre = pre + pwo[k] * xmb;
    }
    outt[1] = dellat[1] * xmb;
    outq[1] = dellaq[1] * xmb;
    outu[1] = dellu[1] * xmb;
    outv[1] = dellv[1] * xmb;
    for (int k = 2; k <= ktop; ++k) {
      outu[k] = R(0.25) * (dellu[k - 1] + R(2.0) * dellu[k] + dellu[k + 1]) * xmb;
      outv[k] = R(0.25) * (dellv[k - 1] + R(2.0) * dellv[k] + dellv[k + 1]) * xmb;
    }
  }

  if (ierr == 0) {
    real dts = R(0.0);
    real fpi = R(0.0);
    for (int k = 1; k <= ktop; ++k) {
      const real dp = (po_cup[k] - po_cup[k + 1]) * R(100.0);
      dts = dts - (outu[k] * us[k] + outv[k] * vs[k]) * dp / C::g;
      fpi = fpi + sqrt(outu[k] * outu[k] + outv[k] * outv[k]) * dp;
    }
    if (fpi > R(0.0)) {
      for (int k = 1; k <= ktop; ++k) {
        const real fp = sqrt(outu[k] * outu[k] + outv[k] * outv[k]) / fpi;
        outt[k] = outt[k] + fp * dts * C::g / C::cp;
      }
    }
  }

  return GflShallowOut{ierr, kbcon, ktop, k22};
}

#undef R

}  // namespace cppwrf::phys::cu_gfl
