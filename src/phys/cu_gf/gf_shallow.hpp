#pragma once
// Grell-Freitas shallow convection: WRF 4.6.1 module_cu_gf_sh.F CUP_gf_sh for one column
// (ichoice_s = 0, MAKE_CALC_FOR_XK = .true.).
//
// get_inversion_layers is not column-local in WRF: its scratch sec_deriv(kts:kte) is reset to 1e9 once
// per call (= per grid row j) and then reused column after column, so the shallow minimum search of
// column i also sees the "mid" values left by the last column processed before it in the same row.
// The C++ therefore runs the shallow scheme in two passes: stage 0 (per column) up to the inversion
// layers (GfInvInfo), the row scan gf_inversion_row (per row, in i order), stage 1 (per column) the rest.
// The first loop of get_inversion_layers reads sec_deriv(2), which WRF never sets (kbcon = 2 only);
// it is taken as `sec2_garbage`.
#include "phys/cu_gf/gf_common.hpp"

namespace cppwrf::phys::cu_gf {

inline constexpr int kInvMax = 12;

template <typename T>
struct GfInvInfo {
  int valid;         // ierr == 0 when get_inversion_layers ran
  int ken;           // maxloc(k_inv_layers) after the merge of the first loop
  int kinv[kInvMax];  // k_inv_layers(1..kInvMax-1) after the first loop (1-based, [0] unused)
  T dp[kInvMax];     // p_cup(k_inv_layers(k)) - p_cup(kstart), k = 1..ken
  int kinv1;         // result of the row scan: k_inv_layers(i,1)
};

// get_inversion_layers, first DO i loop, for one column (ierr == 0)
template <typename T>
CPPWRF_HD inline void gf_inversion_column(const T* p_cup, const T* t_cup, const T* z_cup, int kstart, int kend,
                                          int kte, T sec2_garbage, GfInvInfo<T>& inv) {
  int kinv[KM];
  for (int k = 0; k < KM; ++k) kinv[k] = 1;
  T first[KM] = {}, sec[KM] = {}, dtempdz[KM] = {};
  const int kend_p3 = kend + 3;
  for (int k = 2; k <= kend_p3 + 4; ++k) {
    first[k] = (t_cup[k + 1] - t_cup[k - 1]) / (z_cup[k + 1] - z_cup[k - 1]);
    dtempdz[k] = first[k];
  }
  sec[2] = sec2_garbage;
  for (int k = 3; k <= kend_p3 + 3; ++k) {
    sec[k] = (first[k + 1] - first[k - 1]) / (z_cup[k + 1] - z_cup[k - 1]);
    sec[k] = std::fabs(sec[k]);
  }
  int ilev = 3 > kstart + 1 ? 3 : kstart + 1;
  int ix = 1;
  int k = ilev;
  while (ilev < kend_p3) {
    for (int kk = k; kk <= kend_p3 + 2; ++kk) {
      if (sec[kk] < sec[kk + 1] && sec[kk] < sec[kk - 1]) {
        kinv[ix] = kk;
        ix = ix + 1 < 5 ? ix + 1 : 5;
        ilev = kk + 1;
        break;
      }
      ilev = kk + 1;
    }
    k = ilev;
  }
  int kadd = 0;
  int ken = 1;
  for (int kk = 2; kk <= kte; ++kk)
    if (kinv[kk] > kinv[ken]) ken = kk;
  for (int kq = 1; kq <= ken; ++kq) {
    const int kk = kinv[kq + kadd];
    if (kk == 1) break;
    if (dtempdz[kk] < dtempdz[kk - 1] && dtempdz[kk] < dtempdz[kk + 1]) {
      kadd = kadd + 1;
      for (int kj = kq; kj <= ken; ++kj) {
        if (kinv[kj + kadd] > 1) kinv[kj] = kinv[kj + kadd];
        if (kinv[kj + kadd] == 1) kinv[kj] = 1;
      }
    }
  }
  int ken2 = 1;
  for (int kk = 2; kk <= kte; ++kk)
    if (kinv[kk] > kinv[ken2]) ken2 = kk;
  inv.valid = 1;
  inv.ken = ken2;
  for (int kk = 0; kk < kInvMax; ++kk) inv.kinv[kk] = kinv[kk];
  for (int kk = 0; kk < kInvMax; ++kk) inv.dp[kk] = T(0);
  for (int kk = 1; kk <= ken2 && kk < kInvMax; ++kk) inv.dp[kk] = p_cup[kinv[kk]] - p_cup[kstart];
}

// get_inversion_layers, second DO i loop, over one row (columns in i order); sets inv[i].kinv1
template <typename T>
CPPWRF_HD inline void gf_inversion_row(GfInvInfo<T>* inv, int n, int kte) {
  T sec[KM];
  for (int k = 1; k <= kte; ++k) sec[k] = T(1.e9);
  for (int i = 0; i < n; ++i) {
    GfInvInfo<T>& c = inv[i];
    if (!c.valid) continue;
    for (int k = 1; k <= c.ken; ++k) sec[k] = std::fabs(c.dp[k]) - T(100.);
    int k800 = 1;
    for (int k = 2; k <= kte; ++k)
      if (std::fabs(sec[k]) < std::fabs(sec[k800])) k800 = k;
    for (int k = 1; k <= kte; ++k) sec[k] = T(1.e9);
    for (int k = 1; k <= c.ken; ++k) sec[k] = std::fabs(c.dp[k]) - T(300.);
    c.kinv1 = k800 < kInvMax ? c.kinv[k800] : 1;
  }
}

template <typename T>
struct GfShallowIn {
  int nz;
  const T *t, *q, *tn, *qo, *po, *rho, *zo, *dhdt;
  T z1, psur, hfx, qfx, xland, dtime;
  int kpbl;
};

struct GfShallowOut {
  int ierr, kbcon, ktop, k22;
};

// stage 0: up to get_inversion_layers (fills inv); stage 1: the whole scheme with inv.kinv1
template <typename T>
CPPWRF_HD inline GfShallowOut cup_gf_sh(const GfShallowIn<T>& in, int stage, GfInvInfo<T>& inv, T sec2_garbage,
                                        T* zuo, T& xmb_out, T* outt, T* outq, T* outqc, T* cnvwt, T& pre,
                                        T* cupclw) {
  using K = GfK<T>;
  const int kte = in.nz, ktf = in.nz;
  const T* t = in.t; const T* q = in.q; const T* tn = in.tn; const T* qo = in.qo; const T* po = in.po;
  const T* rho = in.rho; const T* dhdt = in.dhdt;
  const T z1 = in.z1, psur = in.psur, dtime = in.dtime;
  const int kpbl = in.kpbl;
  int ierr = 0, kbcon = 0, ktop = 0, k22 = 0, kbmax = 1, kstabi = 0, ktopx = 0;

  T zo[KM] = {}, z[KM] = {}, xz[KM] = {}, zu[KM] = {}, xzu[KM] = {};
  T qes[KM] = {}, he[KM] = {}, hes[KM] = {}, qeso[KM] = {}, heo[KM] = {}, heso[KM] = {};
  T qes_cup[KM] = {}, q_cup[KM] = {}, he_cup[KM] = {}, hes_cup[KM] = {}, z_cup[KM] = {}, p_cup[KM] = {};
  T gamma_cup[KM] = {}, t_cup[KM] = {};
  T qeso_cup[KM] = {}, qo_cup[KM] = {}, heo_cup[KM] = {}, heso_cup[KM] = {}, zo_cup[KM] = {}, po_cup[KM] = {};
  T gammao_cup[KM] = {}, tn_cup[KM] = {};
  T hc[KM] = {}, dby[KM] = {}, hco[KM] = {}, dbyo[KM] = {}, dbyt[KM] = {}, qco[KM] = {}, qrco[KM] = {}, pwo[KM] = {};
  T cd[KM] = {}, entr_rate_2d[KM] = {};
  T up_massentr[KM] = {}, up_massdetr[KM] = {}, up_massentro[KM] = {}, up_massdetro[KM] = {};
  T dellah[KM] = {}, dellaq[KM] = {}, dellat[KM] = {}, dellaqc[KM] = {};
  T xhe[KM] = {}, xq[KM] = {}, xt[KM] = {}, xqes[KM] = {}, xhes[KM] = {};
  T xqes_cup[KM] = {}, xq_cup[KM] = {}, xhe_cup[KM] = {}, xhes_cup[KM] = {}, xz_cup[KM] = {}, xt_cup[KM] = {};
  T xhc[KM] = {}, xdby[KM] = {};

  for (int k = 1; k <= kte; ++k) zo[k] = in.zo[k];
  int start_level = 0;
  const T flux_tun = K::fluxtune;
  int xland1 = int(in.xland + T(.001));
  if (in.xland > T(1.5) || in.xland < T(0.5)) xland1 = 0;
  (void)xland1;
  pre = T(0.);
  xmb_out = T(0.);
  const T cap_max_increment = T(25.);
  const T entr_rate = T(9.e-5);
  for (int k = 1; k <= ktf; ++k) {
    up_massentro[k] = 0; up_massdetro[k] = 0; z[k] = zo[k]; xz[k] = zo[k]; qrco[k] = 0; pwo[k] = 0;
    cd[k] = T(1.) * entr_rate; dellaqc[k] = 0; cupclw[k] = 0;
  }
  const T cap_maxs = T(125.);
  T aa0 = T(0.), aa1 = T(0.);
  T cap_max = cap_maxs;
  T ztexec = T(0.), zqexec = T(0.), zws = T(0.);
  {
    const T buo_flux = (in.hfx / K::cp + T(0.608) * t[1] * in.qfx / K::xlv) / rho[1];
    zws = fmx(T(0.), flux_tun * T(0.41) * buo_flux * zo[2] * K::g / t[1]);
    if (zws > std::numeric_limits<T>::min()) {
      zws = T(1.2) * std::pow(zws, T(.3333));
      ztexec = fmx(flux_tun * in.hfx / (rho[1] * zws * K::cp), T(0.0));
      zqexec = fmx(flux_tun * in.qfx / K::xlv / (rho[1] * zws), T(0.));
    }
    zws = fmx(T(0.), flux_tun * T(0.41) * buo_flux * zo[kpbl] * K::g / t[kpbl]);
    zws = T(1.2) * std::pow(zws, T(.3333));
    zws = zws * rho[kpbl];
  }
  const T zkbmax = T(3000.);
  cup_env(z, qes, he, hes, t, q, po, ierr, ktf);
  cup_env(zo, qeso, heo, heso, tn, qo, po, ierr, ktf);
  cup_env_clev(t, qes, q, he, hes, z, po, qes_cup, q_cup, he_cup, hes_cup, z_cup, p_cup, gamma_cup, t_cup, psur, ierr,
               z1, ktf);
  cup_env_clev(tn, qeso, qo, heo, heso, zo, po, qeso_cup, qo_cup, heo_cup, heso_cup, zo_cup, po_cup, gammao_cup,
               tn_cup, psur, ierr, z1, ktf);
  if (ierr == 0) {
    for (int k = 1; k <= ktf; ++k)
      if (zo_cup[k] > zkbmax + z1) {
        kbmax = k;
        break;
      }
    kbmax = kbmax < ktf / 2 ? kbmax : ktf / 2;
  }
  if (kpbl > 3) cap_max = po_cup[kpbl];
  if (ierr == 0) {
    k22 = maxloc_abs(heo_cup, 2, kbmax) - 1;  // MAXLOC(HEO_CUP(i,2:kbmax)) is relative to the slice
    k22 = k22 > 2 ? k22 : 2;
    if (k22 > kbmax) {
      ierr = 2;
      ktop = 0;
      k22 = 0;
      kbcon = 0;
    }
  }
  T hkb = T(0), hkbo = T(0);
  if (ierr == 0) {
    const T x_add = K::xlv * zqexec + K::cp * ztexec;
    hkb = get_cloud_bc(he_cup, k22, true, x_add);
    hkbo = get_cloud_bc(heo_cup, k22, true, x_add);
  }
  for (int k = 1; k <= ktf; ++k) dbyo[k] = 0;
  cup_kbcon(cap_max_increment, 5, k22, kbcon, heo_cup, heso_cup, hkbo, ierr, kbmax, po_cup, cap_max, ztexec, zqexec,
            z_cup, entr_rate, heo, 0);
  kstabi = cup_minimi(heso_cup, kbcon, kbmax, ierr);
  if (stage == 0) {
    inv.valid = 0;
    inv.ken = 1;
    inv.kinv1 = 1;
    if (ierr == 0) gf_inversion_column(p_cup, t_cup, z_cup, kbcon, kstabi, kte, sec2_garbage, inv);
    return GfShallowOut{ierr, kbcon, ktop, k22};
  }
  const int kinv1 = ierr == 0 ? inv.kinv1 : 1;
  for (int k = 1; k <= kte; ++k) entr_rate_2d[k] = entr_rate;
  if (ierr == 0) {
    start_level = k22;
    const T x_add = K::xlv * zqexec + K::cp * ztexec;
    hkb = get_cloud_bc(he_cup, k22, true, x_add);
    if (kbcon > ktf - 4) ierr = 231;
    for (int k = 1; k <= ktf; ++k) {
      const T frh = T(2.) * fmn(qo_cup[k] / qeso_cup[k], T(1.));
      entr_rate_2d[k] = entr_rate * (T(2.3) - frh);
      cd[k] = entr_rate_2d[k];
    }
    ktop = 1;
    if (kinv1 > 0 && (po_cup[kbcon] - po_cup[kinv1]) < T(200.)) {
      ktop = kinv1;
    } else {
      for (int k = kbcon + 1; k <= ktf; ++k)
        if ((po_cup[kbcon] - po_cup[k]) > T(200.)) {
          ktop = k;
          break;
        }
    }
  }
  rates_up_pdf(true, ktop, ierr, po_cup, entr_rate_2d, hkbo, heo, heso_cup, zo_cup, kstabi, k22, kbcon, kte, ktf, zuo,
               kpbl, ktopx, ktopx);
  if (ierr == 0) {
    if (k22 > 1)
      for (int k = 1; k <= k22 - 1; ++k) {
        zuo[k] = 0; zu[k] = 0; xzu[k] = 0;
      }
    for (int k = maxloc_abs(zuo, 1, kte); k <= ktop; ++k)
      if (zuo[k] < T(1.e-6)) {
        ktop = k - 1;
        break;
      }
    for (int k = k22; k <= ktop; ++k) {
      xzu[k] = zuo[k];
      zu[k] = zuo[k];
    }
    for (int k = ktop + 1; k <= ktf; ++k) {
      zuo[k] = 0; zu[k] = 0; xzu[k] = 0;
    }
    k22 = k22 > 2 ? k22 : 2;
  }
  get_lateral_massflux(ktf, kte, ierr, ktop, zo_cup, zuo, cd, entr_rate_2d, up_massentro, up_massdetro, up_massentr,
                       up_massdetr, k22, false, (T*)nullptr, (T*)nullptr, T(0));
  for (int k = 1; k <= ktf; ++k) {
    hc[k] = 0; qco[k] = 0; qrco[k] = 0; dby[k] = 0; hco[k] = 0; dbyo[k] = 0;
  }
  if (ierr == 0) {
    for (int k = 1; k <= start_level - 1; ++k) {
      hc[k] = he_cup[k];
      hco[k] = heo_cup[k];
    }
    hc[start_level] = hkb;
    hco[start_level] = hkbo;
  }
  for (int k = 1; k <= kte; ++k) dbyt[k] = 0;
  if (ierr == 0) {
    for (int k = start_level + 1; k <= ktop; ++k) {
      hc[k] = (hc[k - 1] * zu[k - 1] - T(.5) * up_massdetr[k - 1] * hc[k - 1] + up_massentr[k - 1] * he[k - 1]) /
              (zu[k - 1] - T(.5) * up_massdetr[k - 1] + up_massentr[k - 1]);
      dby[k] = fmx(T(0.), hc[k] - hes_cup[k]);
      hco[k] = (hco[k - 1] * zuo[k - 1] - T(.5) * up_massdetro[k - 1] * hco[k - 1] +
                up_massentro[k - 1] * heo[k - 1]) /
               (zuo[k - 1] - T(.5) * up_massdetro[k - 1] + up_massentro[k - 1]);
      dbyo[k] = hco[k] - heso_cup[k];
      const T dz = zo_cup[k + 1] - zo_cup[k];
      dbyt[k] = dbyt[k - 1] + dbyo[k] * dz;
    }
    const int ki = maxloc_abs(dbyt, 1, kte);
    if (ktop > ki + 1) {
      ktop = ki + 1;
      for (int k = ktop + 1; k <= ktf; ++k) {
        zuo[k] = 0; zu[k] = 0; cd[k] = 0;
      }
      up_massdetro[ktop] = zuo[ktop];
      for (int k = ktop; k <= ktf; ++k) up_massentro[k] = 0;
      for (int k = ktop + 1; k <= ktf; ++k) {
        up_massdetro[k] = 0;
        entr_rate_2d[k] = 0;
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
      T qaver = get_cloud_bc(qo_cup, k22, false, T(0));
      qaver = qaver + zqexec;
      for (int k = 1; k <= start_level - 1; ++k) qco[k] = qo_cup[k];
      qco[start_level] = qaver;
      for (int k = start_level + 1; k <= ktop; ++k) {
        const T trash = qeso_cup[k] + (T(1.) / K::xlv) * (gammao_cup[k] / (T(1.) + gammao_cup[k])) * dbyo[k];
        const T trash2 = qco[k - 1];
        qco[k] = (trash2 * (zuo[k - 1] - T(0.5) * up_massdetr[k - 1]) + up_massentr[k - 1] * qo[k - 1]) /
                 (zuo[k - 1] - T(.5) * up_massdetr[k - 1] + up_massentr[k - 1]);
        if (qco[k] >= trash) {
          const T dz = z_cup[k] - z_cup[k - 1];
          qrco[k] = (qco[k] - trash) / (T(1.) + (K::c0_shal + K::c1_shal) * dz);
          pwo[k] = K::c0_shal * dz * qrco[k] * zuo[k];
          qco[k] = trash + qrco[k];
        } else {
          qrco[k] = T(0.0);
        }
        cupclw[k] = qrco[k];
      }
      for (int k = k22 + 1; k <= ktop; ++k) {
        const T dp = T(100.) * (po_cup[k] - po_cup[k + 1]);
        cnvwt[k] = zuo[k] * cupclw[k] * K::g / dp;
        qco[k] = qco[k] - qrco[k];
      }
      for (int k = ktop + 1; k <= ktf - 1; ++k) {
        hc[k] = hes_cup[k]; hco[k] = heso_cup[k]; qco[k] = qeso_cup[k]; qrco[k] = 0;
        dby[k] = 0; dbyo[k] = 0; zu[k] = 0; xzu[k] = 0; zuo[k] = 0;
      }
    }
  }
  aa0 = cup_up_aa0(z, zu, dby, gamma_cup, t_cup, kbcon, ktop, ierr, ktf);
  aa1 = cup_up_aa0(zo, zuo, dbyo, gammao_cup, tn_cup, kbcon, ktop, ierr, ktf);
  if (ierr == 0 && aa1 <= T(0.)) ierr = 17;
  for (int k = 1; k <= kte; ++k) {
    dellah[k] = 0;
    dellaq[k] = 0;
  }
  if (ierr == 0)
    for (int k = k22; k <= ktop; ++k) {
      const T detup = up_massdetro[k];
      const T dp = T(100.) * (po_cup[k] - po_cup[k + 1]);
      dellah[k] = -((zuo[k + 1] * (hco[k + 1] - heo_cup[k + 1]) - zuo[k] * (hco[k] - heo_cup[k])) * K::g / dp);
      const T dz = zo_cup[k + 1] - zo_cup[k];
      if (k < ktop) dellaqc[k] = zuo[k] * K::c1_shal * qrco[k] * dz / dp * K::g;
      else dellaqc[k] = detup * qrco[k] * K::g / dp;
      const T c_up = dellaqc[k] + (zuo[k + 1] * qrco[k + 1] - zuo[k] * qrco[k]) * K::g / dp;
      dellaq[k] = -((zuo[k + 1] * (qco[k + 1] - qo_cup[k + 1]) - zuo[k] * (qco[k] - qo_cup[k])) * K::g / dp) - c_up -
                  T(0.5) * (pwo[k] + pwo[k + 1]) * K::g / dp;
    }
  const T mbdt = T(.5);
  for (int k = 1; k <= ktf; ++k) {
    dellat[k] = T(0.);
    if (ierr != 0) continue;
    xhe[k] = dellah[k] * mbdt + heo[k];
    xq[k] = fmx(T(1.e-16), (dellaq[k] + dellaqc[k]) * mbdt + qo[k]);
    dellat[k] = (T(1.) / K::cp) * (dellah[k] - K::xlv * (dellaq[k]));
    xt[k] = (-(dellaqc[k] * K::xlv / K::cp) + dellat[k]) * mbdt + tn[k];
    xt[k] = fmx(T(190.), xt[k]);
  }
  if (ierr == 0) {
    xhe[ktf] = heo[ktf];
    xq[ktf] = qo[ktf];
    xt[ktf] = tn[ktf];
  }
  cup_env(xz, xqes, xhe, xhes, xt, xq, po, ierr, ktf);
  cup_env_clev(xt, xqes, xq, xhe, xhes, xz, po, xqes_cup, xq_cup, xhe_cup, xhes_cup, xz_cup, po_cup, gamma_cup, xt_cup,
               psur, ierr, z1, ktf);
  for (int k = 1; k <= ktf; ++k) {
    xhc[k] = 0;
    xdby[k] = 0;
  }
  if (ierr == 0) {
    const T x_add = K::xlv * zqexec + K::cp * ztexec;
    const T xhkb = get_cloud_bc(xhe_cup, k22, true, x_add);
    for (int k = 1; k <= start_level - 1; ++k) xhc[k] = xhe_cup[k];
    xhc[start_level] = xhkb;
  }
  if (ierr == 0) {
    for (int k = 1; k <= ktf; ++k) xzu[k] = zuo[k];
    for (int k = start_level + 1; k <= ktop; ++k) {
      xhc[k] = (xhc[k - 1] * xzu[k - 1] - T(.5) * up_massdetro[k - 1] * xhc[k - 1] +
                up_massentro[k - 1] * xhe[k - 1]) /
               (xzu[k - 1] - T(.5) * up_massdetro[k - 1] + up_massentro[k - 1]);
      xdby[k] = xhc[k] - xhes_cup[k];
    }
    for (int k = ktop + 1; k <= ktf; ++k) {
      xhc[k] = xhes_cup[k];
      xdby[k] = 0;
      xzu[k] = 0;
    }
  }
  const T xaa0 = cup_up_aa0(xz, xzu, xdby, gamma_cup, xt_cup, kbcon, ktop, ierr, ktf);
  T xmb = T(0.);
  if (ierr == 0) {
    const T xmbmax = T(1.0);
    T xkshal = (xaa0 - aa1) / mbdt;
    if (xkshal <= T(0) && xkshal > T(-.01) * mbdt) xkshal = T(-.01) * mbdt;
    if (xkshal > T(0) && xkshal < T(1.e-2)) xkshal = T(1.e-2);
    T xff[4];
    xff[1] = fmx(T(0.), -((aa1 - aa0) / (xkshal * dtime)));
    xff[2] = T(.03) * zws;
    T blqe = T(0.);
    for (int k = 1; k <= kpbl; ++k) blqe = blqe + T(100.) * dhdt[k] * (po_cup[k] - po_cup[k + 1]) / K::g;
    const T trash = fmx((hc[kbcon] - he_cup[kbcon]), T(1.e1));
    xff[3] = fmx(T(0.), blqe / trash);
    xff[3] = fmn(xmbmax, xff[3]);
    xmb = (xff[1] + xff[2] + xff[3]) / T(3.);
    xmb = fmn(xmbmax, xmb);
    if (xmb <= T(0.)) ierr = 21;
  }
  if (ierr != 0) {
    k22 = 0;
    kbcon = 0;
    ktop = 0;
    xmb = T(0.);
    for (int k = 1; k <= kte; ++k) {
      outt[k] = 0; outq[k] = 0; outqc[k] = 0;
    }
  } else {
    xmb_out = xmb;
    pre = T(0.);
    for (int k = 2; k <= ktop; ++k) {
      outt[k] = dellat[k] * xmb;
      outq[k] = dellaq[k] * xmb;
      outqc[k] = dellaqc[k] * xmb;
      pre = pre + pwo[k] * xmb;
    }
  }
  return GfShallowOut{ierr, kbcon, ktop, k22};
}

}  // namespace cppwrf::phys::cu_gf
