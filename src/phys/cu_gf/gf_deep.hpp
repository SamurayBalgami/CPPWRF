#pragma once
// Grell-Freitas deep convection: WRF 4.6.1 module_cu_gf_deep.F CUP_gf for one column, as called by
// module_cu_gf_wrfdrv.F GFDRV (imid = 0, dicycle = 1, iversion = 1, nranflag = 0, csum = 0,
// irainevap = 0, autoconv = 1, aeroevap = 1). The second and third cup_kbcon calls (ierr2/ierr3) only
// feed ens_adj = 1 and are not ported.
#include <cstdio>

#include "phys/cu_gf/gf_common.hpp"

namespace cppwrf::phys::cu_gf {

template <typename T>
CPPWRF_HD inline void cup_dd_edt(int ierr, const T* us, const T* vs, const T* z, int ktop, int kbcon, T& edt,
                                 const T* p, T pwav, T pwev, T edtmax, T edtmin, T& edtc, int ktf) {
  edt = T(0.);
  T vws = T(0.), sdp = T(0.), vshear = T(0.);
  edtc = T(0.);
  for (int kk = 1; kk <= ktf - 1; ++kk) {
    if (ierr != 0) continue;
    if (kk <= (ktop < ktf ? ktop : ktf) && kk >= kbcon) {
      vws = vws + (std::fabs((us[kk + 1] - us[kk]) / (z[kk + 1] - z[kk])) +
                   std::fabs((vs[kk + 1] - vs[kk]) / (z[kk + 1] - z[kk]))) *
                      (p[kk] - p[kk + 1]);
      sdp = sdp + p[kk] - p[kk + 1];
    }
    if (kk == ktf - 1) vshear = T(1.e3) * vws / sdp;
  }
  if (ierr == 0) {
    T pef = (T(1.591) - T(.639) * vshear + T(.0953) * (vshear * vshear) - T(.00496) * (vshear * vshear * vshear));
    if (pef > T(0.9)) pef = T(0.9);
    if (pef < T(0.1)) pef = T(0.1);
    const T zkbc = z[kbcon] * T(3.281e-3);
    T prezk = T(.02);
    if (zkbc > T(3.))
      prezk = T(.96729352) +
              zkbc * (T(-.70034167) +
                      zkbc * (T(.162179896) + zkbc * (T(-1.2569798E-2) + zkbc * (T(4.2772E-4) - zkbc * T(5.44E-6)))));
    if (zkbc > T(25)) prezk = T(2.4);
    T pefb = T(1.) / (T(1.) + prezk);
    if (pefb > T(0.9)) pefb = T(0.9);
    if (pefb < T(0.1)) pefb = T(0.1);
    edt = T(1.) - T(.5) * (pefb + pef);
    const T einc = T(.2) * edt;
    edtc = edt - einc;
    edtc = -edtc * pwav / pwev;
    if (edtc > edtmax) edtc = edtmax;
    if (edtc < edtmin) edtc = edtmin;
  }
}

template <typename T>
CPPWRF_HD inline void cup_forcing_ens_3d(T& closure_n, int xland, T aa0, T aa1, T xaa0, T mbdt, T dtime, int ierr,
                                         T* xf_ens, T axx, const T* mconv_p, const T* p_cup, int ktop,
                                         const T* omeg, const T* zd, const T* zu, const T* pr_ens, T edt, int kbcon,
                                         int ichoice, int kte, T tau_ecmwf, T aa1_bl, T& xf_dicycle) {
  using K = GfK<T>;
  (void)p_cup; (void)ktop; (void)xland;
  const T mconv = *mconv_p;
  const T ens_adj = T(1.);
  T xff_dicycle = T(0.);
  T xff_ens3[maxens3 + 1] = {};
  if (ierr == 0) {
    const int kloc = maxloc_abs(zu, 1, kte);
    T a_ave = axx;
    a_ave = fmx(T(0.), a_ave);
    a_ave = fmn(a_ave, aa1);
    a_ave = fmx(T(0.), a_ave);
    (void)a_ave;
    const T xff0 = (aa1 - aa0) / dtime;
    xff_ens3[1] = fmx(T(0.), (aa1 - aa0) / dtime);
    xff_ens3[2] = fmx(T(0.), (aa1 - aa0) / dtime);
    xff_ens3[3] = fmx(T(0.), (aa1 - aa0) / dtime);
    xff_ens3[16] = fmx(T(0.), (aa1 - aa0) / dtime);
    T xomg = T(0.);
    int kk = 0;
    for (int k = kbcon - 1; k <= kbcon + 1; ++k)
      if (zu[k] > T(0.)) {
        xomg = xomg - omeg[k] / T(9.81) / fmx(T(0.5), (T(1.) - edt * zd[k] / zu[k]));
        kk = kk + 1;
      }
    if (kk > 0) xff_ens3[4] = xomg / T(kk);
    xff_ens3[4] = K::betajb * xff_ens3[4];
    xff_ens3[5] = xff_ens3[4];
    xff_ens3[6] = xff_ens3[4];
    if (xff_ens3[4] < T(0.)) xff_ens3[4] = T(0.);
    if (xff_ens3[5] < T(0.)) xff_ens3[5] = T(0.);
    if (xff_ens3[6] < T(0.)) xff_ens3[6] = T(0.);
    xff_ens3[14] = K::betajb * xff_ens3[4];
    const T mc = mconv / fmx(T(0.5), (T(1.) - edt * zd[kbcon] / zu[kloc]));
    xff_ens3[7] = mc;
    xff_ens3[8] = mc;
    xff_ens3[9] = mc;
    xff_ens3[15] = mc;
    xff_ens3[10] = aa1 / tau_ecmwf;
    xff_ens3[11] = aa1 / tau_ecmwf;
    xff_ens3[12] = aa1 / tau_ecmwf;
    xff_ens3[13] = (aa1) / tau_ecmwf;
    xff_dicycle = fmx(T(0.), aa1_bl / tau_ecmwf);  // dicycle == 1
    if (ichoice == 0) {
      if (xff0 < T(0.)) {
        xff_ens3[1] = 0; xff_ens3[2] = 0; xff_ens3[3] = 0; xff_ens3[10] = 0;
        xff_ens3[11] = 0; xff_ens3[12] = 0; xff_ens3[13] = 0; xff_ens3[16] = 0;
        closure_n = T(12.);
      }
    }
    T xk = (xaa0 - aa1) / mbdt;
    if (xk <= T(0) && xk > T(-.01) * mbdt) xk = T(-.01) * mbdt;
    if (xk > T(0) && xk < T(1.e-2)) xk = T(1.e-2);
    (void)ens_adj;  // xland < 0.1 with ierr2/ierr3 > 0 multiplies by ens_adj = 1
    if (xk < T(0.)) {
      if (xff_ens3[1] > T(0)) xf_ens[1] = fmx(T(0.), -xff_ens3[1] / xk);
      if (xff_ens3[2] > T(0)) xf_ens[2] = fmx(T(0.), -xff_ens3[2] / xk);
      if (xff_ens3[3] > T(0)) xf_ens[3] = fmx(T(0.), -xff_ens3[3] / xk);
      if (xff_ens3[16] > T(0)) xf_ens[16] = fmx(T(0.), -xff_ens3[16] / xk);
      for (int n : {1, 2, 3, 16}) xf_ens[n] = xf_ens[n] + xf_ens[n] * T(0.);  // rand_clos = 0
    } else {
      xff_ens3[1] = 0; xff_ens3[2] = 0; xff_ens3[3] = 0; xff_ens3[16] = 0;
    }
    xf_ens[4] = fmx(T(0.), xff_ens3[4]);
    xf_ens[5] = fmx(T(0.), xff_ens3[5]);
    xf_ens[6] = fmx(T(0.), xff_ens3[6]);
    xf_ens[14] = fmx(T(0.), xff_ens3[14]);
    T a1 = fmx(T(1.e-5), pr_ens[7]);
    xf_ens[7] = fmx(T(0.), xff_ens3[7] / a1);
    a1 = fmx(T(1.e-5), pr_ens[8]);
    xf_ens[8] = fmx(T(0.), xff_ens3[8] / a1);
    a1 = fmx(T(1.e-5), pr_ens[9]);
    xf_ens[9] = fmx(T(0.), xff_ens3[9] / a1);
    a1 = fmx(T(1.e-3), pr_ens[15]);
    xf_ens[15] = fmx(T(0.), xff_ens3[15] / a1);
    for (int n : {4, 5, 6, 14, 7, 8, 9, 15}) xf_ens[n] = xf_ens[n] + xf_ens[n] * T(0.);
    if (xk < T(0.)) {
      for (int n = 10; n <= 13; ++n) xf_ens[n] = fmx(T(0.), -xff_ens3[n] / xk);
      for (int n = 10; n <= 13; ++n) xf_ens[n] = xf_ens[n] + xf_ens[n] * T(0.);
    } else {
      for (int n = 10; n <= 13; ++n) xf_ens[n] = T(0.);
    }
    if (xk < T(0.)) xf_dicycle = fmx(T(0.), -xff_dicycle / xk);
    else xf_dicycle = T(0.);
    if (ichoice >= 1) {
      const T v = xf_ens[ichoice];
      for (int n = 1; n <= 16; ++n) xf_ens[n] = v;
    }
  } else if (ierr != 20 && ierr != 0) {
    for (int n = 1; n <= maxens3; ++n) {
      xf_ens[n] = T(0.);
      xf_dicycle = T(0.);
    }
  }
}

// cup_output_ens_3d, imid = 0, dicycle = 1
template <typename T>
CPPWRF_HD inline void cup_output_ens_3d(T* xf_ens, int& ierr, const T* dellat, const T* dellaq, const T* dellaqc,
                                        T* outtem, T* outq, T* outqc, T& pre, const T* pw, T& xmb, int ktop, T edt,
                                        const T* pwd, const T* p_cup, const T* pr_ens, T sig, T closure_n,
                                        int kte, T xf_dicycle) {
  using K = GfK<T>;
  for (int k = 1; k <= kte; ++k) {
    outtem[k] = 0;
    outq[k] = 0;
    outqc[k] = 0;
  }
  pre = T(0.);
  xmb = T(0.);
  if (ierr == 0)
    for (int n = 1; n <= maxens3; ++n)
      if (pr_ens[n] <= T(0.)) xf_ens[n] = T(0.);
  if (ierr == 0) {
    int k = 0;
    T xmb_ave = T(0.);
    for (int n = 1; n <= maxens3; ++n) {
      k = k + 1;
      xmb_ave = xmb_ave + xf_ens[n];
    }
    xmb_ave = xmb_ave / T(k);
    xmb_ave = fmn(xmb_ave, xmb_ave - xf_dicycle);
    xmb_ave = fmx(T(0.), xmb_ave);
    const T clos_wei = T(16.) / fmx(T(1.), closure_n);
    xmb_ave = fmn(xmb_ave, T(100.));
    xmb = clos_wei * sig * xmb_ave;
    if (xmb < T(1.e-16)) ierr = 19;
  }
  if (ierr == 0) {
    T pwtot = T(0.);
    for (int k = 1; k <= ktop; ++k) pwtot = pwtot + pw[k];
    for (int k = 1; k <= ktop; ++k) {
      const T dp = T(100.) * (p_cup[k] - p_cup[k + 1]) / K::g;
      const T dtt = dellat[k];
      const T dtq = dellaq[k];
      T dtpwd = -pwd[k] * edt;
      T dtqc = dellaqc[k] * dp - dtpwd;
      if (dtqc < T(0.)) {
        dtpwd = dtpwd - dellaqc[k] * dp;
        dtqc = T(0.);
      } else {
        dtpwd = T(0.);
        dtqc = dtqc / dp;
      }
      outtem[k] = xmb * dtt;
      outq[k] = xmb * dtq;
      outqc[k] = xmb * dtqc;
      pre = pre - xmb * dtpwd;
    }
    pre = -pre + xmb * pwtot;
  }
}

// neg_check (the moisture part after its RETURN is dead code)
template <typename T>
CPPWRF_HD inline void neg_check(bool shallow, T* outq, T* outt, T* outu, T* outv, T* outqc, T& pret, int ktf) {
  T thresh = T(300.01), names = T(1.);
  if (shallow) {
    thresh = T(148.01);
    names = T(2.);
  }
  T qmemf = T(1.);
  for (int k = 1; k <= ktf; ++k) {
    const T qmem = (outt[k]) * T(86400.);
    if (qmem > thresh) {
      const T qmem2 = thresh / qmem;
      qmemf = fmn(qmemf, qmem2);
    }
    if (qmem < T(-.5) * thresh * names) {
      const T qmem2 = T(-.5) * names * thresh / qmem;
      qmemf = fmn(qmemf, qmem2);
    }
  }
  for (int k = 1; k <= ktf; ++k) {
    outq[k] = outq[k] * qmemf;
    outt[k] = outt[k] * qmemf;
    outu[k] = outu[k] * qmemf;
    outv[k] = outv[k] * qmemf;
    outqc[k] = outqc[k] * qmemf;
  }
  pret = pret * qmemf;
}

// Inputs of CUP_gf for one column (1-based level arrays, levels 1..nz)
template <typename T>
struct GfColumnIn {
  int nz;
  const T *t, *q, *tn, *qo, *po, *us, *vs, *rho, *omeg, *zo, *dhdt;
  T z1, psur, hfx, qfx, xland, dx, dtime;
  int kpbl, ichoice;
  T xmbs_in;  // shallow mass flux (dicycle 2 only; unused)
  int dbg = 0;  // debug print of intermediate values (host only)
};

struct GfDeepOut {
  int ierr, kbcon, ktop, k22, jmin;
};

// CUP_gf deep (imid = 0). outt, outq, outqc, outu, outv, cupclw, cnvwt zeroed by the caller; pre = pret.
template <typename T>
CPPWRF_HD inline GfDeepOut cup_gf_deep(const GfColumnIn<T>& in, T* outu, T* outv, T* outt, T* outq, T* outqc,
                                       T* cupclw, T* cnvwt, T& pre, T& xmb_out) {
  using K = GfK<T>;
  const int kte = in.nz, ktf = in.nz;
  const T* t = in.t; const T* q = in.q; const T* tn = in.tn; const T* qo = in.qo; const T* po = in.po;
  const T* us = in.us; const T* vs = in.vs; const T* rho = in.rho; const T* omeg = in.omeg;
  const T z1 = in.z1, psur = in.psur, dx = in.dx, dtime = in.dtime;
  const int kpbl = in.kpbl;
  const int imid = 0, csum = 0;
  int ierr = 0, kbcon = 0, ktop = 0, k22 = 0, jmin = 0, kdet = 0, kbmax = 1, kstabi = 0, kzdown = 0;
  int ktopdby = 0, pmin_lev = 0;

  T zo[KM] = {}, z[KM] = {}, xz[KM] = {}, zuo[KM] = {}, zdo[KM] = {}, zu[KM] = {}, xzu[KM] = {};
  T qes[KM] = {}, he[KM] = {}, hes[KM] = {}, qeso[KM] = {}, heo[KM] = {}, heso[KM] = {};
  T qes_cup[KM] = {}, q_cup[KM] = {}, he_cup[KM] = {}, hes_cup[KM] = {}, z_cup[KM] = {}, p_cup[KM] = {};
  T gamma_cup[KM] = {}, t_cup[KM] = {};
  T qeso_cup[KM] = {}, qo_cup[KM] = {}, heo_cup[KM] = {}, heso_cup[KM] = {}, zo_cup[KM] = {}, po_cup[KM] = {};
  T gammao_cup[KM] = {}, tn_cup[KM] = {};
  T u_cup[KM] = {}, v_cup[KM] = {}, uc[KM] = {}, vc[KM] = {}, ucd[KM] = {}, vcd[KM] = {};
  T hc[KM] = {}, dby[KM] = {}, hco[KM] = {}, dbyo[KM] = {}, dbyt[KM] = {};
  T cd[KM] = {}, cdd[KM] = {}, c1d[KM] = {}, entr_rate_2d[KM] = {}, mentrd_rate_2d[KM] = {};
  T up_massentr[KM] = {}, up_massdetr[KM] = {}, up_massentro[KM] = {}, up_massdetro[KM] = {};
  T up_massentru[KM] = {}, up_massdetru[KM] = {};
  T dd_massentro[KM] = {}, dd_massdetro[KM] = {}, dd_massentru[KM] = {}, dd_massdetru[KM] = {};
  T hcdo[KM] = {}, qcdo[KM] = {}, qrcdo[KM] = {}, pwdo[KM] = {}, dbydo[KM] = {};
  T qco[KM] = {}, qrco[KM] = {}, pwo[KM] = {}, clw_all[KM] = {};
  T dellu[KM] = {}, dellv[KM] = {}, dellah[KM] = {}, dellaq[KM] = {}, dellat[KM] = {}, dellaqc[KM] = {};
  T xhe[KM] = {}, xq[KM] = {}, xt[KM] = {}, xqes[KM] = {}, xhes[KM] = {};
  T xqes_cup[KM] = {}, xq_cup[KM] = {}, xhe_cup[KM] = {}, xhes_cup[KM] = {}, xz_cup[KM] = {}, xt_cup[KM] = {};
  T xhc[KM] = {}, xdby[KM] = {};
  T xf_ens[maxens3 + 1] = {}, pr_ens[maxens3 + 1] = {};

  for (int k = 1; k <= kte; ++k) zo[k] = in.zo[k];
  const T flux_tun = K::fluxtune;
  const T pmin = T(150.);
  const T lambau = T(2.);
  const T pgcon = T(0.);
  T ztexec = T(0.), zqexec = T(0.), zws = T(0.);
  {
    const T buo_flux = (in.hfx / K::cp + T(0.608) * t[1] * in.qfx / K::xlv) / rho[1];
    const T pgeoh = zo[2] * K::g;
    zws = fmx(T(0.), flux_tun * T(0.41) * buo_flux * zo[2] * K::g / t[1]);
    if (zws > std::numeric_limits<T>::min()) {
      (void)pgeoh;
      zws = T(1.2) * std::pow(zws, T(.3333));
      ztexec = fmx(flux_tun * in.hfx / (rho[1] * zws * K::cp), T(0.0));
      zqexec = fmx(flux_tun * in.qfx / K::xlv / (rho[1] * zws), T(0.));
    }
    zws = fmx(T(0.), flux_tun * T(0.41) * buo_flux * zo[kpbl] * K::g / t[kpbl]);
    zws = T(1.2) * std::pow(zws, T(.3333));
    zws = zws * rho[kpbl];
  }
  const T cap_maxs = T(75.);
  T edto = T(0.);
  T closure_n = T(16.);
  xmb_out = T(0.);
  T cap_max = cap_maxs;
  const T cap_max_increment = T(20.);
  int xland1 = int(in.xland + T(.0001));
  if (in.xland > T(1.5) || in.xland < T(0.5)) {
    xland1 = 0;
  } else {
    if (ztexec > T(0.)) cap_max = cap_max + T(25.);
    if (ztexec < T(0.)) cap_max = cap_max - T(25.);
  }
  int start_level = kte;
  T entr_rate = T(7.e-5) - fmn(T(20.), T(csum)) * T(3.e-6);
  if (xland1 == 0) entr_rate = T(7.e-5);
  T radius = T(.2) / entr_rate;
  T frh = fmn(T(1.), T(3.14) * radius * radius / dx / dx);
  if (frh > K::frh_thresh) {
    frh = K::frh_thresh;
    radius = std::sqrt(frh * dx * dx / T(3.14));
    entr_rate = T(.2) / radius;
  }
  const T sig = (T(1.) - frh) * (T(1.) - frh);
  const T sig_thresh = (T(1.) - K::frh_thresh) * (T(1.) - K::frh_thresh);
  for (int k = 1; k <= ktf; ++k) {
    cnvwt[k] = 0; zuo[k] = 0; zdo[k] = 0; z[k] = zo[k]; xz[k] = zo[k]; cupclw[k] = 0;
    cd[k] = T(1.e-9); cdd[k] = T(1.e-9); hcdo[k] = 0; qrcdo[k] = 0; dellaqc[k] = 0;
  }
  T edtmax = T(1.);
  const T edtmin = T(.1);
  const T depth_min = T(1000.);
  T aa0 = T(0.), aa1 = T(0.), edt = T(0.);
  const int kstabm = ktf - 1;
  const T zkbmax = T(4000.), zcutdown = T(4000.), z_detr = T(1000.);

  cup_env(z, qes, he, hes, t, q, po, ierr, ktf);
  cup_env(zo, qeso, heo, heso, tn, qo, po, ierr, ktf);
  cup_env_clev(t, qes, q, he, hes, z, po, qes_cup, q_cup, he_cup, hes_cup, z_cup, p_cup, gamma_cup, t_cup, psur, ierr,
               z1, ktf);
  cup_env_clev(tn, qeso, qo, heo, heso, zo, po, qeso_cup, qo_cup, heo_cup, heso_cup, zo_cup, po_cup, gammao_cup,
               tn_cup, psur, ierr, z1, ktf);
  if (ierr == 0) {
    u_cup[1] = us[1];
    v_cup[1] = vs[1];
    for (int k = 2; k <= ktf; ++k) {
      u_cup[k] = T(.5) * (us[k - 1] + us[k]);
      v_cup[k] = T(.5) * (vs[k - 1] + vs[k]);
    }
  }
  if (ierr == 0) {
    for (int k = 1; k <= ktf; ++k)
      if (zo_cup[k] > zkbmax + z1) {
        kbmax = k;
        break;
      }
    for (int k = 1; k <= ktf; ++k)
      if (zo_cup[k] > z_detr + z1) {
        kdet = k;
        break;
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
  T hkb = T(0), hkbo = T(0);
  if (ierr == 0) {
    const T x_add = K::xlv * zqexec + K::cp * ztexec;
    hkb = get_cloud_bc(he_cup, k22, true, x_add);
    hkbo = get_cloud_bc(heo_cup, k22, true, x_add);
  }
  cup_kbcon(cap_max_increment, 1, k22, kbcon, heo_cup, heso_cup, hkbo, ierr, kbmax, po_cup, cap_max, ztexec, zqexec,
            z_cup, entr_rate, heo, imid);
  kstabi = cup_minimi(heso_cup, kbcon, kstabm, ierr);
  if (ierr == 0) {
    frh = fmn(qo_cup[kbcon] / qeso_cup[kbcon], T(1.));
    if (frh >= K::rh_thresh && sig <= sig_thresh) {
      ierr = 231;
    } else {
      const T x_add0 = T(0.);
      for (int k = kbcon + 1; k <= ktf; ++k)
        if (po[kbcon] - po[k] > pmin + x_add0) {
          pmin_lev = k;
          break;
        }
      start_level = k22;
      const T x_add = K::xlv * zqexec + K::cp * ztexec;
      hkb = get_cloud_bc(he_cup, k22, true, x_add);
    }
  }
  (void)pmin_lev;
  if (kstabi < kbcon) {
    kbcon = 1;
    ierr = 42;
  }
  for (int k = 1; k <= ktf; ++k) entr_rate_2d[k] = entr_rate;
  if (ierr == 0) {
    kbcon = kbcon > 2 ? kbcon : 2;
    for (int k = 1; k <= ktf; ++k) {
      frh = fmn(qo_cup[k] / qeso_cup[k], T(1.));
      entr_rate_2d[k] = entr_rate * (T(1.3) - frh);
    }
  }
  // rates_up_pdf 'deep' (its kpbl argument is kbcon)
  rates_up_pdf(false, ktop, ierr, po_cup, entr_rate_2d, hkbo, heo, heso_cup, zo_cup, kstabi, k22, kbcon, kte, ktf,
               zuo, kbcon, ktopdby, csum);
  if (ierr == 0) {
    if (k22 > 1)
      for (int k = 1; k <= k22 - 1; ++k) {
        zuo[k] = 0; zu[k] = 0; xzu[k] = 0;
      }
    for (int k = k22; k <= ktop; ++k) {
      xzu[k] = zuo[k];
      zu[k] = zuo[k];
    }
    for (int k = ktop + 1; k <= kte; ++k) {
      zuo[k] = 0; zu[k] = 0; xzu[k] = 0;
    }
  }
  get_lateral_massflux(ktf, kte, ierr, ktop, zo_cup, zuo, cd, entr_rate_2d, up_massentro, up_massdetro, up_massentr,
                       up_massdetr, k22, true, up_massentru, up_massdetru, lambau);
  for (int k = 1; k <= ktf; ++k) {
    uc[k] = 0; vc[k] = 0; hc[k] = 0; dby[k] = 0; hco[k] = 0; dbyo[k] = 0;
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
  int ktopkeep = 0;
  for (int k = 1; k <= kte; ++k) dbyt[k] = 0;
  if (ierr == 0) {
    ktopkeep = ktop;
    for (int k = start_level + 1; k <= ktop; ++k) {
      const T denom = zuo[k - 1] - T(.5) * up_massdetro[k - 1] + up_massentro[k - 1];
      if (denom < T(1.e-8)) {
        ierr = 51;
        break;
      }
      hc[k] = (hc[k - 1] * zu[k - 1] - T(.5) * up_massdetr[k - 1] * hc[k - 1] + up_massentr[k - 1] * he[k - 1]) /
              (zu[k - 1] - T(.5) * up_massdetr[k - 1] + up_massentr[k - 1]);
      uc[k] = (uc[k - 1] * zu[k - 1] - T(.5) * up_massdetru[k - 1] * uc[k - 1] + up_massentru[k - 1] * us[k - 1] -
               pgcon * T(.5) * (zu[k] + zu[k - 1]) * (u_cup[k] - u_cup[k - 1])) /
              (zu[k - 1] - T(.5) * up_massdetru[k - 1] + up_massentru[k - 1]);
      vc[k] = (vc[k - 1] * zu[k - 1] - T(.5) * up_massdetru[k - 1] * vc[k - 1] + up_massentru[k - 1] * vs[k - 1] -
               pgcon * T(.5) * (zu[k] + zu[k - 1]) * (v_cup[k] - v_cup[k - 1])) /
              (zu[k - 1] - T(.5) * up_massdetru[k - 1] + up_massentru[k - 1]);
      dby[k] = hc[k] - hes_cup[k];
      hco[k] = (hco[k - 1] * zuo[k - 1] - T(.5) * up_massdetro[k - 1] * hco[k - 1] +
                up_massentro[k - 1] * heo[k - 1]) /
               (zuo[k - 1] - T(.5) * up_massdetro[k - 1] + up_massentro[k - 1]);
      dbyo[k] = hco[k] - heso_cup[k];
      const T dz = zo_cup[k + 1] - zo_cup[k];
      dbyt[k] = dbyt[k - 1] + dbyo[k] * dz;
    }
    for (int k = ktop - 1; k >= kbcon; --k)
      if (dbyo[k] > T(0.)) {
        ktopkeep = k + 1;
        break;
      }
    ktop = ktopkeep;
  }
  if (ierr == 0)
    for (int k = ktop + 1; k <= ktf; ++k) {
      hc[k] = hes_cup[k]; uc[k] = u_cup[k]; vc[k] = v_cup[k]; hco[k] = heso_cup[k];
      dby[k] = 0; dbyo[k] = 0; zu[k] = 0; zuo[k] = 0; cd[k] = 0; entr_rate_2d[k] = 0;
      up_massentr[k] = 0; up_massdetr[k] = 0; up_massentro[k] = 0; up_massdetro[k] = 0;
    }
  if (ierr == 0)
    if (ktop < kbcon + 2) {
      ierr = 5;
      ktop = 0;
    }
  kzdown = 0;
  if (ierr == 0) {
    T zktop = (zo_cup[ktop] - z1) * T(.6);
    zktop = fmn(zktop + z1, zcutdown + z1);
    for (int k = 1; k <= ktf; ++k)
      if (zo_cup[k] > zktop) {
        kzdown = k;
        kzdown = kzdown < kstabi - 1 ? kzdown : kstabi - 1;
        break;
      }
  }
  jmin = cup_minimi(heso_cup, k22, kzdown, ierr);
  if (ierr == 0) {
    int jmini = jmin;
    bool keep_going = true;
    while (keep_going) {
      keep_going = false;
      if (jmini - 1 < kdet) kdet = jmini - 1;
      if (jmini >= ktop - 1) jmini = ktop - 2;
      const int ki = jmini;
      hcdo[ki] = heso_cup[ki];
      T dh = T(0.);
      for (int k = ki - 1; k >= 1; --k) {
        hcdo[k] = heso_cup[jmini];
        const T dz = zo_cup[k + 1] - zo_cup[k];
        dh = dh + dz * (hcdo[k] - heso_cup[k]);
        if (dh > T(0.)) {
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
    if (jmin - 1 < kdet) kdet = jmin - 1;
    if (-zo_cup[kbcon] + zo_cup[ktop] < depth_min) ierr = 6;
  }
  for (int k = 1; k <= ktf; ++k) {
    zdo[k] = 0; cdd[k] = 0; dd_massentro[k] = 0; dd_massdetro[k] = 0; dd_massentru[k] = 0; dd_massdetru[k] = 0;
    hcdo[k] = heso_cup[k]; ucd[k] = u_cup[k]; vcd[k] = v_cup[k]; dbydo[k] = 0; mentrd_rate_2d[k] = entr_rate;
  }
  T bud = T(0.);
  {
    const T beta = fmx(T(.02), T(.05) - T(csum) * T(.0015));
    (void)beta;
    if (xland1 == 0) edtmax = fmx(T(0.1), T(.4) - T(csum) * T(.015));
    bud = T(0.);
    bool skip = false;
    if (ierr == 0) {
      for (int k = 1; k <= jmin; ++k) cdd[k] = T(1.e-9);
      cdd[jmin] = T(0.);
      for (int k = 1; k <= kte; ++k) {
        dd_massdetro[k] = 0;
        dd_massentro[k] = 0;
      }
      get_zu_zd_pdf_fim(po_cup, T(0.), 2, kdet, jmin, zdo, kte, ktf, kpbl, csum);
      if (zdo[jmin] < T(1.e-8)) {
        zdo[jmin] = T(0.);
        jmin = jmin - 1;
        if (zdo[jmin] < T(1.e-8)) {
          ierr = 876;
          skip = true;
        }
      }
      if (!skip) {
        const int kmz = maxloc_abs(zdo, 1, kte);
        for (int ki = jmin; ki >= kmz; --ki) {
          const T dzo = zo_cup[ki + 1] - zo_cup[ki];
          dd_massdetro[ki] = cdd[ki] * dzo * zdo[ki + 1];
          dd_massentro[ki] = zdo[ki] - zdo[ki + 1] + dd_massdetro[ki];
          if (dd_massentro[ki] < T(0.)) {
            dd_massentro[ki] = T(0.);
            dd_massdetro[ki] = zdo[ki + 1] - zdo[ki];
            if (zdo[ki + 1] > T(0.)) cdd[ki] = dd_massdetro[ki] / (dzo * zdo[ki + 1]);
          }
          if (zdo[ki + 1] > T(0.)) mentrd_rate_2d[ki] = dd_massentro[ki] / (dzo * zdo[ki + 1]);
        }
        mentrd_rate_2d[1] = T(0.);
        for (int ki = kmz - 1; ki >= 1; --ki) {
          const T dzo = zo_cup[ki + 1] - zo_cup[ki];
          dd_massentro[ki] = mentrd_rate_2d[ki] * dzo * zdo[ki + 1];
          dd_massdetro[ki] = zdo[ki + 1] + dd_massentro[ki] - zdo[ki];
          if (dd_massdetro[ki] < T(0.)) {
            dd_massdetro[ki] = T(0.);
            dd_massentro[ki] = zdo[ki] - zdo[ki + 1];
            if (zdo[ki + 1] > T(0.)) mentrd_rate_2d[ki] = dd_massentro[ki] / (dzo * zdo[ki + 1]);
          }
          if (zdo[ki + 1] > T(0.)) cdd[ki] = dd_massdetro[ki] / (dzo * zdo[ki + 1]);
        }
        for (int k = kbcon + 1; k <= ktop - 1; ++k) c1d[k] = K::c1;  // the parabola is overwritten by c1
        for (int k = 2; k <= jmin + 1; ++k) {
          dd_massentru[k - 1] = dd_massentro[k - 1] + lambau * dd_massdetro[k - 1];
          dd_massdetru[k - 1] = dd_massdetro[k - 1] + lambau * dd_massdetro[k - 1];
        }
        dbydo[jmin] = hcdo[jmin] - heso_cup[jmin];
        bud = dbydo[jmin] * (zo_cup[jmin + 1] - zo_cup[jmin]);
        for (int ki = jmin; ki >= 1; --ki) {
          const T dzo = zo_cup[ki + 1] - zo_cup[ki];
          const T h_entr = T(.5) * (heo[ki] + T(.5) * (hco[ki] + hco[ki + 1]));
          ucd[ki] = (ucd[ki + 1] * zdo[ki + 1] - T(.5) * dd_massdetru[ki] * ucd[ki + 1] + dd_massentru[ki] * us[ki] -
                     pgcon * zdo[ki + 1] * (us[ki + 1] - us[ki])) /
                    (zdo[ki + 1] - T(.5) * dd_massdetru[ki] + dd_massentru[ki]);
          vcd[ki] = (vcd[ki + 1] * zdo[ki + 1] - T(.5) * dd_massdetru[ki] * vcd[ki + 1] + dd_massentru[ki] * vs[ki] -
                     pgcon * zdo[ki + 1] * (vs[ki + 1] - vs[ki])) /
                    (zdo[ki + 1] - T(.5) * dd_massdetru[ki] + dd_massentru[ki]);
          hcdo[ki] = (hcdo[ki + 1] * zdo[ki + 1] - T(.5) * dd_massdetro[ki] * hcdo[ki + 1] +
                      dd_massentro[ki] * h_entr) /
                     (zdo[ki + 1] - T(.5) * dd_massdetro[ki] + dd_massentro[ki]);
          dbydo[ki] = hcdo[ki] - heso_cup[ki];
          bud = bud + dbydo[ki] * dzo;
        }
      }
    }
    if (!skip && bud > T(0)) ierr = 7;
  }
  T pwevo = T(0), bu = T(0);
  cup_dd_moisture(zdo, hcdo, heso_cup, qcdo, qeso_cup, pwdo, qo_cup, zo_cup, dd_massentro, dd_massdetro, jmin, ierr,
                  gammao_cup, pwevo, bu, qrcdo, qo, 1, ktf);
  T pwavo = T(0), psum = T(0), psumh = T(0);
  cup_up_moisture(ierr, zo_cup, qco, qrco, pwo, pwavo, kbcon, ktop, dbyo, clw_all, qo, gammao_cup, zuo, qeso_cup, k22,
                  qo_cup, c1d, tn_cup, up_massentr, up_massdetr, psum, psumh, ktf);
  if (ierr == 0)
    for (int k = 2; k <= ktop; ++k) {
      const T dp = T(100.) * (po_cup[1] - po_cup[2]);
      cupclw[k] = qrco[k];
      cnvwt[k] = zuo[k] * cupclw[k] * K::g / dp;
    }
  aa0 = cup_up_aa0(z, zu, dby, gamma_cup, t_cup, kbcon, ktop, ierr, ktf);
  aa1 = cup_up_aa0(zo, zuo, dbyo, gammao_cup, tn_cup, kbcon, ktop, ierr, ktf);
  if (ierr == 0 && aa1 == T(0.)) ierr = 17;
  T aa1_bl = T(0.0), tau_ecmwf = T(0.), wmean = T(0), tau_bl = T(0.);
  if (ierr == 0) {
    wmean = T(7.0);
    tau_ecmwf = (zo_cup[ktopdby] - zo_cup[kbcon]) / wmean;
    tau_ecmwf = tau_ecmwf * (T(1.0061) + T(1.23E-2) * (dx / T(1000.)));
  }
  // dicycle = 1, iversion = 1
  if (ierr == 0) {
    if (xland1 == 0) {
      const T umean = T(2.0) + std::sqrt(T(2.0) * (us[1] * us[1] + vs[1] * vs[1] + us[kbcon] * us[kbcon] +
                                                   vs[kbcon] * vs[kbcon]));
      tau_bl = (zo_cup[kbcon] - z1) / umean;
    } else {
      tau_bl = (zo_cup[ktopdby] - zo_cup[kbcon]) / wmean;
    }
  }
  {
    const T t_star = T(4.);
    // cup_up_aa1bl(aa1_bl, t, tn, q, qo, dtime, zo_cup, ...)
    aa1_bl = T(0.);
    for (int k = 2; k <= ktf; ++k) {
      if (ierr != 0) continue;
      if (k > kbcon) continue;
      const T dz = zo_cup[k] - zo_cup[k - 1];
      const T da = dz * T(9.81) * (tn[k] - t[k] + T(0.608) * (qo[k] - q[k])) / dtime;
      aa1_bl = aa1_bl + da;
    }
    if (ierr == 0) {
      const int kp = kte < kpbl + 1 ? kte : kpbl + 1;
      if (zo_cup[kbcon] - z1 > zo[kp]) aa1_bl = T(0.0);
      else aa1_bl = fmx(T(0.), aa1_bl / t_star * tau_bl);
    }
  }
  const T axx = aa1;
  T edtc = T(0);
  cup_dd_edt(ierr, us, vs, zo, ktop, kbcon, edt, po, pwavo, pwevo, edtmax, edtmin, edtc, ktf);
  if (ierr == 0) edto = edtc;
  for (int k = 1; k <= kte; ++k) {
    dellu[k] = 0; dellv[k] = 0; dellah[k] = 0; dellat[k] = 0; dellaq[k] = 0; dellaqc[k] = 0;
  }
  if (ierr == 0) {
    T dp = T(100.) * (po_cup[1] - po_cup[2]);
    dellu[1] = K::pgcd * (edto * zdo[2] * ucd[2] - edto * zdo[2] * u_cup[2]) * K::g / dp;
    dellv[1] = K::pgcd * (edto * zdo[2] * vcd[2] - edto * zdo[2] * v_cup[2]) * K::g / dp;
    for (int k = 2; k <= ktop; ++k) {
      dp = T(100.) * (po_cup[k] - po_cup[k + 1]);
      dellu[k] = -((zuo[k + 1] * (uc[k + 1] - u_cup[k + 1]) - zuo[k] * (uc[k] - u_cup[k])) * K::g / dp) +
                 (zdo[k + 1] * (ucd[k + 1] - u_cup[k + 1]) - zdo[k] * (ucd[k] - u_cup[k])) * K::g / dp * edto *
                     K::pgcd;
      dellv[k] = -((zuo[k + 1] * (vc[k + 1] - v_cup[k + 1]) - zuo[k] * (vc[k] - v_cup[k])) * K::g / dp) +
                 (zdo[k + 1] * (vcd[k + 1] - v_cup[k + 1]) - zdo[k] * (vcd[k] - v_cup[k])) * K::g / dp * edto *
                     K::pgcd;
    }
  }
  if (ierr == 0) {
    T dp = T(100.) * (po_cup[1] - po_cup[2]);
    dellah[1] = (edto * zdo[2] * hcdo[2] - edto * zdo[2] * heo_cup[2]) * K::g / dp;
    dellaq[1] = (edto * zdo[2] * qcdo[2] - edto * zdo[2] * qo_cup[2]) * K::g / dp;
    T g_rain = T(0.5) * (pwo[1] + pwo[2]) * K::g / dp;
    T e_dn = -(T(0.5) * (pwdo[1] + pwdo[2]) * K::g / dp * edto);
    dellaq[1] = dellaq[1] + e_dn - g_rain;
    for (int k = 2; k <= ktop; ++k) {
      dp = T(100.) * (po_cup[k] - po_cup[k + 1]);
      dellah[k] = -((zuo[k + 1] * (hco[k + 1] - heo_cup[k + 1]) - zuo[k] * (hco[k] - heo_cup[k])) * K::g / dp) +
                  (zdo[k + 1] * (hcdo[k + 1] - heo_cup[k + 1]) - zdo[k] * (hcdo[k] - heo_cup[k])) * K::g / dp * edto;
      const T detup = up_massdetro[k];
      const T dz = zo_cup[k] - zo_cup[k - 1];
      if (k < ktop) dellaqc[k] = zuo[k] * c1d[k] * qrco[k] * dz / dp * K::g;
      if (k == ktop) dellaqc[k] = detup * T(0.5) * (qrco[k + 1] + qrco[k]) * K::g / dp;
      g_rain = T(0.5) * (pwo[k] + pwo[k + 1]) * K::g / dp;
      e_dn = -(T(0.5) * (pwdo[k] + pwdo[k + 1]) * K::g / dp * edto);
      const T c_up = dellaqc[k] + (zuo[k + 1] * qrco[k + 1] - zuo[k] * qrco[k]) * K::g / dp + g_rain;
      dellaq[k] = -((zuo[k + 1] * (qco[k + 1] - qo_cup[k + 1]) - zuo[k] * (qco[k] - qo_cup[k])) * K::g / dp) +
                  (zdo[k + 1] * (qcdo[k + 1] - qo_cup[k + 1]) - zdo[k] * (qcdo[k] - qo_cup[k])) * K::g / dp * edto -
                  c_up + e_dn;
    }
  }
  const T mbdt = T(.1);
  T xaa0_ens = T(0.);
  if (ierr == 0)
    for (int k = 1; k <= ktf; ++k) {
      xhe[k] = dellah[k] * mbdt + heo[k];
      xq[k] = fmx(T(1.e-16), dellaq[k] * mbdt + qo[k]);
      dellat[k] = (T(1.) / K::cp) * (dellah[k] - K::xlv * dellaq[k]);
      xt[k] = dellat[k] * mbdt + tn[k];
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
    for (int k = start_level + 1; k <= ktop; ++k) {
      xhc[k] = (xhc[k - 1] * xzu[k - 1] - T(.5) * up_massdetro[k - 1] * xhc[k - 1] +
                up_massentro[k - 1] * xhe[k - 1]) /
               (xzu[k - 1] - T(.5) * up_massdetro[k - 1] + up_massentro[k - 1]);
      xdby[k] = xhc[k] - xhes_cup[k];
    }
    for (int k = ktop + 1; k <= ktf; ++k) {
      xhc[k] = xhes_cup[k];
      xdby[k] = 0;
    }
  }
  const T xaa0 = cup_up_aa0(xz, xzu, xdby, gamma_cup, xt_cup, kbcon, ktop, ierr, ktf);
  if (ierr == 0) {
    xaa0_ens = xaa0;
    for (int k = 1; k <= ktop; ++k)
      for (int n = 1; n <= maxens3; ++n) pr_ens[n] = pr_ens[n] + pwo[k] + edto * pwdo[k];
    if (pr_ens[7] < T(1.e-6)) {
      ierr = 18;
      for (int n = 1; n <= maxens3; ++n) pr_ens[n] = T(0.);
    }
    for (int n = 1; n <= maxens3; ++n)
      if (pr_ens[n] < T(1.e-5)) pr_ens[n] = T(0.);
  }
  T mconv = T(0);
  if (ierr == 0)
    for (int k = 1; k <= ktop; ++k) {
      const T dq = (qo_cup[k + 1] - qo_cup[k]);
      mconv = mconv + omeg[k] * dq / K::g;
    }
  T xf_dicycle = T(0);
#if !defined(__CUDA_ARCH__)
  if (in.dbg) {
    auto sum = [&](const T* a) { T r = 0; for (int k = 1; k <= kte; ++k) r = r + a[k]; return double(r); };
    std::printf("C ierr k22 kbcon ktop jmin kdet %5d%5d%5d%5d%5d%5d\n", ierr, k22, kbcon, ktop, jmin, kdet);
    std::printf("C aa0 aa1 xaa0 edto mconv %24.16E%24.16E%24.16E%24.16E%24.16E\n", double(aa0), double(aa1),
                double(xaa0), double(edto), double(mconv));
    std::printf("C pwavo pwevo tau aa1bl %24.16E%24.16E%24.16E%24.16E\n", double(pwavo), double(pwevo),
                double(tau_ecmwf), double(aa1_bl));
    std::printf("C sumzuo sumzdo sumdellat %24.16E%24.16E%24.16E\n", sum(zuo), sum(zdo), sum(dellat));
    std::printf("C sumpwo sumpwdo sumdellaq %24.16E%24.16E%24.16E\n", sum(pwo), sum(pwdo), sum(dellaq));
    std::printf("C sumhco sumhcdo sumqco %24.16E%24.16E%24.16E\n", sum(hco), sum(hcdo), sum(qco));
    std::printf("C sumdbyo sumqrco sumcd %24.16E%24.16E%24.16E\n", sum(dbyo), sum(qrco), sum(cd));
    std::printf("C sumentro sumdetro sumentr2d %24.16E%24.16E%24.16E\n", sum(up_massentro), sum(up_massdetro),
                sum(entr_rate_2d));
    std::printf("C hkbo hkb pr_ens7 %24.16E%24.16E%24.16E\n", double(hkbo), double(hkb), double(pr_ens[7]));
    for (int k = 1; k <= ktf; ++k)
      std::printf("CL%3d%24.16E%24.16E%24.16E%24.16E%24.16E%24.16E\n", k, double(heso[k]), double(heso_cup[k]),
                  double(hco[k]), double(dbyo[k]), double(zuo[k]), double(qeso[k]));
  }
#endif
  cup_forcing_ens_3d(closure_n, xland1, aa0, aa1, xaa0_ens, mbdt, dtime, ierr, xf_ens, axx, &mconv, po_cup, ktop, omeg,
                     zdo, zuo, pr_ens, edto, kbcon, in.ichoice, kte, tau_ecmwf, aa1_bl, xf_dicycle);
  // dellat_ens = dellat (ierr = 0) else 0: cup_output_ens_3d reads them only for ierr = 0
  T xmb = T(0);
  cup_output_ens_3d(xf_ens, ierr, dellat, dellaq, dellaqc, outt, outq, outqc, pre, pwo, xmb, ktop, edto, pwdo, po_cup,
                    pr_ens, sig, closure_n, kte, xf_dicycle);
  if (ierr == 0 && pre > T(0.)) {
    pre = fmx(pre, T(0.));
    xmb_out = xmb;
    for (int k = 1; k <= ktop; ++k) {
      outu[k] = dellu[k] * xmb;
      outv[k] = dellv[k] * xmb;
    }
  } else if (ierr != 0 || pre == T(0.)) {
    ktop = 0;
    for (int k = 1; k <= kte; ++k) {
      outt[k] = 0; outq[k] = 0; outqc[k] = 0; outu[k] = 0; outv[k] = 0;
    }
  }
  if (ierr == 0) {
    T dts = T(0.), fpi = T(0.);
    for (int k = 1; k <= ktop; ++k) {
      const T dp = (po_cup[k] - po_cup[k + 1]) * T(100.);
      dts = dts - (outu[k] * us[k] + outv[k] * vs[k]) * dp / K::g;
      fpi = fpi + std::sqrt(outu[k] * outu[k] + outv[k] * outv[k]) * dp;
    }
    if (fpi > T(0.))
      for (int k = 1; k <= ktop; ++k) {
        const T fp = std::sqrt((outu[k] * outu[k] + outv[k] * outv[k])) / fpi;
        outt[k] = outt[k] + fp * dts * K::g / K::cp;
      }
  }
  return GfDeepOut{ierr, kbcon, ktop, k22, jmin};
}

}  // namespace cppwrf::phys::cu_gf
