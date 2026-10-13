#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "phys/phys_constants.hpp"

namespace cppwrf::phys::ysu {

static constexpr int kMaxLevels = 128;

template <typename T>
CPPWRF_HD inline void tridin_ysu(
    const T* cl,
    const T* cm,
    const T* cu,
    const T* r2,
    T* au,
    T* f2,
    int nz) {
  T aul[kMaxLevels];
  for (int k = 0; k < nz; ++k) {
    aul[k] = static_cast<T>(0.0);
  }

  T fk = static_cast<T>(1.0) / cm[0];
  aul[0] = fk * cu[0];
  f2[0] = fk * r2[0];

  for (int k = 1; k < nz - 1; ++k) {
    fk = static_cast<T>(1.0) / (cm[k] - cl[k - 1] * aul[k - 1]);
    aul[k] = fk * cu[k];
    f2[k] = fk * (r2[k] - cl[k - 1] * f2[k - 1]);
  }

  fk = static_cast<T>(1.0) / (cm[nz - 1] - cl[nz - 2] * aul[nz - 2]);
  f2[nz - 1] = fk * (r2[nz - 1] - cl[nz - 2] * f2[nz - 2]);

  for (int k = nz - 2; k >= 0; --k) {
    f2[k] = f2[k] - aul[k] * f2[k + 1];
  }
}

template <typename T>
CPPWRF_HD inline void tridi2n(
    const T* cl,
    const T* cm,
    const T* cm1,
    const T* cu,
    const T* r1,
    const T* r2,
    T* au,
    T* f1,
    T* f2,
    int nz) {
  T fk = static_cast<T>(1.0) / cm[0];
  au[0] = fk * cu[0];
  f1[0] = fk * r1[0];

  fk = static_cast<T>(1.0) / cm1[0];
  f2[0] = fk * r2[0];

  for (int k = 1; k < nz - 1; ++k) {
    fk = static_cast<T>(1.0) / (cm[k] - cl[k - 1] * au[k - 1]);
    au[k] = fk * cu[k];
    f1[k] = fk * (r1[k] - cl[k - 1] * f1[k - 1]);
  }

  for (int k = 1; k < nz - 1; ++k) {
    fk = static_cast<T>(1.0) / (cm1[k] - cl[k - 1] * au[k - 1]);
    f2[k] = fk * (r2[k] - cl[k - 1] * f2[k - 1]);
  }

  fk = static_cast<T>(1.0) / (cm[nz - 1] - cl[nz - 2] * au[nz - 2]);
  f1[nz - 1] = fk * (r1[nz - 1] - cl[nz - 2] * f1[nz - 2]);

  fk = static_cast<T>(1.0) / (cm1[nz - 1] - cl[nz - 2] * au[nz - 2]);
  f2[nz - 1] = fk * (r2[nz - 1] - cl[nz - 2] * f2[nz - 2]);

  for (int k = nz - 2; k >= 0; --k) {
    f1[k] = f1[k] - au[k] * f1[k + 1];
    f2[k] = f2[k] - au[k] * f2[k + 1];
  }
}

template <typename T>
CPPWRF_HD inline void get_pblh(
    int nz,
    T& zi,
    const T* thetav1d,
    const T* qke1d,
    const T* zw1d,
    const T* dz1d,
    T landsea) {
  using std::max;
  using std::min;
  using std::tanh;

  constexpr T sbl_lim = static_cast<T>(200.0);
  constexpr T sbl_damp = static_cast<T>(400.0);

  int k = 1;
  int kthv = 0;
  int ktke = 0;
  T maxqke = static_cast<T>(0.0);
  T minthv = static_cast<T>(9.e9);

  while (k < nz && zw1d[k] <= static_cast<T>(500.0)) {
    T qtke = max(qke1d[k], static_cast<T>(0.0));
    if (maxqke < qtke) {
      maxqke = qtke;
      ktke = k;
    }
    if (minthv > thetav1d[k]) {
      minthv = thetav1d[k];
      kthv = k;
    }
    k++;
  }

  T tkeeps = maxqke / static_cast<T>(40.0);
  tkeeps = max(tkeeps, static_cast<T>(0.025));
  tkeeps = min(tkeeps, static_cast<T>(0.25));

  T delt_thv = ((landsea - static_cast<T>(1.5)) >= static_cast<T>(0.0)) ? static_cast<T>(0.75) : static_cast<T>(1.5);

  zi = static_cast<T>(0.0);
  k = kthv + 1;
  while (zi == static_cast<T>(0.0) && k < nz) {
    if (thetav1d[k] >= (minthv + delt_thv)) {
      zi = zw1d[k] - dz1d[k - 1] *
           min((thetav1d[k] - (minthv + delt_thv)) / max(thetav1d[k] - thetav1d[k - 1], static_cast<T>(1e-6)), static_cast<T>(1.0));
    }
    k++;
    if (k == nz - 1 && zi == static_cast<T>(0.0)) {
      zi = zw1d[1];
    }
  }
  if (zi == static_cast<T>(0.0)) zi = zw1d[1];

  T pblh_tke = static_cast<T>(0.0);
  k = ktke + 1;
  while (pblh_tke == static_cast<T>(0.0) && k < nz) {
    T qtke = max(qke1d[k] / static_cast<T>(2.0), static_cast<T>(0.0));
    T qtkem1 = max(qke1d[k - 1] / static_cast<T>(2.0), static_cast<T>(0.0));
    if (qtke <= tkeeps) {
      pblh_tke = zw1d[k] - dz1d[k - 1] *
                 min((tkeeps - qtke) / max(qtkem1 - qtke, static_cast<T>(1e-6)), static_cast<T>(1.0));
      pblh_tke = max(pblh_tke, zw1d[1]);
    }
    k++;
    if (k == nz - 1 && pblh_tke == static_cast<T>(0.0)) {
      pblh_tke = zw1d[1];
    }
  }
  if (pblh_tke == static_cast<T>(0.0)) pblh_tke = zw1d[1];

  T wt = static_cast<T>(0.5) * tanh((zi - sbl_lim) / sbl_damp) + static_cast<T>(0.5);
  zi = pblh_tke * (static_cast<T>(1.0) - wt) + zi * wt;
}

struct YSUKernel {
  // 3D atmosphere inputs
  View3<const real> ux;
  View3<const real> vx;
  View3<const real> tx;
  View3<const real> qvx;
  View3<const real> qcx;
  View3<const real> qix;
  View3<const real> p2d;
  View3<const real> p2di;
  View3<const real> pi2d;
  View3<const real> dz8w2d;
  View3<const real> rthraten;

  // 2D surface state
  View2<const real> psfcpa;
  View2<const real> znt;
  View2<const real> ust;
  View2<const real> psim;
  View2<const real> psih;
  View2<const real> xland;
  View2<const real> hfx;
  View2<const real> qfx;
  View2<const real> wspd;
  View2<const real> br;
  View2<real> u10;
  View2<real> v10;
  View2<const real> uox;
  View2<const real> vox;
  View2<const real> ctopo;
  View2<const real> ctopo2;

  // 3D tendencies and diagnostics output
  View3<real> utnp;
  View3<real> vtnp;
  View3<real> ttnp;
  View3<real> qvtnp;
  View3<real> qctnp;
  View3<real> qitnp;
  View3<real> exch_hx;
  View3<real> exch_mx;

  // 2D diagnostics output
  View2<real> hpbl;
  View2<int> kpbl1d;
  View2<real> wstar;
  View2<real> delta;
  View2<real> dusfc;
  View2<real> dvsfc;
  View2<real> dtsfc;
  View2<real> dqsfc;

  // Parameters
  real dt;
  int nz;
  bool f_qc;
  bool f_qi;
  bool ysu_topdown_pblmix;
  bool have_ctopo;
  bool have_surface_diag;

  real cp;
  real g;
  real rovcp;
  real rd;
  real rovg;
  real ep1;
  real ep2;
  real karman;
  real xlv;
  real rv;

  CPPWRF_HD void operator()(int i, int j) const {
    using std::abs;
    using std::exp;
    using std::max;
    using std::min;
    using std::pow;
    using std::sqrt;

    // Physical constants from WRF bl_ysu.F90 lines 141-155
    constexpr real xkzminm = real(0.1);
    constexpr real xkzminh = real(0.01);
    constexpr real xkzmax = real(1000.0);
    constexpr real rimin = real(-100.0);
    constexpr real rlam = real(30.0);
    constexpr real prmin = real(0.25);
    constexpr real prmax = real(4.0);
    constexpr real brcr_ub = real(0.0);
    constexpr real brcr_sb = real(0.25);
    constexpr real cori = real(1.e-4);
    constexpr real afac = real(6.8);
    constexpr real bfac = real(6.8);
    constexpr real pfac = real(2.0);
    constexpr real pfac_q = real(2.0);
    constexpr real phifac = real(8.0);
    constexpr real sfcfrac = real(0.1);
    constexpr real d1 = real(0.02);
    constexpr real d2 = real(0.05);
    constexpr real d3 = real(0.001);
    constexpr real h1 = real(0.33333335);
    constexpr real h2 = real(0.6666667);
    constexpr real zfmin = real(1.e-8);
    constexpr real aphi5 = real(5.0);
    constexpr real aphi16 = real(16.0);
    constexpr real tmin = real(1.e-2);
    constexpr real gamcrt = real(3.0);
    constexpr real gamcrq = real(2.e-3);
    constexpr int imvdif = 1;
    constexpr real rcl = real(1.0);

    const real cont = cp / g;
    const real conq = xlv / g;
    const real conw = real(1.0) / g;
    const real conwrc = conw * sqrt(rcl);
    const real conpr = bfac * karman * sfcfrac;
    const real dtstep = dt;
    const real dt2 = real(2.0) * dtstep;
    const real rdt = real(1.0) / dt2;

    // Local column scratch arrays (kMaxLevels = 128)
    real zq[kMaxLevels + 1];
    real za[kMaxLevels];
    real dzq[kMaxLevels];
    real del[kMaxLevels];
    real dza[kMaxLevels];
    real rhox2[kMaxLevels];

    real thx[kMaxLevels];
    real thvx[kMaxLevels];
    real thlix[kMaxLevels];
    real qcxl[kMaxLevels];
    real qixl[kMaxLevels];

    real xkzm[kMaxLevels];
    real xkzh[kMaxLevels];
    real xkzq[kMaxLevels];
    real xkzml[kMaxLevels];
    real xkzhl[kMaxLevels];
    real xkzom[kMaxLevels];
    real xkzoh[kMaxLevels];

    real wscalek[kMaxLevels];
    real wscalek2[kMaxLevels];
    real zfac[kMaxLevels];
    real zfacent[kMaxLevels];
    real entfac[kMaxLevels];

    real au[kMaxLevels];
    real al[kMaxLevels];
    real ad[kMaxLevels];
    real f1[kMaxLevels];
    real f2[kMaxLevels];
    real r1[kMaxLevels];
    real r2[kMaxLevels];
    real cu[kMaxLevels];
    real adv[kMaxLevels];
    real ad1[kMaxLevels];
    real hgamt2[kMaxLevels];

    real tke_ysu[kMaxLevels];

    // Hydrometeor flags (WRF lines 390-416)
    for (int k = 0; k < nz; ++k) {
      qcxl[k] = f_qc ? qcx(i, j, k) : real(0.0);
      qixl[k] = f_qi ? qix(i, j, k) : real(0.0);
      thx[k] = tx(i, j, k) / pi2d(i, j, k);
      thlix[k] = (tx(i, j, k) - xlv * qcxl[k] / cp - real(2.834e6) * qixl[k] / cp) / pi2d(i, j, k);
      real tvcon = real(1.0) + ep1 * qvx(i, j, k);
      thvx[k] = thx[k] * tvcon;
    }

    real uoxl = uox.p ? uox(i, j) : real(0.0);
    real voxl = vox.p ? vox(i, j) : real(0.0);

    real tvcon0 = real(1.0) + ep1 * qvx(i, j, 0);
    real rhox = psfcpa(i, j) / (rd * tx(i, j, 0) * tvcon0);
    real govrth = g / thx[0];

    // Height of levels (WRF lines 507-535)
    zq[0] = real(0.0);
    for (int k = 0; k < nz; ++k) {
      zq[k + 1] = dz8w2d(i, j, k) + zq[k];
      real tvcon = real(1.0) + ep1 * qvx(i, j, k);
      rhox2[k] = p2d(i, j, k) / (rd * tx(i, j, k) * tvcon);
      za[k] = real(0.5) * (zq[k] + zq[k + 1]);
      dzq[k] = zq[k + 1] - zq[k];
      del[k] = p2di(i, j, k) - p2di(i, j, k + 1);
    }

    dza[0] = za[0];
    for (int k = 1; k < nz; ++k) {
      dza[k] = za[k] - za[k - 1];
    }

    // Initialize exchange coefficients and scratch arrays (WRF lines 538-592)
    for (int k = 0; k < nz; ++k) {
      exch_hx(i, j, k) = real(0.0);
      exch_mx(i, j, k) = real(0.0);
      xkzh[k] = real(0.0);
      xkzhl[k] = real(0.0);
      xkzm[k] = real(0.0);
      xkzml[k] = real(0.0);
      xkzq[k] = real(0.0);
      wscalek[k] = real(0.0);
      wscalek2[k] = real(0.0);
      zfac[k] = real(0.0);
      xkzom[k] = xkzminm;
      xkzoh[k] = xkzminh;
    }

    real wspd1 = sqrt((ux(i, j, 0) - uoxl) * (ux(i, j, 0) - uoxl) + (vx(i, j, 0) - voxl) * (vx(i, j, 0) - voxl)) + real(1.e-9);

    real bfxpbl = real(0.0);
    real hfxpbl = real(0.0);
    real qfxpbl = real(0.0);
    real ufxpbl = real(0.0);
    real vfxpbl = real(0.0);
    real hgamu = real(0.0);
    real hgamv = real(0.0);
    delta(i, j) = real(0.0);
    real wstar3_2 = real(0.0);

    if (have_surface_diag) {
      dusfc(i, j) = real(0.0);
      dvsfc(i, j) = real(0.0);
      dtsfc(i, j) = real(0.0);
      dqsfc(i, j) = real(0.0);
    }

    real hgamt = real(0.0);
    real hgamq = real(0.0);
    real wscale = real(0.0);
    int kpbl = 0;  // 0-based in C++ (Fortran kpbl = 1)
    hpbl(i, j) = zq[0];
    real zl1 = za[0];
    real thermal = thvx[0];
    real thermalli = thlix[0];
    bool pblflg = true;
    bool sfcflg = true;
    real sflux = hfx(i, j) / rhox / cp + qfx(i, j) / rhox * ep1 * thx[0];
    if (br(i, j) > real(0.0)) sfcflg = false;

    // First guess of PBL height (WRF lines 618-649)
    bool stable = false;
    real brup = br(i, j);
    real brdn = brup;
    real brcr = brcr_ub;

    for (int k = 1; k < nz; ++k) {
      if (!stable) {
        brdn = brup;
        real spdk2 = max(ux(i, j, k) * ux(i, j, k) + vx(i, j, k) * vx(i, j, k), real(1.0));
        brup = (thvx[k] - thermal) * (g * za[k] / thvx[0]) / spdk2;
        kpbl = k;
        stable = (brup > brcr);
      }
    }

    real brint = real(0.0);
    if (brdn >= brcr) {
      brint = real(0.0);
    } else if (brup <= brcr) {
      brint = real(1.0);
    } else {
      brint = (brcr - brdn) / (brup - brdn);
    }
    hpbl(i, j) = za[kpbl - 1] + brint * (za[kpbl] - za[kpbl - 1]);
    if (hpbl(i, j) < zq[1]) kpbl = 0;
    if (kpbl <= 0) pblflg = false;

    // Similarity functions and scales (WRF lines 650-679)
    real fm_val = psim(i, j);
    real fh_val = psih(i, j);
    real zol1 = max(br(i, j) * fm_val * fm_val / fh_val, rimin);
    if (sfcflg) {
      zol1 = min(zol1, -zfmin);
    } else {
      zol1 = max(zol1, zfmin);
    }
    real hol1 = zol1 * hpbl(i, j) / zl1 * sfcfrac;
    real phim, phih, wstar3;
    if (sfcflg) {
      phim = pow(real(1.0) - aphi16 * hol1, real(-1.0 / 4.0));
      phih = pow(real(1.0) - aphi16 * hol1, real(-1.0 / 2.0));
      real bfx0 = max(sflux, real(0.0));
      wstar3 = govrth * bfx0 * hpbl(i, j);
      wstar(i, j) = pow(wstar3, h1);
    } else {
      phim = real(1.0) + aphi5 * hol1;
      phih = phim;
      wstar(i, j) = real(0.0);
      wstar3 = real(0.0);
    }
    real ust3 = pow(ust(i, j), real(3.0));
    wscale = pow(ust3 + phifac * karman * wstar3 * real(0.5), h1);
    wscale = min(wscale, ust(i, j) * aphi16);
    wscale = max(wscale, ust(i, j) / aphi5);

    // Thermal estimation under unstable conditions (WRF lines 683-700)
    if (sfcflg && sflux > real(0.0)) {
      real gamfac = bfac / rhox / wscale;
      hgamt = min(gamfac * hfx(i, j) / cp, gamcrt);
      hgamq = min(gamfac * qfx(i, j), gamcrq);
      real vpert = (hgamt + ep1 * thx[0] * hgamq) / bfac * afac;
      thermal = thermal + max(vpert, real(0.0)) * min(za[0] / (sfcfrac * hpbl(i, j)), real(1.0));
      thermalli = thermalli + max(vpert, real(0.0)) * min(za[0] / (sfcfrac * hpbl(i, j)), real(1.0));
      hgamt = max(hgamt, real(0.0));
      hgamq = max(hgamq, real(0.0));
      brint = real(-15.9) * ust(i, j) * ust(i, j) / wspd(i, j) * wstar3 / pow(wscale, real(4.0));
      hgamu = brint * ux(i, j, 0);
      hgamv = brint * vx(i, j, 0);
    } else {
      pblflg = false;
    }

    // Enhance PBL height by considering thermal (WRF lines 703-729)
    if (pblflg) {
      kpbl = 0;
      hpbl(i, j) = zq[0];
      stable = false;
      brup = br(i, j);
      brcr = brcr_ub;

      for (int k = 1; k < nz; ++k) {
        if (!stable && pblflg) {
          brdn = brup;
          real spdk2 = max(ux(i, j, k) * ux(i, j, k) + vx(i, j, k) * vx(i, j, k), real(1.0));
          brup = (thvx[k] - thermal) * (g * za[k] / thvx[0]) / spdk2;
          kpbl = k;
          stable = (brup > brcr);
        }
      }

      // Enhance PBL by theta-li (WRF lines 732-753)
      if (ysu_topdown_pblmix) {
        int kpblold = kpbl;
        bool definebrup = false;
        for (int k = kpblold; k < nz - 1; ++k) {
          real spdk2 = max(ux(i, j, k) * ux(i, j, k) + vx(i, j, k) * vx(i, j, k), real(1.0));
          real bruptmp = (thlix[k] - thermalli) * (g * za[k] / thlix[0]) / spdk2;
          stable = (bruptmp >= brcr);
          if (definebrup) {
            kpbl = k;
            brup = bruptmp;
            definebrup = false;
          }
          if (!stable) {
            brdn = bruptmp;
            definebrup = true;
            pblflg = true;
          }
        }
      }

      if (pblflg) {
        if (brdn >= brcr) {
          brint = real(0.0);
        } else if (brup <= brcr) {
          brint = real(1.0);
        } else {
          brint = (brcr - brdn) / (brup - brdn);
        }
        hpbl(i, j) = za[kpbl - 1] + brint * (za[kpbl] - za[kpbl - 1]);
        if (hpbl(i, j) < zq[1]) kpbl = 0;
        if (kpbl <= 0) pblflg = false;
      }
    }

    // Stable boundary layer (WRF lines 772-826)
    if (!sfcflg && hpbl(i, j) < zq[1]) {
      brup = br(i, j);
      stable = false;
    } else {
      stable = true;
    }

    real brcr_sbro = real(0.0);
    if (!stable && (xland(i, j) - real(1.5)) >= real(0.0)) {
      real wspd10 = sqrt(u10(i, j) * u10(i, j) + v10(i, j) * v10(i, j));
      real ross = wspd10 / (cori * znt(i, j));
      brcr_sbro = min(real(0.16) * pow(real(1.e-7) * ross, real(-0.18)), real(0.3));
    }

    if (!stable) {
      if ((xland(i, j) - real(1.5)) >= real(0.0)) {
        brcr = brcr_sbro;
      } else {
        brcr = brcr_sb;
      }

      for (int k = 1; k < nz; ++k) {
        if (!stable) {
          brdn = brup;
          real spdk2 = max(ux(i, j, k) * ux(i, j, k) + vx(i, j, k) * vx(i, j, k), real(1.0));
          brup = (thvx[k] - thermal) * (g * za[k] / thvx[0]) / spdk2;
          kpbl = k;
          stable = (brup > brcr);
        }
      }

      if (!sfcflg && hpbl(i, j) < zq[1]) {
        if (brdn >= brcr) {
          brint = real(0.0);
        } else if (brup <= brcr) {
          brint = real(1.0);
        } else {
          brint = (brcr - brdn) / (brup - brdn);
        }
        hpbl(i, j) = za[kpbl - 1] + brint * (za[kpbl] - za[kpbl - 1]);
        if (hpbl(i, j) < zq[1]) kpbl = 0;
        if (kpbl <= 0) pblflg = false;
      }
    }

    // Entrainment parameters (WRF lines 830-925)
    bool cloudflg = false;
    real wm2 = real(0.0);
    real we = real(0.0);
    real dthvx = real(0.0);

    if (pblflg) {
      int k = kpbl - 1;
      real wm3 = wstar3 + real(5.0) * ust3;
      wm2 = pow(wm3, h2);
      bfxpbl = real(-0.15) * thvx[0] / g * wm3 / hpbl(i, j);
      dthvx = max(thvx[k + 1] - thvx[k], tmin);
      we = max(bfxpbl / dthvx, -sqrt(wm2));

      if ((qcxl[k] + qixl[k]) > real(0.01e-3) && ysu_topdown_pblmix) {
        if (kpbl >= 1) {
          cloudflg = true;
          real templ = thlix[k] * pow(p2di(i, j, k + 1) / real(100000.0), rovcp);
          real rvls = real(100.0) * real(6.112) * exp(real(17.67) * (templ - real(273.16)) / (templ - real(29.65))) * (ep2 / p2di(i, j, k + 1));
          real temps = templ + ((qvx(i, j, k) + qcxl[k]) - rvls) / (cp / xlv + ep2 * xlv * rvls / (rd * templ * templ));
          rvls = real(100.0) * real(6.112) * exp(real(17.67) * (temps - real(273.15)) / (temps - real(29.65))) * (ep2 / p2di(i, j, k + 1));
          real rcldb = max((qvx(i, j, k) + qcxl[k]) - rvls, real(0.0));

          real dthvx_cloud = (thlix[k + 2] + thx[k + 2] * ep1 * (qvx(i, j, k + 2) + qcxl[k + 2])) -
                             (thlix[k] + thx[k] * ep1 * (qvx(i, j, k) + qcxl[k]));
          dthvx_cloud = max(dthvx_cloud, real(0.1));
          real tmp1 = xlv / cp * rcldb / (pi2d(i, j, k) * dthvx_cloud);
          real ent_eff = real(0.2) * real(8.0) * tmp1 + real(0.2);

          real radsum = real(0.0);
          for (int kk = 0; kk <= kpbl - 1; ++kk) {
            real radflux = rthraten(i, j, kk) * pi2d(i, j, kk);
            radflux = radflux * cp / g * (p2di(i, j, kk) - p2di(i, j, kk + 1));
            if (radflux < real(0.0)) radsum = abs(radflux) + radsum;
          }
          radsum = max(radsum, real(0.0));

          real bfx0 = max(sflux, real(0.0));
          wm3 = (govrth * bfx0 * hpbl(i, j)) + real(5.0) * ust3;
          wm2 = pow(wm3, h2);
          bfxpbl = real(-0.15) * thvx[0] / g * wm3 / hpbl(i, j);
          dthvx = max(thvx[k + 1] - thvx[k], tmin);
          we = max(bfxpbl / dthvx, -sqrt(wm2));

          bfx0 = max(radsum / rhox2[k] / cp, real(0.0));
          wm3 = (g / thvx[k] * bfx0 * hpbl(i, j));
          wm2 = wm2 + pow(wm3, h2);
          bfxpbl = -ent_eff * bfx0;
          dthvx = max(thvx[k + 1] - thvx[k], real(0.1));
          we = we + max(bfxpbl / dthvx, -sqrt(pow(wm3, h2)));

          wstar3_2 = (g / thvx[k] * bfx0 * hpbl(i, j));
          wscale = pow(ust3 + phifac * karman * (wstar3 + wstar3_2) * real(0.5), h1);
          wscale = min(wscale, ust(i, j) * aphi16);
          wscale = max(wscale, ust(i, j) / aphi5);
          real gamfac = bfac / rhox / wscale;
          hgamt = min(gamfac * hfx(i, j) / cp, gamcrt);
          hgamq = min(gamfac * qfx(i, j), gamcrq);
          gamfac = bfac / rhox2[k] / wscale;
          hgamt2[k] = min(gamfac * radsum / cp, gamcrt);
          hgamt = max(hgamt, real(0.0)) + max(hgamt2[k], real(0.0));
          brint = real(-15.9) * ust(i, j) * ust(i, j) / wspd(i, j) * (wstar3 + wstar3_2) / pow(wscale, real(4.0));
          hgamu = brint * ux(i, j, 0);
          hgamv = brint * vx(i, j, 0);
        }
      }

      real prpbl = real(1.0);
      real dthx = max(thx[k + 1] - thx[k], tmin);
      real dqx = min(qvx(i, j, k + 1) - qvx(i, j, k), real(0.0));
      hfxpbl = we * dthx;
      qfxpbl = we * dqx;

      real dux = ux(i, j, k + 1) - ux(i, j, k);
      real dvx = vx(i, j, k + 1) - vx(i, j, k);
      if (dux > tmin) {
        ufxpbl = max(prpbl * we * dux, -ust(i, j) * ust(i, j));
      } else if (dux < -tmin) {
        ufxpbl = min(prpbl * we * dux, ust(i, j) * ust(i, j));
      } else {
        ufxpbl = real(0.0);
      }

      if (dvx > tmin) {
        vfxpbl = max(prpbl * we * dvx, -ust(i, j) * ust(i, j));
      } else if (dvx < -tmin) {
        vfxpbl = min(prpbl * we * dvx, ust(i, j) * ust(i, j));
      } else {
        vfxpbl = real(0.0);
      }

      real delb = govrth * d3 * hpbl(i, j);
      delta(i, j) = min(d1 * hpbl(i, j) + d2 * wm2 / delb, real(100.0));
    }

    // Entrainment factor entfac (WRF lines 927-935)
    for (int k = 0; k < nz; ++k) {
      if (pblflg && k >= kpbl) {
        real ratio = (zq[k + 1] - hpbl(i, j)) / delta(i, j);
        entfac[k] = ratio * ratio;
      } else {
        entfac[k] = real(1.e30);
      }
    }

    // Diffusion coefficients below PBL (WRF lines 939-979)
    for (int k = 0; k < nz; ++k) {
      if (k < kpbl) {
        zfac[k] = min(max(real(1.0) - (zq[k + 1] - zl1) / (hpbl(i, j) - zl1), zfmin), real(1.0));
        zfacent[k] = pow(real(1.0) - zfac[k], real(3.0));
        wscalek[k] = pow(ust3 + phifac * karman * wstar3 * (real(1.0) - zfac[k]), h1);
        wscalek2[k] = pow(phifac * karman * wstar3_2 * zfac[k], h1);

        real prfac, prfac2, prnumfac;
        if (sfcflg) {
          prfac = conpr;
          prfac2 = real(15.9) * (wstar3 + wstar3_2) / ust3 / (real(1.0) + real(4.0) * karman * (wstar3 + wstar3_2) / ust3);
          real term = max(zq[k + 1] - sfcfrac * hpbl(i, j), real(0.0));
          prnumfac = real(-3.0) * term * term / (hpbl(i, j) * hpbl(i, j));
        } else {
          prfac = real(0.0);
          prfac2 = real(0.0);
          prnumfac = real(0.0);
          real phim8z = real(1.0) + aphi5 * zol1 * zq[k + 1] / zl1;
          wscalek[k] = ust(i, j) / phim8z;
          wscalek[k] = max(wscalek[k], real(0.001));
        }

        real prnum0 = (phih / phim + prfac);
        prnum0 = max(min(prnum0, prmax), prmin);
        xkzm[k] = wscalek[k] * karman * zq[k + 1] * pow(zfac[k], pfac) +
                  wscalek2[k] * karman * (hpbl(i, j) - zq[k + 1]) * pow(real(1.0) - zfac[k], pfac);

        if (k == kpbl - 1 && cloudflg && we < real(0.0)) {
          xkzm[k] = real(0.0);
        }

        real prnum = real(1.0) + (prnum0 - real(1.0)) * exp(prnumfac);
        xkzq[k] = xkzm[k] / prnum * pow(zfac[k], pfac_q - pfac);
        prnum0 = prnum0 / (real(1.0) + prfac2 * karman * sfcfrac);
        prnum = real(1.0) + (prnum0 - real(1.0)) * exp(prnumfac);
        xkzh[k] = xkzm[k] / prnum;
        xkzm[k] = xkzm[k] + xkzom[k];
        xkzh[k] = xkzh[k] + xkzoh[k];
        xkzq[k] = xkzq[k] + xkzoh[k];
        xkzm[k] = min(xkzm[k], xkzmax);
        xkzh[k] = min(xkzh[k], xkzmax);
        xkzq[k] = min(xkzq[k], xkzmax);
      }
    }

    // Diffusion coefficients in free atmosphere over PBL (WRF lines 983-1029)
    for (int k = 0; k < nz - 1; ++k) {
      if (k >= kpbl) {
        real du = ux(i, j, k + 1) - ux(i, j, k);
        real dv = vx(i, j, k + 1) - vx(i, j, k);
        real ss = (du * du + dv * dv) / (dza[k + 1] * dza[k + 1]) + real(1.e-9);
        real govrthv = g / (real(0.5) * (thvx[k + 1] + thvx[k]));
        real ri = govrthv * (thvx[k + 1] - thvx[k]) / (ss * dza[k + 1]);

        if (imvdif == 1) {
          if ((qcxl[k] + qixl[k]) > real(0.01e-3) && (qcxl[k + 1] + qixl[k + 1]) > real(0.01e-3)) {
            real qmean = real(0.5) * (qvx(i, j, k) + qvx(i, j, k + 1));
            real tmean = real(0.5) * (tx(i, j, k) + tx(i, j, k + 1));
            real alph = xlv * qmean / rd / tmean;
            real chi = xlv * xlv * qmean / cp / rv / tmean / tmean;
            ri = (real(1.0) + alph) * (ri - g * g / ss / tmean / cp * ((chi - alph) / (real(1.0) + chi)));
          }
        }

        real zk = karman * zq[k + 1];
        real rlamdz = min(max(real(0.1) * dza[k + 1], rlam), real(300.0));
        rlamdz = min(dza[k + 1], rlamdz);
        real rl2 = (zk * rlamdz / (rlamdz + zk)) * (zk * rlamdz / (rlamdz + zk));
        real dk = rl2 * sqrt(ss);

        if (ri < real(0.0)) {
          ri = max(ri, rimin);
          real sri = sqrt(-ri);
          xkzm[k] = dk * (real(1.0) + real(8.0) * (-ri) / (real(1.0) + real(1.746) * sri));
          xkzh[k] = dk * (real(1.0) + real(8.0) * (-ri) / (real(1.0) + real(1.286) * sri));
        } else {
          xkzh[k] = dk / ((real(1.0) + real(5.0) * ri) * (real(1.0) + real(5.0) * ri));
          real prnum = real(1.0) + real(2.1) * ri;
          prnum = min(prnum, prmax);
          xkzm[k] = xkzh[k] * prnum;
        }

        xkzm[k] = xkzm[k] + xkzom[k];
        xkzh[k] = xkzh[k] + xkzoh[k];
        xkzm[k] = min(xkzm[k], xkzmax);
        xkzh[k] = min(xkzh[k], xkzmax);
        xkzml[k] = xkzm[k];
        xkzhl[k] = xkzh[k];
      }
    }

    // =========================================================================
    // Heat tridiagonal system (WRF lines 1033-1112)
    // =========================================================================
    for (int k = 0; k < nz; ++k) {
      au[k] = real(0.0);
      al[k] = real(0.0);
      ad[k] = real(0.0);
      f1[k] = real(0.0);
    }

    ad[0] = real(1.0);
    f1[0] = thx[0] - real(300.0) + hfx(i, j) / cont / del[0] * dt2;

    for (int k = 0; k < nz - 1; ++k) {
      real dtodsd = dt2 / del[k];
      real dtodsu = dt2 / del[k + 1];
      real dsig = p2d(i, j, k) - p2d(i, j, k + 1);
      real rdz = real(1.0) / dza[k + 1];
      real tem1 = dsig * xkzh[k] * rdz;

      if (pblflg && k < kpbl) {
        real dsdzt = tem1 * (-hgamt / hpbl(i, j) - hfxpbl * zfacent[k] / xkzh[k]);
        f1[k] = f1[k] + dtodsd * dsdzt;
        f1[k + 1] = thx[k + 1] - real(300.0) - dtodsu * dsdzt;
      } else if (pblflg && k >= kpbl && entfac[k] < real(4.6)) {
        xkzh[k] = -we * dza[kpbl] * exp(-entfac[k]);
        xkzh[k] = sqrt(xkzh[k] * xkzhl[k]);
        xkzh[k] = max(xkzh[k], xkzoh[k]);
        xkzh[k] = min(xkzh[k], xkzmax);
        f1[k + 1] = thx[k + 1] - real(300.0);
      } else {
        f1[k + 1] = thx[k + 1] - real(300.0);
      }

      tem1 = dsig * xkzh[k] * rdz;
      real dsdz2 = tem1 * rdz;
      au[k] = -dtodsd * dsdz2;
      al[k] = -dtodsu * dsdz2;
      ad[k] = ad[k] - au[k];
      ad[k + 1] = real(1.0) - al[k];
      exch_hx(i, j, k + 1) = xkzh[k];
    }

    for (int k = 0; k < nz; ++k) {
      cu[k] = au[k];
      r1[k] = f1[k];
    }

    tridin_ysu(al, ad, cu, r1, au, f1, nz);

    for (int k = nz - 1; k >= 0; --k) {
      real ttend = (f1[k] - thx[k] + real(300.0)) * rdt;
      ttnp(i, j, k) = ttend;
      if (have_surface_diag) {
        dtsfc(i, j) = dtsfc(i, j) + ttend * cont * del[k];
      }
    }

    // =========================================================================
    // Moisture tridiagonal systems (WRF lines 1115-1239)
    // =========================================================================
    for (int k = 0; k < nz - 1; ++k) {
      if (k >= kpbl) xkzq[k] = xkzh[k];
    }

    for (int k = 0; k < nz; ++k) {
      au[k] = real(0.0);
      al[k] = real(0.0);
      ad[k] = real(0.0);
      f1[k] = real(0.0);
      r1[k] = real(0.0);
    }

    ad[0] = real(1.0);
    f1[0] = qvx(i, j, 0) + qfx(i, j) * g / del[0] * dt2;

    for (int k = 0; k < nz - 1; ++k) {
      real dtodsd = dt2 / del[k];
      real dtodsu = dt2 / del[k + 1];
      real dsig = p2d(i, j, k) - p2d(i, j, k + 1);
      real rdz = real(1.0) / dza[k + 1];
      real tem1 = dsig * xkzq[k] * rdz;

      if (pblflg && k < kpbl) {
        real dsdzq = tem1 * (-qfxpbl * zfacent[k] / xkzq[k]);
        f1[k] = f1[k] + dtodsd * dsdzq;
        f1[k + 1] = qvx(i, j, k + 1) - dtodsu * dsdzq;
      } else if (pblflg && k >= kpbl && entfac[k] < real(4.6)) {
        xkzq[k] = -we * dza[kpbl] * exp(-entfac[k]);
        xkzq[k] = sqrt(xkzq[k] * xkzhl[k]);
        xkzq[k] = max(xkzq[k], xkzoh[k]);
        xkzq[k] = min(xkzq[k], xkzmax);
        f1[k + 1] = qvx(i, j, k + 1);
      } else {
        f1[k + 1] = qvx(i, j, k + 1);
      }

      tem1 = dsig * xkzq[k] * rdz;
      real dsdz2 = tem1 * rdz;
      au[k] = -dtodsd * dsdz2;
      al[k] = -dtodsu * dsdz2;
      ad[k] = ad[k] - au[k];
      ad[k + 1] = real(1.0) - al[k];
    }

    for (int k = 0; k < nz; ++k) {
      adv[k] = ad[k];
      cu[k] = au[k];
      r1[k] = f1[k];
    }

    tridin_ysu(al, adv, cu, r1, au, f1, nz);

    for (int k = nz - 1; k >= 0; --k) {
      real qtend = (f1[k] - qvx(i, j, k)) * rdt;
      qvtnp(i, j, k) = qtend;
      if (have_surface_diag) {
        dqsfc(i, j) = dqsfc(i, j) + qtend * conq * del[k];
      }
    }

    // Cloud water (qc)
    if (f_qc) {
      for (int k = 0; k < nz; ++k) {
        f1[k] = qcxl[k];
        r1[k] = f1[k];
      }
      tridin_ysu(al, ad, cu, r1, au, f1, nz);
      for (int k = nz - 1; k >= 0; --k) {
        qctnp(i, j, k) = (f1[k] - qcxl[k]) * rdt;
      }
    }

    // Cloud ice (qi)
    if (f_qi) {
      for (int k = 0; k < nz; ++k) {
        f1[k] = qixl[k];
        r1[k] = f1[k];
      }
      tridin_ysu(al, ad, cu, r1, au, f1, nz);
      for (int k = nz - 1; k >= 0; --k) {
        qitnp(i, j, k) = (f1[k] - qixl[k]) * rdt;
      }
    }

    // =========================================================================
    // Momentum tridiagonal system (WRF lines 1243-1408)
    // =========================================================================
    for (int k = 0; k < nz; ++k) {
      au[k] = real(0.0);
      al[k] = real(0.0);
      ad[k] = real(0.0);
      f1[k] = real(0.0);
      f2[k] = real(0.0);
    }

    for (int k = 0; k < nz - 1; ++k) {
      real du_dz = (ux(i, j, k + 1) - ux(i, j, k)) / dza[k + 1];
      real dv_dz = (vx(i, j, k + 1) - vx(i, j, k)) / dza[k + 1];
      real term_u = -hgamu / hpbl(i, j) + du_dz;
      real term_v = -hgamv / hpbl(i, j) + dv_dz;
      real shear_ysu = xkzm[k] * (term_u * du_dz + term_v * dv_dz);
      real buoy_ysu = xkzh[k] * g * (real(1.0) / thx[k]) * (-hgamt / hpbl(i, j) + (thx[k + 1] - thx[k]) / dza[k + 1]);

      real zk = karman * zq[k + 1];
      real rlamdz;
      if (k >= kpbl) {
        rlamdz = min(max(real(0.1) * dza[k + 1], rlam), real(300.0));
        rlamdz = min(dza[k + 1], rlamdz);
      } else {
        rlamdz = real(150.0);
      }
      real el_ysu = zk * rlamdz / (rlamdz + zk);
      real tke_val = real(16.6) * el_ysu * (shear_ysu - buoy_ysu);
      if (tke_val <= real(0.0)) {
        tke_ysu[k] = real(0.0);
      } else {
        tke_ysu[k] = pow(tke_val, real(0.66));
      }
    }
    tke_ysu[nz - 1] = real(0.0);

    real pblh_ysu = real(0.0);
    get_pblh(nz, pblh_ysu, thvx, tke_ysu, zq, dzq, xland(i, j));

    real vconv = real(0.0);
    if (xland(i, j) < real(1.5)) {
      real fluxc = max(sflux, real(0.0));
      vconv = pow(g / thvx[0] * pblh_ysu * fluxc, real(0.33));
    }
    real vconvfx = vconv;

    real wspd_ratio = (wspd1 / wspd(i, j));
    real fric = ust(i, j) * ust(i, j) / wspd1 * rhox * g / del[0] * dt2 * (wspd_ratio * wspd_ratio);

    if (have_ctopo) {
      real vconvnew = real(0.9) * vconvfx + real(1.5) * (max((pblh_ysu - real(500.0)) / real(1000.0), real(0.0)));
      real vconvlim = min(vconvnew, real(1.0));
      ad[0] = real(1.0) + fric * vconvlim + ctopo(i, j) * fric * (real(1.0) - vconvlim);
    } else {
      ad[0] = real(1.0) + fric;
    }

    f1[0] = ux(i, j, 0) + uoxl * ust(i, j) * ust(i, j) * rhox * g / del[0] * dt2 / wspd1 * (wspd_ratio * wspd_ratio);
    f2[0] = vx(i, j, 0) + voxl * ust(i, j) * ust(i, j) * rhox * g / del[0] * dt2 / wspd1 * (wspd_ratio * wspd_ratio);

    for (int k = 0; k < nz - 1; ++k) {
      real dtodsd = dt2 / del[k];
      real dtodsu = dt2 / del[k + 1];
      real dsig = p2d(i, j, k) - p2d(i, j, k + 1);
      real rdz = real(1.0) / dza[k + 1];
      real tem1 = dsig * xkzm[k] * rdz;

      if (pblflg && k < kpbl) {
        real dsdzu = tem1 * (-hgamu / hpbl(i, j) - ufxpbl * zfacent[k] / xkzm[k]);
        real dsdzv = tem1 * (-hgamv / hpbl(i, j) - vfxpbl * zfacent[k] / xkzm[k]);
        f1[k] = f1[k] + dtodsd * dsdzu;
        f1[k + 1] = ux(i, j, k + 1) - dtodsu * dsdzu;
        f2[k] = f2[k] + dtodsd * dsdzv;
        f2[k + 1] = vx(i, j, k + 1) - dtodsu * dsdzv;
      } else if (pblflg && k >= kpbl && entfac[k] < real(4.6)) {
        xkzm[k] = real(1.0) * xkzh[k];
        xkzm[k] = sqrt(xkzm[k] * xkzml[k]);
        xkzm[k] = max(xkzm[k], xkzom[k]);
        xkzm[k] = min(xkzm[k], xkzmax);
        f1[k + 1] = ux(i, j, k + 1);
        f2[k + 1] = vx(i, j, k + 1);
      } else {
        f1[k + 1] = ux(i, j, k + 1);
        f2[k + 1] = vx(i, j, k + 1);
      }

      tem1 = dsig * xkzm[k] * rdz;
      real dsdz2 = tem1 * rdz;
      au[k] = -dtodsd * dsdz2;
      al[k] = -dtodsu * dsdz2;
      ad[k] = ad[k] - au[k];
      ad[k + 1] = real(1.0) - al[k];
      exch_mx(i, j, k + 1) = xkzm[k];
    }

    for (int k = 0; k < nz; ++k) {
      ad1[k] = ad[k];
      cu[k] = au[k];
      r1[k] = f1[k];
      r2[k] = f2[k];
    }

    tridi2n(al, ad, ad1, cu, r1, r2, au, f1, f2, nz);

    for (int k = nz - 1; k >= 0; --k) {
      real utend = (f1[k] - ux(i, j, k)) * rdt;
      real vtend = (f2[k] - vx(i, j, k)) * rdt;
      utnp(i, j, k) = utend;
      vtnp(i, j, k) = vtend;
      if (have_surface_diag) {
        dusfc(i, j) = dusfc(i, j) + utend * conwrc * del[k];
        dvsfc(i, j) = dvsfc(i, j) + vtend * conwrc * del[k];
      }
    }

    if (have_ctopo) {
      u10(i, j) = ctopo2(i, j) * u10(i, j) + (real(1.0) - ctopo2(i, j)) * ux(i, j, 0);
      v10(i, j) = ctopo2(i, j) * v10(i, j) + (real(1.0) - ctopo2(i, j)) * vx(i, j, 0);
    }

    kpbl1d(i, j) = kpbl + 1;  // Fortran 1-based index
  }
};

void run_ysu(
    View3<const real> ux,
    View3<const real> vx,
    View3<const real> tx,
    View3<const real> qvx,
    View3<const real> qcx,
    View3<const real> qix,
    View3<const real> p2d,
    View3<const real> p2di,
    View3<const real> pi2d,
    View3<const real> dz8w2d,
    View3<const real> rthraten,
    View2<const real> psfcpa,
    View2<const real> znt,
    View2<const real> ust,
    View2<const real> psim,
    View2<const real> psih,
    View2<const real> xland,
    View2<const real> hfx,
    View2<const real> qfx,
    View2<const real> wspd,
    View2<const real> br,
    View2<real> u10,
    View2<real> v10,
    View2<const real> uox,
    View2<const real> vox,
    View2<const real> ctopo,
    View2<const real> ctopo2,
    View3<real> utnp,
    View3<real> vtnp,
    View3<real> ttnp,
    View3<real> qvtnp,
    View3<real> qctnp,
    View3<real> qitnp,
    View3<real> exch_hx,
    View3<real> exch_mx,
    View2<real> hpbl,
    View2<int> kpbl1d,
    View2<real> wstar,
    View2<real> delta,
    View2<real> dusfc,
    View2<real> dvsfc,
    View2<real> dtsfc,
    View2<real> dqsfc,
    real dt,
    const Range2D& range,
    int nz,
    bool f_qc = true,
    bool f_qi = true,
    bool ysu_topdown_pblmix = false,
    bool have_ctopo = false,
    bool have_surface_diag = true,
    real cp = real(cst::cp),
    real g = real(cst::g),
    real rovcp = real(cst::rcp),
    real rd = real(cst::r_d),
    real rovg = real(cst::r_d / cst::g),
    real ep1 = real(cst::ep_1),
    real ep2 = real(cst::ep_2),
    real karman = real(cst::karman),
    real xlv = real(cst::xlv),
    real rv = real(cst::r_v));

}  // namespace cppwrf::phys::ysu
