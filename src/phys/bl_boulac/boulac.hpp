// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Bougeault-Lacarrere PBL scheme (bl_pbl_physics = 8): WRF v4.8.0 phys/module_bl_boulac.F (BOULAC -> pbl_height,
// boulac1D -> dissip_bougeault, length_bougeault, cdtur_bougeault, tke_bougeault (shear, buoy), diff, invert) for one
// column, without the BEP urban terms (flag_bep = .false.: CPPWRF has no multi-layer urban model) and with idiff = 0
// (module_pbl_driver.F sets it). Every statement follows the Fortran in order and parenthesisation (real*8: the literal
// constants are doubles; x**2. = x*x as gfortran folds it, (u**2.+v**2.)**.5 and (...)**(1./3.) are libm pow calls).
// The variances t2_1D / w2_1D of BOULAC (the only loop gfortran vectorises with libmvec's pow) and tstar / tstar_w
// feed only the igamma = 2 / 3 branches; BOULAC fixes igamma = 1, so they are not computed.
// The model's levels k = 0 .. nz-1 are WRF's iz = 1 .. kte; the arrays below use WRF's indices 1 .. kte+1 (kme).
#include <cmath>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "core/libm.hpp"

namespace cppwrf::phys::boulac {

inline constexpr int kMaxLevels = 128;

// module_bl_boulac.F parameters (real*8 under -fdefault-real-8): ceps_b = 1/1.4, ck_b = 0.4, temin = 0.0001
inline constexpr real ceps_b = real(1.0 / 1.4);
inline constexpr real ck_b = real(0.4);
inline constexpr real temin = real(0.0001);

// invert: the tridiagonal system A X = C (a1 lower, a2 principal, a3 upper diagonal; a2 and c are overwritten)
CPPWRF_HD inline void invert(int kte, real* a1, real* a2, const real* a3, real* c, real* x) {
  for (int in = kte - 1; in >= 1; --in) {
    c[in] = c[in] - a3[in] * c[in + 1] / a2[in + 1];
    a2[in] = a2[in] - a3[in] * a1[in + 1] / a2[in + 1];
  }
  for (int in = 2; in <= kte; ++in) c[in] = c[in] - a1[in] * c[in - 1] / a2[in - 1];
  for (int in = 1; in <= kte; ++in) x[in] = c[in] / a2[in];
}

// diff with iz1 = 1, izf = 1 (every call of the scheme): co updated in place, fc the turbulent flux
template <int K>
CPPWRF_HD inline void diff(int kte, real dt, real* co, const real* rho, const real* rhoz, const real* cd,
                           const real* aa, const real* bb, const real* sf, const real* vl, const real* dz, real* fc) {
  real cddz[K + 3], a1[K + 2], a2[K + 2], a3[K + 2], c[K + 2];
  cddz[1] = sf[1] * rhoz[1] * cd[1] / dz[1];
  for (int iz = 2; iz <= kte; ++iz) cddz[iz] = real(2.) * sf[iz] * rhoz[iz] * cd[iz] / (dz[iz] + dz[iz - 1]);
  cddz[kte + 1] = sf[kte + 1] * rhoz[kte + 1] * cd[kte + 1] / dz[kte];
  // do iz = kts, iz1-1: empty
  for (int iz = 1; iz <= kte - 1; ++iz) {
    const real dzv = vl[iz] * dz[iz];
    a1[iz] = -cddz[iz] * dt / dzv / rho[iz];
    a2[iz] = real(1) + dt * (cddz[iz] + cddz[iz + 1]) / dzv / rho[iz] - aa[iz] * dt;
    a3[iz] = -cddz[iz + 1] * dt / dzv / rho[iz];
    c[iz] = co[iz] + bb[iz] * dt;
  }
  a1[kte] = real(0.);
  a2[kte] = real(1);
  a3[kte] = real(0.);
  c[kte] = co[kte];
  invert(kte, a1, a2, a3, c, co);
  fc[1] = real(0.);
  for (int iz = 2; iz <= kte; ++iz) fc[iz] = -(cddz[iz] * (co[iz] - co[iz - 1])) / rho[iz];
}

// pbl_height: the height where thv first exceeds thv(kts) + 0.5 (interpolated between the cell centres)
CPPWRF_HD inline real pbl_height(int kte, const real* dz, const real* z, const real* th, const real* q) {
  real thv[kMaxLevels + 2], zc[kMaxLevels + 2];
  for (int iz = 1; iz <= kte; ++iz) zc[iz] = z[iz] + dz[iz] / real(2.);
  for (int iz = 1; iz <= kte; ++iz) thv[iz] = th[iz] * (real(1.) + real(0.61) * q[iz]);
  real pblh = real(0.);
  const real thsfc = thv[1] + real(0.5);
  for (int iz = 2; iz <= kte; ++iz)
    if (pblh == real(0) && thv[iz] > thsfc)
      pblh = zc[iz - 1] + (thsfc - thv[iz - 1]) / (c_max(real(0.01), thv[iz] - thv[iz - 1])) * (zc[iz] - zc[iz - 1]);
  return pblh;
}

// dissip_bougeault: the length scales up (dlu) and down (dld)
CPPWRF_HD inline void dissip_bougeault(real g, int kte, const real* z, const real* dz, const real* te, real* dlu,
                                       real* dld, const real* th, const real* th0) {
  using std::sqrt;
  for (int iz = 1; iz <= kte; ++iz) {
    real zup = real(0.);
    dlu[iz] = z[kte + 1] - z[iz] - dz[iz] / real(2.);
    real zzz = real(0.);
    real zup_inf = real(0.);
    const real beta = g / th0[iz];  // Buoyancy coefficient
    for (int izz = iz; izz <= kte - 1; ++izz) {
      const real dzt = (dz[izz + 1] + dz[izz]) / real(2.);
      zup = zup - beta * th[iz] * dzt;
      zup = zup + beta * (th[izz + 1] + th[izz]) * dzt / real(2.);
      zzz = zzz + dzt;
      if (te[iz] < zup && te[iz] >= zup_inf) {
        const real bbb = (th[izz + 1] - th[izz]) / dzt;
        real tl;
        if (bbb != real(0)) {
          const real x = beta * (th[izz] - th[iz]);
          tl = (-(beta * (th[izz] - th[iz])) + sqrt(c_max(real(0.), x * x + real(2.) * bbb * beta * (te[iz] - zup_inf)))) /
               bbb / beta;
        } else {
          if (th[izz] != th[iz])
            tl = (te[iz] - zup_inf) / (beta * (th[izz] - th[iz]));
          else
            tl = real(0.);
        }
        dlu[iz] = zzz - dzt + tl;
      }
      zup_inf = zup;
    }
    real zdo = real(0.);
    real zdo_sup = real(0.);
    dld[iz] = z[iz] + dz[iz] / real(2.);
    zzz = real(0.);
    for (int izz = iz; izz >= 2; --izz) {
      const real dzt = (dz[izz - 1] + dz[izz]) / real(2.);
      zdo = zdo + beta * th[iz] * dzt;
      zdo = zdo - beta * (th[izz - 1] + th[izz]) * dzt / real(2.);
      zzz = zzz + dzt;
      if (te[iz] < zdo && te[iz] >= zdo_sup) {
        const real bbb = (th[izz] - th[izz - 1]) / dzt;
        real tl;
        if (bbb != real(0.)) {
          const real x = beta * (th[izz] - th[iz]);
          tl = (beta * (th[izz] - th[iz]) + sqrt(c_max(real(0.), x * x + real(2.) * bbb * beta * (te[iz] - zdo_sup)))) /
               bbb / beta;
        } else {
          if (th[izz] != th[iz])
            tl = (te[iz] - zdo_sup) / (beta * (th[izz] - th[iz]));
          else
            tl = real(0.);
        }
        dld[iz] = zzz - dzt + tl;
      }
      zdo_sup = zdo;
    }
  }
}

// One column of BOULAC (flag_bep = .false., idiff = 0). The View3s are (i, j, k) with k = 0 .. nz-1.
struct BoulacKernel {
  View3<const real> u_phy, v_phy, th_phy, rho, qv, qc, dz8w;  // u_phytmp, v_phytmp, th_phy, rho, qv_curr, qc_curr
  View2<const real> ust, hfx, qfx;
  View2<real> pblh;
  View3<real> tke, dlk, wu, wv, wt, wq, exch_h, exch_m;
  View3<real> rublten, rvblten, rthblten, rqvblten, rqcblten;
  real dt, cp, g;
  int nz;

  template <int K>
  CPPWRF_HD void column(int i, int j) const {
    using std::sqrt;
    const int kte = nz;
    // BOULAC's 1D arrays (kms:kme), zero at the start of the call as WRF sets them before the column loop (and every
    // value a column reads is set by that column)
    real z1D[K + 2], dz1D[K + 2], u1D[K + 2], v1D[K + 2], th1D[K + 2], q1D[K + 2], qc1D[K + 2], rho1D[K + 2];
    real rhoz1D[K + 2], tke1D[K + 2], th01D[K + 2], dlk1D[K + 2], dls1D[K + 2], exch1D[K + 2], sf1D[K + 2],
        vl1D[K + 2];
    real a_u1D[K + 2], a_v1D[K + 2], a_t1D[K + 2], a_q1D[K + 2], a_qc1D[K + 2], a_e1D[K + 2];
    real b_u1D[K + 2], b_v1D[K + 2], b_t1D[K + 2], b_q1D[K + 2], b_qc1D[K + 2], b_e1D[K + 2];
    real dlg1D[K + 2], dl_u1D[K + 2], sh1D[K + 2], bu1D[K + 2], wu1D[K + 2], wv1D[K + 2], wt1D[K + 2], wq1D[K + 2],
        wqc1D[K + 2], gamma1D[K + 2];
    for (int k = 0; k <= kte + 1; ++k) {
      z1D[k] = dz1D[k] = u1D[k] = v1D[k] = th1D[k] = q1D[k] = qc1D[k] = rho1D[k] = real(0.);
      rhoz1D[k] = tke1D[k] = th01D[k] = dlk1D[k] = dls1D[k] = exch1D[k] = real(0.);
      sf1D[k] = vl1D[k] = real(1.);
      a_u1D[k] = a_v1D[k] = a_t1D[k] = a_q1D[k] = a_qc1D[k] = a_e1D[k] = real(0.);
      b_u1D[k] = b_v1D[k] = b_t1D[k] = b_q1D[k] = b_qc1D[k] = b_e1D[k] = real(0.);
      dlg1D[k] = dl_u1D[k] = sh1D[k] = bu1D[k] = wu1D[k] = wv1D[k] = wt1D[k] = wq1D[k] = real(0.);
    }
    z1D[1] = real(0.);
    for (int iz = 1; iz <= kte; ++iz) {
      const int k = iz - 1;
      u1D[iz] = u_phy(i, j, k);
      v1D[iz] = v_phy(i, j, k);
      th1D[iz] = th_phy(i, j, k);
      q1D[iz] = qv(i, j, k);
      qc1D[iz] = qc(i, j, k);
      tke1D[iz] = tke(i, j, k);
      rho1D[iz] = rho(i, j, k);
      th01D[iz] = real(300.);  // th_0 = 300 K (BETA = g/th0)
      dz1D[iz] = dz8w(i, j, k);
      z1D[iz + 1] = z1D[iz] + dz1D[iz];
    }
    rhoz1D[1] = rho1D[1];
    for (int iz = 2; iz <= kte; ++iz)
      rhoz1D[iz] = (rho1D[iz] * dz1D[iz - 1] + rho1D[iz - 1] * dz1D[iz]) / (dz1D[iz - 1] + dz1D[iz]);
    rhoz1D[kte + 1] = rho1D[kte];
    for (int iz = 1; iz <= kte; ++iz) {  // .not. flag_bep
      a_e1D[iz] = real(0.);
      b_e1D[iz] = real(0.);
      dlg1D[iz] = (z1D[iz] + z1D[iz + 1]) / real(2.);
      dl_u1D[iz] = real(0.);
      vl1D[iz] = real(1.);
      sf1D[iz] = real(1.);
    }
    sf1D[kte + 1] = real(1.);

    // pbl_height
    const real pbl = pbl_height(kte, dz1D, z1D, th1D, q1D);
    pblh(i, j) = pbl;

    // wstar; gamma (igamma = 1: Troen and Mahrt)
    const real wts = c_max(real(0.), hfx(i, j) / rho1D[1] / cp);
    const real wstar = libm::pow(g * wts * pbl / th01D[1], real(1. / 3.));
    for (int k = 0; k <= kte + 1; ++k) gamma1D[k] = real(0.);
    for (int iz = 2; iz <= kte; ++iz) {
      if (z1D[iz] <= real(1.0) * pbl && wts > real(0.))
        gamma1D[iz] = real(10.) * wts / wstar / pbl;
      else
        gamma1D[iz] = real(0.);
    }

    // boulac1D: dissip_bougeault, length_bougeault, cdtur_bougeault, tke_bougeault, diff for the TKE
    {
      real dlu[K + 2], dld[K + 2], we[K + 2];
      dissip_bougeault(g, kte, z1D, dz1D, tke1D, dlu, dld, th1D, th01D);
      for (int iz = 1; iz <= kte; ++iz) {  // length_bougeault (dl_u = 0)
        dld[iz] = c_min(dld[iz], dlg1D[iz]);
        dls1D[iz] = sqrt(dlu[iz] * dld[iz]);
        dlk1D[iz] = c_min(dlu[iz], dld[iz]);
        if (dl_u1D[iz] > real(0.)) {
          dls1D[iz] = real(1.) / (real(1.) / dls1D[iz] + real(1.) / dl_u1D[iz]);
          dlk1D[iz] = real(1.) / (real(1.) / dlk1D[iz] + real(1.) / dl_u1D[iz]);
        }
      }
      // cdtur_bougeault
      exch1D[1] = real(0.);
      for (int iz = 2; iz <= kte; ++iz) {
        const real te_m = (tke1D[iz - 1] * dz1D[iz] + tke1D[iz] * dz1D[iz - 1]) / (dz1D[iz] + dz1D[iz - 1]);
        const real dlk_m = (dlk1D[iz - 1] * dz1D[iz] + dlk1D[iz] * dz1D[iz - 1]) / (dz1D[iz] + dz1D[iz - 1]);
        exch1D[iz] = ck_b * dlk_m * sqrt(te_m);
        exch1D[iz] = c_max(exch1D[iz], real(0.1));
      }
      exch1D[kte + 1] = real(0.1);
      // tke_bougeault: shear
      sh1D[1] = real(0.);
      for (int iz = 2; iz <= kte - 1; ++iz) {
        const real u2 = (dz1D[iz + 1] * u1D[iz] + dz1D[iz] * u1D[iz + 1]) / (dz1D[iz] + dz1D[iz + 1]);
        const real u1 = (dz1D[iz] * u1D[iz - 1] + dz1D[iz - 1] * u1D[iz]) / (dz1D[iz - 1] + dz1D[iz]);
        const real v2 = (dz1D[iz + 1] * v1D[iz] + dz1D[iz] * v1D[iz + 1]) / (dz1D[iz] + dz1D[iz + 1]);
        const real v1 = (dz1D[iz] * v1D[iz - 1] + dz1D[iz - 1] * v1D[iz]) / (dz1D[iz - 1] + dz1D[iz]);
        const real cdm = real(0.5) * (exch1D[iz] + exch1D[iz + 1]);
        const real du = (u2 - u1) / dz1D[iz], dv = (v2 - v1) / dz1D[iz];
        const real dumdz = du * du + dv * dv;
        sh1D[iz] = cdm * dumdz;
      }
      sh1D[kte] = real(0.);
      // buoy
      bu1D[1] = real(0.);
      for (int iz = 2; iz <= kte - 1; ++iz) {
        const real dtdz1 = real(2.) * (th1D[iz] - th1D[iz - 1]) / (dz1D[iz - 1] + dz1D[iz]);
        const real dtdz2 = real(2.) * (th1D[iz + 1] - th1D[iz]) / (dz1D[iz + 1] + dz1D[iz]);
        const real dtmdz = real(0.5) * (dtdz1 + dtdz2);
        const real cdm = real(0.5) * (exch1D[iz + 1] + exch1D[iz]);
        const real gammam = real(0.5) * (gamma1D[iz + 1] + gamma1D[iz]);
        bu1D[iz] = -cdm * (dtmdz - gammam) * g / th01D[iz];
      }
      bu1D[kte] = real(0.);
      for (int iz = 1; iz <= kte; ++iz) {
        const real te1 = c_max(tke1D[iz], temin);
        const real dl1 = c_max(dls1D[iz], real(0.1));
        const real td = -ceps_b * sqrt(te1) / dl1;
        sh1D[iz] = sh1D[iz] * sf1D[iz];
        bu1D[iz] = bu1D[iz] * sf1D[iz];
        a_e1D[iz] = a_e1D[iz] + td;
        b_e1D[iz] = b_e1D[iz] + sh1D[iz] + bu1D[iz];
      }
      // solve for tke
      diff<K>(kte, dt, tke1D, rho1D, rhoz1D, exch1D, a_e1D, b_e1D, sf1D, vl1D, dz1D, we);
      for (int iz = 1; iz <= kte; ++iz)
        if (tke1D[iz] < temin) tke1D[iz] = temin;
    }

    for (int iz = 1; iz <= kte; ++iz) {
      const int k = iz - 1;
      tke(i, j, k) = tke1D[iz];
      dlk(i, j, k) = dlk1D[iz];
      exch_h(i, j, k) = exch1D[iz];
      exch_m(i, j, k) = exch1D[iz];
    }

    // idiff /= 1: the tendencies (flag_bep = .false.)
    for (int iz = 1; iz <= kte; ++iz) {
      a_t1D[iz] = real(0.); b_t1D[iz] = real(0.); a_u1D[iz] = real(0.); b_u1D[iz] = real(0.);
      a_v1D[iz] = real(0.); b_v1D[iz] = real(0.); a_q1D[iz] = real(0.); b_q1D[iz] = real(0.);
    }
    b_t1D[1] = hfx(i, j) / dz1D[1] / rho1D[1] / cp;
    b_q1D[1] = qfx(i, j) / dz1D[1] / rho1D[1];
    const real ws1 = libm::pow(u1D[1] * u1D[1] + v1D[1] * v1D[1], real(.5));
    a_u1D[1] = (-ust(i, j) * ust(i, j) / dz1D[1] / ws1);
    a_v1D[1] = (-ust(i, j) * ust(i, j) / dz1D[1] / ws1);
    for (int iz = 2; iz <= kte; ++iz)
      if (z1D[iz] <= real(1.0) * pbl && wts > real(0.))
        b_t1D[iz] = b_t1D[iz] - (exch1D[iz + 1] * gamma1D[iz + 1] - exch1D[iz] * gamma1D[iz]) / dz1D[iz];

    diff<K>(kte, dt, u1D, rho1D, rhoz1D, exch1D, a_u1D, b_u1D, sf1D, vl1D, dz1D, wu1D);
    diff<K>(kte, dt, v1D, rho1D, rhoz1D, exch1D, a_v1D, b_v1D, sf1D, vl1D, dz1D, wv1D);
    diff<K>(kte, dt, th1D, rho1D, rhoz1D, exch1D, a_t1D, b_t1D, sf1D, vl1D, dz1D, wt1D);
    diff<K>(kte, dt, q1D, rho1D, rhoz1D, exch1D, a_q1D, b_q1D, sf1D, vl1D, dz1D, wq1D);
    diff<K>(kte, dt, qc1D, rho1D, rhoz1D, exch1D, a_qc1D, b_qc1D, sf1D, vl1D, dz1D, wqc1D);

    for (int iz = 1; iz <= kte; ++iz) {
      const int k = iz - 1;
      rthblten(i, j, k) = (th1D[iz] - th_phy(i, j, k)) / dt;
      rqvblten(i, j, k) = (q1D[iz] - qv(i, j, k)) / dt;
      rqcblten(i, j, k) = (qc1D[iz] - qc(i, j, k)) / dt;
      rublten(i, j, k) = (u1D[iz] - u_phy(i, j, k)) / dt;
      rvblten(i, j, k) = (v1D[iz] - v_phy(i, j, k)) / dt;
      wu(i, j, k) = wu1D[iz];
      wv(i, j, k) = wv1D[iz];
      wt(i, j, k) = wt1D[iz];
      wq(i, j, k) = wq1D[iz];
    }
  }
};

template <int K>
struct BoulacRunner {
  BoulacKernel k;
  CPPWRF_HD void operator()(int i, int j) const { k.template column<K>(i, j); }
};

}  // namespace cppwrf::phys::boulac
