#pragma once

#include "core/backend.hpp"
#include "core/field.hpp"
#include "dyn/diffusion/diffusion_params.hpp"

namespace cppwrf::dyn::diffusion {

// Ported from WRF v4.6.1 dyn_em/module_diffusion_em.F: cal_deform_and_div (lines 17-1190)
// Calculates deformation tensor components and 3-d divergence:
//   defor11 = 2 * du/dx (mass points)
//   defor22 = 2 * dv/dy (mass points)
//   defor33 = 2 * dw/dz (mass points)
//   div     = du/dx + dv/dy + dw/dz (mass points)
//   defor12 = du/dy + dv/dx (vorticity points: i=0..nx, j=0..ny)
//   defor13 = dw/dx + du/dz (u-w points: i=0..nx, j=0..ny-1, k=0..nz)
//   defor23 = dw/dy + dv/dz (v-w points: i=0..nx-1, j=0..ny, k=0..nz)

struct KDeformDiagAndDiv {
  View3<const real> u, v, w, zx, zy, rdzw;
  BView2<const real> msftx, msfty, msfuy, msfvx;
  View1<const real> fnm, fnp, dn, dnw;
  real cf1, cf2, cf3;
  View3<real> defor11, defor22, defor33, div;
  real rdx, rdy;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    const real mm = msftx(i, j) * msfty(i, j);
    const int ktes1 = nz - 2;
    const real cft2 = -real(0.5) * dnw(ktes1) / dn(ktes1);
    const real cft1 = real(1.0) - cft2;

    // --- du/dx ---
    real hat_u0[128], hat_u1[128], hatavg_u[128];
    for (int k = 0; k < nz; ++k) {
      hat_u0[k] = u(i, j, k) / msfuy(i, j);
      hat_u1[k] = u(i + 1, j, k) / msfuy(i + 1, j);
    }
    for (int k = 1; k < nz; ++k) {
      hatavg_u[k] = real(0.5) * (fnm(k) * (hat_u0[k] + hat_u1[k]) +
                                 fnp(k) * (hat_u0[k - 1] + hat_u1[k - 1]));
    }
    hatavg_u[0] = real(0.5) * (cf1 * hat_u0[0] + cf2 * hat_u0[1] + cf3 * hat_u0[2] +
                               cf1 * hat_u1[0] + cf2 * hat_u1[1] + cf3 * hat_u1[2]);
    hatavg_u[nz] = real(0.5) * (cft1 * (hat_u0[nz - 1] + hat_u1[nz - 1]) +
                                cft2 * (hat_u0[nz - 2] + hat_u1[nz - 2]));

    // --- dv/dy ---
    real hat_v0[128], hat_v1[128], hatavg_v[128];
    for (int k = 0; k < nz; ++k) {
      hat_v0[k] = v(i, j, k) / msfvx(i, j);
      hat_v1[k] = v(i, j + 1, k) / msfvx(i, j + 1);
    }
    for (int k = 1; k < nz; ++k) {
      hatavg_v[k] = real(0.5) * (fnm(k) * (hat_v0[k] + hat_v1[k]) +
                                 fnp(k) * (hat_v0[k - 1] + hat_v1[k - 1]));
    }
    hatavg_v[0] = real(0.5) * (cf1 * hat_v0[0] + cf2 * hat_v0[1] + cf3 * hat_v0[2] +
                               cf1 * hat_v1[0] + cf2 * hat_v1[1] + cf3 * hat_v1[2]);
    hatavg_v[nz] = real(0.5) * (cft1 * (hat_v0[nz - 1] + hat_v1[nz - 1]) +
                                cft2 * (hat_v0[nz - 2] + hat_v1[nz - 2]));

    for (int k = 0; k < nz; ++k) {
      real zx_kp1 = (k + 1 < nz) ? zx(i, j, k + 1) : real(0);
      real zxp1_kp1 = (k + 1 < nz) ? zx(i + 1, j, k + 1) : real(0);
      real tmpzx = real(0.25) * (zx(i, j, k) + zx(i + 1, j, k) +
                                 zx_kp1 + zxp1_kp1);
      real tmp1_u = (hatavg_u[k + 1] - hatavg_u[k]) * tmpzx * rdzw(i, j, k);
      real dudx = mm * (rdx * (hat_u1[k] - hat_u0[k]) - tmp1_u);
      defor11(i, j, k) = real(2.0) * dudx;

      real zy_kp1 = (k + 1 < nz) ? zy(i, j, k + 1) : real(0);
      real zyp1_kp1 = (k + 1 < nz) ? zy(i, j + 1, k + 1) : real(0);
      real tmpzy = real(0.25) * (zy(i, j, k) + zy(i, j + 1, k) +
                                 zy_kp1 + zyp1_kp1);
      real tmp1_v = (hatavg_v[k + 1] - hatavg_v[k]) * tmpzy * rdzw(i, j, k);
      real dvdy = mm * (rdy * (hat_v1[k] - hat_v0[k]) - tmp1_v);
      defor22(i, j, k) = real(2.0) * dvdy;

      real dwdz = (w(i, j, k + 1) - w(i, j, k)) * rdzw(i, j, k);
      defor33(i, j, k) = real(2.0) * dwdz;

      div(i, j, k) = dudx + dvdy + dwdz;
    }
  }
};

struct KDeform12 {
  View3<const real> u, v, zx, zy, rdzw;
  BView2<const real> msfux, msfvy;
  View1<const real> fnm, fnp, dn, dnw;
  real cf1, cf2, cf3;
  View3<real> defor12;
  real rdx, rdy;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    const real mm = real(0.25) * (msfux(i, j - 1) + msfux(i, j)) * (msfvy(i - 1, j) + msfvy(i, j));
    const int ktes1 = nz - 2;
    const real cft2 = -real(0.5) * dnw(ktes1) / dn(ktes1);
    const real cft1 = real(1.0) - cft2;

    // --- u / msfux at u points (j-1 and j) ---
    real hat_um[128], hat_up[128], hatavg_u[128];
    for (int k = 0; k < nz; ++k) {
      hat_um[k] = u(i, j - 1, k) / msfux(i, j - 1);
      hat_up[k] = u(i, j, k) / msfux(i, j);
    }
    for (int k = 1; k < nz; ++k) {
      hatavg_u[k] = real(0.5) * (fnm(k) * (hat_um[k] + hat_up[k]) +
                                 fnp(k) * (hat_um[k - 1] + hat_up[k - 1]));
    }
    hatavg_u[0] = real(0.5) * (cf1 * hat_um[0] + cf2 * hat_um[1] + cf3 * hat_um[2] +
                               cf1 * hat_up[0] + cf2 * hat_up[1] + cf3 * hat_up[2]);
    hatavg_u[nz] = real(0.5) * (cft1 * (hat_um[nz - 1] + hat_up[nz - 1]) +
                                cft2 * (hat_um[nz - 2] + hat_up[nz - 2]));

    // --- v / msfvy at v points (i-1 and i) ---
    real hat_vm[128], hat_vp[128], hatavg_v[128];
    for (int k = 0; k < nz; ++k) {
      hat_vm[k] = v(i - 1, j, k) / msfvy(i - 1, j);
      hat_vp[k] = v(i, j, k) / msfvy(i, j);
    }
    for (int k = 1; k < nz; ++k) {
      hatavg_v[k] = real(0.5) * (fnm(k) * (hat_vm[k] + hat_vp[k]) +
                                 fnp(k) * (hat_vm[k - 1] + hat_vp[k - 1]));
    }
    hatavg_v[0] = real(0.5) * (cf1 * hat_vm[0] + cf2 * hat_vm[1] + cf3 * hat_vm[2] +
                               cf1 * hat_vp[0] + cf2 * hat_vp[1] + cf3 * hat_vp[2]);
    hatavg_v[nz] = real(0.5) * (cft1 * (hat_vm[nz - 1] + hat_vp[nz - 1]) +
                                cft2 * (hat_vm[nz - 2] + hat_vp[nz - 2]));

    for (int k = 0; k < nz; ++k) {
      real zy_kp1 = (k + 1 < nz) ? zy(i - 1, j, k + 1) : real(0);
      real zyp1_kp1 = (k + 1 < nz) ? zy(i, j, k + 1) : real(0);
      real tmpzy = real(0.25) * (zy(i - 1, j, k) + zy(i, j, k) +
                                 zy_kp1 + zyp1_kp1);
      real rdzw_avg = real(0.25) * (rdzw(i, j, k) + rdzw(i - 1, j, k) +
                                    rdzw(i - 1, j - 1, k) + rdzw(i, j - 1, k));
      real tmp1_u = (hatavg_u[k + 1] - hatavg_u[k]) * tmpzy * rdzw_avg;
      real term_u = mm * (rdy * (hat_up[k] - hat_um[k]) - tmp1_u);

      real zx_kp1 = (k + 1 < nz) ? zx(i, j - 1, k + 1) : real(0);
      real zxp1_kp1 = (k + 1 < nz) ? zx(i, j, k + 1) : real(0);
      real tmpzx = real(0.25) * (zx(i, j - 1, k) + zx(i, j, k) +
                                 zx_kp1 + zxp1_kp1);
      real tmp1_v = (hatavg_v[k + 1] - hatavg_v[k]) * tmpzx * rdzw_avg;
      real term_v = mm * (rdx * (hat_vp[k] - hat_vm[k]) - tmp1_v);

      defor12(i, j, k) = term_u + term_v;
    }
  }
};

struct KDeform13 {
  View3<const real> u, w, zx, rdz;
  BView2<const real> msfux, msfuy, msfty;
  View3<real> defor13;
  real rdx;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    const real mm = msfux(i, j) * msfuy(i, j);

    real hat_wm[128], hat_wp[128], hatavg[128];
    for (int k = 0; k <= nz; ++k) {
      hat_wm[k] = w(i - 1, j, k) / msfty(i - 1, j);
      hat_wp[k] = w(i, j, k) / msfty(i, j);
    }
    for (int k = 0; k < nz; ++k) {
      hatavg[k] = real(0.25) * (hat_wp[k] + hat_wp[k + 1] + hat_wm[k] + hat_wm[k + 1]);
    }

    defor13(i, j, 0) = real(0);
    for (int k = 1; k < nz; ++k) {
      real rdz_avg = real(0.5) * (rdz(i, j, k) + rdz(i - 1, j, k));
      real tmp1_w = (hatavg[k] - hatavg[k - 1]) * zx(i, j, k) * rdz_avg;
      real dwdx = mm * (rdx * (hat_wp[k] - hat_wm[k]) - tmp1_w);
      real dudz = (u(i, j, k) - u(i, j, k - 1)) * rdz_avg;
      defor13(i, j, k) = dwdx + dudz;
    }
    defor13(i, j, nz) = real(0);
  }
};

struct KDeform23 {
  View3<const real> v, w, zy, rdz;
  BView2<const real> msfvx, msfvy, msftx;
  View3<real> defor23;
  real rdy;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    const real mm = msfvx(i, j) * msfvy(i, j);

    real hat_wm[128], hat_wp[128], hatavg[128];
    for (int k = 0; k <= nz; ++k) {
      hat_wm[k] = w(i, j - 1, k) / msftx(i, j - 1);
      hat_wp[k] = w(i, j, k) / msftx(i, j);
    }
    for (int k = 0; k < nz; ++k) {
      hatavg[k] = real(0.25) * (hat_wp[k] + hat_wp[k + 1] + hat_wm[k] + hat_wm[k + 1]);
    }

    defor23(i, j, 0) = real(0);
    for (int k = 1; k < nz; ++k) {
      real rdz_avg = real(0.5) * (rdz(i, j, k) + rdz(i, j - 1, k));
      real tmp1_w = (hatavg[k] - hatavg[k - 1]) * zy(i, j, k) * rdz_avg;
      real dwdy = mm * (rdy * (hat_wp[k] - hat_wm[k]) - tmp1_w);
      real dvdz = (v(i, j, k) - v(i, j, k - 1)) * rdz_avg;
      defor23(i, j, k) = dwdy + dvdz;
    }
    defor23(i, j, nz) = real(0);
  }
};

inline void cal_deform_and_div(const DiffusionParams& p,
                               View3<const real> u,
                               View3<const real> v,
                               View3<const real> w,
                               View3<const real> zx,
                               View3<const real> zy,
                               View3<const real> rdz,
                               View3<const real> rdzw,
                               BView2<const real> msftx,
                               BView2<const real> msfty,
                               BView2<const real> msfux,
                               BView2<const real> msfuy,
                               BView2<const real> msfvx,
                               BView2<const real> msfvy,
                               View1<const real> fnm,
                               View1<const real> fnp,
                               View1<const real> dn,
                               View1<const real> dnw,
                               real cf1, real cf2, real cf3,
                               real rdx, real rdy,
                               View3<real> defor11,
                               View3<real> defor22,
                               View3<real> defor33,
                               View3<real> defor12,
                               View3<real> defor13,
                               View3<real> defor23,
                               View3<real> div,
                               int nx, int ny, int nz) {
  // 1. Diagonal deformations and 3D divergence at mass points (0..nx-1, 0..ny-1)
  column_for(Range2D{0, nx, 0, ny},
             KDeformDiagAndDiv{u, v, w, zx, zy, rdzw, msftx, msfty, msfuy, msfvx,
                               fnm, fnp, dn, dnw, cf1, cf2, cf3, defor11, defor22, defor33, div,
                               rdx, rdy, nz});

  // 2. defor12 at vorticity points (staggered in x and y)
  int i12_s = p.specified ? 1 : 0;
  int i12_e = p.specified ? nx : nx + 1;
  int j12_s = p.specified ? 1 : 0;
  int j12_e = p.specified ? ny : ny + 1;

  column_for(Range2D{i12_s, i12_e, j12_s, j12_e},
             KDeform12{u, v, zx, zy, rdzw, msfux, msfvy,
                       fnm, fnp, dn, dnw, cf1, cf2, cf3, defor12,
                       rdx, rdy, nz});

  if (p.specified) {
    // Copy edges as WRF does in cal_deform_and_div
    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j <= ny; ++j) {
        defor12(0, j, k) = defor12(1, j, k);
        defor12(nx, j, k) = defor12(nx - 1, j, k);
      }
      for (int i = 0; i <= nx; ++i) {
        defor12(i, 0, k) = defor12(i, 1, k);
        defor12(i, ny, k) = defor12(i, ny - 1, k);
      }
    }
  }

  // 3. defor13 at u-w points (staggered in x and z: 0..nx, 0..ny-1)
  int i13_s = p.specified ? 1 : 0;
  int i13_e = p.specified ? nx : nx + 1;
  int j13_s = p.specified ? 1 : 0;
  int j13_e = p.specified ? ny - 1 : ny;

  column_for(Range2D{i13_s, i13_e, j13_s, j13_e},
             KDeform13{u, w, zx, rdz, msfux, msfuy, msfty, defor13, rdx, nz});

  if (p.specified) {
    for (int k = 0; k <= nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        defor13(0, j, k) = defor13(1, j, k);
        defor13(nx, j, k) = defor13(nx - 1, j, k);
      }
      for (int i = 0; i <= nx; ++i) {
        defor13(i, 0, k) = defor13(i, 1, k);
        defor13(i, ny - 1, k) = defor13(i, ny - 2, k);
      }
    }
  }

  // 4. defor23 at v-w points (staggered in y and z: 0..nx-1, 0..ny)
  int i23_s = p.specified ? 1 : 0;
  int i23_e = p.specified ? nx - 1 : nx;
  int j23_s = p.specified ? 1 : 0;
  int j23_e = p.specified ? ny : ny + 1;

  column_for(Range2D{i23_s, i23_e, j23_s, j23_e},
             KDeform23{v, w, zy, rdz, msfvx, msfvy, msftx, defor23, rdy, nz});

  if (p.specified) {
    for (int k = 0; k <= nz; ++k) {
      for (int j = 0; j <= ny; ++j) {
        defor23(0, j, k) = defor23(1, j, k);
        defor23(nx - 1, j, k) = defor23(nx - 2, j, k);
      }
      for (int i = 0; i < nx; ++i) {
        defor23(i, 0, k) = defor23(i, 1, k);
        defor23(i, ny, k) = defor23(i, ny - 1, k);
      }
    }
  }
}

}  // namespace cppwrf::dyn::diffusion
