#pragma once

#include "core/backend.hpp"
#include "core/constants.hpp"
#include "core/field.hpp"
#include "dyn/diffusion/diffusion_params.hpp"

namespace cppwrf::dyn::diffusion {

// Ported from WRF v4.8.0 dyn_em/module_diffusion_em.F:
//   horizontal_diffusion_2 (lines 2864-3113)
//   horizontal_diffusion_u_2 (lines 3118-3318)
//   horizontal_diffusion_v_2 (lines 3323-3514)
//   horizontal_diffusion_w_2 (lines 3519-3706)
//   horizontal_diffusion_s   (lines 3711-3999)
//   vertical_diffusion_2   (lines 4004-4458, bl_pbl_physics = 0)
//   vertical_diffusion_u_2 (lines 4463-4571)
//   vertical_diffusion_v_2 (lines 4576-4683)
//   vertical_diffusion_w_2 (lines 4688-4784)
//   vertical_diffusion_s   (lines 4789-4907)
//   cal_titau_11_22_33     (lines 5331-5451)
//   cal_titau_12_21        (lines 5456-5594)
//   cal_titau_13_31        (lines 5598-5744)
//   cal_titau_23_32        (lines 5749-5896)

constexpr real Gr = real(cst::g);

// -----------------------------------------------------------------------------
// Stress tensor calculations
// -----------------------------------------------------------------------------

struct KCalTitauDiag {
  View3<const real> rho, xkmh, defor;
  View3<real> titau;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    titau(i, j, k) = -rho(i, j, k) * xkmh(i, j, k) * defor(i, j, k);
  }
};

struct KCalTitau12 {
  View3<const real> rho, xkmh, defor12;
  View3<real> titau12;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    real rhoavg = real(0.25) * (rho(i - 1, j, k) + rho(i, j, k) +
                                rho(i - 1, j - 1, k) + rho(i, j - 1, k));
    real xkmhavg = real(0.25) * (xkmh(i - 1, j, k) + xkmh(i, j, k) +
                                 xkmh(i - 1, j - 1, k) + xkmh(i, j - 1, k));
    titau12(i, j, k) = -rhoavg * xkmhavg * defor12(i, j, k);
  }
};

struct KCalTitau13 {
  View3<const real> rho, xkmv, defor13;
  View1<const real> fnm, fnp;
  View3<real> titau13;
  int nz;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    real rhoavg, xkmvavg;
    if (k == 0) {
      rhoavg = real(0.5) * (rho(i - 1, j, 0) + rho(i, j, 0));
      xkmvavg = real(0.5) * (xkmv(i - 1, j, 0) + xkmv(i, j, 0));
    } else if (k == nz) {
      rhoavg = real(0.5) * (rho(i - 1, j, nz - 1) + rho(i, j, nz - 1));
      xkmvavg = real(0.5) * (xkmv(i - 1, j, nz - 1) + xkmv(i, j, nz - 1));
    } else {
      rhoavg = real(0.5) * (fnm(k) * (rho(i - 1, j, k) + rho(i, j, k)) +
                            fnp(k) * (rho(i - 1, j, k - 1) + rho(i, j, k - 1)));
      xkmvavg = real(0.5) * (fnm(k) * (xkmv(i - 1, j, k) + xkmv(i, j, k)) +
                             fnp(k) * (xkmv(i - 1, j, k - 1) + xkmv(i, j, k - 1)));
    }
    titau13(i, j, k) = -rhoavg * xkmvavg * defor13(i, j, k);
  }
};

struct KCalTitau23 {
  View3<const real> rho, xkmv, defor23;
  View1<const real> fnm, fnp;
  View3<real> titau23;
  int nz;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    real rhoavg, xkmvavg;
    if (k == 0) {
      rhoavg = real(0.5) * (rho(i, j - 1, 0) + rho(i, j, 0));
      xkmvavg = real(0.5) * (xkmv(i, j - 1, 0) + xkmv(i, j, 0));
    } else if (k == nz) {
      rhoavg = real(0.5) * (rho(i, j - 1, nz - 1) + rho(i, j, nz - 1));
      xkmvavg = real(0.5) * (xkmv(i, j - 1, nz - 1) + xkmv(i, j, nz - 1));
    } else {
      rhoavg = real(0.5) * (fnm(k) * (rho(i, j - 1, k) + rho(i, j, k)) +
                            fnp(k) * (rho(i, j - 1, k - 1) + rho(i, j, k - 1)));
      xkmvavg = real(0.5) * (fnm(k) * (xkmv(i, j - 1, k) + xkmv(i, j, k)) +
                             fnp(k) * (xkmv(i, j - 1, k - 1) + xkmv(i, j, k - 1)));
    }
    titau23(i, j, k) = -rhoavg * xkmvavg * defor23(i, j, k);
  }
};

// -----------------------------------------------------------------------------
// Horizontal diffusion 2
// -----------------------------------------------------------------------------

struct KHorizDiffU2 {
  View3<const real> titau11, titau12, zx, zy, rdzw;
  BView2<const real> msfux, msfuy;
  View1<const real> fnm, fnp, dnw;
  View3<real> ru_tendf;
  real rdx, rdy;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    real titau1avg[128], titau2avg[128], zx_at_u[128], zy_at_u[128];
    for (int k = 1; k < nz; ++k) {
      titau1avg[k] = real(0.5) * (fnm(k) * (titau11(i - 1, j, k) + titau11(i, j, k)) +
                                  fnp(k) * (titau11(i - 1, j, k - 1) + titau11(i, j, k - 1)));
      titau2avg[k] = real(0.5) * (fnm(k) * (titau12(i, j + 1, k) + titau12(i, j, k)) +
                                  fnp(k) * (titau12(i, j + 1, k - 1) + titau12(i, j, k - 1)));
      zx_at_u[k] = real(0.5) * (zx(i, j, k) + zx(i, j, k + 1));
      zy_at_u[k] = real(0.125) * (zy(i - 1, j, k) + zy(i, j, k) +
                                  zy(i - 1, j + 1, k) + zy(i, j + 1, k) +
                                  zy(i - 1, j, k + 1) + zy(i, j, k + 1) +
                                  zy(i - 1, j + 1, k + 1) + zy(i, j + 1, k + 1));
    }
    titau1avg[0] = real(0);
    titau1avg[nz] = real(0);
    titau2avg[0] = real(0);
    titau2avg[nz] = real(0);
    zx_at_u[0] = real(0.5) * (zx(i, j, 0) + zx(i, j, 1));
    zy_at_u[0] = real(0.125) * (zy(i - 1, j, 0) + zy(i, j, 0) +
                                zy(i - 1, j + 1, 0) + zy(i, j + 1, 0) +
                                zy(i - 1, j, 1) + zy(i, j, 1) +
                                zy(i - 1, j + 1, 1) + zy(i, j + 1, 1));

    real mrdx = msfux(i, j) * rdx;
    real mrdy = msfuy(i, j) * rdy;

    for (int k = 0; k < nz; ++k) {
      real tmpdz = real(0.5) * (real(1.0) / rdzw(i, j, k) + real(1.0) / rdzw(i - 1, j, k));
      real term_x = mrdx * (titau11(i, j, k) - titau11(i - 1, j, k));
      real term_y = mrdy * (titau12(i, j + 1, k) - titau12(i, j, k));
      real term_zx = msfux(i, j) * zx_at_u[k] * (titau1avg[k + 1] - titau1avg[k]) / tmpdz;
      real term_zy = msfuy(i, j) * zy_at_u[k] * (titau2avg[k + 1] - titau2avg[k]) / tmpdz;

      ru_tendf(i, j, k) += (Gr * tmpdz / dnw(k)) * (term_x + term_y - term_zx - term_zy);
    }
  }
};

struct KHorizDiffV2 {
  View3<const real> titau12, titau22, zx, zy, rdzw;
  BView2<const real> msfvx, msfvy;
  View1<const real> fnm, fnp, dnw;
  View3<real> rv_tendf;
  real rdx, rdy;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    real titau1avg[128], titau2avg[128], zx_at_v[128], zy_at_v[128];
    for (int k = 1; k < nz; ++k) {
      titau1avg[k] = real(0.5) * (fnm(k) * (titau12(i + 1, j, k) + titau12(i, j, k)) +
                                  fnp(k) * (titau12(i + 1, j, k - 1) + titau12(i, j, k - 1)));
      titau2avg[k] = real(0.5) * (fnm(k) * (titau22(i, j - 1, k) + titau22(i, j, k)) +
                                  fnp(k) * (titau22(i, j - 1, k - 1) + titau22(i, j, k - 1)));
      zx_at_v[k] = real(0.125) * (zx(i, j, k) + zx(i + 1, j, k) +
                                  zx(i, j - 1, k) + zx(i + 1, j - 1, k) +
                                  zx(i, j, k + 1) + zx(i + 1, j, k + 1) +
                                  zx(i, j - 1, k + 1) + zx(i + 1, j - 1, k + 1));
      zy_at_v[k] = real(0.5) * (zy(i, j, k) + zy(i, j, k + 1));
    }
    titau1avg[0] = real(0);
    titau1avg[nz] = real(0);
    titau2avg[0] = real(0);
    titau2avg[nz] = real(0);
    zx_at_v[0] = real(0.125) * (zx(i, j, 0) + zx(i + 1, j, 0) +
                                zx(i, j - 1, 0) + zx(i + 1, j - 1, 0) +
                                zx(i, j, 1) + zx(i + 1, j, 1) +
                                zx(i, j - 1, 1) + zx(i + 1, j - 1, 1));
    zy_at_v[0] = real(0.5) * (zy(i, j, 0) + zy(i, j, 1));

    real mrdx = msfvx(i, j) * rdx;
    real mrdy = msfvy(i, j) * rdy;

    for (int k = 0; k < nz; ++k) {
      real tmpdz = real(0.5) * (real(1.0) / rdzw(i, j, k) + real(1.0) / rdzw(i, j - 1, k));
      real term_y = mrdy * (titau22(i, j, k) - titau22(i, j - 1, k));
      real term_x = mrdx * (titau12(i + 1, j, k) - titau12(i, j, k));
      real term_zx = msfvx(i, j) * zx_at_v[k] * (titau1avg[k + 1] - titau1avg[k]) / tmpdz;
      real term_zy = msfvy(i, j) * zy_at_v[k] * (titau2avg[k + 1] - titau2avg[k]) / tmpdz;

      rv_tendf(i, j, k) += (Gr * tmpdz / dnw(k)) * (term_y + term_x - term_zx - term_zy);
    }
  }
};

struct KHorizDiffW2 {
  View3<const real> titau13, titau23, zx, zy, rdz;
  BView2<const real> msftx, msfty;
  View1<const real> dn;
  View3<real> rw_tendf;
  real rdx, rdy;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    real titau1avg[128], titau2avg[128], zx_at_w[128], zy_at_w[128];
    for (int k = 0; k < nz; ++k) {
      titau1avg[k] = real(0.25) * (titau13(i + 1, j, k + 1) + titau13(i, j, k + 1) +
                                   titau13(i + 1, j, k) + titau13(i, j, k));
      titau2avg[k] = real(0.25) * (titau23(i, j + 1, k + 1) + titau23(i, j, k + 1) +
                                   titau23(i, j + 1, k) + titau23(i, j, k));
      zx_at_w[k] = real(0.5) * (zx(i, j, k) + zx(i + 1, j, k));
      zy_at_w[k] = real(0.5) * (zy(i, j, k) + zy(i, j + 1, k));
    }
    titau1avg[nz] = real(0);
    titau2avg[nz] = real(0);

    real mrdx = msftx(i, j) * rdx;
    real mrdy = msfty(i, j) * rdy;

    for (int k = 1; k < nz; ++k) {
      real term_x = mrdx * (titau13(i + 1, j, k) - titau13(i, j, k));
      real term_y = mrdy * (titau23(i, j + 1, k) - titau23(i, j, k));
      real term_z = msfty(i, j) * rdz(i, j, k) *
                    (zx_at_w[k] * (titau1avg[k] - titau1avg[k - 1]) +
                     zy_at_w[k] * (titau2avg[k] - titau2avg[k - 1]));

      rw_tendf(i, j, k) += (Gr / (dn(k) * rdz(i, j, k))) * (term_x + term_y - term_z);
    }
  }
};

struct KHorizDiffS2 {
  View3<const real> var, xkhh, rho, zx, zy, rdz, rdzw;
  BView2<const real> msftx, msfty, msfux, msfvy;
  View1<const real> fnm, fnp, dn, dnw;
  real cf1, cf2, cf3;
  View3<real> tendf;
  real rdx, rdy;
  int nz;
  bool doing_tke = false;

  CPPWRF_HD void operator()(int i, int j) const {
    const int ktes1 = nz - 1;
    const real dnw_top = dnw(ktes1), dn_top = dn(ktes1);

    // H1 along x: need at u points (i and i+1)
    real H1_0[128], H1_1[128];
    for (int ii : {i, i + 1}) {
      real H1avg[128];
      for (int k = 1; k < nz; ++k) {
        H1avg[k] = real(0.5) * (fnm(k) * (var(ii - 1, j, k) + var(ii, j, k)) +
                                fnp(k) * (var(ii - 1, j, k - 1) + var(ii, j, k - 1)));
      }
      H1avg[0] = real(0.5) * (cf1 * var(ii, j, 0) + cf2 * var(ii, j, 1) + cf3 * var(ii, j, 2) +
                              cf1 * var(ii - 1, j, 0) + cf2 * var(ii - 1, j, 1) + cf3 * var(ii - 1, j, 2));
      H1avg[nz] = real(0.5) * (var(ii, j, nz - 1) + (var(ii, j, nz - 1) - var(ii, j, nz - 2)) * real(0.5) * dnw_top / dn_top +
                               var(ii - 1, j, nz - 1) + (var(ii - 1, j, nz - 1) - var(ii - 1, j, nz - 2)) * real(0.5) * dnw_top / dn_top);

      for (int k = 0; k < nz; ++k) {
        real xkxavg = real(0.5) * (xkhh(ii - 1, j, k) + xkhh(ii, j, k)) * real(0.5) * (rho(ii - 1, j, k) + rho(ii, j, k));
        real tmpzx = real(0.5) * (zx(ii, j, k) + zx(ii, j, k + 1));
        real rdzu = real(2.0) / (real(1.0) / rdzw(ii, j, k) + real(1.0) / rdzw(ii - 1, j, k));
        real h = -msfux(ii, j) * xkxavg * (rdx * (var(ii, j, k) - var(ii - 1, j, k)) -
                                           tmpzx * (H1avg[k + 1] - H1avg[k]) * rdzu);
        if (ii == i) H1_0[k] = h;
        else H1_1[k] = h;
      }
    }

    // H2 along y: need at v points (j and j+1)
    real H2_0[128], H2_1[128];
    for (int jj : {j, j + 1}) {
      real H2avg[128];
      for (int k = 1; k < nz; ++k) {
        H2avg[k] = real(0.5) * (fnm(k) * (var(i, jj - 1, k) + var(i, jj, k)) +
                                fnp(k) * (var(i, jj - 1, k - 1) + var(i, jj, k - 1)));
      }
      H2avg[0] = real(0.5) * (cf1 * var(i, jj, 0) + cf2 * var(i, jj, 1) + cf3 * var(i, jj, 2) +
                              cf1 * var(i, jj - 1, 0) + cf2 * var(i, jj - 1, 1) + cf3 * var(i, jj - 1, 2));
      H2avg[nz] = real(0.5) * (var(i, jj, nz - 1) + (var(i, jj, nz - 1) - var(i, jj, nz - 2)) * real(0.5) * dnw_top / dn_top +
                               var(i, jj - 1, nz - 1) + (var(i, jj - 1, nz - 1) - var(i, jj - 1, nz - 2)) * real(0.5) * dnw_top / dn_top);

      for (int k = 0; k < nz; ++k) {
        real xkxavg = real(0.5) * (xkhh(i, jj - 1, k) + xkhh(i, jj, k)) * real(0.5) * (rho(i, jj - 1, k) + rho(i, jj, k));
        real tmpzy = real(0.5) * (zy(i, jj, k) + zy(i, jj, k + 1));
        real rdzv = real(2.0) / (real(1.0) / rdzw(i, jj, k) + real(1.0) / rdzw(i, jj - 1, k));
        real h = -msfvy(i, jj) * xkxavg * (rdy * (var(i, jj, k) - var(i, jj - 1, k)) -
                                           tmpzy * (H2avg[k + 1] - H2avg[k]) * rdzv);
        if (jj == j) H2_0[k] = h;
        else H2_1[k] = h;
      }
    }

    real H1avg[128], H2avg[128];
    for (int k = 1; k < nz; ++k) {
      H1avg[k] = real(0.5) * (fnm(k) * (H1_1[k] + H1_0[k]) + fnp(k) * (H1_1[k - 1] + H1_0[k - 1]));
      H2avg[k] = real(0.5) * (fnm(k) * (H2_1[k] + H2_0[k]) + fnp(k) * (H2_1[k - 1] + H2_0[k - 1]));
    }
    H1avg[0] = real(0);
    H1avg[nz] = real(0);
    H2avg[0] = real(0);
    H2avg[nz] = real(0);

    real mrdx = msftx(i, j) * rdx;
    real mrdy = msfty(i, j) * rdy;

    for (int k = 0; k < nz; ++k) {
      real tmptendf = tendf(i, j, k);
      real zx_at_m = real(0.25) * (zx(i, j, k) + zx(i + 1, j, k) + zx(i, j, k + 1) + zx(i + 1, j, k + 1));
      real zy_at_m = real(0.25) * (zy(i, j, k) + zy(i, j + 1, k) + zy(i, j, k + 1) + zy(i, j + 1, k + 1));

      real term_x = mrdx * (H1_1[k] - H1_0[k]);
      real term_y = mrdy * (H2_1[k] - H2_0[k]);
      real term_zx = msftx(i, j) * zx_at_m * (H1avg[k + 1] - H1avg[k]) * rdzw(i, j, k);
      real term_zy = msfty(i, j) * zy_at_m * (H2avg[k + 1] - H2avg[k]) * rdzw(i, j, k);

      tendf(i, j, k) += (Gr / (dnw(k) * rdzw(i, j, k))) * (term_x + term_y - term_zx - term_zy);
      if (doing_tke) {
        tendf(i, j, k) = tmptendf + real(2.0) * (tendf(i, j, k) - tmptendf);
      }
    }
  }
};

// -----------------------------------------------------------------------------
// Vertical diffusion 2 (bl_pbl_physics = 0)
// -----------------------------------------------------------------------------

struct KVertDiffU2 {
  View3<const real> titau13;
  View1<const real> dnw;
  View3<real> ru_tendf;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    ru_tendf(i, j, 0) += (Gr / dnw(0)) * titau13(i, j, 1);
    for (int k = 1; k < nz; ++k) {
      ru_tendf(i, j, k) += (Gr / dnw(k)) * (titau13(i, j, k + 1) - titau13(i, j, k));
    }
  }
};

struct KVertDiffV2 {
  View3<const real> titau23;
  View1<const real> dnw;
  View3<real> rv_tendf;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    rv_tendf(i, j, 0) += (Gr / dnw(0)) * titau23(i, j, 1);
    for (int k = 1; k < nz; ++k) {
      rv_tendf(i, j, k) += (Gr / dnw(k)) * (titau23(i, j, k + 1) - titau23(i, j, k));
    }
  }
};

struct KVertDiffW2 {
  View3<const real> titau33;
  View1<const real> dn;
  View3<real> rw_tendf;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    for (int k = 1; k < nz; ++k) {
      rw_tendf(i, j, k) += Gr * (titau33(i, j, k) - titau33(i, j, k - 1)) / dn(k);
    }
  }
};

struct KVertDiffS2 {
  View3<const real> var, xkhv, rho, rdz;
  View1<const real> fnm, fnp, dnw;
  View3<real> tendf;
  int nz;
  bool doing_tke = false;

  CPPWRF_HD void operator()(int i, int j) const {
    real H3[128];
    H3[0] = real(0);
    H3[nz] = real(0);
    for (int k = 1; k < nz; ++k) {
      real xkxavg = (fnm(k) * xkhv(i, j, k) + fnp(k) * xkhv(i, j, k - 1)) *
                    (fnm(k) * rho(i, j, k) + fnp(k) * rho(i, j, k - 1));
      H3[k] = -xkxavg * (var(i, j, k) - var(i, j, k - 1)) * rdz(i, j, k);
    }
    for (int k = 0; k < nz; ++k) {
      real tmptendf = tendf(i, j, k);
      tendf(i, j, k) += Gr * (H3[k + 1] - H3[k]) / dnw(k);
      if (doing_tke) {
        tendf(i, j, k) = tmptendf + real(2.0) * (tendf(i, j, k) - tmptendf);
      }
    }
  }
};

struct KVertDiffS2Pert {
  View3<const real> var, xkhv, rho, rdz;
  View1<const real> base, fnm, fnp, dnw;
  View3<real> tendf;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    real H3[128];
    H3[0] = real(0);
    H3[nz] = real(0);
    for (int k = 1; k < nz; ++k) {
      real xkxavg = (fnm(k) * xkhv(i, j, k) + fnp(k) * xkhv(i, j, k - 1)) *
                    (fnm(k) * rho(i, j, k) + fnp(k) * rho(i, j, k - 1));
      // WRF: var_mix(k) = var(k) - base(k), then var_mix(k) - var_mix(k-1)
      real dpert = (var(i, j, k) - base(k)) - (var(i, j, k - 1) - base(k - 1));
      H3[k] = -xkxavg * dpert * rdz(i, j, k);
    }
    for (int k = 0; k < nz; ++k) {
      tendf(i, j, k) += Gr * (H3[k + 1] - H3[k]) / dnw(k);
    }
  }
};

// Surface momentum flux on U at k = 0 (module_diffusion_em.F:4161-4244)
struct KVertDiffSfcFluxU2 {
  View3<const real> u_2, v_2, rho;
  View2<const real> ust;
  View1<const real> dnw;
  bool has_ust, periodic_x;
  int isfflx, nx;
  real tke_drag_coefficient;
  View3<real> ru_tendf;

  CPPWRF_HD void operator()(int i, int j) const {
    constexpr real epsilon = real(1.0e-15);
    const int k   = 0;
    const int i0  = (i == nx) ? (periodic_x ? 0 : nx - 1) : i;
    const int im1 = (i == 0) ? (periodic_x ? nx - 1 : 0) : (i - 1);
    real v_avg  = (v_2(i0, j, k) + v_2(i0, j + 1, k) + v_2(im1, j, k) + v_2(im1, j + 1, k)) / real(4.0);
    real V0_u   = std::sqrt((u_2(i, j, k) * u_2(i, j, k)) + (v_avg * v_avg)) + epsilon;
    real tao_xz;
    if (isfflx == 0) {
      tao_xz = tke_drag_coefficient * V0_u * u_2(i, j, k);
    } else {
      real ustar = has_ust ? real(0.5) * (ust(i0, j) + ust(im1, j)) : real(0.0);
      tao_xz     = ustar * ustar * u_2(i, j, k) / V0_u;
    }
    ru_tendf(i, j, k) = ru_tendf(i, j, k) +
                        Gr * tao_xz * real(0.5) * (rho(i0, j, k) + rho(im1, j, k)) / dnw(k);
  }
};

// Surface momentum flux on V at k = 0 (module_diffusion_em.F:4185-4241)
struct KVertDiffSfcFluxV2 {
  View3<const real> u_2, v_2, rho;
  View2<const real> ust;
  View1<const real> dnw;
  bool has_ust, periodic_y;
  int isfflx, ny;
  real tke_drag_coefficient;
  View3<real> rv_tendf;

  CPPWRF_HD void operator()(int i, int j) const {
    constexpr real epsilon = real(1.0e-15);
    const int k   = 0;
    const int j0  = (j == ny) ? (periodic_y ? 0 : ny - 1) : j;
    const int jm1 = (j == 0) ? (periodic_y ? ny - 1 : 0) : (j - 1);
    real u_avg  = (u_2(i, j0, k) + u_2(i, jm1, k) + u_2(i + 1, j0, k) + u_2(i + 1, jm1, k)) / real(4.0);
    real V0_v   = std::sqrt((v_2(i, j, k) * v_2(i, j, k)) + (u_avg * u_avg)) + epsilon;
    real tao_yz;
    if (isfflx == 0) {
      tao_yz = tke_drag_coefficient * V0_v * v_2(i, j, k);
    } else {
      real ustar = has_ust ? real(0.5) * (ust(i, j0) + ust(i, jm1)) : real(0.0);
      tao_yz     = ustar * ustar * v_2(i, j, k) / V0_v;
    }
    rv_tendf(i, j, k) = rv_tendf(i, j, k) +
                        Gr * tao_yz * real(0.5) * (rho(i, j0, k) + rho(i, jm1, k)) / dnw(k);
  }
};

// Surface heat flux on theta at k = 0 (module_diffusion_em.F:4286-4325)
struct KVertDiffSfcFluxH2 {
  View3<const real> theta, rho, qv;
  View2<const real> hfx, qfx;
  View1<const real> dnw;
  bool has_qv, has_hfx, has_qfx;
  int isfflx, use_theta_m;
  real tke_heat_flux;
  View3<real> rt_tendf;

  CPPWRF_HD void operator()(int i, int j) const {
    constexpr real r_d    = real(287.0);
    constexpr real cp     = real(7.0) * r_d / real(2.0);
    constexpr real r_v    = real(461.6);
    constexpr real rvovrd = r_v / r_d;
    const int k = 0;
    real qv_0   = has_qv ? qv(i, j, k) : real(0.0);
    real qfx_ij = has_qfx ? qfx(i, j) : real(0.0);
    real cpm    = cp * (real(1.0) + real(0.8) * qv_0);

    if (isfflx == 0 || isfflx == 2) {
      real heat_flux = tke_heat_flux;
      if (hfx_out.p != nullptr) hfx_out(i, j) = heat_flux * cpm * rho(i, j, k);  // module_diffusion_em.F:4293
      if (use_theta_m == 1) {
        rt_tendf(i, j, k) = rt_tendf(i, j, k) -
                            Gr * heat_flux * (real(1.0) + rvovrd * qv_0) * rho(i, j, k) / dnw(k) -
                            Gr * real(1.61) * theta(i, j, k) * qfx_ij / dnw(k);
      } else {
        rt_tendf(i, j, k) = rt_tendf(i, j, k) -
                            Gr * heat_flux * rho(i, j, k) / dnw(k);
      }
    } else {
      real hfx_ij    = has_hfx ? hfx(i, j) : real(0.0);
      real heat_flux = hfx_ij / cpm;
      if (use_theta_m == 1) {
        rt_tendf(i, j, k) = rt_tendf(i, j, k) -
                            Gr * heat_flux * (real(1.0) + rvovrd * qv_0) / dnw(k) -
                            Gr * real(1.61) * theta(i, j, k) * qfx_ij / dnw(k);
      } else {
        rt_tendf(i, j, k) = rt_tendf(i, j, k) -
                            Gr * heat_flux / dnw(k);
      }
    }
  }
};

// Surface moisture flux on qv at k = 0 (module_diffusion_em.F:4384-4400)
struct KVertDiffSfcFluxQ2 {
  View2<const real> qfx;
  View1<const real> dnw;
  bool has_qfx;
  int isfflx;
  View3<real> qv_tendf;

  CPPWRF_HD void operator()(int i, int j) const {
    if (isfflx == 1 || isfflx == 2) {
      real moist_flux = has_qfx ? qfx(i, j) : real(0.0);
      qv_tendf(i, j, 0) = qv_tendf(i, j, 0) - Gr * moist_flux / dnw(0);
    }
  }
};

inline void horizontal_diffusion_2(const DiffusionParams& p,
                                   View3<const real> defor11,
                                   View3<const real> defor22,
                                   View3<const real> defor33,
                                   View3<const real> defor12,
                                   View3<const real> defor13,
                                   View3<const real> defor23,
                                   View3<const real> rho,
                                   View3<const real> xkmh,
                                   View3<const real> xkmv,
                                   View3<const real> xkhh,
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
                                   View3<real> titau11,
                                   View3<real> titau22,
                                   View3<real> titau33,
                                   View3<real> titau12,
                                   View3<real> titau13,
                                   View3<real> titau23,
                                   View3<real> ru_tendf,
                                   View3<real> rv_tendf,
                                   View3<real> rw_tendf,
                                   View3<const real> thp,
                                   View3<real> t_tendf,
                                   View3<const real> moist,
                                   View3<real> moist_tendf,
                                   bool has_moist,
                                   int nx, int ny, int nz,
                                   View3<const real> tke = {},
                                   View3<real> tke_tendf = {},
                                   bool has_tke = false) {
  const bool bdy_xs = p.open_xs || p.specified || p.nested;
  const bool bdy_xe = p.open_xe || p.specified || p.nested;
  const bool bdy_ys = p.open_ys || p.specified || p.nested;
  const bool bdy_ye = p.open_ye || p.specified || p.nested;

  // 1. Calculate stress tensor components
  parallel_for(Range3D{0, nx, 0, ny, 0, nz}, KCalTitauDiag{rho, xkmh, defor11, titau11});
  parallel_for(Range3D{0, nx, 0, ny, 0, nz}, KCalTitauDiag{rho, xkmh, defor22, titau22});
  parallel_for(Range3D{0, nx, 0, ny, 0, nz}, KCalTitauDiag{rho, xkmh, defor33, titau33});

  int i12_s = (!p.periodic_x && bdy_xs) ? 1 : 0;
  int i12_e = (!p.periodic_x && bdy_xe) ? nx : nx + 1;
  int j12_s = bdy_ys ? 1 : 0;
  int j12_e = bdy_ye ? ny : ny + 1;
  parallel_for(Range3D{i12_s, i12_e, j12_s, j12_e, 0, nz},
               KCalTitau12{rho, xkmh, defor12, titau12});

  int i13_s = (!p.periodic_x && bdy_xs) ? 1 : 0;
  int i13_e = (!p.periodic_x && bdy_xe) ? nx : nx + 1;
  int j13_s = bdy_ys ? 1 : 0;
  int j13_e = bdy_ye ? ny - 1 : ny;
  parallel_for(Range3D{i13_s, i13_e, j13_s, j13_e, 0, nz + 1},
               KCalTitau13{rho, xkmv, defor13, fnm, fnp, titau13, nz});

  int i23_s = (!p.periodic_x && bdy_xs) ? 1 : 0;
  int i23_e = (!p.periodic_x && bdy_xe) ? nx - 1 : nx;
  int j23_s = bdy_ys ? 1 : 0;
  int j23_e = bdy_ye ? ny : ny + 1;
  parallel_for(Range3D{i23_s, i23_e, j23_s, j23_e, 0, nz + 1},
               KCalTitau23{rho, xkmv, defor23, fnm, fnp, titau23, nz});

  // 2. Horizontal diffusion of U
  int iu_s = (!p.periodic_x && bdy_xs) ? 1 : 0;
  int iu_e = (!p.periodic_x && bdy_xe) ? nx : nx + 1;
  int ju_s = bdy_ys ? 1 : 0;
  int ju_e = bdy_ye ? ny - 1 : ny;
  column_for(Range2D{iu_s, iu_e, ju_s, ju_e},
             KHorizDiffU2{View3<const real>{titau11.p, titau11.sj, titau11.sk},
                          View3<const real>{titau12.p, titau12.sj, titau12.sk},
                          zx, zy, rdzw, msfux, msfuy, fnm, fnp, dnw,
                          ru_tendf, rdx, rdy, nz});

  // 3. Horizontal diffusion of V
  int iv_s = (!p.periodic_x && bdy_xs) ? 1 : 0;
  int iv_e = (!p.periodic_x && bdy_xe) ? nx - 1 : nx;
  int jv_s = bdy_ys ? 1 : 0;
  int jv_e = bdy_ye ? ny : ny + 1;
  column_for(Range2D{iv_s, iv_e, jv_s, jv_e},
             KHorizDiffV2{View3<const real>{titau12.p, titau12.sj, titau12.sk},
                          View3<const real>{titau22.p, titau22.sj, titau22.sk},
                          zx, zy, rdzw, msfvx, msfvy, fnm, fnp, dnw,
                          rv_tendf, rdx, rdy, nz});

  // 4. Horizontal diffusion of W
  int iw_s = (!p.periodic_x && bdy_xs) ? 1 : 0;
  int iw_e = (!p.periodic_x && bdy_xe) ? nx - 1 : nx;
  int jw_s = bdy_ys ? 1 : 0;
  int jw_e = bdy_ye ? ny - 1 : ny;
  column_for(Range2D{iw_s, iw_e, jw_s, jw_e},
             KHorizDiffW2{View3<const real>{titau13.p, titau13.sj, titau13.sk},
                          View3<const real>{titau23.p, titau23.sj, titau23.sk},
                          zx, zy, rdz, msftx, msfty, dn,
                          rw_tendf, rdx, rdy, nz});

  // 5. Horizontal diffusion of scalars (theta, TKE, and moisture)
  int is_s = (!p.periodic_x && bdy_xs) ? 1 : 0;
  int is_e = (!p.periodic_x && bdy_xe) ? nx - 1 : nx;
  int js_s = bdy_ys ? 1 : 0;
  int js_e = bdy_ye ? ny - 1 : ny;
  column_for(Range2D{is_s, is_e, js_s, js_e},
             KHorizDiffS2{thp, xkhh, rho, zx, zy, rdz, rdzw, msftx, msfty, msfux, msfvy,
                          fnm, fnp, dn, dnw, cf1, cf2, cf3, t_tendf, rdx, rdy, nz, false});

  if (has_tke && ((p.km_opt == 2 && !p.tke_mix2_off) || p.km_opt == 5)) {
    // Note: WRF passes xkmh (not xkhh) and doing_tke = .true. for TKE horizontal diffusion
    column_for(Range2D{is_s, is_e, js_s, js_e},
               KHorizDiffS2{tke, xkmh, rho, zx, zy, rdz, rdzw, msftx, msfty, msfux, msfvy,
                            fnm, fnp, dn, dnw, cf1, cf2, cf3, tke_tendf, rdx, rdy, nz, true});
  }

  if (has_moist && !p.moist_mix2_off) {
    column_for(Range2D{is_s, is_e, js_s, js_e},
               KHorizDiffS2{moist, xkhh, rho, zx, zy, rdz, rdzw, msftx, msfty, msfux, msfvy,
                            fnm, fnp, dn, dnw, cf1, cf2, cf3, moist_tendf, rdx, rdy, nz, false});
  }
}

inline void cal_titau_vert(View3<const real> rho,
                           View3<const real> xkmh,
                           View3<const real> xkmv,
                           View3<const real> defor13,
                           View3<const real> defor23,
                           View3<const real> defor33,
                           View1<const real> fnm,
                           View1<const real> fnp,
                           View3<real> titau13,
                           View3<real> titau23,
                           View3<real> titau33,
                           bool specified, int nx, int ny, int nz,
                           const DiffusionParams& p = {}) {
  const bool bdy_xs = p.open_xs || p.specified || p.nested || specified;
  const bool bdy_xe = p.open_xe || p.specified || p.nested || specified;
  const bool bdy_ys = p.open_ys || p.specified || p.nested || specified;
  const bool bdy_ye = p.open_ye || p.specified || p.nested || specified;
  parallel_for(Range3D{0, nx, 0, ny, 0, nz}, KCalTitauDiag{rho, xkmh, defor33, titau33});
  int i13_s = (!p.periodic_x && bdy_xs) ? 1 : 0;
  int i13_e = (!p.periodic_x && bdy_xe) ? nx : nx + 1;
  int j13_s = bdy_ys ? 1 : 0;
  int j13_e = bdy_ye ? ny - 1 : ny;
  parallel_for(Range3D{i13_s, i13_e, j13_s, j13_e, 0, nz + 1},
               KCalTitau13{rho, xkmv, defor13, fnm, fnp, titau13, nz});
  int i23_s = (!p.periodic_x && bdy_xs) ? 1 : 0;
  int i23_e = (!p.periodic_x && bdy_xe) ? nx - 1 : nx;
  int j23_s = bdy_ys ? 1 : 0;
  int j23_e = bdy_ye ? ny : ny + 1;
  parallel_for(Range3D{i23_s, i23_e, j23_s, j23_e, 0, nz + 1},
               KCalTitau23{rho, xkmv, defor23, fnm, fnp, titau23, nz});
}

inline void vertical_diffusion_2(const DiffusionParams& p,
                                 View3<const real> titau13,
                                 View3<const real> titau23,
                                 View3<const real> titau33,
                                 View3<const real> var_th,
                                 View1<const real> t_base,
                                 View3<const real> xkhv,
                                 View3<const real> rho,
                                 View3<const real> rdz,
                                 View1<const real> fnm,
                                 View1<const real> fnp,
                                 View1<const real> dn,
                                 View1<const real> dnw,
                                 View3<real> ru_tendf,
                                 View3<real> rv_tendf,
                                 View3<real> rw_tendf,
                                 View3<real> t_tendf,
                                 View3<const real> moist,
                                 View3<real> moist_tendf,
                                 bool has_moist,
                                 int nx, int ny, int nz,
                                 View3<const real> tke = {},
                                 View3<real> tke_tendf = {},
                                 View3<const real> xkmv = {},
                                 bool has_tke = false,
                                 View3<const real> u_2 = {},
                                 View3<const real> v_2 = {},
                                 View3<const real> theta = {},
                                 View1<const real> qv_base = {},
                                 bool has_qv_base = false,
                                 View2<const real> ust = {},
                                 bool has_ust = false,
                                 View2<const real> hfx = {},
                                 bool has_hfx = false,
                                 View2<const real> qfx = {},
                                 bool has_qfx = false) {
  const bool bdy_xs = p.open_xs || p.specified || p.nested;
  const bool bdy_xe = p.open_xe || p.specified || p.nested;
  const bool bdy_ys = p.open_ys || p.specified || p.nested;
  const bool bdy_ye = p.open_ye || p.specified || p.nested;

  // 1. Vertical diffusion of U
  int iu_s = (!p.periodic_x && bdy_xs) ? 1 : 0;
  int iu_e = (!p.periodic_x && bdy_xe) ? nx : nx + 1;
  int ju_s = bdy_ys ? 1 : 0;
  int ju_e = bdy_ye ? ny - 1 : ny;
  column_for(Range2D{iu_s, iu_e, ju_s, ju_e},
             KVertDiffU2{titau13, dnw, ru_tendf, nz});

  // 2. Vertical diffusion of V
  int iv_s = (!p.periodic_x && bdy_xs) ? 1 : 0;
  int iv_e = (!p.periodic_x && bdy_xe) ? nx - 1 : nx;
  int jv_s = bdy_ys ? 1 : 0;
  int jv_e = bdy_ye ? ny : ny + 1;
  column_for(Range2D{iv_s, iv_e, jv_s, jv_e},
             KVertDiffV2{titau23, dnw, rv_tendf, nz});

  // 3. Vertical diffusion of W
  int iw_s = (!p.periodic_x && bdy_xs) ? 1 : 0;
  int iw_e = (!p.periodic_x && bdy_xe) ? nx - 1 : nx;
  int jw_s = bdy_ys ? 1 : 0;
  int jw_e = bdy_ye ? ny - 1 : ny;
  column_for(Range2D{iw_s, iw_e, jw_s, jw_e},
             KVertDiffW2{titau33, dn, rw_tendf, nz});

  // 3b. Surface momentum fluxes (vflux) if u_2 and v_2 are provided
  if (u_2.p != nullptr && v_2.p != nullptr) {
    column_for(Range2D{0, nx + 1, 0, ny},
               KVertDiffSfcFluxU2{u_2, v_2, rho, ust, dnw, has_ust, p.periodic_x,
                                  p.isfflx, nx, p.tke_drag_coefficient, ru_tendf});
    column_for(Range2D{0, nx, 0, ny + 1},
               KVertDiffSfcFluxV2{u_2, v_2, rho, ust, dnw, has_ust, p.periodic_y,
                                  p.isfflx, ny, p.tke_drag_coefficient, rv_tendf});
  }

  // 4. Vertical diffusion of scalars (theta, TKE, and moisture)
  int is_s = (!p.periodic_x && bdy_xs) ? 1 : 0;
  int is_e = (!p.periodic_x && bdy_xe) ? nx - 1 : nx;
  int js_s = bdy_ys ? 1 : 0;
  int js_e = bdy_ye ? ny - 1 : ny;
  if (p.mix_full_fields) {
    column_for(Range2D{is_s, is_e, js_s, js_e},
               KVertDiffS2{var_th, xkhv, rho, rdz, fnm, fnp, dnw, t_tendf, nz, false});
  } else {
    column_for(Range2D{is_s, is_e, js_s, js_e},
               KVertDiffS2Pert{var_th, xkhv, rho, rdz, t_base, fnm, fnp, dnw, t_tendf, nz});
  }

  // 4b. Surface heat flux (hflux) if theta is provided (or when isfflx == 0/2 or has_hfx)
  if (theta.p != nullptr || p.isfflx == 0 || p.isfflx == 2 || has_hfx) {
    View3<const real> th_view = (theta.p != nullptr) ? theta : var_th;
    column_for(Range2D{0, nx, 0, ny},
               KVertDiffSfcFluxH2{th_view, rho, moist, hfx, qfx, dnw,
                                  has_moist, has_hfx, has_qfx,
                                  p.isfflx, p.use_theta_m, p.tke_heat_flux, t_tendf});
  }

  // 4c. Vertical diffusion of TKE (when km_opt == 2)
  if (has_tke && p.km_opt == 2) {
    // Note: WRF passes xkmv (not xkhv) and doing_tke = .true. for TKE vertical diffusion
    column_for(Range2D{is_s, is_e, js_s, js_e},
               KVertDiffS2{tke, xkmv, rho, rdz, fnm, fnp, dnw, tke_tendf, nz, true});
  }

  if (has_moist) {
    if (!p.mix_full_fields && has_qv_base) {
      column_for(Range2D{is_s, is_e, js_s, js_e},
                 KVertDiffS2Pert{moist, xkhv, rho, rdz, qv_base, fnm, fnp, dnw, moist_tendf, nz});
    } else {
      column_for(Range2D{is_s, is_e, js_s, js_e},
                 KVertDiffS2{moist, xkhv, rho, rdz, fnm, fnp, dnw, moist_tendf, nz, false});
    }
    if (has_qfx) {
      column_for(Range2D{0, nx, 0, ny},
                 KVertDiffSfcFluxQ2{qfx, dnw, has_qfx, p.isfflx, moist_tendf});
    }
  }
}

}  // namespace cppwrf::dyn::diffusion
