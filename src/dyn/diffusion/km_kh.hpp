#pragma once

#include <cmath>
#include "core/backend.hpp"
#include "core/field.hpp"
#include "dyn/diffusion/diffusion_params.hpp"

namespace cppwrf::dyn::diffusion {

// Ported from WRF v4.6.1 dyn_em/module_diffusion_em.F: calculate_km_kh (lines 1195-1349)
// Supports:
//   km_opt = 1: isotropic_km (lines 1718-1772)
//   km_opt = 4: smag2d_km (lines 1934-2044) with diff_opt = 2 terrain slope alpha reduction

struct KIsotropicKm {
  View3<real> xkmh, xkmv, xkhh, xkhv;
  real khdif, kvdif, prandtl;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    xkmh(i, j, k) = khdif;
    xkmv(i, j, k) = kvdif;
    xkhh(i, j, k) = khdif / prandtl;
    xkhv(i, j, k) = kvdif / prandtl;
  }
};

struct KSmag2dKm {
  View3<const real> defor11, defor22, defor12, rdzw, zx, zy;
  BView2<const real> msftx, msfty;
  View3<real> xkmh, xkmv, xkhh, xkhv;
  real dx, dy, c_s, prandtl;
  int diff_opt, nz;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    // def2 = 0.25 * (defor11 - defor22)^2 + (0.25 * sum(defor12))^2
    real d11_m_d22 = defor11(i, j, k) - defor22(i, j, k);
    real def2 = real(0.25) * (d11_m_d22 * d11_m_d22);
    real d12_avg = real(0.25) * (defor12(i, j, k) + defor12(i, j + 1, k) +
                                 defor12(i + 1, j, k) + defor12(i + 1, j + 1, k));
    def2 += d12_avg * d12_avg;

    real mlen_h = std::sqrt((dx / msftx(i, j)) * (dy / msfty(i, j)));
    real tmp = std::sqrt(def2);
    real kmh = c_s * c_s * mlen_h * mlen_h * tmp;
    kmh = kmh < real(10.0) * mlen_h ? kmh : real(10.0) * mlen_h;

    real kmv = kmh;
    real khh = kmh / prandtl;
    real khv = real(0.0);

    if (diff_opt == 2) {
      real dxm = dx / msftx(i, j);
      real dym = dy / msfty(i, j);
      real zx_kp1 = (k + 1 < nz) ? zx(i, j, k + 1) : real(0);
      real zxp1_kp1 = (k + 1 < nz) ? zx(i + 1, j, k + 1) : real(0);
      real tmpzx = real(0.25) * (std::abs(zx(i, j, k)) + std::abs(zx(i + 1, j, k)) +
                                 std::abs(zx_kp1) + std::abs(zxp1_kp1)) *
                   rdzw(i, j, k) * dxm;
      real zy_kp1 = (k + 1 < nz) ? zy(i, j, k + 1) : real(0);
      real zyp1_kp1 = (k + 1 < nz) ? zy(i, j + 1, k + 1) : real(0);
      real tmpzy = real(0.25) * (std::abs(zy(i, j, k)) + std::abs(zy(i, j + 1, k)) +
                                 std::abs(zy_kp1) + std::abs(zyp1_kp1)) *
                   rdzw(i, j, k) * dym;
      real slope = std::sqrt(tmpzx * tmpzx + tmpzy * tmpzy);
      real alpha = slope > real(1.0) ? slope : real(1.0);
      real def_limit = (real(10.0) / mlen_h) > real(1e-3) ? (real(10.0) / mlen_h) : real(1e-3);

      if (tmp > def_limit) {
        kmh = kmh / (alpha * alpha);
      } else {
        kmh = kmh / alpha;
      }
      khh = kmh / prandtl;
      kmv = kmh;
    }

    xkmh(i, j, k) = kmh;
    xkmv(i, j, k) = kmv;
    xkhh(i, j, k) = khh;
    xkhv(i, j, k) = khv;
  }
};

inline void calculate_km_kh(const DiffusionParams& p,
                            real dx, real dy,
                            View3<const real> defor11,
                            View3<const real> defor22,
                            View3<const real> defor12,
                            View3<const real> rdzw,
                            View3<const real> zx,
                            View3<const real> zy,
                            BView2<const real> msftx,
                            BView2<const real> msfty,
                            View3<real> xkmh,
                            View3<real> xkmv,
                            View3<real> xkhh,
                            View3<real> xkhv,
                            int nx, int ny, int nz) {
  if (p.km_opt == 1) {
    parallel_for(Range3D{0, nx, 0, ny, 0, nz},
                 KIsotropicKm{xkmh, xkmv, xkhh, xkhv, p.khdif, p.kvdif, p.prandtl});
  } else if (p.km_opt == 4) {
    int i_s = p.specified ? 1 : 0;
    int i_e = p.specified ? nx - 1 : nx;
    int j_s = p.specified ? 1 : 0;
    int j_e = p.specified ? ny - 1 : ny;

    parallel_for(Range3D{i_s, i_e, j_s, j_e, 0, nz},
                 KSmag2dKm{defor11, defor22, defor12, rdzw, zx, zy, msftx, msfty,
                           xkmh, xkmv, xkhh, xkhv, dx, dy, p.c_s, p.prandtl, p.diff_opt, nz});

    if (p.specified) {
      // Copy borders as needed for specified boundary conditions
      for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
          xkmh(0, j, k) = xkmh(1, j, k);
          xkmh(nx - 1, j, k) = xkmh(nx - 2, j, k);
          xkmv(0, j, k) = xkmv(1, j, k);
          xkmv(nx - 1, j, k) = xkmv(nx - 2, j, k);
          xkhh(0, j, k) = xkhh(1, j, k);
          xkhh(nx - 1, j, k) = xkhh(nx - 2, j, k);
          xkhv(0, j, k) = xkhv(1, j, k);
          xkhv(nx - 1, j, k) = xkhv(nx - 2, j, k);
        }
        for (int i = 0; i < nx; ++i) {
          xkmh(i, 0, k) = xkmh(i, 1, k);
          xkmh(i, ny - 1, k) = xkmh(i, ny - 2, k);
          xkmv(i, 0, k) = xkmv(i, 1, k);
          xkmv(i, ny - 1, k) = xkmv(i, ny - 2, k);
          xkhh(i, 0, k) = xkhh(i, 1, k);
          xkhh(i, ny - 1, k) = xkhh(i, ny - 2, k);
          xkhv(i, 0, k) = xkhv(i, 1, k);
          xkhv(i, ny - 1, k) = xkhv(i, ny - 2, k);
        }
      }
    }
  }
}

}  // namespace cppwrf::dyn::diffusion
