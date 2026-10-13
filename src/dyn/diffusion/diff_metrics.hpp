#pragma once

#include "core/backend.hpp"
#include "core/constants.hpp"
#include "core/field.hpp"
#include "dyn/diffusion/diffusion_params.hpp"

namespace cppwrf::dyn::diffusion {

// Ported from WRF v4.6.1 dyn_em/module_diffusion_em.F: compute_diff_metrics (lines 6882-7130)
// Computes physical vertical grid distances and metric slopes:
//   z_at_w = (ph + phb) / g
//   rdzw   = 1 / (z_at_w(k+1) - z_at_w(k)) at mass points
//   rdz    = 2 / (z_at_w(k+1) - z_at_w(k-1)) at w levels (and 2/(z1 - z0) at k=0)
//   z      = 0.5 * (ph(k) + phb(k) + ph(k+1) + phb(k+1)) / g at mass points
//   zx     = rdx * (phb(i) - phb(i-1) + ph(i) - ph(i-1)) / g at u points
//   zy     = rdy * (phb(j) - phb(j-1) + ph(j) - ph(j-1)) / g at v points

struct KDiffMetricsZ {
  View3<const real> ph;
  BView3<const real> phb;
  View3<real> z, rdz, rdzw;
  int nz;
  static constexpr real Gr = real(cst::g);

  CPPWRF_HD void operator()(int i, int j) const {
    // 128 max levels fits in stack registers
    real zw[128];
    for (int k = 0; k <= nz; ++k) {
      zw[k] = (ph(i, j, k) + phb(i, j, k)) / Gr;
    }
    for (int k = 0; k < nz; ++k) {
      rdzw(i, j, k) = real(1.0) / (zw[k + 1] - zw[k]);
      z(i, j, k) = real(0.5) * (ph(i, j, k) + phb(i, j, k) + ph(i, j, k + 1) + phb(i, j, k + 1)) / Gr;
    }
    rdz(i, j, 0) = real(2.0) / (zw[1] - zw[0]);
    for (int k = 1; k < nz; ++k) {
      rdz(i, j, k) = real(2.0) / (zw[k + 1] - zw[k - 1]);
    }
  }
};

struct KDiffMetricsZX {
  View3<const real> ph;
  BView3<const real> phb;
  View3<real> zx;
  real rdx;
  int nx;
  bool periodic_x;
  static constexpr real Gr = real(cst::g);

  CPPWRF_HD void operator()(int i, int j, int k) const {
    if (!periodic_x && (i == 0 || i == nx)) {
      zx(i, j, k) = real(0);
    } else {
      int i0 = (i == nx) ? 0 : i;
      int im1 = (i == 0) ? nx - 1 : (i - 1);
      real v = rdx * (phb(i0, j, k) - phb(im1, j, k)) / Gr;
      v += rdx * (ph(i0, j, k) - ph(im1, j, k)) / Gr;
      zx(i, j, k) = v;
    }
  }
};

struct KDiffMetricsZY {
  View3<const real> ph;
  BView3<const real> phb;
  View3<real> zy;
  real rdy;
  int ny;
  bool periodic_y;
  static constexpr real Gr = real(cst::g);

  CPPWRF_HD void operator()(int i, int j, int k) const {
    if (!periodic_y && (j == 0 || j == ny)) {
      zy(i, j, k) = real(0);
    } else {
      int j0 = (j == ny) ? 0 : j;
      int jm1 = (j == 0) ? ny - 1 : (j - 1);
      real v = rdy * (phb(i, j0, k) - phb(i, jm1, k)) / Gr;
      v += rdy * (ph(i, j0, k) - ph(i, jm1, k)) / Gr;
      zy(i, j, k) = v;
    }
  }
};

inline void compute_diff_metrics(const DiffusionParams& p,
                                 View3<const real> ph,
                                 BView3<const real> phb,
                                 View3<real> z,
                                 View3<real> rdz,
                                 View3<real> rdzw,
                                 View3<real> zx,
                                 View3<real> zy,
                                 real rdx, real rdy,
                                 int nx, int ny, int nz) {
  // Column-wise vertical metrics
  column_for(Range2D{0, nx, 0, ny}, KDiffMetricsZ{ph, phb, z, rdz, rdzw, nz});
  // zx at U points: i = 0..nx, j = 0..ny-1, k = 0..nz
  parallel_for(Range3D{0, nx + 1, 0, ny, 0, nz + 1},
               KDiffMetricsZX{ph, phb, zx, rdx, nx, p.periodic_x});
  // zy at V points: i = 0..nx-1, j = 0..ny, k = 0..nz
  parallel_for(Range3D{0, nx, 0, ny + 1, 0, nz + 1},
               KDiffMetricsZY{ph, phb, zy, rdy, ny, p.periodic_y});
}

}  // namespace cppwrf::dyn::diffusion
