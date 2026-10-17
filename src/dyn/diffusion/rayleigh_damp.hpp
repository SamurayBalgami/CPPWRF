#pragma once

#include <cmath>
#include <type_traits>
#include "core/backend.hpp"
#include "core/constants.hpp"
#include "core/field.hpp"
#include "dyn/diffusion/diffusion_params.hpp"

namespace cppwrf::dyn::diffusion {

// Ported from WRF v4.6.1 dyn_em/module_big_step_utilities_em.F:
//   rk_rayleigh_damp (lines 5836-6110, damp_opt = 2)

CPPWRF_HD inline real vsin_rayleigh(real x, int lane_k, int nk) {
#if defined(CPPWRF_HAVE_LIBMVEC)
  if constexpr (std::is_same_v<real, double>) {
    return cppwrf::vm::sin(x, lane_k, nk);
  }
#endif
  (void)lane_k;
  (void)nk;
  return std::sin(x);
}

struct KRayleighDampU {
  View3<const real> u, ph;
  BView3<const real> phb;
  View2<const real> muu;
  View1<const real> c1h, c2h, u_base, z_base;
  View3<real> ru_tendf;
  real dampcoef, zdamp, pii;
  int nx, nz;
  bool periodic_x;

  CPPWRF_HD void operator()(int i, int j) const {
    constexpr real g = real(9.81);
    const int ktf    = nz - 1;
    const int i0     = (i == nx) ? (periodic_x ? 0 : nx - 1) : i;
    const int im1    = (i == 0)  ? (periodic_x ? nx - 1 : 0) : (i - 1);

    real ztop = real(0.5) * (phb(i0, j, nz) + phb(im1, j, nz) +
                             ph(i0, j, nz) + ph(im1, j, nz)) / g;

    real z00[128], u00[128];
    int k1 = ktf;
    real z = ztop;
    while (k1 >= 0 && z >= (ztop - zdamp)) {
      z = real(0.25) * (phb(i0, j, k1) + phb(i0, j, k1 + 1) +
                        phb(im1, j, k1) + phb(im1, j, k1 + 1) +
                        ph(i0, j, k1) + ph(i0, j, k1 + 1) +
                        ph(im1, j, k1) + ph(im1, j, k1 + 1)) / g;
      z00[k1] = z;
      k1 = k1 - 1;
    }
    k1 = k1 + 2;

    for (int k = k1; k <= ktf; ++k) {
      int k2 = ktf;
      while (k2 > 0 && z_base(k2) > z00[k]) {
        k2 = k2 - 1;
      }
      if (k2 + 1 > ktf) {
        u00[k] = u_base(k2) + (u_base(k2) - u_base(k2 - 1)) *
                              (z00[k] - z_base(k2)) /
                              (z_base(k2) - z_base(k2 - 1));
      } else {
        u00[k] = u_base(k2) + (u_base(k2 + 1) - u_base(k2)) *
                              (z00[k] - z_base(k2)) /
                              (z_base(k2 + 1) - z_base(k2));
      }
    }

    const int nk = ktf - k1 + 1;
    for (int k = k1; k <= ktf; ++k) {
      real ratio = (ztop - z00[k]) / zdamp;
      real m     = ratio < real(1.0) ? ratio : real(1.0);
      real dcoef = real(1.0) - m;
      real s     = vsin_rayleigh(real(0.5) * pii * dcoef, k - k1, nk);
      dcoef      = s * s;
      ru_tendf(i, j, k) = ru_tendf(i, j, k) -
                          (c1h(k) * muu(i, j) + c2h(k)) * (dcoef * dampcoef) *
                          (u(i, j, k) - u00[k]);
    }
  }
};

struct KRayleighDampV {
  View3<const real> v, ph;
  BView3<const real> phb;
  View2<const real> muv;
  View1<const real> c1h, c2h, v_base, z_base;
  View3<real> rv_tendf;
  real dampcoef, zdamp, pii;
  int ny, nz;
  bool periodic_y;

  CPPWRF_HD void operator()(int i, int j) const {
    constexpr real g = real(9.81);
    const int ktf    = nz - 1;
    const int j0     = (j == ny) ? (periodic_y ? 0 : ny - 1) : j;
    const int jm1    = (j == 0)  ? (periodic_y ? ny - 1 : 0) : (j - 1);

    real ztop = real(0.5) * (phb(i, j0, nz) + phb(i, jm1, nz) +
                             ph(i, j0, nz) + ph(i, jm1, nz)) / g;

    real z00[128], v00[128];
    int k1 = ktf;
    real z = ztop;
    while (k1 >= 0 && z >= (ztop - zdamp)) {
      z = real(0.25) * (phb(i, j0, k1) + phb(i, j0, k1 + 1) +
                        phb(i, jm1, k1) + phb(i, jm1, k1 + 1) +
                        ph(i, j0, k1) + ph(i, j0, k1 + 1) +
                        ph(i, jm1, k1) + ph(i, jm1, k1 + 1)) / g;
      z00[k1] = z;
      k1 = k1 - 1;
    }
    k1 = k1 + 2;

    for (int k = k1; k <= ktf; ++k) {
      int k2 = ktf;
      while (k2 > 0 && z_base(k2) > z00[k]) {
        k2 = k2 - 1;
      }
      if (k2 + 1 > ktf) {
        v00[k] = v_base(k2) + (v_base(k2) - v_base(k2 - 1)) *
                              (z00[k] - z_base(k2)) /
                              (z_base(k2) - z_base(k2 - 1));
      } else {
        v00[k] = v_base(k2) + (v_base(k2 + 1) - v_base(k2)) *
                              (z00[k] - z_base(k2)) /
                              (z_base(k2 + 1) - z_base(k2));
      }
    }

    const int nk = ktf - k1 + 1;
    for (int k = k1; k <= ktf; ++k) {
      real ratio = (ztop - z00[k]) / zdamp;
      real m     = ratio < real(1.0) ? ratio : real(1.0);
      real dcoef = real(1.0) - m;
      real s     = vsin_rayleigh(real(0.5) * pii * dcoef, k - k1, nk);
      dcoef      = s * s;
      rv_tendf(i, j, k) = rv_tendf(i, j, k) -
                          (c1h(k) * muv(i, j) + c2h(k)) * (dcoef * dampcoef) *
                          (v(i, j, k) - v00[k]);
    }
  }
};

struct KRayleighDampWT {
  View3<const real> w, t, ph;
  BView3<const real> phb;
  View2<const real> mut;
  View1<const real> c1h, c2h, c1f, c2f, t_base, z_base;
  View3<real> rw_tendf, t_tendf;
  real dampcoef, zdamp, pii;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    constexpr real g = real(9.81);
    const int ktf    = nz - 1;

    // 1. Adjust w to base state (scalar sin in WRF)
    real ztop = (phb(i, j, nz) + ph(i, j, nz)) / g;
    for (int k = 0; k <= nz; ++k) {
      real z = (phb(i, j, k) + ph(i, j, k)) / g;
      if (z >= (ztop - zdamp)) {
        real ratio = (ztop - z) / zdamp;
        real m     = ratio < real(1.0) ? ratio : real(1.0);
        real dcoef = real(1.0) - m;
        real s     = std::sin(real(0.5) * pii * dcoef);
        dcoef      = s * s;
        rw_tendf(i, j, k) = rw_tendf(i, j, k) -
                            (c1f(k) * mut(i, j) + c2f(k)) * (dcoef * dampcoef) * w(i, j, k);
      }
    }

    // 2. Adjust potential temperature to base state (vectorized sin over k = k1..ktf)
    real z00[128], t00[128];
    int k1 = ktf;
    real z = ztop;
    while (k1 >= 0 && z >= (ztop - zdamp)) {
      z = real(0.5) * (phb(i, j, k1) + phb(i, j, k1 + 1) +
                       ph(i, j, k1) + ph(i, j, k1 + 1)) / g;
      z00[k1] = z;
      k1 = k1 - 1;
    }
    k1 = k1 + 2;

    for (int k = k1; k <= ktf; ++k) {
      int k2 = ktf;
      while (k2 > 0 && z_base(k2) > z00[k]) {
        k2 = k2 - 1;
      }
      if (k2 + 1 > ktf) {
        t00[k] = t_base(k2) + (t_base(k2) - t_base(k2 - 1)) *
                              (z00[k] - z_base(k2)) /
                              (z_base(k2) - z_base(k2 - 1));
      } else {
        t00[k] = t_base(k2) + (t_base(k2 + 1) - t_base(k2)) *
                              (z00[k] - z_base(k2)) /
                              (z_base(k2 + 1) - z_base(k2));
      }
    }

    const int nk = ktf - k1 + 1;
    for (int k = k1; k <= ktf; ++k) {
      real ratio = (ztop - z00[k]) / zdamp;
      real m     = ratio < real(1.0) ? ratio : real(1.0);
      real dcoef = real(1.0) - m;
      real s     = vsin_rayleigh(real(0.5) * pii * dcoef, k - k1, nk);
      dcoef      = s * s;
      t_tendf(i, j, k) = t_tendf(i, j, k) -
                         (c1h(k) * mut(i, j) + c2h(k)) * (dcoef * dampcoef) *
                         (t(i, j, k) - t00[k]);
    }
  }
};

inline void rk_rayleigh_damp(const DiffusionParams& p,
                             View3<const real> u,
                             View3<const real> v,
                             View3<const real> w,
                             View3<const real> t,
                             View3<const real> ph,
                             BView3<const real> phb,
                             View2<const real> mut,
                             View2<const real> muu,
                             View2<const real> muv,
                             View1<const real> c1h,
                             View1<const real> c2h,
                             View1<const real> c1f,
                             View1<const real> c2f,
                             View1<const real> u_base,
                             View1<const real> v_base,
                             View1<const real> t_base,
                             View1<const real> z_base,
                             View3<real> ru_tendf,
                             View3<real> rv_tendf,
                             View3<real> rw_tendf,
                             View3<real> t_tendf,
                             int nx, int ny, int nz) {
  real pii = real(2.0) * std::asin(real(1.0));

  column_for(Range2D{0, nx + 1, 0, ny},
             KRayleighDampU{u, ph, phb, muu, c1h, c2h, u_base, z_base,
                            ru_tendf, p.dampcoef, p.zdamp, pii, nx, nz, p.periodic_x});

  column_for(Range2D{0, nx, 0, ny + 1},
             KRayleighDampV{v, ph, phb, muv, c1h, c2h, v_base, z_base,
                            rv_tendf, p.dampcoef, p.zdamp, pii, ny, nz, p.periodic_y});

  column_for(Range2D{0, nx, 0, ny},
             KRayleighDampWT{w, t, ph, phb, mut, c1h, c2h, c1f, c2f, t_base, z_base,
                             rw_tendf, t_tendf, p.dampcoef, p.zdamp, pii, nz});
}

}  // namespace cppwrf::dyn::diffusion
