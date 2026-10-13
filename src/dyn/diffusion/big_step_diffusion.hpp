#pragma once

#include "core/backend.hpp"
#include "core/constants.hpp"
#include "core/field.hpp"
#include "dyn/diffusion/diffusion_params.hpp"

namespace cppwrf::dyn::diffusion {

// Ported from WRF v4.6.1 dyn_em/module_big_step_utilities_em.F:
//   horizontal_diffusion      (lines 2715-2950)
//   horizontal_diffusion_3dmp (lines 2954-3060)
//   vertical_diffusion        (lines 3064-3187)
//   vertical_diffusion_mp     (lines 3192-3281)
//   vertical_diffusion_3dmp   (lines 3286-3376)
//   vertical_diffusion_u      (lines 3382-3480)
//   vertical_diffusion_v      (lines 3485-3581)

constexpr real Gr_bs = real(cst::g);

// =============================================================================
// Horizontal diffusion (diff_opt = 1)
// =============================================================================

struct KHorizDiffU1 {
  View3<const real> field, xkmhd;
  View2<const real> mut;
  BView2<const real> msftx, msfty, msfux, msfuy;
  View1<const real> c1, c2;
  View3<real> tendency;
  real rdx, rdy;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    real mkrdxm = (msftx(i - 1, j) / msfty(i - 1, j)) * (c1(k) * mut(i - 1, j) + c2(k)) *
                  xkmhd(i - 1, j, k) * rdx;
    real mkrdxp = (msftx(i, j) / msfty(i, j)) * (c1(k) * mut(i, j) + c2(k)) *
                  xkmhd(i, j, k) * rdx;
    real mrdx = msfux(i, j) * msfuy(i, j) * rdx;

    real mkrdym = ((msfuy(i, j) + msfuy(i, j - 1)) / (msfux(i, j) + msfux(i, j - 1))) *
                  real(0.25) * ((c1(k) * mut(i, j) + c2(k)) + (c1(k) * mut(i, j - 1) + c2(k)) +
                                (c1(k) * mut(i - 1, j - 1) + c2(k)) + (c1(k) * mut(i - 1, j) + c2(k))) *
                  real(0.25) * (xkmhd(i, j, k) + xkmhd(i, j - 1, k) +
                                xkmhd(i - 1, j - 1, k) + xkmhd(i - 1, j, k)) * rdy;

    real mkrdyp = ((msfuy(i, j) + msfuy(i, j + 1)) / (msfux(i, j) + msfux(i, j + 1))) *
                  real(0.25) * ((c1(k) * mut(i, j) + c2(k)) + (c1(k) * mut(i, j + 1) + c2(k)) +
                                (c1(k) * mut(i - 1, j + 1) + c2(k)) + (c1(k) * mut(i - 1, j) + c2(k))) *
                  real(0.25) * (xkmhd(i, j, k) + xkmhd(i, j + 1, k) +
                                xkmhd(i - 1, j + 1, k) + xkmhd(i - 1, j, k)) * rdy;
    real mrdy = msfux(i, j) * msfuy(i, j) * rdy;

    tendency(i, j, k) += mrdx * (mkrdxp * (field(i + 1, j, k) - field(i, j, k)) -
                                 mkrdxm * (field(i, j, k) - field(i - 1, j, k))) +
                         mrdy * (mkrdyp * (field(i, j + 1, k) - field(i, j, k)) -
                                 mkrdym * (field(i, j, k) - field(i, j - 1, k)));
  }
};

struct KHorizDiffV1 {
  View3<const real> field, xkmhd;
  View2<const real> mut;
  BView2<const real> msftx, msfty, msfvx, msfvy;
  View1<const real> c1, c2;
  View3<real> tendency;
  real rdx, rdy;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    real mkrdxm = ((msfvx(i, j) + msfvx(i - 1, j)) / (msfvy(i, j) + msfvy(i - 1, j))) *
                  real(0.25) * ((c1(k) * mut(i, j) + c2(k)) + (c1(k) * mut(i, j - 1) + c2(k)) +
                                (c1(k) * mut(i - 1, j - 1) + c2(k)) + (c1(k) * mut(i - 1, j) + c2(k))) *
                  real(0.25) * (xkmhd(i, j, k) + xkmhd(i, j - 1, k) +
                                xkmhd(i - 1, j - 1, k) + xkmhd(i - 1, j, k)) * rdx;

    real mkrdxp = ((msfvx(i, j) + msfvx(i + 1, j)) / (msfvy(i, j) + msfvy(i + 1, j))) *
                  real(0.25) * ((c1(k) * mut(i, j) + c2(k)) + (c1(k) * mut(i, j - 1) + c2(k)) +
                                (c1(k) * mut(i + 1, j - 1) + c2(k)) + (c1(k) * mut(i + 1, j) + c2(k))) *
                  real(0.25) * (xkmhd(i, j, k) + xkmhd(i, j - 1, k) +
                                xkmhd(i + 1, j - 1, k) + xkmhd(i + 1, j, k)) * rdx;
    real mrdx = msfvx(i, j) * msfvy(i, j) * rdx;

    // Follows WRF v4.6.1 line 2854-2856 exactly (no c1*mut+c2 in y direction for v)
    real mkrdym = (msfty(i, j - 1) / msftx(i, j - 1)) * xkmhd(i, j - 1, k) * rdy;
    real mkrdyp = (msfty(i, j) / msftx(i, j)) * xkmhd(i, j, k) * rdy;
    real mrdy = msfvx(i, j) * msfvy(i, j) * rdy;

    tendency(i, j, k) += mrdx * (mkrdxp * (field(i + 1, j, k) - field(i, j, k)) -
                                 mkrdxm * (field(i, j, k) - field(i - 1, j, k))) +
                         mrdy * (mkrdyp * (field(i, j + 1, k) - field(i, j, k)) -
                                 mkrdym * (field(i, j, k) - field(i, j - 1, k)));
  }
};

struct KHorizDiffW1 {
  View3<const real> field, xkmhd;
  View2<const real> mut;
  BView2<const real> msftx, msfty, msfux, msfuy, msfvy, msfvx_inv;
  View1<const real> c1, c2;
  View3<real> tendency;
  real rdx, rdy;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    real mkrdxm = (msfux(i, j) / msfuy(i, j)) *
                  real(0.25) * ((c1(k) * mut(i, j) + c2(k)) + (c1(k) * mut(i - 1, j) + c2(k)) +
                                (c1(k) * mut(i, j) + c2(k)) + (c1(k) * mut(i - 1, j) + c2(k))) *
                  real(0.25) * (xkmhd(i, j, k) + xkmhd(i - 1, j, k) +
                                xkmhd(i, j, k - 1) + xkmhd(i - 1, j, k - 1)) * rdx;

    real mkrdxp = (msfux(i + 1, j) / msfuy(i + 1, j)) *
                  real(0.25) * ((c1(k) * mut(i + 1, j) + c2(k)) + (c1(k) * mut(i, j) + c2(k)) +
                                (c1(k) * mut(i + 1, j) + c2(k)) + (c1(k) * mut(i, j) + c2(k))) *
                  real(0.25) * (xkmhd(i + 1, j, k) + xkmhd(i, j, k) +
                                xkmhd(i + 1, j, k - 1) + xkmhd(i, j, k - 1)) * rdx;
    real mrdx = msftx(i, j) * msfty(i, j) * rdx;

    real mkrdym = (msfvy(i, j) * msfvx_inv(i, j)) *
                  real(0.25) * ((c1(k) * mut(i, j) + c2(k)) + (c1(k) * mut(i, j - 1) + c2(k)) +
                                (c1(k) * mut(i, j) + c2(k)) + (c1(k) * mut(i, j - 1) + c2(k))) *
                  real(0.25) * (xkmhd(i, j, k) + xkmhd(i, j - 1, k) +
                                xkmhd(i, j, k - 1) + xkmhd(i, j - 1, k - 1)) * rdy;

    real mkrdyp = (msfvy(i, j + 1) * msfvx_inv(i, j + 1)) *
                  real(0.25) * ((c1(k) * mut(i, j + 1) + c2(k)) + (c1(k) * mut(i, j) + c2(k)) +
                                (c1(k) * mut(i, j + 1) + c2(k)) + (c1(k) * mut(i, j) + c2(k))) *
                  real(0.25) * (xkmhd(i, j + 1, k) + xkmhd(i, j, k) +
                                xkmhd(i, j + 1, k - 1) + xkmhd(i, j, k - 1)) * rdy;
    real mrdy = msftx(i, j) * msfty(i, j) * rdy;

    tendency(i, j, k) += mrdx * (mkrdxp * (field(i + 1, j, k) - field(i, j, k)) -
                                 mkrdxm * (field(i, j, k) - field(i - 1, j, k))) +
                         mrdy * (mkrdyp * (field(i, j + 1, k) - field(i, j, k)) -
                                 mkrdym * (field(i, j, k) - field(i, j - 1, k)));
  }
};

struct KHorizDiffM1 {
  View3<const real> field, xkmhd;
  View2<const real> mut;
  BView2<const real> msftx, msfty, msfux, msfuy, msfvy, msfvx_inv;
  View1<const real> c1, c2;
  View3<real> tendency;
  real rdx, rdy;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    real mkrdxm = (msfux(i, j) / msfuy(i, j)) * real(0.5) * (xkmhd(i, j, k) + xkmhd(i - 1, j, k)) *
                  real(0.5) * ((c1(k) * mut(i, j) + c2(k)) + (c1(k) * mut(i - 1, j) + c2(k))) * rdx;
    real mkrdxp = (msfux(i + 1, j) / msfuy(i + 1, j)) * real(0.5) * (xkmhd(i + 1, j, k) + xkmhd(i, j, k)) *
                  real(0.5) * ((c1(k) * mut(i + 1, j) + c2(k)) + (c1(k) * mut(i, j) + c2(k))) * rdx;
    real mrdx = msftx(i, j) * msfty(i, j) * rdx;

    real mkrdym = (msfvy(i, j) * msfvx_inv(i, j)) * real(0.5) * (xkmhd(i, j, k) + xkmhd(i, j - 1, k)) *
                  real(0.5) * ((c1(k) * mut(i, j) + c2(k)) + (c1(k) * mut(i, j - 1) + c2(k))) * rdy;
    real mkrdyp = (msfvy(i, j + 1) * msfvx_inv(i, j + 1)) * real(0.5) * (xkmhd(i, j + 1, k) + xkmhd(i, j, k)) *
                  real(0.5) * ((c1(k) * mut(i, j + 1) + c2(k)) + (c1(k) * mut(i, j) + c2(k))) * rdy;
    real mrdy = msftx(i, j) * msfty(i, j) * rdy;

    tendency(i, j, k) += mrdx * (mkrdxp * (field(i + 1, j, k) - field(i, j, k)) -
                                 mkrdxm * (field(i, j, k) - field(i - 1, j, k))) +
                         mrdy * (mkrdyp * (field(i, j + 1, k) - field(i, j, k)) -
                                 mkrdym * (field(i, j, k) - field(i, j - 1, k)));
  }
};

struct KHorizDiff3DMP {
  View3<const real> field;
  BView3<const real> base_3d;  // base states may be horizontally uniform (BaseField3)
  View3<const real> xkmhd;
  View2<const real> mut;
  BView2<const real> msftx, msfty, msfux, msfuy, msfvy, msfvx_inv;
  View1<const real> c1, c2;
  View3<real> tendency;
  real rdx, rdy;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    real mkrdxm = (msfux(i, j) / msfuy(i, j)) * real(0.5) * (xkmhd(i, j, k) + xkmhd(i - 1, j, k)) *
                  real(0.5) * ((c1(k) * mut(i, j) + c2(k)) + (c1(k) * mut(i - 1, j) + c2(k))) * rdx;
    real mkrdxp = (msfux(i + 1, j) / msfuy(i + 1, j)) * real(0.5) * (xkmhd(i + 1, j, k) + xkmhd(i, j, k)) *
                  real(0.5) * ((c1(k) * mut(i + 1, j) + c2(k)) + (c1(k) * mut(i, j) + c2(k))) * rdx;
    real mrdx = msftx(i, j) * msfty(i, j) * rdx;

    real mkrdym = (msfvy(i, j) * msfvx_inv(i, j)) * real(0.5) * (xkmhd(i, j, k) + xkmhd(i, j - 1, k)) *
                  real(0.5) * ((c1(k) * mut(i, j) + c2(k)) + (c1(k) * mut(i, j - 1) + c2(k))) * rdy;
    real mkrdyp = (msfvy(i, j + 1) * msfvx_inv(i, j + 1)) * real(0.5) * (xkmhd(i, j + 1, k) + xkmhd(i, j, k)) *
                  real(0.5) * ((c1(k) * mut(i, j + 1) + c2(k)) + (c1(k) * mut(i, j) + c2(k))) * rdy;
    real mrdy = msftx(i, j) * msfty(i, j) * rdy;

    real dfield_xp = (field(i + 1, j, k) - field(i, j, k)) - (base_3d(i + 1, j, k) - base_3d(i, j, k));
    real dfield_xm = (field(i, j, k) - field(i - 1, j, k)) - (base_3d(i, j, k) - base_3d(i - 1, j, k));
    real dfield_yp = (field(i, j + 1, k) - field(i, j, k)) - (base_3d(i, j + 1, k) - base_3d(i, j, k));
    real dfield_ym = (field(i, j, k) - field(i, j - 1, k)) - (base_3d(i, j, k) - base_3d(i, j - 1, k));

    tendency(i, j, k) += mrdx * (mkrdxp * dfield_xp - mkrdxm * dfield_xm) +
                         mrdy * (mkrdyp * dfield_yp - mkrdym * dfield_ym);
  }
};

// =============================================================================
// Vertical diffusion (diff_opt = 1)
// =============================================================================

struct KVertDiffW1 {
  View3<const real> field, alt;
  View2<const real> mut;
  View1<const real> c1f, c2f, rdn, rdnw;
  View3<real> tendency;
  real kvdif;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    real vflux[128];
    for (int k = 0; k < nz - 1; ++k) {
      vflux[k] = (kvdif / alt(i, j, k)) * rdnw(k) * (field(i, j, k + 1) - field(i, j, k));
    }
    vflux[nz - 1] = real(0);

    for (int k = 1; k < nz; ++k) {
      real mu = c1f(k) * mut(i, j) + c2f(k);
      real alt_avg = real(0.5) * (alt(i, j, k) + alt(i, j, k - 1));
      tendency(i, j, k) += rdn(k) * Gr_bs * Gr_bs / mu / alt_avg * (vflux[k] - vflux[k - 1]);
    }
  }
};

struct KVertDiffM1 {
  View3<const real> field, alt;
  View2<const real> mut;
  View1<const real> c1h, c2h, rdn, rdnw;
  View3<real> tendency;
  real kvdif;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    real vflux[128];
    for (int k = 0; k < nz - 1; ++k) {
      vflux[k] = kvdif * rdn(k + 1) / (real(0.5) * (alt(i, j, k) + alt(i, j, k + 1))) *
                 (field(i, j, k + 1) - field(i, j, k));
    }
    real vf0 = vflux[0];
    vflux[nz - 1] = real(0);

    tendency(i, j, 0) += Gr_bs * Gr_bs / (c1h(0) * mut(i, j) + c2h(0)) / alt(i, j, 0) *
                         rdnw(0) * (vflux[0] - vf0);  // vflux(0) - vflux(-1) = 0
    for (int k = 1; k < nz; ++k) {
      real mu = c1h(k) * mut(i, j) + c2h(k);
      tendency(i, j, k) += Gr_bs * Gr_bs / mu / alt(i, j, k) * rdnw(k) * (vflux[k] - vflux[k - 1]);
    }
  }
};

struct KVertDiffMP1 {
  View3<const real> field, alt;
  BView3<const real> base;
  View2<const real> mut;
  View1<const real> c1h, c2h, rdn, rdnw;
  View3<real> tendency;
  real kvdif;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    real vflux[128];
    for (int k = 0; k < nz - 1; ++k) {
      real dpert = (field(i, j, k + 1) - field(i, j, k)) - (base(i, j, k + 1) - base(i, j, k));
      vflux[k] = kvdif * rdn(k + 1) / (real(0.5) * (alt(i, j, k) + alt(i, j, k + 1))) * dpert;
    }
    real vf0 = vflux[0];
    vflux[nz - 1] = real(0);

    tendency(i, j, 0) += Gr_bs * Gr_bs / (c1h(0) * mut(i, j) + c2h(0)) / alt(i, j, 0) *
                         rdnw(0) * (vflux[0] - vf0);
    for (int k = 1; k < nz; ++k) {
      real mu = c1h(k) * mut(i, j) + c2h(k);
      tendency(i, j, k) += Gr_bs * Gr_bs / mu / alt(i, j, k) * rdnw(k) * (vflux[k] - vflux[k - 1]);
    }
  }
};

struct KVertDiffU1 {
  View3<const real> field, alt;
  View1<const real> u_base;
  View2<const real> muu;
  View1<const real> c1h, c2h, rdn, rdnw;
  View3<real> tendency;
  real kvdif;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    real vflux[128];
    for (int k = 0; k < nz - 1; ++k) {
      real a = real(0.25) * (alt(i, j, k) + alt(i - 1, j, k) + alt(i, j, k + 1) + alt(i - 1, j, k + 1));
      real dpert = (field(i, j, k + 1) - field(i, j, k)) - (u_base(k + 1) - u_base(k));
      vflux[k] = kvdif * rdn(k + 1) / a * dpert;
    }
    real vf0 = vflux[0];
    vflux[nz - 1] = real(0);

    // WRF vertical_diffusion_u: loops k=0..ktf-1 (top level ktf-1 in 0-based gets no vertical diffusion)
    tendency(i, j, 0) += Gr_bs * Gr_bs * rdnw(0) / (c1h(0) * muu(i, j) + c2h(0)) /
                         (real(0.5) * (alt(i - 1, j, 0) + alt(i, j, 0))) * (vflux[0] - vf0);
    for (int k = 1; k < nz - 1; ++k) {
      real mu = c1h(k) * muu(i, j) + c2h(k);
      real a = real(0.5) * (alt(i - 1, j, k) + alt(i, j, k));
      tendency(i, j, k) += Gr_bs * Gr_bs * rdnw(k) / mu / a * (vflux[k] - vflux[k - 1]);
    }
  }
};

struct KVertDiffV1 {
  View3<const real> field, alt;
  View1<const real> v_base;
  View2<const real> muv;
  View1<const real> c1h, c2h, rdn, rdnw;
  View3<real> tendency;
  real kvdif;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    real vflux[128];
    for (int k = 0; k < nz - 1; ++k) {
      real a = real(0.25) * (alt(i, j, k) + alt(i, j - 1, k) + alt(i, j, k + 1) + alt(i, j - 1, k + 1));
      real dpert = (field(i, j, k + 1) - field(i, j, k)) - (v_base(k + 1) - v_base(k));
      vflux[k] = kvdif * rdn(k + 1) / a * dpert;
    }
    real vf0 = vflux[0];
    vflux[nz - 1] = real(0);

    tendency(i, j, 0) += Gr_bs * Gr_bs * rdnw(0) / (c1h(0) * muv(i, j) + c2h(0)) /
                         (real(0.5) * (alt(i, j - 1, 0) + alt(i, j, 0))) * (vflux[0] - vf0);
    for (int k = 1; k < nz - 1; ++k) {
      real mu = c1h(k) * muv(i, j) + c2h(k);
      real a = real(0.5) * (alt(i, j - 1, k) + alt(i, j, k));
      tendency(i, j, k) += Gr_bs * Gr_bs * rdnw(k) / mu / a * (vflux[k] - vflux[k - 1]);
    }
  }
};

// =============================================================================
// Launchers
// =============================================================================

inline void horizontal_diffusion(char name,
                                 const DiffusionParams& p,
                                 View3<const real> field,
                                 View3<real> tendency,
                                 View2<const real> mut,
                                 BView2<const real> msftx,
                                 BView2<const real> msfty,
                                 BView2<const real> msfux,
                                 BView2<const real> msfuy,
                                 BView2<const real> msfvx,
                                 BView2<const real> msfvy,
                                 BView2<const real> msfvx_inv,
                                 View1<const real> c1,
                                 View1<const real> c2,
                                 View3<const real> xkmhd,
                                 real rdx, real rdy,
                                 int nx, int ny, int nz) {
  if (name == 'u') {
    int is = (p.open_xs || p.specified) ? 1 : 0;
    int ie = (p.open_xe || p.specified) ? nx : nx + 1;
    int js = (p.open_ys || p.specified) ? 1 : 0;
    int je = (p.open_ye || p.specified) ? ny - 1 : ny;
    if (p.periodic_x) { is = 0; ie = nx + 1; }

    parallel_for(Range3D{is, ie, js, je, 0, nz},
                 KHorizDiffU1{field, xkmhd, mut, msftx, msfty, msfux, msfuy, c1, c2, tendency, rdx, rdy});
  } else if (name == 'v') {
    int is = (p.open_xs || p.specified) ? 1 : 0;
    int ie = (p.open_xe || p.specified) ? nx - 1 : nx;
    int js = (p.open_ys || p.specified) ? 1 : 0;
    int je = (p.open_ye || p.specified) ? ny : ny + 1;
    if (p.periodic_x) { is = 0; ie = nx; }

    parallel_for(Range3D{is, ie, js, je, 0, nz},
                 KHorizDiffV1{field, xkmhd, mut, msftx, msfty, msfvx, msfvy, c1, c2, tendency, rdx, rdy});
  } else if (name == 'w') {
    int is = (p.open_xs || p.specified) ? 1 : 0;
    int ie = (p.open_xe || p.specified) ? nx - 1 : nx;
    int js = (p.open_ys || p.specified) ? 1 : 0;
    int je = (p.open_ye || p.specified) ? ny - 1 : ny;
    if (p.periodic_x) { is = 0; ie = nx; }

    parallel_for(Range3D{is, ie, js, je, 1, nz},
                 KHorizDiffW1{field, xkmhd, mut, msftx, msfty, msfux, msfuy, msfvy, msfvx_inv, c1, c2, tendency, rdx, rdy});
  } else {  // 'm' (scalar / theta / moisture)
    int is = (p.open_xs || p.specified) ? 1 : 0;
    int ie = (p.open_xe || p.specified) ? nx - 1 : nx;
    int js = (p.open_ys || p.specified) ? 1 : 0;
    int je = (p.open_ye || p.specified) ? ny - 1 : ny;
    if (p.periodic_x) { is = 0; ie = nx; }

    parallel_for(Range3D{is, ie, js, je, 0, nz},
                 KHorizDiffM1{field, xkmhd, mut, msftx, msfty, msfux, msfuy, msfvy, msfvx_inv, c1, c2, tendency, rdx, rdy});
  }
}

inline void horizontal_diffusion_3dmp(const DiffusionParams& p,
                                      View3<const real> field,
                                      BView3<const real> base_3d,
                                      View3<real> tendency,
                                      View2<const real> mut,
                                      BView2<const real> msftx,
                                      BView2<const real> msfty,
                                      BView2<const real> msfux,
                                      BView2<const real> msfuy,
                                      BView2<const real> msfvy,
                                      BView2<const real> msfvx_inv,
                                      View1<const real> c1,
                                      View1<const real> c2,
                                      View3<const real> xkmhd,
                                      real rdx, real rdy,
                                      int nx, int ny, int nz) {
  int is = (p.open_xs || p.specified) ? 1 : 0;
  int ie = (p.open_xe || p.specified) ? nx - 1 : nx;
  int js = (p.open_ys || p.specified) ? 1 : 0;
  int je = (p.open_ye || p.specified) ? ny - 1 : ny;
  if (p.periodic_x) { is = 0; ie = nx; }

  parallel_for(Range3D{is, ie, js, je, 0, nz},
               KHorizDiff3DMP{field, base_3d, xkmhd, mut, msftx, msfty, msfux, msfuy, msfvy, msfvx_inv, c1, c2, tendency, rdx, rdy});
}

inline void vertical_diffusion(char name,
                               View3<const real> field,
                               View3<real> tendency,
                               View2<const real> mut,
                               View1<const real> c1,
                               View1<const real> c2,
                               View3<const real> alt,
                               View1<const real> rdn,
                               View1<const real> rdnw,
                               real kvdif,
                               int nx, int ny, int nz) {
  if (name == 'w') {
    column_for(Range2D{0, nx, 0, ny},
               KVertDiffW1{field, alt, mut, c1, c2, rdn, rdnw, tendency, kvdif, nz});
  } else {  // 'm'
    column_for(Range2D{0, nx, 0, ny},
               KVertDiffM1{field, alt, mut, c1, c2, rdn, rdnw, tendency, kvdif, nz});
  }
}

inline void vertical_diffusion_mp(View3<const real> field,
                                  View3<real> tendency,
                                  BView3<const real> base,
                                  View2<const real> mut,
                                  View1<const real> c1,
                                  View1<const real> c2,
                                  View3<const real> alt,
                                  View1<const real> rdn,
                                  View1<const real> rdnw,
                                  real kvdif,
                                  int nx, int ny, int nz) {
  column_for(Range2D{0, nx, 0, ny},
             KVertDiffMP1{field, alt, base, mut, c1, c2, rdn, rdnw, tendency, kvdif, nz});
}

inline void vertical_diffusion_u(const DiffusionParams& p,
                                 View3<const real> field,
                                 View3<real> tendency,
                                 View1<const real> u_base,
                                 View2<const real> muu,
                                 View1<const real> c1h,
                                 View1<const real> c2h,
                                 View3<const real> alt,
                                 View1<const real> rdn,
                                 View1<const real> rdnw,
                                 real kvdif,
                                 int nx, int ny, int nz) {
  int is = (p.open_xs || p.specified) ? 1 : 0;
  int ie = (p.open_xe || p.specified) ? nx : nx + 1;
  int js = 0;
  int je = ny;
  if (p.periodic_x) { is = 0; ie = nx + 1; }

  column_for(Range2D{is, ie, js, je},
             KVertDiffU1{field, alt, u_base, muu, c1h, c2h, rdn, rdnw, tendency, kvdif, nz});
}

inline void vertical_diffusion_v(const DiffusionParams& p,
                                 View3<const real> field,
                                 View3<real> tendency,
                                 View1<const real> v_base,
                                 View2<const real> muv,
                                 View1<const real> c1h,
                                 View1<const real> c2h,
                                 View3<const real> alt,
                                 View1<const real> rdn,
                                 View1<const real> rdnw,
                                 real kvdif,
                                 int nx, int ny, int nz) {
  int is = 0;
  int ie = nx;
  int js = (p.open_ys || p.specified) ? 1 : 0;
  int je = (p.open_ye || p.specified) ? ny : ny + 1;

  column_for(Range2D{is, ie, js, je},
             KVertDiffV1{field, alt, v_base, muv, c1h, c2h, rdn, rdnw, tendency, kvdif, nz});
}

}  // namespace cppwrf::dyn::diffusion
