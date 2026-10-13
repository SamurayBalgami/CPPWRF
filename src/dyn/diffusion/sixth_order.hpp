#pragma once

#include "core/backend.hpp"
#include "core/field.hpp"
#include "dyn/diffusion/diffusion_params.hpp"

namespace cppwrf::dyn::diffusion {

// Ported from WRF v4.6.1 dyn_em/module_big_step_utilities_em.F: sixth_order_diffusion (lines 6220-6636)
// Supports diff_6th_opt = 1 (linear) and diff_6th_opt = 2 (monotonic flux-limited), diff_6th_slopeopt = 0.

struct KSixthOrderU {
  View3<const real> field;
  View2<const real> mut;
  BView2<const real> msfux, msfuy;
  View1<const real> c1, c2;
  View3<real> tendency;
  real diff_6th_coef;
  int diff_6th_opt;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    // x-flux
    real dflux_x_p0 = real(10.0) * (field(i, j, k) - field(i - 1, j, k)) -
                      real(5.0) * (field(i + 1, j, k) - field(i - 2, j, k)) +
                      (field(i + 2, j, k) - field(i - 3, j, k));
    real dflux_x_p1 = real(10.0) * (field(i + 1, j, k) - field(i, j, k)) -
                      real(5.0) * (field(i + 2, j, k) - field(i - 1, j, k)) +
                      (field(i + 3, j, k) - field(i - 2, j, k));

    if (diff_6th_opt == 2) {
      if (dflux_x_p0 * (field(i, j, k) - field(i - 1, j, k)) <= real(0)) dflux_x_p0 = real(0);
      if (dflux_x_p1 * (field(i + 1, j, k) - field(i, j, k)) <= real(0)) dflux_x_p1 = real(0);
    }
    real mu_x_p0 = c1(k) * mut(i - 1, j) + c2(k);
    real mu_x_p1 = c1(k) * mut(i, j) + c2(k);
    real tend_x = diff_6th_coef * msfux(i, j) * (mu_x_p1 * dflux_x_p1 - mu_x_p0 * dflux_x_p0);

    // y-flux
    real dflux_y_p0 = real(10.0) * (field(i, j, k) - field(i, j - 1, k)) -
                      real(5.0) * (field(i, j + 1, k) - field(i, j - 2, k)) +
                      (field(i, j + 2, k) - field(i, j - 3, k));
    real dflux_y_p1 = real(10.0) * (field(i, j + 1, k) - field(i, j, k)) -
                      real(5.0) * (field(i, j + 2, k) - field(i, j - 1, k)) +
                      (field(i, j + 3, k) - field(i, j - 2, k));

    if (diff_6th_opt == 2) {
      if (dflux_y_p0 * (field(i, j, k) - field(i, j - 1, k)) <= real(0)) dflux_y_p0 = real(0);
      if (dflux_y_p1 * (field(i, j + 1, k) - field(i, j, k)) <= real(0)) dflux_y_p1 = real(0);
    }
    real mu_y_p0 = real(0.25) * ((c1(k) * mut(i - 1, j - 1) + c2(k)) + (c1(k) * mut(i, j - 1) + c2(k)) +
                                 (c1(k) * mut(i - 1, j) + c2(k)) + (c1(k) * mut(i, j) + c2(k)));
    real mu_y_p1 = real(0.25) * ((c1(k) * mut(i - 1, j) + c2(k)) + (c1(k) * mut(i, j) + c2(k)) +
                                 (c1(k) * mut(i - 1, j + 1) + c2(k)) + (c1(k) * mut(i, j + 1) + c2(k)));
    real tend_y = diff_6th_coef * msfuy(i, j) * (mu_y_p1 * dflux_y_p1 - mu_y_p0 * dflux_y_p0);

    tendency(i, j, k) += tend_x + tend_y;
  }
};

struct KSixthOrderV {
  View3<const real> field;
  View2<const real> mut;
  BView2<const real> msfvx, msfvy;
  View1<const real> c1, c2;
  View3<real> tendency;
  real diff_6th_coef;
  int diff_6th_opt;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    // x-flux
    real dflux_x_p0 = real(10.0) * (field(i, j, k) - field(i - 1, j, k)) -
                      real(5.0) * (field(i + 1, j, k) - field(i - 2, j, k)) +
                      (field(i + 2, j, k) - field(i - 3, j, k));
    real dflux_x_p1 = real(10.0) * (field(i + 1, j, k) - field(i, j, k)) -
                      real(5.0) * (field(i + 2, j, k) - field(i - 1, j, k)) +
                      (field(i + 3, j, k) - field(i - 2, j, k));

    if (diff_6th_opt == 2) {
      if (dflux_x_p0 * (field(i, j, k) - field(i - 1, j, k)) <= real(0)) dflux_x_p0 = real(0);
      if (dflux_x_p1 * (field(i + 1, j, k) - field(i, j, k)) <= real(0)) dflux_x_p1 = real(0);
    }
    real mu_x_p0 = real(0.25) * ((c1(k) * mut(i - 1, j - 1) + c2(k)) + (c1(k) * mut(i, j - 1) + c2(k)) +
                                 (c1(k) * mut(i - 1, j) + c2(k)) + (c1(k) * mut(i, j) + c2(k)));
    real mu_x_p1 = real(0.25) * ((c1(k) * mut(i, j - 1) + c2(k)) + (c1(k) * mut(i + 1, j - 1) + c2(k)) +
                                 (c1(k) * mut(i, j) + c2(k)) + (c1(k) * mut(i + 1, j) + c2(k)));
    real tend_x = diff_6th_coef * msfvx(i, j) * (mu_x_p1 * dflux_x_p1 - mu_x_p0 * dflux_x_p0);

    // y-flux
    real dflux_y_p0 = real(10.0) * (field(i, j, k) - field(i, j - 1, k)) -
                      real(5.0) * (field(i, j + 1, k) - field(i, j - 2, k)) +
                      (field(i, j + 2, k) - field(i, j - 3, k));
    real dflux_y_p1 = real(10.0) * (field(i, j + 1, k) - field(i, j, k)) -
                      real(5.0) * (field(i, j + 2, k) - field(i, j - 1, k)) +
                      (field(i, j + 3, k) - field(i, j - 2, k));

    if (diff_6th_opt == 2) {
      if (dflux_y_p0 * (field(i, j, k) - field(i, j - 1, k)) <= real(0)) dflux_y_p0 = real(0);
      if (dflux_y_p1 * (field(i, j + 1, k) - field(i, j, k)) <= real(0)) dflux_y_p1 = real(0);
    }
    real mu_y_p0 = c1(k) * mut(i, j - 1) + c2(k);
    real mu_y_p1 = c1(k) * mut(i, j) + c2(k);
    real tend_y = diff_6th_coef * msfvy(i, j) * (mu_y_p1 * dflux_y_p1 - mu_y_p0 * dflux_y_p0);

    tendency(i, j, k) += tend_x + tend_y;
  }
};

struct KSixthOrderM {
  View3<const real> field;
  View2<const real> mut;
  BView2<const real> msftx, msfty;
  View1<const real> c1, c2;
  View3<real> tendency;
  real diff_6th_coef;
  int diff_6th_opt;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    // x-flux
    real dflux_x_p0 = real(10.0) * (field(i, j, k) - field(i - 1, j, k)) -
                      real(5.0) * (field(i + 1, j, k) - field(i - 2, j, k)) +
                      (field(i + 2, j, k) - field(i - 3, j, k));
    real dflux_x_p1 = real(10.0) * (field(i + 1, j, k) - field(i, j, k)) -
                      real(5.0) * (field(i + 2, j, k) - field(i - 1, j, k)) +
                      (field(i + 3, j, k) - field(i - 2, j, k));

    if (diff_6th_opt == 2) {
      if (dflux_x_p0 * (field(i, j, k) - field(i - 1, j, k)) <= real(0)) dflux_x_p0 = real(0);
      if (dflux_x_p1 * (field(i + 1, j, k) - field(i, j, k)) <= real(0)) dflux_x_p1 = real(0);
    }
    real mu_x_p0 = real(0.5) * ((c1(k) * mut(i - 1, j) + c2(k)) + (c1(k) * mut(i, j) + c2(k)));
    real mu_x_p1 = real(0.5) * ((c1(k) * mut(i, j) + c2(k)) + (c1(k) * mut(i + 1, j) + c2(k)));
    real tend_x = diff_6th_coef * msftx(i, j) * (mu_x_p1 * dflux_x_p1 - mu_x_p0 * dflux_x_p0);

    // y-flux
    real dflux_y_p0 = real(10.0) * (field(i, j, k) - field(i, j - 1, k)) -
                      real(5.0) * (field(i, j + 1, k) - field(i, j - 2, k)) +
                      (field(i, j + 2, k) - field(i, j - 3, k));
    real dflux_y_p1 = real(10.0) * (field(i, j + 1, k) - field(i, j, k)) -
                      real(5.0) * (field(i, j + 2, k) - field(i, j - 1, k)) +
                      (field(i, j + 3, k) - field(i, j - 2, k));

    if (diff_6th_opt == 2) {
      if (dflux_y_p0 * (field(i, j, k) - field(i, j - 1, k)) <= real(0)) dflux_y_p0 = real(0);
      if (dflux_y_p1 * (field(i, j + 1, k) - field(i, j, k)) <= real(0)) dflux_y_p1 = real(0);
    }
    real mu_y_p0 = real(0.5) * ((c1(k) * mut(i, j - 1) + c2(k)) + (c1(k) * mut(i, j) + c2(k)));
    real mu_y_p1 = real(0.5) * ((c1(k) * mut(i, j) + c2(k)) + (c1(k) * mut(i, j + 1) + c2(k)));
    real tend_y = diff_6th_coef * msfty(i, j) * (mu_y_p1 * dflux_y_p1 - mu_y_p0 * dflux_y_p0);

    tendency(i, j, k) += tend_x + tend_y;
  }
};

inline void sixth_order_diffusion(char name,
                                  const DiffusionParams& p,
                                  View3<const real> field,
                                  View3<real> tendency,
                                  View2<const real> mut,
                                  real dt,
                                  View1<const real> c1,
                                  View1<const real> c2,
                                  BView2<const real> msftx,
                                  BView2<const real> msfty,
                                  BView2<const real> msfux,
                                  BView2<const real> msfuy,
                                  BView2<const real> msfvx,
                                  BView2<const real> msfvy,
                                  int nx, int ny, int nz) {
  if (p.diff_6th_opt == 0 || p.diff_6th_factor == 0) return;

  real diff_6th_coef = p.diff_6th_factor * real(0.015625) / (real(2.0) * dt);

  if (name == 'u') {
    int is = p.specified ? 3 : 0;
    int ie = p.specified ? nx - 2 : nx + 1;
    int js = p.specified ? 3 : 0;
    int je = p.specified ? ny - 3 : ny;
    if (p.periodic_x) { is = 0; ie = nx + 1; }

    parallel_for(Range3D{is, ie, js, je, 0, nz},
                 KSixthOrderU{field, mut, msfux, msfuy, c1, c2, tendency, diff_6th_coef, p.diff_6th_opt});
  } else if (name == 'v') {
    int is = p.specified ? 3 : 0;
    int ie = p.specified ? nx - 3 : nx;
    int js = p.specified ? 3 : 0;
    int je = p.specified ? ny - 2 : ny + 1;
    if (p.periodic_x) { is = 0; ie = nx; }

    parallel_for(Range3D{is, ie, js, je, 0, nz},
                 KSixthOrderV{field, mut, msfvx, msfvy, c1, c2, tendency, diff_6th_coef, p.diff_6th_opt});
  } else if (name == 'w') {
    int is = p.specified ? 3 : 0;
    int ie = p.specified ? nx - 3 : nx;
    int js = p.specified ? 3 : 0;
    int je = p.specified ? ny - 3 : ny;
    if (p.periodic_x) { is = 0; ie = nx; }

    parallel_for(Range3D{is, ie, js, je, 1, nz},
                 KSixthOrderM{field, mut, msftx, msfty, c1, c2, tendency, diff_6th_coef, p.diff_6th_opt});
  } else {  // 'm'
    int is = p.specified ? 3 : 0;
    int ie = p.specified ? nx - 3 : nx;
    int js = p.specified ? 3 : 0;
    int je = p.specified ? ny - 3 : ny;
    if (p.periodic_x) { is = 0; ie = nx; }

    parallel_for(Range3D{is, ie, js, je, 0, nz},
                 KSixthOrderM{field, mut, msftx, msfty, c1, c2, tendency, diff_6th_coef, p.diff_6th_opt});
  }
}

}  // namespace cppwrf::dyn::diffusion
