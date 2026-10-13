#pragma once

#include <algorithm>
#include <cmath>
#include "core/backend.hpp"
#include "core/constants.hpp"
#include "core/field.hpp"
#include "dyn/diffusion/diffusion_params.hpp"

namespace cppwrf::dyn::diffusion {

// Ported from WRF v4.6.1 dyn_em/module_big_step_utilities_em.F:
//   w_damp (lines 2503-2711)

struct KWdamp {
  View3<real> rw_tend;
  View3<const real> ww, w;
  View2<const real> mut;
  View1<const real> c1f, c2f, rdnw;
  real dt;
  real w_alpha, w_crit_cfl, w_damp_on;
  int w_damping;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    real c1c2 = c1f(k) * mut(i, j) + c2f(k);
    real vert_cfl = std::abs(ww(i, j, k) / c1c2 * rdnw(k) * dt);

    if (w_damping == 1 && vert_cfl > w_damp_on) {
      real sgn_w = (w(i, j, k) >= real(0)) ? real(1) : real(-1);
      rw_tend(i, j, k) -= sgn_w * w_alpha * (vert_cfl - w_crit_cfl) * c1c2;
    }
  }
};

inline void w_damp(
    View3<real> rw_tend,
    real& max_vert_cfl,
    real& max_horiz_cfl,
    View3<const real> u,
    View3<const real> v,
    View3<const real> ww,
    View3<const real> w,
    View2<const real> mut,
    View1<const real> c1f,
    View1<const real> c2f,
    View1<const real> rdnw,
    real rdx, real rdy,
    BView2<const real> msfux,
    BView2<const real> msfuy,
    BView2<const real> msfvx,
    BView2<const real> msfvy,
    real dt,
    const DiffusionParams& params,
    int nx, int ny, int nz) {

  const real w_crit_cfl = params.w_crit_cfl;
  const real w_damp_on = (params.zadvect_implicit > 0) ? w_crit_cfl : params.w_beta;
  constexpr double pi = 3.1415926535897932384626433;
  constexpr double degrad = pi / 180.0;
  const real msfxffl = params.polar ? real(1.0 / std::cos(params.fft_filter_lat * degrad)) : real(0.0);

  // Kernel for tendency update: WRF loops i=1..nx, j=1..ny, k=2..kde-1 (0-based: i=0..nx-1, j=0..ny-1, k=1..nz-1)
  if (params.w_damping == 1) {
    KWdamp kw{rw_tend, ww, w, mut, c1f, c2f, rdnw, dt, params.w_alpha, w_crit_cfl, w_damp_on, params.w_damping};
    parallel_for(Range3D{0, nx, 0, ny, 1, nz}, kw);
  }

  // Diagnostics: CFL max calculation
  Reducer red(nz);
  max_vert_cfl = red.max(Range3D{0, nx, 0, ny, 1, nz}, [=] CPPWRF_HD (int i, int j, int k) -> double {
    using std::abs;
    real c1c2 = c1f(k) * mut(i, j) + c2f(k);
    return abs(ww(i, j, k) / c1c2 * rdnw(k) * dt);
  });
  max_horiz_cfl = red.max(Range3D{0, nx, 0, ny, 1, nz}, [=] CPPWRF_HD (int i, int j, int k) -> double {
    using std::abs; using std::max; using std::min;
    real msfuxt = params.polar ? min(msfux(i, j), msfxffl) : msfux(i, j);
    return max(abs(u(i, j, k) * rdx * msfuxt * dt),
               abs(v(i, j, k) * rdy * msfvy(i, j) * dt));
  });
}

}  // namespace cppwrf::dyn::diffusion
