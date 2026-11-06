// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once

// ============================================================================
// CPPWRF Coupling and Uncoupling Kernels (couple_or_uncouple_em)
// Ports WRF 4.8.0 dyn_em/couple_or_uncouple_em.F
// Couples/uncouples prognostic state variables in-place before and after
// nesting forcing/feedback, using WRF's exact reciprocal multiplication on
// uncoupling so that interior roundoff matches WRF bit-for-bit:
//   PH  coupled with mutf = (c1f * mub + c2f) + c1f * mu_2
//   W   coupled with muwt = mutf / msfty
//   T   coupled with muth = (c1h * mub + c2h) + c1h * mu_2
//   Moist / Scalars coupled with muth
//   U   coupled with muut = 0.5 * ((c1h*mub(i)+c2h) + (c1h*mub(i-1)+c2h) + c1h*mu_2(i) + c1h*mu_2(i-1)) / msfuy
//   V   likewise in j with msfvx (WRF's summation order, couple_or_uncouple_em.F:143-144)
// ============================================================================

#include <cmath>
#include <vector>
#include <cstdint>
#include "core/backend.hpp"
#include "core/field.hpp"
#include "dyn/state.hpp"

namespace cppwrf::nest {

// ============================================================================
// Single-Source Device / Host Coupling & Uncoupling on CPPWRF State
// Matches WRF dyn_em/couple_or_uncouple_em.F in-place on u_2, v_2, w_2, ph_2,
// t_2, and all scalars q_2[0..n_scalars-1].
// ============================================================================
struct KCoupleFullLevels {
  View3<real> ph_2, w_2;
  BView2<const real> mub, msfty;
  View2<const real> mu_2;
  View1<const real> c1f, c2f;
  bool couple;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    const real denom = (c1f(k) * mub(i, j) + c2f(k)) + (c1f(k) * mu_2(i, j));
    if (couple) {
      const real mutf = denom;
      const real muwt = denom / msfty(i, j);
      ph_2(i, j, k) *= mutf;
      w_2(i, j, k)  *= muwt;
    } else {
      const real mutf = real(1) / denom;
      const real muwt = msfty(i, j) / denom;
      ph_2(i, j, k) *= mutf;
      w_2(i, j, k)  *= muwt;
    }
  }
};

struct KCoupleHalfLevels {
  View3<real> f;
  BView2<const real> mub;
  View2<const real> mu_2;
  View1<const real> c1h, c2h;
  bool couple;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    const real denom = (c1h(k) * mub(i, j) + c2h(k)) + (c1h(k) * mu_2(i, j));
    const real muth = couple ? denom : (real(1) / denom);
    f(i, j, k) *= muth;
  }
};

struct KCoupleU {
  View3<real> u_2;
  BView2<const real> mub, msfuy;
  View2<const real> mu_2;
  View1<const real> c1h, c2h;
  int nx;
  bool couple;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    real muut;
    if (i == nx) {
      const real denom = (c1h(k) * mub(nx - 1, j) + c2h(k)) + (c1h(k) * mu_2(nx - 1, j));
      muut = couple ? (denom / msfuy(nx, j)) : (msfuy(nx, j) / denom);
    } else {
      const int im1 = i - 1;  // WRF reads the halo column at i = ids (mu as left by its halo schedule)
      const real denom = (c1h(k) * mub(i, j) + c2h(k)) + (c1h(k) * mub(im1, j) + c2h(k)) +
                         (c1h(k) * mu_2(i, j)) + (c1h(k) * mu_2(im1, j));
      muut = couple ? (real(0.5) * denom / msfuy(i, j)) : (real(2) * msfuy(i, j) / denom);
    }
    u_2(i, j, k) *= muut;
  }
};

struct KCoupleV {
  View3<real> v_2;
  BView2<const real> mub, msfvx;
  View2<const real> mu_2;
  View1<const real> c1h, c2h;
  int ny;
  bool couple;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    real muvt;
    if (j == ny) {
      const real denom = (c1h(k) * mub(i, ny - 1) + c2h(k)) + (c1h(k) * mu_2(i, ny - 1));
      muvt = couple ? (denom / msfvx(i, ny)) : (msfvx(i, ny) / denom);
    } else {
      const int jm1 = j - 1;  // WRF reads the halo row at j = jds
      const real denom = (c1h(k) * mub(i, j) + c2h(k)) + (c1h(k) * mub(i, jm1) + c2h(k)) +
                         (c1h(k) * mu_2(i, j)) + (c1h(k) * mu_2(i, jm1));
      muvt = couple ? (real(0.5) * denom / msfvx(i, j)) : (real(2) * msfvx(i, j) / denom);
    }
    v_2(i, j, k) *= muvt;
  }
};

inline void couple_or_uncouple_state(
    const Grid& g,
    dyn::State& s,
    const dyn::BaseState& b,
    bool couple) {

  const int nx = g.nx, ny = g.ny, nz = g.nz;
  parallel_for(Range3D{0, nx, 0, ny, 0, nz + 1},
               KCoupleFullLevels{s.ph_2.view(), s.w_2.view(), b.mub.view(), b.msfty.view(),
                                 s.mu_2.cview(), b.c1f.view(), b.c2f.view(), couple});
  parallel_for(Range3D{0, nx, 0, ny, 0, nz},
               KCoupleHalfLevels{s.t_2.view(), b.mub.view(), s.mu_2.cview(),
                                 b.c1h.view(), b.c2h.view(), couple});
  for (int n = 0; n < s.n_scalars(); ++n) {
    if (n == s.P_TKE) continue;  // couple_or_uncouple_em: moist, chem, scalar, tracer; tke_2 is not coupled
    parallel_for(Range3D{0, nx, 0, ny, 0, nz},
                 KCoupleHalfLevels{s.q_2[n].view(), b.mub.view(), s.mu_2.cview(),
                                   b.c1h.view(), b.c2h.view(), couple});
  }
  parallel_for(Range3D{0, nx + 1, 0, ny, 0, nz},
               KCoupleU{s.u_2.view(), b.mub.view(), b.msfuy.view(), s.mu_2.cview(),
                        b.c1h.view(), b.c2h.view(), nx, couple});
  parallel_for(Range3D{0, nx, 0, ny + 1, 0, nz},
               KCoupleV{s.v_2.view(), b.mub.view(), b.msfvx.view(), s.mu_2.cview(),
                        b.c1h.view(), b.c2h.view(), ny, couple});
}

}  // namespace cppwrf::nest
