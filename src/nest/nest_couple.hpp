#pragma once

// ============================================================================
// CPPWRF Coupling and Uncoupling Kernels (couple_or_uncouple_em)
// Ports WRF 4.6.1 dyn_em/couple_or_uncouple_em.F
// Couples/uncouples prognostic state variables in-place before and after
// nesting forcing/feedback, using WRF's exact reciprocal multiplication on
// uncoupling so that interior roundoff matches WRF bit-for-bit:
//   PH  coupled with mutf = (c1f * mub + c2f) + c1f * mu_2
//   W   coupled with muwt = mutf / msfty
//   T   coupled with muth = (c1h * mub + c2h) + c1h * mu_2
//   Moist / Scalars coupled with muth
//   U   coupled with muut = 0.5 * (muth(i) + muth(i-1)) / msfuy
//   V   coupled with muvt = 0.5 * (muth(j) + muth(j-1)) / msfvx
// ============================================================================

#include <cmath>
#include <vector>
#include <cstdint>
#include "core/backend.hpp"
#include "core/field.hpp"
#include "dyn/state.hpp"
#include "nest/nest_types.hpp"

namespace cppwrf::nest {

// Compute column mass terms mutf (full levels) and muth (half levels)
CPPWRF_HD inline void compute_column_mut(
    double mub, double mu_2,
    const double* c1f, const double* c2f,
    const double* c1h, const double* c2h,
    int nz, double* mutf, double* muth) {

  for (int k = 0; k <= nz; ++k) {
    mutf[k] = (c1f[k] * mub + c2f[k]) + (c1f[k] * mu_2);
  }
  for (int k = 0; k < nz; ++k) {
    muth[k] = (c1h[k] * mub + c2h[k]) + (c1h[k] * mu_2);
  }
}

// ============================================================================
// Routine couple_or_uncouple_em matching WRF dyn_em/couple_or_uncouple_em.F
// Modifies fields in-place using exact reciprocal multiplication on uncoupling.
// ============================================================================
inline void couple_or_uncouple_em(
    double* u_2, double* v_2, double* w_2, double* ph_2, double* t_2,
    double* qv,
    const double* mub, const double* mu_2,
    const double* msfuy, const double* msfvx, const double* msfty,
    const double* c1f, const double* c2f,
    const double* c1h, const double* c2h,
    int nx, int ny, int nz,
    bool couple,
    bool /*is_nested_or_specified*/ = true) {

  const int ni = nx;
  const int nj = ny;
  const int nk = nz;

  auto idx2d = [&](int i, int j) -> std::size_t {
    return static_cast<std::size_t>(i + j * ni);
  };
  auto idx2d_u = [&](int i, int j) -> std::size_t {
    return static_cast<std::size_t>(i + j * (ni + 1));
  };
  auto idx2d_v = [&](int i, int j) -> std::size_t {
    return static_cast<std::size_t>(i + j * ni);
  };
  auto idx3d = [&](int i, int k, int j, int n_i, int n_k) -> std::size_t {
    return static_cast<std::size_t>(i + k * n_i + j * n_i * n_k);
  };

  // 1. Mass-point full-level fields: PH_2 and W_2
  for (int j = 0; j < nj; ++j) {
    for (int k = 0; k <= nk; ++k) {
      for (int i = 0; i < ni; ++i) {
        double mb = mub[idx2d(i, j)];
        double m2 = mu_2[idx2d(i, j)];
        double denom = (c1f[k] * mb + c2f[k]) + (c1f[k] * m2);
        std::size_t p_idx = idx3d(i, k, j, ni, nk + 1);

        if (couple) {
          double mutf = denom;
          double muwt = denom / msfty[idx2d(i, j)];
          if (ph_2) ph_2[p_idx] *= mutf;
          if (w_2)  w_2[p_idx]  *= muwt;
        } else {
          double mutf = 1.0 / denom;
          double muwt = msfty[idx2d(i, j)] / denom;
          if (ph_2) ph_2[p_idx] *= mutf;
          if (w_2)  w_2[p_idx]  *= muwt;
        }
      }
    }
  }

  // 2. Mass-point half-level fields: T_2 and QV
  for (int j = 0; j < nj; ++j) {
    for (int k = 0; k < nk; ++k) {
      for (int i = 0; i < ni; ++i) {
        double mb = mub[idx2d(i, j)];
        double m2 = mu_2[idx2d(i, j)];
        double denom = (c1h[k] * mb + c2h[k]) + (c1h[k] * m2);
        double muth = couple ? denom : (1.0 / denom);
        std::size_t p_idx = idx3d(i, k, j, ni, nk);

        if (t_2) t_2[p_idx] *= muth;
        if (qv)  qv[p_idx]  *= muth;
      }
    }
  }

  // 3. U_2 (x-momentum, ni + 1 points in x)
  if (u_2) {
    for (int j = 0; j < nj; ++j) {
      for (int k = 0; k < nk; ++k) {
        for (int i = 0; i <= ni; ++i) {
          double muut;
          if (couple) {
            if (i == 0) {
              double muth0 = (c1h[k] * mub[idx2d(0, j)] + c2h[k]) + (c1h[k] * mu_2[idx2d(0, j)]);
              muut = muth0 / msfuy[idx2d_u(0, j)];
            } else if (i == ni) {
              double muth_n = (c1h[k] * mub[idx2d(ni - 1, j)] + c2h[k]) + (c1h[k] * mu_2[idx2d(ni - 1, j)]);
              muut = muth_n / msfuy[idx2d_u(ni, j)];
            } else {
              double muth_0  = (c1h[k] * mub[idx2d(i, j)]     + c2h[k]) + (c1h[k] * mu_2[idx2d(i, j)]);
              double muth_m1 = (c1h[k] * mub[idx2d(i - 1, j)] + c2h[k]) + (c1h[k] * mu_2[idx2d(i - 1, j)]);
              muut = 0.5 * (muth_0 + muth_m1) / msfuy[idx2d_u(i, j)];
            }
          } else {
            if (i == 0) {
              double muth0_inv = 1.0 / ((c1h[k] * mub[idx2d(0, j)] + c2h[k]) + (c1h[k] * mu_2[idx2d(0, j)]));
              muut = msfuy[idx2d_u(0, j)] * muth0_inv;
            } else if (i == ni) {
              muut = msfuy[idx2d_u(ni, j)] / ((c1h[k] * mub[idx2d(ni - 1, j)] + c2h[k]) + (c1h[k] * mu_2[idx2d(ni - 1, j)]));
            } else {
              double denom = (c1h[k] * mub[idx2d(i, j)] + c2h[k]) + (c1h[k] * mub[idx2d(i - 1, j)] + c2h[k]) +
                             (c1h[k] * mu_2[idx2d(i, j)]) + (c1h[k] * mu_2[idx2d(i - 1, j)]);
              muut = 2.0 * msfuy[idx2d_u(i, j)] / denom;
            }
          }
          std::size_t u_idx = idx3d(i, k, j, ni + 1, nk);
          u_2[u_idx] *= muut;
        }
      }
    }
  }

  // 4. V_2 (y-momentum, nj + 1 points in y)
  if (v_2) {
    for (int j = 0; j <= nj; ++j) {
      for (int k = 0; k < nk; ++k) {
        for (int i = 0; i < ni; ++i) {
          double muvt;
          if (couple) {
            if (j == 0) {
              double muth0 = (c1h[k] * mub[idx2d(i, 0)] + c2h[k]) + (c1h[k] * mu_2[idx2d(i, 0)]);
              muvt = muth0 / msfvx[idx2d_v(i, 0)];
            } else if (j == nj) {
              double muth_n = (c1h[k] * mub[idx2d(i, nj - 1)] + c2h[k]) + (c1h[k] * mu_2[idx2d(i, nj - 1)]);
              muvt = muth_n / msfvx[idx2d_v(i, nj)];
            } else {
              double muth_0  = (c1h[k] * mub[idx2d(i, j)]     + c2h[k]) + (c1h[k] * mu_2[idx2d(i, j)]);
              double muth_m1 = (c1h[k] * mub[idx2d(i, j - 1)] + c2h[k]) + (c1h[k] * mu_2[idx2d(i, j - 1)]);
              muvt = 0.5 * (muth_0 + muth_m1) / msfvx[idx2d_v(i, j)];
            }
          } else {
            if (j == 0) {
              double muth0_inv = 1.0 / ((c1h[k] * mub[idx2d(i, 0)] + c2h[k]) + (c1h[k] * mu_2[idx2d(i, 0)]));
              muvt = msfvx[idx2d_v(i, 0)] * muth0_inv;
            } else if (j == nj) {
              muvt = msfvx[idx2d_v(i, nj)] / ((c1h[k] * mub[idx2d(i, nj - 1)] + c2h[k]) + (c1h[k] * mu_2[idx2d(i, nj - 1)]));
            } else {
              double denom = (c1h[k] * mub[idx2d(i, j)] + c2h[k]) + (c1h[k] * mub[idx2d(i, j - 1)] + c2h[k]) +
                             (c1h[k] * mu_2[idx2d(i, j)]) + (c1h[k] * mu_2[idx2d(i, j - 1)]);
              muvt = 2.0 * msfvx[idx2d_v(i, j)] / denom;
            }
          }
          std::size_t v_idx = idx3d(i, k, j, ni, nk);
          v_2[v_idx] *= muvt;
        }
      }
    }
  }
}

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
      const int im1 = (i == 0) ? 0 : (i - 1);
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
      const int jm1 = (j == 0) ? 0 : (j - 1);
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
