// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "phys/ra_rrtmg_common/rrtmg_real.hpp"
// Strided per-column array views for the RRTMG ports. One code path serves both layouts:
//  - a single column in local arrays (stride 1, the CPU column driver);
//  - a batch of columns in device memory, column index fastest (stride = number of columns), so a
//    warp of threads working on neighbouring columns reads and writes coalesced.
// SD/SI: x[i] = p[i*s]; SD2: x[a][b] = p[(a*n + b)*s] (row length n).
#include <cstddef>

#include "core/backend.hpp"

namespace cppwrf::phys {

struct SD {
  rr* p;
  int s;
  CPPWRF_FORCEINLINE rr& operator[](int i) const { return p[std::ptrdiff_t(i) * s]; }
};

struct SI {
  int* p;
  int s;
  CPPWRF_FORCEINLINE int& operator[](int i) const { return p[std::ptrdiff_t(i) * s]; }
};

// read-only arrays seen through the (non-const) views: one column (stride 1) and column col of a
// batch of ncol columns
CPPWRF_FORCEINLINE SD sd_const(const rr* p) { return SD{const_cast<rr*>(p), 1}; }
CPPWRF_FORCEINLINE SD col_view(const rr* a, int col, int ncol) { return SD{const_cast<rr*>(a) + col, ncol}; }

// grows a workspace's rr and int buffers to a layout's sizes (lay.nd, lay.ni)
template <class Workspace, class Layout>
inline void ensure_workspace(Workspace& ws, const Layout& lay) {
  if (ws.d.size() < lay.nd) ws.d.assign(lay.nd, rr(0.0));
  if (ws.i.size() < lay.ni) ws.i.assign(lay.ni, 0);
}

struct SD2 {
  rr* p;
  int s;
  int n;
  CPPWRF_FORCEINLINE SD operator[](int a) const { return SD{p + std::ptrdiff_t(a) * n * s, s}; }
};

struct CSD {
  const rr* p;
  int s;
  CPPWRF_FORCEINLINE CSD() : p(nullptr), s(0) {}
  CPPWRF_FORCEINLINE CSD(const rr* p_, int s_) : p(p_), s(s_) {}
  CPPWRF_FORCEINLINE CSD(SD v) : p(v.p), s(v.s) {}
  CPPWRF_FORCEINLINE rr operator[](int i) const {
#if defined(__CUDA_ARCH__)
    return __ldg(p + std::ptrdiff_t(i) * s);
#else
    return p[std::ptrdiff_t(i) * s];
#endif
  }
};

struct CSI {
  const int* p;
  int s;
  CPPWRF_FORCEINLINE CSI() : p(nullptr), s(0) {}
  CPPWRF_FORCEINLINE CSI(const int* p_, int s_) : p(p_), s(s_) {}
  CPPWRF_FORCEINLINE CSI(SI v) : p(v.p), s(v.s) {}
  CPPWRF_FORCEINLINE int operator[](int i) const {
#if defined(__CUDA_ARCH__)
    return __ldg(p + std::ptrdiff_t(i) * s);
#else
    return p[std::ptrdiff_t(i) * s];
#endif
  }
};

struct CSD2 {
  const rr* p;
  int s;
  int n;
  CPPWRF_FORCEINLINE CSD2() : p(nullptr), s(0), n(0) {}
  CPPWRF_FORCEINLINE CSD2(const rr* p_, int s_, int n_) : p(p_), s(s_), n(n_) {}
  CPPWRF_FORCEINLINE CSD2(SD2 v) : p(v.p), s(v.s), n(v.n) {}
  CPPWRF_FORCEINLINE CSD operator[](int a) const { return CSD{p + std::ptrdiff_t(a) * n * s, s}; }
};

// view of a local stride-1 array
CPPWRF_FORCEINLINE SD sd(rr* p) { return SD{p, 1}; }

}  // namespace cppwrf::phys
