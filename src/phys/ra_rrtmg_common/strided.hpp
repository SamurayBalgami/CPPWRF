#pragma once
// Strided per-column array views for the RRTMG ports. One code path serves both layouts:
//  - a single column in local arrays (stride 1, the CPU column driver);
//  - a batch of columns in device memory, column index fastest (stride = number of columns), so a
//    warp of threads working on neighbouring columns reads and writes coalesced.
// SD/SI: x[i] = p[i*s]; SD2: x[a][b] = p[(a*n + b)*s] (row length n).
#include <cstddef>

#include "core/backend.hpp"

namespace cppwrf::phys {

struct SD {
  double* p;
  int s;
  CPPWRF_FORCEINLINE double& operator[](int i) const { return p[std::ptrdiff_t(i) * s]; }
};

struct SI {
  int* p;
  int s;
  CPPWRF_FORCEINLINE int& operator[](int i) const { return p[std::ptrdiff_t(i) * s]; }
};

struct SD2 {
  double* p;
  int s;
  int n;
  CPPWRF_FORCEINLINE SD operator[](int a) const { return SD{p + std::ptrdiff_t(a) * n * s, s}; }
};

struct CSD {
  const double* p;
  int s;
  CPPWRF_FORCEINLINE CSD() : p(nullptr), s(0) {}
  CPPWRF_FORCEINLINE CSD(const double* p_, int s_) : p(p_), s(s_) {}
  CPPWRF_FORCEINLINE CSD(SD v) : p(v.p), s(v.s) {}
  CPPWRF_FORCEINLINE double operator[](int i) const {
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
  const double* p;
  int s;
  int n;
  CPPWRF_FORCEINLINE CSD2() : p(nullptr), s(0), n(0) {}
  CPPWRF_FORCEINLINE CSD2(const double* p_, int s_, int n_) : p(p_), s(s_), n(n_) {}
  CPPWRF_FORCEINLINE CSD2(SD2 v) : p(v.p), s(v.s), n(v.n) {}
  CPPWRF_FORCEINLINE CSD operator[](int a) const { return CSD{p + std::ptrdiff_t(a) * n * s, s}; }
};

// views of local stride-1 arrays
CPPWRF_FORCEINLINE SD sd(double* p) { return SD{p, 1}; }
CPPWRF_FORCEINLINE SI si(int* p) { return SI{p, 1}; }
CPPWRF_FORCEINLINE SD2 sd2(double* p, int n) { return SD2{p, 1, n}; }

}  // namespace cppwrf::phys
