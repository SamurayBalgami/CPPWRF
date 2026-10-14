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
  CPPWRF_HD double& operator[](int i) const { return p[std::ptrdiff_t(i) * s]; }
};

struct SI {
  int* p;
  int s;
  CPPWRF_HD int& operator[](int i) const { return p[std::ptrdiff_t(i) * s]; }
};

struct SD2 {
  double* p;
  int s;
  int n;
  CPPWRF_HD SD operator[](int a) const { return SD{p + std::ptrdiff_t(a) * n * s, s}; }
};

// views of local stride-1 arrays
CPPWRF_HD inline SD sd(double* p) { return SD{p, 1}; }
CPPWRF_HD inline SI si(int* p) { return SI{p, 1}; }
CPPWRF_HD inline SD2 sd2(double* p, int n) { return SD2{p, 1, n}; }

}  // namespace cppwrf::phys
