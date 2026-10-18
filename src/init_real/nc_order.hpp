#pragma once
// Fortran-ordered real.exe arrays (farray.hpp) to NetCDF order.
#include <vector>

#include "init_real/farray.hpp"

namespace cppwrf::init_real {

// reorder Fortran (i,k,j) storage to NetCDF (k, j, i) for i in [1,ni], k in [1,nk], j in [1,nj]
template <typename T>
std::vector<double> nc_order(const F3<T>& f, int ni, int nk, int nj) {
  std::vector<double> v(static_cast<std::size_t>(ni) * nk * nj);
  for (int k = 1; k <= nk; ++k)
    for (int j = 1; j <= nj; ++j)
      for (int i = 1; i <= ni; ++i) v[(std::size_t(k - 1) * nj + (j - 1)) * ni + (i - 1)] = double(f(i, k, j));
  return v;
}
template <typename T>
std::vector<double> nc_order(const F2<T>& f, int ni, int nj) {
  std::vector<double> v(static_cast<std::size_t>(ni) * nj);
  for (int j = 1; j <= nj; ++j)
    for (int i = 1; i <= ni; ++i) v[std::size_t(j - 1) * ni + (i - 1)] = double(f(i, j));
  return v;
}
template <typename T>
std::vector<double> nc_order(const F1<T>& f, int n) {
  std::vector<double> v(static_cast<std::size_t>(n));
  for (int k = 1; k <= n; ++k) v[std::size_t(k - 1)] = double(f(k));
  return v;
}

}  // namespace cppwrf::init_real
