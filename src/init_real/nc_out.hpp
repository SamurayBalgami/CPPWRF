#pragma once
// Minimal NetCDF writer for WRF-format files (one Time record): dimensions, global attributes and
// variables with WRF's per-variable attributes.
#include <map>
#include <string>
#include <vector>

#include "init_real/farray.hpp"

namespace cppwrf::init_real {

class NcOut {
 public:
  explicit NcOut(const std::string& path);
  ~NcOut();
  NcOut(const NcOut&) = delete;
  NcOut& operator=(const NcOut&) = delete;

  void dim(const std::string& name, std::size_t len);  // len 0: unlimited
  void att(const std::string& name, const std::string& v);
  void att(const std::string& name, int v);
  void att_f(const std::string& name, float v);
  void att_f(const std::string& name, const std::vector<float>& v);
  // define a variable; dims slowest first, e.g. {"Time","bottom_top","south_north","west_east"}
  void var(const std::string& name, char type /* 'f','d','i','c' */, const std::vector<std::string>& dims,
           const std::map<std::string, std::string>& satts, int field_type);
  void enddef();
  void put(const std::string& name, const std::vector<double>& v, std::size_t rec = 0);  // converted to the type
  void put_int(const std::string& name, const std::vector<int>& v);
  void put_text(const std::string& name, const std::string& s, std::size_t rec = 0);

 private:
  int nc_ = -1;
  std::map<std::string, int> dims_, vars_;
  std::map<std::string, char> types_;
};

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
