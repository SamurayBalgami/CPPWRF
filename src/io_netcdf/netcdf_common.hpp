#pragma once

#include <cstddef>
#include <netcdf.h>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cppwrf::io_netcdf {

#define NC_CHECK(call)                                                                     \
  do {                                                                                     \
    int _nc_err = (call);                                                                  \
    if (_nc_err != NC_NOERR) {                                                             \
      throw std::runtime_error(std::string("NetCDF error: ") + nc_strerror(_nc_err) +       \
                               " at " + __FILE__ + ":" + std::to_string(__LINE__));        \
    }                                                                                      \
  } while (0)

enum class DType { F32, F64 };

struct VarView {
  std::string name;
  char stagger = '\0';   // 'M' or '\0' (unstaggered), 'U' (X), 'V' (Y), 'W' (Z)
  int nx = 0, ny = 0, nz = 0; // logical extents (1 for 1D/2D if dimension not used)
  const void* data = nullptr; // pointer to host memory (contiguous or strided)
  DType dtype = DType::F32;
  std::ptrdiff_t sj = 0, sk = 0; // element strides for j and k; i stride is always 1

  VarView() = default;
  VarView(std::string n, char stag, int x, int y, int z, const void* d, DType dt,
          std::ptrdiff_t s_j = 0, std::ptrdiff_t s_k = 0)
      : name(std::move(n)), stagger(stag), nx(x), ny(y), nz(z), data(d), dtype(dt), sj(s_j), sk(s_k) {}
};

class ScopedNcFile {
 public:
  explicit ScopedNcFile(int id = -1) : id_(id) {}
  ~ScopedNcFile() { close(); }
  ScopedNcFile(const ScopedNcFile&) = delete;
  ScopedNcFile& operator=(const ScopedNcFile&) = delete;
  ScopedNcFile(ScopedNcFile&& o) noexcept : id_(o.id_) { o.id_ = -1; }
  ScopedNcFile& operator=(ScopedNcFile&& o) noexcept {
    if (this != &o) {
      close();
      id_ = o.id_;
      o.id_ = -1;
    }
    return *this;
  }
  operator int() const { return id_; }
  int* ptr() { return &id_; }
  void close() {
    if (id_ >= 0) {
      nc_close(id_);
      id_ = -1;
    }
  }
  int release() {
    int ret = id_;
    id_ = -1;
    return ret;
  }
 private:
  int id_ = -1;
};

}  // namespace cppwrf::io_netcdf
