#pragma once
// Reader-independent lateral boundary data structures for specified boundaries.
// Supports West (XS), East (XE), South (YS), North (YE) boundary slabs of width 5,
// storing both base boundary values and boundary tendencies.
// Memory layout matches Fortran column-major storage for direct zero-copy interchange.

#include <cstddef>
#include <vector>
#include "core/backend.hpp"
#include "core/grid.hpp"
#include "core/precision.hpp"

namespace cppwrf::bdy {

// View for a single 3D lateral boundary slab (position along boundary, vertical k, depth b)
template <class T>
struct BdySlabView {
  T* p = nullptr;
  int n_pos = 0;
  int n_k = 0;
  int width = 5;
  int s_pos = 1;
  int s_k = 0;
  int s_b = 0;

  CPPWRF_HD T& operator()(int pos, int k, int b) const {
    return p[pos * s_pos + k * s_k + b * s_b];
  }
};

// Container for a lateral boundary slab (flat allocation, row/column pitch contiguous)
template <class T>
class BdySlab {
 public:
  BdySlab() = default;
  BdySlab(int n_pos, int n_k, int width = 5)
      : n_pos_(n_pos), n_k_(n_k), width_(width), data_(n_pos * n_k * width, T(0)) {}

  void allocate(int n_pos, int n_k, int width = 5) {
    n_pos_ = n_pos;
    n_k_ = n_k;
    width_ = width;
    data_.assign(std::size_t(n_pos) * std::size_t(n_k) * std::size_t(width), T(0));
  }

  int n_pos() const { return n_pos_; }
  int n_k() const { return n_k_; }
  int width() const { return width_; }
  std::size_t size() const { return data_.size(); }
  T* data() { return data_.data(); }
  const T* data() const { return data_.data(); }

  BdySlabView<T> view() {
    return BdySlabView<T>{data_.data(), n_pos_, n_k_, width_, 1, n_pos_, n_pos_ * n_k_};
  }
  BdySlabView<const T> cview() const {
    return BdySlabView<const T>{data_.data(), n_pos_, n_k_, width_, 1, n_pos_, n_pos_ * n_k_};
  }

  T& operator()(int pos, int k, int b) {
    return data_[pos + k * n_pos_ + b * (n_pos_ * n_k_)];
  }
  const T& operator()(int pos, int k, int b) const {
    return data_[pos + k * n_pos_ + b * (n_pos_ * n_k_)];
  }

 private:
  int n_pos_ = 0;
  int n_k_ = 0;
  int width_ = 5;
  std::vector<T> data_;
};

// Holds the 4 boundaries (XS, XE, YS, YE) for a single physical variable (values + tendencies)
template <class T>
struct BdyVar {
  BdySlab<T> bxs;   // West values (j, k, b)
  BdySlab<T> bxe;   // East values (j, k, b)
  BdySlab<T> bys;   // South values (i, k, b)
  BdySlab<T> bye;   // North values (i, k, b)

  BdySlab<T> btxs;  // West tendencies (j, k, b)
  BdySlab<T> btxe;  // East tendencies (j, k, b)
  BdySlab<T> btys;  // South tendencies (i, k, b)
  BdySlab<T> btye;  // North tendencies (i, k, b)

  void allocate(int nx, int ny, int nz, int width = 5, Stagger stag = Stagger::M) {
    int nx_s = (stag == Stagger::U) ? nx + 1 : nx;
    int ny_s = (stag == Stagger::V) ? ny + 1 : ny;
    int nz_s = (stag == Stagger::W) ? nz + 1 : nz;

    bxs.allocate(ny_s, nz_s, width);
    bxe.allocate(ny_s, nz_s, width);
    bys.allocate(nx_s, nz_s, width);
    bye.allocate(nx_s, nz_s, width);

    btxs.allocate(ny_s, nz_s, width);
    btxe.allocate(ny_s, nz_s, width);
    btys.allocate(nx_s, nz_s, width);
    btye.allocate(nx_s, nz_s, width);
  }
};

// Reader-independent lateral boundary condition dataset for a domain
template <class T = real>
struct BdyData {
  int nx = 0, ny = 0, nz = 0;
  int width = 5;

  BdyVar<T> u;   // Stagger::U
  BdyVar<T> v;   // Stagger::V
  BdyVar<T> w;   // Stagger::W
  BdyVar<T> t;   // Stagger::M
  BdyVar<T> ph;  // Stagger::W
  BdyVar<T> mu;  // Stagger::M, nz = 1 (2D dry-air column mass)

  std::vector<BdyVar<T>> scalars; // optional scalars (qv, qc, ...)

  void init(int nx_, int ny_, int nz_, int width_ = 5) {
    nx = nx_;
    ny = ny_;
    nz = nz_;
    width = width_;

    u.allocate(nx, ny, nz, width, Stagger::U);
    v.allocate(nx, ny, nz, width, Stagger::V);
    w.allocate(nx, ny, nz, width, Stagger::W);
    t.allocate(nx, ny, nz, width, Stagger::M);
    ph.allocate(nx, ny, nz, width, Stagger::W);
    mu.allocate(nx, ny, 1, width, Stagger::M);
  }

  void add_scalar() {
    scalars.emplace_back();
    scalars.back().allocate(nx, ny, nz, width, Stagger::M);
  }
};

}  // namespace cppwrf::bdy
