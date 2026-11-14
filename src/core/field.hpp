// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Fields: SoA, one allocation per field, i fastest, k slowest, halo 3 in x (and y if ny>1),
// row pitch padded to a multiple of 32 elements, i=0 aligned to 8 elements.
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "core/backend.hpp"
#include "core/grid.hpp"
#include "core/precision.hpp"

namespace cppwrf {

// View strides are 32-bit: index arithmetic stays in 32 bits and is sign-extended once (all
// fields are < 2^31 elements). 64-bit strides triggered a ptxas (CUDA 12.4, sm_75) miscompile
// of negative 64-bit offsets ([R.64+UR] addressing) in KTendPh: see docs/STATUS.md.
template <class T> struct View3 {
  T* p = nullptr;
  int sj = 0, sk = 0;
  CPPWRF_HD T& operator()(int i, int j, int k) const { return p[i + j * sj + k * sk]; }
};
template <class T> struct View2 {
  T* p = nullptr;
  int sj = 0;
  CPPWRF_HD T& operator()(int i, int j) const { return p[i + j * sj]; }
};
template <class T> struct View1 {
  T* p = nullptr;
  CPPWRF_HD T& operator()(int k) const { return p[k]; }
};
// Base-state views carry an i-stride too: si = sj = 0 when horizontally uniform (flat
// terrain), so 1D storage is read through 3D indexing and terrain needs no kernel rewrite.
template <class T> struct BView3 {
  T* p = nullptr;
  int si = 0, sj = 0, sk = 0;
  CPPWRF_HD T& operator()(int i, int j, int k) const { return p[i * si + j * sj + k * sk]; }
};
template <class T> struct BView2 {
  T* p = nullptr;
  int si = 0, sj = 0;
  CPPWRF_HD T& operator()(int i, int j) const { return p[i * si + j * sj]; }
};

struct Layout {
  int nx = 0, ny = 0, nz = 0;  // logical extents (staggered, no halo)
  int hy = 0;
  std::ptrdiff_t si = 1, sj = 0, sk = 0, offset = 0;
  std::size_t total = 0;
  // the indices the allocation holds: wi0 <= i < wi1, wj0 <= j < wj1 (the halo included). The whole domain:
  // [-halo, nx + halo) x [-hy, ny + hy); a patch of a decomposed run (dom >= 0): its points and halo
  int wi0 = 0, wi1 = 0, wj0 = 0, wj1 = 0;
  int dom = -1;
  // the field's points (no lateral halo) in the allocation: host loops over the domain run over these
  int i0() const { return wi0 > 0 ? wi0 : 0; }
  int i1() const { return wi1 < nx ? wi1 : nx; }
  int j0() const { return wj0 > 0 ? wj0 : 0; }
  int j1() const { return wj1 < ny ? wj1 : ny; }
  bool holds(int i, int j) const { return si == 0 || (i >= wi0 && i < wi1 && (sj == 0 || (j >= wj0 && j < wj1))); }
  static constexpr int lpad = 8;
  static Layout make(const Grid& g, Stagger s, int rank, bool uniform = false);
  std::size_t index(int i, int j, int k) const { return std::size_t(offset + i * si + j * sj + k * sk); }
};

namespace tile {
// The registry of a tile thread's fields (src/core/tiles.cpp): every field the thread allocates, in allocation order,
// so that the n-th field of one tile is the n-th field of every other tile (they run the same code). Null without tiles.
struct Track {
  void (*add)(const void* base, const Layout& L, int esz, const std::string& name);
  void (*remove)(const void* base);
};
inline thread_local const Track* track = nullptr;
}  // namespace tile

class Buffer {
 public:
  Buffer() = default;
  explicit Buffer(std::size_t n) : n_(n), p_(static_cast<real*>(dev_alloc(n * sizeof(real)))) {}
  ~Buffer() {
    if (p_ && tile::track) tile::track->remove(p_);
    dev_free(p_);
  }
  Buffer(Buffer&& o) noexcept : n_(o.n_), p_(o.p_) { o.p_ = nullptr; o.n_ = 0; }
  Buffer& operator=(Buffer&& o) noexcept { std::swap(n_, o.n_); std::swap(p_, o.p_); return *this; }
  Buffer(const Buffer&) = delete;
  Buffer& operator=(const Buffer&) = delete;
  real* get() const { return p_; }
  std::size_t size() const { return n_; }

 private:
  std::size_t n_ = 0;
  real* p_ = nullptr;
};

// Prognostic/diagnostic field (rank 3 or 2). Move-only; swap() swaps storage, not names.
template <int R> class FieldT {
 public:
  FieldT() = default;
  FieldT(const Grid& g, Stagger s, std::string name)
      : name_(std::move(name)), st_(s), L_(Layout::make(g, s, R)), buf_(L_.total) {
    if (tile::track) tile::track->add(buf_.get(), L_, int(sizeof(real)), name_);
  }
  const std::string& name() const { return name_; }
  const void* base() const { return buf_.get(); }  // the allocation (halo exchange between tiles)
  Stagger stagger() const { return st_; }
  const Layout& layout() const { return L_; }
  real* data() const { return buf_.get() + L_.offset; }  // element (0,0,0)
  std::size_t bytes() const { return L_.total * sizeof(real); }
  auto view() const {
    if constexpr (R == 3) return View3<real>{data(), int(L_.sj), int(L_.sk)};
    else return View2<real>{data(), int(L_.sj)};
  }
  auto cview() const {
    if constexpr (R == 3) return View3<const real>{data(), int(L_.sj), int(L_.sk)};
    else return View2<const real>{data(), int(L_.sj)};
  }
  // Host-side initialisation (setup only): f(i,j,k) for logical points; halos zero.
  template <class F> void fill_host(F f) {
    std::vector<real> h(L_.total, real(0));
    for (int k = 0; k < L_.nz; ++k)
      for (int j = L_.j0(); j < L_.j1(); ++j)
        for (int i = L_.i0(); i < L_.i1(); ++i) h[L_.index(i, j, k)] = real(f(i, j, k));
    copy_to_dev(buf_.get(), h.data(), bytes());
  }
  void download(std::vector<real>& h) const {  // whole buffer incl. halos
    h.resize(L_.total);
    copy_to_host(h.data(), buf_.get(), bytes());
  }
  void copy_from(const FieldT& o) { copy_dev(buf_.get(), o.buf_.get(), bytes()); }
  // Swaps the storage (data, layout, stagger); each field keeps its name.
  friend void swap(FieldT& a, FieldT& b) noexcept {
    std::swap(a.st_, b.st_); std::swap(a.L_, b.L_); std::swap(a.buf_, b.buf_);
  }

 private:
  std::string name_;
  Stagger st_ = Stagger::M;
  Layout L_;
  Buffer buf_;
};
using Field3 = FieldT<3>;
using Field2 = FieldT<2>;

// Base-state field: stored 1D (k only) / scalar when horizontally uniform.
template <int R> class BaseFieldT {
 public:
  BaseFieldT() = default;
  BaseFieldT(const Grid& g, Stagger s, std::string name, bool uniform)
      : name_(std::move(name)), st_(s), uniform_(uniform), L_(Layout::make(g, s, R, uniform)), buf_(L_.total) {
    if (tile::track && !uniform_) tile::track->add(buf_.get(), L_, int(sizeof(real)), name_);
  }
  bool uniform() const { return uniform_; }
  const void* base() const { return buf_.get(); }
  const Layout& layout() const { return L_; }
  real* data() const { return buf_.get() + L_.offset; }
  auto view() const {
    if constexpr (R == 3) return BView3<const real>{data(), int(L_.si), int(L_.sj), int(L_.sk)};
    else return BView2<const real>{data(), int(L_.si), int(L_.sj)};
  }
  // f(i,j,k); for uniform fields only i=j=0 is evaluated.
  template <class F> void fill_host(F f) {
    std::vector<real> h(L_.total, real(0));
    const int ia = uniform_ ? 0 : L_.i0(), ib = uniform_ ? 1 : L_.i1(), ja = uniform_ ? 0 : L_.j0(), jb = uniform_ ? 1 : L_.j1();
    for (int k = 0; k < L_.nz; ++k)
      for (int j = ja; j < jb; ++j)
        for (int i = ia; i < ib; ++i) h[L_.index(i, j, k)] = real(f(i, j, k));
    copy_to_dev(buf_.get(), h.data(), L_.total * sizeof(real));
  }
  // As fill_host, but a non-uniform field is also filled in its lateral halo:
  // i = -halo..nx+halo-1 (and j likewise when ny > 1); f must accept those indices.
  template <class F> void fill_host_halo(F f) {
    if (uniform_) { fill_host(f); return; }
    std::vector<real> h(L_.total, real(0));
    for (int k = 0; k < L_.nz; ++k)
      for (int j = L_.wj0; j < L_.wj1; ++j)
        for (int i = L_.wi0; i < L_.wi1; ++i) h[L_.index(i, j, k)] = real(f(i, j, k));
    copy_to_dev(buf_.get(), h.data(), L_.total * sizeof(real));
  }
  void download(std::vector<real>& h) const {
    h.resize(L_.total);
    copy_to_host(h.data(), buf_.get(), L_.total * sizeof(real));
  }

 private:
  std::string name_;
  Stagger st_ = Stagger::M;
  bool uniform_ = true;
  Layout L_;
  Buffer buf_;
};
using BaseField3 = BaseFieldT<3>;
using BaseField2 = BaseFieldT<2>;

// 1D vertical array with a host mirror (host copy is authoritative during setup).
class Array1 {
 public:
  Array1() = default;
  explicit Array1(int n) : h_(n, 0.0), buf_(n) {}
  int size() const { return int(h_.size()); }
  double& operator[](int k) { return h_[k]; }
  double operator[](int k) const { return h_[k]; }
  void upload() {
    std::vector<real> t(h_.begin(), h_.end());
    copy_to_dev(buf_.get(), t.data(), t.size() * sizeof(real));
  }
  View1<const real> view() const { return {buf_.get()}; }

 private:
  std::vector<double> h_;
  Buffer buf_;
};

namespace tile {
// Halo exchange of fields between the tiles inside a tiled time step (no-op otherwise): every tile's copy of each field
// gets the values its neighbours own within Grid::halo points of its own patch (and the periodic wrap's sources).
// Collective: every tile calls it with the same fields.
template <class... F>
void exchange(const F&... f) {
  if (!ctx.on || !hooks.exchange) return;
  const void* const b[] = {f.base()...};
  hooks.exchange(b, int(sizeof...(F)), false);
}
template <class F>
void exchange(const std::vector<F>& fs) {
  if (!ctx.on || !hooks.exchange || fs.empty()) return;
  std::vector<const void*> b;
  for (const F& f : fs) b.push_back(f.base());
  hooks.exchange(b.data(), int(b.size()), false);
}
}  // namespace tile

}  // namespace cppwrf
