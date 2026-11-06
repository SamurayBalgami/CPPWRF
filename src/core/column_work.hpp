// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Work arrays of column kernels (column_for: one thread per column) in device memory instead of the thread's stack.
//
// The CUDA driver reserves a kernel's stack frame for every thread the GPU can hold at once (SMs x threads per SM,
// e.g. 188 x 1536 on the RTX PRO 6000), whatever the launch: a scheme that keeps ~170 column arrays of 48 levels
// in FP64 (Thompson, 72 KB per thread) reserves ~20 GB there. With CPPWRF_COLUMN_ARRAY the arrays of a column live
// in one shared device buffer sized for the launch's threads only, laid out as the driver lays out local memory:
// per warp of the launch a block of its 32 threads' arrays, level k of an array of the 32 threads side by side (so
// a warp's access is one contiguous 128 or 256 bytes and a warp's arrays are one contiguous block). The values
// stored and the arithmetic on them are the same as with the local arrays; the CPU build (and the host pass of a
// CUDA build) declares the local arrays as before.
//
// A kernel lists its arrays in an X-macro LIST(X) of X(type, name, levels) entries (X(type, name, levels, n) and
// CPPWRF_COLUMN_ARRAY2 / CPPWRF_COLUMN_BYTES2 for a [levels][n] array), has a member
// `ColumnWork work;`, declares them with
//     CPPWRF_COLUMN_SLOTS(work, i, j)          // once, in operator()(i, j)
//     LIST(CPPWRF_COLUMN_ARRAY)                // T name[n]; on the host, a ColArray<T> into the buffer on the GPU
// and gives its bytes per column as `0 LIST(CPPWRF_COLUMN_BYTES)`. A helper with arrays of its own takes a copy
// of the cursor (`ColumnSlots cppwrf_cw_`) and continues after the caller's arrays. The launcher reserves the
// buffer outside the time loop (column_work_reserve) and launches with column_for_work.
#include <algorithm>
#include <cstddef>

#include "core/backend.hpp"

namespace cppwrf {

// level k of a work array: p[32 k] (the warp's 32 threads interleaved)
template <class T>
struct ColArray {
  T* p;
  CPPWRF_FORCEINLINE T& operator[](int k) const { return p[k * 32]; }
  CPPWRF_FORCEINLINE explicit operator bool() const { return p != nullptr; }  // an optional array (null: absent)
};

// a [n1][N2] work array: element (k, n) of a thread at p[32 (k N2 + n)]
template <class T, int N2>
struct ColArray2 {
  T* p;
  CPPWRF_FORCEINLINE ColArray<T> operator[](int k) const { return ColArray<T>{p + k * N2 * 32}; }
};

// One launch's buffer: warp w of the launch (thread t = w * 32 + lane in column_for's numbering: blocks of
// column_block() threads, x fastest) owns 32 * bpc bytes at base + w * 32 * bpc; in it the arrays one after
// another, each n * sizeof(T) * 32 bytes
struct ColumnWork {
  char* base = nullptr;
  std::size_t bpc = 0;  // bytes per thread (a multiple of 8)
};

class ColumnSlots {
 public:
  CPPWRF_HD ColumnSlots() : base_(nullptr), lane_(0) {}  // a host caller's (local arrays only)
  CPPWRF_HD ColumnSlots(const ColumnWork& w, int, int) : base_(w.base), lane_(0) {
#if defined(__CUDA_ARCH__)
    const unsigned t = ((blockIdx.y * gridDim.x + blockIdx.x) * blockDim.y + threadIdx.y) * blockDim.x + threadIdx.x;
    base_ = w.base + std::size_t(t / 32) * 32 * w.bpc;
    lane_ = int(t % 32);
#endif
  }
  // column (i, j)'s slots in column_for(r, ...) with blocks of bx x by threads (from any kernel)
  CPPWRF_HD ColumnSlots(const ColumnWork& w, int i, int j, const Range2D& r, int bx, int by)
      : base_(w.base), lane_(0) {
    const int ti = i - r.i0, tj = j - r.j0;
    const unsigned gx = unsigned((r.i1 - r.i0 + bx - 1) / bx);
    const unsigned t = ((unsigned(tj / by) * gx + unsigned(ti / bx)) * unsigned(by) + unsigned(tj % by)) * unsigned(bx) +
                       unsigned(ti % bx);
    base_ = w.base + std::size_t(t / 32) * 32 * w.bpc;
    lane_ = int(t % 32);
  }
  template <class T>
  CPPWRF_FORCEINLINE ColArray<T> take(int n) {
    T* p = reinterpret_cast<T*>(base_ + off_ * 32) + lane_;
    off_ += std::size_t(n) * sizeof(T);
    return ColArray<T>{p};
  }
  template <class T, int N2>
  CPPWRF_FORCEINLINE ColArray2<T, N2> take2(int n1) {
    T* p = reinterpret_cast<T*>(base_ + off_ * 32) + lane_;
    off_ += std::size_t(n1) * N2 * sizeof(T);
    return ColArray2<T, N2>{p};
  }

 private:
  char* base_;
  int lane_;
  std::size_t off_ = 0;  // bytes per thread taken so far
};

#if defined(__CUDA_ARCH__)
#define CPPWRF_COLUMN_ARRAY(T, name, n) const ::cppwrf::ColArray<T> name = cppwrf_cw_.take<T>(n);
#define CPPWRF_COLUMN_ARRAY2(T, name, n1, n2) const ::cppwrf::ColArray2<T, n2> name = cppwrf_cw_.take2<T, n2>(n1);
#else
#define CPPWRF_COLUMN_ARRAY(T, name, n) T name[n];
#define CPPWRF_COLUMN_ARRAY2(T, name, n1, n2) T name[n1][n2];
#endif
#define CPPWRF_COLUMN_BYTES(T, name, n) +sizeof(T) * std::size_t(n)
#define CPPWRF_COLUMN_BYTES2(T, name, n1, n2) +sizeof(T) * std::size_t(n1) * std::size_t(n2)
#define CPPWRF_COLUMN_SLOTS(work, i, j) [[maybe_unused]] ::cppwrf::ColumnSlots cppwrf_cw_((work), (i), (j));

// Grows the process's work buffer (outside the time loop) for a kernel of bytes_per_col bytes per column over the
// columns of r: all of them while they fit an eighth of the device memory and a quarter of what is free, else rows
// of them (at least two). The kernels of every scheme and domain run one after another on the one launch stream
// and share it. CPU build: no-op.
void column_work_reserve(const Range2D& r, std::size_t bytes_per_col);

namespace detail {
// the buffer and its size (CUDA; grown to the request when it is smaller: then outside the time loop only if the
// launcher reserved it)
char* column_work_buffer(std::size_t bytes, std::size_t* have);
// bytes per thread: a multiple of 8, so every array starts aligned
constexpr std::size_t column_work_bpc(std::size_t bytes_per_col) { return (bytes_per_col + 7) / 8 * 8; }
// the threads of column_for over ni x nj columns
inline std::size_t column_for_threads(int ni, int nj) {
  const ColumnBlock b = column_block();
  return std::size_t((ni + b.x - 1) / b.x * b.x) * std::size_t((nj + b.y - 1) / b.y * b.y);
}
}  // namespace detail

#if defined(CPPWRF_USE_CUDA)
// f(chunk, work) for chunks of rows of r whose columns' work arrays (bytes_per_col) fit the buffer (one chunk when
// it holds every column); work is the buffer of column_for(chunk, ...)'s threads
template <class F>
void column_work_chunks(const Range2D& r, std::size_t bytes_per_col, F f) {
  if (r.i1 <= r.i0 || r.j1 <= r.j0) return;
  const int ni = r.i1 - r.i0, nj = r.j1 - r.j0;
  const std::size_t bpc = detail::column_work_bpc(bytes_per_col > 0 ? bytes_per_col : 1);
  std::size_t have = 0;
  char* buf = detail::column_work_buffer(0, &have);
  // the most block rows that fit, then chunks of equal size (whole block rows: a block covers by rows)
  const int by = column_block().y;
  const std::size_t pair = detail::column_for_threads(ni, by) * bpc;
  std::size_t fit = have / pair;
  if (fit < 1) {
    buf = detail::column_work_buffer(pair, &have);
    fit = 1;
  }
  const int npairs = (nj + by - 1) / by;
  int pairs = fit >= std::size_t(npairs) ? npairs : int(fit);
  const int nchunks = (npairs + pairs - 1) / pairs;
  pairs = (npairs + nchunks - 1) / nchunks;
  const int rows = pairs * by;
  for (int j0 = r.j0; j0 < r.j1; j0 += rows) {
    const int j1 = j0 + rows < r.j1 ? j0 + rows : r.j1;
    f(Range2D{r.i0, r.i1, j0, j1}, ColumnWork{buf, bpc});
  }
}
#endif

// column_for(r, k) with k.work set (or k.set_work): on the GPU in chunks of rows that fit the buffer (one launch when
// it holds every column; every column is computed by itself, so the chunks do not change any value). K::kWorkBytes
// is the bytes per column.
template <class K>
void column_for_work(const Range2D& r, K k) {
#if defined(CPPWRF_USE_CUDA)
  column_work_chunks(r, K::kWorkBytes, [&](const Range2D& c, const ColumnWork& w) {
    if constexpr (requires { k.set_work(w); }) k.set_work(w);
    else k.work = w;
    column_for(c, k);
  });
#else
  column_for(r, k);
#endif
}

}  // namespace cppwrf
