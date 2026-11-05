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
};

// a [n1][N2] work array: element (k, n) of a thread at p[32 (k N2 + n)]
template <class T, int N2>
struct ColArray2 {
  T* p;
  CPPWRF_FORCEINLINE ColArray<T> operator[](int k) const { return ColArray<T>{p + k * N2 * 32}; }
};

// One launch's buffer: warp w of the launch (thread t = w * 32 + lane in column_for's numbering: blocks of
// kColumnBlockX x kColumnBlockY threads, x fastest) owns 32 * bpc bytes at base + w * 32 * bpc; in it the arrays
// one after another, each n * sizeof(T) * 32 bytes
inline constexpr int kColumnBlockX = 64, kColumnBlockY = 2;  // column_for's block (core/backend.hpp)
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
constexpr std::size_t column_for_threads(int ni, int nj) {
  return std::size_t((ni + kColumnBlockX - 1) / kColumnBlockX * kColumnBlockX) *
         std::size_t((nj + kColumnBlockY - 1) / kColumnBlockY * kColumnBlockY);
}
}  // namespace detail

// column_for(r, k) with k.work set (or k.set_work): on the GPU in chunks of rows that fit the buffer (one launch when
// it holds every column; every column is computed by itself, so the chunks do not change any value). K::kWorkBytes
// is the bytes per column.
template <class K>
void column_for_work(const Range2D& r, K k) {
#if defined(CPPWRF_USE_CUDA)
  if (r.i1 <= r.i0 || r.j1 <= r.j0) return;
  const int ni = r.i1 - r.i0, nj = r.j1 - r.j0;
  const std::size_t bpc = detail::column_work_bpc(K::kWorkBytes > 0 ? K::kWorkBytes : 1);
  std::size_t have = 0;
  char* buf = detail::column_work_buffer(0, &have);
  // the most pairs of rows that fit, then chunks of equal size (whole pairs: a block covers two rows)
  const std::size_t pair = detail::column_for_threads(ni, kColumnBlockY) * bpc;
  std::size_t fit = have / pair;
  if (fit < 1) {
    buf = detail::column_work_buffer(pair, &have);
    fit = 1;
  }
  const int npairs = (nj + kColumnBlockY - 1) / kColumnBlockY;
  int pairs = fit >= std::size_t(npairs) ? npairs : int(fit);
  const int nchunks = (npairs + pairs - 1) / pairs;
  pairs = (npairs + nchunks - 1) / nchunks;
  const int rows = pairs * kColumnBlockY;
  for (int j0 = r.j0; j0 < r.j1; j0 += rows) {
    const int j1 = j0 + rows < r.j1 ? j0 + rows : r.j1;
    const ColumnWork w{buf, bpc};
    if constexpr (requires { k.set_work(w); }) k.set_work(w);
    else k.work = w;
    column_for(Range2D{r.i0, r.i1, j0, j1}, k);
  }
#else
  column_for(r, k);
#endif
}

}  // namespace cppwrf
