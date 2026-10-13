#pragma once
#include <cmath>
// Minimal single-source execution backend (CPU/OpenMP or CUDA).
//   parallel_for(Range3D, k)  : k(i,j,k) for every cell, cells independent.
//   column_for(Range2D, k)    : k(i,j) once per column (vertical sweeps in registers).
//   Reducer::max/sum(Range3D, f): f(i,j,k)->double; deterministic (per-k partials,
//                                 combined on the host in fixed order).
// Ranges are half-open. Kernels are trivially-copyable structs, no virtuals/exceptions.
#include <cstddef>
#include <cstdint>
#include <vector>

#if defined(__CUDACC__)
#define CPPWRF_HD __host__ __device__
#define CPPWRF_DEV __device__
#else
#define CPPWRF_HD
#define CPPWRF_DEV
#endif

// pow in the working precision on host and device (powf for float, as Fortran REAL ** REAL);
// unqualified pow(float, float) would resolve to the double C function on the host.
CPPWRF_HD inline float rpow(float a, float b) {
#if defined(__CUDA_ARCH__)
  return powf(a, b);
#else
  return std::pow(a, b);
#endif
}
CPPWRF_HD inline double rpow(double a, double b) {
#if defined(__CUDA_ARCH__)
  return ::pow(a, b);
#else
  return std::pow(a, b);
#endif
}

namespace cppwrf {

struct Range3D { int i0, i1, j0, j1, k0, k1; };
struct Range2D { int i0, i1, j0, j1; };

// ---- memory (implemented in backend.cpp) ----
void* dev_alloc(std::size_t bytes);   // zero-initialised
void dev_free(void* p);
void copy_to_host(void* dst, const void* src, std::size_t bytes);
void copy_to_dev(void* dst, const void* src, std::size_t bytes);
void copy_dev(void* dst, const void* src, std::size_t bytes);
void dev_sync();
const char* backend_name();

#if defined(CPPWRF_USE_CUDA)
void cuda_check_launch(const char* where);
#if defined(__CUDACC__)
namespace detail {
template <class K> __global__ void k3d(K k, Range3D r) {
  int i = r.i0 + blockIdx.x * blockDim.x + threadIdx.x;
  int j = r.j0 + blockIdx.y * blockDim.y + threadIdx.y;
  int kk = r.k0 + blockIdx.z;
  if (i < r.i1 && j < r.j1) k(i, j, kk);
}
template <class K> __global__ void k2d(K k, Range2D r) {
  int i = r.i0 + blockIdx.x * blockDim.x + threadIdx.x;
  int j = r.j0 + blockIdx.y * blockDim.y + threadIdx.y;
  if (i < r.i1 && j < r.j1) k(i, j);
}
// One block per k plane; fixed-order tree => deterministic.
template <class F, bool MAX> __global__ void kred(F f, Range3D r, double* part) {
  __shared__ double sh[256];
  int kk = r.k0 + blockIdx.x;
  double acc = MAX ? -1e308 : 0.0;
  int ni = r.i1 - r.i0, n = ni * (r.j1 - r.j0);
  for (int t = threadIdx.x; t < n; t += 256) {
    double v = f(r.i0 + t % ni, r.j0 + t / ni, kk);
    acc = MAX ? (v > acc ? v : acc) : acc + v;
  }
  sh[threadIdx.x] = acc;
  __syncthreads();
  for (int s = 128; s > 0; s >>= 1) {
    if (threadIdx.x < s) {
      double a = sh[threadIdx.x], b = sh[threadIdx.x + s];
      sh[threadIdx.x] = MAX ? (a > b ? a : b) : a + b;
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) part[blockIdx.x] = sh[0];
}
}  // namespace detail

template <class K> void parallel_for(const Range3D& r, const K& k) {
  if (r.i1 <= r.i0 || r.j1 <= r.j0 || r.k1 <= r.k0) return;
  dim3 b(64, 2, 1);
  dim3 g((r.i1 - r.i0 + 63) / 64, (r.j1 - r.j0 + 1) / 2, r.k1 - r.k0);
  detail::k3d<<<g, b>>>(k, r);
  cuda_check_launch("parallel_for");
}
template <class K> void column_for(const Range2D& r, const K& k) {
  if (r.i1 <= r.i0 || r.j1 <= r.j0) return;
  dim3 b(64, 2, 1);
  dim3 g((r.i1 - r.i0 + 63) / 64, (r.j1 - r.j0 + 1) / 2, 1);
  detail::k2d<<<g, b>>>(k, r);
  cuda_check_launch("column_for");
}
#endif  // __CUDACC__
#else   // ---- CPU / OpenMP ----
template <class K> void parallel_for(const Range3D& r, const K& k) {
#pragma omp parallel for collapse(2) schedule(static)
  for (int kk = r.k0; kk < r.k1; ++kk)
    for (int j = r.j0; j < r.j1; ++j)
#pragma omp simd
      for (int i = r.i0; i < r.i1; ++i) k(i, j, kk);
}
template <class K> void column_for(const Range2D& r, const K& k) {
#pragma omp parallel for collapse(2) schedule(static)
  for (int j = r.j0; j < r.j1; ++j)
    for (int i = r.i0; i < r.i1; ++i) k(i, j);
}
#endif

// Deterministic reductions (diagnostics only; not in the hot path).
class Reducer {
 public:
  explicit Reducer(int max_k);
  ~Reducer();
  Reducer(const Reducer&) = delete;
  Reducer& operator=(const Reducer&) = delete;
  template <class F> double max(const Range3D& r, const F& f) { return run<F, true>(r, f); }
  template <class F> double sum(const Range3D& r, const F& f) { return run<F, false>(r, f); }

 private:
  int max_k_;
  double* dpart_ = nullptr;  // device partials (CUDA) or unused
  std::vector<double> hpart_;
  template <class F, bool MAX> double run(const Range3D& r, const F& f) {
    int nk = r.k1 - r.k0;
    if (nk <= 0) return MAX ? -1e308 : 0.0;
#if defined(CPPWRF_USE_CUDA)
#if defined(__CUDACC__)
    detail::kred<F, MAX><<<nk, 256>>>(f, r, dpart_);
    cuda_check_launch("Reducer");
    copy_to_host(hpart_.data(), dpart_, sizeof(double) * nk);
#endif
#else
#pragma omp parallel for schedule(static)
    for (int kk = r.k0; kk < r.k1; ++kk) {
      double acc = MAX ? -1e308 : 0.0;
      for (int j = r.j0; j < r.j1; ++j)
        for (int i = r.i0; i < r.i1; ++i) {
          double v = f(i, j, kk);
          acc = MAX ? (v > acc ? v : acc) : acc + v;
        }
      hpart_[kk - r.k0] = acc;
    }
#endif
    double acc = MAX ? -1e308 : 0.0;
    for (int t = 0; t < nk; ++t) acc = MAX ? (hpart_[t] > acc ? hpart_[t] : acc) : acc + hpart_[t];
    return acc;
  }
};

}  // namespace cppwrf
