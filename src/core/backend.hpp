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
#include <type_traits>
#include <typeinfo>
#include <stdexcept>
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
CPPWRF_HD inline float rlog(float a) {
#if defined(__CUDA_ARCH__)
  return logf(a);
#else
  return std::log(a);
#endif
}
CPPWRF_HD inline double rlog(double a) {
#if defined(__CUDA_ARCH__)
  return ::log(a);
#else
  return std::log(a);
#endif
}
CPPWRF_HD inline double rpow(double a, double b) {
#if defined(__CUDA_ARCH__)
  return ::pow(a, b);
#else
  return std::pow(a, b);
#endif
}

// ---- WRF's vectorised libm ----
// WRF (gfortran -O2 -ftree-vectorize on x86-64) computes exp/log/pow/sin/... inside vectorised loops
// with glibc's libmvec (2-lane SSE2 ABI for REAL*8): pairs of iterations from the start of the loop,
// a trailing odd iteration (and a 1-trip loop) with the scalar libm. The vector results are lane
// independent and differ from the scalar libm in the last bit for ~30% of arguments. On the host the
// v* functions below give WRF's value for iteration i (0-based) of an n-trip loop; double only (float
// and the device use the scalar functions).
#if defined(__x86_64__) && !defined(__CUDA_ARCH__) && !defined(CPPWRF_NO_LIBMVEC)
#include <emmintrin.h>
extern "C" __m128d _ZGVbN2vv_pow(__m128d, __m128d);
extern "C" __m128d _ZGVbN2v_log(__m128d);
extern "C" __m128d _ZGVbN2v_exp(__m128d);
extern "C" __m128d _ZGVbN2v_sin(__m128d);
extern "C" __m128d _ZGVbN2v_cos(__m128d);
#define CPPWRF_HAVE_LIBMVEC 1
#endif
namespace cppwrf {
// x**n for a variable INTEGER n as gfortran computes it: libgcc's __powidf2 (square-and-multiply,
// 1/y for n < 0); differs from pow(x, double(n)) in the last bit for some arguments.
template <typename T>
CPPWRF_HD inline T powi(T x, int m) {
  unsigned int n = m < 0 ? 0u - unsigned(m) : unsigned(m);
  T y = (n % 2) ? x : T(1);
  while (n >>= 1) {
    x = x * x;
    if (n % 2) y = y * x;
  }
  return m < 0 ? T(1) / y : y;
}
}  // namespace cppwrf
namespace cppwrf::vm {
CPPWRF_HD inline bool lane(int i, int n) { return n > 1 && !((n & 1) && i == n - 1); }
CPPWRF_HD inline double pow(double a, double b, int i, int n) {
#if defined(CPPWRF_HAVE_LIBMVEC)
  if (lane(i, n)) return _mm_cvtsd_f64(_ZGVbN2vv_pow(_mm_set1_pd(a), _mm_set1_pd(b)));
#endif
  (void)i; (void)n;
  return rpow(a, b);
}
CPPWRF_HD inline double log(double a, int i, int n) {
#if defined(CPPWRF_HAVE_LIBMVEC)
  if (lane(i, n)) return _mm_cvtsd_f64(_ZGVbN2v_log(_mm_set1_pd(a)));
#endif
  (void)i; (void)n;
  return rlog(a);
}
CPPWRF_HD inline double exp(double a, int i, int n) {
#if defined(CPPWRF_HAVE_LIBMVEC)
  if (lane(i, n)) return _mm_cvtsd_f64(_ZGVbN2v_exp(_mm_set1_pd(a)));
#endif
  (void)i; (void)n;
#if defined(__CUDA_ARCH__)
  return ::exp(a);
#else
  return std::exp(a);
#endif
}
CPPWRF_HD inline double sin(double a, int i, int n) {
#if defined(CPPWRF_HAVE_LIBMVEC)
  if (lane(i, n)) return _mm_cvtsd_f64(_ZGVbN2v_sin(_mm_set1_pd(a)));
#endif
  (void)i; (void)n;
#if defined(__CUDA_ARCH__)
  return ::sin(a);
#else
  return std::sin(a);
#endif
}
CPPWRF_HD inline double cos(double a, int i, int n) {
#if defined(CPPWRF_HAVE_LIBMVEC)
  if (lane(i, n)) return _mm_cvtsd_f64(_ZGVbN2v_cos(_mm_set1_pd(a)));
#endif
  (void)i; (void)n;
#if defined(__CUDA_ARCH__)
  return ::cos(a);
#else
  return std::cos(a);
#endif
}
CPPWRF_HD inline float pow(float a, float b, int, int) { return rpow(a, b); }
CPPWRF_HD inline float log(float a, int, int) { return rlog(a); }
CPPWRF_HD inline float exp(float a, int, int) {
#if defined(__CUDA_ARCH__)
  return expf(a);
#else
  return std::exp(a);
#endif
}
}  // namespace cppwrf::vm

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
// Developer kernel timer: CPPWRF_KTIME=1 synchronises around every parallel_for/column_for and
// prints the time per kernel type at exit (no performance counters needed). Off: one branch.
namespace ktime {
bool on();
void begin();
void end(const char* kernel_type);
}  // namespace ktime
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
  static_assert(__nv_is_extended_host_device_lambda_closure_type(K) || std::is_trivially_copyable_v<K>,
                "kernels are copied to the device: trivially copyable structs/lambdas only (no std::string, std::vector, ...)");
  if (r.i1 <= r.i0 || r.j1 <= r.j0 || r.k1 <= r.k0) return;
  dim3 b(64, 2, 1);
  dim3 g((r.i1 - r.i0 + 63) / 64, (r.j1 - r.j0 + 1) / 2, r.k1 - r.k0);
  const bool kt = ktime::on();
  if (kt) ktime::begin();
  detail::k3d<<<g, b>>>(k, r);
  cuda_check_launch("parallel_for");
  if (kt) ktime::end(typeid(K).name());
}
template <class K> void column_for(const Range2D& r, const K& k) {
  static_assert(__nv_is_extended_host_device_lambda_closure_type(K) || std::is_trivially_copyable_v<K>,
                "kernels are copied to the device: trivially copyable structs/lambdas only (no std::string, std::vector, ...)");
  if (r.i1 <= r.i0 || r.j1 <= r.j0) return;
  dim3 b(64, 2, 1);
  dim3 g((r.i1 - r.i0 + 63) / 64, (r.j1 - r.j0 + 1) / 2, 1);
  const bool kt = ktime::on();
  if (kt) ktime::begin();
  detail::k2d<<<g, b>>>(k, r);
  cuda_check_launch("column_for");
  if (kt) ktime::end(typeid(K).name());
}
#endif  // __CUDACC__
#else   // ---- CPU / OpenMP ----
template <class K> void parallel_for(const Range3D& r, const K& k) {
  static_assert(std::is_trivially_copyable_v<K>, "kernels are copied to the device: trivially copyable structs/lambdas only (no std::string, std::vector, ...)");
#pragma omp parallel for collapse(2) schedule(static)
  for (int kk = r.k0; kk < r.k1; ++kk)
    for (int j = r.j0; j < r.j1; ++j)
#pragma omp simd
      for (int i = r.i0; i < r.i1; ++i) k(i, j, kk);
}
template <class K> void column_for(const Range2D& r, const K& k) {
  static_assert(std::is_trivially_copyable_v<K>, "kernels are copied to the device: trivially copyable structs/lambdas only (no std::string, std::vector, ...)");
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
    if (nk > max_k_) throw std::runtime_error("Reducer: k range exceeds max_k");
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
