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
#define CPPWRF_FORCEINLINE __host__ __device__ __forceinline__
#else
#define CPPWRF_HD
#define CPPWRF_DEV
#define CPPWRF_FORCEINLINE inline
#endif

#if defined(__GNUC__) || defined(__clang__) || defined(__CUDACC__)
#define CPPWRF_RESTRICT __restrict__
#else
#define CPPWRF_RESTRICT
#endif

// CUDA FP64 build: glibc's libm and libmvec algorithms on the device and the host (core/glibc_math.hpp), the
// results of WRF's real*8 build and of the CPU build; CMake option CPPWRF_CUDA_GLIBC_MATH (default ON; OFF:
// CPPWRF_CUDA_LIBM, CUDA's math library). The kernels call the math library as cppwrf::libm (core/libm.hpp).
#if defined(CPPWRF_USE_CUDA) && defined(CPPWRF_DOUBLE) && !defined(CPPWRF_CUDA_LIBM)
#define CPPWRF_GLIBC_MATH 1
#endif
#include "core/libm.hpp"

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
#if defined(CPPWRF_GLIBC_MATH)
  return cppwrf::libm::log(a);
#elif defined(__CUDA_ARCH__)
  return ::log(a);
#else
  return std::log(a);
#endif
}
CPPWRF_HD inline double rpow(double a, double b) {
#if defined(CPPWRF_GLIBC_MATH)
  return cppwrf::libm::pow_d(a, b);
#elif defined(__CUDA_ARCH__)
  return ::pow(a, b);
#else
  return std::pow(a, b);
#endif
}
// x**2.0 with the constant exponent 2.0: gcc folds pow(x, 2.0) into x * x (WRF's gfortran build and the CPU build);
// the device cannot test for a constant, so the CUDA FP64 build writes the product out
template <typename T>
CPPWRF_HD inline T rpow2c(T a) {
#if defined(CPPWRF_GLIBC_MATH)
  if constexpr (sizeof(T) == sizeof(double)) return a * a;
#endif
  return rpow(a, T(2.0));
}
CPPWRF_HD inline float rexp(float a) {
#if defined(__CUDA_ARCH__)
  return expf(a);
#else
  return std::exp(a);
#endif
}
CPPWRF_HD inline double rexp(double a) {
#if defined(CPPWRF_GLIBC_MATH)
  return cppwrf::libm::exp(a);
#elif defined(__CUDA_ARCH__)
  return ::exp(a);
#else
  return std::exp(a);
#endif
}

// Fortran MIN / MAX of two values (b on ties and when either is NaN)
template <class T>
CPPWRF_HD inline T c_min(T a, T b) { return (a < b) ? a : b; }
template <class T>
CPPWRF_HD inline T c_max(T a, T b) { return (a > b) ? a : b; }

// ---- WRF's vectorised libm ----
// WRF (gfortran -O2 -ftree-vectorize on x86-64) computes exp/log/pow/sin/... inside vectorised loops
// with glibc's libmvec (2-lane SSE2 ABI for REAL*8): pairs of iterations from the start of the loop,
// a trailing odd iteration (and a 1-trip loop) with the scalar libm. The vector results are lane
// independent and differ from the scalar libm in the last bit for ~30% of arguments. On the host the
// v* functions below give WRF's value for iteration i (0-based) of an n-trip loop; double only (float
// and the device use the scalar functions).
#if defined(__x86_64__) && !defined(__CUDA_ARCH__) && !defined(CPPWRF_NO_LIBMVEC) && !defined(CPPWRF_GLIBC_MATH)
#include <emmintrin.h>
extern "C" __m128d _ZGVbN2vv_pow(__m128d, __m128d);
extern "C" __m128d _ZGVbN2v_log(__m128d);
extern "C" __m128d _ZGVbN2v_exp(__m128d);
extern "C" __m128d _ZGVbN2v_sin(__m128d);
extern "C" __m128d _ZGVbN2v_atan(__m128d);
extern "C" __m128d _ZGVbN2v_cos(__m128d);
extern "C" __m128d _ZGVbN2v_log10(__m128d);
#define CPPWRF_HAVE_LIBMVEC 1
#endif
// log and log10: glibc's SSE4 variants with Intel's rcpps (core/libmvec_x86.hpp), the same on any x86-64 CPU
#include "core/libmvec_x86.hpp"
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
#if defined(CPPWRF_GLIBC_MATH)
// CUDA FP64: libmvec's lanes and the scalar libm as glibc computes them, on the device and the host
CPPWRF_HD inline double pow(double a, double b, int i, int n) {
  return lane(i, n) ? glibc::mvec::pow(a, b) : cppwrf::libm::pow_d(a, b);
}
CPPWRF_HD inline double log(double a, int i, int n) {
  return lane(i, n) ? mvec::log(a) : cppwrf::libm::log(a);
}
CPPWRF_HD inline double exp(double a, int i, int n) {
  return lane(i, n) ? glibc::mvec::exp(a) : cppwrf::libm::exp(a);
}
CPPWRF_HD inline double atan(double a, int i, int n) {
  return lane(i, n) ? glibc::mvec::atan(a) : cppwrf::libm::atan(a);
}
CPPWRF_HD inline double sin(double a, int i, int n) {
  return lane(i, n) ? glibc::mvec::sin(a) : cppwrf::libm::sin(a);
}
CPPWRF_HD inline double cos(double a, int i, int n) {
  return lane(i, n) ? glibc::mvec::cos(a) : cppwrf::libm::cos(a);
}
CPPWRF_HD inline double log10(double a, int i, int n) {
  return lane(i, n) ? mvec::log10(a) : cppwrf::libm::log10(a);
}
#else
CPPWRF_HD inline double pow(double a, double b, int i, int n) {
#if defined(CPPWRF_HAVE_LIBMVEC)
  if (lane(i, n)) return _mm_cvtsd_f64(_ZGVbN2vv_pow(_mm_set1_pd(a), _mm_set1_pd(b)));
#endif
  (void)i; (void)n;
  return rpow(a, b);
}
CPPWRF_HD inline double log(double a, int i, int n) {
#if defined(CPPWRF_HAVE_LIBMVEC) && defined(CPPWRF_PORTABLE_LOG)
  if (lane(i, n)) return mvec::log(a);
#elif defined(CPPWRF_HAVE_LIBMVEC)
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
CPPWRF_HD inline double atan(double a, int i, int n) {
#if defined(CPPWRF_HAVE_LIBMVEC)
  if (lane(i, n)) return _mm_cvtsd_f64(_ZGVbN2v_atan(_mm_set1_pd(a)));
#endif
  (void)i; (void)n;
#if defined(__CUDA_ARCH__)
  return ::atan(a);
#else
  return std::atan(a);
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
CPPWRF_HD inline double log10(double a, int i, int n) {
#if defined(CPPWRF_HAVE_LIBMVEC) && defined(CPPWRF_PORTABLE_LOG)
  if (lane(i, n)) return mvec::log10(a);
#elif defined(CPPWRF_HAVE_LIBMVEC)
  if (lane(i, n)) return _mm_cvtsd_f64(_ZGVbN2v_log10(_mm_set1_pd(a)));
#endif
  (void)i; (void)n;
#if defined(__CUDA_ARCH__)
  return ::log10(a);
#else
  return std::log10(a);
#endif
}
#endif
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

// One model time step (dyn::RK3::step). CUDA: the step's launches and device copies are captured into CUDA graphs
// and replayed (an instantiated graph per domain and segment, updated with the new launch parameters), so the
// device runs the same kernels in the same order with the same arguments as without graphs. An operation that
// waits for the device (copy_to_host, copy_to_dev, dev_alloc, dev_free, dev_sync, Reducer) ends the current graph
// segment, launches it, and the capture resumes after it. Inside the step fatal_check() only copies the device
// flag asynchronously; the flags are read once, when the host next waits for the step (end() waits for the
// previous step before it launches this one). CPPWRF_GRAPHS=0 (or CPPWRF_KTIME) launches without graphs;
// CPPWRF_GRAPHS=stats prints the graph counts at exit. CPU build: no-op.
class StepScope {
 public:
  explicit StepScope(const void* key);
  ~StepScope();
  StepScope(const StepScope&) = delete;
  StepScope& operator=(const StepScope&) = delete;
  void end();  // launches the step's last segment (does not wait for it)

 private:
  bool open_ = false;
};
bool step_graphs_on();  // the steps run as CUDA graphs (the caller need not wait for each step)

#if defined(CPPWRF_USE_CUDA)
void cuda_check_launch(const char* where);
void* cuda_stream_raw();  // the stream of every launch and copy (created on first use; blocking w.r.t. stream 0)
namespace detail {
int* fatal_poll_slot(void (*reset)());  // inside a step: pinned host slot for this check's flag; else nullptr
void fatal_flush();                     // reads the flags of the launched, unchecked checks (waits for them)
}  // namespace detail
// Developer kernel timer: CPPWRF_KTIME=1 synchronises around every parallel_for/column_for and
// prints the time per kernel type at exit (no performance counters needed). Off: one branch.
namespace ktime {
bool on();
void begin();
void end(const char* kernel_type);
}  // namespace ktime
#if defined(__CUDACC__)
inline cudaStream_t cuda_launch_stream() {
  return static_cast<cudaStream_t>(cuda_stream_raw());
}
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
  const unsigned int by = (r.j1 - r.j0 == 1) ? 1u : 2u;
  dim3 b(64, by, 1);
  dim3 g((r.i1 - r.i0 + 63) / 64, (r.j1 - r.j0 + by - 1) / by, r.k1 - r.k0);
  const bool kt = ktime::on();
  if (kt) ktime::begin();
  detail::k3d<<<g, b, 0, cuda_launch_stream()>>>(k, r);
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
  detail::k2d<<<g, b, 0, cuda_launch_stream()>>>(k, r);
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

// ---- launch fusion: independent kernels in one launch ----
// LaunchBatch<K, N> collects up to N kernels of the type K, each with its own range, and launches them as one CUDA
// kernel (a full batch launches by itself, run() launches the rest); parallel_fused(part(r1, k1), part(r2, k2), ...)
// does the same for kernels of different types. Every part computes what its own parallel_for computes: part m
// takes the blocks z0[m] <= blockIdx.z < z0[m + 1] (its k range) with parallel_for's thread mapping in x and y.
// The parts must be independent (no part reads or writes what another part writes): on the GPU their cells run
// in any order. The CPU build runs each part at once as its own parallel_for, in the order given.
// CUDA FP32 (nvcc's default --fmad=true): ptxas contracts a*b + c into an fma depending on the code around it, so
// a kernel with such arithmetic can round differently inside a batched kernel (the FP64 build compiles with
// --fmad=false: every operation is rounded on its own). There only kernels that declare
// `static constexpr bool kNoMulAdd = true` (copies, zeros, comparisons) are batched; the others are launched one by
// one with parallel_for, the kernel they had before.
#if !defined(CPPWRF_USE_CUDA) || defined(__CUDACC__)
#if defined(CPPWRF_USE_CUDA) && !defined(CPPWRF_DOUBLE)
inline constexpr bool kFuseMulAdd = false;
#else
inline constexpr bool kFuseMulAdd = true;
#endif
template <class K> concept Fusable = kFuseMulAdd || requires { requires K::kNoMulAdd; };
template <class K> struct Part {
  Range3D r;
  K k;
};
template <class K> Part<K> part(const Range3D& r, const K& k) { return Part<K>{r, k}; }
namespace detail {
inline bool range_empty(const Range3D& r) { return r.i1 <= r.i0 || r.j1 <= r.j0 || r.k1 <= r.k0; }
#if defined(CPPWRF_USE_CUDA)
struct LaunchShape {  // the grid over every part: x and y as the largest part, z = the parts' k ranges in turn
  unsigned int gx = 0, nj = 0, nz = 0;
  void add(const Range3D& r) {
    if (range_empty(r)) return;
    const unsigned int x = unsigned(r.i1 - r.i0 + 63) / 64u, y = unsigned(r.j1 - r.j0);
    gx = x > gx ? x : gx;
    nj = y > nj ? y : nj;
    nz += unsigned(r.k1 - r.k0);
  }
  unsigned int by() const { return nj == 1 ? 1u : 2u; }
};
template <class K, int N> struct BatchArgs {
  K k[N];
  Range3D r[N];
  int z0[N + 1];
  int n;
};
// compile-time part indices only: a run-time index into the kernel parameters would copy them to local memory
template <int M, class K, int N> __device__ __forceinline__ void batch_call(const BatchArgs<K, N>& a, int z) {
  if (M + 1 == N || M + 1 >= a.n || z < a.z0[M + 1]) {
    const Range3D& r = a.r[M];
    const int i = r.i0 + blockIdx.x * blockDim.x + threadIdx.x;
    const int j = r.j0 + blockIdx.y * blockDim.y + threadIdx.y;
    if (i < r.i1 && j < r.j1) a.k[M](i, j, r.k0 + (z - a.z0[M]));
  } else if constexpr (M + 1 < N) {
    batch_call<M + 1>(a, z);
  }
}
template <class K, int N> __global__ void kbatch(const BatchArgs<K, N> a) { batch_call<0>(a, int(blockIdx.z)); }
template <class... K> struct FusedArgs {};
template <class K0, class... R> struct FusedArgs<K0, R...> {
  K0 k;
  Range3D r;
  int z0, z1;
  FusedArgs<R...> rest;
};
template <class K0, class... R> __device__ __forceinline__ void fused_call(const FusedArgs<K0, R...>& a, int z) {
  if (z >= a.z0 && z < a.z1) {
    const int i = a.r.i0 + blockIdx.x * blockDim.x + threadIdx.x;
    const int j = a.r.j0 + blockIdx.y * blockDim.y + threadIdx.y;
    if (i < a.r.i1 && j < a.r.j1) a.k(i, j, a.r.k0 + (z - a.z0));
  } else if constexpr (sizeof...(R) > 0) {
    fused_call(a.rest, z);
  }
}
template <class... K> __global__ void kfused(const FusedArgs<K...> a) { fused_call(a, int(blockIdx.z)); }
inline void fused_fill(FusedArgs<>&, LaunchShape&) {}
template <class K0, class... R>
void fused_fill(FusedArgs<K0, R...>& a, LaunchShape& s, const Part<K0>& p, const Part<R>&... rest) {
  a.k = p.k;
  a.r = p.r;
  a.z0 = int(s.nz);
  s.add(p.r);
  a.z1 = int(s.nz);
  fused_fill(a.rest, s, rest...);
}
#endif
}  // namespace detail

template <class K, int N = 8> class LaunchBatch {
  static_assert(std::is_trivially_copyable_v<K>, "kernels are copied to the device: trivially copyable structs only");

 public:
  LaunchBatch() = default;
  LaunchBatch(const LaunchBatch&) = delete;
  LaunchBatch& operator=(const LaunchBatch&) = delete;
#if defined(CPPWRF_USE_CUDA)
  static_assert(sizeof(detail::BatchArgs<K, N>) <= 32000, "kernel parameters too large: use a smaller batch");
  ~LaunchBatch() { run(); }
  void add(const Range3D& r, const K& k) {
    if constexpr (!Fusable<K>) {
      parallel_for(r, k);
      return;
    }
    if (detail::range_empty(r)) return;
    a_.k[a_.n] = k;
    a_.r[a_.n] = r;
    if (++a_.n == N) run();
  }
  void run() {
    if (a_.n == 0) return;
    detail::LaunchShape s;
    for (int m = 0; m < a_.n; ++m) {
      a_.z0[m] = int(s.nz);
      s.add(a_.r[m]);
    }
    a_.z0[a_.n] = int(s.nz);
    const unsigned int by = s.by();
    const dim3 b(64, by, 1), g(s.gx, (s.nj + by - 1) / by, s.nz);
    const bool kt = ktime::on();
    if (kt) ktime::begin();
    detail::kbatch<K, N><<<g, b, 0, cuda_launch_stream()>>>(a_);
    a_.n = 0;
    cuda_check_launch("LaunchBatch");
    if (kt) ktime::end(typeid(K).name());
  }

 private:
  detail::BatchArgs<K, N> a_{};
#else
  void add(const Range3D& r, const K& k) { parallel_for(r, k); }
  void run() {}
#endif
};

template <class... K> void parallel_fused(const Part<K>&... parts) {
#if defined(CPPWRF_USE_CUDA)
  if constexpr (!(Fusable<K> && ...)) {
    (parallel_for(parts.r, parts.k), ...);
    return;
  }
  using Args = detail::FusedArgs<K...>;
  static_assert(sizeof(Args) <= 32000, "kernel parameters too large");
  Args a{};
  detail::LaunchShape s;
  detail::fused_fill(a, s, parts...);
  if (s.nz == 0) return;
  const unsigned int by = s.by();
  const dim3 b(64, by, 1), g(s.gx, (s.nj + by - 1) / by, s.nz);
  const bool kt = ktime::on();
  if (kt) ktime::begin();
  detail::kfused<K...><<<g, b, 0, cuda_launch_stream()>>>(a);
  cuda_check_launch("parallel_fused");
  if (kt) ktime::end(typeid(Args).name());
#else
  (parallel_for(parts.r, parts.k), ...);
#endif
}
#endif

// WRF's wrf_error_fatal inside a physics kernel. A kernel cannot stop the run: it records the
// condition (raise_fatal; the first code wins, atomic on host and device) and the host calls
// fatal_check() right after the physics call, which prints WRF's FATAL block and throws FatalError.
// The device flag is local to each CUDA translation unit: call fatal_check() in the file that
// launches the kernel.
enum FatalCode : int {
  kFatalNone = 0,
  kFatalSnowDensity = 1,  // Noah land / sea ice / glacial ice: SNEQV / SNOWH > 1
  kFatalKfQgSurface = 2,  // Kain-Fritsch: QG(1) < 0 after the feedback
  kFatalKfIstop = 3,      // Kain-Fritsch: mass does not balance at LTOP (ISTOP = 1)
  kFatalNoahSoilTypes = 4,     // Noah REDPRM: SOILTYP > SLCATS
  kFatalNoahLanduseTypes = 5,  // VEGTYP > LUCATS
  kFatalNoahSlopeTypes = 6,    // SLOPETYP > SLPCATS
  kFatalNoahRootLayers = 7,    // NROOT > NSOIL
  kFatalNoahmpWater = 8,       // Noah-MP (and glacier) water balance error > 0.1 kg m-2
  kFatalNoahmpSolar = 9,       // shortwave balance error > 0.01 W m-2
  kFatalNoahmpEnergy = 10,     // land energy balance error > 0.01 W m-2
  kFatalNoahmpLongwave = 11,   // emitted longwave <= 0
  kFatalNoahmpZeroPlane = 12,  // canopy height <= zero-plane displacement
  kFatalNoahmpRefHeightCanopy = 13,  // reference height <= zero-plane displacement (above canopy)
  kFatalNoahmpRefHeightBare = 14,    // reference height <= zero-plane displacement (bare ground)
  kFatalNoahmpGlacierEnergy = 15,    // glacier surface energy balance error > 0.01 W m-2
  kFatalNoahUrbanLai = 16,           // Noah rdlai2d: urban point with VEGFRA = 0 and LAI <= 0 (not a WRF stop)
};
struct FatalError : std::runtime_error {
  int code;
  FatalError(int c, const char* msg) : std::runtime_error(msg), code(c) {}
};
extern int fatal_code_host;
const char* fatal_message(int code);
[[noreturn]] void fatal_stop(int code);  // prints WRF's FATAL block, throws FatalError
#if defined(__CUDACC__)
namespace detail {
static __device__ int fatal_code_dev;
}
#endif
CPPWRF_HD inline void raise_fatal(int code) {
#if defined(__CUDA_ARCH__)
  atomicCAS(&detail::fatal_code_dev, 0, code);
#else
  int expected = 0;
  __atomic_compare_exchange_n(&fatal_code_host, &expected, code, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED);
#endif
}
#if defined(CPPWRF_USE_CUDA) && defined(__CUDACC__)
namespace detail {
static inline void fatal_reset_dev() {  // clears this translation unit's device flag after a stop
  const int zero = 0;
  cudaMemcpyToSymbol(fatal_code_dev, &zero, sizeof(int));
}
}  // namespace detail
#endif
static inline void fatal_check() {
#if defined(CPPWRF_USE_CUDA) && defined(__CUDACC__)
  // inside a time step (StepScope): the flag is copied to a pinned slot on the stream and read at the step's end
  if (int* slot = detail::fatal_poll_slot(&detail::fatal_reset_dev)) {
    cudaMemcpyFromSymbolAsync(slot, detail::fatal_code_dev, sizeof(int), 0, cudaMemcpyDeviceToHost,
                              cuda_launch_stream());
    return;
  }
  detail::fatal_flush();
#endif
  int code = __atomic_exchange_n(&fatal_code_host, 0, __ATOMIC_RELAXED);
#if defined(CPPWRF_USE_CUDA) && defined(__CUDACC__)
  int dcode = 0;
  const int zero = 0;
  cudaStream_t s = cuda_launch_stream();
  cudaMemcpyFromSymbolAsync(&dcode, detail::fatal_code_dev, sizeof(int), 0, cudaMemcpyDeviceToHost, s);
  cudaStreamSynchronize(s);
  if (dcode != 0) {
    cudaMemcpyToSymbolAsync(detail::fatal_code_dev, &zero, sizeof(int), 0, cudaMemcpyHostToDevice, s);
    cudaStreamSynchronize(s);
    if (code == 0) code = dcode;
  }
#endif
  if (code != 0) fatal_stop(code);
}

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
    detail::kred<F, MAX><<<nk, 256, 0, cuda_launch_stream()>>>(f, r, dpart_);
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
