#pragma once
// glibc 2.41 libmvec's 2-lane vector log and log10 (_ZGVbN2v_log, _ZGVbN2v_log10; the SSE4 variants that glibc
// selects on every CPU with SSE4.1) re-implemented statement by statement, one lane, from
// sysdeps/x86_64/fpu/multiarch/svml_d_log2_core_sse4.S and svml_d_log102_core_sse4.S. Both take an approximate
// reciprocal of the mantissa with `rcpps`, whose result the hardware defines (Intel and AMD differ); here it is
// Intel's, from the table tools/libmvec_tables.py dumped (src/core/libmvec_tables.hpp), so that the results equal
// libmvec's on an Intel CPU (where WRF's references are made) on any x86-64 CPU. Every other operation is an
// IEEE double operation in the assembly's order (the build has -ffp-contract=off: no FMA contraction).
// Lanes outside the normal range (x < DBL_MIN incl. 0 and negatives, x > DBL_MAX, NaN) call the scalar libm, as
// the assembly does. CPPWRF_PORTABLE_LOG: cppwrf::vm::log / log10 and cppwrf_real use these instead of libmvec's
// log / log10 (CMake option CPPWRF_PORTABLE_LIBMVEC, default ON; OFF defines CPPWRF_LIBMVEC_HW_LOG).
// Host and device functions (the device reads copies of the tables; the CUDA FP64 build computes the lanes of WRF's
// vectorised log / log10 with them, core/backend.hpp CPPWRF_GLIBC_MATH, and the scalar fallback with glibc's log /
// log10 of core/glibc_math.hpp).
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>

#include "core/libmvec_tables.hpp"
#if defined(CPPWRF_GLIBC_MATH)
#include "core/glibc_math.hpp"
#endif

#if defined(__x86_64__) && !defined(CPPWRF_NO_LIBMVEC) && !defined(CPPWRF_LIBMVEC_HW_LOG)
#define CPPWRF_PORTABLE_LOG 1
#endif

#if defined(__CUDACC__)
#define CPPWRF_MV __host__ __device__ __forceinline__
#else
#define CPPWRF_MV inline
#endif

namespace cppwrf::mvec {
namespace detail {
// k(t) = roundeven(rcpps(f)) on Intel, t = 0 ... 2048 (f = 2^-10 * (1 + t/2048 + ...), t = 2048: f = 2^-9)
inline constexpr std::array<double, 2049> kRcp = [] {
  std::array<double, 2049> a{};
  int n = 0;
  for (int t = 0; t <= 2048; ++t) {
    while (n < 512 && mvec_tables::kRcpSteps[n] <= t) ++n;
    a[std::size_t(t)] = double(1024 - n);
  }
  return a;
}();
#if defined(__CUDACC__)
// device copies of the tables (constant-initialised from the host ones)
template <class T, int N>
struct DevTable {
  T v[N];
};
template <class T, std::size_t N>
constexpr DevTable<T, int(N)> dev_table(const std::array<T, N>& a) {
  DevTable<T, int(N)> t{};
  for (std::size_t i = 0; i < N; ++i) t.v[i] = a[i];
  return t;
}
template <class T, int N>
constexpr DevTable<T, N> dev_table(const T (&a)[N]) {
  DevTable<T, N> t{};
  for (int i = 0; i < N; ++i) t.v[i] = a[i];
  return t;
}
static __device__ const DevTable<double, 2049> kRcpDev = dev_table(kRcp);
static __device__ const DevTable<std::uint64_t, 513> kLogLADev = dev_table(mvec_tables::kLogLA);
static __device__ const DevTable<std::uint64_t, 513> kLog10LADev = dev_table(mvec_tables::kLog10LA);
#endif
#if defined(__CUDA_ARCH__)
#define CPPWRF_MVT(host, dev) ::cppwrf::mvec::detail::dev.v
#else
#define CPPWRF_MVT(host, dev) host
#endif
// bit casts (memcpy: std::bit_cast is not available to every nvcc host pass)
CPPWRF_MV double d(std::uint64_t bits) {
#if defined(__CUDA_ARCH__)
  return __longlong_as_double(static_cast<long long>(bits));
#else
  double v;
  std::memcpy(&v, &bits, sizeof v);
  return v;
#endif
}
CPPWRF_MV std::uint64_t bits(double v) {
#if defined(__CUDA_ARCH__)
  return static_cast<std::uint64_t>(__double_as_longlong(v));
#else
  std::uint64_t b;
  std::memcpy(&b, &v, sizeof b);
  return b;
#endif
}
// the scalar log / log10 of the lanes outside the normal range
CPPWRF_MV double scalar_log(double x) {
#if defined(CPPWRF_GLIBC_MATH)
  return ::cppwrf::glibc::log(x);
#else
  return std::log(x);
#endif
}
CPPWRF_MV double scalar_log10(double x) {
#if defined(CPPWRF_GLIBC_MATH)
  return ::cppwrf::glibc::log10(x);
#else
  return std::log10(x);
#endif
}
// cvtpd2ps m; rcpps; cvtps2pd; roundpd $0 (log) or + / - 1.5*2^52 (log10): m in [2^-10, 2^-9)
CPPWRF_MV double rcp_rounded(double m) {
  const float f = static_cast<float>(m);
#if defined(__CUDA_ARCH__)
  const std::uint32_t fb = static_cast<std::uint32_t>(__float_as_uint(f));
#else
  std::uint32_t fb;
  std::memcpy(&fb, &f, sizeof fb);
#endif
  return CPPWRF_MVT(kRcp, kRcpDev)[(fb >> 12) - (117u << 11)];
}
}  // namespace detail

// _ZGVbN2v_log_sse4 (svml_d_log2_core_sse4.S, data svml_d_log_data.S)
CPPWRF_MV double log(double x) {
  using namespace mvec_tables;
  using detail::d;
  if (x < d(kLog_MinNorm) || !(x <= d(kLog_MaxNorm))) return detail::scalar_log(x);  // range mask: scalar log
  const std::uint64_t xb = detail::bits(x);
  // preserve mantissa, set input exponent to 2^(-10)
  double m = d((xb & kLog_ExpMask) | kLog_Two10);
  const double rcp = detail::rcp_rounded(m);
  // argument reduction: R = Mantissa*Rcp - 1
  const double mh = d(detail::bits(m) & kLog_HalfMask);
  m = m - mh;
  double hi = mh * rcp;
  m = m * rcp;
  hi = hi - d(kLog_One);
  const double r = m + hi;
  const bool big = d(kLog_Threshold) < rcp;
  double e = double(std::int32_t(xb >> 52));  // biased exponent (x > 0)
  double p = d(kLog_poly_coeff_1) * r;
  e = e - (big ? d(kLog_Bias) : d(kLog_Bias1));  // the table stores -log(0.5*mantissa) for larger mantissas
  p = p + d(kLog_poly_coeff_2);
  e = e * d(kLog_L2);
  double r2 = r * r;
  double q = d(kLog_poly_coeff_3) * r;
  p = p * r2;
  q = q + d(kLog_poly_coeff_4);
  q = q + p;
  r2 = r2 * q;
  const double rp = r + r2;
  // reconstruction: (exponent*log(2)) + (LogRcp + (R+poly)); LogRcp index = Rcp - 512
  const double l = d(CPPWRF_MVT(kLogLA, kLogLADev)[int(rcp) - 512]) + rp;
  return e + l;
}

// _ZGVbN2v_log10_sse4 (svml_d_log102_core_sse4.S)
CPPWRF_MV double log10(double x) {
  using namespace mvec_tables;
  using detail::d;
  if (x < d(kLog10_MinNorm) || !(x <= d(kLog10_MaxNorm))) return detail::scalar_log10(x);
  const std::uint64_t xb = detail::bits(x);
  double m = d((xb & kLog10_ExpMask) | kLog10_Two10);
  const double rcp = detail::rcp_rounded(m);
  const double mh = d(detail::bits(m) & kLog10_HalfMask);
  m = m - mh;
  double hi = mh * rcp;
  m = m * rcp;
  hi = hi - d(kLog10_One);
  const double r = m + hi;
  double p14 = d(kLog10_poly_coeff_0) * r;
  double r2 = r * r;
  p14 = p14 + d(kLog10_poly_coeff_1);
  double p15 = d(kLog10_poly_coeff_2) * r;
  p14 = p14 * r2;
  p15 = p15 + d(kLog10_poly_coeff_3);
  p15 = p15 + p14;
  const bool big = d(kLog10_Threshold) < rcp;
  r2 = r2 * p15;
  double e = double(std::int32_t(xb >> 52));
  double p11 = d(kLog10_poly_coeff_4) * r;
  e = e - (big ? d(kLog10_Bias) : d(kLog10_Bias1));
  p11 = p11 + r2;
  e = e * d(kLog10_L2);
  const double l = d(CPPWRF_MVT(kLog10LA, kLog10LADev)[int(rcp) - 512]) + p11;
  return e + l;
}
}  // namespace cppwrf::mvec
