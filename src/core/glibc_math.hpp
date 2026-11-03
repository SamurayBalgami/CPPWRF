#pragma once
// glibc's libm and libmvec (Debian 13, glibc 2.41-12, x86-64) on host and device, bit for bit.
//
// WRF's references and CPPWRF's CPU build call glibc's libm (the variants its ifuncs select on CPUs with FMA and
// AVX2: __ieee754_exp_fma, __ieee754_log_fma, ...) and, in the loops gfortran vectorises, libmvec's 2-lane SSE4
// functions (_ZGVbN2v_exp_sse4, ...). CUDA's math library rounds differently, so the CUDA FP64 build computes these
// functions here instead (core/backend.hpp: CPPWRF_GLIBC_MATH), on the device and on the host, the same on any machine.
// Each function is transcribed from the machine code of the libraries (objdump -d with the symbols of their debug
// files; Debian's 2.41 differs from the 2.41 tarball in places, e.g. the layout of __exp_data), statement by
// statement: every fused multiply-add of the FMA variants is an fma() here (gcc contracted a*b+c in them), every
// other operation a separate IEEE double operation in the machine code's order (the CUDA FP64 build has
// --fmad=false, the host -ffp-contract=off). The tables are the libraries' objects bit for bit
// (core/glibc_math_tables.hpp, tools/glibc_math_tables.py), read at the machine code's byte offsets. errno and the
// floating-point exception flags are not set. NaN results may differ from glibc's in sign and payload (CUDA's NaN).
// Proof: tools/glibc_math_proof.cpp (CPU, against the libraries) and tools/glibc_math_proof_cuda.cu (device
// against host).
#include <cmath>
#include <cstdint>
#include <cstring>

#include "core/glibc_math_tables.hpp"

#if defined(__CUDACC__)
// libm's deliberate divisions by zero and overflowing / underflowing constants (math_err.c, e_gamma_r.c)
#pragma nv_diagnostic push
#pragma nv_diag_suppress 39
#pragma nv_diag_suppress 222
#define CPPWRF_GM __host__ __device__ __forceinline__
#else
#define CPPWRF_GM inline
#endif

namespace cppwrf::glibc {
namespace detail {
#if defined(__CUDA_ARCH__)
#define CPPWRF_GT(name) ::cppwrf::glibc_tables::dev::name.v
#else
#define CPPWRF_GT(name) ::cppwrf::glibc_tables::name
#endif
CPPWRF_GM double asd(std::uint64_t b) {
#if defined(__CUDA_ARCH__)
  return __longlong_as_double(static_cast<long long>(b));
#else
  double v;
  std::memcpy(&v, &b, sizeof v);
  return v;
#endif
}
CPPWRF_GM std::uint64_t asu(double v) {
#if defined(__CUDA_ARCH__)
  return static_cast<std::uint64_t>(__double_as_longlong(v));
#else
  std::uint64_t b;
  std::memcpy(&b, &v, sizeof b);
  return b;
#endif
}
// a + k for a double's bits a and a constant k. CUDA 13.0's ptxas for sm_120 (-O1 and above) compiles an integer add of
// an immediate to bits moved from a double register (mov.b64; add.s64) as a DADD of the immediate read as a double;
// two 32-bit adds with carry on the device keep the integer add
CPPWRF_GM std::uint64_t addk(std::uint64_t a, std::uint64_t k) {
#if defined(__CUDA_ARCH__)
  unsigned lo, hi;
  asm("{\n\t.reg .u32 al, ah;\n\tmov.b64 {al, ah}, %2;\n\tadd.cc.u32 %0, al, %3;\n\taddc.u32 %1, ah, %4;\n\t}"
      : "=r"(lo), "=r"(hi)
      : "l"(a), "r"(unsigned(k)), "r"(unsigned(k >> 32)));
  return (std::uint64_t(hi) << 32) | lo;
#else
  return a + k;
#endif
}
CPPWRF_GM double fma(double a, double b, double c) {
#if defined(__CUDA_ARCH__)
  return __fma_rn(a, b, c);
#else
  return std::fma(a, b, c);
#endif
}
// roundpd $0 / roundsd $0: round to nearest, ties to even
CPPWRF_GM double rint(double a) {
#if defined(__CUDA_ARCH__)
  return ::rint(a);
#else
  return std::nearbyint(a);
#endif
}
// libm's error helpers (sysdeps/ieee754/dbl-64/math_err.c)
CPPWRF_GM double math_oflow(std::uint32_t sign) {
  return (sign ? -0x1p769 : 0x1p769) * 0x1p769;
}
CPPWRF_GM double math_uflow(std::uint32_t sign) {
  return (sign ? -0x1p-767 : 0x1p-767) * 0x1p-767;
}
CPPWRF_GM double math_invalid(double x) {
  const double y = x - x;
  return y / y;
}
CPPWRF_GM double math_divzero(std::uint32_t sign) {
  const double zero = 0.0;
  return (sign ? -1.0 : 1.0) / zero;
}
// word w (byte offset 8 w) of a table as a double / as bits
#define CPPWRF_GD(name, off) ::cppwrf::glibc::detail::asd(CPPWRF_GT(name)[(off) / 8])
#define CPPWRF_GU(name, off) (CPPWRF_GT(name)[(off) / 8])
}  // namespace detail

// ---- exp: __ieee754_exp_fma (libm 0x7cfc0) ----
// __exp_data: invln2N +0x0, negln2hiN +0x8, negln2loN +0x10, C2 +0x18, C3 +0x20, C4 +0x28, C5 +0x30, shift +0x38,
// tab +0xb0 (tail, sbits pairs)
CPPWRF_GM double exp(double x) {
  using namespace detail;
  const std::uint64_t ix = asu(x);
  std::uint32_t abstop = std::uint32_t(ix >> 52) & 0x7ff;
  if (abstop - 0x3c9u > 0x3eu) {
    if (std::int32_t(abstop - 0x3c9u) < 0) return x + 1.0;  // tiny x (and 0): 1 + x
    if (abstop > 0x408u) {
      if (ix == 0xfff0000000000000ull) return 0.0;
      if (abstop == 0x7ffu) return x + 1.0;
      return (ix >> 63) ? math_uflow(0) : math_oflow(0);
    }
    abstop = 0;  // large x: special case below
  }
  const double shift = CPPWRF_GD(kExpData, 0x38);
  double kd = fma(x, CPPWRF_GD(kExpData, 0x0), shift);
  const std::uint64_t ki = asu(kd);
  kd = kd - shift;
  double r = fma(kd, CPPWRF_GD(kExpData, 0x8), x);
  r = fma(kd, CPPWRF_GD(kExpData, 0x10), r);
  const std::uint64_t idx = 2 * (ki & 0x7f);
  const std::uint64_t top = ki << 45;
  const double p23 = fma(r, CPPWRF_GD(kExpData, 0x20), CPPWRF_GD(kExpData, 0x18));  // C2 + r C3
  const double tr = r + asd(CPPWRF_GT(kExpData)[0xb0 / 8 + idx]);                  // tail + r
  std::uint64_t sbits = CPPWRF_GT(kExpData)[0xb8 / 8 + idx] + top;
  const double r2 = r * r;
  const double p45 = fma(r, CPPWRF_GD(kExpData, 0x30), CPPWRF_GD(kExpData, 0x28));  // C4 + r C5
  const double t = fma(p23, r2, tr);
  const double tmp = fma(r2 * r2, p45, t);
  if (abstop == 0) {  // specialcase (0x7d0c0)
    if ((ki & 0x80000000u) == 0) {
      sbits -= 1009ull << 52;
      const double scale = asd(sbits);
      return fma(scale, tmp, scale) * 0x1p1009;
    }
    sbits += 1022ull << 52;
    const double scale = asd(sbits);
    const double st = tmp * scale;
    double y = scale + st;
    if (1.0 > y) {
      const double hi = y + 1.0;
      double lo = scale - y;
      lo = lo + st;
      double v = 1.0 - hi;
      v = v + y;
      v = v + lo;
      v = v + hi;
      y = v - 1.0;
      if (y == 0.0) y = 0.0;
    }
    return y * 0x1p-1022;
  }
  const double scale = asd(sbits);
  return fma(scale, tmp, scale);
}

// ---- log: __ieee754_log_fma (libm 0x7d660) ----
// __log_data: ln2hi +0x0, ln2lo +0x8, A[0..4] +0x10..+0x30, B[0..10] +0x38..+0x88, tab +0x90 ({invc, logc})
CPPWRF_GM double log(double x) {
  using namespace detail;
  std::uint64_t ix = asu(x);
  if (addk(ix, 0xc012000000000000ull) <= 0x308ffffffffffull) {  // |x - 1| small: ix - LO < HI - LO
    if (ix == 0x3ff0000000000000ull) return 0.0;
    const double r = x - 1.0;
    double b12 = fma(r, CPPWRF_GD(kLogData, 0x48), CPPWRF_GD(kLogData, 0x40));   // B1 + r B2
    double b45 = fma(r, CPPWRF_GD(kLogData, 0x60), CPPWRF_GD(kLogData, 0x58));   // B4 + r B5
    const double r2 = r * r;
    double b78 = fma(r, CPPWRF_GD(kLogData, 0x78), CPPWRF_GD(kLogData, 0x70));   // B7 + r B8
    b12 = fma(r2, CPPWRF_GD(kLogData, 0x50), b12);                                // + r2 B3
    b45 = fma(r2, CPPWRF_GD(kLogData, 0x68), b45);                                // + r2 B6
    const double r3 = r * r2;
    double p = fma(r2, CPPWRF_GD(kLogData, 0x80), b78);                           // + r2 B9
    p = fma(r3, CPPWRF_GD(kLogData, 0x88), p);                                    // + r3 B10
    p = fma(p, r3, b45);
    p = fma(p, r3, b12);
    const double rw = fma(r, 0x1p27, r);   // r + w, w = r 2^27
    const double rhi = fma(-0x1p27, r, rw);  // r + w - w
    const double b0 = CPPWRF_GD(kLogData, 0x38);
    const double rhi2 = rhi * rhi;
    const double rlo = r - rhi;
    const double hi = fma(rhi2, b0, r);
    const double rmhi = r - hi;
    const double rrhi = r + rhi;
    double lo = fma(rhi2, b0, rmhi);
    lo = fma(b0 * rlo, rrhi, lo);
    const double y = fma(p, r3, lo);
    return hi + y;
  }
  std::uint32_t top = std::uint32_t(ix >> 48);
  if (top - 0x10u > 0x7fdfu) {
    if (2 * ix == 0) return math_divzero(1);
    if (ix == 0x7ff0000000000000ull) return x;
    if ((top & 0x8000u) || (~top & 0x7ff0u) == 0) return math_invalid(x);
    ix = addk(asu(x * 0x1p52), 0xfcc0000000000000ull);  // subnormal: normalise, - 52 << 52
  }
  const std::uint64_t tmp = addk(ix, 0xc01a000000000000ull);  // ix - OFF
  const int i = int((tmp >> 45) & 0x7f);
  const std::int64_t k = std::int64_t(tmp) >> 52;
  const std::uint64_t iz = ix - (tmp & 0xfff0000000000000ull);
  const double kd = double(std::int32_t(k));
  const double z = asd(iz);
  const double invc = CPPWRF_GD(kLogData, 0x90 + 16 * i), logc = CPPWRF_GD(kLogData, 0x98 + 16 * i);
  const double w = fma(kd, CPPWRF_GD(kLogData, 0x0), logc);
  const double r = fma(z, invc, -1.0);
  const double a12 = fma(r, CPPWRF_GD(kLogData, 0x20), CPPWRF_GD(kLogData, 0x18));  // A1 + r A2
  const double hi = r + w;
  const double r2 = r * r;
  double lo = w - hi;
  lo = lo + r;
  lo = fma(kd, CPPWRF_GD(kLogData, 0x8), lo);
  const double rr2 = r * r2;
  double a34 = fma(r, CPPWRF_GD(kLogData, 0x30), CPPWRF_GD(kLogData, 0x28));        // A3 + r A4
  lo = fma(r2, CPPWRF_GD(kLogData, 0x10), lo);                                      // + r2 A0
  a34 = fma(a34, r2, a12);
  const double y = fma(rr2, a34, lo);
  return y + hi;
}

// ---- log10: __ieee754_log10 (libm 0x56440; sysdeps/ieee754/dbl-64/e_log10.c, no FMA variant; calls log) ----
CPPWRF_GM double log10(double x) {
  using namespace detail;
  std::int64_t hx = std::int64_t(asu(x));
  std::int32_t k = 0;
  if (hx < 0x0010000000000000ll) {
    if ((hx & 0x7fffffffffffffffll) == 0) return -0x1p54 / asd(asu(x) & 0x7fffffffffffffffull);  // log(+-0)
    if (hx < 0) return (x - x) / (x - x);                                     // log(-x) = NaN
    k = -54;
    x = x * 0x1p54;  // subnormal: scale up
    hx = std::int64_t(asu(x));
  }
  if (std::uint64_t(hx) > 0x7fefffffffffffffull) return x + x;
  const std::int64_t kk = std::int64_t(k + std::int32_t(hx >> 52) - 1023);
  const std::int64_t i = std::int64_t(std::uint64_t(kk) >> 63);
  const double y = double(kk + i);
  const double m = asd((std::uint64_t(hx) & 0x000fffffffffffffull) | (std::uint64_t(0x3ff - i) << 52));
  const double ylo = y * 0x1.9fef311f12b36p-42;   // log10_2lo
  const double z = glibc::log(m) * 0x1.bcb7b1526e50ep-2 + ylo;  // ivln10
  return z + y * 0x1.34413509f6p-2;               // log10_2hi
}

// ---- pow: __ieee754_pow_fma (libm 0x7dcf0) ----
namespace detail {
// checkint: 0 not an integer, 1 odd, 2 even (iy: a non-zero finite value)
CPPWRF_GM int checkint(std::uint64_t iy) {
  const int e = int(iy >> 52 & 0x7ff);
  if (e < 0x3ff) return 0;
  if (e > 0x3ff + 52) return 2;
  if (iy & ((1ull << (0x3ff + 52 - e)) - 1)) return 0;
  if (iy & (1ull << (0x3ff + 52 - e))) return 1;
  return 2;
}
CPPWRF_GM bool zeroinfnan(std::uint64_t i) { return 2 * i - 1 >= 2 * 0x7ff0000000000000ull - 1; }
CPPWRF_GM bool issignaling(double x) {
  return 2 * (asu(x) ^ 0x0008000000000000ull) > 2 * 0x7ff8000000000000ull;
}
}  // namespace detail

// __pow_log_data: ln2hi +0x0, ln2lo +0x8, A[0..6] +0x10..+0x40, tab +0x48 ({invc, pad, logc, logctail});
// __exp_data as exp()
CPPWRF_GM double pow(double x, double y) {
  using namespace detail;
  std::uint32_t sign_bias = 0;
  std::uint64_t ix = asu(x);
  const std::uint64_t iy = asu(y);
  std::uint32_t topx = std::uint32_t(ix >> 52);
  const std::uint32_t topy = std::uint32_t(iy >> 52);
  if (topx - 0x001u >= 0x7ffu - 0x001u || (topy & 0x7ff) - 0x3beu >= 0x43eu - 0x3beu) {
    if (zeroinfnan(iy)) {
      if (2 * iy == 0) return issignaling(x) ? x + y : 1.0;
      if (ix == 0x3ff0000000000000ull) return issignaling(y) ? x + y : 1.0;
      if (2 * ix > 2 * 0x7ff0000000000000ull || 2 * iy > 2 * 0x7ff0000000000000ull) return x + y;
      if (2 * ix == 2 * 0x3ff0000000000000ull) return 1.0;
      if ((2 * ix < 2 * 0x3ff0000000000000ull) == !(iy >> 63)) return 0.0;
      return y * y;
    }
    if (zeroinfnan(ix)) {
      double x2 = x * x;
      std::uint32_t sb = 0;
      if ((ix >> 63) && checkint(iy) == 1) {
        x2 = -x2;
        sb = 1;
      }
      if (2 * ix == 0 && (iy >> 63)) return math_divzero(sb);
      return (iy >> 63) ? 1 / x2 : x2;
    }
    if (ix >> 63) {
      const int yint = checkint(iy);
      if (yint == 0) return math_invalid(x);
      if (yint == 1) sign_bias = 0x800u << 7;
      ix &= 0x7fffffffffffffffull;
      topx &= 0x7ff;
    }
    if ((topy & 0x7ff) - 0x3beu >= 0x43eu - 0x3beu) {
      if (ix == 0x3ff0000000000000ull) return 1.0;
      if ((topy & 0x7ff) < 0x3be) return ix > 0x3ff0000000000000ull ? 1.0 + y : 1.0 - y;
      return (ix > 0x3ff0000000000000ull) == (topy < 0x800) ? math_oflow(0) : math_uflow(0);
    }
    if (topx == 0) {  // subnormal x
      ix = asu(x * 0x1p52);
      ix &= 0x7fffffffffffffffull;
      ix = addk(ix, 0ull - (52ull << 52));
    }
  }
  // log_inline
  const std::uint64_t tmp = addk(ix, 0xc0196aab00000000ull);  // ix - OFF
  const int i = int((tmp >> 45) & 0x7f);
  const double kd = double(std::int32_t(std::int64_t(tmp) >> 52));
  const int e = 0x48 + 32 * i;
  const double t1 = fma(kd, CPPWRF_GD(kPowLogData, 0x0), CPPWRF_GD(kPowLogData, e + 0x10));  // kd Ln2hi + logc
  const double lo1 = fma(kd, CPPWRF_GD(kPowLogData, 0x8), CPPWRF_GD(kPowLogData, e + 0x18)); // kd Ln2lo + logctail
  const double z = asd(ix - (tmp & 0xfff0000000000000ull));
  const double r = fma(z, CPPWRF_GD(kPowLogData, e), -1.0);
  const double ar = r * CPPWRF_GD(kPowLogData, 0x10);                                      // A0 r
  const double a12 = fma(r, CPPWRF_GD(kPowLogData, 0x20), CPPWRF_GD(kPowLogData, 0x18));   // A1 + r A2
  const double a34 = fma(r, CPPWRF_GD(kPowLogData, 0x30), CPPWRF_GD(kPowLogData, 0x28));   // A3 + r A4
  const double t2 = r + t1;
  double lo2 = t1 - t2;
  lo2 = lo2 + r;
  const double ar2 = r * ar;
  const double ar3 = r * ar2;
  const double lo3 = fma(ar, r, -ar2);
  const double hi = t2 + ar2;
  double q = fma(r, CPPWRF_GD(kPowLogData, 0x40), CPPWRF_GD(kPowLogData, 0x38));           // A5 + r A6
  double lo4 = t2 - hi;
  q = fma(q, ar2, a34);
  lo4 = lo4 + ar2;
  q = fma(ar2, q, a12);
  double lo = lo1 + lo2;
  lo = lo + lo3;
  lo = lo + lo4;
  lo = fma(ar3, q, lo);
  const double ly = hi + lo;
  double ltail = hi - ly;
  ltail = ltail + lo;
  const double ehi = y * ly;
  const double elo = fma(y, ltail, fma(ly, y, -ehi));
  // exp_inline (ehi, elo, sign_bias)
  const std::uint64_t ie = asu(ehi);
  std::uint32_t abstop = std::uint32_t(ie >> 52) & 0x7ff;
  if (abstop - 0x3c9u > 0x3eu) {
    if (std::int32_t(abstop - 0x3c9u) < 0) {
      const double one = ehi + 1.0;
      return sign_bias ? -one : one;
    }
    if (abstop > 0x408u) return (ie >> 63) ? math_uflow(sign_bias) : math_oflow(sign_bias);
    abstop = 0;
  }
  const double shift = CPPWRF_GD(kExpData, 0x38);
  double ek = fma(ehi, CPPWRF_GD(kExpData, 0x0), shift);
  const std::uint64_t ki = asu(ek);
  ek = ek - shift;
  double er = fma(ek, CPPWRF_GD(kExpData, 0x8), ehi);
  er = fma(ek, CPPWRF_GD(kExpData, 0x10), er);
  const std::uint64_t idx = 2 * (ki & 0x7f);
  std::uint64_t sbits = ((ki + sign_bias) << 45) + CPPWRF_GT(kExpData)[0xb8 / 8 + idx];
  er = elo + er;
  const double p23 = fma(er, CPPWRF_GD(kExpData, 0x20), CPPWRF_GD(kExpData, 0x18));
  const double tr = er + asd(CPPWRF_GT(kExpData)[0xb0 / 8 + idx]);
  const double r2 = er * er;
  const double p45 = fma(er, CPPWRF_GD(kExpData, 0x30), CPPWRF_GD(kExpData, 0x28));
  const double t = fma(p23, r2, tr);
  const double etmp = fma(r2 * r2, p45, t);
  if (abstop == 0) {  // specialcase (0x7e070)
    if ((ki & 0x80000000u) == 0) {
      sbits -= 1009ull << 52;
      const double scale = asd(sbits);
      return fma(scale, etmp, scale) * 0x1p1009;
    }
    sbits += 1022ull << 52;
    const double scale = asd(sbits);
    const double st = etmp * scale;
    double v = scale + st;
    if (1.0 > (v < 0 ? -v : v)) {
      const double one = v < 0.0 ? -1.0 : 1.0;
      double l = scale - v;
      l = l + st;
      const double h = v + one;
      double u = one - h;
      u = u + v;
      u = u + l;
      u = u + h;
      v = u - one;
      if (v == 0.0) v = asd(sbits & 0x8000000000000000ull);
    }
    return v * 0x1p-1022;
  }
  const double scale = asd(sbits);
  return fma(scale, etmp, scale);
}

// ---- atan: __atan_fma (libm 0x7ee30; sysdeps/ieee754/dbl-64/s_atan.c, constants of atnat.h, cij of uatan.tbl) ----
namespace detail {
CPPWRF_GM double copysign(double y, double x) {  // andnpd / andpd / orpd with the sign mask
  return asd((asu(y) & 0x7fffffffffffffffull) | (asu(x) & 0x8000000000000000ull));
}
CPPWRF_GM double fabs_(double x) { return asd(asu(x) & 0x7fffffffffffffffull); }  // andpd with the abs mask
}  // namespace detail
CPPWRF_GM double atan(double x) {
  using namespace detail;
  constexpr double TWO52 = 0x1p52, HPI = 0x1.921fb54442d18p0, HPI1 = 0x1.1a62633145c07p-54;
  constexpr double d3 = -0x1.5555555555555p-2, d5 = 0x1.99999999997fdp-3, d7 = -0x1.24924923f7603p-3,
                   d9 = 0x1.c71c6e5129a3bp-4, d11 = -0x1.7458022b13c25p-4, d13 = 0x1.375f08b31cbcep-4;
  const std::uint64_t ix = asu(x);
  const std::uint32_t ux = std::uint32_t(ix >> 32);
  if ((~ux & 0x7ff00000u) == 0 && ((ux & 0xfffffu) | std::uint32_t(ix)) != 0) return x + x;  // NaN
  const double u = x < 0 ? -x : x;
  if (u < 1.0) {
    if (!(u >= 0x1p-4)) {
      if (!(u >= 0x1.bb67ap-27)) return x;  // u < A
      const double v = x * x;               // A <= u < B
      double yy = fma(v, d13, d11);
      yy = fma(v, yy, d9);
      yy = fma(v, yy, d7);
      yy = fma(v, yy, d5);
      yy = fma(v, yy, d3);
      return fma(x * v, yy, x);
    }
    // B <= u < C
    const int i = int(fma(u, 256.0, TWO52) - TWO52) - 16;
    const int c = 7 * i;
    const double z = u - CPPWRF_GD(kAtanCij, 8 * c);
    double yy = fma(z, CPPWRF_GD(kAtanCij, 8 * (c + 6)), CPPWRF_GD(kAtanCij, 8 * (c + 5)));
    yy = fma(z, yy, CPPWRF_GD(kAtanCij, 8 * (c + 4)));
    yy = fma(z, yy, CPPWRF_GD(kAtanCij, 8 * (c + 3)));
    yy = fma(z, yy, CPPWRF_GD(kAtanCij, 8 * (c + 2)));
    return copysign(fma(yy, z, CPPWRF_GD(kAtanCij, 8 * (c + 1))), x);
  }
  if (u < 16.0) {  // C <= u < D
    const double w = 1.0 / u;
    const double t1 = w * u;
    const double t2 = fma(u, w, -t1);
    double ww = 1.0 - t1;
    ww = ww - t2;
    const int i = int(fma(w, 256.0, TWO52) - TWO52) - 16;
    const int c = 7 * i;
    const double z = fma(ww, w, w - CPPWRF_GD(kAtanCij, 8 * c));
    double yy = fma(z, CPPWRF_GD(kAtanCij, 8 * (c + 6)), CPPWRF_GD(kAtanCij, 8 * (c + 5)));
    yy = fma(z, yy, CPPWRF_GD(kAtanCij, 8 * (c + 4)));
    yy = fma(z, yy, CPPWRF_GD(kAtanCij, 8 * (c + 3)));
    yy = fma(z, yy, CPPWRF_GD(kAtanCij, 8 * (c + 2)));
    yy = fma(-yy, z, HPI1);
    const double tt = HPI - CPPWRF_GD(kAtanCij, 8 * (c + 1));
    return copysign(tt + yy, x);
  }
  if (u < 0x1.49ff2p52) {  // D <= u < E
    const double w = 1.0 / u;
    const double t1 = w * u;
    const double t3 = HPI - w;
    const double v = w * w;
    double yy = fma(v, d13, d11);
    yy = fma(v, yy, d9);
    double cor = HPI - t3;
    yy = fma(v, yy, d7);
    const double t2 = fma(u, w, -t1);
    double ww = 1.0 - t1;
    yy = fma(v, yy, d5);
    cor = cor - w;
    cor = cor + HPI1;
    yy = fma(v, yy, d3);
    const double wv = w * v;
    ww = ww - t2;
    const double a = fma(-ww, w, cor);
    const double b = fma(-wv, yy, a);
    return copysign(t3 + b, x);
  }
  return x > 0 ? HPI : -HPI;
}

// ---- expm1: __expm1_fma (libm 0x7feb0; sysdeps/ieee754/dbl-64/s_expm1.c) ----
namespace detail {
CPPWRF_GM double with_high(double y, std::uint32_t hi) {  // SET_HIGH_WORD
  return asd((std::uint64_t(hi) << 32) | (asu(y) & 0xffffffffull));
}
}  // namespace detail
CPPWRF_GM double expm1(double x) {
  using namespace detail;
  constexpr double ln2_hi = 0x1.62e42feep-1, ln2_lo = 0x1.a39ef35793c76p-33, invln2 = 0x1.71547652b82fep0;
  constexpr double Q1 = -0x1.11111111110f4p-5, Q2 = 0x1.a01a019fe5585p-10, Q3 = -0x1.4ce199eaadbb7p-14,
                   Q4 = 0x1.0cfca86e65239p-18, Q5 = -0x1.afdb76e09c32dp-23;
  const std::uint32_t hxw = std::uint32_t(asu(x) >> 32), hx = hxw & 0x7fffffffu, xsb = hxw & 0x80000000u;
  if (hx > 0x40862e41u) {  // |x| >= 709.78, inf, NaN
    if (hx > 0x7fefffffu) {
      if (((hxw & 0xfffffu) | std::uint32_t(asu(x))) != 0) return x + x;  // NaN
      return xsb ? -1.0 : x;
    }
    if (x > 0x1.62e42fefa39efp9) return 0x1.7e43c8800759cp996 * 0x1.7e43c8800759cp996;  // overflow
  }
  if (hx > 0x40436879u && xsb) return 0x1.56e1fc2f8f359p-997 - 1.0;  // x < -56 ln2: tiny - 1
  int k = 0;
  double c = 0.0;
  if (hx > 0x3fd62e42u) {  // |x| > 0.5 ln2: argument reduction
    double hi, lo;
    if (hx <= 0x3ff0a2b1u) {  // |x| < 1.5 ln2
      if (xsb == 0) {
        hi = x - ln2_hi;
        lo = ln2_lo;
        k = 1;
      } else {
        hi = x + ln2_hi;
        lo = -ln2_lo;
        k = -1;
      }
    } else {
      const double kf = x * invln2;
      k = int(kf + (xsb == 0 ? 0.5 : -0.5));
      const double t = double(k);
      hi = fma(-t, ln2_hi, x);
      lo = t * ln2_lo;
    }
    x = hi - lo;
    c = hi - x;
    c = c - lo;
  } else if (hx <= 0x3c8fffffu) {  // |x| < 2^-54
    const double t = x + 0x1.7e43c8800759cp996;
    return x - (t - t);
  }
  // x in the primary range
  const double hfx = x * 0.5;
  const double hxs = x * hfx;
  const double R2 = fma(hxs, Q3, Q2), R3 = fma(hxs, Q5, Q4);
  const double h2 = hxs * hxs;
  const double R1 = fma(hxs, Q1, 1.0);
  const double h4 = h2 * h2;
  const double r1 = fma(h4, R3, fma(h2, R2, R1));
  const double t = fma(-r1, hfx, 3.0);
  double e = (r1 - t) / fma(-x, t, 6.0);
  e = e * hxs;
  if (k == 0) return x - fma(e, x, -hxs);
  e = fma(e - c, x, -c);
  e = e - hxs;
  if (k == -1) return fma(0.5, x - e, -0.5);
  if (k == 1) {
    if (-0.25 > x) return (e - (x + 0.5)) * -2.0;
    return fma(x - e, 2.0, 1.0);
  }
  const std::uint32_t kk = std::uint32_t(k) << 20;
  if (std::uint32_t(k + 1) > 57u) {  // k <= -2 or k > 56
    double y = 1.0 - (e - x);
    y = with_high(y, std::uint32_t(asu(y) >> 32) + kk);
    return y - 1.0;
  }
  if (k > 19) {
    const double t2 = asd(std::uint64_t(std::uint32_t(0x3ff - k) << 20) << 32);  // 2^-k
    double y = x - (e + t2);
    y = y + 1.0;
    return with_high(y, std::uint32_t(asu(y) >> 32) + kk);
  }
  const double t1 = asd(std::uint64_t(0x3ff00000u - (0x200000u >> k)) << 32);  // 1 - 2^-k
  const double y = t1 - (e - x);
  return with_high(y, std::uint32_t(asu(y) >> 32) + kk);
}

// ---- tanh: __tanh_sse2 (libm 0x5e560, the function tanh resolves to; sysdeps/ieee754/dbl-64/s_tanh.c; calls expm1) --
CPPWRF_GM double tanh(double x) {
  using namespace detail;
  const std::uint64_t b = asu(x);
  const std::int32_t jx = std::int32_t(b >> 32);
  const std::uint32_t ix = std::uint32_t(jx) & 0x7fffffffu;
  if (ix > 0x7fefffffu) {  // inf or NaN
    const double r = 1.0 / x;
    return jx >= 0 ? r + 1.0 : r - 1.0;
  }
  double z;
  if (ix > 0x4035ffffu) {
    z = 1.0 - 0x1.56e1fc2f8f359p-997;  // |x| >= 22: 1 - tiny
  } else {
    if ((ix | std::uint32_t(b)) == 0) return x;  // +-0
    const double ax = asd(b & 0x7fffffffffffffffull);
    if (ix <= 0x3c7fffffu) return (x + 1.0) * x;  // |x| < 2^-55: x * (1 + x)
    if (ix > 0x3fefffffu) {                       // |x| >= 1
      const double t = expm1(ax + ax);
      z = 1.0 - 2.0 / (t + 2.0);
    } else {
      const double t = expm1(ax * -2.0);
      z = -t / (t + 2.0);
    }
  }
  return jx >= 0 ? z : -z;
}

// ---- sin, cos: __sin_fma (libm 0x81740), __cos_fma (libm 0x81f50); sysdeps/ieee754/dbl-64/s_sin.c (constants of
// usncs.h, __sincostab of sincostab.c), __branred of branred.c (libm 0x70a40, compiled without FMA) ----
// The rounding-mode switch of SET_RESTORE_ROUND_53BIT is left out (round to nearest throughout).
namespace detail {
// __branred: x (|x| >= 105414350) - n pi/2 as a + aa, returns n mod 4 (branred.c; plain double operations)
CPPWRF_GM int branred(double x, double* a, double* aa) {
  constexpr double tm600 = 0x1p-600, split = 134217729.0, tm24 = 0x1p-24, big = 0x1.8p52, big1 = 0x1.8p54,
                   hp0 = 0x1.921fb54442d18p0, hp1 = 0x1.1a62633145c07p-54, mp1 = 0x1.921fb58p0,
                   mp2 = -0x1.dde974p-27;
  double r[6], s, t, sum, b, bb, sum1, sum2, b1, bb1, b2, bb2, x1, x2, t1, t2;
  x *= tm600;
  t = x * split;  // split x to two numbers
  x1 = t - (t - x);
  x2 = x - x1;
  sum = 0;
  for (int part = 0; part < 2; ++part) {
    const double xp = part == 0 ? x1 : x2;
    int k = int((asu(xp) >> 52) & 2047);
    k = (k - 450) / 24;
    if (k < 0) k = 0;
    double gor = asd((std::uint64_t(0x63f00000u - std::uint32_t((k * 24) << 20)) << 32));  // t576 with exponent - 24 k
    for (int i = 0; i < 6; i++) {
      r[i] = xp * CPPWRF_GD(kBranredToverp, 8 * (k + i)) * gor;
      gor *= tm24;
    }
    for (int i = 0; i < 3; i++) {
      s = (r[i] + big) - big;
      sum += s;
      r[i] -= s;
    }
    t = 0;
    for (int i = 0; i < 6; i++) t += r[5 - i];
    bb = (((((r[0] - t) + r[1]) + r[2]) + r[3]) + r[4]) + r[5];
    s = (t + big) - big;
    sum += s;
    t -= s;
    b = t + bb;
    bb = (t - b) + bb;
    s = (sum + big1) - big1;
    sum -= s;
    if (part == 0) {
      b1 = b;
      bb1 = bb;
      sum1 = sum;
      sum = 0;
    } else {
      b2 = b;
      bb2 = bb;
      sum2 = sum;
    }
  }
  sum = sum1 + sum2;
  b = b1 + b2;
  bb = (fabs_(b1) > fabs_(b2)) ? (b1 - b) + b2 : (b2 - b) + b1;
  if (b > 0.5) {
    b -= 1.0;
    sum += 1.0;
  } else if (b < -0.5) {
    b += 1.0;
    sum -= 1.0;
  }
  s = b + (bb + bb1 + bb2);
  t = ((b - s) + bb) + (bb1 + bb2);
  b = s * split;
  t1 = b - (b - s);
  t2 = s - t1;
  b = s * hp0;
  bb = (((t1 * mp1 - b) + t1 * mp2) + t2 * mp1) + (t2 * mp2 + s * hp1 + t * hp0);
  s = b + bb;
  t = (b - s) + bb;
  *a = s;
  *aa = t;
  return int(sum) & 3;
}
constexpr double kSinBig = 0x1.8p45, kSn3 = -0x1.5555555555515p-3, kSn5 = 0x1.11110e829872fp-7,
                 kCs2 = 0.5, kCs4 = -0x1.5555555555535p-5, kCs6 = 0x1.6c16bedd9e239p-10;
// TAYLOR_SIN and do_sin of s_sin.c
CPPWRF_GM double do_sin(double x, double dx) {
  const double xold = x;
  if (fabs_(x) < 0.126) {
    constexpr double s1 = -0x1.5555555555555p-3, s2 = 0x1.1111111110ecep-7, s3 = -0x1.a01a019db08b8p-13,
                     s4 = 0x1.71de27b9a7ed9p-19, s5 = -0x1.addffc2fcdf59p-26;
    const double xx = x * x;
    double p = fma(xx, s5, s4);
    p = fma(xx, p, s3);
    p = fma(xx, p, s2);
    p = fma(xx, p, s1);
    const double t = fma(xx, fma(p, x, -(0.5 * dx)), dx);
    return x + t;
  }
  if (x <= 0) dx = -dx;
  const double u = kSinBig + fabs_(x);
  x = fabs_(x) - (u - kSinBig);
  const double xx = x * x;
  const double s = x + fma(x * xx, fma(xx, kSn5, kSn3), dx);
  const double c = fma(x, dx, xx * fma(xx, fma(xx, kCs6, kCs4), kCs2));
  const int k = int(std::uint32_t(asu(u)) << 2);
  const double sn = CPPWRF_GD(kSinCosTab, 8 * k), ssn = CPPWRF_GD(kSinCosTab, 8 * (k + 1)),
               cs = CPPWRF_GD(kSinCosTab, 8 * (k + 2)), ccs = CPPWRF_GD(kSinCosTab, 8 * (k + 3));
  const double cor = fma(s, cs, fma(-c, sn, fma(s, ccs, ssn)));
  return copysign(sn + cor, xold);
}
// do_cos of s_sin.c
CPPWRF_GM double do_cos(double x, double dx) {
  if (x < 0) dx = -dx;
  const double u = kSinBig + fabs_(x);
  x = fabs_(x) - (u - kSinBig) + dx;
  const double xx = x * x;
  const double s = fma(x * xx, fma(xx, kSn5, kSn3), x);
  const double c = xx * fma(xx, fma(xx, kCs6, kCs4), kCs2);
  const int k = int(std::uint32_t(asu(u)) << 2);
  const double sn = CPPWRF_GD(kSinCosTab, 8 * k), ssn = CPPWRF_GD(kSinCosTab, 8 * (k + 1)),
               cs = CPPWRF_GD(kSinCosTab, 8 * (k + 2)), ccs = CPPWRF_GD(kSinCosTab, 8 * (k + 3));
  const double cor = fma(-s, sn, fma(-c, cs, fma(-s, ssn, ccs)));
  return cs + cor;
}
// reduce_sincos of s_sin.c (2.426265 < |x| < 105414350)
CPPWRF_GM int reduce_sincos(double x, double* a, double* da) {
  constexpr double toint = 0x1.8p52, hpinv = 0x1.45f306dc9c883p-1, mp1 = 0x1.921fb58p0, mp2 = -0x1.dde973cp-27,
                   pp3 = -0x1.cb3b398p-55, pp4 = -0x1.d747f23e32ed7p-83;
  const double t = fma(x, hpinv, toint);
  const double xn = t - toint;
  double y = fma(-xn, mp1, x);
  y = fma(-xn, mp2, y);
  const int n = int(std::uint32_t(asu(t)) & 3);
  const double t2 = fma(-xn, pp3, y);
  double db = fma(-xn, pp3, y - t2);
  const double b = fma(-xn, pp4, t2);
  db = db + fma(-xn, pp4, t2 - b);
  *a = b;
  *da = db;
  return n;
}
CPPWRF_GM double do_sincos(double a, double da, int n) {
  const double r = (n & 1) ? do_cos(a, da) : do_sin(a, da);
  return (n & 2) ? -r : r;
}
}  // namespace detail
CPPWRF_GM double sin(double x) {
  using namespace detail;
  constexpr double hp0 = 0x1.921fb54442d18p0, hp1 = 0x1.1a62633145c07p-54;
  const std::uint64_t b = asu(x);
  const std::uint32_t k = std::uint32_t(b >> 32) & 0x7fffffffu;
  if (k < 0x3e500000u) return x;
  if (k < 0x3feb6000u) return do_sin(x, 0);
  if (k < 0x400368fdu) return copysign(do_cos(hp0 - fabs_(x), hp1), x);
  double a, da;
  if (k < 0x419921FBu) {
    const int n = reduce_sincos(x, &a, &da);
    return do_sincos(a, da, n);
  }
  if (k < 0x7ff00000u) {
    const int n = branred(x, &a, &da);
    return do_sincos(a, da, n);
  }
  return x / x;
}
CPPWRF_GM double cos(double x) {
  using namespace detail;
  constexpr double hp0 = 0x1.921fb54442d18p0, hp1 = 0x1.1a62633145c07p-54;
  const std::uint64_t b = asu(x);
  const std::uint32_t k = std::uint32_t(b >> 32) & 0x7fffffffu;
  if (k < 0x3e400000u) return 1.0;
  if (k < 0x3feb6000u) return do_cos(x, 0);
  double a, da;
  if (k < 0x400368fdu) {
    const double y = hp0 - fabs_(x);
    a = y + hp1;
    da = (y - a) + hp1;
    return do_sin(a, da);
  }
  if (k < 0x419921FBu) {
    const int n = reduce_sincos(x, &a, &da);
    return do_sincos(a, da, n + 1);
  }
  if (k < 0x7ff00000u) {
    const int n = branred(x, &a, &da);
    return do_sincos(a, da, n + 1);
  }
  return x / x;
}

// ---- sinh: __ieee754_sinh_fma (libm 0x7e740; sysdeps/ieee754/dbl-64/e_sinh.c; calls __expm1_fma, __ieee754_exp_fma) --
CPPWRF_GM double sinh(double x) {
  using namespace detail;
  const std::uint64_t b = asu(x);
  const std::int32_t jx = std::int32_t(b >> 32);
  const std::uint32_t ix = std::uint32_t(jx) & 0x7fffffffu;
  if (ix >= 0x7ff00000u) return x + x;
  const double h = jx < 0 ? -0.5 : 0.5;
  const double ax = fabs_(x);
  if (ix < 0x40360000u) {    // |x| < 22
    if (ix < 0x3e300000u) {  // |x| < 2^-28
      if (1.0e307 + x > 1.0) return x;
    }
    const double t = expm1(ax);
    if (ix < 0x3ff00000u) return h * fma(t, 2.0, -(t * t / (t + 1.0)));
    return h * (t + t / (t + 1.0));
  }
  if (ix < 0x40862e42u) return h * exp(ax);
  const std::uint32_t lx = std::uint32_t(b);
  if (ix < 0x408633ceu || (ix == 0x408633ceu && lx <= 0x8fb9f87du)) {
    const double w = exp(0.5 * ax);
    const double t = h * w;
    return t * w;
  }
  return x * 1.0e307;
}

// ---- erf: __erf (libm 0x58d20; sysdeps/ieee754/dbl-64/s_erf.c, compiled without FMA; calls __ieee754_exp_fma) ----
namespace detail {
struct ErfConst {
  static constexpr double
    tiny = 1e-300,
    half = 5.00000000000000000000e-01,
    one = 1.00000000000000000000e+00,
    two = 2.00000000000000000000e+00,
    erx = 8.45062911510467529297e-01,
    efx = 1.28379167095512586316e-01,
    pp0 = 1.28379167095512558561e-01,
    pp1 = -3.25042107247001499370e-01,
    pp2 = -2.84817495755985104766e-02,
    pp3 = -5.77027029648944159157e-03,
    pp4 = -2.37630166566501626084e-05,
    qq0 = 0.0,
    qq1 = 3.97917223959155352819e-01,
    qq2 = 6.50222499887672944485e-02,
    qq3 = 5.08130628187576562776e-03,
    qq4 = 1.32494738004321644526e-04,
    qq5 = -3.96022827877536812320e-06,
    pa0 = -2.36211856075265944077e-03,
    pa1 = 4.14856118683748331666e-01,
    pa2 = -3.72207876035701323847e-01,
    pa3 = 3.18346619901161753674e-01,
    pa4 = -1.10894694282396677476e-01,
    pa5 = 3.54783043256182359371e-02,
    pa6 = -2.16637559486879084300e-03,
    qa0 = 0.0,
    qa1 = 1.06420880400844228286e-01,
    qa2 = 5.40397917702171048937e-01,
    qa3 = 7.18286544141962662868e-02,
    qa4 = 1.26171219808761642112e-01,
    qa5 = 1.36370839120290507362e-02,
    qa6 = 1.19844998467991074170e-02,
    ra0 = -9.86494403484714822705e-03,
    ra1 = -6.93858572707181764372e-01,
    ra2 = -1.05586262253232909814e+01,
    ra3 = -6.23753324503260060396e+01,
    ra4 = -1.62396669462573470355e+02,
    ra5 = -1.84605092906711035994e+02,
    ra6 = -8.12874355063065934246e+01,
    ra7 = -9.81432934416914548592e+00,
    sa0 = 0.0,
    sa1 = 1.96512716674392571292e+01,
    sa2 = 1.37657754143519042600e+02,
    sa3 = 4.34565877475229228821e+02,
    sa4 = 6.45387271733267880336e+02,
    sa5 = 4.29008140027567833386e+02,
    sa6 = 1.08635005541779435134e+02,
    sa7 = 6.57024977031928170135e+00,
    sa8 = -6.04244152148580987438e-02,
    rb0 = -9.86494292470009928597e-03,
    rb1 = -7.99283237680523006574e-01,
    rb2 = -1.77579549177547519889e+01,
    rb3 = -1.60636384855821916062e+02,
    rb4 = -6.37566443368389627722e+02,
    rb5 = -1.02509513161107724954e+03,
    rb6 = -4.83519191608651397019e+02,
    sb0 = 0.0,
    sb1 = 3.03380607434824582924e+01,
    sb2 = 3.25792512996573918826e+02,
    sb3 = 1.53672958608443695994e+03,
    sb4 = 3.19985821950859553908e+03,
    sb5 = 2.55305040643316442583e+03,
    sb6 = 4.74528541206955367215e+02,
    sb7 = -2.24409524465858183362e+01;
};
}  // namespace detail
CPPWRF_GM double erf(double x) {
  using namespace detail;
  using C = ErfConst;
  constexpr double tiny = C::tiny, one = C::one, erx = C::erx, efx = C::efx;
  const std::int32_t hx = std::int32_t(asu(x) >> 32);
  const std::int32_t ix = hx & 0x7fffffff;
  double R, S, P, Q, s, y, z, r;
  if (ix >= 0x7ff00000) {
    const int i = int((std::uint32_t(hx) >> 31) << 1);
    return double(1 - i) + one / x;
  }
  if (ix < 0x3feb0000) {
    double r1, r2, s1, s2, s3, z2, z4;
    if (ix < 0x3e300000) {
      if (ix < 0x00800000) return 0.0625 * (16.0 * x + (16.0 * efx) * x);
      return x + efx * x;
    }
    z = x * x;
    r1 = C::pp0 + z * C::pp1;
    z2 = z * z;
    r2 = C::pp2 + z * C::pp3;
    z4 = z2 * z2;
    s1 = one + z * C::qq1;
    s2 = C::qq2 + z * C::qq3;
    s3 = C::qq4 + z * C::qq5;
    r = r1 + z2 * r2 + z4 * C::pp4;
    s = s1 + z2 * s2 + z4 * s3;
    y = r / s;
    return x + x * y;
  }
  if (ix < 0x3ff40000) {
    double s2, s4, s6, P1, P2, P3, P4, Q1, Q2, Q3, Q4;
    s = fabs_(x) - one;
    P1 = C::pa0 + s * C::pa1;
    s2 = s * s;
    Q1 = one + s * C::qa1;
    s4 = s2 * s2;
    P2 = C::pa2 + s * C::pa3;
    s6 = s4 * s2;
    Q2 = C::qa2 + s * C::qa3;
    P3 = C::pa4 + s * C::pa5;
    Q3 = C::qa4 + s * C::qa5;
    P4 = C::pa6;
    Q4 = C::qa6;
    P = P1 + s2 * P2 + s4 * P3 + s6 * P4;
    Q = Q1 + s2 * Q2 + s4 * Q3 + s6 * Q4;
    if (hx >= 0) return erx + P / Q;
    return -erx - P / Q;
  }
  if (ix >= 0x40180000) {
    if (hx >= 0) return one - tiny;
    return tiny - one;
  }
  x = fabs_(x);
  s = one / (x * x);
  if (ix < 0x4006DB6E) {
    double R1, R2, R3, R4, S1, S2, S3, S4, s2, s4, s6, s8;
    R1 = C::ra0 + s * C::ra1;
    s2 = s * s;
    S1 = one + s * C::sa1;
    s4 = s2 * s2;
    R2 = C::ra2 + s * C::ra3;
    s6 = s4 * s2;
    S2 = C::sa2 + s * C::sa3;
    s8 = s4 * s4;
    R3 = C::ra4 + s * C::ra5;
    S3 = C::sa4 + s * C::sa5;
    R4 = C::ra6 + s * C::ra7;
    S4 = C::sa6 + s * C::sa7;
    R = R1 + s2 * R2 + s4 * R3 + s6 * R4;
    S = S1 + s2 * S2 + s4 * S3 + s6 * S4 + s8 * C::sa8;
  } else {
    double R1, R2, R3, S1, S2, S3, S4, s2, s4, s6;
    R1 = C::rb0 + s * C::rb1;
    s2 = s * s;
    S1 = one + s * C::sb1;
    s4 = s2 * s2;
    R2 = C::rb2 + s * C::rb3;
    s6 = s4 * s2;
    S2 = C::sb2 + s * C::sb3;
    R3 = C::rb4 + s * C::rb5;
    S3 = C::sb4 + s * C::sb5;
    S4 = C::sb6 + s * C::sb7;
    R = R1 + s2 * R2 + s4 * R3 + s6 * C::rb6;
    S = S1 + s2 * S2 + s4 * S3 + s6 * S4;
  }
  z = asd(asu(x) & 0xffffffff00000000ull);  // SET_LOW_WORD (z, 0)
  r = exp(-z * z - 0.5625) * exp((z - x) * (z + x) + R / S);
  if (hx >= 0) return one - r / x;
  return r / x - one;
}

// ---- asin, acos: __ieee754_asin_fma (libm 0x7a650), __ieee754_acos_fma (libm 0x7ad30); sysdeps/ieee754/dbl-64/
// e_asin.c (constants of uasncs.h, tables asncs, inroot, powtwo) ----
namespace detail {
constexpr double kAsF1 = 1.66666666666664110590506577996662E-01, kAsF2 = 7.50000000026122686814431784722623E-02,
                 kAsF3 = 4.46428561421059750978517350006940E-02, kAsF4 = 3.03821268582119319911193410625235E-02,
                 kAsF5 = 2.23551211026525610742786300334557E-02, kAsF6 = 1.81382903404565056280372531963613E-02,
                 kAsRt0 = 9.99999999859990725855365213134618E-01, kAsRt1 = 4.99999999495955425917856814202739E-01,
                 kAsRt2 = 3.75017500867345182581453026130850E-01, kAsRt3 = 3.12523626554518656309172508769531E-01,
                 kAsHp0 = 0x1.921fb54442d18p0, kAsHp1 = 0x1.1a62633145c07p-54;
#define CPPWRF_ASN(i) CPPWRF_GD(kAsnCs, 8 * (i))
// the table branches: xx = |x| - asncs[n]; t = asncs[n+1] xx + xx^2 (asncs[n+2] + ... + xx^(d-1) asncs[n+d+1]) +
// asncs[n+d+2] (Horner, each step a fused multiply-add); the result's table value asncs[n+d+3] is added by the caller
CPPWRF_GM double asn_poly(double xx, int n, int d) {
  double q = CPPWRF_ASN(n + d + 1);
  for (int i = n + d; i >= n + 2; --i) q = fma(xx, q, CPPWRF_ASN(i));
  const double p = fma(xx * xx, q, CPPWRF_ASN(n + d + 2));
  return fma(xx, CPPWRF_ASN(n + 1), p);
}
// the branch of 0.125 <= |x| < 0.96875: index n and polynomial degree d
CPPWRF_GM void asn_branch(std::int32_t k, int& n, int& d) {
  if (k < 0x3fe00000) {
    n = k < 0x3fd00000 ? 11 * ((k & 0x000fffff) >> 15) : 11 * ((k & 0x000fffff) >> 14) + 352;
    d = 5;
  } else if (k < 0x3fe80000) {
    n = 1056 + ((k & 0x000fe000) >> 11) * 3;
    d = 6;
  } else if (k < 0x3fed8000) {
    n = 992 + ((k & 0x000fe000) >> 13) * 13;
    d = 7;
  } else if (k < 0x3fee8000) {
    n = 884 + ((k & 0x000fe000) >> 13) * 14;
    d = 8;
  } else {
    n = 768 + ((k & 0x000fe000) >> 13) * 15;
    d = 9;
  }
}
// 0.96875 <= |x| < 1: z = (1 - |x|) / 2, the root sqrt(z) as y + cc and the polynomial p
CPPWRF_GM double asn_root(double z, double& y, double& cc, bool acos_split) {
  const std::uint32_t k = std::uint32_t(asu(z) >> 32);
  double t = CPPWRF_GD(kInRoot, 8 * ((k & 0x001fffffu) >> 14)) * CPPWRF_GD(kPowTwo, 8 * (511 - int(k >> 21)));
  const double r = fma(-(t * t), z, 1.0);
  t = t * fma(r, fma(r, fma(r, kAsRt3, kAsRt2), kAsRt1), kAsRt0);
  const double c = t * z;
  const double a = fma(-c, t * 0.5, 1.5);
  if (acos_split) {
    y = fma(-0x1p27, c, fma(c, 0x1p27, c));  // (t27 c + c) - t27 c
  } else {
    y = (c + 0x1p24) - 0x1p24;
  }
  cc = fma(-y, y, z) / fma(a, c, y);
  double p = fma(z, kAsF6, kAsF5);
  p = fma(z, p, kAsF4);
  p = fma(z, p, kAsF3);
  p = fma(z, p, kAsF2);
  p = fma(z, p, kAsF1);
  return p * z;
}
}  // namespace detail
CPPWRF_GM double asin(double x) {
  using namespace detail;
  const std::int32_t m = std::int32_t(asu(x) >> 32);
  const std::int32_t k = m & 0x7fffffff;
  if (k < 0x3e500000) return x;
  if (k < 0x3fc00000) {
    const double x2 = x * x;
    double q = fma(x2, kAsF6, kAsF5);
    q = fma(x2, q, kAsF4);
    q = fma(x2, q, kAsF3);
    q = fma(x2, q, kAsF2);
    q = fma(x2, q, kAsF1);
    return fma(q, x2 * x, x);
  }
  if (k < 0x3fef0000) {
    int n, d;
    asn_branch(k, n, d);
    const double xx = (m > 0 ? x : -x) - CPPWRF_ASN(n);
    const double res = CPPWRF_ASN(n + d + 3) + asn_poly(xx, n, d);
    return m > 0 ? res : -res;
  }
  if (k < 0x3ff00000) {
    const double z = (m > 0 ? 1.0 - x : x + 1.0) * 0.5;
    double y, cc;
    const double p = asn_root(z, y, cc, false);
    const double cor = fma(-((y + cc) + (y + cc)), p, fma(-2.0, cc, kAsHp1));
    const double res = cor + fma(-2.0, y, kAsHp0);
    return m > 0 ? res : -res;
  }
  if (k == 0x3ff00000 && std::uint32_t(asu(x)) == 0) return m > 0 ? kAsHp0 : -kAsHp0;
  return (x - x) / (x - x);
}
CPPWRF_GM double acos(double x) {
  using namespace detail;
  const std::int32_t m = std::int32_t(asu(x) >> 32);
  const std::int32_t k = m & 0x7fffffff;
  if (k < 0x3c880000) return kAsHp0;
  if (k < 0x3fc00000) {
    const double x2 = x * x;
    double q = fma(x2, kAsF6, kAsF5);
    q = fma(x2, q, kAsF4);
    q = fma(x2, q, kAsF3);
    q = fma(x2, q, kAsF2);
    q = fma(x2, q, kAsF1);
    const double r = kAsHp0 - x;
    const double cor = fma(-q, x * x2, (((kAsHp0 - r) - x) + kAsHp1));
    return r + cor;
  }
  if (k < 0x3fef0000) {
    int n, d;
    asn_branch(k, n, d);
    const double xx = (m > 0 ? x : -x) - CPPWRF_ASN(n);
    const double t = asn_poly(xx, n, d);
    const double yv = m > 0 ? kAsHp0 - CPPWRF_ASN(n + d + 3) : kAsHp0 + CPPWRF_ASN(n + d + 3);
    const double tv = m > 0 ? kAsHp1 - t : kAsHp1 + t;
    return tv + yv;
  }
  if (k < 0x3ff00000) {
    const double z = (m > 0 ? 1.0 - x : x + 1.0) * 0.5;
    double y, cc;
    const double p = asn_root(z, y, cc, true);
    const double py = p * (y + cc);
    if (m < 0) {
      const double res = ((kAsHp1 - cc) - py) + (kAsHp0 - y);
      return res + res;
    }
    const double res = (cc + py) + y;
    return res + res;
  }
  if (k == 0x3ff00000 && std::uint32_t(asu(x)) == 0) return m > 0 ? 0.0 : 0x1.921fb54442d18p1;
  return (x - x) / (x - x);
}
#undef CPPWRF_ASN

// ---- tan: __tan_fma (libm 0x84db0; sysdeps/ieee754/dbl-64/s_tan.c, constants of utan.h, xfg of utan.tbl; EADD,
// DIV2 and MUL12 of dla.h with FMA) ----
namespace detail {
struct TanC {
  static constexpr double d3 = 0x1.5555555555555p-2, d5 = 0x1.11111111107c6p-3, d7 = 0x1.ba1ba1cdb8745p-5,
                          d9 = 0x1.664ed49cfc666p-6, d11 = 0x1.2385a3cf2e4eap-7, e0 = 0x1.5555555554dbdp-2,
                          e1 = 0x1.11112e0a6b45fp-3, g1 = 0x1.b096cp-27, g2 = 0x1.f212dp-5, g3 = 0x1.92f1ap-1,
                          g4 = 25.0, g5 = 1e8, mp1 = 0x1.921fb58p0, mp2 = -0x1.dde973cp-27,
                          mp3 = -0x1.cb3b399d747f2p-55, pp3 = -0x1.cb3b398p-55, pp4 = -0x1.d747f23e32ed7p-83,
                          hpinv = 0x1.45f306dc9c883p-1, toint = 0x1.8p52;
};
// the reduced argument a + da (quadrant parity n): tan (n = 0) or -cot (n = 1)
CPPWRF_GM double tan_reduced(double a, double da, int n) {
  using C = TanC;
  double ya, yya, sy;
  if (a < 0.0) {
    ya = -a;
    yya = -da;
    sy = -1;
  } else {
    ya = a;
    yya = da;
    sy = 1;
  }
  if (ya <= C::g2) {  // gy2
    const double a2 = a * a;
    double t2 = fma(a2, C::d11, C::d9);
    t2 = fma(a2, t2, C::d7);
    t2 = fma(a2, t2, C::d5);
    t2 = fma(a2, t2, C::d3);
    t2 = fma(a * a2, t2, da);
    if (n) {
      const double b = a + t2;  // EADD (a, t2, b, db)
      const double db = fabs_(a) > fabs_(t2) ? (a - b) + t2 : (t2 - b) + a;
      const double c = 1.0 / b;  // DIV2 (1.0, 0.0, b, db, c, dc, ...)
      const double u = c * b;
      const double uu = fma(c, b, -u);
      const double cc = fma(-db, c, ((1.0 - u) - uu) + 0.0) / b;
      const double z = c + cc;
      const double zz = (c - z) + cc;
      return -(zz + z);
    }
    return a + t2;
  }
  const int i = int(fma(ya, 256.0, -15.5));
  const double z = (ya - CPPWRF_GD(kTanXfg, 32 * i)) + yya;
  const double z2 = z * z;
  const double pz = fma(z * z2, fma(z2, C::e1, C::e0), z);
  const double fi = CPPWRF_GD(kTanXfg, 32 * i + 8), gi = CPPWRF_GD(kTanXfg, 32 * i + 16);
  const double num = (fi + gi) * pz;
  if (n) return (gi - num / (pz + fi)) * -sy;
  return (num / (gi - pz) + fi) * sy;
}
}  // namespace detail
CPPWRF_GM double tan(double x) {
  using namespace detail;
  using C = TanC;
  const std::uint32_t ux = std::uint32_t(asu(x) >> 32);
  if ((ux & 0x7ff00000u) == 0x7ff00000u) return x - x;
  const double w = x < 0.0 ? -x : x;
  if (w <= C::g1) return x;
  if (w <= C::g2) {
    const double x2 = x * x;
    double t2 = fma(x2, C::d11, C::d9);
    t2 = fma(x2, t2, C::d7);
    t2 = fma(x2, t2, C::d5);
    t2 = fma(x2, t2, C::d3);
    return fma(x * x2, t2, x);
  }
  if (w <= C::g3) {
    const int i = int(fma(w, 256.0, -15.5));
    const double z = w - CPPWRF_GD(kTanXfg, 32 * i);
    const double z2 = z * z;
    const double s = x < 0.0 ? -1.0 : 1.0;
    const double pz = fma(z * z2, fma(z2, C::e1, C::e0), z);
    const double fi = CPPWRF_GD(kTanXfg, 32 * i + 8), gi = CPPWRF_GD(kTanXfg, 32 * i + 16);
    return ((fi + gi) * pz / (gi - pz) + fi) * s;
  }
  double a, da;
  int n;
  if (w <= C::g4) {
    const double t = fma(x, C::hpinv, C::toint);
    const double xn = t - C::toint;
    double t1 = fma(-xn, C::mp1, x);
    t1 = fma(-xn, C::mp2, t1);
    n = int(asu(t) & 1);
    a = fma(-xn, C::mp3, t1);
    da = fma(-xn, C::mp3, t1 - a);
    return tan_reduced(a, da, n);
  }
  if (w <= C::g5) {
    const double t = fma(x, C::hpinv, C::toint);
    const double xn = t - C::toint;
    double t1 = fma(-xn, C::mp1, x);
    t1 = fma(-xn, C::mp2, t1);
    n = int(asu(t) & 1);
    const double tt = fma(-xn, C::pp3, t1);
    da = fma(-xn, C::pp3, t1 - tt);
    a = fma(-xn, C::pp4, tt);
    da = da + fma(-xn, C::pp4, tt - a);
  } else {
    n = branred(x, &a, &da) & 1;
  }
  const double z = a + da;  // EADD (a, da, t1, t2)
  const double zz = fabs_(a) > fabs_(da) ? (a - z) + da : (da - z) + a;
  return tan_reduced(z, zz, n);
}

// ---- atan2: __ieee754_atan2_fma (libm 0x7c170; sysdeps/ieee754/dbl-64/e_atan2.c, constants of atnat2.h, cij of
// uatan.tbl (the same as atan's); EMULV, EADD and ESUB of dla.h with FMA) ----
CPPWRF_GM double atan2(double y, double x) {
  using namespace detail;
  constexpr double d3 = -0x1.5555555555555p-2, d5 = 0x1.99999999997fdp-3, d7 = -0x1.24924923f7603p-3,
                   d9 = 0x1.c71c6e5129a3bp-4, d11 = -0x1.7458022b13c25p-4, d13 = 0x1.375f08b31cbcep-4;
  constexpr double hpi = 0x1.921fb54442d18p0, hpi1 = 0x1.1a62633145c07p-54, opi = 0x1.921fb54442d18p1,
                   opi1 = 0x1.1a62633145c07p-53, qpi = 0x1.921fb54442d18p-1, tqpi = 0x1.2d97c7f3321d2p1,
                   two500 = 0x1p500, twom500 = 0x1p-500, inv16 = 0.0625, TWO52 = 0x1p52;
#define CPPWRF_CIJ(i, j) CPPWRF_GD(kAtanCij, 8 * (7 * (i) + (j)))
  const std::uint64_t bx = asu(x), by = asu(y);
  const std::uint32_t ux = std::uint32_t(bx >> 32), dx = std::uint32_t(bx), uy = std::uint32_t(by >> 32),
                      dy = std::uint32_t(by);
  if ((ux & 0x7ff00000u) == 0x7ff00000u && ((ux & 0x000fffffu) | dx) != 0) return x + y;
  if ((uy & 0x7ff00000u) == 0x7ff00000u && ((uy & 0x000fffffu) | dy) != 0) return y + y;
  if (uy == 0x00000000u && dy == 0) return (ux & 0x80000000u) == 0 ? 0.0 : opi;
  if (uy == 0x80000000u && dy == 0) return (ux & 0x80000000u) == 0 ? -0.0 : -opi;
  if (x == 0) return (uy & 0x80000000u) == 0 ? hpi : -hpi;
  if (ux == 0x7ff00000u && dx == 0) {
    if (uy == 0x7ff00000u) {
      if (dy == 0) return qpi;
    } else if (uy == 0xfff00000u) {
      if (dy == 0) return -qpi;
    } else {
      return (uy & 0x80000000u) == 0 ? 0.0 : -0.0;
    }
  } else if (ux == 0xfff00000u && dx == 0) {
    if (uy == 0x7ff00000u) {
      if (dy == 0) return tqpi;
    } else if (uy == 0xfff00000u) {
      if (dy == 0) return -tqpi;
    } else {
      return (uy & 0x80000000u) == 0 ? opi : -opi;
    }
  }
  if (uy == 0x7ff00000u) {
    if (dy == 0) return hpi;
  } else if (uy == 0xfff00000u) {
    if (dy == 0) return -hpi;
  }
  double ax = x < 0 ? -x : x, ay = y < 0 ? -y : y;
  const std::int32_t de = std::int32_t(uy & 0x7ff00000u) - std::int32_t(ux & 0x7ff00000u);
  if (de >= 59768832) return y > 0 ? hpi : -hpi;
  if (de <= -59768832) {
    if (x > 0) return copysign(ay / ax, y);
    return y > 0 ? opi : -opi;
  }
  if (ax < twom500 || ay < twom500) {
    ax *= two500;
    ay *= two500;
  }
  if (ax > two500 || ay > two500) {
    ax *= twom500;
    ay *= twom500;
  }
  double u, du;
  if (ay < ax) {
    u = ay / ax;
    const double v = ax * u, vv = fma(ax, u, -v);
    du = ((ay - v) - vv) / ax;
  } else {
    u = ax / ay;
    const double v = ay * u, vv = fma(ay, u, -v);
    du = ((ax - v) - vv) / ay;
  }
  auto small_poly = [&](double& zz) {  // u v (d3 + v (d5 + ... + v d13)), v = u^2
    const double v = u * u;
    double q = fma(v, d13, d11);
    q = fma(v, q, d9);
    q = fma(v, q, d7);
    q = fma(v, q, d5);
    q = fma(v, q, d3);
    zz = (u * v) * q;
    return q;
  };
  auto index = [&]() { return int(fma(u, 256.0, TWO52) - TWO52) - 16; };
  auto tbl_poly = [&](int i, double v) {  // cij[i][2] + v (cij[i][3] + ... + v cij[i][6])
    double q = fma(v, CPPWRF_CIJ(i, 6), CPPWRF_CIJ(i, 5));
    q = fma(v, q, CPPWRF_CIJ(i, 4));
    q = fma(v, q, CPPWRF_CIJ(i, 3));
    return fma(v, q, CPPWRF_CIJ(i, 2));
  };
  if (x > 0) {
    if (ay < ax) {  // (i) atan(ay / ax)
      if (u < inv16) {
        const double v = u * u;
        double q = fma(v, d13, d11);
        q = fma(v, q, d9);
        q = fma(v, q, d7);
        q = fma(v, q, d5);
        q = fma(v, q, d3);
        return copysign(u + fma(u * v, q, du), y);
      }
      const int i = index();
      const double t3 = u - CPPWRF_CIJ(i, 0);
      const double v = du + t3;  // EADD (t3, du, v, dv)
      const double dv = fabs_(t3) > fabs_(du) ? (t3 - v) + du : (du - v) + t3;
      const double t2 = CPPWRF_CIJ(i, 2);
      double q = fma(v, CPPWRF_CIJ(i, 6), CPPWRF_CIJ(i, 5));
      q = fma(v, q, CPPWRF_CIJ(i, 4));
      q = fma(v, q, CPPWRF_CIJ(i, 3));
      const double zz = fma(v, t2, fma(dv, t2, (v * v) * q));
      return copysign(zz + CPPWRF_CIJ(i, 1), y);
    }
    if (u < inv16) {  // (ii) pi/2 - atan(ax / ay)
      double zz;
      small_poly(zz);
      const double t2 = hpi - u;  // ESUB (hpi, u, t2, cor)
      const double cor = hpi > fabs_(u) ? (hpi - t2) - u : hpi - (u + t2);
      const double t3 = ((cor + hpi1) - du) - zz;
      return copysign(t3 + t2, y);
    }
    const int i = index();
    const double v = (u - CPPWRF_CIJ(i, 0)) + du;
    const double zz = fma(-v, tbl_poly(i, v), hpi1);
    return copysign((hpi - CPPWRF_CIJ(i, 1)) + zz, y);
  }
  if (ax < ay) {  // (iii) pi/2 + atan(ax / ay)
    if (u < inv16) {
      double zz;
      small_poly(zz);
      const double t2 = u + hpi;  // EADD (hpi, u, t2, cor)
      const double cor = hpi > fabs_(u) ? (hpi - t2) + u : (u - t2) + hpi;
      const double t3 = ((cor + hpi1) + du) + zz;
      return copysign(t3 + t2, y);
    }
    const int i = index();
    const double v = (u - CPPWRF_CIJ(i, 0)) + du;
    const double zz = fma(v, tbl_poly(i, v), hpi1);
    return copysign((hpi + CPPWRF_CIJ(i, 1)) + zz, y);
  }
  if (u < inv16) {  // (iv) pi - atan(ax / ay)
    double zz;
    small_poly(zz);
    const double t2 = opi - u;  // ESUB (opi, u, t2, cor)
    const double cor = opi > fabs_(u) ? (opi - t2) - u : opi - (u + t2);
    const double t3 = ((cor + opi1) - du) - zz;
    return copysign(t3 + t2, y);
  }
  const int i = index();
  const double v = (u - CPPWRF_CIJ(i, 0)) + du;
  const double zz = fma(-v, tbl_poly(i, v), opi1);
  return copysign((opi - CPPWRF_CIJ(i, 1)) + zz, y);
#undef CPPWRF_CIJ
}

// ---- tgamma: __tgamma (libm 0x16f80, w_tgamma_compat.c) -> __ieee754_gamma_r (0x52b00, e_gamma_r.c) with
// gamma_positive, __ieee754_lgamma_r (0x559a0, e_lgamma_r.c), __ieee754_exp2 (0x495b0, e_exp2.c), __scalbn and the
// x86-64 __gamma_product (0x57b00, x87 long double), all compiled without FMA; they call the FMA variants of exp,
// log, pow, expm1, sin and cos. mul_split etc. of math/mul_split.h without FMA (Dekker). ----
namespace detail {
struct LgC {
  static constexpr double
    two52=  4.50359962737049600000e+15,
    half=  5.00000000000000000000e-01,
    one =  1.00000000000000000000e+00,
    pi  =  3.14159265358979311600e+00,
    a0  =  7.72156649015328655494e-02,
    a1  =  3.22467033424113591611e-01,
    a2  =  6.73523010531292681824e-02,
    a3  =  2.05808084325167332806e-02,
    a4  =  7.38555086081402883957e-03,
    a5  =  2.89051383673415629091e-03,
    a6  =  1.19270763183362067845e-03,
    a7  =  5.10069792153511336608e-04,
    a8  =  2.20862790713908385557e-04,
    a9  =  1.08011567247583939954e-04,
    a10 =  2.52144565451257326939e-05,
    a11 =  4.48640949618915160150e-05,
    tc  =  1.46163214496836224576e+00,
    tf  = -1.21486290535849611461e-01,
    tt  = -3.63867699703950536541e-18,
    t0  =  4.83836122723810047042e-01,
    t1  = -1.47587722994593911752e-01,
    t2  =  6.46249402391333854778e-02,
    t3  = -3.27885410759859649565e-02,
    t4  =  1.79706750811820387126e-02,
    t5  = -1.03142241298341437450e-02,
    t6  =  6.10053870246291332635e-03,
    t7  = -3.68452016781138256760e-03,
    t8  =  2.25964780900612472250e-03,
    t9  = -1.40346469989232843813e-03,
    t10 =  8.81081882437654011382e-04,
    t11 = -5.38595305356740546715e-04,
    t12 =  3.15632070903625950361e-04,
    t13 = -3.12754168375120860518e-04,
    t14 =  3.35529192635519073543e-04,
    u0  = -7.72156649015328655494e-02,
    u1  =  6.32827064025093366517e-01,
    u2  =  1.45492250137234768737e+00,
    u3  =  9.77717527963372745603e-01,
    u4  =  2.28963728064692451092e-01,
    u5  =  1.33810918536787660377e-02,
    v1  =  2.45597793713041134822e+00,
    v2  =  2.12848976379893395361e+00,
    v3  =  7.69285150456672783825e-01,
    v4  =  1.04222645593369134254e-01,
    v5  =  3.21709242282423911810e-03,
    s0  = -7.72156649015328655494e-02,
    s1  =  2.14982415960608852501e-01,
    s2  =  3.25778796408930981787e-01,
    s3  =  1.46350472652464452805e-01,
    s4  =  2.66422703033638609560e-02,
    s5  =  1.84028451407337715652e-03,
    s6  =  3.19475326584100867617e-05,
    r1  =  1.39200533467621045958e+00,
    r2  =  7.21935547567138069525e-01,
    r3  =  1.71933865632803078993e-01,
    r4  =  1.86459191715652901344e-02,
    r5  =  7.77942496381893596434e-04,
    r6  =  7.32668430744625636189e-06,
    w0  =  4.18938533204672725052e-01,
    w1  =  8.33333333333329678849e-02,
    w2  = -2.77777777728775536470e-03,
    w3  =  7.93650558643019558500e-04,
    w4  = -5.95187557450339963135e-04,
    w5  =  8.36339918996282139126e-04,
    w6  = -1.63092934096575273989e-03;
};
// __ieee754_lgamma_r for x > 0 (the arguments gamma_positive passes)
CPPWRF_GM double lgamma_pos(double x) {
  using C = LgC;
  const std::uint64_t b = asu(x);
  const std::int32_t hx = std::int32_t(b >> 32), lx = std::int32_t(std::uint32_t(b));
  const std::int32_t ix = hx & 0x7fffffff;
  double t, y, z, p, p1, p2, p3, q, r, w;
  int i;
  if (ix >= 0x7ff00000) return x * x;
  if ((ix | lx) == 0) return C::one / fabs_(x);
  if (ix < 0x3b900000) return -log(x);
  if ((((ix - 0x3ff00000) | lx) == 0) || (((ix - 0x40000000) | lx) == 0)) {
    r = 0;
  } else if (ix < 0x40000000) {
    if (ix <= 0x3feccccc) {
      r = -log(x);
      if (ix >= 0x3FE76944) {
        y = C::one - x;
        i = 0;
      } else if (ix >= 0x3FCDA661) {
        y = x - (C::tc - C::one);
        i = 1;
      } else {
        y = x;
        i = 2;
      }
    } else {
      r = 0.0;
      if (ix >= 0x3FFBB4C3) {
        y = 2.0 - x;
        i = 0;
      } else if (ix >= 0x3FF3B4C4) {
        y = x - C::tc;
        i = 1;
      } else {
        y = x - C::one;
        i = 2;
      }
    }
    switch (i) {
      case 0:
        z = y * y;
        p1 = C::a0 + z * (C::a2 + z * (C::a4 + z * (C::a6 + z * (C::a8 + z * C::a10))));
        p2 = z * (C::a1 + z * (C::a3 + z * (C::a5 + z * (C::a7 + z * (C::a9 + z * C::a11)))));
        p = y * p1 + p2;
        r += (p - 0.5 * y);
        break;
      case 1:
        z = y * y;
        w = z * y;
        p1 = C::t0 + w * (C::t3 + w * (C::t6 + w * (C::t9 + w * C::t12)));
        p2 = C::t1 + w * (C::t4 + w * (C::t7 + w * (C::t10 + w * C::t13)));
        p3 = C::t2 + w * (C::t5 + w * (C::t8 + w * (C::t11 + w * C::t14)));
        p = z * p1 - (C::tt - w * (p2 + y * p3));
        r += (C::tf + p);
        break;
      default:
        p1 = y * (C::u0 + y * (C::u1 + y * (C::u2 + y * (C::u3 + y * (C::u4 + y * C::u5)))));
        p2 = C::one + y * (C::v1 + y * (C::v2 + y * (C::v3 + y * (C::v4 + y * C::v5))));
        r += (-0.5 * y + p1 / p2);
    }
  } else if (ix < 0x40200000) {
    i = int(x);
    y = x - double(i);
    p = y * (C::s0 + y * (C::s1 + y * (C::s2 + y * (C::s3 + y * (C::s4 + y * (C::s5 + y * C::s6))))));
    q = C::one + y * (C::r1 + y * (C::r2 + y * (C::r3 + y * (C::r4 + y * (C::r5 + y * C::r6)))));
    r = C::half * y + p / q;
    z = C::one;
    switch (i) {
      case 7: z *= (y + 6.0); [[fallthrough]];
      case 6: z *= (y + 5.0); [[fallthrough]];
      case 5: z *= (y + 4.0); [[fallthrough]];
      case 4: z *= (y + 3.0); [[fallthrough]];
      case 3: z *= (y + 2.0); r += log(z); break;
      default: break;
    }
  } else if (ix < 0x43900000) {
    t = log(x);
    z = C::one / x;
    y = z * z;
    w = C::w0 + z * (C::w1 + y * (C::w2 + y * (C::w3 + y * (C::w4 + y * (C::w5 + y * C::w6)))));
    r = (x - C::half) * (t - C::one) + w;
  } else {
    r = x * (log(x) - C::one);
  }
  return r;
}
// __ieee754_exp2 (e_exp2.c; __exp_data: exp2_shift +0x40, exp2_poly C1 ... C5 +0x48 ... +0x68, tab +0xb0)
CPPWRF_GM double exp2(double x) {
  std::uint32_t abstop = std::uint32_t(asu(x) >> 52) & 0x7ff;
  if (abstop - 0x3c9u >= 0x408u - 0x3c9u) {
    if (abstop - 0x3c9u >= 0x80000000u) return 1.0 + x;
    if (abstop >= 0x409u) {
      if (asu(x) == 0xfff0000000000000ull) return 0.0;
      if (abstop >= 0x7ffu) return 1.0 + x;
      if (!(asu(x) >> 63)) return math_oflow(0);
      if (asu(x) >= asu(-1075.0)) return math_uflow(0);
    }
    if (2 * asu(x) > 2 * asu(928.0)) abstop = 0;
  }
  const double shift = CPPWRF_GD(kExpData, 0x40);
  double kd = x + shift;
  const std::uint64_t ki = asu(kd);
  kd -= shift;
  const double r = x - kd;
  const std::uint64_t idx = 2 * (ki % 128);
  const std::uint64_t top = ki << (52 - 7);
  const double tail = asd(CPPWRF_GU(kExpData, 0xb0 + 8 * idx));
  std::uint64_t sbits = CPPWRF_GU(kExpData, 0xb0 + 8 * (idx + 1)) + top;
  const double r2 = r * r;
  const double tmp = tail + r * CPPWRF_GD(kExpData, 0x48) + r2 * (CPPWRF_GD(kExpData, 0x50) + r * CPPWRF_GD(kExpData, 0x58)) +
                     r2 * r2 * (CPPWRF_GD(kExpData, 0x60) + r * CPPWRF_GD(kExpData, 0x68));
  if (abstop == 0) {  // specialcase
    if ((ki & 0x80000000u) == 0) {
      sbits -= 1ull << 52;
      const double scale = asd(sbits);
      return 2 * (scale + scale * tmp);
    }
    sbits += 1022ull << 52;
    const double scale = asd(sbits);
    double y = scale + scale * tmp;
    if (y < 1.0) {
      double lo = scale - y + scale * tmp;
      const double hi = 1.0 + y;
      lo = 1.0 - hi + y + lo;
      y = (hi + lo) - 1.0;
      if (y == 0.0) y = 0.0;
    }
    return 0x1p-1022 * y;
  }
  const double scale = asd(sbits);
  return scale + scale * tmp;
}
CPPWRF_GM double scalbn(double x, int n) {
  std::int64_t ix = std::int64_t(asu(x));
  std::int64_t k = (ix >> 52) & 0x7ff;
  if (k == 0) {
    if ((ix & 0xfffffffffffffll) == 0) return x;
    x *= 0x1p54;
    ix = std::int64_t(asu(x));
    k = ((ix >> 52) & 0x7ff) - 54;
  }
  if (k == 0x7ff) return x + x;
  if (n < -50000) return 1.0e-300 * copysign(1.0e-300, x);
  if (n > 50000 || k + n > 0x7fe) return 1.0e+300 * copysign(1.0e+300, x);
  k = k + n;
  if (k > 0) return asd(std::uint64_t((ix & std::int64_t(0x800fffffffffffffull)) | (k << 52)));
  if (k <= -54) return 1.0e-300 * copysign(1.0e-300, x);
  k += 54;
  x = asd(std::uint64_t((ix & std::int64_t(0x800fffffffffffffull)) | (k << 52)));
  return x * 0x1p-54;
}
// x87 extended precision (64-bit significand, round to nearest even; Linux's default precision control) as
// __gamma_product computes in it: value (-1)^neg m 2^(e - 63), m normalised (bit 63 set) or zero
struct X87 {
  std::uint64_t m;
  int e;
  bool neg;
};
CPPWRF_GM int clz128(unsigned __int128 v) {
  const std::uint64_t hi = std::uint64_t(v >> 64), lo = std::uint64_t(v);
#if defined(__CUDA_ARCH__)
  return hi ? __clzll((long long)hi) : 64 + __clzll((long long)lo);
#else
  return hi ? __builtin_clzll(hi) : 64 + __builtin_clzll(lo);
#endif
}
// v 2^E (v > 0, any bits; bit 0 may be a sticky bit) rounded to 64 bits
CPPWRF_GM X87 x87_round(bool neg, unsigned __int128 v, int E) {
  const int s = clz128(v);
  v <<= s;
  E -= s;
  std::uint64_t m = std::uint64_t(v >> 64);
  const std::uint64_t rest = std::uint64_t(v);
  if (rest > 0x8000000000000000ull || (rest == 0x8000000000000000ull && (m & 1))) {
    if (++m == 0) {
      m = 0x8000000000000000ull;
      ++E;
    }
  }
  return X87{m, E + 127, neg};
}
CPPWRF_GM X87 x87_from(double d) {
  const std::uint64_t b = asu(d);
  const bool neg = b >> 63;
  const int ex = int((b >> 52) & 0x7ff);
  std::uint64_t f = b & 0xfffffffffffffull;
  if (ex == 0 && f == 0) return X87{0, 0, neg};
  int e = ex - 1023;
  if (ex == 0) e = -1022;
  else f |= 1ull << 52;
  return x87_round(neg, (unsigned __int128)f, e - 52);
}
CPPWRF_GM X87 x87_add(X87 a, X87 b) {
  if (b.m == 0) return a.m == 0 ? X87{0, 0, a.neg && b.neg} : a;
  if (a.m == 0) return b;
  if (a.e < b.e || (a.e == b.e && a.m < b.m)) {
    const X87 t = a;
    a = b;
    b = t;
  }
  const unsigned __int128 A = (unsigned __int128)a.m << 62;
  const int d = a.e - b.e;
  unsigned __int128 B;
  if (d >= 126) {
    B = 1;  // sticky
  } else {
    const unsigned __int128 full = (unsigned __int128)b.m << 62;
    B = full >> d;
    if ((B << d) != full) B |= 1;
  }
  const unsigned __int128 S = a.neg == b.neg ? A + B : A - B;
  if (S == 0) return X87{0, 0, false};
  return x87_round(a.neg, S, a.e - 63 - 62);
}
CPPWRF_GM X87 x87_mul(X87 a, X87 b) {
  if (a.m == 0 || b.m == 0) return X87{0, 0, a.neg != b.neg};
  return x87_round(a.neg != b.neg, (unsigned __int128)a.m * b.m, a.e + b.e - 126);
}
CPPWRF_GM X87 x87_div(X87 a, X87 b) {
  if (a.m == 0) return X87{0, 0, a.neg != b.neg};
  const unsigned __int128 n1 = (unsigned __int128)a.m << 63;
  const std::uint64_t q1 = std::uint64_t(n1 / b.m), r1 = std::uint64_t(n1 % b.m);
  const unsigned __int128 n2 = (unsigned __int128)r1 << 64;
  const std::uint64_t q2 = std::uint64_t(n2 / b.m), r2 = std::uint64_t(n2 % b.m);
  unsigned __int128 Q = ((unsigned __int128)q1 << 64) | q2;
  if (r2) Q |= 1;
  return x87_round(a.neg != b.neg, Q, a.e - b.e - 127);
}
// fstl: to double, round to nearest even (normal results)
CPPWRF_GM double x87_to_double(X87 a) {
  if (a.m == 0) return a.neg ? -0.0 : 0.0;
  std::uint64_t m = a.m >> 11;
  const std::uint64_t rest = a.m & 0x7ff;
  int e = a.e;
  if (rest > 0x400 || (rest == 0x400 && (m & 1))) {
    if (++m == (1ull << 53)) {
      m >>= 1;
      ++e;
    }
  }
  return asd((std::uint64_t(a.neg) << 63) | (std::uint64_t(e + 1023) << 52) | (m & 0xfffffffffffffull));
}
// __gamma_product (sysdeps/x86_64/fpu/gamma_product.c: long double): x + x_eps times x + x_eps + i, i < n
CPPWRF_GM double gamma_product(double x, double x_eps, int n, double* eps) {
  const X87 s = x87_add(x87_from(x), x87_from(x_eps));
  X87 ret = s;
  for (int i = 1; i < n; i++) ret = x87_mul(ret, x87_add(x87_from(double(i)), s));
  const double d = x87_to_double(ret);
  const X87 dd = x87_from(d);
  *eps = x87_to_double(x87_div(x87_add(ret, X87{dd.m, dd.e, !dd.neg}), dd));
  return d;
}
CPPWRF_GM void mul_split(double* hi, double* lo, double x, double y) {
  *hi = x * y;
  double x1 = x * 134217729.0;
  double y1 = y * 134217729.0;
  x1 = (x - x1) + x1;
  y1 = (y - y1) + y1;
  const double x2 = x - x1;
  const double y2 = y - y1;
  *lo = (((x1 * y1 - *hi) + x1 * y2) + x2 * y1) + x2 * y2;
}
CPPWRF_GM void fast_two_sum(double* hi, double* lo, double a, double b) {
  *hi = a + b;
  const double e = *hi - a;
  *lo = b - e;
}
CPPWRF_GM void mul_expansion(double* hi, double* lo, double h1, double l1, double h2, double l2) {
  double e;
  mul_split(hi, lo, h1, h2);
  const double r = h1 * l2 + h2 * l1;
  fast_two_sum(hi, &e, *hi, r);
  *lo -= e;
}
CPPWRF_GM void div_expansion(double* hi, double* lo, double h1, double l1, double h2, double l2) {
  double e, a, b;
  *hi = h1 / h2;  // div_split
  mul_split(&a, &b, *hi, h2);
  a = h1 - a;
  a = a - b;
  *lo = a / h2;
  const double r = (l1 * h2 - l2 * h1) / (h2 * h2);
  fast_two_sum(hi, &e, *hi, r);
  *lo += e;
  fast_two_sum(hi, lo, *hi, *lo);
}
CPPWRF_GM double gamma_positive(double x, int* exp2_adj) {
  if (x < 0.5) {
    *exp2_adj = 0;
    return exp(lgamma_pos(x + 1)) / x;
  }
  if (x <= 1.5) {
    *exp2_adj = 0;
    return exp(lgamma_pos(x));
  }
  if (x < 6.5) {
    *exp2_adj = 0;
    const double n = std::ceil(x - 1.5);
    const double x_adj = x - n;
    double eps;
    const double prod = gamma_product(x_adj, 0, int(n), &eps);
    return exp(lgamma_pos(x_adj)) * prod * (1.0 + eps);
  }
  double eps = 0, x_eps = 0, x_adj = x, prod = 1;
  if (x < 12.0) {
    const double n = std::ceil(12.0 - x);
    x_adj = x + n;
    x_eps = (x - (x_adj - n));
    prod = gamma_product(x_adj - n, x_eps, int(n), &eps);
  }
  const double x_adj_int = std::round(x_adj);
  const double x_adj_frac = x_adj - x_adj_int;
  int x_adj_log2;
  double x_adj_mant = std::frexp(x_adj, &x_adj_log2);
  if (x_adj_mant < 0.70710678118654752440) {
    x_adj_log2--;
    x_adj_mant *= 2.0;
  }
  *exp2_adj = x_adj_log2 * int(x_adj_int);
  double h1, l1, h2, l2;
  mul_split(&h1, &l1, pow(x_adj_mant, x_adj), exp2(x_adj_log2 * x_adj_frac));
  mul_split(&h2, &l2, exp(-x_adj), std::sqrt(2 * 3.14159265358979323846 / x_adj));
  mul_expansion(&h1, &l1, h1, l1, h2, l2);
  div_expansion(&h1, &l1, h1, l1, prod, prod * eps);
  double exp_adj = x_eps * log(x_adj);
  constexpr double gc[6] = {0x1.5555555555555p-4, -0xb.60b60b60b60b8p-12, 0x3.4034034034034p-12,
                            -0x2.7027027027028p-12, 0x3.72a3c5631fe46p-12, -0x7.daac36664f1f4p-12};
  double bsum = gc[5];
  const double x_adj2 = x_adj * x_adj;
  for (int i = 1; i <= 5; i++) bsum = bsum / x_adj2 + gc[5 - i];
  exp_adj += bsum / x_adj;
  l1 += h1 * expm1(exp_adj);
  return h1 + l1;
}
}  // namespace detail
CPPWRF_GM double tgamma(double x) {
  using namespace detail;
  const std::uint64_t b = asu(x);
  const std::int32_t hx = std::int32_t(b >> 32);
  const std::uint32_t lx = std::uint32_t(b);
  int signgam;
  double ret;
  if (((hx & 0x7fffffff) | lx) == 0) return 1.0 / x;
  if (hx < 0 && std::uint32_t(hx) < 0xfff00000u && rint(x) == x) return (x - x) / (x - x);
  if (std::uint32_t(hx) == 0xfff00000u && lx == 0) return x - x;
  if ((hx & 0x7ff00000) == 0x7ff00000) return x + x;
  if (x >= 172.0) return 0x1.fffffffffffffp1023 * 0x1.fffffffffffffp1023;
  if (x > 0.0) {
    signgam = 0;
    int exp2_adj;
    const double tret = gamma_positive(x, &exp2_adj);
    ret = scalbn(tret, exp2_adj);
  } else if (x >= -0x1p-52 / 4.0) {
    signgam = 0;
    ret = 1.0 / x;
  } else {
    const double tx = std::trunc(x);
    signgam = (tx == 2.0 * std::trunc(tx / 2.0)) ? -1 : 1;
    if (x <= -184.0) {
      ret = 0x1p-1022 * 0x1p-1022;
    } else {
      double frac = tx - x;
      if (frac > 0.5) frac = 1.0 - frac;
      const double sinpix = frac <= 0.25 ? sin(3.14159265358979323846 * frac) : cos(3.14159265358979323846 * (0.5 - frac));
      int exp2_adj;
      double h1, l1, h2, l2;
      h2 = gamma_positive(-x, &exp2_adj);
      mul_split(&h1, &l1, sinpix, h2);
      mul_split(&h2, &l2, h1, x);
      l2 += l1 * x;
      h1 = 0x3.243f6a8885a3p+0;
      l1 = 0x8.d313198a2e038p-56;
      div_expansion(&h1, &l1, h1, l1, h2, l2);
      ret = scalbn(-h1, -exp2_adj);
    }
  }
  if ((asu(ret) & 0x7fffffffffffffffull) == 0x7ff0000000000000ull && x != 0) {
    if (signgam < 0) ret = -(-copysign(0x1.fffffffffffffp1023, ret) * 0x1.fffffffffffffp1023);
    else ret = copysign(0x1.fffffffffffffp1023, ret) * 0x1.fffffffffffffp1023;
  } else if (ret == 0) {
    if (signgam < 0) ret = -(-copysign(0x1p-1022, ret) * 0x1p-1022);
    else ret = copysign(0x1p-1022, ret) * 0x1p-1022;
  }
  return signgam < 0 ? -ret : ret;
}

namespace detail {
// fpclassify / isfinite / isinf / signbit from the bits (device and host alike)
enum Fp { kFpNan, kFpInf, kFpZero, kFpSub, kFpNormal };
CPPWRF_GM int fpcls(double x) {
  const std::uint64_t a = asu(x) & 0x7fffffffffffffffull;
  if (a > 0x7ff0000000000000ull) return kFpNan;
  if (a == 0x7ff0000000000000ull) return kFpInf;
  if (a == 0) return kFpZero;
  return a < 0x0010000000000000ull ? kFpSub : kFpNormal;
}
CPPWRF_GM bool finite_(double x) { return (asu(x) & 0x7fffffffffffffffull) < 0x7ff0000000000000ull; }
CPPWRF_GM bool inf_(double x) { return (asu(x) & 0x7fffffffffffffffull) == 0x7ff0000000000000ull; }
CPPWRF_GM bool signbit_(double x) { return (asu(x) >> 63) != 0; }
}  // namespace detail

// ---- hypot: __hypot (libm 0x532b0; sysdeps/ieee754/dbl-64/e_hypot.c without FMA: no ifunc variant) ----
// ax = maxsd, ay = minsd of |x|, |y|; Borges' correction of sqrt(ax * ax + ay * ay)
namespace detail {
CPPWRF_GM double hypot_kernel(double ax, double ay) {
  double h = std::sqrt(ax * ax + ay * ay);
  double t1, t2;
  if (h <= 2.0 * ay) {
    const double delta = h - ay;
    t1 = ax * (2.0 * delta - ax);
    t2 = (delta - 2.0 * (ax - ay)) * delta;
  } else {
    const double delta = h - ax;
    t1 = 2.0 * delta * (ax - 2.0 * ay);
    t2 = (4.0 * delta - ay) * ay + delta * delta;
  }
  h -= (t1 + t2) / (2.0 * h);
  return h;
}
}  // namespace detail
CPPWRF_GM double hypot(double x, double y) {
  using namespace detail;
  if (!finite_(x) || !finite_(y)) {
    if ((inf_(x) || inf_(y)) && !issignaling(x) && !issignaling(y)) return asd(0x7ff0000000000000ull);
    return x + y;
  }
  x = fabs_(x);
  y = fabs_(y);
  const double ax = y > x ? y : x;  // maxsd |x| into |y|
  const double ay = x < y ? x : y;  // minsd |y| into |x|
  if (ax > 0x1p+511) {
    if (ay <= ax * 0x1p-54) return ax + ay;
    return hypot_kernel(ax * 0x1p-600, ay * 0x1p-600) * 0x1p+600;
  }
  if (ay < 0x1p-459) {
    if (ax >= ay * 0x1p+54) return ax + ay;
    return hypot_kernel(ax * 0x1p+600, ay * 0x1p+600) * 0x1p-600;
  }
  if (ay <= ax * 0x1p-54) return ax + ay;
  return hypot_kernel(ax, ay);
}

// ---- log1p: __log1p_fma (libm 0x80290; sysdeps/ieee754/dbl-64/s_log1p.c built with FMA) ----
CPPWRF_GM double log1p(double x) {
  using namespace detail;
  constexpr double ln2_hi = 0x1.62e42fee00000p-1, ln2_lo = 0x1.a39ef35793c76p-33, two54 = 0x1p54;
  constexpr double Lp1 = 0x1.5555555555593p-1, Lp2 = 0x1.999999997fa04p-2, Lp3 = 0x1.2492494229359p-2,
                   Lp4 = 0x1.c71c51d8e78afp-3, Lp5 = 0x1.7466496cb03dep-3, Lp6 = 0x1.39a09d078c69fp-3,
                   Lp7 = 0x1.2f112df3e5244p-3;
  const std::int32_t hx = std::int32_t(asu(x) >> 32);
  const std::int32_t ax = hx & 0x7fffffff;
  int k = 1;
  double f = 0, c = 0;
  std::int32_t hu = 0;
  if (hx < 0x3FDA827A) {
    if (ax >= 0x3ff00000) {
      if (x == -1.0) return -two54 / 0.0;
      return (x - x) / (x - x);
    }
    if (ax < 0x3e200000) {
      if (ax < 0x3c900000) return x;
      return fma(-(x * x), 0.5, x);  // vfnmadd231sd
    }
    if (hx > 0 || hx <= std::int32_t(0xbfd2bec3)) {
      k = 0;
      f = x;
      hu = 1;
    }
  } else if (hx >= 0x7ff00000) {
    return x + x;
  }
  double u;
  if (k != 0) {
    if (hx < 0x43400000) {
      u = x + 1.0;
      hu = std::int32_t(asu(u) >> 32);
      k = (hu >> 20) - 1023;
      c = (k > 0) ? 1.0 - (u - x) : x - (u - 1.0);
      c /= u;
    } else {
      u = x;
      hu = hx;
      k = (hu >> 20) - 1023;
      c = 0;
    }
    hu &= 0x000fffff;
    if (hu < 0x6a09e) {
      u = with_high(u, std::uint32_t(hu | 0x3ff00000));
    } else {
      k += 1;
      u = with_high(u, std::uint32_t(hu | 0x3fe00000));
      hu = (0x00100000 - hu) >> 2;
    }
    f = u - 1.0;
  }
  const double hfsq = (f * 0.5) * f;
  if (hu == 0) {
    if (f == 0.0) {
      if (k == 0) return 0.0;
      const double dk = double(k);
      c = fma(dk, ln2_lo, c);
      return fma(dk, ln2_hi, c);
    }
    const double R = fma(-f, 0x1.5555555555555p-1, 1.0) * hfsq;
    if (k == 0) return f - R;
    const double dk = double(k);
    return fma(dk, ln2_hi, -((R - fma(dk, ln2_lo, c)) - f));
  }
  const double s = f / (f + 2.0);
  const double z = s * s;
  const double R2 = fma(z, Lp3, Lp2), R3 = fma(z, Lp5, Lp4), R4 = fma(z, Lp7, Lp6);
  const double z2 = z * z, z4 = z2 * z2, z6 = z2 * z4;
  double R = fma(z, Lp1, z2 * R2);
  R = fma(z4, R3, R);
  R = fma(z6, R4, R);
  const double t = (R + hfsq) * s;
  if (k == 0) return f - (hfsq - t);
  const double dk = double(k);
  return fma(dk, ln2_hi, -((hfsq - (fma(dk, ln2_lo, c) + t)) - f));
}

// ---- x2y2m1: __x2y2m1 (libm 0x5edc0; sysdeps/ieee754/dbl-64/x2y2m1.c, mul_split without FMA) ----
// x^2 + y^2 - 1 for 1 > x >= y >= epsilon / 2, x^2 + y^2 >= 0.5. qsort is glibc's stable merge sort: any stable
// sort by |value| gives its order.
namespace detail {
CPPWRF_GM void sort_abs(double* v, int n) {  // stable insertion sort by |v| ascending (qsort with x2y2m1's compare)
  for (int i = 1; i < n; ++i) {
    const double t = v[i];
    int j = i - 1;
    while (j >= 0 && fabs_(v[j]) > fabs_(t)) {
      v[j + 1] = v[j];
      --j;
    }
    v[j + 1] = t;
  }
}
}  // namespace detail
CPPWRF_GM double x2y2m1(double x, double y) {
  using namespace detail;
  double vals[5];
  mul_split(&vals[1], &vals[0], x, x);
  mul_split(&vals[3], &vals[2], y, y);
  vals[4] = -1.0;
  sort_abs(vals, 5);
  for (int i = 0; i <= 3; i++) {
    const double a = vals[i + 1], b = vals[i];
    const double hi = b + a;
    vals[i] = (a - hi) + b;
    vals[i + 1] = hi;
    sort_abs(vals + i + 1, 4 - i);
  }
  return vals[4] + vals[3] + vals[2] + vals[1] + vals[0];
}

// ---- clog: __clog (libm 0x4ca70; math/s_clog_template.c), real part to *re, imaginary part to *im ----
CPPWRF_GM void clog(double xr, double xi, double* re, double* im) {
  using namespace detail;
  constexpr double kMax = 0x1.fffffffffffffp1023, kMin = 0x1p-1022, kEps = 0x1p-52;
  const int rcls = fpcls(xr), icls = fpcls(xi);
  if (rcls == kFpZero && icls == kFpZero) {
    *im = copysign(signbit_(xr) ? 0x1.921fb54442d18p+1 : 0.0, xi);
    *re = -1.0 / fabs_(xr);
    return;
  }
  if (rcls == kFpNan || icls == kFpNan) {
    *im = asd(0x7ff8000000000000ull);
    *re = (rcls == kFpInf || icls == kFpInf) ? asd(0x7ff0000000000000ull) : asd(0x7ff8000000000000ull);
    return;
  }
  double absx = fabs_(xr), absy = fabs_(xi);
  int scale = 0;
  if (absx < absy) {
    const double t = absx;
    absx = absy;
    absy = t;
  }
  if (absx > kMax / 2) {
    scale = -1;
    absx = scalbn(absx, scale);
    absy = (absy >= kMin * 2 ? scalbn(absy, scale) : 0);
  } else if (absx < kMin && absy < kMin) {
    scale = 53;
    absx = scalbn(absx, scale);
    absy = scalbn(absy, scale);
  }
  if (absx == 1 && scale == 0) {
    *re = log1p(absy * absy) * 0.5;
  } else if (absx > 1 && absx < 2 && absy < 1 && scale == 0) {
    double d2m1 = (absx - 1) * (absx + 1);
    if (absy >= kEps) d2m1 += absy * absy;
    *re = log1p(d2m1) * 0.5;
  } else if (absx < 1 && absx >= 0.5 && absy < kEps / 2 && scale == 0) {
    const double d2m1 = (absx - 1) * (absx + 1);
    *re = log1p(d2m1) * 0.5;
  } else if (absx < 1 && absx >= 0.5 && scale == 0 && absx * absx + absy * absy >= 0.5) {
    *re = log1p(x2y2m1(absx, absy)) * 0.5;
  } else {
    const double d = hypot(absx, absy);
    *re = log(d) - double(scale) * 0x1.62e42fefa39efp-1;
  }
  *im = atan2(xi, xr);
}

// ---- csqrt: __csqrt (libm 0x4daa0; math/s_csqrt_template.c), real part to *re, imaginary part to *im ----
CPPWRF_GM void csqrt(double xr, double xi, double* re, double* im) {
  using namespace detail;
  constexpr double kMax = 0x1.fffffffffffffp1023, kMin = 0x1p-1022;
  const double kInf = asd(0x7ff0000000000000ull), kNan = asd(0x7ff8000000000000ull);
  const int rcls = fpcls(xr), icls = fpcls(xi);
  if (rcls == kFpNan || rcls == kFpInf || icls == kFpNan || icls == kFpInf) {
    if (icls == kFpInf) {
      *re = kInf;
      *im = xi;
    } else if (rcls == kFpInf) {
      if (xr < 0) {
        *re = icls == kFpNan ? kNan : 0;
        *im = copysign(kInf, xi);
      } else {
        *re = xr;
        *im = icls == kFpNan ? kNan : copysign(0.0, xi);
      }
    } else {
      *re = kNan;
      *im = kNan;
    }
    return;
  }
  if (icls == kFpZero) {
    if (xr < 0) {
      *re = 0;
      *im = copysign(std::sqrt(-xr), xi);
    } else {
      *re = fabs_(std::sqrt(xr));
      *im = copysign(0.0, xi);
    }
    return;
  }
  if (rcls == kFpZero) {
    double r;
    if (fabs_(xi) >= 2 * kMin) r = std::sqrt(0.5 * fabs_(xi));
    else r = 0.5 * std::sqrt(2 * fabs_(xi));
    *re = r;
    *im = copysign(r, xi);
    return;
  }
  double d, r, s;
  int scale = 0;
  if (fabs_(xr) > kMax / 4) {
    scale = 1;
    xr = scalbn(xr, -2 * scale);
    xi = scalbn(xi, -2 * scale);
  } else if (fabs_(xi) > kMax / 4) {
    scale = 1;
    if (fabs_(xr) >= 4 * kMin) xr = scalbn(xr, -2 * scale);
    else xr = 0;
    xi = scalbn(xi, -2 * scale);
  } else if (fabs_(xr) < 2 * kMin && fabs_(xi) < 2 * kMin) {
    scale = -((53 + 1) / 2);
    xr = scalbn(xr, -2 * scale);
    xi = scalbn(xi, -2 * scale);
  }
  d = hypot(xr, xi);
  if (xr > 0) {
    r = std::sqrt(0.5 * (d + xr));
    if (scale == 1 && fabs_(xi) < 1) {
      s = xi / r;
      r = scalbn(r, scale);
      scale = 0;
    } else {
      s = 0.5 * (xi / r);
    }
  } else {
    s = std::sqrt(0.5 * (d - xr));
    if (scale == 1 && fabs_(xi) < 1) {
      r = fabs_(xi / s);
      s = scalbn(s, scale);
      scale = 0;
    } else {
      r = fabs_(0.5 * (xi / s));
    }
  }
  if (scale) {
    r = scalbn(r, scale);
    s = scalbn(s, scale);
  }
  *re = r;
  *im = copysign(s, xi);
}

// ---- libmvec, one lane of the 2-lane SSE4 functions ----
namespace mvec {
// _ZGVbN2v_exp_sse4 (libmvec 0xba40; __svml_dexp_data: dbT +0, InvLn2 +8192, Shifter +8256, Ln2hi +8320,
// Ln2lo +8384, PC0 +8448, PC1 +8512, PC2 +8576, lIndexMask +8640, iAbsMask +8704, iDomainRange +8768).
// Lanes with |x| beyond the domain (and NaN) call the scalar exp.
CPPWRF_GM double exp(double x) {
  using namespace detail;
  const std::int32_t hx = std::int32_t(std::uint32_t(asu(x) >> 32) & std::uint32_t(CPPWRF_GU(kSvmlExpData, 8704)));
  if (hx > std::int32_t(std::uint32_t(CPPWRF_GU(kSvmlExpData, 8768)))) return glibc::exp(x);
  const double dk = x * CPPWRF_GD(kSvmlExpData, 8192);
  const double dn = rint(dk);
  const double nh = CPPWRF_GD(kSvmlExpData, 8320) * dn;
  const double nl = dn * CPPWRF_GD(kSvmlExpData, 8384);
  const double dm = CPPWRF_GD(kSvmlExpData, 8256) + dk;
  double r = x - nh;
  r = r - nl;
  double p = CPPWRF_GD(kSvmlExpData, 8576) * r;
  p = p + CPPWRF_GD(kSvmlExpData, 8512);
  p = p * r;
  p = p + CPPWRF_GD(kSvmlExpData, 8448);
  const double rp = r * p;
  const std::uint64_t mask = CPPWRF_GU(kSvmlExpData, 8640);
  const std::uint64_t lindex = asu(dm) & mask;
  const std::uint64_t lm = (asu(dm) & ~mask) << 42;
  const double e = CPPWRF_GD(kSvmlExpData, 8448) + rp;
  const double res = CPPWRF_GD(kSvmlExpData, 8 * lindex) * e;
  return asd(asu(res) + lm);
}

// _ZGVbN2vv_pow_sse4 (libmvec 0xca50; svml_d_pow2_core_sse4.S, offsets of svml_d_pow_data.h). Lanes outside the
// ranges (x subnormal / huge / non-positive / NaN, |y| huge / NaN, |y log2 x| beyond the domain) call the scalar pow.
CPPWRF_GM double pow(double x, double y) {
  using namespace detail;
  constexpr int iIndexMask = 40384, iIndexAdd = 40448, i3fe7fe00 = 40512, i2p20_2p19 = 40576, iOne = 40640,
                jIndexMask = 40704, iMantissaMask = 39232, dbOne = 39360, iffffffff00000000 = 39424,
                db2p20_2p19 = 39488, iHighMask = 39552, LHN = 39616, ifff0000000000000 = 39680, db2p45_2p44 = 39744,
                HIDELTA = 40064, LORANGE = 40128, ABSMASK = 40192, INF = 40256, DOMAINRANGE = 40320, rcp_t1 = 11712,
                log2_t1 = 19968, exp2_tbl = 36416, clv_2 = 38528, cev_1 = 38912;
  auto lo32 = [](std::uint64_t w) { return std::uint32_t(w); };
  const std::uint32_t hx = std::uint32_t(asu(x) >> 32);
  // i = (((Hi(x) & 0x000ffe00) + 0x00000200) >> 10)
  const std::uint32_t ii = ((hx & lo32(CPPWRF_GU(kSvmlPowData, iIndexMask))) + lo32(CPPWRF_GU(kSvmlPowData, iIndexAdd))) >> 10;
  const int ircp = std::int32_t(ii << 3), ilog = std::int32_t(ii << 4);
  // x1 = mantissa of x with exponent 0
  const double x1 = asd((asu(x) & CPPWRF_GU(kSvmlPowData, iMantissaMask)) | CPPWRF_GU(kSvmlPowData, dbOne));
  // k = (Hi(x) - 0x3fe7fe00) >> 20 (arithmetic), + 0x41380000 as the high word of a double, - 2^20 - 2^19
  std::int32_t k = std::int32_t(hx - lo32(CPPWRF_GU(kSvmlPowData, i3fe7fe00))) >> 20;
  k = std::int32_t(std::uint32_t(k) + lo32(CPPWRF_GU(kSvmlPowData, i2p20_2p19)));
  const std::uint64_t hm = CPPWRF_GU(kSvmlPowData, iHighMask);
  const double rcp1 = CPPWRF_GD(kSvmlPowData, rcp_t1 + ircp);
  const double x1hi = asd(asu(x1) & hm);
  const double x1lo = x1 - x1hi;
  const double r1 = x1 * rcp1;
  double E = x1hi * rcp1;
  const double xlr = x1lo * rcp1;
  E = E - r1;
  double kd = asd((std::uint64_t(std::uint32_t(k)) << 32) & CPPWRF_GU(kSvmlPowData, iffffffff00000000));
  kd = kd - CPPWRF_GD(kSvmlPowData, db2p20_2p19);
  E = E + xlr;
  const double L1hi = CPPWRF_GD(kSvmlPowData, log2_t1 + ilog), L1lo = CPPWRF_GD(kSvmlPowData, log2_t1 + ilog + 8);
  double cq = CPPWRF_GD(kSvmlPowData, LHN) + r1;
  const double T = L1hi + kd;
  const double T_Rh = cq + T;
  const double T_Rh_Eh = E + T_Rh;
  double Rl = T - T_Rh;
  double HLL = T_Rh - T_Rh_Eh;
  Rl = Rl + cq;
  cq = cq + E;
  HLL = HLL + E;
  const double X = Rl + HLL;
  double hll = L1lo + X;
  double p = CPPWRF_GD(kSvmlPowData, clv_2) * cq;
  for (int c = 1; c < 6; ++c) {  // + clv_3 ... clv_7
    p = p + CPPWRF_GD(kSvmlPowData, clv_2 + 64 * c);
    if (c < 5) p = p * cq;
  }
  const double cqp = cq * p;
  hll = hll + cqp;
  const double T_Rh_Eh_HLLhi = T_Rh_Eh + hll;
  // 2^(y*(HH+HL+HLL))
  const double yH = asd(asu(y) & hm);
  const double yL = y - yH;
  const double HH = asd(asu(T_Rh_Eh_HLLhi) & hm);
  const double HLLhi = T_Rh_Eh_HLLhi - T_Rh_Eh;
  const double HL = T_Rh_Eh_HLLhi - HH;
  const std::int32_t hxd = std::int32_t(hx + lo32(CPPWRF_GU(kSvmlPowData, HIDELTA)));
  bool special = std::int32_t(lo32(CPPWRF_GU(kSvmlPowData, LORANGE))) > hxd;
  const double pH = yH * HH;
  double pL = yL * HL;
  const double yhhl = yH * HL;
  const double hhyl = HH * yL;
  pL = pL + yhhl;
  hll = hll - HLLhi;
  pL = pL + hhyl;
  const double pLL = hll * y;
  double pHH = CPPWRF_GD(kSvmlPowData, db2p45_2p44) + pH;
  const double t = pL + pLL;
  const std::uint32_t absm = lo32(CPPWRF_GU(kSvmlPowData, ABSMASK));
  const std::int32_t hy = std::int32_t(std::uint32_t(asu(y) >> 32) & absm), inf = std::int32_t(lo32(CPPWRF_GU(kSvmlPowData, INF)));
  special = special || hy > inf || hy == inf;
  const std::uint32_t jlo = lo32(asu(pHH));
  pHH = pHH - CPPWRF_GD(kSvmlPowData, db2p45_2p44);
  const std::int32_t hph = std::int32_t(std::uint32_t(asu(pH) >> 32) & absm),
                     dr = std::int32_t(lo32(CPPWRF_GU(kSvmlPowData, DOMAINRANGE)));
  double Z = pH - pHH;
  Z = Z + t;
  special = special || hph > dr || hph == dr;
  if (special) return glibc::pow(x, y);
  const int j = std::int32_t((jlo & lo32(CPPWRF_GU(kSvmlPowData, jIndexMask))) << 4);
  const std::uint32_t n = (jlo << 13) + lo32(CPPWRF_GU(kSvmlPowData, iOne));
  const double two_n = asd((std::uint64_t(n) << 32) & CPPWRF_GU(kSvmlPowData, ifff0000000000000));
  double q = CPPWRF_GD(kSvmlPowData, cev_1) * Z;
  q = q + CPPWRF_GD(kSvmlPowData, cev_1 + 64);
  q = q * Z;
  q = q + CPPWRF_GD(kSvmlPowData, cev_1 + 128);
  q = q * Z;
  q = q + CPPWRF_GD(kSvmlPowData, cev_1 + 192);
  const double res_hi = two_n * CPPWRF_GD(kSvmlPowData, exp2_tbl + j);
  double zq = Z * q;
  const double zr = Z * res_hi;
  zq = zq + CPPWRF_GD(kSvmlPowData, cev_1 + 256);
  zq = zq * zr;
  return res_hi + zq;
}

// _ZGVbN2v_sin_sse4 (libmvec 0xcf20; svml_d_sin2_core_sse4.S, __svml_d_trig_data of svml_d_trig_data.h: AbsMask +0,
// RangeVal +64, InvPI +256, RShifter +320, PI1 ... PI4 +576 ... +768, C1_sin ... C7_sin +1728 ... +2112). Lanes with
// |x| > RangeVal (and NaN) call the scalar sin.
CPPWRF_GM double sin(double x) {
  using namespace detail;
#define CPPWRF_T(off) CPPWRF_GD(kSvmlTrigData, off)
  const double ax = fabs_(x);
  if (!(ax <= CPPWRF_T(64))) return ::cppwrf::glibc::sin(x);
  const double y = ax * CPPWRF_T(256) + CPPWRF_T(320);
  const double n = y - CPPWRF_T(320);
  const std::uint64_t sign_res = asu(y) << 63;
  double r = ax - n * CPPWRF_T(576);
  r = r - n * CPPWRF_T(640);
  r = r - n * CPPWRF_T(704);
  r = r - n * CPPWRF_T(768);
  double r2 = r * r;
  r = asd(asu(r) ^ sign_res);
  double p = CPPWRF_T(2112) * r2 + CPPWRF_T(2048);
  p = p * r2 + CPPWRF_T(1984);
  p = p * r2 + CPPWRF_T(1920);
  p = p * r2 + CPPWRF_T(1856);
  p = p * r2 + CPPWRF_T(1792);
  p = p * r2 + CPPWRF_T(1728);
  r2 = r2 * p;
  r2 = r2 * r;
  r = r + r2;
  return asd(asu(r) ^ (asu(x) & 0x8000000000000000ull));
}
// _ZGVbN2v_cos_sse4 (libmvec 0xb260; svml_d_cos2_core_sse4.S: HalfPI +192, OneHalf +512, C1 ... C7 +1280 ... +1664).
// Lanes with |x + pi/2| > RangeVal (and NaN) call the scalar cos.
CPPWRF_GM double cos(double x) {
  using namespace detail;
  const double xp = CPPWRF_T(192) + x;
  if (!(fabs_(xp) <= CPPWRF_T(64))) return ::cppwrf::glibc::cos(x);
  const double y = xp * CPPWRF_T(256) + CPPWRF_T(320);
  const std::uint64_t sign_res = asu(y) << 63;
  double n = y - CPPWRF_T(320);
  n = n - CPPWRF_T(512);
  double r = x - n * CPPWRF_T(576);
  r = r - n * CPPWRF_T(640);
  r = r - n * CPPWRF_T(704);
  r = r - n * CPPWRF_T(768);
  double r2 = r * r;
  double p = CPPWRF_T(1664) * r2 + CPPWRF_T(1600);
  p = p * r2 + CPPWRF_T(1536);
  p = p * r2 + CPPWRF_T(1472);
  p = p * r2 + CPPWRF_T(1408);
  p = p * r2 + CPPWRF_T(1344);
  p = p * r2 + CPPWRF_T(1280);
  r2 = r2 * p;
  r2 = r2 * r;
  r = r + r2;
  return asd(asu(r) ^ sign_res);
#undef CPPWRF_T
}
// _ZGVbN2v_atan_sse4 (libmvec 0xa8e0; svml_d_atan2_core_sse4.S, __svml_datan_data_internal_avx512: Shifter +16,
// MaxThreshold +32, MOne +48, One +64, LargeX +80, Tbl_H +112, dIndexMed +624, Pi2 +640, coeff +672 ... +752,
// .FLT_11 index mask +768). No scalar path.
CPPWRF_GM double atan(double x) {
  using namespace detail;
#define CPPWRF_T(off) CPPWRF_GD(kSvmlAtanData, off)
  const double ax = fabs_(x);
  const std::uint64_t sign = asu(x) ^ asu(ax);
  const double y = CPPWRF_T(16) + ax;
  double xr = y - CPPWRF_T(16);
  const double diff = ax - xr;
  xr = xr * ax;
  const bool large = CPPWRF_T(32) <= ax;
  xr = xr + CPPWRF_T(64);
  const double lx = CPPWRF_T(80) < ax ? CPPWRF_T(80) : ax;  // minpd: the second operand unless the first is less
  const double num = large ? CPPWRF_T(48) : diff;
  const double den = large ? lx : xr;
  double r = num / den;
  const int idx = int((asu(y) << 3) & CPPWRF_GU(kSvmlAtanData, 768));
  const double r2 = r * r;
  const double r4 = r2 * r2, r3 = r2 * r;
  double p5 = CPPWRF_T(672) * r2 + CPPWRF_T(688);
  double p6 = CPPWRF_T(704) * r2;
  p5 = p5 * r4;
  p6 = p6 + CPPWRF_T(720);
  p6 = p6 + p5;
  double p8 = CPPWRF_T(736) * r2;
  p6 = r4 * p6;
  p8 = p8 + CPPWRF_T(752);
  p8 = p8 + p6;
  r = r + r3 * p8;
  const bool med = CPPWRF_T(624) <= y;
  const double tbl = med ? CPPWRF_T(112 + 128 + idx) : CPPWRF_T(112 + idx);
  const double t = large ? CPPWRF_T(640) : tbl;
  return asd(asu(t + r) ^ sign);
#undef CPPWRF_T
}
}  // namespace mvec
}  // namespace cppwrf::glibc
#if defined(__CUDACC__)
#pragma nv_diagnostic pop
#endif
