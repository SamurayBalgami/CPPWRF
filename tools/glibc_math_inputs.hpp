#pragma once
// The functions and inputs of tools/glibc_math_proof.cpp (CPU: port against glibc) and tools/glibc_math_proof_cuda.cu
// (device against host). Per function the edge values (zeros, every power of two and its neighbours, subnormals, the
// largest values, infinities, NaNs of both signs with payloads, 1 and its neighbours, the thresholds of the
// algorithms; two-argument functions: every edge value with 64 partners, edge values and random ones) and four
// random sets, a quarter of the inputs each: raw 64-bit patterns, log-uniform values over the whole range, the
// function's typical range, and the values WRF computes.
#include <cmath>
#include <cstdint>
#include <cstring>
#include <utility>
#include <vector>

#include "core/glibc_math.hpp"
#include "core/libmvec_x86.hpp"

namespace glibc_math_inputs {
inline std::uint64_t B(double v) {
  std::uint64_t b;
  std::memcpy(&b, &v, 8);
  return b;
}
inline double D(std::uint64_t b) {
  double v;
  std::memcpy(&v, &b, 8);
  return v;
}
inline bool same(double a, double b) { return B(a) == B(b) || (std::isnan(a) && std::isnan(b)); }
inline std::uint64_t splitmix(std::uint64_t& s) {
  std::uint64_t z = (s += 0x9e3779b97f4a7c15ull);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
  return z ^ (z >> 31);
}
inline double unif(std::uint64_t& s, double a, double b) { return a + (b - a) * (double(splitmix(s) >> 11) * 0x1p-53); }
inline double logu(std::uint64_t& s, double a, double b) { return std::exp(unif(s, std::log(a), std::log(b))); }
inline double sgn(std::uint64_t& s) { return (splitmix(s) & 1) ? -1.0 : 1.0; }

// the functions (index F): name, arity, libmvec lane function
enum : int { kExp, kLog, kPow, kVExp, kVPow, kLog10, kAtan, kExpm1, kTanh, kVLog, kVLog10, kSin, kCos, kVSin, kVCos, kVAtan, kSinh, kErf, kAsin, kAcos, kTan, kAtan2, kTgamma, kHypot, kLog1p, kClogRe, kClogIm, kCsqrtRe, kCsqrtIm, kNumFuncs };
inline constexpr const char* kNames[kNumFuncs] = {"exp",   "log",  "pow",   "mvec_exp",  "mvec_pow",  "log10",
                                                  "atan",  "expm1", "tanh", "mvec_log", "mvec_log10", "sin",
                                                  "cos",   "mvec_sin", "mvec_cos", "mvec_atan", "sinh",
                                                  "erf",   "asin", "acos", "tan",   "atan2", "tgamma",
                                                  "hypot", "log1p", "clog_re", "clog_im", "csqrt_re", "csqrt_im"};
inline constexpr int kArity[kNumFuncs] = {1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1,
                                          2, 1, 2, 2, 2, 2};

template <int F>
CPPWRF_GM double dispatch(double x, double y) {
  namespace g = cppwrf::glibc;
  (void)y;
  if constexpr (F == kExp) return g::exp(x);
  else if constexpr (F == kLog) return g::log(x);
  else if constexpr (F == kPow) return g::pow(x, y);
  else if constexpr (F == kVExp) return g::mvec::exp(x);
  else if constexpr (F == kVPow) return g::mvec::pow(x, y);
  else if constexpr (F == kLog10) return g::log10(x);
  else if constexpr (F == kAtan) return g::atan(x);
  else if constexpr (F == kExpm1) return g::expm1(x);
  else if constexpr (F == kTanh) return g::tanh(x);
  else if constexpr (F == kVLog) return cppwrf::mvec::log(x);
  else if constexpr (F == kVLog10) return cppwrf::mvec::log10(x);
  else if constexpr (F == kSin) return g::sin(x);
  else if constexpr (F == kCos) return g::cos(x);
  else if constexpr (F == kVSin) return g::mvec::sin(x);
  else if constexpr (F == kVCos) return g::mvec::cos(x);
  else if constexpr (F == kVAtan) return g::mvec::atan(x);
  else if constexpr (F == kSinh) return g::sinh(x);
  else if constexpr (F == kErf) return g::erf(x);
  else if constexpr (F == kAsin) return g::asin(x);
  else if constexpr (F == kAcos) return g::acos(x);
  else if constexpr (F == kTan) return g::tan(x);
  else if constexpr (F == kAtan2) return g::atan2(x, y);
  else if constexpr (F == kTgamma) return g::tgamma(x);
  else if constexpr (F == kHypot) return g::hypot(x, y);
  else if constexpr (F == kLog1p) return g::log1p(x);
  else if constexpr (F == kClogRe || F == kClogIm || F == kCsqrtRe || F == kCsqrtIm) {  // z = x + i y
    double re, im;
    if constexpr (F == kClogRe || F == kClogIm) g::clog(x, y, &re, &im);
    else g::csqrt(x, y, &re, &im);
    return (F == kClogRe || F == kCsqrtRe) ? re : im;
  } else return 0;
}

inline constexpr double kWrfExp[] = {0.2857142857142857, 287.0 / 1004.0, 1004.0 / 717.0, 0.5, 1.5, 1.0 / 3.0,
                                     0.33333333, 2.0 / 3.0, 0.25, 0.75, 1.6, -1.6364, 0.16, 0.19, 4.0 / 3.0, 2.5, 0.6,
                                     0.86, -0.5, 1.0 / 7.0, 0.1, 0.4, 2.0, 3.0, -2.0, 1.0, -1.0, 0.0, 0.81, 2.2};

// a random input of set `set` (0 ... 3) for function f
inline void gen(int f, std::uint64_t& s, int set, double& x, double& y) {
  y = 0.0;
  if (set == 0) {  // raw bit patterns
    x = D(splitmix(s));
    y = D(splitmix(s));
    return;
  }
  switch (f) {
    case kExp:
    case kVExp:
      if (set == 1) x = sgn(s) * logu(s, 1e-300, 800.);
      else if (set == 2) x = unif(s, -746., 710.);
      else x = unif(s, -40., 40.);
      return;
    case kLog:
    case kLog10:
    case kVLog:
    case kVLog10:
      if (set == 1) x = logu(s, 1e-320, 1e308);
      else if (set == 2) x = unif(s, 0.9, 1.1);
      else x = unif(s, 1e-3, 3.);
      return;
    case kPow:
    case kVPow:
      if (set == 1) {
        x = logu(s, 1e-320, 1e308);
        y = sgn(s) * logu(s, 1e-30, 1e4);
      } else if (set == 2) {
        if (splitmix(s) & 1) {  // negative x, integer y
          x = -logu(s, 1e-5, 1e5);
          y = double(std::int64_t(splitmix(s) % 61) - 30);
        } else {
          x = unif(s, 0.5, 2.0);
          y = unif(s, -4., 4.);
        }
      } else {
        x = unif(s, 0.05, 2.0);
        y = kWrfExp[splitmix(s) % (sizeof kWrfExp / sizeof kWrfExp[0])];
        if (splitmix(s) % 4 == 0) y = unif(s, 0., 3.);
      }
      return;
    case kAtan:
    case kVAtan:
      if (set == 1) x = sgn(s) * logu(s, 1e-300, 1e300);
      else if (set == 2) x = unif(s, -20., 20.);
      else x = sgn(s) * logu(s, 1e-3, 1e3);  // sfclayrev / MYNN stability functions, Noah-MP
      return;
    case kExpm1:
    case kTanh:
    case kSinh:
      if (set == 1) x = sgn(s) * logu(s, 1e-300, 800.);
      else if (set == 2) x = unif(s, -2., 2.);
      else x = unif(s, -25., 25.);
      return;
    case kSin:
    case kCos:
    case kVSin:
    case kVCos:
    case kTan:
      if (set == 1) x = sgn(s) * logu(s, 1e-300, 1e300);
      else if (set == 2) x = unif(s, -4., 4.);  // the polynomial, the pi/2 - |x| and the reduced ranges
      else x = sgn(s) * logu(s, 1e-3, 2e8);     // reduce_sincos and __branred (|x| >= 105414350)
      return;
    case kErf:
      if (set == 1) x = sgn(s) * logu(s, 1e-300, 1e3);
      else if (set == 2) x = unif(s, -3., 3.);
      else x = unif(s, -7., 7.);
      return;
    case kAsin:
    case kAcos:
      if (set == 1) x = sgn(s) * logu(s, 1e-300, 1.0);
      else if (set == 2) x = unif(s, -1., 1.);
      else x = sgn(s) * unif(s, 0.9, 1.0);  // the near-1 branches
      return;
    case kAtan2:  // atan2(x, y)
      if (set == 1) {
        x = sgn(s) * logu(s, 1e-300, 1e300);
        y = sgn(s) * logu(s, 1e-300, 1e300);
      } else if (set == 2) {
        x = unif(s, -10., 10.);
        y = unif(s, -10., 10.);
      } else {
        x = sgn(s) * logu(s, 1e-3, 1e3);
        y = sgn(s) * logu(s, 1e-3, 1e3);
      }
      return;
    case kTgamma:
      if (set == 1) x = sgn(s) * logu(s, 1e-300, 200.);
      else if (set == 2) x = unif(s, -10., 13.);  // every branch of gamma_positive, the negative side
      else x = unif(s, 0.3, 7.);                  // Grell-Freitas' alpha, beta and their sums
      return;
    case kHypot:
    case kClogRe:
    case kClogIm:
    case kCsqrtRe:
    case kCsqrtIm:  // x, y (z = x + i y); set 3: |z| near 1 (clog's log1p / x2y2m1 branches), the radar's values
      if (set == 1) {
        x = sgn(s) * logu(s, 1e-300, 1e300);
        y = sgn(s) * logu(s, 1e-300, 1e300);
      } else if (set == 2) {
        x = unif(s, -100., 100.);
        y = unif(s, -100., 100.);
      } else if (splitmix(s) & 1) {
        const double r = 1.0 + unif(s, -0.6, 0.6) * (splitmix(s) & 1 ? 1.0 : 1e-8), t = unif(s, -3.2, 3.2);
        x = r * std::cos(t);
        y = r * std::sin(t);
        if (splitmix(s) % 4 == 0) y = sgn(s) * logu(s, 1e-30, 1e-8);
      } else {
        x = sgn(s) * logu(s, 1e-3, 1e2);
        y = sgn(s) * logu(s, 1e-6, 1e2);
      }
      return;
    case kLog1p:
      if (set == 1) x = logu(s, 1e-300, 1e308) * (splitmix(s) % 4 == 0 ? -1e-300 : 1.0);
      else if (set == 2) x = unif(s, -0.999, 1.0);
      else x = sgn(s) * logu(s, 1e-30, 0.5);
      return;
    default: x = 0.0; return;
  }
}

inline std::vector<double> edge_values() {
  std::vector<double> e;
  for (int k = -1074; k <= 1023; ++k) {
    const double p = std::ldexp(1.0, k);
    for (double v : {p, std::nextafter(p, 0.0), std::nextafter(p, INFINITY)}) {
      e.push_back(v);
      e.push_back(-v);
    }
  }
  for (std::uint64_t b : {0x0ull, 0x8000000000000000ull, 0x000fffffffffffffull, 0x7fefffffffffffffull,
                          0x7ff0000000000000ull, 0xfff0000000000000ull, 0x7ff8000000000000ull, 0xfff8000000000000ull,
                          0x7ff0000000000001ull, 0x7ff4000000000000ull, 0xfff0000000000001ull, 0x7fffffffffffffffull,
                          0xffffffffffffffffull, 0x7ff8000000000123ull})
    e.push_back(D(b));
  for (double v : {1.0, 0.9375, 1.0 + 0x1.09p-4, 0x1p-54, 512.0, 708.0, 709.78, 709.79, -708.4, -708.5, -744.4,
                   -745.1, -745.2, 0.5, 2.0, 10.0, 0.1, 3.0, 0.2857142857142857, 1e-300, 1e300, 0.7853981633974483,
                   1.5707963267948966, 3.141592653589793, 6.283185307179586, 1e22, 1e6, 105414350.0})
    for (double w : {v, std::nextafter(v, 0.0), std::nextafter(v, INFINITY), -v}) e.push_back(w);
  for (int k = -1100; k <= 1100; k += 7) e.push_back(double(k));
  return e;
}
// edge input u (0 ... n_edges - 1; n_edges = size for one argument, size * 64 for two)
inline void edge_input(int f, const std::vector<double>& e, long long u, std::uint64_t& s, double& x, double& y) {
  const long long ne = (long long)e.size();
  if (kArity[f] == 1) {
    x = e[u];
    y = 0.0;
    return;
  }
  x = e[u % ne];
  const long long v = u / ne;
  if (v < 32) y = e[splitmix(s) % ne];
  else {
    double dx;
    gen(f, s, int(v % 4), dx, y);
  }
  if (v % 3 == 0) std::swap(x, y);
}
}  // namespace glibc_math_inputs
