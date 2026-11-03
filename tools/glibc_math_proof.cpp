// cppwrf::glibc (src/core/glibc_math.hpp) against this machine's glibc libm and libmvec, bit for bit. Run where WRF's
// references are made (Debian 13, glibc 2.41, x86-64 with FMA and AVX2: libm's FMA variants), e.g. in the Debian 13
// chroot of a Colab runtime (tools/colab/remote_chroot_job.sh). The inputs (edge values and four random sets per function) are
// those of tools/glibc_math_inputs.hpp; the vector functions get two inputs per call (one per lane). NaN results are
// equal when both are NaN.
//
//   g++ -std=c++20 -O2 -march=x86-64-v3 -ffp-contract=off -fopenmp -Isrc -Itools tools/glibc_math_proof.cpp
//       -o glibc_math_proof -lmvec -lm && ./glibc_math_proof [N] [FUNCTION,...]
// N random inputs per function (default 10^9). Prints per function the inputs compared and those that differ (the
// first ones with their bits); exit 1 when any differs.
#include <emmintrin.h>
#include <omp.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "glibc_math_inputs.hpp"

extern "C" __m128d _ZGVbN2v_exp(__m128d);
extern "C" __m128d _ZGVbN2vv_pow(__m128d, __m128d);
extern "C" __m128d _ZGVbN2v_log(__m128d);
extern "C" __m128d _ZGVbN2v_log10(__m128d);
extern "C" __m128d _ZGVbN2v_sin(__m128d);
extern "C" __m128d _ZGVbN2v_cos(__m128d);
extern "C" __m128d _ZGVbN2v_atan(__m128d);
extern "C" double _Complex clog(double _Complex);  // C99 (g++ extension _Complex)
extern "C" double _Complex csqrt(double _Complex);

using namespace glibc_math_inputs;

namespace {
// libm / libmvec through pointers (no compile-time folding)
double (*volatile c_exp)(double) = ::exp;
double (*volatile c_log)(double) = ::log;
double (*volatile c_pow)(double, double) = ::pow;
double (*volatile c_log10)(double) = ::log10;
double (*volatile c_atan)(double) = ::atan;
double (*volatile c_expm1)(double) = ::expm1;
double (*volatile c_tanh)(double) = ::tanh;
double (*volatile c_sin)(double) = ::sin;
double (*volatile c_cos)(double) = ::cos;
double (*volatile c_sinh)(double) = ::sinh;
double (*volatile c_erf)(double) = ::erf;
double (*volatile c_asin)(double) = ::asin;
double (*volatile c_acos)(double) = ::acos;
double (*volatile c_tan)(double) = ::tan;
double (*volatile c_atan2)(double, double) = ::atan2;
double (*volatile c_tgamma)(double) = ::tgamma;
double (*volatile c_hypot)(double, double) = ::hypot;
double (*volatile c_log1p)(double) = ::log1p;
double _Complex (*volatile c_clog)(double _Complex) = ::clog;
double _Complex (*volatile c_csqrt)(double _Complex) = ::csqrt;
double _Complex cx(double x, double y) {
  double _Complex z;
  __real__ z = x;
  __imag__ z = y;
  return z;
}

// glibc's result for function f on the lane pair (x[0], y[0]), (x[1], y[1])
void ref(int f, const double* x, const double* y, double* r) {
  switch (f) {
    case kExp: for (int l = 0; l < 2; ++l) r[l] = c_exp(x[l]); break;
    case kLog: for (int l = 0; l < 2; ++l) r[l] = c_log(x[l]); break;
    case kPow: for (int l = 0; l < 2; ++l) r[l] = c_pow(x[l], y[l]); break;
    case kVExp: _mm_storeu_pd(r, _ZGVbN2v_exp(_mm_loadu_pd(x))); break;
    case kVPow: _mm_storeu_pd(r, _ZGVbN2vv_pow(_mm_loadu_pd(x), _mm_loadu_pd(y))); break;
    case kLog10: for (int l = 0; l < 2; ++l) r[l] = c_log10(x[l]); break;
    case kAtan: for (int l = 0; l < 2; ++l) r[l] = c_atan(x[l]); break;
    case kExpm1: for (int l = 0; l < 2; ++l) r[l] = c_expm1(x[l]); break;
    case kTanh: for (int l = 0; l < 2; ++l) r[l] = c_tanh(x[l]); break;
    case kVLog: _mm_storeu_pd(r, _ZGVbN2v_log(_mm_loadu_pd(x))); break;
    case kVLog10: _mm_storeu_pd(r, _ZGVbN2v_log10(_mm_loadu_pd(x))); break;
    case kSin: for (int l = 0; l < 2; ++l) r[l] = c_sin(x[l]); break;
    case kCos: for (int l = 0; l < 2; ++l) r[l] = c_cos(x[l]); break;
    case kVSin: _mm_storeu_pd(r, _ZGVbN2v_sin(_mm_loadu_pd(x))); break;
    case kVCos: _mm_storeu_pd(r, _ZGVbN2v_cos(_mm_loadu_pd(x))); break;
    case kVAtan: _mm_storeu_pd(r, _ZGVbN2v_atan(_mm_loadu_pd(x))); break;
    case kSinh: for (int l = 0; l < 2; ++l) r[l] = c_sinh(x[l]); break;
    case kErf: for (int l = 0; l < 2; ++l) r[l] = c_erf(x[l]); break;
    case kAsin: for (int l = 0; l < 2; ++l) r[l] = c_asin(x[l]); break;
    case kAcos: for (int l = 0; l < 2; ++l) r[l] = c_acos(x[l]); break;
    case kTan: for (int l = 0; l < 2; ++l) r[l] = c_tan(x[l]); break;
    case kAtan2: for (int l = 0; l < 2; ++l) r[l] = c_atan2(x[l], y[l]); break;
    case kTgamma: for (int l = 0; l < 2; ++l) r[l] = c_tgamma(x[l]); break;
    case kHypot: for (int l = 0; l < 2; ++l) r[l] = c_hypot(x[l], y[l]); break;
    case kLog1p: for (int l = 0; l < 2; ++l) r[l] = c_log1p(x[l]); break;
    case kClogRe: for (int l = 0; l < 2; ++l) r[l] = __real__ c_clog(cx(x[l], y[l])); break;
    case kClogIm: for (int l = 0; l < 2; ++l) r[l] = __imag__ c_clog(cx(x[l], y[l])); break;
    case kCsqrtRe: for (int l = 0; l < 2; ++l) r[l] = __real__ c_csqrt(cx(x[l], y[l])); break;
    case kCsqrtIm: for (int l = 0; l < 2; ++l) r[l] = __imag__ c_csqrt(cx(x[l], y[l])); break;
  }
}
template <int... F>
double port(int f, double x, double y, std::integer_sequence<int, F...>) {
  double r = 0;
  ((f == F ? (void)(r = dispatch<F>(x, y)) : (void)0), ...);
  return r;
}

struct Count {
  long long n = 0, bad = 0;
};

void pair(int f, const double* x, const double* y, Count& c, int& shown) {
  double r[2];
  ref(f, x, y, r);
  for (int l = 0; l < 2; ++l) {
    const double p = port(f, x[l], y[l], std::make_integer_sequence<int, kNumFuncs>{});
    c.n++;
    if (!same(r[l], p)) {
      c.bad++;
#pragma omp critical
      if (shown++ < 10)
        std::printf("  %s(%016llx, %016llx) = %016llx (%.17g), port %016llx (%.17g)\n", kNames[f],
                    (unsigned long long)B(x[l]), (unsigned long long)B(y[l]), (unsigned long long)B(r[l]), r[l],
                    (unsigned long long)B(p), p);
    }
  }
}
}  // namespace

int main(int argc, char** argv) {
  const long long n_random = argc > 1 ? std::atoll(argv[1]) : 1000000000ll;
  const std::string only = argc > 2 ? argv[2] : "";
  std::printf("glibc_math_proof: %lld random inputs per function, %d threads\n", n_random, omp_get_max_threads());
  const std::vector<double> edges = edge_values();
  long long total_bad = 0;
  for (int f = 0; f < kNumFuncs; ++f) {
    if (!only.empty() && ("," + only + ",").find("," + std::string(kNames[f]) + ",") == std::string::npos) continue;
    int shown = 0;
    const auto t0 = std::chrono::steady_clock::now();
    const long long ne = (long long)edges.size(), n_edges = kArity[f] == 1 ? ne : ne * 64;
    Count ce, cr;
#pragma omp parallel
    {
      Count le, lr;
      std::uint64_t s = 0x1234567ull + 977ull * omp_get_thread_num() + 131ull * f;
#pragma omp for schedule(static)
      for (long long t = 0; t < (n_edges + 1) / 2; ++t) {
        double x[2], y[2];
        for (int l = 0; l < 2; ++l) edge_input(f, edges, std::min(2 * t + l, n_edges - 1), s, x[l], y[l]);
        pair(f, x, y, le, shown);
      }
#pragma omp for schedule(static)
      for (long long t = 0; t < n_random / 2; ++t) {
        double x[2], y[2];
        for (int l = 0; l < 2; ++l) gen(f, s, int((2 * t + l) % 4), x[l], y[l]);
        pair(f, x, y, lr, shown);
      }
#pragma omp critical
      {
        ce.n += le.n; ce.bad += le.bad;
        cr.n += lr.n; cr.bad += lr.bad;
      }
    }
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("%-10s edges %10lld inputs, %lld differ; random %12lld inputs, %lld differ (%.0f s)\n", kNames[f], ce.n,
                ce.bad, cr.n, cr.bad, sec);
    std::fflush(stdout);
    total_bad += ce.bad + cr.bad;
  }
  std::printf(total_bad ? "glibc_math_proof: %lld DIFFER\n" : "glibc_math_proof: identical\n", total_bad);
  return total_bad ? 1 : 0;
}
