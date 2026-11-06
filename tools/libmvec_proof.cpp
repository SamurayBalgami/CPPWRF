// SPDX-License-Identifier: GPL-3.0-or-later
// cppwrf::mvec::log / log10 (src/core/libmvec_x86.hpp) against this machine's glibc libmvec _ZGVbN2v_log /
// _ZGVbN2v_log10, bit for bit (run on Intel: there they must agree everywhere). Three input sets per function:
//   edges       0, -0, every power of two (subnormal to the largest), its neighbours, DBL_MIN / DBL_MAX and theirs,
//               infinities, NaNs (quiet, signalling, negative, payloads), negatives, 1 and its neighbours, 10^k;
//   classes     every biased exponent (0 ... 2047) and sign x every rcpps input class t (the top 11 mantissa bits)
//               x low mantissa bits 0, all ones, the cvtpd2ps tie and its neighbours, two random values;
//   random      N inputs (default 10^10): a quarter each raw 64-bit patterns, positive normal patterns, uniform in
//               [0.5, 2) and log-uniform in [1e-30, 1e30]; two inputs per libmvec call (one per lane).
//
//   g++ -std=c++20 -O3 -march=skylake -ffp-contract=off -fopenmp -DCPPWRF_DOUBLE -Isrc tools/libmvec_proof.cpp
//       -o libmvec_proof -lmvec -lm && ./libmvec_proof [N]
// Prints per set and function the inputs compared and those that differ (the first few with their bits), the
// time per call of libmvec and of the C++ functions; exit 1 when any differs.
#include <emmintrin.h>
#include <omp.h>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "core/libmvec_x86.hpp"

extern "C" __m128d _ZGVbN2v_log(__m128d);
extern "C" __m128d _ZGVbN2v_log10(__m128d);

static std::uint64_t B(double v) {
  std::uint64_t b;
  std::memcpy(&b, &v, 8);
  return b;
}
static double D(std::uint64_t b) {
  double v;
  std::memcpy(&v, &b, 8);
  return v;
}
static std::uint64_t splitmix(std::uint64_t& s) {
  std::uint64_t z = (s += 0x9e3779b97f4a7c15ull);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
  return z ^ (z >> 31);
}

struct Count {
  long long n = 0, bad = 0;
};
static int shown = 0;

// compares both functions on the pair (x0, x1) (one libmvec call each)
static inline void pair(double x0, double x1, Count& cl, Count& c10) {
  alignas(16) double v[2], w[2];
  _mm_store_pd(v, _ZGVbN2v_log(_mm_set_pd(x1, x0)));
  _mm_store_pd(w, _ZGVbN2v_log10(_mm_set_pd(x1, x0)));
  const double xs[2] = {x0, x1};
  for (int l = 0; l < 2; ++l) {
    const double a = cppwrf::mvec::log(xs[l]), b = cppwrf::mvec::log10(xs[l]);
    cl.n++;
    c10.n++;
    if (B(a) != B(v[l])) {
      cl.bad++;
#pragma omp critical
      if (shown++ < 20) std::printf("  log   x %016llx: libmvec %016llx, C++ %016llx\n", (unsigned long long)B(xs[l]), (unsigned long long)B(v[l]), (unsigned long long)B(a));
    }
    if (B(b) != B(w[l])) {
      c10.bad++;
#pragma omp critical
      if (shown++ < 20) std::printf("  log10 x %016llx: libmvec %016llx, C++ %016llx\n", (unsigned long long)B(xs[l]), (unsigned long long)B(w[l]), (unsigned long long)B(b));
    }
  }
}

static void report(const char* set, const Count& cl, const Count& c10, double sec) {
  std::printf("%-8s log: %lld inputs, %lld differ; log10: %lld inputs, %lld differ (%.1f s)\n", set, cl.n, cl.bad, c10.n,
              c10.bad, sec);
}

int main(int argc, char** argv) {
  const long long n_random = argc > 1 ? std::atoll(argv[1]) : 10000000000ll;
  char model[64] = "?";
  if (FILE* f = std::fopen("/proc/cpuinfo", "r")) {
    char line[256];
    while (std::fgets(line, sizeof line, f))
      if (!std::strncmp(line, "model name", 10)) {
        std::snprintf(model, sizeof model, "%s", std::strchr(line, ':') + 2);
        model[std::strcspn(model, "\n")] = 0;
        break;
      }
    std::fclose(f);
  }
  std::printf("CPU %s, %d threads\n", model, omp_get_max_threads());
  long long total_bad = 0;
  using clk = std::chrono::steady_clock;

  // edges
  {
    auto t0 = clk::now();
    std::vector<double> e;
    for (int k = -1074; k <= 1023; ++k) {
      const double p = std::ldexp(1.0, k);
      for (double v : {p, std::nextafter(p, 0.0), std::nextafter(p, INFINITY), -p}) e.push_back(v);
    }
    for (std::uint64_t b : {0x0ull, 0x8000000000000000ull, 0x000fffffffffffffull, 0x0010000000000000ull, 0x7fefffffffffffffull,
                            0x7ff0000000000000ull, 0xfff0000000000000ull, 0x7ff8000000000000ull, 0xfff8000000000000ull,
                            0x7ff0000000000001ull, 0x7ff4000000000000ull, 0xfff0000000000001ull, 0x7fffffffffffffffull,
                            0xffffffffffffffffull, 0x7ff8000000000123ull, 0x3fefffffffffffffull, 0x3ff0000000000001ull})
      e.push_back(D(b));
    for (int k = -323; k <= 308; ++k) {
      const double p = std::pow(10.0, k);
      for (double v : {p, std::nextafter(p, 0.0), std::nextafter(p, INFINITY)}) e.push_back(v);
    }
    if (e.size() % 2) e.push_back(1.0);
    Count cl, c10;
    for (std::size_t i = 0; i < e.size(); i += 2) pair(e[i], e[i + 1], cl, c10);
    report("edges", cl, c10, std::chrono::duration<double>(clk::now() - t0).count());
    total_bad += cl.bad + c10.bad;
  }
  // classes
  {
    auto t0 = clk::now();
    long long nl = 0, bl = 0, n10 = 0, b10 = 0;
    const std::uint64_t low_mask = (1ull << 41) - 1, tie = (0xfffull << 29) | (1ull << 28);
#pragma omp parallel for schedule(dynamic) reduction(+ : nl, bl, n10, b10)
    for (int se = 0; se < 4096; ++se) {
      Count cl, c10;
      std::uint64_t s = 0x1234567ull + std::uint64_t(se);
      const std::uint64_t hi = std::uint64_t(se) << 52;  // sign and biased exponent
      for (std::uint64_t t = 0; t < 2048; ++t) {
        const std::uint64_t lows[8] = {0, low_mask, tie, tie - 1, tie + 1, splitmix(s) & low_mask, splitmix(s) & low_mask,
                                       (0xfffull << 29)};
        for (int l = 0; l < 8; l += 2) pair(D(hi | t << 41 | lows[l]), D(hi | t << 41 | lows[l + 1]), cl, c10);
      }
      nl += cl.n;
      bl += cl.bad;
      n10 += c10.n;
      b10 += c10.bad;
    }
    Count cl{nl, bl}, c10{n10, b10};
    report("classes", cl, c10, std::chrono::duration<double>(clk::now() - t0).count());
    total_bad += bl + b10;
  }
  // random
  {
    auto t0 = clk::now();
    long long nl = 0, bl = 0, n10 = 0, b10 = 0;
    const long long chunks = 4096, per = (n_random / chunks + 7) / 8 * 8;
#pragma omp parallel for schedule(dynamic) reduction(+ : nl, bl, n10, b10)
    for (long long c = 0; c < chunks; ++c) {
      Count cl, c10;
      std::uint64_t s = 0xabcdef0123ull * std::uint64_t(c + 1);
      auto draw = [&s](int kind) {
        const std::uint64_t r = splitmix(s);
        switch (kind) {
          case 0: return D(r);
          case 1: return D(0x0010000000000000ull + r % 0x7fe0000000000000ull);
          case 2: return 0.5 + 1.5 * double(r >> 11) * 0x1p-53;
          default: return std::pow(10.0, -30.0 + 60.0 * double(r >> 11) * 0x1p-53);
        }
      };
      for (long long i = 0; i < per; i += 2) {
        const int kind = int((i / 2) % 4);
        pair(draw(kind), draw(kind), cl, c10);
      }
      nl += cl.n;
      bl += cl.bad;
      n10 += c10.n;
      b10 += c10.bad;
    }
    Count cl{nl, bl}, c10{n10, b10};
    report("random", cl, c10, std::chrono::duration<double>(clk::now() - t0).count());
    total_bad += bl + b10;
  }
  // cost per call (one thread): libmvec's 2-lane call per input vs the C++ function
  {
    const int n = 1 << 24;
    std::vector<double> x(n);
    std::uint64_t s = 99;
    for (auto& v : x) v = std::pow(10.0, -5.0 + 10.0 * double(splitmix(s) >> 11) * 0x1p-53);
    double acc = 0;
    auto t0 = clk::now();
    for (int i = 0; i < n; ++i) acc += _mm_cvtsd_f64(_ZGVbN2v_log(_mm_set1_pd(x[i])));
    auto t1 = clk::now();
    for (int i = 0; i < n; ++i) acc += cppwrf::mvec::log(x[i]);
    auto t2 = clk::now();
    for (int i = 0; i < n; ++i) acc += std::log(x[i]);
    auto t3 = clk::now();
    auto ns = [n](auto a, auto b) { return std::chrono::duration<double, std::nano>(b - a).count() / n; };
    std::printf("cost per log: libmvec %.2f ns, C++ (portable) %.2f ns, scalar libm %.2f ns (checksum %g)\n", ns(t0, t1),
                ns(t1, t2), ns(t2, t3), acc);
  }
  std::printf("%s\n", total_bad ? "libmvec_proof: DIFFERENCES" : "libmvec_proof: identical");
  return total_bad ? 1 : 0;
}
