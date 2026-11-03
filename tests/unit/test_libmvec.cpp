// The vector log / log10 of src/core/libmvec_x86.hpp (glibc 2.41 libmvec's SSE4 variants with Intel's rcpps):
//  (a) on any x86-64 CPU their results over a fixed input set (every rcpps input class t at 13 exponents with the
//      float-rounding boundaries, edge values, 2^16 random bit patterns, random normals, random values in [0.5, 2))
//      hash to the values glibc's libmvec gives on Intel (`python3 tools/libmvec_tables.py --hash` on an Intel CPU);
//  (b) on Intel also: equal to libmvec's _ZGVbN2v_log / _ZGVbN2v_log10 bit for bit on that set, and the table equals
//      the CPU's rcpps for every input the algorithms give it (2^23 + 1 floats);
//  (c) elsewhere (AMD) the test reports how many of libmvec's own results differ (information, no failure).
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

#include "core/backend.hpp"

#if defined(CPPWRF_HAVE_LIBMVEC)
#include <cpuid.h>
#include <xmmintrin.h>

static std::uint64_t bits_of(double v) {
  std::uint64_t b;
  std::memcpy(&b, &v, sizeof b);
  return b;
}
static double from_bits(std::uint64_t b) {
  double v;
  std::memcpy(&v, &b, sizeof v);
  return v;
}

// the input set (tools/libmvec_tables.py builds the same)
static std::vector<double> inputs() {
  std::uint64_t s = 0;
  auto sm = [&s] {  // splitmix64
    std::uint64_t z = (s += 0x9e3779b97f4a7c15ull);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
  };
  std::vector<double> x;
  for (std::uint64_t b : {0x0ull, 0x8000000000000000ull, 0x1ull, 0x000fffffffffffffull, 0x8000000000000001ull,
                          0x0010000000000000ull, 0x0010000000000001ull, 0x7fefffffffffffffull, 0x7ff0000000000000ull,
                          0xfff0000000000000ull, 0x7ff8000000000000ull, 0xfff8000000000000ull, 0x7ff0000000000001ull,
                          0x7ff4000000000000ull, 0xbff0000000000000ull, 0x3ff0000000000000ull, 0x3fefffffffffffffull,
                          0x3ff0000000000001ull, 0x4024000000000000ull})
    x.push_back(from_bits(b));
  const std::uint64_t low_mask = (1ull << 41) - 1, tie = (0xfffull << 29) | (1ull << 28);
  for (std::uint64_t e : {1, 2, 512, 1000, 1021, 1022, 1023, 1024, 1025, 1100, 1536, 2045, 2046})
    for (std::uint64_t t = 0; t < 2048; ++t)
      for (std::uint64_t low : {std::uint64_t(0), low_mask, tie, tie - 1, sm() & low_mask})
        x.push_back(from_bits(e << 52 | t << 41 | low));
  for (int i = 0; i < 65536; ++i) x.push_back(from_bits(sm()));
  for (int i = 0; i < 65536; ++i) x.push_back(from_bits(0x0010000000000000ull + sm() % 0x7fe0000000000000ull));
  for (int i = 0; i < 65536; ++i) x.push_back(0.5 + 1.5 * double(sm() >> 11) * 0x1p-53);
  return x;
}

static std::uint64_t hash(const std::vector<double>& y) {
  std::uint64_t h = 0xcbf29ce484222325ull;
  for (double v : y) h = (h ^ bits_of(v)) * 0x100000001b3ull;
  return h;
}

static bool is_intel() {
  unsigned a, b, c, d;
  char v[13] = {};
  if (!__get_cpuid(0, &a, &b, &c, &d)) return false;
  std::memcpy(v, &b, 4);
  std::memcpy(v + 4, &d, 4);
  std::memcpy(v + 8, &c, 4);
  return std::strcmp(v, "GenuineIntel") == 0;
}

// glibc libmvec's results on Intel (python3 tools/libmvec_tables.py --hash)
constexpr std::uint64_t kHashLog = 0x0df23b47a0a5826bull, kHashLog10 = 0x6031508ce600671full;

int main() {
  int fails = 0;
  const bool intel = is_intel();
  const std::vector<double> x = inputs();
  std::vector<double> mlog(x.size()), mlog10(x.size()), vlog(x.size()), vlog10(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    mlog[i] = cppwrf::mvec::log(x[i]);
    mlog10[i] = cppwrf::mvec::log10(x[i]);
    vlog[i] = _mm_cvtsd_f64(_ZGVbN2v_log(_mm_set1_pd(x[i])));
    vlog10[i] = _mm_cvtsd_f64(_ZGVbN2v_log10(_mm_set1_pd(x[i])));
  }
  long dlog = 0, dlog10 = 0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    dlog += bits_of(mlog[i]) != bits_of(vlog[i]);
    dlog10 += bits_of(mlog10[i]) != bits_of(vlog10[i]);
  }
  std::printf("%s CPU, %zu inputs: hash log %016llx (Intel %016llx), log10 %016llx (Intel %016llx); this CPU's "
              "libmvec differs for %ld (log) and %ld (log10)\n",
              intel ? "Intel" : "non-Intel", x.size(), (unsigned long long)hash(mlog), (unsigned long long)kHashLog,
              (unsigned long long)hash(mlog10), (unsigned long long)kHashLog10, dlog, dlog10);
  if (hash(mlog) != kHashLog || hash(mlog10) != kHashLog10) {
    std::printf("FAIL: the portable log / log10 differ from Intel's libmvec\n");
    ++fails;
  }
  if (intel) {
    if (dlog || dlog10) {
      std::printf("FAIL: differs from this Intel CPU's libmvec\n");
      ++fails;
    }
    // the table against this CPU's rcpps: every f = (float)m the algorithms give it (biased exponent 117, and 2^-9)
    long dk = 0;
    for (std::uint32_t fb = 117u << 23; fb <= 118u << 23; ++fb) {
      float f, r;
      std::memcpy(&f, &fb, sizeof f);
      r = _mm_cvtss_f32(_mm_rcp_ps(_mm_set1_ps(f)));
      const double k = __builtin_roundeven(double(r));
      dk += k != cppwrf::mvec::detail::kRcp[(fb >> 12) - (117u << 11)];
    }
    std::printf("rcpps table: %ld of %u inputs differ from this CPU\n", dk, (1u << 23) + 1);
    if (dk) ++fails;
  }
  std::printf("%s\n", fails ? "test_libmvec FAILED" : "test_libmvec passed");
  return fails ? 1 : 0;
}
#else
int main() {
  std::printf("test_libmvec: no libmvec in this build, skipped\n");
  return 0;
}
#endif
