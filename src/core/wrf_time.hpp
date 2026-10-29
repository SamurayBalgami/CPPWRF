#pragma once
// WRF's time intervals (external/esmf_time_f90: ESMF_TimeInterval.F90, ESMF_BaseTime.F90, Meat.F90): seconds as
// S + Sn/Sd in 8-byte integers, normalised after every operation. The adaptive time step (adapt_timestep_em.F)
// computes dt and the model time with this arithmetic; the results (sums, rounded fractions, comparisons) are exact
// and the conversion to a floating-point time (real_time) depends only on the value.
#include <cstdint>
#include <cstdlib>
#include <stdexcept>

namespace cppwrf::wrf_time {

struct Interval {
  std::int64_t S = 0, Sn = 0, Sd = 0;
};

namespace detail {
// Meat.F90 simplify / simplify_numprimes: divides out the first nprimes primes (simplify: 9, up to 23)
inline void simplify(std::int64_t ni, std::int64_t di, std::int64_t& no, std::int64_t& d_o, int nprimes) {
  static constexpr std::int64_t primes[62] = {2,   3,   5,   7,   11,  13,  17,  19,  23,  29,  31,  37,  41,
                                              43,  47,  53,  59,  61,  67,  71,  73,  79,  83,  89,  97,  101,
                                              103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
                                              173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239,
                                              241, 251, 257, 263, 269, 271, 277, 281, 283, 293};
  if (ni == 0) { d_o = 1; no = 0; return; }
  if (di % ni == 0) { d_o = di / ni; no = 1; return; }
  std::int64_t d = di, n = ni;
  for (int np = 0; np < nprimes; ++np) {
    const std::int64_t pr = primes[np];
    bool keepgoing = true;
    while (keepgoing) {
      keepgoing = false;
      if (d / pr != 0 && n / pr != 0 && d % pr == 0 && n % pr == 0) {
        d = d / pr;
        n = n / pr;
        keepgoing = true;
      }
    }
  }
  d_o = d;
  no = n;
}

// Meat.F90 compute_lcd (a common multiple: the product with the primes 2..23 divided out)
inline std::int64_t compute_lcd(std::int64_t e1, std::int64_t e2) {
  static constexpr std::int64_t primes[9] = {2, 3, 5, 7, 11, 13, 17, 19, 23};
  std::int64_t d1 = e1, d2 = e2;
  if (d1 == 0 && d2 == 0) return 1;
  if (d1 == 0) d1 = d2;
  if (d2 == 0) d2 = d1;
  if (d1 == d2) return d1;
  std::int64_t lcd = d1 * d2;
  for (const std::int64_t p : primes)
    while (lcd / p != 0 && (lcd / p) % d1 == 0 && (lcd / p) % d2 == 0) lcd = lcd / p;
  return lcd;
}
}  // namespace detail

// Meat.F90 normalize_basetime
inline void normalize(Interval& t) {
  if (t.Sd < 0) throw std::runtime_error("normalize_basetime:  denominator of seconds cannot be negative");
  if (t.Sd == 0 && t.Sn != 0)
    throw std::runtime_error("normalize_basetime:  denominator of seconds cannot be zero when numerator is non-zero");
  if (t.Sd > 0) {
    if (std::llabs(t.Sn) >= t.Sd) {
      t.S = t.S + t.Sn / t.Sd;
      t.Sn = t.Sn % t.Sd;
    }
    if (t.S > 0 && t.Sn < 0) { t.S = t.S - 1; t.Sn = t.Sn + t.Sd; }
    if (t.S < 0 && t.Sn > 0) { t.S = t.S + 1; t.Sn = t.Sn - t.Sd; }
  }
  if (t.Sd > 0 && t.Sn > 0) {
    std::int64_t n, d;
    detail::simplify(t.Sn, t.Sd, n, d, 9);
    t.Sn = n; t.Sd = d;
    if (t.Sd > 10000 && t.Sn > 10000) {
      detail::simplify(t.Sn, t.Sd, n, d, 62);
      t.Sn = n; t.Sd = d;
    }
  }
}

// ESMF_TimeIntervalSet(S=s, Sn=sn, Sd=sd)
inline Interval make(std::int64_t s, std::int64_t sn = 0, std::int64_t sd = 0) {
  Interval t{s, sd != 0 ? sn : 0, sd};
  normalize(t);
  return t;
}

// ESMF_BaseTimeSum (the interval sum and the time + interval)
inline Interval operator+(const Interval& a, const Interval& b) {
  Interval r = a;
  r.S = r.S + b.S;
  if (a.Sd == 0 && b.Sd == 0) { r.Sn = 0; r.Sd = 0; }
  else if (a.Sd != 0 && b.Sd == 0) { r.Sn = a.Sn; r.Sd = a.Sd; }
  else if (a.Sd == 0 && b.Sd != 0) { r.Sn = b.Sn; r.Sd = b.Sd; }
  else {
    const std::int64_t lcd = detail::compute_lcd(a.Sd, b.Sd);
    r.Sd = lcd;
    r.Sn = (a.Sn * lcd / a.Sd) + (b.Sn * lcd / b.Sd);
  }
  normalize(r);
  return r;
}

// ESMF_BaseTimeDifference: a + (-b)
inline Interval operator-(const Interval& a, const Interval& b) { return a + Interval{-b.S, -b.Sn, b.Sd}; }

// ESMF_TimeIntervalProdI
inline Interval operator*(const Interval& a, int m) {
  Interval r{a.S * std::int64_t(m), a.Sn * std::int64_t(m), a.Sd};
  normalize(r);
  return r;
}

// ESMF_TimeIntervalQuotI (ESMF_BaseTimeQuotI8)
inline Interval operator/(const Interval& a, int divisor) {
  if (divisor == 0) throw std::runtime_error("ESMF_TimeIntervalQuotI:  divide by zero");
  const std::int64_t dinit = a.Sd == 0 ? 1 : a.Sd;
  const std::int64_t n = a.S * dinit + a.Sn, d = dinit * std::int64_t(divisor);
  Interval r;
  detail::simplify(n, d, r.Sn, r.Sd, 9);
  normalize(r);
  return r;
}

// Meat.F90 seccmp: -1, 0, 1
inline int cmp(const Interval& a, const Interval& b) {
  std::int64_t n1 = a.Sn, n2 = b.Sn;
  if (n1 != 0 || n2 != 0) {
    const std::int64_t lcd = detail::compute_lcd(a.Sd, b.Sd);
    if (a.Sd != 0) n1 = n1 * (lcd / a.Sd);
    if (b.Sd != 0) n2 = n2 * (lcd / b.Sd);
  }
  if (a.S > b.S) return 1;
  if (a.S < b.S) return -1;
  return n1 > n2 ? 1 : (n1 < n2 ? -1 : 0);
}
inline bool operator<(const Interval& a, const Interval& b) { return cmp(a, b) < 0; }
inline bool operator>(const Interval& a, const Interval& b) { return cmp(a, b) > 0; }
inline bool operator<=(const Interval& a, const Interval& b) { return cmp(a, b) <= 0; }
inline bool operator>=(const Interval& a, const Interval& b) { return cmp(a, b) >= 0; }
inline bool operator==(const Interval& a, const Interval& b) { return cmp(a, b) == 0; }
inline bool operator!=(const Interval& a, const Interval& b) { return cmp(a, b) != 0; }

// adapt_timestep_em.F real_time: REAL(S) + Sn / REAL(Sd) (R = WRF's REAL)
template <class R>
R real_time(const Interval& t) {
  if (std::llabs(t.Sd) < 1) return R(t.S);
  return R(t.S) + R(t.Sn) / R(t.Sd);
}

// domain_clock_get minutesSinceSimulationStart (WRF's XTIME): days*24*60 + seconds/60 + (Sn/Sd)/60
template <class R>
R minutes(const Interval& t) {
  const std::int64_t days = t.S / 86400, sec = t.S % 86400;
  R m = (R(days) * R(24.) * R(60.)) + (R(sec) / R(60.));
  if (t.Sd != 0) m = m + ((R(t.Sn) / R(t.Sd)) / R(60.));
  return m;
}

}  // namespace cppwrf::wrf_time
