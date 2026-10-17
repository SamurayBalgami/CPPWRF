#pragma once
// Column and level routines of WRF's real.exe (dyn_em/module_initialize_real.F,
// dyn_em/nest_init_utils.F compute_vcoord_1d_coeffs), ported statement by statement in WRF's
// default REAL kind T (float for the single-precision build, double for -r8). Literals are written
// as T(...) so that they have Fortran's default kind; REAL(KIND=8) locals stay double.
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "init_real/farray.hpp"

#if defined(__x86_64__) && !defined(CPPWRF_NO_LIBMVEC)
#include <emmintrin.h>
extern "C" __m128d _ZGVbN2v_log(__m128d);
extern "C" __m128d _ZGVbN2v_exp(__m128d);
extern "C" __m128 _ZGVbN4v_logf(__m128);
extern "C" __m128 _ZGVbN4v_expf(__m128);
#define CPPWRF_INIT_LIBMVEC 1
#endif

namespace cppwrf::init_real {

// module_model_constants.F in the default REAL kind
template <typename T>
struct Const {
  static constexpr T g = T(9.81);
  static constexpr T r_d = T(287.);
  static constexpr T cp = T(7.) * r_d / T(2.);
  static constexpr T r_v = T(461.6);
  static constexpr T cv = cp - r_d;
  static constexpr T cvpm = -cv / cp;
  static constexpr T rcp = r_d / cp;
  static constexpr T p1000mb = T(100000.);
  static constexpr T t0 = T(300.);
  static constexpr T cpovcv = cp / (cp - r_d);
  static constexpr T rvovrd = r_v / r_d;
};

// gfortran -O2 -ftree-vectorize evaluates elementwise LOG/EXP of an array section with glibc's
// libmvec: groups of VL elements from the start (2 for REAL*8, 4 for REAL*4), the remainder with
// the scalar function. vlog(x, i, n): element i (0-based) of an n-element array expression.
template <typename T>
inline bool vec_lane(int i, int n) {
  constexpr int vl = sizeof(T) == 8 ? 2 : 4;
  return i < n - n % vl;
}
inline double vlog(double x, int i, int n) {
#if defined(CPPWRF_INIT_LIBMVEC)
  if (vec_lane<double>(i, n)) return _mm_cvtsd_f64(_ZGVbN2v_log(_mm_set1_pd(x)));
#endif
  (void)i; (void)n;
  return std::log(x);
}
inline double vexp(double x, int i, int n) {
#if defined(CPPWRF_INIT_LIBMVEC)
  if (vec_lane<double>(i, n)) return _mm_cvtsd_f64(_ZGVbN2v_exp(_mm_set1_pd(x)));
#endif
  (void)i; (void)n;
  return std::exp(x);
}
inline float vlog(float x, int i, int n) {
#if defined(CPPWRF_INIT_LIBMVEC)
  if (vec_lane<float>(i, n)) return _mm_cvtss_f32(_ZGVbN4v_logf(_mm_set1_ps(x)));
#endif
  (void)i; (void)n;
  return std::log(x);
}
inline float vexp(float x, int i, int n) {
#if defined(CPPWRF_INIT_LIBMVEC)
  if (vec_lane<float>(i, n)) return _mm_cvtss_f32(_ZGVbN4v_expf(_mm_set1_ps(x)));
#endif
  (void)i; (void)n;
  return std::exp(x);
}

// Fortran x**n with a constant integer n (gfortran expands it into multiplications)
template <typename T> inline T sq(T x) { return x * x; }
template <typename T> inline T cube(T x) { return x * x * x; }
template <typename T> inline T pow4(T x) { T y = x * x; return y * y; }

struct Dims {
  int ids, ide, jds, jde, kds, kde;  // domain (staggered ends), serial: memory and tile = domain
  int its() const { return ids; }
  int ite() const { return ide; }
  int jts() const { return jds; }
  int jte() const { return jde; }
  int kts() const { return kds; }
  int kte() const { return kde; }
};

// share/module_soil_pre.F skip_middle_points_t
inline bool skip_middle(const Dims& d, int i, int j, int width, bool hold_ups) {
  if (hold_ups || d.ide == 3 || d.jde == 3) return false;
  return i >= d.ids + width && i <= d.ide - 1 - width && j >= d.jds + width && j <= d.jde - 1 - width;
}

struct SkipCtx {
  int width = 5;
  bool hold_ups = true;
};

// ---------------------------------------------------------------- lagrange_interp / setup
template <typename T>
void lagrange_interp(const T* x, const T* y, int n, T target_x, T& target_y) {  // x, y: 0..n
  T Px = T(0.);
  for (int i = 0; i <= n; ++i) {
    T numer = T(1.), denom = T(1.);
    for (int k = 0; k <= n; ++k) {
      if (k == i) continue;
      numer = numer * (target_x - x[k]);
      denom = denom * (x[i] - x[k]);
    }
    if (denom != T(0.)) {
      T Ln = y[i] * numer / denom;
      Px = Px + Ln;
    }
  }
  target_y = Px;
}

// all_x, all_y: 1..all_dim (passed 0-based), target_x/y: 1..target_dim
template <typename T>
void lagrange_setup(char var_type, int interp_type, const T* all_x, const T* all_y, int all_dim, int n,
                    int extrap_type, const T* target_x, T* target_y, int target_dim, int ii, int jj) {
  if (all_dim < n + 1)
    throw std::runtime_error("real: vertical interpolation order too large at (" + std::to_string(ii) + "," +
                             std::to_string(jj) + ")");
  if (n < 1) throw std::runtime_error("real: lagrange order < 1");
  constexpr T RovCp = Const<T>::rcp;
  const T CRC_const1 = T(11880.516), CRC_const2 = T(0.1902632), CRC_const3 = T(0.0065);
  const int vboundb = 4, vboundt = 0;
  std::vector<T> all_x_full(static_cast<std::size_t>(all_dim)), target_x_full(static_cast<std::size_t>(target_dim));
  for (int tl = 1; tl <= target_dim; ++tl) {
    bool found_loc = false;
    int lcl = 0, lcr = 0;
    for (int loop = 1; loop <= all_dim - 1; ++loop) {
      T a = target_x[tl - 1] - all_x[loop - 1];
      T b = target_x[tl - 1] - all_x[loop];
      if (a * b <= T(0.0)) {
        lcl = loop;
        lcr = loop + 1;
        found_loc = true;
        break;
      }
    }
    if (!found_loc && target_x[tl - 1] > all_x[0]) {
      if (interp_type == 1) {
        for (int m = 0; m < all_dim; ++m) all_x_full[std::size_t(m)] = all_x[m];
        for (int m = 0; m < target_dim; ++m) target_x_full[std::size_t(m)] = target_x[m];
      } else {
        for (int m = 0; m < all_dim; ++m) all_x_full[std::size_t(m)] = vexp(all_x[m], m, all_dim);
        for (int m = 0; m < target_dim; ++m) target_x_full[std::size_t(m)] = vexp(target_x[m], m, target_dim);
      }
      const T txf = target_x_full[std::size_t(tl - 1)];
      if (extrap_type == 1 && var_type == 'T') {
        T temp_1 = all_y[0] * std::pow(all_x_full[0] / T(100000.), RovCp);
        target_y[tl - 1] = temp_1 * std::pow(T(100000.) / txf, RovCp);
      } else if (extrap_type == 2 && var_type == 'T') {
        T depth_of_extrap_in_p = txf - all_x_full[0];
        T avg_of_extrap_p = (txf + all_x_full[0]) * T(0.5);
        T temp_extrap_starting_point = all_y[0] * std::pow(all_x_full[0] / T(100000.), RovCp);
        T dhdp = CRC_const1 * CRC_const2 * std::pow(avg_of_extrap_p / T(100.), CRC_const2 - T(1.));
        T dh = dhdp * (depth_of_extrap_in_p / T(100.));
        T dt = dh * CRC_const3;
        target_y[tl - 1] = (temp_extrap_starting_point + dt) * std::pow(T(100000.) / txf, RovCp);
      } else if (extrap_type == 3 && var_type == 'T') {
        target_y[tl - 1] = all_y[0];
      } else if (extrap_type == 1) {
        target_y[tl - 1] = (all_y[1] * (target_x[tl - 1] - all_x[2]) + all_y[2] * (all_x[1] - target_x[tl - 1])) /
                           (all_x[1] - all_x[2]);
      } else if (extrap_type == 2) {
        target_y[tl - 1] = all_y[0];
      } else if (extrap_type == 3) {
        throw std::runtime_error("real: extrap_option 3 only for theta");
      }
      continue;
    } else if (!found_loc) {
      throw std::runtime_error("real: could not find trapping x locations at (" + std::to_string(ii) + "," +
                               std::to_string(jj) + ")");
    }
    if (n == 9) {
      throw std::runtime_error("real: lagrange_order 9 (cubic spline) is not ported");
    } else if (n % 2 != 0) {
      if (lcl - (((n + 1) / 2) - 1) >= 1 && lcr + (((n + 1) / 2) - 1) <= all_dim) {
        int ist = lcl - (((n + 1) / 2) - 1);
        lagrange_interp(all_x + (ist - 1), all_y + (ist - 1), n, target_x[tl - 1], target_y[tl - 1]);
      } else {
        // WRF leaves target_y unset here (only 2nd order is handled)
      }
    } else if (n % 2 == 0 && tl >= 1 + vboundb && tl <= target_dim - vboundt) {
      if (lcl - ((n / 2) - 1) >= 1 && lcr + (n / 2) <= all_dim && lcl - (n / 2) >= 1 && lcr + ((n / 2) - 1) <= all_dim) {
        T y1, y2;
        int ist = lcl - ((n / 2) - 1);
        lagrange_interp(all_x + (ist - 1), all_y + (ist - 1), n, target_x[tl - 1], y1);
        ist = lcl - (n / 2);
        lagrange_interp(all_x + (ist - 1), all_y + (ist - 1), n, target_x[tl - 1], y2);
        target_y[tl - 1] = (y1 + y2) * T(0.5);
      } else if (lcl - ((n / 2) - 1) >= 1 && lcr + (n / 2) <= all_dim) {
        int ist = lcl - ((n / 2) - 1);
        lagrange_interp(all_x + (ist - 1), all_y + (ist - 1), n, target_x[tl - 1], target_y[tl - 1]);
      } else if (lcl - (n / 2) >= 1 && lcr + ((n / 2) - 1) <= all_dim) {
        int ist = lcl - (n / 2);
        lagrange_interp(all_x + (ist - 1), all_y + (ist - 1), n, target_x[tl - 1], target_y[tl - 1]);
      } else {
        throw std::runtime_error("real: unauthorized area in lagrange_setup");
      }
    } else if (n % 2 == 0) {
      lagrange_interp(all_x + (lcl - 1), all_y + (lcl - 1), 1, target_x[tl - 1], target_y[tl - 1]);
    }
  }
}

// ---------------------------------------------------------------- vert_interp
struct VInterpOpts {
  int interp_type = 2, lagrange_order = 2, extrap_type = 2;
  bool lowest_lev_from_sfc = false, use_levels_below_ground = true, use_surface = true;
  double zap_close_levels = 500;
  int force_sfc_in_vinterp = 1;
};

// fo, po: (i, 1..generic, j); fnew, pnu: (i, kds..kde, j). flag_maxw / flag_trop are 0 in our use
// (use_maxw_level = use_trop_level = 0).
template <typename T>
void vert_interp(const F3<T>& fo, const F3<T>& po, F3<T>& fnew, const F3<T>& pnu, int generic, char var_type,
                 const VInterpOpts& o, const Dims& d, const SkipCtx& sk) {
  const int ids = d.ids, ide = d.ide, jds = d.jds, jde = d.jde, its = d.its(), ite = d.ite(), jts = d.jts(),
            jte = d.jte(), kts = d.kts(), kte = d.kte();
  const T zap_close_levels = T(o.zap_close_levels);
  F3<T> porig(ids, ide, 1, generic, jds, jde), forig(ids, ide, 1, generic, jds, jde), pnew(ids, ide, kts, kte, jds, jde);
  int istart, iend, jstart, jend, kstart, kend;
  auto skip = [&](int i, int j) { return skip_middle(d, i, j, sk.width, sk.hold_ups); };
  if (var_type == 'U') {
    istart = its; iend = ite; jstart = std::max(jds, jts - 1); jend = std::min(jde - 1, jte + 1);
    kstart = kts; kend = kte - 1;
    for (int j = jstart; j <= jend; ++j) {
      for (int k = 1; k <= generic; ++k)
        for (int i = std::max(ids + 1, its - 1); i <= std::min(ide - 1, ite + 1); ++i) {
          if (skip(i, j)) continue;
          porig(i, k, j) = (po(i, k, j) + po(i - 1, k, j)) * T(0.5);
        }
      if (ids == its)
        for (int k = 1; k <= generic; ++k) porig(its, k, j) = po(its, k, j);
      if (ide == ite)
        for (int k = 1; k <= generic; ++k) porig(ite, k, j) = po(ite - 1, k, j);
      for (int k = kstart; k <= kend; ++k)
        for (int i = std::max(ids + 1, its - 1); i <= std::min(ide - 1, ite + 1); ++i) {
          if (skip(i, j)) continue;
          pnew(i, k, j) = (pnu(i, k, j) + pnu(i - 1, k, j)) * T(0.5);
        }
      if (ids == its)
        for (int k = kstart; k <= kend; ++k) pnew(its, k, j) = pnu(its, k, j);
      if (ide == ite)
        for (int k = kstart; k <= kend; ++k) pnew(ite, k, j) = pnu(ite - 1, k, j);
    }
  } else if (var_type == 'V') {
    istart = std::max(ids, its - 1); iend = std::min(ide - 1, ite + 1); jstart = jts; jend = jte;
    kstart = kts; kend = kte - 1;
    for (int i = istart; i <= iend; ++i) {
      for (int k = 1; k <= generic; ++k)
        for (int j = std::max(jds + 1, jts - 1); j <= std::min(jde - 1, jte + 1); ++j) {
          if (skip(i, j)) continue;
          porig(i, k, j) = (po(i, k, j) + po(i, k, j - 1)) * T(0.5);
        }
      if (jds == jts)
        for (int k = 1; k <= generic; ++k) porig(i, k, jts) = po(i, k, jts);
      if (jde == jte)
        for (int k = 1; k <= generic; ++k) porig(i, k, jte) = po(i, k, jte - 1);
      for (int k = kstart; k <= kend; ++k)
        for (int j = std::max(jds + 1, jts - 1); j <= std::min(jde - 1, jte + 1); ++j) {
          if (skip(i, j)) continue;
          pnew(i, k, j) = (pnu(i, k, j) + pnu(i, k, j - 1)) * T(0.5);
        }
      if (jds == jts)
        for (int k = kstart; k <= kend; ++k) pnew(i, k, jts) = pnu(i, k, jts);
      if (jde == jte)
        for (int k = kstart; k <= kend; ++k) pnew(i, k, jte) = pnu(i, k, jte - 1);
    }
  } else {
    istart = its; iend = std::min(ide - 1, ite); jstart = jts; jend = std::min(jde - 1, jte);
    kstart = kts; kend = (var_type == 'W' || var_type == 'Z') ? kte : kte - 1;
    for (int j = jstart; j <= jend; ++j) {
      for (int k = 1; k <= generic; ++k)
        for (int i = istart; i <= iend; ++i) {
          if (skip(i, j)) continue;
          porig(i, k, j) = po(i, k, j);
        }
      for (int k = kstart; k <= kend; ++k)
        for (int i = istart; i <= iend; ++i) {
          if (skip(i, j)) continue;
          pnew(i, k, j) = pnu(i, k, j);
        }
    }
  }
  bool any_valid = false;
  int i_valid = 0, j_valid = 0;
  for (int j = jstart; j <= jend && !any_valid; ++j)
    for (int i = istart; i <= iend; ++i) {
      if (skip(i, j)) continue;
      any_valid = true; i_valid = i; j_valid = j;
      break;
    }
  if (!any_valid) return;
  const bool flip = porig(i_valid, 2, j_valid) < porig(i_valid, generic, j_valid);
  std::vector<int> ko_above_sfc(static_cast<std::size_t>(ide - ids + 1));
  auto kas = [&](int i) -> int& { return ko_above_sfc[std::size_t(i - ids)]; };
  const int nord = generic + 2;  // flag_maxw + flag_trop = 0; room for one insertion
  std::vector<T> op(static_cast<std::size_t>(nord + 1)), of(static_cast<std::size_t>(nord + 1));  // 1-based
  std::vector<T> opn(static_cast<std::size_t>(kte + 1)), ofn(static_cast<std::size_t>(kte + 1));  // ordered_pnew/fnew (kts..kte)
  std::vector<T> lx(static_cast<std::size_t>(nord)), lt(static_cast<std::size_t>(kte + 1));
  for (int j = jstart; j <= jend; ++j) {
    if (flip) {
      for (int kn = 2; kn <= (generic + 1) / 2; ++kn)
        for (int i = istart; i <= iend; ++i) {
          if (skip(i, j)) continue;
          T hold = porig(i, kn, j);
          porig(i, kn, j) = porig(i, generic + 2 - kn, j);
          porig(i, generic + 2 - kn, j) = hold;
          forig(i, kn, j) = fo(i, generic + 2 - kn, j);
          forig(i, generic + 2 - kn, j) = fo(i, kn, j);
        }
      for (int i = istart; i <= iend; ++i) forig(i, 1, j) = fo(i, 1, j);
      if (generic % 2 == 0) {
        int k = generic / 2 + 1;
        for (int i = istart; i <= iend; ++i) forig(i, k, j) = fo(i, k, j);
      }
    } else {
      for (int kn = 1; kn <= generic; ++kn)
        for (int i = istart; i <= iend; ++i) forig(i, kn, j) = fo(i, kn, j);
    }
    for (int i = istart; i <= iend; ++i) kas(i) = -1;
    for (int ko = kstart + 1; ko <= generic; ++ko)
      for (int i = istart; i <= iend; ++i) {
        if (skip(i, j)) continue;
        if (kas(i) == -1 && porig(i, 1, j) > porig(i, ko, j)) kas(i) = ko;
      }
    for (int i = istart; i <= iend; ++i) {
      if (skip(i, j)) continue;
      int zap_below = 0, count, knext, kst;
      if (kas(i) > 2) {
        count = 1;
        for (int ko = 2; ko <= kas(i) - 1; ++ko) {
          op[std::size_t(count)] = porig(i, ko, j);
          of[std::size_t(count)] = forig(i, ko, j);
          count = count + 1;
        }
        if (op[std::size_t(count - 1)] - porig(i, 1, j) < zap_close_levels) {
          count = count - 1;
          zap_below = 1;
        }
        op[std::size_t(count)] = porig(i, 1, j);
        of[std::size_t(count)] = forig(i, 1, j);
        count = count + 1;
        if (o.force_sfc_in_vinterp > 0) {
          knext = kas(i);
          for (int ko = kas(i); ko <= generic; ++ko)
            if (porig(i, ko, j) <= pnew(i, o.force_sfc_in_vinterp, j)) {
              knext = ko;
              break;
            }
        } else {
          knext = kas(i);
        }
        if (op[std::size_t(count - 1)] - porig(i, knext, j) < zap_close_levels) kst = knext + 1;
        else kst = knext;
        for (int ko = kst; ko <= generic; ++ko) {
          op[std::size_t(count)] = porig(i, ko, j);
          of[std::size_t(count)] = forig(i, ko, j);
          count = count + 1;
        }
      } else {
        op[1] = porig(i, 1, j);
        of[1] = forig(i, 1, j);
        count = 2;
        if (o.force_sfc_in_vinterp > 0) {
          knext = 2;
          for (int ko = 2; ko <= generic; ++ko)
            if (porig(i, ko, j) <= pnew(i, o.force_sfc_in_vinterp, j)) {
              knext = ko;
              break;
            }
        } else {
          knext = 2;
        }
        for (int ko = knext; ko <= generic; ++ko) {
          if (op[std::size_t(count - 1)] - porig(i, ko, j) < zap_close_levels && ko < generic) continue;
          op[std::size_t(count)] = porig(i, ko, j);
          of[std::size_t(count)] = forig(i, ko, j);
          count = count + 1;
        }
      }
      for (int kn = kstart; kn <= kend; ++kn) opn[std::size_t(kn)] = pnew(i, kn, j);
      int kinterp_start = 1, kinterp_end = 1;
      if (o.use_levels_below_ground && o.use_surface) {
        count = 0;
        for (int ko = 1; ko <= generic; ++ko) {
          count = count + 1;
          if (porig(i, generic, j) == op[std::size_t(ko)]) break;
        }
        kinterp_start = 1;
        kinterp_end = kinterp_start + count - 1;
      } else if (o.use_levels_below_ground && !o.use_surface) {
        int sfc_level = 0;
        for (int ko = 1; ko <= generic; ++ko)
          if (porig(i, 1, j) == op[std::size_t(ko)]) {
            sfc_level = ko;
            break;
          }
        for (int ko = sfc_level; ko <= generic - 1; ++ko) {
          op[std::size_t(ko)] = op[std::size_t(ko + 1)];
          of[std::size_t(ko)] = of[std::size_t(ko + 1)];
        }
        op[std::size_t(generic)] = T(1.E-5);
        of[std::size_t(generic)] = T(1.E10);
        count = 0;
        for (int ko = 1; ko <= generic; ++ko) {
          count = count + 1;
          if (porig(i, generic, j) == op[std::size_t(ko)]) break;
        }
        kinterp_start = 1;
        kinterp_end = kinterp_start + count - 1;
      } else if (!o.use_levels_below_ground && o.use_surface) {
        int kcount = kas(i) - 1 - zap_below;
        count = 0;
        for (int ko = 1; ko <= generic; ++ko)
          if (porig(i, ko, j) == op[std::size_t(kcount)]) {
            kcount = kcount + 1;
            count = count + 1;
          }
        kinterp_start = kas(i) - 1 - zap_below;
        kinterp_end = kinterp_start + count - 1;
      }
      const int na = kinterp_end - kinterp_start + 1, nt = kend - kstart + 1;
      if (o.interp_type == 1) {
        for (int m = 0; m < na; ++m) lx[std::size_t(m)] = op[std::size_t(kinterp_start + m)];
        for (int m = 0; m < nt; ++m) lt[std::size_t(m)] = opn[std::size_t(kstart + m)];
      } else {
        for (int m = 0; m < na; ++m) lx[std::size_t(m)] = vlog(op[std::size_t(kinterp_start + m)], m, na);
        for (int m = 0; m < nt; ++m) lt[std::size_t(m)] = vlog(opn[std::size_t(kstart + m)], m, nt);
      }
      lagrange_setup<T>(var_type, o.interp_type, lx.data(), of.data() + kinterp_start, na, o.lagrange_order,
                        o.extrap_type, lt.data(), ofn.data() + kstart, nt, i, j);
      for (int kn = kstart; kn <= kend; ++kn) fnew(i, kn, j) = ofn[std::size_t(kn)];
      if (o.lowest_lev_from_sfc) fnew(i, 1, j) = forig(i, 1, j);
    }
  }
}

// ---------------------------------------------------------------- simple column routines
template <typename T>
void p_dry(const F2<T>& mu0, const F1<T>& eta, T pdht, F3<T>& pdry, bool full_levs, const F1<T>& c3f,
           const F1<T>& c3h, const F1<T>& c4f, const F1<T>& c4h, const Dims& d, const SkipCtx& sk) {
  (void)eta;
  const int kte = d.kte();
  for (int j = d.jts(); j <= std::min(d.jde - 1, d.jte()); ++j)
    for (int k = d.kts(); k <= (full_levs ? kte : kte - 1); ++k)
      for (int i = d.its(); i <= std::min(d.ide - 1, d.ite()); ++i) {
        if (skip_middle(d, i, j, sk.width, sk.hold_ups)) continue;
        pdry(i, k, j) = full_levs ? c3f(k) * mu0(i, j) + c4f(k) + pdht : c3h(k) * mu0(i, j) + c4h(k) + pdht;
      }
}

template <typename T>
void p_dts(F2<T>& pdts, const F2<T>& intq, const F2<T>& psfc, T p_top, const Dims& d, const SkipCtx& sk) {
  for (int j = d.jts(); j <= std::min(d.jde - 1, d.jte()); ++j)
    for (int i = d.its(); i <= std::min(d.ide - 1, d.ite()); ++i) {
      if (skip_middle(d, i, j, sk.width, sk.hold_ups)) continue;
      pdts(i, j) = psfc(i, j) - intq(i, j) - p_top;
    }
}

template <typename T>
void p_dhs(F2<T>& pdhs, const F2<T>& ht, T p0, T t0, T a, const Dims& d, const SkipCtx& sk) {
  const T Rd = Const<T>::r_d, g = Const<T>::g;
  for (int j = d.jts(); j <= std::min(d.jde - 1, d.jte()); ++j)
    for (int i = d.its(); i <= std::min(d.ide - 1, d.ite()); ++i) {
      if (skip_middle(d, i, j, sk.width, sk.hold_ups)) continue;
      pdhs(i, j) = p0 * std::exp(-t0 / a + std::sqrt(sq(t0 / a) - T(2.) * g * ht(i, j) / (a * Rd)));
    }
}

// p: (i, 1..nlev, j)
template <typename T>
T find_p_top(const F3<T>& p, int nlev, const Dims& d, const SkipCtx& sk) {
  int i = d.its(), j = d.jts();
  T p_top = p(i, 2, j);
  int min_lev = 2;
  for (int k = 2; k <= nlev; ++k)
    if (p_top > p(i, k, j)) {
      p_top = p(i, k, j);
      min_lev = k;
    }
  int k = min_lev;
  p_top = p(d.its(), k, d.jts());
  for (j = d.jts(); j <= std::min(d.jde - 1, d.jte()); ++j)
    for (i = d.its(); i <= std::min(d.ide - 1, d.ite()); ++i) {
      if (skip_middle(d, i, j, sk.width, sk.hold_ups)) continue;
      p_top = std::max(p_top, p(i, k, j));
    }
  return p_top;
}

// t_to_theta / theta_to_t over k = kts..kte (theta_to_t: kte-1), Fortran's (p00/p)**(Rd/Cp)
template <typename T>
void t_to_theta(F3<T>& t, const F3<T>& p, T p00, int kte, const Dims& d, const SkipCtx& sk) {
  const T e = Const<T>::r_d / Const<T>::cp;
  for (int j = d.jts(); j <= std::min(d.jde - 1, d.jte()); ++j)
    for (int k = d.kts(); k <= kte; ++k)
      for (int i = d.its(); i <= std::min(d.ide - 1, d.ite()); ++i) {
        if (skip_middle(d, i, j, sk.width, sk.hold_ups)) continue;
        t(i, k, j) = t(i, k, j) * std::pow(p00 / p(i, k, j), e);
      }
}

template <typename T>
void theta_to_t(F3<T>& t, const F3<T>& p, T p00, const Dims& d, const SkipCtx& sk) {
  const T e = Const<T>::r_d / Const<T>::cp;
  for (int j = d.jts(); j <= std::min(d.jde - 1, d.jte()); ++j)
    for (int k = d.kts(); k <= d.kte() - 1; ++k)
      for (int i = d.its(); i <= std::min(d.ide - 1, d.ite()); ++i) {
        if (skip_middle(d, i, j, sk.width, sk.hold_ups)) continue;
        if (p(i, k, j) != T(0.)) t(i, k, j) = t(i, k, j) / std::pow(p00 / p(i, k, j), e);
      }
}

// rh_to_mxrat1 over k = kts..kte_
struct QvLimits {
  bool wrt_liquid = true;
  double qv_max_p_safe = 10000, qv_max_flag = 1.E-5, qv_max_value = 3.E-6;
  double qv_min_p_safe = 110000, qv_min_flag = 1.E-6, qv_min_value = 1.E-6;
};

template <typename T>
void rh_to_mxrat1(F3<T>& rh, const F3<T>& t, const F3<T>& p, F3<T>& q, const QvLimits& L, int kte_, const Dims& d,
                  const SkipCtx& sk) {
  const T T_REF = T(0.0), MW_AIR = T(28.966), MW_VAP = T(18.0152);
  const T A0 = T(6.107799961), A1 = T(4.436518521e-01), A2 = T(1.428945805e-02), A3 = T(2.650648471e-04),
          A4 = T(3.031240396e-06), A5 = T(2.034080948e-08), A6 = T(6.136820929e-11);
  const T ES0 = T(6.1121), C1 = T(9.09718), C2 = T(3.56654), C3 = T(0.876793), EIS = T(6.1071), TF = T(273.16);
  const T EPS = T(0.622), SVP1 = T(0.6112), SVP2 = T(17.67), SVP3 = T(29.65), SVPT0 = T(273.15);
  const int i1 = std::min(d.ide - 1, d.ite()), j1 = std::min(d.jde - 1, d.jte());
  auto skip = [&](int i, int j) { return skip_middle(d, i, j, sk.width, sk.hold_ups); };
  for (int j = d.jts(); j <= j1; ++j)
    for (int k = d.kts(); k <= kte_; ++k)
      for (int i = d.its(); i <= i1; ++i) {
        if (skip(i, j)) continue;
        rh(i, k, j) = std::min(std::max(rh(i, k, j), T(0.)), T(100.));
      }
  if (L.wrt_liquid) {
    for (int j = d.jts(); j <= j1; ++j)
      for (int k = d.kts(); k <= kte_; ++k)
        for (int i = d.its(); i <= i1; ++i) {
          if (skip(i, j)) continue;
          if (t(i, k, j) != T(0.)) {
            T es = T(.01) * rh(i, k, j) * SVP1 * T(10.) * std::exp(SVP2 * (t(i, k, j) - SVPT0) / (t(i, k, j) - SVP3));
            if (es >= p(i, k, j) / T(100.)) q(i, k, j) = T(1.E-6);
            else q(i, k, j) = std::max(EPS * es / (p(i, k, j) / T(100.) - es), T(1.E-6));
          } else {
            q(i, k, j) = T(1.E-6);
          }
        }
  } else {
    for (int j = d.jts(); j <= j1; ++j)
      for (int k = d.kts(); k <= kte_; ++k)
        for (int i = d.its(); i <= i1; ++i) {
          if (skip(i, j)) continue;
          T t1 = t(i, k, j) - T(273.16);
          if (t1 < T(-200.)) {
            q(i, k, j) = T(0);
          } else {
            T ew = T(0);
            if (t1 >= T_REF) {
              ew = A0 + t1 * (A1 + t1 * (A2 + t1 * (A3 + t1 * (A4 + t1 * (A5 + t1 * A6)))));
            } else if (t1 < T_REF && t1 >= T(-47.)) {
              ew = ES0 * std::exp(T(17.67) * t1 / (t1 + T(243.5)));
            } else if (t1 < T(-47.)) {
              T tk = t(i, k, j);
              T rhs = -C1 * (TF / tk - T(1.)) - C2 * std::log10(TF / tk) + C3 * (T(1.) - tk / TF) + std::log10(EIS);
              ew = std::pow(T(10.), rhs);
            }
            ew = std::max(ew, T(0.)) * rh(i, k, j) * T(0.01);
            T q1 = MW_VAP * ew;
            q1 = q1 / (q1 + MW_AIR * (p(i, k, j) / T(100.) - ew));
            q(i, k, j) = q1 / (T(1.) - q1);
          }
        }
  }
  const T maxps = T(L.qv_max_p_safe), maxf = T(L.qv_max_flag), maxv = T(L.qv_max_value);
  const T minps = T(L.qv_min_p_safe), minf = T(L.qv_min_flag), minv = T(L.qv_min_value);
  for (int j = d.jts(); j <= j1; ++j)
    for (int k = d.kts(); k <= kte_; ++k)
      for (int i = d.its(); i <= i1; ++i) {
        if (skip(i, j)) continue;
        if (p(i, k, j) < maxps && q(i, k, j) > maxf) q(i, k, j) = maxv;
        if (p(i, k, j) < minps && q(i, k, j) < minf) q(i, k, j) = minv;
      }
}

// integ_moist over metgrid levels 1..nlev
template <typename T>
void integ_moist(const F3<T>& q_in, const F3<T>& p_in, F3<T>& pd_out, const F3<T>& t_in, const F3<T>& ght_in,
                 F2<T>& intq, int nlev, const Dims& d, const SkipCtx& sk) {
  const int kts = 1, kte = nlev;
  const T Rd = Const<T>::r_d, g = Const<T>::g;
  const int i1 = std::min(d.ide - 1, d.ite()), j1 = std::min(d.jde - 1, d.jte());
  auto skip = [&](int i, int j) { return skip_middle(d, i, j, sk.width, sk.hold_ups); };
  bool upside_down = false, assigned = false;
  for (int j = d.jts(); j <= j1 && !assigned; ++j)
    for (int i = d.its(); i <= i1; ++i) {
      if (skip(i, j)) continue;
      upside_down = p_in(i, kts + 1, j) < p_in(i, kte, j);
      assigned = true;
      break;
    }
  F2<T> psfc(d.ids, d.ide, d.jds, d.jde), tsfc(d.ids, d.ide, d.jds, d.jde), qsfc(d.ids, d.ide, d.jds, d.jde),
      zsfc(d.ids, d.ide, d.jds, d.jde);
  for (int j = d.jts(); j <= j1; ++j)
    for (int i = d.its(); i <= i1; ++i) {
      if (skip(i, j)) continue;
      psfc(i, j) = p_in(i, kts, j); tsfc(i, j) = t_in(i, kts, j); qsfc(i, j) = q_in(i, kts, j);
      zsfc(i, j) = ght_in(i, kts, j);
    }
  std::vector<T> q(static_cast<std::size_t>(kte + 1)), p(static_cast<std::size_t>(kte + 1)), t(static_cast<std::size_t>(kte + 1)), ght(static_cast<std::size_t>(kte + 1)),
      pd(std::size_t(kte + 1));
  for (int j = d.jts(); j <= j1; ++j) {
    for (int i = d.its(); i <= i1; ++i) intq(i, j) = T(0.);
    for (int i = d.its(); i <= i1; ++i) {
      if (skip(i, j)) continue;
      if (upside_down) {
        p[1] = p_in(i, kts, j); t[1] = t_in(i, kts, j); q[1] = q_in(i, kts, j); ght[1] = ght_in(i, kts, j);
        for (int k = kts + 1; k <= kte; ++k) {
          p[std::size_t(k)] = p_in(i, kte + 2 - k, j); t[std::size_t(k)] = t_in(i, kte + 2 - k, j);
          q[std::size_t(k)] = q_in(i, kte + 2 - k, j); ght[std::size_t(k)] = ght_in(i, kte + 2 - k, j);
        }
      } else {
        for (int k = kts; k <= kte; ++k) {
          p[std::size_t(k)] = p_in(i, k, j); t[std::size_t(k)] = t_in(i, k, j);
          q[std::size_t(k)] = q_in(i, k, j); ght[std::size_t(k)] = ght_in(i, k, j);
        }
      }
      int las = -1;
      if (p[kts + 1] < psfc(i, j)) {
        las = kts + 1;
      } else {
        for (int k = kts + 1; k <= kte - 1; ++k)
          if (p[std::size_t(k)] - psfc(i, j) >= T(0.) && p[std::size_t(k + 1)] - psfc(i, j) < T(0.)) {
            las = k + 1;
            break;
          }
        if (las == -1) throw std::runtime_error("real: integ_moist could not find level above ground");
      }
      pd[std::size_t(kte)] = p[std::size_t(kte)];
      for (int k = kte - 1; k >= las; --k) {
        T rhobar = (p[std::size_t(k)] / (Rd * t[std::size_t(k)]) + p[std::size_t(k + 1)] / (Rd * t[std::size_t(k + 1)])) * T(0.5);
        T qbar = (q[std::size_t(k)] + q[std::size_t(k + 1)]) * T(0.5);
        T dz = ght[std::size_t(k + 1)] - ght[std::size_t(k)];
        intq(i, j) = intq(i, j) + g * qbar * rhobar / (T(1.) + qbar) * dz;
        pd[std::size_t(k)] = p[std::size_t(k)] - intq(i, j);
      }
      if (p[std::size_t(las - 1)] - psfc(i, j) >= T(0.) && p[std::size_t(las)] - psfc(i, j) < T(0.) && las > kts) {
        T p1 = psfc(i, j), p2 = p[std::size_t(las)], t1 = tsfc(i, j), t2 = t[std::size_t(las)];
        T q1 = qsfc(i, j), q2 = q[std::size_t(las)], z1 = zsfc(i, j), z2 = ght[std::size_t(las)];
        T rhobar = (p1 / (Rd * t1) + p2 / (Rd * t2)) * T(0.5);
        T qbar = (q1 + q2) * T(0.5);
        T dz = z2 - z1;
        if (dz > T(0.1)) intq(i, j) = intq(i, j) + g * qbar * rhobar / (T(1.) + qbar) * dz;
        for (int k = las - 1; k >= kts + 1; --k) pd[std::size_t(k)] = p[std::size_t(k)] - intq(i, j);
      }
      pd[std::size_t(kts)] = psfc(i, j) - intq(i, j);
      if (upside_down) {
        pd_out(i, kts, j) = pd[std::size_t(kts)];
        for (int k = kts + 1; k <= kte; ++k) pd_out(i, kte + 2 - k, j) = pd[std::size_t(k)];
      } else {
        for (int k = kts; k <= kte; ++k) pd_out(i, k, j) = pd[std::size_t(k)];
      }
    }
  }
}

// sfcprs3 over metgrid levels (kte = nlev)
template <typename T>
void sfcprs3(const F3<T>& height, const F3<T>& p, const F2<T>& ter, const F2<T>& slp, F2<T>& psfc, int nlev,
             const Dims& d, const SkipCtx& sk) {
  const int kts = 1, kte = nlev;
  for (int j = d.jts(); j <= std::min(d.jde - 1, d.jte()); ++j)
    for (int i = d.its(); i <= std::min(d.ide - 1, d.ite()); ++i) {
      if (skip_middle(d, i, j, sk.width, sk.hold_ups)) continue;
      if (ter(i, j) < T(50)) {
        psfc(i, j) = slp(i, j) + (p(i, 2, j) - p(i, 3, j)) / (height(i, 2, j) - height(i, 3, j)) * ter(i, j);
        continue;
      }
      bool found = false;
      T zl, zu, pl, pu, zm;
      for (int k = kts + 1; k <= kte - 2; ++k)
        if (height(i, k, j) <= ter(i, j) && height(i, k + 1, j) > ter(i, j)) {
          zl = height(i, k, j); zu = height(i, k + 1, j); zm = ter(i, j); pl = p(i, k, j); pu = p(i, k + 1, j);
          psfc(i, j) = std::exp((std::log(pl) * (zm - zu) + std::log(pu) * (zl - zm)) / (zl - zu));
          found = true;
          break;
        }
      if (!found) {
        if (slp(i, j) >= p(i, 2, j)) {
          zl = T(0.); zu = height(i, 3, j); zm = ter(i, j); pl = slp(i, j); pu = p(i, 3, j);
          psfc(i, j) = std::exp((std::log(pl) * (zm - zu) + std::log(pu) * (zl - zm)) / (zl - zu));
          found = true;
        } else {
          for (int k = kts + 1; k <= kte - 3; ++k)
            if (slp(i, j) >= p(i, k + 1, j) && slp(i, j) < p(i, k, j)) {
              zl = T(0.); zu = height(i, k + 1, j); zm = ter(i, j); pl = slp(i, j); pu = p(i, k + 1, j);
              psfc(i, j) = std::exp((std::log(pl) * (zm - zu) + std::log(pu) * (zl - zm)) / (zl - zu));
              found = true;
              break;
            }
        }
      }
      if (!found) throw std::runtime_error("real: problems with the surface pressure computation");
    }
}

// ---------------------------------------------------------------- levels / compute_eta
// levels: eta(0:nlev) in the default REAL kind
template <typename T>
void levels(int nlev, T ptop, T* eta /*0..nlev*/, T dzmax, T dzbot, T dzstretch_s, T dzstretch_u, T r_d, T g) {
  std::vector<T> zup(static_cast<std::size_t>(nlev + 1)), pup(static_cast<std::size_t>(nlev + 1));
  const T tt = T(290.);
  T ztop = r_d * tt / g * std::log(T(1.e5) / ptop);
  T dz = dzbot;
  zup[1] = dz;
  pup[1] = T(1.e5) * std::exp(-g * zup[1] / r_d / tt);
  eta[0] = T(1.0);
  eta[1] = (pup[1] - ptop) / (T(1.e5) - ptop);
  int isave = 1;
  for (int i = 1; i <= nlev - 1; ++i) {
    T a = dzstretch_u + (dzstretch_s - dzstretch_u) * std::max((dzmax * T(0.5) - dz) / (dzmax * T(0.5)), T(0.));
    dz = a * dz;
    T dztest = (ztop - zup[std::size_t(isave)]) / T(nlev - isave);
    if (dztest < dz) break;
    isave = i + 1;
    zup[std::size_t(i + 1)] = zup[std::size_t(i)] + dz;
    pup[std::size_t(i + 1)] = T(1.e5) * std::exp(-g * zup[std::size_t(i + 1)] / r_d / tt);
    eta[i + 1] = (pup[std::size_t(i + 1)] - ptop) / (T(1.e5) - ptop);
    if (i == nlev - 1) throw std::runtime_error("real: not enough eta levels to reach p_top");
  }
  dz = (ztop - zup[std::size_t(isave)]) / T(nlev - isave);
  if (dz > T(1.5) * dzmax) throw std::runtime_error("real: upper levels may be too thick");
  for (int i = isave; i <= nlev - 1; ++i) {
    zup[std::size_t(i + 1)] = zup[std::size_t(i)] + dz;
    pup[std::size_t(i + 1)] = T(1.e5) * std::exp(-g * zup[std::size_t(i + 1)] / r_d / tt);
    eta[i + 1] = (pup[std::size_t(i + 1)] - ptop) / (T(1.e5) - ptop);
  }
  eta[nlev] = T(0.);
}

struct EtaOpts {
  std::vector<double> eta_levels;  // namelist eta_levels (empty: computed)
  int auto_levels_opt = 2;
  double max_dz = 1000., dzbot = 50., dzstretch_s = 1.3, dzstretch_u = 1.1;
};

// compute_eta: znw(kts..kte). Only the namelist and auto_levels_opt = 2 branches (WRF's default).
template <typename T>
void compute_eta(F1<T>& znw, const EtaOpts& eo, T p_top, int kde) {
  const int kds = 1, kte = kde;
  if (!eo.eta_levels.empty()) {
    std::vector<T> el(static_cast<std::size_t>(kde + 1));
    for (int k = 1; k <= kde && k <= int(eo.eta_levels.size()); ++k) el[std::size_t(k)] = T(eo.eta_levels[std::size_t(k - 1)]);
    if (std::abs(el[1] - T(1.)) < T(0.0000001) && std::abs(el[std::size_t(kde)] - T(0.)) < T(0.0000001)) {
      for (int k = kds + 1; k <= kde - 1; ++k) znw(k) = el[std::size_t(k)];
      znw(1) = T(1.);
      znw(kde) = T(0.);
    } else if (std::abs(el[std::size_t(kde)] - T(1.)) < T(0.0000001) && std::abs(el[1] - T(0.)) < T(0.0000001)) {
      for (int k = kds + 1; k <= kde - 1; ++k) znw(k) = el[std::size_t(kde + 1 - k)];
      znw(1) = T(1.);
      znw(kde) = T(0.);
    } else {
      throw std::runtime_error("real: first eta level should be 1.0 and the last 0.0 in namelist");
    }
    for (int k = kds; k <= kde - 1; ++k)
      if (znw(k) <= znw(k + 1)) throw std::runtime_error("real: eta_levels not monotonic");
    return;
  }
  if (eo.auto_levels_opt != 2) throw std::runtime_error("real: only auto_levels_opt = 2 is ported");
  std::vector<T> eta(static_cast<std::size_t>(kte));  // levels writes znw(1..kte) as eta(0..kte-1)
  levels<T>(kte - 1, p_top, eta.data(), T(eo.max_dz), T(eo.dzbot), T(eo.dzstretch_s), T(eo.dzstretch_u),
            Const<T>::r_d, Const<T>::g);
  for (int k = 1; k <= kte; ++k) znw(k) = eta[std::size_t(k - 1)];
  // the rest of the auto_levels_opt = 2 branch only prints diagnostics
}

// nest_init_utils.F compute_vcoord_1d_coeffs (without the per-column monotonicity check)
template <typename T>
void compute_vcoord_1d_coeffs(T etac, const F1<T>& znw, int hybrid_opt, T p_top, int kde, F1<T>& znu, F1<T>& c1f,
                              F1<T>& c2f, F1<T>& c3f, F1<T>& c4f, F1<T>& c1h, F1<T>& c2h, F1<T>& c3h, F1<T>& c4h) {
  const int kds = 1, kts = 1, kte = kde;
  const T p1000mb = Const<T>::p1000mb;
  for (int k = kts; k <= kte; ++k) {
    if (hybrid_opt == 0 || hybrid_opt == 1) {
      c3f(k) = znw(k);
    } else if (hybrid_opt == 2) {
      T B1 = T(2.) * sq(etac) * (T(1.) - etac);
      T B2 = -etac * (T(4.) - T(3.) * etac - cube(etac));
      T B3 = T(2.) * (T(1.) - cube(etac));
      T B4 = -(T(1.) - sq(etac));
      T B5 = pow4(T(1.) - etac);
      c3f(k) = (B1 + B2 * znw(k) + B3 * sq(znw(k)) + B4 * cube(znw(k))) / B5;
      if (znw(k) < etac) c3f(k) = T(0.);
      if (k == kds) c3f(k) = T(1.);
      else if (k == kde) c3f(k) = T(0.);
    } else {
      throw std::runtime_error("real: hybrid_opt must be 0, 1 or 2");
    }
  }
  for (int k = 1; k <= kde; ++k) c4f(k) = (znw(k) - c3f(k)) * (p1000mb - p_top);
  for (int k = 1; k <= kde - 1; ++k) {
    znu(k) = (znw(k + 1) + znw(k)) * T(0.5);
    c3h(k) = (c3f(k + 1) + c3f(k)) * T(0.5);
    c4h(k) = (znu(k) - c3h(k)) * (p1000mb - p_top);
  }
  for (int k = kds + 1; k <= kde - 1; ++k) c1f(k) = (c3h(k) - c3h(k - 1)) / (znu(k) - znu(k - 1));
  c1f(kds) = T(1.);
  c1f(kde) = (hybrid_opt == 0 || hybrid_opt == 1) ? T(1.) : T(0.);
  for (int k = kds; k <= kde; ++k) c2f(k) = (T(1.) - c1f(k)) * (p1000mb - p_top);
  for (int k = 1; k <= kde - 1; ++k) {
    c1h(k) = (c3f(k + 1) - c3f(k)) / (znw(k + 1) - znw(k));
    c2h(k) = (T(1.) - c1h(k)) * (p1000mb - p_top);
  }
}

}  // namespace cppwrf::init_real
