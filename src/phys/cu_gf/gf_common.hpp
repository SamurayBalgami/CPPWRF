#pragma once
// Grell-Freitas cumulus (cu_physics = 3), WRF 4.6.1 phys/module_cu_gf_deep.F: the routines shared by
// the deep (CUP_gf) and shallow (CUP_gf_sh) schemes, for one column. Arrays keep the Fortran 1-based
// level index (k = 1 is the lowest level), sized KM and zero-initialised by the callers. Every
// statement follows the Fortran in order: x**2 is x*x, x**real is pow, Fortran MAX/MIN are fmax/fmin.
// The module constants (g = 9.81, cp = 1004, xlv = 2.5e6, r_v = 461) are the scheme's own, not WRF's.
#include <cmath>
#include <limits>
#include <type_traits>

#include "core/backend.hpp"

#if defined(__x86_64__) && !defined(CPPWRF_GF_NO_LIBMVEC)
#include <emmintrin.h>
extern "C" __m128d _ZGVbN2vv_pow(__m128d, __m128d);  // glibc libmvec, SSE2 ABI
#define CPPWRF_GF_HAVE_LIBMVEC 1
#endif

namespace cppwrf::phys::cu_gf {

// pow for two loop iterations that gfortran (-O2 -ftree-vectorize) computes with glibc's 2-lane vector
// pow _ZGVbN2vv_pow (it may differ from the scalar pow in the last bit). Host, double only.
template <typename T>
CPPWRF_HD inline void pow_pair(T x0, T x1, T y, T& r0, T& r1) {
#if defined(CPPWRF_GF_HAVE_LIBMVEC) && !defined(__CUDA_ARCH__)
  if constexpr (std::is_same_v<T, double>) {
    alignas(16) double out[2];
    _mm_store_pd(out, _ZGVbN2vv_pow(_mm_set_pd(x1, x0), _mm_set1_pd(y)));
    r0 = out[0];
    r1 = out[1];
    return;
  }
#endif
  r0 = std::pow(x0, y);
  r1 = std::pow(x1, y);
}

// zu(k) = zubeg + fzu * kratio**(alpha-1) * (1-kratio)**(beta-1) for k = lo..hi, kratio = (p(k)-p0)/(p1-p0):
// pairs of iterations through the vector pow, a trailing odd iteration through the scalar pow (as WRF's
// object code does)
template <typename T>
CPPWRF_HD inline void zu_profile(const T* p, T p0, T p1, int lo, int hi, T zubeg, T fzu, T alpha, T beta, T* zu,
                                 bool add_zubeg) {
  const int n = hi - lo + 1;
  if (n <= 0) return;
  const int nv = n - n % 2;
  for (int i = 0; i < nv; i += 2) {
    const int k = lo + i;
    const T k0 = (p[k] - p0) / (p1 - p0), k1 = (p[k + 1] - p0) / (p1 - p0);
    T a0, a1, b0, b1;
    pow_pair(k0, k1, alpha - T(1.0), a0, a1);
    pow_pair(T(1.0) - k0, T(1.0) - k1, beta - T(1.0), b0, b1);
    zu[k] = add_zubeg ? zubeg + fzu * a0 * b0 : fzu * a0 * b0;
    zu[k + 1] = add_zubeg ? zubeg + fzu * a1 * b1 : fzu * a1 * b1;
  }
  for (int k = lo + nv; k <= hi; ++k) {
    const T kratio = (p[k] - p0) / (p1 - p0);
    const T a = std::pow(kratio, alpha - T(1.0)), b = std::pow(T(1.0) - kratio, beta - T(1.0));
    zu[k] = add_zubeg ? zubeg + fzu * a * b : fzu * a * b;
  }
}

inline constexpr int kMaxLevels = 100;
inline constexpr int KM = kMaxLevels + 12;  // room for the k+1..k+8 reads of the Fortran
inline constexpr int maxens3 = 16;

template <typename T>
struct GfK {  // module_cu_gf_deep / module_cu_gf_sh parameters
  static constexpr T g = T(9.81), cp = T(1004.), xlv = T(2.5e6), r_v = T(461.);
  static constexpr T tcrit = T(258.), c1 = T(.001), frh_thresh = T(.9), rh_thresh = T(.97), betajb = T(1.5);
  static constexpr T fluxtune = T(1.5), pgcd = T(1.), ccnclean = T(250.);
  static constexpr T c1_shal = T(0.), c0_shal = T(.001);
};

template <typename T> CPPWRF_HD inline T fmx(T a, T b) { return std::fmax(a, b); }
template <typename T> CPPWRF_HD inline T fmn(T a, T b) { return std::fmin(a, b); }

// Fortran MAXLOC(a(lo:hi),1) as an absolute index (first maximum)
template <typename T>
CPPWRF_HD inline int maxloc_abs(const T* a, int lo, int hi) {
  int m = lo;
  for (int k = lo + 1; k <= hi; ++k)
    if (a[k] > a[m]) m = k;
  return m;
}
template <typename T>
CPPWRF_HD inline T maxval(const T* a, int lo, int hi) {
  T m = a[lo];
  for (int k = lo + 1; k <= hi; ++k)
    if (a[k] > m) m = a[k];
  return m;
}

// satvap (hPa)
template <typename T>
CPPWRF_HD inline T satvap(T temp2) {
  const T temp = temp2 - T(273.155);
  if (temp < T(-20.)) {
    const T toot = T(273.16) / temp2;
    const T toto = T(1) / toot;
    const T eilog = -(T(9.09718) * (toot - T(1))) - T(3.56654) * (std::log(toot) / std::log(T(10.))) +
                    T(.876793) * (T(1) - toto) + (std::log(T(6.1071)) / std::log(T(10.)));
    return std::pow(T(10), eilog);
  }
  const T tsot = T(373.16) / temp2;
  const T ewlog = -(T(7.90298) * (tsot - T(1))) + T(5.02808) * (std::log(tsot) / std::log(T(10.)));
  const T ewlog2 = ewlog - T(1.3816e-07) * (std::pow(T(10), T(11.344) * (T(1) - (T(1) / tsot))) - T(1));
  const T ewlog3 = ewlog2 + T(.0081328) * (std::pow(T(10), -(T(3.49149) * (tsot - T(1)))) - T(1));
  const T ewlog4 = ewlog3 + (std::log(T(1013.246)) / std::log(T(10.)));
  return std::pow(T(10), ewlog4);
}

// cup_env with itest = -1 (z given; he, hes, qes from t, q, p)
template <typename T>
CPPWRF_HD inline void cup_env(const T* z, T* qes, T* he, T* hes, const T* t, const T* q, const T* p, int ierr,
                              int ktf) {
  if (ierr != 0) return;
  for (int k = 1; k <= ktf; ++k) {
    const T e = satvap(t[k]);
    qes[k] = T(0.622) * e / fmx(T(1.e-8), (p[k] - e));
    if (qes[k] <= T(1.E-16)) qes[k] = T(1.E-16);
    if (qes[k] < q[k]) qes[k] = q[k];
  }
  for (int k = 1; k <= ktf; ++k) {
    he[k] = T(9.81) * z[k] + T(1004.) * t[k] + T(2.5E06) * q[k];
    hes[k] = T(9.81) * z[k] + T(1004.) * t[k] + T(2.5E06) * qes[k];
    if (he[k] >= hes[k]) he[k] = hes[k];
  }
}

template <typename T>
CPPWRF_HD inline void cup_env_clev(const T* t, const T* qes, const T* q, const T* he, const T* hes, const T* z,
                                   const T* p, T* qes_cup, T* q_cup, T* he_cup, T* hes_cup, T* z_cup, T* p_cup,
                                   T* gamma_cup, T* t_cup, T psur, int ierr, T z1, int ktf) {
  using K = GfK<T>;
  for (int k = 1; k <= ktf; ++k) {
    qes_cup[k] = 0; q_cup[k] = 0; hes_cup[k] = 0; he_cup[k] = 0;
    z_cup[k] = 0; p_cup[k] = 0; t_cup[k] = 0; gamma_cup[k] = 0;
  }
  if (ierr != 0) return;
  for (int k = 2; k <= ktf; ++k) {
    qes_cup[k] = T(.5) * (qes[k - 1] + qes[k]);
    q_cup[k] = T(.5) * (q[k - 1] + q[k]);
    hes_cup[k] = T(.5) * (hes[k - 1] + hes[k]);
    he_cup[k] = T(.5) * (he[k - 1] + he[k]);
    if (he_cup[k] > hes_cup[k]) he_cup[k] = hes_cup[k];
    z_cup[k] = T(.5) * (z[k - 1] + z[k]);
    p_cup[k] = T(.5) * (p[k - 1] + p[k]);
    t_cup[k] = T(.5) * (t[k - 1] + t[k]);
    gamma_cup[k] = (K::xlv / K::cp) * (K::xlv / (K::r_v * t_cup[k] * t_cup[k])) * qes_cup[k];
  }
  qes_cup[1] = qes[1];
  q_cup[1] = q[1];
  hes_cup[1] = T(9.81) * z1 + T(1004.) * t[1] + T(2.5e6) * qes[1];
  he_cup[1] = T(9.81) * z1 + T(1004.) * t[1] + T(2.5e6) * q[1];
  z_cup[1] = z1;
  p_cup[1] = psur;
  t_cup[1] = t[1];
  gamma_cup[1] = K::xlv / K::cp * (K::xlv / (K::r_v * t_cup[1] * t_cup[1])) * qes_cup[1];
}

// get_cloud_bc: mean of array(k22-2..k22) (fewer near the ground) plus add
template <typename T>
CPPWRF_HD inline T get_cloud_bc(const T* a, int k22, bool has_add, T add) {
  const int lo = k22 < 3 ? k22 : 3;
  T x = T(0.);
  for (int i = 1; i <= lo; ++i) x = x + a[k22 - i + 1];
  x = x / T(lo);
  if (has_add) x = x + add;
  return x;
}

template <typename T>
CPPWRF_HD inline void cup_kbcon(T cap_inc, int iloop_in, int& k22, int& kbcon, const T* he_cup, const T* hes_cup,
                                T& hkb, int& ierr, int kbmax, const T* p_cup, T cap_max, T ztexec, T zqexec,
                                const T* z_cup, T entr_rate, const T* heo, int imid) {
  using K = GfK<T>;
  T hcot[KM] = {};
  int iloop = iloop_in;
  kbcon = 1;
  if (cap_max > T(200) && imid == 1) iloop = 5;
  if (ierr != 0) return;
  int start_level = k22;
  kbcon = k22 + 1;
  if (iloop == 5) kbcon = k22;
  auto fill = [&]() {
    for (int k = 1; k <= start_level; ++k) hcot[k] = hkb;
    for (int k = start_level + 1; k <= kbmax + 3; ++k) {
      const T dz = z_cup[k] - z_cup[k - 1];
      hcot[k] = ((T(1.) - T(0.5) * entr_rate * dz) * hcot[k - 1] + entr_rate * dz * heo[k - 1]) /
                (T(1.) + T(0.5) * entr_rate * dz);
    }
  };
  fill();
  for (;;) {  // label 32
    const T hetest = hcot[kbcon];
    if (hetest < hes_cup[kbcon]) {  // label 31
      kbcon = kbcon + 1;
      if (kbcon > kbmax + 2) {
        if (iloop != 4) ierr = 3;
        return;
      }
      continue;
    }
    if (kbcon - k22 == 1) return;
    if (iloop == 5 && (kbcon - k22) <= 2) return;
    T pbcdif = -p_cup[kbcon] + p_cup[k22];
    T plus = fmx(T(25.), cap_max - T(iloop - 1) * cap_inc);
    if (iloop == 4) plus = cap_max;
    if (iloop == 5) plus = T(150.);
    if (iloop == 5 && cap_max > T(200)) pbcdif = -p_cup[kbcon] + cap_max;
    if (pbcdif <= plus) return;
    k22 = k22 + 1;
    kbcon = k22 + 1;
    const T x_add = K::xlv * zqexec + K::cp * ztexec;
    hkb = get_cloud_bc(he_cup, k22, true, x_add);
    start_level = k22;
    fill();
    if (iloop == 5) kbcon = k22;
    if (kbcon > kbmax + 2) {
      if (iloop != 4) ierr = 3;
      return;
    }
  }
}

template <typename T>
CPPWRF_HD inline int cup_minimi(const T* a, int ks, int kend, int ierr) {
  int kt = ks;
  if (ierr != 0) return kt;
  T x = a[ks];
  const int kstop = ks + 1 > kend ? ks + 1 : kend;
  for (int k = ks + 1; k <= kstop; ++k)
    if (a[k] < x) {
      x = a[k];
      kt = k;
    }
  return kt;
}

template <typename T>
CPPWRF_HD inline T cup_up_aa0(const T* z, const T* zu, const T* dby, const T* gamma_cup, const T* t_cup, int kbcon,
                              int ktop, int ierr, int ktf) {
  T aa0 = T(0.);
  for (int k = 2; k <= ktf; ++k) {
    if (ierr != 0) continue;
    if (k < kbcon) continue;
    if (k > ktop) continue;
    const T dz = z[k] - z[k - 1];
    const T da = zu[k] * dz * (T(9.81) / (T(1004.) * (t_cup[k]))) * dby[k - 1] / (T(1.) + gamma_cup[k]);
    aa0 = aa0 + fmx(T(0.), da);
    if (aa0 < T(0.)) aa0 = T(0.);
  }
  return aa0;
}

// get_zu_zd_pdf_fim for the drafts GF uses here: 0 = "UP", 1 = "SH2", 2 = "DOWN"
template <typename T>
CPPWRF_HD inline void get_zu_zd_pdf_fim(const T* p, T zubeg, int draft, int kb, int kt, T* zu, int kte, int ktf,
                                        int kpbli, int csum) {
  for (int k = 1; k <= kte; ++k) zu[k] = T(0.0);
  int kb_adj = kb > 2 ? kb : 2;
  if (draft == 0) {
    const T lev_start = fmn(T(.9), T(.4) + T(csum) * T(.013));
    kb_adj = kb > 2 ? kb : 2;
    T tunning = p[kt] + (p[kpbli] - p[kt]) * lev_start;
    tunning = fmn(T(0.9), (tunning - p[kb_adj]) / (p[kt] - p[kb_adj]));
    tunning = fmx(T(0.2), tunning);
    const T beta = T(1.3);
    const T alpha = (tunning * (beta - T(2.)) + T(1.)) / (T(1.) - tunning);
    const T fzu = std::tgamma(alpha + beta) / (std::tgamma(alpha) * std::tgamma(beta));
    zu_profile(p, p[kb_adj], p[kt], kb_adj, (kte < kt ? kte : kt), zubeg, fzu, alpha, beta, zu, true);
    const int hi = ktf < kt + 1 ? ktf : kt + 1;
    const T m = maxval(zu, 1, hi);
    if (m > T(0.))
      for (int k = 1; k <= hi; ++k) zu[k] = zu[k] / m;
    for (int k = maxloc_abs(zu, 1, kte); k >= 1; --k)
      if (zu[k] < T(1.e-6)) {
        kb_adj = k + 1;
        break;
      }
    kb_adj = kb_adj > 2 ? kb_adj : 2;
    for (int k = 1; k <= kb_adj - 1; ++k) zu[k] = T(0.);
  } else if (draft == 1) {
    T tunning = fmn(T(0.8), (p[kpbli] - p[kb_adj]) / (p[kt] - p[kb_adj]));
    tunning = fmx(T(0.2), tunning);
    const T beta = T(2.5);
    const T alpha = (tunning * (beta - T(2.)) + T(1.)) / (T(1.) - tunning);
    const T fzu = std::tgamma(alpha + beta) / (std::tgamma(alpha) * std::tgamma(beta));
    zu_profile(p, p[kb_adj], p[kt], kb_adj, (kte < kt ? kte : kt), zubeg, fzu, alpha, beta, zu, true);
    const int hi = ktf < kt + 1 ? ktf : kt + 1;
    const T m = maxval(zu, 1, hi);
    if (m > T(0.))
      for (int k = 1; k <= hi; ++k) zu[k] = zu[k] / m;
  } else {
    T tunning = p[kb];
    tunning = fmn(T(0.9), (tunning - p[1]) / (p[kt] - p[1]));
    tunning = fmx(T(0.2), tunning);
    const T beta = T(4.);
    const T alpha = (tunning * (beta - T(2.)) + T(1.)) / (T(1.) - tunning);
    T fzu = std::tgamma(alpha + beta) / (std::tgamma(alpha) * std::tgamma(beta));
    for (int k = 1; k <= kte; ++k) zu[k] = T(0.);
    zu_profile(p, p[1], p[kt], 2, (kt < ktf ? kt : ktf), T(0), fzu, alpha, beta, zu, false);
    const int hi = ktf < kt + 1 ? ktf : kt + 1;
    fzu = maxval(zu, 1, hi);
    if (fzu > T(0.))
      for (int k = 1; k <= hi; ++k) zu[k] = zu[k] / fzu;
    zu[1] = T(0.);
  }
}

// get_lateral_massflux; u arrays and lambau only for the deep scheme (has_u)
template <typename T>
CPPWRF_HD inline void get_lateral_massflux(int ktf, int kte, int ierr, int ktop, const T* zo_cup, const T* zuo,
                                           T* cd, T* entr_rate_2d, T* up_massentro, T* up_massdetro,
                                           T* up_massentr, T* up_massdetr, int k22, bool has_u, T* up_massentru,
                                           T* up_massdetru, T lambau) {
  for (int k = 1; k <= kte; ++k) {
    up_massentro[k] = 0; up_massdetro[k] = 0; up_massentr[k] = 0; up_massdetr[k] = 0;
  }
  if (has_u)
    for (int k = 1; k <= kte; ++k) {
      up_massentru[k] = 0;
      up_massdetru[k] = 0;
    }
  if (ierr != 0) return;
  const int kmx = maxloc_abs(zuo, 1, kte);
  for (int k = (2 > k22 + 1 ? 2 : k22 + 1); k <= kmx; ++k) {
    const T dz = zo_cup[k] - zo_cup[k - 1];
    up_massdetro[k - 1] = cd[k - 1] * dz * zuo[k - 1];
    up_massentro[k - 1] = zuo[k] - zuo[k - 1] + up_massdetro[k - 1];
    if (up_massentro[k - 1] < T(0.)) {
      up_massentro[k - 1] = T(0.);
      up_massdetro[k - 1] = zuo[k - 1] - zuo[k];
      if (zuo[k - 1] > T(0.)) cd[k - 1] = up_massdetro[k - 1] / (dz * zuo[k - 1]);
    }
    if (zuo[k - 1] > T(0.)) entr_rate_2d[k - 1] = (up_massentro[k - 1]) / (dz * zuo[k - 1]);
  }
  for (int k = kmx + 1; k <= ktop; ++k) {
    const T dz = zo_cup[k] - zo_cup[k - 1];
    up_massentro[k - 1] = entr_rate_2d[k - 1] * dz * zuo[k - 1];
    up_massdetro[k - 1] = zuo[k - 1] + up_massentro[k - 1] - zuo[k];
    if (up_massdetro[k - 1] < T(0.)) {
      up_massdetro[k - 1] = T(0.);
      up_massentro[k - 1] = zuo[k] - zuo[k - 1];
      if (zuo[k - 1] > T(0.)) entr_rate_2d[k - 1] = (up_massentro[k - 1]) / (dz * zuo[k - 1]);
    }
    if (zuo[k - 1] > T(0.)) cd[k - 1] = up_massdetro[k - 1] / (dz * zuo[k - 1]);
  }
  up_massdetro[ktop] = zuo[ktop];
  up_massentro[ktop] = T(0.);
  for (int k = ktop + 1; k <= ktf; ++k) {
    cd[k] = 0; entr_rate_2d[k] = 0; up_massentro[k] = 0; up_massdetro[k] = 0;
  }
  for (int k = 2; k <= ktf - 1; ++k) {
    up_massentr[k - 1] = up_massentro[k - 1];
    up_massdetr[k - 1] = up_massdetro[k - 1];
  }
  if (has_u)
    for (int k = 2; k <= ktf - 1; ++k) {
      up_massentru[k - 1] = up_massentro[k - 1] + lambau * up_massdetro[k - 1];
      up_massdetru[k - 1] = up_massdetro[k - 1] + lambau * up_massdetro[k - 1];
    }
}

// rates_up_pdf for 'deep' (shallow = false) and 'shallow' (zuo, ktop, kbcon, ktopdby, ierr updated)
template <typename T>
CPPWRF_HD inline void rates_up_pdf(bool shallow, int& ktop, int& ierr, const T* p_cup, const T* entr_rate_2d, T hkbo,
                                   const T* heo, const T* heso_cup, const T* z_cup, int kstabi, int k22, int& kbcon,
                                   int kte, int ktf, T* zuo, int kpbl, int& ktopdby, int csum) {
  const T zustart = T(.1), dbythresh = T(1.);
  T dby[KM] = {}, dbm[KM] = {}, hcot[KM] = {};
  const T beta_u = fmx(T(.1), T(.2) - T(csum) * T(.01));
  (void)beta_u;
  for (int k = 1; k <= kte; ++k) zuo[k] = T(0.);
  kbcon = kbcon > 2 ? kbcon : 2;
  if (ierr != 0) return;
  const int start_level = k22;
  zuo[start_level] = zustart;
  for (int k = start_level + 1; k <= kbcon; ++k) {
    const T dz = z_cup[k] - z_cup[k - 1];
    const T massent = dz * entr_rate_2d[k - 1] * zuo[k - 1];
    const T massdetr = dz * T(1.e-9) * zuo[k - 1];
    zuo[k] = zuo[k - 1] + massent - massdetr;
  }
  const T zubeg = zustart;
  if (!shallow) {
    ktop = 0;
    hcot[start_level] = hkbo;
    for (int k = start_level + 1; k <= ktf - 2; ++k) {
      const T dz = z_cup[k] - z_cup[k - 1];
      hcot[k] = ((T(1.) - T(0.5) * entr_rate_2d[k - 1] * dz) * hcot[k - 1] + entr_rate_2d[k - 1] * dz * heo[k - 1]) /
                (T(1.) + T(0.5) * entr_rate_2d[k - 1] * dz);
      if (k >= kbcon) dby[k] = dby[k - 1] + (hcot[k] - heso_cup[k]) * dz;
      if (k >= kbcon) dbm[k] = hcot[k] - heso_cup[k];
    }
    ktopdby = maxloc_abs(dby, 1, kte);
    int kfinalzu = ktf - 2;
    bool found = false;
    for (int k = maxloc_abs(dby, 1, kte) + 1; k <= ktf - 2; ++k)
      if (dby[k] < dbythresh * maxval(dby, 1, kte)) {
        kfinalzu = k - 1;
        found = true;
        break;
      }
    (void)found;
    ktop = kfinalzu;
    if (kfinalzu <= kbcon + 2) {
      ierr = 41;
      ktop = 0;
    } else {
      get_zu_zd_pdf_fim(p_cup, zubeg, 0, k22, kfinalzu, zuo, kte, ktf, kstabi, csum);
    }
  } else {
    if (ktop <= kbcon + 2) {
      ierr = 41;
      ktop = 0;
    } else {
      const int kfinalzu = ktop;
      ktopdby = ktop;
      get_zu_zd_pdf_fim(p_cup, zubeg, 1, k22, kfinalzu, zuo, kte, ktf, kpbl, ktopdby);
    }
  }
}

template <typename T>
CPPWRF_HD inline void cup_dd_moisture(const T* zd, const T* hcd, const T* hes_cup, T* qcd, const T* qes_cup, T* pwd,
                                      const T* q_cup, const T* z_cup, const T* dd_massentr, const T* dd_massdetr,
                                      int jmin, int& ierr, const T* gamma_cup, T& pwev, T& bu, T* qrcd, const T* q,
                                      int iloop, int ktf) {
  using K = GfK<T>;
  bu = T(0.);
  pwev = T(0.);
  for (int k = 1; k <= ktf; ++k) {
    qcd[k] = 0;
    qrcd[k] = 0;
    pwd[k] = 0;
  }
  if (ierr != 0) return;
  int k = jmin;
  T dz = z_cup[k + 1] - z_cup[k];
  qcd[k] = q_cup[k];
  T dh = hcd[k] - hes_cup[k];
  if (dh < T(0)) {
    qrcd[k] = (qes_cup[k] + (T(1.) / K::xlv) * (gamma_cup[k] / (T(1.) + gamma_cup[k])) * dh);
  } else {
    qrcd[k] = qes_cup[k];
  }
  pwd[jmin] = zd[jmin] * fmn(T(0.), qcd[k] - qrcd[k]);
  qcd[k] = qrcd[k];
  pwev = pwev + pwd[jmin];
  bu = dz * dh;
  for (int ki = jmin - 1; ki >= 1; --ki) {
    dz = z_cup[ki + 1] - z_cup[ki];
    const T denom = zd[ki + 1] - T(.5) * dd_massdetr[ki] + dd_massentr[ki];
    if (denom < T(1.e-8)) {
      ierr = 51;
      break;
    }
    qcd[ki] = (qcd[ki + 1] * zd[ki + 1] - T(.5) * dd_massdetr[ki] * qcd[ki + 1] + dd_massentr[ki] * q[ki]) /
              (zd[ki + 1] - T(.5) * dd_massdetr[ki] + dd_massentr[ki]);
    dh = hcd[ki] - hes_cup[ki];
    bu = bu + dz * dh;
    qrcd[ki] = qes_cup[ki] + (T(1.) / K::xlv) * (gamma_cup[ki] / (T(1.) + gamma_cup[ki])) * dh;
    T dqeva = qcd[ki] - qrcd[ki];
    if (dqeva > T(0.)) {
      dqeva = T(0.);
      qrcd[ki] = qcd[ki];
    }
    pwd[ki] = zd[ki] * dqeva;
    qcd[ki] = qrcd[ki];
    pwev = pwev + pwd[ki];
  }
  if ((pwev == T(0.)) && (iloop == 1)) ierr = 7;
  if (bu >= T(0) && iloop == 1) ierr = 7;
}

// cup_up_moisture with autoconv = 1 (iall = 0)
template <typename T>
CPPWRF_HD inline void cup_up_moisture(int& ierr, const T* z_cup, T* qc, T* qrc, T* pw, T& pwav, int kbcon, int ktop,
                                      const T* dby, T* clw_all, const T* q, const T* gamma_cup, const T* zu,
                                      const T* qes_cup, int k22, const T* qe_cup, const T* c1d, const T* t,
                                      const T* up_massentr, const T* up_massdetr, T& psum, T& psumh, int ktf) {
  using K = GfK<T>;
  T qch[KM] = {};
  pwav = T(0.);
  psum = T(0.);
  psumh = T(0.);
  for (int k = 1; k <= ktf; ++k) {
    pw[k] = 0;
    qc[k] = 0;
    if (ierr == 0) qc[k] = qe_cup[k];
    if (ierr == 0) qch[k] = qe_cup[k];
    clw_all[k] = 0;
    qrc[k] = 0;
  }
  if (ierr != 0) return;
  const int start_level = k22;
  const T qaver = get_cloud_bc(qe_cup, k22, false, T(0));
  qc[start_level] = qaver;
  qch[start_level] = qaver;
  for (int k = 1; k <= start_level - 1; ++k) {
    qc[k] = qe_cup[k];
    qch[k] = qe_cup[k];
  }
  T c0 = T(.004);
  for (int k = k22 + 1; k <= kbcon; ++k) {
    if (t[k] < T(273.15)) c0 = c0 * std::exp(T(0.07) * (t[k] - T(273.15)));
    qc[k] = (qc[k - 1] * zu[k - 1] - T(.5) * up_massdetr[k - 1] * qc[k - 1] + up_massentr[k - 1] * q[k - 1]) /
            (zu[k - 1] - T(.5) * up_massdetr[k - 1] + up_massentr[k - 1]);
    T qrch = qes_cup[k] + (T(1.) / K::xlv) * (gamma_cup[k] / (T(1.) + gamma_cup[k])) * dby[k];
    if (k < kbcon) qrch = qc[k];
    if (qc[k] > qrch) {
      const T dz = z_cup[k] - z_cup[k - 1];
      qrc[k] = (qc[k] - qrch) / (T(1.) + c0 * dz);
      pw[k] = c0 * dz * qrc[k] * zu[k];
      qc[k] = qrch + qrc[k];
      clw_all[k] = qrc[k];
    }
  }
  for (int k = kbcon + 1; k <= ktop; ++k) {
    c0 = T(.004);
    if (t[k] < T(273.15)) c0 = c0 * std::exp(T(0.07) * (t[k] - T(273.15)));
    const T denom = zu[k - 1] - T(.5) * up_massdetr[k - 1] + up_massentr[k - 1];
    if (denom < T(1.e-8)) {
      ierr = 51;
      break;
    }
    const T dz = z_cup[k] - z_cup[k - 1];
    const T qrch = qes_cup[k] + (T(1.) / K::xlv) * (gamma_cup[k] / (T(1.) + gamma_cup[k])) * dby[k];
    qc[k] = (qc[k - 1] * zu[k - 1] - T(.5) * up_massdetr[k - 1] * qc[k - 1] + up_massentr[k - 1] * q[k - 1]) /
            (zu[k - 1] - T(.5) * up_massdetr[k - 1] + up_massentr[k - 1]);
    qch[k] = (qch[k - 1] * zu[k - 1] - T(.5) * up_massdetr[k - 1] * qch[k - 1] + up_massentr[k - 1] * q[k - 1]) /
             (zu[k - 1] - T(.5) * up_massdetr[k - 1] + up_massentr[k - 1]);
    if (qc[k] <= qrch) qc[k] = qrch;
    if (qch[k] <= qrch) qch[k] = qrch;
    clw_all[k] = fmx(T(0.), qc[k] - qrch);
    qrc[k] = fmx(T(0.), (qc[k] - qrch));
    qrc[k] = (qc[k] - qrch) / (T(1.) + (c1d[k] + c0) * dz);
    pw[k] = c0 * dz * qrc[k] * zu[k];
    if (qrc[k] < T(0)) {
      qrc[k] = T(0.);
      pw[k] = T(0.);
    }
    qc[k] = qrc[k] + qrch;
    pwav = pwav + pw[k];
    psum = psum + clw_all[k] * zu[k] * dz;
  }
  for (int k = k22 + 1; k <= ktop; ++k) qc[k] = qc[k] - qrc[k];
}

}  // namespace cppwrf::phys::cu_gf
