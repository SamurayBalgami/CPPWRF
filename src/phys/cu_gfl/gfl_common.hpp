#pragma once

// Grell-Freitas-Li (GFL) shared constants and helper routines.
// Ported from WRF v4.8.0 phys/GFL/module_cu_gfl_deep.F90,
// phys/GFL/module_cu_gfl_sh.F90, and phys/GFL/WRF/module_cu_gfl_common.F90.

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "core/backend.hpp"

#ifndef CPPWRF_HOST_DEVICE
#define CPPWRF_HOST_DEVICE CPPWRF_HD
#endif

namespace cppwrf::phys::cu_gfl {

using std::abs;
using std::exp;
using std::log;
using std::max;
using std::min;
using std::pow;
using std::sin;
using std::sqrt;

#define R(x) static_cast<real>(x)

inline constexpr int kMaxLevels = 100;
inline constexpr int KM = kMaxLevels + 12;

template <typename real>
struct GflConstants {
  // Constants inside module_cu_gfl_deep / module_cu_gfl_sh (kind_phys literals)
  static constexpr real g = R(9.81);
  static constexpr real cp = R(1004.0);
  static constexpr real xlv = R(2.5e6);
  static constexpr real r_v = R(461.0);
  static constexpr real tcrit = R(258.0);
  static constexpr real fluxtune = R(1.5);
  static constexpr real pgcd = R(0.1);
  static constexpr real rad_deep = R(1.0);
  static constexpr real lambau_sh = R(2.0);
  static constexpr real ccnclean = R(250.0);

  // Deep module parameters
  static constexpr real c1 = R(0.003);
  static constexpr int irainevap = 1;
  static constexpr int use_fct = 1;
  static constexpr real frh_thresh = R(0.9);
  static constexpr real rh_thresh = R(0.97);
  static constexpr real betajb = R(1.2);
  static constexpr int use_excess = 0;
  static constexpr int autoconv = 1;
  static constexpr int aeroevap = 1;
  static constexpr real dx_thresh = R(6500.0);
  static constexpr real c0_deep_ocean = R(0.002);
  static constexpr real c0_deep_land = R(0.004);
  static constexpr real entr_rate_deep = R(7.e-5);
  static constexpr real entr_rate_deep_hr = R(1.e-4);
  static constexpr int maxens3 = 16;
  static constexpr int maxens_sh = 3;
  static constexpr bool melt_glac = true;
  static constexpr real t_0 = R(273.16);
  static constexpr real t_ice = R(250.16);
  static constexpr real xlf = R(0.333e6);

  // Shallow module parameters
  static constexpr real c1_shal = R(0.0005);
  static constexpr real c0_shal = R(0.001);
  static constexpr real entr_rate_sh = R(1.e-3);
};

// Match gfortran -O2 vectorised pow in get_zu_zd_pdf_fim: gfortran emits
// libmvec _ZGVbN2vv_pow for pairs of iterations and scalar pow for the odd
// remainder element.
template <typename real>
CPPWRF_HOST_DEVICE inline void pow_pair(real x0, real x1, real y, real& r0,
                                        real& r1) {
#if !defined(__CUDA_ARCH__) && defined(__x86_64__) && defined(__GLIBC__)
  if constexpr (std::is_same_v<real, double>) {
    alignas(16) double out[2];
    _mm_store_pd(out, _ZGVbN2vv_pow(_mm_set_pd(x1, x0), _mm_set1_pd(y)));
    r0 = out[0];
    r1 = out[1];
    return;
  }
#endif
  r0 = pow(x0, y);
  r1 = pow(x1, y);
}

template <typename real>
CPPWRF_HOST_DEVICE inline real tgamma_fn(real x) {
  if constexpr (std::is_same_v<real, float>) {
    return ::tgammaf(x);
  } else {
    return ::tgamma(x);
  }
}

template <typename real>
CPPWRF_HOST_DEVICE inline void zu_profile(real* zu, const real* p, int kstart,
                                          int kt_minus_1, real p_kstart,
                                          real dp, real zubeg, real fzu,
                                          real e1, real e2) {
  int k = kstart;
  for (; k + 1 <= kt_minus_1; k += 2) {
    real r0 = (p[k] - p_kstart) / dp;
    real r1 = (p[k + 1] - p_kstart) / dp;
    real a0, a1, b0, b1;
    pow_pair(r0, r1, e1, a0, a1);
    pow_pair(R(1.0) - r0, R(1.0) - r1, e2, b0, b1);
    zu[k] = zubeg + (fzu * a0) * b0;
    zu[k + 1] = zubeg + (fzu * a1) * b1;
  }
  if (k <= kt_minus_1) {
    real r = (p[k] - p_kstart) / dp;
    real a = pow(r, e1);
    real b = pow(R(1.0) - r, e2);
    zu[k] = zubeg + (fzu * a) * b;
  }
}

// ---------------------------------------------------------------------------
// satvap (lines 4963-4986 of module_cu_gfl_deep.F90)
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline real satvap(real temp2) {
  const real temp = temp2 - R(273.155);
  if (temp < R(-20.0)) {
    const real toot = R(273.16) / temp2;
    const real toto = R(1.0) / toot;
    const real eilog =
        -(R(9.09718) * (toot - R(1.0))) -
        R(3.56654) * (log(toot) / log(R(10.0))) +
        R(0.876793) * (R(1.0) - toto) + (log(R(6.1071)) / log(R(10.0)));
    return pow(R(10.0), eilog);
  }
  const real tsot = R(373.16) / temp2;
  const real ewlog = -(R(7.90298) * (tsot - R(1.0))) +
                     R(5.02808) * (log(tsot) / log(R(10.0)));
  const real ewlog2 =
      ewlog - R(1.3816e-07) *
                  (pow(R(10.0), R(11.344) * (R(1.0) - (R(1.0) / tsot))) -
                   R(1.0));
  const real ewlog3 =
      ewlog2 +
      R(0.0081328) * (pow(R(10.0), -(R(3.49149) * (tsot - R(1.0)))) - R(1.0));
  const real ewlog4 = ewlog3 + (log(R(1013.246)) / log(R(10.0)));
  return pow(R(10.0), ewlog4);
}

// ---------------------------------------------------------------------------
// cup_env (lines 3180-3324 of module_cu_gfl_deep.F90)
// 1-indexed arrays [1..kte]
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void cup_env(int ktf, int kte, int ierr, int itest,
                                       real psur, real z1, const real* p,
                                       const real* t, real* z, real* he,
                                       real* hes, real* qes, const real* q) {
  if (ierr != 0) return;
  real tv[KM];
  for (int k = 1; k <= ktf; ++k) {
    real e = satvap(t[k]);
    qes[k] = R(0.622) * e / max(R(1.e-8), (p[k] - e));
    if (qes[k] <= R(1.e-16)) qes[k] = R(1.e-16);
    if (qes[k] < q[k]) qes[k] = q[k];
    tv[k] = t[k] + R(0.608) * q[k] * t[k];
  }

  if (itest == 1 || itest == 0) {
    z[1] =
        max(R(0.0), z1) - (log(p[1]) - log(psur)) * R(287.0) * tv[1] / R(9.81);
    for (int k = 2; k <= ktf; ++k) {
      real tvbar = R(0.5) * tv[k] + R(0.5) * tv[k - 1];
      z[k] =
          z[k - 1] - (log(p[k]) - log(p[k - 1])) * R(287.0) * tvbar / R(9.81);
    }
  } else if (itest == 2) {
    for (int k = 1; k <= ktf; ++k) {
      z[k] = (he[k] - R(1004.0) * t[k] - R(2.5e6) * q[k]) / R(9.81);
      z[k] = max(R(1.e-3), z[k]);
    }
  }

  for (int k = 1; k <= ktf; ++k) {
    if (itest <= 0)
      he[k] = R(9.81) * z[k] + R(1004.0) * t[k] + R(2.5e6) * q[k];
    hes[k] = R(9.81) * z[k] + R(1004.0) * t[k] + R(2.5e6) * qes[k];
    if (he[k] >= hes[k]) he[k] = hes[k];
  }
  (void)kte;
}

// ---------------------------------------------------------------------------
// cup_env_clev (lines 3346-3429 of module_cu_gfl_deep.F90)
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void cup_env_clev(
    int ktf, int kte, int ierr, real psur, real z1, const real* t,
    const real* qes, const real* q, const real* he, const real* hes,
    const real* z, const real* p, real* qes_cup, real* q_cup, real* he_cup,
    real* hes_cup, real* z_cup, real* p_cup, real* gamma_cup, real* t_cup) {
  using C = GflConstants<real>;
  for (int k = 1; k <= ktf; ++k) {
    qes_cup[k] = R(0.0);
    q_cup[k] = R(0.0);
    hes_cup[k] = R(0.0);
    he_cup[k] = R(0.0);
    z_cup[k] = R(0.0);
    p_cup[k] = R(0.0);
    t_cup[k] = R(0.0);
    gamma_cup[k] = R(0.0);
  }
  if (ierr != 0) return;

  for (int k = 2; k <= ktf; ++k) {
    qes_cup[k] = R(0.5) * (qes[k - 1] + qes[k]);
    q_cup[k] = R(0.5) * (q[k - 1] + q[k]);
    hes_cup[k] = R(0.5) * (hes[k - 1] + hes[k]);
    he_cup[k] = R(0.5) * (he[k - 1] + he[k]);
    if (he_cup[k] > hes_cup[k]) he_cup[k] = hes_cup[k];
    z_cup[k] = R(0.5) * (z[k - 1] + z[k]);
    p_cup[k] = R(0.5) * (p[k - 1] + p[k]);
    t_cup[k] = R(0.5) * (t[k - 1] + t[k]);
    gamma_cup[k] = (C::xlv / C::cp) *
                   (C::xlv / (C::r_v * t_cup[k] * t_cup[k])) * qes_cup[k];
  }
  qes_cup[1] = qes[1];
  q_cup[1] = q[1];
  hes_cup[1] = R(9.81) * z1 + R(1004.0) * t[1] + R(2.5e6) * qes[1];
  he_cup[1] = R(9.81) * z1 + R(1004.0) * t[1] + R(2.5e6) * q[1];
  z_cup[1] = z1;
  p_cup[1] = psur;
  t_cup[1] = t[1];
  gamma_cup[1] = (C::xlv / C::cp) *
                 (C::xlv / (C::r_v * t_cup[1] * t_cup[1])) * qes_cup[1];
  (void)kte;
}

// ---------------------------------------------------------------------------
// get_cloud_bc (lines 4989-5012 of module_cu_gfl_deep.F90)
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void get_cloud_bc(int kte, const real* array,
                                            real& x_aver, int k22, real add) {
  constexpr int order_aver = 3;
  int local_order_aver = min(k22, order_aver);
  x_aver = R(0.0);
  for (int i = 1; i <= local_order_aver; ++i) {
    x_aver += array[k22 - i + 1];
  }
  x_aver = x_aver / real(local_order_aver);
  x_aver = x_aver + add;
  (void)kte;
}

// ---------------------------------------------------------------------------
// cup_direction2 (lines 2609-2688 of module_cu_gfl_deep.F90)
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void cup_direction2(int ktf, int kte, int ierr,
                                              int imass, const real* us,
                                              const real* vs, real& maxhs,
                                              real& dir) {
  if (imass == 1) maxhs = R(1.0);
  dir = R(0.0);
  if (ierr != 0) return;

  constexpr real pi = R(3.141592653589793);
  constexpr real rad = pi / R(180.0);

  real diff = R(22.5);
  int id = 0;
  real uust = R(0.0), vvst = R(0.0);
  for (int k = 2; k <= 10; ++k) {
    uust += us[k];
    vvst += vs[k];
  }
  uust /= R(9.0);
  vvst /= R(9.0);
  vvst = max(abs(vvst), R(1.e-10));
  real datan = std::atan(uust / vvst);
  real ddir = datan / pi * R(180.0);
  if (ddir < R(0.0)) ddir = -ddir;
  if (uust > R(0.0) && vvst < R(0.0)) ddir = R(180.0) - ddir;
  if (uust < R(0.0) && vvst < R(0.0)) ddir = R(180.0) + ddir;
  if (uust < R(0.0) && vvst > R(0.0)) ddir = R(360.0) - ddir;
  real dir_work = ddir;

  for (int ia = 1; ia <= 16; ++ia) {
    real ddira = real(ia - 1) * diff;
    real ddirb = ddira - R(0.5) * diff;
    if (ddirb < R(0.0)) ddirb = R(360.0) - R(0.5) * diff;
    real ddirc = ddira + R(0.5) * diff;
    if (ia == 1) {
      if (dir_work >= ddirb && dir_work <= R(360.0)) {
        dir_work = ddira * rad;
        id = ia;
      } else if (dir_work <= ddirc && dir_work >= R(0.0)) {
        dir_work = ddira * rad;
        id = ia;
      }
    } else {
      if (dir_work > ddirb && dir_work <= ddirc) {
        dir_work = ddira * rad;
        id = ia;
      }
    }
  }
  dir = dir_work;
  (void)id;
  (void)ktf;
  (void)kte;
}

// ---------------------------------------------------------------------------
// cup_maximi (lines 4002-4062 of module_cu_gfl_deep.F90)
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void cup_maximi(const real* array, int ks, int ke,
                                          int ierr, int& maxx) {
  maxx = ks;
  if (ierr != 0) return;
  real x = array[ks];
  for (int k = ks; k <= ke; ++k) {
    real xar = array[k];
    if (xar >= x) {
      x = xar;
      maxx = k;
    }
  }
}

// ---------------------------------------------------------------------------
// cup_minimi (lines 4065-4120 of module_cu_gfl_deep.F90)
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void cup_minimi(const real* array, int ks, int kend,
                                          int ierr, int& kt) {
  kt = ks;
  if (ierr != 0) return;
  real x = array[ks];
  int kstop = max(ks + 1, kend);
  for (int k = ks + 1; k <= kstop; ++k) {
    if (array[k] < x) {
      x = array[k];
      kt = k;
    }
  }
}

// ---------------------------------------------------------------------------
// cup_kbcon (lines 3842-3998 of module_cu_gfl_deep.F90)
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void cup_kbcon(
    real cap_inc, int iloop_in, int& k22, int& kbcon, const real* he_cup,
    const real* hes_cup, real& hkb, int& ierr, int kbmax, const real* p_cup,
    real cap_max, real ztexec, real zqexec, int kte, const real* z_cup,
    real entr_rate, const real* heo, int imid) {
  using C = GflConstants<real>;
  real hcot[KM] = {};
  int iloop = iloop_in;
  kbcon = 1;
  if (cap_max > R(200.0) && imid == 1) iloop = 5;
  if (ierr != 0) return;

  int start_level = k22;
  kbcon = k22 + 1;
  if (iloop == 5) kbcon = k22;

  auto fill_hcot = [&]() {
    for (int k = 1; k <= start_level; ++k) hcot[k] = hkb;
    for (int k = start_level + 1; k <= kbmax + 3; ++k) {
      real dz = z_cup[k] - z_cup[k - 1];
      hcot[k] = ((R(1.0) - R(0.5) * entr_rate * dz) * hcot[k - 1] +
                 entr_rate * dz * heo[k - 1]) /
                (R(1.0) + R(0.5) * entr_rate * dz);
    }
  };
  fill_hcot();

  for (;;) {
    real hetest = hcot[kbcon];
    if (hetest < hes_cup[kbcon]) {
      kbcon = kbcon + 1;
      if (kbcon > kbmax + 2) {
        if (iloop != 4) ierr = 3;
        return;
      }
      continue;
    }

    if (kbcon - k22 == 1) return;
    if (iloop == 5 && (kbcon - k22) <= 2) return;
    real pbcdif = -p_cup[kbcon] + p_cup[k22];
    real plus = max(R(25.0), cap_max - real(iloop - 1) * cap_inc);
    if (iloop == 4) plus = cap_max;
    if (iloop == 5) plus = R(150.0);
    if (iloop == 5 && cap_max > R(200.0)) pbcdif = -p_cup[kbcon] + cap_max;
    if (pbcdif <= plus) {
      return;
    } else if (pbcdif > plus) {
      k22 = k22 + 1;
      kbcon = k22 + 1;
      real x_add = C::xlv * zqexec + C::cp * ztexec;
      get_cloud_bc(kte, he_cup, hkb, k22, x_add);
      start_level = k22;
      fill_hcot();
      if (iloop == 5) kbcon = k22;
      if (kbcon > kbmax + 2) {
        if (iloop != 4) ierr = 3;
        return;
      }
    }
  }
}

// ---------------------------------------------------------------------------
// cup_up_he (lines 4520-4582 in 4.6.1 / 4470 in GFL if used; GFL deep/sh inline)
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void cup_up_he(int ktf, int kte, int k22, real hkb,
                                         const real* z_cup, const real* cd,
                                         const real* entr_rate_2d,
                                         const real* he, const real* he_cup,
                                         real* hc, int kbcon, int ierr,
                                         real* dby, const real* hes_cup,
                                         const real* up_massentr,
                                         const real* up_massdetr) {
  for (int k = 1; k <= kte; ++k) {
    hc[k] = R(0.0);
    dby[k] = R(0.0);
  }
  if (ierr != 0) return;
  for (int k = 1; k <= k22; ++k) {
    hc[k] = he_cup[k];
  }
  hc[k22] = hkb;
  for (int k = 2; k <= k22; ++k) {
    dby[k] = hc[k] - hes_cup[k];
  }
  for (int k = k22 + 1; k <= ktf; ++k) {
    hc[k] = (hc[k - 1] * R(1.0) - R(0.5) * up_massdetr[k - 1] * hc[k - 1] +
             up_massentr[k - 1] * he[k - 1]) /
            (R(1.0) - R(0.5) * up_massdetr[k - 1] + up_massentr[k - 1]);
    dby[k] = hc[k] - hes_cup[k];
  }
  (void)z_cup;
  (void)cd;
  (void)entr_rate_2d;
  (void)kbcon;
}

// ---------------------------------------------------------------------------
// cup_up_aa0 (lines 4123-4196 of module_cu_gfl_deep.F90)
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void cup_up_aa0(int ktf, int kte, real& aa0,
                                          const real* z, const real* zu,
                                          const real* dby,
                                          const real* gamma_cup,
                                          const real* t_cup, int kbcon,
                                          int ktop, int ierr) {
  aa0 = R(0.0);
  if (ierr != 0) return;
  for (int k = 2; k <= ktf; ++k) {
    if (k < kbcon) continue;
    if (k > ktop) continue;
    real dz = z[k] - z[k - 1];
    real da = zu[k] * dz * (R(9.81) / (R(1004.0) * t_cup[k])) * dby[k - 1] /
              (R(1.0) + gamma_cup[k]);
    aa0 = aa0 + max(R(0.0), da);
    if (aa0 < R(0.0)) aa0 = R(0.0);
  }
  (void)kte;
}

template <typename T>
CPPWRF_HOST_DEVICE inline int maxloc_abs(const T* a, int lo, int hi) {
  int m = lo;
  for (int k = lo + 1; k <= hi; ++k) {
    if (a[k] > a[m]) m = k;
  }
  return m;
}

template <typename real>
CPPWRF_HOST_DEVICE inline real maxval_range(const real* a, int lo, int hi) {
  real m = a[lo];
  for (int k = lo + 1; k <= hi; ++k) {
    if (a[k] > m) m = a[k];
  }
  return m;
}

// ---------------------------------------------------------------------------
// get_inversion_layers (lines 5486-5594 of module_cu_gfl_deep.F90)
// In GFL v4.8.0: sec_deriv(:) = 0. at entry, ilev = max(4, kstart + 1),
// and sec_deriv(:) = 1.e9 is inside the column loop.
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void get_inversion_layers(
    int ierr, const real* p_cup, const real* t_cup, const real* z_cup,
    int kstart, int kend, int kte, int& kinv1, int& kinv2) {
  int kinv[KM] = {};
  for (int k = 1; k <= kte; ++k) kinv[k] = 1;
  kinv1 = 1;
  kinv2 = 1;
  if (ierr != 0) return;

  real first_deriv[KM] = {}, sec_deriv[KM] = {}, dtempdz[KM] = {};
  int kend_p3 = kend + 3;
  for (int k = 2; k <= kend_p3 + 4; ++k) {
    first_deriv[k] =
        (t_cup[k + 1] - t_cup[k - 1]) / (z_cup[k + 1] - z_cup[k - 1]);
    dtempdz[k] = first_deriv[k];
  }
  for (int k = 3; k <= kend_p3 + 3; ++k) {
    sec_deriv[k] = (first_deriv[k + 1] - first_deriv[k - 1]) /
                   (z_cup[k + 1] - z_cup[k - 1]);
    sec_deriv[k] = abs(sec_deriv[k]);
  }

  int ilev = max(4, kstart + 1);
  int ix = 1;
  int k = ilev;
  while (ilev < kend_p3) {
    for (int kk = k; kk <= kend_p3 + 2; ++kk) {
      if (sec_deriv[kk] < sec_deriv[kk + 1] &&
          sec_deriv[kk] < sec_deriv[kk - 1]) {
        kinv[ix] = kk;
        ix = min(5, ix + 1);
        ilev = kk + 1;
        break;
      }
      ilev = kk + 1;
    }
    k = ilev;
  }

  int kadd = 0;
  int ken = maxloc_abs(kinv, 1, kte);
  for (int kq = 1; kq <= ken; ++kq) {
    int kk = kinv[kq + kadd];
    if (kk == 1) break;
    if (dtempdz[kk] < dtempdz[kk - 1] && dtempdz[kk] < dtempdz[kk + 1]) {
      kadd = kadd + 1;
      for (int kj = kq; kj <= ken; ++kj) {
        if (kinv[kj + kadd] > 1) kinv[kj] = kinv[kj + kadd];
        if (kinv[kj + kadd] == 1) kinv[kj] = 1;
      }
    }
  }

  int ken2 = maxloc_abs(kinv, 1, kte);
  for (int kk = 1; kk <= kte; ++kk) sec_deriv[kk] = R(1.e9);
  for (int kk = 1; kk <= ken2; ++kk) {
    real dp = p_cup[kinv[kk]] - p_cup[kstart];
    sec_deriv[kk] = abs(dp) - R(100.0);
  }
  int k800 = 1;
  for (int kk = 2; kk <= kte; ++kk) {
    if (abs(sec_deriv[kk]) < abs(sec_deriv[k800])) k800 = kk;
  }

  for (int kk = 1; kk <= kte; ++kk) sec_deriv[kk] = R(1.e9);
  for (int kk = 1; kk <= ken2; ++kk) {
    real dp = p_cup[kinv[kk]] - p_cup[kstart];
    sec_deriv[kk] = abs(dp) - R(300.0);
  }
  int k550 = 1;
  for (int kk = 2; kk <= kte; ++kk) {
    if (abs(sec_deriv[kk]) < abs(sec_deriv[k550])) k550 = kk;
  }

  kinv[1] = kinv[k800];
  kinv[2] = kinv[k550];
  kinv1 = kinv[1];
  kinv2 = kinv[2];
}

// ---------------------------------------------------------------------------
// get_zu_zd_pdf_fim (lines 5154-5412 of module_cu_gfl_deep.F90)
// draft: 1 = deep, 2 = shallow, 3 = mid, 4 = downdraft
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void get_zu_zd_pdf_fim(int kklev, const real* p,
                                                 real rand_vmas, real zubeg,
                                                 int ipr, int xland,
                                                 const real* zuh2, int draft,
                                                 int& ierr, int kb, int kt,
                                                 real* zu, int kts, int kte,
                                                 int ktf, real max_mass,
                                                 int kpbli, int csum,
                                                 int pmin_lev) {
  static constexpr real alpha[31] = {
      R(0.0),      R(3.699999), R(3.699999), R(3.699999), R(3.699999),
      R(3.024999), R(2.559999), R(2.249999), R(2.028571), R(1.862500),
      R(1.733333), R(1.630000), R(1.545454), R(1.475000), R(1.415385),
      R(1.364286), R(1.320000), R(1.281250), R(1.247059), R(1.216667),
      R(1.189474), R(1.165000), R(1.142857), R(1.122727), R(1.104348),
      R(1.087500), R(1.075000), R(1.075000), R(1.075000), R(1.075000),
      R(1.075000)};
  static constexpr real g_alpha[31] = {
      R(0.0),       R(4.170645),  R(4.170645),  R(4.170645),  R(4.170645),
      R(2.046925),  R(1.387837),  R(1.133003),  R(1.012418),  R(0.9494680),
      R(0.9153771), R(0.8972442), R(0.8885444), R(0.8856795), R(0.8865333),
      R(0.8897996), R(0.8946404), R(0.9005030), R(0.9070138), R(0.9139161),
      R(0.9210315), R(0.9282347), R(0.9354376), R(0.9425780), R(0.9496124),
      R(0.9565111), R(0.9619183), R(0.9619183), R(0.9619183), R(0.9619183),
      R(0.9619183)};

  constexpr real beta_sh = R(2.2), g_beta_sh = R(0.8974707);
  constexpr real beta_mid = R(1.3), g_beta_mid = R(0.8974707);
  constexpr real beta_dd = R(4.0), g_beta_dd = R(6.0);

  for (int k = 1; k <= kte; ++k) zu[k] = R(0.0);
  int kb_adj = max(kb, 2);

  auto interp_g_alpha = [&](real alpha2) -> real {
    int k = 27;
    for (; k >= 3; --k) {
      if (alpha[k] >= alpha2) break;
    }
    int k1 = k + 1;
    if (alpha[k1] != alpha[k1 - 1]) {
      real a = alpha[k1] - alpha[k1 - 1];
      real b = alpha[k1 - 1] * real(k1) - real(k1 - 1) * alpha[k1];
      real x1 = (alpha2 - b) / a;
      real g_a = g_alpha[k1] - g_alpha[k1 - 1];
      real g_b = g_alpha[k1 - 1] * real(k1) - real(k1 - 1) * g_alpha[k1];
      return g_a * x1 + g_b;
    }
    return g_alpha[k1];
  };

  if (draft == 1) {
    kb_adj = max(kb, 2);
    real trash = -p[kt] + p[kb_adj];
    real tunning = p[kklev];
    if (rand_vmas != R(0.0)) tunning = p[kklev - 1] + R(0.1) * rand_vmas * trash;
    real beta_deep = R(1.3) + (R(1.0) - trash / R(1200.0));
    tunning = min(R(0.95), (tunning - p[kb_adj]) / (p[kt] - p[kb_adj]));
    tunning = max(R(0.02), tunning);
    real alpha2 = (tunning * (beta_deep - R(2.0)) + R(1.0)) / (R(1.0) - tunning);
    (void)interp_g_alpha(alpha2);
    real fzu = tgamma_fn(alpha2 + beta_deep) /
               (tgamma_fn(alpha2) * tgamma_fn(beta_deep));
    zu[kb_adj] = zubeg;
    zu_profile(zu, p, kb_adj + 1, min(kte, kt - 1), p[kb_adj],
               p[kt] - p[kb_adj], zubeg, fzu, alpha2 - R(1.0),
               beta_deep - R(1.0));
    if (zu[kpbli] > R(0.0)) {
      real norm = zu[kpbli];
      int kend = min(ktf, kt - 1);
      for (int k = kts; k <= kend; ++k) zu[k] = zu[k] / norm;
    }
    for (int k = maxloc_abs(zu, 1, kte); k >= 1; --k) {
      if (zu[k] < R(1.e-6)) {
        kb_adj = k + 1;
        break;
      }
    }
    kb_adj = max(2, kb_adj);
    for (int k = kts; k <= kb_adj - 1; ++k) zu[k] = R(0.0);
    constexpr real maxlim = R(1.2);
    real a = maxval_range(zu, 1, kte) - zu[kb_adj];
    for (int k = kb_adj; k <= kt; ++k) {
      if (a > maxlim) {
        zu[k] = (zu[k] - zu[kb_adj]) * maxlim / a + zu[kb_adj];
      }
    }
  } else if (draft == 2) {
    real tunning = p[kklev];
    tunning = min(R(0.95), (tunning - p[kb_adj]) / (p[kt] - p[kb_adj]));
    tunning = max(R(0.02), tunning);
    real alpha2 = (tunning * (beta_sh - R(2.0)) + R(1.0)) / (R(1.0) - tunning);
    real g_alpha2 = interp_g_alpha(alpha2);
    real fzu = tgamma_fn(alpha2 + beta_sh) / (g_alpha2 * g_beta_sh);
    zu[kb_adj] = zubeg;
    zu_profile(zu, p, kb_adj + 1, min(kte, kt - 1), p[kb_adj],
               p[kt] - p[kb_adj], zubeg, fzu, alpha2 - R(1.0),
               beta_sh - R(1.0));
    if (zu[kpbli] > R(0.0)) {
      real norm = zu[kpbli];
      int kend = min(ktf, kt - 1);
      for (int k = kts; k <= kend; ++k) zu[k] = zu[k] / norm;
    }
    for (int k = maxloc_abs(zu, 1, kte); k >= 1; --k) {
      if (zu[k] < R(1.e-6)) {
        kb_adj = k + 1;
        break;
      }
    }
    constexpr real maxlim = R(1.0);
    real a = maxval_range(zu, 1, kte) - zu[kb_adj];
    for (int k = kts; k <= kt; ++k) {
      if (a > maxlim) {
        zu[k] = (zu[k] - zu[kb_adj]) * maxlim / a + zu[kb_adj];
      }
    }
  } else if (draft == 3) {
    kb_adj = max(kb, 2);
    real tunning = R(0.5) * (p[kt] + p[kpbli]);
    tunning = min(R(0.95), (tunning - p[kb_adj]) / (p[kt] - p[kb_adj]));
    tunning = max(R(0.02), tunning);
    real alpha2 = (tunning * (beta_mid - R(2.0)) + R(1.0)) / (R(1.0) - tunning);
    (void)interp_g_alpha(alpha2);
    (void)g_beta_mid;
    real fzu = tgamma_fn(alpha2 + beta_mid) /
               (tgamma_fn(alpha2) * tgamma_fn(beta_mid));
    zu[kb_adj] = zubeg;
    zu_profile(zu, p, kb_adj + 1, min(kte, kt - 1), p[kb_adj],
               p[kt] - p[kb_adj], zubeg, fzu, alpha2 - R(1.0),
               beta_mid - R(1.0));
    if (zu[kpbli] > R(0.0)) {
      real norm = zu[kpbli];
      int kend = min(ktf, kt - 1);
      for (int k = kts; k <= kend; ++k) zu[k] = zu[k] / norm;
    }
    for (int k = maxloc_abs(zu, 1, kte); k >= 1; --k) {
      if (zu[k] < R(1.e-6)) {
        kb_adj = k + 1;
        break;
      }
    }
    kb_adj = max(2, kb_adj);
    for (int k = kts; k <= kb_adj - 1; ++k) zu[k] = R(0.0);
    constexpr real maxlim = R(1.5);
    real a = maxval_range(zu, 1, kte) - zu[kb_adj];
    for (int k = kts; k <= kt; ++k) {
      if (a > maxlim) {
        zu[k] = (zu[k] - zu[kb_adj]) * maxlim / a + zu[kb_adj];
      }
    }
  } else if (draft == 4 || draft == 5) {
    real tunning = p[kb];
    tunning = min(R(0.95), (tunning - p[1]) / (p[kt] - p[1]));
    tunning = max(R(0.02), tunning);
    real alpha2 = (tunning * (beta_dd - R(2.0)) + R(1.0)) / (R(1.0) - tunning);
    real g_alpha2 = interp_g_alpha(alpha2);
    real fzu = tgamma_fn(alpha2 + beta_dd) / (g_alpha2 * g_beta_dd);
    for (int k = 1; k <= kte; ++k) zu[k] = R(0.0);
    zu_profile(zu, p, 2, min(kte, kt - 1), p[1], p[kt] - p[1], R(0.0), fzu,
               alpha2 - R(1.0), beta_dd - R(1.0));
    int kend = min(ktf, kt - 1);
    fzu = maxval_range(zu, kts, kend);
    if (fzu > R(0.0)) {
      for (int k = kts; k <= kend; ++k) zu[k] = zu[k] / fzu;
    }
    zu[1] = R(0.0);
    for (int k = 1; k <= kb - 2; ++k) {
      zu[kb - k] =
          zu[kb - k + 1] - zu[kb] * (p[kb - k] - p[kb - k + 1]) / (p[1] - p[kb]);
    }
    zu[1] = R(0.0);
  }
  (void)ipr;
  (void)xland;
  (void)zuh2;
  (void)ierr;
  (void)max_mass;
  (void)csum;
  (void)pmin_lev;
}

// ---------------------------------------------------------------------------
// rates_up_pdf (lines 5016-5151 of module_cu_gfl_deep.F90)
// mode: 1 = deep, 2 = shallow, 3 = mid
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void rates_up_pdf(
    int mode, int ktf, int kte, int& ierr, const real* p_cup, int& ktop,
    int& k22, const real* z_cup, const real* entr_rate_2d, const real* heo,
    const real* heso_cup, real hkbo, int xland, int& kbcon, real* zuo,
    int& ktopdby, int csum, int pmin_lev) {
  if (ierr > 0) return;
  for (int k = 1; k <= kte; ++k) zuo[k] = R(0.0);

  constexpr real zustart = R(0.1);
  const real dbythresh = (mode == 2 || mode == 3) ? R(1.0) : R(0.8);
  const real beta_u = max(R(0.1), R(0.2) - real(csum) * R(0.01));

  real dby[KM] = {};
  real dbm[KM] = {};
  real hcot[KM] = {};
  real zu_local[KM] = {};
  real zuh2[40] = {};

  kbcon = max(kbcon, 2);
  int start_level = max(1, k22);
  zuo[start_level] = zustart;
  real entr_init = entr_rate_2d[1];
  for (int k = start_level + 1; k <= kbcon; ++k) {
    real dz = z_cup[k] - z_cup[k - 1];
    real massent = dz * entr_rate_2d[k - 1] * zuo[k - 1];
    real massdetr = dz * R(0.1) * entr_init * zuo[k - 1];
    zuo[k] = zuo[k - 1] + massent - massdetr;
  }
  real zubeg = zustart;

  if (mode == 1) {
    ktop = 0;
    hcot[start_level] = hkbo;
    for (int k = start_level + 1; k <= ktf - 2; ++k) {
      real dz = z_cup[k] - z_cup[k - 1];
      hcot[k] = ((R(1.0) - R(0.5) * entr_rate_2d[k - 1] * dz) * hcot[k - 1] +
                 entr_rate_2d[k - 1] * dz * heo[k - 1]) /
                (R(1.0) + R(0.5) * entr_rate_2d[k - 1] * dz);
      if (k >= kbcon) dby[k] = dby[k - 1] + (hcot[k] - heso_cup[k]) * dz;
      if (k >= kbcon) dbm[k] = hcot[k] - heso_cup[k];
    }
    int max_dby_k = 1;
    real max_dby = dby[1];
    int kklev = 1;
    real max_dbm = dbm[1];
    for (int k = 2; k <= kte; ++k) {
      if (dby[k] > max_dby) {
        max_dby = dby[k];
        max_dby_k = k;
      }
      if (dbm[k] > max_dbm) {
        max_dbm = dbm[k];
        kklev = k;
      }
    }
    ktopdby = max_dby_k;
    int kfinalzu = ktf - 2;
    ktop = kfinalzu;
    for (int k = max_dby_k + 1; k <= ktf - 2; ++k) {
      if (dby[k] < dbythresh * max_dby) {
        kfinalzu = k - 1;
        ktop = kfinalzu;
        break;
      }
    }
    ktop = ktopdby;
    kklev = min(kklev + 3, ktop - 2);
    if (kfinalzu <= kbcon + 2) {
      ierr = 41;
      ktop = 0;
    } else {
      for (int k = 1; k <= kte; ++k) zu_local[k] = zuo[k];
      get_zu_zd_pdf_fim(kklev, p_cup, R(0.0), zubeg, 0, xland, zuh2, 1, ierr,
                        k22, kfinalzu + 1, zu_local, 1, kte, ktf, beta_u, kbcon,
                        csum, pmin_lev);
      for (int k = 1; k <= ktf; ++k) zuo[k] = zu_local[k];
    }
  } else if (mode == 3) {
    if (ktop <= kbcon + 2) {
      ierr = 41;
      ktop = 0;
    } else {
      ktopdby = ktop + 1;
      int kklev = 1;
      for (int k = 1; k <= kte; ++k) zu_local[k] = zuo[k];
      get_zu_zd_pdf_fim(kklev, p_cup, R(0.0), zubeg, 0, xland, zuh2, 3, ierr,
                        k22, ktopdby + 1, zu_local, 1, kte, ktf, beta_u, kbcon,
                        csum, pmin_lev);
      for (int k = 1; k <= ktf; ++k) zuo[k] = zu_local[k];
    }
  } else if (mode == 2) {
    if (ktop <= kbcon + 2) {
      ierr = 41;
      ktop = 0;
    } else {
      ktopdby = ktop + 1;
      for (int k = 1; k <= kte; ++k) zu_local[k] = zuo[k];
      get_zu_zd_pdf_fim(kbcon, p_cup, R(0.0), zubeg, 0, xland, zuh2, 2, ierr,
                        k22, ktopdby + 1, zu_local, 1, kte, ktf, beta_u, kbcon,
                        csum, pmin_lev);
      for (int k = 1; k <= ktf; ++k) zuo[k] = zu_local[k];
    }
  }
}

// ---------------------------------------------------------------------------
// get_lateral_massflux (lines 5684-5812 of module_cu_gfl_deep.F90)
// draft: 1 = deep, 2 = shallow, 3 = mid
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void get_lateral_massflux(
    int ktf, int kte, int ierr, int ktop, const real* zo_cup, const real* zuo,
    real* cd, real* entr_rate_2d, real* up_massentro, real* up_massdetro,
    real* up_massentr, real* up_massdetr, int draft, int kbcon, int k22,
    real* up_massentru, real* up_massdetru, real& lambau) {
  for (int k = 1; k <= kte; ++k) {
    up_massentro[k] = R(0.0);
    up_massdetro[k] = R(0.0);
    up_massentr[k] = R(0.0);
    up_massdetr[k] = R(0.0);
    if (up_massentru && up_massdetru) {
      up_massentru[k] = R(0.0);
      up_massdetru[k] = R(0.0);
    }
  }
  if (ierr != 0) return;

  int kmx = maxloc_abs(zuo, 1, kte);
  for (int k = max(2, k22 + 1); k <= kmx; ++k) {
    real dz = zo_cup[k] - zo_cup[k - 1];
    up_massdetro[k - 1] = cd[k - 1] * dz * zuo[k - 1];
    up_massentro[k - 1] = zuo[k] - zuo[k - 1] + up_massdetro[k - 1];
    if (up_massentro[k - 1] < R(0.0)) {
      up_massentro[k - 1] = R(0.0);
      up_massdetro[k - 1] = zuo[k - 1] - zuo[k];
      if (zuo[k - 1] > R(0.0)) {
        cd[k - 1] = up_massdetro[k - 1] / (dz * zuo[k - 1]);
      }
    }
    if (zuo[k - 1] > R(0.0)) {
      entr_rate_2d[k - 1] = up_massentro[k - 1] / (dz * zuo[k - 1]);
    }
  }

  for (int k = kmx + 1; k <= ktop; ++k) {
    real dz = zo_cup[k] - zo_cup[k - 1];
    up_massentro[k - 1] = entr_rate_2d[k - 1] * dz * zuo[k - 1];
    up_massdetro[k - 1] = zuo[k - 1] + up_massentro[k - 1] - zuo[k];
    if (up_massdetro[k - 1] < R(0.0)) {
      up_massdetro[k - 1] = R(0.0);
      up_massentro[k - 1] = zuo[k] - zuo[k - 1];
      if (zuo[k - 1] > R(0.0)) {
        entr_rate_2d[k - 1] = up_massentro[k - 1] / (dz * zuo[k - 1]);
      }
    }
    if (zuo[k - 1] > R(0.0)) {
      cd[k - 1] = up_massdetro[k - 1] / (dz * zuo[k - 1]);
    }
  }

  up_massdetro[ktop] = zuo[ktop];
  up_massentro[ktop] = R(0.0);
  for (int k = ktop + 1; k <= ktf; ++k) {
    cd[k] = R(0.0);
    entr_rate_2d[k] = R(0.0);
    up_massentro[k] = R(0.0);
    up_massdetro[k] = R(0.0);
  }
  for (int k = 2; k <= ktf - 1; ++k) {
    up_massentr[k - 1] = up_massentro[k - 1];
    up_massdetr[k - 1] = up_massdetro[k - 1];
  }

  if (up_massentru && up_massdetru) {
    if (draft == 3) lambau = R(0.0);
    if (draft == 1 || draft == 2 || draft == 3) {
      for (int k = 2; k <= ktf - 1; ++k) {
        up_massentru[k - 1] =
            up_massentro[k - 1] + lambau * up_massdetro[k - 1];
        up_massdetru[k - 1] =
            up_massdetro[k - 1] + lambau * up_massdetro[k - 1];
      }
    }
  }
  (void)kbcon;
}

// ---------------------------------------------------------------------------
// cup_up_moisture (lines 4585-4959 of module_cu_gfl_deep.F90)
// ---------------------------------------------------------------------------
template <typename real>
CPPWRF_HOST_DEVICE inline void cup_up_moisture(
    bool is_mid, int& ierr, const real* z_cup, real* qc, real* qrc, real* pw,
    real& pwav, const real* p_cup, int kbcon, int ktop, const real* dby,
    real* clw_all, int xland1, const real* q, const real* gamma_cup,
    const real* zu, const real* qes_cup, int k22, const real* qe_cup, real c0,
    real* c0t3d, real zqexec, real ccn, real ccnclean, const real* rho,
    real* c1d, const real* t, int autoconv, const real* up_massentr,
    const real* up_massdetr, real& psum, real& psumh, int itest, int ktf,
    int kte) {
  using C = GflConstants<real>;

  constexpr real clwdet = R(0.1);
  constexpr real c0_iceconv = R(0.01);

  for (int k = 1; k <= kte; ++k) c0t3d[k] = R(0.0);
  pwav = R(0.0);
  psum = R(0.0);
  psumh = R(0.0);

  for (int k = 1; k <= ktf; ++k) {
    pw[k] = R(0.0);
    qc[k] = (ierr == 0) ? qe_cup[k] : R(0.0);
    clw_all[k] = R(0.0);
    qrc[k] = R(0.0);
  }
  if (ierr != 0) return;

  int start_level = k22;
  real qaver = R(0.0);
  get_cloud_bc(kte, qe_cup, qaver, k22, R(0.0));
  qc[start_level] = qaver;
  for (int k = 1; k <= start_level - 1; ++k) {
    qc[k] = qe_cup[k];
  }

  for (int k = k22 + 1; k <= kbcon; ++k) {
    real c0t =
        (t[k] > R(273.16)) ? c0 : c0 * exp(c0_iceconv * (t[k] - R(273.16)));
    c0t3d[k] = c0t;
    qc[k] = (qc[k - 1] * zu[k - 1] - R(0.5) * up_massdetr[k - 1] * qc[k - 1] +
             up_massentr[k - 1] * q[k - 1]) /
            (zu[k - 1] - R(0.5) * up_massdetr[k - 1] + up_massentr[k - 1]);
    real qrch = qes_cup[k] + (R(1.0) / C::xlv) *
                                 (gamma_cup[k] / (R(1.0) + gamma_cup[k])) *
                                 dby[k];
    if (k < kbcon) qrch = qc[k];
    if (qc[k] > qrch) {
      real dz = z_cup[k] - z_cup[k - 1];
      qrc[k] = (qc[k] - qrch) / (R(1.0) + c0t * dz);
      pw[k] = c0t * dz * qrc[k] * zu[k];
      qc[k] = qrch + qrc[k];
      clw_all[k] = qrc[k];
    }
  }

  for (int k = kbcon + 1; k <= ktop; ++k) {
    real c0t =
        (t[k] > R(273.16)) ? c0 : c0 * exp(c0_iceconv * (t[k] - R(273.16)));
    if (is_mid) c0t = R(0.004);
    c0t3d[k] = c0t;
    if (autoconv > 1) c0t = c0;

    real denom = zu[k - 1] - R(0.5) * up_massdetr[k - 1] + up_massentr[k - 1];
    if (denom < R(1.e-16)) {
      ierr = 51;
      break;
    }

    real dz = z_cup[k] - z_cup[k - 1];
    real dp = -R(100.0) * (p_cup[k] - p_cup[k - 1]);
    real qrch = qes_cup[k] + (R(1.0) / C::xlv) *
                                 (gamma_cup[k] / (R(1.0) + gamma_cup[k])) *
                                 dby[k];
    qc[k] = (qc[k - 1] * zu[k - 1] - R(0.5) * up_massdetr[k - 1] * qc[k - 1] +
             up_massentr[k - 1] * q[k - 1]) /
            (zu[k - 1] - R(0.5) * up_massdetr[k - 1] + up_massentr[k - 1]);
    if (qc[k] <= qrch) qc[k] = qrch + R(1.e-8);

    clw_all[k] = max(R(0.0), qc[k] - qrch);
    qrc[k] = max(R(0.0), qc[k] - qrch);
    if (k > kbcon + 1) c1d[k] = clwdet * up_massdetr[k - 1];

    qrc[k] = (qc[k] - qrch) / (R(1.0) + (c1d[k] + c0t) * dz);
    if (qrc[k] < R(0.0)) qrc[k] = R(0.0);
    pw[k] = c0t * dz * qrc[k] * zu[k];
    if (qrc[k] < R(0.0)) {
      qrc[k] = R(0.0);
      pw[k] = R(0.0);
    }
    qc[k] = qrc[k] + qrch;
    pwav += pw[k];
    psum += pw[k] * C::g / dp;
  }

  for (int k = k22 + 1; k <= ktop; ++k) {
    qc[k] = qc[k] - qrc[k];
  }
  (void)xland1;
  (void)zqexec;
  (void)ccn;
  (void)ccnclean;
  (void)rho;
  (void)itest;
}

#undef R

}  // namespace cppwrf::phys::cu_gfl
