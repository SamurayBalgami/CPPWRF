#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>
#include "core/backend.hpp"
#include "core/constants.hpp"
#include "core/field.hpp"
#include "dyn/diffusion/diffusion_params.hpp"

namespace cppwrf::dyn::diffusion {

// Ported from WRF v4.8.0 dyn_em/module_diffusion_em.F:
//   calculate_km_kh (lines 1195-1349)
//   cal_dampkm      (lines 1353-1480, damp_opt = 1)
//   calculate_N2    (lines 1485-1713)
//   isotropic_km    (lines 1718-1772, km_opt = 1)
//   smag_km         (lines 1777-1929, km_opt = 3)
//   smag2d_km       (lines 1934-2044, km_opt = 4)
//   tke_km          (lines 2049-2336, km_opt = 2)
//   calc_l_scale    (lines 2341-2406)

CPPWRF_HD inline real vcos_diff(real x, int lane_i, int ni) {
#if defined(CPPWRF_HAVE_LIBMVEC)
  if constexpr (std::is_same_v<real, double>) {
    return cppwrf::vm::cos(x, lane_i, ni);
  }
#endif
  (void)lane_i;
  (void)ni;
  return std::cos(x);
}

struct KIsotropicKm {
  View3<real> xkmh, xkmv, xkhh, xkhv;
  real khdif, kvdif, prandtl;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    xkmh(i, j, k) = khdif;
    xkmv(i, j, k) = kvdif;
    xkhh(i, j, k) = khdif / prandtl;
    xkhv(i, j, k) = kvdif / prandtl;
  }
};

struct KSmag2dKm {
  View3<const real> defor11, defor22, defor12, rdzw, zx, zy;
  BView2<const real> msftx, msfty;
  View3<real> xkmh, xkmv, xkhh, xkhv;
  real dx, dy, c_s, prandtl;
  int diff_opt, nz;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    // def2 = 0.25 * (defor11 - defor22)^2 + (0.25 * sum(defor12))^2
    real d11_m_d22 = defor11(i, j, k) - defor22(i, j, k);
    real def2 = real(0.25) * (d11_m_d22 * d11_m_d22);
    real d12_avg = real(0.25) * (defor12(i, j, k) + defor12(i, j + 1, k) +
                                 defor12(i + 1, j, k) + defor12(i + 1, j + 1, k));
    def2 += d12_avg * d12_avg;

    real mlen_h = std::sqrt(dx / msftx(i, j) * dy / msfty(i, j));  // WRF: ((dx/msftx)*dy)/msfty
    real tmp = std::sqrt(def2);
    real kmh = c_s * c_s * mlen_h * mlen_h * tmp;
    kmh = kmh < real(10.0) * mlen_h ? kmh : real(10.0) * mlen_h;

    real kmv = kmh;
    real khh = kmh / prandtl;
    real khv = real(0.0);

    if (diff_opt == 2) {
      real dxm = dx / msftx(i, j);
      real dym = dy / msfty(i, j);
      real tmpzx = real(0.25) * (std::abs(zx(i, j, k)) + std::abs(zx(i + 1, j, k)) +
                                 std::abs(zx(i, j, k + 1)) + std::abs(zx(i + 1, j, k + 1))) *
                   rdzw(i, j, k) * dxm;
      real tmpzy = real(0.25) * (std::abs(zy(i, j, k)) + std::abs(zy(i, j + 1, k)) +
                                 std::abs(zy(i, j, k + 1)) + std::abs(zy(i, j + 1, k + 1))) *
                   rdzw(i, j, k) * dym;
      real slope = std::sqrt(tmpzx * tmpzx + tmpzy * tmpzy);
      real alpha = slope > real(1.0) ? slope : real(1.0);
      real def_limit = (real(10.0) / mlen_h) > real(1e-3) ? (real(10.0) / mlen_h) : real(1e-3);

      if (tmp > def_limit) {
        kmh = kmh / (alpha * alpha);
      } else {
        kmh = kmh / alpha;
      }
      khh = kmh / prandtl;
      kmv = kmh;
    }

    xkmh(i, j, k) = kmh;
    xkmv(i, j, k) = kmv;
    xkhh(i, j, k) = khh;
    xkhv(i, j, k) = khv;
  }
};

// -----------------------------------------------------------------------------
// calculate_N2 (lines 1485-1713)
// -----------------------------------------------------------------------------
struct KCalculateN2 {
  View3<const real> theta, t_phy, p_phy, p8w, t8w, rdz, rdzw;
  View3<const real> qv, qc, qi;
  bool has_qv, has_qc, has_qi;
  real cf1, cf2, cf3;
  int i_start, ni, nz;
  View3<real> bn2;

  CPPWRF_HD void operator()(int i, int j) const {
    constexpr real g       = real(9.81);
    constexpr real r_d     = real(287.0);
    constexpr real cp      = real(7.0) * r_d / real(2.0);
    constexpr real r_v     = real(461.6);
    constexpr real p1000mb = real(100000.0);
    constexpr real xlv     = real(2.5e6);
    constexpr real ep_2    = r_d / r_v;
    constexpr real svp1    = real(0.6112);
    constexpr real svp2    = real(17.67);
    constexpr real svp3    = real(29.65);
    constexpr real svpt0   = real(273.15);
    constexpr real qc_cr   = real(0.00001);

    const int lane_i = i - i_start;
    real tmp1[128], qvs[128], qctmp[128];
    for (int k = 0; k < nz; ++k) {
      qctmp[k] = has_qc ? qc(i, j, k) : real(0.0);
      tmp1[k]  = real(0.0);
    }
    real tmp1sfc = real(0.0);

    if (has_qv) {
      for (int k = 0; k < nz; ++k) tmp1[k] = tmp1[k] + qv(i, j, k);
      tmp1sfc = tmp1sfc + cf1 * qv(i, j, 0) + cf2 * qv(i, j, 1) + cf3 * qv(i, j, 2);
    }
    if (has_qc) {
      for (int k = 0; k < nz; ++k) tmp1[k] = tmp1[k] + qc(i, j, k);
      tmp1sfc = tmp1sfc + cf1 * qc(i, j, 0) + cf2 * qc(i, j, 1) + cf3 * qc(i, j, 2);
    }
    if (has_qi) {
      for (int k = 0; k < nz; ++k) tmp1[k] = tmp1[k] + qi(i, j, k);
      tmp1sfc = tmp1sfc + cf1 * qi(i, j, 0) + cf2 * qi(i, j, 1) + cf3 * qi(i, j, 2);
    }

    for (int k = 0; k < nz; ++k) {
      real tc = t_phy(i, j, k) - svpt0;
      real es = real(1000.0) * svp1 * cppwrf::vm::exp(svp2 * tc / (t_phy(i, j, k) - svp3), lane_i, ni);
      qvs[k]  = ep_2 * es / (p_phy(i, j, k) - es);
    }

    for (int k = 1; k < nz - 1; ++k) {
      real tmpdz = real(1.0) / rdz(i, j, k) + real(1.0) / rdz(i, j, k + 1);
      real qv_k  = has_qv ? qv(i, j, k) : real(0.0);
      if (qv_k >= qvs[k] || qctmp[k] >= qc_cr) {
        real xlvqv    = xlv * qv_k;
        real coefa    = (real(1.0) + xlvqv / r_d / t_phy(i, j, k)) /
                        (real(1.0) + xlv * xlvqv / cp / r_v / t_phy(i, j, k) / t_phy(i, j, k)) /
                        theta(i, j, k);
        real thetaep1 = theta(i, j, k + 1) *
                        (real(1.0) + xlv * qvs[k + 1] / cp / t_phy(i, j, k + 1));
        real thetaem1 = theta(i, j, k - 1) *
                        (real(1.0) + xlv * qvs[k - 1] / cp / t_phy(i, j, k - 1));
        bn2(i, j, k)  = g * (coefa * (thetaep1 - thetaem1) / tmpdz -
                             (tmp1[k + 1] - tmp1[k - 1]) / tmpdz);
      } else {
        real qv_kp1  = has_qv ? qv(i, j, k + 1) : real(0.0);
        real qv_km1  = has_qv ? qv(i, j, k - 1) : real(0.0);
        bn2(i, j, k) = g * ((theta(i, j, k + 1) - theta(i, j, k - 1)) /
                            theta(i, j, k) / tmpdz +
                            real(1.61) * (qv_kp1 - qv_km1) / tmpdz -
                            (tmp1[k + 1] - tmp1[k - 1]) / tmpdz);
      }
    }

    {
      const int k   = 0;
      real tmpdz    = real(1.0) / rdz(i, j, k + 1) + real(0.5) / rdzw(i, j, k);
      real thetasfc = t8w(i, j, 0) / rpow(p8w(i, j, k) / p1000mb, r_d / cp);
      real qv_0     = has_qv ? qv(i, j, 0) : real(0.0);
      if (qv_0 >= qvs[0] || qctmp[0] >= qc_cr) {
        real qvsfc    = cf1 * qvs[0] + cf2 * qvs[1] + cf3 * qvs[2];
        real xlvqv    = xlv * qv_0;
        real coefa    = (real(1.0) + xlvqv / r_d / t_phy(i, j, 0)) /
                        (real(1.0) + xlv * xlvqv / cp / r_v / t_phy(i, j, 0) / t_phy(i, j, 0)) /
                        theta(i, j, 0);
        real thetaep1  = theta(i, j, 1) *
                         (real(1.0) + xlv * qvs[1] / cp / t_phy(i, j, 1));
        real thetaesfc = thetasfc *
                         (real(1.0) + xlv * qvsfc / cp / t8w(i, j, 0));
        bn2(i, j, 0)   = g * (coefa * (thetaep1 - thetaesfc) / tmpdz -
                              (tmp1[1] - tmp1sfc) / tmpdz);
      } else {
        real qv_1    = has_qv ? qv(i, j, 1) : real(0.0);
        real qv_2    = has_qv ? qv(i, j, 2) : real(0.0);
        real qvsfc   = cf1 * qv_0 + cf2 * qv_1 + cf3 * qv_2;
        tmpdz        = real(1.0) / rdzw(i, j, 0);
        bn2(i, j, 0) = g * ((theta(i, j, 1) - theta(i, j, 0)) /
                            theta(i, j, 0) / tmpdz +
                            real(1.61) * (qv_1 - qvsfc) / tmpdz -
                            (tmp1[1] - tmp1sfc) / tmpdz);
      }
      bn2(i, j, nz - 1) = bn2(i, j, nz - 2);
    }
  }
};

inline void calculate_N2(const DiffusionParams& p,
                         View3<const real> theta,
                         View3<const real> t_phy,
                         View3<const real> p_phy,
                         View3<const real> p8w,
                         View3<const real> t8w,
                         View3<const real> rdz,
                         View3<const real> rdzw,
                         View3<const real> qv, bool has_qv,
                         View3<const real> qc, bool has_qc,
                         View3<const real> qi, bool has_qi,
                         real cf1, real cf2, real cf3,
                         View3<real> bn2,
                         int nx, int ny, int nz) {
  int i_start = (p.open_xs || p.specified || p.nested) ? 1 : 0;
  int i_end   = (p.open_xe || p.specified || p.nested) ? nx - 1 : nx;
  int j_start = (p.open_ys || p.specified || p.nested) ? 1 : 0;
  int j_end   = (p.open_ye || p.specified || p.nested) ? ny - 1 : ny;
  if (p.periodic_x) {
    i_start = 0;
    i_end   = nx;
  }
  const int ni = i_end - i_start;
  column_for(Range2D{i_start, i_end, j_start, j_end},
             KCalculateN2{theta, t_phy, p_phy, p8w, t8w, rdz, rdzw,
                          qv, qc, qi, has_qv, has_qc, has_qi,
                          cf1, cf2, cf3, i_start, ni, nz, bn2});
}

// -----------------------------------------------------------------------------
// calc_l_scale (lines 2341-2406)
// -----------------------------------------------------------------------------
CPPWRF_HD inline void calc_l_scale_col(int i, int j, int nz,
                                       real dx, real dy,
                                       BView2<const real> msftx,
                                       BView2<const real> msfty,
                                       View3<const real> tke,
                                       View3<const real> bn2,
                                       View3<const real> rdzw,
                                       real* l_scale) {
  for (int k = 0; k < nz; ++k) {
    real deltas = rpow(dx / msftx(i, j) * dy / msfty(i, j) / rdzw(i, j, k), real(0.33333333));
    l_scale[k]  = deltas;
    if (bn2(i, j, k) > real(1.0e-6)) {
      real tketmp = tke(i, j, k) > real(1.0e-6) ? tke(i, j, k) : real(1.0e-6);
      real tmp    = std::sqrt(tketmp);
      l_scale[k]  = real(0.76) * tmp / std::sqrt(bn2(i, j, k));
      if (l_scale[k] > deltas) l_scale[k] = deltas;
      real min_l = real(0.001) * deltas;
      if (l_scale[k] < min_l) l_scale[k] = min_l;
    }
  }
}

// -----------------------------------------------------------------------------
// tke_km (lines 2049-2336, km_opt = 2)
// -----------------------------------------------------------------------------
struct KTkeKm {
  View3<const real> bn2, tke, p8w, t8w, theta, rdz, rdzw;
  BView2<const real> msftx, msfty;
  View3<real> xkmh, xkmv, xkhh, xkhv;
  real dx, dy, dt, c_k, prandtl, mix_upper_bound, tke_seed;
  int isotropic, i_start, ni, nz;

  CPPWRF_HD void operator()(int i, int j) const {
    constexpr real g       = real(9.81);
    constexpr real r_d     = real(287.0);
    constexpr real cp      = real(7.0) * r_d / real(2.0);
    constexpr real p1000mb = real(100000.0);

    const int lane_i = i - i_start;
    real dthrdn[128];

    for (int k = 1; k < nz - 1; ++k) {
      real tmpdz = real(1.0) / rdz(i, j, k + 1) + real(1.0) / rdz(i, j, k);
      dthrdn[k]  = (theta(i, j, k + 1) - theta(i, j, k - 1)) / tmpdz;
    }
    {
      const int k   = 0;
      real tmpdz    = real(1.0) / rdzw(i, j, k + 1) + real(1.0) / rdzw(i, j, k);
      real thetasfc = t8w(i, j, 0) / cppwrf::vm::pow(p8w(i, j, k) / p1000mb, r_d / cp, lane_i, ni);
      dthrdn[0]     = (theta(i, j, 1) - thetasfc) / tmpdz;
    }
    {
      const int k   = nz - 1;
      real tmpdz    = real(1.0) / rdz(i, j, k) + real(0.5) / rdzw(i, j, k);
      real thetatop = t8w(i, j, nz) / cppwrf::vm::pow(p8w(i, j, nz) / p1000mb, r_d / cp, lane_i, ni);
      dthrdn[k]     = (thetatop - theta(i, j, k - 1)) / tmpdz;
    }

    if (isotropic == 0) {
      for (int k = 0; k < nz; ++k) {
        real mlen_h  = std::sqrt(dx / msftx(i, j) * dy / msfty(i, j));
        real tke_val = tke(i, j, k) > tke_seed ? tke(i, j, k) : tke_seed;
        real tmp     = std::sqrt(tke_val);
        real deltas  = real(1.0) / rdzw(i, j, k);
        real mlen_v  = deltas;
        if (dthrdn[k] > real(0.0)) {
          real mlen_s = real(0.76) * tmp / rpow(std::abs(g / theta(i, j, k) * dthrdn[k]), real(0.5));
          if (mlen_s < mlen_v) mlen_v = mlen_s;
        }
        real kmh     = c_k * tmp * mlen_h;
        real kmh_min = real(1.0e-6) * mlen_h * mlen_h;
        if (kmh < kmh_min) kmh = kmh_min;
        real kmh_max = mix_upper_bound * mlen_h * mlen_h / dt;
        if (kmh > kmh_max) kmh = kmh_max;

        real kmv     = c_k * tmp * mlen_v;
        real kmv_min = real(1.0e-6) * deltas * deltas;
        if (kmv < kmv_min) kmv = kmv_min;
        real kmv_max = mix_upper_bound * deltas * deltas / dt;
        if (kmv > kmv_max) kmv = kmv_max;

        real pr_inv_h = real(1.0) / prandtl;
        real pr_inv_v = real(1.0) + real(2.0) * mlen_v / deltas;
        xkmh(i, j, k) = kmh;
        xkmv(i, j, k) = kmv;
        xkhh(i, j, k) = kmh * pr_inv_h;
        xkhv(i, j, k) = kmv * pr_inv_v;
      }
    } else {
      real l_scale[128];
      calc_l_scale_col(i, j, nz, dx, dy, msftx, msfty, tke, bn2, rdzw, l_scale);
      for (int k = 0; k < nz; ++k) {
        real tke_val = tke(i, j, k) > tke_seed ? tke(i, j, k) : tke_seed;
        real tmp     = std::sqrt(tke_val);
        real deltas  = cppwrf::vm::pow(dx / msftx(i, j) * dy / msfty(i, j) / rdzw(i, j, k),
                                       real(0.33333333), lane_i, ni);
        real kmh       = c_k * tmp * l_scale[k];
        real kmh_bound = mix_upper_bound * dx / msftx(i, j) * dy / msfty(i, j) / dt;
        if (kmh_bound < kmh) kmh = kmh_bound;

        real kmv       = c_k * tmp * l_scale[k];
        real kmv_bound = mix_upper_bound / rdzw(i, j, k) / rdzw(i, j, k) / dt;
        if (kmv_bound < kmv) kmv = kmv_bound;

        real pr_inv = real(1.0) + real(2.0) * l_scale[k] / deltas;
        real khh    = kmh * pr_inv;
        if (kmh_bound < khh) khh = kmh_bound;
        real khv    = kmv * pr_inv;
        if (kmv_bound < khv) khv = kmv_bound;

        xkmh(i, j, k) = kmh;
        xkmv(i, j, k) = kmv;
        xkhh(i, j, k) = khh;
        xkhv(i, j, k) = khv;
      }
    }
  }
};

inline void tke_km(const DiffusionParams& p,
                   real dx, real dy, real dt,
                   View3<const real> bn2,
                   View3<const real> tke,
                   View3<const real> p8w,
                   View3<const real> t8w,
                   View3<const real> theta,
                   View3<const real> rdz,
                   View3<const real> rdzw,
                   BView2<const real> msftx,
                   BView2<const real> msfty,
                   View3<real> xkmh,
                   View3<real> xkmv,
                   View3<real> xkhh,
                   View3<real> xkhv,
                   int nx, int ny, int nz) {
  int i_start = (p.open_xs || p.specified || p.nested) ? 1 : 0;
  int i_end   = (p.open_xe || p.specified || p.nested) ? nx - 1 : nx;
  int j_start = (p.open_ys || p.specified || p.nested) ? 1 : 0;
  int j_end   = (p.open_ye || p.specified || p.nested) ? ny - 1 : ny;
  if (p.periodic_x) {
    i_start = 0;
    i_end   = nx;
  }
  const int ni = i_end - i_start;

  real tke_seed = real(0.0);
  if (p.isfflx == 0) {
    if (p.diff_opt == 2 && p.bl_pbl_physics == 0) {
      if (p.tke_drag_coefficient < real(1.0e-10) && p.tke_heat_flux < real(1.0e-10)) {
        tke_seed = real(1.0e-6);
      }
    } else {
      tke_seed = real(1.0e-6);
    }
  }

  column_for(Range2D{i_start, i_end, j_start, j_end},
             KTkeKm{bn2, tke, p8w, t8w, theta, rdz, rdzw, msftx, msfty,
                    xkmh, xkmv, xkhh, xkhv, dx, dy, dt, p.c_k, p.prandtl,
                    p.mix_upper_bound, tke_seed, p.isotropic, i_start, ni, nz});
}

// -----------------------------------------------------------------------------
// smag_km (lines 1777-1929, km_opt = 3)
// -----------------------------------------------------------------------------
struct KSmagKm {
  View3<const real> defor11, defor22, defor33, defor12, defor13, defor23, bn2, rdzw;
  BView2<const real> msftx, msfty;
  View3<real> xkmh, xkmv, xkhh, xkhv;
  real dx, dy, dt, c_s, prandtl, mix_upper_bound;
  int isotropic, i_start, ni, nz;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    const int lane_i = i - i_start;
    real def2 = real(0.5) * (defor11(i, j, k) * defor11(i, j, k) +
                             defor22(i, j, k) * defor22(i, j, k) +
                             defor33(i, j, k) * defor33(i, j, k));
    real tmp12 = real(0.25) * (defor12(i, j, k) + defor12(i, j + 1, k) +
                               defor12(i + 1, j, k) + defor12(i + 1, j + 1, k));
    def2 = def2 + tmp12 * tmp12;
    real tmp13 = real(0.25) * (defor13(i, j, k + 1) + defor13(i, j, k) +
                               defor13(i + 1, j, k + 1) + defor13(i + 1, j, k));
    def2 = def2 + tmp13 * tmp13;
    real tmp23 = real(0.25) * (defor23(i, j, k + 1) + defor23(i, j, k) +
                               defor23(i, j + 1, k + 1) + defor23(i, j + 1, k));
    def2 = def2 + tmp23 * tmp23;

    real pr  = prandtl;
    real arg = def2 - bn2(i, j, k) / pr;
    // tmp**0.5 in smag_km's vectorised i loops: gfortran turns it into sqrtpd for the iteration
    // pairs and keeps pow for an odd remainder iteration (module_diffusion_em.o)
    const real tmp0 = arg > real(0.0) ? arg : real(0.0);
    real tmp = cppwrf::vm::lane(lane_i, ni) ? std::sqrt(tmp0) : rpow(tmp0, real(0.5));

    if (isotropic == 0) {
      real mlen_h  = std::sqrt(dx / msftx(i, j) * dy / msfty(i, j));
      real mlen_v  = real(1.0) / rdzw(i, j, k);
      real kmh     = c_s * c_s * mlen_h * mlen_h * tmp;
      real kmh_min = real(1.0e-6) * mlen_h * mlen_h;
      if (kmh < kmh_min) kmh = kmh_min;
      real kmh_max = mix_upper_bound * mlen_h * mlen_h / dt;
      if (kmh > kmh_max) kmh = kmh_max;

      real kmv     = c_s * c_s * mlen_v * mlen_v * tmp;
      real kmv_min = real(1.0e-6) * mlen_v * mlen_v;
      if (kmv < kmv_min) kmv = kmv_min;
      real kmv_max = mix_upper_bound * mlen_v * mlen_v / dt;
      if (kmv > kmv_max) kmv = kmv_max;

      real khh = kmh / pr;
      if (khh > kmh_max) khh = kmh_max;
      real khv = kmv / pr;
      if (khv > kmv_max) khv = kmv_max;

      xkmh(i, j, k) = kmh;
      xkmv(i, j, k) = kmv;
      xkhh(i, j, k) = khh;
      xkhv(i, j, k) = khv;
    } else {
      real deltas  = cppwrf::vm::pow(dx / msftx(i, j) * dy / msfty(i, j) / rdzw(i, j, k),
                                     real(0.33333333), lane_i, ni);
      real kmh     = c_s * c_s * deltas * deltas * tmp;
      real kmh_min = real(1.0e-6) * deltas * deltas;
      if (kmh < kmh_min) kmh = kmh_min;
      real kmh_max = mix_upper_bound * dx / msftx(i, j) * dy / msfty(i, j) / dt;
      if (kmh > kmh_max) kmh = kmh_max;

      real kmv     = kmh;
      real kmv_max = mix_upper_bound / rdzw(i, j, k) / rdzw(i, j, k) / dt;
      if (kmv > kmv_max) kmv = kmv_max;

      real khh = kmh / pr;
      if (khh > kmh_max) khh = kmh_max;
      real khv = kmv / pr;
      if (khv > kmv_max) khv = kmv_max;

      xkmh(i, j, k) = kmh;
      xkmv(i, j, k) = kmv;
      xkhh(i, j, k) = khh;
      xkhv(i, j, k) = khv;
    }
  }
};

// -----------------------------------------------------------------------------
// cal_dampkm (lines 1353-1480, damp_opt = 1)
// -----------------------------------------------------------------------------
struct KCalDampKm {
  View3<const real> rdz, rdzw;
  BView2<const real> msftx, msfty;
  View3<real> xkmh, xkhh, xkmv, xkhv;
  real dx, dy, dt, dampcoef, zdamp;
  int i_start, ni, nz;

  CPPWRF_HD void operator()(int i, int j) const {
    constexpr real piconst  = real(3.1415926535897932384626433);
    constexpr real degrad   = piconst / real(180.0);
    constexpr real degrad90 = degrad * real(90.0);
    const int lane_i = i - i_start;

    real dampk[128], dampkv[128];
    real deltaz;
    {
      const int k   = nz - 1;
      real dx_m     = dx / msftx(i, j);
      real dy_m     = dy / msfty(i, j);
      real ds       = dx_m < dy_m ? dx_m : dy_m;
      real kmmax    = ds * ds / dt;
      real dz       = real(1.0) / rdzw(i, j, k);
      deltaz        = real(0.5) * dz;
      real kmmvmax  = dz * dz / dt;
      real ratio    = deltaz / zdamp;
      real tmp      = ratio < real(1.0) ? ratio : real(1.0);
      real c        = vcos_diff(degrad90 * tmp, lane_i, ni);
      dampk[k]      = c * c * kmmax * dampcoef;
      dampkv[k]     = c * c * kmmvmax * dampcoef;
      if (dampk[k] < dampkv[k]) dampkv[k] = dampk[k];
    }
    for (int k = nz - 2; k >= 0; --k) {
      real dx_m     = dx / msftx(i, j);
      real dy_m     = dy / msfty(i, j);
      real ds       = dx_m < dy_m ? dx_m : dy_m;
      real kmmax    = ds * ds / dt;
      real dz       = real(1.0) / rdz(i, j, k);
      deltaz        = deltaz + dz;
      dz            = real(1.0) / rdzw(i, j, k);
      real kmmvmax  = dz * dz / dt;
      real ratio    = deltaz / zdamp;
      real tmp      = ratio < real(1.0) ? ratio : real(1.0);
      real c        = vcos_diff(degrad90 * tmp, lane_i, ni);
      dampk[k]      = c * c * kmmax * dampcoef;
      dampkv[k]     = c * c * kmmvmax * dampcoef;
      if (dampk[k] < dampkv[k]) dampkv[k] = dampk[k];
    }
    for (int k = 0; k < nz; ++k) {
      if (dampk[k] > xkmh(i, j, k))  xkmh(i, j, k) = dampk[k];
      if (dampk[k] > xkhh(i, j, k))  xkhh(i, j, k) = dampk[k];
      if (dampkv[k] > xkmv(i, j, k)) xkmv(i, j, k) = dampkv[k];
      if (dampkv[k] > xkhv(i, j, k)) xkhv(i, j, k) = dampkv[k];
    }
  }
};

inline void cal_dampkm(const DiffusionParams& p,
                       real dx, real dy, real dt,
                       View3<const real> rdz,
                       View3<const real> rdzw,
                       BView2<const real> msftx,
                       BView2<const real> msfty,
                       View3<real> xkmh,
                       View3<real> xkhh,
                       View3<real> xkmv,
                       View3<real> xkhv,
                       int nx, int ny, int nz) {
  int i_start = 0;
  int i_end   = nx;
  int j_start = 0;
  int j_end   = ny;
  if (p.specified || p.nested) {
    i_start = std::max(i_start, p.spec_bdy_width - 1);
    i_end   = std::min(i_end, nx - p.spec_bdy_width + 1);
    j_start = std::max(j_start, p.spec_bdy_width - 1);
    j_end   = std::min(j_end, ny - p.spec_bdy_width + 1);
  }
  const int ni = i_end - i_start;
  column_for(Range2D{i_start, i_end, j_start, j_end},
             KCalDampKm{rdz, rdzw, msftx, msfty, xkmh, xkhh, xkmv, xkhv,
                        dx, dy, dt, p.dampcoef, p.zdamp, i_start, ni, nz});
}

// Existing calculate_km_kh overload (km_opt = 1, 4 without TKE/damp_opt = 1)
inline void calculate_km_kh(const DiffusionParams& p,
                            real dx, real dy,
                            View3<const real> defor11,
                            View3<const real> defor22,
                            View3<const real> defor12,
                            View3<const real> rdzw,
                            View3<const real> zx,
                            View3<const real> zy,
                            BView2<const real> msftx,
                            BView2<const real> msfty,
                            View3<real> xkmh,
                            View3<real> xkmv,
                            View3<real> xkhh,
                            View3<real> xkhv,
                            int nx, int ny, int nz) {
  if (p.km_opt == 1) {
    parallel_for(Range3D{0, nx, 0, ny, 0, nz},
                 KIsotropicKm{xkmh, xkmv, xkhh, xkhv, p.khdif, p.kvdif, p.prandtl});
  } else if (p.km_opt == 4) {
    int i_s = (p.open_xs || p.specified || p.nested) ? 1 : 0;
    int i_e = (p.open_xe || p.specified || p.nested) ? nx - 1 : nx;
    int j_s = (p.open_ys || p.specified || p.nested) ? 1 : 0;
    int j_e = (p.open_ye || p.specified || p.nested) ? ny - 1 : ny;
    if (p.periodic_x) {
      i_s = 0;
      i_e = nx;
    }

    parallel_for(Range3D{i_s, i_e, j_s, j_e, 0, nz},
                 KSmag2dKm{defor11, defor22, defor12, rdzw, zx, zy, msftx, msfty,
                           xkmh, xkmv, xkhh, xkhv, dx, dy, p.c_s, p.prandtl, p.diff_opt, nz});

    if (p.specified) {
      column_for(Range2D{0, 1, 0, ny}, [=] CPPWRF_HD (int, int j) {
        for (int k = 0; k < nz; ++k) {
          xkmh(0, j, k) = real(0);
          xkmh(nx - 1, j, k) = real(0);
          xkmv(0, j, k) = real(0);
          xkmv(nx - 1, j, k) = real(0);
          xkhh(0, j, k) = real(0);
          xkhh(nx - 1, j, k) = real(0);
          xkhv(0, j, k) = real(0);
          xkhv(nx - 1, j, k) = real(0);
        }
      });
      column_for(Range2D{0, nx, 0, 1}, [=] CPPWRF_HD (int i, int) {
        for (int k = 0; k < nz; ++k) {
          xkmh(i, 0, k) = real(0);
          xkmh(i, ny - 1, k) = real(0);
          xkmv(i, 0, k) = real(0);
          xkmv(i, ny - 1, k) = real(0);
          xkhh(i, 0, k) = real(0);
          xkhh(i, ny - 1, k) = real(0);
          xkhv(i, 0, k) = real(0);
          xkhv(i, ny - 1, k) = real(0);
        }
      });
    }
  }
}

// Full calculate_km_kh (supports km_opt = 1, 2, 3, 4 and damp_opt = 1)
inline void calculate_km_kh(const DiffusionParams& p,
                            real dx, real dy, real dt,
                            real cf1, real cf2, real cf3,
                            View3<const real> defor11,
                            View3<const real> defor22,
                            View3<const real> defor33,
                            View3<const real> defor12,
                            View3<const real> defor13,
                            View3<const real> defor23,
                            View3<const real> tke,
                            View3<const real> p8w,
                            View3<const real> t8w,
                            View3<const real> theta,
                            View3<const real> t_phy,
                            View3<const real> p_phy,
                            View3<const real> qv, bool has_qv,
                            View3<const real> qc, bool has_qc,
                            View3<const real> qi, bool has_qi,
                            View3<const real> rdz,
                            View3<const real> rdzw,
                            View3<const real> zx,
                            View3<const real> zy,
                            BView2<const real> msftx,
                            BView2<const real> msfty,
                            View3<real> bn2,
                            View3<real> xkmh,
                            View3<real> xkmv,
                            View3<real> xkhh,
                            View3<real> xkhv,
                            int nx, int ny, int nz) {
  calculate_N2(p, theta, t_phy, p_phy, p8w, t8w, rdz, rdzw,
               qv, has_qv, qc, has_qc, qi, has_qi, cf1, cf2, cf3,
               bn2, nx, ny, nz);

  View3<const real> bn2_c{bn2.p, bn2.sj, bn2.sk};
  if (p.km_opt == 1 || p.km_opt == 4) {
    calculate_km_kh(p, dx, dy, defor11, defor22, defor12, rdzw, zx, zy,
                    msftx, msfty, xkmh, xkmv, xkhh, xkhv, nx, ny, nz);
  } else if (p.km_opt == 2) {
    tke_km(p, dx, dy, dt, bn2_c, tke, p8w, t8w, theta, rdz, rdzw,
           msftx, msfty, xkmh, xkmv, xkhh, xkhv, nx, ny, nz);
  } else if (p.km_opt == 3) {
    int i_s = (p.open_xs || p.specified || p.nested) ? 1 : 0;
    int i_e = (p.open_xe || p.specified || p.nested) ? nx - 1 : nx;
    int j_s = (p.open_ys || p.specified || p.nested) ? 1 : 0;
    int j_e = (p.open_ye || p.specified || p.nested) ? ny - 1 : ny;
    if (p.periodic_x) {
      i_s = 0;
      i_e = nx;
    }
    const int ni = i_e - i_s;
    parallel_for(Range3D{i_s, i_e, j_s, j_e, 0, nz},
                 KSmagKm{defor11, defor22, defor33, defor12, defor13, defor23, bn2_c, rdzw,
                         msftx, msfty, xkmh, xkmv, xkhh, xkhv, dx, dy, dt,
                         p.c_s, p.prandtl, p.mix_upper_bound, p.isotropic, i_s, ni, nz});
  }

  if (p.damp_opt == 1) {
    cal_dampkm(p, dx, dy, dt, rdz, rdzw, msftx, msfty,
               xkmh, xkhh, xkmv, xkhv, nx, ny, nz);
  }
}

}  // namespace cppwrf::dyn::diffusion
