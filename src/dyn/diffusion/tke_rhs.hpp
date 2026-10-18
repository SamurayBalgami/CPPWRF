#pragma once

#include <cmath>
#include "core/backend.hpp"
#include "core/constants.hpp"
#include "core/field.hpp"
#include "dyn/diffusion/diffusion_params.hpp"
#include "dyn/diffusion/km_kh.hpp"

namespace cppwrf::dyn::diffusion {

// Ported from WRF v4.8.0 dyn_em/module_diffusion_em.F:
//   tke_rhs      (lines 6099-6229)
//   tke_buoyancy (lines 6234-6379)
//   tke_dissip   (lines 6384-6524)
//   tke_shear    (lines 6529-6877)

struct KTkeRhs {
  View3<const real> defor11, defor22, defor33, defor12, defor13, defor23;
  View3<const real> u, v, tke, bn2, theta, rho, rdzw, xkmh, xkmv, xkhv;
  View3<const real> qv;
  View2<const real> mu, ust, hfx;
  BView2<const real> msftx, msfty;
  View1<const real> c1, c2;
  bool has_qv, has_ust, has_hfx;
  real dx, dy, dt, c_k, tke_drag_coefficient, tke_heat_flux;
  int isfflx, nz;
  View3<real> tendency;

  CPPWRF_HD void operator()(int i, int j) const {
    constexpr real g       = real(9.81);
    constexpr real r_d     = real(287.0);
    constexpr real cp      = real(7.0) * r_d / real(2.0);
    constexpr real epsilon = real(1.0e-15);

    // -------------------------------------------------------------------------
    // 1. tke_shear (lines 6529-6877)
    // -------------------------------------------------------------------------
    // defor11
    for (int k = 0; k < nz; ++k) {
      tendency(i, j, k) = tendency(i, j, k) +
                          real(0.5) * (c1(k) * mu(i, j) + c2(k)) * xkmh(i, j, k) *
                          (defor11(i, j, k) * defor11(i, j, k));
    }
    // defor22
    for (int k = 0; k < nz; ++k) {
      tendency(i, j, k) = tendency(i, j, k) +
                          real(0.5) * (c1(k) * mu(i, j) + c2(k)) * xkmh(i, j, k) *
                          (defor22(i, j, k) * defor22(i, j, k));
    }
    // defor33
    for (int k = 0; k < nz; ++k) {
      tendency(i, j, k) = tendency(i, j, k) +
                          real(0.5) * (c1(k) * mu(i, j) + c2(k)) * xkmv(i, j, k) *
                          (defor33(i, j, k) * defor33(i, j, k));
    }
    // defor12
    for (int k = 0; k < nz; ++k) {
      real avg12 = real(0.25) * ((defor12(i, j, k) * defor12(i, j, k)) +
                                 (defor12(i, j + 1, k) * defor12(i, j + 1, k)) +
                                 (defor12(i + 1, j, k) * defor12(i + 1, j, k)) +
                                 (defor12(i + 1, j + 1, k) * defor12(i + 1, j + 1, k)));
      tendency(i, j, k) = tendency(i, j, k) + (c1(k) * mu(i, j) + c2(k)) * xkmh(i, j, k) * avg12;
    }
    // defor13 (tmp2 has k=0 and k=nz zeroed)
    for (int k = 0; k < nz; ++k) {
      real d13_0_kp1 = (k + 1 == nz) ? real(0.0) : defor13(i, j, k + 1);
      real d13_0_k   = (k == 0)      ? real(0.0) : defor13(i, j, k);
      real d13_1_kp1 = (k + 1 == nz) ? real(0.0) : defor13(i + 1, j, k + 1);
      real d13_1_k   = (k == 0)      ? real(0.0) : defor13(i + 1, j, k);
      real avg13 = real(0.25) * ((d13_0_kp1 * d13_0_kp1) + (d13_0_k * d13_0_k) +
                                 (d13_1_kp1 * d13_1_kp1) + (d13_1_k * d13_1_k));
      tendency(i, j, k) = tendency(i, j, k) + (c1(k) * mu(i, j) + c2(k)) * xkmv(i, j, k) * avg13;
    }
    // uflux at k = 0
    {
      const int k = 0;
      real u_sum = u(i, j, k) + u(i + 1, j, k);
      real v_sum = v(i, j, k) + v(i, j + 1, k);
      if (isfflx == 0) {
        real absU = real(0.5) * std::sqrt(u_sum * u_sum + v_sum * v_sum);
        real Cd   = tke_drag_coefficient;
        tendency(i, j, k) = tendency(i, j, k) +
            (c1(k) * mu(i, j) + c2(k)) *
            (u_sum * real(0.5) * Cd * absU * (defor13(i, j, 1) + defor13(i + 1, j, 1)) * real(0.5));
      } else {
        real absU   = real(0.5) * std::sqrt(u_sum * u_sum + v_sum * v_sum) + epsilon;
        real ust_ij = has_ust ? ust(i, j) : real(0.0);
        real Cd     = (ust_ij * ust_ij) / (absU * absU);
        tendency(i, j, k) = tendency(i, j, k) +
            (c1(k) * mu(i, j) + c2(k)) *
            (u_sum * real(0.5) * Cd * absU * (defor13(i, j, 1) + defor13(i + 1, j, 1)) * real(0.5));
      }
    }
    // defor23 (tmp2 has k=0 and k=nz zeroed)
    for (int k = 0; k < nz; ++k) {
      real d23_0_kp1 = (k + 1 == nz) ? real(0.0) : defor23(i, j, k + 1);
      real d23_0_k   = (k == 0)      ? real(0.0) : defor23(i, j, k);
      real d23_1_kp1 = (k + 1 == nz) ? real(0.0) : defor23(i, j + 1, k + 1);
      real d23_1_k   = (k == 0)      ? real(0.0) : defor23(i, j + 1, k);
      real avg23 = real(0.25) * ((d23_0_kp1 * d23_0_kp1) + (d23_0_k * d23_0_k) +
                                 (d23_1_kp1 * d23_1_kp1) + (d23_1_k * d23_1_k));
      tendency(i, j, k) = tendency(i, j, k) + (c1(k) * mu(i, j) + c2(k)) * xkmv(i, j, k) * avg23;
    }
    // vflux at k = 0
    {
      const int k = 0;
      real u_sum = u(i, j, k) + u(i + 1, j, k);
      real v_sum = v(i, j, k) + v(i, j + 1, k);
      if (isfflx == 0) {
        real absU = real(0.5) * std::sqrt(u_sum * u_sum + v_sum * v_sum);
        real Cd   = tke_drag_coefficient;
        tendency(i, j, k) = tendency(i, j, k) +
            (c1(k) * mu(i, j) + c2(k)) *
            (v_sum * real(0.5) * Cd * absU * (defor23(i, j, 1) + defor23(i, j + 1, 1)) * real(0.5));
      } else {
        real absU   = real(0.5) * std::sqrt(u_sum * u_sum + v_sum * v_sum) + epsilon;
        real ust_ij = has_ust ? ust(i, j) : real(0.0);
        real Cd     = (ust_ij * ust_ij) / (absU * absU);
        tendency(i, j, k) = tendency(i, j, k) +
            (c1(k) * mu(i, j) + c2(k)) *
            (v_sum * real(0.5) * Cd * absU * (defor23(i, j, 1) + defor23(i, j + 1, 1)) * real(0.5));
      }
    }

    // -------------------------------------------------------------------------
    // 2. tke_buoyancy (lines 6234-6379)
    // -------------------------------------------------------------------------
    for (int k = 1; k < nz; ++k) {
      tendency(i, j, k) = tendency(i, j, k) - (c1(k) * mu(i, j) + c2(k)) * xkhv(i, j, k) * bn2(i, j, k);
    }
    {
      const int k = 0;
      real heat_flux;
      if (isfflx == 0 || isfflx == 2) {
        heat_flux = tke_heat_flux;
      } else {
        real qv_0   = has_qv ? qv(i, j, 0) : real(0.0);
        real hfx_ij = has_hfx ? hfx(i, j) : real(0.0);
        real cpm    = cp * (real(1.0) + real(0.8) * qv_0);
        heat_flux   = (hfx_ij / cpm) / rho(i, j, 0);
      }
      tendency(i, j, k) = tendency(i, j, k) -
          (c1(k) * mu(i, j) + c2(k)) *
          ((xkhv(i, j, k) * bn2(i, j, k)) - (g / theta(i, j, k)) * heat_flux) / real(2.0);
    }

    // -------------------------------------------------------------------------
    // 3. tke_dissip (lines 6384-6524)
    // -------------------------------------------------------------------------
    real ce1 = (c_k / real(0.10)) * real(0.19);
    real ce2 = (real(0.93) - ce1) > real(0.0) ? (real(0.93) - ce1) : real(0.0);
    real l_scale[128];
    calc_l_scale_col(i, j, nz, dx, dy, msftx, msfty, tke, bn2, rdzw, l_scale);

    for (int k = 0; k < nz; ++k) {
      real deltas = rpow(dx / msftx(i, j) * dy / msfty(i, j) / rdzw(i, j, k), real(0.33333333));
      real tketmp = tke(i, j, k) > real(1.0e-6) ? tke(i, j, k) : real(1.0e-6);
      real coefc  = (k == 0 || k == nz - 1) ? real(3.9) : (ce1 + ce2 * l_scale[k] / deltas);
      tendency(i, j, k) = tendency(i, j, k) -
                          (c1(k) * mu(i, j) + c2(k)) * coefc * rpow(tketmp, real(1.5)) / l_scale[k];
    }

    // -------------------------------------------------------------------------
    // 4. Lower limit on TKE tendency (lines 6221-6227)
    // -------------------------------------------------------------------------
    for (int k = 0; k < nz; ++k) {
      real tke_pos = tke(i, j, k) > real(0.0) ? tke(i, j, k) : real(0.0);
      real lim     = -(c1(k) * mu(i, j) + c2(k)) * tke_pos / dt;
      if (tendency(i, j, k) < lim) tendency(i, j, k) = lim;
    }
  }
};

inline void tke_rhs(const DiffusionParams& p,
                    real dx, real dy, real dt,
                    View3<const real> defor11,
                    View3<const real> defor22,
                    View3<const real> defor33,
                    View3<const real> defor12,
                    View3<const real> defor13,
                    View3<const real> defor23,
                    View3<const real> u,
                    View3<const real> v,
                    View3<const real> tke,
                    View3<const real> bn2,
                    View3<const real> theta,
                    View3<const real> rho,
                    View3<const real> rdzw,
                    View3<const real> xkmh,
                    View3<const real> xkmv,
                    View3<const real> xkhv,
                    View3<const real> qv, bool has_qv,
                    View2<const real> mu,
                    View1<const real> c1,
                    View1<const real> c2,
                    BView2<const real> msftx,
                    BView2<const real> msfty,
                    View2<const real> ust, bool has_ust,
                    View2<const real> hfx, bool has_hfx,
                    View3<real> tendency,
                    int nx, int ny, int nz) {
  int i_start = (p.open_xs || p.specified || p.nested) ? 1 : 0;
  int i_end   = (p.open_xe || p.specified || p.nested) ? nx - 1 : nx;
  int j_start = (p.open_ys || p.specified || p.nested) ? 1 : 0;
  int j_end   = (p.open_ye || p.specified || p.nested) ? ny - 1 : ny;
  if (p.periodic_x) {
    i_start = 0;
    i_end   = nx;
  }

  column_for(Range2D{i_start, i_end, j_start, j_end},
             KTkeRhs{defor11, defor22, defor33, defor12, defor13, defor23,
                     u, v, tke, bn2, theta, rho, rdzw, xkmh, xkmv, xkhv,
                     qv, mu, ust, hfx, msftx, msfty, c1, c2,
                     has_qv, has_ust, has_hfx,
                     dx, dy, dt, p.c_k, p.tke_drag_coefficient, p.tke_heat_flux,
                     p.isfflx, nz, tendency});
}

}  // namespace cppwrf::dyn::diffusion
