#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "phys/phys_constants.hpp"

namespace cppwrf::phys::kessler {

static constexpr int kMaxLevels = 128;

struct KesslerKernel {
  View3<real> t;
  View3<real> qv;
  View3<real> qc;
  View3<real> qr;
  View3<const real> rho;
  View3<const real> pii;
  View3<const real> dz8w;
  View3<const real> z;
  View2<real> rainnc;
  View2<real> rainncv;
  real dt;
  real xlv;
  real cp;
  real ep2;
  real svp1;
  real svp2;
  real svp3;
  real svpt0;
  real rhowater;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    using std::max;
    using std::min;
    using std::sqrt;
    using std::pow;
    using std::exp;
    using std::round;

    // Fixed-size local arrays per column (Fortran module_mp_kessler.F lines 65)
    real vt[kMaxLevels];
    real prodk[kMaxLevels];
    real vtden[kMaxLevels];
    real rdzk[kMaxLevels];
    real rhok[kMaxLevels];
    real factor[kMaxLevels];
    real rdzw[kMaxLevels];
    real prod[kMaxLevels];

    // Fortran line 81: f5 = svp2*(svpt0-svp3)*xlv/cp
    const real f5 = svp2 * (svpt0 - svp3) * xlv / cp;

    // Fortran line 23-28 parameters:
    constexpr real c1 = real(0.001);
    constexpr real c2 = real(0.001);
    constexpr real c3 = real(2.2);
    constexpr real c4 = real(0.875);
    constexpr real max_cr_sedimentation = real(0.75);

    // Fortran lines 102-124: Terminal velocity calculation and advection setup
    real crmax = real(0);
    for (int k = 0; k < nz; ++k) {
      prodk[k] = qr(i, j, k);
      rhok[k] = rho(i, j, k);
      real qrr = max(real(0), qr(i, j, k) * real(0.001) * rhok[k]);
      vtden[k] = sqrt(rhok[0] / rhok[k]);
      vt[k] = real(36.34) * pow(qrr, real(0.1364)) * vtden[k];
      rdzw[k] = real(1) / dz8w(i, j, k);
      crmax = max(vt[k] * dt * rdzw[k], crmax);
    }
    for (int k = 0; k < nz - 1; ++k) {
      rdzk[k] = real(1) / (z(i, j, k + 1) - z(i, j, k));
    }
    rdzk[nz - 1] = real(1) / (z(i, j, nz - 1) - z(i, j, nz - 2));

    // Fortran lines 125-128: Courant number and sedimentation timestep splitting
    int nfall = max(1, int(round(real(0.5) + crmax / max_cr_sedimentation)));
    real dtfall = dt / real(nfall);
    real time_sediment = dt;

    // Fortran lines 134-193: Column sedimentation loop
    while (nfall > 0) {
      time_sediment = time_sediment - dtfall;
      for (int k = 0; k < nz - 1; ++k) {
        factor[k] = dtfall * rdzk[k] / rhok[k];
      }
      factor[nz - 1] = dtfall * rdzk[nz - 1];

      // Surface rain fallout
      real ppt = rhok[0] * prodk[0] * vt[0] * dtfall / rhowater;
      rainncv(i, j) = ppt * real(1000.0);
      rainnc(i, j) = rainnc(i, j) + ppt * real(1000.0);

      // Fallout flux upstream
      for (int k = 0; k < nz - 1; ++k) {
        prodk[k] = prodk[k] - factor[k] * (rhok[k] * prodk[k] * vt[k] - rhok[k + 1] * prodk[k + 1] * vt[k + 1]);
      }
      prodk[nz - 1] = prodk[nz - 1] - factor[nz - 1] * prodk[nz - 1] * vt[nz - 1];

      if (nfall > 1) {
        nfall = nfall - 1;
        crmax = real(0);
        for (int k = 0; k < nz; ++k) {
          real qrr = max(real(0), prodk[k] * real(0.001) * rhok[k]);
          vt[k] = real(36.34) * pow(qrr, real(0.1364)) * vtden[k];
          crmax = max(vt[k] * time_sediment * rdzw[k], crmax);
        }
        int nfall_new = max(1, int(round(real(0.5) + crmax / max_cr_sedimentation)));
        if (nfall_new != nfall) {
          nfall = nfall_new;
          dtfall = time_sediment / real(nfall);
        }
      } else {
        for (int k = 0; k < nz; ++k) {
          prod[k] = prodk[k];
        }
        nfall = 0;
      }
    }

    // Fortran lines 203-236: Microphysics production & evaporation
    for (int k = 0; k < nz; ++k) {
      real factorn = real(1.0) / (real(1.0) + c3 * dt * pow(max(real(0.0), qr(i, j, k)), c4));
      real qrprod = qc(i, j, k) * (real(1.0) - factorn) + factorn * c1 * dt * max(qc(i, j, k) - c2, real(0.0));
      real rcgs = real(0.001) * rho(i, j, k);

      qc(i, j, k) = max(qc(i, j, k) - qrprod, real(0.0));
      qr(i, j, k) = (qr(i, j, k) + prod[k] - qr(i, j, k));
      qr(i, j, k) = max(qr(i, j, k) + qrprod, real(0.0));

      real temp = pii(i, j, k) * t(i, j, k);
      real pressure = real(1.000e+05) * pow(pii(i, j, k), real(1004.0 / 287.0));
      real gam = real(2.5e+06) / (real(1004.0) * pii(i, j, k));
      real es = real(1000.0) * svp1 * exp(svp2 * (temp - svpt0) / (temp - svp3));
      real qvs = ep2 * es / (pressure - es);
      prod[k] = (qv(i, j, k) - qvs) / (real(1.0) + pressure / (pressure - es) * qvs * f5 / ((temp - svp3) * (temp - svp3)));
      real dim_qvs_qv = max(real(0.0), qvs - qv(i, j, k));
      real ern_num = (real(1.6) + real(124.9) * pow(rcgs * qr(i, j, k), real(0.2046))) * pow(rcgs * qr(i, j, k), real(0.525));
      real ern_den = real(2.55e8) / (pressure * qvs) + real(5.4e5);
      real ern_term1 = dt * (ern_num / ern_den) * (dim_qvs_qv / (rcgs * qvs));
      real ern_term2 = max(-prod[k] - qc(i, j, k), real(0.0));
      real ern = min(ern_term1, min(ern_term2, qr(i, j, k)));

      real product = max(prod[k], -qc(i, j, k));
      t(i, j, k) = t(i, j, k) + gam * (product - ern);
      qv(i, j, k) = max(qv(i, j, k) - product + ern, real(0.0));
      qc(i, j, k) = qc(i, j, k) + product;
      qr(i, j, k) = qr(i, j, k) - ern;
    }
  }
};

void run_kessler(
    View3<real> t,
    View3<real> qv,
    View3<real> qc,
    View3<real> qr,
    View3<const real> rho,
    View3<const real> pii,
    View3<const real> dz8w,
    View3<const real> z,
    View2<real> rainnc,
    View2<real> rainncv,
    real dt,
    const Range2D& range,
    int nz,
    real xlv = real(cst::xlv),
    real cp = real(cst::cp),
    real ep2 = real(cst::ep_2),
    real svp1 = real(cst::svp1),
    real svp2 = real(cst::svp2),
    real svp3 = real(cst::svp3),
    real svpt0 = real(cst::svpt0),
    real rhowater = real(cst::rhowater));

}  // namespace cppwrf::phys::kessler
