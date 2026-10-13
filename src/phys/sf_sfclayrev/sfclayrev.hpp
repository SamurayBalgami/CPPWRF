#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "phys/phys_constants.hpp"

namespace cppwrf::phys::sfclayrev {

template <typename T>
struct SFClayRevTables {
  T psim_stab[1001];
  T psim_unstab[1001];
  T psih_stab[1001];
  T psih_unstab[1001];
};

template <typename T>
CPPWRF_HD inline T psim_stable_full(T zolf) {
  using std::log;
  using std::pow;
  return static_cast<T>(-6.1) * log(zolf + pow(static_cast<T>(1.0) + pow(zolf, static_cast<T>(2.5)), static_cast<T>(1.0 / 2.5)));
}

template <typename T>
CPPWRF_HD inline T psih_stable_full(T zolf) {
  using std::log;
  using std::pow;
  return static_cast<T>(-5.3) * log(zolf + pow(static_cast<T>(1.0) + pow(zolf, static_cast<T>(1.1)), static_cast<T>(1.0 / 1.1)));
}

template <typename T>
CPPWRF_HD inline T psim_unstable_full(T zolf) {
  using std::atan;
  using std::log;
  using std::pow;
  using std::sqrt;
  T x = pow(static_cast<T>(1.0) - static_cast<T>(16.0) * zolf, static_cast<T>(0.25));
  T psimk = static_cast<T>(2.0) * log(static_cast<T>(0.5) * (static_cast<T>(1.0) + x)) +
            log(static_cast<T>(0.5) * (static_cast<T>(1.0) + x * x)) -
            static_cast<T>(2.0) * atan(x) + static_cast<T>(2.0) * atan(static_cast<T>(1.0));

  T ym = pow(static_cast<T>(1.0) - static_cast<T>(10.0) * zolf, static_cast<T>(0.33));
  T sqrt3 = sqrt(static_cast<T>(3.0));
  T psimc = (static_cast<T>(3.0) / static_cast<T>(2.0)) * log((ym * ym + ym + static_cast<T>(1.0)) / static_cast<T>(3.0)) -
            sqrt3 * atan((static_cast<T>(2.0) * ym + static_cast<T>(1.0)) / sqrt3) +
            static_cast<T>(4.0) * atan(static_cast<T>(1.0)) / sqrt3;

  return (psimk + zolf * zolf * psimc) / (static_cast<T>(1.0) + zolf * zolf);
}

template <typename T>
CPPWRF_HD inline T psih_unstable_full(T zolf) {
  using std::atan;
  using std::log;
  using std::pow;
  using std::sqrt;
  T y = sqrt(static_cast<T>(1.0) - static_cast<T>(16.0) * zolf);
  T psihk = static_cast<T>(2.0) * log((static_cast<T>(1.0) + y) * static_cast<T>(0.5));

  T yh = pow(static_cast<T>(1.0) - static_cast<T>(34.0) * zolf, static_cast<T>(0.33));
  T sqrt3 = sqrt(static_cast<T>(3.0));
  T psihc = (static_cast<T>(3.0) / static_cast<T>(2.0)) * log((yh * yh + yh + static_cast<T>(1.0)) / static_cast<T>(3.0)) -
            sqrt3 * atan((static_cast<T>(2.0) * yh + static_cast<T>(1.0)) / sqrt3) +
            static_cast<T>(4.0) * atan(static_cast<T>(1.0)) / sqrt3;

  return (psihk + zolf * zolf * psihc) / (static_cast<T>(1.0) + zolf * zolf);
}

template <typename T>
CPPWRF_HD inline T psim_stable(T zolf, const SFClayRevTables<T>& tbl) {
  int nzol = static_cast<int>(zolf * static_cast<T>(100.0));
  T rzol = zolf * static_cast<T>(100.0) - static_cast<T>(nzol);
  if (nzol + 1 < 1000) {
    return tbl.psim_stab[nzol] + rzol * (tbl.psim_stab[nzol + 1] - tbl.psim_stab[nzol]);
  } else {
    return psim_stable_full(zolf);
  }
}

template <typename T>
CPPWRF_HD inline T psih_stable(T zolf, const SFClayRevTables<T>& tbl) {
  int nzol = static_cast<int>(zolf * static_cast<T>(100.0));
  T rzol = zolf * static_cast<T>(100.0) - static_cast<T>(nzol);
  if (nzol + 1 < 1000) {
    return tbl.psih_stab[nzol] + rzol * (tbl.psih_stab[nzol + 1] - tbl.psih_stab[nzol]);
  } else {
    return psih_stable_full(zolf);
  }
}

template <typename T>
CPPWRF_HD inline T psim_unstable(T zolf, const SFClayRevTables<T>& tbl) {
  int nzol = static_cast<int>(-zolf * static_cast<T>(100.0));
  T rzol = -zolf * static_cast<T>(100.0) - static_cast<T>(nzol);
  if (nzol + 1 < 1000) {
    return tbl.psim_unstab[nzol] + rzol * (tbl.psim_unstab[nzol + 1] - tbl.psim_unstab[nzol]);
  } else {
    return psim_unstable_full(zolf);
  }
}

template <typename T>
CPPWRF_HD inline T psih_unstable(T zolf, const SFClayRevTables<T>& tbl) {
  int nzol = static_cast<int>(-zolf * static_cast<T>(100.0));
  T rzol = -zolf * static_cast<T>(100.0) - static_cast<T>(nzol);
  if (nzol + 1 < 1000) {
    return tbl.psih_unstab[nzol] + rzol * (tbl.psih_unstab[nzol + 1] - tbl.psih_unstab[nzol]);
  } else {
    return psih_unstable_full(zolf);
  }
}

template <typename T>
CPPWRF_HD inline T zolri2(T& zol2, T ri2, T z, T z0, const SFClayRevTables<T>& tbl) {
  using std::log;
  if (zol2 * ri2 < static_cast<T>(0.0)) {
    zol2 = static_cast<T>(0.0);
  }
  T zol20 = zol2 * z0 / z;
  T zol3 = zol2 + zol20;

  T psix2, psih2;
  if (ri2 < static_cast<T>(0.0)) {
    psix2 = log((z + z0) / z0) - (psim_unstable(zol3, tbl) - psim_unstable(zol20, tbl));
    psih2 = log((z + z0) / z0) - (psih_unstable(zol3, tbl) - psih_unstable(zol20, tbl));
  } else {
    psix2 = log((z + z0) / z0) - (psim_stable(zol3, tbl) - psim_stable(zol20, tbl));
    psih2 = log((z + z0) / z0) - (psih_stable(zol3, tbl) - psih_stable(zol20, tbl));
  }
  return zol2 * psih2 / (psix2 * psix2) - ri2;
}

template <typename T>
CPPWRF_HD inline T zolri(T ri, T z, T z0, const SFClayRevTables<T>& tbl) {
  using std::abs;
  T x1, x2;
  if (ri < static_cast<T>(0.0)) {
    x1 = static_cast<T>(-5.0);
    x2 = static_cast<T>(0.0);
  } else {
    x1 = static_cast<T>(0.0);
    x2 = static_cast<T>(5.0);
  }

  T fx1 = zolri2(x1, ri, z, z0, tbl);
  T fx2 = zolri2(x2, ri, z, z0, tbl);
  T zolri_val = (abs(fx2) < abs(fx1)) ? x1 : x2;
  int iter = 0;
  while (abs(x1 - x2) > static_cast<T>(0.01)) {
    if (iter == 10) return zolri_val;
    if (fx1 == fx2) return zolri_val;
    if (abs(fx2) < abs(fx1)) {
      x1 = x1 - fx1 / (fx2 - fx1) * (x2 - x1);
      fx1 = zolri2(x1, ri, z, z0, tbl);
      zolri_val = x1;
    } else {
      x2 = x2 - fx2 / (fx2 - fx1) * (x2 - x1);
      fx2 = zolri2(x2, ri, z, z0, tbl);
      zolri_val = x2;
    }
    iter++;
  }
  return zolri_val;
}

template <typename T>
CPPWRF_HD inline T depth_dependent_z0(T water_depth, T z0, T ust) {
  using std::exp;
  using std::log;
  using std::min;
  T effective_depth;
  if (water_depth < static_cast<T>(10.0)) {
    effective_depth = static_cast<T>(10.0);
  } else if (water_depth > static_cast<T>(100.0)) {
    effective_depth = static_cast<T>(100.0);
  } else {
    effective_depth = water_depth;
  }
  T depth_b = static_cast<T>(1.0 / 30.0) * log(static_cast<T>(1260.0) / effective_depth);
  T res = exp((static_cast<T>(2.7) * ust - static_cast<T>(1.8) / depth_b) /
              (ust + static_cast<T>(0.17) / depth_b));
  return min(res, static_cast<T>(0.1));
}

template <typename T>
void init_sfclayrev_tables(SFClayRevTables<T>& tbl) {
  for (int n = 0; n <= 1000; ++n) {
    T zolf = static_cast<T>(n) * static_cast<T>(0.01);
    tbl.psim_stab[n] = psim_stable_full(zolf);
    tbl.psih_stab[n] = psih_stable_full(zolf);

    zolf = -static_cast<T>(n) * static_cast<T>(0.01);
    tbl.psim_unstab[n] = psim_unstable_full(zolf);
    tbl.psih_unstab[n] = psih_unstable_full(zolf);
  }
}

struct SFClayRevKernel {
  // Inputs (level 0 slices of 3D or 2D)
  View2<const real> u1d;
  View2<const real> v1d;
  View2<const real> t1d;
  View2<const real> qv1d;
  View2<const real> p1d;
  View2<const real> dz8w1d;
  View2<const real> psfcpa;
  View2<const real> tsk;
  View2<const real> xland;
  View2<const real> dx;
  View2<const real> mavail;
  View2<const real> pblh;
  View2<const real> lakemask;
  View2<const real> water_depth;

  // Inout
  View2<real> regime;
  View2<real> hfx;
  View2<real> qfx;
  View2<real> qsfc;
  View2<real> mol;
  View2<real> rmol;
  View2<real> gz1oz0;
  View2<real> wspd;
  View2<real> br;
  View2<real> psim;
  View2<real> psih;
  View2<real> fm;
  View2<real> fh;
  View2<real> znt;
  View2<real> zol;
  View2<real> ust;
  View2<real> cpm;
  View2<real> chs2;
  View2<real> cqs2;
  View2<real> chs;
  View2<real> flhc;
  View2<real> flqc;
  View2<real> qgh;
  View2<real> ustm;

  // Outputs
  View2<real> lh;
  View2<real> u10;
  View2<real> v10;
  View2<real> th2;
  View2<real> t2;
  View2<real> q2;
  View2<real> ck;
  View2<real> cka;
  View2<real> cd;
  View2<real> cda;

  // Constants & options
  const SFClayRevTables<real>* tables;
  bool isfflx;
  bool shalwater_z0;
  int isftcflx;
  int iz0tlnd;
  bool scm_force_flux;

  real cp;
  real g;
  real rovcp;
  real r;
  real xlv;
  real svp1;
  real svp2;
  real svp3;
  real svpt0;
  real ep1;
  real ep2;
  real karman;
  real p1000mb;

  CPPWRF_HD void operator()(int i, int j) const {
    using std::abs;
    using std::exp;
    using std::log;
    using std::max;
    using std::min;
    using std::pow;
    using std::sqrt;

    const SFClayRevTables<real>& tbl = *tables;

    constexpr real vconvc = real(1.0);
    constexpr real czo = real(0.0185);
    constexpr real ozo = real(1.59e-5);
    constexpr real xka = real(2.4e-5);
    constexpr real prt = real(1.0);
    constexpr real salinity_factor = real(0.98);

    // WRF sf_sfclayrev.F90 lines 219-222
    real psfc_cb = psfcpa(i, j) / real(1000.0);

    // WRF lines 226-231: Convert ground temp to potential temp
    real tgdsa = tsk(i, j);
    real thgb = tsk(i, j) * pow(p1000mb / psfcpa(i, j), rovcp);

    // WRF lines 253-263
    real pl = p1d(i, j) / real(1000.0);
    real scr3 = t1d(i, j);
    real thcon = pow(p1000mb * real(0.001) / pl, rovcp);
    real thx = scr3 * thcon;
    real scr4 = scr3;
    real thvx = thx;
    real qx = real(0.0);

    // WRF lines 265-270
    qgh(i, j) = real(0.0);
    flhc(i, j) = real(0.0);
    flqc(i, j) = real(0.0);
    cpm(i, j) = cp;

    // WRF lines 273-278
    qx = qv1d(i, j);
    real tvcon = real(1.0) + ep1 * qx;
    thvx = thx * tvcon;
    scr4 = scr3 * tvcon;

    // WRF lines 280-292
    real e1 = svp1 * exp(svp2 * (tgdsa - svpt0) / (tgdsa - svp3));
    if (xland(i, j) > real(1.5) && lakemask(i, j) == real(0.0)) {
      e1 = e1 * salinity_factor;
    }
    if (xland(i, j) > real(1.5) || qsfc(i, j) <= real(0.0)) {
      qsfc(i, j) = ep2 * e1 / (psfc_cb - e1);
    }
    e1 = svp1 * exp(svp2 * (t1d(i, j) - svpt0) / (t1d(i, j) - svp3));
    pl = p1d(i, j) / real(1000.0);
    qgh(i, j) = ep2 * e1 / (pl - e1);
    cpm(i, j) = cp * (real(1.0) + real(0.8) * qx);

    // WRF lines 298-313
    real zqklp1 = real(0.0);
    real rhox = psfc_cb * real(1000.0) / (r * scr4);
    real zqkl = dz8w1d(i, j) + zqklp1;
    real za = real(0.5) * (zqkl + zqklp1);
    real govrth = g / thx;

    // WRF lines 317-358
    gz1oz0(i, j) = log((za + znt(i, j)) / znt(i, j));
    real gz2oz0 = log((real(2.0) + znt(i, j)) / znt(i, j));
    real gz10oz0 = log((real(10.0) + znt(i, j)) / znt(i, j));

    real zl;
    if ((xland(i, j) - real(1.5)) >= real(0.0)) {
      zl = znt(i, j);
    } else {
      zl = real(0.01);
    }

    wspd(i, j) = sqrt(u1d(i, j) * u1d(i, j) + v1d(i, j) * v1d(i, j));
    real tskv = thgb * (real(1.0) + ep1 * qsfc(i, j));
    real dthvdz = thvx - tskv;

    real vconv;
    if (xland(i, j) < real(1.5)) {
      real fluxc = max(hfx(i, j) / rhox / cp + ep1 * tskv * qfx(i, j) / rhox, real(0.0));
      vconv = vconvc * pow(g / tgdsa * pblh(i, j) * fluxc, real(0.33));
    } else {
      real dthvm = (-dthvdz >= real(0.0)) ? -dthvdz : real(0.0);
      vconv = sqrt(dthvm);
    }

    real vsgd = real(0.32) * pow(max(dx(i, j) / real(5000.0) - real(1.0), real(0.0)), real(0.33));
    wspd(i, j) = sqrt(wspd(i, j) * wspd(i, j) + vconv * vconv + vsgd * vsgd);
    wspd(i, j) = max(wspd(i, j), real(0.1));
    br(i, j) = govrth * za * dthvdz / (wspd(i, j) * wspd(i, j));
    if (mol(i, j) < real(0.0)) {
      br(i, j) = min(br(i, j), real(0.0));
    }
    rmol(i, j) = -govrth * dthvdz * za * karman;

    // WRF lines 379-401
    zol(i, j) = real(0.0);
    if (br(i, j) > real(0.0)) {
      if (br(i, j) > real(250.0)) {
        zol(i, j) = zolri(real(250.0), za, znt(i, j), tbl);
      } else {
        zol(i, j) = zolri(br(i, j), za, znt(i, j), tbl);
      }
    }
    if (br(i, j) < real(0.0)) {
      if (ust(i, j) < real(0.001)) {
        zol(i, j) = br(i, j) * gz1oz0(i, j);
      } else {
        if (br(i, j) < real(-250.0)) {
          zol(i, j) = zolri(real(-250.0), za, znt(i, j), tbl);
        } else {
          zol(i, j) = zolri(br(i, j), za, znt(i, j), tbl);
        }
      }
    }

    // WRF lines 405-416
    real zolzz = zol(i, j) * (za + znt(i, j)) / za;
    real zol10 = zol(i, j) * (real(10.0) + znt(i, j)) / za;
    real zol2 = zol(i, j) * (real(2.0) + znt(i, j)) / za;
    real zol0 = zol(i, j) * znt(i, j) / za;
    real zl2 = (real(2.0)) / za * zol(i, j);
    real zl10 = (real(10.0)) / za * zol(i, j);

    if ((xland(i, j) - real(1.5)) < real(0.0)) {
      zl = (real(0.01)) / za * zol(i, j);
    } else {
      zl = zol0;
    }

    real psim10, psih10, psim2, psih2, pq, pq2, pq10;

    // WRF lines 418-500: Regimes
    if (br(i, j) > real(0.0)) {
      // Class 1: Stable
      regime(i, j) = real(1.0);
      psim(i, j) = psim_stable(zolzz, tbl) - psim_stable(zol0, tbl);
      psih(i, j) = psih_stable(zolzz, tbl) - psih_stable(zol0, tbl);
      psim10 = psim_stable(zol10, tbl) - psim_stable(zol0, tbl);
      psih10 = psih_stable(zol10, tbl) - psih_stable(zol0, tbl);
      psim2 = psim_stable(zol2, tbl) - psim_stable(zol0, tbl);
      psih2 = psih_stable(zol2, tbl) - psih_stable(zol0, tbl);

      pq = psih_stable(zol(i, j), tbl) - psih_stable(zl, tbl);
      pq2 = psih_stable(zl2, tbl) - psih_stable(zl, tbl);
      pq10 = psih_stable(zl10, tbl) - psih_stable(zl, tbl);
      rmol(i, j) = zol(i, j) / za;
    } else if (br(i, j) == real(0.0)) {
      // Class 3: Forced convection
      regime(i, j) = real(3.0);
      psim(i, j) = real(0.0);
      psih(i, j) = real(0.0);
      psim10 = real(0.0);
      psih10 = real(0.0);
      psim2 = real(0.0);
      psih2 = real(0.0);

      pq = real(0.0);
      pq2 = real(0.0);
      pq10 = real(0.0);
      zol(i, j) = real(0.0);
      rmol(i, j) = real(0.0);
    } else {
      // Class 4: Free convection
      regime(i, j) = real(4.0);
      psim(i, j) = psim_unstable(zolzz, tbl) - psim_unstable(zol0, tbl);
      psih(i, j) = psih_unstable(zolzz, tbl) - psih_unstable(zol0, tbl);
      psim10 = psim_unstable(zol10, tbl) - psim_unstable(zol0, tbl);
      psih10 = psih_unstable(zol10, tbl) - psih_unstable(zol0, tbl);
      psim2 = psim_unstable(zol2, tbl) - psim_unstable(zol0, tbl);
      psih2 = psih_unstable(zol2, tbl) - psih_unstable(zol0, tbl);

      pq = psih_unstable(zol(i, j), tbl) - psih_unstable(zl, tbl);
      pq2 = psih_unstable(zl2, tbl) - psih_unstable(zl, tbl);
      pq10 = psih_unstable(zl10, tbl) - psih_unstable(zl, tbl);

      psih(i, j) = min(psih(i, j), real(0.9) * gz1oz0(i, j));
      psim(i, j) = min(psim(i, j), real(0.9) * gz1oz0(i, j));
      psih2 = min(psih2, real(0.9) * gz2oz0);
      psim10 = min(psim10, real(0.9) * gz10oz0);
      psih10 = min(psih10, real(0.9) * gz10oz0);
      rmol(i, j) = zol(i, j) / za;
    }

    // WRF lines 506-528
    real dtg = thx - thgb;
    real psix = gz1oz0(i, j) - psim(i, j);
    real psix10 = gz10oz0 - psim10;
    real psit = gz1oz0(i, j) - psih(i, j);
    real psit2 = gz2oz0 - psih2;

    if ((xland(i, j) - real(1.5)) >= real(0.0)) {
      zl = znt(i, j);
    } else {
      zl = real(0.01);
    }

    real psiq = log(karman * ust(i, j) * za / xka + za / zl) - pq;
    real psiq2 = log(karman * ust(i, j) * real(2.0) / xka + real(2.0) / zl) - pq2;
    real psiq10 = log(karman * ust(i, j) * real(10.0) / xka + real(10.0) / zl) - pq10;

    // WRF lines 531-587: Fairall 2003 over water
    if ((xland(i, j) - real(1.5)) >= real(0.0)) {
      real visc = (real(1.32) + real(0.009) * (scr3 - real(273.15))) * real(1.0e-5);
      real restar = ust(i, j) * znt(i, j) / visc;
      real z0t = (real(5.5e-5)) * pow(restar, real(-0.60));
      z0t = min(z0t, real(1.0e-4));
      z0t = max(z0t, real(2.0e-9));
      real z0q = z0t;

      zolzz = zol(i, j) * (za + z0t) / za;
      zol10 = zol(i, j) * (real(10.0) + z0t) / za;
      zol2 = zol(i, j) * (real(2.0) + z0t) / za;
      zol0 = zol(i, j) * z0t / za;

      if (zol(i, j) > real(0.0)) {
        psih(i, j) = psih_stable(zolzz, tbl) - psih_stable(zol0, tbl);
        psih10 = psih_stable(zol10, tbl) - psih_stable(zol0, tbl);
        psih2 = psih_stable(zol2, tbl) - psih_stable(zol0, tbl);
      } else {
        if (zol(i, j) == real(0.0)) {
          psih(i, j) = real(0.0);
          psih10 = real(0.0);
          psih2 = real(0.0);
        } else {
          psih(i, j) = psih_unstable(zolzz, tbl) - psih_unstable(zol0, tbl);
          psih10 = psih_unstable(zol10, tbl) - psih_unstable(zol0, tbl);
          psih2 = psih_unstable(zol2, tbl) - psih_unstable(zol0, tbl);
        }
      }
      psit = log((za + z0t) / z0t) - psih(i, j);
      psit2 = log((real(2.0) + z0t) / z0t) - psih2;

      zolzz = zol(i, j) * (za + z0q) / za;
      zol10 = zol(i, j) * (real(10.0) + z0q) / za;
      zol2 = zol(i, j) * (real(2.0) + z0q) / za;
      zol0 = zol(i, j) * z0q / za;

      if (zol(i, j) > real(0.0)) {
        psih(i, j) = psih_stable(zolzz, tbl) - psih_stable(zol0, tbl);
        psih10 = psih_stable(zol10, tbl) - psih_stable(zol0, tbl);
        psih2 = psih_stable(zol2, tbl) - psih_stable(zol0, tbl);
      } else {
        if (zol(i, j) == real(0.0)) {
          psih(i, j) = real(0.0);
          psih10 = real(0.0);
          psih2 = real(0.0);
        } else {
          psih(i, j) = psih_unstable(zolzz, tbl) - psih_unstable(zol0, tbl);
          psih10 = psih_unstable(zol10, tbl) - psih_unstable(zol0, tbl);
          psih2 = psih_unstable(zol2, tbl) - psih_unstable(zol0, tbl);
        }
      }

      psiq = log((za + z0q) / z0q) - psih(i, j);
      psiq2 = log((real(2.0) + z0q) / z0q) - psih2;
      psiq10 = log((real(10.0) + z0q) / z0q) - psih10;

      // WRF lines 589-626: isftcflx == 1
      if (isftcflx == 1) {
        z0q = real(1.e-4);
        zolzz = zol(i, j) * (za + z0q) / za;
        zol10 = zol(i, j) * (real(10.0) + z0q) / za;
        zol2 = zol(i, j) * (real(2.0) + z0q) / za;
        zol0 = zol(i, j) * z0q / za;

        if (zol(i, j) > real(0.0)) {
          psih(i, j) = psih_stable(zolzz, tbl) - psih_stable(zol0, tbl);
          psih10 = psih_stable(zol10, tbl) - psih_stable(zol0, tbl);
          psih2 = psih_stable(zol2, tbl) - psih_stable(zol0, tbl);
        } else {
          if (zol(i, j) == real(0.0)) {
            psih(i, j) = real(0.0);
            psih10 = real(0.0);
            psih2 = real(0.0);
          } else {
            psih(i, j) = psih_unstable(zolzz, tbl) - psih_unstable(zol0, tbl);
            psih10 = psih_unstable(zol10, tbl) - psih_unstable(zol0, tbl);
            psih2 = psih_unstable(zol2, tbl) - psih_unstable(zol0, tbl);
          }
        }
        psiq = log((za + z0q) / z0q) - psih(i, j);
        psit = psiq;
        psiq2 = log((real(2.0) + z0q) / z0q) - psih2;
        psiq10 = log((real(10.0) + z0q) / z0q) - psih10;
        psit2 = psiq2;
      }

      // WRF lines 627-699: isftcflx == 2
      if (isftcflx == 2) {
        visc = (real(1.32) + real(0.009) * (scr3 - real(273.15))) * real(1.0e-5);
        restar = ust(i, j) * znt(i, j) / visc;
        real gz0ozt = real(0.40) * (real(7.3) * sqrt(sqrt(restar)) * sqrt(real(0.71)) - real(5.0));
        z0t = znt(i, j) / exp(gz0ozt);

        zolzz = zol(i, j) * (za + z0t) / za;
        zol10 = zol(i, j) * (real(10.0) + z0t) / za;
        zol2 = zol(i, j) * (real(2.0) + z0t) / za;
        zol0 = zol(i, j) * z0t / za;

        if (zol(i, j) > real(0.0)) {
          psih(i, j) = psih_stable(zolzz, tbl) - psih_stable(zol0, tbl);
          psih10 = psih_stable(zol10, tbl) - psih_stable(zol0, tbl);
          psih2 = psih_stable(zol2, tbl) - psih_stable(zol0, tbl);
        } else {
          if (zol(i, j) == real(0.0)) {
            psih(i, j) = real(0.0);
            psih10 = real(0.0);
            psih2 = real(0.0);
          } else {
            psih(i, j) = psih_unstable(zolzz, tbl) - psih_unstable(zol0, tbl);
            psih10 = psih_unstable(zol10, tbl) - psih_unstable(zol0, tbl);
            psih2 = psih_unstable(zol2, tbl) - psih_unstable(zol0, tbl);
          }
        }

        psit = log((za + z0t) / z0t) - psih(i, j);
        psit2 = log((real(2.0) + z0t) / z0t) - psih2;

        real gz0ozq = real(0.40) * (real(7.3) * sqrt(sqrt(restar)) * sqrt(real(0.60)) - real(5.0));
        z0q = znt(i, j) / exp(gz0ozq);

        zolzz = zol(i, j) * (za + z0q) / za;
        zol10 = zol(i, j) * (real(10.0) + z0q) / za;
        zol2 = zol(i, j) * (real(2.0) + z0q) / za;
        zol0 = zol(i, j) * z0q / za;

        if (zol(i, j) > real(0.0)) {
          psih(i, j) = psih_stable(zolzz, tbl) - psih_stable(zol0, tbl);
          psih10 = psih_stable(zol10, tbl) - psih_stable(zol0, tbl);
          psih2 = psih_stable(zol2, tbl) - psih_stable(zol0, tbl);
        } else {
          if (zol(i, j) == real(0.0)) {
            psih(i, j) = real(0.0);
            psih10 = real(0.0);
            psih2 = real(0.0);
          } else {
            psih(i, j) = psih_unstable(zolzz, tbl) - psih_unstable(zol0, tbl);
            psih10 = psih_unstable(zol10, tbl) - psih_unstable(zol0, tbl);
            psih2 = psih_unstable(zol2, tbl) - psih_unstable(zol0, tbl);
          }
        }

        psiq = log((za + z0q) / z0q) - psih(i, j);
        psiq2 = log((real(2.0) + z0q) / z0q) - psih2;
        psiq10 = log((real(10.0) + z0q) / z0q) - psih10;
      }
    }

    // WRF lines 700-705: ck, cd, cka, cda
    ck(i, j) = (karman / psix10) * (karman / psiq10);
    cd(i, j) = (karman / psix10) * (karman / psix10);
    cka(i, j) = (karman / psix) * (karman / psiq);
    cda(i, j) = (karman / psix) * (karman / psix);

    // WRF lines 706-756: iz0tlnd over land
    if (iz0tlnd >= 1 && (xland(i, j) - real(1.5)) <= real(0.0)) {
      zl = znt(i, j);
      real visc = (real(1.32) + real(0.009) * (scr3 - real(273.15))) * real(1.0e-5);
      real restar = ust(i, j) * zl / visc;
      real czil = real(0.0);
      if (iz0tlnd == 1) {
        czil = pow(real(10.0), real(-0.40) * (zl / real(0.07)));
      } else if (iz0tlnd == 2) {
        czil = real(0.1);
      }

      real z0t = znt(i, j) / exp(czil * karman * sqrt(restar));
      zolzz = zol(i, j) * (za + z0t) / za;
      zol10 = zol(i, j) * (real(10.0) + z0t) / za;
      zol2 = zol(i, j) * (real(2.0) + z0t) / za;
      zol0 = zol(i, j) * z0t / za;

      if (zol(i, j) > real(0.0)) {
        psih(i, j) = psih_stable(zolzz, tbl) - psih_stable(zol0, tbl);
        psih10 = psih_stable(zol10, tbl) - psih_stable(zol0, tbl);
        psih2 = psih_stable(zol2, tbl) - psih_stable(zol0, tbl);
      } else {
        if (zol(i, j) == real(0.0)) {
          psih(i, j) = real(0.0);
          psih10 = real(0.0);
          psih2 = real(0.0);
        } else {
          psih(i, j) = psih_unstable(zolzz, tbl) - psih_unstable(zol0, tbl);
          psih10 = psih_unstable(zol10, tbl) - psih_unstable(zol0, tbl);
          psih2 = psih_unstable(zol2, tbl) - psih_unstable(zol0, tbl);
        }
      }

      psiq = log((za + z0t) / z0t) - psih(i, j);
      psiq2 = log((real(2.0) + z0t) / z0t) - psih2;
      psit = psiq;
      psit2 = psiq2;
    }

    // WRF lines 757-780: ust update and diagnostics
    ust(i, j) = real(0.5) * ust(i, j) + real(0.5) * karman * wspd(i, j) / psix;
    real wspdi = sqrt(u1d(i, j) * u1d(i, j) + v1d(i, j) * v1d(i, j));
    ustm(i, j) = real(0.5) * ustm(i, j) + real(0.5) * karman * wspdi / psix;

    u10(i, j) = u1d(i, j) * psix10 / psix;
    v10(i, j) = v1d(i, j) * psix10 / psix;
    th2(i, j) = thgb + dtg * psit2 / psit;
    q2(i, j) = qsfc(i, j) + (qx - qsfc(i, j)) * psiq2 / psiq;
    t2(i, j) = th2(i, j) * pow(psfcpa(i, j) / p1000mb, rovcp);

    if ((xland(i, j) - real(1.5)) < real(0.0)) {
      ust(i, j) = max(ust(i, j), real(0.001));
    }
    mol(i, j) = karman * dtg / psit / prt;
    real denomq = psiq;
    real denomq2 = psiq2;
    real denomt2 = psit2;
    fm(i, j) = psix;
    fh(i, j) = psit;

    // WRF lines 784-792: Surface fluxes initialization
    if (!scm_force_flux) {
      qfx(i, j) = real(0.0);
      hfx(i, j) = real(0.0);
    }

    if (isfflx) {
      // WRF lines 797-847
      if ((xland(i, j) - real(1.5)) >= real(0.0)) {
        if (shalwater_z0) {
          znt(i, j) = depth_dependent_z0(water_depth(i, j), znt(i, j), ust(i, j));
        } else {
          znt(i, j) = czo * ust(i, j) * ust(i, j) / g + real(0.11) * real(1.5e-5) / ust(i, j);
          znt(i, j) = min(znt(i, j), real(2.85e-3));
        }

        if (isftcflx != 0) {
          real zw = min(pow(ust(i, j) / real(1.06), real(0.3)), real(1.0));
          real zn1 = real(0.011) * ust(i, j) * ust(i, j) / g + ozo;
          real zn2 = real(10.0) * exp(real(-9.5) * pow(ust(i, j), real(-0.3333))) +
                     real(0.11) * real(1.5e-5) / max(ust(i, j), real(0.01));
          znt(i, j) = (real(1.0) - zw) * zn1 + zw * zn2;
          znt(i, j) = min(znt(i, j), real(2.85e-3));
          znt(i, j) = max(znt(i, j), real(1.27e-7));
        }
        zl = znt(i, j);
      } else {
        zl = real(0.01);
      }

      flqc(i, j) = rhox * mavail(i, j) * ust(i, j) * karman / denomq;
      real dtthx = abs(thx - thgb);
      if (dtthx > real(1.e-5)) {
        flhc(i, j) = cpm(i, j) * rhox * ust(i, j) * mol(i, j) / (thx - thgb);
      } else {
        flhc(i, j) = real(0.0);
      }

      if (!scm_force_flux) {
        qfx(i, j) = flqc(i, j) * (qsfc(i, j) - qx);
        lh(i, j) = xlv * qfx(i, j);
        hfx(i, j) = flhc(i, j) * (thgb - thx);
      }

      chs(i, j) = ust(i, j) * karman / denomq;
      cqs2(i, j) = ust(i, j) * karman / denomq2;
      chs2(i, j) = ust(i, j) * karman / denomt2;
    }
  }
};

void run_sfclayrev(
    View2<const real> u1d,
    View2<const real> v1d,
    View2<const real> t1d,
    View2<const real> qv1d,
    View2<const real> p1d,
    View2<const real> dz8w1d,
    View2<const real> psfcpa,
    View2<const real> tsk,
    View2<const real> xland,
    View2<const real> dx,
    View2<const real> mavail,
    View2<const real> pblh,
    View2<const real> lakemask,
    View2<const real> water_depth,
    View2<real> regime,
    View2<real> hfx,
    View2<real> qfx,
    View2<real> qsfc,
    View2<real> mol,
    View2<real> rmol,
    View2<real> gz1oz0,
    View2<real> wspd,
    View2<real> br,
    View2<real> psim,
    View2<real> psih,
    View2<real> fm,
    View2<real> fh,
    View2<real> znt,
    View2<real> zol,
    View2<real> ust,
    View2<real> cpm,
    View2<real> chs2,
    View2<real> cqs2,
    View2<real> chs,
    View2<real> flhc,
    View2<real> flqc,
    View2<real> qgh,
    View2<real> ustm,
    View2<real> lh,
    View2<real> u10,
    View2<real> v10,
    View2<real> th2,
    View2<real> t2,
    View2<real> q2,
    View2<real> ck,
    View2<real> cka,
    View2<real> cd,
    View2<real> cda,
    const SFClayRevTables<real>* tables,
    const Range2D& range,
    bool isfflx = true,
    bool shalwater_z0 = false,
    int isftcflx = 0,
    int iz0tlnd = 0,
    bool scm_force_flux = false,
    real cp = real(cst::cp),
    real g = real(cst::g),
    real rovcp = real(cst::rcp),
    real r = real(cst::r_d),
    real xlv = real(cst::xlv),
    real svp1 = real(cst::svp1),
    real svp2 = real(cst::svp2),
    real svp3 = real(cst::svp3),
    real svpt0 = real(cst::svpt0),
    real ep1 = real(cst::ep_1),
    real ep2 = real(cst::ep_2),
    real karman = real(cst::karman),
    real p1000mb = real(cst::p1000mb));

}  // namespace cppwrf::phys::sfclayrev
