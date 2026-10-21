#pragma once
// MYNN surface layer (WRF v4.8.0 phys/MYNN-SFC, sf_sfclay_physics = 5): constants, the psi tables of
// mynnsfc_driver (psi_init) and the helpers shared by module_sf_mynnsfc_land / _water / _ice (their copies
// are identical: zolrib, psi*_stable / _unstable, li_etal_2010). Default options only: psi_opt = 0.
// REAL is REAL(kind_phys) = REAL*8 in WRF's -r8 build; x**2 is x*x, x**real is pow.
#include <cmath>

#include "core/backend.hpp"
#include "core/precision.hpp"
#include "phys/driver/wrf_constants.hpp"

namespace cppwrf::phys::mynnsfc {

namespace wc = cppwrf::phys::wrfc;

// module_sf_mynnsfc_common / the modules' parameters
inline constexpr real grav = wc::g;
inline constexpr real cp = wc::cp;
inline constexpr real rd = wc::r_d;
inline constexpr real rovcp = wc::rcp;
inline constexpr real xlv = wc::xlv;
inline constexpr real ep1 = wc::r_v / wc::r_d - real(1.);  // p608 = rvovrd - 1.
inline constexpr real ep2 = wc::ep_2;
inline constexpr real ep3 = real(1.) - wc::ep_2;
inline constexpr real karman = wc::karman;
inline constexpr real wmin = real(0.1);
inline constexpr real svp1 = real(0.6112), svp2 = real(17.67), svp3 = real(29.65), svpt0 = real(273.15);
inline constexpr real vconvc = real(1.25);
inline constexpr real sqrt3 = real(1.7320508075688773);
inline constexpr real atan1 = real(0.785398163397);
inline constexpr real zero = 0, one = 1, two = 2, four = 4, ten = 10, twenty = 20;
inline constexpr real p01 = real(0.01), p1 = real(0.1), p25 = real(0.25), p5 = real(0.5), p75 = real(0.75);
inline constexpr real p9 = real(0.9), p99 = real(0.99);
inline constexpr real p333 = real(1.) / real(3.);

struct PsiTables {
  real psim_stab[1001], psim_unstab[1001], psih_stab[1001], psih_unstab[1001];
};

// mynnsfc_driver's psim/psih *_full (psi_opt = 0). (i, n): iteration and trip count of psi_init's loop, which
// gfortran vectorises (glibc libmvec pow / log / atan); i = 0, n = 1 is the scalar libm. ne: the trip count
// of the call sites that WRF's vectorised loop also takes through libmvec in the last (odd) iteration.
CPPWRF_HD inline real psim_stable_full(real zolf, int i = 0, int n = 1) {
  return real(-6.1) * vm::log(zolf + vm::pow(one + vm::pow(zolf, real(2.5), i, n), real(0.4), i, n), i, n);
}
CPPWRF_HD inline real psih_stable_full(real zolf, int i = 0, int n = 1, int ne = 1) {
  return real(-5.3) *
         vm::log(zolf + vm::pow(one + vm::pow(zolf, real(1.1), i, n), real(0.9090909090909090909), i, n), i, ne);
}
CPPWRF_HD inline real psim_unstable_full(real zolf, int i = 0, int n = 1, int ne = 1) {
  const real x = vm::pow(one - real(16.) * zolf, real(.25), i, n);
  const real psimk = two * vm::log(real(0.5) * (one + x), i, n) + vm::log(real(0.5) * (one + x * x), i, n) -
                     two * vm::atan(x, i, n) + two * atan1;
  const real ym = vm::pow(one - ten * zolf, p333, i, n);
  const real psimc = real(1.5) * vm::log((ym * ym + ym + one) * p333, i, n) -
                     sqrt3 * vm::atan((two * ym + 1) / sqrt3, i, ne) + four * atan1 / sqrt3;
  return (psimk + zolf * zolf * (psimc)) / (1 + zolf * zolf);
}
CPPWRF_HD inline real psih_unstable_full(real zolf, int i = 0, int n = 1) {
  const real y = std::sqrt(one - real(16.) * zolf);  // gfortran: x**.5 is sqrt(x)
  const real psihk = two * vm::log((one + y) * real(0.5), i, n);
  const real yh = vm::pow(one - real(34.) * zolf, p333, i, n);
  const real psihc = real(1.5) * vm::log((yh * yh + yh + one) * p333, i, n) -
                     sqrt3 * vm::atan((two * yh + one) / sqrt3, i, n) + four * atan1 / sqrt3;
  return (psihk + zolf * zolf * (psihc)) / (one + zolf * zolf);
}

// psi_init (psi_opt = 0): n = 0..1000, one vectorised loop of 1001 trips
inline void psi_init(PsiTables& t) {
  constexpr int nt = 1001;  // odd: the last iteration is scalar, except two call sites (nt + 1: always libmvec)
  for (int n = 0; n <= 1000; ++n) {
    real zolf = real(n) * real(0.01);
    t.psim_stab[n] = psim_stable_full(zolf, n, nt);
    t.psih_stab[n] = psih_stable_full(zolf, n, nt, nt + 1);
    zolf = -real(n) * real(0.01);
    t.psim_unstab[n] = psim_unstable_full(zolf, n, nt, nt + 1);
    t.psih_unstab[n] = psih_unstable_full(zolf, n, nt);
  }
}

// table lookups of the land / water / ice modules (scalar libm beyond the table)
CPPWRF_HD inline real psim_stable(real zolf, const PsiTables& t) {
  const int nzol = int(zolf * real(100.));
  const real rzol = zolf * real(100.) - nzol;
  if (nzol + 1 < 1000) return t.psim_stab[nzol] + rzol * (t.psim_stab[nzol + 1] - t.psim_stab[nzol]);
  return psim_stable_full(zolf);
}
CPPWRF_HD inline real psih_stable(real zolf, const PsiTables& t) {
  const int nzol = int(zolf * real(100.));
  const real rzol = zolf * real(100.) - nzol;
  if (nzol + 1 < 1000) return t.psih_stab[nzol] + rzol * (t.psih_stab[nzol + 1] - t.psih_stab[nzol]);
  return psih_stable_full(zolf);
}
CPPWRF_HD inline real psim_unstable(real zolf, const PsiTables& t) {
  const int nzol = int(-zolf * real(100.));
  const real rzol = -zolf * real(100.) - nzol;
  if (nzol + 1 < 1000) return t.psim_unstab[nzol] + rzol * (t.psim_unstab[nzol + 1] - t.psim_unstab[nzol]);
  return psim_unstable_full(zolf);
}
CPPWRF_HD inline real psih_unstable(real zolf, const PsiTables& t) {
  const int nzol = int(-zolf * real(100.));
  const real rzol = -zolf * real(100.) - nzol;
  if (nzol + 1 < 1000) return t.psih_unstab[nzol] + rzol * (t.psih_unstab[nzol + 1] - t.psih_unstab[nzol]);
  return psih_unstable_full(zolf);
}

// li_etal_2010
CPPWRF_HD inline void li_etal_2010(real& zl, real rib, real zaz0, real z0zt) {
  constexpr real au11 = real(0.045), bu11 = real(0.003), bu12 = real(0.0059), bu21 = real(-0.0828),
                 bu22 = real(0.8845), bu31 = real(0.1739), bu32 = real(-0.9213), bu33 = real(-0.1057);
  constexpr real aw11 = real(0.5738), aw12 = real(-0.4399), aw21 = real(-4.901), aw22 = real(52.50),
                 bw11 = real(-0.0539), bw12 = real(1.540), bw21 = real(-0.669), bw22 = real(-3.282);
  constexpr real as11 = real(0.7529), as21 = real(14.94), bs11 = real(0.1569), bs21 = real(-0.3091),
                 bs22 = real(-1.303);
  real zaz02 = zaz0;
  if (zaz0 < real(100.0)) zaz02 = real(100.);
  if (zaz0 > real(100000.0)) zaz02 = real(100000.);
  real z0zt2 = z0zt;
  if (z0zt < real(0.5)) z0zt2 = real(0.5);
  if (z0zt > real(100.0)) z0zt2 = real(100.);
  const real alfa = std::log(zaz02);
  const real beta = std::log(z0zt2);
  if (rib <= real(0.0)) {
    zl = au11 * alfa * (rib * rib) +
         ((bu11 * beta + bu12) * (alfa * alfa) + (bu21 * beta + bu22) * alfa + (bu31 * (beta * beta) + bu32 * beta + bu33)) * rib;
    zl = c_max(zl, real(-15.));
    zl = c_min(zl, real(0.));
  } else if (rib > real(0.0) && rib <= real(0.2)) {
    zl = ((aw11 * beta + aw12) * alfa + (aw21 * beta + aw22)) * (rib * rib) +
         ((bw11 * beta + bw12) * alfa + (bw21 * beta + bw22)) * rib;
    zl = c_min(zl, real(4.));
    zl = c_max(zl, real(0.));
  } else {
    zl = (as11 * alfa + as21) * rib + bs11 * alfa + bs21 * beta + bs22;
    zl = c_min(zl, real(20.));
    zl = c_max(zl, one);
  }
}

// zolrib (Pedro Jimenez's iteration for z/L from the bulk Richardson number)
CPPWRF_HD inline real zolrib(real ri, real za, real z0, real zt, real logz0, real logzt, real& zol1,
                             const PsiTables& t) {
  constexpr int nmax = 20;
  if (zol1 * ri < real(0.)) zol1 = zero;
  real zolold, res;
  if (ri < real(0.)) {
    zolold = real(-99999.);
    res = real(-66666.);
  } else {
    zolold = real(99999.);
    res = real(66666.);
  }
  int n = 1;
  while (std::fabs(zolold - res) > real(0.01) && n < nmax) {
    zolold = (n == 1) ? zol1 : res;
    const real zol20 = zolold * z0 / za;
    const real zol3 = zolold + zol20;
    const real zolt = zolold * zt / za;
    real psit2, psix2;
    if (ri < 0) {
      psit2 = c_max(logzt - (psih_unstable(zol3, t) - psih_unstable(zolt, t)), one);
      psix2 = c_max(logz0 - (psim_unstable(zol3, t) - psim_unstable(zol20, t)), one);
    } else {
      psit2 = c_max(logzt - (psih_stable(zol3, t) - psih_stable(zolt, t)), one);
      psix2 = c_max(logz0 - (psim_stable(zol3, t) - psim_stable(zol20, t)), one);
    }
    res = ri * (psix2 * psix2) / psit2;
    n = n + 1;
  }
  if (n == nmax && std::fabs(zolold - res) > real(0.01)) li_etal_2010(res, ri, za / z0, z0 / zt);
  return res;
}

}  // namespace cppwrf::phys::mynnsfc
