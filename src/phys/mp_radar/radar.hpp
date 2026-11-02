#pragma once
// WRF v4.8.0 radar reflectivity library: phys/module_mp_radar.F (Thompson) and phys/physics_mmm/mp_radar.F90
// (WSM6), whose routines are the same: radar_init's size bins and constants, rayleigh_soak_wetgraupel (Blahak's
// Rayleigh backscattering of melting snow / graupel) with get_m_mix_nested, get_m_mix and
// m_complex_maxwellgarnett for the only mixing WRF sets up (host 'air', matrix 'water', hostmatrix 'icewater',
// inclusions 'spheroidal', rule 'maxwellgarnett'). DOUBLE PRECISION and COMPLEX*16 are double here, as in WRF's
// -r8 -fdefault-double-8 build (a plain -r8 build computes them in REAL*16).
//
// Complex arithmetic as gfortran (-fcx-fortran-rules) lowers it: products (ar*br - ai*bi, ar*bi + ai*br), a real
// times a complex part by part, quotients by GCC's expand_complex_div_wide (Smith's method), SQRT / LOG / ABS by
// glibc's csqrt / clog / cabs on the host (std::complex calls them); the device uses its own formulas.
#include <cmath>
#include <memory>

#include "core/backend.hpp"
#if !defined(__CUDA_ARCH__)
#include <complex>
#endif

namespace cppwrf::phys::radar {

inline constexpr int nrbins = 50;
inline constexpr double lamda_radar = 0.10;
inline constexpr double melt_outside = 0.9;  // melt_outside_s, melt_outside_g
inline constexpr double PIx = 3.1415926535897932384626434;

struct Cplx {
  double re, im;
};

CPPWRF_HD inline Cplx cmul(Cplx a, Cplx b) { return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re}; }
CPPWRF_HD inline Cplx cscale(double r, Cplx a) { return {r * a.re, r * a.im}; }  // REAL * COMPLEX
CPPWRF_HD inline Cplx cadd(Cplx a, Cplx b) { return {a.re + b.re, a.im + b.im}; }
CPPWRF_HD inline Cplx csub(Cplx a, Cplx b) { return {a.re - b.re, a.im - b.im}; }
// GCC expand_complex_div_wide: the branch on |br| < |bi|
CPPWRF_HD inline Cplx cdiv(Cplx a, Cplx b) {
  double ratio, div, tr, ti;
  if (std::fabs(b.re) < std::fabs(b.im)) {
    ratio = b.re / b.im;
    div = b.re * ratio + b.im;
    tr = a.re * ratio + a.im;
    ti = a.im * ratio - a.re;
  } else {
    ratio = b.im / b.re;
    div = b.im * ratio + b.re;
    tr = a.im * ratio + a.re;
    ti = a.im - a.re * ratio;
  }
  return {tr / div, ti / div};
}
CPPWRF_HD inline double cabs_(Cplx z) {
#if defined(__CUDA_ARCH__)
  return ::hypot(z.re, z.im);
#else
  return std::abs(std::complex<double>(z.re, z.im));
#endif
}
CPPWRF_HD inline Cplx clog_(Cplx z) {
#if defined(__CUDA_ARCH__)
  return {::log(::hypot(z.re, z.im)), ::atan2(z.im, z.re)};
#else
  const std::complex<double> r = std::log(std::complex<double>(z.re, z.im));
  return {r.real(), r.imag()};
#endif
}
CPPWRF_HD inline Cplx csqrt_(Cplx z) {
#if defined(__CUDA_ARCH__)
  if (z.re == 0.0 && z.im == 0.0) return {0.0, z.im};
  const double t = ::sqrt(0.5 * (::hypot(z.re, z.im) + ::fabs(z.re)));
  if (z.re >= 0.0) return {t, z.im / (2.0 * t)};
  return {::fabs(z.im) / (2.0 * t), ::copysign(t, z.im)};
#else
  const std::complex<double> r = std::sqrt(std::complex<double>(z.re, z.im));
  return {r.real(), r.imag()};
#endif
}

// radar_init's results. REAL members are Real (WRF's default REAL), the DOUBLE PRECISION / COMPLEX*16 ones double.
template <typename Real>
struct RadarTables {
  double xxDs[nrbins], xdts[nrbins], xxDg[nrbins], xdtg[nrbins], simpson[nrbins + 1];
  double k_w, pi5, lamda4;
  Cplx m_w_0, m_i_0;
  Real xcre[4], xcse[4], xcge[4], xcrg[4], xcsg[4], xcgg[4];
  Real xam_r, xbm_r, xmu_r, xobmr;
  Real xam_s, xbm_s, xmu_s, xoams, xobms, xocms;
  Real xam_g, xbm_g, xmu_g, xoamg, xobmg, xocmg;
  Real xorg2, xosg2, xogg2;
};

// GAMMLN / WGAMMA (Numerical Recipes): DOUBLE PRECISION inside, REAL result
template <typename Real>
inline Real gammln(Real xx) {
  constexpr double stp = 2.5066282746310005;
  constexpr double cof[6] = {76.18009172947146, -86.50532032941677, 24.01409824083091,
                             -1.231739572450155, .1208650973866179e-2, -.5395239384953e-5};
  const double x = double(xx);
  double y = x;
  double tmp = x + 5.5;
  tmp = (x + 0.5) * std::log(tmp) - tmp;
  double ser = 1.000000000190015;
  for (int j = 0; j < 6; ++j) {
    y = y + 1.0;
    ser = ser + cof[j] / y;
  }
  return Real(tmp + std::log(stp * ser / x));
}
// WGAMMA = EXP(GAMMLN(y)); radar_init's 4-value loops are vectorised: libmvec exp on the pairs 1-2, 3-4
template <typename Real>
inline Real wgamma(Real y) { return vm::exp(gammln(y), 0, 2); }

// radar_init (host). The caller's m(D) = am D**bm and shape parameters mu of rain, snow and graupel set the
// x* constants (mp_radar.F90's form: xcre(3) = 4 + mu, xcre(4) = 7 + mu, i.e. bm = 3).
// gfortran folds the constant parts at compile time, correctly rounded (MPFR / MPC): PI5, lamda4, m_w_0 =
// m_complex_water_ray(0.1 m, 0 C), m_i_0 = m_complex_ice_maetzler(0.1 m, 0 C), K_w, xxDx(2) and the exp
// arguments of the bin loops (log(xxDx(51) / xxDx(1)), log(xxDx(1))); those are the constants of WRF's
// module_mp_radar.o / mp_radar.o. Bins 3 ... 50 come from libmvec's 2-lane exp (pairs 3-4, ..., 49-50).
template <typename Real>
inline RadarTables<Real> radar_init(Real xam_r, Real xbm_r, Real xmu_r, Real xam_s, Real xbm_s, Real xmu_s,
                                    Real xam_g, Real xbm_g, Real xmu_g) {
  RadarTables<Real> t{};
  t.pi5 = 3.14159 * 3.14159 * 3.14159 * 3.14159 * 3.14159;
  t.lamda4 = lamda_radar * lamda_radar * lamda_radar * lamda_radar;
  t.m_w_0 = {0x1.213d1a5047315p+3, -0x1.634d6f1c00d6dp+0};
  t.m_i_0 = {0x1.c91dae5095597p+0, -0x1.1acb617ef252cp-13};
  t.k_w = 0x1.de4b5a94b4103p-1;  // (ABS((m_w_0**2 - 1) / (m_w_0**2 + 2)))**2
  constexpr double basis[3] = {1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0};
  for (int n = 0; n <= nrbins; ++n) t.simpson[n] = 0.0;
  for (int n = 0; n < nrbins - 1; n += 2) {
    t.simpson[n] = t.simpson[n] + basis[0];
    t.simpson[n + 1] = t.simpson[n + 1] + basis[1];
    t.simpson[n + 2] = t.simpson[n + 2] + basis[2];
  }
  // snow bins 100 microns ... 2 cm, graupel 100 microns ... 5 cm: xxDx(n) = exp((n-1)/nrbins * log(xxDx(51) /
  // xxDx(1)) + log(xxDx(1)))
  auto bins = [](double top, double log_ratio, double x2, double* xxD, double* xdt) {
    constexpr double log_x1 = -0x1.26bb1bbb55515p+3;  // log(100.d-6)
    double xxDx[nrbins + 1];
    xxDx[0] = 100.e-6;
    xxDx[nrbins] = top;
    xxDx[1] = x2;
    for (int n = 3; n <= nrbins; ++n)
      xxDx[n - 1] = vm::exp(double(n - 1) / double(nrbins) * log_ratio + log_x1, 0, 2);
    for (int n = 0; n < nrbins; ++n) {
      xxD[n] = std::sqrt(xxDx[n] * xxDx[n + 1]);
      xdt[n] = xxDx[n + 1] - xxDx[n];
    }
  };
  bins(0.02, 0x1.5317a1b949c53p+2, 0x1.d250f2923b4dap-14, t.xxDs, t.xdts);
  bins(0.05, 0x1.8dbc239b0b862p+2, 0x1.daf0cbbe98755p-14, t.xxDg, t.xdtg);

  t.xam_r = xam_r; t.xbm_r = xbm_r; t.xmu_r = xmu_r;
  t.xam_s = xam_s; t.xbm_s = xbm_s; t.xmu_s = xmu_s;
  t.xam_g = xam_g; t.xbm_g = xbm_g; t.xmu_g = xmu_g;
  t.xcre[0] = Real(1.) + xbm_r;
  t.xcre[1] = Real(1.) + xmu_r;
  t.xcre[2] = Real(4.) + xmu_r;
  t.xcre[3] = Real(7.) + xmu_r;
  for (int n = 0; n < 4; ++n) t.xcrg[n] = wgamma(t.xcre[n]);
  t.xorg2 = Real(1.) / t.xcrg[1];
  t.xcse[0] = Real(1.) + xbm_s;
  t.xcse[1] = Real(1.) + xmu_s;
  t.xcse[2] = Real(4.) + xmu_s;
  t.xcse[3] = Real(7.) + xmu_s;
  for (int n = 0; n < 4; ++n) t.xcsg[n] = wgamma(t.xcse[n]);
  t.xosg2 = Real(1.) / t.xcsg[1];
  t.xcge[0] = Real(1.) + xbm_g;
  t.xcge[1] = Real(1.) + xmu_g;
  t.xcge[2] = Real(4.) + xmu_g;
  t.xcge[3] = Real(7.) + xmu_g;
  for (int n = 0; n < 4; ++n) t.xcgg[n] = wgamma(t.xcge[n]);
  t.xogg2 = Real(1.) / t.xcgg[1];
  t.xobmr = Real(1.) / xbm_r;
  t.xoams = Real(1.) / xam_s;
  t.xobms = Real(1.) / xbm_s;
  t.xocms = std::pow(t.xoams, t.xobms);
  t.xoamg = Real(1.) / xam_g;
  t.xobmg = Real(1.) / xbm_g;
  t.xocmg = std::pow(t.xoamg, t.xobmg);
  return t;
}

// The tables in device memory (one copy per domain's microphysics)
template <typename Real>
class RadarDeviceTables {
 public:
  explicit RadarDeviceTables(const RadarTables<Real>& h) {
    p_ = static_cast<RadarTables<Real>*>(dev_alloc(sizeof(RadarTables<Real>)));
    copy_to_dev(p_, &h, sizeof(RadarTables<Real>));
  }
  ~RadarDeviceTables() { dev_free(p_); }
  RadarDeviceTables(const RadarDeviceTables&) = delete;
  RadarDeviceTables& operator=(const RadarDeviceTables&) = delete;
  const RadarTables<Real>* get() const { return p_; }

 private:
  RadarTables<Real>* p_ = nullptr;
};

// m_complex_maxwellgarnett, inclusion 'spheroidal': m1 the matrix, m2 / m3 the inclusions. Returns false on
// WRF's error (the volume fractions do not sum to 1)
CPPWRF_HD inline bool maxwellgarnett(double vol1, double vol2, double vol3, Cplx m1, Cplx m2, Cplx m3, Cplx& m) {
  if (std::fabs(vol1 + vol2 + vol3 - 1.0) > 1e-6) return false;
  const Cplx m1t = cmul(m1, m1), m2t = cmul(m2, m2), m3t = cmul(m3, m3);
  const Cplx d2 = csub(m2t, m1t), d3 = csub(m3t, m1t);
  Cplx l2 = cmul(cdiv(m2t, d2), clog_(cdiv(m2t, m1t)));
  l2.re = l2.re - 1.0;
  Cplx l3 = cmul(cdiv(m3t, d3), clog_(cdiv(m3t, m1t)));
  l3.re = l3.re - 1.0;
  const Cplx beta2 = cmul(cdiv(cscale(2.0, m1t), d2), l2);
  const Cplx beta3 = cmul(cdiv(cscale(2.0, m1t), d3), l3);
  const double v1 = 1.0 - vol2 - vol3;
  const Cplx num = cadd(cadd(cscale(v1, m1t), cmul(cscale(vol2, beta2), m2t)), cmul(cscale(vol3, beta3), m3t));
  Cplx den = cscale(vol2, beta2);
  den.re = v1 + den.re;  // REAL + COMPLEX
  den = cadd(den, cscale(vol3, beta3));
  m = csqrt_(cdiv(num, den));
  return true;
}

// rayleigh_soak_wetgraupel: backscattering cross section C_back of a melting particle of mass x_g with the
// geometry D = a_geo * x_g**b_geo, melted fraction fmelt (90 % of the meltwater outside)
CPPWRF_HD inline double rayleigh_soak_wetgraupel(double x_g, double a_geo, double b_geo, double fmelt, Cplx m_w,
                                                 Cplx m_i, double pi5, double lamda4) {
  const Cplx m_air{1.0, 0.0};
  const double fm = c_max(c_min(fmelt, 1.0), 0.0);
  double mra = c_max(c_min(melt_outside, 1.0), 0.0);
  mra = mra + (1.0 - mra) * fm;
  const double x_w = x_g * fm;
  const double D_g = a_geo * std::pow(x_g, b_geo);
  if (!(D_g >= 1e-12)) return 0.0;
  double vg = PIx / 6. * ((D_g * D_g) * D_g);
  const double rhog = c_max(c_min(x_g / vg, 900.0), 10.0);
  vg = x_g / rhog;
  const double meltratio_outside_grenz = 1.0 - rhog / 1000.;
  double volg;
  if (mra <= meltratio_outside_grenz) {
    volg = vg * (1.0 - mra * fm);
  } else {
    const double fmgrenz = (900.0 - rhog) / (mra * 900.0 - rhog + 900.0 * rhog / 1000.);
    if (fm <= fmgrenz) volg = (1.0 - mra * fm) * vg;
    else volg = (x_g - x_w) / 900.0 + x_w / 1000.;
  }
  const double D_large = std::pow(6.0 / PIx * volg, 1. / 3.);
  const double volice = (x_g - x_w) / (volg * 900.0);
  const double volwater = x_w / (1000. * volg);
  const double volair = 1.0 - volice - volwater;
  // get_m_mix_nested, host 'air': water matrix with ice inclusions, then air inclusions in that mixture
  const double vol1 = volice / c_max(volice + volwater, 1e-10);
  const double vol2 = 1.0 - vol1;
  Cplx mtmp, m_core;
  const bool ok1 = maxwellgarnett(vol2, 0.0, vol1, m_w, m_air, m_i, mtmp);
  if (!ok1) mtmp = {-999.99, -999.99};
  const bool ok2 = maxwellgarnett(1.0 - volair, volair, 0.0, mtmp, m_air, cscale(2.0, m_air), m_core);
  if (!ok1 || !ok2) return 0.0;
  const Cplx m2 = cmul(m_core, m_core);
  const double a = cabs_(cdiv({m2.re - 1.0, m2.im}, {m2.re + 2.0, m2.im}));
  const double d3 = (D_large * D_large) * D_large;  // D_large**6 = (D**3)**2 (GCC's powi chain)
  return a * a * pi5 * (d3 * d3) / lamda4;
}

}  // namespace cppwrf::phys::radar
