// SPDX-License-Identifier: GPL-3.0-or-later
#include "phys/mp_p3/p3.hpp"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "core/limits.hpp"

namespace cppwrf::phys::mp_p3 {

// p3_init (module_mp_p3.F) for nCat = 1, trplMomI = .false.: the module constants, lookupTable_1
// p3_lookupTable_1.dat-v5.4_2momI (WRF reads it from the run directory, '.') and the rain fall-speed and ventilation
// tables (vn_table, vm_table, revap_table; mu_r = mu_r_constant = 0 for every one of the 10 mu_r columns, so the
// columns are equal and computed once)
template <typename Real>
P3Host<Real> p3_init(const std::string& dir) {
  P3Host<Real> h;
  P3Params<Real>& c = h.c;
  c.pi = Real(3.14159265);
  c.thrd = Real(1.) / Real(3.);
  const Real sxth = Real(1.) / Real(6.);
  const Real piov3 = c.pi * c.thrd;
  const Real piov6 = c.pi * sxth;
  c.max_total_Ni = Real(2000.e+3);
  c.nccnst = Real(200.e+6);
  c.cp = Real(1005.);
  c.inv_cp = Real(1.) / c.cp;
  c.g = Real(9.816);
  c.rd = Real(287.15);
  c.rv = Real(461.51);
  c.ep_2 = Real(0.622);
  c.rhosur = Real(100000.) / (c.rd * Real(273.15));
  c.rhosui = Real(60000.) / (c.rd * Real(253.15));
  c.f1r = Real(0.78);
  c.f2r = Real(0.32);
  c.rhow = Real(1000.);
  c.cpw = Real(4218.);
  c.inv_rhow = Real(1.) / c.rhow;
  c.mu_r_constant = Real(0.);
  c.inv_Drmax = Real(1.) / Real(0.002);
  c.rho_rimeMin = Real(50.);
  c.rho_rimeMax = Real(900.);
  c.inv_rho_rimeMax = Real(1.) / c.rho_rimeMax;
  c.qsmall = Real(1.e-14);
  c.nsmall = Real(1.e-16);
  c.bsmall = c.qsmall * c.inv_rho_rimeMax;
  const Real bimm = Real(2.);
  c.aimm = Real(0.65);
  c.mi0 = Real(4.) * piov3 * Real(900.) * Real(1.e-18);
  c.eci = Real(0.5);
  c.eri = Real(1.);
  c.bcn = Real(2.);
  c.nmltratio = Real(1.);
  c.e0 = polysvp1(Real(273.15), 0);
  c.cons1 = piov6 * c.rhow;
  const Real cons2 = Real(4.) * piov3 * c.rhow;
  // (25.e-6)**3 and (1.e-6)**3: constants gfortran folds correctly rounded
  c.cons3 = Real(1.) / (cons2 * Real(0x1.19799812dea12p-46));
  c.cons5 = piov6 * bimm;
  c.cons6 = piov6 * piov6 * c.rhow * bimm;
  c.cons7 = Real(4.) * piov3 * c.rhow * Real(0x1.2725dd1d243abp-60);
  c.mw = Real(0.018);
  const Real osm = Real(1.), vi = Real(3.), epsm = Real(0.9), map = Real(0.132);
  const Real rhoa = Real(1777.);
  c.rr = Real(8.3187);
  c.bact = vi * osm * epsm * c.mw * rhoa / (map * c.rhow);
  c.inv_rm1 = Real(2.e+7);
  c.sig1 = Real(2.0);
  c.nanew1 = Real(300.e6);
  c.inv_rm2 = Real(7.6923076e+5);
  c.sig2 = Real(2.5);
  c.nanew2 = Real(0.);

  // lookupTable_1 (the version line, a blank line, then per density and rime-fraction index the 50 ice lines and the
  // 50 x 30 ice-rain collection lines)
  const std::string path = dir.empty() ? std::string("p3_lookupTable_1.dat-v5.4_2momI")
                                       : dir + "/p3_lookupTable_1.dat-v5.4_2momI";
  std::ifstream f(path);
  if (!f) throw std::runtime_error("P3: cannot open " + path + " (WRF's run directory file)");
  std::string line;
  std::getline(f, line);
  {
    std::istringstream s(line);
    std::string dumstr, version;
    s >> dumstr >> version;
    if (version != "5.4_2momI")
      throw std::runtime_error("P3: " + path + " is version '" + version + "', P3 v4.5.2 needs 5.4_2momI (WRF stops)");
  }
  std::getline(f, line);
  h.itab.assign(kItab, Real(0));
  h.itabcoll.assign(kItabColl, Real(0));
  auto values = [&](int nskip, int nval, Real* out[]) {
    if (!std::getline(f, line)) throw std::runtime_error("P3: " + path + " ends early");
    const char* p = line.c_str();
    char* e = nullptr;
    for (int n = 0; n < nskip + nval; ++n) {
      const double v = std::strtod(p, &e);
      if (e == p) throw std::runtime_error("P3: unreadable line in " + path + ": " + line);
      if (n >= nskip) *out[n - nskip] = Real(v);
      p = e;
    }
  };
  auto ix = [](int jj, int ii, int i, int idx) {
    return (jj - 1) + densize * ((ii - 1) + rimsize * ((i - 1) + isize * (idx - 1)));
  };
  auto ixc = [](int jj, int ii, int i, int j, int idx) {
    return (jj - 1) + densize * ((ii - 1) + rimsize * ((i - 1) + isize * ((j - 1) + rcollsize * (idx - 1))));
  };
  for (int jj = 1; jj <= densize; ++jj)
    for (int ii = 1; ii <= rimsize; ++ii) {
      for (int i = 1; i <= isize; ++i) {
        Real* out[tabsize];
        for (int n = 0; n < tabsize; ++n) out[n] = &h.itab[std::size_t(ix(jj, ii, i, n + 1))];
        values(3, tabsize, out);
      }
      for (int i = 1; i <= isize; ++i)
        for (int j = 1; j <= rcollsize; ++j) {
          Real* out[2] = {&h.itabcoll[std::size_t(ixc(jj, ii, i, j, 1))], &h.itabcoll[std::size_t(ixc(jj, ii, i, j, 2))]};
          values(3, 2, out);
        }
    }

  // rain fall speeds and ventilation by numerical integration over the PSD
  h.vn.assign(kRainTab, Real(0));
  h.vm.assign(kRainTab, Real(0));
  h.revap.assign(kRainTab, Real(0));
  const Real mu_r = c.mu_r_constant;
  for (int jj = 1; jj <= 300; ++jj) {
    Real dm;
    if (jj <= 20) dm = (Real(jj) * Real(10.) - Real(5.)) * Real(1.e-6);
    else dm = (Real(jj - 20) * Real(30.) + Real(195.)) * Real(1.e-6);
    const Real lamr = (mu_r + Real(1)) / dm;
    Real dum1 = 0, dum2 = 0, dum3 = 0, dum4 = 0, dum5 = 0;
    const Real dd = Real(2.);
    for (int kk = 1; kk <= 10000; ++kk) {
      const Real dia = (Real(kk) * dd - dd * Real(0.5)) * Real(1.e-6);
      Real amg = piov6 * Real(997.) * ((dia * dia) * dia);
      amg = amg * Real(1000.);
      Real vt;
      if (dia * Real(1.e+6) <= Real(134.43)) vt = Real(4.5795e+3) * libm::pow(amg, Real(2.) * c.thrd);
      else if (dia * Real(1.e+6) < Real(1511.64)) vt = Real(4.962e+1) * libm::pow(amg, c.thrd);
      else if (dia * Real(1.e+6) < Real(3477.84)) vt = Real(1.732e+1) * libm::pow(amg, sxth);
      else vt = Real(9.17);
      const Real ld = libm::log10(dia);
      const Real ex = libm::exp(-(lamr * dia));
      dum1 = dum1 + vt * libm::pow(Real(10.), mu_r * ld + Real(4.) * mu_r) * ex * dd * Real(1.e-6);
      dum2 = dum2 + libm::pow(Real(10.), mu_r * ld + Real(4.) * mu_r) * ex * dd * Real(1.e-6);
      dum3 = dum3 + vt * libm::pow(Real(10.), (mu_r + Real(3.)) * ld + Real(4.) * mu_r) * ex * dd * Real(1.e-6);
      dum4 = dum4 + libm::pow(Real(10.), (mu_r + Real(3.)) * ld + Real(4.) * mu_r) * ex * dd * Real(1.e-6);
      dum5 = dum5 + libm::pow(vt * dia, Real(0.5)) * libm::pow(Real(10.), (mu_r + Real(1.)) * ld + Real(3.) * mu_r) *
                        ex * dd * Real(1.e-6);
    }
    dum2 = c_max(dum2, Real(1.e-30));
    dum4 = c_max(dum4, Real(1.e-30));
    dum5 = c_max(dum5, Real(1.e-30));
    for (int ii = 1; ii <= 10; ++ii) {
      const int n = P3Tables<Real>::r2(jj, ii);
      h.vn[std::size_t(n)] = dum1 / dum2;
      h.vm[std::size_t(n)] = dum3 / dum4;
      h.revap[std::size_t(n)] =
          libm::pow(Real(10.), libm::log10(dum5) + (mu_r + Real(1.)) * libm::log10(lamr) - (Real(3.) * mu_r));
    }
  }
  return h;
}

template <typename Real>
void run_p3(const P3Params<Real>& c, const P3Tables<Real>& tab, Real dt, int itimestep, bool predict_nc,
            const Range2D& range, int nz, View3<Real> th, View3<Real> qv, View3<Real> qc, View3<Real> qr,
            View3<Real> qnr, View3<Real> qi, View3<Real> qni, View3<Real> qir, View3<Real> qib, View3<Real> qnc,
            View3<Real> th_old, View3<Real> qv_old, View3<const Real> p, View3<const Real> dz, View3<Real> refl,
            View3<Real> re_cloud, View3<Real> re_ice, View3<Real> vmi, View3<Real> di, View3<Real> rhopo,
            View2<Real> rainnc, View2<Real> rainncv, View2<Real> sr, View2<Real> snownc, View2<Real> snowncv) {
  check_limit("P3", "nz", nz, kMaxLevels);
  auto launch = [&](auto nl) {
    P3Kernel<Real, decltype(nl)::value> k{c,      tab,    dt,       nz,     itimestep, predict_nc, range.i0,
                                          range.i1 - range.i0, th, qv,    qc,        qr,         qnr,
                                          qi,     qni,    qir,      qib,    qnc,       th_old,     qv_old,
                                          p,      dz,     refl,     re_cloud, re_ice,  vmi,        di,
                                          rhopo,  rainnc, rainncv,  sr,     snownc,    snowncv};
    column_for_work(range, k);
  };
  if (nz <= 48) launch(std::integral_constant<int, 48>{});
  else if (nz <= 64) launch(std::integral_constant<int, 64>{});
  else launch(std::integral_constant<int, kMaxLevels>{});
}

template <typename Real>
void p3_reserve(const Range2D& range, int nz) {
  const std::size_t b = nz <= 48   ? P3Kernel<Real, 48>::kWorkBytes
                        : nz <= 64 ? P3Kernel<Real, 64>::kWorkBytes
                                   : P3Kernel<Real, kMaxLevels>::kWorkBytes;
  column_work_reserve(range, b);
}

#define CPPWRF_P3_INST(R)                                                                                         \
  template P3Host<R> p3_init<R>(const std::string&);                                                              \
  template void p3_reserve<R>(const Range2D&, int);                                                               \
  template void run_p3<R>(const P3Params<R>&, const P3Tables<R>&, R, int, bool, const Range2D&, int, View3<R>,    \
                          View3<R>, View3<R>, View3<R>, View3<R>, View3<R>, View3<R>, View3<R>, View3<R>, View3<R>, \
                          View3<R>, View3<R>, View3<const R>, View3<const R>, View3<R>, View3<R>, View3<R>,        \
                          View3<R>, View3<R>, View3<R>, View2<R>, View2<R>, View2<R>, View2<R>, View2<R>);
CPPWRF_P3_INST(double)
CPPWRF_P3_INST(float)
#undef CPPWRF_P3_INST

}  // namespace cppwrf::phys::mp_p3
