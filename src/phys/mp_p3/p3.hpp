// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// WRF v4.8.0 phys/module_mp_p3.F (P3 v4.5.2): the Predicted Particle Properties microphysics with one free ice
// category and two-moment ice, mp_physics = 50 (P3_1CATEGORY: specified droplet number nccnst = 200 cm-3) and 51
// (P3_1CATEGORY_NC: predicted droplet number QNCLOUD), as WRF's mp_p3_wrapper_wrf runs P3_MAIN on (i, k) slabs.
// p3_init: the constants, lookupTable_1 (p3_lookupTable_1.dat-v5.4_2momI from the run directory) and the rain
// fall-speed / ventilation tables computed by numerical integration. Prognostic: theta, qv, qc, qr, qnr, qi, qni,
// qir (rime mass), qib (rime volume), qnc (51), TH_OLD / QV_OLD (the values after the previous call; zero before
// the first). Diagnostic: REFL_10CM (diag_ze), RE_CLOUD / RE_ICE (diag_effc / diag_effi, RRTMG's effective radii:
// has_reqc = has_reqi = 1, has_reqs = 0), V_ICE, D_ICE, RHO_ICE, the precipitation RAINNC, RAINNCV, SNOWNC,
// SNOWNCV, SR.
//
// WRF's settings in the WRF interface: SCPF off (cloud and precipitation fractions 1, the in-cloud vapour Qv_cld =
// the grid value of the first compute_SCPF call), log_predictSsat = .false., iparam = 3 (Khairoutdinov and Kogan
// 2000), nCat = 1 (log_hmossopOn = .false., no ice-ice category interactions, no icecat_destination),
// typeDiags_ON = .false., debug_on = .false., no visibility diagnostics, clbfact_dep = clbfact_sub = 1. Factors of
// those fractions (x*SCF, x*iSCF, x*(SPF-SPF_clr) = x*1) are left out: multiplying by 1 is exact.
//
// Arithmetic as gfortran -O2 -ftree-vectorize compiles the -r8 -fdefault-double-8 build (DOUBLE PRECISION is 8
// bytes: dexp(dble(x)) = exp(x)): x**2 = x*x, x**3 = (x*x)*x, x**6 = ((x*x)*x)*((x*x)*x) (the front end's
// addition chain), a REAL exponent a pow call (x**0.5 too), constant**constant folded correctly rounded
// ((25.e-6)**3), Fortran's left-to-right order. The only libmvec call is P3_MAIN's whole-slab
// tmparr1 = (pres*1.e-5)**(rd*inv_cp): the i loop of each k in pairs from its (cppwrf::vm::pow).
// C++ k = Fortran k - 1 (kts = 1, kbot = kts, ktop = kte, kdir = 1).
#include <cstddef>
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "core/column_work.hpp"
#include "core/field.hpp"
#include "core/libm.hpp"

namespace cppwrf::phys::mp_p3 {

constexpr int kMaxLevels = 128;
// lookupTable_1 dimensions (module parameters)
constexpr int isize = 50, densize = 5, rimsize = 4, rcollsize = 30, tabsize = 14, colltabsize = 2;
constexpr int kItab = densize * rimsize * isize * tabsize;
constexpr int kItabColl = densize * rimsize * isize * rcollsize * colltabsize;
constexpr int kRainTab = 300 * 10;

// p3_init's module variables that P3_MAIN reads
template <typename Real>
struct P3Params {
  Real pi, thrd, max_total_Ni, nccnst, cp, inv_cp, g, rd, rv, ep_2, rhosur, rhosui, f1r, f2r, rhow, cpw, inv_rhow,
      mu_r_constant, inv_Drmax, rho_rimeMin, rho_rimeMax, inv_rho_rimeMax, qsmall, nsmall, bsmall, aimm, mi0, eci,
      eri, bcn, nmltratio, e0, cons1, cons3, cons5, cons6, cons7, mw, rr, bact, inv_rm1, sig1, nanew1, inv_rm2, sig2,
      nanew2;
};

// the tables as P3_MAIN indexes them (Fortran column-major, 1-based indices in the accessors)
template <typename Real>
struct P3Tables {
  const Real* itab = nullptr;      // itab(densize, rimsize, isize, tabsize)
  const Real* itabcoll = nullptr;  // itabcoll(densize, rimsize, isize, rcollsize, colltabsize)
  const Real* vn = nullptr;        // vn_table(300, 10)
  const Real* vm = nullptr;        // vm_table(300, 10)
  const Real* revap = nullptr;     // revap_table(300, 10)
  CPPWRF_HD Real it(int jj, int ii, int i, int idx) const {
    return itab[(jj - 1) + densize * ((ii - 1) + rimsize * ((i - 1) + isize * (idx - 1)))];
  }
  CPPWRF_HD Real ic(int jj, int ii, int i, int j, int idx) const {
    return itabcoll[(jj - 1) + densize * ((ii - 1) + rimsize * ((i - 1) + isize * ((j - 1) + rcollsize * (idx - 1))))];
  }
  CPPWRF_HD static int r2(int jj, int ii) { return (jj - 1) + 300 * (ii - 1); }
};

// real function polysvp1(T, i_type) [Pa]: Flatau et al. (1992) at and above 195.8 K (ice) / 202.0 K (liquid),
// Goff-Gratch below; alog10(6.1071) and alog10(1013.246) are constants gfortran folds (correctly rounded)
template <typename Real>
CPPWRF_HD inline Real polysvp1(Real t, int i_type) {
  Real p = 0;
  if (i_type == 1 && t < Real(273.15)) {
    if (t >= Real(195.8)) {
      const Real dt = t - Real(273.15);
      p = Real(6.11147274) +
          dt * (Real(0.503160820) +
                dt * (Real(0.188439774e-1) +
                      dt * (Real(0.420895665e-3) +
                            dt * (Real(0.615021634e-5) +
                                  dt * (Real(0.602588177e-7) +
                                        dt * (Real(0.385852041e-9) +
                                              dt * (Real(0.146898966e-11) + Real(0.252751365e-14) * dt)))))));
      p = p * Real(100.);
    } else {
      p = libm::pow(Real(10.), -Real(9.09718) * (Real(273.16) / t - Real(1.)) -
                                   Real(3.56654) * libm::log10(Real(273.16) / t) +
                                   Real(0.876793) * (Real(1.) - t / Real(273.16)) + Real(0x1.9258f81f79246p-1)) *
          Real(100.);
    }
  } else if (i_type == 0 || t >= Real(273.15)) {
    if (t >= Real(202.0)) {
      const Real dt = t - Real(273.15);
      p = Real(6.11239921) +
          dt * (Real(0.443987641) +
                dt * (Real(0.142986287e-1) +
                      dt * (Real(0.264847430e-3) +
                            dt * (Real(0.302950461e-5) +
                                  dt * (Real(0.206739458e-7) +
                                        dt * (Real(0.640689451e-10) +
                                              dt * (Real(-0.952447341e-13) + Real(-0.976195544e-15) * dt)))))));
      p = p * Real(100.);
    } else {
      p = libm::pow(Real(10.),
                    -Real(7.90298) * (Real(373.16) / t - Real(1.)) + Real(5.02808) * libm::log10(Real(373.16) / t) -
                        Real(1.3816e-7) * (libm::pow(Real(10.), Real(11.344) * (Real(1.) - t / Real(373.16))) - Real(1.)) +
                        Real(8.1328e-3) *
                            (libm::pow(Real(10.), -Real(3.49149) * (Real(373.16) / t - Real(1.))) - Real(1.)) +
                        Real(0x1.80bb4409e541bp+1)) *
          Real(100.);
    }
  }
  return p;
}

// real function qv_sat(t_atm, p_atm, i_wrt): saturation mixing ratio over liquid (0) or ice (1)
template <typename Real>
CPPWRF_HD inline Real qv_sat(const P3Params<Real>& c, Real t, Real p, int i_wrt) {
  const Real e = polysvp1(t, i_wrt);
  return c.ep_2 * e / c_max(Real(1.e-3), p - e);
}

// real function DERF(X) of the module
template <typename Real>
CPPWRF_HD inline Real derf(Real x) {
  constexpr double A[65] = {
      0.00000000005958930743E0, -0.00000000113739022964E0, 0.00000001466005199839E0, -0.00000016350354461960E0,
      0.00000164610044809620E0, -0.00001492559551950604E0, 0.00012055331122299265E0, -0.00085483269811296660E0,
      0.00522397762482322257E0, -0.02686617064507733420E0, 0.11283791670954881569E0, -0.37612638903183748117E0,
      1.12837916709551257377E0, 0.00000000002372510631E0, -0.00000000045493253732E0, 0.00000000590362766598E0,
      -0.00000006642090827576E0, 0.00000067595634268133E0, -0.00000621188515924000E0, 0.00005103883009709690E0,
      -0.00037015410692956173E0, 0.00233307631218880978E0, -0.01254988477182192210E0, 0.05657061146827041994E0,
      -0.21379664776456006580E0, 0.84270079294971486929E0, 0.00000000000949905026E0, -0.00000000018310229805E0,
      0.00000000239463074000E0, -0.00000002721444369609E0, 0.00000028045522331686E0, -0.00000261830022482897E0,
      0.00002195455056768781E0, -0.00016358986921372656E0, 0.00107052153564110318E0, -0.00608284718113590151E0,
      0.02986978465246258244E0, -0.13055593046562267625E0, 0.67493323603965504676E0, 0.00000000000382722073E0,
      -0.00000000007421598602E0, 0.00000000097930574080E0, -0.00000001126008898854E0, 0.00000011775134830784E0,
      -0.00000111992758382650E0, 0.00000962023443095201E0, -0.00007404402135070773E0, 0.00050689993654144881E0,
      -0.00307553051439272889E0, 0.01668977892553165586E0, -0.08548534594781312114E0, 0.56909076642393639985E0,
      0.00000000000155296588E0, -0.00000000003032205868E0, 0.00000000040424830707E0, -0.00000000471135111493E0,
      0.00000005011915876293E0, -0.00000048722516178974E0, 0.00000430683284629395E0, -0.00003445026145385764E0,
      0.00024879276133931664E0, -0.00162940941748079288E0, 0.00988786373932350462E0, -0.05962426839442303805E0,
      0.49766113250947636708E0};
  constexpr double B[65] = {
      -0.00000000029734388465E0, 0.00000000269776334046E0, -0.00000000640788827665E0, -0.00000001667820132100E0,
      -0.00000021854388148686E0, 0.00000266246030457984E0, 0.00001612722157047886E0, -0.00025616361025506629E0,
      0.00015380842432375365E0, 0.00815533022524927908E0, -0.01402283663896319337E0, -0.19746892495383021487E0,
      0.71511720328842845913E0, -0.00000000001951073787E0, -0.00000000032302692214E0, 0.00000000522461866919E0,
      0.00000000342940918551E0, -0.00000035772874310272E0, 0.00000019999935792654E0, 0.00002687044575042908E0,
      -0.00011843240273775776E0, -0.00080991728956032271E0, 0.00661062970502241174E0, 0.00909530922354827295E0,
      -0.20160072778491013140E0, 0.51169696718727644908E0, 0.00000000003147682272E0, -0.00000000048465972408E0,
      0.00000000063675740242E0, 0.00000003377623323271E0, -0.00000015451139637086E0, -0.00000203340624738438E0,
      0.00001947204525295057E0, 0.00002854147231653228E0, -0.00101565063152200272E0, 0.00271187003520095655E0,
      0.02328095035422810727E0, -0.16725021123116877197E0, 0.32490054966649436974E0, 0.00000000002319363370E0,
      -0.00000000006303206648E0, -0.00000000264888267434E0, 0.00000002050708040581E0, 0.00000011371857327578E0,
      -0.00000211211337219663E0, 0.00000368797328322935E0, 0.00009823686253424796E0, -0.00065860243990455368E0,
      -0.00075285814895230877E0, 0.02585434424202960464E0, -0.11637092784486193258E0, 0.18267336775296612024E0,
      -0.00000000000367789363E0, 0.00000000020876046746E0, -0.00000000193319027226E0, -0.00000000435953392472E0,
      0.00000018006992266137E0, -0.00000078441223763969E0, -0.00000675407647949153E0, 0.00008428418334440096E0,
      -0.00017604388937031815E0, -0.00239729611435071610E0, 0.02064129023876022970E0, -0.06905562880005864105E0,
      0.09084526782065478489E0};
  const Real w = x < Real(0) ? -x : x;
  Real y;
  if (w < Real(2.2)) {
    Real t = w * w;
    int k = int(t);
    t = t - Real(k);
    k = k * 13;
    y = ((((((((((((Real(A[k]) * t + Real(A[k + 1])) * t + Real(A[k + 2])) * t + Real(A[k + 3])) * t +
                 Real(A[k + 4])) * t + Real(A[k + 5])) * t + Real(A[k + 6])) * t + Real(A[k + 7])) * t +
              Real(A[k + 8])) * t + Real(A[k + 9])) * t + Real(A[k + 10])) * t + Real(A[k + 11])) * t +
         Real(A[k + 12])) * w;
  } else if (w < Real(6.9)) {
    int k = int(w);
    const Real t = w - Real(k);
    k = 13 * (k - 2);
    y = (((((((((((Real(B[k]) * t + Real(B[k + 1])) * t + Real(B[k + 2])) * t + Real(B[k + 3])) * t +
                Real(B[k + 4])) * t + Real(B[k + 5])) * t + Real(B[k + 6])) * t + Real(B[k + 7])) * t +
             Real(B[k + 8])) * t + Real(B[k + 9])) * t + Real(B[k + 10])) * t + Real(B[k + 11])) * t +
        Real(B[k + 12]);
    y = y * y;
    y = y * y;
    y = y * y;
    y = Real(1) - y * y;
  } else {
    y = Real(1);
  }
  if (x < Real(0)) y = -y;
  return y;
}

// access_lookup_table: trilinear interpolation of itab(:,:,:,index)
template <typename Real>
CPPWRF_HD inline Real access_lookup_table(const P3Tables<Real>& T, int dumjj, int dumii, int dumi, int index,
                                          Real dum1, Real dum4, Real dum5) {
  Real iproc1 = T.it(dumjj, dumii, dumi, index) +
                (dum1 - Real(dumi)) * (T.it(dumjj, dumii, dumi + 1, index) - T.it(dumjj, dumii, dumi, index));
  Real gproc1 = T.it(dumjj, dumii + 1, dumi, index) +
                (dum1 - Real(dumi)) * (T.it(dumjj, dumii + 1, dumi + 1, index) - T.it(dumjj, dumii + 1, dumi, index));
  const Real tmp1 = iproc1 + (dum4 - Real(dumii)) * (gproc1 - iproc1);
  iproc1 = T.it(dumjj + 1, dumii, dumi, index) +
           (dum1 - Real(dumi)) * (T.it(dumjj + 1, dumii, dumi + 1, index) - T.it(dumjj + 1, dumii, dumi, index));
  gproc1 = T.it(dumjj + 1, dumii + 1, dumi, index) +
           (dum1 - Real(dumi)) *
               (T.it(dumjj + 1, dumii + 1, dumi + 1, index) - T.it(dumjj + 1, dumii + 1, dumi, index));
  const Real tmp2 = iproc1 + (dum4 - Real(dumii)) * (gproc1 - iproc1);
  return tmp1 + (dum5 - Real(dumjj)) * (tmp2 - tmp1);
}

// access_lookup_table_coll: ice-rain collection, interpolated in ice size, rain size, rime fraction and density
template <typename Real>
CPPWRF_HD inline Real access_lookup_table_coll(const P3Tables<Real>& T, int dumjj, int dumii, int dumj, int dumi,
                                               int index, Real dum1, Real dum3, Real dum4, Real dum5) {
  auto d = [&](int jj, int ii, int j) {
    return T.ic(jj, ii, dumi, j, index) +
           (dum1 - Real(dumi)) * (T.ic(jj, ii, dumi + 1, j, index) - T.ic(jj, ii, dumi, j, index));
  };
  Real dproc1 = d(dumjj, dumii, dumj), dproc2 = d(dumjj, dumii, dumj + 1);
  const Real iproc1 = dproc1 + (dum3 - Real(dumj)) * (dproc2 - dproc1);
  dproc1 = d(dumjj, dumii + 1, dumj);
  dproc2 = d(dumjj, dumii + 1, dumj + 1);
  Real gproc1 = dproc1 + (dum3 - Real(dumj)) * (dproc2 - dproc1);
  const Real tmp1 = iproc1 + (dum4 - Real(dumii)) * (gproc1 - iproc1);
  dproc1 = d(dumjj + 1, dumii, dumj);
  dproc2 = d(dumjj + 1, dumii, dumj + 1);
  const Real iproc2 = dproc1 + (dum3 - Real(dumj)) * (dproc2 - dproc1);
  dproc1 = d(dumjj + 1, dumii + 1, dumj);
  dproc2 = d(dumjj + 1, dumii + 1, dumj + 1);
  gproc1 = dproc1 + (dum3 - Real(dumj)) * (dproc2 - dproc1);
  const Real tmp2 = iproc2 + (dum4 - Real(dumii)) * (gproc1 - iproc2);
  return tmp1 + (dum5 - Real(dumjj)) * (tmp2 - tmp1);
}

// find_lookupTable_indices_1a: indices of normalized ice mass, rime fraction and bulk rime density
template <typename Real>
CPPWRF_HD inline void find_indices_1a(int& dumi, int& dumjj, int& dumii, Real& dum1, Real& dum4, Real& dum5,
                                      Real qitot, Real nitot, Real qirim, Real rhop) {
  dum1 = (libm::log10(qitot / nitot) + Real(18.)) * Real(3.444606) - Real(10.);
  dumi = int(dum1);
  dum1 = c_min(dum1, Real(isize));
  dum1 = c_max(dum1, Real(1.));
  dumi = dumi > 1 ? dumi : 1;
  dumi = isize - 1 < dumi ? isize - 1 : dumi;
  dum4 = (qirim / qitot) * Real(3.) + Real(1.);
  dumii = int(dum4);
  dum4 = c_min(dum4, Real(rimsize));
  dum4 = c_max(dum4, Real(1.));
  dumii = dumii > 1 ? dumii : 1;
  dumii = rimsize - 1 < dumii ? rimsize - 1 : dumii;
  if (rhop <= Real(650.)) dum5 = (rhop - Real(50.)) * Real(0.005) + Real(1.);
  else dum5 = (rhop - Real(650.)) * Real(0.004) + Real(4.);
  dumjj = int(dum5);
  dum5 = c_min(dum5, Real(densize));
  dum5 = c_max(dum5, Real(1.));
  dumjj = dumjj > 1 ? dumjj : 1;
  dumjj = densize - 1 < dumjj ? densize - 1 : dumjj;
}

// find_lookupTable_indices_1b: index of the scaled mean rain size
template <typename Real>
CPPWRF_HD inline void find_indices_1b(const P3Params<Real>& c, int& dumj, Real& dum3, Real qr, Real nr) {
  if (qr >= c.qsmall && nr > Real(0.)) {
    const Real dumlr = libm::pow(qr / (c.pi * c.rhow * nr), c.thrd);
    dum3 = (libm::log10(Real(1.) * dumlr) + Real(5.)) * Real(10.70415);
    dumj = int(dum3);
    dum3 = c_min(dum3, Real(rcollsize));
    dum3 = c_max(dum3, Real(1.));
    dumj = dumj > 1 ? dumj : 1;
    dumj = rcollsize - 1 < dumj ? rcollsize - 1 : dumj;
  } else {
    dumj = 1;
    dum3 = Real(1.);
  }
}

// find_lookupTable_indices_3: indices in the rain tables (scaled mean size, mu_r)
template <typename Real>
CPPWRF_HD inline void find_indices_3(const P3Params<Real>& c, int& dumii, int& dumjj, Real& rdumii, Real& rdumjj,
                                     Real mu_r, Real lamr) {
  const Real dum1 = (mu_r + Real(1.)) / lamr;
  if (dum1 <= Real(195.e-6)) {
    const Real inv_dum3 = Real(0.1);
    rdumii = (dum1 * Real(1.e6) + Real(5.)) * inv_dum3;
    rdumii = c_max(rdumii, Real(1.));
    rdumii = c_min(rdumii, Real(20.));
    dumii = int(rdumii);
    dumii = dumii > 1 ? dumii : 1;
    dumii = dumii < 20 ? dumii : 20;
  } else if (dum1 > Real(195.e-6)) {
    const Real inv_dum3 = c.thrd * Real(0.1);
    rdumii = (dum1 * Real(1.e+6) - Real(195.)) * inv_dum3 + Real(20.);
    rdumii = c_max(rdumii, Real(20.));
    rdumii = c_min(rdumii, Real(300.));
    dumii = int(rdumii);
    dumii = dumii > 20 ? dumii : 20;
    dumii = dumii < 299 ? dumii : 299;
  }
  rdumjj = mu_r + Real(1.);
  rdumjj = c_max(rdumjj, Real(1.));
  rdumjj = c_min(rdumjj, Real(10.));
  dumjj = int(rdumjj);
  dumjj = dumjj > 1 ? dumjj : 1;
  dumjj = dumjj < 9 ? dumjj : 9;
}

// get_cloud_dsd2 (iSCF = 1; iparam = 3: nu is not computed)
template <typename Real>
CPPWRF_HD inline void get_cloud_dsd2(const P3Params<Real>& c, Real qc, Real& nc_grd, Real& mu_c, Real rho, Real& nu,
                                     Real& lamc, Real& cdist, Real& cdist1) {
  if (qc >= c.qsmall) {
    Real nc = nc_grd;
    nc = c_max(nc, c.nsmall);
    mu_c = Real(0.0005714) * (nc * Real(1.e-6) * rho) + Real(0.2714);
    mu_c = Real(1.) / (mu_c * mu_c) - Real(1.);
    mu_c = c_max(mu_c, Real(2.));
    mu_c = c_min(mu_c, Real(15.));
    lamc = libm::pow(c.cons1 * nc * (mu_c + Real(3.)) * (mu_c + Real(2.)) * (mu_c + Real(1.)) / qc, c.thrd);
    const Real lammin = (mu_c + Real(1.)) * Real(2.5e+4);
    const Real lammax = (mu_c + Real(1.)) * Real(1.e+6);
    if (lamc < lammin) {
      lamc = lammin;
      nc = Real(6.) * (lamc * lamc * lamc) * qc /
           (c.pi * c.rhow * (mu_c + Real(3.)) * (mu_c + Real(2.)) * (mu_c + Real(1.)));
    } else if (lamc > lammax) {
      lamc = lammax;
      nc = Real(6.) * (lamc * lamc * lamc) * qc /
           (c.pi * c.rhow * (mu_c + Real(3.)) * (mu_c + Real(2.)) * (mu_c + Real(1.)));
    }
    cdist = nc * (mu_c + Real(1.)) / lamc;
    cdist1 = nc / libm::tgamma(mu_c + Real(1.));
    nc_grd = nc;
  } else {
    mu_c = 0;
    lamc = 0;
    cdist = 0;
    cdist1 = 0;
    nu = 0;
  }
}

// get_rain_dsd2 (iSPF = 1; mu_r = mu_r_constant)
template <typename Real>
CPPWRF_HD inline void get_rain_dsd2(const P3Params<Real>& c, Real qr, Real& nr_grd, Real& mu_r, Real& lamr,
                                    Real& cdistr, Real& logn0r) {
  if (qr >= c.qsmall) {
    Real nr = nr_grd;
    nr = c_max(nr, c.nsmall);
    mu_r = c.mu_r_constant;
    lamr = libm::pow(c.cons1 * nr * (mu_r + Real(3.)) * (mu_r + Real(2)) * (mu_r + Real(1.)) / (qr), c.thrd);
    const Real lammax = (mu_r + Real(1.)) * Real(1.e+5);
    const Real lammin = (mu_r + Real(1.)) * c.inv_Drmax;
    if (lamr < lammin) {
      lamr = lammin;
      nr = libm::exp(Real(3.) * libm::log(lamr) + libm::log(qr) + libm::log(libm::tgamma(mu_r + Real(1.))) -
                     libm::log(libm::tgamma(mu_r + Real(4.)))) /
           (c.cons1);
    } else if (lamr > lammax) {
      lamr = lammax;
      nr = libm::exp(Real(3.) * libm::log(lamr) + libm::log(qr) + libm::log(libm::tgamma(mu_r + Real(1.))) -
                     libm::log(libm::tgamma(mu_r + Real(4.)))) /
           (c.cons1);
    }
    logn0r = libm::log10(nr) + (mu_r + Real(1.)) * libm::log10(lamr) - libm::log10(libm::tgamma(mu_r + Real(1)));
    cdistr = nr / libm::tgamma(mu_r + Real(1.));
    nr_grd = nr;
  } else {
    lamr = 0;
    cdistr = 0;
    logn0r = 0;
  }
}

// calc_bulkRhoRime: bulk rime density, qirim and birim made consistent
template <typename Real>
CPPWRF_HD inline Real calc_bulk_rho_rime(const P3Params<Real>& c, Real qi_tot, Real& qi_rim, Real& bi_rim) {
  Real rho_rime;
  if (bi_rim >= Real(1.e-15)) {
    rho_rime = qi_rim / bi_rim;
    if (rho_rime < c.rho_rimeMin) {
      rho_rime = c.rho_rimeMin;
      bi_rim = qi_rim / rho_rime;
    } else if (rho_rime > c.rho_rimeMax) {
      rho_rime = c.rho_rimeMax;
      bi_rim = qi_rim / rho_rime;
    }
  } else {
    qi_rim = 0;
    bi_rim = 0;
    rho_rime = 0;
  }
  if (qi_rim > qi_tot && rho_rime > Real(0.)) {
    qi_rim = qi_tot;
    bi_rim = qi_rim / rho_rime;
  }
  if (qi_rim < c.qsmall) {
    qi_rim = 0;
    bi_rim = 0;
  }
  return rho_rime;
}

// impose_max_total_Ni (one category: sum(nitot) = 0 + nitot)
template <typename Real>
CPPWRF_HD inline void impose_max_total_ni(const P3Params<Real>& c, Real& nitot, Real inv_rho) {
  const Real s = Real(0) + nitot;
  if (s >= Real(1.e-20)) {
    const Real dum = c.max_total_Ni * inv_rho / s;
    nitot = nitot * c_min(dum, Real(1.));
  }
}

#define CPPWRF_P3_COLUMN(X)                                                                                          \
  X(Real, t, NL) X(Real, t_old, NL) X(Real, invexn, NL) X(Real, rho, NL) X(Real, inv_rho, NL) X(Real, rhofacr, NL)  \
  X(Real, rhofaci, NL) X(Real, acn, NL) X(Real, qvs, NL) X(Real, qvi, NL) X(Real, ssat, NL) X(Real, sup, NL)        \
  X(Real, supi, NL) X(Real, inv_dzq, NL) X(Real, mu_r, NL) X(Real, lamc, NL) X(Real, lamr, NL) X(Real, logn0r, NL) \
  X(Real, mu_c, NL) X(Real, nu, NL) X(Real, cdist, NL) X(Real, cdist1, NL) X(Real, cdistr, NL) X(Real, ze_ice, NL) \
  X(Real, ze_rain, NL) X(Real, qv_cld, NL) X(Real, nc, NL) X(Real, va, NL) X(Real, vb, NL) X(Real, fa, NL)          \
  X(Real, fb, NL) X(Real, fc, NL) X(Real, fd, NL)

// mp_p3_wrapper_wrf + P3_MAIN for one column (i, j), one ice category; i0, ni: the microphysics tile's i range
// (the lanes of the vectorised exner loop)
template <typename Real, int NL = kMaxLevels>
struct P3Kernel {
  P3Params<Real> c;
  P3Tables<Real> tab;
  Real dt;
  int nz, it;  // it: the time step counter (1 for the first step of a cold start)
  bool predict_nc;
  int i0, ni;
  View3<Real> th, qv, qc, qr, qnr, qi, qni, qir, qib, qnc, th_old, qv_old;
  View3<const Real> pres, dzq;
  View3<Real> diag_ze, diag_effc, diag_effi, diag_vmi, diag_di, diag_rhoi;
  View2<Real> rainnc, rainncv, sr, snownc, snowncv;
  ColumnWork work{};
  static constexpr std::size_t kWorkBytes = 0 CPPWRF_P3_COLUMN(CPPWRF_COLUMN_BYTES);

  CPPWRF_HD void operator()(int i, int j) const {
    if (nz > NL) return;
    CPPWRF_COLUMN_SLOTS(work, i, j)
    CPPWRF_P3_COLUMN(CPPWRF_COLUMN_ARRAY)
    const P3Params<Real>& C = c;
    const Real qsmall = C.qsmall, inv_cp = C.inv_cp;
    const Real odt = Real(1.) / dt;
    const int kbot = 0, ktop = nz - 1;
    // the wrapper: nc of the slab (predicted) or 0 (specified; P3_MAIN sets it from nccnst)
    for (int k = 0; k < nz; ++k) nc[k] = predict_nc ? qnc(i, j, k) : Real(0.);
    Real prt_liq = 0, prt_sol = 0;
    for (int k = 0; k < nz; ++k) {
      inv_dzq[k] = Real(1.) / dzq(i, j, k);
      mu_r[k] = 0;
      diag_ze(i, j, k) = Real(-99.);
      ze_ice[k] = Real(1.e-22);
      ze_rain[k] = Real(1.e-22);
      diag_effc(i, j, k) = Real(10.e-6);
      diag_effi(i, j, k) = Real(25.e-6);
      diag_vmi(i, j, k) = 0;
      diag_di(i, j, k) = 0;
      diag_rhoi(i, j, k) = 0;
      // tmparr1 = (pres*1.e-5)**(rd*inv_cp): libmvec's pow in pairs of the i loop
      const Real tmparr1 = vm::pow(pres(i, j, k) * Real(1.e-5), C.rd * inv_cp, i - i0, ni);
      invexn[k] = Real(1.) / tmparr1;
      t[k] = th(i, j, k) * tmparr1;
      t_old[k] = th_old(i, j, k) * tmparr1;
      qv(i, j, k) = c_max(qv(i, j, k), Real(0.));
    }
    Real rhorime_c = Real(400.);
    const Real xxlv = Real(3.1484e6) - Real(2370.) * Real(273.15);  // use constant Lv
    const Real xxls = xxlv + Real(0.3337e6);
    const Real xlf = xxls - xxlv;

    bool log_hydrometeorsPresent = false, log_nucleationPossible = false;
    for (int k = kbot; k <= ktop; ++k) {  // k_loop_1
      const Real p = pres(i, j, k);
      rho[k] = p / (C.rd * t[k]);
      inv_rho[k] = Real(1.) / rho[k];
      qvs[k] = qv_sat(C, c_max(t_old[k], Real(1.)), p, 0);
      qvi[k] = qv_sat(C, c_max(t_old[k], Real(1.)), p, 1);
      ssat[k] = qv_old(i, j, k) - qvs[k];
      sup[k] = qv_old(i, j, k) / qvs[k] - Real(1.);
      supi[k] = qv_old(i, j, k) / qvi[k] - Real(1.);
      rhofacr[k] = libm::pow(C.rhosur * inv_rho[k], Real(0.54));
      rhofaci[k] = libm::pow(C.rhosui * inv_rho[k], Real(0.54));
      const Real dum = Real(1.496e-6) * libm::pow(t[k], Real(1.5)) / (t[k] + Real(120.));
      acn[k] = C.g * C.rhow / (Real(18.) * dum);
      if (!predict_nc) nc[k] = C.nccnst * inv_rho[k];
      // .and. binds tighter than .or. (the WRF condition as written; SCPF_on = .false.)
      if ((t[k] < Real(273.15) && supi[k] >= Real(-0.05)) || (t[k] >= Real(273.15) && sup[k] >= Real(-0.05)))
        log_nucleationPossible = true;
      Real& QC = qc(i, j, k);
      Real& QR = qr(i, j, k);
      Real& QI = qi(i, j, k);
      Real& QV = qv(i, j, k);
      Real& TH = th(i, j, k);
      if (QC < qsmall || (QC < Real(1.e-8) && sup[k] < Real(-0.1))) {
        QV = QV + QC;
        TH = TH - invexn[k] * QC * xxlv * inv_cp;
        QC = 0;
        nc[k] = 0;
      } else {
        log_hydrometeorsPresent = true;
      }
      if (QR < qsmall || (QR < Real(1.e-8) && sup[k] < Real(-0.1))) {
        QV = QV + QR;
        TH = TH - invexn[k] * QR * xxlv * inv_cp;
        QR = 0;
        qnr(i, j, k) = 0;
      } else {
        log_hydrometeorsPresent = true;
      }
      if (QI < qsmall || (QI < Real(1.e-8) && supi[k] < Real(-0.1))) {
        QV = QV + QI;
        TH = TH - invexn[k] * QI * xxls * inv_cp;
        QI = 0;
        qni(i, j, k) = 0;
        qir(i, j, k) = 0;
        qib(i, j, k) = 0;
      } else {
        log_hydrometeorsPresent = true;
      }
      if (QI >= qsmall && QI < Real(1.e-8) && t[k] >= Real(273.15)) {
        QR = QR + QI;
        qnr(i, j, k) = qnr(i, j, k) + qni(i, j, k);
        TH = TH - invexn[k] * QI * xlf * inv_cp;
        QI = 0;
        qni(i, j, k) = 0;
        qir(i, j, k) = 0;
        qib(i, j, k) = 0;
      }
    }
    // first compute_SCPF (SCPF off): Qv_cld = Qv
    for (int k = 0; k < nz; ++k) qv_cld[k] = qv(i, j, k);

    if (log_nucleationPossible || log_hydrometeorsPresent) {
      log_hydrometeorsPresent = false;
      process_rates(i, j, log_hydrometeorsPresent, t, t_old, invexn, rho, inv_rho, rhofaci, acn, qvs, qvi, ssat,
                    sup, supi, mu_r, lamc, lamr, logn0r, mu_c, nu, cdist, cdist1, cdistr, qv_cld, nc, xxlv, xxls,
                    xlf, odt, rhorime_c);
      if (log_hydrometeorsPresent)
        sedimentation_and_final(i, j, prt_liq, prt_sol, t, invexn, rho, inv_rho, rhofacr, rhofaci, acn, inv_dzq,
                                mu_r, lamc, lamr, logn0r, mu_c, nu, cdistr, ze_ice, ze_rain, nc, va, vb,
                                fa, fb, fc, fd, xxlv, xxls, xlf, odt);
    }
    // 333: (WRF interface) th_old, qv_old for the next call; the precipitation of mp_p3_wrapper_wrf
    for (int k = 0; k < nz; ++k) {
      th_old(i, j, k) = th(i, j, k);
      qv_old(i, j, k) = qv(i, j, k);
      if (predict_nc) qnc(i, j, k) = nc[k];
    }
    const Real dum1 = Real(1000.) * dt;
    rainnc(i, j) = rainnc(i, j) + (prt_liq + prt_sol) * dum1;
    rainncv(i, j) = (prt_liq + prt_sol) * dum1;
    snownc(i, j) = snownc(i, j) + prt_sol * dum1;
    snowncv(i, j) = prt_sol * dum1;
    sr(i, j) = prt_sol / (prt_liq + prt_sol + Real(1.E-12));
  }

  // P3_MAIN's main k loop (process rates, conservation, updates)
  template <class A, class CA>
  CPPWRF_HD void process_rates(int i, int j, bool& log_hydrometeorsPresent, CA t, CA t_old, CA invexn, CA rho,
                               CA inv_rho, CA rhofaci, CA acn, CA qvs, CA qvi, CA ssat, CA sup, CA supi, A mu_r,
                               A lamc, A lamr, A logn0r, A mu_c, A nu, A cdist, A cdist1, A cdistr, CA qv_cld, A nc,
                               Real xxlv, Real xxls, Real xlf, Real odt, Real& rhorime_c) const {
    const P3Params<Real>& C = c;
    const Real qsmall = C.qsmall, nsmall = C.nsmall, pi = C.pi, thrd = C.thrd, inv_cp = C.inv_cp, cp = C.cp,
               rv = C.rv, g = C.g;
    Real eii = 0, Eii_fact = 0, epsi = 0, f1pr02 = 0, f1pr03 = 0, f1pr04 = 0, f1pr05 = 0, f1pr07 = 0, f1pr08 = 0,
         f1pr09 = 0, f1pr10 = 0, f1pr14 = 0;
    Real mu = 0, dv = 0, sc = 0, dqsdT = 0, dqsidT = 0, ab = 0, abi = 0, kap = 0, oabi = 0;
    for (int k = 0; k < nz; ++k) {  // k_loop_main
      Real& QC = qc(i, j, k);
      Real& QR = qr(i, j, k);
      Real& NR = qnr(i, j, k);
      Real& QI = qi(i, j, k);
      Real& NI = qni(i, j, k);
      Real& QIR = qir(i, j, k);
      Real& QIB = qib(i, j, k);
      Real& QV = qv(i, j, k);
      Real& TH = th(i, j, k);
      Real& NC = nc[k];
      const Real T = t[k], P = pres(i, j, k);
      bool log_exitlevel = !(QC >= qsmall || QR >= qsmall || QI >= qsmall);
      if (log_exitlevel && ((T < Real(273.15) && supi[k] < Real(-0.05)) || (T >= Real(273.15) && sup[k] < Real(-0.05))))
        continue;  // 555
      Real qcacc = 0, qrevp = 0, qccon = 0, qcaut = 0, qcevp = 0, qrcon = 0, ncacc = 0, ncnuc = 0, ncslf = 0,
           ncautc = 0, qcnuc = 0, nrslf = 0, nrevp = 0, ncautr = 0;
      Real qchetc = 0, qisub = 0, nrshdr = 0, qcheti = 0, qrcol = 0, qcshd = 0, qrhetc = 0, qimlt = 0, qccol = 0,
           qrheti = 0, qinuc = 0, nimlt = 0, nchetc = 0, nccol = 0, ncshdc = 0, ncheti = 0, nrcol = 0, nislf = 0,
           nrhetc = 0, ninuc = 0, qidep = 0, nrheti = 0, nisub = 0, qwgrth = 0, qrmul = 0, nimul = 0, qicol = 0,
           nicol = 0;
      bool log_wetgrowth = false;
      Real ssat_cld, ssat_r, sup_cld, sup_r, supi_cld, Q_nuc, N_nuc, dum, dum1, dum2, dumqvi, dumqvs, ratio, tmp1;
      log_exitlevel = !(QC >= qsmall || QR >= qsmall || QI >= qsmall);
      if (!log_exitlevel) {
        mu = Real(1.496e-6) * libm::pow(T, Real(1.5)) / (T + Real(120.));
        dv = Real(8.794e-5) * libm::pow(T, Real(1.81)) / P;
        sc = mu / (rho[k] * dv);
        dum = Real(1.) / (rv * (T * T));
        dqsdT = xxlv * qvs[k] * dum;
        dqsidT = xxls * qvi[k] * dum;
        ab = Real(1.) + dqsdT * xxlv * inv_cp;
        abi = Real(1.) + dqsidT * xxls * inv_cp;
        kap = Real(1.414e+3) * mu;
        if (T < Real(253.15)) eii = Real(0.001);
        else if (T >= Real(253.15) && T < Real(273.15))
          eii = Real(0.001) + (T - Real(253.15)) * (Real(0.3) - Real(0.001)) / Real(20.);
        else if (T >= Real(273.15)) eii = Real(0.3);
        get_cloud_dsd2(C, QC, NC, mu_c[k], rho[k], nu[k], lamc[k], cdist[k], cdist1[k]);
        get_rain_dsd2(C, QR, NR, mu_r[k], lamr[k], cdistr[k], logn0r[k]);
        Real epsi_tot = 0;
        impose_max_total_ni(C, NI, inv_rho[k]);
        if (QI >= qsmall) {  // qitot_notsmall_1
          NI = c_max(NI, nsmall);
          NR = c_max(NR, nsmall);
          const Real rhop = calc_bulk_rho_rime(C, QI, QIR, QIB);
          int dumi, dumjj, dumii, dumj;
          Real d1, d3, d4, d5;
          find_indices_1a(dumi, dumjj, dumii, d1, d4, d5, QI, NI, QIR, rhop);
          find_indices_1b(C, dumj, d3, QR, NR);
          f1pr02 = access_lookup_table(tab, dumjj, dumii, dumi, 2, d1, d4, d5);
          f1pr03 = access_lookup_table(tab, dumjj, dumii, dumi, 3, d1, d4, d5);
          f1pr04 = access_lookup_table(tab, dumjj, dumii, dumi, 4, d1, d4, d5);
          f1pr05 = access_lookup_table(tab, dumjj, dumii, dumi, 5, d1, d4, d5);
          f1pr09 = access_lookup_table(tab, dumjj, dumii, dumi, 7, d1, d4, d5);
          f1pr10 = access_lookup_table(tab, dumjj, dumii, dumi, 8, d1, d4, d5);
          f1pr14 = access_lookup_table(tab, dumjj, dumii, dumi, 10, d1, d4, d5);
          if (QR >= qsmall) {
            f1pr07 = access_lookup_table_coll(tab, dumjj, dumii, dumj, dumi, 1, d1, d3, d4, d5);
            f1pr08 = access_lookup_table_coll(tab, dumjj, dumii, dumj, dumi, 2, d1, d3, d4, d5);
          } else {
            f1pr07 = Real(-99.);
            f1pr08 = Real(-99.);
          }
          NI = c_min(NI, f1pr09 * QI);
          NI = c_max(NI, f1pr10 * QI);
          if (QIR > Real(0.)) {
            tmp1 = QIR / QI;
            if (tmp1 < Real(0.6)) Eii_fact = Real(1.);
            else if (tmp1 >= Real(0.6) && tmp1 < Real(0.9)) Eii_fact = Real(1.) - (tmp1 - Real(0.6)) / Real(0.3);
            else if (tmp1 >= Real(0.9)) Eii_fact = Real(0.);
          } else {
            Eii_fact = Real(1.);
          }
        }
        // collection of droplets
        if (QI >= qsmall && QC >= qsmall && T <= Real(273.15)) {
          qccol = rhofaci[k] * f1pr04 * QC * C.eci * rho[k] * NI;
          nccol = rhofaci[k] * f1pr04 * NC * C.eci * rho[k] * NI;
        }
        if (QI >= qsmall && QC >= qsmall && T > Real(273.15)) {
          qcshd = rhofaci[k] * f1pr04 * QC * C.eci * rho[k] * NI;
          nccol = rhofaci[k] * f1pr04 * NC * C.eci * rho[k] * NI;
          ncshdc = qcshd * Real(1.923e+6);
        }
        // collection of rain (f1pr07, f1pr08 and logn0r are log10 values)
        if (QI >= qsmall && QR >= qsmall && T <= Real(273.15)) {
          qrcol = libm::pow(Real(10.), f1pr08 + logn0r[k]) * rho[k] * rhofaci[k] * C.eri * NI;
          nrcol = libm::pow(Real(10.), f1pr07 + logn0r[k]) * rho[k] * rhofaci[k] * C.eri * NI;
        }
        if (QI >= qsmall && QR >= qsmall && T > Real(273.15)) {
          nrcol = libm::pow(Real(10.), f1pr07 + logn0r[k]) * rho[k] * rhofaci[k] * C.eri * NI;
        }
        // self-collection of ice
        if (QI >= qsmall) nislf = f1pr03 * rho[k] * eii * Eii_fact * rhofaci[k] * NI * NI;
        // melting (enhanced melting from collision neglected: dum = 0)
        if (QI >= qsmall && T > Real(273.15)) {
          const Real qsat0 = Real(0.622) * C.e0 / (P - C.e0);
          dum = Real(0.);
          qimlt = ((f1pr05 + f1pr14 * libm::pow(sc, thrd) * libm::pow(rhofaci[k] * rho[k] / mu, Real(0.5))) *
                       ((T - Real(273.15)) * kap - rho[k] * xxlv * dv * (qsat0 - qv_cld[k])) * Real(2.) * pi / xlf +
                   dum) *
                  NI;
          qimlt = c_max(qimlt, Real(0.));
          nimlt = qimlt * (NI / QI);
        }
        // wet growth
        if (QI >= qsmall && (QC + QR) >= Real(1.e-6) && T < Real(273.15)) {
          const Real qsat0 = Real(0.622) * C.e0 / (P - C.e0);
          qwgrth = ((f1pr05 + f1pr14 * libm::pow(sc, thrd) * libm::pow(rhofaci[k] * rho[k] / mu, Real(0.5))) *
                    Real(2.) * pi * (rho[k] * xxlv * dv * (qsat0 - qv_cld[k]) - (T - Real(273.15)) * kap) /
                    (xlf + C.cpw * (T - Real(273.15)))) *
                   NI;
          qwgrth = c_max(qwgrth, Real(0.));
          dum = c_max(Real(0.), (qccol + qrcol) - qwgrth);
          if (dum >= Real(1.e-10)) {
            nrshdr = nrshdr + dum * Real(1.923e+6);
            if ((qccol + qrcol) >= Real(1.e-10)) {
              dum1 = Real(1.) / (qccol + qrcol);
              qcshd = qcshd + dum * qccol * dum1;
              qccol = qccol - dum * qccol * dum1;
              qrcol = qrcol - dum * qrcol * dum1;
            }
            log_wetgrowth = true;
          }
        }
        // inverse ice relaxation timescale
        if (QI >= qsmall && T < Real(273.15)) {
          epsi = ((f1pr05 + f1pr14 * libm::pow(sc, thrd) * libm::pow(rhofaci[k] * rho[k] / mu, Real(0.5))) *
                  Real(2.) * pi * rho[k] * dv) *
                 NI;
          epsi_tot = epsi_tot + epsi;
        } else {
          epsi = Real(0.);
        }
        // rime density (Cober and List 1993)
        if (qccol >= qsmall && T < Real(273.15)) {
          const Real vtrmi1 = f1pr02 * rhofaci[k];
          const Real iTc = Real(1.) / c_min(Real(-0.001), T - Real(273.15));
          if (QC >= qsmall) {
            const Real Vt_qc = acn[k] * libm::tgamma(Real(4.) + C.bcn + mu_c[k]) /
                               (libm::pow(lamc[k], C.bcn) * libm::tgamma(mu_c[k] + Real(4.)));
            const Real D_c = (mu_c[k] + Real(4.)) / lamc[k];
            const Real dv_i = vtrmi1 - Vt_qc;
            const Real V_impact = dv_i < Real(0) ? -dv_i : dv_i;
            Real Ri = -(Real(0.5e+6) * D_c) * V_impact * iTc;
            Ri = c_max(Real(1.), c_min(Ri, Real(12.)));
            if (Ri <= Real(8.)) rhorime_c = (Real(0.051) + Real(0.114) * Ri - Real(0.0055) * (Ri * Ri)) * Real(1000.);
            else rhorime_c = Real(611.) + Real(72.25) * (Ri - Real(8.));
          }
        } else {
          rhorime_c = Real(400.);
        }
        // contact and immersion freezing of droplets
        if (QC >= qsmall && T <= Real(269.15)) {
          const Real d = Real(1.) / lamc[k];
          dum = (d * d) * d;
          tmp1 = cdist1[k] * libm::exp(C.aimm * (Real(273.15) - T));
          Q_nuc = C.cons6 * libm::tgamma(Real(7.) + mu_c[k]) * tmp1 * (dum * dum);
          N_nuc = C.cons5 * libm::tgamma(mu_c[k] + Real(4.)) * tmp1 * dum;
          qcheti = Q_nuc;
          ncheti = N_nuc;
        }
        // immersion freezing of rain
        if (QR >= qsmall && T <= Real(269.15)) {
          const Real tmpdbl1 = libm::exp(libm::log(cdistr[k]) + libm::log(libm::tgamma(Real(7.) + mu_r[k])) -
                                         Real(6.) * libm::log(lamr[k]));
          const Real tmpdbl2 = libm::exp(libm::log(cdistr[k]) + libm::log(libm::tgamma(mu_r[k] + Real(4.))) -
                                         Real(3.) * libm::log(lamr[k]));
          const Real tmpdbl3 = libm::exp(C.aimm * (Real(273.15) - T));
          Q_nuc = C.cons6 * (tmpdbl1 * tmpdbl3);
          N_nuc = C.cons5 * (tmpdbl2 * tmpdbl3);
          qrheti = Q_nuc;
          nrheti = N_nuc;
        }
        // condensation / evaporation, deposition / sublimation: rain evaporation with ventilation
        Real epsr, epsc, xx;
        if (QR >= qsmall) {
          int dumii, dumjj;
          Real rdumii, rdumjj;
          find_indices_3(C, dumii, dumjj, rdumii, rdumjj, mu_r[k], lamr[k]);
          const Real* rv_t = tab.revap;
          dum1 = rv_t[tab.r2(dumii, dumjj)] +
                 (rdumii - Real(dumii)) * (rv_t[tab.r2(dumii + 1, dumjj)] - rv_t[tab.r2(dumii, dumjj)]);
          dum2 = rv_t[tab.r2(dumii, dumjj + 1)] +
                 (rdumii - Real(dumii)) * (rv_t[tab.r2(dumii + 1, dumjj + 1)] - rv_t[tab.r2(dumii, dumjj + 1)]);
          dum = dum1 + (rdumjj - Real(dumjj)) * (dum2 - dum1);
          epsr = Real(2.) * pi * cdistr[k] * rho[k] * dv *
                 (C.f1r * libm::tgamma(mu_r[k] + Real(2.)) / (lamr[k]) +
                  C.f2r * libm::pow(rho[k] / mu, Real(0.5)) * libm::pow(sc, thrd) * dum);
        } else {
          epsr = Real(0.);
        }
        if (QC >= qsmall) epsc = Real(2.) * pi * rho[k] * dv * cdist[k];
        else epsc = Real(0.);
        if (T < Real(273.15)) {
          oabi = Real(1.) / abi;
          xx = epsc + epsr + epsi_tot * (Real(1.) + xxls * inv_cp * dqsdT) * oabi;
        } else {
          xx = epsc + epsr;
        }
        dumqvi = qvi[k];
        dum = -cp / g * (T - t_old[k]) * odt;
        Real aaa;
        if (T < Real(273.15))
          aaa = (QV - qv_old(i, j, k)) * odt - dqsdT * (-dum * g * inv_cp) -
                (qvs[k] - dumqvi) * (Real(1.) + xxls * inv_cp * dqsdT) * oabi * epsi_tot;
        else
          aaa = (QV - qv_old(i, j, k)) * odt - dqsdT * (-dum * g * inv_cp);
        xx = c_max(Real(1.e-20), xx);
        const Real oxx = Real(1.) / xx;
        ssat_cld = ssat[k];
        ssat_r = ssat[k];
        sup_cld = sup[k];
        sup_r = sup[k];
        supi_cld = supi[k];
        if (QC >= qsmall)
          qccon = (aaa * epsc * oxx + (ssat_cld - aaa * oxx) * odt * epsc * oxx * (Real(1.) - libm::exp(-(xx * dt)))) / ab;
        if (QR >= qsmall)
          qrcon = (aaa * epsr * oxx + (ssat_r - aaa * oxx) * odt * epsr * oxx * (Real(1.) - libm::exp(-(xx * dt)))) / ab;
        if (sup_cld < Real(-0.001) && QC < Real(1.e-12)) qccon = -QC * odt;
        if (sup_r < Real(-0.001) && QR < Real(1.e-12)) qrcon = -QR * odt;
        if (qccon < Real(0.)) {
          qcevp = -qccon;
          qccon = Real(0.);
        } else {
          qccon = c_min(qccon, QV * odt);
        }
        if (qrcon < Real(0.)) {
          qrevp = -qrcon;
          nrevp = qrevp * (NR / QR);
          qrcon = Real(0.);
        } else {
          qrcon = c_min(qrcon, QV * odt);
        }
        if (QI >= qsmall && T < Real(273.15))
          qidep = (aaa * epsi * oxx + (ssat_cld - aaa * oxx) * odt * epsi * oxx * (Real(1.) - libm::exp(-(xx * dt)))) *
                      oabi +
                  (qvs[k] - dumqvi) * epsi * oabi;
        if (supi_cld < Real(-0.001) && QI < Real(1.e-12)) qidep = -QI * odt;
        if (qidep < Real(0.)) {
          qisub = -qidep;
          qisub = c_min(qisub, QI * odt);
          nisub = qisub * (NI / QI);
          qidep = Real(0.);
        } else {
          qidep = c_min(qidep, QV * odt);
        }
      }
      // 444: deposition / condensation-freezing nucleation
      sup_cld = sup[k];
      supi_cld = supi[k];
      if (T < Real(258.15) && supi_cld >= Real(0.05)) {
        dum = Real(0.005) * libm::exp(Real(0.304) * (Real(273.15) - T)) * Real(1000.) * inv_rho[k];  // Cooper (1986)
        dum = c_min(dum, Real(100.e3) * inv_rho[k]);
        N_nuc = c_max(Real(0.), (dum - (Real(0) + NI)) * odt);
        if (N_nuc >= Real(1.e-20)) {
          Q_nuc = c_max(Real(0.), (dum - (Real(0) + NI)) * C.mi0 * odt);
          qinuc = Q_nuc;
          ninuc = N_nuc;
        }
      }
      // droplet activation
      if (!predict_nc && sup_cld > Real(1.e-6) && it > 1) {
        dum = C.nccnst * inv_rho[k] * C.cons7 - QC;
        dum = c_max(Real(0.), dum);
        dumqvs = qv_sat(C, T, P, 0);
        dqsdT = xxlv * dumqvs / (rv * T * T);
        ab = Real(1.) + dqsdT * xxlv * inv_cp;
        dum = c_max(Real(0.), c_min(dum, (qv_cld[k] - dumqvs) / ab));
        qcnuc = dum * odt;
      }
      if (predict_nc) {
        if (sup_cld > Real(1.e-6)) {
          dum1 = Real(1.) / libm::pow(C.bact, Real(0.5));
          const Real sigvl = Real(0.0761) - Real(1.55e-4) * (T - Real(273.15));
          const Real aact = Real(2.) * C.mw / (C.rhow * C.rr * T) * sigvl;
          const Real sm1 = Real(2.) * dum1 * libm::pow(aact * thrd * C.inv_rm1, Real(1.5));
          const Real sm2 = Real(2.) * dum1 * libm::pow(aact * thrd * C.inv_rm2, Real(1.5));
          const Real uu1 = Real(2.) * libm::log(sm1 / sup_cld) / (Real(4.242) * libm::log(C.sig1));
          const Real uu2 = Real(2.) * libm::log(sm2 / sup_cld) / (Real(4.242) * libm::log(C.sig2));
          dum1 = C.nanew1 * Real(0.5) * (Real(1.) - derf(uu1));
          dum2 = C.nanew2 * Real(0.5) * (Real(1.) - derf(uu2));
          dum2 = c_min((C.nanew1 + C.nanew2), dum1 + dum2);
          dum2 = (dum2 - NC) * odt;
          dum2 = c_max(Real(0.), dum2);
          ncnuc = dum2;
          if (it <= 1) qcnuc = Real(0.);
          else qcnuc = ncnuc * C.cons7;
        }
      }
      // saturation adjustment at the first step
      if (it <= 1) {
        const Real dumt = TH * libm::pow(P * Real(1.e-5), C.rd * inv_cp);
        const Real dumqv = qv_cld[k];
        dumqvs = qv_sat(C, dumt, P, 0);
        const Real dums = dumqv - dumqvs;
        qccon = dums / (Real(1.) + xxlv * xxlv * dumqvs / (cp * rv * (dumt * dumt))) * odt;
        qccon = c_max(Real(0.), qccon);
        if (qccon <= Real(1.e-7)) qccon = Real(0.);
      }
      // autoconversion (iparam = 3: Khairoutdinov and Kogan 2000)
      if (QC >= Real(1.e-8)) {
        dum = QC;
        qcaut = Real(1350.) * libm::pow(dum, Real(2.47)) * libm::pow(NC * Real(1.e-6) * rho[k], Real(-1.79));
        ncautr = qcaut * C.cons3;
        ncautc = qcaut * NC / QC;
        if (qcaut == Real(0.)) ncautc = Real(0.);
        if (ncautc == Real(0.)) qcaut = Real(0.);
      }
      // self-collection of droplets (iparam = 3: none)
      if (QC >= qsmall) ncslf = Real(0.);
      // accretion of cloud by rain (iparam = 3)
      if (QR >= qsmall && QC >= qsmall) {
        qcacc = Real(67.) * libm::pow(QC * QR, Real(1.15));
        ncacc = qcacc * NC / QC;
        if (qcacc == Real(0.)) ncacc = Real(0.);
        if (ncacc == Real(0.)) qcacc = Real(0.);
      }
      // self-collection and breakup of rain
      if (QR >= qsmall) {
        dum1 = Real(280.e-6);
        NR = c_max(NR, nsmall);
        dum2 = libm::pow(QR / (pi * C.rhow * NR), thrd);
        if (dum2 < dum1) dum = Real(1.);
        else if (dum2 >= dum1) dum = Real(2.) - libm::exp(Real(2300.) * (dum2 - dum1));
        nrslf = dum * Real(5.78) * NR * QR * rho[k];
      }
      // conservation: condensation and evaporation limited to saturation adjustment
      dumqvs = qv_sat(C, T, P, 0);
      const Real qcon_satadj = (qv_cld[k] - dumqvs) / (Real(1.) + xxlv * xxlv * dumqvs / (cp * rv * (T * T))) * odt;
      tmp1 = qccon + qrcon + qcnuc;
      if (tmp1 > Real(0.) && qcon_satadj < Real(0.)) {
        qccon = Real(0.);
        qrcon = Real(0.);
        qcnuc = Real(0.);
        ncnuc = Real(0.);
      } else {
        if (tmp1 > Real(0.) && tmp1 > qcon_satadj) {
          ratio = c_max(Real(0.), qcon_satadj) / tmp1;
          ratio = c_min(Real(1.), ratio);
          qccon = qccon * ratio;
          qrcon = qrcon * ratio;
          qcnuc = qcnuc * ratio;
          ncnuc = ncnuc * ratio;
        } else if (qcevp + qrevp > Real(0.)) {
          ratio = c_max(Real(0.), -qcon_satadj) / (qcevp + qrevp);
          ratio = c_min(Real(1.), ratio);
          qcevp = qcevp * ratio;
          qrevp = qrevp * ratio;
          nrevp = nrevp * ratio;
        }
      }
      // deposition (with nucleation) and sublimation limited to saturation adjustment
      const Real qv_tmp = qv_cld[k] + (-qcnuc - qccon - qrcon + qcevp + qrevp) * dt;
      const Real t_tmp = T + (qcnuc + qccon + qrcon - qcevp - qrevp) * xxlv * inv_cp * dt;
      dumqvi = qv_sat(C, t_tmp, P, 1);
      const Real qdep_satadj =
          (qv_tmp - dumqvi) / (Real(1.) + xxls * xxls * dumqvi / (cp * rv * (t_tmp * t_tmp))) * odt;
      tmp1 = (Real(0) + qidep) + (Real(0) + qinuc);
      if (tmp1 > Real(0.) && qdep_satadj < Real(0.)) {
        qidep = Real(0.);
        qinuc = Real(0.);
        ninuc = Real(0.);
      } else {
        if (tmp1 > Real(0.) && tmp1 > qdep_satadj) {
          ratio = c_max(Real(0.), qdep_satadj) / tmp1;
          ratio = c_min(Real(1.), ratio);
          qidep = qidep * ratio;
          qinuc = qinuc * ratio;
          ninuc = ninuc * ratio;
        }
        dum = c_max(qisub, Real(1.e-20));
        qisub = qisub * c_min(Real(1.), c_max(Real(0.), -qdep_satadj) / c_max(Real(0) + qisub, Real(1.e-20)));
        nisub = nisub * c_min(Real(1.), qisub / dum);
      }
      // cloud
      Real sinks = (qcaut + qcacc + (Real(0) + qccol) + qcevp + (Real(0) + qchetc) + (Real(0) + qcheti) +
                    (Real(0) + qcshd)) *
                   dt;
      Real sources = QC + (qccon + qcnuc) * dt;
      if (sinks > sources && sinks >= Real(1.e-20)) {
        ratio = sources / sinks;
        qcaut = qcaut * ratio;
        qcacc = qcacc * ratio;
        qcevp = qcevp * ratio;
        qccol = qccol * ratio;
        qcheti = qcheti * ratio;
        qcshd = qcshd * ratio;
        ncautc = ncautc * ratio;
        ncacc = ncacc * ratio;
        nccol = nccol * ratio;
        ncheti = ncheti * ratio;
      }
      // rain
      sinks = (qrevp + (Real(0) + qrcol) + (Real(0) + qrhetc) + (Real(0) + qrheti) + (Real(0) + qrmul)) * dt;
      sources = QR + (qrcon + qcaut + qcacc + (Real(0) + qimlt) + (Real(0) + qcshd)) * dt;
      if (sinks > sources && sinks >= Real(1.e-20)) {
        ratio = sources / sinks;
        qrevp = qrevp * ratio;
        qrcol = qrcol * ratio;
        qrheti = qrheti * ratio;
        qrmul = qrmul * ratio;
        nrevp = nrevp * ratio;
        nrcol = nrcol * ratio;
        nrheti = nrheti * ratio;
      }
      // ice (one category: qicol(1,1) = 0)
      sinks = (qisub + qimlt) * dt;
      sources = QI + (qidep + qinuc + qrcol + qccol + qrhetc + qrheti + qchetc + qcheti + qrmul) * dt;
      sources = sources + qicol * dt;
      sinks = sinks + qicol * dt;
      if (sinks > sources && sinks >= Real(1.e-20)) {
        ratio = sources / sinks;
        qisub = qisub * ratio;
        qimlt = qimlt * ratio;
        nisub = nisub * ratio;
        nimlt = nimlt * ratio;
        qicol = qicol * ratio;
        nicol = nicol * ratio;
      }
      // vapor
      sinks = (qccon + qrcon + qcnuc + (Real(0) + qidep) + (Real(0) + qinuc)) * dt;
      sources = QV + (qcevp + qrevp + (Real(0) + qisub)) * dt;
      if (sinks > sources && sinks >= Real(1.e-20)) {
        ratio = sources / sinks;
        qccon = qccon * ratio;
        qrcon = qrcon * ratio;
        qcnuc = qcnuc * ratio;
        qidep = qidep * ratio;
        qinuc = qinuc * ratio;
        ninuc = ninuc * ratio;
        ncnuc = ncnuc * ratio;
      }
      // update of the prognostic variables: ice-phase processes
      Real rimevolume = 0, rimefraction = 0;
      if (QI >= qsmall) {
        tmp1 = Real(1.) / QI;
        rimevolume = QIB * tmp1;
        rimefraction = QIR * tmp1;
      }
      QC = QC + (-qchetc - qcheti - qccol - qcshd) * dt;
      if (predict_nc) NC = NC + (-nccol - nchetc - ncheti) * dt;
      QR = QR + (-qrcol + qimlt - qrhetc - qrheti + qcshd - qrmul) * dt;
      NR = NR + (-nrcol - nrhetc - nrheti + C.nmltratio * nimlt + nrshdr + ncshdc) * dt;
      QIB = QIB - ((qisub + qimlt) * dt * rimevolume);
      QIR = QIR - ((qisub + qimlt) * dt * rimefraction);
      QI = QI - (qisub + qimlt) * dt;
      dum = (qrcol + qccol + qrhetc + qrheti + qchetc + qcheti + qrmul) * dt;
      QI = QI + (qidep + qinuc) * dt + dum;
      QIR = QIR + dum;
      QIB = QIB + (qrcol * C.inv_rho_rimeMax + qccol / rhorime_c +
                   (qrhetc + qrheti + qchetc + qcheti + qrmul) * C.inv_rho_rimeMax) *
                      dt;
      NI = NI + (ninuc - nimlt - nisub - nislf + nrhetc + nrheti + nchetc + ncheti + nimul) * dt;
      if (QIR < Real(0.)) {
        QIR = Real(0.);
        QIB = Real(0.);
      }
      if (log_wetgrowth) {  // densify ice during wet growth
        QIR = QI;
        QIB = QIR * C.inv_rho_rimeMax;
      }
      if (QI >= qsmall && QIB >= C.bsmall && qimlt > Real(0.)) {  // densify rimed ice during melting
        tmp1 = QIR / QIB;
        const Real tmp2 = QI + qimlt * dt;
        QIB = QIR / (tmp1 + (Real(917.) - tmp1) * qimlt * dt / tmp2);
      }
      QV = QV + (-qidep + qisub - qinuc) * dt;
      TH = TH + invexn[k] *
                    ((qidep - qisub + qinuc) * xxls * inv_cp +
                     (qrcol + qccol + qchetc + qcheti + qrhetc + qrheti + qrmul - qimlt) * xlf * inv_cp) *
                    dt;
      // warm-phase processes
      QC = QC + (-qcacc - qcaut + qcnuc + qccon - qcevp) * dt;
      QR = QR + (qcacc + qcaut + qrcon - qrevp) * dt;
      if (predict_nc) NC = NC + (-ncacc - ncautc + ncslf + ncnuc) * dt;
      else NC = C.nccnst * inv_rho[k];
      NR = NR + (ncautr - nrslf - nrevp) * dt;
      QV = QV + (-qcnuc - qccon - qrcon + qcevp + qrevp) * dt;
      TH = TH + invexn[k] * ((qcnuc + qccon + qrcon - qcevp - qrevp) * xxlv * inv_cp) * dt;
      // clipping
      if (QC < qsmall) {
        QV = QV + QC;
        TH = TH - invexn[k] * QC * xxlv * inv_cp;
        QC = Real(0.);
        NC = Real(0.);
      } else {
        log_hydrometeorsPresent = true;
      }
      if (QR < qsmall) {
        QV = QV + QR;
        TH = TH - invexn[k] * QR * xxlv * inv_cp;
        QR = Real(0.);
        NR = Real(0.);
      } else {
        log_hydrometeorsPresent = true;
      }
      if (QI < qsmall) {
        QV = QV + QI;
        TH = TH - invexn[k] * QI * xxls * inv_cp;
        QI = Real(0.);
        NI = Real(0.);
        QIR = Real(0.);
        QIB = Real(0.);
      } else {
        log_hydrometeorsPresent = true;
      }
      QV = c_max(Real(0.), QV);
      impose_max_total_ni(C, NI, inv_rho[k]);
      (void)ssat_r;
    }
  }

  // sedimentation (cloud, rain, ice; adaptive sub-steps), homogeneous freezing, final checks and diagnostics
  template <class A, class CA>
  CPPWRF_HD void sedimentation_and_final(int i, int j, Real& prt_liq, Real& prt_sol, CA t, CA invexn, CA rho,
                                         CA inv_rho, CA rhofacr, CA rhofaci, CA acn, CA inv_dzq, A mu_r, A lamc, A lamr,
                                         A logn0r, A mu_c, A nu, A cdistr, A ze_ice, A ze_rain, A nc,
                                         A va, A vb, A fa, A fb, A fc, A fd, Real xxlv, Real xxls, Real xlf,
                                         Real odt) const {
    const P3Params<Real>& C = c;
    const Real qsmall = C.qsmall, nsmall = C.nsmall, inv_cp = C.inv_cp;
    const int kbot = 0, ktop = nz - 1;
    Real tmp1, tmp2;
    // shared sub-stepping: fall speeds by `speeds`, then the fluxes of up to four quantities
    auto substeps = [&](int k_qxtop, int k_qxbot, Real& prt_accum, auto&& speeds, auto&& fluxes, auto&& update) {
      Real dt_left = dt;
      while (dt_left > Real(1.e-4)) {
        Real Co_max = Real(0.);
        for (int k = 0; k < nz; ++k) { va[k] = Real(0.); vb[k] = Real(0.); }
        for (int k = k_qxtop; k >= k_qxbot; --k) {
          speeds(k);
          Co_max = c_max(Co_max, va[k] * dt_left * inv_dzq[k]);
        }
        const int tmpint1 = int(Co_max + Real(1.));
        const Real dt_sub = c_min(dt_left, dt_left / Real(tmpint1));
        const int k_temp = k_qxbot == kbot ? k_qxbot : k_qxbot - 1;
        for (int k = k_temp; k <= k_qxtop; ++k) fluxes(k);
        if (k_qxbot == kbot) prt_accum = prt_accum + fa[kbot] * dt_sub;
        update(k_qxtop, true, dt_sub);
        for (int k = k_qxtop - 1; k >= k_temp; --k) update(k, false, dt_sub);
        dt_left = dt_left - dt_sub;
        if (k_qxbot != kbot) k_qxbot = k_qxbot - 1;
      }
    };
    // cloud
    {
      bool present = false;
      int k_qxtop = kbot, k_qxbot = kbot;
      for (int k = ktop; k >= kbot; --k)
        if (qc(i, j, k) >= qsmall) { present = true; k_qxtop = k; break; }
      if (present) {
        Real prt_accum = Real(0.);
        for (int k = kbot; k <= k_qxtop; ++k)
          if (qc(i, j, k) >= qsmall) { k_qxbot = k; break; }
        substeps(
            k_qxtop, k_qxbot, prt_accum,
            [&](int k) {
              if (qc(i, j, k) > qsmall) {
                Real t1, t2;
                get_cloud_dsd2(C, qc(i, j, k), nc[k], mu_c[k], rho[k], nu[k], lamc[k], t1, t2);
                const Real dum = Real(1.) / libm::pow(lamc[k], C.bcn);
                va[k] = acn[k] * libm::tgamma(Real(4.) + C.bcn + mu_c[k]) * dum / (libm::tgamma(mu_c[k] + Real(4.)));
                if (predict_nc)
                  vb[k] = acn[k] * libm::tgamma(Real(1.) + C.bcn + mu_c[k]) * dum / (libm::tgamma(mu_c[k] + Real(1.)));
              }
            },
            [&](int k) {
              fa[k] = va[k] * qc(i, j, k) * rho[k];
              if (predict_nc) fb[k] = vb[k] * nc[k] * rho[k];
            },
            [&](int k, bool top, Real dt_sub) {
              const Real dq = top ? -fa[k] * inv_dzq[k] : (fa[k + 1] - fa[k]) * inv_dzq[k];
              qc(i, j, k) = qc(i, j, k) + dq * dt_sub * inv_rho[k];
              if (predict_nc) {
                const Real dn = top ? -fb[k] * inv_dzq[k] : (fb[k + 1] - fb[k]) * inv_dzq[k];
                nc[k] = nc[k] + dn * dt_sub * inv_rho[k];
              }
            });
        prt_liq = prt_accum * C.inv_rhow * odt;
      }
    }
    // rain
    {
      bool present = false;
      int k_qxtop = kbot, k_qxbot = kbot;
      for (int k = ktop; k >= kbot; --k)
        if (qr(i, j, k) >= qsmall) { present = true; k_qxtop = k; break; }
      if (present) {
        Real prt_accum = Real(0.);
        for (int k = kbot; k <= k_qxtop; ++k)
          if (qr(i, j, k) >= qsmall) { k_qxbot = k; break; }
        substeps(
            k_qxtop, k_qxbot, prt_accum,
            [&](int k) {
              if (qr(i, j, k) > qsmall) {
                qnr(i, j, k) = c_max(qnr(i, j, k), nsmall);
                get_rain_dsd2(C, qr(i, j, k), qnr(i, j, k), mu_r[k], lamr[k], cdistr[k], logn0r[k]);
                int dumii, dumjj;
                Real rdumii, rdumjj;
                find_indices_3(C, dumii, dumjj, rdumii, rdumjj, mu_r[k], lamr[k]);
                const Real* vm = tab.vm;
                const Real* vn = tab.vn;
                Real dum1 = vm[tab.r2(dumii, dumjj)] +
                            (rdumii - Real(dumii)) * (vm[tab.r2(dumii + 1, dumjj)] - vm[tab.r2(dumii, dumjj)]);
                Real dum2 = vm[tab.r2(dumii, dumjj + 1)] +
                            (rdumii - Real(dumii)) * (vm[tab.r2(dumii + 1, dumjj + 1)] - vm[tab.r2(dumii, dumjj + 1)]);
                va[k] = dum1 + (rdumjj - Real(dumjj)) * (dum2 - dum1);
                va[k] = va[k] * rhofacr[k];
                dum1 = vn[tab.r2(dumii, dumjj)] +
                       (rdumii - Real(dumii)) * (vn[tab.r2(dumii + 1, dumjj)] - vn[tab.r2(dumii, dumjj)]);
                dum2 = vn[tab.r2(dumii, dumjj + 1)] +
                       (rdumii - Real(dumii)) * (vn[tab.r2(dumii + 1, dumjj + 1)] - vn[tab.r2(dumii, dumjj + 1)]);
                vb[k] = dum1 + (rdumjj - Real(dumjj)) * (dum2 - dum1);
                vb[k] = vb[k] * rhofacr[k];
              }
            },
            [&](int k) {
              fa[k] = va[k] * qr(i, j, k) * rho[k];
              fb[k] = vb[k] * qnr(i, j, k) * rho[k];
            },
            [&](int k, bool top, Real dt_sub) {
              const Real dq = top ? -fa[k] * inv_dzq[k] : (fa[k + 1] - fa[k]) * inv_dzq[k];
              const Real dn = top ? -fb[k] * inv_dzq[k] : (fb[k + 1] - fb[k]) * inv_dzq[k];
              qr(i, j, k) = qr(i, j, k) + dq * dt_sub * inv_rho[k];
              qnr(i, j, k) = qnr(i, j, k) + dn * dt_sub * inv_rho[k];
            });
        prt_liq = prt_liq + prt_accum * C.inv_rhow * odt;
      }
    }
    // ice (two-moment)
    {
      bool present = false;
      int k_qxtop = kbot, k_qxbot = kbot;
      for (int k = ktop; k >= kbot; --k)
        if (qi(i, j, k) >= qsmall) { present = true; k_qxtop = k; break; }
      if (present) {
        Real prt_accum = Real(0.);
        for (int k = kbot; k <= k_qxtop; ++k)
          if (qi(i, j, k) >= qsmall) { k_qxbot = k; break; }
        substeps(
            k_qxtop, k_qxbot, prt_accum,
            [&](int k) {
              if (qi(i, j, k) >= qsmall) {
                qni(i, j, k) = c_max(qni(i, j, k), nsmall);
                const Real rhop = calc_bulk_rho_rime(C, qi(i, j, k), qir(i, j, k), qib(i, j, k));
                int dumi, dumjj, dumii;
                Real d1, d4, d5;
                find_indices_1a(dumi, dumjj, dumii, d1, d4, d5, qi(i, j, k), qni(i, j, k), qir(i, j, k), rhop);
                const Real f1pr01 = access_lookup_table(tab, dumjj, dumii, dumi, 1, d1, d4, d5);
                const Real f1pr02 = access_lookup_table(tab, dumjj, dumii, dumi, 2, d1, d4, d5);
                const Real f1pr09 = access_lookup_table(tab, dumjj, dumii, dumi, 7, d1, d4, d5);
                const Real f1pr10 = access_lookup_table(tab, dumjj, dumii, dumi, 8, d1, d4, d5);
                qni(i, j, k) = c_min(qni(i, j, k), f1pr09 * qi(i, j, k));
                qni(i, j, k) = c_max(qni(i, j, k), f1pr10 * qi(i, j, k));
                va[k] = f1pr02 * rhofaci[k];
                vb[k] = f1pr01 * rhofaci[k];
              }
            },
            [&](int k) {
              fa[k] = va[k] * qi(i, j, k) * rho[k];
              fb[k] = vb[k] * qni(i, j, k) * rho[k];
              fc[k] = va[k] * qir(i, j, k) * rho[k];
              fd[k] = va[k] * qib(i, j, k) * rho[k];
            },
            [&](int k, bool top, Real dt_sub) {
              const Real dq = top ? -fa[k] * inv_dzq[k] : (fa[k + 1] - fa[k]) * inv_dzq[k];
              const Real dr = top ? -fc[k] * inv_dzq[k] : (fc[k + 1] - fc[k]) * inv_dzq[k];
              const Real db = top ? -fd[k] * inv_dzq[k] : (fd[k + 1] - fd[k]) * inv_dzq[k];
              const Real dn = top ? -fb[k] * inv_dzq[k] : (fb[k + 1] - fb[k]) * inv_dzq[k];
              qi(i, j, k) = qi(i, j, k) + dq * dt_sub * inv_rho[k];
              qir(i, j, k) = qir(i, j, k) + dr * dt_sub * inv_rho[k];
              qib(i, j, k) = qib(i, j, k) + db * dt_sub * inv_rho[k];
              qni(i, j, k) = qni(i, j, k) + dn * dt_sub * inv_rho[k];
            });
        prt_sol = prt_sol + prt_accum * C.inv_rhow * odt;
      }
    }
    // homogeneous freezing of cloud and rain (t: the temperature at the start of the scheme)
    for (int k = kbot; k <= ktop; ++k) {
      Real& TH = th(i, j, k);
      if (qc(i, j, k) >= qsmall && t[k] < Real(233.15)) {
        const Real Q_nuc = qc(i, j, k);
        nc[k] = c_max(nc[k], nsmall);
        const Real N_nuc = nc[k];
        qir(i, j, k) = qir(i, j, k) + Q_nuc;
        qi(i, j, k) = qi(i, j, k) + Q_nuc;
        qib(i, j, k) = qib(i, j, k) + Q_nuc * C.inv_rho_rimeMax;
        qni(i, j, k) = qni(i, j, k) + N_nuc;
        TH = TH + invexn[k] * Q_nuc * xlf * inv_cp;
        qc(i, j, k) = Real(0.);
        nc[k] = Real(0.);
      }
      if (qr(i, j, k) >= qsmall && t[k] < Real(233.15)) {
        const Real Q_nuc = qr(i, j, k);
        qnr(i, j, k) = c_max(qnr(i, j, k), nsmall);
        const Real N_nuc = qnr(i, j, k);
        qir(i, j, k) = qir(i, j, k) + Q_nuc;
        qi(i, j, k) = qi(i, j, k) + Q_nuc;
        qib(i, j, k) = qib(i, j, k) + Q_nuc * C.inv_rho_rimeMax;
        qni(i, j, k) = qni(i, j, k) + N_nuc;
        TH = TH + invexn[k] * Q_nuc * xlf * inv_cp;
        qr(i, j, k) = Real(0.);
        qnr(i, j, k) = Real(0.);
      }
    }
    // final checks and diagnostics
    for (int k = kbot; k <= ktop; ++k) {
      Real& QC = qc(i, j, k);
      Real& QR = qr(i, j, k);
      Real& NR = qnr(i, j, k);
      Real& QI = qi(i, j, k);
      Real& NI = qni(i, j, k);
      Real& QIR = qir(i, j, k);
      Real& QIB = qib(i, j, k);
      Real& QV = qv(i, j, k);
      Real& TH = th(i, j, k);
      if (QC >= qsmall) {
        get_cloud_dsd2(C, QC, nc[k], mu_c[k], rho[k], nu[k], lamc[k], tmp1, tmp2);
        diag_effc(i, j, k) = Real(0.5) * (mu_c[k] + Real(3.)) / lamc[k];
      } else {
        QV = QV + QC;
        TH = TH - invexn[k] * QC * xxlv * inv_cp;
        QC = Real(0.);
        nc[k] = Real(0.);
      }
      if (QR >= qsmall) {
        get_rain_dsd2(C, QR, NR, mu_r[k], lamr[k], tmp1, tmp2);
        const Real l = lamr[k];
        const Real l3 = (l * l) * l;
        ze_rain[k] = rho[k] * NR * (mu_r[k] + Real(6.)) * (mu_r[k] + Real(5.)) * (mu_r[k] + Real(4.)) *
                     (mu_r[k] + Real(3.)) * (mu_r[k] + Real(2.)) * (mu_r[k] + Real(1.)) / (l3 * l3);
        ze_rain[k] = c_max(ze_rain[k], Real(1.e-22));
      } else {
        QV = QV + QR;
        TH = TH - invexn[k] * QR * xxlv * inv_cp;
        QR = Real(0.);
        NR = Real(0.);
      }
      impose_max_total_ni(C, NI, inv_rho[k]);
      if (QI >= qsmall) {
        NI = c_max(NI, nsmall);
        NR = c_max(NR, nsmall);
        const Real rhop = calc_bulk_rho_rime(C, QI, QIR, QIB);
        int dumi, dumjj, dumii;
        Real d1, d4, d5;
        find_indices_1a(dumi, dumjj, dumii, d1, d4, d5, QI, NI, QIR, rhop);
        const Real f1pr02 = access_lookup_table(tab, dumjj, dumii, dumi, 2, d1, d4, d5);
        const Real f1pr06 = access_lookup_table(tab, dumjj, dumii, dumi, 6, d1, d4, d5);
        const Real f1pr09 = access_lookup_table(tab, dumjj, dumii, dumi, 7, d1, d4, d5);
        const Real f1pr10 = access_lookup_table(tab, dumjj, dumii, dumi, 8, d1, d4, d5);
        const Real f1pr13 = access_lookup_table(tab, dumjj, dumii, dumi, 9, d1, d4, d5);
        const Real f1pr15 = access_lookup_table(tab, dumjj, dumii, dumi, 11, d1, d4, d5);
        const Real f1pr16 = access_lookup_table(tab, dumjj, dumii, dumi, 12, d1, d4, d5);
        NI = c_min(NI, f1pr09 * QI);
        NI = c_max(NI, f1pr10 * QI);
        if (QIR < qsmall) {
          QIR = Real(0.);
          QIB = Real(0.);
        }
        diag_vmi(i, j, k) = f1pr02 * rhofaci[k];
        diag_effi(i, j, k) = f1pr06;
        diag_di(i, j, k) = f1pr15;
        diag_rhoi(i, j, k) = f1pr16;
        ze_ice[k] = ze_ice[k] + Real(0.1892) * f1pr13 * NI * rho[k];
        ze_ice[k] = c_max(ze_ice[k], Real(1.e-22));
      } else {
        QV = QV + QI;
        TH = TH - invexn[k] * QI * xxls * inv_cp;
        QI = Real(0.);
        NI = Real(0.);
        QIR = Real(0.);
        QIB = Real(0.);
        diag_di(i, j, k) = Real(0.);
      }
      diag_ze(i, j, k) = Real(10.) * libm::log10((ze_rain[k] + ze_ice[k]) * Real(1.e+18));
      if (QR < qsmall) NR = Real(0.);
    }
  }
};

// p3_init: the constants, lookupTable_1 and the rain tables on the host
template <typename Real>
struct P3Host {
  P3Params<Real> c{};
  std::vector<Real> itab, itabcoll, vn, vm, revap;
};
template <typename Real>
P3Host<Real> p3_init(const std::string& dir);

// the GPU work buffer of run_p3 over range (outside the time loop; CPU build: no-op)
template <typename Real>
void p3_reserve(const Range2D& range, int nz);

// mp_p3_wrapper_wrf over the microphysics columns `range` (i0 = range.i0 is WRF's its)
template <typename Real>
void run_p3(const P3Params<Real>& c, const P3Tables<Real>& tab, Real dt, int itimestep, bool predict_nc,
            const Range2D& range, int nz, View3<Real> th, View3<Real> qv, View3<Real> qc, View3<Real> qr,
            View3<Real> qnr, View3<Real> qi, View3<Real> qni, View3<Real> qir, View3<Real> qib, View3<Real> qnc,
            View3<Real> th_old, View3<Real> qv_old, View3<const Real> p, View3<const Real> dz, View3<Real> refl,
            View3<Real> re_cloud, View3<Real> re_ice, View3<Real> vmi, View3<Real> di, View3<Real> rhopo,
            View2<Real> rainnc, View2<Real> rainncv, View2<Real> sr, View2<Real> snownc, View2<Real> snowncv);

}  // namespace cppwrf::phys::mp_p3
