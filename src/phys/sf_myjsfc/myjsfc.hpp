#pragma once
// Port of WRF 4.8.0 phys/module_sf_myjsfc.F (myjsfc, sfcdif)
// Mellor-Yamada-Janjic (Eta) surface layer scheme (sf_sfclay_physics = 2).

#include <cmath>
#include <vector>
#include "core/backend.hpp"
#include "core/field.hpp"
#include "core/precision.hpp"
#include "phys/driver/wrf_constants.hpp"

namespace cppwrf::phys::sf_myjsfc {

// Table size (module_sf_myjsfc.F:55)
static constexpr int KZTM = 10001;
static constexpr int KZTM2 = KZTM - 2;
static constexpr int kMaxLevels = 128;

template <typename T>
CPPWRF_HD inline T c_abs(T a) { return (a < T(0)) ? -a : a; }

// Lookup tables and bounds initialized by MYJSFCINIT (module_sf_myjsfc.F:57-59, 1177-1300)
struct MYJSFCTables {
  View1<const real> psim1;
  View1<const real> psih1;
  View1<const real> psim2;
  View1<const real> psih2;
  real ztmin1;
  real ztmax1;
  real dzeta1;
  real fh01;
  real ztmin2;
  real ztmax2;
  real dzeta2;
  real fh02;
};

// Host helper to allocate and initialize MYJSFC lookup tables on device/host
class MYJSFCTableManager {
 public:
  MYJSFCTableManager();
  ~MYJSFCTableManager();
  MYJSFCTableManager(const MYJSFCTableManager&) = delete;
  MYJSFCTableManager& operator=(const MYJSFCTableManager&) = delete;

  MYJSFCTables get_tables() const { return tables_; }

 private:
  void* dev_psim1_ = nullptr;
  void* dev_psih1_ = nullptr;
  void* dev_psim2_ = nullptr;
  void* dev_psih2_ = nullptr;
  MYJSFCTables tables_{};
};

// C++ column kernel for MYJSFC (module_sf_myjsfc.F:66-1056)
struct MYJSFCKernel {
  // 2D/3D Atmospheric and Surface Inputs
  View2<const real> ht;
  View3<const real> dz;
  View3<const real> pmid;
  View3<const real> pint;
  View3<const real> th;
  View3<const real> t;
  View3<const real> qv;
  View3<const real> qc;
  View3<const real> u;
  View3<const real> v;
  View3<const real> q2;       // TKE = 0.5 * q^2
  View2<const real> tsk;
  View2<const real> xland;
  View2<const int> ivgtyp;
  View2<const real> z0base;
  View2<const real> mavail;
  View2<const int> lowlyr;

  // Inouts
  View2<real> qsfc;
  View2<real> thz0;
  View2<real> qz0;
  View2<real> uz0;
  View2<real> vz0;
  View2<real> ustar;
  View2<real> znt;
  View2<real> pblh;
  View2<real> rmol;
  View2<real> akhs;
  View2<real> akms;

  // Outputs
  View2<real> rib;
  View2<real> chs;
  View2<real> chs2;
  View2<real> cqs2;
  View2<real> hfx;
  View2<real> qfx;
  View2<real> flx_lh;
  View2<real> flhc;
  View2<real> flqc;
  View2<real> qgh;
  View2<real> cpm;
  View2<real> ct;
  View2<real> u10;
  View2<real> v10;
  View2<real> t02;
  View2<real> th02;
  View2<real> tshltr;
  View2<real> th10;
  View2<real> q02;
  View2<real> qshltr;
  View2<real> q10;
  View2<real> pshltr;
  View2<real> u10e;
  View2<real> v10e;

  // Tables and parameters
  MYJSFCTables tables;
  int itimestep;
  int isurban;
  int iz0tlnd;
  real p1000mb;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    // Physical constants from WRF module_model_constants.F and module_sf_myjsfc.F:25-53
    constexpr real G = wrfc::g;
    constexpr real R_D = wrfc::r_d;
    constexpr real CP = wrfc::cp;
    constexpr real R_V = wrfc::r_v;
    constexpr real XLV = wrfc::xlv;
    constexpr real VKARMAN = real(0.4);
    constexpr real CAPA = R_D / CP;
    constexpr real ELOCP = real(2.72e6) / CP;
    constexpr real RCAP = real(1.) / CAPA;
    constexpr real GOCP02 = G / CP * real(2.);
    constexpr real GOCP10 = G / CP * real(10.);
    constexpr real EPSU2 = real(1.e-6);
    constexpr real EPSUST = real(1.e-9);
    constexpr real EPSZT = real(1.e-28);
    constexpr real A2S = real(17.2693882);
    constexpr real A3S = real(273.16);
    constexpr real A4S = real(35.86);
    constexpr real PQ0 = real(379.90516);
    constexpr real SEAFC = real(0.98);
    constexpr real PQ0SEA = PQ0 * SEAFC;
    constexpr real EXCML = real(0.0001);
    constexpr real EXCMS = real(0.0001);
    constexpr real GLKBR = real(10.);
    constexpr real GLKBS = real(30.);
    constexpr real QVISC = real(2.1e-5);
    constexpr real RIC = real(0.505);
    constexpr real SMALL = real(0.35);
    constexpr real SQPR = real(0.84);
    constexpr real SQSC = real(0.84);
    constexpr real SQVISC = real(258.2);
    constexpr real TVISC = real(2.1e-5);
    constexpr real USTC = real(0.7);
    constexpr real USTR = real(0.225);
    constexpr real VISC = real(1.5e-5);
    constexpr real FH = real(1.01);
    constexpr real WWST = real(1.2);
    constexpr real ZTFC = real(1.);
    constexpr real P608 = R_V / R_D - real(1.);
    constexpr real EPSQ2 = real(0.2);

    constexpr real CZIV = SMALL * GLKBS;
    constexpr real GRRS = GLKBR / GLKBS;
    constexpr real RTVISC = real(1.) / TVISC;
    constexpr real RVISC = real(1.) / VISC;
    constexpr real ZQRZT = SQSC / SQPR;
    constexpr real FZQ1 = RTVISC * QVISC * ZQRZT;
    constexpr real FZQ2 = RTVISC * QVISC * ZQRZT;
    constexpr real FZT1 = RVISC * TVISC * SQPR;
    constexpr real FZT2 = CZIV * GRRS * TVISC * SQPR;
    constexpr real FZU1 = CZIV * VISC;
    constexpr real RQVISC = real(1.) / QVISC;
    constexpr real USTFC = real(0.018) / G;
    constexpr real WWST2 = WWST * WWST;
    constexpr int ITRMX = 5;

    const real ZTMIN1 = tables.ztmin1;
    const real ZTMAX1 = tables.ztmax1;
    const real DZETA1 = tables.dzeta1;
    const real FH01 = tables.fh01;

    const real ZTMIN2 = tables.ztmin2;
    const real ZTMAX2 = tables.ztmax2;
    const real DZETA2 = tables.dzeta2;
    const real FH02 = tables.fh02;

    const int ntsd = itimestep;

    // module_sf_myjsfc.F:191-210
    if (ntsd == 1) {
      ustar(i, j) = real(0.1);
    }
    ct(i, j) = real(0.);

    // Surface layer index (1-based LOWLYR in WRF)
    int llow_f = lowlyr.p ? lowlyr(i, j) : 1;
    int lmh = nz - llow_f + 1;
    int k_sfc = llow_f - 1;

    real psfc_val = pint(i, j, k_sfc);
    real apesfc = rpow(psfc_val / p1000mb, CAPA);
    real thsk_val = tsk(i, j) / apesfc;
    real seamask = xland(i, j) - real(1.);

    // Layer interface heights ZHK (module_sf_myjsfc.F:165-184, 250-253)
    real zhk[kMaxLevels + 2];  // 1..nz+1
    zhk[nz + 1] = ht(i, j);
    for (int k_f = nz; k_f >= 1; --k_f) {
      int k_wrf = nz - k_f;
      zhk[k_f] = zhk[k_f + 1] + dz(i, j, k_wrf);
    }

    // Find PBL height (module_sf_myjsfc.F:263-277)
    int lpbl = 1;
    for (int k_f = lmh - 1; k_f >= 1; --k_f) {
      int k_wrf = nz - k_f;
      real q2_k = real(2.) * q2(i, j, k_wrf);
      if (q2_k <= EPSQ2 * FH) {
        lpbl = k_f;
        break;
      }
    }
    real pblh_val = zhk[lpbl] - zhk[lmh + 1];
    pblh(i, j) = pblh_val;

    // Atmospheric state at lowest layer (module_sf_myjsfc.F:284-305)
    real plow = pmid(i, j, k_sfc);
    real tlow = t(i, j, k_sfc);
    real thlow = th(i, j, k_sfc);
    real ratiomx = qv(i, j, k_sfc);
    real qlow = ratiomx / (real(1.) + ratiomx);
    real cwmlow = qc(i, j, k_sfc);
    real thelow = (cwmlow * (-ELOCP / tlow) + real(1.)) * thlow;
    real ulow = u(i, j, k_sfc);
    real vlow = v(i, j, k_sfc);
    real zsl = (zhk[lmh] - zhk[lmh + 1]) * real(0.5);

    real tz0 = (ntsd == 1) ? tsk(i, j) : thz0(i, j) * apesfc;

    // ------------------------------------------------------------------
    // SFCDIF computation (module_sf_myjsfc.F:361-1056)
    // ------------------------------------------------------------------
    real rdz = real(1.) / zsl;
    real cxchl = EXCML * rdz;
    real cxchs = EXCMS * rdz;

    real btgx = G / thlow;
    real elfc = VKARMAN * btgx;

    real btgh = (pblh_val > real(1000.)) ? (btgx * pblh_val) : (btgx * real(1000.));

    real ustar_val = ustar(i, j);
    real z0_val = znt(i, j);
    real z0base_val = z0base(i, j);
    real uz0_val = uz0(i, j);
    real vz0_val = vz0(i, j);
    real thz0_val = thz0(i, j);
    real qz0_val = qz0(i, j);
    real qs_val = qsfc(i, j);
    real akms_val = akms(i, j);
    real akhs_val = akhs(i, j);
    real wetm = mavail(i, j);

    real zu = real(0.);
    real zt = real(0.);
    real zq = real(0.);
    real rib_val = real(0.);
    real rlmo_val = real(0.);
    real ustark = ustar_val * VKARMAN;
    real wstar2 = real(0.);
    real psmz = real(0.);
    real pshz = real(0.);

    if (seamask > real(0.5)) {
      // SEA POINTS (module_sf_myjsfc.F:443-636)
      for (int itr = 1; itr <= ITRMX; ++itr) {
        z0_val = c_max(USTFC * ustar_val * ustar_val, real(1.59e-5));
        if (ustar_val < USTC) {
          if (ustar_val < USTR) {
            if (ntsd == 1) {
              akms_val = cxchs;
              akhs_val = cxchs;
              qs_val = qlow;
            }
            zu = FZU1 * std::sqrt(std::sqrt(z0_val * ustar_val * RVISC)) / ustar_val;
            real wght = akms_val * zu * RVISC;
            real rwgh = wght / (wght + real(1.));
            uz0_val = (ulow * rwgh + uz0_val) * real(0.5);
            vz0_val = (vlow * rwgh + vz0_val) * real(0.5);

            zt = FZT1 * zu;
            zq = FZQ1 * zt;
            real wghtt = akhs_val * zt * RTVISC;
            real wghtq = akhs_val * zq * RQVISC;

            if (ntsd > 1) {
              thz0_val = ((wghtt * thlow + thsk_val) / (wghtt + real(1.)) + thz0_val) * real(0.5);
              qz0_val = ((wghtq * qlow + qs_val) / (wghtq + real(1.)) + qz0_val) * real(0.5);
            } else {
              thz0_val = (wghtt * thlow + thsk_val) / (wghtt + real(1.));
              qz0_val = (wghtq * qlow + qs_val) / (wghtq + real(1.));
            }
          }
          if (ustar_val >= USTR && ustar_val < USTC) {
            zu = z0_val;
            uz0_val = real(0.);
            vz0_val = real(0.);

            zt = FZT2 * std::sqrt(std::sqrt(z0_val * ustar_val * RVISC)) / ustar_val;
            zq = FZQ2 * zt;
            real wghtt = akhs_val * zt * RTVISC;
            real wghtq = akhs_val * zq * RQVISC;

            if (ntsd > 1) {
              thz0_val = ((wghtt * thlow + thsk_val) / (wghtt + real(1.)) + thz0_val) * real(0.5);
              qz0_val = ((wghtq * qlow + qs_val) / (wghtq + real(1.)) + qz0_val) * real(0.5);
            } else {
              thz0_val = (wghtt * thlow + thsk_val) / (wghtt + real(1.));
              qz0_val = (wghtq * qlow + qs_val) / (wghtq + real(1.));
            }
          }
        } else {
          zu = z0_val;
          uz0_val = real(0.);
          vz0_val = real(0.);
          zt = z0_val;
          thz0_val = thsk_val;
          zq = z0_val;
          qz0_val = qs_val;
        }

        real tem = (tlow + tz0) * real(0.5);
        real thm = (thelow + thz0_val) * real(0.5);
        real a_term = thm * P608;
        real b_term = (ELOCP / tem - real(1.) - P608) * thm;

        real dthv = ((thelow - thz0_val) * ((qlow + qz0_val + cwmlow) * (real(0.5) * P608) + real(1.)) +
                     (qlow - qz0_val + cwmlow) * a_term + cwmlow * b_term);

        real du2 = c_max((ulow - uz0_val) * (ulow - uz0_val) + (vlow - vz0_val) * (vlow - vz0_val), EPSU2);
        rib_val = btgx * dthv * zsl / du2;

        real zslu = zsl + zu;
        real zslt = zsl + zt;
        real rzsu = zslu / zu;
        real rzst = zslt / zt;
        real rlogu = rlog(rzsu);
        real rlogt = rlog(rzst);

        rlmo_val = elfc * akhs_val * dthv / (ustar_val * ustar_val * ustar_val);

        real zetalu = zslu * rlmo_val;
        real zetalt = zslt * rlmo_val;
        real zetau = zu * rlmo_val;
        real zetat = zt * rlmo_val;

        zetalu = c_min(c_max(zetalu, ZTMIN1), ZTMAX1);
        zetalt = c_min(c_max(zetalt, ZTMIN1), ZTMAX1);
        zetau = c_min(c_max(zetau, ZTMIN1 / rzsu), ZTMAX1 / rzsu);
        zetat = c_min(c_max(zetat, ZTMIN1 / rzst), ZTMAX1 / rzst);

        // Water functions: table 1
        real rz = (zetau - ZTMIN1) / DZETA1;
        int k_idx = static_cast<int>(rz);
        real rdzt = rz - real(k_idx);
        k_idx = c_min(k_idx, KZTM2);
        k_idx = c_max(k_idx, 0);
        psmz = (tables.psim1(k_idx + 1) - tables.psim1(k_idx)) * rdzt + tables.psim1(k_idx);

        rz = (zetalu - ZTMIN1) / DZETA1;
        k_idx = static_cast<int>(rz);
        rdzt = rz - real(k_idx);
        k_idx = c_min(k_idx, KZTM2);
        k_idx = c_max(k_idx, 0);
        real psmzl = (tables.psim1(k_idx + 1) - tables.psim1(k_idx)) * rdzt + tables.psim1(k_idx);

        real simm = psmzl - psmz + rlogu;

        rz = (zetat - ZTMIN1) / DZETA1;
        k_idx = static_cast<int>(rz);
        rdzt = rz - real(k_idx);
        k_idx = c_min(k_idx, KZTM2);
        k_idx = c_max(k_idx, 0);
        pshz = (tables.psih1(k_idx + 1) - tables.psih1(k_idx)) * rdzt + tables.psih1(k_idx);

        rz = (zetalt - ZTMIN1) / DZETA1;
        k_idx = static_cast<int>(rz);
        rdzt = rz - real(k_idx);
        k_idx = c_min(k_idx, KZTM2);
        k_idx = c_max(k_idx, 0);
        real pshzl = (tables.psih1(k_idx + 1) - tables.psih1(k_idx)) * rdzt + tables.psih1(k_idx);

        real simh = (pshzl - pshz + rlogt) * FH01;

        ustark = ustar_val * VKARMAN;
        akms_val = c_max(ustark / simm, cxchs);
        akhs_val = c_max(ustark / simh, cxchs);

        if (dthv <= real(0.)) {
          wstar2 = WWST2 * rpow(c_abs(btgh * akhs_val * dthv), real(2.) / real(3.));
        } else {
          wstar2 = real(0.);
        }
        ustar_val = c_max(std::sqrt(akms_val * std::sqrt(du2 + wstar2)), EPSUST);
      }
    } else {
      // LAND POINTS (module_sf_myjsfc.F:644-807)
      if (ntsd == 1) {
        qs_val = qlow;
      }
      zu = z0_val;
      uz0_val = real(0.);
      vz0_val = real(0.);
      zt = zu * ZTFC;
      thz0_val = thsk_val;
      zq = zt;
      qz0_val = qs_val;

      real tem = (tlow + tz0) * real(0.5);
      real thm = (thelow + thz0_val) * real(0.5);
      real a_term = thm * P608;
      real b_term = (ELOCP / tem - real(1.) - P608) * thm;

      real dthv = ((thelow - thz0_val) * ((qlow + qz0_val + cwmlow) * (real(0.5) * P608) + real(1.)) +
                   (qlow - qz0_val + cwmlow) * a_term + cwmlow * b_term);

      real du2 = c_max((ulow - uz0_val) * (ulow - uz0_val) + (vlow - vz0_val) * (vlow - vz0_val), EPSU2);
      rib_val = btgx * dthv * zsl / du2;

      real zslu = zsl + zu;
      real rzsu = zslu / zu;
      real rlogu = rlog(rzsu);
      real zslt = zsl + zu;

      real czil = real(0.1);
      real zilfc = -czil * VKARMAN * SQVISC;
      real czetmax = real(10.);
      real zzil;
      if (dthv > real(0.)) {
        if (rib_val < RIC) {
          zzil = zilfc * (real(1.0) + (rib_val / RIC) * (rib_val / RIC) * czetmax);
        } else {
          zzil = zilfc * (real(1.0) + czetmax);
        }
      } else {
        zzil = zilfc;
      }

      for (int itr = 1; itr <= ITRMX; ++itr) {
        zt = c_max(std::exp(zzil * std::sqrt(ustar_val * z0base_val)) * z0base_val, EPSZT);
        real rzst = zslt / zt;
        real rlogt = rlog(rzst);

        rlmo_val = elfc * akhs_val * dthv / (ustar_val * ustar_val * ustar_val);

        real zetalu = zslu * rlmo_val;
        real zetalt = zslt * rlmo_val;
        real zetau = zu * rlmo_val;
        real zetat = zt * rlmo_val;

        zetalu = c_min(c_max(zetalu, ZTMIN2), ZTMAX2);
        zetalt = c_min(c_max(zetalt, ZTMIN2), ZTMAX2);
        zetau = c_min(c_max(zetau, ZTMIN2 / rzsu), ZTMAX2 / rzsu);
        zetat = c_min(c_max(zetat, ZTMIN2 / rzst), ZTMAX2 / rzst);

        // Land functions: table 2
        real rz = (zetau - ZTMIN2) / DZETA2;
        int k_idx = static_cast<int>(rz);
        real rdzt = rz - real(k_idx);
        k_idx = c_min(k_idx, KZTM2);
        k_idx = c_max(k_idx, 0);
        psmz = (tables.psim2(k_idx + 1) - tables.psim2(k_idx)) * rdzt + tables.psim2(k_idx);

        rz = (zetalu - ZTMIN2) / DZETA2;
        k_idx = static_cast<int>(rz);
        rdzt = rz - real(k_idx);
        k_idx = c_min(k_idx, KZTM2);
        k_idx = c_max(k_idx, 0);
        real psmzl = (tables.psim2(k_idx + 1) - tables.psim2(k_idx)) * rdzt + tables.psim2(k_idx);

        real simm = psmzl - psmz + rlogu;

        rz = (zetat - ZTMIN2) / DZETA2;
        k_idx = static_cast<int>(rz);
        rdzt = rz - real(k_idx);
        k_idx = c_min(k_idx, KZTM2);
        k_idx = c_max(k_idx, 0);
        pshz = (tables.psih2(k_idx + 1) - tables.psih2(k_idx)) * rdzt + tables.psih2(k_idx);

        rz = (zetalt - ZTMIN2) / DZETA2;
        k_idx = static_cast<int>(rz);
        rdzt = rz - real(k_idx);
        k_idx = c_min(k_idx, KZTM2);
        k_idx = c_max(k_idx, 0);
        real pshzl = (tables.psih2(k_idx + 1) - tables.psih2(k_idx)) * rdzt + tables.psih2(k_idx);

        real simh = (pshzl - pshz + rlogt) * FH02;

        ustark = ustar_val * VKARMAN;
        akms_val = c_max(ustark / simm, cxchl);
        akhs_val = c_max(ustark / simh, cxchl);

        if (dthv <= real(0.)) {
          wstar2 = WWST2 * rpow(c_abs(btgh * akhs_val * dthv), real(2.) / real(3.));
        } else {
          wstar2 = real(0.);
        }
        ustar_val = c_max(std::sqrt(akms_val * std::sqrt(du2 + wstar2)), EPSUST);
      }
    }

    // Diagnostics at 2m and 10m (module_sf_myjsfc.F:836-1020)
    real umflx = akms_val * (ulow - uz0_val);
    real vmflx = akms_val * (vlow - vz0_val);
    real hsflx = akhs_val * (thlow - thz0_val);
    real hlflx = akhs_val * (qlow - qz0_val);

    real zu10 = zu + real(10.);
    real zt02 = zt + real(2.);
    real zt10 = zt + real(10.);

    real rlnu10 = rlog(zu10 / zu);
    real rlnt02 = rlog(zt02 / zt);
    real rlnt10 = rlog(zt10 / zt);

    real ztau10 = zu10 * rlmo_val;
    real ztat02 = zt02 * rlmo_val;
    real ztat10 = zt10 * rlmo_val;

    real akms10, akhs02, akhs10;
    if (seamask > real(0.5)) {
      ztau10 = c_min(c_max(ztau10, ZTMIN1), ZTMAX1);
      ztat02 = c_min(c_max(ztat02, ZTMIN1), ZTMAX1);
      ztat10 = c_min(c_max(ztat10, ZTMIN1), ZTMAX1);

      real rz = (ztau10 - ZTMIN1) / DZETA1;
      int k_idx = static_cast<int>(rz);
      real rdzt = rz - real(k_idx);
      k_idx = c_min(k_idx, KZTM2);
      k_idx = c_max(k_idx, 0);
      real psm10 = (tables.psim1(k_idx + 1) - tables.psim1(k_idx)) * rdzt + tables.psim1(k_idx);
      real simm10 = psm10 - psmz + rlnu10;

      rz = (ztat02 - ZTMIN1) / DZETA1;
      k_idx = static_cast<int>(rz);
      rdzt = rz - real(k_idx);
      k_idx = c_min(k_idx, KZTM2);
      k_idx = c_max(k_idx, 0);
      real psh02 = (tables.psih1(k_idx + 1) - tables.psih1(k_idx)) * rdzt + tables.psih1(k_idx);
      real simh02 = (psh02 - pshz + rlnt02) * FH01;

      rz = (ztat10 - ZTMIN1) / DZETA1;
      k_idx = static_cast<int>(rz);
      rdzt = rz - real(k_idx);
      k_idx = c_min(k_idx, KZTM2);
      k_idx = c_max(k_idx, 0);
      real psh10 = (tables.psih1(k_idx + 1) - tables.psih1(k_idx)) * rdzt + tables.psih1(k_idx);
      real simh10 = (psh10 - pshz + rlnt10) * FH01;

      akms10 = c_max(ustark / simm10, cxchs);
      akhs02 = c_max(ustark / simh02, cxchs);
      akhs10 = c_max(ustark / simh10, cxchs);
    } else {
      ztau10 = c_min(c_max(ztau10, ZTMIN2), ZTMAX2);
      ztat02 = c_min(c_max(ztat02, ZTMIN2), ZTMAX2);
      ztat10 = c_min(c_max(ztat10, ZTMIN2), ZTMAX2);

      real rz = (ztau10 - ZTMIN2) / DZETA2;
      int k_idx = static_cast<int>(rz);
      real rdzt = rz - real(k_idx);
      k_idx = c_min(k_idx, KZTM2);
      k_idx = c_max(k_idx, 0);
      real psm10 = (tables.psim2(k_idx + 1) - tables.psim2(k_idx)) * rdzt + tables.psim2(k_idx);
      real simm10 = psm10 - psmz + rlnu10;

      rz = (ztat02 - ZTMIN2) / DZETA2;
      k_idx = static_cast<int>(rz);
      rdzt = rz - real(k_idx);
      k_idx = c_min(k_idx, KZTM2);
      k_idx = c_max(k_idx, 0);
      real psh02 = (tables.psih2(k_idx + 1) - tables.psih2(k_idx)) * rdzt + tables.psih2(k_idx);
      real simh02 = (psh02 - pshz + rlnt02) * FH02;

      rz = (ztat10 - ZTMIN2) / DZETA2;
      k_idx = static_cast<int>(rz);
      rdzt = rz - real(k_idx);
      k_idx = c_min(k_idx, KZTM2);
      k_idx = c_max(k_idx, 0);
      real psh10 = (tables.psih2(k_idx + 1) - tables.psih2(k_idx)) * rdzt + tables.psih2(k_idx);
      real simh10 = (psh10 - pshz + rlnt10) * FH02;

      akms10 = c_max(ustark / simm10, cxchl);
      akhs02 = c_max(ustark / simh02, cxchl);
      akhs10 = c_max(ustark / simh10, cxchl);
    }

    real u10_val = umflx / akms10 + uz0_val;
    real v10_val = vmflx / akms10 + vz0_val;
    real tshltr_val = hsflx / akhs02 + thz0_val;

    if ((thlow > thz0_val && (tshltr_val < thz0_val || tshltr_val > thlow)) ||
        (thlow < thz0_val && (tshltr_val > thz0_val || tshltr_val < thlow))) {
      tshltr_val = thz0_val + real(2.) * rdz * (thlow - thz0_val);
    }

    real th10_val = hsflx / akhs10 + thz0_val;
    if ((thlow > thz0_val && (th10_val < thz0_val || th10_val > thlow)) ||
        (thlow < thz0_val && (th10_val > thz0_val || th10_val < thlow))) {
      th10_val = thz0_val + real(10.) * rdz * (thlow - thz0_val);
    }

    real qshltr_val = hlflx / akhs02 + qz0_val;
    real q10_val = hlflx / akhs10 + qz0_val;
    real term1 = real(-0.068283) / tlow;
    real pshltr_val = psfc_val * std::exp(term1);

    real u10e_val = u10_val;
    real v10e_val = v10_val;

    if (seamask < real(0.5)) {
      real zuuz = c_min(zu * real(0.50), real(0.18));
      real zu_e = c_max(zu * real(0.35), zuuz);
      real zu10_e = zu_e + real(10.);
      real rzsu_e = zu10_e / zu_e;
      real rlnu10_e = rlog(rzsu_e);

      real zetau_e = zu_e * rlmo_val;
      real ztau10_e = zu10_e * rlmo_val;
      ztau10_e = c_min(c_max(ztau10_e, ZTMIN2), ZTMAX2);
      zetau_e = c_min(c_max(zetau_e, ZTMIN2 / rzsu_e), ZTMAX2 / rzsu_e);

      real rz = (ztau10_e - ZTMIN2) / DZETA2;
      int k_idx = static_cast<int>(rz);
      real rdzt = rz - real(k_idx);
      k_idx = c_min(k_idx, KZTM2);
      k_idx = c_max(k_idx, 0);
      real psm10_e = (tables.psim2(k_idx + 1) - tables.psim2(k_idx)) * rdzt + tables.psim2(k_idx);
      real simm10_e = psm10_e - psmz + rlnu10_e;
      real ekms10 = c_max(ustark / simm10_e, cxchl);

      u10e_val = umflx / ekms10 + uz0_val;
      v10e_val = vmflx / ekms10 + vz0_val;
    }

    // Set other WRF driver arrays (module_sf_myjsfc.F:1029-1053)
    real rlow = plow / (R_D * tlow);
    real hfx_val = -rlow * CP * hsflx;
    real qfx_val = -rlow * hlflx * wetm;
    real flx_lh_val = XLV * qfx_val;
    real flhc_val = rlow * CP * akhs_val;
    real flqc_val = rlow * akhs_val * wetm;

    real qgh_val = ((real(1.) - seamask) * PQ0 + seamask * PQ0SEA) / plow *
                   std::exp(A2S * (tlow - A3S) / (tlow - A4S));
    qgh_val = qgh_val / (real(1.) - qgh_val);
    real cpm_val = CP * (real(1.) + real(0.8) * qlow);

    if (seamask > real(0.5)) {
      qs_val = PQ0SEA / psfc_val * std::exp(A2S * (tsk(i, j) - A3S) / (tsk(i, j) - A4S));
      qs_val = qs_val / (real(1.) - qs_val);
    }

    // Remove supersaturation at 2m and 10m (module_sf_myjsfc.F:324-350)
    real rapa = apesfc;
    real th02p = tshltr_val;
    real th10p = th10_val;
    real th02_val = tshltr_val;

    real rapa02 = rapa - GOCP02 / th02p;
    real rapa10 = rapa - GOCP10 / th10p;

    real t02p = th02p * rapa02;
    real t10p = th10p * rapa10;
    real t02_val = th02_val * apesfc;

    real p02p = rpow(rapa02, RCAP) * p1000mb;
    real p10p = rpow(rapa10, RCAP) * p1000mb;

    constexpr real A2 = real(17.2693882);
    constexpr real A3 = real(273.16);
    constexpr real A4 = real(35.86);

    real qs02 = PQ0 / p02p * std::exp(A2 * (t02p - A3) / (t02p - A4));
    real qs10 = PQ0 / p10p * std::exp(A2 * (t10p - A3) / (t10p - A4));

    if (qshltr_val > qs02) qshltr_val = qs02;
    if (q10_val > qs10) q10_val = qs10;
    real q02_val = qshltr_val / (real(1.) - qshltr_val);

    // Write back inouts and outputs
    znt(i, j) = z0_val;
    ustar(i, j) = ustar_val;
    uz0(i, j) = uz0_val;
    vz0(i, j) = vz0_val;
    thz0(i, j) = thz0_val;
    qz0(i, j) = qz0_val;
    qsfc(i, j) = qs_val;
    akms(i, j) = akms_val;
    akhs(i, j) = akhs_val;
    rmol(i, j) = rlmo_val;
    rib(i, j) = rib_val;

    chs(i, j) = akhs_val;
    chs2(i, j) = akhs02;
    cqs2(i, j) = akhs02;
    hfx(i, j) = hfx_val;
    qfx(i, j) = qfx_val;
    flx_lh(i, j) = flx_lh_val;
    flhc(i, j) = flhc_val;
    flqc(i, j) = flqc_val;
    qgh(i, j) = qgh_val;
    cpm(i, j) = cpm_val;

    u10(i, j) = u10_val;
    v10(i, j) = v10_val;
    t02(i, j) = t02_val;
    th02(i, j) = th02_val;
    tshltr(i, j) = tshltr_val;
    th10(i, j) = th10_val;
    q02(i, j) = q02_val;
    qshltr(i, j) = qshltr_val;
    q10(i, j) = q10_val;
    pshltr(i, j) = pshltr_val;
    u10e(i, j) = u10e_val;
    v10e(i, j) = v10e_val;
  }
};

// Dispatcher functions
void run_myjsfc(
    View2<const real> ht,
    View3<const real> dz,
    View3<const real> pmid,
    View3<const real> pint,
    View3<const real> th,
    View3<const real> t,
    View3<const real> qv,
    View3<const real> qc,
    View3<const real> u,
    View3<const real> v,
    View3<const real> q2,
    View2<const real> tsk,
    View2<const real> xland,
    View2<const int> ivgtyp,
    View2<const real> z0base,
    View2<const real> mavail,
    View2<const int> lowlyr,
    View2<real> qsfc,
    View2<real> thz0,
    View2<real> qz0,
    View2<real> uz0,
    View2<real> vz0,
    View2<real> ustar,
    View2<real> znt,
    View2<real> pblh,
    View2<real> rmol,
    View2<real> akhs,
    View2<real> akms,
    View2<real> rib,
    View2<real> chs,
    View2<real> chs2,
    View2<real> cqs2,
    View2<real> hfx,
    View2<real> qfx,
    View2<real> flx_lh,
    View2<real> flhc,
    View2<real> flqc,
    View2<real> qgh,
    View2<real> cpm,
    View2<real> ct,
    View2<real> u10,
    View2<real> v10,
    View2<real> t02,
    View2<real> th02,
    View2<real> tshltr,
    View2<real> th10,
    View2<real> q02,
    View2<real> qshltr,
    View2<real> q10,
    View2<real> pshltr,
    View2<real> u10e,
    View2<real> v10e,
    const MYJSFCTables& tables,
    const Range2D& range,
    int itimestep,
    int isurban,
    int iz0tlnd,
    real p1000mb,
    int nz);

}  // namespace cppwrf::phys::sf_myjsfc
