#pragma once
// Port of WRF 4.8.0 phys/module_bl_myjpbl.F (myjpbl, mixlen, prodq2, difcof, vdifq, vdifh, vdifv, myjpblinit)
// Mellor-Yamada-Janjic (MYJ) PBL scheme (bl_pbl_physics = 2).

#include <cmath>
#include "core/backend.hpp"
#include "core/field.hpp"
#include "core/precision.hpp"
#include "phys/driver/wrf_constants.hpp"

namespace cppwrf::phys::myjpbl {

template <typename T>
CPPWRF_HD inline T c_min(T a, T b) { return (a < b) ? a : b; }

template <typename T>
CPPWRF_HD inline T c_max(T a, T b) { return (a > b) ? a : b; }

template <typename T>
CPPWRF_HD inline T c_abs(T a) { return (a < T(0)) ? -a : a; }

static constexpr int kMaxLevels = 128;

// Physical constants from module_model_constants.F and module_bl_myjpbl.F:29-125
static constexpr real G = wrfc::g;
static constexpr real R_D = wrfc::r_d;
static constexpr real CP = wrfc::cp;
static constexpr real R_V = wrfc::r_v;
static constexpr real XLV = wrfc::xlv;
static constexpr real XLS = wrfc::xls;
static constexpr real EP_1 = wrfc::ep_1;

static constexpr real PI = real(3.1415926);
static constexpr real VKARMAN = real(0.4);
static constexpr real CAPA = R_D / CP;
static constexpr real RLIVWV = XLS / XLV;
static constexpr real ELOCP = real(2.72e6) / CP;
static constexpr real EPS1 = real(1.e-12);
static constexpr real EPS2 = real(0.);
static constexpr real EPSL = real(0.32);
static constexpr real EPSRU = real(1.e-7);
static constexpr real EPSRS = real(1.e-7);
static constexpr real EPSTRB = real(1.e-24);
static constexpr real EPSA = real(1.e-8);
static constexpr real EPSIT = real(1.e-4);
static constexpr real EPSU2 = real(1.e-4);
static constexpr real EPSUST = real(0.07);
static constexpr real FH = real(1.01);
static constexpr real ALPH = real(0.30);
static constexpr real BETA = real(1.) / real(273.);
static constexpr real EL0MAX = real(1000.);
static constexpr real EL0MIN = real(1.);
static constexpr real ELFC = real(0.23) * real(0.5);
static constexpr real GAM1 = real(0.2222222222222222222);
static constexpr real PRT = real(1.);

static constexpr real A1 = real(0.659888514560862645);
static constexpr real A2x = real(0.6574209922667784586);
static constexpr real B1 = real(11.87799326209552761);
static constexpr real B2 = real(7.226971804046074028);
static constexpr real C1 = real(0.000830955950095854396);

static constexpr real A2S = real(17.2693882);
static constexpr real A3S = real(273.16);
static constexpr real A4S = real(35.86);

static constexpr real ELZ0 = real(0.);
static constexpr real ESQ = real(5.0);
static constexpr real EXCM = real(0.001);
static constexpr real FHNEU = real(0.8);
static constexpr real GLKBR = real(10.);
static constexpr real GLKBS = real(30.);
static constexpr real QVISC = real(2.1e-5);
static constexpr real RFC = real(0.191);
static constexpr real RIC = real(0.505);
static constexpr real SMALL = real(0.35);
static constexpr real SQPR = real(0.84);
static constexpr real SQSC = real(0.84);
static constexpr real SQVISC = real(258.2);
static constexpr real TVISC = real(2.1e-5);
static constexpr real USTC = real(0.7);
static constexpr real USTR = real(0.225);
static constexpr real VISC = real(1.5e-5);
static constexpr real WOLD = real(0.15);
static constexpr real WWST = real(1.2);
static constexpr real ZTMAX = real(1.);
static constexpr real ZTFC = real(1.);
static constexpr real ZTMIN = real(-5.);

static constexpr real SEAFC = real(0.98);
static constexpr real PQ0 = real(379.90516);
static constexpr real PQ0SEA = PQ0 * SEAFC;

static constexpr real BTG = BETA * G;
static constexpr real CZIV = SMALL * GLKBS;
static constexpr real ESQHF = real(0.5) * real(5.0);
static constexpr real GRRS = GLKBR / GLKBS;
static constexpr real RB1 = real(1.) / B1;
static constexpr real RTVISC = real(1.) / TVISC;
static constexpr real RVISC = real(1.) / VISC;
static constexpr real ZQRZT = SQSC / SQPR;

static constexpr real ADNH = real(9.) * A1 * A2x * A2x * (real(12.) * A1 + real(3.) * B2) * BTG * BTG;
static constexpr real ADNM = real(18.) * A1 * A1 * A2x * (B2 - real(3.) * A2x) * BTG;
static constexpr real ANMH = -real(9.) * A1 * A2x * A2x * BTG * BTG;
static constexpr real ANMM = -real(3.) * A1 * A2x * (real(3.) * A2x + real(3.) * B2 * C1 + real(18.) * A1 * C1 - B2) * BTG;
static constexpr real BDNH = real(3.) * A2x * (real(7.) * A1 + B2) * BTG;
static constexpr real BDNM = real(6.) * A1 * A1;
static constexpr real BEQH = A2x * B1 * BTG + real(3.) * A2x * (real(7.) * A1 + B2) * BTG;
static constexpr real BEQM = -A1 * B1 * (real(1.) - real(3.) * C1) + real(6.) * A1 * A1;
static constexpr real BNMH = -A2x * BTG;
static constexpr real BNMM = A1 * (real(1.) - real(3.) * C1);
static constexpr real BSHH = real(9.) * A1 * A2x * A2x * BTG;
static constexpr real BSHM = real(18.) * A1 * A1 * A2x * C1;
static constexpr real BSMH = -real(3.) * A1 * A2x * (real(3.) * A2x + real(3.) * B2 * C1 + real(12.) * A1 * C1 - B2) * BTG;
static constexpr real CESH = A2x;
static constexpr real CESM = A1 * (real(1.) - real(3.) * C1);
static constexpr real CNV = EP_1 * G / BTG;
static constexpr real ELFCS = VKARMAN * BTG;
static constexpr real FZQ1 = RTVISC * QVISC * ZQRZT;
static constexpr real FZQ2 = RTVISC * QVISC * ZQRZT;
static constexpr real FZT1 = RVISC * TVISC * SQPR;
static constexpr real FZT2 = CZIV * GRRS * TVISC * SQPR;
static constexpr real FZU1 = CZIV * VISC;
static constexpr real PIHF = real(0.5) * PI;
static constexpr real RFAC = RIC / (FHNEU * RFC * RFC);
static constexpr real RQVISC = real(1.) / QVISC;
static constexpr real RRIC = real(1.) / RIC;
static constexpr real USTFC = real(0.018) / G;
static constexpr real WNEW = real(1.) - WOLD;
static constexpr real WWST2 = WWST * WWST;

static constexpr real AEQH = real(9.) * A1 * A2x * A2x * B1 * BTG * BTG +
                             real(9.) * A1 * A2x * A2x * (real(12.) * A1 + real(3.) * B2) * BTG * BTG;
static constexpr real AEQM = real(3.) * A1 * A2x * B1 * (real(3.) * A2x + real(3.) * B2 * C1 + real(18.) * A1 * C1 - B2) * BTG +
                             real(18.) * A1 * A1 * A2x * (B2 - real(3.) * A2x) * BTG;

static constexpr real REQU = -AEQH / AEQM;
static constexpr real EPSGH = real(1.e-9);
static constexpr real EPSGM = REQU * EPSGH;

static constexpr real UBRYL = (real(18.) * REQU * A1 * A1 * A2x * B2 * C1 * BTG +
                               real(9.) * A1 * A2x * A2x * B2 * BTG * BTG) /
                              (REQU * ADNM + ADNH);
static constexpr real UBRY = (real(1.) + EPSRS) * UBRYL;
static constexpr real UBRY3 = real(3.) * UBRY;

static constexpr real AUBH = real(27.) * A1 * A2x * A2x * B2 * BTG * BTG - ADNH * UBRY3;
static constexpr real AUBM = real(54.) * A1 * A1 * A2x * B2 * C1 * BTG - ADNM * UBRY3;
static constexpr real BUBH = (real(9.) * A1 * A2x + real(3.) * A2x * B2) * BTG - BDNH * UBRY3;
static constexpr real BUBM = real(18.) * A1 * A1 * C1 - BDNM * UBRY3;
static constexpr real CUBR = real(1.) - UBRY3;
static constexpr real RCUBR = real(1.) / CUBR;

static constexpr real P608 = R_V / R_D - real(1.);
static constexpr real EPSQ2 = real(0.2);

static constexpr real A2 = real(17.2693882);
static constexpr real A3 = real(273.16);
static constexpr real A4 = real(35.86);

// ----------------------------------------------------------------------
// Subroutine MIXLEN (Level 2.5 Mixing Length, module_bl_myjpbl.F:771-963)
// ----------------------------------------------------------------------
CPPWRF_HD inline void mixlen(
    int lmh,
    const real* u,
    const real* v,
    const real* t,
    const real* the,
    const real* q,
    const real* cwm,
    const real* q2,
    const real* z,
    real* gm,
    real* gh,
    real* el,
    real& pblh,
    int& lpbl,
    int& lmxl,
    real& ct,
    real& mixht,
    real* q1,
    real* dth,
    real* elm,
    real* rel) {
  lpbl = 1;
  for (int k = lmh - 1; k >= 1; --k) {
    if (q2[k] <= EPSQ2 * FH) {
      lpbl = k;
      break;
    }
  }

  pblh = z[lpbl + 1] - z[lmh + 1];

  for (int k = 1; k <= lmh; ++k) {
    q1[k] = real(0.);
  }

  for (int k = 1; k <= lmh - 1; ++k) {
    dth[k] = the[k] - the[k + 1];
  }

  for (int k = lmh - 2; k >= 1; --k) {
    if (dth[k] > real(0.) && dth[k + 1] <= real(0.)) {
      dth[k] = dth[k] + ct;
      break;
    }
  }

  ct = real(0.);

  for (int k = 1; k <= lmh - 1; ++k) {
    real rdz = real(2.) / (z[k] - z[k + 2]);
    real du = u[k] - u[k + 1];
    real dv = v[k] - v[k + 1];
    real gml = (du * du + dv * dv) * rdz * rdz;
    gm[k] = c_max(gml, EPSGM);

    real tem = (t[k] + t[k + 1]) * real(0.5);
    real thm = (the[k] + the[k + 1]) * real(0.5);

    real a_term = thm * P608;
    real b_term = (ELOCP / tem - real(1.) - P608) * thm;

    real ghl = (dth[k] * ((q[k] + q[k + 1] + cwm[k] + cwm[k + 1]) * (real(0.5) * P608) + real(1.)) +
                (q[k] - q[k + 1] + cwm[k] - cwm[k + 1]) * a_term +
                (cwm[k] - cwm[k + 1]) * b_term) * rdz;

    if (c_abs(ghl) <= EPSGH) ghl = EPSGH;
    gh[k] = ghl;
  }

  lmxl = lmh;
  for (int k = 1; k <= lmh - 1; ++k) {
    real gml = gm[k];
    real ghl = gh[k];

    if (ghl >= EPSGH) {
      if (gml / ghl <= REQU) {
        elm[k] = EPSL;
        lmxl = k;
      } else {
        real aubr = (AUBM * gml + AUBH * ghl) * ghl;
        real bubr = BUBM * gml + BUBH * ghl;
        real qol2st = (-real(0.5) * bubr + std::sqrt(bubr * bubr * real(0.25) - aubr * CUBR)) * RCUBR;
        real eloq2x = real(1.) / qol2st;
        elm[k] = c_max(std::sqrt(eloq2x * q2[k]), EPSL);
      }
    } else {
      real aden = (ADNM * gml + ADNH * ghl) * ghl;
      real bden = BDNM * gml + BDNH * ghl;
      real qol2un = -real(0.5) * bden + std::sqrt(bden * bden * real(0.25) - aden);
      real eloq2x = real(1.) / (qol2un + EPSRU);
      elm[k] = c_max(std::sqrt(eloq2x * q2[k]), EPSL);
    }
  }

  if (elm[lmh - 1] == EPSL) lmxl = lmh;

  real blmx = z[lmxl] - z[lmh + 1];
  mixht = blmx;

  for (int k = lpbl; k <= lmh; ++k) {
    q1[k] = std::sqrt(q2[k]);
  }

  real szq = real(0.);
  real sq = real(0.);
  for (int k = 1; k <= lmh - 1; ++k) {
    real qdzl = (q1[k] + q1[k + 1]) * (z[k + 1] - z[k + 2]);
    szq = (z[k + 1] + z[k + 2] - z[lmh + 1] - z[lmh + 1]) * qdzl + szq;
    sq = qdzl + sq;
  }

  real el0 = c_min(ALPH * szq * real(0.5) / sq, EL0MAX);
  el0 = c_max(el0, EL0MIN);

  int lpblm = c_max(lpbl - 1, 1);
  for (int k = 1; k <= lpblm; ++k) {
    el[k] = c_min((z[k] - z[k + 2]) * ELFC, elm[k]);
    rel[k] = el[k] / elm[k];
  }

  if (lpbl < lmh) {
    for (int k = lpbl; k <= lmh - 1; ++k) {
      real vkrmz = (z[k + 1] - z[lmh + 1]) * VKARMAN;
      el[k] = c_min(vkrmz / (vkrmz / el0 + real(1.)), elm[k]);
      rel[k] = el[k] / elm[k];
    }
  }

  for (int k = lpbl + 1; k <= lmh - 2; ++k) {
    real srel = c_min(((rel[k - 1] + rel[k + 1]) * real(0.5) + rel[k]) * real(0.5), rel[k]);
    el[k] = c_max(srel * elm[k], EPSL);
  }
}

// ----------------------------------------------------------------------
// Subroutine PRODQ2 (Level 2.5 Q2 Production/Dissipation, module_bl_myjpbl.F:967-1167)
// ----------------------------------------------------------------------
CPPWRF_HD inline void prodq2(
    int lmh,
    real dtturbl,
    real ustar,
    const real* gm,
    const real* gh,
    real* el,
    real* q2) {
  for (int k = 1; k <= lmh - 1; ++k) {
    real gml = gm[k];
    real ghl = gh[k];

    real aequ = (AEQM * gml + AEQH * ghl) * ghl;
    real bequ = BEQM * gml + BEQH * ghl;

    real eqol2 = -real(0.5) * bequ + std::sqrt(bequ * bequ * real(0.25) - aequ);

    if ((gml + ghl * ghl <= EPSTRB) ||
        (ghl >= EPSGH && gml / ghl <= REQU) ||
        (eqol2 <= EPS2)) {
      q2[k] = EPSQ2;
      el[k] = EPSL;
    } else {
      real anum = (ANMM * gml + ANMH * ghl) * ghl;
      real bnum = BNMM * gml + BNMH * ghl;
      real aden = (ADNM * gml + ADNH * ghl) * ghl;
      real bden = BDNM * gml + BDNH * ghl;
      real cden = real(1.);

      real arhs = -(anum * bden - bnum * aden) * real(2.);
      real brhs = -anum * real(4.);
      real crhs = -bnum * real(2.);

      real dloq1 = el[k] / std::sqrt(q2[k]);

      real eloq21 = real(1.) / eqol2;
      real eloq11 = std::sqrt(eloq21);
      real eloq31 = eloq21 * eloq11;
      real eloq41 = eloq21 * eloq21;
      real eloq51 = eloq21 * eloq31;

      real rden1 = real(1.) / (aden * eloq41 + bden * eloq21 + cden);
      real rhsp1 = (arhs * eloq51 + brhs * eloq31 + crhs * eloq11) * rden1 * rden1;

      real eloq12 = eloq11 + (dloq1 - eloq11) * std::exp(rhsp1 * dtturbl);
      eloq12 = c_max(eloq12, EPS1);

      real eloq22 = eloq12 * eloq12;
      real eloq32 = eloq22 * eloq12;
      real eloq42 = eloq22 * eloq22;
      real eloq52 = eloq22 * eloq32;

      real rden2 = real(1.) / (aden * eloq42 + bden * eloq22 + cden);
      real rhs2 = -(anum * eloq42 + bnum * eloq22) * rden2 + RB1;
      real rhsp2 = (arhs * eloq52 + brhs * eloq32 + crhs * eloq12) * rden2 * rden2;
      real rhst2 = rhs2 / rhsp2;

      real eloq13 = eloq12 - rhst2 + (rhst2 + dloq1 - eloq12) * std::exp(rhsp2 * dtturbl);
      eloq13 = c_max(eloq13, EPS1);

      real eloqn = eloq13;
      if (eloqn > EPS1) {
        q2[k] = el[k] * el[k] / (eloqn * eloqn);
        q2[k] = c_max(q2[k], EPSQ2);
        if (q2[k] == EPSQ2) {
          el[k] = EPSL;
        }
      } else {
        q2[k] = EPSQ2;
        el[k] = EPSL;
      }
    }
  }

  q2[lmh] = c_max(rpow(B1, real(2.) / real(3.)) * ustar * ustar, EPSQ2);
}

// ----------------------------------------------------------------------
// Subroutine DIFCOF (Level 2.5 Diffusion Coefficients, module_bl_myjpbl.F:1172-1325)
// ----------------------------------------------------------------------
CPPWRF_HD inline void difcof(
    int lmh,
    const real* gm,
    const real* gh,
    const real* el,
    const real* q2,
    const real* z,
    real* akm,
    real* akh) {
  for (int k = 1; k <= lmh - 1; ++k) {
    real ell = el[k];
    real eloq2 = ell * ell / q2[k];
    real eloq4 = eloq2 * eloq2;

    real gml = gm[k];
    real ghl = gh[k];

    real aden = (ADNM * gml + ADNH * ghl) * ghl;
    real bden = BDNM * gml + BDNH * ghl;
    real cden = real(1.);

    real besm = BSMH * ghl;
    real besh = BSHM * gml + BSHH * ghl;

    real rden = real(1.) / (aden * eloq4 + bden * eloq2 + cden);

    real esm = (besm * eloq2 + CESM) * rden;
    real esh = (besh * eloq2 + CESH) * rden;

    real rdz = real(2.) / (z[k] - z[k + 2]);
    real q1l = std::sqrt(q2[k]);
    real elqdz = ell * q1l * rdz;

    akm[k] = elqdz * esm;
    akh[k] = elqdz * esh;
  }
}

// ----------------------------------------------------------------------
// Subroutine VDIFQ (Vertical Diffusion of Q2, module_bl_myjpbl.F:1330-1406)
// ----------------------------------------------------------------------
CPPWRF_HD inline void vdifq(
    int lmh,
    real dtdif,
    real* q2,
    const real* el,
    const real* z,
    real* akq,
    real* rsq2,
    real* cm,
    real* cr,
    real* dtoz) {
  const real esqhf = real(0.5) * ESQ;

  for (int k = 1; k <= lmh - 2; ++k) {
    dtoz[k] = (dtdif + dtdif) / (z[k] - z[k + 2]);
    akq[k] = std::sqrt((q2[k] + q2[k + 1]) * real(0.5)) * (el[k] + el[k + 1]) * esqhf / (z[k + 1] - z[k + 2]);
    cr[k] = -dtoz[k] * akq[k];
  }

  cm[1] = dtoz[1] * akq[1] + real(1.);
  rsq2[1] = q2[1];

  for (int k = 2; k <= lmh - 2; ++k) {
    real cf = -dtoz[k] * akq[k - 1] / cm[k - 1];
    cm[k] = -cr[k - 1] * cf + (akq[k - 1] + akq[k]) * dtoz[k] + real(1.);
    rsq2[k] = -rsq2[k - 1] * cf + q2[k];
  }

  real dtozs = (dtdif + dtdif) / (z[lmh - 1] - z[lmh + 1]);
  real akqs = std::sqrt((q2[lmh - 1] + q2[lmh]) * real(0.5)) * (el[lmh - 1] + ELZ0) * esqhf / (z[lmh] - z[lmh + 1]);

  real cf = -dtozs * akq[lmh - 2] / cm[lmh - 2];

  q2[lmh - 1] = (dtozs * akqs * q2[lmh] - rsq2[lmh - 2] * cf + q2[lmh - 1]) /
                ((akq[lmh - 2] + akqs) * dtozs - cr[lmh - 2] * cf + real(1.));

  for (int k = lmh - 2; k >= 1; --k) {
    q2[k] = (-cr[k] * q2[k + 1] + rsq2[k]) / cm[k];
  }
}

// ----------------------------------------------------------------------
// Subroutine VDIFH (Vertical Diffusion of Mass Species, module_bl_myjpbl.F:1411-1519)
// Solved per species in-place using caller-allocated 1D scratch arrays to avoid
// multi-KB local arrays inside inline device functions.
// ----------------------------------------------------------------------
CPPWRF_HD inline void vdifh_1spec(
    real dtdif,
    int lmh,
    int lpbl,
    real sz0_m,
    real rkhs_dens,
    real clow_m,
    real cts_m,
    real* sp,
    const real* rkh,
    const real* zhk,
    const real* rho,
    real* rkct,
    real* rss,
    real* cm,
    real* cr,
    real* dtoz) {
  for (int k = 1; k <= lmh - 1; ++k) {
    dtoz[k] = dtdif / (zhk[k] - zhk[k + 1]);
    cr[k] = -dtoz[k] * rkh[k];

    if (k < lpbl) {
      rkct[k] = real(0.);
    } else {
      real rkhz = rkh[k] * (zhk[k] - zhk[k + 2]);
      rkct[k] = rkhz * cts_m * real(0.5);
    }
  }

  real rhok_val = rho[1];
  cm[1] = dtoz[1] * rkh[1] + rhok_val;
  rss[1] = -rkct[1] * dtoz[1] + sp[1] * rhok_val;

  for (int k = 2; k <= lmh - 1; ++k) {
    real dtozl = dtoz[k];
    real cf = -dtozl * rkh[k - 1] / cm[k - 1];
    rhok_val = rho[k];
    cm[k] = -cr[k - 1] * cf + (rkh[k - 1] + rkh[k]) * dtozl + rhok_val;
    rss[k] = -rss[k - 1] * cf + (rkct[k - 1] - rkct[k]) * dtozl + sp[k] * rhok_val;
  }

  real dtozs = dtdif / (zhk[lmh] - zhk[lmh + 1]);
  real rkhh = rkh[lmh - 1];
  real cf = -dtozs * rkhh / cm[lmh - 1];
  real cmb = cr[lmh - 1] * cf;
  rhok_val = rho[lmh];

  real rkss = rkhs_dens * clow_m;
  real cmsb = -cmb + (rkhh + rkss) * dtozs + rhok_val;
  real rssb = -rss[lmh - 1] * cf + rkct[lmh - 1] * dtozs + sp[lmh] * rhok_val;
  sp[lmh] = (dtozs * rkss * sz0_m + rssb) / cmsb;

  for (int k = lmh - 1; k >= 1; --k) {
    real rcml = real(1.) / cm[k];
    sp[k] = (-cr[k] * sp[k + 1] + rss[k]) * rcml;
  }
}

// ----------------------------------------------------------------------
// Subroutine VDIFV (Vertical Diffusion of Velocity Components, module_bl_myjpbl.F:1613-1689)
// ----------------------------------------------------------------------
CPPWRF_HD inline void vdifv(
    int lmh,
    real dtdif,
    real uz0,
    real vz0,
    real rkms_dens,
    real* u,
    real* v,
    const real* rkm,
    const real* z,
    const real* rho,
    real* rsu,
    real* rsv,
    real* cm,
    real* cr,
    real* dtoz) {
  for (int k = 1; k <= lmh - 1; ++k) {
    dtoz[k] = dtdif / (z[k] - z[k + 1]);
    cr[k] = -dtoz[k] * rkm[k];
  }

  real rhok_val = rho[1];
  cm[1] = dtoz[1] * rkm[1] + rhok_val;
  rsu[1] = u[1] * rhok_val;
  rsv[1] = v[1] * rhok_val;

  for (int k = 2; k <= lmh - 1; ++k) {
    real dtozl = dtoz[k];
    real cf = -dtozl * rkm[k - 1] / cm[k - 1];
    rhok_val = rho[k];
    cm[k] = -cr[k - 1] * cf + (rkm[k - 1] + rkm[k]) * dtozl + rhok_val;
    rsu[k] = -rsu[k - 1] * cf + u[k] * rhok_val;
    rsv[k] = -rsv[k - 1] * cf + v[k] * rhok_val;
  }

  real dtozs = dtdif / (z[lmh] - z[lmh + 1]);
  real rkmh = rkm[lmh - 1];
  real cf = -dtozs * rkmh / cm[lmh - 1];
  rhok_val = rho[lmh];
  real rcmvb = real(1.) / ((rkmh + rkms_dens) * dtozs - cr[lmh - 1] * cf + rhok_val);
  real dtozak = dtozs * rkms_dens;

  u[lmh] = (dtozak * uz0 - rsu[lmh - 1] * cf + u[lmh] * rhok_val) * rcmvb;
  v[lmh] = (dtozak * vz0 - rsv[lmh - 1] * cf + v[lmh] * rhok_val) * rcmvb;

  for (int k = lmh - 1; k >= 1; --k) {
    real rcml = real(1.) / cm[k];
    u[k] = (-cr[k] * u[k + 1] + rsu[k]) * rcml;
    v[k] = (-cr[k] * v[k + 1] + rsv[k]) * rcml;
  }
}

// ----------------------------------------------------------------------
// Core Column Kernel: MYJPBLKernel (module_bl_myjpbl.F:131-766)
// ----------------------------------------------------------------------
struct MYJPBLKernel {
  // Atmospheric and surface inputs
  View2<const real> ht;
  View3<const real> dz;
  View3<const real> pmid;
  View3<const real> pint;
  View3<const real> th;
  View3<const real> t;
  View3<const real> exner;
  View3<const real> qv;
  View3<const real> qc;
  View3<const real> qi;
  View3<const real> u;
  View3<const real> v;
  View3<const real> rho;
  View2<const real> tsk;
  View2<const real> chklowq;
  View2<const real> thz0_in;
  View2<const real> qz0_in;
  View2<const real> uz0;
  View2<const real> vz0;
  View2<const int> lowlyr;
  View2<const real> xland;
  View2<const real> sice;
  View2<const real> snow;
  View2<const real> ustar;
  View2<const real> znt;
  View2<const real> akhs;
  View2<const real> akms;
  View2<const real> elflx;

  // Inouts
  View2<real> qsfc;
  View2<real> thz0;
  View2<real> qz0;
  View2<real> ct;
  View3<real> tke_myj;
  View3<real> exch_h;

  // Outputs
  View3<real> el_myj;
  View2<real> pblh;
  View2<int> kpbl;
  View2<real> mixht;
  View3<real> rublten;
  View3<real> rvblten;
  View3<real> rthblten;
  View3<real> rqvblten;
  View3<real> rqcblten;
  View3<real> rqiblten;

  real dt;
  int stepbl;
  int nz;
  bool has_qi;
  int i0;
  int nx;

  CPPWRF_HD void operator()(int i, int j) const {
    real dtturbl = dt * real(stepbl);
    real rdtturbl = real(1.) / dtturbl;
    real dtdif = dtturbl;

    int lowlyr_val = lowlyr(i, j);
    int lmh = nz - lowlyr_val + 1;

    // Local column arrays (1-indexed matching Fortran: 1 to nz)
    real tk[kMaxLevels];
    real thek[kMaxLevels];
    real qk[kMaxLevels];
    real cwmk[kMaxLevels];
    real qcwk[kMaxLevels];
    real qcik[kMaxLevels];
    real uk[kMaxLevels];
    real vk[kMaxLevels];
    real q2k[kMaxLevels];
    real zhk[kMaxLevels + 1];

    // Tridiagonal and level-2.5 workspace buffers (shared across MIXLEN, VDIFQ, VDIFH, VDIFV)
    real cm[kMaxLevels];
    real cr[kMaxLevels];
    real dtoz[kMaxLevels];

    real gm[kMaxLevels];
    real gh[kMaxLevels];
    real el[kMaxLevels];
    real rhok[kMaxLevels];
    real akhk[kMaxLevels];
    real akmk[kMaxLevels];

    // Interface heights (module_bl_myjpbl.F:312-328, 381-384)
    zhk[nz + 1] = ht(i, j);
    for (int k_f = nz; k_f >= 1; --k_f) {
      int k_wrf = nz - k_f;
      zhk[k_f] = zhk[k_f + 1] + dz(i, j, k_wrf);
    }

    // Atmospheric profiles (module_bl_myjpbl.F:328-337, 364-377)
    for (int k_f = 1; k_f <= nz; ++k_f) {
      int k_wrf = nz - k_f;
      tk[k_f] = t(i, j, k_wrf);
      real ratiomx = qv(i, j, k_wrf);
      qk[k_f] = ratiomx / (real(1.) + ratiomx);
      qcwk[k_f] = qc(i, j, k_wrf);
      qcik[k_f] = (has_qi && qi.p) ? qi(i, j, k_wrf) : real(0.);
      cwmk[k_f] = qcwk[k_f];
      if (has_qi) cwmk[k_f] = cwmk[k_f] + qcik[k_f];
      thek[k_f] = (cwmk[k_f] * (-ELOCP / tk[k_f]) + real(1.)) * th(i, j, k_wrf);
      uk[k_f] = u(i, j, k_wrf);
      vk[k_f] = v(i, j, k_wrf);
      q2k[k_f] = real(2.) * tke_myj(i, j, k_wrf);
    }

    int lpbl = lmh;
    int lmxl = lmh;
    real pblh_val = real(0.);
    real mixht_val = real(0.);
    real ct_val = ct(i, j);

    // 1. Mixing length (module_bl_myjpbl.F:399-404; zeros CT(I,J) at line 852)
    mixlen(lmh, uk, vk, tk, thek, qk, cwmk, q2k, zhk, gm, gh, el,
           pblh_val, lpbl, lmxl, ct_val, mixht_val, cm, cr, dtoz, rhok);
    ct(i, j) = ct_val;

    pblh(i, j) = pblh_val;
    mixht(i, j) = mixht_val;
    kpbl(i, j) = nz - lpbl + 1;

    // 2. Production/dissipation of TKE (module_bl_myjpbl.F:412-415)
    real ustar_val = ustar(i, j);
    prodq2(lmh, dtturbl, ustar_val, gm, gh, el, q2k);

    // 3. Exchange coefficients in free atmosphere (module_bl_myjpbl.F:427-444)
    difcof(lmh, gm, gh, el, q2k, zhk, akmk, akhk);

    // Save exch_h: defined on bottoms of layers 0 to nz-2 (upward)
    for (int k_wrf = 0; k_wrf < nz - 1; ++k_wrf) {
      int k_f = k_wrf + 1;
      int k_flip = nz - k_f;  // Fortran KFLIP = KTE - K
      real deltaz = real(0.5) * (zhk[k_flip] - zhk[k_flip + 2]);
      exch_h(i, j, k_wrf) = akhk[k_flip] * deltaz;
    }

    // 4. Vertical diffusion of TKE (module_bl_myjpbl.F:452-464; gm/gh reused as akq/rsq2)
    vdifq(lmh, dtdif, q2k, el, zhk, gm, gh, cm, cr, dtoz);

    for (int k_wrf = 0; k_wrf < nz; ++k_wrf) {
      int k_f = nz - k_wrf;
      q2k[k_f] = c_max(q2k[k_f], EPSQ2);
      tke_myj(i, j, k_wrf) = real(0.5) * q2k[k_f];
      if (k_f < nz) {
        el_myj(i, j, k_wrf) = el[k_f];
      } else {
        el_myj(i, j, k_wrf) = real(0.);
      }
    }

    // 5. Surface preparations and boundary conditions (module_bl_myjpbl.F:474-585)
    // In WRF's module_bl_myjpbl.o, the THSK loop (lines 474-479) is vectorized with VF=4
    // (because LOWLYR is 32-bit int and THSK is 64-bit real*8), calling _ZGVbN2vv_pow
    // for the first (nx & ~3) columns and scalar pow for the remaining (nx & 3) columns.
    real psfc = pint(i, j, lowlyr_val - 1);
    int ii = i - i0;
    real thsk_val = tsk(i, j) * cppwrf::vm::pow(real(1.e5) / psfc, CAPA, ii, (ii < (nx & ~3)) ? 2 : 1);

    for (int k_f = 1; k_f <= nz; ++k_f) {
      int k_wrf = nz - k_f;
      rhok[k_f] = pmid(i, j, k_wrf) / (R_D * tk[k_f] * (real(1.) + P608 * qk[k_f] - cwmk[k_f]));
    }

    for (int k = 1; k <= nz - 1; ++k) {
      akhk[k] = akhk[k] * real(0.5) * (rhok[k] + rhok[k + 1]);
      akmk[k] = akmk[k] * (rhok[k] + rhok[k + 1]) * real(0.5);
    }

    real seamask = xland(i, j) - real(1.);
    real thz0_val = thz0_in(i, j);
    thz0_val = (real(1.) - seamask) * thsk_val + seamask * thz0_val;
    thz0(i, j) = thz0_val;

    real akhs_dens = akhs(i, j) * rhok[lmh];
    real akms_dens = akms(i, j) * rhok[lmh];

    real qsfc_val = qsfc(i, j);
    if (seamask < real(0.5)) {
      real qfc1 = XLV * chklowq(i, j) * akhs_dens;
      if (snow(i, j) > real(0.) || sice(i, j) > real(0.5)) {
        qfc1 = qfc1 * RLIVWV;
      }
      if (qfc1 > real(0.)) {
        real qlow = qk[lmh];
        qsfc_val = qlow + elflx(i, j) / qfc1;
      }
    } else {
      real exnsfc = rpow(real(1.e5) / psfc, CAPA);
      qsfc_val = PQ0SEA / psfc * std::exp(A2 * (thsk_val - A3 * exnsfc) / (thsk_val - A4 * exnsfc));
    }
    qsfc(i, j) = qsfc_val;

    real qz0_val = qz0_in(i, j);
    qz0_val = (real(1.) - seamask) * qsfc_val + seamask * qz0_val;
    qz0(i, j) = qz0_val;

    // 6. Vertical diffusion of heat, moisture, condensate (module_bl_myjpbl.F:501-652; gm/gh reused as rkct/rss)
    vdifh_1spec(dtdif, lmh, lpbl, thz0_val, akhs_dens, real(1.), ct_val,
                thek, akhk, zhk, rhok, gm, gh, cm, cr, dtoz);
    vdifh_1spec(dtdif, lmh, lpbl, qz0_val, akhs_dens, chklowq(i, j), real(0.),
                qk, akhk, zhk, rhok, gm, gh, cm, cr, dtoz);
    vdifh_1spec(dtdif, lmh, lpbl, real(0.), akhs_dens, real(0.), real(0.),
                qcwk, akhk, zhk, rhok, gm, gh, cm, cr, dtoz);
    if (has_qi) {
      vdifh_1spec(dtdif, lmh, lpbl, real(0.), akhs_dens, real(0.), real(0.),
                  qcik, akhk, zhk, rhok, gm, gh, cm, cr, dtoz);
    }

    // Compute thermal and moisture tendencies
    for (int k_wrf = 0; k_wrf < nz; ++k_wrf) {
      int k_f = nz - k_wrf;
      real the_new = thek[k_f];
      real q_new = qk[k_f];
      real qcw_new = qcwk[k_f];
      real qci_new = has_qi ? qcik[k_f] : real(0.);
      real cwm_new = qcw_new;
      if (has_qi) cwm_new = cwm_new + qci_new;

      real apex = real(1.) / exner(i, j, k_wrf);
      real th_old = th(i, j, k_wrf);
      real th_new = the_new + cwm_new * ELOCP * apex;
      rthblten(i, j, k_wrf) = (th_new - th_old) * rdtturbl;

      real q_old = qv(i, j, k_wrf) / (real(1.) + qv(i, j, k_wrf));
      real dqdt = (q_new - q_old) * rdtturbl;
      rqvblten(i, j, k_wrf) = dqdt / ((real(1.) - q_new) * (real(1.) - q_new));

      rqcblten(i, j, k_wrf) = (qcw_new - qc(i, j, k_wrf)) * rdtturbl;
      if (has_qi && rqiblten.p) {
        rqiblten(i, j, k_wrf) = (qci_new - (qi.p ? qi(i, j, k_wrf) : real(0.))) * rdtturbl;
      }
    }

    // 7. Vertical diffusion of velocity components (module_bl_myjpbl.F:720-758; gm/gh reused as rsu/rsv)
    vdifv(lmh, dtdif, uz0(i, j), vz0(i, j), akms_dens,
          uk, vk, akmk, zhk, rhok, gm, gh, cm, cr, dtoz);

    for (int k_wrf = 0; k_wrf < nz; ++k_wrf) {
      int k_f = nz - k_wrf;
      rublten(i, j, k_wrf) = (uk[k_f] - u(i, j, k_wrf)) * rdtturbl;
      rvblten(i, j, k_wrf) = (vk[k_f] - v(i, j, k_wrf)) * rdtturbl;
    }
  }
};

// ----------------------------------------------------------------------
// Kernel: MYJPBLInitKernel (module_bl_myjpbl.F:1694-1736)
// ----------------------------------------------------------------------
struct MYJPBLInitKernel {
  View3<real> rublten;
  View3<real> rvblten;
  View3<real> rthblten;
  View3<real> rqvblten;
  View3<real> tke_myj;
  View3<real> exch_h;
  int nz;

  CPPWRF_HD void operator()(int i, int j) const {
    for (int k = 0; k < nz; ++k) {
      tke_myj(i, j, k) = EPSQ2;
      rublten(i, j, k) = real(0.);
      rvblten(i, j, k) = real(0.);
      rthblten(i, j, k) = real(0.);
      rqvblten(i, j, k) = real(0.);
      exch_h(i, j, k) = real(0.);
    }
  }
};

// Host launcher functions
void run_myjpbl(
    View2<const real> ht,
    View3<const real> dz,
    View3<const real> pmid,
    View3<const real> pint,
    View3<const real> th,
    View3<const real> t,
    View3<const real> exner,
    View3<const real> qv,
    View3<const real> qc,
    View3<const real> qi,
    View3<const real> u,
    View3<const real> v,
    View3<const real> rho,
    View2<const real> tsk,
    View2<const real> chklowq,
    View2<const real> thz0_in,
    View2<const real> qz0_in,
    View2<const real> uz0,
    View2<const real> vz0,
    View2<const int> lowlyr,
    View2<const real> xland,
    View2<const real> sice,
    View2<const real> snow,
    View2<const real> ustar,
    View2<const real> znt,
    View2<const real> akhs,
    View2<const real> akms,
    View2<const real> elflx,
    View2<real> qsfc,
    View2<real> thz0,
    View2<real> qz0,
    View2<real> ct,
    View3<real> tke_myj,
    View3<real> exch_h,
    View3<real> el_myj,
    View2<real> pblh,
    View2<int> kpbl,
    View2<real> mixht,
    View3<real> rublten,
    View3<real> rvblten,
    View3<real> rthblten,
    View3<real> rqvblten,
    View3<real> rqcblten,
    View3<real> rqiblten,
    const Range2D& range,
    real dt,
    int stepbl,
    int nz,
    bool has_qi);

void run_myjpblinit(
    View3<real> rublten,
    View3<real> rvblten,
    View3<real> rthblten,
    View3<real> rqvblten,
    View3<real> tke_myj,
    View3<real> exch_h,
    const Range2D& range,
    int nz);

}  // namespace cppwrf::phys::myjpbl
