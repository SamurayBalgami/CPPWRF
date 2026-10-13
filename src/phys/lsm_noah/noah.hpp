#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "phys/lsm_noah/noah_tables.hpp"
#include "phys/phys_constants.hpp"

namespace cppwrf::phys::noah {

static constexpr int kMaxSoilLevels = 8;

struct NoahKernel {
  // Atmospheric Forcing (2D views)
  View2<const real> lwdn;
  View2<const real> soldn;
  View2<const real> solnet;
  View2<const real> sfcprs;
  View2<const real> prcp;
  View2<const real> sfctmp;
  View2<const real> th2;
  View2<const real> q2;
  View2<const real> sfcspd;
  View2<const real> q2sat;
  View2<const real> dqsdt2;
  View2<const real> ffrozp;

  // Surface characteristics (2D views)
  View2<const int> vegtyp;
  View2<const int> soiltyp;
  View2<const int> slopetyp;
  View2<const real> shdfac;
  View2<const real> shdmin;
  View2<const real> shdmax;
  View2<const real> alb;
  View2<const real> snoalb;
  View2<const real> tbot;
  View2<real> z0brd;
  View2<real> embrd;

  // History / State variables (INOUT)
  View2<real> cmc;
  View2<real> t1;
  View3<real> stc;   // (i, j, k) where k is soil layer
  View3<real> smc;   // (i, j, k)
  View3<real> sh2o;  // (i, j, k)
  View2<real> snowh;
  View2<real> sneqv;
  View2<real> albedo;
  View2<real> ch;
  View2<real> cm;
  View2<real> emissi;
  View2<real> z0;
  View2<real> xlai;
  View2<real> snotime1;
  View2<real> ribb;

  // Outputs (2D views)
  View2<real> eta;
  View2<real> sheat;
  View2<real> eta_kinematic;
  View2<real> fdown;
  View2<real> ec;
  View2<real> edir;
  View3<real> et;    // (i, j, k)
  View2<real> ett;
  View2<real> esnow;
  View2<real> drip;
  View2<real> dew;
  View2<real> beta;
  View2<real> etp;
  View2<real> ssoil;
  View2<real> flx1;
  View2<real> flx2;
  View2<real> flx3;
  View2<real> snomlt;
  View2<real> sncovr;
  View2<real> runoff1;
  View2<real> runoff2;
  View2<real> runoff3;
  View2<real> rc;
  View2<real> pc;
  View2<real> rsmin;
  View2<real> soilw;
  View2<real> soilm;
  View2<real> q1;
  View3<real> smav;  // (i, j, k)

  // Configuration
  real dt;
  real zlvl;
  int nsoil;
  real sldpth[kMaxSoilLevels];
  int isurban;
  NoahTables tables;

  // WRF namelist rdlai2d / usemonalb (SFLX/REDPRM: table LAI and background albedo when false).
  // The standalone SFLX test runs with both true.
  bool rdlai2d = true;
  bool usemonalb = true;
  // Called from the Noah driver (lsm): the driver selects the points (active != 0: land, not
  // ice) and SFLX has no water test of its own; ALB is INOUT (written to alb_out).
  View2<const int> active;
  View2<real> alb_out;

  // Helper static functions mirroring module_sf_noahlsm.F
  CPPWRF_HD static void csnow(real& sncond, real dsnow) {
    constexpr real unit = real(0.11631);
    real c = real(0.328) * std::pow(real(10.0), real(2.25) * dsnow);
    sncond = real(2.0) * unit * c;
  }

  CPPWRF_HD static void snow_new(real temp, real newsn, real& snowh, real& sndens) {
    real snowhc = snowh * real(100.0);
    real newsnc = newsn * real(100.0);
    real tempc = temp - real(273.15);
    real dsnew;
    if (tempc <= real(-15.0)) {
      dsnew = real(0.05);
    } else {
      dsnew = real(0.05) + real(0.0017) * std::pow(tempc + real(15.0), real(1.5));
    }
    real hnewc = newsnc / dsnew;
    if (snowhc + hnewc < real(1.0e-3)) {
      sndens = std::max(dsnew, sndens);
    } else {
      sndens = (snowhc * sndens + hnewc * dsnew) / (snowhc + hnewc);
    }
    snowhc += hnewc;
    snowh = snowhc * real(0.01);
  }

  CPPWRF_HD static void snfrac(real sneqv, real snup, real salp, real snowh, real& sncovr) {
    if (sneqv < snup) {
      real rsnow = sneqv / snup;
      sncovr = real(1.0) - (std::exp(-salp * rsnow) - rsnow * std::exp(-salp));
    } else {
      sncovr = real(1.0);
    }
  }

  CPPWRF_HD static void alcalc(real alb, real snoalb, real embrd, real shdfac, real shdmin,
                               real sncovr, real tsnow, real& albedo, real& emissi,
                               real dt, bool snowng, real& snotime1, real lvcoef) {
    constexpr real emissi_s = real(0.95);
    constexpr real snacca = real(0.94);
    constexpr real snaccb = real(0.58);

    emissi = embrd + sncovr * (emissi_s - embrd);
    real snoalb1 = snoalb + lvcoef * (real(0.85) - snoalb);
    real snoalb2 = snoalb1;
    if (snowng) {
      snotime1 = real(0.0);
    } else {
      snotime1 += dt;
      snoalb2 = snoalb1 * std::pow(snacca, std::pow(snotime1 / real(86400.0), snaccb));
    }
    snoalb2 = std::max(snoalb2, alb);
    albedo = alb + sncovr * (snoalb2 - alb);
    if (albedo > snoalb2) albedo = snoalb2;
  }

  CPPWRF_HD static void snowz0(real sncovr, real& z0, real z0brd, real snowh) {
    constexpr real z0s = real(0.001);
    real burial = real(7.0) * z0brd - snowh;
    real z0eff;
    if (burial <= real(0.0007)) {
      z0eff = z0s;
    } else {
      z0eff = burial / real(7.0);
    }
    z0 = (real(1.0) - sncovr) * z0brd + sncovr * z0eff;
  }

  CPPWRF_HD static void tdfcnd(real& df, real smc, real qz, real smcmax, real sh2o,
                               real bexp, real psisat, int soiltyp, int opt_thcnd) {
    real satratio = smc / smcmax;
    constexpr real thkice = real(2.2);
    constexpr real thkw = real(0.57);
    constexpr real thko = real(2.0);
    constexpr real thkqtz = real(7.7);

    real thks = std::pow(thkqtz, qz) * std::pow(thko, real(1.0) - qz);
    real xunfroz = sh2o / smc;
    real xu = xunfroz * smcmax;
    real thksat = std::pow(thks, real(1.0) - smcmax) *
                  std::pow(thkice, smcmax - xu) *
                  std::pow(thkw, xu);
    real gammd = (real(1.0) - smcmax) * real(2700.0);
    real thkdry = (real(0.135) * gammd + real(64.7)) / (real(2700.0) - real(0.947) * gammd);
    real akei = satratio;
    real akel;
    if (satratio > real(0.1)) {
      akel = std::log10(satratio) + real(1.0);
    } else {
      akel = real(0.0);
    }
    real ake = ((smc - sh2o) * akei + sh2o * akel) / smc;
    df = ake * (thksat - thkdry) + thkdry;
  }

  CPPWRF_HD static void penman(real sfctmp, real sfcprs, real ch, real t2v, real th2,
                               real prcp, real fdown, real& t24, real ssoil,
                               real q2, real q2sat, real& etp, real& rch, real& epsca,
                               real& rr, bool snowng, bool frzgra, real dqsdt2,
                               real& flx2, real emissi, real sneqv, real t1, real sncovr,
                               real aoasis, real albedo, real soldn, real stc1) {
    constexpr real elcp = real(2.4888e3);
    constexpr real subc = real(2.501000e6);
    constexpr real lsubs = real(2.83e6);
    constexpr real cp_val = real(1004.6);
    constexpr real rd_val = real(287.04);
    constexpr real sigma = real(5.67e-8);
    constexpr real cph2o = real(4.218e3);
    constexpr real cpice = real(2.106e3);
    constexpr real lsubf = real(3.335e5);

    real elcp1 = (real(1.0) - sncovr) * elcp + sncovr * elcp * lsubs / subc;
    real lvs = (real(1.0) - sncovr) * subc + sncovr * lsubs;
    flx2 = real(0.0);
    real delta = elcp1 * dqsdt2;
    t24 = sfctmp * sfctmp * sfctmp * sfctmp;
    rr = emissi * t24 * real(6.48e-8) / (sfcprs * ch) + real(1.0);
    real rho = sfcprs / (rd_val * t2v);
    rch = rho * cp_val * ch;
    if (!snowng) {
      if (prcp > real(0.0)) rr += cph2o * prcp / rch;
    } else {
      rr += cpice * prcp / rch;
    }
    real fnet = fdown - emissi * sigma * t24 - ssoil;
    if (frzgra) {
      flx2 = - lsubf * prcp;
      fnet -= flx2;
    }
    real rad = fnet / rch + th2 - sfctmp;
    real a = elcp1 * (q2sat - q2);
    epsca = (a * rr + rad * delta) / (delta + rr);
    if (epsca > real(0.0)) epsca *= aoasis;
    etp = epsca * rch / lvs;
  }

  CPPWRF_HD static void canres(real solar, real ch, real sfctmp, real q2, real sfcprs,
                               const real* sh2o, const real* zsoil, int nsoil,
                               real smcwlt, real smcref, real rsmin, real& rc, real& pc,
                               int nroot, real q2sat, real dqsdt2, real topt, real rsmax,
                               real rgl, real hs, real xlai, real& rcs, real& rct,
                               real& rcq, real& rcsoil, real emissi) {
    constexpr real slv = real(2.501000e6);
    constexpr real cp_val = real(1004.5);
    constexpr real rd_val = real(287.04);
    constexpr real sigma = real(5.67e-8);

    rcs = real(0.0); rct = real(0.0); rcq = real(0.0); rcsoil = real(0.0); rc = real(0.0);
    real ff = real(0.55) * real(2.0) * solar / (rgl * xlai);
    rcs = (ff + rsmin / rsmax) / (real(1.0) + ff);
    rcs = std::max(rcs, real(0.0001));

    rct = real(1.0) - real(0.0016) * ((topt - sfctmp) * (topt - sfctmp));
    rct = std::max(rct, real(0.0001));

    rcq = real(1.0) / (real(1.0) + hs * (q2sat - q2));
    rcq = std::max(rcq, real(0.01));

    real part[kMaxSoilLevels];
    real gx = (sh2o[0] - smcwlt) / (smcref - smcwlt);
    if (gx > real(1.0)) gx = real(1.0);
    if (gx < real(0.0)) gx = real(0.0);
    part[0] = (zsoil[0] / zsoil[nroot - 1]) * gx;
    for (int k = 1; k < nroot; ++k) {
      gx = (sh2o[k] - smcwlt) / (smcref - smcwlt);
      if (gx > real(1.0)) gx = real(1.0);
      if (gx < real(0.0)) gx = real(0.0);
      part[k] = ((zsoil[k] - zsoil[k - 1]) / zsoil[nroot - 1]) * gx;
    }
    for (int k = 0; k < nroot; ++k) rcsoil += part[k];
    rcsoil = std::max(rcsoil, real(0.0001));

    rc = rsmin / (xlai * rcs * rct * rcq * rcsoil);
    real rr = (real(4.0) * emissi * sigma * rd_val / cp_val) *
              std::pow(sfctmp, real(4.0)) / (sfcprs * ch) + real(1.0);
    real delta = (slv / cp_val) * dqsdt2;
    pc = (rr + delta) / (rr * (real(1.0) + rc * ch) + delta);
  }

  CPPWRF_HD static void devap(real& edir, real etp1, real smc, real zsoil, real shdfac,
                              real smcmax, real bexp, real dksat, real dwsat,
                              real smcdry, real smcref, real smcwlt, real fxexp) {
    real sratio = (smc - smcdry) / (smcmax - smcdry);
    real fx;
    if (sratio > real(0.0)) {
      fx = std::pow(sratio, fxexp);
      fx = std::max(std::min(fx, real(1.0)), real(0.0));
    } else {
      fx = real(0.0);
    }
    edir = fx * (real(1.0) - shdfac) * etp1;
  }

  CPPWRF_HD static void transp(real* et, int nsoil, real etp1, const real* smc, real cmc,
                               const real* zsoil, real shdfac, real smcwlt, real cmcmax,
                               real pc, real cfactr, real smcref, real sfctmp, real q2,
                               int nroot, const real* rtdis) {
    for (int k = 0; k < nsoil; ++k) et[k] = real(0.0);
    real etp1a;
    if (cmc != real(0.0)) {
      etp1a = shdfac * pc * etp1 * (real(1.0) - std::pow(cmc / cmcmax, cfactr));
    } else {
      etp1a = shdfac * pc * etp1;
    }
    real gx[kMaxSoilLevels];
    real sgx = real(0.0);
    for (int i = 0; i < nroot; ++i) {
      gx[i] = (smc[i] - smcwlt) / (smcref - smcwlt);
      gx[i] = std::max(std::min(gx[i], real(1.0)), real(0.0));
      sgx += gx[i];
    }
    sgx /= real(nroot);
    real denom = real(0.0);
    for (int i = 0; i < nroot; ++i) {
      real rtx = rtdis[i] + gx[i] - sgx;
      gx[i] = gx[i] * std::max(rtx, real(0.0));
      denom += gx[i];
    }
    if (denom <= real(0.0)) denom = real(1.0);
    for (int i = 0; i < nroot; ++i) {
      et[i] = etp1a * gx[i] / denom;
    }
  }

  CPPWRF_HD static void evapo(real& eta1, const real* smc, int nsoil, real cmc, real etp1,
                              real dt, const real* zsoil, const real* sh2o, real smcmax,
                              real bexp, real pc, real smcwlt, real dksat, real dwsat,
                              real smcref, real shdfac, real cmcmax, real smcdry,
                              real cfactr, real& edir, real& ec, real* et, real& ett,
                              real sfctmp, real q2, int nroot, const real* rtdis, real fxexp) {
    edir = real(0.0);
    ec = real(0.0);
    ett = real(0.0);
    for (int k = 0; k < nsoil; ++k) et[k] = real(0.0);

    if (etp1 > real(0.0)) {
      if (shdfac < real(1.0)) {
        devap(edir, etp1, smc[0], zsoil[0], shdfac, smcmax, bexp, dksat, dwsat, smcdry, smcref, smcwlt, fxexp);
      }
      if (shdfac > real(0.0)) {
        transp(et, nsoil, etp1, sh2o, cmc, zsoil, shdfac, smcwlt, cmcmax, pc, cfactr, smcref, sfctmp, q2, nroot, rtdis);
        for (int k = 0; k < nsoil; ++k) ett += et[k];
        if (cmc > real(0.0)) {
          ec = shdfac * std::pow(cmc / cmcmax, cfactr) * etp1;
        } else {
          ec = real(0.0);
        }
        real cmc2ms = cmc / dt;
        ec = std::min(cmc2ms, ec);
      }
    }
    eta1 = edir + ett + ec;
  }

  CPPWRF_HD static void rosr12(real* p, const real* a, const real* b, real* c,
                               const real* d, real* delta, int nsoil) {
    c[nsoil - 1] = real(0.0);
    p[0] = - c[0] / b[0];
    delta[0] = d[0] / b[0];
    for (int k = 1; k < nsoil; ++k) {
      p[k] = - c[k] * (real(1.0) / (b[k] + a[k] * p[k - 1]));
      delta[k] = (d[k] - a[k] * delta[k - 1]) * (real(1.0) / (b[k] + a[k] * p[k - 1]));
    }
    p[nsoil - 1] = delta[nsoil - 1];
    for (int k = 1; k < nsoil; ++k) {
      int kk = nsoil - 1 - k;
      p[kk] = p[kk] * p[kk + 1] + delta[kk];
    }
  }

  CPPWRF_HD static void wdfcnd(real& wdf, real& wcnd, real smc, real smcmax,
                               real bexp, real dksat, real dwsat, real sicemax) {
    real factr1 = real(0.05) / smcmax;
    real factr2 = smc / smcmax;
    factr1 = std::min(factr1, factr2);
    real expon = bexp + real(2.0);
    wdf = dwsat * std::pow(factr2, expon);
    if (sicemax > real(0.0)) {
      real term = real(500.0) * sicemax;
      real vkwgt = real(1.0) / (real(1.0) + term * term * term);
      wdf = vkwgt * wdf + (real(1.0) - vkwgt) * dwsat * std::pow(factr1, expon);
    }
    expon = real(2.0) * bexp + real(3.0);
    wcnd = dksat * std::pow(factr2, expon);
  }

  CPPWRF_HD static void srt(real* rhstt, real edir, const real* et, const real* sh2o,
                            const real* sh2oa, int nsoil, real pcpdrp, const real* zsoil,
                            real dwsat, real dksat, real smcmax, real bexp,
                            real& runoff1, real& runoff2, real dt, real smcwlt,
                            real slope, real kdt, real frzx, const real* sice,
                            real* ai, real* bi, real* ci) {
    constexpr int cvfrz = 3;
    real sicemax = real(0.0);
    for (int ks = 0; ks < nsoil; ++ks) {
      if (sice[ks] > sicemax) sicemax = sice[ks];
    }
    real pddum = pcpdrp;
    runoff1 = real(0.0);

    if (pcpdrp != real(0.0)) {
      real dt1 = dt / real(86400.0);
      real smcav = smcmax - smcwlt;
      real dmax[kMaxSoilLevels];
      dmax[0] = - zsoil[0] * smcav;
      real dice = - zsoil[0] * sice[0];
      dmax[0] = dmax[0] * (real(1.0) - (sh2oa[0] + sice[0] - smcwlt) / smcav);
      real dd = dmax[0];
      for (int ks = 1; ks < nsoil; ++ks) {
        dice += (zsoil[ks - 1] - zsoil[ks]) * sice[ks];
        dmax[ks] = (zsoil[ks - 1] - zsoil[ks]) * smcav;
        dmax[ks] = dmax[ks] * (real(1.0) - (sh2oa[ks] + sice[ks] - smcwlt) / smcav);
        dd += dmax[ks];
      }
      real val = real(1.0) - std::exp(- kdt * dt1);
      real ddt = dd * val;
      real px = pcpdrp * dt;
      if (px < real(0.0)) px = real(0.0);
      real infmax = (px * (ddt / (px + ddt))) / dt;
      real fcr = real(1.0);
      if (dice > real(1.0e-2)) {
        real acrt = real(cvfrz) * frzx / dice;
        real sum = real(1.0);
        int ialp1 = cvfrz - 1;
        for (int j = 1; j <= ialp1; ++j) {
          int k = 1;
          for (int jj = j + 1; jj <= ialp1; ++jj) k *= jj;
          sum += std::pow(acrt, real(cvfrz - j)) / real(k);
        }
        fcr = real(1.0) - std::exp(- acrt) * sum;
      }
      infmax *= fcr;
      real mxsmc = sh2oa[0];
      real wdf, wcnd;
      wdfcnd(wdf, wcnd, mxsmc, smcmax, bexp, dksat, dwsat, sicemax);
      infmax = std::max(infmax, wcnd);
      infmax = std::min(infmax, px / dt);
      if (pcpdrp > infmax) {
        runoff1 = pcpdrp - infmax;
        pddum = infmax;
      }
    }

    real mxsmc = sh2oa[0];
    real wdf, wcnd;
    wdfcnd(wdf, wcnd, mxsmc, smcmax, bexp, dksat, dwsat, sicemax);
    real ddz = real(1.0) / (- real(0.5) * zsoil[1]);
    ai[0] = real(0.0);
    bi[0] = wdf * ddz / (- zsoil[0]);
    ci[0] = - bi[0];
    real dsmdz = (sh2o[0] - sh2o[1]) / (- real(0.5) * zsoil[1]);
    rhstt[0] = (wdf * dsmdz + wcnd - pddum + edir + et[0]) / zsoil[0];

    real ddz2 = real(0.0);
    for (int k = 1; k < nsoil; ++k) {
      real denom2 = zsoil[k - 1] - zsoil[k];
      real slopx;
      real wdf2, wcnd2, dsmdz2;
      if (k != nsoil - 1) {
        slopx = real(1.0);
        real mxsmc2 = sh2oa[k];
        wdfcnd(wdf2, wcnd2, mxsmc2, smcmax, bexp, dksat, dwsat, sicemax);
        real denom = zsoil[k - 1] - zsoil[k + 1];
        dsmdz2 = (sh2o[k] - sh2o[k + 1]) / (denom * real(0.5));
        ddz2 = real(2.0) / denom;
        ci[k] = - wdf2 * ddz2 / denom2;
      } else {
        slopx = slope;
        wdfcnd(wdf2, wcnd2, sh2oa[nsoil - 1], smcmax, bexp, dksat, dwsat, sicemax);
        dsmdz2 = real(0.0);
        ci[k] = real(0.0);
      }
      real numer = (wdf2 * dsmdz2) + slopx * wcnd2 - (wdf * dsmdz) - wcnd + et[k];
      rhstt[k] = numer / (- denom2);
      ai[k] = - wdf * ddz / denom2;
      bi[k] = - (ai[k] + ci[k]);
      if (k == nsoil - 1) {
        runoff2 = slopx * wcnd2;
      } else {
        wdf = wdf2;
        wcnd = wcnd2;
        dsmdz = dsmdz2;
        ddz = ddz2;
      }
    }
  }

  CPPWRF_HD static void sstep(real* sh2oout, const real* sh2oin, real& cmc, real* rhstt,
                             real rhsct, real dt, int nsoil, real smcmax, real cmcmax,
                             real& runoff3, const real* zsoil, real* smc, const real* sice,
                             real* ai, real* bi, real* ci) {
    real rhsttin[kMaxSoilLevels];
    real ciin[kMaxSoilLevels];
    for (int k = 0; k < nsoil; ++k) {
      rhstt[k] *= dt;
      ai[k] *= dt;
      bi[k] = real(1.0) + bi[k] * dt;
      ci[k] *= dt;
      rhsttin[k] = rhstt[k];
      ciin[k] = ci[k];
    }
    rosr12(ci, ai, bi, ciin, rhsttin, rhstt, nsoil);
    real wplus = real(0.0);
    runoff3 = real(0.0);
    real ddz = - zsoil[0];
    for (int k = 0; k < nsoil; ++k) {
      if (k != 0) ddz = zsoil[k - 1] - zsoil[k];
      sh2oout[k] = sh2oin[k] + ci[k] + wplus / ddz;
      real stot = sh2oout[k] + sice[k];
      if (stot > smcmax) {
        if (k == 0) {
          ddz = - zsoil[0];
        } else {
          ddz = - zsoil[k] + zsoil[k - 1];
        }
        wplus = (stot - smcmax) * ddz;
      } else {
        wplus = real(0.0);
      }
      smc[k] = std::max(std::min(stot, smcmax), real(0.02));
      sh2oout[k] = std::max(smc[k] - sice[k], real(0.0));
    }
    runoff3 = wplus;
    cmc += dt * rhsct;
    if (cmc < real(1.0e-20)) cmc = real(0.0);
    cmc = std::min(cmc, cmcmax);
  }

  CPPWRF_HD static void smflx(real* smc, int nsoil, real& cmc, real dt, real prcp1,
                              const real* zsoil, real* sh2o, real slope, real kdt,
                              real frzfact, real smcmax, real bexp, real smcwlt,
                              real dksat, real dwsat, real shdfac, real cmcmax,
                              real& runoff1, real& runoff2, real& runoff3,
                              real edir, real ec, const real* et, real& drip) {
    real rhsct = shdfac * prcp1 - ec;
    drip = real(0.0);
    real trhsct = dt * rhsct;
    real excess = cmc + trhsct;
    if (excess > cmcmax) drip = excess - cmcmax;
    real pcpdrp = (real(1.0) - shdfac) * prcp1 + drip / dt;

    real sice[kMaxSoilLevels];
    for (int i = 0; i < nsoil; ++i) sice[i] = smc[i] - sh2o[i];

    real fac2 = real(0.0);
    for (int i = 0; i < nsoil; ++i) fac2 = std::max(fac2, sh2o[i] / smcmax);
    real flimit;
    fac2mit(smcmax, flimit);

    real ai[kMaxSoilLevels], bi[kMaxSoilLevels], ci[kMaxSoilLevels];
    real rhstt[kMaxSoilLevels];
    real dummy = real(0.0);

    if (((pcpdrp * dt) > (real(0.0001) * real(1000.0) * (- zsoil[0]) * smcmax)) || (fac2 > flimit)) {
      real sh2ofg[kMaxSoilLevels], sh2oa[kMaxSoilLevels], sh2oin[kMaxSoilLevels];
      srt(rhstt, edir, et, sh2o, sh2o, nsoil, pcpdrp, zsoil, dwsat, dksat, smcmax, bexp,
          runoff1, runoff2, dt, smcwlt, slope, kdt, frzfact, sice, ai, bi, ci);
      sstep(sh2ofg, sh2o, dummy, rhstt, rhsct, dt, nsoil, smcmax, cmcmax, runoff3, zsoil, smc, sice, ai, bi, ci);
      for (int k = 0; k < nsoil; ++k) sh2oa[k] = (sh2o[k] + sh2ofg[k]) * real(0.5);
      srt(rhstt, edir, et, sh2o, sh2oa, nsoil, pcpdrp, zsoil, dwsat, dksat, smcmax, bexp,
          runoff1, runoff2, dt, smcwlt, slope, kdt, frzfact, sice, ai, bi, ci);
      for (int k = 0; k < nsoil; ++k) sh2oin[k] = sh2o[k];
      sstep(sh2o, sh2oin, cmc, rhstt, rhsct, dt, nsoil, smcmax, cmcmax, runoff3, zsoil, smc, sice, ai, bi, ci);
    } else {
      real sh2oin[kMaxSoilLevels];
      srt(rhstt, edir, et, sh2o, sh2o, nsoil, pcpdrp, zsoil, dwsat, dksat, smcmax, bexp,
          runoff1, runoff2, dt, smcwlt, slope, kdt, frzfact, sice, ai, bi, ci);
      for (int k = 0; k < nsoil; ++k) sh2oin[k] = sh2o[k];
      sstep(sh2o, sh2oin, cmc, rhstt, rhsct, dt, nsoil, smcmax, cmcmax, runoff3, zsoil, smc, sice, ai, bi, ci);
    }
  }

  CPPWRF_HD static void fac2mit(real smcmax, real& flimit) {
    flimit = real(0.90);
    if (smcmax == real(0.395)) flimit = real(0.59);
    else if (smcmax == real(0.434) || smcmax == real(0.404)) flimit = real(0.85);
    else if (smcmax == real(0.465) || smcmax == real(0.406)) flimit = real(0.86);
    else if (smcmax == real(0.476) || smcmax == real(0.439)) flimit = real(0.74);
    else if (smcmax == real(0.200) || smcmax == real(0.464)) flimit = real(0.80);
  }

  CPPWRF_HD static void frh2o(real& free_w, real tkelv, real smc, real sh2o,
                              real smcmax, real bexp, real psis) {
    constexpr real ck = real(8.0);
    constexpr real blim = real(5.5);
    constexpr real error_tol = real(0.005);
    constexpr real hlice = real(3.335e5);
    constexpr real gs = real(9.81);
    constexpr real t0 = real(273.15);

    real bx = bexp;
    if (bexp > blim) bx = blim;
    int nlog = 0;
    int kcount = 0;
    if (tkelv > (t0 - real(1.0e-3))) {
      free_w = smc;
    } else {
      real swl = smc - sh2o;
      if (swl > (smc - real(0.02))) swl = smc - real(0.02);
      if (swl < real(0.0)) swl = real(0.0);
      while (nlog < 10 && kcount == 0) {
        nlog++;
        real term1 = (psis * gs / hlice) * ((real(1.0) + ck * swl) * (real(1.0) + ck * swl)) *
                     std::pow(smcmax / (smc - swl), bx);
        real term2 = - (tkelv - t0) / tkelv;
        real df = std::log(term1) - std::log(term2);
        real denom = real(2.0) * ck / (real(1.0) + ck * swl) + bx / (smc - swl);
        real swlk = swl - df / denom;
        if (swlk > (smc - real(0.02))) swlk = smc - real(0.02);
        if (swlk < real(0.0)) swlk = real(0.0);
        real dswl = std::abs(swlk - swl);
        swl = swlk;
        if (dswl <= error_tol) kcount++;
      }
      if (kcount != 0) {
        free_w = smc - swl;
      } else {
        real fk = std::pow((hlice / (gs * (- psis))) * ((tkelv - t0) / tkelv), - real(1.0) / bx) * smcmax;
        if (fk < real(0.02)) fk = real(0.02);
        free_w = std::min(fk, smc);
      }
    }
  }

  CPPWRF_HD static void tbnd(real tu, real tb, const real* zsoil, real zbot, int k, int nsoil, real& tbnd1) {
    real zup = (k == 0) ? real(0.0) : zsoil[k - 1];
    real zb = (k == nsoil - 1) ? (real(2.0) * zbot - zsoil[k]) : zsoil[k + 1];
    tbnd1 = tu + (tb - tu) * (zup - zsoil[k]) / (zup - zb);
  }

  CPPWRF_HD static void tmpavg(real& tavg, real tup, real tm, real tdn, const real* zsoil, int nsoil, int k) {
    constexpr real t0 = real(273.15);
    real dz = (k == 0) ? (- zsoil[0]) : (zsoil[k - 1] - zsoil[k]);
    real dzh = dz * real(0.5);
    if (tup < t0) {
      if (tm < t0) {
        if (tdn < t0) {
          tavg = (tup + real(2.0) * tm + tdn) / real(4.0);
        } else {
          real x0 = (t0 - tm) * dzh / (tdn - tm);
          tavg = real(0.5) * (tup * dzh + tm * (dzh + x0) + t0 * (real(2.0) * dzh - x0)) / dz;
        }
      } else {
        if (tdn < t0) {
          real xup = (t0 - tup) * dzh / (tm - tup);
          real xdn = dzh - (t0 - tm) * dzh / (tdn - tm);
          tavg = real(0.5) * (tup * xup + t0 * (real(2.0) * dz - xup - xdn) + tdn * xdn) / dz;
        } else {
          real xup = (t0 - tup) * dzh / (tm - tup);
          tavg = real(0.5) * (tup * xup + t0 * (real(2.0) * dz - xup)) / dz;
        }
      }
    } else {
      if (tm < t0) {
        if (tdn < t0) {
          real xup = dzh - (t0 - tup) * dzh / (tm - tup);
          tavg = real(0.5) * (t0 * (dz - xup) + tm * (dzh + xup) + tdn * dzh) / dz;
        } else {
          real xup = dzh - (t0 - tup) * dzh / (tm - tup);
          real xdn = (t0 - tm) * dzh / (tdn - tm);
          tavg = real(0.5) * (t0 * (real(2.0) * dz - xup - xdn) + tm * (xup + xdn)) / dz;
        }
      } else {
        if (tdn < t0) {
          real xdn = dzh - (t0 - tm) * dzh / (tdn - tm);
          tavg = (t0 * (dz - xdn) + real(0.5) * (t0 + tdn) * xdn) / dz;
        } else {
          tavg = (tup + real(2.0) * tm + tdn) / real(4.0);
        }
      }
    }
  }

  CPPWRF_HD static void snksrc(real& tsnsr, real tavg, real smc, real& sh2o, const real* zsoil,
                               int nsoil, real smcmax, real psisat, real bexp, real dt,
                               int k, real qtot) {
    constexpr real dh2o = real(1000.0);
    constexpr real hlice = real(3.3350e5);
    real dz = (k == 0) ? (- zsoil[0]) : (zsoil[k - 1] - zsoil[k]);
    real free_w;
    frh2o(free_w, tavg, smc, sh2o, smcmax, bexp, psisat);
    real xh2o = sh2o + qtot * dt / (dh2o * hlice * dz);
    if (xh2o < sh2o && xh2o < free_w) {
      if (free_w > sh2o) xh2o = sh2o;
      else xh2o = free_w;
    }
    if (xh2o > sh2o && xh2o > free_w) {
      if (free_w < sh2o) xh2o = sh2o;
      else xh2o = free_w;
    }
    if (xh2o < real(0.0)) xh2o = real(0.0);
    if (xh2o > smc) xh2o = smc;
    tsnsr = - dh2o * hlice * dz * (xh2o - sh2o) / dt;
    sh2o = xh2o;
  }

  CPPWRF_HD static void hrt(real* rhsts, const real* stc, const real* smc, real smcmax,
                            int nsoil, const real* zsoil, real yy, real zz1, real tbot,
                            real zbot, real psisat, real* sh2o, real dt, real bexp,
                            int soiltyp, int opt_thcnd, real f1, real df1, real quartz,
                            real csoil, real* ai, real* bi, real* ci, int vegtyp, int isurban) {
    constexpr real t0 = real(273.15);
    constexpr real cair = real(1004.0);
    constexpr real cice = real(2.106e6);
    constexpr real ch2o = real(4.2e6);

    real csoil_loc = (vegtyp == isurban) ? real(3.0e6) : csoil;
    real hcpct = sh2o[0] * ch2o + (real(1.0) - smcmax) * csoil_loc + (smcmax - smc[0]) * cair +
                 (smc[0] - sh2o[0]) * cice;

    real ddz = real(1.0) / (- real(0.5) * zsoil[1]);
    ai[0] = real(0.0);
    ci[0] = (df1 * ddz) / (zsoil[0] * hcpct);
    bi[0] = - ci[0] + df1 / (real(0.5) * zsoil[0] * zsoil[0] * hcpct * zz1);
    real dtsdz = (stc[0] - stc[1]) / (- real(0.5) * zsoil[1]);
    real ssoil_h = df1 * (stc[0] - yy) / (real(0.5) * zsoil[0] * zz1);
    real denom = zsoil[0] * hcpct;
    rhsts[0] = (df1 * dtsdz - ssoil_h) / denom;
    real qtot = - real(1.0) * rhsts[0] * denom;

    real sice = smc[0] - sh2o[0];
    real tsurf = (yy + (zz1 - real(1.0)) * stc[0]) / zz1;
    real tbk;
    tbnd(stc[0], stc[1], zsoil, zbot, 0, nsoil, tbk);
    if ((sice > real(0.0)) || (stc[0] < t0) || (tsurf < t0) || (tbk < t0)) {
      real tavg, tsnsr;
      tmpavg(tavg, tsurf, stc[0], tbk, zsoil, nsoil, 0);
      snksrc(tsnsr, tavg, smc[0], sh2o[0], zsoil, nsoil, smcmax, psisat, bexp, dt, 0, qtot);
      rhsts[0] -= tsnsr / denom;
    }

    real ddz2 = real(0.0);
    real df1k = df1;
    for (int k = 1; k < nsoil; ++k) {
      hcpct = sh2o[k] * ch2o + (real(1.0) - smcmax) * csoil_loc + (smcmax - smc[k]) * cair +
              (smc[k] - sh2o[k]) * cice;
      real df1n;
      tdfcnd(df1n, smc[k], quartz, smcmax, sh2o[k], bexp, psisat, soiltyp, opt_thcnd);
      if (vegtyp == isurban) df1n = real(3.24);

      real dtsdz2, tbk1;
      if (k != nsoil - 1) {
        real denom_k = real(0.5) * (zsoil[k - 1] - zsoil[k + 1]);
        dtsdz2 = (stc[k] - stc[k + 1]) / denom_k;
        ddz2 = real(2.0) / (zsoil[k - 1] - zsoil[k + 1]);
        ci[k] = - df1n * ddz2 / ((zsoil[k - 1] - zsoil[k]) * hcpct);
        tbnd(stc[k], stc[k + 1], zsoil, zbot, k, nsoil, tbk1);
      } else {
        real denom_k = real(0.5) * (zsoil[k - 1] + zsoil[k]) - zbot;
        dtsdz2 = (stc[k] - tbot) / denom_k;
        ci[k] = real(0.0);
        tbnd(stc[k], tbot, zsoil, zbot, k, nsoil, tbk1);
      }
      denom = (zsoil[k] - zsoil[k - 1]) * hcpct;
      rhsts[k] = (df1n * dtsdz2 - df1k * dtsdz) / denom;
      qtot = - real(1.0) * denom * rhsts[k];

      sice = smc[k] - sh2o[k];
      real tavg, tsnsr;
      tmpavg(tavg, tbk, stc[k], tbk1, zsoil, nsoil, k);
      if ((sice > real(0.0)) || (stc[k] < t0) || (tbk < t0) || (tbk1 < t0)) {
        snksrc(tsnsr, tavg, smc[k], sh2o[k], zsoil, nsoil, smcmax, psisat, bexp, dt, k, qtot);
        rhsts[k] -= tsnsr / denom;
      }
      ai[k] = - df1k * ddz / ((zsoil[k - 1] - zsoil[k]) * hcpct);
      bi[k] = - (ai[k] + ci[k]);
      tbk = tbk1;
      df1k = df1n;
      dtsdz = dtsdz2;
      ddz = ddz2;
    }
  }

  CPPWRF_HD static void hstep(real* stcout, const real* stcin, real* rhsts,
                              real dt, int nsoil, real* ai, real* bi, real* ci) {
    real rhstsin[kMaxSoilLevels];
    real ciin[kMaxSoilLevels];
    for (int k = 0; k < nsoil; ++k) {
      rhsts[k] *= dt;
      ai[k] *= dt;
      bi[k] = real(1.0) + bi[k] * dt;
      ci[k] *= dt;
      rhstsin[k] = rhsts[k];
      ciin[k] = ci[k];
    }
    rosr12(ci, ai, bi, ciin, rhstsin, rhsts, nsoil);
    for (int k = 0; k < nsoil; ++k) {
      stcout[k] = stcin[k] + ci[k];
    }
  }

  CPPWRF_HD static void shflx(real& ssoil, real* stc, const real* smc, real smcmax,
                              int nsoil, real& t1, real dt, real yy, real zz1,
                              const real* zsoil, real tbot, real zbot, real smcwlt,
                              real psisat, real* sh2o, real bexp, real f1, real df1,
                              real quartz, real csoil, int vegtyp, int isurban,
                              int soiltyp, int opt_thcnd) {
    real ai[kMaxSoilLevels], bi[kMaxSoilLevels], ci[kMaxSoilLevels];
    real stcf[kMaxSoilLevels], rhsts[kMaxSoilLevels];

    hrt(rhsts, stc, smc, smcmax, nsoil, zsoil, yy, zz1, tbot, zbot, psisat, sh2o, dt, bexp,
        soiltyp, opt_thcnd, f1, df1, quartz, csoil, ai, bi, ci, vegtyp, isurban);

    hstep(stcf, stc, rhsts, dt, nsoil, ai, bi, ci);
    for (int i = 0; i < nsoil; ++i) stc[i] = stcf[i];

    t1 = (yy + (zz1 - real(1.0)) * stc[0]) / zz1;
    ssoil = df1 * (stc[0] - t1) / (real(0.5) * zsoil[0]);
  }

  CPPWRF_HD static void snowpack(real esd, real dtsec, real& snowh, real& sndens,
                                real tsnow, real tsoil, real snomlt) {
    constexpr real c1 = real(0.01);
    constexpr real c2 = real(21.0);
    real snowhc = snowh * real(100.0);
    real esdc = esd * real(100.0);
    real dthr = dtsec / real(3600.0);
    real tsnowc = tsnow - real(273.15);
    real tsoilc = tsoil - real(273.15);
    real tavgc = real(0.5) * (tsnowc + tsoilc);
    real esdcx = (esdc > real(1.0e-2)) ? esdc : real(1.0e-2);
    real bfac = dthr * c1 * std::exp(real(0.08) * tavgc - c2 * sndens);
    real pexp = real(0.0);
    for (int j = 4; j >= 1; --j) {
      pexp = (real(1.0) + pexp) * bfac * esdcx / real(j + 1);
    }
    pexp += real(1.0);
    real dsx = sndens * pexp;
    if (dsx > real(0.40)) dsx = real(0.40);
    if (dsx < real(0.05)) dsx = real(0.05);
    sndens = dsx;
    if (tsnowc >= real(0.0)) {
      real dw = real(0.13) * dthr / real(24.0);
      sndens = sndens * (real(1.0) - dw) + dw;
      if (sndens >= real(0.40)) sndens = real(0.40);
    }
    snowhc = esdc / sndens;
    snowh = snowhc * real(0.01);
  }

  CPPWRF_HD static void nopac(real etp, real& eta, real prcp, real* smc, real smcmax,
                              real smcwlt, real smcref, real smcdry, real& cmc, real cmcmax,
                              int nsoil, real dt, real shdfac, real sbeta, real q2,
                              real& t1, real sfctmp, real t24, real th2, real fdown,
                              real f1, real emissi, real& ssoil, real* stc, real epsca,
                              real bexp, real pc, real rch, real rr, real cfactr,
                              real* sh2o, real slope, real kdt, real frzx, real psisat,
                              const real* zsoil, real dksat, real dwsat, real tbot, real zbot,
                              real& runoff1, real& runoff2, real& runoff3, real& edir,
                              real& ec, real* et, real& ett, int nroot, const real* rtdis,
                              real quartz, real fxexp, real csoil, real& beta, real& drip,
                              real& dew, real& flx1, real& flx3, int vegtyp, int isurban,
                              int soiltyp, int opt_thcnd) {
    constexpr real sigma = real(5.67e-8);
    constexpr real cph2o = real(4.218e3);
    real prcp1 = prcp * real(0.001);
    real etp1 = etp * real(0.001);
    dew = real(0.0);
    real edir1 = real(0.0);
    real ec1 = real(0.0);
    real et1[kMaxSoilLevels] = {0};
    real ett1 = real(0.0);
    real eta1 = real(0.0);

    if (etp > real(0.0)) {
      evapo(eta1, smc, nsoil, cmc, etp1, dt, zsoil, sh2o, smcmax, bexp, pc, smcwlt, dksat, dwsat,
            smcref, shdfac, cmcmax, smcdry, cfactr, edir1, ec1, et1, ett1, sfctmp, q2, nroot, rtdis, fxexp);
      smflx(smc, nsoil, cmc, dt, prcp1, zsoil, sh2o, slope, kdt, frzx, smcmax, bexp, smcwlt,
            dksat, dwsat, shdfac, cmcmax, runoff1, runoff2, runoff3, edir1, ec1, et1, drip);
      eta = eta1 * real(1000.0);
    } else {
      dew = - etp1;
      prcp1 += dew;
      smflx(smc, nsoil, cmc, dt, prcp1, zsoil, sh2o, slope, kdt, frzx, smcmax, bexp, smcwlt,
            dksat, dwsat, shdfac, cmcmax, runoff1, runoff2, runoff3, edir1, ec1, et1, drip);
    }

    if (etp <= real(0.0)) {
      beta = (etp < real(0.0)) ? real(1.0) : real(0.0);
      eta = etp;
    } else {
      beta = eta / etp;
    }

    edir = edir1 * real(1000.0);
    ec = ec1 * real(1000.0);
    for (int k = 0; k < nsoil; ++k) et[k] = et1[k] * real(1000.0);
    ett = ett1 * real(1000.0);

    real df1;
    tdfcnd(df1, smc[0], quartz, smcmax, sh2o[0], bexp, psisat, soiltyp, opt_thcnd);
    if (vegtyp == isurban) df1 = real(3.24);
    df1 = df1 * std::exp(sbeta * shdfac);

    real yynum = fdown - emissi * sigma * t24;
    real yy = sfctmp + (yynum / rch + th2 - sfctmp - beta * epsca) / rr;
    real zz1 = df1 / (- real(0.5) * zsoil[0] * rch * rr) + real(1.0);

    shflx(ssoil, stc, smc, smcmax, nsoil, t1, dt, yy, zz1, zsoil, tbot, zbot, smcwlt,
          psisat, sh2o, bexp, f1, df1, quartz, csoil, vegtyp, isurban, soiltyp, opt_thcnd);

    flx1 = cph2o * prcp * (t1 - sfctmp);
    flx3 = real(0.0);
  }

  CPPWRF_HD static void snopac(real& etp, real& eta, real prcp, real prcpf, bool snowng,
                               real* smc, real smcmax, real smcwlt, real smcref, real smcdry,
                               real& cmc, real cmcmax, int nsoil, real dt, real sbeta,
                               real df1, real q2, real& t1, real sfctmp, real t24, real th2,
                               real fdown, real f1, real& ssoil, real* stc, real epsca,
                               real sfcprs, real bexp, real pc, real rch, real rr,
                               real cfactr, real& sncovr, real& esd, real& sndens,
                               real& snowh, real* sh2o, real slope, real kdt, real frzfact,
                               real psisat, const real* zsoil, real dwsat, real dksat,
                               real tbot, real zbot, real shdfac, real& runoff1, real& runoff2,
                               real& runoff3, real& edir, real& ec, real* et, real& ett,
                               int nroot, real& snomlt, const real* rtdis, real quartz,
                               real fxexp, real csoil, real& beta, real& drip, real& dew,
                               real& flx1, real& flx2, real& flx3, real& esnow, real& etns,
                               real emissi, real& ribb, real soldn, int isurban, int vegtyp,
                               int soiltyp, int opt_thcnd) {
    constexpr real esdmin = real(1.0e-6);
    constexpr real subc = real(2.501000e6);
    constexpr real lsubs = real(2.83e6);
    constexpr real tfreez = real(273.15);
    constexpr real snoexp = real(2.0);
    constexpr real sigma = real(5.67e-8);
    constexpr real cpice = real(2.106e3);
    constexpr real cph2o = real(4.218e3);
    constexpr real lsubf = real(3.335e5);

    dew = real(0.0);
    edir = real(0.0);
    real edir1 = real(0.0);
    ec = real(0.0);
    real ec1 = real(0.0);
    ett = real(0.0);
    real ett1 = real(0.0);
    for (int k = 0; k < nsoil; ++k) {
      et[k] = real(0.0);
    }
    real et1[kMaxSoilLevels] = {0};
    etns = real(0.0);
    real etns1 = real(0.0);
    esnow = real(0.0);
    real esnow1 = real(0.0);
    real esnow2 = real(0.0);
    real prcp1 = prcpf * real(0.001);
    beta = real(1.0);
    real etanrg = real(0.0);

    if (etp <= real(0.0)) {
      if (ribb >= real(0.1) && fdown > real(150.0)) {
        etp = (std::min(etp * (real(1.0) - ribb), real(0.0)) * sncovr / real(0.980) +
               etp * (real(0.980) - sncovr)) / real(0.980);
      }
      if (etp == real(0.0)) beta = real(0.0);
      real etp1 = etp * real(0.001);
      dew = - etp1;
      esnow2 = etp1 * dt;
      etanrg = etp * ((real(1.0) - sncovr) * subc + sncovr * lsubs);
    } else {
      real etp1 = etp * real(0.001);
      if (sncovr < real(1.0)) {
        evapo(etns1, smc, nsoil, cmc, etp1, dt, zsoil, sh2o, smcmax, bexp, pc, smcwlt,
              dksat, dwsat, smcref, shdfac, cmcmax, smcdry, cfactr, edir1, ec1, et1,
              ett1, sfctmp, q2, nroot, rtdis, fxexp);
        edir1 *= (real(1.0) - sncovr);
        ec1 *= (real(1.0) - sncovr);
        for (int k = 0; k < nsoil; ++k) et1[k] *= (real(1.0) - sncovr);
        ett1 *= (real(1.0) - sncovr);
        etns1 *= (real(1.0) - sncovr);

        edir = edir1 * real(1000.0);
        ec = ec1 * real(1000.0);
        for (int k = 0; k < nsoil; ++k) et[k] = et1[k] * real(1000.0);
        ett = ett1 * real(1000.0);
        etns = etns1 * real(1000.0);
      }
      esnow = etp * sncovr;
      esnow1 = esnow * real(0.001);
      esnow2 = esnow1 * dt;
      etanrg = esnow * lsubs + etns * subc;
    }

    flx1 = real(0.0);
    if (snowng) {
      flx1 = cpice * prcp * (t1 - sfctmp);
    } else {
      if (prcp > real(0.0)) flx1 = cph2o * prcp * (t1 - sfctmp);
    }

    real dsoil = - (real(0.5) * zsoil[0]);
    real dtot = snowh + dsoil;
    real denom = real(1.0) + df1 / (dtot * rr * rch);
    real t12a = ((fdown - flx1 - flx2 - emissi * sigma * t24) / rch + th2 - sfctmp - etanrg / rch) / rr;
    real t12b = df1 * stc[0] / (dtot * rr * rch);
    real t12 = (sfctmp + t12a + t12b) / denom;
    real ex = real(0.0);

    if (t12 <= tfreez) {
      t1 = t12;
      ssoil = df1 * (t1 - stc[0]) / dtot;
      esd = std::max(real(0.0), esd - esnow2);
      flx3 = real(0.0);
      ex = real(0.0);
      snomlt = real(0.0);
    } else {
      real sn_factor = std::max(real(0.01), std::pow(sncovr, snoexp));
      t1 = tfreez * sn_factor + t12 * (real(1.0) - sn_factor);
      beta = real(1.0);
      ssoil = df1 * (t1 - stc[0]) / dtot;

      if (esd - esnow2 <= esdmin) {
        esd = real(0.0);
        ex = real(0.0);
        snomlt = real(0.0);
        flx3 = real(0.0);
      } else {
        esd -= esnow2;
        real seh = rch * (t1 - th2);
        real t14 = t1 * t1;  // T14 = T1*T1; T14 = T14*T14
        t14 = t14 * t14;
        flx3 = fdown - flx1 - flx2 - emissi * sigma * t14 - ssoil - seh - etanrg;
        if (flx3 <= real(0.0)) flx3 = real(0.0);
        ex = flx3 * real(0.001) / lsubf;
        snomlt = ex * dt;
        if (esd - snomlt >= esdmin) {
          esd -= snomlt;
        } else {
          ex = esd / dt;
          flx3 = ex * real(1000.0) * lsubf;
          snomlt = esd;
          esd = real(0.0);
        }
      }
      prcp1 += ex;
    }

    smflx(smc, nsoil, cmc, dt, prcp1, zsoil, sh2o, slope, kdt, frzfact, smcmax, bexp, smcwlt,
          dksat, dwsat, shdfac, cmcmax, runoff1, runoff2, runoff3, edir1, ec1, et1, drip);

    real zz1 = real(1.0);
    real yy = stc[0] - real(0.5) * ssoil * zsoil[0] * zz1 / df1;
    real ssoil1;
    real t11 = t1;
    shflx(ssoil1, stc, smc, smcmax, nsoil, t11, dt, yy, zz1, zsoil, tbot, zbot, smcwlt,
          psisat, sh2o, bexp, f1, df1, quartz, csoil, vegtyp, isurban, soiltyp, opt_thcnd);

    if (esd > real(0.0)) {
      snowpack(esd, dt, snowh, sndens, t1, yy, snomlt);
    } else {
      esd = real(0.0);
      snowh = real(0.0);
      sndens = real(0.0);
      sncovr = real(0.0);
    }
  }

  CPPWRF_HD void operator()(int i, int j) const {
    int vtyp = vegtyp(i, j);
    int styp = soiltyp(i, j);
    int slop = slopetyp(i, j);

    if (active.p) {
      if (active(i, j) == 0) return;
    } else if (vtyp == 16 || styp == 14) {
      // standalone test: skip water points (USGS category 16, STATSGO soil 14)
      return;
    }

    // Initialize outputs
    runoff1(i, j) = real(0.0);
    runoff2(i, j) = real(0.0);
    runoff3(i, j) = real(0.0);
    snomlt(i, j) = real(0.0);

    // Compute layer depth ZSOIL (negative depth below ground)
    real zsoil[kMaxSoilLevels];
    zsoil[0] = - sldpth[0];
    for (int kz = 1; kz < nsoil; ++kz) {
      zsoil[kz] = - sldpth[kz] + zsoil[kz - 1];
    }

    // Retrieve parameters from tables via REDPRM
    real csoil = tables.csoil_data;
    real bexp = tables.bb[styp - 1];
    real dksat = tables.satdk[styp - 1];
    real dwsat = tables.satdw[styp - 1];
    real f1 = tables.f11[styp - 1];
    real psisat = tables.satpsi[styp - 1];
    real quartz = tables.qtz[styp - 1];
    real smcdry = tables.drysmc[styp - 1];
    real smcmax = tables.maxsmc[styp - 1];
    real smcref = tables.refsmc[styp - 1];
    real smcwlt = tables.wltsmc[styp - 1];

    real zbot = tables.zbot_data;
    real salp = tables.salp_data;
    real sbeta = tables.sbeta_data;
    real refdk = tables.refdk_data;
    real frzk = tables.frzk_data;
    real fxexp = tables.fxexp_data;
    real refkdt = tables.refkdt_data;
    real kdt = refkdt * dksat / refdk;
    real czil = tables.czil_data;
    real slope = tables.slope_data[slop - 1];
    real lvcoef = tables.lvcoef_data;

    real frzfact = (smcmax / smcref) * (real(0.412) / real(0.468));
    real frzx = frzk * frzfact;

    real topt = tables.topt_data;
    real cmcmax = tables.cmcmax_data;
    real cfactr = tables.cfactr_data;
    real rsmax = tables.rsmax_data;
    int nroot = tables.nrotbl[vtyp - 1];
    real snup = tables.snuptbl[vtyp - 1];
    real rsmin_val = tables.rstbl[vtyp - 1];
    real rgl = tables.rgltbl[vtyp - 1];
    real hs = tables.hstbl[vtyp - 1];
    real emissmin = tables.emissmintbl[vtyp - 1];
    real emissmax = tables.emissmaxtbl[vtyp - 1];
    real laimin = tables.laimintbl[vtyp - 1];
    real laimax = tables.laimaxtbl[vtyp - 1];
    real z0min = tables.z0mintbl[vtyp - 1];
    real z0max = tables.z0maxtbl[vtyp - 1];
    real albedomin = tables.albedomintbl[vtyp - 1];
    real albedomax = tables.albedomaxtbl[vtyp - 1];

    real cur_shdfac = shdfac(i, j);
    if (vtyp == tables.bare) cur_shdfac = real(0.0);

    real rtdis[kMaxSoilLevels];
    for (int k = 0; k < nroot; ++k) {
      rtdis[k] = - sldpth[k] / zsoil[nroot - 1];
    }

    // Urban vegetation override
    if (vtyp == isurban) {
      cur_shdfac = real(0.05);
      rsmin_val = real(400.0);
      smcmax = real(0.45);
      smcref = real(0.42);
      smcwlt = real(0.40);
      smcdry = real(0.40);
    }
    rsmin(i, j) = rsmin_val;

    real cur_embrd;
    real cur_alb = alb(i, j);
    real cur_z0brd;
    real cur_xlai = xlai(i, j);
    real cur_shdmax = shdmax(i, j);
    real cur_shdmin = shdmin(i, j);

    if (cur_shdfac >= cur_shdmax) {
      cur_embrd = emissmax;
      if (!rdlai2d) cur_xlai = laimax;
      if (!usemonalb) cur_alb = albedomin;
      cur_z0brd = z0max;
    } else if (cur_shdfac <= cur_shdmin) {
      cur_embrd = emissmin;
      if (!rdlai2d) cur_xlai = laimin;
      if (!usemonalb) cur_alb = albedomax;
      cur_z0brd = z0min;
    } else {
      if (cur_shdmax > cur_shdmin) {
        real interp_fraction = (cur_shdfac - cur_shdmin) / (cur_shdmax - cur_shdmin);
        interp_fraction = std::min(interp_fraction, real(1.0));
        interp_fraction = std::max(interp_fraction, real(0.0));
        cur_embrd = ((real(1.0) - interp_fraction) * emissmin) + (interp_fraction * emissmax);
        if (!rdlai2d) cur_xlai = ((real(1.0) - interp_fraction) * laimin) + (interp_fraction * laimax);
        if (!usemonalb) cur_alb = ((real(1.0) - interp_fraction) * albedomax) + (interp_fraction * albedomin);
        cur_z0brd = ((real(1.0) - interp_fraction) * z0min) + (interp_fraction * z0max);
      } else {
        cur_embrd = real(0.5) * emissmin + real(0.5) * emissmax;
        if (!rdlai2d) cur_xlai = real(0.5) * laimin + real(0.5) * laimax;
        if (!usemonalb) cur_alb = real(0.5) * albedomin + real(0.5) * albedomax;
        cur_z0brd = real(0.5) * z0min + real(0.5) * z0max;
      }
    }
    embrd(i, j) = cur_embrd;
    xlai(i, j) = cur_xlai;
    z0brd(i, j) = cur_z0brd;
    if (alb_out.p) alb_out(i, j) = cur_alb;

    // Precipitation checks and snowpack update
    bool snowng = false;
    bool frzgra = false;
    constexpr real tfreez = real(273.15);

    real cur_sneqv = sneqv(i, j);
    real cur_snowh = snowh(i, j);
    real sndens = real(0.0);
    real sncond = real(1.0);

    if (cur_sneqv <= real(1.0e-7)) {
      cur_sneqv = real(0.0);
      sndens = real(0.0);
      cur_snowh = real(0.0);
      sncond = real(1.0);
    } else {
      sndens = cur_sneqv / cur_snowh;
      csnow(sncond, sndens);
    }

    real cur_prcp = prcp(i, j);
    real cur_sfctmp = sfctmp(i, j);
    real cur_t1 = t1(i, j);

    if (cur_prcp > real(0.0)) {
      if (ffrozp(i, j) > real(0.5)) {
        snowng = true;
      } else {
        if (cur_t1 <= tfreez) frzgra = true;
      }
    }

    real prcpf;
    if (snowng || frzgra) {
      real sn_new = cur_prcp * dt * real(0.001);
      cur_sneqv += sn_new;
      prcpf = real(0.0);
      snow_new(cur_sfctmp, sn_new, cur_snowh, sndens);
      csnow(sncond, sndens);
    } else {
      prcpf = cur_prcp;
    }

    // Snow cover and albedo
    real cur_sncovr;
    real cur_albedo;
    real cur_emissi;
    real cur_snotime1 = snotime1(i, j);

    if (cur_sneqv == real(0.0)) {
      cur_sncovr = real(0.0);
      cur_albedo = cur_alb;
      cur_emissi = cur_embrd;
    } else {
      snfrac(cur_sneqv, snup, salp, cur_snowh, cur_sncovr);
      cur_sncovr = std::min(cur_sncovr, real(0.98));
      alcalc(cur_alb, snoalb(i, j), cur_embrd, cur_shdfac, cur_shdmin, cur_sncovr,
             cur_t1, cur_albedo, cur_emissi, dt, snowng, cur_snotime1, lvcoef);
    }
    snotime1(i, j) = cur_snotime1;
    sncovr(i, j) = cur_sncovr;
    albedo(i, j) = cur_albedo;
    emissi(i, j) = cur_emissi;

    // Soil thermal conductivity & subsurface heat flux
    real stc_col[kMaxSoilLevels];
    real smc_col[kMaxSoilLevels];
    real sh2o_col[kMaxSoilLevels];
    for (int k = 0; k < nsoil; ++k) {
      stc_col[k] = stc(i, j, k);
      smc_col[k] = smc(i, j, k);
      sh2o_col[k] = sh2o(i, j, k);
    }

    int opt_thcnd = 1;
    real df1;
    tdfcnd(df1, smc_col[0], quartz, smcmax, sh2o_col[0], bexp, psisat, styp, opt_thcnd);
    if (vtyp == isurban) df1 = real(3.24);
    df1 = df1 * std::exp(sbeta * cur_shdfac);
    if (cur_sncovr > real(0.97)) df1 = sncond;

    real dsoil = - (real(0.5) * zsoil[0]);
    real ssoil_val;
    if (cur_sneqv == real(0.0)) {
      ssoil_val = df1 * (cur_t1 - stc_col[0]) / dsoil;
    } else {
      real dtot = cur_snowh + dsoil;
      real frcsno = cur_snowh / dtot;
      real frcsoi = dsoil / dtot;
      real df1a = frcsno * sncond + frcsoi * df1;
      df1 = df1a * cur_sncovr + df1 * (real(1.0) - cur_sncovr);
      ssoil_val = df1 * (cur_t1 - stc_col[0]) / dtot;
    }

    real cur_z0;
    if (cur_sncovr > real(0.0)) {
      snowz0(cur_sncovr, cur_z0, cur_z0brd, cur_snowh);
    } else {
      cur_z0 = cur_z0brd;
    }
    z0(i, j) = cur_z0;

    // Radiative and atmospheric forcing
    real cur_fdown = solnet(i, j) + lwdn(i, j);
    fdown(i, j) = cur_fdown;
    real cur_q2 = q2(i, j);
    real t2v = cur_sfctmp * (real(1.0) + real(0.61) * cur_q2);

    real t24, etp_val, rch_val, epsca_val, rr_val, flx2_val;
    real aoasis = real(1.0);
    penman(cur_sfctmp, sfcprs(i, j), ch(i, j), t2v, th2(i, j), cur_prcp, cur_fdown, t24,
           ssoil_val, cur_q2, q2sat(i, j), etp_val, rch_val, epsca_val, rr_val, snowng,
           frzgra, dqsdt2(i, j), flx2_val, cur_emissi, cur_sneqv, cur_t1, cur_sncovr,
           aoasis, cur_albedo, soldn(i, j), stc_col[0]);

    // Canopy resistance
    real rc_val = real(0.0);
    real pc_val = real(0.0);
    real rcs_val, rct_val, rcq_val, rcsoil_val;
    if (cur_shdfac > real(0.0) && cur_xlai > real(0.0)) {
      canres(soldn(i, j), ch(i, j), cur_sfctmp, cur_q2, sfcprs(i, j), sh2o_col, zsoil, nsoil,
             smcwlt, smcref, rsmin_val, rc_val, pc_val, nroot, q2sat(i, j), dqsdt2(i, j),
             topt, rsmax, rgl, hs, cur_xlai, rcs_val, rct_val, rcq_val, rcsoil_val, cur_emissi);
    }
    rc(i, j) = rc_val;
    pc(i, j) = pc_val;

    // Energy balance & hydrological solvers (NOPAC vs SNOPAC)
    real eta_val = real(0.0);
    real ec_val = real(0.0);
    real edir_val = real(0.0);
    real et_col[kMaxSoilLevels] = {0};
    real ett_val = real(0.0);
    real esnow_val = real(0.0);
    real drip_val = real(0.0);
    real dew_val = real(0.0);
    real beta_val = real(0.0);
    real flx1_val = real(0.0);
    real flx3_val = real(0.0);
    real snomlt_val = real(0.0);
    real r1_val = real(0.0);
    real r2_val = real(0.0);
    real r3_val = real(0.0);
    real eta_kin = real(0.0);

    real cur_cmc = cmc(i, j);
    real cur_tbot = tbot(i, j);
    real cur_th2 = th2(i, j);
    real cur_ribb = ribb(i, j);

    if (cur_sneqv == real(0.0)) {
      nopac(etp_val, eta_val, cur_prcp, smc_col, smcmax, smcwlt, smcref, smcdry, cur_cmc, cmcmax,
            nsoil, dt, cur_shdfac, sbeta, cur_q2, cur_t1, cur_sfctmp, t24, cur_th2, cur_fdown,
            f1, cur_emissi, ssoil_val, stc_col, epsca_val, bexp, pc_val, rch_val, rr_val,
            cfactr, sh2o_col, slope, kdt, frzx, psisat, zsoil, dksat, dwsat, cur_tbot, zbot,
            r1_val, r2_val, r3_val, edir_val, ec_val, et_col, ett_val, nroot, rtdis, quartz,
            fxexp, csoil, beta_val, drip_val, dew_val, flx1_val, flx3_val, vtyp, isurban, styp, opt_thcnd);
      eta_kin = eta_val;
    } else {
      real etns_val = real(0.0);
      snopac(etp_val, eta_val, cur_prcp, prcpf, snowng, smc_col, smcmax, smcwlt, smcref, smcdry,
             cur_cmc, cmcmax, nsoil, dt, sbeta, df1, cur_q2, cur_t1, cur_sfctmp, t24, cur_th2,
             cur_fdown, f1, ssoil_val, stc_col, epsca_val, sfcprs(i, j), bexp, pc_val, rch_val,
             rr_val, cfactr, cur_sncovr, cur_sneqv, sndens, cur_snowh, sh2o_col, slope, kdt,
             frzx, psisat, zsoil, dwsat, dksat, cur_tbot, zbot, cur_shdfac, r1_val, r2_val,
             r3_val, edir_val, ec_val, et_col, ett_val, nroot, snomlt_val, rtdis, quartz,
             fxexp, csoil, beta_val, drip_val, dew_val, flx1_val, flx2_val, flx3_val, esnow_val,
             etns_val, cur_emissi, cur_ribb, soldn(i, j), isurban, vtyp, styp, opt_thcnd);
      eta_kin = esnow_val + etns_val - real(1000.0) * dew_val;
    }

    // Diagnostics and unit conversions
    q1(i, j) = cur_q2 + eta_kin * real(cst::cp) / rch_val;
    constexpr real r_air = real(287.04);
    real sheat_val = - (ch(i, j) * real(cst::cp) * sfcprs(i, j)) / (r_air * t2v) * (cur_th2 - cur_t1);
    sheat(i, j) = sheat_val;

    constexpr real lvh2o = real(2.501e6);
    constexpr real lsubs = real(2.83e6);

    edir_val *= lvh2o;
    ec_val *= lvh2o;
    for (int k = 0; k < nsoil; ++k) {
      et_col[k] *= lvh2o;
      et(i, j, k) = et_col[k];
    }
    ett_val *= lvh2o;
    esnow_val *= lsubs;
    etp_val = etp_val * ((real(1.0) - cur_sncovr) * lvh2o + cur_sncovr * lsubs);
    if (etp_val > real(0.0)) {
      eta_val = edir_val + ec_val + ett_val + esnow_val;
    } else {
      eta_val = etp_val;
    }
    if (etp_val == real(0.0)) {
      beta_val = real(0.0);
    } else {
      beta_val = eta_val / etp_val;
    }
    ssoil_val = - real(1.0) * ssoil_val;

    r3_val /= dt;
    r2_val += r3_val;

    real cur_soilm = - real(1.0) * smc_col[0] * zsoil[0];
    for (int k = 1; k < nsoil; ++k) {
      cur_soilm += smc_col[k] * (zsoil[k - 1] - zsoil[k]);
    }
    real soilwm = - real(1.0) * (smcmax - smcwlt) * zsoil[0];
    real soilww = - real(1.0) * (smc_col[0] - smcwlt) * zsoil[0];

    for (int k = 0; k < nsoil; ++k) {
      smav(i, j, k) = (smc_col[k] - smcwlt) / (smcmax - smcwlt);
    }

    if (nroot >= 2) {
      for (int k = 1; k < nroot; ++k) {
        soilwm += (smcmax - smcwlt) * (zsoil[k - 1] - zsoil[k]);
        soilww += (smc_col[k] - smcwlt) * (zsoil[k - 1] - zsoil[k]);
      }
    }
    real cur_soilw;
    if (soilwm < real(1.0e-6)) {
      soilwm = real(0.0);
      cur_soilw = real(0.0);
      cur_soilm = real(0.0);
    } else {
      cur_soilw = soilww / soilwm;
    }

    // Write back updated states and outputs
    cmc(i, j) = cur_cmc;
    t1(i, j) = cur_t1;
    snowh(i, j) = cur_snowh;
    sneqv(i, j) = cur_sneqv;
    ribb(i, j) = cur_ribb;

    for (int k = 0; k < nsoil; ++k) {
      stc(i, j, k) = stc_col[k];
      smc(i, j, k) = smc_col[k];
      sh2o(i, j, k) = sh2o_col[k];
    }

    eta(i, j) = eta_val;
    eta_kinematic(i, j) = eta_kin;
    ec(i, j) = ec_val;
    edir(i, j) = edir_val;
    ett(i, j) = ett_val;
    esnow(i, j) = esnow_val;
    drip(i, j) = drip_val;
    dew(i, j) = dew_val;
    beta(i, j) = beta_val;
    etp(i, j) = etp_val;
    ssoil(i, j) = ssoil_val;
    flx1(i, j) = flx1_val;
    flx2(i, j) = flx2_val;
    flx3(i, j) = flx3_val;
    snomlt(i, j) = snomlt_val;
    runoff1(i, j) = r1_val;
    runoff2(i, j) = r2_val;
    runoff3(i, j) = r3_val;
    soilw(i, j) = cur_soilw;
    soilm(i, j) = cur_soilm;
  }
};

void run_noah(NoahKernel kernel, const Range2D& range);

}  // namespace cppwrf::phys::noah
