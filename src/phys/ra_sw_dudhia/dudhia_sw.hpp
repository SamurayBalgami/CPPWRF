// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// WRF phys/module_ra_sw.F (Dudhia shortwave, ra_sw_physics = 1): SWRAD's column set-up and SWPARA for one
// column, with WRF's module constants (r_d = 287, cp = 7*r_d/2) and without the WRF-Chem aerosols (pm2_5 = 0)
// and the eclipse (obscur = 0). SWRAD passes no slope_rad arguments, so SWPARA's slope code never runs.
// WRF's 1-D arrays run top-down (K = 1 is the model top, NK = kme-1-K+kms); every one of them but SDOWN
// depends on its own level only, so they are formed in the layer loop.
#include <algorithm>
#include <cmath>

#include "phys/ra_rrtmg_common/strided.hpp"

namespace cppwrf::phys {

// swinit: CSSCA = swrad_scat * 1.e-5
CPPWRF_HD inline double dudhia_cssca(double swrad_scat) { return swrad_scat * 1.e-5; }

// Model levels k = 0 (lowest) .. nz-1: t (K), p (Pa), dz (dz8w, m); moist species qv..qg, a null pointer for
// a species the scheme lacks (WRF F_Q*). warm_rain (Kessler): without qi the cloud and rain below 273.15 K are
// not turned into ice and snow. tten: the heating TTEN1D(NK) at level k (K/s, 0 at night); returns GSW.
CPPWRF_HD inline double dudhia_swpara(int nz, CSD t, CSD p, CSD dz, CSD qv, CSD qc, CSD qr, CSD qi, CSD qs, CSD qg,
                                      double coszen, double albedo, double solcon, int icloud, double cssca,
                                      bool warm_rain, SD tten) {
  constexpr double R = 287., CP = 7. * 287. / 2.;
  constexpr double ALBTAB[5][4] = {{0., 0., 0., 0.}, {69., 58., 40., 15.}, {90., 80., 70., 60.},
                                   {94., 90., 82., 78.}, {96., 92., 85., 80.}};  // ALBTAB(i,j) = ALBTAB[j-1][i-1]
  constexpr double ABSTAB[5][4] = {{0., 0., 0., 0.}, {0., 2.5, 4., 5.}, {0., 2.6, 7., 10.},
                                   {0., 3.3, 10., 14.}, {0., 3.7, 10., 15.}};
  constexpr double XMUVAL[4] = {0., 0.2, 0.5, 1.0};
  for (int k = 0; k < nz; ++k) tten[k] = 0.;
  double GSW = 0.0;
  const double bext340 = 5.E-6, bexth2o = 5.E-6;
  const double SOLTOP = solcon * (1 - 0.);
  const double CSZA = coszen;
  if (CSZA <= 1.E-9) return GSW;  // night
  if (SOLTOP <= 1.E-9) return GSW;
  const double XMU = CSZA;
  const double SDOWN1 = SOLTOP * XMU;
  double SDOWN = SDOWN1;  // SDOWN(K)
  double WW = 0., UV = 0., OLDALB = 0., OLDABC = 0., TOTABS = 0., DSCA = 0., DABS = 0., DSCLD = 0., DABSA = 0.;
  for (int K = 1; K <= nz; ++K) {
    const int m = nz - K;  // NK
    const double T = t[m], P = p[m], DZ = dz[m];
    // SWRAD: max(0, q) of the species present; ice and snow from cloud and rain without F_QI
    double QV = qv.p ? std::max(0., qv[m]) : 0., QC = qc.p ? std::max(0., qc[m]) : 0.;
    double QR = qr.p ? std::max(0., qr[m]) : 0., QI = 0., QS = 0., QG = 0.;
    if (qi.p) {
      QI = std::max(0., qi[m]);
    } else if (!warm_rain && T < 273.15) {
      QI = QC;
      QC = 0.;
      QS = QR;
      QR = 0.;
    }
    if (qs.p) QS = std::max(0., qs[m]);
    if (qg.p) QG = std::max(0., qg[m]);
    const double RO = P / (R * T);
    const double XWVP = RO * QV * DZ * 1000.;
    const double XATP = RO * DZ;
    const double XLWP = icloud == 0 ? 0. : RO * 1000. * DZ * (QC + 0.1 * QI + 0.05 * QR + 0.02 * QS + 0.05 * QG);
    WW = WW + XLWP;
    UV = UV + XWVP;
    const double WGM = WW / XMU;
    const double UGCM = UV * 0.0001 / XMU;
    const double OLDABS = TOTABS;
    TOTABS = 2.9 * UGCM / (libm::pow(1. + 141.5 * UGCM, 0.635) + 5.925 * UGCM);
    const double beta = 0.4 * (1.0 - XMU) + 0.1;
    const double aer_dry1 = 0., aer_water1 = 0.;
    double XSCA = (cssca * XATP + beta * aer_dry1 * bext340 * DZ + beta * aer_water1 * bexth2o * DZ) / XMU;
    double XABS = (TOTABS - OLDABS) * (SDOWN1 - DSCLD - DSCA - DABSA) / SDOWN;
    double XABSA = 0.;
    if (XABS < 0.) XABS = 0.;
    double ALW = libm::log10(WGM + 1.);
    if (ALW > 3.999) ALW = 3.999;
    int IIL = 0, IU = 0;
    double XI = 0.;
    for (int II = 1; II <= 3; ++II) {
      if (XMU > XMUVAL[II - 1]) {
        IIL = II;
        IU = II + 1;
        XI = (XMU - XMUVAL[II - 1]) / (XMUVAL[II] - XMUVAL[II - 1]) + double(IIL);
      }
    }
    const int JJL = int(ALW) + 1;
    const int JU = JJL + 1;
    const double YJ = ALW + 1.;
    const double dI1 = XI - double(IIL), dI2 = double(IU) - XI, dJ1 = YJ - double(JJL), dJ2 = double(JU) - YJ;
    const double den = double((IU - IIL) * (JU - JJL));
    const double ALBA = (ALBTAB[JU - 1][IU - 1] * dI1 * dJ1 + ALBTAB[JU - 1][IIL - 1] * dI2 * dJ1 +
                         ALBTAB[JJL - 1][IU - 1] * dI1 * dJ2 + ALBTAB[JJL - 1][IIL - 1] * dI2 * dJ2) / den;
    const double ABSC = (ABSTAB[JU - 1][IU - 1] * dI1 * dJ1 + ABSTAB[JU - 1][IIL - 1] * dI2 * dJ1 +
                         ABSTAB[JJL - 1][IU - 1] * dI1 * dJ2 + ABSTAB[JJL - 1][IIL - 1] * dI2 * dJ2) / den;
    double XALB = (ALBA - OLDALB) * (SDOWN1 - DSCA - DABS) / SDOWN;
    double XABSC = (ABSC - OLDABC) * (SDOWN1 - DSCA - DABS) / SDOWN;
    if (XALB < 0.) XALB = 0.;
    if (XABSC < 0.) XABSC = 0.;
    DSCLD = DSCLD + (XALB + XABSC) * SDOWN * 0.01;
    DSCA = DSCA + XSCA * SDOWN;
    DABS = DABS + XABS * SDOWN;
    DABSA = DABSA + XABSA * SDOWN;
    OLDALB = ALBA;
    OLDABC = ABSC;
    double TRANS0 = 100. - XALB - XABSC - XABS * 100. - XSCA * 100.;
    if (TRANS0 < 1.) {
      const double FF = 99. / (XALB + XABSC + XABS * 100. + XSCA * 100.);
      XALB = XALB * FF;
      XABSC = XABSC * FF;
      XABS = XABS * FF;
      XSCA = XSCA * FF;
      TRANS0 = 1.;
    }
    const double SDOWN_NEXT = std::max(1.E-9, SDOWN * TRANS0 * 0.01);
    tten[m] = SDOWN * (XABSC + XABS * 100. + XABSA * 100.) * 0.01 / (RO * CP * DZ);
    SDOWN = SDOWN_NEXT;
  }
  GSW = (1. - albedo) * SDOWN;
  return GSW;
}

// radiation_driver's surface direct / diffuse split for SW schemes without their own (Ruiz-Arias et al. 2010),
// where coszen > 1e-3: swddif and swddir from SWDOWN, the solar constant and the terrain height ht (m); the
// caller forms swddni = swddir / max(coszen, 1e-4) from the stored swddir. Returns false where coszen <= 1e-3
CPPWRF_HD inline bool sw_diffuse_backup(double coszen, double swdown, double solcon, double ht, double& swddif,
                                        double& swddir) {
  if (!(coszen > 1e-3)) return false;
  const double ioh = solcon * coszen;
  double kt = swdown / std::max(ioh, 1e-3);
  const double airmass =
      libm::exp(-ht / 8434.5) / (coszen + 0.50572 * libm::pow(libm::asin(coszen) * 57.295779513082323 + 6.07995, -1.6364));
  kt = kt / (0.1 + 1.031 * libm::exp(-1.4 / (0.9 + (9.4 / std::max(airmass, 1e-3)))));
  const double kd = 0.952 - 1.041 * libm::exp(-libm::exp(2.300 - 4.702 * kt));
  swddif = kd * swdown;
  swddir = (1. - kd) * swdown;
  return true;
}

}  // namespace cppwrf::phys
