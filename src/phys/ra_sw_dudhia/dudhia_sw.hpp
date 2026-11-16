// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "phys/ra_rrtmg_common/rrtmg_real.hpp"
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
CPPWRF_HD inline rr dudhia_cssca(rr swrad_scat) { return swrad_scat * rr(1.e-5); }

// Model levels k = 0 (lowest) .. nz-1: t (K), p (Pa), dz (dz8w, m); moist species qv..qg, a null pointer for
// a species the scheme lacks (WRF F_Q*). warm_rain (Kessler): without qi the cloud and rain below 273.15 K are
// not turned into ice and snow. tten: the heating TTEN1D(NK) at level k (K/s, 0 at night); returns GSW.
CPPWRF_HD inline rr dudhia_swpara(int nz, CSD t, CSD p, CSD dz, CSD qv, CSD qc, CSD qr, CSD qi, CSD qs, CSD qg,
                                      rr coszen, rr albedo, rr solcon, int icloud, rr cssca,
                                      bool warm_rain, SD tten) {
  constexpr rr R = rr(287.), CP = rr(7.) * rr(287.) / rr(2.);
  constexpr rr ALBTAB[5][4] = {{rr(0.), rr(0.), rr(0.), rr(0.)}, {rr(69.), rr(58.), rr(40.), rr(15.)}, {rr(90.), rr(80.), rr(70.), rr(60.)},
                                   {rr(94.), rr(90.), rr(82.), rr(78.)}, {rr(96.), rr(92.), rr(85.), rr(80.)}};  // ALBTAB(i,j) = ALBTAB[j-1][i-1]
  constexpr rr ABSTAB[5][4] = {{rr(0.), rr(0.), rr(0.), rr(0.)}, {rr(0.), rr(2.5), rr(4.), rr(5.)}, {rr(0.), rr(2.6), rr(7.), rr(10.)},
                                   {rr(0.), rr(3.3), rr(10.), rr(14.)}, {rr(0.), rr(3.7), rr(10.), rr(15.)}};
  constexpr rr XMUVAL[4] = {rr(0.), rr(0.2), rr(0.5), rr(1.0)};
  for (int k = 0; k < nz; ++k) tten[k] = rr(0.);
  rr GSW = rr(0.0);
  const rr bext340 = rr(5.E-6), bexth2o = rr(5.E-6);
  const rr SOLTOP = solcon * (1 - rr(0.));
  const rr CSZA = coszen;
  if (CSZA <= rr(1.E-9)) return GSW;  // night
  if (SOLTOP <= rr(1.E-9)) return GSW;
  const rr XMU = CSZA;
  const rr SDOWN1 = SOLTOP * XMU;
  rr SDOWN = SDOWN1;  // SDOWN(K)
  rr WW = rr(0.), UV = rr(0.), OLDALB = rr(0.), OLDABC = rr(0.), TOTABS = rr(0.), DSCA = rr(0.), DABS = rr(0.), DSCLD = rr(0.), DABSA = rr(0.);
  for (int K = 1; K <= nz; ++K) {
    const int m = nz - K;  // NK
    const rr T = t[m], P = p[m], DZ = dz[m];
    // SWRAD: max(0, q) of the species present; ice and snow from cloud and rain without F_QI
    rr QV = qv.p ? std::max(rr(0.), qv[m]) : rr(0.), QC = qc.p ? std::max(rr(0.), qc[m]) : rr(0.);
    rr QR = qr.p ? std::max(rr(0.), qr[m]) : rr(0.), QI = rr(0.), QS = rr(0.), QG = rr(0.);
    if (qi.p) {
      QI = std::max(rr(0.), qi[m]);
    } else if (!warm_rain && T < rr(273.15)) {
      QI = QC;
      QC = rr(0.);
      QS = QR;
      QR = rr(0.);
    }
    if (qs.p) QS = std::max(rr(0.), qs[m]);
    if (qg.p) QG = std::max(rr(0.), qg[m]);
    const rr RO = P / (R * T);
    const rr XWVP = RO * QV * DZ * rr(1000.);
    const rr XATP = RO * DZ;
    const rr XLWP = icloud == 0 ? rr(0.) : RO * rr(1000.) * DZ * (QC + rr(0.1) * QI + rr(0.05) * QR + rr(0.02) * QS + rr(0.05) * QG);
    WW = WW + XLWP;
    UV = UV + XWVP;
    const rr WGM = WW / XMU;
    const rr UGCM = UV * rr(0.0001) / XMU;
    const rr OLDABS = TOTABS;
    TOTABS = rr(2.9) * UGCM / (libm::pow(rr(1.) + rr(141.5) * UGCM, rr(0.635)) + rr(5.925) * UGCM);
    const rr beta = rr(0.4) * (rr(1.0) - XMU) + rr(0.1);
    const rr aer_dry1 = rr(0.), aer_water1 = rr(0.);
    rr XSCA = (cssca * XATP + beta * aer_dry1 * bext340 * DZ + beta * aer_water1 * bexth2o * DZ) / XMU;
    rr XABS = (TOTABS - OLDABS) * (SDOWN1 - DSCLD - DSCA - DABSA) / SDOWN;
    rr XABSA = rr(0.);
    if (XABS < rr(0.)) XABS = rr(0.);
    rr ALW = libm::log10(WGM + rr(1.));
    if (ALW > rr(3.999)) ALW = rr(3.999);
    int IIL = 0, IU = 0;
    rr XI = rr(0.);
    for (int II = 1; II <= 3; ++II) {
      if (XMU > XMUVAL[II - 1]) {
        IIL = II;
        IU = II + 1;
        XI = (XMU - XMUVAL[II - 1]) / (XMUVAL[II] - XMUVAL[II - 1]) + rr(IIL);
      }
    }
    const int JJL = int(ALW) + 1;
    const int JU = JJL + 1;
    const rr YJ = ALW + rr(1.);
    const rr dI1 = XI - rr(IIL), dI2 = rr(IU) - XI, dJ1 = YJ - rr(JJL), dJ2 = rr(JU) - YJ;
    const rr den = rr((IU - IIL) * (JU - JJL));
    const rr ALBA = (ALBTAB[JU - 1][IU - 1] * dI1 * dJ1 + ALBTAB[JU - 1][IIL - 1] * dI2 * dJ1 +
                         ALBTAB[JJL - 1][IU - 1] * dI1 * dJ2 + ALBTAB[JJL - 1][IIL - 1] * dI2 * dJ2) / den;
    const rr ABSC = (ABSTAB[JU - 1][IU - 1] * dI1 * dJ1 + ABSTAB[JU - 1][IIL - 1] * dI2 * dJ1 +
                         ABSTAB[JJL - 1][IU - 1] * dI1 * dJ2 + ABSTAB[JJL - 1][IIL - 1] * dI2 * dJ2) / den;
    rr XALB = (ALBA - OLDALB) * (SDOWN1 - DSCA - DABS) / SDOWN;
    rr XABSC = (ABSC - OLDABC) * (SDOWN1 - DSCA - DABS) / SDOWN;
    if (XALB < rr(0.)) XALB = rr(0.);
    if (XABSC < rr(0.)) XABSC = rr(0.);
    DSCLD = DSCLD + (XALB + XABSC) * SDOWN * rr(0.01);
    DSCA = DSCA + XSCA * SDOWN;
    DABS = DABS + XABS * SDOWN;
    DABSA = DABSA + XABSA * SDOWN;
    OLDALB = ALBA;
    OLDABC = ABSC;
    rr TRANS0 = rr(100.) - XALB - XABSC - XABS * rr(100.) - XSCA * rr(100.);
    if (TRANS0 < rr(1.)) {
      const rr FF = rr(99.) / (XALB + XABSC + XABS * rr(100.) + XSCA * rr(100.));
      XALB = XALB * FF;
      XABSC = XABSC * FF;
      XABS = XABS * FF;
      XSCA = XSCA * FF;
      TRANS0 = rr(1.);
    }
    const rr SDOWN_NEXT = std::max(rr(1.E-9), SDOWN * TRANS0 * rr(0.01));
    tten[m] = SDOWN * (XABSC + XABS * rr(100.) + XABSA * rr(100.)) * rr(0.01) / (RO * CP * DZ);
    SDOWN = SDOWN_NEXT;
  }
  GSW = (rr(1.) - albedo) * SDOWN;
  return GSW;
}

// radiation_driver's surface direct / diffuse split for SW schemes without their own (Ruiz-Arias et al. 2010),
// where coszen > 1e-3: swddif and swddir from SWDOWN, the solar constant and the terrain height ht (m); the
// caller forms swddni = swddir / max(coszen, 1e-4) from the stored swddir. Returns false where coszen <= 1e-3
CPPWRF_HD inline bool sw_diffuse_backup(rr coszen, rr swdown, rr solcon, rr ht, rr& swddif,
                                        rr& swddir) {
  if (!(coszen > rr(1e-3))) return false;
  const rr ioh = solcon * coszen;
  rr kt = swdown / std::max(ioh, rr(1e-3));
  const rr airmass =
      libm::exp(-ht / rr(8434.5)) / (coszen + rr(0.50572) * libm::pow(libm::asin(coszen) * rr(57.295779513082323) + rr(6.07995), -rr(1.6364)));
  kt = kt / (rr(0.1) + rr(1.031) * libm::exp(-rr(1.4) / (rr(0.9) + (rr(9.4) / std::max(airmass, rr(1e-3))))));
  const rr kd = rr(0.952) - rr(1.041) * libm::exp(-libm::exp(rr(2.300) - rr(4.702) * kt));
  swddif = kd * swdown;
  swddir = (rr(1.) - kd) * swdown;
  return true;
}

}  // namespace cppwrf::phys
