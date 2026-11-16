// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "phys/ra_rrtmg_common/rrtmg_real.hpp"

#include "core/backend.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_tables.hpp"
#include <algorithm>
#include <cmath>

namespace cppwrf::phys {

CPPWRF_FORCEINLINE void reftra_sw_elem(
    bool chk,
    rr zg,
    rr prmuz,
    rr zto1,
    rr zw,
    rr& pref,
    rr& prefd,
    rr& ptra,
    rr& ptrad,
    const RrtmgSwTables& CPPWRF_RESTRICT tables
) {
  constexpr rr tblint = rr(10000.0);
  constexpr rr od_lo = rr(0.06);
  constexpr rr pade = rr(0.278);
  constexpr rr bpade = rr(1.0) / pade;
  constexpr rr zwcrit = rr(0.9999995);
  constexpr rr eps = rr(1.0e-8);

  if (!chk) {
    pref = rr(0.0);
    ptra = rr(1.0);
    prefd = rr(0.0);
    ptrad = rr(1.0);
    return;
  }

  rr zg3 = rr(3.0) * zg;
  // kmodts == 2
  rr zgamma1 = (rr(8.0) - zw * (rr(5.0) + zg3)) * rr(0.25);
  rr zgamma2 = rr(3.0) * (zw * (rr(1.0) - zg)) * rr(0.25);
  rr zgamma3 = (rr(2.0) - zg3 * prmuz) * rr(0.25);
  rr zgamma4 = rr(1.0) - zgamma3;

  rr zwo = rr(0.0);
  rr denom = rr(1.0);
  if (zg != rr(1.0)) {
    rr ratio = zg / (rr(1.0) - zg);
    denom = rr(1.0) - (rr(1.0) - zw) * (ratio * ratio);  // (zg/(1-zg))**2
  }
  if (zw > rr(0.0) && denom != rr(0.0)) {
    zwo = zw / denom;
  }

  if (zwo >= zwcrit) {
    // Conservative scattering
    rr za = zgamma1 * prmuz;
    rr za1 = za - zgamma3;
    rr zgt = zgamma1 * zto1;

    rr ze1 = std::min(zto1 / prmuz, rr(500.0));
    rr ze2;
    if (ze1 <= od_lo) {
      ze2 = rr(1.0) - ze1 + rr(0.5) * ze1 * ze1;
    } else {
      rr tblind = ze1 / (bpade + ze1);
      int itind = static_cast<int>(tblint * tblind + rr(0.5));
      ze2 = tables.exp_tbl[itind];
    }

    pref = (zgt - za1 * (rr(1.0) - ze2)) / (rr(1.0) + zgt);
    ptra = rr(1.0) - pref;

    prefd = zgt / (rr(1.0) + zgt);
    ptrad = rr(1.0) - prefd;

    if (ze2 == rr(1.0)) {
      pref = rr(0.0);
      ptra = rr(1.0);
      prefd = rr(0.0);
      ptrad = rr(1.0);
    }
  } else {
    // Non-conservative scattering
    rr za1 = zgamma1 * zgamma4 + zgamma2 * zgamma3;
    rr za2 = zgamma1 * zgamma3 + zgamma2 * zgamma4;
    rr zrk = std::sqrt(zgamma1 * zgamma1 - zgamma2 * zgamma2);
    rr zrp = zrk * prmuz;
    rr zrp1 = rr(1.0) + zrp;
    rr zrm1 = rr(1.0) - zrp;
    rr zrk2 = rr(2.0) * zrk;
    rr zrpp = rr(1.0) - zrp * zrp;
    rr zrkg = zrk + zgamma1;
    rr zr1 = zrm1 * (za2 + zrk * zgamma3);
    rr zr2 = zrp1 * (za2 - zrk * zgamma3);
    rr zr3 = zrk2 * (zgamma3 - za2 * prmuz);
    rr zr4 = zrpp * zrkg;
    rr zr5 = zrpp * (zrk - zgamma1);
    rr zt1 = zrp1 * (za1 + zrk * zgamma4);
    rr zt2 = zrm1 * (za1 - zrk * zgamma4);
    rr zt3 = zrk2 * (zgamma4 + za1 * prmuz);
    rr zt4 = zr4;
    rr zt5 = zr5;

    rr zbeta = (zgamma1 - zrk) / zrkg;

    rr ze1 = std::min(zrk * zto1, rr(500.0));
    rr ze2 = std::min(zto1 / prmuz, rr(500.0));

    rr zem1, zep1, zem2, zep2;
    if (ze1 <= od_lo) {
      zem1 = rr(1.0) - ze1 + rr(0.5) * ze1 * ze1;
      zep1 = rr(1.0) / zem1;
    } else {
      rr tblind = ze1 / (bpade + ze1);
      int itind = static_cast<int>(tblint * tblind + rr(0.5));
      zem1 = tables.exp_tbl[itind];
      zep1 = rr(1.0) / zem1;
    }

    if (ze2 <= od_lo) {
      zem2 = rr(1.0) - ze2 + rr(0.5) * ze2 * ze2;
      zep2 = rr(1.0) / zem2;
    } else {
      rr tblind = ze2 / (bpade + ze2);
      int itind = static_cast<int>(tblint * tblind + rr(0.5));
      zem2 = tables.exp_tbl[itind];
      zep2 = rr(1.0) / zem2;
    }

    rr zdenr = zr4 * zep1 + zr5 * zem1;
    rr zdent = zt4 * zep1 + zt5 * zem1;
    if (zdenr >= -eps && zdenr <= eps) {
      pref = eps;
      ptra = zem2;
    } else {
      pref = zw * (zr1 * zep1 - zr2 * zem1 - zr3 * zem2) / zdenr;
      ptra = zem2 - zem2 * zw * (zt1 * zep1 - zt2 * zem1 - zt3 * zep2) / zdent;
    }

    rr zemm = zem1 * zem1;
    rr zdend = rr(1.0) / ((rr(1.0) - zbeta * zemm) * zrkg);
    prefd = zgamma2 * (rr(1.0) - zemm) * zdend;
    ptrad = zrk2 * zem1 * zdend;
  }
}

}  // namespace cppwrf::phys
