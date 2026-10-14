#pragma once

#include "core/backend.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_tables.hpp"
#include <algorithm>
#include <cmath>

namespace cppwrf::phys {

CPPWRF_HD inline void reftra_sw(
    int nlayers,
    const bool* lrtchk,       // (nlayers)
    const double* pgg,        // (nlayers)
    double prmuz,
    const double* ptau,       // (nlayers)
    const double* pw,         // (nlayers)
    double* pref,             // (nlayers+1)
    double* prefd,            // (nlayers+1)
    double* ptra,             // (nlayers+1)
    double* ptrad,            // (nlayers+1)
    const RrtmgSwTables& tables
) {
  constexpr double tblint = 10000.0;
  constexpr double od_lo = 0.06;
  constexpr double pade = 0.278;
  constexpr double bpade = 1.0 / pade;
  constexpr double zwcrit = 0.9999995;
  constexpr double eps = 1.0e-8;

  // kmodts = 2 (PIFM)
  for (int jk = 0; jk < nlayers; ++jk) {
    if (!lrtchk[jk]) {
      pref[jk] = 0.0;
      ptra[jk] = 1.0;
      prefd[jk] = 0.0;
      ptrad[jk] = 1.0;
    } else {
      double zto1 = ptau[jk];
      double zw = pw[jk];
      double zg = pgg[jk];

      double zg3 = 3.0 * zg;
      // kmodts == 2
      double zgamma1 = (8.0 - zw * (5.0 + zg3)) * 0.25;
      double zgamma2 = 3.0 * (zw * (1.0 - zg)) * 0.25;
      double zgamma3 = (2.0 - zg3 * prmuz) * 0.25;
      double zgamma4 = 1.0 - zgamma3;

      double zwo = 0.0;
      double denom = 1.0;
      if (zg != 1.0) {
        double ratio = zg / (1.0 - zg);
        denom = 1.0 - (1.0 - zw) * ratio * ratio;
      }
      if (zw > 0.0 && denom != 0.0) {
        zwo = zw / denom;
      }

      if (zwo >= zwcrit) {
        // Conservative scattering
        double za = zgamma1 * prmuz;
        double za1 = za - zgamma3;
        double zgt = zgamma1 * zto1;

        double ze1 = std::min(zto1 / prmuz, 500.0);
        double ze2;
        if (ze1 <= od_lo) {
          ze2 = 1.0 - ze1 + 0.5 * ze1 * ze1;
        } else {
          double tblind = ze1 / (bpade + ze1);
          int itind = static_cast<int>(tblint * tblind + 0.5);
          ze2 = tables.exp_tbl[itind];
        }

        pref[jk] = (zgt - za1 * (1.0 - ze2)) / (1.0 + zgt);
        ptra[jk] = 1.0 - pref[jk];

        prefd[jk] = zgt / (1.0 + zgt);
        ptrad[jk] = 1.0 - prefd[jk];

        if (ze2 == 1.0) {
          pref[jk] = 0.0;
          ptra[jk] = 1.0;
          prefd[jk] = 0.0;
          ptrad[jk] = 1.0;
        }
      } else {
        // Non-conservative scattering
        double za1 = zgamma1 * zgamma4 + zgamma2 * zgamma3;
        double za2 = zgamma1 * zgamma3 + zgamma2 * zgamma4;
        double zrk = std::sqrt(zgamma1 * zgamma1 - zgamma2 * zgamma2);
        double zrp = zrk * prmuz;
        double zrp1 = 1.0 + zrp;
        double zrm1 = 1.0 - zrp;
        double zrk2 = 2.0 * zrk;
        double zrpp = 1.0 - zrp * zrp;
        double zrkg = zrk + zgamma1;
        double zr1 = zrm1 * (za2 + zrk * zgamma3);
        double zr2 = zrp1 * (za2 - zrk * zgamma3);
        double zr3 = zrk2 * (zgamma3 - za2 * prmuz);
        double zr4 = zrpp * zrkg;
        double zr5 = zrpp * (zrk - zgamma1);
        double zt1 = zrp1 * (za1 + zrk * zgamma4);
        double zt2 = zrm1 * (za1 - zrk * zgamma4);
        double zt3 = zrk2 * (zgamma4 + za1 * prmuz);
        double zt4 = zr4;
        double zt5 = zr5;

        double zbeta = (zgamma1 - zrk) / zrkg;

        double ze1 = std::min(zrk * zto1, 500.0);
        double ze2 = std::min(zto1 / prmuz, 500.0);

        double zem1, zep1, zem2, zep2;
        if (ze1 <= od_lo) {
          zem1 = 1.0 - ze1 + 0.5 * ze1 * ze1;
          zep1 = 1.0 / zem1;
        } else {
          double tblind = ze1 / (bpade + ze1);
          int itind = static_cast<int>(tblint * tblind + 0.5);
          zem1 = tables.exp_tbl[itind];
          zep1 = 1.0 / zem1;
        }

        if (ze2 <= od_lo) {
          zem2 = 1.0 - ze2 + 0.5 * ze2 * ze2;
          zep2 = 1.0 / zem2;
        } else {
          double tblind = ze2 / (bpade + ze2);
          int itind = static_cast<int>(tblint * tblind + 0.5);
          zem2 = tables.exp_tbl[itind];
          zep2 = 1.0 / zem2;
        }

        double zdenr = zr4 * zep1 + zr5 * zem1;
        double zdent = zt4 * zep1 + zt5 * zem1;
        if (zdenr >= -eps && zdenr <= eps) {
          pref[jk] = eps;
          ptra[jk] = zem2;
        } else {
          pref[jk] = zw * (zr1 * zep1 - zr2 * zem1 - zr3 * zem2) / zdenr;
          ptra[jk] = zem2 - zem2 * zw * (zt1 * zep1 - zt2 * zem1 - zt3 * zep2) / zdent;
        }

        double zemm = zem1 * zem1;
        double zdend = 1.0 / ((1.0 - zbeta * zemm) * zrkg);
        prefd[jk] = zgamma2 * (1.0 - zemm) * zdend;
        ptrad[jk] = zrk2 * zem1 * zdend;
      }
    }
  }
}

}  // namespace cppwrf::phys
