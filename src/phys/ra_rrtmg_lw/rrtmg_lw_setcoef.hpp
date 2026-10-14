#pragma once

#include "core/backend.hpp"
#include "phys/ra_rrtmg_common/strided.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_tables.hpp"
#include <algorithm>
#include <cmath>

namespace cppwrf::phys {

CPPWRF_HD inline double get_chi_mls(const RrtmgLwTables& tables, int ispec, int ilev) {
  // ispec: 1..7, ilev: 1..59
  return tables.chi_mls[(ilev - 1) * 7 + (ispec - 1)];
}

CPPWRF_HD inline void rrtmg_lw_setcoef(
    int nlayers,
    int istart,
    SD pavel,
    SD tavel,
    SD tz,
    double tbound,
    SD semiss,
    SD coldry,
    SD2 wkl, // [7][nlayers]
    SD wbroad,
    int& laytrop,
    SI jp,
    SI jt,
    SI jt1,
    SD2 planklay,  // [nlayers][16]
    SD2 planklev,  // [nlayers + 1][16]
    SD plankbnd,         // [16]
    SD colh2o,
    SD colco2,
    SD colo3,
    SD coln2o,
    SD colco,
    SD colch4,
    SD colo2,
    SD colbrd,
    SD fac00,
    SD fac01,
    SD fac10,
    SD fac11,
    SD rat_h2oco2,
    SD rat_h2oco2_1,
    SD rat_h2oo3,
    SD rat_h2oo3_1,
    SD rat_h2on2o,
    SD rat_h2on2o_1,
    SD rat_h2och4,
    SD rat_h2och4_1,
    SD rat_n2oco2,
    SD rat_n2oco2_1,
    SD rat_o3co2,
    SD rat_o3co2_1,
    SD selffac,
    SD selffrac,
    SI indself,
    SD forfac,
    SD forfrac,
    SI indfor,
    SD minorfrac,
    SD scaleminor,
    SD scaleminorn2,
    SI indminor,
    const RrtmgLwTables& tables
) {
  constexpr double stpfac = 296.0 / 1013.0;

  int indbound = static_cast<int>(tbound - 159.0);
  if (indbound < 1) indbound = 1;
  else if (indbound > 180) indbound = 180;
  double tbndfrac = tbound - 159.0 - static_cast<double>(indbound);

  int indlev0 = static_cast<int>(tz[0] - 159.0);
  if (indlev0 < 1) indlev0 = 1;
  else if (indlev0 > 180) indlev0 = 180;
  double t0frac = tz[0] - 159.0 - static_cast<double>(indlev0);

  laytrop = 0;

  for (int lay = 0; lay < nlayers; ++lay) {
    int indlay = static_cast<int>(tavel[lay] - 159.0);
    if (indlay < 1) indlay = 1;
    else if (indlay > 180) indlay = 180;
    double tlayfrac = tavel[lay] - 159.0 - static_cast<double>(indlay);

    int indlev = static_cast<int>(tz[lay + 1] - 159.0);
    if (indlev < 1) indlev = 1;
    else if (indlev > 180) indlev = 180;
    double tlevfrac = tz[lay + 1] - 159.0 - static_cast<double>(indlev);

    // Planck function evaluations across the 16 bands
    for (int iband = 0; iband < 15; ++iband) {
      if (lay == 0) {
        double dbdtlev = tables.totplnk[iband * 181 + (indbound)] - tables.totplnk[iband * 181 + (indbound - 1)];
        plankbnd[iband] = semiss[iband] * (tables.totplnk[iband * 181 + (indbound - 1)] + tbndfrac * dbdtlev);

        dbdtlev = tables.totplnk[iband * 181 + (indlev0)] - tables.totplnk[iband * 181 + (indlev0 - 1)];
        planklev[0][iband] = tables.totplnk[iband * 181 + (indlev0 - 1)] + t0frac * dbdtlev;
      }
      double dbdtlev = tables.totplnk[iband * 181 + (indlev)] - tables.totplnk[iband * 181 + (indlev - 1)];
      double dbdtlay = tables.totplnk[iband * 181 + (indlay)] - tables.totplnk[iband * 181 + (indlay - 1)];
      planklay[lay][iband] = tables.totplnk[iband * 181 + (indlay - 1)] + tlayfrac * dbdtlay;
      planklev[lay + 1][iband] = tables.totplnk[iband * 181 + (indlev - 1)] + tlevfrac * dbdtlev;
    }

    // Band 16 (index 15)
    int iband = 15;
    if (istart == 16) {
      if (lay == 0) {
        double dbdtlev = tables.totplk16[indbound] - tables.totplk16[indbound - 1];
        plankbnd[iband] = semiss[iband] * (tables.totplk16[indbound - 1] + tbndfrac * dbdtlev);

        dbdtlev = tables.totplnk[iband * 181 + (indlev0)] - tables.totplnk[iband * 181 + (indlev0 - 1)];
        planklev[0][iband] = tables.totplk16[indlev0 - 1] + t0frac * dbdtlev;
      }
      double dbdtlev = tables.totplk16[indlev] - tables.totplk16[indlev - 1];
      double dbdtlay = tables.totplk16[indlay] - tables.totplk16[indlay - 1];
      planklay[lay][iband] = tables.totplk16[indlay - 1] + tlayfrac * dbdtlay;
      planklev[lay + 1][iband] = tables.totplk16[indlev - 1] + tlevfrac * dbdtlev;
    } else {
      if (lay == 0) {
        double dbdtlev = tables.totplnk[iband * 181 + (indbound)] - tables.totplnk[iband * 181 + (indbound - 1)];
        plankbnd[iband] = semiss[iband] * (tables.totplnk[iband * 181 + (indbound - 1)] + tbndfrac * dbdtlev);

        dbdtlev = tables.totplnk[iband * 181 + (indlev0)] - tables.totplnk[iband * 181 + (indlev0 - 1)];
        planklev[0][iband] = tables.totplnk[iband * 181 + (indlev0 - 1)] + t0frac * dbdtlev;
      }
      double dbdtlev = tables.totplnk[iband * 181 + (indlev)] - tables.totplnk[iband * 181 + (indlev - 1)];
      double dbdtlay = tables.totplnk[iband * 181 + (indlay)] - tables.totplnk[iband * 181 + (indlay - 1)];
      planklay[lay][iband] = tables.totplnk[iband * 181 + (indlay - 1)] + tlayfrac * dbdtlay;
      planklev[lay + 1][iband] = tables.totplnk[iband * 181 + (indlev - 1)] + tlevfrac * dbdtlev;
    }

    // Pressure and temperature interpolation factors
    double plog = std::log(pavel[lay]);
    jp[lay] = static_cast<int>(36.0 - 5.0 * (plog + 0.04));
    if (jp[lay] < 1) jp[lay] = 1;
    else if (jp[lay] > 58) jp[lay] = 58;
    int jp1 = jp[lay] + 1;

    double fp = 5.0 * (tables.preflog[jp[lay] - 1] - plog);

    jt[lay] = static_cast<int>(3.0 + (tavel[lay] - tables.tref[jp[lay] - 1]) / 15.0);
    if (jt[lay] < 1) jt[lay] = 1;
    else if (jt[lay] > 4) jt[lay] = 4;
    double ft = ((tavel[lay] - tables.tref[jp[lay] - 1]) / 15.0) - static_cast<double>(jt[lay] - 3);

    jt1[lay] = static_cast<int>(3.0 + (tavel[lay] - tables.tref[jp1 - 1]) / 15.0);
    if (jt1[lay] < 1) jt1[lay] = 1;
    else if (jt1[lay] > 4) jt1[lay] = 4;
    double ft1 = ((tavel[lay] - tables.tref[jp1 - 1]) / 15.0) - static_cast<double>(jt1[lay] - 3);

    double water = wkl[0][lay] / coldry[lay];
    double scalefac = pavel[lay] * stpfac / tavel[lay];

    if (plog > 4.56) {
      // Lower atmosphere
      laytrop++;

      forfac[lay] = scalefac / (1.0 + water);
      double factor = (332.0 - tavel[lay]) / 36.0;
      indfor[lay] = std::min(2, std::max(1, static_cast<int>(factor)));
      forfrac[lay] = factor - static_cast<double>(indfor[lay]);

      selffac[lay] = water * forfac[lay];
      factor = (tavel[lay] - 188.0) / 7.2;
      indself[lay] = std::min(9, std::max(1, static_cast<int>(factor) - 7));
      selffrac[lay] = factor - static_cast<double>(indself[lay] + 7);

      scaleminor[lay] = pavel[lay] / tavel[lay];
      scaleminorn2[lay] = (pavel[lay] / tavel[lay]) * (wbroad[lay] / (coldry[lay] + wkl[0][lay]));
      factor = (tavel[lay] - 180.8) / 7.2;
      indminor[lay] = std::min(18, std::max(1, static_cast<int>(factor)));
      minorfrac[lay] = factor - static_cast<double>(indminor[lay]);

      rat_h2oco2[lay] = get_chi_mls(tables, 1, jp[lay]) / get_chi_mls(tables, 2, jp[lay]);
      rat_h2oco2_1[lay] = get_chi_mls(tables, 1, jp[lay] + 1) / get_chi_mls(tables, 2, jp[lay] + 1);

      rat_h2oo3[lay] = get_chi_mls(tables, 1, jp[lay]) / get_chi_mls(tables, 3, jp[lay]);
      rat_h2oo3_1[lay] = get_chi_mls(tables, 1, jp[lay] + 1) / get_chi_mls(tables, 3, jp[lay] + 1);

      rat_h2on2o[lay] = get_chi_mls(tables, 1, jp[lay]) / get_chi_mls(tables, 4, jp[lay]);
      rat_h2on2o_1[lay] = get_chi_mls(tables, 1, jp[lay] + 1) / get_chi_mls(tables, 4, jp[lay] + 1);

      rat_h2och4[lay] = get_chi_mls(tables, 1, jp[lay]) / get_chi_mls(tables, 6, jp[lay]);
      rat_h2och4_1[lay] = get_chi_mls(tables, 1, jp[lay] + 1) / get_chi_mls(tables, 6, jp[lay] + 1);

      rat_n2oco2[lay] = get_chi_mls(tables, 4, jp[lay]) / get_chi_mls(tables, 2, jp[lay]);
      rat_n2oco2_1[lay] = get_chi_mls(tables, 4, jp[lay] + 1) / get_chi_mls(tables, 2, jp[lay] + 1);

      colh2o[lay] = 1.0e-20 * wkl[0][lay];
      colco2[lay] = 1.0e-20 * wkl[1][lay];
      colo3[lay]  = 1.0e-20 * wkl[2][lay];
      coln2o[lay] = 1.0e-20 * wkl[3][lay];
      colco[lay]  = 1.0e-20 * wkl[4][lay];
      colch4[lay] = 1.0e-20 * wkl[5][lay];
      colo2[lay]  = 1.0e-20 * wkl[6][lay];
      if (colco2[lay] == 0.0) colco2[lay] = 1.0e-32 * coldry[lay];
      if (colo3[lay]  == 0.0) colo3[lay]  = 1.0e-32 * coldry[lay];
      if (coln2o[lay] == 0.0) coln2o[lay] = 1.0e-32 * coldry[lay];
      if (colco[lay]  == 0.0) colco[lay]  = 1.0e-32 * coldry[lay];
      if (colch4[lay] == 0.0) colch4[lay] = 1.0e-32 * coldry[lay];
      colbrd[lay] = 1.0e-20 * wbroad[lay];
    } else {
      // Upper atmosphere
      forfac[lay] = scalefac / (1.0 + water);
      double factor = (tavel[lay] - 188.0) / 36.0;
      indfor[lay] = 3;
      forfrac[lay] = factor - 1.0;

      selffac[lay] = water * forfac[lay];

      scaleminor[lay] = pavel[lay] / tavel[lay];
      scaleminorn2[lay] = (pavel[lay] / tavel[lay]) * (wbroad[lay] / (coldry[lay] + wkl[0][lay]));
      factor = (tavel[lay] - 180.8) / 7.2;
      indminor[lay] = std::min(18, std::max(1, static_cast<int>(factor)));
      minorfrac[lay] = factor - static_cast<double>(indminor[lay]);

      rat_h2oco2[lay] = get_chi_mls(tables, 1, jp[lay]) / get_chi_mls(tables, 2, jp[lay]);
      rat_h2oco2_1[lay] = get_chi_mls(tables, 1, jp[lay] + 1) / get_chi_mls(tables, 2, jp[lay] + 1);

      rat_o3co2[lay] = get_chi_mls(tables, 3, jp[lay]) / get_chi_mls(tables, 2, jp[lay]);
      rat_o3co2_1[lay] = get_chi_mls(tables, 3, jp[lay] + 1) / get_chi_mls(tables, 2, jp[lay] + 1);

      colh2o[lay] = 1.0e-20 * wkl[0][lay];
      colco2[lay] = 1.0e-20 * wkl[1][lay];
      colo3[lay]  = 1.0e-20 * wkl[2][lay];
      coln2o[lay] = 1.0e-20 * wkl[3][lay];
      colco[lay]  = 1.0e-20 * wkl[4][lay];
      colch4[lay] = 1.0e-20 * wkl[5][lay];
      colo2[lay]  = 1.0e-20 * wkl[6][lay];
      if (colco2[lay] == 0.0) colco2[lay] = 1.0e-32 * coldry[lay];
      if (colo3[lay]  == 0.0) colo3[lay]  = 1.0e-32 * coldry[lay];
      if (coln2o[lay] == 0.0) coln2o[lay] = 1.0e-32 * coldry[lay];
      if (colco[lay]  == 0.0) colco[lay]  = 1.0e-32 * coldry[lay];
      if (colch4[lay] == 0.0) colch4[lay] = 1.0e-32 * coldry[lay];
      colbrd[lay] = 1.0e-20 * wbroad[lay];
    }

    double compfp = 1.0 - fp;
    fac10[lay] = compfp * ft;
    fac00[lay] = compfp * (1.0 - ft);
    fac11[lay] = fp * ft1;
    fac01[lay] = fp * (1.0 - ft1);

    selffac[lay] *= colh2o[lay];
    forfac[lay] *= colh2o[lay];
  }
}

}  // namespace cppwrf::phys
