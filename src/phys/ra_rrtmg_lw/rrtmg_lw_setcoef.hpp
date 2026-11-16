// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "phys/ra_rrtmg_common/rrtmg_real.hpp"

#include "core/backend.hpp"
#include "phys/ra_rrtmg_common/strided.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_tables.hpp"
#include <algorithm>
#include <cmath>

namespace cppwrf::phys {

CPPWRF_HD inline rr get_chi_mls(const RrtmgLwTables& tables, int ispec, int ilev) {
  // ispec: 1..7, ilev: 1..59
  return tables.chi_mls[(ilev - 1) * 7 + (ispec - 1)];
}

// setcoef's statements for layer lay (WRF module_ra_rrtmg_lw.F setcoef, the body of its layer loop); returns
// whether the layer counts into laytrop (plog > 4.56)
CPPWRF_HD inline bool rrtmg_lw_setcoef_layer(
    int lay,
    int istart,
    SD pavel,
    SD tavel,
    SD tz,
    rr tbound,
    SD semiss,
    SD coldry,
    SD2 wkl, // [7][nlayers]
    SD wbroad,
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
  constexpr rr stpfac = rr(296.0) / rr(1013.0);

  int indbound = static_cast<int>(tbound - rr(159.0));
  if (indbound < 1) indbound = 1;
  else if (indbound > 180) indbound = 180;
  rr tbndfrac = tbound - rr(159.0) - static_cast<rr>(indbound);

  int indlev0 = static_cast<int>(tz[0] - rr(159.0));
  if (indlev0 < 1) indlev0 = 1;
  else if (indlev0 > 180) indlev0 = 180;
  rr t0frac = tz[0] - rr(159.0) - static_cast<rr>(indlev0);

  int indlay = static_cast<int>(tavel[lay] - rr(159.0));
  if (indlay < 1) indlay = 1;
  else if (indlay > 180) indlay = 180;
  rr tlayfrac = tavel[lay] - rr(159.0) - static_cast<rr>(indlay);

  int indlev = static_cast<int>(tz[lay + 1] - rr(159.0));
  if (indlev < 1) indlev = 1;
  else if (indlev > 180) indlev = 180;
  rr tlevfrac = tz[lay + 1] - rr(159.0) - static_cast<rr>(indlev);

  // Planck function evaluations across the 16 bands
  for (int iband = 0; iband < 15; ++iband) {
    if (lay == 0) {
      rr dbdtlev = tables.totplnk[iband * 181 + (indbound)] - tables.totplnk[iband * 181 + (indbound - 1)];
      plankbnd[iband] = semiss[iband] * (tables.totplnk[iband * 181 + (indbound - 1)] + tbndfrac * dbdtlev);

      dbdtlev = tables.totplnk[iband * 181 + (indlev0)] - tables.totplnk[iband * 181 + (indlev0 - 1)];
      planklev[0][iband] = tables.totplnk[iband * 181 + (indlev0 - 1)] + t0frac * dbdtlev;
    }
    rr dbdtlev = tables.totplnk[iband * 181 + (indlev)] - tables.totplnk[iband * 181 + (indlev - 1)];
    rr dbdtlay = tables.totplnk[iband * 181 + (indlay)] - tables.totplnk[iband * 181 + (indlay - 1)];
    planklay[lay][iband] = tables.totplnk[iband * 181 + (indlay - 1)] + tlayfrac * dbdtlay;
    planklev[lay + 1][iband] = tables.totplnk[iband * 181 + (indlev - 1)] + tlevfrac * dbdtlev;
  }

  // Band 16 (index 15)
  int iband = 15;
  if (istart == 16) {
    if (lay == 0) {
      rr dbdtlev = tables.totplk16[indbound] - tables.totplk16[indbound - 1];
      plankbnd[iband] = semiss[iband] * (tables.totplk16[indbound - 1] + tbndfrac * dbdtlev);

      dbdtlev = tables.totplnk[iband * 181 + (indlev0)] - tables.totplnk[iband * 181 + (indlev0 - 1)];
      planklev[0][iband] = tables.totplk16[indlev0 - 1] + t0frac * dbdtlev;
    }
    rr dbdtlev = tables.totplk16[indlev] - tables.totplk16[indlev - 1];
    rr dbdtlay = tables.totplk16[indlay] - tables.totplk16[indlay - 1];
    planklay[lay][iband] = tables.totplk16[indlay - 1] + tlayfrac * dbdtlay;
    planklev[lay + 1][iband] = tables.totplk16[indlev - 1] + tlevfrac * dbdtlev;
  } else {
    if (lay == 0) {
      rr dbdtlev = tables.totplnk[iband * 181 + (indbound)] - tables.totplnk[iband * 181 + (indbound - 1)];
      plankbnd[iband] = semiss[iband] * (tables.totplnk[iband * 181 + (indbound - 1)] + tbndfrac * dbdtlev);

      dbdtlev = tables.totplnk[iband * 181 + (indlev0)] - tables.totplnk[iband * 181 + (indlev0 - 1)];
      planklev[0][iband] = tables.totplnk[iband * 181 + (indlev0 - 1)] + t0frac * dbdtlev;
    }
    rr dbdtlev = tables.totplnk[iband * 181 + (indlev)] - tables.totplnk[iband * 181 + (indlev - 1)];
    rr dbdtlay = tables.totplnk[iband * 181 + (indlay)] - tables.totplnk[iband * 181 + (indlay - 1)];
    planklay[lay][iband] = tables.totplnk[iband * 181 + (indlay - 1)] + tlayfrac * dbdtlay;
    planklev[lay + 1][iband] = tables.totplnk[iband * 181 + (indlev - 1)] + tlevfrac * dbdtlev;
  }

  // Pressure and temperature interpolation factors
  rr plog = libm::log(pavel[lay]);
  jp[lay] = static_cast<int>(rr(36.0) - rr(5.0) * (plog + rr(0.04)));
  if (jp[lay] < 1) jp[lay] = 1;
  else if (jp[lay] > 58) jp[lay] = 58;
  int jp1 = jp[lay] + 1;

  rr fp = rr(5.0) * (tables.preflog[jp[lay] - 1] - plog);

  jt[lay] = static_cast<int>(rr(3.0) + (tavel[lay] - tables.tref[jp[lay] - 1]) / rr(15.0));
  if (jt[lay] < 1) jt[lay] = 1;
  else if (jt[lay] > 4) jt[lay] = 4;
  rr ft = ((tavel[lay] - tables.tref[jp[lay] - 1]) / rr(15.0)) - static_cast<rr>(jt[lay] - 3);

  jt1[lay] = static_cast<int>(rr(3.0) + (tavel[lay] - tables.tref[jp1 - 1]) / rr(15.0));
  if (jt1[lay] < 1) jt1[lay] = 1;
  else if (jt1[lay] > 4) jt1[lay] = 4;
  rr ft1 = ((tavel[lay] - tables.tref[jp1 - 1]) / rr(15.0)) - static_cast<rr>(jt1[lay] - 3);

  rr water = wkl[0][lay] / coldry[lay];
  rr scalefac = pavel[lay] * stpfac / tavel[lay];

  if (plog > rr(4.56)) {
    // Lower atmosphere

    forfac[lay] = scalefac / (rr(1.0) + water);
    rr factor = (rr(332.0) - tavel[lay]) / rr(36.0);
    indfor[lay] = std::min(2, std::max(1, static_cast<int>(factor)));
    forfrac[lay] = factor - static_cast<rr>(indfor[lay]);

    selffac[lay] = water * forfac[lay];
    factor = (tavel[lay] - rr(188.0)) / rr(7.2);
    indself[lay] = std::min(9, std::max(1, static_cast<int>(factor) - 7));
    selffrac[lay] = factor - static_cast<rr>(indself[lay] + 7);

    scaleminor[lay] = pavel[lay] / tavel[lay];
    scaleminorn2[lay] = (pavel[lay] / tavel[lay]) * (wbroad[lay] / (coldry[lay] + wkl[0][lay]));
    factor = (tavel[lay] - rr(180.8)) / rr(7.2);
    indminor[lay] = std::min(18, std::max(1, static_cast<int>(factor)));
    minorfrac[lay] = factor - static_cast<rr>(indminor[lay]);

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

    colh2o[lay] = rr(1.0e-20) * wkl[0][lay];
    colco2[lay] = rr(1.0e-20) * wkl[1][lay];
    colo3[lay]  = rr(1.0e-20) * wkl[2][lay];
    coln2o[lay] = rr(1.0e-20) * wkl[3][lay];
    colco[lay]  = rr(1.0e-20) * wkl[4][lay];
    colch4[lay] = rr(1.0e-20) * wkl[5][lay];
    colo2[lay]  = rr(1.0e-20) * wkl[6][lay];
    if (colco2[lay] == rr(0.0)) colco2[lay] = rr(1.0e-32) * coldry[lay];
    if (colo3[lay]  == rr(0.0)) colo3[lay]  = rr(1.0e-32) * coldry[lay];
    if (coln2o[lay] == rr(0.0)) coln2o[lay] = rr(1.0e-32) * coldry[lay];
    if (colco[lay]  == rr(0.0)) colco[lay]  = rr(1.0e-32) * coldry[lay];
    if (colch4[lay] == rr(0.0)) colch4[lay] = rr(1.0e-32) * coldry[lay];
    colbrd[lay] = rr(1.0e-20) * wbroad[lay];
  } else {
    // Upper atmosphere
    forfac[lay] = scalefac / (rr(1.0) + water);
    rr factor = (tavel[lay] - rr(188.0)) / rr(36.0);
    indfor[lay] = 3;
    forfrac[lay] = factor - rr(1.0);

    selffac[lay] = water * forfac[lay];

    scaleminor[lay] = pavel[lay] / tavel[lay];
    scaleminorn2[lay] = (pavel[lay] / tavel[lay]) * (wbroad[lay] / (coldry[lay] + wkl[0][lay]));
    factor = (tavel[lay] - rr(180.8)) / rr(7.2);
    indminor[lay] = std::min(18, std::max(1, static_cast<int>(factor)));
    minorfrac[lay] = factor - static_cast<rr>(indminor[lay]);

    rat_h2oco2[lay] = get_chi_mls(tables, 1, jp[lay]) / get_chi_mls(tables, 2, jp[lay]);
    rat_h2oco2_1[lay] = get_chi_mls(tables, 1, jp[lay] + 1) / get_chi_mls(tables, 2, jp[lay] + 1);

    rat_o3co2[lay] = get_chi_mls(tables, 3, jp[lay]) / get_chi_mls(tables, 2, jp[lay]);
    rat_o3co2_1[lay] = get_chi_mls(tables, 3, jp[lay] + 1) / get_chi_mls(tables, 2, jp[lay] + 1);

    colh2o[lay] = rr(1.0e-20) * wkl[0][lay];
    colco2[lay] = rr(1.0e-20) * wkl[1][lay];
    colo3[lay]  = rr(1.0e-20) * wkl[2][lay];
    coln2o[lay] = rr(1.0e-20) * wkl[3][lay];
    colco[lay]  = rr(1.0e-20) * wkl[4][lay];
    colch4[lay] = rr(1.0e-20) * wkl[5][lay];
    colo2[lay]  = rr(1.0e-20) * wkl[6][lay];
    if (colco2[lay] == rr(0.0)) colco2[lay] = rr(1.0e-32) * coldry[lay];
    if (colo3[lay]  == rr(0.0)) colo3[lay]  = rr(1.0e-32) * coldry[lay];
    if (coln2o[lay] == rr(0.0)) coln2o[lay] = rr(1.0e-32) * coldry[lay];
    if (colco[lay]  == rr(0.0)) colco[lay]  = rr(1.0e-32) * coldry[lay];
    if (colch4[lay] == rr(0.0)) colch4[lay] = rr(1.0e-32) * coldry[lay];
    colbrd[lay] = rr(1.0e-20) * wbroad[lay];
  }

  rr compfp = rr(1.0) - fp;
  fac10[lay] = compfp * ft;
  fac00[lay] = compfp * (rr(1.0) - ft);
  fac11[lay] = fp * ft1;
  fac01[lay] = fp * (rr(1.0) - ft1);

  selffac[lay] *= colh2o[lay];
  forfac[lay] *= colh2o[lay];
  return plog > rr(4.56);
}

CPPWRF_HD inline void rrtmg_lw_setcoef(
    int nlayers,
    int istart,
    SD pavel,
    SD tavel,
    SD tz,
    rr tbound,
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
  laytrop = 0;
  for (int lay = 0; lay < nlayers; ++lay) {
    if (rrtmg_lw_setcoef_layer(lay, istart, pavel, tavel, tz, tbound, semiss, coldry, wkl, wbroad, jp, jt, jt1, planklay,
                               planklev, plankbnd, colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00,
                               fac01, fac10, fac11, rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o,
                               rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1,
                               selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2,
                               indminor, tables))
      ++laytrop;
  }
}

}  // namespace cppwrf::phys
