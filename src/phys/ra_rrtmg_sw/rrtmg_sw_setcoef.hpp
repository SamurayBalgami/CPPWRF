#pragma once

#include "core/backend.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_tables.hpp"
#include <algorithm>
#include <cmath>

namespace cppwrf::phys {

CPPWRF_HD inline void setcoef_sw(
    int nlayers,
    const double* pavel,      // (nlayers)
    const double* tavel,      // (nlayers)
    const double* coldry,     // (nlayers)
    const double* const* wkl, // wkl[imol][lay], 1-based in WRF: wkl[0]=h2o, wkl[1]=co2, wkl[2]=o3, wkl[3]=n2o, wkl[5]=ch4, wkl[6]=o2
    int& laytrop,
    int& layswtch,
    int& laylow,
    int* jp,                  // (nlayers)
    int* jt,                  // (nlayers)
    int* jt1,                 // (nlayers)
    double* colh2o,           // (nlayers)
    double* colco2,           // (nlayers)
    double* colo3,            // (nlayers)
    double* coln2o,           // (nlayers)
    double* colch4,           // (nlayers)
    double* colo2,            // (nlayers)
    double* colmol,           // (nlayers)
    double* co2mult,          // (nlayers)
    int* indself,             // (nlayers)
    int* indfor,              // (nlayers)
    double* selffac,          // (nlayers)
    double* selffrac,         // (nlayers)
    double* forfac,           // (nlayers)
    double* forfrac,          // (nlayers)
    double* fac00,            // (nlayers)
    double* fac01,            // (nlayers)
    double* fac10,            // (nlayers)
    double* fac11,            // (nlayers)
    const RrtmgSwTables& tables
) {
  constexpr double stpfac = 296.0 / 1013.0;

  laytrop = 0;
  layswtch = 0;
  laylow = 0;

  for (int lay = 0; lay < nlayers; ++lay) {
    double plog = std::log(pavel[lay]);
    jp[lay] = static_cast<int>(36.0 - 5.0 * (plog + 0.04));
    if (jp[lay] < 1) {
      jp[lay] = 1;
    } else if (jp[lay] > 58) {
      jp[lay] = 58;
    }
    int jp1 = jp[lay] + 1;
    double fp = 5.0 * (tables.preflog[jp[lay] - 1] - plog);

    jt[lay] = static_cast<int>(3.0 + (tavel[lay] - tables.tref[jp[lay] - 1]) / 15.0);
    if (jt[lay] < 1) {
      jt[lay] = 1;
    } else if (jt[lay] > 4) {
      jt[lay] = 4;
    }
    double ft = ((tavel[lay] - tables.tref[jp[lay] - 1]) / 15.0) - static_cast<double>(jt[lay] - 3);

    jt1[lay] = static_cast<int>(3.0 + (tavel[lay] - tables.tref[jp1 - 1]) / 15.0);
    if (jt1[lay] < 1) {
      jt1[lay] = 1;
    } else if (jt1[lay] > 4) {
      jt1[lay] = 4;
    }
    double ft1 = ((tavel[lay] - tables.tref[jp1 - 1]) / 15.0) - static_cast<double>(jt1[lay] - 3);

    double water = wkl[0][lay] / coldry[lay];
    double scalefac = pavel[lay] * stpfac / tavel[lay];

    colh2o[lay] = 1.0e-20 * wkl[0][lay];
    colco2[lay] = 1.0e-20 * wkl[1][lay];
    colo3[lay]  = 1.0e-20 * wkl[2][lay];
    coln2o[lay] = 1.0e-20 * wkl[3][lay];
    colch4[lay] = 1.0e-20 * wkl[5][lay];
    colo2[lay]  = 1.0e-20 * wkl[6][lay];
    colmol[lay] = 1.0e-20 * coldry[lay] + colh2o[lay];

    if (colco2[lay] == 0.0) colco2[lay] = 1.0e-32 * coldry[lay];
    if (coln2o[lay] == 0.0) coln2o[lay] = 1.0e-32 * coldry[lay];
    if (colch4[lay] == 0.0) colch4[lay] = 1.0e-32 * coldry[lay];
    if (colo2[lay] == 0.0)  colo2[lay]  = 1.0e-32 * coldry[lay];

    double co2reg = 3.55e-24 * coldry[lay];
    co2mult[lay] = (colco2[lay] - co2reg) * 272.63 * std::exp(-1919.4 / tavel[lay]) / (8.7604e-4 * tavel[lay]);

    if (plog > 4.56) {
      // Lower atmosphere
      laytrop++;
      if (plog >= 6.62) laylow++;

      forfac[lay] = scalefac / (1.0 + water);
      double factor = (332.0 - tavel[lay]) / 36.0;
      indfor[lay] = std::min(2, std::max(1, static_cast<int>(factor)));
      forfrac[lay] = factor - static_cast<double>(indfor[lay]);

      selffac[lay] = water * forfac[lay];
      factor = (tavel[lay] - 188.0) / 7.2;
      indself[lay] = std::min(9, std::max(1, static_cast<int>(factor) - 7));
      selffrac[lay] = factor - static_cast<double>(indself[lay] + 7);
    } else {
      // Upper atmosphere
      forfac[lay] = scalefac / (1.0 + water);
      double factor = (tavel[lay] - 188.0) / 36.0;
      indfor[lay] = 3;
      forfrac[lay] = factor - 1.0;

      selffac[lay] = 0.0;
      selffrac[lay] = 0.0;
      indself[lay] = 0;
    }

    double compfp = 1.0 - fp;
    fac10[lay] = compfp * ft;
    fac00[lay] = compfp * (1.0 - ft);
    fac11[lay] = fp * ft1;
    fac01[lay] = fp * (1.0 - ft1);
  }
}

}  // namespace cppwrf::phys
