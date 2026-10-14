#pragma once

#include "core/backend.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_tables.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_reftra.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_vrtqdr.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_taumol.hpp"
#include <algorithm>
#include <cmath>

namespace cppwrf::phys {

CPPWRF_HD inline void spcvmc_sw(
    int nlayers,
    int icpr,
    const double* palbd,      // (14)
    const double* palbp,      // (14)
    const double* pcldfmc,    // (nlayers, 112)
    const double* ptaucmc,    // (nlayers, 112)
    const double* pasycmc,    // (nlayers, 112)
    const double* pomgcmc,    // (nlayers, 112)
    const double* ptaormc,    // (nlayers, 112)
    const double* ptaua,      // (nlayers, 14)
    const double* pasya,      // (nlayers, 14)
    const double* pomga,      // (nlayers, 14)
    double prmu0,
    const double* adjflux,    // (30)
    int laytrop,
    int layswtch,
    int laylow,
    const int* jp,
    const int* jt,
    const int* jt1,
    const double* colh2o,
    const double* colco2,
    const double* colch4,
    const double* colo2,
    const double* colo3,
    const double* colmol,
    const double* fac00,
    const double* fac01,
    const double* fac10,
    const double* fac11,
    const double* selffac,
    const double* selffrac,
    const int* indself,
    const double* forfac,
    const double* forfrac,
    const int* indfor,
    double* pbbfd,            // (nlayers+1)
    double* pbbfu,            // (nlayers+1)
    double* pbbcd,            // (nlayers+1)
    double* pbbcu,            // (nlayers+1)
    double* puvfd,            // (nlayers+1)
    double* puvcd,            // (nlayers+1)
    double* pnifd,            // (nlayers+1)
    double* pnicd,            // (nlayers+1)
    double* pbbfddir,         // (nlayers+1)
    double* pbbcddir,         // (nlayers+1)
    double* puvfddir,         // (nlayers+1)
    double* puvcddir,         // (nlayers+1)
    double* pnifddir,         // (nlayers+1)
    double* pnicddir,         // (nlayers+1)
    const RrtmgSwTables& tables
) {
  constexpr double tblint = 10000.0;
  constexpr double od_lo = 0.06;
  constexpr double pade = 0.278;
  constexpr double bpade = 1.0 / pade;
  constexpr double repclc = 1.0e-12;

  int klev = nlayers;

  // Working buffers for 2-stream calculations
  double ztaug[MAX_SW_LAYERS][112];
  double ztaur[MAX_SW_LAYERS][112];
  double zsflxzen[112];

  taumol_sw(
      nlayers, laytrop, jp, jt, jt1,
      colh2o, colco2, colch4, colo2, colo3, colmol,
      fac00, fac01, fac10, fac11,
      selffac, selffrac, indself, forfac, forfrac, indfor,
      ztaug, ztaur, zsflxzen, tables
  );

  for (int i = 0; i <= klev; ++i) {
    pbbfd[i] = 0.0; pbbfu[i] = 0.0;
    pbbcd[i] = 0.0; pbbcu[i] = 0.0;
    puvfd[i] = 0.0; puvcd[i] = 0.0;
    pnifd[i] = 0.0; pnicd[i] = 0.0;
    pbbfddir[i] = 0.0; pbbcddir[i] = 0.0;
    puvfddir[i] = 0.0; puvcddir[i] = 0.0;
    pnifddir[i] = 0.0; pnicddir[i] = 0.0;
  }

  // Intermediate layer arrays
  double ztauc[MAX_SW_LAYERS], zomcc[MAX_SW_LAYERS], zgcc[MAX_SW_LAYERS];
  double ztauo[MAX_SW_LAYERS], zomco[MAX_SW_LAYERS], zgco[MAX_SW_LAYERS];
  bool lrtchkclr[MAX_SW_LAYERS], lrtchkcld[MAX_SW_LAYERS];

  double zrefc[MAX_SW_LAYERS + 1], zrefdc[MAX_SW_LAYERS + 1];
  double ztrac[MAX_SW_LAYERS + 1], ztradc[MAX_SW_LAYERS + 1];
  double zrefo[MAX_SW_LAYERS + 1], zrefdo[MAX_SW_LAYERS + 1];
  double ztrao[MAX_SW_LAYERS + 1], ztrado[MAX_SW_LAYERS + 1];
  double zdbtc[MAX_SW_LAYERS + 1], zdbtc_nodel[MAX_SW_LAYERS + 1];
  double ztdbtc[MAX_SW_LAYERS + 1], ztdbtc_nodel[MAX_SW_LAYERS + 1];
  double zrdndc[MAX_SW_LAYERS + 1], zrupc[MAX_SW_LAYERS + 1], zrupdc[MAX_SW_LAYERS + 1];

  double zref[MAX_SW_LAYERS + 1], zrefd[MAX_SW_LAYERS + 1];
  double ztra[MAX_SW_LAYERS + 1], ztrad[MAX_SW_LAYERS + 1];
  double zdbt[MAX_SW_LAYERS + 1], zdbt_nodel[MAX_SW_LAYERS + 1];
  double ztdbt[MAX_SW_LAYERS + 1], ztdbt_nodel[MAX_SW_LAYERS + 1];
  double zrdnd[MAX_SW_LAYERS + 1], zrup[MAX_SW_LAYERS + 1], zrupd[MAX_SW_LAYERS + 1];

  double zcd[MAX_SW_LAYERS + 1], zcu[MAX_SW_LAYERS + 1];
  double zfd[MAX_SW_LAYERS + 1], zfu[MAX_SW_LAYERS + 1];

  int iw = 0;
  // Top of shortwave spectral band loop (jb = 16 to 29; ibm = 0 to 13)
  for (int jb = 16; jb <= 29; ++jb) {
    int ibm = jb - 16;
    int igt = tables.ngc[ibm];

    for (int jg = 0; jg < igt; ++jg) {
      double zincflx = adjflux[jb] * zsflxzen[iw] * prmu0;

      // Clear-sky TOA direct beam & surface values
      ztdbtc[0] = 1.0;
      ztdbtc_nodel[0] = 1.0;
      zdbtc[klev] = 0.0;
      ztrac[klev] = 0.0;
      ztradc[klev] = 0.0;
      zrefc[klev] = palbp[ibm];
      zrefdc[klev] = palbd[ibm];
      zrupc[klev] = palbp[ibm];
      zrupdc[klev] = palbd[ibm];

      // Total-sky TOA direct beam & surface values
      ztdbt[0] = 1.0;
      ztdbt_nodel[0] = 1.0;
      zdbt[klev] = 0.0;
      ztra[klev] = 0.0;
      ztrad[klev] = 0.0;
      zref[klev] = palbp[ibm];
      zrefd[klev] = palbd[ibm];
      zrup[klev] = palbp[ibm];
      zrupd[klev] = palbd[ibm];

      // Top of layer loop
      for (int jk = 0; jk < klev; ++jk) {
        int ikl = klev - 1 - jk; // reverse bottom-to-top into top-to-bottom

        lrtchkclr[jk] = true;
        lrtchkcld[jk] = (pcldfmc[ikl * 112 + iw] > repclc);

        ztauc[jk] = ztaur[ikl][iw] + ztaug[ikl][iw] + ptaua[ikl * 14 + ibm];
        zomcc[jk] = ztaur[ikl][iw] * 1.0 + ptaua[ikl * 14 + ibm] * pomga[ikl * 14 + ibm];
        zgcc[jk]  = (zomcc[jk] != 0.0) ? pasya[ikl * 14 + ibm] * pomga[ikl * 14 + ibm] * ptaua[ikl * 14 + ibm] / zomcc[jk] : 0.0;
        zomcc[jk] = (ztauc[jk] != 0.0) ? zomcc[jk] / ztauc[jk] : 0.0;

        double zclear = 1.0 - pcldfmc[ikl * 112 + iw];
        double zcloud = pcldfmc[ikl * 112 + iw];

        // Clear direct beam
        double ze1 = ztauc[jk] / prmu0;
        double zdbtmc;
        if (ze1 <= od_lo) {
          zdbtmc = 1.0 - ze1 + 0.5 * ze1 * ze1;
        } else {
          double tblind = ze1 / (bpade + ze1);
          int itind = static_cast<int>(tblint * tblind + 0.5);
          zdbtmc = tables.exp_tbl[itind];
        }
        zdbtc_nodel[jk] = zdbtmc;
        ztdbtc_nodel[jk + 1] = zdbtc_nodel[jk] * ztdbtc_nodel[jk];

        // Clear + cloud direct beam
        double tauorig = ztauc[jk] + ptaormc[ikl * 112 + iw];
        ze1 = tauorig / prmu0;
        double zdbtmo;
        if (ze1 <= od_lo) {
          zdbtmo = 1.0 - ze1 + 0.5 * ze1 * ze1;
        } else {
          double tblind = ze1 / (bpade + ze1);
          int itind = static_cast<int>(tblint * tblind + 0.5);
          zdbtmo = tables.exp_tbl[itind];
        }
        zdbt_nodel[jk] = zclear * zdbtmc + zcloud * zdbtmo;
        ztdbt_nodel[jk + 1] = zdbt_nodel[jk] * ztdbt_nodel[jk];
      }

      // Delta scaling - clear
      for (int jk = 0; jk < klev; ++jk) {
        double zf = zgcc[jk] * zgcc[jk];
        double zwf = zomcc[jk] * zf;
        ztauc[jk] = (1.0 - zwf) * ztauc[jk];
        zomcc[jk] = (zomcc[jk] - zwf) / (1.0 - zwf);
        zgcc[jk]  = (zgcc[jk] - zf) / (1.0 - zf);
      }

      // Total sky optical parameters
      if (icpr >= 1) {
        for (int jk = 0; jk < klev; ++jk) {
          int ikl = klev - 1 - jk;
          ztauo[jk] = ztauc[jk] + ptaucmc[ikl * 112 + iw];
          zomco[jk] = ztauc[jk] * zomcc[jk] + ptaucmc[ikl * 112 + iw] * pomgcmc[ikl * 112 + iw];
          zgco[jk]  = (ptaucmc[ikl * 112 + iw] * pomgcmc[ikl * 112 + iw] * pasycmc[ikl * 112 + iw] +
                       ztauc[jk] * zomcc[jk] * zgcc[jk]) / zomco[jk];
          zomco[jk] = zomco[jk] / ztauo[jk];
        }
      }

      // Clear sky reflectivities
      reftra_sw(klev, lrtchkclr, zgcc, prmu0, ztauc, zomcc, zrefc, zrefdc, ztrac, ztradc, tables);

      // Total sky reflectivities
      reftra_sw(klev, lrtchkcld, zgco, prmu0, ztauo, zomco, zrefo, zrefdo, ztrao, ztrado, tables);

      for (int jk = 0; jk < klev; ++jk) {
        int ikl = klev - 1 - jk;
        double zclear = 1.0 - pcldfmc[ikl * 112 + iw];
        double zcloud = pcldfmc[ikl * 112 + iw];

        zref[jk]  = zclear * zrefc[jk]  + zcloud * zrefo[jk];
        zrefd[jk] = zclear * zrefdc[jk] + zcloud * zrefdo[jk];
        ztra[jk]  = zclear * ztrac[jk]  + zcloud * ztrao[jk];
        ztrad[jk] = zclear * ztradc[jk] + zcloud * ztrado[jk];

        // Direct beam transmittance
        double ze1 = ztauc[jk] / prmu0;
        double zdbtmc;
        if (ze1 <= od_lo) {
          zdbtmc = 1.0 - ze1 + 0.5 * ze1 * ze1;
        } else {
          double tblind = ze1 / (bpade + ze1);
          int itind = static_cast<int>(tblint * tblind + 0.5);
          zdbtmc = tables.exp_tbl[itind];
        }
        zdbtc[jk] = zdbtmc;
        ztdbtc[jk + 1] = zdbtc[jk] * ztdbtc[jk];

        ze1 = ztauo[jk] / prmu0;
        double zdbtmo;
        if (ze1 <= od_lo) {
          zdbtmo = 1.0 - ze1 + 0.5 * ze1 * ze1;
        } else {
          double tblind = ze1 / (bpade + ze1);
          int itind = static_cast<int>(tblint * tblind + 0.5);
          zdbtmo = tables.exp_tbl[itind];
        }
        zdbt[jk] = zclear * zdbtmc + zcloud * zdbtmo;
        ztdbt[jk + 1] = zdbt[jk] * ztdbt[jk];
      }

      // Vertical quadrature for clear-sky fluxes
      vrtqdr_sw(klev, zrefc, zrefdc, ztrac, ztradc, zdbtc, ztdbtc, zrdndc, zrupc, zrupdc, zcd, zcu);

      // Vertical quadrature for total-sky fluxes
      vrtqdr_sw(klev, zref, zrefd, ztra, ztrad, zdbt, ztdbt, zrdnd, zrup, zrupd, zfd, zfu);

      // Accumulate fluxes (bottom-to-top output indexing)
      for (int jk = 0; jk <= klev; ++jk) {
        int ikl = klev - jk;
        pbbfu[ikl]    += zincflx * zfu[jk];
        pbbfd[ikl]    += zincflx * zfd[jk];
        pbbcu[ikl]    += zincflx * zcu[jk];
        pbbcd[ikl]    += zincflx * zcd[jk];
        pbbfddir[ikl] += zincflx * ztdbt_nodel[jk];
        pbbcddir[ikl] += zincflx * ztdbtc_nodel[jk];
      }

      // UV/visible direct/diffuse fluxes (bands 25-28, ibm 9-12 0-based)
      if (ibm >= 9 && ibm <= 12) {
        for (int jk = 0; jk <= klev; ++jk) {
          int ikl = klev - jk;
          puvcd[ikl]    += zincflx * zcd[jk];
          puvfd[ikl]    += zincflx * zfd[jk];
          puvcddir[ikl] += zincflx * ztdbtc_nodel[jk];
          puvfddir[ikl] += zincflx * ztdbt_nodel[jk];
        }
      } else if (ibm == 13 || ibm <= 8) { // Near-IR direct/diffuse fluxes
        for (int jk = 0; jk <= klev; ++jk) {
          int ikl = klev - jk;
          pnicd[ikl]    += zincflx * zcd[jk];
          pnifd[ikl]    += zincflx * zfd[jk];
          pnicddir[ikl] += zincflx * ztdbtc_nodel[jk];
          pnifddir[ikl] += zincflx * ztdbt_nodel[jk];
        }
      }

      iw++;
    }
  }
}

}  // namespace cppwrf::phys
