#pragma once
// WRF rrtmg_lw rtrnmc (McICA radiative transfer), split into independent pieces so that a GPU can
// run every band and every level of a column in parallel with the same arithmetic as WRF:
//   rrtmg_lw_secdiff        per column: diffusivity angles (16 bands) from pwvcm
//   rrtmg_lw_icldlyr        per layer : cloudy layer flag (any McICA subcolumn cloudy)
//   rrtmg_lw_band           per band  : rtrnmc's g-point loop of one band (sweeps, band sums per level)
//   rrtmg_lw_accum          per level : the band sums in band order
//   rrtmg_lw_heating        per column: heating rates from the net fluxes
// Composed in order for one column they reproduce rtrnmc bit for bit.
#include "core/backend.hpp"
#include "phys/ra_rrtmg_common/strided.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_tables.hpp"
#include <algorithm>
#include <cmath>

namespace cppwrf::phys {

constexpr int RRTMG_LW_SWEEP_MAX = 128;  // layers of one column (stack arrays of rrtmg_lw_sweep)

CPPWRF_HD inline void rrtmg_lw_secdiff(double pwvcm, SD secdiff) {
  const double a0[16] = {1.66, 1.55, 1.58, 1.66, 1.54, 1.454, 1.89, 1.33,
                         1.668, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66};
  const double a1[16] = {0.00, 0.25, 0.22, 0.00, 0.13, 0.446, -0.10, 0.40,
                         -0.006, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
  const double a2[16] = {0.00, -12.0, -11.7, 0.00, -0.72, -0.243, 0.19, -0.062,
                         0.414, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
  for (int ibnd = 0; ibnd < 16; ++ibnd) {
    if (ibnd == 0 || ibnd == 3 || ibnd >= 9) {
      secdiff[ibnd] = 1.66;
    } else {
      double s = a0[ibnd] + a1[ibnd] * std::exp(a2[ibnd] * pwvcm);
      if (s > 1.80) s = 1.80;
      if (s < 1.50) s = 1.50;
      secdiff[ibnd] = s;
    }
  }
}

// icldlyr for layers [lay0, lay1)
CPPWRF_HD inline void rrtmg_lw_icldlyr(SD2 cldfmc, SI icldlyr, int lay0, int lay1) {
  for (int lay = lay0; lay < lay1; ++lay) {
    int c = 0;
    for (int ig = 0; ig < 140; ++ig)
      if (cldfmc[ig][lay] == 1.0) c = 1;
    icldlyr[lay] = c;
  }
}

// One band of rtrnmc's band loop (iband 0..15): its g-points in order, down and up sweeps, the
// band sums urad/drad/clrurad/clrdrad of every level (bu/bd/bcu/bcd [16][nlayers + 1]) with WRF's
// clear-sky assignments (clrdrad = drad, clrurad = urad below the first cloud).
CPPWRF_HD inline void rrtmg_lw_band(int nlayers, int iband, SD secdiff, SD semiss, SD2 cldfmc, SD2 taucmc,
                                    SD2 planklay, SD2 planklev, SD plankbnd, SD2 fracs, SD2 taut, SI icldlyr,
                                    SD2 bu, SD2 bd, SD2 bcu, SD2 bcd, const RrtmgLwTables& tables) {
  constexpr double rec_6 = 0.166667;
  constexpr double bpade = 1.0 / 0.278;
  constexpr double tblint = 10000.0;
  const SD urad = bu[iband], drad = bd[iband], clrurad = bcu[iband], clrdrad = bcd[iband];
  for (int lev = 0; lev <= nlayers; ++lev) {
    urad[lev] = 0.0; drad[lev] = 0.0; clrurad[lev] = 0.0; clrdrad[lev] = 0.0;
  }
  const double sd = secdiff[iband];
  double atrans[RRTMG_LW_SWEEP_MAX], atot[RRTMG_LW_SWEEP_MAX], bbugas[RRTMG_LW_SWEEP_MAX],
      bbutot[RRTMG_LW_SWEEP_MAX], efcl[RRTMG_LW_SWEEP_MAX];
  const int igc0 = iband == 0 ? 0 : tables.ngs[iband - 1];
  for (int igc = igc0; igc < tables.ngs[iband]; ++igc) {
    double radld = 0.0;
    double radclrd = 0.0;
    int iclddn = 0;

    for (int lev = nlayers; lev >= 1; --lev) {
      int lay = lev - 1;
      // odcld, abscld, efclfrac of rtrnmc for (lay, igc)
      double odcld = 0.0;
      efcl[lay] = 0.0;
      const double cf = cldfmc[igc][lay];
      if (cf == 1.0) {
        odcld = secdiff[tables.ngb[igc] - 1] * taucmc[igc][lay];
        double transcld = std::exp(-odcld);
        double abscld = 1.0 - transcld;
        efcl[lay] = abscld * cf;
      }
      double plfrac = fracs[lay][igc];
      double blay = planklay[lay][iband];
      double dplankup = planklev[lev][iband] - blay;
      double dplankdn = planklev[lev - 1][iband] - blay;
      double odepth = sd * taut[lay][igc];
      if (odepth < 0.0) odepth = 0.0;

      double bbd = 0.0;

      if (icldlyr[lay] == 1) {
        iclddn = 1;
        double odtot = odepth + odcld;
        if (odtot < 0.06) {
          atrans[lay] = odepth - 0.5 * odepth * odepth;
          double odepth_rec = rec_6 * odepth;
          double gassrc = plfrac * (blay + dplankdn * odepth_rec) * atrans[lay];

          atot[lay] = odtot - 0.5 * odtot * odtot;
          double odtot_rec = rec_6 * odtot;
          double bbdtot = plfrac * (blay + dplankdn * odtot_rec);
          bbd = plfrac * (blay + dplankdn * odepth_rec);

          radld = radld - radld * (atrans[lay] + efcl[lay] * (1.0 - atrans[lay])) +
                  gassrc + cf * (bbdtot * atot[lay] - gassrc);

          bbugas[lay] = plfrac * (blay + dplankup * odepth_rec);
          bbutot[lay] = plfrac * (blay + dplankup * odtot_rec);
        } else if (odepth <= 0.06) {
          atrans[lay] = odepth - 0.5 * odepth * odepth;
          double odepth_rec = rec_6 * odepth;
          double gassrc = plfrac * (blay + dplankdn * odepth_rec) * atrans[lay];

          double tblind = odtot / (bpade + odtot);
          int ittot = static_cast<int>(tblint * tblind + 0.5);
          double tfactot = tables.tfn_tbl[ittot];
          double bbdtot = plfrac * (blay + tfactot * dplankdn);
          bbd = plfrac * (blay + dplankdn * odepth_rec);
          atot[lay] = 1.0 - tables.exp_tbl[ittot];

          radld = radld - radld * (atrans[lay] + efcl[lay] * (1.0 - atrans[lay])) +
                  gassrc + cf * (bbdtot * atot[lay] - gassrc);

          bbugas[lay] = plfrac * (blay + dplankup * odepth_rec);
          bbutot[lay] = plfrac * (blay + tfactot * dplankup);
        } else {
          double tblind = odepth / (bpade + odepth);
          int itgas = static_cast<int>(tblint * tblind + 0.5);
          odepth = tables.tau_tbl[itgas];
          atrans[lay] = 1.0 - tables.exp_tbl[itgas];
          double tfacgas = tables.tfn_tbl[itgas];
          double gassrc = atrans[lay] * plfrac * (blay + tfacgas * dplankdn);

          odtot = odepth + odcld;
          tblind = odtot / (bpade + odtot);
          int ittot = static_cast<int>(tblint * tblind + 0.5);
          double tfactot = tables.tfn_tbl[ittot];
          double bbdtot = plfrac * (blay + tfactot * dplankdn);
          bbd = plfrac * (blay + tfacgas * dplankdn);
          atot[lay] = 1.0 - tables.exp_tbl[ittot];

          radld = radld - radld * (atrans[lay] + efcl[lay] * (1.0 - atrans[lay])) +
                  gassrc + cf * (bbdtot * atot[lay] - gassrc);

          bbugas[lay] = plfrac * (blay + tfacgas * dplankup);
          bbutot[lay] = plfrac * (blay + tfactot * dplankup);
        }
      } else {
        if (odepth <= 0.06) {
          atrans[lay] = odepth - 0.5 * odepth * odepth;
          double odepth_rec = rec_6 * odepth;
          bbd = plfrac * (blay + dplankdn * odepth_rec);
          bbugas[lay] = plfrac * (blay + dplankup * odepth_rec);
        } else {
          double tblind = odepth / (bpade + odepth);
          int itr = static_cast<int>(tblint * tblind + 0.5);
          double transc = tables.exp_tbl[itr];
          atrans[lay] = 1.0 - transc;
          double tausfac = tables.tfn_tbl[itr];
          bbd = plfrac * (blay + tausfac * dplankdn);
          bbugas[lay] = plfrac * (blay + tausfac * dplankup);
        }
        radld = radld + (bbd - radld) * atrans[lay];
      }
      drad[lev - 1] += radld;

      if (iclddn == 1) {
        radclrd = radclrd + (bbd - radclrd) * atrans[lay];
        clrdrad[lev - 1] += radclrd;
      } else {
        radclrd = radld;
        clrdrad[lev - 1] = drad[lev - 1];
      }
    }

    double rad0 = fracs[0][igc] * plankbnd[iband];
    double reflect = 1.0 - semiss[iband];
    double radlu = rad0 + reflect * radld;
    double radclru = rad0 + reflect * radclrd;

    urad[0] += radlu;
    clrurad[0] += radclru;

    for (int lev = 1; lev <= nlayers; ++lev) {
      int lay = lev - 1;
      if (icldlyr[lay] == 1) {
        double gassrc = bbugas[lay] * atrans[lay];
        radlu = radlu - radlu * (atrans[lay] + efcl[lay] * (1.0 - atrans[lay])) +
                gassrc + cldfmc[igc][lay] * (bbutot[lay] * atot[lay] - gassrc);
        urad[lev] += radlu;
      } else {
        radlu = radlu + (bbugas[lay] - radlu) * atrans[lay];
        urad[lev] += radlu;
      }

      if (iclddn == 1) {
        radclru = radclru + (bbugas[lay] - radclru) * atrans[lay];
        clrurad[lev] += radclru;
      } else {
        radclru = radlu;
        clrurad[lev] = urad[lev];
      }
    }
  }
}

// Level lev (0..nlayers): rtrnmc's band sums in band order, scaled by fluxfac; fnet / fnetc.
CPPWRF_HD inline void rrtmg_lw_accum(int lev, SD2 bu, SD2 bd, SD2 bcu, SD2 bcd, const RrtmgLwTables& tables,
                                     double& totuflux, double& totdflux, double& totuclfl, double& totdclfl,
                                     double& fnet, double& fnetc) {
  constexpr double pi = 3.14159265358979323846;
  constexpr double fluxfac = pi * 2.0e4;
  constexpr double wtdiff = 0.5;
  double tu = 0.0, td = 0.0, tuc = 0.0, tdc = 0.0;
  for (int iband = 0; iband < 16; ++iband) {
    const double dw = tables.delwave[iband];
    double uflux = bu[iband][lev] * wtdiff;
    double dflux = bd[iband][lev] * wtdiff;
    tu += uflux * dw;
    td += dflux * dw;
    double uclfl = bcu[iband][lev] * wtdiff;
    double dclfl = bcd[iband][lev] * wtdiff;
    tuc += uclfl * dw;
    tdc += dclfl * dw;
  }
  tu *= fluxfac;
  td *= fluxfac;
  fnet = tu - td;
  tuc *= fluxfac;
  tdc *= fluxfac;
  fnetc = tuc - tdc;
  totuflux = tu; totdflux = td; totuclfl = tuc; totdclfl = tdc;
}

CPPWRF_HD inline void rrtmg_lw_heating(int nlayers, SD pz, SD fnet, SD fnetc, SD htr, SD htrc, double heatfac) {
  for (int lev = 1; lev <= nlayers; ++lev) {
    int l = lev - 1;
    htr[l] = heatfac * (fnet[l] - fnet[lev]) / (pz[l] - pz[lev]);
    htrc[l] = heatfac * (fnetc[l] - fnetc[lev]) / (pz[l] - pz[lev]);
  }
  htr[nlayers] = 0.0;
  htrc[nlayers] = 0.0;
}

}  // namespace cppwrf::phys
