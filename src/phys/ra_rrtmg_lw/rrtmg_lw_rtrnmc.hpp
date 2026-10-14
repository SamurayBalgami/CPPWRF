#pragma once

#include "core/backend.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_tables.hpp"
#include <algorithm>
#include <cmath>

namespace cppwrf::phys {

CPPWRF_HD inline void rrtmg_lw_rtrnmc(
    int nlayers,
    int istart,
    int iend,
    int iout,
    const double* pz,               // [nlayers + 1] (interface pressures)
    const double* semiss,          // [16]
    int ncbands,
    const double* const* cldfmc,   // [140][nlayers]
    const double* const* taucmc,   // [140][nlayers]
    const double* const* planklay, // [nlayers][16]
    const double* const* planklev, // [nlayers + 1][16]
    const double* plankbnd,        // [16]
    double pwvcm,
    const double* const* fracs,    // [nlayers][140]
    const double* const* taut,     // [nlayers][140]
    double* totuflux,              // [nlayers + 1]
    double* totdflux,              // [nlayers + 1]
    double* fnet,                  // [nlayers + 1]
    double* htr,                   // [nlayers + 1]
    double* totuclfl,              // [nlayers + 1]
    double* totdclfl,              // [nlayers + 1]
    double* fnetc,                 // [nlayers + 1]
    double* htrc,                  // [nlayers + 1]
    double heatfac,
    const RrtmgLwTables& tables
) {
  constexpr double pi = 3.14159265358979323846;
  constexpr double fluxfac = pi * 2.0e4;
  constexpr double wtdiff = 0.5;
  constexpr double rec_6 = 0.166667;
  constexpr double bpade = 1.0 / 0.278;
  constexpr double tblint = 10000.0;

  // Diffusivity angle parameters
  const double a0[16] = {1.66, 1.55, 1.58, 1.66, 1.54, 1.454, 1.89, 1.33,
                         1.668, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66};
  const double a1[16] = {0.00, 0.25, 0.22, 0.00, 0.13, 0.446, -0.10, 0.40,
                         -0.006, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
  const double a2[16] = {0.00, -12.0, -11.7, 0.00, -0.72, -0.243, 0.19, -0.062,
                         0.414, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};

  double secdiff[16];
  for (int ibnd = 0; ibnd < 16; ++ibnd) {
    if (ibnd == 0 || ibnd == 3 || ibnd >= 9) {
      secdiff[ibnd] = 1.66;
    } else {
      secdiff[ibnd] = a0[ibnd] + a1[ibnd] * std::exp(a2[ibnd] * pwvcm);
      if (secdiff[ibnd] > 1.80) secdiff[ibnd] = 1.80;
      if (secdiff[ibnd] < 1.50) secdiff[ibnd] = 1.50;
    }
  }

  // Initialize total flux arrays
  for (int lev = 0; lev <= nlayers; ++lev) {
    totuflux[lev] = 0.0;
    totdflux[lev] = 0.0;
    totuclfl[lev] = 0.0;
    totdclfl[lev] = 0.0;
  }

  // Work arrays (1D vertical vectors)
  constexpr int MAX_K = 128;
  double urad[MAX_K] = {0.0};
  double drad[MAX_K] = {0.0};
  double clrurad[MAX_K] = {0.0};
  double clrdrad[MAX_K] = {0.0};
  double atrans[MAX_K];
  double atot[MAX_K];
  double bbugas[MAX_K];
  double bbutot[MAX_K];
  // Precalculate cloud properties across all layers and g-channels
  double odcld[MAX_K][140];
  double abscld[MAX_K][140];
  double efclfrac[MAX_K][140];
  int icldlyr[MAX_K] = {0};

  for (int lay = 0; lay < nlayers; ++lay) {
    icldlyr[lay] = 0;
    for (int ig = 0; ig < 140; ++ig) {
      if (cldfmc[ig][lay] == 1.0) {
        int ib = tables.ngb[ig] - 1;
        odcld[lay][ig] = secdiff[ib] * taucmc[ig][lay];
        double transcld = std::exp(-odcld[lay][ig]);
        abscld[lay][ig] = 1.0 - transcld;
        efclfrac[lay][ig] = abscld[lay][ig] * cldfmc[ig][lay];
        icldlyr[lay] = 1;
      } else {
        odcld[lay][ig] = 0.0;
        abscld[lay][ig] = 0.0;
        efclfrac[lay][ig] = 0.0;
      }
    }
  }

  int igc = 0;

  // Loop over frequency bands
  for (int iband = istart - 1; iband < iend; ++iband) {
    if (iout > 0 && iband >= 1) {
      igc = tables.ngs[iband - 1];
    }

    int igc_end = tables.ngs[iband];

    // Loop over g-channels
    while (igc < igc_end) {
      // Downward radiative transfer
      double radld = 0.0;
      double radclrd = 0.0;
      int iclddn = 0;

      for (int lev = nlayers; lev >= 1; --lev) {
        int lay = lev - 1;
        double plfrac = fracs[lay][igc];
        double blay = planklay[lay][iband];
        double dplankup = planklev[lev][iband] - blay;
        double dplankdn = planklev[lev - 1][iband] - blay;
        double odepth = secdiff[iband] * taut[lay][igc];
        if (odepth < 0.0) odepth = 0.0;

        double bbd = 0.0;

        if (icldlyr[lay] == 1) {
          iclddn = 1;
          double odtot = odepth + odcld[lay][igc];
          if (odtot < 0.06) {
            atrans[lay] = odepth - 0.5 * odepth * odepth;
            double odepth_rec = rec_6 * odepth;
            double gassrc = plfrac * (blay + dplankdn * odepth_rec) * atrans[lay];

            atot[lay] = odtot - 0.5 * odtot * odtot;
            double odtot_rec = rec_6 * odtot;
            double bbdtot = plfrac * (blay + dplankdn * odtot_rec);
            bbd = plfrac * (blay + dplankdn * odepth_rec);

            radld = radld - radld * (atrans[lay] + efclfrac[lay][igc] * (1.0 - atrans[lay])) +
                    gassrc + cldfmc[igc][lay] * (bbdtot * atot[lay] - gassrc);
            drad[lev - 1] += radld;

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

            radld = radld - radld * (atrans[lay] + efclfrac[lay][igc] * (1.0 - atrans[lay])) +
                    gassrc + cldfmc[igc][lay] * (bbdtot * atot[lay] - gassrc);
            drad[lev - 1] += radld;

            bbugas[lay] = plfrac * (blay + dplankup * odepth_rec);
            bbutot[lay] = plfrac * (blay + tfactot * dplankup);
          } else {
            double tblind = odepth / (bpade + odepth);
            int itgas = static_cast<int>(tblint * tblind + 0.5);
            odepth = tables.tau_tbl[itgas];
            atrans[lay] = 1.0 - tables.exp_tbl[itgas];
            double tfacgas = tables.tfn_tbl[itgas];
            double gassrc = atrans[lay] * plfrac * (blay + tfacgas * dplankdn);

            odtot = odepth + odcld[lay][igc];
            tblind = odtot / (bpade + odtot);
            int ittot = static_cast<int>(tblint * tblind + 0.5);
            double tfactot = tables.tfn_tbl[ittot];
            double bbdtot = plfrac * (blay + tfactot * dplankdn);
            bbd = plfrac * (blay + tfacgas * dplankdn);
            atot[lay] = 1.0 - tables.exp_tbl[ittot];

            radld = radld - radld * (atrans[lay] + efclfrac[lay][igc] * (1.0 - atrans[lay])) +
                    gassrc + cldfmc[igc][lay] * (bbdtot * atot[lay] - gassrc);
            drad[lev - 1] += radld;

            bbugas[lay] = plfrac * (blay + tfacgas * dplankup);
            bbutot[lay] = plfrac * (blay + tfactot * dplankup);
          }
        } else {
          // Clear layer
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
          drad[lev - 1] += radld;
        }

        if (iclddn == 1) {
          radclrd = radclrd + (bbd - radclrd) * atrans[lay];
          clrdrad[lev - 1] += radclrd;
        } else {
          radclrd = radld;
          clrdrad[lev - 1] = drad[lev - 1];
        }
      }

      // 3. Surface reflection and upward radiance
      double rad0 = fracs[0][igc] * plankbnd[iband];
      double reflect = 1.0 - semiss[iband];
      double radlu = rad0 + reflect * radld;
      double radclru = rad0 + reflect * radclrd;

      urad[0] += radlu;
      clrurad[0] += radclru;

      // 4. Upward radiative transfer
      for (int lev = 1; lev <= nlayers; ++lev) {
        int lay = lev - 1;
        if (icldlyr[lay] == 1) {
          double gassrc = bbugas[lay] * atrans[lay];
          radlu = radlu - radlu * (atrans[lay] + efclfrac[lay][igc] * (1.0 - atrans[lay])) +
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

      // Increment g-channel
      igc++;
    }

    // Accumulate band fluxes
    double dw = tables.delwave[iband];
    for (int lev = nlayers; lev >= 0; --lev) {
      double uflux = urad[lev] * wtdiff;
      double dflux = drad[lev] * wtdiff;
      urad[lev] = 0.0;
      drad[lev] = 0.0;
      totuflux[lev] += uflux * dw;
      totdflux[lev] += dflux * dw;

      double uclfl = clrurad[lev] * wtdiff;
      double dclfl = clrdrad[lev] * wtdiff;
      clrurad[lev] = 0.0;
      clrdrad[lev] = 0.0;
      totuclfl[lev] += uclfl * dw;
      totdclfl[lev] += dclfl * dw;
    }
  }

  // 5. Final scale and heating rate calculations
  totuflux[0] *= fluxfac;
  totdflux[0] *= fluxfac;
  fnet[0] = totuflux[0] - totdflux[0];
  totuclfl[0] *= fluxfac;
  totdclfl[0] *= fluxfac;
  fnetc[0] = totuclfl[0] - totdclfl[0];

  for (int lev = 1; lev <= nlayers; ++lev) {
    totuflux[lev] *= fluxfac;
    totdflux[lev] *= fluxfac;
    fnet[lev] = totuflux[lev] - totdflux[lev];
    totuclfl[lev] *= fluxfac;
    totdclfl[lev] *= fluxfac;
    fnetc[lev] = totuclfl[lev] - totdclfl[lev];

    int l = lev - 1;
    htr[l] = heatfac * (fnet[l] - fnet[lev]) / (pz[l] - pz[lev]);
    htrc[l] = heatfac * (fnetc[l] - fnetc[lev]) / (pz[l] - pz[lev]);
  }

  htr[nlayers] = 0.0;
  htrc[nlayers] = 0.0;
}

}  // namespace cppwrf::phys
