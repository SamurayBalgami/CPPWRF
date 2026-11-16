// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "phys/ra_rrtmg_common/rrtmg_real.hpp"
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

CPPWRF_HD inline void rrtmg_lw_secdiff(rr pwvcm, SD secdiff) {
  const rr a0[16] = {rr(1.66), rr(1.55), rr(1.58), rr(1.66), rr(1.54), rr(1.454), rr(1.89), rr(1.33),
                         rr(1.668), rr(1.66), rr(1.66), rr(1.66), rr(1.66), rr(1.66), rr(1.66), rr(1.66)};
  const rr a1[16] = {rr(0.00), rr(0.25), rr(0.22), rr(0.00), rr(0.13), rr(0.446), -rr(0.10), rr(0.40),
                         -rr(0.006), rr(0.00), rr(0.00), rr(0.00), rr(0.00), rr(0.00), rr(0.00), rr(0.00)};
  const rr a2[16] = {rr(0.00), -rr(12.0), -rr(11.7), rr(0.00), -rr(0.72), -rr(0.243), rr(0.19), -rr(0.062),
                         rr(0.414), rr(0.00), rr(0.00), rr(0.00), rr(0.00), rr(0.00), rr(0.00), rr(0.00)};
  for (int ibnd = 0; ibnd < 16; ++ibnd) {
    if (ibnd == 0 || ibnd == 3 || ibnd >= 9) {
      secdiff[ibnd] = rr(1.66);
    } else {
      rr s = a0[ibnd] + a1[ibnd] * libm::exp(a2[ibnd] * pwvcm);
      if (s > rr(1.80)) s = rr(1.80);
      if (s < rr(1.50)) s = rr(1.50);
      secdiff[ibnd] = s;
    }
  }
}

// icldlyr for layers [lay0, lay1)
CPPWRF_FORCEINLINE void rrtmg_lw_icldlyr(CSD2 cldfmc, SI icldlyr, int lay0, int lay1) {
  for (int lay = lay0; lay < lay1; ++lay) {
    int c = 0;
    for (int ig = 0; ig < 140; ++ig)
      if (cldfmc[ig][lay] == rr(1.0)) c = 1;
    icldlyr[lay] = c;
  }
}

// One g-point igc (of band iband) of rtrnmc: its down and up sweeps. The radiances the sweeps add to the band sums
// go to `out` per level, in the order of the sweeps:
//   out.down(l, radld, radclrd, iclddn)  the bottom of layer l (rtrnmc's lev - 1), iclddn after layer l
//   out.up(lev, radlu, radclru, iclddn)  level lev = 0..nlayers, iclddn at the end of the down sweep
// CLEAR_COL: a column without cloudy layers (rtrnmc's iclddn stays 0; radclrd = radld, radclru = radlu).
// A2 (CSD2 or SD2): the g-point's arrays; each element of fracs, taut is read once, at its layer of the down sweep,
// cldfmc's at its layer of each sweep, so `out` may overwrite them after (the GPU batch's in-place store).
template <int MAX_LAY, bool CLEAR_COL, class A2, class Out>
CPPWRF_FORCEINLINE void rrtmg_lw_gpoint(int nlayers, int iband, int igc, CSD secdiff, rr reflect, rr pbnd,
                                        A2 cldfmc, A2 taucmc, CSD2 planklay, CSD2 planklev, A2 fracs, A2 taut,
                                        CSI icldlyr, const RrtmgLwTables& CPPWRF_RESTRICT tables, Out& out) {
  constexpr rr rec_6 = rr(0.166667);
  constexpr rr bpade = rr(1.0) / rr(0.278);
  constexpr rr tblint = rr(10000.0);
  const rr sd = secdiff[iband];

  if constexpr (CLEAR_COL) {
    rr atrans[MAX_LAY], bbugas[MAX_LAY];
    {
      rr radld = rr(0.0);
      rr plfrac0 = rr(0.0);
      for (int lev = nlayers; lev >= 1; --lev) {
        int lay = lev - 1;
        rr plfrac = fracs[lay][igc];
        rr blay = planklay[lay][iband];
        rr dplankup = planklev[lev][iband] - blay;
        rr dplankdn = planklev[lev - 1][iband] - blay;
        rr odepth = sd * taut[lay][igc];
        if (odepth < rr(0.0)) odepth = rr(0.0);
        rr bbd = rr(0.0);
        if (odepth <= rr(0.06)) {
          atrans[lay] = odepth - rr(0.5) * odepth * odepth;
          rr odepth_rec = rec_6 * odepth;
          bbd = plfrac * (blay + dplankdn * odepth_rec);
          bbugas[lay] = plfrac * (blay + dplankup * odepth_rec);
        } else {
          rr tblind = odepth / (bpade + odepth);
          int itr = static_cast<int>(tblint * tblind + rr(0.5));
          rr transc = tables.exp_tbl[itr];
          atrans[lay] = rr(1.0) - transc;
          rr tausfac = tables.tfn_tbl[itr];
          bbd = plfrac * (blay + tausfac * dplankdn);
          bbugas[lay] = plfrac * (blay + tausfac * dplankup);
        }
        radld = radld + (bbd - radld) * atrans[lay];
        if (lay == 0) plfrac0 = plfrac;
        out.down(lev - 1, radld, radld, 0);
      }

      rr rad0 = plfrac0 * pbnd;  // fracs[0][igc] * pbnd
      rr radlu = rad0 + reflect * radld;
      out.up(0, radlu, radlu, 0);

      for (int lev = 1; lev <= nlayers; ++lev) {
        int lay = lev - 1;
        radlu = radlu + (bbugas[lay] - radlu) * atrans[lay];
        out.up(lev, radlu, radlu, 0);
      }
    }
  } else {
    rr atrans[MAX_LAY], bbugas[MAX_LAY], bbutot_atot[MAX_LAY], efcl[MAX_LAY];
    {
      rr radld = rr(0.0);
      rr radclrd = rr(0.0);
      rr plfrac0 = rr(0.0);
      int iclddn = 0;

      for (int lev = nlayers; lev >= 1; --lev) {
        int lay = lev - 1;
        rr plfrac = fracs[lay][igc];
        rr blay = planklay[lay][iband];
        rr dplankup = planklev[lev][iband] - blay;
        rr dplankdn = planklev[lev - 1][iband] - blay;
        rr odepth = sd * taut[lay][igc];
        if (odepth < rr(0.0)) odepth = rr(0.0);

        rr bbd = rr(0.0);

        if (icldlyr[lay] == 1) {
          iclddn = 1;
          rr odcld = rr(0.0);
          rr efcl_lay = rr(0.0);
          const rr cf = cldfmc[igc][lay];
          if (cf == rr(1.0)) {
            odcld = secdiff[tables.ngb[igc] - 1] * taucmc[igc][lay];
            rr transcld = libm::exp(-odcld);
            rr abscld = rr(1.0) - transcld;
            efcl_lay = abscld * cf;
          }
          efcl[lay] = efcl_lay;
          rr odtot = odepth + odcld;
          rr atot_lay, bbutot_lay;
          if (odtot < rr(0.06)) {
            atrans[lay] = odepth - rr(0.5) * odepth * odepth;
            rr odepth_rec = rec_6 * odepth;
            rr gassrc = plfrac * (blay + dplankdn * odepth_rec) * atrans[lay];

            atot_lay = odtot - rr(0.5) * odtot * odtot;
            rr odtot_rec = rec_6 * odtot;
            rr bbdtot = plfrac * (blay + dplankdn * odtot_rec);
            bbd = plfrac * (blay + dplankdn * odepth_rec);

            radld = radld - radld * (atrans[lay] + efcl_lay * (rr(1.0) - atrans[lay])) +
                    gassrc + cf * (bbdtot * atot_lay - gassrc);

            bbugas[lay] = plfrac * (blay + dplankup * odepth_rec);
            bbutot_lay = plfrac * (blay + dplankup * odtot_rec);
          } else if (odepth <= rr(0.06)) {
            atrans[lay] = odepth - rr(0.5) * odepth * odepth;
            rr odepth_rec = rec_6 * odepth;
            rr gassrc = plfrac * (blay + dplankdn * odepth_rec) * atrans[lay];

            rr tblind = odtot / (bpade + odtot);
            int ittot = static_cast<int>(tblint * tblind + rr(0.5));
            rr tfactot = tables.tfn_tbl[ittot];
            rr bbdtot = plfrac * (blay + tfactot * dplankdn);
            bbd = plfrac * (blay + dplankdn * odepth_rec);
            atot_lay = rr(1.0) - tables.exp_tbl[ittot];

            radld = radld - radld * (atrans[lay] + efcl_lay * (rr(1.0) - atrans[lay])) +
                    gassrc + cf * (bbdtot * atot_lay - gassrc);

            bbugas[lay] = plfrac * (blay + dplankup * odepth_rec);
            bbutot_lay = plfrac * (blay + tfactot * dplankup);
          } else {
            rr tblind = odepth / (bpade + odepth);
            int itgas = static_cast<int>(tblint * tblind + rr(0.5));
            odepth = tables.tau_tbl[itgas];
            atrans[lay] = rr(1.0) - tables.exp_tbl[itgas];
            rr tfacgas = tables.tfn_tbl[itgas];
            rr gassrc = atrans[lay] * plfrac * (blay + tfacgas * dplankdn);

            odtot = odepth + odcld;
            tblind = odtot / (bpade + odtot);
            int ittot = static_cast<int>(tblint * tblind + rr(0.5));
            rr tfactot = tables.tfn_tbl[ittot];
            rr bbdtot = plfrac * (blay + tfactot * dplankdn);
            bbd = plfrac * (blay + tfacgas * dplankdn);
            atot_lay = rr(1.0) - tables.exp_tbl[ittot];

            radld = radld - radld * (atrans[lay] + efcl_lay * (rr(1.0) - atrans[lay])) +
                    gassrc + cf * (bbdtot * atot_lay - gassrc);

            bbugas[lay] = plfrac * (blay + tfacgas * dplankup);
            bbutot_lay = plfrac * (blay + tfactot * dplankup);
          }
          bbutot_atot[lay] = bbutot_lay * atot_lay;
        } else {
          if (odepth <= rr(0.06)) {
            atrans[lay] = odepth - rr(0.5) * odepth * odepth;
            rr odepth_rec = rec_6 * odepth;
            bbd = plfrac * (blay + dplankdn * odepth_rec);
            bbugas[lay] = plfrac * (blay + dplankup * odepth_rec);
          } else {
            rr tblind = odepth / (bpade + odepth);
            int itr = static_cast<int>(tblint * tblind + rr(0.5));
            rr transc = tables.exp_tbl[itr];
            atrans[lay] = rr(1.0) - transc;
            rr tausfac = tables.tfn_tbl[itr];
            bbd = plfrac * (blay + tausfac * dplankdn);
            bbugas[lay] = plfrac * (blay + tausfac * dplankup);
          }
          radld = radld + (bbd - radld) * atrans[lay];
        }

        if (iclddn == 1) radclrd = radclrd + (bbd - radclrd) * atrans[lay];
        else radclrd = radld;
        if (lay == 0) plfrac0 = plfrac;
        out.down(lev - 1, radld, radclrd, iclddn);
      }

      rr rad0 = plfrac0 * pbnd;  // fracs[0][igc] * pbnd
      rr radlu = rad0 + reflect * radld;
      rr radclru = rad0 + reflect * radclrd;
      out.up(0, radlu, radclru, iclddn);

      for (int lev = 1; lev <= nlayers; ++lev) {
        int lay = lev - 1;
        if (icldlyr[lay] == 1) {
          rr gassrc = bbugas[lay] * atrans[lay];
          radlu = radlu - radlu * (atrans[lay] + efcl[lay] * (rr(1.0) - atrans[lay])) +
                  gassrc + cldfmc[igc][lay] * (bbutot_atot[lay] - gassrc);
        } else {
          radlu = radlu + (bbugas[lay] - radlu) * atrans[lay];
        }

        if (iclddn == 1) radclru = radclru + (bbugas[lay] - radclru) * atrans[lay];
        else radclru = radlu;
        out.up(lev, radlu, radclru, iclddn);
      }
    }
  }
}

// rtrnmc's band sums of one column (urad, drad, clrurad, clrdrad of one band): what each g-point's sweeps add,
// and below the first cloud (iclddn = 0) clrdrad = drad, clrurad = urad (level 0: clrurad adds radclru)
struct LwBandAcc {
  SD urad, drad, clrurad, clrdrad;
  CPPWRF_FORCEINLINE void down(int l, rr radld, rr radclrd, int iclddn) const {
    drad[l] += radld;
    if (iclddn == 1) clrdrad[l] += radclrd;
    else clrdrad[l] = drad[l];
  }
  CPPWRF_FORCEINLINE void up(int lev, rr radlu, rr radclru, int iclddn) const {
    urad[lev] += radlu;
    if (lev == 0 || iclddn == 1) clrurad[lev] += radclru;
    else clrurad[lev] = urad[lev];
  }
};
// the clear column's sums: the clear-sky sums are copied at the end (rtrnmc's clrdrad = drad, clrurad = urad)
struct LwBandAccClear {
  SD urad, drad;
  CPPWRF_FORCEINLINE void down(int l, rr radld, rr, int) const { drad[l] += radld; }
  CPPWRF_FORCEINLINE void up(int lev, rr radlu, rr, int) const { urad[lev] += radlu; }
};

template <int MAX_LAY, bool CLEAR_COL>
CPPWRF_FORCEINLINE void rrtmg_lw_band_impl(int nlayers, int iband, CSD secdiff, CSD semiss, CSD2 cldfmc, CSD2 taucmc,
                                           CSD2 planklay, CSD2 planklev, CSD plankbnd, CSD2 fracs, CSD2 taut,
                                           CSI icldlyr, SD2 bu, SD2 bd, SD2 bcu, SD2 bcd,
                                           const RrtmgLwTables& CPPWRF_RESTRICT tables) {
  const SD urad = bu[iband], drad = bd[iband], clrurad = bcu[iband], clrdrad = bcd[iband];
  for (int lev = 0; lev <= nlayers; ++lev) {
    urad[lev] = rr(0.0); drad[lev] = rr(0.0);
    if constexpr (!CLEAR_COL) {
      clrurad[lev] = rr(0.0); clrdrad[lev] = rr(0.0);
    }
  }
  const rr reflect = rr(1.0) - semiss[iband];
  const rr pbnd = plankbnd[iband];
  const int igc0 = iband == 0 ? 0 : tables.ngs[iband - 1];
  const int igc1 = tables.ngs[iband];

  if constexpr (CLEAR_COL) {
    LwBandAccClear acc{urad, drad};
    for (int igc = igc0; igc < igc1; ++igc)
      rrtmg_lw_gpoint<MAX_LAY, true, CSD2>(nlayers, iband, igc, secdiff, reflect, pbnd, cldfmc, taucmc, planklay, planklev,
                                     fracs, taut, icldlyr, tables, acc);
    for (int lev = 0; lev <= nlayers; ++lev) {
      clrdrad[lev] = drad[lev];
      clrurad[lev] = urad[lev];
    }
  } else {
    LwBandAcc acc{urad, drad, clrurad, clrdrad};
    for (int igc = igc0; igc < igc1; ++igc)
      rrtmg_lw_gpoint<MAX_LAY, false, CSD2>(nlayers, iband, igc, secdiff, reflect, pbnd, cldfmc, taucmc, planklay,
                                      planklev, fracs, taut, icldlyr, tables, acc);
  }
}

// One band of rtrnmc's band loop (iband 0..15): its g-points in order, down and up sweeps, the
// band sums urad/drad/clrurad/clrdrad of every level (bu/bd/bcu/bcd [16][nlayers + 1]) with WRF's
// clear-sky assignments (clrdrad = drad, clrurad = urad below the first cloud).
template <int MAX_LAY = RRTMG_LW_SWEEP_MAX>
CPPWRF_FORCEINLINE void rrtmg_lw_band(int nlayers, int iband, CSD secdiff, CSD semiss, CSD2 cldfmc, CSD2 taucmc,
                                      CSD2 planklay, CSD2 planklev, CSD plankbnd, CSD2 fracs, CSD2 taut, CSI icldlyr,
                                      SD2 bu, SD2 bd, SD2 bcu, SD2 bcd, const RrtmgLwTables& CPPWRF_RESTRICT tables,
                                      bool clear_col = false) {
#if defined(__CUDA_ARCH__)
  const bool fast_clear = __all_sync(__activemask(), clear_col);
#else
  const bool fast_clear = clear_col;
#endif
  if (fast_clear)
    rrtmg_lw_band_impl<MAX_LAY, true>(nlayers, iband, secdiff, semiss, cldfmc, taucmc, planklay, planklev, plankbnd,
                                      fracs, taut, icldlyr, bu, bd, bcu, bcd, tables);
  else
    rrtmg_lw_band_impl<MAX_LAY, false>(nlayers, iband, secdiff, semiss, cldfmc, taucmc, planklay, planklev, plankbnd,
                                       fracs, taut, icldlyr, bu, bd, bcu, bcd, tables);
}

// Level lev (0..nlayers): rtrnmc's band sums in band order, scaled by fluxfac; fnet / fnetc.
CPPWRF_FORCEINLINE void rrtmg_lw_accum(int lev, CSD2 bu, CSD2 bd, CSD2 bcu, CSD2 bcd,
                                       const RrtmgLwTables& CPPWRF_RESTRICT tables, rr& totuflux, rr& totdflux,
                                       rr& totuclfl, rr& totdclfl, rr& fnet, rr& fnetc) {
  constexpr rr pi = rr(3.14159265358979323846);
  constexpr rr fluxfac = pi * rr(2.0e4);
  constexpr rr wtdiff = rr(0.5);
  rr tu = rr(0.0), td = rr(0.0), tuc = rr(0.0), tdc = rr(0.0);
  for (int iband = 0; iband < 16; ++iband) {
    const rr dw = tables.delwave[iband];
    rr uflux = bu[iband][lev] * wtdiff;
    rr dflux = bd[iband][lev] * wtdiff;
    tu += uflux * dw;
    td += dflux * dw;
    rr uclfl = bcu[iband][lev] * wtdiff;
    rr dclfl = bcd[iband][lev] * wtdiff;
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

CPPWRF_FORCEINLINE void rrtmg_lw_heating(int nlayers, CSD pz, CSD fnet, CSD fnetc, SD htr, SD htrc, rr heatfac) {
  for (int lev = 1; lev <= nlayers; ++lev) {
    int l = lev - 1;
    htr[l] = heatfac * (fnet[l] - fnet[lev]) / (pz[l] - pz[lev]);
    htrc[l] = heatfac * (fnetc[l] - fnetc[lev]) / (pz[l] - pz[lev]);
  }
  htr[nlayers] = rr(0.0);
  htrc[nlayers] = rr(0.0);
}

}  // namespace cppwrf::phys
