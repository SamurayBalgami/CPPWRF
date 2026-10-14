#include "phys/ra_rrtmg_sw/rrtmg_sw.hpp"
#include "phys/ra_rrtmg_common/rrtmg_common.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace cppwrf::phys {

CPPWRF_HD // WRF rrtmg_sw (cldprmc_sw, inatm_sw, setcoef_sw, spcvmc_sw, fluxes, heating rates) for one column whose
// McICA subcolumn arrays are in ws. Level arrays: nlayers + 1 (0 = surface); heating rates K/day.
void rrtmg_sw_core(int nlayers, double heatfac, const double* play, const double* plev, const double* tlay,
                   const double* tlev, const double* h2ovmr, const double* o3vmr, const double* co2vmr,
                   const double* ch4vmr, const double* n2ovmr, const double* o2vmr, double asdir, double asdif,
                   double aldir, double aldif, double coszen, double adjes, double scon, const double* reicmcl,
                   const double* relqmcl, const double* resnmcl, double* swuflx, double* swdflx, double* swhr,
                   double* swuflxc, double* swdflxc, double* swhrc, double* swdkdir, double* swdkdif,
                   double* sibvisdir, double* sibvisdif, double* sibnirdir, double* sibnirdif, RrtmgSwWorkspace& ws,
                   const RrtmgSwTables& tables, int inflag, int iceflag, int liqflag) {
  cldprmc_sw(nlayers, inflag, iceflag, liqflag, ws.cldfmc, ws.ciwpmc, ws.clwpmc, ws.cswpmc,
             reicmcl, relqmcl, resnmcl, ws.fsfcmc, ws.taormc, ws.taucmc, ws.ssacmc, ws.asmcmc, tables);
  // 7. Atmospheric preparation (inatm)
  constexpr double amd = 28.9660;
  constexpr double amw = 18.0160;
  constexpr double avogad = 6.02214199e23;
  constexpr double grav = 9.8066;

  double pavel[RRTMG_SW_MAX_LAYERS];
  double tavel[RRTMG_SW_MAX_LAYERS];
  double pz[RRTMG_SW_MAX_LAYERS + 1];
  double tz[RRTMG_SW_MAX_LAYERS + 1];
  double coldry[RRTMG_SW_MAX_LAYERS];

  pz[0] = plev[0];
  tz[0] = tlev[0];
  ws.init_pointers();

  double pdp[RRTMG_SW_MAX_LAYERS];
  for (int l = 0; l < nlayers; ++l) {
    pavel[l] = play[l];
    tavel[l] = tlay[l];
    pz[l + 1] = plev[l + 1];
    tz[l + 1] = tlev[l + 1];
    pdp[l] = pz[l] - pz[l + 1];

    ws.wkl[0][l] = h2ovmr[l];
    ws.wkl[1][l] = co2vmr[l];
    ws.wkl[2][l] = o3vmr[l];
    ws.wkl[3][l] = n2ovmr[l];
    ws.wkl[4][l] = 0.0;
    ws.wkl[5][l] = ch4vmr[l];
    ws.wkl[6][l] = o2vmr[l];

    double amm = (1.0 - ws.wkl[0][l]) * amd + ws.wkl[0][l] * amw;
    coldry[l] = (pz[l] - pz[l + 1]) * 1.0e3 * avogad / (1.0e2 * grav * amm * (1.0 + ws.wkl[0][l]));
    for (int imol = 0; imol < 7; ++imol) {
      ws.wkl[imol][l] *= coldry[l];
    }
  }

  // 8. Molecular absorption coefficients interpolation (setcoef)
  int laytrop, layswtch, laylow;
  int jp[RRTMG_SW_MAX_LAYERS];
  int jt[RRTMG_SW_MAX_LAYERS];
  int jt1[RRTMG_SW_MAX_LAYERS];
  double colh2o[RRTMG_SW_MAX_LAYERS];
  double colco2[RRTMG_SW_MAX_LAYERS];
  double colo3[RRTMG_SW_MAX_LAYERS];
  double coln2o[RRTMG_SW_MAX_LAYERS];
  double colch4[RRTMG_SW_MAX_LAYERS];
  double colo2[RRTMG_SW_MAX_LAYERS];
  double colmol[RRTMG_SW_MAX_LAYERS];
  double co2mult[RRTMG_SW_MAX_LAYERS];
  int indself[RRTMG_SW_MAX_LAYERS];
  int indfor[RRTMG_SW_MAX_LAYERS];
  double selffac[RRTMG_SW_MAX_LAYERS];
  double selffrac[RRTMG_SW_MAX_LAYERS];
  double forfac[RRTMG_SW_MAX_LAYERS];
  double forfrac[RRTMG_SW_MAX_LAYERS];
  double fac00[RRTMG_SW_MAX_LAYERS];
  double fac01[RRTMG_SW_MAX_LAYERS];
  double fac10[RRTMG_SW_MAX_LAYERS];
  double fac11[RRTMG_SW_MAX_LAYERS];

  setcoef_sw(
      nlayers, pavel, tavel, coldry, ws.wkl_ptrs,
      laytrop, layswtch, laylow, jp, jt, jt1,
      colh2o, colco2, colo3, coln2o, colch4, colo2, colmol, co2mult,
      indself, indfor, selffac, selffrac, forfac, forfrac,
      fac00, fac01, fac10, fac11, tables
  );

  // 9. Surface albedo per spectral band: WRF bands 1-9 and 14 near-IR (al*), 10-13 UV/vis (as*)
  double albdir[14];
  double albdif[14];
  for (int ib = 0; ib < 9; ++ib) { albdir[ib] = aldir; albdif[ib] = aldif; }
  albdir[13] = aldir;
  albdif[13] = aldif;
  for (int ib = 9; ib < 13; ++ib) { albdir[ib] = asdir; albdif[ib] = asdif; }

  // 10. Solar flux scaling (inatm_sw: adjflux = adjflx * solvar, dyofyr = 0 -> adjflx = adjes)
  double adjflux[30];
  constexpr double rrsw_scon = 1.36822e+03;
  for (int ib = 0; ib < 30; ++ib) adjflux[ib] = adjes * (scon / rrsw_scon);
  double prmu0 = coszen;
  if (prmu0 <= 1.e-10) prmu0 = 1.e-10;  // zepzen

  // Zero aerosols
  double taua[RRTMG_SW_MAX_LAYERS * 14] = {0.0};
  double asma[RRTMG_SW_MAX_LAYERS * 14] = {0.0};
  double ssaa[RRTMG_SW_MAX_LAYERS * 14];
  for (int i = 0; i < nlayers * 14; ++i) ssaa[i] = 1.0;

  // 11. Call 2-stream spectral radiative transfer solver (spcvmc_sw)
  double pbbfd[RRTMG_SW_MAX_LAYERS + 1];
  double pbbfu[RRTMG_SW_MAX_LAYERS + 1];
  double pbbcd[RRTMG_SW_MAX_LAYERS + 1];
  double pbbcu[RRTMG_SW_MAX_LAYERS + 1];
  double puvfd[RRTMG_SW_MAX_LAYERS + 1];
  double puvcd[RRTMG_SW_MAX_LAYERS + 1];
  double pnifd[RRTMG_SW_MAX_LAYERS + 1];
  double pnicd[RRTMG_SW_MAX_LAYERS + 1];
  double pbbfddir[RRTMG_SW_MAX_LAYERS + 1];
  double pbbcddir[RRTMG_SW_MAX_LAYERS + 1];
  double puvfddir[RRTMG_SW_MAX_LAYERS + 1];
  double puvcddir[RRTMG_SW_MAX_LAYERS + 1];
  double pnifddir[RRTMG_SW_MAX_LAYERS + 1];
  double pnicddir[RRTMG_SW_MAX_LAYERS + 1];

  spcvmc_sw(
      nlayers,
      /*icpr=*/1,
      albdif, albdir,
      ws.cldfmc, ws.taucmc, ws.asmcmc, ws.ssacmc, ws.taormc,
      taua, asma, ssaa,
      prmu0, adjflux,
      laytrop, layswtch, laylow, jp, jt, jt1,
      colh2o, colco2, colch4, colo2, colo3, colmol,
      fac00, fac01, fac10, fac11,
      selffac, selffrac, indself, forfac, forfrac, indfor,
      pbbfd, pbbfu, pbbcd, pbbcu,
      puvfd, puvcd, pnifd, pnicd,
      pbbfddir, pbbcddir, puvfddir, puvcddir, pnifddir, pnicddir,
      tables
  );

  for (int i = 0; i <= nlayers; ++i) {
    swuflxc[i] = pbbcu[i];
    swdflxc[i] = pbbcd[i];
    swuflx[i] = pbbfu[i];
    swdflx[i] = pbbfd[i];
    swdkdir[i] = pbbfddir[i];
    swdkdif[i] = swdflx[i] - pbbfddir[i];
    sibvisdir[i] = puvfddir[i];
    sibvisdif[i] = puvfd[i] - puvfddir[i];
    sibnirdir[i] = pnifddir[i];
    sibnirdif[i] = pnifd[i] - pnifddir[i];
  }
  for (int i = 0; i < nlayers; ++i) {
    const double zdpgcp = heatfac / pdp[i];
    swhrc[i] = ((swdflxc[i + 1] - swuflxc[i + 1]) - (swdflxc[i] - swuflxc[i])) * zdpgcp;
    swhr[i] = ((swdflx[i + 1] - swuflx[i + 1]) - (swdflx[i] - swuflx[i])) * zdpgcp;
  }
  swhrc[nlayers - 1] = 0.;
  swhr[nlayers - 1] = 0.;
}

void rrtmg_sw_column(
    int num_levels,
    const RrtmgSwConfig& config,
    const RrtmgSwTables& tables,
    double coszen,
    double albedo,
    double xland,
    double xice,
    double snow,
    double xlat,
    const double* p8w,
    const double* t8w,
    const double* p3d,
    const double* pi3d,
    const double* dz8w,
    const double* t3d,
    const double* rho3d,
    const double* qv3d,
    const double* qc3d,
    const double* qr3d,
    const double* qi3d,
    const double* qs3d,
    const double* qg3d,
    const double* cldfra3d,
    const double* o33d,
    double& gsw,
    double& swcf,
    double& swupt,
    double& swuptc,
    double& swdnt,
    double& swdntc,
    double& swupb,
    double& swupbc,
    double& swdnb,
    double& swdnbc,
    double& swvisdir,
    double& swvisdif,
    double& swnirdir,
    double& swnirdif,
    double* rthratensw,
    double* rthratenswc,
    RrtmgSwWorkspace& ws,
    const double* re_cloud,
    const double* re_ice,
    const double* re_snow
) {
  // If sun is below horizon, deactivate shortwave
  if (coszen <= 0.0) {
    gsw = 0.0;
    swcf = 0.0;
    swupt = 0.0;
    swuptc = 0.0;
    swdnt = 0.0;
    swdntc = 0.0;
    swupb = 0.0;
    swupbc = 0.0;
    swdnb = 0.0;
    swdnbc = 0.0;
    swvisdir = 0.0;
    swvisdif = 0.0;
    swnirdir = 0.0;
    swnirdif = 0.0;
    for (int k = 0; k < num_levels; ++k) {
      rthratensw[k] = 0.0;
      rthratenswc[k] = 0.0;
    }
    return;
  }

  int nlayers = num_levels + 1; // SW adds 1 buffer layer to TOA (1.0e-5 mb)

  // 1. Atmosphere base arrays
  double play[RRTMG_SW_MAX_LAYERS];
  double plev[RRTMG_SW_MAX_LAYERS + 1];
  double tlay[RRTMG_SW_MAX_LAYERS];
  double tlev[RRTMG_SW_MAX_LAYERS + 1];
  double pdel[RRTMG_SW_MAX_LAYERS];
  double hgt[RRTMG_SW_MAX_LAYERS];

  double h2ovmr[RRTMG_SW_MAX_LAYERS];
  double o3vmr[RRTMG_SW_MAX_LAYERS];
  double co2vmr[RRTMG_SW_MAX_LAYERS];
  double ch4vmr[RRTMG_SW_MAX_LAYERS];
  double n2ovmr[RRTMG_SW_MAX_LAYERS];
  double o2vmr[RRTMG_SW_MAX_LAYERS];

  double cldfrac[RRTMG_SW_MAX_LAYERS];
  double clwpth[RRTMG_SW_MAX_LAYERS];
  double ciwpth[RRTMG_SW_MAX_LAYERS];
  double cswpth[RRTMG_SW_MAX_LAYERS];
  double rel[RRTMG_SW_MAX_LAYERS];
  double rei[RRTMG_SW_MAX_LAYERS];
  double res[RRTMG_SW_MAX_LAYERS];

  double relqmcl[RRTMG_SW_MAX_LAYERS];
  double reicmcl[RRTMG_SW_MAX_LAYERS];
  double resnmcl[RRTMG_SW_MAX_LAYERS];

  constexpr double amdw = 1.607793;
  constexpr double amdo = 0.603461;
  constexpr double o2_const = 0.209488;
  const double ch4_const = config.ch4;  // ghg_input = 1: read_CAMgases, else the WRF defaults
  const double n2o_const = config.n2o;
  const double co2_const = config.co2;

  plev[0] = p8w[0] / 100.0; // Pa -> mb
  tlev[0] = t8w[0];

  double dzsum = 0.0;
  for (int k = 0; k < num_levels; ++k) {
    play[k] = p3d[k] / 100.0;
    plev[k + 1] = p8w[k + 1] / 100.0;
    tlay[k] = t3d[k];
    tlev[k + 1] = t8w[k + 1];
    pdel[k] = plev[k] - plev[k + 1];

    double dz = dz8w[k];
    hgt[k] = dzsum + 0.5 * dz;
    dzsum += dz;

    h2ovmr[k] = std::max(qv3d[k], 1.0e-12) * amdw;
    co2vmr[k] = co2_const;
    o2vmr[k] = o2_const;
    ch4vmr[k] = ch4_const;
    n2ovmr[k] = n2o_const;
  }

  // 2. Extra layer from model top to TOA (1.0e-5 mb)
  play[num_levels] = 0.5 * plev[num_levels];
  tlay[num_levels] = tlev[num_levels];
  plev[num_levels + 1] = 1.0e-5;
  tlev[num_levels + 1] = tlev[num_levels];
  pdel[num_levels] = plev[num_levels] - plev[num_levels + 1];
  hgt[num_levels] = dzsum + 0.5 * dz8w[num_levels - 1];

  h2ovmr[num_levels] = h2ovmr[num_levels - 1];
  co2vmr[num_levels] = co2vmr[num_levels - 1];
  o2vmr[num_levels] = o2vmr[num_levels - 1];
  ch4vmr[num_levels] = ch4vmr[num_levels - 1];
  n2ovmr[num_levels] = n2ovmr[num_levels - 1];

  // 3. Ozone profile calculation via inirad
  double o3mmr[RRTMG_SW_MAX_LAYERS + 1];
  rrtmg_inirad(o3mmr, plev, 0, num_levels);
  for (int k = 0; k < nlayers; ++k) {
    o3vmr[k] = o3mmr[k] * amdo;
    if (config.o3input == 2) {  // CAM ozone on model levels, the climatology shifted to it above
      if (k < num_levels) {
        o3vmr[k] = o33d[k];
      } else {
        o3vmr[k] = o33d[num_levels - 1] - o3mmr[num_levels - 1] * amdo + o3mmr[k] * amdo;
        if (o3vmr[k] <= 0.) o3vmr[k] = o3mmr[k] * amdo;
      }
    }
  }

  // 4. Cloud physical properties (RRTMG_SWRAD, icloud != 0; EM_CORE branches)
  int inflgsw = 2, iceflgsw = 3;
  const int liqflgsw = 1;
  for (int k = 0; k < num_levels; ++k) {
    const double cldf = (cldfra3d && config.icloud != 0) ? cldfra3d[k] : 0.0;
    double qc1 = std::max(0., qc3d[k]), qi1 = std::max(0., qi3d[k]), qs1 = std::max(0., qs3d[k]);
    double recloud1d = 5.0, reice1d = 10., resnow1d = 10.0;
    if (re_cloud) {
      recloud1d = std::max(2.5, re_cloud[k] * 1.E6);
      if (recloud1d <= 2.5 && cldf > 0. && (xland - 1.5) > 0.) recloud1d = 10.5;       // ocean
      else if (recloud1d <= 2.5 && cldf > 0. && (xland - 1.5) < 0.) recloud1d = 7.5;   // land
      inflgsw = 3;
    }
    if (re_ice) {
      inflgsw = 4; iceflgsw = 4;
      reice1d = std::max(5., re_ice[k] * 1.E6);
      if (reice1d <= 5. && cldf > 0.) {
        int idx_rei = int(t3d[k] - 179.);
        idx_rei = std::min(std::max(idx_rei, 1), 75);
        const double corr = t3d[k] - int(t3d[k]);
        reice1d = RETAB[idx_rei - 1] * (1. - corr) + RETAB[idx_rei] * corr;
        reice1d = std::max(reice1d, 5.0);
      }
    }
    if (re_snow) {
      inflgsw = 5; iceflgsw = 5;
      resnow1d = std::max(10., re_snow[k] * 1.E6);
    }
    if (!re_snow && re_ice && re_cloud) {  // P3: ice into the snow category for the optics
      inflgsw = 5; iceflgsw = 5;
      resnow1d = std::max(10., re_ice[k] * 1.E6);
      qs1 = qi3d[k];
      qi1 = 0.;
      reice1d = 10.;
    }
    cldfrac[k] = cldf;
    double gicewp = (qi1 + qs1) * pdel[k] * 100.0 / config.g * 1000.0; // g/m2
    const double gliqwp = qc1 * pdel[k] * 100.0 / config.g * 1000.0;   // g/m2
    ciwpth[k] = gicewp / std::max(0.01, cldf);
    clwpth[k] = gliqwp / std::max(0.01, cldf);
    if (iceflgsw >= 4) {  // ice path without the snow
      gicewp = qi1 * pdel[k] * 100.0 / config.g * 1000.0;
      ciwpth[k] = gicewp / std::max(0.01, cldf);
    }
    double csnowp = 0.0;
    if (iceflgsw == 5) {
      double snow_mass_factor = 0.99;
      if (resnow1d > 130.) {
        snow_mass_factor = std::min(snow_mass_factor, (130.0 / resnow1d) * (130.0 / resnow1d));
        resnow1d = 130.0;
      }
      const double gsnowp = qs1 * snow_mass_factor * pdel[k] * 100.0 / config.g * 1000.0;
      csnowp = gsnowp / std::max(0.01, cldf);
    }
    rel[k] = inflgsw >= 3 ? recloud1d : rrtmg_relcalc_single(tlay[k], xland, snow, xice);
    if (iceflgsw >= 4) rei[k] = reice1d;
    else rei[k] = std::min(140.0, rrtmg_reicalc_single(tlay[k]) * 1.0315);  // Fu ice size conversion
    if (inflgsw == 5) {
      cswpth[k] = csnowp;
      res[k] = resnow1d;
    } else {
      cswpth[k] = 0.0;
      res[k] = 10.0;
    }
  }

  cldfrac[num_levels] = 0.0;
  clwpth[num_levels] = 0.0;
  ciwpth[num_levels] = 0.0;
  cswpth[num_levels] = 0.0;
  rel[num_levels] = 10.0;
  rei[num_levels] = 10.0;
  res[num_levels] = 10.0;

  // 5. McICA stochastic subcolumn generator (KISSVEC PRNG, permuteseed = 1)
  KissRng rng;
  rng.init_from_pressure(play[0] * 100.0, play[1] * 100.0, play[2] * 100.0, play[3] * 100.0);
  rng.permute(1);

  double cldf[RRTMG_SW_MAX_LAYERS];
  for (int ilev = 0; ilev < nlayers; ++ilev) {
    cldf[ilev] = (cldfrac[ilev] < 1.0e-20) ? 0.0 : cldfrac[ilev];
  }

  double cdf[112][RRTMG_SW_MAX_LAYERS];
  for (int isubcol = 0; isubcol < 112; ++isubcol) {
    for (int ilev = 0; ilev < nlayers; ++ilev) {
      cdf[isubcol][ilev] = rng.step<double>();
    }
  }

  // Max/random overlap (cldovrlp == 2)
  for (int ilev = 1; ilev < nlayers; ++ilev) {
    for (int isubcol = 0; isubcol < 112; ++isubcol) {
      if (cdf[isubcol][ilev - 1] > 1.0 - cldf[ilev - 1]) {
        cdf[isubcol][ilev] = cdf[isubcol][ilev - 1];
      } else {
        cdf[isubcol][ilev] = cdf[isubcol][ilev] * (1.0 - cldf[ilev - 1]);
      }
    }
  }

  for (int ilev = 0; ilev < nlayers; ++ilev) {
    reicmcl[ilev] = rei[ilev];
    relqmcl[ilev] = rel[ilev];
    resnmcl[ilev] = res[ilev];

    for (int isubcol = 0; isubcol < 112; ++isubcol) {
      int idx = ilev * 112 + isubcol;
      bool is_cloudy = (cdf[isubcol][ilev] >= 1.0 - cldf[ilev]);
      if (is_cloudy) {
        ws.cldfmc[idx] = 1.0;
        ws.clwpmc[idx] = clwpth[ilev];
        ws.ciwpmc[idx] = ciwpth[ilev];
        ws.cswpmc[idx] = cswpth[ilev];
        ws.taucmc[idx] = 0.0;
        ws.ssacmc[idx] = 1.0;
        ws.asmcmc[idx] = 0.0;
        ws.fsfcmc[idx] = 0.0;
      } else {
        ws.cldfmc[idx] = 0.0;
        ws.clwpmc[idx] = 0.0;
        ws.ciwpmc[idx] = 0.0;
        ws.cswpmc[idx] = 0.0;
        ws.taucmc[idx] = 0.0;
        ws.ssacmc[idx] = 1.0;
        ws.asmcmc[idx] = 0.0;
        ws.fsfcmc[idx] = 0.0;
      }
    }
  }

  // 6-11. rrtmg_sw core (cldprmc_sw, inatm_sw, setcoef_sw, spcvmc_sw), WRF dyofyr = 0, adjes = 1
  double swuflx[RRTMG_SW_MAX_LAYERS + 1], swdflx[RRTMG_SW_MAX_LAYERS + 1], swhr[RRTMG_SW_MAX_LAYERS];
  double swuflxc[RRTMG_SW_MAX_LAYERS + 1], swdflxc[RRTMG_SW_MAX_LAYERS + 1], swhrc[RRTMG_SW_MAX_LAYERS];
  double swdkdir[RRTMG_SW_MAX_LAYERS + 1], swdkdif[RRTMG_SW_MAX_LAYERS + 1];
  double uvdir[RRTMG_SW_MAX_LAYERS + 1], uvdif[RRTMG_SW_MAX_LAYERS + 1], nidir[RRTMG_SW_MAX_LAYERS + 1],
      nidif[RRTMG_SW_MAX_LAYERS + 1];
  rrtmg_sw_core(nlayers, config.heatfac, play, plev, tlay, tlev, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,
                albedo, albedo, albedo, albedo, coszen, 1.0, config.solcon, reicmcl, relqmcl, resnmcl, swuflx, swdflx,
                swhr, swuflxc, swdflxc, swhrc, swdkdir, swdkdif, uvdir, uvdif, nidir, nidif, ws, tables, inflgsw,
                iceflgsw, liqflgsw);

  // RRTMG_SWRAD outputs (level index 0 = surface, nlayers = top)
  gsw = swdflx[0] - swuflx[0];
  swcf = (swdflx[nlayers] - swuflx[nlayers]) - (swdflxc[nlayers] - swuflxc[nlayers]);
  swupt = swuflx[nlayers];
  swuptc = swuflxc[nlayers];
  swdnt = swdflx[nlayers];
  swdntc = swdflxc[nlayers];
  swupb = swuflx[0];
  swupbc = swuflxc[0];
  swdnb = swdflx[0];
  swdnbc = swdflxc[0];
  swvisdir = uvdir[0];
  swvisdif = uvdif[0];
  swnirdir = nidir[0];
  swnirdif = nidif[0];
  for (int k = 0; k < num_levels; ++k) {
    rthratensw[k] = (swhr[k] / 86400.) / pi3d[k];
    rthratenswc[k] = (swhrc[k] / 86400.) / pi3d[k];
  }
}

void rrtmg_swrad(
    int num_cols,
    int nz,
    const RrtmgSwConfig& config,
    const RrtmgSwTables& tables,
    const double* coszen,
    const double* albedo,
    const double* xland,
    const double* xice,
    const double* snow,
    const double* xlat,
    const double* p8w,
    const double* t8w,
    const double* p3d,
    const double* pi3d,
    const double* dz8w,
    const double* t3d,
    const double* rho3d,
    const double* qv3d,
    const double* qc3d,
    const double* qr3d,
    const double* qi3d,
    const double* qs3d,
    const double* qg3d,
    const double* cldfra3d,
    double* gsw,
    double* swcf,
    double* swupt,
    double* swuptc,
    double* swdnt,
    double* swdntc,
    double* swupb,
    double* swupbc,
    double* swdnb,
    double* swdnbc,
    double* swvisdir,
    double* swvisdif,
    double* swnirdir,
    double* swnirdif,
    double* rthratensw,
    double* rthratenswc
) {
#pragma omp parallel
  {
    RrtmgSwWorkspace ws;

#pragma omp for schedule(static)
    for (int col = 0; col < num_cols; ++col) {
      double col_p8w[RRTMG_SW_MAX_LAYERS + 1];
      double col_t8w[RRTMG_SW_MAX_LAYERS + 1];

      for (int k = 0; k < nz + 1; ++k) {
        int idx = k * num_cols + col;
        col_p8w[k] = p8w[idx];
        col_t8w[k] = t8w[idx];
      }

      double col_p3d[RRTMG_SW_MAX_LAYERS];
      double col_pi3d[RRTMG_SW_MAX_LAYERS];
      double col_dz8w[RRTMG_SW_MAX_LAYERS];
      double col_t3d[RRTMG_SW_MAX_LAYERS];
      double col_rho3d[RRTMG_SW_MAX_LAYERS];
      double col_qv3d[RRTMG_SW_MAX_LAYERS];
      double col_qc3d[RRTMG_SW_MAX_LAYERS];
      double col_qr3d[RRTMG_SW_MAX_LAYERS];
      double col_qi3d[RRTMG_SW_MAX_LAYERS];
      double col_qs3d[RRTMG_SW_MAX_LAYERS];
      double col_qg3d[RRTMG_SW_MAX_LAYERS];
      double col_cldfra3d[RRTMG_SW_MAX_LAYERS];

      for (int k = 0; k < nz; ++k) {
        int idx = k * num_cols + col;
        col_p3d[k] = p3d[idx];
        col_pi3d[k] = pi3d[idx];
        col_dz8w[k] = dz8w[idx];
        col_t3d[k] = t3d[idx];
        col_rho3d[k] = rho3d[idx];
        col_qv3d[k] = qv3d[idx];
        col_qc3d[k] = qc3d[idx];
        col_qr3d[k] = qr3d[idx];
        col_qi3d[k] = qi3d[idx];
        col_qs3d[k] = qs3d[idx];
        col_qg3d[k] = qg3d[idx];
        col_cldfra3d[k] = cldfra3d ? cldfra3d[idx] : 0.0;
      }

      double col_rthratensw[RRTMG_SW_MAX_LAYERS];
      double col_rthratenswc[RRTMG_SW_MAX_LAYERS];

      rrtmg_sw_column(
          nz, config, tables,
          coszen[col], albedo[col], xland[col], xice[col], snow[col], xlat[col],
          col_p8w, col_t8w, col_p3d, col_pi3d, col_dz8w, col_t3d, col_rho3d,
          col_qv3d, col_qc3d, col_qr3d, col_qi3d, col_qs3d, col_qg3d, col_cldfra3d, nullptr,
          gsw[col], swcf[col], swupt[col], swuptc[col], swdnt[col], swdntc[col],
          swupb[col], swupbc[col], swdnb[col], swdnbc[col],
          swvisdir[col], swvisdif[col], swnirdir[col], swnirdif[col],
          col_rthratensw, col_rthratenswc, ws
      );

      for (int k = 0; k < nz; ++k) {
        int idx = k * num_cols + col;
        rthratensw[idx] = col_rthratensw[k];
        rthratenswc[idx] = col_rthratenswc[k];
      }
    }
  }
}

} // namespace cppwrf::phys
