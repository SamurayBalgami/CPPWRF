#include <cstdio>
#include "phys/ra_rrtmg_lw/rrtmg_lw.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace cppwrf::phys {

static constexpr double PPROF[60] = {
    1000.00, 855.47, 731.82, 626.05, 535.57, 458.16, 391.94, 335.29, 286.83, 245.38,
    209.91, 179.57, 153.62, 131.41, 112.42, 96.17,  82.27,  70.38,  60.21,  51.51,
    44.06,   37.69,  32.25,  27.59,  23.60,  20.19,  17.27,  14.77,  12.64,  10.81,
    9.25,    7.91,   6.77,   5.79,   4.95,   4.24,   3.63,   3.10,   2.65,   2.27,
    1.94,    1.66,   1.42,   1.22,   1.04,   0.89,   0.76,   0.65,   0.56,   0.48,
    0.41,    0.35,   0.30,   0.26,   0.22,   0.19,   0.16,   0.14,   0.12,   0.10};

static constexpr double TPROF[60] = {
    286.96, 281.07, 275.16, 268.11, 260.56, 253.02, 245.62, 238.41, 231.57, 225.91,
    221.72, 217.79, 215.06, 212.74, 210.25, 210.16, 210.69, 212.14, 213.74, 215.37,
    216.82, 217.94, 219.03, 220.18, 221.37, 222.64, 224.16, 225.88, 227.63, 229.51,
    231.50, 233.73, 236.18, 238.78, 241.60, 244.44, 247.35, 250.33, 253.32, 256.30,
    259.22, 262.12, 264.80, 266.50, 267.59, 268.44, 268.69, 267.76, 266.13, 263.96,
    261.54, 258.93, 256.15, 253.23, 249.89, 246.67, 243.48, 240.25, 236.66, 233.86};

static constexpr double retab[95] = {
    5.92779, 6.26422, 6.61973, 6.99539, 7.39234, 7.81177, 8.25496, 8.72323, 9.21800, 9.74075,
    10.2930, 10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319, 15.0351, 15.8799, 16.7674,
    17.6986, 18.6744, 19.6955, 20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125, 27.7895,
    28.6450, 29.4167, 30.1088, 30.7306, 31.2943, 31.8151, 32.3077, 32.7870, 33.2657, 33.7540,
    34.2601, 34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078, 38.4720, 39.1508, 39.8442,
    40.5552, 41.2912, 42.0635, 42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221, 50.4710,
    52.4980, 54.8315, 57.4898, 60.4785, 63.7898, 65.5604, 71.2885, 75.4113, 79.7368, 84.2351,
    88.8833, 93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424, 124.954, 130.630, 136.457,
    142.446, 148.608, 154.956, 161.503, 168.262, 175.248, 182.473, 189.952, 197.699, 205.728,
    214.055, 222.694, 231.661, 240.971, 250.639};

void rrtmg_lw_inirad(double* o3mmr, const double* plev, int kts, int kte) {
  static constexpr double o3sum[31] = {
      5.297e-8, 5.852e-8, 6.579e-8, 7.505e-8, 8.577e-8, 9.895e-8, 1.175e-7, 1.399e-7,
      1.677e-7, 2.003e-7, 2.571e-7, 3.325e-7, 4.438e-7, 6.255e-7, 8.168e-7, 1.036e-6,
      1.366e-6, 1.855e-6, 2.514e-6, 3.240e-6, 4.033e-6, 4.854e-6, 5.517e-6, 6.089e-6,
      6.689e-6, 1.106e-5, 1.462e-5, 1.321e-5, 9.856e-6, 5.960e-6, 5.960e-6};
  static constexpr double ppsum[31] = {
      955.890, 850.532, 754.599, 667.742, 589.841, 519.421, 455.480, 398.085,
      347.171, 301.735, 261.310, 225.360, 193.419, 165.490, 141.032, 120.125,
      102.689, 87.829,  75.123,  64.306,  55.086,  47.209,  40.535,  34.795,
      29.865,  19.122,  9.277,   4.660,   2.421,   1.294,   0.647};
  static constexpr double o3win[31] = {
      4.629e-8, 4.686e-8, 5.017e-8, 5.613e-8, 6.871e-8, 8.751e-8, 1.138e-7, 1.516e-7,
      2.161e-7, 3.264e-7, 4.968e-7, 7.338e-7, 1.017e-6, 1.308e-6, 1.625e-6, 2.011e-6,
      2.516e-6, 3.130e-6, 3.840e-6, 4.703e-6, 5.486e-6, 6.289e-6, 6.993e-6, 7.494e-6,
      8.197e-6, 9.632e-6, 1.113e-5, 1.146e-5, 9.389e-6, 6.135e-6, 6.135e-6};
  static constexpr double ppwin[31] = {
      955.747, 841.783, 740.199, 649.538, 568.404, 495.815, 431.069, 373.464,
      322.354, 277.190, 237.635, 203.433, 174.070, 148.949, 127.408, 108.915,
      93.114,  79.551,  67.940,  58.072,  49.593,  42.318,  36.138,  30.907,
      26.362,  16.423,  7.583,   3.620,   1.807,   0.938,   0.469};

  double o3ann[31];
  o3ann[0] = 0.5 * (o3sum[0] + o3win[0]);
  for (int k = 1; k < 31; ++k) {
    o3ann[k] = o3win[k - 1] + (o3win[k] - o3win[k - 1]) / (ppwin[k] - ppwin[k - 1]) * (ppsum[k] - ppwin[k - 1]);
  }
  for (int k = 1; k < 31; ++k) {
    o3ann[k] = 0.5 * (o3ann[k] + o3sum[k]);
  }

  double ppwrkh[32];
  ppwrkh[0] = 1100.0;
  for (int k = 1; k < 31; ++k) {
    ppwrkh[k] = 0.5 * (ppsum[k] + ppsum[k - 1]);
  }
  ppwrkh[31] = 0.0;

  for (int k = kts; k <= kte + 1; ++k) {
    o3mmr[k] = 0.0;
    double pr_k = plev[k];
    double pr_k1 = plev[k + 1];
    for (int jj = 0; jj < 31; ++jj) {
      double pb1 = (pr_k > ppwrkh[jj]) ? (pr_k - ppwrkh[jj]) : 0.0;
      double pb2 = (pr_k > ppwrkh[jj + 1]) ? (pr_k - ppwrkh[jj + 1]) : 0.0;
      double pt1 = (pr_k1 > ppwrkh[jj]) ? (pr_k1 - ppwrkh[jj]) : 0.0;
      double pt2 = (pr_k1 > ppwrkh[jj + 1]) ? (pr_k1 - ppwrkh[jj + 1]) : 0.0;
      o3mmr[k] += (pb2 - pb1 - pt2 + pt1) * o3ann[jj];
    }
    o3mmr[k] /= (pr_k - pr_k1);
  }
}

// WRF rrtmg_lw (inatm, cldprmc, setcoef, taumol, rtrnmc) for one column whose McICA subcolumn
// cloud arrays are already in ws (cldfmc, ciwpmc, clwpmc, cswpmc). Layer arrays have nlayers
// entries, level arrays nlayers + 1 (index 0 = surface); pressures in mb, emis per band.
void rrtmg_default_gases(RrtmgLwConfig& c) {
  c.co2 = (280. + 90. * std::exp(0.02 * (c.yr - 2000))) * 1.e-6;
  c.ch4 = 1774.e-9;
  c.n2o = 319.e-9;
  c.cfc11 = 0.251e-9;
  c.cfc12 = 0.538e-9;
}

std::FILE* g_rrtmg_lw_dbg = nullptr;  // DEBUG (not for main)

void rrtmg_lw_core(int nlayers, double cp, const double* play, const double* plev, const double* tlay,
                   const double* tlev, double tsk, const double* h2ovmr, const double* o3vmr, const double* co2vmr,
                   const double* ch4vmr, const double* n2ovmr, const double* o2vmr, const double* cfc11vmr,
                   const double* cfc12vmr, const double* cfc22vmr, const double* ccl4vmr, const double* semiss,
                   const double* rei, const double* rel, const double* res, double* totuflux, double* totdflux,
                   double* htr, double* totuclfl, double* totdclfl, double* htrc, RrtmgLwWorkspace& ws,
                   const RrtmgLwTables& tables) {
  // 7. Atmospheric preparation (inatm)
  constexpr double amd = 28.9660;
  constexpr double amw = 18.0160;
  constexpr double avogad = 6.02214199e23;
  constexpr double grav = 9.8066;

  double pavel[RRTMG_LW_MAX_LAYERS];
  double tavel[RRTMG_LW_MAX_LAYERS];
  double pz[RRTMG_LW_MAX_LAYERS + 1];
  double tz[RRTMG_LW_MAX_LAYERS + 1];
  double coldry[RRTMG_LW_MAX_LAYERS];
  double wbrodl[RRTMG_LW_MAX_LAYERS];

  pz[0] = plev[0];
  tz[0] = tlev[0];
  for (int l = 0; l < nlayers; ++l) {
    pavel[l] = play[l];
    tavel[l] = tlay[l];
    pz[l + 1] = plev[l + 1];
    tz[l + 1] = tlev[l + 1];

    ws.wkl[0][l] = h2ovmr[l];
    ws.wkl[1][l] = co2vmr[l];
    ws.wkl[2][l] = o3vmr[l];
    ws.wkl[3][l] = n2ovmr[l];
    ws.wkl[4][l] = 0.0; // co
    ws.wkl[5][l] = ch4vmr[l];
    ws.wkl[6][l] = o2vmr[l];

    double amm = (1.0 - ws.wkl[0][l]) * amd + ws.wkl[0][l] * amw;
    coldry[l] = (pz[l] - pz[l + 1]) * 1.0e3 * avogad / (1.0e2 * grav * amm * (1.0 + ws.wkl[0][l]));
  }

  for (int l = 0; l < nlayers; ++l) {
    ws.wx[0][l] = ccl4vmr[l];
    ws.wx[1][l] = cfc11vmr[l];
    ws.wx[2][l] = cfc12vmr[l];
    ws.wx[3][l] = cfc22vmr[l];
  }

  double amttl = 0.0;
  double wvttl = 0.0;
  for (int l = 0; l < nlayers; ++l) {
    double summol = 0.0;
    for (int imol = 1; imol < 7; ++imol) {
      summol += ws.wkl[imol][l];
    }
    wbrodl[l] = coldry[l] * (1.0 - summol);
    for (int imol = 0; imol < 7; ++imol) {
      ws.wkl[imol][l] = coldry[l] * ws.wkl[imol][l];
    }
    amttl = amttl + coldry[l] + ws.wkl[0][l];  // WRF order: (amttl + coldry) + wkl
    wvttl += ws.wkl[0][l];

    for (int ix = 0; ix < 4; ++ix) {
      ws.wx[ix][l] = coldry[l] * ws.wx[ix][l] * 1.0e-20;
    }
  }

  double wvsh = (amw * wvttl) / (amd * amttl);
  double pwvcm = wvsh * (1.0e3 * pz[0]) / (1.0e2 * grav);

  // 8. In-cloud optical depths calculation (cldprmc)
  rrtmg_lw_cldprmc(nlayers, 2, 3, 1, ws.cldfmc_ptrs, ws.ciwpmc_ptrs, ws.clwpmc_ptrs, ws.cswpmc_ptrs,
                   rei, rel, res, ws.taucmc_ptrs, tables);

  // 9. Radiative transfer coefficients and Planck functions (setcoef)
  int laytrop;
  int jp[RRTMG_LW_MAX_LAYERS];
  int jt[RRTMG_LW_MAX_LAYERS];
  int jt1[RRTMG_LW_MAX_LAYERS];
  double plankbnd[16];
  double colh2o[RRTMG_LW_MAX_LAYERS];
  double colco2[RRTMG_LW_MAX_LAYERS];
  double colo3[RRTMG_LW_MAX_LAYERS];
  double coln2o[RRTMG_LW_MAX_LAYERS];
  double colco[RRTMG_LW_MAX_LAYERS];
  double colch4[RRTMG_LW_MAX_LAYERS];
  double colo2[RRTMG_LW_MAX_LAYERS];
  double colbrd[RRTMG_LW_MAX_LAYERS];
  double fac00[RRTMG_LW_MAX_LAYERS];
  double fac01[RRTMG_LW_MAX_LAYERS];
  double fac10[RRTMG_LW_MAX_LAYERS];
  double fac11[RRTMG_LW_MAX_LAYERS];
  double rat_h2oco2[RRTMG_LW_MAX_LAYERS];
  double rat_h2oco2_1[RRTMG_LW_MAX_LAYERS];
  double rat_h2oo3[RRTMG_LW_MAX_LAYERS];
  double rat_h2oo3_1[RRTMG_LW_MAX_LAYERS];
  double rat_h2on2o[RRTMG_LW_MAX_LAYERS];
  double rat_h2on2o_1[RRTMG_LW_MAX_LAYERS];
  double rat_h2och4[RRTMG_LW_MAX_LAYERS];
  double rat_h2och4_1[RRTMG_LW_MAX_LAYERS];
  double rat_n2oco2[RRTMG_LW_MAX_LAYERS];
  double rat_n2oco2_1[RRTMG_LW_MAX_LAYERS];
  double rat_o3co2[RRTMG_LW_MAX_LAYERS];
  double rat_o3co2_1[RRTMG_LW_MAX_LAYERS];
  double selffac[RRTMG_LW_MAX_LAYERS];
  double selffrac[RRTMG_LW_MAX_LAYERS];
  int indself[RRTMG_LW_MAX_LAYERS];
  double forfac[RRTMG_LW_MAX_LAYERS];
  double forfrac[RRTMG_LW_MAX_LAYERS];
  int indfor[RRTMG_LW_MAX_LAYERS];
  double minorfrac[RRTMG_LW_MAX_LAYERS];
  double scaleminor[RRTMG_LW_MAX_LAYERS];
  double scaleminorn2[RRTMG_LW_MAX_LAYERS];
  int indminor[RRTMG_LW_MAX_LAYERS];

  rrtmg_lw_setcoef(nlayers, 1, pavel, tavel, tz, tsk, semiss, coldry,
                   const_cast<const double* const*>(ws.wkl_ptrs), wbrodl,
                   laytrop, jp, jt, jt1, ws.planklay_ptrs, ws.planklev_ptrs, plankbnd,
                   colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd,
                   fac00, fac01, fac10, fac11,
                   rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1,
                   rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1,
                   rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1,
                   selffac, selffrac, indself, forfac, forfrac, indfor,
                   minorfrac, scaleminor, scaleminorn2, indminor, tables);

  // 10. Gaseous absorption optical depths (taumol)
  rrtmg_lw_taumol(nlayers, pavel, const_cast<const double* const*>(ws.wx_ptrs), coldry,
                  laytrop, jp, jt, jt1,
                  colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd,
                  fac00, fac01, fac10, fac11,
                  rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1,
                  rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1,
                  rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1,
                  selffac, selffrac, indself, forfac, forfrac, indfor,
                  minorfrac, scaleminor, scaleminorn2, indminor,
                  ws.fracs_ptrs, ws.taug_ptrs, tables);

  if (g_rrtmg_lw_dbg) {  // DEBUG: same order as the instrumented Fortran
    std::FILE* f = g_rrtmg_lw_dbg;
    auto wi = [&](int v) { std::fwrite(&v, 4, 1, f); };
    auto wd = [&](double v) { std::fwrite(&v, 8, 1, f); };
    wi(nlayers); wi(laytrop + 0); wd(pwvcm);
    for (int k = 0; k < nlayers; ++k) wd(coldry[k]);
    for (int k = 0; k < nlayers; ++k) wd(wbrodl[k]);
    for (int k = 0; k < nlayers; ++k) for (int m = 0; m < 7; ++m) wd(ws.wkl[m][k]);
    for (int k = 0; k < nlayers; ++k) for (int m = 0; m < 4; ++m) wd(ws.wx[m][k]);
    for (int k = 0; k < nlayers; ++k) wi(jp[k]);
    for (int k = 0; k < nlayers; ++k) wi(jt[k]);
    for (int k = 0; k < nlayers; ++k) wi(jt1[k]);
    for (double* a : {fac00, fac01, fac10, fac11, selffac, forfac, minorfrac, scaleminor})
      for (int k = 0; k < nlayers; ++k) wd(a[k]);
    for (int k = 0; k < nlayers; ++k) for (int m = 0; m < 16; ++m) wd(ws.planklay[k][m]);
    for (int k = 0; k <= nlayers; ++k) for (int m = 0; m < 16; ++m) wd(ws.planklev[k][m]);
    for (int m = 0; m < 16; ++m) wd(plankbnd[m]);
    for (int k = 0; k < nlayers; ++k) for (int m = 0; m < 140; ++m) wd(ws.fracs[k][m]);
    for (int k = 0; k < nlayers; ++k) for (int m = 0; m < 140; ++m) wd(ws.taug[k][m]);
  }
  // Combine gaseous and aerosol optical depths (aerosol optical depth is 0)
  for (int l = 0; l < nlayers; ++l) {
    for (int ig = 0; ig < 140; ++ig) {
      ws.taut[l][ig] = ws.taug[l][ig];
    }
  }

  // 11. McICA Radiative transfer solver (rtrnmc)
  double fnet[RRTMG_LW_MAX_LAYERS + 1];
  double fnetc[RRTMG_LW_MAX_LAYERS + 1];

  double heatfac = grav * 86400.0 / (cp * 100.0);

  rrtmg_lw_rtrnmc(nlayers, 1, 16, 0, pz, semiss, 16,
                  const_cast<const double* const*>(ws.cldfmc_ptrs),
                  const_cast<const double* const*>(ws.taucmc_ptrs),
                  const_cast<const double* const*>(ws.planklay_ptrs),
                  const_cast<const double* const*>(ws.planklev_ptrs),
                  plankbnd, pwvcm,
                  const_cast<const double* const*>(ws.fracs_ptrs),
                  const_cast<const double* const*>(ws.taut_ptrs),
                  totuflux, totdflux, fnet, htr,
                  totuclfl, totdclfl, fnetc, htrc,
                  heatfac, tables);

  if (g_rrtmg_lw_dbg) {
    for (int k = 0; k <= nlayers; ++k) std::fwrite(&totuflux[k], 8, 1, g_rrtmg_lw_dbg);
    for (int k = 0; k <= nlayers; ++k) std::fwrite(&totdflux[k], 8, 1, g_rrtmg_lw_dbg);
    for (int k = 0; k < nlayers; ++k) std::fwrite(&htr[k], 8, 1, g_rrtmg_lw_dbg);
    std::fclose(g_rrtmg_lw_dbg);
    g_rrtmg_lw_dbg = nullptr;
  }
}

void rrtmg_lw_column(
    int num_levels,
    const RrtmgLwConfig& config,
    const RrtmgLwTables& tables,
    double tsk,
    double emiss,
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
    double& glw,
    double& olr,
    double& lwcf,
    double& lwupt,
    double& lwuptc,
    double& lwdnt,
    double& lwdntc,
    double& lwupb,
    double& lwupbc,
    double& lwdnb,
    double& lwdnbc,
    double* rthratenlw,
    double* rthratenlwc,
    RrtmgLwWorkspace& ws
) {
  constexpr double deltap = 4.0; // mb
  // rrtmg_lwinit: NLAYERS = kme + nint(p_top*0.01/deltap) - 1 with kme = nz + 1 (Fortran NINT)
  const int nlayers = num_levels + int(std::lround(config.p_top * 0.01 / deltap));

  // 1. Compute effective droplet and ice sizes (CAM relcalc / reicalc)
  double reliq[RRTMG_LW_MAX_LAYERS];
  double reice[RRTMG_LW_MAX_LAYERS];

  double landfrac = 2.0 - xland;
  double landm = landfrac;
  double snowh = 0.001 * snow;
  double icefrac = xice;
  constexpr double rliqland = 8.0, rliqocean = 14.0, rliqice = 14.0, tmelt = 273.16;

  for (int k = 0; k < num_levels; ++k) {
    double r = rliqland + (rliqocean - rliqland) * std::min(1.0, std::max(0.0, (tmelt - t3d[k]) * 0.05));
    r += (rliqocean - r) * std::min(1.0, std::max(0.0, snowh * 10.0));
    r += (rliqocean - r) * std::min(1.0, std::max(0.0, 1.0 - landm));
    r += (rliqice - r) * std::min(1.0, std::max(0.0, icefrac));
    reliq[k] = r;

    int idx = static_cast<int>(t3d[k] - 179.0);
    idx = std::min(std::max(idx, 1), 94);
    double corr = t3d[k] - static_cast<int>(t3d[k]);
    double ri = retab[idx - 1] * (1.0 - corr) + retab[idx] * corr;
    // Fu generalized effective size conversion for iceflglw = 3
    ri = ri * 1.0315;
    if (ri > 140.0) ri = 140.0;
    reice[k] = ri;
  }

  // 2. Base layer vertical variables
  double play[RRTMG_LW_MAX_LAYERS];
  double plev[RRTMG_LW_MAX_LAYERS + 1];
  double tlay[RRTMG_LW_MAX_LAYERS];
  double tlev[RRTMG_LW_MAX_LAYERS + 1];
  double pdel[RRTMG_LW_MAX_LAYERS];
  double hgt[RRTMG_LW_MAX_LAYERS];

  double h2ovmr[RRTMG_LW_MAX_LAYERS];
  double co2vmr[RRTMG_LW_MAX_LAYERS];
  double o2vmr[RRTMG_LW_MAX_LAYERS];
  double ch4vmr[RRTMG_LW_MAX_LAYERS];
  double n2ovmr[RRTMG_LW_MAX_LAYERS];
  double cfc11vmr[RRTMG_LW_MAX_LAYERS];
  double cfc12vmr[RRTMG_LW_MAX_LAYERS];
  double cfc22vmr[RRTMG_LW_MAX_LAYERS];
  double ccl4vmr[RRTMG_LW_MAX_LAYERS];

  double clwpth[RRTMG_LW_MAX_LAYERS];
  double ciwpth[RRTMG_LW_MAX_LAYERS];
  double cswpth[RRTMG_LW_MAX_LAYERS];
  double rel[RRTMG_LW_MAX_LAYERS];
  double rei[RRTMG_LW_MAX_LAYERS];
  double res[RRTMG_LW_MAX_LAYERS];
  double cldfrac[RRTMG_LW_MAX_LAYERS];

  plev[0] = p8w[0] / 100.0;
  tlev[0] = t8w[0];
  double dzsum = 0.0;

  for (int k = 0; k < num_levels; ++k) {
    play[k] = p3d[k] / 100.0;
    plev[k + 1] = p8w[k + 1] / 100.0;
    pdel[k] = plev[k] - plev[k + 1];
    tlay[k] = t3d[k];
    tlev[k + 1] = t8w[k + 1];

    double qv_val = std::max(qv3d[k], 1.0e-12);
    h2ovmr[k] = qv_val * 1.607793;
    co2vmr[k] = config.co2;
    o2vmr[k] = 0.209488;
    ch4vmr[k] = config.ch4;
    n2ovmr[k] = config.n2o;
    cfc11vmr[k] = config.cfc11;
    cfc12vmr[k] = config.cfc12;
    cfc22vmr[k] = 0.169e-9;
    ccl4vmr[k] = 0.093e-9;

    hgt[k] = dzsum + 0.5 * dz8w[k];
    dzsum += dz8w[k];

    // Grid box water paths (g/m2)
    double qc_val = std::max(0.0, qc3d[k]);
    double qi_val = std::max(0.0, qi3d[k]);
    double qs_val = std::max(0.0, qs3d[k]);
    double cf_val = cldfra3d[k];

    double gicewp = (qi_val + qs_val) * pdel[k] * 100.0 / config.g * 1000.0;
    double gliqwp = qc_val * pdel[k] * 100.0 / config.g * 1000.0;

    double denom = std::max(0.01, cf_val);
    clwpth[k] = gliqwp / denom;
    ciwpth[k] = gicewp / denom;
    cswpth[k] = 0.0;
    rel[k] = reliq[k];
    rei[k] = reice[k];
    res[k] = 10.0;
    cldfrac[k] = cf_val;
  }

  // 3. Buffer layer expansion above model top (to 0 mb)
  double dz_top = dz8w[num_levels - 1];
  for (int L = num_levels; L < nlayers; ++L) {
    plev[L + 1] = plev[L] - deltap;
    play[L] = 0.5 * (plev[L] + plev[L + 1]);
    hgt[L] = dzsum + 0.5 * dz_top;
    dzsum += dz_top;
  }
  plev[nlayers] = 0.00;
  play[nlayers - 1] = 0.5 * (plev[nlayers - 1] + plev[nlayers]);

  // Interpolate standard profile temperatures to column pressure levels
  double varint[RRTMG_LW_MAX_LAYERS + 1];
  for (int L = 0; L <= nlayers; ++L) {
    double p_val = plev[L];
    int klev = 59;
    if (PPROF[59] < p_val) {
      for (int LL = 1; LL < 60; ++LL) {
        if (PPROF[LL] < p_val) {
          klev = LL - 1;
          break;
        }
      }
    }
    double wght = 0.0;
    double vark = TPROF[klev];
    double vark1 = (klev < 59) ? TPROF[klev + 1] : TPROF[klev];
    if (klev < 59) {
      wght = (p_val - PPROF[klev]) / (PPROF[klev + 1] - PPROF[klev]);
    }
    varint[L] = wght * (vark1 - vark) + vark;
  }

  for (int L = num_levels; L <= nlayers; ++L) {
    tlev[L] = varint[L] + (tlev[num_levels - 1] - varint[num_levels - 1]);
    tlay[L - 1] = 0.5 * (tlev[L] + tlev[L - 1]);
  }

  for (int L = num_levels; L < nlayers; ++L) {
    h2ovmr[L] = h2ovmr[num_levels - 1];
    co2vmr[L] = co2vmr[num_levels - 1];
    o2vmr[L] = o2vmr[num_levels - 1];
    ch4vmr[L] = ch4vmr[num_levels - 1];
    n2ovmr[L] = n2ovmr[num_levels - 1];
    cfc11vmr[L] = cfc11vmr[num_levels - 1];
    cfc12vmr[L] = cfc12vmr[num_levels - 1];
    cfc22vmr[L] = cfc22vmr[num_levels - 1];
    ccl4vmr[L] = ccl4vmr[num_levels - 1];
    clwpth[L] = 0.0;
    ciwpth[L] = 0.0;
    cswpth[L] = 0.0;
    rel[L] = 10.0;
    rei[L] = 10.0;
    res[L] = 10.0;
    cldfrac[L] = 0.0;
  }

  // 4. Ozone profile calculation via inirad
  double o3mmr[RRTMG_LW_MAX_LAYERS + 1];
  rrtmg_lw_inirad(o3mmr, plev, 0, nlayers - 1);
  double o3vmr[RRTMG_LW_MAX_LAYERS];
  constexpr double amdo = 0.603461;
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

  // 5. Surface emissivity
  double semiss[16];
  for (int nb = 0; nb < 16; ++nb) {
    semiss[nb] = emiss;
  }

  // 6. McICA stochastic subcolumn generator (KISSVEC PRNG, changeSeed = 150)
  KissRng rng;
  rng.init_from_pressure(play[0] * 100.0, play[1] * 100.0, play[2] * 100.0, play[3] * 100.0);
  rng.permute(150);

  double cldf[RRTMG_LW_MAX_LAYERS];
  for (int ilev = 0; ilev < nlayers; ++ilev) {
    cldf[ilev] = (cldfrac[ilev] < 1.0e-20) ? 0.0 : cldfrac[ilev];
  }

  double cdf[140][RRTMG_LW_MAX_LAYERS];
  for (int isubcol = 0; isubcol < 140; ++isubcol) {
    for (int ilev = 0; ilev < nlayers; ++ilev) {
      cdf[isubcol][ilev] = rng.step<double>();
    }
  }

  // Max/random overlap (cldovrlp == 2)
  for (int ilev = 1; ilev < nlayers; ++ilev) {
    for (int isubcol = 0; isubcol < 140; ++isubcol) {
      if (cdf[isubcol][ilev - 1] > 1.0 - cldf[ilev - 1]) {
        cdf[isubcol][ilev] = cdf[isubcol][ilev - 1];
      } else {
        cdf[isubcol][ilev] = cdf[isubcol][ilev] * (1.0 - cldf[ilev - 1]);
      }
    }
  }

  for (int ilev = 0; ilev < nlayers; ++ilev) {
    for (int isubcol = 0; isubcol < 140; ++isubcol) {
      bool is_cloudy = (cdf[isubcol][ilev] >= 1.0 - cldf[ilev]);
      if (is_cloudy) {
        ws.cldfmc[isubcol][ilev] = 1.0;
        ws.clwpmc[isubcol][ilev] = clwpth[ilev];
        ws.ciwpmc[isubcol][ilev] = ciwpth[ilev];
        ws.cswpmc[isubcol][ilev] = cswpth[ilev];
        ws.taucmc[isubcol][ilev] = 0.0;
      } else {
        ws.cldfmc[isubcol][ilev] = 0.0;
        ws.clwpmc[isubcol][ilev] = 0.0;
        ws.ciwpmc[isubcol][ilev] = 0.0;
        ws.cswpmc[isubcol][ilev] = 0.0;
        ws.taucmc[isubcol][ilev] = 0.0;
      }
    }
  }

  double totuflux[RRTMG_LW_MAX_LAYERS + 1];
  double totdflux[RRTMG_LW_MAX_LAYERS + 1];
  double htr[RRTMG_LW_MAX_LAYERS + 1];
  double totuclfl[RRTMG_LW_MAX_LAYERS + 1];
  double totdclfl[RRTMG_LW_MAX_LAYERS + 1];
  double htrc[RRTMG_LW_MAX_LAYERS + 1];
  rrtmg_lw_core(nlayers, config.cp, play, plev, tlay, tlev, tsk, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,
                cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, semiss, rei, rel, res, totuflux, totdflux, htr, totuclfl,
                totdclfl, htrc, ws, tables);

  // 12. Map outputs
  glw = totdflux[0];
  olr = totuflux[nlayers];
  lwcf = totuclfl[nlayers] - totuflux[nlayers];

  lwupt = totuflux[nlayers];
  lwuptc = totuclfl[nlayers];
  lwdnt = totdflux[nlayers];
  lwdntc = totdclfl[nlayers];
  lwupb = totuflux[0];
  lwupbc = totuclfl[0];
  lwdnb = totdflux[0];
  lwdnbc = totdclfl[0];

  for (int k = 0; k < num_levels; ++k) {
    rthratenlw[k] = (htr[k] / 86400.0) / pi3d[k];
    rthratenlwc[k] = (htrc[k] / 86400.0) / pi3d[k];
  }
}

void rrtmg_lwrad(
    int num_cols,
    int nz,
    const RrtmgLwConfig& config,
    const RrtmgLwTables& tables,
    const double* tsk,
    const double* emiss,
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
    double* glw,
    double* olr,
    double* lwcf,
    double* lwupt,
    double* lwuptc,
    double* lwdnt,
    double* lwdntc,
    double* lwupb,
    double* lwupbc,
    double* lwdnb,
    double* lwdnbc,
    double* rthratenlw,
    double* rthratenlwc
) {
  int num_threads = 1;
#ifdef _OPENMP
  num_threads = omp_get_max_threads();
#endif

  std::vector<RrtmgLwWorkspace> workspaces(num_threads);
  for (auto& ws : workspaces) {
    ws.init_pointers();
  }

#pragma omp parallel for schedule(dynamic, 1)
  for (int col = 0; col < num_cols; ++col) {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    RrtmgLwWorkspace& ws = workspaces[tid];

    // Gather 1D profiles for this column
    double col_p8w[RRTMG_LW_MAX_LAYERS + 1];
    double col_t8w[RRTMG_LW_MAX_LAYERS + 1];
    for (int k = 0; k <= nz; ++k) {
      col_p8w[k] = p8w[k * num_cols + col];
      col_t8w[k] = t8w[k * num_cols + col];
    }

    double col_p3d[RRTMG_LW_MAX_LAYERS];
    double col_pi3d[RRTMG_LW_MAX_LAYERS];
    double col_dz8w[RRTMG_LW_MAX_LAYERS];
    double col_t3d[RRTMG_LW_MAX_LAYERS];
    double col_rho3d[RRTMG_LW_MAX_LAYERS];
    double col_qv3d[RRTMG_LW_MAX_LAYERS];
    double col_qc3d[RRTMG_LW_MAX_LAYERS];
    double col_qr3d[RRTMG_LW_MAX_LAYERS];
    double col_qi3d[RRTMG_LW_MAX_LAYERS];
    double col_qs3d[RRTMG_LW_MAX_LAYERS];
    double col_qg3d[RRTMG_LW_MAX_LAYERS];
    double col_cldfra3d[RRTMG_LW_MAX_LAYERS];

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
      col_cldfra3d[k] = cldfra3d[idx];
    }

    double col_rthratenlw[RRTMG_LW_MAX_LAYERS];
    double col_rthratenlwc[RRTMG_LW_MAX_LAYERS];

    rrtmg_lw_column(
        nz, config, tables,
        tsk[col], emiss[col], xland[col], xice[col], snow[col], xlat[col],
        col_p8w, col_t8w, col_p3d, col_pi3d, col_dz8w, col_t3d, col_rho3d,
        col_qv3d, col_qc3d, col_qr3d, col_qi3d, col_qs3d, col_qg3d, col_cldfra3d, nullptr,
        glw[col], olr[col], lwcf[col],
        lwupt[col], lwuptc[col], lwdnt[col], lwdntc[col],
        lwupb[col], lwupbc[col], lwdnb[col], lwdnbc[col],
        col_rthratenlw, col_rthratenlwc,
        ws
    );

    for (int k = 0; k < nz; ++k) {
      int idx = k * num_cols + col;
      rthratenlw[idx] = col_rthratenlw[k];
      rthratenlwc[idx] = col_rthratenlwc[k];
    }
  }
}

} // namespace cppwrf::phys
