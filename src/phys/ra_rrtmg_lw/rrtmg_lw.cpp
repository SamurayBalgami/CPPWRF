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

namespace {
// inirad (module_ra_rrtmg_lw.F) climatology
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


}  // namespace

void rrtmg_lw_fill_consts(RrtmgLwConsts& k) {
  for (int i = 0; i < 60; ++i) { k.pprof[i] = PPROF[i]; k.tprof[i] = TPROF[i]; }
  for (int i = 0; i < 95; ++i) k.retab[i] = retab[i];
  // inirad: annual mean ozone (same arithmetic as WRF, once instead of every column)
  double* o3ann = k.o3ann;
  o3ann[0] = 0.5 * (o3sum[0] + o3win[0]);
  for (int kk = 1; kk < 31; ++kk) {
    o3ann[kk] = o3win[kk - 1] + (o3win[kk] - o3win[kk - 1]) / (ppwin[kk] - ppwin[kk - 1]) * (ppsum[kk] - ppwin[kk - 1]);
  }
  for (int kk = 1; kk < 31; ++kk) {
    o3ann[kk] = 0.5 * (o3ann[kk] + o3sum[kk]);
  }
  double* ppwrkh = k.ppwrkh;
  ppwrkh[0] = 1100.0;
  for (int kk = 1; kk < 31; ++kk) {
    ppwrkh[kk] = 0.5 * (ppsum[kk] + ppsum[kk - 1]);
  }
  ppwrkh[31] = 0.0;
}

const RrtmgLwConsts& rrtmg_lw_consts() {
  static const RrtmgLwConsts k = [] {
    RrtmgLwConsts c{};
    rrtmg_lw_fill_consts(c);
    return c;
  }();
  return k;
}

void rrtmg_default_gases(RrtmgLwConfig& c) {
  c.co2 = (280. + 90. * std::exp(0.02 * (c.yr - 2000))) * 1.e-6;
  c.ch4 = 1774.e-9;
  c.n2o = 319.e-9;
  c.cfc11 = 0.251e-9;
  c.cfc12 = 0.538e-9;
}

RrtmgLwCallConfig rrtmg_lw_call_config(int num_levels, const RrtmgLwConfig& config, bool hrc, bool hri, bool hrs) {
  constexpr double deltap = 4.0;  // mb
  constexpr double grav = 9.8066;
  RrtmgLwCallConfig c{};
  c.num_levels = num_levels;
  // rrtmg_lwinit: NLAYERS = kme + nint(p_top*0.01/deltap) - 1 with kme = nz + 1 (Fortran NINT)
  c.nlayers = num_levels + int(std::lround(config.p_top * 0.01 / deltap));
  c.o3input = config.o3input;
  c.co2 = config.co2; c.n2o = config.n2o; c.ch4 = config.ch4; c.cfc11 = config.cfc11; c.cfc12 = config.cfc12;
  c.g = config.g;
  c.heatfac = grav * 86400.0 / (config.cp * 100.0);
  c.has_reqc = hrc; c.has_reqi = hri; c.has_reqs = hrs;
  return c;
}

namespace {
void ensure(RrtmgLwWorkspace& ws, const LwLayout& lay) {
  if (ws.d.size() < lay.nd) ws.d.assign(lay.nd, 0.0);
  if (ws.i.size() < lay.ni) ws.i.assign(lay.ni, 0);
}
SD cv(const double* p) { return SD{const_cast<double*>(p), 1}; }
}  // namespace

// WRF rrtmg_lw (inatm, cldprmc, setcoef, taumol, rtrnmc) for one column whose McICA subcolumn
// cloud arrays are already in ws (cldfmc, ciwpmc, clwpmc, cswpmc, taucmc). Layer arrays have
// nlayers entries, level arrays nlayers + 1 (index 0 = surface); pressures in mb, emis per band.
void rrtmg_lw_core(int nlayers, double cp, const double* play, const double* plev, const double* tlay,
                   const double* tlev, double tsk, const double* h2ovmr, const double* o3vmr, const double* co2vmr,
                   const double* ch4vmr, const double* n2ovmr, const double* o2vmr, const double* cfc11vmr,
                   const double* cfc12vmr, const double* cfc22vmr, const double* ccl4vmr, const double* semiss,
                   const double* rei, const double* rel, const double* res, double* totuflux, double* totdflux,
                   double* htr, double* totuclfl, double* totdclfl, double* htrc, RrtmgLwWorkspace& ws,
                   const RrtmgLwTables& tables, int inflag, int iceflag, int liqflag) {
  constexpr double grav = 9.8066;
  LwLayout lay(nlayers);
  ensure(ws, lay);
  const LwViews v = lay.views(ws.d.data(), ws.i.data(), 0, 1);
  v.pz[0] = plev[0];
  v.tz[0] = tlev[0];
  for (int l = 0; l < nlayers; ++l) {
    v.pavel[l] = play[l];
    v.tavel[l] = tlay[l];
    v.pz[l + 1] = plev[l + 1];
    v.tz[l + 1] = tlev[l + 1];
    v.wkl[0][l] = h2ovmr[l]; v.wkl[1][l] = co2vmr[l]; v.wkl[2][l] = o3vmr[l]; v.wkl[3][l] = n2ovmr[l];
    v.wkl[4][l] = 0.0; v.wkl[5][l] = ch4vmr[l]; v.wkl[6][l] = o2vmr[l];
    v.wx[0][l] = ccl4vmr[l]; v.wx[1][l] = cfc11vmr[l]; v.wx[2][l] = cfc12vmr[l]; v.wx[3][l] = cfc22vmr[l];
    v.reicmc[l] = rei[l]; v.relqmc[l] = rel[l]; v.resnmc[l] = res[l];
    // McICA flags from ws; the subcolumn water paths are not used (the core check has no clouds)
    for (int g = 0; g < 140; ++g) { v.cldfmc[g][l] = ws.cldfmc[g][l]; v.taucmc[g][l] = ws.taucmc[g][l]; }
    v.ciwpth[l] = 0.0; v.clwpth[l] = 0.0; v.cswpth[l] = 0.0;
  }
  for (int nb = 0; nb < 16; ++nb) v.semiss[nb] = semiss[nb];
  v.iscal[1] = inflag; v.iscal[2] = iceflag; v.iscal[3] = liqflag;
  rrtmg_lw_inatm(nlayers, v);
  int laytrop = 0;
  rrtmg_lw_setcoef(nlayers, 1, v.pavel, v.tavel, v.tz, tsk, v.semiss, v.coldry, v.wkl, v.wbrodl, laytrop, v.jp,
                   v.jt, v.jt1, v.planklay, v.planklev, v.plankbnd, v.colh2o, v.colco2, v.colo3, v.coln2o, v.colco,
                   v.colch4, v.colo2, v.colbrd, v.fac00, v.fac01, v.fac10, v.fac11, v.rat_h2oco2, v.rat_h2oco2_1,
                   v.rat_h2oo3, v.rat_h2oo3_1, v.rat_h2on2o, v.rat_h2on2o_1, v.rat_h2och4, v.rat_h2och4_1,
                   v.rat_n2oco2, v.rat_n2oco2_1, v.rat_o3co2, v.rat_o3co2_1, v.selffac, v.selffrac, v.indself,
                   v.forfac, v.forfrac, v.indfor, v.minorfrac, v.scaleminor, v.scaleminorn2, v.indminor, tables);
  v.iscal[0] = laytrop;
  rrtmg_lw_secdiff(v.dscal[0], v.secdiff);
  for (int l = 0; l < nlayers; ++l) rrtmg_lw_layer(nlayers, l, v, tables);
  for (int b = 0; b < 16; ++b) rrtmg_lw_bandstage(nlayers, b, v, tables);
  for (int lev = 0; lev <= nlayers; ++lev) rrtmg_lw_level(lev, v, tables);
  rrtmg_lw_heating(nlayers, v.pz, v.fnet, v.fnetc, v.htr, v.htrc, grav * 86400.0 / (cp * 100.0));
  for (int lev = 0; lev <= nlayers; ++lev) {
    totuflux[lev] = v.totuflux[lev]; totdflux[lev] = v.totdflux[lev];
    totuclfl[lev] = v.totuclfl[lev]; totdclfl[lev] = v.totdclfl[lev];
    htr[lev] = v.htr[lev]; htrc[lev] = v.htrc[lev];
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
    RrtmgLwWorkspace& ws,
    const double* re_cloud,
    const double* re_ice,
    const double* re_snow
) {
  (void)qr3d; (void)qg3d; (void)xlat; (void)rho3d; (void)dz8w;
  const RrtmgLwCallConfig c = rrtmg_lw_call_config(num_levels, config, re_cloud, re_ice, re_snow);
  LwLayout lay(c.nlayers);
  ensure(ws, lay);
  const LwViews v = lay.views(ws.d.data(), ws.i.data(), 0, 1);
  LwColumnIn in{};
  in.tsk = tsk; in.emiss = emiss; in.xland = xland; in.xice = xice; in.snow = snow; in.xlat = xlat;
  in.p8w = cv(p8w); in.t8w = cv(t8w); in.p3d = cv(p3d); in.pi3d = cv(pi3d); in.dz8w = cv(dz8w); in.t3d = cv(t3d);
  in.rho3d = cv(rho3d); in.qv3d = cv(qv3d); in.qc3d = cv(qc3d); in.qi3d = cv(qi3d); in.qs3d = cv(qs3d);
  in.cldfra3d = cv(cldfra3d); in.o33d = cv(o33d ? o33d : qv3d);
  in.re_cloud = cv(re_cloud); in.re_ice = cv(re_ice); in.re_snow = cv(re_snow);
  in.has_qc = in.has_qi = in.has_qs = true;  // absent species are passed as zero arrays
  rrtmg_lw_prep(c, in, v, tables, rrtmg_lw_consts());
  for (int l = 0; l < c.nlayers; ++l) rrtmg_lw_layer(c.nlayers, l, v, tables);
  for (int b = 0; b < 16; ++b) rrtmg_lw_bandstage(c.nlayers, b, v, tables);
  for (int lev = 0; lev <= c.nlayers; ++lev) rrtmg_lw_level(lev, v, tables);
  const LwColumnOut o = rrtmg_lw_finish(c, v, cv(pi3d), SD{rthratenlw, 1}, SD{rthratenlwc, 1});
  glw = o.glw; olr = o.olr; lwcf = o.lwcf; lwupt = o.lwupt; lwuptc = o.lwuptc; lwdnt = o.lwdnt;
  lwdntc = o.lwdntc; lwupb = o.lwupb; lwupbc = o.lwupbc; lwdnb = o.lwdnb; lwdnbc = o.lwdnbc;
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
