// SPDX-License-Identifier: GPL-3.0-or-later
#include "phys/ra_rrtmg_common/rrtmg_real.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw.hpp"

#include <algorithm>
#include <cmath>

namespace cppwrf::phys {

static constexpr rr PPROF[60] = {
    rr(1000.00), rr(855.47), rr(731.82), rr(626.05), rr(535.57), rr(458.16), rr(391.94), rr(335.29), rr(286.83), rr(245.38),
    rr(209.91), rr(179.57), rr(153.62), rr(131.41), rr(112.42), rr(96.17),  rr(82.27),  rr(70.38),  rr(60.21),  rr(51.51),
    rr(44.06),   rr(37.69),  rr(32.25),  rr(27.59),  rr(23.60),  rr(20.19),  rr(17.27),  rr(14.77),  rr(12.64),  rr(10.81),
    rr(9.25),    rr(7.91),   rr(6.77),   rr(5.79),   rr(4.95),   rr(4.24),   rr(3.63),   rr(3.10),   rr(2.65),   rr(2.27),
    rr(1.94),    rr(1.66),   rr(1.42),   rr(1.22),   rr(1.04),   rr(0.89),   rr(0.76),   rr(0.65),   rr(0.56),   rr(0.48),
    rr(0.41),    rr(0.35),   rr(0.30),   rr(0.26),   rr(0.22),   rr(0.19),   rr(0.16),   rr(0.14),   rr(0.12),   rr(0.10)};

static constexpr rr TPROF[60] = {
    rr(286.96), rr(281.07), rr(275.16), rr(268.11), rr(260.56), rr(253.02), rr(245.62), rr(238.41), rr(231.57), rr(225.91),
    rr(221.72), rr(217.79), rr(215.06), rr(212.74), rr(210.25), rr(210.16), rr(210.69), rr(212.14), rr(213.74), rr(215.37),
    rr(216.82), rr(217.94), rr(219.03), rr(220.18), rr(221.37), rr(222.64), rr(224.16), rr(225.88), rr(227.63), rr(229.51),
    rr(231.50), rr(233.73), rr(236.18), rr(238.78), rr(241.60), rr(244.44), rr(247.35), rr(250.33), rr(253.32), rr(256.30),
    rr(259.22), rr(262.12), rr(264.80), rr(266.50), rr(267.59), rr(268.44), rr(268.69), rr(267.76), rr(266.13), rr(263.96),
    rr(261.54), rr(258.93), rr(256.15), rr(253.23), rr(249.89), rr(246.67), rr(243.48), rr(240.25), rr(236.66), rr(233.86)};

static constexpr rr retab[95] = {
    rr(5.92779), rr(6.26422), rr(6.61973), rr(6.99539), rr(7.39234), rr(7.81177), rr(8.25496), rr(8.72323), rr(9.21800), rr(9.74075),
    rr(10.2930), rr(10.8765), rr(11.4929), rr(12.1440), rr(12.8317), rr(13.5581), rr(14.2319), rr(15.0351), rr(15.8799), rr(16.7674),
    rr(17.6986), rr(18.6744), rr(19.6955), rr(20.7623), rr(21.8757), rr(23.0364), rr(24.2452), rr(25.5034), rr(26.8125), rr(27.7895),
    rr(28.6450), rr(29.4167), rr(30.1088), rr(30.7306), rr(31.2943), rr(31.8151), rr(32.3077), rr(32.7870), rr(33.2657), rr(33.7540),
    rr(34.2601), rr(34.7892), rr(35.3442), rr(35.9255), rr(36.5316), rr(37.1602), rr(37.8078), rr(38.4720), rr(39.1508), rr(39.8442),
    rr(40.5552), rr(41.2912), rr(42.0635), rr(42.8876), rr(43.7863), rr(44.7853), rr(45.9170), rr(47.2165), rr(48.7221), rr(50.4710),
    rr(52.4980), rr(54.8315), rr(57.4898), rr(60.4785), rr(63.7898), rr(65.5604), rr(71.2885), rr(75.4113), rr(79.7368), rr(84.2351),
    rr(88.8833), rr(93.6658), rr(98.5739), rr(103.603), rr(108.752), rr(114.025), rr(119.424), rr(124.954), rr(130.630), rr(136.457),
    rr(142.446), rr(148.608), rr(154.956), rr(161.503), rr(168.262), rr(175.248), rr(182.473), rr(189.952), rr(197.699), rr(205.728),
    rr(214.055), rr(222.694), rr(231.661), rr(240.971), rr(250.639)};

namespace {
// inirad (module_ra_rrtmg_lw.F) climatology
static constexpr rr o3sum[31] = {
    rr(5.297e-8), rr(5.852e-8), rr(6.579e-8), rr(7.505e-8), rr(8.577e-8), rr(9.895e-8), rr(1.175e-7), rr(1.399e-7),
    rr(1.677e-7), rr(2.003e-7), rr(2.571e-7), rr(3.325e-7), rr(4.438e-7), rr(6.255e-7), rr(8.168e-7), rr(1.036e-6),
    rr(1.366e-6), rr(1.855e-6), rr(2.514e-6), rr(3.240e-6), rr(4.033e-6), rr(4.854e-6), rr(5.517e-6), rr(6.089e-6),
    rr(6.689e-6), rr(1.106e-5), rr(1.462e-5), rr(1.321e-5), rr(9.856e-6), rr(5.960e-6), rr(5.960e-6)};
static constexpr rr ppsum[31] = {
    rr(955.890), rr(850.532), rr(754.599), rr(667.742), rr(589.841), rr(519.421), rr(455.480), rr(398.085),
    rr(347.171), rr(301.735), rr(261.310), rr(225.360), rr(193.419), rr(165.490), rr(141.032), rr(120.125),
    rr(102.689), rr(87.829),  rr(75.123),  rr(64.306),  rr(55.086),  rr(47.209),  rr(40.535),  rr(34.795),
    rr(29.865),  rr(19.122),  rr(9.277),   rr(4.660),   rr(2.421),   rr(1.294),   rr(0.647)};
static constexpr rr o3win[31] = {
    rr(4.629e-8), rr(4.686e-8), rr(5.017e-8), rr(5.613e-8), rr(6.871e-8), rr(8.751e-8), rr(1.138e-7), rr(1.516e-7),
    rr(2.161e-7), rr(3.264e-7), rr(4.968e-7), rr(7.338e-7), rr(1.017e-6), rr(1.308e-6), rr(1.625e-6), rr(2.011e-6),
    rr(2.516e-6), rr(3.130e-6), rr(3.840e-6), rr(4.703e-6), rr(5.486e-6), rr(6.289e-6), rr(6.993e-6), rr(7.494e-6),
    rr(8.197e-6), rr(9.632e-6), rr(1.113e-5), rr(1.146e-5), rr(9.389e-6), rr(6.135e-6), rr(6.135e-6)};
static constexpr rr ppwin[31] = {
    rr(955.747), rr(841.783), rr(740.199), rr(649.538), rr(568.404), rr(495.815), rr(431.069), rr(373.464),
    rr(322.354), rr(277.190), rr(237.635), rr(203.433), rr(174.070), rr(148.949), rr(127.408), rr(108.915),
    rr(93.114),  rr(79.551),  rr(67.940),  rr(58.072),  rr(49.593),  rr(42.318),  rr(36.138),  rr(30.907),
    rr(26.362),  rr(16.423),  rr(7.583),   rr(3.620),   rr(1.807),   rr(0.938),   rr(0.469)};


}  // namespace

void rrtmg_lw_fill_consts(RrtmgLwConsts& k) {
  for (int i = 0; i < 60; ++i) { k.pprof[i] = PPROF[i]; k.tprof[i] = TPROF[i]; }
  for (int i = 0; i < 95; ++i) k.retab[i] = retab[i];
  // inirad: annual mean ozone (same arithmetic as WRF, once instead of every column)
  rr* o3ann = k.o3ann;
  o3ann[0] = rr(0.5) * (o3sum[0] + o3win[0]);
  for (int kk = 1; kk < 31; ++kk) {
    o3ann[kk] = o3win[kk - 1] + (o3win[kk] - o3win[kk - 1]) / (ppwin[kk] - ppwin[kk - 1]) * (ppsum[kk] - ppwin[kk - 1]);
  }
  for (int kk = 1; kk < 31; ++kk) {
    o3ann[kk] = rr(0.5) * (o3ann[kk] + o3sum[kk]);
  }
  rr* ppwrkh = k.ppwrkh;
  ppwrkh[0] = rr(1100.0);
  for (int kk = 1; kk < 31; ++kk) {
    ppwrkh[kk] = rr(0.5) * (ppsum[kk] + ppsum[kk - 1]);
  }
  ppwrkh[31] = rr(0.0);
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
  c.co2 = (rr(280.) + rr(90.) * libm::exp(rr(0.02) * (c.yr - 2000))) * rr(1.e-6);
  c.ch4 = rr(1774.e-9);
  c.n2o = rr(319.e-9);
  c.cfc11 = rr(0.251e-9);
  c.cfc12 = rr(0.538e-9);
}

RrtmgLwCallConfig rrtmg_lw_call_config(int num_levels, const RrtmgLwConfig& config, bool hrc, bool hri, bool hrs) {
  constexpr rr deltap = rr(4.0);  // mb
  constexpr rr grav = rr(9.8066);
  RrtmgLwCallConfig c{};
  c.num_levels = num_levels;
  // rrtmg_lwinit: NLAYERS = kme + nint(p_top*0.01/deltap) - 1 with kme = nz + 1 (Fortran NINT)
  c.nlayers = num_levels + int(std::lround(config.p_top * rr(0.01) / deltap));
  c.o3input = config.o3input;
  c.co2 = config.co2; c.n2o = config.n2o; c.ch4 = config.ch4; c.cfc11 = config.cfc11; c.cfc12 = config.cfc12;
  c.g = config.g;
  c.heatfac = grav * rr(86400.0) / (config.cp * rr(100.0));
  c.has_reqc = hrc; c.has_reqi = hri; c.has_reqs = hrs;
  c.icld = config.cldovrlp; c.idcor = config.idcor;
  c.juldat = int(config.julian);  // juldat = julian (INTEGER)
  return c;
}


// WRF rrtmg_lw (inatm, cldprmc, setcoef, taumol, rtrnmc) for one column whose McICA subcolumn
// cloud arrays are already in ws (cldfmc, taucmc). Layer arrays have
// nlayers entries, level arrays nlayers + 1 (index 0 = surface); pressures in mb, emis per band.
void rrtmg_lw_core(int nlayers, rr cp, const rr* play, const rr* plev, const rr* tlay,
                   const rr* tlev, rr tsk, const rr* h2ovmr, const rr* o3vmr, const rr* co2vmr,
                   const rr* ch4vmr, const rr* n2ovmr, const rr* o2vmr, const rr* cfc11vmr,
                   const rr* cfc12vmr, const rr* cfc22vmr, const rr* ccl4vmr, const rr* semiss,
                   const rr* rei, const rr* rel, const rr* res, rr* totuflux, rr* totdflux,
                   rr* htr, rr* totuclfl, rr* totdclfl, rr* htrc, RrtmgLwWorkspace& ws,
                   const RrtmgLwTables& tables, int inflag, int iceflag, int liqflag) {
  constexpr rr grav = rr(9.8066);
  LwLayout lay(nlayers);
  ensure_workspace(ws, lay);
  const LwViews v = lay.views(ws.d.data(), ws.i.data(), 0, 1);
  v.pz[0] = plev[0];
  v.tz[0] = tlev[0];
  for (int l = 0; l < nlayers; ++l) {
    v.pavel[l] = play[l];
    v.tavel[l] = tlay[l];
    v.pz[l + 1] = plev[l + 1];
    v.tz[l + 1] = tlev[l + 1];
    v.wkl[0][l] = h2ovmr[l]; v.wkl[1][l] = co2vmr[l]; v.wkl[2][l] = o3vmr[l]; v.wkl[3][l] = n2ovmr[l];
    v.wkl[4][l] = rr(0.0); v.wkl[5][l] = ch4vmr[l]; v.wkl[6][l] = o2vmr[l];
    v.wx[0][l] = ccl4vmr[l]; v.wx[1][l] = cfc11vmr[l]; v.wx[2][l] = cfc12vmr[l]; v.wx[3][l] = cfc22vmr[l];
    v.reicmc[l] = rei[l]; v.relqmc[l] = rel[l]; v.resnmc[l] = res[l];
    // McICA flags from ws; the subcolumn water paths are not used (the core check has no clouds)
    for (int g = 0; g < 140; ++g) { v.cldfmc[g][l] = ws.cldfmc[g][l]; v.taucmc[g][l] = ws.taucmc[g][l]; }
    v.ciwpth[l] = rr(0.0); v.clwpth[l] = rr(0.0); v.cswpth[l] = rr(0.0);
  }
  for (int nb = 0; nb < 16; ++nb) v.semiss[nb] = semiss[nb];
  v.iscal[1] = inflag; v.iscal[2] = iceflag; v.iscal[3] = liqflag;
  rrtmg_lw_icldlyr(v.cldfmc, v.icldlyr, 0, nlayers);
  int any_cld = 0;
  for (int l = 0; l < nlayers; ++l) if (v.icldlyr[l]) any_cld = 1;
  v.iscal[4] = any_cld;
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
  rrtmg_lw_heating(nlayers, v.pz, v.fnet, v.fnetc, v.htr, v.htrc, grav * rr(86400.0) / (cp * rr(100.0)));
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
    rr tsk,
    rr emiss,
    rr xland,
    rr xice,
    rr snow,
    rr xlat,
    const rr* p8w,
    const rr* t8w,
    const rr* p3d,
    const rr* pi3d,
    const rr* dz8w,
    const rr* t3d,
    const rr* rho3d,
    const rr* qv3d,
    const rr* qc3d,
    const rr* qr3d,
    const rr* qi3d,
    const rr* qs3d,
    const rr* qg3d,
    const rr* cldfra3d,
    const rr* o33d,
    rr& glw,
    rr& olr,
    rr& lwcf,
    rr& lwupt,
    rr& lwuptc,
    rr& lwdnt,
    rr& lwdntc,
    rr& lwupb,
    rr& lwupbc,
    rr& lwdnb,
    rr& lwdnbc,
    rr* rthratenlw,
    rr* rthratenlwc,
    RrtmgLwWorkspace& ws,
    const rr* re_cloud,
    const rr* re_ice,
    const rr* re_snow,
    rr* o3_out
) {
  (void)qr3d; (void)qg3d; (void)xlat; (void)rho3d; (void)dz8w;
  const RrtmgLwCallConfig c = rrtmg_lw_call_config(num_levels, config, re_cloud, re_ice, re_snow);
  LwLayout lay(c.nlayers);
  ensure_workspace(ws, lay);
  const LwViews v = lay.views(ws.d.data(), ws.i.data(), 0, 1);
  LwColumnIn in{};
  in.tsk = tsk; in.emiss = emiss; in.xland = xland; in.xice = xice; in.snow = snow; in.xlat = xlat;
  in.p8w = sd_const(p8w); in.t8w = sd_const(t8w); in.p3d = sd_const(p3d); in.pi3d = sd_const(pi3d); in.dz8w = sd_const(dz8w); in.t3d = sd_const(t3d);
  in.rho3d = sd_const(rho3d); in.qv3d = sd_const(qv3d); in.qc3d = sd_const(qc3d); in.qi3d = sd_const(qi3d); in.qs3d = sd_const(qs3d);
  in.cldfra3d = sd_const(cldfra3d); in.o33d = sd_const(o33d ? o33d : qv3d);
  in.re_cloud = sd_const(re_cloud); in.re_ice = sd_const(re_ice); in.re_snow = sd_const(re_snow);
  in.has_qc = in.has_qi = in.has_qs = true;  // absent species are passed as zero arrays
  in.o3_out = SD{o3_out, 1};
  rrtmg_lw_prep(c, in, v, tables, rrtmg_lw_consts());
  for (int l = 0; l < c.nlayers; ++l) rrtmg_lw_layer(c.nlayers, l, v, tables);
  for (int b = 0; b < 16; ++b) rrtmg_lw_bandstage(c.nlayers, b, v, tables);
  for (int lev = 0; lev <= c.nlayers; ++lev) rrtmg_lw_level(lev, v, tables);
  const LwColumnOut o = rrtmg_lw_finish(c, v, sd_const(pi3d), SD{rthratenlw, 1}, SD{rthratenlwc, 1});
  glw = o.glw; olr = o.olr; lwcf = o.lwcf; lwupt = o.lwupt; lwuptc = o.lwuptc; lwdnt = o.lwdnt;
  lwdntc = o.lwdntc; lwupb = o.lwupb; lwupbc = o.lwupbc; lwdnb = o.lwdnb; lwdnbc = o.lwdnbc;
}

} // namespace cppwrf::phys
