#pragma once
// RRTMG_LWRAD (WRF module_ra_rrtmg_lw.F) for one column as independent stages, every per-column
// array a strided view (src/phys/ra_rrtmg_common/strided.hpp):
//   rrtmg_lw_prep     per column : RRTMG_LWRAD profiles, cloud properties, McICA subcolumns
//                                  (mcica_subcol_lw, max-random overlap), inatm, setcoef, secdiff
//   rrtmg_lw_layer    per layer  : taumol, cldprmc, cloudy-layer flag
//   rrtmg_lw_bandstage per band  : rtrnmc g-point loop of the band (rrtmg_lw_rtrnmc.hpp)
//   rrtmg_lw_level    per level  : rtrnmc band sums
//   rrtmg_lw_finish   per column : heating rates and the RRTMG_LWRAD outputs
// The CPU column driver (rrtmg_lw_column) runs them in this order on stride-1 arrays; the GPU
// batch driver runs each stage for all columns (layers, bands, levels) at once. The arithmetic
// of every value is WRF's (checked bit for bit against WRF dumps on the CPU path).
#include <algorithm>
#include <cmath>
#include <cstddef>

#include "core/backend.hpp"
#include "phys/ra_rrtmg_common/kissvec.hpp"
#include "phys/ra_rrtmg_common/strided.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_cldprmc.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_rtrnmc.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_setcoef.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_tables.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_taumol.hpp"

namespace cppwrf::phys {

// Constants of RRTMG_LWRAD / inirad (module data in WRF), device-copyable.
struct RrtmgLwConsts {
  double pprof[60], tprof[60];  // standard profile above the model top
  double retab[95];             // ice effective radius table (reicalc)
  double o3ann[31], ppwrkh[32]; // inirad annual ozone and its pressure edges
};
void rrtmg_lw_fill_consts(RrtmgLwConsts& k);  // rrtmg_lw.cpp

// Configuration of one radiation call (all columns).
struct RrtmgLwCallConfig {
  int num_levels, nlayers;
  int o3input;
  double co2, n2o, ch4, cfc11, cfc12;
  double g;                    // gravmks of the water paths
  double heatfac;              // rrtmg_lw heatfac = grav*86400/(cp*100)
  int has_reqc, has_reqi, has_reqs;
  int icld, idcor, juldat;     // McICA overlap (cldovrlp), decorrelation length type, day of year
};

// Per-column views. nl = layer dimension (>= nlayers), all level arrays nl + 1.
struct LwViews {
  // inputs of the stages after prep / outputs of prep (layer index 0 = surface layer)
  SD pz, tz, pavel, tavel, coldry, wbrodl, semiss, plankbnd, secdiff, reicmc, relqmc, resnmc;
  SD cldfrac, clwpth, ciwpth, cswpth, o3mmr;
  SD2 wkl, wx, planklay, planklev;
  SI jp, jt, jt1, indself, indfor, indminor, icldlyr;
  SD colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11;
  SD rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1;
  SD rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, selffac, selffrac, forfac, forfrac, minorfrac, scaleminor,
      scaleminorn2;
  SD2 cldfmc, taucmc;                          // [140][nl]
  SD2 fracs, taug;                             // [nl][140]
  SD2 bu, bd, bcu, bcd;                        // [16][nl + 1] band sums of rtrnmc
  SD totuflux, totdflux, totuclfl, totdclfl, fnet, fnetc, htr, htrc;  // [nl + 1]
  SI iscal;  // [0] laytrop, [1] inflag, [2] iceflag, [3] liqflag, [4] any_cloud
  SD dscal;  // [0] pwvcm
};

// Storage of the views: element e of column col of a batch of ncol columns at base[e*ncol + col]
// (doubles and ints separately). nl = layer dimension (>= nlayers).
struct LwLayout {
  int nl;
  std::size_t nd = 0, ni = 0;  // doubles and ints per column
  CPPWRF_HD explicit LwLayout(int nl_) : nl(nl_) { views(nullptr, nullptr, 0, 1, &nd, &ni); }
  CPPWRF_FORCEINLINE LwViews views(double* d, int* ip, int col, int ncol, std::size_t* nd_out = nullptr,
                                   std::size_t* ni_out = nullptr) const {
    LwViews v{};
    const int L = nl, L1 = nl + 1, s = ncol;
    std::size_t od = 0, oi = 0;
    auto take = [&](std::size_t n) {
      double* q = d ? d + od * std::size_t(s) + col : nullptr;
      od += n;
      return q;
    };
    auto v1 = [&](int n) { return SD{take(n), s}; };
    auto v2 = [&](int a, int b) { return SD2{take(std::size_t(a) * b), s, b}; };
    auto i1 = [&](int n) {
      int* q = ip ? ip + oi * std::size_t(s) + col : nullptr;
      oi += n;
      return SI{q, s};
    };
    v.pz = v1(L1); v.tz = v1(L1); v.pavel = v1(L1); v.tavel = v1(L1); v.coldry = v1(L1); v.wbrodl = v1(L1);
    v.semiss = v1(16); v.plankbnd = v1(16); v.secdiff = v1(16);
    v.reicmc = v1(L1); v.relqmc = v1(L1); v.resnmc = v1(L1);
    v.cldfrac = v1(L1); v.clwpth = v1(L1); v.ciwpth = v1(L1); v.cswpth = v1(L1); v.o3mmr = v1(L1);
    v.wkl = v2(7, L); v.wx = v2(4, L); v.planklay = v2(L, 16); v.planklev = v2(L1, 16);
    v.colh2o = v1(L); v.colco2 = v1(L); v.colo3 = v1(L); v.coln2o = v1(L); v.colco = v1(L); v.colch4 = v1(L);
    v.colo2 = v1(L); v.colbrd = v1(L);
    v.fac00 = v1(L); v.fac01 = v1(L); v.fac10 = v1(L); v.fac11 = v1(L);
    v.rat_h2oco2 = v1(L); v.rat_h2oco2_1 = v1(L); v.rat_h2oo3 = v1(L); v.rat_h2oo3_1 = v1(L);
    v.rat_h2on2o = v1(L); v.rat_h2on2o_1 = v1(L); v.rat_h2och4 = v1(L); v.rat_h2och4_1 = v1(L);
    v.rat_n2oco2 = v1(L); v.rat_n2oco2_1 = v1(L); v.rat_o3co2 = v1(L); v.rat_o3co2_1 = v1(L);
    v.selffac = v1(L); v.selffrac = v1(L); v.forfac = v1(L); v.forfrac = v1(L); v.minorfrac = v1(L);
    v.scaleminor = v1(L); v.scaleminorn2 = v1(L);
    v.cldfmc = v2(140, L); v.taucmc = v2(140, L);
    v.fracs = v2(L, 140); v.taug = v2(L, 140);
    v.bu = v2(16, L1); v.bd = v2(16, L1); v.bcu = v2(16, L1); v.bcd = v2(16, L1);
    v.totuflux = v1(L1); v.totdflux = v1(L1); v.totuclfl = v1(L1); v.totdclfl = v1(L1);
    v.fnet = v1(L1); v.fnetc = v1(L1); v.htr = v1(L1); v.htrc = v1(L1);
    v.dscal = v1(1);
    v.jp = i1(L); v.jt = i1(L); v.jt1 = i1(L); v.indself = i1(L); v.indfor = i1(L); v.indminor = i1(L);
    v.icldlyr = i1(L); v.iscal = i1(5);
    if (nd_out) *nd_out = od;
    if (ni_out) *ni_out = oi;
    return v;
  }
};

// Per-column inputs of RRTMG_LWRAD (k = 0 at the surface), model memory views.
struct LwColumnIn {
  double tsk, emiss, xland, xice, snow, xlat;
  SD p8w, t8w, p3d, pi3d, dz8w, t3d, rho3d, qv3d, qc3d, qi3d, qs3d, cldfra3d, o33d;
  SD re_cloud, re_ice, re_snow;  // p == nullptr: has_req* = 0
  bool has_qc, has_qi, has_qs;   // F_QC, F_QI, F_QS (absent species are zero)
  SD o3_out;  // o3input != 2, p != nullptr: receives the ozone vmr of the model levels (WRF's O33D = O31D)
};

// inatm of rrtmg_lw: wkl holds the volume mixing ratios on entry, the column amounts on exit.
CPPWRF_HD inline void rrtmg_lw_inatm(int nlayers, const LwViews& v) {
  constexpr double amd = 28.9660;
  constexpr double amw = 18.0160;
  constexpr double avogad = 6.02214199e23;
  constexpr double grav = 9.8066;
  for (int l = 0; l < nlayers; ++l) {
    double amm = (1.0 - v.wkl[0][l]) * amd + v.wkl[0][l] * amw;
    v.coldry[l] = (v.pz[l] - v.pz[l + 1]) * 1.0e3 * avogad / (1.0e2 * grav * amm * (1.0 + v.wkl[0][l]));
  }
  double amttl = 0.0;
  double wvttl = 0.0;
  for (int l = 0; l < nlayers; ++l) {
    double summol = 0.0;
    for (int imol = 1; imol < 7; ++imol) summol += v.wkl[imol][l];
    v.wbrodl[l] = v.coldry[l] * (1.0 - summol);
    for (int imol = 0; imol < 7; ++imol) v.wkl[imol][l] = v.coldry[l] * v.wkl[imol][l];
    amttl = amttl + v.coldry[l] + v.wkl[0][l];  // WRF order: (amttl + coldry) + wkl
    wvttl += v.wkl[0][l];
    for (int ix = 0; ix < 4; ++ix) v.wx[ix][l] = v.coldry[l] * v.wx[ix][l] * 1.0e-20;
  }
  double wvsh = (amw * wvttl) / (amd * amttl);
  v.dscal[0] = wvsh * (1.0e3 * v.pz[0]) / (1.0e2 * grav);  // pwvcm
}

// ---- stage 1: per column ------------------------------------------------------------------------
CPPWRF_FORCEINLINE void rrtmg_lw_prep(const RrtmgLwCallConfig& c, const LwColumnIn& in, const LwViews& v,
                                      const RrtmgLwTables& CPPWRF_RESTRICT tables,
                                      const RrtmgLwConsts& CPPWRF_RESTRICT kc) {
  const int num_levels = c.num_levels, nlayers = c.nlayers;
  constexpr double deltap = 4.0;  // mb
  const double xland = in.xland;
  auto qc_of = [&](int k) { return in.has_qc ? in.qc3d[k] : 0.0; };
  auto qi_of = [&](int k) { return in.has_qi ? in.qi3d[k] : 0.0; };
  auto qs_of = [&](int k) { return in.has_qs ? in.qs3d[k] : 0.0; };

  // cloud optics flags (RRTMG_LWRAD, icloud != 0)
  int inflglw = 2, iceflglw = 3;
  const int liqflglw = 1;
  const bool hrc = in.re_cloud.p != nullptr, hri = in.re_ice.p != nullptr, hrs = in.re_snow.p != nullptr;
  if (hrc) inflglw = 3;
  if (hri) { inflglw = 4; iceflglw = 4; }
  if (hrs) { inflglw = 5; iceflglw = 5; }
  const bool p3case = !hrs && hri && hrc;  // P3: ice into the snow category for the optics
  if (p3case) { inflglw = 5; iceflglw = 5; }

  // CAM relcalc inputs
  const double landfrac = 2.0 - xland;
  const double landm = landfrac;
  const double snowh = 0.001 * in.snow;
  const double icefrac = in.xice;
  constexpr double rliqland = 8.0, rliqocean = 14.0, rliqice = 14.0, tmelt = 273.16;
  // relcalc / reicalc at the layer temperature t (RRTMG_LWRAD passes them tlay)
  auto relcalc = [&](double t) {
    double reliq = rliqland + (rliqocean - rliqland) * std::min(1.0, std::max(0.0, (tmelt - t) * 0.05));
    reliq += (rliqocean - reliq) * std::min(1.0, std::max(0.0, snowh * 10.0));
    reliq += (rliqocean - reliq) * std::min(1.0, std::max(0.0, 1.0 - landm));
    reliq += (rliqice - reliq) * std::min(1.0, std::max(0.0, icefrac));
    return reliq;
  };
  auto reicalc = [&](double t) {
    int idx = static_cast<int>(t - 179.0);
    idx = std::min(std::max(idx, 1), 94);
    const double corr = t - static_cast<int>(t);
    return kc.retab[idx - 1] * (1.0 - corr) + kc.retab[idx] * corr;
  };

  // base layer vertical variables, gases, cloud physical properties
  v.pz[0] = in.p8w[0] / 100.0;
  v.tz[0] = in.t8w[0];
  for (int k = 0; k < num_levels; ++k) {
    const double t3 = in.t3d[k];
    v.pavel[k] = in.p3d[k] / 100.0;
    v.pz[k + 1] = in.p8w[k + 1] / 100.0;
    const double pdel = v.pz[k] - v.pz[k + 1];
    v.tavel[k] = t3;
    v.tz[k + 1] = in.t8w[k + 1];

    double qv_val = std::max(in.qv3d[k], 1.0e-12);
    v.wkl[0][k] = qv_val * 1.607793;  // h2ovmr
    v.wkl[1][k] = c.co2;
    v.wkl[2][k] = 0.0;                // o3vmr below
    v.wkl[3][k] = c.n2o;
    v.wkl[4][k] = 0.0;                // co
    v.wkl[5][k] = c.ch4;
    v.wkl[6][k] = 0.209488;           // o2
    v.wx[0][k] = 0.093e-9;            // ccl4
    v.wx[1][k] = c.cfc11;
    v.wx[2][k] = c.cfc12;
    v.wx[3][k] = 0.169e-9;            // cfc22

    // hydrometeors (max(0, q)) and effective radii
    double qc1 = std::max(0.0, qc_of(k)), qi1 = std::max(0.0, qi_of(k)), qs1 = std::max(0.0, qs_of(k));
    const double cf3 = in.cldfra3d[k];
    double recloud1d = 5.0, reice1d = 10.0, resnow1d = 10.0;
    if (hrc) {
      double r = std::max(2.5, in.re_cloud[k] * 1.E6);
      if (r <= 2.5 && cf3 > 0. && (xland - 1.5) > 0.) r = 10.5;       // ocean
      else if (r <= 2.5 && cf3 > 0. && (xland - 1.5) < 0.) r = 7.5;   // land
      recloud1d = r;
    }
    if (hri) {
      double r = std::max(5., in.re_ice[k] * 1.E6);
      if (r <= 5. && cf3 > 0.) {
        int idx_rei = int(t3 - 179.);
        idx_rei = std::min(std::max(idx_rei, 1), 75);
        const double corr = t3 - int(t3);
        r = kc.retab[idx_rei - 1] * (1. - corr) + kc.retab[idx_rei] * corr;
        r = std::max(r, 5.0);
      }
      reice1d = r;
    }
    if (hrs) resnow1d = std::max(10., in.re_snow[k] * 1.E6);
    if (p3case) {
      resnow1d = std::max(10., in.re_ice[k] * 1.E6);
      qs1 = qi_of(k);
      qi1 = 0.;
      reice1d = 10.;
    }
    // relcalc / reicalc (the top model layer again below, after the buffer changed its temperature)
    double reliq = inflglw >= 3 ? recloud1d : relcalc(t3);
    double reice = iceflglw >= 4 ? reice1d : reicalc(t3);
    if (iceflglw == 3) {  // Fu (1996) generalized effective size
      reice = reice * 1.0315;
      reice = std::min(140.0, reice);
    }
    // grid box and in-cloud water paths (g/m2), CAM method
    double gicewp = (qi1 + qs1) * pdel * 100.0 / c.g * 1000.0;
    const double gliqwp = qc1 * pdel * 100.0 / c.g * 1000.0;
    double cicewp = gicewp / std::max(0.01, cf3);
    const double cliqwp = gliqwp / std::max(0.01, cf3);
    if (iceflglw >= 4) {
      gicewp = qi1 * pdel * 100.0 / c.g * 1000.0;
      cicewp = gicewp / std::max(0.01, cf3);
    }
    double csnowp = 0.0;
    if (iceflglw == 5) {
      double snow_mass_factor = 0.99;
      if (resnow1d > 130.) {
        snow_mass_factor = std::min(snow_mass_factor, (130.0 / resnow1d) * (130.0 / resnow1d));
        resnow1d = 130.0;
      }
      const double gsnowp = qs1 * snow_mass_factor * pdel * 100.0 / c.g * 1000.0;
      csnowp = gsnowp / std::max(0.01, cf3);
    }
    v.clwpth[k] = cliqwp;
    v.ciwpth[k] = cicewp;
    v.relqmc[k] = reliq;
    v.reicmc[k] = reice;
    if (inflglw == 5) {
      v.cswpth[k] = csnowp;
      v.resnmc[k] = resnow1d;
    } else {
      v.cswpth[k] = 0.;
      v.resnmc[k] = 10.;
    }
    v.cldfrac[k] = cf3;
  }

  // buffer layers above the model top (to 0 mb), Cavallo (2010)
  for (int L = num_levels; L < nlayers; ++L) {
    v.pz[L + 1] = v.pz[L] - deltap;
    v.pavel[L] = 0.5 * (v.pz[L] + v.pz[L + 1]);
  }
  v.pz[nlayers] = 0.00;
  v.pavel[nlayers - 1] = 0.5 * (v.pz[nlayers - 1] + v.pz[nlayers]);
  // standard profile temperatures at the column pressure levels, matched to the model top
  auto varint = [&](int L) {
    const double p_val = v.pz[L];
    int klev = 59;
    if (kc.pprof[59] < p_val) {
      for (int LL = 1; LL < 60; ++LL) {
        if (kc.pprof[LL] < p_val) {
          klev = LL - 1;
          break;
        }
      }
    }
    double wght = 0.0;
    double vark = kc.tprof[klev];
    double vark1 = (klev < 59) ? kc.tprof[klev + 1] : kc.tprof[klev];
    if (klev < 59) wght = (p_val - kc.pprof[klev]) / (kc.pprof[klev + 1] - kc.pprof[klev]);
    return wght * (vark1 - vark) + vark;
  };
  const double vtop = varint(num_levels - 1);
  for (int L = num_levels; L <= nlayers; ++L) {
    v.tz[L] = varint(L) + (v.tz[num_levels - 1] - vtop);
    v.tavel[L - 1] = 0.5 * (v.tz[L] + v.tz[L - 1]);
  }
  // RRTMG_LWRAD calls relcalc / reicalc after this adjustment: the top model layer's sizes use its new tlay
  {
    const int k = num_levels - 1;
    if (inflglw < 3) v.relqmc[k] = relcalc(v.tavel[k]);
    if (iceflglw == 3) v.reicmc[k] = std::min(140.0, reicalc(v.tavel[k]) * 1.0315);
  }
  for (int L = num_levels; L < nlayers; ++L) {
    for (int m = 0; m < 7; ++m) v.wkl[m][L] = v.wkl[m][num_levels - 1];
    for (int m = 0; m < 4; ++m) v.wx[m][L] = v.wx[m][num_levels - 1];
    v.clwpth[L] = 0.0;
    v.ciwpth[L] = 0.0;
    v.cswpth[L] = 0.0;
    v.relqmc[L] = 10.0;
    v.reicmc[L] = 10.0;
    v.resnmc[L] = 10.0;
    v.cldfrac[L] = 0.0;
  }

  // ozone: inirad climatology, o3input = 2 CAM ozone on the model levels, shifted above
  for (int k = 0; k <= nlayers - 1; ++k) {
    double o3 = 0.0;
    const double pr_k = v.pz[k], pr_k1 = v.pz[k + 1];
    for (int jj = 0; jj < 31; ++jj) {
      double pb1 = (pr_k > kc.ppwrkh[jj]) ? (pr_k - kc.ppwrkh[jj]) : 0.0;
      double pb2 = (pr_k > kc.ppwrkh[jj + 1]) ? (pr_k - kc.ppwrkh[jj + 1]) : 0.0;
      double pt1 = (pr_k1 > kc.ppwrkh[jj]) ? (pr_k1 - kc.ppwrkh[jj]) : 0.0;
      double pt2 = (pr_k1 > kc.ppwrkh[jj + 1]) ? (pr_k1 - kc.ppwrkh[jj + 1]) : 0.0;
      o3 += (pb2 - pb1 - pt2 + pt1) * kc.o3ann[jj];
    }
    v.o3mmr[k] = o3 / (pr_k - pr_k1);
  }
  constexpr double amdo = 0.603461;
  for (int k = 0; k < nlayers; ++k) {
    double o3vmr = v.o3mmr[k] * amdo;
    if (c.o3input == 2) {
      if (k < num_levels) {
        o3vmr = in.o33d[k];
      } else {
        o3vmr = in.o33d[num_levels - 1] - v.o3mmr[num_levels - 1] * amdo + v.o3mmr[k] * amdo;
        if (o3vmr <= 0.) o3vmr = v.o3mmr[k] * amdo;
      }
    } else if (in.o3_out.p && k < num_levels) {
      in.o3_out[k] = o3vmr;  // o31d(k) = o3vmr(ncol,k), then O33D(I,K,J) = O31D(K) ("output o3 for o3input=0")
    }
    v.wkl[2][k] = o3vmr;
  }

  for (int nb = 0; nb < 16; ++nb) v.semiss[nb] = in.emiss;

  // McICA subcolumns (mcica_subcol_lw / generate_stochastic_clouds, KISS, changeSeed = 150, overlap icld)
  bool any_cloud = false;
  for (int ilev = 0; ilev < nlayers; ++ilev) {
    if (CSD(v.cldfrac)[ilev] >= 1.0e-20) any_cloud = true;
  }
#if defined(__CUDA_ARCH__)
  const bool run_mcica = __any_sync(__activemask(), any_cloud);
#else
  const bool run_mcica = any_cloud;
#endif
  std::uint64_t icld_mask_lo = 0, icld_mask_hi = 0;
  if (run_mcica && any_cloud) {
    KissRng rng;
    rng.init_from_pressure(CSD(v.pavel)[0] * 100.0, CSD(v.pavel)[1] * 100.0, CSD(v.pavel)[2] * 100.0,
                           CSD(v.pavel)[3] * 100.0);
    rng.permute(150);
    // layer mid-point heights from dz8w (exponential overlaps): above the model top the top layer's dz
    double hgt[kMcicaMaxLayers];
    double zo_inv = 0.0;
    if (c.icld == 4 || c.icld == 5) {
      double dzsum = 0.0, dz = 0.0;
      for (int k = 0; k < nlayers; ++k) {
        if (k < num_levels) dz = in.dz8w[k];
        hgt[k] = dzsum + 0.5 * dz;
        dzsum = dzsum + dz;
      }
      zo_inv = 1.0 / mcica_decorr_len(c.idcor, c.juldat, in.xlat);
    }
    auto cldf_of = [&](int l) { const double cf = CSD(v.cldfrac)[l]; return cf < 1.0e-20 ? 0.0 : cf; };
    auto alpha_of = [&](int l) {
      if (c.icld == 5 && cldf_of(l) == 0.0 && cldf_of(l - 1) > 0.0) return 0.0;
      return exp(-(hgt[l] - hgt[l - 1]) * zo_inv);
    };
    // cloudy subcolumns get the layer's in-cloud paths (cldprmc reads them through cldfmc)
    mcica_subcolumns(c.icld, rng, 140, nlayers, cldf_of, alpha_of, [&](int isubcol, int ilev, bool cld) {
      if (cldf_of(ilev) > 0.0) {
        v.cldfmc[isubcol][ilev] = cld ? 1.0 : 0.0;
        if (cld) {
          if (ilev < 64) icld_mask_lo |= (std::uint64_t(1) << ilev);
          else icld_mask_hi |= (std::uint64_t(1) << (ilev - 64));
        }
      }
    });
  }
  for (int ilev = 0; ilev < nlayers; ++ilev) {
    const int bit = (ilev < 64) ? int((icld_mask_lo >> ilev) & 1u) : int((icld_mask_hi >> (ilev - 64)) & 1u);
    v.icldlyr[ilev] = bit;
  }
  v.iscal[1] = inflglw;
  v.iscal[2] = iceflglw;
  v.iscal[3] = liqflglw;
  v.iscal[4] = (icld_mask_lo != 0 || icld_mask_hi != 0) ? 1 : 0;

  // inatm
  rrtmg_lw_inatm(nlayers, v);
  // setcoef
  int laytrop = 0;
  rrtmg_lw_setcoef(nlayers, 1, v.pavel, v.tavel, v.tz, in.tsk, v.semiss, v.coldry, v.wkl, v.wbrodl, laytrop, v.jp,
                   v.jt, v.jt1, v.planklay, v.planklev, v.plankbnd, v.colh2o, v.colco2, v.colo3, v.coln2o, v.colco,
                   v.colch4, v.colo2, v.colbrd, v.fac00, v.fac01, v.fac10, v.fac11, v.rat_h2oco2, v.rat_h2oco2_1,
                   v.rat_h2oo3, v.rat_h2oo3_1, v.rat_h2on2o, v.rat_h2on2o_1, v.rat_h2och4, v.rat_h2och4_1,
                   v.rat_n2oco2, v.rat_n2oco2_1, v.rat_o3co2, v.rat_o3co2_1, v.selffac, v.selffrac, v.indself,
                   v.forfac, v.forfrac, v.indfor, v.minorfrac, v.scaleminor, v.scaleminorn2, v.indminor, tables);
  v.iscal[0] = laytrop;
  rrtmg_lw_secdiff(v.dscal[0], v.secdiff);
}

// ---- stage 2: per layer -------------------------------------------------------------------------
CPPWRF_FORCEINLINE void rrtmg_lw_layer(int nlayers, int lay, const LwViews& v,
                                       const RrtmgLwTables& CPPWRF_RESTRICT tables) {
  const int laytrop = CSI(v.iscal)[0];
  if (CSI(v.icldlyr)[lay] != 0) {
    rrtmg_lw_cldprmc(nlayers, CSI(v.iscal)[1], CSI(v.iscal)[2], CSI(v.iscal)[3], v.cldfmc, v.ciwpth, v.clwpth,
                     v.cswpth, v.reicmc, v.relqmc, v.resnmc, v.taucmc, tables, lay, lay + 1);
  }
  rrtmg_lw_taumol(nlayers, v.pavel, v.wx, v.coldry, laytrop, v.jp, v.jt, v.jt1, v.colh2o, v.colco2, v.colo3,
                  v.coln2o, v.colco, v.colch4, v.colo2, v.colbrd, v.fac00, v.fac01, v.fac10, v.fac11, v.rat_h2oco2,
                  v.rat_h2oco2_1, v.rat_h2oo3, v.rat_h2oo3_1, v.rat_h2on2o, v.rat_h2on2o_1, v.rat_h2och4,
                  v.rat_h2och4_1, v.rat_n2oco2, v.rat_n2oco2_1, v.rat_o3co2, v.rat_o3co2_1, v.selffac, v.selffrac,
                  v.indself, v.forfac, v.forfrac, v.indfor, v.minorfrac, v.scaleminor, v.scaleminorn2, v.indminor,
                  v.fracs, v.taug, tables, lay, lay + 1);
}

// ---- stage 3: per band (taut = taug: no aerosols) ----------------------------------------------
template <int MAX_LAY = RRTMG_LW_SWEEP_MAX>
CPPWRF_FORCEINLINE void rrtmg_lw_bandstage(int nlayers, int iband, const LwViews& v,
                                           const RrtmgLwTables& CPPWRF_RESTRICT tables) {
  rrtmg_lw_band<MAX_LAY>(nlayers, iband, v.secdiff, v.semiss, v.cldfmc, v.taucmc, v.planklay, v.planklev, v.plankbnd,
                         v.fracs, v.taug, v.icldlyr, v.bu, v.bd, v.bcu, v.bcd, tables, CSI(v.iscal)[4] == 0);
}

// ---- stage 4: per level -------------------------------------------------------------------------
CPPWRF_FORCEINLINE void rrtmg_lw_level(int lev, const LwViews& v, const RrtmgLwTables& CPPWRF_RESTRICT tables) {
  double tu, td, tuc, tdc, fn, fnc;
  rrtmg_lw_accum(lev, v.bu, v.bd, v.bcu, v.bcd, tables, tu, td, tuc, tdc, fn, fnc);
  v.totuflux[lev] = tu; v.totdflux[lev] = td; v.totuclfl[lev] = tuc; v.totdclfl[lev] = tdc;
  v.fnet[lev] = fn; v.fnetc[lev] = fnc;
}

// ---- stage 5: per column ------------------------------------------------------------------------
struct LwColumnOut {
  double glw, olr, lwcf, lwupt, lwuptc, lwdnt, lwdntc, lwupb, lwupbc, lwdnb, lwdnbc;
};
CPPWRF_FORCEINLINE LwColumnOut rrtmg_lw_finish(const RrtmgLwCallConfig& c, const LwViews& v, CSD pi3d,
                                               SD rthratenlw, SD rthratenlwc) {
  const int nlayers = c.nlayers;
  rrtmg_lw_heating(nlayers, v.pz, v.fnet, v.fnetc, v.htr, v.htrc, c.heatfac);
  LwColumnOut o;
  o.glw = v.totdflux[0];
  o.olr = v.totuflux[nlayers];
  o.lwcf = v.totuclfl[nlayers] - v.totuflux[nlayers];
  o.lwupt = v.totuflux[nlayers];
  o.lwuptc = v.totuclfl[nlayers];
  o.lwdnt = v.totdflux[nlayers];
  o.lwdntc = v.totdclfl[nlayers];
  o.lwupb = v.totuflux[0];
  o.lwupbc = v.totuclfl[0];
  o.lwdnb = v.totdflux[0];
  o.lwdnbc = v.totdclfl[0];
  for (int k = 0; k < c.num_levels; ++k) {
    rthratenlw[k] = (v.htr[k] / 86400.0) / pi3d[k];
    if (rthratenlwc.p) rthratenlwc[k] = (v.htrc[k] / 86400.0) / pi3d[k];
  }
  return o;
}

}  // namespace cppwrf::phys
