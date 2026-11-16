// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "phys/ra_rrtmg_common/rrtmg_real.hpp"
// RRTMG_SWRAD (WRF module_ra_rrtmg_sw.F) for one column as independent stages on strided views
// (src/phys/ra_rrtmg_common/strided.hpp), as rrtmg_lw_stages.hpp for the longwave:
//   rrtmg_sw_prep     per column : RRTMG_SWRAD profiles, cloud properties, McICA subcolumns
//                                  (changeSeed = 1, maximum-random), inatm_sw, setcoef_sw, taumol_sw
//   rrtmg_sw_gpoint   per g-point: spcvmc_sw for the g-point (cldprmc_sw of its subcolumn, delta
//                                  scaling, reftra_sw, vrtqdr_sw), its contributions to every level
//   rrtmg_sw_level    per level  : spcvmc_sw's flux sums over the g-points in WRF's order
//   rrtmg_sw_finish   per column : heating rates and the RRTMG_SWRAD outputs
// Night columns (coszen <= 0) skip the work and return zeros as RRTMG_SWRAD.
#include <algorithm>
#include <cmath>
#include <cstddef>

#include "core/backend.hpp"
#include "phys/ra_rrtmg_common/kissvec.hpp"
#include "phys/ra_rrtmg_common/strided.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_stages.hpp"  // RrtmgLwConsts (retab, inirad ozone)
#include "phys/ra_rrtmg_sw/rrtmg_sw_cldprmc.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_reftra.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_setcoef.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_tables.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_taumol.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_vrtqdr.hpp"

namespace cppwrf::phys {

constexpr int RRTMG_SW_COL_MAX = 128;  // layers of one column (stack arrays of rrtmg_sw_gpoint)

struct RrtmgSwCallConfig {
  int num_levels, nlayers;  // nlayers = num_levels + 1
  int o3input;
  rr co2, n2o, ch4;
  rr g;                 // gravmks of the water paths
  rr heatfac;           // rrsw_con heatfac
  rr adjflux;           // adjes*scon/rrsw_scon (dyofyr = 0, adjes = 1)
  int icld, idcor, juldat;  // McICA overlap (cldovrlp), decorrelation length type, day of year
};

struct SwViews {
  SD pz, tz, pavel, tavel, coldry, pdp;
  SD cldfrac, clwpth, ciwpth, cswpth, reicmc, relqmc, resnmc;
  SD2 wkl;  // [7][nl]
  SI jp, jt, jt1, indself, indfor;
  SD colh2o, colco2, colo3, coln2o, colch4, colo2, colmol, co2mult, selffac, selffrac, forfac, forfrac;
  SD fac00, fac01, fac10, fac11;
  SD2 cldfmc;          // [nl][112]
  SD2 taug, taur;      // [nl][112]
  SD sfluxzen;         // [112]
  SD2 cfu, cfd, ccu, ccd, cfdd, ccdd;  // [nl + 1][112] g-point contributions (bottom-to-top levels)
  SD fu, fd, cu, cd, fddir, cddir, uvfd, uvfddir, nifd, nifddir;  // [nl + 1] fluxes
  SD swhr, swhrc;      // [nl]
  SI iscal;  // [0] day (coszen > 0), [1] laytrop, [2] inflag, [3] iceflag, [4] liqflag
  SD dscal;  // [0] prmu0, [1..4] asdir, asdif, aldir, aldif
};

struct SwLayout {
  int nl;
  std::size_t nd = 0, ni = 0;
  CPPWRF_HD explicit SwLayout(int nl_) : nl(nl_) { views(nullptr, nullptr, 0, 1, &nd, &ni); }
  CPPWRF_FORCEINLINE SwViews views(rr* d, int* ip, int col, int ncol, std::size_t* nd_out = nullptr,
                                   std::size_t* ni_out = nullptr) const {
    SwViews v{};
    const int L = nl, L1 = nl + 1, s = ncol;
    std::size_t od = 0, oi = 0;
    auto take = [&](std::size_t n) {
      rr* q = d ? d + od * std::size_t(s) + col : nullptr;
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
    v.pz = v1(L1); v.tz = v1(L1); v.pavel = v1(L1); v.tavel = v1(L1); v.coldry = v1(L1); v.pdp = v1(L1);
    v.cldfrac = v1(L1); v.clwpth = v1(L1); v.ciwpth = v1(L1); v.cswpth = v1(L1);
    v.reicmc = v1(L1); v.relqmc = v1(L1); v.resnmc = v1(L1);
    v.wkl = v2(7, L);
    v.colh2o = v1(L); v.colco2 = v1(L); v.colo3 = v1(L); v.coln2o = v1(L); v.colch4 = v1(L); v.colo2 = v1(L);
    v.colmol = v1(L); v.co2mult = v1(L); v.selffac = v1(L); v.selffrac = v1(L); v.forfac = v1(L);
    v.forfrac = v1(L); v.fac00 = v1(L); v.fac01 = v1(L); v.fac10 = v1(L); v.fac11 = v1(L);
    v.cldfmc = v2(L, 112); v.taug = v2(L, 112); v.taur = v2(L, 112); v.sfluxzen = v1(112);
    v.cfu = v2(L1, 112); v.cfd = v2(L1, 112); v.ccu = v2(L1, 112); v.ccd = v2(L1, 112);
    v.cfdd = v2(L1, 112); v.ccdd = v2(L1, 112);
    v.fu = v1(L1); v.fd = v1(L1); v.cu = v1(L1); v.cd = v1(L1); v.fddir = v1(L1); v.cddir = v1(L1);
    v.uvfd = v1(L1); v.uvfddir = v1(L1); v.nifd = v1(L1); v.nifddir = v1(L1);
    v.swhr = v1(L1); v.swhrc = v1(L1);
    v.dscal = v1(5);
    v.jp = i1(L); v.jt = i1(L); v.jt1 = i1(L); v.indself = i1(L); v.indfor = i1(L); v.iscal = i1(6);
    if (nd_out) *nd_out = od;
    if (ni_out) *ni_out = oi;
    return v;
  }
};

struct SwColumnIn {
  rr coszen, albedo, xland, xice, snow, xlat;
  SD p8w, t8w, p3d, pi3d, dz8w, t3d, qv3d, qc3d, qi3d, qs3d, cldfra3d, o33d;
  SD re_cloud, re_ice, re_snow;  // p == nullptr: has_req* = 0
  bool has_qc, has_qi, has_qs;
  // RRTMG_LWRAD / RRTMG_SWRAD "for MP option 3": without F_QI and warm_rain, below 273.15 K the cloud water is
  // ice and the rain snow (QI1D = QC1D, QS1D = QR1D, QC1D = 0)
  SD qr3d;
  bool split3;
};

// CAM relcalc (module_ra_rrtmg_sw.F) for one layer
CPPWRF_HD inline rr rrtmg_sw_relcalc(rr t_lay, rr xland, rr snow, rr xice) {
  constexpr rr tmelt = rr(273.16), rliqocean = rr(14.0), rliqice = rr(14.0), rliqland = rr(8.0);
  rr landfrac = rr(2.0) - xland;
  rr landm = landfrac;
  rr snowh = rr(0.001) * snow;
  rr icefrac = xice;
  rr r = rliqland + (rliqocean - rliqland) * std::min(rr(1.0), std::max(rr(0.0), (tmelt - t_lay) * rr(0.05)));
  r += (rliqocean - r) * std::min(rr(1.0), std::max(rr(0.0), snowh * rr(10.0)));
  r += (rliqocean - r) * std::min(rr(1.0), std::max(rr(0.0), rr(1.0) - landm));
  r += (rliqice - r) * std::min(rr(1.0), std::max(rr(0.0), icefrac));
  return r;
}

// rrtmg_sw core preparation: inatm_sw (wkl: vmr -> column amounts), setcoef_sw, taumol_sw, the
// cosine of the zenith angle and the band albedos (RRTMG_SWRAD: all four are the albedo).
CPPWRF_HD inline void rrtmg_sw_core_prep(int nlayers, const SwViews& v, const RrtmgSwTables& tables, rr coszen,
                                         rr asdir, rr asdif, rr aldir, rr aldif) {
  v.iscal[0] = 1;
  // 6. inatm_sw
  {
    constexpr rr amd = rr(28.9660);
    constexpr rr amw = rr(18.0160);
    constexpr rr avogad = rr(6.02214199e23);
    constexpr rr grav = rr(9.8066);
    for (int l = 0; l < nlayers; ++l) {
      v.pdp[l] = v.pz[l] - v.pz[l + 1];
      rr amm = (rr(1.0) - v.wkl[0][l]) * amd + v.wkl[0][l] * amw;
      v.coldry[l] = (v.pz[l] - v.pz[l + 1]) * rr(1.0e3) * avogad / (rr(1.0e2) * grav * amm * (rr(1.0) + v.wkl[0][l]));
      for (int imol = 0; imol < 7; ++imol) v.wkl[imol][l] *= v.coldry[l];
    }
  }
  // 7. setcoef_sw, taumol_sw
  int laytrop, layswtch, laylow;
  setcoef_sw(nlayers, v.pavel, v.tavel, v.coldry, v.wkl, laytrop, layswtch, laylow, v.jp, v.jt, v.jt1, v.colh2o,
             v.colco2, v.colo3, v.coln2o, v.colch4, v.colo2, v.colmol, v.co2mult, v.indself, v.indfor, v.selffac,
             v.selffrac, v.forfac, v.forfrac, v.fac00, v.fac01, v.fac10, v.fac11, tables);
  v.iscal[1] = laytrop;
  taumol_sw(nlayers, laytrop, v.jp, v.jt, v.jt1, v.colh2o, v.colco2, v.colch4, v.colo2, v.colo3, v.colmol, v.fac00,
            v.fac01, v.fac10, v.fac11, v.selffac, v.selffrac, v.indself, v.forfac, v.forfrac, v.indfor, v.taug,
            v.taur, v.sfluxzen, tables);
  rr prmu0 = coszen;
  if (prmu0 <= rr(1.e-10)) prmu0 = rr(1.e-10);  // zepzen
  v.dscal[0] = prmu0;
  v.dscal[1] = asdir;
  v.dscal[2] = asdif;
  v.dscal[3] = aldir;
  v.dscal[4] = aldif;
}


// ---- stage 1: per column ------------------------------------------------------------------------
CPPWRF_HD inline void rrtmg_sw_prep(const RrtmgSwCallConfig& c, const SwColumnIn& in, const SwViews& v,
                                    const RrtmgSwTables& tables, const RrtmgLwConsts& kc) {
  v.iscal[0] = in.coszen > rr(0.0) ? 1 : 0;
  if (in.coszen <= rr(0.0)) return;  // the sun below the horizon
  const int num_levels = c.num_levels, nlayers = c.nlayers;
  const rr xland = in.xland;
  auto qc_of = [&](int k) { return in.has_qc ? in.qc3d[k] : rr(0.0); };
  auto qi_of = [&](int k) { return in.has_qi ? in.qi3d[k] : rr(0.0); };
  auto qs_of = [&](int k) { return in.has_qs ? in.qs3d[k] : rr(0.0); };
  constexpr rr amdw = rr(1.607793);
  constexpr rr amdo = rr(0.603461);

  // 1. atmosphere (wkl = vmr until inatm_sw)
  v.pz[0] = in.p8w[0] / rr(100.0);  // Pa -> mb
  v.tz[0] = in.t8w[0];
  for (int k = 0; k < num_levels; ++k) {
    v.pavel[k] = in.p3d[k] / rr(100.0);
    v.pz[k + 1] = in.p8w[k + 1] / rr(100.0);
    v.tavel[k] = in.t3d[k];
    v.tz[k + 1] = in.t8w[k + 1];
    v.wkl[0][k] = std::max(in.qv3d[k], rr(1.0e-12)) * amdw;
    v.wkl[1][k] = c.co2;
    v.wkl[3][k] = c.n2o;
    v.wkl[4][k] = rr(0.0);
    v.wkl[5][k] = c.ch4;
    v.wkl[6][k] = rr(0.209488);
  }
  // 2. extra layer from the model top to the TOA (1e-5 mb)
  v.pavel[num_levels] = rr(0.5) * v.pz[num_levels];
  v.tavel[num_levels] = v.tz[num_levels];
  v.pz[num_levels + 1] = rr(1.0e-5);
  v.tz[num_levels + 1] = v.tz[num_levels];
  for (int m = 0; m < 7; ++m) v.wkl[m][num_levels] = v.wkl[m][num_levels - 1];
  // 3. ozone (inirad, CAM ozone on the model levels for o3input = 2)
  for (int k = 0; k < nlayers; ++k) {
    rr o3 = rr(0.0);
    const rr pr_k = v.pz[k], pr_k1 = v.pz[k + 1];
    for (int jj = 0; jj < 31; ++jj) {
      rr pb1 = (pr_k > kc.ppwrkh[jj]) ? (pr_k - kc.ppwrkh[jj]) : rr(0.0);
      rr pb2 = (pr_k > kc.ppwrkh[jj + 1]) ? (pr_k - kc.ppwrkh[jj + 1]) : rr(0.0);
      rr pt1 = (pr_k1 > kc.ppwrkh[jj]) ? (pr_k1 - kc.ppwrkh[jj]) : rr(0.0);
      rr pt2 = (pr_k1 > kc.ppwrkh[jj + 1]) ? (pr_k1 - kc.ppwrkh[jj + 1]) : rr(0.0);
      o3 += (pb2 - pb1 - pt2 + pt1) * kc.o3ann[jj];
    }
    v.coldry[k] = o3 / (pr_k - pr_k1);  // o3mmr (coldry is set by inatm_sw below)
  }
  for (int k = 0; k < nlayers; ++k) {
    rr o3vmr = v.coldry[k] * amdo;
    if (c.o3input == 2) {
      if (k < num_levels) {
        o3vmr = in.o33d[k];
      } else {
        o3vmr = in.o33d[num_levels - 1] - v.coldry[num_levels - 1] * amdo + v.coldry[k] * amdo;
        if (o3vmr <= rr(0.)) o3vmr = v.coldry[k] * amdo;
      }
    }
    v.wkl[2][k] = o3vmr;
  }
  // 4. cloud physical properties (RRTMG_SWRAD, icloud != 0; EM_CORE branches)
  int inflgsw = 2, iceflgsw = 3;
  const int liqflgsw = 1;
  const bool hrc = in.re_cloud.p != nullptr, hri = in.re_ice.p != nullptr, hrs = in.re_snow.p != nullptr;
  for (int k = 0; k < num_levels; ++k) {
    const rr pdel = v.pz[k] - v.pz[k + 1];
    const rr cldf = in.cldfra3d[k];
    rr qc1 = std::max(rr(0.), qc_of(k)), qi1 = std::max(rr(0.), qi_of(k)), qs1 = std::max(rr(0.), qs_of(k));
    if (in.split3 && in.t3d[k] < rr(273.15)) {  // module_ra_rrtmg_sw.F "For MP option 3" (WSM3)
      qi1 = qc1;
      qs1 = std::max(rr(0.), in.qr3d[k]);
      qc1 = rr(0.);
    }
    rr recloud1d = rr(5.0), reice1d = rr(10.), resnow1d = rr(10.0);
    if (hrc) {
      recloud1d = std::max(rr(2.5), in.re_cloud[k] * rr(1.E6));
      if (recloud1d <= rr(2.5) && cldf > rr(0.) && (xland - rr(1.5)) > rr(0.)) recloud1d = rr(10.5);       // ocean
      else if (recloud1d <= rr(2.5) && cldf > rr(0.) && (xland - rr(1.5)) < rr(0.)) recloud1d = rr(7.5);   // land
      inflgsw = 3;
    }
    if (hri) {
      inflgsw = 4; iceflgsw = 4;
      reice1d = std::max(rr(5.), in.re_ice[k] * rr(1.E6));
      if (reice1d <= rr(5.) && cldf > rr(0.)) {
        int idx_rei = int(in.t3d[k] - rr(179.));
        idx_rei = std::min(std::max(idx_rei, 1), 75);
        const rr corr = in.t3d[k] - int(in.t3d[k]);
        reice1d = kc.retab[idx_rei - 1] * (rr(1.) - corr) + kc.retab[idx_rei] * corr;
        reice1d = std::max(reice1d, rr(5.0));
      }
    }
    if (hrs) {
      inflgsw = 5; iceflgsw = 5;
      resnow1d = std::max(rr(10.), in.re_snow[k] * rr(1.E6));
    }
    if (!hrs && hri && hrc) {  // P3: ice into the snow category for the optics
      inflgsw = 5; iceflgsw = 5;
      resnow1d = std::max(rr(10.), in.re_ice[k] * rr(1.E6));
      qs1 = qi_of(k);
      qi1 = rr(0.);
      reice1d = rr(10.);
    }
    v.cldfrac[k] = cldf;
    rr gicewp = (qi1 + qs1) * pdel * rr(100.0) / c.g * rr(1000.0);  // g/m2
    const rr gliqwp = qc1 * pdel * rr(100.0) / c.g * rr(1000.0);    // g/m2
    v.ciwpth[k] = gicewp / std::max(rr(0.01), cldf);
    v.clwpth[k] = gliqwp / std::max(rr(0.01), cldf);
    if (iceflgsw >= 4) {  // ice path without the snow
      gicewp = qi1 * pdel * rr(100.0) / c.g * rr(1000.0);
      v.ciwpth[k] = gicewp / std::max(rr(0.01), cldf);
    }
    rr csnowp = rr(0.0);
    if (iceflgsw == 5) {
      rr snow_mass_factor = rr(0.99);
      if (resnow1d > rr(130.)) {
        snow_mass_factor = std::min(snow_mass_factor, (rr(130.0) / resnow1d) * (rr(130.0) / resnow1d));
        resnow1d = rr(130.0);
      }
      const rr gsnowp = qs1 * snow_mass_factor * pdel * rr(100.0) / c.g * rr(1000.0);
      csnowp = gsnowp / std::max(rr(0.01), cldf);
    }
    v.relqmc[k] = inflgsw >= 3 ? recloud1d : rrtmg_sw_relcalc(v.tavel[k], xland, in.snow, in.xice);
    if (iceflgsw >= 4) {
      v.reicmc[k] = reice1d;
    } else {  // reicalc, Fu ice size conversion
      const rr t_lay = v.tavel[k];
      int idx = static_cast<int>(t_lay - rr(179.0));
      idx = std::min(std::max(idx, 1), 94);
      rr corr = t_lay - static_cast<int>(t_lay);
      const rr ri = kc.retab[idx - 1] * (rr(1.0) - corr) + kc.retab[idx] * corr;
      v.reicmc[k] = std::min(rr(140.0), ri * rr(1.0315));
    }
    if (inflgsw == 5) {
      v.cswpth[k] = csnowp;
      v.resnmc[k] = resnow1d;
    } else {
      v.cswpth[k] = rr(0.0);
      v.resnmc[k] = rr(10.0);
    }
  }
  v.cldfrac[num_levels] = rr(0.0);
  v.clwpth[num_levels] = rr(0.0);
  v.ciwpth[num_levels] = rr(0.0);
  v.cswpth[num_levels] = rr(0.0);
  v.relqmc[num_levels] = rr(10.0);
  v.reicmc[num_levels] = rr(10.0);
  v.resnmc[num_levels] = rr(10.0);
  v.iscal[2] = inflgsw;
  v.iscal[3] = iceflgsw;
  v.iscal[4] = liqflgsw;

  // 5. McICA subcolumns (KISS, permuteseed = 1, overlap icld), subcolumn by subcolumn; icld = 0: mcica_subcol_sw
  // returns at once and spcvmc_sw gets zcldfmc = 0, zomgcmc = 1 (rrtmg_sw), a clear column
  bool any_cloud = false;
  for (int ilev = 0; ilev < nlayers; ++ilev) {
    if (c.icld != 0 && CSD(v.cldfrac)[ilev] >= rr(1.0e-20)) any_cloud = true;
  }
#if defined(__CUDA_ARCH__)
  const bool run_mcica = __any_sync(__activemask(), any_cloud);
#else
  const bool run_mcica = any_cloud;
#endif
  int any_cld_subcol = 0;
  if (run_mcica && any_cloud) {
    KissRng rng;
    rng.init_from_pressure(CSD(v.pavel)[0] * rr(100.0), CSD(v.pavel)[1] * rr(100.0), CSD(v.pavel)[2] * rr(100.0),
                           CSD(v.pavel)[3] * rr(100.0));
    rng.permute(1);
    // layer mid-point heights from dz8w (exponential overlaps): the buffer layer above the top gets the top dz
    rr hgt[kMcicaMaxLayers];
    rr zo_inv = rr(0.0);
    if (c.icld == 4 || c.icld == 5) {
      rr dzsum = rr(0.0), dz = rr(0.0);
      for (int k = 0; k < nlayers; ++k) {
        if (k < c.num_levels) dz = in.dz8w[k];
        hgt[k] = dzsum + rr(0.5) * dz;
        dzsum = dzsum + dz;
      }
      zo_inv = rr(1.0) / mcica_decorr_len(c.idcor, c.juldat, in.xlat);
    }
    auto cldf_of = [&](int l) { const rr cf = CSD(v.cldfrac)[l]; return cf < rr(1.0e-20) ? rr(0.0) : cf; };
    auto alpha_of = [&](int l) {
      if (c.icld == 5 && cldf_of(l) == rr(0.0) && cldf_of(l - 1) > rr(0.0)) return rr(0.0);
      return libm::exp(-(hgt[l] - hgt[l - 1]) * zo_inv);
    };
    mcica_subcolumns(c.icld, rng, 112, nlayers, cldf_of, alpha_of, [&](int isubcol, int ilev, bool cld) {
      v.cldfmc[ilev][isubcol] = cld ? rr(1.0) : rr(0.0);
      if (cld) any_cld_subcol = 1;
    });
  }
  v.iscal[5] = any_cld_subcol;

  rrtmg_sw_core_prep(nlayers, v, tables, in.coszen, in.albedo, in.albedo, in.albedo, in.albedo);
}

// ---- stage 2: per g-point (spcvmc_sw with icpr = 1, no aerosols) --------------------------------
template <int ML>
CPPWRF_FORCEINLINE void rrtmg_sw_gpoint_impl(const RrtmgSwCallConfig& c, int iw, const SwViews& v,
                                             const RrtmgSwTables& CPPWRF_RESTRICT tables) {
  if (CSI(v.iscal)[0] == 0) return;
  constexpr rr tblint = rr(10000.0);
  constexpr rr od_lo = rr(0.06);
  constexpr rr pade = rr(0.278);
  constexpr rr bpade = rr(1.0) / pade;
  constexpr rr ptaua = rr(0.0), pasya = rr(0.0), pomga = rr(1.0);  // aer_opt = 0
  const int klev = c.nlayers;
  const rr prmu0 = CSD(v.dscal)[0];
  const int inflag = CSI(v.iscal)[2], iceflag = CSI(v.iscal)[3], liqflag = CSI(v.iscal)[4];
  const bool has_col_cld = CSI(v.iscal)[5] != 0;
  int ibm = 0, cum = tables.ngc[0];
  while (iw >= cum) cum += tables.ngc[++ibm];
  // band albedos: WRF bands 1-9 and 14 near-IR (al*), 10-13 UV/visible (as*)
  const bool uv = ibm >= 9 && ibm <= 12;
  const rr palbp = uv ? CSD(v.dscal)[1] : CSD(v.dscal)[3], palbd = uv ? CSD(v.dscal)[2] : CSD(v.dscal)[4];

  rr zref[ML + 1], zrefd[ML + 1], ztra[ML + 1], ztrad[ML + 1];
  rr zdbt[ML + 1], ztdbt[ML + 1];

  const rr zincflx = c.adjflux * CSD(v.sfluxzen)[iw] * prmu0;

  // Clear-sky TOA direct beam and surface values
  ztdbt[0] = rr(1.0);
  zdbt[klev] = rr(0.0);
  ztra[klev] = rr(0.0);
  ztrad[klev] = rr(0.0);
  zref[klev] = palbp;
  zrefd[klev] = palbd;

  rr ztdbtc_nodel = rr(1.0);
  v.ccdd[klev][iw] = zincflx;
  bool any_cld = false;

  for (int jk = 0; jk < klev; ++jk) {
    int ikl = klev - 1 - jk;  // top to bottom
    if (has_col_cld && CSD2(v.cldfmc)[ikl][iw] == rr(1.0)) any_cld = true;
    const rr taur_val = CSD2(v.taur)[ikl][iw];
    const rr taug_val = CSD2(v.taug)[ikl][iw];
    rr ztauc = taur_val + taug_val + ptaua;
    rr zomcc = taur_val * rr(1.0) + ptaua * pomga;
    rr zgcc = (zomcc != rr(0.0)) ? pasya * pomga * ptaua / zomcc : rr(0.0);
    zomcc = (ztauc != rr(0.0)) ? zomcc / ztauc : rr(0.0);

    rr ze1 = ztauc / prmu0;
    rr zdbtmc;
    if (ze1 <= od_lo) {
      zdbtmc = rr(1.0) - ze1 + rr(0.5) * ze1 * ze1;
    } else {
      rr tblind = ze1 / (bpade + ze1);
      int itind = static_cast<int>(tblint * tblind + rr(0.5));
      zdbtmc = tables.exp_tbl[itind];
    }
    ztdbtc_nodel = zdbtmc * ztdbtc_nodel;
    v.ccdd[ikl][iw] = zincflx * ztdbtc_nodel;

    // delta scaling, clear sky
    rr zf = zgcc * zgcc;
    rr zwf = zomcc * zf;
    ztauc = (rr(1.0) - zwf) * ztauc;
    zomcc = (zomcc - zwf) / (rr(1.0) - zwf);
    zgcc = (zgcc - zf) / (rr(1.0) - zf);

    reftra_sw_elem(true, zgcc, prmu0, ztauc, zomcc, zref[jk], zrefd[jk], ztra[jk], ztrad[jk], tables);

    ze1 = ztauc / prmu0;
    if (ze1 <= od_lo) {
      zdbtmc = rr(1.0) - ze1 + rr(0.5) * ze1 * ze1;
    } else {
      rr tblind = ze1 / (bpade + ze1);
      int itind = static_cast<int>(tblint * tblind + rr(0.5));
      zdbtmc = tables.exp_tbl[itind];
    }
    zdbt[jk] = zdbtmc;
    ztdbt[jk + 1] = zdbt[jk] * ztdbt[jk];
  }

  vrtqdr_sw_direct<ML>(klev, zref, zrefd, ztra, ztrad, zdbt, ztdbt, zincflx, v.ccd, v.ccu, iw);

  if (!any_cld) {
    for (int ikl = 0; ikl <= klev; ++ikl) {
      v.cfu[ikl][iw] = v.ccu[ikl][iw];
      v.cfd[ikl][iw] = v.ccd[ikl][iw];
      v.cfdd[ikl][iw] = v.ccdd[ikl][iw];
    }
    return;
  }

  // Cloudy subcolumn: update only cloudy layers in zref/zrefd/ztra/ztrad/zdbt and run total-sky vrtqdr_sw
  rr ztdbt_nodel = rr(1.0);
  v.cfdd[klev][iw] = zincflx;
  ztdbt[0] = rr(1.0);

  for (int jk = 0; jk < klev; ++jk) {
    int ikl = klev - 1 - jk;
    const rr cf = CSD2(v.cldfmc)[ikl][iw];
    const rr taur_val = CSD2(v.taur)[ikl][iw];
    const rr taug_val = CSD2(v.taug)[ikl][iw];
    rr ztauc = taur_val + taug_val + ptaua;

    if (cf == rr(1.0)) {
      rr taorm = rr(0.0), taucm = rr(0.0), ssacm = rr(1.0), asmcm = rr(0.0);
      cldprmc_sw_elem(inflag, iceflag, liqflag, iw, cf, CSD(v.ciwpth)[ikl], CSD(v.clwpth)[ikl],
                      CSD(v.cswpth)[ikl], CSD(v.reicmc)[ikl], CSD(v.relqmc)[ikl], CSD(v.resnmc)[ikl], rr(0.0),
                      taorm, taucm, ssacm, asmcm, tables);
      rr tauorig = ztauc + taorm;
      rr ze1 = tauorig / prmu0;
      rr zdbtmo;
      if (ze1 <= od_lo) {
        zdbtmo = rr(1.0) - ze1 + rr(0.5) * ze1 * ze1;
      } else {
        rr tblind = ze1 / (bpade + ze1);
        int itind = static_cast<int>(tblint * tblind + rr(0.5));
        zdbtmo = tables.exp_tbl[itind];
      }
      ztdbt_nodel = zdbtmo * ztdbt_nodel;
      v.cfdd[ikl][iw] = zincflx * ztdbt_nodel;

      rr zomcc = taur_val * rr(1.0) + ptaua * pomga;
      rr zgcc = (zomcc != rr(0.0)) ? pasya * pomga * ptaua / zomcc : rr(0.0);
      zomcc = (ztauc != rr(0.0)) ? zomcc / ztauc : rr(0.0);
      rr zf = zgcc * zgcc;
      rr zwf = zomcc * zf;
      ztauc = (rr(1.0) - zwf) * ztauc;
      zomcc = (zomcc - zwf) / (rr(1.0) - zwf);
      zgcc = (zgcc - zf) / (rr(1.0) - zf);

      rr ztauo = ztauc + taucm;
      rr zomco = ztauc * zomcc + taucm * ssacm;
      rr zgco = (taucm * ssacm * asmcm + ztauc * zomcc * zgcc) / zomco;
      zomco = zomco / ztauo;

      reftra_sw_elem(true, zgco, prmu0, ztauo, zomco, zref[jk], zrefd[jk], ztra[jk], ztrad[jk], tables);

      ze1 = ztauo / prmu0;
      if (ze1 <= od_lo) {
        zdbtmo = rr(1.0) - ze1 + rr(0.5) * ze1 * ze1;
      } else {
        rr tblind = ze1 / (bpade + ze1);
        int itind = static_cast<int>(tblint * tblind + rr(0.5));
        zdbtmo = tables.exp_tbl[itind];
      }
      zdbt[jk] = zdbtmo;
    } else {
      rr ze1 = ztauc / prmu0;
      rr zdbtmc;
      if (ze1 <= od_lo) {
        zdbtmc = rr(1.0) - ze1 + rr(0.5) * ze1 * ze1;
      } else {
        rr tblind = ze1 / (bpade + ze1);
        int itind = static_cast<int>(tblint * tblind + rr(0.5));
        zdbtmc = tables.exp_tbl[itind];
      }
      ztdbt_nodel = zdbtmc * ztdbt_nodel;
      v.cfdd[ikl][iw] = zincflx * ztdbt_nodel;
    }
    ztdbt[jk + 1] = zdbt[jk] * ztdbt[jk];
  }

  vrtqdr_sw_direct<ML>(klev, zref, zrefd, ztra, ztrad, zdbt, ztdbt, zincflx, v.cfd, v.cfu, iw);
}

template <int ML = RRTMG_SW_COL_MAX>
CPPWRF_FORCEINLINE void rrtmg_sw_gpoint(const RrtmgSwCallConfig& c, int iw, const SwViews& v,
                                        const RrtmgSwTables& CPPWRF_RESTRICT tables) {
  rrtmg_sw_gpoint_impl<ML>(c, iw, v, tables);
}

// ---- stage 3: per level: spcvmc_sw's sums, g-points in order ------------------------------------
CPPWRF_FORCEINLINE void rrtmg_sw_level(int lev, const SwViews& v, const RrtmgSwTables& CPPWRF_RESTRICT tables) {
  if (CSI(v.iscal)[0] == 0) return;
  rr fu = rr(0.0), fd = rr(0.0), cu = rr(0.0), cd = rr(0.0), fddir = rr(0.0), cddir = rr(0.0), uvfd = rr(0.0), uvfddir = rr(0.0), nifd = rr(0.0),
         nifddir = rr(0.0);
  int iw = 0;
  for (int ibm = 0; ibm < 14; ++ibm) {
    const int igt = tables.ngc[ibm];
    for (int jg = 0; jg < igt; ++jg, ++iw) {
      const rr cfu = CSD2(v.cfu)[lev][iw], cfd = CSD2(v.cfd)[lev][iw];
      const rr ccu = CSD2(v.ccu)[lev][iw], ccd = CSD2(v.ccd)[lev][iw];
      const rr cfdd = CSD2(v.cfdd)[lev][iw];
      fu += cfu;
      fd += cfd;
      cu += ccu;
      cd += ccd;
      fddir += cfdd;
      cddir += CSD2(v.ccdd)[lev][iw];
      if (ibm >= 9 && ibm <= 12) {
        uvfd += cfd;
        uvfddir += cfdd;
      } else if (ibm == 13 || ibm <= 8) {
        nifd += cfd;
        nifddir += cfdd;
      }
    }
  }
  v.fu[lev] = fu; v.fd[lev] = fd; v.cu[lev] = cu; v.cd[lev] = cd; v.fddir[lev] = fddir; v.cddir[lev] = cddir;
  v.uvfd[lev] = uvfd; v.uvfddir[lev] = uvfddir; v.nifd[lev] = nifd; v.nifddir[lev] = nifddir;
}

// ---- stage 4: per column ------------------------------------------------------------------------
// surface direct / diffuse fluxes: rrtmg_sw's swdkdir, swdkdif, swdkdirc at level 1 (SWDDIR, SWDDIF, SWDDIRC)
struct SwColumnOut {
  rr gsw, swcf, swupt, swuptc, swdnt, swdntc, swupb, swupbc, swdnb, swdnbc, swvisdir, swvisdif, swnirdir,
      swnirdif, swddir, swddif, swddirc;
};
CPPWRF_FORCEINLINE SwColumnOut rrtmg_sw_finish(const RrtmgSwCallConfig& c, const SwViews& v, CSD pi3d,
                                               SD rthratensw, SD rthratenswc) {
  SwColumnOut o{};
  if (v.iscal[0] == 0) {
    for (int k = 0; k < c.num_levels; ++k) {
      rthratensw[k] = rr(0.0);
      if (rthratenswc.p) rthratenswc[k] = rr(0.0);
    }
    return o;
  }
  const int nlayers = c.nlayers;
  for (int i = 0; i < nlayers; ++i) {
    const rr zdpgcp = c.heatfac / v.pdp[i];
    v.swhrc[i] = ((v.cd[i + 1] - v.cu[i + 1]) - (v.cd[i] - v.cu[i])) * zdpgcp;
    v.swhr[i] = ((v.fd[i + 1] - v.fu[i + 1]) - (v.fd[i] - v.fu[i])) * zdpgcp;
  }
  v.swhrc[nlayers - 1] = rr(0.);
  v.swhr[nlayers - 1] = rr(0.);
  o.gsw = v.fd[0] - v.fu[0];
  o.swcf = (v.fd[nlayers] - v.fu[nlayers]) - (v.cd[nlayers] - v.cu[nlayers]);
  o.swupt = v.fu[nlayers];
  o.swuptc = v.cu[nlayers];
  o.swdnt = v.fd[nlayers];
  o.swdntc = v.cd[nlayers];
  o.swupb = v.fu[0];
  o.swupbc = v.cu[0];
  o.swdnb = v.fd[0];
  o.swdnbc = v.cd[0];
  o.swvisdir = v.uvfddir[0];
  o.swvisdif = v.uvfd[0] - v.uvfddir[0];
  o.swnirdir = v.nifddir[0];
  o.swnirdif = v.nifd[0] - v.nifddir[0];
  o.swddir = v.fddir[0];
  o.swddif = v.fd[0] - v.fddir[0];
  o.swddirc = v.cddir[0];
  for (int k = 0; k < c.num_levels; ++k) {
    rthratensw[k] = (v.swhr[k] / rr(86400.)) / pi3d[k];
    if (rthratenswc.p) rthratenswc[k] = (v.swhrc[k] / rr(86400.)) / pi3d[k];
  }
  return o;
}

}  // namespace cppwrf::phys
