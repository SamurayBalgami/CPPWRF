#pragma once
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
  double co2, n2o, ch4;
  double g;                 // gravmks of the water paths
  double heatfac;           // rrsw_con heatfac
  double adjflux;           // adjes*scon/rrsw_scon (dyofyr = 0, adjes = 1)
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
  SD fu, fd, cu, cd, fddir, uvfd, uvfddir, nifd, nifddir;  // [nl + 1] fluxes
  SD swhr, swhrc;      // [nl]
  SI iscal;  // [0] day (coszen > 0), [1] laytrop, [2] inflag, [3] iceflag, [4] liqflag
  SD dscal;  // [0] prmu0, [1..4] asdir, asdif, aldir, aldif
};

struct SwLayout {
  int nl;
  std::size_t nd = 0, ni = 0;
  CPPWRF_HD explicit SwLayout(int nl_) : nl(nl_) { views(nullptr, nullptr, 0, 1, &nd, &ni); }
  CPPWRF_FORCEINLINE SwViews views(double* d, int* ip, int col, int ncol, std::size_t* nd_out = nullptr,
                                   std::size_t* ni_out = nullptr) const {
    SwViews v{};
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
    v.fu = v1(L1); v.fd = v1(L1); v.cu = v1(L1); v.cd = v1(L1); v.fddir = v1(L1); v.uvfd = v1(L1);
    v.uvfddir = v1(L1); v.nifd = v1(L1); v.nifddir = v1(L1);
    v.swhr = v1(L1); v.swhrc = v1(L1);
    v.dscal = v1(5);
    v.jp = i1(L); v.jt = i1(L); v.jt1 = i1(L); v.indself = i1(L); v.indfor = i1(L); v.iscal = i1(6);
    if (nd_out) *nd_out = od;
    if (ni_out) *ni_out = oi;
    return v;
  }
};

struct SwColumnIn {
  double coszen, albedo, xland, xice, snow;
  SD p8w, t8w, p3d, pi3d, t3d, qv3d, qc3d, qi3d, qs3d, cldfra3d, o33d;
  SD re_cloud, re_ice, re_snow;  // p == nullptr: has_req* = 0
  bool has_qc, has_qi, has_qs;
};

// CAM relcalc for one layer (rrtmg_relcalc_single of rrtmg_common.hpp)
CPPWRF_HD inline double rrtmg_sw_relcalc(double t_lay, double xland, double snow, double xice) {
  constexpr double tmelt = 273.16, rliqocean = 14.0, rliqice = 14.0, rliqland = 8.0;
  double landfrac = 2.0 - xland;
  double landm = landfrac;
  double snowh = 0.001 * snow;
  double icefrac = xice;
  double r = rliqland + (rliqocean - rliqland) * std::min(1.0, std::max(0.0, (tmelt - t_lay) * 0.05));
  r += (rliqocean - r) * std::min(1.0, std::max(0.0, snowh * 10.0));
  r += (rliqocean - r) * std::min(1.0, std::max(0.0, 1.0 - landm));
  r += (rliqice - r) * std::min(1.0, std::max(0.0, icefrac));
  return r;
}

// rrtmg_sw core preparation: inatm_sw (wkl: vmr -> column amounts), setcoef_sw, taumol_sw, the
// cosine of the zenith angle and the band albedos (RRTMG_SWRAD: all four are the albedo).
CPPWRF_HD inline void rrtmg_sw_core_prep(int nlayers, const SwViews& v, const RrtmgSwTables& tables, double coszen,
                                         double asdir, double asdif, double aldir, double aldif) {
  v.iscal[0] = 1;
  // 6. inatm_sw
  {
    constexpr double amd = 28.9660;
    constexpr double amw = 18.0160;
    constexpr double avogad = 6.02214199e23;
    constexpr double grav = 9.8066;
    for (int l = 0; l < nlayers; ++l) {
      v.pdp[l] = v.pz[l] - v.pz[l + 1];
      double amm = (1.0 - v.wkl[0][l]) * amd + v.wkl[0][l] * amw;
      v.coldry[l] = (v.pz[l] - v.pz[l + 1]) * 1.0e3 * avogad / (1.0e2 * grav * amm * (1.0 + v.wkl[0][l]));
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
  double prmu0 = coszen;
  if (prmu0 <= 1.e-10) prmu0 = 1.e-10;  // zepzen
  v.dscal[0] = prmu0;
  v.dscal[1] = asdir;
  v.dscal[2] = asdif;
  v.dscal[3] = aldir;
  v.dscal[4] = aldif;
}


// ---- stage 1: per column ------------------------------------------------------------------------
CPPWRF_HD inline void rrtmg_sw_prep(const RrtmgSwCallConfig& c, const SwColumnIn& in, const SwViews& v,
                                    const RrtmgSwTables& tables, const RrtmgLwConsts& kc) {
  v.iscal[0] = in.coszen > 0.0 ? 1 : 0;
  if (in.coszen <= 0.0) return;  // the sun below the horizon
  const int num_levels = c.num_levels, nlayers = c.nlayers;
  const double xland = in.xland;
  auto qc_of = [&](int k) { return in.has_qc ? in.qc3d[k] : 0.0; };
  auto qi_of = [&](int k) { return in.has_qi ? in.qi3d[k] : 0.0; };
  auto qs_of = [&](int k) { return in.has_qs ? in.qs3d[k] : 0.0; };
  constexpr double amdw = 1.607793;
  constexpr double amdo = 0.603461;

  // 1. atmosphere (wkl = vmr until inatm_sw)
  v.pz[0] = in.p8w[0] / 100.0;  // Pa -> mb
  v.tz[0] = in.t8w[0];
  for (int k = 0; k < num_levels; ++k) {
    v.pavel[k] = in.p3d[k] / 100.0;
    v.pz[k + 1] = in.p8w[k + 1] / 100.0;
    v.tavel[k] = in.t3d[k];
    v.tz[k + 1] = in.t8w[k + 1];
    v.wkl[0][k] = std::max(in.qv3d[k], 1.0e-12) * amdw;
    v.wkl[1][k] = c.co2;
    v.wkl[3][k] = c.n2o;
    v.wkl[4][k] = 0.0;
    v.wkl[5][k] = c.ch4;
    v.wkl[6][k] = 0.209488;
  }
  // 2. extra layer from the model top to the TOA (1e-5 mb)
  v.pavel[num_levels] = 0.5 * v.pz[num_levels];
  v.tavel[num_levels] = v.tz[num_levels];
  v.pz[num_levels + 1] = 1.0e-5;
  v.tz[num_levels + 1] = v.tz[num_levels];
  for (int m = 0; m < 7; ++m) v.wkl[m][num_levels] = v.wkl[m][num_levels - 1];
  // 3. ozone (inirad, CAM ozone on the model levels for o3input = 2)
  for (int k = 0; k < nlayers; ++k) {
    double o3 = 0.0;
    const double pr_k = v.pz[k], pr_k1 = v.pz[k + 1];
    for (int jj = 0; jj < 31; ++jj) {
      double pb1 = (pr_k > kc.ppwrkh[jj]) ? (pr_k - kc.ppwrkh[jj]) : 0.0;
      double pb2 = (pr_k > kc.ppwrkh[jj + 1]) ? (pr_k - kc.ppwrkh[jj + 1]) : 0.0;
      double pt1 = (pr_k1 > kc.ppwrkh[jj]) ? (pr_k1 - kc.ppwrkh[jj]) : 0.0;
      double pt2 = (pr_k1 > kc.ppwrkh[jj + 1]) ? (pr_k1 - kc.ppwrkh[jj + 1]) : 0.0;
      o3 += (pb2 - pb1 - pt2 + pt1) * kc.o3ann[jj];
    }
    v.coldry[k] = o3 / (pr_k - pr_k1);  // o3mmr (coldry is set by inatm_sw below)
  }
  for (int k = 0; k < nlayers; ++k) {
    double o3vmr = v.coldry[k] * amdo;
    if (c.o3input == 2) {
      if (k < num_levels) {
        o3vmr = in.o33d[k];
      } else {
        o3vmr = in.o33d[num_levels - 1] - v.coldry[num_levels - 1] * amdo + v.coldry[k] * amdo;
        if (o3vmr <= 0.) o3vmr = v.coldry[k] * amdo;
      }
    }
    v.wkl[2][k] = o3vmr;
  }
  // 4. cloud physical properties (RRTMG_SWRAD, icloud != 0; EM_CORE branches)
  int inflgsw = 2, iceflgsw = 3;
  const int liqflgsw = 1;
  const bool hrc = in.re_cloud.p != nullptr, hri = in.re_ice.p != nullptr, hrs = in.re_snow.p != nullptr;
  for (int k = 0; k < num_levels; ++k) {
    const double pdel = v.pz[k] - v.pz[k + 1];
    const double cldf = in.cldfra3d[k];
    double qc1 = std::max(0., qc_of(k)), qi1 = std::max(0., qi_of(k)), qs1 = std::max(0., qs_of(k));
    double recloud1d = 5.0, reice1d = 10., resnow1d = 10.0;
    if (hrc) {
      recloud1d = std::max(2.5, in.re_cloud[k] * 1.E6);
      if (recloud1d <= 2.5 && cldf > 0. && (xland - 1.5) > 0.) recloud1d = 10.5;       // ocean
      else if (recloud1d <= 2.5 && cldf > 0. && (xland - 1.5) < 0.) recloud1d = 7.5;   // land
      inflgsw = 3;
    }
    if (hri) {
      inflgsw = 4; iceflgsw = 4;
      reice1d = std::max(5., in.re_ice[k] * 1.E6);
      if (reice1d <= 5. && cldf > 0.) {
        int idx_rei = int(in.t3d[k] - 179.);
        idx_rei = std::min(std::max(idx_rei, 1), 75);
        const double corr = in.t3d[k] - int(in.t3d[k]);
        reice1d = kc.retab[idx_rei - 1] * (1. - corr) + kc.retab[idx_rei] * corr;
        reice1d = std::max(reice1d, 5.0);
      }
    }
    if (hrs) {
      inflgsw = 5; iceflgsw = 5;
      resnow1d = std::max(10., in.re_snow[k] * 1.E6);
    }
    if (!hrs && hri && hrc) {  // P3: ice into the snow category for the optics
      inflgsw = 5; iceflgsw = 5;
      resnow1d = std::max(10., in.re_ice[k] * 1.E6);
      qs1 = qi_of(k);
      qi1 = 0.;
      reice1d = 10.;
    }
    v.cldfrac[k] = cldf;
    double gicewp = (qi1 + qs1) * pdel * 100.0 / c.g * 1000.0;  // g/m2
    const double gliqwp = qc1 * pdel * 100.0 / c.g * 1000.0;    // g/m2
    v.ciwpth[k] = gicewp / std::max(0.01, cldf);
    v.clwpth[k] = gliqwp / std::max(0.01, cldf);
    if (iceflgsw >= 4) {  // ice path without the snow
      gicewp = qi1 * pdel * 100.0 / c.g * 1000.0;
      v.ciwpth[k] = gicewp / std::max(0.01, cldf);
    }
    double csnowp = 0.0;
    if (iceflgsw == 5) {
      double snow_mass_factor = 0.99;
      if (resnow1d > 130.) {
        snow_mass_factor = std::min(snow_mass_factor, (130.0 / resnow1d) * (130.0 / resnow1d));
        resnow1d = 130.0;
      }
      const double gsnowp = qs1 * snow_mass_factor * pdel * 100.0 / c.g * 1000.0;
      csnowp = gsnowp / std::max(0.01, cldf);
    }
    v.relqmc[k] = inflgsw >= 3 ? recloud1d : rrtmg_sw_relcalc(v.tavel[k], xland, in.snow, in.xice);
    if (iceflgsw >= 4) {
      v.reicmc[k] = reice1d;
    } else {  // reicalc, Fu ice size conversion
      const double t_lay = v.tavel[k];
      int idx = static_cast<int>(t_lay - 179.0);
      idx = std::min(std::max(idx, 1), 94);
      double corr = t_lay - static_cast<int>(t_lay);
      const double ri = kc.retab[idx - 1] * (1.0 - corr) + kc.retab[idx] * corr;
      v.reicmc[k] = std::min(140.0, ri * 1.0315);
    }
    if (inflgsw == 5) {
      v.cswpth[k] = csnowp;
      v.resnmc[k] = resnow1d;
    } else {
      v.cswpth[k] = 0.0;
      v.resnmc[k] = 10.0;
    }
  }
  v.cldfrac[num_levels] = 0.0;
  v.clwpth[num_levels] = 0.0;
  v.ciwpth[num_levels] = 0.0;
  v.cswpth[num_levels] = 0.0;
  v.relqmc[num_levels] = 10.0;
  v.reicmc[num_levels] = 10.0;
  v.resnmc[num_levels] = 10.0;
  v.iscal[2] = inflgsw;
  v.iscal[3] = iceflgsw;
  v.iscal[4] = liqflgsw;

  // 5. McICA subcolumns (KISS, permuteseed = 1, maximum-random overlap), subcolumn by subcolumn
  bool any_cloud = false;
  for (int ilev = 0; ilev < nlayers; ++ilev) {
    if (CSD(v.cldfrac)[ilev] >= 1.0e-20) any_cloud = true;
  }
#if defined(__CUDA_ARCH__)
  const bool run_mcica = __any_sync(__activemask(), any_cloud);
#else
  const bool run_mcica = any_cloud;
#endif
  int any_cld_subcol = 0;
  if (run_mcica && any_cloud) {
    KissRng rng;
    rng.init_from_pressure(CSD(v.pavel)[0] * 100.0, CSD(v.pavel)[1] * 100.0, CSD(v.pavel)[2] * 100.0,
                           CSD(v.pavel)[3] * 100.0);
    rng.permute(1);
    for (int isubcol = 0; isubcol < 112; ++isubcol) {
      double cdf_prev = 0.0, cldf_prev = 0.0;
      for (int ilev = 0; ilev < nlayers; ++ilev) {
        const double cf = CSD(v.cldfrac)[ilev];
        const double cldf = (cf < 1.0e-20) ? 0.0 : cf;
        double cdf = rng.step<double>();
        if (ilev > 0) {
          if (cdf_prev > 1.0 - cldf_prev) cdf = cdf_prev;
          else cdf = cdf * (1.0 - cldf_prev);
        }
        const bool cld = cdf >= 1.0 - cldf;
        v.cldfmc[ilev][isubcol] = cld ? 1.0 : 0.0;
        if (cld) any_cld_subcol = 1;
        cdf_prev = cdf;
        cldf_prev = cldf;
      }
    }
  }
  v.iscal[5] = any_cld_subcol;

  rrtmg_sw_core_prep(nlayers, v, tables, in.coszen, in.albedo, in.albedo, in.albedo, in.albedo);
}

// ---- stage 2: per g-point (spcvmc_sw with icpr = 1, no aerosols) --------------------------------
template <int ML>
CPPWRF_FORCEINLINE void rrtmg_sw_gpoint_impl(const RrtmgSwCallConfig& c, int iw, const SwViews& v,
                                             const RrtmgSwTables& CPPWRF_RESTRICT tables) {
  if (CSI(v.iscal)[0] == 0) return;
  constexpr double tblint = 10000.0;
  constexpr double od_lo = 0.06;
  constexpr double pade = 0.278;
  constexpr double bpade = 1.0 / pade;
  constexpr double ptaua = 0.0, pasya = 0.0, pomga = 1.0;  // aer_opt = 0
  const int klev = c.nlayers;
  const double prmu0 = CSD(v.dscal)[0];
  const int inflag = CSI(v.iscal)[2], iceflag = CSI(v.iscal)[3], liqflag = CSI(v.iscal)[4];
  const bool has_col_cld = CSI(v.iscal)[5] != 0;
  int ibm = 0, cum = tables.ngc[0];
  while (iw >= cum) cum += tables.ngc[++ibm];
  // band albedos: WRF bands 1-9 and 14 near-IR (al*), 10-13 UV/visible (as*)
  const bool uv = ibm >= 9 && ibm <= 12;
  const double palbp = uv ? CSD(v.dscal)[1] : CSD(v.dscal)[3], palbd = uv ? CSD(v.dscal)[2] : CSD(v.dscal)[4];

  double zref[ML + 1], zrefd[ML + 1], ztra[ML + 1], ztrad[ML + 1];
  double zdbt[ML + 1], ztdbt[ML + 1];

  const double zincflx = c.adjflux * CSD(v.sfluxzen)[iw] * prmu0;

  // Clear-sky TOA direct beam and surface values
  ztdbt[0] = 1.0;
  zdbt[klev] = 0.0;
  ztra[klev] = 0.0;
  ztrad[klev] = 0.0;
  zref[klev] = palbp;
  zrefd[klev] = palbd;

  double ztdbtc_nodel = 1.0;
  v.ccdd[klev][iw] = zincflx;
  bool any_cld = false;

  for (int jk = 0; jk < klev; ++jk) {
    int ikl = klev - 1 - jk;  // top to bottom
    if (has_col_cld && CSD2(v.cldfmc)[ikl][iw] == 1.0) any_cld = true;
    const double taur_val = CSD2(v.taur)[ikl][iw];
    const double taug_val = CSD2(v.taug)[ikl][iw];
    double ztauc = taur_val + taug_val + ptaua;
    double zomcc = taur_val * 1.0 + ptaua * pomga;
    double zgcc = (zomcc != 0.0) ? pasya * pomga * ptaua / zomcc : 0.0;
    zomcc = (ztauc != 0.0) ? zomcc / ztauc : 0.0;

    double ze1 = ztauc / prmu0;
    double zdbtmc;
    if (ze1 <= od_lo) {
      zdbtmc = 1.0 - ze1 + 0.5 * ze1 * ze1;
    } else {
      double tblind = ze1 / (bpade + ze1);
      int itind = static_cast<int>(tblint * tblind + 0.5);
      zdbtmc = tables.exp_tbl[itind];
    }
    ztdbtc_nodel = zdbtmc * ztdbtc_nodel;
    v.ccdd[ikl][iw] = zincflx * ztdbtc_nodel;

    // delta scaling, clear sky
    double zf = zgcc * zgcc;
    double zwf = zomcc * zf;
    ztauc = (1.0 - zwf) * ztauc;
    zomcc = (zomcc - zwf) / (1.0 - zwf);
    zgcc = (zgcc - zf) / (1.0 - zf);

    reftra_sw_elem(true, zgcc, prmu0, ztauc, zomcc, zref[jk], zrefd[jk], ztra[jk], ztrad[jk], tables);

    ze1 = ztauc / prmu0;
    if (ze1 <= od_lo) {
      zdbtmc = 1.0 - ze1 + 0.5 * ze1 * ze1;
    } else {
      double tblind = ze1 / (bpade + ze1);
      int itind = static_cast<int>(tblint * tblind + 0.5);
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
  double ztdbt_nodel = 1.0;
  v.cfdd[klev][iw] = zincflx;
  ztdbt[0] = 1.0;

  for (int jk = 0; jk < klev; ++jk) {
    int ikl = klev - 1 - jk;
    const double cf = CSD2(v.cldfmc)[ikl][iw];
    const double taur_val = CSD2(v.taur)[ikl][iw];
    const double taug_val = CSD2(v.taug)[ikl][iw];
    double ztauc = taur_val + taug_val + ptaua;

    if (cf == 1.0) {
      double taorm = 0.0, taucm = 0.0, ssacm = 1.0, asmcm = 0.0;
      cldprmc_sw_elem(inflag, iceflag, liqflag, iw, cf, CSD(v.ciwpth)[ikl], CSD(v.clwpth)[ikl],
                      CSD(v.cswpth)[ikl], CSD(v.reicmc)[ikl], CSD(v.relqmc)[ikl], CSD(v.resnmc)[ikl], 0.0,
                      taorm, taucm, ssacm, asmcm, tables);
      double tauorig = ztauc + taorm;
      double ze1 = tauorig / prmu0;
      double zdbtmo;
      if (ze1 <= od_lo) {
        zdbtmo = 1.0 - ze1 + 0.5 * ze1 * ze1;
      } else {
        double tblind = ze1 / (bpade + ze1);
        int itind = static_cast<int>(tblint * tblind + 0.5);
        zdbtmo = tables.exp_tbl[itind];
      }
      ztdbt_nodel = zdbtmo * ztdbt_nodel;
      v.cfdd[ikl][iw] = zincflx * ztdbt_nodel;

      double zomcc = taur_val * 1.0 + ptaua * pomga;
      double zgcc = (zomcc != 0.0) ? pasya * pomga * ptaua / zomcc : 0.0;
      zomcc = (ztauc != 0.0) ? zomcc / ztauc : 0.0;
      double zf = zgcc * zgcc;
      double zwf = zomcc * zf;
      ztauc = (1.0 - zwf) * ztauc;
      zomcc = (zomcc - zwf) / (1.0 - zwf);
      zgcc = (zgcc - zf) / (1.0 - zf);

      double ztauo = ztauc + taucm;
      double zomco = ztauc * zomcc + taucm * ssacm;
      double zgco = (taucm * ssacm * asmcm + ztauc * zomcc * zgcc) / zomco;
      zomco = zomco / ztauo;

      reftra_sw_elem(true, zgco, prmu0, ztauo, zomco, zref[jk], zrefd[jk], ztra[jk], ztrad[jk], tables);

      ze1 = ztauo / prmu0;
      if (ze1 <= od_lo) {
        zdbtmo = 1.0 - ze1 + 0.5 * ze1 * ze1;
      } else {
        double tblind = ze1 / (bpade + ze1);
        int itind = static_cast<int>(tblint * tblind + 0.5);
        zdbtmo = tables.exp_tbl[itind];
      }
      zdbt[jk] = zdbtmo;
    } else {
      double ze1 = ztauc / prmu0;
      double zdbtmc;
      if (ze1 <= od_lo) {
        zdbtmc = 1.0 - ze1 + 0.5 * ze1 * ze1;
      } else {
        double tblind = ze1 / (bpade + ze1);
        int itind = static_cast<int>(tblint * tblind + 0.5);
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
  double fu = 0.0, fd = 0.0, cu = 0.0, cd = 0.0, fddir = 0.0, uvfd = 0.0, uvfddir = 0.0, nifd = 0.0, nifddir = 0.0;
  int iw = 0;
  for (int ibm = 0; ibm < 14; ++ibm) {
    const int igt = tables.ngc[ibm];
    for (int jg = 0; jg < igt; ++jg, ++iw) {
      const double cfu = CSD2(v.cfu)[lev][iw], cfd = CSD2(v.cfd)[lev][iw];
      const double ccu = CSD2(v.ccu)[lev][iw], ccd = CSD2(v.ccd)[lev][iw];
      const double cfdd = CSD2(v.cfdd)[lev][iw];
      fu += cfu;
      fd += cfd;
      cu += ccu;
      cd += ccd;
      fddir += cfdd;
      if (ibm >= 9 && ibm <= 12) {
        uvfd += cfd;
        uvfddir += cfdd;
      } else if (ibm == 13 || ibm <= 8) {
        nifd += cfd;
        nifddir += cfdd;
      }
    }
  }
  v.fu[lev] = fu; v.fd[lev] = fd; v.cu[lev] = cu; v.cd[lev] = cd; v.fddir[lev] = fddir;
  v.uvfd[lev] = uvfd; v.uvfddir[lev] = uvfddir; v.nifd[lev] = nifd; v.nifddir[lev] = nifddir;
}

// ---- stage 4: per column ------------------------------------------------------------------------
struct SwColumnOut {
  double gsw, swcf, swupt, swuptc, swdnt, swdntc, swupb, swupbc, swdnb, swdnbc, swvisdir, swvisdif, swnirdir,
      swnirdif;
};
CPPWRF_FORCEINLINE SwColumnOut rrtmg_sw_finish(const RrtmgSwCallConfig& c, const SwViews& v, CSD pi3d,
                                               SD rthratensw, SD rthratenswc) {
  SwColumnOut o{};
  if (v.iscal[0] == 0) {
    for (int k = 0; k < c.num_levels; ++k) {
      rthratensw[k] = 0.0;
      if (rthratenswc.p) rthratenswc[k] = 0.0;
    }
    return o;
  }
  const int nlayers = c.nlayers;
  for (int i = 0; i < nlayers; ++i) {
    const double zdpgcp = c.heatfac / v.pdp[i];
    v.swhrc[i] = ((v.cd[i + 1] - v.cu[i + 1]) - (v.cd[i] - v.cu[i])) * zdpgcp;
    v.swhr[i] = ((v.fd[i + 1] - v.fu[i + 1]) - (v.fd[i] - v.fu[i])) * zdpgcp;
  }
  v.swhrc[nlayers - 1] = 0.;
  v.swhr[nlayers - 1] = 0.;
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
  for (int k = 0; k < c.num_levels; ++k) {
    rthratensw[k] = (v.swhr[k] / 86400.) / pi3d[k];
    if (rthratenswc.p) rthratenswc[k] = (v.swhrc[k] / 86400.) / pi3d[k];
  }
  return o;
}

}  // namespace cppwrf::phys
