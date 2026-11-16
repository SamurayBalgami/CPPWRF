// SPDX-License-Identifier: GPL-3.0-or-later
#include "phys/ra_rrtmg_common/rrtmg_real.hpp"
#include "phys/driver/radiation.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>

#include "core/data_dir.hpp"
#include "core/tiles.hpp"
#include "phys/phys_constants.hpp"
#include "phys/ra_rrtm_lw/rrtm_lw.hpp"
#include "phys/ra_rrtm_lw/rrtm_lw_tables.hpp"
#include "phys/ra_rrtmg_common/cam_gases.hpp"
#include "phys/ra_rrtmg_common/cam_ozone.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_tables.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_tables.hpp"
#include "phys/ra_sw_dudhia/dudhia_sw.hpp"
#if defined(CPPWRF_USE_CUDA)
#include "phys/ra_rrtmg_lw/rrtmg_lw_batch.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_batch.hpp"
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

namespace cppwrf::phys::driver {

namespace {
constexpr rr piconst = rr(3.1415926535897932384626433);
constexpr rr degrad = piconst / rr(180.);
constexpr rr dpd = rr(360.) / rr(365.);

// radconst: declination (rad) and solar constant from the julian day
void radconst(rr julian, rr& declin, rr& solcon) {
  const rr obecl = rr(23.5) * degrad, sinob = libm::sin(obecl);
  rr sxlong = julian >= rr(80.) ? dpd * (julian - rr(80.)) : dpd * (julian + rr(285.));
  sxlong = sxlong * degrad;
  const rr arg = sinob * libm::sin(sxlong);
  declin = libm::asin(arg);
  const rr djul = julian * rr(360.) / rr(365.), rjul = djul * degrad;
  const rr eccfac = rr(1.000110) + rr(0.034221) * libm::cos(rjul) + rr(0.001280) * libm::sin(rjul) +
                        rr(0.000719) * libm::cos(2 * rjul) + rr(0.000077) * libm::sin(2 * rjul);
  solcon = rr(1370.) * eccfac;
}


#if defined(CPPWRF_USE_CUDA)
// device path: model fields -> rr arrays with the column index fastest: the columns the worker owns (a
// decomposed run: its patch [ci0, ci0 + nx) x [cj0, ...); else the domain), col = (i - ci0) + nx (j - cj0)
struct KStageAll {
  int nx, nz, ncol, n3;
  View3<const real> f3w[2];
  rr* d3w[2];
  View3<const real> f3[13];
  rr* d3[13];
  View2<const real> f2[6];
  rr* d2[6];
  rr* dht = nullptr;  // Dudhia: the terrain height
  BView2<const real> ht;
  // a tile stages the whole rows of its patch (the batch's columns i + nx j are contiguous by row); a column it does
  // not own gets the inputs of its nearest own column in the row (valid values; its results are not used)
  int oi0 = -2147483647 - 1, oi1 = 2147483647;
  int ci0 = 0, cj0 = 0;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    const int col = (i - ci0) + nx * (j - cj0);
    const std::size_t idx = std::size_t(k) * ncol + col;
    i = i < oi0 ? oi0 : (i >= oi1 ? oi1 - 1 : i);
    d3w[0][idx] = rr(f3w[0](i, j, k));
    d3w[1][idx] = rr(f3w[1](i, j, k));
    if (k < nz) {
      for (int m = 0; m < n3; ++m) d3[m][idx] = rr(f3[m](i, j, k));
    }
    if (k == 0) {
      for (int m = 0; m < 6; ++m) d2[m][col] = rr(f2[m](i, j));
      if (dht) dht[col] = rr(ht(i, j));
    }
  }
};
// per column: calc_coszen at xtime + radt/2, CAM ozone (ozn_time_int, ozn_p_int), cal_cldfra1
struct KRadColumn {
  int nx, nz, ncol;
  rr julian, xtime_mid, gmt, declin;
  const rr *xlat, *xlong;
  int o3input;
  const rr *ozmixm, *pin;
  int nm, np;
  rr f1, f2;
  bool clouds;
  const rr *p, *t, *qv, *qc, *qi, *qs;
  rr *o3, *cf, *cz;
  bool has_cldfra;
  View3<real> cldfra_out;
  View2<real> coszen_out;
  bool has_diag;
  View2<real> hrang_out;
  View3<real> o3_out;
  bool o3_field;  // a nest: o3rad as the parent forced it (WRF computes it on domain 1 only)
  View3<const real> o3_in;
  bool bl_clouds;  // icloud_bl: MYNN subgrid clouds into cf and the staged qc, qi
  long itimestep;
  View3<const real> cfbl, qcbl, qibl;
  rr *qcw, *qiw;
  int icloud_cu;  // cu_rad_feedback: the cumulus clouds into cf and the staged qc (qcw), qi (qiw)
  View3<const real> qccu, qicu, cfdp, cfsh;
  int ci0 = 0, cj0 = 0, gnx = 0;  // the patch's first column; the domain's nx (calc_coszen's lanes)
  bool bmj;  // bmj_rad_feedback: BMJ's clouds into the staged qc, qi and cf
  View3<const real> qcconv, qiconv, ccldfra;
  CPPWRF_HD void operator()(int i, int j) const {
    const int col = (i - ci0) + nx * (j - cj0);
    rr hr = 0;
    // lane i of the nx-trip i loop, as the host path: gfortran's 2-lane sin / cos of calc_coszen
    const rr c = calc_coszen(julian, xtime_mid, gmt, declin, xlong[col], xlat[col], i, gnx, &hr);
    if (has_diag) hrang_out(i, j) = real(hr);
    cz[col] = c;
    coszen_out(i, j) = real(c);
    constexpr int per = CamOzone::nmon * CamOzone::levsiz;
    auto v = [&](const rr* a) { return SD{const_cast<rr*>(a) + col, ncol}; };
    if (o3input == 2 && o3_field)
      for (int k = 0; k < nz; ++k) o3[std::size_t(k) * ncol + col] = rr(o3_in(i, j, k));
    else if (o3input == 2)
      ozone_column(&ozmixm[std::size_t(col) * per + nm * CamOzone::levsiz],
                   &ozmixm[std::size_t(col) * per + np * CamOzone::levsiz], f1, f2, pin, nz, v(p), v(o3));
    if (icloud_cu == 1) {  // Grell schemes: qc = qc + qc_cu, qi = qi + qi_cu (F_QC / F_QI)
      for (int k = 0; k < nz; ++k) {
        const std::size_t q = std::size_t(k) * ncol + col;
        if (qcw) qcw[q] = qcw[q] + rr(qccu(i, j, k));
        if (qiw) qiw[q] = qiw[q] + rr(qicu(i, j, k));
      }
    }
    if (clouds) {
      cal_cldfra1(nz, v(qv), v(qc), qi ? v(qi) : SD{nullptr, 0}, qs ? v(qs) : SD{nullptr, 0}, v(t), v(p), v(cf));
    } else {
      for (int k = 0; k < nz; ++k) cf[std::size_t(k) * ncol + col] = rr(0.0);
    }
    if (clouds && icloud_cu == 2) {  // KF: the convective cloud fraction over the resolved one
      for (int k = 0; k < nz; ++k) {
        const std::size_t q = std::size_t(k) * ncol + col;
        const rr cu = rr(cfdp(i, j, k) + cfsh(i, j, k));  // cldfra_cu
        cf[q] = (rr(1.) - cu) * cf[q];
        cf[q] = cf[q] + cu;
        cf[q] = std::min(rr(1.0), cf[q]);
        if (qcw) qcw[q] = qcw[q] + rr(qccu(i, j, k)) * cu;
        if (qiw) qiw[q] = qiw[q] + rr(qicu(i, j, k)) * cu;
      }
    }
    if (clouds && bl_clouds) {
      for (int k = 0; k < nz; ++k) {
        const std::size_t q = std::size_t(k) * ncol + col;
        const rr b = rr(cfbl(i, j, k));
        if (itimestep != 1) cf[q] = b;
        if (qcw[q] < rr(1.E-6) && b > rr(0.001)) qcw[q] = qcw[q] + rr(qcbl(i, j, k));
        if (qiw[q] < rr(1.E-8) && b > rr(0.001)) qiw[q] = qiw[q] + rr(qibl(i, j, k));
      }
    }
    if (bmj) {  // module_radiation_driver.F:1578-1591 (PCC scheme)
      for (int k = 0; k < nz; ++k) {
        const std::size_t q = std::size_t(k) * ncol + col;
        qcw[q] = qcw[q] + rr(qcconv(i, j, k));
        qiw[q] = qiw[q] + rr(qiconv(i, j, k));
        const rr cc = rr(ccldfra(i, j, k));
        const rr itrmx = c_min(cf[q], cc);
        const rr itrmn = c_max(rr(0.), cf[q] + cc - rr(1.));
        cf[q] = cf[q] + cc - rr(0.5) * (itrmx + itrmn);
      }
    }
    if (has_cldfra)
      for (int k = 0; k < nz; ++k) cldfra_out(i, j, k) = real(cf[std::size_t(k) * ncol + col]);
    if (has_diag && o3input == 2 && !o3_field)
      for (int k = 0; k < nz; ++k) o3_out(i, j, k) = real(o3[std::size_t(k) * ncol + col]);
  }
};
// Dudhia SWRAD per column (module_ra_sw.F) on the staged inputs: rsw = TTEN1D(NK)/pi3D, gsw; rswc = 0
struct KDudhia {
  int nx, nz, ncol;
  const rr *t, *p, *dz, *pi, *qv, *qc, *qr, *qi, *qs, *qg, *cz, *alb;
  rr solcon, cssca;
  int icloud;
  bool warm_rain;
  rr *rsw, *rswc, *gsw;
  int ci0 = 0, cj0 = 0;
  CPPWRF_HD void operator()(int i, int j) const {
    const int col = (i - ci0) + nx * (j - cj0);
    auto v = [&](const rr* a) { return a ? CSD(a + col, ncol) : CSD(); };
    gsw[col] = dudhia_swpara(nz, v(t), v(p), v(dz), v(qv), v(qc), v(qr), v(qi), v(qs), v(qg), rr(real(cz[col])),
                             alb[col], solcon, icloud, cssca, warm_rain, SD{rsw + col, ncol});
    for (int k = 0; k < nz; ++k) {
      const std::size_t q = std::size_t(k) * ncol + col;
      rsw[q] = rsw[q] / pi[q];
      rswc[q] = rr(0.0);
    }
  }
};
// RTHRATEN = LW (+ SW in daytime), GLW, GSW, SWDOWN = GSW/(1 - ALBEDO)
struct KRadOut {
  int nx, ncol;
  const rr *rlw, *rsw, *glw, *gsw, *cz, *olr, *lwd, *swd, *rlwc, *rswc, *lwcf, *swcf;
  View3<real> rth;
  View2<real> glw_out, gsw_out, swdown_out;
  View2<const real> alb;
  bool diag;  // RTHRATLW/SW and the LW/SW flux diagnostics
  View3<real> rthlw, rthsw, rthlwc, rthswc;
  View2<real> lw[9], sw[8], lwcf_out, swcf_out;
  const rr* o3lw = nullptr;  // o3input != 2: RRTMG_LWRAD's ozone into O3RAD
  View3<real> o3_out;
  const rr* sws = nullptr;  // SwBatchOutput::sfc into RadiationOutputs::swsfc
  View2<real> swsfc[10];
  bool dudhia = false;  // Dudhia SW: added at night too, no SW diagnostics; the direct / diffuse split from ht
  const rr* ht = nullptr;
  rr solcon = 0;
  int ci0 = 0, cj0 = 0;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    const int col = (i - ci0) + nx * (j - cj0);
    const std::size_t q = std::size_t(k) * ncol + col;
    const bool day = cz[col] > rr(0.0);
    rr r = rlw[q];
    if (dudhia || day) r = r + rsw[q];
    rth(i, j, k) = real(r);
    if (o3lw) o3_out(i, j, k) = real(o3lw[q]);
    if (diag) {
      rthlw(i, j, k) = real(rlw[q]);
      rthsw(i, j, k) = real(r) - real(rlw[q]);  // RTHRATENSW = RTHRATEN - RTHRATENLW
      rthlwc(i, j, k) = real(rlwc[q]);
      rthswc(i, j, k) = real(day && !dudhia ? rswc[q] : rr(0.0));
    }
    if (k == 0) {
      glw_out(i, j) = real(glw[col]);
      const real gs = real(day ? gsw[col] : rr(0.0));
      gsw_out(i, j) = gs;
      swdown_out(i, j) = real(rr(gs) / (rr(1.) - rr(alb(i, j))));
      if (diag) {
        for (int d = 0; d < 8; ++d) lw[d](i, j) = real(lwd[std::size_t(d) * ncol + col]);
        lw[8](i, j) = real(olr[col]);
        if (!dudhia)
          for (int d = 0; d < 8; ++d) sw[d](i, j) = real(day ? swd[std::size_t(d) * ncol + col] : rr(0.0));
        lwcf_out(i, j) = real(lwcf[col]);
        if (!dudhia) swcf_out(i, j) = real(day ? swcf[col] : rr(0.0));
      }
      if (sws) {  // RRTMG_SWRAD (0 at night): SWDDNI = SWDDIR / coszen, SWDDNIC = SWDDIRC / coszen
        for (int d = 0; d < 4; ++d) swsfc[d](i, j) = real(day ? sws[std::size_t(d) * ncol + col] : rr(0.0));
        const real dir = real(day ? sws[std::size_t(4) * ncol + col] : rr(0.0));
        const real dirc = real(day ? sws[std::size_t(6) * ncol + col] : rr(0.0));
        swsfc[4](i, j) = dir;
        swsfc[5](i, j) = real(day ? rr(dir) / cz[col] : rr(0.0));
        swsfc[6](i, j) = real(day ? sws[std::size_t(5) * ncol + col] : rr(0.0));
        swsfc[7](i, j) = dirc;
        swsfc[8](i, j) = real(day ? rr(dirc) / cz[col] : rr(0.0));
        const real dif = swsfc[6](i, j), sd = swdown_out(i, j);
        swsfc[9](i, j) = rr(sd) > rr(0.001) ? real(std::min(rr(dif) / rr(sd), rr(1.0))) : real(0);
      }
      if (dudhia && swsfc[0].p) {  // radiation_driver's direct / diffuse split (0 where coszen <= 1e-3)
        const real czr = real(cz[col]), sd = swdown_out(i, j);
        rr dif = 0, dir = 0;
        swsfc[4](i, j) = real(0);
        swsfc[5](i, j) = real(0);
        swsfc[6](i, j) = real(0);
        if (sw_diffuse_backup(rr(czr), rr(sd), solcon, ht[col], dif, dir)) {
          swsfc[6](i, j) = real(dif);
          swsfc[4](i, j) = real(dir);
          swsfc[5](i, j) = real(rr(swsfc[4](i, j)) / std::max(rr(czr), rr(1e-4)));
        }
        const real dsw = swsfc[6](i, j);
        swsfc[9](i, j) = rr(sd) > rr(0.001) ? real(std::min(rr(dsw) / rr(sd), rr(1.0))) : real(0);
      }
    }
  }
};
#endif

}  // namespace

struct RadiationDriver::Impl {
  std::unique_ptr<RrtmgLwTables> lw;
  std::unique_ptr<RrtmLwTables> rrtm;  // RRTM longwave (ra_lw_physics = 1)
  int rrtm_nlayers = 0;
  std::unique_ptr<RrtmgSwTables> sw;
  CamGases gases;
  CamOzone ozone;
  std::vector<rr> xlat, xlong, ozmixm;  // ozmixm: [the worker's column][month][level] (oznini)
#if defined(CPPWRF_USE_CUDA)
  std::unique_ptr<RrtmgLwBatch> lwb;
  std::unique_ptr<RrtmgSwBatch> swb;
  std::vector<void*> bufs;
  rr* dev(std::size_t n) {
    void* q = dev_alloc(n * sizeof(rr));
    bufs.push_back(q);
    return static_cast<rr*>(q);
  }
  rr *xlat_d = nullptr, *xlong_d = nullptr, *oz_d = nullptr, *pin_d = nullptr;
  // staged inputs, outputs (column fastest)
  rr *dz = nullptr;
  rr *p8w, *t8w, *p3d, *pi3d, *t3d, *qv, *qc, *qi, *qs, *rec, *rei, *res, *cf, *o3;
  rr *qr = nullptr, *qg = nullptr, *ht = nullptr;  // Dudhia SW
  rr *tsk, *emiss, *xland, *xice, *snow, *albedo, *cz;
  rr *rlw, *rsw, *glw, *gsw, *olr, *lwd, *swd, *rlwc, *rswc, *lwcf, *swcf, *sws;
  ~Impl() {
    for (void* q : bufs) dev_free(q);
  }
#endif
};

RadiationDriver::RadiationDriver(const Grid& g, const RadiationOptions& o, const std::vector<rr>& xlat,
                                 const std::vector<rr>& xlong)
    : p_(std::make_unique<Impl>()), g_(g), o_(o) {
  if ((o.ra_lw_physics != 4 && o.ra_lw_physics != 1) || (o.ra_sw_physics != 4 && o.ra_sw_physics != 1))
    throw std::runtime_error("radiation: only ra_lw_physics = 4 (RRTMG) or 1 (RRTM) with ra_sw_physics = 4 (RRTMG) or "
                             "1 (Dudhia) are ported");
  if (o.aer_opt != 0) throw std::runtime_error("radiation: aer_opt != 0 is not ported");
  if (o.o3input != 2 && o.o3input != 0) throw std::runtime_error("radiation: o3input must be 0 or 2");
  // cldovrlp 0 (clear sky only) and the McICA overlaps 1-5 (mcica_subcol_lw/_sw stop for any other icld)
  if (o.cldovrlp < 0 || o.cldovrlp > 5) throw std::runtime_error("radiation: cldovrlp must be 0 to 5");
  if (o.icloud_cu < 0 || o.icloud_cu > 2) throw std::runtime_error("radiation: icloud_cu must be 0, 1 or 2");
  if (o.bmj_feedback && o.icloud_cu != 0) throw std::runtime_error("radiation: bmj_rad_feedback with icloud_cu");
  // wrf.exe writes NaN cloud fractions with the Dudhia shortwave and the Grell clouds (V1 inputs): not ported
  if (o.ra_sw_physics == 1 && o.icloud_cu == 1)
    throw std::runtime_error("radiation: cu_rad_feedback with a Grell scheme and the Dudhia shortwave is not supported");
  const bool rrtmg_lw = o.ra_lw_physics == 4;
  if (rrtmg_lw) {  // the LW column arrays hold RRTMG_LW_SWEEP_MAX layers: the model's and the buffer layers up to 0 hPa
    RrtmgLwConfig lc;
    lc.p_top = o.p_top;
    if (rrtmg_lw_call_config(g.nz, lc, false, false, false).nlayers > RRTMG_LW_SWEEP_MAX)
      throw std::runtime_error("radiation: RRTMG LW layers (nz + nint(p_top/400 Pa)) exceed " +
                               std::to_string(RRTMG_LW_SWEEP_MAX));
  }
  stepra_ = o.radt == 0 ? 1 : std::max(1, int(std::lround(o.radt * rr(60.) / g.dt)));
  if (rrtmg_lw) p_->lw = std::make_unique<RrtmgLwTables>();
  const bool rrtmg_sw = o.ra_sw_physics == 4;
  if (rrtmg_sw) p_->sw = std::make_unique<RrtmgSwTables>();
  // RRTMG tables: phys.table_dir, then CPPWRF's copy (data/rrtmg of the data directory, core/data_dir.hpp)
  auto table_file = [&](const std::string& name) {
    const std::string a = o.table_dir + "/" + name;
    if (std::ifstream(a).good()) return a;
    const std::string b = data_dir() + "/rrtmg/" + name;
    if (std::ifstream(b).good()) return b;
    throw std::runtime_error("radiation: cannot find " + name + " in " + o.table_dir + " or the data directory");
  };
  if (rrtmg_lw && !load_rrtmg_lw_tables(table_file("rrtmg_lw_tables.bin"), *p_->lw))
    throw std::runtime_error("radiation: cannot read rrtmg_lw_tables.bin");
  if (!rrtmg_lw) {  // rrtminit (phy_init): rrtm_lookuptable reads RRTM_DATA of the run directory
    p_->rrtm = std::make_unique<RrtmLwTables>();
    std::string err;
    if (!load_rrtm_lw_tables(o.table_dir + "/RRTM_DATA", *p_->rrtm, err)) throw std::runtime_error("radiation: " + err);
    p_->rrtm_nlayers = rrtm_lw_nlayers(g.nz, o.p_top);
  }
  if (rrtmg_sw && !load_rrtmg_sw_tables(table_file("rrtmg_sw_tables.bin"), *p_->sw))
    throw std::runtime_error("radiation: cannot read rrtmg_sw_tables.bin");
  // the SW column arrays (CPU and GPU) hold RRTMG_SW_COL_MAX layers: the model's and the one above the top (WRF:
  // mxlay = 203)
  if (rrtmg_sw && g.nz + 1 > RRTMG_SW_COL_MAX)
    throw std::runtime_error("radiation: RRTMG SW layers (nz + 1) exceed " + std::to_string(RRTMG_SW_COL_MAX));
  if (o.ghg_input == 1) p_->gases.load(o.table_dir + "/CAMtr_volume_mixing_ratio");
  p_->xlat = xlat;
  p_->xlong = xlong;
  if (o.o3input == 2) {  // oznini (phy_init, domain 1 only)
    // the columns the worker owns (a decomposed run: its patch), row by row
    const tile::PatchBox ob = tile::own_box(g);
    const std::size_t onx = std::size_t(ob.i1 - ob.i0), nc = onx * std::size_t(ob.j1 - ob.j0);
    const std::size_t per = std::size_t(CamOzone::nmon) * CamOzone::levsiz;
    p_->ozmixm.assign(nc * per, rr(0.0));
    if (o.oznini) {
      p_->ozone.load(o.table_dir);
      for (int j = ob.j0; j < ob.j1; ++j)
        for (int i = ob.i0; i < ob.i1; ++i)
          p_->ozone.column(xlat[std::size_t(i) + std::size_t(g.nx) * j],
                           &p_->ozmixm[(std::size_t(i - ob.i0) + onx * std::size_t(j - ob.j0)) * per]);
    } else {
      p_->ozone.pin.assign(CamOzone::levsiz, rr(0.0));
    }
  }
#if defined(CPPWRF_USE_CUDA)
  {
    Impl& q = *p_;
    // the columns the worker owns (a decomposed run: its patch), row by row
    const tile::PatchBox ob = tile::own_box(g);
    const int cnx = ob.i1 - ob.i0;
    const std::size_t nc = std::size_t(cnx) * std::size_t(ob.j1 - ob.j0), n3 = nc * g.nz, n3w = nc * (g.nz + 1);
    auto up = [&](const std::vector<rr>& h) {
      rr* d = q.dev(std::max<std::size_t>(h.size(), 1));
      if (!h.empty()) copy_to_dev(d, h.data(), h.size() * sizeof(rr));
      return d;
    };
    auto own = [&](const std::vector<rr>& h, std::size_t per) {  // per values per column of the domain
      if (h.empty() || !g.patched()) return up(h);
      std::vector<rr> l(nc * per);
      for (int j = ob.j0; j < ob.j1; ++j)
        for (int i = ob.i0; i < ob.i1; ++i)
          std::copy_n(&h[(std::size_t(i) + std::size_t(g.nx) * j) * per], per,
                      &l[(std::size_t(i - ob.i0) + std::size_t(cnx) * (j - ob.j0)) * per]);
      return up(l);
    };
    q.xlat_d = own(xlat, 1);
    q.xlong_d = own(xlong, 1);
    q.oz_d = up(q.ozmixm);  // the worker's columns already
    q.pin_d = up(q.ozone.pin);
    q.p8w = q.dev(n3w); q.t8w = q.dev(n3w);
    for (rr** a : {&q.p3d, &q.pi3d, &q.t3d, &q.qv, &q.qc, &q.qi, &q.qs, &q.rec, &q.rei, &q.res, &q.cf, &q.o3, &q.dz,
                       &q.rlw, &q.rsw, &q.rlwc, &q.rswc})
      *a = q.dev(n3);
    for (rr** a : {&q.tsk, &q.emiss, &q.xland, &q.xice, &q.snow, &q.albedo, &q.cz, &q.glw, &q.gsw, &q.olr, &q.lwcf, &q.swcf}) *a = q.dev(nc);
    q.lwd = q.dev(8 * nc);
    q.swd = q.dev(8 * nc);
    q.sws = q.dev(7 * nc);
    RrtmgLwConfig lc;
    lc.p_top = o.p_top;
    if (rrtmg_lw)
      q.lwb = std::make_unique<RrtmgLwBatch>(*q.lw, rrtmg_lw_call_config(g.nz, lc, false, false, false).nlayers,
                                             int(nc));
    if (rrtmg_sw) q.swb = std::make_unique<RrtmgSwBatch>(*q.sw, g.nz + 1, int(nc));
    for (rr** a : {&q.qr, &q.qg}) *a = q.dev(n3);
    q.ht = q.dev(nc);
  }
#endif
}

RadiationDriver::~RadiationDriver() = default;

bool RadiationDriver::run(long itimestep, rr xtime, rr julian, int yr, const RadiationInputs& in,
                          const RadiationOutputs& out, bool adapt_step_flag, real curr_secs, real* radtacttime) {
  // run decision (module_radiation_driver.F:1083-1138); the adaptive time step's first call sets the activation time
  const bool doing_adapt_dt = adapt_step_flag;
  if (doing_adapt_dt && *radtacttime == real(rr(0.))) *radtacttime = curr_secs + real(o_.radt) * real(rr(60.));
  bool run_param = itimestep == 1 || o_.radt == rr(0.) || stepra_ == 1;
  if (!run_param && !doing_adapt_dt) run_param = int(itimestep % stepra_) == 1 + o_.ra_call_offset;
  if (!run_param && doing_adapt_dt && curr_secs >= *radtacttime) {
    run_param = true;
    *radtacttime = curr_secs + real(o_.radt) * real(60);
  }
  if (!run_param) return false;
  rr declin, solcon;
  radconst(julian, declin, solcon);  // Solar_step (radconst before calc_coszen)
  declin_ = declin;
  solcon_ = solcon;
  RrtmgLwConfig lc;
  lc.icloud = o_.icloud; lc.cldovrlp = o_.cldovrlp; lc.idcor = o_.idcor; lc.p_top = o_.p_top;
  lc.yr = yr; lc.julian = julian; lc.o3input = o_.o3input;
  if (o_.ghg_input == 1) p_->gases.interpolate(yr, julian, lc.co2, lc.n2o, lc.ch4, lc.cfc11, lc.cfc12);
  else rrtmg_default_gases(lc);
  RrtmgSwConfig sc;
  sc.icloud = o_.icloud; sc.cldovrlp = o_.cldovrlp; sc.idcor = o_.idcor; sc.p_top = o_.p_top; sc.yr = yr;
  sc.julian = julian; sc.solcon = solcon; sc.o3input = o_.o3input; sc.co2 = lc.co2; sc.n2o = lc.n2o; sc.ch4 = lc.ch4;
  int nm = 0, np = 0;
  rr f1 = 0, f2 = 0;
  if (o_.o3input == 2) ozone_time_weights(julian, nm, np, f1, f2);
  const std::size_t per = std::size_t(CamOzone::nmon) * CamOzone::levsiz;
#if defined(CPPWRF_USE_CUDA)
  {
    (void)per;
    Impl& q = *p_;
    const int nz = g_.nz;
    const tile::PatchBox ob = tile::own_box(g_);  // the worker's columns (the domain without a decomposition)
    const int cnx = ob.i1 - ob.i0, ncol = cnx * (ob.j1 - ob.j0);
    const bool dudhia = o_.ra_sw_physics == 1, rrtm = o_.ra_lw_physics == 1;
    if (o_.icloud != 0 && o_.icloud != 1) throw std::runtime_error("radiation: only icloud = 0 or 1 are ported");
    // bmj_rad_feedback adds QCCONV / QICONV to qc / qi: WRF's moist(P_QC) / moist(P_QI) without F_QC / F_QI is its
    // scalar placeholder, not ported
    if (o_.bmj_feedback && (!in.qc || !in.qi))
      throw std::runtime_error("radiation: bmj_rad_feedback needs a microphysics scheme with cloud water and ice");
    if (in.qc && !in.qi && in.qs) throw std::runtime_error("radiation: F_QC, F_QS without F_QI is not ported");
    const bool clouds = o_.icloud != 0 && in.qc != nullptr;
    KStageAll ks{};
    ks.nx = cnx; ks.nz = nz; ks.ncol = ncol; ks.n3 = 0;
    ks.ci0 = ob.i0; ks.cj0 = ob.j0;
    ks.f3w[0] = in.p_hyd_w->cview(); ks.d3w[0] = q.p8w;
    ks.f3w[1] = in.t8w->cview();     ks.d3w[1] = q.t8w;
    auto add3 = [&](const Field3* f, rr* d) {
      if (f) {
        ks.f3[ks.n3] = f->cview();
        ks.d3[ks.n3] = d;
        ++ks.n3;
      }
    };
    add3(in.p_hyd, q.p3d); add3(in.pi_phy, q.pi3d); add3(in.t_phy, q.t3d); add3(in.qv, q.qv);
    add3(in.qc, q.qc); add3(in.qi, q.qi); add3(in.qs, q.qs);
    add3(in.re_cloud, q.rec); add3(in.re_ice, q.rei); add3(in.re_snow, q.res);
    add3(in.dz8w, q.dz);
    // RRTMG's "MP option 3" (WSM3: qc without qi): below 273.15 K the cloud water is ice and the rain snow
    const bool split3 = in.qc && !in.qi && !o_.warm_rain;
    if (dudhia || rrtm || split3) add3(in.qr, q.qr);
    if (dudhia) {
      add3(in.qg, q.qg);
      ks.ht = in.ht->view();
      ks.dht = q.ht;
    }
    ks.f2[0] = in.tsk->cview();    ks.d2[0] = q.tsk;
    ks.f2[1] = in.emiss->cview();  ks.d2[1] = q.emiss;
    ks.f2[2] = in.xland->cview();  ks.d2[2] = q.xland;
    ks.f2[3] = in.xice->cview();   ks.d2[3] = q.xice;
    ks.f2[4] = in.snow->cview();   ks.d2[4] = q.snow;
    ks.f2[5] = in.albedo->cview(); ks.d2[5] = q.albedo;
    parallel_for(Range3D{ob.i0, ob.i1, ob.j0, ob.j1, 0, nz + 1}, ks);
    KRadColumn kc{};
    kc.nx = cnx; kc.nz = nz; kc.ncol = ncol;
    kc.ci0 = ob.i0; kc.cj0 = ob.j0; kc.gnx = g_.nx;
    kc.julian = julian; kc.xtime_mid = xtime + o_.radt * rr(0.5); kc.gmt = o_.gmt; kc.declin = declin;
    kc.xlat = q.xlat_d; kc.xlong = q.xlong_d;
    kc.o3input = o_.o3input; kc.ozmixm = q.oz_d; kc.pin = q.pin_d; kc.nm = nm; kc.np = np; kc.f1 = f1; kc.f2 = f2;
    kc.clouds = clouds;
    kc.p = q.p3d; kc.t = q.t3d; kc.qv = q.qv; kc.qc = in.qc ? q.qc : nullptr; kc.qi = in.qi ? q.qi : nullptr;
    kc.qs = in.qs ? q.qs : nullptr;
    kc.o3 = q.o3; kc.cf = q.cf; kc.cz = q.cz;
    kc.has_cldfra = out.cldfra != nullptr;
    if (out.cldfra) kc.cldfra_out = out.cldfra->view();
    kc.coszen_out = out.coszen->view();
    kc.has_diag = out.hrang != nullptr;
    if (kc.has_diag) { kc.hrang_out = out.hrang->view(); kc.o3_out = out.o3rad->view(); }
    kc.o3_field = !o_.oznini && out.o3rad;
    if (kc.o3_field) kc.o3_in = out.o3rad->cview();
    kc.bl_clouds = in.cldfra_bl != nullptr && in.qc && in.qi;
    kc.itimestep = itimestep;
    kc.qcw = in.qc ? q.qc : nullptr; kc.qiw = in.qi ? q.qi : nullptr;
    if (kc.bl_clouds) {
      kc.cfbl = in.cldfra_bl->cview(); kc.qcbl = in.qc_bl->cview(); kc.qibl = in.qi_bl->cview();
    }
    kc.icloud_cu = o_.icloud_cu;
    if (o_.icloud_cu != 0) {
      kc.qccu = in.qc_cu->cview(); kc.qicu = in.qi_cu->cview();
      kc.cfdp = in.cldfra_dp->cview(); kc.cfsh = in.cldfra_sh->cview();
    }
    kc.bmj = o_.bmj_feedback;
    if (kc.bmj) { kc.qcconv = in.qcconv->cview(); kc.qiconv = in.qiconv->cview(); kc.ccldfra = in.ccldfra->cview(); }
    column_for(Range2D{ob.i0, ob.i1, ob.j0, ob.j1}, kc);
    LwBatchInput li{};
    li.ncol = ncol;
    li.col_begin = 0;
    li.col_end = ncol;
    li.tsk = q.tsk; li.emiss = q.emiss; li.xland = q.xland; li.xice = q.xice; li.snow = q.snow;
    li.p8w = q.p8w; li.t8w = q.t8w; li.p3d = q.p3d; li.pi3d = q.pi3d; li.t3d = q.t3d; li.qv = q.qv; li.cldfra = q.cf;
    li.o3 = q.o3; li.xlat = q.xlat_d; li.dz8w = q.dz;
    li.qc = in.qc ? q.qc : nullptr; li.qi = in.qi ? q.qi : nullptr; li.qs = in.qs ? q.qs : nullptr;
    li.re_cloud = in.re_cloud ? q.rec : nullptr; li.re_ice = in.re_ice ? q.rei : nullptr;
    li.re_snow = in.re_snow ? q.res : nullptr;
    li.split3 = split3; li.qr = in.qr ? q.qr : nullptr;
    if (!in.qr) li.split3 = false;
    LwBatchOutput lo{q.rlw, q.glw, q.olr, q.lwd, q.rlwc, q.lwcf};
    // O3RAD from RRTMG_LWRAD (o33d is not read); RRTM leaves it to the driver's CAM ozone (o3input = 2)
    const bool o3lw = !rrtm && o_.o3input != 2 && out.o3rad != nullptr;
    if (o3lw) lo.o3 = q.o3;
    if (rrtm) {
#ifdef CPPWRF_DOUBLE
      // RRTMLWRAD on the host from the staged columns (cloud fractions, BL and cumulus clouds included); no LW flux
      // diagnostics and LWCF (RRTMG only): 0
      const std::size_t n3 = std::size_t(ncol) * nz, n3w = std::size_t(ncol) * (nz + 1), nc = std::size_t(ncol);
      auto host = [](const rr* d, std::size_t n) {
        std::vector<rr> h(n);
        copy_to_host(h.data(), d, n * sizeof(rr));
        return h;
      };
      const std::vector<rr> hp8w = host(q.p8w, n3w), ht8w = host(q.t8w, n3w), hp = host(q.p3d, n3),
                                ht3 = host(q.t3d, n3), hdz = host(q.dz, n3), hpi = host(q.pi3d, n3),
                                hqv = host(q.qv, n3), hcf = host(q.cf, n3), htsk = host(q.tsk, nc),
                                hem = host(q.emiss, nc);
      const std::vector<rr> hqc = in.qc ? host(q.qc, n3) : std::vector<rr>(),
                                hqr = in.qr ? host(q.qr, n3) : std::vector<rr>(),
                                hqi = in.qi ? host(q.qi, n3) : std::vector<rr>(),
                                hqs = in.qs ? host(q.qs, n3) : std::vector<rr>();
      const RrtmLwBatchInput bi{ncol, nz, hp8w.data(), ht8w.data(), hp.data(), ht3.data(), hdz.data(), hpi.data(),
                                hqv.data(), in.qc ? hqc.data() : nullptr, in.qr ? hqr.data() : nullptr,
                                in.qi ? hqi.data() : nullptr, in.qs ? hqs.data() : nullptr, hcf.data(), htsk.data(),
                                hem.data(), in.qi != nullptr, lc.co2, lc.n2o, lc.ch4, cst::r_d, cst::g, o_.icloud,
                                o_.warm_rain, li.col_begin, li.col_end};
      std::vector<rr> hr(n3), hrc(n3), hg(nc), ho(nc), zero(8 * nc, rr(0.0));
      rrtm_lw_batch(*q.rrtm, bi, q.rrtm_nlayers, hr.data(), hrc.data(), hg.data(), ho.data());
      copy_to_dev(q.rlw, hr.data(), n3 * sizeof(rr));
      copy_to_dev(q.rlwc, hrc.data(), n3 * sizeof(rr));
      copy_to_dev(q.glw, hg.data(), nc * sizeof(rr));
      copy_to_dev(q.olr, ho.data(), nc * sizeof(rr));
      copy_to_dev(q.lwd, zero.data(), 8 * nc * sizeof(rr));
      copy_to_dev(q.lwcf, zero.data(), nc * sizeof(rr));
#else
      // FP32 build: RRTM LW computes in double (rrtm_lw.cpp is not converted); the staged columns are widened and
      // the tendencies and fluxes rounded once
      const std::size_t n3 = std::size_t(ncol) * nz, n3w = std::size_t(ncol) * (nz + 1), nc = std::size_t(ncol);
      auto host = [](const rr* d, std::size_t n) {
        std::vector<rr> h(n);
        copy_to_host(h.data(), d, n * sizeof(rr));
        return std::vector<double>(h.begin(), h.end());
      };
      auto dev = [](rr* d, const std::vector<double>& h) {
        const std::vector<rr> f(h.begin(), h.end());
        copy_to_dev(d, f.data(), f.size() * sizeof(rr));
      };
      const std::vector<double> hp8w = host(q.p8w, n3w), ht8w = host(q.t8w, n3w), hp = host(q.p3d, n3),
                                ht3 = host(q.t3d, n3), hdz = host(q.dz, n3), hpi = host(q.pi3d, n3),
                                hqv = host(q.qv, n3), hcf = host(q.cf, n3), htsk = host(q.tsk, nc),
                                hem = host(q.emiss, nc);
      const std::vector<double> hqc = in.qc ? host(q.qc, n3) : std::vector<double>(),
                                hqr = in.qr ? host(q.qr, n3) : std::vector<double>(),
                                hqi = in.qi ? host(q.qi, n3) : std::vector<double>(),
                                hqs = in.qs ? host(q.qs, n3) : std::vector<double>();
      const RrtmLwBatchInput bi{ncol, nz, hp8w.data(), ht8w.data(), hp.data(), ht3.data(), hdz.data(), hpi.data(),
                                hqv.data(), in.qc ? hqc.data() : nullptr, in.qr ? hqr.data() : nullptr,
                                in.qi ? hqi.data() : nullptr, in.qs ? hqs.data() : nullptr, hcf.data(), htsk.data(),
                                hem.data(), in.qi != nullptr, lc.co2, lc.n2o, lc.ch4, cst::r_d, cst::g, o_.icloud,
                                o_.warm_rain, li.col_begin, li.col_end};
      std::vector<double> hr(n3), hrc(n3), hg(nc), ho(nc);
      rrtm_lw_batch(*q.rrtm, bi, q.rrtm_nlayers, hr.data(), hrc.data(), hg.data(), ho.data());
      dev(q.rlw, hr);
      dev(q.rlwc, hrc);
      dev(q.glw, hg);
      dev(q.olr, ho);
      const std::vector<rr> zero(8 * nc, rr(0.0));
      copy_to_dev(q.lwd, zero.data(), 8 * nc * sizeof(rr));
      copy_to_dev(q.lwcf, zero.data(), nc * sizeof(rr));
#endif
    } else {
      q.lwb->run(rrtmg_lw_call_config(nz, lc, in.re_cloud, in.re_ice, in.re_snow), li, lo);
    }
    if (dudhia) {
      KDudhia kd{cnx, nz, ncol, q.t3d, q.p3d, q.dz, q.pi3d, q.qv, in.qc ? q.qc : nullptr, in.qr ? q.qr : nullptr,
                 in.qi ? q.qi : nullptr, in.qs ? q.qs : nullptr, in.qg ? q.qg : nullptr, q.cz, q.albedo, solcon,
                 dudhia_cssca(o_.swrad_scat), o_.icloud, o_.warm_rain, q.rsw, q.rswc, q.gsw};
      kd.ci0 = ob.i0; kd.cj0 = ob.j0;
      column_for(Range2D{ob.i0, ob.i1, ob.j0, ob.j1}, kd);
    } else {
      SwBatchInput si{};
      si.ncol = ncol; si.col_begin = li.col_begin; si.col_end = li.col_end; si.coszen = q.cz; si.albedo = q.albedo; si.xland = q.xland; si.xice = q.xice; si.snow = q.snow;
      si.p8w = q.p8w; si.t8w = q.t8w; si.p3d = q.p3d; si.pi3d = q.pi3d; si.t3d = q.t3d; si.qv = q.qv; si.cldfra = q.cf;
      si.o3 = q.o3; si.xlat = q.xlat_d; si.dz8w = q.dz; si.qc = li.qc; si.qi = li.qi; si.qs = li.qs;
      si.split3 = li.split3; si.qr = li.qr;
      si.re_cloud = li.re_cloud; si.re_ice = li.re_ice; si.re_snow = li.re_snow;
      q.swb->run(rrtmg_sw_call_config(nz, sc), si,
                 SwBatchOutput{q.rsw, q.gsw, q.swd, q.rswc, q.swcf, out.swsfc[0] ? q.sws : nullptr});
    }
    KRadOut ko{cnx, ncol, q.rlw, q.rsw, q.glw, q.gsw, q.cz, q.olr, q.lwd, q.swd, q.rlwc, q.rswc, q.lwcf, q.swcf,
               out.rthraten->view(), out.glw->view(), out.gsw->view(), out.swdown->view(), in.albedo->cview(),
               out.rthratlw != nullptr};
    if (o3lw) {
      ko.o3lw = q.o3;
      ko.o3_out = out.o3rad->view();
    }
    if (ko.diag) {
      ko.rthlwc = out.rthratlwc->view();
      ko.rthswc = out.rthratswc->view();
      ko.lwcf_out = out.lwcf->view();
      ko.swcf_out = out.swcf->view();
      ko.rthlw = out.rthratlw->view();
      ko.rthsw = out.rthratsw->view();
      for (int d = 0; d < 9; ++d) ko.lw[d] = out.lw[d]->view();
      for (int d = 0; d < 8; ++d) ko.sw[d] = out.sw[d]->view();
    }
    if (out.swsfc[0]) {
      if (!dudhia) ko.sws = q.sws;
      for (int d = 0; d < 10; ++d) ko.swsfc[d] = out.swsfc[d]->view();
    }
    ko.dudhia = dudhia;
    ko.ht = q.ht;
    ko.solcon = solcon;
    ko.ci0 = ob.i0; ko.cj0 = ob.j0;
    parallel_for(Range3D{ob.i0, ob.i1, ob.j0, ob.j1, 0, nz}, ko);
    return true;
  }
#else
  const int nx = g_.nx, ny = g_.ny, nz = g_.nz;
  // inputs to the host
  std::vector<real> t, p, pw, pi, rho, dz, tw, qv, qc, qi, qs, rec, rei, res, tsk, alb, ems, xland, xice, snow, cfbl,
      qcbl, qibl, qccu, qicu, cfdp, cfsh, qr, qg, ht;
  in.t_phy->download(t); in.p_hyd->download(p); in.p_hyd_w->download(pw); in.pi_phy->download(pi);
  in.rho->download(rho); in.dz8w->download(dz); in.t8w->download(tw); in.qv->download(qv);
  if (in.qc) in.qc->download(qc);
  if (in.qi) in.qi->download(qi);
  if (in.qs) in.qs->download(qs);
  const bool dudhia = o_.ra_sw_physics == 1, rrtm = o_.ra_lw_physics == 1;
  if ((dudhia || rrtm || (in.qc && !in.qi && !o_.warm_rain)) && in.qr) in.qr->download(qr);
  if (dudhia) {
    if (in.qg) in.qg->download(qg);
    in.ht->download(ht);
  }
  if (in.re_cloud) in.re_cloud->download(rec);
  if (in.re_ice) in.re_ice->download(rei);
  if (in.re_snow) in.re_snow->download(res);
  const bool bl_clouds = in.cldfra_bl != nullptr;
  if (bl_clouds) { in.cldfra_bl->download(cfbl); in.qc_bl->download(qcbl); in.qi_bl->download(qibl); }
  if (o_.icloud_cu != 0) {
    in.qc_cu->download(qccu); in.qi_cu->download(qicu); in.cldfra_dp->download(cfdp); in.cldfra_sh->download(cfsh);
  }
  std::vector<real> qcconv, qiconv, ccldfra;
  if (o_.bmj_feedback) { in.qcconv->download(qcconv); in.qiconv->download(qiconv); in.ccldfra->download(ccldfra); }
  // icloud = 1 with F_QC or F_QI: cal_cldfra1 (else CLDFRA = 0); qc is required with qi (WSM/Kessler)
  const bool clouds = o_.icloud != 0 && in.qc != nullptr;
  if (o_.icloud != 0 && o_.icloud != 1) throw std::runtime_error("radiation: only icloud = 0 or 1 are ported");
  // bmj_rad_feedback adds QCCONV / QICONV to qc / qi: WRF's moist(P_QC) / moist(P_QI) without F_QC / F_QI is its
  // scalar placeholder, not ported
  if (o_.bmj_feedback && (!in.qc || !in.qi))
    throw std::runtime_error("radiation: bmj_rad_feedback needs a microphysics scheme with cloud water and ice");
  if (in.qc && !in.qi && in.qs) throw std::runtime_error("radiation: F_QC, F_QS without F_QI is not ported");
  in.tsk->download(tsk); in.albedo->download(alb); in.emiss->download(ems); in.xland->download(xland);
  in.xice->download(xice); in.snow->download(snow);
  const Layout &Lm = in.t_phy->layout(), &Lw = in.p_hyd_w->layout(), &Lt8 = in.t8w->layout(),
               &L2 = in.tsk->layout();
  std::vector<real> rth(out.rthraten->layout().total, real(0)), glw(out.glw->layout().total, real(0)),
      gsw(out.gsw->layout().total, real(0)), swdown(out.swdown->layout().total, real(0)),
      coszen(out.coszen->layout().total, real(0)), cldfra(out.cldfra ? out.cldfra->layout().total : 0, real(0));
  const Layout &Lo3 = out.rthraten->layout(), &Lo2 = out.glw->layout();
  const bool diag = out.rthratlw != nullptr;
  const bool o3nest = !o_.oznini && out.o3rad;  // a nest: o3rad as the parent forced it
  std::vector<real> o3h;
  if (o3nest) out.o3rad->download(o3h);
  std::vector<real> rthlw(diag ? Lo3.total : 0, real(0)), rthsw(diag ? Lo3.total : 0, real(0)),
      rthlwc(diag ? Lo3.total : 0, real(0)), rthswc(diag ? Lo3.total : 0, real(0)),
      o3r(diag ? Lo3.total : 0, real(0));
  std::vector<real> lwcf2(diag ? Lo2.total : 0, real(0)), swcf2(diag ? Lo2.total : 0, real(0)),
      hrang2(diag ? Lo2.total : 0, real(0));
  std::vector<std::vector<real>> lwd(diag ? 9 : 0, std::vector<real>(Lo2.total, real(0))),
      swd(diag ? 8 : 0, std::vector<real>(Lo2.total, real(0))),
      swsfc(out.swsfc[0] ? 10 : 0, std::vector<real>(Lo2.total, real(0)));

  Range2D own{0, nx, 0, ny};
  tile::clip(own);  // a tile: its own columns (the others' inputs are not set)
  const int own_nx = own.i1 - own.i0, own_cols = own_nx * (own.j1 - own.j0);
#pragma omp parallel
  {
    auto lws = std::make_unique<RrtmgLwWorkspace>();
    std::unique_ptr<RrtmLwWorkspace> rws;
    if (rrtm) rws = std::make_unique<RrtmLwWorkspace>(nz, p_->rrtm_nlayers);
    std::string rmsg;
    auto sws = std::make_unique<RrtmgSwWorkspace>();
    std::vector<rr> cp8w(nz + 1), ct8w(nz + 1), cp(nz), cpi(nz), cdz(nz), ct(nz), crho(nz), cqv(nz),
        zero(nz, rr(0.0)), ccf(nz, rr(0.0)), co3(nz, rr(0.0)), rlw(nz), rlwc(nz), rsw(nz), rswc(nz), cqc(nz, rr(0.0)),
        cqi(nz, rr(0.0)), cqs(nz, rr(0.0)), crec(nz), crei(nz), cres(nz), cqr(nz, rr(0.0)), cqg(nz, rr(0.0)),
        sqc(nz), sqi(nz), sqs(nz);
    // RRTMG_LWRAD / RRTMG_SWRAD "for MP option 3" (no F_QI, not warm_rain): QI1D = QC1D, QS1D = QR1D, QC1D = 0
    // below 273.15 K, on copies (cqc stays as cal_cldfra1 and Dudhia see it)
    const bool split3 = in.qc && !in.qi && !o_.warm_rain;
    auto split_mp3 = [&] {
      for (int k = 0; k < nz; ++k) {
        sqc[k] = cqc[k]; sqi[k] = cqi[k]; sqs[k] = cqs[k];
        if (split3 && ct[k] < rr(273.15)) {
          sqi[k] = std::max(rr(0.0), cqc[k]);
          sqs[k] = std::max(rr(0.0), cqr[k]);
          sqc[k] = rr(0.0);
        }
      }
    };
#pragma omp for schedule(dynamic, 16)
    for (int oc = 0; oc < own_cols; ++oc) {
      const int i = own.i0 + oc % own_nx, j = own.j0 + oc / own_nx, c = i + nx * j;
      for (int k = 0; k <= nz; ++k) { cp8w[k] = pw[Lw.index(i, j, k)]; ct8w[k] = tw[Lt8.index(i, j, k)]; }
      for (int k = 0; k < nz; ++k) {
        const std::size_t q = Lm.index(i, j, k);
        cp[k] = p[q]; cpi[k] = pi[q]; cdz[k] = dz[q]; ct[k] = t[q]; crho[k] = rho[q]; cqv[k] = qv[q];
        if (!qc.empty()) cqc[k] = qc[q];
        if (!qi.empty()) cqi[k] = qi[q];
        if (!qs.empty()) cqs[k] = qs[q];
        if (!rec.empty()) crec[k] = rec[q];
        if (!rei.empty()) crei[k] = rei[q];
        if (!res.empty()) cres[k] = res[q];
        if (!qr.empty()) cqr[k] = qr[q];
        if (!qg.empty()) cqg[k] = qg[q];
      }
      if (o_.icloud_cu == 1) {  // Grell schemes (module_radiation_driver.F:1264-1281): qc + qc_cu, qi + qi_cu
        for (int k = 0; k < nz; ++k) {
          const std::size_t q = Lm.index(i, j, k);
          if (!qc.empty()) cqc[k] = cqc[k] + qccu[q];
          if (!qi.empty()) cqi[k] = cqi[k] + qicu[q];
        }
      }
      if (clouds)
        cal_cldfra1(nz, sd(cqv.data()), sd(cqc.data()), SD{qi.empty() ? nullptr : cqi.data(), 1},
                    SD{qs.empty() ? nullptr : cqs.data(), 1}, sd(ct.data()), sd(cp.data()), sd(ccf.data()));
      if (clouds && o_.icloud_cu == 2) {  // KF (:1333-1347): cldfra_cu = cldfra_dp + cldfra_sh over CLDFRA
        for (int k = 0; k < nz; ++k) {
          const std::size_t q = Lm.index(i, j, k);
          const rr cu = rr(real(cfdp[q] + cfsh[q]));
          ccf[k] = (rr(1.) - cu) * ccf[k];
          ccf[k] = ccf[k] + cu;
          ccf[k] = std::min(rr(1.0), ccf[k]);
          if (!qc.empty()) cqc[k] = cqc[k] + qccu[q] * cu;
          if (!qi.empty()) cqi[k] = cqi[k] + qicu[q] * cu;
        }
      }
      if (clouds && bl_clouds) {  // BL clouds (icloud_bl > 0)
        for (int k = 0; k < nz; ++k) {
          const std::size_t q = Lm.index(i, j, k);
          if (itimestep != 1) ccf[k] = cfbl[q];
          if (cqc[k] < rr(1.E-6) && cfbl[q] > rr(0.001)) cqc[k] = cqc[k] + qcbl[q];
          if (cqi[k] < rr(1.E-8) && cfbl[q] > rr(0.001)) cqi[k] = cqi[k] + qibl[q];
        }
      }
      if (o_.bmj_feedback) {  // BMJ (module_radiation_driver.F:1578-1591, PCC scheme)
        if (!clouds) std::fill(ccf.begin(), ccf.end(), rr(0.0));  // CLDFRA = 0 (the column buffer is reused)
        for (int k = 0; k < nz; ++k) {
          const std::size_t q = Lm.index(i, j, k);
          cqc[k] = cqc[k] + rr(qcconv[q]);
          cqi[k] = cqi[k] + rr(qiconv[q]);
          const rr cc = rr(ccldfra[q]);
          const rr itrmx = c_min(ccf[k], cc);
          const rr itrmn = c_max(rr(0.), ccf[k] + cc - rr(1.));
          ccf[k] = ccf[k] + cc - rr(0.5) * (itrmx + itrmn);
        }
      }
      if (out.cldfra)
        for (int k = 0; k < nz; ++k) cldfra[out.cldfra->layout().index(i, j, k)] = real(ccf[k]);
      const rr* prc = rec.empty() ? nullptr : crec.data();
      const rr* pri = rei.empty() ? nullptr : crei.data();
      const rr* prs = res.empty() ? nullptr : cres.data();
      const std::size_t q2 = L2.index(i, j, 0), o2i = Lo2.index(i, j, 0);
      rr hr = 0;
      const rr cz = calc_coszen(julian, xtime + o_.radt * rr(0.5), o_.gmt, declin, p_->xlong[c], p_->xlat[c], i, nx, &hr);
      if (diag) hrang2[o2i] = real(hr);
      coszen[o2i] = real(cz);
      if (o_.o3input == 2 && o3nest)
        for (int k = 0; k < nz; ++k) co3[k] = rr(o3h[Lo3.index(i, j, k)]);
      else if (o_.o3input == 2)
        ozone_column(&p_->ozmixm[std::size_t(oc) * per + nm * CamOzone::levsiz],
                     &p_->ozmixm[std::size_t(oc) * per + np * CamOzone::levsiz], f1, f2,
                     p_->ozone.pin.data(), nz, sd(cp.data()), sd(co3.data()));
      rr g_lw, olr, lwcf = 0, a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0, a6 = 0, a7 = 0, a8 = 0;
      if (rrtm) {  // RRTMLWRAD: no LW flux diagnostics, LWCF untouched (0)
#ifdef CPPWRF_DOUBLE
        const RrtmLwColumn rc{cp8w.data(), ct8w.data(), cp.data(), ct.data(), cdz.data(), cpi.data(), cqv.data(),
                              in.qc ? cqc.data() : nullptr, in.qr ? cqr.data() : nullptr,
                              in.qi ? cqi.data() : nullptr, in.qs ? cqs.data() : nullptr, ccf.data(),
                              in.qi != nullptr, rr(tsk[q2]), rr(ems[q2]), lc.co2, lc.n2o, lc.ch4, cst::r_d,
                              cst::g, o_.icloud, o_.warm_rain};
        rrtm_lw_column(*p_->rrtm, rc, *rws, g_lw, olr, rlw.data(), rlwc.data(), rmsg);
#else
        // FP32 build: RRTM LW computes in double (rrtm_lw.cpp is not converted); the column is widened and the
        // results rounded once
        auto dv = [](const std::vector<rr>& v) { return std::vector<double>(v.begin(), v.end()); };
        const std::vector<double> dp8w = dv(cp8w), dt8w = dv(ct8w), dp = dv(cp), dt = dv(ct), ddz = dv(cdz),
                                  dpi = dv(cpi), dqv = dv(cqv), dqc = dv(cqc), dqr = dv(cqr), dqi = dv(cqi),
                                  dqs = dv(cqs), dcf = dv(ccf);
        const RrtmLwColumn rc{dp8w.data(), dt8w.data(), dp.data(), dt.data(), ddz.data(), dpi.data(), dqv.data(),
                              in.qc ? dqc.data() : nullptr, in.qr ? dqr.data() : nullptr,
                              in.qi ? dqi.data() : nullptr, in.qs ? dqs.data() : nullptr, dcf.data(),
                              in.qi != nullptr, tsk[q2], ems[q2], lc.co2, lc.n2o, lc.ch4, cst::r_d,
                              cst::g, o_.icloud, o_.warm_rain};
        std::vector<double> drlw(nz), drlwc(nz);
        double dglw = 0, dolr = 0;
        rrtm_lw_column(*p_->rrtm, rc, *rws, dglw, dolr, drlw.data(), drlwc.data(), rmsg);
        g_lw = rr(dglw);
        olr = rr(dolr);
        for (int k = 0; k < nz; ++k) {
          rlw[k] = rr(drlw[k]);
          rlwc[k] = rr(drlwc[k]);
        }
#endif
        if (!rmsg.empty()) {
#pragma omp critical(rrtm_lw_message)
          std::fprintf(stderr, "%s\n", rmsg.c_str());
        }
      } else {
        split_mp3();
        rrtmg_lw_column(nz, lc, *p_->lw, tsk[q2], ems[q2], xland[q2], xice[q2], snow[q2], p_->xlat[c], cp8w.data(),
                        ct8w.data(), cp.data(), cpi.data(), cdz.data(), ct.data(), crho.data(), cqv.data(), sqc.data(),
                        zero.data(), sqi.data(), sqs.data(), zero.data(), ccf.data(), co3.data(), g_lw, olr, lwcf, a1,
                        a2, a3, a4, a5, a6, a7, a8, rlw.data(), rlwc.data(), *lws, prc, pri, prs,
                        o_.o3input != 2 && diag ? co3.data() : nullptr);
      }
      glw[o2i] = real(g_lw);
      if (diag) {
        const rr d[9] = {a1, a2, a3, a4, a5, a6, a7, a8, olr};
        for (int q = 0; q < 9; ++q) lwd[std::size_t(q)][o2i] = real(d[q]);
        lwcf2[o2i] = real(lwcf);
        // O3RAD: the CAM ozone (o3input = 2), else RRTMG_LWRAD's climatology ("output o3 for o3input=0"); RRTM
        // does not write it
        if (!rrtm || o_.o3input == 2)
          for (int k = 0; k < nz; ++k) o3r[Lo3.index(i, j, k)] = real(co3[k]);
      }
      rr g_sw = 0;
      std::fill(rsw.begin(), rsw.end(), rr(0.0));
      if (dudhia) {
        // SWRAD on every column (0 at night): RTHRATEN = RTHRATEN + TTEN1D(NK)/pi3D; the qc / qi the radiation
        // sees (BL clouds added); RTHRATENSWC stays 0
        auto opt = [](bool has, const std::vector<rr>& a) { return has ? CSD(a.data(), 1) : CSD(); };
        g_sw = dudhia_swpara(nz, CSD(ct.data(), 1), CSD(cp.data(), 1), CSD(cdz.data(), 1), CSD(cqv.data(), 1),
                             opt(in.qc, cqc), opt(in.qr, cqr), opt(in.qi, cqi), opt(in.qs, cqs), opt(in.qg, cqg),
                             rr(coszen[o2i]), rr(alb[q2]), solcon, o_.icloud, dudhia_cssca(o_.swrad_scat),
                             o_.warm_rain, sd(rsw.data()));
        for (int k = 0; k < nz; ++k) rsw[k] = rsw[k] / cpi[k];
        std::fill(rswc.begin(), rswc.end(), rr(0.0));
      } else if (cz > rr(0.0)) {
        SwColumnOut so{};
        split_mp3();
        rrtmg_sw_column(nz, sc, *p_->sw, cz, alb[q2], xland[q2], xice[q2], snow[q2], p_->xlat[c], cp8w.data(),
                        ct8w.data(), cp.data(), cpi.data(), cdz.data(), ct.data(), crho.data(), cqv.data(), sqc.data(),
                        zero.data(), sqi.data(), sqs.data(), zero.data(), ccf.data(), co3.data(), so, rsw.data(),
                        rswc.data(), *sws, prc, pri, prs);
        g_sw = so.gsw;
        if (diag) {
          const rr d[8] = {so.swupt, so.swuptc, so.swdnt, so.swdntc, so.swupb, so.swupbc, so.swdnb, so.swdnbc};
          for (int q = 0; q < 8; ++q) swd[std::size_t(q)][o2i] = real(d[q]);
          swcf2[o2i] = real(so.swcf);
        }
        if (!swsfc.empty()) {  // RRTMG_SWRAD: SWDDNI = SWDDIR / coszen, SWDDNIC = SWDDIRC / coszen
          const real dir = real(so.swddir), dirc = real(so.swddirc);
          const real v[9] = {real(so.swvisdir), real(so.swvisdif), real(so.swnirdir), real(so.swnirdif), dir,
                             real(rr(dir) / cz), real(so.swddif), dirc, real(rr(dirc) / cz)};
          for (int q = 0; q < 9; ++q) swsfc[std::size_t(q)][o2i] = v[q];
        }
      }
      gsw[o2i] = real(g_sw);
      // SWDOWN = GSW/(1-ALBEDO) (all columns)
      swdown[o2i] = real(rr(gsw[o2i]) / (rr(1.) - rr(alb[q2])));
      // Dudhia: radiation_driver's direct / diffuse split (swddir, swddni, swddif: 0 where coszen <= 1e-3)
      if (dudhia && !swsfc.empty()) {
        rr dif = 0, dir = 0;
        const real czr = coszen[o2i];
        if (sw_diffuse_backup(rr(czr), rr(swdown[o2i]), solcon, rr(ht[in.ht->layout().index(i, j, 0)]),
                              dif, dir)) {
          swsfc[6][o2i] = real(dif);
          swsfc[4][o2i] = real(dir);
          swsfc[5][o2i] = real(rr(swsfc[4][o2i]) / std::max(rr(czr), rr(1e-4)));
        }
      }
      if (!swsfc.empty()) {  // DIFFUSE_FRAC (radiation_driver, after the schemes)
        const real dif = swsfc[6][o2i], sd = swdown[o2i];
        swsfc[9][o2i] = rr(sd) > rr(0.001) ? real(std::min(rr(dif) / rr(sd), rr(1.0))) : real(0);
      }
      for (int k = 0; k < nz; ++k) {
        rr r = rlw[k];  // RTHRATEN = LW, then RTHRATEN = RTHRATEN + RTHRATENSW (daytime; Dudhia: always)
        if (dudhia || cz > rr(0.0)) r = r + rsw[k];
        rth[Lo3.index(i, j, k)] = real(r);
        if (diag) {
          rthlw[Lo3.index(i, j, k)] = real(rlw[k]);
          rthsw[Lo3.index(i, j, k)] = real(r) - real(rlw[k]);  // RTHRATENSW = RTHRATEN - RTHRATENLW
          rthlwc[Lo3.index(i, j, k)] = real(rlwc[k]);
          rthswc[Lo3.index(i, j, k)] = real(!dudhia && cz > rr(0.0) ? rswc[k] : rr(0.0));
        }
      }
    }
  }
  copy_to_dev(out.rthraten->data() - Lo3.offset, rth.data(), rth.size() * sizeof(real));
  copy_to_dev(out.glw->data() - Lo2.offset, glw.data(), glw.size() * sizeof(real));
  copy_to_dev(out.gsw->data() - out.gsw->layout().offset, gsw.data(), gsw.size() * sizeof(real));
  copy_to_dev(out.swdown->data() - out.swdown->layout().offset, swdown.data(), swdown.size() * sizeof(real));
  copy_to_dev(out.coszen->data() - out.coszen->layout().offset, coszen.data(), coszen.size() * sizeof(real));
  if (out.cldfra)
    copy_to_dev(out.cldfra->data() - out.cldfra->layout().offset, cldfra.data(), cldfra.size() * sizeof(real));
  if (diag) {
    copy_to_dev(out.rthratlw->data() - Lo3.offset, rthlw.data(), rthlw.size() * sizeof(real));
    copy_to_dev(out.rthratsw->data() - Lo3.offset, rthsw.data(), rthsw.size() * sizeof(real));
    for (int q = 0; q < 9; ++q) copy_to_dev(out.lw[q]->data() - Lo2.offset, lwd[std::size_t(q)].data(), Lo2.total * sizeof(real));
    if (!dudhia)
      for (int q = 0; q < 8; ++q) copy_to_dev(out.sw[q]->data() - Lo2.offset, swd[std::size_t(q)].data(), Lo2.total * sizeof(real));
    copy_to_dev(out.rthratlwc->data() - Lo3.offset, rthlwc.data(), rthlwc.size() * sizeof(real));
    copy_to_dev(out.rthratswc->data() - Lo3.offset, rthswc.data(), rthswc.size() * sizeof(real));
    if (o_.o3input != 2 || !o3nest) copy_to_dev(out.o3rad->data() - Lo3.offset, o3r.data(), o3r.size() * sizeof(real));
    copy_to_dev(out.lwcf->data() - Lo2.offset, lwcf2.data(), Lo2.total * sizeof(real));
    if (!dudhia) copy_to_dev(out.swcf->data() - Lo2.offset, swcf2.data(), Lo2.total * sizeof(real));
    copy_to_dev(out.hrang->data() - Lo2.offset, hrang2.data(), Lo2.total * sizeof(real));
  }
  for (std::size_t q = 0; q < swsfc.size(); ++q)
    if (!dudhia || q == 4 || q == 5 || q == 6 || q == 9)
      copy_to_dev(out.swsfc[q]->data() - Lo2.offset, swsfc[q].data(), Lo2.total * sizeof(real));
  return true;
#endif
}

}  // namespace cppwrf::phys::driver
