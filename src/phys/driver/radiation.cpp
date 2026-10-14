#include "phys/driver/radiation.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>

#include "phys/ra_rrtmg_common/cam_gases.hpp"
#include "phys/ra_rrtmg_common/cam_ozone.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_tables.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_tables.hpp"
#if defined(CPPWRF_USE_CUDA)
#include "phys/ra_rrtmg_lw/rrtmg_lw_batch.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_batch.hpp"
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

namespace cppwrf::phys::driver {

namespace {
constexpr double piconst = 3.1415926535897932384626433;
constexpr double degrad = piconst / 180.;
constexpr double dpd = 360. / 365.;

// radconst: declination (rad) and solar constant from the julian day
void radconst(double julian, double& declin, double& solcon) {
  const double obecl = 23.5 * degrad, sinob = std::sin(obecl);
  double sxlong = julian >= 80. ? dpd * (julian - 80.) : dpd * (julian + 285.);
  sxlong = sxlong * degrad;
  const double arg = sinob * std::sin(sxlong);
  declin = std::asin(arg);
  const double djul = julian * 360. / 365., rjul = djul * degrad;
  const double eccfac = 1.000110 + 0.034221 * std::cos(rjul) + 0.001280 * std::sin(rjul) +
                        0.000719 * std::cos(2 * rjul) + 0.000077 * std::sin(2 * rjul);
  solcon = 1370. * eccfac;
}


#if defined(CPPWRF_USE_CUDA)
// device path: model fields -> double arrays with the column index fastest (col = i + nx*j)
struct KStage3 {
  View3<const real> f;
  double* d;
  int nx, ncol;
  CPPWRF_HD void operator()(int i, int j, int k) const { d[std::size_t(k) * ncol + i + nx * j] = double(f(i, j, k)); }
};
struct KStage2 {
  View2<const real> f;
  double* d;
  int nx;
  CPPWRF_HD void operator()(int i, int j, int) const { d[i + nx * j] = double(f(i, j)); }
};
// per column: calc_coszen at xtime + radt/2, CAM ozone (ozn_time_int, ozn_p_int), cal_cldfra1
struct KRadColumn {
  int nx, nz, ncol;
  double julian, xtime_mid, gmt, declin;
  const double *xlat, *xlong;
  int o3input;
  const double *ozmixm, *pin;
  int nm, np;
  double f1, f2;
  bool clouds;
  const double *p, *t, *qv, *qc, *qi, *qs;
  double *o3, *cf, *cz;
  bool has_cldfra;
  View3<real> cldfra_out;
  View2<real> coszen_out;
  CPPWRF_HD void operator()(int i, int j) const {
    const int col = i + nx * j;
    const double c = calc_coszen(julian, xtime_mid, gmt, declin, xlong[col], xlat[col]);
    cz[col] = c;
    coszen_out(i, j) = real(c);
    constexpr int per = CamOzone::nmon * CamOzone::levsiz;
    auto v = [&](const double* a) { return SD{const_cast<double*>(a) + col, ncol}; };
    if (o3input == 2)
      ozone_column(&ozmixm[std::size_t(col) * per + nm * CamOzone::levsiz],
                   &ozmixm[std::size_t(col) * per + np * CamOzone::levsiz], f1, f2, pin, nz, v(p), v(o3));
    if (clouds) {
      cal_cldfra1(nz, v(qv), v(qc), qi ? v(qi) : SD{nullptr, 0}, qs ? v(qs) : SD{nullptr, 0}, v(t), v(p), v(cf));
    } else {
      for (int k = 0; k < nz; ++k) cf[std::size_t(k) * ncol + col] = 0.0;
    }
    if (has_cldfra)
      for (int k = 0; k < nz; ++k) cldfra_out(i, j, k) = real(cf[std::size_t(k) * ncol + col]);
  }
};
// RTHRATEN = LW (+ SW in daytime), GLW, GSW, SWDOWN = GSW/(1 - ALBEDO)
struct KRadOut {
  int nx, ncol;
  const double *rlw, *rsw, *glw, *gsw, *cz;
  View3<real> rth;
  View2<real> glw_out, gsw_out, swdown_out;
  View2<const real> alb;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    const int col = i + nx * j;
    const std::size_t q = std::size_t(k) * ncol + col;
    double r = rlw[q];
    if (cz[col] > 0.0) r = r + rsw[q];
    rth(i, j, k) = real(r);
    if (k == 0) {
      glw_out(i, j) = real(glw[col]);
      const real gs = real(cz[col] > 0.0 ? gsw[col] : 0.0);
      gsw_out(i, j) = gs;
      swdown_out(i, j) = real(double(gs) / (1. - double(alb(i, j))));
    }
  }
};
#endif

}  // namespace

struct RadiationDriver::Impl {
  std::unique_ptr<RrtmgLwTables> lw;
  std::unique_ptr<RrtmgSwTables> sw;
  CamGases gases;
  CamOzone ozone;
  std::vector<double> xlat, xlong, ozmixm;  // ozmixm: [col][month][level] (oznini)
#if defined(CPPWRF_USE_CUDA)
  std::unique_ptr<RrtmgLwBatch> lwb;
  std::unique_ptr<RrtmgSwBatch> swb;
  std::vector<void*> bufs;
  double* dev(std::size_t n) {
    void* q = dev_alloc(n * sizeof(double));
    bufs.push_back(q);
    return static_cast<double*>(q);
  }
  double *xlat_d = nullptr, *xlong_d = nullptr, *oz_d = nullptr, *pin_d = nullptr;
  // staged inputs, outputs (column fastest)
  double *p8w, *t8w, *p3d, *pi3d, *t3d, *qv, *qc, *qi, *qs, *rec, *rei, *res, *cf, *o3;
  double *tsk, *emiss, *xland, *xice, *snow, *albedo, *cz;
  double *rlw, *rsw, *glw, *gsw;
  ~Impl() {
    for (void* q : bufs) dev_free(q);
  }
#endif
};

RadiationDriver::RadiationDriver(const Grid& g, const RadiationOptions& o, const std::vector<double>& xlat,
                                 const std::vector<double>& xlong)
    : p_(std::make_unique<Impl>()), g_(g), o_(o) {
  if (o.ra_lw_physics != 4 || o.ra_sw_physics != 4)
    throw std::runtime_error("radiation: only ra_lw_physics = ra_sw_physics = 4 (RRTMG) are ported");
  if (o.aer_opt != 0) throw std::runtime_error("radiation: aer_opt != 0 is not ported");
  if (o.o3input != 2 && o.o3input != 0) throw std::runtime_error("radiation: o3input must be 0 or 2");
  stepra_ = o.radt == 0 ? 1 : std::max(1, int(std::lround(o.radt * 60. / g.dt)));
  p_->lw = std::make_unique<RrtmgLwTables>();
  p_->sw = std::make_unique<RrtmgSwTables>();
  // RRTMG tables: phys.rrtmg_tables, then the copy in the source tree (data/rrtmg, CPPWRF_DATA_DIR)
  auto table_file = [&](const std::string& name) {
    const std::string a = o.table_dir + "/" + name;
    if (std::ifstream(a).good()) return a;
#ifdef CPPWRF_DATA_DIR
    const std::string b = std::string(CPPWRF_DATA_DIR) + "/rrtmg/" + name;
    if (std::ifstream(b).good()) return b;
#endif
    throw std::runtime_error("radiation: cannot find " + name + " in " + o.table_dir + " or the data directory");
  };
  if (!load_rrtmg_lw_tables(table_file("rrtmg_lw_tables.bin"), *p_->lw))
    throw std::runtime_error("radiation: cannot read rrtmg_lw_tables.bin");
  if (!load_rrtmg_sw_tables(table_file("rrtmg_sw_tables.bin"), *p_->sw))
    throw std::runtime_error("radiation: cannot read rrtmg_sw_tables.bin");
  if (o.ghg_input == 1) p_->gases.load(o.table_dir + "/CAMtr_volume_mixing_ratio");
  p_->xlat = xlat;
  p_->xlong = xlong;
  if (o.o3input == 2) {  // oznini (phy_init, domain 1)
    p_->ozone.load(o.table_dir);
    const std::size_t nc = std::size_t(g.nx) * g.ny, per = std::size_t(CamOzone::nmon) * CamOzone::levsiz;
    p_->ozmixm.resize(nc * per);
    for (std::size_t c = 0; c < nc; ++c) p_->ozone.column(xlat[c], &p_->ozmixm[c * per]);
  }
#if defined(CPPWRF_USE_CUDA)
  {
    Impl& q = *p_;
    const std::size_t nc = std::size_t(g.nx) * g.ny, n3 = nc * g.nz, n3w = nc * (g.nz + 1);
    auto up = [&](const std::vector<double>& h) {
      double* d = q.dev(std::max<std::size_t>(h.size(), 1));
      if (!h.empty()) copy_to_dev(d, h.data(), h.size() * sizeof(double));
      return d;
    };
    q.xlat_d = up(xlat);
    q.xlong_d = up(xlong);
    q.oz_d = up(q.ozmixm);
    q.pin_d = up(q.ozone.pin);
    q.p8w = q.dev(n3w); q.t8w = q.dev(n3w);
    for (double** a : {&q.p3d, &q.pi3d, &q.t3d, &q.qv, &q.qc, &q.qi, &q.qs, &q.rec, &q.rei, &q.res, &q.cf, &q.o3,
                       &q.rlw, &q.rsw})
      *a = q.dev(n3);
    for (double** a : {&q.tsk, &q.emiss, &q.xland, &q.xice, &q.snow, &q.albedo, &q.cz, &q.glw, &q.gsw}) *a = q.dev(nc);
    RrtmgLwConfig lc;
    lc.p_top = o.p_top;
    q.lwb = std::make_unique<RrtmgLwBatch>(*q.lw, rrtmg_lw_call_config(g.nz, lc, false, false, false).nlayers, 2048);
    q.swb = std::make_unique<RrtmgSwBatch>(*q.sw, g.nz + 1, 2048);
  }
#endif
}

RadiationDriver::~RadiationDriver() = default;

bool RadiationDriver::run(long itimestep, double xtime, double julian, const RadiationInputs& in,
                          const RadiationOutputs& out) {
  // run decision (radiation_driver, no adaptive time step)
  const bool run_param = itimestep == 1 || o_.radt == 0. || stepra_ == 1 ||
                         int(itimestep % stepra_) == 1 + o_.ra_call_offset;
  if (!run_param) return false;
  double declin, solcon;
  radconst(julian, declin, solcon);  // Solar_step (radconst before calc_coszen)
  RrtmgLwConfig lc;
  lc.icloud = o_.icloud; lc.cldovrlp = o_.cldovrlp; lc.idcor = o_.idcor; lc.p_top = o_.p_top;
  lc.yr = o_.julyr; lc.julian = julian; lc.o3input = o_.o3input;
  if (o_.ghg_input == 1) p_->gases.interpolate(o_.julyr, julian, lc.co2, lc.n2o, lc.ch4, lc.cfc11, lc.cfc12);
  else rrtmg_default_gases(lc);
  RrtmgSwConfig sc;
  sc.icloud = o_.icloud; sc.cldovrlp = o_.cldovrlp; sc.idcor = o_.idcor; sc.p_top = o_.p_top; sc.yr = o_.julyr;
  sc.julian = julian; sc.solcon = solcon; sc.o3input = o_.o3input; sc.co2 = lc.co2; sc.n2o = lc.n2o; sc.ch4 = lc.ch4;
  int nm = 0, np = 0;
  double f1 = 0, f2 = 0;
  if (o_.o3input == 2) ozone_time_weights(julian, nm, np, f1, f2);
  const std::size_t per = std::size_t(CamOzone::nmon) * CamOzone::levsiz;
#if defined(CPPWRF_USE_CUDA)
  {
    (void)per;
    Impl& q = *p_;
    const int nx = g_.nx, ny = g_.ny, nz = g_.nz, ncol = nx * ny;
    if (o_.icloud != 0 && o_.icloud != 1) throw std::runtime_error("radiation: only icloud = 0 or 1 are ported");
    if (in.qc && !in.qi && in.qs) throw std::runtime_error("radiation: F_QC, F_QS without F_QI is not ported");
    const bool clouds = o_.icloud != 0 && in.qc != nullptr;
    auto st3 = [&](const Field3* f, double* d, int nlev) {
      parallel_for(Range3D{0, nx, 0, ny, 0, nlev}, KStage3{f->cview(), d, nx, ncol});
    };
    auto st2 = [&](const Field2* f, double* d) { parallel_for(Range3D{0, nx, 0, ny, 0, 1}, KStage2{f->cview(), d, nx}); };
    st3(in.p_hyd_w, q.p8w, nz + 1); st3(in.t8w, q.t8w, nz + 1); st3(in.p_hyd, q.p3d, nz); st3(in.pi_phy, q.pi3d, nz);
    st3(in.t_phy, q.t3d, nz); st3(in.qv, q.qv, nz);
    if (in.qc) st3(in.qc, q.qc, nz);
    if (in.qi) st3(in.qi, q.qi, nz);
    if (in.qs) st3(in.qs, q.qs, nz);
    if (in.re_cloud) st3(in.re_cloud, q.rec, nz);
    if (in.re_ice) st3(in.re_ice, q.rei, nz);
    if (in.re_snow) st3(in.re_snow, q.res, nz);
    st2(in.tsk, q.tsk); st2(in.emiss, q.emiss); st2(in.xland, q.xland); st2(in.xice, q.xice); st2(in.snow, q.snow);
    st2(in.albedo, q.albedo);
    KRadColumn kc{};
    kc.nx = nx; kc.nz = nz; kc.ncol = ncol;
    kc.julian = julian; kc.xtime_mid = xtime + o_.radt * 0.5; kc.gmt = o_.gmt; kc.declin = declin;
    kc.xlat = q.xlat_d; kc.xlong = q.xlong_d;
    kc.o3input = o_.o3input; kc.ozmixm = q.oz_d; kc.pin = q.pin_d; kc.nm = nm; kc.np = np; kc.f1 = f1; kc.f2 = f2;
    kc.clouds = clouds;
    kc.p = q.p3d; kc.t = q.t3d; kc.qv = q.qv; kc.qc = in.qc ? q.qc : nullptr; kc.qi = in.qi ? q.qi : nullptr;
    kc.qs = in.qs ? q.qs : nullptr;
    kc.o3 = q.o3; kc.cf = q.cf; kc.cz = q.cz;
    kc.has_cldfra = out.cldfra != nullptr;
    if (out.cldfra) kc.cldfra_out = out.cldfra->view();
    kc.coszen_out = out.coszen->view();
    column_for(Range2D{0, nx, 0, ny}, kc);
    LwBatchInput li{};
    li.ncol = ncol; li.tsk = q.tsk; li.emiss = q.emiss; li.xland = q.xland; li.xice = q.xice; li.snow = q.snow;
    li.p8w = q.p8w; li.t8w = q.t8w; li.p3d = q.p3d; li.pi3d = q.pi3d; li.t3d = q.t3d; li.qv = q.qv; li.cldfra = q.cf;
    li.o3 = q.o3;
    li.qc = in.qc ? q.qc : nullptr; li.qi = in.qi ? q.qi : nullptr; li.qs = in.qs ? q.qs : nullptr;
    li.re_cloud = in.re_cloud ? q.rec : nullptr; li.re_ice = in.re_ice ? q.rei : nullptr;
    li.re_snow = in.re_snow ? q.res : nullptr;
    q.lwb->run(rrtmg_lw_call_config(nz, lc, in.re_cloud, in.re_ice, in.re_snow), li, LwBatchOutput{q.rlw, q.glw, nullptr});
    SwBatchInput si{};
    si.ncol = ncol; si.coszen = q.cz; si.albedo = q.albedo; si.xland = q.xland; si.xice = q.xice; si.snow = q.snow;
    si.p8w = q.p8w; si.t8w = q.t8w; si.p3d = q.p3d; si.pi3d = q.pi3d; si.t3d = q.t3d; si.qv = q.qv; si.cldfra = q.cf;
    si.o3 = q.o3; si.qc = li.qc; si.qi = li.qi; si.qs = li.qs;
    si.re_cloud = li.re_cloud; si.re_ice = li.re_ice; si.re_snow = li.re_snow;
    q.swb->run(rrtmg_sw_call_config(nz, sc), si, SwBatchOutput{q.rsw, q.gsw});
    KRadOut ko{nx, ncol, q.rlw, q.rsw, q.glw, q.gsw, q.cz, out.rthraten->view(), out.glw->view(), out.gsw->view(),
               out.swdown->view(), in.albedo->cview()};
    parallel_for(Range3D{0, nx, 0, ny, 0, nz}, ko);
    return true;
  }
#else
  const int nx = g_.nx, ny = g_.ny, nz = g_.nz;
  // inputs to the host
  std::vector<real> t, p, pw, pi, rho, dz, tw, qv, qc, qi, qs, rec, rei, res, tsk, alb, ems, xland, xice, snow;
  in.t_phy->download(t); in.p_hyd->download(p); in.p_hyd_w->download(pw); in.pi_phy->download(pi);
  in.rho->download(rho); in.dz8w->download(dz); in.t8w->download(tw); in.qv->download(qv);
  if (in.qc) in.qc->download(qc);
  if (in.qi) in.qi->download(qi);
  if (in.qs) in.qs->download(qs);
  if (in.re_cloud) in.re_cloud->download(rec);
  if (in.re_ice) in.re_ice->download(rei);
  if (in.re_snow) in.re_snow->download(res);
  // icloud = 1 with F_QC or F_QI: cal_cldfra1 (else CLDFRA = 0); qc is required with qi (WSM/Kessler)
  const bool clouds = o_.icloud != 0 && in.qc != nullptr;
  if (o_.icloud != 0 && o_.icloud != 1) throw std::runtime_error("radiation: only icloud = 0 or 1 are ported");
  if (in.qc && !in.qi && in.qs) throw std::runtime_error("radiation: F_QC, F_QS without F_QI is not ported");
  in.tsk->download(tsk); in.albedo->download(alb); in.emiss->download(ems); in.xland->download(xland);
  in.xice->download(xice); in.snow->download(snow);
  const Layout &Lm = in.t_phy->layout(), &Lw = in.p_hyd_w->layout(), &Lt8 = in.t8w->layout(),
               &L2 = in.tsk->layout();
  std::vector<real> rth(out.rthraten->layout().total, real(0)), glw(out.glw->layout().total, real(0)),
      gsw(out.gsw->layout().total, real(0)), swdown(out.swdown->layout().total, real(0)),
      coszen(out.coszen->layout().total, real(0)), cldfra(out.cldfra ? out.cldfra->layout().total : 0, real(0));
  const Layout &Lo3 = out.rthraten->layout(), &Lo2 = out.glw->layout();


#pragma omp parallel
  {
    auto lws = std::make_unique<RrtmgLwWorkspace>();
    lws->init_pointers();
    auto sws = std::make_unique<RrtmgSwWorkspace>();
    std::vector<double> cp8w(nz + 1), ct8w(nz + 1), cp(nz), cpi(nz), cdz(nz), ct(nz), crho(nz), cqv(nz),
        zero(nz, 0.0), ccf(nz, 0.0), co3(nz, 0.0), rlw(nz), rlwc(nz), rsw(nz), rswc(nz), cqc(nz, 0.0),
        cqi(nz, 0.0), cqs(nz, 0.0), crec(nz), crei(nz), cres(nz);
#pragma omp for schedule(dynamic, 16)
    for (int c = 0; c < nx * ny; ++c) {
      const int i = c % nx, j = c / nx;
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
      }
      if (clouds)
        cal_cldfra1(nz, sd(cqv.data()), sd(cqc.data()), SD{qi.empty() ? nullptr : cqi.data(), 1},
                    SD{qs.empty() ? nullptr : cqs.data(), 1}, sd(ct.data()), sd(cp.data()), sd(ccf.data()));
      if (out.cldfra)
        for (int k = 0; k < nz; ++k) cldfra[out.cldfra->layout().index(i, j, k)] = real(ccf[k]);
      const double* prc = rec.empty() ? nullptr : crec.data();
      const double* pri = rei.empty() ? nullptr : crei.data();
      const double* prs = res.empty() ? nullptr : cres.data();
      const std::size_t q2 = L2.index(i, j, 0), o2i = Lo2.index(i, j, 0);
      const double cz = calc_coszen(julian, xtime + o_.radt * 0.5, o_.gmt, declin, p_->xlong[c], p_->xlat[c]);
      coszen[o2i] = real(cz);
      if (o_.o3input == 2)
        ozone_column(&p_->ozmixm[c * per + nm * CamOzone::levsiz], &p_->ozmixm[c * per + np * CamOzone::levsiz], f1, f2,
                     p_->ozone.pin.data(), nz, sd(cp.data()), sd(co3.data()));
      double g_lw, olr, lwcf, a1, a2, a3, a4, a5, a6, a7, a8;
      rrtmg_lw_column(nz, lc, *p_->lw, tsk[q2], ems[q2], xland[q2], xice[q2], snow[q2], p_->xlat[c], cp8w.data(),
                      ct8w.data(), cp.data(), cpi.data(), cdz.data(), ct.data(), crho.data(), cqv.data(), cqc.data(),
                      zero.data(), cqi.data(), cqs.data(), zero.data(), ccf.data(), co3.data(), g_lw, olr, lwcf, a1,
                      a2, a3, a4, a5, a6, a7, a8, rlw.data(), rlwc.data(), *lws, prc, pri, prs);
      glw[o2i] = real(g_lw);
      double g_sw = 0;
      std::fill(rsw.begin(), rsw.end(), 0.0);
      if (cz > 0.0) {
        double s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13;
        rrtmg_sw_column(nz, sc, *p_->sw, cz, alb[q2], xland[q2], xice[q2], snow[q2], p_->xlat[c], cp8w.data(),
                        ct8w.data(), cp.data(), cpi.data(), cdz.data(), ct.data(), crho.data(), cqv.data(), cqc.data(),
                        zero.data(), cqi.data(), cqs.data(), zero.data(), ccf.data(), co3.data(), g_sw, s1, s2, s3, s4,
                        s5, s6, s7, s8, s9, s10, s11, s12, s13, rsw.data(), rswc.data(), *sws, prc, pri, prs);
      }
      gsw[o2i] = real(g_sw);
      // SWDOWN = GSW/(1-ALBEDO) (all columns)
      swdown[o2i] = real(double(gsw[o2i]) / (1. - double(alb[q2])));
      for (int k = 0; k < nz; ++k) {
        double r = rlw[k];  // RTHRATEN = LW, then RTHRATEN = RTHRATEN + RTHRATENSW (daytime)
        if (cz > 0.0) r = r + rsw[k];
        rth[Lo3.index(i, j, k)] = real(r);
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
  return true;
#endif
}

}  // namespace cppwrf::phys::driver
