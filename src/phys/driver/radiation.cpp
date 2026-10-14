#include "phys/driver/radiation.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "phys/ra_rrtmg_common/cam_gases.hpp"
#include "phys/ra_rrtmg_common/cam_ozone.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_tables.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_tables.hpp"
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

// calc_coszen for one column
double calc_coszen(double julian, double xtime, double gmt, double declin, double xlon, double xlat) {
  const double da = 6.2831853071795862 * (julian - 1) / 365.;
  const double eot = (0.000075 + 0.001868 * std::cos(da) - 0.032077 * std::sin(da) - 0.014615 * std::cos(2 * da) -
                      0.04089 * std::sin(2 * da)) * (229.18);
  const double xt24 = std::fmod(xtime, 1440.) + eot;
  const double tloctm = gmt + xt24 / 60. + xlon / 15.;
  const double hrang = 15. * (tloctm - 12.) * degrad;
  const double xxlat = xlat * degrad;
  double c = std::sin(xxlat) * std::sin(declin) + std::cos(xxlat) * std::cos(declin) * std::cos(hrang);
  return std::min(std::max(c, -1.0), 1.0);
}

}  // namespace

// cal_cldfra1 (Randall 1994 / Hong-Pan 1998) for one column, mp_physics != FER_MP_HIRES; qi/qs may be
// nullptr (F_QI / F_QS false). F_QC and not F_QI but F_QS needs F_ICE_PHY (not ported).
void cal_cldfra1(int nz, const double* qv, const double* qc, const double* qi, const double* qs, const double* t,
                 const double* p, double* cldfra) {
  constexpr double ALPHA0 = 100., GAMMA = 0.49, QCLDMIN = 1.E-12, PEXP = 0.25, RHGRID = 1.0;
  constexpr double SVP1 = 0.61078, SVP2 = 17.2693882, SVPI2 = 21.8745584, SVP3 = 35.86, SVPI3 = 7.66,
                   SVPT0 = 273.15, r_d = 287., r_v = 461.6, ep_2 = r_d / r_v;
  for (int k = 0; k < nz; ++k) {
    const double tc = t[k] - SVPT0;
    const double esw = 1000.0 * SVP1 * std::exp(SVP2 * tc / (t[k] - SVP3));
    const double esi = 1000.0 * SVP1 * std::exp(SVPI2 * tc / (t[k] - SVPI3));
    const double QVSW = ep_2 * esw / (p[k] - esw);
    const double QVSI = ep_2 * esi / (p[k] - esi);
    double QCLD = 0, weight = 0;
    if (qi && qs) {
      QCLD = qi[k] + qc[k] + qs[k];
      weight = QCLD < QCLDMIN ? 0. : (qi[k] + qs[k]) / QCLD;
    } else if (qi) {
      QCLD = qi[k] + qc[k];
      weight = QCLD < QCLDMIN ? 0. : qi[k] / QCLD;
    } else {  // F_QC only
      QCLD = qc[k];
      if (QCLD < QCLDMIN) weight = 0.;
      else weight = t[k] > 273.15 ? 0. : 1.;
    }
    const double QVS_WEIGHT = (1 - weight) * QVSW + weight * QVSI;
    double RHUM = qv[k] / QVS_WEIGHT;
    if (QCLD < QCLDMIN) {
      cldfra[k] = 0.;
    } else if (RHUM >= RHGRID) {
      cldfra[k] = 1.;
    } else {
      const double SUBSAT = std::max(1.E-10, RHGRID * QVS_WEIGHT - qv[k]);
      const double DENOM = std::pow(SUBSAT, GAMMA);
      const double ARG = std::max(-6.9, -ALPHA0 * QCLD / DENOM);
      RHUM = std::max(1.E-10, RHUM);
      cldfra[k] = std::pow(RHUM / RHGRID, PEXP) * (1. - std::exp(ARG));
      if (cldfra[k] < .01) cldfra[k] = 0.;
    }
  }
}

struct RadiationDriver::Impl {
  std::unique_ptr<RrtmgLwTables> lw;
  std::unique_ptr<RrtmgSwTables> sw;
  CamGases gases;
  CamOzone ozone;
  std::vector<double> xlat, xlong, ozmixm;  // ozmixm: [col][month][level] (oznini)
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
  if (!load_rrtmg_lw_tables(o.table_dir + "/rrtmg_lw_tables.bin", *p_->lw))
    throw std::runtime_error("radiation: cannot read " + o.table_dir + "/rrtmg_lw_tables.bin");
  if (!load_rrtmg_sw_tables(o.table_dir + "/rrtmg_sw_tables.bin", *p_->sw))
    throw std::runtime_error("radiation: cannot read " + o.table_dir + "/rrtmg_sw_tables.bin");
  if (o.ghg_input == 1) p_->gases.load(o.table_dir + "/CAMtr_volume_mixing_ratio");
  p_->xlat = xlat;
  p_->xlong = xlong;
  if (o.o3input == 2) {  // oznini (phy_init, domain 1)
    p_->ozone.load(o.table_dir);
    const std::size_t nc = std::size_t(g.nx) * g.ny, per = std::size_t(CamOzone::nmon) * CamOzone::levsiz;
    p_->ozmixm.resize(nc * per);
    for (std::size_t c = 0; c < nc; ++c) p_->ozone.column(xlat[c], &p_->ozmixm[c * per]);
  }
}

RadiationDriver::~RadiationDriver() = default;

bool RadiationDriver::run(long itimestep, double xtime, double julian, const RadiationInputs& in,
                          const RadiationOutputs& out) {
  // run decision (radiation_driver, no adaptive time step)
  const bool run_param = itimestep == 1 || o_.radt == 0. || stepra_ == 1 ||
                         int(itimestep % stepra_) == 1 + o_.ra_call_offset;
  if (!run_param) return false;
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
        cal_cldfra1(nz, cqv.data(), cqc.data(), qi.empty() ? nullptr : cqi.data(), qs.empty() ? nullptr : cqs.data(),
                    ct.data(), cp.data(), ccf.data());
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
                     p_->ozone.pin.data(), nz, cp.data(), co3.data());
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
}

}  // namespace cppwrf::phys::driver
