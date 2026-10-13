#include "phys/driver/suite.hpp"

#include <cmath>
#include <stdexcept>

#include "phys/bl_ysu/ysu.hpp"
#include "phys/driver/noah_ice.hpp"
#include "phys/driver/wrf_constants.hpp"
#include "phys/lsm_noah/noah.hpp"

namespace cppwrf::phys::driver {

// ------------------------------------------------------------------ IntField2
IntField2::IntField2(const Grid& g, std::string name)
    : name_(std::move(name)), L_(Layout::make(g, Stagger::M, 2)),
      buf_(static_cast<int*>(dev_alloc(L_.total * sizeof(int)))) {}
IntField2::~IntField2() { dev_free(buf_); }
IntField2::IntField2(IntField2&& o) noexcept : name_(std::move(o.name_)), L_(o.L_), buf_(o.buf_) { o.buf_ = nullptr; }
IntField2& IntField2::operator=(IntField2&& o) noexcept {
  std::swap(name_, o.name_);
  std::swap(L_, o.L_);
  std::swap(buf_, o.buf_);
  return *this;
}
void IntField2::upload(const std::vector<int>& h) { copy_to_dev(buf_, h.data(), L_.total * sizeof(int)); }
void IntField2::download(std::vector<int>& h) const {
  h.resize(L_.total);
  copy_to_host(h.data(), buf_, L_.total * sizeof(int));
}

namespace {

using V2 = View2<real>;
using C2 = View2<const real>;
using V3 = View3<real>;
using C3 = View3<const real>;
using I2 = View2<int>;
using CI2 = View2<const int>;

constexpr int kMaxSoil = noah::kMaxSoilLevels;

C2 slice0(const Field3& f) { return C2{f.data(), int(f.layout().sj)}; }  // level k = 0 (WRF kts)

// ---------------------------------------------------------------- surface_driver, every call
// RAINBL accumulation (RAINCV + RAINNCV [+ RAINSHV]), zeroed locals QGH, CHS, CPM, CHS2.
struct KSfcAlways {
  V2 rainbl, qgh, chs, cpm, chs2, cqs2;
  C2 raincv, rainncv, rainshv;
  V3 u_phytmp, v_phytmp;
  int nz;
  CPPWRF_HD void operator()(int i, int j) const {
    for (int k = 0; k < nz; ++k) {
      v_phytmp(i, j, k) = real(0.);
      u_phytmp(i, j, k) = real(0.);
    }
    qgh(i, j) = real(0.);
    chs(i, j) = real(0.);
    cpm(i, j) = real(0.);
    chs2(i, j) = real(0.);
    cqs2(i, j) = real(0.);  // a local WRF does not zero; SFCLAYREV sets it before any use
    rainbl(i, j) = rainbl(i, j) + raincv(i, j) + rainncv(i, j);
    rainbl(i, j) = rainbl(i, j) + rainshv(i, j);
    rainbl(i, j) = rainbl(i, j) > real(0.0) ? rainbl(i, j) : real(0.0);
  }
};

// PSFC = p8w(kts) (p8w = grid%p_hyd_w); u_phytmp = u_phy + u_frame, surface currents removed.
struct KSfcPrep {
  V2 psfc;
  C3 p8w, u_phy, v_phy;
  C2 uoce, voce;
  V3 u_phytmp, v_phytmp;
  int nz;
  CPPWRF_HD void operator()(int i, int j) const {
    constexpr real u_frame = real(0.), v_frame = real(0.);
    psfc(i, j) = p8w(i, j, 0);
    for (int k = 0; k < nz; ++k) {
      v_phytmp(i, j, k) = v_phy(i, j, k) + v_frame;
      u_phytmp(i, j, k) = u_phy(i, j, k) + u_frame;
    }
    u_phytmp(i, j, 0) = u_phytmp(i, j, 0) - uoce(i, j);
    v_phytmp(i, j, 0) = v_phytmp(i, j, 0) - voce(i, j);
  }
};

struct KCopy2 {
  V2 dst;
  C2 src;
  CPPWRF_HD void operator()(int i, int j) const { dst(i, j) = src(i, j); }
};

// ---------------------------------------------------------------- Noah driver (lsm), before SFLX
struct LsmConst {
  // module_sf_noahdrv.F lsm parameters
  static constexpr real capa = wrfc::r_d / wrfc::cp;
  static constexpr real tresh = real(.95e0), a2 = real(17.67), a3 = real(273.15), a4 = real(29.65);
  static constexpr real a23m4 = a2 * (a3 - a4);
  static constexpr real row = real(1.e3), eliw = wrfc::xlf, rowliw = row * eliw;
};

struct KLsmPre {
  // inputs
  C3 dz8w, qv3d, p8w3d, t3d;
  C2 qgh, glw, swdown, rainbl, vegfra, shdmin, shdmax, snoalb, tmn, xice, chs, rib, sr;
  CI2 ivgtyp, isltyp;
  // in/out
  V2 tsk, albedo, albbck, z0, emiss, embck, snow, snowh, snowc, canwat, snotime, lai, smstav, smstot, chklowq,
      xland;
  V3 smois, tslb, sh2o, smcrel;
  // SFLX / SFLX_GLACIAL inputs (scratch)
  V2 lwdn, soldn, solnet, sfcprs, prcp, sfctmp, th2, q2k, q2sat, dqsdt2, ffrozp, shdfac, shmin, shmax, albbrd,
      snoalb1, tbot, z0brd, emissi, embrd, cmc, t1, snowhk, sneqv, chk, snotime1, ribb, xlai, soilw, sncovr;
  V3 stc, smc, swc;
  I2 vegtyp, soiltyp, slopetyp, ice, active;
  int nsoil, itimestep, isurban, isice, sf_urban_physics;
  int lcz[11];
  bool myj, frpcpn, rdlai2d;
  real dt, xice_threshold;
  CPPWRF_HD void operator()(int i, int j) const {
    using C = LsmConst;
    ice(i, j) = 0;
    active(i, j) = 0;
    if (itimestep == 1) {
      if ((xland(i, j) - real(1.5)) >= real(0.)) {
        smstav(i, j) = real(1.0);
        smstot(i, j) = real(1.0);
        for (int ns = 0; ns < nsoil; ++ns) {
          smois(i, j, ns) = real(1.0);
          tslb(i, j, ns) = real(273.16);
          smcrel(i, j, ns) = real(1.0);
        }
      } else if (xice(i, j) >= xice_threshold) {
        smstav(i, j) = real(1.0);
        smstot(i, j) = real(1.0);
        for (int ns = 0; ns < nsoil; ++ns) {
          smois(i, j, ns) = real(1.0);
          smcrel(i, j, ns) = real(1.0);
        }
      }
    }
    real psfc = p8w3d(i, j, 0);
    sfcprs(i, j) = (p8w3d(i, j, 1) + p8w3d(i, j, 0)) * real(0.5);
    q2k(i, j) = qv3d(i, j, 0) / (real(1.0) + qv3d(i, j, 0));
    soilw(i, j) = real(0.);
    real q2s = qgh(i, j) / (real(1.0) + qgh(i, j));
    if (myj && (q2k(i, j) >= q2s * C::tresh)) {  // QZ0 test of MYJ (not used with sfclayrev)
      chklowq(i, j) = real(0.);
    } else {
      chklowq(i, j) = real(1.);
    }
    real sfct = t3d(i, j, 0);
    sfctmp(i, j) = sfct;
    real apes = rpow(real(1.e5) / psfc, C::capa);
    real apelm = rpow(real(1.e5) / sfcprs(i, j), C::capa);
    real sfcth2 = sfct * apelm;
    th2(i, j) = sfcth2 / apes;
    real emis = emiss(i, j);
    emissi(i, j) = emis;
    lwdn(i, j) = glw(i, j) * emis;
    real sold = swdown(i, j);
    soldn(i, j) = sold;
    solnet(i, j) = sold * (real(1.) - albedo(i, j));
    prcp(i, j) = rainbl(i, j) / dt;
    int vtyp = ivgtyp(i, j);
    int styp = isltyp(i, j);
    real shdf = vegfra(i, j) / real(100.);
    t1(i, j) = tsk(i, j);
    chk(i, j) = chs(i, j);
    shmin(i, j) = shdmin(i, j) / real(100.);
    shmax(i, j) = shdmax(i, j) / real(100.);
    real sneq = snow(i, j) * real(0.001);
    real snwh = snowh(i, j);
    sncovr(i, j) = snowc(i, j);
    if (frpcpn) {
      ffrozp(i, j) = sr(i, j);
    } else {
      ffrozp(i, j) = sfct <= real(273.15) ? real(1.0) : real(0.0);
    }
    slopetyp(i, j) = 1;
    if ((xland(i, j) - real(1.5)) >= real(0.)) {  // open water: no land-surface model
      vegtyp(i, j) = vtyp;
      soiltyp(i, j) = styp;
      return;
    }
    int ic;
    if (xice(i, j) >= xice_threshold) {
      ic = 1;
    } else if (vtyp == isice) {
      ic = -1;
    } else {
      ic = 0;
    }
    real dq = q2s * C::a23m4 / ((sfct - C::a4) * (sfct - C::a4));
    if (snow(i, j) > real(0.0)) {
      real sfctsno = sfct;
      real e2sat = real(611.2) * std::exp(real(6174.) * (real(1.) / real(273.15) - real(1.) / sfctsno));
      real q2sati = real(0.622) * e2sat / (sfcprs(i, j) - e2sat);
      q2sati = q2sati / (real(1.0) + q2sati);
      if (t1(i, j) > real(273.14)) {
        q2s = q2s * (real(1.) - snowc(i, j)) + q2sati * snowc(i, j);
        dq = dq * (real(1.) - snowc(i, j)) + q2sati * real(6174.) / (sfctsno * sfctsno) * snowc(i, j);
      } else {
        q2s = q2sati;
        dq = q2sati * real(6174.) / (sfctsno * sfctsno);
      }
      if (t1(i, j) > real(273.) && snowc(i, j) > real(0.) && swdown(i, j) > real(10.)) dq = dq * (real(1.) - snowc(i, j));
    }
    q2sat(i, j) = q2s;
    dqsdt2(i, j) = dq;
    tbot(i, j) = tmn(i, j);
    if (isurban == 1) {  // USGS land data
      if (vtyp == 25) shdf = real(0.0000);
      if (vtyp == 26) shdf = real(0.0000);
      if (vtyp == 27) shdf = real(0.0000);
    }
    if (styp == 14 && xice(i, j) == real(0.)) styp = 7;
    snoalb1(i, j) = snoalb(i, j);
    cmc(i, j) = canwat(i, j) / real(1000.);
    albbrd(i, j) = albbck(i, j);
    z0brd(i, j) = z0(i, j);
    embrd(i, j) = embck(i, j);
    snotime1(i, j) = snotime(i, j);
    ribb(i, j) = rib(i, j);
    for (int ns = 0; ns < nsoil; ++ns) {
      smc(i, j, ns) = smois(i, j, ns);
      stc(i, j, ns) = tslb(i, j, ns);
      swc(i, j, ns) = sh2o(i, j, ns);
    }
    if ((sneq != real(0.) && snwh == real(0.)) || (snwh <= sneq)) snwh = real(5.) * sneq;
    // sf_urban_physics = 0: urban (and LCZ) points use the urban vegetation category
    bool urb = vtyp == isurban;
    for (int n = 0; n < 11; ++n) urb = urb || vtyp == lcz[n];
    if (sf_urban_physics == 0 && urb) vtyp = isurban;
    if (rdlai2d) {
      if (shdf > real(0.0) && lai(i, j) <= real(0.0)) lai(i, j) = real(0.01);
      xlai(i, j) = lai(i, j);
    }
    vegtyp(i, j) = vtyp;
    soiltyp(i, j) = styp;
    shdfac(i, j) = shdf;
    sneqv(i, j) = sneq;
    snowhk(i, j) = snwh;
    ice(i, j) = ic;
    if (ic == 1) {  // sea ice: handled by seaice_noah
      for (int ns = 0; ns < nsoil; ++ns) sh2o(i, j, ns) = real(1.0);
      lai(i, j) = real(0.01);
      ice(i, j) = 1;
      return;
    }
    if (ic == 0) active(i, j) = 1;
  }
};

// SFLX_GLACIAL for land-ice points (ICE == -1); outputs in the SFLX scratch of the driver.
struct KLsmGlacial {
  CI2 ice;
  V2 lwdn, solnet, sfcprs, prcp, sfctmp, th2, q2k, q2sat, dqsdt2, ffrozp, albbrd, snoalb1, tbot, z0brd, emissi,
      embrd, t1, snowhk, sneqv, chk, snotime1, ribb, xlai, soilw, soilm, runoff1, runoff2, runoff3, eta, sheat,
      eta_kin, etp, ssoil, flx1, flx2, flx3, snomlt, sncovr, q1, albedok, z0k;
  V3 stc, smc, swc, smav;
  real dt, lvcoef;
  real sldpth[kMaxSoil];
  int nsoil;
  CPPWRF_HD void operator()(int i, int j) const {
    if (ice(i, j) != -1) return;
    soilm(i, j) = real(0.0);
    xlai(i, j) = real(0.01);
    runoff2(i, j) = real(0.0);
    runoff3(i, j) = real(0.0);
    for (int ns = 0; ns < nsoil; ++ns) {
      swc(i, j, ns) = real(1.0);
      smc(i, j, ns) = real(1.0);
      smav(i, j, ns) = real(1.0);
    }
    real stcl[kMaxSoil];
    for (int ns = 0; ns < nsoil; ++ns) stcl[ns] = stc(i, j, ns);
    real em = emissi(i, j), tt = t1(i, j), sh = snowhk(i, j), se = sneqv(i, j), sc = sncovr(i, j),
         st = snotime1(i, j), rb = ribb(i, j);
    noah_ice::glacial::Out o{};
    noah_ice::glacial::sflx_glacial(ffrozp(i, j), dt, nsoil, sldpth, lwdn(i, j), solnet(i, j), sfcprs(i, j),
                                    prcp(i, j), sfctmp(i, j), q2k(i, j), th2(i, j), q2sat(i, j), dqsdt2(i, j),
                                    albbrd(i, j), snoalb1(i, j), tbot(i, j), z0brd(i, j), em, embrd(i, j), tt, stcl,
                                    sh, se, chk(i, j), sc, st, rb, lvcoef, o);
    for (int ns = 0; ns < nsoil; ++ns) stc(i, j, ns) = stcl[ns];
    emissi(i, j) = em;
    t1(i, j) = tt;
    snowhk(i, j) = sh;
    sneqv(i, j) = se;
    sncovr(i, j) = sc;
    snotime1(i, j) = st;
    ribb(i, j) = rb;
    albedok(i, j) = o.albedo;
    z0k(i, j) = o.z0;
    eta(i, j) = o.eta;
    sheat(i, j) = o.sheat;
    eta_kin(i, j) = o.eta_kinematic;
    etp(i, j) = o.etp;
    ssoil(i, j) = o.ssoil;
    flx1(i, j) = o.flx1;
    flx2(i, j) = o.flx2;
    flx3(i, j) = o.flx3;
    snomlt(i, j) = o.snomlt;
    runoff1(i, j) = o.runoff1;
    q1(i, j) = o.q1;
  }
};

// lsm after SFLX / SFLX_GLACIAL (land points, ICE /= 1)
struct KLsmPost {
  CI2 ice;
  C2 xland;
  C2 lwdn, solnet, prcp, ffrozp, cmc, sneqv, snowhk, albedok, albbrd, z0brd, emissi, z0k, t1, sheat, etp, eta_kin,
      eta, ssoil, sncovr, snotime1, q1, flx1, flx2, flx3, soilw, soilm, runoff1, runoff2, snomlt, xlai, cqs2;
  C3 smc, stc, swc, smav;
  V2 lai, canwat, snow, snowh, albedo, albbck, z0, emiss, znt, tsk, hfx, potevp, qfx, lh, grdflx, snowc, chs2,
      snotime, qsfc, noahres, smstav, smstot, sfcrunoff, udrunoff, acsnow, acsnom, snopcx;
  V3 smois, tslb, sh2o, smcrel;
  real dt;
  int nsoil;
  CPPWRF_HD void operator()(int i, int j) const {
    using C = LsmConst;
    if ((xland(i, j) - real(1.5)) >= real(0.) || ice(i, j) == 1) return;
    const real fdtliw = dt / C::rowliw;
    const real fdtw = dt / (wrfc::xlv * wrfc::rhowater);
    lai(i, j) = xlai(i, j);
    canwat(i, j) = cmc(i, j) * real(1000.);
    snow(i, j) = sneqv(i, j) * real(1000.);
    snowh(i, j) = snowhk(i, j);
    albedo(i, j) = albedok(i, j);
    albbck(i, j) = albbrd(i, j);
    z0(i, j) = z0brd(i, j);
    emiss(i, j) = emissi(i, j);
    znt(i, j) = z0k(i, j);
    real tt = t1(i, j);
    tsk(i, j) = tt;
    hfx(i, j) = sheat(i, j);
    potevp(i, j) = potevp(i, j) + etp(i, j) * fdtw;
    qfx(i, j) = eta_kin(i, j);
    lh(i, j) = eta(i, j);
    grdflx(i, j) = ssoil(i, j);
    snowc(i, j) = sncovr(i, j);
    chs2(i, j) = cqs2(i, j);
    snotime(i, j) = snotime1(i, j);
    qsfc(i, j) = q1(i, j) / (real(1.0) - q1(i, j));
    for (int ns = 0; ns < nsoil; ++ns) {
      smois(i, j, ns) = smc(i, j, ns);
      tslb(i, j, ns) = stc(i, j, ns);
      sh2o(i, j, ns) = swc(i, j, ns);
    }
    noahres(i, j) = (solnet(i, j) + lwdn(i, j)) - sheat(i, j) + ssoil(i, j) - eta(i, j) -
                    (emissi(i, j) * wrfc::stbolt * ((tt * tt) * (tt * tt))) - flx1(i, j) - flx2(i, j) - flx3(i, j);
    smstav(i, j) = soilw(i, j);
    smstot(i, j) = soilm(i, j) * real(1000.);
    for (int ns = 0; ns < nsoil; ++ns) smcrel(i, j, ns) = smav(i, j, ns);
    sfcrunoff(i, j) = sfcrunoff(i, j) + runoff1(i, j) * dt * real(1000.0);
    udrunoff(i, j) = udrunoff(i, j) + runoff2(i, j) * dt * real(1000.0);
    if (ffrozp(i, j) > real(0.5)) acsnow(i, j) = acsnow(i, j) + prcp(i, j) * dt;
    if (snow(i, j) > real(0.)) {
      acsnom(i, j) = acsnom(i, j) + snomlt(i, j) * real(1000.);
      snopcx(i, j) = snopcx(i, j) - snomlt(i, j) / fdtliw;
    }
  }
};

// module_sf_noah_seaice_drv.F seaice_noah
struct KSeaiceNoah {
  C3 t3d, qv3d, p8w3d, dz8w;
  C2 glw, swdown, rainbl, snoalb2d, qgh, xice, rib, sr;
  V2 emiss, albedo, z02d, tsk, snow, snowc, snowh2d, chs, chs2, cqs2, znt, lh, hfx, qfx, potevp, grdflx, qsfc,
      acsnow, acsnom, snopcx, sfcrunoff, noahres;
  V3 tslb;
  noah_ice::seaice::Opts op;
  int thickness_opt;
  real thickness_default;
  real dt, xice_threshold;
  int nsoil;
  bool frpcpn;
  CPPWRF_HD void operator()(int i, int j) const {
    using C = LsmConst;
    const real fdtliw = dt / C::rowliw;
    const real fdtw = dt / (wrfc::xlv * wrfc::rhowater);
    if (xice(i, j) < xice_threshold) return;  // (ICEDEPTH/SNOWSI resets only for opt 1)
    real sithick = thickness_default;  // SEAICE_THICKNESS_OPT = 0 (opt 1 needs ICEDEPTH input)
    real sfctmp = t3d(i, j, 0);
    real t1 = tsk(i, j);
    real snoalb = snoalb2d(i, j);
    real albedosi = real(0), snowonsi = real(0);
    real emissi = emiss(i, j);
    real lwdn = glw(i, j) * emissi;
    real sneqv = snow(i, j) * real(0.001);
    real snowh = snowh2d(i, j);
    real sncovr = snowc(i, j);
    real solnet = swdown(i, j) * (real(1.) - albedo(i, j));
    real sfcprs = (p8w3d(i, j, 1) + p8w3d(i, j, 0)) * real(0.5);
    real psfc = p8w3d(i, j, 0);
    real q2 = qv3d(i, j, 0) / (real(1.0) + qv3d(i, j, 0));
    real apes = rpow(real(1.e5) / psfc, C::capa);
    real apelm = rpow(real(1.e5) / sfcprs, C::capa);
    real th2 = (sfctmp * apelm) / apes;
    real q2sat = qgh(i, j) / (real(1.0) + qgh(i, j));
    real dqsdt2 = q2sat * C::a23m4 / ((sfctmp - C::a4) * (sfctmp - C::a4));
    if (snow(i, j) > real(0.0)) {
      real sfctsno = sfctmp;
      real e2sat = real(611.2) * std::exp(real(6174.) * (real(1.) / real(273.15) - real(1.) / sfctsno));
      real q2sati = real(0.622) * e2sat / (sfcprs - e2sat);
      q2sati = q2sati / (real(1.0) + q2sati);
      if (t1 > real(273.14)) {
        q2sat = q2sat * (real(1.) - snowc(i, j)) + q2sati * snowc(i, j);
        dqsdt2 = dqsdt2 * (real(1.) - snowc(i, j)) + q2sati * real(6174.) / (sfctsno * sfctsno) * snowc(i, j);
      } else {
        q2sat = q2sati;
        dqsdt2 = q2sati * real(6174.) / (sfctsno * sfctsno);
      }
      if ((t1 > real(273.)) && (snowc(i, j) > real(0.0))) dqsdt2 = dqsdt2 * (real(1.) - snowc(i, j));
    }
    real prcp = rainbl(i, j) / dt;
    real ffrozp;
    if (frpcpn) {
      ffrozp = sr(i, j);
    } else {
      ffrozp = sfctmp <= real(273.15) ? real(1.0) : real(0.0);
    }
    real tbot = real(271.36);
    real z0brd = z02d(i, j);
    real stc[kMaxSoil];
    for (int ns = 0; ns < nsoil; ++ns) stc[ns] = tslb(i, j, ns);
    real ch = chs(i, j);
    real ribb = rib(i, j);
    noah_ice::seaice::Out o{};
    noah_ice::seaice::sflx_seaice(op, ffrozp, dt, nsoil, sithick, lwdn, solnet, sfcprs, prcp, sfctmp, q2, th2, q2sat,
                                  dqsdt2, snoalb, tbot, z0brd, emissi, t1, stc, snowh, sneqv, ch, albedosi,
                                  snowonsi, sncovr, ribb, o);
    albedo(i, j) = o.albedo;
    emiss(i, j) = emissi;
    tsk(i, j) = t1;
    z02d(i, j) = z0brd;
    snowh2d(i, j) = snowh;
    snowc(i, j) = sncovr;
    snow(i, j) = sneqv * real(1000.);
    for (int ns = 0; ns < nsoil; ++ns) tslb(i, j, ns) = stc[ns];
    znt(i, j) = o.z0;
    lh(i, j) = o.eta;
    hfx(i, j) = o.sheat;
    qfx(i, j) = o.eta_kinematic;
    potevp(i, j) = potevp(i, j) + o.etp * fdtw;
    grdflx(i, j) = o.ssoil;
    chs2(i, j) = cqs2(i, j);
    if (o.q1 > qsfc(i, j)) cqs2(i, j) = chs(i, j);
    qsfc(i, j) = o.q1 / (real(1.0) - o.q1);
    if (ffrozp > real(0.5)) acsnow(i, j) = acsnow(i, j) + prcp * dt;
    acsnom(i, j) = acsnom(i, j) + o.snomlt * real(1000.);
    snopcx(i, j) = snopcx(i, j) - o.snomlt / fdtliw;
    sfcrunoff(i, j) = sfcrunoff(i, j) + o.runoff1 * dt * real(1000.0);
    noahres(i, j) = (solnet + lwdn) - o.sheat + o.ssoil - o.eta - (emissi * wrfc::stbolt * ((t1 * t1) * (t1 * t1))) -
                    o.flx1 - o.flx2 - o.flx3;
    (void)thickness_opt;
  }
};

// LSMSCHEME tail of surface_driver: accumulations and SFCDIAGS (module_sf_sfcdiags.F)
struct KLsmDiag {
  C2 qfx, chs, hfx, lh, grdflx, tsk, qsfc, chs2, cqs2, psfc;
  C3 t3d, qv3d;
  V2 sfcevp, sfcexc, achfx, aclhf, acgrdflx, q2, t2, th2;
  real dtbl, dt;
  CPPWRF_HD void operator()(int i, int j) const {
    using namespace wrfc;
    sfcevp(i, j) = sfcevp(i, j) + qfx(i, j) * dtbl;
    sfcexc(i, j) = chs(i, j);
    achfx(i, j) = achfx(i, j) + hfx(i, j) * dt;
    aclhf(i, j) = aclhf(i, j) + lh(i, j) * dt;
    acgrdflx(i, j) = acgrdflx(i, j) + grdflx(i, j) * dt;
    // SFCDIAGS (ua_phys = .false.)
    real rho = psfc(i, j) / (r_d * tsk(i, j));
    if (cqs2(i, j) < real(1.e-5)) {
      q2(i, j) = qsfc(i, j);
    } else {
      q2(i, j) = qsfc(i, j) - qfx(i, j) / (rho * cqs2(i, j));
    }
    if (chs2(i, j) < real(1.e-5)) {
      t2(i, j) = tsk(i, j);
    } else {
      t2(i, j) = tsk(i, j) - hfx(i, j) / (rho * cp * chs2(i, j));
    }
    th2(i, j) = t2(i, j) * rpow(real(1.e5) / psfc(i, j), rcp);
    (void)t3d; (void)qv3d;
  }
};

// end of surface_driver: RAINBL reset, Q2 limited over land
struct KSfcEnd {
  V2 rainbl, q2;
  C2 xland;
  C3 qv;
  CPPWRF_HD void operator()(int i, int j) const {
    rainbl(i, j) = real(0.);
    if (xland(i, j) < real(1.5)) {
      real q = real(1.05) * qv(i, j, 0);
      q2(i, j) = q2(i, j) < q ? q2(i, j) : q;
    }
  }
};

// ---------------------------------------------------------------- pbl_driver
struct KPblPrep {
  C3 u_phy, v_phy, p8w;
  V3 u_phytmp, v_phytmp, rthblten, rublten, rvblten, rqcblten, rqvblten, rqiblten;
  V2 psfc;
  int nz;
  bool flag_qi;
  CPPWRF_HD void operator()(int i, int j) const {
    constexpr real u_frame = real(0.), v_frame = real(0.);
    for (int k = 0; k < nz; ++k) {
      v_phytmp(i, j, k) = v_phy(i, j, k) + v_frame;
      u_phytmp(i, j, k) = u_phy(i, j, k) + u_frame;
    }
    psfc(i, j) = p8w(i, j, 0);
    for (int k = 0; k < nz; ++k) {  // kts..min(kte+1,kde): the field has nz levels, level kde is never read
      rthblten(i, j, k) = real(0.);
      rublten(i, j, k) = real(0.);
      rvblten(i, j, k) = real(0.);
      rqcblten(i, j, k) = real(0.);
      rqvblten(i, j, k) = real(0.);
    }
    if (flag_qi)
      for (int k = 0; k < nz; ++k) rqiblten(i, j, k) = real(0.);
  }
};

}  // namespace

// ------------------------------------------------------------------ PhysicsSuite
PhysicsSuite::PhysicsSuite(const Grid& g, const SuiteOptions& o)
    : tsk(g, Stagger::M, "tsk"), xland(g, Stagger::M, "xland"), xice(g, Stagger::M, "xice"),
      lakemask(g, Stagger::M, "lakemask"), water_depth(g, Stagger::M, "water_depth"), mavail(g, Stagger::M, "mavail"),
      tmn(g, Stagger::M, "tmn"), vegfra(g, Stagger::M, "vegfra"), shdmin(g, Stagger::M, "shdmin"),
      shdmax(g, Stagger::M, "shdmax"), snoalb(g, Stagger::M, "snoalb"), albbck(g, Stagger::M, "albbck"),
      albedo(g, Stagger::M, "albedo"), emiss(g, Stagger::M, "emiss"), embck(g, Stagger::M, "embck"),
      znt(g, Stagger::M, "znt"), z0(g, Stagger::M, "z0"), ust(g, Stagger::M, "ust"), ustm(g, Stagger::M, "ustm"),
      mol(g, Stagger::M, "mol"), rmol(g, Stagger::M, "rmol"), pblh(g, Stagger::M, "pblh"), hfx(g, Stagger::M, "hfx"),
      qfx(g, Stagger::M, "qfx"), lh(g, Stagger::M, "lh"), grdflx(g, Stagger::M, "grdflx"), qsfc(g, Stagger::M, "qsfc"),
      flhc(g, Stagger::M, "flhc"), flqc(g, Stagger::M, "flqc"), wspd(g, Stagger::M, "wspd"), br(g, Stagger::M, "br"),
      fm(g, Stagger::M, "fm"), fh(g, Stagger::M, "fh"), zol(g, Stagger::M, "zol"), regime(g, Stagger::M, "regime"),
      u10(g, Stagger::M, "u10"), v10(g, Stagger::M, "v10"), t2(g, Stagger::M, "t2"), th2(g, Stagger::M, "th2"),
      q2(g, Stagger::M, "q2"), psfc(g, Stagger::M, "psfc"), ch(g, Stagger::M, "ch"), snow(g, Stagger::M, "snow"),
      snowh(g, Stagger::M, "snowh"), snowc(g, Stagger::M, "snowc"), canwat(g, Stagger::M, "canwat"),
      lai(g, Stagger::M, "lai"), snotime(g, Stagger::M, "snotime"), qz0(g, Stagger::M, "qz0"),
      rainbl(g, Stagger::M, "rainbl"), raincv(g, Stagger::M, "raincv"), rainncv(g, Stagger::M, "rainncv"),
      rainshv(g, Stagger::M, "rainshv"), sr(g, Stagger::M, "sr"), glw(g, Stagger::M, "glw"),
      swdown(g, Stagger::M, "swdown"), gsw(g, Stagger::M, "gsw"), smstav(g, Stagger::M, "smstav"),
      smstot(g, Stagger::M, "smstot"), sfcrunoff(g, Stagger::M, "sfcrunoff"), udrunoff(g, Stagger::M, "udrunoff"),
      acsnom(g, Stagger::M, "acsnom"), acsnow(g, Stagger::M, "acsnow"), snopcx(g, Stagger::M, "snopcx"),
      potevp(g, Stagger::M, "potevp"), noahres(g, Stagger::M, "noahres"), sfcevp(g, Stagger::M, "sfcevp"),
      sfcexc(g, Stagger::M, "sfcexc"), achfx(g, Stagger::M, "achfx"), aclhf(g, Stagger::M, "aclhf"),
      acgrdflx(g, Stagger::M, "acgrdflx"), dx2d(g, Stagger::M, "dx2d"), uoce(g, Stagger::M, "uoce"),
      voce(g, Stagger::M, "voce"), ctopo(g, Stagger::M, "ctopo"), ctopo2(g, Stagger::M, "ctopo2"),
      ck(g, Stagger::M, "ck"), cka(g, Stagger::M, "cka"), cd(g, Stagger::M, "cd"), cda(g, Stagger::M, "cda"),
      wstar(g, Stagger::M, "wstar"), delta(g, Stagger::M, "delta"), hol(g, Stagger::M, "hol"),
      ivgtyp(g, "ivgtyp"), isltyp(g, "isltyp"), kpbl(g, "kpbl"), pt(g), prep(g),
      psim(g, Stagger::M, "psim"), psih(g, Stagger::M, "psih"), gz1oz0(g, Stagger::M, "gz1oz0"),
      chklowq(g, Stagger::M, "chklowq"), qgh(g, Stagger::M, "qgh"), chs(g, Stagger::M, "chs"),
      chs2(g, Stagger::M, "chs2"), cqs2(g, Stagger::M, "cqs2"), cpm(g, Stagger::M, "cpm"),
      u_phytmp(g, Stagger::M, "u_phytmp"), v_phytmp(g, Stagger::M, "v_phytmp"), g_(g), gs_(g), o_(o) {
  if (o.sf_sfclay_physics != 1 || o.sf_surface_physics != 2 || o.bl_pbl_physics != 1)
    throw std::runtime_error("PhysicsSuite: only sf_sfclay_physics=1, sf_surface_physics=2, bl_pbl_physics=1");
  if (o.fractional_seaice != 0) throw std::runtime_error("PhysicsSuite: fractional_seaice=1 is not ported");
  if (o.seaice_thickness_opt != 0 || o.seaice_snowdepth_opt != 0)
    throw std::runtime_error("PhysicsSuite: seaice_thickness_opt/seaice_snowdepth_opt 1 need ICEDEPTH/SNOWSI input");
  if (o.num_soil_layers > noah::kMaxSoilLevels) throw std::runtime_error("PhysicsSuite: too many soil layers");
  if (g.nz > ysu::kMaxLevels) throw std::runtime_error("PhysicsSuite: nz exceeds the YSU column limit");
  gs_.nz = o.num_soil_layers;
  auto mk3 = [&](Field3& f, const Grid& gg, const char* n) { f = Field3(gg, Stagger::M, n); };
  mk3(tslb, gs_, "tslb"); mk3(smois, gs_, "smois"); mk3(sh2o, gs_, "sh2o"); mk3(smcrel, gs_, "smcrel");
  mk3(exch_h, g, "exch_h"); mk3(exch_m, g, "exch_m");
  mk3(n_stc, gs_, "n_stc"); mk3(n_smc, gs_, "n_smc"); mk3(n_swc, gs_, "n_swc"); mk3(n_et, gs_, "n_et");
  mk3(n_smav, gs_, "n_smav"); mk3(zero3_, g, "zero3");
  for (Field2* f : {&n_lwdn, &n_soldn, &n_solnet, &n_sfcprs, &n_prcp, &n_sfctmp, &n_th2, &n_q2, &n_q2sat, &n_dqsdt2,
                    &n_ffrozp, &n_shdfac, &n_shmin, &n_shmax, &n_albbrd, &n_snoalb1, &n_tbot, &n_z0brd, &n_emissi,
                    &n_embrd, &n_cmc, &n_t1, &n_snowhk, &n_sneqv, &n_chk, &n_cm, &n_albedok, &n_z0k, &n_xlai,
                    &n_snotime1, &n_ribb, &n_eta, &n_sheat, &n_eta_kin, &n_fdown, &n_ec, &n_edir, &n_ett, &n_esnow,
                    &n_drip, &n_dew, &n_beta, &n_etp, &n_ssoil, &n_flx1, &n_flx2, &n_flx3, &n_snomlt, &n_sncovr,
                    &n_runoff1, &n_runoff2, &n_runoff3, &n_rc, &n_pc, &n_rsmin, &n_soilw, &n_soilm, &n_q1})
    *f = Field2(g, Stagger::M, "noah_scratch");
  n_vegtyp = IntField2(g, "n_vegtyp"); n_soiltyp = IntField2(g, "n_soiltyp"); n_slopetyp = IntField2(g, "n_slopetyp");
  n_ice = IntField2(g, "n_ice"); n_active = IntField2(g, "n_active");
  dzs.assign(o.num_soil_layers, 0.0);
  // tables: read once (SOIL_VEG_GEN_PARM), sfclayrev psi tables built once (sf_sfclayrev_init)
  if (!noah::read_noah_tables(o.table_dir, o.mminlu, o.mminsl, noah_tables_))
    throw std::runtime_error("PhysicsSuite: cannot read the Noah tables from " + o.table_dir);
  sfclayrev::SFClayRevTables<real> h{};
  sfclayrev::init_sfclayrev_tables(h);
  d_sfclay_tables_ = static_cast<sfclayrev::SFClayRevTables<real>*>(dev_alloc(sizeof(h)));
  copy_to_dev(d_sfclay_tables_, &h, sizeof(h));
  // bldt = 0: every step; otherwise WRF stepbl = nint(bldt*60/dt), at least 1
  stepbl = o.bldt == 0 ? 1 : int(std::lround(o.bldt * 60. / g.dt));
  if (stepbl < 1) stepbl = 1;
}

PhysicsSuite::~PhysicsSuite() { dev_free(d_sfclay_tables_); }

void PhysicsSuite::part1(const dyn::State& s, const dyn::BaseState& b, dyn::Tendencies&) {
  ++itimestep;  // solve_em: grid%itimestep = grid%itimestep + 1
  run_phy_prep(s, b, prep, PhyPrepParams{o_.use_theta_m});
  surface_driver(s);
  pbl_driver(s);
}

void PhysicsSuite::part2(const dyn::State& s, const dyn::BaseState& b, dyn::Tendencies& t) {
  const CouplingFlags f{false, true};
  calculate_phy_tend(pt, s, b, g_, f);
  update_phy_ten(t, pt, s, g_, f);
  if (o_.use_theta_m == 1) conv_t_tendf_to_moist(t, s, g_);
}

void PhysicsSuite::after_rk(const dyn::State& s, const dyn::BaseState& b) {
  phy_prep_part2(pt, s, b, g_, CouplingFlags{false, true});
}

// module_surface_driver.F surface_driver for SFCLAYREVSCHEME + LSMSCHEME
void PhysicsSuite::surface_driver(const dyn::State& s) {
  const Range2D cols = g_.columns();
  const int nz = g_.nz;
  const real xice_threshold = o_.fractional_seaice == 0 ? real(0.5) : real(0.02);
  column_for(cols, KSfcAlways{rainbl.view(), qgh.view(), chs.view(), cpm.view(), chs2.view(), cqs2.view(),
                              raincv.cview(), rainncv.cview(), rainshv.cview(), u_phytmp.view(), v_phytmp.view(), nz});
  // run_param: itimestep == 1, bldt == 0 or stepbl == 1, mod(itimestep, stepbl) == 0
  const bool run_param = itimestep == 1 || o_.bldt == 0 || stepbl == 1 || itimestep % stepbl == 0;
  if (!run_param) return;
  const real dt = real(g_.dt);
  const real dtbl = o_.bldt == 0 ? dt : dt * real(stepbl);

  column_for(cols, KSfcPrep{psfc.view(), prep.p_hyd_w.cview(), prep.u_phy.cview(), prep.v_phy.cview(),
                            uoce.cview(), voce.cview(), u_phytmp.view(), v_phytmp.view(), nz});
  // SFCLAYREV (module_sf_sfclayrev.F wrapper: level kts of u_phytmp, v_phytmp, t_phy, qv,
  // P_PHY = grid%p_hyd, dz8w)
  {
    using namespace wrfc;
    sfclayrev::run_sfclayrev(
        slice0(u_phytmp), slice0(v_phytmp), slice0(prep.t_phy), slice0(s.q_2[s.P_QV]), slice0(prep.p_hyd),
        slice0(prep.dz8w), psfc.cview(), tsk.cview(), xland.cview(), dx2d.cview(), mavail.cview(), pblh.cview(),
        lakemask.cview(), water_depth.cview(), regime.view(), hfx.view(), qfx.view(), qsfc.view(), mol.view(),
        rmol.view(), gz1oz0.view(), wspd.view(), br.view(), psim.view(), psih.view(), fm.view(), fh.view(),
        znt.view(), zol.view(), ust.view(), cpm.view(), chs2.view(), cqs2.view(), chs.view(), flhc.view(),
        flqc.view(), qgh.view(), ustm.view(), lh.view(), u10.view(), v10.view(), th2.view(), t2.view(), q2.view(),
        ck.view(), cka.view(), cd.view(), cda.view(), d_sfclay_tables_, cols, o_.isfflx != 0, o_.shalwater_z0 != 0,
        o_.isftcflx, o_.iz0tlnd, o_.scm_force_flux != 0, cp, g, rcp, r_d, xlv, svp1, svp2, svp3, svpt0, ep_1,
        ep_2, karman, p1000mb);
  }
  column_for(cols, KCopy2{ch.view(), chs.cview()});  // ch = chs (EM_CORE)

  if (o_.isfflx != 0) {
    const int nsoil = o_.num_soil_layers;
    // ---- lsm (module_sf_noahdrv.F) ----
    KLsmPre pre{};
    pre.dz8w = prep.dz8w.cview(); pre.qv3d = s.q_2[s.P_QV].cview(); pre.p8w3d = prep.p_hyd_w.cview();
    pre.t3d = prep.t_phy.cview();
    pre.qgh = qgh.cview(); pre.glw = glw.cview(); pre.swdown = swdown.cview(); pre.rainbl = rainbl.cview();
    pre.vegfra = vegfra.cview(); pre.shdmin = shdmin.cview(); pre.shdmax = shdmax.cview(); pre.snoalb = snoalb.cview();
    pre.tmn = tmn.cview(); pre.xice = xice.cview(); pre.chs = chs.cview(); pre.rib = br.cview(); pre.sr = sr.cview();
    pre.ivgtyp = ivgtyp.cview(); pre.isltyp = isltyp.cview();
    pre.tsk = tsk.view(); pre.albedo = albedo.view(); pre.albbck = albbck.view(); pre.z0 = z0.view();
    pre.emiss = emiss.view(); pre.embck = embck.view(); pre.snow = snow.view(); pre.snowh = snowh.view();
    pre.snowc = snowc.view(); pre.canwat = canwat.view(); pre.snotime = snotime.view(); pre.lai = lai.view();
    pre.smstav = smstav.view(); pre.smstot = smstot.view(); pre.chklowq = chklowq.view(); pre.xland = xland.view();
    pre.smois = smois.view(); pre.tslb = tslb.view(); pre.sh2o = sh2o.view(); pre.smcrel = smcrel.view();
    pre.lwdn = n_lwdn.view(); pre.soldn = n_soldn.view(); pre.solnet = n_solnet.view(); pre.sfcprs = n_sfcprs.view();
    pre.prcp = n_prcp.view(); pre.sfctmp = n_sfctmp.view(); pre.th2 = n_th2.view(); pre.q2k = n_q2.view();
    pre.q2sat = n_q2sat.view(); pre.dqsdt2 = n_dqsdt2.view(); pre.ffrozp = n_ffrozp.view();
    pre.shdfac = n_shdfac.view(); pre.shmin = n_shmin.view(); pre.shmax = n_shmax.view();
    pre.albbrd = n_albbrd.view(); pre.snoalb1 = n_snoalb1.view(); pre.tbot = n_tbot.view();
    pre.z0brd = n_z0brd.view(); pre.emissi = n_emissi.view(); pre.embrd = n_embrd.view(); pre.cmc = n_cmc.view();
    pre.t1 = n_t1.view(); pre.snowhk = n_snowhk.view(); pre.sneqv = n_sneqv.view(); pre.chk = n_chk.view();
    pre.snotime1 = n_snotime1.view(); pre.ribb = n_ribb.view(); pre.xlai = n_xlai.view(); pre.soilw = n_soilw.view();
    pre.sncovr = n_sncovr.view();
    pre.stc = n_stc.view(); pre.smc = n_smc.view(); pre.swc = n_swc.view();
    pre.vegtyp = n_vegtyp.view(); pre.soiltyp = n_soiltyp.view(); pre.slopetyp = n_slopetyp.view();
    pre.ice = n_ice.view(); pre.active = n_active.view();
    pre.nsoil = nsoil; pre.itimestep = int(itimestep); pre.isurban = o_.isurban; pre.isice = o_.isice;
    pre.sf_urban_physics = 0;
    for (int n = 0; n < 11; ++n) pre.lcz[n] = noah_tables_.lcz[n];
    pre.myj = false; pre.frpcpn = true;  // SR is passed by first_rk_step_part1
    pre.rdlai2d = o_.rdlai2d; pre.dt = dtbl; pre.xice_threshold = xice_threshold;
    column_for(cols, pre);

    // SFLX on the land points (ICE == 0)
    noah::NoahKernel nk{};
    nk.lwdn = n_lwdn.cview(); nk.soldn = n_soldn.cview(); nk.solnet = n_solnet.cview(); nk.sfcprs = n_sfcprs.cview();
    nk.prcp = n_prcp.cview(); nk.sfctmp = n_sfctmp.cview(); nk.th2 = n_th2.cview(); nk.q2 = n_q2.cview();
    nk.sfcspd = n_q2.cview();  // SFCSPD: DUMMY in lsm (not used by SFLX)
    nk.q2sat = n_q2sat.cview(); nk.dqsdt2 = n_dqsdt2.cview(); nk.ffrozp = n_ffrozp.cview();
    nk.vegtyp = n_vegtyp.cview(); nk.soiltyp = n_soiltyp.cview(); nk.slopetyp = n_slopetyp.cview();
    nk.shdfac = n_shdfac.cview(); nk.shdmin = n_shmin.cview(); nk.shdmax = n_shmax.cview();
    nk.alb = n_albbrd.cview(); nk.snoalb = n_snoalb1.cview(); nk.tbot = n_tbot.cview();
    nk.z0brd = n_z0brd.view(); nk.embrd = n_embrd.view();
    nk.cmc = n_cmc.view(); nk.t1 = n_t1.view(); nk.stc = n_stc.view(); nk.smc = n_smc.view(); nk.sh2o = n_swc.view();
    nk.snowh = n_snowhk.view(); nk.sneqv = n_sneqv.view(); nk.albedo = n_albedok.view(); nk.ch = n_chk.view();
    nk.cm = n_cm.view(); nk.emissi = n_emissi.view(); nk.z0 = n_z0k.view(); nk.xlai = n_xlai.view();
    nk.snotime1 = n_snotime1.view(); nk.ribb = n_ribb.view();
    nk.eta = n_eta.view(); nk.sheat = n_sheat.view(); nk.eta_kinematic = n_eta_kin.view(); nk.fdown = n_fdown.view();
    nk.ec = n_ec.view(); nk.edir = n_edir.view(); nk.et = n_et.view(); nk.ett = n_ett.view(); nk.esnow = n_esnow.view();
    nk.drip = n_drip.view(); nk.dew = n_dew.view(); nk.beta = n_beta.view(); nk.etp = n_etp.view();
    nk.ssoil = n_ssoil.view(); nk.flx1 = n_flx1.view(); nk.flx2 = n_flx2.view(); nk.flx3 = n_flx3.view();
    nk.snomlt = n_snomlt.view(); nk.sncovr = n_sncovr.view(); nk.runoff1 = n_runoff1.view();
    nk.runoff2 = n_runoff2.view(); nk.runoff3 = n_runoff3.view(); nk.rc = n_rc.view(); nk.pc = n_pc.view();
    nk.rsmin = n_rsmin.view(); nk.soilw = n_soilw.view(); nk.soilm = n_soilm.view(); nk.q1 = n_q1.view();
    nk.smav = n_smav.view();
    nk.dt = dtbl; nk.zlvl = real(0); nk.nsoil = nsoil;
    for (int k = 0; k < nsoil; ++k) nk.sldpth[k] = real(dzs[k]);
    nk.isurban = o_.isurban;
    nk.tables = noah_tables_;
    nk.rdlai2d = o_.rdlai2d; nk.usemonalb = o_.usemonalb;
    nk.active = n_active.cview(); nk.alb_out = n_albbrd.view();
    noah::run_noah(nk, cols);

    // SFLX_GLACIAL on the land-ice points (ICE == -1)
    KLsmGlacial gl{};
    gl.ice = n_ice.cview(); gl.lwdn = n_lwdn.view(); gl.solnet = n_solnet.view(); gl.sfcprs = n_sfcprs.view();
    gl.prcp = n_prcp.view(); gl.sfctmp = n_sfctmp.view(); gl.th2 = n_th2.view(); gl.q2k = n_q2.view();
    gl.q2sat = n_q2sat.view(); gl.dqsdt2 = n_dqsdt2.view(); gl.ffrozp = n_ffrozp.view(); gl.albbrd = n_albbrd.view();
    gl.snoalb1 = n_snoalb1.view(); gl.tbot = n_tbot.view(); gl.z0brd = n_z0brd.view(); gl.emissi = n_emissi.view();
    gl.embrd = n_embrd.view(); gl.t1 = n_t1.view(); gl.snowhk = n_snowhk.view(); gl.sneqv = n_sneqv.view();
    gl.chk = n_chk.view(); gl.snotime1 = n_snotime1.view(); gl.ribb = n_ribb.view(); gl.xlai = n_xlai.view();
    gl.soilw = n_soilw.view(); gl.soilm = n_soilm.view(); gl.runoff1 = n_runoff1.view(); gl.runoff2 = n_runoff2.view();
    gl.runoff3 = n_runoff3.view(); gl.eta = n_eta.view(); gl.sheat = n_sheat.view(); gl.eta_kin = n_eta_kin.view();
    gl.etp = n_etp.view(); gl.ssoil = n_ssoil.view(); gl.flx1 = n_flx1.view(); gl.flx2 = n_flx2.view();
    gl.flx3 = n_flx3.view(); gl.snomlt = n_snomlt.view(); gl.sncovr = n_sncovr.view(); gl.q1 = n_q1.view();
    gl.albedok = n_albedok.view(); gl.z0k = n_z0k.view();
    gl.stc = n_stc.view(); gl.smc = n_smc.view(); gl.swc = n_swc.view(); gl.smav = n_smav.view();
    gl.dt = dtbl; gl.lvcoef = noah_tables_.lvcoef_data; gl.nsoil = nsoil;
    for (int k = 0; k < nsoil; ++k) gl.sldpth[k] = real(dzs[k]);
    column_for(cols, gl);

    KLsmPost post{};
    post.ice = n_ice.cview(); post.xland = xland.cview();
    post.lwdn = n_lwdn.cview(); post.solnet = n_solnet.cview(); post.prcp = n_prcp.cview();
    post.ffrozp = n_ffrozp.cview(); post.cmc = n_cmc.cview(); post.sneqv = n_sneqv.cview();
    post.snowhk = n_snowhk.cview(); post.albedok = n_albedok.cview(); post.albbrd = n_albbrd.cview();
    post.z0brd = n_z0brd.cview(); post.emissi = n_emissi.cview(); post.z0k = n_z0k.cview(); post.t1 = n_t1.cview();
    post.sheat = n_sheat.cview(); post.etp = n_etp.cview(); post.eta_kin = n_eta_kin.cview(); post.eta = n_eta.cview();
    post.ssoil = n_ssoil.cview(); post.sncovr = n_sncovr.cview(); post.snotime1 = n_snotime1.cview();
    post.q1 = n_q1.cview(); post.flx1 = n_flx1.cview(); post.flx2 = n_flx2.cview(); post.flx3 = n_flx3.cview();
    post.soilw = n_soilw.cview(); post.soilm = n_soilm.cview(); post.runoff1 = n_runoff1.cview();
    post.runoff2 = n_runoff2.cview(); post.snomlt = n_snomlt.cview(); post.xlai = n_xlai.cview();
    post.cqs2 = cqs2.cview();
    post.smc = n_smc.cview(); post.stc = n_stc.cview(); post.swc = n_swc.cview(); post.smav = n_smav.cview();
    post.lai = lai.view(); post.canwat = canwat.view(); post.snow = snow.view(); post.snowh = snowh.view();
    post.albedo = albedo.view(); post.albbck = albbck.view(); post.z0 = z0.view(); post.emiss = emiss.view();
    post.znt = znt.view(); post.tsk = tsk.view(); post.hfx = hfx.view(); post.potevp = potevp.view();
    post.qfx = qfx.view(); post.lh = lh.view(); post.grdflx = grdflx.view(); post.snowc = snowc.view();
    post.chs2 = chs2.view(); post.snotime = snotime.view(); post.qsfc = qsfc.view(); post.noahres = noahres.view();
    post.smstav = smstav.view(); post.smstot = smstot.view(); post.sfcrunoff = sfcrunoff.view();
    post.udrunoff = udrunoff.view(); post.acsnow = acsnow.view(); post.acsnom = acsnom.view();
    post.snopcx = snopcx.view();
    post.smois = smois.view(); post.tslb = tslb.view(); post.sh2o = sh2o.view(); post.smcrel = smcrel.view();
    post.dt = dtbl; post.nsoil = nsoil;
    column_for(cols, post);

    // ---- seaice_noah (module_sf_noah_seaice_drv.F), DT = model time step ----
    KSeaiceNoah si{};
    si.t3d = prep.t_phy.cview(); si.qv3d = s.q_2[s.P_QV].cview(); si.p8w3d = prep.p_hyd_w.cview();
    si.dz8w = prep.dz8w.cview();
    si.glw = glw.cview(); si.swdown = swdown.cview(); si.rainbl = rainbl.cview(); si.snoalb2d = snoalb.cview();
    si.qgh = qgh.cview(); si.xice = xice.cview(); si.rib = br.cview(); si.sr = sr.cview();
    si.emiss = emiss.view(); si.albedo = albedo.view(); si.z02d = z0.view(); si.tsk = tsk.view();
    si.snow = snow.view(); si.snowc = snowc.view(); si.snowh2d = snowh.view(); si.chs = chs.view();
    si.chs2 = chs2.view(); si.cqs2 = cqs2.view(); si.znt = znt.view(); si.lh = lh.view(); si.hfx = hfx.view();
    si.qfx = qfx.view(); si.potevp = potevp.view(); si.grdflx = grdflx.view(); si.qsfc = qsfc.view();
    si.acsnow = acsnow.view(); si.acsnom = acsnom.view(); si.snopcx = snopcx.view();
    si.sfcrunoff = sfcrunoff.view(); si.noahres = noahres.view(); si.tslb = tslb.view();
    si.op.albedo_opt = o_.seaice_albedo_opt; si.op.albedo_default = real(o_.seaice_albedo_default);
    si.op.snowdepth_opt = o_.seaice_snowdepth_opt; si.op.snowdepth_max = real(o_.seaice_snowdepth_max);
    si.op.snowdepth_min = real(o_.seaice_snowdepth_min);
    si.thickness_opt = o_.seaice_thickness_opt; si.thickness_default = real(o_.seaice_thickness_default);
    si.dt = dt; si.xice_threshold = xice_threshold; si.nsoil = nsoil; si.frpcpn = true;
    column_for(cols, si);

    column_for(cols, KLsmDiag{qfx.cview(), chs.cview(), hfx.cview(), lh.cview(), grdflx.cview(), tsk.cview(),
                              qsfc.cview(), chs2.cview(), cqs2.cview(), psfc.cview(), prep.t_phy.cview(),
                              s.q_2[s.P_QV].cview(), sfcevp.view(), sfcexc.view(), achfx.view(), aclhf.view(),
                              acgrdflx.view(), q2.view(), t2.view(), th2.view(), dtbl, dt});
  }
  column_for(cols, KSfcEnd{rainbl.view(), q2.view(), xland.cview(), s.q_2[s.P_QV].cview()});
}

// module_pbl_driver.F pbl_driver for YSUSCHEME
void PhysicsSuite::pbl_driver(const dyn::State& s) {
  const Range2D cols = g_.columns();
  const bool run_param = itimestep == 1 || o_.bldt == 0 || stepbl == 1 || itimestep % stepbl == 0;
  if (!run_param) return;
  const real dt = real(g_.dt);
  const real dtbl = o_.bldt == 0 ? dt : dt * real(stepbl);
  const bool flag_qc = s.P_QC >= 0, flag_qi = s.P_QI >= 0;
  column_for(cols, KPblPrep{prep.u_phy.cview(), prep.v_phy.cview(), prep.p_hyd_w.cview(), u_phytmp.view(),
                            v_phytmp.view(), pt.rthblten.view(), pt.rublten.view(), pt.rvblten.view(),
                            pt.rqcblten.view(), pt.rqvblten.view(), pt.rqiblten.view(), psfc.view(), g_.nz, flag_qi});
  using namespace wrfc;
  // YSU: U3D=u_phytmp, T3D=t_phy, P3D=p_phy (grid%p_hyd), P3DI=p8w (grid%p_hyd_w), PI3D=pi_phy,
  // PSIM=fm, PSIH=fhh, DT=dtbl
  ysu::run_ysu(u_phytmp.cview(), v_phytmp.cview(), prep.t_phy.cview(), s.q_2[s.P_QV].cview(),
               flag_qc ? s.q_2[s.P_QC].cview() : zero3_.cview(), flag_qi ? s.q_2[s.P_QI].cview() : zero3_.cview(),
               prep.p_hyd.cview(), prep.p_hyd_w.cview(), prep.pi_phy.cview(), prep.dz8w.cview(),
               pt.rthraten.cview(), psfc.cview(), znt.cview(), ust.cview(), fm.cview(), fh.cview(), xland.cview(),
               hfx.cview(), qfx.cview(), wspd.cview(), br.cview(), u10.view(), v10.view(), uoce.cview(),
               voce.cview(), ctopo.cview(), ctopo2.cview(), pt.rublten.view(), pt.rvblten.view(),
               pt.rthblten.view(), pt.rqvblten.view(), pt.rqcblten.view(), pt.rqiblten.view(), exch_h.view(),
               exch_m.view(), pblh.view(), kpbl.view(), wstar.view(), delta.view(), hol.view(), hol.view(),
               hol.view(), hol.view(), dtbl, cols, g_.nz, flag_qc, flag_qi, o_.ysu_topdown_pblmix == 1,
               /*have_ctopo=*/true, /*have_surface_diag=*/false, cp, g, rcp, r_d, rovg, ep_1, ep_2, karman,
               xlv, r_v);
}

}  // namespace cppwrf::phys::driver
