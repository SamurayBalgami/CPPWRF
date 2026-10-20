// PhysicsSuite::surface_driver (WRF phys/module_surface_driver.F)
#include "phys/driver/suite.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

#include "core/vm_dump.hpp"
#include "phys/driver/noah_ice.hpp"
#include "phys/driver/wrf_constants.hpp"
#include "phys/lsm_noah/noah.hpp"

#include "phys/driver/suite_internal.hpp"

namespace cppwrf::phys::driver {

using namespace suite_detail;

namespace {

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

// module_surface_driver.F: obs-nudging ratios after the surface layer (P_PHY = p_hyd)
struct KSfcRatios {
  C3 u_phytmp, v_phytmp, t_phy, p_phy;
  C2 u10, v10, th2;
  V2 uratx, vratx, tratx;
  int nx;
  CPPWRF_HD void operator()(int i, int j) const {
    uratx(i, j) = abs(u10(i, j)) > real(1.e-10) ? u_phytmp(i, j, 0) / u10(i, j) : real(1.2);
    vratx(i, j) = abs(v10(i, j)) > real(1.e-10) ? v_phytmp(i, j, 0) / v10(i, j) : real(1.2);
    // scalar libm: gfortran does not vectorise this loop of surface_driver
    tratx(i, j) = (t_phy(i, j, 0) * rpow(real(wrfc::p1000mb) * real(0.001) / (p_phy(i, j, 0) / real(1000.)),
                                     real(wrfc::rcp))) / th2(i, j);
  }
};
// RA = WSPD/UST**2 before the Noah call (Noah-CMAQ consistency)
struct KNoahRa {
  C2 wspd, ust;
  V2 ra;
  CPPWRF_HD void operator()(int i, int j) const { ra(i, j) = wspd(i, j) / (ust(i, j) * ust(i, j)); }
};
// rc2 = RC (RS) after SFLX and after SFLX_GLACIAL: SFLX_GLACIAL leaves the local RC alone, so a land-ice point
// gets the RC of the last SFLX point before it in lsm's loop order (j, then i); one thread walks the domain
struct KNoahRs {
  C2 rc;
  CI2 active, ice;
  V2 rs;
  int i0, i1, j0, j1;
  CPPWRF_HD void operator()(int, int) const {
    real last = real(0.);  // WRF: uninitialised before the first SFLX point
    for (int j = j0; j < j1; ++j)
      for (int i = i0; i < i1; ++i) {
        if (active(i, j) == 1) {
          last = rc(i, j);
          rs(i, j) = last;
        } else if (ice(i, j) == -1) {
          rs(i, j) = last;
        }
      }
  }
};

struct KCopy2 {
  V2 dst;
  C2 src;
  CPPWRF_HD void operator()(int i, int j) const { dst(i, j) = src(i, j); }
};

// module_surface_driver.F:2226-2232 after MYJSFC
struct KMyjsfcPost {
  V2 ustm, wspd, ch;
  C2 ust, chs;
  C3 u_phy, v_phy;
  CPPWRF_HD void operator()(int i, int j) const {
    ustm(i, j) = ust(i, j);
    real w = std::sqrt(u_phy(i, j, 0) * u_phy(i, j, 0) + v_phy(i, j, 0) * v_phy(i, j, 0));
    wspd(i, j) = w > real(0.001) ? w : real(0.001);
    ch(i, j) = chs(i, j);
  }
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
  C2 qgh, glw, swdown, rainbl, vegfra, shdmin, shdmax, snoalb, tmn, xice, chs, rib, sr, qz0;
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
    if (myj && (q2k(i, j) >= q2s * C::tresh) && (q2k(i, j) < qz0(i, j))) {  // QZ0 test of MYJ (not used with sfclayrev)
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
                                  dqsdt2, snoalb, tbot, z0brd, emissi, t1, stc, snowh, sneqv, ch, sncovr,
                                  ribb, o);
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

// NOAHMPSCHEME tail of surface_driver (module_surface_driver.F:3431-3481)
struct KNoahmpDiag {
  CI2 ivgtyp;
  C2 xice, qfx, chs, hfx, lh, grdflx, tsk, qsfc, chs2, cqs2, psfc, fvegxy, t2mvxy, t2mbxy, q2mvxy, q2mbxy;
  V2 chklowq, sfcevp, sfcexc, achfx, aclhf, acgrdflx, q2, t2, th2;
  real dtbl, dt, xice_threshold;
  int iswater, isice, isurban;
  int lcz[11];
  CPPWRF_HD void operator()(int i, int j) const {
    using namespace wrfc;
    chklowq(i, j) = real(1.0);
    sfcevp(i, j) = sfcevp(i, j) + qfx(i, j) * dtbl;
    sfcexc(i, j) = chs(i, j);
    achfx(i, j) = achfx(i, j) + hfx(i, j) * dt;
    aclhf(i, j) = aclhf(i, j) + lh(i, j) * dt;
    acgrdflx(i, j) = acgrdflx(i, j) + grdflx(i, j) * dt;
    const int vt = ivgtyp(i, j);
    bool is_lcz = false;
    for (int n = 0; n < 11; ++n) is_lcz = is_lcz || (vt == lcz[n]);
    if (vt == iswater || (vt == isice && xice(i, j) >= xice_threshold)) {
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
    } else if (vt == isurban || is_lcz || (vt == isice && xice(i, j) < xice_threshold)) {
      q2(i, j) = q2mbxy(i, j);
      t2(i, j) = t2mbxy(i, j);
      th2(i, j) = t2(i, j) * rpow(real(1.e5) / psfc(i, j), rcp);
    } else {
      t2(i, j) = fvegxy(i, j) * t2mvxy(i, j) + (real(1.) - fvegxy(i, j)) * t2mbxy(i, j);
      q2(i, j) = fvegxy(i, j) * q2mvxy(i, j) + (real(1.) - fvegxy(i, j)) * q2mbxy(i, j);
      th2(i, j) = t2(i, j) * rpow(real(1.e5) / psfc(i, j), rcp);
    }
  }
};

// sst_update = 1, every surface_driver call (module_surface_driver.F:1697-1850, lakemodel = 0): ice and water
// points follow the SEAICE read from wrflowinp; TSK and the top soil level take the SST over water
struct KSstUpdate {
  C2 sst, xice;
  V2 xicem, tsk_save, tsk, xland, vegfra, tmn, albedo, albbck, emiss, embck, snow, snowc, snowh;
  I2 ivgtyp, isltyp;
  V3 tslb, smois, sh2o;
  real xice_threshold, seaice_albedo_default;
  int nsoil, isice, iswater;
  bool fractional;
  CPPWRF_HD void operator()(int i, int j) const {
    const real x = xice(i, j);
    if (fractional) {
      if ((x != xicem(i, j)) && (xicem(i, j) > xice_threshold)) {
        // ALBEDO and EMISS hold the averages of the previous ice fraction XICEM
        albedo(i, j) = real(0.08) + x / xicem(i, j) * (albedo(i, j) - real(0.08));
        emiss(i, j) = real(0.98) + x / xicem(i, j) * (emiss(i, j) - real(0.98));
        tsk(i, j) = tsk_save(i, j) * x + (real(1.) - x) * sst(i, j);  // the ice part of the last TSK
      }
    }
    if (xland(i, j) > real(1.5) && x >= xice_threshold && xicem(i, j) < xice_threshold) {  // water turns to ice
      xicem(i, j) = x;
      xland(i, j) = real(1.);
      ivgtyp(i, j) = isice;
      isltyp(i, j) = 16;
      vegfra(i, j) = real(0.);
      tmn(i, j) = real(271.4);
      albedo(i, j) = seaice_albedo_default * x + real(0.08) * (real(1.0) - x);  // seaice_albedo_opt 0
      albbck(i, j) = seaice_albedo_default;
      emiss(i, j) = real(0.98) * x + real(0.98) * (real(1.0) - x);
      embck(i, j) = real(0.98);
      for (int nk = 0; nk < nsoil; ++nk) {
        tslb(i, j, nk) = tsk(i, j);
        smois(i, j, nk) = real(1.0);
        sh2o(i, j, nk) = real(0.0);
      }
    }
    if (xland(i, j) > real(1.5)) {
      if (sst(i, j) < real(350.) && sst(i, j) > real(250.)) {
        tsk(i, j) = sst(i, j);
        tslb(i, j, 0) = sst(i, j);
      }
    }
    if (xland(i, j) < real(1.5) && xicem(i, j) >= xice_threshold && x < xice_threshold) {  // ice turns to water
      xicem(i, j) = x;
      xland(i, j) = real(2.);
      ivgtyp(i, j) = iswater;
      isltyp(i, j) = 14;
      vegfra(i, j) = real(0.);
      snow(i, j) = real(0.);
      snowc(i, j) = real(0.);
      snowh(i, j) = real(0.);
      tmn(i, j) = sst(i, j);
      albedo(i, j) = real(0.08);
      albbck(i, j) = real(0.08);
      emiss(i, j) = real(0.98);
      embck(i, j) = real(0.98);
      for (int nk = 0; nk < nsoil; ++nk) {
        tslb(i, j, nk) = sst(i, j);
        smois(i, j, nk) = real(1.0);
        sh2o(i, j, nk) = real(1.0);
      }
    }
    xicem(i, j) = x;
    tsk_save(i, j) = tsk(i, j);
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

// ---------------------------------------------------------------- fractional_seaice = 1
// the sea-ice points of the fractional treatment: XICE_THRESHOLD <= XICE <= 1
CPPWRF_HD inline bool frac_ice(real xice, real xice_threshold) { return xice >= xice_threshold && xice <= real(1.0); }

// get_local_ice_tsk: SST under ice raised to 271.4 (and capped in the first steps), TSK_SEA = SST and the ice
// portion's skin temperature from the blended TSK; land and open-water points keep TSK
CPPWRF_HD inline void local_ice_tsk(real xice, real xice_threshold, real tsk, long itimestep, bool sfc_layer_values,
                                    bool tice2tsk_if2cold, real& sst, real& tsk_sea, real& tsk_ice) {
  const real tice_min = real(221.4);
  if (xice >= xice_threshold && xice <= real(1.0)) {
    if (sst < real(271.4)) sst = real(271.4);
    if (sfc_layer_values) {
      if (sst > real(273.) && itimestep <= 3) {
        if (xice >= real(0.6)) {
          sst = real(271.4);
        } else if (xice >= real(0.4)) {
          sst = real(273.);
        } else if (xice >= real(0.2) && sst > real(275.)) {
          sst = real(275.);
        } else if (sst > real(278.)) {
          sst = real(278.);
        }
      }
    }
    tsk_sea = sst;
    if (tice2tsk_if2cold) {
      tsk_ice = tsk < real(273.15) ? tsk : real(273.15);
    } else {
      tsk_ice = (tsk - (real(1.0) - xice) * sst) / xice;
      if (tsk_ice < tice_min) tsk_ice = tice_min;
    }
    if ((xice < real(0.2)) && (tsk < real(253.15))) tsk_ice = real(253.15);
    if ((xice < real(0.1)) && (tsk < real(263.15))) tsk_ice = real(263.15);
  } else {
    tsk_sea = tsk;
    tsk_ice = tsk;
  }
}

// sfclayrev_seaice_wrapper before the two SFCLAYREV calls: local ice TSK, the *_HOLD copies (the open-water
// call starts from them) and the open-water surface (XLAND 2, MAVAIL 1, TSK = SST)
struct KSfclayrevSeaPre {
  C2 xice, tsk, xland, mavail, br, chs2, chs, cpm, cqs2, flhc, flqc, gz1oz0, hfx, lh, mol, psih, psim, fh, fm, qfx, qgh,
      regime, rmol, ust, wspd, zol, znt, qsfc;
  V2 sst, tsk_local, tsk_sea, xland_sea, mavail_sea, br_sea, chs2_sea, chs_sea, cpm_sea, cqs2_sea, flhc_sea, flqc_sea,
      gz1oz0_sea, hfx_sea, lh_sea, mol_sea, psih_sea, psim_sea, fh_sea, fm_sea, qfx_sea, qgh_sea, regime_sea,
      regime_hold, rmol_sea, ust_sea, wspd_sea, zol_sea, znt_sea, qsfc_sea;
  real xice_threshold;
  long itimestep;
  bool tice2tsk_if2cold;
  CPPWRF_HD void operator()(int i, int j) const {
    real ss = sst(i, j), ts, tl;
    local_ice_tsk(xice(i, j), xice_threshold, tsk(i, j), itimestep, true, tice2tsk_if2cold, ss, ts, tl);
    br_sea(i, j) = br(i, j); chs2_sea(i, j) = chs2(i, j); chs_sea(i, j) = chs(i, j); cpm_sea(i, j) = cpm(i, j);
    cqs2_sea(i, j) = cqs2(i, j); flhc_sea(i, j) = flhc(i, j); flqc_sea(i, j) = flqc(i, j);
    gz1oz0_sea(i, j) = gz1oz0(i, j); hfx_sea(i, j) = hfx(i, j); lh_sea(i, j) = lh(i, j); mol_sea(i, j) = mol(i, j);
    psih_sea(i, j) = psih(i, j); psim_sea(i, j) = psim(i, j); fh_sea(i, j) = fh(i, j); fm_sea(i, j) = fm(i, j);
    qfx_sea(i, j) = qfx(i, j); qgh_sea(i, j) = qgh(i, j); regime_sea(i, j) = regime(i, j);
    regime_hold(i, j) = regime(i, j); rmol_sea(i, j) = rmol(i, j); ust_sea(i, j) = ust(i, j);
    wspd_sea(i, j) = wspd(i, j); zol_sea(i, j) = zol(i, j); qsfc_sea(i, j) = qsfc(i, j);
    znt_sea(i, j) = znt(i, j);  // the wrapper's ZNT_SEA = 0.0001 is replaced by ZNT_HOLD before the call
    if (frac_ice(xice(i, j), xice_threshold)) {
      xland_sea(i, j) = real(2.);
      mavail_sea(i, j) = real(1.);
      ts = ss;
      if (ss < real(271.4)) {
        ss = real(271.4);
        ts = ss;
      }
    } else {
      xland_sea(i, j) = xland(i, j);
      mavail_sea(i, j) = mavail(i, j);
      ts = tl;
    }
    sst(i, j) = ss;
    tsk_sea(i, j) = ts;
    tsk_local(i, j) = tl;
  }
};

// sfclayrev_seaice_wrapper after the open-water call: ice/open-water averages (REGIME of the ice call when
// XICE >= 0.5, else of the open-water call)
struct KSfclayrevSeaBlend {
  C2 xice, br_sea, gz1oz0_sea, mol_sea, psih_sea, psim_sea, fh_sea, fm_sea, regime_hold, rmol_sea, ust_sea, wspd_sea,
      zol_sea, cd_sea, cda_sea, ck_sea, cka_sea, q2_sea, t2_sea, th2_sea, u10_sea, ustm_sea, v10_sea;
  V2 br, gz1oz0, mol, psih, psim, fh, fm, regime, rmol, ust, wspd, zol, cd, cda, ck, cka, q2, t2, th2, u10, ustm, v10;
  real xice_threshold;
  CPPWRF_HD void operator()(int i, int j) const {
    const real x = xice(i, j);
    if (!frac_ice(x, xice_threshold)) return;
    auto bl = [&](V2 f, C2 fs) { f(i, j) = (f(i, j) * x) + ((real(1.0) - x) * fs(i, j)); };
    bl(br, br_sea); bl(gz1oz0, gz1oz0_sea); bl(mol, mol_sea); bl(psih, psih_sea); bl(psim, psim_sea);
    bl(fh, fh_sea); bl(fm, fm_sea);
    if (x >= real(0.5)) regime(i, j) = regime_hold(i, j);
    bl(rmol, rmol_sea); bl(ust, ust_sea); bl(wspd, wspd_sea); bl(zol, zol_sea);
    bl(cd, cd_sea); bl(cda, cda_sea); bl(ck, ck_sea); bl(cka, cka_sea);
    bl(q2, q2_sea); bl(t2, t2_sea); bl(th2, th2_sea); bl(u10, u10_sea);
    bl(ustm, ustm_sea);  // WRF: USTM_SEA is an uninitialised local (zero here)
    bl(v10, v10_sea);
  }
};

// myjsfc_seaice_wrapper before the two MYJSFC calls: TSK becomes the local ice TSK, QSFC the ice portion's,
// *_HOLD copies and the open-water surface (XLAND 2, MAVAIL 1, Z0 = 0.0001, TSK = SST)
struct KMyjsfcSeaPre {
  C2 xice, xland, mavail, z0, znt, qz0, thz0, uz0, vz0, ust, pblh, rmol, akhs, akms, hfx, qfx, lh;
  C3 pint;
  CI2 lowlyr;
  V2 tsk, qsfc, sst, tsk_sea, xland_sea, mavail_sea, znt_sea, z0base_sea, qsfc_sea, qz0_sea, thz0_sea, uz0_sea,
      vz0_sea, ust_sea, pblh_sea, rmol_sea, akhs_sea, akms_sea, hfx_sea, qfx_sea, lh_sea;
  real xice_threshold;
  long itimestep;
  bool tice2tsk_if2cold;
  CPPWRF_HD void operator()(int i, int j) const {
    // module_sf_myjsfc.F parameters
    constexpr real a2s = real(17.2693882), a3s = real(273.16), a4s = real(35.86);
    constexpr real pq0sea = real(379.90516) * real(0.98);
    real ss = sst(i, j), ts, tl;
    local_ice_tsk(xice(i, j), xice_threshold, tsk(i, j), itimestep, true, tice2tsk_if2cold, ss, ts, tl);
    tsk(i, j) = tl;
    const bool ice = frac_ice(xice(i, j), xice_threshold);
    const real psfc = pint(i, j, lowlyr(i, j) - 1);
    if (ice) {
      const real qs = pq0sea / psfc * std::exp(a2s * (tsk(i, j) - a3s) / (tsk(i, j) - a4s));
      qsfc(i, j) = (qsfc(i, j) - (real(1.0) - xice(i, j)) * qs) / xice(i, j);
      hfx_sea(i, j) = hfx(i, j);
      qfx_sea(i, j) = qfx(i, j);
      lh_sea(i, j) = lh(i, j);
    }
    qz0_sea(i, j) = qz0(i, j); thz0_sea(i, j) = thz0(i, j); uz0_sea(i, j) = uz0(i, j); vz0_sea(i, j) = vz0(i, j);
    ust_sea(i, j) = ust(i, j); pblh_sea(i, j) = pblh(i, j); rmol_sea(i, j) = rmol(i, j);
    akhs_sea(i, j) = akhs(i, j); akms_sea(i, j) = akms(i, j);
    if (ice) {
      xland_sea(i, j) = real(2.);
      mavail_sea(i, j) = real(1.);
      znt_sea(i, j) = real(0.0001);
      z0base_sea(i, j) = znt_sea(i, j);
      if (ss < real(271.4)) ss = real(271.4);
      ts = ss;
      qsfc_sea(i, j) = pq0sea / psfc * std::exp(a2s * (ts - a3s) / (ts - a4s));
    } else {
      xland_sea(i, j) = xland(i, j);
      mavail_sea(i, j) = mavail(i, j);
      znt_sea(i, j) = znt(i, j);
      z0base_sea(i, j) = z0(i, j);
      ts = tsk(i, j);
      qsfc_sea(i, j) = qsfc(i, j);
    }
    sst(i, j) = ss;
    tsk_sea(i, j) = ts;
  }
};

// myjsfc_seaice_wrapper after the open-water call
struct KMyjsfcSeaBlend {
  C2 xice, ct_sea, pshltr_sea, qshltr_sea, q2_sea, q10_sea, th2_sea, th10_sea, tshltr_sea, t2_sea, u10_sea, v10_sea,
      thz0_sea, uz0_sea, vz0_sea, ust_sea, pblh_sea, rmol_sea, akhs_sea, akms_sea;
  V2 ct, pshltr, qshltr, q2, q10, th2, th10, tshltr, t2, u10, v10, u10e, v10e, thz0, uz0, vz0, ust, pblh, rmol, akhs,
      akms;
  real xice_threshold;
  CPPWRF_HD void operator()(int i, int j) const {
    const real x = xice(i, j);
    if (!frac_ice(x, xice_threshold)) return;
    auto bl = [&](V2 f, C2 fs) { f(i, j) = f(i, j) * x + (real(1.0) - x) * fs(i, j); };
    bl(ct, ct_sea); bl(pshltr, pshltr_sea); bl(qshltr, qshltr_sea); bl(q2, q2_sea); bl(q10, q10_sea);
    bl(th2, th2_sea); bl(th10, th10_sea); bl(tshltr, tshltr_sea); bl(t2, t2_sea); bl(u10, u10_sea);
    bl(v10, v10_sea);
    u10e(i, j) = u10(i, j);
    v10e(i, j) = v10(i, j);
    bl(thz0, thz0_sea); bl(uz0, uz0_sea); bl(vz0, vz0_sea); bl(ust, ust_sea); bl(pblh, pblh_sea);
    bl(rmol, rmol_sea); bl(akhs, akhs_sea); bl(akms, akms_sea);
  }
};

// LSM / Noah-MP entry: ALBEDO and EMISS of the ice portion (open water 0.08 and 0.98)
struct KSeaiceLsmPre {
  C2 xice;
  V2 albedo, emiss;
  real xice_threshold;
  CPPWRF_HD void operator()(int i, int j) const {
    const real x = xice(i, j);
    if (!(x >= xice_threshold && x <= real(1))) return;
    albedo(i, j) = (albedo(i, j) - (real(1.) - x) * real(0.08)) / x;
    emiss(i, j) = (emiss(i, j) - (real(1.) - x) * real(0.98)) / x;
  }
};

// after seaice_noah: back to the ice/open-water averages (the surface layer ran for both portions)
struct KSeaiceLsmPost {
  C2 xice, flhc_sea, flqc_sea, cpm_sea, cqs2_sea, chs2_sea, chs_sea, qsfc_sea, qgh_sea, qz0_sea, hfx_sea, qfx_sea,
      lh_sea, tsk_sea;
  V2 albedo, emiss, flhc, flqc, cpm, cqs2, chs2, chs, qsfc, qgh, qz0, hfx, qfx, lh, tsk_save, tsk;
  real xice_threshold;
  CPPWRF_HD void operator()(int i, int j) const {
    const real x = xice(i, j);
    if (!frac_ice(x, xice_threshold)) return;
    albedo(i, j) = (albedo(i, j) * x) + ((real(1.0) - x) * real(0.08));
    emiss(i, j) = (emiss(i, j) * x) + ((real(1.0) - x) * real(0.98));
    auto bl = [&](V2 f, C2 fs) { f(i, j) = (f(i, j) * x) + ((real(1.0) - x) * fs(i, j)); };
    bl(flhc, flhc_sea); bl(flqc, flqc_sea); bl(cpm, cpm_sea); bl(cqs2, cqs2_sea); bl(chs2, chs2_sea);
    bl(chs, chs_sea); bl(qsfc, qsfc_sea); bl(qgh, qgh_sea); bl(qz0, qz0_sea); bl(hfx, hfx_sea); bl(qfx, qfx_sea);
    bl(lh, lh_sea);
    tsk_save(i, j) = tsk(i, j);
    bl(tsk, tsk_sea);
  }
};

}  // namespace

// module_surface_driver.F surface_driver for SFCLAYREVSCHEME / MYJSFCSCHEME + LSMSCHEME / NOAHMPSCHEME
void PhysicsSuite::surface_driver(const dyn::State& s, const dyn::BaseState& b) {
  const Range2D cols = g_.columns();
  const int nz = g_.nz;
  const real xice_threshold = o_.fractional_seaice == 0 ? real(0.5) : real(0.02);
  column_for(cols, KSfcAlways{rainbl.view(), qgh.view(), chs.view(), cpm.view(), chs2.view(), cqs2.view(),
                              raincv.cview(), rainncv.cview(), rainshv.cview(), u_phytmp.view(), v_phytmp.view(), nz});
  if (o_.sst_update == 1)
    column_for(cols, KSstUpdate{sst.cview(), xice.cview(), xicem.view(), tsk_save.view(), tsk.view(), xland.view(),
                                vegfra.view(), tmn.view(), albedo.view(), albbck.view(), emiss.view(), embck.view(),
                                snow.view(), snowc.view(), snowh.view(), ivgtyp.view(), isltyp.view(), tslb.view(),
                                smois.view(), sh2o.view(), xice_threshold, real(o_.seaice_albedo_default),
                                o_.num_soil_layers, o_.isice, o_.iswater, o_.fractional_seaice == 1});
  // run_param: itimestep == 1, bldt == 0 or stepbl == 1, mod(itimestep, stepbl) == 0
  const bool run_param = itimestep == 1 || o_.bldt == 0 || stepbl == 1 || itimestep % stepbl == 0;
  if (!run_param) return;
  const real dt = real(g_.dt);
  const real dtbl = o_.bldt == 0 ? dt : dt * real(stepbl);

  column_for(cols, KSfcPrep{psfc.view(), prep.p_hyd_w.cview(), prep.u_phy.cview(), prep.v_phy.cview(),
                            uoce.cview(), voce.cview(), u_phytmp.view(), v_phytmp.view(), nz});
  const bool frac = o_.fractional_seaice == 1;
  if (o_.sf_sfclay_physics == 1) {
    // SFCLAYREV (module_sf_sfclayrev.F wrapper: level kts of u_phytmp, v_phytmp, t_phy, qv,
    // P_PHY = grid%p_hyd, dz8w)
    using namespace wrfc;
    auto sfclayrev_call = [&](C2 tsk_v, C2 xland_v, C2 mavail_v, V2 regime_v, V2 hfx_v, V2 qfx_v, V2 qsfc_v, V2 mol_v,
                              V2 rmol_v, V2 gz1oz0_v, V2 wspd_v, V2 br_v, V2 psim_v, V2 psih_v, V2 fm_v, V2 fh_v,
                              V2 znt_v, V2 zol_v, V2 ust_v, V2 cpm_v, V2 chs2_v, V2 cqs2_v, V2 chs_v, V2 flhc_v,
                              V2 flqc_v, V2 qgh_v, V2 ustm_v, V2 lh_v, V2 u10_v, V2 v10_v, V2 th2_v, V2 t2_v, V2 q2_v,
                              V2 ck_v, V2 cka_v, V2 cd_v, V2 cda_v) {
      sfclayrev::run_sfclayrev(
          slice0(u_phytmp), slice0(v_phytmp), slice0(prep.t_phy), slice0(s.q_2[s.P_QV]), slice0(prep.p_hyd),
          slice0(prep.dz8w), psfc.cview(), tsk_v, xland_v, dx2d.cview(), mavail_v, pblh.cview(), lakemask.cview(),
          water_depth.cview(), regime_v, hfx_v, qfx_v, qsfc_v, mol_v, rmol_v, gz1oz0_v, wspd_v, br_v, psim_v,
          psih_v, fm_v, fh_v, znt_v, zol_v, ust_v, cpm_v, chs2_v, cqs2_v, chs_v, flhc_v, flqc_v, qgh_v, ustm_v, lh_v,
          u10_v, v10_v, th2_v, t2_v, q2_v, ck_v, cka_v, cd_v, cda_v, d_sfclay_tables_, cols, o_.isfflx != 0,
          o_.shalwater_z0 != 0, o_.isftcflx, o_.iz0tlnd, o_.scm_force_flux != 0, cp, g, rcp, r_d, xlv, svp1, svp2,
          svp3, svpt0, ep_1, ep_2, karman, p1000mb, var2d.cview(), o_.kim_tofd, real(o_.tofd_factor));
    };
    if (!frac) {
      sfclayrev_call(tsk.cview(), xland.cview(), mavail.cview(), regime.view(), hfx.view(), qfx.view(), qsfc.view(),
                     mol.view(), rmol.view(), gz1oz0.view(), wspd.view(), br.view(), psim.view(), psih.view(),
                     fm.view(), fh.view(), znt.view(), zol.view(), ust.view(), cpm.view(), chs2.view(), cqs2.view(),
                     chs.view(), flhc.view(), flqc.view(), qgh.view(), ustm.view(), lh.view(), u10.view(), v10.view(),
                     th2.view(), t2.view(), q2.view(), ck.view(), cka.view(), cd.view(), cda.view());
    } else {  // sfclayrev_seaice_wrapper: the ice portion, then the open-water portion of the sea-ice points
      SeaIceFields& w = *sea_;
      column_for(cols, KSfclayrevSeaPre{
          xice.cview(), tsk.cview(), xland.cview(), mavail.cview(), br.cview(), chs2.cview(), chs.cview(),
          cpm.cview(), cqs2.cview(), flhc.cview(), flqc.cview(), gz1oz0.cview(), hfx.cview(), lh.cview(),
          mol.cview(), psih.cview(), psim.cview(), fh.cview(), fm.cview(), qfx.cview(), qgh.cview(), regime.cview(),
          rmol.cview(), ust.cview(), wspd.cview(), zol.cview(), znt.cview(), qsfc.cview(), sst.view(),
          w.tsk_local.view(), w.tsk_sea.view(), w.xland_sea.view(), w.mavail_sea.view(), w.br_sea.view(),
          w.chs2_sea.view(), w.chs_sea.view(), w.cpm_sea.view(), w.cqs2_sea.view(), w.flhc_sea.view(),
          w.flqc_sea.view(), w.gz1oz0_sea.view(), w.hfx_sea.view(), w.lh_sea.view(), w.mol_sea.view(),
          w.psih_sea.view(), w.psim_sea.view(), w.fh_sea.view(), w.fm_sea.view(), w.qfx_sea.view(),
          w.qgh_sea.view(), w.regime_sea.view(), w.regime_hold.view(), w.rmol_sea.view(), w.ust_sea.view(),
          w.wspd_sea.view(), w.zol_sea.view(), w.znt_sea.view(), w.qsfc_sea.view(), xice_threshold, itimestep,
          o_.tice2tsk_if2cold});
      sfclayrev_call(w.tsk_local.cview(), xland.cview(), mavail.cview(), regime.view(), hfx.view(), qfx.view(),
                     qsfc.view(), mol.view(), rmol.view(), gz1oz0.view(), wspd.view(), br.view(), psim.view(),
                     psih.view(), fm.view(), fh.view(), znt.view(), zol.view(), ust.view(), cpm.view(), chs2.view(),
                     cqs2.view(), chs.view(), flhc.view(), flqc.view(), qgh.view(), ustm.view(), lh.view(), u10.view(),
                     v10.view(), th2.view(), t2.view(), q2.view(), ck.view(), cka.view(), cd.view(), cda.view());
      sfclayrev_call(w.tsk_sea.cview(), w.xland_sea.cview(), w.mavail_sea.cview(), w.regime_sea.view(),
                     w.hfx_sea.view(), w.qfx_sea.view(), w.qsfc_sea.view(), w.mol_sea.view(), w.rmol_sea.view(),
                     w.gz1oz0_sea.view(), w.wspd_sea.view(), w.br_sea.view(), w.psim_sea.view(), w.psih_sea.view(),
                     w.fm_sea.view(), w.fh_sea.view(), w.znt_sea.view(), w.zol_sea.view(), w.ust_sea.view(),
                     w.cpm_sea.view(), w.chs2_sea.view(), w.cqs2_sea.view(), w.chs_sea.view(), w.flhc_sea.view(),
                     w.flqc_sea.view(), w.qgh_sea.view(), w.ustm_sea.view(), w.lh_sea.view(), w.u10_sea.view(),
                     w.v10_sea.view(), w.th2_sea.view(), w.t2_sea.view(), w.q2_sea.view(), w.ck_sea.view(),
                     w.cka_sea.view(), w.cd_sea.view(), w.cda_sea.view());
      column_for(cols, KSfclayrevSeaBlend{
          xice.cview(), w.br_sea.cview(), w.gz1oz0_sea.cview(), w.mol_sea.cview(), w.psih_sea.cview(),
          w.psim_sea.cview(), w.fh_sea.cview(), w.fm_sea.cview(), w.regime_hold.cview(), w.rmol_sea.cview(),
          w.ust_sea.cview(), w.wspd_sea.cview(), w.zol_sea.cview(), w.cd_sea.cview(), w.cda_sea.cview(),
          w.ck_sea.cview(), w.cka_sea.cview(), w.q2_sea.cview(), w.t2_sea.cview(), w.th2_sea.cview(),
          w.u10_sea.cview(), w.ustm_sea.cview(), w.v10_sea.cview(), br.view(), gz1oz0.view(), mol.view(),
          psih.view(), psim.view(), fh.view(), fm.view(), regime.view(), rmol.view(), ust.view(), wspd.view(),
          zol.view(), cd.view(), cda.view(), ck.view(), cka.view(), q2.view(), t2.view(), th2.view(), u10.view(),
          ustm.view(), v10.view(), xice_threshold});
    }
    column_for(cols, KCopy2{ch.view(), chs.cview()});  // ch = chs (EM_CORE)
  } else if (o_.sf_sfclay_physics == 2) {
    // MYJSFC (module_sf_myjsfc.F: u_phy, v_phy, p_phy=grid%p_hyd, p8w=grid%p_hyd_w)
    const C2 ht_v = b.ht.uniform() ? slice0(zero3_) : C2{b.ht.data(), int(b.ht.layout().sj)};
    auto myjsfc_call = [&](C2 tsk_v, C2 xland_v, C2 z0_v, C2 mavail_v, V2 qsfc_v, V2 thz0_v, V2 qz0_v, V2 uz0_v,
                           V2 vz0_v, V2 ust_v, V2 znt_v, V2 pblh_v, V2 rmol_v, V2 akhs_v, V2 akms_v, V2 br_v, V2 chs_v,
                           V2 chs2_v, V2 cqs2_v, V2 hfx_v, V2 qfx_v, V2 lh_v, V2 flhc_v, V2 flqc_v, V2 qgh_v,
                           V2 cpm_v, V2 ct_v, V2 u10_v, V2 v10_v, V2 t2_v, V2 th2_v, V2 tshltr_v, V2 th10_v, V2 q2_v,
                           V2 qshltr_v, V2 q10_v, V2 pshltr_v, V2 u10e_v, V2 v10e_v) {
      sf_myjsfc::run_myjsfc(
          ht_v, prep.dz8w.cview(), prep.p_hyd.cview(), prep.p_hyd_w.cview(), prep.th_phy.cview(),
          prep.t_phy.cview(), s.q_2[s.P_QV].cview(), s.P_QC >= 0 ? s.q_2[s.P_QC].cview() : zero3_.cview(),
          prep.u_phy.cview(), prep.v_phy.cview(), tke_pbl.cview(), tsk_v, xland_v, ivgtyp.cview(), z0_v, mavail_v,
          lowlyr.cview(), qsfc_v, thz0_v, qz0_v, uz0_v, vz0_v, ust_v, znt_v, pblh_v, rmol_v, akhs_v, akms_v, br_v,
          chs_v, chs2_v, cqs2_v, hfx_v, qfx_v, lh_v, flhc_v, flqc_v, qgh_v, cpm_v, ct_v, u10_v, v10_v, t2_v, th2_v,
          tshltr_v, th10_v, q2_v, qshltr_v, q10_v, pshltr_v, u10e_v, v10e_v, myjsfc_tables_.get_tables(), cols,
          int(itimestep), o_.isurban, o_.iz0tlnd, wrfc::p1000mb, nz);
    };
    if (frac) {  // myjsfc_seaice_wrapper: TSK and QSFC become the ice portion's before the ice call
      SeaIceFields& w = *sea_;
      column_for(cols, KMyjsfcSeaPre{
          xice.cview(), xland.cview(), mavail.cview(), z0.cview(), znt.cview(), qz0.cview(), thz0.cview(),
          uz0.cview(), vz0.cview(), ust.cview(), pblh.cview(), rmol.cview(), akhs.cview(), akms.cview(),
          hfx.cview(), qfx.cview(), lh.cview(), prep.p_hyd_w.cview(), lowlyr.cview(), tsk.view(), qsfc.view(),
          sst.view(), w.tsk_sea.view(), w.xland_sea.view(), w.mavail_sea.view(), w.znt_sea.view(),
          w.z0base_sea.view(), w.qsfc_sea.view(), w.qz0_sea.view(), w.thz0_sea.view(), w.uz0_sea.view(),
          w.vz0_sea.view(), w.ust_sea.view(), w.pblh_sea.view(), w.rmol_sea.view(), w.akhs_sea.view(),
          w.akms_sea.view(), w.hfx_sea.view(), w.qfx_sea.view(), w.lh_sea.view(), xice_threshold, itimestep,
          o_.tice2tsk_if2cold});
    }
    myjsfc_call(tsk.cview(), xland.cview(), z0.cview(), mavail.cview(), qsfc.view(), thz0.view(), qz0.view(),
                uz0.view(), vz0.view(), ust.view(), znt.view(), pblh.view(), rmol.view(), akhs.view(), akms.view(),
                br.view(), chs.view(), chs2.view(), cqs2.view(), hfx.view(), qfx.view(), lh.view(), flhc.view(),
                flqc.view(), qgh.view(), cpm.view(), ct.view(), u10.view(), v10.view(), t2.view(), th2.view(),
                tshltr.view(), th10.view(), q2.view(), qshltr.view(), q10.view(), pshltr.view(), u10e.view(),
                v10e.view());
    if (frac) {
      SeaIceFields& w = *sea_;
      myjsfc_call(w.tsk_sea.cview(), w.xland_sea.cview(), w.z0base_sea.cview(), w.mavail_sea.cview(),
                  w.qsfc_sea.view(), w.thz0_sea.view(), w.qz0_sea.view(), w.uz0_sea.view(), w.vz0_sea.view(),
                  w.ust_sea.view(), w.znt_sea.view(), w.pblh_sea.view(), w.rmol_sea.view(), w.akhs_sea.view(),
                  w.akms_sea.view(), w.br_sea.view(), w.chs_sea.view(), w.chs2_sea.view(), w.cqs2_sea.view(),
                  w.hfx_sea.view(), w.qfx_sea.view(), w.lh_sea.view(), w.flhc_sea.view(), w.flqc_sea.view(),
                  w.qgh_sea.view(), w.cpm_sea.view(), w.ct_sea.view(), w.u10_sea.view(), w.v10_sea.view(),
                  w.t2_sea.view(), w.th2_sea.view(), w.tshltr_sea.view(), w.th10_sea.view(), w.q2_sea.view(),
                  w.qshltr_sea.view(), w.q10_sea.view(), w.pshltr_sea.view(), w.u10e_sea.view(), w.v10e_sea.view());
      column_for(cols, KMyjsfcSeaBlend{
          xice.cview(), w.ct_sea.cview(), w.pshltr_sea.cview(), w.qshltr_sea.cview(), w.q2_sea.cview(),
          w.q10_sea.cview(), w.th2_sea.cview(), w.th10_sea.cview(), w.tshltr_sea.cview(), w.t2_sea.cview(),
          w.u10_sea.cview(), w.v10_sea.cview(), w.thz0_sea.cview(), w.uz0_sea.cview(), w.vz0_sea.cview(),
          w.ust_sea.cview(), w.pblh_sea.cview(), w.rmol_sea.cview(), w.akhs_sea.cview(), w.akms_sea.cview(),
          ct.view(), pshltr.view(), qshltr.view(), q2.view(), q10.view(), th2.view(), th10.view(), tshltr.view(),
          t2.view(), u10.view(), v10.view(), u10e.view(), v10e.view(), thz0.view(), uz0.view(), vz0.view(),
          ust.view(), pblh.view(), rmol.view(), akhs.view(), akms.view(), xice_threshold});
    }
    column_for(cols, KMyjsfcPost{ustm.view(), wspd.view(), ch.view(), ust.cview(), chs.cview(),
                                 prep.u_phy.cview(), prep.v_phy.cview()});
  }

  column_for(cols, KSfcRatios{u_phytmp.cview(), v_phytmp.cview(), prep.t_phy.cview(), prep.p_hyd.cview(), u10.cview(),
                              v10.cview(), th2.cview(), uratx.view(), vratx.view(), tratx.view(), g_.nx});
  if (o_.isfflx != 0) {
    const int nsoil = o_.num_soil_layers;
    column_for(cols, KNoahRa{wspd.cview(), ust.cview(), ra.view()});
    auto run_seaice = [&]() {
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
      si.op.albedo_default = real(o_.seaice_albedo_default);
      si.op.snowdepth_max = real(o_.seaice_snowdepth_max);
      si.op.snowdepth_min = real(o_.seaice_snowdepth_min);
      si.thickness_opt = o_.seaice_thickness_opt; si.thickness_default = real(o_.seaice_thickness_default);
      si.dt = dt; si.xice_threshold = xice_threshold; si.nsoil = nsoil; si.frpcpn = true;
      column_for(cols, si);
      fatal_check();
    };

    auto seaice_lsm_pre = [&]() {
      if (frac) column_for(cols, KSeaiceLsmPre{xice.cview(), albedo.view(), emiss.view(), xice_threshold});
    };
    auto seaice_lsm_post = [&]() {
      if (!frac) return;
      SeaIceFields& w = *sea_;
      column_for(cols, KSeaiceLsmPost{xice.cview(), w.flhc_sea.cview(), w.flqc_sea.cview(), w.cpm_sea.cview(),
                                      w.cqs2_sea.cview(), w.chs2_sea.cview(), w.chs_sea.cview(), w.qsfc_sea.cview(),
                                      w.qgh_sea.cview(), w.qz0_sea.cview(), w.hfx_sea.cview(), w.qfx_sea.cview(),
                                      w.lh_sea.cview(), w.tsk_sea.cview(), albedo.view(), emiss.view(), flhc.view(),
                                      flqc.view(), cpm.view(), cqs2.view(), chs2.view(), chs.view(), qsfc.view(),
                                      qgh.view(), qz0.view(), hfx.view(), qfx.view(), lh.view(), tsk_save.view(),
                                      tsk.view(), xice_threshold});
    };

    if (o_.sf_surface_physics == 2) {
      // ---- lsm (module_sf_noahdrv.F) ----
      seaice_lsm_pre();
      KLsmPre pre{};
      pre.dz8w = prep.dz8w.cview(); pre.qv3d = s.q_2[s.P_QV].cview(); pre.p8w3d = prep.p_hyd_w.cview();
      pre.t3d = prep.t_phy.cview();
      pre.qgh = qgh.cview(); pre.glw = glw.cview(); pre.swdown = swdown.cview(); pre.rainbl = rainbl.cview();
      pre.vegfra = vegfra.cview(); pre.shdmin = shdmin.cview(); pre.shdmax = shdmax.cview(); pre.snoalb = snoalb.cview();
      pre.tmn = tmn.cview(); pre.xice = xice.cview(); pre.chs = chs.cview(); pre.rib = br.cview(); pre.sr = sr.cview();
      pre.qz0 = qz0.cview();
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
      pre.myj = (o_.sf_sfclay_physics == 2); pre.frpcpn = true;  // SR is passed by first_rk_step_part1
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
      nk.rdlai2d = o_.rdlai2d; nk.usemonalb = o_.usemonalb; nk.opt_thcnd = o_.opt_thcnd;
      nk.active = n_active.cview(); nk.alb_out = n_albbrd.view();
      // debug: CPPWRF_NOAH_DUMP=dir, CPPWRF_NOAH_DUMP_STEP=n: SFLX inputs/outputs of the land columns at step n
      // in the stream format of tests/phys/lsm_noah/fortran/driver_noah.f90 (replay through WRF's SFLX)
      const char* ndump = std::getenv("CPPWRF_NOAH_DUMP");
      const char* nstep = std::getenv("CPPWRF_NOAH_DUMP_STEP");
      const bool do_ndump = ndump && nstep && itimestep == std::atol(nstep);
      std::vector<std::pair<int, int>> ncol;
      auto L2 = n_lwdn.layout();
      auto h2 = [&](const Field2& f) { std::vector<real> h; f.download(h); return h; };
      auto hi = [&](const IntField2& f) { std::vector<int> h; f.download(h); return h; };
      auto h3 = [&](const Field3& f) { std::vector<real> h; f.download(h); return h; };
      auto write_cols = [&](std::FILE* fp, bool out) {
        auto wr2 = [&](const Field2& f) { auto h = h2(f); for (auto [i, j] : ncol) { double v = h[L2.index(i, j, 0)]; std::fwrite(&v, 8, 1, fp); } };
        auto wr3 = [&](const Field3& f) { auto h = h3(f); const Layout& L = f.layout();
          for (int k = 0; k < nsoil; ++k) for (auto [i, j] : ncol) { double v = h[L.index(i, j, k)]; std::fwrite(&v, 8, 1, fp); } };
        if (!out) {
          for (const IntField2* f : {&n_vegtyp, &n_soiltyp, &n_slopetyp}) { auto h = hi(*f); for (auto [i, j] : ncol) { int v = h[L2.index(i, j, 0)]; std::fwrite(&v, 4, 1, fp); } }
          for (const Field2* f : {&n_shdfac, &n_shmin, &n_shmax, &n_albbrd, &n_snoalb1, &n_tbot, &n_z0brd, &n_embrd, &n_lwdn, &n_soldn,
                                  &n_solnet, &n_sfcprs, &n_prcp, &n_sfctmp, &n_th2, &n_q2, &n_q2, &n_q2sat, &n_dqsdt2, &n_ffrozp,
                                  &n_cmc, &n_t1, &n_snowhk, &n_sneqv, &n_albedok, &n_chk, &n_cm, &n_emissi, &n_z0k, &n_xlai,
                                  &n_snotime1, &n_ribb}) wr2(*f);
          wr3(n_stc); wr3(n_smc); wr3(n_swc);
        } else {
          for (const Field2* f : {&n_eta, &n_sheat, &n_eta_kin, &n_fdown, &n_ec, &n_edir, &n_ett, &n_esnow, &n_drip, &n_dew,
                                  &n_beta, &n_etp, &n_ssoil, &n_flx1, &n_flx2, &n_flx3, &n_snomlt, &n_sncovr, &n_runoff1,
                                  &n_runoff2, &n_runoff3, &n_rc, &n_pc, &n_rsmin, &n_soilw, &n_soilm, &n_q1, &n_cmc, &n_t1,
                                  &n_snowhk, &n_sneqv, &n_albedok, &n_emissi, &n_z0k, &n_xlai, &n_snotime1, &n_ribb}) wr2(*f);
          wr3(n_stc); wr3(n_smc); wr3(n_swc); wr3(n_et); wr3(n_smav);
        }
      };
      if (do_ndump) {
        auto act = hi(n_active);
        for (int j = cols.j0; j < cols.j1; ++j)
          for (int i = cols.i0; i < cols.i1; ++i)
            if (act[L2.index(i, j, 0)]) ncol.emplace_back(i, j);
        std::FILE* fp = std::fopen((std::string(ndump) + "/noah_in.bin").c_str(), "wb");
        const int hdr[2] = {int(ncol.size()), nsoil};
        std::fwrite(hdr, 4, 2, fp);
        const double dd[2] = {double(dtbl), 0.0};
        std::fwrite(dd, 8, 2, fp);
        std::fwrite(&o_.isurban, 4, 1, fp);
        for (int k = 0; k < nsoil; ++k) { double v = dzs[k]; std::fwrite(&v, 8, 1, fp); }
        write_cols(fp, false);
        std::fclose(fp);
        fp = std::fopen((std::string(ndump) + "/noah_cols.txt").c_str(), "w");
        for (auto [i, j] : ncol) std::fprintf(fp, "%d %d\n", i, j);
        std::fclose(fp);
      }
      noah::run_noah(nk, cols);
      column_for(Range2D{0, 1, 0, 1}, KNoahRs{n_rc.cview(), n_active.cview(), n_ice.cview(), rs.view(), cols.i0, cols.i1,
                                             cols.j0, cols.j1});
      if (do_ndump) {
        std::FILE* fp = std::fopen((std::string(ndump) + "/noah_out_c.bin").c_str(), "wb");
        const int hdr[2] = {int(ncol.size()), nsoil};
        std::fwrite(hdr, 4, 2, fp);
        write_cols(fp, true);
        std::fclose(fp);
      }

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
      fatal_check();

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

      run_seaice();
      seaice_lsm_post();

      column_for(cols, KLsmDiag{qfx.cview(), chs.cview(), hfx.cview(), lh.cview(), grdflx.cview(), tsk.cview(),
                                qsfc.cview(), chs2.cview(), cqs2.cview(), psfc.cview(), prep.t_phy.cview(),
                                s.q_2[s.P_QV].cview(), sfcevp.view(), sfcexc.view(), achfx.view(), aclhf.view(),
                                acgrdflx.view(), q2.view(), t2.view(), th2.view(), dtbl, dt});
    } else if (o_.sf_surface_physics == 4) {
      // ---- Noah-MP (NoahmpWRFmain -> NoahmpDriverMain, module_surface_driver.F:3134-3481) ----
      const double rsec = double(o_.rad.julday - 1) * 86400.0 + o_.rad.gmt * 3600.0 + double(itimestep - 1) * g_.dt;
      const real julian_in = real(rsec / 86400.0 + 1.0) - real(1.0);
      noahmp::NoahmpKernel nmp{};
      nmp.ivgtyp = ivgtyp.cview();
      nmp.isltyp = isltyp.cview();
      nmp.cropcat = cropcat.cview();
      nmp.coszen = coszen.cview();
      nmp.xlat = xlat_sfc.cview();
      nmp.xlong = xlong_sfc.cview();
      nmp.dz8w = prep.dz8w.cview();
      nmp.vegfra = vegfra.cview();
      nmp.shdmax = shdmax.cview();
      nmp.tmn = tmn.cview();
      nmp.xland = xland.cview();
      nmp.xice = xice.cview();
      nmp.t_phy = prep.t_phy.cview();
      nmp.qv_curr = s.q_2[s.P_QV].cview();
      nmp.u_phy = prep.u_phy.cview();
      nmp.v_phy = prep.v_phy.cview();
      nmp.swdown = swdown.cview();
      nmp.glw = glw.cview();
      nmp.p8w = prep.p_hyd_w.cview();
      nmp.rainbl = rainbl.cview();
      nmp.sr = sr.cview();
      nmp.raincv = raincv.cview();
      nmp.rainncv = rainncv.cview();
      nmp.rainshv = rainshv.cview();
      nmp.snowncv = snowncv_mp_ ? snowncv_mp_->cview() : C2{};
      nmp.graupelncv = graupelncv_mp_ ? graupelncv_mp_->cview() : C2{};
      nmp.tsk = tsk.view();
      nmp.hfx = hfx.view();
      nmp.qfx = qfx.view();
      nmp.lh = lh.view();
      nmp.grdflx = grdflx.view();
      nmp.smstav = smstav.view();
      nmp.smstot = smstot.view();
      nmp.sfcrunoff = sfcrunoff.view();
      nmp.udrunoff = udrunoff.view();
      nmp.albedo = albedo.view();
      nmp.snowc = snowc.view();
      nmp.snow = snow.view();
      nmp.snowh = snowh.view();
      nmp.canwat = canwat.view();
      nmp.acsnom = acsnom.view();
      nmp.acsnow = acsnow.view();
      nmp.emiss = emiss.view();
      nmp.qsfc = qsfc.view();
      nmp.z0 = z0.view();
      nmp.znt = znt.view();
      nmp.isnowxy = isnowxy.view();
      nmp.tvxy = tvxy.view();
      nmp.tgxy = tgxy.view();
      nmp.canicexy = canicexy.view();
      nmp.canliqxy = canliqxy.view();
      nmp.eahxy = eahxy.view();
      nmp.tahxy = tahxy.view();
      nmp.cmxy = cmxy.view();
      nmp.chxy = chxy.view();
      nmp.fwetxy = fwetxy.view();
      nmp.sneqvoxy = sneqvoxy.view();
      nmp.alboldxy = alboldxy.view();
      nmp.qsnowxy = qsnowxy.view();
      nmp.qrainxy = qrainxy.view();
      nmp.wslakexy = wslakexy.view();
      nmp.zwtxy = zwtxy.view();
      nmp.waxy = waxy.view();
      nmp.wtxy = wtxy.view();
      nmp.lfmassxy = lfmassxy.view();
      nmp.rtmassxy = rtmassxy.view();
      nmp.stmassxy = stmassxy.view();
      nmp.woodxy = woodxy.view();
      nmp.stblcpxy = stblcpxy.view();
      nmp.fastcpxy = fastcpxy.view();
      nmp.lai = lai.view();
      nmp.xsaixy = xsaixy.view();
      nmp.taussxy = taussxy.view();
      nmp.smcwtdxy = smcwtdxy.view();
      nmp.deeprechxy = deeprechxy.view();
      nmp.rechxy = rechxy.view();
      nmp.grainxy = grainxy.view();
      nmp.gddxy = gddxy.view();
      nmp.pgsxy = pgsxy.view();
      nmp.qtdrain = qtdrain.view();
      nmp.smois = smois.view();
      nmp.sh2o = sh2o.view();
      nmp.tslb = tslb.view();
      nmp.smoiseq = smoiseq.view();
      nmp.tsnoxy = tsnoxy.view();
      nmp.zsnsoxy = zsnsoxy.view();
      nmp.snicexy = snicexy.view();
      nmp.snliqxy = snliqxy.view();
      nmp.albsoildirxy = albsoildirxy.view();
      nmp.albsoildifxy = albsoildifxy.view();
      nmp.acc_etranixy = acc_etranixy.view();
      nmp.t2mvxy = t2mvxy.view();
      nmp.t2mbxy = t2mbxy.view();
      nmp.q2mvxy = q2mvxy.view();
      nmp.q2mbxy = q2mbxy.view();
      nmp.tradxy = tradxy.view();
      nmp.neexy = neexy.view();
      nmp.gppxy = gppxy.view();
      nmp.nppxy = nppxy.view();
      nmp.fvegxy = fvegxy.view();
      nmp.runsfxy = runsfxy.view();
      nmp.runsbxy = runsbxy.view();
      nmp.ecanxy = ecanxy.view();
      nmp.edirxy = edirxy.view();
      nmp.etranxy = etranxy.view();
      nmp.fsaxy = fsaxy.view();
      nmp.firaxy = firaxy.view();
      nmp.aparxy = aparxy.view();
      nmp.psnxy = psnxy.view();
      nmp.savxy = savxy.view();
      nmp.sagxy = sagxy.view();
      nmp.rssunxy = rssunxy.view();
      nmp.rsshaxy = rsshaxy.view();
      nmp.bgapxy = bgapxy.view();
      nmp.wgapxy = wgapxy.view();
      nmp.tgvxy = tgvxy.view();
      nmp.tgbxy = tgbxy.view();
      nmp.chvxy = chvxy.view();
      nmp.chbxy = chbxy.view();
      nmp.shgxy = shgxy.view();
      nmp.shcxy = shcxy.view();
      nmp.shbxy = shbxy.view();
      nmp.evgxy = evgxy.view();
      nmp.evbxy = evbxy.view();
      nmp.ghvxy = ghvxy.view();
      nmp.ghbxy = ghbxy.view();
      nmp.irgxy = irgxy.view();
      nmp.ircxy = ircxy.view();
      nmp.irbxy = irbxy.view();
      nmp.trxy = trxy.view();
      nmp.evcxy = evcxy.view();
      nmp.chleafxy = chleafxy.view();
      nmp.chucxy = chucxy.view();
      nmp.chv2xy = chv2xy.view();
      nmp.chb2xy = chb2xy.view();
      nmp.rs = rs.view();
      nmp.acc_ssoilxy = acc_ssoilxy.view();
      nmp.acc_qinsurxy = acc_qinsurxy.view();
      nmp.acc_qsevaxy = acc_qsevaxy.view();
      nmp.eflxbxy = eflxbxy.view();
      nmp.acc_dwaterxy = acc_dwaterxy.view();
      nmp.acc_prcpxy = acc_prcpxy.view();
      nmp.acc_ecanxy = acc_ecanxy.view();
      nmp.acc_etranxy = acc_etranxy.view();
      nmp.acc_edirxy = acc_edirxy.view();
      nmp.acc_glaflwxy = acc_glaflwxy.view();
      nmp.tables = d_noahmp_tables_;
      nmp.opts = o_.noahmp;
      nmp.itimestep = int(itimestep);
      nmp.yr = o_.rad.julyr > 0 ? o_.rad.julyr : 2000;
      nmp.julian_in = julian_in;
      nmp.dtbl = dtbl;
      nmp.dx = real(g_.dx);
      nmp.xice_threshold = xice_threshold;
      nmp.nsoil = nsoil;
      nmp.nsnow = noahmp::kMaxSnowLayers;
      for (int k = 0; k < nsoil; ++k) nmp.dzs[k] = real(dzs[k]);
      seaice_lsm_pre();
      noahmp::run_noahmp(nmp, cols);

      run_seaice();
      seaice_lsm_post();

      KNoahmpDiag nd{};
      nd.ivgtyp = ivgtyp.cview();
      nd.xice = xice.cview();
      nd.qfx = qfx.cview();
      nd.chs = chs.cview();
      nd.hfx = hfx.cview();
      nd.lh = lh.cview();
      nd.grdflx = grdflx.cview();
      nd.tsk = tsk.cview();
      nd.qsfc = qsfc.cview();
      nd.chs2 = chs2.cview();
      nd.cqs2 = cqs2.cview();
      nd.psfc = psfc.cview();
      nd.fvegxy = fvegxy.cview();
      nd.t2mvxy = t2mvxy.cview();
      nd.t2mbxy = t2mbxy.cview();
      nd.q2mvxy = q2mvxy.cview();
      nd.q2mbxy = q2mbxy.cview();
      nd.chklowq = chklowq.view();
      nd.sfcevp = sfcevp.view();
      nd.sfcexc = sfcexc.view();
      nd.achfx = achfx.view();
      nd.aclhf = aclhf.view();
      nd.acgrdflx = acgrdflx.view();
      nd.q2 = q2.view();
      nd.t2 = t2.view();
      nd.th2 = th2.view();
      nd.dtbl = dtbl;
      nd.dt = dt;
      nd.xice_threshold = xice_threshold;
      nd.iswater = o_.iswater;
      nd.isice = o_.isice;
      nd.isurban = o_.isurban;
      const int lcz_arr[11] = {noahmp_tables_.lcz_1, noahmp_tables_.lcz_2, noahmp_tables_.lcz_3, noahmp_tables_.lcz_4,
                               noahmp_tables_.lcz_5, noahmp_tables_.lcz_6, noahmp_tables_.lcz_7, noahmp_tables_.lcz_8,
                               noahmp_tables_.lcz_9, noahmp_tables_.lcz_10, noahmp_tables_.lcz_11};
      for (int n = 0; n < 11; ++n) nd.lcz[n] = lcz_arr[n];
      column_for(cols, nd);
    }
  }
  column_for(cols, KSfcEnd{rainbl.view(), q2.view(), xland.cview(), s.q_2[s.P_QV].cview()});
}

}  // namespace cppwrf::phys::driver
