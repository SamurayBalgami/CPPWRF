#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "phys/lsm_noahmp/noahmp_energy.hpp"
#include "phys/lsm_noahmp/noahmp_glacier.hpp"
#include "phys/lsm_noahmp/noahmp_rad_therm.hpp"
#include "phys/lsm_noahmp/noahmp_tables.hpp"
#include "phys/lsm_noahmp/noahmp_var.hpp"
#include "phys/lsm_noahmp/noahmp_water.hpp"

namespace cppwrf::phys::noahmp {

// ============================================================================
// Calendar utility (from module_ra_gfdleta: cal_mon_day)
// ============================================================================
CPPWRF_HD inline void cal_mon_day(int julday, int yr, int& month, int& day) {
  int day_of_mon[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  if (yr % 4 == 0) day_of_mon[1] = 29;
  for (int m = 1; m <= 12; ++m) {
    if (julday <= day_of_mon[m - 1]) {
      month = m;
      day = julday;
      break;
    } else {
      julday = julday - day_of_mon[m - 1];
    }
  }
}

// ============================================================================
// Initialization Kernel (NoahmpWRFinit -> NoahmpInitMain + NoahmpSnowInitMain)
// ============================================================================
struct NoahmpInitKernel {
  View2<const int> ivgtyp;
  View2<const int> isltyp;
  View2<const real> tsk;
  View2<const real> xice;
  View2<real> tmn;
  View2<real> snow;
  View2<real> snowh;
  View2<real> canwat;
  View3<real> tslb;
  View3<real> smois;
  View3<real> sh2o;
  View2<int> isnowxy;
  View2<real> tvxy;
  View2<real> tgxy;
  View2<real> canicexy;
  View2<real> canliqxy;
  View2<real> eahxy;
  View2<real> tahxy;
  View2<real> cmxy;
  View2<real> chxy;
  View2<real> fwetxy;
  View2<real> sneqvoxy;
  View2<real> alboldxy;
  View2<real> qsnowxy;
  View2<real> qrainxy;
  View2<real> wslakexy;
  View2<real> zwtxy;
  View2<real> waxy;
  View2<real> wtxy;
  View2<real> lfmassxy;
  View2<real> rtmassxy;
  View2<real> stmassxy;
  View2<real> woodxy;
  View2<real> stblcpxy;
  View2<real> fastcpxy;
  View2<real> xsaixy;
  View2<real> lai;
  View2<real> grainxy;
  View2<real> gddxy;
  View2<int> cropcat;
  View2<real> t2mvxy;
  View2<real> t2mbxy;
  View2<real> chstarxy;
  View2<real> qtdrain;
  View3<real> tsnoxy;
  View3<real> snicexy;
  View3<real> snliqxy;
  View3<real> zsnsoxy;
  View3<real> albsoildirxy;
  View3<real> albsoildifxy;

  const NoahmpTables* tables = nullptr;
  int nsoil = 4;
  int nsnow = 3;
  real dzs[kMaxSoilLayers] = {0};
  bool fndsnowh = true;
  bool restart = false;
  int iopt_runsub = 3;
  int iopt_crop = 0;
  int sf_urban_physics = 0;

  CPPWRF_HD void operator()(int i, int j) const {
    if (restart) {
      if (chstarxy.p) chstarxy(i, j) = real(0.1);
      return;
    }

    const NoahmpTables& tb = *tables;
    const int st = isltyp(i, j);
    const int vt = ivgtyp(i, j);

    // 1. Initialize physical snow height SNOWH and check consistency (NoahmpInitMainMod lines 48-73)
    if (!fndsnowh) {
      snowh(i, j) = snow(i, j) * real(0.005);
    }
    if (snow(i, j) < real(0.0)) snow(i, j) = real(0.0);
    if (snowh(i, j) < real(0.0)) snowh(i, j) = real(0.0);
    if ((snow(i, j) > real(0.0)) && (snowh(i, j) == real(0.0))) {
      snowh(i, j) = snow(i, j) * real(0.005);
    }
    if ((snowh(i, j) > real(0.0)) && (snow(i, j) == real(0.0))) {
      snow(i, j) = snowh(i, j) / real(0.005);
    }
    if (snow(i, j) > real(2000.0)) {
      snowh(i, j) = snowh(i, j) * real(2000.0) / snow(i, j);
      snow(i, j) = real(2000.0);
    }

    // 2. Soil layer depths ZSOIL
    FArray1D<real, 1, kMaxSoilLayers> zsoil;
    zsoil(1) = -dzs[0];
    for (int ns = 2; ns <= nsoil; ++ns) {
      zsoil(ns) = zsoil(ns - 1) - dzs[ns - 1];
    }

    // 3. Initialize soil liquid water content SH2O (NoahmpInitMainMod lines 95-132)
    constexpr real hlice = real(3.335e5);
    constexpr real grav0 = real(9.81);
    constexpr real t0 = real(273.15);
    if ((vt == tb.isice) && (xice(i, j) <= real(0.0))) {
      for (int ns = 0; ns < nsoil; ++ns) {
        smois(i, j, ns) = real(1.0);
        sh2o(i, j, ns) = real(0.0);
        tslb(i, j, ns) = std::min(tslb(i, j, ns), real(263.15));
      }
      snow(i, j) = std::max(snow(i, j), real(10.0));
      snowh(i, j) = snow(i, j) * real(0.005);
    } else {
      real bexp = tb.bexp[st];
      real smcmax = tb.smcmax[st];
      real psisat = tb.psisat[st];
      for (int ns = 0; ns < nsoil; ++ns) {
        if (smois(i, j, ns) > smcmax) smois(i, j, ns) = smcmax;
      }
      if ((bexp > real(0.0)) && (smcmax > real(0.0)) && (psisat > real(0.0))) {
        for (int ns = 0; ns < nsoil; ++ns) {
          if (tslb(i, j, ns) < real(273.149)) {
            real fk = rpow(((hlice / (grav0 * (-psisat))) *
                            ((tslb(i, j, ns) - t0) / tslb(i, j, ns))),
                           -(real(1.0) / bexp)) *
                      smcmax;
            fk = std::max(fk, real(0.02));
            sh2o(i, j, ns) = std::min(fk, smois(i, j, ns));
          } else {
            sh2o(i, j, ns) = smois(i, j, ns);
          }
        }
      } else {
        for (int ns = 0; ns < nsoil; ++ns) {
          sh2o(i, j, ns) = smois(i, j, ns);
        }
      }
    }

    // 4. Initialize other quantities (NoahmpInitMainMod lines 135-266)
    if (qtdrain.p) qtdrain(i, j) = real(0.0);
    tvxy(i, j) = tsk(i, j);
    tgxy(i, j) = tsk(i, j);
    if ((snow(i, j) > real(0.0)) && (tsk(i, j) > real(273.15))) tvxy(i, j) = real(273.15);
    if ((snow(i, j) > real(0.0)) && (tsk(i, j) > real(273.15))) tgxy(i, j) = real(273.15);
    canwat(i, j) = real(0.0);
    canliqxy(i, j) = canwat(i, j);
    canicexy(i, j) = real(0.0);
    eahxy(i, j) = real(2000.0);
    tahxy(i, j) = tsk(i, j);
    t2mvxy(i, j) = tsk(i, j);
    t2mbxy(i, j) = tsk(i, j);
    if ((snow(i, j) > real(0.0)) && (tsk(i, j) > real(273.15))) tahxy(i, j) = real(273.15);
    if ((snow(i, j) > real(0.0)) && (tsk(i, j) > real(273.15))) t2mvxy(i, j) = real(273.15);
    if ((snow(i, j) > real(0.0)) && (tsk(i, j) > real(273.15))) t2mbxy(i, j) = real(273.15);
    cmxy(i, j) = real(0.0);
    chxy(i, j) = real(0.0);
    fwetxy(i, j) = real(0.0);
    sneqvoxy(i, j) = real(0.0);
    alboldxy(i, j) = real(0.65);
    qsnowxy(i, j) = real(0.0);
    qrainxy(i, j) = real(0.0);
    wslakexy(i, j) = real(0.0);

    if (iopt_runsub != 5) {
      waxy(i, j) = real(4900.0);
      wtxy(i, j) = waxy(i, j);
      zwtxy(i, j) = (real(25.0) - zsoil(nsoil)) - waxy(i, j) / real(1000.0) / real(0.2);
    } else {
      waxy(i, j) = real(0.0);
      wtxy(i, j) = real(0.0);
    }

    bool urbanpt_flag = (vt == tb.isurban) || (vt > tb.urbtype_beg);
    if ((vt == tb.isbarren) || (vt == tb.isice) ||
        ((sf_urban_physics == 0) && urbanpt_flag) ||
        (vt == tb.iswater)) {
      lai(i, j) = real(0.0);
      xsaixy(i, j) = real(0.0);
      lfmassxy(i, j) = real(0.0);
      stmassxy(i, j) = real(0.0);
      rtmassxy(i, j) = real(0.0);
      woodxy(i, j) = real(0.0);
      stblcpxy(i, j) = real(0.0);
      fastcpxy(i, j) = real(0.0);
      if (grainxy.p) grainxy(i, j) = real(1.0e-10);
      if (gddxy.p) gddxy(i, j) = real(0.0);
      if (cropcat.p) cropcat(i, j) = 0;
    } else {
      if ((lai(i, j) > real(100.0)) || (lai(i, j) < real(0.0))) {
        lai(i, j) = real(0.0);
      }
      lai(i, j) = std::max(lai(i, j), real(0.05));
      xsaixy(i, j) = std::max(real(0.1) * lai(i, j), real(0.05));
      if (urbanpt_flag) {
        lfmassxy(i, j) = lai(i, j) * real(1000.0) / std::max(tb.sla[tb.natural], real(1.0));
      } else {
        lfmassxy(i, j) = lai(i, j) * real(1000.0) / std::max(tb.sla[vt], real(1.0));
      }
      stmassxy(i, j) = xsaixy(i, j) * real(1000.0) / real(3.0);
      rtmassxy(i, j) = real(500.0);
      woodxy(i, j) = real(500.0);
      stblcpxy(i, j) = real(1000.0);
      fastcpxy(i, j) = real(1000.0);
      if (grainxy.p) grainxy(i, j) = real(1.0e-10);
      if (gddxy.p) gddxy(i, j) = real(0.0);
      if (cropcat.p) {
        if (iopt_crop == 1) {
          cropcat(i, j) = tb.default_crop;
        } else {
          cropcat(i, j) = undefined_int;
        }
      }
    }

    if (albsoildirxy.p && albsoildifxy.p) {
      for (int ib = 0; ib < kNumRadBands; ++ib) {
        albsoildirxy(i, j, ib) = real(0.0);
        albsoildifxy(i, j, ib) = real(0.0);
      }
    }
    if (chstarxy.p) chstarxy(i, j) = real(0.1);

    // 5. NoahmpSnowInitMain
    FArray1D<real, -kMaxSnowLayers + 1, 0> dzsno;
    FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> dzsnso;
    FArray1D<real, -kMaxSnowLayers + 1, 0> tsno;
    FArray1D<real, -kMaxSnowLayers + 1, 0> snice;
    FArray1D<real, -kMaxSnowLayers + 1, 0> snliq;
    FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> zsnso;
    dzsno.fill(real(0.0));
    dzsnso.fill(real(0.0));
    tsno.fill(real(0.0));
    snice.fill(real(0.0));
    snliq.fill(real(0.0));
    zsnso.fill(real(0.0));

    int isnow = 0;
    real sh = snowh(i, j);
    if (sh < real(0.025)) {
      isnow = 0;
      dzsno.fill(real(0.0));
    } else {
      if ((sh >= real(0.025)) && (sh <= real(0.05))) {
        isnow = -1;
        dzsno(0) = sh;
      } else if ((sh > real(0.05)) && (sh <= real(0.10))) {
        isnow = -2;
        dzsno(-1) = sh / real(2.0);
        dzsno(0) = sh / real(2.0);
      } else if ((sh > real(0.10)) && (sh <= real(0.25))) {
        isnow = -2;
        dzsno(-1) = real(0.05);
        dzsno(0) = sh - dzsno(-1);
      } else if ((sh > real(0.25)) && (sh <= real(0.45))) {
        isnow = -3;
        dzsno(-2) = real(0.05);
        dzsno(-1) = real(0.5) * (sh - dzsno(-2));
        dzsno(0) = real(0.5) * (sh - dzsno(-2));
      } else if (sh > real(0.45)) {
        isnow = -3;
        dzsno(-2) = real(0.05);
        dzsno(-1) = real(0.20);
        dzsno(0) = sh - dzsno(-1) - dzsno(-2);
      }
    }

    for (int iz = isnow + 1; iz <= 0; ++iz) {
      tsno(iz) = tgxy(i, j);
      snliq(iz) = real(0.0);
      snice(iz) = real(1.0) * dzsno(iz) * (snow(i, j) / sh);
      dzsnso(iz) = -dzsno(iz);
    }
    dzsnso(1) = zsoil(1);
    for (int iz = 2; iz <= nsoil; ++iz) {
      dzsnso(iz) = zsoil(iz) - zsoil(iz - 1);
    }
    zsnso(isnow + 1) = dzsnso(isnow + 1);
    for (int iz = isnow + 2; iz <= nsoil; ++iz) {
      zsnso(iz) = zsnso(iz - 1) + dzsnso(iz);
    }

    isnowxy(i, j) = isnow;
    for (int iz = -nsnow + 1; iz <= 0; ++iz) {
      tsnoxy(i, j, iz + nsnow - 1) = tsno(iz);
      snicexy(i, j, iz + nsnow - 1) = snice(iz);
      snliqxy(i, j, iz + nsnow - 1) = snliq(iz);
    }
    for (int iz = -nsnow + 1; iz <= nsoil; ++iz) {
      zsnsoxy(i, j, iz + nsnow - 1) = zsnso(iz);
    }
  }
};

// ============================================================================
// Per-Step Column Kernel (NoahmpWRFmain -> NoahmpDriverMain)
// ============================================================================
struct NoahmpKernel {
  // 2D/3D const inputs
  View2<const int> ivgtyp;
  View2<const int> isltyp;
  View2<const int> cropcat;
  View2<const real> coszen;
  View2<const real> xlat;
  View2<const real> xlong;
  View3<const real> dz8w;
  View2<const real> dz8w_0;
  View2<const real> vegfra;
  View2<const real> shdmax;
  View2<const real> tmn;
  View2<const real> xland;
  View2<const real> xice;
  View3<const real> t_phy;
  View2<const real> t_phy_0;
  View3<const real> qv_curr;
  View2<const real> qv_curr_0;
  View3<const real> u_phy;
  View2<const real> u_phy_0;
  View3<const real> v_phy;
  View2<const real> v_phy_0;
  View2<const real> swdown;
  View2<const real> swddir;
  View2<const real> swddif;
  View2<const real> glw;
  View3<const real> p8w;
  View2<const real> p8w_0;
  View2<const real> p8w_1;
  View2<const real> rainbl;
  View2<const real> sr;
  View2<const real> raincv;
  View2<const real> rainncv;
  View2<const real> rainshv;
  View2<const real> snowncv;
  View2<const real> graupelncv;
  View2<const real> hailncv;

  // 2D INOUT state & fluxes
  View2<real> tsk;
  View2<real> hfx;
  View2<real> qfx;
  View2<real> lh;
  View2<real> grdflx;
  View2<real> smstav;
  View2<real> smstot;
  View2<real> sfcrunoff;
  View2<real> udrunoff;
  View2<real> albedo;
  View2<real> snowc;
  View2<real> snow;
  View2<real> snowh;
  View2<real> canwat;
  View2<real> acsnom;
  View2<real> acsnow;
  View2<real> emiss;
  View2<real> qsfc;
  View2<real> z0;
  View2<real> znt;

  // 2D INOUT Noah-MP state
  View2<int> isnowxy;
  View2<real> tvxy;
  View2<real> tgxy;
  View2<real> canicexy;
  View2<real> canliqxy;
  View2<real> eahxy;
  View2<real> tahxy;
  View2<real> cmxy;
  View2<real> chxy;
  View2<real> fwetxy;
  View2<real> sneqvoxy;
  View2<real> alboldxy;
  View2<real> qsnowxy;
  View2<real> qrainxy;
  View2<real> wslakexy;
  View2<real> zwtxy;
  View2<real> waxy;
  View2<real> wtxy;
  View2<real> lfmassxy;
  View2<real> rtmassxy;
  View2<real> stmassxy;
  View2<real> woodxy;
  View2<real> stblcpxy;
  View2<real> fastcpxy;
  View2<real> lai;
  View2<real> xsaixy;
  View2<real> taussxy;
  View2<real> smcwtdxy;
  View2<real> deeprechxy;
  View2<real> rechxy;
  View2<real> grainxy;
  View2<real> gddxy;
  View2<int> pgsxy;
  View2<real> qtdrain;

  // 3D INOUT state
  View3<real> smois;
  View3<real> sh2o;
  View3<real> tslb;
  View3<real> smoiseq;
  View3<real> tsnoxy;
  View3<real> zsnsoxy;
  View3<real> snicexy;
  View3<real> snliqxy;
  View3<real> albsoildirxy;
  View3<real> albsoildifxy;
  View3<real> acc_etranixy;

  // 2D OUT / INOUT diagnostics
  View2<real> t2mvxy;
  View2<real> t2mbxy;
  View2<real> q2mvxy;
  View2<real> q2mbxy;
  View2<real> tradxy;
  View2<real> neexy;
  View2<real> gppxy;
  View2<real> nppxy;
  View2<real> fvegxy;
  View2<real> runsfxy;
  View2<real> runsbxy;
  View2<real> ecanxy;
  View2<real> edirxy;
  View2<real> etranxy;
  View2<real> fsaxy;
  View2<real> firaxy;
  View2<real> aparxy;
  View2<real> psnxy;
  View2<real> savxy;
  View2<real> sagxy;
  View2<real> rssunxy;
  View2<real> rsshaxy;
  View2<real> bgapxy;
  View2<real> wgapxy;
  View2<real> tgvxy;
  View2<real> tgbxy;
  View2<real> chvxy;
  View2<real> chbxy;
  View2<real> shgxy;
  View2<real> shcxy;
  View2<real> shbxy;
  View2<real> evgxy;
  View2<real> evbxy;
  View2<real> ghvxy;
  View2<real> ghbxy;
  View2<real> irgxy;
  View2<real> ircxy;
  View2<real> irbxy;
  View2<real> trxy;
  View2<real> evcxy;
  View2<real> chleafxy;
  View2<real> chucxy;
  View2<real> chv2xy;
  View2<real> chb2xy;
  View2<real> rs;
  View2<real> qintsxy;
  View2<real> qintrxy;
  View2<real> qdripsxy;
  View2<real> qdriprxy;
  View2<real> qthrosxy;
  View2<real> qthrorxy;
  View2<real> qsnsubxy;
  View2<real> qsnfroxy;
  View2<real> qsubcxy;
  View2<real> qfrocxy;
  View2<real> qevacxy;
  View2<real> qdewcxy;
  View2<real> qfrzcxy;
  View2<real> qmeltcxy;
  View2<real> qsnbotxy;
  View2<real> qmeltxy;
  View2<real> pondingxy;
  View2<real> pahxy;
  View2<real> pahgxy;
  View2<real> pahvxy;
  View2<real> pahbxy;
  View2<real> fpicexy;
  View2<real> rainlsm;
  View2<real> snowlsm;
  View2<real> forctlsm;
  View2<real> forcqlsm;
  View2<real> forcplsm;
  View2<real> forczlsm;
  View2<real> forcwlsm;
  View2<real> acc_ssoilxy;
  View2<real> acc_qinsurxy;
  View2<real> acc_qsevaxy;
  View2<real> eflxbxy;
  View2<real> soilenergy;
  View2<real> snowenergy;
  View2<real> canhsxy;
  View2<real> acc_dwaterxy;
  View2<real> acc_prcpxy;
  View2<real> acc_ecanxy;
  View2<real> acc_etranxy;
  View2<real> acc_edirxy;
  View2<real> acc_glaflwxy;

  // Configuration
  const NoahmpTables* tables = nullptr;
  NoahmpOptions opts;
  int itimestep = 1;
  int yr = 2000;
  real julian_in = 0.0;
  real dtbl = 0.0;
  real dx = 0.0;
  real xice_threshold = 0.5;
  int nsoil = 4;
  int nsnow = 3;
  real dzs[kMaxSoilLayers] = {0};

  CPPWRF_HD void operator()(int i, int j) const {
    const NoahmpTables& tb = *tables;

    int soil_update_steps = nint_f(real(opts.soiltstep) / dtbl);
    soil_update_steps = std::max(soil_update_steps, 1);
    bool calculate_soil = (itimestep % soil_update_steps == 0);

    if ((soil_update_steps == 1) || ((soil_update_steps > 1) && (itimestep % soil_update_steps == 1))) {
      if (acc_ssoilxy.p) acc_ssoilxy(i, j) = real(0.0);
      if (acc_qinsurxy.p) acc_qinsurxy(i, j) = real(0.0);
      if (acc_qsevaxy.p) acc_qsevaxy(i, j) = real(0.0);
      if (acc_etranixy.p) {
        for (int k = 0; k < nsoil; ++k) acc_etranixy(i, j, k) = real(0.0);
      }
      if (acc_dwaterxy.p) acc_dwaterxy(i, j) = real(0.0);
      if (acc_prcpxy.p) acc_prcpxy(i, j) = real(0.0);
      if (acc_ecanxy.p) acc_ecanxy(i, j) = real(0.0);
      if (acc_etranxy.p) acc_etranxy(i, j) = real(0.0);
      if (acc_edirxy.p) acc_edirxy(i, j) = real(0.0);
      if (acc_glaflwxy.p) acc_glaflwxy(i, j) = real(0.0);
    }

    int yearlen = 365;
    if (yr % 4 == 0) {
      yearlen = 366;
      if (yr % 100 == 0) {
        yearlen = 365;
        if (yr % 400 == 0) {
          yearlen = 366;
        }
      }
    }

    FArray1D<real, 1, kMaxSoilLayers> zsoil;
    zsoil(1) = -dzs[0];
    for (int k = 2; k <= nsoil; ++k) {
      zsoil(k) = -dzs[k - 1] + zsoil(k - 1);
    }

    // Initialization over ocean / sea-ice at itimestep == 1 (NoahmpDriverMainMod lines 125-145)
    if (itimestep == 1) {
      if ((xland(i, j) - real(1.5)) >= real(0.0)) {
        smstav(i, j) = real(1.0);
        smstot(i, j) = real(1.0);
        for (int k = 0; k < nsoil; ++k) {
          smois(i, j, k) = real(1.0);
          tslb(i, j, k) = real(273.16);
        }
      } else {
        if (xice(i, j) == real(1.0)) {
          smstav(i, j) = real(1.0);
          smstot(i, j) = real(1.0);
          for (int k = 0; k < nsoil; ++k) {
            smois(i, j, k) = real(1.0);
          }
        }
      }
    }

    // Skip sea-ice and open water points (NoahmpDriverMainMod lines 150-156)
    auto set_skipped_outputs = [&]() {
      q2mvxy(i, j) = undefined_real;
      q2mbxy(i, j) = undefined_real;
      fvegxy(i, j) = undefined_real;
      if (tradxy.p) tradxy(i, j) = undefined_real;
      if (neexy.p) neexy(i, j) = undefined_real;
      if (gppxy.p) gppxy(i, j) = undefined_real;
      if (nppxy.p) nppxy(i, j) = undefined_real;
      if (runsfxy.p) runsfxy(i, j) = undefined_real;
      if (runsbxy.p) runsbxy(i, j) = undefined_real;
      if (ecanxy.p) ecanxy(i, j) = undefined_real;
      if (edirxy.p) edirxy(i, j) = undefined_real;
      if (etranxy.p) etranxy(i, j) = undefined_real;
      if (fsaxy.p) fsaxy(i, j) = undefined_real;
      if (firaxy.p) firaxy(i, j) = undefined_real;
      if (aparxy.p) aparxy(i, j) = undefined_real;
      if (psnxy.p) psnxy(i, j) = undefined_real;
      if (savxy.p) savxy(i, j) = undefined_real;
      if (sagxy.p) sagxy(i, j) = undefined_real;
      if (rssunxy.p) rssunxy(i, j) = undefined_real;
      if (rsshaxy.p) rsshaxy(i, j) = undefined_real;
      if (bgapxy.p) bgapxy(i, j) = undefined_real;
      if (wgapxy.p) wgapxy(i, j) = undefined_real;
      if (tgvxy.p) tgvxy(i, j) = undefined_real;
      if (tgbxy.p) tgbxy(i, j) = undefined_real;
      if (chvxy.p) chvxy(i, j) = undefined_real;
      if (chbxy.p) chbxy(i, j) = undefined_real;
      if (shgxy.p) shgxy(i, j) = undefined_real;
      if (shcxy.p) shcxy(i, j) = undefined_real;
      if (shbxy.p) shbxy(i, j) = undefined_real;
      if (evgxy.p) evgxy(i, j) = undefined_real;
      if (evbxy.p) evbxy(i, j) = undefined_real;
      if (ghvxy.p) ghvxy(i, j) = undefined_real;
      if (ghbxy.p) ghbxy(i, j) = undefined_real;
      if (irgxy.p) irgxy(i, j) = undefined_real;
      if (ircxy.p) ircxy(i, j) = undefined_real;
      if (irbxy.p) irbxy(i, j) = undefined_real;
      if (trxy.p) trxy(i, j) = undefined_real;
      if (evcxy.p) evcxy(i, j) = undefined_real;
      if (chleafxy.p) chleafxy(i, j) = undefined_real;
      if (chucxy.p) chucxy(i, j) = undefined_real;
      if (chv2xy.p) chv2xy(i, j) = undefined_real;
      if (chb2xy.p) chb2xy(i, j) = undefined_real;
      if (qintsxy.p) qintsxy(i, j) = undefined_real;
      if (qintrxy.p) qintrxy(i, j) = undefined_real;
      if (qdripsxy.p) qdripsxy(i, j) = undefined_real;
      if (qdriprxy.p) qdriprxy(i, j) = undefined_real;
      if (qthrosxy.p) qthrosxy(i, j) = undefined_real;
      if (qthrorxy.p) qthrorxy(i, j) = undefined_real;
      if (qsnsubxy.p) qsnsubxy(i, j) = undefined_real;
      if (qsnfroxy.p) qsnfroxy(i, j) = undefined_real;
      if (qsubcxy.p) qsubcxy(i, j) = undefined_real;
      if (qfrocxy.p) qfrocxy(i, j) = undefined_real;
      if (qevacxy.p) qevacxy(i, j) = undefined_real;
      if (qdewcxy.p) qdewcxy(i, j) = undefined_real;
      if (qfrzcxy.p) qfrzcxy(i, j) = undefined_real;
      if (qmeltcxy.p) qmeltcxy(i, j) = undefined_real;
      if (qsnbotxy.p) qsnbotxy(i, j) = undefined_real;
      if (qmeltxy.p) qmeltxy(i, j) = undefined_real;
      if (pondingxy.p) pondingxy(i, j) = real(0.0);
      if (pahxy.p) pahxy(i, j) = undefined_real;
      if (pahgxy.p) pahgxy(i, j) = undefined_real;
      if (pahvxy.p) pahvxy(i, j) = undefined_real;
      if (pahbxy.p) pahbxy(i, j) = undefined_real;
      if (fpicexy.p) fpicexy(i, j) = undefined_real;
      if (rainlsm.p) rainlsm(i, j) = undefined_real;
      if (snowlsm.p) snowlsm(i, j) = undefined_real;
      if (forctlsm.p) forctlsm(i, j) = undefined_real;
      if (forcqlsm.p) forcqlsm(i, j) = undefined_real;
      if (forcplsm.p) forcplsm(i, j) = undefined_real;
      if (forczlsm.p) forczlsm(i, j) = undefined_real;
      if (forcwlsm.p) forcwlsm(i, j) = undefined_real;
      if (eflxbxy.p) eflxbxy(i, j) = undefined_real;
      if (soilenergy.p) soilenergy(i, j) = undefined_real;
      if (snowenergy.p) snowenergy(i, j) = undefined_real;
      if (canhsxy.p) canhsxy(i, j) = undefined_real;
    };
    if (xice(i, j) >= xice_threshold) {
      for (int k = 0; k < nsoil; ++k) {
        sh2o(i, j, k) = real(1.0);
      }
      lai(i, j) = real(0.01);
      set_skipped_outputs();
      return;
    }
    if ((xland(i, j) - real(1.5)) >= real(0.0)) {
      set_skipped_outputs();
      return;
    }

    NoahmpType noahmp;
    noahmp.InitDefault();

    // ------------------------------------------------------------------------
    // 1. ConfigVarInTransfer (ConfigVarInTransferMod.F90)
    // ------------------------------------------------------------------------
    auto& dom = noahmp.config.domain;
    auto& nml = noahmp.config.nmlist;
    nml.OptDynamicVeg = opts.dveg;
    nml.OptRainSnowPartition = opts.opt_snf;
    nml.OptSoilWaterTranspiration = opts.opt_btr;
    nml.OptGroundResistanceEvap = opts.opt_rsf;
    nml.OptSurfaceDrag = opts.opt_sfc;
    nml.OptStomataResistance = opts.opt_crs;
    nml.OptSnowAlbedo = opts.opt_alb;
    nml.OptCanopyRadiationTransfer = opts.opt_rad;
    nml.OptSnowSoilTempTime = opts.opt_stc;
    nml.OptSnowThermConduct = opts.opt_tksno;
    nml.OptSoilTemperatureBottom = opts.opt_tbot;
    nml.OptSoilSupercoolWater = opts.opt_frz;
    nml.OptSoilPermeabilityFrozen = opts.opt_inf;
    nml.OptDynVicInfiltration = opts.opt_infdv;
    nml.OptTileDrainage = opts.opt_tdrn;
    nml.OptIrrigation = opts.opt_irr;
    nml.OptIrrigationMethod = opts.opt_irrm;
    nml.OptCropModel = opts.opt_crop;
    nml.OptSoilProperty = opts.opt_soil;
    nml.OptPedotransfer = opts.opt_pedo;
    nml.OptRunoffSurface = opts.opt_runsrf;
    nml.OptRunoffSubsurface = opts.opt_runsub;
    nml.OptGlacierTreatment = opts.opt_gla;
    nml.OptSnowCompaction = opts.opt_compact;
    nml.OptWetlandModel = opts.opt_wetland;
    nml.OptSnowCoverGround = opts.opt_scf;

    real dz8w_val = dz8w.p ? dz8w(i, j, 0) : dz8w_0(i, j);
    real dx_max = std::max(real(10.0), dx);

    dom.SurfaceType = 1;
    dom.NumSwRadBand = kNumRadBands;
    dom.SoilColor = 4;
    dom.NumCropGrowStage = 8;
    dom.FlagSoilProcess = calculate_soil;
    dom.NumSoilTimeStep = soil_update_steps;
    dom.NumSnowLayerMax = nsnow;
    dom.NumSnowLayerNeg = isnowxy(i, j);
    dom.NumSoilLayer = nsoil;
    dom.GridIndexI = i + 1;
    dom.GridIndexJ = j + 1;
    dom.MainTimeStep = dtbl;
    dom.SoilTimeStep = dtbl * static_cast<real>(soil_update_steps);
    dom.GridSize = std::sqrt(dx_max * dx_max);
    dom.VegType = ivgtyp(i, j);
    dom.CropType = cropcat.p ? cropcat(i, j) : 0;
    dom.IndicatorIceSfc = 0;
    dom.DayJulianInYear = julian_in;
    dom.NumDayInYear = yearlen;
    dom.Latitude = xlat(i, j);
    dom.RefHeightAboveSfc = dz8w_val * real(0.5);
    dom.ThicknessAtmosBotLayer = dz8w_val;
    dom.CosSolarZenithAngle = coszen(i, j);
    dom.IndexWaterPoint = tb.iswater;
    dom.IndexBarrenPoint = tb.isbarren;
    dom.IndexIcePoint = tb.isice;
    dom.IndexCropPoint = tb.iscrop;
    dom.IndexEBLForest = tb.eblforest;
    dom.RunoffSlopeType = 1;
    dom.DepthSoilTempBottom = tb.zbot;

    for (int k = 1; k <= nsoil; ++k) {
      dom.SoilType(k) = isltyp(i, j);
      dom.DepthSoilLayer(k) = zsoil(k);
    }
    for (int iz = -nsnow + 1; iz <= nsoil; ++iz) {
      dom.DepthSnowSoilLayer(iz) = zsnsoxy(i, j, iz + nsnow - 1);
    }

    real gvfmax_val = shdmax(i, j);
    real vegfra_val = vegfra(i, j);
    if ((ivgtyp(i, j) == tb.isurban) || (ivgtyp(i, j) > tb.urbtype_beg)) {
      if (opts.sf_urban_physics == 0) {
        dom.VegType = tb.isurban;
        dom.FlagUrban = true;
      } else {
        dom.VegType = tb.natural;
        gvfmax_val = real(0.96) * real(100.0);
      }
    }

    dom.CropType = 0;
    if ((opts.opt_crop > 0) && (ivgtyp(i, j) == tb.iscrop)) {
      dom.CropType = tb.default_crop;
    }
    if ((opts.opt_crop > 0) && cropcat.p && (cropcat(i, j) > 0)) {
      dom.CropType = cropcat(i, j);
      dom.VegType = tb.iscrop;
      vegfra_val = real(0.95) * real(100.0);
      gvfmax_val = real(0.95) * real(100.0);
    }

    bool any_water_soil = false;
    for (int k = 1; k <= nsoil; ++k) {
      if (dom.SoilType(k) == 14) any_water_soil = true;
    }
    if (any_water_soil && (xice(i, j) == real(0.0))) {
      for (int k = 1; k <= nsoil; ++k) dom.SoilType(k) = 7;
    }

    // ------------------------------------------------------------------------
    // 2. ForcingVarInTransfer (ForcingVarInTransferMod.F90)
    // ------------------------------------------------------------------------
    auto& frc = noahmp.forcing;
    real t_val = t_phy.p ? t_phy(i, j, 0) : t_phy_0(i, j);
    real u_val = u_phy.p ? u_phy(i, j, 0) : u_phy_0(i, j);
    real v_val = v_phy.p ? v_phy(i, j, 0) : v_phy_0(i, j);
    real qv_val = qv_curr.p ? qv_curr(i, j, 0) : qv_curr_0(i, j);
    real p8w_k0 = p8w.p ? p8w(i, j, 0) : p8w_0(i, j);
    real p8w_k1 = p8w.p ? p8w(i, j, 1) : (p8w_1.p ? p8w_1(i, j) : p8w_0(i, j));

    frc.TemperatureAirRefHeight = t_val;
    frc.WindEastwardRefHeight = u_val;
    frc.WindNorthwardRefHeight = v_val;
    frc.SpecHumidityRefHeight = qv_val / (real(1.0) + qv_val);
    frc.PressureAirRefHeight = (p8w_k0 + p8w_k1) * real(0.5);
    frc.PressureAirSurface = p8w_k0;
    frc.RadLwDownRefHeight = glw(i, j);
    frc.RadSwDownRefHeight = swdown(i, j);
    frc.TemperatureSoilBottom = tmn(i, j);

    real precip_total = rainbl(i, j) / dtbl;
    frc.PrecipConvRefHeight = raincv.p ? (raincv(i, j) / dtbl) : real(0.0);
    frc.PrecipNonConvRefHeight = rainncv.p ? (rainncv(i, j) / dtbl) : real(0.0);
    frc.PrecipShConvRefHeight = rainshv.p ? (rainshv(i, j) / dtbl) : real(0.0);
    frc.PrecipSnowRefHeight = snowncv.p ? (snowncv(i, j) / dtbl) : real(0.0);
    frc.PrecipGraupelRefHeight = graupelncv.p ? (graupelncv(i, j) / dtbl) : real(0.0);
    frc.PrecipHailRefHeight = hailncv.p ? (hailncv(i, j) / dtbl) : real(0.0);

    real precip_other = precip_total - frc.PrecipConvRefHeight -
                        frc.PrecipNonConvRefHeight - frc.PrecipShConvRefHeight;
    precip_other = std::max(real(0.0), precip_other);
    frc.PrecipNonConvRefHeight = frc.PrecipNonConvRefHeight + precip_other;
    frc.PrecipSnowRefHeight = frc.PrecipSnowRefHeight + precip_other * sr(i, j);

    frc.RadSwDirFrac = real(0.7);
    frc.RadSwVisFrac = real(0.5);

    // ------------------------------------------------------------------------
    // 3. EnergyVarInTransfer (EnergyVarInTransferMod.F90)
    // ------------------------------------------------------------------------
    auto& es = noahmp.energy.state;
    auto& ef = noahmp.energy.flux;
    auto& ep = noahmp.energy.param;
    const int vt = dom.VegType;

    es.LeafAreaIndex = lai(i, j);
    es.StemAreaIndex = xsaixy(i, j);
    es.SpecHumiditySfcMean = qsfc(i, j);
    es.TemperatureGrd = tgxy(i, j);
    es.TemperatureCanopy = tvxy(i, j);
    es.SnowAgeNondim = taussxy(i, j);
    es.AlbedoSnowPrev = alboldxy(i, j);
    es.PressureVaporCanAir = eahxy(i, j);
    es.TemperatureCanopyAir = tahxy(i, j);
    es.ExchCoeffShSfc = chxy(i, j);
    es.ExchCoeffMomSfc = cmxy(i, j);
    for (int iz = -nsnow + 1; iz <= 0; ++iz) {
      es.TemperatureSoilSnow(iz) = tsnoxy(i, j, iz + nsnow - 1);
    }
    for (int k = 1; k <= nsoil; ++k) {
      es.TemperatureSoilSnow(k) = tslb(i, j, k - 1);
    }
    es.PressureAtmosCO2 = tb.co2 * frc.PressureAirRefHeight;
    es.PressureAtmosO2 = tb.o2 * frc.PressureAirRefHeight;
    for (int ib = 1; ib <= kNumRadBands; ++ib) {
      es.AlbedoSoilDir(ib) = albsoildirxy.p ? albsoildirxy(i, j, ib - 1) : real(0.0);
      es.AlbedoSoilDif(ib) = albsoildifxy.p ? albsoildifxy(i, j, ib - 1) : real(0.0);
    }
    if ((vt == 25) || (vt == 26) || (vt == 27)) {
      es.VegFrac = real(0.0);
      es.LeafAreaIndex = real(0.0);
    }
    ef.HeatGroundTotAcc = acc_ssoilxy.p ? acc_ssoilxy(i, j) : real(0.0);

    ep.SoilHeatCapacity = tb.csoil;
    ep.SnowAgeFacBats = tb.tau0;
    ep.SnowGrowVapFacBats = tb.grain_growth;
    ep.SnowSootFacBats = tb.dirt_soot;
    ep.SnowGrowFrzFacBats = tb.extra_growth;
    ep.SolarZenithAdjBats = tb.bats_cosz;
    ep.FreshSnoAlbVisBats = tb.bats_vis_new;
    ep.FreshSnoAlbNirBats = tb.bats_nir_new;
    ep.SnoAgeFacDifVisBats = tb.bats_vis_age;
    ep.SnoAgeFacDifNirBats = tb.bats_nir_age;
    ep.SzaFacDirVisBats = tb.bats_vis_dir;
    ep.SzaFacDirNirBats = tb.bats_nir_dir;
    ep.SnowAlbRefClass = tb.class_alb_ref;
    ep.SnowAgeFacClass = tb.class_sno_age;
    ep.SnowAlbFreshClass = tb.class_alb_new;
    ep.UpscatterCoeffSnowDir = tb.betads;
    ep.UpscatterCoeffSnowDif = tb.betais;
    ep.ZilitinkevichCoeff = tb.czil;
    ep.EmissivitySnow = tb.snow_emis;
    for (int ib = 1; ib <= 2; ++ib) {
      ep.EmissivitySoilLake(ib) = tb.eg[ib];
    }
    for (int ib = 1; ib <= kNumRadBands; ++ib) {
      ep.AlbedoLandIce(ib) = tb.albice[ib];
    }
    ep.RoughLenMomSnow = tb.z0sno;
    ep.RoughLenMomSoil = tb.z0soil;
    ep.RoughLenMomLake = tb.z0lake;
    ep.EmissivityIceSfc = tb.eice;
    ep.ResistanceSoilExp = tb.rsurf_exp;
    ep.ResistanceSnowSfc = tb.rsurf_snow;
    ep.VegFracAnnMax = gvfmax_val / real(100.0);
    ep.VegFracGreen = vegfra_val / real(100.0);
    ep.TreeCrownRadius = tb.rc[vt];
    ep.HeightCanopyTop = tb.hvt[vt];
    ep.HeightCanopyBot = tb.hvb[vt];
    ep.RoughLenMomVeg = tb.z0mvt[vt];
    ep.CanopyWindExtFac = tb.cwpvt[vt];
    ep.TreeDensity = tb.den[vt];
    ep.CanopyOrientIndex = tb.xl[vt];
    ep.ConductanceLeafMin = tb.bp[vt];
    ep.Co2MmConst25C = tb.kc25[vt];
    ep.O2MmConst25C = tb.ko25[vt];
    ep.Co2MmConstQ10 = tb.akc[vt];
    ep.O2MmConstQ10 = tb.ako[vt];
    ep.RadiationStressFac = tb.rgl[vt];
    ep.ResistanceStomataMin = tb.rs[vt];
    ep.ResistanceStomataMax = tb.rsmax[vt];
    ep.AirTempOptimTransp = tb.topt[vt];
    ep.VaporPresDeficitFac = tb.hs[vt];
    ep.LeafDimLength = tb.dleaf[vt];
    ep.HeatCapacCanFac = tb.cbiom[vt];
    for (int m = 1; m <= 12; ++m) {
      ep.LeafAreaIndexMon(m) = tb.laim[vt][m];
      ep.StemAreaIndexMon(m) = tb.saim[vt][m];
    }
    for (int ib = 1; ib <= kNumRadBands; ++ib) {
      ep.ReflectanceLeaf(ib) = tb.rhol[vt][ib];
      ep.ReflectanceStem(ib) = tb.rhos[vt][ib];
      ep.TransmittanceLeaf(ib) = tb.taul[vt][ib];
      ep.TransmittanceStem(ib) = tb.taus[vt][ib];
      ep.AlbedoSoilSat(ib) = tb.albsat[dom.SoilColor][ib];
      ep.AlbedoSoilDry(ib) = tb.albdry[dom.SoilColor][ib];
      ep.AlbedoLakeFrz(ib) = tb.alblak[ib];
      ep.ScatterCoeffSnow(ib) = tb.omegas[ib];
    }
    for (int k = 1; k <= nsoil; ++k) {
      ep.SoilQuartzFrac(k) = tb.quartz[dom.SoilType(k)];
    }
    if (dom.FlagUrban) ep.SoilHeatCapacity = real(3.0e6);
    if (dom.CropType > 0) {
      ep.ConductanceLeafMin = tb.bpi[dom.CropType];
      ep.Co2MmConst25C = tb.kc25i[dom.CropType];
      ep.O2MmConst25C = tb.ko25i[dom.CropType];
      ep.Co2MmConstQ10 = tb.akci[dom.CropType];
      ep.O2MmConstQ10 = tb.akoi[dom.CropType];
    }

    // ------------------------------------------------------------------------
    // 4. WaterVarInTransfer (WaterVarInTransferMod.F90)
    // ------------------------------------------------------------------------
    auto& ws = noahmp.water.state;
    auto& wf = noahmp.water.flux;
    auto& wp = noahmp.water.param;
    ws.CanopyLiqWater = canliqxy(i, j);
    ws.CanopyIce = canicexy(i, j);
    ws.CanopyWetFrac = fwetxy(i, j);
    ws.SnowWaterEquiv = snow(i, j);
    ws.SnowWaterEquivPrev = sneqvoxy(i, j);
    ws.SnowDepth = snowh(i, j);
    ws.IrrigationFracFlood = real(0.0);
    ws.IrrigationAmtFlood = real(0.0);
    ws.IrrigationFracMicro = real(0.0);
    ws.IrrigationAmtMicro = real(0.0);
    ws.IrrigationFracSprinkler = real(0.0);
    ws.IrrigationAmtSprinkler = real(0.0);
    ws.WaterTableDepth = zwtxy(i, j);
    ws.SoilMoistureToWT = smcwtdxy(i, j);
    ws.TileDrainFrac = undefined_real;
    ws.WaterStorageAquifer = waxy(i, j);
    ws.WaterStorageSoilAqf = wtxy(i, j);
    ws.WaterStorageLake = wslakexy(i, j);
    ws.IrrigationFracGrid = real(0.0);
    ws.IrrigationCntSprinkler = 0;
    ws.IrrigationCntMicro = 0;
    ws.IrrigationCntFlood = 0;
    for (int iz = -nsnow + 1; iz <= 0; ++iz) {
      ws.SnowIce(iz) = snicexy(i, j, iz + nsnow - 1);
      ws.SnowLiqWater(iz) = snliqxy(i, j, iz + nsnow - 1);
    }
    for (int k = 1; k <= nsoil; ++k) {
      ws.SoilLiqWater(k) = sh2o(i, j, k - 1);
      ws.SoilMoisture(k) = smois(i, j, k - 1);
      ws.SoilMoistureEqui(k) = smoiseq.p ? smoiseq(i, j, k - 1) : undefined_real;
    }
    ws.RechargeGwDeepWT = real(0.0);
    ws.RechargeGwShallowWT = real(0.0);

    wf.EvapSoilSfcLiqAcc = acc_qsevaxy.p ? acc_qsevaxy(i, j) : real(0.0);
    wf.SoilSfcInflowAcc = acc_qinsurxy.p ? acc_qinsurxy(i, j) : real(0.0);
    wf.SfcWaterTotChgAcc = acc_dwaterxy.p ? acc_dwaterxy(i, j) : real(0.0);
    wf.PrecipTotAcc = acc_prcpxy.p ? acc_prcpxy(i, j) : real(0.0);
    wf.EvapCanopyNetAcc = acc_ecanxy.p ? acc_ecanxy(i, j) : real(0.0);
    wf.TranspirationAcc = acc_etranxy.p ? acc_etranxy(i, j) : real(0.0);
    wf.EvapGroundNetAcc = acc_edirxy.p ? acc_edirxy(i, j) : real(0.0);
    for (int k = 1; k <= nsoil; ++k) {
      wf.TranspWatLossSoilAcc(k) = acc_etranixy.p ? acc_etranixy(i, j, k - 1) : real(0.0);
    }
    wf.GlacierExcessFlowAcc = acc_glaflwxy.p ? acc_glaflwxy(i, j) : real(0.0);

    wp.DrainSoilLayerInd = tb.drain_layer_opt;
    wp.CanopyLiqHoldCap = tb.ch2op[vt];
    wp.SnowCompactBurdenFac = tb.c2_snowcompact;
    wp.SnowCompactAgingFac1 = tb.c3_snowcompact;
    wp.SnowCompactAgingFac2 = tb.c4_snowcompact;
    wp.SnowCompactAgingFac3 = tb.c5_snowcompact;
    wp.SnowCompactAgingMax = tb.dm_snowcompact;
    wp.SnowViscosityCoeff = tb.eta0_snowcompact;
    wp.SnowCompactmAR24 = tb.snowcompactm_ar24;
    wp.SnowCompactbAR24 = tb.snowcompactb_ar24;
    wp.SnowCompactP1AR24 = tb.snowcompact_p1_ar24;
    wp.SnowCompactP2AR24 = tb.snowcompact_p2_ar24;
    wp.SnowCompactP3AR24 = tb.snowcompact_p3_ar24;
    wp.SnowCoverM1AR25 = tb.scfm1_ar25;
    wp.SnowCoverM2AR25 = tb.scfm2_ar25;
    wp.SnowCoverFac1AR25 = tb.scfac1_ar25;
    wp.SnowCoverFac2AR25 = tb.scfac2_ar25;
    wp.BurdenFacUpAR24 = tb.snowcompact_up_ar24;
    wp.SnowLiqFracMax = tb.snliqmaxfrac;
    wp.SnowLiqHoldCap = tb.ssi;
    wp.SnowLiqReleaseFac = tb.snow_ret_fac;
    wp.IrriFloodRateFac = tb.firtfac;
    wp.IrriMicroRate = tb.micir_rate;
    wp.SoilConductivityRef = tb.refdk;
    wp.SoilInfilFacRef = tb.refkdt;
    wp.GroundFrzCoeff = tb.frzk;
    wp.GridTopoIndex = tb.timean;
    wp.SoilSfcSatFracMax = tb.fsatmx;
    wp.SpecYieldGw = tb.rous;
    wp.MicroPoreContent = tb.cmic;
    wp.WaterStorageLakeMax = tb.wslmax;
    wp.SnoWatEqvMaxGlacier = tb.swemaxgla;
    wp.IrriStopDayBfHarvest = tb.irr_har;
    wp.IrriTriggerLaiMin = tb.irr_lai;
    wp.SoilWatDeficitAllow = tb.irr_mad;
    wp.IrriFloodLossFrac = tb.filoss;
    wp.IrriSprinklerRate = tb.sprir_rate;
    wp.IrriFracThreshold = tb.irr_frac;
    wp.IrriStopPrecipThr = tb.ir_rain;
    wp.SnowfallDensityMax = tb.snowden_max;
    wp.SnowMassFullCoverOld = tb.swemx;
    wp.SoilMatPotentialWilt = tb.psiwlt;
    wp.SnowMeltFac = tb.mfsno[vt];
    wp.SnowCoverFac = tb.scffac[vt];
    const int st1 = dom.SoilType(1);
    wp.InfilFacVic = tb.bvic[st1];
    wp.TensionWatDistrInfl = tb.axaj[st1];
    wp.TensionWatDistrShp = tb.bxaj[st1];
    wp.FreeWatDistrShp = tb.xxaj[st1];
    wp.InfilHeteroDynVic = tb.bbvic[st1];
    wp.InfilCapillaryDynVic = tb.gdvic[st1];
    wp.InfilFacDynVic = tb.bdvic[st1];
    wp.TileDrainCoeffSp = tb.td_dc[st1];
    wp.TileDrainTubeDepth = tb.td_depth[st1];
    wp.DrainFacSoilWat = tb.tdsmc_fac[st1];
    wp.TileDrainCoeff = tb.td_dcoef[st1];
    wp.DrainDepthToImperv = tb.td_adepth[st1];
    wp.LateralWatCondFac = tb.klat_fac[st1];
    wp.TileDrainDepth = tb.td_ddrain[st1];
    wp.DrainTubeDist = tb.td_spac[st1];
    wp.DrainTubeRadius = tb.td_radi[st1];
    wp.DrainWatDepToImperv = tb.td_d[st1];
    wp.NumSoilLayerRoot = static_cast<int>(tb.nroot[vt]);
    wp.SoilDrainSlope = tb.slope[dom.RunoffSlopeType];
    wp.WetlandCapMax = tb.wcap;

    for (int k = 1; k <= nsoil; ++k) {
      int st_k = dom.SoilType(k);
      wp.SoilMoistureSat(k) = tb.smcmax[st_k];
      wp.SoilMoistureWilt(k) = tb.smcwlt[st_k];
      wp.SoilMoistureFieldCap(k) = tb.smcref[st_k];
      wp.SoilMoistureDry(k) = tb.smcdry[st_k];
      wp.SoilWatDiffusivitySat(k) = tb.dwsat[st_k];
      wp.SoilWatConductivitySat(k) = tb.dksat[st_k];
      wp.SoilExpCoeffB(k) = tb.bexp[st_k];
      wp.SoilMatPotentialSat(k) = tb.psisat[st_k];
    }

    wp.SoilInfilMaxCoeff = wp.SoilInfilFacRef * wp.SoilWatConductivitySat(1) / wp.SoilConductivityRef;
    if (dom.FlagUrban) {
      wp.SoilMoistureSat.fill(real(0.45));
      wp.SoilMoistureFieldCap.fill(real(0.42));
      wp.SoilMoistureWilt.fill(real(0.40));
      wp.SoilMoistureDry.fill(real(0.40));
    }
    if (st1 != 14) {
      wp.SoilImpervFracCoeff = wp.GroundFrzCoeff *
                               ((wp.SoilMoistureSat(1) / wp.SoilMoistureFieldCap(1)) *
                                (real(0.412) / real(0.468)));
    }
    ws.SnowIceFracPrev.fill(real(0.0));
    for (int iz = dom.NumSnowLayerNeg + 1; iz <= 0; ++iz) {
      real sni = snicexy(i, j, iz + nsnow - 1);
      real snl = snliqxy(i, j, iz + nsnow - 1);
      ws.SnowIceFracPrev(iz) = sni / (sni + snl);
    }

    // ------------------------------------------------------------------------
    // 5. BiochemVarInTransfer (BiochemVarInTransferMod.F90)
    // ------------------------------------------------------------------------
    auto& bs = noahmp.biochem.state;
    auto& bp = noahmp.biochem.param;
    bs.PlantGrowStage = pgsxy.p ? pgsxy(i, j) : undefined_int;
    bs.LeafMass = lfmassxy(i, j);
    bs.RootMass = rtmassxy(i, j);
    bs.StemMass = stmassxy(i, j);
    bs.WoodMass = woodxy(i, j);
    bs.CarbonMassDeepSoil = stblcpxy(i, j);
    bs.CarbonMassShallowSoil = fastcpxy(i, j);
    bs.GrainMass = grainxy.p ? grainxy(i, j) : undefined_real;
    bs.GrowDegreeDay = gddxy.p ? gddxy(i, j) : undefined_real;
    bs.NitrogenConcFoliage = real(1.0);

    bp.NitrogenConcFoliageMax = tb.folnmx[vt];
    bp.QuantumEfficiency25C = tb.qe25[vt];
    bp.CarboxylRateMax25C = tb.vcmx25[vt];
    bp.CarboxylRateMaxQ10 = tb.avcmx[vt];
    bp.PhotosynPathC3 = tb.c3psn[vt];
    bp.SlopeConductToPhotosyn = tb.mp[vt];
    bp.RespMaintQ10 = tb.arm[vt];
    bp.RespMaintLeaf25C = tb.rmf25[vt];
    bp.RespMaintStem25C = tb.rms25[vt];
    bp.RespMaintRoot25C = tb.rmr25[vt];
    bp.WoodToRootRatio = tb.wrrat[vt];
    bp.WoodPoolIndex = tb.wdpool[vt];
    bp.TurnoverCoeffLeafVeg = tb.ltovrc[vt];
    bp.TemperaureLeafFreeze = tb.tdlef[vt];
    bp.LeafDeathWaterCoeffVeg = tb.dilefw[vt];
    bp.LeafDeathTempCoeffVeg = tb.dilefc[vt];
    bp.GrowthRespFrac = tb.fragr[vt];
    bp.MicroRespCoeff = tb.mrp[vt];
    bp.TemperatureMinPhotosyn = tb.tmin[vt];
    bp.LeafAreaPerMass1side = tb.sla[vt];
    bp.StemAreaIndexMin = tb.xsamin[vt];
    bp.WoodAllocFac = tb.bf[vt];
    bp.WaterStressCoeff = tb.wstrc[vt];
    bp.LeafAreaIndexMin = tb.laimin[vt];
    bp.TurnoverCoeffRootVeg = tb.rtovrc[vt];
    bp.WoodRespCoeff = tb.rswoodc[vt];

    // ------------------------------------------------------------------------
    // 6. Execute NoahmpMain or NoahmpMainGlacier (NoahmpDriverMainMod lines 198-206)
    // ------------------------------------------------------------------------
    if (dom.VegType == dom.IndexIcePoint) {
      dom.IndicatorIceSfc = -1;
      frc.TemperatureSoilBottom = std::min(frc.TemperatureSoilBottom, real(263.15));
      NoahmpMainGlacier(noahmp);
    } else {
      dom.IndicatorIceSfc = 0;
      NoahmpMain(noahmp);
    }

    // ------------------------------------------------------------------------
    // 7. Transfer Outputs Back to Views (*VarOutTransferMod.F90)
    // ------------------------------------------------------------------------
    const int ice = dom.IndicatorIceSfc;

    // ConfigVarOutTransfer
    isnowxy(i, j) = dom.NumSnowLayerNeg;
    for (int iz = -nsnow + 1; iz <= nsoil; ++iz) {
      zsnsoxy(i, j, iz + nsnow - 1) = dom.DepthSnowSoilLayer(iz);
    }
    if (forczlsm.p) forczlsm(i, j) = dom.RefHeightAboveSfc;

    // ForcingVarOutTransfer
    if (forctlsm.p) forctlsm(i, j) = frc.TemperatureAirRefHeight;
    if (forcqlsm.p) forcqlsm(i, j) = frc.SpecHumidityRefHeight;
    if (forcplsm.p) forcplsm(i, j) = frc.PressureAirRefHeight;
    if (forcwlsm.p) {
      forcwlsm(i, j) = std::sqrt(frc.WindEastwardRefHeight * frc.WindEastwardRefHeight +
                                 frc.WindNorthwardRefHeight * frc.WindNorthwardRefHeight);
    }

    // EnergyVarOutTransfer
    if (ice == -1) {
      es.VegFrac = real(0.0);
      es.RoughLenMomSfcToAtm = real(0.002);
      ef.RadSwAbsVeg = real(0.0);
      ef.RadLwNetCanopy = real(0.0);
      ef.RadLwNetVegGrd = real(0.0);
      ef.HeatSensibleCanopy = real(0.0);
      ef.HeatSensibleVegGrd = real(0.0);
      ef.HeatLatentVegGrd = real(0.0);
      ef.HeatGroundVegGrd = real(0.0);
      ef.HeatCanStorageChg = real(0.0);
      ef.HeatLatentCanTransp = real(0.0);
      ef.HeatLatentCanEvap = real(0.0);
      ef.HeatPrecipAdvCanopy = real(0.0);
      ef.HeatPrecipAdvVegGrd = real(0.0);
      ef.HeatLatentCanopy = real(0.0);
      ef.HeatLatentTransp = real(0.0);
      ef.RadLwNetBareGrd = ef.RadLwNetSfc;
      ef.HeatSensibleBareGrd = ef.HeatSensibleSfc;
      ef.HeatLatentBareGrd = ef.HeatLatentGrd;
      ef.HeatGroundBareGrd = ef.HeatGroundTot;
      es.TemperatureGrdBare = es.TemperatureGrd;
      es.ExchCoeffShBare = es.ExchCoeffShSfc;
      lh(i, j) = ef.HeatLatentGrd;
    } else if (ice == 0) {
      lh(i, j) = ef.HeatLatentGrd + ef.HeatLatentCanopy + ef.HeatLatentTransp + ef.HeatLatentIrriEvap;
    }

    hfx(i, j) = ef.HeatSensibleSfc;
    grdflx(i, j) = ef.HeatGroundTot;
    if (fsaxy.p) fsaxy(i, j) = ef.RadSwAbsSfc;
    if (firaxy.p) firaxy(i, j) = ef.RadLwNetSfc;
    if (aparxy.p) aparxy(i, j) = ef.RadPhotoActAbsCan;
    if (savxy.p) savxy(i, j) = ef.RadSwAbsVeg;
    if (sagxy.p) sagxy(i, j) = ef.RadSwAbsGrd;
    if (ircxy.p) ircxy(i, j) = ef.RadLwNetCanopy;
    if (irgxy.p) irgxy(i, j) = ef.RadLwNetVegGrd;
    if (shcxy.p) shcxy(i, j) = ef.HeatSensibleCanopy;
    if (shgxy.p) shgxy(i, j) = ef.HeatSensibleVegGrd;
    if (evgxy.p) evgxy(i, j) = ef.HeatLatentVegGrd;
    if (ghvxy.p) ghvxy(i, j) = ef.HeatGroundVegGrd;
    if (irbxy.p) irbxy(i, j) = ef.RadLwNetBareGrd;
    if (shbxy.p) shbxy(i, j) = ef.HeatSensibleBareGrd;
    if (evbxy.p) evbxy(i, j) = ef.HeatLatentBareGrd;
    if (ghbxy.p) ghbxy(i, j) = ef.HeatGroundBareGrd;
    if (trxy.p) trxy(i, j) = ef.HeatLatentCanTransp;
    if (evcxy.p) evcxy(i, j) = ef.HeatLatentCanEvap;
    if (canhsxy.p) canhsxy(i, j) = ef.HeatCanStorageChg;
    if (pahxy.p) pahxy(i, j) = ef.HeatPrecipAdvSfc;
    if (pahgxy.p) pahgxy(i, j) = ef.HeatPrecipAdvVegGrd;
    if (pahvxy.p) pahvxy(i, j) = ef.HeatPrecipAdvCanopy;
    if (pahbxy.p) pahbxy(i, j) = ef.HeatPrecipAdvBareGrd;
    if (acc_ssoilxy.p) acc_ssoilxy(i, j) = ef.HeatGroundTotAcc;
    if (eflxbxy.p) eflxbxy(i, j) = ef.HeatFromSoilBot;

    tsk(i, j) = es.TemperatureRadSfc;
    emiss(i, j) = es.EmissivitySfc;
    qsfc(i, j) = es.SpecHumiditySfcMean;
    tvxy(i, j) = es.TemperatureCanopy;
    tgxy(i, j) = es.TemperatureGrd;
    eahxy(i, j) = es.PressureVaporCanAir;
    tahxy(i, j) = es.TemperatureCanopyAir;
    cmxy(i, j) = es.ExchCoeffMomSfc;
    chxy(i, j) = es.ExchCoeffShSfc;
    alboldxy(i, j) = es.AlbedoSnowPrev;
    lai(i, j) = es.LeafAreaIndex;
    xsaixy(i, j) = es.StemAreaIndex;
    taussxy(i, j) = es.SnowAgeNondim;
    z0(i, j) = es.RoughLenMomSfcToAtm;
    znt(i, j) = es.RoughLenMomSfcToAtm;
    t2mvxy(i, j) = es.TemperatureAir2mVeg;
    t2mbxy(i, j) = es.TemperatureAir2mBare;
    if (tradxy.p) tradxy(i, j) = es.TemperatureRadSfc;
    fvegxy(i, j) = es.VegFrac;
    if (rssunxy.p) rssunxy(i, j) = es.ResistanceStomataSunlit;
    if (rsshaxy.p) rsshaxy(i, j) = es.ResistanceStomataShade;
    if (bgapxy.p) bgapxy(i, j) = es.GapBtwCanopy;
    if (wgapxy.p) wgapxy(i, j) = es.GapInCanopy;
    if (tgvxy.p) tgvxy(i, j) = es.TemperatureGrdVeg;
    if (tgbxy.p) tgbxy(i, j) = es.TemperatureGrdBare;
    if (chvxy.p) chvxy(i, j) = es.ExchCoeffShAbvCan;
    if (chbxy.p) chbxy(i, j) = es.ExchCoeffShBare;
    if (chleafxy.p) chleafxy(i, j) = es.ExchCoeffShLeaf;
    if (chucxy.p) chucxy(i, j) = es.ExchCoeffShUndCan;
    if (chv2xy.p) chv2xy(i, j) = es.ExchCoeffSh2mVeg;
    if (chb2xy.p) chb2xy(i, j) = es.ExchCoeffSh2mBare;
    q2mvxy(i, j) = es.SpecHumidity2mVeg / (real(1.0) - es.SpecHumidity2mVeg);
    q2mbxy(i, j) = es.SpecHumidity2mBare / (real(1.0) - es.SpecHumidity2mBare);
    albedo(i, j) = es.AlbedoSfc;

    for (int k = 1; k <= nsoil; ++k) {
      tslb(i, j, k - 1) = es.TemperatureSoilSnow(k);
    }
    for (int iz = -nsnow + 1; iz <= 0; ++iz) {
      tsnoxy(i, j, iz + nsnow - 1) = es.TemperatureSoilSnow(iz);
    }
    if (albsoildirxy.p && albsoildifxy.p) {
      for (int ib = 1; ib <= kNumRadBands; ++ib) {
        albsoildirxy(i, j, ib - 1) = es.AlbedoSoilDir(ib);
        albsoildifxy(i, j, ib - 1) = es.AlbedoSoilDif(ib);
      }
    }

    if (rs.p) {
      real lai_sun = std::max(es.LeafAreaIndSunlit, real(0.0));
      real lai_sha = std::max(es.LeafAreaIndShade, real(0.0));
      real rb_leaf = std::max(es.ResistanceLeafBoundary, real(0.0));
      if ((es.ResistanceStomataSunlit <= real(0.0)) || (es.ResistanceStomataShade <= real(0.0)) ||
          (lai_sun == real(0.0)) || (lai_sha == real(0.0)) ||
          (es.ResistanceStomataSunlit == undefined_real) ||
          (es.ResistanceStomataShade == undefined_real)) {
        rs(i, j) = real(0.0);
      } else {
        real cond = ((real(1.0) / (es.ResistanceStomataSunlit + rb_leaf)) * es.LeafAreaIndSunlit) +
                    ((real(1.0) / (es.ResistanceStomataShade + rb_leaf)) * es.LeafAreaIndShade);
        rs(i, j) = real(1.0) / cond;
      }
    }

    if (snowenergy.p && soilenergy.p) {
      snowenergy(i, j) = real(0.0);
      soilenergy(i, j) = real(0.0);
      for (int iz = dom.NumSnowLayerNeg + 1; iz <= nsoil; ++iz) {
        real dz_iz = (iz == dom.NumSnowLayerNeg + 1)
                         ? (-dom.DepthSnowSoilLayer(iz))
                         : (dom.DepthSnowSoilLayer(iz - 1) - dom.DepthSnowSoilLayer(iz));
        if (iz >= 1) {
          soilenergy(i, j) = soilenergy(i, j) + dz_iz * es.HeatCapacSoilSnow(iz) *
                             (es.TemperatureSoilSnow(iz) - real(273.16)) * real(0.001);
        } else {
          snowenergy(i, j) = snowenergy(i, j) + dz_iz * es.HeatCapacSoilSnow(iz) *
                             (es.TemperatureSoilSnow(iz) - real(273.16)) * real(0.001);
        }
      }
    }

    // WaterVarOutTransfer
    if (ice == -1) {
      ws.SnowCoverFrac = real(1.0);
      wf.EvapCanopyNet = real(0.0);
      wf.Transpiration = real(0.0);
      wf.InterceptCanopySnow = real(0.0);
      wf.InterceptCanopyRain = real(0.0);
      wf.DripCanopySnow = real(0.0);
      wf.DripCanopyRain = real(0.0);
      wf.ThroughfallSnow = wf.SnowfallRefHeight;
      wf.ThroughfallRain = wf.RainfallRefHeight;
      wf.SublimCanopyIce = real(0.0);
      wf.FrostCanopyIce = real(0.0);
      wf.FreezeCanopyLiq = real(0.0);
      wf.MeltCanopyIce = real(0.0);
      wf.EvapCanopyLiq = real(0.0);
      wf.DewCanopyLiq = real(0.0);
      ws.CanopyIce = real(0.0);
      ws.CanopyLiqWater = real(0.0);
      wf.TileDrain = real(0.0);
      wf.RunoffSurface = wf.RunoffSurface * dom.MainTimeStep;
      wf.RunoffSubsurface = wf.RunoffSubsurface * dom.MainTimeStep;
      qfx(i, j) = wf.EvapGroundNet;
    } else if (ice == 0) {
      qfx(i, j) = wf.EvapCanopyNet + wf.EvapGroundNet + wf.Transpiration + wf.EvapIrriSprinkler;
    }

    smstav(i, j) = real(0.0);
    smstot(i, j) = real(0.0);
    sfcrunoff(i, j) = sfcrunoff(i, j) + wf.RunoffSurface;
    udrunoff(i, j) = udrunoff(i, j) + wf.RunoffSubsurface;
    if (qtdrain.p) qtdrain(i, j) = qtdrain(i, j) + wf.TileDrain;
    snowc(i, j) = ws.SnowCoverFrac;
    snow(i, j) = ws.SnowWaterEquiv;
    snowh(i, j) = ws.SnowDepth;
    canwat(i, j) = ws.CanopyLiqWater + ws.CanopyIce;
    acsnow(i, j) = acsnow(i, j) + rainbl(i, j) * ws.FrozenPrecipFrac;
    acsnom(i, j) = acsnom(i, j) + wf.MeltGroundSnow * dtbl;
    canliqxy(i, j) = ws.CanopyLiqWater;
    canicexy(i, j) = ws.CanopyIce;
    fwetxy(i, j) = ws.CanopyWetFrac;
    sneqvoxy(i, j) = ws.SnowWaterEquivPrev;
    qsnowxy(i, j) = wf.SnowfallGround;
    qrainxy(i, j) = wf.RainfallGround;
    wslakexy(i, j) = ws.WaterStorageLake;
    zwtxy(i, j) = ws.WaterTableDepth;
    waxy(i, j) = ws.WaterStorageAquifer;
    wtxy(i, j) = ws.WaterStorageSoilAqf;
    if (runsfxy.p) runsfxy(i, j) = wf.RunoffSurface;
    if (runsbxy.p) runsbxy(i, j) = wf.RunoffSubsurface;
    if (ecanxy.p) ecanxy(i, j) = wf.EvapCanopyNet;
    if (edirxy.p) edirxy(i, j) = wf.EvapGroundNet;
    if (etranxy.p) etranxy(i, j) = wf.Transpiration;
    if (qintsxy.p) qintsxy(i, j) = wf.InterceptCanopySnow;
    if (qintrxy.p) qintrxy(i, j) = wf.InterceptCanopyRain;
    if (qdripsxy.p) qdripsxy(i, j) = wf.DripCanopySnow;
    if (qdriprxy.p) qdriprxy(i, j) = wf.DripCanopyRain;
    if (qthrosxy.p) qthrosxy(i, j) = wf.ThroughfallSnow;
    if (qthrorxy.p) qthrorxy(i, j) = wf.ThroughfallRain;
    if (qsnsubxy.p) qsnsubxy(i, j) = wf.SublimSnowSfcIce;
    if (qsnfroxy.p) qsnfroxy(i, j) = wf.FrostSnowSfcIce;
    if (qsubcxy.p) qsubcxy(i, j) = wf.SublimCanopyIce;
    if (qfrocxy.p) qfrocxy(i, j) = wf.FrostCanopyIce;
    if (qevacxy.p) qevacxy(i, j) = wf.EvapCanopyLiq;
    if (qdewcxy.p) qdewcxy(i, j) = wf.DewCanopyLiq;
    if (qfrzcxy.p) qfrzcxy(i, j) = wf.FreezeCanopyLiq;
    if (qmeltcxy.p) qmeltcxy(i, j) = wf.MeltCanopyIce;
    if (qsnbotxy.p) qsnbotxy(i, j) = wf.SnowBotOutflow;
    if (qmeltxy.p) qmeltxy(i, j) = wf.MeltGroundSnow;
    if (pondingxy.p) pondingxy(i, j) = ws.PondSfcThinSnwTrans + ws.PondSfcThinSnwComb + ws.PondSfcThinSnwMelt;
    if (fpicexy.p) fpicexy(i, j) = ws.FrozenPrecipFrac;
    if (rainlsm.p) rainlsm(i, j) = wf.RainfallRefHeight;
    if (snowlsm.p) snowlsm(i, j) = wf.SnowfallRefHeight;
    if (acc_qinsurxy.p) acc_qinsurxy(i, j) = wf.SoilSfcInflowAcc;
    if (acc_qsevaxy.p) acc_qsevaxy(i, j) = wf.EvapSoilSfcLiqAcc;
    if (acc_dwaterxy.p) acc_dwaterxy(i, j) = wf.SfcWaterTotChgAcc;
    if (acc_prcpxy.p) acc_prcpxy(i, j) = wf.PrecipTotAcc;
    if (acc_ecanxy.p) acc_ecanxy(i, j) = wf.EvapCanopyNetAcc;
    if (acc_etranxy.p) acc_etranxy(i, j) = wf.TranspirationAcc;
    if (acc_edirxy.p) acc_edirxy(i, j) = wf.EvapGroundNetAcc;
    if (acc_glaflwxy.p) acc_glaflwxy(i, j) = wf.GlacierExcessFlowAcc;
    rechxy(i, j) = rechxy(i, j) + (ws.RechargeGwShallowWT * real(1.0e3));
    deeprechxy(i, j) = deeprechxy(i, j) + ws.RechargeGwDeepWT;
    smcwtdxy(i, j) = ws.SoilMoistureToWT;

    for (int k = 1; k <= nsoil; ++k) {
      smois(i, j, k - 1) = ws.SoilMoisture(k);
      sh2o(i, j, k - 1) = ws.SoilLiqWater(k);
      if (acc_etranixy.p) acc_etranixy(i, j, k - 1) = wf.TranspWatLossSoilAcc(k);
    }
    for (int iz = -nsnow + 1; iz <= 0; ++iz) {
      snicexy(i, j, iz + nsnow - 1) = ws.SnowIce(iz);
      snliqxy(i, j, iz + nsnow - 1) = ws.SnowLiqWater(iz);
    }

    // BiochemVarOutTransfer (BiochemVarOutTransferMod.F90)
    auto& bf = noahmp.biochem.flux;
    lfmassxy(i, j) = bs.LeafMass;
    rtmassxy(i, j) = bs.RootMass;
    stmassxy(i, j) = bs.StemMass;
    woodxy(i, j) = bs.WoodMass;
    stblcpxy(i, j) = bs.CarbonMassDeepSoil;
    fastcpxy(i, j) = bs.CarbonMassShallowSoil;
    if (gddxy.p) gddxy(i, j) = bs.GrowDegreeDay;
    if (pgsxy.p) pgsxy(i, j) = bs.PlantGrowStage;
    if (grainxy.p) grainxy(i, j) = bs.GrainMass;
    if (neexy.p) neexy(i, j) = bf.NetEcoExchange;
    if (gppxy.p) gppxy(i, j) = bf.GrossPriProduction;
    if (nppxy.p) nppxy(i, j) = bf.NetPriProductionTot;
    if (psnxy.p) psnxy(i, j) = bf.PhotosynTotal;
  }
};

void run_noahmp_init(NoahmpInitKernel kernel, const Range2D& range);
void run_noahmp(NoahmpKernel kernel, const Range2D& range);

}  // namespace cppwrf::phys::noahmp
