#pragma once
// Noah LSM ice surfaces, statement-by-statement ports of WRF v4.6.1:
//   SFLX_GLACIAL  phys/module_sf_noahlsm_glacial_only.F  (land ice: IVGTYP == ISICE on land)
//   SFLX_SEAICE   phys/module_sf_noah_seaice.F           (sea ice: XICE >= xice_threshold)
// Each module has its own private CSNOW, SNOW_NEW, SNOWZ0, PENMAN, HRTICE, SNOPAC, SNOWPACK
// (they differ from the module_sf_noahlsm versions); HSTEP/ROSR12 are module_sf_noahlsm's.
// Single source: device functions, column arrays in registers, literals as real(...).
#include <cmath>

#include "core/backend.hpp"
#include "core/precision.hpp"
#include "phys/lsm_noah/noah.hpp"

namespace cppwrf::phys::noah_ice {

using noah::kMaxSoilLevels;

// module_sf_noahlsm parameters used by both modules
struct IceConst {
  static constexpr real rd = real(287.04), sigma = real(5.67e-8), cph2o = real(4.218e+3),
                        cpice = real(2.106e+3), lsubf = real(3.335e+5), emissi_s = real(0.95);
  static constexpr real cp = real(7.) * real(287.) / real(2.);  // module_model_constants cp
};

// ------------------------------------------------------------------ glacial ice
namespace glacial {

CPPWRF_HD inline void csnow(real& sncond, real dsnow) {
  constexpr real unit = real(0.11631);
  real c = real(0.328) * rpow(real(10), real(2.25) * dsnow);
  sncond = real(2.0) * unit * c;
}

CPPWRF_HD inline void snow_new(real temp, real newsn, real& snowh, real& sndens) {
  real tempc = temp - real(273.15);
  real dsnew;
  if (tempc <= real(-15.)) {
    dsnew = real(0.05);
  } else {
    dsnew = real(0.05) + real(0.0017) * rpow(tempc + real(15.), real(1.5));
  }
  real snowhc = snowh * real(100.);
  real newsnc = newsn * real(100.);
  real hnewc = newsnc / dsnew;
  if (snowhc + hnewc < real(1.0e-3)) {
    sndens = dsnew > sndens ? dsnew : sndens;
  } else {
    sndens = (snowhc * sndens + hnewc * dsnew) / (snowhc + hnewc);
  }
  snowhc = snowhc + hnewc;
  snowh = snowhc * real(0.01);
}

CPPWRF_HD inline void alcalc(real alb, real snoalb, real embrd, real tsnow, real& albedo, real& emissi,
                             real dt, bool snowng, real& snotime1, real lvcoef) {
  constexpr real snacca = real(0.94), snaccb = real(0.58);
  (void)tsnow;
  albedo = alb + (snoalb - alb);
  emissi = embrd + (IceConst::emissi_s - embrd);
  real snoalb1 = snoalb + lvcoef * (real(0.85) - snoalb);
  real snoalb2 = snoalb1;
  if (snowng) {
    snotime1 = real(0.);
  } else {
    snotime1 = snotime1 + dt;
    snoalb2 = snoalb1 * rpow(snacca, rpow(snotime1 / real(86400.0), snaccb));
  }
  snoalb2 = snoalb2 > alb ? snoalb2 : alb;
  albedo = alb + (snoalb2 - alb);
  if (albedo > snoalb2) albedo = snoalb2;
}

CPPWRF_HD inline void snowz0(real& z0, real z0brd, real snowh) {
  constexpr real z0s = real(0.001);
  real burial = real(7.0) * z0brd - snowh;
  real z0eff;
  if (burial <= real(0.0007)) {
    z0eff = z0s;
  } else {
    z0eff = burial / real(7.0);
  }
  z0 = z0eff;
}

CPPWRF_HD inline void penman(real sfctmp, real sfcprs, real ch, real th2, real prcp, real fdown, real t24,
                             real ssoil, real q2, real q2sat, real& etp, real rch, real& rr, bool snowng,
                             bool frzgra, real dqsdt2, real& flx2, real emissi, real t1) {
  constexpr real elcp = real(2.4888e+3), lsubc = real(2.501000e+6), lsubs = real(2.83e+6);
  real elcp1, lvs;
  if (t1 > real(273.15)) {
    elcp1 = elcp;
    lvs = lsubc;
  } else {
    elcp1 = elcp * lsubs / lsubc;
    lvs = lsubs;
  }
  real delta = elcp1 * dqsdt2;
  real a = elcp1 * (q2sat - q2);
  rr = emissi * t24 * real(6.48e-8) / (sfcprs * ch) + real(1.0);
  if (!snowng) {
    if (prcp > real(0.0)) rr = rr + IceConst::cph2o * prcp / rch;
  } else {
    rr = rr + IceConst::cpice * prcp / rch;
  }
  if (frzgra) {
    flx2 = -IceConst::lsubf * prcp;
  } else {
    flx2 = real(0.0);
  }
  real fnet = fdown - (emissi * IceConst::sigma * t24) - ssoil - flx2;
  real rad = fnet / rch + th2 - sfctmp;
  real epsca = (a * rr + rad * delta) / (delta + rr);
  etp = epsca * rch / lvs;
}

CPPWRF_HD inline void hrtice(real* rhsts, const real* stc, real tbot, int nsoil, const real* zsoil, real yy,
                             real zz1, real df1, real* ai, real* bi, real* ci) {
  constexpr real zbot = real(-25.0);
  real hcpct = real(1.e6) * (real(0.8194) - real(0.1309) * real(0.5) * zsoil[0]);
  real df1k = df1;
  real ddz = real(1.0) / (real(-0.5) * zsoil[1]);
  ai[0] = real(0.0);
  ci[0] = (df1 * ddz) / (zsoil[0] * hcpct);
  bi[0] = -ci[0] + df1 / (real(0.5) * zsoil[0] * zsoil[0] * hcpct * zz1);
  real dtsdz = (stc[0] - stc[1]) / (real(-0.5) * zsoil[1]);
  real ssoil = df1 * (stc[0] - yy) / (real(0.5) * zsoil[0] * zz1);
  rhsts[0] = (df1 * dtsdz - ssoil) / (zsoil[0] * hcpct);
  real ddz2 = real(0.0);
  df1k = df1;
  real df1n = df1;
  for (int k = 1; k < nsoil; ++k) {
    real zmd = real(0.5) * (zsoil[k] + zsoil[k - 1]);
    hcpct = real(1.e6) * (real(0.8194) - real(0.1309) * zmd);
    df1n = real(0.32333) - (real(0.10073) * zmd);
    real dtsdz2;
    if (k != nsoil - 1) {
      real denom = real(0.5) * (zsoil[k - 1] - zsoil[k + 1]);
      dtsdz2 = (stc[k] - stc[k + 1]) / denom;
      ddz2 = real(2.) / (zsoil[k - 1] - zsoil[k + 1]);
      ci[k] = -df1n * ddz2 / ((zsoil[k - 1] - zsoil[k]) * hcpct);
    } else {
      dtsdz2 = (stc[k] - tbot) / (real(.5) * (zsoil[k - 1] + zsoil[k]) - zbot);
      ci[k] = real(0.);
    }
    real denom = (zsoil[k] - zsoil[k - 1]) * hcpct;
    rhsts[k] = (df1n * dtsdz2 - df1k * dtsdz) / denom;
    ai[k] = -df1k * ddz / ((zsoil[k - 1] - zsoil[k]) * hcpct);
    bi[k] = -(ai[k] + ci[k]);
    df1k = df1n;
    dtsdz = dtsdz2;
    ddz = ddz2;
  }
}

CPPWRF_HD inline void shflx(real* stc, int nsoil, real dt, real yy, real zz1, const real* zsoil, real tbot,
                            real df1) {
  real ai[kMaxSoilLevels], bi[kMaxSoilLevels], ci[kMaxSoilLevels], stcf[kMaxSoilLevels],
      rhsts[kMaxSoilLevels];
  hrtice(rhsts, stc, tbot, nsoil, zsoil, yy, zz1, df1, ai, bi, ci);
  noah::NoahKernel::hstep(stcf, stc, rhsts, dt, nsoil, ai, bi, ci);  // identical HSTEP/ROSR12
  for (int i = 0; i < nsoil; ++i) stc[i] = stcf[i];
}

CPPWRF_HD inline void snowpack(real sneqv, real dtsec, real& snowh, real& sndens, real tsnow, real tsoil) {
  constexpr real c1 = real(0.01), c2 = real(21.0);
  real esdc = sneqv * real(100.);
  real dthr = dtsec / real(3600.);
  real tsnowc = tsnow - real(273.15);
  real tsoilc = tsoil - real(273.15);
  real tavgc = real(0.5) * (tsnowc + tsoilc);
  real esdcx = esdc > real(1.e-2) ? esdc : real(1.e-2);
  real bfac = dthr * c1 * std::exp(real(0.08) * tavgc - c2 * sndens);
  real pexp = real(0.);
  for (int j = 4; j >= 1; --j) pexp = (real(1.) + pexp) * bfac * esdcx / real(j + 1);
  pexp = pexp + real(1.);
  real dsx = sndens * (pexp);
  if (dsx > real(0.40)) dsx = real(0.40);
  if (dsx < real(0.05)) dsx = real(0.05);
  sndens = dsx;
  if (tsnowc >= real(0.)) {
    real dw = real(0.13) * dthr / real(24.);
    sndens = sndens * (real(1.) - dw) + dw;
    if (sndens >= real(0.40)) sndens = real(0.40);
  }
  real snowhc = esdc / sndens;
  snowh = snowhc * real(0.01);
}

CPPWRF_HD inline void snopac(real& etp, real& eta, real prcp, real& prcpf, bool snowng, int nsoil, real dt,
                             real df1, real q2, real& t1, real sfctmp, real t24, real th2, real fdown,
                             real& ssoil, real* stc, real sfcprs, real rch, real rr, real& sneqv, real& sndens,
                             real& snowh, const real* zsoil, real tbot, real& snomlt, real& dew, real& flx1,
                             real& flx2, real& flx3, real& esnow, real emissi, real& ribb) {
  constexpr real esdmin = real(1.e-6), lsubc = real(2.501000e+6), lsubs = real(2.83e+6), tfreez = real(273.15);
  (void)eta; (void)q2; (void)sfcprs;
  snomlt = real(0.0);
  dew = real(0.);
  esnow = real(0.);
  real esnow1 = real(0.);
  real esnow2 = real(0.);
  real prcp1 = prcpf * real(0.001);
  (void)prcp1;
  real etanrg;
  if (etp <= real(0.0)) {
    if ((ribb >= real(0.1)) && (fdown > real(150.0))) {
      real m = etp * (real(1.0) - ribb);
      m = m < real(0.) ? m : real(0.);
      etp = (m / real(0.980) + etp * (real(0.980) - real(1.0))) / real(0.980);
    }
    real etp1 = etp * real(0.001);
    dew = -etp1;
    esnow2 = etp1 * dt;
    etanrg = etp * lsubs;
  } else {
    esnow = etp;
    esnow1 = esnow * real(0.001);
    esnow2 = esnow1 * dt;
    etanrg = esnow * lsubs;
  }
  flx1 = real(0.0);
  if (snowng) {
    flx1 = IceConst::cpice * prcp * (t1 - sfctmp);
  } else {
    if (prcp > real(0.0)) flx1 = IceConst::cph2o * prcp * (t1 - sfctmp);
  }
  real dsoil = -(real(0.5) * zsoil[0]);
  real dtot = snowh + dsoil;
  real denom = real(1.0) + df1 / (dtot * rr * rch);
  real t12a = ((fdown - flx1 - flx2 - emissi * IceConst::sigma * t24) / rch + th2 - sfctmp - etanrg / rch) / rr;
  real t12b = df1 * stc[0] / (dtot * rr * rch);
  real t12 = (sfctmp + t12a + t12b) / denom;
  if (t12 <= tfreez) {
    t1 = t12;
    ssoil = df1 * (t1 - stc[0]) / dtot;
    real d = sneqv - esnow2;
    sneqv = real(0.0) > d ? real(0.0) : d;
    flx3 = real(0.0);
    snomlt = real(0.0);
  } else {
    t1 = tfreez;
    if (dtot > real(2.0) * dsoil) dtot = real(2.0) * dsoil;
    ssoil = df1 * (t1 - stc[0]) / dtot;
    if (sneqv - esnow2 <= esdmin) {
      sneqv = real(0.0);
      snomlt = real(0.0);
      flx3 = real(0.0);
    } else {
      sneqv = sneqv - esnow2;
      real seh = rch * (t1 - th2);
      real t14 = (t1 * t1) * (t1 * t1);
      flx3 = fdown - flx1 - flx2 - emissi * IceConst::sigma * t14 - ssoil - seh - etanrg;
      if (flx3 <= real(0.0)) flx3 = real(0.0);
      real ex = flx3 * real(0.001) / IceConst::lsubf;
      snomlt = ex * dt;
      if (sneqv - snomlt >= esdmin) {
        sneqv = sneqv - snomlt;
      } else {
        ex = sneqv / dt;
        flx3 = ex * real(1000.0) * IceConst::lsubf;
        snomlt = sneqv;
        sneqv = real(0.0);
      }
    }
  }
  (void)lsubc;
  real zz1 = real(1.0);
  real yy = stc[0] - real(0.5) * ssoil * zsoil[0] * zz1 / df1;
  shflx(stc, nsoil, dt, yy, zz1, zsoil, tbot, df1);
  if (sneqv >= real(0.10)) {
    snowpack(sneqv, dt, snowh, sndens, t1, yy);
  } else {
    sneqv = real(0.10);
    snowh = real(0.50);
  }
}

// SFLX_GLACIAL. ALB in/out (unchanged), EMBRD is read only (INTENT(OUT) in WRF but never set;
// the driver's value is used by ALCALC).
struct Out {
  real eta, sheat, eta_kinematic, fdown, esnow, dew, etp, ssoil, flx1, flx2, flx3, snomlt, runoff1, q1, z0,
      albedo;
};

CPPWRF_HD inline void sflx_glacial(real ffrozp, real dt, int nsoil, const real* sldpth, real lwdn, real solnet,
                                   real sfcprs, real prcp, real sfctmp, real q2, real th2, real q2sat,
                                   real dqsdt2, real alb, real snoalb, real tbot, real z0brd, real& emissi,
                                   real embrd, real& t1, real* stc, real& snowh, real& sneqv, real ch,
                                   real& sncovr, real& snotime1, real& ribb, real lvcoef, Out& o) {
  constexpr real tfreez = real(273.15), lvh2o = real(2.501e+6), lsubs = real(2.83e+6), r = real(287.04);
  (void)lvh2o;
  real zsoil[kMaxSoilLevels];
  zsoil[0] = -sldpth[0];
  for (int kz = 1; kz < nsoil; ++kz) zsoil[kz] = -sldpth[kz] + zsoil[kz - 1];
  if (sneqv < real(0.10)) {
    sneqv = real(0.10);
    snowh = real(0.50);
  }
  real sndens = sneqv / snowh;
  real sncond;
  csnow(sncond, sndens);
  bool snowng = false, frzgra = false;
  if (prcp > real(0.0)) {
    if (ffrozp > real(0.5)) {
      snowng = true;
    } else {
      if (t1 <= tfreez) frzgra = true;
    }
  }
  real prcpf;
  if (snowng || frzgra) {
    real sn_new = prcp * dt * real(0.001);
    sneqv = sneqv + sn_new;
    prcpf = real(0.0);
    snow_new(sfctmp, sn_new, snowh, sndens);
    if (sncovr > real(0.99)) {
      if (stc[0] < (tfreez - real(5.))) sndens = real(0.2);
      if (snowng && (t1 < real(273.)) && (sfctmp < real(273.))) sndens = real(0.2);
    }
    csnow(sncond, sndens);
  } else {
    prcpf = prcp;
  }
  sncovr = real(1.0);
  alcalc(alb, snoalb, embrd, t1, o.albedo, emissi, dt, snowng, snotime1, lvcoef);
  real df1 = sncond;
  real dsoil = -(real(0.5) * zsoil[0]);
  real dtot = snowh + dsoil;
  real frcsno = snowh / dtot;
  real frcsoi = dsoil / dtot;
  df1 = frcsno * sncond + frcsoi * df1;
  if (dtot > real(2.) * dsoil) dtot = real(2.) * dsoil;
  real ssoil = df1 * (t1 - stc[0]) / dtot;
  snowz0(o.z0, z0brd, snowh);
  real fdown = solnet + lwdn;
  real t2v = sfctmp * (real(1.0) + real(0.61) * q2);
  real rho = sfcprs / (IceConst::rd * t2v);
  real rch = rho * real(1004.6) * ch;
  real t24 = sfctmp * sfctmp * sfctmp * sfctmp;
  real etp, rr, flx2;
  penman(sfctmp, sfcprs, ch, th2, prcp, fdown, t24, ssoil, q2, q2sat, etp, rch, rr, snowng, frzgra, dqsdt2, flx2,
         emissi, t1);
  real eta = real(0), snomlt, dew, flx1, flx3, esnow;
  snopac(etp, eta, prcp, prcpf, snowng, nsoil, dt, df1, q2, t1, sfctmp, t24, th2, fdown, ssoil, stc, sfcprs, rch,
         rr, sneqv, sndens, snowh, zsoil, tbot, snomlt, dew, flx1, flx2, flx3, esnow, emissi, ribb);
  real eta_kinematic = etp;
  o.q1 = q2 + eta_kinematic * IceConst::cp / rch;
  o.sheat = -(ch * IceConst::cp * sfcprs) / (r * t2v) * (th2 - t1);
  esnow = esnow * lsubs;
  etp = etp * lsubs;
  if (etp > real(0.)) {
    eta = esnow;
  } else {
    eta = etp;
  }
  o.ssoil = real(-1.0) * ssoil;
  o.runoff1 = snomlt / dt;
  o.eta = eta;
  o.eta_kinematic = eta_kinematic;
  o.fdown = fdown;
  o.esnow = esnow;
  o.dew = dew;
  o.etp = etp;
  o.flx1 = flx1;
  o.flx2 = flx2;
  o.flx3 = flx3;
  o.snomlt = snomlt;
}

}  // namespace glacial

// ------------------------------------------------------------------ sea ice
namespace seaice {

constexpr real tfreez = real(273.15);

CPPWRF_HD inline void csnow(real& sncond, real dsnow) {
  constexpr real unit = real(0.11631);
  real c = real(0.328) * rpow(real(10), real(2.25) * dsnow);
  sncond = real(2.0) * unit * c;
}

CPPWRF_HD inline void snow_new(real temp, real newsn, real& snowh, real& sndens) {
  real snowhc = snowh * real(100.);
  real newsnc = newsn * real(100.);
  real tempc = temp - real(273.15);
  real dsnew;
  if (tempc <= real(-15.)) {
    dsnew = real(0.05);
  } else {
    dsnew = real(0.05) + real(0.0017) * rpow(tempc + real(15.), real(1.5));
  }
  real hnewc = newsnc / dsnew;
  if (snowhc + hnewc < real(1.0e-3)) {
    sndens = dsnew > sndens ? dsnew : sndens;
  } else {
    sndens = (snowhc * sndens + hnewc * dsnew) / (snowhc + hnewc);
  }
  snowhc = snowhc + hnewc;
  snowh = snowhc * real(0.01);
}

CPPWRF_HD inline void snowz0(real sncovr, real& z0, real z0brd, real snowh) {
  constexpr real z0s = real(0.001);
  real burial = real(7.0) * z0brd - snowh;
  real z0eff;
  if (burial <= real(0.0007)) {
    z0eff = z0s;
  } else {
    z0eff = burial / real(7.0);
  }
  z0 = (real(1.) - sncovr) * z0brd + sncovr * z0eff;
}

CPPWRF_HD inline void penman(real sfctmp, real sfcprs, real ch, real th2, real prcp, real fdown, real t24,
                             real ssoil, real q2, real q2sat, real& etp, real rch, real& rr, bool snowng,
                             bool frzgra, real dqsdt2, real& flx2, real emissi, real t1) {
  constexpr real elcp = real(2.4888e+3), lsubc = real(2.501000e+6), lsubs = real(2.83e+6);
  real elcp1, lvs;
  if (t1 > real(273.15)) {
    elcp1 = elcp;
    lvs = lsubc;
  } else {
    elcp1 = elcp * lsubs / lsubc;
    lvs = lsubs;
  }
  flx2 = real(0.0);
  real delta = elcp1 * dqsdt2;
  rr = emissi * t24 * real(6.48e-8) / (sfcprs * ch) + real(1.0);
  if (prcp > real(0.0)) {
    if (!snowng) {
      rr = rr + IceConst::cph2o * prcp / rch;
    } else {
      rr = rr + IceConst::cpice * prcp / rch;
    }
  }
  real fnet = fdown - emissi * IceConst::sigma * t24 - ssoil;
  if (frzgra) {
    flx2 = -IceConst::lsubf * prcp;
    fnet = fnet - flx2;
  }
  real rad = fnet / rch + th2 - sfctmp;
  real a = elcp1 * (q2sat - q2);
  real epsca = (a * rr + rad * delta) / (delta + rr);
  etp = epsca * rch / lvs;
}

CPPWRF_HD inline void hrtice(real* rhsts, const real* stc, real tbot, int nsoil, const real* zsoil, real yy,
                             real zz1, real df1, real* ai, real* bi, real* ci) {
  real hcpct = real(1.72396e+6);
  real zbot = zsoil[nsoil - 1];
  real ddz = real(1.0) / (real(-0.5) * zsoil[1]);
  ai[0] = real(0.0);
  ci[0] = (df1 * ddz) / (zsoil[0] * hcpct);
  bi[0] = -ci[0] + df1 / (real(0.5) * zsoil[0] * zsoil[0] * hcpct * zz1);
  real dtsdz = (stc[0] - stc[1]) / (real(-0.5) * zsoil[1]);
  real ssoil = df1 * (stc[0] - yy) / (real(0.5) * zsoil[0] * zz1);
  rhsts[0] = (df1 * dtsdz - ssoil) / (zsoil[0] * hcpct);
  real ddz2 = real(0.0);
  real df1k = df1;
  real df1n = df1;
  for (int k = 1; k < nsoil; ++k) {
    real dtsdz2;
    if (k != nsoil - 1) {
      real denom = real(0.5) * (zsoil[k - 1] - zsoil[k + 1]);
      dtsdz2 = (stc[k] - stc[k + 1]) / denom;
      ddz2 = real(2.) / (zsoil[k - 1] - zsoil[k + 1]);
      ci[k] = -df1n * ddz2 / ((zsoil[k - 1] - zsoil[k]) * hcpct);
    } else {
      dtsdz2 = (stc[k] - tbot) / (real(.5) * (zsoil[k - 1] + zsoil[k]) - zbot);
      ci[k] = real(0.);
    }
    real denom = (zsoil[k] - zsoil[k - 1]) * hcpct;
    rhsts[k] = (df1n * dtsdz2 - df1k * dtsdz) / denom;
    ai[k] = -df1k * ddz / ((zsoil[k - 1] - zsoil[k]) * hcpct);
    bi[k] = -(ai[k] + ci[k]);
    df1k = df1n;
    dtsdz = dtsdz2;
    ddz = ddz2;
  }
}

CPPWRF_HD inline void shflx(real* stc, int nsoil, real dt, real yy, real zz1, const real* zsoil, real tbot,
                            real df1) {
  real ai[kMaxSoilLevels], bi[kMaxSoilLevels], ci[kMaxSoilLevels], stcf[kMaxSoilLevels],
      rhsts[kMaxSoilLevels];
  hrtice(rhsts, stc, tbot, nsoil, zsoil, yy, zz1, df1, ai, bi, ci);
  noah::NoahKernel::hstep(stcf, stc, rhsts, dt, nsoil, ai, bi, ci);  // module_sf_noahlsm HSTEP
  for (int i = 0; i < nsoil; ++i) stc[i] = stcf[i];
}

CPPWRF_HD inline void snowpack(real esd, real dtsec, real& snowh, real& sndens, real tsnow, real tsoil) {
  constexpr real c1 = real(0.01), c2 = real(21.0);
  real esdc = esd * real(100.);
  real dthr = dtsec / real(3600.);
  real tsnowc = tsnow - real(273.15);
  real tsoilc = tsoil - real(273.15);
  real tavgc = real(0.5) * (tsnowc + tsoilc);
  real esdcx = esdc > real(1.e-2) ? esdc : real(1.e-2);
  real bfac = dthr * c1 * std::exp(real(0.08) * tavgc - c2 * sndens);
  real pexp = real(0.);
  for (int j = 4; j >= 1; --j) pexp = (real(1.) + pexp) * bfac * esdcx / real(j + 1);
  pexp = pexp + real(1.);
  real dsx = sndens * (pexp);
  if (dsx > real(0.40)) dsx = real(0.40);
  if (dsx < real(0.05)) dsx = real(0.05);
  sndens = dsx;
  if (tsnowc >= real(0.)) {
    real dw = real(0.13) * dthr / real(24.);
    sndens = sndens * (real(1.) - dw) + dw;
    if (sndens >= real(0.40)) sndens = real(0.40);
  }
  real snowhc = esdc / sndens;
  snowh = snowhc * real(0.01);
}

CPPWRF_HD inline void snopac(real& etp, real prcp, bool snowng, int nsoil, real dt, real df1, real& t1,
                             real sfctmp, real t24, real th2, real fdown, real& ssoil, real* stc, real rch, real rr,
                             real& sncovr, real& esd, real& sndens, real& snowh, const real* zsoil, real tbot,
                             real& snomlt, real& dew, real& flx1, real& flx2, real& flx3, real& esnow,
                             real emissi, real& ribb, int seaice_albedo_opt) {
  constexpr real esdmin = real(1.e-6), lsubc = real(2.501000e+6), lsubs = real(2.83e+6);
  dew = real(0.);
  esnow = real(0.);
  real esnow1 = real(0.);
  real esnow2 = real(0.);
  real etanrg;
  if (etp <= real(0.0)) {
    if ((ribb >= real(0.1)) && (fdown > real(150.0))) {
      real m = etp * (real(1.0) - ribb);
      m = m < real(0.) ? m : real(0.);
      etp = (m * sncovr / real(0.980) + etp * (real(0.980) - sncovr)) / real(0.980);
    }
    real etp1 = etp * real(0.001);
    dew = -etp1;
    esnow2 = etp1 * dt;
    etanrg = etp * ((real(1.) - sncovr) * lsubc + sncovr * lsubs);
  } else {
    esnow = etp;
    esnow1 = esnow * real(0.001);
    esnow2 = esnow1 * dt;
    etanrg = esnow * lsubs;
    esnow = etp * sncovr;
    esnow1 = esnow * real(0.001);
    esnow2 = esnow1 * dt;
    etanrg = esnow * lsubs;
  }
  flx1 = real(0.0);
  if (snowng) {
    flx1 = IceConst::cpice * prcp * (t1 - sfctmp);
  } else {
    if (prcp > real(0.0)) flx1 = IceConst::cph2o * prcp * (t1 - sfctmp);
  }
  real dsoil = -(real(0.5) * zsoil[0]);
  real dtot = snowh + dsoil;
  real denom = real(1.0) + df1 / (dtot * rr * rch);
  real t12a = ((fdown - flx1 - flx2 - emissi * IceConst::sigma * t24) / rch + th2 - sfctmp - etanrg / rch) / rr;
  real t12b = df1 * stc[0] / (dtot * rr * rch);
  real t12 = (sfctmp + t12a + t12b) / denom;
  if (t12 <= tfreez) {
    t1 = t12;
    ssoil = df1 * (t1 - stc[0]) / dtot;
    real d = esd - esnow2;
    esd = real(0.0) > d ? real(0.0) : d;
    flx3 = real(0.0);
    snomlt = real(0.0);
  } else {
    t1 = tfreez;
    ssoil = df1 * (t1 - stc[0]) / dtot;
    if (esd - esnow2 <= esdmin) {
      esd = real(0.0);
      snomlt = real(0.0);
      flx3 = real(0.0);
    } else {
      esd = esd - esnow2;
      real seh = rch * (t1 - th2);
      real t14 = (t1 * t1) * (t1 * t1);
      flx3 = fdown - flx1 - flx2 - emissi * IceConst::sigma * t14 - ssoil - seh - etanrg;
      if (flx3 <= real(0.0)) flx3 = real(0.0);
      real ex = flx3 * real(0.001) / IceConst::lsubf;
      snomlt = ex * dt;
      if (esd - snomlt >= esdmin) {
        esd = esd - snomlt;
      } else {
        ex = esd / dt;
        flx3 = ex * real(1000.0) * IceConst::lsubf;
        snomlt = esd;
        esd = real(0.0);
      }
    }
  }
  real zz1 = real(1.0);
  real yy = stc[0] - real(0.5) * ssoil * zsoil[0] * zz1 / df1;
  shflx(stc, nsoil, dt, yy, zz1, zsoil, tbot, df1);
  if (seaice_albedo_opt == 1) {
    if (esd >= real(0.0001)) {
      snowpack(esd, dt, snowh, sndens, t1, yy);
    } else {
      esd = real(0.0001);
      snowh = real(0.0005);
      sncovr = real(0.005);
    }
  } else {
    if (esd >= real(0.01)) {
      snowpack(esd, dt, snowh, sndens, t1, yy);
    } else {
      esd = real(0.01);
      snowh = real(0.05);
      sncovr = real(1.0);
    }
  }
}

struct Opts {
  int albedo_opt = 0;
  real albedo_default = real(0.65);
  int snowdepth_opt = 0;
  real snowdepth_max = real(1.e10);
  real snowdepth_min = real(0.001);
};

struct Out {
  real albedo, eta, sheat, eta_kinematic, fdown, esnow, dew, etp, ssoil, flx1, flx2, flx3, snomlt, runoff1, q1,
      z0;
};

// SFLX_SEAICE. Returns false where WRF stops with "Physical snow depth is less than snow water equiv."
CPPWRF_HD inline bool sflx_seaice(const Opts& op, real ffrozp, real dt, int nsoil, real sithick, real lwdn,
                                  real solnet, real sfcprs, real prcp, real sfctmp, real q2, real th2, real q2sat,
                                  real dqsdt2, real snoalb, real tbot, real& z0brd, real& emissi, real& t1,
                                  real* stc, real& snowh, real& sneqv, real ch, real albedosi, real snowonsi,
                                  real& sncovr, real& ribb, Out& o) {
  constexpr real lvh2o = real(2.501e+6), lsubs = real(2.83e+6), r = real(287.04);
  (void)snoalb;
  o.runoff1 = real(0.0);
  o.snomlt = real(0.0);
  real zsoil[kMaxSoilLevels];
  for (int kz = 0; kz < nsoil; ++kz) zsoil[kz] = -sithick * real(kz + 1) / real(nsoil);
  z0brd = real(0.001);
  bool snowng = false, frzgra = false;
  if (op.albedo_opt == 1) {
    if (sneqv < real(0.0001)) {
      sneqv = real(0.0001);
      snowh = real(0.0005);
    }
  } else {
    if (sneqv < real(0.01)) {
      sneqv = real(0.01);
      snowh = real(0.05);
    }
  }
  real sndens = real(0);
  if (op.snowdepth_opt == 0) {
    sndens = sneqv / snowh;
    real m = snowh < op.snowdepth_max ? snowh : op.snowdepth_max;
    snowh = op.snowdepth_min > m ? op.snowdepth_min : m;
    sneqv = snowh * sndens;
  } else if (op.snowdepth_opt == 1) {
    sndens = real(0.3);
    snowh = snowonsi;
    sneqv = snowh * sndens;
  }
  sndens = sneqv / snowh;
  if (sndens > real(1.0)) return false;
  real sncond;
  csnow(sncond, sndens);
  if (prcp > real(0.0)) {
    if (ffrozp > real(0.5)) {
      snowng = true;
    } else {
      if (t1 <= tfreez) frzgra = true;
    }
  }
  if (snowng || frzgra) {
    real sn_new = prcp * dt * real(0.001);
    sneqv = sneqv + sn_new;
    snow_new(sfctmp, sn_new, snowh, sndens);
    if (sncovr > real(0.99)) {
      if (stc[0] < (tfreez - real(5.))) sndens = real(0.2);
      if (snowng && (t1 < real(273.)) && (sfctmp < real(273.))) sndens = real(0.2);
    }
    csnow(sncond, sndens);
  }
  real albedo;
  if (op.albedo_opt == 1) {
    real alb_snow, alb_ice;
    if (t1 < real(268.15)) {
      alb_snow = real(0.8);
    } else if ((t1 >= real(268.15)) && (t1 < real(273.15))) {
      alb_snow = real(0.65) - (real(0.03) * (t1 - real(273.15)));
    } else {
      alb_snow = real(0.65);
    }
    if (sfctmp <= real(273.15)) {
      alb_ice = real(0.65);
    } else if ((sfctmp > real(273.15)) && (sfctmp < real(278.15))) {
      alb_ice = real(0.65) - (real(0.04) * (sfctmp - real(273.15)));
    } else {
      alb_ice = real(0.45);
    }
    real z0n = real(0.10);
    real sncovrr = snowh / (snowh + z0n);
    albedo = (sncovrr * alb_snow) + ((real(1.0) - sncovrr) * alb_ice);
  } else if (op.albedo_opt == 2) {
    sncovr = real(1.0);
    emissi = real(0.98);
    albedo = albedosi;
  } else {
    sncovr = real(1.0);
    emissi = real(0.98);
    albedo = op.albedo_default;
  }
  real df1 = real(2.2);
  real dsoil = -(real(0.5) * zsoil[0]);
  real dtot = snowh + dsoil;
  real frcsno = snowh / dtot;
  real frcsoi = dsoil / dtot;
  real df1a = frcsno * sncond + frcsoi * df1;
  df1 = df1a * sncovr + df1 * (real(1.0) - sncovr);
  real ssoil = df1 * (t1 - stc[0]) / dtot;
  snowz0(sncovr, o.z0, z0brd, snowh);
  real fdown = solnet + lwdn;
  real t2v = sfctmp * (real(1.0) + real(0.61) * q2);
  real t24 = sfctmp * sfctmp * sfctmp * sfctmp;
  real rho = sfcprs / (IceConst::rd * t2v);
  real rch = rho * real(1004.6) * ch;
  real etp, rr, flx2;
  penman(sfctmp, sfcprs, ch, th2, prcp, fdown, t24, ssoil, q2, q2sat, etp, rch, rr, snowng, frzgra, dqsdt2, flx2,
         emissi, t1);
  real esnow = real(0.0), snomlt, dew, flx1, flx3;
  snopac(etp, prcp, snowng, nsoil, dt, df1, t1, sfctmp, t24, th2, fdown, ssoil, stc, rch, rr, sncovr, sneqv,
         sndens, snowh, zsoil, tbot, snomlt, dew, flx1, flx2, flx3, esnow, emissi, ribb, op.albedo_opt);
  real eta_kinematic = etp;
  if (op.snowdepth_opt == 0) {
    sndens = sneqv / snowh;
    real m = snowh < op.snowdepth_max ? snowh : op.snowdepth_max;
    snowh = op.snowdepth_min > m ? op.snowdepth_min : m;
    sneqv = snowh * sndens;
  } else if (op.snowdepth_opt == 1) {
    sndens = real(0.3);
    snowh = snowonsi;
    sneqv = snowh * sndens;
  }
  o.q1 = q2 + eta_kinematic * IceConst::cp / rch;
  o.sheat = -(ch * IceConst::cp * sfcprs) / (r * t2v) * (th2 - t1);
  esnow = esnow * lsubs;
  etp = etp * ((real(1.) - sncovr) * lvh2o + sncovr * lsubs);
  real eta;
  if (etp > real(0.)) {
    eta = esnow;
  } else {
    eta = etp;
  }
  o.ssoil = real(-1.0) * ssoil;
  o.runoff1 = snomlt / dt;
  o.albedo = albedo;
  o.eta = eta;
  o.eta_kinematic = eta_kinematic;
  o.fdown = fdown;
  o.esnow = esnow;
  o.dew = dew;
  o.etp = etp;
  o.flx1 = flx1;
  o.flx2 = flx2;
  o.flx3 = flx3;
  o.snomlt = snomlt;
  return true;
}

}  // namespace seaice
}  // namespace cppwrf::phys::noah_ice
