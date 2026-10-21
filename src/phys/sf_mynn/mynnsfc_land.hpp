#pragma once
// module_sf_mynnsfc_land.F90 (WRF v4.8.0 phys/MYNN-SFC), subroutine mynnsfc_land, for WRF's call from
// mynnsfc_driver: flag_iter = .true., spp_sfc = 0, psi_opt = 0, compute_flux = compute_diag = .true.,
// sf_mynn_sfcflux_land = 0 (Zilitinkevich 1995, czil = 0.095), lsm /= ruclsm. Other options are rejected
// by the configuration.
#include "phys/sf_mynn/mynnsfc_common.hpp"

namespace cppwrf::phys::mynnsfc {

// zilitinkevich_1995 with zt_opt = 0
CPPWRF_HD inline void zilitinkevich_1995(real z_0, real& zt, real& zq, real restar) {
  const real czil = real(0.095);
  const real znt = z_0;
  zt = znt * std::exp(-karman * czil * std::sqrt(restar));
  zt = c_min(zt, p75 * znt);
  zt = c_max(zt, real(0.0001));
  zq = zt;
}

// Point state of mynnsfc_land (the INTENT(INOUT) / OUT scalars of the call in mynnsfc_driver)
struct SfcPoint {
  real chs, chs2, cqs2, cqs, ust, ustm, stress, qsfc, rmol, zol, mol, psim, psih, hfx, qfx, u10, v10, th2, t2, q2,
      flhc, flqc, lh, gz1oz0, wspd, rb, cpm, ch, cm, wstar, qstar, znt, qgh, psix, psit, psix10, psit2, ck, cka, cd,
      cda;
};

struct SfcIn {
  int itimestep;
  bool flag_restart;
  real dx, u_1, v_1, t_1, qv_1, p_1, dz8w_1, rho_1, u_2, v_2, dz8w_2, pblh, psfcpa, mavail, tskin, tsurf, snowh;
  int isfflx;
};

CPPWRF_HD inline void mynnsfc_land(const SfcIn& in, SfcPoint& o, const PsiTables& tb) {
  constexpr real prt = real(1.);
  real zol = o.zol;
  const real psfc = in.psfcpa / real(1000.);
  real qsfcmr;
  (void)qsfcmr;
  // flag_iter
  {
    const real tabs = (in.itimestep == 1) ? p99 * in.tskin + p01 * in.t_1 : in.tskin;
    const real e1 = svp1 * std::exp(svp2 * (tabs - svpt0) / (tabs - svp3));
    o.qsfc = ep2 * e1 / (psfc - ep3 * e1);
    qsfcmr = o.qsfc / (one - o.qsfc);
  }
  {
    const real e1 = svp1 * std::exp(svp2 * (in.t_1 - svpt0) / (in.t_1 - svp3));
    const real pl = in.p_1 / real(1000.);
    o.qgh = ep2 * e1 / (pl - e1);
  }
  const real thcon = std::pow(real(100000.) / in.psfcpa, rovcp);
  const real tsk = in.tskin;
  const real thsk = tsk * thcon;
  const real thvsk = thsk * (one + ep1 * o.qsfc);

  const real th_1 = in.t_1 * std::pow(real(100000.) / in.p_1, rovcp);
  const real tc_1 = in.t_1 - real(273.15);
  const real thv_1 = th_1 * (one + ep1 * in.qv_1);
  const real tv_1 = in.t_1 * (one + ep1 * in.qv_1);
  const real rho_1 = in.p_1 / (rd * tv_1);
  const real za = p5 * in.dz8w_1;
  const real za2 = in.dz8w_1 + p5 * in.dz8w_2;
  const real govrth = grav / th_1;
  o.cpm = cp * (one + real(0.84) * in.qv_1);

  o.wspd = std::sqrt(in.u_1 * in.u_1 + in.v_1 * in.v_1);
  const real dthvdz = (thv_1 - thvsk);
  const real fluxc = c_max(o.hfx / rho_1 / cp + ep1 * thvsk * o.qfx / rho_1, zero);
  o.wstar = vconvc * std::pow(grav / tsk * c_min(real(1.5) * in.pblh, real(4000.)) * fluxc, p333);
  const real vsgd = c_min(real(0.32) * std::pow(c_max(in.dx / real(5000.) - one, zero), p333), p5);
  o.wspd = std::sqrt(o.wspd * o.wspd + o.wstar * o.wstar + vsgd * vsgd);
  o.wspd = c_max(o.wspd, wmin);
  o.rb = govrth * za * dthvdz / (o.wspd * o.wspd);
  o.rb = c_max(o.rb, real(-2.0));
  o.rb = c_min(o.rb, two);

  const real visc = real(1.326e-5) * (one + real(6.542e-3) * tc_1 + real(8.301e-6) * tc_1 * tc_1 -
                                      real(4.84e-9) * tc_1 * tc_1 * tc_1);
  real zntstoch = o.znt;
  zntstoch = c_min(zntstoch, in.dz8w_1 * real(0.0666));
  const real restar = c_max(o.ust * zntstoch / visc, p1);
  real zt, zq;
  zilitinkevich_1995(zntstoch, zt, zq, restar);

  o.gz1oz0 = std::log((za + zntstoch) / zntstoch);
  const real gz1ozt = std::log((za + zntstoch) / zt);
  const real gz2oz0 = std::log((two + zntstoch) / zntstoch);
  (void)gz2oz0;
  const real gz2ozt = std::log((two + zntstoch) / zt);
  const real gz10oz0 = std::log((ten + zntstoch) / zntstoch);
  const real gz10ozt = std::log((ten + zntstoch) / zt);
  const real zratio = zntstoch / zt;

  real psim10, psih10, psih2;
  if (o.rb > zero) {
    if (!in.flag_restart || (in.flag_restart && in.itimestep > 1)) {
      li_etal_2010(zol, o.rb, za / zntstoch, zratio);
      zol = c_max(zol, zero);
      zol = c_min(zol, twenty);
      zol = zolrib(o.rb, za, zntstoch, zt, o.gz1oz0, gz1ozt, zol, tb);
    }
    zol = c_max(zol, zero);
    zol = c_min(zol, twenty);
    const real zolzt = zol * zt / za;
    const real zolz0 = zol * zntstoch / za;
    const real zolza = zol * (za + zntstoch) / za;
    const real zol10 = zol * (ten + zntstoch) / za;
    const real zol2 = zol * (two + zntstoch) / za;
    o.psim = psim_stable(zolza, tb) - psim_stable(zolz0, tb);
    o.psih = psih_stable(zolza, tb) - psih_stable(zolzt, tb);
    psim10 = psim_stable(zol10, tb) - psim_stable(zolz0, tb);
    psih10 = psih_stable(zol10, tb) - psih_stable(zolz0, tb);
    psih2 = psih_stable(zol2, tb) - psih_stable(zolzt, tb);
    o.rmol = zol / za;
  } else if (o.rb == zero) {
    o.psim = zero;
    o.psih = o.psim;
    psim10 = zero;
    psih10 = zero;
    psih2 = zero;
    zol = zero;
    o.rmol = zero;
  } else {
    if (!in.flag_restart || (in.flag_restart && in.itimestep > 1)) {
      li_etal_2010(zol, o.rb, za / zntstoch, zratio);
      zol = c_max(zol, real(-20.0));
      zol = c_min(zol, zero);
      zol = zolrib(o.rb, za, zntstoch, zt, o.gz1oz0, gz1ozt, zol, tb);
    }
    zol = c_max(zol, real(-20.0));
    zol = c_min(zol, zero);
    const real zolzt = zol * zt / za;
    const real zolz0 = zol * zntstoch / za;
    const real zolza = zol * (za + zntstoch) / za;
    const real zol10 = zol * (ten + zntstoch) / za;
    const real zol2 = zol * (two + zntstoch) / za;
    o.psim = psim_unstable(zolza, tb) - psim_unstable(zolz0, tb);
    o.psih = psih_unstable(zolza, tb) - psih_unstable(zolzt, tb);
    psim10 = psim_unstable(zol10, tb) - psim_unstable(zolz0, tb);
    psih10 = psih_unstable(zol10, tb) - psih_unstable(zolz0, tb);
    psih2 = psih_unstable(zol2, tb) - psih_unstable(zolzt, tb);
    o.psih = c_min(o.psih, p9 * gz1ozt);
    o.psim = c_min(o.psim, p9 * o.gz1oz0);
    psih2 = c_min(psih2, p9 * gz2ozt);
    psim10 = c_min(psim10, p9 * gz10oz0);
    psih10 = c_min(psih10, p9 * gz10ozt);
    o.rmol = zol / za;
  }
  o.zol = zol;
  o.psix = c_max(o.gz1oz0 - o.psim, one);
  o.psix10 = c_max(gz10oz0 - psim10, one);
  o.psit = c_max(gz1ozt - o.psih, one);
  o.psit2 = c_max(gz2ozt - psih2, one);
  const real psiq = c_max(std::log((za + zq) / zq) - o.psih, one);
  const real psiq2 = c_max(std::log((two + zq) / zq) - psih2, one);
  const real psiq10 = c_max(std::log((ten + zq) / zq) - psih10, one);

  o.ust = p5 * (o.ust + karman * o.wspd / o.psix);
  o.ust = c_max(o.ust, real(0.005));
  o.stress = o.ust * o.ust;
  o.ustm = o.ust;
  real dtg = thv_1 - thvsk;
  o.mol = karman * dtg / o.psit / prt;
  const real dqg = (in.qv_1 - o.qsfc) * real(1000.);
  o.qstar = karman * dqg / psiq / prt;

  if (in.isfflx < 1) {
    o.qfx = zero; o.hfx = zero; o.flhc = zero; o.flqc = zero; o.lh = zero; o.chs = zero; o.ch = zero;
    o.chs2 = zero; o.cqs = zero; o.cqs2 = zero; o.cm = zero;
    o.ck = zero; o.cd = zero; o.cka = zero; o.cda = zero;
  } else {
    o.flqc = rho_1 * in.mavail * o.ust * karman / psiq;
    o.flhc = rho_1 * o.cpm * o.ust * karman / o.psit;
    o.qfx = o.flqc * (o.qsfc - in.qv_1);
    o.qfx = c_max(o.qfx, real(-0.02));
    o.lh = xlv * o.qfx;
    o.hfx = rho_1 * o.cpm * karman * o.wspd / o.psix * karman / o.psit * (thsk - th_1);
    o.hfx = c_max(o.hfx, real(-250.));
    o.chs = o.ust * karman / o.psit;
    o.cqs = o.ust * karman / psiq;
    o.cqs2 = o.ust * karman / psiq2;
    o.chs2 = o.ust * karman / o.psit2;
    o.ck = (karman / o.psix10) * (karman / psiq10);
    o.cd = (karman / o.psix10) * (karman / o.psix10);
    o.cka = (karman / o.psix) * (karman / psiq);
    o.cda = (karman / o.psix) * (karman / o.psix);
    o.ch = (karman / o.psix) * (karman / o.psit);
    o.cm = (karman / o.psix) * (karman / o.psix);
  }

  // compute_diag
  if (za <= real(7.0)) {
    if (za2 > real(7.0) && za2 < real(13.0)) {
      o.u10 = in.u_2;
      o.v10 = in.v_2;
    } else {
      o.u10 = in.u_1 * std::log(ten / zntstoch) / std::log(za / zntstoch);
      o.v10 = in.v_1 * std::log(ten / zntstoch) / std::log(za / zntstoch);
    }
  } else if (za > real(7.0) && za < real(13.0)) {
    o.u10 = in.u_1 * std::log(ten / zntstoch) / std::log(za / zntstoch);
    o.v10 = in.v_1 * std::log(ten / zntstoch) / std::log(za / zntstoch);
  } else {
    o.u10 = in.u_1 * o.psix10 / o.psix;
    o.v10 = in.v_1 * o.psix10 / o.psix;
  }
  dtg = th_1 - thsk;
  o.th2 = thsk + dtg * o.psit2 / o.psit;
  if ((th_1 > thsk && (o.th2 < thsk || o.th2 > th_1)) || (th_1 < thsk && (o.th2 > thsk || o.th2 < th_1)))
    o.th2 = thsk + two * (th_1 - thsk) / za;
  o.t2 = o.th2 * std::pow(in.psfcpa / real(100000.), rovcp);
  o.q2 = o.qsfc + (in.qv_1 - o.qsfc) * psiq2 / psiq;
  o.q2 = c_max(o.q2, c_min(o.qsfc, in.qv_1));
  o.q2 = c_min(o.q2, real(1.05) * in.qv_1);
}

}  // namespace cppwrf::phys::mynnsfc
