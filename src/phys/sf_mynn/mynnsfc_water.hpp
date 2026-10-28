#pragma once
// module_sf_mynnsfc_water.F90 (WRF v4.8.0 phys/MYNN-SFC), subroutine mynnsfc_water, for WRF's call from
// mynnsfc_driver: flag_iter = .true., spp_sfc = 0, psi_opt = 0, sfc_z0_type = 0, compute_flux = compute_diag =
// .true., shalwater_z0 = 0, sf_mynn_sfcflux_water 0 (COARE 3.0), 1 (COARE 3.5: Edson et al. 2013 z0, Fairall et
// al. 2014 zt / zq), 2 / 3 (Davis et al. 2008 z0; Fairall 2014 / Garratt 1992 zt) or 4 (Taylor and Yelland
// 2001 z0); 5 (GFS) is rejected by the configuration.
#include "phys/sf_mynn/mynnsfc_land.hpp"

namespace cppwrf::phys::mynnsfc {

inline constexpr real g_inv = real(1.) / grav;

CPPWRF_HD inline void edson_etal_2013(real& z_0, real ustar, real wsp10, real visc, real zu) {
  constexpr real m = real(0.0017), b = real(-0.005);
  real wsp10m = wsp10 * std::log(ten / real(1e-4)) / std::log(zu / real(1e-4));
  wsp10m = c_min(real(19.), wsp10m);
  real czc = m * wsp10m + b;
  czc = c_max(czc, zero);
  z_0 = czc * ustar * ustar * g_inv + (real(0.11) * visc / c_max(ustar, real(0.05)));
  z_0 = c_max(z_0, real(1.27e-7));
  z_0 = c_min(z_0, real(2.85e-3));
}

CPPWRF_HD inline void coare30_znt(real& z_0, real ustar, real wsp10, real visc, real zu) {
  constexpr real czo2 = real(0.011);
  const real wsp10m = wsp10 * std::log(ten / real(1e-4)) / std::log(zu / real(1e-4));
  const real czc = czo2 + real(0.007) * c_min(c_max((wsp10m - ten) / real(8.), zero), one);
  z_0 = czc * ustar * ustar * g_inv + (real(0.11) * visc / c_max(ustar, real(0.05)));
  z_0 = c_max(z_0, real(1.27e-7));
  z_0 = c_min(z_0, real(2.85e-3));
}

CPPWRF_HD inline void davis_etal_2008(real& z_0, real ustar) {
  constexpr real ozo = real(1.59e-5);
  const real zw = c_min(std::pow(ustar / real(1.06), real(0.3)), one);
  const real zn1 = real(0.011) * ustar * ustar * g_inv + ozo;
  const real zn2 = ten * std::exp(-(real(9.5) * std::pow(ustar, -p333))) +
                   real(0.11) * real(1.5e-5) / c_max(ustar, real(0.01));
  z_0 = (one - zw) * zn1 + zw * zn2;
  z_0 = c_max(z_0, real(1.27e-7));
  z_0 = c_min(z_0, real(2.85e-3));
}

CPPWRF_HD inline void taylor_yelland_2001(real& z_0, real ustar, real wsp10) {
  (void)ustar;
  constexpr real pi = real(3.14159265);
  const real hs = real(0.0248) * (wsp10 * wsp10);
  const real tp = real(0.729) * c_max(wsp10, real(0.1));
  const real lp = grav * (tp * tp) / (two * pi);
  z_0 = real(1200.) * hs * std::pow(hs / lp, real(4.5));
  z_0 = c_max(z_0, real(1.27e-7));
  z_0 = c_min(z_0, real(2.85e-3));
}

// garratt_1992 of the water module with landsea = 2 (the water branch)
CPPWRF_HD inline void garratt_1992_water(real& zt, real& zq, real z_0, real ren) {
  zt = z_0 * std::exp(two - (real(2.48) * std::pow(ren, real(0.25))));
  zq = z_0 * std::exp(two - (real(2.28) * std::pow(ren, real(0.25))));
  zq = c_min(zq, real(5.5e-5));
  zq = c_max(zq, real(2.0e-9));
  zt = c_min(zt, real(5.5e-5));
  zt = c_max(zt, real(2.0e-9));
}

CPPWRF_HD inline void fairall_etal_2003(real& zt, real& zq, real ren) {
  zt = (real(5.5e-5)) * std::pow(ren, real(-0.60));  // both branches of WRF's ren <= 2 test are the same
  zq = zt;
  zt = c_min(zt, real(1.0e-4));
  zt = c_max(zt, real(2.0e-9));
  zq = c_min(zq, real(1.0e-4));
  zq = c_max(zq, real(2.0e-9));
}

CPPWRF_HD inline void fairall_etal_2014(real& zt, real& zq, real ren) {
  zt = c_min(real(1.6e-4), real(5.8e-5) / std::pow(ren, real(0.72)));
  zq = zt;
  zt = c_min(zt, real(1.0e-4));
  zt = c_max(zt, real(2.0e-9));
  zq = c_min(zq, real(1.0e-4));
  zq = c_max(zq, real(2.0e-9));
}

CPPWRF_HD inline void mynnsfc_water(const SfcIn& in, real lakemask, SfcPoint& o, const PsiTables& tb) {
  constexpr real prt = real(1.);
  real zol = o.zol;
  const real psfc = in.psfcpa / real(1000.);
  const real tsk = in.tskin;
  {
    real e1 = svp1 * std::exp(svp2 * (tsk - svpt0) / (tsk - svp3));
    if (lakemask < p5) e1 = e1 * real(0.9812);
    o.qsfc = ep2 * e1 / (psfc - ep3 * e1);
  }
  {
    const real e1 = svp1 * std::exp(svp2 * (in.t_1 - svpt0) / (in.t_1 - svp3));
    const real pl = in.p_1 / real(1000.);
    o.qgh = ep2 * e1 / (pl - e1);
  }
  const real thcon = std::pow(real(100000.) / in.psfcpa, rovcp);
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
  o.wstar = vconvc * std::pow(grav / tsk * in.pblh * fluxc, p333);
  const real vsgd = c_min(p25 * std::pow(c_max(in.dx / real(5000.) - one, zero), p333), p5);
  o.wspd = std::sqrt(o.wspd * o.wspd + o.wstar * o.wstar + vsgd * vsgd);
  o.wspd = c_max(o.wspd, wmin);
  o.rb = govrth * za * dthvdz / (o.wspd * o.wspd);
  o.rb = c_max(o.rb, real(-2.0));
  o.rb = c_min(o.rb, two);

  const real visc = real(1.326e-5) * (one + real(6.542e-3) * tc_1 + real(8.301e-6) * tc_1 * tc_1 -
                                      real(4.84e-9) * tc_1 * tc_1 * tc_1);
  // sfc_z0_type >= 0
  if (in.sfcflux_water == 0) coare30_znt(o.znt, o.ust, o.wspd, visc, za);
  else if (in.sfcflux_water == 1) edson_etal_2013(o.znt, o.ust, o.wspd, visc, za);
  else if (in.sfcflux_water == 2 || in.sfcflux_water == 3) davis_etal_2008(o.znt, o.ust);
  else taylor_yelland_2001(o.znt, o.ust, o.wspd);  // 4
  const real zntstoch = o.znt;
  const real restar = c_max(o.ust * zntstoch / visc, p1);
  real zt, zq;
  if (in.sfcflux_water == 0) fairall_etal_2003(zt, zq, restar);
  else if (in.sfcflux_water == 3) garratt_1992_water(zt, zq, zntstoch, restar);
  else fairall_etal_2014(zt, zq, restar);  // 1, 2, 4

  o.gz1oz0 = std::log((za + zntstoch) / zntstoch);
  const real gz1ozt = std::log((za + zntstoch) / zt);
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
  } else if (o.rb == real(0.)) {
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

  o.ust = p5 * o.ust + p5 * karman * o.wspd / o.psix;
  o.stress = o.ust * o.ust;
  const real wspdi = c_max(std::sqrt(in.u_1 * in.u_1 + in.v_1 * in.v_1), wmin);
  o.ustm = p5 * o.ustm + p5 * karman * wspdi / o.psix;
  real dtg = thv_1 - thvsk;
  o.mol = karman * dtg / o.psit / prt;
  const real dqg = (in.qv_1 - o.qsfc) * real(1000.);
  o.qstar = karman * dqg / psiq / prt;

  if (in.isfflx < 1) {
    o.qfx = zero; o.hfx = zero; o.flhc = zero; o.flqc = zero; o.lh = zero; o.chs = zero; o.ch = zero;
    o.chs2 = zero; o.cqs2 = zero; o.cqs = zero; o.cm = zero;
    o.ck = zero; o.cd = zero; o.cka = zero; o.cda = zero;
  } else {
    o.flqc = rho_1 * in.mavail * o.ust * karman / psiq;
    o.flhc = rho_1 * o.cpm * o.ust * karman / o.psit;
    o.qfx = o.flqc * (o.qsfc - in.qv_1);
    o.qfx = c_max(o.qfx, real(-0.02));
    o.lh = xlv * o.qfx;
    o.hfx = rho_1 * o.cpm * karman * o.wspd / o.psix * karman / o.psit * (thsk - th_1);
    if (in.sfcflux_water > 1) o.hfx = o.hfx + rho_1 * o.ustm * o.ustm * wspdi;  // ahw: dissipative heating
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
