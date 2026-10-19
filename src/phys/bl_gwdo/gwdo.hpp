#pragma once
// KIM orographic gravity wave drag (gwd_opt = 1): WRF 4.8.0 phys/physics_mmm/bl_gwdo.F90 bl_gwdo_run for one
// column (module_bl_gwdo.F gwdo passes the columns of a row; they are independent). Level arrays use the
// Fortran 1-based index k = 1 (lowest) .. kte, interfaces 1 .. kte + 1. Every statement follows the Fortran
// in order; x**2 and x**2. are x*x, x**real is pow, nint is lround.
#include <cmath>

#include "core/backend.hpp"

namespace cppwrf::phys::bl_gwdo {

inline constexpr int kMaxLevels = 100;
inline constexpr int KM = kMaxLevels + 2;

template <typename T>
struct GwdoColumn {
  T uproj[KM], vproj[KM], t1[KM], q1[KM], prslk[KM], zl[KM], prsl[KM], prsi[KM];
  T rublten[KM], rvblten[KM];  // in: the PBL tendencies; out: + the drag
  T dtaux3d[KM], dtauy3d[KM];
  T var, oc1, omax, oa2d1, oa2d2, oa2d3, oa2d4, ol2d1, ol2d2, ol2d3, ol2d4, sina, cosa, dxmeter;
  T dusfcg, dvsfcg;
};

template <typename T>
CPPWRF_HD inline void bl_gwdo_run(GwdoColumn<T>& c, int kte, T dx_factor, bool if_nonhyd, T g_, T cp_, T rd_, T fv_,
                                  T pi_, T deltim) {
  constexpr int kts = 1, kgwdmin = 2, mdir = 8;
  const T qmin = T(1.0e-30), ric = T(0.25), velmin = T(1.), rimin = T(-100.), bnv2min = T(1.0e-5), efacmin = T(0.);
  const T efacmax = T(10.), frc = T(1.0), frmax = T(10.), ce = T(0.8), cg = T(1.), hmt_min = T(50.);
  const T oc_min = T(1.), oc_max = T(10.), olmin = T(1.0e-5), odmin = T(0.1), odmax = T(10.);
  const T pcutoff = T(7.5e2), pcutoff_den = T(0.01);
  constexpr int nwdir[mdir] = {6, 7, 5, 8, 2, 3, 1, 4};

  const int lcap = kte;
  const T fdir = T(mdir) / (T(2.0) * pi_);
  const int kgwdmax = kte / 2;  // maximum height for gwd stress
  // grid lengths for the flow-blocking drag
  const T delx = c.dxmeter, dely = c.dxmeter;
  const T area = delx * dely;
  T dxy4[5], dxy4p[5];
  dxy4[1] = delx;
  dxy4[2] = dely;
  dxy4[3] = std::sqrt(delx * delx + dely * dely);
  dxy4[4] = dxy4[3];
  dxy4p[1] = dxy4[2];
  dxy4p[2] = dxy4[1];
  dxy4p[3] = dxy4[4];
  dxy4p[4] = dxy4[3];
  const T dx_eff = dx_factor * (delx + dely) * T(0.5);
  // initialisation
  bool ldrag = false, icrilv = false, flag = true;
  int kref = 0, komax = 0, kbomax = 0;
  T dtaux = T(0.), dtauy = T(0.), xn = T(0.), yn = T(0.), ubar = T(0.), vbar = T(0.), rhobar = T(0.), ulow = T(0.);
  T oa = T(0.), ol = T(0.), oc = T(0.), taub = T(0.), xlinv = T(1.0), denfac = T(1.0), nhd_effect = T(1.0);
  T coefm = T(0.), bnv = T(0.), fr = T(0.), brvf = T(0.), rulow = T(0.), olp = T(0.), od = T(0.), dxy = T(0.), dxyp = T(0.);
  T usqj[KM], bnv2[KM], vtj[KM], vtk[KM], taup[KM], taud[KM], dtfac[KM], rho[KM], del[KM], u1[KM], v1[KM];
  T velco[KM], taufb[KM], dudt[KM], dvdt[KM], dtaux2d[KM], dtauy2d[KM];
  for (int k = 0; k < KM; ++k) {
    taufb[k] = T(0.); usqj[k] = T(0.); bnv2[k] = T(0.); vtj[k] = T(0.); vtk[k] = T(0.); taup[k] = T(0.);
    taud[k] = T(0.); dtaux2d[k] = T(0.); dtauy2d[k] = T(0.); dtfac[k] = T(1.0);
  }
  for (int k = kts; k <= kte; ++k) {
    vtj[k] = c.t1[k] * (T(1.) + fv_ * c.q1[k]);
    vtk[k] = vtj[k] / c.prslk[k];
    rho[k] = T(1.) / rd_ * c.prsl[k] / vtj[k];
    del[k] = c.prsi[k] - c.prsi[k + 1];
    u1[k] = c.uproj[k] * c.cosa - c.vproj[k] * c.sina;
    v1[k] = c.uproj[k] * c.sina + c.vproj[k] * c.cosa;
  }
  const T zref = T(2.) * c.var;
  flag = true;
  for (int k = kts + 1; k <= kte; ++k) {
    if (flag && c.zl[k] - c.zl[1] >= zref) {
      kref = k + 1;
      flag = false;
    }
  }
  kref = c_max(c_min(kref, kgwdmax), kgwdmin);
  flag = true;
  for (int k = kts + 1; k <= kte; ++k) {
    if (flag && c.zl[k] >= c.omax) {
      komax = k + 1;
      flag = false;
    }
  }
  komax = c_max(komax, kref);
  flag = true;
  for (int k = kts + 1; k <= kte; ++k) {
    if (flag && c.zl[k] >= (c.omax + zref)) {
      kbomax = k + 1;
      flag = false;
    }
  }
  const T delks = T(1.0) / (c.prsi[1] - c.prsi[kref]);
  const T delks1 = T(1.0) / (c.prsl[1] - c.prsl[kref]);
  // low level average below the reference level
  for (int k = kts; k <= kgwdmax; ++k) {
    if (k < kref) {
      const T rcsks = del[k] * delks;
      const T rdelks = del[k] * delks;
      ubar = ubar + rcsks * u1[k];
      vbar = vbar + rcsks * v1[k];
      rhobar = rhobar + rdelks * rho[k];
    }
  }
  // low-level horizontal wind direction
  T oa4[5], ol4[5], ol4p[5];
  oa4[1] = c.oa2d1; oa4[2] = c.oa2d2; oa4[3] = c.oa2d3; oa4[4] = c.oa2d4;
  ol4[1] = c.ol2d1; ol4[2] = c.ol2d2; ol4[3] = c.ol2d3; ol4[4] = c.ol2d4;
  {
    const T wdir = std::atan2(ubar, vbar) + pi_;
    const int idir = int(std::lround(fdir * wdir) % mdir) + 1;
    const int nwd = nwdir[idir - 1];
    oa = T(1 - 2 * ((nwd - 1) / 4)) * oa4[(nwd - 1) % 4 + 1];
    ol = c_max(ol4[(nwd - 1) % 4 + 1], olmin);
    oc = c_min(c_max(c.oc1, oc_min), oc_max);
    ol4p[1] = ol4[2]; ol4p[2] = ol4[1]; ol4p[3] = ol4[4]; ol4p[4] = ol4[3];
    olp = c_max(ol4p[(nwd - 1) % 4 + 1], olmin);
    od = olp / ol;
    od = c_min(od, odmax);
    od = c_max(od, odmin);
    dxy = dxy4[(nwd - 1) % 4 + 1];
    dxyp = dxy4p[(nwd - 1) % 4 + 1];
  }
  // richardson number in usqj
  for (int k = kts; k <= kte - 1; ++k) {
    const T ti = T(2.0) / (c.t1[k] + c.t1[k + 1]);
    const T rdz = T(1.) / (c.zl[k + 1] - c.zl[k]);
    const T tem1 = u1[k] - u1[k + 1];
    const T tem2 = v1[k] - v1[k + 1];
    const T dw2 = tem1 * tem1 + tem2 * tem2;
    const T shr2 = c_max(dw2, velmin) * rdz * rdz;
    const T bvf2 = g_ * (g_ / cp_ + rdz * (vtj[k + 1] - vtj[k])) * ti;
    usqj[k] = c_max(bvf2 / shr2, rimin);
    bnv2[k] = T(2.0) * g_ * rdz * (vtk[k + 1] - vtk[k]) / (vtk[k + 1] + vtk[k]);
  }
  ulow = c_max(std::sqrt(ubar * ubar + vbar * vbar), T(1.0));
  rulow = T(1.) / ulow;
  // horizontal wind speed projected to the low-level wind
  for (int k = kts; k <= kte - 1; ++k) {
    velco[k] = T(0.5) * ((u1[k] + u1[k + 1]) * ubar + (v1[k] + v1[k + 1]) * vbar);
    velco[k] = velco[k] * rulow;
    if ((velco[k] < velmin) && (velco[k] > T(0.))) velco[k] = velmin;
  }
  ldrag = c.omax <= hmt_min;
  ldrag = ldrag || c.var <= T(0.);
  ldrag = ldrag || velco[1] <= T(0.);
  ldrag = ldrag || ulow == T(1.0);
  for (int k = kgwdmin; k <= kgwdmax; ++k)
    if (k < kref) ldrag = ldrag || velco[k] <= T(0.) || bnv2[k] < T(0.);
  // low level weighted average ri (usqj) and n**2 (bnv2)
  {
    const T wtkbj = (c.prsl[1] - c.prsl[2]) * delks1;
    bnv2[1] = wtkbj * bnv2[1];
    usqj[1] = wtkbj * usqj[1];
  }
  for (int k = kgwdmin; k <= kgwdmax; ++k) {
    if (k < kref) {
      const T rdelks = (c.prsl[k] - c.prsl[k + 1]) * delks1;
      bnv2[1] = bnv2[1] + bnv2[k] * rdelks;
      usqj[1] = usqj[1] + usqj[k] * rdelks;
    }
  }
  ldrag = ldrag || bnv2[1] <= T(0.);
  for (int k = kgwdmin; k <= kgwdmax; ++k)
    if (k < kref) usqj[k] = usqj[1];
  if (!ldrag) {
    bnv = std::sqrt(bnv2[1]);
    fr = bnv * rulow * c.var * od;
    fr = c_min(fr, frmax);
    xn = ubar * rulow;
    yn = vbar * rulow;
  }
  // base level stress
  if (!ldrag) {
    T efac = std::pow(oa + T(2.), ce * fr / frc);
    efac = c_min(c_max(efac, efacmin), efacmax);
    coefm = std::pow(T(1.) + ol, oa + T(1.));
    xlinv = coefm / dx_eff;
    T tem = fr * fr * oc;
    const T gfac = tem / (tem + cg);
    taub = rhobar * efac * xlinv * gfac * ulow * ulow * ulow / bnv;
    if (if_nonhyd) {
      tem = fr * fr;
      nhd_effect = -(T(9.) / T(8.)) * tem +
                   std::exp(-(T(2.) / fr)) * (-(T(5.) / T(4.) / tem) - T(0.5) / fr + T(5.) / T(4.) + T(9.) / T(4.) * fr + T(9.) / T(8.) * tem);
      taub = taub * (T(1.) + nhd_effect);
    }
  } else {
    taub = T(0.);
    xn = T(0.);
    yn = T(0.);
  }
  // vertical structure of the stress
  for (int k = kts; k <= kgwdmax; ++k)
    if (k <= kref) taup[k] = taub;
  for (int k = kgwdmin; k <= kte - 1; ++k) {
    const int kp1 = k + 1;
    if (k >= kref) {
      icrilv = icrilv || (usqj[k] < ric) || (velco[k] <= T(0.));
      brvf = c_max(bnv2[k], bnv2min);
      brvf = std::sqrt(brvf);
    }
    if (k >= kref && !ldrag) {
      if (!icrilv && taup[k] > T(0.)) {
        const T temv = T(1.0) / velco[k];
        const T tem1 = coefm / dxy * (rho[kp1] + rho[k]) * brvf * velco[k] * T(0.5);
        T hd = std::sqrt(taup[k] / tem1);
        const T fro = brvf * hd * temv;
        const T tem2 = std::sqrt(usqj[k]);
        const T tem = T(1.) + tem2 * fro;
        const T rim = usqj[k] * (T(1.) - fro) / (tem * tem);
        if (rim <= ric) {
          const T temc = T(2.0) + T(1.0) / tem2;
          hd = velco[k] * (T(2.) * std::sqrt(temc) - temc) / brvf;
          taup[kp1] = tem1 * hd * hd;
        } else {
          taup[kp1] = taup[k];
        }
      }
    }
  }
  if (lcap < kte) {
    for (int klcap = lcap + 1; klcap <= kte; ++klcap) taup[klcap] = c.prsi[klcap] / c.prsi[lcap] * taup[lcap];
  }
  if (!ldrag) {
    // height of the flow-blocking layer
    int kblk = 0;
    T fbdpe = T(0.), fbdke = T(0.), zblk = T(0.);
    for (int k = kte; k >= kgwdmin; --k) {
      if (kblk == 0 && k <= kbomax) {
        fbdpe = fbdpe + bnv2[k] * (c.zl[komax] - c.zl[k]) * del[k] / g_ / rho[k];
        fbdke = T(0.5) * (u1[k] * u1[k] + v1[k] * v1[k]);
        if (fbdpe >= fbdke && k <= komax) {
          kblk = k;
          zblk = c.zl[kblk] - c.zl[kts];
        }
      }
    }
    if (kblk != 0) {
      const T fbdcd = c_max(T(2.0) - T(1.0) / od, T(0.));
      taufb[kts] = T(0.5) * rhobar * coefm / area * fbdcd * dxyp * olp * zblk * (ulow * ulow);
      const T tautem = taufb[kts] / T(kblk - kts);
      for (int k = kts + 1; k <= kblk; ++k) taufb[k] = taufb[k - 1] - tautem;
      for (int k = kts; k <= kte + 1; ++k) taup[k] = taup[k] + taufb[k];
    }
  }
  // - g d(tau)/d(pressure), deceleration
  for (int k = kts; k <= kte; ++k) taud[k] = T(1.) * (taup[k + 1] - taup[k]) * g_ / del[k];
  // damping factor against wind reversal
  for (int k = kts; k <= kgwdmax - 1; ++k)
    if (std::fabs(taud[k]) >= qmin) dtfac[k] = c_min(dtfac[k], std::fabs(velco[k] / (deltim * taud[k])));
  // mesosphere drag limiter
  for (int k = kgwdmax; k <= kte - 1; ++k) {
    if (std::fabs(taud[k]) >= qmin && c.prsl[k] <= pcutoff) {
      denfac = c_min(rho[k] / pcutoff_den, T(1.));
      dtfac[k] = c_min(dtfac[k], denfac * std::fabs(velco[k] / (deltim * taud[k])));
    }
  }
  dtfac[kte] = dtfac[kte - 1];
  T dusfc = T(0.), dvsfc = T(0.);
  for (int k = kts; k <= kte; ++k) {
    taud[k] = taud[k] * dtfac[k];
    dtaux = taud[k] * xn;
    dtauy = taud[k] * yn;
    dtaux2d[k] = dtaux;
    dtauy2d[k] = dtauy;
    dudt[k] = dtaux;
    dvdt[k] = dtauy;
    dusfc = dusfc + dtaux * del[k];
    dvsfc = dvsfc + dtauy * del[k];
  }
  dusfc = (-(T(1.) / g_)) * dusfc;
  dvsfc = (-(T(1.) / g_)) * dvsfc;
  // rotate the tendencies back to the model grid
  for (int k = kts; k <= kte; ++k) {
    c.rublten[k] = c.rublten[k] + dudt[k] * c.cosa + dvdt[k] * c.sina;
    c.rvblten[k] = c.rvblten[k] - dudt[k] * c.sina + dvdt[k] * c.cosa;
    c.dtaux3d[k] = dtaux2d[k] * c.cosa + dtauy2d[k] * c.sina;
    c.dtauy3d[k] = -dtaux2d[k] * c.sina + dtauy2d[k] * c.cosa;
  }
  c.dusfcg = dusfc * c.cosa + dvsfc * c.sina;
  c.dvsfcg = -dusfc * c.sina + dvsfc * c.cosa;
}

}  // namespace cppwrf::phys::bl_gwdo
