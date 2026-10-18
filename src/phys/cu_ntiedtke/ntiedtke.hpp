#pragma once
// New Tiedtke cumulus scheme (cu_physics = 16): WRF 4.8.0 phys/physics_mmm/cu_ntiedtke.F90 for one column
// (klon = 1). Arrays keep the Fortran 1-based level index (k = 1 is the model TOP, k = klev the lowest
// level; klev + 1 is the surface for half levels), sized kMaxLevels + 2. Every statement follows the
// Fortran in order; x**2 is x*x (gfortran powi), x**real is pow.
#include <cmath>

#include "core/backend.hpp"

namespace cppwrf::phys::cu_ntiedtke {

inline constexpr int kMaxLevels = 100;
inline constexpr int KM = kMaxLevels + 2;

// cu_ntiedtke_common (module data and cu_ntiedtke_init)
template <typename T>
struct NtConst {
  T alf, als, alv, cpd, g, rd, rv;
  T rcpd, c2es, c5les, c5ies, r5alvcp, r5alscp, ralvdcp, ralsdcp, ralfdcp, vtmpc1, zrg;
  static constexpr double t13 = 1.0 / 3.0;
  static constexpr double tmelt = 273.16, c1es = 610.78, c3les = 17.2693882, c3ies = 21.875;
  static constexpr double c4les = 35.86, c4ies = 7.66;
  static constexpr double rtwat = tmelt, rtber = tmelt - 5., rtice = tmelt - 23.;
  static constexpr int momtrans = 2;
  static constexpr double entrdd = 2.0e-4, cmfcmax = 1.0, cmfcmin = 1.e-10, cmfdeps = 0.30, zdnoprc = 2.0e4;
  static constexpr double cprcon = 1.4e-3, pgcoef = 0.7;
  static constexpr bool nonequil = true, lmfpen = true, lmfmid = true, lmfscv = true, lmfdd = true, lmfdudv = true;
};

template <typename T>
inline NtConst<T> make_ntconst(T con_cp, T con_rd, T con_rv, T con_xlv, T con_xls, T con_xlf, T con_grav) {
  using C = NtConst<T>;
  NtConst<T> c{};
  c.alf = con_xlf; c.als = con_xls; c.alv = con_xlv; c.cpd = con_cp; c.g = con_grav; c.rd = con_rd; c.rv = con_rv;
  c.rcpd = T(1.0) / con_cp;
  c.c2es = T(C::c1es) * con_rd / con_rv;
  c.c5les = T(C::c3les) * (T(C::tmelt) - T(C::c4les));
  c.c5ies = T(C::c3ies) * (T(C::tmelt) - T(C::c4ies));
  c.r5alvcp = c.c5les * con_xlv * c.rcpd;
  c.r5alscp = c.c5ies * con_xls * c.rcpd;
  c.ralvdcp = con_xlv * c.rcpd;
  c.ralsdcp = con_xls * c.rcpd;
  c.ralfdcp = con_xlf * c.rcpd;
  c.vtmpc1 = con_rv / con_rd - T(1.0);
  c.zrg = T(1.0) / con_grav;
  return c;
}

// statement functions (end of cu_ntiedtke.F90)
template <typename T>
CPPWRF_HD inline T foealfa(T tt) {
  using C = NtConst<T>;
  const T x = (std::fmax(T(C::rtice), std::fmin(T(C::rtwat), tt)) - T(C::rtice)) / (T(C::rtwat) - T(C::rtice));
  return std::fmin(T(1.), x * x);
}
template <typename T>
CPPWRF_HD inline T foelhm(const NtConst<T>& c, T tt) {
  return foealfa(tt) * c.alv + (T(1.) - foealfa(tt)) * c.als;
}
template <typename T>
CPPWRF_HD inline T foeewm(const NtConst<T>& c, T tt) {
  using C = NtConst<T>;
  return c.c2es * (foealfa(tt) * std::exp(T(C::c3les) * (tt - T(C::tmelt)) / (tt - T(C::c4les))) +
                   (T(1.) - foealfa(tt)) * std::exp(T(C::c3ies) * (tt - T(C::tmelt)) / (tt - T(C::c4ies))));
}
template <typename T>
CPPWRF_HD inline T foedem(const NtConst<T>& c, T tt) {
  using C = NtConst<T>;
  const T a = tt - T(C::c4les), b = tt - T(C::c4ies);
  return foealfa(tt) * c.r5alvcp * (T(1.) / (a * a)) + (T(1.) - foealfa(tt)) * c.r5alscp * (T(1.) / (b * b));
}
template <typename T>
CPPWRF_HD inline T foeldcpm(const NtConst<T>& c, T tt) {
  return foealfa(tt) * c.ralvdcp + (T(1.) - foealfa(tt)) * c.ralsdcp;
}

// cuadjtqn: saturation adjustment of pt(kk), pq(kk) at pressure psp
template <typename T>
CPPWRF_HD inline void cuadjtqn(const NtConst<T>& c, int kk, T psp, T* pt, T* pq, bool ldflag, int kcall) {
  using C = NtConst<T>;
  const T c3les = T(C::c3les), c3ies = T(C::c3ies), c4les = T(C::c4les), c4ies = T(C::c4ies), tmelt = T(C::tmelt);
  if (kcall == 1) {
    if (ldflag) {
      const T zqp = T(1.) / psp;
      T zl = T(1.) / (pt[kk] - c4les);
      T zi = T(1.) / (pt[kk] - c4ies);
      T zqsat = c.c2es * (foealfa(pt[kk]) * std::exp(c3les * (pt[kk] - tmelt) * zl) +
                          (T(1.) - foealfa(pt[kk])) * std::exp(c3ies * (pt[kk] - tmelt) * zi));
      zqsat = zqsat * zqp;
      zqsat = std::fmin(T(0.5), zqsat);
      T zcor = T(1.) - c.vtmpc1 * zqsat;
      T zf = foealfa(pt[kk]) * c.r5alvcp * (zl * zl) + (T(1.) - foealfa(pt[kk])) * c.r5alscp * (zi * zi);
      const T zcond = (pq[kk] * (zcor * zcor) - zqsat * zcor) / ((zcor * zcor) + zqsat * zf);
      if (zcond > T(0.)) {
        pt[kk] = pt[kk] + foeldcpm(c, pt[kk]) * zcond;
        pq[kk] = pq[kk] - zcond;
        zl = T(1.) / (pt[kk] - c4les);
        zi = T(1.) / (pt[kk] - c4ies);
        zqsat = c.c2es * (foealfa(pt[kk]) * std::exp(c3les * (pt[kk] - tmelt) * zl) +
                          (T(1.) - foealfa(pt[kk])) * std::exp(c3ies * (pt[kk] - tmelt) * zi));
        zqsat = zqsat * zqp;
        zqsat = std::fmin(T(0.5), zqsat);
        zcor = T(1.) - c.vtmpc1 * zqsat;
        zf = foealfa(pt[kk]) * c.r5alvcp * (zl * zl) + (T(1.) - foealfa(pt[kk])) * c.r5alscp * (zi * zi);
        T zcond1 = (pq[kk] * (zcor * zcor) - zqsat * zcor) / ((zcor * zcor) + zqsat * zf);
        if (std::fabs(zcond) < T(1.e-20)) zcond1 = T(0.);
        pt[kk] = pt[kk] + foeldcpm(c, pt[kk]) * zcond1;
        pq[kk] = pq[kk] - zcond1;
      }
    }
  } else if (kcall == 2) {
    if (ldflag) {
      const T zqp = T(1.) / psp;
      T zqsat = foeewm(c, pt[kk]) * zqp;
      zqsat = std::fmin(T(0.5), zqsat);
      T zcor = T(1.) / (T(1.) - c.vtmpc1 * zqsat);
      zqsat = zqsat * zcor;
      T zcond = (pq[kk] - zqsat) / (T(1.) + zqsat * zcor * foedem(c, pt[kk]));
      zcond = std::fmin(zcond, T(0.));
      pt[kk] = pt[kk] + foeldcpm(c, pt[kk]) * zcond;
      pq[kk] = pq[kk] - zcond;
      zqsat = foeewm(c, pt[kk]) * zqp;
      zqsat = std::fmin(T(0.5), zqsat);
      zcor = T(1.) / (T(1.) - c.vtmpc1 * zqsat);
      zqsat = zqsat * zcor;
      T zcond1 = (pq[kk] - zqsat) / (T(1.) + zqsat * zcor * foedem(c, pt[kk]));
      if (std::fabs(zcond) < T(1.e-20)) zcond1 = std::fmin(zcond1, T(0.));
      pt[kk] = pt[kk] + foeldcpm(c, pt[kk]) * zcond1;
      pq[kk] = pq[kk] - zcond1;
    }
  } else if (kcall == 0) {
    const T zqp = T(1.) / psp;
    T zqsat = foeewm(c, pt[kk]) * zqp;
    zqsat = std::fmin(T(0.5), zqsat);
    T zcor = T(1.) / (T(1.) - c.vtmpc1 * zqsat);
    zqsat = zqsat * zcor;
    T zcond1 = (pq[kk] - zqsat) / (T(1.) + zqsat * zcor * foedem(c, pt[kk]));
    pt[kk] = pt[kk] + foeldcpm(c, pt[kk]) * zcond1;
    pq[kk] = pq[kk] - zcond1;
    zqsat = foeewm(c, pt[kk]) * zqp;
    zqsat = std::fmin(T(0.5), zqsat);
    zcor = T(1.) / (T(1.) - c.vtmpc1 * zqsat);
    zqsat = zqsat * zcor;
    zcond1 = (pq[kk] - zqsat) / (T(1.) + zqsat * zcor * foedem(c, pt[kk]));
    pt[kk] = pt[kk] + foeldcpm(c, pt[kk]) * zcond1;
    pq[kk] = pq[kk] - zcond1;
  }
}

// cuinin: half-level environment, initial updraft/downdraft profiles
template <typename T>
CPPWRF_HD inline void cuinin(const NtConst<T>& c, int klev, const T* pten, const T* pqen, const T* pqsen, const T* puen,
                             const T* pven, const T* pverv, const T* pgeo, const T* paph, const T* pgeoh, T* ptenh,
                             T* pqenh, T* pqsenh, int& klwmin, T* ptu, T* pqu, T* ptd, T* pqd, T* puu, T* pvu, T* pud,
                             T* pvd, T* pmfu, T* pmfd, T* pmfus, T* pmfds, T* pmfuq, T* pmfdq, T* pdmfup, T* pdmfdp,
                             T* pdpmel, T* plu, T* plude, int* klab) {
  (void)pmfu; (void)pmfd; (void)pmfus; (void)pmfds; (void)pmfuq; (void)pmfdq; (void)pdmfup; (void)pdmfdp;
  (void)pdpmel; (void)plude;
  for (int jk = 2; jk <= klev; ++jk) {
    ptenh[jk] = (std::fmax(c.cpd * pten[jk - 1] + pgeo[jk - 1], c.cpd * pten[jk] + pgeo[jk]) - pgeoh[jk]) * c.rcpd;
    pqenh[jk] = pqen[jk - 1];
    pqsenh[jk] = pqsen[jk - 1];
    const T zph = paph[jk];
    if (jk >= klev - 1 || jk < 2) continue;
    cuadjtqn(c, jk, zph, ptenh, pqsenh, true, 0);
    pqenh[jk] = std::fmin(pqen[jk - 1], pqsen[jk - 1]) + (pqsenh[jk] - pqsen[jk - 1]);
    pqenh[jk] = std::fmax(pqenh[jk], T(0.));
  }
  ptenh[klev] = (c.cpd * pten[klev] + pgeo[klev] - pgeoh[klev]) * c.rcpd;
  pqenh[klev] = pqen[klev];
  ptenh[1] = pten[1];
  pqenh[1] = pqen[1];
  klwmin = klev;
  T zwmax = T(0.);
  for (int jk = klev - 1; jk >= 2; --jk) {
    const T zzs = std::fmax(c.cpd * ptenh[jk] + pgeoh[jk], c.cpd * ptenh[jk + 1] + pgeoh[jk + 1]);
    ptenh[jk] = (zzs - pgeoh[jk]) * c.rcpd;
  }
  for (int jk = klev; jk >= 3; --jk) {
    if (pverv[jk] < zwmax) {
      zwmax = pverv[jk];
      klwmin = jk;
    }
  }
  for (int jk = 1; jk <= klev; ++jk) {
    int ik = jk - 1;
    if (jk == 1) ik = 1;
    ptu[jk] = ptenh[jk];
    ptd[jk] = ptenh[jk];
    pqu[jk] = pqenh[jk];
    pqd[jk] = pqenh[jk];
    plu[jk] = T(0.);
    puu[jk] = puen[ik];
    pud[jk] = puen[ik];
    pvu[jk] = pven[ik];
    pvd[jk] = pven[ik];
    klab[jk] = 0;
  }
}

// cubasmcn: mid-level convection base
template <typename T>
CPPWRF_HD inline void cubasmcn(const NtConst<T>& c, int kk, const T* pten, const T* pqen, const T* pqsen,
                               const T* pverv, const T* pgeo, const T* pgeoh, bool ldcum, int& ktype, int* klab,
                               T* plrain, T* pmfu, T& pmfub, int& kcbot, T* ptu, T* pqu, T* plu, T* pmfus, T* pmfuq,
                               T* pmful, T* pdmfup) {
  using C = NtConst<T>;
  if (!ldcum && klab[kk + 1] == 0) {
    if (C::lmfmid && pqen[kk] > T(0.80) * pqsen[kk] && pgeo[kk] * c.zrg > T(5.0e2) && pgeo[kk] * c.zrg < T(1.0e4)) {
      ptu[kk + 1] = (c.cpd * pten[kk] + pgeo[kk] - pgeoh[kk + 1]) * c.rcpd;
      pqu[kk + 1] = pqen[kk];
      plu[kk + 1] = T(0.);
      T zzzmb = std::fmax(T(C::cmfcmin), -pverv[kk] * c.zrg);
      zzzmb = std::fmin(zzzmb, T(C::cmfcmax));
      pmfub = zzzmb;
      pmfu[kk + 1] = pmfub;
      pmfus[kk + 1] = pmfub * (c.cpd * ptu[kk + 1] + pgeoh[kk + 1]);
      pmfuq[kk + 1] = pmfub * pqu[kk + 1];
      pmful[kk + 1] = T(0.);
      pdmfup[kk + 1] = T(0.);
      kcbot = kk;
      klab[kk + 1] = 1;
      plrain[kk + 1] = T(0.0);
      ktype = 3;
    }
  }
}

// cuentrn (zentr = 0: no entrainment from this routine, detrainment 0.75e-4 per m)
template <typename T>
CPPWRF_HD inline void cuentrn(const NtConst<T>& c, int kk, int kcbot, bool ldcum, bool ldwork, const T* pgeoh,
                              const T* pmfu, T& pdmfen, T& pdmfde) {
  if (ldwork) {
    pdmfen = T(0.);
    pdmfde = T(0.);
    const T zentr = T(0.);
    if (ldcum) {
      const T zdz = (pgeoh[kk] - pgeoh[kk + 1]) * c.zrg;
      const T zmf = pmfu[kk + 1] * zdz;
      const bool llo1 = kk < kcbot;
      if (llo1) {
        pdmfen = zentr * zmf;
        pdmfde = T(0.75e-4) * zmf;
      }
    }
  }
}

// cutypen: convection type and cloud base/top from test parcels (shallow from the surface layer, then
// deep parcels launched from each level of the lower half)
template <typename T>
CPPWRF_HD inline void cutypen(const NtConst<T>& c, int klev, const T* pqen, const T* ptenh, const T* pqenh,
                              const T* pqsenh, const T* pgeoh, const T* paph, T hfx, T qfx, const T* pgeo,
                              const T* pqsen, const T* pap, const T* pten, int lndj, T* cutu, T* cuqu, int* culab,
                              bool& ldcum, int& cubot, int& cutop, int& ktype, T& wbase, T* culu, int& kdpl) {
  using C = NtConst<T>;
  (void)pqsenh; (void)lndj;
  const int klevm1 = klev - 1;
  const T c4les = T(C::c4les), c4ies = T(C::c4ies);
  int kcbot = klev, kctop = klev;
  kdpl = klev;
  ktype = 0;
  wbase = T(0.);
  ldcum = false;
  T plu[KM], ptu[KM], pqu[KM], dh[KM], dhen[KM], kup[KM], vptu[KM], vten[KM], zbuo[KM], abuoy[KM];
  int klab[KM];
  for (int jk = 1; jk <= klev; ++jk) {
    plu[jk] = culu[jk];
    ptu[jk] = cutu[jk];
    pqu[jk] = cuqu[jk];
    klab[jk] = culab[jk];
    dh[jk] = T(0.0);
    dhen[jk] = T(0.0);
    kup[jk] = T(0.0);
    vptu[jk] = T(0.0);
    vten[jk] = T(0.0);
    zbuo[jk] = T(0.0);
    abuoy[jk] = T(0.0);
  }
  T zqold = T(0.);
  bool lldcum = false;
  bool loflag = true;
  T eta = T(0), dz = T(0), coef = T(0), zph = T(0);
  int zcbase = 0;

  // test parcel shared by the column: depth-dependent entrainment, saturation adjustment, kinetic energy
  auto parcel_step = [&](int jk, bool deep) {
    if (loflag) {
      T zdq = std::fmax((zqold - pqu[jk]), T(0.));
      plu[jk] = plu[jk + 1] + zdq;
      const T zlglac = zdq * ((T(1.) - foealfa(ptu[jk])) - (T(1.) - foealfa(ptu[jk + 1])));
      if (deep) plu[jk] = T(0.5) * plu[jk];
      else plu[jk] = std::fmin(plu[jk], T(5.e-3));
      dh[jk] = pgeoh[jk] + c.cpd * (ptu[jk] + c.ralfdcp * zlglac);
      vptu[jk] = ptu[jk] * (T(1.) + c.vtmpc1 * pqu[jk] - plu[jk]) + c.ralfdcp * zlglac;
      vten[jk] = ptenh[jk] * (T(1.) + c.vtmpc1 * pqenh[jk]);
      zbuo[jk] = (vptu[jk] - vten[jk]) / vten[jk];
      abuoy[jk] = (zbuo[jk] + zbuo[jk + 1]) * T(0.5) * c.g;
      const T atop1 = T(1.0) - T(2.) * coef;
      const T atop2 = T(2.0) * dz * abuoy[jk];
      const T abot = T(1.0) + T(2.) * coef;
      kup[jk] = (atop1 * kup[jk + 1] + atop2) / abot;
      if (plu[jk] > T(0.) && klab[jk + 1] == 1) {
        const int ik = jk + 1;
        T zqsu = foeewm(c, ptu[ik]) / paph[ik];
        zqsu = std::fmin(T(0.5), zqsu);
        T zcor = T(1.) / (T(1.) - c.vtmpc1 * zqsu);
        zqsu = zqsu * zcor;
        zdq = std::fmin(T(0.), pqu[ik] - zqsu);
        const T zalfaw = foealfa(ptu[ik]);
        const T aw = ptu[ik] - c4les, ai = ptu[ik] - c4ies;
        const T zfacw = c.c5les / (aw * aw);
        const T zfaci = c.c5ies / (ai * ai);
        const T zfac = zalfaw * zfacw + (T(1.) - zalfaw) * zfaci;
        const T zesdp = foeewm(c, ptu[ik]) / paph[ik];
        zcor = T(1.) / (T(1.) - c.vtmpc1 * zesdp);
        const T zdqsdt = zfac * zcor * zqsu;
        const T zdtdp = c.rd * ptu[ik] / (c.cpd * paph[ik]);
        const T zdp = zdq / (zdqsdt * zdtdp);
        zcbase = int(paph[ik] + zdp);  // INTEGER zcbase in WRF: truncated
        const T zpdifftop = T(zcbase) - paph[jk];
        const T zpdiffbot = paph[jk + 1] - T(zcbase);
        if (zpdifftop > zpdiffbot && kup[jk + 1] > T(0.)) {
          const int ikb = (klev - 1 < jk + 1) ? klev - 1 : jk + 1;
          klab[ikb] = 2;
          klab[jk] = 2;
          kcbot = ikb;
          plu[jk + 1] = T(1.0e-8);
        } else if (zpdifftop <= zpdiffbot && kup[jk] > T(0.)) {
          klab[jk] = 2;
          kcbot = jk;
        }
      }
      if (kup[jk] < T(0.)) {
        loflag = false;
        if (plu[jk + 1] > T(0.)) {
          kctop = jk;
          lldcum = true;
        } else {
          lldcum = false;
        }
      } else {
        if (plu[jk] > T(0.)) klab[jk] = 2;
        else klab[jk] = 1;
      }
    }
  };

  // shallow convection test from the lowest level
  for (int jk = klevm1; jk >= 2; --jk) {
    if (jk == klevm1) {
      const T rho = pap[klev] / (c.rd * (pten[klev] * (T(1.) + c.vtmpc1 * pqen[klev])));
      const T part1 = T(1.5) * T(0.4) * pgeo[klev] / (rho * pten[klev]);
      const T part2 = -hfx * c.rcpd - c.vtmpc1 * pten[klev] * qfx;
      const T root = T(0.001) - part1 * part2;
      if (part2 < T(0.)) {
        const T conw = T(1.2) * std::pow(root, T(C::t13));
        const T deltt = std::fmax(T(1.5) * hfx / (rho * c.cpd * conw), T(0.));
        const T deltq = std::fmax(T(1.5) * qfx / (rho * conw), T(0.));
        kup[klev] = T(0.5) * (conw * conw);
        pqu[klev] = pqenh[klev] + deltq;
        dhen[klev] = pgeoh[klev] + ptenh[klev] * c.cpd;
        dh[klev] = dhen[klev] + deltt * c.cpd;
        ptu[klev] = (dh[klev] - pgeoh[klev]) * c.rcpd;
        vptu[klev] = ptu[klev] * (T(1.) + c.vtmpc1 * pqu[klev]);
        vten[klev] = ptenh[klev] * (T(1.) + c.vtmpc1 * pqenh[klev]);
        zbuo[klev] = (vptu[klev] - vten[klev]) / vten[klev];
        klab[klev] = 1;
      } else {
        loflag = false;
      }
    }
    if (!loflag) break;
    if (loflag) {
      eta = T(0.8) / (pgeo[jk] * c.zrg) + T(2.e-4);
      dz = (pgeoh[jk] - pgeoh[jk + 1]) * c.zrg;
      coef = T(0.5) * eta * dz;
      dhen[jk] = pgeoh[jk] + c.cpd * ptenh[jk];
      dh[jk] = (coef * (dhen[jk + 1] + dhen[jk]) + (T(1.) - coef) * dh[jk + 1]) / (T(1.) + coef);
      pqu[jk] = (coef * (pqenh[jk + 1] + pqenh[jk]) + (T(1.) - coef) * pqu[jk + 1]) / (T(1.) + coef);
      ptu[jk] = (dh[jk] - pgeoh[jk]) * c.rcpd;
      zqold = pqu[jk];
      zph = paph[jk];
    }
    cuadjtqn(c, jk, zph, ptu, pqu, loflag, 1);
    parcel_step(jk, false);
  }
  {
    const int ikb = kcbot, ikt = kctop;
    if (paph[ikb] - paph[ikt] > T(C::zdnoprc)) lldcum = false;
    if (lldcum) {
      ktype = 2;
      ldcum = true;
      wbase = std::sqrt(std::fmax(T(2.) * kup[ikb], T(0.)));
      cubot = ikb;
      cutop = ikt;
      kdpl = klev;
    } else {
      cutop = -1;
      cubot = -1;
      kdpl = klev - 1;
      ldcum = false;
      wbase = T(0.);
    }
  }
  for (int jk = klev; jk >= 1; --jk) {
    const int ikt = kctop;
    if (jk >= ikt) {
      culab[jk] = klab[jk];
      cutu[jk] = ptu[jk];
      cuqu[jk] = pqu[jk];
      culu[jk] = plu[jk];
    }
  }

  // deep convection tests
  const T deltt = T(0.2);
  const T deltq = T(1.0e-4);
  bool deepflag = false;
  int itoppacel = 0;
  for (int jk = klev; jk >= 1; --jk)
    if ((paph[klev + 1] - paph[jk]) < T(350.e2)) itoppacel = jk;
  for (int levels = klevm1 - 1; levels >= klev / 2 + 1; --levels) {
    for (int jk = 1; jk <= klev; ++jk) {
      plu[jk] = T(0.0);
      ptu[jk] = T(0.0);
      pqu[jk] = T(0.0);
      dh[jk] = T(0.0);
      dhen[jk] = T(0.0);
      kup[jk] = T(0.0);
      vptu[jk] = T(0.0);
      vten[jk] = T(0.0);
      abuoy[jk] = T(0.0);
      zbuo[jk] = T(0.0);
      klab[jk] = 0;
    }
    kcbot = levels;
    kctop = levels;
    zqold = T(0.);
    lldcum = false;
    bool resetflag = false;
    loflag = (!deepflag) && (levels >= itoppacel);
    for (int jk = levels; jk >= 2; --jk) {
      if (!loflag) break;
      if (jk == levels) {
        if (loflag) {
          T tmix, qmix, zmix;
          if ((paph[klev + 1] - paph[jk]) < T(60.e2)) {
            tmix = T(0.);
            qmix = T(0.);
            zmix = T(0.);
            T pmix = T(0.);
            for (int nk = jk + 2; nk >= jk; --nk) {
              if (pmix < T(50.e2)) {
                const T dp = paph[nk] - paph[nk - 1];
                tmix = tmix + dp * ptenh[nk];
                qmix = qmix + dp * pqenh[nk];
                zmix = zmix + dp * pgeoh[nk];
                pmix = pmix + dp;
              }
            }
            tmix = tmix / pmix;
            qmix = qmix / pmix;
            zmix = zmix / pmix;
          } else {
            tmix = ptenh[jk + 1];
            qmix = pqenh[jk + 1];
            zmix = pgeoh[jk + 1];
          }
          pqu[jk + 1] = qmix + deltq;
          dhen[jk + 1] = zmix + tmix * c.cpd;
          dh[jk + 1] = dhen[jk + 1] + deltt * c.cpd;
          ptu[jk + 1] = (dh[jk + 1] - pgeoh[jk + 1]) * c.rcpd;
          kup[jk + 1] = T(0.5);
          klab[jk + 1] = 1;
          vptu[jk + 1] = ptu[jk + 1] * (T(1.) + c.vtmpc1 * pqu[jk + 1]);
          vten[jk + 1] = ptenh[jk + 1] * (T(1.) + c.vtmpc1 * pqenh[jk + 1]);
          zbuo[jk + 1] = (vptu[jk + 1] - vten[jk + 1]) / vten[jk + 1];
        }
      }
      if (loflag) {
        const T r = pqsen[jk] / pqsen[levels];
        const T fscale = std::fmin(T(1.), (r * r) * r);
        eta = T(1.75e-3) * fscale;
        dz = (pgeoh[jk] - pgeoh[jk + 1]) * c.zrg;
        coef = T(0.5) * eta * dz;
        dhen[jk] = pgeoh[jk] + c.cpd * ptenh[jk];
        dh[jk] = (coef * (dhen[jk + 1] + dhen[jk]) + (T(1.) - coef) * dh[jk + 1]) / (T(1.) + coef);
        pqu[jk] = (coef * (pqenh[jk + 1] + pqenh[jk]) + (T(1.) - coef) * pqu[jk + 1]) / (T(1.) + coef);
        ptu[jk] = (dh[jk] - pgeoh[jk]) * c.rcpd;
        zqold = pqu[jk];
        zph = paph[jk];
      }
      cuadjtqn(c, jk, zph, ptu, pqu, loflag, 1);
      parcel_step(jk, true);
    }
    {
      const int ikb = kcbot, ikt = kctop;
      if (paph[ikb] - paph[ikt] < T(C::zdnoprc)) lldcum = false;
      if (lldcum) {
        ktype = 1;
        ldcum = true;
        deepflag = true;
        wbase = std::sqrt(std::fmax(T(2.) * kup[ikb], T(0.)));
        cubot = ikb;
        cutop = ikt;
        kdpl = levels + 1;
        resetflag = true;
      }
    }
    if (resetflag) {
      for (int jk = klev; jk >= 1; --jk) {
        const int ikt = kctop, ikb = kdpl;
        if (jk <= ikb && jk >= ikt) {
          culab[jk] = klab[jk];
          cutu[jk] = ptu[jk];
          cuqu[jk] = pqu[jk];
          culu[jk] = plu[jk];
        } else {
          culab[jk] = 1;
          cutu[jk] = ptenh[jk];
          cuqu[jk] = pqenh[jk];
          culu[jk] = T(0.);
        }
        if (jk < ikt) culab[jk] = 0;
      }
    }
  }
}

// cuascn: updraft (entrainment/detrainment, condensation, precipitation, mid-level convection)
template <typename T>
CPPWRF_HD inline void cuascn(const NtConst<T>& c, int klev, T* ptenh, T* pqenh, const T* puen, const T* pven,
                             const T* pten, const T* pqen, const T* pqsen, const T* pgeo, const T* pgeoh, const T* pap,
                             const T* paph, const T* pqte, const T* pverv, int klwmin, bool& ldcum, T& phcbase,
                             int& ktype, int* klab, T* ptu, T* pqu, T* plu, T* puu, T* pvu, T* pmfu, T& pmfub,
                             T* pmfus, T* pmfuq, T* pmful, T* plude, T* pdmfup, int& kcbot, int& kctop, int& kctop0,
                             int& kcum, T ztmst, T* pqsenh, T* plglac, int lndj, T& wup, T wbase, int kdpl,
                             T* pmfude_rate) {
  using C = NtConst<T>;
  (void)puen; (void)pven; (void)pqte; (void)klwmin; (void)phcbase; (void)puu; (void)pvu; (void)kcum; (void)pqsenh;
  const int klevm1 = klev - 1;
  const T zcons2 = T(3.) / (c.g * ztmst);
  const T zfacbuo = T(0.5) / (T(1.) + T(0.5));
  const T zprcdgw = T(C::cprcon) * c.zrg;
  const T z_cldmax = T(5.e-3);
  const T z_cwifrac = T(0.5);
  const T z_cprc2 = T(0.5);
  const T z_cwdrag = (T(3.0) / T(8.0)) * T(0.506) / T(0.2);
  const T cmfcmin = T(C::cmfcmin);
  bool llo3 = false;
  T zluold = T(0.);
  wup = T(0.);
  T zdpmean = T(0.);
  T zoentr = T(0.);
  if (!ldcum) {
    ktype = 0;
    kcbot = -1;
    pmfub = T(0.);
    pqu[klev] = T(0.);
  }
  T zlrain[KM], zbuo[KM], kup[KM], pdmfen[KM];
  for (int jk = 1; jk <= klev; ++jk) {
    if (jk != kcbot) plu[jk] = T(0.);
    pmfu[jk] = T(0.);
    pmfus[jk] = T(0.);
    pmfuq[jk] = T(0.);
    pmful[jk] = T(0.);
    plude[jk] = T(0.);
    plglac[jk] = T(0.);
    pdmfup[jk] = T(0.);
    zlrain[jk] = T(0.);
    zbuo[jk] = T(0.);
    kup[jk] = T(0.);
    pdmfen[jk] = T(0.);
    pmfude_rate[jk] = T(0.);
    if (!ldcum || ktype == 3) klab[jk] = 0;
    if (!ldcum && paph[jk] < T(4.e4)) kctop0 = jk;
  }
  if (ktype == 3) ldcum = false;
  kctop = kcbot;
  if (ldcum) {
    const int ikb = kcbot;
    kup[ikb] = T(0.5) * (wbase * wbase);
    pmfu[ikb] = pmfub;
    pmfus[ikb] = pmfub * (c.cpd * ptu[ikb] + pgeoh[ikb]);
    pmfuq[ikb] = pmfub * pqu[ikb];
    pmful[ikb] = pmfub * plu[ikb];
  }
  T zdmfen = T(0.), zdmfde = T(0.), zqold = T(0.);
  for (int jk = klevm1; jk >= 3; --jk) {
    cubasmcn(c, jk, pten, pqen, pqsen, pverv, pgeo, pgeoh, ldcum, ktype, klab, zlrain, pmfu, pmfub, kcbot, ptu, pqu,
             plu, pmfus, pmfuq, pmful, pdmfup);
    int is = 0;
    bool loflag = false;
    T zprecip = T(0.);
    bool llo1 = false;
    is = is + klab[jk + 1];
    if (klab[jk + 1] == 0) klab[jk] = 0;
    if ((ldcum && klab[jk + 1] == 2) || (ktype == 3 && klab[jk + 1] == 1)) loflag = true;
    const T zph = paph[jk];
    if (ktype == 3 && jk == kcbot) {
      const T zmfmax = (paph[jk] - paph[jk - 1]) * zcons2;
      if (pmfub > zmfmax) {
        const T zfac = zmfmax / pmfub;
        pmfu[jk + 1] = pmfu[jk + 1] * zfac;
        pmfus[jk + 1] = pmfus[jk + 1] * zfac;
        pmfuq[jk + 1] = pmfuq[jk + 1] * zfac;
        pmfub = zmfmax;
      }
      pmfub = std::fmin(pmfub, zmfmax);
    }
    if (is > 0) llo3 = true;
    cuentrn(c, jk, kcbot, ldcum, llo3, pgeoh, pmfu, zdmfen, zdmfde);
    if (llo3) {
      zqold = T(0.);
      if (loflag) {
        zdmfde = std::fmin(zdmfde, T(0.75) * pmfu[jk + 1]);
        if (jk == kcbot) {
          zoentr = -T(1.75e-3) * (std::fmin(T(1.), pqen[jk] / pqsen[jk]) - T(1.)) * (pgeoh[jk] - pgeoh[jk + 1]) * c.zrg;
          zoentr = std::fmin(T(0.4), zoentr) * pmfu[jk + 1];
        }
        if (jk < kcbot) {
          const T zmfmax = (paph[jk] - paph[jk - 1]) * zcons2;
          const T zxs = std::fmax(pmfu[jk + 1] - zmfmax, T(0.));
          wup = wup + kup[jk + 1] * (pap[jk + 1] - pap[jk]);
          zdpmean = zdpmean + pap[jk + 1] - pap[jk];
          zdmfen = zoentr;
          if (ktype >= 2) {
            zdmfen = T(2.0) * zdmfen;
            zdmfde = zdmfen;
          }
          zdmfde = zdmfde * (T(1.6) - std::fmin(T(1.), pqen[jk] / pqsen[jk]));
          const T zmftest = pmfu[jk + 1] + zdmfen - zdmfde;
          T zchange = std::fmax(zmftest - zmfmax, T(0.));
          const T zxe = std::fmax(zchange - zxs, T(0.));
          zdmfen = zdmfen - zxe;
          zchange = zchange - zxe;
          zdmfde = zdmfde + zchange;
        }
        pdmfen[jk] = zdmfen - zdmfde;
        pmfu[jk] = pmfu[jk + 1] + zdmfen - zdmfde;
        const T zqeen = pqenh[jk + 1] * zdmfen;
        const T zseen = (c.cpd * ptenh[jk + 1] + pgeoh[jk + 1]) * zdmfen;
        const T zscde = (c.cpd * ptu[jk + 1] + pgeoh[jk + 1]) * zdmfde;
        const T zqude = pqu[jk + 1] * zdmfde;
        plude[jk] = plu[jk + 1] * zdmfde;
        const T zmfusk = pmfus[jk + 1] + zseen - zscde;
        const T zmfuqk = pmfuq[jk + 1] + zqeen - zqude;
        const T zmfulk = pmful[jk + 1] - plude[jk];
        plu[jk] = zmfulk * (T(1.) / std::fmax(cmfcmin, pmfu[jk]));
        pqu[jk] = zmfuqk * (T(1.) / std::fmax(cmfcmin, pmfu[jk]));
        ptu[jk] = (zmfusk * (T(1.) / std::fmax(cmfcmin, pmfu[jk])) - pgeoh[jk]) * c.rcpd;
        ptu[jk] = std::fmax(T(100.), ptu[jk]);
        ptu[jk] = std::fmin(T(400.), ptu[jk]);
        zqold = pqu[jk];
        zlrain[jk] = zlrain[jk + 1] * (pmfu[jk + 1] - zdmfde) * (T(1.) / std::fmax(cmfcmin, pmfu[jk]));
        zluold = plu[jk];
      }
      if (jk > kdpl) {
        ptu[jk] = ptenh[jk];
        pqu[jk] = pqenh[jk];
        plu[jk] = T(0.);
        zluold = plu[jk];
      }
      if (loflag) cuadjtqn(c, jk, zph, ptu, pqu, loflag, 1);
      if (loflag) {
        if (pqu[jk] != zqold) {
          plglac[jk] = plu[jk] * ((T(1.) - foealfa(ptu[jk])) - (T(1.) - foealfa(ptu[jk + 1])));
          ptu[jk] = ptu[jk] + c.ralfdcp * plglac[jk];
        }
      }
      if (loflag) {
        if (pqu[jk] != zqold) {
          klab[jk] = 2;
          plu[jk] = plu[jk] + zqold - pqu[jk];
          const T zbc = ptu[jk] * (T(1.) + c.vtmpc1 * pqu[jk] - plu[jk + 1] - zlrain[jk + 1]);
          const T zbe = ptenh[jk] * (T(1.) + c.vtmpc1 * pqenh[jk]);
          zbuo[jk] = zbc - zbe;
          if (ktype == 3 && klab[jk + 1] == 1) {
            if (zbuo[jk] > -T(0.5)) {
              ldcum = true;
              kctop = jk;
              kup[jk] = T(0.5);
            } else {
              klab[jk] = 0;
              pmfu[jk] = T(0.);
              plude[jk] = T(0.);
              plu[jk] = T(0.);
            }
          }
          if (klab[jk + 1] == 2) {
            if (zbuo[jk] < T(0.)) {
              ptenh[jk] = T(0.5) * (pten[jk] + pten[jk - 1]);
              pqenh[jk] = T(0.5) * (pqen[jk] + pqen[jk - 1]);
              zbuo[jk] = zbc - ptenh[jk] * (T(1.) + c.vtmpc1 * pqenh[jk]);
            }
            const T zbuoc = (zbuo[jk] / (ptenh[jk] * (T(1.) + c.vtmpc1 * pqenh[jk])) +
                             zbuo[jk + 1] / (ptenh[jk + 1] * (T(1.) + c.vtmpc1 * pqenh[jk + 1]))) *
                            T(0.5);
            const T zdkbuo = (pgeoh[jk] - pgeoh[jk + 1]) * zfacbuo * zbuoc;
            T zdken;
            if (zdmfen > T(0.)) {
              zdken = std::fmin(T(1.), (T(1.) + z_cwdrag) * zdmfen / std::fmax(cmfcmin, pmfu[jk + 1]));
            } else {
              zdken = std::fmin(T(1.), (T(1.) + z_cwdrag) * zdmfde / std::fmax(cmfcmin, pmfu[jk + 1]));
            }
            kup[jk] = (kup[jk + 1] * (T(1.) - zdken) + zdkbuo) / (T(1.) + zdken);
            if (zbuo[jk] < T(0.)) {
              T zkedke = kup[jk] / std::fmax(T(1.e-10), kup[jk + 1]);
              zkedke = std::fmax(T(0.), std::fmin(T(1.), zkedke));
              const T zmfun = std::sqrt(zkedke) * pmfu[jk + 1];
              zdmfde = std::fmax(zdmfde, pmfu[jk + 1] - zmfun);
              plude[jk] = plu[jk + 1] * zdmfde;
              pmfu[jk] = pmfu[jk + 1] + zdmfen - zdmfde;
            }
            if (zbuo[jk] > -T(0.2)) {
              const int ikb = kcbot;
              const T m = std::fmin(T(1.), pqsen[jk] / pqsen[ikb]);
              zoentr = T(1.75e-3) * (T(0.3) - (std::fmin(T(1.), pqen[jk - 1] / pqsen[jk - 1]) - T(1.))) *
                       (pgeoh[jk - 1] - pgeoh[jk]) * c.zrg * ((m * m) * m);
              zoentr = std::fmin(T(0.4), zoentr) * pmfu[jk];
            } else {
              zoentr = T(0.);
            }
            if (jk > kdpl) {
              pmfu[jk] = pmfu[jk + 1];
              kup[jk] = T(0.5);
            }
            if (kup[jk] > T(0.) && pmfu[jk] > T(0.)) {
              kctop = jk;
              llo1 = true;
            } else {
              klab[jk] = 0;
              pmfu[jk] = T(0.);
              kup[jk] = T(0.);
              zdmfde = pmfu[jk + 1];
              plude[jk] = plu[jk + 1] * zdmfde;
            }
            if (pmfu[jk + 1] > T(0.)) pmfude_rate[jk] = zdmfde;
          }
        } else if (ktype == 2 && pqu[jk] == zqold) {
          klab[jk] = 0;
          pmfu[jk] = T(0.);
          kup[jk] = T(0.);
          zdmfde = pmfu[jk + 1];
          plude[jk] = plu[jk + 1] * zdmfde;
          pmfude_rate[jk] = zdmfde;
        }
      }
      if (llo1) {
        T zdshrd;
        if (lndj == 1) zdshrd = T(5.e-4);
        else zdshrd = T(3.e-4);
        if (plu[jk] > zdshrd) {
          const T zwu = std::fmin(T(15.0), std::sqrt(T(2.) * std::fmax(T(0.1), kup[jk + 1])));
          const T zprcon = zprcdgw / (T(0.75) * zwu);
          const T zdt = std::fmin(T(C::rtber) - T(C::rtice), std::fmax(T(C::rtber) - ptu[jk], T(0.)));
          const T zcbf = T(1.) + z_cprc2 * std::sqrt(zdt);
          const T zzco = zprcon * zcbf;
          const T zlcrit = zdshrd / zcbf;
          const T zdfi = pgeoh[jk] - pgeoh[jk + 1];
          const T zc = (plu[jk] - zluold);
          const T q = plu[jk] / zlcrit;
          const T zarg = q * q;
          T zd;
          if (zarg < T(25.0)) zd = zzco * (T(1.) - std::exp(-zarg)) * zdfi;
          else zd = zzco * zdfi;
          const T zint = std::exp(-zd);
          T zlnew = zluold * zint + zc / zd * (T(1.) - zint);
          zlnew = std::fmax(T(0.), std::fmin(plu[jk], zlnew));
          zlnew = std::fmin(z_cldmax, zlnew);
          zprecip = std::fmax(T(0.), zluold + zc - zlnew);
          pdmfup[jk] = zprecip * pmfu[jk];
          zlrain[jk] = zlrain[jk] + zprecip;
          plu[jk] = zlnew;
        }
      }
      if (llo1) {
        if (zlrain[jk] > T(0.)) {
          const T zvw = T(21.18) * std::pow(zlrain[jk], T(0.2));
          const T zvi = z_cwifrac * zvw;
          const T zalfaw = foealfa(ptu[jk]);
          const T zvv = zalfaw * zvw + (T(1.) - zalfaw) * zvi;
          const T zrold = zlrain[jk] - zprecip;
          const T zc = zprecip;
          const T zwu = std::fmin(T(15.0), std::sqrt(T(2.) * std::fmax(T(0.1), kup[jk])));
          const T zd = zvv / zwu;
          const T zint = std::exp(-zd);
          T zrnew = zrold * zint + zc / zd * (T(1.) - zint);
          zrnew = std::fmax(T(0.), std::fmin(zlrain[jk], zrnew));
          zlrain[jk] = zrnew;
        }
      }
      if (loflag) {
        pmful[jk] = plu[jk] * pmfu[jk];
        pmfus[jk] = (c.cpd * ptu[jk] + pgeoh[jk]) * pmfu[jk];
        pmfuq[jk] = pqu[jk] * pmfu[jk];
      }
    }
  }
  if (kctop == -1) ldcum = false;
  kcbot = kcbot > kctop ? kcbot : kctop;
  if (ldcum) {
    wup = std::fmax(T(1.e-2), wup / std::fmax(T(1.), zdpmean));
    wup = std::sqrt(T(2.) * wup);
  }
}

// cudlfsn: level of free sinking (downdraft top)
template <typename T>
CPPWRF_HD inline void cudlfsn(const NtConst<T>& c, int klev, int kcbot, int kctop, int lndj, bool ldcum,
                              const T* ptenh, const T* pqenh, const T* puen, const T* pven, const T* pten,
                              const T* pqsen, const T* pgeo, const T* pgeoh, const T* paph, const T* ptu, const T* pqu,
                              const T* plu, const T* puu, const T* pvu, T pmfub, T& prfl, T* ptd, T* pqd, T* pud,
                              T* pvd, T* pmfd, T* pmfds, T* pmfdq, T* pdmfdp, int& kdtop, bool& lddraf) {
  using C = NtConst<T>;
  (void)lndj; (void)puen; (void)pven; (void)plu; (void)puu; (void)pvu; (void)pud; (void)pvd;
  lddraf = false;
  kdtop = klev + 1;
  int ikhsmin = klev + 1;
  T zhsmin = T(1.e8);
  for (int jk = 3; jk <= klev - 2; ++jk) {
    const T zhsk = c.cpd * pten[jk] + pgeo[jk] + foelhm(c, pten[jk]) * pqsen[jk];
    if (zhsk < zhsmin) {
      zhsmin = zhsk;
      ikhsmin = jk;
    }
  }
  const int ike = klev - 3;
  T ztenwb[KM], zqenwb[KM];
  for (int jk = 3; jk <= ike; ++jk) {
    ztenwb[jk] = ptenh[jk];
    zqenwb[jk] = pqenh[jk];
    const T zph = paph[jk];
    const bool llo2 = ldcum && prfl > T(0.) && !lddraf && (jk < kcbot && jk > kctop) && jk >= ikhsmin;
    if (!llo2) continue;
    cuadjtqn(c, jk, zph, ztenwb, zqenwb, llo2, 2);
    const T zttest = T(0.5) * (ptu[jk] + ztenwb[jk]);
    const T zqtest = T(0.5) * (pqu[jk] + zqenwb[jk]);
    const T zbuo = zttest * (T(1.) + c.vtmpc1 * zqtest) - ptenh[jk] * (T(1.) + c.vtmpc1 * pqenh[jk]);
    const T zcond = pqenh[jk] - zqenwb[jk];
    const T zmftop = -T(C::cmfdeps) * pmfub;
    if (zbuo < T(0.) && prfl > T(10.) * zmftop * zcond) {
      kdtop = jk;
      lddraf = true;
      ptd[jk] = zttest;
      pqd[jk] = zqtest;
      pmfd[jk] = zmftop;
      pmfds[jk] = pmfd[jk] * (c.cpd * ptd[jk] + pgeoh[jk]);
      pmfdq[jk] = pmfd[jk] * pqd[jk];
      pdmfdp[jk - 1] = -T(0.5) * pmfd[jk] * zcond;
      prfl = prfl + pdmfdp[jk - 1];
    }
  }
}

// cuddrafn: downdraft
template <typename T>
CPPWRF_HD inline void cuddrafn(const NtConst<T>& c, int klev, bool lddraf, const T* ptenh, const T* pqenh,
                               const T* puen, const T* pven, const T* pgeo, const T* pgeoh, const T* paph, T& prfl,
                               T* ptd, T* pqd, T* pud, T* pvd, const T* pmfu, T* pmfd, T* pmfds, T* pmfdq, T* pdmfdp,
                               T* pmfdde_rate) {
  using C = NtConst<T>;
  (void)puen; (void)pven; (void)pud; (void)pvd;
  T zoentr = T(0.), zbuoy = T(0.), zdmfen = T(0.), zdmfde = T(0.);
  int itopde = klev;
  for (int jk = klev; jk >= 1; --jk) {
    pmfdde_rate[jk] = T(0.);
    if ((paph[klev + 1] - paph[jk]) < T(60.e2)) itopde = jk;
  }
  for (int jk = 3; jk <= klev; ++jk) {
    const T zph = paph[jk];
    const bool llo2 = lddraf && pmfd[jk - 1] < T(0.);
    if (!llo2) continue;
    {
      const T zentr = T(C::entrdd) * pmfd[jk - 1] * (pgeoh[jk - 1] - pgeoh[jk]) * c.zrg;
      zdmfen = zentr;
      zdmfde = zentr;
    }
    if (jk > itopde) {
      zdmfen = T(0.);
      zdmfde = pmfd[itopde] * (paph[jk] - paph[jk - 1]) / (paph[klev + 1] - paph[itopde]);
    }
    if (jk <= itopde) {
      const T zdz = -(pgeoh[jk - 1] - pgeoh[jk]) * c.zrg;
      const T zzentr = zoentr * zdz * pmfd[jk - 1];
      zdmfen = zdmfen + zzentr;
      zdmfen = std::fmax(zdmfen, T(0.3) * pmfd[jk - 1]);
      zdmfen = std::fmax(zdmfen, -T(0.75) * pmfu[jk] - (pmfd[jk - 1] - zdmfde));
      zdmfen = std::fmin(zdmfen, T(0.));
    }
    pmfd[jk] = pmfd[jk - 1] + zdmfen - zdmfde;
    const T zseen = (c.cpd * ptenh[jk - 1] + pgeoh[jk - 1]) * zdmfen;
    const T zqeen = pqenh[jk - 1] * zdmfen;
    const T zsdde = (c.cpd * ptd[jk - 1] + pgeoh[jk - 1]) * zdmfde;
    const T zqdde = pqd[jk - 1] * zdmfde;
    const T zmfdsk = pmfds[jk - 1] + zseen - zsdde;
    const T zmfdqk = pmfdq[jk - 1] + zqeen - zqdde;
    pqd[jk] = zmfdqk * (T(1.) / std::fmin(-T(C::cmfcmin), pmfd[jk]));
    ptd[jk] = (zmfdsk * (T(1.) / std::fmin(-T(C::cmfcmin), pmfd[jk])) - pgeoh[jk]) * c.rcpd;
    ptd[jk] = std::fmin(T(400.), ptd[jk]);
    ptd[jk] = std::fmax(T(100.), ptd[jk]);
    T zcond = pqd[jk];
    cuadjtqn(c, jk, zph, ptd, pqd, llo2, 2);
    zcond = zcond - pqd[jk];
    T zbuo = ptd[jk] * (T(1.) + c.vtmpc1 * pqd[jk]) - ptenh[jk] * (T(1.) + c.vtmpc1 * pqenh[jk]);
    if (prfl > T(0.) && pmfu[jk] > T(0.)) {
      const T zrain = prfl / pmfu[jk];
      zbuo = zbuo - ptd[jk] * zrain;
    }
    if (zbuo >= T(0) || prfl <= (pmfd[jk] * zcond)) {
      pmfd[jk] = T(0.);
      zbuo = T(0.);
    }
    pmfds[jk] = (c.cpd * ptd[jk] + pgeoh[jk]) * pmfd[jk];
    pmfdq[jk] = pqd[jk] * pmfd[jk];
    const T zdmfdp = -pmfd[jk] * zcond;
    pdmfdp[jk - 1] = zdmfdp;
    prfl = prfl + zdmfdp;
    T zbuoyz = zbuo / ptenh[jk];
    zbuoyz = std::fmin(zbuoyz, T(0.0));
    const T zdz = -(pgeo[jk - 1] - pgeo[jk]);
    zbuoy = zbuoy + zbuoyz * zdz;
    zoentr = c.g * zbuoyz * T(0.5) / (T(1.) + zbuoy);
    pmfdde_rate[jk] = -zdmfde;
  }
}

// cuflxn: final convective fluxes, precipitation fluxes, melting and evaporation below cloud base
template <typename T>
CPPWRF_HD inline void cuflxn(const NtConst<T>& c, int klev, T ztmst, const T* pten, const T* pqen, T* pqsen,
                             const T* ptenh, const T* pqenh, const T* paph, const T* pap, const T* pgeoh, int lndj,
                             bool ldcum, int kcbot, int kctop, int kdtop, int& ktopm2, int& ktype, bool& lddraf,
                             T* pmfu, T* pmfd, T* pmfus, T* pmfds, T* pmfuq, T* pmfdq, T* pmful, T* plude, T* pdmfup,
                             T* pdmfdp, T* pdpmel, T* plglac, T& prain, T* pmfdde_rate, T* pmflxr, T* pmflxs) {
  using C = NtConst<T>;
  const T tmelt = T(C::tmelt);
  const T ztaumel = T(18000.);
  const T zcons1a = c.cpd / (c.alf * c.g * ztaumel);
  const T zcons2 = T(3.) / (c.g * ztmst);
  const T zcucov = T(0.05);
  const T zcpecons = T(5.44e-4) / c.g;
  prain = T(0.);
  if (!ldcum || kdtop < kctop) lddraf = false;
  if (!ldcum) ktype = 0;
  int idbas = klev;
  T rhevap;
  if (lndj == 1) rhevap = T(0.7);
  else rhevap = T(0.9);
  ktopm2 = 2;
  for (int jk = ktopm2; jk <= klev; ++jk) {
    const int ikb = (jk + 1 < klev) ? jk + 1 : klev;
    pmflxr[jk] = T(0.);
    pmflxs[jk] = T(0.);
    pdpmel[jk] = T(0.);
    if (ldcum && jk >= kctop) {
      pmfus[jk] = pmfus[jk] - pmfu[jk] * (c.cpd * ptenh[jk] + pgeoh[jk]);
      pmfuq[jk] = pmfuq[jk] - pmfu[jk] * pqenh[jk];
      plglac[jk] = pmfu[jk] * plglac[jk];
      const bool llddraf = lddraf && jk >= kdtop;
      if (llddraf && jk >= kdtop) {
        pmfds[jk] = pmfds[jk] - pmfd[jk] * (c.cpd * ptenh[jk] + pgeoh[jk]);
        pmfdq[jk] = pmfdq[jk] - pmfd[jk] * pqenh[jk];
      } else {
        pmfd[jk] = T(0.);
        pmfds[jk] = T(0.);
        pmfdq[jk] = T(0.);
        pdmfdp[jk - 1] = T(0.);
      }
      if (llddraf && pmfd[jk] < T(0.) && std::fabs(pmfd[ikb]) < T(1.e-20)) idbas = jk;
    } else {
      pmfu[jk] = T(0.);
      pmfd[jk] = T(0.);
      pmfus[jk] = T(0.);
      pmfds[jk] = T(0.);
      pmfuq[jk] = T(0.);
      pmfdq[jk] = T(0.);
      pmful[jk] = T(0.);
      plglac[jk] = T(0.);
      pdmfup[jk - 1] = T(0.);
      pdmfdp[jk - 1] = T(0.);
      plude[jk - 1] = T(0.);
    }
  }
  pmflxr[klev + 1] = T(0.);
  pmflxs[klev + 1] = T(0.);
  if (ldcum) {
    const int ikb = kcbot;
    const int ik = ikb + 1;
    T zzp = ((paph[klev + 1] - paph[ik]) / (paph[klev + 1] - paph[ikb]));
    if (ktype == 3) zzp = zzp * zzp;
    pmfu[ik] = pmfu[ikb] * zzp;
    pmfus[ik] = (pmfus[ikb] - foelhm(c, ptenh[ikb]) * pmful[ikb]) * zzp;
    pmfuq[ik] = (pmfuq[ikb] + pmful[ikb]) * zzp;
    pmful[ik] = T(0.);
  }
  for (int jk = ktopm2; jk <= klev; ++jk) {
    if (ldcum && jk > kcbot + 1) {
      const int ikb = kcbot + 1;
      T zzp = ((paph[klev + 1] - paph[jk]) / (paph[klev + 1] - paph[ikb]));
      if (ktype == 3) zzp = zzp * zzp;
      pmfu[jk] = pmfu[ikb] * zzp;
      pmfus[jk] = pmfus[ikb] * zzp;
      pmfuq[jk] = pmfuq[ikb] * zzp;
      pmful[jk] = T(0.);
    }
    const int ik = idbas;
    const bool llddraf = lddraf && jk > ik && ik < klev;
    if (llddraf && ik == kcbot + 1) {
      T zzp = ((paph[klev + 1] - paph[jk]) / (paph[klev + 1] - paph[ik]));
      if (ktype == 3) zzp = zzp * zzp;
      pmfd[jk] = pmfd[ik] * zzp;
      pmfds[jk] = pmfds[ik] * zzp;
      pmfdq[jk] = pmfdq[ik] * zzp;
      pmfdde_rate[jk] = -(pmfd[jk - 1] - pmfd[jk]);
    } else if (llddraf && ik != kcbot + 1 && jk == ik + 1) {
      pmfdde_rate[jk] = -(pmfd[jk - 1] - pmfd[jk]);
    }
  }
  for (int jk = ktopm2; jk <= klev; ++jk) {
    if (ldcum && jk >= kctop - 1) {
      prain = prain + pdmfup[jk];
      if (pmflxs[jk] > T(0.) && pten[jk] > tmelt) {
        const T zcons1 = zcons1a * (T(1.) + T(0.5) * (pten[jk] - tmelt));
        const T zfac = zcons1 * (paph[jk + 1] - paph[jk]);
        const T zsnmlt = std::fmin(pmflxs[jk], zfac * (pten[jk] - tmelt));
        pdpmel[jk] = zsnmlt;
        pqsen[jk] = foeewm(c, pten[jk] - zsnmlt / zfac) / pap[jk];
      }
      T zalfaw = foealfa(pten[jk]);
      if (pten[jk] < tmelt && zalfaw > T(0.)) {
        plglac[jk] = plglac[jk] + zalfaw * (pdmfup[jk] + pdmfdp[jk]);
        zalfaw = T(0.);
      }
      pmflxr[jk + 1] = pmflxr[jk] + zalfaw * (pdmfup[jk] + pdmfdp[jk]) + pdpmel[jk];
      pmflxs[jk + 1] = pmflxs[jk] + (T(1.) - zalfaw) * (pdmfup[jk] + pdmfdp[jk]) - pdpmel[jk];
      if (pmflxr[jk + 1] + pmflxs[jk + 1] < T(0.0)) {
        pdmfdp[jk] = -(pmflxr[jk] + pmflxs[jk] + pdmfup[jk]);
        pmflxr[jk + 1] = T(0.0);
        pmflxs[jk + 1] = T(0.0);
        pdpmel[jk] = T(0.0);
      } else if (pmflxr[jk + 1] < T(0.)) {
        pmflxs[jk + 1] = pmflxs[jk + 1] + pmflxr[jk + 1];
        pmflxr[jk + 1] = T(0.);
      } else if (pmflxs[jk + 1] < T(0.)) {
        pmflxr[jk + 1] = pmflxr[jk + 1] + pmflxs[jk + 1];
        pmflxs[jk + 1] = T(0.);
      }
    }
  }
  for (int jk = ktopm2; jk <= klev; ++jk) {
    if (ldcum && jk >= kcbot) {
      const T zrfl = pmflxr[jk] + pmflxs[jk];
      if (zrfl > T(1.e-20)) {
        const T zdrfl1 = zcpecons * std::fmax(T(0.), pqsen[jk] - pqen[jk]) * zcucov *
                         std::pow(std::sqrt(paph[jk] / paph[klev + 1]) / T(5.09e-3) * zrfl / zcucov, T(0.5777)) *
                         (paph[jk + 1] - paph[jk]);
        T zrnew = zrfl - zdrfl1;
        const T zrmin = zrfl - zcucov * std::fmax(T(0.), rhevap * pqsen[jk] - pqen[jk]) * zcons2 *
                                   (paph[jk + 1] - paph[jk]);
        zrnew = std::fmax(zrnew, zrmin);
        const T zrfln = std::fmax(zrnew, T(0.));
        const T zdrfl = std::fmin(T(0.), zrfln - zrfl);
        const T zdenom = T(1.) / std::fmax(T(1.e-20), pmflxr[jk] + pmflxs[jk]);
        T zalfaw = foealfa(pten[jk]);
        if (pten[jk] < tmelt) zalfaw = T(0.);
        const T zpdr = zalfaw * pdmfdp[jk];
        const T zpds = (T(1.) - zalfaw) * pdmfdp[jk];
        pmflxr[jk + 1] = pmflxr[jk] + zpdr + pdpmel[jk] + zdrfl * pmflxr[jk] * zdenom;
        pmflxs[jk + 1] = pmflxs[jk] + zpds - pdpmel[jk] + zdrfl * pmflxs[jk] * zdenom;
        pdmfup[jk] = pdmfup[jk] + zdrfl;
        if (pmflxr[jk + 1] + pmflxs[jk + 1] < T(0.)) {
          pdmfup[jk] = pdmfup[jk] - (pmflxr[jk + 1] + pmflxs[jk + 1]);
          pmflxr[jk + 1] = T(0.);
          pmflxs[jk + 1] = T(0.);
          pdpmel[jk] = T(0.);
        } else if (pmflxr[jk + 1] < T(0.)) {
          pmflxs[jk + 1] = pmflxs[jk + 1] + pmflxr[jk + 1];
          pmflxr[jk + 1] = T(0.);
        } else if (pmflxs[jk + 1] < T(0.)) {
          pmflxr[jk + 1] = pmflxr[jk + 1] + pmflxs[jk + 1];
          pmflxs[jk + 1] = T(0.);
        }
      } else {
        pmflxr[jk + 1] = T(0.0);
        pmflxs[jk + 1] = T(0.0);
        pdmfdp[jk] = T(0.0);
        pdpmel[jk] = T(0.0);
      }
    }
  }
}

// cudtdqn: temperature and moisture tendencies, detrained condensate
template <typename T>
CPPWRF_HD inline void cudtdqn(const NtConst<T>& c, int klev, int ktopm2, bool ldcum, const T* paph, const T* pten,
                              const T* plglac, const T* plude, const T* pmfu, const T* pmfd, const T* pmfus,
                              const T* pmfds, const T* pmfuq, const T* pmfdq, const T* pmful, const T* pdmfup,
                              const T* pdmfdp, const T* pdpmel, T* ptent, T* ptenq, T* pcte) {
  (void)pmfu; (void)pmfd;
  if (!ldcum) return;
  T zdp[KM], zdtdt[KM], zdqdt[KM];
  for (int jk = 1; jk <= klev; ++jk) zdp[jk] = c.g / (paph[jk + 1] - paph[jk]);
  for (int jk = ktopm2; jk <= klev; ++jk) {
    const T zalv = foelhm(c, pten[jk]);
    if (jk < klev) {
      zdtdt[jk] = zdp[jk] * c.rcpd *
                  (pmfus[jk + 1] - pmfus[jk] + pmfds[jk + 1] - pmfds[jk] + c.alf * plglac[jk] - c.alf * pdpmel[jk] -
                   zalv * (pmful[jk + 1] - pmful[jk] - plude[jk] - pdmfup[jk] - pdmfdp[jk]));
      zdqdt[jk] = zdp[jk] * (pmfuq[jk + 1] - pmfuq[jk] + pmfdq[jk + 1] - pmfdq[jk] + pmful[jk + 1] - pmful[jk] -
                             plude[jk] - pdmfup[jk] - pdmfdp[jk]);
    } else {
      zdtdt[jk] = -zdp[jk] * c.rcpd *
                  (pmfus[jk] + pmfds[jk] + c.alf * pdpmel[jk] - zalv * (pmful[jk] + pdmfup[jk] + pdmfdp[jk] + plude[jk]));
      zdqdt[jk] = -zdp[jk] * (pmfuq[jk] + plude[jk] + pmfdq[jk] + (pmful[jk] + pdmfup[jk] + pdmfdp[jk]));
    }
  }
  for (int jk = ktopm2; jk <= klev; ++jk) {
    ptent[jk] = ptent[jk] + zdtdt[jk];
    ptenq[jk] = ptenq[jk] + zdqdt[jk];
    pcte[jk] = zdp[jk] * plude[jk];
  }
}

// cududvn: momentum tendencies
template <typename T>
CPPWRF_HD inline void cududvn(const NtConst<T>& c, int klev, int ktopm2, int ktype, int kcbot, bool ldcum,
                              const T* paph, const T* puen, const T* pven, const T* pmfu, const T* pmfd, const T* puu,
                              const T* pud, const T* pvu, const T* pvd, T* ptenu, T* ptenv) {
  if (!ldcum) return;
  T zdp[KM], zuen[KM], zven[KM], zmfuu[KM], zmfuv[KM], zmfdu[KM], zmfdv[KM], zdudt[KM], zdvdt[KM];
  for (int jk = 1; jk <= klev; ++jk) {
    zuen[jk] = puen[jk];
    zven[jk] = pven[jk];
    zdp[jk] = c.g / (paph[jk + 1] - paph[jk]);
  }
  for (int jk = ktopm2; jk <= klev; ++jk) {
    const int ik = jk - 1;
    zmfuu[jk] = pmfu[jk] * (puu[jk] - zuen[ik]);
    zmfuv[jk] = pmfu[jk] * (pvu[jk] - zven[ik]);
    zmfdu[jk] = pmfd[jk] * (pud[jk] - zuen[ik]);
    zmfdv[jk] = pmfd[jk] * (pvd[jk] - zven[ik]);
  }
  for (int jk = ktopm2; jk <= klev; ++jk) {
    if (jk > kcbot) {
      const int ikb = kcbot;
      T zzp = ((paph[klev + 1] - paph[jk]) / (paph[klev + 1] - paph[ikb]));
      if (ktype == 3) zzp = zzp * zzp;
      zmfuu[jk] = zmfuu[ikb] * zzp;
      zmfuv[jk] = zmfuv[ikb] * zzp;
      zmfdu[jk] = zmfdu[ikb] * zzp;
      zmfdv[jk] = zmfdv[ikb] * zzp;
    }
  }
  for (int jk = ktopm2; jk <= klev; ++jk) {
    if (jk < klev) {
      const int ik = jk + 1;
      zdudt[jk] = zdp[jk] * (zmfuu[ik] - zmfuu[jk] + zmfdu[ik] - zmfdu[jk]);
      zdvdt[jk] = zdp[jk] * (zmfuv[ik] - zmfuv[jk] + zmfdv[ik] - zmfdv[jk]);
    } else {
      zdudt[jk] = -zdp[jk] * (zmfuu[jk] + zmfdu[jk]);
      zdvdt[jk] = -zdp[jk] * (zmfuv[jk] + zmfdv[jk]);
    }
  }
  for (int jk = ktopm2; jk <= klev; ++jk) {
    ptenu[jk] = ptenu[jk] + zdudt[jk];
    ptenv[jk] = ptenv[jk] + zdvdt[jk];
  }
}

// cumastrn: master routine (one column)
template <typename T>
CPPWRF_HD inline void cumastrn(const NtConst<T>& c, int klev, const T* pten, const T* pqen, const T* puen,
                               const T* pven, const T* pverv, T* pqsen, T pqhfl, T ztmst, const T* pap, const T* paph,
                               const T* pgeo, T* ptte, T* pqte, T* pvom, T* pvol, T& prsfc, T& pssfc, bool& ldcum,
                               int& ktype, int& kcbot, int& kctop, T* ptu, T* pqu, T* plu, T* plude, T* pmfu,
                               T* pmfd, T& prain, T* pcte, T phhfl, int lndj, const T* zgeoh, T dx, T scale_fac,
                               T scale_fac2) {
  using C = NtConst<T>;
  (void)dx;
  const T cmfcmin = T(C::cmfcmin);
  const T zcons = T(1.) / (c.g * ztmst);
  const T zcons2 = T(3.) / (c.g * ztmst);
  int ilwmin = 0, ictop0 = 0, kdpl = 0, idtop = 0, icum = 0, itopm2 = 0;
  int ilab[KM];
  T ztenh[KM], zqenh[KM], zqsenh[KM], ztd[KM], zqd[KM], zmfus[KM], zmfds[KM], zmfuq[KM], zmfdq[KM], zdmfup[KM],
      zdmfdp[KM], zmful[KM], zuu[KM], zvu[KM], zud[KM], zvd[KM], zlglac[KM], zdpmel[KM];
  T pmfude_rate[KM], pmfdde_rate[KM], zmfuus[KM], zmfdus[KM], zuv2[KM], ztenu[KM], ztenv[KM];
  T pmflxr[KM], pmflxs[KM];
  T wbase = T(0), wup = T(0), zhcbase = T(0), zmfub = T(0), zmfub1 = T(0), zdhpbl = T(0), upbl = T(0), zrfl = T(0);
  T zmfs = T(0), zmfuub = T(0), zheat = T(0), zcape = T(0), zcape1 = T(0), zcape2 = T(0), ztauc = T(0),
    ztaubl = T(0), zsum12 = T(0), zsum22 = T(0);
  bool loddraf = false, llo2 = false;
  for (int k = 0; k < KM; ++k) {  // Fortran locals are uninitialised; zero them (only written values are read)
    pmfude_rate[k] = T(0); pmfdde_rate[k] = T(0); zlglac[k] = T(0); pmflxr[k] = T(0); pmflxs[k] = T(0);
  }

  cuinin(c, klev, pten, pqen, pqsen, puen, pven, pverv, pgeo, paph, zgeoh, ztenh, zqenh, zqsenh, ilwmin, ptu, pqu, ztd,
         zqd, zuu, zvu, zud, zvd, pmfu, pmfd, zmfus, zmfds, zmfuq, zmfdq, zdmfup, zdmfdp, zdpmel, plu, plude, ilab);
  cutypen(c, klev, pqen, ztenh, zqenh, zqsenh, zgeoh, paph, phhfl, pqhfl, pgeo, pqsen, pap, pten, lndj, ptu, pqu, ilab,
          ldcum, kcbot, ictop0, ktype, wbase, plu, kdpl);

  zdhpbl = T(0.0);
  upbl = T(0.0);
  idtop = 0;
  for (int jk = 2; jk <= klev; ++jk) {
    if (jk >= kcbot && ldcum) {
      zdhpbl = zdhpbl + (c.alv * pqte[jk] + c.cpd * ptte[jk]) * (paph[jk + 1] - paph[jk]);
      if (lndj == 0) {
        const T wspeed = std::sqrt(puen[jk] * puen[jk] + pven[jk] * pven[jk]);
        upbl = upbl + wspeed * (paph[jk + 1] - paph[jk]);
      }
    }
  }
  if (ldcum) {
    const int ikb = kcbot;
    const T zmfmax = (paph[ikb] - paph[ikb - 1]) * zcons2;
    if (ktype == 1) {
      zmfub = T(0.1) * zmfmax;
    } else if (ktype == 2) {
      const T zqumqe = pqu[ikb] + plu[ikb] - zqenh[ikb];
      const T zdqmin = std::fmax(T(0.01) * zqenh[ikb], T(1.e-10));
      T zdh = c.cpd * (ptu[ikb] - ztenh[ikb]) + c.alv * zqumqe;
      zdh = c.g * std::fmax(zdh, T(1.e5) * zdqmin);
      if (zdhpbl > T(0.)) {
        zmfub = zdhpbl / zdh;
        zmfub = std::fmin(zmfub, zmfmax);
      } else {
        zmfub = T(0.1) * zmfmax;
        ldcum = false;
      }
    }
  } else {
    zmfub = T(0.);
  }

  cuascn(c, klev, ztenh, zqenh, puen, pven, pten, pqen, pqsen, pgeo, zgeoh, pap, paph, pqte, pverv, ilwmin, ldcum,
         zhcbase, ktype, ilab, ptu, pqu, plu, zuu, zvu, pmfu, zmfub, zmfus, zmfuq, zmful, plude, zdmfup, kcbot, kctop,
         ictop0, icum, ztmst, zqsenh, zlglac, lndj, wup, wbase, kdpl, pmfude_rate);

  if (ldcum) {
    const int ikb = kcbot;
    itopm2 = kctop;
    const T zpbmpt = paph[ikb] - paph[itopm2];
    if (ktype == 1 && zpbmpt < T(C::zdnoprc)) ktype = 2;
    if (ktype == 2 && zpbmpt >= T(C::zdnoprc)) ktype = 1;
    ictop0 = kctop;
  }
  zrfl = zdmfup[1];
  for (int jk = 2; jk <= klev; ++jk) zrfl = zrfl + zdmfup[jk];
  for (int jk = 1; jk <= klev; ++jk) {
    pmfd[jk] = T(0.);
    zmfds[jk] = T(0.);
    zmfdq[jk] = T(0.);
    zdmfdp[jk] = T(0.);
    zdpmel[jk] = T(0.);
  }
  if (C::lmfdd) {
    cudlfsn(c, klev, kcbot, kctop, lndj, ldcum, ztenh, zqenh, puen, pven, pten, pqsen, pgeo, zgeoh, paph, ptu, pqu,
            plu, zuu, zvu, zmfub, zrfl, ztd, zqd, zud, zvd, pmfd, zmfds, zmfdq, zdmfdp, idtop, loddraf);
    cuddrafn(c, klev, loddraf, ztenh, zqenh, puen, pven, pgeo, zgeoh, paph, zrfl, ztd, zqd, zud, zvd, pmfu, pmfd,
             zmfds, zmfdq, zdmfdp, pmfdde_rate);
  }
  if (ldcum && ktype == 1) {
    const int ikb = kcbot, ikt = kctop;
    zheat = T(0.0);
    zcape = T(0.0);
    zcape1 = T(0.0);
    zcape2 = T(0.0);
    zmfub1 = zmfub;
    ztauc = (zgeoh[ikt] - zgeoh[ikb]) / ((T(2.) + std::fmin(T(15.0), wup)) * c.g);
    if (lndj == 0) {
      upbl = T(2.) + upbl / (paph[klev + 1] - paph[ikb]);
      ztaubl = (zgeoh[ikb] - zgeoh[klev + 1]) / (c.g * upbl);
      ztaubl = std::fmin(T(300.), ztaubl);
    } else {
      ztaubl = ztauc;
    }
  }
  for (int jk = 1; jk <= klev; ++jk) {
    const bool llo1 = ldcum && ktype == 1;
    if (llo1 && jk <= kcbot && jk > kctop) {
      const T zdz = pgeo[jk - 1] - pgeo[jk];
      const T zdp = pap[jk] - pap[jk - 1];
      zheat = zheat + ((pten[jk - 1] - pten[jk] + zdz * c.rcpd) / ztenh[jk] + c.vtmpc1 * (pqen[jk - 1] - pqen[jk])) *
                          (c.g * (pmfu[jk] + pmfd[jk]));
      zcape1 = zcape1 + ((ptu[jk] - ztenh[jk]) / ztenh[jk] + c.vtmpc1 * (pqu[jk] - zqenh[jk]) - plu[jk]) * zdp;
    }
    if (llo1 && jk >= kcbot) {
      if ((paph[klev + 1] - paph[kdpl]) < T(50.e2)) {
        const T zdp = paph[jk + 1] - paph[jk];
        zcape2 = zcape2 + ztaubl * ((T(1.) + c.vtmpc1 * pqen[jk]) * ptte[jk] + c.vtmpc1 * pten[jk] * pqte[jk]) * zdp;
      }
    }
  }
  if (ldcum && ktype == 1) {
    const int ikb = kcbot;
    ztauc = std::fmax(ztmst, ztauc);
    ztauc = std::fmax(T(360.), ztauc);
    ztauc = std::fmin(T(10800.), ztauc);
    const T ztau = ztauc * scale_fac;
    if (C::nonequil) {
      zcape2 = std::fmax(T(0.), zcape2);
      zcape = std::fmax(T(0.), std::fmin(zcape1 - zcape2, T(5000.)));
    } else {
      zcape = std::fmax(T(0.), std::fmin(zcape1, T(5000.)));
    }
    zheat = std::fmax(T(1.e-4), zheat);
    zmfub1 = (zcape * zmfub) / (zheat * ztau);
    zmfub1 = std::fmax(zmfub1, T(0.001));
    const T zmfmax = (paph[ikb] - paph[ikb - 1]) * zcons2;
    zmfub1 = std::fmin(zmfub1, zmfmax);
  }
  if (ldcum && ktype == 2) {
    const int ikb = kcbot;
    T zeps;
    if (pmfd[ikb] < T(0.0) && loddraf) zeps = -pmfd[ikb] / std::fmax(zmfub, cmfcmin);
    else zeps = T(0.);
    const T zqumqe = pqu[ikb] + plu[ikb] - zeps * zqd[ikb] - (T(1.) - zeps) * zqenh[ikb];
    const T zdqmin = std::fmax(T(0.01) * zqenh[ikb], cmfcmin);
    const T zmfmax = (paph[ikb] - paph[ikb - 1]) * zcons2;
    T zdh = c.cpd * (ptu[ikb] - zeps * ztd[ikb] - (T(1.) - zeps) * ztenh[ikb]) + c.alv * zqumqe;
    zdh = c.g * std::fmax(zdh, T(1.e5) * zdqmin);
    if (zdhpbl > T(0.)) zmfub1 = zdhpbl / zdh;
    else zmfub1 = zmfub;
    zmfub1 = zmfub1 / scale_fac2;
    zmfub1 = std::fmin(zmfub1, zmfmax);
  }
  if (ldcum && ktype == 3) zmfub1 = zmfub;
  if (ldcum) {
    const T zfac = zmfub1 / std::fmax(zmfub, cmfcmin);
    for (int jk = 1; jk <= klev; ++jk) {
      pmfd[jk] = pmfd[jk] * zfac;
      zmfds[jk] = zmfds[jk] * zfac;
      zmfdq[jk] = zmfdq[jk] * zfac;
      zdmfdp[jk] = zdmfdp[jk] * zfac;
      pmfdde_rate[jk] = pmfdde_rate[jk] * zfac;
    }
  }
  if (ldcum) zmfs = zmfub1 / std::fmax(cmfcmin, zmfub);
  for (int jk = 2; jk <= klev; ++jk) {
    if (ldcum && jk >= kctop - 1) {
      const int ikb = kcbot;
      if (jk > ikb) {
        const T zdz = ((paph[klev + 1] - paph[jk]) / (paph[klev + 1] - paph[ikb]));
        pmfu[jk] = pmfu[ikb] * zdz;
      }
      const T zmfmax = (paph[jk] - paph[jk - 1]) * zcons2;
      if (pmfu[jk] * zmfs > zmfmax) zmfs = std::fmin(zmfs, zmfmax / pmfu[jk]);
    }
  }
  for (int jk = 2; jk <= klev; ++jk) {
    if (ldcum && jk <= kcbot && jk >= kctop - 1) {
      pmfu[jk] = pmfu[jk] * zmfs;
      zmfus[jk] = zmfus[jk] * zmfs;
      zmfuq[jk] = zmfuq[jk] * zmfs;
      zmful[jk] = zmful[jk] * zmfs;
      zdmfup[jk] = zdmfup[jk] * zmfs;
      plude[jk] = plude[jk] * zmfs;
      pmfude_rate[jk] = pmfude_rate[jk] * zmfs;
    }
  }
  if (ktype == 2 && kcbot == kctop && kcbot >= klev - 1) {
    ldcum = false;
    ktype = 0;
  }
  if (!C::lmfscv || !C::lmfpen) {
    llo2 = false;
    if ((!C::lmfscv && ktype == 2) || (!C::lmfpen && ktype == 1)) {
      llo2 = true;
      ldcum = false;
    }
  }
  if (loddraf && idtop <= kctop) idtop = kctop + 1;
  for (int jk = 2; jk <= klev; ++jk) {
    if (loddraf) {
      if (jk < idtop) {
        pmfd[jk] = T(0.);
        zmfds[jk] = T(0.);
        zmfdq[jk] = T(0.);
        pmfdde_rate[jk] = T(0.);
        zdmfdp[jk] = T(0.);
      } else if (jk == idtop) {
        pmfdde_rate[jk] = T(0.);
      }
    }
  }

  cuflxn(c, klev, ztmst, pten, pqen, pqsen, ztenh, zqenh, paph, pap, zgeoh, lndj, ldcum, kcbot, kctop, idtop, itopm2,
         ktype, loddraf, pmfu, pmfd, zmfus, zmfds, zmfuq, zmfdq, zmful, plude, zdmfup, zdmfdp, zdpmel, zlglac, prain,
         pmfdde_rate, pmflxr, pmflxs);

  zmfs = T(1.);
  zmfuub = T(0.);
  for (int jk = 2; jk <= klev; ++jk) {
    if (loddraf && jk >= idtop - 1) {
      const T zmfmax = pmfu[jk] * T(0.98);
      if (pmfd[jk] + zmfmax + T(1.e-15) < T(0.)) zmfs = std::fmin(zmfs, -zmfmax / pmfd[jk]);
    }
  }
  for (int jk = 2; jk <= klev; ++jk) {
    if (zmfs < T(1.) && jk >= idtop - 1) {
      pmfd[jk] = pmfd[jk] * zmfs;
      zmfds[jk] = zmfds[jk] * zmfs;
      zmfdq[jk] = zmfdq[jk] * zmfs;
      pmfdde_rate[jk] = pmfdde_rate[jk] * zmfs;
      zmfuub = zmfuub - (T(1.) - zmfs) * zdmfdp[jk];
      pmflxr[jk + 1] = pmflxr[jk + 1] + zmfuub;
      zdmfdp[jk] = zdmfdp[jk] * zmfs;
    }
  }
  for (int jk = 2; jk <= klev - 1; ++jk) {
    if (loddraf && jk >= idtop - 1) {
      const T zerate = -pmfd[jk] + pmfd[jk - 1] + pmfdde_rate[jk];
      if (zerate < T(0.)) pmfdde_rate[jk] = pmfdde_rate[jk] - zerate;
    }
    if (ldcum && jk >= kctop - 1) {
      const T zerate = pmfu[jk] - pmfu[jk + 1] + pmfude_rate[jk];
      if (zerate < T(0.)) pmfude_rate[jk] = pmfude_rate[jk] - zerate;
      zdmfup[jk] = pmflxr[jk + 1] + pmflxs[jk + 1] - pmflxr[jk] - pmflxs[jk];
      zdmfdp[jk] = T(0.);
    }
  }
  if (loddraf) {
    const int jk = idtop;
    const int ik = (jk + 1 < klev) ? jk + 1 : klev;
    if (zmfdq[jk] < T(0.3) * zmfdq[ik]) zmfdq[jk] = T(0.3) * zmfdq[ik];
  }
  for (int jk = 2; jk <= klev; ++jk) {
    if (ldcum && jk >= kctop - 1 && jk < kcbot) {
      const T zdz = ztmst * c.g / (paph[jk + 1] - paph[jk]);
      T zmfa = zmfuq[jk + 1] + zmfdq[jk + 1] - zmfuq[jk] - zmfdq[jk] + zmful[jk + 1] - zmful[jk] + zdmfup[jk];
      zmfa = (zmfa - plude[jk]) * zdz;
      if (pqen[jk] + zmfa < T(0.)) plude[jk] = plude[jk] + T(2.) * (pqen[jk] + zmfa) / zdz;
      if (plude[jk] < T(0.)) plude[jk] = T(0.);
    }
    if (!ldcum) pmfude_rate[jk] = T(0.);
    if (std::fabs(pmfd[jk - 1]) < T(1.0e-20)) pmfdde_rate[jk] = T(0.);
  }
  prsfc = pmflxr[klev + 1];
  pssfc = pmflxs[klev + 1];

  cudtdqn(c, klev, itopm2, ldcum, paph, pten, zlglac, plude, pmfu, pmfd, zmfus, zmfds, zmfuq, zmfdq, zmful, zdmfup,
          zdmfdp, zdpmel, ptte, pqte, pcte);

  if (C::lmfdudv) {
    for (int jk = klev - 1; jk >= 2; --jk) {
      const int ik = jk + 1;
      if (ldcum) {
        if (jk == kcbot && ktype < 3) {
          const int ikb = kdpl;
          zuu[jk] = puen[ikb - 1];
          zvu[jk] = pven[ikb - 1];
        } else if (jk == kcbot && ktype == 3) {
          zuu[jk] = puen[jk - 1];
          zvu[jk] = pven[jk - 1];
        }
        if (jk < kcbot && jk >= kctop) {
          if (C::momtrans == 1) {
            T zfac = T(0.);
            if (ktype == 1 || ktype == 3) zfac = T(2.);
            if (ktype == 1 && jk <= kctop + 2) zfac = T(3.);
            const T zerate = pmfu[jk] - pmfu[ik] + (T(1.) + zfac) * pmfude_rate[jk];
            const T zderate = (T(1.) + zfac) * pmfude_rate[jk];
            const T zmfa = T(1.) / std::fmax(cmfcmin, pmfu[jk]);
            zuu[jk] = (zuu[ik] * pmfu[ik] + zerate * puen[jk] - zderate * zuu[ik]) * zmfa;
            zvu[jk] = (zvu[ik] * pmfu[ik] + zerate * pven[jk] - zderate * zvu[ik]) * zmfa;
          } else {
            const T pgf_u = -T(C::pgcoef) * T(0.5) *
                            (pmfu[ik] * (puen[ik] - puen[jk]) + pmfu[jk] * (puen[jk] - puen[jk - 1]));
            const T pgf_v = -T(C::pgcoef) * T(0.5) *
                            (pmfu[ik] * (pven[ik] - pven[jk]) + pmfu[jk] * (pven[jk] - pven[jk - 1]));
            const T zerate = pmfu[jk] - pmfu[ik] + pmfude_rate[jk];
            const T zderate = pmfude_rate[jk];
            const T zmfa = T(1.) / std::fmax(cmfcmin, pmfu[jk]);
            zuu[jk] = (zuu[ik] * pmfu[ik] + zerate * puen[jk] - zderate * zuu[ik] + pgf_u) * zmfa;
            zvu[jk] = (zvu[ik] * pmfu[ik] + zerate * pven[jk] - zderate * zvu[ik] + pgf_v) * zmfa;
          }
        }
      }
    }
    if (C::lmfdd) {
      for (int jk = 3; jk <= klev; ++jk) {
        const int ik = jk - 1;
        if (ldcum) {
          if (jk == idtop) {
            zud[jk] = T(0.5) * (zuu[jk] + puen[ik]);
            zvd[jk] = T(0.5) * (zvu[jk] + pven[ik]);
          } else if (jk > idtop) {
            const T zerate = -pmfd[jk] + pmfd[ik] + pmfdde_rate[jk];
            const T zmfa = T(1.) / std::fmin(-cmfcmin, pmfd[jk]);
            zud[jk] = (zud[ik] * pmfd[ik] - zerate * puen[ik] + pmfdde_rate[jk] * zud[ik]) * zmfa;
            zvd[jk] = (zvd[ik] * pmfd[ik] - zerate * pven[ik] + pmfdde_rate[jk] * zvd[ik]) * zmfa;
          }
        }
      }
    }
    zmfs = T(1.);
    for (int jk = 2; jk <= klev; ++jk) {
      if (ldcum && jk >= kctop - 1) {
        const T zmfmax = (paph[jk] - paph[jk - 1]) * zcons;
        if (pmfu[jk] > zmfmax && jk >= kctop) zmfs = std::fmin(zmfs, zmfmax / pmfu[jk]);
      }
    }
    for (int jk = 1; jk <= klev; ++jk) {
      zmfuus[jk] = pmfu[jk];
      zmfdus[jk] = pmfd[jk];
      if (ldcum && jk >= kctop - 1) {
        zmfuus[jk] = pmfu[jk] * zmfs;
        zmfdus[jk] = pmfd[jk] * zmfs;
      }
    }
    for (int jk = 1; jk <= klev; ++jk) {
      ztenu[jk] = pvom[jk];
      ztenv[jk] = pvol[jk];
    }
    cududvn(c, klev, itopm2, ktype, kcbot, ldcum, paph, puen, pven, zmfuus, zmfdus, zuu, zud, zvu, zvd, pvom, pvol);
    zsum12 = T(0.);
    zsum22 = T(0.);
    for (int jk = 1; jk <= klev; ++jk) {
      zuv2[jk] = T(0.);
      if (ldcum && jk >= kctop - 1) {
        const T zdz = (paph[jk + 1] - paph[jk]);
        const T zduten = pvom[jk] - ztenu[jk];
        const T zdvten = pvol[jk] - ztenv[jk];
        zuv2[jk] = std::sqrt(zduten * zduten + zdvten * zdvten);
        zsum22 = zsum22 + zuv2[jk] * zdz;
        zsum12 = zsum12 - (puen[jk] * zduten + pven[jk] * zdvten) * zdz;
      }
    }
    for (int jk = 1; jk <= klev; ++jk) {
      if (ldcum && jk >= kctop - 1) {
        const T ztdis = c.rcpd * zsum12 * zuv2[jk] / std::fmax(T(1.e-15), zsum22);
        ptte[jk] = ptte[jk] + ztdis;
      }
    }
  }
  if (!C::lmfscv || !C::lmfpen) {
    for (int jk = 2; jk <= klev; ++jk) {
      if (llo2 && jk >= kctop - 1) {
        ptu[jk] = pten[jk];
        pqu[jk] = pqen[jk];
        plu[jk] = T(0.);
        pmfude_rate[jk] = T(0.);
        pmfdde_rate[jk] = T(0.);
      }
    }
    if (llo2) {
      kctop = klev - 1;
      kcbot = klev - 1;
    }
  }
}

// cu_ntiedtke_run for one column (lq = 1). Arrays 1..km (top to bottom), half levels 1..km+1.
// On exit pu, pv, pt, pqv, pqc, pqi are the updated state and zprecc the convective precipitation (mm).
template <typename T>
CPPWRF_HD inline void cu_ntiedtke_run(const NtConst<T>& c, int km, T* pu, T* pv, T* pt, T* pqv, T* pqc, T* pqi,
                                      const T* pqvf, const T* ptf, const T* poz, const T* pzz, const T* pomg,
                                      const T* pap, const T* paph, T evap, T hfx, T& zprecc, int lndj, T dt, T dx) {
  const T ztmst = dt;
  const T dxref = T(15000.);
  T scale_fac, scale_fac2;
  if (dx < dxref) {
    const T a = T(1.06133) + std::log(dxref / dx);
    scale_fac = (a * a) * a;
    scale_fac2 = std::pow(scale_fac, T(0.5));
  } else {
    scale_fac = T(1.) + T(1.33e-5) * dx;
    scale_fac2 = T(1.);
  }
  T zrain = T(0.0);
  bool locum = false;
  T prsfc = T(0.0), pssfc = T(0.0);
  const T pqhfl = evap, phhfl = hfx;
  T pgeoh[KM], pcte[KM], pvom[KM], pvol[KM], ztp1[KM], zqp1[KM], pum1[KM], pvm1[KM], pverv[KM], pgeo[KM], zqsat[KM];
  T pqte[KM], zqq[KM], ptte[KM], ztt[KM], ztu[KM], zqu[KM], zlu[KM], zlude[KM], zmfu[KM], zmfd[KM];
  pgeoh[km + 1] = c.g * pzz[km + 1];
  for (int k = 1; k <= km; ++k) {
    pcte[k] = T(0.0);
    pvom[k] = T(0.0);
    pvol[k] = T(0.0);
    ztp1[k] = pt[k];
    zqp1[k] = pqv[k] / (T(1.0) + pqv[k]);
    pum1[k] = pu[k];
    pvm1[k] = pv[k];
    pverv[k] = pomg[k];
    pgeo[k] = c.g * poz[k];
    pgeoh[k] = c.g * pzz[k];
    const T tt = ztp1[k];
    const T zew = foeewm(c, tt);
    T zqs = zew / pap[k];
    zqs = std::fmin(T(0.5), zqs);
    const T zcor = T(1.) / (T(1.) - c.vtmpc1 * zqs);
    zqsat[k] = zqs * zcor;
    pqte[k] = pqvf[k];
    zqq[k] = pqte[k];
    ptte[k] = ptf[k];
    ztt[k] = ptte[k];
  }
  int ktype = 0, icbot = 0, ictop = 0;
  cumastrn(c, km, ztp1, zqp1, pum1, pvm1, pverv, zqsat, pqhfl, ztmst, pap, paph, pgeo, ptte, pqte, pvom, pvol, prsfc,
           pssfc, locum, ktype, icbot, ictop, ztu, zqu, zlu, zlude, zmfu, zmfd, zrain, pcte, phhfl, lndj, pgeoh, dx,
           scale_fac, scale_fac2);
  for (int k = 1; k <= km; ++k) {
    if (pcte[k] > T(0.)) {
      const T fliq = foealfa(ztp1[k]);
      const T fice = T(1.0) - fliq;
      pqc[k] = pqc[k] + fliq * pcte[k] * ztmst;
      pqi[k] = pqi[k] + fice * pcte[k] * ztmst;
    }
  }
  for (int k = 1; k <= km; ++k) {
    pt[k] = ztp1[k] + (ptte[k] - ztt[k]) * ztmst;
    zqp1[k] = zqp1[k] + (pqte[k] - zqq[k]) * ztmst;
    pqv[k] = zqp1[k] / (T(1.0) - zqp1[k]);
  }
  zprecc = std::fmax(T(0.0), (prsfc + pssfc) * ztmst);
  if (NtConst<T>::lmfdudv) {
    for (int k = 1; k <= km; ++k) {
      pu[k] = pu[k] + pvom[k] * ztmst;
      pv[k] = pv[k] + pvol[k] * ztmst;
    }
  }
}

}  // namespace cppwrf::phys::cu_ntiedtke
