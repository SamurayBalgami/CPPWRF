#pragma once
// Tiedtke cumulus scheme (cu_physics = 6): WRF 4.8.0 phys/module_cu_tiedtke.F (tiecnv, cumastr_new and the
// level 3/4 routines) for one column (klon = 1). Arrays keep the Fortran 1-based level index (k = 1 is the
// model TOP, k = klev the lowest level; klev + 1 is the surface for half levels), sized kMaxLevels + 2.
// Every statement follows the Fortran in order; x**2 is x*x (gfortran powi).
//
// Row coupling: WRF calls the scheme for a row of columns (klon = ite - its + 1). Columns are independent
// except for one test, cumastr_new's "if (icum == 0) return" after the first cloud ascent: it skips the
// rest of the scheme only when no column of the row convects. A column whose first ascent fails can still
// get midlevel convection in the second ascent when another column of its row convects. The driver
// therefore runs the column twice: tiedtke_first_ascent() gives the column's ldcum after the first ascent,
// the row's OR of it is row_icum for tiecnv().
//
// libmvec: gfortran vectorises tiecnv's initialisation loop over the row (zqsat = tlucua(t)/p, tlucua
// inlined), so its exp is glibc's 2-lane vector exp for column `lane` of `nlane` (cppwrf::vm::exp).
#include <cmath>

#include "core/backend.hpp"

namespace cppwrf::phys::cu_tiedtke {

inline constexpr int kMaxLevels = 100;
inline constexpr int KM = kMaxLevels + 2;

// module_cu_tiedtke parameters (module_model_constants rd, rv, cpd, alv, als, alf, g)
template <typename T>
struct TdConst {
  T rd, rv, cpd, alv, als, alf, g;
  T rcpd, zrg, c2es, c5les, c5ies, vtmpc1, cevapcu1, cevapcu2, cprcon;
  static constexpr double rhoh2o = 1.0e03, tmelt = 273.16, t000 = 273.15, hgfr = 233.15;
  static constexpr double c1es = 610.78, c3les = 17.269, c4les = 35.86, c3ies = 21.875, c4ies = 7.66;
  static constexpr double entrpen = 1.0e-4, entrscv = 1.2e-3, entrmid = 1.0e-4, entrdd = 2.0e-4;
  static constexpr double cmfctop = 0.30, cmfcmax = 1.0, cmfcmin = 1.e-10, cmfdeps = 0.30, zdnoprc = 1.5e4;
  static constexpr double rhc = 0.80, rhm = 1.0, zbuo0 = 0.50, crirh = 0.70, fdbk = 1.0, ztau = 2400.0;
  static constexpr bool lmfpen = true, lmfmid = true, lmfscv = true, lmfdd = true, lmfdudv = true;
};

template <typename T>
inline TdConst<T> make_tdconst(T r_d, T r_v, T cp, T xlv, T xls, T xlf, T g) {
  using C = TdConst<T>;
  TdConst<T> c{};
  c.rd = r_d; c.rv = r_v; c.cpd = cp; c.alv = xlv; c.als = xls; c.alf = xlf; c.g = g;
  c.rcpd = T(1.0) / cp;
  c.zrg = T(1.0) / g;
  c.c2es = T(C::c1es) * r_d / r_v;
  c.c5les = T(C::c3les) * (T(C::tmelt) - T(C::c4les));
  c.c5ies = T(C::c3ies) * (T(C::tmelt) - T(C::c4ies));
  c.vtmpc1 = r_v / r_d - T(1.0);
  c.cevapcu1 = T(1.93e-6) * T(261.0) * T(0.5) / g;
  c.cevapcu2 = T(1.e3) / (T(38.3) * T(0.293));
  c.cprcon = T(1.1e-3) / g;
  return c;
}

// tlucua, tlucub, tlucuc (lane/nlane: libmvec exp of a vectorised caller; nlane = 0: scalar exp)
template <typename T>
CPPWRF_HD inline T tlucua(const TdConst<T>& c, T tt, int lane = 0, int nlane = 0) {
  using C = TdConst<T>;
  T zcvm3, zcvm4;
  if (tt - T(C::tmelt) > T(0.)) {
    zcvm3 = T(C::c3les);
    zcvm4 = T(C::c4les);
  } else {
    zcvm3 = T(C::c3ies);
    zcvm4 = T(C::c4ies);
  }
  return c.c2es * vm::exp(zcvm3 * (tt - T(C::tmelt)) * (T(1.) / (tt - zcvm4)), lane, nlane);
}
template <typename T>
CPPWRF_HD inline T tlucub(const TdConst<T>& c, T tt) {
  using C = TdConst<T>;
  const T z5alvcp = c.c5les * c.alv / c.cpd;
  const T z5alscp = c.c5ies * c.als / c.cpd;
  T zcvm4, zcvm5;
  if (tt - T(C::tmelt) > T(0.)) {
    zcvm4 = T(C::c4les);
    zcvm5 = z5alvcp;
  } else {
    zcvm4 = T(C::c4ies);
    zcvm5 = z5alscp;
  }
  const T r = T(1.) / (tt - zcvm4);
  return zcvm5 * (r * r);
}
template <typename T>
CPPWRF_HD inline T tlucuc(const TdConst<T>& c, T tt) {
  using C = TdConst<T>;
  const T zalvdcp = c.alv / c.cpd;
  const T zalsdcp = c.als / c.cpd;
  return tt - T(C::tmelt) > T(0.) ? zalvdcp : zalsdcp;
}

// cuadjtq: t and q adjustment at level kk (kcall 0: environment, 1: condensation, 2: evaporation)
template <typename T>
CPPWRF_HD inline void cuadjtq(const TdConst<T>& c, int kk, T pp, T* pt, T* pq, bool ldflag, int kcall) {
  if (kcall == 1 || kcall == 2) {
    if (!ldflag) return;
    const T zqp = T(1.) / pp;
    T tt = pt[kk];
    T zqsat = tlucua(c, tt) * zqp;
    zqsat = c_min(T(0.5), zqsat);
    T zcor = T(1.) / (T(1.) - c.vtmpc1 * zqsat);
    zqsat = zqsat * zcor;
    T zcond = (pq[kk] - zqsat) / (T(1.) + zqsat * zcor * tlucub(c, tt));
    zcond = kcall == 1 ? c_max(zcond, T(0.)) : c_min(zcond, T(0.));
    pt[kk] = pt[kk] + tlucuc(c, tt) * zcond;
    pq[kk] = pq[kk] - zcond;
    if (zcond != T(0.)) {
      tt = pt[kk];
      zqsat = tlucua(c, tt) * zqp;
      zqsat = c_min(T(0.5), zqsat);
      zcor = T(1.) / (T(1.) - c.vtmpc1 * zqsat);
      zqsat = zqsat * zcor;
      const T zcond1 = (pq[kk] - zqsat) / (T(1.) + zqsat * zcor * tlucub(c, tt));
      pt[kk] = pt[kk] + tlucuc(c, tt) * zcond1;
      pq[kk] = pq[kk] - zcond1;
    }
  } else if (kcall == 0) {
    // the second pass runs for every column of the row when any column condensed; with zcond = 0 it
    // repeats the first pass exactly (zcond1 = 0), so the column alone gives the same result
    const T zqp = T(1.) / pp;
    T tt = pt[kk];
    T zqsat = tlucua(c, tt) * zqp;
    zqsat = c_min(T(0.5), zqsat);
    T zcor = T(1.) / (T(1.) - c.vtmpc1 * zqsat);
    zqsat = zqsat * zcor;
    const T zcond = (pq[kk] - zqsat) / (T(1.) + zqsat * zcor * tlucub(c, tt));
    pt[kk] = pt[kk] + tlucuc(c, tt) * zcond;
    pq[kk] = pq[kk] - zcond;
    if (zcond != T(0.)) {
      tt = pt[kk];
      zqsat = tlucua(c, tt) * zqp;
      zqsat = c_min(T(0.5), zqsat);
      zcor = T(1.) / (T(1.) - c.vtmpc1 * zqsat);
      zqsat = zqsat * zcor;
      const T zcond1 = (pq[kk] - zqsat) / (T(1.) + zqsat * zcor * tlucub(c, tt));
      pt[kk] = pt[kk] + tlucuc(c, tt) * zcond1;
      pq[kk] = pq[kk] - zcond1;
    }
  }
}

// cuini: half-level environment values, updraft/downdraft initial values
template <typename T>
CPPWRF_HD inline void cuini(const TdConst<T>& c, int klev, const T* pten, const T* pqen, const T* pqsen,
                            const T* puen, const T* pven, const T* pverv, const T* pgeo, const T* paph, T* pgeoh,
                            T* ptenh, T* pqenh, T* pqsenh, int& klwmin, T* ptu, T* pqu, T* ptd, T* pqd, T* puu,
                            T* pvu, T* pud, T* pvd, T* pmfu, T* pmfd, T* pmfus, T* pmfds, T* pmfuq, T* pmfdq,
                            T* pdmfup, T* pdmfdp, T* pdpmel, T* plu, T* plude, int* klab) {
  const int klevm1 = klev - 1;
  const T zdp = T(0.5);
  for (int jk = 2; jk <= klev; ++jk) {
    pgeoh[jk] = pgeo[jk] + (pgeo[jk - 1] - pgeo[jk]) * zdp;
    ptenh[jk] = (c_max(c.cpd * pten[jk - 1] + pgeo[jk - 1], c.cpd * pten[jk] + pgeo[jk]) - pgeoh[jk]) * c.rcpd;
    pqsenh[jk] = pqsen[jk - 1];
    const T zph = paph[jk];
    cuadjtq(c, jk, zph, ptenh, pqsenh, true, 0);
    pqenh[jk] = c_min(pqen[jk - 1], pqsen[jk - 1]) + (pqsenh[jk] - pqsen[jk - 1]);
    pqenh[jk] = c_max(pqenh[jk], T(0.));
  }
  ptenh[klev] = (c.cpd * pten[klev] + pgeo[klev] - pgeoh[klev]) * c.rcpd;
  pqenh[klev] = pqen[klev];
  ptenh[1] = pten[1];
  pqenh[1] = pqen[1];
  pgeoh[1] = pgeo[1];
  klwmin = klev;
  T zwmax = T(0.);
  for (int jk = klevm1; jk >= 2; --jk) {
    const T zzs = c_max(c.cpd * ptenh[jk] + pgeoh[jk], c.cpd * ptenh[jk + 1] + pgeoh[jk + 1]);
    ptenh[jk] = (zzs - pgeoh[jk]) * c.rcpd;
  }
  for (int jk = klev; jk >= 3; --jk) {
    if (pverv[jk] < zwmax) {
      zwmax = pverv[jk];
      klwmin = jk;
    }
  }
  for (int jk = 1; jk <= klev; ++jk) {
    const int ik = jk == 1 ? 1 : jk - 1;
    ptu[jk] = ptenh[jk];
    ptd[jk] = ptenh[jk];
    pqu[jk] = pqenh[jk];
    pqd[jk] = pqenh[jk];
    plu[jk] = T(0.);
    puu[jk] = puen[ik];
    pud[jk] = puen[ik];
    pvu[jk] = pven[ik];
    pvd[jk] = pven[ik];
    pmfu[jk] = T(0.);
    pmfd[jk] = T(0.);
    pmfus[jk] = T(0.);
    pmfds[jk] = T(0.);
    pmfuq[jk] = T(0.);
    pmfdq[jk] = T(0.);
    pdmfup[jk] = T(0.);
    pdmfdp[jk] = T(0.);
    pdpmel[jk] = T(0.);
    plude[jk] = T(0.);
    klab[jk] = 0;
  }
}

// cubase: cloud base (subcloud ascent). The row's "if (is == 0) cycle" only skips levels where this
// column does not ascend either: a column that finds a cloud base ascends at every level below it.
template <typename T>
CPPWRF_HD inline void cubase(const TdConst<T>& c, int klev, const T* ptenh, const T* pqenh, const T* pgeoh,
                             const T* paph, T* ptu, T* pqu, T* plu, const T* puen, const T* pven, T* puu, T* pvu,
                             bool& ldcum, int& kcbot, int* klab) {
  using C = TdConst<T>;
  const int klevm1 = klev - 1, klevp1 = klev + 1;
  T zqold[KM];
  klab[klev] = 1;
  kcbot = klevm1;
  ldcum = false;
  bool ldbase = false;
  puu[klev] = puen[klev] * (paph[klevp1] - paph[klev]);
  pvu[klev] = pven[klev] * (paph[klevp1] - paph[klev]);
  for (int jk = 1; jk <= klev; ++jk) zqold[jk] = T(0.0);
  for (int jk = klevm1; jk >= 2; --jk) {
    const bool loflag = klab[jk + 1] == 1 || (ldcum && kcbot == jk + 1);
    const T zph = paph[jk];
    if (!loflag) continue;
    if (C::lmfdudv) {
      if (!ldbase) {
        puu[klev] = puu[klev] + puen[jk] * (paph[jk + 1] - paph[jk]);
        pvu[klev] = pvu[klev] + pven[jk] * (paph[jk + 1] - paph[jk]);
      }
    }
    pqu[jk] = pqu[jk + 1];
    ptu[jk] = (c.cpd * ptu[jk + 1] + pgeoh[jk + 1] - pgeoh[jk]) * c.rcpd;
    zqold[jk] = pqu[jk];
    cuadjtq(c, jk, zph, ptu, pqu, loflag, 1);
    if (pqu[jk] == zqold[jk]) {
      const T zbuo = ptu[jk] * (T(1.) + c.vtmpc1 * pqu[jk]) - ptenh[jk] * (T(1.) + c.vtmpc1 * pqenh[jk]) + T(C::zbuo0);
      if (zbuo > T(0.)) klab[jk] = 1;
    } else {
      klab[jk] = 2;
      plu[jk] = plu[jk] + zqold[jk] - pqu[jk];
      const T zbuo = ptu[jk] * (T(1.) + c.vtmpc1 * pqu[jk]) - ptenh[jk] * (T(1.) + c.vtmpc1 * pqenh[jk]) + T(C::zbuo0);
      const bool llo1 = zbuo > T(0.) && klab[jk + 1] == 1;
      if (llo1) {
        kcbot = jk;
        ldcum = true;
        ldbase = true;
      }
    }
  }
  if (C::lmfdudv) {
    if (ldcum) {
      const int ikb = kcbot;
      const T zz = T(1.) / (paph[klevp1] - paph[ikb]);
      puu[klev] = puu[klev] * zz;
      pvu[klev] = pvu[klev] * zz;
    } else {
      puu[klev] = puen[klevm1];
      pvu[klev] = pven[klevm1];
    }
  }
}

// cubasmc: cloud base values for midlevel convection at level kk
template <typename T>
CPPWRF_HD inline void cubasmc(const TdConst<T>& c, int kk, const T* pten, const T* pqen, const T* pqsen,
                              const T* puen, const T* pven, const T* pverv, const T* pgeo, const T* pgeoh,
                              bool ldcum, int& ktype, int* klab, T* pmfu, T& pmfub, T& pentr, int& kcbot, T* ptu,
                              T* pqu, T* plu, T* puu, T* pvu, T* pmfus, T* pmfuq, T* pmful, T* pdmfup, T& pmfuu,
                              T& pmfuv) {
  using C = TdConst<T>;
  if (!ldcum && klab[kk + 1] == 0 && pqen[kk] > T(0.80) * pqsen[kk]) {
    ptu[kk + 1] = (c.cpd * pten[kk] + pgeo[kk] - pgeoh[kk + 1]) * c.rcpd;
    pqu[kk + 1] = pqen[kk];
    plu[kk + 1] = T(0.);
    T zzzmb = c_max(T(C::cmfcmin), -pverv[kk] / c.g);
    zzzmb = c_min(zzzmb, T(C::cmfcmax));
    pmfub = zzzmb;
    pmfu[kk + 1] = pmfub;
    pmfus[kk + 1] = pmfub * (c.cpd * ptu[kk + 1] + pgeoh[kk + 1]);
    pmfuq[kk + 1] = pmfub * pqu[kk + 1];
    pmful[kk + 1] = T(0.);
    pdmfup[kk + 1] = T(0.);
    kcbot = kk;
    klab[kk + 1] = 1;
    ktype = 3;
    pentr = T(C::entrmid);
    if (C::lmfdudv) {
      puu[kk + 1] = puen[kk];
      pvu[kk + 1] = pven[kk];
      pmfuu = pmfub * puu[kk + 1];
      pmfuv = pmfub * pvu[kk + 1];
    }
  }
}

// cuentr_new: entrainment / detrainment rates at level kk
template <typename T>
CPPWRF_HD inline void cuentr_new(const TdConst<T>& c, int kk, const T* ptenh, const T* paph, const T* pap,
                                 const T* pgeoh, int klwmin, bool ldcum, int ktype, int kcbot, int kctop0,
                                 T& zpbase, const T* pmfu, T pentr, T& zdmfen, T& zdmfde, T* zodetr, int khmin) {
  using C = TdConst<T>;
  zpbase = paph[kcbot];
  const T zrrho = (c.rd * ptenh[kk + 1]) / paph[kk + 1];
  T zdprho = (paph[kk + 1] - paph[kk]) * c.zrg;
  const T zpmid = T(0.5) * (zpbase + paph[kctop0]);
  const T zentr = pentr * pmfu[kk + 1] * zdprho * zrrho;
  const bool llo1 = kk < kcbot && ldcum;
  zdmfde = llo1 ? zentr : T(0.0);
  bool llo2 = llo1 && ktype == 2 && ((zpbase - paph[kk]) < T(C::zdnoprc) || paph[kk] > zpmid);
  zdmfen = llo2 ? zentr : T(0.0);
  const int iklwmin = c_max(klwmin, kctop0 + 2);
  llo2 = llo1 && ktype == 3 && (kk >= iklwmin || pap[kk] > zpmid);
  if (llo2) zdmfen = zentr;
  llo2 = llo1 && ktype == 1;
  if (llo2) zdmfen = zentr;  // turbulent entrainment
  // organized detrainment, detrainment starts at khmin
  zodetr[kk] = T(0.);
  if (llo2 && kk <= khmin && kk >= kctop0) {
    const int ikt = kctop0;
    const int ikh = khmin;
    if (ikh > ikt) {
      const T zzmzk = -(pgeoh[ikh] - pgeoh[kk]) * c.zrg;
      const T ztmzk = -(pgeoh[ikh] - pgeoh[ikt]) * c.zrg;
      const T arg = T(3.1415) * (zzmzk / ztmzk) * T(0.5);
      const T zorgde = std::tan(arg) * T(3.1415) * T(0.5) / ztmzk;
      zdprho = (paph[kk + 1] - paph[kk]) * (c.zrg * zrrho);
      zodetr[kk] = c_min(zorgde, T(1.e-3)) * pmfu[kk + 1] * zdprho;
    }
  }
}

// cuasc_new: cloud ascent (entraining / detraining plume). The row's "if (is == 0) cycle" skips levels
// where this column does not ascend; what it skips for such a column (cuentr_new rates, the momentum
// mixing of zdmfen / zdmfde) is recomputed before it is read again.
template <typename T>
CPPWRF_HD inline void cuasc_new(const TdConst<T>& c, int klev, const T* ptenh, const T* pqenh, const T* puen,
                                const T* pven, const T* pten, const T* pqen, const T* pqsen, const T* pgeo,
                                const T* pgeoh, const T* pap, const T* paph, const T* pverv, int klwmin, bool& ldcum,
                                T phcbase, int& ktype, int* klab, T* ptu, T* pqu, T* plu, T* puu, T* pvu, T* pmfu,
                                T& pmfub, T& pentr, T* pmfus, T* pmfuq, T* pmful, T* plude, T* pdmfup, int& kcbot,
                                int& kctop, int& kctop0, int& kcum, T ztmst, int khmin, const T* phhatt,
                                const T* pqsenh) {
  using C = TdConst<T>;
  const int klevm1 = klev - 1;
  const T zcons2 = T(1.) / (c.g * ztmst);
  T zodetr[KM], zoentr[KM];
  T zmfuu = T(0.), zmfuv = T(0.), zbuoy = T(0.), zpbase = T(0.), zqold = T(0.), zdmfen = T(0.), zdmfde = T(0.);
  if (!ldcum) ktype = 0;
  for (int jk = 1; jk <= klev; ++jk) {
    plu[jk] = T(0.);
    pmfu[jk] = T(0.);
    pmfus[jk] = T(0.);
    pmfuq[jk] = T(0.);
    pmful[jk] = T(0.);
    plude[jk] = T(0.);
    pdmfup[jk] = T(0.);
    zoentr[jk] = T(0.);
    zodetr[jk] = T(0.);
    if (!ldcum || ktype == 3) klab[jk] = 0;
    if (!ldcum && paph[jk] < T(4.e4)) kctop0 = jk;
  }
  kctop = klevm1;
  if (!ldcum) {
    kcbot = klevm1;
    pmfub = T(0.);
    pqu[klev] = T(0.);
  }
  pmfu[klev] = pmfub;
  pmfus[klev] = pmfub * (c.cpd * ptu[klev] + pgeoh[klev]);
  pmfuq[klev] = pmfub * pqu[klev];
  if (C::lmfdudv) {
    zmfuu = pmfub * puu[klev];
    zmfuv = pmfub * pvu[klev];
  }
  // 3.1 organized entrainment at cloud base
  ldcum = false;
  if (ktype == 1) {
    const int ikb = kcbot;
    zbuoy = c.g * ((ptu[ikb] - ptenh[ikb]) / ptenh[ikb] + T(0.608) * (pqu[ikb] - pqenh[ikb]));
    if (zbuoy > T(0.)) {
      const T zdz = (pgeo[ikb - 1] - pgeo[ikb]) * c.zrg;
      const T zdrodz = -std::log(pten[ikb - 1] / pten[ikb]) / zdz - c.g / (c.rd * ptenh[ikb]);
      zoentr[ikb - 1] = zbuoy * T(0.5) / (T(1.) + zbuoy * zdz) + zdrodz;
      zoentr[ikb - 1] = c_min(zoentr[ikb - 1], T(1.e-3));
      zoentr[ikb - 1] = c_max(zoentr[ikb - 1], T(0.));
    }
  }
  // 4. ascent
  for (int jk = klevm1; jk >= 2; --jk) {
    const int ik = jk;
    if (C::lmfmid && ik < klevm1 && ik > klev - 13)
      cubasmc(c, ik, pten, pqen, pqsen, puen, pven, pverv, pgeo, pgeoh, ldcum, ktype, klab, pmfu, pmfub, pentr, kcbot,
              ptu, pqu, plu, puu, pvu, pmfus, pmfuq, pmful, pdmfup, zmfuu, zmfuv);
    zqold = T(0.0);
    const int is = klab[jk + 1];
    if (klab[jk + 1] == 0) klab[jk] = 0;
    const bool loflag = klab[jk + 1] > 0;
    const T zph = paph[jk];
    if (ktype == 3 && jk == kcbot) {
      const T zmfmax = (paph[jk] - paph[jk - 1]) * zcons2;
      if (pmfub > zmfmax) {
        const T zfac = zmfmax / pmfub;
        pmfu[jk + 1] = pmfu[jk + 1] * zfac;
        pmfus[jk + 1] = pmfus[jk + 1] * zfac;
        pmfuq[jk + 1] = pmfuq[jk + 1] * zfac;
        zmfuu = zmfuu * zfac;
        zmfuv = zmfuv * zfac;
        pmfub = zmfmax;
      }
    }
    if (is == 0) continue;
    cuentr_new(c, ik, ptenh, paph, pap, pgeoh, klwmin, ldcum, ktype, kcbot, kctop0, zpbase, pmfu, pentr, zdmfen,
               zdmfde, zodetr, khmin);
    if (loflag) {
      if (jk < kcbot) {
        const T zmftest = pmfu[jk + 1] + zdmfen - zdmfde;
        const T zmfmax = c_min(zmftest, (paph[jk] - paph[jk - 1]) * zcons2);
        zdmfen = c_max(zdmfen - c_max(zmftest - zmfmax, T(0.)), T(0.));
      }
      zdmfde = c_min(zdmfde, T(0.75) * pmfu[jk + 1]);
      pmfu[jk] = pmfu[jk + 1] + zdmfen - zdmfde;
      if (jk < kcbot) {
        const T zdprho = (pgeoh[jk] - pgeoh[jk + 1]) * c.zrg;
        zoentr[jk] = zoentr[jk] * zdprho * pmfu[jk + 1];
        const T zmftest = pmfu[jk] + zoentr[jk] - zodetr[jk];
        const T zmfmax = c_min(zmftest, (paph[jk] - paph[jk - 1]) * zcons2);
        zoentr[jk] = c_max(zoentr[jk] - c_max(zmftest - zmfmax, T(0.)), T(0.));
      }
      // limit organized detrainment to not allowing for too deep clouds
      if (ktype == 1 && jk < kcbot && jk <= khmin) {
        const T zmse = c.cpd * ptu[jk + 1] + c.alv * pqu[jk + 1] + pgeoh[jk + 1];
        const int ikt = kctop0;
        T znevn = (pgeoh[ikt] - pgeoh[jk + 1]) * (zmse - phhatt[jk + 1]) * c.zrg;
        if (znevn <= T(0.)) znevn = T(1.);
        const T zdprho = (pgeoh[jk] - pgeoh[jk + 1]) * c.zrg;
        T zodmax = ((phcbase - zmse) / znevn) * zdprho * pmfu[jk + 1];
        zodmax = c_max(zodmax, T(0.));
        zodetr[jk] = c_min(zodetr[jk], zodmax);
      }
      zodetr[jk] = c_min(zodetr[jk], T(0.75) * pmfu[jk]);
      pmfu[jk] = pmfu[jk] + zoentr[jk] - zodetr[jk];
      T zqeen = pqenh[jk + 1] * zdmfen;
      zqeen = zqeen + pqenh[jk + 1] * zoentr[jk];
      T zseen = (c.cpd * ptenh[jk + 1] + pgeoh[jk + 1]) * zdmfen;
      zseen = zseen + (c.cpd * ptenh[jk + 1] + pgeoh[jk + 1]) * zoentr[jk];
      T zscde = (c.cpd * ptu[jk + 1] + pgeoh[jk + 1]) * zdmfde;
      // find moist static energy that give nonbuoyant air
      const T zga = c.alv * pqsenh[jk + 1] / (c.rv * (ptenh[jk + 1] * ptenh[jk + 1]));
      const T zdt = (plu[jk + 1] - T(0.608) * (pqsenh[jk + 1] - pqenh[jk + 1])) / (T(1.) / ptenh[jk + 1] + T(0.608) * zga);
      const T zscod = c.cpd * ptenh[jk + 1] + pgeoh[jk + 1] + c.cpd * zdt;
      zscde = zscde + zodetr[jk] * zscod;
      T zqude = pqu[jk + 1] * zdmfde;
      const T zqcod = pqsenh[jk + 1] + zga * zdt;
      zqude = zqude + zodetr[jk] * zqcod;
      plude[jk] = plu[jk + 1] * zdmfde;
      plude[jk] = plude[jk] + plu[jk + 1] * zodetr[jk];
      const T zmfusk = pmfus[jk + 1] + zseen - zscde;
      const T zmfuqk = pmfuq[jk + 1] + zqeen - zqude;
      const T zmfulk = pmful[jk + 1] - plude[jk];
      plu[jk] = zmfulk * (T(1.) / c_max(T(C::cmfcmin), pmfu[jk]));
      pqu[jk] = zmfuqk * (T(1.) / c_max(T(C::cmfcmin), pmfu[jk]));
      ptu[jk] = (zmfusk * (T(1.) / c_max(T(C::cmfcmin), pmfu[jk])) - pgeoh[jk]) * c.rcpd;
      ptu[jk] = c_max(T(100.), ptu[jk]);
      ptu[jk] = c_min(T(400.), ptu[jk]);
      zqold = pqu[jk];
    }
    // corrections for moist ascent
    cuadjtq(c, ik, zph, ptu, pqu, loflag, 1);
    if (loflag && pqu[jk] != zqold) {
      klab[jk] = 2;
      plu[jk] = plu[jk] + zqold - pqu[jk];
      T zbuo = ptu[jk] * (T(1.) + c.vtmpc1 * pqu[jk] - plu[jk]) - ptenh[jk] * (T(1.) + c.vtmpc1 * pqenh[jk]);
      if (klab[jk + 1] == 1) zbuo = zbuo + T(C::zbuo0);
      if (zbuo > T(0.) && pmfu[jk] > T(0.01) * pmfub && jk >= kctop0) {
        kctop = jk;
        ldcum = true;
        const T zprcon = zpbase - paph[jk] >= T(C::zdnoprc) ? c.cprcon : T(0.);
        const T zlnew = plu[jk] / (T(1.) + zprcon * (pgeoh[jk] - pgeoh[jk + 1]));
        pdmfup[jk] = c_max(T(0.), (plu[jk] - zlnew) * pmfu[jk]);
        plu[jk] = zlnew;
      } else {
        klab[jk] = 0;
        pmfu[jk] = T(0.);
      }
    }
    if (loflag) {
      pmful[jk] = plu[jk] * pmfu[jk];
      pmfus[jk] = (c.cpd * ptu[jk] + pgeoh[jk]) * pmfu[jk];
      pmfuq[jk] = pqu[jk] * pmfu[jk];
    }
    if (C::lmfdudv) {
      zdmfen = zdmfen + zoentr[jk];
      zdmfde = zdmfde + zodetr[jk];
      if (loflag) {
        T zz;
        if (ktype == 1 || ktype == 3) {
          zz = zdmfen <= T(1.e-20) ? T(3.) : T(2.);
        } else {
          zz = zdmfen <= T(1.0e-20) ? T(1.) : T(0.);
        }
        const T zdmfeu = zdmfen + zz * zdmfde;
        T zdmfdu = zdmfde + zz * zdmfde;
        zdmfdu = c_min(zdmfdu, T(0.75) * pmfu[jk + 1]);
        zmfuu = zmfuu + zdmfeu * puen[jk] - zdmfdu * puu[jk + 1];
        zmfuv = zmfuv + zdmfeu * pven[jk] - zdmfdu * pvu[jk + 1];
        if (pmfu[jk] > T(0.)) {
          puu[jk] = zmfuu * (T(1.) / pmfu[jk]);
          pvu[jk] = zmfuv * (T(1.) / pmfu[jk]);
        }
      }
    }
    // organized entrainment for use at the next level
    if (loflag && ktype == 1) {
      T zbuoyz = c.g * ((ptu[jk] - ptenh[jk]) / ptenh[jk] + T(0.608) * (pqu[jk] - pqenh[jk]) - plu[jk]);
      zbuoyz = c_max(zbuoyz, T(0.0));
      const T zdz = (pgeo[jk - 1] - pgeo[jk]) * c.zrg;
      const T zdrodz = -std::log(pten[jk - 1] / pten[jk]) / zdz - c.g / (c.rd * ptenh[jk]);
      zbuoy = zbuoy + zbuoyz * zdz;
      zoentr[jk - 1] = zbuoyz * T(0.5) / (T(1.) + zbuoy) + zdrodz;
      zoentr[jk - 1] = c_min(zoentr[jk - 1], T(1.e-3));
      zoentr[jk - 1] = c_max(zoentr[jk - 1], T(0.));
    }
  }
  // 5. convective fluxes above the non-buoyancy level
  if (kctop == klevm1) ldcum = false;
  kcbot = c_max(kcbot, kctop);
  kcum = ldcum ? 1 : 0;
  if (!ldcum) return;
  {
    const int jk = kctop - 1;
    const T zzdmf = T(C::cmfctop);
    zdmfde = (T(1.) - zzdmf) * pmfu[jk + 1];
    plude[jk] = zdmfde * plu[jk + 1];
    pmfu[jk] = pmfu[jk + 1] - zdmfde;
    pmfus[jk] = (c.cpd * ptu[jk] + pgeoh[jk]) * pmfu[jk];
    pmfuq[jk] = pqu[jk] * pmfu[jk];
    pmful[jk] = plu[jk] * pmfu[jk];
    plude[jk - 1] = pmful[jk];
    pdmfup[jk] = T(0.);
  }
  if (C::lmfdudv) {
    const int jk = kctop - 1;
    puu[jk] = puu[jk + 1];
    pvu[jk] = pvu[jk + 1];
  }
}

// cudlfs: level of free sinking for downdrafts
template <typename T>
CPPWRF_HD inline void cudlfs(const TdConst<T>& c, int klev, const T* ptenh, const T* pqenh, const T* puen,
                             const T* pven, const T* pgeoh, const T* paph, const T* ptu, const T* pqu, const T* puu,
                             const T* pvu, bool ldcum, int kcbot, int kctop, T pmfub, T& prfl, T* ptd, T* pqd,
                             T* pud, T* pvd, T* pmfd, T* pmfds, T* pmfdq, T* pdmfdp, int& kdtop, bool& lddraf) {
  using C = TdConst<T>;
  T ztenwb[KM], zqenwb[KM];
  lddraf = false;
  kdtop = klev + 1;
  if (!C::lmfdd) return;
  const int ke = klev - 3;
  for (int jk = 3; jk <= ke; ++jk) {
    ztenwb[jk] = ptenh[jk];
    zqenwb[jk] = pqenh[jk];
    const T zph = paph[jk];
    const bool llo2 = ldcum && prfl > T(0.) && !lddraf && (jk < kcbot && jk > kctop);
    if (!llo2) continue;  // the row's is == 0 cycle; pmfd(jk) < 0 below needs llo2 at this level
    cuadjtq(c, jk, zph, ztenwb, zqenwb, llo2, 2);
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
    if (C::lmfdudv) {
      if (pmfd[jk] < T(0.)) {
        pud[jk] = T(0.5) * (puu[jk] + puen[jk - 1]);
        pvd[jk] = T(0.5) * (pvu[jk] + pven[jk - 1]);
      }
    }
  }
}

// cuddraf: downdraft descent
template <typename T>
CPPWRF_HD inline void cuddraf(const TdConst<T>& c, int klev, const T* ptenh, const T* pqenh, const T* puen,
                              const T* pven, const T* pgeoh, const T* paph, T& prfl, bool lddraf, T* ptd, T* pqd,
                              T* pud, T* pvd, T* pmfd, T* pmfds, T* pmfdq, T* pdmfdp) {
  using C = TdConst<T>;
  const int klevp1 = klev + 1;
  for (int jk = 3; jk <= klev; ++jk) {
    const T zph = paph[jk];
    const bool llo2 = lddraf && pmfd[jk - 1] < T(0.);
    if (!llo2) continue;
    const T zentr = T(C::entrdd) * pmfd[jk - 1] * c.rd * ptenh[jk - 1] / (c.g * paph[jk - 1]) * (paph[jk] - paph[jk - 1]);
    T zdmfen = zentr;
    T zdmfde = zentr;
    const int itopde = klev - 2;
    if (jk > itopde) {
      zdmfen = T(0.);
      zdmfde = pmfd[itopde] * (paph[jk] - paph[jk - 1]) / (paph[klevp1] - paph[itopde]);
    }
    pmfd[jk] = pmfd[jk - 1] + zdmfen - zdmfde;
    const T zseen = (c.cpd * ptenh[jk - 1] + pgeoh[jk - 1]) * zdmfen;
    const T zqeen = pqenh[jk - 1] * zdmfen;
    const T zsdde = (c.cpd * ptd[jk - 1] + pgeoh[jk - 1]) * zdmfde;
    const T zqdde = pqd[jk - 1] * zdmfde;
    const T zmfdsk = pmfds[jk - 1] + zseen - zsdde;
    const T zmfdqk = pmfdq[jk - 1] + zqeen - zqdde;
    pqd[jk] = zmfdqk * (T(1.) / c_min(-T(C::cmfcmin), pmfd[jk]));
    ptd[jk] = (zmfdsk * (T(1.) / c_min(-T(C::cmfcmin), pmfd[jk])) - pgeoh[jk]) * c.rcpd;
    ptd[jk] = c_min(T(400.), ptd[jk]);
    ptd[jk] = c_max(T(100.), ptd[jk]);
    T zcond = pqd[jk];
    cuadjtq(c, jk, zph, ptd, pqd, llo2, 2);
    zcond = zcond - pqd[jk];
    const T zbuo = ptd[jk] * (T(1.) + c.vtmpc1 * pqd[jk]) - ptenh[jk] * (T(1.) + c.vtmpc1 * pqenh[jk]);
    if (zbuo >= T(0.) || prfl <= (pmfd[jk] * zcond)) pmfd[jk] = T(0.);
    pmfds[jk] = (c.cpd * ptd[jk] + pgeoh[jk]) * pmfd[jk];
    pmfdq[jk] = pqd[jk] * pmfd[jk];
    const T zdmfdp = -pmfd[jk] * zcond;
    pdmfdp[jk - 1] = zdmfdp;
    prfl = prfl + zdmfdp;
    if (C::lmfdudv) {
      if (pmfd[jk] < T(0.)) {
        const T zmfduk = pmfd[jk - 1] * pud[jk - 1] + zdmfen * puen[jk - 1] - zdmfde * pud[jk - 1];
        const T zmfdvk = pmfd[jk - 1] * pvd[jk - 1] + zdmfen * pven[jk - 1] - zdmfde * pvd[jk - 1];
        pud[jk] = zmfduk * (T(1.) / c_min(-T(C::cmfcmin), pmfd[jk]));
        pvd[jk] = zmfdvk * (T(1.) / c_min(-T(C::cmfcmin), pmfd[jk]));
      }
    }
  }
}

// cuflx: final convective fluxes, precipitation, melting and evaporation. WRF's ktopm2 is the row's
// minimum cloud top - 2; the levels above this column's own top only add zeros, so its own top is used.
template <typename T>
CPPWRF_HD inline void cuflx(const TdConst<T>& c, int klev, const T* pqen, const T* pqsen, const T* ptenh,
                            const T* pqenh, const T* paph, const T* pgeoh, int kcbot, int kctop, int kdtop,
                            int& ktype, bool& lddraf, bool& ldcum, T* pmfu, T* pmfd, T* pmfus, T* pmfds, T* pmfuq,
                            T* pmfdq, T* pmful, T* plude, T* pdmfup, T* pdmfdp, T& prfl, T& prain, const T* pten,
                            T& psfl, T* pdpmel, int& ktopm2, T ztmst, const T* sig1) {
  using C = TdConst<T>;
  const int klevp1 = klev + 1;
  const T zcons1 = c.cpd / (c.alf * c.g * ztmst);
  const T zcons2 = T(1.) / (c.g * ztmst);
  const T zcucov = T(0.05);
  const T ztmelp2 = T(C::tmelt) + T(2.);
  int itop = klev;
  prfl = T(0.);
  psfl = T(0.);
  prain = T(0.);
  if (!C::lmfscv && ktype == 2) {
    ldcum = false;
    lddraf = false;
  }
  itop = c_min(itop, kctop);
  if (!ldcum || kdtop < kctop) lddraf = false;
  if (!ldcum) ktype = 0;
  ktopm2 = itop - 2;
  for (int jk = ktopm2; jk <= klev; ++jk) {
    if (ldcum && jk >= kctop - 1) {
      pmfus[jk] = pmfus[jk] - pmfu[jk] * (c.cpd * ptenh[jk] + pgeoh[jk]);
      pmfuq[jk] = pmfuq[jk] - pmfu[jk] * pqenh[jk];
      if (lddraf && jk >= kdtop) {
        pmfds[jk] = pmfds[jk] - pmfd[jk] * (c.cpd * ptenh[jk] + pgeoh[jk]);
        pmfdq[jk] = pmfdq[jk] - pmfd[jk] * pqenh[jk];
      } else {
        pmfd[jk] = T(0.);
        pmfds[jk] = T(0.);
        pmfdq[jk] = T(0.);
        pdmfdp[jk - 1] = T(0.);
      }
    } else {
      pmfu[jk] = T(0.);
      pmfd[jk] = T(0.);
      pmfus[jk] = T(0.);
      pmfds[jk] = T(0.);
      pmfuq[jk] = T(0.);
      pmfdq[jk] = T(0.);
      pmful[jk] = T(0.);
      pdmfup[jk - 1] = T(0.);
      pdmfdp[jk - 1] = T(0.);
      plude[jk - 1] = T(0.);
    }
  }
  for (int jk = ktopm2; jk <= klev; ++jk) {
    if (ldcum && jk > kcbot) {
      const int ikb = kcbot;
      T zzp = ((paph[klevp1] - paph[jk]) / (paph[klevp1] - paph[ikb]));
      if (ktype == 3) zzp = zzp * zzp;
      pmfu[jk] = pmfu[ikb] * zzp;
      pmfus[jk] = pmfus[ikb] * zzp;
      pmfuq[jk] = pmfuq[ikb] * zzp;
      pmful[jk] = pmful[ikb] * zzp;
    }
    // rain / snow fall rates, melting, evaporation
    if (ldcum) {
      prain = prain + pdmfup[jk];
      if (pten[jk] > T(C::tmelt)) {
        prfl = prfl + pdmfup[jk] + pdmfdp[jk];
        if (psfl > T(0.) && pten[jk] > ztmelp2) {
          const T zfac = zcons1 * (paph[jk + 1] - paph[jk]);
          const T zsnmlt = c_min(psfl, zfac * (pten[jk] - ztmelp2));
          pdpmel[jk] = zsnmlt;
          psfl = psfl - zsnmlt;
          prfl = prfl + zsnmlt;
        }
      } else {
        psfl = psfl + pdmfup[jk] + pdmfdp[jk];
      }
    }
  }
  prfl = c_max(prfl, T(0.));
  psfl = c_max(psfl, T(0.));
  T zpsubcl = prfl + psfl;
  for (int jk = ktopm2; jk <= klev; ++jk) {
    if (ldcum && jk >= kcbot && zpsubcl > T(1.e-20)) {
      const T zrfl = zpsubcl;
      const T cevapcu = c.cevapcu1 * std::sqrt(c.cevapcu2 * std::sqrt(sig1[jk]));
      const T m = c_max(T(0.), std::sqrt(zrfl / zcucov) - cevapcu * (paph[jk + 1] - paph[jk]) * c_max(T(0.), pqsen[jk] - pqen[jk]));
      T zrnew = (m * m) * zcucov;
      const T zrmin = zrfl - zcucov * c_max(T(0.), T(0.8) * pqsen[jk] - pqen[jk]) * zcons2 * (paph[jk + 1] - paph[jk]);
      zrnew = c_max(zrnew, zrmin);
      const T zrfln = c_max(zrnew, T(0.));
      const T zdrfl = c_min(T(0.), zrfln - zrfl);
      pdmfup[jk] = pdmfup[jk] + zdrfl;
      zpsubcl = zrfln;
    }
  }
  const T zdpevap = zpsubcl - (prfl + psfl);
  prfl = prfl + zdpevap * prfl * (T(1.) / c_max(T(1.e-20), prfl + psfl));
  psfl = psfl + zdpevap * psfl * (T(1.) / c_max(T(1.e-20), prfl + psfl));
}

// cudtdq: t and q tendencies, surface precipitation rates (the global budgets are diagnostics only)
template <typename T>
CPPWRF_HD inline void cudtdq(const TdConst<T>& c, int klev, int ktopm2, const T* paph, bool ldcum, const T* pten,
                             T* ptte, T* pqte, const T* pmfus, const T* pmfds, const T* pmfuq, const T* pmfdq,
                             const T* pmful, const T* pdmfup, const T* pdmfdp, const T* pdpmel, T prfl, T psfl,
                             T& prsfc, T& pssfc, const T* pqen, const T* pqsen, const T* plude, T* pcte) {
  using C = TdConst<T>;
  const T rhc = T(C::rhc), rhm = T(C::rhm), fdbk = T(C::fdbk);
  for (int jk = ktopm2; jk <= klev; ++jk) {
    if (!ldcum) continue;
    const T zalv = pten[jk] > T(C::tmelt) ? c.alv : c.als;
    const T rhk = c_min(T(1.0), pqen[jk] / pqsen[jk]);
    const T rhcoe = c_max(T(0.0), (rhk - rhc) / (rhm - rhc));
    const T pldfd = c_max(T(0.0), rhcoe * fdbk * plude[jk]);
    if (jk < klev) {
      const T zdtdt = (c.g / (paph[jk + 1] - paph[jk])) * c.rcpd *
                      (pmfus[jk + 1] - pmfus[jk] + pmfds[jk + 1] - pmfds[jk] - c.alf * pdpmel[jk] -
                       zalv * (pmful[jk + 1] - pmful[jk] - pldfd - (pdmfup[jk] + pdmfdp[jk])));
      ptte[jk] = ptte[jk] + zdtdt;
      const T zdqdt = (c.g / (paph[jk + 1] - paph[jk])) *
                      (pmfuq[jk + 1] - pmfuq[jk] + pmfdq[jk + 1] - pmfdq[jk] + pmful[jk + 1] - pmful[jk] - pldfd -
                       (pdmfup[jk] + pdmfdp[jk]));
      pqte[jk] = pqte[jk] + zdqdt;
      pcte[jk] = (c.g / (paph[jk + 1] - paph[jk])) * pldfd;
    } else {
      const T zdtdt = -(c.g / (paph[jk + 1] - paph[jk])) * c.rcpd *
                      (pmfus[jk] + pmfds[jk] + c.alf * pdpmel[jk] - zalv * (pmful[jk] + pdmfup[jk] + pdmfdp[jk] + pldfd));
      ptte[jk] = ptte[jk] + zdtdt;
      const T zdqdt = -(c.g / (paph[jk + 1] - paph[jk])) *
                      (pmfuq[jk] + pmfdq[jk] + pldfd + (pmful[jk] + pdmfup[jk] + pdmfdp[jk]));
      pqte[jk] = pqte[jk] + zdqdt;
      pcte[jk] = (c.g / (paph[jk + 1] - paph[jk])) * pldfd;
    }
  }
  prsfc = prfl;
  pssfc = psfl;
}

// cududv: u and v tendencies (cumulus friction)
template <typename T>
CPPWRF_HD inline void cududv(const TdConst<T>& c, int klev, int ktopm2, int ktype, int kcbot, const T* paph,
                             bool ldcum, const T* puen, const T* pven, T* pvom, T* pvol, const T* puu, const T* pud,
                             const T* pvu, const T* pvd, const T* pmfu, const T* pmfd) {
  if (!ldcum) return;
  const int klevp1 = klev + 1;
  T zmfuu[KM], zmfdu[KM], zmfuv[KM], zmfdv[KM];
  for (int jk = ktopm2; jk <= klev; ++jk) {
    const int ik = jk - 1;
    zmfuu[jk] = pmfu[jk] * (puu[jk] - puen[ik]);
    zmfuv[jk] = pmfu[jk] * (pvu[jk] - pven[ik]);
    zmfdu[jk] = pmfd[jk] * (pud[jk] - puen[ik]);
    zmfdv[jk] = pmfd[jk] * (pvd[jk] - pven[ik]);
  }
  for (int jk = ktopm2; jk <= klev; ++jk) {
    if (jk > kcbot) {
      const int ikb = kcbot;
      T zzp = ((paph[klevp1] - paph[jk]) / (paph[klevp1] - paph[ikb]));
      if (ktype == 3) zzp = zzp * zzp;
      zmfuu[jk] = zmfuu[ikb] * zzp;
      zmfuv[jk] = zmfuv[ikb] * zzp;
      zmfdu[jk] = zmfdu[ikb] * zzp;
      zmfdv[jk] = zmfdv[ikb] * zzp;
    }
  }
  for (int jk = ktopm2; jk <= klev; ++jk) {
    if (jk < klev) {
      const T zdudt = (c.g / (paph[jk + 1] - paph[jk])) * (zmfuu[jk + 1] - zmfuu[jk] + zmfdu[jk + 1] - zmfdu[jk]);
      const T zdvdt = (c.g / (paph[jk + 1] - paph[jk])) * (zmfuv[jk + 1] - zmfuv[jk] + zmfdv[jk + 1] - zmfdv[jk]);
      pvom[jk] = pvom[jk] + zdudt;
      pvol[jk] = pvol[jk] + zdvdt;
    } else {
      const T zdudt = -(c.g / (paph[jk + 1] - paph[jk])) * (zmfuu[jk] + zmfdu[jk]);
      const T zdvdt = -(c.g / (paph[jk + 1] - paph[jk])) * (zmfuv[jk] + zmfdv[jk]);
      pvom[jk] = pvom[jk] + zdudt;
      pvol[jk] = pvol[jk] + zdvdt;
    }
  }
}

// cumastr_new state of one column (the level arrays of the master routine)
template <typename T>
struct TdColumn {
  T ztenh[KM], zqenh[KM], zgeoh[KM], zqsenh[KM], ztd[KM], zqd[KM], zmfus[KM], zmfds[KM], zmfuq[KM], zmfdq[KM];
  T zdmfup[KM], zdmfdp[KM], zmful[KM], zuu[KM], zvu[KM], zud[KM], zvd[KM], zdpmel[KM], zhhatt[KM];
  T ptu[KM], pqu[KM], plu[KM], plude[KM], pmfu[KM], pmfd[KM];
  int ilab[KM];
};

// cumastr_new up to and including the first cloud ascent; returns icum (1 if the column convects)
template <typename T>
CPPWRF_HD inline int cumastr_first(const TdConst<T>& c, int klev, const T* pten, const T* pqen, const T* puen,
                                   const T* pven, const T* pverv, const T* pqsen, T pqhfl, T ztmst, const T* pap,
                                   const T* paph, const T* pgeo, const T* pqte, bool& ldcum, int& ktype, int& kcbot,
                                   int& kctop, int lndj, TdColumn<T>& s, int& ilwmin, int& ictop0, int& ihmin,
                                   T& zhcbase, T& zmfub, T& zentr, T& zdqpbl) {
  using C = TdConst<T>;
  const int klevm1 = klev - 1;
  const T zcons2 = T(1.) / (c.g * ztmst);
  cuini(c, klev, pten, pqen, pqsen, puen, pven, pverv, pgeo, paph, s.zgeoh, s.ztenh, s.zqenh, s.zqsenh, ilwmin, s.ptu,
        s.pqu, s.ztd, s.zqd, s.zuu, s.zvu, s.zud, s.zvd, s.pmfu, s.pmfd, s.zmfus, s.zmfds, s.zmfuq, s.zmfdq, s.zdmfup,
        s.zdmfdp, s.zdpmel, s.plu, s.plude, s.ilab);
  cubase(c, klev, s.ztenh, s.zqenh, s.zgeoh, paph, s.ptu, s.pqu, s.plu, puen, pven, s.zuu, s.zvu, ldcum, kcbot, s.ilab);
  // (b) total moisture convergence, type of convection
  int jk = 1;
  T zdqcv = pqte[jk] * (paph[jk + 1] - paph[jk]);
  zdqpbl = T(0.0);
  for (jk = 2; jk <= klev; ++jk) {
    zdqcv = zdqcv + pqte[jk] * (paph[jk + 1] - paph[jk]);
    if (jk >= kcbot) zdqpbl = zdqpbl + pqte[jk] * (paph[jk + 1] - paph[jk]);
  }
  ktype = 0;
  ktype = zdqcv > c_max(T(0.), T(1.1) * pqhfl * c.g) ? 1 : 2;
  // (c) moisture supply of the boundary layer, cloud base massflux without downdrafts
  int ikb = kcbot;
  const T zqumqe = s.pqu[ikb] + s.plu[ikb] - s.zqenh[ikb];
  const T zdqmin = c_max(T(0.01) * s.zqenh[ikb], T(1.e-10));
  if (zdqpbl > T(0.) && zqumqe > zdqmin && ldcum) {
    zmfub = zdqpbl / (c.g * c_max(zqumqe, zdqmin));
  } else {
    zmfub = T(0.01);
    ldcum = false;
  }
  const T zmfmax = (paph[ikb] - paph[ikb - 1]) * zcons2;
  zmfub = c_min(zmfub, zmfmax);
  // 4.0 (a) cloud height estimate for entrainment / detrainment
  ikb = kcbot;
  zhcbase = c.cpd * s.ptu[ikb] + s.zgeoh[ikb] + c.alv * s.pqu[ikb];
  ictop0 = kcbot - 1;
  const T zalvdcp = c.alv / c.cpd;
  const T zqalv = T(1.) / c.alv;
  for (jk = klevm1; jk >= 3; --jk) {
    const T zhsat = c.cpd * s.ztenh[jk] + s.zgeoh[jk] + c.alv * s.zqsenh[jk];
    const T d = s.ztenh[jk] - T(C::c4les);
    const T zgam = c.c5les * zalvdcp * s.zqsenh[jk] / ((T(1.) - c.vtmpc1 * s.zqsenh[jk]) * (d * d));
    const T zzz = c.cpd * s.ztenh[jk] * T(0.608);
    const T zhhat = zhsat - (zzz + zgam * zzz) / (T(1.) + zgam * zzz * zqalv) * c_max(s.zqsenh[jk] - s.zqenh[jk], T(0.));
    s.zhhatt[jk] = zhhat;
    if (jk < ictop0 && zhcbase > zhhat) ictop0 = jk;
  }
  {
    jk = kcbot;
    const T zhsat = c.cpd * s.ztenh[jk] + s.zgeoh[jk] + c.alv * s.zqsenh[jk];
    const T d = s.ztenh[jk] - T(C::c4les);
    const T zgam = c.c5les * zalvdcp * s.zqsenh[jk] / ((T(1.) - c.vtmpc1 * s.zqsenh[jk]) * (d * d));
    const T zzz = c.cpd * s.ztenh[jk] * T(0.608);
    const T zhhat = zhsat - (zzz + zgam * zzz) / (T(1.) + zgam * zzz * zqalv) * c_max(s.zqsenh[jk] - s.zqenh[jk], T(0.));
    s.zhhatt[jk] = zhhat;
  }
  // lowest possible organized detrainment level
  T zhmin = T(0.);
  ihmin = ldcum && ktype == 1 ? kcbot : -1;
  const T zbi = T(1.) / (T(25.) * c.g);
  for (jk = klev; jk >= 1; --jk) {
    const bool llo1 = ldcum && ktype == 1 && ihmin == kcbot;
    if (llo1 && jk < kcbot && jk >= ictop0) {
      ikb = kcbot;
      const T zro = c.rd * s.ztenh[jk] / (c.g * paph[jk]);
      const T zdz = (paph[jk] - paph[jk - 1]) * zro;
      const T zdhdz = (c.cpd * (pten[jk - 1] - pten[jk]) + c.alv * (pqen[jk - 1] - pqen[jk]) + (pgeo[jk - 1] - pgeo[jk])) *
                      c.g / (pgeo[jk - 1] - pgeo[jk]);
      const T zdepth = s.zgeoh[jk] - s.zgeoh[ikb];
      const T zfac = std::sqrt(T(1.) + zdepth * zbi);
      zhmin = zhmin + zdhdz * zfac * zdz;
      const T zrh = -c.alv * (s.zqsenh[jk] - s.zqenh[jk]) * zfac;
      if (zhmin > zrh) ihmin = jk;
    }
  }
  if (ldcum && ktype == 1) {
    if (ihmin < ictop0) ihmin = ictop0;
  }
  zentr = ktype == 1 ? T(C::entrpen) : T(C::entrscv);
  if (lndj == 1) zentr = zentr * T(1.1);
  // (b) ascent without downdrafts
  int icum = 0;
  cuasc_new(c, klev, s.ztenh, s.zqenh, puen, pven, pten, pqen, pqsen, pgeo, s.zgeoh, pap, paph, pverv, ilwmin, ldcum,
            zhcbase, ktype, s.ilab, s.ptu, s.pqu, s.plu, s.zuu, s.zvu, s.pmfu, zmfub, zentr, s.zmfus, s.zmfuq, s.zmful,
            s.plude, s.zdmfup, kcbot, kctop, ictop0, icum, ztmst, ihmin, s.zhhatt, s.zqsenh);
  return icum;
}

// cumastr_new (row_icum: the row's icum after the first ascent)
template <typename T>
CPPWRF_HD inline void cumastr_new(const TdConst<T>& c, int klev, const T* pten, const T* pqen, const T* puen,
                                  const T* pven, const T* pverv, const T* pqsen, T pqhfl, T ztmst, const T* pap,
                                  const T* paph, const T* pgeo, T* ptte, T* pqte, T* pvom, T* pvol, T& prsfc,
                                  T& pssfc, int& ktype, T* pcte, const T* sig1, int lndj, int row_icum) {
  using C = TdConst<T>;
  const T zcons2 = T(1.) / (c.g * ztmst);
  TdColumn<T> s;
  bool ldcum = false;
  int kcbot = 0, kctop = 0, ilwmin = 0, ictop0 = 0, ihmin = 0, idtop = 0, itopm2 = 0, icum = 0;
  T zhcbase = T(0), zmfub = T(0), zentr = T(0), zdqpbl = T(0);
  cumastr_first(c, klev, pten, pqen, puen, pven, pverv, pqsen, pqhfl, ztmst, pap, paph, pgeo, pqte, ldcum, ktype, kcbot,
                kctop, lndj, s, ilwmin, ictop0, ihmin, zhcbase, zmfub, zentr, zdqpbl);
  idtop = 0;
  if (row_icum == 0) return;
  // (c) cloud depth, entrainment rate, precipitation rate for the downdrafts
  const T zpbmpt = paph[kcbot] - paph[kctop];
  if (ldcum) ictop0 = kctop;
  if (ldcum && ktype == 1 && zpbmpt < T(C::zdnoprc)) ktype = 2;
  if (ktype == 2) {
    zentr = T(C::entrscv);
    if (lndj == 1) zentr = zentr * T(1.1);
  }
  T zrfl = s.zdmfup[1];
  for (int jk = 2; jk <= klev; ++jk) zrfl = zrfl + s.zdmfup[jk];
  // 5.0 downdrafts
  bool loddraf = false;
  if (C::lmfdd) {
    cudlfs(c, klev, s.ztenh, s.zqenh, puen, pven, s.zgeoh, paph, s.ptu, s.pqu, s.zuu, s.zvu, ldcum, kcbot, kctop, zmfub,
           zrfl, s.ztd, s.zqd, s.zud, s.zvd, s.pmfd, s.zmfds, s.zmfdq, s.zdmfdp, idtop, loddraf);
    cuddraf(c, klev, s.ztenh, s.zqenh, puen, pven, s.zgeoh, paph, zrfl, loddraf, s.ztd, s.zqd, s.zud, s.zvd, s.pmfd,
            s.zmfds, s.zmfdq, s.zdmfdp);
  }
  // 5.1 cloud base massflux from a CAPE closure (ktype 1)
  T zheat = T(0.0), zcape = T(0.0), zrelh = T(0.0), zmfub1 = zmfub;
  if (ldcum && ktype == 1) {
    const int ktop0 = c_max(12, kctop);
    const int ikb = kcbot;
    for (int jk = 2; jk <= klev; ++jk) {
      if (jk <= kcbot && jk > kctop) {
        const T zro = paph[jk] / (c.rd * s.ztenh[jk]);
        const T zdz = (paph[jk] - paph[jk - 1]) / (c.g * zro);
        zheat = zheat + ((pten[jk - 1] - pten[jk] + c.g * zdz / c.cpd) / s.ztenh[jk] + T(0.608) * (pqen[jk - 1] - pqen[jk])) *
                            (s.pmfu[jk] + s.pmfd[jk]) * c.g / zro;
        zcape = zcape + c.g * ((s.ptu[jk] * (T(1.) + T(.608) * s.pqu[jk] - s.plu[jk])) / (s.ztenh[jk] * (T(1.) + T(.608) * s.zqenh[jk])) - T(1.0)) * zdz;
      }
      if (jk <= kcbot && jk > ktop0) {
        const T dept = (paph[jk + 1] - paph[jk]) / (paph[ikb + 1] - paph[ktop0 + 1]);
        zrelh = zrelh + dept * pqen[jk] / pqsen[jk];
      }
    }
    if (zrelh >= T(C::crirh)) {
      const T zht = c_max(T(0.0), (zcape - T(0.0))) / (T(C::ztau) * zheat);
      zmfub1 = c_max(zmfub * zht, T(0.01));
      const T zmfmax = (paph[ikb] - paph[ikb - 1]) * zcons2;
      zmfub1 = c_min(zmfub1, zmfmax);
    } else {
      zmfub1 = T(0.01);
      zmfub = T(0.01);
      ldcum = false;
    }
  }
  // 5.2 convective fluxes with the downdraft effect on the boundary layer moisture budget
  if (ktype != 1) {
    const int ikb = kcbot;
    const T zeps = s.pmfd[ikb] < T(0.0) && loddraf ? T(C::cmfdeps) : T(0.);
    const T zqumqe = s.pqu[ikb] + s.plu[ikb] - zeps * s.zqd[ikb] - (T(1.) - zeps) * s.zqenh[ikb];
    const T zdqmin = c_max(T(0.01) * s.zqenh[ikb], T(1.e-10));
    const T zmfmax = (paph[ikb] - paph[ikb - 1]) * zcons2;
    if (zdqpbl > T(0.) && zqumqe > zdqmin && ldcum && zmfub < zmfmax) {
      zmfub1 = zdqpbl / (c.g * c_max(zqumqe, zdqmin));
    } else {
      zmfub1 = zmfub;
    }
    const bool llo1 = (ktype == 2) && std::fabs(zmfub1 - zmfub) < T(0.2) * zmfub;
    if (!llo1) zmfub1 = zmfub;
    zmfub1 = c_min(zmfub1, zmfmax);
  }
  for (int jk = 1; jk <= klev; ++jk) {
    if (ldcum) {
      const T zfac = zmfub1 / c_max(zmfub, T(1.e-10));
      s.pmfd[jk] = s.pmfd[jk] * zfac;
      s.zmfds[jk] = s.zmfds[jk] * zfac;
      s.zmfdq[jk] = s.zmfdq[jk] * zfac;
      s.zdmfdp[jk] = s.zdmfdp[jk] * zfac;
    } else {
      s.pmfd[jk] = T(0.0);
      s.zmfds[jk] = T(0.0);
      s.zmfdq[jk] = T(0.0);
      s.zdmfdp[jk] = T(0.0);
    }
  }
  zmfub = ldcum ? zmfub1 : T(0.0);
  // 6.0 final cloud ascent
  cuasc_new(c, klev, s.ztenh, s.zqenh, puen, pven, pten, pqen, pqsen, pgeo, s.zgeoh, pap, paph, pverv, ilwmin, ldcum,
            zhcbase, ktype, s.ilab, s.ptu, s.pqu, s.plu, s.zuu, s.zvu, s.pmfu, zmfub, zentr, s.zmfus, s.zmfuq, s.zmful,
            s.plude, s.zdmfup, kcbot, kctop, ictop0, icum, ztmst, ihmin, s.zhhatt, s.zqsenh);
  // 7.0 final convective fluxes
  T prain = T(0), zsfl = T(0);
  cuflx(c, klev, pqen, pqsen, s.ztenh, s.zqenh, paph, s.zgeoh, kcbot, kctop, idtop, ktype, loddraf, ldcum, s.pmfu, s.pmfd,
        s.zmfus, s.zmfds, s.zmfuq, s.zmfdq, s.zmful, s.plude, s.zdmfup, s.zdmfdp, zrfl, prain, pten, zsfl, s.zdpmel,
        itopm2, ztmst, sig1);
  // 8.0 t and q tendencies
  cudtdq(c, klev, itopm2, paph, ldcum, pten, ptte, pqte, s.zmfus, s.zmfds, s.zmfuq, s.zmfdq, s.zmful, s.zdmfup,
         s.zdmfdp, s.zdpmel, zrfl, zsfl, prsfc, pssfc, pqen, pqsen, s.plude, pcte);
  // 9.0 u and v tendencies
  if (C::lmfdudv)
    cududv(c, klev, itopm2, ktype, kcbot, paph, ldcum, puen, pven, pvom, pvol, s.zuu, s.zud, s.zvu, s.zvd, s.pmfu,
           s.pmfd);
}

// tiecnv inputs of one column: pt, pqv (mixing ratio) ... at levels 1 (top) .. km, paph at 1 .. km + 1
template <typename T>
struct TdInput {
  T pu[KM], pv[KM], pt[KM], pqv[KM], pqc[KM], pqi[KM], pqvf[KM], pqvbl[KM], poz[KM], pomg[KM], pap[KM], paph[KM];
  T sig1[KM];
};

// tiecnv's conversion of the model variables (shared by both passes)
template <typename T>
CPPWRF_HD inline void tiecnv_prep(const TdConst<T>& c, int km, const TdInput<T>& in, int lane, int nlane, T* ztp1,
                                  T* zqp1, T* zqsat, T* pgeo, T* pqte, T* zqq) {
  for (int k = 1; k <= km; ++k) {
    ztp1[k] = in.pt[k];
    zqp1[k] = in.pqv[k] / (T(1.0) + in.pqv[k]);
    pgeo[k] = c.g * in.poz[k];
    const T tt = ztp1[k];
    zqsat[k] = tlucua(c, tt, lane, nlane) / in.pap[k];  // vectorised loop: libmvec exp
    zqsat[k] = c_min(T(0.5), zqsat[k]);
    zqsat[k] = zqsat[k] / (T(1.) - c.vtmpc1 * zqsat[k]);
    pqte[k] = in.pqvf[k] + in.pqvbl[k];
    zqq[k] = pqte[k];
  }
}

// first pass: the column's icum after cumastr_new's first cloud ascent
template <typename T>
CPPWRF_HD inline int tiedtke_first_ascent(const TdConst<T>& c, int km, const TdInput<T>& in, T evap, int lndj, T dt,
                                          int lane, int nlane) {
  T ztp1[KM], zqp1[KM], zqsat[KM], pgeo[KM], pqte[KM], zqq[KM];
  tiecnv_prep(c, km, in, lane, nlane, ztp1, zqp1, zqsat, pgeo, pqte, zqq);
  TdColumn<T> s;
  bool ldcum = false;
  int ktype = 0, kcbot = 0, kctop = 0, ilwmin = 0, ictop0 = 0, ihmin = 0;
  T zhcbase = T(0), zmfub = T(0), zentr = T(0), zdqpbl = T(0);
  return cumastr_first(c, km, ztp1, zqp1, in.pu, in.pv, in.pomg, zqsat, evap, dt, in.pap, in.paph, pgeo, pqte, ldcum,
                       ktype, kcbot, kctop, lndj, s, ilwmin, ictop0, ihmin, zhcbase, zmfub, zentr, zdqpbl);
}

// tiecnv: updates pu, pv, pt, pqv, pqc, pqi of the column (levels 1 .. km) and returns zprecc
template <typename T>
CPPWRF_HD inline T tiecnv(const TdConst<T>& c, int km, const TdInput<T>& in, T* pu, T* pv, T* pt, T* pqv, T* pqc,
                          T* pqi, T evap, int lndj, T dt, int lane, int nlane, int row_icum) {
  using C = TdConst<T>;
  const T ztmst = dt;
  T ztp1[KM], zqp1[KM], zqsat[KM], pgeo[KM], pqte[KM], zqq[KM], ptte[KM], pcte[KM], pvom[KM], pvol[KM];
  T prsfc = T(0.0), pssfc = T(0.0);
  const T pqhfl = evap;
  for (int k = 1; k <= km; ++k) {
    ptte[k] = T(0.0);
    pcte[k] = T(0.0);
    pvom[k] = T(0.0);
    pvol[k] = T(0.0);
  }
  tiecnv_prep(c, km, in, lane, nlane, ztp1, zqp1, zqsat, pgeo, pqte, zqq);
  int ktype = 0;
  cumastr_new(c, km, ztp1, zqp1, in.pu, in.pv, in.pomg, zqsat, pqhfl, ztmst, in.pap, in.paph, pgeo, ptte, pqte, pvom,
              pvol, prsfc, pssfc, ktype, pcte, in.sig1, lndj, row_icum);
  // cloud water and cloud ice detrained from convection
  if (C::fdbk >= 1.0e-9) {
    for (int k = 1; k <= km; ++k) {
      if (pcte[k] > T(0.0)) {
        const T ztpp1 = pt[k] + ptte[k] * ztmst;
        T fliq, zalf;
        if (ztpp1 >= T(C::t000)) {
          fliq = T(1.0);
          zalf = T(0.0);
        } else if (ztpp1 <= T(C::hgfr)) {
          fliq = T(0.0);
          zalf = c.alf;
        } else {
          const T ztc = ztpp1 - T(C::t000);
          fliq = T(0.0059) + T(0.9941) * std::exp(-T(0.003102) * ztc * ztc);
          zalf = c.alf;
        }
        const T fice = T(1.0) - fliq;
        pqc[k] = pqc[k] + fliq * pcte[k] * ztmst;
        pqi[k] = pqi[k] + fice * pcte[k] * ztmst;
        ptte[k] = ptte[k] - zalf * c.rcpd * fliq * pcte[k];
      }
    }
  }
  for (int k = 1; k <= km; ++k) {
    pt[k] = ztp1[k] + ptte[k] * ztmst;
    zqp1[k] = zqp1[k] + (pqte[k] - zqq[k]) * ztmst;
    pqv[k] = zqp1[k] / (T(1.0) - zqp1[k]);
  }
  const T zprecc = c_max(T(0.0), (prsfc + pssfc) * ztmst);
  if (C::lmfdudv) {
    for (int k = 1; k <= km; ++k) {
      pu[k] = pu[k] + pvom[k] * ztmst;
      pv[k] = pv[k] + pvol[k] * ztmst;
    }
  }
  return zprecc;
}

}  // namespace cppwrf::phys::cu_tiedtke
