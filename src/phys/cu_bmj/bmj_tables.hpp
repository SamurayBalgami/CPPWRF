// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// BMJINIT's look-up tables (WRF 4.8.0 phys/module_cu_bmj.F:1813-2077) and SPLINE (:2081-2189), on the host. The
// statements and their order are WRF's (TH = TH + DTH and P = P + DP accumulated as there); gfortran compiles these
// loops without libmvec (scalar exp and pow).
#include <vector>

#include "core/libm.hpp"
#include "phys/cu_bmj/bmj.hpp"

namespace cppwrf::phys::cu_bmj {

// SPLINE(JTBX, NOLD, XOLD, YOLD, Y2, NNEW, XNEW, YNEW, P, Q): the arrays 1-based (index 0 unused)
template <typename T>
inline void spline(int nold, const T* xold, const T* yold, T* y2, int nnew, const T* xnew, T* ynew, T* p, T* q) {
  const int noldm1 = nold - 1;
  T dxl = xold[2] - xold[1];
  T dxr = xold[3] - xold[2];
  T dydxl = (yold[2] - yold[1]) / dxl;
  T dydxr = (yold[3] - yold[2]) / dxr;
  const T rtdxc = T(0.5) / (dxl + dxr);
  p[1] = rtdxc * (T(6.) * (dydxr - dydxl) - dxl * y2[1]);
  q[1] = -rtdxc * dxr;
  int k;
  if (nold != 3) {
    k = 3;
    do {  // 100
      dxl = dxr;
      dydxl = dydxr;
      dxr = xold[k + 1] - xold[k];
      dydxr = (yold[k + 1] - yold[k]) / dxr;
      const T dxc = dxl + dxr;
      const T den = T(1.) / (dxl * q[k - 2] + dxc + dxc);
      p[k - 1] = den * (T(6.) * (dydxr - dydxl) - dxl * p[k - 2]);
      q[k - 1] = -den * dxr;
      k = k + 1;
    } while (k < nold);
  }
  k = noldm1;  // 150
  do {         // 200
    y2[k] = p[k - 1] + q[k - 1] * y2[k + 1];
    k = k - 1;
  } while (k > 1);
  T ak = T(0.), bk = T(0.), ck = T(0.);
  int k1 = 1;
  do {  // 300
    const T xk = xnew[k1];
    int kold = 0;
    bool found = false;
    for (int k2 = 2; k2 <= nold; ++k2) {
      if (xold[k2] > xk) {
        kold = k2 - 1;
        found = true;
        break;
      }
    }
    if (!found) {
      ynew[k1] = yold[nold];
    } else {
      if (k1 == 1 || k != kold) {  // 450 -> 500 (k1 == 1 or k /= kold), else 550
        k = kold;
        const T y2k = y2[k];
        const T y2kp1 = y2[k + 1];
        const T dx = xold[k + 1] - xold[k];
        const T rdx = T(1.) / dx;
        ak = T(.1666667) * rdx * (y2kp1 - y2k);
        bk = T(0.5) * y2k;
        ck = rdx * (yold[k + 1] - yold[k]) - T(.1666667) * dx * (y2kp1 + y2k + y2k);
      }
      const T x = xk - xold[k];  // 550
      const T xsq = x * x;
      ynew[k1] = ak * xsq * x + bk * xsq + ck * x + yold[k];
    }
    k1 = k1 + 1;  // 600
  } while (k1 <= nnew);
}

// BMJINIT's tables for CP, RD (the cumulus driver's cp and r_d)
template <typename T>
inline void init_bmj_tables(BmjTables<T>& tb, T cp, T rd) {
  using C = BmjC<T>;
  constexpr T EPS = T(1.E-9);
  constexpr int N = JTBQ + 1;  // 1-based work arrays (JTB / JTBQ long in WRF)
  std::vector<T> w(std::size_t(21) * N, T(0.));  // the work arrays (one block per array)
  auto a = [&](int n) { return w.data() + std::size_t(n) * N; };
  T *app = a(0), *apt = a(1), *aqp = a(2), *aqt = a(3), *pnew = a(4), *pold = a(5), *qsnew = a(6), *qsold = a(7),
    *thenew = a(8), *theold = a(9), *tnew = a(10), *told = a(11), *y2p = a(12), *y2t = a(13);
  T *aptq = a(14), *aqtq = a(15), *thenewq = a(16), *theoldq = a(17), *tnewq = a(18), *toldq = a(19), *y2tq = a(20);
  const T elocp = C::ELIWV / cp;
  //----------------COARSE LOOK-UP TABLE FOR SATURATION POINT--------------
  int kthm = JTB, kpm = ITB;
  int kthm1 = kthm - 1, kpm1 = kpm - 1;
  T dth = (C::THH - C::THL) / T(kthm - 1);
  T dp = (C::PH - C::PL) / T(kpm - 1);
  T th = C::THL - dth;
  T p, ape, denom, qs;
  for (int kth = 1; kth <= kthm; ++kth) {
    th = th + dth;
    p = C::PL - dp;
    for (int kp = 1; kp <= kpm; ++kp) {
      p = p + dp;
      ape = libm::pow(T(100000.) / p, rd / cp);
      denom = th - C::A4 * ape;
      if (denom > EPS) {
        qsold[kp] = C::PQ0 / p * libm::exp(C::A2 * (th - C::A3 * ape) / denom);
      } else {
        qsold[kp] = T(0.);
      }
      pold[kp] = p;
    }
    const T qs0k = qsold[1];
    const T sqsk = qsold[kpm] - qsold[1];
    qsold[1] = T(0.);
    qsold[kpm] = T(1.);
    for (int kp = 2; kp <= kpm1; ++kp) {
      qsold[kp] = (qsold[kp] - qs0k) / sqsk;
      if ((qsold[kp] - qsold[kp - 1]) < EPS) qsold[kp] = qsold[kp - 1] + EPS;
    }
    tb.qs0[kth - 1] = qs0k;
    tb.sqs[kth - 1] = sqsk;
    qsnew[1] = T(0.);
    qsnew[kpm] = T(1.);
    const T dqs = T(1.) / T(kpm - 1);
    for (int kp = 2; kp <= kpm1; ++kp) qsnew[kp] = qsnew[kp - 1] + dqs;
    y2p[1] = T(0.);
    y2p[kpm] = T(0.);
    spline(kpm, qsold, pold, y2p, kpm, qsnew, pnew, app, aqp);
    for (int kp = 1; kp <= kpm; ++kp) tb.ptbl[kth - 1][kp - 1] = pnew[kp];
  }
  //------------COARSE LOOK-UP TABLE FOR T(P) FROM CONSTANT THE------------
  p = C::PL - dp;
  for (int kp = 1; kp <= kpm; ++kp) {
    p = p + dp;
    th = C::THL - dth;
    for (int kth = 1; kth <= kthm; ++kth) {
      th = th + dth;
      ape = libm::pow(T(1.E5) / p, rd / cp);
      denom = th - C::A4 * ape;
      if (denom > EPS) {
        qs = C::PQ0 / p * libm::exp(C::A2 * (th - C::A3 * ape) / denom);
      } else {
        qs = T(0.);
      }
      told[kth] = th / ape;
      theold[kth] = th * libm::exp(elocp * qs / told[kth]);
    }
    const T the0k = theold[1];
    const T sthek = theold[kthm] - theold[1];
    theold[1] = T(0.);
    theold[kthm] = T(1.);
    for (int kth = 2; kth <= kthm1; ++kth) {
      theold[kth] = (theold[kth] - the0k) / sthek;
      if ((theold[kth] - theold[kth - 1]) < EPS) theold[kth] = theold[kth - 1] + EPS;
    }
    tb.the0[kp - 1] = the0k;
    tb.sthe[kp - 1] = sthek;
    thenew[1] = T(0.);
    thenew[kthm] = T(1.);
    const T dthe = T(1.) / T(kthm - 1);
    for (int kth = 2; kth <= kthm1; ++kth) thenew[kth] = thenew[kth - 1] + dthe;
    y2t[1] = T(0.);
    y2t[kthm] = T(0.);
    spline(kthm, theold, told, y2t, kthm, thenew, tnew, apt, aqt);
    for (int kth = 1; kth <= kthm; ++kth) tb.ttbl[kp - 1][kth - 1] = tnew[kth];
  }
  //---------------FINE LOOK-UP TABLE FOR T(P) FROM CONSTANT THE-----------
  kthm = JTBQ;
  kpm = ITBQ;
  kthm1 = kthm - 1;
  kpm1 = kpm - 1;
  dth = (C::THHQ - C::THL) / T(kthm - 1);
  dp = (C::PH - C::PLQ) / T(kpm - 1);
  th = C::THL - dth;
  p = C::PLQ - dp;
  for (int kp = 1; kp <= kpm; ++kp) {
    p = p + dp;
    th = C::THL - dth;
    for (int kth = 1; kth <= kthm; ++kth) {
      th = th + dth;
      ape = libm::pow(T(1.E5) / p, rd / cp);
      denom = th - C::A4 * ape;
      if (denom > EPS) {
        qs = C::PQ0 / p * libm::exp(C::A2 * (th - C::A3 * ape) / denom);
      } else {
        qs = T(0.);
      }
      toldq[kth] = th / ape;
      theoldq[kth] = th * libm::exp(elocp * qs / toldq[kth]);
    }
    const T the0k = theoldq[1];
    const T sthek = theoldq[kthm] - theoldq[1];
    theoldq[1] = T(0.);
    theoldq[kthm] = T(1.);
    for (int kth = 2; kth <= kthm1; ++kth) {
      theoldq[kth] = (theoldq[kth] - the0k) / sthek;
      if ((theoldq[kth] - theoldq[kth - 1]) < EPS) theoldq[kth] = theoldq[kth - 1] + EPS;
    }
    tb.the0q[kp - 1] = the0k;
    tb.stheq[kp - 1] = sthek;
    thenewq[1] = T(0.);
    thenewq[kthm] = T(1.);
    const T dthe = T(1.) / T(kthm - 1);
    for (int kth = 2; kth <= kthm1; ++kth) thenewq[kth] = thenewq[kth - 1] + dthe;
    y2tq[1] = T(0.);
    y2tq[kthm] = T(0.);
    spline(kthm, theoldq, toldq, y2tq, kthm, thenewq, tnewq, aptq, aqtq);
    for (int kth = 1; kth <= kthm; ++kth) tb.ttblq[kp - 1][kth - 1] = tnewq[kth];
  }
  (void)kpm1;
}

}  // namespace cppwrf::phys::cu_bmj
