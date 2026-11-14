// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Betts-Miller-Janjic cumulus scheme (cu_physics = 2, BMJSCHEME): WRF 4.8.0 phys/module_cu_bmj.F, BMJDRV for one
// column (bmjdrv_column), BMJ (bmj) and TTBLEX (ttblex); the look-up tables of BMJINIT are in bmj_tables.hpp.
// The arrays keep the Fortran 1-based level index: BMJ counts from the model top (L = 1) down to LMH, BMJDRV's
// JPR / PVPR use the model's own levels (K = 1 the lowest); they are sized kMaxLevels + 2. Every statement follows the
// Fortran in order and parenthesisation; x**2 is x*x, x**real is pow.
// gfortran -O2 -ftree-vectorize vectorises four loops with math functions (2-lane libmvec, cppwrf::vm with the
// iteration's index in the loop): APE = (1.E5/P)**CAPA (module_cu_bmj.F:521, pow), the deep convection's
// QREFK below LQM (:1148, exp), the shallow convection's QSATK (:1378, exp) and the Sterling branch of the cloud
// fraction (:311, log, exp; never reached: JPR(BBOT) = 1/PEPS); the other exp / pow / log / GAMMA calls are scalar.
#include <cmath>

#include "core/backend.hpp"
#include "core/column_work.hpp"
#include "core/libm.hpp"

namespace cppwrf::phys::cu_bmj {

inline constexpr int kMaxLevels = 100;
inline constexpr int KM = kMaxLevels + 2;

inline constexpr int ITB = 76, JTB = 134, ITBQ = 152, JTBQ = 440;
inline constexpr int ITREFI_MAX = 3;

// BMJINIT's tables (module variables): PTBL(IQ,IT) = ptbl[IT-1][IQ-1], TTBL(KTH,KP) = ttbl[KP-1][KTH-1],
// TTBLQ(KTH,KP) = ttblq[KP-1][KTH-1] (Fortran's column-major layout)
template <typename T>
struct BmjTables {
  T sthe[ITB], the0[ITB];
  T qs0[JTB], sqs[JTB];
  T stheq[ITBQ], the0q[ITBQ];
  T ptbl[JTB][ITB];
  T ttbl[ITB][JTB];
  T ttblq[ITBQ][JTBQ];
};

// MODULE_CU_BMJ's PARAMETERs and those of MODULE_MODEL_CONSTANTS it uses (evaluated in the working precision as
// gfortran folds them in the default REAL kind)
template <typename R>
struct BmjC {
  static constexpr R DTTOP = R(0.), EFIFC = R(5.0), EFIMN = R(0.20), EFMNT = R(0.70);
  static constexpr R ELIWV = R(2.683E6), ENPLO = R(20000.), ENPUP = R(15000.);
  static constexpr R EPSDN = R(1.05), EPSDT = R(0.);
  static constexpr R EPSNTP = R(.0001), EPSNTT = R(.0001), EPSPR = R(1.E-7);
  static constexpr R EPSUP = R(1.00);
  static constexpr R FR = R(1.00), FSL = R(0.85), FSS = R(0.85), GAM = R(0.5), PEPS = R(1.) / R(2.5);
  static constexpr R FUP = R(0.), FCC = R(5.00), CRMN = R(0.14), CRMX = R(85.0);
  static constexpr R PBM = R(13000.), PFRZ = R(15000.), PNO = R(1000.);
  static constexpr R PONE = R(2500.), PQM = R(20000.);
  static constexpr R PSH = R(20000.), PSHU = R(45000.);
  static constexpr R RHLSC = R(0.00), RHHSC = R(1.10);
  static constexpr R ROW = R(1.E3);
  static constexpr R STABDF = R(0.90), STABDS = R(0.90);
  static constexpr R STABS = R(1.0), STRESH = R(1.10);
  static constexpr R DTSHAL = R(-1.0), TREL = R(2400.);
  static constexpr R DTtrigr = -R(0.0), DTPtrigr = DTtrigr * PONE;
  static constexpr R DSPBFL = R(-3875.) * FR, DSP0FL = R(-5875.) * FR, DSPTFL = R(-1875.) * FR;
  static constexpr R DSPBFS = R(-3875.), DSP0FS = R(-5875.), DSPTFS = R(-1875.);
  static constexpr R PL = R(2500.), PLQ = R(70000.), PH = R(105000.);
  static constexpr R THL = R(210.), THH = R(365.), THHQ = R(325.);
  static constexpr R RDP = (R(ITB) - R(1.)) / (PH - PL), RDPQ = (R(ITBQ) - R(1.)) / (PH - PLQ);
  static constexpr R RDQ = R(ITB - 1), RDTH = (R(JTB) - R(1.)) / (THH - THL);
  static constexpr R RDTHE = R(JTB) - R(1.), RDTHEQ = R(JTBQ) - R(1.);
  static constexpr R RSFCP = R(1.) / R(101300.);
  static constexpr R AVGEFI = (EFIMN + R(1.)) * R(0.5);
  // BMJ's local PARAMETERs
  static constexpr R DSPBSL = DSPBFL * FSL, DSP0SL = DSP0FL * FSL, DSPTSL = DSPTFL * FSL;
  static constexpr R DSPBSS = DSPBFS * FSS, DSP0SS = DSP0FS * FSS, DSPTSS = DSPTFS * FSS;
  static constexpr R ELEVFC = R(0.6), STEFI = R(1.);
  static constexpr R SLOPBL = (DSPBFL - DSPBSL) / (R(1.) - EFIMN), SLOP0L = (DSP0FL - DSP0SL) / (R(1.) - EFIMN);
  static constexpr R SLOPTL = (DSPTFL - DSPTSL) / (R(1.) - EFIMN), SLOPBS = (DSPBFS - DSPBSS) / (R(1.) - EFIMN);
  static constexpr R SLOP0S = (DSP0FS - DSP0SS) / (R(1.) - EFIMN), SLOPTS = (DSPTFS - DSPTSS) / (R(1.) - EFIMN);
  static constexpr R SLOPST = (STABDF - STABDS) / (R(1.) - EFIMN), SLOPE = (R(1.) - EFMNT) / (R(1.) - EFIMN);
  // module_model_constants
  static constexpr R PQ0 = R(379.90516), A2 = R(17.2693882), A3 = R(273.16), A4 = R(35.86), EPSQ = R(1.e-12);
};

// TTBLEX: the temperature of the moist adiabat from TTBL (ITBX x JTBX, ttbl[(IPTB-1) JTBX + ITHTB-1])
template <typename R>
CPPWRF_HD inline R ttblex(int itbx, int jtbx, R plx, R prsmid, R rdpx, R rdthex, const R* sthe, const R* the0,
                          const R* ttbl, R thesp) {
  const R pk = prsmid;
  const R tpk = (pk - plx) * rdpx;
  R qq = tpk - std::trunc(tpk);
  int iptb = int(tpk) + 1;
  if (iptb < 1) {
    iptb = 1;
    qq = R(0.);
  }
  if (iptb >= itbx) {
    iptb = itbx - 1;
    qq = R(0.);
  }
  const R bthe00k = the0[iptb - 1], sthe00k = sthe[iptb - 1];
  const R bthe10k = the0[iptb], sthe10k = sthe[iptb];
  const R bthk = (bthe10k - bthe00k) * qq + bthe00k;
  const R sthk = (sthe10k - sthe00k) * qq + sthe00k;
  const R tthk = (thesp - bthk) / sthk * rdthex;
  R pp = tthk - std::trunc(tthk);
  int ithtb = int(tthk) + 1;
  if (ithtb < 1) {
    ithtb = 1;
    pp = R(0.);
  }
  if (ithtb >= jtbx) {
    ithtb = jtbx - 1;
    pp = R(0.);
  }
  const R t00k = ttbl[std::size_t(iptb - 1) * jtbx + (ithtb - 1)];
  const R t10k = ttbl[std::size_t(iptb - 1) * jtbx + ithtb];
  const R t01k = ttbl[std::size_t(iptb) * jtbx + (ithtb - 1)];
  const R t11k = ttbl[std::size_t(iptb) * jtbx + ithtb];
  return (t00k + (t10k - t00k) * pp + (t01k - t00k) * qq + (t00k - t10k - t01k + t11k) * pp * qq);
}

// the coarse (P < PLQ) or the fine table
template <typename R>
CPPWRF_HD inline R ttblex_p(const BmjTables<R>& tb, R p, R thes) {
  using C = BmjC<R>;
  if (p < C::PLQ) return ttblex<R>(ITB, JTB, C::PL, p, C::RDP, C::RDTHE, tb.sthe, tb.the0, &tb.ttbl[0][0], thes);
  return ttblex<R>(ITBQ, JTBQ, C::PLQ, p, C::RDPQ, C::RDTHEQ, tb.stheq, tb.the0q, &tb.ttblq[0][0], thes);
}

// BMJ's column arrays (core/column_work.hpp)
#define CPPWRF_BMJ_LOCALS(X)                                                                                       \
  X(R, apek, cu_bmj::KM) X(R, apesk, cu_bmj::KM) X(R, el, cu_bmj::KM) X(R, fpk, cu_bmj::KM) X(R, pk, cu_bmj::KM)   \
  X(R, psk, cu_bmj::KM) X(R, qk, cu_bmj::KM) X(R, qrefk, cu_bmj::KM) X(R, qsatk, cu_bmj::KM)                      \
  X(R, therk, cu_bmj::KM) X(R, thsk, cu_bmj::KM) X(R, thvref, cu_bmj::KM) X(R, tk, cu_bmj::KM)                     \
  X(R, trefk, cu_bmj::KM) X(R, ape, cu_bmj::KM) X(R, difq, cu_bmj::KM) X(R, dift, cu_bmj::KM)                     \
  X(R, thes, cu_bmj::KM) X(R, tref, cu_bmj::KM) X(R, cpe, cu_bmj::KM) X(R, cpecnv, cu_bmj::KM)                    \
  X(R, dtv, cu_bmj::KM) X(R, dtvcnv, cu_bmj::KM) X(R, thescnv, cu_bmj::KM)
// BMJDRV's column arrays
#define CPPWRF_BMJDRV_LOCALS(X)                                                                                    \
  X(R, dpcol, cu_bmj::KM) X(R, dqdt, cu_bmj::KM) X(R, dtdt, cu_bmj::KM) X(R, pcol, cu_bmj::KM)                     \
  X(R, qcol, cu_bmj::KM) X(R, tcol, cu_bmj::KM) X(R, pvpr, cu_bmj::KM) X(R, jpr, cu_bmj::KM)

template <typename R>
inline constexpr std::size_t kWorkBytes =
    0 CPPWRF_BMJ_LOCALS(CPPWRF_COLUMN_BYTES) CPPWRF_BMJDRV_LOCALS(CPPWRF_COLUMN_BYTES);

// SUBROUTINE BMJ (module_cu_bmj.F:391-1730) for KTS = 1 .. KTE; dprs, prsmid, q, t, dqdt, dtdt are 1-based columns
// from the top. A: the array type of the columns (R* on the host, ColArray<R> on the GPU)
template <typename R, class A>
CPPWRF_HD inline void bmj(const BmjTables<R>& tb, R dtcnvc, int lmh, R sm, R& cldefi, A dprs, A prsmid, A q, A t,
                          R psfc, A dqdt, A dtdt, R& pcpcol, int& lbot, int& ltop, R& pwcol, R& dqcol, R& dqcolmin,
                          R cp, R r, R elwv, R g, R tfrz, R d608, int kte, [[maybe_unused]] ColumnSlots cppwrf_cw_) {
  using C = BmjC<R>;
  constexpr int kts = 1;
  CPPWRF_BMJ_LOCALS(CPPWRF_COLUMN_ARRAY)
  auto QS0 = [&](int k) { return tb.qs0[k - 1]; };
  auto SQS = [&](int k) { return tb.sqs[k - 1]; };
  auto PTBL = [&](int iq, int it) { return tb.ptbl[it - 1][iq - 1]; };
  R apests, avrgt, avrgtl, bq, bqk, bqs00k, bqs10k, cup, den, dentpy, depth, depwl, dhdt, difql, diftl, dp, dpkl,
      dpmix, dpot, dqref, drheat, dsp, dsp0k, dspbk, dsptk, dsq, dst, dstq, dthem, dtdp, efi, fefi, fptk, hcorr, otsum,
      p, p00k, p01k, p10k, p11k, part1, part2, part3, pbot, pbotfc, pk0, pkb, pkl, pkt, pkxxxx, pkxxxy, plmh, pelevfc,
      pbtmx, plo, potsum, pp1, ppk, preck, prwd, pdqd, psp, psum, ptop, ptpk, pup, qbt, qkl, qnew, qotsum, qq1, qqk,
      qrfkl, qrftp, qsp, qsum, qup, rdp0t, rdpsum, rotsum, rtbar, rhavg, sm1, smix, sq, sqk, sqs00k, sqs10k, stabdl,
      sumde, sumdp, sumdt, tcorr, thbt, therkx, therky, thesp, thskl, thtpk, thvmkl, tkl, tnew, tq, tqk, trefkx, trfkl,
      trmlo, trmup, tsp, tth, tthk, tup, apekl, apekxx, apekxy, apes, capecnv, pspcnv, thbtcnv, capetrigr, cape, tlev2,
      qsat1, qsat2, rhshmax, qwat;
  int iq, iqtb, it, ittb, ittbk, l, l0, l0m1, lb, lbm1, lcor, lqm, ltp1, ltsh, lbotcnv, ltopcnv, lmid, kb;
  bool shallow;
  //---------------------------------------------------------------------
  const R capa = r / cp;
  const R cprlg = cp / (C::ROW * g * elwv);
  const R elocp = C::ELIWV / cp;
  const R rcp = R(1.) / cp;
  const R a23m4l = C::A2 * (C::A3 - C::A4) * elwv;
  const R rdtcnvc = R(1.) / dtcnvc;
  const R depmin = C::PSH * psfc * C::RSFCP;
  shallow = false;
  const R tauk = dtcnvc / C::TREL;
  const R tauksc = dtcnvc / (R(1.0) * C::TREL);
  //-----------------------------PREPARATIONS------------------------------
  lbot = lmh;
  dqcol = R(0.);
  pwcol = R(0.);
  pcpcol = R(0.);
  dqcolmin = R(0.);
  tref[kts] = t[kts];
  for (l = kts; l <= lmh; ++l) {
    apests = prsmid[l];
    ape[l] = vm::pow(R(1.E5) / apests, capa, l - kts, lmh - kts + 1);  // :521 vectorised
    cpecnv[l] = R(0.);
    dtvcnv[l] = R(0.);
    thescnv[l] = R(0.);
  }
  //----------------SEARCH FOR MAXIMUM BUOYANCY LEVEL----------------------
  plmh = prsmid[lmh];
  pelevfc = plmh * C::ELEVFC;
  pbtmx = prsmid[lmh] - C::PONE;
  capecnv = R(0.);
  pspcnv = R(0.);
  thbtcnv = R(0.);
  lbotcnv = lbot;
  ltopcnv = lbot;
  // values of the last trial level (WRF leaves them from max_buoy_loop when no level has CAPE)
  psp = R(0.); thbt = R(0.); pbot = R(0.); ptop = R(0.); ltop = lmh;
  for (kb = lmh; kb >= 1; --kb) {  // max_buoy_loop
    pkl = prsmid[kb];
    if (pkl < pelevfc) break;
    lbot = lmh;
    ltop = lmh;
    qbt = q[kb];
    thbt = t[kb] * ape[kb];
    tth = (thbt - C::THL) * C::RDTH;
    qq1 = tth - std::trunc(tth);
    ittb = int(tth) + 1;
    if (ittb < 1) {
      ittb = 1;
      qq1 = R(0.);
    } else if (ittb >= JTB) {
      ittb = JTB - 1;
      qq1 = R(0.);
    }
    ittbk = ittb;
    bqs00k = QS0(ittbk);
    sqs00k = SQS(ittbk);
    bqs10k = QS0(ittbk + 1);
    sqs10k = SQS(ittbk + 1);
    bq = (bqs10k - bqs00k) * qq1 + bqs00k;
    sq = (sqs10k - sqs00k) * qq1 + sqs00k;
    tq = (qbt - bq) / sq * C::RDQ;
    pp1 = tq - std::trunc(tq);
    iqtb = int(tq) + 1;
    if (iqtb < 1) {
      iqtb = 1;
      pp1 = R(0.);
    } else if (iqtb >= ITB) {
      iqtb = ITB - 1;
      pp1 = R(0.);
    }
    iq = iqtb;
    it = ittb;
    p00k = PTBL(iq, it);
    p10k = PTBL(iq + 1, it);
    p01k = PTBL(iq, it + 1);
    p11k = PTBL(iq + 1, it + 1);
    psp = p00k + (p10k - p00k) * pp1 + (p01k - p00k) * qq1 + (p00k - p10k - p01k + p11k) * pp1 * qq1;
    apes = libm::pow(R(1.E5) / psp, capa);
    thesp = thbt * libm::exp(elocp * qbt * apes / thbt);
    //-----------CHOOSE CLOUD BASE AS MODEL LEVEL JUST BELOW PSP-------------
    for (l = kts; l <= lmh - 1; ++l) {
      p = prsmid[l];
      if (p < psp && p >= C::PQM) lbot = l + 1;
    }
    pbot = prsmid[lbot];
    if (pbot >= pbtmx || lbot >= lmh) {
      for (l = kts; l <= lmh - 1; ++l) {
        p = prsmid[l];
        if (p < pbtmx) lbot = l;
      }
      pbot = prsmid[lbot];
    }
    //----------------CLOUD TOP COMPUTATION----------------------------------
    ltop = lbot;
    ptop = pbot;
    for (l = lmh; l >= kts; --l) thes[l] = thesp;
    //------------------FIRST ENTROPY CHECK----------------------------------
    for (l = kts; l <= lmh; ++l) {
      cpe[l] = R(0.);
      dtv[l] = R(0.);
    }
    dentpy = R(0.);
    l = kb;
    plo = prsmid[l];
    trmlo = R(0.);
    capetrigr = C::DTPtrigr / t[lbot];
    //--- Below cloud base
    if (kb > lbot) {
      for (l = kb - 1; l >= lbot + 1; --l) {
        pup = prsmid[l];
        tup = thbt / ape[l];
        dp = plo - pup;
        trmup = (tup * (qbt * R(0.608) + R(1.)) - t[l] * (q[l] * R(0.608) + R(1.))) * R(0.5) /
                (t[l] * (q[l] * R(0.608) + R(1.)));
        dtv[l] = trmlo + trmup;
        dentpy = dtv[l] * dp + dentpy;
        cpe[l] = dentpy;
        if (dentpy < capetrigr) goto l170;
        plo = pup;
        trmlo = trmup;
      }
    } else {
      l = lbot + 1;
      plo = prsmid[l];
      tup = thbt / ape[l];
      trmlo = (tup * (qbt * R(0.608) + R(1.)) - t[l] * (q[l] * R(0.608) + R(1.))) * R(0.5) /
              (t[l] * (q[l] * R(0.608) + R(1.)));
    }
    //--- At cloud base
    l = lbot;
    pup = psp;
    tup = thbt / apes;
    tsp = (t[l + 1] - t[l]) / (plo - pbot) * (pup - pbot) + t[l];
    qsp = (q[l + 1] - q[l]) / (plo - pbot) * (pup - pbot) + q[l];
    dp = plo - pup;
    trmup = (tup * (qbt * R(0.608) + R(1.)) - tsp * (qsp * R(0.608) + R(1.))) * R(0.5) / (tsp * (qsp * R(0.608) + R(1.)));
    dtv[l] = trmlo + trmup;
    dentpy = dtv[l] * dp + dentpy;
    cpe[l] = dentpy;
    dtv[l] = dtv[l] * dp;
    plo = pup;
    trmlo = trmup;
    pup = prsmid[l];
    //--- Calculate updraft temperature along moist adiabat (TUP)
    tup = ttblex_p(tb, pup, R(thes[l]));
    qup = C::PQ0 / pup * libm::exp(C::A2 * (tup - C::A3) / (tup - C::A4));
    qwat = qbt - qup;
    dp = plo - pup;
    trmup = (tup * (qup * R(0.608) + R(1.) - qwat) - t[l] * (q[l] * R(0.608) + R(1.))) * R(0.5) /
            (t[l] * (q[l] * R(0.608) + R(1.)));
    dentpy = (trmlo + trmup) * dp + dentpy;
    cpe[l] = dentpy;
    dtv[l] = (dtv[l] + (trmlo + trmup) * dp) / (prsmid[lbot + 1] - prsmid[lbot]);
    if (dentpy < capetrigr) goto l170;
    plo = pup;
    trmlo = trmup;
    //--- In cloud above cloud base
    for (l = lbot - 1; l >= kts; --l) {
      pup = prsmid[l];
      tup = ttblex_p(tb, pup, R(thes[l]));
      qup = C::PQ0 / pup * libm::exp(C::A2 * (tup - C::A3) / (tup - C::A4));
      qwat = qbt - qup;
      dp = plo - pup;
      trmup = (tup * (qup * R(0.608) + R(1.) - qwat) - t[l] * (q[l] * R(0.608) + R(1.))) * R(0.5) /
              (t[l] * (q[l] * R(0.608) + R(1.)));
      dtv[l] = trmlo + trmup;
      dentpy = dtv[l] * dp + dentpy;
      cpe[l] = dentpy;
      if (dentpy < capetrigr) goto l170;
      plo = pup;
      trmlo = trmup;
    }
  l170:
    ltp1 = kb;
    cape = R(0.);
    //--- Cloud top level (LTOP) is located where CAPE is a maximum
    for (l = kb; l >= kts; --l) {
      if (cpe[l] < capetrigr) {
        break;
      } else if (cpe[l] > cape) {
        ltp1 = l;
        cape = cpe[l];
      }
    }
    ltop = ltp1 < lbot ? ltp1 : lbot;
    //--------------- CHECK FOR MAXIMUM INSTABILITY  ------------------------
    if (cape > capecnv) {
      capecnv = cape;
      pspcnv = psp;
      thbtcnv = thbt;
      lbotcnv = lbot;
      ltopcnv = ltop;
      for (l = lmh; l >= kts; --l) {
        cpecnv[l] = cpe[l];
        dtvcnv[l] = dtv[l];
        thescnv[l] = thes[l];
      }
    }
  }  // max_buoy_loop
  //------------------------  MAXIMUM INSTABILITY  ------------------------
  if (capecnv > R(0.)) {
    psp = pspcnv;
    thbt = thbtcnv;
    lbot = lbotcnv;
    ltop = ltopcnv;
    pbot = prsmid[lbot];
    ptop = prsmid[ltop];
    for (l = lmh; l >= kts; --l) {
      cpe[l] = cpecnv[l];
      dtv[l] = dtvcnv[l];
      thes[l] = thescnv[l];
    }
  }
  //-----  Quick exit if cloud is too thin or no CAPE is present  ---------
  if (ptop > pbot - C::PNO || ltop > lbot - 2 || capecnv <= R(0.)) {
    lbot = 0;
    ltop = kte;
    pbot = prsmid[lmh];
    ptop = pbot;
    cldefi = C::AVGEFI * sm + C::STEFI * (R(1.) - sm);
    return;  // GO TO 800
  }
  //***  DEPTH OF CLOUD REQUIRED TO MAKE THE POINT A DEEP CONVECTION POINT
  depth = pbot - ptop;
  if (!(depth >= depmin)) {
    shallow = true;
    goto l600;
  }
  //DCDCDC    DEEP CONVECTION   DCDCDC
  lb = lbot;
  efi = cldefi;
  for (l = kts; l <= lmh; ++l) {
    dift[l] = R(0.);
    difq[l] = R(0.);
    tkl = t[l];
    tk[l] = tkl;
    trefk[l] = tkl;
    qkl = q[l];
    qk[l] = qkl;
    qrefk[l] = qkl;
    pkl = prsmid[l];
    pk[l] = pkl;
    psk[l] = pkl;
    apekl = ape[l];
    apek[l] = apekl;
    //--- Calculate temperature along moist adiabat (TREF)
    tref[l] = ttblex_p(tb, pkl, R(thes[l]));
    therk[l] = tref[l] * apekl;
  }
  //------------DEEP CONVECTION REFERENCE TEMPERATURE PROFILE------------
  ltp1 = ltop + 1;
  lbm1 = lb - 1;
  pkb = pk[lb];
  pkt = pk[ltop];
  stabdl = (efi - C::EFIMN) * C::SLOPST + C::STABDS;
  //------------TEMPERATURE REFERENCE PROFILE BELOW FREEZING LEVEL-------
  el[lb] = elwv;
  l0 = lb;
  pk0 = pk[lb];
  trefkx = trefk[lb];
  therkx = therk[lb];
  apekxx = apek[lb];
  therky = therk[lbm1];
  apekxy = apek[lbm1];
  for (l = lbm1; l >= ltop; --l) {
    if (t[l + 1] < tfrz) goto l430;
    trefkx = ((therky - therkx) * stabdl + trefkx * apekxx) / apekxy;
    trefk[l] = trefkx;
    el[l] = elwv;
    apekxx = apekxy;
    therkx = therky;
    apekxy = apek[l - 1];
    therky = therk[l - 1];
    l0 = l;
    pk0 = pk[l0];
  }
  //--------------FREEZING LEVEL AT OR ABOVE THE CLOUD TOP-----------------
  goto l450;
  //--------------TEMPERATURE REFERENCE PROFILE ABOVE FREEZING LEVEL-------
l430:
  l0m1 = l0 - 1;
  rdp0t = R(1.) / (pk0 - pkt);
  dthem = therk[l0] - trefk[l0] * apek[l0];
  for (l = ltop; l <= l0m1; ++l) {
    trefk[l] = (therk[l] - (pk[l] - pkt) * dthem * rdp0t) / apek[l];
    el[l] = elwv;  // ELIV
  }
  //--------------DEEP CONVECTION REFERENCE HUMIDITY PROFILE---------------
l450:
  depwl = pkb - pk0;
  depth = C::PFRZ * psfc * C::RSFCP;
  sm1 = R(1.) - sm;
  pbotfc = R(1.);
  //--------------- ITERATION LOOP FOR CLOUD EFFICIENCY -------------------
  sumdp = R(0.); dentpy = R(0.); preck = R(0.); pdqd = R(0.); prwd = R(0.);
  for (int itrefi = 1; itrefi <= ITREFI_MAX; ++itrefi) {  // cloud_efficiency
    dspbk = ((efi - C::EFIMN) * C::SLOPBS + C::DSPBSS * pbotfc) * sm + ((efi - C::EFIMN) * C::SLOPBL + C::DSPBSL * pbotfc) * sm1;
    dsp0k = ((efi - C::EFIMN) * C::SLOP0S + C::DSP0SS * pbotfc) * sm + ((efi - C::EFIMN) * C::SLOP0L + C::DSP0SL * pbotfc) * sm1;
    dsptk = ((efi - C::EFIMN) * C::SLOPTS + C::DSPTSS * pbotfc) * sm + ((efi - C::EFIMN) * C::SLOPTL + C::DSPTSL * pbotfc) * sm1;
    for (l = ltop; l <= lb; ++l) {
      //***  SATURATION PRESSURE DIFFERENCE
      if (depwl >= depth) {
        if (l < l0) {
          dsp = ((pk0 - pk[l]) * dsptk + (pk[l] - pkt) * dsp0k) / (pk0 - pkt);
        } else {
          dsp = ((pkb - pk[l]) * dsp0k + (pk[l] - pk0) * dspbk) / (pkb - pk0);
        }
      } else {
        dsp = dsp0k;
        if (l < l0) dsp = ((pk0 - pk[l]) * dsptk + (pk[l] - pkt) * dsp0k) / (pk0 - pkt);
      }
      //***  HUMIDITY PROFILE
      psk[l] = pk[l] + dsp;
      apesk[l] = libm::pow(R(1.E5) / psk[l], capa);
      if (pk[l] > C::PQM) {
        thsk[l] = trefk[l] * apek[l];
        qrefk[l] = C::PQ0 / psk[l] * libm::exp(C::A2 * (thsk[l] - C::A3 * apesk[l]) / (thsk[l] - C::A4 * apesk[l]));
      } else {
        qrefk[l] = qk[l];
      }
    }
    //***  ENTHALPY CONSERVATION INTEGRAL
    for (int iter = 1; iter <= 2; ++iter) {  // enthalpy_conservation
      sumde = R(0.);
      sumdp = R(0.);
      dhdt = R(0.);
      for (l = ltop; l <= lb; ++l) {
        sumde = ((tk[l] - trefk[l]) * cp + (qk[l] - qrefk[l]) * el[l]) * dprs[l] + sumde;
        const R x = (trefk[l] * apek[l] / apesk[l]) - C::A4;
        dhdt = (qrefk[l] * a23m4l / (x * x) + cp) * dprs[l] + dhdt;
        sumdp = sumdp + dprs[l];
      }
      hcorr = sumde / (sumdp - dprs[ltop]);
      dhdt = dhdt / (sumdp - dprs[ltop]);
      lcor = ltop + 1;
      //***  FIND LQM
      lqm = 1;
      for (l = kts; l <= lb; ++l) {
        if (pk[l] <= C::PQM) lqm = l;
      }
      //***  ABOVE LQM CORRECT TEMPERATURE ONLY
      if (lcor <= lqm) {
        for (l = lcor; l <= lqm; ++l) trefk[l] = trefk[l] + hcorr * rcp;
        lcor = lqm + 1;
      }
      //***  BELOW LQM CORRECT BOTH TEMPERATURE AND MOISTURE
      for (l = lcor; l <= lb; ++l) {
        trefk[l] = hcorr / dhdt + trefk[l];
        thskl = trefk[l] * apek[l];
        qrefk[l] = C::PQ0 / psk[l] *
                   vm::exp(C::A2 * (thskl - C::A3 * apesk[l]) / (thskl - C::A4 * apesk[l]), l - lcor, lb - lcor + 1);
      }
    }
    //***  HEATING, MOISTENING, PRECIPITATION
    avrgt = R(0.);
    preck = R(0.);
    pdqd = R(0.);
    prwd = R(0.);
    dsq = R(0.);
    dst = R(0.);
    for (l = ltop; l <= lb; ++l) {
      tkl = tk[l];
      diftl = (trefk[l] - tkl) * tauk;
      difql = (qrefk[l] - qk[l]) * tauk;
      avrgtl = (tkl + tkl + diftl);
      dpot = dprs[l] / avrgtl;
      dst = diftl * dpot + dst;
      dsq = difql * el[l] * dpot + dsq;
      avrgt = avrgtl * dprs[l] + avrgt;
      preck = diftl * dprs[l] + preck;
      pdqd = (qk[l] - qrefk[l]) * dprs[l] + pdqd;
      prwd = qk[l] * dprs[l] + prwd;
      dift[l] = diftl;
      difq[l] = difql;
    }
    dst = (dst + dst) * cp;
    dsq = dsq + dsq;
    dentpy = dst + dsq;
    avrgt = avrgt / (sumdp + sumdp);
    drheat = (preck * sm + c_max(R(1.E-7), preck) * (R(1.) - sm)) * cp / avrgt;  // As in Eta!
    drheat = c_max(drheat, R(1.E-20));
    efi = C::EFIFC * dentpy / drheat;
    efi = c_min(efi, R(1.));
    efi = c_max(efi, C::EFIMN);
  }
  //---------------------- DEEP CONVECTION --------------------------------
  if (dentpy >= C::EPSNTP && preck > C::EPSPR) {
    cldefi = efi;
    fefi = C::EFMNT + C::SLOPE * (efi - C::EFIMN);
    fefi = (dentpy - C::EPSNTP) * fefi / dentpy;
    preck = preck * fefi;
    //***  UPDATE PRECIPITATION AND TENDENCIES OF TEMPERATURE AND MOISTURE
    cup = preck * cprlg;
    pcpcol = cup;
    dqcol = pdqd / g;
    pwcol = prwd / g;
    dqcolmin = (C::CRMN * C::TREL) / (fefi * R(86400.));
    for (l = ltop; l <= lb; ++l) {
      dtdt[l] = dift[l] * fefi * rdtcnvc;
      dqdt[l] = difq[l] * fefi * rdtcnvc;
    }
  } else {
    cldefi = C::EFIMN * sm + C::STEFI * (R(1.) - sm);
    //***  SEARCH FOR SHALLOW CLOUD TOP
    ptpk = c_max(C::PSHU, pk[lbot] - depmin);
    //***  CLOUD TOP IS THE LEVEL JUST BELOW PBTK-PSH or JUST BELOW PSHU
    for (l = kts; l <= lmh; ++l) {
      if (pk[l] <= ptpk) ltop = l + 1;
    }
    //-- Make shallow cloud top a function of virtual temperature excess (DTV)
    ltp1 = lbot;
    for (l = lbot - 1; l >= ltop; --l) {
      if (dtv[l] > R(0.)) {
        ltp1 = l;
      } else {
        break;
      }
    }
    ltop = ltp1 < lbot ? ltp1 : lbot;
    ptop = pk[ltop];
    shallow = true;
  }
  //DCDCDC          END OF DEEP CONVECTION            DCDCDC
l600:
  if (!shallow) return;  // GO TO 800
  //SCSCSC         SHALLOW CONVECTION          SCSCSC
  for (l = kts; l <= lmh; ++l) {
    tkl = t[l];
    tk[l] = tkl;
    trefk[l] = tkl;
    qkl = q[l];
    qk[l] = qkl;
    qrefk[l] = qkl;
    pkl = prsmid[l];
    pk[l] = pkl;
    qsatk[l] = C::PQ0 / pk[l] * vm::exp(C::A2 * (tk[l] - C::A3) / (tk[l] - C::A4), l - kts, lmh - kts + 1);  // :1378
    apekl = ape[l];
    apek[l] = apekl;
    thvmkl = tkl * apekl * (qkl * d608 + R(1.));
    thvref[l] = thvmkl;
    el[l] = elwv;
  }
  //-- Begin: Raise cloud top if avg RH>RHSHmax and CAPE>0
  tlev2 = t[lbot] * libm::pow((pk[lbot] - C::PONE) / pk[lbot], capa);
  qsat1 = C::PQ0 / pk[lbot] * libm::exp(C::A2 * (t[lbot] - C::A3) / (tk[lbot] - C::A4));
  qsat2 = C::PQ0 / (pk[lbot] - C::PONE) * libm::exp(C::A2 * (tlev2 - C::A3) / (tlev2 - C::A4));
  rhshmax = qsat2 / qsat1;
  sumdp = R(0.);
  rhavg = R(0.);
  for (l = lbot; l >= ltop; --l) {
    rhavg = rhavg + dprs[l] * qk[l] / qsatk[l];
    sumdp = sumdp + dprs[l];
  }
  if (rhavg / sumdp > rhshmax) {
    ltsh = ltop;
    for (l = ltop - 1; l >= kts; --l) {
      rhavg = rhavg + dprs[l] * qk[l] / qsatk[l];
      sumdp = sumdp + dprs[l];
      if (cpe[l] > R(0.)) {
        ltsh = l;
      } else {
        break;
      }
      if (rhavg / sumdp <= rhshmax) break;
      if (pk[l] <= C::PSHU) break;
    }
    ltop = ltsh;
  }
  //---------------------------SHALLOW CLOUD TOP---------------------------
  lbm1 = lbot - 1;
  ptpk = ptop;
  ltp1 = ltop - 1;
  depth = pbot - ptop;
  if (ptop > pbot - C::PNO || ltop > lbot - 2) {
    lbot = 0;
    ltop = kte;
    ptop = pbot;
    return;
  }
  //--------------SCALING POTENTIAL TEMPERATURE & TABLE INDEX AT TOP-------
  thtpk = t[ltp1] * ape[ltp1];
  tthk = (thtpk - C::THL) * C::RDTH;
  qqk = tthk - std::trunc(tthk);
  it = int(tthk) + 1;
  if (it < 1) {
    it = 1;
    qqk = R(0.);
  }
  if (it >= JTB) {
    it = JTB - 1;
    qqk = R(0.);
  }
  //--------------BASE AND SCALING FACTOR FOR SPEC. HUMIDITY AT TOP--------
  bqs00k = QS0(it);
  sqs00k = SQS(it);
  bqs10k = QS0(it + 1);
  sqs10k = SQS(it + 1);
  //--------------SCALING SPEC. HUMIDITY & TABLE INDEX AT TOP--------------
  bqk = (bqs10k - bqs00k) * qqk + bqs00k;
  sqk = (sqs10k - sqs00k) * qqk + sqs00k;
  tqk = (q[ltp1] - bqk) / sqk * C::RDQ;
  ppk = tqk - std::trunc(tqk);
  iq = int(tqk) + 1;
  if (iq < 1) {
    iq = 1;
    ppk = R(0.);
  }
  if (iq >= ITB) {
    iq = ITB - 1;
    ppk = R(0.);
  }
  //----------------CLOUD TOP SATURATION POINT PRESSURE--------------------
  part1 = (PTBL(iq + 1, it) - PTBL(iq, it)) * ppk;
  part2 = (PTBL(iq, it + 1) - PTBL(iq, it)) * qqk;
  part3 = (PTBL(iq, it) - PTBL(iq + 1, it) - PTBL(iq, it + 1) + PTBL(iq + 1, it + 1)) * ppk * qqk;
  ptpk = PTBL(iq, it) + part1 + part2 + part3;
  dpmix = ptpk - psp;
  if (std::fabs(dpmix) < R(3000.)) dpmix = R(-3000.);
  //----------------TEMPERATURE PROFILE SLOPE------------------------------
  smix = (thtpk - thbt) / dpmix * C::STABS;
  trefkx = trefk[lbot + 1];
  pkxxxx = pk[lbot + 1];
  pkxxxy = pk[lbot];
  apekxx = apek[lbot + 1];
  apekxy = apek[lbot];
  lmid = int(R(.5) * R(lbot + ltop));
  for (l = lbot; l >= ltop; --l) {
    trefkx = ((pkxxxy - pkxxxx) * smix + trefkx * apekxx) / apekxy;
    trefk[l] = trefkx;
    if (l <= lmid) trefk[l] = c_max(R(trefk[l]), R(tk[l] + C::DTSHAL));
    apekxx = apekxy;
    pkxxxx = pkxxxy;
    apekxy = apek[l - 1];
    pkxxxy = pk[l - 1];
  }
  //----------------TEMPERATURE REFERENCE PROFILE CORRECTION---------------
  sumdt = R(0.);
  sumdp = R(0.);
  for (l = ltop; l <= lbot; ++l) {
    sumdt = (tk[l] - trefk[l]) * dprs[l] + sumdt;
    sumdp = sumdp + dprs[l];
  }
  rdpsum = R(1.) / sumdp;
  fpk[lbot] = trefk[lbot];
  tcorr = sumdt * rdpsum;
  for (l = ltop; l <= lbot; ++l) {
    trfkl = trefk[l] + tcorr;
    trefk[l] = trfkl;
    fpk[l] = trfkl;
  }
  //----------------HUMIDITY PROFILE EQUATIONS-----------------------------
  psum = R(0.);
  qsum = R(0.);
  potsum = R(0.);
  qotsum = R(0.);
  otsum = R(0.);
  dst = R(0.);
  fptk = fpk[ltop];
  for (l = ltop; l <= lbot; ++l) {
    dpkl = fpk[l] - fptk;
    psum = dpkl * dprs[l] + psum;
    qsum = qk[l] * dprs[l] + qsum;
    rtbar = R(2.) / (trefk[l] + tk[l]);
    otsum = dprs[l] * rtbar + otsum;
    potsum = dpkl * rtbar * dprs[l] + potsum;
    qotsum = qk[l] * rtbar * dprs[l] + qotsum;
    dst = (trefk[l] - tk[l]) * rtbar * dprs[l] / el[l] + dst;
  }
  psum = psum * rdpsum;
  qsum = qsum * rdpsum;
  rotsum = R(1.) / otsum;
  potsum = potsum * rotsum;
  qotsum = qotsum * rotsum;
  dst = dst * rotsum * cp;
  //----------------ENSURE POSITIVE ENTROPY CHANGE-------------------------
  if (dst > R(0.)) {
    lbot = 0;
    ltop = kte;
    ptop = pbot;
    return;
  } else {
    dstq = dst * C::EPSDN;
  }
  //----------------CHECK FOR ISOTHERMAL ATMOSPHERE------------------------
  den = potsum - psum;
  if (-den / psum < R(5.E-5)) {
    lbot = 0;
    ltop = kte;
    ptop = pbot;
    return;
  } else {
    dqref = (qotsum - dstq - qsum) / den;
  }
  //-------------- HUMIDITY DOES NOT INCREASE WITH HEIGHT------------------
  if (dqref < R(0.)) {
    lbot = 0;
    ltop = kte;
    ptop = pbot;
    return;
  }
  //----------------HUMIDITY AT THE CLOUD TOP------------------------------
  qrftp = qsum - dqref * psum;
  //----------------HUMIDITY PROFILE---------------------------------------
  for (l = ltop; l <= lbot; ++l) {
    qrfkl = (fpk[l] - fptk) * dqref + qrftp;
    //***  TOO DRY CLOUDS NOT ALLOWED
    tnew = (trefk[l] - tk[l]) * tauksc + tk[l];
    qsatk[l] = C::PQ0 / pk[l] * libm::exp(C::A2 * (tnew - C::A3) / (tnew - C::A4));
    qnew = (qrfkl - qk[l]) * tauksc + qk[l];
    if (qnew < qsatk[l] * C::RHLSC) {
      lbot = 0;
      ltop = kte;
      ptop = pbot;
      return;
    }
    //-------------TOO MOIST CLOUDS NOT ALLOWED------------------------------
    if (qnew > qsatk[l] * C::RHHSC) {
      lbot = 0;
      ltop = kte;
      ptop = pbot;
      return;
    }
    thvref[l] = trefk[l] * apek[l] * (qrfkl * d608 + R(1.));
    qrefk[l] = qrfkl;
  }
  //-------------- ELIMINATE IMPOSSIBLE SLOPES (BETTS,DTHETA/DQ)------------
  for (l = ltop; l <= lbot; ++l) {
    dtdp = (thvref[l - 1] - thvref[l]) / (prsmid[l] - prsmid[l - 1]);
    if (dtdp < C::EPSDT) {
      lbot = 0;
      ltop = kte;
      ptop = pbot;
      return;
    }
  }
  //--------------RELAXATION TOWARD REFERENCE PROFILES---------------------
  for (l = ltop; l <= lbot; ++l) {
    dtdt[l] = (trefk[l] - tk[l]) * tauksc * rdtcnvc;
    dqdt[l] = (qrefk[l] - qk[l]) * tauksc * rdtcnvc;
  }
  (void)lbm1; (void)ltp1; (void)depth; (void)pkxxxx;
}

// SUBROUTINE BMJDRV (module_cu_bmj.F:76-387) for column (i, j). F: the column's fields (a functor with the views,
// C++ k = WRF k - 1), see KBmj in phys/driver/cumulus_driver.cpp.
template <typename R, class F>
CPPWRF_HD inline void bmjdrv_column(const F& f, const BmjTables<R>& tb, int i, int j, int kte, R dt, int stepcu,
                                    R cp, R r, R elwv, R g, R tfrz, R d608, bool bmj_rad_feedback,
                                    [[maybe_unused]] ColumnSlots cppwrf_cw_) {
  using C = BmjC<R>;
  constexpr int kts = 1;
  CPPWRF_BMJDRV_LOCALS(CPPWRF_COLUMN_ARRAY)
  const R dtcnvc = dt * R(stepcu);
  for (int k = kts; k <= kte; ++k) {
    dqdt[k] = R(0.);
    dtdt[k] = R(0.);
    jpr[k] = R(0.);
    pvpr[k] = R(0.);
    f.qcconv(i, j, k - 1) = R(0.);
    f.qiconv(i, j, k - 1) = R(0.);
    f.ccldfra(i, j, k - 1) = R(0.);
  }
  R dqcol = R(0.), pwcol = R(0.), pcpcol = R(0.), dqcolmin = R(0.);
  f.raincv(i, j) = R(0.);
  f.pratec(i, j) = R(0.);
  f.convcld(i, j) = R(0.);
  const int lowlyr = f.lowlyr(i, j);
  const R psfc = f.p8w(i, j, lowlyr - 1);
  //***  CONVERT TO BMJ LAND MASK (1.0 FOR SEA; 0.0 FOR LAND)
  const R landmask = f.xland(i, j) - R(1.);
  //***  FILL 1-D VERTICAL ARRAYS AND FLIP DIRECTION SINCE BMJ SCHEME COUNTS DOWNWARD FROM THE DOMAIN'S TOP
  for (int k = kts; k <= kte; ++k) {
    const int kflip = kte + 1 - k;
    const R qv = f.qv(i, j, kflip - 1);
    qcol[k] = c_max(C::EPSQ, qv / (R(1.) + qv));
    tcol[k] = f.t(i, j, kflip - 1);
    pcol[k] = f.p(i, j, kflip - 1);
    dpcol[k] = f.rho(i, j, kflip - 1) * g * f.dz8w(i, j, kflip - 1);
  }
  const int lmh = kte + 1 - lowlyr;
  int lbot = 0, ltop = 0;
  R cldefi = f.cldefi(i, j);
  bmj<R>(tb, dtcnvc, lmh, landmask, cldefi, dpcol, pcol, qcol, tcol, psfc, dqdt, dtdt, pcpcol, lbot, ltop, pwcol, dqcol,
         dqcolmin, cp, r, elwv, g, tfrz, d608, kte, cppwrf_cw_);
  f.cldefi(i, j) = cldefi;
  //***  COMPUTE HEATING AND MOISTENING TENDENCIES
  for (int k = kts; k <= kte; ++k) {
    const int kflip = kte + 1 - k;
    f.rthcuten(i, j, k - 1) = dtdt[kflip] / f.pi(i, j, k - 1);
    //***  CONVERT FROM SPECIFIC HUMIDTY BACK TO MIXING RATIO
    const R d = R(1.) - qcol[kflip];
    f.rqvcuten(i, j, k - 1) = dqdt[kflip] / (d * d);
  }
  //***  PRECIP FROM METERS TO MILLIMETERS PER STEP FOR OUTPUT
  f.raincv(i, j) = pcpcol * R(1.E3) / R(stepcu);
  const R pratec = pcpcol * R(1.E3) / (R(stepcu) * dt);
  f.pratec(i, j) = pratec;
  //***  CONVECTIVE CLOUD TOP AND BOTTOM FROM THIS CALL
  const R cutop = R(kte + 1 - ltop), cubot = R(kte + 1 - lbot);
  f.cutop(i, j) = cutop;
  f.cubot(i, j) = cubot;
  if (bmj_rad_feedback) {
    if (dqcol > dqcolmin) {
      //***  CONVECTIVE CLOUD FRACTION: BASED ON SLINGO (1987) WITH A POISSON VERTICAL PROFILE
      R cumx = R(0.), mpvpr, factl;
      const R prrt = (pratec * R(86400.0)) / C::CRMN;
      const R rrp = R(0.8) / (libm::log(C::CRMX / C::CRMN));  // folded by gfortran; glibc's log gives the same
      if (prrt < C::CRMX / C::CRMN) {
        cumx = rrp * libm::log(prrt);
      } else {
        cumx = R(0.8);
      }
      //***  THE CONVECTIVE CLOUD FRACTION (CCLDFRA) AT EACH MODEL LEVEL; STERLING FOR N >= 17
      const int ttop = int(std::lround(cutop));
      const int bbot = int(std::lround(cubot));
      const R pavg = R(1.) / ((C::PEPS * R(3.)) * (C::PEPS * R(3.)));
      for (int k = kts; k <= kte; ++k) {
        if (k >= bbot && k <= ttop) {
          jpr[k] = (R(1.) / C::PEPS) * ((f.p(i, j, k - 1) - f.p(i, j, ttop - 1)) /
                                        (f.p(i, j, bbot - 1) - f.p(i, j, ttop - 1)));
        } else {
          jpr[k] = R(0.0);
        }
      }
      for (int k = kts; k <= kte; ++k) pvpr[k] = R(0.);
      if (jpr[bbot] < R(17)) {
        for (int k = bbot; k <= ttop; ++k) pvpr[k] = libm::pow(pavg, R(jpr[k])) / libm::tgamma(R(jpr[k] + R(1.)));
      } else {
        const R lpavg = libm::log(R(1.) * pavg);  // folded by gfortran; glibc's log gives the same
        const R pi = R(3.141592653589793238462643383279502884);  // ACOS(-1.)
        const int n = ttop - bbot + 1;
        for (int k = bbot; k <= ttop; ++k) {  // :311 vectorised
          const R x = jpr[k];
          factl = x * vm::log(x, k - bbot, n) - x + R(0.5) * vm::log(R(2.) * x * pi, k - bbot, n) +
                  vm::log(R(1.) + R(1.) / (R(12.) * x) + R(1.) / (R(288.) * (x * x)), k - bbot, n);
          pvpr[k] = vm::exp((x)*lpavg - factl, k - bbot, n);
        }
      }
      mpvpr = pvpr[kts];  // MAXVAL(PVPR)
      for (int k = kts + 1; k <= kte; ++k)
        if (pvpr[k] > mpvpr) mpvpr = pvpr[k];
      for (int k = bbot; k <= ttop; ++k) pvpr[k] = pvpr[k] / mpvpr;
      for (int k = kts; k <= kte; ++k) f.ccldfra(i, j, k - 1) = cumx * pvpr[k];
      //***  THE CONVECTIVE CLOUD CONDENSATES (QCCONV,QICONV)
      R mcol = R(0.);
      for (int k = int(cubot); k <= int(cutop); ++k) {
        const int kflip = kte + 1 - k;
        mcol = f.rho(i, j, k - 1) * f.dz8w(i, j, k - 1) * qcol[kflip] * f.ccldfra(i, j, k - 1) + mcol;
      }
      const R convcld = libm::pow(pwcol, C::GAM) * libm::pow(dqcol - dqcolmin, R(1.) - C::GAM);
      f.convcld(i, j) = convcld;
      const R qcis = convcld / mcol;
      for (int k = kts; k <= kte; ++k) {
        const int kflip = kte + 1 - k;
        if (tcol[kflip] >= tfrz) {
          f.qiconv(i, j, k - 1) = R(0.);
          f.qcconv(i, j, k - 1) = (qcis * qcol[kflip] * f.ccldfra(i, j, k - 1)) / (R(1.) - qcol[kflip]);
        } else {
          f.qiconv(i, j, k - 1) = (qcis * qcol[kflip] * f.ccldfra(i, j, k - 1)) / (R(1.) - qcol[kflip]);
          f.qcconv(i, j, k - 1) = R(0.);
        }
      }
      (void)factl;
    }
  }
}

}  // namespace cppwrf::phys::cu_bmj
