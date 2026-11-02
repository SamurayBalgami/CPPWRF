#include "phys/ra_rrtm_lw/rrtm_lw.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "core/backend.hpp"

namespace cppwrf::phys {

namespace {

constexpr int NGPT = 140, NBANDS = 16, MAXXSEC = 4, NMOL = 6;
constexpr int NG1 = 8, NG2 = 14, NG3 = 16, NG4 = 14, NG5 = 16, NG6 = 8, NG7 = 12, NG8 = 8, NG9 = 12, NG10 = 6,
              NG11 = 8, NG12 = 8, NG13 = 4, NG14 = 2, NG15 = 2, NG16 = 2;
constexpr double ONEMINUS = 1. - 1.E-6;
constexpr double deltap = 4.;  // pressure interval of the buffer layers (mb)

// A(i,j) of a Fortran array with first extent n1 (1-based indices)
inline double F2(const double* a, int n1, int i, int j) { return a[(j - 1) * n1 + (i - 1)]; }

struct Col {  // the shared arrays of RRTM's subroutines, Fortran-indexed
  RrtmLwWorkspace& w;
  int nl;
  double& TAUG(int ig, int lay) { return w.taug[std::size_t((lay - 1) * NGPT + ig - 1)]; }
  double& PFRAC(int ig, int lay) { return w.pfrac[std::size_t((lay - 1) * NGPT + ig - 1)]; }
  int& ITR(int ig, int lay) { return w.itr[std::size_t((lay - 1) * NGPT + ig - 1)]; }
  double& WKL(int m, int lay) { return w.wkl[std::size_t((lay - 1) * 35 + m - 1)]; }
  double& WX(int m, int lay) { return w.wx[std::size_t((lay - 1) * MAXXSEC + m - 1)]; }
};

// O3DATA (after INIRAD's O3PROF = 0): the annual mean of the summer / winter mid-latitude ozone profiles, averaged
// over the layers between Pw levels (mb, kts:kte+1, top-down: Pw(kts) the top)
void o3data(const RrtmLwTables& t, double* O3PROF, const double* Pw, int kts, int kte) {
  double PRLEVH[1024], PPWRKH[33], O3WRK[32], PPWRK[32], O3ANN[32], PPANN[32];
  const double* O3SUM = t.o3sum - 1;
  const double* PPSUM = t.ppsum - 1;
  const double* O3WIN = t.o3win - 1;
  const double* PPWIN = t.ppwin - 1;
  for (int K = 1; K <= 31; ++K) PPANN[K] = PPSUM[K];
  O3ANN[1] = 0.5 * (O3SUM[1] + O3WIN[1]);
  for (int K = 2; K <= 31; ++K)
    O3ANN[K] = O3WIN[K - 1] + (O3WIN[K] - O3WIN[K - 1]) / (PPWIN[K] - PPWIN[K - 1]) * (PPSUM[K] - PPWIN[K - 1]);
  for (int K = 2; K <= 31; ++K) O3ANN[K] = 0.5 * (O3ANN[K] + O3SUM[K]);
  for (int K = 1; K <= 31; ++K) {
    O3WRK[K] = O3ANN[K];
    PPWRK[K] = PPANN[K];
  }
  for (int K = kts; K <= kte + 1; ++K) {
    const int NK = kte + 1 - K + 1;
    PRLEVH[K] = Pw[NK];
  }
  PPWRKH[1] = 1100.;
  for (int K = 2; K <= 31; ++K) PPWRKH[K] = (PPWRK[K] + PPWRK[K - 1]) / 2.;
  PPWRKH[32] = 0.;
  for (int K = kts; K <= kte; ++K) {
    for (int JJ = 1; JJ <= 31; ++JJ) {
      double PB1, PB2, PT1, PT2;
      if ((-(PRLEVH[K] - PPWRKH[JJ])) >= 0.)
        PB1 = 0.;
      else
        PB1 = PRLEVH[K] - PPWRKH[JJ];
      if ((-(PRLEVH[K] - PPWRKH[JJ + 1])) >= 0.)
        PB2 = 0.;
      else
        PB2 = PRLEVH[K] - PPWRKH[JJ + 1];
      if ((-(PRLEVH[K + 1] - PPWRKH[JJ])) >= 0.)
        PT1 = 0.;
      else
        PT1 = PRLEVH[K + 1] - PPWRKH[JJ];
      if ((-(PRLEVH[K + 1] - PPWRKH[JJ + 1])) >= 0.)
        PT2 = 0.;
      else
        PT2 = PRLEVH[K + 1] - PPWRKH[JJ + 1];
      O3PROF[K] = O3PROF[K] + (PB2 - PB1 - PT2 + PT1) * O3WRK[JJ];
    }
    O3PROF[K] = O3PROF[K] / (PRLEVH[K] - PRLEVH[K + 1]);
  }
}

// INIRAD: O3PROF(kts:kte) = 0, then O3DATA
void inirad(const RrtmLwTables& t, double* O3PROF, const double* Pw, int kts, int kte) {
  for (int K = kts; K <= kte; ++K) O3PROF[K] = 0.;
  o3data(t, O3PROF, Pw, kts, kte);
}

// MM5ATM: the RRTM column (bottom-up, L = 1 the lowest layer) from the top-down model arrays, the buffer layers
// above the model top, the gas amounts, the cloud optical depths
void mm5atm(const RrtmLwTables& t, Col& C, int kte, double co2vmr, double n2ovmr, double ch4vmr, double R,
            double G, double EMISS, double TSFC, double& TBOUND, double* SEMISS, std::string& message) {
  RrtmLwWorkspace& w = C.w;
  const int kts = 1, NLAYERS = C.nl;
  const double *CLDFRA = w.cldfra1d.data(), *O3PROF = w.o3prof.data(), *T = w.t1d.data(), *Tw = w.tw1d.data(),
               *QC = w.qc1d.data(), *QR = w.qr1d.data(), *QI = w.qi1d.data(), *QS = w.qs1d.data(),
               *P = w.p1d.data(), *Pw = w.pw1d.data(), *DELZ = w.dz1d.data();
  double* QV = w.qv1d.data();
  double *PAVEL = w.pavel.data(), *TAVEL = w.tavel.data(), *PZ = w.pz.data(), *TZ = w.tz.data(),
         *CLDFRAC = w.cldfrac.data(), *TAUCLOUD = w.taucloud.data(), *COLDRY = w.coldry.data();
  double *CLDFRC = w.cldfrc.data(), *O3 = w.o3.data(), *N2O = w.n2o.data(), *CH4 = w.ch4.data(),
         *CLWP = w.clwp.data(), *CIWP = w.ciwp.data(), *PLWP = w.plwp.data(), *PIWP = w.piwp.data();
  double *O3PROF2 = w.o3prof2.data(), *PZR = w.pzr.data(), *varint = w.varint.data();
  constexpr int nproflevs = 60;
  const double* PPROF = t.pprof - 1;
  const double* TPROF = t.tprof - 1;
  constexpr double amd = 28.9644, amw = 18.0154, avgdro = 6.022E23;
  constexpr double amdw = 1.607758, amdo = 0.603461, amdc = 1.805423, amdn = 0.658090;
  constexpr double h2ovmr = 5.00e-6;
  constexpr double ABCW = 0.144, ABICE = 0.0735, ABRN = 0.330E-3, ABSN = 2.34E-3;
  const double GRAVIT = G * 100.;
  // the K loop is vectorised by gfortran: the pow lanes of libmvec on pairs of levels from kts
  for (int K = kts; K <= kte; ++K) {
    const double RO = P[K] / (R * T[K]) * 100.;
    const double DZ = DELZ[K];
    QV[K] = c_max(QV[K], 1.E-12);
    CLDFRC[K] = CLDFRA[K];
    CLWP[K] = RO * QC[K] * DZ * 1000.;
    CIWP[K] = RO * QI[K] * DZ * 1000.;
    PLWP[K] = vm::pow(RO * QR[K], 0.75, K - kts, kte - kts + 1) * DZ * 1000.;
    PIWP[K] = vm::pow(RO * QS[K], 0.75, K - kts, kte - kts + 1) * DZ * 1000.;
    O3[K] = O3PROF[K];
    N2O[K] = 0.;
    CH4[K] = 0.;
  }
  for (int K = kts; K <= kte; ++K) N2O[K] = n2ovmr / amdn;
  for (int K = kts; K <= kte; ++K) CH4[K] = ch4vmr / amdc;
  for (int ILAY = kts; ILAY <= NLAYERS; ++ILAY) {
    for (int ISP = 1; ISP <= 35; ++ISP) C.WKL(ISP, ILAY) = 0.0;
    for (int ISP = 1; ISP <= MAXXSEC; ++ISP) C.WX(ISP, ILAY) = 0.0;
  }
  TBOUND = TSFC;
  if (TBOUND > 340.) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%10.3f", TBOUND);
    message = std::string("rrtm: TBOUND exceeds table limit: reset ") + buf;
    TBOUND = 339.99;
  }
  PZ[0] = Pw[kte + 1];
  TZ[0] = Tw[kte + 1];
  for (int L = 1; L <= kte; ++L) {
    PAVEL[L] = P[kte + 1 - L];
    TAVEL[L] = T[kte + 1 - L];
    PZ[L] = Pw[kte + 1 - L];
    TZ[L] = Tw[kte + 1 - L];
    C.WKL(1, L) = QV[kte + 1 - L] * amdw;
    C.WKL(2, L) = co2vmr;
    C.WKL(3, L) = O3[kte + 1 - L];
    C.WKL(3, L) = O3[L] * amdo;
    C.WKL(4, L) = N2O[kte + 1 - L] * amdn;
    C.WKL(6, L) = CH4[kte + 1 - L] * amdc;
    const double amm = (1 - C.WKL(1, L)) * amd + C.WKL(1, L) * amw;
    COLDRY[L] = (PZ[L - 1] - PZ[L]) * 1.E3 * avgdro / (GRAVIT * amm * (1 + C.WKL(1, L)));
  }
  for (int L = 1; L <= kte; ++L) {
    C.WX(2, L) = 0.;
    C.WX(3, L) = 0.;
  }
  // buffer layers: deltap mb thick up to 0 mb
  for (int L = kte + 1; L <= NLAYERS - 1; ++L) {
    PZ[L] = PZ[L - 1] - deltap;
    PAVEL[L] = 0.5 * (PZ[L] + PZ[L - 1]);
  }
  PZ[NLAYERS] = 0.00;
  PAVEL[NLAYERS] = 0.5 * (PZ[NLAYERS] + PZ[NLAYERS - 1]);
  // temperature of the buffer layers: the profile TPROF interpolated to PZ, shifted to the model top
  for (int L = 1; L <= NLAYERS; ++L) {
    int klev = 0;
    if (PPROF[nproflevs] < PZ[L]) {
      for (int LL = 2; LL <= nproflevs; ++LL) {
        if (PPROF[LL] < PZ[L]) {
          klev = LL - 1;
          break;
        }
      }
    } else {
      klev = nproflevs;
    }
    double vark, vark1, wght;
    if (klev != nproflevs) {
      vark = TPROF[klev];
      vark1 = TPROF[klev + 1];
      wght = (PZ[L] - PPROF[klev]) / (PPROF[klev + 1] - PPROF[klev]);
    } else {
      vark = TPROF[klev];
      vark1 = TPROF[klev];
      wght = 0.0;
    }
    varint[L] = wght * (vark1 - vark) + vark;
  }
  for (int L = kte + 1; L <= NLAYERS; ++L) {
    TZ[L] = varint[L] + (TZ[kte] - varint[kte]);
    TAVEL[L] = 0.5 * (TZ[L] + TZ[L - 1]);
  }
  // ozone of the extended column
  for (int L = kts; L <= NLAYERS; ++L) {
    const int klev = NLAYERS - L + 1;
    PZR[L] = PZ[klev];
  }
  inirad(t, O3PROF2, PZR, kts, NLAYERS - 1);
  O3PROF2[NLAYERS] = 6.135E-6;
  if (kte != NLAYERS) {
    for (int L = 1; L <= NLAYERS; ++L) {
      C.WKL(3, L) = O3PROF2[L] * amdo;
      if (L > kte) {
        C.WKL(1, L) = h2ovmr;
        C.WKL(2, L) = co2vmr;
        C.WKL(4, L) = C.WKL(4, kte);
        C.WKL(6, L) = C.WKL(6, kte);
        const double amm = (1 - C.WKL(1, L)) * amd + C.WKL(1, L) * amw;
        COLDRY[L] = (PZ[L - 1] - PZ[L]) * 1.E3 * avgdro / (GRAVIT * amm * (1 + C.WKL(1, L)));
        C.WX(2, L) = C.WX(2, kte);
        C.WX(3, L) = C.WX(3, kte);
      }
    }
  }
  for (int L = 1; L <= NLAYERS; ++L) {
    for (int IMOL = 1; IMOL <= NMOL; ++IMOL) C.WKL(IMOL, L) = COLDRY[L] * C.WKL(IMOL, L);
    for (int IX = 1; IX <= MAXXSEC; ++IX)
      if (t.ixindx[IX - 1] != 0) C.WX(t.ixindx[IX - 1], L) = COLDRY[L] * C.WX(IX, L) * 1.E-20;
  }
  for (int N = 1; N <= NBANDS; ++N) SEMISS[N] = EMISS;
  for (int L = 1; L <= kte; ++L) {
    TAUCLOUD[L] = ABCW * CLWP[kte + 1 - L] + ABICE * CIWP[kte + 1 - L] + ABRN * PLWP[kte + 1 - L] +
                  ABSN * PIWP[kte + 1 - L];
    if (TAUCLOUD[L] > 0.01) CLDFRC[kte + 1 - L] = 1.;
    CLDFRAC[L] = CLDFRC[kte + 1 - L];
  }
  for (int L = kte + 1; L <= NLAYERS; ++L) {
    CLDFRAC[L] = 0.0;
    TAUCLOUD[L] = 0.0;
  }
}

// SETCOEF: the reference-pressure / temperature indices and interpolation factors, the column amounts
void setcoef(const RrtmLwTables& t, Col& C, int& LAYTROP, int& LAYSWTCH, int& LAYLOW) {
  RrtmLwWorkspace& w = C.w;
  const int NLAYERS = C.nl;
  const double* PREFLOG = t.preflog - 1;
  const double* TREF = t.tref - 1;
  const double STPFAC = 296. / 1013.;
  LAYTROP = 0;
  LAYSWTCH = 0;
  LAYLOW = 0;
  for (int LAY = 1; LAY <= NLAYERS; ++LAY) {
    const double PLOG = std::log(w.pavel[LAY]);
    w.jp[LAY] = int(36. - 5 * (PLOG + 0.04));
    if (w.jp[LAY] < 1)
      w.jp[LAY] = 1;
    else if (w.jp[LAY] > 58)
      w.jp[LAY] = 58;
    const int JP1 = w.jp[LAY] + 1;
    const double FP = 5. * (PREFLOG[w.jp[LAY]] - PLOG);
    w.jt[LAY] = int(3. + (w.tavel[LAY] - TREF[w.jp[LAY]]) / 15.);
    if (w.jt[LAY] < 1)
      w.jt[LAY] = 1;
    else if (w.jt[LAY] > 4)
      w.jt[LAY] = 4;
    const double FT = ((w.tavel[LAY] - TREF[w.jp[LAY]]) / 15.) - double(w.jt[LAY] - 3);
    w.jt1[LAY] = int(3. + (w.tavel[LAY] - TREF[JP1]) / 15.);
    if (w.jt1[LAY] < 1)
      w.jt1[LAY] = 1;
    else if (w.jt1[LAY] > 4)
      w.jt1[LAY] = 4;
    const double FT1 = ((w.tavel[LAY] - TREF[JP1]) / 15.) - double(w.jt1[LAY] - 3);
    const double WATER = C.WKL(1, LAY) / w.coldry[LAY];
    const double SCALEFAC = w.pavel[LAY] * STPFAC / w.tavel[LAY];
    if (PLOG > 4.56) {  // troposphere (else GO TO 5300)
      LAYTROP = LAYTROP + 1;
      if (PLOG > 5.76) LAYSWTCH = LAYSWTCH + 1;
      if (PLOG >= 6.62) LAYLOW = LAYLOW + 1;
      w.forfac[LAY] = SCALEFAC / (1. + WATER);
      w.selffac[LAY] = WATER * w.forfac[LAY];
      const double FACTOR = (w.tavel[LAY] - 188.0) / 7.2;
      w.indself[LAY] = std::min(9, std::max(1, int(FACTOR) - 7));
      w.selffrac[LAY] = FACTOR - double(w.indself[LAY] + 7);
    } else {  // 5300
      w.forfac[LAY] = SCALEFAC / (1. + WATER);
    }
    w.colh2o[LAY] = 1.E-20 * C.WKL(1, LAY);
    w.colco2[LAY] = 1.E-20 * C.WKL(2, LAY);
    w.colo3[LAY] = 1.E-20 * C.WKL(3, LAY);
    w.coln2o[LAY] = 1.E-20 * C.WKL(4, LAY);
    w.colch4[LAY] = 1.E-20 * C.WKL(6, LAY);
    w.colo2[LAY] = 1.E-20 * C.WKL(7, LAY);
    if (w.colco2[LAY] == 0.) w.colco2[LAY] = 1.E-32 * w.coldry[LAY];
    if (w.coln2o[LAY] == 0.) w.coln2o[LAY] = 1.E-32 * w.coldry[LAY];
    if (w.colch4[LAY] == 0.) w.colch4[LAY] = 1.E-32 * w.coldry[LAY];
    const double CO2REG = 3.55E-24 * w.coldry[LAY];
    w.co2mult[LAY] =
        (w.colco2[LAY] - CO2REG) * 272.63 * std::exp(-1919.4 / w.tavel[LAY]) / (8.7604E-4 * w.tavel[LAY]);
    // 5400
    const double COMPFP = 1. - FP;
    w.fac10[LAY] = COMPFP * FT;
    w.fac00[LAY] = COMPFP * (1. - FT);
    w.fac11[LAY] = FP * FT1;
    w.fac01[LAY] = FP * (1. - FT1);
  }
  if (LAYLOW == 0) LAYLOW = 1;
  if (w.jp[LAYSWTCH + 1] <= 6) LAYSWTCH = LAYSWTCH + 1;
}

// TAUGB1..16: the gaseous optical depths TAUG and Planck fractions PFRAC of the 16 bands (band n's g-points at
// NGS(n-1)+IG). The sums follow the Fortran left to right.
struct Taugb {
  const RrtmLwTables& t;
  Col& C;
  RrtmLwWorkspace& w;
  int NLAYERS, LAYTROP, LAYSWTCH, LAYLOW;

  // the eight-coefficient interpolation of the binary-species bands (JS, JS+1 in the species index)
  struct Fac8 {
    double f000, f010, f100, f110, f001, f011, f101, f111;
  };
  Fac8 fac8(int LAY, double FS) const {
    return {(1. - FS) * w.fac00[LAY], (1. - FS) * w.fac10[LAY], FS * w.fac00[LAY], FS * w.fac10[LAY],
            (1. - FS) * w.fac01[LAY], (1. - FS) * w.fac11[LAY], FS * w.fac01[LAY], FS * w.fac11[LAY]};
  }
  // FAC000*A(IND0) + FAC100*A(IND0+1) + FAC010*A(IND0+s) + FAC110*A(IND0+s+1) + (the same at IND1)
  static double sum8(const Fac8& f, const double* A, int n1, int IND0, int IND1, int s, int IG) {
    return f.f000 * F2(A, n1, IND0, IG) + f.f100 * F2(A, n1, IND0 + 1, IG) + f.f010 * F2(A, n1, IND0 + s, IG) +
           f.f110 * F2(A, n1, IND0 + s + 1, IG) + f.f001 * F2(A, n1, IND1, IG) + f.f101 * F2(A, n1, IND1 + 1, IG) +
           f.f011 * F2(A, n1, IND1 + s, IG) + f.f111 * F2(A, n1, IND1 + s + 1, IG);
  }
  // FAC00*A(IND0) + FAC10*A(IND0+1) + FAC01*A(IND1) + FAC11*A(IND1+1)
  double sum4(int LAY, const double* A, int n1, int IND0, int IND1, int IG) const {
    return w.fac00[LAY] * F2(A, n1, IND0, IG) + w.fac10[LAY] * F2(A, n1, IND0 + 1, IG) +
           w.fac01[LAY] * F2(A, n1, IND1, IG) + w.fac11[LAY] * F2(A, n1, IND1 + 1, IG);
  }
  // SELFREFC(INDS,IG) + SELFFRAC*(SELFREFC(INDS+1,IG) - SELFREFC(INDS,IG))
  double self(int LAY, const double* S, int INDS, int IG) const {
    return F2(S, 10, INDS, IG) + w.selffrac[LAY] * (F2(S, 10, INDS + 1, IG) - F2(S, 10, INDS, IG));
  }
  // PFRACn = FRAC(IG,JS) + FS*(FRAC(IG,JS+1) - FRAC(IG,JS))
  static double frac(const double* F, int ng, int IG, int JS, double FS) {
    return F2(F, ng, IG, JS) + FS * (F2(F, ng, IG, JS + 1) - F2(F, ng, IG, JS));
  }
  int ind0a(int LAY, int b) const { return ((w.jp[LAY] - 1) * 5 + (w.jt[LAY] - 1)) * t.nspa[b - 1]; }
  int ind1a(int LAY, int b) const { return (w.jp[LAY] * 5 + (w.jt1[LAY] - 1)) * t.nspa[b - 1]; }
  int ind0b(int LAY, int b) const { return ((w.jp[LAY] - 13) * 5 + (w.jt[LAY] - 1)) * t.nspb[b - 1]; }
  int ind1b(int LAY, int b) const { return ((w.jp[LAY] - 12) * 5 + (w.jt1[LAY] - 1)) * t.nspb[b - 1]; }
  // SPECPARM = num / (num + STRRAT*den) capped at ONEMINUS; JS = 1 + INT(mult*SPECPARM), FS = MOD(.., 1.0)
  static void spec(double num, double SPECCOMB, double mult, int& JS, double& FS) {
    double SPECPARM = num / SPECCOMB;
    if (SPECPARM >= ONEMINUS) SPECPARM = ONEMINUS;
    const double SPECMULT = mult * (SPECPARM);
    JS = 1 + int(SPECMULT);
    FS = std::fmod(SPECMULT, 1.0);
  }

  void taugb1() {
    for (int LAY = 1; LAY <= LAYTROP; ++LAY) {
      const int IND0 = ind0a(LAY, 1) + 1, IND1 = ind1a(LAY, 1) + 1, INDS = w.indself[LAY];
      for (int IG = 1; IG <= NG1; ++IG) {
        C.TAUG(IG, LAY) = w.colh2o[LAY] * (sum4(LAY, t.absa1, 65, IND0, IND1, IG) +
                                           w.selffac[LAY] * self(LAY, t.selfrefc1, INDS, IG) +
                                           w.forfac[LAY] * t.forrefc1[IG - 1]);
        C.PFRAC(IG, LAY) = t.fracrefac1[IG - 1];
      }
    }
    for (int LAY = LAYTROP + 1; LAY <= NLAYERS; ++LAY) {
      const int IND0 = ind0b(LAY, 1) + 1, IND1 = ind1b(LAY, 1) + 1;
      for (int IG = 1; IG <= NG1; ++IG) {
        C.TAUG(IG, LAY) =
            w.colh2o[LAY] * (sum4(LAY, t.absb1, 235, IND0, IND1, IG) + w.forfac[LAY] * t.forrefc1[IG - 1]);
        C.PFRAC(IG, LAY) = t.fracrefbc1[IG - 1];
      }
    }
  }

  void taugb2() {
    constexpr int NGS1 = 8;
    const double* REFPARAM = t.refparam2 - 1;
    double FC00, FC01, FC10, FC11;
    for (int LAY = 1; LAY <= LAYTROP; ++LAY) {
      const double WATER = 1.E20 * w.colh2o[LAY] / w.coldry[LAY];
      const double H2OPARAM = WATER / (WATER + .002);
      int IFRAC = 2;
      for (; IFRAC <= 12; ++IFRAC)
        if (H2OPARAM >= REFPARAM[IFRAC]) break;  // GO TO 1900; else the loop ends with IFRAC = 13
      const double FRACINT = (H2OPARAM - REFPARAM[IFRAC]) / (REFPARAM[IFRAC - 1] - REFPARAM[IFRAC]);
      const double FP = w.fac11[LAY] + w.fac01[LAY];
      int IFP = int(2.E2 * FP + 0.5);
      if (IFP <= 0) IFP = 0;
      FC00 = w.fac00[LAY] * t.corr2[IFP];
      FC10 = w.fac10[LAY] * t.corr2[IFP];
      FC01 = w.fac01[LAY] * t.corr1[IFP];
      FC11 = w.fac11[LAY] * t.corr1[IFP];
      const int IND0 = ind0a(LAY, 2) + 1, IND1 = ind1a(LAY, 2) + 1, INDS = w.indself[LAY];
      for (int IG = 1; IG <= NG2; ++IG) {
        C.TAUG(NGS1 + IG, LAY) =
            w.colh2o[LAY] * (FC00 * F2(t.absa2, 65, IND0, IG) + FC10 * F2(t.absa2, 65, IND0 + 1, IG) +
                             FC01 * F2(t.absa2, 65, IND1, IG) + FC11 * F2(t.absa2, 65, IND1 + 1, IG) +
                             w.selffac[LAY] * self(LAY, t.selfrefc2, INDS, IG) + w.forfac[LAY] * t.forrefc2[IG - 1]);
        C.PFRAC(NGS1 + IG, LAY) =
            F2(t.fracrefac2, NG2, IG, IFRAC) +
            FRACINT * (F2(t.fracrefac2, NG2, IG, IFRAC - 1) - F2(t.fracrefac2, NG2, IG, IFRAC));
      }
    }
    for (int LAY = LAYTROP + 1; LAY <= NLAYERS; ++LAY) {
      const double FP = w.fac11[LAY] + w.fac01[LAY];
      int IFP = int(2.E2 * FP + 0.5);
      if (IFP <= 0) IFP = 0;
      FC00 = w.fac00[LAY] * t.corr2[IFP];
      FC10 = w.fac10[LAY] * t.corr2[IFP];
      FC01 = w.fac01[LAY] * t.corr1[IFP];
      FC11 = w.fac11[LAY] * t.corr1[IFP];
      const int IND0 = ind0b(LAY, 2) + 1, IND1 = ind1b(LAY, 2) + 1;
      for (int IG = 1; IG <= NG2; ++IG) {
        C.TAUG(NGS1 + IG, LAY) =
            w.colh2o[LAY] * (FC00 * F2(t.absb2, 235, IND0, IG) + FC10 * F2(t.absb2, 235, IND0 + 1, IG) +
                             FC01 * F2(t.absb2, 235, IND1, IG) + FC11 * F2(t.absb2, 235, IND1 + 1, IG) +
                             w.forfac[LAY] * t.forrefc2[IG - 1]);
        C.PFRAC(NGS1 + IG, LAY) = t.fracrefbc2[IG - 1];
      }
    }
  }

  void taugb3() {
    constexpr int NGS2 = 22;
    const double *H2OREF = t.h2oref3 - 1, *N2OREF = t.n2oref3 - 1, *CO2REF = t.co2ref3 - 1,
                 *ETAREF = t.etaref3 - 1;
    const double STRRAT = 1.19268;
    for (int LAY = 1; LAY <= LAYTROP; ++LAY) {
      const double SPECCOMB = w.colh2o[LAY] + STRRAT * w.colco2[LAY];
      int JS;
      double FS;
      spec(w.colh2o[LAY], SPECCOMB, 8., JS, FS);
      if (JS == 8) {
        if (FS >= 0.9) {
          JS = 9;
          FS = 10. * (FS - 0.9);
        } else {
          FS = FS / 0.9;
        }
      }
      const int NS = JS + int(FS + 0.5);
      const double FP = w.fac01[LAY] + w.fac11[LAY];
      const Fac8 f = fac8(LAY, FS);
      const int IND0 = ind0a(LAY, 3) + JS, IND1 = ind1a(LAY, 3) + JS, INDS = w.indself[LAY];
      const double COLREF1 = N2OREF[w.jp[LAY]], COLREF2 = N2OREF[w.jp[LAY] + 1];
      double WCOMB1, WCOMB2;
      if (NS == 10) {
        WCOMB1 = H2OREF[w.jp[LAY]];
        WCOMB2 = H2OREF[w.jp[LAY] + 1];
      } else {
        WCOMB1 = STRRAT * CO2REF[w.jp[LAY]] / (1. - ETAREF[NS]);
        WCOMB2 = STRRAT * CO2REF[w.jp[LAY] + 1] / (1. - ETAREF[NS]);
      }
      const double RATIO = (COLREF1 / WCOMB1) + FP * ((COLREF2 / WCOMB2) - (COLREF1 / WCOMB1));
      const double CURRN2O = SPECCOMB * RATIO;
      const double N2OMULT = w.coln2o[LAY] - CURRN2O;
      for (int IG = 1; IG <= NG3; ++IG) {
        C.TAUG(NGS2 + IG, LAY) = SPECCOMB * sum8(f, t.absa3, 650, IND0, IND1, 10, IG) +
                                 w.colh2o[LAY] * (w.selffac[LAY] * self(LAY, t.selfrefc3, INDS, IG) +
                                                  w.forfac[LAY] * t.forrefc3[IG - 1]) +
                                 N2OMULT * t.absn2oac3[IG - 1];
        C.PFRAC(NGS2 + IG, LAY) = frac(t.fracrefac3, NG3, IG, JS, FS);
      }
    }
    for (int LAY = LAYTROP + 1; LAY <= NLAYERS; ++LAY) {
      const double SPECCOMB = w.colh2o[LAY] + STRRAT * w.colco2[LAY];
      int JS;
      double FS;
      spec(w.colh2o[LAY], SPECCOMB, 4., JS, FS);
      const int NS = JS + int(FS + 0.5);
      const double FP = w.fac01[LAY] + w.fac11[LAY];
      const Fac8 f = fac8(LAY, FS);
      const int IND0 = ind0b(LAY, 3) + JS, IND1 = ind1b(LAY, 3) + JS;
      const double COLREF1 = N2OREF[w.jp[LAY]], COLREF2 = N2OREF[w.jp[LAY] + 1];
      double WCOMB1, WCOMB2;
      if (NS == 5) {
        WCOMB1 = H2OREF[w.jp[LAY]];
        WCOMB2 = H2OREF[w.jp[LAY] + 1];
      } else {
        WCOMB1 = STRRAT * CO2REF[w.jp[LAY]] / (1. - ETAREF[NS]);
        WCOMB2 = STRRAT * CO2REF[w.jp[LAY] + 1] / (1. - ETAREF[NS]);
      }
      const double RATIO = (COLREF1 / WCOMB1) + FP * ((COLREF2 / WCOMB2) - (COLREF1 / WCOMB1));
      const double CURRN2O = SPECCOMB * RATIO;
      const double N2OMULT = w.coln2o[LAY] - CURRN2O;
      for (int IG = 1; IG <= NG3; ++IG) {
        C.TAUG(NGS2 + IG, LAY) = SPECCOMB * sum8(f, t.absb3, 1175, IND0, IND1, 5, IG) +
                                 w.colh2o[LAY] * w.forfac[LAY] * t.forrefc3[IG - 1] +
                                 N2OMULT * t.absn2obc3[IG - 1];
        C.PFRAC(NGS2 + IG, LAY) = frac(t.fracrefbc3, NG3, IG, JS, FS);
      }
    }
  }

  void taugb4() {
    constexpr int NGS3 = 38;
    const double STRRAT1 = 850.577, STRRAT2 = 35.7416;
    for (int LAY = 1; LAY <= LAYTROP; ++LAY) {
      const double SPECCOMB = w.colh2o[LAY] + STRRAT1 * w.colco2[LAY];
      int JS;
      double FS;
      spec(w.colh2o[LAY], SPECCOMB, 8., JS, FS);
      const Fac8 f = fac8(LAY, FS);
      const int IND0 = ind0a(LAY, 4) + JS, IND1 = ind1a(LAY, 4) + JS, INDS = w.indself[LAY];
      for (int IG = 1; IG <= NG4; ++IG) {
        C.TAUG(NGS3 + IG, LAY) = SPECCOMB * sum8(f, t.absa4, 585, IND0, IND1, 9, IG) +
                                 w.colh2o[LAY] * w.selffac[LAY] * self(LAY, t.selfrefc4, INDS, IG);
        C.PFRAC(NGS3 + IG, LAY) = frac(t.fracrefac4, NG4, IG, JS, FS);
      }
    }
    for (int LAY = LAYTROP + 1; LAY <= NLAYERS; ++LAY) {
      const double SPECCOMB = w.colo3[LAY] + STRRAT2 * w.colco2[LAY];
      int JS;
      double FS;
      spec(w.colo3[LAY], SPECCOMB, 4., JS, FS);
      if (JS > 1) {
        JS = JS + 1;
      } else if (FS >= 0.0024) {
        JS = 2;
        FS = (FS - 0.0024) / 0.9976;
      } else {
        JS = 1;
        FS = FS / 0.0024;
      }
      const Fac8 f = fac8(LAY, FS);
      const int IND0 = ind0b(LAY, 4) + JS, IND1 = ind1b(LAY, 4) + JS;
      for (int IG = 1; IG <= NG4; ++IG) {
        C.TAUG(NGS3 + IG, LAY) = SPECCOMB * sum8(f, t.absb4, 1410, IND0, IND1, 6, IG);
        C.PFRAC(NGS3 + IG, LAY) = frac(t.fracrefbc4, NG4, IG, JS, FS);
      }
    }
  }

  void taugb5() {
    constexpr int NGS4 = 52;
    const double STRRAT1 = 90.4894, STRRAT2 = 0.900502;
    for (int LAY = 1; LAY <= LAYTROP; ++LAY) {
      const double SPECCOMB = w.colh2o[LAY] + STRRAT1 * w.colco2[LAY];
      int JS;
      double FS;
      spec(w.colh2o[LAY], SPECCOMB, 8., JS, FS);
      const Fac8 f = fac8(LAY, FS);
      const int IND0 = ind0a(LAY, 5) + JS, IND1 = ind1a(LAY, 5) + JS, INDS = w.indself[LAY];
      for (int IG = 1; IG <= NG5; ++IG) {
        C.TAUG(NGS4 + IG, LAY) = SPECCOMB * sum8(f, t.absa5, 585, IND0, IND1, 9, IG) +
                                 w.colh2o[LAY] * w.selffac[LAY] * self(LAY, t.selfrefc5, INDS, IG) +
                                 C.WX(1, LAY) * t.ccl4c5[IG - 1];
        C.PFRAC(NGS4 + IG, LAY) = frac(t.fracrefac5, NG5, IG, JS, FS);
      }
    }
    for (int LAY = LAYTROP + 1; LAY <= NLAYERS; ++LAY) {
      const double SPECCOMB = w.colo3[LAY] + STRRAT2 * w.colco2[LAY];
      int JS;
      double FS;
      spec(w.colo3[LAY], SPECCOMB, 4., JS, FS);
      const Fac8 f = fac8(LAY, FS);
      const int IND0 = ind0b(LAY, 5) + JS, IND1 = ind1b(LAY, 5) + JS;
      for (int IG = 1; IG <= NG5; ++IG) {
        C.TAUG(NGS4 + IG, LAY) =
            SPECCOMB * sum8(f, t.absb5, 1175, IND0, IND1, 5, IG) + C.WX(1, LAY) * t.ccl4c5[IG - 1];
        C.PFRAC(NGS4 + IG, LAY) = frac(t.fracrefbc5, NG5, IG, JS, FS);
      }
    }
  }

  void taugb6() {
    constexpr int NGS5 = 68;
    for (int LAY = 1; LAY <= LAYTROP; ++LAY) {
      const int IND0 = ind0a(LAY, 6) + 1, IND1 = ind1a(LAY, 6) + 1, INDS = w.indself[LAY];
      for (int IG = 1; IG <= NG6; ++IG) {
        C.TAUG(NGS5 + IG, LAY) =
            w.colh2o[LAY] *
                (sum4(LAY, t.absa6, 65, IND0, IND1, IG) + w.selffac[LAY] * self(LAY, t.selfrefc6, INDS, IG)) +
            C.WX(2, LAY) * t.cfc11adjc6[IG - 1] + C.WX(3, LAY) * t.cfc12c6[IG - 1] +
            w.co2mult[LAY] * t.absco2c6[IG - 1];
        C.PFRAC(NGS5 + IG, LAY) = t.fracrefac6[IG - 1];
      }
    }
    for (int LAY = LAYTROP + 1; LAY <= NLAYERS; ++LAY) {
      for (int IG = 1; IG <= NG6; ++IG) {
        C.TAUG(NGS5 + IG, LAY) = 0.0 + C.WX(2, LAY) * t.cfc11adjc6[IG - 1] + C.WX(3, LAY) * t.cfc12c6[IG - 1];
        C.PFRAC(NGS5 + IG, LAY) = t.fracrefac6[IG - 1];
      }
    }
  }

  void taugb7() {
    constexpr int NGS6 = 76;
    const double STRRAT1 = 8.21104E4;
    for (int LAY = 1; LAY <= LAYTROP; ++LAY) {
      const double SPECCOMB = w.colh2o[LAY] + STRRAT1 * w.colo3[LAY];
      int JS;
      double FS;
      spec(w.colh2o[LAY], SPECCOMB, 8., JS, FS);
      const Fac8 f = fac8(LAY, FS);
      const int IND0 = ind0a(LAY, 7) + JS, IND1 = ind1a(LAY, 7) + JS, INDS = w.indself[LAY];
      for (int IG = 1; IG <= NG7; ++IG) {
        C.TAUG(NGS6 + IG, LAY) = SPECCOMB * sum8(f, t.absa7, 585, IND0, IND1, 9, IG) +
                                 w.colh2o[LAY] * w.selffac[LAY] * self(LAY, t.selfrefc7, INDS, IG) +
                                 w.co2mult[LAY] * t.absco2c7[IG - 1];
        C.PFRAC(NGS6 + IG, LAY) = frac(t.fracrefac7, NG7, IG, JS, FS);
      }
    }
    for (int LAY = LAYTROP + 1; LAY <= NLAYERS; ++LAY) {
      const int IND0 = ind0b(LAY, 7) + 1, IND1 = ind1b(LAY, 7) + 1;
      for (int IG = 1; IG <= NG7; ++IG) {
        C.TAUG(NGS6 + IG, LAY) =
            w.colo3[LAY] * sum4(LAY, t.absb7, 235, IND0, IND1, IG) + w.co2mult[LAY] * t.absco2c7[IG - 1];
        C.PFRAC(NGS6 + IG, LAY) = t.fracrefbc7[IG - 1];
      }
    }
  }

  void taugb8() {
    constexpr int NGS7 = 88;
    const double *H2OREF = t.h2oref8 - 1, *N2OREF = t.n2oref8 - 1, *O3REF = t.o3ref8 - 1;
    for (int LAY = 1; LAY <= LAYSWTCH; ++LAY) {
      const double FP = w.fac01[LAY] + w.fac11[LAY];
      const int IND0 = ind0a(LAY, 8) + 1, IND1 = ind1a(LAY, 8) + 1, INDS = w.indself[LAY];
      const double COLREF1 = N2OREF[w.jp[LAY]], COLREF2 = N2OREF[w.jp[LAY] + 1];
      const double WCOMB1 = H2OREF[w.jp[LAY]], WCOMB2 = H2OREF[w.jp[LAY] + 1];
      const double RATIO = (COLREF1 / WCOMB1) + FP * ((COLREF2 / WCOMB2) - (COLREF1 / WCOMB1));
      const double CURRN2O = w.colh2o[LAY] * RATIO;
      const double N2OMULT = w.coln2o[LAY] - CURRN2O;
      for (int IG = 1; IG <= NG8; ++IG) {
        C.TAUG(NGS7 + IG, LAY) =
            w.colh2o[LAY] *
                (sum4(LAY, t.absa8, 35, IND0, IND1, IG) + w.selffac[LAY] * self(LAY, t.selfrefc8, INDS, IG)) +
            C.WX(3, LAY) * t.cfc12c8[IG - 1] + C.WX(4, LAY) * t.cfc22adjc8[IG - 1] +
            w.co2mult[LAY] * t.absco2ac8[IG - 1] + N2OMULT * t.absn2oac8[IG - 1];
        C.PFRAC(NGS7 + IG, LAY) = t.fracrefac8[IG - 1];
      }
    }
    for (int LAY = LAYSWTCH + 1; LAY <= NLAYERS; ++LAY) {
      const double FP = w.fac01[LAY] + w.fac11[LAY];
      const int IND0 = ((w.jp[LAY] - 7) * 5 + (w.jt[LAY] - 1)) * t.nspb[7] + 1;
      const int IND1 = ((w.jp[LAY] - 6) * 5 + (w.jt1[LAY] - 1)) * t.nspb[7] + 1;
      const double COLREF1 = N2OREF[w.jp[LAY]], COLREF2 = N2OREF[w.jp[LAY] + 1];
      const double WCOMB1 = O3REF[w.jp[LAY]], WCOMB2 = O3REF[w.jp[LAY] + 1];
      const double RATIO = (COLREF1 / WCOMB1) + FP * ((COLREF2 / WCOMB2) - (COLREF1 / WCOMB1));
      const double CURRN2O = w.colo3[LAY] * RATIO;
      const double N2OMULT = w.coln2o[LAY] - CURRN2O;
      for (int IG = 1; IG <= NG8; ++IG) {
        C.TAUG(NGS7 + IG, LAY) = w.colo3[LAY] * sum4(LAY, t.absb8, 265, IND0, IND1, IG) +
                                 C.WX(3, LAY) * t.cfc12c8[IG - 1] + C.WX(4, LAY) * t.cfc22adjc8[IG - 1] +
                                 w.co2mult[LAY] * t.absco2bc8[IG - 1] + N2OMULT * t.absn2obc8[IG - 1];
        C.PFRAC(NGS7 + IG, LAY) = t.fracrefbc8[IG - 1];
      }
    }
  }

  void taugb9() {
    constexpr int NGS8 = 96;
    const double *N2OREF = t.n2oref9 - 1, *H2OREF = t.h2oref9 - 1, *CH4REF = t.ch4ref9 - 1,
                 *ETAREF = t.etaref9 - 1;
    const double STRRAT = 21.6282;
    int IOFF = 0;
    for (int LAY = 1; LAY <= LAYTROP; ++LAY) {
      const double SPECCOMB = w.colh2o[LAY] + STRRAT * w.colch4[LAY];
      int JS;
      double FS;
      spec(w.colh2o[LAY], SPECCOMB, 8., JS, FS);
      int JFRAC = JS;
      double FFRAC = FS;
      if (JS == 8) {
        if (FS <= 0.68) {
          FS = FS / 0.68;
        } else if (FS <= 0.92) {
          JS = JS + 1;
          FS = (FS - 0.68) / 0.24;
        } else {
          JS = JS + 2;
          FS = (FS - 0.92) / 0.08;
        }
      } else if (JS == 9) {
        JS = 10;
        FS = 1.;
        JFRAC = 8;
        FFRAC = 1.;
      }
      const double FP = w.fac01[LAY] + w.fac11[LAY];
      const int NS = JS + int(FS + 0.5);
      const Fac8 f = fac8(LAY, FS);
      const int IND0 = ind0a(LAY, 9) + JS, IND1 = ind1a(LAY, 9) + JS, INDS = w.indself[LAY];
      if (LAY == LAYLOW) IOFF = NG9;
      if (LAY == LAYSWTCH) IOFF = 2 * NG9;
      const double COLREF1 = N2OREF[w.jp[LAY]], COLREF2 = N2OREF[w.jp[LAY] + 1];
      double WCOMB1, WCOMB2;
      if (NS == 11) {
        WCOMB1 = H2OREF[w.jp[LAY]];
        WCOMB2 = H2OREF[w.jp[LAY] + 1];
      } else {
        WCOMB1 = STRRAT * CH4REF[w.jp[LAY]] / (1. - ETAREF[NS]);
        WCOMB2 = STRRAT * CH4REF[w.jp[LAY] + 1] / (1. - ETAREF[NS]);
      }
      const double RATIO = (COLREF1 / WCOMB1) + FP * ((COLREF2 / WCOMB2) - (COLREF1 / WCOMB1));
      const double CURRN2O = SPECCOMB * RATIO;
      const double N2OMULT = w.coln2o[LAY] - CURRN2O;
      for (int IG = 1; IG <= NG9; ++IG) {
        C.TAUG(NGS8 + IG, LAY) = SPECCOMB * sum8(f, t.absa9, 715, IND0, IND1, 11, IG) +
                                 w.colh2o[LAY] * w.selffac[LAY] * self(LAY, t.selfrefc9, INDS, IG) +
                                 N2OMULT * t.absn2oc9[IG + IOFF - 1];
        C.PFRAC(NGS8 + IG, LAY) = frac(t.fracrefac9, NG9, IG, JFRAC, FFRAC);
      }
    }
    for (int LAY = LAYTROP + 1; LAY <= NLAYERS; ++LAY) {
      const int IND0 = ind0b(LAY, 9) + 1, IND1 = ind1b(LAY, 9) + 1;
      for (int IG = 1; IG <= NG9; ++IG) {
        C.TAUG(NGS8 + IG, LAY) = w.colch4[LAY] * sum4(LAY, t.absb9, 235, IND0, IND1, IG);
        C.PFRAC(NGS8 + IG, LAY) = t.fracrefbc9[IG - 1];
      }
    }
  }

  void taugb10() {
    constexpr int NGS9 = 108;
    for (int LAY = 1; LAY <= LAYTROP; ++LAY) {
      const int IND0 = ind0a(LAY, 10) + 1, IND1 = ind1a(LAY, 10) + 1;
      for (int IG = 1; IG <= NG10; ++IG) {
        C.TAUG(NGS9 + IG, LAY) = w.colh2o[LAY] * sum4(LAY, t.absa10, 65, IND0, IND1, IG);
        C.PFRAC(NGS9 + IG, LAY) = t.fracrefac10[IG - 1];
      }
    }
    for (int LAY = LAYTROP + 1; LAY <= NLAYERS; ++LAY) {
      const int IND0 = ind0b(LAY, 10) + 1, IND1 = ind1b(LAY, 10) + 1;
      for (int IG = 1; IG <= NG10; ++IG) {
        C.TAUG(NGS9 + IG, LAY) = w.colh2o[LAY] * sum4(LAY, t.absb10, 235, IND0, IND1, IG);
        C.PFRAC(NGS9 + IG, LAY) = t.fracrefbc10[IG - 1];
      }
    }
  }

  void taugb11() {
    constexpr int NGS10 = 114;
    for (int LAY = 1; LAY <= LAYTROP; ++LAY) {
      const int IND0 = ind0a(LAY, 11) + 1, IND1 = ind1a(LAY, 11) + 1, INDS = w.indself[LAY];
      for (int IG = 1; IG <= NG11; ++IG) {
        C.TAUG(NGS10 + IG, LAY) = w.colh2o[LAY] * (sum4(LAY, t.absa11, 65, IND0, IND1, IG) +
                                                   w.selffac[LAY] * self(LAY, t.selfrefc11, INDS, IG));
        C.PFRAC(NGS10 + IG, LAY) = t.fracrefac11[IG - 1];
      }
    }
    for (int LAY = LAYTROP + 1; LAY <= NLAYERS; ++LAY) {
      const int IND0 = ind0b(LAY, 11) + 1, IND1 = ind1b(LAY, 11) + 1;
      for (int IG = 1; IG <= NG11; ++IG) {
        C.TAUG(NGS10 + IG, LAY) = w.colh2o[LAY] * sum4(LAY, t.absb11, 235, IND0, IND1, IG);
        C.PFRAC(NGS10 + IG, LAY) = t.fracrefbc11[IG - 1];
      }
    }
  }

  // TAUGB12 / 13 / 15 / 16: a binary-species lower atmosphere, 0 above LAYTROP
  void taugb_lower(int NGS, int ng, int b, const double* absa, const double* selfc, const double* fracc,
                   double strrat, const std::vector<double>& num, const std::vector<double>& den,
                   const std::vector<double>& self_col) {
    for (int LAY = 1; LAY <= LAYTROP; ++LAY) {
      const double SPECCOMB = num[LAY] + strrat * den[LAY];
      int JS;
      double FS;
      spec(num[LAY], SPECCOMB, 8., JS, FS);
      const Fac8 f = fac8(LAY, FS);
      const int IND0 = ind0a(LAY, b) + JS, IND1 = ind1a(LAY, b) + JS, INDS = w.indself[LAY];
      for (int IG = 1; IG <= ng; ++IG) {
        C.TAUG(NGS + IG, LAY) = SPECCOMB * sum8(f, absa, 585, IND0, IND1, 9, IG) +
                                self_col[LAY] * w.selffac[LAY] * self(LAY, selfc, INDS, IG);
        C.PFRAC(NGS + IG, LAY) = frac(fracc, ng, IG, JS, FS);
      }
    }
    for (int LAY = LAYTROP + 1; LAY <= NLAYERS; ++LAY) {
      for (int IG = 1; IG <= ng; ++IG) {
        C.TAUG(NGS + IG, LAY) = 0.0;
        C.PFRAC(NGS + IG, LAY) = 0.0;
      }
    }
  }

  void taugb14() {
    constexpr int NGS13 = 134;
    for (int LAY = 1; LAY <= LAYTROP; ++LAY) {
      const int IND0 = ind0a(LAY, 14) + 1, IND1 = ind1a(LAY, 14) + 1, INDS = w.indself[LAY];
      for (int IG = 1; IG <= NG14; ++IG) {
        C.TAUG(NGS13 + IG, LAY) = w.colco2[LAY] * (sum4(LAY, t.absa14, 65, IND0, IND1, IG) +
                                                   w.selffac[LAY] * self(LAY, t.selfrefc14, INDS, IG));
        C.PFRAC(NGS13 + IG, LAY) = t.fracrefac14[IG - 1];
      }
    }
    for (int LAY = LAYTROP + 1; LAY <= NLAYERS; ++LAY) {
      const int IND0 = ind0b(LAY, 14) + 1, IND1 = ind1b(LAY, 14) + 1;
      for (int IG = 1; IG <= NG14; ++IG) {
        C.TAUG(NGS13 + IG, LAY) = w.colco2[LAY] * sum4(LAY, t.absb14, 235, IND0, IND1, IG);
        C.PFRAC(NGS13 + IG, LAY) = t.fracrefbc14[IG - 1];
      }
    }
  }

  void run() {
    taugb1();
    taugb2();
    taugb3();
    taugb4();
    taugb5();
    taugb6();
    taugb7();
    taugb8();
    taugb9();
    taugb10();
    taugb11();
    // TAUGB12: H2O / CO2, TAUGB13: H2O / N2O, TAUGB15: N2O / CO2 (self-continuum times COLH2O), TAUGB16: H2O / CH4
    taugb_lower(122, NG12, 12, t.absa12, t.selfrefc12, t.fracrefac12, 0.009736757, w.colh2o, w.colco2, w.colh2o);
    taugb_lower(130, NG13, 13, t.absa13, t.selfrefc13, t.fracrefac13, 16658.87, w.colh2o, w.coln2o, w.colh2o);
    taugb14();
    taugb_lower(136, NG15, 15, t.absa15, t.selfrefc15, t.fracrefac15, 0.2883201, w.coln2o, w.colco2, w.colh2o);
    taugb_lower(138, NG16, 16, t.absa16, t.selfrefc16, t.fracrefac16, 830.411, w.colh2o, w.colch4, w.colh2o);
  }
};

// GASABS: TAUGB1..16, then the transmittance-table index of the optical depth along the diffusivity angle
void gasabs(const RrtmLwTables& t, Col& C, int LAYTROP, int LAYSWTCH, int LAYLOW) {
  Taugb tg{t, C, C.w, C.nl, LAYTROP, LAYSWTCH, LAYLOW};
  tg.run();
  constexpr double SECANG = 1.66;
  for (int LAY = 1; LAY <= C.nl; ++LAY) {
    for (int IPR = 1; IPR <= NGPT; ++IPR) {
      const double ODEPTH = SECANG * C.TAUG(IPR, LAY);
      double TFF = ODEPTH / (t.bpade + ODEPTH);
      if (ODEPTH <= 0.) TFF = 0.;
      C.ITR(IPR, LAY) = int(5.E3 * TFF + 0.5);
    }
  }
}

// RTRN: the clear and cloudy radiative transfer (one diffusivity angle, the linear-in-tau Planck source), fluxes
// and heating rates on levels 0 (surface) .. NLAYERS
void rtrn(const RrtmLwTables& t, Col& C, double TBOUND, const double* SEMISS) {
  RrtmLwWorkspace& w = C.w;
  const int NLAYERS = C.nl;
  constexpr double SECANG = 1.66, WTNUM = 0.5;
  const double* TAVEL = w.tavel.data();
  const double *PZ = w.pz.data(), *TZ = w.tz.data(), *CLDFRAC = w.cldfrac.data(), *TAUCLOUD = w.taucloud.data();
  const int* ICLDLYR = w.icldlyr.data();
  double *TOTUFLUX = w.totuflux.data(), *TOTDFLUX = w.totdflux.data(), *TOTUCLFL = w.totuclfl.data(),
         *TOTDCLFL = w.totdclfl.data(), *FNET = w.fnet.data(), *FNETC = w.fnetc.data(), *HTR = w.htr.data(),
         *HTRC = w.htrc.data();
  double *TLAYFRAC = w.tlayfrac.data(), *TLEVFRAC = w.tlevfrac.data(), *ABSCLD = w.abscld.data(),
         *ODCLD = w.odcld.data(), *EFCLFRAC = w.efclfrac.data();
  int *INDLAY = w.indlay.data(), *INDLEV = w.indlev.data();
  double *BBU = w.bbu.data() - 1, *BBUTOT = w.bbutot.data() - 1, *ABSS = w.abss.data() - 1,
         *ATOT = w.atot.data() - 1;
  auto ODCLR = [&](int ipr, int lev) -> double& { return w.odclr[std::size_t((lev - 1) * NGPT + ipr - 1)]; };
  auto PLVL = [&](int ib, int lev) -> double& { return w.plvl[std::size_t(lev * NBANDS + ib - 1)]; };
  auto PLAY = [&](int ib, int lev) -> double& { return w.play[std::size_t((lev - 1) * NBANDS + ib - 1)]; };
  auto TOTPLNK = [&](int i, int ib) { return F2(t.totplnk, 181, i, ib); };
  const double* DELWAVE = t.delwave - 1;
  const int* NGB = t.ngb - 1;
  double PLANKBND[NBANDS + 1], PLNKEMIT[NBANDS + 1];
  double RADLU[NGPT + 1], RADLD[NGPT + 1], RADCLRU[NGPT + 1], RADCLRD[NGPT + 1], SEMIS[NGPT + 1],
      RADUEMIT[NGPT + 1], BGLEV[NGPT + 1];

  const int INDBOUND = int(TBOUND - 159.);
  const double TBNDFRAC = TBOUND - double(int(TBOUND));
  for (int LAY = 0; LAY <= NLAYERS; ++LAY) {
    TOTUFLUX[LAY] = 0.0;
    TOTDFLUX[LAY] = 0.0;
    TOTUCLFL[LAY] = 0.0;
    TOTDCLFL[LAY] = 0.0;
    INDLEV[LAY] = int(TZ[LAY] - 159.);
    TLEVFRAC[LAY] = TZ[LAY] - double(int(TZ[LAY]));
  }
  for (int LEV = 1; LEV <= NLAYERS; ++LEV) {
    if (ICLDLYR[LEV] == 1) {
      INDLAY[LEV] = int(TAVEL[LEV] - 159.);
      TLAYFRAC[LEV] = TAVEL[LEV] - double(int(TAVEL[LEV]));
      ODCLD[LEV] = SECANG * TAUCLOUD[LEV];
      const double TRANSCLD = std::exp(-ODCLD[LEV]);
      ABSCLD[LEV] = 1. - TRANSCLD;
      EFCLFRAC[LEV] = ABSCLD[LEV] * CLDFRAC[LEV];
      for (int IPR = 1; IPR <= NGPT; ++IPR) {
        const int IND = C.ITR(IPR, LEV);
        ODCLR(IPR, LEV) = t.tau[IND];
      }
    } else {
      INDLAY[LEV] = int(TAVEL[LEV] - 159.);
      TLAYFRAC[LEV] = TAVEL[LEV] - double(int(TAVEL[LEV]));
    }
  }
  for (int IBAND = 1; IBAND <= NBANDS; ++IBAND) {
    double DBDTLEV = TOTPLNK(INDBOUND + 1, IBAND) - TOTPLNK(INDBOUND, IBAND);
    PLANKBND[IBAND] = DELWAVE[IBAND] * (TOTPLNK(INDBOUND, IBAND) + TBNDFRAC * DBDTLEV);
    DBDTLEV = TOTPLNK(INDLEV[0] + 1, IBAND) - TOTPLNK(INDLEV[0], IBAND);
    PLVL(IBAND, 0) = DELWAVE[IBAND] * (TOTPLNK(INDLEV[0], IBAND) + TLEVFRAC[0] * DBDTLEV);
    PLNKEMIT[IBAND] = SEMISS[IBAND] * PLANKBND[IBAND];
    for (int LEV = 1; LEV <= NLAYERS; ++LEV) {
      DBDTLEV = TOTPLNK(INDLEV[LEV] + 1, IBAND) - TOTPLNK(INDLEV[LEV], IBAND);
      const double DBDTLAY = TOTPLNK(INDLAY[LEV] + 1, IBAND) - TOTPLNK(INDLAY[LEV], IBAND);
      PLAY(IBAND, LEV) = DELWAVE[IBAND] * (TOTPLNK(INDLAY[LEV], IBAND) + TLAYFRAC[LEV] * DBDTLAY);
      PLVL(IBAND, LEV) = DELWAVE[IBAND] * (TOTPLNK(INDLEV[LEV], IBAND) + TLEVFRAC[LEV] * DBDTLEV);
    }
  }
  for (int IPR = 1; IPR <= NGPT; ++IPR) {
    RADCLRD[IPR] = 0.;
    RADLD[IPR] = 0.;
    SEMIS[IPR] = SEMISS[NGB[IPR]];
    RADUEMIT[IPR] = C.PFRAC(IPR, 1) * PLNKEMIT[NGB[IPR]];
    BGLEV[IPR] = C.PFRAC(IPR, NLAYERS) * PLVL(NGB[IPR], NLAYERS);
  }
  // downward radiative transfer
  int ICLDDN = 0;
  for (int LEV = NLAYERS; LEV >= 1; --LEV) {
    double DRAD = 0.0, CLRDRAD = 0.0;
    if (ICLDLYR[LEV] == 1) {
      ICLDDN = 1;
      const int IENT = NGPT * (LEV - 1);
      for (int IPR = 1; IPR <= NGPT; ++IPR) {
        const int INDEX = IENT + IPR;
        const int IND = C.ITR(IPR, LEV);
        const double ODSM = ODCLR(IPR, LEV) + ODCLD[LEV];
        const double FACTOT = ODSM / (t.bpade + ODSM);
        const double BGLAY = C.PFRAC(IPR, LEV) * PLAY(NGB[IPR], LEV);
        const double DELBGUP = BGLEV[IPR] - BGLAY;
        const double TAUF = t.tf[IND];
        BBU[INDEX] = BGLAY + TAUF * DELBGUP;
        BBUTOT[INDEX] = BGLAY + FACTOT * DELBGUP;
        BGLEV[IPR] = C.PFRAC(IPR, LEV) * PLVL(NGB[IPR], LEV - 1);
        const double DELBGDN = BGLEV[IPR] - BGLAY;
        const double BBD = BGLAY + TAUF * DELBGDN;
        const double BBDLEVD = BGLAY + FACTOT * DELBGDN;
        ABSS[INDEX] = 1. - t.trans[IND];
        ATOT[INDEX] = ABSS[INDEX] + ABSCLD[LEV] - ABSS[INDEX] * ABSCLD[LEV];
        const double GASSRC = BBD * ABSS[INDEX];
        RADLD[IPR] = RADLD[IPR] - RADLD[IPR] * (ABSS[INDEX] + EFCLFRAC[LEV] * (1. - ABSS[INDEX])) + GASSRC +
                     CLDFRAC[LEV] * (BBDLEVD * ATOT[INDEX] - GASSRC);
        DRAD = DRAD + RADLD[IPR];
        RADCLRD[IPR] = RADCLRD[IPR] + (BBD - RADCLRD[IPR]) * ABSS[INDEX];
        CLRDRAD = CLRDRAD + RADCLRD[IPR];
      }
    } else {
      const int IENT = NGPT * (LEV - 1);
      for (int IPR = 1; IPR <= NGPT; ++IPR) {
        const int INDEX = IENT + IPR;
        const int IND = C.ITR(IPR, LEV);
        const double BGLAY = C.PFRAC(IPR, LEV) * PLAY(NGB[IPR], LEV);
        const double DELBGUP = BGLEV[IPR] - BGLAY;
        const double TAUF = t.tf[IND];
        BBU[INDEX] = BGLAY + TAUF * DELBGUP;
        BGLEV[IPR] = C.PFRAC(IPR, LEV) * PLVL(NGB[IPR], LEV - 1);
        const double DELBGDN = BGLEV[IPR] - BGLAY;
        const double BBD = BGLAY + TAUF * DELBGDN;
        ABSS[INDEX] = 1. - t.trans[IND];
        RADLD[IPR] = RADLD[IPR] + (BBD - RADLD[IPR]) * ABSS[INDEX];
        DRAD = DRAD + RADLD[IPR];
        if (ICLDDN == 1) {
          RADCLRD[IPR] = RADCLRD[IPR] + (BBD - RADCLRD[IPR]) * ABSS[INDEX];
          CLRDRAD = CLRDRAD + RADCLRD[IPR];
        } else {
          RADCLRD[IPR] = RADLD[IPR];
          CLRDRAD = DRAD;
        }
      }
    }
    TOTDFLUX[LEV - 1] = DRAD * WTNUM;
    TOTDCLFL[LEV - 1] = CLRDRAD * WTNUM;
  }
  // upward radiative transfer
  double URAD = 0.0, CLRURAD = 0.0;
  for (int IPR = 1; IPR <= NGPT; ++IPR) {
    RADLU[IPR] = RADUEMIT[IPR] + (1. - SEMIS[IPR]) * RADLD[IPR];
    URAD = URAD + RADLU[IPR];
    RADCLRU[IPR] = RADUEMIT[IPR] + (1. - SEMIS[IPR]) * RADCLRD[IPR];
    CLRURAD = CLRURAD + RADCLRU[IPR];
  }
  TOTUFLUX[0] = URAD * WTNUM;
  TOTUCLFL[0] = CLRURAD * WTNUM;
  for (int LEV = 1; LEV <= NLAYERS; ++LEV) {
    URAD = 0.0;
    CLRURAD = 0.0;
    const int IENT = NGPT * (LEV - 1);
    if (ICLDLYR[LEV] == 1) {
      for (int IPR = 1; IPR <= NGPT; ++IPR) {
        const int INDEX = IENT + IPR;
        const double GASSRC = BBU[INDEX] * ABSS[INDEX];
        RADLU[IPR] = RADLU[IPR] - RADLU[IPR] * (ABSS[INDEX] + EFCLFRAC[LEV] * (1. - ABSS[INDEX])) + GASSRC +
                     CLDFRAC[LEV] * (BBUTOT[INDEX] * ATOT[INDEX] - GASSRC);
        URAD = URAD + RADLU[IPR];
        RADCLRU[IPR] = RADCLRU[IPR] + (BBU[INDEX] - RADCLRU[IPR]) * ABSS[INDEX];
        CLRURAD = CLRURAD + RADCLRU[IPR];
      }
    } else {
      for (int IPR = 1; IPR <= NGPT; ++IPR) {
        const int INDEX = IENT + IPR;
        RADLU[IPR] = RADLU[IPR] + (BBU[INDEX] - RADLU[IPR]) * ABSS[INDEX];
        URAD = URAD + RADLU[IPR];
        RADCLRU[IPR] = RADCLRU[IPR] + (BBU[INDEX] - RADCLRU[IPR]) * ABSS[INDEX];
        CLRURAD = CLRURAD + RADCLRU[IPR];
      }
    }
    TOTUFLUX[LEV] = URAD * WTNUM;
    TOTUCLFL[LEV] = CLRURAD * WTNUM;
  }
  // fluxes (W/m2) and heating rates (K/day)
  TOTUFLUX[0] = TOTUFLUX[0] * t.fluxfac;
  TOTDFLUX[0] = TOTDFLUX[0] * t.fluxfac;
  FNET[0] = TOTUFLUX[0] - TOTDFLUX[0];
  TOTUCLFL[0] = TOTUCLFL[0] * t.fluxfac;
  TOTDCLFL[0] = TOTDCLFL[0] * t.fluxfac;
  FNETC[0] = TOTUCLFL[0] - TOTDCLFL[0];
  for (int LEV = 1; LEV <= NLAYERS; ++LEV) {
    TOTUFLUX[LEV] = TOTUFLUX[LEV] * t.fluxfac;
    TOTDFLUX[LEV] = TOTDFLUX[LEV] * t.fluxfac;
    FNET[LEV] = TOTUFLUX[LEV] - TOTDFLUX[LEV];
    TOTUCLFL[LEV] = TOTUCLFL[LEV] * t.fluxfac;
    TOTDCLFL[LEV] = TOTDCLFL[LEV] * t.fluxfac;
    FNETC[LEV] = TOTUCLFL[LEV] - TOTDCLFL[LEV];
    const int L = LEV - 1;
    HTR[L] = t.heatfac * (FNET[L] - FNET[LEV]) / (PZ[L] - PZ[LEV]);
    HTRC[L] = t.heatfac * (FNETC[L] - FNETC[LEV]) / (PZ[L] - PZ[LEV]);
  }
  HTR[NLAYERS] = 0.0;
  HTRC[NLAYERS] = 0.0;
}

}  // namespace

int rrtm_lw_nlayers(int nz, double p_top) {
  const int kme = nz + 1;
  return kme + int(std::lround(p_top * 0.01 / deltap)) - 1;
}

RrtmLwWorkspace::RrtmLwWorkspace(int nz_, int nl) : nz(nz_), nlayers(nl) {
  for (auto* v : {&pw1d, &tw1d, &tten1d, &tten1dc, &cldfra1d, &dz1d, &p1d, &t1d, &qv1d, &qc1d, &qr1d, &qi1d, &qs1d,
                  &o3prof, &cldfrc, &o3, &n2o, &ch4, &clwp, &ciwp, &plwp, &piwp})
    v->assign(std::size_t(nz) + 2, 0.0);
  for (auto* v : {&pavel, &tavel, &cldfrac, &taucloud, &coldry, &colh2o, &colco2, &colo3, &coln2o, &colch4, &colo2,
                  &co2mult, &fac00, &fac01, &fac10, &fac11, &forfac, &selffac, &selffrac, &o3prof2, &pzr, &varint,
                  &tlayfrac, &abscld, &odcld, &efclfrac, &pz, &tz, &totdflux, &totuflux, &htr, &htrc, &totuclfl,
                  &totdclfl, &fnet, &fnetc, &tlevfrac})
    v->assign(std::size_t(nl) + 2, 0.0);
  for (auto* v : {&icldlyr, &jp, &jt, &jt1, &indself, &indlay, &indlev}) v->assign(std::size_t(nl) + 2, 0);
  const std::size_t g = std::size_t(NGPT) * std::size_t(nl);
  for (auto* v : {&taug, &pfrac, &odclr, &bbu, &bbutot, &abss, &atot}) v->assign(g, 0.0);
  itr.assign(g, 0);
  wkl.assign(std::size_t(35) * std::size_t(nl), 0.0);
  wx.assign(std::size_t(MAXXSEC) * std::size_t(nl), 0.0);
  plvl.assign(std::size_t(NBANDS) * (std::size_t(nl) + 1), 0.0);
  play.assign(std::size_t(NBANDS) * std::size_t(nl), 0.0);
}

void rrtm_lw_column(const RrtmLwTables& t, const RrtmLwColumn& c, RrtmLwWorkspace& w, double& glw, double& olr,
                    double* rth, double* rthc, std::string& message) {
  const int nz = w.nz, kts = 1, kte = nz, kms = 1, kme = nz + 1;
  message.clear();
  // RRTMLWRAD: the top-down 1-D arrays (K = 1 the model top); model level NK (1-based) is c.*[NK - 1]
  for (int k = kts; k <= kte + 1; ++k) {
    const int NK = kme - k + kms;
    w.pw1d[k] = c.p8w[NK - 1] / 100.;
    w.tw1d[k] = c.t8w[NK - 1];
  }
  for (int K = kts; K <= kte; ++K) {
    w.qv1d[K] = 0.;
    w.qc1d[K] = 0.;
    w.qr1d[K] = 0.;
    w.qi1d[K] = 0.;
    w.qs1d[K] = 0.;
    w.cldfra1d[K] = 0.;
  }
  for (int K = kts; K <= kte; ++K) {
    const int NK = kme - 1 - K + kms;
    w.qv1d[K] = c.qv[NK - 1];
    w.qv1d[K] = std::max(0., w.qv1d[K]);
  }
  for (int K = kts; K <= kte; ++K) {
    const int NK = kme - 1 - K + kms;
    w.tten1d[K] = 0.;
    w.tten1dc[K] = 0.;
    w.t1d[K] = c.t[NK - 1];
    w.p1d[K] = c.p[NK - 1] / 100.;
    w.dz1d[K] = c.dz[NK - 1];
  }
  if (c.icloud != 0) {
    for (int K = kts; K <= kte; ++K) w.cldfra1d[K] = c.cldfra[kme - 1 - K + kms - 1];
    if (c.qc)
      for (int K = kts; K <= kte; ++K) w.qc1d[K] = std::max(0., c.qc[kme - 1 - K + kms - 1]);
    if (c.qr)
      for (int K = kts; K <= kte; ++K) w.qr1d[K] = std::max(0., c.qr[kme - 1 - K + kms - 1]);
    if (!c.has_qi && !c.warm_rain) {
      for (int K = kts; K <= kte; ++K) {
        if (w.t1d[K] < 273.15) {
          w.qi1d[K] = w.qc1d[K];
          w.qs1d[K] = w.qr1d[K];
          w.qc1d[K] = 0.;
          w.qr1d[K] = 0.;
        }
      }
    }
    // QI1D from QI3D whenever F_QI is passed (moist slot 1, 0, for a scheme without ice)
    for (int K = kts; K <= kte; ++K) w.qi1d[K] = std::max(0., c.qi ? c.qi[kme - 1 - K + kms - 1] : 0.);
    if (c.qs)
      for (int K = kts; K <= kte; ++K) w.qs1d[K] = std::max(0., c.qs[kme - 1 - K + kms - 1]);
  }
  const double EMISS0 = c.emiss, TSFC = c.tsk;

  // RRTM
  Col C{w, w.nlayers};
  inirad(t, w.o3prof.data(), w.pw1d.data(), kts, kte);
  double TBOUND = 0, SEMISS[NBANDS + 1];
  mm5atm(t, C, kte, c.co2, c.n2o, c.ch4, c.r_d, c.g, EMISS0, TSFC, TBOUND, SEMISS, message);
  int LAYTROP, LAYSWTCH, LAYLOW;
  setcoef(t, C, LAYTROP, LAYSWTCH, LAYLOW);
  gasabs(t, C, LAYTROP, LAYSWTCH, LAYLOW);
  for (int k = 1; k <= w.nlayers; ++k) w.icldlyr[k] = w.cldfrac[k] > 0. ? 1 : 0;
  rtrn(t, C, TBOUND, SEMISS);
  for (int k = 1; k <= kte; ++k) {
    w.tten1d[k] = w.htr[kte - k] / 86400.;
    w.tten1dc[k] = w.htrc[kte - k] / 86400.;
  }
  glw = w.totdflux[0];
  olr = w.totuflux[kte];
  // RTHRATEN(i,k,j) = RTHRATEN(i,k,j) + TTEN1D(nk)/pi3d(i,k,j), nk = kme-1-k+kms (RTHRATEN zeroed by the driver)
  for (int k = kts; k <= kte; ++k) {
    const int nk = kme - 1 - k + kms;
    rth[k - 1] = 0. + w.tten1d[nk] / c.pi[k - 1];
    rthc[k - 1] = 0. + w.tten1dc[nk] / c.pi[k - 1];
  }
}

void rrtm_lw_batch(const RrtmLwTables& t, const RrtmLwBatchInput& in, int nlayers, double* rth, double* rthc,
                   double* glw, double* olr) {
  const int ncol = in.ncol, nz = in.nz;
#pragma omp parallel
  {
    RrtmLwWorkspace ws(nz, nlayers);
    std::vector<double> p8w(nz + 1), t8w(nz + 1), p(nz), tt(nz), dz(nz), pi(nz), qv(nz), qc(nz), qr(nz), qi(nz),
        qs(nz), cf(nz), r(nz), rc(nz);
    std::string msg;
#pragma omp for schedule(dynamic, 16)
    for (int col = 0; col < ncol; ++col) {
      auto at = [&](int k) { return std::size_t(k) * std::size_t(ncol) + std::size_t(col); };
      for (int k = 0; k <= nz; ++k) {
        p8w[k] = in.p8w[at(k)];
        t8w[k] = in.t8w[at(k)];
      }
      for (int k = 0; k < nz; ++k) {
        p[k] = in.p[at(k)];
        tt[k] = in.t[at(k)];
        dz[k] = in.dz[at(k)];
        pi[k] = in.pi[at(k)];
        qv[k] = in.qv[at(k)];
        cf[k] = in.cldfra[at(k)];
        if (in.qc) qc[k] = in.qc[at(k)];
        if (in.qr) qr[k] = in.qr[at(k)];
        if (in.qi) qi[k] = in.qi[at(k)];
        if (in.qs) qs[k] = in.qs[at(k)];
      }
      const RrtmLwColumn c{p8w.data(), t8w.data(), p.data(), tt.data(), dz.data(), pi.data(), qv.data(),
                           in.qc ? qc.data() : nullptr, in.qr ? qr.data() : nullptr, in.qi ? qi.data() : nullptr,
                           in.qs ? qs.data() : nullptr, cf.data(), in.has_qi, in.tsk[col], in.emiss[col], in.co2,
                           in.n2o, in.ch4, in.r_d, in.g, in.icloud, in.warm_rain};
      rrtm_lw_column(t, c, ws, glw[col], olr[col], r.data(), rc.data(), msg);
      if (!msg.empty()) {
#pragma omp critical(rrtm_lw_message)
        std::fprintf(stderr, "%s\n", msg.c_str());
      }
      for (int k = 0; k < nz; ++k) {
        rth[at(k)] = r[k];
        rthc[at(k)] = rc[k];
      }
    }
  }
}

}  // namespace cppwrf::phys
