// SPDX-License-Identifier: GPL-3.0-or-later
#include "phys/ra_rrtm_lw/rrtm_lw_tables.hpp"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <vector>

#include "core/libm.hpp"
#include "phys/ra_rrtm_lw/rrtm_lw_data.hpp"

namespace cppwrf::phys {

namespace {

namespace rd = rrtm_data;
constexpr int MG = 16;

// one record of a big-endian unformatted sequential file with 4-byte markers (WRF: -fconvert=big-endian
// -frecord-marker=4); its first n reals, false when the record is shorter (Fortran READ: ERR=)
bool read_record(std::ifstream& f, std::size_t n, std::vector<double>& out) {
  auto be32 = [](const unsigned char* b) {
    return std::uint32_t(b[0]) << 24 | std::uint32_t(b[1]) << 16 | std::uint32_t(b[2]) << 8 | std::uint32_t(b[3]);
  };
  unsigned char m[4];
  if (!f.read(reinterpret_cast<char*>(m), 4)) return false;
  const std::uint32_t len = be32(m);
  if (len < n * 8) return false;
  std::vector<unsigned char> buf(len);
  if (!f.read(reinterpret_cast<char*>(buf.data()), len) || !f.read(reinterpret_cast<char*>(m), 4)) return false;
  out.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    std::uint64_t u = 0;
    for (int k = 0; k < 8; ++k) u = u << 8 | buf[i * 8 + std::size_t(k)];
    std::memcpy(&out[i], &u, 8);
  }
  return true;
}

// CMBGBn, band b (1..16): the band's NGC(b) combined g-points take NGN(NGS(b-1)+IGC) of its 16 g-points each
// (NGN(IGC) for band 1), weighted by RWGT(IPRSM + 16*(b-1)). Every CMBGB has this form; the outer loops only order
// the independent sums.
struct Combine {
  const double* rwgt;
  int ngc(int b) const { return rd::NGC[b - 1]; }
  int ngn(int b, int igc) const { return rd::NGN[(b == 1 ? 0 : rd::NGS[b - 2]) + igc - 1]; }
  // src(rows, MG) -> out(rows, NGC(b)), weighted: abscoefL / abscoefH (the leading indices flattened), SELFREF(10,MG),
  // FORREF(MG) and the other 1-D coefficients (rows = 1)
  void weighted(int b, const double* src, int rows, double* out) const {
    for (int r = 0; r < rows; ++r) {
      int iprsm = 0;
      for (int igc = 1; igc <= ngc(b); ++igc) {
        double sumk = 0.;
        for (int ipr = 1; ipr <= ngn(b, igc); ++ipr) {
          iprsm = iprsm + 1;
          sumk = sumk + src[r + rows * (iprsm - 1)] * rwgt[iprsm + 16 * (b - 1) - 1];
        }
        out[r + rows * (igc - 1)] = sumk;
      }
    }
  }
  // CMBGB9's ABSN2O(3*MG) -> ABSN2OC(3*NG9): JND = (JN-1)*16, JNDC = (JN-1)*NGC(9)
  void weighted_gfirst(int b, const double* src, int cols, double* out) const {
    for (int jn = 1; jn <= cols; ++jn) {
      int iprsm = 0;
      for (int igc = 1; igc <= ngc(b); ++igc) {
        double sumk = 0.;
        for (int ipr = 1; ipr <= ngn(b, igc); ++ipr) {
          iprsm = iprsm + 1;
          const int jnd = (jn - 1) * 16;
          sumk = sumk + src[jnd + iprsm - 1] * rwgt[iprsm + 16 * (b - 1) - 1];
        }
        const int jndc = (jn - 1) * ngc(b);
        out[jndc + igc - 1] = sumk;
      }
    }
  }
  // FRACREFA(MG, cols) -> FRACREFAC(NGC(b), cols): plain sums
  void fractions(int b, const double* src, int cols, double* out) const {
    for (int c = 0; c < cols; ++c) {
      int iprsm = 0;
      for (int igc = 1; igc <= ngc(b); ++igc) {
        double sumf = 0.;
        for (int ipr = 1; ipr <= ngn(b, igc); ++ipr) {
          iprsm = iprsm + 1;
          sumf = sumf + src[iprsm - 1 + MG * c];
        }
        out[igc - 1 + ngc(b) * c] = sumf;
      }
    }
  }
};

template <std::size_t N>
void copy(double (&dst)[N], const double (&src)[N]) {
  for (std::size_t i = 0; i < N; ++i) dst[i] = src[i];
}
template <std::size_t N>
void copy(int (&dst)[N], const int (&src)[N]) {
  for (std::size_t i = 0; i < N; ++i) dst[i] = src[i];
}

}  // namespace

bool load_rrtm_lw_tables(const std::string& path, RrtmLwTables& t, std::string& error) {
  std::ifstream f(path, std::ios::binary);
  if (!f) {
    error = "module_ra_rrtm: error opening RRTM_DATA (" + path + ")";
    return false;
  }
  // the 16 READs: abscoefLn (rows of the low-pressure table), abscoefHn, SELFREFn (10)
  constexpr int rowsL[16] = {65, 65, 650, 585, 585, 65, 585, 35, 715, 65, 65, 585, 585, 65, 585, 585};
  constexpr int rowsH[16] = {235, 235, 1175, 1410, 1175, 0, 235, 265, 235, 235, 235, 0, 0, 235, 0, 0};
  constexpr bool self[16] = {true, true, true, true, true, true, true, true, true, false, true, true, true, true,
                             true, true};
  std::vector<double> L[16], H[16], S[16];
  for (int b = 0; b < 16; ++b) {
    const std::size_t nl = std::size_t(rowsL[b]) * MG, nh = std::size_t(rowsH[b]) * MG, ns = self[b] ? 10 * MG : 0;
    std::vector<double> rec;
    if (!read_record(f, nl + nh + ns, rec)) {
      error = "module_ra_rrtm: error reading RRTM_DATA (" + path + ")";
      return false;
    }
    L[b].assign(rec.begin(), rec.begin() + long(nl));
    H[b].assign(rec.begin() + long(nl), rec.begin() + long(nl + nh));
    S[b].assign(rec.begin() + long(nl + nh), rec.end());
  }

  // rrtm_lookuptable: the Pade approximation tables
  t.tau[0] = 0.0;
  t.tau[5000] = 1.E10;
  t.trans[0] = 1.0;
  t.trans[5000] = 0.0;
  t.tf[0] = 0.0;
  t.tf[5000] = 1.0;
  t.bpade = 1. / 0.278;
  for (int itre = 1; itre <= 4999; ++itre) {
    const double tfn = double(itre) / 5.E3;
    t.tau[itre] = t.bpade * tfn / (1. - tfn);
    t.trans[itre] = libm::exp(-t.tau[itre]);
    if (t.tau[itre] < 0.1)
      t.tf[itre] = t.tau[itre] / 6.;
    else
      t.tf[itre] = 1. - 2. * ((1. / t.tau[itre]) - (t.trans[itre] / (1. - t.trans[itre])));
  }
  t.corr1[0] = 1.;
  t.corr1[200] = 1.;
  t.corr2[0] = 1.;
  t.corr2[200] = 1.;
  for (int i = 1; i <= 199; ++i) {
    const double fp = 0.005 * double(i);
    const double rtfp = std::sqrt(fp);
    t.corr1[i] = rtfp / fp;
    t.corr2[i] = (1. - rtfp) / (1. - fp);
  }
  // the g-point weights RWGT
  double rwgt[MG * 16];
  int igcsm = 0;
  for (int ibnd = 1; ibnd <= 16; ++ibnd) {
    int iprsm = 0;
    if (rd::NGC[ibnd - 1] < 16) {
      double wtsm[MG];
      for (int igc = 1; igc <= rd::NGC[ibnd - 1]; ++igc) {
        igcsm = igcsm + 1;
        double wtsum = 0.;
        for (int ipr = 1; ipr <= rd::NGN[igcsm - 1]; ++ipr) {
          iprsm = iprsm + 1;
          wtsum = wtsum + rd::WT[iprsm - 1];
        }
        wtsm[igc - 1] = wtsum;
      }
      for (int ig = 1; ig <= rd::NG[ibnd - 1]; ++ig) {
        const int ind = (ibnd - 1) * 16 + ig;
        rwgt[ind - 1] = rd::WT[ig - 1] / wtsm[rd::NGM[ind - 1] - 1];
      }
    } else {
      for (int ig = 1; ig <= rd::NG[ibnd - 1]; ++ig) {
        igcsm = igcsm + 1;
        const int ind = (ibnd - 1) * 16 + ig;
        rwgt[ind - 1] = 1.0;
      }
    }
  }

  // CMBGB1..16
  const Combine c{rwgt};
  auto lhs = [&](int b, double* absa, double* absb, double* selfc) {
    c.weighted(b, L[b - 1].data(), rowsL[b - 1], absa);
    if (absb) c.weighted(b, H[b - 1].data(), rowsH[b - 1], absb);
    if (selfc) c.weighted(b, S[b - 1].data(), 10, selfc);
  };
  lhs(1, t.absa1, t.absb1, t.selfrefc1);
  c.weighted(1, rd::FORREF1, 1, t.forrefc1);
  c.fractions(1, rd::FRACREFA1, 1, t.fracrefac1);
  c.fractions(1, rd::FRACREFB1, 1, t.fracrefbc1);
  lhs(2, t.absa2, t.absb2, t.selfrefc2);
  c.fractions(2, rd::FRACREFA2, 13, t.fracrefac2);
  c.weighted(2, rd::FORREF2, 1, t.forrefc2);
  c.fractions(2, rd::FRACREFB2, 1, t.fracrefbc2);
  lhs(3, t.absa3, t.absb3, t.selfrefc3);
  c.fractions(3, rd::FRACREFA3, 10, t.fracrefac3);
  c.fractions(3, rd::FRACREFB3, 5, t.fracrefbc3);
  c.weighted(3, rd::FORREF3, 1, t.forrefc3);
  c.weighted(3, rd::ABSN2OA3, 1, t.absn2oac3);
  c.weighted(3, rd::ABSN2OB3, 1, t.absn2obc3);
  lhs(4, t.absa4, t.absb4, t.selfrefc4);
  c.fractions(4, rd::FRACREFA4, 9, t.fracrefac4);
  c.fractions(4, rd::FRACREFB4, 6, t.fracrefbc4);
  lhs(5, t.absa5, t.absb5, t.selfrefc5);
  c.fractions(5, rd::FRACREFA5, 9, t.fracrefac5);
  c.fractions(5, rd::FRACREFB5, 5, t.fracrefbc5);
  c.weighted(5, rd::CCL45, 1, t.ccl4c5);
  lhs(6, t.absa6, nullptr, t.selfrefc6);
  c.fractions(6, rd::FRACREFA6, 1, t.fracrefac6);
  c.weighted(6, rd::ABSCO26, 1, t.absco2c6);
  c.weighted(6, rd::CFC11ADJ6, 1, t.cfc11adjc6);
  c.weighted(6, rd::CFC126, 1, t.cfc12c6);
  lhs(7, t.absa7, t.absb7, t.selfrefc7);
  c.fractions(7, rd::FRACREFA7, 9, t.fracrefac7);
  c.fractions(7, rd::FRACREFB7, 1, t.fracrefbc7);
  c.weighted(7, rd::ABSCO27, 1, t.absco2c7);
  lhs(8, t.absa8, t.absb8, t.selfrefc8);
  c.fractions(8, rd::FRACREFA8, 1, t.fracrefac8);
  c.fractions(8, rd::FRACREFB8, 1, t.fracrefbc8);
  c.weighted(8, rd::ABSCO2A8, 1, t.absco2ac8);
  c.weighted(8, rd::ABSCO2B8, 1, t.absco2bc8);
  c.weighted(8, rd::ABSN2OA8, 1, t.absn2oac8);
  c.weighted(8, rd::ABSN2OB8, 1, t.absn2obc8);
  c.weighted(8, rd::CFC128, 1, t.cfc12c8);
  c.weighted(8, rd::CFC22ADJ8, 1, t.cfc22adjc8);
  lhs(9, t.absa9, t.absb9, t.selfrefc9);
  c.weighted_gfirst(9, rd::ABSN2O9, 3, t.absn2oc9);
  c.fractions(9, rd::FRACREFA9, 9, t.fracrefac9);
  c.fractions(9, rd::FRACREFB9, 1, t.fracrefbc9);
  lhs(10, t.absa10, t.absb10, nullptr);
  c.fractions(10, rd::FRACREFA10, 1, t.fracrefac10);
  c.fractions(10, rd::FRACREFB10, 1, t.fracrefbc10);
  lhs(11, t.absa11, t.absb11, t.selfrefc11);
  c.fractions(11, rd::FRACREFA11, 1, t.fracrefac11);
  c.fractions(11, rd::FRACREFB11, 1, t.fracrefbc11);
  lhs(12, t.absa12, nullptr, t.selfrefc12);
  c.fractions(12, rd::FRACREFA12, 9, t.fracrefac12);
  lhs(13, t.absa13, nullptr, t.selfrefc13);
  c.fractions(13, rd::FRACREFA13, 9, t.fracrefac13);
  lhs(14, t.absa14, t.absb14, t.selfrefc14);
  c.fractions(14, rd::FRACREFA14, 1, t.fracrefac14);
  c.fractions(14, rd::FRACREFB14, 1, t.fracrefbc14);
  lhs(15, t.absa15, nullptr, t.selfrefc15);
  c.fractions(15, rd::FRACREFA15, 9, t.fracrefac15);
  lhs(16, t.absa16, nullptr, t.selfrefc16);
  c.fractions(16, rd::FRACREFA16, 9, t.fracrefac16);

  // rrtminit: FLUXFAC = PI * 2.D4 with PI = 2.*ASIN(1.); the DATA the column routines read
  const double pi = 2. * libm::asin(1.);
  t.fluxfac = pi * 2.E4;
  t.heatfac = rd::HEATFAC;
  copy(t.nspa, rd::NSPA);
  copy(t.nspb, rd::NSPB);
  copy(t.ngb, rd::NGB);
  copy(t.ixindx, rd::IXINDX);
  copy(t.preflog, rd::PREFLOG);
  copy(t.tref, rd::TREF);
  copy(t.totplnk, rd::TOTPLNK);
  copy(t.delwave, rd::DELWAVE);
  copy(t.refparam2, rd::TAUGB2_REFPARAM);
  copy(t.etaref3, rd::TAUGB3_ETAREF);
  copy(t.h2oref3, rd::TAUGB3_H2OREF);
  copy(t.n2oref3, rd::TAUGB3_N2OREF);
  copy(t.co2ref3, rd::TAUGB3_CO2REF);
  copy(t.h2oref8, rd::TAUGB8_H2OREF);
  copy(t.n2oref8, rd::TAUGB8_N2OREF);
  copy(t.o3ref8, rd::TAUGB8_O3REF);
  copy(t.n2oref9, rd::TAUGB9_N2OREF);
  copy(t.h2oref9, rd::TAUGB9_H2OREF);
  copy(t.ch4ref9, rd::TAUGB9_CH4REF);
  copy(t.etaref9, rd::TAUGB9_ETAREF);
  copy(t.o3sum, rd::O3DATA_O3SUM);
  copy(t.ppsum, rd::O3DATA_PPSUM);
  copy(t.o3win, rd::O3DATA_O3WIN);
  copy(t.ppwin, rd::O3DATA_PPWIN);
  copy(t.pprof, rd::MM5ATM_PPROF);
  copy(t.tprof, rd::MM5ATM_TPROF);
  return true;
}

}  // namespace cppwrf::phys
