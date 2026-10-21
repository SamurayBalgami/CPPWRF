// CPPWRF side of the MYNN-EDMF differential test (see driver_bl_mynn.f90): mynnedmf_driver's column work
// (flux limits, mynnedmf_pre_run, mynnedmf, mynnedmf_post_run) on ni columns for nsteps calls, the same records out.
#include <cstdio>
#include <memory>
#include <vector>

#include "phys/bl_mynn/mynnedmf.hpp"

using namespace cppwrf;
using namespace cppwrf::phys::mynnedmf;

int main(int argc, char** argv) {
  if (argc < 3) { std::fprintf(stderr, "usage: %s IN.bin OUT.bin\n", argv[0]); return 2; }
  std::FILE* f = std::fopen(argv[1], "rb");
  if (!f) return 1;
  int hdr[3];
  double rh[2];
  if (std::fread(hdr, sizeof(int), 3, f) != 3 || std::fread(rh, sizeof(double), 2, f) != 2) return 1;
  const int ni = hdr[0], nz = hdr[1], nsteps = hdr[2];
  const double delt = rh[0], dx = rh[1];
  constexpr int n3 = 14, n2 = 12, no3 = 25, no2 = 7;
  constexpr int K = 128;
  if (nz > K) return 3;
  std::vector<double> in3(size_t(n3) * nz * ni), in2(size_t(n2) * ni), o3(size_t(no3) * nz * ni, 0.0),
      o2(size_t(no2) * ni, 0.0);
  if (std::fread(in3.data(), sizeof(double), in3.size(), f) != in3.size()) return 1;
  if (std::fread(in2.data(), sizeof(double), in2.size(), f) != in2.size()) return 1;
  std::fclose(f);
  // Fortran (i, k) order: element (n, i, k) at (n * nz + (k - 1)) * ni + i
  auto I3 = [&](int n, int i, int k) -> double& { return in3[(size_t(n) * nz + (k - 1)) * ni + i]; };
  auto I2 = [&](int n, int i) -> double& { return in2[size_t(n) * ni + i]; };
  auto O3 = [&](int n, int i, int k) -> double& { return o3[(size_t(n) * nz + (k - 1)) * ni + i]; };
  auto O2 = [&](int n, int i) -> double& { return o2[size_t(n) * ni + i]; };
  auto col = std::make_unique<Column<K>>();
  Column<K>& c = *col;
  c.kte = nz;
  for (int s = 1; s <= nsteps; ++s) {
    const int initflag = (s == 1) ? 1 : 0;
    for (int i = 0; i < ni; ++i) {
      double qv[K + 2], qc[K + 2], qi[K + 2], qs[K + 2];
      for (int k = 1; k <= nz; ++k) {
        c.dz[k] = I3(0, i, k); c.u[k] = I3(1, i, k); c.v[k] = I3(2, i, k); c.w[k] = I3(3, i, k);
        c.th[k] = I3(4, i, k); c.tk[k] = I3(5, i, k); c.pres[k] = I3(6, i, k); c.ex[k] = I3(7, i, k);
        c.rho[k] = I3(8, i, k);
        qv[k] = I3(9, i, k); qc[k] = I3(10, i, k); qi[k] = I3(11, i, k); qs[k] = I3(12, i, k);
      }
      double hfx1 = I2(6, i);
      if (hfx1 > 1200.) hfx1 = 1200.;
      if (hfx1 < -600.) hfx1 = -600.;
      double qfx1 = I2(7, i);
      if (qfx1 > 9e-4) qfx1 = 9e-4;
      if (qfx1 < -3e-4) qfx1 = -3e-4;
      if (initflag == 0) {
        for (int k = 1; k <= nz; ++k) {
          c.qc_bl[k] = O3(16, i, k); c.qi_bl[k] = O3(17, i, k); c.cldfra_bl[k] = O3(18, i, k);
          c.qke[k] = O3(0, i, k); c.qsq[k] = O3(7, i, k); c.tsq[k] = O3(6, i, k); c.cov[k] = O3(8, i, k);
          c.sh[k] = O3(2, i, k); c.sm[k] = O3(3, i, k); c.kh[k] = O3(4, i, k); c.km[k] = O3(5, i, k);
          c.el[k] = O3(1, i, k); c.qke_adv[k] = O3(0, i, k);
        }
        c.pblh = O2(0, i);
        c.kpbl = int(O2(6, i));
      }
      for (int k = 1; k <= nz; ++k) {  // mynnedmf_pre_run
        c.sqv[k] = qv[k] / (one + qv[k]);
        c.sqc[k] = qc[k] / (one + qv[k]);
        c.sqi[k] = qi[k] / (one + qv[k]);
        c.sqs[k] = qs[k] / (one + qv[k]);
      }
      mynnedmf<K>(c, initflag, delt, dx, I2(0, i), I2(1, i), I2(3, i), I2(4, i), hfx1, qfx1, I2(8, i), I2(10, i),
                  I2(11, i));
      for (int k = 1; k <= nz; ++k) {  // mynnedmf_post_run
        const double sq = qv[k] / (one + qv[k]);
        const double sqv = sq + c.dqv[k] * delt;
        const double rq = sqv / (one - sqv);
        c.dqv[k] = (rq - qv[k]) / delt;
        double sqx = qc[k] / (one + qv[k]);
        double sqn = sqx + c.dqc[k] * delt;
        c.dqc[k] = (sqn * (one + sqv) - qc[k]) / delt;
        sqx = qi[k] / (one + qv[k]);
        sqn = sqx + c.dqi[k] * delt;
        c.dqi[k] = (sqn * (one + sqv) - qi[k]) / delt;
      }
      for (int k = 1; k <= nz; ++k) {
        O3(0, i, k) = c.qke[k]; O3(1, i, k) = c.el[k]; O3(2, i, k) = c.sh[k]; O3(3, i, k) = c.sm[k];
        O3(4, i, k) = c.kh[k]; O3(5, i, k) = c.km[k]; O3(6, i, k) = c.tsq[k]; O3(7, i, k) = c.qsq[k];
        O3(8, i, k) = c.cov[k]; O3(9, i, k) = c.qke_adv[k];
        O3(10, i, k) = c.du[k]; O3(11, i, k) = c.dv[k]; O3(12, i, k) = c.dth[k]; O3(13, i, k) = c.dqv[k];
        O3(14, i, k) = c.dqc[k]; O3(15, i, k) = c.dqi[k];
        O3(16, i, k) = c.qc_bl[k] / (one - c.sqv[k]);
        O3(17, i, k) = c.qi_bl[k] / (one - c.sqv[k]);
        O3(18, i, k) = c.cldfra_bl[k];
        O3(19, i, k) = c.edmf_a[k]; O3(20, i, k) = c.edmf_w[k]; O3(21, i, k) = c.edmf_qt[k];
        O3(22, i, k) = c.edmf_thl[k]; O3(23, i, k) = c.edmf_ent[k]; O3(24, i, k) = c.edmf_qc[k];
      }
      O2(0, i) = c.pblh; O2(1, i) = c.maxwidth; O2(2, i) = c.maxmf; O2(3, i) = c.ztop_plume;
      O2(4, i) = c.excess_h; O2(5, i) = c.excess_q; O2(6, i) = double(c.kpbl);
    }
    if (s < nsteps)
      for (int i = 0; i < ni; ++i)
        for (int k = 1; k <= nz; ++k) {
          I3(1, i, k) = I3(1, i, k) + O3(10, i, k) * delt;
          I3(2, i, k) = I3(2, i, k) + O3(11, i, k) * delt;
          I3(4, i, k) = I3(4, i, k) + O3(12, i, k) * delt;
          I3(9, i, k) = c_max(I3(9, i, k) + O3(13, i, k) * delt, 0.0);
          I3(10, i, k) = c_max(I3(10, i, k) + O3(14, i, k) * delt, 0.0);
          I3(11, i, k) = c_max(I3(11, i, k) + O3(15, i, k) * delt, 0.0);
          I3(5, i, k) = I3(4, i, k) * I3(7, i, k);
        }
  }
  std::FILE* g = std::fopen(argv[2], "wb");
  if (!g) return 1;
  std::fwrite(o3.data(), sizeof(double), o3.size(), g);
  std::fwrite(o2.data(), sizeof(double), o2.size(), g);
  std::fclose(g);
  return 0;
}
