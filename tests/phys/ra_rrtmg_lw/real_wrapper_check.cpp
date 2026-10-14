// Check: WRF's RRTMG_LWRAD inputs (debug WRF dump, first radiation call, big-endian whole
// memory arrays) -> C++ rrtmg_lw_column per column -> compare with WRF's RRTMG_LWRAD outputs.
//   real_wrapper_check <rrtmg_lw_in.bin> <rrtmg_lw_out.bin> <rrtmg_lw_tables.bin> <CAMtr_volume_mixing_ratio>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <memory>
#include <vector>

#include "phys/driver/radiation.hpp"
#include "phys/ra_rrtmg_common/cam_gases.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_tables.hpp"

using namespace cppwrf::phys;

struct BE {
  std::vector<char> b;
  size_t o = 0;
  explicit BE(const char* p) { std::ifstream f(p, std::ios::binary); b.assign(std::istreambuf_iterator<char>(f), {}); }
  int i4() { uint32_t v; std::memcpy(&v, &b[o], 4); o += 4; v = __builtin_bswap32(v); int r; std::memcpy(&r, &v, 4); return r; }
  double f8() { uint64_t v; std::memcpy(&v, &b[o], 8); o += 8; v = __builtin_bswap64(v); double r; std::memcpy(&r, &v, 8); return r; }
  std::vector<double> arr(size_t n) { std::vector<double> a(n); for (auto& x : a) x = f8(); return a; }
};

int main(int argc, char** argv) {
  if (argc < 5) return 2;
  BE in(argv[1]), out(argv[2]);
  const int ims = in.i4(), ime = in.i4(), kms = in.i4(), kme = in.i4(), jms = in.i4(), jme = in.i4();
  const int its = in.i4(), ite = in.i4(), jts = in.i4(), jte = in.i4(), kts = in.i4(), kte = in.i4();
  const double p_top = in.f8(), julian = in.f8();
  const int yr = in.i4(), icloud = in.i4(), cldovrlp = in.i4(), idcor = in.i4(), o3input = in.i4(), ghg = in.i4();
  const size_t ni = ime - ims + 1, nk = kme - kms + 1, nj = jme - jms + 1, n3 = ni * nk * nj, n2 = ni * nj;
  auto p8w = in.arr(n3), p = in.arr(n3), pi = in.arr(n3), dz8w = in.arr(n3), t = in.arr(n3), t8w = in.arr(n3),
       rho = in.arr(n3), qv = in.arr(n3), o3 = in.arr(n3), cldfra = in.arr(n3);
  auto tsk = in.arr(n2), emiss = in.arr(n2), xlat = in.arr(n2), xland = in.arr(n2), xice = in.arr(n2), snow = in.arr(n2);
  // cloud extension of the dump (09 UTC dumps): has_reqc/i/s, mp_physics, qc qr qi qs qg, re_cloud/ice/snow
  int hrc = 0, hri = 0, hrs = 0, mp = 0;
  std::vector<double> qc(n3, 0.), qi(n3, 0.), qs(n3, 0.), rec(n3, 0.), rei(n3, 0.), res(n3, 0.);
  const bool clouds = in.o < in.b.size();
  if (clouds) {
    hrc = in.i4(); hri = in.i4(); hrs = in.i4(); mp = in.i4();
    qc = in.arr(n3); in.arr(n3); qi = in.arr(n3); qs = in.arr(n3); in.arr(n3);
    rec = in.arr(n3); rei = in.arr(n3); res = in.arr(n3);
    std::printf("clouds: has_reqc %d has_reqi %d has_reqs %d mp_physics %d\n", hrc, hri, hrs, mp);
  }
  auto w_rth = out.arr(n3), w_rthc = out.arr(n3);
  auto w_glw = out.arr(n2), w_olr = out.arr(n2), w_lwcf = out.arr(n2), w_lwupt = out.arr(n2), w_lwupb = out.arr(n2),
       w_lwdnb = out.arr(n2);
  auto I3 = [&](int i, int k, int j) { return (size_t(j - jms) * nk + (k - kms)) * ni + (i - ims); };
  auto I2 = [&](int i, int j) { return size_t(j - jms) * ni + (i - ims); };

  auto tables = std::make_unique<RrtmgLwTables>();
  if (!load_rrtmg_lw_tables(argv[3], *tables)) return 1;
  auto ws = std::make_unique<RrtmgLwWorkspace>();
  ws->init_pointers();
  RrtmgLwConfig c;
  c.icloud = icloud; c.cldovrlp = cldovrlp; c.idcor = idcor; c.p_top = p_top; c.yr = yr; c.julian = julian;
  c.g = 9.81; c.r_d = 287.; c.cp = 1004.5; c.o3input = o3input;
  if (ghg == 1) {
    CamGases cg;
    cg.load(argv[4]);
    cg.interpolate(yr, julian, c.co2, c.n2o, c.ch4, c.cfc11, c.cfc12);
  } else {
    rrtmg_default_gases(c);
  }
  std::printf("ghg %d: co2 %.17g n2o %.17g ch4 %.17g cfc11 %.17g cfc12 %.17g\n", ghg, c.co2, c.n2o, c.ch4, c.cfc11, c.cfc12);
  const int nz = kte - kts + 1;
  std::vector<double> cp8w(nz + 1), ct8w(nz + 1), cp(nz), cpi(nz), cdz(nz), ct(nz), crho(nz), cqv(nz), zero(nz, 0.0),
      ccld(nz), co3(nz), rth(nz), rthc(nz), cqc(nz), cqi(nz), cqs(nz), crec(nz), crei(nz), cres(nz), cfc(nz);
  double cf_mx = 0;
  long cf_ndiff = 0, ncloudy = 0;
  double mx[4] = {0}, sc[4] = {0};
  long ncol = 0, nexact = 0;
  for (int j = jts; j <= jte; ++j)
    for (int i = its; i <= ite; ++i) {
      for (int k = 0; k <= nz; ++k) { cp8w[k] = p8w[I3(i, k + kts, j)]; ct8w[k] = t8w[I3(i, k + kts, j)]; }
      for (int k = 0; k < nz; ++k) {
        const size_t q = I3(i, k + kts, j);
        cp[k] = p[q]; cpi[k] = pi[q]; cdz[k] = dz8w[q]; ct[k] = t[q]; crho[k] = rho[q]; cqv[k] = qv[q];
        ccld[k] = cldfra[q]; co3[k] = o3[q];
        cqc[k] = qc[q]; cqi[k] = qi[q]; cqs[k] = qs[q]; crec[k] = rec[q]; crei[k] = rei[q]; cres[k] = res[q];
      }
      if (clouds) {  // C++ cal_cldfra1 vs WRF's CLDFRA (t = t_phy, p = p_hyd)
        driver::cal_cldfra1(nz, cqv.data(), cqc.data(), cqi.data(), cqs.data(), ct.data(), cp.data(), cfc.data());
        bool cl = false;
        for (int k = 0; k < nz; ++k) {
          cf_mx = std::max(cf_mx, std::abs(cfc[k] - ccld[k]));
          cf_ndiff += cfc[k] != ccld[k];
          cl = cl || ccld[k] > 0;
        }
        ncloudy += cl;
      }
      double glw, olr, lwcf, lwupt, lwuptc, lwdnt, lwdntc, lwupb, lwupbc, lwdnb, lwdnbc;
      const size_t q2 = I2(i, j);
      rrtmg_lw_column(nz, c, *tables, tsk[q2], emiss[q2], xland[q2], xice[q2], snow[q2], xlat[q2], cp8w.data(),
                      ct8w.data(), cp.data(), cpi.data(), cdz.data(), ct.data(), crho.data(), cqv.data(), cqc.data(),
                      zero.data(), cqi.data(), cqs.data(), zero.data(), ccld.data(), co3.data(), glw, olr, lwcf,
                      lwupt, lwuptc, lwdnt, lwdntc, lwupb, lwupbc, lwdnb, lwdnbc, rth.data(), rthc.data(), *ws,
                      hrc ? crec.data() : nullptr, hri ? crei.data() : nullptr, hrs ? cres.data() : nullptr);
      bool ex = true;
      auto upd = [&](int v, double w, double x) {
        mx[v] = std::max(mx[v], std::abs(w - x)); sc[v] = std::max(sc[v], std::abs(w)); if (w != x) ex = false;
      };
      upd(0, w_glw[q2], glw); upd(1, w_olr[q2], olr); upd(2, w_lwupb[q2], lwupb);
      for (int k = 0; k < nz; ++k) upd(3, w_rth[I3(i, k + kts, j)], rth[k]);
      ++ncol; nexact += ex;
    }
  const char* nm[4] = {"glw", "olr", "lwupb", "rthratlw"};
  std::printf("%ld columns, %ld bit-identical\n", ncol, nexact);
  if (clouds) std::printf("cal_cldfra1: %ld cloudy columns, %ld values differ, max|d| %.3e\n", ncloudy, cf_ndiff, cf_mx);
  for (int v = 0; v < 4; ++v) std::printf("%-9s max|d| %.3e  rel %.3e\n", nm[v], mx[v], mx[v] / sc[v]);
  return 0;
}
