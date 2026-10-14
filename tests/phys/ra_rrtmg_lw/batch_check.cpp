// Check of the batch (device) RRTMG LW: WRF's RRTMG_LWRAD inputs (debug WRF dump, big-endian
// whole memory arrays) -> RrtmgLwBatch over all columns -> compare with WRF's RRTMG_LWRAD outputs.
//   batch_check <rrtmg_lw_in.bin> <rrtmg_lw_out.bin> <rrtmg_lw_tables.bin> <CAMtr_volume_mixing_ratio> [chunk]
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
#include "phys/ra_rrtmg_lw/rrtmg_lw_batch.hpp"
#include "core/backend.hpp"
#include <chrono>
#include "phys/ra_rrtmg_lw/rrtmg_lw_tables.hpp"

using namespace cppwrf;
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
  const int chunk = argc > 5 ? std::atoi(argv[5]) : 3072;
  const int ims = in.i4(), ime = in.i4(), kms = in.i4(), kme = in.i4(), jms = in.i4(), jme = in.i4();
  const int its = in.i4(), ite = in.i4(), jts = in.i4(), jte = in.i4(), kts = in.i4(), kte = in.i4();
  const double p_top = in.f8(), julian = in.f8();
  const int yr = in.i4(), icloud = in.i4(), cldovrlp = in.i4(), idcor = in.i4(), o3input = in.i4(), ghg = in.i4();
  (void)icloud; (void)cldovrlp; (void)idcor;
  const size_t ni = ime - ims + 1, nk = kme - kms + 1, nj = jme - jms + 1, n3 = ni * nk * nj, n2 = ni * nj;
  auto p8w = in.arr(n3), p = in.arr(n3), pi = in.arr(n3), dz8w = in.arr(n3), t = in.arr(n3), t8w = in.arr(n3),
       rho = in.arr(n3), qv = in.arr(n3), o3 = in.arr(n3), cldfra = in.arr(n3);
  auto tsk = in.arr(n2), emiss = in.arr(n2), xlat = in.arr(n2), xland = in.arr(n2), xice = in.arr(n2), snow = in.arr(n2);
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
  auto w_glw = out.arr(n2), w_olr = out.arr(n2);
  auto I3 = [&](int i, int k, int j) { return (size_t(j - jms) * nk + (k - kms)) * ni + (i - ims); };
  auto I2 = [&](int i, int j) { return size_t(j - jms) * ni + (i - ims); };

  auto tables = std::make_unique<RrtmgLwTables>();
  if (!load_rrtmg_lw_tables(argv[3], *tables)) return 1;
  RrtmgLwConfig c;
  c.p_top = p_top; c.yr = yr; c.julian = julian; c.g = 9.81; c.cp = 1004.5; c.o3input = o3input;
  if (ghg == 1) {
    CamGases cg;
    cg.load(argv[4]);
    cg.interpolate(yr, julian, c.co2, c.n2o, c.ch4, c.cfc11, c.cfc12);
  } else {
    rrtmg_default_gases(c);
  }
  const int nz = kte - kts + 1, nx = ite - its + 1, ny = jte - jts + 1, ncol = nx * ny;
  // column-fastest host arrays
  auto col3 = [&](const std::vector<double>& a, int nlev) {
    std::vector<double> r(size_t(nlev) * ncol);
    for (int j = jts; j <= jte; ++j)
      for (int i = its; i <= ite; ++i)
        for (int k = 0; k < nlev; ++k) r[size_t(k) * ncol + (i - its) + size_t(j - jts) * nx] = a[I3(i, k + kts, j)];
    return r;
  };
  auto col2 = [&](const std::vector<double>& a) {
    std::vector<double> r(ncol);
    for (int j = jts; j <= jte; ++j)
      for (int i = its; i <= ite; ++i) r[(i - its) + size_t(j - jts) * nx] = a[I2(i, j)];
    return r;
  };
  std::vector<void*> allocs;
  auto dev = [&](const std::vector<double>& h) {
    void* d = dev_alloc(h.size() * sizeof(double));
    copy_to_dev(d, h.data(), h.size() * sizeof(double));
    allocs.push_back(d);
    return static_cast<const double*>(d);
  };
  LwBatchInput bi{};
  bi.ncol = ncol;
  bi.tsk = dev(col2(tsk)); bi.emiss = dev(col2(emiss)); bi.xland = dev(col2(xland)); bi.xice = dev(col2(xice));
  bi.snow = dev(col2(snow));
  bi.p8w = dev(col3(p8w, nz + 1)); bi.t8w = dev(col3(t8w, nz + 1));
  bi.p3d = dev(col3(p, nz)); bi.pi3d = dev(col3(pi, nz)); bi.t3d = dev(col3(t, nz)); bi.qv = dev(col3(qv, nz));
  bi.cldfra = dev(col3(cldfra, nz)); bi.o3 = dev(col3(o3, nz));
  if (clouds) {
    bi.qc = dev(col3(qc, nz)); bi.qi = dev(col3(qi, nz)); bi.qs = dev(col3(qs, nz));
    if (hrc) bi.re_cloud = dev(col3(rec, nz));
    if (hri) bi.re_ice = dev(col3(rei, nz));
    if (hrs) bi.re_snow = dev(col3(res, nz));
  }
  double* d_rth = static_cast<double*>(dev_alloc(sizeof(double) * nz * ncol));
  double* d_glw = static_cast<double*>(dev_alloc(sizeof(double) * ncol));
  double* d_olr = static_cast<double*>(dev_alloc(sizeof(double) * ncol));
  const RrtmgLwCallConfig cc = rrtmg_lw_call_config(nz, c, hrc, hri, hrs);
  RrtmgLwBatch batch(*tables, cc.nlayers, chunk);
  for (int rep = 0; rep < 2; ++rep) {
    const auto t0 = std::chrono::steady_clock::now();
    batch.run(cc, bi, LwBatchOutput{d_rth, d_glw, d_olr});
    std::vector<double> tmp(1);
    copy_to_host(tmp.data(), d_glw, sizeof(double));  // synchronises
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("batch run %d: %.3f s (%d columns, chunk %d)\n", rep, sec, ncol, chunk);
  }
  std::vector<double> rth(size_t(nz) * ncol), glw(ncol), olr(ncol);
  copy_to_host(rth.data(), d_rth, rth.size() * sizeof(double));
  copy_to_host(glw.data(), d_glw, glw.size() * sizeof(double));
  copy_to_host(olr.data(), d_olr, olr.size() * sizeof(double));
  double mx[3] = {0}, sc[3] = {0};
  long nexact = 0;
  for (int j = jts; j <= jte; ++j)
    for (int i = its; i <= ite; ++i) {
      const int col = (i - its) + (j - jts) * nx;
      bool ex = true;
      auto upd = [&](int v, double w, double x) {
        mx[v] = std::max(mx[v], std::abs(w - x)); sc[v] = std::max(sc[v], std::abs(w)); if (w != x) ex = false;
      };
      upd(0, w_glw[I2(i, j)], glw[col]); upd(1, w_olr[I2(i, j)], olr[col]);
      for (int k = 0; k < nz; ++k) upd(2, w_rth[I3(i, k + kts, j)], rth[size_t(k) * ncol + col]);
      nexact += ex;
    }
  const char* nm[3] = {"glw", "olr", "rthratlw"};
  std::printf("%d columns, %ld bit-identical to WRF\n", ncol, nexact);
  for (int v = 0; v < 3; ++v) std::printf("%-9s max|d| %.3e  rel %.3e\n", nm[v], mx[v], mx[v] / sc[v]);
  for (void* d : allocs) dev_free(d);
  dev_free(d_rth); dev_free(d_glw); dev_free(d_olr);
  return 0;
}
