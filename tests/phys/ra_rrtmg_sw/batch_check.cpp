// Check of the batch (device) RRTMG SW: WRF's RRTMG_SWRAD inputs (debug WRF dump) ->
// RrtmgSwBatch over all columns -> compare with WRF's RRTMG_SWRAD outputs.
//   batch_check <rrtmg_sw_in.bin> <rrtmg_sw_out.bin> <rrtmg_sw_tables.bin> <CAMtr_volume_mixing_ratio> [chunk]
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <memory>
#include <vector>

#include "phys/ra_rrtmg_common/cam_gases.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_batch.hpp"
#include "core/backend.hpp"
#include <chrono>
#include <cstdlib>
#include "phys/ra_rrtmg_sw/rrtmg_sw_tables.hpp"

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
  const int chunk = argc > 5 ? std::atoi(argv[5]) : 2048;
  BE in(argv[1]), out(argv[2]);
  const int ims = in.i4(), ime = in.i4(), kms = in.i4(), kme = in.i4(), jms = in.i4(), jme = in.i4();
  const int its = in.i4(), ite = in.i4(), jts = in.i4(), jte = in.i4(), kts = in.i4(), kte = in.i4();
  const double xtime = in.f8(), gmt = in.f8(), radt = in.f8(), degrad = in.f8(), declin = in.f8(), solcon = in.f8(),
               julian = in.f8();
  const int julday = in.i4(), yr = in.i4(), icloud = in.i4(), cldovrlp = in.i4(), idcor = in.i4(), o3input = in.i4(),
            ghg = in.i4(), aer_opt = in.i4();
  std::printf("xtime %g gmt %g radt %g declin %.17g solcon %.17g julian %.17g julday %d yr %d aer_opt %d\n", xtime, gmt,
              radt, declin, solcon, julian, julday, yr, aer_opt);
  (void)degrad;
  const size_t ni = ime - ims + 1, nk = kme - kms + 1, nj = jme - jms + 1, n3 = ni * nk * nj, n2 = ni * nj;
  auto p8w = in.arr(n3), p = in.arr(n3), pi = in.arr(n3), dz8w = in.arr(n3), t = in.arr(n3), t8w = in.arr(n3),
       rho = in.arr(n3), qv = in.arr(n3), o3 = in.arr(n3), cldfra = in.arr(n3);
  auto xlat = in.arr(n2), xlong = in.arr(n2), albedo = in.arr(n2), tsk = in.arr(n2), xland = in.arr(n2),
       xice = in.arr(n2), snow = in.arr(n2), coszen = in.arr(n2);
  int hrc = 0, hri = 0, hrs = 0, mp = 0;  // cloud extension of the dump (09 UTC dumps)
  std::vector<double> qc(n3, 0.), qi(n3, 0.), qs(n3, 0.), rec(n3, 0.), rei(n3, 0.), res(n3, 0.);
  if (in.o < in.b.size()) {
    hrc = in.i4(); hri = in.i4(); hrs = in.i4(); mp = in.i4();
    qc = in.arr(n3); in.arr(n3); qi = in.arr(n3); qs = in.arr(n3); in.arr(n3);
    rec = in.arr(n3); rei = in.arr(n3); res = in.arr(n3);
    std::printf("clouds: has_reqc %d has_reqi %d has_reqs %d mp_physics %d\n", hrc, hri, hrs, mp);
  }
  auto w_rth = out.arr(n3), w_rthc = out.arr(n3);
  auto w_gsw = out.arr(n2), w_swdnb = out.arr(n2), w_swupb = out.arr(n2), w_swdnt = out.arr(n2),
       w_swupt = out.arr(n2);
  auto I3 = [&](int i, int k, int j) { return (size_t(j - jms) * nk + (k - kms)) * ni + (i - ims); };
  auto I2 = [&](int i, int j) { return size_t(j - jms) * ni + (i - ims); };

  auto tables = std::make_unique<RrtmgSwTables>();
  if (!load_rrtmg_sw_tables(argv[3], *tables)) { std::fprintf(stderr, "tables\n"); return 1; }
  RrtmgSwConfig c;
  c.icloud = icloud; c.cldovrlp = cldovrlp; c.idcor = idcor; c.p_top = 5000.; c.yr = yr; c.julian = julian;
  c.solcon = solcon; c.g = 9.81; c.r_d = 287.; c.cp = 1004.5; c.o3input = o3input;
  if (ghg == 1) {
    CamGases cg;
    cg.load(argv[4]);
    double c11, c12;
    cg.interpolate(yr, julian, c.co2, c.n2o, c.ch4, c11, c12);
  } else {
    c.co2 = (280. + 90. * std::exp(0.02 * (yr - 2000))) * 1.e-6; c.ch4 = 1774.e-9; c.n2o = 319.e-9;
  }
  const int nz = kte - kts + 1, nx = ite - its + 1, ny = jte - jts + 1, ncol = nx * ny;
  const bool clouds = hrc || hri || hrs || mp;
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
  SwBatchInput bi{};
  bi.ncol = ncol;
  bi.coszen = dev(col2(coszen)); bi.albedo = dev(col2(albedo)); bi.xland = dev(col2(xland));
  bi.xice = dev(col2(xice)); bi.snow = dev(col2(snow));
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
  double* d_gsw = static_cast<double*>(dev_alloc(sizeof(double) * ncol));
  const RrtmgSwCallConfig cc = rrtmg_sw_call_config(nz, c);
  RrtmgSwBatch batch(*tables, cc.nlayers, chunk);
  for (int rep = 0; rep < 2; ++rep) {
    const auto t0 = std::chrono::steady_clock::now();
    batch.run(cc, bi, SwBatchOutput{d_rth, d_gsw});
    std::vector<double> tmp(1);
    copy_to_host(tmp.data(), d_gsw, sizeof(double));
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("batch run %d: %.3f s (%d columns, chunk %d)\n", rep, sec, ncol, chunk);
  }
  std::vector<double> rth(size_t(nz) * ncol), gsw(ncol);
  copy_to_host(rth.data(), d_rth, rth.size() * sizeof(double));
  copy_to_host(gsw.data(), d_gsw, gsw.size() * sizeof(double));
  double mx[2] = {0}, sc[2] = {0};
  long nexact = 0, nday = 0;
  for (int j = jts; j <= jte; ++j)
    for (int i = its; i <= ite; ++i) {
      if (!(coszen[I2(i, j)] > 0.0)) continue;
      ++nday;
      const int col = (i - its) + (j - jts) * nx;
      bool ex = true;
      auto upd = [&](int v, double w, double x) {
        mx[v] = std::max(mx[v], std::abs(w - x)); sc[v] = std::max(sc[v], std::abs(w)); if (w != x) ex = false;
      };
      upd(0, w_gsw[I2(i, j)], gsw[col]);
      for (int k = 0; k < nz; ++k) upd(1, w_rth[I3(i, k + kts, j)], rth[size_t(k) * ncol + col]);
      nexact += ex;
    }
  const char* nm[2] = {"gsw", "rthratsw"};
  std::printf("%ld daytime columns, %ld bit-identical to WRF\n", nday, nexact);
  for (int v = 0; v < 2; ++v) std::printf("%-9s max|d| %.3e  rel %.3e\n", nm[v], mx[v], mx[v] / std::max(sc[v], 1e-300));
  for (void* d : allocs) dev_free(d);
  dev_free(d_rth); dev_free(d_gsw);
  return 0;
}
