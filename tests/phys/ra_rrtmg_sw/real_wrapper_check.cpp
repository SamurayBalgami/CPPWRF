// Check: WRF's RRTMG_SWRAD inputs (debug WRF dump, big-endian whole memory arrays) ->
// C++ rrtmg_sw_column per daytime column -> compare with WRF's RRTMG_SWRAD outputs.
//   real_wrapper_check <rrtmg_sw_in.bin> <rrtmg_sw_out.bin> <rrtmg_sw_tables.bin> <CAMtr_volume_mixing_ratio>
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
#include "phys/ra_rrtmg_sw/rrtmg_sw_tables.hpp"

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
  auto ws = std::make_unique<RrtmgSwWorkspace>();
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
  const int nz = kte - kts + 1;
  std::vector<double> cp8w(nz + 1), ct8w(nz + 1), cp(nz), cpi(nz), cdz(nz), ct(nz), crho(nz), cqv(nz), zero(nz, 0.0),
      ccld(nz), co3(nz), rth(nz), rthc(nz), cqc(nz), cqi(nz), cqs(nz), crec(nz), crei(nz), cres(nz);
  double mx[4] = {0}, sc[4] = {0};
  long nday = 0, nexact = 0, nnight = 0;
  for (int j = jts; j <= jte; ++j)
    for (int i = its; i <= ite; ++i) {
      const size_t q2 = I2(i, j);
      if (coszen[q2] <= 0.0) { ++nnight; continue; }
      for (int k = 0; k <= nz; ++k) { cp8w[k] = p8w[I3(i, k + kts, j)]; ct8w[k] = t8w[I3(i, k + kts, j)]; }
      for (int k = 0; k < nz; ++k) {
        const size_t q = I3(i, k + kts, j);
        cp[k] = p[q]; cpi[k] = pi[q]; cdz[k] = dz8w[q]; ct[k] = t[q]; crho[k] = rho[q]; cqv[k] = qv[q];
        ccld[k] = cldfra[q]; co3[k] = o3[q];
        cqc[k] = qc[q]; cqi[k] = qi[q]; cqs[k] = qs[q]; crec[k] = rec[q]; crei[k] = rei[q]; cres[k] = res[q];
      }
      SwColumnOut o{};
      rrtmg_sw_column(nz, c, *tables, coszen[q2], albedo[q2], xland[q2], xice[q2], snow[q2], xlat[q2], cp8w.data(),
                      ct8w.data(), cp.data(), cpi.data(), cdz.data(), ct.data(), crho.data(), cqv.data(), cqc.data(),
                      zero.data(), cqi.data(), cqs.data(), zero.data(), ccld.data(), co3.data(), o, rth.data(),
                      rthc.data(), *ws, hrc ? crec.data() : nullptr, hri ? crei.data() : nullptr,
                      hrs ? cres.data() : nullptr);
      bool ex = true;
      auto upd = [&](int v, double w, double x) {
        mx[v] = std::max(mx[v], std::abs(w - x)); sc[v] = std::max(sc[v], std::abs(w)); if (w != x) ex = false;
      };
      upd(0, w_gsw[q2], o.gsw); upd(1, w_swdnb[q2], o.swdnb); upd(2, w_swupt[q2], o.swupt);
      for (int k = 0; k < nz; ++k) upd(3, w_rth[I3(i, k + kts, j)], rth[k]);
      ++nday; nexact += ex;
    }
  const char* nm[4] = {"gsw", "swdnb", "swupt", "rthratsw"};
  std::printf("%ld daytime columns (%ld night), %ld bit-identical\n", nday, nnight, nexact);
  for (int v = 0; v < 4; ++v) std::printf("%-9s max|d| %.3e  rel %.3e\n", nm[v], mx[v], mx[v] / std::max(sc[v], 1e-300));
  return 0;
}
