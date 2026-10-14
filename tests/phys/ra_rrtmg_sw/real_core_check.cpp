// Check: WRF's rrtmg_sw core inputs (debug WRF dump, daytime columns) -> C++ rrtmg_sw_core -> compare
//   real_core_check <rrtmg_sw_core.bin> <rrtmg_sw_tables.bin>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <memory>
#include <vector>
#include "phys/ra_rrtmg_sw/rrtmg_sw.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_tables.hpp"
using namespace cppwrf::phys;
struct BE { std::vector<char> b; size_t o = 0;
  explicit BE(const char* p) { std::ifstream f(p, std::ios::binary); b.assign(std::istreambuf_iterator<char>(f), {}); }
  bool more() const { return o < b.size(); }
  int i4() { uint32_t v; std::memcpy(&v, &b[o], 4); o += 4; v = __builtin_bswap32(v); int r; std::memcpy(&r, &v, 4); return r; }
  double f8() { uint64_t v; std::memcpy(&v, &b[o], 8); o += 8; v = __builtin_bswap64(v); double r; std::memcpy(&r, &v, 8); return r; }
  void arr(double* a, int n) { for (int k = 0; k < n; ++k) a[k] = f8(); } };
int main(int, char** argv) {
  BE in(argv[1]);
  auto t = std::make_unique<RrtmgSwTables>(); if (!load_rrtmg_sw_tables(argv[2], *t)) return 1;
  auto ws = std::make_unique<RrtmgSwWorkspace>(); ws->init_pointers();
  const int M = RRTMG_SW_MAX_LAYERS + 1;
  std::vector<double> play(M), plev(M), tlay(M), tlev(M), h2o(M), o3(M), co2(M), ch4(M), n2o(M), o2(M), rei(M, 10.), rel(M, 5.), res(M, 10.);
  std::vector<double> wu(M), wd(M), whr(M), wuc(M), wdc(M), whrc(M), wdir(M), wdif(M);
  std::vector<double> u(M), d(M), hr(M), uc(M), dc(M), hrc(M), dir(M), dif(M), a1(M), a2(M), a3(M), a4(M);
  double mx[4] = {0}, sc[4] = {0}; long n = 0, ne = 0;
  while (in.more()) {
    in.i4(); in.i4(); const int nlay = in.i4(); for (int q = 0; q < 5; ++q) in.i4();
    const double adjes = in.f8(), scon = in.f8(), cz = in.f8();
    in.arr(play.data(), nlay); in.arr(plev.data(), nlay + 1); in.arr(tlay.data(), nlay); in.arr(tlev.data(), nlay + 1); in.f8();
    for (auto* a : {&h2o, &o3, &co2, &ch4, &n2o, &o2}) in.arr(a->data(), nlay);
    const double asdir = in.f8(), asdif = in.f8(), aldir = in.f8(), aldif = in.f8();
    in.arr(wu.data(), nlay + 1); in.arr(wd.data(), nlay + 1); in.arr(whr.data(), nlay);
    in.arr(wuc.data(), nlay + 1); in.arr(wdc.data(), nlay + 1); in.arr(whrc.data(), nlay);
    in.arr(wdir.data(), nlay + 1); in.arr(wdif.data(), nlay + 1);
    for (int i = 0; i < nlay * 112; ++i) { ws->cldfmc[i] = ws->ciwpmc[i] = ws->clwpmc[i] = ws->cswpmc[i] = ws->taucmc[i] = 0; ws->ssacmc[i] = 1; ws->asmcmc[i] = ws->fsfcmc[i] = ws->taormc[i] = 0; }
    if (n == 0) g_rrtmg_sw_dbg = std::fopen("csw_dbg.bin", "wb");
    rrtmg_sw_core(nlay, 8.4349451468392243, play.data(), plev.data(), tlay.data(), tlev.data(), h2o.data(), o3.data(), co2.data(),
                  ch4.data(), n2o.data(), o2.data(), asdir, asdif, aldir, aldif, cz, adjes, scon, rei.data(), rel.data(), res.data(),
                  u.data(), d.data(), hr.data(), uc.data(), dc.data(), hrc.data(), dir.data(), dif.data(), a1.data(), a2.data(), a3.data(), a4.data(), *ws, *t);
    bool ex = true;
    auto upd = [&](int v, const std::vector<double>& w, const std::vector<double>& c, int k) {
      for (int i = 0; i < k; ++i) { mx[v] = std::max(mx[v], std::abs(w[i] - c[i])); sc[v] = std::max(sc[v], std::abs(w[i])); if (w[i] != c[i]) ex = false; } };
    upd(0, wu, u, nlay + 1); upd(1, wd, d, nlay + 1); upd(2, whr, hr, nlay); upd(3, wdir, dir, nlay + 1);
    ++n; ne += ex;
  }
  const char* nm[4] = {"swuflx", "swdflx", "swhr", "swdkdir"};
  std::printf("%ld columns, %ld bit-identical\n", n, ne);
  for (int v = 0; v < 4; ++v) std::printf("%-8s max|d| %.3e rel %.3e\n", nm[v], mx[v], mx[v] / sc[v]);
}
