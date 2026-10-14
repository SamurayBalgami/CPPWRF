// Check: WRF's rrtmg_lw core inputs (dumped by a debug WRF build at the first radiation
// call, big-endian stream) -> C++ rrtmg_lw_core -> compare with WRF's core outputs.
//   real_core_check <rrtmg_lw_core.bin> <rrtmg_lw_tables.bin>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <memory>
#include <vector>

#include "phys/ra_rrtmg_lw/rrtmg_lw.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_tables.hpp"

using namespace cppwrf::phys;
namespace cppwrf::phys { extern std::FILE* g_rrtmg_lw_dbg; }

struct BE {
  std::vector<char> b;
  size_t o = 0;
  bool more() const { return o < b.size(); }
  int i4() { uint32_t v; std::memcpy(&v, &b[o], 4); o += 4; v = __builtin_bswap32(v); int r; std::memcpy(&r, &v, 4); return r; }
  double f8() { uint64_t v; std::memcpy(&v, &b[o], 8); o += 8; v = __builtin_bswap64(v); double r; std::memcpy(&r, &v, 8); return r; }
  void arr(double* a, int n) { for (int k = 0; k < n; ++k) a[k] = f8(); }
};

int main(int argc, char** argv) {
  BE in;
  {
    std::ifstream f(argv[1], std::ios::binary);
    in.b.assign(std::istreambuf_iterator<char>(f), {});
  }
  auto tables = std::make_unique<RrtmgLwTables>();
  if (!load_rrtmg_lw_tables(argv[2], *tables)) { std::fprintf(stderr, "tables\n"); return 1; }
  auto ws = std::make_unique<RrtmgLwWorkspace>();
  ws->init_pointers();
  const int M = RRTMG_LW_MAX_LAYERS + 1;
  std::vector<double> play(M), plev(M), tlay(M), tlev(M), h2o(M), o3(M), co2(M), ch4(M), n2o(M), o2(M), c11(M), c12(M),
      c22(M), ccl4(M), emis(16), rei(M, 10.), rel(M, 5.), res(M, 10.);
  std::vector<double> wu(M), wd(M), whr(M), wuc(M), wdc(M), whrc(M), cu(M), cd(M), chr(M), cuc(M), cdc(M), chrc(M);
  double mx[6] = {0}, sc[6] = {0};
  int ncol = 0, ncol_exact = 0;
  while (in.more()) {
    const int i = in.i4(), j = in.i4(), nlay = in.i4();
    in.i4(); in.i4(); in.i4(); in.i4();  // icld, inflglw, iceflglw, liqflglw
    (void)i; (void)j;
    in.arr(play.data(), nlay); in.arr(plev.data(), nlay + 1); in.arr(tlay.data(), nlay); in.arr(tlev.data(), nlay + 1);
    const double tsfc = in.f8();
    for (auto* a : {&h2o, &o3, &co2, &ch4, &n2o, &o2, &c11, &c12, &c22, &ccl4}) in.arr(a->data(), nlay);
    in.arr(emis.data(), 16);
    in.arr(wu.data(), nlay + 1); in.arr(wd.data(), nlay + 1); in.arr(whr.data(), nlay);
    in.arr(wuc.data(), nlay + 1); in.arr(wdc.data(), nlay + 1); in.arr(whrc.data(), nlay);
    for (int g = 0; g < 140; ++g)
      for (int l = 0; l < nlay; ++l) ws->cldfmc[g][l] = ws->ciwpmc[g][l] = ws->clwpmc[g][l] = ws->cswpmc[g][l] = ws->taucmc[g][l] = 0;
    rrtmg_lw_core(nlay, 1004.5, play.data(), plev.data(), tlay.data(), tlev.data(), tsfc, h2o.data(), o3.data(),
                  co2.data(), ch4.data(), n2o.data(), o2.data(), c11.data(), c12.data(), c22.data(), ccl4.data(),
                  emis.data(), rei.data(), rel.data(), res.data(), cu.data(), cd.data(), chr.data(), cuc.data(),
                  cdc.data(), chrc.data(), *ws, *tables);
    bool exact = true;
    auto upd = [&](int v, const std::vector<double>& w, const std::vector<double>& c, int n) {
      for (int k = 0; k < n; ++k) {
        mx[v] = std::max(mx[v], std::abs(w[k] - c[k]));
        sc[v] = std::max(sc[v], std::abs(w[k]));
        if (w[k] != c[k]) exact = false;
      }
    };
    upd(0, wu, cu, nlay + 1); upd(1, wd, cd, nlay + 1); upd(2, whr, chr, nlay);
    upd(3, wuc, cuc, nlay + 1); upd(4, wdc, cdc, nlay + 1); upd(5, whrc, chrc, nlay);
    ++ncol; ncol_exact += exact;
  }
  const char* nm[6] = {"uflx", "dflx", "hr", "uflxc", "dflxc", "hrc"};
  std::printf("%d columns, %d bit-identical\n", ncol, ncol_exact);
  for (int v = 0; v < 6; ++v) std::printf("%-6s max|d| %.3e  rel %.3e\n", nm[v], mx[v], mx[v] / sc[v]);
  return 0;
}
