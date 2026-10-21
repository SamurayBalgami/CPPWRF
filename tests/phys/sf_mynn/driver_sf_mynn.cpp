// CPPWRF side of the MYNN surface-layer differential test (see driver_sf_mynn.f90 and fields.txt):
// MynnSfcKernel on ni columns for nsteps calls, the same records out.
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "phys/sf_mynn/mynnsfc.hpp"

using namespace cppwrf;
using namespace cppwrf::phys::mynnsfc;

int main(int argc, char** argv) {
  if (argc < 3) { std::fprintf(stderr, "usage: %s IN.bin OUT.bin\n", argv[0]); return 2; }
  std::FILE* f = std::fopen(argv[1], "rb");
  if (!f) return 1;
  int hdr[4];
  if (std::fread(hdr, sizeof(int), 4, f) != 4) return 1;
  const int ni = hdr[0], itimestep = hdr[1], isfflx = hdr[2], nsteps = hdr[3];
  constexpr int n3 = 8, n2 = 9, nio = 34;
  std::vector<double> in3(size_t(ni) * 2 * n3), in2(size_t(ni) * n2), io2(size_t(ni) * nio);
  for (int n = 0; n < n3; ++n)
    for (int k = 0; k < 2; ++k)
      if (std::fread(&in3[(size_t(n) * 2 + k) * ni], sizeof(double), ni, f) != size_t(ni)) return 1;
  if (std::fread(in2.data(), sizeof(double), in2.size(), f) != in2.size()) return 1;
  if (std::fread(io2.data(), sizeof(double), io2.size(), f) != io2.size()) return 1;
  std::fclose(f);
  // View3 (i, j, k) = p[i + j*sj + k*sk]: one row (j = 0), k stride ni
  auto v3 = [&](int n) { return View3<const real>{&in3[size_t(n) * 2 * ni], ni, ni}; };
  auto v2 = [&](int n) { return View2<const real>{&in2[size_t(n) * ni], ni}; };
  auto w2 = [&](int n) { return View2<real>{&io2[size_t(n) * ni], ni}; };
  MynnSfcKernel k{};
  k.u3d = v3(0); k.v3d = v3(1); k.t3d = v3(2); k.qv3d = v3(3); k.p3d = v3(4); k.dz8w = v3(5); k.rho3d = v3(7);
  k.psfcpa = v2(0); k.mavail = v2(1); k.xland = v2(2); k.tsk = v2(3); k.snowh = v2(4); k.dx = v2(5); k.xice = v2(6);
  k.lakemask = v2(7);
  k.pblh = w2(0); k.chs = w2(1); k.chs2 = w2(2); k.cqs = w2(3); k.cqs2 = w2(4); k.cpm = w2(5); k.znt = w2(6);
  k.ust = w2(7); k.ustm = w2(8); k.zol = w2(9); k.mol = w2(10); k.rmol = w2(11); k.psim = w2(12); k.psih = w2(13);
  k.qgh = w2(14); k.hfx = w2(15); k.qfx = w2(16); k.lh = w2(17); k.flhc = w2(18); k.flqc = w2(19); k.qsfc = w2(20);
  k.u10 = w2(21); k.v10 = w2(22); k.th2 = w2(23); k.t2 = w2(24); k.q2 = w2(25); k.gz1oz0 = w2(26); k.wspd = w2(27);
  k.br = w2(28); k.ch = w2(29); k.ck = w2(30); k.cka = w2(31); k.cd = w2(32); k.cda = w2(33);
  static PsiTables tables;
  psi_init(tables);
  k.tb = &tables;
  k.isfflx = isfflx;
  k.restart = false;
  k.xice_threshold = real(0.5);
  for (int s = 0; s < nsteps; ++s) {
    k.itimestep = itimestep + s;
    for (int i = 0; i < ni; ++i) k(i, 0);
  }
  std::FILE* g = std::fopen(argv[2], "wb");
  if (!g) return 1;
  std::fwrite(io2.data(), sizeof(double), io2.size(), g);
  std::fclose(g);
  return 0;
}
