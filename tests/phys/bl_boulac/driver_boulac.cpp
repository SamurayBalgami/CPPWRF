// SPDX-License-Identifier: GPL-3.0-or-later
// The C++ BouLac PBL (phys/bl_boulac/boulac.hpp, CPU FP64) for the columns of generate_inputs.py; writes the outputs
// in driver_boulac.f90's order.
#include <cstdint>
#include <cstdio>
#include <vector>

#include "phys/bl_boulac/boulac.hpp"

using namespace cppwrf;

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: driver_boulac_cpu IN OUT\n");
    return 2;
  }
  FILE* f = std::fopen(argv[1], "rb");
  if (!f) return 1;
  int32_t dims[2];
  double sc[3];
  if (std::fread(dims, 4, 2, f) != 2 || std::fread(sc, 8, 3, f) != 3) return 1;
  const int n = dims[0], nz = dims[1];
  const std::size_t n3 = std::size_t(n) * nz;
  auto rd = [&](std::size_t m) {
    std::vector<double> v(m);
    if (std::fread(v.data(), 8, m, f) != m) std::fprintf(stderr, "short read\n");
    return v;
  };
  std::vector<double> u = rd(n3), v = rd(n3), th = rd(n3), rho = rd(n3), qv = rd(n3), qc = rd(n3), dz = rd(n3),
                      tke = rd(n3), ust = rd(n), hfx = rd(n), qfx = rd(n);
  std::fclose(f);
  std::vector<std::vector<double>> out(13, std::vector<double>(n3, 0.0));
  std::vector<double> pblh(n, 0.0);
  auto c3 = [&](const std::vector<double>& x) { return View3<const real>{x.data(), n, n}; };
  auto v3 = [&](std::vector<double>& x) { return View3<real>{x.data(), n, n}; };
  auto c2 = [&](const std::vector<double>& x) { return View2<const real>{x.data(), n}; };
  out[0] = tke;
  phys::boulac::BoulacKernel k{};
  k.u_phy = c3(u); k.v_phy = c3(v); k.th_phy = c3(th); k.rho = c3(rho); k.qv = c3(qv); k.qc = c3(qc);
  k.dz8w = c3(dz); k.ust = c2(ust); k.hfx = c2(hfx); k.qfx = c2(qfx); k.pblh = View2<real>{pblh.data(), n};
  k.tke = v3(out[0]); k.dlk = v3(out[1]); k.wu = v3(out[2]); k.wv = v3(out[3]); k.wt = v3(out[4]);
  k.wq = v3(out[5]); k.exch_h = v3(out[6]); k.exch_m = v3(out[7]); k.rublten = v3(out[8]); k.rvblten = v3(out[9]);
  k.rthblten = v3(out[10]); k.rqvblten = v3(out[11]); k.rqcblten = v3(out[12]);
  k.dt = sc[0]; k.cp = sc[1]; k.g = sc[2]; k.nz = nz;
  for (int i = 0; i < n; ++i) {
    if (nz <= 48)
      phys::boulac::BoulacRunner<48>{k}(i, 0);
    else
      phys::boulac::BoulacRunner<phys::boulac::kMaxLevels>{k}(i, 0);
  }
  FILE* o = std::fopen(argv[2], "wb");
  if (!o) return 1;
  for (auto& x : out) std::fwrite(x.data(), 8, n3, o);
  std::fwrite(pblh.data(), 8, n, o);
  std::fclose(o);
  return 0;
}
