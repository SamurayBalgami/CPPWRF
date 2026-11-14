// SPDX-License-Identifier: GPL-3.0-or-later
// The C++ side of the P3 column test: p3_init and run_p3 (CPU FP64) on the slab of generate_inputs.py, the same
// calls and output layout as driver_p3.f90.
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "phys/mp_p3/p3.hpp"

using namespace cppwrf;
using namespace cppwrf::phys::mp_p3;

int main(int argc, char** argv) {
  if (argc < 3) { std::fprintf(stderr, "usage: driver_p3_cpu IN OUT\n"); return 2; }
  FILE* f = std::fopen(argv[1], "rb");
  if (!f) return 1;
  int32_t hdr[4];
  double dt;
  int32_t start_it;
  if (std::fread(hdr, 4, 4, f) != 4 || std::fread(&dt, 8, 1, f) != 1 || std::fread(&start_it, 4, 1, f) != 1) return 1;
  const int n = hdr[0], nz = hdr[1], nsteps = hdr[2];
  const bool predict = hdr[3] == 1;
  const std::size_t m = std::size_t(n) * nz;
  std::vector<double> th(m), qv(m), qc(m), qr(m), qnr(m), qi(m), qni(m), qir(m), qib(m), qnc(m), tho(m), qvo(m), p(m),
      dz(m), refl(m, 0.), effc(m, 0.), effi(m, 0.), vmi(m, 0.), di(m, 0.), rhopo(m, 0.);
  for (auto* v : {&th, &qv, &qc, &qr, &qnr, &qi, &qni, &qir, &qib, &qnc, &tho, &qvo, &p, &dz})
    if (std::fread(v->data(), 8, m, f) != m) return 1;
  std::fclose(f);
  std::vector<double> rainnc(n, 0.), rainncv(n, 0.), sr(n, 0.), snownc(n, 0.), snowncv(n, 0.);
  const P3Host<double> h = p3_init<double>(".");
  P3Tables<double> tab{h.itab.data(), h.itabcoll.data(), h.vn.data(), h.vm.data(), h.revap.data()};
  // (i, k) slab with j = 0: View3 index i + k * n
  auto V = [&](std::vector<double>& a) { return View3<double>{a.data(), 0, n}; };
  auto CV = [&](std::vector<double>& a) { return View3<const double>{a.data(), 0, n}; };
  auto V2 = [&](std::vector<double>& a) { return View2<double>{a.data(), 0}; };
  FILE* o = std::fopen(argv[2], "wb");
  for (int s = 0; s < nsteps; ++s) {
    run_p3<double>(h.c, tab, dt, start_it + s, predict, Range2D{0, n, 0, 1}, nz, V(th), V(qv), V(qc), V(qr), V(qnr),
                   V(qi), V(qni), V(qir), V(qib), V(qnc), V(tho), V(qvo), CV(p), CV(dz), V(refl), V(effc), V(effi),
                   V(vmi), V(di), V(rhopo), V2(rainnc), V2(rainncv), V2(sr), V2(snownc), V2(snowncv));
    for (auto* v : {&th, &qv, &qc, &qr, &qnr, &qi, &qni, &qir, &qib, &qnc, &tho, &qvo, &refl, &effc, &effi, &vmi, &di,
                    &rhopo})
      std::fwrite(v->data(), 8, m, o);
    for (auto* v : {&rainnc, &rainncv, &sr, &snownc, &snowncv}) std::fwrite(v->data(), 8, std::size_t(n), o);
  }
  std::fclose(o);
  return 0;
}
