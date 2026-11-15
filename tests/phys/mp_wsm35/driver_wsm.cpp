// SPDX-License-Identifier: GPL-3.0-or-later
// The C++ side of the WSM5 / WSM3 column test: wsm5_params (wsm5init) and Wsm5Kernel + WSM6's Wsm6EffRadKernel, or
// Wsm3Kernel (its effective radii included), on the slab of generate_inputs.py with the lanes i0 = 0, ni as
// src/dyn/physics_coupling.cpp passes them; the output layout of driver_wsm5.f90.
#include <cstdint>
#include <cstdio>
#include <vector>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "phys/mp_wsm3/wsm3.hpp"
#include "phys/mp_wsm5/wsm5.hpp"

using namespace cppwrf;
using namespace cppwrf::phys;

int main(int argc, char** argv) {
  if (argc < 3) { std::fprintf(stderr, "usage: driver_wsm_cpu IN OUT\n"); return 2; }
  FILE* f = std::fopen(argv[1], "rb");
  if (!f) return 1;
  int32_t hdr[5];
  double dt;
  if (std::fread(hdr, 4, 5, f) != 5 || std::fread(&dt, 8, 1, f) != 1) return 1;
  const int ni = hdr[0], nj = hdr[1], nz = hdr[2], nsteps = hdr[3], scheme = hdr[4];
  const std::size_t m = std::size_t(ni) * nz * nj, m2 = std::size_t(ni) * nj;
  std::vector<double> th(m), q(m), qc(m), qr(m), qi(m), qs(m), den(m), pii(m), p(m), delz(m), w(m), rec(m, 0.),
      rei(m, 0.), res(m, 0.), rain(m2), rainncv(m2, 0.), snow(m2), snowncv(m2, 0.), sr(m2, 0.);
  for (auto* v : {&th, &q, &qc, &qr, &qi, &qs, &den, &pii, &p, &delz, &w})
    if (std::fread(v->data(), 8, m, f) != m) return 1;
  for (auto* v : {&rain, &snow})
    if (std::fread(v->data(), 8, m2, f) != m2) return 1;
  std::fclose(f);
  // Fortran (i, k, j) order: index i + k * ni + j * ni * nz
  const int sj = ni * nz, sk = ni;
  auto V = [&](std::vector<double>& a) { return View3<double>{a.data(), sj, sk}; };
  auto CV = [&](std::vector<double>& a) { return View3<const double>{a.data(), sj, sk}; };
  auto V2 = [&](std::vector<double>& a) { return View2<double>{a.data(), ni}; };
  const wsm5::Wsm5Params<double> c = wsm5::wsm5_params<double>();
  FILE* o = std::fopen(argv[2], "wb");
  const double cv[33] = {c.qc0, c.qck1, c.pidnc, c.g1pbr, c.g3pbr, c.g4pbr, c.g5pbro2, c.pvtr, c.pacrr, c.precr1,
                         c.precr2, c.roqimax, c.g1pbs, c.g3pbs, c.g4pbs, c.g5pbso2, c.pvts, c.pacrs, c.precs1,
                         c.precs2, c.pidn0r, c.pidn0s, c.xlv1, scheme == 5 ? c.pacrc : 0., c.pi, c.rslopermax,
                         c.rslopesmax, c.rsloperbmax, c.rslopesbmax, c.rsloper2max, c.rslopes2max, c.rsloper3max,
                         c.rslopes3max};
  std::fwrite(cv, 8, 33, o);
  const Range2D r{0, ni, 0, nj};
  for (int s = 0; s < nsteps; ++s) {
    if (scheme == 5) {
      column_for(r, wsm5::Wsm5Kernel<double, 128>{c, dt, nz, V(th), V(q), V(qc), V(qr), V(qi), V(qs), CV(den), CV(pii),
                                                  CV(p), CV(delz), V2(rain), V2(rainncv), V2(snow), V2(snowncv),
                                                  V2(sr), 0, ni});
      parallel_for(Range3D{0, ni, 0, nj, 0, nz},
                   wsm6::Wsm6EffRadKernel<double>{c.pidnc, c.pidn0s, c.alpha, c.n0s, c.n0smax, c.t0c, c.qmin, 2.49E-6,
                                                  4.99E-6, 9.99E-6, 50.E-6, 125.E-6, 999.E-6, CV(th), CV(pii), CV(qc),
                                                  CV(qi), CV(qs), CV(den), V(rec), V(rei), V(res)});
    } else {
      column_for(r, wsm3::Wsm3Kernel<double, 128>{c, dt, nz, V(th), V(q), V(qc), V(qr), CV(w), CV(den), CV(pii), CV(p),
                                                  CV(delz), V2(rain), V2(rainncv), V2(snow), V2(snowncv), V2(sr),
                                                  V(rec), V(rei), V(res), 0, ni});
    }
    for (auto* v : {&th, &q, &qc, &qr, &qi, &qs, &rec, &rei, &res}) std::fwrite(v->data(), 8, m, o);
    for (auto* v : {&rain, &rainncv, &snow, &snowncv, &sr}) std::fwrite(v->data(), 8, m2, o);
  }
  std::fclose(o);
  return 0;
}
