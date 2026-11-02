// WSM6's radar reflectivity (do_radar_ref = 1) in C++: radar_init with mp_wsm6_init's m(D) of rain, snow and
// graupel (hail) and Wsm6ReflKernel per column; reads the WSM6 test's input, writes what driver_radar.f90 writes.
//   driver_radar_cpu INPUT OUTPUT
#include <cstdint>
#include <cstdio>
#include <vector>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "phys/mp_wsm6/wsm6.hpp"

using namespace cppwrf;
using namespace cppwrf::phys;

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: driver_radar_cpu INPUT OUTPUT\n");
    return 2;
  }
  FILE* f = std::fopen(argv[1], "rb");
  if (!f) {
    std::fprintf(stderr, "Cannot open input file: %s\n", argv[1]);
    return 1;
  }
  int32_t n_cols = 0, nz = 0, hail_opt = 0;
  double c[19];
  if (std::fread(&n_cols, sizeof(int32_t), 1, f) != 1 || std::fread(&nz, sizeof(int32_t), 1, f) != 1 ||
      std::fread(c, sizeof(double), 19, f) != 19 || std::fread(&hail_opt, sizeof(int32_t), 1, f) != 1) {
    std::fprintf(stderr, "Failed reading the header of %s\n", argv[1]);
    return 1;
  }
  const std::size_t n3d = std::size_t(n_cols) * std::size_t(nz);
  // th, q, qc, qr, qi, qs, qg, den, pii, p
  std::vector<std::vector<real>> a(10, std::vector<real>(n3d));
  std::vector<double> tmp(n3d);
  for (auto& v : a) {
    if (std::fread(tmp.data(), sizeof(double), n3d, f) != n3d) {
      std::fprintf(stderr, "Failed reading the arrays of %s\n", argv[1]);
      return 1;
    }
    for (std::size_t n = 0; n < n3d; ++n) v[n] = real(tmp[n]);
  }
  std::fclose(f);

  // mp_wsm6_init(den0, denr, dens, cliq, cpv, hail_opt)
  wsm6::Wsm6Params<real> wp;
  wp.init(real(c[13]), real(c[14]), real(c[15]), real(c[16]), real(c[3]), hail_opt);
  const radar::RadarTables<real> rt = radar::radar_init<real>(
      wp.pi * wp.denr / real(6.), real(3.), real(0.), wp.pi * wp.dens / real(6.), real(3.), real(0.),
      wp.pi * wp.deng / real(6.), real(3.), real(0.));

  std::vector<real> refl(n3d, real(0));
  const int sj = n_cols * nz, sk = n_cols;
  auto cv = [&](int n) { return View3<const real>{a[n].data(), sj, sk}; };
  const wsm6::Wsm6ReflKernel<real> kernel{&rt, wp.n0r, wp.n0s, wp.n0smax, wp.alpha, wp.n0g, wp.pi, nz,
                                          cv(0), cv(8), cv(9), cv(1), cv(3), cv(5), cv(6),
                                          View3<real>{refl.data(), sj, sk}};
  for (int i = 0; i < n_cols; ++i) kernel(i, 0);

  FILE* fo = std::fopen(argv[2], "wb");
  if (!fo) {
    std::fprintf(stderr, "Cannot open output file: %s\n", argv[2]);
    return 1;
  }
  auto put = [&](double x) { std::fwrite(&x, sizeof(double), 1, fo); };
  for (const double* t : {rt.xxDs, rt.xdts, rt.xxDg, rt.xdtg})
    for (int n = 0; n < radar::nrbins; ++n) put(t[n]);
  for (int n = 0; n <= radar::nrbins; ++n) put(rt.simpson[n]);
  for (double x : {rt.k_w, rt.pi5, rt.lamda4, rt.m_w_0.re, rt.m_w_0.im, rt.m_i_0.re, rt.m_i_0.im}) put(x);
  for (const real* t : {rt.xcre, rt.xcse, rt.xcge, rt.xcrg, rt.xcsg, rt.xcgg})
    for (int n = 0; n < 4; ++n) put(double(t[n]));
  for (real x : {rt.xam_r, rt.xbm_r, rt.xmu_r, rt.xobmr, rt.xam_s, rt.xbm_s, rt.xmu_s, rt.xoams, rt.xobms,
                 rt.xocms, rt.xam_g, rt.xbm_g, rt.xmu_g, rt.xoamg, rt.xobmg, rt.xocmg, rt.xorg2, rt.xosg2,
                 rt.xogg2})
    put(double(x));
  for (real x : refl) put(double(x));
  std::fclose(fo);
  return 0;
}
