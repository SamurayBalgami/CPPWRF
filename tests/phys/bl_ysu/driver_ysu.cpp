#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "phys/bl_ysu/ysu.hpp"

#ifdef __CUDACC__
#include <cuda_runtime.h>
#endif

using namespace cppwrf;

int main(int argc, char** argv) {
#ifdef __CUDACC__
  cudaDeviceSetLimit(cudaLimitStackSize, 64 * 1024);
#endif

  std::string infile = "ysu_input.bin";
  std::string outfile = "ysu_output_cpp.bin";
  bool bench_mode = false;
  int bench_iters = 10;

  bool infile_set = false;
  for (int a = 1; a < argc; ++a) {
    std::string arg = argv[a];
    if (arg == "--bench") {
      bench_mode = true;
      if (a + 1 < argc && argv[a + 1][0] != '-') {
        bench_iters = std::atoi(argv[++a]);
      }
    } else if (!infile_set) {
      infile = arg;
      infile_set = true;
    } else {
      outfile = arg;
    }
  }

  FILE* f = std::fopen(infile.c_str(), "rb");
  if (!f) {
    std::fprintf(stderr, "Cannot open input file: %s\n", infile.c_str());
    return 1;
  }

  int32_t n_cols = 0, nz = 0;
  if (std::fread(&n_cols, sizeof(int32_t), 1, f) != 1 ||
      std::fread(&nz, sizeof(int32_t), 1, f) != 1) {
    std::fprintf(stderr, "Failed reading dimensions from %s\n", infile.c_str());
    std::fclose(f);
    return 1;
  }

  double dt_in = 0.0;
  if (std::fread(&dt_in, sizeof(double), 1, f) != 1) {
    std::fprintf(stderr, "Failed reading dt\n");
    std::fclose(f);
    return 1;
  }

  int32_t ysu_topdown_pblmix_i = 0, flag_qc_i = 1, flag_qi_i = 1;
  if (std::fread(&ysu_topdown_pblmix_i, sizeof(int32_t), 1, f) != 1 ||
      std::fread(&flag_qc_i, sizeof(int32_t), 1, f) != 1 ||
      std::fread(&flag_qi_i, sizeof(int32_t), 1, f) != 1) {
    std::fprintf(stderr, "Failed reading flags\n");
    std::fclose(f);
    return 1;
  }

  double cp, g, rovcp, rd, rovg, ep1, ep2, karman, xlv, rv;
  if (std::fread(&cp, sizeof(double), 1, f) != 1 ||
      std::fread(&g, sizeof(double), 1, f) != 1 ||
      std::fread(&rovcp, sizeof(double), 1, f) != 1 ||
      std::fread(&rd, sizeof(double), 1, f) != 1 ||
      std::fread(&rovg, sizeof(double), 1, f) != 1 ||
      std::fread(&ep1, sizeof(double), 1, f) != 1 ||
      std::fread(&ep2, sizeof(double), 1, f) != 1 ||
      std::fread(&karman, sizeof(double), 1, f) != 1 ||
      std::fread(&xlv, sizeof(double), 1, f) != 1 ||
      std::fread(&rv, sizeof(double), 1, f) != 1) {
    std::fprintf(stderr, "Failed reading constants\n");
    std::fclose(f);
    return 1;
  }

  std::size_t n3d = std::size_t(n_cols) * std::size_t(nz);
  std::size_t n3di = std::size_t(n_cols) * std::size_t(nz + 1);
  std::size_t n2d = std::size_t(n_cols);

  // Read FP64 double arrays from file
  std::vector<double> h_u(n3d), h_v(n3d), h_t(n3d), h_qv(n3d);
  std::vector<double> h_qc(n3d), h_qi(n3d), h_p(n3d), h_pi(n3di), h_pii(n3d);
  std::vector<double> h_dz8w(n3d), h_rthraten(n3d);

  std::vector<double> h_psfcpa(n2d), h_znt(n2d), h_ust(n2d), h_psim(n2d), h_psih(n2d);
  std::vector<double> h_xland(n2d), h_hfx(n2d), h_qfx(n2d), h_wspd(n2d), h_br(n2d);
  std::vector<double> h_u10(n2d), h_v10(n2d), h_uoce(n2d), h_voce(n2d);

  auto read_arr = [&](std::vector<double>& arr, std::size_t count) {
    if (std::fread(arr.data(), sizeof(double), count, f) != count) {
      std::fprintf(stderr, "Failed reading array data\n");
      std::exit(1);
    }
  };

  read_arr(h_u, n3d); read_arr(h_v, n3d); read_arr(h_t, n3d); read_arr(h_qv, n3d);
  read_arr(h_qc, n3d); read_arr(h_qi, n3d); read_arr(h_p, n3d); read_arr(h_pi, n3di);
  read_arr(h_pii, n3d); read_arr(h_dz8w, n3d); read_arr(h_rthraten, n3d);

  read_arr(h_psfcpa, n2d); read_arr(h_znt, n2d); read_arr(h_ust, n2d); read_arr(h_psim, n2d); read_arr(h_psih, n2d);
  read_arr(h_xland, n2d); read_arr(h_hfx, n2d); read_arr(h_qfx, n2d); read_arr(h_wspd, n2d); read_arr(h_br, n2d);
  read_arr(h_u10, n2d); read_arr(h_v10, n2d); read_arr(h_uoce, n2d); read_arr(h_voce, n2d);
  std::fclose(f);

  // Convert to working precision `real`
  auto to_real = [](const std::vector<double>& in) {
    return std::vector<real>(in.begin(), in.end());
  };

  auto r_u = to_real(h_u), r_v = to_real(h_v), r_t = to_real(h_t), r_qv = to_real(h_qv);
  auto r_qc = to_real(h_qc), r_qi = to_real(h_qi), r_p = to_real(h_p), r_pi = to_real(h_pi);
  auto r_pii = to_real(h_pii), r_dz8w = to_real(h_dz8w), r_rthraten = to_real(h_rthraten);

  auto r_psfcpa = to_real(h_psfcpa), r_znt = to_real(h_znt), r_ust = to_real(h_ust);
  auto r_psim = to_real(h_psim), r_psih = to_real(h_psih), r_xland = to_real(h_xland);
  auto r_hfx = to_real(h_hfx), r_qfx = to_real(h_qfx), r_wspd = to_real(h_wspd), r_br = to_real(h_br);
  auto r_u10 = to_real(h_u10), r_v10 = to_real(h_v10), r_uoce = to_real(h_uoce), r_voce = to_real(h_voce);

  // Allocate device memory using Buffer
  Buffer d_u(n3d), d_v(n3d), d_t(n3d), d_qv(n3d);
  Buffer d_qc(n3d), d_qi(n3d), d_p(n3d), d_pi(n3di), d_pii(n3d);
  Buffer d_dz8w(n3d), d_rthraten(n3d);

  Buffer d_psfcpa(n2d), d_znt(n2d), d_ust(n2d), d_psim(n2d), d_psih(n2d);
  Buffer d_xland(n2d), d_hfx(n2d), d_qfx(n2d), d_wspd(n2d), d_br(n2d);
  Buffer d_u10(n2d), d_v10(n2d), d_uoce(n2d), d_voce(n2d);

  Buffer d_utnp(n3d), d_vtnp(n3d), d_ttnp(n3d), d_qvtnp(n3d), d_qctnp(n3d), d_qitnp(n3d);
  Buffer d_exch_h(n3d), d_exch_m(n3d);
  Buffer d_hpbl(n2d), d_kpbl1d(n2d * sizeof(int) / sizeof(real) + 1);
  Buffer d_wstar(n2d), d_delta(n2d);
  Buffer d_dusfc(n2d), d_dvsfc(n2d), d_dtsfc(n2d), d_dqsfc(n2d);

  auto copy_up = [&](Buffer& b, const std::vector<real>& v, std::size_t sz) {
    copy_to_dev(b.get(), v.data(), sz * sizeof(real));
  };

  copy_up(d_u, r_u, n3d); copy_up(d_v, r_v, n3d); copy_up(d_t, r_t, n3d); copy_up(d_qv, r_qv, n3d);
  copy_up(d_qc, r_qc, n3d); copy_up(d_qi, r_qi, n3d); copy_up(d_p, r_p, n3d); copy_up(d_pi, r_pi, n3di);
  copy_up(d_pii, r_pii, n3d); copy_up(d_dz8w, r_dz8w, n3d); copy_up(d_rthraten, r_rthraten, n3d);

  copy_up(d_psfcpa, r_psfcpa, n2d); copy_up(d_znt, r_znt, n2d); copy_up(d_ust, r_ust, n2d);
  copy_up(d_psim, r_psim, n2d); copy_up(d_psih, r_psih, n2d); copy_up(d_xland, r_xland, n2d);
  copy_up(d_hfx, r_hfx, n2d); copy_up(d_qfx, r_qfx, n2d); copy_up(d_wspd, r_wspd, n2d); copy_up(d_br, r_br, n2d);
  copy_up(d_u10, r_u10, n2d); copy_up(d_v10, r_v10, n2d); copy_up(d_uoce, r_uoce, n2d); copy_up(d_voce, r_voce, n2d);

  // Set up views
  // In binary layout, 3D element (i, k) is at i + k * n_cols
  int sj = n_cols * nz;
  int sk = n_cols;
  int sji = n_cols * (nz + 1);

  View3<const real> v_ux{d_u.get(), sj, sk};
  View3<const real> v_vx{d_v.get(), sj, sk};
  View3<const real> v_tx{d_t.get(), sj, sk};
  View3<const real> v_qvx{d_qv.get(), sj, sk};
  View3<const real> v_qcx{d_qc.get(), sj, sk};
  View3<const real> v_qix{d_qi.get(), sj, sk};
  View3<const real> v_p2d{d_p.get(), sj, sk};
  View3<const real> v_p2di{d_pi.get(), sji, sk};
  View3<const real> v_pi2d{d_pii.get(), sj, sk};
  View3<const real> v_dz8w2d{d_dz8w.get(), sj, sk};
  View3<const real> v_rthraten{d_rthraten.get(), sj, sk};

  View2<const real> v_psfcpa{d_psfcpa.get(), n_cols};
  View2<const real> v_znt{d_znt.get(), n_cols};
  View2<const real> v_ust{d_ust.get(), n_cols};
  View2<const real> v_psim{d_psim.get(), n_cols};
  View2<const real> v_psih{d_psih.get(), n_cols};
  View2<const real> v_xland{d_xland.get(), n_cols};
  View2<const real> v_hfx{d_hfx.get(), n_cols};
  View2<const real> v_qfx{d_qfx.get(), n_cols};
  View2<const real> v_wspd{d_wspd.get(), n_cols};
  View2<const real> v_br{d_br.get(), n_cols};
  View2<real> v_u10{d_u10.get(), n_cols};
  View2<real> v_v10{d_v10.get(), n_cols};
  View2<const real> v_uoce{d_uoce.get(), n_cols};
  View2<const real> v_voce{d_voce.get(), n_cols};

  View2<const real> v_ctopo{nullptr, 0};
  View2<const real> v_ctopo2{nullptr, 0};

  View3<real> v_utnp{d_utnp.get(), sj, sk};
  View3<real> v_vtnp{d_vtnp.get(), sj, sk};
  View3<real> v_ttnp{d_ttnp.get(), sj, sk};
  View3<real> v_qvtnp{d_qvtnp.get(), sj, sk};
  View3<real> v_qctnp{d_qctnp.get(), sj, sk};
  View3<real> v_qitnp{d_qitnp.get(), sj, sk};
  View3<real> v_exch_hx{d_exch_h.get(), sj, sk};
  View3<real> v_exch_mx{d_exch_m.get(), sj, sk};

  View2<real> v_hpbl{d_hpbl.get(), n_cols};
  View2<int> v_kpbl1d{reinterpret_cast<int*>(d_kpbl1d.get()), n_cols};
  View2<real> v_wstar{d_wstar.get(), n_cols};
  View2<real> v_delta{d_delta.get(), n_cols};

  View2<real> v_dusfc{d_dusfc.get(), n_cols};
  View2<real> v_dvsfc{d_dvsfc.get(), n_cols};
  View2<real> v_dtsfc{d_dtsfc.get(), n_cols};
  View2<real> v_dqsfc{d_dqsfc.get(), n_cols};

  Range2D range{0, n_cols, 0, 1};

  bool f_qc = (flag_qc_i == 1);
  bool f_qi = (flag_qi_i == 1);
  bool ysu_topdown_pblmix = (ysu_topdown_pblmix_i == 1);

  // Run kernel once
  phys::ysu::run_ysu(
      v_ux, v_vx, v_tx, v_qvx, v_qcx, v_qix, v_p2d, v_p2di, v_pi2d, v_dz8w2d, v_rthraten,
      v_psfcpa, v_znt, v_ust, v_psim, v_psih, v_xland, v_hfx, v_qfx, v_wspd, v_br,
      v_u10, v_v10, v_uoce, v_voce, v_ctopo, v_ctopo2,
      v_utnp, v_vtnp, v_ttnp, v_qvtnp, v_qctnp, v_qitnp, v_exch_hx, v_exch_mx,
      v_hpbl, v_kpbl1d, v_wstar, v_delta, v_dusfc, v_dvsfc, v_dtsfc, v_dqsfc,
      real(dt_in), range, nz,
      f_qc, f_qi, ysu_topdown_pblmix, false, true,
      real(cp), real(g), real(rovcp), real(rd), real(rovg),
      real(ep1), real(ep2), real(karman), real(xlv), real(rv));
  dev_sync();

  if (bench_mode) {
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < bench_iters; ++iter) {
      phys::ysu::run_ysu(
          v_ux, v_vx, v_tx, v_qvx, v_qcx, v_qix, v_p2d, v_p2di, v_pi2d, v_dz8w2d, v_rthraten,
          v_psfcpa, v_znt, v_ust, v_psim, v_psih, v_xland, v_hfx, v_qfx, v_wspd, v_br,
          v_u10, v_v10, v_uoce, v_voce, v_ctopo, v_ctopo2,
          v_utnp, v_vtnp, v_ttnp, v_qvtnp, v_qctnp, v_qitnp, v_exch_hx, v_exch_mx,
          v_hpbl, v_kpbl1d, v_wstar, v_delta, v_dusfc, v_dvsfc, v_dtsfc, v_dqsfc,
          real(dt_in), range, nz,
          f_qc, f_qi, ysu_topdown_pblmix, false, true,
          real(cp), real(g), real(rovcp), real(rd), real(rovg),
          real(ep1), real(ep2), real(karman), real(xlv), real(rv));
    }
    dev_sync();
    auto t1 = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    double avg_ms = total_ms / bench_iters;
    double ns_per_col = (avg_ms * 1e6) / n_cols;
    std::printf("BENCHMARK: %d cols x %d levels | %d iters | avg %.3f ms | %.1f ns/col\n",
                n_cols, nz, bench_iters, avg_ms, ns_per_col);
  }

  // Copy outputs back to host
  std::vector<real> r_utnp(n3d), r_vtnp(n3d), r_ttnp(n3d), r_qvtnp(n3d), r_qctnp(n3d), r_qitnp(n3d);
  std::vector<real> r_exch_h(n3d), r_exch_m(n3d);
  std::vector<real> r_hpbl(n2d), r_wstar(n2d), r_delta(n2d);
  std::vector<int> r_kpbl1d(n2d);

  auto copy_down = [&](std::vector<real>& v, Buffer& b, std::size_t sz) {
    copy_to_host(v.data(), b.get(), sz * sizeof(real));
  };

  copy_down(r_utnp, d_utnp, n3d); copy_down(r_vtnp, d_vtnp, n3d); copy_down(r_ttnp, d_ttnp, n3d);
  copy_down(r_qvtnp, d_qvtnp, n3d); copy_down(r_qctnp, d_qctnp, n3d); copy_down(r_qitnp, d_qitnp, n3d);
  copy_down(r_exch_h, d_exch_h, n3d); copy_down(r_exch_m, d_exch_m, n3d);
  copy_down(r_hpbl, d_hpbl, n2d); copy_down(r_wstar, d_wstar, n2d); copy_down(r_delta, d_delta, n2d);
  copy_to_host(r_kpbl1d.data(), d_kpbl1d.get(), n2d * sizeof(int));

  // Convert back to double for output file
  auto to_double = [](const std::vector<real>& in) {
    return std::vector<double>(in.begin(), in.end());
  };

  auto out_utnp = to_double(r_utnp), out_vtnp = to_double(r_vtnp), out_ttnp = to_double(r_ttnp);
  auto out_qvtnp = to_double(r_qvtnp), out_qctnp = to_double(r_qctnp), out_qitnp = to_double(r_qitnp);
  auto out_exch_h = to_double(r_exch_h), out_exch_m = to_double(r_exch_m);
  auto out_hpbl = to_double(r_hpbl), out_wstar = to_double(r_wstar), out_delta = to_double(r_delta);
  std::vector<double> out_kpbl1d(n2d);
  for (std::size_t c = 0; c < n2d; ++c) out_kpbl1d[c] = static_cast<double>(r_kpbl1d[c]);

  FILE* fo = std::fopen(outfile.c_str(), "wb");
  if (!fo) {
    std::fprintf(stderr, "Cannot open output file: %s\n", outfile.c_str());
    return 1;
  }

  auto write_arr = [&](const std::vector<double>& arr, std::size_t sz) {
    std::fwrite(arr.data(), sizeof(double), sz, fo);
  };

  write_arr(out_utnp, n3d); write_arr(out_vtnp, n3d); write_arr(out_ttnp, n3d);
  write_arr(out_qvtnp, n3d); write_arr(out_qctnp, n3d); write_arr(out_qitnp, n3d);
  write_arr(out_exch_h, n3d); write_arr(out_exch_m, n3d);
  write_arr(out_hpbl, n2d); write_arr(out_wstar, n2d); write_arr(out_delta, n2d);
  write_arr(out_kpbl1d, n2d);

  std::fclose(fo);
  return 0;
}
