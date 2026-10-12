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
#include "phys/mp_wsm6/wsm6.hpp"

using namespace cppwrf;
using namespace cppwrf::phys::wsm6;

int main(int argc, char** argv) {
  std::string infile = "wsm6_input.bin";
  std::string outfile = "wsm6_output_cpp.bin";
  bool bench_mode = false;
  int bench_iters = 10;

  for (int a = 1; a < argc; ++a) {
    std::string arg = argv[a];
    if (arg == "--bench") {
      bench_mode = true;
      if (a + 1 < argc && argv[a + 1][0] != '-') {
        bench_iters = std::atoi(argv[++a]);
      }
    } else if (infile == "wsm6_input.bin") {
      infile = arg;
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

  double delt = 0.0, g = 0.0, cpd = 0.0, cpv = 0.0, rd = 0.0, rv = 0.0, t0c = 0.0;
  double ep1 = 0.0, ep2 = 0.0, qmin = 0.0, xls = 0.0, xlv0 = 0.0, xlf0 = 0.0;
  double den0 = 0.0, denr = 0.0, dens = 0.0, cliq = 0.0, cice = 0.0, psat = 0.0;
  int32_t hail_opt = 0;

  double consts[19];
  if (std::fread(consts, sizeof(double), 19, f) != 19 ||
      std::fread(&hail_opt, sizeof(int32_t), 1, f) != 1) {
    std::fprintf(stderr, "Failed reading constants from %s\n", infile.c_str());
    std::fclose(f);
    return 1;
  }

  delt = consts[0]; g = consts[1]; cpd = consts[2]; cpv = consts[3];
  rd = consts[4]; rv = consts[5]; t0c = consts[6]; ep1 = consts[7];
  ep2 = consts[8]; qmin = consts[9]; xls = consts[10]; xlv0 = consts[11];
  xlf0 = consts[12]; den0 = consts[13]; denr = consts[14]; dens = consts[15];
  cliq = consts[16]; cice = consts[17]; psat = consts[18];

  std::size_t n3d = std::size_t(n_cols) * std::size_t(nz);
  std::size_t n2d = std::size_t(n_cols);

  // Read FP64 double arrays from file
  std::vector<double> h_th(n3d), h_q(n3d), h_qc(n3d), h_qr(n3d), h_qi(n3d), h_qs(n3d), h_qg(n3d);
  std::vector<double> h_den(n3d), h_pii(n3d), h_p(n3d), h_delz(n3d);
  std::vector<double> h_rain(n2d), h_rainncv(n2d), h_snow(n2d), h_snowncv(n2d);
  std::vector<double> h_graupel(n2d), h_graupelncv(n2d), h_sr(n2d);

  if (std::fread(h_th.data(), sizeof(double), n3d, f) != n3d ||
      std::fread(h_q.data(), sizeof(double), n3d, f) != n3d ||
      std::fread(h_qc.data(), sizeof(double), n3d, f) != n3d ||
      std::fread(h_qr.data(), sizeof(double), n3d, f) != n3d ||
      std::fread(h_qi.data(), sizeof(double), n3d, f) != n3d ||
      std::fread(h_qs.data(), sizeof(double), n3d, f) != n3d ||
      std::fread(h_qg.data(), sizeof(double), n3d, f) != n3d ||
      std::fread(h_den.data(), sizeof(double), n3d, f) != n3d ||
      std::fread(h_pii.data(), sizeof(double), n3d, f) != n3d ||
      std::fread(h_p.data(), sizeof(double), n3d, f) != n3d ||
      std::fread(h_delz.data(), sizeof(double), n3d, f) != n3d ||
      std::fread(h_rain.data(), sizeof(double), n2d, f) != n2d ||
      std::fread(h_rainncv.data(), sizeof(double), n2d, f) != n2d ||
      std::fread(h_snow.data(), sizeof(double), n2d, f) != n2d ||
      std::fread(h_snowncv.data(), sizeof(double), n2d, f) != n2d ||
      std::fread(h_graupel.data(), sizeof(double), n2d, f) != n2d ||
      std::fread(h_graupelncv.data(), sizeof(double), n2d, f) != n2d ||
      std::fread(h_sr.data(), sizeof(double), n2d, f) != n2d) {
    std::fprintf(stderr, "Failed reading array data from %s\n", infile.c_str());
    std::fclose(f);
    return 1;
  }
  std::fclose(f);

  // Convert to working precision `real`
  std::vector<real> r_th(h_th.begin(), h_th.end());
  std::vector<real> r_q(h_q.begin(), h_q.end());
  std::vector<real> r_qc(h_qc.begin(), h_qc.end());
  std::vector<real> r_qr(h_qr.begin(), h_qr.end());
  std::vector<real> r_qi(h_qi.begin(), h_qi.end());
  std::vector<real> r_qs(h_qs.begin(), h_qs.end());
  std::vector<real> r_qg(h_qg.begin(), h_qg.end());

  std::vector<real> r_den(h_den.begin(), h_den.end());
  std::vector<real> r_pii(h_pii.begin(), h_pii.end());
  std::vector<real> r_p(h_p.begin(), h_p.end());
  std::vector<real> r_delz(h_delz.begin(), h_delz.end());

  std::vector<real> r_rain(h_rain.begin(), h_rain.end());
  std::vector<real> r_rainncv(h_rainncv.begin(), h_rainncv.end());
  std::vector<real> r_snow(h_snow.begin(), h_snow.end());
  std::vector<real> r_snowncv(h_snowncv.begin(), h_snowncv.end());
  std::vector<real> r_graupel(h_graupel.begin(), h_graupel.end());
  std::vector<real> r_graupelncv(h_graupelncv.begin(), h_graupelncv.end());
  std::vector<real> r_sr(h_sr.begin(), h_sr.end());

  // Allocate device memory using Buffer
  Buffer d_th(n3d), d_q(n3d), d_qc(n3d), d_qr(n3d);
  Buffer d_qi(n3d), d_qs(n3d), d_qg(n3d);
  Buffer d_den(n3d), d_pii(n3d), d_p(n3d), d_delz(n3d);

  Buffer d_rain(n2d), d_rainncv(n2d);
  Buffer d_snow(n2d), d_snowncv(n2d);
  Buffer d_graupel(n2d), d_graupelncv(n2d);
  Buffer d_sr(n2d);

  copy_to_dev(d_th.get(), r_th.data(), n3d * sizeof(real));
  copy_to_dev(d_q.get(), r_q.data(), n3d * sizeof(real));
  copy_to_dev(d_qc.get(), r_qc.data(), n3d * sizeof(real));
  copy_to_dev(d_qr.get(), r_qr.data(), n3d * sizeof(real));
  copy_to_dev(d_qi.get(), r_qi.data(), n3d * sizeof(real));
  copy_to_dev(d_qs.get(), r_qs.data(), n3d * sizeof(real));
  copy_to_dev(d_qg.get(), r_qg.data(), n3d * sizeof(real));

  copy_to_dev(d_den.get(), r_den.data(), n3d * sizeof(real));
  copy_to_dev(d_pii.get(), r_pii.data(), n3d * sizeof(real));
  copy_to_dev(d_p.get(), r_p.data(), n3d * sizeof(real));
  copy_to_dev(d_delz.get(), r_delz.data(), n3d * sizeof(real));

  copy_to_dev(d_rain.get(), r_rain.data(), n2d * sizeof(real));
  copy_to_dev(d_rainncv.get(), r_rainncv.data(), n2d * sizeof(real));
  copy_to_dev(d_snow.get(), r_snow.data(), n2d * sizeof(real));
  copy_to_dev(d_snowncv.get(), r_snowncv.data(), n2d * sizeof(real));
  copy_to_dev(d_graupel.get(), r_graupel.data(), n2d * sizeof(real));
  copy_to_dev(d_graupelncv.get(), r_graupelncv.data(), n2d * sizeof(real));
  copy_to_dev(d_sr.get(), r_sr.data(), n2d * sizeof(real));

  int sj = n_cols;
  int sk = n_cols;

  View3<real> v_th{d_th.get(), sj, sk};
  View3<real> v_q{d_q.get(), sj, sk};
  View3<real> v_qc{d_qc.get(), sj, sk};
  View3<real> v_qr{d_qr.get(), sj, sk};
  View3<real> v_qi{d_qi.get(), sj, sk};
  View3<real> v_qs{d_qs.get(), sj, sk};
  View3<real> v_qg{d_qg.get(), sj, sk};

  View3<const real> v_den{d_den.get(), sj, sk};
  View3<const real> v_pii{d_pii.get(), sj, sk};
  View3<const real> v_p{d_p.get(), sj, sk};
  View3<const real> v_delz{d_delz.get(), sj, sk};

  View2<real> v_rain{d_rain.get(), sj};
  View2<real> v_rainncv{d_rainncv.get(), sj};
  View2<real> v_snow{d_snow.get(), sj};
  View2<real> v_snowncv{d_snowncv.get(), sj};
  View2<real> v_graupel{d_graupel.get(), sj};
  View2<real> v_graupelncv{d_graupelncv.get(), sj};
  View2<real> v_sr{d_sr.get(), sj};

  Range2D range{0, n_cols, 0, 1};

  Wsm6Params<real> params;
  params.g = real(g);
  params.cpd = real(cpd);
  params.cpv = real(cpv);
  params.rd = real(rd);
  params.rv = real(rv);
  params.t0c = real(t0c);
  params.ep1 = real(ep1);
  params.ep2 = real(ep2);
  params.qmin = real(qmin);
  params.xls = real(xls);
  params.xlv0 = real(xlv0);
  params.xlf0 = real(xlf0);
  params.cliq = real(cliq);
  params.cice = real(cice);
  params.psat = real(psat);

  params.init(real(den0), real(denr), real(dens), real(cliq), real(cpv), hail_opt);

#if defined(CPPWRF_USE_CUDA)
  cudaDeviceSetLimit(cudaLimitStackSize, 64 * 1024);
#endif

  // Warmup run
  run_wsm6(v_th, v_q, v_qc, v_qr, v_qi, v_qs, v_qg,
           v_den, v_pii, v_p, v_delz,
           v_rain, v_rainncv, v_snow, v_snowncv,
           v_graupel, v_graupelncv, v_sr,
           real(delt), range, nz, params);
  dev_sync();

#if defined(CPPWRF_USE_CUDA)
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    std::fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(err));
  }
#endif

  if (bench_mode) {
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < bench_iters; ++iter) {
      run_wsm6(v_th, v_q, v_qc, v_qr, v_qi, v_qs, v_qg,
               v_den, v_pii, v_p, v_delz,
               v_rain, v_rainncv, v_snow, v_snowncv,
               v_graupel, v_graupelncv, v_sr,
               real(delt), range, nz, params);
    }
    dev_sync();
    auto t1 = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    double ms_per_iter = total_ms / bench_iters;
    double ns_per_col = (ms_per_iter * 1e6) / n_cols;

    std::printf("WSM6 Benchmark (%s): %d iters, %.3f ms/iter, %.1f ns/col (cols=%d, nz=%d)\n",
                backend_name(), bench_iters, ms_per_iter, ns_per_col, n_cols, nz);
  }

  // Copy results back to host
  copy_to_host(r_th.data(), d_th.get(), n3d * sizeof(real));
  copy_to_host(r_q.data(), d_q.get(), n3d * sizeof(real));
  copy_to_host(r_qc.data(), d_qc.get(), n3d * sizeof(real));
  copy_to_host(r_qr.data(), d_qr.get(), n3d * sizeof(real));
  copy_to_host(r_qi.data(), d_qi.get(), n3d * sizeof(real));
  copy_to_host(r_qs.data(), d_qs.get(), n3d * sizeof(real));
  copy_to_host(r_qg.data(), d_qg.get(), n3d * sizeof(real));

  copy_to_host(r_rain.data(), d_rain.get(), n2d * sizeof(real));
  copy_to_host(r_rainncv.data(), d_rainncv.get(), n2d * sizeof(real));
  copy_to_host(r_snow.data(), d_snow.get(), n2d * sizeof(real));
  copy_to_host(r_snowncv.data(), d_snowncv.get(), n2d * sizeof(real));
  copy_to_host(r_graupel.data(), d_graupel.get(), n2d * sizeof(real));
  copy_to_host(r_graupelncv.data(), d_graupelncv.get(), n2d * sizeof(real));
  copy_to_host(r_sr.data(), d_sr.get(), n2d * sizeof(real));

  // Convert to FP64 for output file
  std::vector<double> out_th(r_th.begin(), r_th.end());
  std::vector<double> out_q(r_q.begin(), r_q.end());
  std::vector<double> out_qc(r_qc.begin(), r_qc.end());
  std::vector<double> out_qr(r_qr.begin(), r_qr.end());
  std::vector<double> out_qi(r_qi.begin(), r_qi.end());
  std::vector<double> out_qs(r_qs.begin(), r_qs.end());
  std::vector<double> out_qg(r_qg.begin(), r_qg.end());

  std::vector<double> out_rain(r_rain.begin(), r_rain.end());
  std::vector<double> out_rainncv(r_rainncv.begin(), r_rainncv.end());
  std::vector<double> out_snow(r_snow.begin(), r_snow.end());
  std::vector<double> out_snowncv(r_snowncv.begin(), r_snowncv.end());
  std::vector<double> out_graupel(r_graupel.begin(), r_graupel.end());
  std::vector<double> out_graupelncv(r_graupelncv.begin(), r_graupelncv.end());
  std::vector<double> out_sr(r_sr.begin(), r_sr.end());

  FILE* fo = std::fopen(outfile.c_str(), "wb");
  if (!fo) {
    std::fprintf(stderr, "Cannot open output file: %s\n", outfile.c_str());
    return 1;
  }

  std::fwrite(&n_cols, sizeof(int32_t), 1, fo);
  std::fwrite(&nz, sizeof(int32_t), 1, fo);

  std::fwrite(out_th.data(), sizeof(double), n3d, fo);
  std::fwrite(out_q.data(), sizeof(double), n3d, fo);
  std::fwrite(out_qc.data(), sizeof(double), n3d, fo);
  std::fwrite(out_qr.data(), sizeof(double), n3d, fo);
  std::fwrite(out_qi.data(), sizeof(double), n3d, fo);
  std::fwrite(out_qs.data(), sizeof(double), n3d, fo);
  std::fwrite(out_qg.data(), sizeof(double), n3d, fo);

  std::fwrite(out_rain.data(), sizeof(double), n2d, fo);
  std::fwrite(out_rainncv.data(), sizeof(double), n2d, fo);
  std::fwrite(out_snow.data(), sizeof(double), n2d, fo);
  std::fwrite(out_snowncv.data(), sizeof(double), n2d, fo);
  std::fwrite(out_graupel.data(), sizeof(double), n2d, fo);
  std::fwrite(out_graupelncv.data(), sizeof(double), n2d, fo);
  std::fwrite(out_sr.data(), sizeof(double), n2d, fo);

  std::fclose(fo);
  return 0;
}
