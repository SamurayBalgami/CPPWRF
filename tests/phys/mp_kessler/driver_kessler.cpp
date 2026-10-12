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
#include "phys/mp_kessler/kessler.hpp"

using namespace cppwrf;

int main(int argc, char** argv) {
  std::string infile = "kessler_input.bin";
  std::string outfile = "kessler_output_cpp.bin";
  bool bench_mode = false;
  int bench_iters = 10;

  for (int a = 1; a < argc; ++a) {
    std::string arg = argv[a];
    if (arg == "--bench") {
      bench_mode = true;
      if (a + 1 < argc && argv[a + 1][0] != '-') {
        bench_iters = std::atoi(argv[++a]);
      }
    } else if (infile == "kessler_input.bin") {
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

  double dt_in = 0.0;
  double c_xlv = 0.0, c_cp = 0.0, c_ep2 = 0.0, c_svp1 = 0.0;
  double c_svp2 = 0.0, c_svp3 = 0.0, c_svpt0 = 0.0, c_rhowater = 0.0;

  if (std::fread(&dt_in, sizeof(double), 1, f) != 1 ||
      std::fread(&c_xlv, sizeof(double), 1, f) != 1 ||
      std::fread(&c_cp, sizeof(double), 1, f) != 1 ||
      std::fread(&c_ep2, sizeof(double), 1, f) != 1 ||
      std::fread(&c_svp1, sizeof(double), 1, f) != 1 ||
      std::fread(&c_svp2, sizeof(double), 1, f) != 1 ||
      std::fread(&c_svp3, sizeof(double), 1, f) != 1 ||
      std::fread(&c_svpt0, sizeof(double), 1, f) != 1 ||
      std::fread(&c_rhowater, sizeof(double), 1, f) != 1) {
    std::fprintf(stderr, "Failed reading constants from %s\n", infile.c_str());
    std::fclose(f);
    return 1;
  }

  std::size_t n3d = std::size_t(n_cols) * std::size_t(nz);
  std::size_t n2d = std::size_t(n_cols);

  // Read FP64 double arrays from file
  std::vector<double> h_t(n3d), h_qv(n3d), h_qc(n3d), h_qr(n3d);
  std::vector<double> h_rho(n3d), h_pii(n3d), h_dz8w(n3d), h_z(n3d);
  std::vector<double> h_rainnc(n2d), h_rainncv(n2d);

  if (std::fread(h_t.data(), sizeof(double), n3d, f) != n3d ||
      std::fread(h_qv.data(), sizeof(double), n3d, f) != n3d ||
      std::fread(h_qc.data(), sizeof(double), n3d, f) != n3d ||
      std::fread(h_qr.data(), sizeof(double), n3d, f) != n3d ||
      std::fread(h_rho.data(), sizeof(double), n3d, f) != n3d ||
      std::fread(h_pii.data(), sizeof(double), n3d, f) != n3d ||
      std::fread(h_dz8w.data(), sizeof(double), n3d, f) != n3d ||
      std::fread(h_z.data(), sizeof(double), n3d, f) != n3d ||
      std::fread(h_rainnc.data(), sizeof(double), n2d, f) != n2d ||
      std::fread(h_rainncv.data(), sizeof(double), n2d, f) != n2d) {
    std::fprintf(stderr, "Failed reading array data from %s\n", infile.c_str());
    std::fclose(f);
    return 1;
  }
  std::fclose(f);

  // Convert to working precision `real`
  std::vector<real> r_t(h_t.begin(), h_t.end());
  std::vector<real> r_qv(h_qv.begin(), h_qv.end());
  std::vector<real> r_qc(h_qc.begin(), h_qc.end());
  std::vector<real> r_qr(h_qr.begin(), h_qr.end());
  std::vector<real> r_rho(h_rho.begin(), h_rho.end());
  std::vector<real> r_pii(h_pii.begin(), h_pii.end());
  std::vector<real> r_dz8w(h_dz8w.begin(), h_dz8w.end());
  std::vector<real> r_z(h_z.begin(), h_z.end());
  std::vector<real> r_rainnc(h_rainnc.begin(), h_rainnc.end());
  std::vector<real> r_rainncv(h_rainncv.begin(), h_rainncv.end());

  // Allocate device memory using Buffer
  Buffer d_t(n3d), d_qv(n3d), d_qc(n3d), d_qr(n3d);
  Buffer d_rho(n3d), d_pii(n3d), d_dz8w(n3d), d_z(n3d);
  Buffer d_rainnc(n2d), d_rainncv(n2d);

  // Copy host to device
  copy_to_dev(d_t.get(), r_t.data(), n3d * sizeof(real));
  copy_to_dev(d_qv.get(), r_qv.data(), n3d * sizeof(real));
  copy_to_dev(d_qc.get(), r_qc.data(), n3d * sizeof(real));
  copy_to_dev(d_qr.get(), r_qr.data(), n3d * sizeof(real));
  copy_to_dev(d_rho.get(), r_rho.data(), n3d * sizeof(real));
  copy_to_dev(d_pii.get(), r_pii.data(), n3d * sizeof(real));
  copy_to_dev(d_dz8w.get(), r_dz8w.data(), n3d * sizeof(real));
  copy_to_dev(d_z.get(), r_z.data(), n3d * sizeof(real));
  copy_to_dev(d_rainnc.get(), r_rainnc.data(), n2d * sizeof(real));
  copy_to_dev(d_rainncv.get(), r_rainncv.data(), n2d * sizeof(real));

  // Set up views
  // Indexing: element for column i (0 <= i < n_cols), level k (0 <= k < nz)
  // is at offset i + k * n_cols.
  int sj = n_cols * nz;
  int sk = n_cols;

  View3<real> v_t{d_t.get(), sj, sk};
  View3<real> v_qv{d_qv.get(), sj, sk};
  View3<real> v_qc{d_qc.get(), sj, sk};
  View3<real> v_qr{d_qr.get(), sj, sk};
  View3<const real> v_rho{d_rho.get(), sj, sk};
  View3<const real> v_pii{d_pii.get(), sj, sk};
  View3<const real> v_dz8w{d_dz8w.get(), sj, sk};
  View3<const real> v_z{d_z.get(), sj, sk};
  View2<real> v_rainnc{d_rainnc.get(), n_cols};
  View2<real> v_rainncv{d_rainncv.get(), n_cols};

  Range2D range{0, n_cols, 0, 1};

  // Run kernel once
  phys::kessler::run_kessler(
      v_t, v_qv, v_qc, v_qr, v_rho, v_pii, v_dz8w, v_z,
      v_rainnc, v_rainncv,
      real(dt_in), range, nz,
      real(c_xlv), real(c_cp), real(c_ep2),
      real(c_svp1), real(c_svp2), real(c_svp3),
      real(c_svpt0), real(c_rhowater));
  dev_sync();

  if (bench_mode) {
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < bench_iters; ++iter) {
      phys::kessler::run_kessler(
          v_t, v_qv, v_qc, v_qr, v_rho, v_pii, v_dz8w, v_z,
          v_rainnc, v_rainncv,
          real(dt_in), range, nz,
          real(c_xlv), real(c_cp), real(c_ep2),
          real(c_svp1), real(c_svp2), real(c_svp3),
          real(c_svpt0), real(c_rhowater));
    }
    dev_sync();
    auto t1 = std::chrono::high_resolution_clock::now();
    double total_sec = std::chrono::duration<double>(t1 - t0).count();
    double ns_per_col = (total_sec * 1e9) / (double(bench_iters) * double(n_cols));
    std::cout << "Benchmark (" << backend_name() << "): "
              << n_cols << " columns x " << nz << " levels, "
              << bench_iters << " iterations: "
              << total_sec * 1e3 / bench_iters << " ms/iter, "
              << ns_per_col << " ns/column" << std::endl;
  }

  // Copy outputs back to host
  copy_to_host(r_t.data(), d_t.get(), n3d * sizeof(real));
  copy_to_host(r_qv.data(), d_qv.get(), n3d * sizeof(real));
  copy_to_host(r_qc.data(), d_qc.get(), n3d * sizeof(real));
  copy_to_host(r_qr.data(), d_qr.get(), n3d * sizeof(real));
  copy_to_host(r_rainnc.data(), d_rainnc.get(), n2d * sizeof(real));
  copy_to_host(r_rainncv.data(), d_rainncv.get(), n2d * sizeof(real));

  // Convert to double for binary output
  for (std::size_t i = 0; i < n3d; ++i) {
    h_t[i] = double(r_t[i]);
    h_qv[i] = double(r_qv[i]);
    h_qc[i] = double(r_qc[i]);
    h_qr[i] = double(r_qr[i]);
  }
  for (std::size_t i = 0; i < n2d; ++i) {
    h_rainnc[i] = double(r_rainnc[i]);
    h_rainncv[i] = double(r_rainncv[i]);
  }

  FILE* fo = std::fopen(outfile.c_str(), "wb");
  if (!fo) {
    std::fprintf(stderr, "Cannot open output file: %s\n", outfile.c_str());
    return 1;
  }
  std::fwrite(&n_cols, sizeof(int32_t), 1, fo);
  std::fwrite(&nz, sizeof(int32_t), 1, fo);
  std::fwrite(h_t.data(), sizeof(double), n3d, fo);
  std::fwrite(h_qv.data(), sizeof(double), n3d, fo);
  std::fwrite(h_qc.data(), sizeof(double), n3d, fo);
  std::fwrite(h_qr.data(), sizeof(double), n3d, fo);
  std::fwrite(h_rainnc.data(), sizeof(double), n2d, fo);
  std::fwrite(h_rainncv.data(), sizeof(double), n2d, fo);
  std::fclose(fo);

  return 0;
}
