#include "phys/ra_rrtmg_sw/rrtmg_sw.hpp"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace cppwrf::phys;

int main(int argc, char** argv) {
  const char* infile = (argc > 1) ? argv[1] : "input_rrtmg_sw.bin";
  const char* outfile = (argc > 2) ? argv[2] : "output_cpp.bin";
  const char* tables_file = (argc > 3) ? argv[3] : "rrtmg_sw_tables.bin";

  FILE* fin = std::fopen(infile, "rb");
  if (!fin) {
    std::cerr << "Failed to open input file: " << infile << std::endl;
    return 1;
  }

  int num_cols = 0;
  int nz = 0;

  if (std::fread(&num_cols, sizeof(int), 1, fin) != 1 ||
      std::fread(&nz, sizeof(int), 1, fin) != 1) {
    std::cerr << "Failed to read header from: " << infile << std::endl;
    std::fclose(fin);
    return 1;
  }

  std::cout << "Loading inputs: num_cols = " << num_cols << ", nz = " << nz << std::endl;

  // 2D inputs
  std::vector<double> coszen(num_cols);
  std::vector<double> albedo(num_cols);
  std::vector<double> xland(num_cols);
  std::vector<double> xice(num_cols);
  std::vector<double> snow(num_cols);
  std::vector<double> xlat(num_cols);

  std::fread(coszen.data(), sizeof(double), num_cols, fin);
  std::fread(albedo.data(), sizeof(double), num_cols, fin);
  std::fread(xland.data(), sizeof(double), num_cols, fin);
  std::fread(xice.data(), sizeof(double), num_cols, fin);
  std::fread(snow.data(), sizeof(double), num_cols, fin);
  std::fread(xlat.data(), sizeof(double), num_cols, fin);

  // 3D layer inputs: (nz, num_cols)
  size_t n_layer = static_cast<size_t>(nz) * num_cols;
  std::vector<double> p3d(n_layer);
  std::vector<double> pi3d(n_layer);
  std::vector<double> dz8w(n_layer);
  std::vector<double> t3d(n_layer);
  std::vector<double> rho3d(n_layer);
  std::vector<double> qv3d(n_layer);
  std::vector<double> qc3d(n_layer);
  std::vector<double> qr3d(n_layer);
  std::vector<double> qi3d(n_layer);
  std::vector<double> qs3d(n_layer);
  std::vector<double> qg3d(n_layer);
  std::vector<double> cldfra3d(n_layer);

  std::fread(p3d.data(), sizeof(double), n_layer, fin);
  std::fread(pi3d.data(), sizeof(double), n_layer, fin);
  std::fread(dz8w.data(), sizeof(double), n_layer, fin);
  std::fread(t3d.data(), sizeof(double), n_layer, fin);
  std::fread(rho3d.data(), sizeof(double), n_layer, fin);
  std::fread(qv3d.data(), sizeof(double), n_layer, fin);
  std::fread(qc3d.data(), sizeof(double), n_layer, fin);
  std::fread(qr3d.data(), sizeof(double), n_layer, fin);
  std::fread(qi3d.data(), sizeof(double), n_layer, fin);
  std::fread(qs3d.data(), sizeof(double), n_layer, fin);
  std::fread(qg3d.data(), sizeof(double), n_layer, fin);
  std::fread(cldfra3d.data(), sizeof(double), n_layer, fin);

  // 3D interface inputs: (nz + 1, num_cols)
  size_t n_iface = static_cast<size_t>(nz + 1) * num_cols;
  std::vector<double> p8w(n_iface);
  std::vector<double> t8w(n_iface);

  std::fread(p8w.data(), sizeof(double), n_iface, fin);
  std::fread(t8w.data(), sizeof(double), n_iface, fin);

  std::fclose(fin);

  // Outputs
  std::vector<double> gsw(num_cols, 0.0);
  std::vector<double> swcf(num_cols, 0.0);
  std::vector<double> swupt(num_cols, 0.0);
  std::vector<double> swuptc(num_cols, 0.0);
  std::vector<double> swdnt(num_cols, 0.0);
  std::vector<double> swdntc(num_cols, 0.0);
  std::vector<double> swupb(num_cols, 0.0);
  std::vector<double> swupbc(num_cols, 0.0);
  std::vector<double> swdnb(num_cols, 0.0);
  std::vector<double> swdnbc(num_cols, 0.0);
  std::vector<double> swvisdir(num_cols, 0.0);
  std::vector<double> swvisdif(num_cols, 0.0);
  std::vector<double> swnirdir(num_cols, 0.0);
  std::vector<double> swnirdif(num_cols, 0.0);
  std::vector<double> rthratensw(n_layer, 0.0);
  std::vector<double> rthratenswc(n_layer, 0.0);

  // Load tables
  RrtmgSwTables tables;
  if (!load_rrtmg_sw_tables(tables_file, tables)) {
    if (!load_rrtmg_sw_tables("../rrtmg_sw_tables.bin", tables)) {
      std::cerr << "Failed to load RRTMG SW tables from: " << tables_file << std::endl;
      return 1;
    }
  }

  RrtmgSwConfig config;

  std::cout << "Running C++ RRTMG SWRAD..." << std::endl;
  auto t_start = std::chrono::high_resolution_clock::now();

  rrtmg_swrad(
      num_cols, nz, config, tables,
      coszen.data(), albedo.data(), xland.data(), xice.data(), snow.data(), xlat.data(),
      p8w.data(), t8w.data(), p3d.data(), pi3d.data(), dz8w.data(), t3d.data(), rho3d.data(),
      qv3d.data(), qc3d.data(), qr3d.data(), qi3d.data(), qs3d.data(), qg3d.data(), cldfra3d.data(),
      gsw.data(), swcf.data(), swupt.data(), swuptc.data(), swdnt.data(), swdntc.data(),
      swupb.data(), swupbc.data(), swdnb.data(), swdnbc.data(),
      swvisdir.data(), swvisdif.data(), swnirdir.data(), swnirdif.data(),
      rthratensw.data(), rthratenswc.data()
  );

  auto t_end = std::chrono::high_resolution_clock::now();
  double elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
  std::cout << "Finished RRTMG SWRAD in " << elapsed_ms << " ms ("
            << (num_cols / (elapsed_ms * 1e-3)) << " col/s)." << std::endl;
  std::cout << "Sample Col 0: GSW = " << gsw[0] << ", SWCF = " << swcf[0]
            << ", SWDNB = " << swdnb[0] << ", SWUPT = " << swupt[0] << std::endl;

  // Write outputs
  FILE* fout = std::fopen(outfile, "wb");
  if (!fout) {
    std::cerr << "Failed to open output file: " << outfile << std::endl;
    return 1;
  }

  std::fwrite(gsw.data(), sizeof(double), num_cols, fout);
  std::fwrite(swcf.data(), sizeof(double), num_cols, fout);
  std::fwrite(swupt.data(), sizeof(double), num_cols, fout);
  std::fwrite(swuptc.data(), sizeof(double), num_cols, fout);
  std::fwrite(swdnt.data(), sizeof(double), num_cols, fout);
  std::fwrite(swdntc.data(), sizeof(double), num_cols, fout);
  std::fwrite(swupb.data(), sizeof(double), num_cols, fout);
  std::fwrite(swupbc.data(), sizeof(double), num_cols, fout);
  std::fwrite(swdnb.data(), sizeof(double), num_cols, fout);
  std::fwrite(swdnbc.data(), sizeof(double), num_cols, fout);
  std::fwrite(swvisdir.data(), sizeof(double), num_cols, fout);
  std::fwrite(swvisdif.data(), sizeof(double), num_cols, fout);
  std::fwrite(swnirdir.data(), sizeof(double), num_cols, fout);
  std::fwrite(swnirdif.data(), sizeof(double), num_cols, fout);

  std::fwrite(rthratensw.data(), sizeof(double), n_layer, fout);
  std::fwrite(rthratenswc.data(), sizeof(double), n_layer, fout);

  std::fclose(fout);
  std::cout << "Successfully wrote C++ RRTMG-SW outputs to: " << outfile << std::endl;

  return 0;
}
