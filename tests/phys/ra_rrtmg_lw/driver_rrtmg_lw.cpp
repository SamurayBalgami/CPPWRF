#include "phys/ra_rrtmg_lw/rrtmg_lw.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

using namespace cppwrf::phys;

int main(int argc, char** argv) {
  const char* infile = "input_rrtmg_lw.bin";
  const char* outfile = "output_cpp.bin";
  const char* tables_file = "rrtmg_lw_tables.bin";

  if (argc > 1) infile = argv[1];
  if (argc > 2) outfile = argv[2];
  if (argc > 3) tables_file = argv[3];

  FILE* fin = std::fopen(infile, "rb");
  if (!fin) {
    std::cerr << "Failed to open input file: " << infile << std::endl;
    return 1;
  }

  uint32_t magic;
  int32_t num_cols, nz;
  double p_top;
  int32_t icloud, cldovrlp, idcor, mp_physics;

  if (std::fread(&magic, sizeof(uint32_t), 1, fin) != 1 ||
      std::fread(&num_cols, sizeof(int32_t), 1, fin) != 1 ||
      std::fread(&nz, sizeof(int32_t), 1, fin) != 1 ||
      std::fread(&p_top, sizeof(double), 1, fin) != 1 ||
      std::fread(&icloud, sizeof(int32_t), 1, fin) != 1 ||
      std::fread(&cldovrlp, sizeof(int32_t), 1, fin) != 1 ||
      std::fread(&idcor, sizeof(int32_t), 1, fin) != 1 ||
      std::fread(&mp_physics, sizeof(int32_t), 1, fin) != 1) {
    std::cerr << "Failed to read header from " << infile << std::endl;
    std::fclose(fin);
    return 1;
  }

  std::cout << "Header: num_cols = " << num_cols << ", nz = " << nz
            << ", p_top = " << p_top << ", icloud = " << icloud
            << ", cldovrlp = " << cldovrlp << std::endl;

  // 2D inputs
  std::vector<double> tsk(num_cols);
  std::vector<double> emiss(num_cols);
  std::vector<double> xland(num_cols);
  std::vector<double> xice(num_cols);
  std::vector<double> snow(num_cols);
  std::vector<double> xlat(num_cols);

  std::fread(tsk.data(), sizeof(double), num_cols, fin);
  std::fread(emiss.data(), sizeof(double), num_cols, fin);
  std::fread(xland.data(), sizeof(double), num_cols, fin);
  std::fread(xice.data(), sizeof(double), num_cols, fin);
  std::fread(snow.data(), sizeof(double), num_cols, fin);
  std::fread(xlat.data(), sizeof(double), num_cols, fin);

  // 3D interface inputs: (nz + 1, num_cols)
  size_t n_iface = static_cast<size_t>(nz + 1) * num_cols;
  std::vector<double> p8w(n_iface);
  std::vector<double> t8w(n_iface);

  std::fread(p8w.data(), sizeof(double), n_iface, fin);
  std::fread(t8w.data(), sizeof(double), n_iface, fin);

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

  std::fclose(fin);

  // Outputs
  std::vector<double> glw(num_cols, 0.0);
  std::vector<double> olr(num_cols, 0.0);
  std::vector<double> lwcf(num_cols, 0.0);
  std::vector<double> lwupt(num_cols, 0.0);
  std::vector<double> lwuptc(num_cols, 0.0);
  std::vector<double> lwdnt(num_cols, 0.0);
  std::vector<double> lwdntc(num_cols, 0.0);
  std::vector<double> lwupb(num_cols, 0.0);
  std::vector<double> lwupbc(num_cols, 0.0);
  std::vector<double> lwdnb(num_cols, 0.0);
  std::vector<double> lwdnbc(num_cols, 0.0);
  std::vector<double> rthratenlw(n_layer, 0.0);
  std::vector<double> rthratenlwc(n_layer, 0.0);

  // Load tables
  RrtmgLwTables tables;
  if (!load_rrtmg_lw_tables(tables_file, tables)) {
    // Try in parent directory if needed
    if (!load_rrtmg_lw_tables("../rrtmg_lw_tables.bin", tables)) {
      std::cerr << "Failed to load RRTMG LW tables from: " << tables_file << std::endl;
      return 1;
    }
  }

  RrtmgLwConfig config;
  config.p_top = p_top;
  config.icloud = icloud;
  config.cldovrlp = cldovrlp;
  config.idcor = idcor;
  config.mp_physics = mp_physics;

  std::cout << "Running C++ RRTMG LWRAD..." << std::endl;
  auto t_start = std::chrono::high_resolution_clock::now();

  rrtmg_lwrad(
      num_cols, nz, config, tables,
      tsk.data(), emiss.data(), xland.data(), xice.data(), snow.data(), xlat.data(),
      p8w.data(), t8w.data(), p3d.data(), pi3d.data(), dz8w.data(), t3d.data(), rho3d.data(),
      qv3d.data(), qc3d.data(), qr3d.data(), qi3d.data(), qs3d.data(), qg3d.data(), cldfra3d.data(),
      glw.data(), olr.data(), lwcf.data(),
      lwupt.data(), lwuptc.data(), lwdnt.data(), lwdntc.data(),
      lwupb.data(), lwupbc.data(), lwdnb.data(), lwdnbc.data(),
      rthratenlw.data(), rthratenlwc.data()
  );

  auto t_end = std::chrono::high_resolution_clock::now();
  double elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
  std::cout << "Finished RRTMG LWRAD in " << elapsed_ms << " ms ("
            << (num_cols / (elapsed_ms * 1e-3)) << " col/s)." << std::endl;
  std::cout << "Sample Col 0: GLW = " << glw[0] << ", OLR = " << olr[0] << ", LWCF = " << lwcf[0] << std::endl;

  // Write outputs
  FILE* fout = std::fopen(outfile, "wb");
  if (!fout) {
    std::cerr << "Failed to open output file: " << outfile << std::endl;
    return 1;
  }

  std::fwrite(&num_cols, sizeof(int32_t), 1, fout);
  std::fwrite(&nz, sizeof(int32_t), 1, fout);

  std::fwrite(rthratenlw.data(), sizeof(double), n_layer, fout);
  std::fwrite(rthratenlwc.data(), sizeof(double), n_layer, fout);

  std::fwrite(glw.data(), sizeof(double), num_cols, fout);
  std::fwrite(olr.data(), sizeof(double), num_cols, fout);
  std::fwrite(lwcf.data(), sizeof(double), num_cols, fout);
  std::fwrite(lwupt.data(), sizeof(double), num_cols, fout);
  std::fwrite(lwuptc.data(), sizeof(double), num_cols, fout);
  std::fwrite(lwdnt.data(), sizeof(double), num_cols, fout);
  std::fwrite(lwdntc.data(), sizeof(double), num_cols, fout);
  std::fwrite(lwupb.data(), sizeof(double), num_cols, fout);
  std::fwrite(lwupbc.data(), sizeof(double), num_cols, fout);
  std::fwrite(lwdnb.data(), sizeof(double), num_cols, fout);
  std::fwrite(lwdnbc.data(), sizeof(double), num_cols, fout);

  std::fclose(fout);
  std::cout << "Output successfully written to " << outfile << std::endl;

  return 0;
}
