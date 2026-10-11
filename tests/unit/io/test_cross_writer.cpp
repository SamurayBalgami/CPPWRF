#include "io/output.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

using namespace cppwrf::io;

int main(int argc, char* argv[]) {
  std::string out_dir = "/tmp/cppwrf_cross_test";
  if (argc > 1) {
    out_dir = argv[1];
  }

  std::filesystem::create_directories(out_dir);

  // 1. Write run.json
  std::string run_json_content = "{\"case\":\"straka\",\"nx\":4,\"ny\":3,\"nz\":2,\"dt\":0.5}";
  write_text_file_atomic((std::filesystem::path(out_dir) / "run.json").string(), run_json_content);

  // 2. Write diag.csv
  {
    DiagWriter diag((std::filesystem::path(out_dir) / "diag.csv").string(), {"step", "time", "energy", "temp"});
    diag.row({0.0, 0.0, 1000.5, 300.0});
    diag.row({100.0, 12.5, 999.8, 300.2});
    diag.row({200.0, 25.0, 998.2, 300.4});
  }

  // 3. Write snapshots with strided data
  const int nx = 4, ny = 3, nz = 2;
  const int halo_x = 2, halo_y = 2, halo_z = 1;
  const int alloc_nx = nx + 2 * halo_x;
  const int alloc_ny = ny + 2 * halo_y;
  const int alloc_nz = nz + 2 * halo_z;
  const std::ptrdiff_t sj = alloc_nx;
  const std::ptrdiff_t sk = alloc_nx * alloc_ny;

  std::vector<float> buf_u(sk * alloc_nz, -999.0f);
  float* origin_u = buf_u.data() + (halo_x + halo_y * sj + halo_z * sk);

  std::vector<double> buf_t(sk * alloc_nz, -999.0);
  double* origin_t = buf_t.data() + (halo_x + halo_y * sj + halo_z * sk);

  std::vector<float> buf_psfc(sj * alloc_ny, -999.0f);
  float* origin_psfc = buf_psfc.data() + (halo_x + halo_y * sj);

  SnapshotWriter snap_writer(out_dir);

  for (long step : {100L, 200L}) {
    double time = (step == 100L ? 12.5 : 25.0);

    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          origin_u[i + j * sj + k * sk] = static_cast<float>(step + i * 1.0f + j * 10.0f + k * 100.0f);
          origin_t[i + j * sj + k * sk] = static_cast<double>(step * 2.0 + i * 0.1 + j * 0.01 + k * 0.001);
        }
      }
    }
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        origin_psfc[i + j * sj] = static_cast<float>(step * 10.0f + i + j * 5.0f);
      }
    }

    std::vector<VarView> vars = {
      {"u", 'U', nx, ny, nz, origin_u, DType::F32, sj, sk},
      {"temp", 'M', nx, ny, nz, origin_t, DType::F64, sj, sk},
      {"psfc", 'M', nx, ny, 1, origin_psfc, DType::F32, sj, 0}
    };
    snap_writer.write(step, time, vars);
  }

  std::cout << "Successfully wrote test cross outputs to: " << out_dir << std::endl;
  return 0;
}
