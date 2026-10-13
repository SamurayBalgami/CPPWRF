#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "io_netcdf/netcdf_common.hpp"
#include "io_netcdf/wrfout_writer.hpp"

using namespace cppwrf::io_netcdf;

int main() {
  std::printf("Running test_wrfout_writer...\n");

  const int nx = 16, ny = 8, nz = 10;
  const double dx = 1000.0, dy = 1000.0, dt = 10.0;

  WrfoutMeta meta;
  meta.nx = nx;
  meta.ny = ny;
  meta.nz = nz;
  meta.dx = dx;
  meta.dy = dy;
  meta.dt = dt;
  meta.ztop = 10000.0;
  meta.p_top = 5000.0;
  meta.start_date = "0001-01-01_00:00:00";
  meta.simulation_start_date = "0001-01-01_00:00:00";
  meta.diff_opt = 1;
  meta.km_opt = 1;
  meta.khdif = 75.0;
  meta.kvdif = 75.0;

  // Allocate host arrays
  std::vector<float> u((nx + 1) * ny * nz, 10.5f);
  std::vector<float> v(nx * (ny + 1) * nz, -2.5f);
  std::vector<float> w(nx * ny * (nz + 1), 0.15f);
  std::vector<float> ph(nx * ny * (nz + 1), 980.0f);
  std::vector<float> phb(nx * ny * (nz + 1), 9800.0f);
  std::vector<float> t(nx * ny * nz, 3.5f); // theta - 300
  std::vector<float> p(nx * ny * nz, -50.0f);
  std::vector<float> pb(nx * ny * nz, 100000.0f);
  std::vector<float> qv(nx * ny * nz, 0.005f);
  std::vector<float> mu(nx * ny, 20.0f);
  std::vector<float> mub(nx * ny, 95000.0f);
  std::vector<float> xlat(nx * ny, 0.0f);
  std::vector<float> xlong(nx * ny, 0.0f);
  std::vector<float> hgt(nx * ny, 0.0f);

  std::vector<float> znu(nz);
  std::vector<float> znw(nz + 1);
  for (int k = 0; k <= nz; ++k) znw[k] = 1.0f - float(k) / float(nz);
  for (int k = 0; k < nz; ++k) znu[k] = 0.5f * (znw[k] + znw[k + 1]);

  std::vector<float> c1h(nz, 1.0f), c2h(nz, 0.0f);
  std::vector<float> c1f(nz + 1, 1.0f), c2f(nz + 1, 0.0f);
  std::vector<float> c3h(nz, 0.5f), c4h(nz, 0.0f);
  std::vector<float> c3f(nz + 1, 0.5f), c4f(nz + 1, 0.0f);

  float rdx = float(1.0 / dx);
  float rdy = float(1.0 / dy);
  float p_top = float(meta.p_top);

  std::vector<VarView> vars;
  vars.emplace_back("U", 'U', nx + 1, ny, nz, u.data(), DType::F32);
  vars.emplace_back("V", 'V', nx, ny + 1, nz, v.data(), DType::F32);
  vars.emplace_back("W", 'W', nx, ny, nz + 1, w.data(), DType::F32);
  vars.emplace_back("PH", 'W', nx, ny, nz + 1, ph.data(), DType::F32);
  vars.emplace_back("PHB", 'W', nx, ny, nz + 1, phb.data(), DType::F32);
  vars.emplace_back("T", 'M', nx, ny, nz, t.data(), DType::F32);
  vars.emplace_back("P", 'M', nx, ny, nz, p.data(), DType::F32);
  vars.emplace_back("PB", 'M', nx, ny, nz, pb.data(), DType::F32);
  vars.emplace_back("QVAPOR", 'M', nx, ny, nz, qv.data(), DType::F32);
  vars.emplace_back("MU", 'M', nx, ny, 1, mu.data(), DType::F32);
  vars.emplace_back("MUB", 'M', nx, ny, 1, mub.data(), DType::F32);
  vars.emplace_back("XLAT", 'M', nx, ny, 1, xlat.data(), DType::F32);
  vars.emplace_back("XLONG", 'M', nx, ny, 1, xlong.data(), DType::F32);
  vars.emplace_back("HGT", 'M', nx, ny, 1, hgt.data(), DType::F32);
  vars.emplace_back("ZNU", '\0', 1, 1, nz, znu.data(), DType::F32);
  vars.emplace_back("ZNW", 'W', 1, 1, nz + 1, znw.data(), DType::F32);
  vars.emplace_back("C1H", '\0', 1, 1, nz, c1h.data(), DType::F32);
  vars.emplace_back("C2H", '\0', 1, 1, nz, c2h.data(), DType::F32);
  vars.emplace_back("C1F", 'W', 1, 1, nz + 1, c1f.data(), DType::F32);
  vars.emplace_back("C2F", 'W', 1, 1, nz + 1, c2f.data(), DType::F32);
  vars.emplace_back("C3H", '\0', 1, 1, nz, c3h.data(), DType::F32);
  vars.emplace_back("C4H", '\0', 1, 1, nz, c4h.data(), DType::F32);
  vars.emplace_back("C3F", 'W', 1, 1, nz + 1, c3f.data(), DType::F32);
  vars.emplace_back("C4F", 'W', 1, 1, nz + 1, c4f.data(), DType::F32);
  vars.emplace_back("RDX", '\0', 1, 1, 1, &rdx, DType::F32);
  vars.emplace_back("RDY", '\0', 1, 1, 1, &rdy, DType::F32);
  vars.emplace_back("P_TOP", '\0', 1, 1, 1, &p_top, DType::F32);

  std::string test_file = "tests/io_netcdf/test_out/wrfout_d01_test.nc";
  std::remove(test_file.c_str());
  WrfoutWriter writer("tests/io_netcdf/test_out", meta);

  // Write frame 0 (t = 0 s)
  writer.write_frame(0, 0.0, vars, test_file);

  // Update a few values and write frame 1 (t = 300 s)
  u[0] = 12.0f;
  t[0] = 4.2f;
  writer.write_frame(30, 300.0, vars, test_file);

  // Verify file with NetCDF C API
  int ncid = -1;
  NC_CHECK(nc_open(test_file.c_str(), NC_NOWRITE, &ncid));

  // Check dimensions
  int dim_time = -1, dim_we = -1, dim_sn = -1, dim_bt = -1;
  NC_CHECK(nc_inq_dimid(ncid, "Time", &dim_time));
  NC_CHECK(nc_inq_dimid(ncid, "west_east", &dim_we));
  NC_CHECK(nc_inq_dimid(ncid, "south_north", &dim_sn));
  NC_CHECK(nc_inq_dimid(ncid, "bottom_top", &dim_bt));

  size_t time_len = 0, we_len = 0, sn_len = 0, bt_len = 0;
  NC_CHECK(nc_inq_dimlen(ncid, dim_time, &time_len));
  NC_CHECK(nc_inq_dimlen(ncid, dim_we, &we_len));
  NC_CHECK(nc_inq_dimlen(ncid, dim_sn, &sn_len));
  NC_CHECK(nc_inq_dimlen(ncid, dim_bt, &bt_len));

  assert(time_len == 2);
  assert(int(we_len) == nx);
  assert(int(sn_len) == ny);
  assert(int(bt_len) == nz);

  // Check variable U
  int var_u = -1;
  NC_CHECK(nc_inq_varid(ncid, "U", &var_u));
  char mem_order[4] = {0};
  NC_CHECK(nc_get_att_text(ncid, var_u, "MemoryOrder", mem_order));
  assert(std::string(mem_order, 3) == "XYZ");

  // Read back U frame 0 and frame 1
  std::vector<float> u_read((nx + 1) * ny * nz);
  size_t start[4] = {0, 0, 0, 0};
  size_t count[4] = {1, size_t(nz), size_t(ny), size_t(nx + 1)};
  NC_CHECK(nc_get_vara_float(ncid, var_u, start, count, u_read.data()));
  assert(std::fabs(u_read[0] - 10.5f) < 1e-6);

  start[0] = 1;
  NC_CHECK(nc_get_vara_float(ncid, var_u, start, count, u_read.data()));
  assert(std::fabs(u_read[0] - 12.0f) < 1e-6);

  // Check Times variable
  int var_times = -1;
  NC_CHECK(nc_inq_varid(ncid, "Times", &var_times));
  char tstr[20] = {0};
  size_t tstart[2] = {0, 0};
  size_t tcount[2] = {1, 19};
  NC_CHECK(nc_get_vara_text(ncid, var_times, tstart, tcount, tstr));
  assert(std::string(tstr, 19) == "0001-01-01_00:00:00");

  tstart[0] = 1;
  NC_CHECK(nc_get_vara_text(ncid, var_times, tstart, tcount, tstr));
  assert(std::string(tstr, 19) == "0001-01-01_00:05:00");

  NC_CHECK(nc_close(ncid));
  std::printf("test_wrfout_writer PASSED!\n");
  return 0;
}
