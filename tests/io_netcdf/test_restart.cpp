#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

#include "io_netcdf/netcdf_common.hpp"
#include "io_netcdf/restart.hpp"

using namespace cppwrf::io_netcdf;

static void test_precision(DType dt) {
  bool is_dp = (dt == DType::F64);
  std::printf("Testing restart bitwise identity for %s...\n", is_dp ? "FP64" : "FP32");

  const int nx = 12, ny = 8, nz = 6;
  const size_t sz_u = (nx + 1) * ny * nz;
  const size_t sz_v = nx * (ny + 1) * nz;
  const size_t sz_w = nx * ny * (nz + 1);
  const size_t sz_m = nx * ny * nz;
  const size_t sz_mu = nx * ny;

  RestartMeta meta_out;
  meta_out.step = 300;
  meta_out.time_sec = 300.0;
  meta_out.nx = nx;
  meta_out.ny = ny;
  meta_out.nz = nz;
  meta_out.dx = 100.0;
  meta_out.dy = 100.0;
  meta_out.dt = 1.0;
  meta_out.case_name = "straka";
  meta_out.date_str = "0001-01-01_00:05:00";

  // Lambda to fill test pattern
  auto fill = [is_dp](void* buf, size_t n, double seed) {
    if (is_dp) {
      double* p = static_cast<double*>(buf);
      for (size_t i = 0; i < n; ++i) {
        p[i] = std::sin(seed + double(i) * 0.123456789) * 1.23456789e4 + std::cos(double(i)) * 1e-7;
      }
    } else {
      float* p = static_cast<float*>(buf);
      for (size_t i = 0; i < n; ++i) {
        p[i] = float(std::sin(seed + double(i) * 0.12345) * 1.234e4 + std::cos(double(i)) * 1e-4);
      }
    }
  };

  const size_t item_size = is_dp ? sizeof(double) : sizeof(float);

  // Allocate original state buffers
  std::vector<char> u1(sz_u * item_size), u2(sz_u * item_size);
  std::vector<char> v1(sz_v * item_size), v2(sz_v * item_size);
  std::vector<char> w1(sz_w * item_size), w2(sz_w * item_size);
  std::vector<char> t1(sz_m * item_size), t2(sz_m * item_size);
  std::vector<char> ph1(sz_w * item_size), ph2(sz_w * item_size);
  std::vector<char> mu1(sz_mu * item_size), mu2(sz_mu * item_size);
  std::vector<char> q1(sz_m * item_size), q2(sz_m * item_size);

  fill(u1.data(), sz_u, 1.1); fill(u2.data(), sz_u, 1.2);
  fill(v1.data(), sz_v, 2.1); fill(v2.data(), sz_v, 2.2);
  fill(w1.data(), sz_w, 3.1); fill(w2.data(), sz_w, 3.2);
  fill(t1.data(), sz_m, 4.1); fill(t2.data(), sz_m, 4.2);
  fill(ph1.data(), sz_w, 5.1); fill(ph2.data(), sz_w, 5.2);
  fill(mu1.data(), sz_mu, 6.1); fill(mu2.data(), sz_mu, 6.2);
  fill(q1.data(), sz_m, 7.1); fill(q2.data(), sz_m, 7.2);

  std::vector<VarView> vars_out;
  vars_out.emplace_back("u_1", 'U', nx + 1, ny, nz, u1.data(), dt);
  vars_out.emplace_back("u_2", 'U', nx + 1, ny, nz, u2.data(), dt);
  vars_out.emplace_back("v_1", 'V', nx, ny + 1, nz, v1.data(), dt);
  vars_out.emplace_back("v_2", 'V', nx, ny + 1, nz, v2.data(), dt);
  vars_out.emplace_back("w_1", 'W', nx, ny, nz + 1, w1.data(), dt);
  vars_out.emplace_back("w_2", 'W', nx, ny, nz + 1, w2.data(), dt);
  vars_out.emplace_back("t_1", 'M', nx, ny, nz, t1.data(), dt);
  vars_out.emplace_back("t_2", 'M', nx, ny, nz, t2.data(), dt);
  vars_out.emplace_back("ph_1", 'W', nx, ny, nz + 1, ph1.data(), dt);
  vars_out.emplace_back("ph_2", 'W', nx, ny, nz + 1, ph2.data(), dt);
  vars_out.emplace_back("mu_1", 'M', nx, ny, 1, mu1.data(), dt);
  vars_out.emplace_back("mu_2", 'M', nx, ny, 1, mu2.data(), dt);
  vars_out.emplace_back("q1_qv", 'M', nx, ny, nz, q1.data(), dt);
  vars_out.emplace_back("q2_qv", 'M', nx, ny, nz, q2.data(), dt);

  std::string filename = is_dp ? "tests/io_netcdf/test_out/restart_dp.nc" : "tests/io_netcdf/test_out/restart_sp.nc";
  write_restart(filename, meta_out, vars_out);

  // Allocate fresh read buffers zeroed out
  std::vector<char> ru1(sz_u * item_size, 0), ru2(sz_u * item_size, 0);
  std::vector<char> rv1(sz_v * item_size, 0), rv2(sz_v * item_size, 0);
  std::vector<char> rw1(sz_w * item_size, 0), rw2(sz_w * item_size, 0);
  std::vector<char> rt1(sz_m * item_size, 0), rt2(sz_m * item_size, 0);
  std::vector<char> rph1(sz_w * item_size, 0), rph2(sz_w * item_size, 0);
  std::vector<char> rmu1(sz_mu * item_size, 0), rmu2(sz_mu * item_size, 0);
  std::vector<char> rq1(sz_m * item_size, 0), rq2(sz_m * item_size, 0);

  std::vector<VarView> vars_in;
  vars_in.emplace_back("u_1", 'U', nx + 1, ny, nz, ru1.data(), dt);
  vars_in.emplace_back("u_2", 'U', nx + 1, ny, nz, ru2.data(), dt);
  vars_in.emplace_back("v_1", 'V', nx, ny + 1, nz, rv1.data(), dt);
  vars_in.emplace_back("v_2", 'V', nx, ny + 1, nz, rv2.data(), dt);
  vars_in.emplace_back("w_1", 'W', nx, ny, nz + 1, rw1.data(), dt);
  vars_in.emplace_back("w_2", 'W', nx, ny, nz + 1, rw2.data(), dt);
  vars_in.emplace_back("t_1", 'M', nx, ny, nz, rt1.data(), dt);
  vars_in.emplace_back("t_2", 'M', nx, ny, nz, rt2.data(), dt);
  vars_in.emplace_back("ph_1", 'W', nx, ny, nz + 1, rph1.data(), dt);
  vars_in.emplace_back("ph_2", 'W', nx, ny, nz + 1, rph2.data(), dt);
  vars_in.emplace_back("mu_1", 'M', nx, ny, 1, rmu1.data(), dt);
  vars_in.emplace_back("mu_2", 'M', nx, ny, 1, rmu2.data(), dt);
  vars_in.emplace_back("q1_qv", 'M', nx, ny, nz, rq1.data(), dt);
  vars_in.emplace_back("q2_qv", 'M', nx, ny, nz, rq2.data(), dt);

  RestartMeta meta_in;
  read_restart(filename, meta_in, vars_in);

  // Validate metadata
  assert(meta_in.step == meta_out.step);
  assert(std::fabs(meta_in.time_sec - meta_out.time_sec) < 1e-12);
  assert(meta_in.nx == meta_out.nx);
  assert(meta_in.ny == meta_out.ny);
  assert(meta_in.nz == meta_out.nz);
  assert(std::fabs(meta_in.dx - meta_out.dx) < 1e-12);
  assert(std::fabs(meta_in.dy - meta_out.dy) < 1e-12);
  assert(std::fabs(meta_in.dt - meta_out.dt) < 1e-12);
  assert(meta_in.case_name == meta_out.case_name);
  assert(meta_in.date_str == meta_out.date_str);

  // Validate bitwise exact identity using memcmp
  assert(std::memcmp(u1.data(), ru1.data(), u1.size()) == 0);
  assert(std::memcmp(u2.data(), ru2.data(), u2.size()) == 0);
  assert(std::memcmp(v1.data(), rv1.data(), v1.size()) == 0);
  assert(std::memcmp(v2.data(), rv2.data(), v2.size()) == 0);
  assert(std::memcmp(w1.data(), rw1.data(), w1.size()) == 0);
  assert(std::memcmp(w2.data(), rw2.data(), w2.size()) == 0);
  assert(std::memcmp(t1.data(), rt1.data(), t1.size()) == 0);
  assert(std::memcmp(t2.data(), rt2.data(), t2.size()) == 0);
  assert(std::memcmp(ph1.data(), rph1.data(), ph1.size()) == 0);
  assert(std::memcmp(ph2.data(), rph2.data(), ph2.size()) == 0);
  assert(std::memcmp(mu1.data(), rmu1.data(), mu1.size()) == 0);
  assert(std::memcmp(mu2.data(), rmu2.data(), mu2.size()) == 0);
  assert(std::memcmp(q1.data(), rq1.data(), q1.size()) == 0);
  assert(std::memcmp(q2.data(), rq2.data(), q2.size()) == 0);

  std::printf("PASS: Bitwise identity verified for %s (memcmp == 0)\n", is_dp ? "FP64" : "FP32");
}

int main() {
  std::printf("Running test_restart acceptance test...\n");
  test_precision(DType::F32);
  test_precision(DType::F64);
  std::printf("RESTART BITWISE TEST PASSED!\n");
  return 0;
}
