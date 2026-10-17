#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

#include "io_netcdf/wrf_readers.hpp"

using namespace cppwrf::io_netcdf;

static void check_near(double val, double expected, double tol, const std::string& name) {
  if (std::abs(val - expected) > tol) {
    std::cerr << "FAIL: " << name << " mismatch: got " << val << ", expected " << expected
              << " (diff=" << std::abs(val - expected) << ", tol=" << tol << ")\n";
    std::exit(1);
  }
}

int main() {
  const char* runs_env = std::getenv("CPPWRF_WRF_RUNS");
  if (!runs_env || !*runs_env) {
    std::cout << "Notice: CPPWRF_WRF_RUNS is not set; skipping real file checks.\n";
    return 0;
  }
  std::string runs = runs_env;
  std::string wrfinput_path = runs + "/wrf_run/wrfinput_d01";
  std::string wrfbdy_path = runs + "/wrf_run/wrfbdy_d01";

  if (!std::filesystem::exists(wrfinput_path)) {
    std::cout << "Notice: " << wrfinput_path << " not found, skipping real file checks.\n";
    return 0;
  }

  std::cout << "Testing read_wrfinput on " << wrfinput_path << "...\n";
  WrfInputData in = read_wrfinput(wrfinput_path);

  // 1. Grid dimensions
  assert(in.nx == 119);
  assert(in.ny == 89);
  assert(in.nz == 44);
  assert(in.nx_stag == 120);
  assert(in.ny_stag == 90);
  assert(in.nz_stag == 45);
  std::cout << "  PASS: Grid dimensions (119x89x44)\n";

  // 2. Global attributes
  check_near(in.dx, 15000.0, 1e-3, "DX");
  check_near(in.dy, 15000.0, 1e-3, "DY");
  check_near(in.dt, 60.0, 1e-3, "DT");
  assert(in.hybrid_opt == 2);
  check_near(in.etac, 0.20, 1e-4, "ETAC");
  assert(in.map_proj == 1);
  check_near(in.cen_lat, 39.0, 0.1, "CEN_LAT");
  check_near(in.cen_lon, 35.0, 0.1, "CEN_LON");
  check_near(in.p_top, 5000.0, 1e-3, "P_TOP");
  std::cout << "  PASS: Global attributes & resolution\n";

  // 3. Finite-difference constants & 1D coordinate arrays
  check_near(in.cf1, 1.828658, 1e-4, "CF1");
  check_near(in.cf2, -1.135213, 1e-4, "CF2");
  check_near(in.cf3, 0.306555, 1e-4, "CF3");
  check_near(in.cfn, 1.478440, 1e-4, "CFN");
  check_near(in.cfn1, -0.478440, 1e-4, "CFN1");

  assert(in.znw.size() == 45);
  check_near(in.znw.front(), 1.0, 1e-6, "ZNW[0]");
  check_near(in.znw.back(), 0.0, 1e-6, "ZNW[top]");

  assert(in.znu.size() == 44);
  check_near(in.znu.front(), 0.996907, 1e-4, "ZNU[0]");
  check_near(in.znu.back(), 0.002371, 1e-4, "ZNU[top]");

  assert(in.dn.size() == 44);
  assert(in.dnw.size() == 44);
  assert(in.rdn.size() == 44);
  assert(in.rdnw.size() == 44);
  assert(in.fnm.size() == 44);
  assert(in.fnp.size() == 44);
  assert(in.c1h.size() == 44);
  assert(in.c2h.size() == 44);
  assert(in.c1f.size() == 45);
  assert(in.c2f.size() == 45);
  check_near(in.c1h.front(), 1.013471, 1e-4, "C1H[0]");
  check_near(in.c2h.front(), -1279.82, 1.0, "C2H[0]");
  std::cout << "  PASS: 1D vertical and hybrid coordinate fields\n";

  // 4. 2D Horizontal fields
  assert(in.mub.size() == 89 * 119);
  assert(in.mu.size() == 89 * 119);
  assert(in.hgt.size() == 89 * 119);
  assert(in.xlat.size() == 89 * 119);
  assert(in.xlong.size() == 89 * 119);
  assert(in.mapfac_mx.size() == 89 * 119);
  assert(in.mapfac_my.size() == 89 * 119);
  assert(in.mapfac_ux.size() == 89 * 120);
  assert(in.mapfac_uy.size() == 89 * 120);
  assert(in.mapfac_vx.size() == 90 * 119);
  assert(in.mapfac_vy.size() == 90 * 119);
  assert(in.f.size() == 89 * 119);
  assert(in.e.size() == 89 * 119);
  assert(in.sinalpha.size() == 89 * 119);
  assert(in.cosalpha.size() == 89 * 119);

  check_near(in.mub[0], 95000.0, 1.0, "MUB[0,0]");
  check_near(in.mu[0], 1093.78, 1.0, "MU[0,0]");
  check_near(in.mapfac_mx[0], 0.99016, 1e-4, "MAPFAC_MX[0,0]");
  check_near(in.mapfac_ux[0], 0.99019, 1e-4, "MAPFAC_UX[0,0]");
  check_near(in.mapfac_vx[0], 0.990408, 1e-4, "MAPFAC_VX[0,0]");
  check_near(in.f[0], 7.82701e-05, 1e-8, "F[0,0]");
  check_near(in.e[0], 0.000123057, 1e-8, "E[0,0]");
  check_near(in.sinalpha[0], 0.118992, 1e-4, "SINALPHA[0,0]");
  check_near(in.cosalpha[0], 0.992895, 1e-4, "COSALPHA[0,0]");
  std::cout << "  PASS: 2D horizontal fields and map factors\n";

  // 5. 3D Fields & C++ indexing
  assert(in.pb.size() == 44 * 89 * 119);
  assert(in.p.size() == 44 * 89 * 119);
  assert(in.t.size() == 44 * 89 * 119);
  assert(in.t_init.size() == 44 * 89 * 119);
  assert(in.alb.size() == 44 * 89 * 119);
  assert(in.al.size() == 44 * 89 * 119);
  assert(in.qvapor.size() == 44 * 89 * 119);
  assert(in.u.size() == 44 * 89 * 120);
  assert(in.v.size() == 44 * 90 * 119);
  assert(in.w.size() == 45 * 89 * 119);
  assert(in.phb.size() == 45 * 89 * 119);
  assert(in.ph.size() == 45 * 89 * 119);

  check_near(in.pb[in.idx_m3(0, 0, 0)], 99706.2, 1.0, "PB[k=0,j=0,i=0]");
  check_near(in.t[in.idx_m3(0, 0, 0)], -3.75198, 1e-3, "T[k=0,j=0,i=0]");
  check_near(in.t_init[in.idx_m3(0, 0, 0)], -9.90335, 1e-3, "T_INIT[k=0,j=0,i=0]");
  check_near(in.alb[in.idx_m3(0, 0, 0)], 0.834329, 1e-4, "ALB[k=0,j=0,i=0]");
  check_near(in.al[in.idx_m3(0, 0, 0)], 0.024764, 1e-4, "AL[k=0,j=0,i=0]");
  check_near(in.p[in.idx_m3(0, 0, 0)], 1300.78, 1.0, "P[k=0,j=0,i=0]");
  check_near(in.u[in.idx_u3(0, 0, 0)], 2.17271, 1e-3, "U[k=0,j=0,i=0]");
  check_near(in.v[in.idx_v3(0, 0, 0)], 0.593197, 1e-3, "V[k=0,j=0,i=0]");
  check_near(in.qvapor[in.idx_m3(0, 0, 0)], 0.0109924, 1e-5, "QVAPOR[k=0,j=0,i=0]");
  std::cout << "  PASS: 3D fields and indexing match reference\n";

  // 6. Test read_wrfbdy
  if (std::filesystem::exists(wrfbdy_path)) {
    std::cout << "Testing read_wrfbdy on " << wrfbdy_path << "...\n";
    WrfBdyData bdy = read_wrfbdy(wrfbdy_path);

    assert(bdy.n_times == 1);
    assert(bdy.bdy_width == 5);
    assert(bdy.nx == 119);
    assert(bdy.ny == 89);
    assert(bdy.nz == 44);

    assert(!bdy.times.empty());
    assert(bdy.times[0] == "2026-10-12_00:00:00");
    std::cout << "  PASS: Boundary dimensions and times (" << bdy.times[0] << ")\n";

    // Check U boundary fields
    assert(bdy.u_bxs.size() == 1 * 5 * 44 * 89);
    assert(bdy.u_btxs.size() == 1 * 5 * 44 * 89);
    assert(bdy.u_bys.size() == 1 * 5 * 44 * 120);
    assert(bdy.u_btys.size() == 1 * 5 * 44 * 120);
    check_near(bdy.u_bxs[0], 210885.0, 10.0, "U_BXS[0]");
    check_near(bdy.u_btxs[0], 4.96951, 1e-3, "U_BTXS[0]");

    // Check V boundary fields
    assert(bdy.v_bxs.size() == 1 * 5 * 44 * 90);
    assert(bdy.v_bys.size() == 1 * 5 * 44 * 119);
    check_near(bdy.v_bxs[0], 57563.4, 10.0, "V_BXS[0]");

    // Check MU boundary fields (2D)
    assert(bdy.mu_bxs.size() == 1 * 5 * 89);
    check_near(bdy.mu_bxs[0], 1093.78, 1.0, "MU_BXS[0]");

    // Check T and QVAPOR
    assert(bdy.t_bxs.size() == 1 * 5 * 44 * 89);
    check_near(bdy.t_bxs[0], 142778.0, 10.0, "T_BXS[0]");
    assert(bdy.qvapor_bxs.size() == 1 * 5 * 44 * 89);
    check_near(bdy.qvapor_bxs[0], 1056.46, 1.0, "QVAPOR_BXS[0]");
    std::cout << "  PASS: Boundary fields spot-checked\n";
  }

  std::cout << "ALL WRF READER TESTS PASSED SUCCESSFULLY!\n";
  return 0;
}
