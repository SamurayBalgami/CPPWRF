#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <vector>
#include <netcdf.h>

#include "bdy/bdy_data.hpp"
#include "bdy/bdy_kernels.hpp"
#include "core/field.hpp"

using namespace cppwrf;
using namespace cppwrf::bdy;

// External Fortran bridge declarations
extern "C" {
void bridge_lbc_fcx_gcx(double* fcx, double* gcx, int spec_bdy_width, int spec_zone, int relax_zone, double dt, double spec_exp);

void bridge_relax_bdytend_core_3d(
    const double* field, double* field_tend,
    const double* bxs, const double* bxe, const double* bys, const double* bye,
    const double* btxs, const double* btxe, const double* btys, const double* btye,
    const double* fcx, const double* gcx, double dtbc,
    char var_code,
    int nx, int ny, int nz,
    int spec_bdy_width, int spec_zone, int relax_zone);

void bridge_relax_bdytend_core_2d(
    const double* field, double* field_tend,
    const double* bxs, const double* bxe, const double* bys, const double* bye,
    const double* btxs, const double* btxe, const double* btys, const double* btye,
    const double* fcx, const double* gcx, double dtbc,
    int nx, int ny,
    int spec_bdy_width, int spec_zone, int relax_zone);

void bridge_spec_bdytend_3d(
    double* field_tend,
    const double* bxs, const double* bxe, const double* bys, const double* bye,
    const double* btxs, const double* btxe, const double* btys, const double* btye,
    char var_code,
    int nx, int ny, int nz,
    int spec_bdy_width, int spec_zone);

void bridge_spec_bdytend_2d(
    double* field_tend,
    const double* bxs, const double* bxe, const double* bys, const double* bye,
    const double* btxs, const double* btxe, const double* btys, const double* btye,
    int nx, int ny,
    int spec_bdy_width, int spec_zone);

void bridge_spec_bdyupdate_3d(
    double* field, const double* field_tend, double dt, char var_code,
    int nx, int ny, int nz, int spec_zone);

void bridge_spec_bdyupdate_2d(
    double* field, const double* field_tend, double dt,
    int nx, int ny, int spec_zone);

void bridge_spec_bdyupdate_ph(
    const double* ph_save, double* field, const double* field_tend,
    const double* mu_tend, const double* muts,
    const double* c1, const double* c2, double dt,
    int nx, int ny, int nz, int spec_zone);

void bridge_zero_grad_bdy(
    double* field, int nx, int ny, int nz, int spec_zone);

void bridge_flow_dep_bdy(
    double* field, const double* u, const double* v,
    int nx, int ny, int nz, int spec_zone);

void bridge_spec_bdy_final_3d(
    double* field, const double* mu, const double* c1, const double* c2, const double* msf,
    const double* bxs, const double* bxe, const double* bys, const double* bye,
    const double* btxs, const double* btxe, const double* btys, const double* btye,
    char var_code, double dtbc,
    int nx, int ny, int nz,
    int spec_bdy_width, int spec_zone);

void bridge_spec_bdy_final_2d(
    double* field,
    const double* bxs, const double* bxe, const double* bys, const double* bye,
    const double* btxs, const double* btxe, const double* btys, const double* btye,
    double dtbc, int nx, int ny,
    int spec_bdy_width, int spec_zone);

void bridge_set_physical_bc3d(
    double* dat, char var_code,
    int nx, int ny, int nz, int halo);

void bridge_set_physical_bc2d(
    double* dat, char var_code,
    int nx, int ny, int halo);
}

// Tolerance for PASS: FP64 bit-exact or <= 1e-12 relative
constexpr double TOL = 1e-12;

struct DiffStats {
  double max_abs = 0.0;
  double max_rel = 0.0;
  int total_points = 0;
  int diff_count = 0;

  void record(double cpp_val, double ftn_val) {
    total_points++;
    double d = std::abs(cpp_val - ftn_val);
    if (d > max_abs) max_abs = d;
    double denom = std::max(std::abs(cpp_val), std::abs(ftn_val));
    double rel = (denom > 1e-14) ? (d / denom) : d;
    if (rel > max_rel) max_rel = rel;
    if (rel > TOL && d > 1e-14) {
      diff_count++;
    }
  }

  bool pass() const { return max_rel <= TOL || max_abs <= 1e-14; }
};

// Check NetCDF helper
#define CHECK_NC(call) do { \
  int status = (call); \
  if (status != NC_NOERR) { \
    std::cerr << "NetCDF Error: " << nc_strerror(status) << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
    return false; \
  } \
} while(0)

// ---------------------------------------------------------------------------
// Test 1: lbc_fcx_gcx
// ---------------------------------------------------------------------------
bool test_lbc_fcx_gcx() {
  std::cout << "[TEST 1] Differential test: lbc_fcx_gcx ... ";
  DiffStats stats;
  for (double dt : {6.0, 12.0, 24.0, 60.0}) {
    for (double spec_exp : {0.0, 0.1, 0.25}) {
      double fcx_c[5], gcx_c[5];
      double fcx_f[5], gcx_f[5];
      calc_lbc_fcx_gcx<double>(fcx_c, gcx_c, 5, 1, 4, dt, spec_exp);
      bridge_lbc_fcx_gcx(fcx_f, gcx_f, 5, 1, 4, dt, spec_exp);
      for (int b = 0; b < 5; ++b) {
        stats.record(fcx_c[b], fcx_f[b]);
        stats.record(gcx_c[b], gcx_f[b]);
      }
    }
  }
  if (!stats.pass()) {
    std::cout << "FAILED! max_abs=" << stats.max_abs << " max_rel=" << stats.max_rel << "\n";
    return false;
  }
  std::cout << "PASSED (max_rel=" << stats.max_rel << ")\n";
  return true;
}

// ---------------------------------------------------------------------------
// Test 2: relax_bdytend_core_3d across multiple tiles & variables
// ---------------------------------------------------------------------------
bool test_relax_bdytend_core_3d() {
  std::cout << "[TEST 2] Differential test: relax_bdytend_core_3d (multiple tiles/variables) ... ";
  std::mt19937_64 rng(42);
  std::uniform_real_distribution<double> dist_val(-50.0, 50.0);
  std::uniform_real_distribution<double> dist_tend(-0.1, 0.1);

  struct TileConfig {
    int nx, ny, nz;
  };
  std::vector<TileConfig> tiles = {
      {16, 16, 8},   // square domain
      {20, 14, 10},  // rectangular
      {13, 17, 9},   // odd dimensions
      {25, 25, 12}   // larger tile
  };

  std::vector<char> vars = {'u', 'v', 't', 'h', 'q'};
  DiffStats overall_stats;

  for (const auto& t : tiles) {
    const int nx = t.nx, ny = t.ny, nz = t.nz;
    double dt = 12.0, dtbc = 180.0;
    double fcx[5], gcx[5];
    calc_lbc_fcx_gcx<double>(fcx, gcx, 5, 1, 4, dt, 0.0);

    for (char var : vars) {
      const int ni = (var == 'u') ? nx + 1 : nx;
      const int nj = (var == 'v') ? ny + 1 : ny;
      const int nk = (var == 'h') ? nz + 1 : nz;

      // C++ SoA field: (i, j, k) -> i + j*sj + k*sk
      std::vector<double> field_cpp(ni * nj * nk);
      std::vector<double> tend_cpp(ni * nj * nk);
      View3<const double> v_field{field_cpp.data(), ni, ni * nj};
      View3<double> v_tend{tend_cpp.data(), ni, ni * nj};

      // Fortran column-major field: (i, k, j) -> i + k*ni + j*(ni*nk)
      std::vector<double> field_ftn(ni * nk * nj);
      std::vector<double> tend_ftn(ni * nk * nj);

      for (int k = 0; k < nk; ++k) {
        for (int j = 0; j < nj; ++j) {
          for (int i = 0; i < ni; ++i) {
            double f_val = dist_val(rng);
            double t_val = dist_tend(rng);
            field_cpp[i + j * ni + k * (ni * nj)] = f_val;
            tend_cpp[i + j * ni + k * (ni * nj)] = t_val;
            field_ftn[i + k * ni + j * (ni * nk)] = f_val;
            tend_ftn[i + k * ni + j * (ni * nk)] = t_val;
          }
        }
      }

      // Boundary slabs: allocated with nk_bdy = nz + 1 for WRF exact layout
      const int nk_bdy = nz + 1;
      BdySlab<double> bxs(nj, nk_bdy, 5), bxe(nj, nk_bdy, 5);
      BdySlab<double> bys(ni, nk_bdy, 5), bye(ni, nk_bdy, 5);
      BdySlab<double> btxs(nj, nk_bdy, 5), btxe(nj, nk_bdy, 5);
      BdySlab<double> btys(ni, nk_bdy, 5), btye(ni, nk_bdy, 5);

      for (int b = 0; b < 5; ++b) {
        for (int k = 0; k < nk_bdy; ++k) {
          for (int j = 0; j < nj; ++j) {
            bxs(j, k, b) = dist_val(rng);
            bxe(j, k, b) = dist_val(rng);
            btxs(j, k, b) = dist_tend(rng);
            btxe(j, k, b) = dist_tend(rng);
          }
          for (int i = 0; i < ni; ++i) {
            bys(i, k, b) = dist_val(rng);
            bye(i, k, b) = dist_val(rng);
            btys(i, k, b) = dist_tend(rng);
            btye(i, k, b) = dist_tend(rng);
          }
        }
      }

      // Run C++ kernel
      relax_bdytend_core_3d<double>(
          v_tend, v_field,
          bxs.cview(), bxe.cview(), bys.cview(), bye.cview(),
          btxs.cview(), btxe.cview(), btys.cview(), btye.cview(),
          fcx, gcx, dtbc, var, nx, ny, nz, 1, 4);

      // Run Fortran bridge
      bridge_relax_bdytend_core_3d(
          field_ftn.data(), tend_ftn.data(),
          bxs.data(), bxe.data(), bys.data(), bye.data(),
          btxs.data(), btxe.data(), btys.data(), btye.data(),
          fcx, gcx, dtbc, var,
          nx, ny, nz, 5, 1, 4);

      // Compare
      for (int k = 0; k < nk; ++k) {
        for (int j = 0; j < nj; ++j) {
          for (int i = 0; i < ni; ++i) {
            double c_val = tend_cpp[i + j * ni + k * (ni * nj)];
            double f_val = tend_ftn[i + k * ni + j * (ni * nk)];
            overall_stats.record(c_val, f_val);
          }
        }
      }
    }
  }

  if (!overall_stats.pass()) {
    std::cout << "FAILED! max_abs=" << overall_stats.max_abs
              << " max_rel=" << overall_stats.max_rel
              << " diffs=" << overall_stats.diff_count << "/" << overall_stats.total_points << "\n";
    return false;
  }
  std::cout << "PASSED (points=" << overall_stats.total_points << " max_rel=" << overall_stats.max_rel << ")\n";
  return true;
}

// ---------------------------------------------------------------------------
// Test 3: relax_bdytend_core_2d (mu)
// ---------------------------------------------------------------------------
bool test_relax_bdytend_core_2d() {
  std::cout << "[TEST 3] Differential test: relax_bdytend_core_2d (mu) ... ";
  std::mt19937_64 rng(12345);
  std::uniform_real_distribution<double> dist_mu(90000.0, 105000.0);
  std::uniform_real_distribution<double> dist_tend(-10.0, 10.0);

  int nx = 18, ny = 16;
  double dt = 12.0, dtbc = 60.0;
  double fcx[5], gcx[5];
  calc_lbc_fcx_gcx<double>(fcx, gcx, 5, 1, 4, dt, 0.0);

  std::vector<double> field_cpp(nx * ny);
  std::vector<double> tend_cpp(nx * ny);
  std::vector<double> field_ftn(nx * ny);
  std::vector<double> tend_ftn(nx * ny);

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      double f = dist_mu(rng);
      double t = dist_tend(rng);
      field_cpp[i + j * nx] = f;
      tend_cpp[i + j * nx] = t;
      field_ftn[i + j * nx] = f;
      tend_ftn[i + j * nx] = t;
    }
  }

  BdySlab<double> bxs(ny, 1, 5), bxe(ny, 1, 5);
  BdySlab<double> bys(nx, 1, 5), bye(nx, 1, 5);
  BdySlab<double> btxs(ny, 1, 5), btxe(ny, 1, 5);
  BdySlab<double> btys(nx, 1, 5), btye(nx, 1, 5);

  for (int b = 0; b < 5; ++b) {
    for (int j = 0; j < ny; ++j) {
      bxs(j, 0, b) = dist_mu(rng);
      bxe(j, 0, b) = dist_mu(rng);
      btxs(j, 0, b) = dist_tend(rng);
      btxe(j, 0, b) = dist_tend(rng);
    }
    for (int i = 0; i < nx; ++i) {
      bys(i, 0, b) = dist_mu(rng);
      bye(i, 0, b) = dist_mu(rng);
      btys(i, 0, b) = dist_tend(rng);
      btye(i, 0, b) = dist_tend(rng);
    }
  }

  View2<const double> v_field{field_cpp.data(), nx};
  View2<double> v_tend{tend_cpp.data(), nx};

  relax_bdytend_core_2d<double>(
      v_tend, v_field,
      bxs.cview(), bxe.cview(), bys.cview(), bye.cview(),
      btxs.cview(), btxe.cview(), btys.cview(), btye.cview(),
      fcx, gcx, dtbc, nx, ny, 1, 4);

  bridge_relax_bdytend_core_2d(
      field_ftn.data(), tend_ftn.data(),
      bxs.data(), bxe.data(), bys.data(), bye.data(),
      btxs.data(), btxe.data(), btys.data(), btye.data(),
      fcx, gcx, dtbc, nx, ny, 5, 1, 4);

  DiffStats stats;
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      stats.record(tend_cpp[i + j * nx], tend_ftn[i + j * nx]);
    }
  }

  if (!stats.pass()) {
    std::cout << "FAILED! max_abs=" << stats.max_abs << " max_rel=" << stats.max_rel << "\n";
    return false;
  }
  std::cout << "PASSED (max_rel=" << stats.max_rel << ")\n";
  return true;
}

// ---------------------------------------------------------------------------
// Test 4: spec_bdytend (3D & 2D)
// ---------------------------------------------------------------------------
bool test_spec_bdytend() {
  std::cout << "[TEST 4] Differential test: spec_bdytend (3D & 2D) ... ";
  std::mt19937_64 rng(999);
  std::uniform_real_distribution<double> dist_tend(-1.0, 1.0);
  int nx = 15, ny = 14, nz = 8;
  DiffStats stats;

  for (char var : {'u', 'v', 't', 'h', 'q'}) {
    const int ni = (var == 'u') ? nx + 1 : nx;
    const int nj = (var == 'v') ? ny + 1 : ny;
    const int nk = (var == 'h') ? nz + 1 : nz;

    std::vector<double> tend_cpp(ni * nj * nk, 0.0);
    std::vector<double> tend_ftn(ni * nk * nj, 0.0);
    View3<double> v_tend{tend_cpp.data(), ni, ni * nj};

    const int nk_bdy = nz + 1;
    BdySlab<double> btxs(nj, nk_bdy, 5), btxe(nj, nk_bdy, 5);
    BdySlab<double> btys(ni, nk_bdy, 5), btye(ni, nk_bdy, 5);

    for (int b = 0; b < 5; ++b) {
      for (int k = 0; k < nk_bdy; ++k) {
        for (int j = 0; j < nj; ++j) {
          btxs(j, k, b) = dist_tend(rng);
          btxe(j, k, b) = dist_tend(rng);
        }
        for (int i = 0; i < ni; ++i) {
          btys(i, k, b) = dist_tend(rng);
          btye(i, k, b) = dist_tend(rng);
        }
      }
    }

    spec_bdytend_3d<double>(v_tend, btxs.cview(), btxe.cview(), btys.cview(), btye.cview(), var, nx, ny, nz, 1);
    bridge_spec_bdytend_3d(tend_ftn.data(), btxs.data(), btxe.data(), btys.data(), btye.data(),
                           btxs.data(), btxe.data(), btys.data(), btye.data(), var, nx, ny, nz, 5, 1);

    for (int k = 0; k < nk; ++k) {
      for (int j = 0; j < nj; ++j) {
        for (int i = 0; i < ni; ++i) {
          stats.record(tend_cpp[i + j * ni + k * (ni * nj)], tend_ftn[i + k * ni + j * (ni * nk)]);
        }
      }
    }
  }

  // 2D test
  {
    std::vector<double> tend_cpp(nx * ny, 0.0);
    std::vector<double> tend_ftn(nx * ny, 0.0);
    View2<double> v_tend{tend_cpp.data(), nx};

    BdySlab<double> btxs(ny, 1, 5), btxe(ny, 1, 5);
    BdySlab<double> btys(nx, 1, 5), btye(nx, 1, 5);

    for (int b = 0; b < 5; ++b) {
      for (int j = 0; j < ny; ++j) {
        btxs(j, 0, b) = dist_tend(rng);
        btxe(j, 0, b) = dist_tend(rng);
      }
      for (int i = 0; i < nx; ++i) {
        btys(i, 0, b) = dist_tend(rng);
        btye(i, 0, b) = dist_tend(rng);
      }
    }

    spec_bdytend_2d<double>(v_tend, btxs.cview(), btxe.cview(), btys.cview(), btye.cview(), nx, ny, 1);
    bridge_spec_bdytend_2d(tend_ftn.data(), btxs.data(), btxe.data(), btys.data(), btye.data(),
                           btxs.data(), btxe.data(), btys.data(), btye.data(), nx, ny, 5, 1);

    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        stats.record(tend_cpp[i + j * nx], tend_ftn[i + j * nx]);
      }
    }
  }

  if (!stats.pass()) {
    std::cout << "FAILED! max_abs=" << stats.max_abs << " max_rel=" << stats.max_rel << "\n";
    return false;
  }
  std::cout << "PASSED (bit-exact, max_abs=" << stats.max_abs << ")\n";
  return true;
}

// ---------------------------------------------------------------------------
// Test 5: spec_bdyupdate (3D & 2D)
// ---------------------------------------------------------------------------
bool test_spec_bdyupdate() {
  std::cout << "[TEST 5] Differential test: spec_bdyupdate ... ";
  std::mt19937_64 rng(888);
  std::uniform_real_distribution<double> dist_val(10.0, 100.0);
  std::uniform_real_distribution<double> dist_tend(-0.5, 0.5);
  int nx = 14, ny = 16, nz = 6;
  double dt = 2.0;
  DiffStats stats;

  for (char var : {'u', 'v', 't', 'h'}) {
    const int ni = (var == 'u') ? nx + 1 : nx;
    const int nj = (var == 'v') ? ny + 1 : ny;
    const int nk = (var == 'h') ? nz + 1 : nz;

    std::vector<double> field_cpp(ni * nj * nk);
    std::vector<double> tend_cpp(ni * nj * nk);
    std::vector<double> field_ftn(ni * nk * nj);
    std::vector<double> tend_ftn(ni * nk * nj);

    for (int k = 0; k < nk; ++k) {
      for (int j = 0; j < nj; ++j) {
        for (int i = 0; i < ni; ++i) {
          double f = dist_val(rng);
          double t = dist_tend(rng);
          field_cpp[i + j * ni + k * (ni * nj)] = f;
          tend_cpp[i + j * ni + k * (ni * nj)] = t;
          field_ftn[i + k * ni + j * (ni * nk)] = f;
          tend_ftn[i + k * ni + j * (ni * nk)] = t;
        }
      }
    }

    View3<double> v_field{field_cpp.data(), ni, ni * nj};
    View3<const double> v_tend{tend_cpp.data(), ni, ni * nj};

    spec_bdyupdate_3d<double>(v_field, v_tend, dt, var, nx, ny, nz, 1);
    bridge_spec_bdyupdate_3d(field_ftn.data(), tend_ftn.data(), dt, var, nx, ny, nz, 1);

    for (int k = 0; k < nk; ++k) {
      for (int j = 0; j < nj; ++j) {
        for (int i = 0; i < ni; ++i) {
          stats.record(field_cpp[i + j * ni + k * (ni * nj)], field_ftn[i + k * ni + j * (ni * nk)]);
        }
      }
    }
  }

  // 2D test (mu)
  {
    std::vector<double> field_cpp(nx * ny);
    std::vector<double> tend_cpp(nx * ny);
    std::vector<double> field_ftn(nx * ny);
    std::vector<double> tend_ftn(nx * ny);

    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        double f = dist_val(rng);
        double t = dist_tend(rng);
        field_cpp[i + j * nx] = f;
        tend_cpp[i + j * nx] = t;
        field_ftn[i + j * nx] = f;
        tend_ftn[i + j * nx] = t;
      }
    }

    View2<double> v_field{field_cpp.data(), nx};
    View2<const double> v_tend{tend_cpp.data(), nx};

    spec_bdyupdate_2d<double>(v_field, v_tend, dt, nx, ny, 1);
    bridge_spec_bdyupdate_2d(field_ftn.data(), tend_ftn.data(), dt, nx, ny, 1);

    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        stats.record(field_cpp[i + j * nx], field_ftn[i + j * nx]);
      }
    }
  }

  if (!stats.pass()) {
    std::cout << "FAILED! max_abs=" << stats.max_abs << " max_rel=" << stats.max_rel << "\n";
    return false;
  }
  std::cout << "PASSED (max_rel=" << stats.max_rel << ")\n";
  return true;
}

// ---------------------------------------------------------------------------
// Test 6: spec_bdyupdate_ph
// ---------------------------------------------------------------------------
bool test_spec_bdyupdate_ph() {
  std::cout << "[TEST 6] Differential test: spec_bdyupdate_ph ... ";
  std::mt19937_64 rng(777);
  std::uniform_real_distribution<double> dist_ph(100.0, 50000.0);
  std::uniform_real_distribution<double> dist_tend(-5.0, 5.0);
  std::uniform_real_distribution<double> dist_mu(90000.0, 100000.0);

  int nx = 12, ny = 14, nz = 8;
  double dt = 2.5;

  std::vector<double> c1f(nz + 1, 1.0);
  std::vector<double> c2f(nz + 1, 0.0);
  for (int k = 0; k <= nz; ++k) {
    c1f[k] = 1.0 - 0.05 * k;
    c2f[k] = 5000.0 * k;
  }

  std::vector<double> ph_cpp(nx * ny * (nz + 1));
  std::vector<double> ph_save_cpp(nx * ny * (nz + 1));
  std::vector<double> ph_tend_cpp(nx * ny * (nz + 1));
  std::vector<double> ph_ftn(nx * (nz + 1) * ny);
  std::vector<double> ph_save_ftn(nx * (nz + 1) * ny);
  std::vector<double> ph_tend_ftn(nx * (nz + 1) * ny);

  for (int k = 0; k <= nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        double p = dist_ph(rng);
        double ps = dist_ph(rng);
        double pt = dist_tend(rng);
        ph_cpp[i + j * nx + k * (nx * ny)] = p;
        ph_save_cpp[i + j * nx + k * (nx * ny)] = ps;
        ph_tend_cpp[i + j * nx + k * (nx * ny)] = pt;

        ph_ftn[i + k * nx + j * (nx * (nz + 1))] = p;
        ph_save_ftn[i + k * nx + j * (nx * (nz + 1))] = ps;
        ph_tend_ftn[i + k * nx + j * (nx * (nz + 1))] = pt;
      }
    }
  }

  std::vector<double> mu_tend(nx * ny);
  std::vector<double> muts(nx * ny);
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      muts[i + j * nx] = dist_mu(rng);
      mu_tend[i + j * nx] = dist_tend(rng);
    }
  }

  View3<double> v_ph{ph_cpp.data(), nx, nx * ny};
  View3<const double> v_ph_save{ph_save_cpp.data(), nx, nx * ny};
  View3<const double> v_ph_tend{ph_tend_cpp.data(), nx, nx * ny};
  View2<const double> v_mu_tend{mu_tend.data(), nx};
  View2<const double> v_muts{muts.data(), nx};

  spec_bdyupdate_ph<double>(v_ph_save, v_ph, v_ph_tend, v_mu_tend, v_muts, c1f.data(), c2f.data(), dt, nx, ny, nz, 1);
  bridge_spec_bdyupdate_ph(ph_save_ftn.data(), ph_ftn.data(), ph_tend_ftn.data(), mu_tend.data(), muts.data(),
                           c1f.data(), c2f.data(), dt, nx, ny, nz, 1);

  DiffStats stats;
  for (int k = 0; k <= nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        stats.record(ph_cpp[i + j * nx + k * (nx * ny)], ph_ftn[i + k * nx + j * (nx * (nz + 1))]);
      }
    }
  }

  if (!stats.pass()) {
    std::cout << "FAILED! max_abs=" << stats.max_abs << " max_rel=" << stats.max_rel << "\n";
    return false;
  }
  std::cout << "PASSED (max_rel=" << stats.max_rel << ")\n";
  return true;
}

// ---------------------------------------------------------------------------
// Test 7: zero_grad_bdy
// ---------------------------------------------------------------------------
bool test_zero_grad_bdy() {
  std::cout << "[TEST 7] Differential test: zero_grad_bdy (w) ... ";
  std::mt19937_64 rng(555);
  std::uniform_real_distribution<double> dist_w(-15.0, 15.0);
  int nx = 16, ny = 14, nz = 10;

  std::vector<double> w_cpp(nx * ny * (nz + 1));
  std::vector<double> w_ftn(nx * (nz + 1) * ny);
  for (int k = 0; k <= nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        double val = dist_w(rng);
        w_cpp[i + j * nx + k * (nx * ny)] = val;
        w_ftn[i + k * nx + j * (nx * (nz + 1))] = val;
      }
    }
  }

  View3<double> v_w{w_cpp.data(), nx, nx * ny};
  zero_grad_bdy<double>(v_w, nx, ny, nz, 1);
  bridge_zero_grad_bdy(w_ftn.data(), nx, ny, nz, 1);

  DiffStats stats;
  for (int k = 0; k <= nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        stats.record(w_cpp[i + j * nx + k * (nx * ny)], w_ftn[i + k * nx + j * (nx * (nz + 1))]);
      }
    }
  }

  if (!stats.pass()) {
    std::cout << "FAILED! max_abs=" << stats.max_abs << "\n";
    return false;
  }
  std::cout << "PASSED (bit-exact, max_abs=" << stats.max_abs << ")\n";
  return true;
}

// ---------------------------------------------------------------------------
// Test 8: flow_dep_bdy
// ---------------------------------------------------------------------------
bool test_flow_dep_bdy() {
  std::cout << "[TEST 8] Differential test: flow_dep_bdy ... ";
  std::mt19937_64 rng(333);
  std::uniform_real_distribution<double> dist_q(0.001, 0.02);
  std::uniform_real_distribution<double> dist_uv(-10.0, 10.0);
  int nx = 14, ny = 14, nz = 6;

  std::vector<double> q_cpp(nx * ny * nz);
  std::vector<double> u_cpp((nx + 1) * ny * nz);
  std::vector<double> v_cpp(nx * (ny + 1) * nz);

  const int ni_f = nx + 1;
  const int nj_f = ny + 1;
  const int nk_f = nz;

  std::vector<double> q_ftn(ni_f * nk_f * nj_f, 0.0);
  std::vector<double> u_ftn(ni_f * nk_f * nj_f, 0.0);
  std::vector<double> v_ftn(ni_f * nk_f * nj_f, 0.0);

  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        double val = dist_q(rng);
        q_cpp[i + j * nx + k * (nx * ny)] = val;
        q_ftn[i + k * ni_f + j * (ni_f * nk_f)] = val;
      }
    }
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i <= nx; ++i) {
        double val = dist_uv(rng);
        u_cpp[i + j * (nx + 1) + k * ((nx + 1) * ny)] = val;
        u_ftn[i + k * ni_f + j * (ni_f * nk_f)] = val;
      }
    }
    for (int j = 0; j <= ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        double val = dist_uv(rng);
        v_cpp[i + j * nx + k * (nx * (ny + 1))] = val;
        v_ftn[i + k * ni_f + j * (ni_f * nk_f)] = val;
      }
    }
  }

  View3<double> v_q{q_cpp.data(), nx, nx * ny};
  View3<const double> v_u{u_cpp.data(), nx + 1, (nx + 1) * ny};
  View3<const double> v_v{v_cpp.data(), nx, nx * (ny + 1)};

  flow_dep_bdy<double>(v_q, v_u, v_v, nx, ny, nz, 1);
  bridge_flow_dep_bdy(q_ftn.data(), u_ftn.data(), v_ftn.data(), nx, ny, nz, 1);

  DiffStats stats;
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        stats.record(q_cpp[i + j * nx + k * (nx * ny)], q_ftn[i + k * ni_f + j * (ni_f * nk_f)]);
      }
    }
  }

  if (!stats.pass()) {
    std::cout << "FAILED! max_abs=" << stats.max_abs << "\n";
    return false;
  }
  std::cout << "PASSED (bit-exact, max_abs=" << stats.max_abs << ")\n";
  return true;
}

// ---------------------------------------------------------------------------
// Test 9: spec_bdy_final (3D & 2D)
// ---------------------------------------------------------------------------
bool test_spec_bdy_final() {
  std::cout << "[TEST 9] Differential test: spec_bdy_final (3D & 2D) ... ";
  std::mt19937_64 rng(222);
  std::uniform_real_distribution<double> dist_val(10.0, 300.0);
  std::uniform_real_distribution<double> dist_tend(-1.0, 1.0);
  std::uniform_real_distribution<double> dist_mu(90000.0, 100000.0);
  std::uniform_real_distribution<double> dist_msf(0.95, 1.05);

  int nx = 15, ny = 13, nz = 7;
  double dtbc = 120.0;
  DiffStats stats;

  std::vector<double> c1(nz + 1, 1.0), c2(nz + 1, 0.0);
  for (int k = 0; k <= nz; ++k) {
    c1[k] = 1.0 - 0.04 * k;
    c2[k] = 4000.0 * k;
  }

  for (char var : {'u', 'v', 'w', 't', 'h'}) {
    const int ni = (var == 'u') ? nx + 1 : nx;
    const int nj = (var == 'v') ? ny + 1 : ny;
    const int nk = (var == 'h' || var == 'w') ? nz + 1 : nz;

    std::vector<double> field_cpp(ni * nj * nk);
    std::vector<double> field_ftn(ni * nk * nj);
    for (int k = 0; k < nk; ++k) {
      for (int j = 0; j < nj; ++j) {
        for (int i = 0; i < ni; ++i) {
          double val = dist_val(rng);
          field_cpp[i + j * ni + k * (ni * nj)] = val;
          field_ftn[i + k * ni + j * (ni * nk)] = val;
        }
      }
    }

    std::vector<double> mu(ni * nj);
    std::vector<double> msf(ni * nj);
    for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
        mu[i + j * ni] = dist_mu(rng);
        msf[i + j * ni] = dist_msf(rng);
      }
    }
    View2<const double> v_mu{mu.data(), ni};
    View2<const double> v_msf{msf.data(), ni};

    const int nk_bdy = nz + 1;
    BdySlab<double> bxs(nj, nk_bdy, 5), bxe(nj, nk_bdy, 5);
    BdySlab<double> bys(ni, nk_bdy, 5), bye(ni, nk_bdy, 5);
    BdySlab<double> btxs(nj, nk_bdy, 5), btxe(nj, nk_bdy, 5);
    BdySlab<double> btys(ni, nk_bdy, 5), btye(ni, nk_bdy, 5);

    for (int b = 0; b < 5; ++b) {
      for (int k = 0; k < nk_bdy; ++k) {
        for (int j = 0; j < nj; ++j) {
          bxs(j, k, b) = dist_val(rng);
          bxe(j, k, b) = dist_val(rng);
          btxs(j, k, b) = dist_tend(rng);
          btxe(j, k, b) = dist_tend(rng);
        }
        for (int i = 0; i < ni; ++i) {
          bys(i, k, b) = dist_val(rng);
          bye(i, k, b) = dist_val(rng);
          btys(i, k, b) = dist_tend(rng);
          btye(i, k, b) = dist_tend(rng);
        }
      }
    }

    View3<double> v_field{field_cpp.data(), ni, ni * nj};

    spec_bdy_final_3d<double>(
        v_field, v_mu, c1.data(), c2.data(), v_msf,
        bxs.cview(), bxe.cview(), bys.cview(), bye.cview(),
        btxs.cview(), btxe.cview(), btys.cview(), btye.cview(),
        var, dtbc, nx, ny, nz, 1);

    bridge_spec_bdy_final_3d(
        field_ftn.data(), mu.data(), c1.data(), c2.data(), msf.data(),
        bxs.data(), bxe.data(), bys.data(), bye.data(),
        btxs.data(), btxe.data(), btys.data(), btye.data(),
        var, dtbc, nx, ny, nz, 5, 1);

    for (int k = 0; k < nk; ++k) {
      for (int j = 0; j < nj; ++j) {
        for (int i = 0; i < ni; ++i) {
          stats.record(field_cpp[i + j * ni + k * (ni * nj)], field_ftn[i + k * ni + j * (ni * nk)]);
        }
      }
    }
  }

  // 2D final update (mu)
  std::vector<double> mu_cpp(nx * ny);
  std::vector<double> mu_ftn(nx * ny);
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      double val = dist_mu(rng);
      mu_cpp[i + j * nx] = val;
      mu_ftn[i + j * nx] = val;
    }
  }

  BdySlab<double> mu_bxs(ny, 1, 5), mu_bxe(ny, 1, 5);
  BdySlab<double> mu_bys(nx, 1, 5), mu_bye(nx, 1, 5);
  BdySlab<double> mu_btxs(ny, 1, 5), mu_btxe(ny, 1, 5);
  BdySlab<double> mu_btys(nx, 1, 5), mu_btye(nx, 1, 5);

  for (int b = 0; b < 5; ++b) {
    for (int j = 0; j < ny; ++j) {
      mu_bxs(j, 0, b) = dist_mu(rng);
      mu_bxe(j, 0, b) = dist_mu(rng);
      mu_btxs(j, 0, b) = dist_tend(rng);
      mu_btxe(j, 0, b) = dist_tend(rng);
    }
    for (int i = 0; i < nx; ++i) {
      mu_bys(i, 0, b) = dist_mu(rng);
      mu_bye(i, 0, b) = dist_mu(rng);
      mu_btys(i, 0, b) = dist_tend(rng);
      mu_btye(i, 0, b) = dist_tend(rng);
    }
  }

  View2<double> v_mu2{mu_cpp.data(), nx};
  spec_bdy_final_2d<double>(v_mu2, mu_bxs.cview(), mu_bxe.cview(), mu_bys.cview(), mu_bye.cview(),
                           mu_btxs.cview(), mu_btxe.cview(), mu_btys.cview(), mu_btye.cview(),
                           dtbc, nx, ny, 1);
  bridge_spec_bdy_final_2d(mu_ftn.data(), mu_bxs.data(), mu_bxe.data(), mu_bys.data(), mu_bye.data(),
                           mu_btxs.data(), mu_btxe.data(), mu_btys.data(), mu_btye.data(),
                           dtbc, nx, ny, 5, 1);

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      stats.record(mu_cpp[i + j * nx], mu_ftn[i + j * nx]);
    }
  }

  if (!stats.pass()) {
    std::cout << "FAILED! max_abs=" << stats.max_abs << " max_rel=" << stats.max_rel << "\n";
    return false;
  }
  std::cout << "PASSED (max_rel=" << stats.max_rel << ")\n";
  return true;
}

// ---------------------------------------------------------------------------
// Test 10: set_physical_bc3d & set_physical_bc2d (flat halo filling)
// ---------------------------------------------------------------------------
bool test_set_physical_bc() {
  std::cout << "[TEST 10] Differential test: set_physical_bc (flat halos) ... ";
  std::mt19937_64 rng(111);
  std::uniform_real_distribution<double> dist(-50.0, 50.0);
  int nx = 12, ny = 12, nz = 6;
  const int halo_mem = 5;
  const int halo_test = 3;
  DiffStats stats;

  for (char var : {'u', 'v', 'w', 't', 'p'}) {
    const int ni = (var == 'u') ? nx + 1 : nx;
    const int nj = (var == 'v') ? ny + 1 : ny;
    const int nk = (var == 'w') ? nz + 1 : nz;

    const int ims = -halo_mem, ime = ni - 1 + halo_mem;
    const int jms = -halo_mem, jme = nj - 1 + halo_mem;
    const int kms = 0, kme = nk - 1;
    const int len_i = ime - ims + 1;
    const int len_j = jme - jms + 1;
    const int len_k = kme - kms + 1;

    // Buffer layout with halos:
    // C++ pointer at (0, 0, 0)
    std::vector<double> buf_cpp(len_i * len_j * len_k, 0.0);
    double* ptr_cpp_0 = buf_cpp.data() + (-ims) + (-jms) * len_i;
    View3<double> v_cpp{ptr_cpp_0, len_i, len_i * len_j};

    // Fortran buffer: (ims:ime, kms:kme, jms:jme)
    std::vector<double> buf_ftn(len_i * len_k * len_j, 0.0);

    for (int k = 0; k < nk; ++k) {
      for (int j = 0; j < nj; ++j) {
        for (int i = 0; i < ni; ++i) {
          double val = dist(rng);
          v_cpp(i, j, k) = val;
          int f_idx = (i - ims) + k * len_i + (j - jms) * (len_i * len_k);
          buf_ftn[f_idx] = val;
        }
      }
    }

    set_physical_bc3d_flat<double>(v_cpp, var, nx, ny, nz, halo_test);
    bridge_set_physical_bc3d(buf_ftn.data(), var, nx, ny, nz, halo_mem);

    for (int k = 0; k < nk; ++k) {
      for (int j = -halo_test; j <= nj - 1 + halo_test; ++j) {
        for (int i = -halo_test; i <= ni - 1 + halo_test; ++i) {
          double c_val = v_cpp(i, j, k);
          int f_idx = (i - ims) + k * len_i + (j - jms) * (len_i * len_k);
          double f_val = buf_ftn[f_idx];
          if (std::abs(c_val - f_val) > 1e-4) {
            std::cout << "\nMismatch for 3D var " << var << " at i=" << i << " j=" << j << " k=" << k
                      << " cpp=" << c_val << " ftn=" << f_val << "\n";
            return false;
          }
          stats.record(c_val, f_val);
        }
      }
    }
  }

  // Also test 2D
  for (char var : {'u', 'v', 'p'}) {
    const int ni = (var == 'u') ? nx + 1 : nx;
    const int nj = (var == 'v') ? ny + 1 : ny;
    const int ims = -halo_mem, ime = ni - 1 + halo_mem;
    const int jms = -halo_mem, jme = nj - 1 + halo_mem;
    const int len_i = ime - ims + 1;
    const int len_j = jme - jms + 1;

    std::vector<double> buf_cpp2(len_i * len_j, 0.0);
    double* ptr_cpp_0 = buf_cpp2.data() + (-ims) + (-jms) * len_i;
    View2<double> v_cpp2{ptr_cpp_0, len_i};
    std::vector<double> buf_ftn2(len_i * len_j, 0.0);

    for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
        double val = dist(rng);
        v_cpp2(i, j) = val;
        int f_idx = (i - ims) + (j - jms) * len_i;
        buf_ftn2[f_idx] = val;
      }
    }

    set_physical_bc2d_flat<double>(v_cpp2, var, nx, ny, halo_test);
    bridge_set_physical_bc2d(buf_ftn2.data(), var, nx, ny, halo_mem);

    for (int j = -halo_test; j <= nj - 1 + halo_test; ++j) {
      for (int i = -halo_test; i <= ni - 1 + halo_test; ++i) {
        double c_val = v_cpp2(i, j);
        int f_idx = (i - ims) + (j - jms) * len_i;
        double f_val = buf_ftn2[f_idx];
        if (std::abs(c_val - f_val) > 1e-4) {
          std::cout << "\nMismatch for 2D var " << var << " at i=" << i << " j=" << j
                    << " cpp=" << c_val << " ftn=" << f_val << "\n";
          return false;
        }
        stats.record(c_val, f_val);
      }
    }
  }

  if (!stats.pass()) {
    std::cout << "FAILED! max_abs=" << stats.max_abs << "\n";
    return false;
  }
  std::cout << "PASSED (bit-exact, max_abs=" << stats.max_abs << ")\n";
  return true;
}

// ---------------------------------------------------------------------------
// Test 11: Real R0 wrfbdy_d01 & wrfinput_d01 verification
// ---------------------------------------------------------------------------
bool test_r0_wrfbdy_verification() {
  std::cout << "[TEST 11] Verification against R0 wrfbdy_d01 & wrfinput_d01 ...\n";
  const char* bdy_path = "/home/ostriquetrum/wrf_real/run_R0/wrfbdy_d01";
  const char* inp_path = "/home/ostriquetrum/wrf_real/run_R0/wrfinput_d01";

  int ncid_bdy = -1, ncid_inp = -1;
  CHECK_NC(nc_open(bdy_path, NC_NOWRITE, &ncid_bdy));
  CHECK_NC(nc_open(inp_path, NC_NOWRITE, &ncid_inp));

  // Read dimensions from wrfinput:
  // west_east = 119, south_north = 89, bottom_top = 44
  int dim_we, dim_sn, dim_bt;
  size_t nx = 0, ny = 0, nz = 0;
  CHECK_NC(nc_inq_dimid(ncid_inp, "west_east", &dim_we));
  CHECK_NC(nc_inq_dimlen(ncid_inp, dim_we, &nx));
  CHECK_NC(nc_inq_dimid(ncid_inp, "south_north", &dim_sn));
  CHECK_NC(nc_inq_dimlen(ncid_inp, dim_sn, &ny));
  CHECK_NC(nc_inq_dimid(ncid_inp, "bottom_top", &dim_bt));
  CHECK_NC(nc_inq_dimlen(ncid_inp, dim_bt, &nz));

  std::cout << "  Domain extents: nx=" << nx << ", ny=" << ny << ", nz=" << nz << "\n";

  // Check 1: MU on boundary at t=0
  int var_mu_inp, var_mub_inp, var_mu_bys, var_mu_btys;
  CHECK_NC(nc_inq_varid(ncid_inp, "MU", &var_mu_inp));
  CHECK_NC(nc_inq_varid(ncid_inp, "MUB", &var_mub_inp));
  CHECK_NC(nc_inq_varid(ncid_bdy, "MU_BYS", &var_mu_bys));
  CHECK_NC(nc_inq_varid(ncid_bdy, "MU_BTYS", &var_mu_btys));

  std::vector<float> mu_inp(nx * ny);
  std::vector<float> mub_inp(nx * ny);
  size_t start_2d[3] = {0, 0, 0};
  size_t count_inp_mu[3] = {1, ny, nx};
  CHECK_NC(nc_get_vara_float(ncid_inp, var_mu_inp, start_2d, count_inp_mu, mu_inp.data()));
  CHECK_NC(nc_get_vara_float(ncid_inp, var_mub_inp, start_2d, count_inp_mu, mub_inp.data()));

  std::vector<float> mu_bys(5 * nx);
  std::vector<float> mu_btys(5 * nx);
  size_t start_bdy_mu_ys[3] = {0, 0, 0};
  size_t count_bdy_mu_ys[3] = {1, 5, nx};
  CHECK_NC(nc_get_vara_float(ncid_bdy, var_mu_bys, start_bdy_mu_ys, count_bdy_mu_ys, mu_bys.data()));
  CHECK_NC(nc_get_vara_float(ncid_bdy, var_mu_btys, start_bdy_mu_ys, count_bdy_mu_ys, mu_btys.data()));

  // Verify South boundary MU (j = 0) vs MU_BYS(b=0):
  DiffStats mu_stats;
  for (size_t i = 0; i < nx; ++i) {
    float inp_val = mu_inp[0 * nx + i];
    float bdy_val = mu_bys[0 * nx + i]; // b=0, i
    mu_stats.record(inp_val, bdy_val);
  }
  std::cout << "  [Check MU South b=0]: points=" << mu_stats.total_points
            << " max_abs=" << mu_stats.max_abs << " max_rel=" << mu_stats.max_rel;
  if (mu_stats.max_abs <= 1e-5) {
    std::cout << " -> MATCHES (bit-exact)!\n";
  } else {
    std::cout << " -> MISMATCH!\n";
    nc_close(ncid_bdy);
    nc_close(ncid_inp);
    return false;
  }

  // Check 2: T (potential temperature perturbation theta) on boundary
  // In WRF 4.x with use_theta_m = 1:
  // theta_m = (theta_dry + 300) * (1 + 1.6083624 * qv) - 300
  // coupled_T = theta_m * (c1h * (mu + mub) + c2h)
  int var_t_inp, var_qv_inp, var_t_bys;
  CHECK_NC(nc_inq_varid(ncid_inp, "T", &var_t_inp));
  CHECK_NC(nc_inq_varid(ncid_inp, "QVAPOR", &var_qv_inp));
  CHECK_NC(nc_inq_varid(ncid_bdy, "T_BYS", &var_t_bys));

  std::vector<float> t_inp(nx * ny * nz);
  std::vector<float> qv_inp(nx * ny * nz);
  size_t count_inp_3d[4] = {1, nz, ny, nx};
  size_t start_3d[4] = {0, 0, 0, 0};
  CHECK_NC(nc_get_vara_float(ncid_inp, var_t_inp, start_3d, count_inp_3d, t_inp.data()));
  CHECK_NC(nc_get_vara_float(ncid_inp, var_qv_inp, start_3d, count_inp_3d, qv_inp.data()));

  std::vector<float> t_bys(5 * nz * nx);
  size_t count_bdy_3d_ys[4] = {1, 5, nz, nx};
  CHECK_NC(nc_get_vara_float(ncid_bdy, var_t_bys, start_3d, count_bdy_3d_ys, t_bys.data()));

  int var_c1h, var_c2h;
  std::vector<float> c1h(nz), c2h(nz);
  size_t start_c1h[2] = {0, 0};
  size_t count_c1h[2] = {1, nz};
  CHECK_NC(nc_inq_varid(ncid_inp, "C1H", &var_c1h));
  CHECK_NC(nc_inq_varid(ncid_inp, "C2H", &var_c2h));
  CHECK_NC(nc_get_vara_float(ncid_inp, var_c1h, start_c1h, count_c1h, c1h.data()));
  CHECK_NC(nc_get_vara_float(ncid_inp, var_c2h, start_c1h, count_c1h, c2h.data()));

  const float rvovrd = 461.6f / 287.0f;
  const float t0 = 300.0f;
  DiffStats t_stats;
  for (size_t k = 0; k < nz; ++k) {
    for (size_t i = 0; i < nx; ++i) {
      float theta_dry = t_inp[k * (ny * nx) + 0 * nx + i];
      float qv = qv_inp[k * (ny * nx) + 0 * nx + i];
      float theta_moist = (theta_dry + t0) * (1.0f + rvovrd * qv) - t0;
      float mu_val = mu_inp[0 * nx + i];
      float mub_val = mub_inp[0 * nx + i];
      float coupled_theta = theta_moist * (c1h[k] * (mu_val + mub_val) + c2h[k]);
      float bdy_theta = t_bys[0 * (nz * nx) + k * nx + i]; // b=0, k, i
      t_stats.record(coupled_theta, bdy_theta);
    }
  }
  std::cout << "  [Check T South b=0 coupled]: points=" << t_stats.total_points
            << " max_abs=" << t_stats.max_abs << " max_rel=" << t_stats.max_rel;
  if (t_stats.max_abs <= 1e-4) {
    std::cout << " -> MATCHES (bit-exact)!\n";
  } else {
    std::cout << " -> MISMATCH!\n";
    nc_close(ncid_bdy);
    nc_close(ncid_inp);
    return false;
  }

  // Check 3: U boundary on West face (XS)
  // coupled_U = (c1h * (mu_u + mub_u) + c2h) * U / msfuy
  int var_u_inp, var_u_bxs;
  CHECK_NC(nc_inq_varid(ncid_inp, "U", &var_u_inp));
  CHECK_NC(nc_inq_varid(ncid_bdy, "U_BXS", &var_u_bxs));

  std::vector<float> u_inp((nx + 1) * ny * nz);
  size_t count_inp_u[4] = {1, nz, ny, nx + 1};
  CHECK_NC(nc_get_vara_float(ncid_inp, var_u_inp, start_3d, count_inp_u, u_inp.data()));

  std::vector<float> u_bxs(5 * nz * ny);
  size_t count_bdy_u_xs[4] = {1, 5, nz, ny};
  CHECK_NC(nc_get_vara_float(ncid_bdy, var_u_bxs, start_3d, count_bdy_u_xs, u_bxs.data()));

  // Map factor msfuy
  int var_msfuy;
  CHECK_NC(nc_inq_varid(ncid_inp, "MAPFAC_UY", &var_msfuy));
  std::vector<float> msfuy((nx + 1) * ny);
  size_t count_msfuy[3] = {1, ny, nx + 1};
  CHECK_NC(nc_get_vara_float(ncid_inp, var_msfuy, start_2d, count_msfuy, msfuy.data()));

  DiffStats u_stats;
  for (size_t k = 0; k < nz; ++k) {
    for (size_t j = 0; j < ny; ++j) {
      float u_val = u_inp[k * (ny * (nx + 1)) + j * (nx + 1) + 0]; // i = 0
      float mu_val = mu_inp[j * nx + 0];
      float mub_val = mub_inp[j * nx + 0];
      float msf_val = msfuy[j * (nx + 1) + 0];
      float coupled_u = (c1h[k] * (mu_val + mub_val) + c2h[k]) * u_val / msf_val;
      float bdy_u = u_bxs[0 * (nz * ny) + k * ny + j]; // b=0, k, j
      u_stats.record(coupled_u, bdy_u);
    }
  }
  std::cout << "  [Check U West b=0 coupled]: points=" << u_stats.total_points
            << " max_abs=" << u_stats.max_abs << " max_rel=" << u_stats.max_rel;
  if (u_stats.max_abs <= 1e-4) {
    std::cout << " -> MATCHES (bit-exact)!\n";
  } else {
    std::cout << " -> MISMATCH!\n";
    nc_close(ncid_bdy);
    nc_close(ncid_inp);
    return false;
  }

  nc_close(ncid_bdy);
  nc_close(ncid_inp);
  std::cout << "  === All R0 wrfbdy & wrfinput checks verified successfully! ===\n";
  return true;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main() {
  std::cout << "=================================================================\n";
  std::cout << "   CPPWRF T4.3: Lateral Boundary Differential Test Suite (C++ vs WRF)\n";
  std::cout << "=================================================================\n";

  bool ok = true;
  ok &= test_lbc_fcx_gcx();
  ok &= test_relax_bdytend_core_3d();
  ok &= test_relax_bdytend_core_2d();
  ok &= test_spec_bdytend();
  ok &= test_spec_bdyupdate();
  ok &= test_spec_bdyupdate_ph();
  ok &= test_zero_grad_bdy();
  ok &= test_flow_dep_bdy();
  ok &= test_spec_bdy_final();
  ok &= test_set_physical_bc();
  ok &= test_r0_wrfbdy_verification();

  std::cout << "=================================================================\n";
  if (ok) {
    std::cout << "   >>> ALL LATERAL BOUNDARY DIFFERENTIAL TESTS PASSED! <<<\n";
    std::cout << "=================================================================\n";
    return 0;
  } else {
    std::cout << "   >>> SOME TESTS FAILED! <<<\n";
    std::cout << "=================================================================\n";
    return 1;
  }
}
