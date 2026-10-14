#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <cuda_runtime.h>

#include "core/backend.hpp"
#include "dyn/lateral_bc.hpp"
#include "nest/nest_types.hpp"
#include "nest/nest_interp.hpp"
#include "nest/nest_bdy.hpp"
#include "nest/nest_feedback.hpp"
#include "nest/nest_couple.hpp"
#include "nest/nest_init.hpp"
#include "tests/nest/test_data_loader.hpp"
#include "tests/nest/fortran_nest_bridge.h"

using namespace cppwrf;
using namespace cppwrf::nest;
using namespace cppwrf::nest::test;

#define CUDA_CHECK(call) do { \
  cudaError_t err = (call); \
  if (err != cudaSuccess) { \
    std::cerr << "CUDA Error at " << __FILE__ << ":" << __LINE__ << " - " \
              << cudaGetErrorString(err) << std::endl; \
    std::exit(1); \
  } \
} while(0)

struct DiffStats {
  double max_abs_diff = 0.0;
  double max_rel_diff = 0.0;
  std::size_t count = 0;

  void update(double a, double b) {
    double diff = std::abs(a - b);
    double denom = std::max(std::abs(a), std::abs(b));
    double rel = (denom > 1.0e-30) ? (diff / denom) : 0.0;
    if (diff > max_abs_diff) max_abs_diff = diff;
    if (rel > max_rel_diff) max_rel_diff = rel;
    count++;
  }

  bool passed(double tol = 1.0e-12) const {
    return max_rel_diff <= tol;
  }
};

// ============================================================================
// CUDA Kernels
// ============================================================================
__global__ void k_run_sint3d(KSint3D k, int nits, int nite, int njts, int njte, int kts, int kte) {
  int ni_raw = nits + blockIdx.x * blockDim.x + threadIdx.x;
  int nj_raw = njts + blockIdx.y * blockDim.y + threadIdx.y;
  int kk     = kts  + blockIdx.z;

  if (ni_raw <= nite && nj_raw <= njte && kk <= kte) {
    k(ni_raw, nj_raw, kk);
  }
}

__global__ void k_run_bdy_slab(KBdySlab k, int npos, int nk, int width) {
  int pos = blockIdx.x * blockDim.x + threadIdx.x;
  int kk  = blockIdx.y * blockDim.y + threadIdx.y;
  int b   = blockIdx.z * blockDim.z + threadIdx.z;

  if (pos < npos && kk < nk && b < width) {
    k(pos, kk, b);
  }
}

__global__ void k_run_copy_fcn(KCopyFcn k, int ci_start, int ci_end, int cj_start, int cj_end, int ck_start, int ck_end) {
  int ci = ci_start + blockIdx.x * blockDim.x + threadIdx.x;
  int cj = cj_start + blockIdx.y * blockDim.y + threadIdx.y;
  int ck = ck_start + blockIdx.z;

  if (ci <= ci_end && cj <= cj_end && ck <= ck_end) {
    k(ci, cj, ck);
  }
}

__global__ void k_run_couple_uncouple(
    double* u_2, double* v_2, double* w_2, double* ph_2, double* t_2, double* qv,
    const double* mub, const double* mu_2,
    const double* msfuy, const double* msfvx, const double* msfty,
    const double* c1f, const double* c2f, const double* c1h, const double* c2h,
    int nx, int ny, int nz, bool couple) {

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  if (i < nx && j < ny) {
    std::size_t idx2 = i + j * nx;
    double mut_b = mub[idx2];
    double mut_2 = mu_2[idx2];

    for (int k = 0; k < nz; ++k) {
      double denom = (c1h[k] * mut_b + c2h[k]) + (c1h[k] * mut_2);
      double muth = couple ? denom : (1.0 / denom);
      std::size_t idx3 = i + k * nx + j * nx * nz;
      t_2[idx3] *= muth;
      qv[idx3]  *= muth;
    }

    for (int k = 0; k <= nz; ++k) {
      double denom = (c1f[k] * mut_b + c2f[k]) + (c1f[k] * mut_2);
      std::size_t idx3_w = i + k * nx + j * nx * (nz + 1);
      if (couple) {
        double mutf = denom;
        double muwt = denom / msfty[idx2];
        ph_2[idx3_w] *= mutf;
        w_2[idx3_w]  *= muwt;
      } else {
        double mutf = 1.0 / denom;
        double muwt = msfty[idx2] / denom;
        ph_2[idx3_w] *= mutf;
        w_2[idx3_w]  *= muwt;
      }
    }
  }

  // U-component
  if (i < nx + 1 && j < ny) {
    for (int k = 0; k < nz; ++k) {
      double muut;
      if (couple) {
        if (i == 0) {
          double muth = (c1h[k] * mub[0 + j * nx] + c2h[k]) + (c1h[k] * mu_2[0 + j * nx]);
          muut = muth / msfuy[0 + j * (nx + 1)];
        } else if (i == nx) {
          double muth = (c1h[k] * mub[nx - 1 + j * nx] + c2h[k]) + (c1h[k] * mu_2[nx - 1 + j * nx]);
          muut = muth / msfuy[nx + j * (nx + 1)];
        } else {
          double muth_m1 = (c1h[k] * mub[i - 1 + j * nx] + c2h[k]) + (c1h[k] * mu_2[i - 1 + j * nx]);
          double muth_0  = (c1h[k] * mub[i + j * nx] + c2h[k]) + (c1h[k] * mu_2[i + j * nx]);
          muut = 0.5 * (muth_0 + muth_m1) / msfuy[i + j * (nx + 1)];
        }
      } else {
        if (i == 0) {
          double muth_inv = 1.0 / ((c1h[k] * mub[0 + j * nx] + c2h[k]) + (c1h[k] * mu_2[0 + j * nx]));
          muut = msfuy[0 + j * (nx + 1)] * muth_inv;
        } else if (i == nx) {
          muut = msfuy[nx + j * (nx + 1)] / ((c1h[k] * mub[nx - 1 + j * nx] + c2h[k]) + (c1h[k] * mu_2[nx - 1 + j * nx]));
        } else {
          double denom = (c1h[k] * mub[i + j * nx] + c2h[k]) + (c1h[k] * mub[i - 1 + j * nx] + c2h[k]) +
                         (c1h[k] * mu_2[i + j * nx]) + (c1h[k] * mu_2[i - 1 + j * nx]);
          muut = 2.0 * msfuy[i + j * (nx + 1)] / denom;
        }
      }
      std::size_t idx3_u = i + k * (nx + 1) + j * (nx + 1) * nz;
      u_2[idx3_u] *= muut;
    }
  }

  // V-component
  if (i < nx && j < ny + 1) {
    for (int k = 0; k < nz; ++k) {
      double muvt;
      if (couple) {
        if (j == 0) {
          double muth = (c1h[k] * mub[i + 0 * nx] + c2h[k]) + (c1h[k] * mu_2[i + 0 * nx]);
          muvt = muth / msfvx[i + 0 * nx];
        } else if (j == ny) {
          double muth = (c1h[k] * mub[i + (ny - 1) * nx] + c2h[k]) + (c1h[k] * mu_2[i + (ny - 1) * nx]);
          muvt = muth / msfvx[i + ny * nx];
        } else {
          double muth_m1 = (c1h[k] * mub[i + (j - 1) * nx] + c2h[k]) + (c1h[k] * mu_2[i + (j - 1) * nx]);
          double muth_0  = (c1h[k] * mub[i + j * nx] + c2h[k]) + (c1h[k] * mu_2[i + j * nx]);
          muvt = 0.5 * (muth_0 + muth_m1) / msfvx[i + j * nx];
        }
      } else {
        if (j == 0) {
          double muth_inv = 1.0 / ((c1h[k] * mub[i + 0 * nx] + c2h[k]) + (c1h[k] * mu_2[i + 0 * nx]));
          muvt = msfvx[i + 0 * nx] * muth_inv;
        } else if (j == ny) {
          muvt = msfvx[i + ny * nx] / ((c1h[k] * mub[i + (ny - 1) * nx] + c2h[k]) + (c1h[k] * mu_2[i + (ny - 1) * nx]));
        } else {
          double denom = (c1h[k] * mub[i + j * nx] + c2h[k]) + (c1h[k] * mub[i + (j - 1) * nx] + c2h[k]) +
                         (c1h[k] * mu_2[i + j * nx]) + (c1h[k] * mu_2[i + (j - 1) * nx]);
          muvt = 2.0 * msfvx[i + j * nx] / denom;
        }
      }
      std::size_t idx3_v = i + k * nx + j * nx * nz;
      v_2[idx3_v] *= muvt;
    }
  }
}

// ============================================================================
// Test 1: SINT & SINTB CUDA Kernel Execution
// ============================================================================
void run_cuda_sint_tests(const CoarseDataset& ds) {
  std::cout << "\n======================================================\n";
  std::cout << "TEST 1 (CUDA): KSint3D Device Kernel\n";
  std::cout << "======================================================\n";

  struct TestCase {
    int ratio;
    int ipos, jpos;
    int n_nx, n_ny, n_nz;
    const char* desc;
  };

  std::vector<TestCase> cases = {
    {3, 10, 10, 30, 30, 4, "Ratio 3, interior, even start (10, 10)"},
    {3, 11, 13, 33, 33, 4, "Ratio 3, interior, odd start (11, 13)"},
    {3, 3,  3,  27, 27, 4, "Ratio 3, parent edge (ipos=3, jpos=3)"},
    {5, 10, 10, 40, 40, 4, "Ratio 5, interior, even start (10, 10)"},
    {5, 15, 17, 45, 45, 4, "Ratio 5, interior, odd start (15, 17)"},
    {5, 3,  3,  35, 35, 4, "Ratio 5, parent edge (ipos=3, jpos=3)"}
  };

  for (const auto& tc : cases) {
    for (int var_type = 0; var_type < 3; ++var_type) {
      bool xstag = (var_type == 1);
      bool ystag = (var_type == 2);
      const char* vname = (var_type == 0 ? "T (unstaggered)" : (var_type == 1 ? "U (xstag)" : "V (ystag)"));

      const int c_ims = 1, c_ime = ds.nx + (xstag ? 1 : 0);
      const int c_kms = 1, c_kme = tc.n_nz;
      const int c_jms = 1, c_jme = ds.ny + (ystag ? 1 : 0);
      const int c_size = (c_ime - c_ims + 1) * (c_kme - c_kms + 1) * (c_jme - c_jms + 1);

      const int n_ims = 1, n_ime = tc.n_nx + (xstag ? 1 : 0);
      const int n_kms = 1, n_kme = tc.n_nz;
      const int n_jms = 1, n_jme = tc.n_ny + (ystag ? 1 : 0);
      const int n_size = (n_ime - n_ims + 1) * (n_kme - n_kms + 1) * (n_jme - n_jms + 1);

      const double* cfld = (var_type == 1) ? ds.u.data() : (var_type == 2 ? ds.v.data() : ds.t.data());

      std::vector<int> imask((n_ime - n_ims + 1) * (n_jme - n_jms + 1), 1);
      std::vector<double> nfld_f(n_size, 0.0);
      std::vector<double> nfld_dev(n_size, 0.0);

      wrf_bridge_interp_fcn_sint(
          cfld,
          1, ds.nx + 1, 1, tc.n_nz + 1, 1, ds.ny + 1,
          c_ims, c_ime, c_kms, c_kme, c_jms, c_jme,
          tc.ipos, tc.ipos + (n_ime - 1) / tc.ratio, 1, tc.n_nz, tc.jpos, tc.jpos + (n_jme - 1) / tc.ratio,
          nfld_f.data(),
          1, tc.n_nx + 1, 1, tc.n_nz + 1, 1, tc.n_ny + 1,
          n_ims, n_ime, n_kms, n_kme, n_jms, n_jme,
          1, n_ime, 1, tc.n_nz, 1, n_jme,
          0, imask.data(), xstag ? 1 : 0, ystag ? 1 : 0,
          tc.ipos, tc.jpos, tc.ratio, tc.ratio);

      double *d_cfld, *d_nfld;
      int *d_imask;
      CUDA_CHECK(cudaMalloc(&d_cfld, c_size * sizeof(double)));
      CUDA_CHECK(cudaMalloc(&d_nfld, n_size * sizeof(double)));
      CUDA_CHECK(cudaMalloc(&d_imask, imask.size() * sizeof(int)));

      CUDA_CHECK(cudaMemcpy(d_cfld, cfld, c_size * sizeof(double), cudaMemcpyHostToDevice));
      CUDA_CHECK(cudaMemset(d_nfld, 0, n_size * sizeof(double)));
      CUDA_CHECK(cudaMemcpy(d_imask, imask.data(), imask.size() * sizeof(int), cudaMemcpyHostToDevice));

      int ioff = xstag ? std::max((tc.ratio - 1) / 2, 1) : 0;
      int joff = ystag ? std::max((tc.ratio - 1) / 2, 1) : 0;

      KSint3D ksint;
      ksint.cfld = d_cfld;
      ksint.nfld = d_nfld;
      ksint.imask = d_imask;
      ksint.c_ni = c_ime - c_ims + 1;
      ksint.c_nk = c_kme - c_kms + 1;
      ksint.c_nj = c_jme - c_jms + 1;
      ksint.cims = c_ims; ksint.ckms = c_kms; ksint.cjms = c_jms;
      ksint.n_ni = n_ime - n_ims + 1;
      ksint.n_nk = n_kme - n_kms + 1;
      ksint.n_nj = n_jme - n_jms + 1;
      ksint.nims = n_ims; ksint.nkms = n_kms; ksint.njms = n_jms;
      ksint.nits = 1; ksint.nite = n_ime;
      ksint.njts = 1; ksint.njte = n_jme;
      ksint.ipos = tc.ipos; ksint.jpos = tc.jpos;
      ksint.nri = tc.ratio; ksint.nrj = tc.ratio;
      ksint.ioff = ioff; ksint.joff = joff;
      ksint.xstag = xstag; ksint.ystag = ystag;

      dim3 block(16, 16, 1);
      dim3 grid((n_ime + ioff + block.x - 1) / block.x,
                (n_jme + joff + block.y - 1) / block.y,
                tc.n_nz);

      k_run_sint3d<<<grid, block>>>(ksint, 1, n_ime + ioff, 1, n_jme + joff, 1, tc.n_nz);
      CUDA_CHECK(cudaDeviceSynchronize());

      CUDA_CHECK(cudaMemcpy(nfld_dev.data(), d_nfld, n_size * sizeof(double), cudaMemcpyDeviceToHost));

      CUDA_CHECK(cudaFree(d_cfld));
      CUDA_CHECK(cudaFree(d_nfld));
      CUDA_CHECK(cudaFree(d_imask));

      DiffStats s_dev;
      for (std::size_t i = 0; i < nfld_f.size(); ++i) {
        s_dev.update(nfld_dev[i], nfld_f[i]);
      }

      std::cout << "  [" << tc.desc << "] " << vname << ":\n"
                << "    -> CUDA vs Fortran: max_rel=" << s_dev.max_rel_diff
                << " max_abs=" << s_dev.max_abs_diff
                << " [" << (s_dev.passed() ? "PASS" : "FAIL") << "]\n";
      assert(s_dev.passed());
    }
  }
}

// ============================================================================
// Test 2: Boundary Interpolation & LateralBC Slabs on Device
// ============================================================================
void run_cuda_bdy_tests(const CoarseDataset& ds) {
  std::cout << "\n======================================================\n";
  std::cout << "TEST 2 (CUDA): LateralBC Device Slab Generation\n";
  std::cout << "======================================================\n";

  const int ratio = 3;
  const int ipos = 10, jpos = 10;
  const int n_nx = 30, n_ny = 30, n_nz = 4;
  const int spec_bdy_width = 5;
  const double cdt = 60.0;
  const double ndt = cdt / ratio;

  for (int var_type = 0; var_type < 3; ++var_type) {
    bool xstag = (var_type == 1);
    bool ystag = (var_type == 2);
    const char* vname = (var_type == 0 ? "T (unstaggered)" : (var_type == 1 ? "U (xstag)" : "V (ystag)"));

    const int c_ims = 1, c_ime = ds.nx + (xstag ? 1 : 0);
    const int c_kms = 1, c_kme = n_nz;
    const int c_jms = 1, c_jme = ds.ny + (ystag ? 1 : 0);

    const int npos_x = xstag ? n_nx + 1 : n_nx;
    const int npos_y = ystag ? n_ny + 1 : n_ny;

    const double* cfld = (var_type == 1) ? ds.u.data() : (var_type == 2 ? ds.v.data() : ds.t.data());

    std::vector<double> nfld(npos_x * n_nz * npos_y, 280.0);
    for (std::size_t i = 0; i < nfld.size(); ++i) nfld[i] += std::sin(static_cast<double>(i) * 0.05);

    std::vector<double> f_xs(npos_y * n_nz * spec_bdy_width, 0.0), f_txs(npos_y * n_nz * spec_bdy_width, 0.0);
    std::vector<double> f_xe(npos_y * n_nz * spec_bdy_width, 0.0), f_txe(npos_y * n_nz * spec_bdy_width, 0.0);
    std::vector<double> f_ys(npos_x * n_nz * spec_bdy_width, 0.0), f_tys(npos_x * n_nz * spec_bdy_width, 0.0);
    std::vector<double> f_ye(npos_x * n_nz * spec_bdy_width, 0.0), f_tye(npos_x * n_nz * spec_bdy_width, 0.0);

    std::vector<int> dummy_mask(npos_x * npos_y, 0);

    wrf_bridge_bdy_interp1(
        cfld,
        1, ds.nx + 1, 1, n_nz + 1, 1, ds.ny + 1,
        c_ims, c_ime, c_kms, c_kme, c_jms, c_jme,
        1, ds.nx, 1, n_nz, 1, ds.ny,
        nfld.data(), spec_bdy_width,
        1, n_nx + 1, 1, n_nz + 1, 1, n_ny + 1,
        1, npos_x, 1, n_nz, 1, npos_y,
        1, npos_x, 1, n_nz, 1, npos_y,
        0, dummy_mask.data(), xstag ? 1 : 0, ystag ? 1 : 0,
        ipos, jpos, ratio, ratio,
        f_xs.data(), f_txs.data(), f_xe.data(), f_txe.data(),
        f_ys.data(), f_tys.data(), f_ye.data(), f_tye.data(),
        cdt, ndt);

    Grid g;
    g.nx = n_nx; g.ny = n_ny; g.nz = n_nz;
    g.bc_x = BC::Specified; g.bc_y = BC::Specified;
    dyn::LateralBC lbc(g, cdt, spec_bdy_width, 1, 4, true);

    dyn::LateralBC::Var lbc_v = (var_type == 1 ? dyn::LateralBC::U : (var_type == 2 ? dyn::LateralBC::V : dyn::LateralBC::T));

    // Also test KBdySlab directly on GPU!
    const int c_size = (c_ime - c_ims + 1) * (c_kme - c_kms + 1) * (c_jme - c_jms + 1);
    double *d_cfld, *d_nfld;
    CUDA_CHECK(cudaMalloc(&d_cfld, c_size * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_nfld, nfld.size() * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_cfld, cfld, c_size * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_nfld, nfld.data(), nfld.size() * sizeof(double), cudaMemcpyHostToDevice));

    for (int side = 0; side < 4; ++side) {
      auto s_enum = static_cast<dyn::LateralBC::Side>(side);
      int npos = lbc.slab_npos(lbc_v, s_enum);
      KBdySlab kbdy{
        d_cfld, d_nfld,
        lbc.slab(lbc_v, s_enum, false),
        lbc.slab(lbc_v, s_enum, true),
        ds.nx, ds.ny, n_nz,
        n_nx, n_ny, n_nz,
        ipos, jpos, ratio, ratio,
        cdt, xstag, ystag,
        spec_bdy_width, 5, side, npos, n_nz
      };
      dim3 block(16, 4, 1);
      dim3 grid((npos + block.x - 1) / block.x, (n_nz + block.y - 1) / block.y, spec_bdy_width);
      k_run_bdy_slab<<<grid, block>>>(kbdy, npos, n_nz, spec_bdy_width);
    }
    CUDA_CHECK(cudaDeviceSynchronize());
    lbc.mark_set(lbc_v);
    CUDA_CHECK(cudaFree(d_cfld));
    CUDA_CHECK(cudaFree(d_nfld));

    assert(lbc.has(lbc_v));

    std::vector<double> dev_xs(f_xs.size()), dev_txs(f_txs.size());
    std::vector<double> dev_xe(f_xe.size()), dev_txe(f_txe.size());
    std::vector<double> dev_ys(f_ys.size()), dev_tys(f_tys.size());
    std::vector<double> dev_ye(f_ye.size()), dev_tye(f_tye.size());

    copy_to_host(dev_xs.data(), lbc.slab(lbc_v, dyn::LateralBC::XS, false), dev_xs.size() * sizeof(double));
    copy_to_host(dev_txs.data(), lbc.slab(lbc_v, dyn::LateralBC::XS, true), dev_txs.size() * sizeof(double));
    copy_to_host(dev_xe.data(), lbc.slab(lbc_v, dyn::LateralBC::XE, false), dev_xe.size() * sizeof(double));
    copy_to_host(dev_txe.data(), lbc.slab(lbc_v, dyn::LateralBC::XE, true), dev_txe.size() * sizeof(double));
    copy_to_host(dev_ys.data(), lbc.slab(lbc_v, dyn::LateralBC::YS, false), dev_ys.size() * sizeof(double));
    copy_to_host(dev_tys.data(), lbc.slab(lbc_v, dyn::LateralBC::YS, true), dev_tys.size() * sizeof(double));
    copy_to_host(dev_ye.data(), lbc.slab(lbc_v, dyn::LateralBC::YE, false), dev_ye.size() * sizeof(double));
    copy_to_host(dev_tye.data(), lbc.slab(lbc_v, dyn::LateralBC::YE, true), dev_tye.size() * sizeof(double));

    DiffStats s_xs, s_xe, s_ys, s_ye;
    for (std::size_t i = 0; i < f_xs.size(); ++i) {
      s_xs.update(dev_xs[i], f_xs[i]);
      s_xs.update(dev_txs[i], f_txs[i]);
      s_xe.update(dev_xe[i], f_xe[i]);
      s_xe.update(dev_txe[i], f_txe[i]);
    }
    for (std::size_t i = 0; i < f_ys.size(); ++i) {
      s_ys.update(dev_ys[i], f_ys[i]);
      s_ys.update(dev_tys[i], f_tys[i]);
      s_ye.update(dev_ye[i], f_ye[i]);
      s_ye.update(dev_tye[i], f_tye[i]);
    }

    std::cout << "  LateralBC device slabs for " << vname << ":\n"
              << "    -> XS: max_rel=" << s_xs.max_rel_diff << " max_abs=" << s_xs.max_abs_diff << " [" << (s_xs.passed() ? "PASS" : "FAIL") << "]\n"
              << "    -> XE: max_rel=" << s_xe.max_rel_diff << " max_abs=" << s_xe.max_abs_diff << " [" << (s_xe.passed() ? "PASS" : "FAIL") << "]\n"
              << "    -> YS: max_rel=" << s_ys.max_rel_diff << " max_abs=" << s_ys.max_abs_diff << " [" << (s_ys.passed() ? "PASS" : "FAIL") << "]\n"
              << "    -> YE: max_rel=" << s_ye.max_rel_diff << " max_abs=" << s_ye.max_abs_diff << " [" << (s_ye.passed() ? "PASS" : "FAIL") << "]\n";

    assert(s_xs.passed());
    assert(s_xe.passed());
    assert(s_ys.passed());
    assert(s_ye.passed());
  }

  // Also test scalar slab (QVAPOR)
  {
    Grid g;
    g.nx = n_nx; g.ny = n_ny; g.nz = n_nz;
    g.bc_x = BC::Specified; g.bc_y = BC::Specified;
    dyn::LateralBC lbc(g, cdt, spec_bdy_width, 1, 4, true);

    std::vector<double> nfld(n_nx * n_nz * n_ny, 0.015);
    force_nest_lateral_bc_scalar(
        lbc, 0,
        ds.qvapor.data(), ds.nx, ds.ny, n_nz,
        nfld.data(), n_nx, n_ny, n_nz,
        ipos, jpos, ratio, ratio,
        cdt, 1, 4);

    assert(lbc.has_scalar(0));
    std::cout << "  LateralBC scalar slab (QVAPOR): [PASS]\n";
  }
}

// ============================================================================
// Test 3: Feedback Averaging on Device
// ============================================================================
void run_cuda_feedback_tests(const CoarseDataset& ds) {
  std::cout << "\n======================================================\n";
  std::cout << "TEST 3 (CUDA): Feedback Averaging (KCopyFcn)\n";
  std::cout << "======================================================\n";

  for (int ratio : {3, 5}) {
    const int ipos = 10, jpos = 10;
    const int n_nx = 15 * ratio, n_ny = 15 * ratio, n_nz = 4;

    for (int var_type = 0; var_type < 3; ++var_type) {
      bool xstag = (var_type == 1);
      bool ystag = (var_type == 2);
      const char* vname = (var_type == 0 ? "T (mass)" : (var_type == 1 ? "U (x-face)" : "V (y-face)"));

      const int c_ims = 1, c_ime = ds.nx + (xstag ? 1 : 0);
      const int c_kms = 1, c_kme = n_nz;
      const int c_jms = 1, c_jme = ds.ny + (ystag ? 1 : 0);
      const int c_size = (c_ime - c_ims + 1) * (c_kme - c_kms + 1) * (c_jme - c_jms + 1);

      const int n_ims = 1, n_ime = n_nx + (xstag ? 1 : 0);
      const int n_kms = 1, n_kme = n_nz;
      const int n_jms = 1, n_jme = n_ny + (ystag ? 1 : 0);
      const int n_size = (n_ime - n_ims + 1) * (n_kme - n_kms + 1) * (n_jme - n_jms + 1);

      std::vector<double> cfld_f(c_size, 0.0);
      std::vector<double> cfld_dev(c_size, 0.0);

      std::vector<double> nfld(n_size);
      for (std::size_t i = 0; i < nfld.size(); ++i) nfld[i] = 300.0 + std::sin(static_cast<double>(i) * 0.03);

      std::vector<int> dummy_mask((n_ime - n_ims + 1) * (n_jme - n_jms + 1), 0);

      wrf_bridge_copy_fcn(
          cfld_f.data(),
          1, ds.nx + 1, 1, n_nz + 1, 1, ds.ny + 1,
          c_ims, c_ime, c_kms, c_kme, c_jms, c_jme,
          1, ds.nx, 1, n_nz, 1, ds.ny,
          nfld.data(),
          1, n_nx + 1, 1, n_nz + 1, 1, n_ny + 1,
          n_ims, n_ime, n_kms, n_kme, n_jms, n_jme,
          1, n_ime, 1, n_nz, 1, n_jme,
          0, dummy_mask.data(), xstag ? 1 : 0, ystag ? 1 : 0,
          ipos, jpos, ratio, ratio);

      double *d_cfld, *d_nfld;
      CUDA_CHECK(cudaMalloc(&d_cfld, c_size * sizeof(double)));
      CUDA_CHECK(cudaMalloc(&d_nfld, n_size * sizeof(double)));
      CUDA_CHECK(cudaMemset(d_cfld, 0, c_size * sizeof(double)));
      CUDA_CHECK(cudaMemcpy(d_nfld, nfld.data(), n_size * sizeof(double), cudaMemcpyHostToDevice));

      const int istag = xstag ? 0 : 1;
      const int jstag = ystag ? 0 : 1;
      const int cj_start = ipos + 1;
      const int cj_end   = jpos + n_ny / ratio - jstag - 1;
      const int ci_start = ipos + 1;
      const int ci_end   = ipos + n_nx / ratio - istag - 1;

      KCopyFcn kcopy;
      kcopy.cfld = d_cfld;
      kcopy.nfld = d_nfld;
      kcopy.c_ni = c_ime - c_ims + 1;
      kcopy.c_nk = c_kme - c_kms + 1;
      kcopy.c_nj = c_jme - c_jms + 1;
      kcopy.cims = c_ims; kcopy.ckms = c_kms; kcopy.cjms = c_jms;
      kcopy.n_ni = n_ime - n_ims + 1;
      kcopy.n_nk = n_kme - n_kms + 1;
      kcopy.n_nj = n_jme - n_jms + 1;
      kcopy.nims = n_ims; kcopy.nkms = n_kms; kcopy.njms = n_jms;
      kcopy.ci_start = ci_start; kcopy.ci_end = ci_end;
      kcopy.cj_start = cj_start; kcopy.cj_end = cj_end;
      kcopy.ck_start = 1; kcopy.ck_end = n_nz;
      kcopy.ipos = ipos; kcopy.jpos = jpos;
      kcopy.nri = ratio; kcopy.nrj = ratio;
      kcopy.xstag = xstag; kcopy.ystag = ystag;

      dim3 block(16, 16, 1);
      dim3 grid((ci_end - ci_start + 1 + block.x - 1) / block.x,
                (cj_end - cj_start + 1 + block.y - 1) / block.y,
                n_nz);

      k_run_copy_fcn<<<grid, block>>>(kcopy, ci_start, ci_end, cj_start, cj_end, 1, n_nz);
      CUDA_CHECK(cudaDeviceSynchronize());

      CUDA_CHECK(cudaMemcpy(cfld_dev.data(), d_cfld, c_size * sizeof(double), cudaMemcpyDeviceToHost));
      CUDA_CHECK(cudaFree(d_cfld));
      CUDA_CHECK(cudaFree(d_nfld));

      DiffStats s_f;
      for (std::size_t i = 0; i < cfld_f.size(); ++i) {
        s_f.update(cfld_dev[i], cfld_f[i]);
      }

      std::cout << "  KCopyFcn ratio=" << ratio << " " << vname << ": max_rel="
                << s_f.max_rel_diff << " max_abs=" << s_f.max_abs_diff
                << " [" << (s_f.passed() ? "PASS" : "FAIL") << "]\n";
      assert(s_f.passed());
    }
  }
}

// ============================================================================
// Test 4: Coupling & Uncoupling on Device
// ============================================================================
void run_cuda_coupling_tests(const CoarseDataset& ds) {
  std::cout << "\n======================================================\n";
  std::cout << "TEST 4 (CUDA): Coupling & Uncoupling on Device\n";
  std::cout << "======================================================\n";

  const int nx = ds.nx, ny = ds.ny, nz = 5;

  std::vector<double> u = ds.u;
  std::vector<double> v = ds.v;
  std::vector<double> w = ds.w;
  std::vector<double> ph = ds.ph;
  std::vector<double> t = ds.t;
  std::vector<double> qv = ds.qvapor;

  double *d_u, *d_v, *d_w, *d_ph, *d_t, *d_qv;
  double *d_mub, *d_mu, *d_msfuy, *d_msfvx, *d_msfty;
  double *d_c1f, *d_c2f, *d_c1h, *d_c2h;

  CUDA_CHECK(cudaMalloc(&d_u, (nx + 1) * nz * ny * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_v, nx * nz * (ny + 1) * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_w, nx * (nz + 1) * ny * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_ph, nx * (nz + 1) * ny * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_t, nx * nz * ny * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_qv, nx * nz * ny * sizeof(double)));

  CUDA_CHECK(cudaMalloc(&d_mub, nx * ny * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_mu, nx * ny * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_msfuy, (nx + 1) * ny * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_msfvx, nx * (ny + 1) * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_msfty, nx * ny * sizeof(double)));

  CUDA_CHECK(cudaMalloc(&d_c1f, (nz + 1) * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_c2f, (nz + 1) * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_c1h, nz * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_c2h, nz * sizeof(double)));

  CUDA_CHECK(cudaMemcpy(d_u, u.data(), (nx + 1) * nz * ny * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_v, v.data(), nx * nz * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_w, w.data(), nx * (nz + 1) * ny * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_ph, ph.data(), nx * (nz + 1) * ny * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_t, t.data(), nx * nz * ny * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_qv, qv.data(), nx * nz * ny * sizeof(double), cudaMemcpyHostToDevice));

  CUDA_CHECK(cudaMemcpy(d_mub, ds.mub.data(), nx * ny * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_mu, ds.mu.data(), nx * ny * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_msfuy, ds.msfuy.data(), (nx + 1) * ny * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_msfvx, ds.msfvx.data(), nx * (ny + 1) * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_msfty, ds.msfty.data(), nx * ny * sizeof(double), cudaMemcpyHostToDevice));

  CUDA_CHECK(cudaMemcpy(d_c1f, ds.c1f.data(), (nz + 1) * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_c2f, ds.c2f.data(), (nz + 1) * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_c1h, ds.c1h.data(), nz * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_c2h, ds.c2h.data(), nz * sizeof(double), cudaMemcpyHostToDevice));

  dim3 block(16, 16);
  dim3 grid((nx + 1 + block.x - 1) / block.x, (ny + 1 + block.y - 1) / block.y);

  // 1. Couple on GPU
  k_run_couple_uncouple<<<grid, block>>>(
      d_u, d_v, d_w, d_ph, d_t, d_qv,
      d_mub, d_mu, d_msfuy, d_msfvx, d_msfty,
      d_c1f, d_c2f, d_c1h, d_c2h,
      nx, ny, nz, true);
  CUDA_CHECK(cudaDeviceSynchronize());

  std::vector<double> u_f = u, v_f = v, w_f = w, ph_f = ph, t_f = t, qv_f = qv;
  wrf_bridge_couple_or_uncouple_em(
      u_f.data(), v_f.data(), w_f.data(), ph_f.data(), t_f.data(), qv_f.data(),
      ds.mub.data(), ds.mu.data(),
      ds.msfuy.data(), ds.msfvx.data(), ds.msfty.data(),
      ds.c1f.data(), ds.c2f.data(), ds.c1h.data(), ds.c2h.data(),
      nx, ny, nz, 1);

  std::vector<double> u_dev((nx + 1) * nz * ny), v_dev(nx * nz * (ny + 1));
  std::vector<double> t_dev(nx * nz * ny), qv_dev(nx * nz * ny);
  CUDA_CHECK(cudaMemcpy(u_dev.data(), d_u, u_dev.size() * sizeof(double), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(v_dev.data(), d_v, v_dev.size() * sizeof(double), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(t_dev.data(), d_t, t_dev.size() * sizeof(double), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(qv_dev.data(), d_qv, qv_dev.size() * sizeof(double), cudaMemcpyDeviceToHost));

  DiffStats su, sv, st, sqv;
  for (std::size_t i = 0; i < u_dev.size(); ++i) su.update(u_dev[i], u_f[i]);
  for (std::size_t i = 0; i < v_dev.size(); ++i) sv.update(v_dev[i], v_f[i]);
  for (std::size_t i = 0; i < t_dev.size(); ++i) {
    st.update(t_dev[i], t_f[i]);
    sqv.update(qv_dev[i], qv_f[i]);
  }

  std::cout << "  CUDA Couple step:\n"
            << "    -> U:  max_rel=" << su.max_rel_diff  << " max_abs=" << su.max_abs_diff  << " [" << (su.passed() ? "PASS" : "FAIL") << "]\n"
            << "    -> V:  max_rel=" << sv.max_rel_diff  << " max_abs=" << sv.max_abs_diff  << " [" << (sv.passed() ? "PASS" : "FAIL") << "]\n"
            << "    -> T:  max_rel=" << st.max_rel_diff  << " max_abs=" << st.max_abs_diff  << " [" << (st.passed() ? "PASS" : "FAIL") << "]\n"
            << "    -> QV: max_rel=" << sqv.max_rel_diff << " max_abs=" << sqv.max_abs_diff << " [" << (sqv.passed() ? "PASS" : "FAIL") << "]\n";
  assert(su.passed());
  assert(sv.passed());
  assert(st.passed());
  assert(sqv.passed());

  // 2. Uncouple on GPU
  k_run_couple_uncouple<<<grid, block>>>(
      d_u, d_v, d_w, d_ph, d_t, d_qv,
      d_mub, d_mu, d_msfuy, d_msfvx, d_msfty,
      d_c1f, d_c2f, d_c1h, d_c2h,
      nx, ny, nz, false);
  CUDA_CHECK(cudaDeviceSynchronize());

  wrf_bridge_couple_or_uncouple_em(
      u_f.data(), v_f.data(), w_f.data(), ph_f.data(), t_f.data(), qv_f.data(),
      ds.mub.data(), ds.mu.data(),
      ds.msfuy.data(), ds.msfvx.data(), ds.msfty.data(),
      ds.c1f.data(), ds.c2f.data(), ds.c1h.data(), ds.c2h.data(),
      nx, ny, nz, 0);

  CUDA_CHECK(cudaMemcpy(u_dev.data(), d_u, u_dev.size() * sizeof(double), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(v_dev.data(), d_v, v_dev.size() * sizeof(double), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(t_dev.data(), d_t, t_dev.size() * sizeof(double), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(qv_dev.data(), d_qv, qv_dev.size() * sizeof(double), cudaMemcpyDeviceToHost));

  DiffStats su2, sv2, st2, sqv2;
  for (std::size_t i = 0; i < u_dev.size(); ++i) su2.update(u_dev[i], u_f[i]);
  for (std::size_t i = 0; i < v_dev.size(); ++i) sv2.update(v_dev[i], v_f[i]);
  for (std::size_t i = 0; i < t_dev.size(); ++i) {
    st2.update(t_dev[i], t_f[i]);
    sqv2.update(qv_dev[i], qv_f[i]);
  }

  std::cout << "  CUDA Uncouple step:\n"
            << "    -> U:  max_rel=" << su2.max_rel_diff  << " max_abs=" << su2.max_abs_diff  << " [" << (su2.passed() ? "PASS" : "FAIL") << "]\n"
            << "    -> V:  max_rel=" << sv2.max_rel_diff  << " max_abs=" << sv2.max_abs_diff  << " [" << (sv2.passed() ? "PASS" : "FAIL") << "]\n"
            << "    -> T:  max_rel=" << st2.max_rel_diff  << " max_abs=" << st2.max_abs_diff  << " [" << (st2.passed() ? "PASS" : "FAIL") << "]\n"
            << "    -> QV: max_rel=" << sqv2.max_rel_diff << " max_abs=" << sqv2.max_abs_diff << " [" << (sqv2.passed() ? "PASS" : "FAIL") << "]\n";
  assert(su2.passed());
  assert(sv2.passed());
  assert(st2.passed());
  assert(sqv2.passed());

  CUDA_CHECK(cudaFree(d_u));
  CUDA_CHECK(cudaFree(d_v));
  CUDA_CHECK(cudaFree(d_w));
  CUDA_CHECK(cudaFree(d_ph));
  CUDA_CHECK(cudaFree(d_t));
  CUDA_CHECK(cudaFree(d_qv));
  CUDA_CHECK(cudaFree(d_mub));
  CUDA_CHECK(cudaFree(d_mu));
  CUDA_CHECK(cudaFree(d_msfuy));
  CUDA_CHECK(cudaFree(d_msfvx));
  CUDA_CHECK(cudaFree(d_msfty));
  CUDA_CHECK(cudaFree(d_c1f));
  CUDA_CHECK(cudaFree(d_c2f));
  CUDA_CHECK(cudaFree(d_c1h));
  CUDA_CHECK(cudaFree(d_c2h));
}

// ============================================================================
// Test 5: Initialisation Kernels on Device
// ============================================================================
void run_cuda_initialisation_tests(const CoarseDataset& ds) {
  std::cout << "\n======================================================\n";
  std::cout << "TEST 5 (CUDA): Initialisation (blend_terrain, adjust_tempqv)\n";
  std::cout << "======================================================\n";

  const int nx = ds.nx, ny = ds.ny, nz = 5;

  // 1. blend_terrain
  std::vector<double> ter_int(nx * ny, 200.0);
  std::vector<double> ter_in_f = ds.hgt;
  std::vector<double> ter_in_dev = ds.hgt;

  wrf_bridge_blend_terrain(ter_int.data(), ter_in_f.data(),
                           1, nx + 1, 1, ny + 1, 1, 2,
                           1, nx, 1, ny, 1, 1,
                           1, nx, 1, ny, 1, 1);

  blend_terrain(ter_int.data(), ter_in_dev.data(),
                1, nx + 1, 1, ny + 1, 1, 2,
                1, nx, 1, ny, 1, 1,
                1, nx, 1, ny, 1, 1,
                5, 5);

  DiffStats s_blend;
  for (std::size_t i = 0; i < ter_in_f.size(); ++i) s_blend.update(ter_in_dev[i], ter_in_f[i]);
  std::cout << "  blend_terrain: max_rel=" << s_blend.max_rel_diff << " [" << (s_blend.passed() ? "PASS" : "FAIL") << "]\n";
  assert(s_blend.passed());

  // 2. adjust_tempqv (use_theta_m = 0 and 1)
  for (int use_theta_m : {0, 1}) {
    std::vector<double> save_mub = ds.mub;
    std::vector<double> new_mub = ds.mub;
    for (std::size_t i = 0; i < new_mub.size(); ++i) new_mub[i] += 50.0 * std::cos(static_cast<double>(i) * 0.01);

    std::vector<double> th_f(ds.t.begin(), ds.t.begin() + nx * nz * ny);
    std::vector<double> th_dev = th_f;
    std::vector<double> qv_f(ds.qvapor.begin(), ds.qvapor.begin() + nx * nz * ny);
    std::vector<double> qv_dev = qv_f;
    std::vector<double> pp(nx * nz * ny, 100.0);

    wrf_bridge_adjust_tempqv(
        new_mub.data(), save_mub.data(),
        ds.c3h.data(), ds.c4h.data(), ds.c1f.data(), ds.p_top,
        th_f.data(), pp.data(), qv_f.data(), use_theta_m,
        1, nx + 1, 1, ny + 1, 1, nz + 1,
        1, nx, 1, ny, 1, nz,
        1, nx, 1, ny, 1, nz);

    adjust_tempqv(
        new_mub.data(), save_mub.data(),
        ds.c3h.data(), ds.c4h.data(), ds.c1f.data(), ds.p_top,
        th_dev.data(), pp.data(), qv_dev.data(), use_theta_m,
        1, nx + 1, 1, ny + 1, 1, nz + 1,
        1, nx, 1, ny, 1, nz,
        1, nx, 1, ny, 1, nz);

    DiffStats s_th, s_qv;
    for (std::size_t i = 0; i < th_f.size(); ++i) {
      s_th.update(th_dev[i], th_f[i]);
      s_qv.update(qv_dev[i], qv_f[i]);
    }

    std::cout << "  adjust_tempqv (use_theta_m=" << use_theta_m << "):\n"
              << "    -> TH: max_rel=" << s_th.max_rel_diff << " [" << (s_th.passed() ? "PASS" : "FAIL") << "]\n"
              << "    -> QV: max_rel=" << s_qv.max_rel_diff << " [" << (s_qv.passed() ? "PASS" : "FAIL") << "]\n";
    assert(s_th.passed());
    assert(s_qv.passed());
  }
}

int main(int argc, char** argv) {
  std::string dataset_path = "/home/ostriquetrum/cppwrf_scratch/nest_test/coarse_fields.bin";
  if (argc > 1) {
    dataset_path = argv[1];
  }

  std::cout << "======================================================\n";
  std::cout << "CPPWRF Nesting Numerics CUDA FP64 Test\n";
  std::cout << "======================================================\n";
  std::cout << "Loading coarse dataset from " << dataset_path << "...\n";

  CoarseDataset ds;
  if (!ds.load(dataset_path)) {
    return 1;
  }
  std::cout << "Coarse grid: nx=" << ds.nx << " ny=" << ds.ny << " nz=" << ds.nz
            << " p_top=" << ds.p_top << "\n";

  int device_count = 0;
  CUDA_CHECK(cudaGetDeviceCount(&device_count));
  if (device_count == 0) {
    std::cerr << "No CUDA devices found!\n";
    return 1;
  }
  cudaDeviceProp prop;
  CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
  std::cout << "Using GPU Device 0: " << prop.name
            << " (Compute Capability " << prop.major << "." << prop.minor
            << ", " << (prop.totalGlobalMem / (1024 * 1024)) << " MB)\n";

  wrf_bridge_init_constants();

  run_cuda_sint_tests(ds);
  run_cuda_bdy_tests(ds);
  run_cuda_feedback_tests(ds);
  run_cuda_coupling_tests(ds);
  run_cuda_initialisation_tests(ds);

  std::cout << "\n======================================================\n";
  std::cout << "ALL NEST NUMERICS CUDA TESTS PASSED! <= 1e-12\n";
  std::cout << "======================================================\n";

  return 0;
}
