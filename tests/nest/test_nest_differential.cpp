// SPDX-License-Identifier: GPL-3.0-or-later
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#undef NDEBUG  // the checks below are asserts
#include <cassert>

#include "nest/nest_interp.hpp"
#include "nest/nest_bdy.hpp"
#include "nest/nest_feedback.hpp"
#include "nest/nest_init.hpp"

#include "tests/nest/fortran_nest_bridge.h"
#include "tests/nest/test_data_loader.hpp"

using namespace cppwrf::nest;
using namespace cppwrf::nest::test;

// Tolerance 0 (CPU FP64): a value passes only with the same value as WRF's (NaN equals only NaN, -0 differs from
// +0, as tests/verify/exact.py); max |diff| and max rel diff are for the report.
std::size_t g_mismatches = 0;

struct DiffStats {
  double max_abs_diff = 0.0;
  double max_rel_diff = 0.0;
  std::size_t count = 0;
  std::size_t mismatches = 0;

  void update(double cpp_val, double f_val) {
    count++;
    double diff = std::abs(cpp_val - f_val);
    if (diff > max_abs_diff) max_abs_diff = diff;
    double denom = std::max(std::abs(f_val), 1.0e-20);
    double rel = diff / denom;
    if (rel > max_rel_diff) max_rel_diff = rel;
    const bool same = (std::isnan(cpp_val) && std::isnan(f_val)) ||
                      (cpp_val == f_val && std::signbit(cpp_val) == std::signbit(f_val));
    if (!same) {
      mismatches++;
      g_mismatches++;
    }
  }

  bool passed() const { return count > 0 && mismatches == 0; }
};

void run_sint_sintb_tests(const CoarseDataset& ds) {
  std::cout << "\n======================================================\n";
  std::cout << "TEST 1: SINT & SINTB Array Kernels (Ratio 3 and Ratio 5)\n";
  std::cout << "======================================================\n";

  for (int ratio : {3, 5}) {
    for (bool xstag : {false, true}) {
      for (bool ystag : {false, true}) {
        const int nf = ratio * ratio;
        const int ims = 1, ime = 20, jms = 1, jme = 20;
        const int its = 4, ite = 16, jts = 4, jte = 16;
        const int ni = ime - ims + 1;
        const int nj = jme - jms + 1;

        std::vector<double> xf_cpp(ni * nj * nf, 0.0);
        std::vector<double> xf_f(ni * nj * nf, 0.0);
        std::vector<int> mask(ni * nj, 1);

        // Fill with realistic slice from ds.t
        for (int m = 0; m < nf; ++m) {
          for (int j = 0; j < nj; ++j) {
            for (int i = 0; i < ni; ++i) {
              double val = ds.t[i + 5 + 10 * ds.nx + (j + 5) * ds.nx * ds.nz] + (m + 1) * 0.1;
              xf_cpp[i + j * ni + m * ni * nj] = val;
              xf_f[i + j * ni + m * ni * nj] = val;
            }
          }
        }

        // Test sint
        sint_fortran_layout(xf_cpp.data(), ims, ime, jms, jme, mask.data(), its, ite, jts, jte, nf, xstag, ystag);
        wrf_bridge_sint(xf_f.data(), ims, ime, jms, jme, mask.data(), its, ite, jts, jte, nf, xstag ? 1 : 0, ystag ? 1 : 0);

        DiffStats stats_sint;
        for (std::size_t idx = 0; idx < xf_cpp.size(); ++idx) {
          stats_sint.update(xf_cpp[idx], xf_f[idx]);
        }
        std::cout << "  sint ratio=" << ratio << " xstag=" << xstag << " ystag=" << ystag
                  << " -> max_rel=" << stats_sint.max_rel_diff << " max_abs=" << stats_sint.max_abs_diff
                  << " [" << (stats_sint.passed() ? "PASS" : "FAIL") << "]\n";
        assert(stats_sint.passed());

        // Test sintb
        std::vector<double> xf_out_cpp(ni * nj * nf, 0.0);
        std::vector<double> xf_out_f(ni * nj * nf, 0.0);
        sintb_fortran_layout(xf_cpp.data(), xf_out_cpp.data(), ims, ime, jms, jme, mask.data(), its, ite, jts, jte, nf, xstag, ystag);
        wrf_bridge_sintb(xf_cpp.data(), xf_out_f.data(), ims, ime, jms, jme, mask.data(), its, ite, jts, jte, nf, xstag ? 1 : 0, ystag ? 1 : 0);

        DiffStats stats_sintb;
        for (std::size_t idx = 0; idx < xf_out_cpp.size(); ++idx) {
          stats_sintb.update(xf_out_cpp[idx], xf_out_f[idx]);
        }
        std::cout << "  sintb ratio=" << ratio << " xstag=" << xstag << " ystag=" << ystag
                  << " -> max_rel=" << stats_sintb.max_rel_diff << " max_abs=" << stats_sintb.max_abs_diff
                  << " [" << (stats_sintb.passed() ? "PASS" : "FAIL") << "]\n";
        assert(stats_sintb.passed());
      }
    }
  }
}

void run_interp_fcn_sint_tests(const CoarseDataset& ds) {
  std::cout << "\n======================================================\n";
  std::cout << "TEST 2: interp_fcn_sint & KSint3D Device Kernel\n";
  std::cout << "  (ratio 3/5, odd/even starts, parent edge-touching)\n";
  std::cout << "======================================================\n";

  struct Config {
    int ratio;
    int ipos, jpos;
    int n_nx_cells, n_ny_cells;
    std::string desc;
  };

  std::vector<Config> configs = {
    {3, 10, 10, 15, 15, "Ratio 3, interior, even start (10, 10)"},
    {3, 11, 13, 15, 15, "Ratio 3, interior, odd start (11, 13)"},
    {3, 3, 3, 12, 12,   "Ratio 3, touching parent edge (ipos=3, jpos=3)"},
    {5, 10, 10, 10, 10, "Ratio 5, interior, even start (10, 10)"},
    {5, 15, 17, 10, 10, "Ratio 5, interior, odd start (15, 17)"},
    {5, 3, 3, 8, 8,     "Ratio 5, touching parent edge (ipos=3, jpos=3)"},
  };

  for (const auto& cfg : configs) {
    const int ratio = cfg.ratio;
    const int ipos = cfg.ipos;
    const int jpos = cfg.jpos;
    const int c_nx = ds.nx;
    const int c_ny = ds.ny;
    const int c_nz = 4;

    const int n_nx = cfg.n_nx_cells * ratio;
    const int n_ny = cfg.n_ny_cells * ratio;
    const int n_nz = c_nz;

    for (int var_type = 0; var_type < 3; ++var_type) {
      bool xstag = (var_type == 1);
      bool ystag = (var_type == 2);

      const int c_ims = 1, c_ime = c_nx + (xstag ? 1 : 0);
      const int c_kms = 1, c_kme = c_nz;
      const int c_jms = 1, c_jme = c_ny + (ystag ? 1 : 0);

      const int n_ims = 1, n_ime = n_nx + (xstag ? 1 : 0);
      const int n_kms = 1, n_kme = n_nz;
      const int n_jms = 1, n_jme = n_ny + (ystag ? 1 : 0);

      const int c_size = (c_ime - c_ims + 1) * (c_kme - c_kms + 1) * (c_jme - c_jms + 1);
      const int n_size = (n_ime - n_ims + 1) * (n_kme - n_kms + 1) * (n_jme - n_jms + 1);

      const double* cfld_src = (var_type == 1) ? ds.u.data() : (var_type == 2 ? ds.v.data() : ds.t.data());
      std::vector<double> cfld(cfld_src, cfld_src + c_size);

      std::vector<double> nfld_cpp(n_size, 0.0);
      std::vector<double> nfld_f(n_size, 0.0);
      std::vector<double> nfld_dev(n_size, 0.0);

      std::vector<int> imask((n_ime - n_ims + 1) * (n_jme - n_jms + 1), 1);

      const int cits = ipos;
      const int cite = ipos + (n_ime - 1) / ratio;
      const int cjts = jpos;
      const int cjte = jpos + (n_jme - 1) / ratio;

      // 1. Fortran reference
      wrf_bridge_interp_fcn_sint(
          cfld.data(),
          1, c_nx + 1, 1, c_nz + 1, 1, c_ny + 1,
          c_ims, c_ime, c_kms, c_kme, c_jms, c_jme,
          cits, cite, 1, c_nz, cjts, cjte,
          nfld_f.data(),
          1, n_nx + 1, 1, n_nz + 1, 1, n_ny + 1,
          n_ims, n_ime, n_kms, n_kme, n_jms, n_jme,
          1, n_ime, 1, n_kme, 1, n_jme,
          0, imask.data(), xstag ? 1 : 0, ystag ? 1 : 0,
          ipos, jpos, ratio, ratio);

      // 2. C++ host interp_fcn_sint
      interp_fcn_sint(
          cfld.data(),
          1, c_nx + 1, 1, c_nz + 1, 1, c_ny + 1,
          c_ims, c_ime, c_kms, c_kme, c_jms, c_jme,
          cits, cite, 1, c_nz, cjts, cjte,
          nfld_cpp.data(),
          1, n_nx + 1, 1, n_nz + 1, 1, n_ny + 1,
          n_ims, n_ime, n_kms, n_kme, n_jms, n_jme,
          1, n_ime, 1, n_kme, 1, n_jme,
          0, imask.data(), xstag, ystag,
          ipos, jpos, ratio, ratio);

      // 3. C++ single-source KSint3D kernel
      int ioff = 0, joff = 0;
      if (xstag) ioff = (ratio - 1) / 2;
      if (ystag) joff = (ratio - 1) / 2;

      KSint3D kernel{
          cfld.data(),
          nfld_dev.data(),
          imask.data(),
          c_ime - c_ims + 1, c_kme - c_kms + 1, c_jme - c_jms + 1,
          c_ims, c_kms, c_jms,
          n_ime - n_ims + 1, n_kme - n_kms + 1, n_jme - n_jms + 1,
          n_ims, n_kms, n_jms,
          1, n_ime, 1, n_jme,
          ipos, jpos, ratio, ratio,
          ioff, joff,
          xstag, ystag
      };

      for (int k = 1; k <= n_nz; ++k) {
        for (int nj_raw = 1; nj_raw <= n_jme + joff; ++nj_raw) {
          for (int ni_raw = 1; ni_raw <= n_ime + ioff; ++ni_raw) {
            kernel(ni_raw, nj_raw, k);
          }
        }
      }

      // 4. the model's kernel (KSintView3D: src/cases/nest_driver.cpp interpolates PHB and O3RAD, unstaggered,
      //    with it) on the same Fortran-layout arrays (i, k, j)
      std::vector<double> nfld_model(n_size, 0.0);
      if (var_type == 0) {
        const int cni = c_ime - c_ims + 1, cnk = c_kme - c_kms + 1;
        const int nni = n_ime - n_ims + 1, nnk = n_kme - n_kms + 1;
        cppwrf::View3<const double> cv{cfld.data(), cni * cnk, cni};
        cppwrf::View3<double> nv{nfld_model.data(), nni * nnk, nni};
        KSintView3D<cppwrf::View3<const double>, cppwrf::View3<double>> kv{cv, nv, ipos, jpos, ratio, ratio, false,
                                                                           false};
        for (int k = 0; k < n_nz; ++k)
          for (int nj = 0; nj < n_jme; ++nj)
            for (int ni = 0; ni < n_ime; ++ni) kv(ni, nj, k);
      }

      DiffStats stats_cpp, stats_dev, stats_model;
      for (std::size_t idx = 0; idx < nfld_f.size(); ++idx) {
        stats_cpp.update(nfld_cpp[idx], nfld_f[idx]);
        stats_dev.update(nfld_dev[idx], nfld_f[idx]);
        if (var_type == 0) stats_model.update(nfld_model[idx], nfld_f[idx]);
      }
      if (var_type == 0) {
        std::cout << "    -> model KSintView3D vs Fortran: max_abs=" << stats_model.max_abs_diff << " ["
                  << (stats_model.passed() ? "PASS" : "FAIL") << "]\n";
        assert(stats_model.passed());
      }

      std::string vname = (var_type == 0 ? "T (unstaggered)" : (var_type == 1 ? "U (xstag)" : "V (ystag)"));
      std::cout << "  [" << cfg.desc << "] " << vname
                << "\n    -> CPU vs Fortran: max_rel=" << stats_cpp.max_rel_diff
                << " max_abs=" << stats_cpp.max_abs_diff << " [" << (stats_cpp.passed() ? "PASS" : "FAIL") << "]"
                << "\n    -> Dev vs Fortran: max_rel=" << stats_dev.max_rel_diff
                << " max_abs=" << stats_dev.max_abs_diff << " [" << (stats_dev.passed() ? "PASS" : "FAIL") << "]\n";

      assert(stats_cpp.passed());
      assert(stats_dev.passed());
    }
  }
}

void run_bdy_interp_tests(const CoarseDataset& ds) {
  std::cout << "\n======================================================\n";
  std::cout << "TEST 3: Boundary Interpolation (bdy_interp1 & LateralBC slabs)\n";
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

    LateralBCSlabs cpp_slabs = compute_lateral_bc_slabs(
        cfld, ds.nx, ds.ny, n_nz,
        nfld.data(), n_nx, n_ny, n_nz,
        ipos, jpos, ratio, ratio,
        cdt, xstag, ystag,
        spec_bdy_width, 1, 4);

    DiffStats s_xs, s_xe, s_ys, s_ye;
    for (std::size_t i = 0; i < f_xs.size(); ++i) {
      s_xs.update(cpp_slabs.val_xs[i], f_xs[i]);
      s_xs.update(cpp_slabs.tend_xs[i], f_txs[i]);
      s_xe.update(cpp_slabs.val_xe[i], f_xe[i]);
      s_xe.update(cpp_slabs.tend_xe[i], f_txe[i]);
    }
    for (std::size_t i = 0; i < f_ys.size(); ++i) {
      s_ys.update(cpp_slabs.val_ys[i], f_ys[i]);
      s_ys.update(cpp_slabs.tend_ys[i], f_tys[i]);
      s_ye.update(cpp_slabs.val_ye[i], f_ye[i]);
      s_ye.update(cpp_slabs.tend_ye[i], f_tye[i]);
    }

    std::string vname = (var_type == 0 ? "T (unstaggered)" : (var_type == 1 ? "U (xstag)" : "V (ystag)"));
    std::cout << "  bdy_interp1 & LateralBC slabs for " << vname << ":\n"
              << "    -> XS: max_rel=" << s_xs.max_rel_diff << " [" << (s_xs.passed() ? "PASS" : "FAIL") << "]\n"
              << "    -> XE: max_rel=" << s_xe.max_rel_diff << " [" << (s_xe.passed() ? "PASS" : "FAIL") << "]\n"
              << "    -> YS: max_rel=" << s_ys.max_rel_diff << " [" << (s_ys.passed() ? "PASS" : "FAIL") << "]\n"
              << "    -> YE: max_rel=" << s_ye.max_rel_diff << " [" << (s_ye.passed() ? "PASS" : "FAIL") << "]\n";

    assert(s_xs.passed());
    assert(s_xe.passed());
    assert(s_ys.passed());
    assert(s_ye.passed());
  }
}

void run_feedback_tests(const CoarseDataset& ds) {
  std::cout << "\n======================================================\n";
  std::cout << "TEST 4: Feedback (copy_fcn, sm121, smdsm)\n";
  std::cout << "======================================================\n";

  for (int ratio : {3, 4, 5}) {
    const int ipos = 10, jpos = 10;
    const int n_nx = 15 * ratio, n_ny = 15 * ratio, n_nz = 4;

    for (int var_type = 0; var_type < 3; ++var_type) {
      bool xstag = (var_type == 1);
      bool ystag = (var_type == 2);

      const int c_ims = 1, c_ime = ds.nx + (xstag ? 1 : 0);
      const int c_kms = 1, c_kme = n_nz;
      const int c_jms = 1, c_jme = ds.ny + (ystag ? 1 : 0);

      const int n_ims = 1, n_ime = n_nx + (xstag ? 1 : 0);
      const int n_kms = 1, n_kme = n_nz;
      const int n_jms = 1, n_jme = n_ny + (ystag ? 1 : 0);

      const int c_size = (c_ime - c_ims + 1) * (c_kme - c_kms + 1) * (c_jme - c_jms + 1);
      const int n_size = (n_ime - n_ims + 1) * (n_kme - n_kms + 1) * (n_jme - n_jms + 1);

      std::vector<double> nfld(n_size);
      for (std::size_t i = 0; i < nfld.size(); ++i) {
        nfld[i] = 300.0 + std::sin(static_cast<double>(i) * 0.02) * 15.0;
      }

      std::vector<double> cfld_cpp(c_size, 0.0);
      std::vector<double> cfld_f(c_size, 0.0);
      std::vector<int> dummy_mask(n_size, 1);

      // copy_fcn
      copy_fcn(
          cfld_cpp.data(),
          1, ds.nx + 1, 1, n_nz + 1, 1, ds.ny + 1,
          c_ims, c_ime, c_kms, c_kme, c_jms, c_jme,
          1, ds.nx, 1, n_nz, 1, ds.ny,
          nfld.data(),
          1, n_nx + 1, 1, n_nz + 1, 1, n_ny + 1,
          n_ims, n_ime, n_kms, n_kme, n_jms, n_jme,
          1, n_ime, 1, n_kme, 1, n_jme,
          0, dummy_mask.data(),
          xstag, ystag,
          ipos, jpos, ratio, ratio, 1);

      wrf_bridge_copy_fcn(
          cfld_f.data(),
          1, ds.nx + 1, 1, n_nz + 1, 1, ds.ny + 1,
          c_ims, c_ime, c_kms, c_kme, c_jms, c_jme,
          1, ds.nx, 1, n_nz, 1, ds.ny,
          nfld.data(),
          1, n_nx + 1, 1, n_nz + 1, 1, n_ny + 1,
          n_ims, n_ime, n_kms, n_kme, n_jms, n_jme,
          1, n_ime, 1, n_kme, 1, n_jme,
          0, dummy_mask.data(),
          xstag ? 1 : 0, ystag ? 1 : 0,
          ipos, jpos, ratio, ratio);

      DiffStats s_copy;
      for (std::size_t i = 0; i < cfld_f.size(); ++i) s_copy.update(cfld_cpp[i], cfld_f[i]);

      std::string vname = (var_type == 0 ? "T (mass 2D avg)" : (var_type == 1 ? "U (x-face 1D avg)" : "V (y-face 1D avg)"));
      std::cout << "  copy_fcn ratio=" << ratio << " " << vname
                << " -> max_rel=" << s_copy.max_rel_diff << " max_abs=" << s_copy.max_abs_diff
                << " [" << (s_copy.passed() ? "PASS" : "FAIL") << "]\n";
      assert(s_copy.passed());

      if (ratio == 4) continue;

      // sm121 (smooth_option = 1)
      std::vector<double> sm1_cpp = cfld_cpp;
      std::vector<double> sm1_f = cfld_f;
      sm121(sm1_cpp.data(),
            1, ds.nx + 1, 1, n_nz + 1, 1, ds.ny + 1,
            c_ims, c_ime, c_kms, c_kme, c_jms, c_jme,
            1, ds.nx, 1, n_nz, 1, ds.ny,
            xstag, ystag,
            1, n_nx + 1, 1, n_nz + 1, 1, n_ny + 1,
            n_ims, n_ime, n_kms, n_kme, n_jms, n_jme,
            1, n_ime, 1, n_kme, 1, n_jme,
            ratio, ratio, ipos, jpos);

      wrf_bridge_sm121(sm1_f.data(),
                       1, ds.nx + 1, 1, n_nz + 1, 1, ds.ny + 1,
                       c_ims, c_ime, c_kms, c_kme, c_jms, c_jme,
                       1, ds.nx, 1, n_nz, 1, ds.ny,
                       xstag ? 1 : 0, ystag ? 1 : 0,
                       1, n_nx + 1, 1, n_nz + 1, 1, n_ny + 1,
                       n_ims, n_ime, n_kms, n_kme, n_jms, n_jme,
                       1, n_ime, 1, n_kme, 1, n_jme,
                       ratio, ratio, ipos, jpos);

      DiffStats s_sm1;
      for (std::size_t i = 0; i < sm1_f.size(); ++i) s_sm1.update(sm1_cpp[i], sm1_f[i]);
      std::cout << "  sm121 ratio=" << ratio << " " << vname
                << " -> max_rel=" << s_sm1.max_rel_diff << " [" << (s_sm1.passed() ? "PASS" : "FAIL") << "]\n";
      assert(s_sm1.passed());

      // smdsm (smooth_option = 2)
      std::vector<double> sm2_cpp = cfld_cpp;
      std::vector<double> sm2_f = cfld_f;
      smdsm(sm2_cpp.data(),
            1, ds.nx + 1, 1, n_nz + 1, 1, ds.ny + 1,
            c_ims, c_ime, c_kms, c_kme, c_jms, c_jme,
            1, ds.nx, 1, n_nz, 1, ds.ny,
            xstag, ystag,
            1, n_nx + 1, 1, n_nz + 1, 1, n_ny + 1,
            n_ims, n_ime, n_kms, n_kme, n_jms, n_jme,
            1, n_ime, 1, n_kme, 1, n_jme,
            ratio, ratio, ipos, jpos);

      wrf_bridge_smdsm(sm2_f.data(),
                       1, ds.nx + 1, 1, n_nz + 1, 1, ds.ny + 1,
                       c_ims, c_ime, c_kms, c_kme, c_jms, c_jme,
                       1, ds.nx, 1, n_nz, 1, ds.ny,
                       xstag ? 1 : 0, ystag ? 1 : 0,
                       1, n_nx + 1, 1, n_nz + 1, 1, n_ny + 1,
                       n_ims, n_ime, n_kms, n_kme, n_jms, n_jme,
                       1, n_ime, 1, n_kme, 1, n_jme,
                       ratio, ratio, ipos, jpos);

      DiffStats s_sm2;
      for (std::size_t i = 0; i < sm2_f.size(); ++i) s_sm2.update(sm2_cpp[i], sm2_f[i]);
      std::cout << "  smdsm ratio=" << ratio << " " << vname
                << " -> max_rel=" << s_sm2.max_rel_diff << " [" << (s_sm2.passed() ? "PASS" : "FAIL") << "]\n";
      assert(s_sm2.passed());
    }
  }
}

void run_initialisation_tests(const CoarseDataset& ds) {
  std::cout << "\n======================================================\n";
  std::cout << "TEST 6: Nest Initialisation (blend_terrain, adjust_tempqv, copy_3d_field)\n";
  std::cout << "======================================================\n";

  const int nx = ds.nx, ny = ds.ny, nz = 5;

  // 1. blend_terrain
  std::vector<double> ter_int(nx * ny, 200.0);
  std::vector<double> ter_in_cpp = ds.hgt;
  std::vector<double> ter_in_f = ds.hgt;

  blend_terrain(ter_int.data(), ter_in_cpp.data(),
                1, nx + 1, 1, ny + 1, 1, 2,
                1, nx, 1, ny, 1, 1,
                1, nx, 1, ny, 1, 1,
                5, 5);

  wrf_bridge_blend_terrain(ter_int.data(), ter_in_f.data(),
                           1, nx + 1, 1, ny + 1, 1, 2,
                           1, nx, 1, ny, 1, 1,
                           1, nx, 1, ny, 1, 1);

  DiffStats s_blend;
  for (std::size_t i = 0; i < ter_in_f.size(); ++i) s_blend.update(ter_in_cpp[i], ter_in_f[i]);
  std::cout << "  blend_terrain: max_rel=" << s_blend.max_rel_diff << " [" << (s_blend.passed() ? "PASS" : "FAIL") << "]\n";
  assert(s_blend.passed());

  // 2. copy_3d_field
  std::vector<double> copy_dst_cpp(nx * nz * ny, 0.0);
  std::vector<double> copy_dst_f(nx * nz * ny, 0.0);

  copy_3d_field(copy_dst_cpp.data(), ds.t.data(),
                1, nx + 1, 1, ny + 1, 1, nz + 1,
                1, nx, 1, ny, 1, nz,
                1, nx, 1, ny, 1, nz);

  wrf_bridge_copy_3d_field(copy_dst_f.data(), ds.t.data(),
                           1, nx + 1, 1, ny + 1, 1, nz + 1,
                           1, nx, 1, ny, 1, nz,
                           1, nx, 1, ny, 1, nz);

  DiffStats s_copy3d;
  for (std::size_t i = 0; i < copy_dst_f.size(); ++i) s_copy3d.update(copy_dst_cpp[i], copy_dst_f[i]);
  std::cout << "  copy_3d_field: max_rel=" << s_copy3d.max_rel_diff << " [" << (s_copy3d.passed() ? "PASS" : "FAIL") << "]\n";
  assert(s_copy3d.passed());

  // 3. adjust_tempqv (both use_theta_m = 0 and 1)
  for (int use_theta_m : {0, 1}) {
    std::vector<double> save_mub = ds.mub;
    std::vector<double> new_mub = ds.mub;
    for (std::size_t i = 0; i < new_mub.size(); ++i) new_mub[i] += 50.0 * std::cos(static_cast<double>(i) * 0.01);

    std::vector<double> th_cpp(ds.t.begin(), ds.t.begin() + nx * nz * ny);
    std::vector<double> th_f = th_cpp;

    std::vector<double> qv_cpp(ds.qvapor.begin(), ds.qvapor.begin() + nx * nz * ny);
    std::vector<double> qv_f = qv_cpp;

    std::vector<double> pp(nx * nz * ny, 100.0);

    adjust_tempqv(
        new_mub.data(), save_mub.data(),
        ds.c3h.data(), ds.c4h.data(), ds.c1f.data(), ds.p_top,
        th_cpp.data(), pp.data(), qv_cpp.data(), use_theta_m,
        1, nx + 1, 1, ny + 1, 1, nz + 1,
        1, nx, 1, ny, 1, nz,
        1, nx, 1, ny, 1, nz);

    wrf_bridge_adjust_tempqv(
        new_mub.data(), save_mub.data(),
        ds.c3h.data(), ds.c4h.data(), ds.c1f.data(), ds.p_top,
        th_f.data(), pp.data(), qv_f.data(), use_theta_m,
        1, nx + 1, 1, ny + 1, 1, nz + 1,
        1, nx, 1, ny, 1, nz,
        1, nx, 1, ny, 1, nz);

    DiffStats s_th, s_qv;
    for (std::size_t i = 0; i < th_f.size(); ++i) {
      s_th.update(th_cpp[i], th_f[i]);
      s_qv.update(qv_cpp[i], qv_f[i]);
    }

    std::cout << "  adjust_tempqv (use_theta_m=" << use_theta_m << "):\n"
              << "    -> TH: max_rel=" << s_th.max_rel_diff << " [" << (s_th.passed() ? "PASS" : "FAIL") << "]\n"
              << "    -> QV: max_rel=" << s_qv.max_rel_diff << " [" << (s_qv.passed() ? "PASS" : "FAIL") << "]\n";

    assert(s_th.passed());
    assert(s_qv.passed());
  }
}

int main(int argc, char** argv) {
  std::string dataset_path = "coarse_fields.bin";  // run_test.sh passes the path
  if (argc > 1) dataset_path = argv[1];

  std::cout << "======================================================\n";
  std::cout << "CPPWRF nesting numerics: CPU differential test\n";
  std::cout << "======================================================\n";
  std::cout << "Loading coarse dataset from " << dataset_path << "...\n";

  CoarseDataset ds;
  if (!ds.load(dataset_path)) {
    return 1;
  }
  std::cout << "Coarse grid: nx=" << ds.nx << " ny=" << ds.ny << " nz=" << ds.nz
            << " p_top=" << ds.p_top << "\n";

  wrf_bridge_init_constants();

  run_sint_sintb_tests(ds);
  run_interp_fcn_sint_tests(ds);
  run_bdy_interp_tests(ds);
  run_feedback_tests(ds);
  run_initialisation_tests(ds);

  if (g_mismatches) {
    std::cout << "NEST NUMERICS: " << g_mismatches << " values differ from WRF's\n";
    return 1;
  }
  std::cout << "\n======================================================\n";
  std::cout << "ALL NEST NUMERICS CPU DIFFERENTIAL TESTS PASSED (bit-identical)\n";
  std::cout << "======================================================\n";
  return 0;
}
