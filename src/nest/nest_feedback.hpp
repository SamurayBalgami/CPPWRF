#pragma once

// ============================================================================
// CPPWRF Feedback Averaging & Smoothing Kernels (copy_fcn, sm121, smdsm)
// Ports WRF 4.6.1 share/interp_fcn.F (copy_fcn, smoother, sm121, smdsm)
// Single-source C++20 / CUDA FP64 bit-exact numerics for odd and even ratios.
// ============================================================================

#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdint>
#include "core/backend.hpp"
#include "core/field.hpp"
#include "nest/nest_types.hpp"

namespace cppwrf::nest {

// ============================================================================
// copy_fcn: Feedback averaging from fine to coarse grid (odd & even ratios)
// Matches WRF share/interp_fcn.F lines 1397-1740 bit-for-bit (term-by-term
// accumulation of inv_pts * nfld(...) in ijpoints loop order).
// Boundary zone (spec_zone) on coarse grid footprint is skipped.
// ============================================================================
inline void copy_fcn(
    double* cfld,
    int cids, int cide, int ckds, int ckde, int cjds, int cjde,
    int cims, int cime, int ckms, int ckme, int cjms, int cjme,
    int cits, int cite, int ckts, int ckte, int cjts, int cjte,
    const double* nfld,
    int nids, int nide, int nkds, int nkde, int njds, int njde,
    int nims, int nime, int nkms, int nkme, int njms, int njme,
    int nits, int nite, int nkts, int nkte, int njts, int njte,
    int /*shw*/,
    const int* /*imask*/,
    bool xstag, bool ystag,
    int ipos, int jpos,
    int nri, int nrj,
    int spec_zone = 1) {

  const int c_ni = cime - cims + 1;
  const int c_nk = ckme - ckms + 1;
  auto c_idx = [&](int i, int k, int j) -> std::size_t {
    return static_cast<std::size_t>((i - cims) + (k - ckms) * c_ni + (j - cjms) * c_ni * c_nk);
  };

  const int n_ni = nime - nims + 1;
  const int n_nk = nkme - nkms + 1;
  auto n_idx = [&](int i, int k, int j) -> std::size_t {
    return static_cast<std::size_t>((i - nims) + (k - nkms) * n_ni + (j - njms) * n_ni * n_nk);
  };

  const int istag = xstag ? 0 : 1;
  const int jstag = ystag ? 0 : 1;

  const int cj_start = std::max(jpos + spec_zone, cjts);
  const int cj_end   = std::min(jpos + (njde - njds) / nrj - jstag - spec_zone, cjte);
  const int ci_start = std::max(ipos + spec_zone, cits);
  const int ci_end   = std::min(ipos + (nide - nids) / nri - istag - spec_zone, cite);

  if (!xstag && !ystag) {
    // Unstaggered (Mass / T / PH / W / scalars): nri x nrj cell average
    const double inv_pts = 1.0 / static_cast<double>(nri * nrj);
    for (int cj = cj_start; cj <= cj_end; ++cj) {
      const int nj_base = (cj - jpos) * nrj + 1;
      for (int ck = ckts; ck <= ckte; ++ck) {
        const int nk = ck;
        for (int ci = ci_start; ci <= ci_end; ++ci) {
          const int ni_base = (ci - ipos) * nri + 1;
          double sum = 0.0;
          for (int dj = 0; dj < nrj; ++dj) {
            for (int di = 0; di < nri; ++di) {
              sum += inv_pts * nfld[n_idx(ni_base + di, nk, nj_base + dj)];
            }
          }
          cfld[c_idx(ci, ck, cj)] = sum;
        }
      }
    }
  } else if (xstag && !ystag) {
    // X-staggered (U momentum): 1 x nrj face average
    const double inv_pts = 1.0 / static_cast<double>(nrj);
    for (int cj = cj_start; cj <= cj_end; ++cj) {
      const int nj_base = (cj - jpos) * nrj + 1;
      for (int ck = ckts; ck <= ckte; ++ck) {
        const int nk = ck;
        for (int ci = ci_start; ci <= ci_end; ++ci) {
          const int ni_face = (ci - ipos) * nri + 1;
          double sum = 0.0;
          for (int dj = 0; dj < nrj; ++dj) {
            sum += inv_pts * nfld[n_idx(ni_face, nk, nj_base + dj)];
          }
          cfld[c_idx(ci, ck, cj)] = sum;
        }
      }
    }
  } else if (!xstag && ystag) {
    // Y-staggered (V momentum): nri x 1 face average
    const double inv_pts = 1.0 / static_cast<double>(nri);
    for (int cj = cj_start; cj <= cj_end; ++cj) {
      const int nj_face = (cj - jpos) * nrj + 1;
      for (int ck = ckts; ck <= ckte; ++ck) {
        const int nk = ck;
        for (int ci = ci_start; ci <= ci_end; ++ci) {
          const int ni_base = (ci - ipos) * nri + 1;
          double sum = 0.0;
          for (int di = 0; di < nri; ++di) {
            sum += inv_pts * nfld[n_idx(ni_base + di, nk, nj_face)];
          }
          cfld[c_idx(ci, ck, cj)] = sum;
        }
      }
    }
  }
}

// ============================================================================
// Single-Source Device Functor KCopyFcn (raw Fortran-layout arrays)
// Direct parallel evaluation of feedback averaging on device/CPU (ci, cj, ck)
// ============================================================================
struct KCopyFcn {
  double* cfld;
  const double* nfld;
  int c_ni, c_nk, c_nj;
  int cims, ckms, cjms;
  int n_ni, n_nk, n_nj;
  int nims, nkms, njms;
  int ci_start, ci_end, cj_start, cj_end, ck_start, ck_end;
  int ipos, jpos, nri, nrj;
  bool xstag, ystag;

  CPPWRF_HD void operator()(int ci, int cj, int ck) const {
    if (ci < ci_start || ci > ci_end || cj < cj_start || cj > cj_end || ck < ck_start || ck > ck_end) return;
    const int nk = ck;

    std::size_t c_idx = static_cast<std::size_t>((ci - cims) + (ck - ckms) * c_ni + (cj - cjms) * c_ni * c_nk);

    if (!xstag && !ystag) {
      const double inv_pts = 1.0 / static_cast<double>(nri * nrj);
      const int nj_base = (cj - jpos) * nrj + 1;
      const int ni_base = (ci - ipos) * nri + 1;
      double sum = 0.0;
      for (int dj = 0; dj < nrj; ++dj) {
        for (int di = 0; di < nri; ++di) {
          std::size_t n_idx = static_cast<std::size_t>(
              (ni_base + di - nims) + (nk - nkms) * n_ni + (nj_base + dj - njms) * n_ni * n_nk);
          sum += inv_pts * nfld[n_idx];
        }
      }
      cfld[c_idx] = sum;
    } else if (xstag && !ystag) {
      const double inv_pts = 1.0 / static_cast<double>(nrj);
      const int nj_base = (cj - jpos) * nrj + 1;
      const int ni_face = (ci - ipos) * nri + 1;
      double sum = 0.0;
      for (int dj = 0; dj < nrj; ++dj) {
        std::size_t n_idx = static_cast<std::size_t>(
            (ni_face - nims) + (nk - nkms) * n_ni + (nj_base + dj - njms) * n_ni * n_nk);
        sum += inv_pts * nfld[n_idx];
      }
      cfld[c_idx] = sum;
    } else if (!xstag && ystag) {
      const double inv_pts = 1.0 / static_cast<double>(nri);
      const int nj_face = (cj - jpos) * nrj + 1;
      const int ni_base = (ci - ipos) * nri + 1;
      double sum = 0.0;
      for (int di = 0; di < nri; ++di) {
        std::size_t n_idx = static_cast<std::size_t>(
            (ni_base + di - nims) + (nk - nkms) * n_ni + (nj_face - njms) * n_ni * n_nk);
        sum += inv_pts * nfld[n_idx];
      }
      cfld[c_idx] = sum;
    }
  }
};

// ============================================================================
// Single-Source View-based KCopyFcnView3D / KCopyFcnView2D for Field3 / Field2
// 0-based coarse indices ci0, cj0, k
// ============================================================================
template <class CView, class NView>
struct KCopyFcnView3D {
  CView cfld;
  NView nfld;
  int ipos, jpos, nri, nrj;
  bool xstag, ystag;

  CPPWRF_HD void operator()(int ci0, int cj0, int k) const {
    const int ci = ci0 + 1;
    const int cj = cj0 + 1;
    const int rx2 = nri / 2;
    const int ry2 = nrj / 2;

    if (!xstag && !ystag) {
      const double inv_pts = 1.0 / static_cast<double>(nri * nrj);
      const int nj_center = (cj - jpos) * nrj + ry2 + 1;
      const int ni_center = (ci - ipos) * nri + rx2 + 1;
      double sum = 0.0;
      for (int dj = -ry2; dj <= ry2; ++dj) {
        for (int di = -rx2; di <= rx2; ++di) {
          sum += inv_pts * static_cast<double>(nfld(ni_center + di - 1, nj_center + dj - 1, k));
        }
      }
      cfld(ci0, cj0, k) = static_cast<real>(sum);
    } else if (xstag && !ystag) {
      const double inv_pts = 1.0 / static_cast<double>(nrj);
      const int nj_center = (cj - jpos) * nrj + ry2 + 1;
      const int ni_face = (ci - ipos) * nri + 1;
      double sum = 0.0;
      for (int dj = -ry2; dj <= ry2; ++dj) {
        sum += inv_pts * static_cast<double>(nfld(ni_face - 1, nj_center + dj - 1, k));
      }
      cfld(ci0, cj0, k) = static_cast<real>(sum);
    } else if (!xstag && ystag) {
      const double inv_pts = 1.0 / static_cast<double>(nri);
      const int nj_face = (cj - jpos) * nrj + 1;
      const int ni_center = (ci - ipos) * nri + rx2 + 1;
      double sum = 0.0;
      for (int di = -rx2; di <= rx2; ++di) {
        sum += inv_pts * static_cast<double>(nfld(ni_center + di - 1, nj_face - 1, k));
      }
      cfld(ci0, cj0, k) = static_cast<real>(sum);
    }
  }
};

template <class CView, class NView>
struct KCopyFcnView2D {
  CView cfld;
  NView nfld;
  int ipos, jpos, nri, nrj;

  CPPWRF_HD void operator()(int ci0, int cj0, int /*k*/) const {
    const int ci = ci0 + 1;
    const int cj = cj0 + 1;
    const int rx2 = nri / 2;
    const int ry2 = nrj / 2;
    const double inv_pts = 1.0 / static_cast<double>(nri * nrj);
    const int nj_center = (cj - jpos) * nrj + ry2 + 1;
    const int ni_center = (ci - ipos) * nri + rx2 + 1;
    double sum = 0.0;
    for (int dj = -ry2; dj <= ry2; ++dj) {
      for (int di = -rx2; di <= rx2; ++di) {
        sum += inv_pts * static_cast<double>(nfld(ni_center + di - 1, nj_center + dj - 1));
      }
    }
    cfld(ci0, cj0) = static_cast<real>(sum);
  }
};

// ============================================================================
// sm121: 1-2-1 smoothing filter on coarse grid (smooth_option = 1)
// Applied in j first, then in i.
// ============================================================================
inline void sm121(
    double* cfld,
    int cids, int cide, int ckds, int ckde, int cjds, int cjde,
    int cims, int cime, int ckms, int ckme, int cjms, int cjme,
    int cits, int cite, int ckts, int ckte, int cjts, int cjte,
    bool xstag, bool ystag,
    int nids, int nide, int nkds, int nkde, int njds, int njde,
    int nims, int nime, int nkms, int nkme, int njms, int njme,
    int nits, int nite, int nkts, int nkte, int njts, int njte,
    int nri, int nrj,
    int ipos, int jpos) {

  const int c_ni = cime - cims + 1;
  const int c_nk = ckme - ckms + 1;
  auto c_idx = [&](int i, int k, int j) -> std::size_t {
    return static_cast<std::size_t>((i - cims) + (k - ckms) * c_ni + (j - cjms) * c_ni * c_nk);
  };
  const int c_nj = cjme - cjms + 1;
  auto new_idx = [&](int i, int j) -> std::size_t {
    return static_cast<std::size_t>((i - cims) + (j - cjms) * c_ni);
  };

  const int istag = xstag ? 0 : 1;
  const int jstag = ystag ? 0 : 1;

  const int i_min = std::max(ipos + 2, cits - 2);
  const int i_max = std::min(ipos + (nide - nids) / nri - 2 - istag, cite + 2);
  const int j_min = std::max(jpos + 2, cjts - 2);
  const int j_max = std::min(jpos + (njde - njds) / nrj - 2 - jstag, cjte + 2);

  std::vector<double> cfldnew(static_cast<std::size_t>(c_ni * c_nj), 0.0);

  for (int k = ckts; k <= ckte; ++k) {
    // Initialize cfldnew from cfld
    for (int i = std::max(ipos, cits - 3); i <= std::min(ipos + (nide - nids) / nri, cite + 3); ++i) {
      for (int j = std::max(jpos, cjts - 3); j <= std::min(jpos + (njde - njds) / nrj, cjte + 3); ++j) {
        cfldnew[new_idx(i, j)] = cfld[c_idx(i, k, j)];
      }
    }

    // 1-2-1 smoothing in j
    for (int i = i_min; i <= i_max; ++i) {
      for (int j = j_min; j <= j_max; ++j) {
        cfldnew[new_idx(i, j)] = 0.25 * (cfld[c_idx(i, k, j + 1)] + 2.0 * cfld[c_idx(i, k, j)] + cfld[c_idx(i, k, j - 1)]);
      }
    }
    // 1-2-1 smoothing in i
    for (int j = j_min; j <= j_max; ++j) {
      for (int i = i_min; i <= i_max; ++i) {
        cfld[c_idx(i, k, j)] = 0.25 * (cfldnew[new_idx(i + 1, j)] + 2.0 * cfldnew[new_idx(i, j)] + cfldnew[new_idx(i - 1, j)]);
      }
    }
  }
}

// ============================================================================
// smdsm: Smoother-desmoother filter on coarse grid (smooth_option = 2)
// Pass 1: nu = 0.50 (smoothing)
// Pass 2: nu = -0.52 (desmoothing)
// ============================================================================
inline void smdsm(
    double* cfld,
    int cids, int cide, int ckds, int ckde, int cjds, int cjde,
    int cims, int cime, int ckms, int ckme, int cjms, int cjme,
    int cits, int cite, int ckts, int ckte, int cjts, int cjte,
    bool xstag, bool ystag,
    int nids, int nide, int nkds, int nkde, int njds, int njde,
    int nims, int nime, int nkms, int nkme, int njms, int njme,
    int nits, int nite, int nkts, int nkte, int njts, int njte,
    int nri, int nrj,
    int ipos, int jpos) {

  const int c_ni = cime - cims + 1;
  const int c_nk = ckme - ckms + 1;
  auto c_idx = [&](int i, int k, int j) -> std::size_t {
    return static_cast<std::size_t>((i - cims) + (k - ckms) * c_ni + (j - cjms) * c_ni * c_nk);
  };
  const int c_nj = cjme - cjms + 1;
  auto new_idx = [&](int i, int j) -> std::size_t {
    return static_cast<std::size_t>((i - cims) + (j - cjms) * c_ni);
  };

  const int istag = xstag ? 0 : 1;
  const int jstag = ystag ? 0 : 1;

  const int i_min = std::max(ipos + 2, cits - 2);
  const int i_max = std::min(ipos + (nide - nids) / nri - 2 - istag, cite + 2);
  const int j_min = std::max(jpos + 2, cjts - 2);
  const int j_max = std::min(jpos + (njde - njds) / nrj - 2 - jstag, cjte + 2);

  const double xnu[2] = {0.50, -0.52};
  std::vector<double> cfldnew(static_cast<std::size_t>(c_ni * c_nj), 0.0);

  for (int loop = 1; loop <= 2; ++loop) {
    const double nu = (loop == 1) ? xnu[0] : xnu[1];

    for (int k = ckts; k <= ckte; ++k) {
      // Sweep in j
      for (int j = j_min; j <= j_max; ++j) {
        for (int i = i_min; i <= i_max; ++i) {
          cfldnew[new_idx(i, j)] = cfld[c_idx(i, k, j)] +
              nu * ((cfld[c_idx(i, k, j + 1)] + cfld[c_idx(i, k, j - 1)]) * 0.5 - cfld[c_idx(i, k, j)]);
        }
      }
      for (int j = j_min; j <= j_max; ++j) {
        for (int i = i_min; i <= i_max; ++i) {
          cfld[c_idx(i, k, j)] = cfldnew[new_idx(i, j)];
        }
      }

      // Sweep in i
      for (int j = j_min; j <= j_max; ++j) {
        for (int i = i_min; i <= i_max; ++i) {
          cfldnew[new_idx(i, j)] = cfld[c_idx(i, k, j)] +
              nu * ((cfld[c_idx(i + 1, k, j)] + cfld[c_idx(i - 1, k, j)]) * 0.5 - cfld[c_idx(i, k, j)]);
        }
      }
      for (int j = j_min; j <= j_max; ++j) {
        for (int i = i_min; i <= i_max; ++i) {
          cfld[c_idx(i, k, j)] = cfldnew[new_idx(i, j)];
        }
      }
    }
  }
}

// Wrapper for smoother matching WRF share/interp_fcn.F smoother subroutine
inline void smoother(
    double* cfld,
    int cids, int cide, int ckds, int ckde, int cjds, int cjde,
    int cims, int cime, int ckms, int ckme, int cjms, int cjme,
    int cits, int cite, int ckts, int ckte, int cjts, int cjte,
    int nids, int nide, int nkds, int nkde, int njds, int njde,
    int nims, int nime, int nkms, int nkme, int njms, int njme,
    int nits, int nite, int nkts, int nkte, int njts, int njte,
    bool xstag, bool ystag,
    int ipos, int jpos,
    int nri, int nrj,
    int smooth_option) {

  if (smooth_option == 1) {
    sm121(cfld, cids, cide, ckds, ckde, cjds, cjde,
          cims, cime, ckms, ckme, cjms, cjme,
          cits, cite, ckts, ckte, cjts, cjte,
          xstag, ystag,
          nids, nide, nkds, nkde, njds, njde,
          nims, nime, nkms, nkme, njms, njme,
          nits, nite, nkts, nkte, njts, njte,
          nri, nrj, ipos, jpos);
  } else if (smooth_option == 2) {
    smdsm(cfld, cids, cide, ckds, ckde, cjds, cjde,
          cims, cime, ckms, ckme, cjms, cjme,
          cits, cite, ckts, ckte, cjts, cjte,
          xstag, ystag,
          nids, nide, nkds, nkde, njds, njde,
          nims, nime, nkms, nkme, njms, njme,
          nits, nite, nkts, nkte, njts, njte,
          nri, nrj, ipos, jpos);
  }
}

}  // namespace cppwrf::nest
