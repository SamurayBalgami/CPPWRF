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
    const int ni0_base = (ci - ipos) * nri;
    const int nj0_base = (cj - jpos) * nrj;

    if (!xstag && !ystag) {
      const double inv_pts = 1.0 / static_cast<double>(nri * nrj);
      double sum = 0.0;
      for (int dj = 0; dj < nrj; ++dj) {
        for (int di = 0; di < nri; ++di) {
          sum += inv_pts * static_cast<double>(nfld(ni0_base + di, nj0_base + dj, k));
        }
      }
      cfld(ci0, cj0, k) = static_cast<real>(sum);
    } else if (xstag && !ystag) {
      const double inv_pts = 1.0 / static_cast<double>(nrj);
      double sum = 0.0;
      for (int dj = 0; dj < nrj; ++dj) {
        sum += inv_pts * static_cast<double>(nfld(ni0_base, nj0_base + dj, k));
      }
      cfld(ci0, cj0, k) = static_cast<real>(sum);
    } else if (!xstag && ystag) {
      const double inv_pts = 1.0 / static_cast<double>(nri);
      double sum = 0.0;
      for (int di = 0; di < nri; ++di) {
        sum += inv_pts * static_cast<double>(nfld(ni0_base + di, nj0_base, k));
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
    const int ni0_base = (ci - ipos) * nri;
    const int nj0_base = (cj - jpos) * nrj;
    const double inv_pts = 1.0 / static_cast<double>(nri * nrj);
    double sum = 0.0;
    for (int dj = 0; dj < nrj; ++dj) {
      for (int di = 0; di < nri; ++di) {
        sum += inv_pts * static_cast<double>(nfld(ni0_base + di, nj0_base + dj));
      }
    }
    cfld(ci0, cj0) = static_cast<real>(sum);
  }
};

// Point-sampling feedback functor matching WRF share/interp_fcn.F copy_fcnm
template <class CView, class NView>
struct KCopyFcnmView2D {
  CView cfld;
  NView nfld;
  int ipos, jpos, nri, nrj;
  bool xstag, ystag;

  CPPWRF_HD void operator()(int ci0, int cj0, int /*k*/) const {
    const int ci = ci0 + 1;
    const int cj = cj0 + 1;
    const int istag = xstag ? 0 : 1;
    const int jstag = ystag ? 0 : 1;
    if (nrj % 2 != 0) {
      const int ni0 = (ci - ipos) * nri + istag;
      const int nj0 = (cj - jpos) * nrj + jstag;
      cfld(ci0, cj0) = nfld(ni0, nj0);
    } else {
      const int ni0 = (ci - ipos) * nri + (nri / 2 - 1);
      const int nj0 = (cj - jpos) * nrj + (nrj / 2 - 1);
      cfld(ci0, cj0) = nfld(ni0, nj0);
    }
  }
};

// Helper to obtain mutable View3<real> / View2<real> and const View3<const real> / View2<const real>
// from either FieldT or non-uniform BaseFieldT
template <class Field3T>
inline View3<real> mut_view3(Field3T& f) {
  return View3<real>{f.data(), int(f.layout().sj), int(f.layout().sk)};
}
template <class Field3T>
inline View3<const real> const_view3(const Field3T& f) {
  return View3<const real>{f.data(), int(f.layout().sj), int(f.layout().sk)};
}
template <class Field2T>
inline View2<real> mut_view2(Field2T& f) {
  return View2<real>{f.data(), int(f.layout().sj)};
}
template <class Field2T>
inline View2<const real> const_view2(const Field2T& f) {
  return View2<const real>{f.data(), int(f.layout().sj)};
}

struct KCopy3D {
  View3<real> dst;
  View3<const real> src;
  CPPWRF_HD void operator()(int i, int j, int k) const { dst(i, j, k) = src(i, j, k); }
};

struct KSm121J3D {
  View3<real> tv;
  View3<const real> ccv;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    tv(i, j, k) = static_cast<real>(0.25 * (static_cast<double>(ccv(i, j + 1, k)) +
                                            2.0 * static_cast<double>(ccv(i, j, k)) +
                                            static_cast<double>(ccv(i, j - 1, k))));
  }
};

struct KSm121I3D {
  View3<real> cv;
  View3<const real> tcv;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    cv(i, j, k) = static_cast<real>(0.25 * (static_cast<double>(tcv(i + 1, j, k)) +
                                            2.0 * static_cast<double>(tcv(i, j, k)) +
                                            static_cast<double>(tcv(i - 1, j, k))));
  }
};

struct KSmdsmJ3D {
  View3<real> tv;
  View3<const real> ccv;
  double nu;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    double c0 = static_cast<double>(ccv(i, j, k));
    double cp = static_cast<double>(ccv(i, j + 1, k));
    double cm = static_cast<double>(ccv(i, j - 1, k));
    tv(i, j, k) = static_cast<real>(c0 + nu * ((cp + cm) * 0.5 - c0));
  }
};

struct KSmdsmI3D {
  View3<real> tv;
  View3<const real> ccv;
  double nu;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    double c0 = static_cast<double>(ccv(i, j, k));
    double cp = static_cast<double>(ccv(i + 1, j, k));
    double cm = static_cast<double>(ccv(i - 1, j, k));
    tv(i, j, k) = static_cast<real>(c0 + nu * ((cp + cm) * 0.5 - c0));
  }
};

struct KCopy2D {
  View2<real> dst;
  View2<const real> src;
  CPPWRF_HD void operator()(int i, int j, int /*k*/) const { dst(i, j) = src(i, j); }
};

struct KSm121J2D {
  View2<real> tv;
  View2<const real> ccv;
  CPPWRF_HD void operator()(int i, int j, int /*k*/) const {
    tv(i, j) = static_cast<real>(0.25 * (static_cast<double>(ccv(i, j + 1)) +
                                         2.0 * static_cast<double>(ccv(i, j)) +
                                         static_cast<double>(ccv(i, j - 1))));
  }
};

struct KSm121I2D {
  View2<real> cv;
  View2<const real> tcv;
  CPPWRF_HD void operator()(int i, int j, int /*k*/) const {
    cv(i, j) = static_cast<real>(0.25 * (static_cast<double>(tcv(i + 1, j)) +
                                         2.0 * static_cast<double>(tcv(i, j)) +
                                         static_cast<double>(tcv(i - 1, j))));
  }
};

struct KSmdsmJ2D {
  View2<real> tv;
  View2<const real> ccv;
  double nu;
  CPPWRF_HD void operator()(int i, int j, int /*k*/) const {
    double c0 = static_cast<double>(ccv(i, j));
    double cp = static_cast<double>(ccv(i, j + 1));
    double cm = static_cast<double>(ccv(i, j - 1));
    tv(i, j) = static_cast<real>(c0 + nu * ((cp + cm) * 0.5 - c0));
  }
};

struct KSmdsmI2D {
  View2<real> tv;
  View2<const real> ccv;
  double nu;
  CPPWRF_HD void operator()(int i, int j, int /*k*/) const {
    double c0 = static_cast<double>(ccv(i, j));
    double cp = static_cast<double>(ccv(i + 1, j));
    double cm = static_cast<double>(ccv(i - 1, j));
    tv(i, j) = static_cast<real>(c0 + nu * ((cp + cm) * 0.5 - c0));
  }
};

// Single-source View-based feedback + optional smoother for Field3 / BaseField3
template <class Field3T>
inline void feedback_field3(
    const Grid& c_grid,
    Field3T& cfld,
    const Field3T& nfld,
    int n_nx, int n_ny, int nk,
    int ipos, int jpos,
    int nri, int nrj,
    bool xstag, bool ystag,
    int smooth_option = 0,
    int spec_zone = 1) {

  const int ci0_begin = (ipos - 1) + spec_zone;
  const int ci0_end   = (ipos - 1) + n_nx / nri - spec_zone + (xstag ? 1 : 0);
  const int cj0_begin = (jpos - 1) + spec_zone;
  const int cj0_end   = (jpos - 1) + n_ny / nrj - spec_zone + (ystag ? 1 : 0);

  auto cv = mut_view3(cfld);
  auto nv = const_view3(nfld);
  parallel_for(Range3D{ci0_begin, ci0_end, cj0_begin, cj0_end, 0, nk},
               KCopyFcnView3D<decltype(cv), decltype(nv)>{cv, nv, ipos, jpos, nri, nrj, xstag, ystag});

  if (smooth_option == 1 || smooth_option == 2) {
    const int s_i0 = ipos + 1;
    const int s_i1 = (ipos - 1) + n_nx / nri - 2 + (xstag ? 1 : 0);
    const int s_j0 = jpos + 1;
    const int s_j1 = (jpos - 1) + n_ny / nrj - 2 + (ystag ? 1 : 0);
    if (s_i0 < s_i1 && s_j0 < s_j1) {
      Stagger stag = xstag ? Stagger::U : (ystag ? Stagger::V : ((nk == c_grid.nz + 1) ? Stagger::W : Stagger::M));
      Field3 tmp(c_grid, stag, "fb_tmp3");
      auto tv = tmp.view();
      auto tcv = tmp.cview();
      auto ccv = const_view3(cfld);
      if (smooth_option == 1) {
        parallel_for(Range3D{s_i0 - 1, s_i1 + 1, s_j0 - 1, s_j1 + 1, 0, nk}, KCopy3D{tv, ccv});
        parallel_for(Range3D{s_i0, s_i1, s_j0, s_j1, 0, nk}, KSm121J3D{tv, ccv});
        parallel_for(Range3D{s_i0, s_i1, s_j0, s_j1, 0, nk}, KSm121I3D{cv, tcv});
      } else {
        const double xnu[2] = {0.50, -0.52};
        for (int pass = 0; pass < 2; ++pass) {
          const double nu = xnu[pass];
          parallel_for(Range3D{s_i0, s_i1, s_j0, s_j1, 0, nk}, KSmdsmJ3D{tv, ccv, nu});
          parallel_for(Range3D{s_i0, s_i1, s_j0, s_j1, 0, nk}, KCopy3D{cv, tcv});
          parallel_for(Range3D{s_i0, s_i1, s_j0, s_j1, 0, nk}, KSmdsmI3D{tv, ccv, nu});
          parallel_for(Range3D{s_i0, s_i1, s_j0, s_j1, 0, nk}, KCopy3D{cv, tcv});
        }
      }
    }
  }
}

// Single-source View-based feedback + optional smoother for Field2 / BaseField2
template <class Field2T>
inline void feedback_field2(
    const Grid& c_grid,
    Field2T& cfld,
    const Field2T& nfld,
    int n_nx, int n_ny,
    int ipos, int jpos,
    int nri, int nrj,
    int smooth_option = 0,
    int spec_zone = 1) {

  const int ci0_begin = (ipos - 1) + spec_zone;
  const int ci0_end   = (ipos - 1) + n_nx / nri - spec_zone;
  const int cj0_begin = (jpos - 1) + spec_zone;
  const int cj0_end   = (jpos - 1) + n_ny / nrj - spec_zone;

  auto cv = mut_view2(cfld);
  auto nv = const_view2(nfld);
  parallel_for(Range3D{ci0_begin, ci0_end, cj0_begin, cj0_end, 0, 1},
               KCopyFcnView2D<decltype(cv), decltype(nv)>{cv, nv, ipos, jpos, nri, nrj});

  if (smooth_option == 1 || smooth_option == 2) {
    const int s_i0 = ipos + 1;
    const int s_i1 = (ipos - 1) + n_nx / nri - 2;
    const int s_j0 = jpos + 1;
    const int s_j1 = (jpos - 1) + n_ny / nrj - 2;
    if (s_i0 < s_i1 && s_j0 < s_j1) {
      Field2 tmp(c_grid, Stagger::M, "fb_tmp2");
      auto tv = tmp.view();
      auto tcv = tmp.cview();
      auto ccv = const_view2(cfld);
      if (smooth_option == 1) {
        parallel_for(Range3D{s_i0 - 1, s_i1 + 1, s_j0 - 1, s_j1 + 1, 0, 1}, KCopy2D{tv, ccv});
        parallel_for(Range3D{s_i0, s_i1, s_j0, s_j1, 0, 1}, KSm121J2D{tv, ccv});
        parallel_for(Range3D{s_i0, s_i1, s_j0, s_j1, 0, 1}, KSm121I2D{cv, tcv});
      } else {
        const double xnu[2] = {0.50, -0.52};
        for (int pass = 0; pass < 2; ++pass) {
          const double nu = xnu[pass];
          parallel_for(Range3D{s_i0, s_i1, s_j0, s_j1, 0, 1}, KSmdsmJ2D{tv, ccv, nu});
          parallel_for(Range3D{s_i0, s_i1, s_j0, s_j1, 0, 1}, KCopy2D{cv, tcv});
          parallel_for(Range3D{s_i0, s_i1, s_j0, s_j1, 0, 1}, KSmdsmI2D{tv, ccv, nu});
          parallel_for(Range3D{s_i0, s_i1, s_j0, s_j1, 0, 1}, KCopy2D{cv, tcv});
        }
      }
    }
  }
}

// Single-source View-based point-sampling feedback (copy_fcnm) for 2D map-factor / Coriolis fields
template <class Field2T>
inline void feedback_field2_fcnm(
    Field2T& cfld,
    const Field2T& nfld,
    int n_nx, int n_ny,
    int ipos, int jpos,
    int nri, int nrj,
    bool xstag = false,
    bool ystag = false,
    int spec_zone = 1) {

  const int ci0_begin = (ipos - 1) + spec_zone;
  const int ci0_end   = (ipos - 1) + n_nx / nri - spec_zone + (xstag ? 1 : 0);
  const int cj0_begin = (jpos - 1) + spec_zone;
  const int cj0_end   = (jpos - 1) + n_ny / nrj - spec_zone + (ystag ? 1 : 0);

  auto cv = mut_view2(cfld);
  auto nv = const_view2(nfld);
  parallel_for(Range3D{ci0_begin, ci0_end, cj0_begin, cj0_end, 0, 1},
               KCopyFcnmView2D<decltype(cv), decltype(nv)>{cv, nv, ipos, jpos, nri, nrj, xstag, ystag});
}

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
