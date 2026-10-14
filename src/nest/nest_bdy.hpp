#pragma once

// ============================================================================
// CPPWRF Lateral Boundary Interpolation Kernels (bdy_interp / bdy_interp1)
// Ports WRF 4.6.1 share/interp_fcn.F (bdy_interp1)
// Computes boundary values and tendencies:
//   value = current nest field
//   tendency = rdt * (sint(parent new) - nest), rdt = 1.D0 / cdt
// In our LateralBC slab layout (pos, k, b), pos fastest (src/dyn/lateral_bc.hpp)
// ============================================================================

#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdint>
#include "core/backend.hpp"
#include "nest/nest_types.hpp"
#include "nest/nest_interp.hpp"
#include "dyn/lateral_bc.hpp"

namespace cppwrf::nest {

// Structure holding the 4 boundary slabs for values and tendencies in LateralBC layout
struct LateralBCSlabs {
  // Slabs: XS (West), XE (East), YS (South), YE (North)
  // Each slab has size npos * nk * width
  std::vector<double> val_xs, tend_xs;
  std::vector<double> val_xe, tend_xe;
  std::vector<double> val_ys, tend_ys;
  std::vector<double> val_ye, tend_ye;
};

// ============================================================================
// Routine bdy_interp1 matching WRF share/interp_fcn.F
// Fills separate Fortran-dimensioned slabs:
//   bdy_xs,  bdy_txs  : (njms:njme, nkms:nkme, spec_bdy_width)
//   bdy_xe,  bdy_txe  : (njms:njme, nkms:nkme, spec_bdy_width)
//   bdy_ys,  bdy_tys  : (nims:nime, nkms:nkme, spec_bdy_width)
//   bdy_ye,  bdy_tye  : (nims:nime, nkms:nkme, spec_bdy_width)
// ============================================================================
inline void bdy_interp1(
    const double* cfld,
    int cids, int cide, int ckds, int ckde, int cjds, int cjde,
    int cims, int cime, int ckms, int ckme, int cjms, int cjme,
    int cits, int cite, int ckts, int ckte, int cjts, int cjte,
    const double* nfld,
    int spec_bdy_width,
    int spec_zone, int relax_zone,
    int nids, int nide, int nkds, int nkde, int njds, int njde,
    int nims, int nime, int nkms, int nkme, int njms, int njme,
    int nits, int nite, int nkts, int nkte, int njts, int njte,
    int /*shw1*/,
    const int* /*imask*/,
    bool xstag, bool ystag,
    int ipos, int jpos,
    int nri, int nrj,
    double* bdy_xs, double* bdy_txs,
    double* bdy_xe, double* bdy_txe,
    double* bdy_ys, double* bdy_tys,
    double* bdy_ye, double* bdy_tye,
    double cdt, double /*ndt*/) {

  const double rdt = 1.0 / cdt;
  int ioff = 0, joff = 0;
  if (xstag) ioff = std::max((nri - 1) / 2, 1);
  if (ystag) joff = std::max((nrj - 1) / 2, 1);

  const int sz = std::min(std::max(spec_zone, relax_zone + 1), spec_bdy_width);
  const int nfx = nri * nrj;

  const int c_ni = cime - cims + 1;
  const int c_nk = ckme - ckms + 1;
  const int c_nj = cjme - cjms + 1;

  auto c_idx = [&](int i, int k, int j) -> std::size_t {
    return static_cast<std::size_t>((i - cims) + (k - ckms) * c_ni + (j - cjms) * c_ni * c_nk);
  };

  const int n_ni = nime - nims + 1;
  const int n_nk = nkme - nkms + 1;
  const int n_nj = njme - njms + 1;

  auto n_idx = [&](int i, int k, int j) -> std::size_t {
    return static_cast<std::size_t>((i - nims) + (k - nkms) * n_ni + (j - njms) * n_ni * n_nk);
  };

  // Slab indexing in Fortran column-major layout: (pos, k, b) where pos is fastest, b 1-based (1..spec_bdy_width)
  auto xs_idx = [&](int nj, int k, int b) -> std::size_t {
    return static_cast<std::size_t>((nj - njms) + (k - nkms) * n_nj + (b - 1) * n_nj * n_nk);
  };
  auto ys_idx = [&](int ni, int k, int b) -> std::size_t {
    return static_cast<std::size_t>((ni - nims) + (k - nkms) * n_ni + (b - 1) * n_ni * n_nk);
  };

  std::vector<double> psca1(static_cast<std::size_t>(c_ni * c_nj * nfx));
  std::vector<double> psca(static_cast<std::size_t>(c_ni * c_nj * nfx));
  std::vector<int> dummy_mask(static_cast<std::size_t>(c_ni * c_nj), 0);

  auto psca_idx = [&](int i, int j, int nf_idx) -> std::size_t {
    return static_cast<std::size_t>((i - cims) + (j - cjms) * c_ni + (nf_idx - 1) * c_ni * c_nj);
  };

  auto n2ci = [&](int n) -> int { return (n + ipos * nri - 1) / nri; };
  auto n2cj = [&](int n) -> int { return (n + jpos * nrj - 1) / nrj; };

  for (int k = ckts; k <= ckte; ++k) {
    // Fill psca1 for level k
    for (int nf = 1; nf <= nfx; ++nf) {
      for (int j = cjms; j <= cjme; ++j) {
        for (int i = cims; i <= cime; ++i) {
          psca1[psca_idx(i, j, nf)] = cfld[c_idx(i, k, j)];
        }
      }
    }

    // Call sintb on the 4 boundaries
    // South BDY
    if (njts >= njds && njts <= njds + sz + joff) {
      sintb_fortran_layout(
          psca1.data(), psca.data(),
          cims, cime, cjms, cjme, dummy_mask.data(),
          n2ci(nits) - 1, n2ci(nite) + 1,
          n2cj(std::max(njts, njds)), n2cj(std::min(njte, njds + sz + joff)),
          nfx, xstag, ystag);
    }
    // North BDY
    if (njte <= njde && njte >= njde - sz - joff) {
      sintb_fortran_layout(
          psca1.data(), psca.data(),
          cims, cime, cjms, cjme, dummy_mask.data(),
          n2ci(nits) - 1, n2ci(nite) + 1,
          n2cj(std::max(njts, njde - sz - joff)), n2cj(std::min(njte, njde - 1 + joff)),
          nfx, xstag, ystag);
    }
    // West BDY
    if (nits >= nids && nits <= nids + sz + ioff) {
      sintb_fortran_layout(
          psca1.data(), psca.data(),
          cims, cime, cjms, cjme, dummy_mask.data(),
          n2ci(std::max(nits, nids)), n2ci(std::min(nite, nids + sz + ioff)),
          n2cj(njts) - 1, n2cj(njte) + 1,
          nfx, xstag, ystag);
    }
    // East BDY
    if (nite <= nide && nite >= nide - sz - ioff) {
      sintb_fortran_layout(
          psca1.data(), psca.data(),
          cims, cime, cjms, cjme, dummy_mask.data(),
          n2ci(std::max(nits, nide - sz - ioff)), n2ci(std::min(nite, nide - 1 + ioff)),
          n2cj(njts) - 1, n2cj(njte) + 1,
          nfx, xstag, ystag);
    }

    // Unpack interpolated boundary values and tendencies
    for (int nj1 = std::max(njds, njts - 1); nj1 <= std::min(njde + joff, njte + joff + 1); ++nj1) {
      int cj = jpos + (nj1 - 1) / nrj;
      int jp = (nj1 - 1) % nrj;
      for (int ni1 = std::max(nids, nits - 1); ni1 <= std::min(nide + ioff, nite + ioff + 1); ++ni1) {
        int ci = ipos + (ni1 - 1) / nri;
        int ip = (ni1 - 1) % nri;
        int ni = ni1 - ioff;
        int nj = nj1 - joff;

        if (ni < nids || ni > nime || nj < njds || nj > njme) continue;

        double psca_val = psca[psca_idx(ci, cj, ip + 1 + jp * nri)];
        double cur_val  = nfld[n_idx(ni, k, nj)];
        double tend_val = rdt * (psca_val - cur_val);

        // West
        if (ni >= nids && ni < nids + sz) {
          int b = ni - nids + 1;
          bdy_txs[xs_idx(nj, k, b)] = tend_val;
          bdy_xs[xs_idx(nj, k, b)]  = cur_val;
        }
        // South
        if (nj >= njds && nj < njds + sz) {
          int b = nj - njds + 1;
          bdy_tys[ys_idx(ni, k, b)] = tend_val;
          bdy_ys[ys_idx(ni, k, b)]  = cur_val;
        }
        // East
        if (xstag) {
          if (ni >= nide - sz + 1 && ni <= nide) {
            int b = nide - ni + 1;
            bdy_txe[xs_idx(nj, k, b)] = tend_val;
            bdy_xe[xs_idx(nj, k, b)]  = cur_val;
          }
        } else {
          if (ni >= nide - sz && ni <= nide - 1) {
            int b = nide - ni;
            bdy_txe[xs_idx(nj, k, b)] = tend_val;
            bdy_xe[xs_idx(nj, k, b)]  = cur_val;
          }
        }
        // North
        if (ystag) {
          if (nj >= njde - sz + 1 && nj <= njde) {
            int b = njde - nj + 1;
            bdy_tye[ys_idx(ni, k, b)] = tend_val;
            bdy_ye[ys_idx(ni, k, b)]  = cur_val;
          }
        } else {
          if (nj >= njde - sz && nj <= njde - 1) {
            int b = njde - nj;
            bdy_tye[ys_idx(ni, k, b)] = tend_val;
            bdy_ye[ys_idx(ni, k, b)]  = cur_val;
          }
        }
      }
    }
  }
}

// ============================================================================
// Direct LateralBC Slab generator
// Generates slabs directly in CPPWRF LateralBC layout:
//   (pos, k, b) with pos fastest, 0 <= b < width
// Slabs returned in LateralBCSlabs:
//   XS, XE: pos is j, npos = ny (or ny+1 for V)
//   YS, YE: pos is i, npos = nx (or nx+1 for U)
// ============================================================================
inline LateralBCSlabs compute_lateral_bc_slabs(
    const double* cfld,
    int c_nx, int c_ny, int c_nz,
    const double* nfld,
    int n_nx, int n_ny, int n_nz,
    int ipos, int jpos,
    int ratio_x, int ratio_y,
    double cdt,
    bool xstag, bool ystag,
    int width = 5, int spec_zone = 1, int relax_zone = 4) {

  LateralBCSlabs slabs;

  const int npos_x = xstag ? n_nx + 1 : n_nx;
  const int npos_y = ystag ? n_ny + 1 : n_ny;
  const int nk = n_nz;

  slabs.val_xs.assign(static_cast<std::size_t>(npos_y * nk * width), 0.0);
  slabs.tend_xs.assign(static_cast<std::size_t>(npos_y * nk * width), 0.0);
  slabs.val_xe.assign(static_cast<std::size_t>(npos_y * nk * width), 0.0);
  slabs.tend_xe.assign(static_cast<std::size_t>(npos_y * nk * width), 0.0);

  slabs.val_ys.assign(static_cast<std::size_t>(npos_x * nk * width), 0.0);
  slabs.tend_ys.assign(static_cast<std::size_t>(npos_x * nk * width), 0.0);
  slabs.val_ye.assign(static_cast<std::size_t>(npos_x * nk * width), 0.0);
  slabs.tend_ye.assign(static_cast<std::size_t>(npos_x * nk * width), 0.0);

  const int cids = 1, cide = c_nx + 1;
  const int ckds = 1, ckde = c_nz + 1;
  const int cjds = 1, cjde = c_ny + 1;

  const int cims = 1, cime = c_nx + (xstag ? 1 : 0);
  const int ckms = 1, ckme = c_nz;
  const int cjms = 1, cjme = c_ny + (ystag ? 1 : 0);

  const int nids = 1, nide = n_nx + 1;
  const int nkds = 1, nkde = n_nz + 1;
  const int njds = 1, njde = n_ny + 1;

  const int nims = 1, nime = npos_x;
  const int nkms = 1, nkme = nk;
  const int njms = 1, njme = npos_y;

  bdy_interp1(
      cfld,
      cids, cide, ckds, ckde, cjds, cjde,
      cims, cime, ckms, ckme, cjms, cjme,
      1, c_nx, 1, c_nz, 1, c_ny,
      nfld,
      width, spec_zone, relax_zone,
      nids, nide, nkds, nkde, njds, njde,
      nims, nime, nkms, nkme, njms, njme,
      1, npos_x, 1, nk, 1, npos_y,
      0, nullptr,
      xstag, ystag,
      ipos, jpos,
      ratio_x, ratio_y,
      slabs.val_xs.data(), slabs.tend_xs.data(),
      slabs.val_xe.data(), slabs.tend_xe.data(),
      slabs.val_ys.data(), slabs.tend_ys.data(),
      slabs.val_ye.data(), slabs.tend_ye.data(),
      cdt, cdt / ratio_x);

  return slabs;
}

// ============================================================================
// Single-Source Device Functor KBdySlab (raw Fortran-layout arrays)
// Direct parallel evaluation of LateralBC slabs on device/CPU (pos, k, b)
// ============================================================================
struct KBdySlab {
  const double* cfld;
  const double* nfld;
  real* val_slab;
  real* tend_slab;
  int c_nx, c_ny, c_nz;
  int n_nx, n_ny, n_nz;
  int ipos, jpos;
  int nri, nrj;
  double cdt;
  bool xstag, ystag;
  int width;
  int sz;
  int side;  // 0: XS, 1: XE, 2: YS, 3: YE
  int npos;
  int nk;

  CPPWRF_HD void operator()(int pos, int k, int b) const {
    if (b >= sz) return;
    int ioff = xstag ? ((nri - 1) / 2) : 0;
    int joff = ystag ? ((nrj - 1) / 2) : 0;
    if (xstag && ioff < 1) ioff = 1;
    if (ystag && joff < 1) joff = 1;

    int ni = 0, nj = 0;
    if (side == 0) {  // XS
      ni = b + 1;
      nj = pos + 1;
    } else if (side == 1) {  // XE
      ni = (xstag ? n_nx + 1 : n_nx) - b;
      nj = pos + 1;
    } else if (side == 2) {  // YS
      ni = pos + 1;
      nj = b + 1;
    } else {  // YE
      ni = pos + 1;
      nj = (ystag ? n_ny + 1 : n_ny) - b;
    }

    int ni_raw = ni + ioff;
    int nj_raw = nj + joff;
    int ci = ipos + (ni_raw - 1) / nri;
    int cj = jpos + (nj_raw - 1) / nrj;
    int ip = (ni_raw - 1) % nri;
    int jp = (nj_raw - 1) % nrj;

    double ax, ay;
    compute_sint_subgrid_offsets(nri, xstag, ystag, ip, jp, ax, ay);

    const int c_stride_x = 1;
    const int c_stride_k = (xstag ? c_nx + 1 : c_nx);
    const int c_stride_j = (xstag ? c_nx + 1 : c_nx) * c_nz;

    double patch[5][5];
    for (int dj = -2; dj <= 2; ++dj) {
      for (int di = -2; di <= 2; ++di) {
        int pci = ci + di;
        int pcj = cj + dj;
        int max_ci = xstag ? c_nx + 1 : c_nx;
        int max_cj = ystag ? c_ny + 1 : c_ny;
        if (pci < 1) pci = 1;
        if (pci > max_ci) pci = max_ci;
        if (pcj < 1) pcj = 1;
        if (pcj > max_cj) pcj = max_cj;
        std::size_t c_idx = static_cast<std::size_t>((pci - 1) * c_stride_x + k * c_stride_k + (pcj - 1) * c_stride_j);
        patch[dj + 2][di + 2] = cfld[c_idx];
      }
    }
    double interp_val = sint_point_2d(patch, ax, ay);

    int n_stride_k = xstag ? n_nx + 1 : n_nx;
    int n_stride_j = (xstag ? n_nx + 1 : n_nx) * nk;
    std::size_t n_idx = static_cast<std::size_t>((ni - 1) + k * n_stride_k + (nj - 1) * n_stride_j);
    double cur_val = nfld[n_idx];
    const double rdt = (cdt > 0.0) ? (1.0 / cdt) : 0.0;
    double tend_val = rdt * (interp_val - cur_val);

    std::size_t s_idx = static_cast<std::size_t>(pos + k * npos + b * npos * nk);
    val_slab[s_idx] = static_cast<real>(cur_val);
    tend_slab[s_idx] = static_cast<real>(tend_val);
  }
};

// ============================================================================
// Single-Source Device Functor KBdySlabView3D / KBdySlabView2D
// Evaluates LateralBC slabs directly on device/CPU from View3<const real> / View2<const real>
// ============================================================================
struct KBdySlabView3D {
  View3<const real> cfld;
  View3<const real> nfld;
  real* val_slab;
  real* tend_slab;
  int n_nx, n_ny;
  int ipos, jpos;
  int nri, nrj;
  double cdt;
  bool xstag, ystag;
  int sz;
  int side;  // 0: XS, 1: XE, 2: YS, 3: YE
  int npos;
  int nk;

  CPPWRF_HD void operator()(int pos, int k, int b) const {
    if (b >= sz) return;
    int ioff = xstag ? ((nri - 1) / 2) : 0;
    int joff = ystag ? ((nrj - 1) / 2) : 0;
    if (xstag && ioff < 1) ioff = 1;
    if (ystag && joff < 1) joff = 1;

    int ni = 0, nj = 0;  // 1-based fine indices
    if (side == 0) {     // XS
      ni = b + 1;
      nj = pos + 1;
    } else if (side == 1) {  // XE
      ni = (xstag ? n_nx + 1 : n_nx) - b;
      nj = pos + 1;
    } else if (side == 2) {  // YS
      ni = pos + 1;
      nj = b + 1;
    } else {  // YE
      ni = pos + 1;
      nj = (ystag ? n_ny + 1 : n_ny) - b;
    }

    int ni_raw = ni + ioff;
    int nj_raw = nj + joff;
    int ci0 = (ipos - 1) + (ni_raw - 1) / nri;
    int cj0 = (jpos - 1) + (nj_raw - 1) / nrj;
    int ip = (ni_raw - 1) % nri;
    int jp = (nj_raw - 1) % nrj;

    double ax, ay;
    compute_sint_subgrid_offsets(nri, xstag, ystag, ip, jp, ax, ay);

    double patch[5][5];
    for (int dj = -2; dj <= 2; ++dj) {
      for (int di = -2; di <= 2; ++di) {
        patch[dj + 2][di + 2] = static_cast<double>(cfld(ci0 + di, cj0 + dj, k));
      }
    }
    double interp_val = sint_point_2d(patch, ax, ay);
    double cur_val = static_cast<double>(nfld(ni - 1, nj - 1, k));
    const double rdt = (cdt > 0.0) ? (1.0 / cdt) : 0.0;
    double tend_val = rdt * (interp_val - cur_val);

    std::size_t s_idx = static_cast<std::size_t>(pos) + static_cast<std::size_t>(k) * npos +
                        static_cast<std::size_t>(b) * npos * nk;
    val_slab[s_idx] = static_cast<real>(cur_val);
    tend_slab[s_idx] = static_cast<real>(tend_val);
  }
};

struct KBdySlabView2D {
  View2<const real> cfld;
  View2<const real> nfld;
  real* val_slab;
  real* tend_slab;
  int n_nx, n_ny;
  int ipos, jpos;
  int nri, nrj;
  double cdt;
  int sz;
  int side;  // 0: XS, 1: XE, 2: YS, 3: YE
  int npos;

  CPPWRF_HD void operator()(int pos, int /*k*/, int b) const {
    if (b >= sz) return;
    int ni = 0, nj = 0;  // 1-based fine indices
    if (side == 0) {     // XS
      ni = b + 1;
      nj = pos + 1;
    } else if (side == 1) {  // XE
      ni = n_nx - b;
      nj = pos + 1;
    } else if (side == 2) {  // YS
      ni = pos + 1;
      nj = b + 1;
    } else {  // YE
      ni = pos + 1;
      nj = n_ny - b;
    }

    int ci0 = (ipos - 1) + (ni - 1) / nri;
    int cj0 = (jpos - 1) + (nj - 1) / nrj;
    int ip = (ni - 1) % nri;
    int jp = (nj - 1) % nrj;

    double ax, ay;
    compute_sint_subgrid_offsets(nri, false, false, ip, jp, ax, ay);

    double patch[5][5];
    for (int dj = -2; dj <= 2; ++dj) {
      for (int di = -2; di <= 2; ++di) {
        patch[dj + 2][di + 2] = static_cast<double>(cfld(ci0 + di, cj0 + dj));
      }
    }
    double interp_val = sint_point_2d(patch, ax, ay);
    double cur_val = static_cast<double>(nfld(ni - 1, nj - 1));
    const double rdt = (cdt > 0.0) ? (1.0 / cdt) : 0.0;
    double tend_val = rdt * (interp_val - cur_val);

    std::size_t s_idx = static_cast<std::size_t>(pos) + static_cast<std::size_t>(b) * npos;
    val_slab[s_idx] = static_cast<real>(cur_val);
    tend_slab[s_idx] = static_cast<real>(tend_val);
  }
};

// Launch single-source boundary forcing for a 3D Field3 variable into LateralBC
inline void force_nest_bdy_field3(
    dyn::LateralBC& lbc,
    dyn::LateralBC::Var v,
    const Field3& cfld,
    const Field3& nfld,
    int n_nx, int n_ny,
    int ipos, int jpos,
    int nri, int nrj,
    double cdt,
    int spec_zone = 1, int relax_zone = 4) {

  const bool xstag = (v == dyn::LateralBC::U);
  const bool ystag = (v == dyn::LateralBC::V);
  const int nk = lbc.slab_nk(v);
  const int width = lbc.width();
  const int sz = std::min(std::max(spec_zone, relax_zone + 1), width);

  for (int s = 0; s < 4; ++s) {
    const auto side = static_cast<dyn::LateralBC::Side>(s);
    const int npos = lbc.slab_npos(v, side);
    real* d_val = lbc.slab(v, side, false);
    real* d_tend = lbc.slab(v, side, true);
    parallel_for(Range3D{0, npos, 0, nk, 0, sz},
                 KBdySlabView3D{cfld.cview(), nfld.cview(), d_val, d_tend,
                                n_nx, n_ny, ipos, jpos, nri, nrj, cdt,
                                xstag, ystag, sz, s, npos, nk});
  }
  lbc.mark_set(v);
}

// Launch single-source boundary forcing for 2D MU into LateralBC
inline void force_nest_bdy_field2(
    dyn::LateralBC& lbc,
    dyn::LateralBC::Var v,
    const Field2& cfld,
    const Field2& nfld,
    int n_nx, int n_ny,
    int ipos, int jpos,
    int nri, int nrj,
    double cdt,
    int spec_zone = 1, int relax_zone = 4) {

  const int width = lbc.width();
  const int sz = std::min(std::max(spec_zone, relax_zone + 1), width);

  for (int s = 0; s < 4; ++s) {
    const auto side = static_cast<dyn::LateralBC::Side>(s);
    const int npos = lbc.slab_npos(v, side);
    real* d_val = lbc.slab(v, side, false);
    real* d_tend = lbc.slab(v, side, true);
    parallel_for(Range3D{0, npos, 0, 1, 0, sz},
                 KBdySlabView2D{cfld.cview(), nfld.cview(), d_val, d_tend,
                                n_nx, n_ny, ipos, jpos, nri, nrj, cdt,
                                sz, s, npos});
  }
  lbc.mark_set(v);
}

// Launch single-source boundary forcing for a 3D scalar into LateralBC
inline void force_nest_bdy_scalar_field3(
    dyn::LateralBC& lbc,
    int scalar_idx,
    const Field3& cfld,
    const Field3& nfld,
    int n_nx, int n_ny, int n_nz,
    int ipos, int jpos,
    int nri, int nrj,
    double cdt,
    int spec_zone = 1, int relax_zone = 4) {

  const int width = lbc.width();
  const int sz = std::min(std::max(spec_zone, relax_zone + 1), width);

  for (int s = 0; s < 4; ++s) {
    const auto side = static_cast<dyn::LateralBC::Side>(s);
    const int npos = (s < 2) ? n_ny : n_nx;
    real* d_val = lbc.scalar_slab(scalar_idx, side, false);
    real* d_tend = lbc.scalar_slab(scalar_idx, side, true);
    parallel_for(Range3D{0, npos, 0, n_nz, 0, sz},
                 KBdySlabView3D{cfld.cview(), nfld.cview(), d_val, d_tend,
                                n_nx, n_ny, ipos, jpos, nri, nrj, cdt,
                                false, false, sz, s, npos, n_nz});
  }
}

// ============================================================================
// Routine to force LateralBC from raw host double arrays (used in unit tests)
// ============================================================================
inline void force_nest_lateral_bc(
    dyn::LateralBC& lbc,
    dyn::LateralBC::Var v,
    const double* cfld,
    int c_nx, int c_ny, int c_nz,
    const double* nfld,
    int n_nx, int n_ny, int n_nz,
    int ipos, int jpos,
    int ratio_x, int ratio_y,
    double cdt,
    int spec_zone = 1, int relax_zone = 4) {

  const bool xstag = (v == dyn::LateralBC::U);
  const bool ystag = (v == dyn::LateralBC::V);
  const int nk = (v == dyn::LateralBC::PH || v == dyn::LateralBC::W) ? (n_nz + 1) : ((v == dyn::LateralBC::MU) ? 1 : n_nz);
  const int c_nk = (v == dyn::LateralBC::PH || v == dyn::LateralBC::W) ? (c_nz + 1) : ((v == dyn::LateralBC::MU) ? 1 : c_nz);
  const int width = lbc.width();

  LateralBCSlabs slabs = compute_lateral_bc_slabs(
      cfld, c_nx, c_ny, c_nk,
      nfld, n_nx, n_ny, nk,
      ipos, jpos,
      ratio_x, ratio_y,
      cdt, xstag, ystag,
      width, spec_zone, relax_zone);

  lbc.set(v, dyn::LateralBC::XS, slabs.val_xs, slabs.tend_xs);
  lbc.set(v, dyn::LateralBC::XE, slabs.val_xe, slabs.tend_xe);
  lbc.set(v, dyn::LateralBC::YS, slabs.val_ys, slabs.tend_ys);
  lbc.set(v, dyn::LateralBC::YE, slabs.val_ye, slabs.tend_ye);
}

inline void force_nest_lateral_bc_scalar(
    dyn::LateralBC& lbc,
    int scalar_idx,
    const double* cfld,
    int c_nx, int c_ny, int c_nz,
    const double* nfld,
    int n_nx, int n_ny, int n_nz,
    int ipos, int jpos,
    int ratio_x, int ratio_y,
    double cdt,
    int spec_zone = 1, int relax_zone = 4) {

  const int width = lbc.width();

  LateralBCSlabs slabs = compute_lateral_bc_slabs(
      cfld, c_nx, c_ny, c_nz,
      nfld, n_nx, n_ny, n_nz,
      ipos, jpos,
      ratio_x, ratio_y,
      cdt, false, false,
      width, spec_zone, relax_zone);

  lbc.set_scalar(scalar_idx, dyn::LateralBC::XS, slabs.val_xs, slabs.tend_xs);
  lbc.set_scalar(scalar_idx, dyn::LateralBC::XE, slabs.val_xe, slabs.tend_xe);
  lbc.set_scalar(scalar_idx, dyn::LateralBC::YS, slabs.val_ys, slabs.tend_ys);
  lbc.set_scalar(scalar_idx, dyn::LateralBC::YE, slabs.val_ye, slabs.tend_ye);
}

}  // namespace cppwrf::nest
