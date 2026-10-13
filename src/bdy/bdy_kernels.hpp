#pragma once
// Port of WRF dyn_em/module_bc_em.F and share/module_bc.F lateral boundary routines.
// Single-source C++20 (FP32/FP64) with parallel_for / column_for support.
// Staggering and indexing convention: C++ k = Fortran k - 1, halo 3, SoA layout.

#include <cmath>
#include <algorithm>
#include "core/backend.hpp"
#include "core/field.hpp"
#include "core/grid.hpp"
#include "core/precision.hpp"
#include "bdy/bdy_data.hpp"

namespace cppwrf::bdy {

// ---------------------------------------------------------------------------
// 1. Boundary relaxation profile computation (lbc_fcx_gcx)
// ---------------------------------------------------------------------------
template <class Real>
void calc_lbc_fcx_gcx(Real* fcx, Real* gcx, int spec_bdy_width, int spec_zone,
                      int relax_zone, Real dt, Real spec_exp = Real(0)) {
  for (int b = 0; b < spec_bdy_width; ++b) {
    fcx[b] = Real(0);
    gcx[b] = Real(0);
  }
  // b in C++ is 0-indexed (b=0 is spec_zone; b=1..relax_zone-1 are relaxation cells)
  for (int b = spec_zone; b < spec_zone + relax_zone; ++b) {
    Real factor = Real(spec_zone + relax_zone - 1 - b) / Real(relax_zone - 1);
    Real spongeweight = std::exp(-Real(b - spec_zone) * spec_exp);
    fcx[b] = (Real(0.1) / dt) * factor * spongeweight;
    gcx[b] = (Real(1.0) / (Real(50.0) * dt)) * factor * spongeweight;
  }
}

// ---------------------------------------------------------------------------
// 2. Core relaxation tendency calculation (relax_bdytend_core, 3D)
// ---------------------------------------------------------------------------
template <class Real>
void relax_bdytend_core_3d(
    View3<Real> field_tend,
    View3<const Real> field,
    BdySlabView<const Real> bxs,
    BdySlabView<const Real> bxe,
    BdySlabView<const Real> bys,
    BdySlabView<const Real> bye,
    BdySlabView<const Real> btxs,
    BdySlabView<const Real> btxe,
    BdySlabView<const Real> btys,
    BdySlabView<const Real> btye,
    const Real* fcx,
    const Real* gcx,
    Real dtbc,
    char var,
    int nx, int ny, int nz,
    int spec_zone = 1,
    int relax_zone = 4) {

  const int ibs = 0;
  const int ibe = (var == 'u') ? nx : nx - 1;
  const int jbs = 0;
  const int jbe = (var == 'v') ? ny : ny - 1;
  const int ktf = (var == 'h') ? nz : nz - 1;

  // South (YS) boundary
  for (int b_dist = spec_zone; b_dist < spec_zone + relax_zone - 1; ++b_dist) {
    const int j = jbs + b_dist;
    const int b_limit = b_dist;
    const Real fc = fcx[b_dist];
    const Real gc = gcx[b_dist];
    for (int k = 0; k <= ktf; ++k) {
      for (int i = b_limit + ibs; i <= ibe - b_limit; ++i) {
        const int im1 = std::max(i - 1, ibs);
        const int ip1 = std::min(i + 1, ibe);
        Real fls0 = (bys(i, k, b_dist) + dtbc * btys(i, k, b_dist)) - field(i, j, k);
        Real fls1 = (bys(im1, k, b_dist) + dtbc * btys(im1, k, b_dist)) - field(im1, j, k);
        Real fls2 = (bys(ip1, k, b_dist) + dtbc * btys(ip1, k, b_dist)) - field(ip1, j, k);
        Real fls3 = (bys(i, k, b_dist - 1) + dtbc * btys(i, k, b_dist - 1)) - field(i, j - 1, k);
        Real fls4 = (bys(i, k, b_dist + 1) + dtbc * btys(i, k, b_dist + 1)) - field(i, j + 1, k);
        field_tend(i, j, k) += fc * fls0 - gc * (fls1 + fls2 + fls3 + fls4 - Real(4) * fls0);
      }
    }
  }

  // North (YE) boundary
  for (int b_dist = spec_zone; b_dist < spec_zone + relax_zone - 1; ++b_dist) {
    const int j = jbe - b_dist;
    const int b_limit = b_dist;
    const Real fc = fcx[b_dist];
    const Real gc = gcx[b_dist];
    for (int k = 0; k <= ktf; ++k) {
      for (int i = b_limit + ibs; i <= ibe - b_limit; ++i) {
        const int im1 = std::max(i - 1, ibs);
        const int ip1 = std::min(i + 1, ibe);
        Real fls0 = (bye(i, k, b_dist) + dtbc * btye(i, k, b_dist)) - field(i, j, k);
        Real fls1 = (bye(im1, k, b_dist) + dtbc * btye(im1, k, b_dist)) - field(im1, j, k);
        Real fls2 = (bye(ip1, k, b_dist) + dtbc * btye(ip1, k, b_dist)) - field(ip1, j, k);
        Real fls3 = (bye(i, k, b_dist - 1) + dtbc * btye(i, k, b_dist - 1)) - field(i, j + 1, k);
        Real fls4 = (bye(i, k, b_dist + 1) + dtbc * btye(i, k, b_dist + 1)) - field(i, j - 1, k);
        field_tend(i, j, k) += fc * fls0 - gc * (fls1 + fls2 + fls3 + fls4 - Real(4) * fls0);
      }
    }
  }

  // West (XS) boundary
  for (int b_dist = spec_zone; b_dist < spec_zone + relax_zone - 1; ++b_dist) {
    const int i = ibs + b_dist;
    const Real fc = fcx[b_dist];
    const Real gc = gcx[b_dist];
    for (int k = 0; k <= ktf; ++k) {
      for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
        Real fls0 = (bxs(j, k, b_dist) + dtbc * btxs(j, k, b_dist)) - field(i, j, k);
        Real fls1 = (bxs(j - 1, k, b_dist) + dtbc * btxs(j - 1, k, b_dist)) - field(i, j - 1, k);
        Real fls2 = (bxs(j + 1, k, b_dist) + dtbc * btxs(j + 1, k, b_dist)) - field(i, j + 1, k);
        Real fls3 = (bxs(j, k, b_dist - 1) + dtbc * btxs(j, k, b_dist - 1)) - field(i - 1, j, k);
        Real fls4 = (bxs(j, k, b_dist + 1) + dtbc * btxs(j, k, b_dist + 1)) - field(i + 1, j, k);
        field_tend(i, j, k) += fc * fls0 - gc * (fls1 + fls2 + fls3 + fls4 - Real(4) * fls0);
      }
    }
  }

  // East (XE) boundary
  for (int b_dist = spec_zone; b_dist < spec_zone + relax_zone - 1; ++b_dist) {
    const int i = ibe - b_dist;
    const Real fc = fcx[b_dist];
    const Real gc = gcx[b_dist];
    for (int k = 0; k <= ktf; ++k) {
      for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
        Real fls0 = (bxe(j, k, b_dist) + dtbc * btxe(j, k, b_dist)) - field(i, j, k);
        Real fls1 = (bxe(j - 1, k, b_dist) + dtbc * btxe(j - 1, k, b_dist)) - field(i, j - 1, k);
        Real fls2 = (bxe(j + 1, k, b_dist) + dtbc * btxe(j + 1, k, b_dist)) - field(i, j + 1, k);
        Real fls3 = (bxe(j, k, b_dist - 1) + dtbc * btxe(j, k, b_dist - 1)) - field(i + 1, j, k);
        Real fls4 = (bxe(j, k, b_dist + 1) + dtbc * btxe(j, k, b_dist + 1)) - field(i - 1, j, k);
        field_tend(i, j, k) += fc * fls0 - gc * (fls1 + fls2 + fls3 + fls4 - Real(4) * fls0);
      }
    }
  }
}

// 2D overload for mu
template <class Real>
void relax_bdytend_core_2d(
    View2<Real> field_tend,
    View2<const Real> field,
    BdySlabView<const Real> bxs,
    BdySlabView<const Real> bxe,
    BdySlabView<const Real> bys,
    BdySlabView<const Real> bye,
    BdySlabView<const Real> btxs,
    BdySlabView<const Real> btxe,
    BdySlabView<const Real> btys,
    BdySlabView<const Real> btye,
    const Real* fcx,
    const Real* gcx,
    Real dtbc,
    int nx, int ny,
    int spec_zone = 1,
    int relax_zone = 4) {

  const int ibs = 0, ibe = nx - 1;
  const int jbs = 0, jbe = ny - 1;
  const int k = 0;

  // South
  for (int b_dist = spec_zone; b_dist < spec_zone + relax_zone - 1; ++b_dist) {
    const int j = jbs + b_dist;
    const int b_limit = b_dist;
    const Real fc = fcx[b_dist];
    const Real gc = gcx[b_dist];
    for (int i = b_limit + ibs; i <= ibe - b_limit; ++i) {
      const int im1 = std::max(i - 1, ibs);
      const int ip1 = std::min(i + 1, ibe);
      Real fls0 = (bys(i, k, b_dist) + dtbc * btys(i, k, b_dist)) - field(i, j);
      Real fls1 = (bys(im1, k, b_dist) + dtbc * btys(im1, k, b_dist)) - field(im1, j);
      Real fls2 = (bys(ip1, k, b_dist) + dtbc * btys(ip1, k, b_dist)) - field(ip1, j);
      Real fls3 = (bys(i, k, b_dist - 1) + dtbc * btys(i, k, b_dist - 1)) - field(i, j - 1);
      Real fls4 = (bys(i, k, b_dist + 1) + dtbc * btys(i, k, b_dist + 1)) - field(i, j + 1);
      field_tend(i, j) += fc * fls0 - gc * (fls1 + fls2 + fls3 + fls4 - Real(4) * fls0);
    }
  }

  // North
  for (int b_dist = spec_zone; b_dist < spec_zone + relax_zone - 1; ++b_dist) {
    const int j = jbe - b_dist;
    const int b_limit = b_dist;
    const Real fc = fcx[b_dist];
    const Real gc = gcx[b_dist];
    for (int i = b_limit + ibs; i <= ibe - b_limit; ++i) {
      const int im1 = std::max(i - 1, ibs);
      const int ip1 = std::min(i + 1, ibe);
      Real fls0 = (bye(i, k, b_dist) + dtbc * btye(i, k, b_dist)) - field(i, j);
      Real fls1 = (bye(im1, k, b_dist) + dtbc * btye(im1, k, b_dist)) - field(im1, j);
      Real fls2 = (bye(ip1, k, b_dist) + dtbc * btye(ip1, k, b_dist)) - field(ip1, j);
      Real fls3 = (bye(i, k, b_dist - 1) + dtbc * btye(i, k, b_dist - 1)) - field(i, j + 1);
      Real fls4 = (bye(i, k, b_dist + 1) + dtbc * btye(i, k, b_dist + 1)) - field(i, j - 1);
      field_tend(i, j) += fc * fls0 - gc * (fls1 + fls2 + fls3 + fls4 - Real(4) * fls0);
    }
  }

  // West
  for (int b_dist = spec_zone; b_dist < spec_zone + relax_zone - 1; ++b_dist) {
    const int i = ibs + b_dist;
    const Real fc = fcx[b_dist];
    const Real gc = gcx[b_dist];
    for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
      Real fls0 = (bxs(j, k, b_dist) + dtbc * btxs(j, k, b_dist)) - field(i, j);
      Real fls1 = (bxs(j - 1, k, b_dist) + dtbc * btxs(j - 1, k, b_dist)) - field(i, j - 1);
      Real fls2 = (bxs(j + 1, k, b_dist) + dtbc * btxs(j + 1, k, b_dist)) - field(i, j + 1);
      Real fls3 = (bxs(j, k, b_dist - 1) + dtbc * btxs(j, k, b_dist - 1)) - field(i - 1, j);
      Real fls4 = (bxs(j, k, b_dist + 1) + dtbc * btxs(j, k, b_dist + 1)) - field(i + 1, j);
      field_tend(i, j) += fc * fls0 - gc * (fls1 + fls2 + fls3 + fls4 - Real(4) * fls0);
    }
  }

  // East
  for (int b_dist = spec_zone; b_dist < spec_zone + relax_zone - 1; ++b_dist) {
    const int i = ibe - b_dist;
    const Real fc = fcx[b_dist];
    const Real gc = gcx[b_dist];
    for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
      Real fls0 = (bxe(j, k, b_dist) + dtbc * btxe(j, k, b_dist)) - field(i, j);
      Real fls1 = (bxe(j - 1, k, b_dist) + dtbc * btxe(j - 1, k, b_dist)) - field(i, j - 1);
      Real fls2 = (bxe(j + 1, k, b_dist) + dtbc * btxe(j + 1, k, b_dist)) - field(i, j + 1);
      Real fls3 = (bxe(j, k, b_dist - 1) + dtbc * btxe(j, k, b_dist - 1)) - field(i + 1, j);
      Real fls4 = (bxe(j, k, b_dist + 1) + dtbc * btxe(j, k, b_dist + 1)) - field(i - 1, j);
      field_tend(i, j) += fc * fls0 - gc * (fls1 + fls2 + fls3 + fls4 - Real(4) * fls0);
    }
  }
}

// ---------------------------------------------------------------------------
// 3. Specified boundary tendencies injection (spec_bdytend)
// ---------------------------------------------------------------------------
template <class Real>
void spec_bdytend_3d(
    View3<Real> field_tend,
    BdySlabView<const Real> btxs,
    BdySlabView<const Real> btxe,
    BdySlabView<const Real> btys,
    BdySlabView<const Real> btye,
    char var,
    int nx, int ny, int nz,
    int spec_zone = 1) {

  const int ibs = 0;
  const int ibe = (var == 'u') ? nx : nx - 1;
  const int jbs = 0;
  const int jbe = (var == 'v') ? ny : ny - 1;
  const int ktf = (var == 'h') ? nz : nz - 1;

  // South (YS)
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int j = jbs + b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int i = b_dist + ibs; i <= ibe - b_dist; ++i) {
        field_tend(i, j, k) = btys(i, k, b_dist);
      }
    }
  }

  // North (YE)
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int j = jbe - b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int i = b_dist + ibs; i <= ibe - b_dist; ++i) {
        field_tend(i, j, k) = btye(i, k, b_dist);
      }
    }
  }

  // West (XS)
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int i = ibs + b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
        field_tend(i, j, k) = btxs(j, k, b_dist);
      }
    }
  }

  // East (XE)
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int i = ibe - b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
        field_tend(i, j, k) = btxe(j, k, b_dist);
      }
    }
  }
}

template <class Real>
void spec_bdytend_2d(
    View2<Real> field_tend,
    BdySlabView<const Real> btxs,
    BdySlabView<const Real> btxe,
    BdySlabView<const Real> btys,
    BdySlabView<const Real> btye,
    int nx, int ny,
    int spec_zone = 1) {

  const int ibs = 0, ibe = nx - 1;
  const int jbs = 0, jbe = ny - 1;
  const int k = 0;

  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int j = jbs + b_dist;
    for (int i = b_dist + ibs; i <= ibe - b_dist; ++i) {
      field_tend(i, j) = btys(i, k, b_dist);
    }
  }
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int j = jbe - b_dist;
    for (int i = b_dist + ibs; i <= ibe - b_dist; ++i) {
      field_tend(i, j) = btye(i, k, b_dist);
    }
  }
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int i = ibs + b_dist;
    for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
      field_tend(i, j) = btxs(j, k, b_dist);
    }
  }
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int i = ibe - b_dist;
    for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
      field_tend(i, j) = btxe(j, k, b_dist);
    }
  }
}

// ---------------------------------------------------------------------------
// 4. Mass weighting helper (mass_weight)
// ---------------------------------------------------------------------------
template <class Real>
void mass_weight(
    View3<Real> rfield,
    View3<const Real> field,
    View2<const Real> mut,
    const Real* c1,
    const Real* c2,
    int nx, int ny, int nz) {
  for (int k = 0; k < nz; ++k) {
    const Real c1k = c1[k];
    const Real c2k = c2[k];
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        rfield(i, j, k) = field(i, j, k) * (c1k * mut(i, j) + c2k);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// 5. Full dry boundary relaxation (relax_bdy_dry)
// ---------------------------------------------------------------------------
template <class Real>
void relax_bdy_dry(
    View3<Real> ru_tend,
    View3<Real> rv_tend,
    View3<Real> ph_tend,
    View3<Real> t_tend,
    View2<Real> mu_tend,
    View3<const Real> ru,
    View3<const Real> rv,
    View3<const Real> ph,
    View3<const Real> t,
    View2<const Real> mu,
    View2<const Real> mut,
    const Real* c1h, const Real* c2h,
    const Real* c1f, const Real* c2f,
    const BdyData<Real>& bdy,
    const Real* fcx, const Real* gcx,
    Real dtbc,
    int nx, int ny, int nz,
    int spec_zone = 1,
    int relax_zone = 4) {

  // ru
  relax_bdytend_core_3d<Real>(ru_tend, ru, bdy.u.bxs.cview(), bdy.u.bxe.cview(),
                             bdy.u.bys.cview(), bdy.u.bye.cview(),
                             bdy.u.btxs.cview(), bdy.u.btxe.cview(),
                             bdy.u.btys.cview(), bdy.u.btye.cview(),
                             fcx, gcx, dtbc, 'u', nx, ny, nz, spec_zone, relax_zone);

  // rv
  relax_bdytend_core_3d<Real>(rv_tend, rv, bdy.v.bxs.cview(), bdy.v.bxe.cview(),
                             bdy.v.bys.cview(), bdy.v.bye.cview(),
                             bdy.v.btxs.cview(), bdy.v.btxe.cview(),
                             bdy.v.btys.cview(), bdy.v.btye.cview(),
                             fcx, gcx, dtbc, 'v', nx, ny, nz, spec_zone, relax_zone);

  // ph (mass-weighted with c1f, c2f, nz+1 levels)
  std::vector<Real> rfield_ph_buf(nx * ny * (nz + 1), Real(0));
  View3<Real> rfield_ph{rfield_ph_buf.data(), nx, nx * ny};
  mass_weight(rfield_ph, ph, mut, c1f, c2f, nx, ny, nz + 1);
  relax_bdytend_core_3d<Real>(ph_tend, rfield_ph, bdy.ph.bxs.cview(), bdy.ph.bxe.cview(),
                             bdy.ph.bys.cview(), bdy.ph.bye.cview(),
                             bdy.ph.btxs.cview(), bdy.ph.btxe.cview(),
                             bdy.ph.btys.cview(), bdy.ph.btye.cview(),
                             fcx, gcx, dtbc, 'h', nx, ny, nz, spec_zone, relax_zone);

  // t (mass-weighted with c1h, c2h, nz levels)
  std::vector<Real> rfield_t_buf(nx * ny * nz, Real(0));
  View3<Real> rfield_t{rfield_t_buf.data(), nx, nx * ny};
  mass_weight(rfield_t, t, mut, c1h, c2h, nx, ny, nz);
  relax_bdytend_core_3d<Real>(t_tend, rfield_t, bdy.t.bxs.cview(), bdy.t.bxe.cview(),
                             bdy.t.bys.cview(), bdy.t.bye.cview(),
                             bdy.t.btxs.cview(), bdy.t.btxe.cview(),
                             bdy.t.btys.cview(), bdy.t.btye.cview(),
                             fcx, gcx, dtbc, 't', nx, ny, nz, spec_zone, relax_zone);

  // mu (2D)
  relax_bdytend_core_2d<Real>(mu_tend, mu, bdy.mu.bxs.cview(), bdy.mu.bxe.cview(),
                             bdy.mu.bys.cview(), bdy.mu.bye.cview(),
                             bdy.mu.btxs.cview(), bdy.mu.btxe.cview(),
                             bdy.mu.btys.cview(), bdy.mu.btye.cview(),
                             fcx, gcx, dtbc, nx, ny, spec_zone, relax_zone);
}

// ---------------------------------------------------------------------------
// 6. Full dry boundary specified tendencies (spec_bdy_dry)
// ---------------------------------------------------------------------------
template <class Real>
void spec_bdy_dry(
    View3<Real> ru_tend,
    View3<Real> rv_tend,
    View3<Real> ph_tend,
    View3<Real> t_tend,
    View2<Real> mu_tend,
    const BdyData<Real>& bdy,
    int nx, int ny, int nz,
    int spec_zone = 1) {

  spec_bdytend_3d<Real>(ru_tend, bdy.u.btxs.cview(), bdy.u.btxe.cview(),
                        bdy.u.btys.cview(), bdy.u.btye.cview(), 'u', nx, ny, nz, spec_zone);
  spec_bdytend_3d<Real>(rv_tend, bdy.v.btxs.cview(), bdy.v.btxe.cview(),
                        bdy.v.btys.cview(), bdy.v.btye.cview(), 'v', nx, ny, nz, spec_zone);
  spec_bdytend_3d<Real>(ph_tend, bdy.ph.btxs.cview(), bdy.ph.btxe.cview(),
                        bdy.ph.btys.cview(), bdy.ph.btye.cview(), 'h', nx, ny, nz, spec_zone);
  spec_bdytend_3d<Real>(t_tend, bdy.t.btxs.cview(), bdy.t.btxe.cview(),
                        bdy.t.btys.cview(), bdy.t.btye.cview(), 't', nx, ny, nz, spec_zone);
  spec_bdytend_2d<Real>(mu_tend, bdy.mu.btxs.cview(), bdy.mu.btxe.cview(),
                        bdy.mu.btys.cview(), bdy.mu.btye.cview(), nx, ny, spec_zone);
}

// ---------------------------------------------------------------------------
// 7. Scalar boundary relaxation & specified tendencies
// ---------------------------------------------------------------------------
template <class Real>
void relax_bdy_scalar(
    View3<Real> scalar_tend,
    View3<const Real> scalar,
    View2<const Real> mu,
    const Real* c1h, const Real* c2h,
    const BdyVar<Real>& bdy_scalar,
    const Real* fcx, const Real* gcx,
    Real dtbc,
    int nx, int ny, int nz,
    int spec_zone = 1,
    int relax_zone = 4) {

  std::vector<Real> rscalar_buf(nx * ny * nz, Real(0));
  View3<Real> rscalar{rscalar_buf.data(), nx, nx * ny};
  mass_weight(rscalar, scalar, mu, c1h, c2h, nx, ny, nz);
  relax_bdytend_core_3d<Real>(scalar_tend, rscalar,
                             bdy_scalar.bxs.cview(), bdy_scalar.bxe.cview(),
                             bdy_scalar.bys.cview(), bdy_scalar.bye.cview(),
                             bdy_scalar.btxs.cview(), bdy_scalar.btxe.cview(),
                             bdy_scalar.btys.cview(), bdy_scalar.btye.cview(),
                             fcx, gcx, dtbc, 'q', nx, ny, nz, spec_zone, relax_zone);
}

template <class Real>
void spec_bdy_scalar(
    View3<Real> scalar_tend,
    const BdyVar<Real>& bdy_scalar,
    int nx, int ny, int nz,
    int spec_zone = 1) {
  spec_bdytend_3d<Real>(scalar_tend,
                        bdy_scalar.btxs.cview(), bdy_scalar.btxe.cview(),
                        bdy_scalar.btys.cview(), bdy_scalar.btye.cview(),
                        'q', nx, ny, nz, spec_zone);
}

// ---------------------------------------------------------------------------
// 8. Acoustic substep boundary state update (spec_bdyupdate)
// ---------------------------------------------------------------------------
template <class Real>
void spec_bdyupdate_3d(
    View3<Real> field,
    View3<const Real> field_tend,
    Real dt,
    char var,
    int nx, int ny, int nz,
    int spec_zone = 1) {

  const int ibs = 0;
  const int ibe = (var == 'u') ? nx : nx - 1;
  const int jbs = 0;
  const int jbe = (var == 'v') ? ny : ny - 1;
  const int ktf = (var == 'h') ? nz : nz - 1;

  // South
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int j = jbs + b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int i = b_dist + ibs; i <= ibe - b_dist; ++i) {
        field(i, j, k) += dt * field_tend(i, j, k);
      }
    }
  }

  // North
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int j = jbe - b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int i = b_dist + ibs; i <= ibe - b_dist; ++i) {
        field(i, j, k) += dt * field_tend(i, j, k);
      }
    }
  }

  // West
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int i = ibs + b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
        field(i, j, k) += dt * field_tend(i, j, k);
      }
    }
  }

  // East
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int i = ibe - b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
        field(i, j, k) += dt * field_tend(i, j, k);
      }
    }
  }
}

template <class Real>
void spec_bdyupdate_2d(
    View2<Real> field,
    View2<const Real> field_tend,
    Real dt,
    int nx, int ny,
    int spec_zone = 1) {

  const int ibs = 0, ibe = nx - 1;
  const int jbs = 0, jbe = ny - 1;

  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int j = jbs + b_dist;
    for (int i = b_dist + ibs; i <= ibe - b_dist; ++i) {
      field(i, j) += dt * field_tend(i, j);
    }
  }
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int j = jbe - b_dist;
    for (int i = b_dist + ibs; i <= ibe - b_dist; ++i) {
      field(i, j) += dt * field_tend(i, j);
    }
  }
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int i = ibs + b_dist;
    for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
      field(i, j) += dt * field_tend(i, j);
    }
  }
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int i = ibe - b_dist;
    for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
      field(i, j) += dt * field_tend(i, j);
    }
  }
}

// ---------------------------------------------------------------------------
// 9. Acoustic geopotential boundary update (spec_bdyupdate_ph)
// ---------------------------------------------------------------------------
template <class Real>
void spec_bdyupdate_ph(
    View3<const Real> ph_save,
    View3<Real> ph,
    View3<const Real> ph_tend,
    View2<const Real> mu_tend,
    View2<const Real> muts,
    const Real* c1f,
    const Real* c2f,
    Real dt,
    int nx, int ny, int nz,
    int spec_zone = 1) {

  const int ibs = 0, ibe = nx - 1;
  const int jbs = 0, jbe = ny - 1;
  const int ktf = nz; // geopotential has nz+1 levels (0..nz)

  auto update_cell = [&](int i, int j) {
    Real mu_old = muts(i, j) - dt * mu_tend(i, j);
    for (int k = 0; k <= ktf; ++k) {
      Real m_old = c1f[k] * mu_old + c2f[k];
      Real m_new = c1f[k] * muts(i, j) + c2f[k];
      Real ratio = m_old / m_new;
      ph(i, j, k) = ph(i, j, k) * ratio + dt * ph_tend(i, j, k) / m_new +
                    ph_save(i, j, k) * (ratio - Real(1));
    }
  };

  // South
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int j = jbs + b_dist;
    for (int i = b_dist + ibs; i <= ibe - b_dist; ++i) update_cell(i, j);
  }
  // North
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int j = jbe - b_dist;
    for (int i = b_dist + ibs; i <= ibe - b_dist; ++i) update_cell(i, j);
  }
  // West
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int i = ibs + b_dist;
    for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) update_cell(i, j);
  }
  // East
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int i = ibe - b_dist;
    for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) update_cell(i, j);
  }
}

// ---------------------------------------------------------------------------
// 10. Vertical velocity zero-gradient condition (zero_grad_bdy)
// ---------------------------------------------------------------------------
template <class Real>
void zero_grad_bdy(
    View3<Real> field,
    int nx, int ny, int nz,
    int spec_zone = 1) {

  const int ibs = 0, ibe = nx - 1;
  const int jbs = 0, jbe = ny - 1;
  const int ktf = nz; // w has nz+1 levels (0..nz)

  // South
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int j = jbs + b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int i = b_dist + ibs; i <= ibe - b_dist; ++i) {
        int i_inner = std::min(std::max(i, ibs + spec_zone), ibe - spec_zone);
        field(i, j, k) = field(i_inner, jbs + spec_zone, k);
      }
    }
  }

  // North
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int j = jbe - b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int i = b_dist + ibs; i <= ibe - b_dist; ++i) {
        int i_inner = std::min(std::max(i, ibs + spec_zone), ibe - spec_zone);
        field(i, j, k) = field(i_inner, jbe - spec_zone, k);
      }
    }
  }

  // West
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int i = ibs + b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
        int j_inner = std::min(std::max(j, jbs + spec_zone), jbe - spec_zone);
        field(i, j, k) = field(ibs + spec_zone, j_inner, k);
      }
    }
  }

  // East
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int i = ibe - b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
        int j_inner = std::min(std::max(j, jbs + spec_zone), jbe - spec_zone);
        field(i, j, k) = field(ibe - spec_zone, j_inner, k);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// 11. Flow-dependent boundary condition (flow_dep_bdy)
// ---------------------------------------------------------------------------
template <class Real>
void flow_dep_bdy(
    View3<Real> field,
    View3<const Real> u,
    View3<const Real> v,
    int nx, int ny, int nz,
    int spec_zone = 1) {

  const int ibs = 0, ibe = nx - 1;
  const int jbs = 0, jbe = ny - 1;
  const int ktf = nz - 1;

  // South: outflow if v < 0
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int j = jbs + b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int i = b_dist + ibs; i <= ibe - b_dist; ++i) {
        int i_inner = std::min(std::max(i, ibs + spec_zone), ibe - spec_zone);
        if (v(i, j, k) < Real(0)) {
          field(i, j, k) = field(i_inner, jbs + spec_zone, k);
        } else {
          field(i, j, k) = Real(0);
        }
      }
    }
  }

  // North: outflow if v > 0 (v at j+1)
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int j = jbe - b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int i = b_dist + ibs; i <= ibe - b_dist; ++i) {
        int i_inner = std::min(std::max(i, ibs + spec_zone), ibe - spec_zone);
        if (v(i, j + 1, k) > Real(0)) {
          field(i, j, k) = field(i_inner, jbe - spec_zone, k);
        } else {
          field(i, j, k) = Real(0);
        }
      }
    }
  }

  // West: outflow if u < 0
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int i = ibs + b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
        int j_inner = std::min(std::max(j, jbs + spec_zone), jbe - spec_zone);
        if (u(i, j, k) < Real(0)) {
          field(i, j, k) = field(ibs + spec_zone, j_inner, k);
        } else {
          field(i, j, k) = Real(0);
        }
      }
    }
  }

  // East: outflow if u > 0 (u at i+1)
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int i = ibe - b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
        int j_inner = std::min(std::max(j, jbs + spec_zone), jbe - spec_zone);
        if (u(i + 1, j, k) > Real(0)) {
          field(i, j, k) = field(ibe - spec_zone, j_inner, k);
        } else {
          field(i, j, k) = Real(0);
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
// 12. Final step boundary overwriting (spec_bdy_final)
// ---------------------------------------------------------------------------
template <class Real>
void spec_bdy_final_3d(
    View3<Real> field,
    View2<const Real> mu,
    const Real* c1,
    const Real* c2,
    View2<const Real> msf,
    BdySlabView<const Real> bxs,
    BdySlabView<const Real> bxe,
    BdySlabView<const Real> bys,
    BdySlabView<const Real> bye,
    BdySlabView<const Real> btxs,
    BdySlabView<const Real> btxe,
    BdySlabView<const Real> btys,
    BdySlabView<const Real> btye,
    char var,
    Real dtbc,
    int nx, int ny, int nz,
    int spec_zone = 1) {

  const int ibs = 0;
  const int ibe = (var == 'u') ? nx : nx - 1;
  const int jbs = 0;
  const int jbe = (var == 'v') ? ny : ny - 1;
  const int ktf = (var == 'h' || var == 'w') ? nz : nz - 1;
  const bool msfcouple = (var == 'u' || var == 'v' || var == 'w');
  const bool mucouple = (var != 'm');

  auto eval = [&](Real bval, int i, int j, int k) {
    Real xmsf = msfcouple ? msf(i, j) : Real(1);
    Real xmu = mucouple ? (c1[k] * mu(i, j) + c2[k]) : Real(1);
    return xmsf * bval / xmu;
  };

  // South
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int j = jbs + b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int i = b_dist + ibs; i <= ibe - b_dist; ++i) {
        Real bfield = bys(i, k, b_dist) + dtbc * btys(i, k, b_dist);
        field(i, j, k) = eval(bfield, i, j, k);
      }
    }
  }

  // North
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int j = jbe - b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int i = b_dist + ibs; i <= ibe - b_dist; ++i) {
        Real bfield = bye(i, k, b_dist) + dtbc * btye(i, k, b_dist);
        field(i, j, k) = eval(bfield, i, j, k);
      }
    }
  }

  // West
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int i = ibs + b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
        Real bfield = bxs(j, k, b_dist) + dtbc * btxs(j, k, b_dist);
        field(i, j, k) = eval(bfield, i, j, k);
      }
    }
  }

  // East
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int i = ibe - b_dist;
    for (int k = 0; k <= ktf; ++k) {
      for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
        Real bfield = bxe(j, k, b_dist) + dtbc * btxe(j, k, b_dist);
        field(i, j, k) = eval(bfield, i, j, k);
      }
    }
  }
}

template <class Real>
void spec_bdy_final_2d(
    View2<Real> field,
    BdySlabView<const Real> bxs,
    BdySlabView<const Real> bxe,
    BdySlabView<const Real> bys,
    BdySlabView<const Real> bye,
    BdySlabView<const Real> btxs,
    BdySlabView<const Real> btxe,
    BdySlabView<const Real> btys,
    BdySlabView<const Real> btye,
    Real dtbc,
    int nx, int ny,
    int spec_zone = 1) {

  const int ibs = 0, ibe = nx - 1;
  const int jbs = 0, jbe = ny - 1;
  const int k = 0;

  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int j = jbs + b_dist;
    for (int i = b_dist + ibs; i <= ibe - b_dist; ++i) {
      field(i, j) = bys(i, k, b_dist) + dtbc * btys(i, k, b_dist);
    }
  }
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int j = jbe - b_dist;
    for (int i = b_dist + ibs; i <= ibe - b_dist; ++i) {
      field(i, j) = bye(i, k, b_dist) + dtbc * btye(i, k, b_dist);
    }
  }
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int i = ibs + b_dist;
    for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
      field(i, j) = bxs(j, k, b_dist) + dtbc * btxs(j, k, b_dist);
    }
  }
  for (int b_dist = 0; b_dist < spec_zone; ++b_dist) {
    const int i = ibe - b_dist;
    for (int j = b_dist + jbs + 1; j <= jbe - b_dist - 1; ++j) {
      field(i, j) = bxe(j, k, b_dist) + dtbc * btxe(j, k, b_dist);
    }
  }
}

// ---------------------------------------------------------------------------
// 13. Halo flat extrapolation (set_physical_bc3d / set_physical_bc2d)
// ---------------------------------------------------------------------------
template <class Real>
void set_physical_bc3d_flat(
    View3<Real> field,
    char var,
    int nx, int ny, int nz,
    int halo = 3) {

  const bool is_u = (var == 'u' || var == 'U' || var == 'x' || var == 'X');
  const bool is_v = (var == 'v' || var == 'V' || var == 'y' || var == 'Y');
  const int ni = is_u ? (nx + 1) : nx;
  const int nj = is_v ? (ny + 1) : ny;
  const int ktf = (var == 'h' || var == 'H' || var == 'w' || var == 'W') ? nz : (nz - 1);

  const int ims = -halo, ime = ni - 1 + halo;
  const int jms = -halo, jme = nj - 1 + halo;

  for (int k = 0; k <= ktf; ++k) {
    // West & East halos
    for (int j = 0; j < nj; ++j) {
      for (int h = 1; h <= halo; ++h) {
        field(-h, j, k) = field(0, j, k);
        field(ni - 1 + h, j, k) = field(ni - 1, j, k);
      }
    }
    // South & North halos (including corners, matching WRF order)
    for (int i = ims; i <= ime; ++i) {
      for (int h = 1; h <= halo; ++h) {
        field(i, -h, k) = field(i, 0, k);
        field(i, nj - 1 + h, k) = field(i, nj - 1, k);
      }
    }
  }
}

template <class Real>
void set_physical_bc2d_flat(
    View2<Real> field,
    char var,
    int nx, int ny,
    int halo = 3) {

  const bool is_u = (var == 'u' || var == 'U' || var == 'x' || var == 'X');
  const bool is_v = (var == 'v' || var == 'V' || var == 'y' || var == 'Y');
  const int ni = is_u ? (nx + 1) : nx;
  const int nj = is_v ? (ny + 1) : ny;

  const int ims = -halo, ime = ni - 1 + halo;
  const int jms = -halo, jme = nj - 1 + halo;

  // West & East halos
  for (int j = 0; j < nj; ++j) {
    for (int h = 1; h <= halo; ++h) {
      field(-h, j) = field(0, j);
      field(ni - 1 + h, j) = field(ni - 1, j);
    }
  }
  // South & North halos
  for (int i = ims; i <= ime; ++i) {
    for (int h = 1; h <= halo; ++h) {
      field(i, -h) = field(i, 0);
      field(i, nj - 1 + h) = field(i, nj - 1);
    }
  }
}

}  // namespace cppwrf::bdy
