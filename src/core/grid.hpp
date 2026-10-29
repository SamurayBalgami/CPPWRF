#pragma once
// Grid geometry and staggering. Arakawa-C horizontally, Lorenz vertically (TN 3.1).
// Indexing: mass point (i,j,k), 0 <= i < nx, 0 <= j < ny, 0 <= k < nz (k=0 lowest).
//   U point (i,j,k) is at the west face of mass cell i  (0 <= i <= nx)
//   V point (i,j,k) is at the south face of mass cell j (0 <= j <= ny; ny==1 -> j=0 only)
//   W point (i,j,k) is at the bottom face of mass cell k (0 <= k <= nz)
// 2D x-z runs: ny == 1, the j-stride is 0 so f(i,j+-1,k) == f(i,j,k) (periodic y for free).
#include "core/backend.hpp"
#include "core/precision.hpp"

namespace cppwrf {

enum class Stagger : char { M = 'M', U = 'U', V = 'V', W = 'W' };
// Specified: real-data lateral boundaries (WRF specified = .true.); halos are filled by copying
// the outermost row outward (set_physical_bc*), boundary rows are set by src/dyn/lateral_bc.
// Open: WRF open_xs/xe (open_ys/ye) on both sides of the direction: the same outward halo copy as
// specified, radiative normal velocity and flow-dependent (inflow zero-gradient) scalars at the edges.
enum class BC { Periodic, Symmetric, Specified, Open };

struct Grid {
  int nx = 1, ny = 1, nz = 1;  // number of mass points
  double dx = 1, dy = 1;       // m
  double dt = 1;               // s (large step)
  double ztop = 0;             // m (nominal, for flat terrain)
  double p_top = 0;            // Pa
  BC bc_x = BC::Periodic, bc_y = BC::Periodic;
  static constexpr int halo = 3;

  bool two_d() const { return ny == 1; }
  int hy() const { return two_d() ? 0 : halo; }
  int nxs(Stagger s) const { return s == Stagger::U ? nx + 1 : nx; }
  int nys(Stagger s) const { return (s == Stagger::V && !two_d()) ? ny + 1 : ny; }
  int nzs(Stagger s) const { return s == Stagger::W ? nz + 1 : nz; }
  Range3D cells(Stagger s) const { return {0, nxs(s), 0, nys(s), 0, nzs(s)}; }
  Range3D mass_cells() const { return {0, nx, 0, ny, 0, nz}; }
  Range2D columns() const { return {0, nx, 0, ny}; }
  // Range of U (V) points updated by the dynamics: boundary faces stay 0 for rigid walls,
  // and u(nx) is a halo copy of u(0) for periodic.
  int u_i0() const { return bc_x != BC::Periodic ? 1 : 0; }
  int u_i1() const { return nx; }
  int v_j0() const { return (bc_y != BC::Periodic && !two_d()) ? 1 : 0; }
  bool specified() const { return bc_x == BC::Specified || bc_y == BC::Specified; }
  bool open_x() const { return bc_x == BC::Open; }
  bool open_y() const { return bc_y == BC::Open && !two_d(); }
  // halo by outward copy of the outermost row (WRF set_physical_bc*: open, specified, nested)
  bool flat_x() const { return bc_x == BC::Specified || bc_x == BC::Open; }
  bool flat_y() const { return bc_y == BC::Specified || bc_y == BC::Open; }
  int v_j1() const { return ny; }
};

}  // namespace cppwrf
