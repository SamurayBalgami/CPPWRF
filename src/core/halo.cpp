#include "core/halo.hpp"

#include <stdexcept>

namespace cppwrf {

Layout Layout::make(const Grid& g, Stagger s, int rank, bool uniform) {
  Layout L;
  L.nx = g.nxs(s);
  L.ny = g.nys(s);
  L.nz = rank == 3 ? g.nzs(s) : 1;
  if (uniform) {  // horizontally uniform base-state storage: one value per level
    L.si = 0; L.sj = 0; L.sk = rank == 3 ? 1 : 0; L.offset = 0;
    L.total = std::size_t(L.nz);
    return L;
  }
  L.hy = g.hy();
  std::ptrdiff_t pitch = (lpad + L.nx + Grid::halo + 31) / 32 * 32;
  L.si = 1;
  L.sj = g.two_d() ? 0 : pitch;
  std::ptrdiff_t plane = pitch * (L.ny + 2 * L.hy);
  L.sk = rank == 3 ? plane : 0;
  L.offset = lpad + L.hy * pitch;
  L.total = std::size_t(plane) * (rank == 3 ? L.nz : 1);
  if (L.total >= (std::size_t(1) << 31)) throw std::runtime_error("field too large for 32-bit view indexing");
  return L;
}

namespace {
// Fills the halo along one direction. a = halo slot (0..nl+nr-1), b = tangential index.
struct HaloLine {
  real* p;
  std::ptrdiff_t sn, st, sk;
  int n;        // number of mass points in this direction
  bool face;    // field staggered on faces in this direction
  bool periodic;
  bool odd;
  CPPWRF_HD void operator()(int a, int b, int k) const {
    const int nl = face ? 4 : 3;
    int i = a < nl ? (face ? -a : -1 - a) : n + (a - nl);
    int src;
    if (periodic) {
      src = i < 0 ? i + n : (i >= n ? i - n : i);
    } else if (face) {
      src = i <= 0 ? -i : 2 * n - i;
    } else {
      src = i < 0 ? -1 - i : 2 * n - 1 - i;
    }
    real* base = p + b * st + k * sk;
    real v = base[src * sn];
    if (!periodic && odd) v = (src == i) ? real(0) : -v;
    base[i * sn] = v;
  }
};

void fill_impl(real* p, const Layout& L, const Grid& g, bool ux, bool vy, Parity px, Parity py) {
  const int nk = L.nz;
  {  // x direction over interior rows
    HaloLine h{p, 1, L.sj, L.sk, g.nx, ux, g.bc_x == BC::Periodic, px == Parity::Odd};
    int na = ux ? 8 : 6;
    parallel_for(Range3D{0, na, 0, L.ny, 0, nk}, h);
  }
  if (!g.two_d()) {  // y direction over full x extent incl. halos -> corners
    HaloLine h{p, L.sj, 1, L.sk, g.ny, vy, g.bc_y == BC::Periodic, py == Parity::Odd};
    int na = vy ? 8 : 6;
    parallel_for(Range3D{0, na, -Grid::halo, L.nx + Grid::halo, 0, nk}, h);
  }
}
}  // namespace

void fill_halo(Field3& f, const Grid& g, Parity px, Parity py) {
  fill_impl(f.data(), f.layout(), g, f.stagger() == Stagger::U, f.stagger() == Stagger::V, px, py);
}
void fill_halo(Field2& f, const Grid& g, Parity px, Parity py) {
  fill_impl(f.data(), f.layout(), g, f.stagger() == Stagger::U, f.stagger() == Stagger::V, px, py);
}
void fill_halo(Field3& f, const Grid& g) {
  fill_halo(f, g, f.stagger() == Stagger::U ? Parity::Odd : Parity::Even,
            f.stagger() == Stagger::V ? Parity::Odd : Parity::Even);
}
void fill_halo(Field2& f, const Grid& g) { fill_halo(f, g, Parity::Even, Parity::Even); }

}  // namespace cppwrf
