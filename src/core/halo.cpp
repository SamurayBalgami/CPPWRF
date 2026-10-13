#include "core/halo.hpp"

#include <algorithm>
#include <initializer_list>
#include <type_traits>
#include <stdexcept>
#include <vector>

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
  bool flat = false;  // specified: copy the outermost row outward (WRF set_physical_bc*)
  CPPWRF_HD void operator()(int a, int b, int k) const {
    const int nl = face ? 4 : 3;
    int i = a < nl ? (face ? -a : -1 - a) : n + (a - nl);
    int src;
    if (flat) {
      const int last = face ? n : n - 1;
      src = i < 0 ? 0 : (i > last ? last : i);
    } else if (periodic) {
      src = i < 0 ? i + n : (i >= n ? i - n : i);
    } else if (face) {
      src = i <= 0 ? -i : 2 * n - i;
    } else {
      src = i < 0 ? -1 - i : 2 * n - 1 - i;
    }
    real* base = p + b * st + k * sk;
    real v = base[src * sn];
    if (!periodic && !flat && odd) v = (src == i) ? real(0) : -v;
    base[i * sn] = v;
  }
};

// Several fields in one launch per direction; field f and level k are folded into the third
// index (kk = f * nkmax + k). In the y pass the x index is the fastest thread index (coalesced).
constexpr int kMaxFields = 8;
struct HaloMulti {
  int nf, nkmax;
  real* p[kMaxFields];
  std::ptrdiff_t st[kMaxFields], sk[kMaxFields];
  int nk[kMaxFields], nt[kMaxFields];
  bool face[kMaxFields], odd[kMaxFields];
  std::ptrdiff_t sn;  // stride along the filled direction (1 for x, pitch for y): per field below
  std::ptrdiff_t snf[kMaxFields];
  int n; bool periodic, ydir, flat;
  CPPWRF_HD void operator()(int x, int y, int kk) const {
    const int f = kk / nkmax, k = kk - f * nkmax;
    if (f >= nf || k >= nk[f]) return;
    const int a = ydir ? y : x, b = ydir ? x : y;  // a = halo slot, b = tangential index
    const int na = face[f] ? 8 : 6;
    if (a >= na) return;
    if (ydir ? (b < -3 || b >= nt[f] + 3) : (b >= nt[f])) return;
    HaloLine h{p[f], snf[f], st[f], sk[f], n, face[f], periodic, odd[f], flat};
    h(a, b, k);
  }
};

void fill_impl(real* p, const Layout& L, const Grid& g, bool ux, bool vy, Parity px, Parity py) {
  const int nk = L.nz;
  {  // x direction over interior rows
    HaloLine h{p, 1, L.sj, L.sk, g.nx, ux, g.bc_x == BC::Periodic, px == Parity::Odd, g.bc_x == BC::Specified};
    int na = ux ? 8 : 6;
    parallel_for(Range3D{0, na, 0, L.ny, 0, nk}, h);
  }
  if (!g.two_d()) {  // y direction over full x extent incl. halos -> corners
    HaloLine h{p, L.sj, 1, L.sk, g.ny, vy, g.bc_y == BC::Periodic, py == Parity::Odd, g.bc_y == BC::Specified};
    int na = vy ? 8 : 6;
    parallel_for(Range3D{0, na, -Grid::halo, L.nx + Grid::halo, 0, nk}, h);
  }
}
}  // namespace

template <class F> static void fill_many(std::initializer_list<F*> fs, const Grid& g) {
  std::vector<F*> v(fs);
  for (std::size_t s0 = 0; s0 < v.size(); s0 += kMaxFields) {
    const int nf = int(std::min<std::size_t>(kMaxFields, v.size() - s0));
    HaloMulti hx{}, hy{};
    int nkmax = 1, ntx = 1, nty = 1;
    for (int f = 0; f < nf; ++f) {
      F& fld = *v[s0 + f];
      const Layout& L = fld.layout();
      const bool ux = fld.stagger() == Stagger::U, vy = fld.stagger() == Stagger::V;
      nkmax = std::max(nkmax, L.nz);
      hx.p[f] = hy.p[f] = fld.data();
      hx.sk[f] = hy.sk[f] = L.sk;
      hx.nk[f] = hy.nk[f] = L.nz;
      // x pass: along i (stride 1), tangential j over interior rows
      // default parity as fill_halo(f, g): Field3 U odd in x, V odd in y; Field2 always even
      constexpr bool is3 = std::is_same_v<F, Field3>;
      hx.snf[f] = 1; hx.st[f] = L.sj; hx.nt[f] = L.ny; hx.face[f] = ux; hx.odd[f] = is3 && ux;
      // y pass: along j (stride sj), tangential i over the full x extent incl. halos
      hy.snf[f] = L.sj; hy.st[f] = 1; hy.nt[f] = L.nx; hy.face[f] = vy; hy.odd[f] = is3 && vy;
      ntx = std::max(ntx, L.ny);
      nty = std::max(nty, L.nx);
    }
    hx.nf = hy.nf = nf;
    hx.nkmax = hy.nkmax = nkmax;
    hx.n = g.nx; hx.periodic = g.bc_x == BC::Periodic; hx.ydir = false; hx.flat = g.bc_x == BC::Specified;
    hy.n = g.ny; hy.periodic = g.bc_y == BC::Periodic; hy.ydir = true; hy.flat = g.bc_y == BC::Specified;
    parallel_for(Range3D{0, 8, 0, ntx, 0, nf * nkmax}, hx);
    if (!g.two_d()) parallel_for(Range3D{-Grid::halo, nty + Grid::halo, 0, 8, 0, nf * nkmax}, hy);
  }
}

void fill_halos(std::initializer_list<Field3*> fs, const Grid& g) { fill_many(fs, g); }
void fill_halos(std::initializer_list<Field2*> fs, const Grid& g) { fill_many(fs, g); }

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
