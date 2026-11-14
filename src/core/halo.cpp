// SPDX-License-Identifier: GPL-3.0-or-later
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
  // the indices held: the whole domain, or the worker's patch (its points, the last face and the lateral halo at
  // the domain's edges) and Grid::halo points around it
  L.wi0 = -Grid::halo;
  L.wi1 = L.nx + Grid::halo;
  L.wj0 = -L.hy;
  L.wj1 = L.ny + L.hy;
  if (g.patched()) {
    L.dom = g.dom;
    L.wi0 = g.pi0 - Grid::halo;
    L.wi1 = (g.pi1 == g.nx ? L.nx : g.pi1) + Grid::halo;
    if (!g.two_d()) {
      L.wj0 = g.pj0 - L.hy;
      L.wj1 = (g.pj1 == g.ny ? L.ny : g.pj1) + L.hy;
    }
  }
  // a row: lpad elements, the points, the halo; the first point after the west halo is aligned (element lpad)
  std::ptrdiff_t pitch = (lpad + (L.wi1 - L.wi0 - Grid::halo) + 31) / 32 * 32;
  L.si = 1;
  L.sj = g.two_d() ? 0 : pitch;
  std::ptrdiff_t plane = pitch * (L.wj1 - L.wj0);
  L.sk = rank == 3 ? plane : 0;
  L.offset = lpad - (L.wi0 + Grid::halo) - L.wj0 * L.sj;
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
  bool odd_face = false;  // symmetric: the wall face itself changes sign too (Parity::OddFace)
  int lo = -(1 << 30), hi = 1 << 30;  // the halo points lo <= i < hi only (a patch: those its allocation holds)
  CPPWRF_HD void operator()(int a, int b, int k) const {
    const int nl = face ? 4 : 3;
    int i = a < nl ? (face ? -a : -1 - a) : n + (a - nl);
    if (i < lo || i >= hi) return;
    int src;
    if (flat) {
      const int last = face ? n : n - 1;
      src = i < 0 ? 0 : (i > last ? last : i);
    } else if (periodic) {
      src = ((i % n) + n) % n;  // modular: also for n < halo width (WRF's sequential copy wraps the same way)
    } else if (face) {
      src = i <= 0 ? -i : 2 * n - i;
    } else {
      src = i < 0 ? -1 - i : 2 * n - 1 - i;
    }
    real* base = p + b * st + k * sk;
    real v = base[src * sn];
    // symmetric: the wall face itself (src == i) keeps its value (WRF set_physical_bc3d fills only
    // ids-1.. / ide+1..: u(ids-i) = -u(ids+i)); odd fields change sign in the mirror
    if (!periodic && !flat && src == i) {
      if (odd_face) base[i * sn] = -v;
      return;
    }
    if (!periodic && !flat && odd) v = -v;
    base[i * sn] = v;
  }
};

struct HaloLineT {  // HaloLine with (tangential, slot) thread order
  HaloLine h;
  CPPWRF_HD void operator()(int b, int a, int k) const { h(a, b, k); }
};

// Several fields in one launch per direction; field f and level k are folded into the third
// index (kk = f * nkmax + k). In the y pass the x index is the fastest thread index (coalesced).
constexpr int kMaxFields = 8;
struct HaloMulti {
  int nf, nkmax;
  real* p[kMaxFields];
  std::ptrdiff_t st[kMaxFields], sk[kMaxFields];
  int nk[kMaxFields];
  int t0[kMaxFields], t1[kMaxFields];  // tangential range of field f
  int lo[kMaxFields], hi[kMaxFields];  // the halo points along the filled direction (HaloLine::lo, hi)
  bool face[kMaxFields], odd[kMaxFields], odd_face[kMaxFields];
  std::ptrdiff_t sn;  // stride along the filled direction (1 for x, pitch for y): per field below
  std::ptrdiff_t snf[kMaxFields];
  int n; bool periodic, ydir, flat;
  // Both passes are launched with the tangential index as the fastest thread index and the
  // halo slot (8) as the second one: all lanes of a warp busy (the slot range is only 8 wide).
  CPPWRF_HD void operator()(int x, int y, int kk) const {
    const int f = kk / nkmax, k = kk - f * nkmax;
    if (f >= nf || k >= nk[f]) return;
    const int a = y, b = x;  // a = halo slot, b = tangential index
    const int na = face[f] ? 8 : 6;
    if (a >= na) return;
    if (b < t0[f] || b >= t1[f]) return;
    HaloLine h{p[f], snf[f], st[f], sk[f], n, face[f], periodic, odd[f], flat, odd_face[f], lo[f], hi[f]};
    h(a, b, k);
  }
};

// A decomposed run (core/tiles.hpp). First the exchange: the points the other patches own inside this patch's
// allocation. Then, x before y: the boundary fill of the direction on the halo points the allocation holds (their
// sources are inside it too), over its tangential range; a periodic direction with several patches along it takes
// the halo from the patch that owns the wrapped point instead (tile::hooks.wrap). The halo lines' index space is not
// (i, j): the launches are not clipped to the patch (tile::Whole). No decomposition: the fill of the whole halo.
void exchange_bases(const void* const* b, int n) {
  if (tile::ctx.patch && tile::hooks.exchange) tile::hooks.exchange(b, n, false);
}
bool wraps(const Layout& L, int axis) { return tile::ctx.patch && tile::hooks.wraps && tile::hooks.wraps(L.dom, axis); }

void fill_impl(real* p, const Layout& L, const Grid& g, bool ux, bool vy, Parity px, Parity py) {
  const void* const base = p - L.offset;
  exchange_bases(&base, 1);
  const tile::Whole whole;
  const int nk = L.nz;
  if (wraps(L, 0)) {
    tile::hooks.wrap(&base, 1, 0);
    tile::after_launch("halo wrap x");  // the check counts the launch of a run with one patch
  } else {  // x direction over interior rows
    HaloLine h{p, 1, L.sj, L.sk, g.nx, ux, g.bc_x == BC::Periodic, px != Parity::Even, g.flat_x(), px == Parity::OddFace,
               L.wi0, L.wi1};
    int na = ux ? 8 : 6;
    parallel_for(Range3D{L.j0(), L.j1(), 0, na, 0, nk}, HaloLineT{h});  // rows as the fastest thread index
  }
  if (g.two_d()) return;
  if (wraps(L, 1)) {
    tile::hooks.wrap(&base, 1, 1);
    tile::after_launch("halo wrap y");
  } else {  // y direction over full x extent incl. halos -> corners
    HaloLine h{p, L.sj, 1, L.sk, g.ny, vy, g.bc_y == BC::Periodic, py != Parity::Even, g.flat_y(), py == Parity::OddFace,
               L.wj0, L.wj1};
    int na = vy ? 8 : 6;
    parallel_for(Range3D{0, na, L.wi0, L.wi1, 0, nk}, h);
  }
}
}  // namespace

static void fill_many(const std::vector<HaloField>& v, const Grid& g) {
  if (v.empty()) return;
  std::vector<const void*> b;
  if (tile::ctx.patch) {
    b.reserve(v.size());
    for (const HaloField& f : v) b.push_back(f.p - f.L.offset);
    exchange_bases(b.data(), int(b.size()));
  }
  const tile::Whole whole;
  const bool wrap_x = wraps(v[0].L, 0), wrap_y = !g.two_d() && wraps(v[0].L, 1);
  // each field's x pass before its y pass: every x pass, then every y pass (the fields are distinct)
  for (int pass = 0; pass < (g.two_d() ? 1 : 2); ++pass) {
    if (pass == 0 ? wrap_x : wrap_y) {
      tile::hooks.wrap(b.data(), int(b.size()), pass);
      for (std::size_t s0 = 0; s0 < v.size(); s0 += kMaxFields) tile::after_launch("halo wrap");  // as the launches below
      continue;
    }
    for (std::size_t s0 = 0; s0 < v.size(); s0 += kMaxFields) {
      const int nf = int(std::min<std::size_t>(kMaxFields, v.size() - s0));
      HaloMulti h{};
      int nkmax = 1, ta = 1 << 30, tb = -(1 << 30);
      for (int f = 0; f < nf; ++f) {
        const HaloField& fld = v[s0 + f];
        const Layout& L = fld.L;
        const bool ux = fld.st == Stagger::U, vy = fld.st == Stagger::V;
        nkmax = std::max(nkmax, L.nz);
        h.p[f] = fld.p;
        h.sk[f] = L.sk;
        h.nk[f] = L.nz;
        if (pass == 0) {  // x pass: along i (stride 1), tangential j over interior rows (as fill_impl)
          h.snf[f] = 1; h.st[f] = L.sj; h.t0[f] = L.j0(); h.t1[f] = L.j1(); h.face[f] = ux;
          h.lo[f] = L.wi0; h.hi[f] = L.wi1;
          h.odd[f] = fld.px != Parity::Even; h.odd_face[f] = fld.px == Parity::OddFace;
        } else {  // y pass: along j (stride sj), tangential i over the full x extent incl. halos
          h.snf[f] = L.sj; h.st[f] = 1; h.t0[f] = L.wi0; h.t1[f] = L.wi1; h.face[f] = vy;
          h.lo[f] = L.wj0; h.hi[f] = L.wj1;
          h.odd[f] = fld.py != Parity::Even; h.odd_face[f] = fld.py == Parity::OddFace;
        }
        ta = std::min(ta, h.t0[f]);
        tb = std::max(tb, h.t1[f]);
      }
      h.nf = nf;
      h.nkmax = nkmax;
      if (pass == 0) { h.n = g.nx; h.periodic = g.bc_x == BC::Periodic; h.ydir = false; h.flat = g.flat_x(); }
      else { h.n = g.ny; h.periodic = g.bc_y == BC::Periodic; h.ydir = true; h.flat = g.flat_y(); }
      parallel_for(Range3D{ta, tb, 0, 8, 0, nf * nkmax}, h);
    }
  }
}

template <class F> static void fill_ptrs(const std::vector<F*>& fs, const Grid& g) {
  std::vector<HaloField> v;
  v.reserve(fs.size());
  for (F* f : fs) v.emplace_back(*f);
  fill_many(v, g);
}

void fill_halos(std::initializer_list<Field3*> fs, const Grid& g) { fill_ptrs(std::vector<Field3*>(fs), g); }
void fill_halos(const std::vector<Field3*>& fs, const Grid& g) { fill_ptrs(fs, g); }
void fill_halos(std::initializer_list<Field2*> fs, const Grid& g) { fill_ptrs(std::vector<Field2*>(fs), g); }
void fill_halos(std::initializer_list<HaloField> fs, const Grid& g) { fill_many(std::vector<HaloField>(fs), g); }
void fill_halos(const std::vector<HaloField>& fs, const Grid& g) { fill_many(fs, g); }

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
