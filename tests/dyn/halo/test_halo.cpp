// Differential test of src/core/halo.cpp (fill_halo, fill_halos) against WRF v4.8.0's
// set_physical_bc3d / set_physical_bc2d (halo_bridge.f90): periodic, specified, open and symmetric
// boundaries, every staggering and parity, from the same random interior and the same marked halo,
// compared bit for bit over CPPWRF's whole halo. Halo points WRF does not set (it leaves them as
// they were) are counted apart: CPPWRF may fill them, nothing reads them as WRF's.
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "core/grid.hpp"
#include "core/halo.hpp"

using namespace cppwrf;

extern "C" void phys_bc_c(int nx, int ny, int nz, int bc, char var, int rank, double* dat);

namespace {

constexpr int H = 4;  // WRF's serial memory halo (bdyzone)

std::uint64_t rng = 0x9E3779B97F4A7C15ull;
double rnd(double lo, double hi) {
  rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
  return lo + (hi - lo) * double(rng >> 11) / double(1ull << 53);
}

int failures = 0;
bool same(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }
const char* bc_name(int bc) { return bc == 0 ? "periodic" : bc == 1 ? "specified" : bc == 2 ? "open" : "symmetric"; }
const char* par_name(Parity p) { return p == Parity::Even ? "even" : p == Parity::Odd ? "odd" : "odd-face"; }

// One field: interior random, every halo point marked with its own value; WRF and CPPWRF fill it.
// multi: through fill_halos (default parity) instead of fill_halo with px, py.
template <int R>
void one(int nx, int ny, int nz, int bc, Stagger st, Parity px, Parity py, char var, bool multi) {
  Grid g;
  g.nx = nx; g.ny = ny; g.nz = nz;
  const BC kinds[] = {BC::Periodic, BC::Specified, BC::Open, BC::Symmetric};
  g.bc_x = g.bc_y = kinds[bc];
  FieldT<R> f(g, st, "f");
  const Layout& L = f.layout();
  const int NI = nx + 1 + 2 * H, NJ = ny + 1 + 2 * H, NK = R == 3 ? nz + 1 : 1;
  auto w = [&](int i, int j, int k) { return std::size_t(i + H) + NI * (std::size_t(k) + NK * std::size_t(j + H)); };
  std::vector<double> a(std::size_t(NI) * NJ * NK);
  for (std::size_t n = 0; n < a.size(); ++n) a[n] = -1e30 - double(n);  // marks: WRF untouched
  std::vector<real> h(L.total, real(0));
  const int hx = Grid::halo, hy = L.hy;
  for (int k = 0; k < L.nz; ++k)
    for (int j = -hy; j < L.ny + hy; ++j)
      for (int i = -hx; i < L.nx + hx; ++i) {
        const bool inner = i >= 0 && i < L.nx && j >= 0 && j < L.ny;
        const double v = inner ? rnd(-1, 1) : a[w(i, j, k)];
        h[L.index(i, j, k)] = real(v);
        a[w(i, j, k)] = v;
      }
  copy_to_dev(f.data() - L.offset, h.data(), f.bytes());
  const std::vector<double> before = a;
  phys_bc_c(nx, ny, nz, bc, var, R, a.data());
  if (multi) {
    if constexpr (R == 3) fill_halos({&f}, g);
    else fill_halos({&f}, g);
  } else {
    fill_halo(f, g, px, py);
  }
  f.download(h);
  long bad = 0, untouched = 0, n = 0;
  for (int k = 0; k < L.nz; ++k)
    for (int j = -hy; j < L.ny + hy; ++j)
      for (int i = -hx; i < L.nx + hx; ++i, ++n) {
        const double c = h[L.index(i, j, k)], r = a[w(i, j, k)];
        if (same(c, r)) continue;
        if (same(r, before[w(i, j, k)])) { ++untouched; continue; }
        if (++bad <= 3) std::printf("    (i=%d, j=%d, k=%d): CPPWRF %a  WRF %a\n", i, j, k, c, r);
      }
  const std::string tag = std::string(bc_name(bc)) + " " + (R == 3 ? "3-D " : "2-D ") + char(st) + " '" + var + "' " +
                          (multi ? "fill_halos" : std::string(par_name(px)) + "/" + par_name(py)) + " " +
                          std::to_string(nx) + "x" + std::to_string(ny);
  std::printf("  %-52s %5ld points  %s", tag.c_str(), n, bad ? (std::to_string(bad) + " differ  FAIL").c_str() : "identical");
  if (untouched) std::printf("  (%ld halo points WRF leaves)", untouched);
  std::printf("\n");
  if (bad) ++failures;
}

// as large_step.cpp fills them: muu with (OddFace, Even) and muv with (Even, Odd), as WRF's
// set_physical_bc2d 'u' flips the wall face of x and its 'v' loop in y starts one row out
void grid_cases(int nx, int ny, int nz, int nbc) {
  using P = Parity;
  using S = Stagger;
  for (int bc = 0; bc < nbc; ++bc) {
    one<3>(nx, ny, nz, bc, S::M, P::Even, P::Even, 't', false);
    one<3>(nx, ny, nz, bc, S::U, P::Odd, P::Even, 'u', false);
    one<3>(nx, ny, nz, bc, S::V, P::Even, P::Odd, 'v', false);
    one<3>(nx, ny, nz, bc, S::W, P::Even, P::Even, 'w', false);
    one<3>(nx, ny, nz, bc, S::U, P::Odd, P::Even, 'u', true);
    one<3>(nx, ny, nz, bc, S::V, P::Even, P::Odd, 'v', true);
    one<3>(nx, ny, nz, bc, S::M, P::Even, P::Even, 't', true);
    one<2>(nx, ny, nz, bc, S::M, P::Even, P::Even, 't', false);
    one<2>(nx, ny, nz, bc, S::M, P::Even, P::Even, 't', true);
    one<2>(nx, ny, nz, bc, S::U, P::OddFace, P::Even, 'u', false);
    one<2>(nx, ny, nz, bc, S::V, P::Even, P::Odd, 'v', false);
  }
}

}  // namespace

int main() {
  std::printf("fill_halo vs WRF v4.8.0 set_physical_bc3d / set_physical_bc2d\n");
  grid_cases(9, 7, 3, 4);
  // fewer points than the halo is wide: periodic only (a mirror of 2 points reads WRF's own unset
  // halo, specified/open behave as on the large grid)
  grid_cases(2, 2, 2, 1);
  std::printf(failures ? "FAIL: %d comparisons differ\n" : "PASS: all comparisons bit-identical\n", failures);
  return failures ? 1 : 0;
}
