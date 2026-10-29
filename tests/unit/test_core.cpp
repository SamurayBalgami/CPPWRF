// Core unit tests: layout/indexing, staggering extents, halo fill (periodic, symmetric,
// parity), 2D y-periodicity, uniform base-state views, vertical coefficients, Reducer.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "core/field.hpp"
#include "core/halo.hpp"
#include "dyn/state.hpp"

using namespace cppwrf;
static int fails = 0;
#define CHECK(c)                                                       \
  do {                                                                 \
    if (!(c)) { std::printf("FAIL %s:%d %s\n", __FILE__, __LINE__, #c); ++fails; } \
  } while (0)

static double val(int i, int j, int k) { return 1 + i + 100 * j + 10000 * k; }

struct Host {  // host copy of a field with logical indexing (halos included)
  std::vector<real> h;
  Layout L;
  double operator()(int i, int j, int k) const { return h[L.index(i, j, k)]; }
};
template <int R> Host get(const FieldT<R>& f) { Host x; f.download(x.h); x.L = f.layout(); return x; }

static Grid mkgrid(int nx, int ny, int nz, BC bx, BC by) {
  Grid g; g.nx = nx; g.ny = ny; g.nz = nz; g.bc_x = bx; g.bc_y = by; return g;
}

static void test_layout() {
  Grid g = mkgrid(37, 5, 4, BC::Periodic, BC::Periodic);
  for (Stagger s : {Stagger::M, Stagger::U, Stagger::V, Stagger::W}) {
    Field3 f(g, s, "f");
    const Layout& L = f.layout();
    CHECK(L.sj % 32 == 0);
    CHECK(L.nx == g.nxs(s) && L.ny == g.nys(s) && L.nz == g.nzs(s));
    CHECK(reinterpret_cast<std::uintptr_t>(f.data()) % (8 * sizeof(real)) == 0 || sizeof(real) == 8);
    f.fill_host(val);
    Host x = get(f);
    for (int k = 0; k < L.nz; ++k)
      for (int j = 0; j < L.ny; ++j)
        for (int i = 0; i < L.nx; ++i) CHECK(x(i, j, k) == real(val(i, j, k)));
  }
  CHECK(g.nxs(Stagger::U) == 38 && g.nys(Stagger::V) == 6 && g.nzs(Stagger::W) == 5);
  Grid g2 = mkgrid(16, 1, 4, BC::Periodic, BC::Periodic);
  Field3 f(g2, Stagger::V, "v");
  CHECK(f.layout().sj == 0 && f.layout().ny == 1);
  f.fill_host(val);
  Host x = get(f);
  CHECK(x(3, 1, 2) == x(3, 0, 2) && x(3, -1, 2) == x(3, 0, 2));  // free y periodicity
}

static void test_halo_periodic() {
  Grid g = mkgrid(10, 6, 3, BC::Periodic, BC::Periodic);
  for (Stagger s : {Stagger::M, Stagger::U, Stagger::V, Stagger::W}) {
    Field3 f(g, s, "f");
    f.fill_host([&](int i, int j, int k) { return val(i % g.nx, j % g.ny, k); });
    fill_halo(f, g);
    Host x = get(f);
    int nzs = g.nzs(s);
    for (int k = 0; k < nzs; ++k)
      for (int j = -3; j < g.nys(s) + 3; ++j)
        for (int i = -3; i < g.nxs(s) + 3; ++i) {
          int ii = ((i % g.nx) + g.nx) % g.nx, jj = ((j % g.ny) + g.ny) % g.ny;
          CHECK(x(i, j, k) == real(val(ii, jj, k)));
        }
  }
}

static void test_halo_symmetric() {
  Grid g = mkgrid(10, 6, 2, BC::Symmetric, BC::Symmetric);
  Field3 m(g, Stagger::M, "m"), u(g, Stagger::U, "u"), v(g, Stagger::V, "v");
  m.fill_host(val); u.fill_host(val); v.fill_host(val);
  fill_halo(m, g); fill_halo(u, g); fill_halo(v, g);
  Host xm = get(m), xu = get(u), xv = get(v);
  for (int k = 0; k < 2; ++k)
    for (int j = 0; j < 6; ++j)
      for (int t = 1; t <= 3; ++t) {
        CHECK(xm(-t, j, k) == xm(t - 1, j, k));
        CHECK(xm(9 + t, j, k) == xm(10 - t, j, k));
        CHECK(xu(-t, j, k) == -xu(t, j, k));
        CHECK(xu(10 + t, j, k) == -xu(10 - t, j, k));
        CHECK(xm(4, -t, k) == xm(4, t - 1, k));       // y mirror, mass
        CHECK(xv(4, -t, k) == -xv(4, t, k));          // y odd, V
        CHECK(xv(4, 6 + t, k) == -xv(4, 6 - t, k));
        CHECK(xu(-t, -t, k) == -xu(t, t - 1, k));     // corner: x odd, y even
      }
  // the wall face of the normal velocity keeps its value (WRF fills only the mirror rows)
  for (int j = 0; j < 6; ++j) { CHECK(xu(0, j, 0) == real(val(0, j, 0))); CHECK(xu(10, j, 1) == real(val(10, j, 1))); }
  for (int i = 0; i < 10; ++i) { CHECK(xv(i, 0, 0) == real(val(i, 0, 0))); CHECK(xv(i, 6, 1) == real(val(i, 6, 1))); }
  // even parity on a U-staggered 2D field (e.g. muu)
  Field2 mu(g, Stagger::U, "muu");
  mu.fill_host(val);
  fill_halo(mu, g);
  Host x = get(mu);
  CHECK(x(0, 2, 0) == real(val(0, 2, 0)) && x(-2, 2, 0) == x(2, 2, 0));
}

static void test_base_and_levels() {
  Grid g = mkgrid(8, 1, 10, BC::Periodic, BC::Periodic);
  g.dx = 100; g.p_top = 1e4;
  dyn::BaseState b(g, true);
  b.pb.fill_host([](int, int, int k) { return 1000.0 * k; });
  std::vector<real> h(b.pb.layout().total);
  copy_to_host(h.data(), b.pb.data(), h.size() * sizeof(real));
  auto v = b.pb.view();
  BView3<const real> hv{h.data(), v.si, v.sj, v.sk};
  CHECK(v.si == 0 && v.sj == 0 && hv(5, 3, 7) == real(7000));
  std::vector<double> e(11);
  for (int k = 0; k <= 10; ++k) e[k] = 1.0 - std::pow(k / 10.0, 1.3);
  e[10] = 0.0;
  b.set_levels(e);
  for (int k = 1; k < 10; ++k) {
    CHECK(std::fabs(b.fnm[k] + b.fnp[k] - 1) < 1e-14);
    // linear-in-eta interpolation is exact
    CHECK(std::fabs(b.fnm[k] * b.znu[k] + b.fnp[k] * b.znu[k - 1] - b.znw[k]) < 1e-14);
  }
  CHECK(std::fabs(b.cf1 + b.cf2 + b.cf3 - 1) < 1e-13);
  CHECK(std::fabs(b.cfn + b.cfn1 - 1) < 1e-13);
  CHECK(std::fabs(b.cfn * b.znu[9] + b.cfn1 * b.znu[8] - b.znw[10]) < 1e-13);  // linear exact
  CHECK(std::fabs(b.cf1 * b.znu[0] + b.cf2 * b.znu[1] + b.cf3 * b.znu[2] - b.znw[0]) < 1e-13);
}

struct SumF {
  View3<const real> f;
  CPPWRF_HD double operator()(int i, int j, int k) const { return double(f(i, j, k)) * 1.000001; }
};

static void test_reducer() {
  Grid g = mkgrid(64, 9, 7, BC::Periodic, BC::Periodic);
  Field3 f(g, Stagger::M, "f");
  f.fill_host([](int i, int j, int k) { return std::sin(0.37 * i + 1.1 * j + 0.3 * k); });
  Reducer r(16);
  double a = r.sum(g.mass_cells(), SumF{f.cview()});
  double b = r.sum(g.mass_cells(), SumF{f.cview()});
  double m = r.max(g.mass_cells(), SumF{f.cview()});
  CHECK(a == b);
  CHECK(m <= 1.000001 && m > 0.99);
}

int main() {
  test_layout();
  test_halo_periodic();
  test_halo_symmetric();
  test_base_and_levels();
  test_reducer();
  std::printf("test_core (%s, %s): %s\n", backend_name(), precision_name, fails ? "FAILED" : "OK");
  return fails ? 1 : 0;
}
