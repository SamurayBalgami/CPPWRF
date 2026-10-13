// fill_halos (several fields, one launch per direction) must give exactly the same halos as
// fill_halo field by field: all staggers, Field2 and Field3, periodic and symmetric, 2D and 3D.
#include <cstdio>
#include <vector>

#include "core/field.hpp"
#include "core/halo.hpp"

using namespace cppwrf;
static int fails = 0;

template <int R> static std::vector<real> host(const FieldT<R>& f) { std::vector<real> h; f.download(h); return h; }

static void run(int nx, int ny, int nz, BC bx, BC by) {
  Grid g; g.nx = nx; g.ny = ny; g.nz = nz; g.bc_x = bx; g.bc_y = by;
  auto val = [](int i, int j, int k) { return 0.37 * i - 1.1 * j + 3.3 * k + 0.01 * i * j * (k + 1); };
  Field3 a1(g, Stagger::M, "m"), b1(g, Stagger::U, "u"), c1(g, Stagger::V, "v"), d1(g, Stagger::W, "w");
  Field3 a2(g, Stagger::M, "m"), b2(g, Stagger::U, "u"), c2(g, Stagger::V, "v"), d2(g, Stagger::W, "w");
  Field2 e1(g, Stagger::M, "mu"), f1(g, Stagger::U, "muu"), e2(g, Stagger::M, "mu"), f2(g, Stagger::U, "muu");
  for (Field3* f : {&a1, &b1, &c1, &d1, &a2, &b2, &c2, &d2}) f->fill_host(val);
  for (Field2* f : {&e1, &f1, &e2, &f2}) f->fill_host(val);
  for (Field3* f : {&a1, &b1, &c1, &d1}) fill_halo(*f, g);
  for (Field2* f : {&e1, &f1}) fill_halo(*f, g);
  fill_halos({&a2, &b2, &c2, &d2}, g);
  fill_halos({&e2, &f2}, g);
  int bad = 0;
  const std::pair<const Field3*, const Field3*> p3[] = {{&a1, &a2}, {&b1, &b2}, {&c1, &c2}, {&d1, &d2}};
  for (auto [x, y] : p3) bad += host(*x) != host(*y);
  bad += host(e1) != host(e2);
  bad += host(f1) != host(f2);
  if (bad) { std::printf("FAIL fill_halos nx=%d ny=%d bx=%d by=%d: %d fields differ\n", nx, ny, int(bx), int(by), bad); ++fails; }
}

int main() {
  for (BC bx : {BC::Periodic, BC::Symmetric})
    for (BC by : {BC::Periodic, BC::Symmetric}) {
      run(13, 1, 5, bx, by);
      run(13, 9, 5, bx, by);
      run(40, 33, 7, bx, by);
    }
  std::printf(fails ? "test_halo_multi: %d FAILED\n" : "test_halo_multi: all passed\n", fails);
  return fails ? 1 : 0;
}
