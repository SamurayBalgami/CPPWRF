// SPDX-License-Identifier: GPL-3.0-or-later
// LaunchBatch and parallel_fused (several independent kernels in one launch) must give exactly what parallel_for of
// each kernel gives: parts with different ranges (staggers, k offsets, one-row and 2-D ranges), a batch that fills
// up and launches by itself, empty parts, and kernels of two types in one launch (the CUDA FP32 build batches only
// kernels without multiply-add: KCopy).
#include <cstdio>
#include <vector>

#include "core/backend.hpp"
#include "core/field.hpp"

using namespace cppwrf;
static int fails = 0;

namespace {
struct KMix {  // out = in * a + (i - 2 j + 3 k) * b, a different value at every point
  View3<real> out;
  View3<const real> in;
  real a, b;
  CPPWRF_HD void operator()(int i, int j, int k) const { out(i, j, k) = in(i, j, k) * a + real(i - 2 * j + 3 * k) * b; }
};
struct KNeg {  // a second kernel type: out = -in at (i, j, k) and its neighbour i - 1
  View3<real> out;
  View3<const real> in;
  CPPWRF_HD void operator()(int i, int j, int k) const { out(i, j, k) = -(in(i, j, k) + real(0.5) * in(i - 1, j, k)); }
};
struct KCopy {  // no multiply-add: batched in the CUDA FP32 build too
  static constexpr bool kNoMulAdd = true;
  View3<real> out;
  View3<const real> in;
  CPPWRF_HD void operator()(int i, int j, int k) const { out(i, j, k) = -in(i, j - 1 < 0 ? 0 : j - 1, k); }
};
std::vector<real> host(const Field3& f) {
  std::vector<real> h;
  f.download(h);
  return h;
}
}  // namespace

static void run(int nx, int ny, int nz) {
  Grid g;
  g.nx = nx; g.ny = ny; g.nz = nz;
  auto val = [](int i, int j, int k) { return 0.37 * i - 1.1 * j + 3.3 * k + 0.01 * i * j * (k + 1); };
  const Stagger st[] = {Stagger::M, Stagger::U, Stagger::V, Stagger::W, Stagger::M};
  const int nf = 5;
  std::vector<Field3> in, ref, bat;
  for (int f = 0; f < nf; ++f) {
    in.emplace_back(g, st[f], "in");
    ref.emplace_back(g, st[f], "ref");
    bat.emplace_back(g, st[f], "bat");
    in[f].fill_host(val);
  }
  // part ranges: the whole staggered field, a k sub-range, one row, a single level (2-D), an interior box
  auto range = [&](int f) {
    Range3D r = g.cells(st[f]);
    if (f == 1) { r.k0 = 1; r.k1 = nz - 1; }
    if (f == 2) { r.j0 = r.j1 - 1; }
    if (f == 3) { r.k0 = nz; r.k1 = nz + 1; }
    if (f == 4) { r.i0 = 1; r.i1 = nx - 1; r.j0 = ny > 2 ? 1 : 0; r.j1 = ny > 2 ? ny - 1 : ny; }
    return r;
  };
  for (int f = 0; f < nf; ++f)
    parallel_for(range(f), KMix{ref[f].view(), in[f].cview(), real(1.5 + f), real(0.25 * f - 1)});
  {
    LaunchBatch<KMix, 3> b;  // 5 parts: the batch launches by itself after 3, run() launches the other 2
    for (int f = 0; f < nf; ++f) {
      b.add(Range3D{0, 0, 0, 1, 0, 1}, KMix{bat[f].view(), in[f].cview(), real(0), real(0)});  // empty: no launch
      b.add(range(f), KMix{bat[f].view(), in[f].cview(), real(1.5 + f), real(0.25 * f - 1)});
    }
    b.run();
  }
  int bad = 0;
  for (int f = 0; f < nf; ++f) bad += host(ref[f]) != host(bat[f]);
  for (int f = 0; f < nf; ++f) parallel_for(range(f), KCopy{ref[f].view(), in[f].cview()});
  {
    LaunchBatch<KCopy, 2> b;
    for (int f = 0; f < nf; ++f) b.add(range(f), KCopy{bat[f].view(), in[f].cview()});
    b.run();
  }
  for (int f = 0; f < nf; ++f) bad += host(ref[f]) != host(bat[f]);
  // two kernel types in one launch
  Field3 r1(g, Stagger::U, "r1"), r2(g, Stagger::M, "r2"), b1(g, Stagger::U, "b1"), b2(g, Stagger::M, "b2");
  const Range3D ru{1, nx + 1, 0, ny, 0, nz}, rm{0, nx, 0, ny, 2, nz};
  parallel_for(ru, KNeg{r1.view(), in[1].cview()});
  parallel_for(rm, KMix{r2.view(), in[0].cview(), real(-0.75), real(2)});
  parallel_fused(part(ru, KNeg{b1.view(), in[1].cview()}), part(rm, KMix{b2.view(), in[0].cview(), real(-0.75), real(2)}));
  bad += host(r1) != host(b1);
  bad += host(r2) != host(b2);
  if (bad) {
    std::printf("FAIL launch batch nx=%d ny=%d nz=%d: %d fields differ\n", nx, ny, nz, bad);
    ++fails;
  }
}

int main() {
  run(13, 1, 5);
  run(13, 9, 5);
  run(70, 33, 7);
  std::printf(fails ? "test_launch_batch: %d FAILED\n" : "test_launch_batch: all passed\n", fails);
  return fails ? 1 : 0;
}
