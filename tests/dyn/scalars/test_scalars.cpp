// Differential test of src/dyn/scalars.cpp (ScalarTransport::run_stage) against WRF v4.8.0's own
// rk_update_scalar_pd, rk_scalar_tend (advect_scalar / advect_scalar_pd) and rk_update_scalar
// (scalar_bridge.f90), from the same random fields, compared bit for bit (-0 differs from +0).
// Periodic, specified, open and symmetric boundaries; every RK stage, original and positive-definite
// advection, 5th/3rd and 6th/4th order. The grid is not square and the map factors vary.
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "core/field.hpp"
#include "core/grid.hpp"
#include "core/halo.hpp"
#include "dyn/scalars.hpp"
#include "dyn/state.hpp"

using namespace cppwrf;
using namespace cppwrf::dyn;

extern "C" {
void scalar_stage_c(int nx, int ny, int nz, int bc, int h_order, int v_order, int adv_opt, int rk_step, double dt_rk,
                    double* q1, double* q2, double* q_tend, const double* ru, const double* rv, const double* ww,
                    const double* muts, const double* mub, const double* mu1, const double* mu2, const double* c1h,
                    const double* c2h, const double* c1f, const double* c2f, const double* fnm, const double* fnp,
                    const double* rdnw, const double* msfux, const double* msfuy, const double* msfvx,
                    const double* msfvy, const double* msftx, const double* msfty, double rdx, double rdy);
}

namespace {

// H: the memory halo of a serial WRF (bdyzone = 4), one point wider than CPPWRF's (Grid::halo = 3):
// advect_scalar_pd forms the limiter one cell outside the tile, the 5th-order flux there reads q(-4)
constexpr int NX = 13, NY = 11, NZ = 7, H = 4;
constexpr int NI = NX + 1 + 2 * H, NJ = NY + 1 + 2 * H, NK = NZ + 1;  // WRF memory dims (ims:ime, kms:kme, jms:jme)
constexpr double DT = 30.0;

// WRF memory index of C++ point (i, j, k): Fortran (i+1, k+1, j+1), ims = jms = 1-H, kms = 1
std::size_t f3(int i, int j, int k) { return std::size_t(i + H) + NI * (std::size_t(k) + NK * std::size_t(j + H)); }
std::size_t f2(int i, int j) { return std::size_t(i + H) + NI * std::size_t(j + H); }

std::uint64_t rng = 0x9E3779B97F4A7C15ull;
double rnd(double lo, double hi) {
  rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
  return lo + (hi - lo) * double(rng >> 11) / double(1ull << 53);
}

void fill(Field3& f, double lo, double hi) { f.fill_host([&](int, int, int) { return rnd(lo, hi); }); }
void fill(Field2& f, double lo, double hi) { f.fill_host([&](int, int, int) { return rnd(lo, hi); }); }

// WRF memory array of a field (halo points CPPWRF does not store are 0; set_physical_bc3d fills
// the scalar's outermost ring in the bridge)
template <int R> std::vector<double> wrf(const FieldT<R>& f) {
  std::vector<real> h;
  f.download(h);
  const Layout& L = f.layout();
  std::vector<double> a(R == 3 ? std::size_t(NI) * NK * NJ : std::size_t(NI) * NJ, 0.0);
  for (int j = -H; j <= NY + H; ++j)
    for (int k = 0; k < (R == 3 ? NK : 1); ++k)
      for (int i = -H; i <= NX + H; ++i) {
        if (k >= L.nz || j < -L.hy || j >= L.ny + L.hy || i < -Grid::halo || i >= L.nx + Grid::halo) continue;
        a[R == 3 ? f3(i, j, k) : f2(i, j)] = h[L.index(i, j, k)];
      }
  return a;
}
std::vector<double> wrf(const BaseField2& f) {
  std::vector<real> h;
  f.download(h);
  const Layout& L = f.layout();
  std::vector<double> a(std::size_t(NI) * NJ, 0.0);
  for (int j = -H; j <= NY + H; ++j)
    for (int i = -H; i <= NX + H; ++i) {
      if (L.si == 0) { a[f2(i, j)] = h[0]; continue; }
      if (j < -L.hy || j >= L.ny + L.hy || i < -Grid::halo || i >= L.nx + Grid::halo) continue;
      a[f2(i, j)] = h[L.index(i, j, 0)];
    }
  return a;
}
std::vector<double> wrf(const Array1& c) {  // (kms:kme)
  std::vector<double> a(NK, 0.0);
  for (int k = 0; k < c.size(); ++k) a[k] = c[k];
  return a;
}

int failures = 0;
bool same(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }
void check(const std::string& what, const Field3& f, const std::vector<double>& ref) {
  std::vector<real> h;
  f.download(h);
  const Layout& L = f.layout();
  long bad = 0, n = 0;
  for (int k = 0; k < L.nz; ++k)
    for (int j = 0; j < L.ny; ++j)
      for (int i = 0; i < L.nx; ++i, ++n) {
        const double c = h[L.index(i, j, k)], w = ref[f3(i, j, k)];
        if (same(c, w)) continue;
        if (++bad <= 3) std::printf("    %s(i=%d, j=%d, k=%d): CPPWRF %a  WRF %a\n", what.c_str(), i, j, k, c, w);
      }
  std::printf("  %-44s %6ld points  %s\n", what.c_str(), n, bad ? (std::to_string(bad) + " differ  FAIL").c_str() : "identical");
  if (bad) ++failures;
}

const char* bc_name(int bc) { return bc == 0 ? "periodic" : bc == 1 ? "specified" : bc == 2 ? "open" : "symmetric"; }

// One stage from fresh random fields: WRF (bridge) and ScalarTransport::run_stage side by side
void stage(int bc, int rk, int adv_opt, int h_order, int v_order) {
  Grid g;
  g.nx = NX; g.ny = NY; g.nz = NZ; g.dt = DT; g.dx = 3000; g.dy = 2500;
  const BC kinds[] = {BC::Periodic, BC::Specified, BC::Open, BC::Symmetric};
  g.bc_x = g.bc_y = kinds[bc];
  State s(g);
  Tendencies t(g);
  BaseState b(g, false, true);  // horizontally varying mub and map factors
  s.add_scalar(g, "QVAPOR", true);
  t.add_scalar(g, "QVAPOR");

  std::vector<double> znw(NZ + 1);
  for (int k = 0; k <= NZ; ++k) znw[k] = k == 0 ? 1.0 : k == NZ ? 0.0 : 1.0 - k / double(NZ) + rnd(-0.03, 0.03);
  b.set_levels(znw);
  for (int k = 0; k <= NZ; ++k) {  // hybrid coordinate: c1, c2 not 1, 0 (mix-ups show)
    if (k < NZ) { b.c1h[k] = rnd(0.3, 1.0); b.c2h[k] = rnd(0.0, 5000.0); }
    b.c1f[k] = rnd(0.3, 1.0); b.c2f[k] = rnd(0.0, 5000.0);
  }
  b.upload_levels();
  b.rdx = 1.0 / g.dx; b.rdy = 1.0 / g.dy;
  // static 2-D fields with halos as a run has them: periodic copies, mirror images (symmetric) or
  // the outermost row copied outward (specified, open); n = interior points, st = staggered
  auto fold = [&](int i, int n, bool st) {
    const int m = st ? n - 1 : n;  // mass-point count of this direction
    if (bc == 0) return ((i % m) + m) % m;
    if (bc == 3) return i < 0 ? (st ? -i : -1 - i) : i >= n ? (st ? 2 * m - i : 2 * n - 1 - i) : i;
    return i < 0 ? 0 : i >= n ? n - 1 : i;
  };
  auto field2 = [&](BaseField2& m, bool sx, bool sy, double lo, double hi) {
    const int ni = NX + (sx ? 1 : 0), nj = NY + (sy ? 1 : 0);
    std::vector<double> v(std::size_t(ni) * nj);
    for (double& x : v) x = rnd(lo, hi);
    m.fill_host_halo([&](int i, int j, int) { return v[fold(i, ni, sx) + std::size_t(ni) * fold(j, nj, sy)]; });
  };
  field2(b.msftx, false, false, 0.9, 1.1); field2(b.msfty, false, false, 0.9, 1.1);
  field2(b.msfux, true, false, 0.9, 1.1); field2(b.msfuy, true, false, 0.9, 1.1);
  field2(b.msfvx, false, true, 0.9, 1.1); field2(b.msfvy, false, true, 0.9, 1.1);
  field2(b.mub, false, false, 6e4, 9e4);

  // mass fluxes ru_m, rv_m, ww_m (coupled, ~ mu u / msf), stage masses, the scalar and its tendency
  fill(s.ru_m, -8e5, 8e5); fill(s.rv_m, -8e5, 8e5); fill(s.ww_m, -2e2, 2e2);
  fill(s.muts, 6e4, 9e4); fill(s.mu_1, -2e3, 2e3); fill(s.mu_2, -2e3, 2e3);
  fill(s.q_1[0], 0.0, 2e-2); fill(s.q_2[0], 0.0, 2e-2);
  // a few near-empty cells, so the positive-definite limiter acts
  s.q_2[0].fill_host([&](int, int, int) { const double x = rnd(0.0, 2e-2); return x < 4e-3 ? x * 1e-4 : x; });
  fill(t.q_tendf[0], -3e-2, 3e-2);
  if (bc == 3) {  // rigid walls: no mass flux through them
    s.ru_m.fill_host([&](int i, int j, int k) { return i == 0 || i == NX ? 0.0 : rnd(-8e5, 8e5); });
    s.rv_m.fill_host([&](int i, int j, int k) { return j == 0 || j == NY ? 0.0 : rnd(-8e5, 8e5); });
  }
  fill_halos({&s.ru_m, &s.rv_m, &s.ww_m}, g);
  fill_halos({&s.muts, &s.mu_1, &s.mu_2}, g);
  fill_halos({&s.q_1[0], &s.q_2[0]}, g);

  const double dt_rk = rk == 1 ? DT / 3 : rk == 2 ? DT / 2 : DT;
  std::vector<double> q1 = wrf(s.q_1[0]), q2 = wrf(s.q_2[0]), qt = wrf(t.q_tendf[0]);
  scalar_stage_c(NX, NY, NZ, bc, h_order, v_order, adv_opt, rk, dt_rk, q1.data(), q2.data(), qt.data(),
                 wrf(s.ru_m).data(), wrf(s.rv_m).data(), wrf(s.ww_m).data(), wrf(s.muts).data(), wrf(b.mub).data(),
                 wrf(s.mu_1).data(), wrf(s.mu_2).data(), wrf(b.c1h).data(), wrf(b.c2h).data(), wrf(b.c1f).data(),
                 wrf(b.c2f).data(), wrf(b.fnm).data(), wrf(b.fnp).data(), wrf(b.rdnw).data(), wrf(b.msfux).data(),
                 wrf(b.msfuy).data(), wrf(b.msfvx).data(), wrf(b.msfvy).data(), wrf(b.msftx).data(),
                 wrf(b.msfty).data(), b.rdx, b.rdy);

  ScalarParams p;
  p.moist_adv_opt = adv_opt; p.h_adv_order = h_order; p.v_adv_order = v_order;
  ScalarTransport st(g, p);
  if (rk == 1) s.q_1[0].copy_from(s.q_2[0]);  // rk3.cpp at the start of the step (rk_update_scalar, rk_step 1)
  st.run_stage(s, t, b, rk, real(dt_rk));

  const std::string tag = std::string(bc_name(bc)) + " rk" + std::to_string(rk) + (adv_opt ? " pd" : " orig") +
                          " " + std::to_string(h_order) + "/" + std::to_string(v_order);
  check(tag + " q_2", s.q_2[0], q2);
  if (rk == 3 && adv_opt == 1) check(tag + " q_1", s.q_1[0], q1);
}

}  // namespace

int main() {
  std::printf("ScalarTransport vs WRF v4.8.0 (nx %d, ny %d, nz %d, dt %g s)\n", NX, NY, NZ, DT);
  for (int bc = 0; bc < 4; ++bc) {
    stage(bc, 1, 1, 5, 3);
    stage(bc, 2, 1, 5, 3);
    stage(bc, 3, 1, 5, 3);
    stage(bc, 3, 0, 5, 3);
    stage(bc, 3, 1, 6, 4);
    stage(bc, 2, 0, 6, 4);
  }
  std::printf(failures ? "FAIL: %d comparisons differ\n" : "PASS: all comparisons bit-identical\n", failures);
  return failures ? 1 : 0;
}
