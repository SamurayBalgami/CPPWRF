// Differential test of src/dyn/physics_coupling.cpp (PhysicsCoupling::step, mp_physics = 1) against
// WRF v4.8.0's moist_physics_prep_em, kessler and moist_physics_finish_em (mp_bridge.f90), from the
// same random but physically plausible columns (near-saturated, with cloud and rain), compared bit
// for bit over the whole domain: the spec zone (left to spec_bdy_final) and the halo (tile_bc_loop_2).
// Periodic, specified, open and symmetric boundaries; the grid is not square.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "core/constants.hpp"
#include "core/field.hpp"
#include "core/grid.hpp"
#include "core/halo.hpp"
#include "dyn/physics_coupling.hpp"
#include "dyn/state.hpp"

using namespace cppwrf;
using namespace cppwrf::dyn;

extern "C" void mp_step_c(int nx, int ny, int nz, int bc, double dt, double mp_tend_lim, double* t2, double* t1,
                          double* qv, double* qc, double* qr, double* h_diabatic, double* qv_diabatic,
                          double* qc_diabatic, const double* al, const double* alb, const double* p, const double* pb,
                          const double* ph, const double* phb, const double* fnm, const double* fnp, double* rainnc,
                          double* rainncv, double* th_phy, double* pii, double* rho, double* z, double* dz8w);

namespace {

constexpr int NX = 11, NY = 9, NZ = 12, H = 4;  // H: WRF's serial memory halo (bdyzone)
constexpr int NI = NX + 1 + 2 * H, NJ = NY + 1 + 2 * H, NK = NZ + 1;
constexpr double DT = 12.0, LIM = 10.0;

std::size_t f3(int i, int j, int k) { return std::size_t(i + H) + NI * (std::size_t(k) + NK * std::size_t(j + H)); }
std::size_t f2(int i, int j) { return std::size_t(i + H) + NI * std::size_t(j + H); }

std::uint64_t rng = 0x9E3779B97F4A7C15ull;
double rnd(double lo, double hi) {
  rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
  return lo + (hi - lo) * double(rng >> 11) / double(1ull << 53);
}

// WRF memory array of a field (WRF's fourth halo row, which CPPWRF does not store, is 0)
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
template <int R> std::vector<double> wrf(const BaseFieldT<R>& f) {
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
std::vector<double> wrf(const Array1& c) {
  std::vector<double> a(NK, 0.0);
  for (int k = 0; k < c.size(); ++k) a[k] = c[k];
  return a;
}

int failures = 0;
bool same(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }
// CPPWRF's field against WRF's array over i0..i1-1, j0..j1-1 (halo: the halo too)
template <int R>
void check(const std::string& what, const FieldT<R>& f, const std::vector<double>& ref, int i0, int i1, int j0, int j1) {
  std::vector<real> h;
  f.download(h);
  const Layout& L = f.layout();
  long bad = 0, n = 0;
  for (int k = 0; k < L.nz; ++k)
    for (int j = j0; j < j1; ++j)
      for (int i = i0; i < i1; ++i, ++n) {
        const double c = h[L.index(i, j, k)], w = ref[R == 3 ? f3(i, j, k) : f2(i, j)];
        if (same(c, w)) continue;
        if (++bad <= 3) std::printf("    %s(i=%d, j=%d, k=%d): CPPWRF %a  WRF %a\n", what.c_str(), i, j, k, c, w);
      }
  std::printf("  %-40s %6ld points  %s\n", what.c_str(), n, bad ? (std::to_string(bad) + " differ  FAIL").c_str() : "identical");
  if (bad) ++failures;
}

const char* bc_name(int bc) { return bc == 0 ? "periodic" : bc == 1 ? "specified" : bc == 2 ? "open" : "symmetric"; }

void run(int bc) {
  Grid g;
  g.nx = NX; g.ny = NY; g.nz = NZ; g.dx = 3000; g.dy = 3000; g.dt = DT;
  const BC kinds[] = {BC::Periodic, BC::Specified, BC::Open, BC::Symmetric};
  g.bc_x = g.bc_y = kinds[bc];
  State s(g);
  BaseState b(g, false);
  const int qv = s.add_scalar(g, "qv", true), qc = s.add_scalar(g, "qc", true), qr = s.add_scalar(g, "qr", true);
  std::vector<double> znw(NZ + 1);
  for (int k = 0; k <= NZ; ++k) znw[k] = 1.0 - double(k) / NZ;
  b.set_levels(znw);
  b.upload_levels();

  // columns: heights of w levels, base pressure falling with height, theta rising, qv near saturation
  std::vector<double> zw(std::size_t(NX + 2 * Grid::halo) * (NY + 2 * Grid::halo) * (NZ + 1));
  auto zat = [&](int i, int j, int k) -> double& {
    return zw[std::size_t(i + Grid::halo) + (NX + 2 * Grid::halo) * (std::size_t(j + Grid::halo) + (NY + 2 * Grid::halo) * k)];
  };
  for (int j = -Grid::halo; j < NY + Grid::halo; ++j)
    for (int i = -Grid::halo; i < NX + Grid::halo; ++i) {
      zat(i, j, 0) = rnd(0, 300);
      for (int k = 1; k <= NZ; ++k) zat(i, j, k) = zat(i, j, k - 1) + rnd(400, 1200);
    }
  b.phb.fill_host_halo([&](int i, int j, int k) { return cst::g * zat(i, j, k); });
  auto zm = [&](int i, int j, int k) { return 0.5 * (zat(i, j, k) + zat(i, j, k + 1)); };
  b.pb.fill_host_halo([&](int i, int j, int k) { return 1e5 * std::exp(-zm(i, j, k) / 8000.0); });
  b.alb.fill_host_halo([&](int i, int j, int k) { return 287.0 * 285.0 / (1e5 * std::exp(-zm(i, j, k) / 8000.0)); });
  s.ph_2.fill_host([](int, int, int) { return rnd(-20, 20); });
  s.p.fill_host([](int, int, int) { return rnd(-80, 80); });
  s.al.fill_host([](int, int, int) { return rnd(-1e-3, 1e-3); });
  s.q_2[qv].fill_host([&](int i, int j, int k) { return 0.018 * std::exp(-zm(i, j, k) / 2500.0) * rnd(0.6, 1.25); });
  s.q_2[qc].fill_host([&](int, int, int k) { return k < NZ - 2 && rnd(0, 1) < 0.5 ? rnd(0, 1.5e-3) : 0.0; });
  s.q_2[qr].fill_host([&](int, int, int k) { return k < NZ - 2 && rnd(0, 1) < 0.5 ? rnd(0, 1.5e-3) : 0.0; });
  s.t_2.fill_host([&](int i, int j, int k) {  // theta_m - t0
    const double th = 295.0 + 3.5e-3 * zm(i, j, k) + rnd(-1, 1);
    return th * (1.0 + cst::rvovrd * 0.018 * std::exp(-zm(i, j, k) / 2500.0)) - cst::t0;
  });
  fill_halos({&s.t_2, &s.ph_2, &s.p, &s.al, &s.q_2[qv], &s.q_2[qc], &s.q_2[qr]}, g);
  s.h_diabatic.fill_host([](int, int, int) { return 0.0; });  // as at the start of a run

  PhysicsParams pp;
  pp.mp_physics = 1;
  pp.mp_tend_lim = real(LIM);
  PhysicsCoupling pc(g, pp);

  std::vector<double> t2 = wrf(s.t_2), t1 = t2, wqv = wrf(s.q_2[qv]), wqc = wrf(s.q_2[qc]), wqr = wrf(s.q_2[qr]),
                      hd = wrf(s.h_diabatic), nul(std::size_t(NI) * NK * NJ, 0.0), qvd = nul, qcd = nul,
                      th_phy = nul, pii = nul, rho = nul, z = nul, dz8w = nul, rainnc(std::size_t(NI) * NJ, 0.0),
                      rainncv = rainnc;
  mp_step_c(NX, NY, NZ, bc, DT, LIM, t2.data(), t1.data(), wqv.data(), wqc.data(), wqr.data(), hd.data(), qvd.data(),
            qcd.data(), wrf(s.al).data(), wrf(b.alb).data(), wrf(s.p).data(), wrf(b.pb).data(), wrf(s.ph_2).data(),
            wrf(b.phb).data(), wrf(b.fnm).data(), wrf(b.fnp).data(), rainnc.data(), rainncv.data(), th_phy.data(),
            pii.data(), rho.data(), z.data(), dz8w.data());
  pc.step(s, b, real(DT));

  const int sz = bc == 1 ? 1 : 0, hx = Grid::halo;
  std::printf("PhysicsCoupling (Kessler) vs WRF v4.8.0, %s boundaries\n", bc_name(bc));
  // the state: whole domain (WRF leaves the spec zone to spec_bdy_final)
  check("t_2 (with halo)", s.t_2, t2, -hx, NX + hx, -hx, NY + hx);
  check("qv (with halo)", s.q_2[qv], wqv, -hx, NX + hx, -hx, NY + hx);
  check("qc (with halo)", s.q_2[qc], wqc, -hx, NX + hx, -hx, NY + hx);
  check("qr (with halo)", s.q_2[qr], wqr, -hx, NX + hx, -hx, NY + hx);
  check("h_diabatic (with halo)", s.h_diabatic, hd, -hx, NX + hx, -hx, NY + hx);
  check("rainnc (mp columns)", pc.rainnc, rainnc, sz, NX - sz, sz, NY - sz);
  check("rainncv (mp columns)", pc.rainncv, rainncv, sz, NX - sz, sz, NY - sz);
  const std::vector<const Field3*> fs = pc.fields();  // rho, qv_diabatic, qc_diabatic
  check("rho (mp columns)", *fs[0], rho, sz, NX - sz, sz, NY - sz);
  check("qv_diabatic", *fs[1], qvd, 0, NX, 0, NY);
  check("qc_diabatic", *fs[2], qcd, 0, NX, 0, NY);
}

}  // namespace

int main() {
  for (int bc = 0; bc < 4; ++bc) run(bc);
  std::printf(failures ? "FAIL: %d comparisons differ\n" : "PASS: all comparisons bit-identical\n", failures);
  return failures ? 1 : 0;
}
