// Diagnostics::compute (diag.csv row) against a plain host loop over the definitions, from random
// fields on a non-square grid with hybrid c1/c2: maxima exact, sums to reduction-order roundoff.
//  (a) dry, (b) moist: t_2 holds theta_m - t0, the energy uses the dry theta (theta_mass keeps
//  theta_m, the quantity flux-form advection conserves), (c) a NaN shows as a non-finite row.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "core/constants.hpp"
#include "core/precision.hpp"
#include "dyn/diagnostics.hpp"
#include "dyn/state.hpp"

using namespace cppwrf;
using namespace cppwrf::dyn;
static int fails = 0;
#define CHECK(c, ...)                                                                                       \
  do {                                                                                                      \
    if (!(c)) { std::printf("FAIL %s:%d %s ", __FILE__, __LINE__, #c); std::printf(__VA_ARGS__); std::printf("\n"); ++fails; } \
  } while (0)

static const double tol = sizeof(real) == 8 ? 1e-12 : 1e-5;  // sums: reduction order differs
static std::uint64_t rng = 0x2545F4914F6CDD1Dull;
static double rnd(double lo, double hi) {
  rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
  return lo + (hi - lo) * double(rng >> 11) / double(1ull << 53);
}
static bool close(double a, double b) { return std::fabs(a - b) <= tol * std::fabs(b); }

struct Host {  // downloaded field with its layout
  std::vector<real> h; Layout L;
  template <class F> explicit Host(const F& f) : L(f.layout()) { f.download(h); }
  double operator()(int i, int j, int k = 0) const { return double(h[L.index(i, j, k)]); }
};

static void run(bool moist) {
  Grid g;
  g.nx = 7; g.ny = 5; g.nz = 6; g.dx = 2000; g.dy = 1500;
  State s(g);
  BaseState b(g, false);
  if (moist) s.add_scalar(g, "qv", true);
  std::vector<double> znw = {1.0, 0.87, 0.71, 0.52, 0.33, 0.15, 0.0};
  b.set_levels(znw);
  for (int k = 0; k < g.nz; ++k) { b.c1h[k] = rnd(0.4, 1.0); b.c2h[k] = rnd(0.0, 4000.0); }
  b.upload_levels();
  auto fill3 = [](Field3& f, double lo, double hi) { f.fill_host([&](int, int, int) { return rnd(lo, hi); }); };
  fill3(s.u_2, -30, 30); fill3(s.v_2, -30, 30); fill3(s.w_2, -5, 5); fill3(s.t_2, -10, 40);
  fill3(s.p, -500, 500); fill3(s.ph_2, -50, 50);
  s.mu_2.fill_host([](int, int, int) { return rnd(-800, 800); });
  b.mub.fill_host_halo([](int, int, int) { return rnd(8e4, 9.5e4); });
  b.pb.fill_host_halo([](int, int, int k) { return 9.5e4 - 1.4e4 * k + rnd(-100, 100); });
  b.phb.fill_host_halo([](int, int, int k) { return 9.81 * (1500.0 * k + rnd(0, 50)); });
  b.t_init.fill_host_halo([](int, int, int) { return rnd(0, 30); });
  if (moist) s.q_2[0].fill_host([](int, int, int) { return rnd(0, 0.02); });

  Diagnostics d(g);
  const DiagRow r = d.compute(s, b);

  const Host u(s.u_2), v(s.v_2), w(s.w_2), t(s.t_2), p(s.p), ph(s.ph_2), mu(s.mu_2), mub(b.mub), pb(b.pb),
      phb(b.phb), ti(b.t_init);
  const Host qv = moist ? Host(s.q_2[0]) : Host(s.t_2);
  double mxu = 0, mxv = 0, mxw = 0, thmax = -1e300, thmin = 1e300, mass = 0, thmass = 0, energy = 0;
  for (int k = 0; k < g.nz; ++k)
    for (int j = 0; j <= g.ny; ++j)
      for (int i = 0; i <= g.nx; ++i) {
        if (j < g.ny) mxu = std::fmax(mxu, std::fabs(u(i, j, k)));
        if (i < g.nx) mxv = std::fmax(mxv, std::fabs(v(i, j, k)));
      }
  for (int k = 0; k <= g.nz; ++k)
    for (int j = 0; j < g.ny; ++j)
      for (int i = 0; i < g.nx; ++i) mxw = std::fmax(mxw, std::fabs(w(i, j, k)));
  for (int j = 0; j < g.ny; ++j)
    for (int i = 0; i < g.nx; ++i) {
      const double mfull = mu(i, j) + mub(i, j);
      mass += mfull;
      for (int k = 0; k < g.nz; ++k) {
        const double dp = -double(b.dnw[k]) * (double(b.c1h[k]) * mfull + double(b.c2h[k]));  // layer dry mass * g
        const double thm = t(i, j, k) + cst::t0;
        const double th = moist ? thm / (1.0 + cst::rvovrd * qv(i, j, k)) : thm;  // dry theta
        thmax = std::fmax(thmax, th - cst::t0 - ti(i, j, k));
        thmin = std::fmin(thmin, th - cst::t0 - ti(i, j, k));
        thmass += dp * thm;
        const double uc = 0.5 * (u(i, j, k) + u(i + 1, j, k)), vc = 0.5 * (v(i, j, k) + v(i, j + 1, k)),
                     wc = 0.5 * (w(i, j, k) + w(i, j, k + 1));
        const double T = th * std::pow((p(i, j, k) + pb(i, j, k)) / cst::p1000mb, cst::rcp);
        const double gz = 0.5 * (ph(i, j, k) + ph(i, j, k + 1) + phb(i, j, k) + phb(i, j, k + 1));
        energy += dp * (0.5 * (uc * uc + vc * vc + wc * wc) + cst::cv * T + gz);
      }
    }
  const double area = g.dx * g.dy / cst::g;
  const char* tag = moist ? "moist" : "dry";
  CHECK(r.max_abs_u == mxu, "%s max|u| %.17g ref %.17g", tag, r.max_abs_u, mxu);
  CHECK(r.max_abs_v == mxv, "%s max|v| %.17g ref %.17g", tag, r.max_abs_v, mxv);
  CHECK(r.max_abs_w == mxw, "%s max|w| %.17g ref %.17g", tag, r.max_abs_w, mxw);
  CHECK(close(r.max_th_pert, thmax), "%s max th' %.17g ref %.17g", tag, r.max_th_pert, thmax);
  CHECK(close(r.min_th_pert, thmin), "%s min th' %.17g ref %.17g", tag, r.min_th_pert, thmin);
  CHECK(close(r.dry_mass, area * mass), "%s dry mass %.17g ref %.17g", tag, r.dry_mass, area * mass);
  CHECK(close(r.theta_mass, area * thmass), "%s theta mass %.17g ref %.17g", tag, r.theta_mass, area * thmass);
  CHECK(close(r.total_energy, area * energy), "%s energy %.17g ref %.17g", tag, r.total_energy, area * energy);
  CHECK(r.finite(), "%s row not finite", tag);

  // a NaN anywhere in a field makes the row non-finite
  s.w_2.fill_host([](int i, int j, int k) { return i == 3 && j == 2 && k == 4 ? std::nan("") : 0.0; });
  CHECK(!d.compute(s, b).finite(), "%s NaN in w not flagged", tag);
}

int main() {
  run(false);
  run(true);
  std::printf(fails ? "test_diagnostics: %d FAILED\n" : "test_diagnostics: all passed\n", fails);
  return fails ? 1 : 0;
}
