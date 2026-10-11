// Large-step tests: (1) uniform flow U advecting a smooth theta' -> t_tend converges to
// -mu U d(theta')/dx at ~5th order; (2) balanced initial state -> rw_tend ~ 0, ww == 0;
// (3) tendencies are bitwise reproducible across repeated evaluations (OpenMP determinism).
#include <cmath>
#include <cstdio>
#include <memory>
#include <numbers>
#include <vector>

#include "cases/cases.hpp"
#include "core/constants.hpp"
#include "dyn/large_step.hpp"

using namespace cppwrf;
static int fails = 0;
#define CHECK(c, ...)                                         \
  do {                                                        \
    if (!(c)) { std::printf("FAIL %s:%d %s ", __FILE__, __LINE__, #c); std::printf(__VA_ARGS__); std::printf("\n"); ++fails; } \
  } while (0)

struct Setup {
  Grid g;
  cases::CaseDef cd;
  std::unique_ptr<cases::Sounding> snd;
  std::unique_ptr<dyn::BaseState> b;
  std::unique_ptr<dyn::State> s;
  std::unique_ptr<dyn::Tendencies> t;
  Setup(int nx, double U, double amp, bool stratified) {
    g.nx = nx; g.ny = 1; g.nz = 10; g.dx = 100000.0 / nx; g.dy = g.dx; g.ztop = 10000; g.dt = 1;
    io::Config cfg = io::Config::from_string("[case]\nu_mean = 0\n");
    cd = cases::make_case("igw", cfg, g);
    double L = g.nx * g.dx;
    if (!stratified) cd.theta_env = [](double) { return 300.0; };
    cd.u_env = [U](double) { return U; };
    cd.theta_pert = [amp, L](double x, double, double z, double) {
      return amp * std::sin(2 * std::numbers::pi * x / L) * std::sin(std::numbers::pi * z / 10000.0);
    };
    snd = std::make_unique<cases::Sounding>(cd.theta_env, cd.p_surf, 13000.0);
    g.p_top = snd->p(g.ztop);
    auto eta = cases::make_eta_levels(cfg, g, *snd);
    b = std::make_unique<dyn::BaseState>(g, true);
    s = std::make_unique<dyn::State>(g);
    t = std::make_unique<dyn::Tendencies>(g);
    cases::initialize(cd, *snd, eta, g, *b, *s);
  }
};

template <int R> static std::vector<real> host(const FieldT<R>& f) { std::vector<real> h; f.download(h); return h; }

static void test_advection_order() {
  // theta' small enough that the induced p'/ph' do not matter for t_tend (pure advection).
  double err[2];
  int nxs[2] = {50, 100};
  for (int r = 0; r < 2; ++r) {
    Setup S(nxs[r], 10.0, 1.0, false);
    dyn::LargeStep ls(S.g, {});
    ls.prep(*S.s, *S.b);
    ls.tendencies(*S.s, *S.t, *S.b, 1);
    auto tt = host(S.t->t_tend), mut = host(S.s->mut);
    const Layout &L = S.t->t_tend.layout(), &L2 = S.s->mut.layout();
    double L_x = S.g.nx * S.g.dx, e = 0;
    // Exact: t_tend = -U * mu * d(theta')/dx with theta'(x,k) = A_k sin(kx x); A_k from t_2.
    auto t2 = host(S.s->t_2);
    auto ti = std::vector<real>();
    S.b->t_init.download(ti);
    const Layout& Lt = S.s->t_2.layout();
    const Layout& Lti = S.b->t_init.layout();
    double kx = 2 * std::numbers::pi / L_x;
    for (int k = 0; k < S.g.nz; ++k) {
      // amplitude from a discrete Fourier projection of the initial field
      double a = 0;
      for (int i = 0; i < S.g.nx; ++i)
        a += (t2[Lt.index(i, 0, k)] - ti[Lti.index(i, 0, k)]) * std::sin(kx * (i + 0.5) * S.g.dx);
      a *= 2.0 / S.g.nx;
      for (int i = 0; i < S.g.nx; ++i) {
        double exact = -10.0 * mut[L2.index(i, 0, 0)] * a * kx * std::cos(kx * (i + 0.5) * S.g.dx);
        e = std::max(e, std::fabs(tt[L.index(i, 0, k)] - exact));
      }
    }
    err[r] = e;
  }
  double order = std::log2(err[0] / err[1]);
  std::printf("advection: err(nx=50)=%.3e err(nx=100)=%.3e order=%.2f\n", err[0], err[1], order);
  if (sizeof(real) == 8) CHECK(order > 4.5, "order %.2f", order);
  else CHECK(err[0] < 1e-3 && err[1] < 1e-3, "fp32 error above roundoff level");
}

static void test_balance_and_determinism() {
  Setup S(64, 0.0, 0.0, true);  // stratified, at rest: everything must vanish
  dyn::LargeStep ls(S.g, {});
  ls.prep(*S.s, *S.b);
  ls.tendencies(*S.s, *S.t, *S.b, 1);
  double mx = 0;
  for (auto* f : {&S.t->ru_tend, &S.t->rw_tend, &S.t->t_tend, &S.t->ph_tend, &S.s->ww})
    for (real v : host(*f)) mx = std::max(mx, double(std::fabs(v)));
  CHECK(mx == 0.0, "rest-state tendencies not zero: %g", mx);

  Setup B(64, 15.0, 3.0, true);  // non-trivial state for determinism
  dyn::LargeStepParams p;
  p.khdif = 50; p.kvdif = 50;
  dyn::LargeStep l2(B.g, p);
  l2.prep(*B.s, *B.b);
  l2.tendencies(*B.s, *B.t, *B.b, 1);
  auto a1 = host(B.t->ru_tend), b1 = host(B.t->rw_tend), c1 = host(B.t->t_tend), d1 = host(B.t->ph_tend);
  for (int rep = 0; rep < 3; ++rep) {
    l2.prep(*B.s, *B.b);
    l2.tendencies(*B.s, *B.t, *B.b, 1);
    CHECK(host(B.t->ru_tend) == a1 && host(B.t->rw_tend) == b1 && host(B.t->t_tend) == c1 && host(B.t->ph_tend) == d1,
          "rep %d differs", rep);
  }
  double m = 0;
  for (real v : c1) m = std::max(m, double(std::fabs(v)));
  CHECK(m > 0, "t_tend unexpectedly zero");
}

int main() {
  test_advection_order();
  test_balance_and_determinism();
  std::printf("test_large_step (%s, %s): %s\n", backend_name(), precision_name, fails ? "FAILED" : "OK");
  return fails ? 1 : 0;
}
