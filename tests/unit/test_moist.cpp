// Moist dynamics tests (M2 item 2, TN 2.2, use_theta_m = 1):
//  (a) a moist (qv0 = 10 g/kg, unsaturated, no microphysics) atmosphere at rest stays at rest:
//      max|w| at roundoff, of the same order as the dry rest state;
//  (b) with qv carried as a species but identically 0, the dynamics is bitwise identical to the
//      dry run (cqu = cqv = cq1 = 1, cq2 = 0 exactly; moist init reduces to the dry init).
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "cases/model.hpp"
#include "core/backend.hpp"
#include "core/precision.hpp"

using namespace cppwrf;
static int fails = 0;
#define CHECK(c, ...)                                                                                       \
  do {                                                                                                      \
    if (!(c)) { std::printf("FAIL %s:%d %s ", __FILE__, __LINE__, #c); std::printf(__VA_ARGS__); std::printf("\n"); ++fails; } \
  } while (0)

static const bool fp64 = sizeof(real) == 8;

static cases::Model run(const std::string& ini) {
  cases::Model M = cases::Model::build(io::Config::from_string(ini));
  for (int n = 0; n < M.nsteps; ++n) M.rk->step(*M.st, *M.tend, *M.base);
  dev_sync();
  return M;
}

// max |w| (m/s): W / (c1f mu_full + c2f) is approximated by W / (mub + mu) (exact for hybrid_opt = 0)
static double max_w(const cases::Model& M) {
  std::vector<real> w, mu, mub;
  M.st->w_2.download(w);
  M.st->mu_2.download(mu);
  M.base->mub.download(mub);
  const Layout &L = M.st->w_2.layout(), &L2 = M.st->mu_2.layout(), &Lb = M.base->mub.layout();
  double m = 0;
  for (int k = 0; k <= M.g.nz; ++k)
    for (int j = 0; j < M.g.ny; ++j)
      for (int i = 0; i < M.g.nx; ++i)
        m = std::max(m, std::fabs(double(w[L.index(i, j, k)])) /
                            (double(mub[Lb.index(i, j, 0)]) + double(mu[L2.index(i, j, 0)])));
  return m;
}

static bool same(const Field3& a, const Field3& b) {
  std::vector<real> x, y;
  a.download(x);
  b.download(y);
  return x == y;
}
static bool same(const Field2& a, const Field2& b) {
  std::vector<real> x, y;
  a.download(x);
  b.download(y);
  return x == y;
}

static const char* kRest =
    "[run]\ncase = rest_state\nrun_time = 1200\n[domain]\nnx = 24\nnz = 40\ndx = 2000\nztop = 20000\n"
    "bc_x = periodic\neta = uniform_z\n[time]\ndt = 10\n[dyn]\nkhdif = 0\nkvdif = 0\n";
static const char* kBubble =
    "[run]\ncase = warm_bubble\nrun_time = 200\n[domain]\nnx = 40\nnz = 40\ndx = 250\nztop = 10000\n"
    "bc_x = symmetric\n[time]\ndt = 2\n[dyn]\nkhdif = 10\nkvdif = 10\n";

int main() {
  // (a) moist rest state vs dry rest state
  const double w_dry = max_w(run(kRest));
  const double w_moist = max_w(run(std::string(kRest) + "[moist]\nspecies = qv\n[case]\nqv0 = 0.01\n"));
  std::printf("rest_state 1200 s: max|w| dry %.3e  moist (qv0 = 10 g/kg) %.3e m/s\n", w_dry, w_moist);
  const double tol_w = fp64 ? 1e-14 : 1e-8;
  CHECK(w_moist < tol_w, "moist rest state max|w| %.3e", w_moist);

  // (b) qv species identically zero: bitwise identical dynamics
  cases::Model D = run(kBubble);
  cases::Model Q = run(std::string(kBubble) + "[moist]\nspecies = qv\n");
  const bool eq = same(D.st->u_2, Q.st->u_2) && same(D.st->w_2, Q.st->w_2) && same(D.st->t_2, Q.st->t_2) &&
                  same(D.st->ph_2, Q.st->ph_2) && same(D.st->mu_2, Q.st->mu_2) && same(D.st->p, Q.st->p);
  std::printf("warm_bubble 200 s, qv = 0 species vs dry: %s\n", eq ? "bitwise identical" : "DIFFER");
  CHECK(eq, "qv = 0 run differs from the dry run");

  std::printf("test_moist (%s, %s): %s\n", backend_name(), precision_name, fails ? "FAILED" : "OK");
  return fails ? 1 : 0;
}
