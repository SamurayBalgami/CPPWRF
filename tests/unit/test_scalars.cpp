// Scalar transport tests (M2 item 1):
//  (a) a tracer that is identically 1 stays 1 to roundoff (straka, warm_bubble, igw setups);
//  (b) with positive-definite advection (scalar_adv_opt = 1) a blob stays >= 0 (to roundoff),
//      without it the 5th-order scheme undershoots (checked in the strongly sheared straka flow);
//  (c) the dry-mass-weighted tracer mass sum((c1 mu + c2) q |dnw|) is conserved to roundoff
//      (periodic, no diffusion).
// Reduced grids and run lengths so the test is quick; the full cases are run by tools/.
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

struct QStats { double qmin, qmax, dev1, mass; };

static QStats stats(const cases::Model& M, int n) {
  std::vector<real> q, mu, mub;
  M.st->q_2[n].download(q);
  M.st->mu_2.download(mu);
  M.base->mub.download(mub);
  const Layout &L = M.st->q_2[n].layout(), &L2 = M.st->mu_2.layout(), &Lb = M.base->mub.layout();
  QStats s{1e300, -1e300, 0, 0};
  for (int k = 0; k < M.g.nz; ++k)
    for (int j = 0; j < M.g.ny; ++j)
      for (int i = 0; i < M.g.nx; ++i) {
        const double v = q[L.index(i, j, k)];
        s.qmin = std::min(s.qmin, v);
        s.qmax = std::max(s.qmax, v);
        s.dev1 = std::max(s.dev1, std::fabs(v - 1.0));
        const double m = M.base->c1h[k] * (double(mub[Lb.index(i, j, 0)]) + double(mu[L2.index(i, j, 0)])) + M.base->c2h[k];
        s.mass += m * v * (-M.base->dnw[k]);
      }
  return s;
}

static cases::Model run(const std::string& ini, int steps) {
  cases::Model M = cases::Model::build(io::Config::from_string(ini));
  for (int n = 0; n < steps; ++n) M.rk->step(*M.st, *M.tend, *M.base);
  dev_sync();
  return M;
}

static const char* kStraka =
    "[run]\ncase = straka\nrun_time = 400\n[domain]\nnx = 128\nnz = 16\ndx = 400\nztop = 6400\nbc_x = periodic\n"
    "[time]\ndt = 4\n[dyn]\nkhdif = 75\nkvdif = 75\npr_inv = 1\n";
static const char* kBubble =
    "[run]\ncase = warm_bubble\nrun_time = 400\n[domain]\nnx = 40\nnz = 40\ndx = 250\nztop = 10000\nbc_x = symmetric\n"
    "[time]\ndt = 2\n";
static const char* kIgw =
    "[run]\ncase = igw\nrun_time = 600\n[domain]\nnx = 100\nnz = 10\ndx = 3000\nztop = 10000\nbc_x = periodic\n"
    "[time]\ndt = 12\n[case]\nxc = 150000\n";

int main() {
  const double tol_one = fp64 ? 1e-12 : 1e-5;
  const double tol_mass = fp64 ? 1e-13 : 1e-5;
  struct C { const char* name; const char* ini; int steps; };
  for (C c : {C{"straka", kStraka, 100}, C{"warm_bubble", kBubble, 200}, C{"igw", kIgw, 50}}) {
    for (int pd = 0; pd <= 1; ++pd) {
      const bool igw = std::string(c.name) == "igw";
      std::string ini = std::string(c.ini) + "scalar_adv_opt = " + std::to_string(pd) +
                        "\n[tracers]\nnames = one, blob\nblob_init = blob\n" +
                        (igw ? "blob_xc = 150000\nblob_radius = 20000\nblob_zc = 5000\n"
                             : "blob_zc = 3000\nblob_radius = 1500\n");
      cases::Model M = run(ini, c.steps);
      QStats one = stats(M, 0), blob = stats(M, 1);
      cases::Model M0 = cases::Model::build(io::Config::from_string(ini));
      QStats blob0 = stats(M0, 1), one0 = stats(M0, 0);
      const double dm_one = (one.mass - one0.mass) / one0.mass, dm_blob = (blob.mass - blob0.mass) / blob0.mass;
      std::printf("%-11s pd=%d  max|one-1| %.2e  blob min %+.3e max %.4f  dmass one %+.2e blob %+.2e\n", c.name, pd,
                  one.dev1, blob.qmin, blob.qmax, dm_one, dm_blob);
      CHECK(one.dev1 < tol_one, "%s: constant tracer deviates by %.3e", c.name, one.dev1);
      if (pd == 1) CHECK(blob.qmin >= -(fp64 ? 1e-12 : 1e-6) * blob0.qmax, "%s: PD blob min %.3e", c.name, blob.qmin);
      const bool diff_v = std::string(c.name) == "straka";  // vertical diffusion: WRF form, not conservative
      if (!diff_v) CHECK(std::fabs(dm_blob) < tol_mass, "%s: tracer mass drift %.3e", c.name, dm_blob);
      if (pd == 0 && std::string(c.name) == "straka")
        CHECK(blob.qmin < -1e-4, "%s: expected undershoot without PD, min %.3e", c.name, blob.qmin);
    }
  }
  std::printf("test_scalars (%s, %s): %s\n", backend_name(), precision_name, fails ? "FAILED" : "OK");
  return fails ? 1 : 0;
}
