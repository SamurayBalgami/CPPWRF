// Acoustic step count of solve_em.F (time_step_sound = 0) against WRF's own arithmetic: the expected
// values are what gfortran gives for solve_em.F:462 with -fdefault-real-8 (FP64 builds) and with the
// default REAL(4) (FP32 builds). 300*dt/dx - 0.01 is exactly an integer for dx 1000 m, dt 6.7 s, so
// the two precisions round to different sides of it.
#include <cstdio>

#include "core/precision.hpp"
#include "dyn/rk3.hpp"

using namespace cppwrf;
static int fails = 0;
#define CHECK(c, ...)                                                                                       \
  do {                                                                                                      \
    if (!(c)) { std::printf("FAIL %s:%d %s ", __FILE__, __LINE__, #c); std::printf(__VA_ARGS__); std::printf("\n"); ++fails; } \
  } while (0)

int main() {
  const bool fp64 = sizeof(real) == 8;
  struct Case { double dt, dx, dy; int r8, r4; };
  const Case cases[] = {
      {6.7, 1000, 1000, 4, 6},  {20.1, 3000, 3000, 4, 6}, {12.0, 1000, 1000, 8, 8},
      {6.0, 1000, 1000, 4, 4},  {6.7, 1000, 4000, 4, 6},  {20.1, 5000, 3000, 4, 6},  // min(dx, dy)
      {1.0, 10000, 10000, 4, 4}, {60.0, 9000, 9000, 4, 4}, {90.0, 9000, 9000, 6, 6},
  };
  for (const Case& c : cases) {
    const int n = dyn::RK3::auto_sound_steps(c.dt, c.dx, c.dy), want = fp64 ? c.r8 : c.r4;
    CHECK(n == want, "dt %g dx %g dy %g: %d sound steps, WRF %d", c.dt, c.dx, c.dy, n, want);
  }
  std::printf(fails ? "test_sound_steps: %d FAILED\n" : "test_sound_steps: all passed\n", fails);
  return fails ? 1 : 0;
}
