#pragma once
// Wicker-Skamarock RK3 with split-explicit acoustic substeps (TN 3.1-3.2; WRF solve_em).
//   stage 1: dt/3, 1 acoustic step of dt/3
//   stage 2: dt/2, n/2 acoustic steps of dts = dt/n
//   stage 3: dt,   n   acoustic steps of dts
#include <memory>

#include "dyn/large_step.hpp"
#include "dyn/state.hpp"

namespace cppwrf::dyn {

struct AcousticParams { real epssm = 0.1, smdiv = 0.1, emdiv = 0.01; };

class SmallStepHolder;  // wraps SmallStep (src/dyn/small_step.hpp) when present

class RK3 {
 public:
  RK3(const Grid& g, const LargeStepParams& lp, const AcousticParams& ap, int time_step_sound);
  ~RK3();
  // WRF solve_em auto rule when time_step_sound == 0 (n even, >= 4).
  static int auto_sound_steps(double dt, double dx, double dy, bool two_d);
  int n_sound() const { return n_sound_; }
  static bool small_step_available();
  void step(State& s, Tendencies& t, const BaseState& b);

 private:
  Grid g_;
  LargeStep large_;
  int n_sound_;
  std::unique_ptr<SmallStepHolder> small_;
};

}  // namespace cppwrf::dyn
