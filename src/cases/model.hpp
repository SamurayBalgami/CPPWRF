#pragma once
// Model assembly from a config (host side): grid, case, vertical levels, base and initial
// state, scalars, numerical parameters and the RK3 integrator. Shared by main and tests.
#include <memory>
#include <string>
#include <vector>

#include "cases/cases.hpp"
#include "dyn/lateral_bc.hpp"
#include "dyn/rk3.hpp"
#include "io/config.hpp"

namespace cppwrf::cases {

struct Model {
  Grid g;
  CaseDef cd;
  std::vector<double> eta;
  std::unique_ptr<dyn::BaseState> base;
  std::unique_ptr<dyn::State> st;
  std::unique_ptr<dyn::Tendencies> tend;
  dyn::LargeStepParams lp;
  dyn::AcousticParams ap;
  dyn::ScalarParams sp;
  double run_time = 0;
  long nsteps = 0;
  int n_sound = 0;
  std::unique_ptr<dyn::RK3> rk;  // null when nsteps == 0
  std::unique_ptr<dyn::LateralBC> bc;  // real data: specified lateral boundaries
  static Model build(const io::Config& cfg);
};

// case = real: grid, base and initial state, boundaries from wrfinput/wrfbdy (src/cases/real_case.cpp)
void init_real_case(const io::Config& cfg, Model& m);

// Comma-separated list (whitespace trimmed, empty entries dropped).
std::vector<std::string> split_list(const std::string& s);

}  // namespace cppwrf::cases
