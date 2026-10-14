#pragma once

#include <cstddef>
#include <cstdint>
#include "core/backend.hpp"
#include "core/field.hpp"

namespace cppwrf::nest {

// Refinement ratio and nest boundary configurations
struct NestParams {
  int ratio_x = 3;            // parent to child grid refinement ratio in x (e.g. 3 or 5)
  int ratio_y = 3;            // parent to child grid refinement ratio in y (e.g. 3 or 5)
  int i_parent_start = 1;     // 1-based index in parent grid where nest starts in x
  int j_parent_start = 1;     // 1-based index in parent grid where nest starts in y
  int spec_bdy_width = 5;     // width of lateral boundary zone in nest
  int spec_zone = 1;          // specified zone width
  int relax_zone = 4;         // relaxation zone width
  int blend_width = 5;        // terrain blending transition zone width
  int smooth_option = 0;      // feedback smoothing: 0=none, 1=sm121, 2=smdsm
  int feedback = 1;           // two-way feedback: 0=one-way, 1=two-way
  double cdt = 1.0;           // coarse grid time step (s)
  double ndt = 1.0 / 3.0;     // fine grid time step (s)
  bool use_theta_m = true;    // whether moist potential temperature is used
  double p_top = 5000.0;      // top of atmosphere pressure (Pa)
};

// Side indicators matching lateral boundary conventions
enum class BdySide : int {
  West = 0,   // XS
  East = 1,   // XE
  South = 2,  // YS
  North = 3   // YE
};

}  // namespace cppwrf::nest
