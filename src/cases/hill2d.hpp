#pragma once

#include <vector>

#include "core/grid.hpp"
#include "dyn/state.hpp"
#include "io/config.hpp"

namespace cppwrf::cases {

// Initializes 2D hill terrain (case hill2d_x, hybrid_opt = 0) from WRF input_sounding.
// Sets non-flat BaseState (ht, mub, pb, alb, t_init, phb) and dry prognostic State.
void init_hill2d(const io::Config& cfg, const std::vector<double>& eta, const Grid& g,
                 dyn::BaseState& b, dyn::State& s);

}  // namespace cppwrf::cases
