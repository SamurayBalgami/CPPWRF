#pragma once

#include <vector>

#include "core/grid.hpp"
#include "dyn/state.hpp"
#include "io/config.hpp"

namespace cppwrf::cases {

// Initializes 2D squall line case matching WRF em_squall2d_x from input_sounding.
void init_squall2d_x(const io::Config& cfg, const std::vector<double>& eta, const Grid& g,
                     dyn::BaseState& b, dyn::State& s);

}  // namespace cppwrf::cases
