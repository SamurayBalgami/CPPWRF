#pragma once

#include <vector>

#include "io/config.hpp"
#include "core/grid.hpp"
#include "dyn/state.hpp"

namespace cppwrf::cases {

// Initializes 3D supercell case (em_quarter_ss, hybrid_opt = 0) from input_sounding.
void init_quarter_ss(const io::Config& cfg, const std::vector<double>& eta, const Grid& g,
                     dyn::BaseState& b, dyn::State& s);

}  // namespace cppwrf::cases
