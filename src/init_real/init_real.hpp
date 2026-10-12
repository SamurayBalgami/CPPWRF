#pragma once

#include "real_types.hpp"
#include "met_em_reader.hpp"
#include "vcoord.hpp"
#include "base_state.hpp"
#include "vert_interp.hpp"
#include "hydrostatic.hpp"
#include "bdy.hpp"
#include <string>
#include <vector>

namespace cppwrf::init_real {

struct RealInitResult {
    WrfInputData initial_state;
    WrfBdyData bdy_data;
};

// Top-level driver: takes met_em file paths (1 or more) and hybrid_opt (0 or 2),
// produces initial state (wrfinput) and lateral boundary data (wrfbdy).
RealInitResult run_real_init(const std::vector<std::string>& met_em_files,
                             int hybrid_opt = 2);

} // namespace cppwrf::init_real
