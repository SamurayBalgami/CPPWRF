#pragma once

#include "real_types.hpp"
#include <vector>

namespace cppwrf::init_real {

// Packs lateral boundaries for a single WrfInputData state (coupled with mu and msf)
BoundaryInterval pack_single_boundary(const WrfInputData& state,
                                      const std::string& this_time,
                                      const std::string& next_time);

// Computes boundary tendencies for consecutive intervals across all 5 dates
WrfBdyData build_wrfbdy_data(const std::vector<WrfInputData>& states,
                             const std::vector<std::string>& date_strings,
                             float dt_interval = 21600.0f);

} // namespace cppwrf::init_real
