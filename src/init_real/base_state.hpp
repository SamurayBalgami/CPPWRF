#pragma once

#include "real_types.hpp"

namespace cppwrf::init_real {

// Computes WRF ideal hydrostatic base state fields:
// mub, pb, t_init, alb, phb
void compute_base_state(WrfInputData& data, int hypsometric_opt = 2);

} // namespace cppwrf::init_real
