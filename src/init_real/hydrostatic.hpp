#pragma once

#include "real_types.hpp"
#include "met_em_reader.hpp"

namespace cppwrf::init_real {

// Performs hydrostatic integration to compute:
// - geopotential perturbation ph
// - perturbation specific volume al
// - perturbation pressure p
// - final qvapor from diagnosed pressure
void integrate_hydrostatic(const MetEmData& met,
                           WrfInputData& data,
                           const Array3D& rh_interp,
                           int hypsometric_opt = 2);

} // namespace cppwrf::init_real
