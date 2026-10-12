#pragma once

#include "real_types.hpp"

namespace cppwrf::init_real {

// Computes 1D eta levels using auto_levels_opt = 2 (SUBROUTINE levels in module_initialize_real.F)
// and derives znu, dnw, dn, rdnw, rdn.
void compute_eta_levels(WrfInputData& data,
                        float max_dz = 1000.0f,
                        float dzbot = 50.0f,
                        float dzstretch_s = 1.3f,
                        float dzstretch_u = 1.1f);

// Computes vertical coordinate 1D coefficients c1f..c4f and c1h..c4h
// for hybrid_opt = 0 (sigma) or hybrid_opt = 2 (Klemp polynomial).
void compute_vcoord_coeffs(WrfInputData& data);

} // namespace cppwrf::init_real
