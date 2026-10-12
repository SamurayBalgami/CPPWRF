#pragma once

#include "real_types.hpp"
#include "met_em_reader.hpp"

namespace cppwrf::init_real {

// Converts relative humidity (%) to water vapor mixing ratio (kg/kg)
// using MM5 formula (rh2qv_method == 1)
float rh_to_mxrat_single(float rh, float t_k, float p_pa,
                         float qv_max_p_safe = 10000.0f,
                         float qv_max_flag = 1.0e-5f,
                         float qv_max_value = 3.0e-6f,
                         float qv_min_p_safe = 110000.0f,
                         float qv_min_flag = 1.0e-6f,
                         float qv_min_value = 1.0e-6f);

// Integrates column moisture to get dry pressure pd_gc and dry column mass mu0 / mu
void compute_dry_pressure_and_mu(const MetEmData& met,
                                 WrfInputData& data,
                                 Array3D& pd_gc,
                                 Array3D& qv_gc);

// Performs vertical interpolation of all fields (U, V, T, RH->QV, PH)
// from met_em isobaric levels to WRF eta levels in dry-pressure space.
void vertical_interpolate_fields(const MetEmData& met,
                                 const Array3D& pd_gc,
                                 const Array3D& qv_gc,
                                 WrfInputData& data,
                                 Array3D& rh_interp);

} // namespace cppwrf::init_real
