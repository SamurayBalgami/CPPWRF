// SPDX-License-Identifier: GPL-3.0-or-later
#include "phys/lsm_noahmp/noahmp.hpp"
#include "core/limits.hpp"

namespace cppwrf::phys::noahmp {

void run_noahmp_init(NoahmpInitKernel kernel, const Range2D& range) {
  check_limit("Noah-MP", "num_soil_layers", kernel.nsoil, kMaxSoilLayers);
  check_limit("Noah-MP", "snow layers", kernel.nsnow, kMaxSnowLayers);
  column_for(range, kernel);
}

void run_noahmp(NoahmpKernel kernel, const Range2D& range) {
  check_limit("Noah-MP", "num_soil_layers", kernel.nsoil, kMaxSoilLayers);
  check_limit("Noah-MP", "snow layers", kernel.nsnow, kMaxSnowLayers);
  column_for(range, kernel);
  fatal_check();  // Noah-MP's stop conditions (the device flag of this translation unit)
}

}  // namespace cppwrf::phys::noahmp
