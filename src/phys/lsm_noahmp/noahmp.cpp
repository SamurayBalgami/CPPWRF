#include "phys/lsm_noahmp/noahmp.hpp"

namespace cppwrf::phys::noahmp {

void run_noahmp_init(NoahmpInitKernel kernel, const Range2D& range) {
  assert(kernel.nsoil <= kMaxSoilLayers);
  assert(kernel.nsnow <= kMaxSnowLayers);
  column_for(range, kernel);
}

void run_noahmp(NoahmpKernel kernel, const Range2D& range) {
  assert(kernel.nsoil <= kMaxSoilLayers);
  assert(kernel.nsnow <= kMaxSnowLayers);
  column_for(range, kernel);
  fatal_check();  // Noah-MP's stop conditions (the device flag of this translation unit)
}

}  // namespace cppwrf::phys::noahmp
