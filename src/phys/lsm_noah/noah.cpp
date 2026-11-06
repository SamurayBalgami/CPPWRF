#include "phys/lsm_noah/noah.hpp"
#include "core/limits.hpp"

namespace cppwrf::phys::noah {

void run_noah(NoahKernel kernel, const Range2D& range) {
  check_limit("Noah", "num_soil_layers", kernel.nsoil, kMaxSoilLevels);
  column_for(range, kernel);
  fatal_check();
}

}  // namespace cppwrf::phys::noah
