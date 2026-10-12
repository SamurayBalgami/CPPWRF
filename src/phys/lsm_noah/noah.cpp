#include "phys/lsm_noah/noah.hpp"

namespace cppwrf::phys::noah {

void run_noah(NoahKernel kernel, const Range2D& range) {
  assert(kernel.nsoil <= kMaxSoilLevels);
  column_for(range, kernel);
}

}  // namespace cppwrf::phys::noah
