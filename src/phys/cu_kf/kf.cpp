#include "kf.hpp"

namespace cppwrf::phys::cu_kf {

template <typename T>
void run_kf(KFKernel<T> kernel, const Range2D& range) {
    assert(kernel.nz <= kMaxLevels);
    if (kernel.nz <= 48) {
        column_for(range, KFRunner<T, 48>{kernel});
    } else {
        column_for(range, KFRunner<T, kMaxLevels>{kernel});
    }
}

template void run_kf<double>(KFKernel<double>, const Range2D&);
template void run_kf<float>(KFKernel<float>, const Range2D&);

} // namespace cppwrf::phys::cu_kf
