#include "kf.hpp"

namespace cppwrf::phys::cu_kf {

template <typename T>
void run_kf(KFKernel<T> kernel, const Range2D& range) {
    assert(kernel.nz <= kMaxLevels);
    column_for(range, kernel);
}

template void run_kf<double>(KFKernel<double>, const Range2D&);
template void run_kf<float>(KFKernel<float>, const Range2D&);

} // namespace cppwrf::phys::cu_kf
