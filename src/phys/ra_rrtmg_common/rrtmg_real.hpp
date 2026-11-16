// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// RRTMG's arithmetic type (src/phys/ra_rrtmg_lw, ra_rrtmg_sw, ra_rrtmg_common, the radiation driver, Dudhia SW).
// CUDA FP32 builds: float, as RRTMG's rb in WRF real*4; the FP64 rate of most GPUs is 1/32 to 1/64 of the FP32 rate,
// and RRTMG in double dominated the FP32 step with radt = 1 (docs/RRTMG_FP32_EXPERIMENT.md). CPU builds and FP64
// builds: double, the same code as with double written out (rr(x) is x).
#include "core/precision.hpp"
namespace cppwrf::phys {
#if defined(CPPWRF_USE_CUDA)
using rr = ::cppwrf::real;
#else
using rr = double;
#endif
}  // namespace cppwrf::phys
