#pragma once
// Shared by the PhysicsSuite translation units (suite.cpp: construction, part1 / part2 and the carried
// state; surface_driver.cpp, pbl_driver.cpp, cumulus_driver.cpp: WRF's module_surface_driver,
// module_pbl_driver and module_cumulus_driver).

#include "core/field.hpp"
#include "phys/lsm_noah/noah.hpp"

namespace cppwrf::phys::driver::suite_detail {

using V2 = View2<real>;
using C2 = View2<const real>;
using V3 = View3<real>;
using C3 = View3<const real>;
using I2 = View2<int>;
using CI2 = View2<const int>;

constexpr int kMaxSoil = noah::kMaxSoilLevels;

inline C2 slice0(const Field3& f) { return C2{f.data(), int(f.layout().sj)}; }  // level k = 0 (WRF kts)

}  // namespace cppwrf::phys::driver::suite_detail
