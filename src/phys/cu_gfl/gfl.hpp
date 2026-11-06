// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Grell-Freitas-Li (cu_physics = 3, GFLSCHEME): gfl_common.hpp, gfl_deep.hpp, gfl_shallow.hpp and gfl_driver.hpp,
// compiled once per level tier. The scheme's column arrays have KM = levels + 12 elements (1-based levels and the
// look-ahead above the top); CUDA reserves a kernel's stack frame for every thread the GPU holds at once (FP64:
// 128 KB per thread with KM = 112), so a domain of up to 48 levels takes the 60-element arrays. Only the arrays'
// zero initialisations reach index KM - 1; every other access stays below levels + 12 (the farthest, levels + 8,
// in get_inversion_layers reads the zeros above the top), so the tier changes no value.
#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

#include "core/backend.hpp"

#ifndef CPPWRF_HOST_DEVICE
#define CPPWRF_HOST_DEVICE CPPWRF_HD
#endif

namespace cppwrf::phys::cu_gfl {
inline constexpr int kMaxLevels = 100;
}  // namespace cppwrf::phys::cu_gfl

namespace cppwrf::phys::cu_gfl::km48 {
inline constexpr int KM = 48 + 12;
#include "phys/cu_gfl/gfl_common.hpp"
#include "phys/cu_gfl/gfl_deep.hpp"
#include "phys/cu_gfl/gfl_shallow.hpp"
#include "phys/cu_gfl/gfl_driver.hpp"
}  // namespace cppwrf::phys::cu_gfl::km48

namespace cppwrf::phys::cu_gfl::km100 {
inline constexpr int KM = kMaxLevels + 12;
#include "phys/cu_gfl/gfl_common.hpp"
#include "phys/cu_gfl/gfl_deep.hpp"
#include "phys/cu_gfl/gfl_shallow.hpp"
#include "phys/cu_gfl/gfl_driver.hpp"
}  // namespace cppwrf::phys::cu_gfl::km100

namespace cppwrf::phys::cu_gfl {
// the tier of up to L levels: its column types and cu_grell_freitas_li
template <int L, typename real>
struct Tier;
template <typename real>
struct Tier<48, real> {
  static constexpr int KM = km48::KM;
  using Column = km48::GflDrvColumn<real>;
  using Out = km48::GflDrvOut<real>;
  CPPWRF_HD static void run(const Column& c, int ishallow, int ichoice, Out& o) {
    km48::gfl_column(c, ishallow, ichoice, o);
  }
};
template <typename real>
struct Tier<kMaxLevels, real> {
  static constexpr int KM = km100::KM;
  using Column = km100::GflDrvColumn<real>;
  using Out = km100::GflDrvOut<real>;
  CPPWRF_HD static void run(const Column& c, int ishallow, int ichoice, Out& o) {
    km100::gfl_column(c, ishallow, ichoice, o);
  }
};
}  // namespace cppwrf::phys::cu_gfl
