#pragma once

#include "core/backend.hpp"
#include <cmath>
#include <cstdint>

namespace cppwrf::phys {

// Marsaglia's 32-bit KISS ("Keep It Simple Stupid") PRNG
// Bit-exact equivalent to WRF module_ra_rrtmg_lw / module_ra_rrtmg_sw kissvec.
struct KissRng {
  uint32_t seed1{0};
  uint32_t seed2{0};
  uint32_t seed3{0};
  uint32_t seed4{0};

  CPPWRF_HD void init(uint32_t s1, uint32_t s2, uint32_t s3, uint32_t s4) {
    seed1 = s1;
    seed2 = s2;
    seed3 = s3;
    seed4 = s4;
  }

  // WRF McICA seed derivation from bottom 4 layer mid-pressures (in Pa)
  // Fortran: seed1 = (pmid(1) - int(pmid(1))) * 1000000000_im
  CPPWRF_HD void init_from_pressure(double pmid0, double pmid1, double pmid2, double pmid3) {
    seed1 = static_cast<uint32_t>((pmid0 - static_cast<int64_t>(pmid0)) * 1.0e9);
    seed2 = static_cast<uint32_t>((pmid1 - static_cast<int64_t>(pmid1)) * 1.0e9);
    seed3 = static_cast<uint32_t>((pmid2 - static_cast<int64_t>(pmid2)) * 1.0e9);
    seed4 = static_cast<uint32_t>((pmid3 - static_cast<int64_t>(pmid3)) * 1.0e9);
  }

  // Advance PRNG by one step and return uniform float in [0, 1)
  template <typename T = double>
  CPPWRF_HD T step() {
    seed1 = 69069u * seed1 + 1327217885u;
    seed2 ^= (seed2 << 13);
    seed2 ^= (seed2 >> 17);
    seed2 ^= (seed2 << 5);
    seed3 = 18000u * (seed3 & 65535u) + (seed3 >> 16);
    seed4 = 30903u * (seed4 & 65535u) + (seed4 >> 16);

    uint32_t kiss_u = seed1 + seed2 + (seed3 << 16) + seed4;
    int32_t kiss = static_cast<int32_t>(kiss_u);
    return static_cast<T>(kiss) * static_cast<T>(2.328306e-10) + static_cast<T>(0.5);
  }

  // Advance PRNG by count steps without returning values (e.g. changeSeed = 150)
  CPPWRF_HD void permute(int count) {
    for (int i = 0; i < count; ++i) {
      step<double>();
    }
  }
};

}  // namespace cppwrf::phys
