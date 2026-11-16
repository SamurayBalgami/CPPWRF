// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "phys/ra_rrtmg_common/rrtmg_real.hpp"

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

  // WRF McICA seed derivation from bottom 4 layer mid-pressures (in Pa)
  // Fortran: seed1 = (pmid(1) - int(pmid(1))) * 1000000000_im
  CPPWRF_HD void init_from_pressure(rr pmid0, rr pmid1, rr pmid2, rr pmid3) {
    seed1 = static_cast<uint32_t>((pmid0 - static_cast<int64_t>(pmid0)) * rr(1.0e9));
    seed2 = static_cast<uint32_t>((pmid1 - static_cast<int64_t>(pmid1)) * rr(1.0e9));
    seed3 = static_cast<uint32_t>((pmid2 - static_cast<int64_t>(pmid2)) * rr(1.0e9));
    seed4 = static_cast<uint32_t>((pmid3 - static_cast<int64_t>(pmid3)) * rr(1.0e9));
  }

  // Advance PRNG by one step and return uniform float in [0, 1)
  template <typename T = rr>
  CPPWRF_HD T step() {
    seed1 = 69069u * seed1 + 1327217885u;
    seed2 ^= (seed2 << 13);
    seed2 ^= (seed2 >> 17);
    seed2 ^= (seed2 << 5);
    seed3 = 18000u * (seed3 & 65535u) + (seed3 >> 16);
    seed4 = 30903u * (seed4 & 65535u) + (seed4 >> 16);

    uint32_t kiss_u = seed1 + seed2 + (seed3 << 16) + seed4;
    int32_t kiss = static_cast<int32_t>(kiss_u);
    return static_cast<T>(kiss) * static_cast<T>(rr(2.328306e-10)) + static_cast<T>(rr(0.5));
  }

  // Advance PRNG by count steps without returning values (e.g. changeSeed = 150)
  CPPWRF_HD void permute(int count) {
    for (int i = 0; i < count; ++i) {
      step<rr>();
    }
  }
};

constexpr int kMcicaMaxLayers = 128;  // RRTMG_LW/SW_MAX_LAYERS (layer heights of the exponential overlaps)

// generate_stochastic_clouds (module_ra_rrtmg_lw.F / _sw.F, irng = 0): the cloud overlap (icld = cldovrlp)
// 1 random, 2 maximum-random, 3 maximum, 4 exponential, 5 exponential-random. cldf(l): layer cloud fraction
// (0 below cldmin), alpha(l): exponential overlap parameter of layer l (icld 4/5, l >= 1), out(isub, l, cloudy).
// The random numbers are drawn subcolumn by subcolumn, level by level (icld 3: one per subcolumn; icld 4/5:
// CDF then CDF2 at every level), and each level's CDF transform reads the previous level's final CDF.
template <class CF, class AL, class OUT>
CPPWRF_HD void mcica_subcolumns(int icld, KissRng& rng, int nsub, int nlay, CF cldf, AL alpha, OUT out) {
  for (int isub = 0; isub < nsub; ++isub) {
    rr cdf_prev = rr(0.0), cldf_prev = rr(0.0);
    const rr r3 = icld == 3 ? rng.step<rr>() : rr(0.0);
    for (int l = 0; l < nlay; ++l) {
      const rr cf = cldf(l);
      rr cdf;
      if (icld == 3) {
        cdf = r3;
      } else if (icld == 4 || icld == 5) {
        cdf = rng.step<rr>();
        const rr cdf2 = rng.step<rr>();
        if (l > 0 && cdf2 < alpha(l)) cdf = cdf_prev;
      } else {
        cdf = rng.step<rr>();
        if (icld == 2 && l > 0) {
          if (cdf_prev > rr(1.0) - cldf_prev) cdf = cdf_prev;
          else cdf = cdf * (rr(1.0) - cldf_prev);
        }
      }
      out(isub, l, cdf >= rr(1.0) - cf);
      cdf_prev = cdf;
      cldf_prev = cf;
    }
  }
}

// decorr_len of mcica_subcol_lw / _sw: idcor 1 from latitude (degrees) and day of year, else 2500 m
CPPWRF_HD inline rr mcica_decorr_len(int idcor, int juldat, rr lat) {
  if (idcor != 1) return rr(2500.0);
  constexpr rr am1 = rr(1.4315), am2 = rr(2.1219), am4 = -rr(25.584), amr = rr(7.0);
  const rr am3 = juldat > 181 ? -rr(4.0) * amr / rr(365.0) * rr(juldat - 272) : rr(4.0) * amr / rr(365.0) * rr(juldat - 91);
  const rr d = lat - am3;
  return (am1 + am2 * libm::exp(-(d * d) / (am4 * am4))) * rr(1.e3);
}

}  // namespace cppwrf::phys
