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

constexpr int kMcicaMaxLayers = 128;  // RRTMG_LW/SW_MAX_LAYERS (layer heights of the exponential overlaps)

// generate_stochastic_clouds (module_ra_rrtmg_lw.F / _sw.F, irng = 0): the cloud overlap (icld = cldovrlp)
// 1 random, 2 maximum-random, 3 maximum, 4 exponential, 5 exponential-random. cldf(l): layer cloud fraction
// (0 below cldmin), alpha(l): exponential overlap parameter of layer l (icld 4/5, l >= 1), out(isub, l, cloudy).
// The random numbers are drawn subcolumn by subcolumn, level by level (icld 3: one per subcolumn; icld 4/5:
// CDF then CDF2 at every level), and each level's CDF transform reads the previous level's final CDF.
template <class CF, class AL, class OUT>
CPPWRF_HD void mcica_subcolumns(int icld, KissRng& rng, int nsub, int nlay, CF cldf, AL alpha, OUT out) {
  for (int isub = 0; isub < nsub; ++isub) {
    double cdf_prev = 0.0, cldf_prev = 0.0;
    const double r3 = icld == 3 ? rng.step<double>() : 0.0;
    for (int l = 0; l < nlay; ++l) {
      const double cf = cldf(l);
      double cdf;
      if (icld == 3) {
        cdf = r3;
      } else if (icld == 4 || icld == 5) {
        cdf = rng.step<double>();
        const double cdf2 = rng.step<double>();
        if (l > 0 && cdf2 < alpha(l)) cdf = cdf_prev;
      } else {
        cdf = rng.step<double>();
        if (icld == 2 && l > 0) {
          if (cdf_prev > 1.0 - cldf_prev) cdf = cdf_prev;
          else cdf = cdf * (1.0 - cldf_prev);
        }
      }
      out(isub, l, cdf >= 1.0 - cf);
      cdf_prev = cdf;
      cldf_prev = cf;
    }
  }
}

// decorr_len of mcica_subcol_lw / _sw: idcor 1 from latitude (degrees) and day of year, else 2500 m
CPPWRF_HD inline double mcica_decorr_len(int idcor, int juldat, double lat) {
  if (idcor != 1) return 2500.0;
  constexpr double am1 = 1.4315, am2 = 2.1219, am4 = -25.584, amr = 7.0;
  const double am3 = juldat > 181 ? -4.0 * amr / 365.0 * double(juldat - 272) : 4.0 * amr / 365.0 * double(juldat - 91);
  const double d = lat - am3;
  return (am1 + am2 * exp(-(d * d) / (am4 * am4))) * 1.e3;
}

}  // namespace cppwrf::phys
