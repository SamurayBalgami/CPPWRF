#pragma once

#include "core/backend.hpp"
#include <algorithm>
#include <cmath>

namespace cppwrf::phys {

static constexpr double RETAB[95] = {
    5.92779, 6.26422, 6.61973, 6.99539, 7.39234, 7.81177, 8.25496, 8.72323, 9.21800, 9.74075,
    10.2930, 10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319, 15.0351, 15.8799, 16.7674,
    17.6986, 18.6744, 19.6955, 20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125, 27.7895,
    28.6450, 29.4167, 30.1088, 30.7306, 31.2943, 31.8151, 32.3077, 32.7870, 33.2657, 33.7540,
    34.2601, 34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078, 38.4720, 39.1508, 39.8442,
    40.5552, 41.2912, 42.0635, 42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221, 50.4710,
    52.4980, 54.8315, 57.4898, 60.4785, 63.7898, 65.5604, 71.2885, 75.4113, 79.7368, 84.2351,
    88.8833, 93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424, 124.954, 130.630, 136.457,
    142.446, 148.608, 154.956, 161.503, 168.262, 175.248, 182.473, 189.952, 197.699, 205.728,
    214.055, 222.694, 231.661, 240.971, 250.639};

// Climatological ozone profile calculation (INIRAD from WRF module_ra_rrtmg_lw)
CPPWRF_HD inline void rrtmg_inirad(
    double* o3mmr,
    const double* plev, // [nlayers + 1] in mb
    int kts,
    int kte
) {
  static constexpr double o3sum[31] = {
      5.297e-8, 5.852e-8, 6.579e-8, 7.505e-8, 8.577e-8, 9.895e-8, 1.175e-7, 1.399e-7,
      1.677e-7, 2.003e-7, 2.571e-7, 3.325e-7, 4.438e-7, 6.255e-7, 8.168e-7, 1.036e-6,
      1.366e-6, 1.855e-6, 2.514e-6, 3.240e-6, 4.033e-6, 4.854e-6, 5.517e-6, 6.089e-6,
      6.689e-6, 1.106e-5, 1.462e-5, 1.321e-5, 9.856e-6, 5.960e-6, 5.960e-6};
  static constexpr double ppsum[31] = {
      955.890, 850.532, 754.599, 667.742, 589.841, 519.421, 455.480, 398.085,
      347.171, 301.735, 261.310, 225.360, 193.419, 165.490, 141.032, 120.125,
      102.689, 87.829,  75.123,  64.306,  55.086,  47.209,  40.535,  34.795,
      29.865,  19.122,  9.277,   4.660,   2.421,   1.294,   0.647};
  static constexpr double o3win[31] = {
      4.629e-8, 4.686e-8, 5.017e-8, 5.613e-8, 6.871e-8, 8.751e-8, 1.138e-7, 1.516e-7,
      2.161e-7, 3.264e-7, 4.968e-7, 7.338e-7, 1.017e-6, 1.308e-6, 1.625e-6, 2.011e-6,
      2.516e-6, 3.130e-6, 3.840e-6, 4.703e-6, 5.486e-6, 6.289e-6, 6.993e-6, 7.494e-6,
      8.197e-6, 9.632e-6, 1.113e-5, 1.146e-5, 9.389e-6, 6.135e-6, 6.135e-6};
  static constexpr double ppwin[31] = {
      955.747, 841.783, 740.199, 649.538, 568.404, 495.815, 431.069, 373.464,
      322.354, 277.190, 237.635, 203.433, 174.070, 148.949, 127.408, 108.915,
      93.114,  79.551,  67.940,  58.072,  49.593,  42.318,  36.138,  30.907,
      26.362,  16.423,  7.583,   3.620,   1.807,   0.938,   0.469};

  double o3ann[31];
  o3ann[0] = 0.5 * (o3sum[0] + o3win[0]);
  for (int k = 1; k < 31; ++k) {
    o3ann[k] = o3win[k - 1] + (o3win[k] - o3win[k - 1]) / (ppwin[k] - ppwin[k - 1]) * (ppsum[k] - ppwin[k - 1]);
  }
  for (int k = 1; k < 31; ++k) {
    o3ann[k] = 0.5 * (o3ann[k] + o3sum[k]);
  }

  double ppwrkh[32];
  ppwrkh[0] = 1100.0;
  for (int k = 1; k < 31; ++k) {
    ppwrkh[k] = 0.5 * (ppsum[k] + ppsum[k - 1]);
  }
  ppwrkh[31] = 0.0;

  for (int k = kts; k <= kte + 1; ++k) {
    o3mmr[k] = 0.0;
    double pr_k = plev[k];
    double pr_k1 = plev[k + 1];
    for (int jj = 0; jj < 31; ++jj) {
      double pb1 = (pr_k > ppwrkh[jj]) ? (pr_k - ppwrkh[jj]) : 0.0;
      double pb2 = (pr_k > ppwrkh[jj + 1]) ? (pr_k - ppwrkh[jj + 1]) : 0.0;
      double pt1 = (pr_k1 > ppwrkh[jj]) ? (pr_k1 - ppwrkh[jj]) : 0.0;
      double pt2 = (pr_k1 > ppwrkh[jj + 1]) ? (pr_k1 - ppwrkh[jj + 1]) : 0.0;
      o3mmr[k] += (pb2 - pb1 - pt2 + pt1) * o3ann[jj];
    }
    o3mmr[k] /= (pr_k - pr_k1);
  }
}

// CAM relcalc: liquid droplet effective radius (microns)
CPPWRF_HD inline double rrtmg_relcalc_single(
    double t_lay,
    double xland,
    double snow,
    double xice
) {
  constexpr double tmelt = 273.16;
  constexpr double rliqocean = 14.0;
  constexpr double rliqice = 14.0;
  constexpr double rliqland = 8.0;

  double landfrac = 2.0 - xland;
  double landm = landfrac;
  double snowh = 0.001 * snow;
  double icefrac = xice;

  double r = rliqland + (rliqocean - rliqland) * std::min(1.0, std::max(0.0, (tmelt - t_lay) * 0.05));
  r += (rliqocean - r) * std::min(1.0, std::max(0.0, snowh * 10.0));
  r += (rliqocean - r) * std::min(1.0, std::max(0.0, 1.0 - landm));
  r += (rliqice - r) * std::min(1.0, std::max(0.0, icefrac));
  return r;
}

// CAM reicalc: ice effective drop size (microns)
CPPWRF_HD inline double rrtmg_reicalc_single(double t_lay) {
  int idx = static_cast<int>(t_lay - 179.0);
  idx = std::min(std::max(idx, 1), 94);
  double corr = t_lay - static_cast<int>(t_lay);
  return RETAB[idx - 1] * (1.0 - corr) + RETAB[idx] * corr;
}

} // namespace cppwrf::phys
