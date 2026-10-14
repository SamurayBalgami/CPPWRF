#pragma once
// WRF CAM ozone climatology for o3input = 2 (module_ra_cam_support.F oznini, module_radiation_driver.F
// ozn_time_int, ozn_p_int): monthly zonal-mean ozone (ozone.formatted, 59 pressure levels x 64
// latitudes), interpolated to each column's latitude once (oznini), then per radiation call in
// time (julian) and to the column's pressures.
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "phys/ra_rrtmg_common/strided.hpp"

namespace cppwrf::phys {

struct CamOzone {
  static constexpr int levsiz = 59, latsiz = 64, nmon = 12;
  std::vector<double> pin;      // [levsiz] Pa
  std::vector<double> lat;      // [latsiz]
  std::vector<double> ozmixin;  // [m][j][k] (month, latitude, level) as read
  void load(const std::string& dir);
  // oznini: ozmixm for one latitude, [m][k]
  void column(double xlat, double* ozmixm) const;
};

// ozn_time_int: weights of the two bracketing months (0-based) for julian (WRF r8 real)
void ozone_time_weights(double julian, int& nm, int& np, double& fact1, double& fact2);

// ozn_time_int + ozn_p_int for one column: p[nz] (Pa, bottom to top), o3[nz] (vmr)
CPPWRF_HD inline void ozone_column(const double* ozm_month_m, const double* ozm_month_p, double fact1, double fact2,
                                   const double* pin, int nz, SD p, SD o3) {
  constexpr int levsiz = CamOzone::levsiz;
  double ozmixt[levsiz];
  for (int k = 0; k < levsiz; ++k) ozmixt[k] = ozm_month_m[k] * fact1 + ozm_month_p[k] * fact2;
  int kupper = 1;  // Fortran 1-based level index into pin
  for (int k = 1; k <= nz; ++k) {  // k = 1 is the model top (pmid reversed)
    const double pmid = p[nz - k];
    const int kout = nz - k + 1;  // o3vmr(kout) = model level kout (1-based)
    bool found = false;
    for (int kk = kupper; kk <= levsiz - 1; ++kk) {
      if (pin[kk - 1] < pmid && pmid <= pin[kk]) { kupper = kk; found = true; }
    }
    double v;
    if (found) {
      const double dpu = pmid - pin[kupper - 1], dpl = pin[kupper] - pmid;
      v = (ozmixt[kupper - 1] * dpl + ozmixt[kupper] * dpu) / (dpl + dpu);
    } else if (pmid < pin[0]) {
      v = ozmixt[0] * pmid / pin[0];
    } else if (pmid > pin[levsiz - 1]) {
      v = ozmixt[levsiz - 1];
    } else {
      const double dpu = pmid - pin[kupper - 1], dpl = pin[kupper] - pmid;
      v = (ozmixt[kupper - 1] * dpl + ozmixt[kupper] * dpu) / (dpl + dpu);
    }
    o3[kout - 1] = v;
  }
}

}  // namespace cppwrf::phys
