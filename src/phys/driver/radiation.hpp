#pragma once
// WRF phys/module_radiation_driver.F radiation_driver for ra_lw_physics = 4 / ra_sw_physics = 4
// (RRTMG) with o3input = 2 (CAM ozone), ghg_input = 1 (CAMtr_volume_mixing_ratio), aer_opt = 0,
// swint_opt = 0, icloud = 1, cldovrlp, idcor, slope_rad = topo_shading = 0, single domain.
// CPU builds: one RRTMG column at a time per OpenMP thread (rrtmg_lw_column / rrtmg_sw_column).
// CUDA builds: everything on the device: inputs staged column-fastest, coszen, CAM ozone and
// cal_cldfra1 per column, the RRTMG stages as kernels (RrtmgLwBatch / RrtmgSwBatch).
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "core/field.hpp"
#include "core/grid.hpp"
#include "phys/ra_rrtmg_common/strided.hpp"

namespace cppwrf::phys::driver {

struct RadiationOptions {
  int ra_lw_physics = 4, ra_sw_physics = 4;
  double radt = 0;          // minutes; stepra = max(nint(radt*60/dt), 1)
  int ra_call_offset = 0;
  int icloud = 1, cldovrlp = 2, idcor = 0, o3input = 2, ghg_input = 1, aer_opt = 0;
  int julyr = 0, julday = 1;  // start (wrfinput JULYR/JULDAY), gmt start hour
  double gmt = 0;
  double p_top = 0;
  std::string table_dir = ".";  // RRTMG_*_DATA tables dump, ozone*.formatted, CAMtr_volume_mixing_ratio
};

struct RadiationInputs {  // device fields (mass points) read by the radiation
  const Field3 *t_phy, *p_hyd, *p_hyd_w, *pi_phy, *rho, *dz8w, *t8w, *qv;
  const Field3 *qc, *qi, *qs;  // moist species (nullptr when the scheme has none: WRF F_QC, F_QI, F_QS)
  const Field3 *re_cloud, *re_ice, *re_snow;  // microphysics effective radii (m): has_reqc/i/s, else nullptr
  const Field2 *tsk, *albedo, *emiss, *xland, *xice, *snow;
};

struct RadiationOutputs {
  Field3* rthraten;  // uncoupled theta tendency (K/s), LW + SW
  Field2 *glw, *gsw, *swdown, *coszen;
  Field3* cldfra;    // cloud fraction (icloud = 1: cal_cldfra1), may be nullptr
};

// calc_coszen (module_radiation_driver.F) for one column; xtime in minutes, angles in degrees
CPPWRF_HD inline double calc_coszen(double julian, double xtime, double gmt, double declin, double xlon,
                                    double xlat) {
  constexpr double piconst = 3.1415926535897932384626433;
  constexpr double degrad = piconst / 180.;
  const double da = 6.2831853071795862 * (julian - 1) / 365.;
  const double eot = (0.000075 + 0.001868 * std::cos(da) - 0.032077 * std::sin(da) - 0.014615 * std::cos(2 * da) -
                      0.04089 * std::sin(2 * da)) * (229.18);
  const double xt24 = std::fmod(xtime, 1440.) + eot;
  const double tloctm = gmt + xt24 / 60. + xlon / 15.;
  const double hrang = 15. * (tloctm - 12.) * degrad;
  const double xxlat = xlat * degrad;
  double c = std::sin(xxlat) * std::sin(declin) + std::cos(xxlat) * std::cos(declin) * std::cos(hrang);
  return std::min(std::max(c, -1.0), 1.0);
}

// module_radiation_driver.F cal_cldfra1 for one column (bottom to top): qi / qs with p == nullptr
// when the scheme has no F_QI / F_QS; t (K) and p (Pa) are t_phy and p_hyd. mp_physics !=
// FER_MP_HIRES; F_QC and not F_QI but F_QS needs F_ICE_PHY (not ported).
CPPWRF_HD inline void cal_cldfra1(int nz, SD qv, SD qc, SD qi, SD qs, SD t, SD p, SD cldfra) {
  constexpr double ALPHA0 = 100., GAMMA = 0.49, QCLDMIN = 1.E-12, PEXP = 0.25, RHGRID = 1.0;
  constexpr double SVP1 = 0.61078, SVP2 = 17.2693882, SVPI2 = 21.8745584, SVP3 = 35.86, SVPI3 = 7.66,
                   SVPT0 = 273.15, r_d = 287., r_v = 461.6, ep_2 = r_d / r_v;
  for (int k = 0; k < nz; ++k) {
    const double tk = t[k], pk = p[k];
    const double tc = tk - SVPT0;
    const double esw = 1000.0 * SVP1 * std::exp(SVP2 * tc / (tk - SVP3));
    const double esi = 1000.0 * SVP1 * std::exp(SVPI2 * tc / (tk - SVPI3));
    const double QVSW = ep_2 * esw / (pk - esw);
    const double QVSI = ep_2 * esi / (pk - esi);
    double QCLD = 0, weight = 0;
    if (qi.p && qs.p) {
      QCLD = qi[k] + qc[k] + qs[k];
      weight = QCLD < QCLDMIN ? 0. : (qi[k] + qs[k]) / QCLD;
    } else if (qi.p) {
      QCLD = qi[k] + qc[k];
      weight = QCLD < QCLDMIN ? 0. : qi[k] / QCLD;
    } else {  // F_QC only
      QCLD = qc[k];
      if (QCLD < QCLDMIN) weight = 0.;
      else weight = tk > 273.15 ? 0. : 1.;
    }
    const double QVS_WEIGHT = (1 - weight) * QVSW + weight * QVSI;
    double RHUM = qv[k] / QVS_WEIGHT;
    double cf;
    if (QCLD < QCLDMIN) {
      cf = 0.;
    } else if (RHUM >= RHGRID) {
      cf = 1.;
    } else {
      const double SUBSAT = std::max(1.E-10, RHGRID * QVS_WEIGHT - qv[k]);
      const double DENOM = std::pow(SUBSAT, GAMMA);
      const double ARG = std::max(-6.9, -ALPHA0 * QCLD / DENOM);
      RHUM = std::max(1.E-10, RHUM);
      cf = std::pow(RHUM / RHGRID, PEXP) * (1. - std::exp(ARG));
      if (cf < .01) cf = 0.;
    }
    cldfra[k] = cf;
  }
}

class RadiationDriver {
 public:
  RadiationDriver(const Grid& g, const RadiationOptions& o, const std::vector<double>& xlat,
                  const std::vector<double>& xlong);
  ~RadiationDriver();
  int stepra() const { return stepra_; }
  // radiation_driver for time step itimestep (1-based) at xtime minutes after the start and
  // julian (day of year - 1 + fraction); returns false when the radiation is not called
  bool run(long itimestep, double xtime, double julian, const RadiationInputs& in, const RadiationOutputs& out);

 private:
  struct Impl;
  std::unique_ptr<Impl> p_;
  Grid g_;
  RadiationOptions o_;
  int stepra_ = 1;
};

}  // namespace cppwrf::phys::driver
