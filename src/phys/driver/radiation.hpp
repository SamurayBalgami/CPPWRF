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
  // oznini and the ozone interpolation run on domain 1 only (module_physics_init.F, module_radiation_driver.F,
  // EM_CORE): a nest's radiation uses o3rad, which med_nest_force interpolates from the parent (Registry f=p2c)
  bool oznini = true;
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
  // icloud_bl = 1 with the MYNN PBL: subgrid clouds (module_radiation_driver.F, ICLOUD == 1): after step 1
  // CLDFRA = CLDFRA_BL; QC_BL / QI_BL are added where qc < 1e-6 / qi < 1e-8 and CLDFRA_BL > 0.001
  const Field3 *cldfra_bl = nullptr, *qc_bl = nullptr, *qi_bl = nullptr;
};

struct RadiationOutputs {
  Field3* rthraten;  // uncoupled theta tendency (K/s), LW + SW
  Field2 *glw, *gsw, *swdown, *coszen;
  Field3* cldfra;    // cloud fraction (icloud = 1: cal_cldfra1), may be nullptr
  // diagnostics (may be nullptr): RTHRATLW, RTHRATSW; LWUPT, LWUPTC, LWDNT, LWDNTC, LWUPB, LWUPBC, LWDNB,
  // LWDNBC, OLR; SWUPT, SWUPTC, SWDNT, SWDNTC, SWUPB, SWUPBC, SWDNB, SWDNBC (0 at night)
  Field3 *rthratlw = nullptr, *rthratsw = nullptr;
  Field3 *rthratlwc = nullptr, *rthratswc = nullptr, *o3rad = nullptr;  // clear-sky heating, ozone (o3input = 2)
  Field2 *lwcf = nullptr, *swcf = nullptr, *hrang = nullptr;
  Field2* lw[9] = {};
  Field2* sw[8] = {};
  // RRTMG_SWRAD's surface components (restart fields, may be nullptr; 0 at night): SWVISDIR, SWVISDIF, SWNIRDIR,
  // SWNIRDIF, SWDDIR, SWDDNI, SWDDIF, SWDDIRC, SWDDNIC (DNI = the direct flux / coszen); radiation_driver's
  // DIFFUSE_FRAC = min(SWDDIF / SWDOWN, 1) where SWDOWN > 0.001, else 0
  Field2* swsfc[10] = {};
};

// grid%julian at the start of a step elapsed_s seconds after the simulation start (julyr, julday, gmt):
// domain_clock_get's REAL(dayOfYear_r8) - 1 with dayOfYear_r8 = (seconds since 1 January of the current year)
// / 86400 + 1 (ESMF_TimeGet); 0 again after a year change
inline double wrf_julian(int julyr, int julday, double gmt, double elapsed_s) {
  double s = double(julday - 1) * 86400.0 + double(std::llround(gmt * 3600.0)) + elapsed_s;
  for (int yr = julyr;; ++yr) {
    const double ylen = (((yr % 4 == 0 && yr % 100 != 0) || yr % 400 == 0) ? 366.0 : 365.0) * 86400.0;
    if (s < ylen) break;
    s -= ylen;
  }
  return (s / 86400.0 + 1.0) - 1.0;
}

// the year of that time (first_rk_step_part1: WRFU_TimeGet(currentTime, YY=yr)), the YR of surface_driver
inline int wrf_year(int julyr, int julday, double gmt, double elapsed_s) {
  double s = double(julday - 1) * 86400.0 + double(std::llround(gmt * 3600.0)) + elapsed_s;
  int yr = julyr;
  for (;; ++yr) {
    const double ylen = (((yr % 4 == 0 && yr % 100 != 0) || yr % 400 == 0) ? 366.0 : 365.0) * 86400.0;
    if (s < ylen) break;
    s -= ylen;
  }
  return yr;
}

// calc_coszen (module_radiation_driver.F) for one column; xtime in minutes, angles in degrees
// li, ln: lane li of WRF's n = ln-trip i loop (its..ite = 0..nx-1), which gfortran vectorises: sin(xxlat),
// cos(xxlat) and cos(hrang) through libmvec's 2-lane functions (cppwrf::vm); ln = 1: scalar libm
CPPWRF_HD inline double calc_coszen(double julian, double xtime, double gmt, double declin, double xlon,
                                    double xlat, int li = 0, int ln = 1, double* hrang_out = nullptr) {
  constexpr double piconst = 3.1415926535897932384626433;
  constexpr double degrad = piconst / 180.;
  const double da = 6.2831853071795862 * (julian - 0.5) / 365.;
  const double eot = (0.000075 + 0.001868 * std::cos(da) - 0.032077 * std::sin(da) - 0.014615 * std::cos(2 * da) -
                      0.040849 * std::sin(2 * da)) * (229.18);
  const double xt24 = std::fmod(xtime, 1440.) + eot;
  const double tloctm = gmt + xt24 / 60. + xlon / 15.;
  const double hrang = 15. * (tloctm - 12.) * degrad;
  if (hrang_out) *hrang_out = hrang;
  const double xxlat = xlat * degrad;
  double c = vm::sin(xxlat, li, ln) * std::sin(declin) + vm::cos(xxlat, li, ln) * std::cos(declin) * vm::cos(hrang, li, ln);
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
  double declin() const { return declin_; }  // radconst of the last call (DECLINX, SOLCONX)
  double solcon() const { return solcon_; }

 private:
  struct Impl;
  std::unique_ptr<Impl> p_;
  Grid g_;
  RadiationOptions o_;
  int stepra_ = 1;
  double declin_ = 0, solcon_ = 0;
};

}  // namespace cppwrf::phys::driver
