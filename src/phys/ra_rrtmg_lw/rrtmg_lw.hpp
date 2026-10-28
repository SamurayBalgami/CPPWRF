#pragma once

#include "core/backend.hpp"
#include "phys/ra_rrtmg_common/kissvec.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_cldprmc.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_rtrnmc.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_setcoef.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_tables.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_taumol.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_stages.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace cppwrf::phys {

constexpr int RRTMG_LW_MAX_LAYERS = 128;

// Per-thread column workspace: the stride-1 storage of the stage views (rrtmg_lw_stages.hpp) of the
// column driver, and the McICA cloud flags / optical depths that rrtmg_lw_core takes as input
struct RrtmgLwWorkspace {
  double cldfmc[140][RRTMG_LW_MAX_LAYERS];
  double taucmc[140][RRTMG_LW_MAX_LAYERS];
  std::vector<double> d;
  std::vector<int> i;
};

struct RrtmgLwConfig {
  int icloud{1};
  int cldovrlp{2};
  int idcor{0};
  int mp_physics{0};
  double p_top{5000.0};
  double yr{2000.0};
  double julian{1.0};
  double g{9.81};
  double r_d{287.0};
  double cp{1004.5};
  int o3input{0};  // 2: o33d (CAM climatology from the radiation driver) on model levels
  // trace gases (volume mixing ratios): ghg_input = 1 -> read_CAMgases, else rrtmg_default_gases
  double co2{0}, n2o{0}, ch4{0}, cfc11{0}, cfc12{0};
};

// call configuration of the stages (nlayers, heatfac, gases) and the RRTMG_LWRAD constants
RrtmgLwCallConfig rrtmg_lw_call_config(int num_levels, const RrtmgLwConfig& config, bool hrc, bool hri, bool hrs);
const RrtmgLwConsts& rrtmg_lw_consts();

// RRTMG_LWRAD trace gases without ghg_input (WRF's fixed values, CO2 from the year)
void rrtmg_default_gases(RrtmgLwConfig& c);

// WRF rrtmg_lw (inatm, cldprmc, setcoef, taumol, rtrnmc) for one column whose McICA subcolumn
// cloud arrays are already in ws (cldfmc, taucmc). Layer arrays have nlayers
// entries, level arrays nlayers + 1 (index 0 = surface); pressures in mb, emis per band.
void rrtmg_lw_core(int nlayers, double cp, const double* play, const double* plev, const double* tlay,
                   const double* tlev, double tsk, const double* h2ovmr, const double* o3vmr, const double* co2vmr,
                   const double* ch4vmr, const double* n2ovmr, const double* o2vmr, const double* cfc11vmr,
                   const double* cfc12vmr, const double* cfc22vmr, const double* ccl4vmr, const double* semiss,
                   const double* rei, const double* rel, const double* res, double* totuflux, double* totdflux,
                   double* htr, double* totuclfl, double* totdclfl, double* htrc, RrtmgLwWorkspace& ws,
                   const RrtmgLwTables& tables, int inflag = 2, int iceflag = 3, int liqflag = 1);

// Single column execution
void rrtmg_lw_column(
    int num_levels, // model nz (e.g. 50)
    const RrtmgLwConfig& config,
    const RrtmgLwTables& tables,
    // 2D surface variables for this column
    double tsk,
    double emiss,
    double xland,
    double xice,
    double snow,
    double xlat,
    // 3D vertical profiles for this column
    const double* p8w,   // [nz + 1] Pa
    const double* t8w,   // [nz + 1] K
    const double* p3d,   // [nz] Pa
    const double* pi3d,  // [nz]
    const double* dz8w,  // [nz] m
    const double* t3d,   // [nz] K
    const double* rho3d, // [nz] kg/m3
    const double* qv3d,  // [nz] kg/kg
    const double* qc3d,  // [nz] kg/kg
    const double* qr3d,  // [nz] kg/kg
    const double* qi3d,  // [nz] kg/kg
    const double* qs3d,  // [nz] kg/kg
    const double* qg3d,  // [nz] kg/kg
    const double* cldfra3d, // [nz]
    const double* o33d,     // [nz] ozone vmr (o3input = 2), else unused
    // 2D outputs
    double& glw,
    double& olr,
    double& lwcf,
    double& lwupt,
    double& lwuptc,
    double& lwdnt,
    double& lwdntc,
    double& lwupb,
    double& lwupbc,
    double& lwdnb,
    double& lwdnbc,
    // 3D outputs
    double* rthratenlw,  // [nz] K/s
    double* rthratenlwc, // [nz] K/s
    // Workspace
    RrtmgLwWorkspace& ws,
    // effective radii from the microphysics [nz] (m): WRF has_reqc / has_reqi / has_reqs != 0
    // when given (nullptr: CAM relcalc / reicalc sizes, no separate snow)
    const double* re_cloud = nullptr,
    const double* re_ice = nullptr,
    const double* re_snow = nullptr,
    // [nz] o3input != 2: the climatology ozone vmr on the model levels, WRF's O33D output (nullptr: not wanted)
    double* o3_out = nullptr
);

} // namespace cppwrf::phys
