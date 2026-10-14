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

// Per-thread column workspace to avoid dynamic allocations or stack overflow
struct RrtmgLwWorkspace {
  double cldfmc[140][RRTMG_LW_MAX_LAYERS];
  double ciwpmc[140][RRTMG_LW_MAX_LAYERS];
  double clwpmc[140][RRTMG_LW_MAX_LAYERS];
  double cswpmc[140][RRTMG_LW_MAX_LAYERS];
  double taucmc[140][RRTMG_LW_MAX_LAYERS];
  double fracs[RRTMG_LW_MAX_LAYERS][140];
  double taug[RRTMG_LW_MAX_LAYERS][140];
  double taut[RRTMG_LW_MAX_LAYERS][140];
  double planklay[RRTMG_LW_MAX_LAYERS][16];
  double planklev[RRTMG_LW_MAX_LAYERS + 1][16];
  double wkl[7][RRTMG_LW_MAX_LAYERS];
  double wx[4][RRTMG_LW_MAX_LAYERS];

  // Pointer arrays for calling module kernels
  double* cldfmc_ptrs[140];
  double* ciwpmc_ptrs[140];
  double* clwpmc_ptrs[140];
  double* cswpmc_ptrs[140];
  double* taucmc_ptrs[140];
  double* fracs_ptrs[RRTMG_LW_MAX_LAYERS];
  double* taug_ptrs[RRTMG_LW_MAX_LAYERS];
  double* taut_ptrs[RRTMG_LW_MAX_LAYERS];
  double* planklay_ptrs[RRTMG_LW_MAX_LAYERS];
  double* planklev_ptrs[RRTMG_LW_MAX_LAYERS + 1];
  double* wkl_ptrs[7];
  double* wx_ptrs[4];

  // stride-1 storage of the stage views (rrtmg_lw_stages.hpp) of the column driver
  std::vector<double> d;
  std::vector<int> i;

  void init_pointers() {
    for (int i = 0; i < 140; ++i) {
      cldfmc_ptrs[i] = cldfmc[i];
      ciwpmc_ptrs[i] = ciwpmc[i];
      clwpmc_ptrs[i] = clwpmc[i];
      cswpmc_ptrs[i] = cswpmc[i];
      taucmc_ptrs[i] = taucmc[i];
    }
    for (int l = 0; l < RRTMG_LW_MAX_LAYERS; ++l) {
      fracs_ptrs[l] = fracs[l];
      taug_ptrs[l] = taug[l];
      taut_ptrs[l] = taut[l];
      planklay_ptrs[l] = planklay[l];
      planklev_ptrs[l] = planklev[l];
    }
    planklev_ptrs[RRTMG_LW_MAX_LAYERS] = planklev[RRTMG_LW_MAX_LAYERS];
    for (int m = 0; m < 7; ++m) {
      wkl_ptrs[m] = wkl[m];
    }
    for (int x = 0; x < 4; ++x) {
      wx_ptrs[x] = wx[x];
    }
  }
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
// cloud arrays are already in ws (cldfmc, ciwpmc, clwpmc, cswpmc). Layer arrays have nlayers
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
    const double* re_snow = nullptr
);

// Multi-column driver (OpenMP CPU or sequential)
void rrtmg_lwrad(
    int num_cols,
    int nz,
    const RrtmgLwConfig& config,
    const RrtmgLwTables& tables,
    const double* tsk,
    const double* emiss,
    const double* xland,
    const double* xice,
    const double* snow,
    const double* xlat,
    const double* p8w,
    const double* t8w,
    const double* p3d,
    const double* pi3d,
    const double* dz8w,
    const double* t3d,
    const double* rho3d,
    const double* qv3d,
    const double* qc3d,
    const double* qr3d,
    const double* qi3d,
    const double* qs3d,
    const double* qg3d,
    const double* cldfra3d,
    double* glw,
    double* olr,
    double* lwcf,
    double* lwupt,
    double* lwuptc,
    double* lwdnt,
    double* lwdntc,
    double* lwupb,
    double* lwupbc,
    double* lwdnb,
    double* lwdnbc,
    double* rthratenlw,
    double* rthratenlwc
);

} // namespace cppwrf::phys
