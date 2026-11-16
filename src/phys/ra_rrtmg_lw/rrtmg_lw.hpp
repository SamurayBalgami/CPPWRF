// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "phys/ra_rrtmg_common/rrtmg_real.hpp"

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
  rr cldfmc[140][RRTMG_LW_MAX_LAYERS];
  rr taucmc[140][RRTMG_LW_MAX_LAYERS];
  std::vector<rr> d;
  std::vector<int> i;
};

struct RrtmgLwConfig {
  int icloud{1};
  int cldovrlp{2};
  int idcor{0};
  int mp_physics{0};
  rr p_top{rr(5000.0)};
  rr yr{rr(2000.0)};
  rr julian{rr(1.0)};
  rr g{rr(9.81)};
  rr r_d{rr(287.0)};
  rr cp{rr(1004.5)};
  int o3input{0};  // 2: o33d (CAM climatology from the radiation driver) on model levels
  // trace gases (volume mixing ratios): ghg_input = 1 -> read_CAMgases, else rrtmg_default_gases
  rr co2{0}, n2o{0}, ch4{0}, cfc11{0}, cfc12{0};
};

// call configuration of the stages (nlayers, heatfac, gases) and the RRTMG_LWRAD constants
RrtmgLwCallConfig rrtmg_lw_call_config(int num_levels, const RrtmgLwConfig& config, bool hrc, bool hri, bool hrs);
const RrtmgLwConsts& rrtmg_lw_consts();

// RRTMG_LWRAD trace gases without ghg_input (WRF's fixed values, CO2 from the year)
void rrtmg_default_gases(RrtmgLwConfig& c);

// WRF rrtmg_lw (inatm, cldprmc, setcoef, taumol, rtrnmc) for one column whose McICA subcolumn
// cloud arrays are already in ws (cldfmc, taucmc). Layer arrays have nlayers
// entries, level arrays nlayers + 1 (index 0 = surface); pressures in mb, emis per band.
void rrtmg_lw_core(int nlayers, rr cp, const rr* play, const rr* plev, const rr* tlay,
                   const rr* tlev, rr tsk, const rr* h2ovmr, const rr* o3vmr, const rr* co2vmr,
                   const rr* ch4vmr, const rr* n2ovmr, const rr* o2vmr, const rr* cfc11vmr,
                   const rr* cfc12vmr, const rr* cfc22vmr, const rr* ccl4vmr, const rr* semiss,
                   const rr* rei, const rr* rel, const rr* res, rr* totuflux, rr* totdflux,
                   rr* htr, rr* totuclfl, rr* totdclfl, rr* htrc, RrtmgLwWorkspace& ws,
                   const RrtmgLwTables& tables, int inflag = 2, int iceflag = 3, int liqflag = 1);

// Single column execution
void rrtmg_lw_column(
    int num_levels, // model nz (e.g. 50)
    const RrtmgLwConfig& config,
    const RrtmgLwTables& tables,
    // 2D surface variables for this column
    rr tsk,
    rr emiss,
    rr xland,
    rr xice,
    rr snow,
    rr xlat,
    // 3D vertical profiles for this column
    const rr* p8w,   // [nz + 1] Pa
    const rr* t8w,   // [nz + 1] K
    const rr* p3d,   // [nz] Pa
    const rr* pi3d,  // [nz]
    const rr* dz8w,  // [nz] m
    const rr* t3d,   // [nz] K
    const rr* rho3d, // [nz] kg/m3
    const rr* qv3d,  // [nz] kg/kg
    const rr* qc3d,  // [nz] kg/kg
    const rr* qr3d,  // [nz] kg/kg
    const rr* qi3d,  // [nz] kg/kg
    const rr* qs3d,  // [nz] kg/kg
    const rr* qg3d,  // [nz] kg/kg
    const rr* cldfra3d, // [nz]
    const rr* o33d,     // [nz] ozone vmr (o3input = 2), else unused
    // 2D outputs
    rr& glw,
    rr& olr,
    rr& lwcf,
    rr& lwupt,
    rr& lwuptc,
    rr& lwdnt,
    rr& lwdntc,
    rr& lwupb,
    rr& lwupbc,
    rr& lwdnb,
    rr& lwdnbc,
    // 3D outputs
    rr* rthratenlw,  // [nz] K/s
    rr* rthratenlwc, // [nz] K/s
    // Workspace
    RrtmgLwWorkspace& ws,
    // effective radii from the microphysics [nz] (m): WRF has_reqc / has_reqi / has_reqs != 0
    // when given (nullptr: CAM relcalc / reicalc sizes, no separate snow)
    const rr* re_cloud = nullptr,
    const rr* re_ice = nullptr,
    const rr* re_snow = nullptr,
    // [nz] o3input != 2: the climatology ozone vmr on the model levels, WRF's O33D output (nullptr: not wanted)
    rr* o3_out = nullptr
);

} // namespace cppwrf::phys
