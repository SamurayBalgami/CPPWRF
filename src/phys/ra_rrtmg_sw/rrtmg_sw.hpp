// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "phys/ra_rrtmg_common/rrtmg_real.hpp"

#include "core/backend.hpp"
#include "phys/ra_rrtmg_common/kissvec.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_tables.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_cldprmc.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_setcoef.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_reftra.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_vrtqdr.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_taumol.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_stages.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace cppwrf::phys {

// Per-thread column workspace: the stage views of rrtmg_sw_stages.hpp (stride 1); cldfmc holds the McICA
// subcolumn flags rrtmg_sw_core takes ([layer][g-point])
struct RrtmgSwWorkspace {
  rr cldfmc[RRTMG_SW_COL_MAX * 112];
  std::vector<rr> d;
  std::vector<int> i;
};

struct RrtmgSwConfig {
  int icloud{1};
  int cldovrlp{2};
  int idcor{0};
  rr p_top{rr(5000.0)};
  rr yr{rr(2000.0)};
  rr julian{rr(1.0)};
  rr solcon{rr(1367.0)};
  rr g{rr(9.81)};
  rr r_d{rr(287.0)};
  rr cp{rr(1004.5)};
  int o3input{0};                     // 2: o33d on model levels
  rr co2{0}, n2o{0}, ch4{0};      // trace gases (vmr), as for RrtmgLwConfig
  rr heatfac{rr(8.4349451468392243)};  // rrsw_con heatfac (rrtmg_sw_ini), checked against WRF
};

RrtmgSwCallConfig rrtmg_sw_call_config(int num_levels, const RrtmgSwConfig& config);

void rrtmg_sw_core(int nlayers, rr heatfac, const rr* play, const rr* plev, const rr* tlay,
                   const rr* tlev, const rr* h2ovmr, const rr* o3vmr, const rr* co2vmr,
                   const rr* ch4vmr, const rr* n2ovmr, const rr* o2vmr, rr asdir, rr asdif,
                   rr aldir, rr aldif, rr coszen, rr adjes, rr scon, const rr* reicmcl,
                   const rr* relqmcl, const rr* resnmcl, rr* swuflx, rr* swdflx, rr* swhr,
                   rr* swuflxc, rr* swdflxc, rr* swhrc, rr* swdkdir, rr* swdkdif,
                   rr* sibvisdir, rr* sibvisdif, rr* sibnirdir, rr* sibnirdif, RrtmgSwWorkspace& ws,
                   const RrtmgSwTables& tables, int inflag = 2, int iceflag = 3, int liqflag = 1);

// Single column execution
void rrtmg_sw_column(
    int num_levels, // model nz (e.g. 50)
    const RrtmgSwConfig& config,
    const RrtmgSwTables& tables,
    // 2D surface / solar variables for this column
    rr coszen,
    rr albedo,
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
    // 2D outputs (gsw, swcf, the TOA / surface fluxes, the surface direct / diffuse components)
    SwColumnOut& out,
    // 3D outputs
    rr* rthratensw,  // [nz] K/s
    rr* rthratenswc, // [nz] K/s
    // Workspace
    RrtmgSwWorkspace& ws,
    // effective radii from the microphysics [nz] (m): WRF has_reqc / has_reqi / has_reqs != 0
    const rr* re_cloud = nullptr,
    const rr* re_ice = nullptr,
    const rr* re_snow = nullptr
);

} // namespace cppwrf::phys
