#pragma once

#include "core/backend.hpp"
#include "phys/ra_rrtmg_common/kissvec.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_tables.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_cldprmc.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_setcoef.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_reftra.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_vrtqdr.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_taumol.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_spcvmc.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace cppwrf::phys {

constexpr int RRTMG_SW_MAX_LAYERS = 128;

// Per-thread column workspace to avoid dynamic allocations
struct RrtmgSwWorkspace {
  double cldfmc[RRTMG_SW_MAX_LAYERS * 112];
  double ciwpmc[RRTMG_SW_MAX_LAYERS * 112];
  double clwpmc[RRTMG_SW_MAX_LAYERS * 112];
  double cswpmc[RRTMG_SW_MAX_LAYERS * 112];
  double taucmc[RRTMG_SW_MAX_LAYERS * 112];
  double taormc[RRTMG_SW_MAX_LAYERS * 112];
  double ssacmc[RRTMG_SW_MAX_LAYERS * 112];
  double asmcmc[RRTMG_SW_MAX_LAYERS * 112];
  double fsfcmc[RRTMG_SW_MAX_LAYERS * 112];
  double wkl[7][RRTMG_SW_MAX_LAYERS];
  double* wkl_ptrs[7];

  CPPWRF_HD void init_pointers() {
    for (int m = 0; m < 7; ++m) {
      wkl_ptrs[m] = wkl[m];
    }
  }
};

struct RrtmgSwConfig {
  int icloud{1};
  int cldovrlp{2};
  int idcor{0};
  double p_top{5000.0};
  double yr{2000.0};
  double julian{1.0};
  double solcon{1367.0};
  double g{9.81};
  double r_d{287.0};
  double cp{1004.5};
  int o3input{0};                     // 2: o33d on model levels
  double co2{0}, n2o{0}, ch4{0};      // trace gases (vmr), as for RrtmgLwConfig
  double heatfac{8.4349451468392243};  // rrsw_con heatfac (rrtmg_sw_ini), checked against WRF
};

void rrtmg_sw_core(int nlayers, double heatfac, const double* play, const double* plev, const double* tlay,
                   const double* tlev, const double* h2ovmr, const double* o3vmr, const double* co2vmr,
                   const double* ch4vmr, const double* n2ovmr, const double* o2vmr, double asdir, double asdif,
                   double aldir, double aldif, double coszen, double adjes, double scon, const double* reicmcl,
                   const double* relqmcl, const double* resnmcl, double* swuflx, double* swdflx, double* swhr,
                   double* swuflxc, double* swdflxc, double* swhrc, double* swdkdir, double* swdkdif,
                   double* sibvisdir, double* sibvisdif, double* sibnirdir, double* sibnirdif, RrtmgSwWorkspace& ws,
                   const RrtmgSwTables& tables);

// Single column execution
CPPWRF_HD void rrtmg_sw_column(
    int num_levels, // model nz (e.g. 50)
    const RrtmgSwConfig& config,
    const RrtmgSwTables& tables,
    // 2D surface / solar variables for this column
    double coszen,
    double albedo,
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
    double& gsw,
    double& swcf,
    double& swupt,
    double& swuptc,
    double& swdnt,
    double& swdntc,
    double& swupb,
    double& swupbc,
    double& swdnb,
    double& swdnbc,
    double& swvisdir,
    double& swvisdif,
    double& swnirdir,
    double& swnirdif,
    // 3D outputs
    double* rthratensw,  // [nz] K/s
    double* rthratenswc, // [nz] K/s
    // Workspace
    RrtmgSwWorkspace& ws
);

// Multi-column driver (OpenMP CPU or sequential)
void rrtmg_swrad(
    int num_cols,
    int nz,
    const RrtmgSwConfig& config,
    const RrtmgSwTables& tables,
    const double* coszen,
    const double* albedo,
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
    double* gsw,
    double* swcf,
    double* swupt,
    double* swuptc,
    double* swdnt,
    double* swdntc,
    double* swupb,
    double* swupbc,
    double* swdnb,
    double* swdnbc,
    double* swvisdir,
    double* swvisdif,
    double* swnirdir,
    double* swnirdif,
    double* rthratensw,
    double* rthratenswc
);

} // namespace cppwrf::phys
