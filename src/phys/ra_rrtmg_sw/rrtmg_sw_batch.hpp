// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// RRTMG_SWRAD for all columns of a domain on the device (CPU builds: OpenMP), the stages of
// rrtmg_sw_stages.hpp run as kernels over (column), (column, g-point) and (column, level), all columns
// in one chunk where the device memory allows (as RrtmgLwBatch). Arrays have the column index fastest: x[k*ncol + col].
#include <cstddef>
#include <memory>

#include "phys/ra_rrtmg_common/batch_memory.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_stages.hpp"

namespace cppwrf::phys {

struct SwBatchInput {
  int ncol;
  const double *coszen, *albedo, *xland, *xice, *snow, *xlat;   // [ncol]
  const double *dz8w;                                           // [nz*ncol] m
  const double *p8w, *t8w;                                      // [(nz+1)*ncol] Pa, K
  const double *p3d, *pi3d, *t3d, *qv, *cldfra, *o3;            // [nz*ncol]
  const double *qc, *qi, *qs;                                   // nullptr: species absent
  const double *re_cloud, *re_ice, *re_snow;                    // nullptr: has_req* = 0
};

struct SwBatchOutput {
  double* rthraten;  // [nz*ncol] K/s (0 at night)
  double* gsw;       // [ncol] W/m2
  // [8*ncol] (may be nullptr): swupt, swuptc, swdnt, swdntc, swupb, swupbc, swdnb, swdnbc (d*ncol + col)
  double* diag = nullptr;
  double* rthratenc = nullptr;  // [nz*ncol] clear-sky heating rate (may be nullptr)
  double* cf = nullptr;         // [ncol] SW cloud forcing (may be nullptr)
  // [7*ncol] (may be nullptr): swvisdir, swvisdif, swnirdir, swnirdif, swddir, swddif, swddirc (d*ncol + col)
  double* sfc = nullptr;
};

class RrtmgSwBatch {
 public:
  // ncol: the domain's columns (they set the chunk size, rrtmg_chunk_columns)
  RrtmgSwBatch(const RrtmgSwTables& tables, int nlayers, int ncol);
  ~RrtmgSwBatch();
  RrtmgSwBatch(const RrtmgSwBatch&) = delete;
  RrtmgSwBatch& operator=(const RrtmgSwBatch&) = delete;
  void run(const RrtmgSwCallConfig& c, const SwBatchInput& in, const SwBatchOutput& out);

 private:
  RrtmgSwTables* d_tables_ = nullptr;  // shared device copies (rrtmg_shared_dev_copy)
  RrtmgLwConsts* d_consts_ = nullptr;
  SwLayout lay_;
  int chunk_;
  std::unique_ptr<RrtmgScratch> work_;  // the shared work arrays
};

}  // namespace cppwrf::phys
