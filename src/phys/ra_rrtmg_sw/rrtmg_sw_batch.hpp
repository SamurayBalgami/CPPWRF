// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "phys/ra_rrtmg_common/rrtmg_real.hpp"
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
  int col_begin = 0, col_end = -1;  // the columns computed: [col_begin, col_end) (-1: ncol); a tile's rows
  const rr *coszen, *albedo, *xland, *xice, *snow, *xlat;   // [ncol]
  const rr *dz8w;                                           // [nz*ncol] m
  const rr *p8w, *t8w;                                      // [(nz+1)*ncol] Pa, K
  const rr *p3d, *pi3d, *t3d, *qv, *cldfra, *o3;            // [nz*ncol]
  const rr *qc, *qi, *qs;                                   // nullptr: species absent
  const rr *re_cloud, *re_ice, *re_snow;                    // nullptr: has_req* = 0
  const rr* qr = nullptr;  // split3: rain (WSM3's rain / snow), [nz*ncol]
  bool split3 = false;         // RRTMG's "MP option 3": no F_QI, not warm_rain (LwColumnIn::split3)
};

struct SwBatchOutput {
  rr* rthraten;  // [nz*ncol] K/s (0 at night)
  rr* gsw;       // [ncol] W/m2
  // [8*ncol] (may be nullptr): swupt, swuptc, swdnt, swdntc, swupb, swupbc, swdnb, swdnbc (d*ncol + col)
  rr* diag = nullptr;
  rr* rthratenc = nullptr;  // [nz*ncol] clear-sky heating rate (may be nullptr)
  rr* cf = nullptr;         // [ncol] SW cloud forcing (may be nullptr)
  // [7*ncol] (may be nullptr): swvisdir, swvisdif, swnirdir, swnirdif, swddir, swddif, swddirc (d*ncol + col)
  rr* sfc = nullptr;
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
