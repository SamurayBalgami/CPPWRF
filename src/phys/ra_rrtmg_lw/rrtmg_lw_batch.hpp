// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "phys/ra_rrtmg_common/rrtmg_real.hpp"
// RRTMG_LWRAD for all columns of a domain on the device (CPU builds: OpenMP), the stages of
// rrtmg_lw_stages.hpp run as kernels over (column), (column, layer), (column, band) and
// (column, level), setcoef over (column, layer). All columns in one chunk where the device memory
// allows (rrtmg_chunk_columns), the work arrays shared with the other batches (batch_memory.hpp).
// Inputs and outputs are device arrays with the column index fastest: x[k*ncol + col].
#include <cstddef>
#include <memory>

#include "phys/ra_rrtmg_common/batch_memory.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_stages.hpp"

namespace cppwrf::phys {

struct LwBatchInput {
  int ncol;
  int col_begin = 0, col_end = -1;  // the columns computed: [col_begin, col_end) (-1: ncol); a tile's rows
  const rr *tsk, *emiss, *xland, *xice, *snow, *xlat;       // [ncol]
  const rr *dz8w;                                           // [nz*ncol] m
  const rr *p8w, *t8w;                                      // [(nz+1)*ncol] Pa, K
  const rr *p3d, *pi3d, *t3d, *qv, *cldfra, *o3;            // [nz*ncol]
  const rr *qc, *qi, *qs;                                   // nullptr: species absent (F_Q* false)
  const rr *re_cloud, *re_ice, *re_snow;                    // nullptr: has_req* = 0
  const rr* qr = nullptr;  // split3: rain (WSM3's rain / snow), [nz*ncol]
  bool split3 = false;         // RRTMG's "MP option 3": no F_QI, not warm_rain (LwColumnIn::split3)
};

struct LwBatchOutput {
  rr* rthraten;  // [nz*ncol] K/s
  rr* glw;       // [ncol] W/m2
  rr* olr;       // [ncol] (may be nullptr)
  // [8*ncol] (may be nullptr): lwupt, lwuptc, lwdnt, lwdntc, lwupb, lwupbc, lwdnb, lwdnbc (d*ncol + col)
  rr* diag = nullptr;
  rr* rthratenc = nullptr;  // [nz*ncol] clear-sky heating rate (may be nullptr)
  rr* cf = nullptr;         // [ncol] LW cloud forcing (may be nullptr)
  rr* o3 = nullptr;         // [nz*ncol] o3input != 2: the climatology ozone vmr (WRF's O33D) (may be nullptr)
};

class RrtmgLwBatch {
 public:
  // ncol: the domain's columns (they set the chunk size, rrtmg_chunk_columns)
  RrtmgLwBatch(const RrtmgLwTables& tables, int nlayers, int ncol);
  ~RrtmgLwBatch();
  RrtmgLwBatch(const RrtmgLwBatch&) = delete;
  RrtmgLwBatch& operator=(const RrtmgLwBatch&) = delete;
  void run(const RrtmgLwCallConfig& c, const LwBatchInput& in, const LwBatchOutput& out);

 private:
  RrtmgLwTables* d_tables_ = nullptr;  // shared device copies (rrtmg_shared_dev_copy)
  RrtmgLwConsts* d_consts_ = nullptr;
  LwLayout lay_;
  int chunk_;
  std::unique_ptr<RrtmgScratch> work_;  // the shared work arrays
};

}  // namespace cppwrf::phys
