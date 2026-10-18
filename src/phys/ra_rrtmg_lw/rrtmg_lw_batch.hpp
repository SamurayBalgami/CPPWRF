#pragma once
// RRTMG_LWRAD for all columns of a domain on the device (CPU builds: OpenMP), the stages of
// rrtmg_lw_stages.hpp run as kernels over (column), (column, layer), (column, band) and
// (column, level). Columns are processed in chunks so that the work arrays stay below ~1 GB.
// Inputs and outputs are device arrays with the column index fastest: x[k*ncol + col].
#include <cstddef>

#include "phys/ra_rrtmg_lw/rrtmg_lw_stages.hpp"

namespace cppwrf::phys {

struct LwBatchInput {
  int ncol;
  const double *tsk, *emiss, *xland, *xice, *snow;              // [ncol]
  const double *p8w, *t8w;                                      // [(nz+1)*ncol] Pa, K
  const double *p3d, *pi3d, *t3d, *qv, *cldfra, *o3;            // [nz*ncol]
  const double *qc, *qi, *qs;                                   // nullptr: species absent (F_Q* false)
  const double *re_cloud, *re_ice, *re_snow;                    // nullptr: has_req* = 0
};

struct LwBatchOutput {
  double* rthraten;  // [nz*ncol] K/s
  double* glw;       // [ncol] W/m2
  double* olr;       // [ncol] (may be nullptr)
  // [8*ncol] (may be nullptr): lwupt, lwuptc, lwdnt, lwdntc, lwupb, lwupbc, lwdnb, lwdnbc (d*ncol + col)
  double* diag = nullptr;
  double* rthratenc = nullptr;  // [nz*ncol] clear-sky heating rate (may be nullptr)
  double* cf = nullptr;         // [ncol] LW cloud forcing (may be nullptr)
};

class RrtmgLwBatch {
 public:
  RrtmgLwBatch(const RrtmgLwTables& tables, int nlayers, int max_chunk = 4096);
  ~RrtmgLwBatch();
  RrtmgLwBatch(const RrtmgLwBatch&) = delete;
  RrtmgLwBatch& operator=(const RrtmgLwBatch&) = delete;
  void run(const RrtmgLwCallConfig& c, const LwBatchInput& in, const LwBatchOutput& out);

 private:
  RrtmgLwTables* d_tables_ = nullptr;
  RrtmgLwConsts* d_consts_ = nullptr;
  LwLayout lay_;
  int chunk_;
  double* d_ = nullptr;
  int* i_ = nullptr;
};

}  // namespace cppwrf::phys
