#pragma once
// RRTMG_SWRAD for all columns of a domain on the device (CPU builds: OpenMP), the stages of
// rrtmg_sw_stages.hpp run as kernels over (column), (column, g-point) and (column, level), in
// chunks of columns. Arrays have the column index fastest: x[k*ncol + col].
#include <cstddef>

#include "phys/ra_rrtmg_sw/rrtmg_sw_stages.hpp"

namespace cppwrf::phys {

struct SwBatchInput {
  int ncol;
  const double *coszen, *albedo, *xland, *xice, *snow;          // [ncol]
  const double *p8w, *t8w;                                      // [(nz+1)*ncol] Pa, K
  const double *p3d, *pi3d, *t3d, *qv, *cldfra, *o3;            // [nz*ncol]
  const double *qc, *qi, *qs;                                   // nullptr: species absent
  const double *re_cloud, *re_ice, *re_snow;                    // nullptr: has_req* = 0
};

struct SwBatchOutput {
  double* rthraten;  // [nz*ncol] K/s (0 at night)
  double* gsw;       // [ncol] W/m2
};

class RrtmgSwBatch {
 public:
  RrtmgSwBatch(const RrtmgSwTables& tables, int nlayers, int max_chunk = 2048);
  ~RrtmgSwBatch();
  RrtmgSwBatch(const RrtmgSwBatch&) = delete;
  RrtmgSwBatch& operator=(const RrtmgSwBatch&) = delete;
  void run(const RrtmgSwCallConfig& c, const SwBatchInput& in, const SwBatchOutput& out);

 private:
  RrtmgSwTables* d_tables_ = nullptr;
  RrtmgLwConsts* d_consts_ = nullptr;
  SwLayout lay_;
  int chunk_;
  double* d_ = nullptr;
  int* i_ = nullptr;
};

}  // namespace cppwrf::phys
