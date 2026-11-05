#pragma once
// Device memory of the RRTMG batches (RrtmgLwBatch, RrtmgSwBatch), shared by every batch of the process: the
// batches of all domains run one after another on the one launch stream (LW, then SW; domain after domain), so
// one work buffer of the largest request serves them all, and a table set with the same contents is copied to
// the device once. The stages read the McICA cloud arrays only under the cloudy-layer / cloudy-column flags that
// the same call set, so what a previous batch left in the buffers does not enter the results (as before, when a
// domain's last, smaller chunk reused the arrays of its full chunks with another stride, and on the CPU, where a
// thread's workspace keeps the previous column's values).
#include <cstddef>

namespace cppwrf::phys {

// Columns per chunk for ncol columns of bytes_per_col work arrays each: all columns when they fit the budget (an
// eighth of the device memory and at most a quarter of what is free, at least 2048 columns or the shared buffer's
// present size), else the fewest chunks of equal size that do.
int rrtmg_chunk_columns(int ncol, std::size_t bytes_per_col);

// One batch's share of the work buffers (doubles and ints apart). The buffers grow to the largest request when a
// batch is constructed (outside the time loop) and are freed with the last batch; d() / i() are the current
// buffers, read at each run (a later batch may have moved them).
class RrtmgScratch {
 public:
  RrtmgScratch(std::size_t doubles, std::size_t ints);
  ~RrtmgScratch();
  RrtmgScratch(const RrtmgScratch&) = delete;
  RrtmgScratch& operator=(const RrtmgScratch&) = delete;
  double* d() const;
  int* i() const;
};

// Device copy of `bytes` of host data, shared by every caller with the same bytes (reference counted).
void* rrtmg_shared_dev_copy(const void* host, std::size_t bytes);
void rrtmg_shared_dev_release(void* dev);

}  // namespace cppwrf::phys
