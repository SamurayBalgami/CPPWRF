// SPDX-License-Identifier: GPL-3.0-or-later
#include "phys/ra_rrtmg_common/batch_memory.hpp"

#include <algorithm>
#include <cstring>
#include <mutex>
#include <vector>

#include "core/backend.hpp"

#if defined(CPPWRF_USE_CUDA)
#include <cuda_runtime.h>
#endif

namespace cppwrf::phys {

namespace {

std::mutex mtx;

struct Buffer {
  void* p = nullptr;
  std::size_t bytes = 0;
  void grow(std::size_t b) {
    if (b <= bytes) return;
    dev_free(p);
    p = dev_alloc(b);
    bytes = b;
  }
  void release() {
    dev_free(p);
    p = nullptr;
    bytes = 0;
  }
};
Buffer dbuf, ibuf;
int users = 0;

struct SharedCopy {
  std::vector<unsigned char> host;
  void* dev;
  int refs;
};
std::vector<SharedCopy> copies;

}  // namespace

int rrtmg_chunk_columns(int ncol, std::size_t bytes_per_col) {
  ncol = std::max(ncol, 1);
  bytes_per_col = std::max<std::size_t>(bytes_per_col, 1);
  std::size_t budget = 0;
#if defined(CPPWRF_USE_CUDA)
  std::size_t free_b = 0, total_b = 0;
  if (cudaMemGetInfo(&free_b, &total_b) == cudaSuccess) budget = std::min(total_b / 8, free_b / 4);
#endif
  {
    std::lock_guard<std::mutex> lk(mtx);
    budget = std::max(budget, dbuf.bytes + ibuf.bytes);
  }
  const std::size_t cap = std::max<std::size_t>(budget / bytes_per_col, 2048);
  if (std::size_t(ncol) <= cap) return ncol;
  const std::size_t nchunks = (std::size_t(ncol) + cap - 1) / cap;
  return int((std::size_t(ncol) + nchunks - 1) / nchunks);
}

RrtmgScratch::RrtmgScratch(std::size_t doubles, std::size_t ints) {
  std::lock_guard<std::mutex> lk(mtx);
  dbuf.grow(std::max<std::size_t>(doubles, 1) * sizeof(double));
  ibuf.grow(std::max<std::size_t>(ints, 1) * sizeof(int));
  ++users;
}

RrtmgScratch::~RrtmgScratch() {
  std::lock_guard<std::mutex> lk(mtx);
  if (--users == 0) {
    dbuf.release();
    ibuf.release();
  }
}

double* RrtmgScratch::d() const {
  std::lock_guard<std::mutex> lk(mtx);
  return static_cast<double*>(dbuf.p);
}

int* RrtmgScratch::i() const {
  std::lock_guard<std::mutex> lk(mtx);
  return static_cast<int*>(ibuf.p);
}

void* rrtmg_shared_dev_copy(const void* host, std::size_t bytes) {
  std::lock_guard<std::mutex> lk(mtx);
  const auto* h = static_cast<const unsigned char*>(host);
  for (SharedCopy& c : copies) {
    if (c.host.size() == bytes && std::memcmp(c.host.data(), h, bytes) == 0) {
      ++c.refs;
      return c.dev;
    }
  }
  void* d = dev_alloc(bytes);
  copy_to_dev(d, host, bytes);
  copies.push_back(SharedCopy{std::vector<unsigned char>(h, h + bytes), d, 1});
  return d;
}

void rrtmg_shared_dev_release(void* dev) {
  if (!dev) return;
  std::lock_guard<std::mutex> lk(mtx);
  for (std::size_t k = 0; k < copies.size(); ++k) {
    if (copies[k].dev != dev) continue;
    if (--copies[k].refs == 0) {
      dev_free(copies[k].dev);
      copies.erase(copies.begin() + std::ptrdiff_t(k));
    }
    return;
  }
}

}  // namespace cppwrf::phys
