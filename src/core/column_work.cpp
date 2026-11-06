#include "core/column_work.hpp"

#include <algorithm>
#include <mutex>
#include <stdexcept>
#include <string>

#if defined(CPPWRF_USE_CUDA)
#include <cuda_runtime.h>
#endif

namespace cppwrf {

#if defined(CPPWRF_USE_CUDA)
namespace {
std::mutex mtx;
char* buf = nullptr;
std::size_t buf_bytes = 0;
// grows the buffer to `bytes` (its contents are not kept: every kernel writes its arrays before it reads them)
void grow(std::size_t bytes) {
  if (bytes <= buf_bytes) return;
  dev_free(buf);
  buf = nullptr;
  buf_bytes = 0;
  try {
    buf = static_cast<char*>(dev_alloc(bytes));
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("the column physics' work arrays: ") + e.what());
  }
  buf_bytes = bytes;
}
}  // namespace

void column_work_reserve(const Range2D& r, std::size_t bytes_per_col) {
  if (r.i1 <= r.i0 || r.j1 <= r.j0 || bytes_per_col == 0) return;
  const int ni = r.i1 - r.i0, nj = r.j1 - r.j0;
  const std::size_t bpc = detail::column_work_bpc(bytes_per_col);
  std::size_t budget = 0, free_b = 0, total_b = 0;
  if (cudaMemGetInfo(&free_b, &total_b) == cudaSuccess) budget = std::min(total_b / 8, free_b / 4);
  std::size_t want = detail::column_for_threads(ni, nj) * bpc;
  if (want > budget) {  // whole block rows (a block covers column_block().y rows), at least one
    const std::size_t pair = detail::column_for_threads(ni, column_block().y) * bpc;
    want = std::max<std::size_t>(budget / pair, 1) * pair;
  }
  std::lock_guard<std::mutex> lk(mtx);
  grow(want);
}

namespace detail {
char* column_work_buffer(std::size_t bytes, std::size_t* have) {
  std::lock_guard<std::mutex> lk(mtx);
  grow(bytes);
  *have = buf_bytes;
  return buf;
}
}  // namespace detail
#else
void column_work_reserve(const Range2D&, std::size_t) {}
namespace detail {
char* column_work_buffer(std::size_t, std::size_t* have) {
  *have = 0;
  return nullptr;
}
}  // namespace detail
#endif

}  // namespace cppwrf
