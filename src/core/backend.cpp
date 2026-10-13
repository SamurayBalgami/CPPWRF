#include "core/backend.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>

#if defined(CPPWRF_USE_CUDA)
#include <cuda_runtime.h>
#endif

namespace cppwrf {

#if defined(CPPWRF_USE_CUDA)
static void ck(cudaError_t e, const char* what) {
  if (e != cudaSuccess) throw std::runtime_error(std::string("CUDA error in ") + what + ": " + cudaGetErrorString(e));
}
void cuda_check_launch(const char* where) {
  ck(cudaGetLastError(), where);
#ifdef CPPWRF_CUDA_SYNC_CHECK
  ck(cudaDeviceSynchronize(), where);
#endif
}
void* dev_alloc(std::size_t bytes) {
  void* p = nullptr;
  ck(cudaMalloc(&p, bytes), "cudaMalloc");
  ck(cudaMemset(p, 0, bytes), "cudaMemset");
  return p;
}
void dev_free(void* p) { if (p) cudaFree(p); }
void copy_to_host(void* d, const void* s, std::size_t n) { ck(cudaMemcpy(d, s, n, cudaMemcpyDeviceToHost), "D2H"); }
void copy_to_dev(void* d, const void* s, std::size_t n) { ck(cudaMemcpy(d, s, n, cudaMemcpyHostToDevice), "H2D"); }
void copy_dev(void* d, const void* s, std::size_t n) { ck(cudaMemcpy(d, s, n, cudaMemcpyDeviceToDevice), "D2D"); }
void dev_sync() { ck(cudaDeviceSynchronize(), "sync"); }
const char* backend_name() { return "cuda"; }
Reducer::Reducer(int max_k) : max_k_(max_k), hpart_(max_k) { dpart_ = static_cast<double*>(dev_alloc(sizeof(double) * max_k)); }
Reducer::~Reducer() { dev_free(dpart_); }
#else
namespace {
// Threaded memset/memcpy in contiguous static chunks: the same thread-to-memory mapping as the
// static parallel_for over (k, j), so first-touch places pages on the NUMA node that uses them
// (review X2), and the per-step field copies of RK3 use all cores.
constexpr std::size_t kChunk = std::size_t(1) << 16;
void par_set0(void* p, std::size_t n) {
  const long nc = long((n + kChunk - 1) / kChunk);
#pragma omp parallel for schedule(static)
  for (long c = 0; c < nc; ++c) {
    const std::size_t o = std::size_t(c) * kChunk;
    std::memset(static_cast<char*>(p) + o, 0, std::min(kChunk, n - o));
  }
}
void par_copy(void* d, const void* s, std::size_t n) {
  const long nc = long((n + kChunk - 1) / kChunk);
#pragma omp parallel for schedule(static)
  for (long c = 0; c < nc; ++c) {
    const std::size_t o = std::size_t(c) * kChunk;
    std::memcpy(static_cast<char*>(d) + o, static_cast<const char*>(s) + o, std::min(kChunk, n - o));
  }
}
}  // namespace
void* dev_alloc(std::size_t bytes) {
  std::size_t n = (bytes + 63) / 64 * 64;
  void* p = std::aligned_alloc(64, n ? n : 64);
  if (!p) throw std::runtime_error("allocation failed");
  par_set0(p, n);
  return p;
}
void dev_free(void* p) { std::free(p); }
void copy_to_host(void* d, const void* s, std::size_t n) { par_copy(d, s, n); }
void copy_to_dev(void* d, const void* s, std::size_t n) { par_copy(d, s, n); }
void copy_dev(void* d, const void* s, std::size_t n) { par_copy(d, s, n); }
void dev_sync() {}
const char* backend_name() { return "cpu"; }
Reducer::Reducer(int max_k) : max_k_(max_k), hpart_(max_k) {}
Reducer::~Reducer() {}
#endif

}  // namespace cppwrf
