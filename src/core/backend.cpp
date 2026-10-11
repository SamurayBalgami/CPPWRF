#include "core/backend.hpp"

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
void* dev_alloc(std::size_t bytes) {
  std::size_t n = (bytes + 63) / 64 * 64;
  void* p = std::aligned_alloc(64, n ? n : 64);
  if (!p) throw std::runtime_error("allocation failed");
  std::memset(p, 0, n);
  return p;
}
void dev_free(void* p) { std::free(p); }
void copy_to_host(void* d, const void* s, std::size_t n) { std::memcpy(d, s, n); }
void copy_to_dev(void* d, const void* s, std::size_t n) { std::memcpy(d, s, n); }
void copy_dev(void* d, const void* s, std::size_t n) { std::memcpy(d, s, n); }
void dev_sync() {}
const char* backend_name() { return "cpu"; }
Reducer::Reducer(int max_k) : max_k_(max_k), hpart_(max_k) {}
Reducer::~Reducer() {}
#endif

}  // namespace cppwrf
