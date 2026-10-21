#include "core/backend.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <map>
#include <string>
#include <vector>

#if defined(CPPWRF_USE_CUDA)
#include <cuda_runtime.h>
#include <cxxabi.h>
#endif

namespace cppwrf {

#if defined(CPPWRF_USE_CUDA)
static void ck(cudaError_t e, const char* what) {
  if (e != cudaSuccess) throw std::runtime_error(std::string("CUDA error in ") + what + ": " + cudaGetErrorString(e));
}
namespace ktime {
namespace {
struct Table {
  std::map<std::string, std::pair<double, long>> t;  // kernel type -> (seconds, launches)
  std::chrono::steady_clock::time_point t0;
  ~Table() {
    if (t.empty()) return;
    std::vector<std::pair<std::string, std::pair<double, long>>> v(t.begin(), t.end());
    std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) { return a.second.first > b.second.first; });
    double tot = 0;
    for (auto& e : v) tot += e.second.first;
    std::fprintf(stderr, "cppwrf ktime: %.3f s in kernels\n", tot);
    for (auto& e : v) {
      int st = 0;
      char* dm = abi::__cxa_demangle(e.first.c_str(), nullptr, nullptr, &st);
      std::string name = (st == 0 && dm) ? dm : e.first;
      std::free(dm);
      if (name.size() > 110) name = name.substr(0, 110);
      std::fprintf(stderr, "%6.2f%% %9.2f ms %7ld  %s\n", 100 * e.second.first / tot, 1e3 * e.second.first,
                   e.second.second, name.c_str());
    }
  }
};
Table& table() { static Table tb; return tb; }
}  // namespace
bool on() {
  static const bool enabled = [] { const char* e = std::getenv("CPPWRF_KTIME"); return e && *e && *e != '0'; }();
  return enabled;
}
void begin() {
  ck(cudaDeviceSynchronize(), "ktime");
  table().t0 = std::chrono::steady_clock::now();
}
void end(const char* kernel_type) {
  ck(cudaDeviceSynchronize(), "ktime");
  auto& e = table().t[kernel_type];
  e.first += std::chrono::duration<double>(std::chrono::steady_clock::now() - table().t0).count();
  e.second += 1;
}
}  // namespace ktime

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
void copy_dev(void* d, const void* s, std::size_t n) { ck(cudaMemcpyAsync(d, s, n, cudaMemcpyDeviceToDevice, 0), "D2D"); }
void dev_sync() { ck(cudaDeviceSynchronize(), "sync"); }
const char* backend_name() { return "cuda"; }
Reducer::Reducer(int max_k) : max_k_(max_k), hpart_(max_k) { dpart_ = static_cast<double*>(dev_alloc(sizeof(double) * max_k)); }
Reducer::~Reducer() { dev_free(dpart_); }
#else
namespace {
// Threaded memset/memcpy in contiguous static chunks: the same thread-to-memory mapping as the
// static parallel_for over (k, j), so first-touch places pages on the NUMA node that uses them
// and the per-step field copies of RK3 use all cores.
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

int fatal_code_host = 0;

const char* fatal_message(int code) {
  switch (code) {
    case kFatalSnowDensity: return "Physical snow depth is less than snow water equiv.";
    case kFatalKfQgSurface: return "QG, QG(NK).LT.0";
    case kFatalKfIstop: return "KAIN-FRITSCH, istop=1, diags";
    case kFatalNoahSoilTypes: return "Warning: too many input soil types";
    case kFatalNoahLanduseTypes: return "Warning: too many input landuse types";
    case kFatalNoahSlopeTypes: return "Warning: too many input slope types";
    case kFatalNoahRootLayers: return "Error: too many root layers";
    case kFatalNoahmpWater: return "Error: Water budget problem in NoahMP LSM";
    case kFatalNoahmpSolar: return "Error: Solar radiation budget problem in NoahMP LSM";
    case kFatalNoahmpEnergy: return "Error: Energy budget problem in NoahMP LSM";
    case kFatalNoahmpLongwave: return "Error: Longwave radiation budget problem in NoahMP LSM";
    case kFatalNoahmpZeroPlane: return "Error: ZeroPlaneDisp problem in NoahMP LSM";
    case kFatalNoahmpRefHeightCanopy: return "Error in ResistanceAboveCanopyMostMod.F90";
    case kFatalNoahmpRefHeightBare: return "Error in ResistanceBareGroundMostMod.F90";
    case kFatalNoahmpGlacierEnergy: return "Error: Surface energy budget problem in NoahMP LSM";
    default: return "unknown fatal condition in a physics kernel";
  }
}

void fatal_stop(int code) {
  // wrf_error_fatal: the block on stdout (wrf_message) and on stderr (wrf_message2)
  const char* msg = fatal_message(code);
  for (std::FILE* f : {stdout, stderr}) {
    std::fprintf(f, "-------------- FATAL CALLED ---------------\n%s\n-------------------------------------------\n", msg);
    std::fflush(f);
  }
  throw FatalError(code, msg);
}

}  // namespace cppwrf
