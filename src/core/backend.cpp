#include "core/backend.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <map>
#include <mutex>
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
void* cuda_stream_raw() {
  static cudaStream_t s = [] {
    cudaStream_t t = nullptr;
    ck(cudaStreamCreate(&t), "cudaStreamCreate");  // blocking: orders with any work on the legacy stream 0
    return t;
  }();
  return s;
}
static cudaStream_t stream() { return static_cast<cudaStream_t>(cuda_stream_raw()); }

// ---- time step graphs (StepScope) and the fatal flags of the step ----
namespace {
int graphs_mode() {  // 0 off, 1 on, 2 on with the counts at exit
  static const int m = [] {
#ifdef CPPWRF_CUDA_SYNC_CHECK
    return 0;  // a device synchronisation after every launch cannot be captured
#endif
    if (ktime::on()) return 0;
    const char* e = std::getenv("CPPWRF_GRAPHS");
    if (!e || !*e) return 1;
    if (std::strcmp(e, "stats") == 0) return 2;
    return *e == '0' ? 0 : 1;
  }();
  return m;
}
struct Poll {
  int* slot;
  void (*reset)();
};
constexpr int kSlots = 256;  // fatal checks per segment
struct Graphs {
  bool in_step = false, capturing = false;
  const void* key = nullptr;
  int seg = 0;
  std::map<std::pair<const void*, int>, std::vector<cudaGraphExec_t>> cache;  // most recently used first
  int* slots = nullptr;  // pinned, 2 banks of kSlots
  int bank = 0;
  std::vector<Poll> cur, pend;  // checks of the open segment; of the launched, unchecked ones
  cudaEvent_t ev = nullptr;
  bool ev_pending = false;
  long n_launch = 0, n_update = 0, n_inst = 0, n_seg_nodes = 0;
  void init() {
    if (slots) return;
    ck(cudaHostAlloc(reinterpret_cast<void**>(&slots), 2 * kSlots * sizeof(int), cudaHostAllocDefault), "cudaHostAlloc");
    ck(cudaEventCreateWithFlags(&ev, cudaEventDisableTiming), "cudaEventCreate");
  }
  void begin_capture() {
    ck(cudaStreamBeginCapture(stream(), cudaStreamCaptureModeThreadLocal), "cudaStreamBeginCapture");
    capturing = true;
  }
  // the instantiated graph of this domain and segment, updated to g's launch parameters (a step of another
  // kind has another topology: up to 4 per segment are kept)
  cudaGraphExec_t exec_for(cudaGraph_t g) {
    auto& v = cache[{key, seg}];
    for (std::size_t i = 0; i < v.size(); ++i) {
      cudaGraphExecUpdateResultInfo info;
      if (cudaGraphExecUpdate(v[i], g, &info) == cudaSuccess) {
        ++n_update;
        std::rotate(v.begin(), v.begin() + long(i), v.begin() + long(i) + 1);
        return v[0];
      }
      (void)cudaGetLastError();  // a failed update is no launch error
    }
    cudaGraphExec_t e = nullptr;
    ck(cudaGraphInstantiate(&e, g, 0), "cudaGraphInstantiate");
    ++n_inst;
    if (v.size() >= 4) {
      cudaGraphExecDestroy(v.back());
      v.pop_back();
    }
    v.insert(v.begin(), e);
    return e;
  }
  // ends the open segment: launches its graph (after the previous segment's flags were read) and records the
  // event that its fatal checks wait for
  void segment_point() {
    bool launched = false;
    if (capturing) {
      cudaGraph_t g = nullptr;
      capturing = false;
      ck(cudaStreamEndCapture(stream(), &g), "cudaStreamEndCapture");
      std::size_t nn = 0;
      ck(cudaGraphGetNodes(g, nullptr, &nn), "cudaGraphGetNodes");
      cudaGraphExec_t e = nn ? exec_for(g) : nullptr;
      cudaGraphDestroy(g);
      flush();  // the previous step or segment is done and its flags are clean before this one runs
      if (e) {
        ck(cudaGraphLaunch(e, stream()), "cudaGraphLaunch");
        launched = true;
        ++n_launch;
        n_seg_nodes += long(nn);
      }
      ++seg;
    } else if (!cur.empty()) {
      flush();
    }
    if (launched || !cur.empty()) {  // at most one step or segment in flight unchecked
      pend.swap(cur);
      cur.clear();
      bank ^= 1;
      ck(cudaEventRecord(ev, stream()), "cudaEventRecord");
      ev_pending = true;
    }
  }
  // reads the flags of the launched checks in the order of the calls; the first raised code stops the run
  void flush() {
    if (!ev_pending) return;
    ck(cudaEventSynchronize(ev), "cudaEventSynchronize");
    ev_pending = false;
    int code = 0;
    for (const Poll& p : pend)
      if (*p.slot != 0) {
        if (code == 0) code = *p.slot;
        p.reset();
      }
    pend.clear();
    if (code != 0) fatal_stop(code);
  }
  void abort() {  // an exception inside the step
    if (capturing) {
      cudaGraph_t g = nullptr;
      cudaStreamEndCapture(stream(), &g);
      if (g) cudaGraphDestroy(g);
      (void)cudaGetLastError();
      capturing = false;
    }
    cur.clear();
    in_step = false;
  }
};
Graphs& graphs() {  // never destroyed: device memory may be freed by static objects at exit
  static Graphs* gs = new Graphs;
  return *gs;
}
struct GraphStats {
  ~GraphStats() {
    const Graphs& g = graphs();
    if (graphs_mode() == 2 && g.n_launch > 0)
      std::fprintf(stderr, "cppwrf graphs: %ld launches, %ld updates, %ld instantiations, %.1f nodes per launch\n",
                   g.n_launch, g.n_update, g.n_inst, double(g.n_seg_nodes) / double(g.n_launch));
  }
} graph_stats;
// an operation that waits for the device: the open segment is launched first, the capture resumes after it
void sync_begin() { graphs().segment_point(); }
void sync_end() {
  Graphs& gs = graphs();
  gs.flush();
  if (gs.in_step && graphs_mode() != 0 && !gs.capturing) gs.begin_capture();
}
}  // namespace

bool step_graphs_on() { return graphs_mode() != 0; }

StepScope::StepScope(const void* key) {
  Graphs& gs = graphs();
  if (gs.in_step) throw std::logic_error("StepScope: a step inside a step");
  gs.init();
  gs.in_step = true;
  gs.key = key;
  gs.seg = 0;
  open_ = true;
  if (graphs_mode() != 0) gs.begin_capture();
}
StepScope::~StepScope() {
  if (open_) graphs().abort();
}
void StepScope::end() {
  Graphs& gs = graphs();
  open_ = false;
  gs.in_step = false;
  gs.segment_point();
}

namespace detail {
int* fatal_poll_slot(void (*reset)()) {
  Graphs& gs = graphs();
  if (!gs.in_step) return nullptr;
  // a stop raised by host code (a kernel run on the host): the checks before this one come first, as before
  if (int hc = __atomic_exchange_n(&fatal_code_host, 0, __ATOMIC_RELAXED)) {
    sync_begin();
    gs.flush();
    fatal_stop(hc);
  }
  if (gs.cur.size() >= std::size_t(kSlots)) throw std::runtime_error("fatal_check: too many checks in one step");
  int* slot = gs.slots + gs.bank * kSlots + int(gs.cur.size());
  gs.cur.push_back({slot, reset});
  return slot;
}
void fatal_flush() { graphs().flush(); }

void reserve_local_memory(const void* kernel, const char* kernel_type) {
  cudaFuncAttributes a{};
  ck(cudaFuncGetAttributes(&a, kernel), "cudaFuncGetAttributes");
  static std::mutex m;
  static std::size_t reserved = 0;  // the largest frame of the process's launches (one device)
  std::lock_guard<std::mutex> lk(m);
  if (a.localSizeBytes <= reserved) return;
  std::size_t limit = 0;
  ck(cudaDeviceGetLimit(&limit, cudaLimitStackSize), "cudaDeviceGetLimit");
  if (a.localSizeBytes <= limit) {
    reserved = a.localSizeBytes;
    return;
  }
  int dev = 0, sms = 0, tpsm = 0;
  ck(cudaGetDevice(&dev), "cudaGetDevice");
  ck(cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, dev), "cudaDeviceGetAttribute");
  ck(cudaDeviceGetAttribute(&tpsm, cudaDevAttrMaxThreadsPerMultiProcessor, dev), "cudaDeviceGetAttribute");
  sync_begin();  // inside a step: the captured part runs first (as for an allocation)
  std::size_t free_b = 0, total_b = 0;
  const cudaError_t mi = cudaMemGetInfo(&free_b, &total_b);
  const cudaError_t e = cudaDeviceSetLimit(cudaLimitStackSize, a.localSizeBytes);
  if (e != cudaSuccess) (void)cudaGetLastError();
  sync_end();
  if (e != cudaSuccess) {
    int st = 0;
    char* dm = abi::__cxa_demangle(kernel_type, nullptr, nullptr, &st);
    const std::string name = (st == 0 && dm) ? dm : kernel_type;
    std::free(dm);
    const double mib = 1024. * 1024.;
    char msg[640];
    std::snprintf(msg, sizeof msg,
                  "device memory too small: the kernel %.200s keeps %zu bytes per thread in local memory, which the GPU "
                  "reserves for each of its %d x %d resident threads (%.0f MiB in all, %.0f MiB of it reserved so far); %.0f MiB of "
                  "%.0f MiB are free (%s). Use a smaller domain, fewer levels or a GPU with more memory.",
                  name.c_str(), std::size_t(a.localSizeBytes), sms, tpsm,
                  double(a.localSizeBytes) * sms * tpsm / mib, double(limit) * sms * tpsm / mib,
                  mi == cudaSuccess ? double(free_b) / mib : 0., mi == cudaSuccess ? double(total_b) / mib : 0.,
                  cudaGetErrorString(e));
    throw std::runtime_error(msg);
  }
  reserved = a.localSizeBytes;
}
}  // namespace detail

void* dev_alloc(std::size_t bytes) {
  sync_begin();
  void* p = nullptr;
  if (const cudaError_t e = cudaMalloc(&p, bytes); e != cudaSuccess) {
    (void)cudaGetLastError();
    std::size_t free_b = 0, total_b = 0;
    cudaMemGetInfo(&free_b, &total_b);
    const double mib = 1024. * 1024.;
    char msg[320];
    std::snprintf(msg, sizeof msg,
                  "device memory too small: an allocation of %.0f MiB failed (%s); %.0f MiB of %.0f MiB are free. Use "
                  "a smaller domain, fewer levels or a GPU with more memory.",
                  double(bytes) / mib, cudaGetErrorString(e), double(free_b) / mib, double(total_b) / mib);
    throw std::runtime_error(msg);
  }
  ck(cudaMemsetAsync(p, 0, bytes, stream()), "cudaMemset");
  ck(cudaStreamSynchronize(stream()), "cudaMemset");
  sync_end();
  return p;
}
void dev_free(void* p) {
  if (!p) return;
  sync_begin();
  cudaFree(p);
  sync_end();
}
void copy_to_host(void* d, const void* s, std::size_t n) {
  sync_begin();
  ck(cudaMemcpyAsync(d, s, n, cudaMemcpyDeviceToHost, stream()), "D2H");
  ck(cudaStreamSynchronize(stream()), "D2H");
  sync_end();
}
void copy_to_dev(void* d, const void* s, std::size_t n) {
  sync_begin();
  ck(cudaMemcpyAsync(d, s, n, cudaMemcpyHostToDevice, stream()), "H2D");
  ck(cudaStreamSynchronize(stream()), "H2D");
  sync_end();
}
void copy_dev(void* d, const void* s, std::size_t n) { ck(cudaMemcpyAsync(d, s, n, cudaMemcpyDeviceToDevice, stream()), "D2D"); }
void dev_sync() {
  sync_begin();
  ck(cudaDeviceSynchronize(), "sync");
  sync_end();
}
const char* backend_name() { return "cuda"; }
void device_memory_report(const char* when) {
  static const bool on = [] { const char* e = std::getenv("CPPWRF_MEMINFO"); return e && *e && *e != '0'; }();
  if (!on) return;
  dev_sync();
  std::size_t free_b = 0, total_b = 0, limit = 0;
  int dev = 0, sms = 0, tpsm = 0;
  ck(cudaMemGetInfo(&free_b, &total_b), "cudaMemGetInfo");
  ck(cudaDeviceGetLimit(&limit, cudaLimitStackSize), "cudaDeviceGetLimit");
  ck(cudaGetDevice(&dev), "cudaGetDevice");
  ck(cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, dev), "cudaDeviceGetAttribute");
  ck(cudaDeviceGetAttribute(&tpsm, cudaDevAttrMaxThreadsPerMultiProcessor, dev), "cudaDeviceGetAttribute");
  const double mib = 1024. * 1024.;
  std::printf("cppwrf: device memory %s: %.0f MiB used of %.0f MiB; stack limit %zu bytes x %d resident threads "
              "(%.0f MiB)\n",
              when, double(total_b - free_b) / mib, double(total_b) / mib, limit, sms * tpsm,
              double(limit) * sms * tpsm / mib);
}
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
void device_memory_report(const char*) {}
StepScope::StepScope(const void*) {}
StepScope::~StepScope() {}
void StepScope::end() {}
bool step_graphs_on() { return false; }
Reducer::Reducer(int max_k) : max_k_(max_k), hpart_(max_k) {}
Reducer::~Reducer() {}
#endif

// One warp per block: a domain of ~10^4 columns gives ~330 blocks, every SM of a large GPU gets some (128-thread
// blocks left more than half of the RTX PRO 6000's 188 SMs idle); one thread per column either way
ColumnBlock column_block() {
  static const ColumnBlock b = [] {
    ColumnBlock r{32, 1};
    if (const char* e = std::getenv("CPPWRF_COLUMN_BLOCK")) {
      int x = 0, y = 0;
      if (std::sscanf(e, "%dx%d", &x, &y) == 2 && x > 0 && x % 32 == 0 && y > 0 && x * y <= 1024) r = {x, y};
      else throw std::runtime_error(std::string("CPPWRF_COLUMN_BLOCK=") + e + ": XxY with X a multiple of 32");
    }
    return r;
  }();
  return b;
}

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
    case kFatalNoahUrbanLai:
      return "rdlai2d: an urban point with VEGFRA = 0 and LAI <= 0 is not supported (WRF's SFLX sets SHDFAC = 0.05, "
             "skips CANRES for XLAI = 0 and uses the plant coefficient PC of the previous land point, an undefined value)";
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
