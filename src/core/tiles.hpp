// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Domain decomposition (docs/MPI.md). A domain runs as P x Q patches; a patch is run by a worker: a tile (a host
// thread of this process, run.tiles / CPPWRF_TILES) or an MPI rank (a CPPWRF_MPI build under mpirun). Each worker runs
// the whole program with its own fields, allocated for its patch and a halo of Grid::halo points and indexed with the
// domain's indices (Grid::dom, pi0 .. pj1; Layout::wi0 .. wj1), on its own stream and device, and computes only the
// points it owns: every launch over (i, j) is clipped to the patch (core/backend.hpp, tile::ctx). The workers meet in
// collective calls, which every worker makes in the same order with the same arguments:
//   - the halo exchange (every core/halo.hpp fill, and tile::exchange where a kernel reads neighbours that no fill
//     refreshes): a worker receives the points the other workers own inside its allocation;
//   - the wrap exchange of a periodic direction with several patches along it (inside the fills);
//   - reductions (the maximum, sums in the order of a run with one patch) and gathers for the one writer of the output.
// The collectives only move values and every point is computed by one worker with the inputs of a run with one patch,
// so the result does not depend on the number of patches.
#include <functional>
#include <mutex>
#include <string>
#include <vector>

#include "core/field.hpp"

namespace cppwrf::tile {

struct Options {
  int px = 1, py = 1;        // patches along x and y
  std::vector<int> devices;  // CUDA device of tile t: devices[t % size] (empty: the current device)
  bool check = false;        // compare every launch of a step with a run with one patch (tiles, CPU build; slow)
  int count() const { return px * py; }
  bool on() const { return count() > 1 || check; }
};
// run.tiles ("P", "PxQ"; the environment CPPWRF_TILES overrides it), CPPWRF_TILE_DEVICES=d0,d1,... and
// CPPWRF_TILE_CHECK=1 (the reference run with one patch is compared with every tile after every launch of a time
// step: the first point that differs names the launch and the field)
Options options(const std::string& run_tiles);

// Tiles: body(t) for t = 0 .. count - 1, each on its own host thread (and CUDA device), and with o.check also
// body(-1), the reference run. The first error stops every tile and is rethrown here; else body(0)'s value.
int run(const Options& o, const std::function<int(int)>& body);

// MPI ranks (a CPPWRF_MPI build; otherwise mpi_start returns false and the others are not called).
// mpi_start: MPI_Init; true when this process is one of several ranks. mpi_run: the calling rank as the worker of
// patch `rank` of the layout (run.ranks / CPPWRF_RANKS = PxQ, default the most square P x Q with P Q = the rank
// count, P <= Q as WRF's): body(rank); an error aborts every rank. mpi_finish: MPI_Finalize.
bool mpi_start(int* argc, char*** argv);
int mpi_run(const std::string& run_ranks, const std::function<int(int)>& body);
void mpi_finish();
void mpi_abort(int code);  // an error outside mpi_run: MPI_Abort (no-op without MPI)

int id();       // the calling worker's patch (0 without a decomposition, -1 for the check's reference run)
int count();    // the workers (1 without a decomposition)
bool tiled();   // the calling thread is a worker of a decomposed run
// the calling worker writes the output files and prints the log: no decomposition, patch 0, the reference run
bool writer();

// The calling worker's patch of a new domain, before its first field: sets g's patch (Grid::dom, pi0 .. pj1) and
// makes it the patch that clips the launches. Patch t = tx + px ty owns the mass points [tx nx / px, (tx + 1) nx / px)
// x [ty ny / py, (ty + 1) ny / py) and the faces with the same indices; the edge patches also own the last faces and
// the lateral halos beyond the edge. No decomposition (and the reference run): g is left as it is.
void set_domain(Grid& g);
// The launches in this scope are clipped to g's patch (nested runs: the domain whose points the launches index).
class Use {
 public:
  explicit Use(const Grid& g);
  ~Use();
  Use(const Use&) = delete;
  Use& operator=(const Use&) = delete;

 private:
  Ctx saved_;
};
// patch t of decomposition dom: i0, i1, j0, j1 (owned mass points)
struct PatchBox {
  int i0, i1, j0, j1;
};
PatchBox patch(int dom, int t);
// The mass points of a domain of nx x ny points the calling worker's fields hold (its patch and Grid::halo points
// around it, inside the domain), before the domain's Grid exists (the input readers); {0, nx, 0, ny} without a
// decomposition.
PatchBox input_window(int nx, int ny);
int patches_x(), patches_y();

// A time step: the check's comparisons on.
class Step {
 public:
  Step();
  ~Step();
  Step(const Step&) = delete;
  Step& operator=(const Step&) = delete;
};

// ---- collectives for the host code (no-ops / the plain operation without a decomposition) ----
// Every field's halo from the other workers (after host code that set fields point by point: the start, an input)
void exchange_all();
// The whole field on the writer, in the layout of a run with one patch (whole; halo points 0): h is empty on the
// other workers. Without a decomposition: the field's download.
void download(const void* base, std::vector<char>& h, Layout& whole);
template <class F>
void download(const F& f, std::vector<real>& h, Layout& whole) {
  if (!tiled() || f.layout().si == 0) {  // no decomposition; a horizontally uniform field
    whole = f.layout();
    f.download(h);
    return;
  }
  std::vector<char> b;
  download(f.base(), b, whole);
  h.assign(reinterpret_cast<const real*>(b.data()), reinterpret_cast<const real*>(b.data()) + b.size() / sizeof(real));
}
// the layout of a field with layout L in a run with one patch
Layout whole_layout(const Layout& L);
// The whole field on every worker (collective), as tile::download gives it to the writer: every point from its owner
// (the lateral halo from the edge patches).
void gather(const void* base, std::vector<char>& h, Layout& whole);
void gather_box(const void* base, const std::function<PatchBox(int)>& box, std::vector<char>& h, Layout& mine);
// A field of a decomposed domain as a whole on this worker's device (nests: the kernels that read one domain and
// write another, whose patches cover other areas). get() is collective; without a decomposition (and for a
// horizontally uniform field) on() stays false and the caller reads the field itself.
class WholeCopy {
 public:
  WholeCopy() = default;
  ~WholeCopy();
  WholeCopy(const WholeCopy&) = delete;
  WholeCopy& operator=(const WholeCopy&) = delete;
  void get(const void* base, const Layout& L);
  // only a box of the field: box(t) = the points (i0, i1, j0, j1 of the field's indices) worker t reads, a function
  // every worker evaluates alike; each worker receives its own box (layout(): that box under the field's indices)
  void get(const void* base, const Layout& L, const std::function<PatchBox(int)>& box);
  bool on() const { return on_; }
  const Layout& layout() const { return L_; }
  template <class T>
  T* data() const { return reinterpret_cast<T*>(buf_) + L_.offset; }  // element (0, 0, 0)

 private:
  char* buf_ = nullptr;
  std::size_t bytes_ = 0, esz_ = 0;
  Layout L_;
  bool on_ = false;
};
// The calling worker's field with allocation `base` as a number that names the same field on every worker, and the
// field of a number: its allocation, layout and element size.
struct FieldRef {
  void* base;
  Layout L;
  int esz;
  template <class T>
  T* data() const { return static_cast<T*>(base) + L.offset; }
};
int field_index(const void* base);
FieldRef field_at(int index);
// the halo exchange of one field by its allocation (as tile::exchange; also outside a time step's clipping)
void exchange_base(const void* base);
// The maximum over the workers; the smallest non-zero value over the workers (0: none)
double max(double v);
int first_code(int code);
// v[q] = the maximum over the workers of their v[q] (n values; every worker the same n)
void max_each(double* v, std::size_t n);
// the mass points the calling worker owns of decomposition dom, inside the domain: i0, i1, j0, j1
PatchBox own_box(const Grid& g);
// v (n bytes) of the worker that has it (have true on exactly one worker) on every worker
void from_owner(bool have, void* v, std::size_t n);
// Raw messages between workers for the collectives of single sites (each a fixed protocol both sides compute)
void send(int to, const void* p, std::size_t n);
void recv(int from, void* p, std::size_t n);

// The writer's host code that reads and sets whole fields (the WRF output, restart files): only the writer runs the
// scope's code; the other workers answer its requests (tile::download, tile::upload) until the scope ends.
//   tile::Service svc;
//   if (svc.served()) return;   // not the writer: the scope is over
// Inside the scope the writer makes no other collective call (no launches over the domain, no halo fills).
class Service {
 public:
  Service();
  ~Service();
  Service(const Service&) = delete;
  Service& operator=(const Service&) = delete;
  bool served() const { return served_; }

 private:
  bool served_ = false, opened_ = false;
};
// Inside a Service scope, on the writer: the field's points (no lateral halo) from h, a whole field in the layout
// `whole` (tile::download's), to the workers that hold them (their own points and their halo's points inside the
// domain). Without a decomposition: the points into the field.
void upload(const void* base, const char* h, const Layout& whole);
template <class F>
void upload(F& f, const std::vector<real>& h, const Layout& whole) {
  upload(f.base(), reinterpret_cast<const char*>(h.data()), whole);
}

// f's download as a whole field: without a decomposition its own buffer and layout; in a decomposed run the whole
// field from the workers (collective; on the writer, h is empty elsewhere). Returns h's layout.
template <class F, class V>
Layout download_any(const F& f, V& h) {
  Layout L = f.layout();
  if (tiled() && L.si != 0) {
    std::vector<char> b;
    download(f.base(), b, L);
    using T = typename V::value_type;
    h.assign(reinterpret_cast<const T*>(b.data()), reinterpret_cast<const T*>(b.data()) + b.size() / sizeof(T));
  } else {
    f.download(h);
  }
  return L;
}
// Not collective: the worker's own buffer of f laid out as a whole field (the points it does not hold are 0), and
// back: the points the worker holds from a whole field into f. Without a decomposition: the buffer itself.
template <class T>
void window_copy(const Layout& L, const Layout& W, std::vector<T>& local, std::vector<T>& whole, bool to_whole) {
  const int nk = L.sk ? (L.nz > 1 ? L.nz : 1) : 1;
  for (int k = 0; k < nk; ++k)
    for (int j = L.sj ? L.wj0 : 0; j < (L.sj ? L.wj1 : 1); ++j)
      for (int i = L.wi0; i < L.wi1; ++i) {
        if (to_whole) whole[W.index(i, j, k)] = local[L.index(i, j, k)];
        else local[L.index(i, j, k)] = whole[W.index(i, j, k)];
      }
}
template <class F, class V>
Layout local_whole(const F& f, V& h) {
  const Layout& L = f.layout();
  if (!tiled() || L.si == 0) {
    f.download(h);
    return L;
  }
  V local;
  f.download(local);
  const Layout W = whole_layout(L);
  h.assign(W.total, typename V::value_type(0));
  window_copy(L, W, local, h, true);
  return W;
}
// the worker's buffer of f (to copy to the device) from the whole field h of local_whole
template <class F, class V>
V local_part(const F& f, V& h) {
  const Layout& L = f.layout();
  if (!tiled() || L.si == 0) return h;
  V local(L.total, typename V::value_type(0));
  window_copy(L, whole_layout(L), local, h, false);
  return local;
}

// Tiles: the host code of this scope one tile at a time (the start, the input and the output are not written for
// several threads of one process: NetCDF, the tables' caches). A collective call inside releases it while the tile
// waits for the others. MPI ranks and runs without a decomposition: nothing.
class Serial {
 public:
  Serial();
  ~Serial();
  Serial(const Serial&) = delete;
  Serial& operator=(const Serial&) = delete;
  void unlock();
  void lock();  // again after unlock()

 private:
  std::unique_lock<std::mutex> lk_;
};

}  // namespace cppwrf::tile
