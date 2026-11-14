// SPDX-License-Identifier: GPL-3.0-or-later
#include "core/tiles.hpp"

#include <algorithm>
#include <atomic>
#include <climits>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <exception>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "core/tiles_internal.hpp"

#if defined(CPPWRF_USE_CUDA)
#include <cuda_runtime.h>
#endif
#if defined(_OPENMP)
#include <omp.h>
#endif

namespace cppwrf::tile {

using detail::Aborted;
using detail::Comm;

namespace {

// A barrier that another tile's error releases (Aborted); the check's meetings of the tiles and the reference run
class Barrier {
 public:
  void init(int n) { n_ = n; }
  void wait() {
    std::unique_lock<std::mutex> lk(m_);
    if (cancel_) throw Aborted();
    const long g = gen_;
    if (++count_ == n_) {
      count_ = 0;
      ++gen_;
      cv_.notify_all();
      return;
    }
    cv_.wait(lk, [&] { return gen_ != g || cancel_; });
    if (gen_ == g) throw Aborted();
  }
  void cancel() {
    std::lock_guard<std::mutex> lk(m_);
    cancel_ = true;
    cv_.notify_all();
  }

 private:
  std::mutex m_;
  std::condition_variable cv_;
  int n_ = 1, count_ = 0;
  long gen_ = 0;
  bool cancel_ = false;
};

struct Entry {
  char* base = nullptr;
  Layout L;
  int esz = 0;
  std::string name;
  bool live = false;
};
struct Registry {
  std::vector<Entry> e;
  std::unordered_map<const void*, int> idx;  // live allocations by base
  std::map<const char*, int> by_addr;        // live allocations by the address of element (0, 0, 0) (a view's data)
};

struct Iv {
  int lo, hi;
  bool empty() const { return hi <= lo; }
  int n() const { return hi - lo; }
};
Iv meet(Iv a, Iv b) { return Iv{std::max(a.lo, b.lo), std::min(a.hi, b.hi)}; }

struct Decomp {
  int nx = 0, ny = 0;
  bool per_x = false, per_y = false;
  std::vector<PatchBox> patch;
};

// A worker: the thread's (process's) patch number, its fields and its domains' decompositions (every worker creates
// the same domains and fields in the same order)
struct Worker {
  int id = 0, n = 1, px = 1, py = 1;
  Comm* comm = nullptr;
  Registry reg;
  std::vector<Decomp> doms;
  std::uint32_t seq = 0;  // collectives so far (checked in every message)
  bool service = false;   // the writer inside a tile::Service scope
  int dev = 0;
#if defined(CPPWRF_USE_CUDA)
  char* dbuf = nullptr;  // device staging buffer of the messages
  std::size_t dbuf_bytes = 0;
#endif
};
thread_local Worker* W = nullptr;

// tiles of one process: what they share
struct Set {
  Options o;
  int n = 1;
  Registry* ref = nullptr;  // the check's reference run
  Barrier all;
  std::mutex serial;
  std::mutex err_m;
  std::exception_ptr err;
  bool err_root = false;
  std::string chk_msg;  // the check's first difference
  std::mutex chk_m;
  std::function<void()> cancel_net;
};
Set* S = nullptr;
thread_local long n_checks = 0;
// the check: fields whose values around the patch differed at the last check
thread_local std::vector<char> stale_before;
thread_local Registry* ref_reg = nullptr;  // the reference run's thread: its registry
thread_local std::unique_lock<std::mutex>* held_serial = nullptr;  // the Serial lock this thread holds

bool worker() { return W != nullptr; }
Registry* my_reg() { return W ? &W->reg : ref_reg; }

void reg_add(const void* base, const Layout& L, int esz, const std::string& name) {
  Registry* reg = my_reg();
  Entry e;
  e.base = static_cast<char*>(const_cast<void*>(base));
  e.L = L;
  e.esz = esz;
  e.name = name;
  e.live = true;
  reg->idx[base] = int(reg->e.size());
  reg->by_addr[e.base + L.offset * std::ptrdiff_t(esz)] = int(reg->e.size());
  reg->e.push_back(std::move(e));
}
void reg_remove(const void* base) {
  Registry* reg = my_reg();
  auto it = reg->idx.find(base);
  if (it == reg->idx.end()) return;
  Entry& e = reg->e[std::size_t(it->second)];
  e.live = false;
  reg->by_addr.erase(e.base + e.L.offset * std::ptrdiff_t(e.esz));
  reg->idx.erase(it);
}
const Track kTrack{&reg_add, &reg_remove};

// a field of a decomposed domain
bool domain_field(const Entry& e) { return e.live && e.L.dom >= 0 && e.L.si == 1; }
bool two_d(const Layout& L) { return L.sj == 0; }
int levels(const Layout& L) { return L.sk ? std::max(L.nz, 1) : 1; }

// ---- geometry: all from the patch table, the same on every worker ----
// the indices worker t's allocation of a field with layout L (any worker's) holds
Iv win_x(const Decomp& D, int t, const Layout& L) {
  const PatchBox& p = D.patch[std::size_t(t)];
  return Iv{p.i0 - Grid::halo, (p.i1 == D.nx ? L.nx : p.i1) + Grid::halo};
}
Iv win_y(const Decomp& D, int t, const Layout& L) {
  if (two_d(L)) return Iv{0, 1};
  const PatchBox& p = D.patch[std::size_t(t)];
  return Iv{p.j0 - L.hy, (p.j1 == D.ny ? L.ny : p.j1) + L.hy};
}
// the points worker s owns of the field: its patch; at the domain's edges up to the end of the halo
Iv own_x(const Decomp& D, int s, const Layout& L) {
  const PatchBox& p = D.patch[std::size_t(s)];
  return Iv{p.i0 == 0 ? -Grid::halo : p.i0, p.i1 == D.nx ? L.nx + Grid::halo : p.i1};
}
Iv own_y(const Decomp& D, int s, const Layout& L) {
  if (two_d(L)) return Iv{0, 1};
  const PatchBox& p = D.patch[std::size_t(s)];
  return Iv{p.j0 == 0 ? -L.hy : p.j0, p.j1 == D.ny ? L.ny + L.hy : p.j1};
}

const Decomp& decomp(const Layout& L) {
  if (L.dom < 0 || L.dom >= int(W->doms.size())) throw std::logic_error("decomposition: a field of no decomposed domain");
  return W->doms[std::size_t(L.dom)];
}

// ---- boxes of a field <-> message bytes (k slowest, i fastest) ----
struct Piece {
  int entry;
  Iv x, y;
  std::size_t off = 0;  // in the message
};
std::size_t piece_bytes(const Entry& e, const Piece& p) {
  return std::size_t(p.x.n()) * std::size_t(p.y.n()) * std::size_t(levels(e.L)) * std::size_t(e.esz);
}
// the messages' offsets (8-byte aligned) and the total
std::size_t layout_pieces(std::vector<Piece>& ps) {
  std::size_t o = 0;
  for (Piece& p : ps) {
    p.off = o;
    o += (piece_bytes(W->reg.e[std::size_t(p.entry)], p) + 7) / 8 * 8;
  }
  return o;
}

#if defined(CPPWRF_USE_CUDA)
template <class U>
struct KBox {
  static constexpr bool kNoMulAdd = true;
  U* f;  // the field's element (0, 0, 0)
  U* b;  // the box in the staging buffer
  int sj, sk, x0, y0, nx, ny;
  bool put;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    U& a = f[(x0 + i) + (y0 + j) * sj + k * sk];
    U& c = b[(std::size_t(k) * std::size_t(ny) + std::size_t(j)) * std::size_t(nx) + std::size_t(i)];
    if (put) a = c;
    else c = a;
  }
};
char* staging(std::size_t bytes) {
  if (bytes > W->dbuf_bytes) {
    dev_free(W->dbuf);
    W->dbuf = nullptr;
    W->dbuf_bytes = 0;
    W->dbuf = static_cast<char*>(dev_alloc(bytes));
    W->dbuf_bytes = bytes;
  }
  return W->dbuf;
}
template <class U>
void launch_box(const Entry& e, const Piece& p, char* dbuf, bool put) {
  const Layout& L = e.L;
  const Whole whole;
  parallel_for(Range3D{0, p.x.n(), 0, p.y.n(), 0, levels(L)},
               KBox<U>{reinterpret_cast<U*>(e.base) + L.offset, reinterpret_cast<U*>(dbuf + p.off), int(L.sj), int(L.sk),
                       p.x.lo, two_d(L) ? 0 : p.y.lo, p.x.n(), p.y.n(), put});
}
void launch_box(const Entry& e, const Piece& p, char* dbuf, bool put) {
  if (e.esz == 8) launch_box<std::uint64_t>(e, p, dbuf, put);
  else launch_box<std::uint32_t>(e, p, dbuf, put);
}
#endif

// the pieces ps of this worker's fields into msg (get) or from msg (put); msg has layout_pieces(ps) bytes
void move_pieces(const std::vector<Piece>& ps, char* msg, std::size_t bytes, bool put) {
  if (ps.empty()) return;
#if defined(CPPWRF_USE_CUDA)
  char* d = staging(bytes);
  if (put) copy_to_dev(d, msg, bytes);
  for (const Piece& p : ps) launch_box(W->reg.e[std::size_t(p.entry)], p, d, put);
  if (!put) copy_to_host(msg, d, bytes);
#else
  (void)bytes;
  for (const Piece& p : ps) {
    const Entry& e = W->reg.e[std::size_t(p.entry)];
    const Layout& L = e.L;
    const std::size_t esz = std::size_t(e.esz), w = std::size_t(p.x.n()) * esz;
    char* m = msg + p.off;
    const int nk = levels(L);
    for (int k = 0; k < nk; ++k)
      for (int j = p.y.lo; j < p.y.hi; ++j) {
        char* f = e.base + (L.offset + p.x.lo + (two_d(L) ? 0 : j) * L.sj + k * L.sk) * std::ptrdiff_t(esz);
        if (put) std::memcpy(f, m, w);
        else std::memcpy(m, f, w);
        m += w;
      }
  }
#endif
}

// ---- messages: a header that names the collective, then the payload ----
struct Hdr {
  std::uint32_t seq, kind;
  std::uint64_t bytes;
};
enum Kind : std::uint32_t { kExchange = 1, kWrap, kMax, kSum, kCode, kOwner, kDownload, kRaw, kService, kUpload, kGather };
enum ServiceOp : std::int32_t { kOpDone = 0, kOpDownload, kOpUpload };

// CPPWRF_TILE_TRACE=1: every message on stderr (worker, collective, kind, peer, bytes)
bool trace() {
  static const bool on = [] {
    const char* e = std::getenv("CPPWRF_TILE_TRACE");
    return e && *e && *e != '0';
  }();
  return on;
}

void send_msg(int to, Kind kind, const char* p, std::size_t n) {
  if (trace()) std::fprintf(stderr, "tile %d: collective %u kind %u send to %d, %zu bytes\n", W->id, W->seq, unsigned(kind), to, n);
  std::vector<char> m(sizeof(Hdr) + n);
  const Hdr h{W->seq, kind, n};
  std::memcpy(m.data(), &h, sizeof h);
  if (n) std::memcpy(m.data() + sizeof h, p, n);
  W->comm->send(to, std::move(m));
}
// the next message from `from`: it must be of this collective (every worker makes the same collective calls in the
// same order: a different one is an error of the program, named here)
std::vector<char> recv_msg(int from, Kind kind, std::size_t n) {
  if (trace()) std::fprintf(stderr, "tile %d: collective %u kind %u recv from %d, %zu bytes\n", W->id, W->seq, unsigned(kind), from, n);
  std::vector<char> m = W->comm->recv(from);
  Hdr h{};
  if (m.size() >= sizeof h) std::memcpy(&h, m.data(), sizeof h);
  if (m.size() < sizeof h || h.seq != W->seq || h.kind != kind || h.bytes != n || m.size() != sizeof h + n) {
    char b[240];
    std::snprintf(b, sizeof b,
                  "decomposition: worker %d expected collective %u of kind %u with %zu bytes from worker %d and got "
                  "collective %u of kind %u with %llu bytes (the workers' calls differ)",
                  W->id, W->seq, unsigned(kind), n, from, h.seq, h.kind, static_cast<unsigned long long>(h.bytes));
    throw std::logic_error(b);
  }
  m.erase(m.begin(), m.begin() + std::ptrdiff_t(sizeof h));
  return m;
}
void end_collective() {
  W->comm->flush();
  ++W->seq;
}

// sends the pieces snd[r] to every worker r and receives rcv[s] from every worker s (empty lists: no message)
void transfer(std::vector<std::vector<Piece>>& snd, std::vector<std::vector<Piece>>& rcv, Kind kind) {
  dev_sync();  // this worker's launches are done
  std::vector<std::size_t> sb(std::size_t(W->n), 0), rb(std::size_t(W->n), 0);
  std::size_t stot = 0, rtot = 0;
  // one staging buffer for all the messages of a direction: the pieces' offsets continue from worker to worker
  std::vector<Piece> all;
  for (int r = 0; r < W->n; ++r) {
    sb[std::size_t(r)] = layout_pieces(snd[std::size_t(r)]);
    for (Piece& p : snd[std::size_t(r)]) {
      p.off += stot;
      all.push_back(p);
    }
    stot += sb[std::size_t(r)];
  }
  std::vector<char> buf(stot);
  move_pieces(all, buf.data(), stot, false);
  std::size_t o = 0;
  for (int r = 0; r < W->n; ++r) {
    if (sb[std::size_t(r)] > 0) send_msg(r, kind, buf.data() + o, sb[std::size_t(r)]);
    o += sb[std::size_t(r)];
  }
  all.clear();
  for (int s = 0; s < W->n; ++s) {
    rb[std::size_t(s)] = layout_pieces(rcv[std::size_t(s)]);
    for (Piece& p : rcv[std::size_t(s)]) {
      p.off += rtot;
      all.push_back(p);
    }
    rtot += rb[std::size_t(s)];
  }
  buf.assign(rtot, 0);
  o = 0;
  for (int s = 0; s < W->n; ++s) {
    if (rb[std::size_t(s)] > 0) {
      const std::vector<char> m = recv_msg(s, kind, rb[std::size_t(s)]);
      std::memcpy(buf.data() + o, m.data(), m.size());
    }
    o += rb[std::size_t(s)];
  }
  move_pieces(all, buf.data(), rtot, true);
  dev_sync();
  end_collective();
}

// the field with the allocation `p` (data false) or with the element (0, 0, 0) at `p` (a view's data pointer; with
// a patch it may lie outside the allocation, inside another field's)
int index_of(const void* p, bool data = false) {
  Registry* reg = &W->reg;
  if (!data) {
    auto it = reg->idx.find(p);
    if (it != reg->idx.end()) return it->second;
  } else {
    auto it = reg->by_addr.find(static_cast<const char*>(p));
    if (it != reg->by_addr.end()) return it->second;
  }
  throw std::logic_error("decomposition: exchange of a field the worker did not allocate");
}

std::vector<int> entries_of(const void* const* bases, int n, bool data = false) {
  std::vector<int> idx;
  idx.reserve(std::size_t(n));
  for (int k = 0; k < n; ++k) {
    const int i = index_of(bases[k], data);
    if (!domain_field(W->reg.e[std::size_t(i)]))
      throw std::logic_error("decomposition: exchange of " + W->reg.e[std::size_t(i)].name + ", not a field of a decomposed domain");
    idx.push_back(i);
  }
  return idx;
}

// the points the other workers own inside each worker's allocation
void exchange_entries(const std::vector<int>& idx) {
  const int t = W->id;
  std::vector<std::vector<Piece>> snd(std::size_t(W->n)), rcv(std::size_t(W->n));
  for (int i : idx) {
    const Entry& e = W->reg.e[std::size_t(i)];
    const Layout& L = e.L;
    const Decomp& D = decomp(L);
    for (int o = 0; o < W->n; ++o) {
      if (o == t) continue;
      const Iv sx = meet(win_x(D, o, L), own_x(D, t, L)), sy = meet(win_y(D, o, L), own_y(D, t, L));
      if (!sx.empty() && !sy.empty()) snd[std::size_t(o)].push_back(Piece{i, sx, sy});
      const Iv rx = meet(win_x(D, t, L), own_x(D, o, L)), ry = meet(win_y(D, t, L), own_y(D, o, L));
      if (!rx.empty() && !ry.empty()) rcv[std::size_t(o)].push_back(Piece{i, rx, ry});
    }
  }
  transfer(snd, rcv, kExchange);
}

void exchange_hook(const void* const* bases, int n, bool data) {
  if (!W) return;
  exchange_entries(entries_of(bases, n, data));
}

bool wraps_hook(int dom, int axis) {
  if (!W || dom < 0) return false;
  const Decomp& D = W->doms[std::size_t(dom)];
  return axis == 0 ? (D.per_x && W->px > 1) : (D.per_y && W->py > 1);
}

// The periodic halo of direction `axis` (several patches along it): the halo point i takes the value of the point
// i mod n (core/halo.cpp HaloLine: -3 .. -1 <- n - 3 .. n - 1; n .. <- 0 ..; faces also n + 3 <- 3) from the worker
// that owns it, over the receiver's tangential range: its rows (x) or its whole x range with the halo (y).
void wrap_hook(const void* const* bases, int n, int axis) {
  if (!W) return;
  const std::vector<int> idx = entries_of(bases, n);
  const int t = W->id;
  std::vector<std::vector<Piece>> snd(std::size_t(W->n)), rcv(std::size_t(W->n));
  for (int i : idx) {
    const Entry& e = W->reg.e[std::size_t(i)];
    const Layout& L = e.L;
    const Decomp& D = decomp(L);
    const int nn = axis == 0 ? D.nx : D.ny, ln = axis == 0 ? L.nx : L.ny;
    auto lo = [&](int w) { return axis == 0 ? D.patch[std::size_t(w)].i0 : D.patch[std::size_t(w)].j0; };
    auto hi = [&](int w) { return axis == 0 ? D.patch[std::size_t(w)].i1 : D.patch[std::size_t(w)].j1; };
    auto tang = [&](int w) { return axis == 0 ? w / W->px : w % W->px; };  // the row (x) or column (y) of patches
    // the tangential range of receiver w
    auto trange = [&](int w) { return axis == 0 ? meet(win_y(D, w, L), Iv{0, L.ny}) : win_x(D, w, L); };
    // from sender s to receiver r: the receiver's halo boxes (dst) and the sender's source boxes (dst shifted)
    auto pieces = [&](int s, int r, bool dst, std::vector<Piece>& out) {
      if (s == r || tang(s) != tang(r)) return;
      const Iv own{lo(s), hi(s)};
      const Iv tr = trange(r);
      for (int side = 0; side < 2; ++side) {
        if (side == 0 ? lo(r) != 0 : hi(r) != nn) continue;
        const Iv halo = side == 0 ? Iv{-Grid::halo, 0} : Iv{nn, ln + Grid::halo};
        const int shift = side == 0 ? nn : -nn;
        const Iv src = meet(Iv{halo.lo + shift, halo.hi + shift}, own);
        if (src.empty() || tr.empty()) continue;
        const Iv a = dst ? Iv{src.lo - shift, src.hi - shift} : src;
        out.push_back(axis == 0 ? Piece{i, a, tr} : Piece{i, tr, a});
      }
    };
    for (int o = 0; o < W->n; ++o) {
      pieces(t, o, false, snd[std::size_t(o)]);
      pieces(o, t, true, rcv[std::size_t(o)]);
    }
  }
  transfer(snd, rcv, kWrap);
}

// One line of a field from the other side of a periodic direction with several patches along it (hooks.line): the
// worker that owns the line `dst` receives, for the tangential positions it owns in [t0, t1), the line `src` from
// the worker that owns it (the same row or column of patches: it owns the same tangential positions).
bool line_hook(const void* data, int axis, int dst, int src, int t0, int t1) {
  if (!W) return false;
  const int i = index_of(data, true);
  const Entry& e = W->reg.e[std::size_t(i)];
  if (!domain_field(e)) throw std::logic_error("decomposition: a line of " + e.name + ", not a field of a decomposed domain");
  const Layout& L = e.L;
  const Decomp& D = decomp(L);
  if (!(axis == 0 ? D.per_x && W->px > 1 : D.per_y && W->py > 1)) return false;
  const int t = W->id;
  auto along = [&](int w) { return axis == 0 ? own_x(D, w, L) : own_y(D, w, L); };
  auto across = [&](int w) { return meet(axis == 0 ? own_y(D, w, L) : own_x(D, w, L), Iv{t0, t1}); };
  auto tang = [&](int w) { return axis == 0 ? w / W->px : w % W->px; };
  auto has = [](Iv v, int p) { return p >= v.lo && p < v.hi; };
  std::vector<std::vector<Piece>> snd(std::size_t(W->n)), rcv(std::size_t(W->n));
  auto piece = [&](int line, Iv tr) { return axis == 0 ? Piece{i, Iv{line, line + 1}, tr} : Piece{i, tr, Iv{line, line + 1}}; };
  for (int o = 0; o < W->n; ++o) {
    if (o == t || tang(o) != tang(t)) continue;
    // this worker sends the line src to o, which owns dst
    if (has(along(t), src) && has(along(o), dst) && !across(o).empty()) snd[std::size_t(o)].push_back(piece(src, across(o)));
    if (has(along(o), src) && has(along(t), dst) && !across(t).empty()) rcv[std::size_t(o)].push_back(piece(dst, across(t)));
  }
  transfer(snd, rcv, kWrap);
  return true;
}

// every worker's value on every worker
template <class T>
std::vector<T> share(T v, Kind kind) {
  for (int o = 0; o < W->n; ++o)
    if (o != W->id) send_msg(o, kind, reinterpret_cast<const char*>(&v), sizeof v);
  std::vector<T> all(std::size_t(W->n), v);
  for (int o = 0; o < W->n; ++o)
    if (o != W->id) {
      const std::vector<char> m = recv_msg(o, kind, sizeof v);
      std::memcpy(&all[std::size_t(o)], m.data(), sizeof v);
    }
  end_collective();
  return all;
}

double max_hook(double v) {
  if (!W) return v;
  const std::vector<double> all = share(v, kMax);
  double m = all[0];
  for (int s = 1; s < W->n; ++s) m = all[std::size_t(s)] > m ? all[std::size_t(s)] : m;
  return m;
}

// Reducer::sum of a run with one patch over the range r of the values v (k slowest, i fastest): the CPU build adds
// each level's values in (j, i) order and then the levels; the CUDA build each level's values in 256 strided
// partial sums combined by a tree (core/backend.hpp kred), then the levels
double ordered_sum(const Range3D& r, const std::vector<double>& v) {
  const std::size_t ni = std::size_t(r.i1 - r.i0), nj = std::size_t(r.j1 - r.j0), n = ni * nj;
  double acc = 0.0;
  for (int k = r.k0; k < r.k1; ++k) {
    const double* p = v.data() + std::size_t(k - r.k0) * n;
#if defined(CPPWRF_USE_CUDA)
    double sh[256];
    for (std::size_t t = 0; t < 256; ++t) {
      double a = 0.0;
      for (std::size_t q = t; q < n; q += 256) a = a + p[q];
      sh[t] = a;
    }
    for (int s = 128; s > 0; s >>= 1)
      for (int t = 0; t < s; ++t) sh[t] = sh[t] + sh[t + s];
    acc = acc + sh[0];
#else
    double a = 0.0;
    for (std::size_t q = 0; q < n; ++q) a = a + p[q];
    acc = acc + a;
#endif
  }
  return acc;
}

double sum_hook(const Range3D& r, const Range3D& own, const double* v) {
  const std::size_t nown = std::size_t(std::max(own.i1 - own.i0, 0)) * std::size_t(std::max(own.j1 - own.j0, 0)) *
                           std::size_t(std::max(own.k1 - own.k0, 0));
  if (!W) return 0.0;
  double result = 0.0;
  if (W->id != 0) {
    std::vector<char> m(sizeof(Range3D) + nown * sizeof(double));
    std::memcpy(m.data(), &own, sizeof own);
    if (nown) std::memcpy(m.data() + sizeof own, v, nown * sizeof(double));
    W->comm->send(0, [&] {
      std::vector<char> f(sizeof(Hdr) + m.size());
      const Hdr h{W->seq, kSum, m.size()};
      std::memcpy(f.data(), &h, sizeof h);
      std::memcpy(f.data() + sizeof h, m.data(), m.size());
      return f;
    }());
    const std::vector<char> back = recv_msg(0, kSum, sizeof(double));
    std::memcpy(&result, back.data(), sizeof result);
  } else {
    const std::size_t ni = std::size_t(std::max(r.i1 - r.i0, 0)), nj = std::size_t(std::max(r.j1 - r.j0, 0));
    std::vector<double> all(ni * nj * std::size_t(std::max(r.k1 - r.k0, 0)), 0.0);
    auto place = [&](const Range3D& o, const double* p) {
      const std::size_t oi = std::size_t(o.i1 - o.i0);
      for (int k = o.k0; k < o.k1; ++k)
        for (int j = o.j0; j < o.j1; ++j)
          std::memcpy(&all[(std::size_t(k - r.k0) * nj + std::size_t(j - r.j0)) * ni + std::size_t(o.i0 - r.i0)],
                      p + (std::size_t(k - o.k0) * std::size_t(o.j1 - o.j0) + std::size_t(j - o.j0)) * oi, oi * sizeof(double));
    };
    if (nown) place(own, v);
    for (int s = 1; s < W->n; ++s) {
      std::vector<char> m = W->comm->recv(s);
      Hdr h{};
      if (m.size() >= sizeof h) std::memcpy(&h, m.data(), sizeof h);
      if (m.size() < sizeof h + sizeof(Range3D) || h.seq != W->seq || h.kind != kSum)
        throw std::logic_error("decomposition: the workers' calls differ (a sum)");
      Range3D o;
      std::memcpy(&o, m.data() + sizeof h, sizeof o);
      if (o.i1 > o.i0 && o.j1 > o.j0 && o.k1 > o.k0)
        place(o, reinterpret_cast<const double*>(m.data() + sizeof h + sizeof o));
    }
    result = ordered_sum(r, all);
    for (int s = 1; s < W->n; ++s) send_msg(s, kSum, reinterpret_cast<const char*>(&result), sizeof result);
  }
  end_collective();
  return result;
}

std::string value_str(const char* p, int esz) {
  char b[64];
  if (esz == 8) {
    double v;
    std::memcpy(&v, p, 8);
    std::snprintf(b, sizeof b, "%.17g", v);
  } else {
    std::int32_t iv;
    float fv;
    std::memcpy(&iv, p, 4);
    std::memcpy(&fv, p, 4);
    std::snprintf(b, sizeof b, "%.9g (int %d)", double(fv), int(iv));
  }
  return b;
}

// the first point of the box x * y of the worker's field `mine` that differs from the reference run's field
// (same field, the whole domain); empty: none
std::string first_difference(const Entry& mine, const Entry& ref, Iv x, Iv y) {
  const Layout& L = mine.L;
  const Layout& R = ref.L;
  const int nk = levels(L);
  const std::size_t esz = std::size_t(mine.esz), bytes = std::size_t(x.n()) * esz;
  for (int k = 0; k < nk; ++k)
    for (int j = y.lo; j < y.hi; ++j) {
      const int jj = two_d(L) ? 0 : j;
      const char* a = mine.base + (L.offset + x.lo + jj * L.sj + k * L.sk) * std::ptrdiff_t(esz);
      const char* b = ref.base + (R.offset + x.lo + jj * R.sj + k * R.sk) * std::ptrdiff_t(esz);
      if (std::memcmp(a, b, bytes) == 0) continue;
      for (int ii = 0; ii < x.n(); ++ii) {
        const char* pa = a + std::size_t(ii) * esz;
        const char* pb = b + std::size_t(ii) * esz;
        if (std::memcmp(pa, pb, esz) == 0) continue;
        char h[200];
        std::snprintf(h, sizeof h, "%s (i, j, k) = (%d, %d, %d) from 0: ", mine.name.c_str(), x.lo + ii, j, k);
        return h + value_str(pa, mine.esz) + " (tile " + std::to_string(W->id) + "), " + value_str(pb, mine.esz) +
               " (reference)";
      }
    }
  return "";
}

// CPPWRF_TILE_CHECK: after every launch of a step, each tile's points against the reference run
void check_hook(const char* what) {
  if (!S || !S->o.check) return;
  ++n_checks;
  S->all.wait();
  if (W) {
    const Registry& r = *S->ref;
    const Registry* reg = &W->reg;
    std::string msg;
    if (r.e.size() != reg->e.size()) {
      msg = "the fields allocated differ from the reference run's";
    } else {
      for (std::size_t i = 0; i < reg->e.size() && msg.empty(); ++i) {
        const Entry& mine = reg->e[i];
        const Entry& ref = r.e[i];
        if (mine.live != ref.live || (mine.live && (mine.L.nx != ref.L.nx || mine.L.ny != ref.L.ny || mine.L.nz != ref.L.nz))) {
          msg = "the fields allocated differ from the reference run's (" + mine.name + ")";
          break;
        }
        if (!domain_field(mine)) continue;
        const Decomp& D = decomp(mine.L);
        msg = first_difference(mine, ref, meet(own_x(D, W->id, mine.L), win_x(D, W->id, mine.L)),
                               meet(own_y(D, W->id, mine.L), win_y(D, W->id, mine.L)));
      }
    }
    // the fields whose values around the patch are stale now (the inputs of the next launch)
    std::vector<char> stale(reg->e.size(), 0);
    for (std::size_t i = 0; i < reg->e.size() && i < r.e.size(); ++i) {
      const Entry& mine = reg->e[i];
      const Entry& ref = r.e[i];
      if (!domain_field(mine) || !ref.live) continue;
      const Decomp& D = decomp(mine.L);
      if (!first_difference(mine, ref, win_x(D, W->id, mine.L), win_y(D, W->id, mine.L)).empty()) stale[i] = 1;
    }
    if (!msg.empty()) {
      std::string names;
      int nn = 0;
      for (std::size_t i = 0; i < stale_before.size() && i < reg->e.size(); ++i)
        if (stale_before[i] && nn++ < 80) names += (names.empty() ? "" : ", ") + reg->e[i].name;
      msg += "; stale around tile " + std::to_string(W->id) + " before the launch: " + (names.empty() ? "none" : names);
    }
    stale_before.swap(stale);
    if (!msg.empty()) {
      std::lock_guard<std::mutex> lk(S->chk_m);
      if (S->chk_msg.empty())
        S->chk_msg = "tile check: " + msg + " after launch " + std::to_string(n_checks) + " of the run (" + what + ")";
    }
  }
  S->all.wait();
  if (!S->chk_msg.empty()) {
    if (W && W->id == 0) {
      std::fprintf(stderr, "cppwrf: %s\n", S->chk_msg.c_str());
      if (const char* f = std::getenv("CPPWRF_TILE_CHECK_LOG"); f && *f)
        if (std::FILE* fp = std::fopen(f, "a")) {
          std::fprintf(fp, "%s\n", S->chk_msg.c_str());
          std::fclose(fp);
        }
    }
    throw std::runtime_error(S->chk_msg);
  }
}

void fail(std::exception_ptr e, bool root) {
  {
    std::lock_guard<std::mutex> lk(S->err_m);
    if (!S->err || (root && !S->err_root)) {
      S->err = e;
      S->err_root = root;
    }
  }
  S->all.cancel();
  if (S->cancel_net) S->cancel_net();
}

std::vector<int> parse_ints(const std::string& s) {
  std::vector<int> v;
  std::size_t p = 0;
  while (p < s.size()) {
    std::size_t q = s.find(',', p);
    if (q == std::string::npos) q = s.size();
    v.push_back(std::stoi(s.substr(p, q - p)));
    p = q + 1;
  }
  return v;
}

// ---- tiles' transport: a mailbox per (receiver, sender) ----
class ThreadNet {
 public:
  explicit ThreadNet(int n) : n_(n), box_(std::size_t(n) * std::size_t(n)) {}
  void send(int from, int to, std::vector<char>&& m) {
    Mail& b = box_[std::size_t(to) * std::size_t(n_) + std::size_t(from)];
    {
      std::lock_guard<std::mutex> lk(b.m);
      b.q.push_back(std::move(m));
    }
    b.cv.notify_one();
  }
  std::vector<char> recv(int me, int from) {
    Mail& b = box_[std::size_t(me) * std::size_t(n_) + std::size_t(from)];
    std::vector<char> m;
    bool got = false;
    {
      std::unique_lock<std::mutex> lk(b.m);
      if (!b.q.empty()) {
        m = std::move(b.q.front());
        b.q.pop_front();
        got = true;
      }
    }
    if (got) return m;
    // wait for the sender without the tiles' Serial lock (the sender may need it to get here)
    std::unique_lock<std::mutex>* ser = held_serial && held_serial->owns_lock() ? held_serial : nullptr;
    if (ser) ser->unlock();
    {
      std::unique_lock<std::mutex> lk(b.m);
      b.cv.wait(lk, [&] { return !b.q.empty() || cancel_.load(); });
      if (!b.q.empty()) {
        m = std::move(b.q.front());
        b.q.pop_front();
        got = true;
      }
    }
    if (ser) ser->lock();
    if (!got) throw Aborted();
    return m;
  }
  void cancel() {
    cancel_.store(true);
    for (Mail& b : box_) {
      std::lock_guard<std::mutex> lk(b.m);
      b.cv.notify_all();
    }
  }

 private:
  struct Mail {
    std::mutex m;
    std::condition_variable cv;
    std::deque<std::vector<char>> q;
  };
  int n_;
  std::vector<Mail> box_;
  std::atomic<bool> cancel_{false};
};
class ThreadComm : public Comm {
 public:
  ThreadComm(ThreadNet& net, int me) : net_(net), me_(me) {}
  void send(int to, std::vector<char>&& m) override { net_.send(me_, to, std::move(m)); }
  std::vector<char> recv(int from) override { return net_.recv(me_, from); }

 private:
  ThreadNet& net_;
  int me_;
};

void set_ctx(const Grid& g) {
  ctx.on = ctx.patch = g.patched();
  if (!ctx.on) return;
  ctx.i0 = g.pi0 == 0 ? INT_MIN : g.pi0;
  ctx.i1 = g.pi1 == g.nx ? INT_MAX : g.pi1;
  ctx.j0 = (g.ny == 1 || g.pj0 == 0) ? INT_MIN : g.pj0;
  ctx.j1 = (g.ny == 1 || g.pj1 == g.ny) ? INT_MAX : g.pj1;
}

void set_hooks(bool multi) {
  hooks.multi = multi;
  hooks.max = &max_hook;
  hooks.exchange = &exchange_hook;
  hooks.wraps = &wraps_hook;
  hooks.wrap = &wrap_hook;
  hooks.sum = &sum_hook;
  hooks.line = &line_hook;
}

}  // namespace

namespace detail {
int run_worker(Comm& c, int id, const Options& o, const std::function<int(int)>& body) {
  Worker w;
  w.id = id;
  w.n = o.count();
  w.px = o.px;
  w.py = o.py;
  w.comm = &c;
  w.dev = current_device();
  W = &w;
  track = &kTrack;
  struct Reset {
    ~Reset() {
#if defined(CPPWRF_USE_CUDA)
      if (W && W->dbuf) dev_free(W->dbuf);
#endif
      track = nullptr;
      W = nullptr;
      ctx = Ctx{};
    }
  } reset;
  return body(id);
}
}  // namespace detail

Options options(const std::string& run_tiles) {
  Options o;
  std::string t = run_tiles;
  if (const char* e = std::getenv("CPPWRF_TILES"); e && *e) t = e;
  int px = 1, py = 1;
  char c = 0;
  if (t.empty() || std::sscanf(t.c_str(), "%d%c%d", &px, &c, &py) < 1) px = 1;
  if (c == 0) py = 1;
  else if (c != 'x' && c != 'X') throw std::runtime_error("run.tiles: P or PxQ (tiles along x and y), not '" + t + "'");
  if (px < 1 || py < 1) throw std::runtime_error("run.tiles: P or PxQ with P, Q >= 1, not '" + t + "'");
  o.px = px;
  o.py = py;
  if (const char* e = std::getenv("CPPWRF_TILE_DEVICES"); e && *e) o.devices = parse_ints(e);
  if (const char* e = std::getenv("CPPWRF_TILE_CHECK"); e && *e && *e != '0') o.check = true;
#if defined(CPPWRF_USE_CUDA)
  if (o.check) throw std::runtime_error("CPPWRF_TILE_CHECK: the CPU build only");
#endif
  return o;
}

int id() { return W ? W->id : (ref_reg ? -1 : 0); }
int count() { return W ? W->n : 1; }
bool tiled() { return W != nullptr; }
bool writer() { return !W || W->id == 0; }
int patches_x() { return W ? W->px : 1; }
int patches_y() { return W ? W->py : 1; }

void set_domain(Grid& g) {
  if (!W) return;
  const int px = W->px, py = W->py;
  const int nx = g.nx, ny = g.ny;
  if (ny == 1 && py > 1) throw std::runtime_error("run.tiles: a 2-D (x-z) domain has one patch along y");
  if ((px > 1 && nx / px < Grid::halo) || (py > 1 && ny / py < Grid::halo))
    throw std::runtime_error("run.tiles: the patches must be at least " + std::to_string(Grid::halo) +
                             " points wide (" + std::to_string(nx) + " x " + std::to_string(ny) + " points in " +
                             std::to_string(px) + " x " + std::to_string(py) + " patches)");
  Decomp D;
  D.nx = nx;
  D.ny = ny;
  D.per_x = g.bc_x == BC::Periodic;
  D.per_y = g.bc_y == BC::Periodic;
  D.patch.resize(std::size_t(W->n));
  for (int t = 0; t < W->n; ++t) {
    const int tx = t % px, ty = t / px;
    D.patch[std::size_t(t)] = PatchBox{int(long(tx) * nx / px), int(long(tx + 1) * nx / px), int(long(ty) * ny / py),
                                       int(long(ty + 1) * ny / py)};
  }
  g.dom = int(W->doms.size());
  const PatchBox& p = D.patch[std::size_t(W->id)];
  g.pi0 = p.i0;
  g.pi1 = p.i1;
  g.pj0 = p.j0;
  g.pj1 = p.j1;
  W->doms.push_back(std::move(D));
  set_ctx(g);
}

PatchBox patch(int dom, int t) {
  if (!W || dom < 0 || dom >= int(W->doms.size())) throw std::logic_error("decomposition: no such domain");
  return W->doms[std::size_t(dom)].patch[std::size_t(t)];
}

Use::Use(const Grid& g) : saved_(ctx) {
  if (W) set_ctx(g);
}
Use::~Use() { ctx = saved_; }

Step::Step() { in_step = true; }
Step::~Step() { in_step = false; }

void exchange_all() {
  if (!W) return;
  std::vector<int> idx;
  for (std::size_t i = 0; i < W->reg.e.size(); ++i)
    if (domain_field(W->reg.e[i])) idx.push_back(int(i));
  exchange_entries(idx);
}

Layout whole_layout(const Layout& L) {
  Layout w = L;
  if (L.si == 0) return w;
  w.dom = -1;
  w.wi0 = -Grid::halo;
  w.wi1 = L.nx + Grid::halo;
  w.wj0 = -L.hy;
  w.wj1 = L.ny + L.hy;
  const std::ptrdiff_t pitch = (Layout::lpad + L.nx + Grid::halo + 31) / 32 * 32;
  const std::ptrdiff_t plane = pitch * (L.ny + 2 * L.hy);
  w.sj = L.sj ? pitch : 0;
  w.sk = L.sk ? plane : 0;
  w.offset = Layout::lpad + L.hy * w.sj;
  w.total = std::size_t(plane) * std::size_t(L.sk ? L.nz : 1);
  return w;
}

namespace {
Piece own_piece(const Entry& e, int i, int s) {
  const Decomp& D = decomp(e.L);
  return Piece{i, meet(own_x(D, s, e.L), win_x(D, s, e.L)), meet(own_y(D, s, e.L), win_y(D, s, e.L))};
}
// the points of the domain (no lateral halo) worker s holds: its own and its halo's
Piece inner_piece(const Entry& e, int i, int s) {
  const Decomp& D = decomp(e.L);
  return Piece{i, meet(win_x(D, s, e.L), Iv{0, e.L.nx}), meet(win_y(D, s, e.L), Iv{0, two_d(e.L) ? 1 : e.L.ny})};
}
void request(ServiceOp op, int entry) {
  const std::int32_t r[2] = {op, entry};
  for (int o = 1; o < W->n; ++o) send_msg(o, kService, reinterpret_cast<const char*>(r), sizeof r);
}
// a worker's part of a download: its points of field i to the writer
void download_send(int i) {
  const Entry& e = W->reg.e[std::size_t(i)];
  dev_sync();
  std::vector<Piece> mine{own_piece(e, i, W->id)};
  const std::size_t nb = layout_pieces(mine);
  std::vector<char> buf(nb);
  move_pieces(mine, buf.data(), nb, false);
  send_msg(0, kDownload, buf.data(), nb);
}
void upload_recv(int i) {
  const Entry& e = W->reg.e[std::size_t(i)];
  std::vector<Piece> mine{inner_piece(e, i, W->id)};
  const std::size_t nb = layout_pieces(mine);
  std::vector<char> m = recv_msg(0, kUpload, nb);
  move_pieces(mine, m.data(), nb, true);
  dev_sync();
}
void check_domain_field(const Entry& e, const char* what) {
  if (!domain_field(e)) throw std::logic_error(std::string("decomposition: ") + what + " of " + e.name + ", not a field of a decomposed domain");
}
}  // namespace

namespace {
// a worker's piece p of a field (message bytes m) into the whole field h (layout whole)
void place_piece(const Layout& L, const Layout& whole, std::size_t esz, const Piece& p, const char* m, std::vector<char>& h) {
  const std::size_t w = std::size_t(p.x.n()) * esz;
  for (int k = 0; k < levels(L); ++k)
    for (int j = p.y.lo; j < p.y.hi; ++j) {
      std::memcpy(h.data() + (whole.offset + p.x.lo + (two_d(L) ? 0 : j) * whole.sj + k * whole.sk) * std::ptrdiff_t(esz), m, w);
      m += w;
    }
}
}  // namespace

void download(const void* base, std::vector<char>& h, Layout& whole) {
  if (!W) throw std::logic_error("decomposition: tile::download without a decomposition");
  const int i = index_of(base);
  const Entry& e = W->reg.e[std::size_t(i)];
  check_domain_field(e, "download");
  const Layout& L = e.L;
  whole = whole_layout(L);
  h.clear();
  if (W->service) request(kOpDownload, i);
  if (W->id != 0) {
    download_send(i);
  } else {
    dev_sync();
    std::vector<Piece> mine{own_piece(e, i, 0)};
    const std::size_t nb = layout_pieces(mine);
    std::vector<char> buf(nb);
    move_pieces(mine, buf.data(), nb, false);
    const std::size_t esz = std::size_t(e.esz);
    h.assign(whole.total * esz, 0);
    place_piece(L, whole, esz, mine[0], buf.data(), h);
    for (int s = 1; s < W->n; ++s) {
      std::vector<Piece> ps{own_piece(e, i, s)};
      const std::size_t n = layout_pieces(ps);
      const std::vector<char> m = recv_msg(s, kDownload, n);
      place_piece(L, whole, esz, ps[0], m.data(), h);
    }
  }
  end_collective();
}

void gather(const void* base, std::vector<char>& h, Layout& whole) {
  if (!W) throw std::logic_error("decomposition: tile::gather without a decomposition");
  if (W->service) throw std::logic_error("decomposition: tile::gather inside a tile::Service scope");
  const int i = index_of(base);
  const Entry& e = W->reg.e[std::size_t(i)];
  check_domain_field(e, "gather");
  const Layout& L = e.L;
  whole = whole_layout(L);
  dev_sync();
  std::vector<Piece> mine{own_piece(e, i, W->id)};
  const std::size_t nb = layout_pieces(mine);
  std::vector<char> buf(nb);
  move_pieces(mine, buf.data(), nb, false);
  for (int o = 0; o < W->n; ++o)
    if (o != W->id) send_msg(o, kGather, buf.data(), nb);
  const std::size_t esz = std::size_t(e.esz);
  h.assign(whole.total * esz, 0);
  place_piece(L, whole, esz, mine[0], buf.data(), h);
  for (int s = 0; s < W->n; ++s) {
    if (s == W->id) continue;
    std::vector<Piece> ps{own_piece(e, i, s)};
    const std::size_t n = layout_pieces(ps);
    const std::vector<char> m = recv_msg(s, kGather, n);
    place_piece(L, whole, esz, ps[0], m.data(), h);
  }
  end_collective();
}

namespace {
// the layout of the box b of a field with layout L (inside the whole field's allocation), under the field's indices
Layout box_layout(const Layout& L, const PatchBox& b) {
  Layout w = L;
  w.dom = -1;
  w.wi0 = std::max(b.i0, -Grid::halo);
  w.wi1 = std::max(w.wi0, std::min(b.i1, L.nx + Grid::halo));
  if (two_d(L)) {
    w.wj0 = 0;
    w.wj1 = 1;
  } else {
    w.wj0 = std::max(b.j0, -L.hy);
    w.wj1 = std::max(w.wj0, std::min(b.j1, L.ny + L.hy));
  }
  const std::ptrdiff_t pitch = (Layout::lpad + (w.wi1 - w.wi0) + 31) / 32 * 32;
  const std::ptrdiff_t plane = pitch * (w.wj1 - w.wj0);
  w.sj = L.sj ? pitch : 0;
  w.sk = L.sk ? plane : 0;
  w.offset = Layout::lpad - w.wi0 - (two_d(L) ? 0 : w.wj0) * w.sj;
  w.total = std::size_t(plane) * std::size_t(L.sk ? L.nz : 1);
  return w;
}
// the part of worker s's own piece of entry e inside the box layout B
Piece box_piece(const Entry& e, int i, int s, const Layout& B) {
  Piece p = own_piece(e, i, s);
  p.x = meet(p.x, Iv{B.wi0, B.wi1});
  if (!two_d(e.L)) p.y = meet(p.y, Iv{B.wj0, B.wj1});
  if (p.x.hi < p.x.lo) p.x.hi = p.x.lo;
  if (p.y.hi < p.y.lo) p.y.hi = p.y.lo;
  return p;
}
}  // namespace

void gather_box(const void* base, const std::function<PatchBox(int)>& box, std::vector<char>& h, Layout& mine_l) {
  if (!W) throw std::logic_error("decomposition: tile::gather_box without a decomposition");
  if (W->service) throw std::logic_error("decomposition: tile::gather_box inside a tile::Service scope");
  const int i = index_of(base);
  const Entry& e = W->reg.e[std::size_t(i)];
  check_domain_field(e, "gather_box");
  const Layout& L = e.L;
  dev_sync();
  for (int o = 0; o < W->n; ++o) {  // this worker's points in every other worker's box
    if (o == W->id) continue;
    std::vector<Piece> ps{box_piece(e, i, W->id, box_layout(L, box(o)))};
    const std::size_t nb = layout_pieces(ps);
    if (ps[0].x.n() == 0 || ps[0].y.n() == 0) continue;
    std::vector<char> buf(nb);
    move_pieces(ps, buf.data(), nb, false);
    send_msg(o, kGather, buf.data(), nb);
  }
  mine_l = box_layout(L, box(W->id));
  const std::size_t esz = std::size_t(e.esz);
  h.assign(mine_l.total * esz, 0);
  for (int s = 0; s < W->n; ++s) {
    std::vector<Piece> ps{box_piece(e, i, s, mine_l)};
    const std::size_t nb = layout_pieces(ps);
    if (ps[0].x.n() == 0 || ps[0].y.n() == 0) continue;
    if (s == W->id) {
      std::vector<char> buf(nb);
      move_pieces(ps, buf.data(), nb, false);
      place_piece(L, mine_l, esz, ps[0], buf.data(), h);
    } else {
      const std::vector<char> m = recv_msg(s, kGather, nb);
      place_piece(L, mine_l, esz, ps[0], m.data(), h);
    }
  }
  end_collective();
}

WholeCopy::~WholeCopy() { dev_free(buf_); }
void WholeCopy::get(const void* base, const Layout& L, const std::function<PatchBox(int)>& box) {
  if (!W || L.si == 0) return;
  std::vector<char> h;
  gather_box(base, box, h, L_);
  if (h.size() != bytes_) {
    dev_free(buf_);
    buf_ = static_cast<char*>(dev_alloc(std::max<std::size_t>(h.size(), 8)));
    bytes_ = h.size();
  }
  if (!h.empty()) copy_to_dev(buf_, h.data(), h.size());
  esz_ = L_.total ? h.size() / L_.total : 0;
  on_ = true;
}
void WholeCopy::get(const void* base, const Layout& L) {
  if (!W || L.si == 0) return;
  std::vector<char> h;
  gather(base, h, L_);
  if (h.size() != bytes_) {
    dev_free(buf_);
    buf_ = static_cast<char*>(dev_alloc(h.size()));
    bytes_ = h.size();
  }
  copy_to_dev(buf_, h.data(), h.size());
  esz_ = L_.total ? h.size() / L_.total : 0;
  on_ = true;
}

int field_index(const void* base) {
  if (!W) throw std::logic_error("decomposition: tile::field_index without a decomposition");
  return index_of(base);
}
FieldRef field_at(int index) {
  if (!W || index < 0 || index >= int(W->reg.e.size()) || !W->reg.e[std::size_t(index)].live)
    throw std::logic_error("decomposition: tile::field_at of no live field");
  const Entry& e = W->reg.e[std::size_t(index)];
  return FieldRef{e.base, e.L, e.esz};
}
void exchange_base(const void* base) {
  if (!W) return;
  exchange_entries({index_of(base)});
}

void upload(const void* base, const char* h, const Layout& whole) {
  if (!W) throw std::logic_error("decomposition: tile::upload without a decomposition");
  if (!W->service) throw std::logic_error("decomposition: tile::upload outside a tile::Service scope");
  const int i = index_of(base);
  const Entry& e = W->reg.e[std::size_t(i)];
  check_domain_field(e, "upload");
  const Layout& L = e.L;
  const std::size_t esz = std::size_t(e.esz);
  request(kOpUpload, i);
  for (int s = 0; s < W->n; ++s) {
    std::vector<Piece> ps{inner_piece(e, i, s)};
    const std::size_t n = layout_pieces(ps);
    std::vector<char> m(n);
    char* q = m.data();
    const Piece& p = ps[0];
    const std::size_t w = std::size_t(p.x.n()) * esz;
    for (int k = 0; k < levels(L); ++k)
      for (int j = p.y.lo; j < p.y.hi; ++j) {
        std::memcpy(q, h + (whole.offset + p.x.lo + (two_d(L) ? 0 : j) * whole.sj + k * whole.sk) * std::ptrdiff_t(esz), w);
        q += w;
      }
    if (s == 0) {
      move_pieces(ps, m.data(), n, true);
      dev_sync();
    } else {
      send_msg(s, kUpload, m.data(), n);
    }
  }
  end_collective();
}

Service::Service() {
  if (!W) return;
  if (W->id == 0) {
    if (W->service) return;  // a nested scope
    W->service = true;
    opened_ = true;
    return;
  }
  served_ = true;
  for (;;) {
    const std::vector<char> m = recv_msg(0, kService, 2 * sizeof(std::int32_t));
    std::int32_t r[2];
    std::memcpy(r, m.data(), sizeof r);
    if (r[0] == kOpDownload) download_send(r[1]);
    else if (r[0] == kOpUpload) upload_recv(r[1]);
    end_collective();
    if (r[0] == kOpDone) break;
  }
}
Service::~Service() {
  if (!opened_ || !W) return;
  W->service = false;
  try {
    request(kOpDone, 0);
    end_collective();
  } catch (...) {  // the run is stopping
  }
}

double max(double v) { return max_hook(v); }

int first_code(int code) {
  if (!W) return code;
  const std::vector<int> all = share(code, kCode);
  int c = 0;
  for (int v : all)
    if (v != 0 && (c == 0 || v < c)) c = v;
  return c;
}

void max_each(double* v, std::size_t n) {
  if (!W) return;
  for (int o = 0; o < W->n; ++o)
    if (o != W->id) send_msg(o, kMax, reinterpret_cast<const char*>(v), n * sizeof(double));
  for (int o = 0; o < W->n; ++o)
    if (o != W->id) {
      const std::vector<char> m = recv_msg(o, kMax, n * sizeof(double));
      const double* p = reinterpret_cast<const double*>(m.data());
      for (std::size_t q = 0; q < n; ++q) v[q] = p[q] > v[q] ? p[q] : v[q];
    }
  end_collective();
}

PatchBox input_window(int nx, int ny) {
  if (!W) return PatchBox{0, nx, 0, ny};
  const int px = W->px, py = W->py, tx = W->id % px, ty = W->id / px;  // set_domain's patches
  const int i0 = int(long(tx) * nx / px), i1 = int(long(tx + 1) * nx / px);
  const int j0 = int(long(ty) * ny / py), j1 = int(long(ty + 1) * ny / py);
  return PatchBox{std::max(0, i0 - Grid::halo), std::min(nx, i1 + Grid::halo), std::max(0, j0 - Grid::halo),
                  std::min(ny, j1 + Grid::halo)};
}

PatchBox own_box(const Grid& g) {
  if (!g.patched()) return PatchBox{0, g.nx, 0, g.ny};
  return PatchBox{g.pi0, g.pi1, g.pj0, g.pj1};
}

void from_owner(bool have, void* v, std::size_t n) {
  if (!W) return;
  std::vector<char> m(1 + n);
  m[0] = have ? 1 : 0;
  std::memcpy(m.data() + 1, v, n);
  for (int o = 0; o < W->n; ++o)
    if (o != W->id) send_msg(o, kOwner, m.data(), m.size());
  int owners = have ? 1 : 0;
  for (int o = 0; o < W->n; ++o)
    if (o != W->id) {
      const std::vector<char> r = recv_msg(o, kOwner, m.size());
      if (r[0]) {
        std::memcpy(v, r.data() + 1, n);
        ++owners;
      }
    }
  end_collective();
  if (owners != 1) throw std::logic_error("decomposition: a value with " + std::to_string(owners) + " owners");
}

void send(int to, const void* p, std::size_t n) {
  if (!W) throw std::logic_error("decomposition: tile::send without a decomposition");
  send_msg(to, kRaw, static_cast<const char*>(p), n);
}
void recv(int from, void* p, std::size_t n) {
  if (!W) throw std::logic_error("decomposition: tile::recv without a decomposition");
  const std::vector<char> m = recv_msg(from, kRaw, n);
  if (n) std::memcpy(p, m.data(), n);
}

Serial::Serial() {
  if (S && !held_serial) {  // a nested scope: the outer one holds the lock
    lk_ = std::unique_lock<std::mutex>(S->serial);
    held_serial = &lk_;
  }
}
Serial::~Serial() {
  if (held_serial == &lk_) held_serial = nullptr;
}
void Serial::unlock() {
  if (held_serial == &lk_) held_serial = nullptr;
  if (lk_.owns_lock()) lk_.unlock();
}
void Serial::lock() {
  if (!S || held_serial || !lk_.mutex()) return;
  lk_.lock();
  held_serial = &lk_;
}

int run(const Options& o, const std::function<int(int)>& body) {
  if (S || W) throw std::logic_error("tiles: a tiled run inside a decomposed run");
  Set set;
  set.o = o;
  set.n = o.count();
  Registry ref;
  set.ref = &ref;
  ThreadNet net(set.n);
  set.cancel_net = [&net] { net.cancel(); };
  std::vector<int> dev(std::size_t(set.n), current_device());
  if (!o.devices.empty())
    for (int t = 0; t < set.n; ++t) dev[std::size_t(t)] = o.devices[std::size_t(t) % o.devices.size()];
  set.all.init(set.n + (o.check ? 1 : 0));
  S = &set;
  set_hooks(true);
  if (o.check) hooks.check = &check_hook;
  int per = 1;
#if defined(_OPENMP)
  per = std::max(1, omp_get_max_threads() / (set.n + (o.check ? 1 : 0)));
#endif
  std::printf("cppwrf: %d x %d tiles%s, %d OpenMP threads each\n", o.px, o.py,
              o.check ? " and the check's reference run" : "", per);
  std::vector<int> rc(std::size_t(set.n), 0);
  auto thread_main = [&](int t) {
    try {
#if defined(CPPWRF_USE_CUDA)
      if (t >= 0) {
        const cudaError_t e = cudaSetDevice(dev[std::size_t(t)]);
        if (e != cudaSuccess) throw std::runtime_error(std::string("tiles: cudaSetDevice: ") + cudaGetErrorString(e));
      }
#endif
#if defined(_OPENMP)
      omp_set_num_threads(per);
#endif
      int r = 0;
      if (t >= 0) {
        ThreadComm comm(net, t);
        r = detail::run_worker(comm, t, o, body);
        rc[std::size_t(t)] = r;
      } else {  // the check's reference run: one patch, its own fields
        ref_reg = &ref;
        track = &kTrack;
        r = body(t);
        track = nullptr;
      }
      if (r != 0) fail(std::make_exception_ptr(std::runtime_error("a tile stopped")), false);
    } catch (const Aborted&) {
      fail(std::current_exception(), false);
    } catch (...) {
      fail(std::current_exception(), true);
    }
    track = nullptr;
    ref_reg = nullptr;
  };
  std::vector<std::thread> th;
  for (int t = 0; t < set.n; ++t) th.emplace_back(thread_main, t);
  if (o.check) th.emplace_back(thread_main, -1);
  for (std::thread& x : th) x.join();
  hooks = Hooks{};
  S = nullptr;
  if (set.err && set.err_root) std::rethrow_exception(set.err);
  for (int r : rc)
    if (r != 0) return r;
  if (set.err) std::rethrow_exception(set.err);
  return rc[0];
}

#if !defined(CPPWRF_MPI)
bool mpi_start(int*, char***) { return false; }
int mpi_run(const std::string&, const std::function<int(int)>&) {
  throw std::logic_error("this build has no MPI (CPPWRF_MPI)");
}
void mpi_finish() {}
void mpi_abort(int) {}
#else
namespace detail {
void install_hooks() { set_hooks(true); }
}  // namespace detail
#endif

}  // namespace cppwrf::tile
