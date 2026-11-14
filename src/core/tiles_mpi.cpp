// SPDX-License-Identifier: GPL-3.0-or-later
// MPI ranks as the workers of a decomposed run (core/tiles.hpp, docs/MPI.md): the transport of core/tiles.cpp's
// collectives over MPI. Built only with CPPWRF_MPI.
#if defined(CPPWRF_MPI)
#include <mpi.h>
#include <sys/resource.h>

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/tiles_internal.hpp"

#if defined(CPPWRF_USE_CUDA)
#include <cuda_runtime.h>
#endif

namespace cppwrf::tile {

namespace detail {
void install_hooks();
}

namespace {

void check(int rc, const char* what) {
  if (rc != MPI_SUCCESS) throw std::runtime_error(std::string("MPI: ") + what + " failed");
}

// Messages in order per (sender, receiver) pair (one tag: MPI does not overtake). A send is MPI_Isend of a buffer
// this object keeps until flush(), so two ranks that send to each other before they receive cannot block.
class MpiComm : public detail::Comm {
 public:
  void send(int to, std::vector<char>&& m) override {
    if (m.size() > std::size_t(INT_MAX)) throw std::runtime_error("MPI: a message of more than 2 GiB");
    pending_.push_back(std::move(m));
    reqs_.push_back(MPI_REQUEST_NULL);
    check(MPI_Isend(pending_.back().data(), int(pending_.back().size()), MPI_BYTE, to, 0, MPI_COMM_WORLD, &reqs_.back()),
          "MPI_Isend");
  }
  std::vector<char> recv(int from) override {
    MPI_Status st;
    check(MPI_Probe(from, 0, MPI_COMM_WORLD, &st), "MPI_Probe");
    int n = 0;
    check(MPI_Get_count(&st, MPI_BYTE, &n), "MPI_Get_count");
    std::vector<char> m(std::size_t(n > 0 ? n : 0));
    check(MPI_Recv(m.data(), n, MPI_BYTE, from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE), "MPI_Recv");
    return m;
  }
  void flush() override {
    if (reqs_.empty()) return;
    check(MPI_Waitall(int(reqs_.size()), reqs_.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");
    reqs_.clear();
    pending_.clear();
  }

 private:
  std::deque<std::vector<char>> pending_;
  std::vector<MPI_Request> reqs_;
};

bool started = false;

}  // namespace

bool mpi_start(int* argc, char*** argv) {
  int provided = 0;
  check(MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &provided), "MPI_Init_thread");
  started = true;
  int n = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &n);
  // a launcher that does not match the MPI library starts every process as a world of its own: each would run the
  // whole case as one process and write the same files
  if (n == 1) {
    const char* e = std::getenv("CPPWRF_RANKS");
    int p = 1, q = 1;
    if (e && std::sscanf(e, "%dx%d", &p, &q) >= 1 && p * q > 1)
      throw std::runtime_error(std::string("CPPWRF_RANKS=") + e + " but this process is the only MPI rank (does mpirun "
                               "belong to the MPI library the build links?)");
  }
  return n > 1;
}

void mpi_finish() {
  if (started) MPI_Finalize();
  started = false;
}

void mpi_abort(int code) {
  if (started) MPI_Abort(MPI_COMM_WORLD, code);
}

int mpi_run(const std::string& run_ranks, const std::function<int(int)>& body) {
  int n = 1, rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &n);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int rc = 2;
  try {
    Options o;
    std::string t = run_ranks;
    if (const char* e = std::getenv("CPPWRF_RANKS"); e && *e) t = e;
    if (t.empty()) {  // the most square P x Q = n with P <= Q (WRF's MPASPECT)
      o.px = 1;
      for (int p = 1; p * p <= n; ++p)
        if (n % p == 0) o.px = p;
      o.py = n / o.px;
    } else {
      char c = 0;
      if (std::sscanf(t.c_str(), "%d%c%d", &o.px, &c, &o.py) != 3 || (c != 'x' && c != 'X') || o.px < 1 || o.py < 1)
        throw std::runtime_error("run.ranks: PxQ (patches along x and y), not '" + t + "'");
      if (o.px * o.py != n)
        throw std::runtime_error("run.ranks = " + t + " needs " + std::to_string(o.px * o.py) + " MPI ranks, this run has " +
                                 std::to_string(n));
    }
#if defined(CPPWRF_USE_CUDA)
    // the rank's device: CPPWRF_RANK_DEVICES=d0,d1,... (rank r takes entry r mod the count); default the current one
    if (const char* e = std::getenv("CPPWRF_RANK_DEVICES"); e && *e) {
      std::vector<int> dev;
      for (const char* p = e; *p;) {
        dev.push_back(std::atoi(p));
        while (*p && *p != ',') ++p;
        if (*p == ',') ++p;
      }
      if (!dev.empty()) {
        const cudaError_t ce = cudaSetDevice(dev[std::size_t(rank) % dev.size()]);
        if (ce != cudaSuccess) throw std::runtime_error(std::string("MPI: cudaSetDevice: ") + cudaGetErrorString(ce));
      }
    }
#endif
    if (rank == 0) std::printf("cppwrf: %d MPI ranks as %d x %d patches\n", n, o.px, o.py);
    // the log is the writer's (rank 0); errors of every rank go to stderr
    else if (!std::freopen("/dev/null", "w", stdout)) throw std::runtime_error("MPI: cannot silence a rank's stdout");
    detail::install_hooks();
    MpiComm comm;
    rc = detail::run_worker(comm, rank, o, body);
    comm.flush();
    // the ranks' peak resident memory: a rank holds its patch, not the domain; the writer also the fields it
    // gathers for the output
    struct rusage ru {};
    getrusage(RUSAGE_SELF, &ru);
    long mine = ru.ru_maxrss / 1024;
    std::vector<long> all(static_cast<std::size_t>(n), 0);
    MPI_Gather(&mine, 1, MPI_LONG, all.data(), 1, MPI_LONG, 0, MPI_COMM_WORLD);
    if (rank == 0) {
      long other = 0, sum = 0;
      for (int r = 0; r < n; ++r) {
        sum += all[static_cast<std::size_t>(r)];
        if (r > 0 && all[static_cast<std::size_t>(r)] > other) other = all[static_cast<std::size_t>(r)];
      }
      std::printf("cppwrf: peak host memory: %ld MiB the writer (rank 0), %ld MiB the largest other rank, %ld MiB in all\n",
                  all[0], other, sum);
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "cppwrf: error (rank %d): %s\n", rank, e.what());
    std::fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, 2);
  }
  return rc;
}

}  // namespace cppwrf::tile
#endif
