#pragma once
// Parallel compression of NetCDF-4 output without changing a byte of the files.
//
// WRF's NetCDF-4 output (external/io_netcdf: shuffle + deflate level 2 per chunk) is compressed by HDF5's deflate
// filter, one chunk after the other in the writing thread. Here the chunks of a record are cut, shuffled and
// compressed beforehand on a pool of host threads, with the same zlib call as HDF5's filter (compress2 at the
// variable's level). HDF5's deflate filter is replaced by one that takes the precomputed result when the chunk it is
// given has the same bytes (compared in full) and otherwise compresses as HDF5 does (H5Zdeflate.c, H5Z__filter_deflate
// of HDF5 1.10 to 1.14), so the file is the same either way; a chunk HDF5 pads differently only loses the speed-up.
// Only zlib runs on the pool threads; every HDF5 / NetCDF call stays in the caller's thread.
// CPPWRF_NC_THREADS=N: pool size (default: the hardware threads, at most 64); 1 switches the precomputation off (the
// filter compresses each chunk itself); 0 leaves HDF5's own filter in place.
#include <cstddef>
#include <vector>

namespace cppwrf::wrf_io::deflate {

// one record of a chunked, deflated variable in memory type = file type, row-major over len (Time excluded)
struct Record {
  const unsigned char* data = nullptr;
  std::vector<std::size_t> len, chunk;  // dimension lengths and chunk lengths, slowest first
  std::size_t elsize = 0;
  std::vector<unsigned char> fill;  // the fill value of the type (elsize bytes): HDF5 pads edge chunks with it
  bool shuffle = true;
  int level = 2;
};

bool active();  // the filter is installed and the pool has more than one thread
void install();  // replaces HDF5's deflate filter (once; no-op without HDF5)
void submit(const Record& r);  // queues the compression of each chunk of r (r.data must live until wait())
void wait();  // until every queued chunk is compressed
void clear();  // drops the results HDF5 did not take

struct Stats {
  long precomputed = 0, hits = 0, misses = 0;  // chunks: compressed on the pool, taken by the filter, compressed by it
  double miss_s = 0, wait_s = 0;  // filter's own compression; waiting for the pool
  double download_s = 0;  // device-to-host copies of the output fields (cases/wrf_history.cpp)
};
Stats& stats();
bool timing();  // CPPWRF_IO_TIMING=1: per-record timings of the output files

}  // namespace cppwrf::wrf_io::deflate
