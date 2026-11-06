// NetCDF-4 output with the chunks compressed on the host threads (src/wrf_io/nc_deflate.hpp) must have the bytes of
// HDF5's own deflate filter. Writes a WRF-like file through NcOut (float, double and int variables in the chunked,
// shuffled, deflate-2 layout of external/io_netcdf, odd dimensions for edge chunks, two records), reads it back and
// checks the values; tests/io_netcdf/nc_deflate_test.sh writes it once with the precomputed chunks and once with
// CPPWRF_NC_THREADS=0 (HDF5's own filter) and compares the files byte for byte.
//   test_nc_deflate FILE        exit 0: written and read back; 77: the filter could not be installed (no speed-up)
#undef NDEBUG  // the checks are asserts: keep them in Release builds
#include <netcdf.h>

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "wrf_io/nc_deflate.hpp"
#include "wrf_io/nc_out.hpp"

using cppwrf::wrf_io::NcOut;
namespace deflate = cppwrf::wrf_io::deflate;

int main(int argc, char** argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: test_nc_deflate FILE\n");
    return 2;
  }
  const std::string path = argv[1];
  const std::size_t nz = 7, ny = 9, nx = 11, ns = 4;
  std::uint32_t seed = 12345;
  auto noise = [&] {  // a small LCG: values that deflate compresses, but not to nothing
    seed = seed * 1664525u + 1013904223u;
    return double(seed >> 8) / double(1u << 24);
  };
  auto field = [&](std::size_t n, double base, double amp) {
    std::vector<double> v(n);
    for (std::size_t i = 0; i < n; ++i) v[i] = base + amp * std::sin(0.1 * double(i)) + 1e-3 * noise();
    return v;
  };
  struct Rec {
    std::vector<double> t, u, ph, soil, hgt;
    std::vector<int> lu;
  };
  std::vector<Rec> recs(2);
  for (Rec& r : recs) {
    r.t = field(nz * ny * nx, 290, 10);
    r.u = field(nz * ny * (nx + 1), 0, 20);
    r.ph = field((nz + 1) * ny * nx, 5000, 3000);
    r.soil = field(ns * ny * nx, 280, 5);
    r.hgt = field(ny * nx, 100, 50);
    r.lu.resize(ny * nx);
    for (int& x : r.lu) x = int(noise() * 20.0) + 1;
  }
  {
    NcOut nc(path, true);
    nc.dim("Time", 0);
    nc.dim("DateStrLen", 19);
    nc.dim("west_east", nx);
    nc.dim("south_north", ny);
    nc.dim("bottom_top", nz);
    nc.dim("west_east_stag", nx + 1);
    nc.dim("bottom_top_stag", nz + 1);
    nc.dim("soil_layers", ns);
    nc.var("Times", 'c', {"Time", "DateStrLen"}, {}, 0);
    nc.var("T", 'f', {"Time", "bottom_top", "south_north", "west_east"}, {{"MemoryOrder", "XYZ"}}, 104);
    nc.var("U", 'd', {"Time", "bottom_top", "south_north", "west_east_stag"}, {{"MemoryOrder", "XYZ"}}, 105);
    nc.var("PH", 'f', {"Time", "bottom_top_stag", "south_north", "west_east"}, {{"MemoryOrder", "XYZ"}}, 104);
    nc.var("TSLB", 'f', {"Time", "soil_layers", "south_north", "west_east"}, {{"MemoryOrder", "XYZ"}}, 104);
    nc.var("HGT", 'd', {"Time", "south_north", "west_east"}, {{"MemoryOrder", "XY "}}, 105);
    nc.var("LU_INDEX", 'i', {"Time", "south_north", "west_east"}, {{"MemoryOrder", "XY "}}, 106);
    nc.var("P_TOP", 'f', {"Time"}, {{"MemoryOrder", "0  "}}, 104);
    nc.enddef();
    for (std::size_t r = 0; r < recs.size(); ++r) {
      nc.put_text("Times", r == 0 ? "2026-10-12_00:00:00" : "2026-10-12_00:10:00", r);
      nc.stage("T", recs[r].t, r);
      nc.stage("U", recs[r].u, r);
      nc.stage("PH", recs[r].ph, r);
      nc.stage("TSLB", recs[r].soil, r);
      nc.stage("HGT", recs[r].hgt, r);
      nc.stage_int("LU_INDEX", recs[r].lu, r);
      nc.stage("P_TOP", {5000.0}, r);
      nc.flush();
      nc.sync();
      deflate::clear();
    }
  }
  const deflate::Stats& s = deflate::stats();
  std::printf("chunks precompressed %ld, taken by the filter %ld, compressed by the filter %ld\n", s.precomputed, s.hits,
              s.misses);
  // read back (the inflate side of the filter)
  int id;
  assert(nc_open(path.c_str(), NC_NOWRITE, &id) == NC_NOERR);
  auto check = [&](const char* name, std::size_t r, const std::vector<double>& want, bool f) {
    int v;
    assert(nc_inq_varid(id, name, &v) == NC_NOERR);
    std::vector<double> got(want.size());
    int nd;
    nc_inq_varndims(id, v, &nd);
    std::vector<int> dims(static_cast<std::size_t>(nd));
    nc_inq_vardimid(id, v, dims.data());
    std::vector<std::size_t> st(static_cast<std::size_t>(nd), 0), ct(static_cast<std::size_t>(nd), 1);
    st[0] = r;
    for (int d = 1; d < nd; ++d) nc_inq_dimlen(id, dims[std::size_t(d)], &ct[std::size_t(d)]);
    assert(nc_get_vara_double(id, v, st.data(), ct.data(), got.data()) == NC_NOERR);
    for (std::size_t i = 0; i < want.size(); ++i) assert(got[i] == (f ? double(float(want[i])) : want[i]));
  };
  for (std::size_t r = 0; r < recs.size(); ++r) {
    check("T", r, recs[r].t, true);
    check("U", r, recs[r].u, false);
    check("PH", r, recs[r].ph, true);
    check("TSLB", r, recs[r].soil, true);
    check("HGT", r, recs[r].hgt, false);
    check("LU_INDEX", r, std::vector<double>(recs[r].lu.begin(), recs[r].lu.end()), false);
  }
  nc_close(id);
  const char* e = std::getenv("CPPWRF_NC_THREADS");
  if (e && std::atoi(e) <= 1) {  // HDF5's own filter, or the filter compressing each chunk itself
    assert(s.hits == 0);
    return 0;
  }
  if (!deflate::active()) {
    std::printf("SKIP: the deflate filter is not installed (HDF5 without H5Z_register)\n");
    return 77;
  }
  // every chunk of the deflated variables (T, U, PH, TSLB: 4 per record; HGT, LU_INDEX: 4) came from the pool
  assert(s.precomputed == 2 * 6 * 4);
  assert(s.hits == s.precomputed && s.misses == 0);
  std::printf("PASS\n");
  return 0;
}
