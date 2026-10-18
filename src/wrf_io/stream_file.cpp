#include "wrf_io/stream_file.hpp"

#include <cstdio>
#include <map>
#include <stdexcept>

namespace cppwrf::wrf_io {

void define_stream_vars(NcOut& nc, std::vector<OutVar>& vars, int real_bytes) {
  nc.dim("Time", 0);
  nc.dim("DateStrLen", 19);
  std::vector<std::pair<std::string, std::size_t>> defined = {{"Time", 0}, {"DateStrLen", 19}};
  for (OutVar& v : vars)
    for (int k = int(v.dimnames.size()) - 1; k >= 0; --k) {  // Fortran order of each variable
      std::string& n = v.dimnames[std::size_t(k)];
      const std::size_t len = v.dimlens[std::size_t(k)];
      if (n.empty()) {
        for (auto& [dn, dl] : defined)
          if (dl == len && dn != "Time") { n = dn; break; }
        if (n.empty()) {
          char b[16];
          std::snprintf(b, sizeof b, "DIM%04d", int(defined.size()));
          n = b;
        }
      }
      bool found = false;
      for (auto& [dn, dl] : defined) found = found || dn == n;
      if (!found) {
        nc.dim(n, len);
        defined.push_back({n, len});
      }
    }
  nc.var("Times", 'c', {"Time", "DateStrLen"}, {}, 0);
  const int ftr = real_bytes == 8 ? 105 : 104;
  for (const OutVar& v : vars) {
    std::vector<std::string> dims = {"Time"};
    for (const std::string& n : v.dimnames) dims.push_back(n);
    std::map<std::string, std::string> at = {{"MemoryOrder", v.memorder_att}, {"description", v.desc},
                                             {"units", v.units}, {"stagger", v.stagger}};
    if (!v.coordinates.empty()) at["coordinates"] = v.coordinates;
    const char t = v.type == 'l' ? 'i' : v.type;
    const int ft = t == 'f' ? 104 : t == 'd' ? 105 : t == 'i' ? 106 : ftr;
    nc.var(v.name, t, dims, at, ft);
  }
}

StreamFile::StreamFile(const std::string& path, std::vector<OutVar> vars, const WrfConfig& c, const DomainAtts& d,
                       Stream s, bool netcdf4)
    : nc_(std::make_unique<NcOut>(path, netcdf4)), vars_(std::move(vars)) {
  define_stream_vars(*nc_, vars_, d.real_bytes);
  write_global_atts(*nc_, c, d, s);
  nc_->enddef();
}

void StreamFile::write(const std::string& date, const StreamValueFn& value) {
  nc_->put_text("Times", date, rec_);
  std::vector<double> data;
  for (const OutVar& v : vars_) {
    std::size_t n = 1;
    for (std::size_t l : v.dimlens) n *= l;
    data.assign(n, 0.0);
    if (value(v, data) && data.size() != n)
      throw std::runtime_error("wrf_io: " + v.name + " has " + std::to_string(data.size()) + " values, expected " +
                               std::to_string(n));
    if (v.type == 'i' || v.type == 'l') nc_->put_int(v.name, std::vector<int>(data.begin(), data.end()), rec_);
    else if (v.type != 'c') nc_->put(v.name, data, rec_);
  }
  nc_->sync();
  ++rec_;
}

}  // namespace cppwrf::wrf_io
