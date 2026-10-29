#pragma once
// Minimal NetCDF writer for WRF-format files (one Time record): dimensions, global attributes and
// variables with WRF's per-variable attributes.
#include <map>
#include <string>
#include <vector>


namespace cppwrf::wrf_io {

class NcOut {
 public:
  // netcdf4: WRF built with NetCDF-4 features (deflate level 2 and chunking); else 64-bit offset classic
  explicit NcOut(const std::string& path, bool netcdf4 = true);
  ~NcOut();
  NcOut(const NcOut&) = delete;
  NcOut& operator=(const NcOut&) = delete;

  void dim(const std::string& name, std::size_t len);  // len 0: unlimited
  void att(const std::string& name, const std::string& v);
  void att(const std::string& name, int v);
  void att_f(const std::string& name, float v);
  void att_d(const std::string& name, double v);
  // a global attribute of a file in data mode (redef, put, enddef: WRF's ext_ncd_put_dom_ti in a written file)
  void reatt_real(const std::string& name, double v, int real_bytes);
  // define a variable; dims slowest first, e.g. {"Time","bottom_top","south_north","west_east"}
  void var(const std::string& name, char type /* 'f','d','i','c' */, const std::vector<std::string>& dims,
           const std::map<std::string, std::string>& satts, int field_type);
  void enddef();
  void put(const std::string& name, const std::vector<double>& v, std::size_t rec = 0);  // converted to the type
  void put_int(const std::string& name, const std::vector<int>& v, std::size_t rec = 0);
  void sync();
  void put_text(const std::string& name, const std::string& s, std::size_t rec = 0);

 private:
  int nc_ = -1;
  bool nc4_ = true;
  std::map<std::string, std::size_t> dimlen_;
  std::map<std::string, int> dims_, vars_;
  std::map<std::string, char> types_;
};

}  // namespace cppwrf::wrf_io
