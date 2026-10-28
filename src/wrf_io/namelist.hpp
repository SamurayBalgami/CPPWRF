#pragma once
// Reader for WRF's namelist.input (Fortran namelist syntax: &group, key = v1, v2, ..., /).
#include <cstdlib>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace cppwrf::wrf_io {

class Namelist {
 public:
  static Namelist read(const std::string& path);
  static Namelist parse(const std::string& text);

  // Value for domain `dom` (1-based); a single value applies to all domains.
  std::optional<std::string> raw(const std::string& group, const std::string& key, int dom = 1) const;
  double get_real(const std::string& group, const std::string& key, double def, int dom = 1) const;
  std::vector<std::string> values(const std::string& group, const std::string& key) const;
  const std::map<std::string, std::map<std::string, std::vector<std::string>>>& groups() const { return g_; }

 private:
  std::map<std::string, std::map<std::string, std::vector<std::string>>> g_;  // lower-case names
};

}  // namespace cppwrf::wrf_io
