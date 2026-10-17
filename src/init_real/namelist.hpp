#pragma once
// Reader for WRF's namelist.input (Fortran namelist syntax: &group, key = v1, v2, ..., /).
#include <cstdlib>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace cppwrf::init_real {

class Namelist {
 public:
  static Namelist read(const std::string& path);
  static Namelist parse(const std::string& text);

  bool has(const std::string& group, const std::string& key) const;
  // Value for domain `dom` (1-based); a single value applies to all domains.
  std::optional<std::string> raw(const std::string& group, const std::string& key, int dom = 1) const;
  int get_int(const std::string& group, const std::string& key, int def, int dom = 1) const;
  double get_real(const std::string& group, const std::string& key, double def, int dom = 1) const;
  float get_real_f(const std::string& group, const std::string& key, float def, int dom = 1) const;
  bool get_bool(const std::string& group, const std::string& key, bool def, int dom = 1) const;
  std::string get_str(const std::string& group, const std::string& key, const std::string& def, int dom = 1) const;
  std::vector<std::string> values(const std::string& group, const std::string& key) const;

 private:
  std::map<std::string, std::map<std::string, std::vector<std::string>>> g_;  // lower-case names
};

}  // namespace cppwrf::init_real
