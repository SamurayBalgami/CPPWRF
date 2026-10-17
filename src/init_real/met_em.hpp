#pragma once
// metgrid output (met_em.d0N.<date>.nc) as read by real.exe's input_auxinput1: every variable of the
// first time, global attributes, and the dimension names. Values are kept in double (met_em is
// float, so this is exact for both real.exe precisions).
#include <map>
#include <string>
#include <vector>

namespace cppwrf::init_real {

struct MetVar {
  std::vector<std::string> dims;  // without Time, slowest first (as in the file)
  std::vector<std::size_t> shape;
  std::vector<double> v;
};

class MetEm {
 public:
  static MetEm read(const std::string& path);

  bool has(const std::string& name) const { return vars_.count(name) != 0; }
  const MetVar& var(const std::string& name) const;
  bool has_att(const std::string& name) const { return iatt_.count(name) || ratt_.count(name) || satt_.count(name); }
  int att_int(const std::string& name, int def = 0) const;
  double att_real(const std::string& name, double def = 0) const;
  std::string att_str(const std::string& name, const std::string& def = "") const;
  int flag(const std::string& name) const { return att_int("FLAG_" + name, 0); }
  const std::map<std::string, MetVar>& vars() const { return vars_; }
  const std::map<std::string, std::vector<int>>& int_atts() const { return iatt_; }
  const std::map<std::string, std::vector<double>>& real_atts() const { return ratt_; }
  const std::map<std::string, std::string>& str_atts() const { return satt_; }
  std::string time;  // Times(0), e.g. 2026-10-12_00:00:00
  std::string path;

 private:
  std::map<std::string, MetVar> vars_;
  std::map<std::string, std::vector<int>> iatt_;
  std::map<std::string, std::vector<double>> ratt_;
  std::map<std::string, std::string> satt_;
};

}  // namespace cppwrf::init_real
