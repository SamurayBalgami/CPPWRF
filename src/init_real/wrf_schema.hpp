#pragma once
// WRF I/O schema (tools/wrf_schema/gen_wrf_schema.py): which state variables WRF writes to a stream for
// a given namelist, with WRF's names, dimensions and attributes.
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "init_real/namelist.hpp"

namespace cppwrf::init_real {

struct SchemaCond {
  std::string opt, op, val;
};

struct SchemaVar {
  bool array4d = false;
  std::string varname, dataname, type, memorder, stagger, desc, units;
  std::string desc_full;  // untruncated (the boundary file uses it)
  int ndim = 0, ntl = 0;
  std::uint32_t s1 = 0, s2 = 0;
  std::vector<std::pair<std::string, std::string>> dims;  // (lo, hi) expressions, Fortran dimension order
  std::vector<std::string> dimnames;
  std::vector<SchemaCond> conds;  // in use if empty or any holds
};

struct SchemaMember {
  std::string array, member, dname, desc, units;
  std::vector<SchemaCond> conds;
  std::uint32_t s1 = 0, s2 = 0;
  bool boundary = false;
};

// Option values: namelist.input, then values set by real.exe itself, then the Registry defaults.
class WrfConfig {
 public:
  // dom: the domain id whose per-domain namelist entries apply (WRF's model_config_rec%x(id))
  WrfConfig(const Namelist& nl, const std::map<std::string, std::pair<std::string, std::string>>& defaults,
            int dom = 1);
  std::string get(const std::string& name, int dom = 0) const;  // "" if unknown; dom 0: this domain
  double num(const std::string& name, int dom = 0) const;
  int integer(const std::string& name, int dom = 0) const { return int(num(name, dom)); }
  bool logical(const std::string& name, int dom = 0) const;
  int domain() const { return dom_; }
  void set(const std::string& name, const std::string& value) { derived_[name] = value; }

 private:
  const Namelist& nl_;
  const std::map<std::string, std::pair<std::string, std::string>>& defaults_;  // name -> (group, default)
  std::map<std::string, std::string> derived_;
  int dom_ = 1;
};

// One variable as written to a file.
struct OutVar {
  std::string name, desc, units, stagger, memorder_att, coordinates;
  char type = 'f';  // 'f', 'd', 'i', 'l', 'c'
  std::vector<std::string> dimnames;  // NetCDF order without Time
  std::vector<std::size_t> dimlens;
  const SchemaVar* var = nullptr;     // state variable (nullptr for a 4-D array member)
  const SchemaMember* member = nullptr;
};

class WrfSchema {
 public:
  static WrfSchema load(const std::string& path);

  const std::map<std::string, std::pair<std::string, std::string>>& defaults() const { return defaults_; }
  // input stream (wrfinput) variables in WRF's write order. real_bytes: 4 or 8 (RWORDSIZE).
  std::vector<OutVar> input_vars(const WrfConfig& c, int ide, int jde, int kde, int real_bytes) const;
  // 4-D array members written to wrfbdy (boundary_table true), e.g. QVAPOR, QCLOUD, ...
  std::vector<const SchemaMember*> boundary_members(const std::string& array, const WrfConfig& c) const;
  bool in_use(const std::vector<SchemaCond>& conds, const WrfConfig& c) const;
  const SchemaVar* find(const std::string& varname) const;

 private:
  std::map<std::string, std::pair<std::string, std::string>> defaults_;
  std::vector<SchemaVar> vars_;
  std::vector<SchemaMember> members_;
};

}  // namespace cppwrf::init_real
