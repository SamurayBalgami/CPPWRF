#pragma once
// WRF I/O schema (tools/wrf_schema/gen_wrf_schema.py): which state variables WRF writes to a stream for
// a given namelist, with WRF's names, dimensions and attributes.
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "wrf_io/namelist.hpp"

namespace cppwrf::wrf_io {

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
  std::string nest_down, nest_up;  // Registry nest interpolation (d) / feedback (u) function, "" if none
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

// Values WRF derives from the namelist before its I/O (setup_physics_suite with physics_suite = 'none',
// set_physics_rconfigs, check_nml_consistency) that the Registry package conditions use. real_exe: as
// real.exe sets them (use_wps_input = 1); otherwise as wrf.exe.
void derive_config(WrfConfig& c, int max_dom, bool real_exe);

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
  std::vector<OutVar> input_vars(const WrfConfig& c, int ide, int jde, int kde, int real_bytes) const {
    return stream_vars(c, ide, jde, kde, real_bytes, 0, 25);
  }
  // history stream (wrfout): word 1, bit 0
  std::vector<OutVar> history_vars(const WrfConfig& c, int ide, int jde, int kde, int real_bytes) const {
    return stream_vars(c, ide, jde, kde, real_bytes, 0, 0);
  }
  // variables of the stream whose bit is `bit` of streams word `word` (0 or 1; 31 bits per word)
  std::vector<OutVar> stream_vars(const WrfConfig& c, int ide, int jde, int kde, int real_bytes, int word,
                                  int bit) const;
  // 4-D array members written to wrfbdy (boundary_table true), e.g. QVAPOR, QCLOUD, ...
  std::vector<const SchemaMember*> boundary_members(const std::string& array, const WrfConfig& c) const;
  bool in_use(const std::vector<SchemaCond>& conds, const WrfConfig& c) const;
  // namelist iofields_filename (frame/module_domain.F modify_io_masks1): lines op:stream:id:fields with op
  // + or -, stream h (history, auxhist id) or i (input, auxinput id); fields by DataName. Returns false
  // when the file cannot be read.
  bool apply_iofields(const std::string& path);
  const SchemaVar* find(const std::string& varname) const;

 private:
  std::map<std::string, std::pair<std::string, std::string>> defaults_;
  std::vector<SchemaVar> vars_;
  std::vector<SchemaMember> members_;
};

}  // namespace cppwrf::wrf_io
