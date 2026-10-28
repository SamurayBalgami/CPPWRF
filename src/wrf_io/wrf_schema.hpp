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
  // dom: the domain id whose per-domain namelist entries apply (WRF's model_config_rec%x(id)). entries (the
  // Registry's entry count per option, WrfSchema::entries()): WRF's namelist read, i.e. a max_domains option
  // without a value for dom takes the Registry default and other options their first value; without it the
  // last value given stands for the missing domains.
  WrfConfig(const Namelist& nl, const std::map<std::string, std::pair<std::string, std::string>>& defaults,
            int dom = 1, const std::map<std::string, std::string>* entries = nullptr);
  std::string get(const std::string& name, int dom = 0) const;  // "" if unknown; dom 0: this domain
  double num(const std::string& name, int dom = 0) const;
  int integer(const std::string& name, int dom = 0) const { return int(num(name, dom)); }
  bool logical(const std::string& name, int dom = 0) const;
  int domain() const { return dom_; }
  void set(const std::string& name, const std::string& value) { derived_[name] = value; }

 private:
  const Namelist& nl_;
  const std::map<std::string, std::pair<std::string, std::string>>& defaults_;  // name -> (group, default)
  const std::map<std::string, std::string>* entries_ = nullptr;                  // name -> 1, max_domains, ...
  std::map<std::string, std::string> derived_;
  int dom_ = 1;
};

// Values WRF derives from the namelist before its I/O (setup_physics_suite,
// set_physics_rconfigs, check_nml_consistency) that the Registry package conditions use. real_exe: as
// real.exe sets them (use_wps_input = 1); otherwise as wrf.exe.
void derive_config(WrfConfig& c, int max_dom, bool real_exe);
// setup_physics_suite (module_check_a_mundo.F): a -1 physics option takes the scheme of physics_suite
// ('CONUS', 'tropical'; any case and quoting); other values, options and suites are returned unchanged.
int physics_suite_option(std::string physics_suite, const std::string& option, int value);

// One variable as written to a file.
struct OutVar {
  std::string name, desc, units, stagger, memorder_att, coordinates;
  char type = 'f';  // 'f', 'd', 'i', 'l', 'c'
  std::vector<std::string> dimnames;  // NetCDF order without Time
  std::vector<std::size_t> dimlens;
  const SchemaVar* var = nullptr;     // state variable (nullptr for a 4-D array member)
  const SchemaMember* member = nullptr;
};

// What iofields_filename did to the stream masks (WrfSchema::apply_iofields).
struct IofieldsResult {
  bool warning = false;  // modify_io_masks1's gavewarning: WRF stops unless ignore_iofields_warning
  std::string message;   // its last W A R N I N G
  std::vector<std::pair<int, std::string>> changed;  // (stream bit, DataName) of every mask the file changed
};

class WrfSchema {
 public:
  static WrfSchema load(const std::string& path);

  const std::map<std::string, std::pair<std::string, std::string>>& defaults() const { return defaults_; }
  const std::map<std::string, std::string>& entries() const { return entries_; }
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
                                  int bit, bool all_levels = false) const;
  // restart stream (wrfrst): word 2, bit 21; every time level under its own name (U_1, U_2)
  std::vector<OutVar> restart_vars(const WrfConfig& c, int ide, int jde, int kde, int real_bytes) const {
    return stream_vars(c, ide, jde, kde, real_bytes, 1, 21, true);
  }
  // 4-D array members written to wrfbdy (boundary_table true), e.g. QVAPOR, QCLOUD, ...
  std::vector<const SchemaMember*> boundary_members(const std::string& array, const WrfConfig& c) const;
  bool in_use(const std::vector<SchemaCond>& conds, const WrfConfig& c) const;
  // namelist iofields_filename (frame/module_domain.F modify_io_masks1): lines op:stream:id:fields with op
  // + or -, stream h (history, auxhist id) or i (input, auxinput id); fields by DataName among the variables
  // and 4-D array members in use for `c` (WRF's state variable list). Reading stops at the first blank line;
  // a line with a '#' anywhere is skipped.
  IofieldsResult apply_iofields(const std::string& path, const WrfConfig& c);
  const SchemaVar* find(const std::string& varname) const;

 private:
  std::map<std::string, std::pair<std::string, std::string>> defaults_;
  std::map<std::string, std::string> entries_;
  std::vector<SchemaVar> vars_;
  std::vector<SchemaMember> members_;
};

// The namelist's iofields_filename of domain c.domain() (a max_domains option: a domain without an entry
// has none) applied to `s` as WRF's modify_io_masks; a relative name is looked up in `dir` (WRF's run
// directory). Errors: a warning without ignore_iofields_warning (WRF stops), and a change of the masks of
// the input streams `input_bits` (25 + n for auxinput n; CPPWRF reads those files whole).
void apply_iofields_option(WrfSchema& s, const WrfConfig& c, const std::string& dir, const std::vector<int>& input_bits);

}  // namespace cppwrf::wrf_io
