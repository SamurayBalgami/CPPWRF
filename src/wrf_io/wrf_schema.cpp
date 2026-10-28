#include "wrf_io/wrf_schema.hpp"
#include "core/strings.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace cppwrf::wrf_io {

namespace {

std::vector<std::string> split(const std::string& s, char sep) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : s) {
    if (c == sep) {
      out.push_back(cur);
      cur.clear();
    } else {
      cur += c;
    }
  }
  out.push_back(cur);
  return out;
}

std::vector<SchemaCond> parse_conds(const std::string& s) {
  std::vector<SchemaCond> out;
  if (s.empty() || s == "-") return out;
  for (const std::string& c : split(s, '|')) {
    std::vector<std::string> p = split(c, '.');
    if (p.size() < 3) throw std::runtime_error("schema: bad condition " + c);
    std::string val = p[2];
    for (std::size_t k = 3; k < p.size(); ++k) val += "." + p[k];  // real values / .TRUE.
    out.push_back({p[0], p[1], val});
  }
  return out;
}

std::string rtrim(std::string s) {
  while (!s.empty() && (s.back() == ' ' || s.back() == '\0')) s.pop_back();
  return s;
}

bool is_true(const std::string& v) {
  std::string t = to_lower(v);
  if (!t.empty() && t[0] == '.') t = t.substr(1);
  return !t.empty() && t[0] == 't';
}

}  // namespace

WrfConfig::WrfConfig(const Namelist& nl, const std::map<std::string, std::pair<std::string, std::string>>& defaults,
                     int dom, const std::map<std::string, std::string>* entries)
    : nl_(nl), defaults_(defaults), entries_(entries), dom_(dom) {}

std::string WrfConfig::get(const std::string& name0, int dom) const {
  if (dom == 0) dom = dom_;
  const std::string name = to_lower(name0);
  auto dv = derived_.find(name);
  if (dv != derived_.end()) return dv->second;
  auto d = defaults_.find(name);
  if (d != defaults_.end() && d->second.first != "derived") {
    if (entries_) {  // WRF's namelist read
      auto e = entries_->find(name);
      const std::vector<std::string> v = nl_.values(d->second.first, name);
      const std::size_t at = e != entries_->end() && e->second == "max_domains" ? std::size_t(dom - 1) : 0;
      return at < v.size() ? v[at] : d->second.second;
    }
    auto r = nl_.raw(d->second.first, name, dom);
    if (r) return *r;
    return d->second.second;
  }
  return d == defaults_.end() ? std::string() : d->second.second;
}

double WrfConfig::num(const std::string& name, int dom) const {
  std::string v = get(name, dom);
  if (v.empty()) return 0;
  std::string t = to_lower(v);
  if (t == ".true." || t == "t" || t == ".t.") return 1;
  if (t == ".false." || t == "f" || t == ".f.") return 0;
  for (char& c : v)
    if (c == 'd' || c == 'D') c = 'e';
  return std::strtod(v.c_str(), nullptr);
}

bool WrfConfig::logical(const std::string& name, int dom) const { return is_true(get(name, dom)); }

int physics_suite_option(std::string suite, const std::string& option, int value) {
  if (value != -1) return value;
  suite.erase(std::remove_if(suite.begin(), suite.end(), [](char ch) { return ch == '\'' || ch == '"' || ch == ' '; }),
              suite.end());
  suite = to_lower(suite);
  static const std::map<std::string, std::map<std::string, int>> suites = {
      {"conus", {{"mp_physics", 8}, {"cu_physics", 6}, {"ra_lw_physics", 4}, {"ra_sw_physics", 4},
                 {"bl_pbl_physics", 2}, {"sf_sfclay_physics", 2}, {"sf_surface_physics", 2}}},
      {"tropical", {{"mp_physics", 6}, {"cu_physics", 16}, {"ra_lw_physics", 4}, {"ra_sw_physics", 4},
                    {"bl_pbl_physics", 1}, {"sf_sfclay_physics", 1}, {"sf_surface_physics", 2}}}};
  auto it = suites.find(suite);
  if (it == suites.end()) return value;
  auto o = it->second.find(option);
  return o == it->second.end() ? value : o->second;
}

void derive_config(WrfConfig& c, int max_dom, bool real_exe) {
  for (const char* p : {"mp_physics", "cu_physics", "ra_lw_physics", "ra_sw_physics", "sf_sfclay_physics",
                        "sf_surface_physics", "bl_pbl_physics"}) {
    const int v = c.integer(p);
    const int w = std::max(0, physics_suite_option(c.get("physics_suite"), p, v));  // no suite: -1 is no package
    if (w != v) c.set(p, std::to_string(w));
  }
  c.set("cu_used", c.integer("cu_physics") != 0 ? "1" : "0");
  c.set("shcu_used", c.integer("shcu_physics") != 0 ? "1" : "0");
  c.set("gwd_used", c.integer("gwd_opt") == 1 ? "1" : c.integer("gwd_opt") == 3 ? "3" : "0");
  c.set("aercu_used", c.integer("aercu_opt") > 0 ? "1" : "0");
  c.set("fft_used", "0");
  c.set("cam_used", "0");
  c.set("use_wps_input", real_exe ? "1" : "0");
  // check_nml_consistency (WRF 4.8): ghg_input is reset to 0 when a domain has no shortwave radiation
  if (c.integer("ghg_input") == 1)
    for (int id = 1; id <= max_dom; ++id)
      if (c.integer("ra_sw_physics", id) == 0) c.set("ghg_input", "0");
  // MYNN EDMF options only with the MYNN PBL (bl_pbl_physics = 5)
  if (c.integer("bl_pbl_physics") != 5) {
    c.set("bl_mynn_edmf", "0");
    c.set("bl_mynn_edmf_dd", "0");
    c.set("bl_mynn_output", "0");
  }
  // cu_diag = 1 exactly for the Grell schemes and KF-CuP (GD 93, GFL 3, KF-CuP 10, G3 5)
  {
    const int cu = c.integer("cu_physics");
    c.set("cu_diag", cu == 93 || cu == 3 || cu == 10 || cu == 5 ? "1" : "0");
  }
  // mfshconv only with the QNSE PBL (bl_pbl_physics = 4)
  if (c.integer("bl_pbl_physics") != 4) c.set("mfshconv", "0");
  // wrf.exe: seed_dim from random_seed(size=...) of gfortran's runtime
  if (!real_exe) c.set("seed_dim", "8");
}

WrfSchema WrfSchema::load(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("schema: cannot open " + path);
  WrfSchema s;
  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::vector<std::string> t = split(line, '\t');
    if (t[0] == "RCONFIG" && t.size() >= 5) {
      s.defaults_[to_lower(t[1])] = {t[3], t[4]};
      if (t.size() >= 6) s.entries_[to_lower(t[1])] = t[5];
    } else if ((t[0] == "VAR" || t[0] == "ARRAY4D") && t.size() >= 15) {
      SchemaVar v;
      v.array4d = t[0] == "ARRAY4D";
      v.varname = t[1]; v.dataname = rtrim(t[2]); v.type = t[3]; v.ndim = std::atoi(t[4].c_str());
      v.memorder = t[5]; v.stagger = t[6]; v.ntl = std::atoi(t[7].c_str());
      v.s1 = std::uint32_t(std::strtoul(t[8].c_str(), nullptr, 10));
      v.s2 = std::uint32_t(std::strtoul(t[9].c_str(), nullptr, 10));
      if (!t[10].empty())
        for (const std::string& d : split(t[10], ';')) {
          std::size_t c = d.find(':');
          v.dims.push_back({d.substr(0, c), c == std::string::npos ? "" : d.substr(c + 1)});
        }
      v.dimnames = split(t[11], ';');
      v.desc = rtrim(t[12].substr(0, 80)); v.units = rtrim(t[13].substr(0, 80));  // CHARACTER*80 in WRF
      v.desc_full = rtrim(t[12]);
      v.conds = parse_conds(t[14]);
      if (t.size() >= 17) { v.nest_down = t[15]; v.nest_up = rtrim(t[16]); }
      s.vars_.push_back(std::move(v));
    } else if (t[0] == "MEMBER" && t.size() >= 10) {
      SchemaMember m;
      m.array = t[1]; m.member = t[2]; m.conds = parse_conds(t[3]); m.dname = rtrim(t[4]);
      m.desc = rtrim(t[5].substr(0, 80)); m.units = rtrim(t[6].substr(0, 80));
      m.s1 = std::uint32_t(std::strtoul(t[7].c_str(), nullptr, 10));
      m.s2 = std::uint32_t(std::strtoul(t[8].c_str(), nullptr, 10));
      m.boundary = is_true(t[9]);
      s.members_.push_back(std::move(m));
    }
  }
  return s;
}

bool WrfSchema::in_use(const std::vector<SchemaCond>& conds, const WrfConfig& c) const {
  if (conds.empty()) return true;
  for (const SchemaCond& k : conds) {
    double a = c.num(k.opt), b;
    std::string v = to_lower(k.val);
    if (v == ".true.") b = 1;
    else if (v == ".false.") b = 0;
    else b = std::strtod(k.val.c_str(), nullptr);
    bool r = k.op == "EQ" ? a == b : k.op == "NE" ? a != b : k.op == "GT" ? a > b : k.op == "GE" ? a >= b
             : k.op == "LT" ? a < b : a <= b;
    if (r) return true;
  }
  return false;
}

namespace {

// dimension bound expressions of the allocs tables
long eval_bound(const std::string& e0, const WrfConfig& c, int ide, int jde, int kde) {
  std::string e;
  for (char ch : e0)
    if (ch != ' ') e += ch;
  if (e.empty()) return 1;
  if (e.front() == '(' && e.back() == ')') e = e.substr(1, e.size() - 2);
  auto atom = [&](const std::string& a) -> long {
    if (a == "ids" || a == "jds" || a == "kds") return 1;
    if (a == "ide") return ide;
    if (a == "jde") return jde;
    if (a == "kde") return kde;
    if (a.rfind("config_flags%", 0) == 0) return long(c.num(a.substr(13)));
    if (a.rfind("model_config_rec%", 0) == 0) return long(c.num(a.substr(17, a.find('(') - 17)));
    if (a == "sr_x" || a == "sr_y") return long(c.num(a));
    return std::strtol(a.c_str(), nullptr, 10);
  };
  if (e.rfind("MAX(", 0) == 0) {  // MAX(1,ide * sr_x)
    std::string in = e.substr(4, e.size() - 5);
    std::size_t comma = in.find(',');
    long a = eval_bound(in.substr(0, comma), c, ide, jde, kde);
    long b = eval_bound(in.substr(comma + 1), c, ide, jde, kde);
    return std::max(a, b);
  }
  std::size_t op = e.find_first_of("-*+", 1);
  if (op != std::string::npos) {
    long a = atom(e.substr(0, op)), b = eval_bound(e.substr(op + 1), c, ide, jde, kde);
    return e[op] == '-' ? a - b : e[op] == '+' ? a + b : a * b;
  }
  return atom(e);
}

// false: a memory order WRF's netCDF output does not take, or a zero-length dimension (not written)
bool fill_out(OutVar& o, const SchemaVar& v, const WrfConfig& c, int ide, int jde, int kde, int real_bytes) {
  o.type = v.type == "r" ? (real_bytes == 8 ? 'd' : 'f') : v.type == "d" ? 'd' : v.type == "i" ? 'i'
           : v.type == "l" ? 'l' : 'c';
  std::vector<std::size_t> lens;
  for (int k = 0; k < v.ndim && k < int(v.dims.size()); ++k) {
    long lo = eval_bound(v.dims[std::size_t(k)].first, c, ide, jde, kde);
    long hi = eval_bound(v.dims[std::size_t(k)].second, c, ide, jde, kde);
    lens.push_back(std::size_t(std::max(0L, hi - lo + 1)));
  }
  // WRF writes X fastest, then Y, then Z: NetCDF order (Z, Y, X)
  const std::string& mo = v.memorder;
  auto pos = [&](char L) { return int(mo.find(L)); };
  std::vector<int> order;
  if (v.ndim == 0) {
    o.memorder_att = "0  ";
  } else if (mo == "XZY" || mo == "XYZ") {
    order = {pos('Z'), pos('Y'), pos('X')};
    o.memorder_att = "XYZ";
  } else if (mo == "XY") {
    order = {pos('Y'), pos('X')};
    o.memorder_att = "XY ";
  } else if (mo == "Z") {
    order = {0};
    o.memorder_att = "Z  ";
  } else if (mo.find('X') == std::string::npos && mo.find('Y') == std::string::npos && !mo.empty() && mo != "ZZ") {
    for (int k = int(mo.size()) - 1; k >= 0; --k) order.push_back(k);  // C, Z, CC, CZ, ZZ, CCZ: slowest first
    o.memorder_att = mo + std::string(3 - std::min<std::size_t>(3, mo.size()), ' ');
  } else {
    return false;
  }
  for (std::size_t l : lens)
    if (l == 0) return false;
  for (int p : order) {
    o.dimnames.push_back(v.dimnames[std::size_t(p)]);
    o.dimlens.push_back(lens[std::size_t(p)]);
  }
  o.stagger = v.stagger;
  if (mo.find('X') != std::string::npos && mo.find('Y') != std::string::npos) {
    const std::string& n = o.name;
    if (n == "XLAT" || n == "XLONG") o.coordinates = "XLONG XLAT";
    else if (n == "XLAT_U" || n == "XLONG_U") o.coordinates = "XLONG_U XLAT_U";
    else if (n == "XLAT_V" || n == "XLONG_V") o.coordinates = "XLONG_V XLAT_V";
    else if (v.stagger == "X") o.coordinates = "XLONG_U XLAT_U XTIME";
    else if (v.stagger == "Y") o.coordinates = "XLONG_V XLAT_V XTIME";
    else o.coordinates = "XLONG XLAT XTIME";
  }
  return true;
}

}  // namespace

std::vector<OutVar> WrfSchema::stream_vars(const WrfConfig& c, int ide, int jde, int kde, int real_bytes, int word,
                                           int bit, bool all_levels) const {
  const std::uint32_t mask = 1u << bit;
  auto on = [&](std::uint32_t s1, std::uint32_t s2) { return ((word == 0 ? s1 : s2) & mask) != 0; };
  std::vector<OutVar> out;
  for (const SchemaVar& v : vars_) {
    if (!on(v.s1, v.s2) || !in_use(v.conds, c)) continue;
    if (!all_levels && v.ntl > 0 && v.ntl / 100 != v.ntl % 100) continue;
    if (!v.array4d) {
      OutVar o;
      o.name = v.dataname;
      if (v.ntl > 0 && !all_levels) o.name = o.name.substr(0, o.name.size() - 2);
      o.desc = v.desc; o.units = v.units; o.var = &v;
      if (fill_out(o, v, c, ide, jde, kde, real_bytes)) out.push_back(std::move(o));
    } else {
      std::set<std::string> seen;
      for (const SchemaMember& m : members_) {
        if (m.array != v.varname || seen.count(m.member) || !in_use(m.conds, c)) continue;
        seen.insert(m.member);
        if (!on(m.s1, m.s2)) continue;
        OutVar o;
        o.name = m.dname; o.desc = m.desc; o.units = m.units; o.var = &v; o.member = &m;
        if (fill_out(o, v, c, ide, jde, kde, real_bytes)) out.push_back(std::move(o));
      }
    }
  }
  return out;
}

namespace {

std::string rtrim_blanks(std::string s) {  // Fortran TRIM: trailing blanks only
  while (!s.empty() && s.back() == ' ') s.pop_back();
  return s;
}

// frame/module_domain.F get_fieldstr: field ix (1-based) of `in` between separators c. A separator search
// starts two characters after the previous separator (an empty field swallows the next separator). false
// (and "") when there is no such field or it is longer than noutstr + 1 characters.
bool get_fieldstr(int ix, char c, const std::string& in, long noutstr, std::string& out) {
  out.clear();
  const long len = long(rtrim_blanks(in).size());
  std::vector<long> dexes = {1};
  long prev = 1;
  while (dexes.size() < 1000) {  // MAX_DEXES
    const std::size_t p = prev <= len ? in.find(c, std::size_t(prev - 1)) : std::string::npos;
    if (p != std::string::npos && long(p) < len) {
      dexes.push_back(long(p) + 2);  // idex + prev
      prev = dexes.back() + 1;
    } else {
      dexes.push_back(len + 2);
    }
  }
  if (ix < 1 || ix >= int(dexes.size())) return false;
  const long a = dexes[std::size_t(ix - 1)], b = dexes[std::size_t(ix)];
  if ((b - 2) - a > noutstr || a == b) return false;
  if (b - 1 > a) out = in.substr(std::size_t(a - 1), std::size_t(b - 1 - a));
  return true;
}

std::string lower_az(std::string s) {  // change_to_lower_case: A-Z only
  for (char& ch : s)
    if (ch >= 'A' && ch <= 'Z') ch = char(ch - 'A' + 'a');
  return s;
}

}  // namespace

IofieldsResult WrfSchema::apply_iofields(const std::string& path, const WrfConfig& c) {
  IofieldsResult r;
  std::ifstream f(path);
  if (!f) {
    r.warning = true;
    r.message = "W A R N I N G : Problem opening " + path;
    return r;
  }
  auto warn = [&](const std::string& m) { r.warning = true; r.message = "W A R N I N G : " + m; };
  // set_mask / reset_mask (frame/pack_utils.c): stream e is bit e % 31 of word e / 31 (two words)
  auto change = [&](std::uint32_t& s1, std::uint32_t& s2, char op, int e, const std::string& dname) {
    if (e < 0 || e / 31 > 1) return;
    std::uint32_t& w = e / 31 == 0 ? s1 : s2;
    const std::uint32_t m = 1u << (e % 31);
    const std::uint32_t old = w;
    if (op == '+') w |= m;
    else if (op == '-') w &= ~m;
    if (w != old) r.changed.push_back({e, dname});
  };
  std::string inln;
  int lineno = 0, strmtyp_int = -1;  // strmtyp_int keeps the previous line's value for an unknown stream type
  while (std::getline(f, inln)) {
    if (!inln.empty() && inln.back() == '\r') inln.pop_back();  // a CRLF line ends the record
    if (inln.size() > 8000) inln.resize(8000);                     // CHARACTER*8000 inln
    if (rtrim_blanks(inln).empty()) break;  // keepgoing = LEN(TRIM(inln)) > 0: a blank line ends the file
    ++lineno;
    const std::string at = " Line " + std::to_string(lineno) + " of " + path;
    if (rtrim_blanks(inln).size() >= 8000) warn("Line too long." + at);
    if (inln.find('#') != std::string::npos) continue;  // a comment: '#' anywhere on the line
    bool noerr = true;
    std::string t1;
    noerr = get_fieldstr(1, ':', inln, 1, t1) && noerr;
    const char op = t1.empty() ? ' ' : t1[0];  // CHARACTER*1 op
    if (op != '+' && op != '-') warn(std::string("unknown operation ") + op + " (should be + or -)." + at);
    noerr = get_fieldstr(2, ':', inln, 1, t1) && noerr;
    const char strmtyp = t1.empty() ? ' ' : lower_az(t1.substr(0, 1))[0];
    if (strmtyp == 'h') strmtyp_int = 1;        // first_history
    else if (strmtyp == 'i') strmtyp_int = 26;  // first_input
    else warn(std::string("unknown stream type ") + strmtyp + "." + at);
    noerr = get_fieldstr(3, ':', inln, 3, t1) && noerr;
    std::string strmid;  // READ(strmid,'(I3)'): CHARACTER*3, blanks ignored
    for (char ch : t1.substr(0, 3))
      if (ch != ' ') strmid += ch;
    int istrm = 0;
    if (!strmid.empty()) {
      std::size_t k = strmid[0] == '+' || strmid[0] == '-' ? 1 : 0;
      if (k == strmid.size() || strmid.find_first_not_of("0123456789", k) != std::string::npos)
        throw std::runtime_error("iofields_filename " + path + ": bad stream id '" + t1 + "' at line " +
                                 std::to_string(lineno) + " (a Fortran read error in WRF)");
      istrm = std::atoi(strmid.c_str());
    }
    if (istrm < 0 || istrm > 25) warn("invalid stream id " + std::to_string(istrm) + " (should be 0 <= id <= 25)." + at);
    std::string fieldlst;
    noerr = get_fieldstr(4, ':', inln, 8000, fieldlst) && noerr;
    if (!noerr) {
      warn("Problem reading " + path + " at line " + std::to_string(lineno));
      continue;
    }
    const int e = strmtyp_int < 0 ? -1 : strmtyp_int + istrm - 1;
    int fieldno = 1;
    noerr = get_fieldstr(fieldno, ',', fieldlst, 8000, t1);
    while (noerr) {
      const std::string lookee = rtrim_blanks(lower_az(t1.substr(0, 256)));  // CHARACTER*256 lookee
      bool found = false;
      for (SchemaVar& v : vars_) {
        if (v.array4d) {  // a scalar array: its members in use
          for (SchemaMember& m : members_)
            if (m.array == v.varname && in_use(m.conds, c) && rtrim_blanks(lower_az(m.dname)) == lookee) {
              found = true;
              change(m.s1, m.s2, op, e, m.dname);
            }
          continue;
        }
        if (!in_use(v.conds, c)) continue;  // not allocated: not on the state variable list
        std::string dn = v.dataname;
        if (v.ntl > 0 && dn.size() >= 2) dn.resize(dn.size() - 2);
        if (rtrim_blanks(lower_az(dn)) == lookee) {
          found = true;
          change(v.s1, v.s2, op, e, v.dataname);
        }
      }
      if (!found) warn("Unable to modify mask for " + lookee + ".  Variable not found. File: " + path + " at line " +
                       std::to_string(lineno));
      ++fieldno;
      noerr = get_fieldstr(fieldno, ',', fieldlst, 256, t1);
    }
  }
  return r;
}

void apply_iofields_option(WrfSchema& s, const WrfConfig& c, const std::string& dir, const std::vector<int>& input_bits) {
  std::string io = c.get("iofields_filename");
  if (io.size() >= 2 && (io.front() == '"' || io.front() == '\'') && io.back() == io.front())
    io = io.substr(1, io.size() - 2);
  if (rtrim_blanks(io).empty() || rtrim_blanks(io) == "NONE_SPECIFIED") return;
  if (io.front() != '/' && !dir.empty()) io = dir + "/" + io;
  const IofieldsResult r = s.apply_iofields(io, c);
  if (r.warning && !c.logical("ignore_iofields_warning"))  // modify_io_masks1: wrf_error_fatal
    throw std::runtime_error("iofields_filename: " + r.message +
                             " (WRF stops: set ignore_iofields_warning to true in namelist to ignore)");
  for (const auto& [bit, name] : r.changed)
    if (std::find(input_bits.begin(), input_bits.end(), bit) != input_bits.end())
      throw std::runtime_error("iofields_filename " + io + ": " + name + " on input stream i:" +
                               std::to_string(bit - 25) + " is not supported (CPPWRF reads its input files whole)");
}

const SchemaVar* WrfSchema::find(const std::string& varname) const {
  for (const SchemaVar& v : vars_)
    if (v.varname == varname) return &v;
  return nullptr;
}

std::vector<const SchemaMember*> WrfSchema::boundary_members(const std::string& array, const WrfConfig& c) const {
  std::vector<const SchemaMember*> out;
  std::set<std::string> seen;
  for (const SchemaMember& m : members_) {
    if (m.array != array || seen.count(m.member) || !in_use(m.conds, c)) continue;
    seen.insert(m.member);
    if (m.boundary) out.push_back(&m);
  }
  return out;
}

}  // namespace cppwrf::wrf_io
