#include "wrf_io/wrf_schema.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
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

std::string lower(std::string s) {
  for (char& c : s) c = char(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

std::string rtrim(std::string s) {
  while (!s.empty() && (s.back() == ' ' || s.back() == '\0')) s.pop_back();
  return s;
}

bool is_true(const std::string& v) {
  std::string t = lower(v);
  if (!t.empty() && t[0] == '.') t = t.substr(1);
  return !t.empty() && t[0] == 't';
}

}  // namespace

WrfConfig::WrfConfig(const Namelist& nl, const std::map<std::string, std::pair<std::string, std::string>>& defaults,
                     int dom)
    : nl_(nl), defaults_(defaults), dom_(dom) {}

std::string WrfConfig::get(const std::string& name0, int dom) const {
  if (dom == 0) dom = dom_;
  const std::string name = lower(name0);
  auto dv = derived_.find(name);
  if (dv != derived_.end()) return dv->second;
  auto d = defaults_.find(name);
  if (d != defaults_.end() && d->second.first != "derived") {
    auto r = nl_.raw(d->second.first, name, dom);
    if (r) return *r;
    return d->second.second;
  }
  return d == defaults_.end() ? std::string() : d->second.second;
}

double WrfConfig::num(const std::string& name, int dom) const {
  std::string v = get(name, dom);
  if (v.empty()) return 0;
  std::string t = lower(v);
  if (t == ".true." || t == "t" || t == ".t.") return 1;
  if (t == ".false." || t == "f" || t == ".f.") return 0;
  for (char& c : v)
    if (c == 'd' || c == 'D') c = 'e';
  return std::strtod(v.c_str(), nullptr);
}

bool WrfConfig::logical(const std::string& name, int dom) const { return is_true(get(name, dom)); }

void derive_config(WrfConfig& c, int max_dom, bool real_exe) {
  for (const char* p : {"mp_physics", "cu_physics", "ra_lw_physics", "ra_sw_physics", "sf_sfclay_physics",
                        "sf_surface_physics", "bl_pbl_physics"})
    if (c.integer(p) < 0) c.set(p, "0");
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
      s.defaults_[lower(t[1])] = {t[3], t[4]};
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
    std::string v = lower(k.val);
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

void fill_out(OutVar& o, const SchemaVar& v, const WrfConfig& c, int ide, int jde, int kde, int real_bytes) {
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
  } else if (mo == "C") {
    order = {0};
    o.memorder_att = "C  ";
  } else {
    throw std::runtime_error("schema: memory order " + mo + " of " + v.varname + " is not supported");
  }
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
}

}  // namespace

std::vector<OutVar> WrfSchema::stream_vars(const WrfConfig& c, int ide, int jde, int kde, int real_bytes, int word,
                                           int bit) const {
  const std::uint32_t mask = 1u << bit;
  auto on = [&](std::uint32_t s1, std::uint32_t s2) { return ((word == 0 ? s1 : s2) & mask) != 0; };
  std::vector<OutVar> out;
  for (const SchemaVar& v : vars_) {
    if (!on(v.s1, v.s2) || !in_use(v.conds, c)) continue;
    if (v.ntl > 0 && v.ntl / 100 != v.ntl % 100) continue;
    if (!v.array4d) {
      OutVar o;
      o.name = v.dataname;
      if (v.ntl > 0) o.name = o.name.substr(0, o.name.size() - 2);
      o.desc = v.desc; o.units = v.units; o.var = &v;
      fill_out(o, v, c, ide, jde, kde, real_bytes);
      out.push_back(std::move(o));
    } else {
      std::set<std::string> seen;
      for (const SchemaMember& m : members_) {
        if (m.array != v.varname || seen.count(m.member) || !in_use(m.conds, c)) continue;
        seen.insert(m.member);
        if (!on(m.s1, m.s2)) continue;
        OutVar o;
        o.name = m.dname; o.desc = m.desc; o.units = m.units; o.var = &v; o.member = &m;
        fill_out(o, v, c, ide, jde, kde, real_bytes);
        out.push_back(std::move(o));
      }
    }
  }
  return out;
}

bool WrfSchema::apply_iofields(const std::string& path) {
  std::ifstream f(path);
  if (!f) return false;
  std::string line;
  while (std::getline(f, line)) {
    std::vector<std::string> t = split(line, ':');
    if (t.size() < 4 || line[0] == '#') continue;
    const std::string op = rtrim(t[0]), st = lower(rtrim(t[1]));
    if ((op != "+" && op != "-") || (st != "h" && st != "i")) continue;
    // stream numbers: history = 1, auxhist n = 1 + n, input = 26, auxinput n = 26 + n; bit = number - 1
    const int idx = (st == "h" ? 0 : 25) + std::atoi(t[2].c_str());
    const int word = idx / 31;
    const std::uint32_t mask = 1u << (idx % 31);
    auto apply = [&](std::uint32_t& s1, std::uint32_t& s2) {
      std::uint32_t& w = word == 0 ? s1 : s2;
      if (op == "+") w |= mask; else w &= ~mask;
    };
    for (std::string fld : split(t[3], ',')) {
      std::string look;
      for (char ch : fld)
        if (ch != ' ' && ch != '\r') look += char(std::tolower(static_cast<unsigned char>(ch)));
      if (look.empty()) continue;
      for (SchemaVar& v : vars_) {
        if (v.array4d) continue;
        std::string dn = lower(v.dataname);
        if (v.ntl > 0 && dn.size() > 2) dn = dn.substr(0, dn.size() - 2);
        if (dn == look) apply(v.s1, v.s2);
      }
      for (SchemaMember& m : members_)
        if (lower(m.dname) == look) apply(m.s1, m.s2);
    }
  }
  return true;
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
