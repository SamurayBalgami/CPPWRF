#include "wrf_io/namelist.hpp"
#include "core/strings.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace cppwrf::wrf_io {

namespace {

// The values of one assignment (Fortran list-directed input): separated by commas and/or blanks outside
// quotes, "r*v" repeated, "r*" r null values; a null value (",," or a leading comma) is "" (the element keeps
// its value). A comma at the end is no null value.
std::vector<std::string> split_values(const std::string& s) {
  std::vector<std::string> out;
  std::string cur;
  bool have = false, comma = false, started = false;
  char quote = 0;
  auto push = [&] {
    std::string t = trim(cur);
    cur.clear();
    have = false;
    std::size_t star = t.find('*');
    if (star != std::string::npos && star > 0 && t[0] != '\'' && t[0] != '"' &&
        std::all_of(t.begin(), t.begin() + long(star), [](char c) { return std::isdigit(static_cast<unsigned char>(c)); })) {
      const int n = std::atoi(t.substr(0, star).c_str());
      const std::string v = t.substr(star + 1);
      for (int i = 0; i < n; ++i) out.push_back(v);
      return;
    }
    out.push_back(t);
  };
  for (char c : s) {
    if (quote) {
      cur += c;
      if (c == quote) quote = 0;
    } else if (c == ',') {
      if (have) push();
      else if (comma || !started) out.push_back("");
      comma = started = true;
    } else if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
      if (have) {
        push();
        comma = false;
      }
    } else {
      if (c == '\'' || c == '"') quote = c;
      cur += c;
      have = started = true;
      comma = false;
    }
  }
  if (have) push();
  return out;
}

}  // namespace

Namelist Namelist::read(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("namelist: cannot open " + path);
  std::stringstream ss;
  ss << f.rdbuf();
  return parse(ss.str());
}

Namelist Namelist::parse(const std::string& text) {
  Namelist nl;
  // strip comments (! outside quotes)
  std::string t;
  {
    char quote = 0;
    bool comment = false;
    for (char c : text) {
      if (c == '\n') comment = false;
      if (comment) continue;
      if (quote) {
        if (c == quote) quote = 0;
      } else if (c == '\'' || c == '"') {
        quote = c;
      } else if (c == '!') {
        comment = true;
        continue;
      }
      t += c;
    }
  }
  std::size_t pos = 0;
  std::set<std::string> seen;
  while (true) {
    std::size_t amp = t.find('&', pos);
    if (amp == std::string::npos) break;
    std::size_t e = amp + 1;
    while (e < t.size() && (std::isalnum(static_cast<unsigned char>(t[e])) || t[e] == '_')) ++e;
    std::string group = to_lower(t.substr(amp + 1, e - amp - 1));
    // body ends at a '/' outside quotes
    std::size_t end = e;
    char quote = 0;
    for (; end < t.size(); ++end) {
      char c = t[end];
      if (quote) {
        if (c == quote) quote = 0;
      } else if (c == '\'' || c == '"') {
        quote = c;
      } else if (c == '/') {
        break;
      }
    }
    std::string body = t.substr(e, end - e);
    // WRF reads a group with REWIND + READ(NML=...): only its first occurrence counts
    const bool first = seen.insert(group).second;
    auto& grp = nl.g_[group];
    // assignments: key = values ... up to the next "key =" (key: a name, optionally with an element designator)
    std::vector<std::pair<std::size_t, std::size_t>> eqs;  // (key start, '=' position)
    quote = 0;
    for (std::size_t i = 0; i < body.size(); ++i) {
      char c = body[i];
      if (quote) {
        if (c == quote) quote = 0;
        continue;
      }
      if (c == '\'' || c == '"') {
        quote = c;
        continue;
      }
      if (c == '=') {
        std::size_t k = i;
        while (k > 0 && std::isspace(static_cast<unsigned char>(body[k - 1]))) --k;
        if (k > 0 && body[k - 1] == ')') {  // name(i) or name(i:j)
          const std::size_t open = body.rfind('(', k - 1);
          if (open != std::string::npos) k = open;
          while (k > 0 && std::isspace(static_cast<unsigned char>(body[k - 1]))) --k;
        }
        std::size_t ks = k;
        while (ks > 0 && (std::isalnum(static_cast<unsigned char>(body[ks - 1])) || body[ks - 1] == '_')) --ks;
        eqs.emplace_back(ks, i);
      }
    }
    for (std::size_t n = 0; first && n < eqs.size(); ++n) {
      std::string key;
      for (char ch : body.substr(eqs[n].first, eqs[n].second - eqs[n].first))
        if (!std::isspace(static_cast<unsigned char>(ch))) key += char(std::tolower(static_cast<unsigned char>(ch)));
      std::size_t base = 0;  // the first element assigned (0-based)
      if (const std::size_t lp = key.find('('); lp != std::string::npos) {
        const std::string d = key.substr(lp + 1, key.size() - lp - 2);
        const std::size_t colon = d.find(':');
        const std::string lo = d.substr(0, colon);
        if (key.back() != ')' || lo.empty() || lo.find_first_not_of("0123456789") != std::string::npos ||
            (colon != std::string::npos && d.substr(colon + 1).find_first_not_of("0123456789") != std::string::npos) ||
            std::atoi(lo.c_str()) < 1)
          throw std::runtime_error("namelist: unsupported element designator in &" + group + " " + key);
        base = std::size_t(std::atoi(lo.c_str()) - 1);
        key = key.substr(0, lp);
      }
      std::size_t vs = eqs[n].second + 1, ve = n + 1 < eqs.size() ? eqs[n + 1].first : body.size();
      // element by element: an element without a (non-null) value keeps its earlier one
      std::vector<std::string>& dst = grp[key];
      const std::vector<std::string> vals = split_values(body.substr(vs, ve - vs));
      if (dst.size() < base + vals.size()) dst.resize(base + vals.size());
      for (std::size_t k = 0; k < vals.size(); ++k)
        if (!vals[k].empty()) dst[base + k] = vals[k];
    }
    pos = end + 1;
  }
  return nl;
}

bool Namelist::has(const std::string& group, const std::string& key) const {
  auto g = g_.find(to_lower(group));
  return g != g_.end() && g->second.count(to_lower(key));
}

std::vector<std::string> Namelist::values(const std::string& group, const std::string& key) const {
  auto g = g_.find(to_lower(group));
  if (g == g_.end()) return {};
  auto k = g->second.find(to_lower(key));
  return k == g->second.end() ? std::vector<std::string>{} : k->second;
}

std::optional<std::string> Namelist::raw(const std::string& group, const std::string& key, int dom) const {
  std::vector<std::string> v = values(group, key);
  if (v.empty()) return std::nullopt;
  const std::string& r = dom - 1 < int(v.size()) ? v[std::size_t(dom - 1)] : v.back();
  if (r.empty()) return std::nullopt;  // a null value: the default
  return r;
}

int Namelist::get_int(const std::string& group, const std::string& key, int def, int dom) const {
  auto r = raw(group, key, dom);
  return r ? std::atoi(r->c_str()) : def;
}

static std::string fortran_real(std::string s) {
  for (char& c : s)
    if (c == 'd' || c == 'D') c = 'e';
  return s;
}

double Namelist::get_real(const std::string& group, const std::string& key, double def, int dom) const {
  auto r = raw(group, key, dom);
  return r ? std::strtod(fortran_real(*r).c_str(), nullptr) : def;
}

float Namelist::get_real_f(const std::string& group, const std::string& key, float def, int dom) const {
  auto r = raw(group, key, dom);
  return r ? std::strtof(fortran_real(*r).c_str(), nullptr) : def;
}

bool Namelist::get_bool(const std::string& group, const std::string& key, bool def, int dom) const {
  auto r = raw(group, key, dom);
  if (!r) return def;
  std::string v = to_lower(*r);
  if (!v.empty() && v[0] == '.') v = v.substr(1);
  return !v.empty() && v[0] == 't';
}

std::string Namelist::get_str(const std::string& group, const std::string& key, const std::string& def, int dom) const {
  auto r = raw(group, key, dom);
  if (!r) return def;
  std::string v = *r;
  if (v.size() >= 2 && (v.front() == '\'' || v.front() == '"')) v = v.substr(1, v.size() - 2);
  return v;
}

}  // namespace cppwrf::wrf_io
