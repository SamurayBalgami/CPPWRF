#include "wrf_io/namelist.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace cppwrf::wrf_io {

namespace {

std::string lower(std::string s) {
  for (char& c : s) c = char(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

std::string trim(const std::string& s) {
  std::size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
  return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
}

// Splits the value list of one assignment at commas/blanks outside quotes; expands n*value.
std::vector<std::string> split_values(const std::string& s) {
  std::vector<std::string> out;
  std::string cur;
  char quote = 0;
  auto flush = [&] {
    std::string t = trim(cur);
    cur.clear();
    if (t.empty()) return;
    std::size_t star = t.find('*');
    if (star != std::string::npos && t[0] != '\'' && t[0] != '"' &&
        std::all_of(t.begin(), t.begin() + long(star), [](char c) { return std::isdigit(static_cast<unsigned char>(c)); })) {
      int n = std::atoi(t.substr(0, star).c_str());
      std::string v = t.substr(star + 1);
      for (int i = 0; i < n; ++i) out.push_back(v);
      return;
    }
    out.push_back(t);
  };
  for (char c : s) {
    if (quote) {
      cur += c;
      if (c == quote) quote = 0;
    } else if (c == '\'' || c == '"') {
      quote = c;
      cur += c;
    } else if (c == ',' || c == ' ' || c == '\t' || c == '\n' || c == '\r') {
      flush();
    } else {
      cur += c;
    }
  }
  flush();
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
  while (true) {
    std::size_t amp = t.find('&', pos);
    if (amp == std::string::npos) break;
    std::size_t e = amp + 1;
    while (e < t.size() && (std::isalnum(static_cast<unsigned char>(t[e])) || t[e] == '_')) ++e;
    std::string group = lower(t.substr(amp + 1, e - amp - 1));
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
    auto& grp = nl.g_[group];
    // assignments: key = values ... up to the next "key ="
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
        std::size_t ks = k;
        while (ks > 0 && (std::isalnum(static_cast<unsigned char>(body[ks - 1])) || body[ks - 1] == '_' ||
                          body[ks - 1] == '(' || body[ks - 1] == ')'))
          --ks;
        eqs.emplace_back(ks, i);
      }
    }
    for (std::size_t n = 0; n < eqs.size(); ++n) {
      std::string key = lower(trim(body.substr(eqs[n].first, eqs[n].second - eqs[n].first)));
      std::size_t vs = eqs[n].second + 1, ve = n + 1 < eqs.size() ? eqs[n + 1].first : body.size();
      grp[key] = split_values(body.substr(vs, ve - vs));
    }
    pos = end + 1;
  }
  return nl;
}

bool Namelist::has(const std::string& group, const std::string& key) const {
  auto g = g_.find(lower(group));
  return g != g_.end() && g->second.count(lower(key));
}

std::vector<std::string> Namelist::values(const std::string& group, const std::string& key) const {
  auto g = g_.find(lower(group));
  if (g == g_.end()) return {};
  auto k = g->second.find(lower(key));
  return k == g->second.end() ? std::vector<std::string>{} : k->second;
}

std::optional<std::string> Namelist::raw(const std::string& group, const std::string& key, int dom) const {
  std::vector<std::string> v = values(group, key);
  if (v.empty()) return std::nullopt;
  if (dom - 1 < int(v.size())) return v[std::size_t(dom - 1)];
  return v.back();
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
  std::string v = lower(*r);
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
