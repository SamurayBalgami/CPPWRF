#include "io/config.hpp"
#include "io/output.hpp"

#include <charconv>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string_view>

namespace cppwrf::io {

namespace {

std::string trim(std::string_view s) {
  size_t start = 0;
  while (start < s.size() && (s[start] == ' ' || s[start] == '\t' || s[start] == '\r' || s[start] == '\n')) {
    ++start;
  }
  size_t end = s.size();
  while (end > start && (s[end - 1] == ' ' || s[end - 1] == '\t' || s[end - 1] == '\r' || s[end - 1] == '\n')) {
    --end;
  }
  return std::string(s.substr(start, end - start));
}

bool is_valid_ident(std::string_view s) {
  if (s.empty()) return false;
  for (char c : s) {
    if (!((c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_')) {
      return false;
    }
  }
  return true;
}

long long parse_int_strict(const std::string& key, const std::string& val) {
  if (val.empty()) {
    throw ConfigError("Key '" + key + "': empty integer value");
  }
  std::string_view sv = val;
  if (sv[0] == '+') sv.remove_prefix(1);
  if (sv.empty()) {
    throw ConfigError("Key '" + key + "': invalid integer value '" + val + "'");
  }
  long long res = 0;
  auto [ptr, ec] = std::from_chars(sv.data(), sv.data() + sv.size(), res);
  if (ec != std::errc{} || ptr != sv.data() + sv.size()) {
    throw ConfigError("Key '" + key + "': invalid integer value '" + val + "'");
  }
  return res;
}

double parse_double_strict(const std::string& key, const std::string& val) {
  if (val.empty()) {
    throw ConfigError("Key '" + key + "': empty double value");
  }
  std::string_view sv = val;
  if (sv[0] == '+') sv.remove_prefix(1);
  if (sv.empty()) {
    throw ConfigError("Key '" + key + "': invalid double value '" + val + "'");
  }
  double res = 0.0;
  auto [ptr, ec] = std::from_chars(sv.data(), sv.data() + sv.size(), res);
  if (ec != std::errc{} || ptr != sv.data() + sv.size() || std::isnan(res) || std::isinf(res)) {
    throw ConfigError("Key '" + key + "': invalid double value '" + val + "'");
  }
  return res;
}

bool parse_bool_strict(const std::string& key, const std::string& val) {
  std::string lower;
  lower.reserve(val.size());
  for (char c : val) {
    if (c >= 'A' && c <= 'Z') {
      lower.push_back(static_cast<char>(c - 'A' + 'a'));
    } else {
      lower.push_back(c);
    }
  }
  if (lower == "true" || lower == "1" || lower == "yes") return true;
  if (lower == "false" || lower == "0" || lower == "no") return false;
  throw ConfigError("Key '" + key + "': invalid boolean value '" + val + "'");
}

std::vector<double> parse_double_list_strict(const std::string& key, const std::string& val) {
  std::vector<double> list;
  std::string trimmed_val = trim(val);
  if (trimmed_val.empty()) return list;

  size_t start = 0;
  while (start <= trimmed_val.size()) {
    size_t comma = trimmed_val.find(',', start);
    std::string token;
    if (comma == std::string::npos) {
      token = trim(trimmed_val.substr(start));
      start = trimmed_val.size() + 1;
    } else {
      token = trim(trimmed_val.substr(start, comma - start));
      start = comma + 1;
    }
    if (token.empty()) {
      throw ConfigError("Key '" + key + "': invalid empty item in double list");
    }
    list.push_back(parse_double_strict(key, token));
  }
  return list;
}

} // namespace

Config Config::parse(const std::string& text, const std::string& source_name, const std::vector<std::string>& overrides) {
  Config cfg;
  std::stringstream ss(text);
  std::string raw_line;
  std::string current_section;
  int lineno = 0;

  while (std::getline(ss, raw_line)) {
    ++lineno;
    std::string line = raw_line;
    // Strip comments (# or ;)
    // Inline comment allowed after whitespace
    for (size_t i = 0; i < line.size(); ++i) {
      if (line[i] == '#' || line[i] == ';') {
        if (i == 0 || line[i - 1] == ' ' || line[i - 1] == '\t') {
          line.erase(i);
          break;
        }
      }
    }

    std::string trimmed = trim(line);
    if (trimmed.empty()) {
      continue;
    }

    if (trimmed.front() == '[') {
      if (trimmed.back() != ']') {
        throw ConfigError(source_name + ":" + std::to_string(lineno) + ": unclosed section bracket in '" + trimmed + "'");
      }
      std::string sec = trimmed.substr(1, trimmed.size() - 2);
      if (!is_valid_ident(sec)) {
        throw ConfigError(source_name + ":" + std::to_string(lineno) + ": invalid section header '" + sec + "'");
      }
      current_section = sec;
    } else {
      auto eq_pos = trimmed.find('=');
      if (eq_pos == std::string::npos) {
        throw ConfigError(source_name + ":" + std::to_string(lineno) + ": syntax error: missing '=' in '" + trimmed + "'");
      }
      if (current_section.empty()) {
        throw ConfigError(source_name + ":" + std::to_string(lineno) + ": key before any section header: '" + trimmed + "'");
      }
      std::string key_part = trim(trimmed.substr(0, eq_pos));
      if (!is_valid_ident(key_part)) {
        throw ConfigError(source_name + ":" + std::to_string(lineno) + ": invalid key name '" + key_part + "'");
      }
      std::string full_key = current_section + "." + key_part;
      if (cfg.data_.find(full_key) != cfg.data_.end()) {
        throw ConfigError(source_name + ":" + std::to_string(lineno) + ": duplicate key '" + full_key + "'");
      }
      std::string val_part = trim(trimmed.substr(eq_pos + 1));
      cfg.data_[full_key] = val_part;
    }
  }

  cfg.apply_overrides(overrides);
  return cfg;
}

void Config::apply_overrides(const std::vector<std::string>& overrides) {
  for (const auto& ov : overrides) {
    auto eq_pos = ov.find('=');
    if (eq_pos == std::string::npos) {
      throw ConfigError("Invalid override syntax (missing '='): '" + ov + "'");
    }
    std::string full_key = trim(ov.substr(0, eq_pos));
    std::string val = trim(ov.substr(eq_pos + 1));
    auto dot_pos = full_key.find('.');
    if (dot_pos == std::string::npos || full_key.find('.', dot_pos + 1) != std::string::npos) {
      throw ConfigError("Invalid override key '" + full_key + "' (must be section.key)");
    }
    std::string sec = full_key.substr(0, dot_pos);
    std::string key = full_key.substr(dot_pos + 1);
    if (!is_valid_ident(sec) || !is_valid_ident(key)) {
      throw ConfigError("Invalid override key '" + full_key + "' (invalid characters)");
    }
    data_[full_key] = val;
  }
}

Config Config::load(const std::string& path, const std::vector<std::string>& overrides) {
  std::ifstream in(path, std::ios::in | std::ios::binary);
  if (!in.is_open()) {
    throw ConfigError("Failed to open file: " + path);
  }
  std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  return parse(text, path, overrides);
}

Config Config::from_string(const std::string& text, const std::vector<std::string>& overrides) {
  return parse(text, "<string>", overrides);
}

bool Config::has(const std::string& key) const {
  return data_.find(key) != data_.end();
}

std::string Config::get_string(const std::string& key) const {
  auto it = data_.find(key);
  if (it == data_.end()) {
    throw ConfigError("Key not found: '" + key + "'");
  }
  return it->second;
}

std::string Config::get_string(const std::string& key, const std::string& def) const {
  auto it = data_.find(key);
  return (it != data_.end()) ? it->second : def;
}

long long Config::get_int(const std::string& key) const {
  auto it = data_.find(key);
  if (it == data_.end()) {
    throw ConfigError("Key not found: '" + key + "'");
  }
  return parse_int_strict(key, it->second);
}

long long Config::get_int(const std::string& key, long long def) const {
  auto it = data_.find(key);
  if (it == data_.end()) {
    return def;
  }
  return parse_int_strict(key, it->second);
}

double Config::get_double(const std::string& key) const {
  auto it = data_.find(key);
  if (it == data_.end()) {
    throw ConfigError("Key not found: '" + key + "'");
  }
  return parse_double_strict(key, it->second);
}

double Config::get_double(const std::string& key, double def) const {
  auto it = data_.find(key);
  if (it == data_.end()) {
    return def;
  }
  return parse_double_strict(key, it->second);
}

bool Config::get_bool(const std::string& key) const {
  auto it = data_.find(key);
  if (it == data_.end()) {
    throw ConfigError("Key not found: '" + key + "'");
  }
  return parse_bool_strict(key, it->second);
}

bool Config::get_bool(const std::string& key, bool def) const {
  auto it = data_.find(key);
  if (it == data_.end()) {
    return def;
  }
  return parse_bool_strict(key, it->second);
}

std::vector<double> Config::get_double_list(const std::string& key) const {
  auto it = data_.find(key);
  if (it == data_.end()) {
    throw ConfigError("Key not found: '" + key + "'");
  }
  return parse_double_list_strict(key, it->second);
}

std::vector<std::string> Config::keys() const {
  std::vector<std::string> res;
  res.reserve(data_.size());
  for (const auto& [k, v] : data_) {
    res.push_back(k);
  }
  return res;
}

std::string Config::to_json() const {
  std::map<std::string, std::map<std::string, std::string>> sections;
  for (const auto& [full_key, val] : data_) {
    auto dot_pos = full_key.find('.');
    if (dot_pos != std::string::npos) {
      sections[full_key.substr(0, dot_pos)][full_key.substr(dot_pos + 1)] = val;
    }
  }

  JsonWriter writer;
  writer.start_object();
  for (const auto& [sec, kvs] : sections) {
    writer.key(sec);
    writer.start_object();
    for (const auto& [k, v] : kvs) {
      writer.key(k);
      writer.write_string(v);
    }
    writer.end_object();
  }
  writer.end_object();
  return writer.str();
}

} // namespace cppwrf::io
