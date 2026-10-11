#ifndef CPPWRF_IO_CONFIG_HPP
#define CPPWRF_IO_CONFIG_HPP

#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace cppwrf::io {

struct ConfigError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

class Config {
 public:
  static Config load(const std::string& path, const std::vector<std::string>& overrides = {});
  static Config from_string(const std::string& text, const std::vector<std::string>& overrides = {});
  bool has(const std::string& key) const;
  std::string get_string(const std::string& key) const;                  // throws if missing
  std::string get_string(const std::string& key, const std::string& def) const;
  long long   get_int(const std::string& key) const;     long long get_int(const std::string& key, long long def) const;
  double      get_double(const std::string& key) const;  double    get_double(const std::string& key, double def) const;
  bool        get_bool(const std::string& key) const;    bool      get_bool(const std::string& key, bool def) const; // true/false/1/0/yes/no
  std::vector<double> get_double_list(const std::string& key) const;     // comma separated
  std::vector<std::string> keys() const;                                 // sorted
  std::string to_json() const;  // {"domain":{"nx":"512",...},...}; values as JSON strings

 private:
  std::map<std::string, std::string> data_;
  static Config parse(const std::string& text, const std::string& source_name, const std::vector<std::string>& overrides);
  void apply_overrides(const std::vector<std::string>& overrides);
};

} // namespace cppwrf::io

#endif // CPPWRF_IO_CONFIG_HPP
