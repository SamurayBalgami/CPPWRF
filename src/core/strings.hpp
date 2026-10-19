#pragma once
// String helpers of the configuration, namelist and table readers.

#include <cctype>
#include <string>
#include <string_view>

namespace cppwrf {

// s without leading and trailing blanks, tabs and line ends
inline std::string trim(std::string_view s) {
  const std::size_t a = s.find_first_not_of(" \t\r\n");
  if (a == std::string_view::npos) return {};
  const std::size_t b = s.find_last_not_of(" \t\r\n");
  return std::string(s.substr(a, b - a + 1));
}

inline std::string to_upper(std::string_view s) {
  std::string r(s);
  for (char& c : r) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return r;
}

inline std::string to_lower(std::string_view s) {
  std::string r(s);
  for (char& c : r) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return r;
}

}  // namespace cppwrf
