#pragma once
// Fortran (gfortran) numeric input of the radiation tables WRF reads with formatted or list-directed READs:
// a value the Fortran read rejects is reported (false / -1) so the caller fails as WRF does.
#include <cstdlib>
#include <string>

namespace cppwrf::phys {

// [sign] digits [. digits] [exponent], exponent = (E|D|Q|e|d|q) [sign] digits or sign digits, without blanks;
// implied_d: the Fw.d edit descriptor's d (digits without a decimal point: the last d are the fraction), 0 in
// list-directed input. Converted as gfortran does (strtod of the digits: correctly rounded).
inline bool fortran_real(const std::string& s, int implied_d, double& v) {
  const std::size_t n = s.size();
  std::size_t i = 0, nd = 0;
  bool dp = false;
  std::string m;
  if (i < n && (s[i] == '+' || s[i] == '-')) m += s[i++];
  for (; i < n; ++i) {
    if (s[i] >= '0' && s[i] <= '9') {
      m += s[i];
      ++nd;
    } else if (s[i] == '.' && !dp) {
      m += '.';
      dp = true;
    } else {
      break;
    }
  }
  if (nd == 0) return false;
  long e = 0;
  if (i < n) {
    switch (s[i]) {
      case 'E': case 'e': case 'D': case 'd': case 'Q': case 'q': ++i; break;
      case '+': case '-': break;
      default: return false;
    }
    bool neg = false;
    if (i < n && (s[i] == '+' || s[i] == '-')) neg = s[i++] == '-';
    if (i == n) return false;
    for (; i < n; ++i) {
      if (s[i] < '0' || s[i] > '9') return false;
      if (e < 100000) e = e * 10 + (s[i] - '0');
    }
    if (neg) e = -e;
  }
  if (!dp) e -= implied_d;
  v = std::strtod((m + "e" + std::to_string(e)).c_str(), nullptr);
  return true;
}

// formatted input (BLANK = 'NULL', PAD = 'YES'): the field's blanks are ignored, a blank field (or one past the
// end of the record) is 0
inline std::string fortran_field(const std::string& rec, std::size_t a, std::size_t w) {
  std::string s;
  for (std::size_t i = a; i < a + w && i < rec.size(); ++i)
    if (rec[i] != ' ') s += rec[i];
  return s;
}

// Iw
inline bool fortran_i_field(const std::string& rec, std::size_t a, std::size_t w, long& v) {
  const std::string s = fortran_field(rec, a, w);
  v = 0;
  if (s.empty()) return true;
  std::size_t i = (s[0] == '+' || s[0] == '-') ? 1 : 0;
  if (i == s.size()) return false;
  for (std::size_t k = i; k < s.size(); ++k)
    if (s[k] < '0' || s[k] > '9') return false;
  v = std::strtol(s.c_str(), nullptr, 10);
  return true;
}

// Fw.d
inline bool fortran_f_field(const std::string& rec, std::size_t a, std::size_t w, int d, double& v) {
  const std::string s = fortran_field(rec, a, w);
  v = 0.0;
  return s.empty() || fortran_real(s, d, v);
}

// READ (unit, *) x of one real: the record's first value (r*c: c; the rest of the record is skipped).
// 1: read; 0: blank record (the read goes on with the next one); -1: not a number, or a null value / slash
// (x left undefined)
inline int fortran_list_real(const std::string& rec, double& v) {
  const std::size_t a = rec.find_first_not_of(" \t\r");
  if (a == std::string::npos) return 0;
  std::string t = rec.substr(a, rec.find_first_of(" \t\r,/", a) - a);
  const std::size_t star = t.find('*');
  if (star != std::string::npos) {
    if (star == 0) return -1;
    for (std::size_t k = 0; k < star; ++k)
      if (t[k] < '0' || t[k] > '9') return -1;
    if (std::strtol(t.substr(0, star).c_str(), nullptr, 10) <= 0) return -1;
    t = t.substr(star + 1);
  }
  if (t.empty()) return -1;
  return fortran_real(t, 0, v) ? 1 : -1;
}

}  // namespace cppwrf::phys
