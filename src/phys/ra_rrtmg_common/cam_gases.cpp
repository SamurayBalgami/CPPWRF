#include "phys/ra_rrtmg_common/cam_gases.hpp"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <stdexcept>

namespace cppwrf::phys {

namespace {
// Fortran F/I edit descriptor read of columns [a, a+w) (blanks ignored, empty = 0)
double field(const std::string& line, std::size_t a, std::size_t w) {
  if (a >= line.size()) return 0.0;
  std::string s = line.substr(a, w);
  s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
  return s.empty() ? 0.0 : std::strtod(s.c_str(), nullptr);
}
bool leap(int y) { return (y % 4 == 0 && y % 100 != 0) || y % 400 == 0; }

// valid_years (1-based indices as in WRF)
void valid_years(int yearIN, const std::vector<double>& gas, int tot_years, int& yr1, int& yr2) {
  auto g = [&](int i) { return gas[std::size_t(i - 1)]; };
  int yr_loc = yearIN;
  yr2 = yr_loc;
  yr1 = yr_loc - 1;
  int below = 0;
  for (int i = 1; i <= yr_loc - 1; ++i) below += g(i) > 0.;
  if (below == 0) {
    for (int i = yr_loc; i <= tot_years - 1; ++i)
      if (g(i) > 0.) { yr1 = i; break; }
    for (int i = yr1 + 1; i <= tot_years; ++i)
      if (g(i) > 0.) { yr2 = i; break; }
  } else {
    if (g(yr_loc) < 0.) {
      bool any_above = false;
      for (int i = yr_loc; i <= tot_years; ++i) any_above |= g(i) > 0;
      if (any_above) {
        for (int i = yr_loc + 1; i <= tot_years; ++i)
          if (g(i) > 0) { yr2 = i; break; }
      } else {
        for (int i = yr_loc - 1; i >= 2; --i)
          if (g(i) > 0) { yr2 = i; break; }
      }
    }
    yr_loc = std::min(yr_loc - 1, yr2 - 1);
    if (g(yr_loc) < 0.) {
      for (int i = yr_loc - 1; i >= 1; --i)
        if (g(i) > 0) { yr1 = i; break; }
    } else {
      yr1 = yr_loc;
    }
  }
}

double interpolate(int yr, double julian, int yeari, int juli, int yr1, int yr2, int yearf, int julf,
                   const std::vector<double>& gas) {
  const int minyear = std::min(yr, yeari), maxyear = std::max(yr, yearf);
  double fact2 = juli, fact1 = julf, x = julian;
  for (int ny = minyear; ny <= maxyear - 1; ++ny) {
    const int nday = leap(ny) ? 366 : 365;
    if (ny < yeari) fact2 = fact2 + nday;
    if (ny < yr) x = x + nday;
    if (ny < yearf) fact1 = fact1 + nday;
  }
  const double deltat = fact1 - fact2;
  fact1 = (fact1 - x) / deltat;
  fact2 = (x - fact2) / deltat;
  double v = gas[std::size_t(yr1 - 1)] * fact1 + gas[std::size_t(yr2 - 1)] * fact2;
  if (v < 0.) v = -99999.;
  return v;
}
}  // namespace

void CamGases::load(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("CLWRF: '" + path + "' does not exist");
  std::string line;
  std::getline(f, line);
  std::getline(f, line);
  // FORMAT (I4, 1x, F8.3, 1x, 4(F10.3, 1x)); WRF stops at the first unreadable record
  while (std::getline(f, line)) {
    // a blank record reads as zeros in Fortran (no error)
    const std::string ys = line.substr(0, std::min<std::size_t>(4, line.size()));
    long y = 0;
    if (ys.find_first_not_of(" \t\r") != std::string::npos) {
      char* end = nullptr;
      y = std::strtol(ys.c_str(), &end, 10);
      if (end == ys.c_str()) break;
    }
    yr.push_back(int(y));
    co2.push_back(field(line, 5, 8));
    n2o.push_back(field(line, 14, 10));
    ch4.push_back(field(line, 25, 10));
    cfc11.push_back(field(line, 36, 10));
    cfc12.push_back(field(line, 47, 10));
    // juldata = sum(mondayi(1:6)) + mondayi(7)/2 - 0.5 (mid-year, mondata = 6)
    const double s = 151.0 + (leap(int(y)) ? 1.0 : 0.0);
    jul.push_back(s + 30.0 / 2. - 0.5);
  }
}

void CamGases::interpolate(int year, double julian, double& co2v, double& n2ov, double& ch4v, double& cfc11v,
                           double& cfc12v) const {
  const int max_years = int(yr.size());
  constexpr int cyr = 1800;
  auto yrd = [&](int i) { return i <= max_years ? yr[std::size_t(i - 1)] : 0; };
  auto jld = [&](int i) { return i <= max_years ? jul[std::size_t(i - 1)] : 0.0; };
  int found = 0, yearIN = 0, iyear = 1;
  while (found == 0 && iyear <= cyr) {
    if (yrd(iyear) > year) { yearIN = iyear; found = 1; }
    else if (yrd(iyear) == year && jld(iyear) > julian) { yearIN = iyear; found = 1; }
    ++iyear;
  }
  if (iyear >= max_years) { yearIN = max_years - 1; found = 1; }
  int nyrm = 0, njulm = 0, nyrp = 0, njulp = 0;
  if (found != 0) {
    if (yearIN == 1) yearIN = yearIN + 1;
    nyrm = yrd(yearIN - 1); njulm = int(jld(yearIN - 1));
    nyrp = yrd(yearIN); njulp = int(jld(yearIN));
  }
  co2v = n2ov = ch4v = cfc11v = cfc12v = -9999.999;
  if (found != 0) {
    int tot_valid = 0;
    for (int i = 0; i < max_years; ++i) tot_valid += co2[std::size_t(i)] > 0;
    if (tot_valid >= 2) {
      int yr1, yr2;
      valid_years(yearIN, co2, max_years, yr1, yr2);
      nyrm = yrd(yr1); njulm = int(jld(yr1)); nyrp = yrd(yr2); njulp = int(jld(yr2));
      co2v = phys::interpolate(year, julian, nyrm, njulm, yr1, yr2, nyrp, njulp, co2);
      n2ov = phys::interpolate(year, julian, nyrm, njulm, yr1, yr2, nyrp, njulp, n2o);
      ch4v = phys::interpolate(year, julian, nyrm, njulm, yr1, yr2, nyrp, njulp, ch4);
      cfc11v = phys::interpolate(year, julian, nyrm, njulm, yr1, yr2, nyrp, njulp, cfc11);
      cfc12v = phys::interpolate(year, julian, nyrm, njulm, yr1, yr2, nyrp, njulp, cfc12);
    }
  }
  if (co2v < 0. || n2ov < 0. || ch4v < 0. || cfc11v < 0. || cfc12v < 0. || found == 0) {  // orig_val, RRTMG
    co2v = 379.e-6; n2ov = 319e-9; ch4v = 1774.e-9; cfc11v = 0.251e-9; cfc12v = 0.538e-9;
  } else {
    if (co2v < 270.) co2v = 270.;
    if (n2ov < 270.) n2ov = 270.;
    if (ch4v < 700.) ch4v = 700.;
    co2v = co2v * 1.e-06;
    n2ov = n2ov * 1.e-09;
    ch4v = ch4v * 1.e-09;
    cfc11v = cfc11v * 1.e-12;
    cfc12v = cfc12v * 1.e-12;
  }
}

}  // namespace cppwrf::phys
