#include "phys/ra_rrtmg_common/cam_gases.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>

#include "phys/ra_rrtmg_common/fortran_read.hpp"

namespace cppwrf::phys {

namespace {
constexpr int cyr = 1800;  // maximum number of records (module_ra_clWRF_support.F)
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
  // two READ (FMT='(1X)') of the header: a missing header record ends WRF with a Fortran read error
  for (int h = 0; h < 2; ++h)
    if (!std::getline(f, line)) throw std::runtime_error("CLWRF: end of file in the header of '" + path + "'");
  // FORMAT (I4, 1x, F8.3, 1x, 4(F10.3, 1x)) up to the end of the file; a record the read rejects is fatal
  // ("Not normal ending of CAMtr_volume_mixing_ratio file")
  static constexpr std::size_t at[5] = {5, 14, 25, 36, 47}, width[5] = {8, 10, 10, 10, 10};
  std::vector<double>* gas[5] = {&co2, &n2o, &ch4, &cfc11, &cfc12};
  while (std::getline(f, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    // yrdata etc. hold cyr = 1800 records, mondata(idata) is set after the end-of-file read too
    if (yr.size() + 1 >= std::size_t(cyr))
      throw std::runtime_error("CLWRF: '" + path + "' has more records than WRF's arrays hold (1799)");
    long y;
    double v[5];
    bool ok = fortran_i_field(line, 0, 4, y);
    for (int n = 0; n < 5 && ok; ++n) ok = fortran_f_field(line, at[n], width[n], 3, v[n]);
    if (!ok) throw std::runtime_error("CLWRF: Not normal ending of '" + path + "' (record '" + line + "')");
    yr.push_back(int(y));
    for (int n = 0; n < 5; ++n) gas[n]->push_back(v[n]);
    // juldata = sum(mondayi(1:6)) + mondayi(7)/2 - 0.5 (mid-year, mondata = 6)
    const double s = 151.0 + (leap(int(y)) ? 1.0 : 0.0);
    jul.push_back(s + 30.0 / 2. - 0.5);
  }
}

void CamGases::interpolate(int year, double julian, double& co2v, double& n2ov, double& ch4v, double& cfc11v,
                           double& cfc12v) const {
  const int max_years = int(yr.size());
  // records past the table are 0 (WRF's SAVE arrays); fewer than 2 records give yearIN < 1, whose values WRF
  // reads before its arrays but does not use (tot_valid < 2: orig_val)
  auto yrd = [&](int i) { return i >= 1 && i <= max_years ? yr[std::size_t(i - 1)] : 0; };
  auto jld = [&](int i) { return i >= 1 && i <= max_years ? jul[std::size_t(i - 1)] : 0.0; };
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
