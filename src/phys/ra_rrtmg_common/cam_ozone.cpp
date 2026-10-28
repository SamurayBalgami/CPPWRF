#include "phys/ra_rrtmg_common/cam_ozone.hpp"

#include <fstream>
#include <stdexcept>

#include "phys/ra_rrtmg_common/fortran_read.hpp"

namespace cppwrf::phys {

namespace {
// oznini's READ (unit,*) x per value: the first value of the next non-blank record, the rest of the record
// skipped; a value the read rejects and the end of the file are Fortran read errors in WRF
std::vector<double> read_list(const std::string& path, std::size_t n) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("CAM ozone: cannot open " + path);
  std::vector<double> v(n);
  std::string line;
  for (std::size_t i = 0; i < n; ++i) {
    int st = 0;
    while (st == 0 && std::getline(f, line))
      if ((st = fortran_list_real(line, v[i])) < 0)
        throw std::runtime_error("CAM ozone: bad value in " + path + ": '" + line + "'");
    if (st == 0) throw std::runtime_error("CAM ozone: short file " + path);
  }
  return v;
}
}  // namespace

void CamOzone::load(const std::string& dir) {
  pin = read_list(dir + "/ozone_plev.formatted", levsiz);
  for (double& p : pin) p = p * 100.;
  // ozn_p_int brackets each model pressure between increasing ozone levels (non-monotonic levels: WRF's
  // row search may differ or end with "Bad ozone data: non-monotonicity suspected")
  for (int k = 1; k < levsiz; ++k)
    if (!(pin[k - 1] < pin[k])) throw std::runtime_error("CAM ozone: ozone_plev.formatted is not increasing");
  lat = read_list(dir + "/ozone_lat.formatted", latsiz);
  // READ order: do m = 2, 13; do j = 1, latsiz; do k = 1, levsiz (lonsiz = 1)
  ozmixin = read_list(dir + "/ozone.formatted", std::size_t(nmon) * latsiz * levsiz);
}

void CamOzone::column(double xlat, double* ozmixm) const {
  // lin_interpol2(lat_ozone, ozmixin(1,k,:,m), xlat)
  const int n = latsiz;
  int k;
  if (xlat <= lat[0]) k = 1;
  else if (xlat >= lat[n - 1]) k = n - 1;
  else {
    k = 1;
    while (xlat > lat[k] && k < n) k = k + 1;
  }
  for (int m = 0; m < nmon; ++m)
    for (int lv = 0; lv < levsiz; ++lv) {
      auto f = [&](int jj) { return ozmixin[(std::size_t(m) * latsiz + (jj - 1)) * levsiz + lv]; };
      const double a = (f(k + 1) - f(k)) / (lat[k] - lat[k - 1]);
      ozmixm[m * levsiz + lv] = f(k) + a * (xlat - lat[k - 1]);
    }
}

void ozone_time_weights(double julian, int& nm, int& np, double& fact1, double& fact2) {
  static const int date_oz[12] = {16, 45, 75, 105, 136, 166, 197, 228, 258, 289, 319, 350};
  constexpr double daysperyear = 365.;
  double intJULIAN = julian + 1.0;
  int IJUL = int(intJULIAN);
  intJULIAN = intJULIAN - double(IJUL);
  IJUL = IJUL % 365;
  if (IJUL == 0) IJUL = 365;
  intJULIAN = intJULIAN + IJUL;
  int np1 = 1;
  bool finddate = false;
  for (int m = 1; m <= 12; ++m)
    if (date_oz[m - 1] > intJULIAN && !finddate) { np1 = m; finddate = true; }
  const double cdayozp = date_oz[np1 - 1];
  double cdayozm;
  if (np1 > 1) { cdayozm = date_oz[np1 - 2]; np = np1; nm = np - 1; }
  else { cdayozm = date_oz[11]; np = np1; nm = 12; }
  if (np1 == 1) {  // Dec-Jan interpolation (ozncyc)
    const double deltat = cdayozp + daysperyear - cdayozm;
    if (intJULIAN > cdayozp) {
      fact1 = (cdayozp + daysperyear - intJULIAN) / deltat;
      fact2 = (intJULIAN - cdayozm) / deltat;
    } else {
      fact1 = (cdayozp - intJULIAN) / deltat;
      fact2 = (intJULIAN + daysperyear - cdayozm) / deltat;
    }
  } else {
    const double deltat = cdayozp - cdayozm;
    fact1 = (cdayozp - intJULIAN) / deltat;
    fact2 = (intJULIAN - cdayozm) / deltat;
  }
  nm -= 1;  // 0-based months for CamOzone
  np -= 1;
}

}  // namespace cppwrf::phys
