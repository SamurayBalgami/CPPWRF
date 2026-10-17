#include "phys/sf_myjsfc/myjsfc.hpp"
#include <cmath>
#include <vector>

namespace cppwrf::phys::sf_myjsfc {

MYJSFCTableManager::MYJSFCTableManager() {
  std::vector<real> host_psim1(KZTM);
  std::vector<real> host_psih1(KZTM);
  std::vector<real> host_psim2(KZTM);
  std::vector<real> host_psih2(KZTM);

  // Exact statement-by-statement port of MYJSFCINIT table setup (module_sf_myjsfc.F:1095, 1177-1300)
  const real pihf = real(3.1415926) / real(2.);
  const real eps = real(1.e-6);

  real fh01 = real(1.);
  real fh02 = real(1.);

  real ztmin1 = real(-5.0);
  real ztmax1 = real(1.0);
  real ztmin2 = real(-5.0);
  real ztmax2 = real(1.0);

  real zrng1 = ztmax1 - ztmin1;
  real zrng2 = ztmax2 - ztmin2;

  real dzeta1 = zrng1 / real(KZTM - 1);
  real dzeta2 = zrng2 / real(KZTM - 1);

  real zeta1 = ztmin1;
  real zeta2 = ztmin2;

  for (int k = 1; k <= KZTM; ++k) {
    if (zeta1 < real(0.)) {
      real x = std::sqrt(std::sqrt(real(1.) - real(16.) * zeta1));
      host_psim1[k - 1] = real(-2.) * std::log((x + real(1.)) / real(2.)) -
                          std::log((x * x + real(1.)) / real(2.)) +
                          real(2.) * std::atan(x) - pihf;
      host_psih1[k - 1] = real(-2.) * std::log((x * x + real(1.)) / real(2.));
    } else {
      host_psim1[k - 1] = real(0.7) * zeta1 + real(0.75) * zeta1 * (real(6.) - real(0.35) * zeta1) * std::exp(real(-0.35) * zeta1);
      host_psih1[k - 1] = real(0.7) * zeta1 + real(0.75) * zeta1 * (real(6.) - real(0.35) * zeta1) * std::exp(real(-0.35) * zeta1);
    }

    if (zeta2 < real(0.)) {
      real x = std::sqrt(std::sqrt(real(1.) - real(16.) * zeta2));
      host_psim2[k - 1] = real(-2.) * std::log((x + real(1.)) / real(2.)) -
                          std::log((x * x + real(1.)) / real(2.)) +
                          real(2.) * std::atan(x) - pihf;
      host_psih2[k - 1] = real(-2.) * std::log((x * x + real(1.)) / real(2.));
    } else {
      host_psim2[k - 1] = real(0.7) * zeta2 + real(0.75) * zeta2 * (real(6.) - real(0.35) * zeta2) * std::exp(real(-0.35) * zeta2);
      host_psih2[k - 1] = real(0.7) * zeta2 + real(0.75) * zeta2 * (real(6.) - real(0.35) * zeta2) * std::exp(real(-0.35) * zeta2);
    }

    if (k == KZTM) {
      ztmax1 = zeta1;
      ztmax2 = zeta2;
    }

    zeta1 = zeta1 + dzeta1;
    zeta2 = zeta2 + dzeta2;
  }

  ztmax1 = ztmax1 - eps;
  ztmax2 = ztmax2 - eps;

  dev_psim1_ = dev_alloc(KZTM * sizeof(real));
  dev_psih1_ = dev_alloc(KZTM * sizeof(real));
  dev_psim2_ = dev_alloc(KZTM * sizeof(real));
  dev_psih2_ = dev_alloc(KZTM * sizeof(real));

  copy_to_dev(dev_psim1_, host_psim1.data(), KZTM * sizeof(real));
  copy_to_dev(dev_psih1_, host_psih1.data(), KZTM * sizeof(real));
  copy_to_dev(dev_psim2_, host_psim2.data(), KZTM * sizeof(real));
  copy_to_dev(dev_psih2_, host_psih2.data(), KZTM * sizeof(real));

  tables_.psim1 = View1<const real>{static_cast<const real*>(dev_psim1_)};
  tables_.psih1 = View1<const real>{static_cast<const real*>(dev_psih1_)};
  tables_.psim2 = View1<const real>{static_cast<const real*>(dev_psim2_)};
  tables_.psih2 = View1<const real>{static_cast<const real*>(dev_psih2_)};
  tables_.ztmin1 = ztmin1;
  tables_.ztmax1 = ztmax1;
  tables_.dzeta1 = dzeta1;
  tables_.fh01 = fh01;
  tables_.ztmin2 = ztmin2;
  tables_.ztmax2 = ztmax2;
  tables_.dzeta2 = dzeta2;
  tables_.fh02 = fh02;
}

MYJSFCTableManager::~MYJSFCTableManager() {
  if (dev_psim1_) { dev_free(dev_psim1_); dev_psim1_ = nullptr; }
  if (dev_psih1_) { dev_free(dev_psih1_); dev_psih1_ = nullptr; }
  if (dev_psim2_) { dev_free(dev_psim2_); dev_psim2_ = nullptr; }
  if (dev_psih2_) { dev_free(dev_psih2_); dev_psih2_ = nullptr; }
}

void run_myjsfc(
    View2<const real> ht,
    View3<const real> dz,
    View3<const real> pmid,
    View3<const real> pint,
    View3<const real> th,
    View3<const real> t,
    View3<const real> qv,
    View3<const real> qc,
    View3<const real> u,
    View3<const real> v,
    View3<const real> q2,
    View2<const real> tsk,
    View2<const real> xland,
    View2<const int> ivgtyp,
    View2<const real> z0base,
    View2<const real> mavail,
    View2<const int> lowlyr,
    View2<real> qsfc,
    View2<real> thz0,
    View2<real> qz0,
    View2<real> uz0,
    View2<real> vz0,
    View2<real> ustar,
    View2<real> znt,
    View2<real> pblh,
    View2<real> rmol,
    View2<real> akhs,
    View2<real> akms,
    View2<real> rib,
    View2<real> chs,
    View2<real> chs2,
    View2<real> cqs2,
    View2<real> hfx,
    View2<real> qfx,
    View2<real> flx_lh,
    View2<real> flhc,
    View2<real> flqc,
    View2<real> qgh,
    View2<real> cpm,
    View2<real> ct,
    View2<real> u10,
    View2<real> v10,
    View2<real> t02,
    View2<real> th02,
    View2<real> tshltr,
    View2<real> th10,
    View2<real> q02,
    View2<real> qshltr,
    View2<real> q10,
    View2<real> pshltr,
    View2<real> u10e,
    View2<real> v10e,
    const MYJSFCTables& tables,
    const Range2D& range,
    int itimestep,
    int isurban,
    int iz0tlnd,
    real p1000mb,
    int nz) {
  MYJSFCKernel kernel{
      ht, dz, pmid, pint, th, t, qv, qc, u, v, q2,
      tsk, xland, ivgtyp, z0base, mavail, lowlyr,
      qsfc, thz0, qz0, uz0, vz0, ustar, znt, pblh, rmol, akhs, akms,
      rib, chs, chs2, cqs2, hfx, qfx, flx_lh, flhc, flqc, qgh, cpm, ct,
      u10, v10, t02, th02, tshltr, th10, q02, qshltr, q10, pshltr, u10e, v10e,
      tables, itimestep, isurban, iz0tlnd, p1000mb, nz};

  column_for(range, kernel);
}

void run_myjsfcinit(
    View2<int> lowlyr,
    View2<real> ustar,
    const Range2D& range,
    bool restart) {
  MYJSFCInitKernel kernel{lowlyr, ustar, restart};
  column_for(range, kernel);
}

}  // namespace cppwrf::phys::sf_myjsfc
