#include "phys/bl_myj/myjpbl.hpp"

namespace cppwrf::phys::myjpbl {

void run_myjpbl(
    View2<const real> ht,
    View3<const real> dz,
    View3<const real> pmid,
    View3<const real> pint,
    View3<const real> th,
    View3<const real> t,
    View3<const real> exner,
    View3<const real> qv,
    View3<const real> qc,
    View3<const real> qi,
    View3<const real> u,
    View3<const real> v,
    View3<const real> rho,
    View2<const real> tsk,
    View2<const real> chklowq,
    View2<const real> thz0_in,
    View2<const real> qz0_in,
    View2<const real> uz0,
    View2<const real> vz0,
    View2<const int> lowlyr,
    View2<const real> xland,
    View2<const real> sice,
    View2<const real> snow,
    View2<const real> ustar,
    View2<const real> znt,
    View2<const real> akhs,
    View2<const real> akms,
    View2<const real> elflx,
    View2<real> qsfc,
    View2<real> thz0,
    View2<real> qz0,
    View2<real> ct,
    View3<real> tke_myj,
    View3<real> exch_h,
    View3<real> el_myj,
    View2<real> pblh,
    View2<int> kpbl,
    View2<real> mixht,
    View3<real> rublten,
    View3<real> rvblten,
    View3<real> rthblten,
    View3<real> rqvblten,
    View3<real> rqcblten,
    View3<real> rqiblten,
    const Range2D& range,
    real dt,
    int stepbl,
    int nz,
    bool has_qi) {
  MYJPBLKernel kernel{
      ht, dz, pmid, pint, th, t, exner, qv, qc, qi, u, v, rho,
      tsk, chklowq, thz0_in, qz0_in, uz0, vz0, lowlyr, xland, sice, snow,
      ustar, znt, akhs, akms, elflx,
      qsfc, thz0, qz0, ct, tke_myj, exch_h,
      el_myj, pblh, kpbl, mixht,
      rublten, rvblten, rthblten, rqvblten, rqcblten, rqiblten,
      dt, stepbl, nz, has_qi,
      range.i0, range.i1 - range.i0};

  column_for(range, kernel);
}

void run_myjpblinit(
    View3<real> rublten,
    View3<real> rvblten,
    View3<real> rthblten,
    View3<real> rqvblten,
    View3<real> tke_myj,
    View3<real> exch_h,
    const Range2D& range,
    int nz) {
  MYJPBLInitKernel kernel{
      rublten, rvblten, rthblten, rqvblten, tke_myj, exch_h, nz};

  column_for(range, kernel);
}

}  // namespace cppwrf::phys::myjpbl
