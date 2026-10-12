#include "phys/mp_kessler/kessler.hpp"

namespace cppwrf::phys::kessler {

void run_kessler(
    View3<real> t,
    View3<real> qv,
    View3<real> qc,
    View3<real> qr,
    View3<const real> rho,
    View3<const real> pii,
    View3<const real> dz8w,
    View3<const real> z,
    View2<real> rainnc,
    View2<real> rainncv,
    real dt,
    const Range2D& range,
    int nz,
    real xlv,
    real cp,
    real ep2,
    real svp1,
    real svp2,
    real svp3,
    real svpt0,
    real rhowater) {
  assert(nz <= kMaxLevels);
  KesslerKernel kernel{
      t, qv, qc, qr, rho, pii, dz8w, z, rainnc, rainncv,
      dt, xlv, cp, ep2, svp1, svp2, svp3, svpt0, rhowater, nz};
  column_for(range, kernel);
}

}  // namespace cppwrf::phys::kessler
