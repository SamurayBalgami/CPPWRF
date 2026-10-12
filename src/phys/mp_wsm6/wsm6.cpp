#include "phys/mp_wsm6/wsm6.hpp"

namespace cppwrf::phys::wsm6 {

void run_wsm6(
    View3<real> th,
    View3<real> q,
    View3<real> qc,
    View3<real> qr,
    View3<real> qi,
    View3<real> qs,
    View3<real> qg,
    View3<const real> den,
    View3<const real> pii,
    View3<const real> p,
    View3<const real> delz,
    View2<real> rain,
    View2<real> rainncv,
    View2<real> snow,
    View2<real> snowncv,
    View2<real> graupel,
    View2<real> graupelncv,
    View2<real> sr,
    real dt,
    const Range2D& range,
    int nz,
    const Wsm6Params<real>& params) {
  assert(nz <= kMaxLevels);
  Wsm6Kernel<real> kernel{
      params, dt, nz,
      th, q, qc, qr, qi, qs, qg,
      den, pii, p, delz,
      rain, rainncv, snow, snowncv, graupel, graupelncv, sr};
  column_for(range, kernel);
}

}  // namespace cppwrf::phys::wsm6
