// SPDX-License-Identifier: GPL-3.0-or-later
#include "phys/mp_wsm6/wsm6.hpp"
#include "core/limits.hpp"

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
  check_limit("WSM6", "nz", nz, kMaxLevels);
#if defined(CPPWRF_DOUBLE)
  // FP64 builds: in column stages, as src/dyn/physics_coupling.cpp
  Buffer work(Wsm6Work<real>::size(range, nz));
  if (nz <= 48)
    wsm6_stages<real, 48>(params, dt, nz, range, work.get(), th, q, qc, qr, qi, qs, qg, den, pii, p, delz, rain,
                          rainncv, snow, snowncv, graupel, graupelncv, sr);
  else
    wsm6_stages<real, kMaxLevels>(params, dt, nz, range, work.get(), th, q, qc, qr, qi, qs, qg, den, pii, p, delz,
                                  rain, rainncv, snow, snowncv, graupel, graupelncv, sr);
  dev_sync();  // the work buffer is freed on return
#else
  if (nz <= 48) {
    Wsm6Kernel<real, 48> kernel{
        params, dt, nz,
        th, q, qc, qr, qi, qs, qg,
        den, pii, p, delz,
        rain, rainncv, snow, snowncv, graupel, graupelncv, sr};
    column_for(range, kernel);
  } else {
    Wsm6Kernel<real, kMaxLevels> kernel{
        params, dt, nz,
        th, q, qc, qr, qi, qs, qg,
        den, pii, p, delz,
        rain, rainncv, snow, snowncv, graupel, graupelncv, sr};
    column_for(range, kernel);
  }
#endif
}

}  // namespace cppwrf::phys::wsm6
