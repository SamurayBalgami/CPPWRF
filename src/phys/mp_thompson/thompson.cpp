#include "thompson.hpp"
#include <cassert>

namespace cppwrf::phys::mp_thompson {

template <typename Real>
void run_thompson(
    View3<Real> th,
    View3<Real> qv,
    View3<Real> qc,
    View3<Real> qr,
    View3<Real> qi,
    View3<Real> qs,
    View3<Real> qg,
    View3<Real> ni,
    View3<Real> nr,
    View3<const Real> pii,
    View3<const Real> p,
    View3<const Real> w,
    View3<const Real> dz,
    View2<Real> rainnc,
    View2<Real> rainncv,
    View2<Real> snownc,
    View2<Real> snowncv,
    View2<Real> graupelnc,
    View2<Real> graupelncv,
    View2<Real> sr,
    Real dt,
    const Range2D& range,
    int nz,
    const ThompsonTablesView<Real>& tables) {
    assert(nz <= kMaxLevels);
    ThompsonKernel<Real> kernel{
        tables,
        dt,
        nz,
        th,
        qv,
        qc,
        qr,
        qi,
        qs,
        qg,
        ni,
        nr,
        pii,
        p,
        w,
        dz,
        rainnc,
        rainncv,
        snownc,
        snowncv,
        graupelnc,
        graupelncv,
        sr
    };
    column_for(range, kernel);
}

template void run_thompson<double>(
    View3<double>, View3<double>, View3<double>, View3<double>, View3<double>,
    View3<double>, View3<double>, View3<double>, View3<double>,
    View3<const double>, View3<const double>, View3<const double>, View3<const double>,
    View2<double>, View2<double>, View2<double>, View2<double>,
    View2<double>, View2<double>, View2<double>,
    double, const Range2D&, int, const ThompsonTablesView<double>&);

template void run_thompson<float>(
    View3<float>, View3<float>, View3<float>, View3<float>, View3<float>,
    View3<float>, View3<float>, View3<float>, View3<float>,
    View3<const float>, View3<const float>, View3<const float>, View3<const float>,
    View2<float>, View2<float>, View2<float>, View2<float>,
    View2<float>, View2<float>, View2<float>,
    float, const Range2D&, int, const ThompsonTablesView<float>&);

} // namespace cppwrf::phys::mp_thompson
