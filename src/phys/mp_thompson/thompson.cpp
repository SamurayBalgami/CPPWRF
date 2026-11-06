#include "thompson.hpp"
#include <type_traits>
#include "core/limits.hpp"

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
    const ThompsonTablesView<Real>& tables,
    View3<Real> re_cloud,
    View3<Real> re_ice,
    View3<Real> re_snow,
    const ThompsonAeroFields<Real>* aero,
    View3<Real> refl,
    const radar::RadarTables<Real>* radar) {
    check_limit("Thompson", "nz", nz, kMaxLevels);
    // the column arrays sized to the next tier above nz (on the GPU in the shared work buffer: thompson_reserve)
    auto launch = [&](auto nl) {
        ThompsonKernel<Real, decltype(nl)::value> kernel{
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
            sr,
            re_cloud.p != nullptr,
            re_cloud,
            re_ice,
            re_snow
        };
        if (aero) {
            kernel.aero = true;
            kernel.nc_f = aero->nc;
            kernel.nwfa_f = aero->nwfa;
            kernel.nifa_f = aero->nifa;
            kernel.nbca_f = aero->nbca;
            kernel.nwfa2d = aero->nwfa2d;
            kernel.nifa2d = aero->nifa2d;
        }
        if (refl.p) {
            kernel.refl = refl;
            kernel.rtab = radar;
        }
#if defined(CPPWRF_USE_CUDA)
        // the GPU: level0 (per (i, j, k)), stage1, level (per (i, j, k)), stage3 of each chunk of columns (ThompsonKernel::operator())
        using K = ThompsonKernel<Real, decltype(nl)::value>;
        constexpr int NL = decltype(nl)::value;
        column_work_chunks(range, K::kWorkBytes, [&](const Range2D& c, const ColumnWork& wk) {
            kernel.work = wk;
            const ColumnBlock cb = column_block();
            parallel_for(Range3D{c.i0, c.i1, c.j0, c.j1, 0, nz}, ThompsonLevel0<Real, NL>{kernel, c, cb.x, cb.y});
            column_for(c, ThompsonStage1<Real, NL>{kernel});
            parallel_for(Range3D{c.i0, c.i1, c.j0, c.j1, 0, nz}, ThompsonLevel<Real, NL>{kernel, c, cb.x, cb.y});
            column_for(c, ThompsonStage3<Real, NL>{kernel});
        });
#else
        column_for_work(range, kernel);
#endif
    };
    if (nz <= 48) launch(std::integral_constant<int, 48>{});
    else if (nz <= 64) launch(std::integral_constant<int, 64>{});
    else launch(std::integral_constant<int, kMaxLevels>{});
}

template <typename Real>
void thompson_reserve(const Range2D& range, int nz) {
    const std::size_t b = nz <= 48 ? ThompsonKernel<Real, 48>::kWorkBytes
                          : nz <= 64 ? ThompsonKernel<Real, 64>::kWorkBytes
                                     : ThompsonKernel<Real, kMaxLevels>::kWorkBytes;
    column_work_reserve(range, b);
}

template void thompson_reserve<double>(const Range2D&, int);
template void thompson_reserve<float>(const Range2D&, int);

template void run_thompson<double>(
    View3<double>, View3<double>, View3<double>, View3<double>, View3<double>,
    View3<double>, View3<double>, View3<double>, View3<double>,
    View3<const double>, View3<const double>, View3<const double>, View3<const double>,
    View2<double>, View2<double>, View2<double>, View2<double>,
    View2<double>, View2<double>, View2<double>,
    double, const Range2D&, int, const ThompsonTablesView<double>&, View3<double>, View3<double>, View3<double>,
    const ThompsonAeroFields<double>*, View3<double>, const radar::RadarTables<double>*);

template void run_thompson<float>(
    View3<float>, View3<float>, View3<float>, View3<float>, View3<float>,
    View3<float>, View3<float>, View3<float>, View3<float>,
    View3<const float>, View3<const float>, View3<const float>, View3<const float>,
    View2<float>, View2<float>, View2<float>, View2<float>,
    View2<float>, View2<float>, View2<float>,
    float, const Range2D&, int, const ThompsonTablesView<float>&, View3<float>, View3<float>, View3<float>,
    const ThompsonAeroFields<float>*, View3<float>, const radar::RadarTables<float>*);

} // namespace cppwrf::phys::mp_thompson
