#pragma once
#include <cstddef>

namespace cppwrf::phys::cu_kf {

inline constexpr int KFNT = 250;
inline constexpr int KFNP = 220;
inline constexpr int KFNALU = 200;

template <typename T>
struct KFTables {
    T rdpr;
    T rdthk;
    T plutop;
    T the0k[KFNP];
    T alu[KFNALU];
    // ttab[kp][it] matches Fortran TTAB(it, kp) in binary memory layout
    T ttab[KFNP][KFNT];
    T qstab[KFNP][KFNT];
};

template <typename T>
void init_kf_tables(KFTables<T>& tables, T svp1 = static_cast<T>(0.6112),
                    T svp2 = static_cast<T>(17.67),
                    T svp3 = static_cast<T>(29.65),
                    T svpt0 = static_cast<T>(273.15));

} // namespace cppwrf::phys::cu_kf
