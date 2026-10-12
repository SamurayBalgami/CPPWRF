#include "kf_tables.hpp"
#include <cmath>

namespace cppwrf::phys::cu_kf {

template <typename T>
void init_kf_tables(KFTables<T>& tables, T svp1, T svp2, T svp3, T svpt0) {
    using std::exp;
    using std::log;
    using std::pow;
    using std::abs;

    const T dth = static_cast<T>(1.0);
    const T tmin = static_cast<T>(150.0);
    const T toler = static_cast<T>(0.001);
    const T plutop = static_cast<T>(5000.0);
    const T pbot = static_cast<T>(110000.0);

    const T aliq = svp1 * static_cast<T>(1000.0);
    const T bliq = svp2;
    const T cliq = svp2 * svpt0;
    const T dliq = svp3;

    tables.rdthk = static_cast<T>(1.0) / dth;
    const T dpr = (pbot - plutop) / static_cast<T>(KFNP - 1);
    tables.rdpr = static_cast<T>(1.0) / dpr;
    tables.plutop = plutop;

    // Starting sat. equiv. theta
    const T temp = tmin;
    T p = plutop - dpr;
    for (int kp = 0; kp < KFNP; ++kp) {
        p += dpr;
        T es = aliq * exp((bliq * temp - cliq) / (temp - dliq));
        T qs = static_cast<T>(0.622) * es / (p - es);
        T pi = pow(static_cast<T>(1.0e5) / p, static_cast<T>(0.2854) * (static_cast<T>(1.0) - static_cast<T>(0.28) * qs));
        tables.the0k[kp] = temp * pi * exp((static_cast<T>(3374.6525) / temp - static_cast<T>(2.5403)) * qs * (static_cast<T>(1.0) + static_cast<T>(0.81) * qs));
    }

    // Compute temperatures for each sat. equiv. potential temp.
    p = plutop - dpr;
    for (int kp = 0; kp < KFNP; ++kp) {
        T thes = tables.the0k[kp] - dth;
        p += dpr;
        for (int it = 0; it < KFNT; ++it) {
            thes += dth;
            T tgues = (it == 0) ? tmin : tables.ttab[kp][it - 1];
            T es = aliq * exp((bliq * tgues - cliq) / (tgues - dliq));
            T qs = static_cast<T>(0.622) * es / (p - es);
            T pi = pow(static_cast<T>(1.0e5) / p, static_cast<T>(0.2854) * (static_cast<T>(1.0) - static_cast<T>(0.28) * qs));
            T thgues = tgues * pi * exp((static_cast<T>(3374.6525) / tgues - static_cast<T>(2.5403)) * qs * (static_cast<T>(1.0) + static_cast<T>(0.81) * qs));
            T f0 = thgues - thes;
            T t1 = tgues - static_cast<T>(0.5) * f0;
            T t0 = tgues;

            for (int itcnt = 1; itcnt <= 11; ++itcnt) {
                es = aliq * exp((bliq * t1 - cliq) / (t1 - dliq));
                qs = static_cast<T>(0.622) * es / (p - es);
                pi = pow(static_cast<T>(1.0e5) / p, static_cast<T>(0.2854) * (static_cast<T>(1.0) - static_cast<T>(0.28) * qs));
                T thtgs = t1 * pi * exp((static_cast<T>(3374.6525) / t1 - static_cast<T>(2.5403)) * qs * (static_cast<T>(1.0) + static_cast<T>(0.81) * qs));
                T f1 = thtgs - thes;
                if (abs(f1) < toler) {
                    break;
                }
                T dt = f1 * (t1 - t0) / (f1 - f0);
                t0 = t1;
                f0 = f1;
                t1 -= dt;
            }
            tables.ttab[kp][it] = t1;
            tables.qstab[kp][it] = qs;
        }
    }

    // Lookup table for tlog(emix/aliq)
    const T astrt = static_cast<T>(1.0e-3);
    const T ainc = static_cast<T>(0.075);
    T a1 = astrt - ainc;
    for (int i = 0; i < KFNALU; ++i) {
        a1 += ainc;
        tables.alu[i] = log(a1);
    }
}

template void init_kf_tables<double>(KFTables<double>&, double, double, double, double);
template void init_kf_tables<float>(KFTables<float>&, float, float, float, float);

} // namespace cppwrf::phys::cu_kf
