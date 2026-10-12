#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "phys/phys_constants.hpp"

namespace cppwrf::phys::wsm6 {

constexpr int kMaxLevels = 128;

template <typename Real>
inline Real rgmma_calc(Real x) {
    if (x == static_cast<Real>(1.0)) return static_cast<Real>(0.0);
    constexpr double euler = 0.577215664901532;
    double r = static_cast<double>(x) * std::exp(euler * static_cast<double>(x));
    for (int i = 1; i <= 10000; ++i) {
        double y = static_cast<double>(i);
        r = r * (1.000 + static_cast<double>(x) / y) * std::exp(-static_cast<double>(x) / y);
    }
    return static_cast<Real>(1.0 / r);
}

template <typename Real>
struct Wsm6Params {
    Real g{9.81};
    Real cpd{1004.5};
    Real cpv{1846.4};
    Real rd{287.0};
    Real rv{461.6};
    Real t0c{273.15};
    Real ep1{0.6083623693379791};
    Real ep2{0.6217504332755632};
    Real qmin{1.0e-15};
    Real xls{2.85e6};
    Real xlv0{2.50e6};
    Real xlf0{0.35e6};
    Real den0{1.28};
    Real denr{1000.0};
    Real dens{100.0};
    Real cliq{4190.0};
    Real cice{2106.0};
    Real psat{610.78};
    int hail_opt{0};

    Real dtcldcr{120.0};
    Real n0r{8.0e6};
    Real avtr{841.9};
    Real bvtr{0.8};
    Real r0{0.8e-5};
    Real peaut{0.55};
    Real xncr{3.0e8};
    Real xmyu{1.718e-5};
    Real avts{11.72};
    Real bvts{0.41};
    Real lamdarmax{8.0e4};
    Real lamdasmax{1.0e5};
    Real dicon{11.9};
    Real dimax{500.0e-6};
    Real pfrz1{100.0};
    Real pfrz2{0.66};
    Real qcrmin{1.0e-9};
    Real eacrc{1.0};
    Real qs0{6.0e-4};
    Real n0smax{1.0e11};
    Real n0s{2.0e6};
    Real alpha{0.12};

    Real n0g{4.0e6};
    Real deng{500.0};
    Real avtg{330.0};
    Real bvtg{0.8};
    Real lamdagmax{6.0e4};

    Real pi{3.14159265358979323846};
    Real xlv1{0.0};
    Real qc0{0.0};
    Real qck1{0.0};
    Real pidnc{0.0};

    Real bvtr1{0.0}, bvtr2{0.0}, bvtr3{0.0}, bvtr4{0.0}, bvtr6{0.0};
    Real g1pbr{0.0}, g3pbr{0.0}, g4pbr{0.0}, g6pbr{0.0}, g5pbro2{0.0};
    Real pvtr{0.0}, eacrr{1.0}, pacrr{0.0}, precr1{0.0}, precr2{0.0}, roqimax{0.0};

    Real bvts1{0.0}, bvts2{0.0}, bvts3{0.0}, bvts4{0.0};
    Real g1pbs{0.0}, g3pbs{0.0}, g4pbs{0.0}, g5pbso2{0.0};
    Real pvts{0.0}, pacrs{0.0}, precs1{0.0}, precs2{0.0}, pidn0r{0.0}, pidn0s{0.0}, pacrc{0.0};

    Real bvtg1{0.0}, bvtg2{0.0}, bvtg3{0.0}, bvtg4{0.0};
    Real g1pbg{0.0}, g3pbg{0.0}, g4pbg{0.0}, g5pbgo2{0.0};
    Real pvtg{0.0}, pacrg{0.0}, precg1{0.0}, precg2{0.0}, pidn0g{0.0};

    Real rslopermax{0.0}, rslopesmax{0.0}, rslopegmax{0.0};
    Real rsloperbmax{0.0}, rslopesbmax{0.0}, rslopegbmax{0.0};
    Real rsloper2max{0.0}, rslopes2max{0.0}, rslopeg2max{0.0};
    Real rsloper3max{0.0}, rslopes3max{0.0}, rslopeg3max{0.0};

    void init(Real p_den0, Real p_denr, Real p_dens, Real p_cliq, Real p_cpv, int p_hail_opt) {
        den0 = p_den0;
        denr = p_denr;
        dens = p_dens;
        cliq = p_cliq;
        cpv = p_cpv;
        hail_opt = p_hail_opt;

        if (hail_opt == 1) {
            n0g = 4.0e4;
            deng = 700.0;
            avtg = 285.0;
            bvtg = 0.8;
            lamdagmax = 2.0e4;
        } else {
            n0g = 4.0e6;
            deng = 500.0;
            avtg = 330.0;
            bvtg = 0.8;
            lamdagmax = 6.0e4;
        }

        pi = 4.0 * std::atan(1.0);
        xlv1 = cliq - cpv;

        qc0 = (4.0 / 3.0) * pi * denr * (r0 * r0 * r0) * xncr / den0;
        qck1 = 0.104 * 9.8 * peaut / std::pow(xncr * denr, 1.0 / 3.0) / xmyu * std::pow(den0, 4.0 / 3.0);
        pidnc = pi * denr / 6.0;

        bvtr1 = 1.0 + bvtr;
        bvtr2 = 2.5 + 0.5 * bvtr;
        bvtr3 = 3.0 + bvtr;
        bvtr4 = 4.0 + bvtr;
        bvtr6 = 6.0 + bvtr;
        g1pbr = rgmma_calc(bvtr1);
        g3pbr = rgmma_calc(bvtr3);
        g4pbr = rgmma_calc(bvtr4);
        g6pbr = rgmma_calc(bvtr6);
        g5pbro2 = rgmma_calc(bvtr2);
        pvtr = avtr * g4pbr / 6.0;
        eacrr = 1.0;
        pacrr = pi * n0r * avtr * g3pbr * 0.25 * eacrr;
        precr1 = 2.0 * pi * n0r * 0.78;
        precr2 = 2.0 * pi * n0r * 0.31 * std::sqrt(avtr) * g5pbro2;
        roqimax = 2.08e22 * std::pow(dimax, 8.0);

        bvts1 = 1.0 + bvts;
        bvts2 = 2.5 + 0.5 * bvts;
        bvts3 = 3.0 + bvts;
        bvts4 = 4.0 + bvts;
        g1pbs = rgmma_calc(bvts1);
        g3pbs = rgmma_calc(bvts3);
        g4pbs = rgmma_calc(bvts4);
        g5pbso2 = rgmma_calc(bvts2);
        pvts = avts * g4pbs / 6.0;
        pacrs = pi * n0s * avts * g3pbs * 0.25;
        precs1 = 4.0 * n0s * 0.65;
        precs2 = 4.0 * n0s * 0.44 * std::sqrt(avts) * g5pbso2;
        pidn0r = pi * denr * n0r;
        pidn0s = pi * dens * n0s;
        pacrc = pi * n0s * avts * g3pbs * 0.25 * eacrc;

        bvtg1 = 1.0 + bvtg;
        bvtg2 = 2.5 + 0.5 * bvtg;
        bvtg3 = 3.0 + bvtg;
        bvtg4 = 4.0 + bvtg;
        g1pbg = rgmma_calc(bvtg1);
        g3pbg = rgmma_calc(bvtg3);
        g4pbg = rgmma_calc(bvtg4);
        g5pbgo2 = rgmma_calc(bvtg2);
        pvtg = avtg * g4pbg / 6.0;
        pacrg = pi * n0g * avtg * g3pbg * 0.25;
        precg1 = 2.0 * pi * n0g * 0.78;
        precg2 = 2.0 * pi * n0g * 0.31 * std::sqrt(avtg) * g5pbgo2;
        pidn0g = pi * deng * n0g;

        rslopermax = 1.0 / lamdarmax;
        rslopesmax = 1.0 / lamdasmax;
        rslopegmax = 1.0 / lamdagmax;
        rsloperbmax = std::pow(rslopermax, bvtr);
        rslopesbmax = std::pow(rslopesmax, bvts);
        rslopegbmax = std::pow(rslopegmax, bvtg);
        rsloper2max = rslopermax * rslopermax;
        rslopes2max = rslopesmax * rslopesmax;
        rslopeg2max = rslopegmax * rslopegmax;
        rsloper3max = rsloper2max * rslopermax;
        rslopes3max = rslopes2max * rslopesmax;
        rslopeg3max = rslopeg2max * rslopegmax;
    }

    CPPWRF_HD inline Real cpmcal(Real x) const {
        return cpd * (1.0 - std::max(x, qmin)) + std::max(x, qmin) * cpv;
    }

    CPPWRF_HD inline Real xlcal(Real x) const {
        return xlv0 - xlv1 * (x - t0c);
    }

    CPPWRF_HD inline Real diffus(Real x, Real y) const {
        return 8.794e-5 * std::exp(std::log(x) * 1.81) / y;
    }

    CPPWRF_HD inline Real viscos(Real x, Real y) const {
        return 1.496e-6 * (x * std::sqrt(x)) / (x + 120.0) / y;
    }

    CPPWRF_HD inline Real xka(Real x, Real y) const {
        return 1.414e3 * viscos(x, y) * y;
    }

    CPPWRF_HD inline Real diffac(Real a, Real b, Real c, Real d, Real e) const {
        return d * a * a / (xka(c, d) * rv * c * c) + 1.0 / (e * diffus(c, b));
    }

    CPPWRF_HD inline Real venfac(Real a, Real b, Real c) const {
        Real v = viscos(b, c);
        return std::exp(std::log(v / diffus(b, a)) * 0.3333333) / std::sqrt(v) * std::sqrt(std::sqrt(den0 / c));
    }

    CPPWRF_HD inline Real conden(Real a, Real b, Real c, Real d, Real e) const {
        return (std::max(b, qmin) - c) / (1.0 + d * d / (rv * e) * c / (a * a));
    }

    CPPWRF_HD inline void slope_rain(int km, const Real* qrs, const Real* den, const Real* denfac, const Real* t, Real* vt) const {
        for (int k = 0; k < km; ++k) {
            Real rslope, rslopeb;
            if (qrs[k] <= qcrmin) {
                rslope = rslopermax;
                rslopeb = rsloperbmax;
            } else {
                Real lamdar = std::sqrt(std::sqrt(pidn0r / (qrs[k] * den[k])));
                rslope = 1.0 / lamdar;
                rslopeb = std::pow(rslope, bvtr);
            }
            vt[k] = pvtr * rslopeb * denfac[k];
            if (qrs[k] <= 0.0) vt[k] = 0.0;
        }
    }

    CPPWRF_HD inline void slope_snow(int km, const Real* qrs, const Real* den, const Real* denfac, const Real* t, Real* vt) const {
        for (int k = 0; k < km; ++k) {
            Real supcol = t0c - t[k];
            Real n0sfac = std::max(std::min(std::exp(alpha * supcol), n0smax / n0s), static_cast<Real>(1.0));
            Real rslope, rslopeb;
            if (qrs[k] <= qcrmin) {
                rslope = rslopesmax;
                rslopeb = rslopesbmax;
            } else {
                Real lamdas = std::sqrt(std::sqrt(pidn0s * n0sfac / (qrs[k] * den[k])));
                rslope = 1.0 / lamdas;
                rslopeb = std::pow(rslope, bvts);
            }
            vt[k] = pvts * rslopeb * denfac[k];
            if (qrs[k] <= 0.0) vt[k] = 0.0;
        }
    }

    CPPWRF_HD inline void slope_graup(int km, const Real* qrs, const Real* den, const Real* denfac, const Real* t, Real* vt) const {
        for (int k = 0; k < km; ++k) {
            Real rslope, rslopeb;
            if (qrs[k] <= qcrmin) {
                rslope = rslopegmax;
                rslopeb = rslopegbmax;
            } else {
                Real lamdag = std::sqrt(std::sqrt(pidn0g / (qrs[k] * den[k])));
                rslope = 1.0 / lamdag;
                rslopeb = std::pow(rslope, bvtg);
            }
            vt[k] = pvtg * rslopeb * denfac[k];
            if (qrs[k] <= 0.0) vt[k] = 0.0;
        }
    }

    CPPWRF_HD inline void slope_wsm6(int km, const Real* qr, const Real* qs, const Real* qg,
                                    const Real* den, const Real* denfac, const Real* t,
                                    Real rslope[][3], Real rslopeb[][3], Real rslope2[][3],
                                    Real rslope3[][3], Real vt[][3]) const {
        for (int k = 0; k < km; ++k) {
            Real supcol = t0c - t[k];
            Real n0sfac = std::max(std::min(std::exp(alpha * supcol), n0smax / n0s), static_cast<Real>(1.0));

            // Rain
            if (qr[k] <= qcrmin) {
                rslope[k][0] = rslopermax;
                rslopeb[k][0] = rsloperbmax;
                rslope2[k][0] = rsloper2max;
                rslope3[k][0] = rsloper3max;
            } else {
                Real lamdar = std::sqrt(std::sqrt(pidn0r / (qr[k] * den[k])));
                rslope[k][0] = 1.0 / lamdar;
                rslopeb[k][0] = std::pow(rslope[k][0], bvtr);
                rslope2[k][0] = rslope[k][0] * rslope[k][0];
                rslope3[k][0] = rslope2[k][0] * rslope[k][0];
            }

            // Snow
            if (qs[k] <= qcrmin) {
                rslope[k][1] = rslopesmax;
                rslopeb[k][1] = rslopesbmax;
                rslope2[k][1] = rslopes2max;
                rslope3[k][1] = rslopes3max;
            } else {
                Real lamdas = std::sqrt(std::sqrt(pidn0s * n0sfac / (qs[k] * den[k])));
                rslope[k][1] = 1.0 / lamdas;
                rslopeb[k][1] = std::pow(rslope[k][1], bvts);
                rslope2[k][1] = rslope[k][1] * rslope[k][1];
                rslope3[k][1] = rslope2[k][1] * rslope[k][1];
            }

            // Graupel
            if (qg[k] <= qcrmin) {
                rslope[k][2] = rslopegmax;
                rslopeb[k][2] = rslopegbmax;
                rslope2[k][2] = rslopeg2max;
                rslope3[k][2] = rslopeg3max;
            } else {
                Real lamdag = std::sqrt(std::sqrt(pidn0g / (qg[k] * den[k])));
                rslope[k][2] = 1.0 / lamdag;
                rslopeb[k][2] = std::pow(rslope[k][2], bvtg);
                rslope2[k][2] = rslope[k][2] * rslope[k][2];
                rslope3[k][2] = rslope2[k][2] * rslope[k][2];
            }

            vt[k][0] = pvtr * rslopeb[k][0] * denfac[k];
            vt[k][1] = pvts * rslopeb[k][1] * denfac[k];
            vt[k][2] = pvtg * rslopeb[k][2] * denfac[k];

            if (qr[k] <= 0.0) vt[k][0] = 0.0;
            if (qs[k] <= 0.0) vt[k][1] = 0.0;
            if (qg[k] <= 0.0) vt[k][2] = 0.0;
        }
    }

    CPPWRF_HD inline void nislfv_rain_plm(int km, const Real* den, const Real* denfac,
                                        const Real* t, const Real* dz,
                                        Real* ww, Real* rql, Real& precip,
                                        Real dt, int id, int iter) const {
        precip = 0.0;
        Real allold = 0.0;
        for (int k = 0; k < km; ++k) {
            allold += rql[k];
        }
        if (allold <= 0.0) {
            return;
        }

        Real zi[kMaxLevels + 1];
        Real za[kMaxLevels + 1];
        Real wi[kMaxLevels + 1];
        Real dza[kMaxLevels + 1];
        Real qa[kMaxLevels + 1];
        Real qmi[kMaxLevels + 1];
        Real qpi[kMaxLevels + 1];
        Real wd[kMaxLevels];
        Real wa[kMaxLevels];
        Real was[kMaxLevels];
        Real qn[kMaxLevels];
        Real qr_work[kMaxLevels];

        zi[0] = 0.0;
        for (int k = 0; k < km; ++k) {
            zi[k + 1] = zi[k] + dz[k];
        }

        for (int k = 0; k < km; ++k) {
            wd[k] = ww[k];
        }

        for (int n = 1; ; ++n) {
            // 2nd order wi
            wi[0] = ww[0];
            wi[km] = ww[km - 1];
            for (int k = 1; k < km; ++k) {
                wi[k] = (ww[k] * dz[k - 1] + ww[k - 1] * dz[k]) / (dz[k - 1] + dz[k]);
            }
            // 3rd order wi
            constexpr Real fa1 = 9.0 / 16.0;
            constexpr Real fa2 = 1.0 / 16.0;
            wi[0] = ww[0];
            wi[1] = 0.5 * (ww[1] + ww[0]);
            for (int k = 2; k < km - 1; ++k) {
                wi[k] = fa1 * (ww[k] + ww[k - 1]) - fa2 * (ww[k + 1] + ww[k - 2]);
            }
            if (km >= 2) {
                wi[km - 1] = 0.5 * (ww[km - 1] + ww[km - 2]);
            }
            wi[km] = ww[km - 1];

            // terminate of top of raingroup
            for (int k = 1; k < km; ++k) {
                if (ww[k] == 0.0) {
                    wi[k] = ww[k - 1];
                }
            }

            // diffusivity of wi
            constexpr Real con1 = 0.05;
            for (int k = km - 1; k >= 0; --k) {
                Real decfl = (wi[k + 1] - wi[k]) * dt / dz[k];
                if (decfl > con1) {
                    wi[k] = wi[k + 1] - con1 * dz[k] / dt;
                }
            }

            // arrival point
            for (int k = 0; k <= km; ++k) {
                za[k] = zi[k] - wi[k] * dt;
            }

            for (int k = 0; k < km; ++k) {
                dza[k] = za[k + 1] - za[k];
            }
            dza[km] = zi[km] - za[km];

            for (int k = 0; k < km; ++k) {
                qa[k] = rql[k] * dz[k] / dza[k];
                qr_work[k] = qa[k] / den[k];
            }
            qa[km] = 0.0;

            if (n <= iter) {
                slope_rain(km, qr_work, den, denfac, t, wa);
                if (n >= 2) {
                    for (int k = 0; k < km; ++k) {
                        wa[k] = 0.5 * (wa[k] + was[k]);
                    }
                }
                for (int k = 0; k < km; ++k) {
                    ww[k] = 0.5 * (wd[k] + wa[k]);
                    was[k] = wa[k];
                }
            } else {
                break;
            }
        }

        // estimate values at arrival cell interface with monotone
        for (int k = 1; k < km; ++k) {
            Real dip = (qa[k + 1] - qa[k]) / (dza[k + 1] + dza[k]);
            Real dim = (qa[k] - qa[k - 1]) / (dza[k - 1] + dza[k]);
            if (dip * dim <= 0.0) {
                qmi[k] = qa[k];
                qpi[k] = qa[k];
            } else {
                qpi[k] = qa[k] + 0.5 * (dip + dim) * dza[k];
                qmi[k] = 2.0 * qa[k] - qpi[k];
                if (qpi[k] < 0.0 || qmi[k] < 0.0) {
                    qpi[k] = qa[k];
                    qmi[k] = qa[k];
                }
            }
        }
        qpi[0] = qa[0];
        qmi[0] = qa[0];
        qmi[km] = qa[km];
        qpi[km] = qa[km];

        // interpolation to regular point
        for (int k = 0; k < km; ++k) {
            qn[k] = 0.0;
        }
        int kb = 0;
        int kt = 0;
        for (int k = 0; k < km; ++k) {
            kb = (kb - 1 > 0) ? kb - 1 : 0;
            kt = (kt - 1 > 0) ? kt - 1 : 0;
            if (zi[k] >= za[km]) {
                break;
            }
            for (int kk = kb; kk < km; ++kk) {
                if (zi[k] <= za[kk + 1]) {
                    kb = kk;
                    break;
                }
            }
            for (int kk = kt; kk < km; ++kk) {
                if (zi[k + 1] <= za[kk]) {
                    kt = kk;
                    break;
                }
            }
            kt = kt - 1;

            if (kt == kb) {
                Real tl = (zi[k] - za[kb]) / dza[kb];
                Real th = (zi[k + 1] - za[kb]) / dza[kb];
                Real tl2 = tl * tl;
                Real th2 = th * th;
                Real qqd = 0.5 * (qpi[kb] - qmi[kb]);
                Real qqh = qqd * th2 + qmi[kb] * th;
                Real qql = qqd * tl2 + qmi[kb] * tl;
                qn[k] = (qqh - qql) / (th - tl);
            } else if (kt > kb) {
                Real tl = (zi[k] - za[kb]) / dza[kb];
                Real tl2 = tl * tl;
                Real qqd = 0.5 * (qpi[kb] - qmi[kb]);
                Real qql = qqd * tl2 + qmi[kb] * tl;
                Real dql = qa[kb] - qql;
                Real zsum = (1.0 - tl) * dza[kb];
                Real qsum = dql * dza[kb];
                if (kt - kb > 1) {
                    for (int m = kb + 1; m <= kt - 1; ++m) {
                        zsum += dza[m];
                        qsum += qa[m] * dza[m];
                    }
                }
                Real th = (zi[k + 1] - za[kt]) / dza[kt];
                Real th2 = th * th;
                qqd = 0.5 * (qpi[kt] - qmi[kt]);
                Real dqh = qqd * th2 + qmi[kt] * th;
                zsum += th * dza[kt];
                qsum += dqh * dza[kt];
                qn[k] = qsum / zsum;
            }
        }

        // rain out
        for (int k = 0; k < km; ++k) {
            if (za[k] < 0.0 && za[k + 1] < 0.0) {
                precip += qa[k] * dza[k];
            } else if (za[k] < 0.0 && za[k + 1] >= 0.0) {
                precip += qa[k] * (0.0 - za[k]);
                break;
            } else {
                break;
            }
        }

        for (int k = 0; k < km; ++k) {
            rql[k] = qn[k];
        }
    }

    CPPWRF_HD inline void nislfv_rain_plm6(int km, const Real* den, const Real* denfac,
                                         const Real* t, const Real* dz,
                                         Real* ww, Real* rql, Real* rql2,
                                         Real& precip1, Real& precip2,
                                         Real dt, int id, int iter) const {
        precip1 = 0.0;
        precip2 = 0.0;
        Real allold = 0.0;
        for (int k = 0; k < km; ++k) {
            allold += rql[k] + rql2[k];
        }
        if (allold <= 0.0) {
            return;
        }

        Real zi[kMaxLevels + 1];
        Real za[kMaxLevels + 1];
        Real wi[kMaxLevels + 1];
        Real dza[kMaxLevels + 1];
        Real qa[kMaxLevels + 1];
        Real qa2[kMaxLevels + 1];
        Real qmi[kMaxLevels + 1];
        Real qpi[kMaxLevels + 1];
        Real wd[kMaxLevels];
        Real wa[kMaxLevels];
        Real wa2[kMaxLevels];
        Real was[kMaxLevels];
        Real qn[kMaxLevels];
        Real qr_work[kMaxLevels];
        Real qr2_work[kMaxLevels];

        zi[0] = 0.0;
        for (int k = 0; k < km; ++k) {
            zi[k + 1] = zi[k] + dz[k];
        }

        for (int k = 0; k < km; ++k) {
            wd[k] = ww[k];
        }

        for (int n = 1; ; ++n) {
            // 2nd order wi
            wi[0] = ww[0];
            wi[km] = ww[km - 1];
            for (int k = 1; k < km; ++k) {
                wi[k] = (ww[k] * dz[k - 1] + ww[k - 1] * dz[k]) / (dz[k - 1] + dz[k]);
            }
            // 3rd order wi
            constexpr Real fa1 = 9.0 / 16.0;
            constexpr Real fa2 = 1.0 / 16.0;
            wi[0] = ww[0];
            wi[1] = 0.5 * (ww[1] + ww[0]);
            for (int k = 2; k < km - 1; ++k) {
                wi[k] = fa1 * (ww[k] + ww[k - 1]) - fa2 * (ww[k + 1] + ww[k - 2]);
            }
            if (km >= 2) {
                wi[km - 1] = 0.5 * (ww[km - 1] + ww[km - 2]);
            }
            wi[km] = ww[km - 1];

            // terminate of top of raingroup
            for (int k = 1; k < km; ++k) {
                if (ww[k] == 0.0) {
                    wi[k] = ww[k - 1];
                }
            }

            // diffusivity of wi
            constexpr Real con1 = 0.05;
            for (int k = km - 1; k >= 0; --k) {
                Real decfl = (wi[k + 1] - wi[k]) * dt / dz[k];
                if (decfl > con1) {
                    wi[k] = wi[k + 1] - con1 * dz[k] / dt;
                }
            }

            // arrival point
            for (int k = 0; k <= km; ++k) {
                za[k] = zi[k] - wi[k] * dt;
            }

            for (int k = 0; k < km; ++k) {
                dza[k] = za[k + 1] - za[k];
            }
            dza[km] = zi[km] - za[km];

            for (int k = 0; k < km; ++k) {
                qa[k] = rql[k] * dz[k] / dza[k];
                qa2[k] = rql2[k] * dz[k] / dza[k];
                qr_work[k] = qa[k] / den[k];
                qr2_work[k] = qa2[k] / den[k];
            }
            qa[km] = 0.0;
            qa2[km] = 0.0;

            if (n <= iter) {
                slope_snow(km, qr_work, den, denfac, t, wa);
                slope_graup(km, qr2_work, den, denfac, t, wa2);
                for (int k = 0; k < km; ++k) {
                    Real tmp_k = std::max(qr_work[k] + qr2_work[k], static_cast<Real>(1.0e-15));
                    if (tmp_k > static_cast<Real>(1.0e-15)) {
                        wa[k] = (wa[k] * qr_work[k] + wa2[k] * qr2_work[k]) / tmp_k;
                    } else {
                        wa[k] = 0.0;
                    }
                }
                if (n >= 2) {
                    for (int k = 0; k < km; ++k) {
                        wa[k] = 0.5 * (wa[k] + was[k]);
                    }
                }
                for (int k = 0; k < km; ++k) {
                    ww[k] = 0.5 * (wd[k] + wa[k]);
                    was[k] = wa[k];
                }
            } else {
                break;
            }
        }

        for (int ist = 1; ist <= 2; ++ist) {
            const Real* cur_qa = (ist == 1) ? qa : qa2;
            Real cur_precip = 0.0;

            for (int k = 1; k < km; ++k) {
                Real dip = (cur_qa[k + 1] - cur_qa[k]) / (dza[k + 1] + dza[k]);
                Real dim = (cur_qa[k] - cur_qa[k - 1]) / (dza[k - 1] + dza[k]);
                if (dip * dim <= 0.0) {
                    qmi[k] = cur_qa[k];
                    qpi[k] = cur_qa[k];
                } else {
                    qpi[k] = cur_qa[k] + 0.5 * (dip + dim) * dza[k];
                    qmi[k] = 2.0 * cur_qa[k] - qpi[k];
                    if (qpi[k] < 0.0 || qmi[k] < 0.0) {
                        qpi[k] = cur_qa[k];
                        qmi[k] = cur_qa[k];
                    }
                }
            }
            qpi[0] = cur_qa[0];
            qmi[0] = cur_qa[0];
            qmi[km] = cur_qa[km];
            qpi[km] = cur_qa[km];

            for (int k = 0; k < km; ++k) {
                qn[k] = 0.0;
            }
            int kb = 0;
            int kt = 0;
            for (int k = 0; k < km; ++k) {
                kb = (kb - 1 > 0) ? kb - 1 : 0;
                kt = (kt - 1 > 0) ? kt - 1 : 0;
                if (zi[k] >= za[km]) {
                    break;
                }
                for (int kk = kb; kk < km; ++kk) {
                    if (zi[k] <= za[kk + 1]) {
                        kb = kk;
                        break;
                    }
                }
                for (int kk = kt; kk < km; ++kk) {
                    if (zi[k + 1] <= za[kk]) {
                        kt = kk;
                        break;
                    }
                }
                kt = kt - 1;

                if (kt == kb) {
                    Real tl = (zi[k] - za[kb]) / dza[kb];
                    Real th = (zi[k + 1] - za[kb]) / dza[kb];
                    Real tl2 = tl * tl;
                    Real th2 = th * th;
                    Real qqd = 0.5 * (qpi[kb] - qmi[kb]);
                    Real qqh = qqd * th2 + qmi[kb] * th;
                    Real qql = qqd * tl2 + qmi[kb] * tl;
                    qn[k] = (qqh - qql) / (th - tl);
                } else if (kt > kb) {
                    Real tl = (zi[k] - za[kb]) / dza[kb];
                    Real tl2 = tl * tl;
                    Real qqd = 0.5 * (qpi[kb] - qmi[kb]);
                    Real qql = qqd * tl2 + qmi[kb] * tl;
                    Real dql = cur_qa[kb] - qql;
                    Real zsum = (1.0 - tl) * dza[kb];
                    Real qsum = dql * dza[kb];
                    if (kt - kb > 1) {
                        for (int m = kb + 1; m <= kt - 1; ++m) {
                            zsum += dza[m];
                            qsum += cur_qa[m] * dza[m];
                        }
                    }
                    Real th = (zi[k + 1] - za[kt]) / dza[kt];
                    Real th2 = th * th;
                    qqd = 0.5 * (qpi[kt] - qmi[kt]);
                    Real dqh = qqd * th2 + qmi[kt] * th;
                    zsum += th * dza[kt];
                    qsum += dqh * dza[kt];
                    qn[k] = qsum / zsum;
                }
            }

            // rain out
            for (int k = 0; k < km; ++k) {
                if (za[k] < 0.0 && za[k + 1] < 0.0) {
                    cur_precip += cur_qa[k] * dza[k];
                } else if (za[k] < 0.0 && za[k + 1] >= 0.0) {
                    cur_precip += cur_qa[k] * (0.0 - za[k]);
                    break;
                } else {
                    break;
                }
            }

            if (ist == 1) {
                for (int k = 0; k < km; ++k) rql[k] = qn[k];
                precip1 = cur_precip;
            } else {
                for (int k = 0; k < km; ++k) rql2[k] = qn[k];
                precip2 = cur_precip;
            }
        }
    }

    CPPWRF_HD inline void operator()(int nz, Real delt,
                                    Real* th, Real* q, Real* qc, Real* qr,
                                    Real* qi, Real* qs, Real* qg,
                                    const Real* den, const Real* pii, const Real* p,
                                    const Real* delz,
                                    Real& rain, Real& rainncv,
                                    Real& snow, Real& snowncv,
                                    Real& graupel, Real& graupelncv,
                                    Real& sr) const {
        for (int k = 0; k < nz; ++k) {
            qc[k] = std::max(qc[k], static_cast<Real>(0.0));
            qr[k] = std::max(qr[k], static_cast<Real>(0.0));
            qi[k] = std::max(qi[k], static_cast<Real>(0.0));
            qs[k] = std::max(qs[k], static_cast<Real>(0.0));
            qg[k] = std::max(qg[k], static_cast<Real>(0.0));
        }

        Real t[kMaxLevels];
        Real cpm[kMaxLevels];
        Real xl[kMaxLevels];
        for (int k = 0; k < nz; ++k) {
            t[k] = th[k] * pii[k];
            cpm[k] = cpmcal(q[k]);
            xl[k] = xlcal(t[k]);
        }

        rainncv = 0.0;
        snowncv = 0.0;
        graupelncv = 0.0;
        sr = 0.0;
        Real tstepsnow = 0.0;
        Real tstepgraup = 0.0;

        int loops = std::max(1, static_cast<int>(std::round(delt / dtcldcr)));
        Real dtcld = delt / static_cast<Real>(loops);
        if (delt <= dtcldcr) dtcld = delt;

        for (int loop = 0; loop < loops; ++loop) {
            Real denfac[kMaxLevels];
            for (int k = 0; k < nz; ++k) {
                denfac[k] = std::sqrt((static_cast<Real>(1.0) / den[k]) * den0);
            }

            Real qsat[kMaxLevels][2];
            Real rh[kMaxLevels][2];
            Real hsub = xls;
            Real hvap = xlv0;
            Real cvap = cpv;
            Real ttp = t0c + 0.01;
            Real dldt = cvap - cliq;
            Real xa = -dldt / rv;
            Real xb = xa + hvap / (rv * ttp);
            Real dldti = cvap - cice;
            Real xai = -dldti / rv;
            Real xbi = xai + hsub / (rv * ttp);

            for (int k = 0; k < nz; ++k) {
                Real tr = ttp / t[k];
                Real qs1 = psat * std::exp(std::log(tr) * xa) * std::exp(xb * (1.0 - tr));
                qs1 = std::min(qs1, static_cast<Real>(0.99) * p[k]);
                qs1 = ep2 * qs1 / (p[k] - qs1);
                qsat[k][0] = std::max(qs1, qmin);
                rh[k][0] = std::max(q[k] / qsat[k][0], qmin);

                Real qs2;
                if (t[k] < ttp) {
                    qs2 = psat * std::exp(std::log(tr) * xai) * std::exp(xbi * (1.0 - tr));
                } else {
                    qs2 = psat * std::exp(std::log(tr) * xa) * std::exp(xb * (1.0 - tr));
                }
                qs2 = std::min(qs2, static_cast<Real>(0.99) * p[k]);
                qs2 = ep2 * qs2 / (p[k] - qs2);
                qsat[k][1] = std::max(qs2, qmin);
                rh[k][1] = std::max(q[k] / qsat[k][1], qmin);
            }

            Real xni[kMaxLevels];
            for (int k = 0; k < nz; ++k) {
                Real temp = den[k] * std::max(qi[k], qmin);
                temp = std::sqrt(std::sqrt(temp * temp * temp));
                xni[k] = std::min(std::max(static_cast<Real>(5.38e7) * temp, static_cast<Real>(1.0e3)), static_cast<Real>(1.0e6));
            }

            Real rslope[kMaxLevels][3];
            Real rslopeb[kMaxLevels][3];
            Real rslope2[kMaxLevels][3];
            Real rslope3[kMaxLevels][3];
            Real vt[kMaxLevels][3];
            slope_wsm6(nz, qr, qs, qg, den, denfac, t, rslope, rslopeb, rslope2, rslope3, vt);

            Real workr[kMaxLevels];
            Real worka[kMaxLevels];
            Real denqrs1[kMaxLevels];
            Real denqrs2[kMaxLevels];
            Real denqrs3[kMaxLevels];
            for (int k = nz - 1; k >= 0; --k) {
                workr[k] = vt[k][0];
                Real qsum_val = std::max(qs[k] + qg[k], static_cast<Real>(1.0e-15));
                if (qsum_val > static_cast<Real>(1.0e-15)) {
                    worka[k] = (vt[k][1] * qs[k] + vt[k][2] * qg[k]) / qsum_val;
                } else {
                    worka[k] = 0.0;
                }
                denqrs1[k] = den[k] * qr[k];
                denqrs2[k] = den[k] * qs[k];
                denqrs3[k] = den[k] * qg[k];
                if (qr[k] <= 0.0) workr[k] = 0.0;
            }

            Real delqrs1 = 0.0, delqrs2 = 0.0, delqrs3 = 0.0;
            nislfv_rain_plm(nz, den, denfac, t, delz, workr, denqrs1, delqrs1, dtcld, 1, 1);
            nislfv_rain_plm6(nz, den, denfac, t, delz, worka, denqrs2, denqrs3, delqrs2, delqrs3, dtcld, 1, 1);

            Real fall[kMaxLevels][3];
            for (int k = 0; k < nz; ++k) {
                qr[k] = std::max(denqrs1[k] / den[k], static_cast<Real>(0.0));
                qs[k] = std::max(denqrs2[k] / den[k], static_cast<Real>(0.0));
                qg[k] = std::max(denqrs3[k] / den[k], static_cast<Real>(0.0));
                fall[k][0] = denqrs1[k] * workr[k] / delz[k];
                fall[k][1] = denqrs2[k] * worka[k] / delz[k];
                fall[k][2] = denqrs3[k] * worka[k] / delz[k];
            }
            fall[0][0] = delqrs1 / delz[0] / dtcld;
            fall[0][1] = delqrs2 / delz[0] / dtcld;
            fall[0][2] = delqrs3 / delz[0] / dtcld;

            slope_wsm6(nz, qr, qs, qg, den, denfac, t, rslope, rslopeb, rslope2, rslope3, vt);

            for (int k = nz - 1; k >= 0; --k) {
                Real supcol = t0c - t[k];
                Real n0sfac_val = std::max(std::min(std::exp(alpha * supcol), n0smax / n0s), static_cast<Real>(1.0));
                if (t[k] > t0c) {
                    Real xlf = xlf0;
                    Real work2_val = venfac(p[k], t[k], den[k]);
                    if (qs[k] > 0.0) {
                        Real coeres = rslope2[k][1] * std::sqrt(rslope[k][1] * rslopeb[k][1]);
                        Real psmlt_val = xka(t[k], den[k]) / xlf * (t0c - t[k]) * (pi / 2.0)
                            * n0sfac_val * (precs1 * rslope2[k][1] + precs2 * work2_val * coeres) / den[k];
                        psmlt_val = std::min(std::max(psmlt_val * dtcld, -qs[k]), static_cast<Real>(0.0));
                        qs[k] += psmlt_val;
                        qr[k] -= psmlt_val;
                        t[k] += (xlf / cpm[k]) * psmlt_val;
                    }
                    if (qg[k] > 0.0) {
                        Real coeres = rslope2[k][2] * std::sqrt(rslope[k][2] * rslopeb[k][2]);
                        Real pgmlt_val = xka(t[k], den[k]) / xlf * (t0c - t[k])
                            * (precg1 * rslope2[k][2] + precg2 * work2_val * coeres) / den[k];
                        pgmlt_val = std::min(std::max(pgmlt_val * dtcld, -qg[k]), static_cast<Real>(0.0));
                        qg[k] += pgmlt_val;
                        qr[k] -= pgmlt_val;
                        t[k] += (xlf / cpm[k]) * pgmlt_val;
                    }
                }
            }

            Real work1c[kMaxLevels];
            Real denqci[kMaxLevels];
            for (int k = nz - 1; k >= 0; --k) {
                if (qi[k] <= 0.0) {
                    work1c[k] = 0.0;
                } else {
                    Real xmi = den[k] * qi[k] / xni[k];
                    Real diameter = std::max(std::min(dicon * std::sqrt(xmi), dimax), static_cast<Real>(1.0e-25));
                    work1c[k] = 1.49e4 * std::exp(std::log(diameter) * 1.31);
                }
                denqci[k] = den[k] * qi[k];
            }
            Real delqi = 0.0;
            nislfv_rain_plm(nz, den, denfac, t, delz, work1c, denqci, delqi, dtcld, 1, 0);
            for (int k = 0; k < nz; ++k) {
                qi[k] = std::max(denqci[k] / den[k], static_cast<Real>(0.0));
            }
            Real fallc = delqi / delz[0] / dtcld;

            Real fallsum = fall[0][0] + fall[0][1] + fall[0][2] + fallc;
            Real fallsum_qsi = fall[0][1] + fallc;
            Real fallsum_qg = fall[0][2];
            if (fallsum > 0.0) {
                Real d_rain = fallsum * delz[0] / denr * dtcld * 1000.0;
                rainncv += d_rain;
                rain += d_rain;
            }
            if (fallsum_qsi > 0.0) {
                Real d_snow = fallsum_qsi * delz[0] / denr * dtcld * 1000.0;
                tstepsnow += d_snow;
                snowncv += d_snow;
                snow += d_snow;
            }
            if (fallsum_qg > 0.0) {
                Real d_graup = fallsum_qg * delz[0] / denr * dtcld * 1000.0;
                tstepgraup += d_graup;
                graupelncv += d_graup;
                graupel += d_graup;
            }
            if (fallsum > 0.0) {
                sr = (snowncv + graupelncv) / (rainncv + 1.0e-12);
            }

            for (int k = 0; k < nz; ++k) {
                Real supcol = t0c - t[k];
                Real xlf = (supcol < 0.0) ? xlf0 : (xls - xl[k]);

                // pimlt
                if (supcol < 0.0 && qi[k] > 0.0) {
                    qc[k] += qi[k];
                    t[k] -= (xlf / cpm[k]) * qi[k];
                    qi[k] = 0.0;
                }
                // pihmf
                if (supcol > 40.0 && qc[k] > 0.0) {
                    qi[k] += qc[k];
                    t[k] += (xlf / cpm[k]) * qc[k];
                    qc[k] = 0.0;
                }
                // pihtf
                if (supcol > 0.0 && qc[k] > qmin) {
                    Real supcolt = std::min(supcol, static_cast<Real>(50.0));
                    Real pfrzdtc = std::min(pfrz1 * (std::exp(pfrz2 * supcolt) - static_cast<Real>(1.0))
                        * den[k] / denr / xncr * qc[k] * qc[k] * dtcld, qc[k]);
                    qi[k] += pfrzdtc;
                    t[k] += (xlf / cpm[k]) * pfrzdtc;
                    qc[k] -= pfrzdtc;
                }
                // pgfrz
                if (supcol > 0.0 && qr[k] > 0.0) {
                    Real temp = rslope3[k][0];
                    temp = temp * temp * rslope[k][0];
                    Real supcolt = std::min(supcol, static_cast<Real>(50.0));
                    Real pfrzdtr = std::min(static_cast<Real>(20.0) * (pi * pi) * pfrz1 * n0r * denr / den[k]
                        * (std::exp(pfrz2 * supcolt) - static_cast<Real>(1.0)) * temp * dtcld, qr[k]);
                    qg[k] += pfrzdtr;
                    t[k] += (xlf / cpm[k]) * pfrzdtr;
                    qr[k] -= pfrzdtr;
                }
            }

            slope_wsm6(nz, qr, qs, qg, den, denfac, t, rslope, rslopeb, rslope2, rslope3, vt);

            Real work1[kMaxLevels][2];
            Real work2[kMaxLevels];
            for (int k = 0; k < nz; ++k) {
                work1[k][0] = diffac(xl[k], p[k], t[k], den[k], qsat[k][0]);
                work1[k][1] = diffac(xls, p[k], t[k], den[k], qsat[k][1]);
                work2[k] = venfac(p[k], t[k], den[k]);
            }

            // Warm and cold rain microphysics
            for (int k = 0; k < nz; ++k) {
                Real supsat = std::max(q[k], qmin) - qsat[k][0];
                Real satdt = supsat / dtcld;

                Real praut = 0.0;
                Real pracw = 0.0;
                Real prevp = 0.0;

                // praut
                if (qc[k] > qc0) {
                    praut = qck1 * std::pow(qc[k], static_cast<Real>(7.0 / 3.0));
                    praut = std::min(praut, qc[k] / dtcld);
                }
                // pracw
                if (qr[k] > qcrmin && qc[k] > qmin) {
                    pracw = std::min(pacrr * rslope3[k][0] * rslopeb[k][0] * qc[k] * denfac[k], qc[k] / dtcld);
                }
                // prevp
                if (qr[k] > 0.0) {
                    Real coeres = rslope2[k][0] * std::sqrt(rslope[k][0] * rslopeb[k][0]);
                    prevp = (rh[k][0] - 1.0) * (precr1 * rslope2[k][0] + precr2 * work2[k] * coeres) / work1[k][0];
                    if (prevp < 0.0) {
                        prevp = std::max(prevp, -qr[k] / dtcld);
                        prevp = std::max(prevp, satdt / static_cast<Real>(2.0));
                    } else {
                        prevp = std::min(prevp, satdt / static_cast<Real>(2.0));
                    }
                }

                Real supcol = t0c - t[k];
                Real n0sfac_val = std::max(std::min(std::exp(alpha * supcol), n0smax / n0s), static_cast<Real>(1.0));
                Real supsat_ice = std::max(q[k], qmin) - qsat[k][1];
                Real satdt_ice = supsat_ice / dtcld;
                int ifsat = 0;

                Real temp = den[k] * std::max(qi[k], qmin);
                temp = std::sqrt(std::sqrt(temp * temp * temp));
                xni[k] = std::min(std::max(static_cast<Real>(5.38e7) * temp, static_cast<Real>(1.0e3)), static_cast<Real>(1.0e6));
                Real eacrs = std::exp(0.07 * (-supcol));

                Real xmi = den[k] * qi[k] / xni[k];
                Real diameter = std::min(dicon * std::sqrt(xmi), dimax);
                Real vt2i = 1.49e4 * std::pow(diameter, static_cast<Real>(1.31));
                Real vt2r = pvtr * rslopeb[k][0] * denfac[k];
                Real vt2s = pvts * rslopeb[k][1] * denfac[k];
                Real vt2g = pvtg * rslopeb[k][2] * denfac[k];

                Real qsum_val = std::max(qs[k] + qg[k], static_cast<Real>(1.0e-15));
                Real vt2ave = (qsum_val > static_cast<Real>(1.0e-15))
                    ? (vt2s * qs[k] + vt2g * qg[k]) / qsum_val : static_cast<Real>(0.0);

                Real praci = 0.0;
                Real piacr = 0.0;
                Real psaci = 0.0;
                Real pgaci = 0.0;

                if (supcol > 0.0 && qi[k] > qmin) {
                    if (qr[k] > qcrmin) {
                        Real acrfac = 2.0 * rslope3[k][0] + 2.0 * diameter * rslope2[k][0]
                            + diameter * diameter * rslope[k][0];
                        praci = pi * qi[k] * n0r * std::abs(vt2r - vt2i) * acrfac / 4.0;
                        Real eff = std::min(std::max(static_cast<Real>(0.0), qr[k] / qi[k]), static_cast<Real>(1.0));
                        praci = praci * (eff * eff);
                        praci = std::min(praci, qi[k] / dtcld);

                        piacr = (pi * pi) * avtr * n0r * denr * xni[k] * denfac[k]
                            * g6pbr * rslope3[k][0] * rslope3[k][0]
                            * rslopeb[k][0] / 24.0 / den[k];
                        Real eff2 = std::min(std::max(static_cast<Real>(0.0), qi[k] / qr[k]), static_cast<Real>(1.0));
                        piacr = piacr * (eff2 * eff2);
                        piacr = std::min(piacr, qr[k] / dtcld);
                    }
                    if (qs[k] > qcrmin) {
                        Real acrfac = 2.0 * rslope3[k][1] + 2.0 * diameter * rslope2[k][1]
                            + diameter * diameter * rslope[k][1];
                        psaci = pi * qi[k] * eacrs * n0s * n0sfac_val * std::abs(vt2ave - vt2i) * acrfac / 4.0;
                        psaci = std::min(psaci, qi[k] / dtcld);
                    }
                    if (qg[k] > qcrmin) {
                        Real egi = std::exp(0.07 * (-supcol));
                        Real acrfac = 2.0 * rslope3[k][2] + 2.0 * diameter * rslope2[k][2]
                            + diameter * diameter * rslope[k][2];
                        pgaci = pi * egi * qi[k] * n0g * std::abs(vt2ave - vt2i) * acrfac / 4.0;
                        pgaci = std::min(pgaci, qi[k] / dtcld);
                    }
                }

                Real psacw = 0.0;
                if (qs[k] > qcrmin && qc[k] > qmin) {
                    Real eff = std::min(std::max(static_cast<Real>(0.0), qs[k] / qc[k]), static_cast<Real>(1.0));
                    psacw = std::min(pacrc * n0sfac_val * rslope3[k][1] * rslopeb[k][1]
                        * (eff * eff) * qc[k] * denfac[k], qc[k] / dtcld);
                }

                Real pgacw = 0.0;
                if (qg[k] > qcrmin && qc[k] > qmin) {
                    Real eff = std::min(std::max(static_cast<Real>(0.0), qg[k] / qc[k]), static_cast<Real>(1.0));
                    pgacw = std::min(pacrg * rslope3[k][2] * rslopeb[k][2]
                        * (eff * eff) * qc[k] * denfac[k], qc[k] / dtcld);
                }

                Real paacw = 0.0;
                if (qsum_val > static_cast<Real>(1.0e-15)) {
                    paacw = (qs[k] * psacw + qg[k] * pgacw) / qsum_val;
                }

                Real pracs = 0.0;
                Real psacr = 0.0;
                if (qs[k] > qcrmin && qr[k] > qcrmin) {
                    if (supcol > 0.0) {
                        Real acrfac = 5.0 * rslope3[k][1] * rslope3[k][1] * rslope[k][0]
                            + 2.0 * rslope3[k][1] * rslope2[k][1] * rslope2[k][0]
                            + 0.5 * rslope2[k][1] * rslope2[k][1] * rslope3[k][0];
                        pracs = (pi * pi) * n0r * n0s * n0sfac_val * std::abs(vt2r - vt2ave)
                            * (dens / den[k]) * acrfac;
                        Real eff = std::min(std::max(static_cast<Real>(0.0), qr[k] / qs[k]), static_cast<Real>(1.0));
                        pracs = pracs * (eff * eff);
                        pracs = std::min(pracs, qs[k] / dtcld);
                    }
                    Real acrfac = 5.0 * rslope3[k][0] * rslope3[k][0] * rslope[k][1]
                        + 2.0 * rslope3[k][0] * rslope2[k][0] * rslope2[k][1]
                        + 0.5 * rslope2[k][0] * rslope2[k][0] * rslope3[k][1];
                    psacr = (pi * pi) * n0r * n0s * n0sfac_val * std::abs(vt2ave - vt2r)
                        * (denr / den[k]) * acrfac;
                    Real eff = std::min(std::max(static_cast<Real>(0.0), qs[k] / qr[k]), static_cast<Real>(1.0));
                    psacr = psacr * (eff * eff);
                    psacr = std::min(psacr, qr[k] / dtcld);
                }

                Real pgacr = 0.0;
                if (qg[k] > qcrmin && qr[k] > qcrmin) {
                    Real acrfac = 5.0 * rslope3[k][0] * rslope3[k][0] * rslope[k][2]
                        + 2.0 * rslope3[k][0] * rslope2[k][0] * rslope2[k][2]
                        + 0.5 * rslope2[k][0] * rslope2[k][0] * rslope3[k][2];
                    pgacr = (pi * pi) * n0r * n0g * std::abs(vt2ave - vt2r) * (denr / den[k]) * acrfac;
                    Real eff = std::min(std::max(static_cast<Real>(0.0), qg[k] / qr[k]), static_cast<Real>(1.0));
                    pgacr = pgacr * (eff * eff);
                    pgacr = std::min(pgacr, qr[k] / dtcld);
                }

                Real pgacs = 0.0;
                Real pseml = 0.0;
                Real pgeml = 0.0;
                if (supcol <= 0.0) {
                    Real xlf = xlf0;
                    if (qs[k] > 0.0) {
                        pseml = std::min(std::max(cliq * supcol * (paacw + psacr) / xlf, -qs[k] / dtcld), static_cast<Real>(0.0));
                    }
                    if (qg[k] > 0.0) {
                        pgeml = std::min(std::max(cliq * supcol * (paacw + pgacr) / xlf, -qg[k] / dtcld), static_cast<Real>(0.0));
                    }
                }

                Real pidep = 0.0;
                Real psdep = 0.0;
                Real pgdep = 0.0;
                Real pigen = 0.0;
                Real psaut = 0.0;
                Real pgaut = 0.0;

                if (supcol > 0.0) {
                    if (qi[k] > 0.0 && ifsat != 1) {
                        pidep = 4.0 * diameter * xni[k] * (rh[k][1] - 1.0) / work1[k][1];
                        Real supice = satdt_ice - prevp;
                        if (pidep < 0.0) {
                            pidep = std::max(std::max(pidep, satdt_ice / static_cast<Real>(2.0)), supice);
                            pidep = std::max(pidep, -qi[k] / dtcld);
                        } else {
                            pidep = std::min(std::min(pidep, satdt_ice / static_cast<Real>(2.0)), supice);
                        }
                        if (std::abs(prevp + pidep) >= std::abs(satdt_ice)) ifsat = 1;
                    }

                    if (qs[k] > 0.0 && ifsat != 1) {
                        Real coeres = rslope2[k][1] * std::sqrt(rslope[k][1] * rslopeb[k][1]);
                        psdep = (rh[k][1] - 1.0) * n0sfac_val * (precs1 * rslope2[k][1]
                            + precs2 * work2[k] * coeres) / work1[k][1];
                        Real supice = satdt_ice - prevp - pidep;
                        if (psdep < 0.0) {
                            psdep = std::max(psdep, -qs[k] / dtcld);
                            psdep = std::max(std::max(psdep, satdt_ice / static_cast<Real>(2.0)), supice);
                        } else {
                            psdep = std::min(std::min(psdep, satdt_ice / static_cast<Real>(2.0)), supice);
                        }
                        if (std::abs(prevp + pidep + psdep) >= std::abs(satdt_ice)) ifsat = 1;
                    }

                    if (qg[k] > 0.0 && ifsat != 1) {
                        Real coeres = rslope2[k][2] * std::sqrt(rslope[k][2] * rslopeb[k][2]);
                        pgdep = (rh[k][1] - 1.0) * (precg1 * rslope2[k][2]
                            + precg2 * work2[k] * coeres) / work1[k][1];
                        Real supice = satdt_ice - prevp - pidep - psdep;
                        if (pgdep < 0.0) {
                            pgdep = std::max(pgdep, -qg[k] / dtcld);
                            pgdep = std::max(std::max(pgdep, satdt_ice / static_cast<Real>(2.0)), supice);
                        } else {
                            pgdep = std::min(std::min(pgdep, satdt_ice / static_cast<Real>(2.0)), supice);
                        }
                        if (std::abs(prevp + pidep + psdep + pgdep) >= std::abs(satdt_ice)) ifsat = 1;
                    }

                    if (supsat_ice > 0.0 && ifsat != 1) {
                        Real supice = satdt_ice - prevp - pidep - psdep - pgdep;
                        Real xni0 = 1.0e3 * std::exp(0.1 * supcol);
                        Real roqi0 = 4.92e-11 * std::pow(xni0, static_cast<Real>(1.33));
                        pigen = std::max(static_cast<Real>(0.0), (roqi0 / den[k] - std::max(qi[k], static_cast<Real>(0.0))) / dtcld);
                        pigen = std::min(std::min(pigen, satdt_ice), supice);
                    }

                    if (qi[k] > 0.0) {
                        Real qimax = roqimax / den[k];
                        psaut = std::max(static_cast<Real>(0.0), (qi[k] - qimax) / dtcld);
                    }

                    if (qs[k] > 0.0) {
                        Real alpha2 = 1.0e-3 * std::exp(0.09 * (-supcol));
                        pgaut = std::min(std::max(static_cast<Real>(0.0), alpha2 * (qs[k] - qs0)), qs[k] / dtcld);
                    }
                }

                Real psevp = 0.0;
                Real pgevp = 0.0;
                if (supcol < 0.0) {
                    if (qs[k] > 0.0 && rh[k][0] < 1.0) {
                        Real coeres = rslope2[k][1] * std::sqrt(rslope[k][1] * rslopeb[k][1]);
                        psevp = (rh[k][0] - 1.0) * n0sfac_val * (precs1 * rslope2[k][1]
                            + precs2 * work2[k] * coeres) / work1[k][0];
                        psevp = std::min(std::max(psevp, -qs[k] / dtcld), static_cast<Real>(0.0));
                    }
                    if (qg[k] > 0.0 && rh[k][0] < 1.0) {
                        Real coeres = rslope2[k][2] * std::sqrt(rslope[k][2] * rslopeb[k][2]);
                        pgevp = (rh[k][0] - 1.0) * (precg1 * rslope2[k][2]
                            + precg2 * work2[k] * coeres) / work1[k][0];
                        pgevp = std::min(std::max(pgevp, -qg[k] / dtcld), static_cast<Real>(0.0));
                    }
                }

                // Mass conservation check and feedback
                Real delta2 = (qr[k] < 1.0e-4 && qs[k] < 1.0e-4) ? static_cast<Real>(1.0) : static_cast<Real>(0.0);
                Real delta3 = (qr[k] < 1.0e-4) ? static_cast<Real>(1.0) : static_cast<Real>(0.0);

                if (t[k] <= t0c) {
                    // cloud water
                    Real value = std::max(qmin, qc[k]);
                    Real source = (praut + pracw + paacw + paacw) * dtcld;
                    if (source > value) {
                        Real factor = value / source;
                        praut *= factor;
                        pracw *= factor;
                        paacw *= factor;
                    }

                    // cloud ice
                    value = std::max(qmin, qi[k]);
                    source = (psaut - pigen - pidep + praci + psaci + pgaci) * dtcld;
                    if (source > value) {
                        Real factor = value / source;
                        psaut *= factor;
                        pigen *= factor;
                        pidep *= factor;
                        praci *= factor;
                        psaci *= factor;
                        pgaci *= factor;
                    }

                    // rain
                    value = std::max(qmin, qr[k]);
                    source = (-praut - prevp - pracw + piacr + psacr + pgacr) * dtcld;
                    if (source > value) {
                        Real factor = value / source;
                        praut *= factor;
                        prevp *= factor;
                        pracw *= factor;
                        piacr *= factor;
                        psacr *= factor;
                        pgacr *= factor;
                    }

                    // snow
                    value = std::max(qmin, qs[k]);
                    source = -(psdep + psaut - pgaut + paacw + piacr * delta3
                        + praci * delta3 - pracs * (1.0 - delta2)
                        + psacr * delta2 + psaci - pgacs) * dtcld;
                    if (source > value) {
                        Real factor = value / source;
                        psdep *= factor;
                        psaut *= factor;
                        pgaut *= factor;
                        paacw *= factor;
                        piacr *= factor;
                        praci *= factor;
                        psaci *= factor;
                        pracs *= factor;
                        psacr *= factor;
                        pgacs *= factor;
                    }

                    // graupel
                    value = std::max(qmin, qg[k]);
                    source = -(pgdep + pgaut + piacr * (1.0 - delta3) + praci * (1.0 - delta3)
                        + psacr * (1.0 - delta2) + pracs * (1.0 - delta2)
                        + pgaci + paacw + pgacr + pgacs) * dtcld;
                    if (source > value) {
                        Real factor = value / source;
                        pgdep *= factor;
                        pgaut *= factor;
                        piacr *= factor;
                        praci *= factor;
                        psacr *= factor;
                        pracs *= factor;
                        paacw *= factor;
                        pgaci *= factor;
                        pgacr *= factor;
                        pgacs *= factor;
                    }

                    Real work2_val = -(prevp + psdep + pgdep + pigen + pidep);
                    q[k] += work2_val * dtcld;
                    qc[k] = std::max(qc[k] - (praut + pracw + paacw + paacw) * dtcld, static_cast<Real>(0.0));
                    qr[k] = std::max(qr[k] + (praut + pracw + prevp - piacr - pgacr - psacr) * dtcld, static_cast<Real>(0.0));
                    qi[k] = std::max(qi[k] - (psaut + praci + psaci + pgaci - pigen - pidep) * dtcld, static_cast<Real>(0.0));
                    qs[k] = std::max(qs[k] + (psdep + psaut + paacw - pgaut + piacr * delta3
                        + praci * delta3 + psaci - pgacs - pracs * (static_cast<Real>(1.0) - delta2) + psacr * delta2) * dtcld, static_cast<Real>(0.0));
                    qg[k] = std::max(qg[k] + (pgdep + pgaut + piacr * (static_cast<Real>(1.0) - delta3)
                        + praci * (static_cast<Real>(1.0) - delta3) + psacr * (static_cast<Real>(1.0) - delta2)
                        + pracs * (static_cast<Real>(1.0) - delta2) + pgaci + paacw + pgacr + pgacs) * dtcld, static_cast<Real>(0.0));

                    Real xlf = xls - xl[k];
                    Real xlwork2 = -xls * (psdep + pgdep + pidep + pigen)
                        - xl[k] * prevp - xlf * (piacr + paacw + paacw + pgacr + psacr);
                    t[k] -= xlwork2 / cpm[k] * dtcld;

                } else {
                    // t[k] > t0c
                    // cloud water
                    Real value = std::max(qmin, qc[k]);
                    Real source = (praut + pracw + paacw + paacw) * dtcld;
                    if (source > value) {
                        Real factor = value / source;
                        praut *= factor;
                        pracw *= factor;
                        paacw *= factor;
                    }

                    // rain
                    value = std::max(qmin, qr[k]);
                    source = (-paacw - praut + pseml + pgeml - pracw - paacw - prevp) * dtcld;
                    if (source > value) {
                        Real factor = value / source;
                        praut *= factor;
                        prevp *= factor;
                        pracw *= factor;
                        paacw *= factor;
                        pseml *= factor;
                        pgeml *= factor;
                    }

                    // snow
                    value = std::max(qcrmin, qs[k]);
                    source = (pgacs - pseml - psevp) * dtcld;
                    if (source > value) {
                        Real factor = value / source;
                        pgacs *= factor;
                        psevp *= factor;
                        pseml *= factor;
                    }

                    // graupel
                    value = std::max(qcrmin, qg[k]);
                    source = -(pgacs + pgevp + pgeml) * dtcld;
                    if (source > value) {
                        Real factor = value / source;
                        pgacs *= factor;
                        pgevp *= factor;
                        pgeml *= factor;
                    }

                    Real work2_val = -(prevp + psevp + pgevp);
                    q[k] += work2_val * dtcld;
                    qc[k] = std::max(qc[k] - (praut + pracw + paacw + paacw) * dtcld, static_cast<Real>(0.0));
                    qr[k] = std::max(qr[k] + (praut + pracw + prevp + paacw + paacw - pseml - pgeml) * dtcld, static_cast<Real>(0.0));
                    qs[k] = std::max(qs[k] + (psevp - pgacs + pseml) * dtcld, static_cast<Real>(0.0));
                    qg[k] = std::max(qg[k] + (pgacs + pgevp + pgeml) * dtcld, static_cast<Real>(0.0));

                    Real xlf = xls - xl[k];
                    Real xlwork2 = -xl[k] * (prevp + psevp + pgevp) - xlf * (pseml + pgeml);
                    t[k] -= xlwork2 / cpm[k] * dtcld;
                }
            }

            // Saturation adjustment (conden)
            for (int k = 0; k < nz; ++k) {
                Real tr = ttp / t[k];
                Real qs1 = psat * std::exp(std::log(tr) * xa) * std::exp(xb * (1.0 - tr));
                qs1 = std::min(qs1, static_cast<Real>(0.99) * p[k]);
                qs1 = ep2 * qs1 / (p[k] - qs1);
                qs1 = std::max(qs1, qmin);

                Real work1_cond = conden(t[k], q[k], qs1, xl[k], cpm[k]);
                Real pcond = std::min(std::max(work1_cond / dtcld, static_cast<Real>(0.0)),
                                      std::max(q[k], static_cast<Real>(0.0)) / dtcld);
                if (qc[k] > 0.0 && work1_cond < 0.0) {
                    pcond = std::max(work1_cond, -qc[k]) / dtcld;
                }
                q[k] -= pcond * dtcld;
                qc[k] = std::max(qc[k] + pcond * dtcld, static_cast<Real>(0.0));
                t[k] += pcond * xl[k] / cpm[k] * dtcld;
            }

            // Small value padding
            for (int k = 0; k < nz; ++k) {
                if (qc[k] <= qmin) qc[k] = 0.0;
                if (qi[k] <= qmin) qi[k] = 0.0;
            }
        } // end loops

        for (int k = 0; k < nz; ++k) {
            th[k] = t[k] / pii[k];
        }
    }
};

template <typename Real>
struct Wsm6Kernel {
    Wsm6Params<Real> params;
    Real dt;
    int nz;

    View3<Real> th;
    View3<Real> q;
    View3<Real> qc;
    View3<Real> qr;
    View3<Real> qi;
    View3<Real> qs;
    View3<Real> qg;

    View3<const Real> den;
    View3<const Real> pii;
    View3<const Real> p;
    View3<const Real> delz;

    View2<Real> rain;
    View2<Real> rainncv;
    View2<Real> snow;
    View2<Real> snowncv;
    View2<Real> graupel;
    View2<Real> graupelncv;
    View2<Real> sr;

    CPPWRF_HD void operator()(int i, int j) const {
        Real th_col[kMaxLevels];
        Real q_col[kMaxLevels];
        Real qc_col[kMaxLevels];
        Real qr_col[kMaxLevels];
        Real qi_col[kMaxLevels];
        Real qs_col[kMaxLevels];
        Real qg_col[kMaxLevels];

        Real den_col[kMaxLevels];
        Real pii_col[kMaxLevels];
        Real p_col[kMaxLevels];
        Real delz_col[kMaxLevels];

        for (int k = 0; k < nz; ++k) {
            th_col[k] = th(i, j, k);
            q_col[k] = q(i, j, k);
            qc_col[k] = qc(i, j, k);
            qr_col[k] = qr(i, j, k);
            qi_col[k] = qi(i, j, k);
            qs_col[k] = qs(i, j, k);
            qg_col[k] = qg(i, j, k);

            den_col[k] = den(i, j, k);
            pii_col[k] = pii(i, j, k);
            p_col[k] = p(i, j, k);
            delz_col[k] = delz(i, j, k);
        }

        Real rain_val = rain(i, j);
        Real rainncv_val = rainncv(i, j);
        Real snow_val = snow(i, j);
        Real snowncv_val = snowncv(i, j);
        Real graupel_val = graupel(i, j);
        Real graupelncv_val = graupelncv(i, j);
        Real sr_val = sr(i, j);

        params(nz, dt, th_col, q_col, qc_col, qr_col, qi_col, qs_col, qg_col,
               den_col, pii_col, p_col, delz_col,
               rain_val, rainncv_val, snow_val, snowncv_val,
               graupel_val, graupelncv_val, sr_val);

        for (int k = 0; k < nz; ++k) {
            th(i, j, k) = th_col[k];
            q(i, j, k) = q_col[k];
            qc(i, j, k) = qc_col[k];
            qr(i, j, k) = qr_col[k];
            qi(i, j, k) = qi_col[k];
            qs(i, j, k) = qs_col[k];
            qg(i, j, k) = qg_col[k];
        }

        rain(i, j) = rain_val;
        rainncv(i, j) = rainncv_val;
        snow(i, j) = snow_val;
        snowncv(i, j) = snowncv_val;
        graupel(i, j) = graupel_val;
        graupelncv(i, j) = graupelncv_val;
        sr(i, j) = sr_val;
    }
};

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
    const Wsm6Params<real>& params);

} // namespace cppwrf::phys::wsm6

