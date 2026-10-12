#pragma once

#include "core/backend.hpp"
#include "core/field.hpp"
#include "phys/phys_constants.hpp"
#include "phys/cu_kf/kf_tables.hpp"
#include <cmath>
#include <cassert>

namespace cppwrf::phys::cu_kf {

inline constexpr int kMaxLevels = 128;

template <typename T>
CPPWRF_HD inline void envirtht(const KFTables<T>& tables, T p1, T t1, T q1, T& tht1,
                               T aliq, T bliq, T cliq, T dliq) {
    using std::exp;
    using std::pow;

    const T t00 = static_cast<T>(273.16);
    const T p00 = static_cast<T>(1.0e5);
    const T c1 = static_cast<T>(3374.6525);
    const T c2 = static_cast<T>(2.5403);

    T ee = q1 * p1 / (static_cast<T>(0.622) + q1);
    const T astrt = static_cast<T>(1.0e-3);
    const T ainc = static_cast<T>(0.075);
    T a1 = ee / aliq;
    T tp = (a1 - astrt) / ainc;
    int indlu = static_cast<int>(tp) + 1; // 1-based Fortran
    if (indlu < 1) indlu = 1;
    if (indlu > 199) indlu = 199;
    T value = static_cast<T>(indlu - 1) * ainc + astrt;
    T aintrp = (a1 - value) / ainc;
    T tlog = aintrp * tables.alu[indlu] + (static_cast<T>(1.0) - aintrp) * tables.alu[indlu - 1];

    T tdpt = (cliq - dliq * tlog) / (bliq - tlog);
    T tsat = tdpt - (static_cast<T>(0.212) + static_cast<T>(1.571e-3) * (tdpt - t00) -
                     static_cast<T>(4.36e-4) * (t1 - t00)) * (t1 - tdpt);
    T tht = t1 * pow(p00 / p1, static_cast<T>(0.2854) * (static_cast<T>(1.0) - static_cast<T>(0.28) * q1));
    tht1 = tht * exp((c1 / tsat - c2) * q1 * (static_cast<T>(1.0) + static_cast<T>(0.81) * q1));
}

template <typename T>
CPPWRF_HD inline void tpmix2(const KFTables<T>& tables, T p, T thes, T& tu, T& qu,
                             T& qliq, T& qice, T& qnewlq, T& qnewic, T xlv1, T xlv0) {
    T tp = (p - tables.plutop) * tables.rdpr;
    int ip = static_cast<int>(tp);
    if (ip < 0) ip = 0;
    if (ip > KFNP - 2) ip = KFNP - 2;
    T qq = tp - static_cast<T>(static_cast<int>(tp));

    T bth = (tables.the0k[ip + 1] - tables.the0k[ip]) * qq + tables.the0k[ip];
    T tth = (thes - bth) * tables.rdthk;
    int it = static_cast<int>(tth);
    if (it < 0) it = 0;
    if (it > KFNT - 2) it = KFNT - 2;
    T pp = tth - static_cast<T>(static_cast<int>(tth));

    T t00 = tables.ttab[ip][it];
    T t10 = tables.ttab[ip][it + 1];
    T t01 = tables.ttab[ip + 1][it];
    T t11 = tables.ttab[ip + 1][it + 1];

    T q00 = tables.qstab[ip][it];
    T q10 = tables.qstab[ip][it + 1];
    T q01 = tables.qstab[ip + 1][it];
    T q11 = tables.qstab[ip + 1][it + 1];

    T temp = (t00 + (t10 - t00) * pp + (t01 - t00) * qq + (t00 - t10 - t01 + t11) * pp * qq);
    T qs = (q00 + (q10 - q00) * pp + (q01 - q00) * qq + (q00 - q10 - q01 + q11) * pp * qq);

    T dq = qs - qu;
    T qnew;
    if (dq <= static_cast<T>(0.0)) {
        qnew = qu - qs;
        qu = qs;
    } else {
        qnew = static_cast<T>(0.0);
        T qtot = qliq + qice;
        if (qtot >= dq) {
            qliq = qliq - dq * qliq / (qtot + static_cast<T>(1.0e-10));
            qice = qice - dq * qice / (qtot + static_cast<T>(1.0e-10));
            qu = qs;
        } else {
            T rll = xlv0 - xlv1 * temp;
            T cpp = static_cast<T>(1004.5) * (static_cast<T>(1.0) + static_cast<T>(0.89) * qu);
            if (qtot < static_cast<T>(1.0e-10)) {
                temp = temp + rll * (dq / (static_cast<T>(1.0) + dq)) / cpp;
            } else {
                temp = temp + rll * ((dq - qtot) / (static_cast<T>(1.0) + dq - qtot)) / cpp;
                qu = qu + qtot;
                qtot = static_cast<T>(0.0);
                qliq = static_cast<T>(0.0);
                qice = static_cast<T>(0.0);
            }
        }
    }
    tu = temp;
    qnewlq = qnew;
    qnewic = static_cast<T>(0.0);
}

template <typename T>
CPPWRF_HD inline void tpmix2dd(const KFTables<T>& tables, T p, T thes, T& ts, T& qs) {
    T tp = (p - tables.plutop) * tables.rdpr;
    int ip = static_cast<int>(tp);
    if (ip < 0) ip = 0;
    if (ip > KFNP - 2) ip = KFNP - 2;
    T qq = tp - static_cast<T>(static_cast<int>(tp));

    T bth = (tables.the0k[ip + 1] - tables.the0k[ip]) * qq + tables.the0k[ip];
    T tth = (thes - bth) * tables.rdthk;
    int it = static_cast<int>(tth);
    if (it < 0) it = 0;
    if (it > KFNT - 2) it = KFNT - 2;
    T pp = tth - static_cast<T>(static_cast<int>(tth));

    T t00 = tables.ttab[ip][it];
    T t10 = tables.ttab[ip][it + 1];
    T t01 = tables.ttab[ip + 1][it];
    T t11 = tables.ttab[ip + 1][it + 1];

    T q00 = tables.qstab[ip][it];
    T q10 = tables.qstab[ip][it + 1];
    T q01 = tables.qstab[ip + 1][it];
    T q11 = tables.qstab[ip + 1][it + 1];

    ts = (t00 + (t10 - t00) * pp + (t01 - t00) * qq + (t00 - t10 - t01 + t11) * pp * qq);
    qs = (q00 + (q10 - q00) * pp + (q01 - q00) * qq + (q00 - q10 - q01 + q11) * pp * qq);
}

template <typename T>
CPPWRF_HD inline void dtfrznew(T& tu, T p, T& thteu, T& qu, T qfrz, T& qice,
                               T aliq, T bliq, T cliq, T dliq) {
    using std::exp;
    using std::pow;

    T rlc = static_cast<T>(2.5e6) - static_cast<T>(2369.276) * (tu - static_cast<T>(273.16));
    T rls = static_cast<T>(2833922.0) - static_cast<T>(259.532) * (tu - static_cast<T>(273.16));
    T rlf = rls - rlc;
    T cpp = static_cast<T>(1004.5) * (static_cast<T>(1.0) + static_cast<T>(0.89) * qu);

    T a = (cliq - bliq * dliq) / ((tu - dliq) * (tu - dliq));
    T dtfrz = rlf * qfrz / (cpp + rls * qu * a);
    tu = tu + dtfrz;

    T es = aliq * exp((bliq * tu - cliq) / (tu - dliq));
    T qs = es * static_cast<T>(0.622) / (p - es);

    T dqevap = qs - qu;
    qice = qice - dqevap;
    qu = qu + dqevap;
    T pii = pow(static_cast<T>(1.0e5) / p, static_cast<T>(0.2854) * (static_cast<T>(1.0) - static_cast<T>(0.28) * qu));
    thteu = tu * pii * exp((static_cast<T>(3374.6525) / tu - static_cast<T>(2.5403)) * qu * (static_cast<T>(1.0) + static_cast<T>(0.81) * qu));
}

template <typename T>
CPPWRF_HD inline void condload(T& qliq, T& qice, T& wtw, T dz, T boterm, T enterm,
                               T rate, T& qnewlq, T& qnewic, T& qlqout, T& qicout, T g) {
    using std::sqrt;
    using std::exp;
    using std::abs;

    T qtot = qliq + qice;
    T qnew = qnewlq + qnewic;
    T qest = static_cast<T>(0.5) * (qtot + qnew);
    T g1 = wtw + boterm - enterm - static_cast<T>(2.0) * g * dz * qest / static_cast<T>(1.5);
    if (g1 < static_cast<T>(0.0)) g1 = static_cast<T>(0.0);
    T wavg = static_cast<T>(0.5) * (sqrt(wtw) + sqrt(g1));
    T conv = rate * dz / wavg;

    T ratio3 = qnewlq / (qnew + static_cast<T>(1.0e-8));
    qtot = qtot + static_cast<T>(0.6) * qnew;
    T oldq = qtot;
    T ratio4 = (static_cast<T>(0.6) * qnewlq + qliq) / (qtot + static_cast<T>(1.0e-8));
    qtot = qtot * exp(-conv);

    T dq = oldq - qtot;
    qlqout = ratio4 * dq;
    qicout = (static_cast<T>(1.0) - ratio4) * dq;

    T pptdrg = static_cast<T>(0.5) * (oldq + qtot - static_cast<T>(0.2) * qnew);
    wtw = wtw + boterm - enterm - static_cast<T>(2.0) * g * dz * pptdrg / static_cast<T>(1.5);
    if (abs(wtw) < static_cast<T>(1.0e-4)) wtw = static_cast<T>(1.0e-4);

    qliq = ratio4 * qtot + ratio3 * static_cast<T>(0.4) * qnew;
    qice = (static_cast<T>(1.0) - ratio4) * qtot + (static_cast<T>(1.0) - ratio3) * static_cast<T>(0.4) * qnew;
    qnewlq = static_cast<T>(0.0);
    qnewic = static_cast<T>(0.0);
}

template <typename T>
CPPWRF_HD inline void prof5(T eq, T& ee, T& ud) {
    using std::exp;
    using std::abs;

    const T sqrt2p = static_cast<T>(2.506628);
    const T a1 = static_cast<T>(0.4361836);
    const T a2 = static_cast<T>(-0.1201676);
    const T a3 = static_cast<T>(0.9372980);
    const T p = static_cast<T>(0.33267);
    const T sigma = static_cast<T>(0.166666667);
    const T fe = static_cast<T>(0.202765151);

    T y = static_cast<T>(6.0) * eq - static_cast<T>(3.0);
    T ey = exp(y * y / static_cast<T>(-2.0));
    T e45 = exp(static_cast<T>(-4.5));
    T t2 = static_cast<T>(1.0) / (static_cast<T>(1.0) + p * abs(y));
    T t1 = static_cast<T>(0.500498);
    T c1 = a1 * t1 + a2 * t1 * t1 + a3 * t1 * t1 * t1;
    T c2 = a1 * t2 + a2 * t2 * t2 + a3 * t2 * t2 * t2;

    if (y >= static_cast<T>(0.0)) {
        ee = sigma * (static_cast<T>(0.5) * (sqrt2p - e45 * c1 - ey * c2) + sigma * (e45 - ey)) -
             e45 * eq * eq / static_cast<T>(2.0);
        ud = sigma * (static_cast<T>(0.5) * (ey * c2 - e45 * c1) + sigma * (e45 - ey)) -
             e45 * (static_cast<T>(0.5) + eq * eq / static_cast<T>(2.0) - eq);
    } else {
        ee = sigma * (static_cast<T>(0.5) * (ey * c2 - e45 * c1) + sigma * (e45 - ey)) -
             e45 * eq * eq / static_cast<T>(2.0);
        ud = sigma * (static_cast<T>(0.5) * (sqrt2p - e45 * c1 - ey * c2) + sigma * (e45 - ey)) -
             e45 * (static_cast<T>(0.5) + eq * eq / static_cast<T>(2.0) - eq);
    }
    ee = ee / fe;
    ud = ud / fe;
}

template <typename T>
struct KFKernel {
    // Grid sizes
    int nx;
    int ny;
    int nz;

    // Scalars
    T dt;
    T dx;
    T cudt;
    int stepcu;
    int ktau;
    int trigger;
    bool warm_rain;
    bool flag_qi;
    bool flag_qs;
    bool flag_qr;

    // Physical constants
    T xlv0;
    T xlv1;
    T xls0;
    T xls1;
    T cp;
    T r;
    T g;
    T ep1;
    T ep2;
    T svp1;
    T svp2;
    T svp3;
    T svpt0;

    // Lookup tables
    const KFTables<T>* tables;

    // 3D Inputs (nx, ny, nz)
    View3<const T> u;
    View3<const T> v;
    View3<const T> w; // (nx, ny, nz+1)
    View3<const T> th;
    View3<const T> t;
    View3<const T> qv;
    View3<const T> dz8w;
    View3<const T> pcps;
    View3<const T> rho;
    View3<const T> pi;

    // 3D Inout
    View3<T> w0avg;

    // 2D Inout
    View2<T> raincv;
    View2<T> pratec;
    View2<T> nca;

    // 2D Out
    View2<T> cutop;
    View2<T> cubot;

    // 3D Tendency Outputs
    View3<T> rthcuten;
    View3<T> rqvcuten;
    View3<T> rqccuten;
    View3<T> rqrcuten;
    View3<T> rqicuten;
    View3<T> rqscuten;

    // Diagnostics
    View3<T> cldfra_dp_kf;
    View3<T> cldfra_sh_kf;
    View3<T> qc_kf;
    View3<T> qi_kf;

    CPPWRF_HD void operator()(int i, int j) const {
        using std::sqrt;
        using std::exp;
        using std::log;
        using std::pow;
        using std::abs;

        const int kx = nz;
        const int kl = nz;
        const T dxsq = dx * dx;

        // 1. Running mean vertical velocity W0AVG
        const T ntst = static_cast<T>(stepcu);
        const T tst = ntst * static_cast<T>(2.0);
        const T w0avgfctr = tst - static_cast<T>(1.0);
        const T w0fctr = static_cast<T>(1.0);
        const T w0den = tst;

        for (int k = 0; k < kx; ++k) {
            T w0 = static_cast<T>(0.5) * (w(i, j, k) + w(i, j, k + 1));
            w0avg(i, j, k) = (w0avg(i, j, k) * w0avgfctr + w0 * w0fctr) / w0den;
        }

        // 2. Check NCA relaxation condition
        if (nca(i, j) >= static_cast<T>(0.5) * dt) {
            // Convection active from prior steps, skip new initiation
            return;
        }

        // Initialize tendencies and outputs
        for (int k = 0; k < kx; ++k) {
            rthcuten(i, j, k) = static_cast<T>(0.0);
            rqvcuten(i, j, k) = static_cast<T>(0.0);
            rqccuten(i, j, k) = static_cast<T>(0.0);
            rqrcuten(i, j, k) = static_cast<T>(0.0);
            if (flag_qi) rqicuten(i, j, k) = static_cast<T>(0.0);
            if (flag_qs) rqscuten(i, j, k) = static_cast<T>(0.0);
            cldfra_dp_kf(i, j, k) = static_cast<T>(0.0);
            cldfra_sh_kf(i, j, k) = static_cast<T>(0.0);
            qc_kf(i, j, k) = static_cast<T>(0.0);
            qi_kf(i, j, k) = static_cast<T>(0.0);
        }
        raincv(i, j) = static_cast<T>(0.0);
        pratec(i, j) = static_cast<T>(0.0);
        cutop(i, j) = static_cast<T>(1.0);          // KTS = 1
        cubot(i, j) = static_cast<T>(kx + 1);       // KTE + 1

        // Local column scratch arrays
        T u0[kMaxLevels], v0[kMaxLevels], t0[kMaxLevels], qv0[kMaxLevels], p0[kMaxLevels];
        T rhoe[kMaxLevels], dzq[kMaxLevels], w0avg1d[kMaxLevels];
        T tpart_h0[kMaxLevels], tpart_v0[kMaxLevels];

        for (int k = 0; k < kx; ++k) {
            u0[k] = u(i, j, k);
            v0[k] = v(i, j, k);
            t0[k] = t(i, j, k);
            qv0[k] = qv(i, j, k);
            p0[k] = pcps(i, j, k);
            rhoe[k] = rho(i, j, k);
            dzq[k] = dz8w(i, j, k);
            w0avg1d[k] = w0avg(i, j, k);
            tpart_h0[k] = static_cast<T>(0.0);
            tpart_v0[k] = static_cast<T>(0.0);
        }

        // --- Call KF_eta_PARA logic ---
        const T p00 = static_cast<T>(1.0e5);
        const T t00 = static_cast<T>(273.16);
        const T rlf = static_cast<T>(3.339e5);
        const T ttfrz = static_cast<T>(268.16);
        const T tbfrz = static_cast<T>(248.16);
        const T rate = static_cast<T>(0.03);

        const T gdry = -g / cp;
        const T rocp = r / cp;
        const T aliq = svp1 * static_cast<T>(1000.0);
        const T bliq = svp2;
        const T cliq = svp2 * svpt0;
        const T dliq = svp3;

        T fbfrc = static_cast<T>(0.0);
        int ishall = 0;
        const T dpmin = static_cast<T>(5.0e3);
        const T p300 = p0[0] - static_cast<T>(30000.0);

        T q0[kMaxLevels], z0[kMaxLevels], tv0[kMaxLevels], tu[kMaxLevels], tvu[kMaxLevels];
        T qu[kMaxLevels], tz[kMaxLevels], tvd[kMaxLevels], qd[kMaxLevels], qes[kMaxLevels];
        T tg[kMaxLevels], tvg[kMaxLevels], qg[kMaxLevels], wu[kMaxLevels], wd[kMaxLevels];
        T ems[kMaxLevels], emsd[kMaxLevels], umf[kMaxLevels], uer[kMaxLevels], udr[kMaxLevels];
        T dmf[kMaxLevels], der[kMaxLevels], ddr[kMaxLevels], umf2[kMaxLevels], uer2[kMaxLevels];
        T udr2[kMaxLevels], dmf2[kMaxLevels], der2[kMaxLevels], ddr2[kMaxLevels], dza[kMaxLevels];
        T thta0[kMaxLevels], thetee[kMaxLevels], thtau[kMaxLevels], theteu[kMaxLevels];
        T thtad[kMaxLevels], theted[kMaxLevels], qliq[kMaxLevels], qice[kMaxLevels];
        T qlqout[kMaxLevels], qicout[kMaxLevels], pptliq[kMaxLevels], pptice[kMaxLevels];
        T detlq[kMaxLevels], detic[kMaxLevels], detlq2[kMaxLevels], detic2[kMaxLevels];
        T ratio2[kMaxLevels], domgdp[kMaxLevels], exn[kMaxLevels], tvqu[kMaxLevels];
        T dp[kMaxLevels], rh[kMaxLevels], eqfrc[kMaxLevels], wspd[kMaxLevels];
        T qdt[kMaxLevels], fxm[kMaxLevels], thtag[kMaxLevels], thpa[kMaxLevels];
        T thfxout[kMaxLevels], thfxin[kMaxLevels], qpa[kMaxLevels], qfxout[kMaxLevels];
        T qfxin[kMaxLevels], qlpa[kMaxLevels], qlfxin[kMaxLevels], qlfxout[kMaxLevels];
        T qipa[kMaxLevels], qifxin[kMaxLevels], qifxout[kMaxLevels], qrpa[kMaxLevels];
        T qrfxin[kMaxLevels], qrfxout[kMaxLevels], qspa[kMaxLevels], qsfxin[kMaxLevels];
        T qsfxout[kMaxLevels], ql0[kMaxLevels], qlg[kMaxLevels], qi0[kMaxLevels];
        T qig[kMaxLevels], qr0[kMaxLevels], qrg[kMaxLevels], qs0[kMaxLevels], qsg[kMaxLevels];
        T rainfb[kMaxLevels], snowfb[kMaxLevels], cldhgt[kMaxLevels], qsd[kMaxLevels];
        T dilfrc[kMaxLevels], ddilfrc[kMaxLevels], tgu[kMaxLevels], qgu[kMaxLevels], thteeg[kMaxLevels];
        T omg[kMaxLevels + 1];
        int kcheck[kMaxLevels];

        int l5 = 0;
        int llfc = 0;
        for (int k = 0; k < kx; ++k) {
            T es = aliq * exp((bliq * t0[k] - cliq) / (t0[k] - dliq));
            qes[k] = static_cast<T>(0.622) * es / (p0[k] - es);
            q0[k] = (qes[k] < qv0[k]) ? qes[k] : qv0[k];
            if (q0[k] < static_cast<T>(0.000001)) q0[k] = static_cast<T>(0.000001);
            ql0[k] = static_cast<T>(0.0);
            qi0[k] = static_cast<T>(0.0);
            qr0[k] = static_cast<T>(0.0);
            qs0[k] = static_cast<T>(0.0);
            rh[k] = q0[k] / qes[k];
            dilfrc[k] = static_cast<T>(1.0);
            tv0[k] = t0[k] * (static_cast<T>(1.0) + static_cast<T>(0.608) * q0[k]);
            dp[k] = rhoe[k] * g * dzq[k];
            cldhgt[k] = static_cast<T>(0.0);
            if (p0[k] >= static_cast<T>(0.5) * p0[0]) l5 = k;
            if (p0[k] >= p300) llfc = k;
        }

        z0[0] = static_cast<T>(0.5) * dzq[0];
        for (int k = 1; k < kl; ++k) {
            z0[k] = z0[k - 1] + static_cast<T>(0.5) * (dzq[k] + dzq[k - 1]);
            dza[k - 1] = z0[k] - z0[k - 1];
        }
        dza[kl - 1] = static_cast<T>(0.0);

        int ncheck = 1;
        kcheck[0] = 0;
        T pm15 = p0[0] - static_cast<T>(15.0e2);
        for (int k = 1; k <= llfc; ++k) {
            if (p0[k] < pm15) {
                kcheck[ncheck] = k;
                ncheck++;
                pm15 -= static_cast<T>(15.0e2);
            }
        }

        int nu = -1;
        int nuchm = -1;
        int lc = 0;
        int kpbl = 0;
        int klcl = 0;
        int ltop = 0;
        int let = 0;
        int lcl = 0;
        T zlcl = static_cast<T>(0.0);
        T tlcl = static_cast<T>(0.0);
        T tvlcl = static_cast<T>(0.0);
        T plcl = static_cast<T>(0.0);
        T wlcl = static_cast<T>(0.0);
        T wkl = static_cast<T>(0.0);
        T dpthmx = static_cast<T>(0.0);
        T vmflcl = static_cast<T>(0.0);
        T upold = static_cast<T>(0.0);
        T upnew = static_cast<T>(0.0);
        T trppt = static_cast<T>(0.0);
        T au0 = static_cast<T>(0.0);
        T tmix = static_cast<T>(0.0);
        T qmix = static_cast<T>(0.0);
        T zmix = static_cast<T>(0.0);
        T pmix = static_cast<T>(0.0);
        T abe = static_cast<T>(0.0);

        // USL search loop
        while (true) {
            nu++;
            if (nu >= ncheck) {
                if (ishall == 1) {
                    T chmax = static_cast<T>(0.0);
                    for (int nk = 0; nk < ncheck; ++nk) {
                        int nnn = kcheck[nk];
                        if (cldhgt[nnn] > chmax) {
                            nuchm = nk;
                            chmax = cldhgt[nnn];
                        }
                    }
                    nu = nuchm - 1;
                    fbfrc = static_cast<T>(1.0);
                    continue;
                } else {
                    return; // No convection occurs
                }
            }

            int kmix = kcheck[nu];
            lc = kmix;

            int nlayrs = 0;
            dpthmx = static_cast<T>(0.0);
            for (int nk = lc; nk < kx; ++nk) {
                dpthmx += dp[nk];
                nlayrs++;
                if (dpthmx > dpmin) break;
            }
            if (dpthmx < dpmin) {
                return;
            }
            kpbl = lc + nlayrs - 1;

            tmix = static_cast<T>(0.0);
            qmix = static_cast<T>(0.0);
            zmix = static_cast<T>(0.0);
            pmix = static_cast<T>(0.0);
            for (int nk = lc; nk <= kpbl; ++nk) {
                tmix += dp[nk] * t0[nk];
                qmix += dp[nk] * q0[nk];
                zmix += dp[nk] * z0[nk];
                pmix += dp[nk] * p0[nk];
            }
            tmix /= dpthmx;
            qmix /= dpthmx;
            zmix /= dpthmx;
            pmix /= dpthmx;
            T emix = qmix * pmix / (static_cast<T>(0.622) + qmix);

            const T astrt = static_cast<T>(1.0e-3);
            const T ainc_lut = static_cast<T>(0.075);
            T a1 = emix / aliq;
            T tp = (a1 - astrt) / ainc_lut;
            int indlu = static_cast<int>(tp) + 1;
            if (indlu < 1) indlu = 1;
            if (indlu > 199) indlu = 199;
            T value = static_cast<T>(indlu - 1) * ainc_lut + astrt;
            T aintrp = (a1 - value) / ainc_lut;
            T tlog = aintrp * tables->alu[indlu] + (static_cast<T>(1.0) - aintrp) * tables->alu[indlu - 1];
            T tdpt = (cliq - dliq * tlog) / (bliq - tlog);
            tlcl = tdpt - (static_cast<T>(0.212) + static_cast<T>(1.571e-3) * (tdpt - t00) -
                           static_cast<T>(4.36e-4) * (tmix - t00)) * (tmix - tdpt);
            if (tlcl > tmix) tlcl = tmix;
            tvlcl = tlcl * (static_cast<T>(1.0) + static_cast<T>(0.608) * qmix);
            zlcl = zmix + (tlcl - tmix) / gdry;

            klcl = lc;
            for (int nk = lc; nk < kl; ++nk) {
                klcl = nk;
                if (zlcl <= z0[nk]) break;
            }
            if (zlcl > z0[kl - 1]) {
                return;
            }

            int k_lcl = klcl - 1;
            T dlp = (zlcl - z0[k_lcl]) / (z0[klcl] - z0[k_lcl]);
            T tenv = t0[k_lcl] + (t0[klcl] - t0[k_lcl]) * dlp;
            T qenv = q0[k_lcl] + (q0[klcl] - q0[k_lcl]) * dlp;
            T tven = tenv * (static_cast<T>(1.0) + static_cast<T>(0.608) * qenv);

            T wklcl = (zlcl < static_cast<T>(2.0e3)) ? static_cast<T>(0.02) * zlcl / static_cast<T>(2.0e3)
                                                     : static_cast<T>(0.02);
            wkl = (w0avg1d[k_lcl] + (w0avg1d[klcl] - w0avg1d[k_lcl]) * dlp) * dx / static_cast<T>(25.0e3) - wklcl;
            T dtlcl = (wkl < static_cast<T>(0.0001)) ? static_cast<T>(0.0) : static_cast<T>(4.64) * pow(wkl, static_cast<T>(0.33));
            if (trigger == 2) {
                T sum = tpart_h0[klcl] + tpart_v0[klcl];
                dtlcl = (sum > static_cast<T>(0.0)) ? sum : static_cast<T>(0.0);
            }
            T dtrh = static_cast<T>(0.0);
            if (trigger == 3) {
                T qslcl = qes[k_lcl] + (qes[klcl] - qes[k_lcl]) * dlp;
                T rhlcl = qenv / qslcl;
                T dqssdt = qmix * (cliq - bliq * dliq) / ((tlcl - dliq) * (tlcl - dliq));
                if (rhlcl >= static_cast<T>(0.75) && rhlcl <= static_cast<T>(0.95)) {
                    dtrh = static_cast<T>(0.25) * (rhlcl - static_cast<T>(0.75)) * qmix / dqssdt;
                } else if (rhlcl > static_cast<T>(0.95)) {
                    dtrh = (static_cast<T>(1.0) / rhlcl - static_cast<T>(1.0)) * qmix / dqssdt;
                }
            }

            // Check trigger condition
            if (tlcl + dtlcl + dtrh < tenv) {
                continue; // Not buoyant
            }

            // Parcel buoyant
            envirtht(*tables, pmix, tmix, qmix, theteu[k_lcl], aliq, bliq, cliq, dliq);
            T dttot = dtlcl + dtrh;
            if (dttot > static_cast<T>(1.0e-4)) {
                T gdt = static_cast<T>(2.0) * g * dttot * static_cast<T>(500.0) / tven;
                wlcl = static_cast<T>(1.0) + static_cast<T>(0.5) * sqrt(gdt);
                if (wlcl > static_cast<T>(3.0)) wlcl = static_cast<T>(3.0);
            } else {
                wlcl = static_cast<T>(1.0);
            }
            plcl = p0[k_lcl] + (p0[klcl] - p0[k_lcl]) * dlp;
            T wtw = wlcl * wlcl;
            tvlcl = tlcl * (static_cast<T>(1.0) + static_cast<T>(0.608) * qmix);
            T rholcl = plcl / (r * tvlcl);

            lcl = klcl;
            let = lcl;
            T rad = (wkl < static_cast<T>(0.0)) ? static_cast<T>(1000.0) :
                    (wkl > static_cast<T>(0.1)) ? static_cast<T>(2000.0) :
                    (static_cast<T>(1000.0) + static_cast<T>(1000.0) * wkl / static_cast<T>(0.1));

            wu[k_lcl] = wlcl;
            au0 = static_cast<T>(0.01) * dxsq;
            umf[k_lcl] = rholcl * au0;
            vmflcl = umf[k_lcl];
            upold = vmflcl;
            upnew = upold;

            ratio2[k_lcl] = static_cast<T>(0.0);
            uer[k_lcl] = static_cast<T>(0.0);
            abe = static_cast<T>(0.0);
            trppt = static_cast<T>(0.0);
            tu[k_lcl] = tlcl;
            tvu[k_lcl] = tvlcl;
            qu[k_lcl] = qmix;
            eqfrc[k_lcl] = static_cast<T>(1.0);
            qliq[k_lcl] = static_cast<T>(0.0);
            qice[k_lcl] = static_cast<T>(0.0);
            qlqout[k_lcl] = static_cast<T>(0.0);
            qicout[k_lcl] = static_cast<T>(0.0);
            detlq[k_lcl] = static_cast<T>(0.0);
            detic[k_lcl] = static_cast<T>(0.0);
            pptliq[k_lcl] = static_cast<T>(0.0);
            pptice[k_lcl] = static_cast<T>(0.0);

            T ttemp = ttfrz;
            T ee1 = static_cast<T>(1.0);
            T ud1 = static_cast<T>(0.0);
            T rei = static_cast<T>(0.0);
            T dilbe = static_cast<T>(0.0);

            // Updraft loop
            ltop = kl - 1;
            for (int nk = k_lcl; nk < kl - 1; ++nk) {
                int nk1 = nk + 1;
                ratio2[nk1] = ratio2[nk];
                T frc1 = static_cast<T>(0.0);
                tu[nk1] = t0[nk1];
                theteu[nk1] = theteu[nk];
                qu[nk1] = qu[nk];
                qliq[nk1] = qliq[nk];
                qice[nk1] = qice[nk];
                T qnewlq, qnewic;
                tpmix2(*tables, p0[nk1], theteu[nk1], tu[nk1], qu[nk1], qliq[nk1],
                       qice[nk1], qnewlq, qnewic, xlv1, xlv0);

                if (tu[nk1] <= ttfrz) {
                    if (tu[nk1] > tbfrz) {
                        if (ttemp > ttfrz) ttemp = ttfrz;
                        frc1 = (ttemp - tu[nk1]) / (ttemp - tbfrz);
                    } else {
                        frc1 = static_cast<T>(1.0);
                    }
                    ttemp = tu[nk1];
                    T qfrz = (qliq[nk1] + qnewlq) * frc1;
                    qnewic = qnewic + qnewlq * frc1;
                    qnewlq = qnewlq - qnewlq * frc1;
                    qice[nk1] = qice[nk1] + qliq[nk1] * frc1;
                    qliq[nk1] = qliq[nk1] - qliq[nk1] * frc1;
                    dtfrznew(tu[nk1], p0[nk1], theteu[nk1], qu[nk1], qfrz, qice[nk1],
                             aliq, bliq, cliq, dliq);
                }
                tvu[nk1] = tu[nk1] * (static_cast<T>(1.0) + static_cast<T>(0.608) * qu[nk1]);

                T be, boterm, dzz;
                if (nk == k_lcl) {
                    be = (tvlcl + tvu[nk1]) / (tven + tv0[nk1]) - static_cast<T>(1.0);
                    boterm = static_cast<T>(2.0) * (z0[nk1] - zlcl) * g * be / static_cast<T>(1.5);
                    dzz = z0[nk1] - zlcl;
                } else {
                    be = (tvu[nk] + tvu[nk1]) / (tv0[nk] + tv0[nk1]) - static_cast<T>(1.0);
                    boterm = static_cast<T>(2.0) * dza[nk] * g * be / static_cast<T>(1.5);
                    dzz = dza[nk];
                }
                T enterm = static_cast<T>(2.0) * rei * wtw / upold;

                condload(qliq[nk1], qice[nk1], wtw, dzz, boterm, enterm, rate,
                         qnewlq, qnewic, qlqout[nk1], qicout[nk1], g);

                if (wtw < static_cast<T>(1.0e-3)) {
                    ltop = nk;
                    break;
                } else {
                    wu[nk1] = sqrt(wtw);
                }

                envirtht(*tables, p0[nk1], t0[nk1], q0[nk1], thetee[nk1], aliq, bliq, cliq, dliq);
                rei = vmflcl * dp[nk1] * static_cast<T>(0.03) / rad;
                tvqu[nk1] = tu[nk1] * (static_cast<T>(1.0) + static_cast<T>(0.608) * qu[nk1] - qliq[nk1] - qice[nk1]);
                if (nk == k_lcl) {
                    dilbe = ((tvlcl + tvqu[nk1]) / (tven + tv0[nk1]) - static_cast<T>(1.0)) * dzz;
                } else {
                    dilbe = ((tvqu[nk] + tvqu[nk1]) / (tv0[nk] + tv0[nk1]) - static_cast<T>(1.0)) * dzz;
                }
                if (dilbe > static_cast<T>(0.0)) abe += dilbe * g;

                T ee2, ud2;
                if (tvqu[nk1] <= tv0[nk1]) {
                    ee2 = static_cast<T>(0.5);
                    ud2 = static_cast<T>(1.0);
                    eqfrc[nk1] = static_cast<T>(0.0);
                } else {
                    let = nk1;
                    T ttmp = tvqu[nk1];
                    T f1 = static_cast<T>(0.95);
                    T f2 = static_cast<T>(1.0) - f1;
                    T thttmp = f1 * thetee[nk1] + f2 * theteu[nk1];
                    T qtmp = f1 * q0[nk1] + f2 * qu[nk1];
                    T tmpliq = f2 * qliq[nk1];
                    T tmpice = f2 * qice[nk1];
                    T dummy1, dummy2;
                    tpmix2(*tables, p0[nk1], thttmp, ttmp, qtmp, tmpliq, tmpice, dummy1, dummy2, xlv1, xlv0);
                    T tu95 = ttmp * (static_cast<T>(1.0) + static_cast<T>(0.608) * qtmp - tmpliq - tmpice);
                    if (tu95 > tv0[nk1]) {
                        ee2 = static_cast<T>(1.0);
                        ud2 = static_cast<T>(0.0);
                        eqfrc[nk1] = static_cast<T>(1.0);
                    } else {
                        f1 = static_cast<T>(0.10);
                        f2 = static_cast<T>(1.0) - f1;
                        thttmp = f1 * thetee[nk1] + f2 * theteu[nk1];
                        qtmp = f1 * q0[nk1] + f2 * qu[nk1];
                        tmpliq = f2 * qliq[nk1];
                        tmpice = f2 * qice[nk1];
                        tpmix2(*tables, p0[nk1], thttmp, ttmp, qtmp, tmpliq, tmpice, dummy1, dummy2, xlv1, xlv0);
                        T tu10 = ttmp * (static_cast<T>(1.0) + static_cast<T>(0.608) * qtmp - tmpliq - tmpice);
                        T tvdiff = abs(tu10 - tvqu[nk1]);
                        if (tvdiff < static_cast<T>(1.0e-3)) {
                            ee2 = static_cast<T>(1.0);
                            ud2 = static_cast<T>(0.0);
                            eqfrc[nk1] = static_cast<T>(1.0);
                        } else {
                            eqfrc[nk1] = (tv0[nk1] - tvqu[nk1]) * f1 / (tu10 - tvqu[nk1]);
                            if (eqfrc[nk1] < static_cast<T>(0.0)) eqfrc[nk1] = static_cast<T>(0.0);
                            if (eqfrc[nk1] > static_cast<T>(1.0)) eqfrc[nk1] = static_cast<T>(1.0);
                            if (eqfrc[nk1] == static_cast<T>(1.0)) {
                                ee2 = static_cast<T>(1.0);
                                ud2 = static_cast<T>(0.0);
                            } else if (eqfrc[nk1] == static_cast<T>(0.0)) {
                                ee2 = static_cast<T>(0.0);
                                ud2 = static_cast<T>(1.0);
                            } else {
                                prof5(eqfrc[nk1], ee2, ud2);
                            }
                        }
                    }
                }

                if (ee2 < static_cast<T>(0.5)) ee2 = static_cast<T>(0.5);
                ud2 = static_cast<T>(1.5) * ud2;
                uer[nk1] = static_cast<T>(0.5) * rei * (ee1 + ee2);
                udr[nk1] = static_cast<T>(0.5) * rei * (ud1 + ud2);

                if (umf[nk] - udr[nk1] < static_cast<T>(10.0)) {
                    if (dilbe > static_cast<T>(0.0)) {
                        abe -= dilbe * g;
                    }
                    let = nk;
                    ltop = nk;
                    break;
                } else {
                    ee1 = ee2;
                    ud1 = ud2;
                    upold = umf[nk] - udr[nk1];
                    upnew = upold + uer[nk1];
                    umf[nk1] = upnew;
                    dilfrc[nk1] = upnew / upold;
                    detlq[nk1] = qliq[nk1] * udr[nk1];
                    detic[nk1] = qice[nk1] * udr[nk1];
                    qdt[nk1] = qu[nk1];
                    qu[nk1] = (upold * qu[nk1] + uer[nk1] * q0[nk1]) / upnew;
                    theteu[nk1] = (theteu[nk1] * upold + thetee[nk1] * uer[nk1]) / upnew;
                    qliq[nk1] = qliq[nk1] * upold / upnew;
                    qice[nk1] = qice[nk1] * upold / upnew;
                    pptliq[nk1] = qlqout[nk1] * umf[nk];
                    pptice[nk1] = qicout[nk1] * umf[nk];
                    trppt += pptliq[nk1] + pptice[nk1];
                    if (nk1 <= kpbl) uer[nk1] += vmflcl * dp[nk1] / dpthmx;
                }
            } // updraft loop

            cldhgt[lc] = z0[ltop] - zlcl;
            T chmin;
            if (tlcl > static_cast<T>(293.0)) {
                chmin = static_cast<T>(4.0e3);
            } else if (tlcl >= static_cast<T>(273.0)) {
                chmin = static_cast<T>(2.0e3) + static_cast<T>(100.0) * (tlcl - static_cast<T>(273.0));
            } else {
                chmin = static_cast<T>(2.0e3);
            }

            for (int nk = k_lcl; nk <= ltop; ++nk) {
                qc_kf(i, j, nk) = qliq[nk];
                qi_kf(i, j, nk) = qice[nk];
            }

            if (ltop <= klcl || ltop <= kpbl || let + 1 <= kpbl) {
                cldhgt[lc] = static_cast<T>(0.0);
                for (int nk = k_lcl; nk <= ltop; ++nk) {
                    umf[nk] = static_cast<T>(0.0);
                    udr[nk] = static_cast<T>(0.0);
                    uer[nk] = static_cast<T>(0.0);
                    detlq[nk] = static_cast<T>(0.0);
                    detic[nk] = static_cast<T>(0.0);
                    pptliq[nk] = static_cast<T>(0.0);
                    pptice[nk] = static_cast<T>(0.0);
                    cldfra_dp_kf(i, j, nk) = static_cast<T>(0.0);
                    cldfra_sh_kf(i, j, nk) = static_cast<T>(0.0);
                    qc_kf(i, j, nk) = static_cast<T>(0.0);
                    qi_kf(i, j, nk) = static_cast<T>(0.0);
                }
            } else if (cldhgt[lc] > chmin && abe > static_cast<T>(1.0)) {
                ishall = 0;
                for (int nk = k_lcl; nk <= ltop; ++nk) {
                    cldfra_sh_kf(i, j, nk) = static_cast<T>(0.0);
                }
                break; // Deep convection found
            } else {
                ishall = 1;
                for (int nk = k_lcl; nk <= ltop; ++nk) {
                    cldfra_dp_kf(i, j, nk) = static_cast<T>(0.0);
                }
                if (nu == nuchm) {
                    break; // Shallow convection found
                } else {
                    for (int nk = k_lcl; nk <= ltop; ++nk) {
                        umf[nk] = static_cast<T>(0.0);
                        udr[nk] = static_cast<T>(0.0);
                        uer[nk] = static_cast<T>(0.0);
                        detlq[nk] = static_cast<T>(0.0);
                        detic[nk] = static_cast<T>(0.0);
                        pptliq[nk] = static_cast<T>(0.0);
                        pptice[nk] = static_cast<T>(0.0);
                        cldfra_dp_kf(i, j, nk) = static_cast<T>(0.0);
                        cldfra_sh_kf(i, j, nk) = static_cast<T>(0.0);
                        qc_kf(i, j, nk) = static_cast<T>(0.0);
                        qi_kf(i, j, nk) = static_cast<T>(0.0);
                    }
                }
            }
        } // end usl while loop

        if (ishall == 1) {
            int kstart = (kpbl > klcl) ? kpbl : klcl;
            let = kstart;
        }

        if (let == ltop) {
            udr[ltop] = umf[ltop] + udr[ltop] - uer[ltop];
            detlq[ltop] = qliq[ltop] * udr[ltop] * upnew / upold;
            detic[ltop] = qice[ltop] * udr[ltop] * upnew / upold;
            uer[ltop] = static_cast<T>(0.0);
            umf[ltop] = static_cast<T>(0.0);
        } else {
            T dptt = static_cast<T>(0.0);
            for (int nj = let + 1; nj <= ltop; ++nj) {
                dptt += dp[nj];
            }
            T dumfdp = umf[let] / dptt;
            for (int nk = let + 1; nk <= ltop; ++nk) {
                if (nk == ltop) {
                    udr[nk] = umf[nk - 1];
                    uer[nk] = static_cast<T>(0.0);
                    detlq[nk] = udr[nk] * qliq[nk] * dilfrc[nk];
                    detic[nk] = udr[nk] * qice[nk] * dilfrc[nk];
                } else {
                    umf[nk] = umf[nk - 1] - dp[nk] * dumfdp;
                    uer[nk] = umf[nk] * (static_cast<T>(1.0) - static_cast<T>(1.0) / dilfrc[nk]);
                    udr[nk] = umf[nk - 1] - umf[nk] + uer[nk];
                    detlq[nk] = udr[nk] * qliq[nk] * dilfrc[nk];
                    detic[nk] = udr[nk] * qice[nk] * dilfrc[nk];
                }
                if (nk >= let + 2) {
                    trppt = trppt - pptliq[nk] - pptice[nk];
                    pptliq[nk] = umf[nk - 1] * qlqout[nk];
                    pptice[nk] = umf[nk - 1] * qicout[nk];
                    trppt = trppt + pptliq[nk] + pptice[nk];
                }
            }
        }

        int ml = 0;
        for (int nk = 0; nk <= ltop; ++nk) {
            if (t0[nk] > t00) ml = nk;
        }

        int k_lcl = klcl - 1;
        for (int nk = 0; nk <= k_lcl; ++nk) {
            if (nk >= lc) {
                if (nk == lc) {
                    umf[nk] = vmflcl * dp[nk] / dpthmx;
                    uer[nk] = vmflcl * dp[nk] / dpthmx;
                } else if (nk <= kpbl) {
                    uer[nk] = vmflcl * dp[nk] / dpthmx;
                    umf[nk] = umf[nk - 1] + uer[nk];
                } else {
                    umf[nk] = vmflcl;
                    uer[nk] = static_cast<T>(0.0);
                }
                tu[nk] = tmix + (z0[nk] - zmix) * gdry;
                qu[nk] = qmix;
                wu[nk] = wlcl;
            } else {
                tu[nk] = static_cast<T>(0.0);
                qu[nk] = static_cast<T>(0.0);
                umf[nk] = static_cast<T>(0.0);
                wu[nk] = static_cast<T>(0.0);
                uer[nk] = static_cast<T>(0.0);
                cldfra_dp_kf(i, j, nk) = static_cast<T>(0.0);
                cldfra_sh_kf(i, j, nk) = static_cast<T>(0.0);
                qc_kf(i, j, nk) = static_cast<T>(0.0);
                qi_kf(i, j, nk) = static_cast<T>(0.0);
            }
            udr[nk] = static_cast<T>(0.0);
            qdt[nk] = static_cast<T>(0.0);
            qliq[nk] = static_cast<T>(0.0);
            qice[nk] = static_cast<T>(0.0);
            qlqout[nk] = static_cast<T>(0.0);
            qicout[nk] = static_cast<T>(0.0);
            pptliq[nk] = static_cast<T>(0.0);
            pptice[nk] = static_cast<T>(0.0);
            detlq[nk] = static_cast<T>(0.0);
            detic[nk] = static_cast<T>(0.0);
            ratio2[nk] = static_cast<T>(0.0);
            envirtht(*tables, p0[nk], t0[nk], q0[nk], thetee[nk], aliq, bliq, cliq, dliq);
            eqfrc[nk] = static_cast<T>(1.0);
        }

        int ltop1 = ltop + 1;
        int ltopm1 = ltop - 1;
        for (int nk = ltop1; nk < kx; ++nk) {
            umf[nk] = static_cast<T>(0.0);
            udr[nk] = static_cast<T>(0.0);
            uer[nk] = static_cast<T>(0.0);
            qdt[nk] = static_cast<T>(0.0);
            qliq[nk] = static_cast<T>(0.0);
            qice[nk] = static_cast<T>(0.0);
            qlqout[nk] = static_cast<T>(0.0);
            qicout[nk] = static_cast<T>(0.0);
            detlq[nk] = static_cast<T>(0.0);
            detic[nk] = static_cast<T>(0.0);
            pptliq[nk] = static_cast<T>(0.0);
            pptice[nk] = static_cast<T>(0.0);
            if (nk > ltop1) {
                tu[nk] = static_cast<T>(0.0);
                qu[nk] = static_cast<T>(0.0);
                wu[nk] = static_cast<T>(0.0);
                cldfra_dp_kf(i, j, nk) = static_cast<T>(0.0);
                cldfra_sh_kf(i, j, nk) = static_cast<T>(0.0);
                qc_kf(i, j, nk) = static_cast<T>(0.0);
                qi_kf(i, j, nk) = static_cast<T>(0.0);
            }
            thta0[nk] = static_cast<T>(0.0);
            thtau[nk] = static_cast<T>(0.0);
            ems[nk] = static_cast<T>(0.0);
            emsd[nk] = static_cast<T>(0.0);
            tg[nk] = t0[nk];
            qg[nk] = q0[nk];
            qlg[nk] = static_cast<T>(0.0);
            qig[nk] = static_cast<T>(0.0);
            qrg[nk] = static_cast<T>(0.0);
            qsg[nk] = static_cast<T>(0.0);
            omg[nk] = static_cast<T>(0.0);
        }
        omg[kx] = static_cast<T>(0.0);

        for (int nk = 0; nk <= ltop; ++nk) {
            ems[nk] = dp[nk] * dxsq / g;
            emsd[nk] = static_cast<T>(1.0) / ems[nk];
            exn[nk] = pow(p00 / p0[nk], static_cast<T>(0.2854) * (static_cast<T>(1.0) - static_cast<T>(0.28) * qdt[nk]));
            thtau[nk] = tu[nk] * exn[nk];
            exn[nk] = pow(p00 / p0[nk], static_cast<T>(0.2854) * (static_cast<T>(1.0) - static_cast<T>(0.28) * q0[nk]));
            thta0[nk] = t0[nk] * exn[nk];
            ddilfrc[nk] = static_cast<T>(1.0) / dilfrc[nk];
            omg[nk] = static_cast<T>(0.0);
        }

        // Convective timescale
        wspd[klcl] = sqrt(u0[klcl] * u0[klcl] + v0[klcl] * v0[klcl]);
        wspd[l5] = sqrt(u0[l5] * u0[l5] + v0[l5] * v0[l5]);
        wspd[ltop] = sqrt(u0[ltop] * u0[ltop] + v0[ltop] * v0[ltop]);
        T vconv = static_cast<T>(0.5) * (wspd[klcl] + wspd[l5]);
        T timec = dx / vconv;
        T tadvec = timec;
        if (timec < static_cast<T>(1800.0)) timec = static_cast<T>(1800.0);
        if (timec > static_cast<T>(3600.0)) timec = static_cast<T>(3600.0);
        if (ishall == 1) timec = static_cast<T>(2400.0);
        int nic = static_cast<int>(std::round(timec / dt));
        timec = static_cast<T>(nic) * dt;

        T shsign = (wspd[ltop] > wspd[klcl]) ? static_cast<T>(1.0) : static_cast<T>(-1.0);
        T vws = (u0[ltop] - u0[klcl]) * (u0[ltop] - u0[klcl]) + (v0[ltop] - v0[klcl]) * (v0[ltop] - v0[klcl]);
        vws = static_cast<T>(1.0e3) * shsign * sqrt(vws) / (z0[ltop] - z0[lcl]);
        T pef = static_cast<T>(1.591) + vws * (static_cast<T>(-0.639) + vws * (static_cast<T>(9.53e-2) - vws * static_cast<T>(4.96e-3)));
        if (pef < static_cast<T>(0.2)) pef = static_cast<T>(0.2);
        if (pef > static_cast<T>(0.9)) pef = static_cast<T>(0.9);

        T cbh = (zlcl - z0[0]) * static_cast<T>(3.281e-3);
        T rcbh;
        if (cbh < static_cast<T>(3.0)) {
            rcbh = static_cast<T>(0.02);
        } else {
            rcbh = static_cast<T>(0.96729352) + cbh * (static_cast<T>(-0.70034167) + cbh * (static_cast<T>(0.162179896) +
                   cbh * (static_cast<T>(-1.2569798e-2) + cbh * (static_cast<T>(4.2772e-4) - cbh * static_cast<T>(5.44e-6)))));
        }
        if (cbh > static_cast<T>(25.0)) rcbh = static_cast<T>(2.4);
        T pefcbh = static_cast<T>(1.0) / (static_cast<T>(1.0) + rcbh);
        if (pefcbh > static_cast<T>(0.9)) pefcbh = static_cast<T>(0.9);
        T peff = static_cast<T>(0.5) * (pef + pefcbh);

        // Downdraft properties
        T tder = static_cast<T>(0.0);
        int lfs = 0;
        int ldt = 0;
        int ldb = 0;
        int kstart = kpbl + 1;
        T dmffrc = static_cast<T>(0.0);

        if (ishall == 1) {
            lfs = 0;
        } else {
            int klfs = let - 1;
            for (int nk = kstart + 1; nk < kl; ++nk) {
                T dppp = p0[kstart] - p0[nk];
                if (dppp > static_cast<T>(150.0e2)) {
                    klfs = nk;
                    break;
                }
            }
            if (klfs > let - 1) klfs = let - 1;
            lfs = klfs;

            if (p0[kstart] - p0[lfs] > static_cast<T>(50.0e2)) {
                theted[lfs] = thetee[lfs];
                qd[lfs] = q0[lfs];
                T qss;
                tpmix2dd(*tables, p0[lfs], theted[lfs], tz[lfs], qss);
                thtad[lfs] = tz[lfs] * pow(p00 / p0[lfs], static_cast<T>(0.2854) * (static_cast<T>(1.0) - static_cast<T>(0.28) * qss));
                tvd[lfs] = tz[lfs] * (static_cast<T>(1.0) + static_cast<T>(0.608) * qss);
                T rdd = p0[lfs] / (r * tvd[lfs]);
                T a1 = (static_cast<T>(1.0) - peff) * au0;
                dmf[lfs] = -a1 * rdd;
                der[lfs] = dmf[lfs];
                ddr[lfs] = static_cast<T>(0.0);
                T rhbar = rh[lfs] * dp[lfs];
                T dptt = dp[lfs];
                for (int nd = lfs - 1; nd >= kstart; --nd) {
                    int nd1 = nd + 1;
                    der[nd] = der[lfs] * ems[nd] / ems[lfs];
                    ddr[nd] = static_cast<T>(0.0);
                    dmf[nd] = dmf[nd1] + der[nd];
                    theted[nd] = (theted[nd1] * dmf[nd1] + thetee[nd] * der[nd]) / dmf[nd];
                    qd[nd] = (qd[nd1] * dmf[nd1] + q0[nd] * der[nd]) / dmf[nd];
                    dptt += dp[nd];
                    rhbar += rh[nd] * dp[nd];
                }
                rhbar /= dptt;
                dmffrc = static_cast<T>(2.0) * (static_cast<T>(1.0) - rhbar);
                T dpdd = static_cast<T>(0.0);
                T pptmlt = static_cast<T>(0.0);
                for (int nk = klcl; nk <= ltop; ++nk) {
                    pptmlt += pptice[nk];
                }
                T dtmelt = static_cast<T>(0.0);
                if (lc < ml) {
                    dtmelt = rlf * pptmlt / (cp * umf[klcl]);
                }
                ldt = (lfs - 1 < kstart - 1) ? (lfs - 1) : (kstart - 1);
                tpmix2dd(*tables, p0[kstart], theted[kstart], tz[kstart], qss);
                tz[kstart] -= dtmelt;
                T es = aliq * exp((bliq * tz[kstart] - cliq) / (tz[kstart] - dliq));
                qss = static_cast<T>(0.622) * es / (p0[kstart] - es);
                theted[kstart] = tz[kstart] * pow(static_cast<T>(1.0e5) / p0[kstart], static_cast<T>(0.2854) * (static_cast<T>(1.0) - static_cast<T>(0.28) * qss)) *
                                exp((static_cast<T>(3374.6525) / tz[kstart] - static_cast<T>(2.5403)) * qss * (static_cast<T>(1.0) + static_cast<T>(0.81) * qss));

                ldt = (lfs - 1 < kstart - 1) ? (lfs - 1) : (kstart - 1);
                ldb = 0;
                for (int nd = ldt; nd >= 0; --nd) {
                    dpdd += dp[nd];
                    theted[nd] = theted[kstart];
                    qd[nd] = qd[kstart];
                    tpmix2dd(*tables, p0[nd], theted[nd], tz[nd], qss);
                    qsd[nd] = qss;
                    T rhh = static_cast<T>(1.0) - static_cast<T>(0.2) / static_cast<T>(1000.0) * (z0[kstart] - z0[nd]);
                    if (rhh < static_cast<T>(1.0)) {
                        T dssdt = (cliq - bliq * dliq) / ((tz[nd] - dliq) * (tz[nd] - dliq));
                        T rl = xlv0 - xlv1 * tz[nd];
                        T dtmp = rl * qss * (static_cast<T>(1.0) - rhh) / (cp + rl * rhh * qss * dssdt);
                        T t1rh = tz[nd] + dtmp;
                        es = rhh * aliq * exp((bliq * t1rh - cliq) / (t1rh - dliq));
                        T qsrh = static_cast<T>(0.622) * es / (p0[nd] - es);
                        if (qsrh < qd[nd]) {
                            qsrh = qd[nd];
                            t1rh = tz[nd] + (qss - qsrh) * rl / cp;
                        }
                        tz[nd] = t1rh;
                        qss = qsrh;
                        qsd[nd] = qss;
                    }
                    tvd[nd] = tz[nd] * (static_cast<T>(1.0) + static_cast<T>(0.608) * qsd[nd]);
                    if (tvd[nd] > tv0[nd] || nd == 0) {
                        ldb = nd;
                        break;
                    }
                }
                if (p0[ldb] - p0[lfs] > static_cast<T>(50.0e2)) {
                    for (int nd = ldt; nd >= ldb; --nd) {
                        int nd1 = nd + 1;
                        ddr[nd] = -dmf[kstart] * dp[nd] / dpdd;
                        der[nd] = static_cast<T>(0.0);
                        dmf[nd] = dmf[nd1] + ddr[nd];
                        tder += (qsd[nd] - qd[nd]) * ddr[nd];
                        qd[nd] = qsd[nd];
                        thtad[nd] = tz[nd] * pow(p00 / p0[nd], static_cast<T>(0.2854) * (static_cast<T>(1.0) - static_cast<T>(0.28) * qd[nd]));
                    }
                }
            }
        }

        T cpr = trppt;
        T pptflx = trppt;
        if (tder < static_cast<T>(1.0)) {
            pptflx = trppt;
            cpr = trppt;
            tder = static_cast<T>(0.0);
            ldb = lfs;
            for (int ndk = 0; ndk <= ltop; ++ndk) {
                dmf[ndk] = static_cast<T>(0.0);
                der[ndk] = static_cast<T>(0.0);
                ddr[ndk] = static_cast<T>(0.0);
                thtad[ndk] = static_cast<T>(0.0);
                wd[ndk] = static_cast<T>(0.0);
                tz[ndk] = static_cast<T>(0.0);
                qd[ndk] = static_cast<T>(0.0);
            }
        } else {
            T ddinc = -dmffrc * umf[klcl] / dmf[kstart];
            if (tder * ddinc > trppt) {
                ddinc = trppt / tder;
            }
            tder *= ddinc;
            for (int nk = ldb; nk <= lfs; ++nk) {
                dmf[nk] *= ddinc;
                der[nk] *= ddinc;
                ddr[nk] *= ddinc;
            }
            cpr = trppt;
            pptflx = trppt - tder;
            if (pptflx < static_cast<T>(0.0)) pptflx = static_cast<T>(0.0);
            peff = pptflx / trppt;

            if (ldb > 0) {
                for (int nk = 0; nk < ldb; ++nk) {
                    dmf[nk] = static_cast<T>(0.0);
                    der[nk] = static_cast<T>(0.0);
                    ddr[nk] = static_cast<T>(0.0);
                    wd[nk] = static_cast<T>(0.0);
                    tz[nk] = static_cast<T>(0.0);
                    qd[nk] = static_cast<T>(0.0);
                    thtad[nk] = static_cast<T>(0.0);
                }
            }
            for (int nk = lfs + 1; nk < kx; ++nk) {
                dmf[nk] = static_cast<T>(0.0);
                der[nk] = static_cast<T>(0.0);
                ddr[nk] = static_cast<T>(0.0);
                wd[nk] = static_cast<T>(0.0);
                tz[nk] = static_cast<T>(0.0);
                qd[nk] = static_cast<T>(0.0);
                thtad[nk] = static_cast<T>(0.0);
            }
            for (int nk = ldt + 1; nk < lfs; ++nk) {
                tz[nk] = static_cast<T>(0.0);
                qd[nk] = static_cast<T>(0.0);
                thtad[nk] = static_cast<T>(0.0);
            }
        }

        T aincmx = static_cast<T>(1000.0);
        int lmax = (klcl > lfs) ? klcl : lfs;
        for (int nk = lc; nk <= lmax; ++nk) {
            if ((uer[nk] - der[nk]) > static_cast<T>(1.0e-3)) {
                T aincm1 = ems[nk] / ((uer[nk] - der[nk]) * timec);
                if (aincm1 < aincmx) aincmx = aincm1;
            }
        }
        T ainc = static_cast<T>(1.0);
        if (aincmx < ainc) ainc = aincmx;

        T tder2 = tder;
        T pptfl2 = pptflx;
        for (int nk = 0; nk <= ltop; ++nk) {
            detlq2[nk] = detlq[nk];
            detic2[nk] = detic[nk];
            udr2[nk] = udr[nk];
            uer2[nk] = uer[nk];
            ddr2[nk] = ddr[nk];
            der2[nk] = der[nk];
            umf2[nk] = umf[nk];
            dmf2[nk] = dmf[nk];
        }
        T fabe = static_cast<T>(1.0);
        const T stab = static_cast<T>(0.95);
        int noitr = 0;

        if (ishall == 1) {
            T tkemax = static_cast<T>(5.0);
            T evac = static_cast<T>(0.5) * tkemax * static_cast<T>(0.1);
            ainc = evac * dpthmx * dxsq / (vmflcl * g * timec);
            tder = tder2 * ainc;
            pptflx = pptfl2 * ainc;
            for (int nk = 0; nk <= ltop; ++nk) {
                umf[nk] = umf2[nk] * ainc;
                dmf[nk] = dmf2[nk] * ainc;
                detlq[nk] = detlq2[nk] * ainc;
                detic[nk] = detic2[nk] * ainc;
                udr[nk] = udr2[nk] * ainc;
                uer[nk] = uer2[nk] * ainc;
                der[nk] = der2[nk] * ainc;
                ddr[nk] = ddr2[nk] * ainc;
            }
        }

        // Iterative closure loop
        T aincold = static_cast<T>(0.0);
        T fabeold = static_cast<T>(0.0);
        int nstep = 0;
        T dtime = static_cast<T>(0.0);

        for (int ncount = 1; ncount <= 10; ++ncount) {
            T dtt = timec;
            for (int nk = 0; nk <= ltop; ++nk) {
                domgdp[nk] = -(uer[nk] - der[nk] - udr[nk] - ddr[nk]) * emsd[nk];
                if (nk > 0) {
                    omg[nk] = omg[nk - 1] - dp[nk - 1] * domgdp[nk - 1];
                    T absomg = abs(omg[nk]);
                    T absomgtc = absomg * timec;
                    T frdp = static_cast<T>(0.75) * dp[nk - 1];
                    if (absomgtc > frdp) {
                        T dtt1 = frdp / absomg;
                        if (dtt1 < dtt) dtt = dtt1;
                    }
                }
            }
            nstep = static_cast<int>(std::round(timec / dtt + static_cast<T>(1.0)));
            dtime = timec / static_cast<T>(nstep);
            for (int nk = 0; nk <= ltop; ++nk) {
                thpa[nk] = thta0[nk];
                qpa[nk] = q0[nk];
                fxm[nk] = omg[nk] * dxsq / g;
            }

            for (int ntc = 1; ntc <= nstep; ++ntc) {
                for (int nk = 0; nk <= ltop; ++nk) {
                    thfxin[nk] = static_cast<T>(0.0);
                    thfxout[nk] = static_cast<T>(0.0);
                    qfxin[nk] = static_cast<T>(0.0);
                    qfxout[nk] = static_cast<T>(0.0);
                }
                for (int nk = 1; nk <= ltop; ++nk) {
                    if (omg[nk] <= static_cast<T>(0.0)) {
                        thfxin[nk] = -fxm[nk] * thpa[nk - 1];
                        qfxin[nk] = -fxm[nk] * qpa[nk - 1];
                        thfxout[nk - 1] += thfxin[nk];
                        qfxout[nk - 1] += qfxin[nk];
                    } else {
                        thfxout[nk] = fxm[nk] * thpa[nk];
                        qfxout[nk] = fxm[nk] * qpa[nk];
                        thfxin[nk - 1] += thfxout[nk];
                        qfxin[nk - 1] += qfxout[nk];
                    }
                }
                for (int nk = 0; nk <= ltop; ++nk) {
                    thpa[nk] += (thfxin[nk] + udr[nk] * thtau[nk] + ddr[nk] * thtad[nk] -
                                 thfxout[nk] - (uer[nk] - der[nk]) * thta0[nk]) * dtime * emsd[nk];
                    qpa[nk] += (qfxin[nk] + udr[nk] * qdt[nk] + ddr[nk] * qd[nk] -
                                qfxout[nk] - (uer[nk] - der[nk]) * q0[nk]) * dtime * emsd[nk];
                }
            }
            for (int nk = 0; nk <= ltop; ++nk) {
                thtag[nk] = thpa[nk];
                qg[nk] = qpa[nk];
            }

            // Check non-negative moisture
            for (int nk = 0; nk <= ltop; ++nk) {
                if (qg[nk] < static_cast<T>(0.0)) {
                    if (nk == 0) {
                        qg[nk] = static_cast<T>(1.0e-9);
                        break;
                    }
                    int nk1 = (nk == ltop) ? klcl : (nk + 1);
                    T tma = qg[nk1] * ems[nk1];
                    T tmb = qg[nk - 1] * ems[nk - 1];
                    T tmm = (qg[nk] - static_cast<T>(1.0e-9)) * ems[nk];
                    T bcoeff = -tmm / ((tma * tma) / tmb + tmb);
                    T acoeff = bcoeff * tma / tmb;
                    tmb *= (static_cast<T>(1.0) - bcoeff);
                    tma *= (static_cast<T>(1.0) - acoeff);
                    qg[nk] = static_cast<T>(1.0e-9);
                    qg[nk1] = tma * emsd[nk1];
                    qg[nk - 1] = tmb * emsd[nk - 1];
                }
            }

            T topomg = (udr[ltop] - uer[ltop]) * dp[ltop] * emsd[ltop];
            if (abs(topomg - omg[ltop]) > static_cast<T>(1.0e-3)) {
                break;
            }

            for (int nk = 0; nk <= ltop; ++nk) {
                exn[nk] = pow(p00 / p0[nk], static_cast<T>(0.2854) * (static_cast<T>(1.0) - static_cast<T>(0.28) * qg[nk]));
                tg[nk] = thtag[nk] / exn[nk];
                tvg[nk] = tg[nk] * (static_cast<T>(1.0) + static_cast<T>(0.608) * qg[nk]);
            }
            if (ishall == 1) {
                break;
            }

            // Adjusted ABE
            tmix = static_cast<T>(0.0);
            qmix = static_cast<T>(0.0);
            for (int nk = lc; nk <= kpbl; ++nk) {
                tmix += dp[nk] * tg[nk];
                qmix += dp[nk] * qg[nk];
            }
            tmix /= dpthmx;
            qmix /= dpthmx;
            T es = aliq * exp((tmix * bliq - cliq) / (tmix - dliq));
            T qss = static_cast<T>(0.622) * es / (pmix - es);
            if (qmix > qss) {
                T rl = xlv0 - xlv1 * tmix;
                T cpm = cp * (static_cast<T>(1.0) + static_cast<T>(0.887) * qmix);
                T dssdt = qss * (cliq - bliq * dliq) / ((tmix - dliq) * (tmix - dliq));
                T dq = (qmix - qss) / (static_cast<T>(1.0) + rl * dssdt / cpm);
                tmix += rl / cp * dq;
                qmix -= dq;
                tlcl = tmix;
            } else {
                if (qmix < static_cast<T>(0.0)) qmix = static_cast<T>(0.0);
                T emix = qmix * pmix / (static_cast<T>(0.622) + qmix);
                const T astrt = static_cast<T>(1.0e-3);
                const T binc = static_cast<T>(0.075);
                T a1 = emix / aliq;
                T tp = (a1 - astrt) / binc;
                int indlu = static_cast<int>(tp) + 1;
                if (indlu < 1) indlu = 1;
                if (indlu > 199) indlu = 199;
                T value = static_cast<T>(indlu - 1) * binc + astrt;
                T aintrp = (a1 - value) / binc;
                T tlog = aintrp * tables->alu[indlu] + (static_cast<T>(1.0) - aintrp) * tables->alu[indlu - 1];
                T tdpt = (cliq - dliq * tlog) / (bliq - tlog);
                tlcl = tdpt - (static_cast<T>(0.212) + static_cast<T>(1.571e-3) * (tdpt - t00) -
                               static_cast<T>(4.36e-4) * (tmix - t00)) * (tmix - tdpt);
                if (tlcl > tmix) tlcl = tmix;
            }
            tvlcl = tlcl * (static_cast<T>(1.0) + static_cast<T>(0.608) * qmix);
            zlcl = zmix + (tlcl - tmix) / gdry;
            for (int nk = lc; nk < kl; ++nk) {
                klcl = nk;
                if (zlcl <= z0[nk]) break;
            }
            int k_lcl_adj = klcl - 1;
            T dlp = (zlcl - z0[k_lcl_adj]) / (z0[klcl] - z0[k_lcl_adj]);
            T tenv = tg[k_lcl_adj] + (tg[klcl] - tg[k_lcl_adj]) * dlp;
            T qenv = qg[k_lcl_adj] + (qg[klcl] - qg[k_lcl_adj]) * dlp;
            T tven = tenv * (static_cast<T>(1.0) + static_cast<T>(0.608) * qenv);
            theteu[k_lcl_adj] = tmix * pow(static_cast<T>(1.0e5) / pmix, static_cast<T>(0.2854) * (static_cast<T>(1.0) - static_cast<T>(0.28) * qmix)) *
                               exp((static_cast<T>(3374.6525) / tlcl - static_cast<T>(2.5403)) * qmix * (static_cast<T>(1.0) + static_cast<T>(0.81) * qmix));

            T abeg = static_cast<T>(0.0);
            for (int nk = k_lcl_adj; nk <= ltopm1; ++nk) {
                int nk1 = nk + 1;
                theteu[nk1] = theteu[nk];
                tpmix2dd(*tables, p0[nk1], theteu[nk1], tgu[nk1], qgu[nk1]);
                tvqu[nk1] = tgu[nk1] * (static_cast<T>(1.0) + static_cast<T>(0.608) * qgu[nk1] - qliq[nk1] - qice[nk1]);
                T dzz = (nk == k_lcl_adj) ? (z0[klcl] - zlcl) : dza[nk];
                T dilbe = ((nk == k_lcl_adj) ? ((tvlcl + tvqu[nk1]) / (tven + tvg[nk1]) - static_cast<T>(1.0))
                                             : ((tvqu[nk] + tvqu[nk1]) / (tvg[nk] + tvg[nk1]) - static_cast<T>(1.0))) * dzz;
                if (dilbe > static_cast<T>(0.0)) abeg += dilbe * g;
                envirtht(*tables, p0[nk1], tg[nk1], qg[nk1], thteeg[nk1], aliq, bliq, cliq, dliq);
                theteu[nk1] = theteu[nk1] * ddilfrc[nk1] + thteeg[nk1] * (static_cast<T>(1.0) - ddilfrc[nk1]);
            }

            if (noitr == 1) {
                break;
            }
            T dabe = (abe - abeg > static_cast<T>(0.1) * abe) ? (abe - abeg) : (static_cast<T>(0.1) * abe);
            fabe = abeg / abe;
            if (fabe > static_cast<T>(1.0) && ishall == 0) {
                return;
            }
            if (ncount != 1) {
                if (abs(ainc - aincold) < static_cast<T>(0.0001)) {
                    noitr = 1;
                    ainc = aincold;
                    continue;
                }
                T dfda = (fabe - fabeold) / (ainc - aincold);
                if (dfda > static_cast<T>(0.0)) {
                    noitr = 1;
                    ainc = aincold;
                    continue;
                }
            }
            aincold = ainc;
            fabeold = fabe;
            if (ainc / aincmx > static_cast<T>(0.999) && fabe > static_cast<T>(1.05) - stab) {
                break;
            }
            if ((fabe <= static_cast<T>(1.05) - stab && fabe >= static_cast<T>(0.95) - stab) || ncount == 10) {
                break;
            } else {
                if (fabe == static_cast<T>(0.0)) {
                    ainc *= static_cast<T>(0.5);
                } else {
                    if (dabe < static_cast<T>(1.0e-4)) {
                        noitr = 1;
                        ainc = aincold;
                        continue;
                    } else {
                        ainc = ainc * stab * abe / dabe;
                    }
                }
                if (ainc > aincmx) ainc = aincmx;
                if (ainc < static_cast<T>(0.05)) {
                    return;
                }
                tder = tder2 * ainc;
                pptflx = pptfl2 * ainc;
                for (int nk = 0; nk <= ltop; ++nk) {
                    umf[nk] = umf2[nk] * ainc;
                    dmf[nk] = dmf2[nk] * ainc;
                    detlq[nk] = detlq2[nk] * ainc;
                    detic[nk] = detic2[nk] * ainc;
                    udr[nk] = udr2[nk] * ainc;
                    uer[nk] = uer2[nk] * ainc;
                    der[nk] = der2[nk] * ainc;
                    ddr[nk] = ddr2[nk] * ainc;
                }
            }
        } // iteration loop

        // Cloud fractions
        int k_start = (klcl - 1 < 0) ? 0 : (klcl - 1);
        int k_end = (ltop1 >= kl) ? (kl - 1) : ltop1;
        if (ishall == 1) {
            for (int nk = k_start; nk <= k_end; ++nk) {
                T umf_new = umf[nk] / dxsq;
                T xcldfra = static_cast<T>(0.07) * log(static_cast<T>(1.0) + (static_cast<T>(500.0) * umf_new));
                if (xcldfra < static_cast<T>(0.01)) xcldfra = static_cast<T>(0.01);
                if (xcldfra > static_cast<T>(0.2)) xcldfra = static_cast<T>(0.2);
                cldfra_sh_kf(i, j, nk) = xcldfra;
            }
        } else {
            for (int nk = k_start; nk <= k_end; ++nk) {
                T umf_new = umf[nk] / dxsq;
                T xcldfra = static_cast<T>(0.14) * log(static_cast<T>(1.0) + (static_cast<T>(500.0) * umf_new));
                if (xcldfra < static_cast<T>(0.01)) xcldfra = static_cast<T>(0.01);
                if (xcldfra > static_cast<T>(0.6)) xcldfra = static_cast<T>(0.6);
                cldfra_dp_kf(i, j, nk) = xcldfra;
            }
        }

        // Hydrometeor tendencies
        T frc2 = (cpr > static_cast<T>(0.0)) ? (pptflx / (cpr * ainc)) : static_cast<T>(0.0);
        for (int nk = 0; nk <= ltop; ++nk) {
            qlpa[nk] = ql0[nk];
            qipa[nk] = qi0[nk];
            qrpa[nk] = qr0[nk];
            qspa[nk] = qs0[nk];
            rainfb[nk] = pptliq[nk] * ainc * fbfrc * frc2;
            snowfb[nk] = pptice[nk] * ainc * fbfrc * frc2;
        }

        for (int ntc = 1; ntc <= nstep; ++ntc) {
            for (int nk = 0; nk <= ltop; ++nk) {
                qlfxin[nk] = static_cast<T>(0.0);
                qlfxout[nk] = static_cast<T>(0.0);
                qifxin[nk] = static_cast<T>(0.0);
                qifxout[nk] = static_cast<T>(0.0);
                qrfxin[nk] = static_cast<T>(0.0);
                qrfxout[nk] = static_cast<T>(0.0);
                qsfxin[nk] = static_cast<T>(0.0);
                qsfxout[nk] = static_cast<T>(0.0);
            }
            for (int nk = 1; nk <= ltop; ++nk) {
                if (omg[nk] <= static_cast<T>(0.0)) {
                    qlfxin[nk] = -fxm[nk] * qlpa[nk - 1];
                    qifxin[nk] = -fxm[nk] * qipa[nk - 1];
                    qrfxin[nk] = -fxm[nk] * qrpa[nk - 1];
                    qsfxin[nk] = -fxm[nk] * qspa[nk - 1];
                    qlfxout[nk - 1] += qlfxin[nk];
                    qifxout[nk - 1] += qifxin[nk];
                    qrfxout[nk - 1] += qrfxin[nk];
                    qsfxout[nk - 1] += qsfxin[nk];
                } else {
                    qlfxout[nk] = fxm[nk] * qlpa[nk];
                    qifxout[nk] = fxm[nk] * qipa[nk];
                    qrfxout[nk] = fxm[nk] * qrpa[nk];
                    qsfxout[nk] = fxm[nk] * qspa[nk];
                    qlfxin[nk - 1] += qlfxout[nk];
                    qifxin[nk - 1] += qifxout[nk];
                    qrfxin[nk - 1] += qrfxout[nk];
                    qsfxin[nk - 1] += qsfxout[nk];
                }
            }
            for (int nk = 0; nk <= ltop; ++nk) {
                qlpa[nk] += (qlfxin[nk] + detlq[nk] - qlfxout[nk]) * dtime * emsd[nk];
                qipa[nk] += (qifxin[nk] + detic[nk] - qifxout[nk]) * dtime * emsd[nk];
                qrpa[nk] += (qrfxin[nk] - qrfxout[nk] + rainfb[nk]) * dtime * emsd[nk];
                qspa[nk] += (qsfxin[nk] - qsfxout[nk] + snowfb[nk]) * dtime * emsd[nk];
            }
        }
        for (int nk = 0; nk <= ltop; ++nk) {
            qlg[nk] = qlpa[nk];
            qig[nk] = qipa[nk];
            qrg[nk] = qrpa[nk];
            qsg[nk] = qspa[nk];
        }

        // Convective timescale & NCA
        if (tadvec < timec) {
            nic = static_cast<int>(std::round(tadvec / dt));
        }
        nca(i, j) = static_cast<T>(nic) * dt;
        if (ishall == 1) {
            timec = static_cast<T>(2400.0);
            nca(i, j) = cudt * static_cast<T>(60.0);
        }

        T dqcdt[kMaxLevels], dqidt[kMaxLevels], dqrdt[kMaxLevels], dqsdt[kMaxLevels];
        T dtdt[kMaxLevels], dqdt[kMaxLevels];

        for (int k_idx = 0; k_idx < kx; ++k_idx) {
            if (warm_rain) {
                T cpm = cp * (static_cast<T>(1.0) + static_cast<T>(0.887) * qg[k_idx]);
                tg[k_idx] -= (qig[k_idx] + qsg[k_idx]) * rlf / cpm;
                dqcdt[k_idx] = (qlg[k_idx] + qig[k_idx] - ql0[k_idx] - qi0[k_idx]) / timec;
                dqidt[k_idx] = static_cast<T>(0.0);
                dqrdt[k_idx] = (qrg[k_idx] + qsg[k_idx] - qr0[k_idx] - qs0[k_idx]) / timec;
                dqsdt[k_idx] = static_cast<T>(0.0);
            } else if (!flag_qs) {
                T cpm = cp * (static_cast<T>(1.0) + static_cast<T>(0.887) * qg[k_idx]);
                if (k_idx <= ml) {
                    tg[k_idx] -= (qig[k_idx] + qsg[k_idx]) * rlf / cpm;
                } else {
                    tg[k_idx] += (qlg[k_idx] + qrg[k_idx]) * rlf / cpm;
                }
                dqcdt[k_idx] = (qlg[k_idx] + qig[k_idx] - ql0[k_idx] - qi0[k_idx]) / timec;
                dqidt[k_idx] = static_cast<T>(0.0);
                dqrdt[k_idx] = (qrg[k_idx] + qsg[k_idx] - qr0[k_idx] - qs0[k_idx]) / timec;
                dqsdt[k_idx] = static_cast<T>(0.0);
            } else {
                dqcdt[k_idx] = (qlg[k_idx] - ql0[k_idx]) / timec;
                dqsdt[k_idx] = (qsg[k_idx] - qs0[k_idx]) / timec;
                dqrdt[k_idx] = (qrg[k_idx] - qr0[k_idx]) / timec;
                if (flag_qi) {
                    dqidt[k_idx] = (qig[k_idx] - qi0[k_idx]) / timec;
                } else {
                    dqsdt[k_idx] += (qig[k_idx] - qi0[k_idx]) / timec;
                }
            }
            dtdt[k_idx] = (tg[k_idx] - t0[k_idx]) / timec;
            dqdt[k_idx] = (qg[k_idx] - q0[k_idx]) / timec;
        }

        pratec(i, j) = pptflx * (static_cast<T>(1.0) - fbfrc) / dxsq;
        if (pratec(i, j) < static_cast<T>(0.0)) pratec(i, j) = static_cast<T>(0.0);
        raincv(i, j) = dt * pratec(i, j);
        if (raincv(i, j) < static_cast<T>(0.0)) raincv(i, j) = static_cast<T>(0.0);
        cutop(i, j) = static_cast<T>(ltop + 1); // 1-based level
        cubot(i, j) = static_cast<T>(lcl + 1);  // 1-based level

        // Feedback to 3D model tendencies
        for (int k_idx = 0; k_idx < kx; ++k_idx) {
            rthcuten(i, j, k_idx) = dtdt[k_idx] / pi(i, j, k_idx);
            rqvcuten(i, j, k_idx) = dqdt[k_idx];
            if (flag_qr) {
                rqrcuten(i, j, k_idx) = dqrdt[k_idx];
                rqccuten(i, j, k_idx) = dqcdt[k_idx];
            } else {
                rqrcuten(i, j, k_idx) = static_cast<T>(0.0);
                rqccuten(i, j, k_idx) = dqrdt[k_idx] + dqcdt[k_idx];
            }
            if (flag_qi) rqicuten(i, j, k_idx) = dqidt[k_idx];
            if (flag_qs) rqscuten(i, j, k_idx) = dqsdt[k_idx];
        }
    }
};

template <typename T>
void run_kf(KFKernel<T> kernel, const Range2D& range);

} // namespace cppwrf::phys::cu_kf
