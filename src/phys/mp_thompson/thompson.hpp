#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "phys/mp_thompson/thompson_tables.hpp"

namespace cppwrf::phys::mp_thompson {

inline constexpr int kMaxLevels = 128;

// Physical constants matching module_mp_thompson.F
inline constexpr double PI_CONST = 3.1415926536;
inline constexpr double rho_w = 1000.0;
inline constexpr double rho_i = 890.0;
inline constexpr double rho_s_default = 100.0;
inline constexpr double rho_g_val = 400.0; // idx_bg1 = 5, rho_g(5) = 400.0
inline constexpr double am_r = PI_CONST * rho_w / 6.0;
inline constexpr double bm_r = 3.0;
inline constexpr double am_s = 0.069;
inline constexpr double bm_s = 2.0;
inline constexpr double am_g5 = PI_CONST * rho_g_val / 6.0;
inline constexpr double bm_g = 3.0;
inline constexpr double am_i = PI_CONST * rho_i / 6.0;
inline constexpr double bm_i = 3.0;

inline constexpr double av_r = 4854.0;
inline constexpr double bv_r = 1.0;
inline constexpr double fv_r = 195.0;
inline constexpr double av_s = 40.0;
inline constexpr double bv_s = 0.55;
inline constexpr double fv_s = 100.0;
inline constexpr double av_g_old = 442.0;
inline constexpr double bv_g_old = 0.89;
inline constexpr double av_g5 = 442.0;
inline constexpr double bv_g5 = 0.89;
inline constexpr double av_i = 1493.9;
inline constexpr double bv_i = 1.0;
inline constexpr double av_c = 0.316946e8;
inline constexpr double bv_c = 2.0;

inline constexpr double C_cube = 0.5;
inline constexpr double C_sqrd = 0.15;
inline constexpr double Ef_si = 0.05;
inline constexpr double Ef_rs = 0.95;
inline constexpr double Ef_rg = 0.75;
inline constexpr double Ef_ri = 0.95;

inline constexpr double R1 = 1.0e-12;
inline constexpr double R2 = 1.0e-6;
inline constexpr double eps_val = 1.0e-15;

inline constexpr double TNO = 5.0;
inline constexpr double ATO = 0.304;
inline constexpr double RHO_NOT = 101325.0 / (287.05 * 298.0);
inline constexpr double HGFR = 235.16;

inline constexpr double Rv_val = 461.5;
inline constexpr double oRv_val = 1.0 / Rv_val;
inline constexpr double R_d = 287.04;
inline constexpr double Cp_val = 1004.0;

inline constexpr double lsub = 2.834e6;
inline constexpr double lvap0 = 2.5e6;
inline constexpr double lfus = lsub - lvap0;
inline constexpr double olfus = 1.0 / lfus;

inline constexpr double xm0i = 1.0e-12;
inline constexpr double D0c = 1.0e-6;
inline constexpr double D0r = 50.0e-6;
inline constexpr double D0s = 300.0e-6;
inline constexpr double D0g = 350.0e-6;
inline constexpr double D0i = 1.289198539077227e-05;
inline constexpr double Sc = 0.632;

inline constexpr double mu_r = 0.0;
inline constexpr double mu_g = 0.0;
inline constexpr double mu_i = 0.0;
inline constexpr double mu_s = 0.6357;
inline constexpr double Kap0 = 490.6;
inline constexpr double Kap1 = 17.46;
inline constexpr double Lam0 = 20.78;
inline constexpr double Lam1 = 3.29;

inline constexpr double Nt_c = 100.0e6;
inline constexpr double Nt_c_max = 1999.0e6;
inline constexpr double naIN1 = 0.5e6;

// Thermodynamic helper functions
template <typename Real>
CPPWRF_HD inline Real rslf(Real p, Real t) {
    constexpr Real C0 = Real(0.611583699e03);
    constexpr Real C1 = Real(0.444606896e02);
    constexpr Real C2 = Real(0.143177157e01);
    constexpr Real C3 = Real(0.264224321e-1);
    constexpr Real C4 = Real(0.299291081e-3);
    constexpr Real C5 = Real(0.203154182e-5);
    constexpr Real C6 = Real(0.702620698e-8);
    constexpr Real C7 = Real(0.379534310e-11);
    constexpr Real C8 = Real(-0.321582393e-13);

    Real x = (t - Real(273.16) > Real(-80.0)) ? t - Real(273.16) : Real(-80.0);
    Real esl = C0 + x * (C1 + x * (C2 + x * (C3 + x * (C4 + x * (C5 + x * (C6 + x * (C7 + x * C8)))))));
    Real esl_max = p * Real(0.15);
    if (esl > esl_max) esl = esl_max;
    return Real(0.622) * esl / (p - esl);
}

template <typename Real>
CPPWRF_HD inline Real rsif(Real p, Real t) {
    constexpr Real C0 = Real(0.609868993e03);
    constexpr Real C1 = Real(0.499320233e02);
    constexpr Real C2 = Real(0.184672631e01);
    constexpr Real C3 = Real(0.402737184e-1);
    constexpr Real C4 = Real(0.565392987e-3);
    constexpr Real C5 = Real(0.521693933e-5);
    constexpr Real C6 = Real(0.307839583e-7);
    constexpr Real C7 = Real(0.105785160e-9);
    constexpr Real C8 = Real(0.161444444e-12);

    Real x = (t - Real(273.16) > Real(-80.0)) ? t - Real(273.16) : Real(-80.0);
    Real esi = C0 + x * (C1 + x * (C2 + x * (C3 + x * (C4 + x * (C5 + x * (C6 + x * (C7 + x * C8)))))));
    Real esi_max = p * Real(0.15);
    if (esi > esi_max) esi = esi_max;
    Real denom = p - esi;
    if (denom < Real(1.0e-4)) denom = Real(1.0e-4);
    return Real(0.622) * esi / denom;
}

template <typename Real>
CPPWRF_HD inline Real t_dew(Real pres_Pa, Real w_non) {
    Real rr = w_non + Real(1.0e-8);
    Real es = pres_Pa * rr / (Real(0.622) + rr);
    Real esln = std::log(es);
    return (Real(35.86) * esln - Real(4947.2325)) / (esln - Real(23.6837));
}

template <typename Real>
CPPWRF_HD inline Real t_lcl(Real temp_K, Real tdew_K) {
    Real denom = (Real(1.0) / (tdew_K - Real(56.0))) + (std::log(temp_K / tdew_K) / Real(800.0));
    return (Real(1.0) / denom) + Real(56.0);
}

template <typename Real>
CPPWRF_HD inline Real theta_e(Real pres_Pa, Real temp_K, Real w_non, Real tlcl_K) {
    Real rr = w_non + Real(1.0e-8);
    Real power = Real(0.2854) * (Real(1.0) - Real(0.28) * rr);
    Real xx = temp_K * std::pow(Real(100000.0) / pres_Pa, power);
    Real p1 = (Real(3.376) / tlcl_K) - Real(0.00254);
    Real p2 = (rr * Real(1000.0)) * (Real(1.0) + Real(0.81) * rr);
    return xx * std::exp(p1 * p2);
}

template <typename Real>
CPPWRF_HD inline Real theta_wetb(Real thetae_K) {
    constexpr double c[7] = {
        -1.00922292e-10, -1.47945344e-8, -1.7303757e-6,
        -0.00012709, 1.15849867e-6, -3.518296861e-9, 3.5741522e-12
    };
    constexpr double d[7] = {
        0.00000000, -3.5223513e-10, -5.7250807e-8,
        -5.83975422e-6, 4.72445163e-8, -1.13402845e-10, 8.729580402e-14
    };
    double x = static_cast<double>(thetae_K);
    if (x > 475.0) x = 475.0;
    double answer = 0.0;
    if (x <= 335.5) {
        answer = c[0] + x * (c[1] + x * (c[2] + x * (c[3] + x * (c[4] + x * (c[5] + x * c[6])))));
    } else {
        answer = d[0] + x * (d[1] + x * (d[2] + x * (d[3] + x * (d[4] + x * (d[5] + x * d[6])))));
    }
    return static_cast<Real>(answer + 273.15);
}

template <typename Real>
CPPWRF_HD inline Real compT_fr_The(Real thelcl_K, Real pres_Pa) {
    Real diff_th = thelcl_K - Real(270.0);
    Real max_diff = (diff_th > Real(0.0)) ? diff_th : Real(0.0);
    Real guess = (thelcl_K - Real(0.5) * std::pow(max_diff, Real(1.05)))
               * std::pow(pres_Pa / Real(100000.0), Real(0.2));
    constexpr Real epsilon = Real(0.01);
    for (int iter = 0; iter < 100; ++iter) {
        Real w1 = rslf(pres_Pa, guess);
        Real w2 = rslf(pres_Pa, guess + Real(1.0));
        Real tenu = theta_e(pres_Pa, guess, w1, guess);
        Real tenup = theta_e(pres_Pa, guess + Real(1.0), w2, guess + Real(1.0));
        Real cor = (thelcl_K - tenu) / (tenup - tenu);
        guess += cor;
        if (std::abs(cor) < epsilon) {
            return guess;
        }
    }
    Real thwlcl_K = theta_wetb(thelcl_K);
    return thwlcl_K * std::pow(pres_Pa / Real(100000.0), Real(0.286));
}

// Single-source Thompson column kernel
template <typename Real>
struct ThompsonKernel {
    ThompsonTablesView<Real> tables;
    Real dt;
    int nz;

    View3<Real> th;
    View3<Real> qv;
    View3<Real> qc;
    View3<Real> qr;
    View3<Real> qi;
    View3<Real> qs;
    View3<Real> qg;
    View3<Real> ni;
    View3<Real> nr;

    View3<const Real> pii;
    View3<const Real> p;
    View3<const Real> w;
    View3<const Real> dz;

    View2<Real> rainnc;
    View2<Real> rainncv;
    View2<Real> snownc;
    View2<Real> snowncv;
    View2<Real> graupelnc;
    View2<Real> graupelncv;
    View2<Real> sr;

    CPPWRF_HD void operator()(int i, int j) const {
        if (nz > kMaxLevels) return;

        Real temp[kMaxLevels], twet[kMaxLevels], pres[kMaxLevels], qv_loc[kMaxLevels];
        Real qc1d[kMaxLevels], qi1d[kMaxLevels], qr1d[kMaxLevels], qs1d[kMaxLevels], qg1d[kMaxLevels], qb1d[kMaxLevels];
        Real nc1d[kMaxLevels], ni1d[kMaxLevels], nr1d[kMaxLevels], ng1d[kMaxLevels];
        Real rc[kMaxLevels], ri[kMaxLevels], rr[kMaxLevels], rs[kMaxLevels], rg[kMaxLevels], rb[kMaxLevels];
        Real ni_loc[kMaxLevels], nr_loc[kMaxLevels], nc[kMaxLevels], ns[kMaxLevels], ng[kMaxLevels];
        Real rho[kMaxLevels], rhof[kMaxLevels], rhof2[kMaxLevels];
        Real qvs[kMaxLevels], qvsi[kMaxLevels], delQvs[kMaxLevels];
        Real satw[kMaxLevels], sati[kMaxLevels], ssatw[kMaxLevels], ssati[kMaxLevels];
        Real diffu[kMaxLevels], visco[kMaxLevels], vsc2[kMaxLevels], tcond[kMaxLevels], lvap[kMaxLevels], ocp[kMaxLevels], lvt2[kMaxLevels];
        Real ilamr[kMaxLevels], ilamg[kMaxLevels], N0_r[kMaxLevels], N0_g[kMaxLevels];
        Real mvd_r[kMaxLevels], mvd_c[kMaxLevels], mvd_g[kMaxLevels];
        Real smo0[kMaxLevels], smo1[kMaxLevels], smo2[kMaxLevels], smob[kMaxLevels], smoc[kMaxLevels], smod[kMaxLevels], smoe[kMaxLevels], smof[kMaxLevels], smog[kMaxLevels];
        Real vts_boost[kMaxLevels];
        Real vtik[kMaxLevels + 1], vtnik[kMaxLevels + 1], vtrk[kMaxLevels + 1], vtnrk[kMaxLevels + 1];
        Real vtsk[kMaxLevels + 1], vtgk[kMaxLevels + 1], vtngk[kMaxLevels + 1], vtck[kMaxLevels + 1], vtnck[kMaxLevels + 1];
        Real sed_r[kMaxLevels], sed_s[kMaxLevels], sed_g[kMaxLevels], sed_i[kMaxLevels], sed_n[kMaxLevels], sed_c[kMaxLevels], sed_b[kMaxLevels];

        Real tten[kMaxLevels], qvten[kMaxLevels], qcten[kMaxLevels], qiten[kMaxLevels], qrten[kMaxLevels];
        Real qsten[kMaxLevels], qgten[kMaxLevels], qbten[kMaxLevels], niten[kMaxLevels], nrten[kMaxLevels], ncten[kMaxLevels], ngten[kMaxLevels];

        Real prw_vcd[kMaxLevels], pnc_wcd[kMaxLevels], pnc_wau[kMaxLevels], pnc_rcw[kMaxLevels], pnc_scw[kMaxLevels], pnc_gcw[kMaxLevels];
        Real prr_wau[kMaxLevels], prr_rcw[kMaxLevels], prr_rcs[kMaxLevels], prr_rcg[kMaxLevels], prr_sml[kMaxLevels], prr_gml[kMaxLevels], prr_rci[kMaxLevels], prv_rev[kMaxLevels];
        Real pnr_wau[kMaxLevels], pnr_rcs[kMaxLevels], pnr_rcg[kMaxLevels], pnr_rci[kMaxLevels], pnr_sml[kMaxLevels], pnr_gml[kMaxLevels], pnr_rev[kMaxLevels], pnr_rcr[kMaxLevels], pnr_rfz[kMaxLevels];
        Real pri_inu[kMaxLevels], pni_inu[kMaxLevels], pri_ihm[kMaxLevels], pni_ihm[kMaxLevels], pri_wfz[kMaxLevels], pni_wfz[kMaxLevels], pri_rfz[kMaxLevels], pni_rfz[kMaxLevels];
        Real pri_ide[kMaxLevels], pni_ide[kMaxLevels], pri_rci[kMaxLevels], pni_rci[kMaxLevels], pni_sci[kMaxLevels], pni_iau[kMaxLevels], pri_iha[kMaxLevels], pni_iha[kMaxLevels];
        Real prs_iau[kMaxLevels], prs_sci[kMaxLevels], prs_rcs[kMaxLevels], prs_scw[kMaxLevels], prs_sde[kMaxLevels], prs_ihm[kMaxLevels], prs_ide[kMaxLevels];
        Real prg_scw[kMaxLevels], prg_rfz[kMaxLevels], prg_gde[kMaxLevels], prg_gcw[kMaxLevels], prg_rci[kMaxLevels], prg_rcs[kMaxLevels], prg_rcg[kMaxLevels], prg_ihm[kMaxLevels];
        Real png_scw[kMaxLevels], png_rcs[kMaxLevels], png_rcg[kMaxLevels], png_gde[kMaxLevels];
        Real pbg_scw[kMaxLevels], pbg_rfz[kMaxLevels], pbg_gcw[kMaxLevels], pbg_rci[kMaxLevels], pbg_rcs[kMaxLevels], pbg_rcg[kMaxLevels], pbg_sml[kMaxLevels], pbg_gml[kMaxLevels];

        bool L_qc[kMaxLevels], L_qi[kMaxLevels], L_qr[kMaxLevels], L_qs[kMaxLevels], L_qg[kMaxLevels];

        Real pptrain = Real(0.0);
        Real pptsnow = Real(0.0);
        Real pptgraul = Real(0.0);
        Real pptice = Real(0.0);

        Real dtsave = dt;
        Real odt = Real(1.0) / dt;
        Real odts = Real(1.0) / dtsave;
        Real obmr = Real(1.0) / Real(bm_r);
        Real obmi = Real(1.0) / Real(bm_i);

        bool no_micro = true;

        for (int k = 0; k < nz; ++k) {
            tten[k] = Real(0.0);
            qvten[k] = Real(0.0);
            qcten[k] = Real(0.0);
            qiten[k] = Real(0.0);
            qrten[k] = Real(0.0);
            qsten[k] = Real(0.0);
            qgten[k] = Real(0.0);
            ngten[k] = Real(0.0);
            qbten[k] = Real(0.0);
            niten[k] = Real(0.0);
            nrten[k] = Real(0.0);
            ncten[k] = Real(0.0);

            prw_vcd[k] = Real(0.0);
            pnc_wcd[k] = Real(0.0);
            pnc_wau[k] = Real(0.0);
            pnc_rcw[k] = Real(0.0);
            pnc_scw[k] = Real(0.0);
            pnc_gcw[k] = Real(0.0);

            prv_rev[k] = Real(0.0);
            prr_wau[k] = Real(0.0);
            prr_rcw[k] = Real(0.0);
            prr_rcs[k] = Real(0.0);
            prr_rcg[k] = Real(0.0);
            prr_sml[k] = Real(0.0);
            prr_gml[k] = Real(0.0);
            prr_rci[k] = Real(0.0);
            pnr_wau[k] = Real(0.0);
            pnr_rcs[k] = Real(0.0);
            pnr_rcg[k] = Real(0.0);
            pnr_rci[k] = Real(0.0);
            pnr_sml[k] = Real(0.0);
            pnr_gml[k] = Real(0.0);
            pnr_rev[k] = Real(0.0);
            pnr_rcr[k] = Real(0.0);
            pnr_rfz[k] = Real(0.0);

            pri_inu[k] = Real(0.0);
            pni_inu[k] = Real(0.0);
            pri_ihm[k] = Real(0.0);
            pni_ihm[k] = Real(0.0);
            pri_wfz[k] = Real(0.0);
            pni_wfz[k] = Real(0.0);
            pri_rfz[k] = Real(0.0);
            pni_rfz[k] = Real(0.0);
            pri_ide[k] = Real(0.0);
            pni_ide[k] = Real(0.0);
            pri_rci[k] = Real(0.0);
            pni_rci[k] = Real(0.0);
            pni_sci[k] = Real(0.0);
            pni_iau[k] = Real(0.0);
            pri_iha[k] = Real(0.0);
            pni_iha[k] = Real(0.0);

            prs_iau[k] = Real(0.0);
            prs_sci[k] = Real(0.0);
            prs_rcs[k] = Real(0.0);
            prs_scw[k] = Real(0.0);
            prs_sde[k] = Real(0.0);
            prs_ihm[k] = Real(0.0);
            prs_ide[k] = Real(0.0);

            prg_scw[k] = Real(0.0);
            prg_rfz[k] = Real(0.0);
            prg_gde[k] = Real(0.0);
            prg_gcw[k] = Real(0.0);
            prg_rci[k] = Real(0.0);
            prg_rcs[k] = Real(0.0);
            prg_rcg[k] = Real(0.0);
            prg_ihm[k] = Real(0.0);
            png_scw[k] = Real(0.0);
            png_rcs[k] = Real(0.0);
            png_rcg[k] = Real(0.0);
            png_gde[k] = Real(0.0);

            pbg_scw[k] = Real(0.0);
            pbg_rfz[k] = Real(0.0);
            pbg_gcw[k] = Real(0.0);
            pbg_rci[k] = Real(0.0);
            pbg_rcs[k] = Real(0.0);
            pbg_rcg[k] = Real(0.0);
            pbg_sml[k] = Real(0.0);
            pbg_gml[k] = Real(0.0);

            smo0[k] = Real(0.0);
            smo1[k] = Real(0.0);
            smo2[k] = Real(0.0);
            smob[k] = Real(0.0);
            smoc[k] = Real(0.0);
            smod[k] = Real(0.0);
            smoe[k] = Real(0.0);
            smof[k] = Real(0.0);
            smog[k] = Real(0.0);
            ns[k]   = Real(0.0);
            mvd_r[k] = Real(0.0);
            mvd_c[k] = Real(0.0);
            mvd_g[k] = Real(0.0);

            temp[k] = th(i, j, k) * pii(i, j, k);
            pres[k] = p(i, j, k);
            Real qv_in = qv(i, j, k);
            qc1d[k] = qc(i, j, k);
            qi1d[k] = qi(i, j, k);
            qr1d[k] = qr(i, j, k);
            qs1d[k] = qs(i, j, k);
            qg1d[k] = qg(i, j, k);
            ni1d[k] = ni(i, j, k);
            nr1d[k] = nr(i, j, k);
            Real rho_drv = Real(0.622) * pres[k] / (Real(R_d) * temp[k] * (qv_in + Real(0.622)));
            nc1d[k] = Real(Nt_c) / rho_drv;
            if (qg1d[k] > Real(R1)) {
                Real rg_drv = qg1d[k] * rho_drv;
                Real ygra1 = std::log10((rg_drv > Real(1.0e-9)) ? rg_drv : Real(1.0e-9));
                Real zans1 = Real(3.0) + Real(2.0 / 7.0) * (ygra1 + Real(8.0));
                if (zans1 < Real(2.0)) zans1 = Real(2.0);
                if (zans1 > Real(6.0)) zans1 = Real(6.0);
                Real N0_exp = std::pow(Real(10.0), zans1);
                Real lam_exp = std::pow(N0_exp * Real(am_g5) * tables.cgg[0] / (rho_drv * qg1d[k]), tables.oge1);
                Real lamg = lam_exp * std::pow(tables.cgg[2] * tables.ogg2 * tables.ogg1, tables.obmg);
                ng1d[k] = tables.cgg[1] * tables.ogg3 * rho_drv * qg1d[k] * std::pow(lamg, Real(bm_g)) / Real(am_g5);
                ng1d[k] = (ng1d[k] / rho_drv > Real(R2)) ? (ng1d[k] / rho_drv) : Real(R2);
                qb1d[k] = qg1d[k] / Real(rho_g_val);
            } else {
                ng1d[k] = Real(0.0);
                qb1d[k] = Real(0.0);
            }

            qv_loc[k] = (qv_in > Real(1.0e-10)) ? qv_in : Real(1.0e-10);
            rho[k] = Real(0.622) * pres[k] / (Real(R_d) * temp[k] * (qv_loc[k] + Real(0.622)));

            // Cloud water
            if (qc1d[k] > Real(R1)) {
                no_micro = false;
                rc[k] = qc1d[k] * rho[k];
                nc[k] = Real(Nt_c);
                L_qc[k] = true;
            } else {
                qc1d[k] = Real(0.0);
                nc1d[k] = Real(0.0);
                rc[k] = Real(R1);
                nc[k] = Real(2.0);
                L_qc[k] = false;
            }

            // Cloud ice
            if (qi1d[k] > Real(R1)) {
                no_micro = false;
                ri[k] = qi1d[k] * rho[k];
                Real ni_in = ni1d[k] * rho[k];
                ni_loc[k] = (ni_in > Real(R2)) ? ni_in : Real(R2);
                if (ni_loc[k] <= Real(R2)) {
                    Real lami = tables.cie[1] / Real(5.0e-6);
                    Real term = tables.cig[0] * tables.oig2 * ri[k] / Real(am_i) * std::pow(lami, Real(bm_i));
                    ni_loc[k] = (term < Real(999.0e3)) ? term : Real(999.0e3);
                }
                L_qi[k] = true;
                Real lami = std::pow(Real(am_i) * tables.cig[1] * tables.oig1 * ni_loc[k] / ri[k], obmi);
                Real ilami = Real(1.0) / lami;
                Real xDi = (Real(bm_i) + Real(mu_i) + Real(1.0)) * ilami;
                if (xDi < Real(5.0e-6)) {
                    lami = tables.cie[1] / Real(5.0e-6);
                    Real term = tables.cig[0] * tables.oig2 * ri[k] / Real(am_i) * std::pow(lami, Real(bm_i));
                    ni_loc[k] = (term < Real(999.0e3)) ? term : Real(999.0e3);
                } else if (xDi > Real(300.0e-6)) {
                    lami = tables.cie[1] / Real(300.0e-6);
                    ni_loc[k] = tables.cig[0] * tables.oig2 * ri[k] / Real(am_i) * std::pow(lami, Real(bm_i));
                }
            } else {
                qi1d[k] = Real(0.0);
                ni1d[k] = Real(0.0);
                ri[k] = Real(R1);
                ni_loc[k] = Real(R2);
                L_qi[k] = false;
            }

            // Rain
            if (qr1d[k] > Real(R1)) {
                no_micro = false;
                rr[k] = qr1d[k] * rho[k];
                Real nr_in = nr1d[k] * rho[k];
                nr_loc[k] = (nr_in > Real(R2)) ? nr_in : Real(R2);
                if (nr_loc[k] <= Real(R2)) {
                    mvd_r[k] = Real(1.0e-3);
                    Real lamr = (Real(3.0) + Real(mu_r) + Real(0.672)) / mvd_r[k];
                    nr_loc[k] = tables.crg[1] * tables.org3 * rr[k] * std::pow(lamr, Real(bm_r)) / Real(am_r);
                }
                L_qr[k] = true;
                Real lamr = std::pow(Real(am_r) * tables.crg[2] * tables.org2 * nr_loc[k] / rr[k], obmr);
                mvd_r[k] = (Real(3.0) + Real(mu_r) + Real(0.672)) / lamr;
                if (mvd_r[k] > Real(2.5e-3)) {
                    mvd_r[k] = Real(2.5e-3);
                    lamr = (Real(3.0) + Real(mu_r) + Real(0.672)) / mvd_r[k];
                    nr_loc[k] = tables.crg[1] * tables.org3 * rr[k] * std::pow(lamr, Real(bm_r)) / Real(am_r);
                } else if (mvd_r[k] < Real(D0r) * Real(0.75)) {
                    mvd_r[k] = Real(D0r) * Real(0.75);
                    lamr = (Real(3.0) + Real(mu_r) + Real(0.672)) / mvd_r[k];
                    nr_loc[k] = tables.crg[1] * tables.org3 * rr[k] * std::pow(lamr, Real(bm_r)) / Real(am_r);
                }
            } else {
                qr1d[k] = Real(0.0);
                nr1d[k] = Real(0.0);
                rr[k] = Real(R1);
                nr_loc[k] = Real(R2);
                L_qr[k] = false;
            }

            // Snow
            if (qs1d[k] > Real(R1)) {
                no_micro = false;
                rs[k] = qs1d[k] * rho[k];
                L_qs[k] = true;
            } else {
                qs1d[k] = Real(0.0);
                rs[k] = Real(R1);
                L_qs[k] = false;
            }

            // Graupel
            if (qg1d[k] > Real(R1)) {
                no_micro = false;
                L_qg[k] = true;
                rg[k] = qg1d[k] * rho[k];
                Real ng_in = ng1d[k] * rho[k];
                ng[k] = (ng_in > Real(R2)) ? ng_in : Real(R2);
                rb[k] = qb1d[k];
                if (ng[k] <= Real(R2)) {
                    mvd_g[k] = Real(1.5e-3);
                    Real lamg = (Real(3.0) + Real(mu_g) + Real(0.672)) / mvd_g[k];
                    ng[k] = tables.cgg[1] * tables.ogg3 * rg[k] * std::pow(lamg, Real(bm_g)) / Real(am_g5);
                }
                Real lamg = std::pow(Real(am_g5) * tables.cgg[2] * tables.ogg2 * ng[k] / rg[k], tables.obmg);
                mvd_g[k] = (Real(3.0) + Real(mu_g) + Real(0.672)) / lamg;
                if (mvd_g[k] > Real(25.4e-3)) {
                    mvd_g[k] = Real(25.4e-3);
                    lamg = (Real(3.0) + Real(mu_g) + Real(0.672)) / mvd_g[k];
                    ng[k] = tables.cgg[1] * tables.ogg3 * rg[k] * std::pow(lamg, Real(bm_g)) / Real(am_g5);
                } else if (mvd_g[k] < Real(D0r)) {
                    mvd_g[k] = Real(D0r);
                    lamg = (Real(3.0) + Real(mu_g) + Real(0.672)) / mvd_g[k];
                    ng[k] = tables.cgg[1] * tables.ogg3 * rg[k] * std::pow(lamg, Real(bm_g)) / Real(am_g5);
                }
            } else {
                qg1d[k] = Real(0.0);
                ng1d[k] = Real(0.0);
                qb1d[k] = Real(0.0);
                rg[k] = Real(R1);
                ng[k] = Real(R2);
                rb[k] = Real(R1) / rho[k] / Real(800.0);
                L_qg[k] = false;
            }
        }

        int k_melting = -1;
        for (int k = 0; k < nz; ++k) {
            Real tempc = temp[k] - Real(273.15);
            rhof[k] = std::sqrt(Real(RHO_NOT) / rho[k]);
            rhof2[k] = std::sqrt(rhof[k]);
            qvs[k] = rslf(pres[k], temp[k]);
            Real delq = rslf(pres[k], Real(273.15)) - qv_loc[k];
            delQvs[k] = (delq > Real(0.0)) ? delq : Real(0.0);
            if (tempc <= Real(0.0)) {
                qvsi[k] = rsif(pres[k], temp[k]);
            } else {
                qvsi[k] = qvs[k];
                if (k > k_melting) k_melting = k;
            }
            satw[k] = qv_loc[k] / qvs[k];
            sati[k] = qv_loc[k] / qvsi[k];
            ssatw[k] = satw[k] - Real(1.0);
            ssati[k] = sati[k] - Real(1.0);
            if (std::abs(ssatw[k]) < Real(eps_val)) ssatw[k] = Real(0.0);
            if (std::abs(ssati[k]) < Real(eps_val)) ssati[k] = Real(0.0);
            if (no_micro && ssati[k] > Real(0.0)) no_micro = false;

            diffu[k] = Real(2.11e-5) * std::pow(temp[k] / Real(273.15), Real(1.94)) * (Real(101325.0) / pres[k]);
            if (tempc >= Real(0.0)) {
                visco[k] = (Real(1.718) + Real(0.0049) * tempc) * Real(1.0e-5);
            } else {
                visco[k] = (Real(1.718) + Real(0.0049) * tempc - Real(1.2e-5) * tempc * tempc) * Real(1.0e-5);
            }
            ocp[k] = Real(1.0) / (Real(Cp_val) * (Real(1.0) + Real(0.887) * qv_loc[k]));
            vsc2[k] = std::sqrt(rho[k] / visco[k]);
            lvap[k] = Real(lvap0) + (Real(2106.0) - Real(4218.0)) * tempc;
            tcond[k] = (Real(5.69) + Real(0.0168) * tempc) * Real(1.0e-5) * Real(418.936);  // WRF: (..)*1.0E-5 * 418.936
            twet[k] = temp[k];
        }

        if (k_melting >= 0) {
            for (int k = 0; k <= k_melting; ++k) {
                if (satw[k] < Real(0.999)) {
                    Real td = t_dew(pres[k], qv_loc[k]);
                    Real dew_t = (temp[k] - Real(0.001) < td) ? temp[k] - Real(0.001) : td;
                    Real Tlcl = t_lcl(temp[k], dew_t);
                    Real The = theta_e(pres[k], temp[k], qv_loc[k], Tlcl);
                    Real t_comp = compT_fr_The(The, pres[k]);
                    twet[k] = (temp[k] < t_comp) ? temp[k] : t_comp;
                }
            }
        }

        if (no_micro) {
            for (int k = 0; k < nz; ++k) {
                qc(i, j, k) = qc1d[k];
                qi(i, j, k) = qi1d[k];
                qr(i, j, k) = qr1d[k];
                qs(i, j, k) = qs1d[k];
                qg(i, j, k) = qg1d[k];
                ni(i, j, k) = ni1d[k];
                nr(i, j, k) = nr1d[k];
                th(i, j, k) = temp[k] / pii(i, j, k);
            }
            rainncv(i, j) = Real(0.0);
            snowncv(i, j) = Real(0.0);
            graupelncv(i, j) = Real(0.0);
            sr(i, j) = Real(0.0);
            return;
        }

        // Snow moments
        for (int k = 0; k < nz; ++k) {
            if (!L_qs[k]) continue;
            Real tc0 = (temp[k] - Real(273.15) < Real(-0.1)) ? temp[k] - Real(273.15) : Real(-0.1);
            smob[k] = rs[k] * tables.oams;
            smo2[k] = smob[k]; // bm_s == 2.0

            Real loga0 = tables.sa[0] + tables.sa[1]*tc0 + tables.sa[4]*tc0*tc0 + tables.sa[8]*tc0*tc0*tc0;
            Real a0 = std::pow(Real(10.0), loga0);
            Real b0 = tables.sb[0] + tables.sb[1]*tc0 + tables.sb[4]*tc0*tc0 + tables.sb[8]*tc0*tc0*tc0;
            smo0[k] = a0 * std::pow(smo2[k], b0);

            Real loga1 = tables.sa[0] + tables.sa[1]*tc0 + tables.sa[2]
                       + tables.sa[3]*tc0 + tables.sa[4]*tc0*tc0
                       + tables.sa[5] + tables.sa[6]*tc0*tc0
                       + tables.sa[7]*tc0 + tables.sa[8]*tc0*tc0*tc0
                       + tables.sa[9];
            Real a1 = std::pow(Real(10.0), loga1);
            Real b1 = tables.sb[0] + tables.sb[1]*tc0 + tables.sb[2] + tables.sb[3]*tc0
                    + tables.sb[4]*tc0*tc0 + tables.sb[5]
                    + tables.sb[6]*tc0*tc0 + tables.sb[7]*tc0
                    + tables.sb[8]*tc0*tc0*tc0 + tables.sb[9];
            smo1[k] = a1 * std::pow(smo2[k], b1);

            Real cse1 = tables.cse[0];
            Real logac = tables.sa[0] + tables.sa[1]*tc0 + tables.sa[2]*cse1
                       + tables.sa[3]*tc0*cse1 + tables.sa[4]*tc0*tc0
                       + tables.sa[5]*cse1*cse1 + tables.sa[6]*tc0*tc0*cse1
                       + tables.sa[7]*tc0*cse1*cse1 + tables.sa[8]*tc0*tc0*tc0
                       + tables.sa[9]*cse1*cse1*cse1;
            Real ac = std::pow(Real(10.0), logac);
            Real bc = tables.sb[0] + tables.sb[1]*tc0 + tables.sb[2]*cse1 + tables.sb[3]*tc0*cse1
                    + tables.sb[4]*tc0*tc0 + tables.sb[5]*cse1*cse1
                    + tables.sb[6]*tc0*tc0*cse1 + tables.sb[7]*tc0*cse1*cse1
                    + tables.sb[8]*tc0*tc0*tc0 + tables.sb[9]*cse1*cse1*cse1;
            smoc[k] = ac * std::pow(smo2[k], bc);

            Real M0 = smob[k] / smoc[k];
            Real Mrat = smob[k] * M0 * M0 * M0;
            Real slam1 = M0 * Real(Lam0);
            Real slam2 = M0 * Real(Lam1);
            ns[k] = Mrat * Real(Kap0) / slam1 + Mrat * Real(Kap1) * std::pow(M0, Real(mu_s)) * tables.csg[14] / std::pow(slam2, tables.cse[14]);

            Real cse13 = tables.cse[12];
            Real logae = tables.sa[0] + tables.sa[1]*tc0 + tables.sa[2]*cse13
                       + tables.sa[3]*tc0*cse13 + tables.sa[4]*tc0*tc0
                       + tables.sa[5]*cse13*cse13 + tables.sa[6]*tc0*tc0*cse13
                       + tables.sa[7]*tc0*cse13*cse13 + tables.sa[8]*tc0*tc0*tc0
                       + tables.sa[9]*cse13*cse13*cse13;
            Real ae = std::pow(Real(10.0), logae);
            Real be = tables.sb[0] + tables.sb[1]*tc0 + tables.sb[2]*cse13 + tables.sb[3]*tc0*cse13
                    + tables.sb[4]*tc0*tc0 + tables.sb[5]*cse13*cse13
                    + tables.sb[6]*tc0*tc0*cse13 + tables.sb[7]*tc0*cse13*cse13
                    + tables.sb[8]*tc0*tc0*tc0 + tables.sb[9]*cse13*cse13*cse13;
            smoe[k] = ae * std::pow(smo2[k], be);

            Real cse16 = tables.cse[15];
            Real logaf = tables.sa[0] + tables.sa[1]*tc0 + tables.sa[2]*cse16
                       + tables.sa[3]*tc0*cse16 + tables.sa[4]*tc0*tc0
                       + tables.sa[5]*cse16*cse16 + tables.sa[6]*tc0*tc0*cse16
                       + tables.sa[7]*tc0*cse16*cse16 + tables.sa[8]*tc0*tc0*tc0
                       + tables.sa[9]*cse16*cse16*cse16;
            Real af = std::pow(Real(10.0), logaf);
            Real bf = tables.sb[0] + tables.sb[1]*tc0 + tables.sb[2]*cse16 + tables.sb[3]*tc0*cse16
                    + tables.sb[4]*tc0*tc0 + tables.sb[5]*cse16*cse16
                    + tables.sb[6]*tc0*tc0*cse16 + tables.sb[7]*tc0*cse16*cse16
                    + tables.sb[8]*tc0*tc0*tc0 + tables.sb[9]*cse16*cse16*cse16;
            smof[k] = af * std::pow(smo2[k], bf);

            Real cse17 = tables.cse[16];
            Real logag = tables.sa[0] + tables.sa[1]*tc0 + tables.sa[2]*cse17
                       + tables.sa[3]*tc0*cse17 + tables.sa[4]*tc0*tc0
                       + tables.sa[5]*cse17*cse17 + tables.sa[6]*tc0*tc0*cse17
                       + tables.sa[7]*tc0*cse17*cse17 + tables.sa[8]*tc0*tc0*tc0
                       + tables.sa[9]*cse17*cse17*cse17;
            Real ag = std::pow(Real(10.0), logag);
            Real bg = tables.sb[0] + tables.sb[1]*tc0 + tables.sb[2]*cse17 + tables.sb[3]*tc0*cse17
                    + tables.sb[4]*tc0*tc0 + tables.sb[5]*cse17*cse17
                    + tables.sb[6]*tc0*tc0*cse17 + tables.sb[7]*tc0*cse17*cse17
                    + tables.sb[8]*tc0*tc0*tc0 + tables.sb[9]*cse17*cse17*cse17;
            smog[k] = ag * std::pow(smo2[k], bg);
        }

        for (int k = nz - 1; k >= 0; --k) {
            Real lamg = std::pow(Real(am_g5) * tables.cgg[2] * tables.ogg2 * ng[k] / rg[k], tables.obmg);
            ilamg[k] = Real(1.0) / lamg;
            N0_g[k] = ng[k] * tables.ogg2 * std::pow(lamg, tables.cge[1]);
        }

        // gfortran vectorises this k loop (kte down to kts): libmvec's 2-lane pow (cppwrf::vm)
        for (int k = nz - 1; k >= 0; --k) {
            Real lamr = vm::pow(Real(am_r) * tables.crg[2] * tables.org2 * nr_loc[k] / rr[k], obmr, nz - 1 - k, nz);
            ilamr[k] = Real(1.0) / lamr;
            mvd_r[k] = (Real(3.0) + Real(mu_r) + Real(0.672)) / lamr;
            N0_r[k] = nr_loc[k] * tables.org2 * vm::pow(lamr, tables.cre[1], nz - 1 - k, nz);
        }

        // Warm-rain processes
        for (int k = 0; k < nz; ++k) {
            if (L_qr[k] && mvd_r[k] > Real(D0r)) {
                Real Ef_rr = Real(1.0) - std::exp(Real(2300.0) * (mvd_r[k] - Real(1950.0e-6)));
                pnr_rcr[k] = Ef_rr * Real(2.0) * nr_loc[k] * rr[k];
            }

            mvd_c[k] = Real(D0c);
            int nu_c = 15;
            Real lamc = Real(0.0);
            Real xDc = Real(0.0);
            if (L_qc[k]) {
                nu_c = static_cast<int>(std::round(Real(1000.0e6) / nc[k])) + 2;
                if (nu_c > 15) nu_c = 15;
                int nu_idx = nu_c - 1;
                Real xDc_calc = std::pow(rc[k] / (Real(am_r) * nc[k]), obmr) * Real(1.0e6);
                xDc = (Real(D0c * 1.0e6) > xDc_calc) ? Real(D0c * 1.0e6) : xDc_calc;
                lamc = std::pow(nc[k] * Real(am_r) * tables.ccg[1][nu_idx] * tables.ocg1[nu_idx] / rc[k], obmr);
                mvd_c[k] = (Real(3.0) + Real(nu_c) + Real(0.672)) / lamc;
                if (mvd_c[k] < Real(D0c)) mvd_c[k] = Real(D0c);
                if (mvd_c[k] > Real(D0r)) mvd_c[k] = Real(D0r);
            }

            if (rc[k] > Real(0.01e-3)) {
                int nu_idx = nu_c - 1;
                Real Dc_g = (std::pow(tables.ccg[2][nu_idx] * tables.ocg2[nu_idx], obmr) / lamc) * Real(1.0e6);
                // WRF's left-to-right products: xDc*xDc*xDc*Dc_g*Dc_g*Dc_g - xDc*xDc*xDc*xDc*xDc*xDc
                Real diff6 = xDc * xDc * xDc * Dc_g * Dc_g * Dc_g - xDc * xDc * xDc * xDc * xDc * xDc;
                Real Dc_b = (diff6 > Real(0.0)) ? std::pow(diff6, Real(1.0 / 6.0)) : Real(0.0);
                Real arg_zeta = Real(6.25e-6) * xDc * Dc_b * Dc_b * Dc_b - Real(0.4);
                Real zeta1 = Real(0.5) * (arg_zeta + std::abs(arg_zeta));
                Real zeta = Real(0.027) * rc[k] * zeta1;
                Real arg_tau = Real(0.5) * Dc_b - Real(7.5);
                Real taud = Real(0.5) * (arg_tau + std::abs(arg_tau)) + Real(R1);
                Real tau = Real(3.72) / (rc[k] * taud);
                prr_wau[k] = zeta / tau;
                Real max_wau = rc[k] * odts;
                if (prr_wau[k] > max_wau) prr_wau[k] = max_wau;
                pnr_wau[k] = prr_wau[k] / (Real(am_r) * Real(nu_c) * Real(10.0) * Real(D0r) * Real(D0r) * Real(D0r));
                Real mvd_c3 = mvd_c[k] * mvd_c[k] * mvd_c[k];
                Real pnc_calc = prr_wau[k] / (Real(am_r) * mvd_c3);
                Real max_pnc = nc[k] * odts;
                pnc_wau[k] = (pnc_calc < max_pnc) ? pnc_calc : max_pnc;
            }

            if (L_qr[k] && mvd_r[k] > Real(D0r) && mvd_c[k] > Real(D0c)) {
                Real lamr = Real(1.0) / ilamr[k];
                int idx = 1 + static_cast<int>(Real(nbr) * std::log(mvd_r[k] / tables.Dr_0) / tables.log_Dr);
                if (idx > nbr) idx = nbr;
                if (idx < 1) idx = 1;
                int j_c = static_cast<int>(mvd_c[k] * Real(1.0e6));
                if (j_c > nbc) j_c = nbc;
                if (j_c < 1) j_c = 1;

                Real Ef_rw = tables.t_Efrw[tables.idx_efrw(idx - 1, j_c - 1)];
                Real rate = rhof[k] * tables.t1_qr_qc * Ef_rw * rc[k] * N0_r[k] * std::pow(lamr + Real(fv_r), -tables.cre[8]);
                prr_rcw[k] = (rate < rc[k] * odts) ? rate : rc[k] * odts;
                Real rate_nc = rhof[k] * tables.t1_qr_qc * Ef_rw * nc[k] * N0_r[k] * std::pow(lamr + Real(fv_r), -tables.cre[8]);
                pnc_rcw[k] = (rate_nc < nc[k] * odts) ? rate_nc : nc[k] * odts;
            }
        }

        // Frozen processes
        for (int k = 0; k < nz; ++k) {
            vts_boost[k] = Real(1.0);
            Real orho = Real(1.0) / rho[k];

            Real xDs = Real(0.0);
            Real rho_s = Real(rho_s_default);
            if (L_qs[k]) {
                xDs = smoc[k] / smob[k];
                Real rho_calc = Real(0.13) / xDs;
                Real max_s = Real(rho_i - 100.0);
                if (rho_calc > max_s) rho_calc = max_s;
                if (rho_calc < Real(50.0)) rho_calc = Real(50.0);
                rho_s = rho_calc;
            }

            Real tempc = temp[k] - Real(273.15);
            int idx_tc = static_cast<int>(std::round(-tempc));
            if (idx_tc < 1) idx_tc = 1;
            if (idx_tc > 45) idx_tc = 45;

            int idx_t = static_cast<int>((tempc - Real(2.5)) / Real(5.0)) - 1;
            idx_t = -idx_t;
            if (idx_t < 1) idx_t = 1;
            if (idx_t > ntb_t) idx_t = ntb_t;

            int idx_c = (rc[k] > tables.r_c[0]) ? tables.get_table_index(rc[k], tables.nic2, ntb_c) : 1;
            int idx_n = static_cast<int>(std::round(Real(1.0) + Real(nbc) * std::log(nc[k] / tables.t_Nc_0) / Real(tables.nic1)));
            if (idx_n < 1) idx_n = 1;
            if (idx_n > nbc) idx_n = nbc;

            int idx_i = (ri[k] > tables.r_i[0]) ? tables.get_table_index(ri[k], tables.nii2, ntb_i) : 1;
            int idx_i1 = (ni_loc[k] > tables.Nt_i[0]) ? tables.get_table_index(ni_loc[k], tables.nii3, ntb_i1) : 1;

            int idx_r = 1;
            int idx_r1 = ntb_r1;
            if (rr[k] > tables.r_r[0]) {
                idx_r = tables.get_table_index(rr[k], tables.nir2, ntb_r);
                Real lamr = Real(1.0) / ilamr[k];
                Real lam_exp = lamr * std::pow(tables.crg[2] * tables.org2 * tables.org1, Real(bm_r));
                Real N0_exp = tables.org1 * rr[k] / Real(am_r) * std::pow(lam_exp, tables.cre[0]);
                idx_r1 = tables.get_table_index(N0_exp, tables.nir3, ntb_r1);
            }

            int idx_s = (rs[k] > tables.r_s[0]) ? tables.get_table_index(rs[k], tables.nis2, ntb_s) : 1;

            int idx_g = 1;
            int idx_g1 = ntb_g1;
            if (rg[k] > tables.r_g[0]) {
                idx_g = tables.get_table_index(rg[k], tables.nig2, ntb_g);
                Real lamg = Real(1.0) / ilamg[k];
                Real lam_exp = lamg * std::pow(tables.cgg[2] * tables.ogg2 * tables.ogg1, Real(bm_g));
                Real N0_exp = tables.ogg1 * rg[k] / Real(am_g5) * std::pow(lam_exp, tables.cge[0]);
                idx_g1 = tables.get_table_index(N0_exp, tables.nig3, ntb_g1);
            }

            Real otemp = Real(1.0) / temp[k];
            Real rvs = rho[k] * qvsi[k];
            Real rvs_p = rvs * otemp * (Real(lsub) * otemp * Real(oRv_val) - Real(1.0));
            Real rvs_pp = rvs * (otemp * (Real(lsub) * otemp * Real(oRv_val) - Real(1.0))
                               * otemp * (Real(lsub) * otemp * Real(oRv_val) - Real(1.0))
                               + (-Real(2.0) * Real(lsub) * otemp * otemp * otemp * Real(oRv_val))
                               + otemp * otemp);
            Real gamsc = Real(lsub) * diffu[k] / tcond[k] * rvs_p;
            Real alphsc = Real(0.5) * (gamsc / (Real(1.0) + gamsc)) * (gamsc / (Real(1.0) + gamsc))
                        * rvs_pp / rvs_p * rvs / rvs_p;
            if (alphsc < Real(1.0e-9)) alphsc = Real(1.0e-9);
            Real xsat = ssati[k];
            if (std::abs(xsat) < Real(1.0e-9)) xsat = Real(0.0);
            Real t1_subl = Real(4.0 * PI_CONST) * (Real(1.0) - alphsc * xsat
                         + Real(2.0) * alphsc * alphsc * xsat * xsat
                         - Real(5.0) * alphsc * alphsc * alphsc * xsat * xsat * xsat)
                         / (Real(1.0) + gamsc);

            if (L_qc[k] && mvd_c[k] > Real(D0c)) {
                if (xDs > Real(D0s)) {
                    int idx = 1 + static_cast<int>(Real(nbs) * std::log(xDs / tables.Ds_0) / tables.log_Ds);
                    if (idx > nbs) idx = nbs;
                    if (idx < 1) idx = 1;
                    int j_c = static_cast<int>(mvd_c[k] * Real(1.0e6));
                    if (j_c > nbc) j_c = nbc;
                    if (j_c < 1) j_c = 1;

                    Real Ef_sw = tables.t_Efsw[tables.idx_efsw(idx - 1, j_c - 1)];
                    Real rate = rhof[k] * tables.t1_qs_qc * Ef_sw * rc[k] * smoe[k];
                    prs_scw[k] = (rate < rc[k] * odts) ? rate : rc[k] * odts;
                    Real rate_nc = rhof[k] * tables.t1_qs_qc * Ef_sw * nc[k] * smoe[k];
                    pnc_scw[k] = (rate_nc < nc[k] * odts) ? rate_nc : nc[k] * odts;
                }

                if (rg[k] >= tables.r_g[0] && mvd_c[k] > Real(D0c)) {
                    Real xDg = (Real(bm_g) + Real(mu_g) + Real(1.0)) * ilamg[k];
                    Real vtg = rhof[k] * Real(av_g_old) * tables.cgg[5] * tables.ogg3 * std::pow(ilamg[k], Real(bv_g_old));
                    Real stoke_g = mvd_c[k] * mvd_c[k] * vtg * Real(rho_w) / (Real(9.0) * visco[k] * xDg);
                    Real tc_min = (tempc < Real(-0.1)) ? tempc : Real(-0.1);
                    Real const_Ri = -Real(1.0) * (mvd_c[k] * Real(0.5e6)) * vtg / tc_min;
                    if (const_Ri < Real(0.1)) const_Ri = Real(0.1);
                    if (const_Ri > Real(10.0)) const_Ri = Real(10.0);
                    Real rime_dens = (Real(0.051) + Real(0.114) * const_Ri - Real(0.0055) * const_Ri * const_Ri) * Real(1000.0);

                    Real Ef_gw = Real(0.0);
                    if (stoke_g >= Real(0.4) && stoke_g <= Real(10.0)) {
                        Ef_gw = Real(0.55) * std::log10(Real(2.51) * stoke_g);
                    } else if (stoke_g > Real(10.0)) {
                        Ef_gw = Real(0.77);
                    }
                    if (twet[k] > Real(273.15)) Ef_gw *= Real(0.1);

                    prg_gcw[k] = rhof[k] * tables.t1_qg_qc * Ef_gw * rc[k] * N0_g[k] * std::pow(ilamg[k], tables.cge[8]);
                    Real rate_nc = rhof[k] * tables.t1_qg_qc * Ef_gw * nc[k] * N0_g[k] * std::pow(ilamg[k], tables.cge[8]);
                    pnc_gcw[k] = (rate_nc < nc[k] * odts) ? rate_nc : nc[k] * odts;
                    if (twet[k] < Real(273.15)) pbg_gcw[k] = prg_gcw[k] / rime_dens;
                }
            }

            // Rain collecting snow
            if (rr[k] >= tables.r_r[0] && rs[k] >= tables.r_s[0]) {
                size_t racs_idx = tables.idx_racs(idx_s - 1, idx_t - 1, idx_r1 - 1, idx_r - 1);
                if (twet[k] < Real(273.15)) {
                    prr_rcs[k] = -(tables.tmr_racs2[racs_idx] + tables.tcr_sacr2[racs_idx]
                                 + tables.tmr_racs1[racs_idx] + tables.tcr_sacr1[racs_idx]);
                    prs_rcs[k] = tables.tmr_racs2[racs_idx] + tables.tcr_sacr2[racs_idx]
                               - tables.tcs_racs1[racs_idx] - tables.tms_sacr1[racs_idx];
                    prg_rcs[k] = tables.tmr_racs1[racs_idx] + tables.tcr_sacr1[racs_idx]
                               + tables.tcs_racs1[racs_idx] + tables.tms_sacr1[racs_idx];

                    if (prr_rcs[k] < -rr[k] * odts) prr_rcs[k] = -rr[k] * odts;
                    if (prs_rcs[k] < -rs[k] * odts) prs_rcs[k] = -rs[k] * odts;
                    if (prg_rcs[k] > (rr[k] + rs[k]) * odts) prg_rcs[k] = (rr[k] + rs[k]) * odts;

                    pnr_rcs[k] = tables.tnr_racs1[racs_idx] + tables.tnr_racs2[racs_idx]
                               + tables.tnr_sacr1[racs_idx] + tables.tnr_sacr2[racs_idx];
                    if (pnr_rcs[k] > nr_loc[k] * odts) pnr_rcs[k] = nr_loc[k] * odts;
                    png_rcs[k] = pnr_rcs[k];
                    pbg_rcs[k] = prg_rcs[k] / Real(rho_i);
                } else {
                    prs_rcs[k] = -tables.tcs_racs1[racs_idx] - tables.tms_sacr1[racs_idx]
                                 + tables.tmr_racs2[racs_idx] + tables.tcr_sacr2[racs_idx];
                    if (prs_rcs[k] < -rs[k] * odts) prs_rcs[k] = -rs[k] * odts;
                    prr_rcs[k] = -prs_rcs[k];
                }
            }

            // Rain collecting graupel
            if (rr[k] >= tables.r_r[0] && rg[k] >= tables.r_g[0]) {
                // WRF (not hail-aware) indexes the (ntb_g1, ntb_g, dimNRHG = 1, ntb_r1, ntb_r) tables with
                // idx_bg = idx_bg1 = 5: past the size-1 dimension, i.e. 4 steps further along idx_r1 (an
                // out-of-bounds read of WRF's; kept for identical results, clamped at the end of the table)
                size_t racg_idx = tables.idx_racg(idx_g1 - 1, idx_g - 1, 4, idx_r1 - 1, idx_r - 1);
                if (racg_idx >= N_RACG) racg_idx = N_RACG - 1;
                if (twet[k] < Real(273.15)) {
                    prg_rcg[k] = tables.tmr_racg[racg_idx] + tables.tcr_gacr[racg_idx];
                    if (prg_rcg[k] > rr[k] * odts) prg_rcg[k] = rr[k] * odts;
                    prr_rcg[k] = -prg_rcg[k];
                    pnr_rcg[k] = tables.tnr_racg[racg_idx] + tables.tnr_gacr[racg_idx];
                    if (pnr_rcg[k] > nr_loc[k] * odts) pnr_rcg[k] = nr_loc[k] * odts;
                    pbg_rcg[k] = prg_rcg[k] / Real(rho_i);
                } else {
                    prr_rcg[k] = tables.tcg_racg[racg_idx];
                    if (prr_rcg[k] > rg[k] * odts) prr_rcg[k] = rg[k] * odts;
                    prg_rcg[k] = -prr_rcg[k];
                    png_rcg[k] = tables.tnr_racg[racg_idx];
                    if (png_rcg[k] > ng[k] * odts) png_rcg[k] = ng[k] * odts;
                    pbg_rcg[k] = prg_rcg[k] / Real(rho_g_val);
                    pnr_rcg[k] = Real(-1.5) * tables.tnr_gacr[racg_idx];
                }
            }

            if (temp[k] < Real(273.15)) {
                vts_boost[k] = Real(1.0);
                Real rate_max = (qv_loc[k] - qvsi[k]) * rho[k] * odts * Real(0.999);

                Real xni = Real(1000.0); // 1.0 per Liter default
                int idx_IN = (xni > tables.Nt_IN[0]) ? tables.get_table_index(xni, tables.niIN2, ntb_IN) : 1;

                if (rr[k] > tables.r_r[0]) {
                    size_t qrfz_idx = tables.idx_qrfz(idx_r - 1, idx_r1 - 1, idx_tc - 1, idx_IN - 1);
                    prg_rfz[k] = tables.tpg_qrfz[qrfz_idx] * odts;
                    pri_rfz[k] = tables.tpi_qrfz[qrfz_idx] * odts;
                    pni_rfz[k] = tables.tni_qrfz[qrfz_idx] * odts;
                    pnr_rfz[k] = tables.tnr_qrfz[qrfz_idx] * odts;
                    if (pnr_rfz[k] > nr_loc[k] * odts) pnr_rfz[k] = nr_loc[k] * odts;
                } else if (rr[k] > Real(R1) && temp[k] < Real(HGFR)) {
                    pri_rfz[k] = rr[k] * odts;
                    pni_rfz[k] = nr_loc[k] * odts;
                }
                pbg_rfz[k] = prg_rfz[k] / Real(rho_i);

                if (rc[k] > tables.r_c[0]) {
                    size_t qcfz_idx = tables.idx_qcfz(idx_c - 1, idx_n - 1, idx_tc - 1, idx_IN - 1);
                    pri_wfz[k] = tables.tpi_qcfz[qcfz_idx] * odts;
                    if (pri_wfz[k] > rc[k] * odts) pri_wfz[k] = rc[k] * odts;
                    Real pni_tbl = tables.tni_qcfz[qcfz_idx] * odts;
                    Real pni_lim = pri_wfz[k] / (Real(2.0 * xm0i));
                    Real pni_eff = (pni_lim < pni_tbl) ? pni_lim : pni_tbl;
                    pni_wfz[k] = (nc[k] * odts < pni_eff) ? nc[k] * odts : pni_eff;
                } else if (rc[k] > Real(R1) && temp[k] < Real(HGFR)) {
                    pri_wfz[k] = rc[k] * odts;
                    pni_wfz[k] = nc[k] * odts;
                }

                if (ssati[k] >= Real(0.25) || (ssatw[k] > Real(eps_val) && temp[k] < Real(253.15))) {
                    Real xnc = Real(TNO) * std::exp(Real(ATO) * (Real(273.15) - temp[k]));
                    if (xnc > Real(250.0e3)) xnc = Real(250.0e3);
                    Real xni_cur = ni_loc[k] + (pni_rfz[k] + pni_wfz[k]) * dtsave;
                    Real diff_ni = xnc - xni_cur;
                    pni_inu[k] = Real(0.5) * (diff_ni + std::abs(diff_ni)) * odts;
                    Real rate_inu = Real(xm0i) * pni_inu[k];
                    pri_inu[k] = (rate_max < rate_inu) ? rate_max : rate_inu;
                    pni_inu[k] = pri_inu[k] / Real(xm0i);
                }

                if (L_qi[k]) {
                    Real lami = std::pow(Real(am_i) * tables.cig[1] * tables.oig1 * ni_loc[k] / ri[k], obmi);
                    Real ilami = Real(1.0) / lami;
                    Real xDi = (Real(bm_i) + Real(mu_i) + Real(1.0)) * ilami;
                    if (xDi < tables.D0i) xDi = tables.D0i;
                    Real xmi = Real(am_i) * std::pow(xDi, Real(bm_i));
                    Real oxmi = Real(1.0) / xmi;

                    pri_ide[k] = Real(C_cube) * t1_subl * diffu[k] * ssati[k] * rvs
                               * tables.oig1 * tables.cig[4] * ni_loc[k] * ilami;
                    if (pri_ide[k] < Real(0.0)) {
                        Real min_ri = -ri[k] * odts;
                        if (pri_ide[k] < min_ri) pri_ide[k] = min_ri;
                        if (pri_ide[k] < rate_max) pri_ide[k] = rate_max;
                        pni_ide[k] = pri_ide[k] * oxmi;
                        Real min_ni = -ni_loc[k] * odts;
                        if (pni_ide[k] < min_ni) pni_ide[k] = min_ni;
                    } else {
                        if (pri_ide[k] > rate_max) pri_ide[k] = rate_max;
                        size_t iaus_idx = tables.idx_iaus(idx_i - 1, idx_i1 - 1);
                        Real factor = tables.tpi_ide[iaus_idx];
                        prs_ide[k] = (Real(1.0) - factor) * pri_ide[k];
                        pri_ide[k] = factor * pri_ide[k];
                    }

                    if (idx_i == ntb_i || xDi > Real(5.0 * D0s)) {
                        prs_iau[k] = ri[k] * Real(0.99) * odts;
                        pni_iau[k] = ni_loc[k] * Real(0.95) * odts;
                    } else if (xDi < Real(0.1 * D0s)) {
                        prs_iau[k] = Real(0.0);
                        pni_iau[k] = Real(0.0);
                    } else {
                        size_t iaus_idx = tables.idx_iaus(idx_i - 1, idx_i1 - 1);
                        prs_iau[k] = tables.tps_iaus[iaus_idx] * odts;
                        if (prs_iau[k] > ri[k] * Real(0.99) * odts) prs_iau[k] = ri[k] * Real(0.99) * odts;
                        pni_iau[k] = tables.tni_iaus[iaus_idx] * odts;
                        if (pni_iau[k] > ni_loc[k] * Real(0.95) * odts) pni_iau[k] = ni_loc[k] * Real(0.95) * odts;
                    }
                }

                if (L_qs[k]) {
                    Real C_snow = Real(C_sqrd) + (tempc + Real(1.5)) * Real(C_cube - C_sqrd) / Real(-30.0 + 1.5);
                    if (C_snow < Real(C_sqrd)) C_snow = Real(C_sqrd);
                    if (C_snow > Real(C_cube)) C_snow = Real(C_cube);
                    prs_sde[k] = C_snow * t1_subl * diffu[k] * ssati[k] * rvs
                               * (tables.t1_qs_sd * smo1[k] + tables.t2_qs_sd * rhof2[k] * vsc2[k] * smof[k]);
                    if (prs_sde[k] < Real(0.0)) {
                        Real min_rs = -rs[k] * odts;
                        if (prs_sde[k] < min_rs) prs_sde[k] = min_rs;
                        if (prs_sde[k] < rate_max) prs_sde[k] = rate_max;
                    } else {
                        if (prs_sde[k] > rate_max) prs_sde[k] = rate_max;
                    }
                }

                if (L_qg[k] && ssati[k] < -Real(eps_val)) {
                    prg_gde[k] = Real(C_cube) * t1_subl * diffu[k] * ssati[k] * rvs * N0_g[k]
                               * (tables.t1_qg_sd * std::pow(ilamg[k], tables.cge[9])
                                + tables.t2_qg_sd * vsc2[k] * rhof2[k] * std::pow(ilamg[k], tables.cge[10]));
                    if (prg_gde[k] < Real(0.0)) {
                        Real min_rg = -rg[k] * odts;
                        if (prg_gde[k] < min_rg) prg_gde[k] = min_rg;
                        if (prg_gde[k] < rate_max) prg_gde[k] = rate_max;
                        png_gde[k] = prg_gde[k] * ng[k] / rg[k];
                    } else {
                        if (prg_gde[k] > rate_max) prg_gde[k] = rate_max;
                    }
                }

                if (L_qi[k]) {
                    Real lami = std::pow(Real(am_i) * tables.cig[1] * tables.oig1 * ni_loc[k] / ri[k], obmi);
                    Real ilami = Real(1.0) / lami;
                    Real xDi = (Real(bm_i) + Real(mu_i) + Real(1.0)) * ilami;
                    if (xDi < tables.D0i) xDi = tables.D0i;
                    Real xmi = Real(am_i) * std::pow(xDi, Real(bm_i));
                    Real oxmi = Real(1.0) / xmi;

                    if (rs[k] >= tables.r_s[0]) {
                        prs_sci[k] = tables.t1_qs_qi * rhof[k] * Real(Ef_si) * ri[k] * smoe[k];
                        pni_sci[k] = prs_sci[k] * oxmi;
                    }

                    if (rr[k] >= tables.r_r[0] && mvd_r[k] > Real(4.0) * xDi) {
                        Real lamr = Real(1.0) / ilamr[k];
                        pri_rci[k] = rhof[k] * tables.t1_qr_qi * Real(Ef_ri) * ri[k] * N0_r[k]
                                   * std::pow(lamr + Real(fv_r), -tables.cre[8]);
                        pnr_rci[k] = rhof[k] * tables.t1_qr_qi * Real(Ef_ri) * ni_loc[k] * N0_r[k]
                                   * std::pow(lamr + Real(fv_r), -tables.cre[8]);
                        if (pnr_rci[k] > nr_loc[k] * odts) pnr_rci[k] = nr_loc[k] * odts;
                        pni_rci[k] = pri_rci[k] * oxmi;
                        prr_rci[k] = rhof[k] * tables.t2_qr_qi * Real(Ef_ri) * ni_loc[k] * N0_r[k]
                                   * std::pow(lamr + Real(fv_r), -tables.cre[7]);
                        if (prr_rci[k] > rr[k] * odts) prr_rci[k] = rr[k] * odts;
                        prg_rci[k] = pri_rci[k] + prr_rci[k];
                        pbg_rci[k] = prg_rci[k] / Real(rho_i);
                    }
                }

                if (prg_gcw[k] > Real(eps_val) && tempc > Real(-8.0)) {
                    Real tf = Real(0.0);
                    if (tempc >= Real(-5.0) && tempc < Real(-3.0)) {
                        tf = Real(0.5) * (Real(-3.0) - tempc);
                    } else if (tempc > Real(-8.0) && tempc < Real(-5.0)) {
                        tf = Real(0.33333333) * (Real(8.0) + tempc);
                    }
                    pni_ihm[k] = Real(3.5e8) * tf * prg_gcw[k];
                    pri_ihm[k] = Real(xm0i) * pni_ihm[k];
                    prs_ihm[k] = prs_scw[k] / (prs_scw[k] + prg_gcw[k]) * pri_ihm[k];
                    prg_ihm[k] = prg_gcw[k] / (prs_scw[k] + prg_gcw[k]) * pri_ihm[k];
                }

                if (prs_scw[k] > Real(2.0) * prs_sde[k] && prs_sde[k] > Real(eps_val)) {
                    Real r_frac = prs_scw[k] / prs_sde[k];
                    if (r_frac > Real(30.0)) r_frac = Real(30.0);
                    Real g_frac = Real(0.15) + (r_frac - Real(2.0)) * Real(0.028);
                    if (g_frac > Real(0.95)) g_frac = Real(0.95);
                    Real v_boost = Real(1.1) + (r_frac - Real(2.0)) * Real(0.014);
                    if (v_boost > Real(1.5)) v_boost = Real(1.5);
                    vts_boost[k] = v_boost;
                    prg_scw[k] = g_frac * prs_scw[k];
                    png_scw[k] = prg_scw[k] * smo0[k] / rs[k];

                    Real vts = Real(av_s) * std::pow(xDs, Real(bv_s)) * std::exp(-Real(fv_s) * xDs);
                    Real tc_min = (tempc < Real(-0.1)) ? tempc : Real(-0.1);
                    Real const_Ri = -Real(1.0) * (mvd_c[k] * Real(0.5e6)) * vts / tc_min;
                    if (const_Ri < Real(0.1)) const_Ri = Real(0.1);
                    if (const_Ri > Real(10.0)) const_Ri = Real(10.0);
                    Real rime_dens = (Real(0.051) + Real(0.114) * const_Ri - Real(0.0055) * const_Ri * const_Ri) * Real(1000.0);
                    if (rime_dens < Real(150.0)) {
                        g_frac = Real(0.0);
                        prg_scw[k] = Real(0.0);
                        png_scw[k] = Real(0.0);
                    }
                    pbg_scw[k] = prg_scw[k] / (Real(0.5) * (rime_dens + rho_s));
                    prs_scw[k] = (Real(1.0) - g_frac) * prs_scw[k];
                }
            } else {
                // Melting and above 0C
                if (L_qs[k]) {
                    prr_sml[k] = (tempc * tcond[k] - Real(lvap0) * diffu[k] * delQvs[k])
                               * (tables.t1_qs_me * smo1[k] + tables.t2_qs_me * rhof2[k] * vsc2[k] * smof[k]);
                    if (prr_sml[k] > Real(0.0)) {
                        prr_sml[k] += Real(4218.0 * olfus) * (twet[k] - Real(273.15)) * (prr_rcs[k] + prs_scw[k]);
                    }
                    if (prr_sml[k] < Real(0.0)) prr_sml[k] = Real(0.0);
                    if (prr_sml[k] > rs[k] * odts) prr_sml[k] = rs[k] * odts;

                    if (prr_sml[k] > Real(0.0)) {
                        pnr_sml[k] = smo0[k] / rs[k] * prr_sml[k] * std::pow(Real(10.0), Real(-0.25) * (twet[k] - Real(273.15)));
                    } else if (ssati[k] < Real(0.0)) {
                        prs_sde[k] = Real(C_sqrd) * t1_subl * diffu[k] * ssati[k] * rvs
                                   * (tables.t1_qs_sd * smo1[k] + tables.t2_qs_sd * rhof2[k] * vsc2[k] * smof[k]);
                        if (prs_sde[k] < -rs[k] * odts) prs_sde[k] = -rs[k] * odts;
                    }
                }

                if (L_qg[k]) {
                    Real N0_melt = N0_g[k];
                    if (rg[k] * ng[k] < Real(1.0e-4)) {
                        Real lamg = Real(1.0) / ilamg[k];
                        N0_melt = (Real(1.0e-4) / rg[k]) * tables.ogg2 * std::pow(lamg, tables.cge[1]);
                    }
                    prr_gml[k] = (tempc * tcond[k] - Real(lvap0) * diffu[k] * delQvs[k])
                               * N0_melt * (tables.t1_qg_me * std::pow(ilamg[k], tables.cge[9])
                                          + tables.t2_qg_me * rhof2[k] * vsc2[k] * std::pow(ilamg[k], tables.cge[10]));
                    if (prr_gml[k] < Real(0.0)) prr_gml[k] = Real(0.0);
                    if (prr_gml[k] > rg[k] * odts) prr_gml[k] = rg[k] * odts;

                    if (prr_gml[k] > Real(0.0)) {
                        Real melt_f = prr_gml[k] * dt / rg[k];
                        if (melt_f < Real(0.05)) melt_f = Real(0.05);
                        if (melt_f > Real(1.0)) melt_f = Real(1.0);
                        Real dens = melt_f * Real(rho_g_val);
                        if (dens < Real(50.0)) dens = Real(50.0);
                        if (dens > Real(1000.0)) dens = Real(1000.0);
                        pbg_gml[k] = prr_gml[k] / dens;
                        pnr_gml[k] = prr_gml[k] * ng[k] / rg[k] * std::pow(Real(10.0), Real(-0.33) * (twet[k] - Real(273.15)));
                    } else if (ssati[k] < Real(0.0)) {
                        prg_gde[k] = Real(C_cube) * t1_subl * diffu[k] * ssati[k] * rvs
                                   * N0_g[k] * (tables.t1_qg_sd * std::pow(ilamg[k], tables.cge[9])
                                              + tables.t2_qg_sd * vsc2[k] * rhof2[k] * std::pow(ilamg[k], tables.cge[10]));
                        if (prg_gde[k] < -rg[k] * odts) prg_gde[k] = -rg[k] * odts;
                        png_gde[k] = prg_gde[k] * ng[k] / rg[k];
                    }
                }

                if (dt > Real(120.0)) {
                    prr_rcw[k] += prs_scw[k] + prg_gcw[k];
                    prs_scw[k] = Real(0.0);
                    prg_gcw[k] = Real(0.0);
                }
            }
        }

        // Conservation
        for (int k = 0; k < nz; ++k) {
            Real sump = pri_inu[k] + pri_ide[k] + prs_ide[k] + prs_sde[k] + prg_gde[k] + pri_iha[k];
            Real rate_max = (qv_loc[k] - qvsi[k]) * rho[k] * odts * Real(0.999);
            if ((sump > Real(eps_val) && sump > rate_max) || (sump < -Real(eps_val) && sump < rate_max)) {
                Real ratio = rate_max / sump;
                pri_inu[k] *= ratio;
                pri_ide[k] *= ratio;
                pni_ide[k] *= ratio;
                prs_ide[k] *= ratio;
                prs_sde[k] *= ratio;
                prg_gde[k] *= ratio;
                pri_iha[k] *= ratio;
            }

            sump = -prr_wau[k] - pri_wfz[k] - prr_rcw[k] - prs_scw[k] - prg_scw[k] - prg_gcw[k];
            rate_max = -rc[k] * odts;
            if (sump < rate_max && L_qc[k]) {
                Real ratio = rate_max / sump;
                prr_wau[k] *= ratio;
                pri_wfz[k] *= ratio;
                prr_rcw[k] *= ratio;
                prs_scw[k] *= ratio;
                prg_scw[k] *= ratio;
                prg_gcw[k] *= ratio;
            }

            sump = pri_ide[k] - prs_iau[k] - prs_sci[k] - pri_rci[k];
            rate_max = -ri[k] * odts;
            if (sump < rate_max && L_qi[k]) {
                Real ratio = rate_max / sump;
                pri_ide[k] *= ratio;
                prs_iau[k] *= ratio;
                prs_sci[k] *= ratio;
                pri_rci[k] *= ratio;
            }

            sump = -prg_rfz[k] - pri_rfz[k] - prr_rci[k] + prr_rcs[k] + prr_rcg[k];
            rate_max = -rr[k] * odts;
            if (sump < rate_max && L_qr[k]) {
                Real ratio = rate_max / sump;
                prg_rfz[k] *= ratio;
                pri_rfz[k] *= ratio;
                prr_rci[k] *= ratio;
                prr_rcs[k] *= ratio;
                prr_rcg[k] *= ratio;
            }

            sump = prs_sde[k] - prs_ihm[k] - prr_sml[k] + prs_rcs[k];
            rate_max = -rs[k] * odts;
            if (sump < rate_max && L_qs[k]) {
                Real ratio = rate_max / sump;
                prs_sde[k] *= ratio;
                prs_ihm[k] *= ratio;
                prr_sml[k] *= ratio;
                prs_rcs[k] *= ratio;
            }

            sump = prg_gde[k] - prg_ihm[k] - prr_gml[k] + prg_rcg[k];
            rate_max = -rg[k] * odts;
            if (sump < rate_max && L_qg[k]) {
                Real ratio = rate_max / sump;
                prg_gde[k] *= ratio;
                prg_ihm[k] *= ratio;
                prr_gml[k] *= ratio;
                prg_rcg[k] *= ratio;
            }

            pri_ihm[k] = prs_ihm[k] + prg_ihm[k];
            Real ratio = (std::abs(prr_rcg[k]) < std::abs(prg_rcg[k])) ? std::abs(prr_rcg[k]) : std::abs(prg_rcg[k]);
            prr_rcg[k] = ratio * ((prr_rcg[k] >= Real(0.0)) ? Real(1.0) : Real(-1.0));
            prg_rcg[k] = -prr_rcg[k];
            if (twet[k] > Real(273.15)) {
                Real ratio_s = (std::abs(prr_rcs[k]) < std::abs(prs_rcs[k])) ? std::abs(prr_rcs[k]) : std::abs(prs_rcs[k]);
                prr_rcs[k] = ratio_s * ((prr_rcs[k] >= Real(0.0)) ? Real(1.0) : Real(-1.0));
                prs_rcs[k] = -prr_rcs[k];
            }
        }

        // Tendencies
        for (int k = 0; k < nz; ++k) {
            Real orho = Real(1.0) / rho[k];
            Real lfus2 = Real(lsub) - lvap[k];

            qvten[k] += (-pri_inu[k] - pri_iha[k] - pri_ide[k] - prs_ide[k] - prs_sde[k] - prg_gde[k]) * orho;
            qcten[k] += (-prr_wau[k] - pri_wfz[k] - prr_rcw[k] - prs_scw[k] - prg_scw[k] - prg_gcw[k]) * orho;
            ncten[k] += (-pnc_wau[k] - pnc_rcw[k] - pni_wfz[k] - pnc_scw[k] - pnc_gcw[k]) * orho;

            Real xrc = (qc1d[k] + qcten[k] * dtsave) * rho[k];
            if (xrc < Real(R1)) xrc = Real(R1);
            Real xnc = (nc1d[k] + ncten[k] * dtsave) * rho[k];
            if (xnc < Real(2.0)) xnc = Real(2.0);
            if (xrc > Real(R1)) {
                int nu_c = static_cast<int>(std::round(Real(1000.0e6) / xnc)) + 2;
                if (nu_c > 15) nu_c = 15;
                int nu_idx = nu_c - 1;
                Real lamc = std::pow(xnc * Real(am_r) * tables.ccg[1][nu_idx] * tables.ocg1[nu_idx] / rc[k], obmr);
                Real xDc = (Real(bm_r) + Real(nu_c) + Real(1.0)) / lamc;
                if (xDc < Real(D0c)) {
                    lamc = tables.cce[1][nu_idx] / Real(D0c);
                    xnc = tables.ccg[0][nu_idx] * tables.ocg2[nu_idx] * xrc / Real(am_r) * std::pow(lamc, Real(bm_r));
                    ncten[k] = (xnc - nc1d[k] * rho[k]) * odts * orho;  // WRF: (x-n1d*rho)*odts*orho
                } else if (xDc > Real(D0r * 2.0)) {
                    lamc = tables.cce[1][nu_idx] / Real(D0r * 2.0);
                    xnc = tables.ccg[0][nu_idx] * tables.ocg2[nu_idx] * xrc / Real(am_r) * std::pow(lamc, Real(bm_r));
                    ncten[k] = (xnc - nc1d[k] * rho[k]) * odts * orho;  // WRF: (x-n1d*rho)*odts*orho
                }
            } else {
                ncten[k] = -nc1d[k] * odts;
            }
            xnc = (nc1d[k] + ncten[k] * dtsave) * rho[k];
            if (xnc > Real(Nt_c_max)) ncten[k] = (Real(Nt_c_max) - nc1d[k] * rho[k]) * odts * orho;  // WRF: (x-n1d*rho)*odts*orho

            qiten[k] += (pri_inu[k] + pri_iha[k] + pri_ihm[k] + pri_wfz[k] + pri_rfz[k] + pri_ide[k]
                       - prs_iau[k] - prs_sci[k] - pri_rci[k]) * orho;
            niten[k] += (pni_inu[k] + pni_iha[k] + pni_ihm[k] + pni_wfz[k] + pni_rfz[k] + pni_ide[k]
                       - pni_iau[k] - pni_sci[k] - pni_rci[k]) * orho;

            Real xri = (qi1d[k] + qiten[k] * dtsave) * rho[k];
            if (xri < Real(R1)) xri = Real(R1);
            Real xni = (ni1d[k] + niten[k] * dtsave) * rho[k];
            if (xni < Real(R2)) xni = Real(R2);
            if (xri > Real(R1)) {
                Real lami = std::pow(Real(am_i) * tables.cig[1] * tables.oig1 * xni / xri, obmi);
                Real ilami = Real(1.0) / lami;
                Real xDi = (Real(bm_i) + Real(mu_i) + Real(1.0)) * ilami;
                if (xDi < Real(5.0e-6)) {
                    lami = tables.cie[1] / Real(5.0e-6);
                    Real term = tables.cig[0] * tables.oig2 * xri / Real(am_i) * std::pow(lami, Real(bm_i));
                    xni = (term < Real(999.0e3)) ? term : Real(999.0e3);
                    niten[k] = (xni - ni1d[k] * rho[k]) * odts * orho;  // WRF: (x-n1d*rho)*odts*orho
                } else if (xDi > Real(300.0e-6)) {
                    lami = tables.cie[1] / Real(300.0e-6);
                    xni = tables.cig[0] * tables.oig2 * xri / Real(am_i) * std::pow(lami, Real(bm_i));
                    niten[k] = (xni - ni1d[k] * rho[k]) * odts * orho;  // WRF: (x-n1d*rho)*odts*orho
                }
            } else {
                niten[k] = -ni1d[k] * odts;
            }
            xni = (ni1d[k] + niten[k] * dtsave) * rho[k];
            if (xni > Real(999.0e3)) niten[k] = (Real(999.0e3) - ni1d[k] * rho[k]) * odts * orho;  // WRF: (x-n1d*rho)*odts*orho

            qrten[k] += (prr_wau[k] + prr_rcw[k] + prr_sml[k] + prr_gml[k] + prr_rcs[k] + prr_rcg[k]
                       - prg_rfz[k] - pri_rfz[k] - prr_rci[k]) * orho;
            nrten[k] += (pnr_wau[k] + pnr_sml[k] + pnr_gml[k]
                       - (pnr_rfz[k] + pnr_rcr[k] + pnr_rcg[k] + pnr_rcs[k] + pnr_rci[k] + pni_rfz[k])) * orho;

            Real xrr = (qr1d[k] + qrten[k] * dtsave) * rho[k];
            if (xrr < Real(R1)) xrr = Real(R1);
            Real xnr = (nr1d[k] + nrten[k] * dtsave) * rho[k];
            if (xnr < Real(R2)) xnr = Real(R2);
            if (xrr > Real(R1)) {
                Real lamr = std::pow(Real(am_r) * tables.crg[2] * tables.org2 * xnr / xrr, obmr);
                mvd_r[k] = (Real(3.0) + Real(mu_r) + Real(0.672)) / lamr;
                if (mvd_r[k] > Real(2.5e-3)) {
                    mvd_r[k] = Real(2.5e-3);
                    lamr = (Real(3.0) + Real(mu_r) + Real(0.672)) / mvd_r[k];
                    xnr = tables.crg[1] * tables.org3 * xrr * std::pow(lamr, Real(bm_r)) / Real(am_r);
                    nrten[k] = (xnr - nr1d[k] * rho[k]) * odts * orho;  // WRF: (x-n1d*rho)*odts*orho
                } else if (mvd_r[k] < Real(D0r) * Real(0.75)) {
                    mvd_r[k] = Real(D0r) * Real(0.75);
                    lamr = (Real(3.0) + Real(mu_r) + Real(0.672)) / mvd_r[k];
                    xnr = tables.crg[1] * tables.org3 * xrr * std::pow(lamr, Real(bm_r)) / Real(am_r);
                    nrten[k] = (xnr - nr1d[k] * rho[k]) * odts * orho;  // WRF: (x-n1d*rho)*odts*orho
                }
            } else {
                qrten[k] = -qr1d[k] * odts;
                nrten[k] = -nr1d[k] * odts;
            }

            qsten[k] += (prs_iau[k] + prs_sde[k] + prs_sci[k] + prs_scw[k] + prs_rcs[k] + prs_ide[k]
                       - prs_ihm[k] - prr_sml[k]) * orho;
            qgten[k] += (prg_scw[k] + prg_rfz[k] + prg_gde[k] + prg_rcg[k] + prg_gcw[k] + prg_rci[k] + prg_rcs[k]
                       - prg_ihm[k] - prr_gml[k]) * orho;
            ngten[k] += (png_scw[k] + pnr_rfz[k] - png_rcg[k] + pnr_rci[k] + png_rcs[k] + png_gde[k]
                       - pnr_gml[k]) * orho;
            qbten[k] += (pbg_scw[k] + pbg_rfz[k] + pbg_gcw[k] + pbg_rci[k] + pbg_rcs[k] + pbg_rcg[k] + pbg_sml[k] - pbg_gml[k]
                       + (prg_gde[k] - prg_ihm[k]) / Real(rho_g_val)) * orho;

            Real xrg = (qg1d[k] + qgten[k] * dtsave) * rho[k];
            if (xrg < Real(R1)) xrg = Real(R1);
            Real xng = (ng1d[k] + ngten[k] * dtsave) * rho[k];
            if (xng < Real(R2)) xng = Real(R2);
            if (xrg > Real(R1)) {
                Real lamg = std::pow(Real(am_g5) * tables.cgg[2] * tables.ogg2 * xng / xrg, tables.obmg);
                mvd_g[k] = (Real(3.0) + Real(mu_g) + Real(0.672)) / lamg;
                if (mvd_g[k] > Real(25.4e-3)) {
                    mvd_g[k] = Real(25.4e-3);
                    lamg = (Real(3.0) + Real(mu_g) + Real(0.672)) / mvd_g[k];
                    xng = tables.cgg[1] * tables.ogg3 * xrg * std::pow(lamg, Real(bm_g)) / Real(am_g5);
                    ngten[k] = (xng - ng1d[k] * rho[k]) * odts * orho;  // WRF: (x-n1d*rho)*odts*orho
                } else if (mvd_g[k] < Real(D0r)) {
                    mvd_g[k] = Real(D0r);
                    lamg = (Real(3.0) + Real(mu_g) + Real(0.672)) / mvd_g[k];
                    xng = tables.cgg[1] * tables.ogg3 * xrg * std::pow(lamg, Real(bm_g)) / Real(am_g5);
                    ngten[k] = (xng - ng1d[k] * rho[k]) * odts * orho;  // WRF: (x-n1d*rho)*odts*orho
                }
            } else {
                qgten[k] = -qg1d[k] * odts;
                ngten[k] = -ng1d[k] * odts;
                qbten[k] = -qb1d[k] * odts;
            }

            if (temp[k] < Real(273.15)) {
                tten[k] += (Real(lsub) * ocp[k] * (pri_inu[k] + pri_ide[k] + prs_ide[k] + prs_sde[k] + prg_gde[k] + pri_iha[k])
                          + lfus2 * ocp[k] * (pri_wfz[k] + pri_rfz[k] + prg_rfz[k] + prs_scw[k] + prg_scw[k] + prg_gcw[k]
                                            + prg_rcs[k] + prs_rcs[k] + prr_rci[k] + prg_rcg[k])) * orho;
            } else {
                tten[k] += (Real(lfus) * ocp[k] * (-prr_sml[k] - prr_gml[k] - prr_rcg[k] - prr_rcs[k])
                          + Real(lsub) * ocp[k] * (prs_sde[k] + prg_gde[k])) * orho;
            }
        }

        // TAU+1 updates
        for (int k = 0; k < nz; ++k) {
            temp[k] = (th(i, j, k) * pii(i, j, k)) + dt * tten[k];
            Real otemp = Real(1.0) / temp[k];
            Real tempc = temp[k] - Real(273.15);
            qv_loc[k] = qv(i, j, k) + dt * qvten[k];
            if (qv_loc[k] < Real(1.0e-10)) qv_loc[k] = Real(1.0e-10);
            rho[k] = Real(0.622) * pres[k] / (Real(R_d) * temp[k] * (qv_loc[k] + Real(0.622)));
            rhof[k] = std::sqrt(Real(RHO_NOT) / rho[k]);
            rhof2[k] = std::sqrt(rhof[k]);
            qvs[k] = rslf(pres[k], temp[k]);
            ssatw[k] = qv_loc[k] / qvs[k] - Real(1.0);
            if (std::abs(ssatw[k]) < Real(eps_val)) ssatw[k] = Real(0.0);
            diffu[k] = Real(2.11e-5) * std::pow(temp[k] / Real(273.15), Real(1.94)) * (Real(101325.0) / pres[k]);
            if (tempc >= Real(0.0)) {
                visco[k] = (Real(1.718) + Real(0.0049) * tempc) * Real(1.0e-5);
            } else {
                visco[k] = (Real(1.718) + Real(0.0049) * tempc - Real(1.2e-5) * tempc * tempc) * Real(1.0e-5);
            }
            vsc2[k] = std::sqrt(rho[k] / visco[k]);
            lvap[k] = Real(lvap0) + (Real(2106.0) - Real(4218.0)) * tempc;
            tcond[k] = (Real(5.69) + Real(0.0168) * tempc) * Real(1.0e-5) * Real(418.936);  // WRF: (..)*1.0E-5 * 418.936
            ocp[k] = Real(1.0) / (Real(Cp_val) * (Real(1.0) + Real(0.887) * qv_loc[k]));
            lvt2[k] = lvap[k] * lvap[k] * ocp[k] * Real(oRv_val) * otemp * otemp;

            Real qc_new = qc1d[k] + qcten[k] * dt;
            if (qc_new > Real(R1)) {
                rc[k] = qc_new * rho[k];
                nc[k] = Real(Nt_c);
                L_qc[k] = true;
            } else {
                rc[k] = Real(R1);
                nc[k] = Real(2.0);
                L_qc[k] = false;
            }

            Real qi_new = qi1d[k] + qiten[k] * dt;
            if (qi_new > Real(R1)) {
                ri[k] = qi_new * rho[k];
                Real ni_calc = (ni1d[k] + niten[k] * dt) * rho[k];
                ni_loc[k] = (ni_calc > Real(R2)) ? ni_calc : Real(R2);
                L_qi[k] = true;
            } else {
                ri[k] = Real(R1);
                ni_loc[k] = Real(R2);
                L_qi[k] = false;
            }

            Real qr_new = qr1d[k] + qrten[k] * dt;
            if (qr_new > Real(R1)) {
                rr[k] = qr_new * rho[k];
                Real nr_calc = (nr1d[k] + nrten[k] * dt) * rho[k];
                nr_loc[k] = (nr_calc > Real(R2)) ? nr_calc : Real(R2);
                L_qr[k] = true;
                Real lamr = std::pow(Real(am_r) * tables.crg[2] * tables.org2 * nr_loc[k] / rr[k], obmr);
                mvd_r[k] = (Real(3.0) + Real(mu_r) + Real(0.672)) / lamr;
                if (mvd_r[k] > Real(2.5e-3)) {
                    mvd_r[k] = Real(2.5e-3);
                    lamr = (Real(3.0) + Real(mu_r) + Real(0.672)) / mvd_r[k];
                    nr_loc[k] = tables.crg[1] * tables.org3 * rr[k] * std::pow(lamr, Real(bm_r)) / Real(am_r);
                } else if (mvd_r[k] < Real(D0r) * Real(0.75)) {
                    mvd_r[k] = Real(D0r) * Real(0.75);
                    lamr = (Real(3.0) + Real(mu_r) + Real(0.672)) / mvd_r[k];
                    nr_loc[k] = tables.crg[1] * tables.org3 * rr[k] * std::pow(lamr, Real(bm_r)) / Real(am_r);
                }
            } else {
                rr[k] = Real(R1);
                nr_loc[k] = Real(R2);
                L_qr[k] = false;
            }

            Real qs_new = qs1d[k] + qsten[k] * dt;
            if (qs_new > Real(R1)) {
                rs[k] = qs_new * rho[k];
                L_qs[k] = true;
            } else {
                rs[k] = Real(R1);
                L_qs[k] = false;
            }

            Real qg_new = qg1d[k] + qgten[k] * dt;
            if (qg_new > Real(R1)) {
                rg[k] = qg_new * rho[k];
                Real ygra1 = std::log10((rg[k] > Real(1.0e-9)) ? rg[k] : Real(1.0e-9));
                Real zans1 = Real(3.0) + Real(2.0 / 7.0) * (ygra1 + Real(8.0));
                if (zans1 < Real(2.0)) zans1 = Real(2.0);
                if (zans1 > Real(6.0)) zans1 = Real(6.0);
                Real N0_exp = std::pow(Real(10.0), zans1);
                Real lam_exp = std::pow(N0_exp * Real(am_g5) * tables.cgg[0] / rg[k], tables.oge1);
                Real lamg = lam_exp * std::pow(tables.cgg[2] * tables.ogg2 * tables.ogg1, tables.obmg);
                ng[k] = tables.cgg[1] * tables.ogg3 * rg[k] * std::pow(lamg, Real(bm_g)) / Real(am_g5);
                rb[k] = rg[k] / rho[k] / Real(rho_g_val);
            } else {
                rg[k] = Real(R1);
                ng[k] = Real(R2);
                rb[k] = Real(R1) / rho[k] / Real(800.0);
                L_qg[k] = false;
            }
        }

        // Snow moments recalculation
        for (int k = 0; k < nz; ++k) {
            smo2[k] = Real(0.0);
            smob[k] = Real(0.0);
            smoc[k] = Real(0.0);
            smod[k] = Real(0.0);
            if (!L_qs[k]) continue;

            Real tc0 = (temp[k] - Real(273.15) < Real(-0.1)) ? temp[k] - Real(273.15) : Real(-0.1);
            smob[k] = rs[k] * tables.oams;
            smo2[k] = smob[k];

            Real cse1 = tables.cse[0];
            Real logac = tables.sa[0] + tables.sa[1]*tc0 + tables.sa[2]*cse1
                       + tables.sa[3]*tc0*cse1 + tables.sa[4]*tc0*tc0
                       + tables.sa[5]*cse1*cse1 + tables.sa[6]*tc0*tc0*cse1
                       + tables.sa[7]*tc0*cse1*cse1 + tables.sa[8]*tc0*tc0*tc0
                       + tables.sa[9]*cse1*cse1*cse1;
            Real ac = std::pow(Real(10.0), logac);
            Real bc = tables.sb[0] + tables.sb[1]*tc0 + tables.sb[2]*cse1 + tables.sb[3]*tc0*cse1
                    + tables.sb[4]*tc0*tc0 + tables.sb[5]*cse1*cse1
                    + tables.sb[6]*tc0*tc0*cse1 + tables.sb[7]*tc0*cse1*cse1
                    + tables.sb[8]*tc0*tc0*tc0 + tables.sb[9]*cse1*cse1*cse1;
            smoc[k] = ac * std::pow(smo2[k], bc);

            Real cse14 = tables.cse[13];
            Real logad = tables.sa[0] + tables.sa[1]*tc0 + tables.sa[2]*cse14
                       + tables.sa[3]*tc0*cse14 + tables.sa[4]*tc0*tc0
                       + tables.sa[5]*cse14*cse14 + tables.sa[6]*tc0*tc0*cse14
                       + tables.sa[7]*tc0*cse14*cse14 + tables.sa[8]*tc0*tc0*tc0
                       + tables.sa[9]*cse14*cse14*cse14;
            Real ad = std::pow(Real(10.0), logad);
            Real bd = tables.sb[0] + tables.sb[1]*tc0 + tables.sb[2]*cse14 + tables.sb[3]*tc0*cse14
                    + tables.sb[4]*tc0*tc0 + tables.sb[5]*cse14*cse14
                    + tables.sb[6]*tc0*tc0*cse14 + tables.sb[7]*tc0*cse14*cse14
                    + tables.sb[8]*tc0*tc0*tc0 + tables.sb[9]*cse14*cse14*cse14;
            smod[k] = ad * std::pow(smo2[k], bd);
        }

        for (int k = nz - 1; k >= 0; --k) {
            Real lamg = std::pow(Real(am_g5) * tables.cgg[2] * tables.ogg2 * ng[k] / rg[k], tables.obmg);
            ilamg[k] = Real(1.0) / lamg;
            N0_g[k] = ng[k] * tables.ogg2 * std::pow(lamg, tables.cge[1]);
        }

        // gfortran vectorises this k loop (kte down to kts): libmvec's 2-lane pow (cppwrf::vm)
        for (int k = nz - 1; k >= 0; --k) {
            Real lamr = vm::pow(Real(am_r) * tables.crg[2] * tables.org2 * nr_loc[k] / rr[k], obmr, nz - 1 - k, nz);
            ilamr[k] = Real(1.0) / lamr;
            mvd_r[k] = (Real(3.0) + Real(mu_r) + Real(0.672)) / lamr;
            N0_r[k] = nr_loc[k] * tables.org2 * vm::pow(lamr, tables.cre[1], nz - 1 - k, nz);
        }

        // Condensation & evaporation of cloud water
        for (int k = 0; k < nz; ++k) {
            Real orho = Real(1.0) / rho[k];
            if (ssatw[k] > Real(eps_val) || (ssatw[k] < -Real(eps_val) && L_qc[k])) {
                Real clap = (qv_loc[k] - qvs[k]) / (Real(1.0) + lvt2[k] * qvs[k]);
                for (int iter = 0; iter < 3; ++iter) {
                    Real fcd = qvs[k] * std::exp(lvt2[k] * clap) - qv_loc[k] + clap;
                    Real dfcd = qvs[k] * lvt2[k] * std::exp(lvt2[k] * clap) + Real(1.0);
                    clap -= fcd / dfcd;
                }
                Real xrc = rc[k] + clap * rho[k];
                if (xrc > Real(R1)) {
                    prw_vcd[k] = clap * odt;
                    if (clap > Real(eps_val)) {
                        Real xnc = Real(Nt_c);
                        Real diff_nc = xnc - nc[k];
                        pnc_wcd[k] = Real(0.5) * (diff_nc + std::abs(diff_nc)) * odts * orho;
                    }
                } else {
                    prw_vcd[k] = -rc[k] * orho * odt;
                    pnc_wcd[k] = -nc[k] * orho * odt;
                }

                qvten[k] -= prw_vcd[k];
                qcten[k] += prw_vcd[k];
                ncten[k] += pnc_wcd[k];
                tten[k] += lvap[k] * ocp[k] * prw_vcd[k];
                rc[k] = (qc1d[k] + dt * qcten[k]) * rho[k];
                if (rc[k] < Real(R1)) rc[k] = Real(R1);
                if (rc[k] == Real(R1)) L_qc[k] = false;
                nc[k] = Real(Nt_c);
                qv_loc[k] = qv(i, j, k) + dt * qvten[k];
                if (qv_loc[k] < Real(1.0e-10)) qv_loc[k] = Real(1.0e-10);
                temp[k] = (th(i, j, k) * pii(i, j, k)) + dt * tten[k];
                rho[k] = Real(0.622) * pres[k] / (Real(R_d) * temp[k] * (qv_loc[k] + Real(0.622)));
                qvs[k] = rslf(pres[k], temp[k]);
                ssatw[k] = qv_loc[k] / qvs[k] - Real(1.0);
            }
        }

        // Rain evaporation
        for (int k = 0; k < nz; ++k) {
            if (ssatw[k] < -Real(eps_val) && L_qr[k] && !(prw_vcd[k] > Real(0.0))) {
                Real tempc = temp[k] - Real(273.15);
                Real otemp = Real(1.0) / temp[k];
                Real orho = Real(1.0) / rho[k];
                rhof[k] = std::sqrt(Real(RHO_NOT) * orho);
                rhof2[k] = std::sqrt(rhof[k]);
                diffu[k] = Real(2.11e-5) * std::pow(temp[k] / Real(273.15), Real(1.94)) * (Real(101325.0) / pres[k]);
                if (tempc >= Real(0.0)) {
                    visco[k] = (Real(1.718) + Real(0.0049) * tempc) * Real(1.0e-5);
                } else {
                    visco[k] = (Real(1.718) + Real(0.0049) * tempc - Real(1.2e-5) * tempc * tempc) * Real(1.0e-5);
                }
                vsc2[k] = std::sqrt(rho[k] / visco[k]);
                lvap[k] = Real(lvap0) + (Real(2106.0) - Real(4218.0)) * tempc;
                tcond[k] = (Real(5.69) + Real(0.0168) * tempc) * Real(1.0e-5) * Real(418.936);  // WRF: (..)*1.0E-5 * 418.936
                ocp[k] = Real(1.0) / (Real(Cp_val) * (Real(1.0) + Real(0.887) * qv_loc[k]));

                Real rvs = rho[k] * qvs[k];
                Real rvs_p = rvs * otemp * (lvap[k] * otemp * Real(oRv_val) - Real(1.0));
                Real rvs_pp = rvs * (otemp * (lvap[k] * otemp * Real(oRv_val) - Real(1.0))
                                   * otemp * (lvap[k] * otemp * Real(oRv_val) - Real(1.0))
                                   + (-Real(2.0) * lvap[k] * otemp * otemp * otemp * Real(oRv_val))
                                   + otemp * otemp);
                Real gamsc = lvap[k] * diffu[k] / tcond[k] * rvs_p;
                Real alphsc = Real(0.5) * (gamsc / (Real(1.0) + gamsc)) * (gamsc / (Real(1.0) + gamsc))
                            * rvs_pp / rvs_p * rvs / rvs_p;
                if (alphsc < Real(1.0e-9)) alphsc = Real(1.0e-9);
                Real xsat = (ssatw[k] < -Real(1.0e-9)) ? ssatw[k] : -Real(1.0e-9);
                Real t1_evap = Real(2.0 * PI_CONST) * (Real(1.0) - alphsc * xsat
                             + Real(2.0) * alphsc * alphsc * xsat * xsat
                             - Real(5.0) * alphsc * alphsc * alphsc * xsat * xsat * xsat)
                             / (Real(1.0) + gamsc);

                Real lamr = Real(1.0) / ilamr[k];
                if (qv_loc[k] / qvs[k] < Real(0.95) && rr[k] * orho <= Real(1.0e-8)) {
                    prv_rev[k] = rr[k] * orho * odts;
                } else {
                    prv_rev[k] = t1_evap * diffu[k] * (-ssatw[k]) * N0_r[k] * rvs
                               * (tables.t1_qr_ev * std::pow(ilamr[k], tables.cre[9])
                                + tables.t2_qr_ev * vsc2[k] * rhof2[k] * std::pow(lamr + Real(0.5 * fv_r), -tables.cre[10]));
                    Real rate_max = (rr[k] * orho * odts < (qvs[k] - qv_loc[k]) * odts) ? rr[k] * orho * odts : (qvs[k] - qv_loc[k]) * odts;
                    Real prv_eff = prv_rev[k] * orho;
                    prv_rev[k] = (rate_max < prv_eff) ? rate_max : prv_eff;

                    if (prr_gml[k] > Real(0.0)) {
                        Real eva_factor = Real(0.01) + (Real(0.99) - Real(0.01)) * (tempc / Real(20.0));
                        if (eva_factor > Real(1.0)) eva_factor = Real(1.0);
                        prv_rev[k] *= eva_factor;
                    }
                }

                Real pnr_lim = nr_loc[k] * Real(0.99) * orho * odts;
                Real pnr_calc = prv_rev[k] * nr_loc[k] / rr[k];
                pnr_rev[k] = (pnr_lim < pnr_calc) ? pnr_lim : pnr_calc;

                qrten[k] -= prv_rev[k];
                qvten[k] += prv_rev[k];
                nrten[k] -= pnr_rev[k];
                tten[k] -= lvap[k] * ocp[k] * prv_rev[k];

                Real rr_new = (qr1d[k] + dt * qrten[k]) * rho[k];
                rr[k] = (rr_new > Real(R1)) ? rr_new : Real(R1);
                Real qv_new = qv(i, j, k) + dt * qvten[k];
                qv_loc[k] = (qv_new > Real(1.0e-10)) ? qv_new : Real(1.0e-10);
                Real nr_new = (nr1d[k] + dt * nrten[k]) * rho[k];
                nr_loc[k] = (nr_new > Real(R2)) ? nr_new : Real(R2);
                temp[k] = (th(i, j, k) * pii(i, j, k)) + dt * tten[k];
                rho[k] = Real(0.622) * pres[k] / (Real(R_d) * temp[k] * (qv_loc[k] + Real(0.622)));
            }
        }

        // Sedimentation velocities
        int nstep_r = 0, nstep_i = 0, nstep_s = 0, nstep_g = 0;
        Real onstep_r = Real(1.0), onstep_i = Real(1.0), onstep_s = Real(1.0), onstep_g = Real(1.0);
        int ksed1_r = 0, ksed1_i = 0, ksed1_s = 0, ksed1_g = 0, ksed1_c = 0;

        for (int k = 0; k <= nz; ++k) {
            vtrk[k] = Real(0.0);
            vtnrk[k] = Real(0.0);
            vtik[k] = Real(0.0);
            vtnik[k] = Real(0.0);
            vtsk[k] = Real(0.0);
            vtgk[k] = Real(0.0);
            vtngk[k] = Real(0.0);
            vtck[k] = Real(0.0);
            vtnck[k] = Real(0.0);
        }

        bool has_qr = false, has_qc = false, has_qi = false, has_qs = false, has_qg = false;
        for (int k = 0; k < nz; ++k) {
            if (L_qr[k]) has_qr = true;
            if (L_qc[k]) has_qc = true;
            if (L_qi[k]) has_qi = true;
            if (L_qs[k]) has_qs = true;
            if (L_qg[k]) has_qg = true;
        }

        // 1. Rain fallspeed
        if (has_qr) {
            for (int k = nz - 1; k >= 0; --k) {
                rhof[k] = std::sqrt(Real(RHO_NOT) / rho[k]);
                if (rr[k] > Real(R1)) {
                    Real lamr = std::pow(Real(am_r) * tables.crg[2] * tables.org2 * nr_loc[k] / rr[k], obmr);
                    // gfortran pairs the two (lamr+fv_r)**(-cre) into one libmvec pow (SLP): vector lanes
                    vtrk[k] = rhof[k] * Real(av_r) * tables.crg[5] * tables.org3 * std::pow(lamr, tables.cre[2])
                            * vm::pow(lamr + Real(fv_r), -tables.cre[5], 0, 2);
                    vtnrk[k] = rhof[k] * Real(av_r) * tables.crg[6] / tables.crg[11] * std::pow(lamr, tables.cre[11])
                             * vm::pow(lamr + Real(fv_r), -tables.cre[6], 0, 2);
                } else {
                    vtrk[k] = vtrk[k + 1];
                    vtnrk[k] = vtnrk[k + 1];
                }

                Real max_v = (vtrk[k] > vtnrk[k]) ? vtrk[k] : vtnrk[k];
                if (max_v > Real(1.0e-3)) {
                    if (k > ksed1_r) ksed1_r = k;
                    Real delta_tp = dz(i, j, k) / max_v;
                    int cur_steps = static_cast<int>(dt / delta_tp + Real(1.0));
                    if (cur_steps > nstep_r) nstep_r = cur_steps;
                }
            }
            if (ksed1_r == nz - 1) ksed1_r = nz - 2;
            if (nstep_r > 0) onstep_r = Real(1.0) / Real(nstep_r);
        }

        // 2. Cloud droplet settling
        if (has_qc) {
            Real hgt_agl = Real(0.0);
            for (int k = 0; k < nz - 1; ++k) {
                if (rc[k] > Real(R2)) ksed1_c = k;
                hgt_agl += dz(i, j, k);
                if (hgt_agl > Real(500.0)) break;
            }
            for (int k = ksed1_c; k >= 0; --k) {
                if (rc[k] > Real(R1) && w(i, j, k) < Real(0.1)) {
                    int nu_c = static_cast<int>(std::round(Real(1000.0e6) / nc[k])) + 2;
                    if (nu_c > 15) nu_c = 15;
                    int nu_idx = nu_c - 1;
                    Real lamc = std::pow(nc[k] * Real(am_r) * tables.ccg[1][nu_idx] * tables.ocg1[nu_idx] / rc[k], obmr);
                    Real ilamc = Real(1.0) / lamc;
                    vtck[k] = rhof[k] * Real(av_c) * tables.ccg[4][nu_idx] * tables.ocg2[nu_idx] * std::pow(ilamc, Real(bv_c));
                    vtnck[k] = rhof[k] * Real(av_c) * tables.ccg[3][nu_idx] * tables.ocg1[nu_idx] * std::pow(ilamc, Real(bv_c));
                }
            }
        }

        // 3. Ice fallspeed
        if (has_qi) {
            for (int k = nz - 1; k >= 0; --k) {
                if (ri[k] > Real(R1)) {
                    Real lami = std::pow(Real(am_i) * tables.cig[1] * tables.oig1 * ni_loc[k] / ri[k], obmi);
                    Real ilami = Real(1.0) / lami;
                    vtik[k] = rhof[k] * Real(av_i) * tables.cig[2] * tables.oig2 * std::pow(ilami, Real(bv_i));
                    vtnik[k] = rhof[k] * Real(av_i) * tables.cig[5] / tables.cig[6] * std::pow(ilami, Real(bv_i));
                } else {
                    vtik[k] = vtik[k + 1];
                    vtnik[k] = vtnik[k + 1];
                }

                if (vtik[k] > Real(1.0e-3)) {
                    if (k > ksed1_i) ksed1_i = k;
                    Real delta_tp = dz(i, j, k) / vtik[k];
                    int cur_steps = static_cast<int>(dt / delta_tp + Real(1.0));
                    if (cur_steps > nstep_i) nstep_i = cur_steps;
                }
            }
            if (ksed1_i == nz - 1) ksed1_i = nz - 2;
            if (nstep_i > 0) onstep_i = Real(1.0) / Real(nstep_i);
        }

        // 4. Snow fallspeed
        if (has_qs) {
            for (int k = nz - 1; k >= 0; --k) {
                if (rs[k] > Real(R1)) {
                    Real xDs = smoc[k] / smob[k];
                    Real Mrat = Real(1.0) / xDs;
                    Real ils1 = Real(1.0) / (Mrat * Real(Lam0) + Real(fv_s));
                    Real ils2 = Real(1.0) / (Mrat * Real(Lam1) + Real(fv_s));
                    Real t1_vts = Real(Kap0) * tables.csg[3] * std::pow(ils1, tables.cse[3]);
                    Real t2_vts = Real(Kap1) * std::pow(Mrat, Real(mu_s)) * tables.csg[9] * std::pow(ils2, tables.cse[9]);
                    ils1 = Real(1.0) / (Mrat * Real(Lam0));
                    ils2 = Real(1.0) / (Mrat * Real(Lam1));
                    Real t3_vts = Real(Kap0) * tables.csg[0] * std::pow(ils1, tables.cse[0]);
                    Real t4_vts = Real(Kap1) * std::pow(Mrat, Real(mu_s)) * tables.csg[6] * std::pow(ils2, tables.cse[6]);
                    Real vts = rhof[k] * Real(av_s) * (t1_vts + t2_vts) / (t3_vts + t4_vts);
                    if (prr_sml[k] > Real(0.0)) {
                        Real SR_k = rs[k] / (rs[k] + rr[k]);
                        vtsk[k] = vts * SR_k + (Real(1.0) - SR_k) * vtrk[k];
                    } else {
                        vtsk[k] = vts * vts_boost[k];
                    }
                } else {
                    vtsk[k] = vtsk[k + 1];
                }

                if (vtsk[k] > Real(1.0e-3)) {
                    if (k > ksed1_s) ksed1_s = k;
                    Real delta_tp = dz(i, j, k) / vtsk[k];
                    int cur_steps = static_cast<int>(dt / delta_tp + Real(1.0));
                    if (cur_steps > nstep_s) nstep_s = cur_steps;
                }
            }
            if (ksed1_s == nz - 1) ksed1_s = nz - 2;
            if (nstep_s > 0) onstep_s = Real(1.0) / Real(nstep_s);
        }

        // 5. Graupel fallspeed
        if (has_qg) {
            for (int k = nz - 1; k >= 0; --k) {
                if (rg[k] > Real(R1)) {
                    vtgk[k] = rhof[k] * Real(av_g_old) * tables.cgg[5] * tables.ogg3 * std::pow(ilamg[k], Real(bv_g_old));
                    vtngk[k] = rhof[k] * Real(av_g_old) * tables.cgg[6] / tables.cgg[11] * std::pow(ilamg[k], Real(bv_g_old));
                } else {
                    vtgk[k] = vtgk[k + 1];
                    vtngk[k] = vtngk[k + 1];
                }

                if (vtgk[k] > Real(1.0e-3)) {
                    if (k > ksed1_g) ksed1_g = k;
                    Real delta_tp = dz(i, j, k) / vtgk[k];
                    int cur_steps = static_cast<int>(dt / delta_tp + Real(1.0));
                    if (cur_steps > nstep_g) nstep_g = cur_steps;
                }
            }
            if (ksed1_g == nz - 1) ksed1_g = nz - 2;
            if (nstep_g > 0) onstep_g = Real(1.0) / Real(nstep_g);
        }

        // Sedimentation execution
        // Rain
        if (has_qr) {
            int nsteps = static_cast<int>(std::round(Real(1.0) / onstep_r));
            for (int n = 0; n < nsteps; ++n) {
                for (int k = nz - 1; k >= 0; --k) {
                    sed_r[k] = vtrk[k] * rr[k];
                    sed_n[k] = vtnrk[k] * nr_loc[k];
                }
                int k = nz - 1;
                Real odzq = Real(1.0) / dz(i, j, k);
                Real orho = Real(1.0) / rho[k];
                qrten[k] -= sed_r[k] * odzq * onstep_r * orho;
                nrten[k] -= sed_n[k] * odzq * onstep_r * orho;
                Real rr_sub = rr[k] - sed_r[k] * odzq * dt * onstep_r;
                rr[k] = (rr_sub > Real(R1)) ? rr_sub : Real(R1);
                Real nr_sub = nr_loc[k] - sed_n[k] * odzq * dt * onstep_r;
                nr_loc[k] = (nr_sub > Real(R2)) ? nr_sub : Real(R2);

                for (int kk = ksed1_r; kk >= 0; --kk) {
                    Real odzq_k = Real(1.0) / dz(i, j, kk);
                    Real orho_k = Real(1.0) / rho[kk];
                    Real diff_r = sed_r[kk + 1] - sed_r[kk];
                    Real diff_n = sed_n[kk + 1] - sed_n[kk];
                    qrten[kk] += diff_r * odzq_k * onstep_r * orho_k;
                    nrten[kk] += diff_n * odzq_k * onstep_r * orho_k;
                    Real rr_add = rr[kk] + diff_r * odzq_k * dt * onstep_r;
                    rr[kk] = (rr_add > Real(R1)) ? rr_add : Real(R1);
                    Real nr_add = nr_loc[kk] + diff_n * odzq_k * dt * onstep_r;
                    nr_loc[kk] = (nr_add > Real(R2)) ? nr_add : Real(R2);
                }

                if (rr[0] > Real(R1 * 1000.0)) {
                    pptrain += sed_r[0] * dt * onstep_r;
                }
            }
        }

        // Cloud water settling
        if (has_qc) {
            for (int k = nz - 1; k >= 0; --k) {
                sed_c[k] = vtck[k] * rc[k];
                sed_n[k] = vtnck[k] * nc[k];
            }
            for (int kk = ksed1_c; kk >= 0; --kk) {
                Real odzq_k = Real(1.0) / dz(i, j, kk);
                Real orho_k = Real(1.0) / rho[kk];
                Real diff_c = sed_c[kk + 1] - sed_c[kk];
                Real diff_n = sed_n[kk + 1] - sed_n[kk];
                qcten[kk] += diff_c * odzq_k * orho_k;
                ncten[kk] += diff_n * odzq_k * orho_k;
                Real rc_add = rc[kk] + diff_c * odzq_k * dt;
                rc[kk] = (rc_add > Real(R1)) ? rc_add : Real(R1);
                Real nc_add = nc[kk] + diff_n * odzq_k * dt;
                nc[kk] = (nc_add > Real(10.0)) ? nc_add : Real(10.0);
            }
        }

        // Ice
        if (has_qi) {
            int nsteps = static_cast<int>(std::round(Real(1.0) / onstep_i));
            for (int n = 0; n < nsteps; ++n) {
                for (int k = nz - 1; k >= 0; --k) {
                    sed_i[k] = vtik[k] * ri[k];
                    sed_n[k] = vtnik[k] * ni_loc[k];
                }
                int k = nz - 1;
                Real odzq = Real(1.0) / dz(i, j, k);
                Real orho = Real(1.0) / rho[k];
                qiten[k] -= sed_i[k] * odzq * onstep_i * orho;
                niten[k] -= sed_n[k] * odzq * onstep_i * orho;
                Real ri_sub = ri[k] - sed_i[k] * odzq * dt * onstep_i;
                ri[k] = (ri_sub > Real(R1)) ? ri_sub : Real(R1);
                Real ni_sub = ni_loc[k] - sed_n[k] * odzq * dt * onstep_i;
                ni_loc[k] = (ni_sub > Real(R2)) ? ni_sub : Real(R2);

                for (int kk = ksed1_i; kk >= 0; --kk) {
                    Real odzq_k = Real(1.0) / dz(i, j, kk);
                    Real orho_k = Real(1.0) / rho[kk];
                    Real diff_i = sed_i[kk + 1] - sed_i[kk];
                    Real diff_n = sed_n[kk + 1] - sed_n[kk];
                    qiten[kk] += diff_i * odzq_k * onstep_i * orho_k;
                    niten[kk] += diff_n * odzq_k * onstep_i * orho_k;
                    Real ri_add = ri[kk] + diff_i * odzq_k * dt * onstep_i;
                    ri[kk] = (ri_add > Real(R1)) ? ri_add : Real(R1);
                    Real ni_add = ni_loc[kk] + diff_n * odzq_k * dt * onstep_i;
                    ni_loc[kk] = (ni_add > Real(R2)) ? ni_add : Real(R2);
                }

                if (ri[0] > Real(R1 * 1000.0)) {
                    pptice += sed_i[0] * dt * onstep_i;
                }
            }
        }

        // Snow
        if (has_qs) {
            int nsteps = static_cast<int>(std::round(Real(1.0) / onstep_s));
            for (int n = 0; n < nsteps; ++n) {
                for (int k = nz - 1; k >= 0; --k) {
                    sed_s[k] = vtsk[k] * rs[k];
                }
                int k = nz - 1;
                Real odzq = Real(1.0) / dz(i, j, k);
                Real orho = Real(1.0) / rho[k];
                qsten[k] -= sed_s[k] * odzq * onstep_s * orho;
                Real rs_sub = rs[k] - sed_s[k] * odzq * dt * onstep_s;
                rs[k] = (rs_sub > Real(R1)) ? rs_sub : Real(R1);

                for (int kk = ksed1_s; kk >= 0; --kk) {
                    Real odzq_k = Real(1.0) / dz(i, j, kk);
                    Real orho_k = Real(1.0) / rho[kk];
                    Real diff_s = sed_s[kk + 1] - sed_s[kk];
                    qsten[kk] += diff_s * odzq_k * onstep_s * orho_k;
                    Real rs_add = rs[kk] + diff_s * odzq_k * dt * onstep_s;
                    rs[kk] = (rs_add > Real(R1)) ? rs_add : Real(R1);
                }

                if (rs[0] > Real(R1 * 1000.0)) {
                    pptsnow += sed_s[0] * dt * onstep_s;
                }
            }
        }

        // Graupel
        if (has_qg) {
            int nsteps = static_cast<int>(std::round(Real(1.0) / onstep_g));
            for (int n = 0; n < nsteps; ++n) {
                for (int k = nz - 1; k >= 0; --k) {
                    sed_g[k] = vtgk[k] * rg[k];
                    sed_n[k] = vtngk[k] * ng[k];
                    sed_b[k] = vtgk[k] * rb[k];
                }
                int k = nz - 1;
                Real odzq = Real(1.0) / dz(i, j, k);
                Real orho = Real(1.0) / rho[k];
                qgten[k] -= sed_g[k] * odzq * onstep_g * orho;
                ngten[k] -= sed_n[k] * odzq * onstep_g * orho;
                qbten[k] -= sed_b[k] * odzq * onstep_g;
                Real rg_sub = rg[k] - sed_g[k] * odzq * dt * onstep_g;
                rg[k] = (rg_sub > Real(R1)) ? rg_sub : Real(R1);
                Real ng_sub = ng[k] - sed_n[k] * odzq * dt * onstep_g;
                ng[k] = (ng_sub > Real(R2)) ? ng_sub : Real(R2);
                Real rb_sub = rb[k] - sed_b[k] * odzq * dt * onstep_g;
                Real rb_min = Real(R1) / (rho[k] * Real(800.0));
                rb[k] = (rb_sub > rb_min) ? rb_sub : rb_min;

                for (int kk = ksed1_g; kk >= 0; --kk) {
                    Real odzq_k = Real(1.0) / dz(i, j, kk);
                    Real orho_k = Real(1.0) / rho[kk];
                    Real diff_g = sed_g[kk + 1] - sed_g[kk];
                    Real diff_n = sed_n[kk + 1] - sed_n[kk];
                    Real diff_b = sed_b[kk + 1] - sed_b[kk];
                    qgten[kk] += diff_g * odzq_k * onstep_g * orho_k;
                    ngten[kk] += diff_n * odzq_k * onstep_g * orho_k;
                    qbten[kk] += diff_b * odzq_k * onstep_g;
                    Real rg_add = rg[kk] + diff_g * odzq_k * dt * onstep_g;
                    rg[kk] = (rg_add > Real(R1)) ? rg_add : Real(R1);
                    Real ng_add = ng[kk] + diff_n * odzq_k * dt * onstep_g;
                    ng[kk] = (ng_add > Real(R2)) ? ng_add : Real(R2);
                    Real rb_add = rb[kk] + diff_b * odzq_k * dt * onstep_g;
                    Real rb_floor = rg[kk] / (rho[kk] * Real(800.0));
                    rb[kk] = (rb_add > rb_floor) ? rb_add : rb_floor;
                }

                if (rg[0] > Real(R1 * 1000.0)) {
                    pptgraul += sed_g[0] * dt * onstep_g;
                }
            }
        }

        // Instantaneous melting of cloud ice > 0C & freezing of cloud water < HGFR
        for (int k = 0; k < nz; ++k) {
            Real xri = qi1d[k] + qiten[k] * dt;
            if (xri < Real(0.0)) xri = Real(0.0);
            if (temp[k] > Real(273.15) && xri > Real(0.0)) {
                qcten[k] += xri * odt;
                ncten[k] += ni1d[k] * odt;
                qiten[k] -= xri * odt;
                niten[k] = -ni1d[k] * odt;
                tten[k] -= Real(lfus) * ocp[k] * xri * odt;
            }

            Real xrc = qc1d[k] + qcten[k] * dt;
            if (xrc < Real(0.0)) xrc = Real(0.0);
            if (temp[k] < Real(HGFR) && xrc > Real(0.0)) {
                Real lfus2 = Real(lsub) - lvap[k];
                Real xnc = nc1d[k] + ncten[k] * dt;
                qiten[k] += xrc * odt;
                niten[k] += xnc * odt;
                qcten[k] -= xrc * odt;
                ncten[k] -= xnc * odt;
                tten[k] += lfus2 * ocp[k] * xrc * odt;
            }
        }

        // Final state updates
        for (int k = 0; k < nz; ++k) {
            Real t_final = (th(i, j, k) * pii(i, j, k)) + tten[k] * dt;
            th(i, j, k) = t_final / pii(i, j, k);

            Real qv_final = qv(i, j, k) + qvten[k] * dt;
            qv(i, j, k) = (qv_final > Real(1.0e-10)) ? qv_final : Real(1.0e-10);

            Real qc_final = qc1d[k] + qcten[k] * dt;
            if (qc_final <= Real(R1)) {
                qc(i, j, k) = Real(0.0);
            } else {
                qc(i, j, k) = qc_final;
            }

            Real qi_final = qi1d[k] + qiten[k] * dt;
            Real ni_final = ni1d[k] + niten[k] * dt;
            Real ni_min = Real(R2) / rho[k];
            if (ni_final < ni_min) ni_final = ni_min;

            if (qi_final <= Real(R1)) {
                qi(i, j, k) = Real(0.0);
                ni(i, j, k) = Real(0.0);
            } else {
                qi(i, j, k) = qi_final;
                Real lami = std::pow(Real(am_i) * tables.cig[1] * tables.oig1 * ni_final / qi_final, obmi);
                Real ilami = Real(1.0) / lami;
                Real xDi = (Real(bm_i) + Real(mu_i) + Real(1.0)) * ilami;
                if (xDi < Real(5.0e-6)) {
                    lami = tables.cie[1] / Real(5.0e-6);
                } else if (xDi > Real(300.0e-6)) {
                    lami = tables.cie[1] / Real(300.0e-6);
                }
                Real term = tables.cig[0] * tables.oig2 * qi_final / Real(am_i) * std::pow(lami, Real(bm_i));
                Real max_ni = Real(999.0e3) / rho[k];
                ni(i, j, k) = (term < max_ni) ? term : max_ni;
            }

            Real qr_final = qr1d[k] + qrten[k] * dt;
            Real nr_final = nr1d[k] + nrten[k] * dt;
            Real nr_min = Real(R2) / rho[k];
            if (nr_final < nr_min) nr_final = nr_min;

            if (qr_final <= Real(R1)) {
                qr(i, j, k) = Real(0.0);
                nr(i, j, k) = Real(0.0);
            } else {
                qr(i, j, k) = qr_final;
                Real lamr = std::pow(Real(am_r) * tables.crg[2] * tables.org2 * nr_final / qr_final, obmr);
                Real mvd = (Real(3.0) + Real(mu_r) + Real(0.672)) / lamr;
                if (mvd > Real(2.5e-3)) {
                    mvd = Real(2.5e-3);
                } else if (mvd < Real(D0r) * Real(0.75)) {
                    mvd = Real(D0r) * Real(0.75);
                }
                lamr = (Real(3.0) + Real(mu_r) + Real(0.672)) / mvd;
                nr(i, j, k) = tables.crg[1] * tables.org3 * qr_final * std::pow(lamr, Real(bm_r)) / Real(am_r);
            }

            Real qs_final = qs1d[k] + qsten[k] * dt;
            qs(i, j, k) = (qs_final > Real(R1)) ? qs_final : Real(0.0);

            Real qg_final = qg1d[k] + qgten[k] * dt;
            if (qg_final <= Real(R1)) {
                qg(i, j, k) = Real(0.0);
            } else {
                qg(i, j, k) = qg_final;
            }
        }

        // Precipitation totals
        Real rain_step = pptrain + pptsnow + pptgraul + pptice;
        rainncv(i, j) = rain_step;
        rainnc(i, j) += rain_step;
        Real snow_step = pptsnow + pptice;
        snowncv(i, j) = snow_step;
        snownc(i, j) += snow_step;
        graupelncv(i, j) = pptgraul;
        graupelnc(i, j) += pptgraul;
        sr(i, j) = (pptsnow + pptgraul + pptice) / (rain_step + Real(1.0e-12));
    }
};

// Host runner function
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
    const ThompsonTablesView<Real>& tables);

} // namespace cppwrf::phys::mp_thompson
