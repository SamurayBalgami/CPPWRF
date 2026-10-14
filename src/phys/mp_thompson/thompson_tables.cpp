#include "phys/mp_thompson/thompson_tables.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <string>
#include <thread>
#if defined(_OPENMP)
#include <omp.h>
#endif

namespace cppwrf::phys::mp_thompson {

namespace {

template <typename F>
void par_for_host(int n_items, F&& fn) {
#if defined(_OPENMP)
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n_items; ++i) fn(i);
#else
    const int nthreads = 4;
    std::atomic<int> next{0};
    std::vector<std::thread> threads;
    threads.reserve(nthreads);
    for (int t = 0; t < nthreads; ++t) {
        threads.emplace_back([&]() {
            for (;;) {
                int i = next.fetch_add(1, std::memory_order_relaxed);
                if (i >= n_items) break;
                fn(i);
            }
        });
    }
    for (auto& th : threads) th.join();
#endif
}

// Exact Numerical Recipes / Lancaster approximations matching module_mp_thompson.F
inline double gammln(double xx) {
    constexpr double stp = 2.5066282746310005;
    constexpr double cof[6] = {
        76.18009172947146, -86.50532032941677,
        24.01409824083091, -1.231739572450155,
        0.1208650973866179e-2, -0.5395239384953e-5
    };
    double x = xx;
    double y = x;
    double tmp = x + 5.5;
    tmp = (x + 0.5) * std::log(tmp) - tmp;
    double ser = 1.000000000190015;
    for (int j = 0; j < 6; ++j) {
        y += 1.0;
        ser += cof[j] / y;
    }
    return tmp + std::log(stp * ser / x);
}

inline double wgamma(double y) {
    return std::exp(gammln(y));
}

inline void gser(double& gamser, double a, double x, double& gln) {
    constexpr int itmax = 100;
    constexpr double geps = 3.0e-7;
    gln = gammln(a);
    if (x <= 0.0) {
        gamser = 0.0;
        return;
    }
    double ap = a;
    double sum = 1.0 / a;
    double del = sum;
    for (int n = 1; n <= itmax; ++n) {
        ap += 1.0;
        del = del * x / ap;
        sum += del;
        if (std::abs(del) < std::abs(sum) * geps) break;
    }
    gamser = sum * std::exp(-x + a * std::log(x) - gln);
}

inline void gcf(double& gammcf, double a, double x, double& gln) {
    constexpr int itmax = 100;
    constexpr double geps = 3.0e-7;
    constexpr double fpmin = 1.0e-30;
    gln = gammln(a);
    double b = x + 1.0 - a;
    double c = 1.0 / fpmin;
    double d = 1.0 / b;
    double h = d;
    for (int i = 1; i <= itmax; ++i) {
        double an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (std::abs(d) < fpmin) d = fpmin;
        c = b + an / c;
        if (std::abs(c) < fpmin) c = fpmin;
        d = 1.0 / d;
        double del = d * c;
        h = h * del;
        if (std::abs(del - 1.0) < geps) break;
    }
    gammcf = std::exp(-x + a * std::log(x) - gln) * h;
}

inline double gammp(double a, double x) {
    if (x < 0.0 || a <= 0.0) return 0.0;
    double gamser = 0, gammcf = 0, gln = 0;
    if (x < a + 1.0) {
        gser(gamser, a, x, gln);
        return gamser;
    } else {
        gcf(gammcf, a, x, gln);
        return 1.0 - gammcf;
    }
}

} // namespace

template <typename Real>
void init_thompson_constants(ThompsonHostTables<Real>& tables) {
    constexpr double PI = 3.1415926536;
    constexpr double rho_w = 1000.0;
    constexpr double rho_i = 890.0;
    constexpr double rho_g_val = 400.0;

    constexpr double am_r = PI * rho_w / 6.0;
    constexpr double bm_r = 3.0;
    constexpr double am_s = 0.069;
    constexpr double bm_s = 2.0;
    constexpr double am_g = PI * rho_g_val / 6.0;
    constexpr double bm_g = 3.0;
    constexpr double am_i = PI * rho_i / 6.0;
    constexpr double bm_i = 3.0;

    constexpr double av_r = 4854.0;
    constexpr double bv_r = 1.0;
    constexpr double av_s = 40.0;
    constexpr double bv_s = 0.55;
    constexpr double av_g = 442.0;
    constexpr double bv_g = 0.89;
    constexpr double bv_i = 1.0;
    constexpr double bv_c = 2.0;

    constexpr double C_cube = 0.5;
    constexpr double C_sqrd = 0.15;

    constexpr double Sc = 0.632;
    constexpr double lsub = 2.834e6;
    constexpr double lvap0 = 2.5e6;
    constexpr double lfus = lsub - lvap0;
    constexpr double olfus = 1.0 / lfus;

    constexpr double xm0i = 1.e-12;
    constexpr double D0r = 50.e-6;
    constexpr double D0s = 300.e-6;
    constexpr double D0g = 350.e-6;

    constexpr double mu_r = 0.0;
    constexpr double mu_g = 0.0;
    constexpr double mu_i = 0.0;
    constexpr double mu_s = 0.6357;

    constexpr double sa[10] = {5.065339, -0.062659, -3.032362, 0.029469, -0.000285,
                               0.31255, 0.000204, 0.003199, 0.0, -0.015952};
    constexpr double sb[10] = {0.476221, -0.015896, 0.165977, 0.007468, -0.000141,
                               0.060366, 0.000079, 0.000594, 0.0, -0.003577};

    const double r_c_init[ntb_c] = {
        1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6,
        1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5,
        1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4,
        1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3,
        1.e-2
    };
    for (int i = 0; i < ntb_c; ++i) tables.r_c[i] = Real(r_c_init[i]);

    const double r_i_init[ntb_i] = {
        1.e-10,2.e-10,3.e-10,4.e-10,
        5.e-10,6.e-10,7.e-10,8.e-10,9.e-10,
        1.e-9,2.e-9,3.e-9,4.e-9,5.e-9,6.e-9,7.e-9,8.e-9,9.e-9,
        1.e-8,2.e-8,3.e-8,4.e-8,5.e-8,6.e-8,7.e-8,8.e-8,9.e-8,
        1.e-7,2.e-7,3.e-7,4.e-7,5.e-7,6.e-7,7.e-7,8.e-7,9.e-7,
        1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6,
        1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5,
        1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4,
        1.e-3
    };
    for (int i = 0; i < ntb_i; ++i) tables.r_i[i] = Real(r_i_init[i]);

    for (int i = 0; i < ntb_r; ++i) tables.r_r[i] = Real(r_c_init[i]);
    for (int i = 0; i < ntb_g; ++i) tables.r_g[i] = Real(r_c_init[i]);
    for (int i = 0; i < ntb_s; ++i) tables.r_s[i] = Real(r_c_init[i]);

    const double N0r_exp_init[ntb_r1] = {
        1.e6,2.e6,3.e6,4.e6,5.e6,6.e6,7.e6,8.e6,9.e6,
        1.e7,2.e7,3.e7,4.e7,5.e7,6.e7,7.e7,8.e7,9.e7,
        1.e8,2.e8,3.e8,4.e8,5.e8,6.e8,7.e8,8.e8,9.e8,
        1.e9,2.e9,3.e9,4.e9,5.e9,6.e9,7.e9,8.e9,9.e9,
        1.e10
    };
    for (int i = 0; i < ntb_r1; ++i) tables.N0r_exp[i] = Real(N0r_exp_init[i]);

    const double N0g_exp_init[ntb_g1] = {
        1.e2,2.e2,3.e2,4.e2,5.e2,6.e2,7.e2,8.e2,9.e2,
        1.e3,2.e3,3.e3,4.e3,5.e3,6.e3,7.e3,8.e3,9.e3,
        1.e4,2.e4,3.e4,4.e4,5.e4,6.e4,7.e4,8.e4,9.e4,
        1.e5,2.e5,3.e5,4.e5,5.e5,6.e5,7.e5,8.e5,9.e5,
        1.e6
    };
    for (int i = 0; i < ntb_g1; ++i) tables.N0g_exp[i] = Real(N0g_exp_init[i]);

    const double Nt_i_init[ntb_i1] = {
        1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,
        1.e1,2.e1,3.e1,4.e1,5.e1,6.e1,7.e1,8.e1,9.e1,
        1.e2,2.e2,3.e2,4.e2,5.e2,6.e2,7.e2,8.e2,9.e2,
        1.e3,2.e3,3.e3,4.e3,5.e3,6.e3,7.e3,8.e3,9.e3,
        1.e4,2.e4,3.e4,4.e4,5.e4,6.e4,7.e4,8.e4,9.e4,
        1.e5,2.e5,3.e5,4.e5,5.e5,6.e5,7.e5,8.e5,9.e5,
        1.e6
    };
    for (int i = 0; i < ntb_i1; ++i) tables.Nt_i[i] = Real(Nt_i_init[i]);
    for (int i = 0; i < ntb_IN; ++i) tables.Nt_IN[i] = Real(Nt_i_init[i]);

    const double Tc_init[ntb_t] = {-0.01, -5., -10., -15., -20., -25., -30., -35., -40.};
    for (int i = 0; i < ntb_t; ++i) tables.Tc[i] = Real(Tc_init[i]);

    double Sc3 = std::pow(Sc, 1.0 / 3.0);
    double D0i = std::pow(xm0i / am_i, 1.0 / bm_i);
    double xm0s = am_s * std::pow(D0s, bm_s);
    double xm0g = (PI * 800.0 / 6.0) * std::pow(D0g, bm_g);

    double cce[5][15], ccg[5][15], ocg1[15], ocg2[15];
    for (int n = 0; n < 15; ++n) {
        int n_1 = n + 1;
        cce[0][n] = n_1 + 1.0;
        cce[1][n] = bm_r + n_1 + 1.0;
        cce[2][n] = bm_r + n_1 + 4.0;
        cce[3][n] = n_1 + bv_c + 1.0;
        cce[4][n] = bm_r + n_1 + bv_c + 1.0;
        for (int k = 0; k < 5; ++k) ccg[k][n] = wgamma(cce[k][n]);
        ocg1[n] = 1.0 / ccg[0][n];
        ocg2[n] = 1.0 / ccg[1][n];
    }

    double cie[7], cig[7];
    cie[0] = mu_i + 1.0;
    cie[1] = bm_i + mu_i + 1.0;
    cie[2] = bm_i + mu_i + bv_i + 1.0;
    cie[3] = mu_i + bv_i + 1.0;
    cie[4] = mu_i + 2.0;
    cie[5] = bm_i * 0.5 + mu_i + bv_i + 1.0;
    cie[6] = bm_i * 0.5 + mu_i + 1.0;
    for (int i = 0; i < 7; ++i) cig[i] = wgamma(cie[i]);
    double oig1 = 1.0 / cig[0];
    double oig2 = 1.0 / cig[1];
    double obmi = 1.0 / bm_i;

    double cre[13], crg[13];
    cre[0] = bm_r + 1.0;
    cre[1] = mu_r + 1.0;
    cre[2] = bm_r + mu_r + 1.0;
    cre[3] = bm_r * 2.0 + mu_r + 1.0;
    cre[4] = mu_r + bv_r + 1.0;
    cre[5] = bm_r + mu_r + bv_r + 1.0;
    cre[6] = bm_r * 0.5 + mu_r + bv_r + 1.0;
    cre[7] = bm_r + mu_r + bv_r + 3.0;
    cre[8] = mu_r + bv_r + 3.0;
    cre[9] = mu_r + 2.0;
    cre[10] = 0.5 * (bv_r + 5.0 + 2.0 * mu_r);
    cre[11] = bm_r * 0.5 + mu_r + 1.0;
    cre[12] = bm_r * 2.0 + mu_r + bv_r + 1.0;
    for (int n = 0; n < 13; ++n) crg[n] = wgamma(cre[n]);
    double obmr = 1.0 / bm_r;
    double ore1 = 1.0 / cre[0];
    double org1 = 1.0 / crg[0];
    double org2 = 1.0 / crg[1];
    double org3 = 1.0 / crg[2];

    double cse[17], csg[17];
    cse[0] = bm_s + 1.0;
    cse[1] = bm_s + 2.0;
    cse[2] = bm_s * 2.0;
    cse[3] = bm_s + bv_s + 1.0;
    cse[4] = bm_s * 2.0 + bv_s + 1.0;
    cse[5] = bm_s * 2.0 + 1.0;
    cse[6] = bm_s + mu_s + 1.0;
    cse[7] = bm_s + mu_s + 2.0;
    cse[8] = bm_s + mu_s + 3.0;
    cse[9] = bm_s + mu_s + bv_s + 1.0;
    cse[10] = bm_s * 2.0 + mu_s + bv_s + 1.0;
    cse[11] = bm_s * 2.0 + mu_s + 1.0;
    cse[12] = bv_s + 2.0;
    cse[13] = bm_s + bv_s;
    cse[14] = mu_s + 1.0;
    cse[15] = 1.0 + (1.0 + bv_s) / 2.0;
    cse[16] = bm_s + bv_s + 2.0;
    for (int n = 0; n < 17; ++n) csg[n] = wgamma(cse[n]);
    double oams = 1.0 / am_s;
    double obms = 1.0 / bm_s;
    double ocms = std::pow(oams, obms);

    double cge[12], cgg[12];
    cge[0] = bm_g + 1.0;
    cge[1] = mu_g + 1.0;
    cge[2] = bm_g + mu_g + 1.0;
    cge[3] = bm_g * 2.0 + mu_g + 1.0;
    cge[4] = bm_g * 2.0 + mu_g + bv_g + 1.0;
    cge[5] = bm_g + mu_g + bv_g + 1.0;
    cge[6] = bm_g * 0.5 + mu_g + bv_g + 1.0;
    cge[7] = mu_g + bv_g + 1.0;
    cge[8] = mu_g + bv_g + 3.0;
    cge[9] = mu_g + 2.0;
    cge[10] = 0.5 * (bv_g + 5.0 + 2.0 * mu_g);
    cge[11] = bm_g * 0.5 + mu_g + 1.0;
    for (int n = 0; n < 12; ++n) cgg[n] = wgamma(cge[n]);
    double oamg = 1.0 / am_g;
    double obmg = 1.0 / bm_g;
    double ocmg = std::pow(oamg, obmg);
    double oge1 = 1.0 / cge[0];
    double ogg1 = 1.0 / cgg[0];
    double ogg2 = 1.0 / cgg[1];
    double ogg3 = 1.0 / cgg[2];

    tables.t1_qr_qc = Real(PI * 0.25 * av_r * crg[8]);
    tables.t1_qr_qi = Real(PI * 0.25 * av_r * crg[8]);
    tables.t2_qr_qi = Real(PI * 0.25 * am_r * av_r * crg[7]);
    tables.t1_qs_qc = Real(PI * 0.25 * av_s);
    tables.t1_qs_qi = Real(PI * 0.25 * av_s);
    tables.t1_qr_ev = Real(0.78 * crg[9]);
    tables.t2_qr_ev = Real(0.308 * Sc3 * std::sqrt(av_r) * crg[10]);
    tables.t1_qs_sd = Real(0.86);
    tables.t2_qs_sd = Real(0.28 * Sc3 * std::sqrt(av_s));
    tables.t1_qs_me = Real(PI * 4.0 * C_sqrd * olfus * 0.86);
    tables.t2_qs_me = Real(PI * 4.0 * C_sqrd * olfus * 0.28 * Sc3 * std::sqrt(av_s));
    tables.t1_qg_sd = Real(0.86 * cgg[9]);
    tables.t2_qg_sd = Real(0.28 * Sc3 * std::sqrt(av_g) * cgg[10]);
    tables.t1_qg_me = Real(PI * 4.0 * C_cube * olfus * 0.86 * cgg[9]);
    tables.t2_qg_me = Real(PI * 4.0 * C_cube * olfus * 0.28 * Sc3 * std::sqrt(av_g) * cgg[10]);

    tables.nic2 = static_cast<int>(std::round(std::log10(tables.r_c[0])));
    tables.nii2 = static_cast<int>(std::round(std::log10(tables.r_i[0])));
    tables.nii3 = static_cast<int>(std::round(std::log10(tables.Nt_i[0])));
    tables.nir2 = static_cast<int>(std::round(std::log10(tables.r_r[0])));
    tables.nir3 = static_cast<int>(std::round(std::log10(tables.N0r_exp[0])));
    tables.nis2 = static_cast<int>(std::round(std::log10(tables.r_s[0])));
    tables.nig2 = static_cast<int>(std::round(std::log10(tables.r_g[0])));
    tables.nig3 = static_cast<int>(std::round(std::log10(tables.N0g_exp[0])));
    tables.niIN2 = static_cast<int>(std::round(std::log10(tables.Nt_IN[0])));

    double xDx[nbins + 1];
    xDx[0] = 1.0;
    xDx[nbc] = 3000.0;
    for (int n = 1; n < nbc; ++n) {
        xDx[n] = std::exp((double(n) / double(nbc)) * std::log(xDx[nbc] / xDx[0]) + std::log(xDx[0]));
    }
    double t_Nc[nbc];
    for (int n = 0; n < nbc; ++n) {
        t_Nc[n] = std::sqrt(xDx[n] * xDx[n + 1]) * 1.0e6;
    }
    tables.nic1 = static_cast<int>(std::log(t_Nc[nbc - 1] / t_Nc[0]));

    xDx[0] = D0r;
    xDx[nbr] = 0.005;
    for (int n = 1; n < nbr; ++n) {
        xDx[n] = std::exp((double(n) / double(nbr)) * std::log(xDx[nbr] / xDx[0]) + std::log(xDx[0]));
    }
    double Dr[nbr];
    for (int n = 0; n < nbr; ++n) {
        Dr[n] = std::sqrt(xDx[n] * xDx[n + 1]);
    }

    xDx[0] = D0s;
    xDx[nbs] = 0.02;
    for (int n = 1; n < nbs; ++n) {
        xDx[n] = std::exp((double(n) / double(nbs)) * std::log(xDx[nbs] / xDx[0]) + std::log(xDx[0]));
    }
    double Ds[nbs];
    for (int n = 0; n < nbs; ++n) {
        Ds[n] = std::sqrt(xDx[n] * xDx[n + 1]);
    }

    for (int k = 0; k < 5; ++k) {
        for (int n = 0; n < 15; ++n) {
            tables.cce[k][n] = Real(cce[k][n]);
            tables.ccg[k][n] = Real(ccg[k][n]);
        }
    }
    for (int n = 0; n < 15; ++n) {
        tables.ocg1[n] = Real(ocg1[n]);
        tables.ocg2[n] = Real(ocg2[n]);
    }
    for (int i = 0; i < 7; ++i) {
        tables.cie[i] = Real(cie[i]);
        tables.cig[i] = Real(cig[i]);
    }
    tables.oig1 = Real(oig1); tables.oig2 = Real(oig2); tables.obmi = Real(obmi);

    for (int i = 0; i < 13; ++i) {
        tables.cre[i] = Real(cre[i]);
        tables.crg[i] = Real(crg[i]);
    }
    tables.ore1 = Real(ore1); tables.org1 = Real(org1); tables.org2 = Real(org2); tables.org3 = Real(org3); tables.obmr = Real(obmr);

    for (int i = 0; i < 17; ++i) {
        tables.cse[i] = Real(cse[i]);
        tables.csg[i] = Real(csg[i]);
    }
    tables.oams = Real(oams); tables.obms = Real(obms); tables.ocms = Real(ocms);

    for (int i = 0; i < 12; ++i) {
        tables.cge[i] = Real(cge[i]);
        tables.cgg[i] = Real(cgg[i]);
    }
    tables.oamg = Real(oamg); tables.obmg = Real(obmg); tables.ocmg = Real(ocmg);
    tables.oge1 = Real(oge1); tables.ogg1 = Real(ogg1); tables.ogg2 = Real(ogg2); tables.ogg3 = Real(ogg3);

    for (int i = 0; i < 10; ++i) {
        tables.sa[i] = Real(sa[i]);
        tables.sb[i] = Real(sb[i]);
    }

    tables.Dr_0 = Real(Dr[0]);
    tables.log_Dr = Real(std::log(Dr[nbr - 1] / Dr[0]));
    tables.inv_log_Dr = Real(1.0 / std::log(Dr[nbr - 1] / Dr[0]));
    tables.Ds_0 = Real(Ds[0]);
    tables.log_Ds = Real(std::log(Ds[nbs - 1] / Ds[0]));
    tables.inv_log_Ds = Real(1.0 / std::log(Ds[nbs - 1] / Ds[0]));
    tables.t_Nc_0 = Real(t_Nc[0]);
    tables.inv_nic1 = Real(1.0 / tables.nic1);

    tables.t1_qg_qc = Real(PI * 0.25 * av_g * cgg[8]);
    tables.D0i = Real(D0i);
    tables.xm0s = Real(xm0s);
    tables.xm0g = Real(xm0g);
}

template <typename Real>
bool load_thompson_tables(ThompsonHostTables<Real>& tables, const std::string& path) {
    std::string search_paths[] = {
        path,
        "tables_fortran.bin",
        "tests/phys/mp_thompson/build/tables_fortran.bin",
        "build/tables_fortran.bin",
        "../build/tables_fortran.bin"
    };
    std::ifstream in;
    for (const auto& sp : search_paths) {
        if (!sp.empty()) {
            in.open(sp, std::ios::binary);
            if (in.is_open()) break;
        }
    }
    if (!in.is_open()) {
        return false;
    }
    init_thompson_constants(tables);
    auto read_table = [&](std::vector<Real>& vec) {
        if constexpr (std::is_same_v<Real, double>) {
            in.read(reinterpret_cast<char*>(vec.data()), vec.size() * sizeof(double));
        } else {
            std::vector<double> buf(vec.size());
            in.read(reinterpret_cast<char*>(buf.data()), buf.size() * sizeof(double));
            for (size_t i = 0; i < vec.size(); ++i) vec[i] = static_cast<float>(buf[i]);
        }
    };
    read_table(tables.tcg_racg);
    read_table(tables.tmr_racg);
    read_table(tables.tcr_gacr);
    read_table(tables.tnr_racg);
    read_table(tables.tnr_gacr);

    read_table(tables.tcs_racs1);
    read_table(tables.tmr_racs1);
    read_table(tables.tcs_racs2);
    read_table(tables.tmr_racs2);
    read_table(tables.tcr_sacr1);
    read_table(tables.tms_sacr1);
    read_table(tables.tcr_sacr2);
    read_table(tables.tms_sacr2);
    read_table(tables.tnr_racs1);
    read_table(tables.tnr_racs2);
    read_table(tables.tnr_sacr1);
    read_table(tables.tnr_sacr2);

    read_table(tables.tpi_qrfz);
    read_table(tables.tpg_qrfz);
    read_table(tables.tni_qrfz);
    read_table(tables.tnr_qrfz);

    read_table(tables.tpi_qcfz);
    read_table(tables.tni_qcfz);

    read_table(tables.tps_iaus);
    read_table(tables.tni_iaus);
    read_table(tables.tpi_ide);

    read_table(tables.t_Efrw);
    read_table(tables.t_Efsw);

    read_table(tables.tpc_wev);
    read_table(tables.tnc_wev);

    return !in.fail();
}

template <typename Real>
void init_thompson_tables(ThompsonHostTables<Real>& tables) {
    // Physical constants matching module_mp_thompson.F
    constexpr double PI = 3.1415926536;
    constexpr double rho_w = 1000.0;
    constexpr double rho_i = 890.0;
    constexpr double rho_g_val = 400.0; // idx_bg1 = 5, rho_g(5) = 400.0

    constexpr double am_r = PI * rho_w / 6.0;
    constexpr double bm_r = 3.0;
    constexpr double am_s = 0.069;
    constexpr double bm_s = 2.0;
    constexpr double am_g = PI * rho_g_val / 6.0;
    constexpr double bm_g = 3.0;
    constexpr double am_i = PI * rho_i / 6.0;
    constexpr double bm_i = 3.0;

    constexpr double av_r = 4854.0;
    constexpr double bv_r = 1.0;
    constexpr double av_s = 40.0;
    constexpr double bv_s = 0.55;
    constexpr double fv_s = 100.0;
    constexpr double av_g = 442.0; // av_g_old
    constexpr double bv_g = 0.89;  // bv_g_old
    constexpr double bv_i = 1.0;
    constexpr double bv_c = 2.0;

    constexpr double C_cube = 0.5;
    constexpr double C_sqrd = 0.15;
    constexpr double Ef_rs = 0.95;
    constexpr double Ef_rg = 0.75;

    constexpr double Sc = 0.632;
    constexpr double lsub = 2.834e6;
    constexpr double lvap0 = 2.5e6;
    constexpr double lfus = lsub - lvap0;
    constexpr double olfus = 1.0 / lfus;

    constexpr double xm0i = 1.e-12;
    constexpr double D0c = 1.e-6;
    constexpr double D0r = 50.e-6;
    constexpr double D0s = 300.e-6;
    constexpr double D0g = 350.e-6;

    constexpr double mu_r = 0.0;
    constexpr double mu_g = 0.0;
    constexpr double mu_i = 0.0;
    constexpr double mu_s = 0.6357;
    constexpr double Kap0 = 490.6;
    constexpr double Kap1 = 17.46;
    constexpr double Lam0 = 20.78;
    constexpr double Lam1 = 3.29;

    constexpr double sa[10] = {5.065339, -0.062659, -3.032362, 0.029469, -0.000285,
                               0.31255, 0.000204, 0.003199, 0.0, -0.015952};
    constexpr double sb[10] = {0.476221, -0.015896, 0.165977, 0.007468, -0.000141,
                               0.060366, 0.000079, 0.000594, 0.0, -0.003577};

    // Table coordinates initialization
    const double r_c_init[ntb_c] = {
        1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6,
        1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5,
        1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4,
        1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3,
        1.e-2
    };
    for (int i = 0; i < ntb_c; ++i) tables.r_c[i] = Real(r_c_init[i]);

    const double r_i_init[ntb_i] = {
        1.e-10,2.e-10,3.e-10,4.e-10,
        5.e-10,6.e-10,7.e-10,8.e-10,9.e-10,
        1.e-9,2.e-9,3.e-9,4.e-9,5.e-9,6.e-9,7.e-9,8.e-9,9.e-9,
        1.e-8,2.e-8,3.e-8,4.e-8,5.e-8,6.e-8,7.e-8,8.e-8,9.e-8,
        1.e-7,2.e-7,3.e-7,4.e-7,5.e-7,6.e-7,7.e-7,8.e-7,9.e-7,
        1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6,
        1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5,
        1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4,
        1.e-3
    };
    for (int i = 0; i < ntb_i; ++i) tables.r_i[i] = Real(r_i_init[i]);

    for (int i = 0; i < ntb_r; ++i) tables.r_r[i] = Real(r_c_init[i]);
    for (int i = 0; i < ntb_g; ++i) tables.r_g[i] = Real(r_c_init[i]);
    for (int i = 0; i < ntb_s; ++i) tables.r_s[i] = Real(r_c_init[i]);

    const double N0r_exp_init[ntb_r1] = {
        1.e6,2.e6,3.e6,4.e6,5.e6,6.e6,7.e6,8.e6,9.e6,
        1.e7,2.e7,3.e7,4.e7,5.e7,6.e7,7.e7,8.e7,9.e7,
        1.e8,2.e8,3.e8,4.e8,5.e8,6.e8,7.e8,8.e8,9.e8,
        1.e9,2.e9,3.e9,4.e9,5.e9,6.e9,7.e9,8.e9,9.e9,
        1.e10
    };
    for (int i = 0; i < ntb_r1; ++i) tables.N0r_exp[i] = Real(N0r_exp_init[i]);

    const double N0g_exp_init[ntb_g1] = {
        1.e2,2.e2,3.e2,4.e2,5.e2,6.e2,7.e2,8.e2,9.e2,
        1.e3,2.e3,3.e3,4.e3,5.e3,6.e3,7.e3,8.e3,9.e3,
        1.e4,2.e4,3.e4,4.e4,5.e4,6.e4,7.e4,8.e4,9.e4,
        1.e5,2.e5,3.e5,4.e5,5.e5,6.e5,7.e5,8.e5,9.e5,
        1.e6
    };
    for (int i = 0; i < ntb_g1; ++i) tables.N0g_exp[i] = Real(N0g_exp_init[i]);

    const double Nt_i_init[ntb_i1] = {
        1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,
        1.e1,2.e1,3.e1,4.e1,5.e1,6.e1,7.e1,8.e1,9.e1,
        1.e2,2.e2,3.e2,4.e2,5.e2,6.e2,7.e2,8.e2,9.e2,
        1.e3,2.e3,3.e3,4.e3,5.e3,6.e3,7.e3,8.e3,9.e3,
        1.e4,2.e4,3.e4,4.e4,5.e4,6.e4,7.e4,8.e4,9.e4,
        1.e5,2.e5,3.e5,4.e5,5.e5,6.e5,7.e5,8.e5,9.e5,
        1.e6
    };
    for (int i = 0; i < ntb_i1; ++i) tables.Nt_i[i] = Real(Nt_i_init[i]);
    for (int i = 0; i < ntb_IN; ++i) tables.Nt_IN[i] = Real(Nt_i_init[i]);

    const double Tc_init[ntb_t] = {-0.01, -5., -10., -15., -20., -25., -30., -35., -40.};
    for (int i = 0; i < ntb_t; ++i) tables.Tc[i] = Real(Tc_init[i]);

    // Exponents and Gamma calculations
    double Sc3 = std::pow(Sc, 1.0 / 3.0);
    double D0i = std::pow(xm0i / am_i, 1.0 / bm_i);
    double xm0s = am_s * std::pow(D0s, bm_s);
    double xm0g = (PI * 800.0 / 6.0) * std::pow(D0g, bm_g); // am_g(NRHG) with rho_g=800

    double cce[5][15], ccg[5][15], ocg1[15], ocg2[15];
    for (int n = 0; n < 15; ++n) {
        int n_1 = n + 1;
        cce[0][n] = n_1 + 1.0;
        cce[1][n] = bm_r + n_1 + 1.0;
        cce[2][n] = bm_r + n_1 + 4.0;
        cce[3][n] = n_1 + bv_c + 1.0;
        cce[4][n] = bm_r + n_1 + bv_c + 1.0;
        for (int k = 0; k < 5; ++k) ccg[k][n] = wgamma(cce[k][n]);
        ocg1[n] = 1.0 / ccg[0][n];
        ocg2[n] = 1.0 / ccg[1][n];
    }

    double cie[7], cig[7];
    cie[0] = mu_i + 1.0;
    cie[1] = bm_i + mu_i + 1.0;
    cie[2] = bm_i + mu_i + bv_i + 1.0;
    cie[3] = mu_i + bv_i + 1.0;
    cie[4] = mu_i + 2.0;
    cie[5] = bm_i * 0.5 + mu_i + bv_i + 1.0;
    cie[6] = bm_i * 0.5 + mu_i + 1.0;
    for (int i = 0; i < 7; ++i) cig[i] = wgamma(cie[i]);
    double oig1 = 1.0 / cig[0];
    double oig2 = 1.0 / cig[1];
    double obmi = 1.0 / bm_i;

    double cre[13], crg[13];
    cre[0] = bm_r + 1.0;
    cre[1] = mu_r + 1.0;
    cre[2] = bm_r + mu_r + 1.0;
    cre[3] = bm_r * 2.0 + mu_r + 1.0;
    cre[4] = mu_r + bv_r + 1.0;
    cre[5] = bm_r + mu_r + bv_r + 1.0;
    cre[6] = bm_r * 0.5 + mu_r + bv_r + 1.0;
    cre[7] = bm_r + mu_r + bv_r + 3.0;
    cre[8] = mu_r + bv_r + 3.0;
    cre[9] = mu_r + 2.0;
    cre[10] = 0.5 * (bv_r + 5.0 + 2.0 * mu_r);
    cre[11] = bm_r * 0.5 + mu_r + 1.0;
    cre[12] = bm_r * 2.0 + mu_r + bv_r + 1.0;
    for (int n = 0; n < 13; ++n) crg[n] = wgamma(cre[n]);
    double obmr = 1.0 / bm_r;
    double ore1 = 1.0 / cre[0];
    double org1 = 1.0 / crg[0];
    double org2 = 1.0 / crg[1];
    double org3 = 1.0 / crg[2];

    double cse[17], csg[17];
    cse[0] = bm_s + 1.0;
    cse[1] = bm_s + 2.0;
    cse[2] = bm_s * 2.0;
    cse[3] = bm_s + bv_s + 1.0;
    cse[4] = bm_s * 2.0 + bv_s + 1.0;
    cse[5] = bm_s * 2.0 + 1.0;
    cse[6] = bm_s + mu_s + 1.0;
    cse[7] = bm_s + mu_s + 2.0;
    cse[8] = bm_s + mu_s + 3.0;
    cse[9] = bm_s + mu_s + bv_s + 1.0;
    cse[10] = bm_s * 2.0 + mu_s + bv_s + 1.0;
    cse[11] = bm_s * 2.0 + mu_s + 1.0;
    cse[12] = bv_s + 2.0;
    cse[13] = bm_s + bv_s;
    cse[14] = mu_s + 1.0;
    cse[15] = 1.0 + (1.0 + bv_s) / 2.0;
    cse[16] = bm_s + bv_s + 2.0;
    for (int n = 0; n < 17; ++n) csg[n] = wgamma(cse[n]);
    double oams = 1.0 / am_s;
    double obms = 1.0 / bm_s;
    double ocms = std::pow(oams, obms);

    double cge[12], cgg[12];
    cge[0] = bm_g + 1.0;
    cge[1] = mu_g + 1.0;
    cge[2] = bm_g + mu_g + 1.0;
    cge[3] = bm_g * 2.0 + mu_g + 1.0;
    cge[4] = bm_g * 2.0 + mu_g + bv_g + 1.0;
    cge[5] = bm_g + mu_g + bv_g + 1.0;
    cge[6] = bm_g * 0.5 + mu_g + bv_g + 1.0;
    cge[7] = mu_g + bv_g + 1.0;
    cge[8] = mu_g + bv_g + 3.0;
    cge[9] = mu_g + 2.0;
    cge[10] = 0.5 * (bv_g + 5.0 + 2.0 * mu_g);
    cge[11] = bm_g * 0.5 + mu_g + 1.0;
    for (int n = 0; n < 12; ++n) cgg[n] = wgamma(cge[n]);
    double oamg = 1.0 / am_g;
    double obmg = 1.0 / bm_g;
    double ocmg = std::pow(oamg, obmg);
    double oge1 = 1.0 / cge[0];
    double ogg1 = 1.0 / cgg[0];
    double ogg2 = 1.0 / cgg[1];
    double ogg3 = 1.0 / cgg[2];

    // Precomputed constants
    tables.t1_qr_qc = Real(PI * 0.25 * av_r * crg[8]);
    tables.t1_qr_qi = Real(PI * 0.25 * av_r * crg[8]);
    tables.t2_qr_qi = Real(PI * 0.25 * am_r * av_r * crg[7]);
    tables.t1_qs_qc = Real(PI * 0.25 * av_s);
    tables.t1_qs_qi = Real(PI * 0.25 * av_s);
    tables.t1_qr_ev = Real(0.78 * crg[9]);
    tables.t2_qr_ev = Real(0.308 * Sc3 * std::sqrt(av_r) * crg[10]);
    tables.t1_qs_sd = Real(0.86);
    tables.t2_qs_sd = Real(0.28 * Sc3 * std::sqrt(av_s));
    tables.t1_qs_me = Real(PI * 4.0 * C_sqrd * olfus * 0.86);
    tables.t2_qs_me = Real(PI * 4.0 * C_sqrd * olfus * 0.28 * Sc3 * std::sqrt(av_s));
    tables.t1_qg_sd = Real(0.86 * cgg[9]);
    tables.t2_qg_sd = Real(0.28 * Sc3 * std::sqrt(av_g) * cgg[10]);
    tables.t1_qg_me = Real(PI * 4.0 * C_cube * olfus * 0.86 * cgg[9]);
    tables.t2_qg_me = Real(PI * 4.0 * C_cube * olfus * 0.28 * Sc3 * std::sqrt(av_g) * cgg[10]);

    // Index rounding constants
    tables.nic2 = static_cast<int>(std::round(std::log10(tables.r_c[0])));
    tables.nii2 = static_cast<int>(std::round(std::log10(tables.r_i[0])));
    tables.nii3 = static_cast<int>(std::round(std::log10(tables.Nt_i[0])));
    tables.nir2 = static_cast<int>(std::round(std::log10(tables.r_r[0])));
    tables.nir3 = static_cast<int>(std::round(std::log10(tables.N0r_exp[0])));
    tables.nis2 = static_cast<int>(std::round(std::log10(tables.r_s[0])));
    tables.nig2 = static_cast<int>(std::round(std::log10(tables.r_g[0])));
    tables.nig3 = static_cast<int>(std::round(std::log10(tables.N0g_exp[0])));
    tables.niIN2 = static_cast<int>(std::round(std::log10(tables.Nt_IN[0])));

    // Bins
    double Dc[nbc], dtc[nbc];
    Dc[0] = D0c;
    dtc[0] = D0c;
    for (int n = 1; n < nbc; ++n) {
        Dc[n] = Dc[n - 1] + 1.0e-6;
        dtc[n] = Dc[n] - Dc[n - 1];
    }

    double xDx[nbins + 1];
    xDx[0] = D0i;
    xDx[nbi] = 2.0 * D0s;
    for (int n = 1; n < nbi; ++n) {
        xDx[n] = std::exp((double(n) / double(nbi)) * std::log(xDx[nbi] / xDx[0]) + std::log(xDx[0]));
    }
    double Di[nbi], dti[nbi];
    for (int n = 0; n < nbi; ++n) {
        Di[n] = std::sqrt(xDx[n] * xDx[n + 1]);
        dti[n] = xDx[n + 1] - xDx[n];
    }

    xDx[0] = D0r;
    xDx[nbr] = 0.005;
    for (int n = 1; n < nbr; ++n) {
        xDx[n] = std::exp((double(n) / double(nbr)) * std::log(xDx[nbr] / xDx[0]) + std::log(xDx[0]));
    }
    double Dr[nbr], dtr[nbr];
    for (int n = 0; n < nbr; ++n) {
        Dr[n] = std::sqrt(xDx[n] * xDx[n + 1]);
        dtr[n] = xDx[n + 1] - xDx[n];
    }

    xDx[0] = D0s;
    xDx[nbs] = 0.02;
    for (int n = 1; n < nbs; ++n) {
        xDx[n] = std::exp((double(n) / double(nbs)) * std::log(xDx[nbs] / xDx[0]) + std::log(xDx[0]));
    }
    double Ds[nbs], dts[nbs];
    for (int n = 0; n < nbs; ++n) {
        Ds[n] = std::sqrt(xDx[n] * xDx[n + 1]);
        dts[n] = xDx[n + 1] - xDx[n];
    }

    xDx[0] = D0g;
    xDx[nbg] = 0.05;
    for (int n = 1; n < nbg; ++n) {
        xDx[n] = std::exp((double(n) / double(nbg)) * std::log(xDx[nbg] / xDx[0]) + std::log(xDx[0]));
    }
    double Dg[nbg], dtg[nbg];
    for (int n = 0; n < nbg; ++n) {
        Dg[n] = std::sqrt(xDx[n] * xDx[n + 1]);
        dtg[n] = xDx[n + 1] - xDx[n];
    }

    xDx[0] = 1.0;
    xDx[nbc] = 3000.0;
    for (int n = 1; n < nbc; ++n) {
        xDx[n] = std::exp((double(n) / double(nbc)) * std::log(xDx[nbc] / xDx[0]) + std::log(xDx[0]));
    }
    double t_Nc[nbc];
    for (int n = 0; n < nbc; ++n) {
        t_Nc[n] = std::sqrt(xDx[n] * xDx[n + 1]) * 1.0e6;
    }
    tables.nic1 = static_cast<int>(std::log(t_Nc[nbc - 1] / t_Nc[0]));

    for (int k = 0; k < 5; ++k) {
        for (int n = 0; n < 15; ++n) {
            tables.cce[k][n] = Real(cce[k][n]);
            tables.ccg[k][n] = Real(ccg[k][n]);
        }
    }
    for (int n = 0; n < 15; ++n) {
        tables.ocg1[n] = Real(ocg1[n]);
        tables.ocg2[n] = Real(ocg2[n]);
    }
    for (int i = 0; i < 7; ++i) {
        tables.cie[i] = Real(cie[i]);
        tables.cig[i] = Real(cig[i]);
    }
    tables.oig1 = Real(oig1); tables.oig2 = Real(oig2); tables.obmi = Real(obmi);

    for (int i = 0; i < 13; ++i) {
        tables.cre[i] = Real(cre[i]);
        tables.crg[i] = Real(crg[i]);
    }
    tables.ore1 = Real(ore1); tables.org1 = Real(org1); tables.org2 = Real(org2); tables.org3 = Real(org3); tables.obmr = Real(obmr);

    for (int i = 0; i < 17; ++i) {
        tables.cse[i] = Real(cse[i]);
        tables.csg[i] = Real(csg[i]);
    }
    tables.oams = Real(oams); tables.obms = Real(obms); tables.ocms = Real(ocms);

    for (int i = 0; i < 12; ++i) {
        tables.cge[i] = Real(cge[i]);
        tables.cgg[i] = Real(cgg[i]);
    }
    tables.oamg = Real(oamg); tables.obmg = Real(obmg); tables.ocmg = Real(ocmg);
    tables.oge1 = Real(oge1); tables.ogg1 = Real(ogg1); tables.ogg2 = Real(ogg2); tables.ogg3 = Real(ogg3);

    for (int i = 0; i < 10; ++i) {
        tables.sa[i] = Real(sa[i]);
        tables.sb[i] = Real(sb[i]);
    }

    tables.Dr_0 = Real(Dr[0]);
    tables.log_Dr = Real(std::log(Dr[nbr - 1] / Dr[0]));
    tables.inv_log_Dr = Real(1.0 / std::log(Dr[nbr - 1] / Dr[0]));
    tables.Ds_0 = Real(Ds[0]);
    tables.log_Ds = Real(std::log(Ds[nbs - 1] / Ds[0]));
    tables.inv_log_Ds = Real(1.0 / std::log(Ds[nbs - 1] / Ds[0]));
    tables.t_Nc_0 = Real(t_Nc[0]);
    tables.inv_nic1 = Real(1.0 / tables.nic1);

    tables.t1_qg_qc = Real(PI * 0.25 * av_g * cgg[8]);
    tables.D0i = Real(D0i);
    tables.xm0s = Real(xm0s);
    tables.xm0g = Real(xm0g);

    auto view = tables.view();

    // 1. table_Efrw
    for (int j = 0; j < nbc; ++j) {
        for (int i = 0; i < nbr; ++i) {
            double Ef_rw = 0.0;
            double p = Dc[j] / Dr[i];
            if (Dr[i] < 50.e-6 || Dc[j] < 3.e-6) {
                Ef_rw = 0.0;
            } else if (p > 0.25) {
                double X = Dc[j] * 1.0e6;
                if (Dr[i] < 75.e-6) {
                    Ef_rw = 0.026794 * X - 0.20604;
                } else if (Dr[i] < 125.e-6) {
                    Ef_rw = -0.00066842 * X * X + 0.061542 * X - 0.37089;
                } else if (Dr[i] < 175.e-6) {
                    Ef_rw = 4.091e-06 * X * X * X * X - 0.00030908 * X * X * X
                          + 0.0066237 * X * X - 0.0013687 * X - 0.073022;
                } else if (Dr[i] < 250.e-6) {
                    Ef_rw = 9.6719e-5 * X * X * X - 0.0068901 * X * X + 0.17305 * X - 0.65988;
                } else if (Dr[i] < 350.e-6) {
                    Ef_rw = 9.0488e-5 * X * X * X - 0.006585 * X * X + 0.16606 * X - 0.56125;
                } else {
                    Ef_rw = 0.00010721 * X * X * X - 0.0072962 * X * X + 0.1704 * X - 0.46929;
                }
            } else {
                double vtr = -0.1021 + 4.932e3 * Dr[i] - 0.9551e6 * Dr[i] * Dr[i]
                           + 0.07934e9 * Dr[i] * Dr[i] * Dr[i]
                           - 0.002362e12 * Dr[i] * Dr[i] * Dr[i] * Dr[i];
                double stokes = Dc[j] * Dc[j] * vtr * rho_w / (9.0 * 1.718e-5 * Dr[i]);
                double reynolds = 9.0 * stokes / (p * p * rho_w);
                double F = std::log(reynolds);
                double G = -0.1007 - 0.358 * F + 0.0261 * F * F;
                double K0 = std::exp(G);
                double z = std::log(stokes / (K0 + 1.0e-15));
                double H = 0.1465 + 1.302 * z - 0.607 * z * z + 0.293 * z * z * z;
                double yc0 = 2.0 / PI * std::atan(H);
                Ef_rw = (yc0 + p) * (yc0 + p) / ((1.0 + p) * (1.0 + p));
            }
            tables.t_Efrw[view.idx_efrw(i, j)] = Real(std::max(0.0, std::min(Ef_rw, 0.95)));
        }
    }

    // 2. table_Efsw
    for (int j = 0; j < nbc; ++j) {
        double vtc = 1.19e4 * (1.0e4 * Dc[j] * Dc[j] * 0.25);
        for (int i = 0; i < nbs; ++i) {
            double vts = av_s * std::pow(Ds[i], bv_s) * std::exp(-fv_s * Ds[i]) - vtc;
            double Ds_m = std::pow(am_s * std::pow(Ds[i], bm_s) / am_r, obmr);
            double p = Dc[j] / Ds_m;
            double Ef_sw = 0.0;
            if (p > 0.25 || Ds[i] < D0s || Dc[j] < 6.e-6 || vts < 1.e-3) {
                Ef_sw = 0.0;
            } else {
                double stokes = Dc[j] * Dc[j] * vts * rho_w / (9.0 * 1.718e-5 * Ds_m);
                double reynolds = 9.0 * stokes / (p * p * rho_w);
                double F = std::log(reynolds);
                double G = -0.1007 - 0.358 * F + 0.0261 * F * F;
                double K0 = std::exp(G);
                double z = std::log(stokes / (K0 + 1.0e-15));
                double H = 0.1465 + 1.302 * z - 0.607 * z * z + 0.293 * z * z * z;
                double yc0 = 2.0 / PI * std::atan(H);
                Ef_sw = (yc0 + p) * (yc0 + p) / ((1.0 + p) * (1.0 + p));
            }
            tables.t_Efsw[view.idx_efsw(i, j)] = Real(std::max(0.0, std::min(Ef_sw, 0.95)));
        }
    }

    // 3. table_dropEvap
    double massc[nbc];
    for (int n = 0; n < nbc; ++n) massc[n] = am_r * std::pow(Dc[n], bm_r);

    for (int k = 0; k < nbc; ++k) {
        int nu_c = std::min(15, static_cast<int>(std::round(1000.e6 / t_Nc[k])) + 2);
        int nu_c_idx = nu_c - 1; // 0-based for cce/ccg
        for (int j = 0; j < ntb_c; ++j) {
            double lamc = std::pow(t_Nc[k] * am_r * ccg[1][nu_c_idx] * ocg1[nu_c_idx] / tables.r_c[j], obmr);
            double N0_c = t_Nc[k] * ocg1[nu_c_idx] * std::pow(lamc, cce[0][nu_c_idx]);
            for (int i = 0; i < nbc; ++i) {
                double summ = 0.0;
                double summ2 = 0.0;
                for (int n = 0; n <= i; ++n) {
                    double N_c_n = N0_c * std::pow(Dc[n], nu_c) * std::exp(-lamc * Dc[n]) * dtc[n];
                    summ += massc[n] * N_c_n;
                    summ2 += N_c_n;
                }
                tables.tpc_wev[view.idx_wev(i, j, k)] = Real(summ);
                tables.tnc_wev[view.idx_wev(i, j, k)] = Real(summ2);
            }
        }
    }

    // 4. qr_acr_qg
    double vr[nbr];
    for (int n2 = 0; n2 < nbr; ++n2) {
        vr[n2] = -0.1021 + 4.932e3 * Dr[n2] - 0.9551e6 * Dr[n2] * Dr[n2]
               + 0.07934e9 * Dr[n2] * Dr[n2] * Dr[n2]
               - 0.002362e12 * Dr[n2] * Dr[n2] * Dr[n2] * Dr[n2];
    }
    double vg[nbg];
    for (int n = 0; n < nbg; ++n) {
        vg[n] = av_g * std::pow(Dg[n], bv_g);
    }

    std::vector<double> N_g_pre(ntb_g * ntb_g1 * nbg);
    auto idx_ng = [](int j, int i, int n) {
        return static_cast<size_t>(n) + nbg * (static_cast<size_t>(i) + ntb_g1 * static_cast<size_t>(j));
    };
    for (int j = 0; j < ntb_g; ++j) {
        for (int i = 0; i < ntb_g1; ++i) {
            double lam_exp_g = std::pow(tables.N0g_exp[i] * am_g * cgg[0] / tables.r_g[j], oge1);
            double lamg = lam_exp_g * std::pow(cgg[2] * ogg2 * ogg1, obmg);
            double N0_g = tables.N0g_exp[i] / (cgg[1] * lam_exp_g) * std::pow(lamg, cge[1]);
            for (int n = 0; n < nbg; ++n) {
                N_g_pre[idx_ng(j, i, n)] = N0_g * std::pow(Dg[n], mu_g) * std::exp(-lamg * Dg[n]) * dtg[n];
            }
        }
    }

    par_for_host(ntb_r * ntb_r1, [&](int km) {
        int m = km / ntb_r1;
        int k = km % ntb_r1;

        double lam_exp_r = std::pow(tables.N0r_exp[k] * am_r * crg[0] / tables.r_r[m], ore1);
        double lamr = lam_exp_r * std::pow(crg[2] * org2 * org1, obmr);
        double N0_r = tables.N0r_exp[k] / (crg[1] * lam_exp_r) * std::pow(lamr, cre[1]);
        double N_r[nbr];
        for (int n2 = 0; n2 < nbr; ++n2) {
            N_r[n2] = N0_r * std::pow(Dr[n2], mu_r) * std::exp(-lamr * Dr[n2]) * dtr[n2];
        }

        constexpr int n3 = 0; // dimNRHG1 = 1
        for (int j = 0; j < ntb_g; ++j) {
            for (int i = 0; i < ntb_g1; ++i) {
                const double* N_g = &N_g_pre[idx_ng(j, i, 0)];

                double t1 = 0.0, t2 = 0.0, z1 = 0.0, y1 = 0.0, y2 = 0.0;
                for (int n2 = 0; n2 < nbr; ++n2) {
                    double massr = am_r * std::pow(Dr[n2], bm_r);
                    for (int n = 0; n < nbg; ++n) {
                        double massg = am_g * std::pow(Dg[n], bm_g);
                        double dvg = 0.5 * ((vr[n2] - vg[n]) + std::abs(vr[n2] - vg[n]));
                        double dvr = 0.5 * ((vg[n] - vr[n2]) + std::abs(vg[n] - vr[n2]));
                        double geom = PI * 0.25 * Ef_rg * (Dg[n] + Dr[n2]) * (Dg[n] + Dr[n2]);

                        t1 += geom * dvg * massg * N_g[n] * N_r[n2];
                        z1 += geom * dvg * massr * N_g[n] * N_r[n2];
                        y1 += geom * dvg * N_g[n] * N_r[n2];

                        t2 += geom * dvr * massr * N_g[n] * N_r[n2];
                        y2 += geom * dvr * N_g[n] * N_r[n2];
                    }
                }
                size_t idx = view.idx_racg(i, j, n3, k, m);
                tables.tcg_racg[idx] = Real(t1);
                tables.tmr_racg[idx] = Real(std::min(z1, static_cast<double>(tables.r_r[m])));
                tables.tcr_gacr[idx] = Real(t2);
                tables.tnr_racg[idx] = Real(y1);
                tables.tnr_gacr[idx] = Real(y2);
            }
        }
    });

    // 5. qr_acr_qs
    double D1[nbr];
    for (int n2 = 0; n2 < nbr; ++n2) {
        D1[n2] = std::pow(vr[n2] / av_s, 1.0 / bv_s);
    }
    double vs[nbs];
    for (int n = 0; n < nbs; ++n) {
        vs[n] = 1.5 * av_s * std::pow(Ds[n], bv_s) * std::exp(-fv_s * Ds[n]);
    }

    std::vector<double> N_s_pre(ntb_t * ntb_s * nbs);
    auto idx_ns = [](int j, int i, int n) {
        return static_cast<size_t>(n) + nbs * (static_cast<size_t>(i) + ntb_s * static_cast<size_t>(j));
    };
    for (int j = 0; j < ntb_t; ++j) {
        for (int i = 0; i < ntb_s; ++i) {
            double M2 = tables.r_s[i] * oams;
            double second = 0.0;
            if (bm_s > 2.0 - 1.e-3 && bm_s < 2.0 + 1.e-3) {
                double loga_ = sa[0] + sa[1] * tables.Tc[j] + sa[2] * bm_s
                             + sa[3] * tables.Tc[j] * bm_s + sa[4] * tables.Tc[j] * tables.Tc[j]
                             + sa[5] * bm_s * bm_s + sa[6] * tables.Tc[j] * tables.Tc[j] * bm_s
                             + sa[7] * tables.Tc[j] * bm_s * bm_s + sa[8] * tables.Tc[j] * tables.Tc[j] * tables.Tc[j]
                             + sa[9] * bm_s * bm_s * bm_s;
                double a_ = std::pow(10.0, loga_);
                double b_ = sb[0] + sb[1] * tables.Tc[j] + sb[2] * bm_s
                          + sb[3] * tables.Tc[j] * bm_s + sb[4] * tables.Tc[j] * tables.Tc[j]
                          + sb[5] * bm_s * bm_s + sb[6] * tables.Tc[j] * tables.Tc[j] * bm_s
                          + sb[7] * tables.Tc[j] * bm_s * bm_s + sb[8] * tables.Tc[j] * tables.Tc[j] * tables.Tc[j]
                          + sb[9] * bm_s * bm_s * bm_s;
                second = std::pow(M2 / a_, 1.0 / b_);
            } else {
                second = M2;
            }

            double loga_ = sa[0] + sa[1] * tables.Tc[j] + sa[2] * cse[0]
                         + sa[3] * tables.Tc[j] * cse[0] + sa[4] * tables.Tc[j] * tables.Tc[j]
                         + sa[5] * cse[0] * cse[0] + sa[6] * tables.Tc[j] * tables.Tc[j] * cse[0]
                         + sa[7] * tables.Tc[j] * cse[0] * cse[0] + sa[8] * tables.Tc[j] * tables.Tc[j] * tables.Tc[j]
                         + sa[9] * cse[0] * cse[0] * cse[0];
            double a_ = std::pow(10.0, loga_);
            double b_ = sb[0] + sb[1] * tables.Tc[j] + sb[2] * cse[0]
                      + sb[3] * tables.Tc[j] * cse[0] + sb[4] * tables.Tc[j] * tables.Tc[j]
                      + sb[5] * cse[0] * cse[0] + sb[6] * tables.Tc[j] * tables.Tc[j] * cse[0]
                      + sb[7] * tables.Tc[j] * cse[0] * cse[0] + sb[8] * tables.Tc[j] * tables.Tc[j] * tables.Tc[j]
                      + sb[9] * cse[0] * cse[0] * cse[0];
            double M3 = a_ * std::pow(second, b_);
            double oM3 = 1.0 / M3;
            double Mrat = M2 * (M2 * oM3) * (M2 * oM3) * (M2 * oM3);
            double M0 = std::pow(M2 * oM3, mu_s);
            double slam1 = M2 * oM3 * Lam0;
            double slam2 = M2 * oM3 * Lam1;

            for (int n = 0; n < nbs; ++n) {
                N_s_pre[idx_ns(j, i, n)] = Mrat * (Kap0 * std::exp(-slam1 * Ds[n])
                       + Kap1 * M0 * std::pow(Ds[n], mu_s) * std::exp(-slam2 * Ds[n])) * dts[n];
            }
        }
    }

    par_for_host(ntb_r * ntb_r1, [&](int km) {
        int m = km / ntb_r1;
        int k = km % ntb_r1;

        double lam_exp_r = std::pow(tables.N0r_exp[k] * am_r * crg[0] / tables.r_r[m], ore1);
        double lamr = lam_exp_r * std::pow(crg[2] * org2 * org1, obmr);
        double N0_r = tables.N0r_exp[k] / (crg[1] * lam_exp_r) * std::pow(lamr, cre[1]);
        double N_r[nbr];
        for (int n2 = 0; n2 < nbr; ++n2) {
            N_r[n2] = N0_r * std::pow(Dr[n2], mu_r) * std::exp(-lamr * Dr[n2]) * dtr[n2];
        }

        for (int j = 0; j < ntb_t; ++j) {
            for (int i = 0; i < ntb_s; ++i) {
                const double* N_s = &N_s_pre[idx_ns(j, i, 0)];

                double t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0;
                double z1 = 0.0, z2 = 0.0, z3 = 0.0, z4 = 0.0;
                double y1 = 0.0, y2 = 0.0, y3 = 0.0, y4 = 0.0;

                for (int n2 = 0; n2 < nbr; ++n2) {
                    double massr = am_r * std::pow(Dr[n2], bm_r);
                    for (int n = 0; n < nbs; ++n) {
                        double masss = am_s * std::pow(Ds[n], bm_s);
                        double dvs = 0.5 * ((vr[n2] - vs[n]) + std::abs(vr[n2] - vs[n]));
                        double dvr = 0.5 * ((vs[n] - vr[n2]) + std::abs(vs[n] - vr[n2]));
                        double geom = PI * 0.25 * Ef_rs * (Ds[n] + Dr[n2]) * (Ds[n] + Dr[n2]);

                        if (massr > 1.5 * masss) {
                            t1 += geom * dvs * masss * N_s[n] * N_r[n2];
                            z1 += geom * dvs * massr * N_s[n] * N_r[n2];
                            y1 += geom * dvs * N_s[n] * N_r[n2];
                        } else {
                            t3 += geom * dvs * masss * N_s[n] * N_r[n2];
                            z3 += geom * dvs * massr * N_s[n] * N_r[n2];
                            y3 += geom * dvs * N_s[n] * N_r[n2];
                        }

                        if (massr > 1.5 * masss) {
                            t2 += geom * dvr * massr * N_s[n] * N_r[n2];
                            y2 += geom * dvr * N_s[n] * N_r[n2];
                            z2 += geom * dvr * masss * N_s[n] * N_r[n2];
                        } else {
                            t4 += geom * dvr * massr * N_s[n] * N_r[n2];
                            y4 += geom * dvr * N_s[n] * N_r[n2];
                            z4 += geom * dvr * masss * N_s[n] * N_r[n2];
                        }
                    }
                }
                size_t idx = view.idx_racs(i, j, k, m);
                tables.tcs_racs1[idx] = Real(t1);
                tables.tmr_racs1[idx] = Real(std::min(z1, static_cast<double>(tables.r_r[m])));
                tables.tcs_racs2[idx] = Real(t3);
                tables.tmr_racs2[idx] = Real(z3);
                tables.tcr_sacr1[idx] = Real(t2);
                tables.tms_sacr1[idx] = Real(z2);
                tables.tcr_sacr2[idx] = Real(t4);
                tables.tms_sacr2[idx] = Real(z4);
                tables.tnr_racs1[idx] = Real(y1);
                tables.tnr_racs2[idx] = Real(y3);
                tables.tnr_sacr1[idx] = Real(y2);
                tables.tnr_sacr2[idx] = Real(y4);
            }
        }
    });

    // 6. freezeH2O
    double orho_w = 1.0 / rho_w;
    double massr_all[nbr];
    for (int n2 = 0; n2 < nbr; ++n2) massr_all[n2] = am_r * std::pow(Dr[n2], bm_r);
    double massc_all[nbc];
    for (int n = 0; n < nbc; ++n) massc_all[n] = am_r * std::pow(Dc[n], bm_r);

    par_for_host(ntb_IN * ntb_freeze_temp, [&](int km) {
        int m = km / ntb_freeze_temp;
        int k = km % ntb_freeze_temp; // 0..44 -> 1..45 in Fortran (temp = -1 to -45 C)
        int k_f = k + 1;

        double T_adjust = std::max(-3.0, std::min(3.0 - std::log10(tables.Nt_IN[m]), 3.0));
        double Texp = std::exp(double(k_f) - T_adjust * 1.0) - 1.0;

        for (int j = 0; j < ntb_r1; ++j) {
            for (int i = 0; i < ntb_r; ++i) {
                double lam_exp = std::pow(tables.N0r_exp[j] * am_r * crg[0] / tables.r_r[i], ore1);
                double lamr = lam_exp * std::pow(crg[2] * org2 * org1, obmr);
                double N0_r = tables.N0r_exp[j] / (crg[1] * lam_exp) * std::pow(lamr, cre[1]);
                double sum1 = 0.0, sum2 = 0.0, sumn1 = 0.0, sumn2 = 0.0;
                for (int n2 = nbr - 1; n2 >= 0; --n2) {
                    double N_r = N0_r * std::pow(Dr[n2], mu_r) * std::exp(-lamr * Dr[n2]) * dtr[n2];
                    double vol = massr_all[n2] * orho_w;
                    double prob = std::max(0.0, 1.0 - std::exp(-120.0 * vol * 5.2e-4 * Texp));
                    if (massr_all[n2] < xm0g) {
                        sumn1 += prob * N_r;
                        sum1 += prob * N_r * massr_all[n2];
                    } else {
                        sumn2 += prob * N_r;
                        sum2 += prob * N_r * massr_all[n2];
                    }
                    if ((sum1 + sum2) >= tables.r_r[i]) break;
                }
                size_t idx = view.idx_qrfz(i, j, k, m);
                tables.tpi_qrfz[idx] = Real(sum1);
                tables.tni_qrfz[idx] = Real(sumn1);
                tables.tpg_qrfz[idx] = Real(sum2);
                tables.tnr_qrfz[idx] = Real(sumn2);
            }
        }

        for (int j = 0; j < nbc; ++j) {
            int nu_c = std::min(15, static_cast<int>(std::round(1000.e6 / t_Nc[j])) + 2);
            int nu_c_idx = nu_c - 1;
            for (int i = 0; i < ntb_c; ++i) {
                double lamc = std::pow(t_Nc[j] * am_r * ccg[1][nu_c_idx] * ocg1[nu_c_idx] / tables.r_c[i], obmr);
                double N0_c = t_Nc[j] * ocg1[nu_c_idx] * std::pow(lamc, cce[0][nu_c_idx]);
                double sum1 = 0.0, sumn2 = 0.0;
                for (int n = nbc - 1; n >= 0; --n) {
                    double vol = massc_all[n] * orho_w;
                    double prob = std::max(0.0, 1.0 - std::exp(-120.0 * vol * 5.2e-4 * Texp));
                    double N_c = N0_c * std::pow(Dc[n], nu_c) * std::exp(-lamc * Dc[n]) * dtc[n];
                    sumn2 = std::min(t_Nc[j], sumn2 + prob * N_c);
                    sum1 += prob * N_c * massc_all[n];
                    if (sum1 >= tables.r_c[i]) break;
                }
                size_t idx = view.idx_qcfz(i, j, k, m);
                tables.tpi_qcfz[idx] = Real(sum1);
                tables.tni_qcfz[idx] = Real(sumn2);
            }
        }
    });

    // 7. qi_aut_qs
    for (int j = 0; j < ntb_i1; ++j) {
        for (int i = 0; i < ntb_i; ++i) {
            double lami = std::pow(am_i * cig[1] * oig1 * tables.Nt_i[j] / tables.r_i[i], obmi);
            double Di_mean = (bm_i + mu_i + 1.0) / lami;
            double N0_i = tables.Nt_i[j] * oig1 * std::pow(lami, cie[0]);
            double t1 = 0.0, t2 = 0.0;
            size_t idx = view.idx_iaus(i, j);
            if (Di_mean > 5.0 * D0s) {
                t1 = tables.r_i[i];
                t2 = tables.Nt_i[j];
                tables.tpi_ide[idx] = Real(0);
            } else if (Di_mean < D0i) {
                t1 = 0.0;
                t2 = 0.0;
                tables.tpi_ide[idx] = Real(1);
            } else {
                double xlimit_intg = lami * D0s;
                tables.tpi_ide[idx] = Real(gammp(mu_i + 2.0, xlimit_intg));
                for (int n2 = 0; n2 < nbi; ++n2) {
                    double N_i_n2 = N0_i * std::pow(Di[n2], mu_i) * std::exp(-lami * Di[n2]) * dti[n2];
                    if (Di[n2] >= D0s) {
                        t1 += N_i_n2 * am_i * std::pow(Di[n2], bm_i);
                        t2 += N_i_n2;
                    }
                }
            }
            tables.tps_iaus[idx] = Real(t1);
            tables.tni_iaus[idx] = Real(t2);
        }
    }
}

template void init_thompson_constants<double>(ThompsonHostTables<double>& tables);
template void init_thompson_constants<float>(ThompsonHostTables<float>& tables);
template void init_thompson_tables<double>(ThompsonHostTables<double>& tables);
template void init_thompson_tables<float>(ThompsonHostTables<float>& tables);
template bool load_thompson_tables<double>(ThompsonHostTables<double>& tables, const std::string& path);
template bool load_thompson_tables<float>(ThompsonHostTables<float>& tables, const std::string& path);

} // namespace cppwrf::phys::mp_thompson
