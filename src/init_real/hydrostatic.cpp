#include "hydrostatic.hpp"
#include "vert_interp.hpp"
#include <cmath>
#include <vector>
#include <algorithm>

namespace cppwrf::init_real {

void integrate_hydrostatic(const MetEmData& met,
                           WrfInputData& data,
                           const Array3D& rh_interp,
                           int hypsometric_opt) {
    const int nx = data.dims.nx;
    const int ny = data.dims.ny;
    const int nz = data.dims.nz;
    const int nz_stag = data.dims.nz_stag;

    const float r_d = 287.0f;
    const float r_v = 461.6f;
    const float rvovrd = r_v / r_d;
    const float cp = 7.0f * r_d / 2.0f; // 1004.5
    const float cv = cp - r_d;          // 717.5
    const float cvpm = -cv / cp;
    const float cpovcv = cp / cv;
    const float p1000mb = 100000.0f;
    const float t0 = 300.0f;
    const float p00 = data.p00;
    const float rcp = r_d / cp;

    data.ph = Array3D(nz_stag, ny, nx);
    data.al = Array3D(nz, ny, nx);
    data.w = Array3D(nz_stag, ny, nx, 0.0f); // w is 0 initially

    Array3D alt(nz, ny, nx);

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            float mu_val = data.mu(j, i);
            float mub_val = data.mub(j, i);

            // Step 1: Top-most level pressure perturbation (kk = nz - 1)
            int kk = nz - 1;
            int k = kk + 1; // full level index = nz = 44 (top full level)
            float qtot = data.qvapor(kk, j, i);
            float qvf2 = 1.0f / (1.0f + qtot);
            float qvf1 = qtot * qvf2;

            data.p(kk, j, i) = -0.5f * (data.c1f[k] * mu_val + qvf1 * (data.c1f[k] * mub_val + data.c2f[k])) / data.rdnw[kk] / qvf2;
            float qvf = 1.0f + rvovrd * data.qvapor(kk, j, i);
            alt(kk, j, i) = (r_d / p1000mb) * (data.t(kk, j, i) + t0) * qvf *
                            std::pow((data.p(kk, j, i) + data.pb(kk, j, i)) / p1000mb, cvpm);
            data.al(kk, j, i) = alt(kk, j, i) - data.alb(kk, j, i);

            // Step 2: Downward integration
            for (kk = nz - 2; kk >= 0; --kk) {
                k = kk + 1;
                qtot = 0.5f * (data.qvapor(kk, j, i) + data.qvapor(kk + 1, j, i));
                qvf2 = 1.0f / (1.0f + qtot);
                qvf1 = qtot * qvf2;
                data.p(kk, j, i) = data.p(kk + 1, j, i) -
                                   (data.c1f[k] * mu_val + qvf1 * (data.c1f[k] * mub_val + data.c2f[k])) / qvf2 / data.rdn[k];
                qvf = 1.0f + rvovrd * data.qvapor(kk, j, i);
                alt(kk, j, i) = (r_d / p1000mb) * (data.t(kk, j, i) + t0) * qvf *
                                std::pow((data.p(kk, j, i) + data.pb(kk, j, i)) / p1000mb, cvpm);
                data.al(kk, j, i) = alt(kk, j, i) - data.alb(kk, j, i);
            }

            // Step 3: Hypsometric integration
            std::vector<float> ph_total(nz_stag, 0.0f);
            ph_total[0] = data.phb(0, j, i);

            float mu0 = mu_val + mub_val;
            if (hypsometric_opt == 2) {
                for (int lev = 1; lev < nz_stag; ++lev) {
                    float pfu = data.c3f[lev] * mu0 + data.c4f[lev] + data.p_top;
                    float pfd = data.c3f[lev - 1] * mu0 + data.c4f[lev - 1] + data.p_top;
                    float phm = data.c3h[lev - 1] * mu0 + data.c4h[lev - 1] + data.p_top;
                    ph_total[lev] = ph_total[lev - 1] + alt(lev - 1, j, i) * phm * std::log(pfd / pfu);
                }
            } else {
                for (int lev = 1; lev < nz_stag; ++lev) {
                    int lev_m1 = lev - 1;
                    float c1_term = (data.c1h[lev_m1] * mub_val + data.c2h[lev_m1]) + (data.c1h[lev_m1] * mu_val);
                    ph_total[lev] = ph_total[lev - 1] - data.dnw[lev_m1] * (c1_term * data.al(lev_m1, j, i) + (data.c1h[lev_m1] * mu_val) * data.alb(lev_m1, j, i));
                }
            }

            for (int lev = 0; lev < nz_stag; ++lev) {
                data.ph(lev, j, i) = ph_total[lev] - data.phb(lev, j, i);
            }

            // Step 4: Recompute al and perturbation pressure p from ph
            if (hypsometric_opt == 2) {
                for (int lev = 0; lev < nz; ++lev) {
                    float pfu = data.c3f[lev + 1] * mu0 + data.c4f[lev + 1] + data.p_top;
                    float pfd = data.c3f[lev] * mu0 + data.c4f[lev] + data.p_top;
                    float phm = data.c3h[lev] * mu0 + data.c4h[lev] + data.p_top;
                    float delta_ph = data.ph(lev + 1, j, i) - data.ph(lev, j, i) +
                                     data.phb(lev + 1, j, i) - data.phb(lev, j, i);
                    data.al(lev, j, i) = delta_ph / (phm * std::log(pfd / pfu)) - data.alb(lev, j, i);

                    qvf = 1.0f + rvovrd * data.qvapor(lev, j, i);
                    float density_term = (r_d * (t0 + data.t(lev, j, i)) * qvf) /
                                         (p1000mb * (data.al(lev, j, i) + data.alb(lev, j, i)));
                    data.p(lev, j, i) = p1000mb * std::pow(density_term, cpovcv) - data.pb(lev, j, i);
                }
            }

            // Step 5: Final diagnosis of qvapor from hydrostatic full pressure p_hyd
            for (int lev = 0; lev < nz; ++lev) {
                float p_hyd = data.p(lev, j, i) + data.pb(lev, j, i);
                float t_k = (data.t(lev, j, i) + t0) / std::pow(p00 / p_hyd, rcp);
                data.qvapor(lev, j, i) = rh_to_mxrat_single(rh_interp(lev, j, i), t_k, p_hyd);
            }
        }
    }
}

} // namespace cppwrf::init_real
