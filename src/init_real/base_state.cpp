#include "base_state.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace cppwrf::init_real {

void compute_base_state(WrfInputData& data, int hypsometric_opt) {
    const int nx = data.dims.nx;
    const int ny = data.dims.ny;
    const int nz = data.dims.nz;
    const int nz_stag = data.dims.nz_stag;

    const float p00 = data.p00;
    const float t00 = data.t00;
    const float a = data.tlp;
    const float tiso = data.tiso;
    const float p_top = data.p_top;
    const float p_strat = data.p_strat;
    const float a_strat = data.tlp_strat;

    const float r_d = 287.0f;
    const float g = 9.81f;
    const float cp = 7.0f * r_d / 2.0f; // 1004.5
    const float cv = cp - r_d;          // 717.5
    const float cvpm = -cv / cp;
    const float p1000mb = 100000.0f;
    const float t0 = 300.0f;

    data.mub = Array2D(ny, nx);
    data.pb = Array3D(nz, ny, nx);
    data.t_init = Array3D(nz, ny, nx);
    data.alb = Array3D(nz, ny, nx);
    data.phb = Array3D(nz_stag, ny, nx);

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            float ht = data.hgt(j, i);
            float arg = (t00 / a) * (t00 / a) - 2.0f * g * ht / (a * r_d);
            float p_surf = p00 * std::exp(-t00 / a + std::sqrt(arg));
            float mub_val = p_surf - p_top;
            data.mub(j, i) = mub_val;

            for (int k = 0; k < nz; ++k) {
                float pb_val = data.c3h[k] * mub_val + data.c4h[k] + p_top;
                data.pb(k, j, i) = pb_val;

                float temp = std::max(tiso, t00 + a * std::log(pb_val / p00));
                if (p_strat > 0.0f && pb_val < p_strat) {
                    temp = tiso + a_strat * std::log(pb_val / p_strat);
                }
                float t_init_val = temp * std::pow(p00 / pb_val, r_d / cp) - t0;
                data.t_init(k, j, i) = t_init_val;

                float alb_val = (r_d / p1000mb) * (t_init_val + t0) * std::pow(pb_val / p1000mb, cvpm);
                data.alb(k, j, i) = alb_val;
            }

            // Integrate base geopotential on full levels
            data.phb(0, j, i) = ht * g;

            if (hypsometric_opt == 1) {
                for (int k = 1; k < nz_stag; ++k) {
                    int km1 = k - 1;
                    float c1_term = data.c1h[km1] * mub_val + data.c2h[km1];
                    data.phb(k, j, i) = data.phb(km1, j, i) - data.dnw[km1] * c1_term * data.alb(km1, j, i);
                }
            } else if (hypsometric_opt == 2) {
                for (int k = 1; k < nz_stag; ++k) {
                    int km1 = k - 1;
                    float pfu = data.c3f[k] * mub_val + data.c4f[k] + p_top;
                    float pfd = data.c3f[km1] * mub_val + data.c4f[km1] + p_top;
                    float phm = data.c3h[km1] * mub_val + data.c4h[km1] + p_top;
                    data.phb(k, j, i) = data.phb(km1, j, i) + data.alb(km1, j, i) * phm * std::log(pfd / pfu);
                }
            } else {
                throw std::runtime_error("Invalid hypsometric_opt: " + std::to_string(hypsometric_opt));
            }
        }
    }
}

} // namespace cppwrf::init_real
