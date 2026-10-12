#include "vcoord.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

namespace cppwrf::init_real {

void compute_eta_levels(WrfInputData& data,
                        float max_dz,
                        float dzbot,
                        float dzstretch_s,
                        float dzstretch_u) {
    const int nlev = data.dims.nz; // 44
    const float ptop = data.p_top; // 5000.0f
    const float r_d = 287.0f;
    const float g = 9.81f;

    data.znw = Array1D(data.dims.nz_stag); // 45
    data.znu = Array1D(data.dims.nz);      // 44
    data.dnw = Array1D(data.dims.nz);      // 44
    data.dn = Array1D(data.dims.nz, 0.0f);   // 44
    data.rdnw = Array1D(data.dims.nz);     // 44
    data.rdn = Array1D(data.dims.nz, 0.0f);  // 44

    std::vector<float> zup(nlev + 1, 0.0f);
    std::vector<float> pup(nlev + 1, 0.0f);

    float tt = 290.0f;
    float ztop = r_d * tt / g * std::log(1.0e5f / ptop);
    float dz = dzbot;
    zup[1] = dz;
    pup[1] = 1.0e5f * std::exp(-g * zup[1] / r_d / tt);
    data.znw[0] = 1.0f;
    data.znw[1] = (pup[1] - ptop) / (1.0e5f - ptop);

    int isave = 1;
    for (int i = 1; i <= nlev - 1; ++i) {
        float a = dzstretch_u + (dzstretch_s - dzstretch_u) * std::max((max_dz * 0.5f - dz) / (max_dz * 0.5f), 0.0f);
        dz = a * dz;
        float dztest = (ztop - zup[isave]) / static_cast<float>(nlev - isave);
        if (dztest < dz) {
            break;
        }
        isave = i + 1;
        zup[i + 1] = zup[i] + dz;
        pup[i + 1] = 1.0e5f * std::exp(-g * zup[i + 1] / r_d / tt);
        data.znw[i + 1] = (pup[i + 1] - ptop) / (1.0e5f - ptop);
    }

    dz = (ztop - zup[isave]) / static_cast<float>(nlev - isave);
    for (int i = isave; i <= nlev - 1; ++i) {
        zup[i + 1] = zup[i] + dz;
        pup[i + 1] = 1.0e5f * std::exp(-g * zup[i + 1] / r_d / tt);
        data.znw[i + 1] = (pup[i + 1] - ptop) / (1.0e5f - ptop);
    }
    data.znw[nlev] = 0.0f;

    // Derived eta metrics: znu, dnw, rdnw, dn, rdn
    for (int k = 0; k < nlev; ++k) {
        data.dnw[k] = data.znw[k + 1] - data.znw[k];
        data.rdnw[k] = 1.0f / data.dnw[k];
        data.znu[k] = 0.5f * (data.znw[k + 1] + data.znw[k]);
    }

    for (int k = 1; k < nlev; ++k) {
        data.dn[k] = 0.5f * (data.dnw[k] + data.dnw[k - 1]);
        data.rdn[k] = 1.0f / data.dn[k];
    }
}

void compute_vcoord_coeffs(WrfInputData& data) {
    const int kde = data.dims.nz_stag; // 45
    const int nlev = data.dims.nz;     // 44
    const float p1000mb = 100000.0f;
    const float p_top = data.p_top;
    const float etac = data.etac;

    data.c1f = Array1D(kde);
    data.c2f = Array1D(kde);
    data.c3f = Array1D(kde);
    data.c4f = Array1D(kde);

    data.c1h = Array1D(nlev);
    data.c2h = Array1D(nlev);
    data.c3h = Array1D(nlev);
    data.c4h = Array1D(nlev);

    if (data.hybrid_opt == 0) {
        for (int k = 0; k < kde; ++k) {
            data.c3f[k] = data.znw[k];
            data.c1f[k] = 1.0f;
            data.c2f[k] = 0.0f;
            data.c4f[k] = 0.0f;
        }
        for (int k = 0; k < nlev; ++k) {
            data.c3h[k] = data.znu[k];
            data.c1h[k] = 1.0f;
            data.c2h[k] = 0.0f;
            data.c4h[k] = 0.0f;
        }
    } else if (data.hybrid_opt == 2) {
        const float b1 = 2.0f * etac * etac * (1.0f - etac);
        const float b2 = -etac * (4.0f - 3.0f * etac - etac * etac * etac);
        const float b3 = 2.0f * (1.0f - etac * etac * etac);
        const float b4 = -(1.0f - etac * etac);
        const float b5 = (1.0f - etac) * (1.0f - etac) * (1.0f - etac) * (1.0f - etac);

        for (int k = 0; k < kde; ++k) {
            float eta = data.znw[k];
            data.c3f[k] = (b1 + b2 * eta + b3 * eta * eta + b4 * eta * eta * eta) / b5;
            if (eta < etac) {
                data.c3f[k] = 0.0f;
            }
        }
        data.c3f[0] = 1.0f;
        data.c3f[kde - 1] = 0.0f;

        for (int k = 0; k < kde; ++k) {
            data.c4f[k] = (data.znw[k] - data.c3f[k]) * (p1000mb - p_top);
        }

        for (int k = 0; k < nlev; ++k) {
            data.c3h[k] = 0.5f * (data.c3f[k + 1] + data.c3f[k]);
            data.c4h[k] = (data.znu[k] - data.c3h[k]) * (p1000mb - p_top);
        }

        // c1f on full levels
        data.c1f[0] = 1.0f;
        data.c1f[kde - 1] = 0.0f;
        for (int k = 1; k < kde - 1; ++k) {
            data.c1f[k] = (data.c3h[k] - data.c3h[k - 1]) / (data.znu[k] - data.znu[k - 1]);
        }

        for (int k = 0; k < kde; ++k) {
            data.c2f[k] = (1.0f - data.c1f[k]) * (p1000mb - p_top);
        }

        for (int k = 0; k < nlev; ++k) {
            data.c1h[k] = (data.c3f[k + 1] - data.c3f[k]) / (data.znw[k + 1] - data.znw[k]);
            data.c2h[k] = (1.0f - data.c1h[k]) * (p1000mb - p_top);
        }
    } else {
        throw std::runtime_error("Unsupported hybrid_opt: " + std::to_string(data.hybrid_opt));
    }
}

} // namespace cppwrf::init_real
