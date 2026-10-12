#include "bdy.hpp"
#include <stdexcept>

namespace cppwrf::init_real {

static void compute_coupled_fields(const WrfInputData& state,
                                   Array3D& u_cp,
                                   Array3D& v_cp,
                                   Array3D& t_cp,
                                   Array3D& ph_cp,
                                   Array3D& q_cp) {
    const int nx = state.dims.nx;
    const int ny = state.dims.ny;
    const int nz = state.dims.nz;
    const int nz_stag = state.dims.nz_stag;
    const int nx_stag = state.dims.nx_stag;
    const int ny_stag = state.dims.ny_stag;

    u_cp = Array3D(nz, ny, nx_stag);
    v_cp = Array3D(nz, ny_stag, nx);
    t_cp = Array3D(nz, ny, nx);
    ph_cp = Array3D(nz_stag, ny, nx);
    q_cp = Array3D(nz, ny, nx);

    // Couple U: muu on staggered x-faces, divide by msfuy
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx_stag; ++i) {
            float mu_l = (i == 0) ? state.mu(j, 0) : state.mu(j, i - 1);
            float mu_r = (i == nx) ? state.mu(j, nx - 1) : state.mu(j, i);
            float mub_l = (i == 0) ? state.mub(j, 0) : state.mub(j, i - 1);
            float mub_r = (i == nx) ? state.mub(j, nx - 1) : state.mub(j, i);
            float muu = 0.5f * (mu_l + mu_r + mub_l + mub_r);
            float msf = state.msfuy(j, i);

            for (int k = 0; k < nz; ++k) {
                u_cp(k, j, i) = state.u(k, j, i) * (state.c1h[k] * muu + state.c2h[k]) / msf;
            }
        }
    }

    // Couple V: muv on staggered y-faces, divide by msfvx
    for (int j = 0; j < ny_stag; ++j) {
        for (int i = 0; i < nx; ++i) {
            float mu_d = (j == 0) ? state.mu(0, i) : state.mu(j - 1, i);
            float mu_u = (j == ny) ? state.mu(ny - 1, i) : state.mu(j, i);
            float mub_d = (j == 0) ? state.mub(0, i) : state.mub(j - 1, i);
            float mub_u = (j == ny) ? state.mub(ny - 1, i) : state.mub(j, i);
            float muv = 0.5f * (mu_d + mu_u + mub_d + mub_u);
            float msf = state.msfvx(j, i);

            for (int k = 0; k < nz; ++k) {
                v_cp(k, j, i) = state.v(k, j, i) * (state.c1h[k] * muv + state.c2h[k]) / msf;
            }
        }
    }

    // Couple T, PH, Q on mass points
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            float mu_val = state.mu(j, i);
            float mub_val = state.mub(j, i);

            const float rvovrd = 461.6f / 287.0f;
            for (int k = 0; k < nz; ++k) {
                float mut = (state.c1h[k] * mu_val) + (state.c1h[k] * mub_val + state.c2h[k]);
                float theta = state.t(k, j, i) + 300.0f;
                float thm = theta * (1.0f + rvovrd * state.qvapor(k, j, i)) - 300.0f;
                t_cp(k, j, i) = thm * mut;
                q_cp(k, j, i) = state.qvapor(k, j, i) * mut;
            }

            for (int k = 0; k < nz_stag; ++k) {
                float mut_f = (state.c1f[k] * mu_val) + (state.c1f[k] * mub_val + state.c2f[k]);
                ph_cp(k, j, i) = state.ph(k, j, i) * mut_f;
            }
        }
    }
}

BoundaryInterval pack_single_boundary(const WrfInputData& state,
                                      const std::string& this_time,
                                      const std::string& next_time) {
    BoundaryInterval bdy;
    bdy.this_bdy_time = this_time;
    bdy.next_bdy_time = next_time;

    const int nx = state.dims.nx;
    const int ny = state.dims.ny;
    const int nz = state.dims.nz;
    const int nz_stag = state.dims.nz_stag;
    const int nx_stag = state.dims.nx_stag;
    const int ny_stag = state.dims.ny_stag;
    const int bw = state.dims.bdy_width; // 5

    Array3D u_cp, v_cp, t_cp, ph_cp, q_cp;
    compute_coupled_fields(state, u_cp, v_cp, t_cp, ph_cp, q_cp);

    bdy.u_bxs = Array3D(bw, nz, ny);
    bdy.u_bxe = Array3D(bw, nz, ny);
    bdy.u_bys = Array3D(bw, nz, nx_stag);
    bdy.u_bye = Array3D(bw, nz, nx_stag);

    bdy.v_bxs = Array3D(bw, nz, ny_stag);
    bdy.v_bxe = Array3D(bw, nz, ny_stag);
    bdy.v_bys = Array3D(bw, nz, nx);
    bdy.v_bye = Array3D(bw, nz, nx);

    bdy.t_bxs = Array3D(bw, nz, ny);
    bdy.t_bxe = Array3D(bw, nz, ny);
    bdy.t_bys = Array3D(bw, nz, nx);
    bdy.t_bye = Array3D(bw, nz, nx);

    bdy.ph_bxs = Array3D(bw, nz_stag, ny);
    bdy.ph_bxe = Array3D(bw, nz_stag, ny);
    bdy.ph_bys = Array3D(bw, nz_stag, nx);
    bdy.ph_bye = Array3D(bw, nz_stag, nx);

    bdy.q_bxs = Array3D(bw, nz, ny);
    bdy.q_bxe = Array3D(bw, nz, ny);
    bdy.q_bys = Array3D(bw, nz, nx);
    bdy.q_bye = Array3D(bw, nz, nx);

    bdy.mu_bxs = Array3D(bw, 1, ny);
    bdy.mu_bxe = Array3D(bw, 1, ny);
    bdy.mu_bys = Array3D(bw, 1, nx);
    bdy.mu_bye = Array3D(bw, 1, nx);

    // XS (West boundary)
    for (int w = 0; w < bw; ++w) {
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                bdy.u_bxs(w, k, j) = u_cp(k, j, w);
                bdy.t_bxs(w, k, j) = t_cp(k, j, w);
                bdy.q_bxs(w, k, j) = q_cp(k, j, w);
            }
            for (int j = 0; j < ny_stag; ++j) {
                bdy.v_bxs(w, k, j) = v_cp(k, j, w);
            }
        }
        for (int k = 0; k < nz_stag; ++k) {
            for (int j = 0; j < ny; ++j) {
                bdy.ph_bxs(w, k, j) = ph_cp(k, j, w);
            }
        }
        for (int j = 0; j < ny; ++j) {
            bdy.mu_bxs(w, 0, j) = state.mu(j, w);
        }
    }

    // XE (East boundary)
    for (int w = 0; w < bw; ++w) {
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                bdy.u_bxe(w, k, j) = u_cp(k, j, nx_stag - 1 - w);
                bdy.t_bxe(w, k, j) = t_cp(k, j, nx - 1 - w);
                bdy.q_bxe(w, k, j) = q_cp(k, j, nx - 1 - w);
            }
            for (int j = 0; j < ny_stag; ++j) {
                bdy.v_bxe(w, k, j) = v_cp(k, j, nx - 1 - w);
            }
        }
        for (int k = 0; k < nz_stag; ++k) {
            for (int j = 0; j < ny; ++j) {
                bdy.ph_bxe(w, k, j) = ph_cp(k, j, nx - 1 - w);
            }
        }
        for (int j = 0; j < ny; ++j) {
            bdy.mu_bxe(w, 0, j) = state.mu(j, nx - 1 - w);
        }
    }

    // YS (South boundary)
    for (int w = 0; w < bw; ++w) {
        for (int k = 0; k < nz; ++k) {
            for (int i = 0; i < nx_stag; ++i) {
                bdy.u_bys(w, k, i) = u_cp(k, w, i);
            }
            for (int i = 0; i < nx; ++i) {
                bdy.v_bys(w, k, i) = v_cp(k, w, i);
                bdy.t_bys(w, k, i) = t_cp(k, w, i);
                bdy.q_bys(w, k, i) = q_cp(k, w, i);
            }
        }
        for (int k = 0; k < nz_stag; ++k) {
            for (int i = 0; i < nx; ++i) {
                bdy.ph_bys(w, k, i) = ph_cp(k, w, i);
            }
        }
        for (int i = 0; i < nx; ++i) {
            bdy.mu_bys(w, 0, i) = state.mu(w, i);
        }
    }

    // YE (North boundary)
    for (int w = 0; w < bw; ++w) {
        for (int k = 0; k < nz; ++k) {
            for (int i = 0; i < nx_stag; ++i) {
                bdy.u_bye(w, k, i) = u_cp(k, ny - 1 - w, i);
            }
            for (int i = 0; i < nx; ++i) {
                bdy.v_bye(w, k, i) = v_cp(k, ny_stag - 1 - w, i);
                bdy.t_bye(w, k, i) = t_cp(k, ny - 1 - w, i);
                bdy.q_bye(w, k, i) = q_cp(k, ny - 1 - w, i);
            }
        }
        for (int k = 0; k < nz_stag; ++k) {
            for (int i = 0; i < nx; ++i) {
                bdy.ph_bye(w, k, i) = ph_cp(k, ny - 1 - w, i);
            }
        }
        for (int i = 0; i < nx; ++i) {
            bdy.mu_bye(w, 0, i) = state.mu(ny - 1 - w, i);
        }
    }

    return bdy;
}

static void compute_tendency(const Array3D& b1, const Array3D& b2, float dt, Array3D& bt) {
    bt = Array3D(b1.d1, b1.d2, b1.d3);
    for (size_t i = 0; i < b1.size(); ++i) {
        bt.data[i] = (b2.data[i] - b1.data[i]) / dt;
    }
}

WrfBdyData build_wrfbdy_data(const std::vector<WrfInputData>& states,
                             const std::vector<std::string>& date_strings,
                             float dt_interval) {
    if (states.size() < 2 || states.size() != date_strings.size()) {
        throw std::runtime_error("Need at least 2 states and matching date strings for wrfbdy");
    }

    WrfBdyData bdy_data;
    bdy_data.dims = states[0].dims;
    bdy_data.dt_interval = dt_interval;

    std::vector<BoundaryInterval> raw_intervals;
    for (size_t t = 0; t < states.size(); ++t) {
        std::string next_t = (t + 1 < states.size()) ? date_strings[t + 1] : date_strings[t];
        raw_intervals.push_back(pack_single_boundary(states[t], date_strings[t], next_t));
    }

    // Now compute tendencies for intervals 0 .. states.size()-2
    for (size_t t = 0; t < states.size() - 1; ++t) {
        BoundaryInterval interval = raw_intervals[t];
        const BoundaryInterval& next_interval = raw_intervals[t + 1];

        compute_tendency(interval.u_bxs, next_interval.u_bxs, dt_interval, interval.u_btxs);
        compute_tendency(interval.u_bxe, next_interval.u_bxe, dt_interval, interval.u_btxe);
        compute_tendency(interval.u_bys, next_interval.u_bys, dt_interval, interval.u_btys);
        compute_tendency(interval.u_bye, next_interval.u_bye, dt_interval, interval.u_btye);

        compute_tendency(interval.v_bxs, next_interval.v_bxs, dt_interval, interval.v_btxs);
        compute_tendency(interval.v_bxe, next_interval.v_bxe, dt_interval, interval.v_btxe);
        compute_tendency(interval.v_bys, next_interval.v_bys, dt_interval, interval.v_btys);
        compute_tendency(interval.v_bye, next_interval.v_bye, dt_interval, interval.v_btye);

        compute_tendency(interval.t_bxs, next_interval.t_bxs, dt_interval, interval.t_btxs);
        compute_tendency(interval.t_bxe, next_interval.t_bxe, dt_interval, interval.t_btxe);
        compute_tendency(interval.t_bys, next_interval.t_bys, dt_interval, interval.t_btys);
        compute_tendency(interval.t_bye, next_interval.t_bye, dt_interval, interval.t_btye);

        compute_tendency(interval.ph_bxs, next_interval.ph_bxs, dt_interval, interval.ph_btxs);
        compute_tendency(interval.ph_bxe, next_interval.ph_bxe, dt_interval, interval.ph_btxe);
        compute_tendency(interval.ph_bys, next_interval.ph_bys, dt_interval, interval.ph_btys);
        compute_tendency(interval.ph_bye, next_interval.ph_bye, dt_interval, interval.ph_btye);

        compute_tendency(interval.q_bxs, next_interval.q_bxs, dt_interval, interval.q_btxs);
        compute_tendency(interval.q_bxe, next_interval.q_bxe, dt_interval, interval.q_btxe);
        compute_tendency(interval.q_bys, next_interval.q_bys, dt_interval, interval.q_btys);
        compute_tendency(interval.q_bye, next_interval.q_bye, dt_interval, interval.q_btye);

        compute_tendency(interval.mu_bxs, next_interval.mu_bxs, dt_interval, interval.mu_btxs);
        compute_tendency(interval.mu_bxe, next_interval.mu_bxe, dt_interval, interval.mu_btxe);
        compute_tendency(interval.mu_bys, next_interval.mu_bys, dt_interval, interval.mu_btys);
        compute_tendency(interval.mu_bye, next_interval.mu_bye, dt_interval, interval.mu_btye);

        bdy_data.intervals.push_back(interval);
    }

    return bdy_data;
}

} // namespace cppwrf::init_real
