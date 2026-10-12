#include "vert_interp.hpp"
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>

namespace cppwrf::init_real {

float rh_to_mxrat_single(float rh, float t_k, float p_pa,
                         float qv_max_p_safe,
                         float qv_max_flag,
                         float qv_max_value,
                         float qv_min_p_safe,
                         float qv_min_flag,
                         float qv_min_value) {
    const float svp1 = 0.6112f;
    const float svp2 = 17.67f;
    const float svp3 = 29.65f;
    const float svpt0 = 273.15f;
    const float eps = 0.622f;

    float rh_clamped = std::clamp(rh, 0.0f, 100.0f);
    float q = 1.0e-6f;

    if (t_k != 0.0f) {
        float es = 0.01f * rh_clamped * svp1 * 10.0f * std::exp(svp2 * (t_k - svpt0) / (t_k - svp3));
        float p_mb = p_pa / 100.0f;
        if (es >= p_mb) {
            q = 1.0e-6f;
        } else {
            q = std::max(eps * es / (p_mb - es), 1.0e-6f);
        }
    }

    if (p_pa < qv_max_p_safe && q > qv_max_flag) {
        q = qv_max_value;
    }
    if (p_pa < qv_min_p_safe && q < qv_min_flag) {
        q = qv_min_value;
    }
    return q;
}

static float compute_single_sfcprs3(float ter, float slp, const std::vector<float>& p, const std::vector<float>& h) {
    const int nlev = static_cast<int>(p.size()); // 34
    // If near sea level (ter < 50 m)
    if (ter < 50.0f) {
        return slp + (p[1] - p[2]) / (h[1] - h[2]) * ter;
    }

    // Trapping isobaric levels (1 .. nlev - 3, comparing with k + 1)
    for (int k = 1; k <= nlev - 3; ++k) {
        if (h[k] <= ter && h[k + 1] > ter) {
            float zl = h[k];
            float zu = h[k + 1];
            float zm = ter;
            float pl = p[k];
            float pu = p[k + 1];
            return std::exp((std::log(pl) * (zm - zu) + std::log(pu) * (zl - zm)) / (zl - zu));
        }
    }

    // Extrapolate between slp and first isobaric level
    if (slp >= p[1]) {
        float zl = 0.0f;
        float zu = h[2];
        float zm = ter;
        float pl = slp;
        float pu = p[2];
        return std::exp((std::log(pl) * (zm - zu) + std::log(pu) * (zl - zm)) / (zl - zu));
    } else {
        for (int k = 1; k <= nlev - 4; ++k) {
            if (slp >= p[k + 1] && slp < p[k]) {
                float zl = 0.0f;
                float zu = h[k + 1];
                float zm = ter;
                float pl = slp;
                float pu = p[k + 1];
                return std::exp((std::log(pl) * (zm - zu) + std::log(pu) * (zl - zm)) / (zl - zu));
            }
        }
    }

    return p[1];
}

void compute_dry_pressure_and_mu(const MetEmData& met,
                                 WrfInputData& data,
                                 Array3D& pd_gc,
                                 Array3D& qv_gc) {
    const int nx = data.dims.nx;
    const int ny = data.dims.ny;
    const int nlev = met.dims.num_metgrid_levels; // 34
    const float r_d = 287.0f;
    const float g = 9.81f;
    const float p_top = data.p_top;

    data.mu = Array2D(ny, nx);
    pd_gc = Array3D(nlev, ny, nx);
    qv_gc = Array3D(nlev, ny, nx);

    // Compute qv on metgrid levels
    for (int k = 0; k < nlev; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                qv_gc(k, j, i) = rh_to_mxrat_single(met.rh(k, j, i), met.tt(k, j, i), met.pres(k, j, i));
            }
        }
    }

    // Exact WRF integ_moist + sfcprs3
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            float psfc_raw = met.psfc(j, i);
            float tsfc = met.tt(0, j, i);
            float qsfc = qv_gc(0, j, i);
            float zsfc = met.ght(0, j, i); // GHT(0) is SOILHGT

            std::vector<float> p(nlev), t(nlev), q(nlev), ght(nlev), pd(nlev);
            for (int k = 0; k < nlev; ++k) {
                p[k] = met.pres(k, j, i);
                t[k] = met.tt(k, j, i);
                q[k] = qv_gc(k, j, i);
                ght[k] = met.ght(k, j, i);
            }

            // Find first isobaric level above ground (levels 1 .. nlev-1)
            int level_above_sfc = -1;
            if (p[1] < psfc_raw) {
                level_above_sfc = 1;
            } else {
                for (int k = 1; k < nlev - 1; ++k) {
                    if ((p[k] - psfc_raw >= 0.0f) && (p[k + 1] - psfc_raw < 0.0f)) {
                        level_above_sfc = k + 1;
                        break;
                    }
                }
                if (level_above_sfc == -1) {
                    level_above_sfc = 1;
                }
            }

            float intq = 0.0f;
            pd[nlev - 1] = p[nlev - 1];
            for (int k = nlev - 2; k >= level_above_sfc; --k) {
                float rhobar = 0.5f * (p[k] / (r_d * t[k]) + p[k + 1] / (r_d * t[k + 1]));
                float qbar = 0.5f * (q[k] + q[k + 1]);
                float dz = ght[k + 1] - ght[k];
                intq += g * qbar * rhobar / (1.0f + qbar) * dz;
                pd[k] = p[k] - intq;
            }

            if (level_above_sfc > 0 && (p[level_above_sfc - 1] - psfc_raw >= 0.0f) && (p[level_above_sfc] - psfc_raw < 0.0f)) {
                float p1 = psfc_raw;
                float p2 = p[level_above_sfc];
                float t1 = tsfc;
                float t2 = t[level_above_sfc];
                float q1 = qsfc;
                float q2 = q[level_above_sfc];
                float z1 = zsfc;
                float z2 = ght[level_above_sfc];
                float rhobar = 0.5f * (p1 / (r_d * t1) + p2 / (r_d * t2));
                float qbar = 0.5f * (q1 + q2);
                float dz = z2 - z1;
                if (dz > 0.1f) {
                    intq += g * qbar * rhobar / (1.0f + qbar) * dz;
                }
                for (int k = level_above_sfc - 1; k >= 1; --k) {
                    pd[k] = p[k] - intq;
                }
            }
            pd[0] = psfc_raw - intq;

            for (int k = 0; k < nlev; ++k) {
                pd_gc(k, j, i) = pd[k];
            }

            // High-res surface pressure adjusted to model terrain via sfcprs3
            float psfc_hires = compute_single_sfcprs3(met.hgt_m(j, i), met.pmsl(j, i), p, ght);

            float mu0 = psfc_hires - intq - p_top;
            data.mu(j, i) = mu0 - data.mub(j, i);
        }
    }
}

static float lagrange_eval(const float* x, const float* y, int n, float target_x) {
    float px = 0.0f;
    for (int i = 0; i <= n; ++i) {
        float numer = 1.0f;
        float denom = 1.0f;
        for (int k = 0; k <= n; ++k) {
            if (k == i) continue;
            numer *= (target_x - x[k]);
            denom *= (x[i] - x[k]);
        }
        if (denom != 0.0f) {
            px += y[i] * numer / denom;
        }
    }
    return px;
}

static void interpolate_column(char var_type,
                               int interp_type, // 1: p, 2: log(p)
                               const std::vector<float>& all_x_in,
                               const std::vector<float>& all_y_in,
                               const std::vector<float>& target_x_in,
                               std::vector<float>& target_y_out) {
    const int all_dim = static_cast<int>(all_x_in.size());
    const int target_dim = static_cast<int>(target_x_in.size());
    target_y_out.resize(target_dim);

    std::vector<float> all_x = all_x_in;
    std::vector<float> target_x = target_x_in;
    if (interp_type == 2) {
        for (auto& val : all_x) val = std::log(val);
        for (auto& val : target_x) val = std::log(val);
    }

    const float rovcp = 287.0f / 1004.5f;
    const float crc1 = 11880.516f;
    const float crc2 = 0.1902632f;
    const float crc3 = 0.0065f;
    const int vboundb = 4;

    for (int t = 0; t < target_dim; ++t) {
        float tx = target_x[t];
        int loc_l = -1, loc_r = -1;
        bool found = false;

        for (int l = 0; l < all_dim - 1; ++l) {
            float a = tx - all_x[l];
            float b = tx - all_x[l + 1];
            if (a * b <= 0.0f) {
                loc_l = l;
                loc_r = l + 1;
                found = true;
                break;
            }
        }

        if (!found && target_x_in[t] > all_x_in[0]) {
            // Extrapolation below ground
            float p_target = target_x_in[t];
            float p_sfc = all_x_in[0];
            if (var_type == 'T') {
                float depth_p = p_target - p_sfc;
                float avg_p = 0.5f * (p_target + p_sfc);
                float temp_start = all_y_in[0] * std::pow(p_sfc / 100000.0f, rovcp);
                float dhdp = crc1 * crc2 * std::pow(avg_p / 100.0f, crc2 - 1.0f);
                float dh = dhdp * (depth_p / 100.0f);
                float dt = dh * crc3;
                target_y_out[t] = (temp_start + dt) * std::pow(100000.0f / p_target, rovcp);
            } else {
                target_y_out[t] = all_y_in[0];
            }
            continue;
        }

        if (!found) {
            // Above top of atmosphere, clamp to top
            target_y_out[t] = all_y_in.back();
            continue;
        }

        // Within bounds
        if (t < vboundb) {
            target_y_out[t] = lagrange_eval(&all_x[loc_l], &all_y_in[loc_l], 1, tx);
        } else {
            bool can_left = (loc_l - 1 >= 0 && loc_r + 1 <= all_dim);
            bool can_right = (loc_l >= 0 && loc_r + 2 <= all_dim);

            if (can_left && can_right) {
                float y1 = lagrange_eval(&all_x[loc_l], &all_y_in[loc_l], 2, tx);
                float y2 = lagrange_eval(&all_x[loc_l - 1], &all_y_in[loc_l - 1], 2, tx);
                target_y_out[t] = 0.5f * (y1 + y2);
            } else if (can_left) {
                target_y_out[t] = lagrange_eval(&all_x[loc_l - 1], &all_y_in[loc_l - 1], 2, tx);
            } else if (can_right) {
                target_y_out[t] = lagrange_eval(&all_x[loc_l], &all_y_in[loc_l], 2, tx);
            } else {
                target_y_out[t] = lagrange_eval(&all_x[loc_l], &all_y_in[loc_l], 1, tx);
            }
        }
    }
}

// Builds the sorted input column matching WRF vert_interp with force_sfc_in_vinterp = 1
static void build_ordered_column(const std::vector<float>& p_in,
                                 const std::vector<float>& f_in,
                                 float p_target_bottom,
                                 float zap_close,
                                 std::vector<float>& p_out,
                                 std::vector<float>& f_out) {
    const int generic = static_cast<int>(p_in.size()); // 34
    int ko_above_sfc = -1;
    for (int k = 1; k < generic; ++k) {
        if (p_in[0] > p_in[k]) {
            ko_above_sfc = k;
            break;
        }
    }

    p_out.clear();
    f_out.clear();

    if (ko_above_sfc > 1) {
        // Underground levels 1 .. ko_above_sfc - 1
        for (int k = 1; k < ko_above_sfc; ++k) {
            p_out.push_back(p_in[k]);
            f_out.push_back(f_in[k]);
        }
        if (!p_out.empty() && (p_out.back() - p_in[0] < zap_close)) {
            p_out.pop_back();
            f_out.pop_back();
        }
        p_out.push_back(p_in[0]);
        f_out.push_back(f_in[0]);

        // force_sfc_in_vinterp = 1
        int knext = ko_above_sfc;
        for (int k = ko_above_sfc; k < generic; ++k) {
            if (p_in[k] <= p_target_bottom) {
                knext = k;
                break;
            }
        }

        int kst = knext;
        if (kst < generic && (p_out.back() - p_in[kst] < zap_close)) {
            kst++;
        }
        for (int k = kst; k < generic; ++k) {
            p_out.push_back(p_in[k]);
            f_out.push_back(f_in[k]);
        }
    } else {
        // Surface is lowest level
        p_out.push_back(p_in[0]);
        f_out.push_back(f_in[0]);

        int knext = 1;
        for (int k = 1; k < generic; ++k) {
            if (p_in[k] <= p_target_bottom) {
                knext = k;
                break;
            }
        }
        int kst = knext;
        if (kst < generic && (p_out.back() - p_in[kst] < zap_close)) {
            kst++;
        }
        for (int k = kst; k < generic; ++k) {
            p_out.push_back(p_in[k]);
            f_out.push_back(f_in[k]);
        }
    }
}

void vertical_interpolate_fields(const MetEmData& met,
                                 const Array3D& pd_gc,
                                 const Array3D& qv_gc,
                                 WrfInputData& data,
                                 Array3D& rh_interp) {
    const int nx = data.dims.nx;
    const int ny = data.dims.ny;
    const int nz = data.dims.nz;
    const int nx_stag = data.dims.nx_stag;
    const int ny_stag = data.dims.ny_stag;
    const int nlev = met.dims.num_metgrid_levels;
    const float zap_close = 500.0f; // Pa

    data.u = Array3D(nz, ny, nx_stag);
    data.v = Array3D(nz, ny_stag, nx);
    data.t = Array3D(nz, ny, nx);
    data.qvapor = Array3D(nz, ny, nx);
    data.p = Array3D(nz, ny, nx);
    rh_interp = Array3D(nz, ny, nx);

    // Compute input theta from met.tt and met.pres
    Array3D theta_in(nlev, ny, nx);
    const float p00 = data.p00;
    const float rcp = 287.0f / 1004.5f;

    for (int k = 0; k < nlev; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                theta_in(k, j, i) = met.tt(k, j, i) * std::pow(p00 / met.pres(k, j, i), rcp);
            }
        }
    }

    // 1. Interpolate T and RH on mass points
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            float mu0 = data.mu(j, i) + data.mub(j, i);
            std::vector<float> p_target(nz);
            for (int k = 0; k < nz; ++k) {
                p_target[k] = data.c3h[k] * mu0 + data.c4h[k] + data.p_top;
            }

            std::vector<float> p_col(nlev), t_col(nlev), rh_col(nlev);
            for (int k = 0; k < nlev; ++k) {
                p_col[k] = pd_gc(k, j, i);
                t_col[k] = theta_in(k, j, i);
                rh_col[k] = met.rh(k, j, i);
            }

            std::vector<float> ord_p, ord_t, ord_rh;
            build_ordered_column(p_col, t_col, p_target[0], zap_close, ord_p, ord_t);
            build_ordered_column(p_col, rh_col, p_target[0], zap_close, ord_p, ord_rh);

            std::vector<float> t_out, rh_out;
            interpolate_column('T', 2, ord_p, ord_t, p_target, t_out);
            interpolate_column('Q', 2, ord_p, ord_rh, p_target, rh_out);

            for (int k = 0; k < nz; ++k) {
                data.t(k, j, i) = t_out[k] - 300.0f;
                rh_interp(k, j, i) = rh_out[k];
                float t_k = t_out[k] * std::pow(p_target[k] / p00, rcp);
                data.qvapor(k, j, i) = rh_to_mxrat_single(rh_out[k], t_k, p_target[k]);
            }
        }
    }

    // 2. Interpolate U on staggered west_east grid [nz, ny, nx_stag]
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx_stag; ++i) {
            float mu0_l = (i == 0) ? (data.mu(j, 0) + data.mub(j, 0)) : (data.mu(j, i - 1) + data.mub(j, i - 1));
            float mu0_r = (i == nx) ? (data.mu(j, nx - 1) + data.mub(j, nx - 1)) : (data.mu(j, i) + data.mub(j, i));
            float mu0_u = (i == 0) ? mu0_l : (i == nx) ? mu0_r : 0.5f * (mu0_l + mu0_r);

            std::vector<float> p_target(nz);
            for (int k = 0; k < nz; ++k) {
                p_target[k] = data.c3h[k] * mu0_u + data.c4h[k] + data.p_top;
            }

            std::vector<float> p_col(nlev), u_col(nlev);
            for (int k = 0; k < nlev; ++k) {
                if (i == 0) {
                    p_col[k] = pd_gc(k, j, 0);
                } else if (i == nx) {
                    p_col[k] = pd_gc(k, j, nx - 1);
                } else {
                    p_col[k] = 0.5f * (pd_gc(k, j, i) + pd_gc(k, j, i - 1));
                }
                u_col[k] = met.uu(k, j, i);
            }

            std::vector<float> ord_p, ord_u;
            build_ordered_column(p_col, u_col, p_target[0], zap_close, ord_p, ord_u);

            std::vector<float> u_out;
            interpolate_column('U', 1, ord_p, ord_u, p_target, u_out);
            for (int k = 0; k < nz; ++k) {
                data.u(k, j, i) = u_out[k];
            }
        }
    }

    // 3. Interpolate V on staggered south_north grid [nz, ny_stag, nx]
    for (int j = 0; j < ny_stag; ++j) {
        for (int i = 0; i < nx; ++i) {
            float mu0_d = (j == 0) ? (data.mu(0, i) + data.mub(0, i)) : (data.mu(j - 1, i) + data.mub(j - 1, i));
            float mu0_u = (j == ny) ? (data.mu(ny - 1, i) + data.mub(ny - 1, i)) : (data.mu(j, i) + data.mub(j, i));
            float mu0_v = (j == 0) ? mu0_d : (j == ny) ? mu0_u : 0.5f * (mu0_d + mu0_u);

            std::vector<float> p_target(nz);
            for (int k = 0; k < nz; ++k) {
                p_target[k] = data.c3h[k] * mu0_v + data.c4h[k] + data.p_top;
            }

            std::vector<float> p_col(nlev), v_col(nlev);
            for (int k = 0; k < nlev; ++k) {
                if (j == 0) {
                    p_col[k] = pd_gc(k, 0, i);
                } else if (j == ny) {
                    p_col[k] = pd_gc(k, ny - 1, i);
                } else {
                    p_col[k] = 0.5f * (pd_gc(k, j, i) + pd_gc(k, j - 1, i));
                }
                v_col[k] = met.vv(k, j, i);
            }

            std::vector<float> ord_p, ord_v;
            build_ordered_column(p_col, v_col, p_target[0], zap_close, ord_p, ord_v);

            std::vector<float> v_out;
            interpolate_column('V', 1, ord_p, ord_v, p_target, v_out);
            for (int k = 0; k < nz; ++k) {
                data.v(k, j, i) = v_out[k];
            }
        }
    }
}

} // namespace cppwrf::init_real
