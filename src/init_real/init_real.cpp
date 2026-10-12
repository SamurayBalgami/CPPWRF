#include "init_real.hpp"
#include <iostream>
#include <stdexcept>

namespace cppwrf::init_real {

static WrfInputData process_single_time(const MetEmData& met, int hybrid_opt) {
    WrfInputData data;
    data.dims = met.dims;
    data.dims.nz = 44;
    data.dims.nz_stag = 45;
    data.hybrid_opt = hybrid_opt;
    data.p_top = 5000.0f;
    data.etac = 0.2f;

    // 1. Vertical coordinate eta levels
    compute_eta_levels(data);

    // 2. Hybrid coefficients c1..c4
    compute_vcoord_coeffs(data);

    // 3. Copy 2D map factors and static grid fields
    data.hgt = met.hgt_m;
    data.msft = met.mapfac_m;
    data.msfu = met.mapfac_u;
    data.msfv = met.mapfac_v;
    data.msftx = met.mapfac_mx;
    data.msfty = met.mapfac_my;
    data.msfux = met.mapfac_ux;
    data.msfuy = met.mapfac_uy;
    data.msfvx = met.mapfac_vx;
    data.msfvy = met.mapfac_vy;

    data.f = met.f;
    data.e = met.e;
    data.sina = met.sinalpha;
    data.cosa = met.cosalpha;

    data.xlat = met.xlat_m;
    data.xlong = met.xlong_m;
    data.xlat_u = met.xlat_u;
    data.xlong_u = met.xlong_u;
    data.xlat_v = met.xlat_v;
    data.xlong_v = met.xlong_v;
    data.clat = met.clat;
    data.clong = met.clong;

    // 4. Base state
    compute_base_state(data, 2);

    // 5. Dry surface pressure and column mass mu
    Array3D pd_gc, qv_gc;
    compute_dry_pressure_and_mu(met, data, pd_gc, qv_gc);

    // 6. Vertical interpolation of U, V, T, QV
    Array3D rh_interp;
    vertical_interpolate_fields(met, pd_gc, qv_gc, data, rh_interp);

    // 7. Hydrostatic integration of geopotential ph, diagnosis of al, p, and final qvapor
    integrate_hydrostatic(met, data, rh_interp, 2);

    return data;
}

RealInitResult run_real_init(const std::vector<std::string>& met_em_files,
                             int hybrid_opt) {
    if (met_em_files.empty()) {
        throw std::runtime_error("No met_em files provided to run_real_init");
    }

    RealInitResult res;
    std::vector<WrfInputData> all_states;
    std::vector<std::string> date_strings;

    for (const auto& file : met_em_files) {
        MetEmData met = read_met_em(file);
        WrfInputData state = process_single_time(met, hybrid_opt);
        all_states.push_back(state);
        date_strings.push_back(met.time_str);
    }

    res.initial_state = all_states[0];

    if (all_states.size() >= 2) {
        res.bdy_data = build_wrfbdy_data(all_states, date_strings, 21600.0f);
    }

    return res;
}

} // namespace cppwrf::init_real
