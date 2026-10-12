#pragma once

#include "real_types.hpp"
#include <string>

namespace cppwrf::init_real {

struct MetEmData {
    GridDims dims;
    std::string time_str;

    // 2D fields [ny, nx] or staggered
    Array2D xlat_m;      // [ny, nx]
    Array2D xlong_m;     // [ny, nx]
    Array2D xlat_u;      // [ny, nx_stag]
    Array2D xlong_u;     // [ny, nx_stag]
    Array2D xlat_v;      // [ny_stag, nx]
    Array2D xlong_v;     // [ny_stag, nx]
    Array2D clat;        // [ny, nx]
    Array2D clong;       // [ny, nx]

    Array2D mapfac_m;    // [ny, nx]
    Array2D mapfac_u;    // [ny, nx_stag]
    Array2D mapfac_v;    // [ny_stag, nx]
    Array2D mapfac_mx;   // [ny, nx]
    Array2D mapfac_my;   // [ny, nx]
    Array2D mapfac_ux;   // [ny, nx_stag]
    Array2D mapfac_uy;   // [ny, nx_stag]
    Array2D mapfac_vx;   // [ny_stag, nx]
    Array2D mapfac_vy;   // [ny_stag, nx]

    Array2D f;           // [ny, nx]
    Array2D e;           // [ny, nx]
    Array2D sinalpha;    // [ny, nx]
    Array2D cosalpha;    // [ny, nx]

    Array2D hgt_m;       // [ny, nx]
    Array2D psfc;        // [ny, nx]
    Array2D pmsl;        // [ny, nx]

    // 3D fields [num_metgrid_levels, ny, nx]
    Array3D pres;        // [num_levels, ny, nx] Pa
    Array3D ght;         // [num_levels, ny, nx] m
    Array3D tt;          // [num_levels, ny, nx] K
    Array3D rh;          // [num_levels, ny, nx] %
    Array3D uu;          // [num_levels, ny, nx_stag] m/s
    Array3D vv;          // [num_levels, ny_stag, nx] m/s
};

MetEmData read_met_em(const std::string& filepath);

} // namespace cppwrf::init_real
