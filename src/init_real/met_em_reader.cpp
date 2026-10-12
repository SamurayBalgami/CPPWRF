#include "met_em_reader.hpp"
#include <netcdf.h>
#include <stdexcept>
#include <iostream>
#include <vector>

namespace cppwrf::init_real {

#define NC_CHECK(cmd) do { \
    int _e = (cmd); \
    if (_e != NC_NOERR) { \
        throw std::runtime_error(std::string("NetCDF error: ") + nc_strerror(_e) + " at " + __FILE__ + ":" + std::to_string(__LINE__)); \
    } \
} while(0)

static void read_var_2d(int ncid, const char* varname, Array2D& arr) {
    int varid;
    NC_CHECK(nc_inq_varid(ncid, varname, &varid));
    size_t start[3] = {0, 0, 0};
    size_t count[3] = {1, static_cast<size_t>(arr.d1), static_cast<size_t>(arr.d2)};
    NC_CHECK(nc_get_vara_float(ncid, varid, start, count, arr.data_ptr()));
}

static void read_var_3d(int ncid, const char* varname, Array3D& arr) {
    int varid;
    NC_CHECK(nc_inq_varid(ncid, varname, &varid));
    size_t start[4] = {0, 0, 0, 0};
    size_t count[4] = {1, static_cast<size_t>(arr.d1), static_cast<size_t>(arr.d2), static_cast<size_t>(arr.d3)};
    NC_CHECK(nc_get_vara_float(ncid, varid, start, count, arr.data_ptr()));
}

MetEmData read_met_em(const std::string& filepath) {
    int ncid;
    NC_CHECK(nc_open(filepath.c_str(), NC_NOWRITE, &ncid));

    MetEmData met;
    int dimid_we, dimid_sn, dimid_lev, dimid_we_stag, dimid_sn_stag;
    size_t len_we, len_sn, len_lev, len_we_stag, len_sn_stag;

    NC_CHECK(nc_inq_dimid(ncid, "west_east", &dimid_we));
    NC_CHECK(nc_inq_dimlen(ncid, dimid_we, &len_we));
    NC_CHECK(nc_inq_dimid(ncid, "south_north", &dimid_sn));
    NC_CHECK(nc_inq_dimlen(ncid, dimid_sn, &len_sn));
    NC_CHECK(nc_inq_dimid(ncid, "num_metgrid_levels", &dimid_lev));
    NC_CHECK(nc_inq_dimlen(ncid, dimid_lev, &len_lev));
    NC_CHECK(nc_inq_dimid(ncid, "west_east_stag", &dimid_we_stag));
    NC_CHECK(nc_inq_dimlen(ncid, dimid_we_stag, &len_we_stag));
    NC_CHECK(nc_inq_dimid(ncid, "south_north_stag", &dimid_sn_stag));
    NC_CHECK(nc_inq_dimlen(ncid, dimid_sn_stag, &len_sn_stag));

    met.dims.nx = static_cast<int>(len_we);
    met.dims.ny = static_cast<int>(len_sn);
    met.dims.num_metgrid_levels = static_cast<int>(len_lev);
    met.dims.nx_stag = static_cast<int>(len_we_stag);
    met.dims.ny_stag = static_cast<int>(len_sn_stag);

    // Read date string
    int varid_times;
    if (nc_inq_varid(ncid, "Times", &varid_times) == NC_NOERR) {
        char time_buf[64] = {0};
        size_t t_start[2] = {0, 0};
        size_t t_count[2] = {1, 19};
        nc_get_vara_text(ncid, varid_times, t_start, t_count, time_buf);
        met.time_str = std::string(time_buf, 19);
    }

    const int nx = met.dims.nx;
    const int ny = met.dims.ny;
    const int nx_stag = met.dims.nx_stag;
    const int ny_stag = met.dims.ny_stag;
    const int nlev = met.dims.num_metgrid_levels;

    met.xlat_m = Array2D(ny, nx);
    met.xlong_m = Array2D(ny, nx);
    met.xlat_u = Array2D(ny, nx_stag);
    met.xlong_u = Array2D(ny, nx_stag);
    met.xlat_v = Array2D(ny_stag, nx);
    met.xlong_v = Array2D(ny_stag, nx);
    met.clat = Array2D(ny, nx);
    met.clong = Array2D(ny, nx);

    met.mapfac_m = Array2D(ny, nx);
    met.mapfac_u = Array2D(ny, nx_stag);
    met.mapfac_v = Array2D(ny_stag, nx);
    met.mapfac_mx = Array2D(ny, nx);
    met.mapfac_my = Array2D(ny, nx);
    met.mapfac_ux = Array2D(ny, nx_stag);
    met.mapfac_uy = Array2D(ny, nx_stag);
    met.mapfac_vx = Array2D(ny_stag, nx);
    met.mapfac_vy = Array2D(ny_stag, nx);

    met.f = Array2D(ny, nx);
    met.e = Array2D(ny, nx);
    met.sinalpha = Array2D(ny, nx);
    met.cosalpha = Array2D(ny, nx);

    met.hgt_m = Array2D(ny, nx);
    met.psfc = Array2D(ny, nx);
    met.pmsl = Array2D(ny, nx);

    read_var_2d(ncid, "XLAT_M", met.xlat_m);
    read_var_2d(ncid, "XLONG_M", met.xlong_m);
    read_var_2d(ncid, "XLAT_U", met.xlat_u);
    read_var_2d(ncid, "XLONG_U", met.xlong_u);
    read_var_2d(ncid, "XLAT_V", met.xlat_v);
    read_var_2d(ncid, "XLONG_V", met.xlong_v);
    read_var_2d(ncid, "CLAT", met.clat);
    read_var_2d(ncid, "CLONG", met.clong);

    read_var_2d(ncid, "MAPFAC_M", met.mapfac_m);
    read_var_2d(ncid, "MAPFAC_U", met.mapfac_u);
    read_var_2d(ncid, "MAPFAC_V", met.mapfac_v);
    read_var_2d(ncid, "MAPFAC_MX", met.mapfac_mx);
    read_var_2d(ncid, "MAPFAC_MY", met.mapfac_my);
    read_var_2d(ncid, "MAPFAC_UX", met.mapfac_ux);
    read_var_2d(ncid, "MAPFAC_UY", met.mapfac_uy);
    read_var_2d(ncid, "MAPFAC_VX", met.mapfac_vx);
    read_var_2d(ncid, "MAPFAC_VY", met.mapfac_vy);

    read_var_2d(ncid, "F", met.f);
    read_var_2d(ncid, "E", met.e);
    read_var_2d(ncid, "SINALPHA", met.sinalpha);
    read_var_2d(ncid, "COSALPHA", met.cosalpha);

    read_var_2d(ncid, "HGT_M", met.hgt_m);
    read_var_2d(ncid, "PSFC", met.psfc);
    read_var_2d(ncid, "PMSL", met.pmsl);

    // 3D variables
    met.pres = Array3D(nlev, ny, nx);
    met.ght = Array3D(nlev, ny, nx);
    met.tt = Array3D(nlev, ny, nx);
    met.rh = Array3D(nlev, ny, nx);
    met.uu = Array3D(nlev, ny, nx_stag);
    met.vv = Array3D(nlev, ny_stag, nx);

    read_var_3d(ncid, "PRES", met.pres);
    read_var_3d(ncid, "GHT", met.ght);
    read_var_3d(ncid, "TT", met.tt);
    read_var_3d(ncid, "RH", met.rh);
    read_var_3d(ncid, "UU", met.uu);
    read_var_3d(ncid, "VV", met.vv);

    NC_CHECK(nc_close(ncid));
    return met;
}

} // namespace cppwrf::init_real
