#include "io_netcdf/restart.hpp"

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <vector>

namespace cppwrf::io_netcdf {

namespace {

void get_var_dims(int stagger, int nz, int dim_time, int dim_bt, int dim_bt_stag,
                  int dim_sn, int dim_sn_stag, int dim_we, int dim_we_stag,
                  std::vector<int>& dimids) {
  dimids.clear();
  dimids.push_back(dim_time);
  if (nz > 1) {
    if (stagger == 'W' || stagger == 'Z') dimids.push_back(dim_bt_stag);
    else dimids.push_back(dim_bt);
  }
  if (stagger == 'V' || stagger == 'Y') dimids.push_back(dim_sn_stag);
  else dimids.push_back(dim_sn);

  if (stagger == 'U' || stagger == 'X') dimids.push_back(dim_we_stag);
  else dimids.push_back(dim_we);
}

}  // namespace

void write_restart(const std::string& path, const RestartMeta& meta,
                   const std::vector<VarView>& vars) {
  auto p = std::filesystem::path(path).parent_path();
  if (!p.empty()) std::filesystem::create_directories(p);

  int ncid = -1;
  NC_CHECK(nc_create(path.c_str(), NC_CLOBBER | NC_64BIT_OFFSET, &ncid));

  // Global attributes
  NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "TITLE", 25, "RESTART FILE FROM CPPWRF"));
  int flag_restart = 1;
  NC_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "FLAG_RESTART", NC_INT, 1, &flag_restart));
  int step_i = int(meta.step);
  NC_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "STEP", NC_INT, 1, &step_i));
  NC_CHECK(nc_put_att_double(ncid, NC_GLOBAL, "TIME", NC_DOUBLE, 1, &meta.time_sec));

  NC_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "NX", NC_INT, 1, &meta.nx));
  NC_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "NY", NC_INT, 1, &meta.ny));
  NC_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "NZ", NC_INT, 1, &meta.nz));

  NC_CHECK(nc_put_att_double(ncid, NC_GLOBAL, "DX", NC_DOUBLE, 1, &meta.dx));
  NC_CHECK(nc_put_att_double(ncid, NC_GLOBAL, "DY", NC_DOUBLE, 1, &meta.dy));
  NC_CHECK(nc_put_att_double(ncid, NC_GLOBAL, "DT", NC_DOUBLE, 1, &meta.dt));

  if (!meta.case_name.empty()) {
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "CASE", meta.case_name.size(), meta.case_name.data()));
  }
  if (!meta.date_str.empty()) {
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "CURRENT_DATE", meta.date_str.size(), meta.date_str.data()));
  }

  // Dimensions
  int dim_time = -1, dim_bt = -1, dim_bt_stag = -1;
  int dim_sn = -1, dim_sn_stag = -1, dim_we = -1, dim_we_stag = -1;

  NC_CHECK(nc_def_dim(ncid, "Time", 1, &dim_time));
  NC_CHECK(nc_def_dim(ncid, "bottom_top", meta.nz, &dim_bt));
  NC_CHECK(nc_def_dim(ncid, "bottom_top_stag", meta.nz + 1, &dim_bt_stag));
  NC_CHECK(nc_def_dim(ncid, "south_north", meta.ny, &dim_sn));
  NC_CHECK(nc_def_dim(ncid, "south_north_stag", meta.ny + 1, &dim_sn_stag));
  NC_CHECK(nc_def_dim(ncid, "west_east", meta.nx, &dim_we));
  NC_CHECK(nc_def_dim(ncid, "west_east_stag", meta.nx + 1, &dim_we_stag));

  // Variable definitions
  std::vector<int> varids(vars.size(), -1);
  for (size_t i = 0; i < vars.size(); ++i) {
    const auto& v = vars[i];
    std::vector<int> dimids;
    get_var_dims(v.stagger, v.nz, dim_time, dim_bt, dim_bt_stag,
                 dim_sn, dim_sn_stag, dim_we, dim_we_stag, dimids);
    nc_type xtype = (v.dtype == DType::F64) ? NC_DOUBLE : NC_FLOAT;
    NC_CHECK(nc_def_var(ncid, v.name.c_str(), xtype, int(dimids.size()), dimids.data(), &varids[i]));
    char stag_str[2] = {v.stagger ? v.stagger : '\0', '\0'};
    NC_CHECK(nc_put_att_text(ncid, varids[i], "stagger", 1, stag_str));
    int shape[3] = {v.nx, v.ny, v.nz};
    NC_CHECK(nc_put_att_int(ncid, varids[i], "shape", NC_INT, 3, shape));
  }

  NC_CHECK(nc_enddef(ncid));

  // Write data
  for (size_t i = 0; i < vars.size(); ++i) {
    const auto& v = vars[i];
    int varid = varids[i];
    const int nx = v.nx > 0 ? v.nx : 1;
    const int ny = v.ny > 0 ? v.ny : 1;
    const int nz = v.nz > 0 ? v.nz : 1;
    size_t count = size_t(nx) * ny * nz;

    const std::ptrdiff_t sj = v.sj != 0 ? v.sj : nx;
    const std::ptrdiff_t sk = v.sk != 0 ? v.sk : (std::ptrdiff_t)nx * ny;
    bool is_contiguous = (sj == nx && sk == (std::ptrdiff_t)nx * ny);

    if (v.dtype == DType::F64) {
      const double* src = static_cast<const double*>(v.data);
      if (is_contiguous) {
        NC_CHECK(nc_put_var_double(ncid, varid, src));
      } else {
        std::vector<double> buf(count);
        for (int k = 0; k < nz; ++k)
          for (int j = 0; j < ny; ++j)
            for (int ii = 0; ii < nx; ++ii)
              buf[(size_t(k) * ny + j) * nx + ii] = src[ii + j * sj + k * sk];
        NC_CHECK(nc_put_var_double(ncid, varid, buf.data()));
      }
    } else {
      const float* src = static_cast<const float*>(v.data);
      if (is_contiguous) {
        NC_CHECK(nc_put_var_float(ncid, varid, src));
      } else {
        std::vector<float> buf(count);
        for (int k = 0; k < nz; ++k)
          for (int j = 0; j < ny; ++j)
            for (int ii = 0; ii < nx; ++ii)
              buf[(size_t(k) * ny + j) * nx + ii] = src[ii + j * sj + k * sk];
        NC_CHECK(nc_put_var_float(ncid, varid, buf.data()));
      }
    }
  }

  NC_CHECK(nc_close(ncid));
}

void read_restart(const std::string& path, RestartMeta& meta,
                  std::vector<VarView>& vars) {
  int ncid = -1;
  NC_CHECK(nc_open(path.c_str(), NC_NOWRITE, &ncid));

  // Global metadata
  int step_i = 0;
  if (nc_get_att_int(ncid, NC_GLOBAL, "STEP", &step_i) == NC_NOERR) {
    meta.step = step_i;
  }
  NC_CHECK(nc_get_att_double(ncid, NC_GLOBAL, "TIME", &meta.time_sec));
  NC_CHECK(nc_get_att_int(ncid, NC_GLOBAL, "NX", &meta.nx));
  NC_CHECK(nc_get_att_int(ncid, NC_GLOBAL, "NY", &meta.ny));
  NC_CHECK(nc_get_att_int(ncid, NC_GLOBAL, "NZ", &meta.nz));
  NC_CHECK(nc_get_att_double(ncid, NC_GLOBAL, "DX", &meta.dx));
  NC_CHECK(nc_get_att_double(ncid, NC_GLOBAL, "DY", &meta.dy));
  NC_CHECK(nc_get_att_double(ncid, NC_GLOBAL, "DT", &meta.dt));

  size_t att_len = 0;
  if (nc_inq_attlen(ncid, NC_GLOBAL, "CASE", &att_len) == NC_NOERR && att_len > 0) {
    meta.case_name.resize(att_len);
    nc_get_att_text(ncid, NC_GLOBAL, "CASE", meta.case_name.data());
  }
  if (nc_inq_attlen(ncid, NC_GLOBAL, "CURRENT_DATE", &att_len) == NC_NOERR && att_len > 0) {
    meta.date_str.resize(att_len);
    nc_get_att_text(ncid, NC_GLOBAL, "CURRENT_DATE", meta.date_str.data());
  }

  // Read each requested variable
  for (auto& v : vars) {
    int varid = -1;
    NC_CHECK(nc_inq_varid(ncid, v.name.c_str(), &varid));

    const int nx = v.nx > 0 ? v.nx : 1;
    const int ny = v.ny > 0 ? v.ny : 1;
    const int nz = v.nz > 0 ? v.nz : 1;
    size_t count = size_t(nx) * ny * nz;

    const std::ptrdiff_t sj = v.sj != 0 ? v.sj : nx;
    const std::ptrdiff_t sk = v.sk != 0 ? v.sk : (std::ptrdiff_t)nx * ny;
    bool is_contiguous = (sj == nx && sk == (std::ptrdiff_t)nx * ny);

    if (v.dtype == DType::F64) {
      double* dst = static_cast<double*>(const_cast<void*>(v.data));
      if (is_contiguous) {
        NC_CHECK(nc_get_var_double(ncid, varid, dst));
      } else {
        std::vector<double> buf(count);
        NC_CHECK(nc_get_var_double(ncid, varid, buf.data()));
        for (int k = 0; k < nz; ++k)
          for (int j = 0; j < ny; ++j)
            for (int ii = 0; ii < nx; ++ii)
              dst[ii + j * sj + k * sk] = buf[(size_t(k) * ny + j) * nx + ii];
      }
    } else {
      float* dst = static_cast<float*>(const_cast<void*>(v.data));
      if (is_contiguous) {
        NC_CHECK(nc_get_var_float(ncid, varid, dst));
      } else {
        std::vector<float> buf(count);
        NC_CHECK(nc_get_var_float(ncid, varid, buf.data()));
        for (int k = 0; k < nz; ++k)
          for (int j = 0; j < ny; ++j)
            for (int ii = 0; ii < nx; ++ii)
              dst[ii + j * sj + k * sk] = buf[(size_t(k) * ny + j) * nx + ii];
      }
    }
  }

  NC_CHECK(nc_close(ncid));
}

}  // namespace cppwrf::io_netcdf
