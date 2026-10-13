#pragma once

#include <cstddef>
#include <netcdf.h>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace cppwrf::io_netcdf {

#define NC_CHECK(call)                                                                     \
  do {                                                                                     \
    int _nc_err = (call);                                                                  \
    if (_nc_err != NC_NOERR) {                                                             \
      throw std::runtime_error(std::string("NetCDF error: ") + nc_strerror(_nc_err) +       \
                               " at " + __FILE__ + ":" + std::to_string(__LINE__));        \
    }                                                                                      \
  } while (0)

enum class DType { F32, F64 };

struct VarView {
  std::string name;
  char stagger = '\0';   // 'M' or '\0' (unstaggered), 'U' (X), 'V' (Y), 'W' (Z)
  int nx = 0, ny = 0, nz = 0; // logical extents (1 for 1D/2D if dimension not used)
  const void* data = nullptr; // pointer to host memory (contiguous or strided)
  DType dtype = DType::F32;
  std::ptrdiff_t sj = 0, sk = 0; // element strides for j and k; i stride is always 1
  std::string description;
  std::string units;

  VarView() = default;
  VarView(std::string n, char stag, int x, int y, int z, const void* d, DType dt,
          std::ptrdiff_t s_j = 0, std::ptrdiff_t s_k = 0, std::string desc = "", std::string u = "")
      : name(std::move(n)), stagger(stag), nx(x), ny(y), nz(z), data(d), dtype(dt),
        sj(s_j), sk(s_k), description(std::move(desc)), units(std::move(u)) {}
};

struct WrfoutMeta {
  int nx = 0, ny = 0, nz = 0;
  double dx = 0.0, dy = 0.0, dt = 0.0;
  double ztop = 0.0, p_top = 0.0;
  int map_proj = 0; // 0 = Cartesian, 1 = Lambert Conformal, etc.
  double cen_lat = 0.0, cen_lon = 0.0;
  double truelat1 = 0.0, truelat2 = 0.0;
  double moad_cen_lat = 0.0, stand_lon = 0.0;
  double pole_lat = 90.0, pole_lon = 0.0;
  std::string start_date = "2000-01-01_00:00:00";
  std::string simulation_start_date = "2000-01-01_00:00:00";
  // physics/dynamics params
  int diff_opt = 1;
  int km_opt = 1;
  int damp_opt = 0;
  double dampcoef = 0.0;
  double khdif = 0.0;
  double kvdif = 0.0;
  int mp_physics = 0;
  int ra_lw_physics = 0;
  int ra_sw_physics = 0;
  int sf_sfclay_physics = 0;
  int sf_surface_physics = 0;
  int bl_pbl_physics = 0;
  int cu_physics = 0;
  int hybrid_opt = 0;
  int hypsometric_opt = 1;
  double etac = 0.2;
};

// Formats seconds into WRF 19-char timestamp "YYYY-MM-DD_hh:mm:ss" starting from base_date
std::string format_wrf_date(const std::string& base_date, double time_sec);

}  // namespace cppwrf::io_netcdf
