#pragma once

#include <map>
#include <string>
#include <vector>

namespace cppwrf::io_netcdf {

/**
 * Data structure holding all required dycore initialization fields from wrfinput.
 *
 * C++ Index Order Documentation:
 * In WRF NetCDF files, multidimensional arrays are stored in row-major order:
 * - 3D Mass fields: (Time, bottom_top, south_north, west_east) -> shape (nz, ny, nx).
 *   In C++ memory, index i (west_east) varies fastest (stride 1),
 *   index j (south_north) varies next (stride nx),
 *   index k (bottom_top) varies slowest (stride nx * ny).
 *   Flat index: ((k * ny) + j) * nx + i
 *   C++ vertical index k in [0, nz-1] corresponds to Fortran level k_fortran - 1 (1 .. nz).
 *
 * - 3D Staggered U: (Time, bottom_top, south_north, west_east_stag) -> shape (nz, ny, nx + 1).
 *   Flat index: ((k * ny) + j) * (nx + 1) + i
 *
 * - 3D Staggered V: (Time, bottom_top, south_north_stag, west_east) -> shape (nz, ny + 1, nx).
 *   Flat index: ((k * (ny + 1)) + j) * nx + i
 *
 * - 3D Staggered W and Geopotential (PH, PHB):
 *   (Time, bottom_top_stag, south_north, west_east) -> shape (nz + 1, ny, nx).
 *   Flat index: ((k * ny) + j) * nx + i
 *
 * - 2D Mass fields: (Time, south_north, west_east) -> shape (ny, nx).
 *   Flat index: j * nx + i
 */
struct WrfInputData {
  // Grid dimensions
  int nx = 0;         // west_east (mass grid)
  int ny = 0;         // south_north (mass grid)
  int nz = 0;         // bottom_top (mass grid)
  int nx_stag = 0;    // west_east_stag (nx + 1)
  int ny_stag = 0;    // south_north_stag (ny + 1)
  int nz_stag = 0;    // bottom_top_stag (nz + 1)

  // Spatial resolution and time step
  double dx = 0.0;
  double dy = 0.0;
  double dt = 0.0;

  // Global attributes
  int hybrid_opt = 0;
  double etac = 0.0;
  double p_top = 0.0;
  int map_proj = 0;
  double cen_lat = 0.0;
  double cen_lon = 0.0;
  double truelat1 = 0.0;
  double truelat2 = 0.0;
  double moad_cen_lat = 0.0;
  double stand_lon = 0.0;
  std::string start_date;
  std::string simulation_start_date;
  std::map<std::string, double> attrs_real;
  std::map<std::string, int> attrs_int;
  std::map<std::string, std::string> attrs_str;

  // 1D Vertical fields
  std::vector<double> znw;   // [nz + 1] eta on full (w) levels
  std::vector<double> znu;   // [nz]     eta on half (mass) levels
  std::vector<double> dn;    // [nz]
  std::vector<double> dnw;   // [nz]
  std::vector<double> rdn;   // [nz]
  std::vector<double> rdnw;  // [nz]
  std::vector<double> fnm;   // [nz]
  std::vector<double> fnp;   // [nz]

  // Hybrid coordinate weights
  std::vector<double> c1h, c2h, c3h, c4h;  // [nz]
  std::vector<double> c1f, c2f, c3f, c4f;  // [nz + 1]

  // Coordinate / finite-difference constants (from variables)
  double cf1 = 0.0;
  double cf2 = 0.0;
  double cf3 = 0.0;
  double cfn = 0.0;
  double cfn1 = 0.0;

  // 2D Horizontal fields (mass grid: size ny * nx, flat index: j * nx + i)
  std::vector<double> mub;        // [ny * nx]
  std::vector<double> mu;         // [ny * nx]
  std::vector<double> hgt;        // [ny * nx]
  std::vector<double> xlat;       // [ny * nx]
  std::vector<double> xlong;      // [ny * nx]
  std::vector<double> mapfac_mx;  // [ny * nx]
  std::vector<double> mapfac_my;  // [ny * nx]
  std::vector<double> f;          // [ny * nx] Coriolis f
  std::vector<double> e;          // [ny * nx] Coriolis e
  std::vector<double> sinalpha;   // [ny * nx]
  std::vector<double> cosalpha;   // [ny * nx]

  // Staggered horizontal map factors
  std::vector<double> mapfac_ux;  // [ny * nx_stag] (j * nx_stag + i)
  std::vector<double> mapfac_uy;  // [ny * nx_stag]
  std::vector<double> mapfac_vx;  // [ny_stag * nx] (j * nx + i)
  std::vector<double> mapfac_vy;  // [ny_stag * nx]

  // 3D Mass grid fields: size nz * ny * nx
  std::vector<double> pb;       // Base pressure [Pa]
  std::vector<double> p;        // Perturbation pressure [Pa]
  std::vector<double> t;        // Perturbation potential temperature (theta - 300 K)
  std::vector<double> t_init;   // Initial potential temperature
  std::vector<double> alb;      // Base inverse density
  std::vector<double> al;       // Perturbation inverse density
  std::vector<double> qvapor;   // Water vapor mixing ratio [kg/kg]

  // 3D Staggered fields
  std::vector<double> u;        // [nz * ny * nx_stag]
  std::vector<double> v;        // [nz * ny_stag * nx]
  std::vector<double> w;        // [nz_stag * ny * nx]
  std::vector<double> phb;      // [nz_stag * ny * nx]
  std::vector<double> ph;       // [nz_stag * ny * nx]

  // Other moist species if present
  std::vector<double> qcloud;
  std::vector<double> qrain;
  std::vector<double> qice;
  std::vector<double> qsnow;
  std::vector<double> qgraup;
  std::map<std::string, std::vector<double>> moist_species;

  // Indexing helper functions
  inline size_t idx_m2(int i, int j) const { return size_t(j) * nx + i; }
  inline size_t idx_u2(int i, int j) const { return size_t(j) * nx_stag + i; }
  inline size_t idx_v2(int i, int j) const { return size_t(j) * nx + i; }
  inline size_t idx_m3(int i, int j, int k) const { return (size_t(k) * ny + j) * nx + i; }
  inline size_t idx_u3(int i, int j, int k) const { return (size_t(k) * ny + j) * nx_stag + i; }
  inline size_t idx_v3(int i, int j, int k) const { return (size_t(k) * ny_stag + j) * nx + i; }
  inline size_t idx_w3(int i, int j, int k) const { return (size_t(k) * ny + j) * nx + i; }
};

/**
 * Data structure holding lateral boundary condition data from wrfbdy.
 *
 * Boundary Dimensions:
 * - Time: n_times
 * - bdy_width: width of boundary zone (typically 5)
 * - West (XS) and East (XE):
 *     3D: (Time, bdy_width, nz, ny)       -> flat: ((t * bdy_width + b) * nz + k) * ny + j
 *     2D: (Time, bdy_width, ny)           -> flat: (t * bdy_width + b) * ny + j
 *     Note for V: ny is replaced by ny_stag.
 *     Note for W, PH: nz is replaced by nz_stag.
 * - South (YS) and North (YE):
 *     3D: (Time, bdy_width, nz, nx)       -> flat: ((t * bdy_width + b) * nz + k) * nx + i
 *     2D: (Time, bdy_width, nx)           -> flat: (t * bdy_width + b) * nx + i
 *     Note for U: nx is replaced by nx_stag.
 *     Note for W, PH: nz is replaced by nz_stag.
 */
struct WrfBdyData {
  int n_times = 0;
  int bdy_width = 0;
  int nx = 0, ny = 0, nz = 0;
  int nx_stag = 0, ny_stag = 0, nz_stag = 0;

  std::vector<std::string> times;
  std::vector<std::string> this_bdy_times;
  std::vector<std::string> next_bdy_times;

  // 3D Boundary fields and tendencies
  // U
  std::vector<double> u_bxs, u_bxe, u_bys, u_bye;
  std::vector<double> u_btxs, u_btxe, u_btys, u_btye;
  // V
  std::vector<double> v_bxs, v_bxe, v_bys, v_bye;
  std::vector<double> v_btxs, v_btxe, v_btys, v_btye;
  // W
  std::vector<double> w_bxs, w_bxe, w_bys, w_bye;
  std::vector<double> w_btxs, w_btxe, w_btys, w_btye;
  // PH
  std::vector<double> ph_bxs, ph_bxe, ph_bys, ph_bye;
  std::vector<double> ph_btxs, ph_btxe, ph_btys, ph_btye;
  // T
  std::vector<double> t_bxs, t_bxe, t_bys, t_bye;
  std::vector<double> t_btxs, t_btxe, t_btys, t_btye;
  // MU (2D)
  std::vector<double> mu_bxs, mu_bxe, mu_bys, mu_bye;
  std::vector<double> mu_btxs, mu_btxe, mu_btys, mu_btye;
  // QVAPOR
  std::vector<double> qvapor_bxs, qvapor_bxe, qvapor_bys, qvapor_bye;
  std::vector<double> qvapor_btxs, qvapor_btxe, qvapor_btys, qvapor_btye;

  // Optional moist species if present
  std::vector<double> qcloud_bxs, qcloud_bxe, qcloud_bys, qcloud_bye;
  std::vector<double> qcloud_btxs, qcloud_btxe, qcloud_btys, qcloud_btye;
  std::vector<double> qrain_bxs, qrain_bxe, qrain_bys, qrain_bye;
  std::vector<double> qrain_btxs, qrain_btxe, qrain_btys, qrain_btye;
  std::vector<double> qice_bxs, qice_bxe, qice_bys, qice_bye;
  std::vector<double> qice_btxs, qice_btxe, qice_btys, qice_btye;

  // Generic map for all variables present in file
  std::map<std::string, std::vector<double>> raw_fields;
};

// Reads wrfinput file at specified path.
WrfInputData read_wrfinput(const std::string& path);

// Reads wrfbdy file at specified path.
WrfBdyData read_wrfbdy(const std::string& path);

}  // namespace cppwrf::io_netcdf
