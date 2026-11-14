// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once

#include <cstddef>
#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace cppwrf::io_netcdf {

// The part of the horizontal grid a reader keeps of the 3-D variables: the mass points [i0, i1) x [j0, j1) (and
// the staggered point after the last); whole(): everything. A decomposed run (docs/MPI.md) reads only the
// worker's patch and its halo, so the start does not hold the domain's 3-D arrays on every worker.
struct Window {
  int i0 = 0, i1 = -1, j0 = 0, j1 = -1;
  bool whole() const { return i1 < 0; }
};

// A (nk, ny, nx) variable of a WRF file (a 2-D one: nk = 1), whole or only a Window of it, indexed with the whole
// array's flat index (k ny + j) nx + i as a std::vector of the whole variable is; an index outside the window is an
// error. size() is the whole variable's.
class WArray {
 public:
  WArray() = default;
  // whole: n values (the flat vector)
  explicit WArray(std::size_t n) : v_(n), n_(n) {}
  // (nk, ny, nx) with the columns [i0, i1) x [j0, j1)
  WArray(std::size_t nx, std::size_t ny, std::size_t nk, std::size_t i0, std::size_t i1, std::size_t j0, std::size_t j1)
      : v_((i1 - i0) * (j1 - j0) * nk), n_(nx * ny * nk), nx_(nx), ny_(ny), i0_(i0), i1_(i1), j0_(j0), j1_(j1),
        win_(!(i0 == 0 && i1 == nx && j0 == 0 && j1 == ny)) {}
  // (nk, ny, nx) with the window w of the mass grid: the columns [w.i0, min(nx, w.i1 + 1)) x [w.j0, min(ny, w.j1 + 1))
  static WArray make(std::size_t nx, std::size_t ny, std::size_t nk, const Window& w) {
    if (w.whole()) return WArray(nx, ny, nk, 0, nx, 0, ny);
    const std::size_t i1 = std::size_t(w.i1) + 1 < nx ? std::size_t(w.i1) + 1 : nx;
    const std::size_t j1 = std::size_t(w.j1) + 1 < ny ? std::size_t(w.j1) + 1 : ny;
    return WArray(nx, ny, nk, std::size_t(w.i0), i1, std::size_t(w.j0), j1);
  }
  std::size_t size() const { return n_; }
  bool empty() const { return n_ == 0; }
  void clear() { *this = WArray(); }
  bool windowed() const { return win_; }
  // the columns held: i0, i1, j0, j1
  std::size_t i0() const { return i0_; }
  std::size_t i1() const { return i1_; }
  std::size_t j0() const { return j0_; }
  std::size_t j1() const { return j1_; }
  double operator[](std::size_t n) const { return v_[at(n)]; }
  double& operator[](std::size_t n) { return v_[at(n)]; }
  // the values held, (nk, j1 - j0, i1 - i0) (the reader's buffer)
  std::vector<double>& held() { return v_; }
  const std::vector<double>& held() const { return v_; }

 private:
  std::size_t at(std::size_t n) const {
    if (!win_) return n;
    const std::size_t nxy = nx_ * ny_, k = n / nxy, r = n - k * nxy, j = r / nx_, i = r - j * nx_;
    if (i < i0_ || i >= i1_ || j < j0_ || j >= j1_ || n >= n_)
      throw std::out_of_range("WArray: point (" + std::to_string(i) + ", " + std::to_string(j) + ") is outside the window [" +
                              std::to_string(i0_) + ", " + std::to_string(i1_) + ") x [" + std::to_string(j0_) + ", " +
                              std::to_string(j1_) + ")");
    return (k * (j1_ - j0_) + (j - j0_)) * (i1_ - i0_) + (i - i0_);
  }
  std::vector<double> v_;
  std::size_t n_ = 0, nx_ = 0, ny_ = 0, i0_ = 0, i1_ = 0, j0_ = 0, j1_ = 0;
  bool win_ = false;
};
// The window of a file with nx x ny mass points (the reader calls it once the dimensions are known)
using WindowFn = std::function<Window(int nx, int ny)>;
// A (Time = 1) variable with horizontal dimensions as the last two, the window w of it; empty if absent
WArray read_window_var(const std::string& path, const char* name, const Window& w);

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
  WArray pb;       // Base pressure [Pa]
  WArray p;        // Perturbation pressure [Pa]
  WArray t;        // Perturbation potential temperature (theta - 300 K)
  WArray t_init;   // Initial potential temperature
  WArray alb;      // Base inverse density
  WArray al;       // Perturbation inverse density
  WArray qvapor;   // Water vapor mixing ratio [kg/kg]

  // 3D Staggered fields
  WArray u;        // [nz * ny * nx_stag]
  WArray v;        // [nz * ny_stag * nx]
  WArray w;        // [nz_stag * ny * nx]
  WArray phb;      // [nz_stag * ny * nx]
  WArray ph;       // [nz_stag * ny * nx]

  // Other moist species if present
  WArray qcloud;
  WArray qrain;
  WArray qice;
  WArray qsnow;
  WArray qgraup;
  std::map<std::string, WArray> moist_species;

  // Surface and soil fields
  int num_soil_layers = 0;
  std::vector<double> zs;
  std::vector<double> dzs;
  std::vector<double> tsk;
  std::vector<double> tslb;
  std::vector<double> smois;
  std::vector<double> sh2o;
  std::vector<double> snow;
  std::vector<double> snowh;
  std::vector<double> snowc;
  std::vector<double> canwat;
  std::vector<double> albedo;
  std::vector<double> albbck;
  std::vector<double> emiss;
  std::vector<double> embck;
  std::vector<double> snoalb;
  std::vector<double> tmn;
  std::vector<double> xland;
  std::vector<double> seaice;
  std::vector<double> lakemask;
  std::vector<double> water_depth;
  std::vector<double> vegfra;
  std::vector<double> shdmin;
  std::vector<double> shdmax;
  std::vector<double> lai;
  std::vector<double> lu_index;
  std::vector<double> ivgtyp;
  std::vector<double> isltyp;
  std::vector<double> u10;
  std::vector<double> v10;
  std::vector<double> th2;
  std::vector<double> t2;
  std::vector<double> q2;
  std::vector<double> ust;
  std::vector<double> pblh;
  std::vector<double> hfx;
  std::vector<double> qfx;
  std::vector<double> lh;
  std::vector<double> znt;
  std::vector<double> qsfc;

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
// window: the part of the 3-D variables to keep (none: everything)
WrfInputData read_wrfinput(const std::string& path, const WindowFn& window = {});

// Reads wrfbdy file at specified path.
WrfBdyData read_wrfbdy(const std::string& path);

}  // namespace cppwrf::io_netcdf
