#pragma once

#include <vector>
#include <string>
#include <array>
#include <cstddef>

namespace cppwrf::init_real {

struct GridDims {
    int nx = 119;         // west_east (unstaggered)
    int ny = 89;          // south_north (unstaggered)
    int nz = 44;          // bottom_top (half levels)
    int nx_stag = 120;    // west_east_stag
    int ny_stag = 90;     // south_north_stag
    int nz_stag = 45;     // bottom_top_stag
    int bdy_width = 5;    // boundary zone width
    int num_metgrid_levels = 34; // input isobaric levels
};

// Plain flat vector container with 2D/3D index helpers
struct Array1D {
    std::vector<float> data;
    explicit Array1D(size_t n = 0, float val = 0.0f) : data(n, val) {}
    size_t size() const { return data.size(); }
    float& operator[](size_t i) { return data[i]; }
    float operator[](size_t i) const { return data[i]; }
    float* data_ptr() { return data.data(); }
    const float* data_ptr() const { return data.data(); }
};

struct Array2D {
    int d1 = 0, d2 = 0;
    std::vector<float> data;
    Array2D() = default;
    Array2D(int dim1, int dim2, float val = 0.0f) : d1(dim1), d2(dim2), data(dim1 * dim2, val) {}
    size_t size() const { return data.size(); }
    float& operator()(int i1, int i2) { return data[i1 * d2 + i2]; }
    float operator()(int i1, int i2) const { return data[i1 * d2 + i2]; }
    float* data_ptr() { return data.data(); }
    const float* data_ptr() const { return data.data(); }
};

struct Array3D {
    int d1 = 0, d2 = 0, d3 = 0;
    std::vector<float> data;
    Array3D() = default;
    Array3D(int dim1, int dim2, int dim3, float val = 0.0f) : d1(dim1), d2(dim2), d3(dim3), data(dim1 * dim2 * dim3, val) {}
    size_t size() const { return data.size(); }
    float& operator()(int i1, int i2, int i3) { return data[(i1 * d2 + i2) * d3 + i3]; }
    float operator()(int i1, int i2, int i3) const { return data[(i1 * d2 + i2) * d3 + i3]; }
    float* data_ptr() { return data.data(); }
    const float* data_ptr() const { return data.data(); }
};

// Initial state matching wrfinput_d01
struct WrfInputData {
    GridDims dims;

    // Base state constants
    float p00 = 100000.0f;
    float t00 = 290.0f;
    float tlp = 50.0f;
    float tiso = 200.0f;
    float p_top = 5000.0f;
    float p_strat = 0.0f;
    float tlp_strat = -11.0f;
    int hybrid_opt = 2;
    float etac = 0.2f;

    // 1D vertical coordinate arrays
    Array1D znw;  // [nz_stag] full eta levels
    Array1D znu;  // [nz] half eta levels
    Array1D dnw;  // [nz] full level eta spacing
    Array1D dn;   // [nz] half level eta spacing (dn[0] = 0)
    Array1D rdnw; // [nz] 1 / dnw
    Array1D rdn;  // [nz] 1 / dn (rdn[0] = 0)

    // Hybrid coordinate 1D coefficients
    Array1D c1f, c2f, c3f, c4f; // [nz_stag]
    Array1D c1h, c2h, c3h, c4h; // [nz]

    // 2D Mass / Grid fields [ny, nx]
    Array2D mub;       // [ny, nx] base state dry air mass in column (Pa)
    Array2D mu;        // [ny, nx] perturbation dry air mass (Pa)
    Array2D msft;      // [ny, nx] mapfac_m
    Array2D msfu;      // [ny, nx_stag] mapfac_u
    Array2D msfv;      // [ny_stag, nx] mapfac_v
    Array2D msftx;     // [ny, nx] mapfac_mx
    Array2D msfty;     // [ny, nx] mapfac_my
    Array2D msfux;     // [ny, nx_stag] mapfac_ux
    Array2D msfuy;     // [ny, nx_stag] mapfac_uy
    Array2D msfvx;     // [ny_stag, nx] mapfac_vx
    Array2D msfvy;     // [ny_stag, nx] mapfac_vy
    Array2D f;         // [ny, nx] Coriolis parameter
    Array2D e;         // [ny, nx] 2*Omega*cos(lat)
    Array2D sina;      // [ny, nx] sin of rotation angle
    Array2D cosa;      // [ny, nx] cos of rotation angle
    Array2D hgt;       // [ny, nx] terrain elevation (m)
    Array2D xlat;      // [ny, nx]
    Array2D xlong;     // [ny, nx]
    Array2D xlat_u;    // [ny, nx_stag]
    Array2D xlong_u;   // [ny, nx_stag]
    Array2D xlat_v;    // [ny_stag, nx]
    Array2D xlong_v;   // [ny_stag, nx]
    Array2D clat;      // [ny, nx]
    Array2D clong;     // [ny, nx]

    // 3D Prognostic and diagnostic state
    Array3D u;         // [nz, ny, nx_stag] u wind (m/s)
    Array3D v;         // [nz, ny_stag, nx] v wind (m/s)
    Array3D w;         // [nz_stag, ny, nx] w wind (m/s)
    Array3D ph;        // [nz_stag, ny, nx] perturbation geopotential (m^2/s^2)
    Array3D phb;       // [nz_stag, ny, nx] base geopotential (m^2/s^2)
    Array3D t;         // [nz, ny, nx] perturbation potential temperature (theta - 300) (K)
    Array3D t_init;    // [nz, ny, nx] base state potential temperature (K)
    Array3D p;         // [nz, ny, nx] perturbation pressure (Pa)
    Array3D pb;        // [nz, ny, nx] base pressure (Pa)
    Array3D al;        // [nz, ny, nx] perturbation inverse density (m^3/kg)
    Array3D alb;       // [nz, ny, nx] base inverse density (m^3/kg)
    Array3D qvapor;    // [nz, ny, nx] water vapor mixing ratio (kg/kg)
};

// Lateral boundary state for a single time interval
struct BoundaryInterval {
    // Current boundary values at start of interval (coupled with dry column mass mu)
    Array3D u_bxs, u_bxe; // [bdy_width, nz, ny]
    Array3D u_bys, u_bye; // [bdy_width, nz, nx_stag]
    Array3D v_bxs, v_bxe; // [bdy_width, nz, ny_stag]
    Array3D v_bys, v_bye; // [bdy_width, nz, nx]
    Array3D t_bxs, t_bxe; // [bdy_width, nz, ny]
    Array3D t_bys, t_bye; // [bdy_width, nz, nx]
    Array3D ph_bxs, ph_bxe; // [bdy_width, nz_stag, ny]
    Array3D ph_bys, ph_bye; // [bdy_width, nz_stag, nx]
    Array3D q_bxs, q_bxe; // [bdy_width, nz, ny]
    Array3D q_bys, q_bye; // [bdy_width, nz, nx]
    Array3D mu_bxs, mu_bxe; // [bdy_width, 1, ny]
    Array3D mu_bys, mu_bye; // [bdy_width, 1, nx]

    // Boundary tendencies over interval (per second: / dt)
    Array3D u_btxs, u_btxe;
    Array3D u_btys, u_btye;
    Array3D v_btxs, v_btxe;
    Array3D v_btys, v_btye;
    Array3D t_btxs, t_btxe;
    Array3D t_btys, t_btye;
    Array3D ph_btxs, ph_btxe;
    Array3D ph_btys, ph_btye;
    Array3D q_btxs, q_btxe;
    Array3D q_btys, q_btye;
    Array3D mu_btxs, mu_btxe;
    Array3D mu_btys, mu_btye;

    std::string this_bdy_time;
    std::string next_bdy_time;
};

struct WrfBdyData {
    GridDims dims;
    float dt_interval = 21600.0f; // 6 hours
    std::vector<BoundaryInterval> intervals; // 4 intervals across 5 input dates
};

} // namespace cppwrf::init_real
