#pragma once

// ============================================================================
// CPPWRF Nest Initialisation Kernels (blend_terrain, adjust_tempqv, copy_3d_field)
// Ports WRF 4.8.0 dyn_em/nest_init_utils.F
// Implements nest terrain blending and thermodynamic lapse-rate adjustments
// for input_from_file initialisation.
// ============================================================================

#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdint>
#include "core/backend.hpp"
#include "core/field.hpp"
#include "nest/nest_types.hpp"

namespace cppwrf::nest {

// ============================================================================
// blend_terrain: Smoothly transitions terrain height, mub, and phb between
// interpolated parent values along nest boundary and fine input in the interior.
//
// Boundary zone (1..spec_bdy_width): 100% interpolated parent terrain
// Transition zone (spec_bdy_width+1 .. spec_bdy_width+blend_width): linear blend
// Interior: 100% fine input terrain
// ============================================================================
inline void blend_terrain(
    const double* ter_interpolated,
    double* ter_input,
    int ids, int ide, int jds, int jde, int kds, int kde,
    int ims, int ime, int jms, int jme, int kms, int kme,
    int ips, int ipe, int jps, int jpe, int kps, int kpe,
    int spec_bdy_width = 5,
    int blend_width = 5) {

  const int ni = ime - ims + 1;
  const int nk = kme - kms + 1;
  auto idx = [&](int i, int k, int j) -> std::size_t {
    return static_cast<std::size_t>((i - ims) + (k - kms) * ni + (j - jms) * ni * nk);
  };

  const double r_blend_zones = 1.0 / static_cast<double>(blend_width + 1);

  std::vector<double> ter_temp(static_cast<std::size_t>(ni * nk * (jme - jms + 1)));

  for (int j = jps; j <= std::min(jpe, jde - 1); ++j) {
    for (int k = kps; k <= kpe; ++k) {
      for (int i = ips; i <= std::min(ipe, ide - 1); ++i) {
        ter_temp[idx(i, k, j)] = ter_input[idx(i, k, j)];
      }
    }
  }

  for (int j = jps; j <= std::min(jpe, jde - 1); ++j) {
    for (int k = kps; k <= kpe; ++k) {
      for (int i = ips; i <= std::min(ipe, ide - 1); ++i) {
        for (int blend_cell = blend_width; blend_cell >= 1; --blend_cell) {
          if ((i == spec_bdy_width + blend_cell) || (j == spec_bdy_width + blend_cell) ||
              (i == ide - spec_bdy_width - blend_cell) || (j == jde - spec_bdy_width - blend_cell)) {
            ter_temp[idx(i, k, j)] = (static_cast<double>(blend_cell) * ter_input[idx(i, k, j)] +
                                      static_cast<double>(blend_width + 1 - blend_cell) * ter_interpolated[idx(i, k, j)]) *
                                     r_blend_zones;
          }
        }
        if ((i <= spec_bdy_width) || (j <= spec_bdy_width) ||
            (i >= ide - spec_bdy_width) || (j >= jde - spec_bdy_width)) {
          ter_temp[idx(i, k, j)] = ter_interpolated[idx(i, k, j)];
        }
      }
    }
  }

  for (int j = jps; j <= std::min(jpe, jde - 1); ++j) {
    for (int k = kps; k <= kpe; ++k) {
      for (int i = ips; i <= std::min(ipe, ide - 1); ++i) {
        ter_input[idx(i, k, j)] = ter_temp[idx(i, k, j)];
      }
    }
  }
}

// ============================================================================
// copy_3d_field: Copies 3D field from source to destination
// ============================================================================
inline void copy_3d_field(
    double* dest,
    const double* src,
    int ids, int ide, int jds, int jde, int kds, int kde,
    int ims, int ime, int jms, int jme, int kms, int kme,
    int ips, int ipe, int jps, int jpe, int kps, int kpe) {

  const int ni = ime - ims + 1;
  const int nk = kme - kms + 1;
  auto idx = [&](int i, int k, int j) -> std::size_t {
    return static_cast<std::size_t>((i - ims) + (k - kms) * ni + (j - jms) * ni * nk);
  };

  for (int j = jps; j <= std::min(jpe, jde - 1); ++j) {
    for (int k = kps; k <= kpe; ++k) {
      for (int i = ips; i <= std::min(ipe, ide - 1); ++i) {
        dest[idx(i, k, j)] = src[idx(i, k, j)];
      }
    }
  }
}

// ============================================================================
// adjust_tempqv: Adjusts potential temperature (t_2) and water vapor (qv)
// following terrain blending so that relative humidity (RH) is conserved
// across the hydrostatic pressure change dp = p_new - p_old.
// ============================================================================
inline void adjust_tempqv(
    const double* mub,
    const double* save_mub,
    const double* c3h,
    const double* c4h,
    const double* /*znw*/,
    double p_top,
    double* th,
    const double* pp,
    double* qv,
    int use_theta_m,
    int ids, int ide, int jds, int jde, int kds, int kde,
    int ims, int ime, int jms, int jme, int kms, int kme,
    int ips, int ipe, int jps, int jpe, int kps, int kpe) {

  constexpr double R_d = 287.0;
  constexpr double R_v = 461.6;
  constexpr double R_ratio = R_v / R_d;
  constexpr double p00 = 1.0e5;
  constexpr double rcp = 2.0 / 7.0;
  constexpr double gamma_coeff = -191.86e-3;

  const int ni = ime - ims + 1;
  const int nk = kme - kms + 1;

  auto idx2d = [&](int i, int j) -> std::size_t {
    return static_cast<std::size_t>((i - ims) + (j - jms) * ni);
  };
  auto idx3d = [&](int i, int k, int j) -> std::size_t {
    return static_cast<std::size_t>((i - ims) + (k - kms) * ni + (j - jms) * ni * nk);
  };

  const int j_end = std::min(jpe, jde - 1);
  const int i_end = std::min(ipe, ide - 1);

  std::vector<double> rh(static_cast<std::size_t>(ni * nk * (jme - jms + 1)));
  std::vector<double> p_old_arr(static_cast<std::size_t>(ni * nk * (jme - jms + 1)));

  // Pass 1: compute p_old, temperature, and initial RH (conserved during adjustment)
  for (int j = jps; j <= j_end; ++j) {
    for (int k = kps; k <= kpe - 1; ++k) {
      for (int i = ips; i <= i_end; ++i) {
        std::size_t p3 = idx3d(i, k, j);
        std::size_t p2 = idx2d(i, j);

        double p_old = c4h[k - kms] + c3h[k - kms] * save_mub[p2] + p_top + pp[p3];
        p_old_arr[p3] = p_old;

        double tc;
        if (use_theta_m == 1) {
          tc = (th[p3] + 300.0) * std::pow(p_old / p00, rcp) / (1.0 + R_ratio * qv[p3]) - 273.15;
        } else {
          tc = (th[p3] + 300.0) * std::pow(p_old / p00, rcp) - 273.15;
        }

        double es = 610.78 * std::exp(17.0809 * tc / (234.175 + tc));
        double e = qv[p3] * p_old / (0.622 + qv[p3]);
        rh[p3] = e / es;
      }
    }
  }

  // Pass 2: compute p_new, temperature correction dth, and new qv
  for (int j = jps; j <= j_end; ++j) {
    for (int k = kps; k <= kpe - 1; ++k) {
      for (int i = ips; i <= i_end; ++i) {
        std::size_t p3 = idx3d(i, k, j);
        std::size_t p2 = idx2d(i, j);

        double p_old = p_old_arr[p3];
        double p_new = c4h[k - kms] + c3h[k - kms] * mub[p2] + p_top + pp[p3];

        double thloc;
        if (use_theta_m == 1) {
          thloc = (th[p3] + 300.0) / (1.0 + R_ratio * qv[p3]);
        } else {
          thloc = (th[p3] + 300.0);
        }

        double dp = p_new - p_old;
        double p_sum = p_new + p_old;
        double dth1 = gamma_coeff * (thloc) / p_sum * dp;
        double dth  = gamma_coeff * (thloc + 0.5 * dth1) / p_sum * dp;

        if (use_theta_m == 1) {
          th[p3] = (thloc + dth) * (1.0 + R_ratio * qv[p3]) - 300.0;
        } else {
          th[p3] = (thloc + dth) - 300.0;
        }

        double tc = (thloc + dth) * std::pow(p_new / p00, rcp) - 273.15;
        double es = 610.78 * std::exp(17.0809 * tc / (234.175 + tc));
        double e = rh[p3] * es;
        qv[p3] = 0.622 * e / (p_new - e);
      }
    }
  }
}

// ============================================================================
// Single-Source View-based blend_terrain for 2D/3D fields (0-based i0, j0, k)
// ============================================================================
struct KBlendTerrain2D {
  View2<const real> interp_v;
  View2<real> in_v;
  int ide, jde, spec_bdy_width, blend_width;
  double r_blend_zones;

  CPPWRF_HD void operator()(int i0, int j0, int /*k*/) const {
    const int i = i0 + 1;
    const int j = j0 + 1;
    const double v_in = static_cast<double>(in_v(i0, j0));
    const double v_int = static_cast<double>(interp_v(i0, j0));
    double v_out = v_in;
    for (int blend_cell = blend_width; blend_cell >= 1; --blend_cell) {
      if ((i == spec_bdy_width + blend_cell) || (j == spec_bdy_width + blend_cell) ||
          (i == ide - spec_bdy_width - blend_cell) || (j == jde - spec_bdy_width - blend_cell)) {
        v_out = (static_cast<double>(blend_cell) * v_in +
                 static_cast<double>(blend_width + 1 - blend_cell) * v_int) * r_blend_zones;
      }
    }
    if ((i <= spec_bdy_width) || (j <= spec_bdy_width) ||
        (i >= ide - spec_bdy_width) || (j >= jde - spec_bdy_width)) {
      v_out = v_int;
    }
    in_v(i0, j0) = static_cast<real>(v_out);
  }
};

template <class Field2T>
inline void blend_terrain_field2(
    const Field2T& ter_interpolated,
    Field2T& ter_input,
    int nx, int ny,
    int spec_bdy_width = 5,
    int blend_width = 5) {

  const int ide = nx + 1;
  const int jde = ny + 1;
  const double r_blend_zones = 1.0 / static_cast<double>(blend_width + 1);
  View2<const real> interp_v{ter_interpolated.data(), int(ter_interpolated.layout().sj)};
  View2<real> in_v{ter_input.data(), int(ter_input.layout().sj)};

  parallel_for(Range3D{0, nx, 0, ny, 0, 1},
               KBlendTerrain2D{interp_v, in_v, ide, jde, spec_bdy_width, blend_width, r_blend_zones});
}

struct KBlendTerrain3D {
  View3<const real> interp_v;
  View3<real> in_v;
  int ide, jde, spec_bdy_width, blend_width;
  double r_blend_zones;

  CPPWRF_HD void operator()(int i0, int j0, int k) const {
    const int i = i0 + 1;
    const int j = j0 + 1;
    const double v_in = static_cast<double>(in_v(i0, j0, k));
    const double v_int = static_cast<double>(interp_v(i0, j0, k));
    double v_out = v_in;
    for (int blend_cell = blend_width; blend_cell >= 1; --blend_cell) {
      if ((i == spec_bdy_width + blend_cell) || (j == spec_bdy_width + blend_cell) ||
          (i == ide - spec_bdy_width - blend_cell) || (j == jde - spec_bdy_width - blend_cell)) {
        v_out = (static_cast<double>(blend_cell) * v_in +
                 static_cast<double>(blend_width + 1 - blend_cell) * v_int) * r_blend_zones;
      }
    }
    if ((i <= spec_bdy_width) || (j <= spec_bdy_width) ||
        (i >= ide - spec_bdy_width) || (j >= jde - spec_bdy_width)) {
      v_out = v_int;
    }
    in_v(i0, j0, k) = static_cast<real>(v_out);
  }
};

template <class Field3T>
inline void blend_terrain_field3(
    const Field3T& ter_interpolated,
    Field3T& ter_input,
    int nx, int ny, int nk,
    int spec_bdy_width = 5,
    int blend_width = 5) {

  const int ide = nx + 1;
  const int jde = ny + 1;
  const double r_blend_zones = 1.0 / static_cast<double>(blend_width + 1);
  View3<const real> interp_v{ter_interpolated.data(), int(ter_interpolated.layout().sj), int(ter_interpolated.layout().sk)};
  View3<real> in_v{ter_input.data(), int(ter_input.layout().sj), int(ter_input.layout().sk)};

  parallel_for(Range3D{0, nx, 0, ny, 0, nk},
               KBlendTerrain3D{interp_v, in_v, ide, jde, spec_bdy_width, blend_width, r_blend_zones});
}

// ============================================================================
// Single-Source View-based adjust_tempqv for CPPWRF fields
// ============================================================================
struct KAdjustTempQv {
  View2<const real> mub_v;
  View2<const real> smub_v;
  View1<const real> c3h_v;
  View1<const real> c4h_v;
  double p_top;
  View3<real> th_v;
  View3<const real> pp_v;
  View3<real> qv_v;
  int use_theta_m;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    constexpr double R_d = 287.0;
    constexpr double R_v = 461.6;
    constexpr double R_ratio = R_v / R_d;
    constexpr double p00 = 1.0e5;
    constexpr double rcp = 2.0 / 7.0;
    constexpr double gamma_coeff = -191.86e-3;

    const double p_old = static_cast<double>(c4h_v(k)) +
                         static_cast<double>(c3h_v(k)) * static_cast<double>(smub_v(i, j)) +
                         p_top + static_cast<double>(pp_v(i, j, k));
    const double th_in = static_cast<double>(th_v(i, j, k));
    const double qv_in = static_cast<double>(qv_v(i, j, k));
    const double moist_fac = (use_theta_m == 1) ? (1.0 + R_ratio * qv_in) : 1.0;

    const double tc_old = (th_in + 300.0) * rpow(p_old / p00, rcp) / moist_fac - 273.15;
    const double es_old = 610.78 * std::exp(17.0809 * tc_old / (234.175 + tc_old));
    const double e_old = qv_in * p_old / (0.622 + qv_in);
    const double rh = e_old / es_old;

    const double p_new = static_cast<double>(c4h_v(k)) +
                         static_cast<double>(c3h_v(k)) * static_cast<double>(mub_v(i, j)) +
                         p_top + static_cast<double>(pp_v(i, j, k));
    const double thloc = (th_in + 300.0) / moist_fac;
    const double dp = p_new - p_old;
    const double p_sum = p_new + p_old;
    const double dth1 = gamma_coeff * thloc / p_sum * dp;
    const double dth = gamma_coeff * (thloc + 0.5 * dth1) / p_sum * dp;

    th_v(i, j, k) = static_cast<real>((thloc + dth) * moist_fac - 300.0);

    const double tc_new = (thloc + dth) * rpow(p_new / p00, rcp) - 273.15;
    const double es_new = 610.78 * std::exp(17.0809 * tc_new / (234.175 + tc_new));
    const double e_new = rh * es_new;
    qv_v(i, j, k) = static_cast<real>(0.622 * e_new / (p_new - e_new));
  }
};

template <class BaseField2T, class Array1T, class Field3T>
inline void adjust_tempqv_fields(
    const BaseField2T& mub,
    const BaseField2T& save_mub,
    const Array1T& c3h,
    const Array1T& c4h,
    double p_top,
    Field3T& th,
    const Field3T& pp,
    Field3T& qv,
    int use_theta_m,
    int nx, int ny, int nz) {

  View2<const real> mub_v{mub.data(), int(mub.layout().sj)};
  View2<const real> smub_v{save_mub.data(), int(save_mub.layout().sj)};
  View1<const real> c3h_v = c3h.view();
  View1<const real> c4h_v = c4h.view();
  View3<real> th_v = th.view();
  View3<const real> pp_v = pp.cview();
  View3<real> qv_v = qv.view();

  parallel_for(Range3D{0, nx, 0, ny, 0, nz},
               KAdjustTempQv{mub_v, smub_v, c3h_v, c4h_v, p_top, th_v, pp_v, qv_v, use_theta_m});
}

}  // namespace cppwrf::nest
