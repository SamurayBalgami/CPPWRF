#pragma once
// Specified lateral boundary conditions of a real-data run (WRF specified = .true.,
// spec_bdy_width = 5, spec_zone = 1, relax_zone = 4), device-resident and single source.
// The point formulas and index ranges are those of src/bdy/bdy_kernels.hpp (port of
// share/module_bc.F and dyn_em/module_bc_em.F, differential-tested against the Fortran); here
// they are restructured as parallel kernels over the boundary ring.
//
// Call order inside a time step (WRF solve_em), done by the dycore when a LateralBC is attached:
//   start of step: advance_time(dt)                       (dtbc = dtbc + dt)
//   each RK stage, after the tendencies: rk 1: relax_dry  (into *_tendf / mu_tend, the WRF
//     *_save -> rk_addtend_dry path), every stage: spec_dry (spec-zone tendencies)
//   small steps: update_uv after advance_uv, update_mu_t after advance_mu_t,
//     update_ph_w after advance_w (spec_bdyupdate_ph, zero_grad_bdy on w)
//   scalars: relax_scalar (rk 1) and spec_scalar into the held scalar tendency
//   end of step: final (spec_bdy_final for u, v, t, ph, mu)
// Boundary data layout per side and variable: (pos, k, b), pos fastest; xs/xe use pos = j,
// ys/ye use pos = i (staggered extents for U in x and V in y), b = 0 is the outermost row.
#include <array>
#include <vector>

#include "core/field.hpp"
#include "dyn/state.hpp"

namespace cppwrf::dyn {

class LateralBC {
 public:
  enum Var { U = 0, V, T, PH, MU, QV, NVAR };
  enum Side { XS = 0, XE, YS, YE };

  LateralBC(const Grid& g, double dt, int width = 5, int spec_zone = 1, int relax_zone = 4);

  // Host data of one variable and side for the current boundary interval: values and
  // tendencies, (pos, k, b) with pos fastest.
  void set(Var v, Side s, const std::vector<double>& values, const std::vector<double>& tends);
  bool has(Var v) const { return have_[v]; }
  int width() const { return width_; }

  double dtbc = 0;  // seconds since the start of the current boundary interval
  void advance_time(double dt) { dtbc += dt; }

  void relax_dry(const State& s, Tendencies& t, const BaseState& b) const;
  void spec_dry(Tendencies& t) const;
  void update_uv(State& s, const Tendencies& t, real dts) const;
  void update_mu_t(State& s, const Tendencies& t, real dts) const;
  void update_ph_w(State& s, const Field3& ph_save, const Tendencies& t, const BaseState& b, real dts) const;
  void relax_scalar(Field3& q_tend, const Field3& q, const State& s, const BaseState& b) const;
  void spec_scalar(Field3& q_tend) const;
  void final(State& s, const BaseState& b) const;

 private:
  Grid g_;
  int width_, spec_zone_, relax_zone_;
  std::vector<double> fcx_h_, gcx_h_;
  Buffer coef_;  // fcx[width], gcx[width] on the device
  std::array<bool, NVAR> have_{};
  // per variable: 8 slabs (values xs, xe, ys, ye, then tendencies xs, xe, ys, ye)
  std::array<Buffer, NVAR> data_;
  std::array<std::array<std::size_t, 8>, NVAR> off_{};
  std::array<int, NVAR> nk_{};
  int npos(Var v, Side s) const;
};

}  // namespace cppwrf::dyn
