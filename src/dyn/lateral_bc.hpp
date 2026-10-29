#pragma once
// Specified lateral boundary conditions of a real-data run (WRF specified = .true.,
// spec_bdy_width = 5, spec_zone = 1, relax_zone = 4), device-resident and single source.
// The port of share/module_bc.F and dyn_em/module_bc_em.F as parallel kernels over the boundary
// ring; tests/dyn/lateral_bc compares every call with WRF's own routines (specified and nested).
//
// Call order inside a time step (WRF solve_em), done by the dycore when a LateralBC is attached:
//   start of step: advance_time(dt)                       (dtbc = dtbc + dt)
//   each RK stage, after the tendencies: rk 1: relax_dry  (into *_tendf / mu_tend, the WRF
//     *_save -> rk_addtend_dry path), every stage: spec_dry (spec-zone tendencies)
//   small steps: update_uv after advance_uv, update_mu_t after advance_mu_t,
//     update_ph_w after advance_w (spec_bdyupdate_ph, zero_grad_bdy on w)
//   scalars: relax_scalar (rk 1) and spec_scalar into the held scalar tendency
//   end of step: final (spec_bdy_final for u, v, t, ph, mu; and w in nests)
// Nests (nested = true): W also relaxed/specified/updated/finalised, all scalars with data too.
// Boundary data layout per side and variable: (pos, k, b), pos fastest; xs/xe use pos = j,
// ys/ye use pos = i (staggered extents for U in x and V in y), b = 0 is the outermost row.
#include <array>
#include <vector>

#include "core/field.hpp"
#include "dyn/state.hpp"

namespace cppwrf::dyn {

class LateralBC {
 public:
  enum Var { U = 0, V, T, PH, MU, W, NVAR };  // W: nests only
  enum Side { XS = 0, XE, YS, YE };

  // nested = WRF config_flags%nested: W is relaxed, specified and updated from boundary data
  // (instead of zero_grad_bdy) and every scalar with data gets relax/spec (no flow_dep).
  LateralBC(const Grid& g, double dt, int width = 5, int spec_zone = 1, int relax_zone = 4, bool nested = false);
  int spec_zone() const { return spec_zone_; }
  // the relaxation-zone weights fcx, gcx for the large step dt
  void set_dt(double dt);

  // Host data of one variable and side for the current boundary interval: values and
  // tendencies, (pos, k, b) with pos fastest.
  void set(Var v, Side s, const std::vector<double>& values, const std::vector<double>& tends);
  // The same for scalar n (State order: moist species first), e.g. QVAPOR of wrfbdy.
  void set_scalar(int n, Side s, const std::vector<double>& values, const std::vector<double>& tends);
  // Device-side access for nest forcing (src/nest): the slab of variable v / scalar n, side s,
  // values (tend = false) or tendencies, layout (pos, k, b) with pos fastest; npos = slab_npos,
  // nk = slab_nk. mark_set() declares the data valid (has() / has_scalar() true).
  real* slab(Var v, Side s, bool tend) { return data_[v].get() + off_[v][(tend ? 4 : 0) + s]; }
  real* scalar_slab(int n, Side s, bool tend);  // allocates scalar n on first use
  int slab_npos(Var v, Side s) const { return npos(v, s); }
  int slab_nk(Var v) const { return nk_[v]; }
  void mark_set(Var v) { have_[v] = true; }
  bool has(Var v) const { return have_[v]; }
  bool has_scalar(int n) const { return n >= 0 && n < int(sdata_.size()) && !sdata_[n].off.empty(); }
  int width() const { return width_; }

  double dtbc = 0;  // seconds since the start of the current boundary interval
  void advance_time(double dt) { dtbc += dt; }

  void relax_dry(const State& s, Tendencies& t, const BaseState& b) const;
  void spec_dry(Tendencies& t) const;
  void update_uv(State& s, const Tendencies& t, real dts) const;
  void update_mu_t(State& s, const Tendencies& t, real dts) const;
  void update_ph_w(State& s, const Field3& ph_save, const Tendencies& t, const BaseState& b, real dts) const;
  void relax_scalar(int n, Field3& q_tend, const Field3& q, const State& s, const BaseState& b) const;
  void spec_scalar(int n, Field3& q_tend) const;
  // flow_dep_bdy for a moist species without boundary data (have_bcs_moist = .false.): spec zone
  // = inner neighbour at outflow, 0 at inflow, by the sign of the averaged fluxes ru_m, rv_m.
  void flow_dep(Field3& q, const State& s) const;
  void final(State& s, const BaseState& b) const;

 private:
  struct ScalarData {
    Buffer data;
    std::vector<std::size_t> off;  // 8 slabs; empty = no data for this scalar
  };
  Grid g_;
  int width_, spec_zone_, relax_zone_;
  bool nested_;
  std::vector<double> fcx_h_, gcx_h_;
  Buffer coef_;  // fcx[width], gcx[width] on the device
  std::array<bool, NVAR> have_{};
  // per variable: 8 slabs (values xs, xe, ys, ye, then tendencies xs, xe, ys, ye)
  std::array<Buffer, NVAR> data_;
  std::array<std::array<std::size_t, 8>, NVAR> off_{};
  std::array<int, NVAR> nk_{};
  std::vector<ScalarData> sdata_;
  int npos(Var v, Side s) const;
};

}  // namespace cppwrf::dyn
