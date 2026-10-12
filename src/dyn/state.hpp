#pragma once
// Core dynamical-core containers. Member names follow the WRF Registry so that the TN and
// the WRF source map 1:1 onto this code.
//
// Index conventions (C++ k=0 is the lowest level; WRF Fortran k = C++ k + 1):
//   mass levels k = 0..nz-1, w levels k = 0..nz (k=0 surface, k=nz model top).
//   znw(k)  eta at w level k (znw(0)=1, znw(nz)=0); znu(k) = (znw(k)+znw(k+1))/2.
//   dnw(k)  = znw(k+1)-znw(k) (<0), k=0..nz-1;  rdnw = 1/dnw.
//   dn(k)   = znu(k)-znu(k-1) = (dnw(k)+dnw(k-1))/2 (<0), k=1..nz-1; rdn = 1/dn; dn(0)=dn(nz)=0.
//   fnm(k), fnp(k), k=1..nz-1: f_w(k) = fnm(k) f(k) + fnp(k) f(k-1)   (linear in eta).
//   c1h,c2h (mass levels), c1f,c2f (w levels): column mass at a level = c1*mu + c2 (TN 2.1).
//   c3h,c4h,c3f,c4f: dry hydrostatic pressure p_d = c3*mu_d + c4 + p_top (TN eq. 2.2).
//   hybrid_opt = 0: c1 = 1, c2 = 0, c3 = eta, c4 = 0.
//
// Prognostic convention as in WRF: *_2 hold the uncoupled stage values
// (u, v, w [m/s], t = theta - t0 [K], ph = phi' [m2/s2], mu = mu_d' [Pa]);
// *_1 hold the time-t values (the RK3 driver copies *_2 -> *_1 at the start of each step).
// Tendencies are mass-coupled, as in WRF.
// Map factors: idealized cases use m = 1 everywhere. When map factors are added, they go
// into BaseState (msfux, msfuy, msfvx, msfvy, msftx, msfty) -- nowhere else.
#include <string>
#include <vector>

#include "core/field.hpp"
#include "core/grid.hpp"

namespace cppwrf::dyn {

struct State {
  explicit State(const Grid& g);
  // time-t (_1) and stage (_2) prognostics
  Field3 u_1, u_2;    // U
  Field3 v_1, v_2;    // V
  Field3 w_1, w_2;    // W
  Field3 t_1, t_2;    // M, theta - t0
  Field3 ph_1, ph_2;  // W, perturbation geopotential
  Field2 mu_1, mu_2;  // M, perturbation dry column mass
  // diagnostics
  Field3 p;    // M, perturbation pressure
  Field3 al;   // M, perturbation inverse density
  Field3 alt;  // M, full inverse density (al + alb)
  Field3 php;  // M, full geopotential at mass levels
  Field3 ww;   // W, mass-coupled eta velocity (mu * deta/dt)
  // column masses (full = mub + mu): at t (mut, muu, muv) and stage (muts, muus, muvs)
  Field2 mut, muu, muv;     // M, U, V
  Field2 muts, muus, muvs;  // M, U, V
  // mass-coupled velocities of the stage state (rk_step_prep: couple_momentum, calc_ww_cp)
  Field3 ru, rv, rw;  // U, V, W
  // time-averaged mass fluxes from the acoustic steps (for scalar advection)
  Field3 ru_m, rv_m, ww_m;  // U, V, W

  // Mass-point scalars, uncoupled like the WRF moist/scalar arrays: entries 0..n_moist-1
  // are moist species [kg/kg] (WRF moist: qv, qc, qr, ...), the rest passive tracers (WRF
  // tracer). q_1 = time t (WRF moist_old), q_2 = stage value (WRF moist).
  std::vector<std::string> q_names;
  std::vector<Field3> q_1, q_2;
  int n_moist = 0;
  int P_QV = -1, P_QC = -1, P_QR = -1, P_QI = -1, P_QS = -1, P_QG = -1;  // indices into q_* (-1: absent)
  Field3 qtot;  // M, sum of the moist species of the stage state (calc_cq input); allocated if n_moist > 0
  Field3 h_diabatic;  // M, microphysics heating of the last step [theta_m K/s] (WRF h_diabatic); with qtot
  // Adds a scalar (moist species before passive tracers). Returns its index.
  int add_scalar(const Grid& g, const std::string& name, bool moist);
  int n_scalars() const { return int(q_2.size()); }
};

struct Tendencies {
  explicit Tendencies(const Grid& g);
  // large-step tendencies (all RK stages), mass-coupled
  Field3 ru_tend, rv_tend, rw_tend, t_tend, ph_tend;  // U V W M W
  Field2 mu_tend;                                     // M
  // tendencies computed in RK stage 1 only and held fixed (diffusion/physics, WRF *_tendf)
  Field3 ru_tendf, rv_tendf, rw_tendf, t_tendf, ph_tendf;
  // scalar tendencies computed in stage 1 (diffusion, physics), mass-coupled (WRF moist_tend)
  std::vector<Field3> q_tendf;
  void add_scalar(const Grid& g, const std::string& name);
};

struct BaseState {
  BaseState(const Grid& g, bool flat);
  bool flat = true;  // horizontally uniform base state -> 1D storage
  BaseField2 mub;                // M, base dry column mass (p_surf_base - p_top)
  BaseField3 pb, alb, t_init;    // M: base pressure, base inverse density, initial theta - t0
  BaseField3 qv_base;            // M: base-state vapour profile (WRF qv_base; 0 if dry)
  BaseField3 phb;                // W: base geopotential (incl. terrain)
  double t0 = 300.0;
  double p_top = 0.0;
  Array1 znu, znw, dnw, dn, rdnw, rdn, fnm, fnp;
  Array1 c1h, c2h, c3h, c4h, c1f, c2f, c3f, c4f;
  double cf1 = 0, cf2 = 0, cf3 = 0;  // surface extrapolation from the lowest 3 mass levels
  double cfn = 0, cfn1 = 0;          // top extrapolation from the highest 2 mass levels
  double rdx = 0, rdy = 0;
  // Fill the 1D vertical arrays and coefficients from znw (size nz+1), hybrid_opt = 0.
  void set_levels(const std::vector<double>& znw_in);
  void upload_levels();
};

}  // namespace cppwrf::dyn
