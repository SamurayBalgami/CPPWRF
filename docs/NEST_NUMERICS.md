# CPPWRF Nesting Numerics (`src/nest/`)

This document specifies the exact mathematical formulations, index conventions, staggering rules, coupling/uncoupling reciprocals, and differential test results for the CPPWRF static nesting numerics module (`src/nest/`), verified against WRF 4.6.1 double-precision Fortran object code (`~/wrf_real/WRF-4.6.1-r8`).

---

## 1. Architecture & Files (`src/nest/`)

All kernels are implemented in single-source C++17 / CUDA (`KOKKOS_FUNCTION`) in `src/nest/`:

- [`nest_types.hpp`](../src/nest/nest_types.hpp): `NestParams` (`i_parent_start`, `j_parent_start`, `parent_grid_ratio`, `parent_time_step_ratio`, `spec_bdy_width`, `smooth_option`, coarse/nest `ids..jde`, `kds..kde`) and `BdySide` (`XS=0`, `XE=1`, `YS=2`, `YE=3`).
- [`nest_interp.hpp`](../src/nest/nest_interp.hpp): 2D 4-point Smolarkiewicz (`sint` / `sintb` / `interp_fcn_sint`) horizontal interpolation (`share/sint.F`, `share/interp_fcn.F:42`).
- [`nest_bdy.hpp`](../src/nest/nest_bdy.hpp): Lateral boundary value (`*B`) and tendency (`*BT`) slab filling (`bdy_interp1`, `share/interp_fcn.F:2423`) and `dyn::LateralBC` helpers (`force_nest_lateral_bc`, `force_nest_lateral_bc_scalar`).
- [`nest_feedback.hpp`](../src/nest/nest_feedback.hpp): Two-way cell-averaging feedback (`copy_fcn`, `share/interp_fcn.F:1397`) for odd and even refinement ratios, plus post-feedback parent smoothers `sm121` (`smooth_option = 1`, 1-2-1 filter) and `smdsm` (`smooth_option = 2`, smoother-desmoother).
- [`nest_couple.hpp`](../src/nest/nest_couple.hpp): In-place column-mass / map-scale coupling and reciprocal-multiply uncoupling (`couple_or_uncouple_em`, `dyn_em/couple_or_uncouple_em.F`).
- [`nest_init.hpp`](../src/nest/nest_init.hpp): Initialization utilities (`blend_terrain`, `copy_3d_field`, `adjust_tempqv`, `dyn_em/nest_init_utils.F`, `share/mediation_integrate.F:509`).

---

## 2. Horizontal Prolongation: `interp_fcn_sint` (`nest_interp.hpp`)

Matches `share/sint.F` and `share/interp_fcn.F:42-357`:

1. **1D 4-point Smolarkiewicz stencil (`smolarkiewicz_1d`):**
   For normalized subgrid offset $x \in [-0.5, 0.5]$ and 4 stencil values $(f_0, f_1, f_2, f_3)$ centered between $f_1$ and $f_2$:
   - Compute 3rd-order Lagrange baseline $y_0$ and forward/backward slope limiters preventing overshoot/undershoot beyond $\max(f_0, f_1, f_2, f_3)$ and $\min(f_0, f_1, f_2, f_3)$ with $\varepsilon = 10^{-15}$.
   - Apply non-oscillatory Flux-Corrected Transport (FCT) fluxes across the 3 intervals when $f_3 - f_0$ and $f_2 - f_1$ share sign.
   - Falls back to 4-point Lagrange interpolation (`fl4`) if $\max - \min \le 10^{-10}$, or linear (`fl2`) if $|x| < 10^{-4}$.

2. **2D Dimension-Split Evaluation (`sint_point_2d`):**
   - **Pass 1 (X sweep):** Interpolate 4 rows $m \in \{0, 1, 2, 3\}$ in $X$ using `smolarkiewicz_1d` at offset $x_{\text{off}}$ to get $g_m$.
   - **Pass 2 (Y sweep):** Interpolate $(g_0, g_1, g_2, g_3)$ in $Y$ using `smolarkiewicz_1d` at offset $y_{\text{off}}$.

3. **Staggering & Subgrid Offset Rules (`compute_sint_subgrid_offsets`):**
   Let $r_x = \text{nri}, r_y = \text{nrj}$ be the refinement ratios and $(i_p, j_p) = (\text{ipos}, \text{jpos})$ the 1-based parent start index:
   - **X-staggered (`xstag = true`):** Coarse parent cell index and offset:
     $$rx = \frac{(i_{\text{nest}} - 1) + 0.5\,((r_x + 1) \bmod 2)}{r_x}, \quad i_c = i_p + \lfloor rx + \tfrac{1}{2} \rfloor, \quad x_{\text{off}} = rx - \lfloor rx + \tfrac{1}{2} \rfloor$$
   - **X-unstaggered (`xstag = false`):**
     $$rx = \frac{i_{\text{nest}} - 1 + 0.5}{r_x} - 0.5, \quad i_c = i_p + \lfloor rx + \tfrac{1}{2} \rfloor, \quad x_{\text{off}} = rx - \lfloor rx + \tfrac{1}{2} \rfloor$$
   - Identical formulas hold in $Y$ with `ystag`, $r_y$, and $j_p$. The 4x4 coarse stencil uses indices $(i_c - 1 + l, j_c - 1 + m)$ for $l, m \in \{0, 1, 2, 3\}$.

---

## 3. Lateral Boundary Slab Forcing: `bdy_interp1` (`nest_bdy.hpp`)

Matches `share/interp_fcn.F:2423`:

- Given `ndt = parent_time_step_ratio` and nest step $\text{cdt} = \text{ndt} \cdot dt_{\text{nest}}$ (coarse parent time step):
  $$\text{rdt} = \frac{1}{\text{cdt}}$$
- For each lateral boundary zone of width `spec_bdy_width = 5` on `XS` ($i = 1..5$), `XE` ($i = n_{ie}, \dots, n_{ie}-4$), `YS` ($j = 1..5$), and `YE` ($j = n_{je}, \dots, n_{je}-4$):
  1. `sint_point_2d` interpolates the parent field onto `nrepl` at the start of the parent step (`ncall == 1`, `spec_zone == 0`), storing `B = nrepl` and initializing `BT = 0`.
  2. After the parent advances by $\Delta t_{\text{parent}} = \text{cdt}$, `bdy_interp1` is called again (`ncall == 2`, `spec_zone == spec_bdy_width`): `sint_point_2d` interpolates the updated parent field onto `nrepl`, and sets:
     $$\text{BT}(i, k, b) = \text{rdt} \cdot \left(\text{nrepl} - \text{B}(i, k, b)\right)$$
     leaving `B(i, k, b)` unchanged as the value at the start of the parent step.

---

## 4. Two-Way Feedback & Smoothing: `copy_fcn`, `sm121`, `smdsm` (`nest_feedback.hpp`)

Matches `share/interp_fcn.F:1397-2258`:

1. **Cell-Averaging Feedback (`copy_fcn`):**
   - Operates over the interior parent feedback window excluding `spec_bdy_width + 1` coarse rows/columns inside the nest footprint:
     $$i_s = i_p + W_{\text{bdy}}, \quad i_e = i_p + \frac{n_{de} - n_{ds}}{r_x} - W_{\text{bdy}} - \begin{cases} 0 & \text{if } \text{xstag} \\ 1 & \text{if } \neg\text{xstag} \end{cases}$$
     (and analogously for $j_s, j_e$).
   - **Staggered points** (`xstag` or `ystag`) use the centered $(2\lfloor r/2 \rfloor + 1)$ cell window with weights $\frac{1}{\text{npts}}$ accumulated term-by-term in the exact Fortran loop order (`dj` outer, `di` inner).
   - **Unstaggered points** (for both odd and even ratios $r_x, r_y$) average over the $r_y \times r_x$ nest cells starting at $n_{j,\text{base}} = (j_c - j_p)\,r_y + 1$, $n_{i,\text{base}} = (i_c - i_p)\,r_x + 1$ with term-by-term multiplication `sum += inv_pts * nfld(...)`.

2. **Point-Sampling Map-Scale & Coriolis Feedback (`copy_fcnm`):**
   - In `inc/nest_feedbackup_interp.inc` (included by `share/mediation_feedback_domain.F`), WRF also feeds back the 11 static 2D map-scale and Coriolis fields (`msftx, msfty, msfux, msfuy, msfvx, msfvx_inv, msfvy, f, e, sina, cosa`) using `copy_fcnm` (`share/interp_fcn.F:1747`) with `.FALSE.` (no post-smoothing).
   - Unlike cell-averaging `copy_fcn`, `copy_fcnm` samples the fine nest grid at the exact center of the coarse cell (`ioff = (nri - 1) / 2`, `joff = (nrj - 1) / 2` for unstaggered; `ioff = 0` for X-staggered; `joff = 0` for Y-staggered):
     $$c(i_c, j_c) = n\big((i_c - i_p)\,r_x + 1 + i_{\text{off}},\, (j_c - j_p)\,r_y + 1 + j_{\text{off}}\big)$$
   - Because `wrfinput_d01` and `wrfinput_d02` store map factors in FP32 (`~1e-7` precision), omitting `copy_fcnm` leaves a `2.4e-7` discrepancy in `d01`'s `msfuy`/`msfvx` inside the feedback window after `t = 0`, which changes `d01`'s `W` diagnostic by `2.6e-8` at `t = 0`. Including `feedback_field2_fcnm` (`KCopyFcnmView2D`) eliminates this discrepancy completely (`0.000e+00` bit difference on `U, V, W` at `t = 0` on `d01`).

3. **Parent Smoothers (`sm121` and `smdsm`):**
   - Applied over $i \in [i_p + 3, i_{\text{end}} - 3]$, $j \in [j_p + 3, j_{\text{end}} - 3]$.
   - `sm121` (`smooth_option = 1`): `npasses = (parent_grid_ratio + 1) / 2` passes of a 2D separable 1-2-1 filter ($Y$-pass into scratch, then $X$-pass, with odd passes sweeping forward in $i$ and even passes sweeping backward in $i$).
   - `smdsm` (`smooth_option = 2`): Two-pass smoother-desmoother with coefficients $x_{\nu} = \{0.50, -0.52\}$ and `npasses = parent_grid_ratio - 1` pairs ($Y$-pass then $X$-pass per $\nu$).

---

## 5. In-Place Coupling & Reciprocal Uncoupling: `couple_or_uncouple_em` (`nest_couple.hpp`)

Matches `dyn_em/couple_or_uncouple_em.F` bit-for-bit:

1. **Column Dry-Air Mass Factors (`mutf`, `mut`, `muu`, `muv`):**
   - **Full levels (`W`, `PH_2`):**
     $$\text{denom}_f(i,k,j) = (c_{1f}(k)\,\text{MUB}(i,j) + c_{2f}(k)) + (c_{1f}(k)\,\text{MU\_2}(i,j))$$
   - **Half levels (`T_2`, `QV`, `QC`, `QR`, `QI`, `QS`, `QG`):**
     $$\text{denom}_h(i,k,j) = (c_{1h}(k)\,\text{MUB}(i,j) + c_{2h}(k)) + (c_{1h}(k)\,\text{MU\_2}(i,j))$$
   - **U-staggered half levels (`U_2`, with west halo clamped $i_{-1} = \max(i-1, 0)$ via `set_physical_bc2d`):**
     $$\text{denom}_u(i,k,j) = \Big(\big((c_{1h}(k)\,\text{MUB}(i,j) + c_{2h}(k)) + (c_{1h}(k)\,\text{MUB}(i_{-1},j) + c_{2h}(k))\big) + (c_{1h}(k)\,\text{MU\_2}(i,j))\Big) + (c_{1h}(k)\,\text{MU\_2}(i_{-1},j))$$
   - **V-staggered half levels (`V_2`, with south halo clamped $j_{-1} = \max(j-1, 0)$ via `set_physical_bc2d`):**
     $$\text{denom}_v(i,k,j) = \Big(\big((c_{1h}(k)\,\text{MUB}(i,j) + c_{2h}(k)) + (c_{1h}(k)\,\text{MUB}(i,j_{-1}) + c_{2h}(k))\big) + (c_{1h}(k)\,\text{MU\_2}(i,j))\Big) + (c_{1h}(k)\,\text{MU\_2}(i,j_{-1}))$$

2. **Coupling vs. Uncoupling Reciprocal Multiplication:**
   WRF does **not** divide on uncoupling; instead it computes a reciprocal multiplier and multiplies in-place:
   - **Mass (`MU_2`):**
     $$\text{factor}_{\text{couple}} = \frac{1}{\text{msftx}\cdot\text{msfty}}, \qquad \text{factor}_{\text{uncouple}} = \text{msftx}\cdot\text{msfty}$$
   - **Full/Half scalars (`W`, `PH_2`, `T_2`, moisture):**
     $$\text{factor}_{\text{couple}} = \frac{\text{denom}}{\text{msftx}\cdot\text{msfty}}, \qquad \text{factor}_{\text{uncouple}} = \frac{\text{msftx}\cdot\text{msfty}}{\text{denom}}$$
   - **Horizontal velocities (`U_2`, `V_2`):**
     $$\text{factor}_{u,\text{couple}} = \frac{0.5 \cdot \text{denom}_u}{\text{msfuy}}, \qquad \text{factor}_{u,\text{uncouple}} = \frac{2.0 \cdot \text{msfuy}}{\text{denom}_u}$$
     $$\text{factor}_{v,\text{couple}} = \frac{0.5 \cdot \text{denom}_v}{\text{msfvx}}, \qquad \text{factor}_{v,\text{uncouple}} = \frac{2.0 \cdot \text{msfvx}}{\text{denom}_v}$$
   Because $(x \cdot \text{factor}_{\text{couple}}) \cdot \text{factor}_{\text{uncouple}}$ incurs up to 1 ULP IEEE-754 rounding drift on both parent and nest every time `med_nest_force` or `med_nest_feedback` runs, performing the exact couple/uncouple round-trip in place is mandatory for bit-reproducibility. Verified against `~/wrf_real/run_N1_r8/steps/wrfout_d02_2026-10-12_00:00:00` frame 0 (written after `med_nest_force`): **0 bit mismatches across all 475,200 `U` grid points**.

---

## 6. Nest Initialization (`input_from_file = true`): `nest_init.hpp` & `nest_driver.cpp`

Matches `share/mediation_integrate.F:509` (`med_nest_initial`) and `dyn_em/nest_init_utils.F`:

1. `med_interp_domain`: Interpolates parent `ht`, `mub`, `phb`, `t_2`, and `qv` onto nest scratch arrays (`ht_int`, `mub_int`, `phb_int`, `t_int`, `qv_int`) via `interp_fcn_sint`.
2. `blend_terrain`: Blends `ht_int`, `mub_int`, and `phb_int` into the nest's `ht`, `mub`, and `phb` (loaded from `wrfinput_d02`) in the outer `spec_bdy_width` zone (setting the outermost 2 rows/cols equal to the parent-interpolated field, and linear blending weights $\frac{5 - b}{3}$ for $b \in \{2, 3, 4\}$).
3. `copy_3d_field`: Copies parent-interpolated `t_int` and `qv_int` into the outer `spec_bdy_width` zone of `t_2` and `qv`, leaving the interior fine `wrfinput_d02` values untouched, and saves the pre-blended `wrfinput_d02` terrain height in `ht_fine`.
4. **Optional Two-Way Feedback (`feedback == 1`):** Calls `med_feedback_domain(parent, child)`, which couples parent and child in-place, copies interior nest prognostic fields (`u_2, v_2, w, ph_2, t_2, mu_2, moist`), base fields (`phb, mub, alb, pb, ht`), and point-sampled map/Coriolis factors (`msftx..cosa`) onto the parent, smooths `u_2, v_2, w, ph_2, t_2, mu_2, mub, phb, ht` if `smooth_option > 0`, uncouples both grids, recomputes parent reference profiles (`p_hyd`, `pb`, `alb`, `t_init`, `al`, `p`) inside the feedback window via `recompute_parent_feedback_ref`, and updates parent kinematic `W` via `init_em_recompute_w`.
5. `adjust_tempqv`: Hydrostatically adjusts interior nest potential temperature `t_2` using the height difference $\Delta z = \text{ht}_{\text{fine}} - \text{ht}$ and lapse rate $-0.0065\,\text{K/m}$, converting between moist potential temperature $\theta_m = (\theta + T_0)(1 + 1.61\,q_v) - T_0$ and dry $\theta$ when `use_theta_m == 1`.
6. `init_domain_constants_em` (`recompute_nest_after_blend`): Recomputes the nest's hydrostatic base state (`phb`, `pb`, `alb`, `t_init`), perturbation `ph_2`, dry air density `alt`/`al`, perturbation pressure `p`, hydrostatic pressure `p_hyd`/`p_hyd_w`, geopotential `z`, and surface diagnostics `p_top`, `psfc`, `mu0` from the blended `ht` and `mub`, followed by `init_em_recompute_w` on the nest.
7. `med_nest_force` (`ncall = 1` at $t = 0$): Couples parent and child in-place, fills the child's lateral boundary value slabs (`*B`, with `*BT = 0`) via `bdy_interp1`, and uncouples both parent and child in-place.

---

## 7. Differential Test Results vs WRF 4.6.1 Fortran (`tests/nest/`)

Command:
```bash
tools/heavy.sh bash tests/nest/run_test.sh
```

### 7.1 CPU Differential Tests (`test_nest_differential` vs `~/wrf_real/WRF-4.6.1-r8` Fortran `.o`)

| Test Group | Sub-tests | Max Rel Error | Max Abs Error | Status |
|---|---|---|---|---|
| 1. `sint` / `sintb` | `sint` (`xstag=0,1`), `sintb` (`xstag=0,1`) | `0.000000e+00` | `0.000000e+00` | PASS (Bit-exact) |
| 2. `interp_fcn_sint` / `KSint3D` | `(F,F)`, `(T,F)`, `(F,T)`, `(T,T)` + `KSint3D` | `0.000000e+00` | `0.000000e+00` | PASS (Bit-exact) |
| 3. `bdy_interp1` / `LateralBC` | `call1` (`B`), `call2` (`BT`) on `XS, XE, YS, YE` for `(F,F), (T,F), (F,T)` | `0.000000e+00` | `0.000000e+00` | PASS (Bit-exact) |
| 4. `copy_fcn` + Smoothers | `copy_fcn` ratios `3, 4, 5` (`(F,F), (T,F), (F,T)`), `sm121`, `smdsm` | `0.000000e+00` | `0.000000e+00` | PASS (Bit-exact) |
| 5. `couple_or_uncouple_em` | `couple` and `uncouple` (`mu_2, u_2, v_2, w_2, ph_2, t_2, moist`) | `0.000000e+00` | `0.000000e+00` | PASS (Bit-exact) |
| 6. `nest_init` | `blend_terrain`, `copy_3d_field`, `adjust_tempqv` (`use_theta_m=0,1`) | `0.000000e+00` | `0.000000e+00` | PASS (Bit-exact) |

### 7.2 CUDA FP64 GPU Kernel Tests (`test_nest_cuda` on NVIDIA GeForce RTX 2060)

| GPU Kernel Test Group | Max Rel Error (GPU vs Fortran/CPU) | Status |
|---|---|---|
| 1. `KSint3D` (`(F,F), (T,F), (F,T)`) | `<= 7.868e-13` (FMA in 4-pt Smolarkiewicz limiter) | PASS |
| 2. `KBdySlab` (`B` and `BT` on `XS, XE, YS, YE`) | `<= 6.469e-16` | PASS |
| 3. `KCopyFcn` (`(F,F), (T,F), (F,T)`) | `<= 3.802e-16` | PASS |
| 4. `couple_or_uncouple_state` (`couple` & `uncouple` all fields) | `0.000000e+00` (Bit-exact) | PASS |
| 5. `nest_init` (`blend_terrain`, `copy_3d_field`, `adjust_tempqv`) | `0.000000e+00` (Bit-exact) | PASS |

---

## 8. Nest Driver (`src/cases/nest_driver.{hpp,cpp}`) & End-to-End Acceptance

### 8.1 Recursive Time-Stepping Order (`frame/module_integrate.F:integrate`)
 implemented in [`src/cases/nest_driver.cpp`](../src/cases/nest_driver.cpp) (`cases::run_nested`) against `cases::Model` and `cases::DomainRunner`:
1. Construct root (`d01`) and child (`d02..d0N`) `DomainRunner` instances in `<out_dir>/d01`, `<out_dir>/d02`, etc.
2. At $t = 0$ on `d01`: write `d01`'s initial output ONLY if `feedback == 0`; otherwise defer `d01`'s initial output until after `med_nest_initial(d01, d02)` has fed back `d02`'s $t = 0$ fields onto `d01` (`share/mediation_integrate.F:med_Before_Solve_io`).
3. For each parent step $s_p = 0, \dots, N_{\text{steps},p}-1$:
   - Advance parent by one time step $\Delta t_p$ (`parent.step_one()`), deferring parent output when child nests are present (`has_children = true`).
   - For each child nest:
     - Call `med_nest_force(parent, child)` (`ncall = 2` at $t_p + \Delta t_p$) to compute lateral boundary tendencies `*BT = (F(t_p + \Delta t_p) - *B) / (r_t \Delta t_c)` while leaving `*B` at $t_p$.
     - Advance child by $r_t = \text{parent\_time\_step\_ratio}$ sub-steps (`integrate_domain(child)`), writing child output at its configured `out_interval`.
     - Call `med_nest_feedback(parent, child)` (if `feedback == 1`), updating parent fields, reference profiles, and kinematic `W` in-place.
     - Refresh child's start-of-window lateral boundary values `*B` at $t_p + \Delta t_p$ (`med_nest_force` with `ncall = 1` at the start of the next parent step).
   - Write parent output (`parent.write_if_due()`) after all child sub-steps and two-way feedback complete.

### 8.2 Acceptance Summary vs WRF 4.6.1 (`run_N1_r8`, `run_N2_r8`, `run_N1`, `run_N2`)

| Case | Mode | Precision / Backend | Domain | Window | Max `relL2` (`U, V`) | Max `relL2` (`W`) | Max `relL2` (`PH, MU, THM`) | Max `relL2` (`P, QV`) | Status |
|---|---|---|---|---|---|---|---|---|---|
| `run_N1_r8/steps` | One-way (`feedback=0`) | FP64 CPU | `d01` | 3 steps (`0..180 s`) | `2.42e-14` | `4.59e-12` | `4.43e-16` | `7.05e-15` | PASS (`<= 1e-10`) |
| `run_N1_r8/steps` | One-way (`feedback=0`) | FP64 CPU | `d02` | 9 steps (`0..180 s`) | `7.83e-14` | `8.06e-12` | `7.00e-16` | `9.60e-15` | PASS (`<= 1e-10`) |
| `run_N2_r8/steps` | Two-way (`feedback=1`) | FP64 CPU | `d01` | 3 steps (`0..180 s`) | `2.44e-14` | `4.48e-12` | `4.29e-16` | `7.09e-15` | PASS (`<= 1e-10`) |
| `run_N2_r8/steps` | Two-way (`feedback=1`) | FP64 CPU | `d02` | 9 steps (`0..180 s`) | `7.70e-14` | `8.05e-12` | `7.07e-16` | `9.71e-15` | PASS (`<= 1e-10`) |
| `run_N1_r8/steps` | One-way (`feedback=0`) | FP64 CUDA | `d01` | 3 steps (`0..180 s`) | `1.80e-13` | `6.30e-12` | `5.95e-16` | `1.50e-14` | PASS (`<= 1e-10`) |
| `run_N1_r8/steps` | One-way (`feedback=0`) | FP64 CUDA | `d02` | 9 steps (`0..180 s`) | `2.09e-13` | `1.04e-11` | `8.14e-16` | `1.62e-14` | PASS (`<= 1e-10`) |
| `run_N2_r8/steps` | Two-way (`feedback=1`) | FP64 CUDA | `d01` | 3 steps (`0..180 s`) | `1.79e-13` | `6.11e-12` | `5.75e-16` | `1.49e-14` | PASS (`<= 1e-10`) |
| `run_N2_r8/steps` | Two-way (`feedback=1`) | FP64 CUDA | `d02` | 9 steps (`0..180 s`) | `1.99e-13` | `1.04e-11` | `8.11e-16` | `1.62e-14` | PASS (`<= 1e-10`) |
| `run_N1_r8` (1 h) | One-way (`feedback=0`) | FP64 CUDA | `d01` | 60 steps (`3600 s`) | `4.77e-13` | `7.10e-12` | `1.02e-15` | `4.06e-14` | PASS (`<= 1e-10`) |
| `run_N1_r8` (1 h) | One-way (`feedback=0`) | FP64 CUDA | `d02` | 180 steps (`3600 s`) | `2.23e-13` | `1.01e-11` | `9.26e-16` | `2.49e-14` | PASS (`<= 1e-10`) |
| `run_N2_r8` (1 h) | Two-way (`feedback=1`) | FP64 CUDA | `d01` | 60 steps (`3600 s`) | `4.78e-13` | `6.78e-12` | `9.75e-16` | `4.03e-14` | PASS (`<= 1e-10`) |
| `run_N2_r8` (1 h) | Two-way (`feedback=1`) | FP64 CUDA | `d02` | 180 steps (`3600 s`) | `2.40e-13` | `1.01e-11` | `9.33e-16` | `2.45e-14` | PASS (`<= 1e-10`) |
| `run_N1` (1 h) | One-way (`feedback=0`) | FP32 CUDA | `d01` | 60 steps (`3600 s`) | `1.91e-04` (`in: 5.09e-5`) | `3.62e-03` | `6.48e-07` | `2.00e-05` | PASS (FP32) |
| `run_N1` (1 h) | One-way (`feedback=0`) | FP32 CUDA | `d02` | 180 steps (`3600 s`) | `2.82e-04` (`in: 9.46e-5`) | `4.72e-03` | `9.67e-07` | `3.14e-05` | PASS (FP32) |
| `run_N2` (1 h) | Two-way (`feedback=1`) | FP32 CUDA | `d01` | 60 steps (`3600 s`) | `1.92e-04` (`in: 5.78e-5`) | `3.46e-03` | `6.36e-07` | `1.96e-05` | PASS (FP32) |
| `run_N2` (1 h) | Two-way (`feedback=1`) | FP32 CUDA | `d02` | 180 steps (`3600 s`) | `2.83e-04` (`in: 9.82e-5`) | `4.76e-03` | `1.01e-06` | `3.20e-05` | PASS (FP32) |

