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

2. **Parent Smoothers (`sm121` and `smdsm`):**
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

## 6. Nest Initialization (`input_from_file = true`): `nest_init.hpp`

Matches `share/mediation_integrate.F:509` (`med_nest_initial`) and `dyn_em/nest_init_utils.F`:

1. `med_interp_domain`: Interpolates parent coupled fields to nest scratch arrays (`ht_fine`, `t_save`, `qv_save`, etc.) via `interp_fcn_sint`.
2. `blend_terrain`: Blends `ht_fine` (parent-interpolated) into fine `ht` (from `wrfinput_d02`) in the outer `spec_bdy_width` zone (setting the outermost 2 rows/cols equal to `ht_fine`, and linear blending weights $\frac{5 - b}{3}$ for $b \in \{2, 3, 4\}$).
3. `copy_3d_field`: Copies parent-interpolated `t_save` and `qv_save` into the outer `spec_bdy_width` zone of `t_2` and `qv`, leaving the interior fine `wrfinput_d02` values untouched.
4. `adjust_tempqv`: Hydrostatically adjusts interior nest potential temperature `t_2` using the height difference $\Delta z = \text{ht}_{\text{fine}} - \text{ht}$ and lapse rate $-0.0065\,\text{K/m}$, converting between moist potential temperature $\theta_m = (\theta + T_0)(1 + 1.61\,q_v) - T_0$ and dry $\theta$ when `use_theta_m == 1`.

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
