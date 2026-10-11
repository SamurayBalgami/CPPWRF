# CPPWRF memory-traffic budget

Model: each 3D field touched by a kernel is read (or written) once from DRAM per kernel,
i.e. stencil neighbours are assumed to hit cache; 1D base-state arrays and 2D fields are
negligible (O(1/nz)). Bytes per cell for FP32 (x2 for FP64). These are **counts from the
code, not measurements**; measured ns/cell/step will be added once M1 runs.

## Large step, per RK stage
| Kernel | reads (3D) | writes (3D) | B/cell FP32 |
|---|---|---|---|
| halo fills (8 prognostics/diagnostics + 6 derived) | halo only | halo only | ~0 (O(6/nx)) |
| KMu (mut, muu, muv) | 2D | 2D | ~0 |
| KPrep3 (ru, rv, rw, alt, php) | u v w al ph | ru rv rw alt php | 40 |
| KWW (calc_ww_cp, column) | ru rv (x2 passes) | ww | 12–20 |
| KTendU (advect_u + PGF [+tendf]) | u ru rv ww p ph alt al php [tendf] | ru_tend | 40–44 |
| KTendV (3D only) | v ru rv ww p ph alt al php [tendf] | rv_tend | 40–44 |
| KTendW (advect_w + pg_buoy_w) | w ru rv ww p [tendf] | rw_tend | 24–28 |
| KTendT (advect_scalar theta) | t ru rv ww [tendf] | t_tend | 20–24 |
| KTendPh (rhs_ph) | u v w ww ph | ph_tend | 24 |

Per stage: x-z (no KTendV): ~176 B/cell; 3D: ~220 B/cell.

## Once per time step
| Item | B/cell FP32 |
|---|---|
| copy *_2 -> *_1 (u v w t ph; mu is 2D) | 40 |
| diffusion, stage 1 only (KDiffU/V/W/T: field + alt -> tendf) | 36 (x-z) / 48 (3D) |

## Total large-step traffic per full time step
x-z: 3 x 176 + 40 + 36 ≈ **604 B/cell** (FP32); 3D: 3 x 220 + 40 + 48 ≈ **748 B/cell**.
At 336 GB/s this is a lower bound of ~1.8 ns/cell/step (x-z) for the large step alone on the
RTX 2060; the acoustic substeps (n_sound/2 + n_sound + 1 substeps) come on top.

## Known inefficiencies (to revisit after M1 is correct)
- KTendU/V/W/T recompute both face fluxes per cell (2x flux arithmetic, no scratch arrays).
- Separate halo-fill launches per field (14 per stage x up to 2 launches).
- `copy_dev` on CPU is a single-threaded memcpy.
