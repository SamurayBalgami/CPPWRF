# cppwrf_real: WRF real.exe in C++

`cppwrf_real` makes `wrfinput_d0N` and `wrfbdy_d01` from WPS `met_em` files, as WRF v4.8.0's
`real.exe` does. Its output files are identical to real.exe's: every value (tolerance 0, for the
default single-precision build and the `-r8` build), the variable list and order, dimensions, types,
attributes and global attributes. The one exception is `TITLE`, which names the program.

## Usage

```bash
build/src/init_real/cppwrf_real --namelist namelist.input --met-dir /path/to/met_em \
    --out-dir . --precision float
```

| Option | Default | Meaning |
|---|---|---|
| `--namelist` | `namelist.input` | WRF namelist (the same file real.exe reads) |
| `--met-dir` | `.` | Directory of the `met_em.d0N.<date>.nc` files |
| `--out-dir` | `.` | Where `wrfinput_d0N` and `wrfbdy_d01` are written |
| `--precision` | `double` | WRF's default REAL kind: `float` (standard build) or `double` (`-r8`) |
| `--schema` | `${CPPWRF_DATA_DIR}/wrf_schema/wrf-4.8.0.tsv` | WRF I/O schema (below) |

Domain 1 is processed for every time from `start_date` to `end_date` with `wrfbdy_d01`. Nests are
processed for their first time only, as real.exe does.

## What is ported

`src/init_real` follows `main/real_em.F` and `dyn_em/module_initialize_real.F` statement by statement:

- `real_core.hpp`: vertical interpolation (`vert_interp`, `lagrange_setup`, `lagrange_interp`), `p_dry`,
  `p_dts`, `p_dhs`, `integ_moist`, `sfcprs3`, `rh_to_mxrat1`, `compute_eta` (including
  `auto_levels_opt = 2`), and `compute_vcoord_1d_coeffs`.
- `real_domain.cpp`: the met_em input of `input_auxinput1` and `optional_input`, and `init_domain_rk`.
  That covers the base state, the hydrostatic dry-pressure loop, the soil levels and soil temperatures,
  the sea-ice adjustments and the land/water consistency checks.
- `real_bdy.hpp`: `couple` (with `calc_mu_staggered`), `stuff_bdy_new` and `stuff_bdytend_new`.

The arithmetic matches gfortran's. That includes the libmvec vector lanes of `log`/`exp`, the expansion
of constant powers, and the kind of each Fortran literal.

## Output schema

Which variables WRF writes to the input stream, and with which names, dimensions and attributes,
depends on the namelist through the Registry's package conditions. `tools/wrf_schema/gen_wrf_schema.py`
extracts this from a configured WRF tree into `data/wrf_schema/wrf-<version>.tsv`. The extracted data
are:

- the state variables (`inc/allocs_*.f90`) and their stream bits;
- the package conditions;
- the 4-D array members (`inc/scalar_indices.inc`);
- the Registry defaults.

`wrf_schema.cpp` evaluates the conditions for a namelist. `real_output.hpp` writes the files the way
`share/output_wrf.F` and `share/wrf_bdyout` do. That includes the NetCDF-4 chunking, shuffle and
deflate, and the boundary file's order: a state variable's 8 arrays together; for each 4-D array, all
members' values, then all members' tendencies.

## Tests

`tools/init_real/make_real_refs.sh WRF_SP_DIR WRF_R8_DIR` runs WRF's `real.exe` for every
`${CPPWRF_WRF_RUNS}/run_<case>`. Cases named `*_r8*` use the `-r8` build. The results go to
`${CPPWRF_WRF_RUNS}/real480/<case>`. `tests/init_real/run_real_test.sh [BUILD_DIR]` then runs
`cppwrf_real` on each case and compares the files with `tests/init_real/compare_real.py`. `ctest -R real_vs_wrf` runs the same test. The test is skipped when `${CPPWRF_WRF_RUNS}/real480` does not exist.

The reference set covers 41 cases:

- single domains with Kessler, WSM6 and Thompson microphysics;
- Noah and no LSM;
- `hybrid_opt` 0 and 2;
- 6 h and 12 h boundaries;
- two-domain nests.

All 41 are identical.

## Not covered

The following real.exe paths are not ported. `cppwrf_real` stops with an error when the namelist
selects one:

- `sf_surface_physics` other than 0 and 2 (RUC, Noah-MP and CLM soil handling);
- `use_tavg_for_tsk`;
- `chem_opt` and `tracer_opt` (the 4-D chemistry and tracer arrays);
- `vert_refine_method`.
