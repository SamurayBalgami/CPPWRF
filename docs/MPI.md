# Domain decomposition: patches, tiles and MPI ranks

Sections 1 to 3 are the design, written before the code (with the reasons and the rejected alternatives); section
4 is how to build and run, section 5 what has been proved and how to repeat it, section 6 what has not.

## 1. Goal

One domain (and its nests) runs as N *patches*. A patch is a rectangle of the horizontal grid with its own memory
(its points plus a halo), its own device and its own execution context. The output is bit-identical to WRF real*8
for any patch count and layout, and so identical to the run with one patch.

A patch is run by a *worker*. There are two kinds of worker, and they share every line of the decomposition code:

| Worker | Selected by | Transport between workers |
|---|---|---|
| tile: a host thread of one process | `run.tiles = PxQ` / `CPPWRF_TILES` | in-process mailboxes |
| rank: an MPI process | a `CPPWRF_MPI` build under `mpirun -n N` | `MPI_Isend` / `MPI_Recv` |

Tiles and ranks are not combined (an MPI run with `run.tiles` is an error).

## 2. Why exact decomposition is possible here

WRF dmpar is not bit-identical to serial WRF: its loop bounds are the patch's, so gfortran's vectorised loops put
different points in the SIMD body and in the scalar remainder, and glibc's libmvec lanes differ from the scalar
libm in the last bit (docs/PERF.md, "dmpar"). CPPWRF's kernels compute every point from its global index:
`cppwrf::vm::pow(x, y, i, n)` takes the lane from the position in *WRF's serial loop*, never from the launch
range. So a point computed on any patch gets the operations of the serial run, provided its inputs are the serial
run's. The decomposition therefore only has to (a) give every kernel the inputs of the serial run and (b) keep
every domain-wide operation (sums, maxima, scans in WRF's loop order, I/O order) independent of the layout.

## 3. Architecture

### 3.1 Patches and ownership

`P x Q` patches. Patch `t = tx + P ty` owns the mass points `[tx nx / P, (tx + 1) nx / P) x [ty ny / Q,
(ty + 1) ny / Q)` and the faces with the same indices. The patches at the domain's edge also own the last face
(`nx`, `ny`) and the lateral (physical-boundary) halo beyond the edge. Every point of every array of the serial run
has exactly one owner.

### 3.2 Memory: a window with global indices

A field of a patch is allocated for its *window*: the owned points plus `Grid::halo` (3) points on each side.
The view's origin is shifted so that `f(i, j, k)` takes the **global** indices: kernels, loop bounds, boundary
distances, `vm::` lanes and WRF's index-dependent branches are unchanged. `Grid` carries the patch (`Grid::win`),
`Layout::make` allocates the window, and `Layout` records the window so that host code can loop over it.

Rejected: *local indices* (a patch as a small domain with `nx = patch width`). Every use of an index as a position
(boundary zones, `i == ide - 1` branches, libmvec lanes, staggering at the last face) would need a translation, in
90 000 lines of verified kernels. Rejected: *whole-domain arrays in every worker* (P03 phase 1): memory and I/O
not divided, which is the point of this task. Rejected: *whole-domain virtual memory with untouched pages*: not
possible on a GPU and it hides out-of-patch writes instead of preventing them.

### 3.3 Launches

`tile::ctx` is on for the whole life of a patch (P03: only inside the step). Every launch over `(i, j)` is clipped
to the points the patch owns (edge patches: open-ended, as the serial range). Launches whose index space is not
`(i, j)` (halo lines, boundary rings, edge copies, column batches) clip themselves to the window; `tile::Whole`
keeps meaning "this launch is not clipped to the owned points", and each of its uses states which points it
writes. Host loops over the domain (`fill_host`, the start-up code) loop over the window.

### 3.4 Halo exchange

Semantics (unchanged from P03, which the suite proves): *a patch receives the points other patches own inside its
window*. It is collective: every worker calls it with the same fields in the same order, so messages need no
field identifiers. For each (receiver, sender) pair both sides compute the same boxes from the patch table; the
sender packs them (a kernel on its device), the transport moves bytes, the receiver unpacks.

A boundary fill (`fill_halo*`) is: exchange; then per direction, x first, then y (WRF's order, which fixes the
corners):

- direction not periodic, or periodic with one patch along it: the serial halo-line kernel, on the window's
  tangential range, writing only halo points inside the window;
- direction periodic with several patches along it: a *wrap exchange*: halo point `i` receives the value of
  point `i mod n` from its owner, over the window's whole tangential range (so the y pass carries the x halos, as
  the serial y pass does).

Proof that this equals the serial fill: the serial fill is `value(i, j) = s_x s_y v(m_x(i), m_y(j))` with index
maps `m` and signs `s` per direction (the wall face keeps its value). The exchange delivers every interior point of
the window; each pass then applies one direction's map to rows/columns that already hold the serial values.

GPU: pack and unpack are kernels. The bytes are staged through pinned host buffers by default (works with any
MPI); a CUDA-aware MPI may take the device buffers (`CPPWRF_MPI_CUDA_AWARE=1`, **unproved**: no such MPI and no
second GPU are available to us).

Rejected: one-sided pulls from the neighbours' arrays (P03's `cudaMemcpy3D` between tiles' copies): needs every
worker's memory in one address space.

### 3.5 Domain-wide operations

| Operation | How it stays layout-independent |
|---|---|
| maximum (CFL, adaptive time step, map factor maxima) | maximum over own points, then over the workers (exact) |
| sums (diagnostics) | the summands go to the root in global order; the root adds them in the serial order |
| "any point in the row" / scans in WRF's loop order (Tiedtke rows, Noah's carried RC) | a collective per-row reduction / west-to-east hand-over |
| a value at one point (ideal `u_base` = column (0, 0)) | the owner sends it to all |
| fatal stops | each worker's code; the smallest non-zero code wins everywhere (WRF's first stop) |

### 3.6 Input and output

Output (wrfout, wrfrst, auxhist, native snapshots and restart, diag.csv, run.json): **one writer**, patch 0. The
code that builds WRF's output (`WrfHistory`, the restart writers) runs on the writer only, inside a
`tile::Service` scope: the other workers wait in that scope and answer the writer's requests, each a whole field
gathered from its owners (`tile::download`) or the writer's values spread to the workers that hold the points
(`tile::upload`). The file is the serial run's because the same writer code gets the same arrays in the same
order. The writer holds one whole variable at a time, as WRF's monitor process does (and the few whole arrays
`WrfHistory` keeps between calls). Parallel netCDF was rejected: it would need the netCDF-4 parallel layout and a
proof of identical bytes for no gain in correctness. Running `WrfHistory` on every worker with gathers inside was
rejected too: every branch of that code would have to be symmetric on workers that hold no whole arrays.

Input: a real-data start reads only the worker's window (patch and halo) of wrfinput's 3-D variables
(`io_netcdf::WArray`: a window of a variable under the whole variable's flat index, so the start's code is the
one-patch code) and computes the base state on it. What is still whole on every worker: wrfinput's 2-D variables
during the start (1 / nz of a 3-D field each), the wrfbdy records and the lateral boundary slabs (whole edges:
perimeter, not area), a restart file's fields one at a time while they are read, and the ideal cases' initial
arrays during the start. The writer also holds the output frame it stages for the NetCDF file (every variable of
the frame, as a run with one patch does). None of this changes values.

### 3.7 Nests

Each domain has its own decomposition, the same `P x Q` over its own points; `tile::Use` selects the domain whose
patch clips the launches. Force, feedback and the nest's initialisation read one grid while writing the other, and
a worker's patch of the parent and its patch of the nest cover different areas. The kernel that writes a domain
runs on the worker's patch of that domain and reads a **copy of a box** of the other domain's field
(`tile::WholeCopy` with a box, `nest::WholeField`): the points of that field the worker's kernel reads, each from
its owner (`tile::gather_box`: every worker computes every worker's box, so a worker sends each other worker the
part of its own points inside that worker's box, and nothing where they do not meet). The copy is allocated for
the box and indexed with the domain's indices, so the interpolation and averaging kernels are unchanged.

- Forcing: of the nest's boundary slabs (whole edges on every worker) a worker computes the entries of the nest's
  points it holds, its patch and halo (the relaxation reads the entries of a point's neighbours). The nest's own
  field is the worker's, its halo exchanged after the coupling; the parent's box is the area under those points
  with the 5 x 5 stencil.
- Feedback: the worker's patch of the parent reads the nest's box over that patch (two parent cells around it); a
  worker whose patch does not meet the nest receives nothing. The smoother reads neighbours on the parent: the
  parent field's halo is exchanged after every pass.
- The nest's initialisation (once) and the history fields' feedback (the writer finds the pairs of fields and names
  them to the other workers, `tile::field_index`) read whole copies (`tile::gather`).
- The coupling of the nest forcing averages `mu_2` across the patches' edges, so `mu_2`'s points of the other
  workers are exchanged before it; its lateral halo stays as the dycore's schedule left it, as in WRF.

Rejected: whole copies of both fields for every forced and fed-back field (the first implementation: exact, but
every worker received every field of both domains at every nest step, and a nested run slowed down with more
ranks); one-sided reads of another worker's memory (tiles only).

### 3.8 What a build without MPI keeps

`CPPWRF_MPI` is off by default. Without it no MPI header is included and the transport is the in-process one only.
With one patch nothing is exchanged, the window is the whole domain and every range is the serial range.

## 4. Building and running

```bash
# CPU, MPI ranks (any MPI with the C API: MPICH, Open MPI); CPPWRF_MPI is OFF by default
cmake -B build-mpi -DCMAKE_BUILD_TYPE=Release -DCPPWRF_DOUBLE=ON -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_MPI=ON
cmake --build build-mpi -j --target cppwrf
OMP_NUM_THREADS=4 mpirun -np 6 build-mpi/cppwrf CONFIG.ini                    # 2 x 3 patches (P <= Q, most square)
CPPWRF_RANKS=3x2 OMP_NUM_THREADS=4 mpirun -np 6 build-mpi/cppwrf CONFIG.ini   # or run.ranks = 3x2
# CUDA with MPI ranks: all ranks on one GPU, or CPPWRF_RANK_DEVICES=d0,d1,... (rank r on device d[r mod n])
cmake -B build-cuda-mpi -DCMAKE_BUILD_TYPE=Release -DCPPWRF_DOUBLE=ON -DCPPWRF_ENABLE_CUDA=ON -DCPPWRF_MPI=ON
mpirun -np 2 build-cuda-mpi/cppwrf CONFIG.ini
# tiles: the same decomposition as host threads of one process, any build
CPPWRF_TILES=2x2 build/cppwrf CONFIG.ini
```

- Ranks times `OMP_NUM_THREADS` should not exceed the cores (CPU build). The result does not depend on either.
- Rank 0 writes every output file and the log; the other ranks' standard output is discarded, their errors go to
  standard error, and an error on any rank ends the run (`MPI_Abort`).
- A patch must be at least 3 points wide. A 2-D (x-z) domain needs `Nx1`, a y-z domain `1xN`.
- A nested run uses the same `P x Q` for every domain.
- The run prints `cppwrf: peak host memory: ... MiB the writer (rank 0), ... MiB the largest other rank, ... MiB in all`.
- The MPI binary without `mpirun`, or with one rank, is the one-process program.
- Debugging: `CPPWRF_TILE_TRACE=1` prints every message of the collectives; `CPPWRF_TILE_CHECK=1` (tiles, CPU
  build, no nests) compares every launch of every step with a run with one patch and names the first field and
  point that differ.

## 5. What has been proved

All runs are on one Colab G4 runtime (AMD EPYC 9B45, 48 vCPUs, 176 GB; NVIDIA RTX PRO 6000 Blackwell), MPICH
4.2.1 of Debian 13 (CPU) and MPICH 4.2.0 of the runtime's Ubuntu (GPU), with `tools/colab/colab_mpi.sh` and
`tools/bench/colab_bench.sh --scaling`. "Identical" means: the suite's comparison of every variable of every
wrfout frame with the WRF real*8 reference at tolerance 0 passes, as for a run without MPI, and every output file
of the layout (wrfout, wrfrst, auxhist, native snapshots and restarts; not the logs and the files that hold wall
times) has the sha256 of the run with one rank.

**CPU FP64, MPI ranks, every layout identical:**

| Cases | Layouts (P x Q) | Rank counts |
|---|---|---|
| R1 (dry dynamics), P4 (WSM6, KF, YSU, Noah, RRTMG), MYNN, N3 (two-way nest, feedback, smooth_option 2), em_quarter_ss (periodic x and y, TKE), em_squall2d_x (2-D: N x 1), the native restart check of MYNN, Morrison 1 h restarted from WRF's wrfrst in the middle | 1, 1x2, 1x3, 3x1, 2x2, 1x5, 5x1, 2x3, 3x2, 2x4, 4x2, 8x1, 1x8 | 1, 2, 3, 4, 5, 6, 8 |
| 53 suite cases: the physics suites P1 - P9, MYNN, Morrison, ACM2, Noah-MP, the slab LSM, Thompson / Tiedtke / MYJ / GWD (WRF's default namelist, 6 h, also on the domain with land ice and sea ice), the adaptive time step, grid nudging (also nested), sst_update, the skin-temperature options, auxhist3, restarts from WRF's wrfrst (also nested), five nested cases, six idealised cases (open, periodic, 2-D in x and in y), w_damping, the drop-in runs | 1, 2x3, 4x2, 1x5 | 1, 5, 6, 8 |
| P4, L (359 x 269 x 44, a domain the suite has never seen), N3, 1 h with history: every layout's wrfout against the one-process run's | 1x2, 2x2, 2x4, 3x4, 4x6, 6x8 | 2 to 48 |

**CUDA FP64.** On the workstation's GPU, thread tiles (the same engine and the same pack / unpack kernels as MPI
ranks): squall2d_x with 1, 2, 3 tiles, quarter_ss 2x2, 1x3, 3x1 and P4 2x1, 2x2 are byte-identical to the CUDA run
with one patch, whose snapshots and wrfout are the CPU build's. The CUDA build with `CPPWRF_MPI` compiles and
links on the Colab runtime, and under `mpirun -np 1` 8 cuda_dp suite cases pass against the WRF references; the
CUDA FP32 build as one process writes the base revision's files (5 of 6 cases every file; MYNN's native restart
differs, a file that carries the build's commit). **Several MPI ranks on one GPU are not proved**: the runtime's
Ubuntu MPICH package (4.2.0, built for PMIx) starts every process as a world of one under its own `mpirun`, so the
"2 rank" runs there were two one-process runs writing the same files. A CPPWRF_MPI build now stops when
`CPPWRF_RANKS` asks for several ranks and the process is alone, and `tools/colab/mpi_gpu.sh` uses the runtime's
Open MPI (whose `mpirun` starts two ranks there); that run has not been made.

**Deliberate one-line bugs** (a test counts only if a bug makes it fail):

| Bug | Caught by |
|---|---|
| halo exchange between workers off | R1 2x2: 4 of 5 output files differ |
| one exchange of the feedback smoother removed | N3 2x2: 9 of 18 files differ |
| `mu_2` not exchanged before the nest coupling (found by the first MPI run, not planted) | N3 1x5: d01 V differs in the last bit at one point after 2 min |
| Noah's carried RC not handed across the workers | the land-ice domain (D0V) 3x2 with wrfrst files: RS differs (suite case `real_d0v_r8_fp64_cpu_wrst` under `colab_mpi.sh`) |
| Tiedtke's row flag not combined over the workers | D0 and D0V, 6 h, 8x1 and 4x2: U, V, W differ from 1 - 2 h on (1x8, with whole rows per rank: identical, as it must be) |

Thread tiles (the same engine) run the suite's `*_tPxQ` cases; AddressSanitizer runs of P4 2x2 and N3 2x2, 4x2,
1x5 are clean.

**Strong scaling against WRF v4.8.0 real*8 dmpar** (same runtime, same MPICH, WRF built by its GNU configure;
stepping seconds per simulated hour, history off; `N x 1`: N ranks with one thread each, as WRF dmpar):

| Case | Ranks | WRF dmpar | CPPWRF N x 1 | CPPWRF / WRF speed | writer MiB | largest other rank MiB |
|---|---|---|---|---|---|---|
| P4 (466 k cells) | 1 | 168.7 | 164.0 | 1.03x | - | - |
| | 2 | 85.4 | 82.5 | 1.03x | 1649 | 861 |
| | 4 | 45.2 | 44.5 | 1.02x | 1375 | 550 |
| | 8 | 27.2 | 23.5 | 1.16x | 1156 | 354 |
| | 12 | 16.6 | 17.3 | 0.96x | 1056 | 273 |
| | 24 | 10.3 | 11.2 | 0.92x | 1040 | 225 |
| | 48 | 13.8 | 11.2 | 1.24x | 1021 | 155 |
| L (4.25 M cells) | 6 | 87.2 | 142.2 | 0.61x | 6740 | 1953 |
| | 12 | 53.9 | 71.1 | 0.76x | 5945 | 1288 |
| | 24 | 44.9 | 49.4 | 0.91x | 5660 | 784 |
| | 48 | 50.0 | 53.4 | 0.94x | 5404 | 717 |
| N3 (two-way nest) | 4 | 123.9 | 139.6 | 0.89x | 1850 | 990 |
| | 12 | 49.3 | 58.0 | 0.85x | 1365 | 530 |
| | 24 | 35.8 | 41.8 | 0.86x | 1323 | 428 |
| | 48 | 44.6 | 47.4 | 0.94x | 1236 | 380 |

One 48-vCPU runtime (AMD EPYC 9B45), MPICH 4.2.1 over shared memory. One process with 48 OpenMP threads: P4
11.0 s, L 52.6 s, N3 37.6 s. The memory of a rank other than the writer falls with the rank count (L: 1953 MiB
with 6 ranks, 717 MiB with 48); the writer holds the output frame it stages and whole 2-D inputs. With whole
copies in place of box fetches the nest ran at 0.49x of WRF with 48 ranks (89.6 s).

WRF's own results depend on its rank layout (docs/PERF.md); CPPWRF's do not, which the comparisons above check.

## 6. What has not been proved

Several GPUs and several nodes were not available; every run above is on one machine. A CUDA build has not run
as several MPI ranks at all (section 5): its decomposition is proved with thread tiles on one GPU only, and MPI
itself only with MPICH on the CPU build. Open MPI has not run.
The design should hold there because (a) a rank's GPU is its own CUDA context and stream as a tile's is
(`CPPWRF_RANK_DEVICES` only changes which device the context is on), (b) the transport
carries bytes and nothing in the numerical path depends on where a rank runs, (c) the result is independent of
message timing: every exchange is a fixed sequence of blocking receives of fully determined boxes. Not run:
`CPPWRF_RANK_DEVICES` with more than one device, an MPI run across nodes (a network transport, a file system
shared by the ranks at the start), Open MPI (the tools take `--mpi openmpi`; every run here used MPICH 4.2),
more than 48 ranks.
