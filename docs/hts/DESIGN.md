# Hierarchical time stepping for MPAS-Atmosphere

Design for MPAS-Atmosphere **v8.4.2**
(`f34984b2c89353cf31c235b08d3b4acd96b37bb3`, upstream tag `v8.4.2`,
4 September 2026). The fork's `master` is left at the 2024 v8.2.2 history.
Phase 2 implements the dry scheme in `src/core_atmosphere/dynamics/mpas_atm_hts.F`
and in masked loops inside `mpas_atm_time_integration.F`. Baseline numbers
from the unmodified dynamical core are in [BASELINE.md](BASELINE.md).
Measured HTS gates are in [PHASE2.md](PHASE2.md).

This scheme is designed for MPAS-A's split-explicit RK integrator. It is
informed by ClusterTech CPAS (Sze et al. 2022, doi:10.1029/2022EA002342;
Cheung et al. PASC'22, doi:10.1145/3539781.3539790) and by MPAS-Ocean local
time stepping (Capodaglio & Petersen 2022, doi:10.1016/j.jcp.2021.110818;
Lilly et al. 2024, doi:10.1016/j.jcp.2024.113511; E3SM PRs 5460 and 6224).
Ocean LTS is not present in MPAS-Model v8.4.2. The atmosphere algorithm
below is not a port of either code base.

## 1. What the current integrator does

One namelist value, `config_dt`, is the step for every cell. The simulation
clock is built from that single interval in `atm_simulation_clock_init`
(`src/core_atmosphere/mpas_atm_core.F`). On a variable-resolution mesh the
step is the one required by the smallest cells.

Each clock tick calls `atm_do_timestep`, which:

1. If `moist_physics` is true, calls `physics_timetracker` and `physics_driver` for the whole block, using the same `dt`.
2. Calls `atm_timestep` → `atm_srk3` (the only implemented `config_time_integration` value is `SRK3`).

`atm_srk3` (`src/core_atmosphere/dynamics/mpas_atm_time_integration.F`) is a
split-explicit Runge–Kutta scheme (Wicker & Skamarock / Klemp et al.):

- **Dynamics / transport split.** If `config_split_dynamics_transport` is true (the registry default), the acoustic/RK dynamics are repeated `config_dynamics_split_steps` times (default 3) with `dt_dynamics = dt / config_dynamics_split_steps`, and scalar transport runs once on the full `dt` using mass fluxes accumulated over those substeps (`ruAvg`, `wwAvg` via `atm_rk_dynamics_substep_finish`). The Jablonowski–Williamson namelist turns the split off and also turns scalar advection off, so one dynamics pass uses the full `dt`.
- **RK stages.** `config_time_integration_order` is 2 by default (3 is supported). For order 2 the stage advances are `dt_dynamics/2`, `dt_dynamics/2`, `dt_dynamics`. For order 3 they are `dt_dynamics/3`, `dt_dynamics/2`, `dt_dynamics`.
- **Acoustic substeps.** Inside each RK stage, `atm_advance_acoustic_step` integrates the fast horizontally explicit, vertically implicit acoustic system. With `config_number_of_sub_steps = N` and order 2, the three stages take `N/2`, `N/2`, and `N` substeps, and the acoustic step is `dt_dynamics/N` on every stage. The JW baseline uses `N = 6`, so the acoustic step is `dt/6`. Slow tendencies (advection, Coriolis, pressure-gradient, physics) are frozen over the acoustic substeps of a stage.
- **Physics tendencies** (`tend_ru_physics`, `tend_rtheta_physics`, `tend_rho_physics`) are taken once at the start of `atm_srk3` and held for the whole large step. Microphysics, if on, runs at the end of `atm_srk3`, after transport. With `config_physics_suite = 'none'` every scheme is set to `off` and `moist_physics` is false, so the driver call is skipped and the tendency arrays are zeroed.
- **Regional lateral boundaries**, if `config_apply_lbcs` is true, time-interpolate a driving state with `mpas_atm_get_bdy_state` at the stage time. That interpolant is the closest existing pattern to what a coarse-to-fine HTS boundary needs. The global idealized baselines do not use it.

### Halo exchange

Groups are defined in `src/core_atmosphere/mpas_atm_halos.F` and invoked as
`exchange_halo_group(domain, 'dynamics:<fields>', ...)`. Acoustic fields
(`rho_pp`, `rtheta_pp`) are registered with one halo layer. The hot path,
executed on every acoustic substep, is:

1. exchange `rho_pp`
2. `atm_advance_acoustic_step` (horizontally explicit, so it needs the neighbor `rho_pp`)
3. exchange `rtheta_pp`
4. `atm_divergence_damping_3d`

Around that inner loop each RK stage also exchanges `tend_u` before the acoustic loop and `ru_p`, `rw_p`, `rho_pp`, `rtheta_pp`, then `u`, then `w` / `pv_edge` / `rho_edge` after the large-step recovery. There is no overlap of communication with the acoustic update. Every rank calls every exchange; a rank that skipped an exchange would deadlock.

Decomposition is a METIS partition file (`config_block_decomp_file_prefix` + task count), one block per rank. The graph is unweighted, so a fine cell and a coarse cell cost the same in the partition even though a fine cell will do `2^L` times more work under HTS. This build does not link Scotch. OpenMP and OpenACC are off.

### How the mesh describes resolution

Horizontal grids from NCAR are SCVTs on the unit sphere. `nominalMinDc` is the nominal minimum `dcEdge` where `meshDensity == 1`, in units of the sphere radius. Physical spacing is `nominalMinDc * sphere_radius` with the default radius 6371229 m. On the tutorial variable mesh `x5.30210`, `nominalMinDc * R = 48000` m and `(max meshDensity / min meshDensity)^0.25 = 5` exactly, i.e. mesh density is stored as `(dx_min / dx)^4`. That matches the del4 scaling in `atm_compute_mesh_scaling`:

```
meshScalingDel4 = 1 / ((meshDensity(cell1) + meshDensity(cell2))/2)^0.75
```

which is `(dx / dx_min)^3` when `config_len_disp` is the minimum spacing. `config_h_ScaleWithMesh` defaults to true. The baseline therefore sets `config_len_disp` to the nominal finest spacing (120000 m on `x1.40962`, 48000 m on `x5.30210`).

## 2. Proposed HTS scheme

`config_dt` stays the fine-level step. Coarser levels take `dt * 2^L`. The outer clock advances by the fine step so that physics alarms and history output, which are clock-based, do not need a new time type. Coarse columns are simply not updated on fine ticks that are not their tick.

### Level assignment

At dynamics init, for each cell:

```
dx = nominalMinDc * sphere_radius / meshDensity^(1/4)
L  = nint(log2(max(dx / dx_min, 1)))
```

clamped to `config_hts_max_level` (registry default 1, hard cap 8). The default is 1 so the first variable-mesh test is two levels; raise it to 2 for the three-level run. The 2019 `x1.40962` grid used for the uniform baseline has no `nominalMinDc` variable (init warns and continues). `meshDensity` is present and is identically 1, so that mesh is one level either way. Level assignment keys off `meshDensity` and treats a missing `nominalMinDc`, or a mesh whose density does not vary, as "all cells level 0." An edge's level is the finer (smaller) of its two `cellsOnEdge` levels: the acoustic CFL at that edge is set by the shorter adjacent spacing. An edge whose cells differ in `L` is an interface edge.

Power-of-two bins are required so a coarse step contains an integer number of fine steps. A mesh whose refinement ratio is not a power of two (the 5× tutorial mesh) is binned conservatively: a 230 km cell on a 48 km mesh lands in `L = 2` (`4 dt`, not `5 dt`) and is integrated slightly more often than its own CFL limit. On `x5.30210` that assignment is:

| level | cells | fraction | mean dx | steps relative to fine |
| --- | ---: | ---: | ---: | ---: |
| 0 | 16652 | 55.1% | 51 km | 1 |
| 1 | 3729 | 12.3% | 91 km | 1/2 |
| 2 | 9829 | 32.5% | 230 km | 1/4 |

Interface edges (the two cells have different `L`) are 1223 of 90624 (1.3%). The work if every level takes its own step is 0.69 of an all-fine integration, a ceiling of about **1.44×** on this mesh. Most cells are fine, so this mesh is a correctness test, not a speed showcase. A mesh whose fine region is a small fraction of the cells (the published 60–3 km mesh is the scientific example; it is too large for this VM) is the right speed test later.

Uniform meshes have `meshDensity = 1` everywhere, so every cell is level 0. That is the regression test: HTS enabled with one level must reproduce the baseline.

### Time march

`config_hts = false` calls the unmodified `atm_srk3`. `config_hts = true` calls `atm_hts_advance` once per fine clock tick. The outer clock still advances by `config_dt`.

On a tick, level L is at the start of its interval when `mod(tick, 2^L) = 0` and is due to advance when `mod(tick, 2^L) = 2^L - 1`. The advance is `dt * 2^L` through one call to the existing RK3 integrator, masked to that level. With one level (`hts_max_level_used = 0`) there is a single masked call and the coarse machinery is not allocated.

Order on a tick with more than one level:

1. Snapshot any level whose interval is opening (`mod(tick, 2^L) = 0`). The snapshot is the time-level-1 state at the start of that level's step. The previous-step increment is stored but not applied (see below).
2. Write that snapshot into every cell and edge coarser than the level about to run, on both time levels and on the diagnostics the acoustic solver reads. Reconstruct the horizontal wind, then call `atm_srk3` masked to the active level.
3. Restore every coarser cell and edge to its snapshot, so the fine call cannot commit the extrapolant.
4. Add `ruAvg * dt` on level-0 interface edges into a double-precision accumulator.
5. For each coarser level that is due, from fine to coarse: copy the finer time-level-2 state back onto time level 1 and onto the saved acoustic base state; set the imposed `ruAvg` target to (accumulated flux) / `dt_L` on interface edges whose coarser cell is this level; hold levels coarser than L at their snapshot; call `atm_srk3` with `dt_L`; restore those coarser levels; patch any flux residual (below); accumulate this level's own outgoing interface flux; clear the accumulator entries this level consumed.

A rejected alternative was a predictor: integrate the coarse region once to get the interval endpoint, interpolate, then integrate it again with the fine region. That integrates the coarse region twice and spends the speedup. Linear extrapolation of the previous coarse increment (`alpha = (mod(tick, 2^L) + 0.5) / 2^L`) was implemented and rejected after the two-level test: the coarse acoustic transient was extrapolated into the fine flux, and the next coarse step amplified it until `w` became NaN within about two hours. The fine region with the coarse state held at the start-of-interval snapshot is stable. The coarse state is therefore held fixed for the fine substeps and advanced once, at the end of its interval. The lag is one coarse step. It is the phase-2 coupling.

### Interface state and the imposed flux

The held snapshot is written into both time levels and into the diagnostics the acoustic pressure gradient reads (`rho_p`, `rtheta_p`, `exner`, `pressure_p`, `ru`, `rw`). It covers every coarser cell, not only the interface band, because the slow tendency loop is unmasked and would otherwise differentiate a stale neighbor. After the call, those cells are put back.

On an imposed edge the acoustic loop does not compute a pressure gradient. Each acoustic substep sets `ru_p = target − ru_save` and folds that into `ruAvg`, so the recovered `ruAvg = ru_save + mean(ru_p)` equals the target. Divergence damping is skipped on edges that are not the active level, including imposed edges, so damping cannot move `ru_p` off the target. Recovery updates `ru`, `u`, and `ruAvg` on imposed edges and skips other inactive edges.

On level L, `config_number_of_sub_steps` is multiplied by `2^L`, so the acoustic step stays the fine acoustic step. A coarse step at the nominal substep count (acoustic step `2^L` times longer) put the `w` error on the interface cells; that configuration was still growing at 2 hours and had become NaN within a few hours. Matching the acoustic step removed that growth. It also spends the acoustic-work part of the speedup: a coarse cell does as many acoustic updates per unit time as a fine cell. The tendency loop is still unmasked and the halo exchange is still global, so the two-level run is slower than all-fine. The profile, the acoustic-CFL check, and the interface-versus-interior `w` breakdown are in [PHASE2.md](PHASE2.md).

Inactive cells zero `rho_pp` and `rtheta_pp` on the first acoustic substep and then skip the column update. The active edge therefore sees a neighbor acoustic perturbation of zero and a neighbor exner from the held snapshot. That is a known source of interface acoustic noise. Do not hide it by adding diffusion. A fine-region run with the coarse state frozen and the coarse advance skipped stays at the all-fine `w`. The noise that remains is in the coarse advance.

Halo exchanges stay on the existing full-domain groups. Every rank participates. Per-level groups are phase 3.

### Where the mask goes

The mask is a cycle at the top of the loops that commit a level's state:

- `atm_rk_integration_setup` (so a second `atm_srk3` on the same tick does not copy an inactive time level 1 over the state just advanced)
- both acoustic edge loops and the acoustic cell loop
- `atm_divergence_damping_3d`
- `atm_recover_large_step_variables_work`
- `atm_rk_dynamics_substep_finish`, including the `rho_zz` time-level-1 restore

`atm_compute_dyn_tend` and `atm_compute_solve_diagnostics` are not masked. They recompute tendencies and diagnostics from whatever prognostics are in the pools. On the measured two-level run those full-domain tendency evaluations, plus global halo exchanges around a masked acoustic loop, outweigh the cells that skip the acoustic update. The 1.44× ceiling counts cell updates, not this split. The six-hour profile is in [PHASE2.md](PHASE2.md): with this loop structure the variable mesh is slower than all-fine.

The mask is inside `#ifndef MPAS_OPENACC`. An OpenACC build that sets `config_hts` aborts in setup. `config_hts = false` does not enter the masked bodies.

Rejected for phase 2:

- A second copy of `atm_srk3`.
- Separate fine and coarse domain blocks. That is phase 3, with the halo groups and the partition.

### Conservation

Dry air mass lives in `rho` over `areaCell * (zgrid(k+1) − zgrid(k))`. The acoustic step's horizontal density change is `−dt * div(ruAvg)`, because the horizontal part of `tend_rho` is frozen from the start-of-step `ru` and the acoustic update adds `div(ru_p)`, with `ruAvg = ru_save + mean(ru_p)`. Vertical fluxes move mass inside the column. Matching `ruAvg` on an interface edge therefore matches the two cells' horizontal mass budgets. The Jablonowski–Williamson mesh has flat `zgrid`, so that is the dry-mass integral.

The coarse acoustic step is asked to land on the accumulated fine flux. If it does not, `atm_hts_correct_mass` adds the missing divergence to the coarse column's `rho_zz` (and the matching `theta_m` / pressure diagnostics) and sets `ru`, `u`, and `ruAvg` to the target. The sign is the same as the acoustic divergence: `drho_zz = −edgesOnCell_sign * dvEdge * invAreaCell * dt * (target − ruAvg)`. The first three coarse steps, and any step whose residual exceeds `1e-3`, log the max `|ruAvg − target|`.

Checks:

- On a uniform mesh (one level) the accumulator is idle and the result must match the baseline to within a small multiple of single-precision roundoff. Bit identity is not required; the mask branches can change reduction order.
- On a variable mesh, the relative mass drift of HTS over 10 days must be in the same range as the all-fine baseline recorded in BASELINE.md. A drift that grows with the number of interface edges is a failed flux match, not a CFL issue.
- Accumulate the interface flux in double precision. Single-precision running sums over `2^L * N_acoustic` terms are a plausible way to lose the conservation the scheme is supposed to protect.

Scalar mass (when `config_scalar_advection` is true) uses the same `ruAvg`. Phase 2 runs with scalar advection off, matching the JW namelist. Turning it on is the test that the accumulated flux, not a separately diagnosed velocity, is what transport sees.

### Physics

Phase 2 leaves physics off (`config_physics_suite = 'none'`), which is the JW configuration.

The coupling that landed with phase 3, and the phase 3b column mask:

- `physics_timetracker` still runs once per fine clock tick, so radiation and convection alarms stay on the fine clock. On one level and with HTS off, `physics_driver` is still that one full-block call.
- On a multi-level step the driver moves inside the level loop. It sees only the owned cells of the active level (`column_mask_on`, contiguous runs of that mask). `dt_dyn` for that call is the level dt. A radiation or convection alarm that rings on a tick when a coarser level does not step stays latched until that level steps, so the coarse columns are not dropped. One level never sets the mask, and the driver call is the original its:ite.
- `physics_get_tend` still runs at the start of every `atm_srk3`. `hts_restrict_physics_tendencies` then zeros `tend_ru_physics`, `tend_rtheta_physics`, and `tend_rho_physics` off the active level. All-fine and a single level do not enter the zeroing.
- Microphysics still sits at the bottom of `atm_srk3` and uses that level's dt (Thompson keeps its 90 s split). The copy loops skip inactive cells. Kessler skips inside its column loops. WSM6, Thompson, and Thompson-aerosols are called through `column_mask_next_segment`, so with the mask off each kernel still sees the original `its:ite` in one call. Radar reflectivity uses the same segments. The 1D ocean mixed layer, the surface diagnostics, the skin SST, and the deep-soil temperature read `column_is_active`. Skin SST and deep-soil temperature move onto the level step only when more than one level is in use; both namelist flags default to false. Ozone, the monthly albedo, and the SST file update stay on the fine-clock alarm for every owned column. The mask is the owned active list and not the dynamics work ring. One level and all-fine leave `column_mask_on` false.

Measurements are in [PHASE3.md](PHASE3.md) and [PHASE3B.md](PHASE3B.md).

### MPI and halos

Phase 2 keeps the existing unweighted `graph.info.part.*` files and the global exchanges described above. That is load-imbalanced: on `x5.30210` about 55% of cells are level 0, but a rank can easily own a different fraction.

Phase 3 partition: METIS vertex weights proportional to the number of substeps (`weight = 2^(L_max - L)`). `gpmetis` is not bundled; a small preprocessing script writes a weighted graph and the model keeps reading a plain partition file. No framework change is required for the first weighted partition. Per-level halo groups come after the weights, and only if profiles show the full-domain acoustic exchange dominating.

Regional LBCs and HTS together are out of scope until both work separately. The LBC interpolant should be reused, not reimplemented.

## 3. Validation

All comparisons use the same compiler, flags, core count, and namelist as [BASELINE.md](BASELINE.md) except for the HTS switches. The baseline is single precision, which is the MPAS default. One double-precision (`PRECISION=double`) pair is required for the conservation test, because single-precision drift can hide a flux bug. It is not required for the speed test.

Metrics, uniform mesh (`x1.40962`), HTS forced to one level:

- Relative L2 and Linf of `u`, `w`, `theta`, `rho` against the baseline history at day 1 and day 10.
- Target: relative L2 below `1e-5` (single) and below `1e-10` (double). Larger than that means the mask changed the numerics.
- Relative dry-mass drift no worse than the baseline drift in BASELINE.md.

Metrics, variable mesh (`x5.30210`), HTS versus an all-fine reference at `dt = 180 s`:

- Same L2 / Linf at day 10, over all cells and again over interface cells only (cells with a neighbor at a different level).
- Dry-mass drift compared with the reference.
- Interface noise: rms of `w` in the interface band versus rms of `w` in the fine interior (cells whose neighbors are all the same level). A persistent spike at the interface that is absent from the reference is a failed interpolant or a failed flux match, even if the global L2 looks acceptable.
- Surface-pressure departure from the initial state, rms over the sphere, at each day. The JW wave should grow on the same schedule as the reference. Day 1 is too early to see the wave; day 10 is the comparison that matters.

Speed, same binary family, 4 ranks, `x5.30210`, 10 days:

- Wall time of the all-fine reference versus HTS. On this mesh the ceiling estimated from the level counts is about 1.44×. A result near that, with the accuracy targets met, is the phase-2 success condition. A much larger claimed speedup would mean the fine region was not actually subcycled at `dt`.
- A second mesh with a small fine patch is the phase-4 speed case. It is not this mesh.

## 4. Risks

- The vertically implicit acoustic solve is per column, but the horizontal pressure gradient reads neighbors. A mask that skips the neighbor update without refreshing the interpolated halo will either deadlock, drift, or radiate acoustic noise from the interface.
- The coarse boundary seen by the fine step is the start-of-interval state, held fixed, formally lower order than RK2/RK3. Extrapolating the previous increment was unstable on `x5.30210`. Expect a local error at the interface. Do not "fix" it by adding diffusion.
- Single-precision flux accumulation.
- `config_split_dynamics_transport = true` (the registry default, off in the JW namelist) accumulates transport fluxes across several dynamics substeps. HTS has to define that average on interface edges or the scalar and the continuity equation will use different fluxes. Phase 2 keeps the split off until the unsplit case is clean.
- Physics packages that key off the step index rather than the clock will fire at the wrong rate on coarse columns.
- OpenACC branches inside `atm_srk3` must keep compiling. Phase 2 is a CPU path; the mask has to be a no-op the accelerator build can ignore, not a second solver.
- Load imbalance will make a correct HTS look slow on equal-cell partitions. Do not tune the numerics to hide that. Fix it with weights in phase 3.
- The 5× tutorial mesh does not have a 4× or 8× interface. Binning to powers of two is slightly conservative and must be stated when quoting a CFL.

## 5. Meshes and real-data inputs

Scripts and URLs are in `docs/hts/scripts/`. Meshes are downloaded, not committed.

### Idealized meshes used for the baseline

| mesh | cells | spacing | file | source | server mtime | size |
| --- | ---: | --- | --- | --- | --- | ---: |
| `x1.40962` | 40962 | 120 km quasi-uniform (`dcEdge` about 97–127 km on the unit sphere, scaled by 6371229 m) | `x1.40962.grid.nc` inside the JW tarball | https://www2.mmm.ucar.edu/projects/mpas/test_cases/v7.0/jw_baroclinic_wave.tar.gz | 2019-06-08 | tarball 26 MB; grid 58 MB |
| `x1.10242` | 10242 | 240 km, `nominalMinDc * R = 240000` m exactly | `x1.10242.grid.nc` | https://www2.mmm.ucar.edu/projects/mpas/atmosphere_meshes/x1.10242.tar.gz | 2023-06-11 | tarball 6.1 MB |
| `x5.30210` | 30210 | 240 km to 48 km, ratio 5, refinement centered at (0°N, 0°E). `nominalMinDc * R = 48000` m. `dcEdge` spans about 39–287 km. | `x5.30210.grid.nc` | https://www2.mmm.ucar.edu/projects/mpas/atmosphere_meshes/x5.30210.tar.gz | 2023-06-13 | tarball 22 MB; grid 34 MB |

`x1.10242` was downloaded and checked, and is not the baseline integration. The baseline uniform case is `x1.40962` because that is the mesh shipped with the official JW namelist (`config_dt = 450` s, `config_len_disp = 120000` m, 26 levels, case 2, physics suite `none`).

`x5.30210` is not listed on the mesh table at https://www2.mmm.ucar.edu/projects/mpas/site/downloads/meshes.html. It is in the same directory, returns HTTP 200, and is the mesh the current MPAS tutorial uses for its 240–48 km exercise (tutorial text: 30210 cells, factor-of-five ellipse centered at 0°N, 0°E). The smallest variable mesh that *is* on the public table is `x4.163842` (92–25 km, 163842 cells, 135 MB). That cell count is the same as a global 60 km mesh and is a tight fit for this VM (about 16 GB RAM, 4 cores). `x5.30210` is the variable mesh the baseline actually integrates.

The published table's other variable meshes (46–12 km and finer, hundreds of thousands to millions of cells) are too large for this VM.

### Real-data files that can be downloaded without WPS

NCAR's real-data example page (https://mpas-dev.github.io/atmosphere/real_data.html and https://www2.mmm.ucar.edu/projects/mpas/site/downloads/real_data.html) does **not** ship a pre-built `init.nc`. `https://www2.mmm.ucar.edu/projects/mpas/real/v7.0/x1.40962.init.nc` returns 404. What it does ship, all for the uniform 120 km mesh:

| file | URL | server mtime | size | what it is |
| --- | --- | --- | --- | --- |
| `x1.40962.static.nc` | https://www2.mmm.ucar.edu/projects/mpas/real/v7.0/x1.40962.static.nc | 2019-06-08 | 76 MB | pre-interpolated static fields (terrain, land use, …) for `x1.40962` only |
| `cfsr.2010102300.tar.gz` | https://www2.mmm.ucar.edu/projects/mpas/real/cfsr.2010102300.tar.gz | 2014-05-28 | 81 MB | ungrib intermediate `CFSR:2010-10-23_00` |
| `namelist_files.tar.gz` | https://www2.mmm.ucar.edu/projects/mpas/real/v7.0/namelist_files.tar.gz | 2019-06-08 | 2.5 KB | v7 namelists and streams. `config_static_interp = false`, `config_met_interp = true`, `config_physics_suite = 'mesoscale_reference'`, `dt = 720` s, 55 levels, 5-day run |

Quasi-uniform static tarballs also exist for every resolution on the mesh table (the 120 km static tarball is listed as 16 MB compressed). Variable-resolution meshes on that table are grid-only. `x5.30210_static.tar.gz` returns 404. The tutorial's precomputed `SouthAmerica.static.nc` is on NCAR Glade (`/glade/campaign/mmm/wmr/mpas_tutorial/precomputed/`), not on a public URL.

The geographic dataset needed to *create* a static file is public:

- https://www2.mmm.ucar.edu/projects/mpas/mpas_static.tar.bz2
- 2.2 GB, Last-Modified Thu, 27 Jun 2024 20:17:25 GMT

It was not downloaded in phase 1.

### Cheapest path to a variable-resolution real-data case

There is no public pre-built `init.nc` for any mesh, and no public pre-built `static.nc` for a variable mesh.

1. Keep `x5.30210` (30210 cells). Do not start from `x4.163842` on this VM.
2. Download `mpas_static.tar.bz2` (2.2 GB) once.
3. Run `init_atmosphere_model` serially with `config_init_case = 7`, `config_static_interp = true`, GWDO static on, met interpolation off. This is the slow serial step. The users' guide warns that the GWDO stage wants on the order of 4 GB for the 30-arc-second terrain field, which this VM can spare if nothing else is running.
4. Run `init_atmosphere_model` again with the already-downloaded `CFSR:2010-10-23_00` intermediate (`config_met_interp = true`, static stages off). That is `init.nc` without WPS and without ungrib.
5. Integrate with the sample physics suite only after phase 3. A one-hour dry dynamics check on that `init.nc` can happen earlier and would not test physics coupling.

The uniform 120 km real-data case can skip step 2 and 3: the 76 MB `static.nc` and the CFSR intermediate are enough to build `init.nc`. That run was not done in phase 1; the dry JW baselines were. The v7 real-data namelist still sets `config_epssm`, which v8.4.2 ignores (it warns and the height-dependent `config_epssm_minimum` / `config_epssm_maximum` are what the solver uses). Any real-data run on v8.4.2 needs that namelist updated the same way the JW baseline was.

## 6. Phase 2 plan

Phase 2 is dry dynamics only. Physics stays `none`. Scalar transport stays off until the unsplit dynamics conserve mass.

Order of work:

1. Registry switches: `config_hts` (logical, default false), `config_hts_max_level` (integer, default 1). Default false is a no-op so phase-1 namelists do not change. Done.
2. Level arrays filled in `mpas_atm_dynamics_init`, with a log line of the cell count per level. Done.
3. Mask the committing loops. Uniform mesh, one level, comparison against BASELINE.md. Done: day 1 and day 10 are bit-identical to the all-fine history of the same binary and to the phase-1 baseline.
4. Two-level flux accumulation and the held snapshot above. Test on `x5.30210` with `config_hts_max_level = 1` (230 km cells run at `2 dt`). Compare mass, interface `w`, and fields against an all-fine run of the same binary. Done for 6 hours and 1 day. The profile and the `w` breakdown are in [PHASE2.md](PHASE2.md). Wall time is about 1.65× the all-fine run. A third level was not started: the same acoustic-step rule and the same full-domain tendency and halo calls would add cost.
5. Scalar transport and `config_split_dynamics_transport` stay off until a design exists that can be faster than all-fine. The uniform regression for that configuration was not run.

The driver lives in `mpas_atm_time_integration.F` (`atm_hts_advance`) so it can call `atm_srk3` without a module cycle. Level assignment, the held snapshot, the flux accumulator, and the mass patch are in `mpas_atm_hts.F`. Halo-group definitions are unchanged.

Out of phase 2 as first written: physics coupling, weighted partitions, per-level halo groups, regional LBCs, GPU. Phase 2b did the weighted partition and the index lists. It left the halo groups as they are. The measurements are in [PHASE2B.md](PHASE2B.md).

## 7. Phase 2b

Phase 2b keeps the held state and the matched acoustic step. It changes where the work is done.

- Tendency, dissipation, acoustic, and damping loops walk a per-level index list (`src/core_atmosphere/dynamics/hts_loops.h`). The work list is the active level plus two cell rings. Diagnostics walk a separate list at four rings, wide enough for the APVM and kinetic-energy stencil the tendency reads, including the stick-out edges of those vertices. A diagnostic that stops on the two-ring work list leaves that fringe stale, and the next level's first RK tendency consumes it. On the six-hour dry `sp8.24755` history the four-ring list keeps the all-fine errors at the same three significant figures (`u` relative L2 7.89×10⁻⁴). Against the previous full-domain diagnostic binary the same frame is `u` relative L2 1.21×10⁻⁵.
- Halo exchanges stay on the existing full-domain groups. Every rank has to enter every collective. After the partition below, every rank owns every level, so a per-level group would still include every rank. Fusing the two acoustic messages and packing only the active-edge halo were both measured on the 6-hour dry run and removed; the numbers are in [FINAL.md](FINAL.md). The seconds that remain after that are the coarse-state hold and the fine acoustic step, not the message count.
- The partition is a METIS multi-constraint graph, one constraint per level. `docs/hts/scripts/hts_partition.py` writes it. MPAS reads the plain `*.part.N` file through `config_block_decomp_file_prefix`.

```
docs/hts/scripts/hts_partition.py GRID.nc GRAPH.info NPARTS --max-level 1
```

The script writes `GRAPH.info.htsw` with header `n m 10 ncon` and, on each line, `ncon` positive weights (1000 on the cell's own level, 1 on the others) followed by the original 1-based neighbors. It then runs `gpmetis -ptype=kway GRAPH.info.htsw NPARTS`. Copies of the 4-rank files used for the phase 2b timings are in `docs/hts/partitions/`. The weighted graphs are regenerable and are not committed.

A longer acoustic step on long coarse edges was not added. The interface and the short transition edges do not have the CFL room, and the earlier 60 s acoustic step on `x5.30210` put the `w` peak on the interface.
