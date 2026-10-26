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

1. Snapshot any level whose interval is opening (`mod(tick, 2^L) = 0`). The snapshot is the time-level-1 state at the start of that level's step. The previous-step increment is stored.
2. On the diagnostic ring, rewrite coarser `rho`, `theta`, and `w` from the last realized snapshot and rebuild the pressure diagnostics from that `rho` and `theta`. The snapshot taken at the open is already the end of the previous interval, so `snapshot + alpha * previous_increment` leads that state. That lead is not applied: a led `rho`/`theta` beside a snapshot `w` and a snapshot face flux is the boundary the fine column advects. An interface edge whose edge level is the active level is locked to that snapshot wind, `ru` is that wind times the two cells' densities, and the acoustic step holds `ru` there. The face is not moved to cancel the column divergence. Coarser levels are not rewritten. Before the fine step, scalar time level 1 is copied onto time level 2 for every coarser column, because the flux reads time level 2 and the RK setup copies scalars only on the active level. After the tick the same copy is repeated for columns whose level did not close, so the core shift cannot publish the previous scalar buffer. The scalar commit is masked to the active level. The first interval has no previous increment, so the edge stays a free fine prognostic until the coarse level has stepped once. Interior columns stay on the snapshot. A level whose interior already matches the snapshot is not rewritten. Reconstruct the horizontal wind, then call `atm_srk3` masked to the active level.
3. If the hold was rewritten, put those coarser cells and edges back on the snapshot before the next coarser RK. The fine step does not commit into a held level.
4. Add `ruAvg * dt` on level-0 interface edges into a double-precision accumulator.
5. For each coarser level that is due, from fine to coarse: copy the finer time-level-2 state back onto time level 1 and onto the saved acoustic base state; set the imposed `ruAvg` target to (accumulated flux) / `dt_L` on interface edges whose coarser cell is this level; hold levels coarser than L at their snapshot; call `atm_srk3` with `dt_L`; restore those coarser levels; patch any flux residual (below); accumulate this level's own outgoing interface flux; clear the accumulator entries this level consumed.

A rejected alternative was a predictor: integrate the coarse region once to get the interval endpoint, interpolate, then integrate it again with the fine region. That integrates the coarse region twice and spends the speedup. Linear extrapolation of the full previous increment (`alpha = (mod(tick, 2^L) + 0.5) / 2^L` on `rho`, `theta`, `u`, and `w`) was implemented and rejected after the two-level test: the coarse acoustic transient was extrapolated into the fine flux, and the next coarse step amplified it until `w` became NaN within about two hours. Holding every coarse field at the snapshot is stable for dry dynamics and is not stable for land. The fine interface edge then integrates `tend_u_euler` against a neighbor whose pressure does not move while the fine Noah/YSU column does, and on the eds land VR that edge wind ran away by about two hours. Job 199900, the fine-dt PBL/Noah subcycle, reached the same cell sooner: cell 22181 (31.84°N, 0.14°E, land, `meshDensity` 0.610) saw `|u|` jump to about ±269 m/s on the level-3 step at 00:39:30, and `u` was NaN at 00:45. That subcycle is removed. `bb575a1b2` advanced only the open-interface parent's acoustic `rho_pp` and `rtheta_pp`. Eds job 199953, rebuild 199950, had `u` about ±7.2×10⁵ m/s and `w` about −1×10⁴ / +8×10⁴ m/s at 00:17:30 and NaNed `u` at 00:18. `cf1dc652d` ramped only held-ring `rho` and `theta` and left `u` and `w` pinned. Eds job 199977, rebuild 199974, GWDO off, stayed near ±88 m/s through 00:09, reached `u` about ±330–532 m/s and `w` about ±17–26 m/s at 00:11:30, and NaNed `u` at 00:12:30 (`w` still finite; HTS wall 57.63 s, exit 11). CPAS interpolates a full parent state. Neither of those partial boundaries is that contract, and both are removed. This tree holds every coarser field on the snapshot until that level steps. Job 199996 on the pure hold `20b4f3a34` reproduced that floor: all-fine finished the 2 h land pair, and HTS NaNed `w` at 01:57:30 after a `u` spike, the same clock as `d3fe6fb75`. `45dfe9329` is that hold plus the 199977 writeup. It is not the `615836256` thermo-only ramp. This tree keeps the close-step of that hold (same-interval imposed flux, one physics call at `dt_level`, acoustic parents at zero perturbation). Once a level has completed an interval, the diagnostic ring blends `rho` and `theta` toward that realized increment. Ring `u` and `w`, and an interface edge the active level would otherwise integrate, stay on the last realized snapshot. The snapshot at the open is already the previous endpoint, so adding `alpha * delta` there leads the wind by up to one coarse step and the impose commits it. Scalars are not rewritten by that snapshot. Time level 2 is copied from time level 1 on columns that are not advancing, so the flux and the following shift both see the committed moisture, and the scalar commit is masked to the active level. Integrating that edge from the parent's own pressure gradient NaNed dry `w` at 00:42 and is not in this tree. Job 200072 on `95ddadd85` (the lead) finished the 2 h land pair with no NaN and max `|u|` about 160–171 m/s. A local dry 6 h without that lead finished with 0 critical errors and no NaN: at 01:57 `|u|` is −35.0813 / 35.0435 m/s and `|w|` about −0.016 / 0.022 m/s, and at 6 h `|u|` is −35.026 / 35.018 m/s. Integration on that sitting was 63.35 s; quote the published 62.95 s. Eds job 200101 (`efe084b70`, `hts_land_snap2h`) reran that no-lead lock. Survival was unchanged and the lead was not the cause: max `|w|` was already 18.8 / 18.6 / 17.5 / 14.5 m/s on the first four frames (all-fine at most 0.59), `w` relative L2 was 3.86 / 6.48 / 6.87 / 7.69, and max `|u|` stayed with the jet through 01:00 (87.7 / 88.6 / 89.1) then 111.3 at 01:30 and 171.1 at 02:00. At 2 h, `u` rms / relL2 1.27 / 0.098, `theta` rms 0.89 K, `p` rms 39 Pa, land `t2m` / skin 1.24 / 1.61 K, dry-mass drift −2.6×10⁻⁸ against all-fine 1.2×10⁻¹¹, time integration 286.11 s versus 352.62 s (1.23×). Job 200144 (`1e3255beb`) then moved the locked level-0 face by up to 20 m/s so that column's horizontal divergence was zero. That is a regression against `efe084b70`: `w` NaN at 01:29:30, `u` NaN at 01:30:00, SIGSEGV at 01:30:30 (182/240 steps). At 00:30 max `|w|` was 18.80 m/s, the same cell and value as job 200101 (18.84 m/s, cell 22181, k=33 from the bottom, 31.84°N, 0.14°E). That cell's six neighbors are level 0; the nearest level-1 cell is 133 km away, so the face move does not act on the column where `|w|` peaks. The vertical implicit coefficients are column-local. Forcing every dry correction (about 0.05 m/s) NaNed `u` near 04:35, and matching coarser levels wrote `|u|` about 2700 m/s at 00:02:30. The face move is not in this tree. The ring `rho`/`theta` lead is not either: those fields now stay on the snapshot with `w` and the locked wind, and the face mass flux is that wind times the two cells' densities. Job 200185 (`4c880bf40`, `hts_land_snapth2h`) reran that snapshot ring. It finished 240/240 with no NaN and 0 critical errors, and the numbers match job 200101: max `|w|` 0.328/18.84/18.89/17.32/14.53 against all-fine 0.328/0.559/0.434/0.502/0.594, `w` relL2 0/3.86/6.49/6.85/7.68, max `|u|` 87.73/88.60/89.06/111.35/171.10, frame-4 `u` relL2 0.098, `theta` and `rho` relL2 0.0022, `u` rms 1.275, `theta` rms 0.893 K, `p` rms 39.2 Pa, land `t2m`/`skintemp` 1.244/1.610 K, dry-mass drift −2.6×10⁻⁸, time integration 286.25 s versus 351.7 s (1.23×). The 00:30 max `|w|` is again cell 22181, k=33, 31.84°N 0.14°E. Holding ring `rho` and `theta` on the snapshot did not move it. The `|u|` overshoot from 01:30 is not that interface: edges 1321 and 2319 have both cells at level 3, near 11–12°N, 100.5–100.9°W, k=41–45. The 00:00 and 00:30 `|u|` maxima are the level-3 southern-ocean jet at k=52. Across 200101, 200144, and 200185, `rho`, `theta`, `u`, and `w` are written into both time levels for a column that does not step, and the scalars are not. The shift then publishes the previous scalar buffer while `theta_m` stays on the snapshot, and the fine flux reads that buffer. This tree copies scalar time level 1 onto time level 2 before the fine step and again, for columns that did not close, before the shift. A local dry 6 h of that copy finished with exit 0, 0 critical errors, and no NaN (`|u|` −35.0255 / 35.0183 and `|w|` −0.0273 / 0.0345 at 05:59:30; run extrema −35.0828 / 35.0612 and −0.0475 / 0.0429; integration 62.59 s; quote the published 62.95 s). Dry scalar advection is off, so this sitting does not show that the land `w` moved. Job 200237 (`350e1fb3f`) passed the 2 h land pair: 1.61×, max `|u|` 89.1 matching all-fine, final `w` relL2 0.433. Job 200248 ran that binary for 6 h (720/720, no NaN, critical 0, 1.55×, time integration 723.6 s versus 1124.6 s) and grew one interior level-3 column. Cell 3177 (16.16°S, 137.11°E, area about 21900 km², spacing 8.85) and its six neighbors are level 3. The nearest other level is 97 cells away, so this is not an interface. From 02:30, `w` there doubled about every 30 min and the maximum rose from k=22 to k=36, while all-fine stayed near 0.05 m/s. A level-3 dynamics step is 240 s. Kessler removes the whole supersaturation in that one call, and the Smagorinsky cap is proportional to `1/dt`, eight times tighter than on the all-fine step. Job 201055 (`cb274ed09`) ran a `config_dt` dynamics step on every level. Accuracy at 2 h matched job 200237 (max `|u|` 89.3, `w` relL2 0.427, no NaN) and the speedup was 0.756× (time integration 492.1 s versus 372.0 s): 960 RK steps versus 240, scalar advection 123 s versus 32 s. The coarsest level is now four coupled steps, so level 3 is 60 s. The acoustic scale is `dt_sub/config_dt`, which keeps `dts` at 5 s, and `ruAvg` is the mean of those steps. Levels below the coarsest stay one long step. PBL and Noah stay one call at the level dt. On a moist coarse step the Smagorinsky cap uses the fine dt. Dry does not enter the split and does not change the cap. A local dry 6 h of this split finished with exit 0, 0 critical errors, and no NaN (`|u|` −35.0255 / 35.0183 and `|w|` −0.0273 / 0.0345 at 05:59:30; run extrema −35.0828 / 35.0612 and −0.0475 / 0.0429; integration 63.80 s; quote the published 62.95 s). Eds has not run it. Expected land speedup is about 1.3–1.4×. The 8 h pair is not accepted yet.

### Interface state and the imposed flux

Interior coarse columns stay on the snapshot in both time levels. During the fine step the diagnostic ring is rewritten to that same snapshot (`rho`, `theta`, `w`, and the pressure diagnostics), and a held interface edge holds the last realized coarse wind with `ru` from the two cells' densities. The restore writes that snapshot again before the time-level shift. The coarse RK therefore still starts from the snapshot.

On an imposed edge the acoustic loop does not compute a pressure gradient. Each acoustic substep sets `ru_p = target − ru_save` and folds that into `ruAvg`, so the recovered `ruAvg = ru_save + mean(ru_p)` equals the target. The same replacement is used on an interface edge the active level is holding at the last realized coarse wind. Adding the previous `u` increment on that edge is a lead: the fine-step flux is the overshoot, and the close-step impose commits it. Job 200072 (`95ddadd85`) survived the 2 h land pair with that lead (no NaN, all-fine wall 355.5 s, HTS wall 289.1 s) and finished with max `|u|` about 160–171 m/s against an all-fine jet near 89 m/s, `u` rms 1.27 m/s, `theta` rms 0.89 K, land `t2m`/`skintemp` rms 1.25 / 1.61 K, and `w` relative L2 about 7.6. Job 200101 removed the wind lead and the numbers barely moved, including max `|w|` about 19 m/s already at 00:30. Job 200144's level-0 divergence match (up to 20 m/s on the locked face) left that `|w|` at 18.80 m/s and NaNed the run at 01:29. The acoustic step holds the snapshot flux (`ru_p = 0`). Parent acoustic pressure stays zero. Divergence damping is skipped on those held edges and on edges that are not the active level, including imposed edges, so damping cannot move `ru_p` off the target. Recovery updates `ru`, `u`, and `ruAvg` on imposed edges and on held interface edges.

On level L, `config_number_of_sub_steps` is multiplied by `2^L`, so the acoustic step stays the fine acoustic step. A coarse step at the nominal substep count (acoustic step `2^L` times longer) put the `w` error on the interface cells; that configuration was still growing at 2 hours and had become NaN within a few hours. Matching the acoustic step removed that growth. It also spends the acoustic-work part of the speedup: a coarse cell does as many acoustic updates per unit time as a fine cell. The tendency loop is still unmasked and the halo exchange is still global, so the two-level run is slower than all-fine. The profile, the acoustic-CFL check, and the interface-versus-interior `w` breakdown are in [PHASE2.md](PHASE2.md).

Interior coarse cells, including open-interface parents, zero `rho_pp` and `rtheta_pp` on the first acoustic substep and then skip the column update. The active edge therefore sees a parent whose acoustic perturbation is zero. On the diagnostic ring the exner is rebuilt from the snapshot `rho` and `theta`; interior exner stays the snapshot. Eds job 199953 showed that advancing the parent's acoustic pressure alone drives `w` to about 10⁴ m/s and moves the land-VR NaN from about 01:57 to 00:18. Eds job 199977 (rebuild 199974, GWDO off) showed that ramping only `rho` and `theta` moves it to 00:12:30, with `w` still finite on the first NaN print and HTS wall 57.63 s. Do not hide a remaining interface signal by adding diffusion.

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
- `physics_get_tend` still runs at the start of every `atm_srk3`. `hts_restrict_physics_tendencies` then zeros `tend_ru_physics`, `tend_rtheta_physics`, `tend_rho_physics`, and `scalars_tend` off the active level. All-fine and a single level do not enter the zeroing. The monotonic final RK stage also skips that source update on an idle column, because it writes `scalars_tend` into time level 1.
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
- The per-substep layer-1 exchanges of `rho_pp` and `rtheta_pp` use a compact list of halo cells the active level touches. A mismatch in the all-to-all count check keeps the full group. Every other halo stays on the existing full-domain group, and every rank still enters every collective. After the partition below, every rank owns every level, so a per-level group would still include every rank. Fusing the two acoustic messages, and an earlier ownership-filtered pack, were measured and removed; the numbers are in [FINAL.md](FINAL.md). The coarse-state hold skips a level that already matches its snapshot and reconstructs wind only on cells whose edges changed. What remains after the compact layer-1 list is the fine acoustic step and the wider halo groups.
- The partition is a METIS multi-constraint graph, one constraint per level. `docs/hts/scripts/hts_partition.py` writes it. MPAS reads the plain `*.part.N` file through `config_block_decomp_file_prefix`.

```
docs/hts/scripts/hts_partition.py GRID.nc GRAPH.info NPARTS --max-level 1
```

The script writes `GRAPH.info.htsw` with header `n m 10 ncon` and, on each line, `ncon` positive weights (1000 on the cell's own level, 1 on the others) followed by the original 1-based neighbors. It then runs `gpmetis -ptype=kway GRAPH.info.htsw NPARTS`. Copies of the 4-rank files used for the phase 2b timings are in `docs/hts/partitions/`. The weighted graphs are regenerable and are not committed.

A longer acoustic step on long coarse edges was not added. The interface and the short transition edges do not have the CFL room, and the earlier 60 s acoustic step on `x5.30210` put the `w` peak on the interface.

Moist runs cap the active level at 2. At `config_dt = 30` s the coarse dynamics, Kessler, and surface step is 120 s. `config_hts_pbl_substeps = 0` is the config_dt surface fold and is not the default. The runs through job 201388 still used one surface call of the level dt. Dry runs keep `config_hts_max_level`. A single 240 s level-3 step (`350e1fb3f`) matched all-fine at 2 h and then grew `w` at interior cell 3177 (job 200248). Eight 30 s dynamics steps on every level, with PBL still once per 240 s (`cb274ed09`, job 201055), moved that `w` closer at 2 h and cost 0.756×. Four 60 s dynamics steps on the coarsest level only, PBL still once per 240 s (`0fa07c436`), cost 1.17× on the 2 h pair and failed the 8 h pair (job 201240): interior edge 64121, 63 cells from any other level, reached `|u|` 221 m/s by 03:12 on a two-step oscillation, and cell 3177 stayed on the 240 s trajectory. The acoustic lock and the mass correction only touch interface edges, so they do not set that wind. A local Kessler dump with PBL and Noah off moves the same edge by about 2×10⁻⁴ m/s per 60 s substep. The land failure is the restarted dynamics under one long surface call. YSU stores `(u_mixed - u) / (2 * dt_pbl)` and does not change `u`. The split calls that once at 240 s, then runs four 60 s SRK3s. Each SRK3 multiplies the stored rate by its own dt, so the integrated surface impulse is `0.5 * (u_mixed - u)`, the same total one 240 s SRK3 applies. The paths diverge because every SRK3 zeros `rho_pp` and `rtheta_pp` on its first acoustic step and rebuilds the pressure response, while the surface rate stays the one computed for the whole 240 s. The promote between substeps copies `u`, `w`, `theta_m`, `rho_zz`, and the scalars, and does not carry that acoustic pressure. Edge 64121 is interior level 3, between land cell 7793 and water cell 8033, so the interface lock and the mass correction never set it. The printed global min alternating every 240 s is that physics cadence: the wind one group leaves is the wind the next YSU/Noah call sees. Smagorinsky is not the driver. On the split `invDt = 2 / config_dt`, twice the coupled level-3 cap, so the deformation ceiling is looser, not tighter. Cell 3177 stayed on the 240 s `w` path because the surface call both the split and `350e1fb3f` share is still one 240 s YSU/Noah step. Calling PBL on each 60 s substep with `dt_pbl` equal to that substep would retune the rate, add about three extra coarse surface calls on top of a path that was already 1.17×, and is the same class of split the 120 s cap removes. Eight coupled 30 s steps on the coarse level would match the dt that moved cell 3177 at 2 h and, on the split2h timers, would spend the speedup. The 120 s cap is the run eds measured next.

Eds job 201325 (2 h) passed at 1.48×: max `|u|` 88.93/88.99, `w` relative L2 0.428, edge 64121 k2 −11.9 versus −13.4, no alternation, dry-mass drift about 1e-9. Job 201338 (8 h, level-2 cap, one 120 s surface call) was 1.40×, all 17 frames, no NaN, critical 0, and edge 64121 stayed near all-fine for the whole 8 h. Cell 3177 did not. It is an interior level-2 water column (every neighbor level 2, nearest finer cell 99 hops). `w` grew from 0.076 at 03:00 to 4.19 at 07:00 while all-fine stayed 0.04–0.11, then fell to 2.93 by 08:00. The same column at level 3 (240 s) started that growth near 02:30. The upper-level `u` burst (edge 41342, and edge 26872 at k40) rose and fell with that `w`. Low-level `|u|` stayed within a few m/s. This is not the 64121 surface-stress restart, and it is not an interface.

Job 201388 reran `2022d65c7` with `config_hts_max_level = 1` (885 cells at 30 s, 23870 at 60 s) for 8 h at 1.22×. Cell 3177 stayed within 0.01 of all-fine through 04:00, then `w` was 0.086 at 04:30, 0.294 at 06:00, 0.827 at 07:00, and 3.780 at 08:00. All-fine stayed 0.04–0.11. Once started, the column grew by ×1.5–2.3 per 30 min, the same rate as the 120 s run (onset 03:30) and the 240 s run (onset 02:30). Edge 41342 `|u|` was 26.3 at 07:00 and 58.5 at 08:00 versus all-fine 24.9. The largest edge difference at 08:00 was 48.7 at edge 26872 k35 (15.4°S 137.1°E). Edge 64121 stayed within 0.1–2 m/s, `k ≤ 5` max `|u|` within 1.4, global max `|u|` 89.45/89.51, final relative L2 `u` 0.085, `theta` 1.8×10⁻³, `qv` 0.094. Timers versus all-fine: tendencies 251/447, scalar advection 291/216, acoustic 261.7+46.3 halo/262.3, microphysics 26.7/46.7, physics driver 154.6/217.1.

Jobs 200248, 201338, and 201388 made the onset look like a function of the coarse step: about 02:30 at 240 s, 03:30 at 120 s, and 04:00–04:30 at 60 s, with the same ×1.5–2.3 per 30 min once it started. Two `c12f3c629` ablations at `config_hts_max_level = 1` removed that reading. Kessler at `config_dt` (`pbl_substeps = 1`, `kessler_substeps = 0`, log `Kessler uses 2 calls of 30 s`) and every process at `config_dt` (`micro_substeps = 0`, `pbl_substeps = 0`, log `coarse level 1: micro substeps 2, pbl substeps 2`) both tracked the 60 s baseline. Cell 3177 `w` at 04:00, 05:00, 06:00, 07:00, and 08:00 was 0.055 / 0.126 / 0.294 / 0.827 / 3.78 on the baseline, 0.057 / 0.127 / 0.293 / 0.828 / 3.82 with Kessler at 30 s, and 0.066 / 0.123 / 0.294 / 0.904 / 4.07 with every process at 30 s. All-fine was 0.047 at 04:00 and 0.112 at 08:00. Edge 41342 `|u|` was 58.5, 58.9, and 60.6 against all-fine 24.9. Global max `|u|` stayed within 0.05 m/s, `w` relative L2 was 1.39–1.40 at 08:00, edge 64121 was quiet, and dry-mass drift was −2.5×10⁻⁹. The largest edge difference grew from 6.8 to 33 m/s on the level-0 Sahara patch (28°N, 1–3°E, k30–40). About 305 edges at the level-0/1 interface near 27–28°N had k0 `|Δu|` > 5. Run B does twice the dynamics-tendency and scalar-advection calls of all-fine and was 0.83×. The Kessler-only run was 1.155×.

`config_hts_debug_force_path` with every cell at level 0 (job 201473) was bit-identical to all-fine for 8 h, about 2% slower. The 3177 growth turns on when a level ≥ 1 exists. It does not depend on the coarse dt or on the Kessler, PBL, or dynamics cadence. Run B takes two 30 s dynamics steps and two 30 s surface calls through the level-1 path and still matches the single 60 s step. What those runs share, and the level-0 force path does not exercise on a coarse column, is the hold and restore, the time-level promote between substeps, the `ruAvg` mean and the single mass correction, the per-level scalar commit and the idle-scalar copy, the latched radiation alarm (`hts_note_physics_alarms`), diagnostics that physics reads on a coarse column, and the acoustic `*_pp` restart on inactive neighbors. Cell 3177 (16°S, 137°E) is about 99 hops from a finer cell. A 20 m/s signal from the interface does not reach it in 4 h, and the visible growth starts near 13:30 local. The mass correction and the acoustic lock only touch interface edges. They stay in the switch list because run B still executes them. The 8db2f20c9 bisect below turned each of those off in turn. The idle copy is the one that moves the column. The open-tick shift after the switch list is the mechanism, and the tip copies the prognostics with the scalars.

The diabatic rate is not long by 2^L. `hts_push_microp_dt` sets `dt_dyn` to the SRK3 dt before Kessler, `post_microphysics` divides the theta change by that dt, and the RK3 subtract removes the same heat from `rtheta_p`. `w` is not subtracted. That is the same one-step lag all-fine has. Extra Kessler calls inside one dynamics step do not move the onset, and the 30 s Kessler ablation above did not either.

The fold is a different failure. `config_hts_pbl_substeps = 0` with one dynamics step of the level dt (jobs 201458 and 201494) warmed polar land skin by 3–7 K, raised land `t2m` rms from 0.23 to 0.43 at 2 h (1.348×), and blew up cell 23237 / edge 1018 (26.9°N, 0.9°W, level 1 beside level 0) by 03:08. The YSU increment `rublten_Edge * dt_sub` matches `physics_get_tend` over one short step. The failure is the operator split: several explicit PBL increments, rates zeroed, then one friction-free SRK3 of `dt_level`. An edge takes the finer cell, so edge 1018 is level 0 and the level-1 fold never updates it. The default is again `config_hts_pbl_substeps = 1`, one surface call of the level dt. Explicit 0 still selects the fold. `micro_substeps > 1` with `pbl_substeps <= 1` warns once: that is the edge 64121 restart under one long YSU rate.

`config_hts_max_level = 0` does not turn `config_hts` off and does not run the hold. Setup clamps the mesh to a single level and returns before the level lists. The time step is one SRK3 of `config_dt`. `config_hts_debug_force_path = true` with that cap keeps every cell at level 0 and still runs the lists, the acoustic halo, snapshot and restore, flux accumulation, the idle-scalar copy, and HTS column physics. Job 201473 matched all-fine. A mesh that already has a coarse level logs that the switch is ignored.

The integers:

- `config_hts_moist_max_level` (default 2). Dry runs ignore it.
- `config_hts_micro_substeps` (default 1). `0` steps dynamics, scalars, and Kessler at `config_dt` on every moist level above 0 and averages `ruAvg` before the single mass correction. A positive integer is that many steps on the coarsest moist level only. The acoustic dt stays `config_dt / config_number_of_sub_steps`.
- `config_hts_pbl_substeps` (default 1). One call of the level dt. `0` folds a `config_dt` surface increment and then takes one dynamics step of the level dt. When `micro_substeps` is also 0, the surface count matches the short dynamics step and the fold is not used. A positive integer other than 1 is that many calls on the coarsest moist level only. An explicit mismatch with both counts greater than 1 is raised to the shorter dt.
- `config_hts_kessler_substeps` (default 1). `0` calls Kessler at `config_dt` inside the one dynamics step on every moist level above 0. Ignored when `micro_substeps` is 0 or greater than 1.

Debug switches, all default off. `config_hts_debug_match_fine` turns on the first two only. It does not skip the idle copy, the mass correction, the promote, or the `ruAvg` mean, and it does not keep acoustic `*_pp`.

- `config_hts_debug_refresh_diag`. Before each column-physics call on a level above 0, rebuild `exner`, `pressure_p`, `rho_p`, and `rtheta_p` from time level 1 and reconstruct the cell wind after copying `u` time level 2 from time level 1 on that level's edges.
- `config_hts_debug_rad_every_substep`. While a coarse step takes more than one substep, restore a latched radiation or convection alarm before substeps 2..N. A single level-dt call is unchanged, so on the 60 s baseline this flag does nothing.
- `config_hts_debug_skip_idle_copy`. Skip the scalars inside the idle-state copy. `rho_zz`, `theta_m`, `w`, and `u` are still copied. The 8db2f20c9 switch skipped the whole scalar copy and left every prognostic on the previous buffer.
- `config_hts_debug_skip_mass_correct`. Skip the interface mass correction.
- `config_hts_debug_ruavg_last`. Leave `ruAvg` at the last substep. `wwAvg` is already rebuilt inside each SRK3. The mass correction then sees that last flux against the full-interval target.
- `config_hts_debug_skip_promote`. Do not copy time level 2 onto time level 1 between dynamics substeps.
- `config_hts_debug_keep_pp`. Do not zero `rho_pp`, `rtheta_pp`, and `rw_p` on inactive acoustic neighbors. That is the job 199953 configuration. It can blow the interface and is not part of `match_fine`.
- `config_hts_debug_probe_cell`. Fortran cell index, printed every fine step on rank 0. 0 disables it. The eds history cell 3177 is 3178. The line is `HTS probe tick cell own L closed k w th qv rho t2 tsk usfc gmax|w| th1 qv1 td1`. `k` is the vertical index of that `w` (1 is the surface). `L` is the mesh level. `closed` is 1 when that level stepped on this tick. `th` and `qv` are `theta_m` and water vapor at that `k`. `th1` and `qv1` are `theta_m` and water vapor at k=1. `td1` is `theta_m / (1 + Rv/Rd * qv)` at k=1, the dry potential temperature history writes. `qv1` is kg/kg.

Run B is the config_dt path through level-1 machinery: `config_hts_max_level = 1`, `config_hts_micro_substeps = 0`, `config_hts_pbl_substeps = 0`. The log must say `coarse level 1: micro substeps 2, pbl substeps 2, dt 60` and must not warn about an acoustic restart.

The 8db2f20c9 bisect ran that namelist for 8 h with one switch at a time (probe cell 3178, all-fine reused). Cell 3177 `|w|`:

| t | allfine | B | mf | rad | refresh | mass | ruavg | promote | idle |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 04:00 | 0.047 | 0.066 | 0.068 | 0.061 | 0.066 | 0.066 | 0.066 | 0.066 | 0.044 |
| 05:00 | 0.048 | 0.123 | 0.126 | 0.125 | 0.123 | 0.123 | 0.123 | 0.122 | 0.066 |
| 05:30 | 0.051 | 0.192 | 0.197 | 0.195 | 0.192 | 0.193 | 0.190 | 0.181 | 0.079 |
| 06:00 | 0.059 | 0.294 | 0.299 | 0.298 | 0.294 | | | | 0.102 |
| 07:00 | 0.077 | 0.904 | | | | | | | 0.184 |
| 08:00 | 0.112 | 4.069 | | | | | | | 0.388 |

`match_fine`, radiation, and diagnostic refresh tracked B and were cancelled around 06:00–06:30. Mass correction, last `ruAvg`, and skip-promote tracked B and were cancelled around 05:30. Skip-idle ran the full 8 h, critical 0, at 0.834×. Edge 41342 `|u|` followed B on every leg except skip-idle, which stayed near 25.0 against all-fine 24.9 (B was 60.6). Global max `|u|` stayed within 0.05 m/s on every leg. `k ≤ 5` and edge 64121 stayed quiet. The mf probe at 00:30 already had dry theta about −0.5 K and `qv` about +1.05 g/kg at k0 against all-fine, on ocean with SST and skin temperature frozen at 320.5 K. `|w|` growth starts near 04:00 as the modal level drops from about k=36 to k=11 (about 13:00 local). Growth ratios on radiation ticks and on other ticks were about 1.00 and 1.01.

No switch removed the growth. Skip-idle delayed the onset by about 1–1.5 h and cut 08:00 `|w|` from 4.07 to 0.388, and the column was still growing by about ×1.4 per 30 min after 05:00.

`95db9939c` copied `rho_zz`, `theta_m`, `w`, and `u` with the scalars before the shift. Job 201673 passed the 2 h safety gate (max `|u|` 89.33/89.34, `w` relative L2 0.19/0.27/0.34/0.40, T2/TSK rms 0.15/0.47 K, dry-mass drift near 1e-9, 0.88×). The 8 h tracked old B to three digits (0.841×): 00:30 k0 was still −0.495 K and +1.045 g/kg, and 08:00 `|w|` was 4.076 against old B 4.069 and all-fine 0.112. Edge 41342 `|u|` was 60.7/24.9. `w` relative L2 at 08:00 was 1.401. Open ticks near 00:30 repeated the previous close in `th1` and in `w`/`th`/`qv` at the level of max `|w|`, while `qv1` and `td1` still moved.

That probe is the write. With `config_monotonic` true (the default), RK stages 1 and 2 use the masked scalar update and the final stage uses `atm_advance_scalars_mono_work`. The final stage adds `scalars_tend * dt / rho_zz` to time level 1 for every owned cell, then transports, then copies the result back to time level 2 only on the active level. `physics_get_tend` fills `scalars_tend` from `rqvblten` on every cell, including columns whose physics did not run on this tick, and `hts_restrict_physics_tendencies` was zeroing `tend_rtheta_physics` but not `scalars_tend`. On an open tick the level-1 column therefore gains the previous surface moisture flux in `qv` and keeps the old `theta_m`. The same addition runs again on the level-0 step that precedes the close, and the coarse RK setup then copies that time level 1 into the start of the level-1 step. History theta is `theta_m / (1 + Rv/Rd * qv)`. At about 300 K, an extra 1.05 g/kg with no matching `theta_m` increment is about −0.5 K. Two extra 30 s applications per 60 s match one full surface flux, which is the 00:30 bias. The prognostic copy published that time level 1, so it could not remove the split. Skip-idle left the extra `qv` in time level 1 and shifted the previous time level 2, which is why the growth only slowed.

The final RK stage now skips that source update on a column whose level is not active, and `hts_restrict_physics_tendencies` zeros `scalars_tend` there too. Active columns, a single level, and all-fine still take the stock update. The idle prognostic copy stays: it keeps time level 2 equal to the committed time level 1 on a column that did not step.

The next eds pair is the same run B, probe on, other debug switches off:

```
config_hts = true
config_hts_max_level = 1
config_hts_micro_substeps = 0
config_hts_pbl_substeps = 0
config_hts_debug_probe_cell = 3178
```

Rebuild the dynamics objects. `gen_includes` is not required when the tree already has the 8db2f20c9 namelist. The probe line must contain `th1`, `qv1`, and `td1`. On an open tick, `qv1` and `td1` should repeat the previous close, as `th1` already did. 2 h is the safety check (no polar fold, edge 64121 quiet, dry-mass drift near 1e-9, `k ≤ 5` `|u|` within a few m/s). 8 h is the cell 3177 test: 00:30 k0 dry theta and `qv` against the old −0.5 K / +1.05 g/kg, and 08:00 `|w|` against all-fine 0.112 and the old B value 4.069. A quiet 2 h is not a pass. The local Kessler mesh has PBL, surface, and radiation off, so it cannot see this column.

`pbl_substeps = 1` with `micro_substeps = 1` is the reverted 60 s baseline and should reproduce job 201388. An explicit `pbl_substeps = 0` with `micro_substeps = 1` is the fold that already failed. The 2 h gate still includes the max per-edge `|u|` difference and max `|u|` for `k ≤ 5`.

Acoustic work on the 120 s cap does not change relative to the 240 s split: both take 48 acoustic steps of 5 s per 240 s. The measured 120 s pair was 1.48× at 2 h and 1.40× at 8 h (jobs 201325 and 201338). Calling PBL at a fine dt and then multiplying the stored rate by a long dynamics step is the `902a6084d` failure and is not the default.
