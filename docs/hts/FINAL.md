# Hierarchical time stepping: where it landed

Branch `hts` of this fork, base `upstream-base` (`f34984b2`, MPAS-Atmosphere v8.4.2). `config_hts` defaults to false. The all-fine path of the HTS binary matches the non-HTS history bit for bit on the gates below. Writeups:

- [DESIGN.md](DESIGN.md) — scheme and the physics contract
- [BASELINE.md](BASELINE.md) — unmodified v8.4.2 JW baselines
- [PHASE2.md](PHASE2.md) — held coarse state, matched acoustic step
- [PHASE2B.md](PHASE2B.md) — index lists and the level-balanced partition
- [MESH_SMALL_PATCH.md](MESH_SMALL_PATCH.md) — the 24755-cell circular patch
- [PHASE2C.md](PHASE2C.md) — dry speed and accuracy on that mesh
- [PHASE3.md](PHASE3.md) — physics coupling and the real-data run
- [PHASE3B.md](PHASE3B.md) — column physics skips inactive columns

All of the timings are 4 MPI ranks, OpenMP off, on this VM, single precision.

## Accuracy against all-fine

Relative L2 is `rms(test − reference) / rms(reference)` on the last history frame. Mass drift is the relative change in dry mass from the first frame.

| run | `u` | `w` | `theta` | `rho` | HTS mass drift | all-fine mass drift |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 120 km, 1 level, dry, 10 steps, and again after the phase 3 binary | 0 | 0 | 0 | 0 | −1.98×10⁻¹⁰ | −1.98×10⁻¹⁰ |
| 120 km, day 10, phase 2, dry | 0 | 0 | 0 | 0 | −2.12×10⁻¹⁰ | −2.12×10⁻¹⁰ |
| `x5.30210`, 6 h, 2 levels, dry | 2.74×10⁻⁴ | 0.180 | 3.61×10⁻⁶ | 2.41×10⁻⁶ | −3.23×10⁻¹⁰ | −2.39×10⁻¹⁰ |
| `x4.10242`, 6 h, 3 levels, dry | 1.09×10⁻³ | 0.444 | 1.02×10⁻⁵ | 1.61×10⁻⁵ | −8.04×10⁻¹⁰ | −6.98×10⁻¹¹ |
| `sp8.24755`, 6 h, 4 levels, dry, `dt = 30` s | 7.89×10⁻⁴ | 0.438 | 3.68×10⁻⁵ | 4.39×10⁻⁶ | −1.27×10⁻⁹ | −7.95×10⁻¹¹ |
| `sp8.24755`, 1 h, 4 levels, Kessler | 6.16×10⁻⁴ | 0.293 | 8.41×10⁻⁶ | 1.28×10⁻⁵ | −6.89×10⁻¹¹ | −9.60×10⁻¹¹ |
| `sp8.24755`, 30 min, 4 levels, WSM6 | 3.63×10⁻⁴ | 0.486 | 3.67×10⁻⁶ | 1.03×10⁻⁵ | −7.91×10⁻⁹ | −1.31×10⁻¹⁰ |
| `sp8.24755`, 15 min, 4 levels, Thompson | 2.12×10⁻⁴ | 0.300 | 2.34×10⁻⁶ | 6.20×10⁻⁶ | +4.03×10⁻⁹ | −8.70×10⁻¹¹ |
| `sp8.24755`, 8 min, 4 levels, RRTMG + YSU + Noah + Kessler | 3.66×10⁻⁴ | 0.267 | 2.24×10⁻⁵ | 1.67×10⁻⁴ | −1.44×10⁻¹⁰ | −5.28×10⁻¹¹ |
| `x1.40962` real, 6 h, 1 level, `mesoscale_reference` | 0 | 0 | 0 | 0 | −4.31×10⁻¹¹ | −4.31×10⁻¹¹ |

`w` relative L2 is large because the reference `w` is small (about 10⁻⁴ to 10⁻³ m/s on the dry JW runs). The absolute rms difference on `sp8.24755` at 6 hours is 4.6×10⁻⁴ m/s. The largest piece sits on the interface. `u`, `theta`, and `rho` stay at relative 10⁻⁴, 10⁻⁵, and 10⁻⁶. Dry-mass drift on the fast four-level dry run is 1.27×10⁻⁹. The previous full-domain-diagnostic binary on that same frame was 6.96×10⁻¹⁰, and `u` relative L2 between those two HTS histories is 1.21×10⁻⁵.

The held-state skip below does not move the dry 6-hour or Kessler 1-hour history: both match the previous HTS binary bit for bit on `u`, `w`, `theta`, `rho`, and `pressure`. WSM6, Thompson, and the radiation smoke move against that previous HTS history by `u` relative L2 1.96×10⁻⁶, 1.22×10⁻⁶, and 6.93×10⁻⁷. Against all-fine the three significant figures in the table are unchanged. WSM6 dry-mass drift moves from −7.93×10⁻⁹ to −7.91×10⁻⁹, and the radiation smoke from −1.47×10⁻¹⁰ to −1.44×10⁻¹⁰. Thompson stays at +4.03×10⁻⁹. The all-fine side of every pair matches its previous history bit for bit.

The layer-1 acoustic halo below does not move any of those histories. Dry, Kessler, WSM6, Thompson, and the radiation smoke match the held-state binary bit for bit on the same five fields, and each all-fine run matches its previous history bit for bit. Against the paired all-fine history the dry 6-hour figures stay at 7.89×10⁻⁴ / 0.438 / 3.68×10⁻⁵ / 4.39×10⁻⁶.

`dt = 60` s on `sp8.24755` is omitted from the table. All-fine finished. Every HTS clamp signaled invalid and overflow. The acoustic step stays the fine step. The comparisons below use `dt = 30` s. The CFL numbers are in the next section.

## Speed against all-fine

Ratio is all-fine / HTS on `time integration`. Higher than 1 is faster. The `sp8.24755` rows are 4 MPI ranks, OpenMP off, `dt = 30` s. The bold rows are this binary. Earlier binaries are labeled. An all-fine history that matches bit for bit can still move by about 2 s between sittings, so each bold ratio uses the paired run from the same sitting.

| case | all-fine | HTS | ratio |
| --- | ---: | ---: | ---: |
| `x5.30210`, 6 h, 2 levels, 55% of cells fine | 25.06 s | 27.34 s | 0.92× |
| `x4.10242`, 6 h, 3 levels, 49% of cells fine | 3.62 s | 4.25 s | 0.85× |
| `sp8.24755`, 6 h, dry, phase 2c binary | 127.10 s | 104.05 s | 1.221× |
| `sp8.24755`, 6 h, dry, diagnostic-ring binary | 120.66 s | 80.66 s | 1.496× |
| `sp8.24755`, 6 h, dry, held-state binary | 122.81 s | 66.50 s | 1.847× |
| `sp8.24755`, 6 h, dry | 125.98 s | 62.95 s | **2.001×** |
| `sp8.24755`, 1 h, Kessler, phase 3 restore | 23.72 s | 23.04 s | 1.029× |
| `sp8.24755`, 1 h, Kessler, phase 3b column mask | 23.35 s | 19.02 s | 1.227× |
| `sp8.24755`, 1 h, Kessler, diagnostic-ring binary | 22.98 s | 14.35 s | 1.601× |
| `sp8.24755`, 1 h, Kessler, held-state binary | 23.73 s | 11.69 s | 2.031× |
| `sp8.24755`, 1 h, Kessler | 22.86 s | 11.17 s | **2.047×** |
| `sp8.24755`, 30 min, WSM6, phase 3b binary | 15.91 s | 10.65 s | 1.494× |
| `sp8.24755`, 30 min, WSM6, diagnostic-ring binary | 16.15 s | 7.69 s | 2.100× |
| `sp8.24755`, 30 min, WSM6, held-state binary | 16.47 s | 6.42 s | 2.567× |
| `sp8.24755`, 30 min, WSM6 | 16.92 s | 6.22 s | **2.720×** |
| `sp8.24755`, 15 min, Thompson, diagnostic-ring binary | 7.32 s | 3.71 s | 1.971× |
| `sp8.24755`, 15 min, Thompson, held-state binary | 7.30 s | 3.04 s | 2.398× |
| `sp8.24755`, 15 min, Thompson | 7.22 s | 2.96 s | **2.439×** |
| `sp8.24755`, 8 min, RRTMG + YSU + Noah + Kessler, diagnostic-ring binary | 6.49 s | 4.81 s | 1.349× |
| `sp8.24755`, 8 min, RRTMG + YSU + Noah + Kessler, held-state binary | 6.64 s | 4.26 s | 1.560× |
| `sp8.24755`, 8 min, RRTMG + YSU + Noah + Kessler | 6.50 s | 4.17 s | **1.559×** |
| `x1.40962` real, 6 h, 1 level, full suite, phase 3 binary | 226.16 s | 218.99 s | bit-identical; the 3% is timer scatter |

Wall (`/usr/bin/time`, all-fine / HTS) on the bold rows: dry 130.73 s / 67.81 s, Kessler 24.50 s / 12.95 s, WSM6 18.28 s / 7.66 s, Thompson 8.72 s / 5.11 s, radiation 7.63 s / 5.37 s.

## What the timers paid for

Level 0 is 885 of 24755 cells. The spacing ratio is 10×. The large-step tendency ceiling is 6.10×. On the 6-hour dry pair the loops that were walking the whole domain, and that the timers showed growing under HTS, are now on a list:

| 6 h dry | all-fine | HTS | calls |
| --- | ---: | ---: | ---: |
| `atm_compute_dyn_tend` | 58.24 s | 11.95 s | 2160 → 4050 |
| `atm_compute_solve_diagnostics` | 10.89 s | 2.32 s | 2160 → 4050 |
| `atm_compute_moist_coefficients` | 0.902 s | 0.176 s | 720 → 1350 |
| `atm_compute_vert_imp_coefs` | 1.120 s | 0.199 s | 720 → 1350 |
| `small_step_prep` | 3.064 s | 0.504 s | 2160 → 4050 |
| `hts_acoustic_halo` | — | 4.79 s | 69120 |
| `atm_advance_acoustic_step` | 28.03 s | 28.24 s | 8640 → 34560 |
| `atm_divergence_damping_3d` | 8.28 s | 8.87 s | 8640 → 34560 |
| `physics_get_tend` | 0.601 s | 0.618 s | 720 → 1350 |

Diagnostics use four rings (rank 0, level 3: 6127 cells of the local domain; levels 0–2 are a few hundred cells). The previous HTS binary spent 16.56 s in diagnostics on this run. Moist coefficients and the vertical implicit coefficients follow the dynamics work list and the active columns. `physics_get_tend` still zeros the full scalar-tendency array and then fills it. Its 6-hour cost is 0.67 s, so that zero stays full-domain: an off-level leftover would persist into the next level.

The acoustic and damping gaps against all-fine are the extra calls at `dts = 5` s, not edges that are still updated when they are inactive. Both loops walk the active-edge list and cycle an edge whose level is not the active one, so an inactive edge is not updated. Each of the four levels does 8640 acoustic calls. That is the same 8640 calls all-fine spends on the whole mesh, split across levels, because every level keeps the fine acoustic step. A four-level call therefore does about a quarter of an all-fine call:

| | all-fine per call | quarter of that | HTS per call | HTS total − all-fine total |
| --- | ---: | ---: | ---: | ---: |
| `atm_advance_acoustic_step` | 28.03 s / 8640 = 3.244 ms | 0.811 ms | 28.24 s / 34560 = 0.817 ms | 0.21 s |
| `atm_divergence_damping_3d` | 8.28 s / 8640 = 0.958 ms | 0.239 ms | 8.87 s / 34560 = 0.257 ms | 0.59 s |

The 0.006 ms and 0.017 ms left on each HTS call, times 34560 calls, are the 0.21 s and 0.59 s gaps. That is loop overhead on four times as many entries. Dropping the update on an active edge would change the fine-edge flux, and `dts` is already 5 s, so that compute stays.

`hts_acoustic_halo` is the per-substep layer-1 exchange, 4.79 s over 69120 calls. The untied remainder of `time integration`, after every direct child timer, is 6.15 s all-fine and 1.95 s with four levels. The all-fine remainder still contains its acoustic halo. On the HTS side that 1.95 s is what the large-step groups look like before they are timed on their own; the split is below. Every rank enters every collective. After the level-balanced partition every rank owns every level, so a per-level group would still include every rank. The acoustic-step parallel efficiency on this partition is 0.96. The compact halo's is 0.52, because level 3 still sends most of the layer-1 halo and the coarse levels send almost none. Level weights are already in the graph. A rank that skipped a collective would deadlock.

`mp_kessler` on the hour is 1.044 s / 0.151 s = 6.91×. `mp_wsm6` on the half hour is 3.799 s / 0.643 s = 5.91×. Both sit on the 6.10× cell-count ceiling. `mp_thompson` on 15 minutes is 1.052 s / 0.350 s = 3.00×. Thompson keeps its 90 s internal substep, and a 15-minute run is short relative to that split.

## Thompson

`build_tables` (`src/core_atmosphere/utils/build_tables.F`) wrote the four lookup files in 176.5 s:

- `MP_THOMPSON_QRacrQG_DATA.DBL`
- `MP_THOMPSON_QRacrQS_DATA.DBL`
- `MP_THOMPSON_freezeH2O_DATA.DBL`
- `MP_THOMPSON_QIautQS_DATA.DBL`

They live outside the repo. The 15-minute smoke is `mp_thompson` with convection, PBL, radiation, the surface layer, and the land model off, the same shape as the Kessler namelist. `qi`, `qs`, and `qg` start at zero. IEEE flags are underflow and denormal. `u` stays about ±35 m/s. Integration is **2.439×**. Dry-mass frame difference is +4.12×10⁻⁹.

## Radiation and PBL

The JW init for `sp8.24755` has blank land use, `ivgtyp = 0`, `isltyp = 0`, `xland = 0`, `ter = 0`, and `smois = 0`. A copy used for this smoke sets USGS water (`ivgtyp = 16`, `isltyp = 14`, `xland = 2`, `mminlu = 'USGS'`), which matches `iswater_lu = 16` on that file. GWDO stays off because `var2d` is zero. Convection stays off. The schemes that run are `rrtmg_lw`, `rrtmg_sw`, `cld_fraction`, `sf_monin_obukhov_rev`, `sf_noah`, `bl_ysu`, and `mp_kessler`. Radiation interval 30 minutes. Duration 8 minutes, so the t = 0 alarm is the only radiation firing. It is latched until each level steps, then cleared. Shortwave time is 1.592 s in one all-fine call and 1.602 s across the masked segments. `bl_ysu` is 0.563 s / 0.066 s = 8.56×.

The lowest-layer `theta` on levels 0–2 matches all-fine to about 0.001 K. On level 3 the same layer warms by a mean of 0.266 K (HTS) and 0.299 K (all-fine). The mean difference is −0.034 K and the maximum is 0.38 K, on cells thousands of kilometres from the fine patch. Both runs warm that layer because the filled water static has zero soil moisture. Setting `dt_pbl` to the level `dt` (240 s on level 3) raised that mean difference to 0.093 K and the maximum to 0.44 K on the same pair, so `dt_pbl` stays at the fine step stored at initialization. YSU, Noah, sea ice, and GWDO return a rate normalized by that `dt`, and the level RK multiplies the stored rate by `dt_level`.

Ozone, monthly albedo, and the SST and sea-ice file updates stay on the fine clock for every owned column. They are file-driven fields, not a per-column rate the level step can rescale. GWDO and Tiedtke were not timed here. The JW static has blank land use and zero orography, `var2d` is zero, and a land variable-resolution forecast cannot be initialized on this mesh. The GWDO driver reads the precomputed `var2d` array and does not index neighbor cells, so the column mask is the same skip used for YSU once a real static exists. That case was not run.

A pass over the column drivers found no further stencil bug on a variable mesh. Radiation, PBL, convection, the surface layer, and microphysics already step through `column_mask_next_segment`. Skin SST, the ocean mixed layer, and the surface diagnostics use `column_is_active`. YSU, RRTMG, and GWDO are handed an `its:ite` segment and do not read `cellsOnCell`, so a segment boundary is not a mesh neighbor. The ozone climatology, the monthly albedo and greenness interpolation, and the SST and sea-ice file update still walk every owned column when the fine-clock alarm rings. Those fields follow the calendar, and the value written on the alarm is the value a later coarse step would interpolate. Skin SST and the deep-soil temperature integrate a rate; they are deferred and use the level `dt`. Scalar advection stays off in these namelists. `x1.40962` is one level, so it does not exercise a variable-mesh radiation path.

## Acoustic halo and the coarse-state hold

Two changes to the per-substep acoustic halo were timed on the 6-hour dry pair and removed. The rows above the bold ones keep the separate full-domain `dynamics:rho_pp` and `dynamics:rtheta_pp` groups. The bold rows add the compact list at the end of this section.

Fusing the two layer-1 fields into one `mpas_halo` group is bit-identical on the 10-step 120 km gate, the 2-step `x1.40962` `mesoscale_reference` smoke, and the 6-hour dry all-fine and four-level histories. It is slower. Integration went from 120.66 s to 127.14 s all-fine and from 80.66 s to 83.21 s with four levels (walls 125.47 s → 132.02 s and 85.62 s → 88.46 s). Default `config_halo_exch_method` is already `mpas_halo`, which packs a group into one message per neighbor. Putting both fields in that message does not reduce the bytes, and the pack dominates the wait.

Packing only the halo cells an active-level edge reads does cut the bytes. On this 4-rank partition the layer-1 halo is 312–355 cells per rank. The active-edge subset, summed across ranks, is 149, 50, 31, and 1096 cells on levels 0–3. Each level does the same number of acoustic substeps, so the packed volume is a quarter of exchanging the full halo every substep. A timer on that exchange recorded 4.08 s over 69120 calls. Integration was 79.68 s (wall 84.76 s) and the paired all-fine run was 122.82 s (wall 127.68 s). The previous four-level integration was 80.66 s, and all-fine moved from 120.66 s to 122.82 s between sittings, so the 1 s change sits in the timer scatter. The 6-hour state is not the previous HTS history: `u` relative L2 1.29×10⁻⁵, `w` 1.17×10⁻², `theta` 5.92×10⁻⁷, `rho` 3.07×10⁻⁷. Against all-fine on the same binary the dry errors stay at the published values (`u` 7.89×10⁻⁴). The subset exchange was removed with the fused group.

Splitting that timer out left 18.4 s of `time integration` with no child timer, against 24.3 s when the acoustic exchanges were inside the remainder. About 6 s of the old remainder was the acoustic halo. The other 18 s was the coarse-state hold. `atm_hts_apply_interpolation` and `atm_hts_restore_above` rewrote every coarser column on each fine step, and `hts_refresh_wind` called `mpas_reconstruct` on the owned domain. Over 6 hours that reconstruction ran once before each of the 2070 level steps and once at the end of every `atm_srk3`.

The value written by the hold is the start-of-interval snapshot (`alpha = 0`). Extrapolating the previous increment fed a coarse acoustic transient back into the interface flux, and that stays out. After the snapshot is written, both time levels match it and stay matched across the time-level shift until that level steps. A finer step does not write those prognostic fields. Copying the snapshot in again is a no-op, so a level with `hts_hold_valid` set is skipped. The flag is cleared only when that level steps, and the next fine step copies the new snapshot into both time levels once per level interval. Cells are written before edges so an edge density sees the updated columns.

The cell wind follows the same rule. `hts_refresh_wind` reconstructs only the owned cells incident to an edge that level just rewrote. The end-of-`atm_srk3` reconstruct on the HTS path does the same for the active level. Cells whose edges did not change keep the previous reconstruction, which is what a full-domain call would have written. The non-HTS path, and a one-level run, still call `mpas_reconstruct` on the owned domain.

On the 6-hour dry run the new timers are `hts_begin_intervals` 0.151 s (720 calls), `hts_apply_interpolation` 0.518 s (1260), `hts_refresh_wind` 0.259 s (360), `hts_restore_above` 0.000 s (1260), and `hts_reconstruct_wind` 0.374 s (1350). Integration is 66.50 s. The dry and Kessler histories match the previous HTS binary bit for bit.

The acoustic step still had one full-domain copy: `rtheta_pp_old` on every cell, and a first-small-step zero of `rho_pp`, `rtheta_pp`, and `rw_p` on every inactive owned cell. Divergence damping reads `rtheta_pp_old` only on the two cells of each active edge, so both loops now walk that cell list. Acoustic time went from 29.42 s to 28.43 s. Damping stayed at 8.98 s. Dropping the damping or the acoustic update on those edges would change the fine-edge flux, and `dts` is already the fine step, so that compute stays. On that binary the untied remainder is 9.51 s, down from 24.3 s, and what remained of it was the full-domain halo.

The per-substep exchanges are layer 1, and the acoustic step reads `rho_pp` and `rtheta_pp` only on the two cells of each edge. Those two groups now use a compact list. A halo cell is on the list for level L when its cell level is L or an incident edge has edge level L. That test is the mesh, not which rank owns the cell, so both ends of a message list the same cells. Send and receive counts are checked with an all-to-all. A mismatch leaves the compact path off and the run keeps `exchange_halo_group`. Non-HTS and one-level runs never enter it. On this 4-rank partition rank 0 sends 41, 12, 9, and 256 of 311 layer-1 halo cells on levels 0–3, and the lists matched. The 6-hour dry integration is 62.95 s against a paired all-fine 125.98 s. `hts_acoustic_halo` is 4.79 s. The histories listed above match the held-state binary bit for bit.

Level 3 is 256 of those 311 cells and does the most substeps, so the message that runs most often barely shrinks. The groups that stay on the full halo are wider than the two cells of an active edge. The tendency reads two rings, diagnostics read four, and `u` is already exchanged at layer 3 only. Cutting those lists to the layer-1 edge endpoints would drop halo cells the next operator still reads.

Those groups were timed on one 6-hour dry run by wrapping the existing `exchange_halo_group` calls. The history matched the compact-halo binary bit for bit. The printed timer is the slowest rank; rank-average is that value times the parallel efficiency. Integration on that sitting was 64.00 s, against 62.95 s for the uninstrumented binary. The same all-fine history moves by about 2 s between sittings.

| group | calls | slowest rank | rank average |
| --- | ---: | ---: | ---: |
| `theta_m`, scalars, `pressure_p`, `rtheta_p` | 1350 | 0.79 s | 0.63 s |
| `exner` | 1350 | 0.35 s | 0.22 s |
| `tend_u` | 4050 | 1.46 s | 0.87 s |
| post-acoustic `rw_p`, `ru_p`, `rho_pp`, `rtheta_pp` | 4050 | 0.97 s | 0.92 s |
| `u` layer 3 | 4050 | 0.65 s | 0.56 s |
| `w`, `pv_edge`, `rho_edge` | 4050 | 1.39 s | 1.36 s |

The six averages sum to 4.56 s. That is the time a rank spends in the collectives, not time that can be deleted. The active level still has to exchange its stencil, and every rank still enters every call. The layer-1 acoustic shrink saved about 3.5 s across 69120 calls, about 50 μs a call, and those messages were a single field. The same 50 μs on these 18900 calls is about 1.0 s, 1.5% of the 62.95 s integration. The wide groups are already one message per neighbor. A shorter payload does not remove the collective, and splitting a group into per-field messages adds one. The byte cut is also smaller than the layer-1 edge-endpoint cut, because the safe stencil is two to four rings and level 3 still owns most of the halo on the calls where the message is large. One second is inside the sitting scatter, so the filters were not added and the timer wrappers were not kept. Nothing left in the dynamics is both above that 1–2% bar and inside the published error table.

## Safe step

`config_dt = 30` s, acoustic `dts = 5` s, `config_number_of_sub_steps = 6`. A 6-hour attempt at `config_dt = 60` s (`dts = 10` s) is written up in [PHASE2C.md](PHASE2C.md). All-fine finished (wall 64.47 s, underflow and denormal). Every HTS clamp signaled `IEEE_INVALID` and `IEEE_OVERFLOW`. With `u ~ 35` m/s and `dt_level = 60 * 2^L`, the advective CFL on the shortest edge of each level stays under 0.21. The fine-edge acoustic CFL at `dts = 10` s is about 0.22. The overflow is the coarse and interface acoustic update. The acoustic step was left at the fine step.

## Turning it on

`config_hts` defaults to false, and false does not enter the masked loops. Set `config_hts = true` and `config_hts_max_level` to the coarsest level to run (an integer from 0 through 8). `config_dt` is the level-0 step. Level L advances by `config_dt * 2^L`. Leave `config_number_of_sub_steps` at the fine acoustic count, so `dts = config_dt / config_number_of_sub_steps`. The integrator multiplies that count by `2^L` on coarser levels. Do not lengthen `dts` to the coarse step.

A cell's level comes from `meshDensity`. Spacing relative to the finest cells is `meshDensity^(-1/4)`, and `L = nint(log2(max(that, 1)))`, clamped to `config_hts_max_level`. An edge takes the finer of its two cells. A uniform mesh (`meshDensity` identically 1, or no `nominalMinDc`) is all level 0, and HTS on with one level matches HTS off bit for bit.

Every rank should own cells of every level. `docs/hts/scripts/hts_partition.py GRID.nc GRAPH.info NPARTS --max-level N` writes the partition file named by `config_block_decomp_file_prefix`. A default partition that puts the whole fine patch on one rank leaves the other ranks idle inside the halo exchange.

An OpenACC build aborts if `config_hts` is true. The masked loops are compiled out of that path.

A variable-resolution forecast with radiation, PBL, and the land model needs a real static: land use, soil type, orography, and the GWDO `var2d` field. The JW init for `sp8.24755` has none of those; the 8-minute cloud smoke filled USGS water by hand. On eds, personal WPS_GEOG now supplies a real static on this mesh (see [EDS_BUILD.md](EDS_BUILD.md)); the land-VR short forecast is documented there. Remaining known limits are the acoustic fine-`dts` floor, wide-halo residual below the noise floor, and that this remains a research prototype—not production NWP.

## Gates

`config_hts` defaults to false in `Registry.xml`. One level returns before the lists are built, and the masked loops require `hts_max_level_used > 0`, so they keep the original index order.

Uniform 120 km, physics suite `none`, 10 steps, `dt = 450` s. HTS on, HTS off, and each against its previous history: relative L2 of `u`, `w`, `theta`, and `rho` is 0. Drift −1.977×10⁻¹⁰.

`x1.40962`, CFSR `2010-10-23_00`, suite `mesoscale_reference`, `dt = 720` s, 2 steps. HTS on against HTS off, and each against the previous smoke history: relative L2 0 on `u`, `w`, `theta`, `rho`, and `pressure`. Drift −1.799×10⁻¹¹. The 6-hour forecast that matches bit for bit was run on the phase 3 binary.

## Bottom line

On `sp8.24755` at `dt = 30` s, four levels, this binary is **2.001×** dry (6 h), **2.047×** Kessler (1 h), **2.720×** WSM6 (30 min), **2.439×** Thompson (15 min), and **1.559×** for an 8-minute RRTMG/YSU/Noah/Kessler smoke. The one-level and HTS-off histories match bit for bit. Dry, Kessler, WSM6, Thompson, and the radiation smoke also match the held-state HTS binary bit for bit, and the vs-all-fine figures stay at the values in the table. The acoustic step stays the fine step, because `dt = 60` s overflows on this mesh while the linear CFL stays small. Of the 62.95 s dry integration, the acoustic update is 28.24 s and the divergence damping is 8.87 s, against 28.03 s and 8.28 s all-fine. Both already walk the active-edge list, and the gap is the quadrupled call count at `dts = 5` s. The coarse-state hold is about 1.3 s. The per-substep layer-1 halo is 4.79 s on the cells the active level touches. The other halo groups were timed at 4.56 s rank-average and left on the full stencil: a safe shrink of those messages is about a second, inside the sitting scatter. Fusing the two acoustic messages made an earlier binary slower, and an ownership-filtered pack of the same halo moved the 6-hour state by `u` relative L2 1.3×10⁻⁵; both stayed out. The compact list that landed uses the cell level and the incident edge level, and the all-to-all count check has to match before it is used. A land-surface variable-resolution forecast on this mesh is now unblocked on eds with personal WPS_GEOG (`/eds/energy/dev01/lvhao/WPS_GEOG/`). Real static (`work/hts_land_static/sp8.24755.static.nc`, job 199304) has landmask/terrain/MODIS land use and native GWD `var2d`. Real GFS IC + 15-minute Noah/YSU/RRTMG/Kessler HTS vs all-fine (job 199311) finished with 0 critical errors and **1.57×** on `time integration` (46.37 / 29.55 s). Land `skintemp`/`t2m` rms differences at 15 min are ~0.4 / 0.3 K. That does not change the acoustic fine-`dts` floor or the wide-halo residual documented above; higher overall acceleration than the published dry/moist tables is still a known limit, not a remaining task. GWDO was left off in the land-VR short run (optional). Noah-MP static (`soilgrids/`) was not downloaded; `config_noahmp_static=false` and the atmosphere used `sf_noah`. The earlier JW radiation smoke (USGS-water fill) remains the cloud reference for the 1.559× RRTMG row. See [EDS_BUILD.md](EDS_BUILD.md) and `work/README_LAND_VR.md`.

`src/core_atmosphere/dynamics/Makefile` depends `mpas_atm_dissipation_models.o` on `mpas_atm_hts.o`, so a parallel `make` builds `mpas_atm_hts.mod` first. `manage_externals` skips `git remote add` when that remote name is already present, so a second `make` does not stop on "remote already exists".
