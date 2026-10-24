# Phase 3b: column physics skips the inactive columns

Phase 3 ran microphysics on every owned column on every level step and copied the off-level result back. On `sp8.24755` that raised the Kessler call count from 120 to 225 (`1 + 1/2 + 1/4 + 1/8 = 1.875`) and left `time integration` at 1.029× all-fine ([PHASE3.md](PHASE3.md)). This phase skips those columns before the work, and advances each column with the `dt` of the level that owns it.

The binary for every number below is md5 `c1bb56d80d0b58f9c88bf62b86333e0a`. Four MPI ranks, OpenMP off, `--oversubscribe`, single precision, the same `dt = 30` s Kessler namelist and the same level-balanced partition as phase 3. `config_hts` still defaults to false. The acoustic step was not lengthened.

## What changed

`physics_timetracker` still runs once per fine clock tick, so the radiation and convection alarms stay on the fine clock. On one level, and with HTS off, `physics_driver` is still that one full-block call.

On a multi-level step the driver moves inside the level loop in `atm_hts_advance`, after the interpolation and the wind refresh and before `atm_srk3`. `mpas_atmphys_column_mask` holds `column_mask_on` and `column_active(:)`. The active set is the owned cells of the active level (`hts_active_own`), the same cells whose state this step is allowed to update. It is not the two-ring dynamics work list: a column scheme has no horizontal stencil, and writing the ring would advance a neighbor level with the wrong `dt`. `column_mask_next_segment` walks that mask as contiguous `its:ite` runs. With the mask off it returns the original range in one call.

`dt_dyn` and `dt_microp` are pushed to the level `dt` for the driver and popped before the dynamics RK. An alarm that rings on a tick when a coarser level does not step is latched in `hts_radlw_due`, `hts_radsw_due`, and `hts_conv_due` until that level steps, then cleared. The delay is at most one coarse interval.

Microphysics no longer snapshots the full state and restores the off-level columns. The copy loops in `microphysics_from_MPAS`, `microphysics_to_MPAS`, `precip_from_MPAS`, `precip_to_MPAS`, and `compute_relhum`, and both Kessler `i` loops, cycle a cell when the mask is on and that cell is inactive. The scratch arrays stay allocated across masked calls. `hts_restrict_physics_tendencies` is unchanged and still zeros `tend_ru_physics`, `tend_rtheta_physics`, and `tend_rho_physics` off the active level. Diagnostics stay full-domain.

The mask check is a nested `if`, so `column_active` is not read when the mask is off. One level returns before the level lists are built, so this path is not entered and the driver call is the original block.

## Kessler, 1 hour, sp8.24755

Suite `mesoscale_reference` with only `mp_kessler` left on. Scalar advection off. 120 fine steps. IEEE flags are underflow and denormal. `u` stays about ±35 m/s. The log line is `column physics runs on the active level; microphysics skips the other columns`. One non-fatal error is the missing `diagnostics` stream, the same line as phase 3.

| | all-fine | HTS, max level 3 |
| --- | ---: | ---: |
| `/usr/bin/time` wall | 25.18 s | 20.84 s |
| MPAS total | 24.840 s | 20.500 s |
| `time integration` | 23.349 s / 120 | 19.024 s / 120 |
| `physics driver` | 0.00067 s / 120 | 0.00135 s / 225 |
| `atm_rk_integration_setup` | 0.226 s / 120 | 0.061 s / 225 |
| `atm_compute_moist_coefficients` | 0.224 s / 120 | 0.433 s / 225 |
| `physics_get_tend` | 0.115 s / 120 | 0.237 s / 225 |
| `atm_compute_vert_imp_coefs` | 0.189 s / 120 | 0.381 s / 225 |
| `atm_compute_dyn_tend` | 9.616 s / 360 | 2.152 s / 675 |
| `small_step_prep` | 0.435 s / 360 | 0.797 s / 675 |
| `atm_advance_acoustic_step` | 4.622 s / 1440 | 5.042 s / 5760 |
| `atm_divergence_damping_3d` | 1.393 s / 1440 | 1.512 s / 5760 |
| `atm_recover_large_step_variables` | 0.907 s / 360 | 0.234 s / 675 |
| `atm_compute_solve_diagnostics` | 1.554 s / 360 | 2.780 s / 675 |
| `atm_rk_dynamics_substep_finish` | 0.134 s / 120 | 0.046 s / 225 |
| `microphysics` | 2.171 s / 120 | 0.529 s / 225 |
| `mp_kessler` | 0.908 s / 120 | 0.165 s / 225 |
| gap after every child of `time integration` | 1.76 s | 4.82 s |

Integration ratio 23.349 / 19.024 = **1.227×**. Wall ratio 25.18 / 20.84 = 1.208×. Phase 3's restore coupling on the same case was 23.717 / 23.044 = 1.029×. The 1.15× bar on `time integration` is met.

`mp_kessler` is 0.90828 / 0.16530 = 5.49×. The cell-count ceiling for four levels on this mesh is 6.10×: active-cell work per fine step is `(885 + 452/2 + 150/4 + 23268/8) / 24755 = 0.164` of an all-fine sweep. The kernel is most of the way there. The `microphysics` wrapper is 4.10× (2.171 / 0.529). It still enters 225 times and the copy loops still walk the index to decide which cells to skip.

The rest of the 4.3 s is the dry tendency saving, partly given back. `atm_compute_dyn_tend` fell 9.62 → 2.15 s. Diagnostics, moist coefficients, `physics_get_tend`, and the vertical-implicit coefficients still walk every owned column and are entered once per level step, so they grow by about the 1.875 call factor (diagnostics 1.55 → 2.78 s). Acoustic time stays near 5 s while the call count goes from 1440 to 5760, because that loop was already on the per-level work list. The untimed gap, halo exchanges plus the HTS driver, grew from 1.76 s to 4.82 s. The net is 1.23×, not the 6.10× tendency ceiling and not a 6× physics ceiling. `physics driver` itself is 0.001 s: with convection, PBL, radiation, the surface layer, and the land model off, the driver returns before those schemes.

## Accuracy

Last frame, HTS against the Kessler all-fine history of this binary. The same frame is bit-identical to the phase 3 Kessler histories (`history.phase3.nc` in each run directory, not in git): relative L2 of `u`, `w`, `theta`, `rho`, and near-surface pressure is 0, and the dry-mass frame difference is 0. Skipping the inactive columns leaves the same state the restore used to keep.

| field | rel L2 | Linf |
| --- | ---: | ---: |
| `u` | 6.159×10⁻⁴ | 0.185 m/s |
| `w` | 0.293 | 6.56×10⁻³ m/s |
| `theta` | 8.409×10⁻⁶ | 0.140 K |
| `rho` | 1.282×10⁻⁵ | 8.345×10⁻⁵ kg/m³ |
| near-surface `pressure` | 1.586×10⁻⁵ | 9.98 Pa (rms 1.56 Pa) |

Dry-mass drift −2.862×10⁻¹¹ (HTS) and −9.604×10⁻¹¹ (all-fine). Frame difference +6.741×10⁻¹¹. The initial sounding is dry, so Kessler has little water to convert. The test is that the per-level `dt` and the skipped columns do not move the mass. They do not, and they match the phase 3 frame exactly.

## One level stays the all-fine history

Uniform 120 km, physics suite `none`, 10 steps, `dt = 450` s. HTS on against HTS off, and each against the phase 2b histories `hts_off_120km_10step` and `hts_on_120km_10step`: relative L2 of `u`, `w`, `theta`, and `rho` is 0, dry-mass frame difference 0, drift −1.977×10⁻¹⁰. The log line is `HTS enabled with a single level; masked integrator, no subcycling`. The column-physics line is absent, because `hts_max_level_used` is 0.

`x1.40962`, CFSR `2010-10-23_00`, suite `mesoscale_reference`, `dt = 720` s, 2 steps (24 minutes). HTS on against HTS off: relative L2 is 0 on `u`, `w`, `theta`, and `rho`, dry-mass frame difference 0, drift −1.799×10⁻¹¹. `physics driver` is called twice, once per step. Radiation fired once (`rrtmg_swrad` and `rrtmg_lwrad` each have one call), and WSM6, YSU, Noah, and Tiedtke ran on both steps. Walls 19.53 s (off) and 19.98 s (on) on a bit-identical pair. The 6-hour forecast was not repeated on this binary. Its bit-identical result in [PHASE3.md](PHASE3.md) is from the phase 3 binary, which took the same one-level branch.

## What is still open

- Thompson and WSM6 are not masked inside the WRF kernel. A multi-level step with the column mask on and any scheme other than `mp_kessler` aborts in `driver_microphysics` with `HTS column mask is implemented for mp_kessler only`. Skipping the copies and then calling an unmasked kernel would read uninitialized scratch. Kessler is the scheme this phase timed.
- Radiation, PBL, the land model, convection, the surface layer, GWDO, and sea ice are segmented onto the owned active cells and use the level `dt`, with the radiation and convection alarms latched until that level steps. That path was not timed on a variable-resolution mesh. There is no `sp8.24755` static here, and the Kessler namelist turns those schemes off. A level's cells are scattered, so one RRTMG tile can become many short calls. That overhead is unknown. `driver_oml1d` is still one full-domain call. `compute_radar_reflectivity` is not masked; Kessler does not call it.
- Diagnostics, `physics_get_tend`, and the moist-coefficient loops still walk every owned column once per level step. Together with the halo gap they are why the net is 1.23×.
- The acoustic step is still the fine step on every level. `dt = 60` s was not retested. The phase 2c overflow at that `dt` still stands.
- GWDO reads neighboring columns. The mask uses the owned active list, so a neighbor on another level is not updated on this step. That dependence was not tested.
- OpenACC was not built. The mask is host-side.
