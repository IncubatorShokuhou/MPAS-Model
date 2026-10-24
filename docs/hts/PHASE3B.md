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

## Schemes beyond Kessler

The Kessler-only stop is gone. `column_is_active` is the shared test: with the mask off it is true and does not read `column_active`. WSM6, Thompson, and Thompson-aerosols are column-local, so `driver_microphysics` calls them through `column_mask_next_segment`. With the mask off that is one call on the original `its:ite`. Radar reflectivity uses the same segments. The 1D ocean mixed layer, the Noah/Noah-MP 2 m diagnostics, the skin SST, and the deep-soil temperature skip inactive columns. Skin SST and deep-soil temperature run from the level loop, with that level's `dt`, only when `hts_max_level_used > 0`. Both namelist flags default to false, so the one-level time tracker still makes the original calls.

Radiation, PBL, the land model, convection, the surface layer, GWDO, sea ice, and cloudiness were already called on those segments from `physics_driver`. No microphysics scheme in this source aborts because the mask is on.

Ozone on the radiation alarm, the monthly albedo, and the SST/sea-ice file update stay on the fine clock for every owned column. They are not integrated with the level `dt`. A coarser level's later radiation reads the ozone written on the alarm tick. `physics_get_tend`, the moist-coefficient loops, and the diagnostics stay full-domain.

Binary md5 `fd5bfbe02ba0a5cc4705f48689bfeab4`. The Kessler hour on it matches the history behind the table above bit for bit (relative L2 0 on `u`, `w`, `theta`, `rho`, and near-surface pressure). A repeat of that hour timed `time integration` at 22.912 s / 18.820 s = **1.217×**, wall 24.96 s / 20.60 s. The earlier 1.227× is the same state; the tenth of a second is timer scatter. The 10-step dry gate and the 2-step `x1.40962` `mesoscale_reference` smoke match the previous histories of those runs bit for bit.

WSM6 on the same mesh, `dt = 30` s, 30 minutes (60 steps), 4 ranks, the same partition. `qi`, `qs`, and `qg` are absent from the JW init and start at zero. Convection, PBL, radiation, the surface layer, and the land model stay off, as in the Kessler namelist. IEEE flags are underflow and denormal. `u` stays about ±35 m/s.

| | all-fine | HTS, max level 3 |
| --- | ---: | ---: |
| `/usr/bin/time` wall | 17.33 s | 12.01 s |
| MPAS total | 17.011 s | 11.677 s |
| `time integration` | 15.914 s / 60 | 10.654 s / 60 |
| `atm_compute_dyn_tend` | 4.966 s / 180 | 1.088 s / 336 |
| `microphysics` | 5.212 s / 60 | 1.065 s / 112 |
| `mp_wsm6` | 3.592 s / 60 | 0.662 s / 112 |

Integration ratio 15.914 / 10.654 = **1.494×**. Wall 17.33 / 12.01 = 1.443×. `mp_wsm6` is 3.592 / 0.662 = 5.43×, again near the 6.10× cell-count ceiling. The wrapper is 4.90×. The call count is 60 against 112, the same 1.875 level factor as the Kessler hour. WSM6 is a larger share of the step than Kessler, so skipping the coarse columns moves the integration ratio from 1.22× to 1.49× on this half hour.

Last frame, HTS against the WSM6 all-fine history.

| field | rel L2 | Linf |
| --- | ---: | ---: |
| `u` | 3.626×10⁻⁴ | 0.132 m/s |
| `w` | 0.486 | 8.58×10⁻³ m/s |
| `theta` | 3.671×10⁻⁶ | 0.136 K |
| `rho` | 1.032×10⁻⁵ | 1.17×10⁻⁴ kg/m³ |
| near-surface `pressure` | 1.267×10⁻⁵ | 13.3 Pa (rms 1.25 Pa) |

Dry-mass drift −7.88×10⁻⁹ (HTS) and −1.31×10⁻¹⁰ (all-fine). Frame difference −7.75×10⁻⁹. That is larger than the Kessler hour. The initial sounding is dry, WSM6 still condenses, and the level `dt` is not the fine `dt`. It is not a NaN, and `u` does not grow.

Thompson and Thompson-aerosols take the same segmented `mp_gt_driver` call. They were not run: `MP_THOMPSON_QRacrQG_DATA.DBL` and the other Thompson lookup files are not on this machine. A variable-resolution hour with RRTMG and PBL was not run. There is no `sp8.24755` static here, and a level's cells are scattered, so one RRTMG tile becomes many short calls. That overhead is still unknown. GWDO's static orography is per column; a neighbor on another level is not updated on this step, and that dependence was not tested. The acoustic step was not lengthened. OpenACC was not built.
