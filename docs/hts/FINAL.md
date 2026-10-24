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
| `sp8.24755`, 30 min, 4 levels, WSM6 | 3.63×10⁻⁴ | 0.486 | 3.67×10⁻⁶ | 1.03×10⁻⁵ | −7.93×10⁻⁹ | −1.31×10⁻¹⁰ |
| `sp8.24755`, 15 min, 4 levels, Thompson | 2.12×10⁻⁴ | 0.300 | 2.34×10⁻⁶ | 6.20×10⁻⁶ | +4.03×10⁻⁹ | −8.70×10⁻¹¹ |
| `sp8.24755`, 8 min, 4 levels, RRTMG + YSU + Noah + Kessler | 3.66×10⁻⁴ | 0.267 | 2.24×10⁻⁵ | 1.67×10⁻⁴ | −1.47×10⁻¹⁰ | −5.28×10⁻¹¹ |
| `x1.40962` real, 6 h, 1 level, `mesoscale_reference` | 0 | 0 | 0 | 0 | −4.31×10⁻¹¹ | −4.31×10⁻¹¹ |

`w` relative L2 is large because the reference `w` is small (about 10⁻⁴ to 10⁻³ m/s on the dry JW runs). The absolute rms difference on `sp8.24755` at 6 hours is 4.6×10⁻⁴ m/s. The largest piece sits on the interface. `u`, `theta`, and `rho` stay at relative 10⁻⁴, 10⁻⁵, and 10⁻⁶. Dry-mass drift on the fast four-level dry run is 1.27×10⁻⁹. The previous full-domain-diagnostic binary on that same frame was 6.96×10⁻¹⁰, and `u` relative L2 between those two HTS histories is 1.21×10⁻⁵.

`dt = 60` s on `sp8.24755` is omitted from the table. All-fine finished. Every HTS clamp signaled invalid and overflow. The acoustic step stays the fine step. The comparisons below use `dt = 30` s. The CFL numbers are in the next section.

## Speed against all-fine

Ratio is all-fine / HTS on `time integration`. Higher than 1 is faster. The `sp8.24755` rows are 4 MPI ranks, OpenMP off, `dt = 30` s, on the diagnostic-ring binary. Earlier binaries are labeled.

| case | all-fine | HTS | ratio |
| --- | ---: | ---: | ---: |
| `x5.30210`, 6 h, 2 levels, 55% of cells fine | 25.06 s | 27.34 s | 0.92× |
| `x4.10242`, 6 h, 3 levels, 49% of cells fine | 3.62 s | 4.25 s | 0.85× |
| `sp8.24755`, 6 h, dry, phase 2c binary | 127.10 s | 104.05 s | 1.221× |
| `sp8.24755`, 6 h, dry | 120.66 s | 80.66 s | **1.496×** |
| `sp8.24755`, 1 h, Kessler, phase 3 restore | 23.72 s | 23.04 s | 1.029× |
| `sp8.24755`, 1 h, Kessler, phase 3b column mask | 23.35 s | 19.02 s | 1.227× |
| `sp8.24755`, 1 h, Kessler | 22.98 s | 14.35 s | **1.601×** |
| `sp8.24755`, 30 min, WSM6, phase 3b binary | 15.91 s | 10.65 s | 1.494× |
| `sp8.24755`, 30 min, WSM6 | 16.15 s | 7.69 s | **2.100×** |
| `sp8.24755`, 15 min, Thompson | 7.32 s | 3.71 s | **1.971×** |
| `sp8.24755`, 8 min, RRTMG + YSU + Noah + Kessler | 6.49 s | 4.81 s | **1.349×** |
| `x1.40962` real, 6 h, 1 level, full suite, phase 3 binary | 226.16 s | 218.99 s | bit-identical; the 3% is timer scatter |

Wall (`/usr/bin/time`) on the new `sp8` rows: dry 125.47 s / 85.62 s, Kessler 25.34 s / 16.10 s, WSM6 17.53 s / 9.08 s, Thompson 9.12 s / 5.23 s, radiation 7.58 s / 5.95 s.

## What the timers paid for

Level 0 is 885 of 24755 cells. The spacing ratio is 10×. The large-step tendency ceiling is 6.10×. On the 6-hour dry pair the loops that were walking the whole domain, and that the timers showed growing under HTS, are now on a list:

| 6 h dry | all-fine | HTS | calls |
| --- | ---: | ---: | ---: |
| `atm_compute_dyn_tend` | 56.92 s | 12.14 s | 2160 → 4050 |
| `atm_compute_solve_diagnostics` | 10.07 s | 2.37 s | 2160 → 4050 |
| `atm_compute_moist_coefficients` | 0.845 s | 0.181 s | 720 → 1350 |
| `atm_compute_vert_imp_coefs` | 1.102 s | 0.200 s | 720 → 1350 |
| `small_step_prep` | 2.775 s | 0.527 s | 2160 → 4050 |
| `atm_advance_acoustic_step` | 27.15 s | 29.42 s | 8640 → 34560 |
| `physics_get_tend` | 0.542 s | 0.670 s | 720 → 1350 |

Diagnostics use four rings (rank 0, level 3: 6127 cells of the local domain; levels 0–2 are a few hundred cells). The previous HTS binary spent 16.56 s in diagnostics on this run. Moist coefficients and the vertical implicit coefficients follow the dynamics work list and the active columns. `physics_get_tend` still zeros the full scalar-tendency array and then fills it. Its 6-hour cost is 0.67 s, so that zero stays full-domain: an off-level leftover would persist into the next level.

The acoustic timer stays near 29 s while the call count quadruples, because the acoustic loop was already on the per-level list and the substep count grows by `2^L` to keep `dts = 5` s. The untied remainder of `time integration`, after every direct child timer, is 5.39 s all-fine and 24.34 s with four levels. Those seconds are the halo exchanges and the HTS driver. Every rank enters every collective on the existing full-domain groups. After the level-balanced partition every rank owns every level, so a per-level group would still include every rank.

`mp_kessler` on the hour is 1.045 s / 0.158 s = 6.61×. `mp_wsm6` on the half hour is 3.734 s / 0.646 s = 5.78×. Both sit on the 6.10× cell-count ceiling. `mp_thompson` on 15 minutes is 1.060 s / 0.350 s = 3.03×. Thompson keeps its 90 s internal substep, and a 15-minute run is short relative to that split.

## Thompson

`build_tables` (`src/core_atmosphere/utils/build_tables.F`) wrote the four lookup files in 176.5 s:

- `MP_THOMPSON_QRacrQG_DATA.DBL`
- `MP_THOMPSON_QRacrQS_DATA.DBL`
- `MP_THOMPSON_freezeH2O_DATA.DBL`
- `MP_THOMPSON_QIautQS_DATA.DBL`

They live outside the repo. The 15-minute smoke is `mp_thompson` with convection, PBL, radiation, the surface layer, and the land model off, the same shape as the Kessler namelist. `qi`, `qs`, and `qg` start at zero. IEEE flags are underflow and denormal. `u` stays about ±35 m/s. Integration is **1.971×**. Dry-mass frame difference is +4.12×10⁻⁹.

## Radiation and PBL

The JW init for `sp8.24755` has blank land use, `ivgtyp = 0`, `isltyp = 0`, `xland = 0`, `ter = 0`, and `smois = 0`. A copy used for this smoke sets USGS water (`ivgtyp = 16`, `isltyp = 14`, `xland = 2`, `mminlu = 'USGS'`), which matches `iswater_lu = 16` on that file. GWDO stays off because `var2d` is zero. Convection stays off. The schemes that run are `rrtmg_lw`, `rrtmg_sw`, `cld_fraction`, `sf_monin_obukhov_rev`, `sf_noah`, `bl_ysu`, and `mp_kessler`. Radiation interval 30 minutes. Duration 8 minutes, so the t = 0 alarm is the only radiation firing. It is latched until each level steps, then cleared. Shortwave time is 1.585 s in one all-fine call and 1.600 s across the masked segments. `bl_ysu` is 0.570 s / 0.068 s = 8.36×.

The lowest-layer `theta` on levels 0–2 matches all-fine to about 0.001 K. On level 3 the same layer warms by a mean of 0.266 K (HTS) and 0.299 K (all-fine). The mean difference is −0.034 K and the maximum is 0.38 K, on cells thousands of kilometres from the fine patch. Both runs warm that layer because the filled water static has zero soil moisture. Setting `dt_pbl` to the level `dt` (240 s on level 3) raised that mean difference to 0.093 K on the same pair, so `dt_pbl` stays at the fine step stored at initialization. YSU and Noah return a rate normalized by that `dt`, and the level RK multiplies the stored rate by `dt_level`.

## Safe step

`config_dt = 30` s, acoustic `dts = 5` s, `config_number_of_sub_steps = 6`. A 6-hour attempt at `config_dt = 60` s (`dts = 10` s) is written up in [PHASE2C.md](PHASE2C.md). All-fine finished (wall 64.47 s, underflow and denormal). Every HTS clamp signaled `IEEE_INVALID` and `IEEE_OVERFLOW`. With `u ~ 35` m/s and `dt_level = 60 * 2^L`, the advective CFL on the shortest edge of each level stays under 0.21. The fine-edge acoustic CFL at `dts = 10` s is about 0.22. The overflow is the coarse and interface acoustic update. The acoustic step was left at the fine step.

## Gates

`config_hts` defaults to false in `Registry.xml`. One level returns before the lists are built, and the masked loops require `hts_max_level_used > 0`, so they keep the original index order.

Uniform 120 km, physics suite `none`, 10 steps, `dt = 450` s. HTS on, HTS off, and each against its previous history: relative L2 of `u`, `w`, `theta`, and `rho` is 0. Drift −1.977×10⁻¹⁰.

`x1.40962`, CFSR `2010-10-23_00`, suite `mesoscale_reference`, `dt = 720` s, 2 steps. HTS on against HTS off, and each against the previous smoke history: relative L2 0 on `u`, `w`, `theta`, `rho`, and `pressure`. Drift −1.799×10⁻¹¹. The 6-hour forecast that matches bit for bit was run on the phase 3 binary.

## Bottom line

On `sp8.24755` at `dt = 30` s, four levels, this binary is **1.496×** dry (6 h), **1.601×** Kessler (1 h), **2.100×** WSM6 (30 min), **1.971×** Thompson (15 min), and **1.349×** for an 8-minute RRTMG/YSU/Noah/Kessler smoke. The one-level and HTS-off histories match bit for bit. The acoustic step stays the fine step, because `dt = 60` s overflows on this mesh while the linear CFL stays small. The seconds that remain are the extra halo calls around that fine acoustic step: 24.3 s of the 80.7 s dry integration. A land-surface variable-resolution forecast is out of reach on this machine. The JW static for this mesh has blank land use and zero orography, so the radiation smoke used a USGS-water fill and left GWDO and convection off.
