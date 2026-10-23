# Phase 3: physics on the level that owns the column

Phase 2c's four-level dry run on `sp8.24755` was 1.221× faster than all-fine on `time integration` ([PHASE2C.md](PHASE2C.md)). That was the gate for physics. The coupling is in `atm_srk3`. `config_hts` still defaults to false.

The binary for every number in this file is md5 `3a4a65d0b6f7504cf7d445fa3578ed36`. The one-level dry gate was re-run on it: 10 steps, uniform 120 km, HTS on and HTS off, against each other and against the phase 2b histories. Relative L2 of `u`, `w`, `theta`, and `rho` is 0, and the dry-mass frame difference is 0. Drift −1.977×10⁻¹⁰.

## What the coupling does

`atm_do_timestep` still calls `physics_timetracker` and `physics_driver` once per fine clock tick, on the whole block, before `atm_hts_advance`. Radiation alarms stay on `config_dt`. That call was not multiplied by the number of levels.

`physics_get_tend` still runs at the start of every `atm_srk3`. On a multi-level step, `hts_restrict_physics_tendencies` then zeros `tend_ru_physics`, `tend_rtheta_physics`, and `tend_rho_physics` off the active level. A coarse column is not forced on a fine tick, and a fine column is not integrated again over `dt * 2^L`. The rate a coarse step uses is the one `physics_driver` wrote on the latest fine tick. All-fine and a single level return before any of those zeros.

Microphysics sits at the end of `atm_srk3` and updates the state, so a zeroed tendency is not enough. On a multi-level step the call uses that level's `dt` (`dt_dyn` and `dt_microp` are pushed and popped; Thompson keeps the 90 s split it uses at initialization). Columns that are not on the active level are copied back, including `theta_m`, the scalars, `pressure_p`, `exner`, `rtheta_p`, the diabatic rate, surface pressure, and the precipitation accumulators. The driver itself still walks every owned column and the restore throws the off-level writes away. Call counts match that: Kessler microphysics ran 120 times all-fine and 225 times at four levels, and `1 + 1/2 + 1/4 + 1/8 = 1.875`.

`docs/hts/scripts/hts_partition.py --physics-weight N` adds a uniform per-cell physics cost on every METIS constraint and keeps the dynamics weight on the cell's own level. `N = 0` is the phase 2b graph. It was not used for the runs below. On the uniform real mesh every cell is one level, so a physics weight does not move the partition. On `sp8.24755` the Kessler cost is paid on every column, not only on the active level, and `dyn_tend` parallel efficiency was already 0.96. A uniform addend does not remove the extra microphysics calls.

## Real data, uniform 120 km

Ready-made inputs, both already on this machine:

- Static: `/home/ubuntu/mpas-data/x1.40962.static.nc`, 73 MB, file date 2019-06-08, `NETCDF3_64BIT_OFFSET`, 40962 cells, terrain from −130 m to 5380 m, land use `MODIFIED_IGBP_MODIS_NOAH`. This is the standard 120 km static. It is not in git.
- First guess: `cfsr.2010102300.tar.gz`, member `CFSR:2010-10-23_00`, member date 2014-05-28, 230 MB uncompressed. The default init namelist's start time is that date.

`init_atmosphere` case 7, `config_static_interp = false`, `config_native_gwd_static = false`, `config_met_interp = true`, `config_vertical_grid = true`, 55 levels, start `2010-10-23_00:00:00`. This MPAS version's met-interp input stream also reads `soilcomp`, `soilcl*`, and the vertical-metric fields. The 2019 static does not have them (`nSoilComps` is absent, and the reader aborts). A local copy adds those fields as zeros, with `soilcomp` component 0 set to 100 so the array is not empty. `sf_noah` reads `isltyp` from the original static, not `soilcomp` (that array is for Noah-MP). The vertical stage recomputes `zgrid` and `zz`. Init was one rank, 11 s, and wrote `x1.40962.init.nc` (356 MB, not in git). `u` is −90 to +89 m/s, `qv` up to 0.022 kg/kg, with a few −8×10⁻⁵ kg/kg undershoots from the interpolator.

Forecast namelist: the v7 sample in `namelist_files.tar.gz`. `config_dt = 720` s, `config_len_disp = 120000`, physics suite `mesoscale_reference` (`mp_wsm6`, `cu_ntiedtke`, `bl_ysu`, `bl_ysu_gwdo`, RRTMG, Noah), scalar advection on, split dynamics and transport on, 4 ranks, the stock `x1.40962.graph.info.part.4`. Duration 6 hours (30 steps). `qi`, `qs`, and `qg` are absent from the init file and start at zero; the log warns and continues. There is no variable-resolution static for `sp8.24755` here (`mpas_static` was not downloaded), so this case does not exercise subcycling. The mesh is uniform, the log line is `HTS enabled with a single level; masked integrator, no subcycling`, and the new zeroing and restore do not run.

Against all-fine of the same binary, 6-hour frame: relative L2 of `u`, `w`, `theta`, and `rho` is 0. Dry-mass drift is −4.31×10⁻¹¹ for both, frame difference 0. Final `u` is about −91 to +92 m/s. No NaN lines.

| | all-fine | HTS, one level |
| --- | ---: | ---: |
| `/usr/bin/time` wall | 229.11 s | 221.68 s |
| MPAS total | 228.71 s | 221.25 s |
| `time integration` | 226.16 s / 30 | 218.99 s / 30 |
| `physics driver` | 133.21 s / 30 | 130.48 s / 30 |
| `rrtmg_swrad` | 88.64 s / 12 | 87.75 s / 12 |
| `rrtmg_lwrad` | 26.90 s / 12 | 26.69 s / 12 |
| `atm_compute_dyn_tend` | 22.45 s / 270 | 22.16 s / 270 |
| `microphysics` | 37.23 s / 30 | 35.67 s / 30 |

The histories match bit for bit, so the 3% wall gap is not a scheme effect. Shortwave radiation is 0.54 parallel-efficient on this 4-rank oversubscribe, and that timer moves by about a second between the two runs. Physics is 133 s of the 226 s integration (59%). Microphysics is another 37 s (16%). The slow tendency is 22 s (10%). On a variable mesh the tendency is the only one of those three that HTS shrinks.

## Kessler on the small patch

`sp8.24755`, JW case 2, `qv = qc = qr = 0` at init, suite `mesoscale_reference` with only `mp_kessler` left on (convection, PBL, gravity-wave drag, radiation, surface layer, and the land model off). Scalar advection off, `dt = 30` s, the same 4-rank level-balanced partition as phase 2c, 1 hour (120 steps). The ozone tables are read at physics init even with radiation off; they were linked into the run directory. A missing `diagnostics` stream prints one error line and does not abort.

| | all-fine | HTS, max level 3 |
| --- | ---: | ---: |
| wall | 25.70 s | 24.82 s |
| `time integration` | 23.717 s / 120 | 23.044 s / 120 |
| `atm_compute_dyn_tend` | 9.837 s / 360 | 2.219 s / 675 |
| `atm_compute_solve_diagnostics` | 1.757 s / 360 | 2.886 s / 675 |
| `microphysics` | 2.356 s / 120 | 4.087 s / 225 |
| `mp_kessler` | 0.962 s / 120 | 1.657 s / 225 |

Integration ratio 23.717 / 23.044 = **1.029×**. The tendency saving is the dry one (9.84 → 2.22 s). Microphysics grew by the 1.875 call factor and put 1.7 s back. Diagnostics grew with the RK count. The net is a few percent, not the dry 1.22×.

Accuracy against the Kessler all-fine history, 1 hour. `u` stayed ±35 m/s.

| field | rel L2 | Linf |
| --- | ---: | ---: |
| `u` | 6.16×10⁻⁴ | 0.185 m/s |
| `w` | 0.293 | 6.56×10⁻³ m/s |
| `theta` | 8.41×10⁻⁶ | 0.140 K |
| `rho` | 1.28×10⁻⁵ | 8.34×10⁻⁵ kg/m³ |

Dry-mass drift −2.86×10⁻¹¹ (HTS) and −9.60×10⁻¹¹ (all-fine). Frame difference +6.7×10⁻¹¹. The initial sounding is dry, so Kessler has little water to convert. The test is that the save/restore and the per-level `dt` do not move the mass or blow up. They do not.

## What this says about speed

A full mesoscale step on the 120 km mesh spends about three quarters of `time integration` in `physics driver` and microphysics. Both are column-local. The driver already runs once per fine tick, which is what the alarms require, and HTS does not reduce it. Microphysics, because it updates the state, runs on each level step and currently walks every column on that step. On `sp8.24755` that turned a 1.22× dry integration speedup into a 1.03× Kessler speedup. A variable-resolution mesh with RRTMG on top of that would be slower than all-fine unless microphysics is limited to the active columns before the call, not after it. That limit was not added: the microphysics drivers take a contiguous `its:ite` range, and a level is not a contiguous range.
