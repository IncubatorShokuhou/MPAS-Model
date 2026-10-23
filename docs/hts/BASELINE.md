# Phase 1 baseline: unmodified MPAS-Atmosphere v8.4.2

Dry Jablonowski–Williamson baroclinic wave, init case 2 (perturbation on),
26 levels, physics suite `none`, scalar advection off. The binary was built
from tag `v8.4.2` (`f34984b2`, `git describe` reports `v8.4.2`) before any
HTS documentation commits. No model source was changed.

Raw NetCDF (`init.nc`, `history.nc`) is not in git. A 10-day history is
456 MB (uniform) and 337 MB (variable). Reproduce with
`docs/hts/scripts/run_jw_baseline.sh both 10_00:00:00`. The runs below are
in `$HOME/mpas-runs/jw_uniform_120km` and
`$HOME/mpas-runs/jw_variable_240_48km` on the machine that produced them.
`docs/hts/scripts/summarize_baseline.py` recomputes mass from those files.

## Build

| item | value |
| --- | --- |
| compiler | gfortran 13.3.0, `make -j2 gfortran` |
| MPI | Open MPI 4.1.6 |
| NetCDF | C 4.9.2, Fortran 4.5.4 (Ubuntu packages) |
| PnetCDF | 1.12.3 (Ubuntu package) |
| I/O library | SMIOL (bundled). `PIO` was not set |
| precision | single (MPAS default) |
| OpenMP / OpenACC | off |
| flags | `-O3 -ffree-line-length-none -fconvert=big-endian -ffree-form` |

`docs/hts/scripts/build_mpas.sh` recreates the symlink prefix the Makefile
expects (`$HOME/mpas-deps/{netcdf,pnetcdf}/{include,lib}`) and builds
`atmosphere_model` then `init_atmosphere_model`.

## Hardware and launch

Cloud agent VM, 4 cores, about 16 GB RAM. `mpirun -np 4`, `OMP_NUM_THREADS=1`.
`/usr/bin/time` `%M` is the max RSS of the `mpirun` process it timed, which
matched a single rank, not the sum of ranks. Per-rank RSS was sampled from
`/proc` during the run.

## Case setup

Both cases use the v7.0 JW namelist dynamics, with `config_epssm` omitted.
v8.4.2 ignores that option and warns. Constant off-centering of 0.1 is set
with `config_epssm_minimum = config_epssm_maximum = 0.1`. `config_xnutr = 0`
as in the published namelist (no upper-level w damping). History is written
daily, including the initial time, so each `history.nc` has 11 frames
(day 0 through day 10).

`surface_pressure` in the history file stays at 100000 Pa. It is not the
evolving surface pressure. The numbers below use the lowest-level
`pressure` (vertical index 0, which is the near-surface level: day-0 values
are about 97900–98506 Pa).

| | uniform | variable |
| --- | --- | --- |
| mesh | `x1.40962` | `x5.30210` |
| cells | 40962 | 30210 |
| spacing | 120 km quasi-uniform | 240 km to 48 km, ratio 5, patch centered at 0°N, 0°E |
| `config_dt` | 450 s | 180 s |
| `config_len_disp` | 120000 m | 48000 m |
| steps in 10 days | 1920 | 4800 |
| acoustic substeps | 6 (`dt_acoust = dt/6`) | 6 |
| `meshScalingDel4` min/max | 1.000 / 1.000 | 1.000 / 121.782 |

The del4 max of 121.8 on the variable mesh is `(dx/dx_min)^3` at the coarse
end (5³ = 125). `config_len_disp = 48000` is the right finest-spacing choice.
Init warned that the 2019 `x1.40962` grid has no `nominalMinDc`. The 2023
`x5.30210` grid has it. Both inits finished. Atmosphere logs reported 0
errors and 0 warnings.

## Timings

| | uniform 120 km | variable 240–48 km |
| --- | ---: | ---: |
| init wall | 2.94 s | 2.08 s |
| init max RSS (`/usr/bin/time`) | 547 MB | 429 MB |
| model wall, 10 days, 4 ranks | 535.27 s | 930.31 s |
| model user time (all ranks) | 2133.46 s | 3712.61 s |
| integration timer | 522.19 s over 1920 calls | 903.03 s over 4800 calls |
| mean step | 0.272 s (min 0.258, max 0.458) | 0.188 s (min 0.176, max 0.454) |
| wall per simulated day | 53.5 s | 93.0 s |
| sampled RSS per rank | 439–475 MB | 324–355 MB |
| sampled RSS, 4 ranks | about 1.80 GB | about 1.33 GB |

Integration is 98% of model wall time. A 16-day run matching the published
namelist duration would be about 1.6× these wall times (about 14 min and
25 min) on this VM.

## Mass

Dry mass is `sum(rho * areaCell * dz)` with `dz` from `zgrid`. Global
atmospheric mass is 5.20463878×10¹⁸ kg (uniform) and 5.20463862×10¹⁸ kg
(variable). `sum(surface_pressure * areaCell) / g` with `g = 9.80616` is
5.20184382×10¹⁸ kg on both, and it does not move because that field is not
updated.

Relative drift of the `rho` integral from day 0 to day 10:

| day | uniform | variable |
| ---: | ---: | ---: |
| 1 | −6.0×10⁻¹¹ | −2.0×10⁻¹⁰ |
| 5 | −4.6×10⁻¹¹ | +2.8×10⁻¹¹ |
| 10 | −2.1×10⁻¹⁰ | −2.8×10⁻¹⁰ |

The drift stays a few parts in 10¹⁰ over 10 days. That is the number HTS
has to match, not a looser single-precision budget. It is still a
single-precision run; the global sum is just very well conserved by the
flux-form integrator.

## Wave growth

Lowest-level pressure anomaly against day 0. The perturbation is small
through day 4 and then grows, which is the JW wave.

| day | uniform rms (Pa) | uniform max \|Δp\| (Pa) | variable rms (Pa) | variable max \|Δp\| (Pa) |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 2.8 | 33 | 3.1 | 30 |
| 4 | 10.2 | 145 | 6.4 | 131 |
| 7 | 78.4 | 1128 | 44.0 | 958 |
| 9 | 275 | 5085 | 158 | 4344 |
| 10 | 401 | 6257 | 239 | 6110 |

Day-10 extrema of normal velocity `u` (m/s), from the model log:

| | day 0 min/max | day 10 min/max |
| --- | --- | --- |
| uniform | −35.14 / 35.16 | −52.14 / 51.94 |
| variable | −34.99 / 34.99 | −48.46 / 50.22 |

Day-10 3-D rms difference from the initial state: `u` 1.40 m/s (uniform)
and 0.90 m/s (variable); `theta` 1.31 K and 0.83 K; `w` 1.22×10⁻³ m/s and
1.11×10⁻³ m/s. Peak `|w|` at day 10 is 0.044 m/s (uniform) and 0.039 m/s
(variable). Nothing in the logs indicates instability or interface noise;
this is the reference the HTS run has to stay close to.

On the variable mesh at day 10, rms of `w` binned with the level rule in
DESIGN.md (`L = nint(log2(dx/dx_min))`) is 9.8×10⁻⁴ m/s on level 0
(16652 cells), 4.3×10⁻⁴ m/s on level 1 (3729 cells), and 1.5×10⁻³ m/s on
level 2 (9829 cells).

## What was not run

No real-data forecast. A pre-built `init.nc` is not published. The uniform
120 km `static.nc` and the CFSR intermediate file were downloaded and are
described in DESIGN.md; `init_atmosphere` was not run on them in this phase.
The 240 km uniform mesh `x1.10242` was downloaded and checked
(`nominalMinDc * R = 240000` m, `meshDensity = 1`) and was not integrated.
