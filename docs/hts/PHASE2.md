# Phase 2 results: dry hierarchical time stepping

Phase 2 stops here. The two-level scheme on `x5.30210` is stable for a day once the coarse level uses the fine acoustic step, mass is conserved at the owned interface edges, and the one-level uniform run is bit-identical to all-fine. It is also slower than all-fine. A third level was not run. Physics stays off.

All runs below are single precision, gfortran, 4 MPI ranks, OpenMP off, physics suite `none`, `config_split_dynamics_transport = false`, `config_number_of_sub_steps = 6`. The variable mesh is `x5.30210` (30210 cells, `dt = 180 s`). The uniform mesh is `x1.40962` (`dt = 450 s`). Relative L2 is `rms(test − reference) / rms(reference)`.

With `config_hts_max_level = 1` the 5× mesh is two levels: 16652 cells at level 0 (55.1%, about 51 km) and 13558 cells at level 1 (44.9%, the 91 km and 230 km cells clamped together). Interface cells, those with a neighbor at the other level, are 739. Fine interior cells are 16279. Coarse interior cells are 13192. The cell-count ceiling for two levels is `0.551 + 0.449/2 = 0.776` of an all-fine step, about 1.29×. The 1.44× figure in [DESIGN.md](DESIGN.md) is the three-level ceiling and does not apply to these runs.

## What the scheme does

`config_hts = false` calls the unmodified `atm_srk3`. `config_hts = true` calls `atm_hts_advance`.

- The coarse state is held at the start-of-interval snapshot while the fine level steps. Linear extrapolation of the previous coarse increment was tried and rejected: it fed the coarse acoustic transient back into the fine flux and `w` became NaN in about two hours.
- Interface mass flux is accumulated in double precision and imposed on the coarse acoustic step (`ru_p = target − ru_save` on interface edges). Owned-edge flux residual is 0. The halo residual printed in the log is a stale halo edge and is not a mass error.
- On level L the acoustic substep count is `config_number_of_sub_steps * 2^L`, so the acoustic `dts` matches the fine acoustic `dts` (30 s here). The slow tendency loop (`atm_compute_dyn_tend`) and `atm_compute_solve_diagnostics` are not masked. Acoustic, damping, recovery, setup, and finish skip inactive cells with a `cycle`, still walking every cell or edge. Halo exchanges stay on the full-domain groups.

## Gate 1: one level, 120 km

`config_hts = true` and a uniform mesh (every cell level 0). Relative L2 of `u`, `w`, `theta`, and `rho` is 0 against the HTS-off history of the same binary at 1 hour, 1 day, and 10 days, and against the phase-1 10-day history. Dry-mass drift at day 10 is −2.115×10⁻¹⁰, the same as the baseline. The 10-day HTS-on wall time is 546 s against the phase-1 535 s. That case has no coarse level, so it is a regression, not a speed test.

## Gate 2: two levels, speed

Wall clock, same binary family, HTS versus all-fine of that binary. The acoustic step is matched (30 s on both levels).

| run | all-fine wall | HTS wall | HTS / all-fine |
| --- | ---: | ---: | ---: |
| 6 h (120 fine steps) | 24.13 s | 39.80 s | 1.65× |
| 1 day (480 fine steps) | 94.89 s | 154.09 s | 1.62× |

The integration timer for the 6 h pair is 22.30 s all-fine and 37.84 s HTS (1.70×). A net speedup is not reachable with the current loop and halo structure on this mesh.

### Where the 6 hours go

Child timers from the uninstrumented 6 h pair. Halo exchanges are compiled out of the timer unless OpenACC is on, so they sit in the gap between `time integration` and the sum of its children.

| timer | all-fine | calls | HTS | calls | per-call ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| `atm_compute_dyn_tend` | 10.45 s | 360 | 15.84 s | 540 | 1.01 |
| `atm_advance_acoustic_step` | 4.98 s | 1440 | 5.49 s | 2880 | 0.55 |
| `atm_divergence_damping_3d` | 1.43 s | 1440 | 1.53 s | 2880 | 0.53 |
| `atm_compute_solve_diagnostics` | 1.98 s | 360 | 2.94 s | 540 | 0.99 |
| `atm_recover_large_step_variables` | 1.12 s | 360 | 1.14 s | 540 | 0.68 |
| remainder of `time integration` | 2.34 s | — | 10.90 s | — | — |
| `time integration` | 22.30 s | 120 | 37.84 s | 120 |  |

HTS performs 180 RK3 calls in those 120 clock steps: 120 fine and 60 coarse. Each coarse call is a full three-stage step.

A second 6 h pair compiled the existing halo-exchange timer in, and timed the HTS driver copies and the wind reconstruction. Integration was 23.51 s all-fine and 38.08 s HTS. The split of that pair:

| piece | all-fine | HTS | HTS − all-fine |
| --- | ---: | ---: | ---: |
| dynamics tendencies | 10.71 s | 15.86 s | +5.15 s |
| halo exchanges | 1.59 s (4560 calls) | 8.16 s (8280 calls) | +6.56 s |
| acoustic step | 5.16 s | 5.54 s | +0.38 s |
| divergence damping | 1.46 s | 1.52 s | +0.06 s |
| solve diagnostics | 2.14 s | 2.94 s | +0.80 s |
| HTS snapshot / flux / mass patch | 0 | 1.55 s | +1.55 s |
| wind reconstruction | 0.14 s | 0.55 s | +0.40 s |

Three facts follow.

1. **Tendencies dominate and do not get cheaper.** `atm_compute_dyn_tend` is 44–47% of the all-fine step and costs the same per call with HTS on (0.029 s). The coarse tick adds another full-domain evaluation, 180 extra calls, about +5.4 s. Masking that loop with an index list, so a coarse call touched only the 44.9% coarse cells, would bring the 6 h tendency time to about 8.1 s, a 2.3 s saving against all-fine. That saving is smaller than the current penalty.

2. **The acoustic kernel is not the place a speedup was lost, and it is not a source of one either.** Matching `dts` doubles the acoustic call count (1440 to 2880). The `cycle` mask cuts the per-call time to 0.55×, close to the average active fraction of 0.50, so the mask overhead inside the acoustic loop is small. Net acoustic time is about 10% higher than all-fine. Damping behaves the same way. Leaving the coarse acoustic step at 60 s would be the only way for this kernel to cost less, and that configuration is the unstable one below.

3. **Global halos on an unweighted partition are a large fraction of the slowdown.** Call count rises by 1.82× (4560 to 8280), which at the all-fine cost per exchange would be about 2.9 s of halo time. The measured halo time is 8.16 s. The extra ~5 s is wait inside the exchange. The METIS partition is unweighted: rank 2 owns 7552 cells and 0 of them are level 0, while rank 0 is 91% level 0. Every rank still enters every acoustic halo. The fast rank pays the wait, and the timer records it as halo time.

Bookkeeping (snapshot, restore, flux sum, mass patch) is 1.55 s that all-fine does not pay. Wind reconstruction is called on the way into each level and again at the end of the tick.

### Is a net speedup reachable on this design?

On this design, on this mesh, no. The measured result is a 1.65× slowdown at 6 hours and 1.62× at 1 day. The pieces that would have to change before a speedup is even arithmetically available:

- Stop visiting inactive cells in the tendency and diagnostic loops. Necessary, and on this mesh worth about 2 s per 6 hours once the extra coarse calls are included. It does not cover the deficit by itself.
- Stop doing a full-domain halo on every acoustic substep of a masked step, and weight the partition by level. The halo line is the largest single addition in the profile (+6.6 s), and most of that is imbalance rather than message count.
- Take a longer acoustic step on long edges. Matched `dts` keeps the acoustic work flat. Edges with `dcEdge ≥ 150 km` are the ones that have room; the interface does not. See the next section. This was not implemented.
- Shrink the driver copies. 1.55 s of bookkeeping is already a 7% tax on the all-fine integration.

Even if all four landed, the two-level cell-count ceiling on `x5.30210` is 1.29×, because 55% of the cells are fine and must take every step. This mesh is a correctness case. A mesh whose fine region is a small patch is the speed case, and it still needs the halo and tendency changes or the same tax shows up there.

A third level under the current rule would multiply the level-2 acoustic count by 4 and add another full-domain tendency evaluation. That spends more of the time the coarse region was supposed to save. It was not run.

## Acoustic step and the interface

The acoustic step is horizontally explicit and vertically implicit, so the CFL that matters is `cs * dts / dcEdge`, with `cs = sqrt(1.4 * p / ρ)` from the all-fine history. The fine namelist runs at `dts = 30 s`.

| edges | count | median `dcEdge` | CFL at 30 s, max | CFL at 60 s, median / max |
| --- | ---: | ---: | ---: | --- |
| both cells level 0 | 49580 | 50 km | 0.270 | 0.42 / 0.54 |
| interface | 739 | 68 km | 0.186 | 0.31 / 0.37 |
| both cells level 1, `dc < 80 km` | 3500 | 74 km | 0.182 | 0.29 / 0.37 |
| both cells level 1, `dc ≥ 150 km` | 28586 | 233 km | 0.070 | 0.09 / 0.14 |

The fine mesh operates at a maximum horizontal acoustic CFL of 0.27. At a doubled acoustic step (`dts = 60 s`, which is what level 1 does if the substep count is left at 6):

- Edges longer than 150 km, 71% of the coarse edges, stay at or below 0.14. At `dts = 120 s` their maximum is 0.28, the fine operating point. Those columns have CFL room.
- Interface edges and the short coarse edges in the transition (`dc < 80 km`) reach 0.37, above the fine operating point.

So the horizontal-CFL hypothesis holds for the long coarse edges. It does not hold for the interface or for the transition edges that the level-1 bin also contains. The imposed interface edge does not integrate a pressure gradient; the cells beside it do, at whatever `dts` the coarse step is using.

The check that was actually integrated is the 2-hour run with the held snapshot and the nominal substep count (`dts = 60 s` on level 1), against all-fine. `w` at 2 hours:

| region | rel L2 | rms(difference) | Linf |
| --- | ---: | ---: | ---: |
| all cells | 1.84 | 2.90×10⁻⁴ m/s | 3.44×10⁻² m/s |
| interface | 7.82 | 1.56×10⁻³ m/s | 3.44×10⁻² m/s |
| coarse, one cell from the interface | 3.26 | 7.15×10⁻⁴ m/s | 1.23×10⁻² m/s |
| coarse, two or more cells in | 1.08 | 1.96×10⁻⁴ m/s | 4.35×10⁻³ m/s |
| fine, two or more cells in | 0.33 | 4.32×10⁻⁵ m/s | 1.04×10⁻³ m/s |

The peak sits on an interface cell. The adjacent coarse ring is next. The deep coarse interior is quieter than the interface by about 8× in rms and 8× in the peak, and it is already larger than its own reference `w` (reference rms in that region is 1.8×10⁻⁴ m/s). A sound wave crosses the coarse region in well under two hours, so that deep-interior signal can be noise that traveled, rather than a local CFL failure. This run was not continued as a split-step experiment (fine `dts` on the interface and the short transition edges, longer `dts` only on the long edges). An earlier continuation of the same 60 s configuration, since overwritten by the matched-dts 1-day run, reached NaN before 6 hours. The 2-hour file is the one that remains.

The hypothesis is supported as a statement about where the error is and about which edges have CFL room. It is not a demonstration that the coarse interior would stay stable for a day at `dts = 60 s`, and the profile says that even a successful split would not by itself make this mesh faster than all-fine.

Matching `dts` removes the spike. The 6-hour and 1-day numbers in the next section are that configuration. Owned flux residual stays 0.

## `w` relative L2 of 0.18

Reference `w` at 6 hours has rms 2.06×10⁻⁴ m/s. The difference has rms 3.72×10⁻⁵ m/s. The ratio is 0.180. The same comparison at 1 day is 5.66×10⁻⁵ / 2.69×10⁻⁴ = 0.210. Both are the JW vertical velocity before the baroclinic wave has much vertical motion (day-10 surface pressure anomalies in [BASELINE.md](BASELINE.md) are the part of the test where `w` is large; that pair was not run). Absolute differences are a few ×10⁻⁵ m/s, and the largest pointwise `|Δw|` at 6 hours is 9.9×10⁻⁴ m/s.

Region breakdown, matched acoustic step, against all-fine:

|  | 6 h `w` rel L2 | 6 h rms(diff) | 1 day `w` rel L2 | 1 day rms(diff) |
| --- | ---: | ---: | ---: | ---: |
| all | 0.180 | 3.72×10⁻⁵ | 0.210 | 5.66×10⁻⁵ |
| interface (739 cells) | 0.299 | 7.95×10⁻⁵ | 0.207 | 5.36×10⁻⁵ |
| fine interior | 0.177 | 4.02×10⁻⁵ | 0.221 | 7.36×10⁻⁵ |
| coarse interior | 0.167 | 2.86×10⁻⁵ | 0.140 | 2.22×10⁻⁵ |
| fine, distance ≥ 2 | 0.175 | 3.96×10⁻⁵ | 0.222 | 7.42×10⁻⁵ |
| coarse, distance ≥ 2 | 0.162 | 2.73×10⁻⁵ | 0.139 | 2.15×10⁻⁵ |

At 6 hours the interface difference is about twice the interior difference. It is a mild elevation spread over the 739 interface cells, with an interface Linf of 7.6×10⁻⁴ m/s, below the domain Linf of 9.9×10⁻⁴ m/s (that point is in the coarse interior). At 1 day the fine interior holds the largest `w` difference, and the coarse interior is the smallest. The 0.18 is the small denominator. It is not a localized interface spike of the kind the 60 s acoustic step produced.

The other fields at the same times:

| field | 6 h rel L2 | 6 h Linf | 1 day rel L2 | 1 day Linf |
| --- | ---: | ---: | ---: | ---: |
| `u` | 2.74×10⁻⁴ | 0.108 m/s | 4.60×10⁻⁴ | 0.089 m/s |
| `theta` | 3.61×10⁻⁶ | 0.054 K | 1.14×10⁻⁵ | 0.091 K |
| `rho` | 2.41×10⁻⁶ | 3.7×10⁻⁵ kg/m³ | 4.05×10⁻⁶ | 1.1×10⁻⁴ kg/m³ |

`u` rms of the reference is 9.0 m/s, so the 6-hour difference rms is 2.5×10⁻³ m/s. Interface-edge `u` rel L2 at 6 hours is 1.75×10⁻³. `theta` and `rho` stay near or under the 10⁻⁵ uniform target at 6 hours; `theta` at 1 day is 1.1×10⁻⁵. Dry-mass drift of the HTS run is −3.2×10⁻¹⁰ at 6 hours and +1.3×10⁻⁹ at 1 day. The all-fine drifts are −2.4×10⁻¹⁰ and −2.0×10⁻¹⁰. The 1-day mass difference between the two runs is +1.5×10⁻⁹. That is larger than the baseline drift and still a relative 10⁻⁹, with owned flux residual 0.

## What was not done

- Three levels (`config_hts_max_level = 2`). The two-level profile already shows a slowdown, and level 2 would take four times the fine acoustic count.
- A 10-day two-level comparison. The longest stable multi-level run is 1 day.
- A split acoustic step (long `dts` only where `dcEdge` allows it).
- Scalar transport and `config_split_dynamics_transport = true`.
- Physics, weighted partitions, per-level halo groups.

`docs/hts/scripts/compare_interface.py` prints the region table above (`--max-level 1`).
