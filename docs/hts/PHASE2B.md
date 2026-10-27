# Phase 2b: index lists and a level-balanced partition

Phase 2b keeps the phase 2 scheme (held coarse state, matched acoustic step, imposed interface flux) and changes which cells a call touches. The goal was a dry HTS that stays as close to all-fine as phase 2 and is faster than all-fine on the same ranks. On the meshes this VM can run, it is not faster.

All new timings are one binary, single precision, gfortran, 4 MPI ranks, OpenMP off, `--oversubscribe`, physics suite `none`, scalar advection off, `config_number_of_sub_steps = 6`. Relative L2 is `rms(test − reference) / rms(reference)`. The "before" column is the phase 2 pair in [PHASE2.md](PHASE2.md). The "after" column is this binary.

## What changed

Tendency, dissipation, acoustic, and damping loops use the index lists in `hts_loops.h`. A work-list call iterates the active level plus two cell rings. An acoustic or damping call iterates the edges and cells of the active level, including the interface edges that level owns. Recovery of the large-step edge state still walks every local edge, because the halo exchange of `u` needs the pure-halo edges of the active level. When HTS is off, or the mesh has a single level, the macros walk the original `start..end` span in the original order, and the lists are not built.

Diagnostics stay full-domain. `solve_diagnostics_work` sets `hts_lists_work = .false.` for that subroutine. A masked diagnostic leaves a fringe. The next level's first RK tendency reads those diagnostics before that level recomputes them (diagnostics run at the end of each RK stage). Widening the diagnostic ring from 2 to 8 moved the `|Δu|` peak outward with the ring. Writing every local cell was bit-identical to the previous HTS binary, so the diagnostic loop is the straight full-domain loop and does not pay indirect addressing.

`delsq_*` and `dpdz` are zeroed at the start of a masked step. Hyperdiffusion reads `delsq` on neighbors of the active level. Without the zero, the fringe produced IEEE invalid and overflow. With the zero, the 6-hour history stays bit-identical to the previous HTS binary.

Halo groups are unchanged. See below.

The partition is level-balanced. `docs/hts/scripts/hts_partition.py` builds a METIS multi-constraint graph and runs `gpmetis`. Usage and the weight convention are in [DESIGN.md](DESIGN.md) section 7. The 4-rank files used below are in `docs/hts/partitions/`.

A longer acoustic step on long edges was not added. The CFL tables below are why.

## Correctness gates

One level, uniform 120 km (`x1.40962`, `dt = 450 s`), 10 steps (01:15:00), this binary, `config_hts_max_level = 1`. The log line is `HTS enabled with a single level; masked integrator, no subcycling`. Against HTS-off of the same binary, relative L2 of `u`, `w`, `theta`, and `rho` is 0, and the dry-mass frame difference is 0. Both drifts are −1.977×10⁻¹⁰.

Six hours on `x5.30210`, original unweighted partition, this binary versus the phase 2 binary: all-fine against all-fine and HTS against HTS are both bit-identical (relative L2 0 on `u`, `w`, `theta`, `rho`, mass frame difference 0). Index lists reproduce the committed scheme. The all-fine path, which is the sequential branch of the same macros, did not move.

## Profile, x5.30210, 6 hours

`dt = 180 s`, 120 steps, acoustic `dts = 30 s` on every level. Two levels: 16652 cells at level 0 (55.1%) and 13558 at level 1 (44.9%). Large-step work over one coarse interval is `0.551 + 0.449/2 = 0.776` of an all-fine step, a ceiling of 1.29×. That ceiling is for the slow tendency only. With a matched acoustic step, every cell takes the same number of acoustic updates per unit time, so the acoustic work does not shrink.

The phase 2 pair used cycle masks and the unweighted partition. This binary's original-partition pair is index lists without rebalancing. The balanced pair uses `x5.30210.graph.info.htsw.part.4`.

| | phase 2 all-fine | phase 2 HTS | index lists, original part, all-fine | index lists, original part, HTS | index lists, balanced, all-fine | index lists, balanced, HTS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `/usr/bin/time` wall | 24.13 s | 39.80 s | 26.05 s | 41.11 s | 26.95 s | 29.28 s |
| MPAS total | — | — | 25.72 s | 40.76 s | 26.61 s | 28.96 s |
| `time integration` | 22.30 s | 37.84 s | 24.17 s | 39.16 s | 25.06 s | 27.34 s |
| `atm_compute_dyn_tend` | 10.45 s / 360 | 15.84 s / 540 | 10.92 s / 360 | 11.20 s / 540 | 11.43 s / 360 | 9.88 s / 540 |
| `atm_advance_acoustic_step` | 4.98 s / 1440 | 5.49 s / 2880 | 5.24 s / 1440 | 5.65 s / 2880 | 5.44 s / 1440 | 5.80 s / 2880 |
| `atm_divergence_damping_3d` | 1.43 s / 1440 | 1.53 s / 2880 | 1.59 s / 1440 | 1.73 s / 2880 | 1.64 s / 1440 | 1.79 s / 2880 |
| `atm_compute_solve_diagnostics` | 1.98 s / 360 | 2.94 s / 540 | 2.14 s / 360 | 3.04 s / 540 | 2.31 s / 360 | 3.22 s / 540 |
| `atm_recover_large_step_variables` | 1.12 s / 360 | 1.14 s / 540 | 1.21 s / 360 | 1.22 s / 540 | 1.35 s / 360 | 1.12 s / 540 |
| gap after those five timers | 2.34 s | 10.90 s | 3.07 s | 16.32 s | 2.89 s | 5.53 s |

`dyn_tend` parallel efficiency on the HTS run is 0.83 with the original partition and 0.97 with the balanced partition. Rank 0 of the original partition owns 6869 of its 7559 cells at level 0. The balanced file gives each rank both levels: rank 0 has 7014 cells (3927 fine, 3087 coarse), rank 1 has 7779 (4288, 3491), rank 2 has 7647 (4158, 3489), rank 3 has 7770 (4279, 3491). Wall ratio, balanced, same binary, same 4 ranks: `29.28 / 26.95 = 1.086`. Integration ratio `27.34 / 25.06 = 1.091`. Index lists without the balanced partition stay at `41.11 / 26.05 = 1.58`, the same slowdown as phase 2, because the critical rank is still almost all fine cells and still waits in the collective.

The last row is `time integration` minus the five timers above it, which is how [PHASE2.md](PHASE2.md) defined the gap. It still contains the smaller child timers (setup, moist, physics, vertical coefficients, `small_step_prep`, finish) plus the halo exchanges and the HTS driver copies. Halos are compiled out of the timer unless OpenACC is on. Subtracting every child timer instead, the balanced pair's gap is 1.34 s all-fine and 3.43 s HTS.

Where the balanced 2.29 s of integration went. Tendency fell by 1.55 s (11.43 → 9.88). The large-step ceiling on that timer is `11.43 × 0.776 = 8.86 s`, so the lists recovered 1.55 of the 2.57 s that ceiling allows. Diagnostics grew 0.91 s (full domain, 360 → 540 calls). The all-child gap grew 2.08 s (1.34 → 3.43). Acoustic grew 0.37 s, damping 0.15 s, `small_step_prep` 0.29 s, and the moist, physics, and vertical-coefficient timers together about 0.28 s. Recovery fell 0.24 s.

## Halo exchanges

Phase 2 timed them once: 1.59 s and 4560 calls all-fine, 8.16 s and 8280 calls with HTS. About 5 s of that gap was wait. The unweighted partition had given one rank 7552 cells and zero fine cells, so that rank entered every fine exchange and waited. Balancing removes that wait. It does not remove the calls. HTS still performs a fine RK3 and a coarse RK3, and every acoustic substep still exchanges.

Skipping an exchange on a rank that has nothing to send deadlocks. The exchange is a collective, and every rank must call it. With the balanced partition every rank owns every level, so a per-level group would still include every rank. Packing only the edges that changed would cut message size. On 4 ranks and 30k cells the exchange is latency-bound (phase 2, about 0.35 ms per call). The extra cost is the extra call, and that call stays until the acoustic substep count drops. No sparse exchange was added.

## Speed-case mesh

A variable mesh with a few percent of its cells in a 16× patch, small enough for this VM, is not in the public set, and this machine has neither MPAS-Tools nor JIGSAW (`import mpas_tools` fails, `jigsaw` is not on `PATH`). Cell count scales as `1/dx²`. A patch that is 3% of the cells at a 16× ratio covers `0.03 / (0.03 + 0.97 × 256) = 0.012%` of the sphere. The public NCAR meshes with a large ratio are the other way around: `x4.163842` is 92–25 km and 163842 cells, and the 60–3 km mesh is 835586 cells. Neither fits this VM. `meshDensity` was not edited by hand.

The mesh that was integrated is the experimental 4× mesh `x4.10242` (10242 cells):

https://www2.mmm.ucar.edu/projects/mpas/atmosphere_meshes/experimental/x4.10242.tar.gz

Last-Modified Fri, 26 Mar 2021, Content-Length 6858237. `meshDensity` runs from 1/256 to 0.999, ratio `(md_max/md_min)^0.25 = 4.00`. `dcEdge` is 76–398 km. Unclamped levels are 4976, 1515, and 3751 cells. The fine fraction is 48.6% at every `config_hts_max_level`.

| clamp | cells by level | large-step work / all-fine | ceiling |
| --- | --- | ---: | ---: |
| max level 1 | 4976 + 5266 | 0.743 | 1.35× |
| max level 2 | 4976 + 1515 + 3751 | 0.651 | 1.54× |

The same ratio on `x4.2562` (2562 cells, same URL directory, tarball 1.6 MB) is 1236 / 385 / 941. The fine fraction is again about 48%, and 2562 cells is latency-dominated. It was not integrated.

`x4.10242` was initialized as JW case 2, 26 levels, and integrated 6 hours at `dt = 360 s` (`config_len_disp = 100000`, acoustic `dts = 60 s`). All three runs used the 3-level balanced partition `x4.10242.graph.info.htsw.L2.part.4` (rank cells 2636, 2633, 2359, 2614; level counts on rank 0 are 1281, 390, 965). A 2-level partition is also in `docs/hts/partitions/` and was not the file these runs read.

| | all-fine | HTS, max level 1 | HTS, max level 2 |
| --- | ---: | ---: | ---: |
| `/usr/bin/time` wall | 4.64 s | 5.23 s | 5.28 s |
| MPAS total | 4.33 s | 4.90 s | 4.96 s |
| `time integration` | 3.62 s / 60 | 4.18 s / 60 | 4.25 s / 60 |
| `atm_compute_dyn_tend` | 1.71 s / 180 | 1.50 s / 270 | 1.36 s / 315 |
| `atm_advance_acoustic_step` | 0.79 s / 720 | 0.89 s / 1440 | 0.92 s / 2160 |
| `atm_divergence_damping_3d` | 0.26 s / 720 | 0.31 s / 1440 | 0.31 s / 2160 |
| `atm_compute_solve_diagnostics` | 0.27 s / 180 | 0.43 s / 270 | 0.48 s / 315 |
| gap after every child timer | 0.24 s | 0.65 s | 0.76 s |

Integration ratios: `4.18 / 3.62 = 1.15` at two levels and `4.25 / 3.62 = 1.18` at three levels. Wall ratios are 1.13 and 1.14, and on a 5-second run the startup share is large, so the integration timer is the comparison. Tendency moved toward the ceiling (1.71 → 1.50 against a 2-level ceiling of 1.27 s, and → 1.36 against a 3-level ceiling of 1.12 s) and the diagnostic and remainder lines grew by more than that saving. Diagnostics track the call count: `0.27 × 270/180 = 0.41` versus a measured 0.43, and `0.27 × 315/180 = 0.48` versus a measured 0.48. All three runs finished, and `u`, `w`, `theta`, and `rho` are finite.

## Accuracy against all-fine

`x5.30210`, 6 hours, balanced partition, this binary. The original-partition pair is bit-identical to phase 2, and the balanced pair matches these digits.

| field | rel L2 | Linf | region notes for `w` |
| --- | ---: | ---: | --- |
| `u` | 2.740×10⁻⁴ | 0.108 m/s | interface-edge rel L2 1.75×10⁻³ |
| `w` | 0.180 | 9.92×10⁻⁴ m/s | rms(diff) 3.72×10⁻⁵, rms(ref) 2.06×10⁻⁴ |
| `theta` | 3.609×10⁻⁶ | 0.054 K | |
| `rho` | 2.411×10⁻⁶ | 3.71×10⁻⁵ kg/m³ | |

`w` by region: interface 0.299 (rms diff 7.95×10⁻⁵), fine interior 0.177 (4.02×10⁻⁵), coarse interior 0.167 (2.86×10⁻⁵). The domain Linf is in the coarse interior. Dry-mass drift is −3.230×10⁻¹⁰ (HTS) and −2.387×10⁻¹⁰ (all-fine). Frame difference −8.43×10⁻¹¹.

`x4.10242`, 6 hours, same partition for all three runs.

| | 2 levels vs all-fine | 3 levels vs all-fine |
| --- | ---: | ---: |
| `u` rel L2 / Linf | 4.49×10⁻⁴ / 0.157 m/s | 1.09×10⁻³ / 0.158 m/s |
| `w` rel L2 / Linf | 0.368 / 1.24×10⁻³ m/s | 0.444 / 1.57×10⁻³ m/s |
| `w` rms(diff) / rms(ref) | 7.08×10⁻⁵ / 1.92×10⁻⁴ | 8.55×10⁻⁵ / 1.92×10⁻⁴ |
| `theta` rel L2 | 5.99×10⁻⁶ | 1.02×10⁻⁵ |
| `rho` rel L2 | 5.30×10⁻⁶ | 1.61×10⁻⁵ |
| dry-mass drift, HTS | −1.34×10⁻⁹ | −8.04×10⁻¹⁰ |
| dry-mass drift, all-fine | −6.98×10⁻¹¹ | −6.98×10⁻¹¹ |
| mass frame difference | −1.27×10⁻⁹ | −7.34×10⁻¹⁰ |

`w` rms difference at two levels: interface 8.58×10⁻⁵, fine interior 6.51×10⁻⁵, coarse interior 7.47×10⁻⁵. The Linf is in the fine interior (interface Linf 6.2×10⁻⁴). At three levels the interface rms difference is 1.13×10⁻⁴, the fine interior 7.69×10⁻⁵, and the coarse interior 8.95×10⁻⁵. The Linf (1.57×10⁻³) is in the coarse interior, two or more cells in from the interface. The elevation at the interface is mild. It is the same shape as the `x5` result, on a smaller `w`.

The owned-edge flux residual printed on the first coarse steps is 1×10⁻¹⁰ to 2×10⁻¹⁰. The 6-hour mass drift on `x4` is about 1×10⁻⁹, larger than the all-fine drift and still a relative 10⁻⁹. On `x5` the drift stays at 3×10⁻¹⁰.

## Longer acoustic steps

Not implemented.

On `x5`, phase 2 already measured the horizontal acoustic CFL at the fine operating point of 0.27 (`dts = 30 s`). At 60 s the interface max is 0.37 and the short coarse edges (`dc < 80 km`) max is 0.37. Edges longer than 150 km max out at 0.14. The 2-hour run at a uniform 60 s coarse step put the `w` peak on the interface, and a longer continuation of that configuration became NaN before 6 hours.

On `x4` the fine acoustic step is already 60 s. From the all-fine history (`cs = sqrt(1.4 p / ρ)`):

| edges | count | median `dcEdge` | CFL at 60 s, max | CFL at 120 s, max |
| --- | ---: | ---: | ---: | ---: |
| fine | 14745 | 99.5 km | 0.277 | — |
| interface | 355 (2 levels) | 131 km | 0.175 | 0.35 |
| coarse | 15620 | 356 km | 0.167 | 0.33 |

The coarse median would be comfortable at 120 s (CFL near 0.12). The interface and the shortest coarse edges would sit above the fine operating point. That is the same split the `x5` run rejected. Acoustic time is also not the line that loses the comparison: on the balanced `x5` pair it is 5.80 s against 5.44 s, while diagnostics and the untimed remainder grow by 0.91 s and 2.08 s.

## Is the approach viable?

On these meshes, no. The measured speedups are ratios below 1.

- `x5.30210`, 4 ranks, 6 hours, index lists plus a level-balanced partition: wall 29.28 s versus 26.95 s (1.09×), integration 27.34 s versus 25.06 s (1.09×).
- `x4.10242`, 4 ranks, 6 hours, same partition: integration 4.18 s versus 3.62 s at two levels (1.15×) and 4.25 s versus 3.62 s at three levels (1.18×).

The arithmetic that would have to close, and does not:

1. The large-step ceiling is 1.29× on `x5` and 1.35× or 1.54× on `x4`, because 55% and 49% of the cells are fine. Index lists recover about 1.55 s of a 2.57 s tendency ceiling on `x5`, and about 0.22 s of a 0.44 s ceiling on `x4` at two levels. They recover none of the acoustic work. Matched `dts` keeps the acoustic cell-updates equal to all-fine.
2. Diagnostics have to stay full-domain to keep the next level's tendency exact. They scale with the extra RK calls (1.5× at two levels, 1.75× at three).
3. Halo call count stays near 1.8×. Balancing removed the idle-rank wait (`dyn_tend` parallel efficiency 0.83 → 0.97) and left the calls. A per-level exchange that any rank skips will deadlock, and after balancing no rank can skip.
4. A public mesh with a few-percent fine patch at 16× is not available at a size this VM can integrate, and it cannot be generated here. That mesh would still pay items 2 and 3.
5. The remaining lever is fewer acoustic substeps on long edges. The interface CFL and the earlier NaN say that step is not stable if the interface takes it. It was not added.

The scheme matches the previous HTS solution bit for bit, the one-level uniform run is bit-identical, and the multi-level mass drift stays near 10⁻¹⁰ on `x5` and near 10⁻⁹ on `x4`. The `w` difference is a small absolute error (a few ×10⁻⁵ m/s) with a mild interface elevation, spread over the domain. Those gates hold. They do not make the step faster than taking every cell at the fine `dt`.

Physics was not started.
