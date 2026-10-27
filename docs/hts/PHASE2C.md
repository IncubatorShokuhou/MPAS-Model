# Phase 2c: a few percent of the cells are fine

Phase 2b's dry HTS was slower than all-fine on `x5.30210` (fine fraction 55%) and on `x4.10242` (fine fraction 49%). Both ceilings were under 1.6×, and full-domain diagnostics plus the extra halo calls used the tendency saving. This phase asks the same question on a mesh whose fine patch is a few percent of the cells. The mesh is [MESH_SMALL_PATCH.md](MESH_SMALL_PATCH.md): 24755 cells, level 0 is 885 cells (3.58%), realized spacing ratio 10.1×, four levels, large-step tendency ceiling 6.10× at `config_hts_max_level = 3`.

The binary is the phase 2b binary, single precision, gfortran, md5 `e8b4ea38e91c3100c3054d5f64acaf89`. No dynamics source changed for these runs. Four MPI ranks, OpenMP off, `--oversubscribe`, physics suite `none`, scalar advection off, `config_number_of_sub_steps = 6`, `config_len_disp = 16000`. The partition is `docs/hts/partitions/sp8.24755.graph.info.htsw.part.4` for every run, including all-fine. Relative L2 is `rms(test − reference) / rms(reference)` on the 6-hour history frame.

## dt = 60 s does not stay finite

The first 6-hour attempt used `config_dt = 60` (acoustic `dts = 10` s). All-fine finished with only underflow/denormal flags and wall 64.47 s. Every HTS clamp (max level 1, 2, and 3) exited 0 and signaled `IEEE_INVALID` and `IEEE_OVERFLOW` on the ranks. Those walls (64.16 s, 56.05 s, 55.63 s) are not a speed measurement. The atmosphere logs from that attempt were replaced by the runs below.

It is not the initial linear CFL. With `u ~ 35` m/s and `dt_level = 60 * 2^L`, the advective CFL on the shortest edge of each level stays under 0.21. The fine-edge acoustic CFL at `dts = 10` s is about 0.22 (the reference-frame sound speed gives a max of 0.65 at `dts = 30` s, and 10/30 of that is 0.22). Phase 2 already diverged when the interface took a 60 s acoustic step. The acoustic step was not lengthened here. The comparison uses `config_dt = 30` (`dts = 5` s), where all four runs finish with no NaN lines in the log and `u` stays near ±35 m/s through 05:59.

## Speed, 6 hours, dt = 30 s

720 fine steps. Rank 0 of the 4-level run owns 226 level-0 cells (work list 273, 845 edges), 116 level-1, 39 level-2, and 5521 level-3. `dyn_tend` parallel efficiency is 0.97–0.98.

| | all-fine | HTS, max level 1 | HTS, max level 2 | HTS, max level 3 |
| --- | ---: | ---: | ---: | ---: |
| `/usr/bin/time` wall | 132.48 s | 119.35 s | 113.80 s | 109.01 s |
| MPAS total | 132.124 s | 119.001 s | 113.474 s | 108.677 s |
| `time integration` | 127.102 s / 720 | 114.381 s / 720 | 108.829 s / 720 | 104.054 s / 720 |
| `atm_compute_dyn_tend` | 58.857 s / 2160 | 32.719 s / 3240 | 19.497 s / 3780 | 12.522 s / 4050 |
| `atm_advance_acoustic_step` | 27.856 s / 8640 | 28.747 s / 17280 | 29.501 s / 25920 | 29.768 s / 34560 |
| `atm_divergence_damping_3d` | 8.330 s / 8640 | 8.952 s / 17280 | 8.985 s / 25920 | 8.927 s / 34560 |
| `atm_compute_solve_diagnostics` | 10.536 s / 2160 | 13.730 s / 3240 | 15.655 s / 3780 | 16.252 s / 4050 |
| `atm_recover_large_step_variables` | 5.761 s / 2160 | 3.458 s / 3240 | 2.065 s / 3780 | 1.352 s / 4050 |
| `small_step_prep` | 2.732 s / 2160 | 4.033 s / 3240 | 4.613 s / 3780 | 4.686 s / 4050 |
| gap after every child timer | 8.29 s | 17.81 s | 23.37 s | 25.42 s |

Integration ratios, all-fine / HTS, higher is faster: **1.111×** at two levels, **1.168×** at three levels, **1.221×** at four levels. Wall ratios are 1.110, 1.164, and 1.215. The four-level integration timer is past the 1.2× bar. Two and three levels are faster than all-fine and short of that bar.

Where the four-level 23.0 s went. Tendency fell 58.86 → 12.52 s. The ceiling on that timer is `58.86 / 6.10 = 9.65` s, so the lists got most of the way (46.3 of the 49.2 s the ceiling allows). Acoustic stayed near 28–30 s while the call count went from 8640 to 34560. Diagnostics rose with the call count (2160 → 4050 is 1.875×; `10.54 × 1.875 = 19.8` s against a measured 16.25 s). Recovery fell 5.76 → 1.35 s. The untimed gap, which is the halo exchanges and the HTS driver, grew 8.3 → 25.4 s. That gap is why the net is 1.22× and not the 6× tendency ceiling. The same pattern, smaller, is the two- and three-level columns.

## Accuracy against all-fine

Same histories, last frame (6 hours). `w` reference rms is 1.056×10⁻³ m/s.

| | 2 levels | 3 levels | 4 levels |
| --- | ---: | ---: | ---: |
| `u` rel L2 / Linf | 2.867×10⁻⁴ / 0.340 m/s | 6.128×10⁻⁴ / 0.454 m/s | 7.890×10⁻⁴ / 0.452 m/s |
| `w` rel L2 / Linf | 0.214 / 1.33×10⁻² m/s | 0.377 / 2.28×10⁻² m/s | 0.438 / 2.15×10⁻² m/s |
| `w` rms(diff) / rms(ref) | 2.26×10⁻⁴ / 1.06×10⁻³ | 3.98×10⁻⁴ / 1.06×10⁻³ | 4.62×10⁻⁴ / 1.06×10⁻³ |
| `theta` rel L2 | 1.15×10⁻⁵ | 2.78×10⁻⁵ | 3.68×10⁻⁵ |
| `rho` rel L2 | 1.75×10⁻⁶ | 3.09×10⁻⁶ | 4.39×10⁻⁶ |
| dry-mass drift, HTS | −1.665×10⁻⁸ | +7.83×10⁻¹⁰ | −6.96×10⁻¹⁰ |
| dry-mass drift, all-fine | −7.95×10⁻¹¹ | −7.95×10⁻¹¹ | −7.95×10⁻¹¹ |
| mass frame difference | −1.66×10⁻⁸ | +8.63×10⁻¹⁰ | −6.16×10⁻¹⁰ |

`w` by region at four levels (the fast configuration): interface rel L2 0.516 (rms diff 1.72×10⁻³, ref rms 3.32×10⁻³), fine interior 0.410 (1.99×10⁻³), coarse interior 0.457 (1.96×10⁻⁴). The domain Linf, 2.15×10⁻² m/s, is on the interface. Deep coarse cells (`dist >= 2`) have rms diff 1.21×10⁻⁴. The interface error is larger than the `x5` 6-hour result. It is not a NaN spike.

`theta` at four levels: interface rms diff 0.108 K, coarse interior 0.013 K, domain Linf 1.09 K on the interface. `rho` interface rel L2 is 2.18×10⁻⁵ against 4.39×10⁻⁶ on the whole mesh.

At two levels the interface has 191 cells (fine interior 786, coarse interior 23778). `w` rms diff is 1.24×10⁻³ on the interface, 7.33×10⁻⁴ in the fine interior, and 1.52×10⁻⁴ in the coarse interior. The domain Linf is on the interface. The two-level mass drift, −1.7×10⁻⁸, is about 200× the all-fine drift. The four-level drift is −7×10⁻¹⁰, the same order as all-fine.

## One-level gate

The one-level uniform 120 km, 10-step gate was run on this binary in phase 2b: relative L2 of `u`, `w`, `theta`, and `rho` is 0, and the dry-mass frame difference is 0. Both drifts are −1.977×10⁻¹⁰. These phase 2c runs did not rebuild the binary, so that result still belongs to this executable. The source change in phase 3 is guarded by `hts_max_level_used > 0` and is checked again after that rebuild.

## Decision

Four-level `time integration` is 1.221× the all-fine timer on the same binary, the same partition, and a finite 6-hour history. That meets the bar for starting physics. Two- and three-level HTS are faster than all-fine and below 1.2×. The speedup is at `dt = 30` s. `dt = 60` s is not a usable HTS step on this mesh.
