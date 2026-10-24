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
| `sp8.24755`, 6 h, 4 levels, dry, `dt = 30` s | 7.89×10⁻⁴ | 0.438 | 3.68×10⁻⁵ | 4.39×10⁻⁶ | −6.96×10⁻¹⁰ | −7.95×10⁻¹¹ |
| `sp8.24755`, 1 h, 4 levels, Kessler (phase 3 restore, phase 3b, and the scheme-mask rebuild, same frame) | 6.16×10⁻⁴ | 0.293 | 8.41×10⁻⁶ | 1.28×10⁻⁵ | −2.86×10⁻¹¹ | −9.60×10⁻¹¹ |
| `sp8.24755`, 30 min, 4 levels, WSM6 | 3.63×10⁻⁴ | 0.486 | 3.67×10⁻⁶ | 1.03×10⁻⁵ | −7.88×10⁻⁹ | −1.31×10⁻¹⁰ |
| `x1.40962` real, 6 h, 1 level, `mesoscale_reference` | 0 | 0 | 0 | 0 | −4.31×10⁻¹¹ | −4.31×10⁻¹¹ |

`w` relative L2 is large because the reference `w` is small (about 10⁻⁴ to 10⁻³ m/s on the dry JW runs). The absolute rms difference on `sp8.24755` at 6 hours is 4.6×10⁻⁴ m/s. The largest piece sits on the interface. It is not a NaN. `u`, `theta`, and `rho` stay at relative 10⁻⁴, 10⁻⁵, and 10⁻⁶. Dry-mass drift stays at or below about 10⁻⁸; the fast four-level dry run and the Kessler run are 10⁻¹⁰ and 10⁻¹¹.

`dt = 60` s on `sp8.24755` is not in the table. All-fine finished. Every HTS clamp signaled invalid and overflow. The acoustic step was not lengthened. The dry speedup is at `dt = 30` s.

## Speed against all-fine

Ratio is all-fine / HTS on `time integration`. Higher than 1 is faster.

| case | all-fine | HTS | ratio |
| --- | ---: | ---: | ---: |
| `x5.30210`, 6 h, 2 levels, 55% of cells fine | 25.06 s | 27.34 s | 0.92× (slower) |
| `x4.10242`, 6 h, 3 levels, 49% of cells fine | 3.62 s | 4.25 s | 0.85× (slower) |
| `sp8.24755`, 6 h, 4 levels, 3.6% of cells fine, dry | 127.10 s | 104.05 s | **1.221×** |
| `sp8.24755`, 1 h, 4 levels, Kessler, phase 3 restore | 23.72 s | 23.04 s | 1.029× |
| `sp8.24755`, 1 h, 4 levels, Kessler, phase 3b column mask | 23.35 s | 19.02 s | **1.227×** |
| same Kessler hour, rebuilt with WSM6/Thompson segments (bit-identical history) | 22.91 s | 18.82 s | 1.217× |
| `sp8.24755`, 30 min, 4 levels, WSM6 | 15.91 s | 10.65 s | **1.494×** |
| `x1.40962` real, 6 h, 1 level, full suite | 226.16 s | 218.99 s | bit-identical; the 3% is timer scatter |

The dry win is real and narrow. Level 0 is 885 of 24755 cells, the spacing ratio is 10×, and the large-step tendency ceiling is 6.10×. Tendency time fell from 58.9 s to 12.5 s. Acoustic time stayed near 29 s, diagnostics rose from 10.5 s to 16.3 s, and the untied halo-and-driver gap rose from 8.3 s to 25.4 s. The net is 1.22×, not 6×.

Phase 3's restore coupling did not keep that net. Microphysics ran once per level step and still visited every column, so 225 calls against 120 put the Kessler integration timer at 1.029×. Phase 3b skips the inactive columns. On the same 1-hour Kessler case, `time integration` is 23.349 s all-fine and 19.024 s with four levels, **1.227×**. Wall is 25.18 s against 20.84 s. `mp_kessler` itself is 5.49× (0.90828 s / 0.16530 s), close to the 6.10× cell-count ceiling. The wrapper around it is 4.10× (2.171 s / 0.529 s), because the copy loops still scan the index. Diagnostics, `physics_get_tend`, and the moist-coefficient loops still walk every column once per level step, and the halo gap grew from 1.76 s to 4.82 s. That is why the net is 1.23×.

On the 120 km mesoscale step the physics driver is 59% of integration and microphysics is 16%. That mesh is one level, so the mask stays off. WSM6, Thompson, and Thompson-aerosols now run on contiguous active columns instead of aborting. On `sp8.24755` for 30 minutes, WSM6 integration is 15.914 s all-fine and 10.654 s with four levels, **1.494×**. `mp_wsm6` itself is 5.43× (3.592 s / 0.662 s). Dry-mass frame difference is −7.75×10⁻⁹ and `u` relative L2 is 3.63×10⁻⁴. Thompson was not timed: its lookup tables are not on this machine. Radiation and PBL on a variable-resolution mesh were not timed.

## Bottom line

Dry HTS on a small fine patch is faster than all-fine, by 1.22× on the integration timer, and it stays within a few parts in 10⁴ on `u` with dry-mass drift near 10⁻¹⁰. Kessler on the same mesh is 1.23× on `time integration` (1.029× before the columns were skipped). WSM6 on a 30-minute run of the same mesh is 1.49×. The one-level path matches the all-fine history bit for bit, including a 10-step dry gate and a 2-step CFSR `mesoscale_reference` smoke, and that smoke still matches after WSM6 and Thompson were put on the column mask. The 6-hour CFSR forecast that matches bit for bit was run on the phase 3 binary. The acoustic step was not lengthened. The open limit is a full physics suite on a variable mesh: Thompson's lookup tables were not here, and radiation on scattered active cells has no timing yet.
