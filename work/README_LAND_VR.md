# Land VR research prototype (eds)

Personal WPS GEOG (do not use shared ACL path):

- `/eds/energy/dev01/lvhao/WPS_GEOG/` (trailing slash in namelist)
- Downloaded 2026-10-25, Slurm jobs 199006 + 199019 (~50G)

## Outputs

| path | what |
| --- | --- |
| `work/hts_land_static/sp8.24755.static.nc` | real static on VR mesh (job 199304) |
| `work/hts_land_init/sp8.24755.init.nc` | GFS 2026-10-15_00 real IC (job 199309 init) |
| `work/hts_land_run/{allfine,hts}/` | 15 min Noah+YSU+RRTMG+Kessler HTS vs all-fine (job 199311) |
| `work/slurm/hts_land_static.sbatch` | static generator |
| `work/slurm/hts_land_init_run.sbatch` | init + paired atmosphere runs |

## Notes

- `config_noahmp_static = false` (no `soilgrids/` in personal GEOG; atmosphere uses `sf_noah`).
- Ungrib intermediates reused from tutorial: `MPAS_tutorial/work/beijing_240km_2026101500/FILE:2026-10-15_00`.
- Benign atmosphere ERROR: diagnostics stream not in simplified `streams.atmosphere` (critical = 0).
- Multi-level (`config_hts_max_level = 3`) NaNed `u` near 00:25 with the fine `dt_pbl` (jobs 199635, 199659, 199660). GWDO was not the cause. `dt_pbl` follows the level step. `max_level = 0` (job 199662, 45 min, bit-identical) was the workaround and is not required.
- Jobs 199818 (2 h, GWDO off) and 199819 (6 h, GWDO on) are `d3fe6fb75`. Both NaN `w` at `2026-10-15_01:57:30` on land cells near 32°N, 0°E (level-0 interface, `|u|` about 243 m/s).
- Job 199900 (rebuild 199893, `902a6084d` fine-dt PBL/Noah subcycle, GWDO off) regressed versus `d3fe6fb75`. `|u|` stays near ±88 m/s through 00:32, jumps to −269 / +272 m/s at 00:39:30 (level-3 step, `hts_tick` 79), stays there until 00:43:30, and the first NaN is `u` at `2026-10-15_00:45:00` (`w` after). all-fine wall 373.60 s (critical 0); HTS wall 165.30 s, exit 11. 0 pressure inversions. At 0.5 h the state is still finite (`w` relL2 3.89, land `skintemp` rms 0.27 K, land `t2m` rms 0.43 K) and the hotspot is cell 22181 (31.84°N, 0.14°E, land, `meshDensity` 0.610, max `|Δw|` 18.85 m/s), the same neighborhood as the 01:57 cells.
- Job 199953 (`bb575a1b2` open-interface acoustic parent, `hts_ifix2h`) regressed further. all-fine wall 374.34 s (critical 0). HTS `|u|` was about 1000 m/s by 00:17:30, first `u` NaN at 00:18:00, SIGSEGV exit 11, wall 71 s. Updating `rho_pp` and `rtheta_pp` on the parent without a prognostic step is worse than the pure hold. The diagnostic-ring thermodynamic ramp was not run on eds and is removed with that parent update. This tree is the `d3fe6fb75` coupling again: level `dt_pbl`, coarse state on the snapshot. eds has not re-run ≥2 h on the restore. The 01:57 NaN is not claimed gone. A 6 h multi-level land run waits until that ≥2 h pair is finite.
