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
- Job 199953 (`bb575a1b2` open-interface acoustic parent, `hts_ifix2h`, rebuild 199950, GWDO off) regressed further. all-fine wall 374.34 s (critical 0). HTS stayed near ±88 m/s through 00:11, grew to about ±250 by 00:15:30, then at 00:17:30 `u` was about ±7.2×10⁵ m/s and `w` about −1×10⁴ / +8×10⁴ m/s. First NaN is `u` at 00:18:00. SIGSEGV exit 11, wall 71.02 s. Radiation off.
- Job 199977 (`cf1dc652d` held-ring `rho`/`theta` ramp, `hts_ramp2h`, rebuild 199974, GWDO off, convection off) regressed again. all-fine wall 357.09 s, TI 353.46 s (critical 0, one non-critical error). HTS stays near ±88 m/s through 00:09, spikes to ±98–108 at 00:09:30, `|u|` about 300 then ~134 (`|w|` about 12) at 00:10:30, `|u|` about 183 at 00:11:00, `u` ±330–532 and `w` ±17–26 at 00:11:30, `u` ±227–274 at 00:12:00. First NaN is `u` at 00:12:30 (`w` still finite). Both NaN from 00:13:30. SIGSEGV exit 11 / mpirun 139, wall 57.63 s. Radiation off. Logs: `work/hts_land_ramp2h/` and `work/slurm/hts_ramp2h.199977.out`. The writeup does not name a new hotspot cell. CPAS interpolates a full parent state. Ramping only `rho` and `theta` while pinning `u` and `w` is not that coupling, and dry 6 h does not validate it. The acoustic-parent update is out. `45dfe9329` is the pure hold, not the `615836256` thermo-only ramp.
- Job 199996 (pure hold `20b4f3a34`, 2 h) reproduced the `d3fe6fb75` ceiling: all-fine OK, HTS first `w` NaN at 01:57:30 after a `u` spike. Relative to `d3fe6fb75`, `cf1dc652d` is three stacked changes: the `902a6084d` fine-dt subcycle (job 199900, `u` NaN at 00:45; removed again inside `cf1dc652d`), the `bb575a1b2` acoustic-parent update (job 199953, `u` NaN at 00:18), and a held-ring `rho`/`theta` ramp on top of that update (job 199977, `u` NaN at 00:12:30). The current tree has none of those. It keeps the pure-hold close-step, blends the diagnostic ring across the previous completed `rho`, `theta`, `u`, and `w` increment, and holds an active-level interface edge on the coarser side of that increment so `tend_u_euler` does not integrate it. A parent integration of the same edge NaNed dry `w` at 00:42 and is not in this tree. A local dry 6 h of the held edge had 0 critical errors and no NaN (`|u|` −35.081 / 35.044 m/s at 01:57, −35.026 / 35.018 m/s at 6 h). Eds has not run it. The 01:57 land NaN is not claimed gone. A 6 h multi-level land run waits until that ≥2 h pair is finite.
