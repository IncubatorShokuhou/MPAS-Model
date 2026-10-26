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
- Multi-level (`config_hts_max_level = 3`) NaNed `u` near 00:25 with the fine `dt_pbl` (jobs 199635, 199659, 199660). GWDO was not the cause. `dt_pbl` now follows the level step. `max_level = 0` (job 199662, 45 min, bit-identical) was the workaround and is not required. The ≥2 h / ≥6 h confirmation is a later eds rerun.
