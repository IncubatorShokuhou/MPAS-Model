# EDS build notes (MPAS-A HTS)

Host: `eds` (`en_op01@116.114.29.62:12223`), path `/eds/energy/dev01/lvhao/MPAS-Model`, branch `hts`.

## Compiler / MPI / I/O stack (reuse tutorial local install)

No environment modules on this host (`module: command not found`). The working stack is the one built by `MPAS_tutorial/install_io_stack.sh` + `install_mpas.sh` (not system MPI).

| item | path / value |
| --- | --- |
| GCC / GFortran | system `/usr/bin/gcc`, `/usr/bin/gfortran` |
| OpenMPI | `/eds/energy/dev01/lvhao/MPAS_tutorial/local/io-stack` (4.1.8) |
| NetCDF-C / Fortran | same prefix (C ~4.9.x via nc-config; Fortran 4.6.2) |
| PnetCDF | 1.15.0, same prefix |
| PIO | static `libpioc.a` / `libpiof.a` in same prefix (used by tutorial CMake MPAS 8.4.2 install) |
| Baseline MPAS 8.4.2 (no HTS) | `/eds/energy/dev01/lvhao/MPAS_tutorial/local/mpas-8.4.2/bin/{mpas_atmosphere,mpas_init_atmosphere}` |

Environment for HTS Makefile build (SMIOL + PnetCDF, matches cloud Phase-1 `build_mpas.sh` style; single precision default; OpenMP off):

```bash
IO=/eds/energy/dev01/lvhao/MPAS_tutorial/local/io-stack
export PATH="$IO/bin:$PATH"
export LD_LIBRARY_PATH="$IO/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /eds/energy/dev01/lvhao/MPAS-Model
make -j$(nproc) gfortran CORE=atmosphere NETCDF=$IO PNETCDF=$IO
make -j$(nproc) gfortran CORE=init_atmosphere NETCDF=$IO PNETCDF=$IO
```

Alternate (tutorial CMake + PIO), for comparison only:

```bash
cmake -S . -B work/hts_build/cmake-pio \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=$IO/bin/mpicc \
  -DCMAKE_Fortran_COMPILER=$IO/bin/mpifort \
  -DCMAKE_PREFIX_PATH=$IO -DPnetCDF_ROOT=$IO -DPIO_ROOT=$IO -DNetCDF_ROOT=$IO \
  -DMPAS_CORES="atmosphere;init_atmosphere" \
  -DMPAS_USE_PIO=ON -DMPAS_OPENMP=OFF
cmake --build work/hts_build/cmake-pio --parallel $(nproc)
```

Cloud HTS timings used the **Makefile / SMIOL / single-precision** path. Prefer that for apples-to-apples.

## Mesh tools (read-only tutorial)

- Script: `docs/hts/scripts/build_small_patch_mesh.py` → `sp8.24755` style VR patch.
- Conda env (already present): `/eds/energy/dev01/lvhao/MPAS_tutorial/work/mesh_lesson/root/local/miniforge3/envs/mpas_mesh` (`mpas_tools`, `jigsaw`, `MpasMeshConverter.x`).
- Existing VR grid (copy OK): `MPAS_tutorial/work/beijing_92-25km/beijing.x4.163842.grid.nc` (large; prefer regenerating small `sp8` under `work/hts_mesh/`).
- WPS_GEOG ACL often blocks static interp; JW idealized init does not need GEOG.

## Work dirs (writable, under this repo only)

- `work/hts_build/` — build logs / CMake tree
- `work/hts_mesh/` — grid / graph / partition copies
- `work/hts_smoke/` — short JW / HTS vs all-fine runs

## Safety

Write/delete only under `/eds/energy/dev01/lvhao/MPAS-Model`. Tutorial, GFS, WPS_GEOG are read-only sources.

## Measured on eds (2026-10-24, America/Los_Angeles session; host clock CST+8)

| item | value |
| --- | --- |
| gcc / gfortran | 11.4.0 (Ubuntu 22.04) |
| OpenMPI | 4.1.8 (`MPAS_tutorial/local/io-stack`) |
| NetCDF-C | 4.9.3 |
| NetCDF-Fortran | 4.6.2 |
| PnetCDF | 1.15.0 |
| PIO | present (static); HTS Makefile build uses **SMIOL** (`-DMPAS_SMIOL_SUPPORT`), not PIO |
| CPUs | 64 logical |
| RAM | ~1 TiB |
| Mesh env | `MPAS_tutorial/work/mesh_lesson/.../envs/mpas_mesh` |

## Mesh status

Regenerated under `work/hts_mesh/`:

- `sp8.24755.grid.nc` (30 MB, NETCDF3_64BIT, unit sphere)
- `sp8.24755.graph.info`
- `sp8.24755.graph.info.htsw.part.4` (METIS 4-constraint, max-level 3)
- `sp8.24755.meta.json`

Same cell counts as cloud docs: 24755 cells, levels [885, 452, 150, 23268], large-step ceiling ~6.10×.

## Slurm (mandatory for heavy work)

The login node is only for light commands (ls, git, editing). Builds, mesh generation, init and model runs go through `sbatch`.

| item | value |
| --- | --- |
| account | `energy` (`sacctmgr show assoc`) |
| QOS | `normal_qos` (also `high_priority`) |
| partitions | `normal` (default; node27, node32–37; 64 or 50 CPUs, ~350–500 GB), `gpu4090` (node53), `a100` (node54) |
| MPI launch | `mpirun` from `MPAS_tutorial/local/io-stack` inside the allocation (project OpenMPI 4.1.8) |
| job scripts / logs | `work/slurm/*.sbatch`, `work/slurm/%x.%j.out` |

On 2026-10-24 the `normal` partition was nearly full with other `en_op01` array jobs (≈2 free CPUs per node), so small jobs wait in `Priority`/`Resources`. `gpu4090` and `a100` were idle; using them for CPU-only runs needs the user's OK.

Note: the `atmosphere_model` and `init_atmosphere_model` builds on 2026-10-24 ran on the login node before this rule was set. Rebuild through Slurm from now on (template: `work/slurm/hts_smoke.sbatch`; swap the body for the `make` lines above).

## Build fixes made on eds

1. `src/core_atmosphere/dynamics/Makefile`: added `mpas_atm_dissipation_models.o: mpas_atm_hts.o`. Without it, a parallel `make -j32` races and fails with `Cannot open module file 'mpas_atm_hts.mod'`.
2. `src/core_atmosphere/tools/manage_externals/manic/repository_git.py`: `_git_remote_add` now skips a remote that already exists, so re-running `make` does not fail with `remote NCAR_MMM-physics.git already exists`.
3. `src/core_atmosphere/Externals.cfg`: GitHub URLs go through `https://gh-proxy.org/` because eds has no direct GitHub access. This change is for eds only; do not merge it upstream.

## Build result (2026-10-24)

- `atmosphere_model` (8.2 MB; HTS; SMIOL; single precision; mpi_f08; no OpenMP); copy kept at `work/hts_build/atmosphere_model`
- `init_atmosphere_model` (2.5 MB)
- Logs: `work/hts_build/atmosphere_make.log`, `work/hts_build/init_make.log`

## Dry 6 h on eds (sp8.24755, 4 MPI, dt=30 s)

Job `work/slurm/hts_dry6h.sbatch` on `gpu4090` / node53, finished `2026-10-25T00:10:49+08:00`. Binary from Makefile SMIOL single-precision build at commit `1b0264ca` (doc-only `84b4956f` was tip of branch afterward). Critical errors 0.

| case | wall | time integration | ratio (all-fine / HTS on TI) |
| --- | ---: | ---: | ---: |
| all-fine | 133.96 s | 127.41 s | — |
| HTS | 71.80 s | 65.70 s | **1.94×** |

Cloud Phase-1 dry 6 h on the same mesh was **2.00×** (125.98 / 62.95 s TI).

History compare (all time levels max abs diff HTS vs all-fine):

| field | max\|Δ\| | rms |
| --- | ---: | ---: |
| u | 4.29×10⁻¹ | 8.37×10⁻³ |
| w | 1.71×10⁻² | 2.88×10⁻⁴ |
| theta | 9.44×10⁻¹ | 1.50×10⁻² |
| rho | 9.35×10⁻⁵ | 3.96×10⁻⁶ |
| pressure | 9.95 | 3.98×10⁻¹ |

Last time level alone: u 3.94×10⁻¹, w 1.71×10⁻², theta 9.44×10⁻¹, rho 9.35×10⁻⁵, pressure 1.31. Outputs under `work/hts_dry6h/{init,allfine,hts}/`. Short smoke (20 steps) also finished earlier with 0 critical errors (`work/hts_smoke/`).

## Land VR static (personal WPS_GEOG)

Shared `/eds/shared/dm/WPS_GEOG` is not writable/readable for init (ACL). Use the personal tree:

- Path: `/eds/energy/dev01/lvhao/WPS_GEOG/` (trailing slash required by MPAS)
- Contents: `topo_gmted2010_30s`, `modis_landuse_20class_30s_with_lakes`, `soiltype_{top,bot}_30s`, `albedo_modis`, `maxsnowalb_modis`, `greenfrac_fpar_modis`, `soiltemp_1deg`, native GWD (`varsso*`, `orogwd_*`)
- Download jobs: 199006, 199019 (2026-10-25)

Static on `sp8.24755` (job **199304**, ~18 min wall on `gpu4090`):

- Work dir: `work/hts_land_static/`
- Output: `work/hts_land_static/sp8.24755.static.nc` (39 MB)
- Namelist: `config_static_interp=true`, `config_native_gwd_static=true`, MODIS/STATSGO/GMTED2010, `config_noahmp_static=false` (no `soilgrids/` dataset here)
- Sanity: 0 critical / 0 error; `landmask` mean ≈0.32; `ter` ≈[-27, 5336] m; `var2d` present for GWDO

Real-land IC + short forecast (jobs **199309** init, **199311** atmosphere):

- IC: GFS ungrib `FILE:2026-10-15_00` (symlink from tutorial Beijing case) → `work/hts_land_init/sp8.24755.init.nc` (210 MB)
- Physics: Noah + YSU + RRTMG + Kessler, convection/GWDO off, `dt=30` s, 15 min, 4 MPI
- `time integration`: all-fine **46.37 s**, HTS **29.55 s** → **1.57×** (wall 49.05 / 31.59 s)
- Accuracy (t=15 min): t0 identical; land `skintemp` rms≈0.40 K, land `t2m` rms≈0.34 K; domain `u` max|Δ|≈3.8 m/s. Critical errors 0.
- Pointer: `work/README_LAND_VR.md`

Do not write outside `/eds/energy/dev01/lvhao/MPAS-Model` and the personal GEOG path. All heavy work via Slurm (`energy` / `normal_qos`; prefer `gpu4090` when `normal` is full).

## Multi-level land stability

The 15-minute pair above ends before the failure. Same mesh, IC, `dt = 30` s, Noah+YSU+RRTMG+Kessler, multi-level (`config_hts_max_level = 3`):

| job | max_level | GWDO | duration | HTS |
| --- | ---: | --- | --- | --- |
| 199311 | 3 | off | 15 min | completed (before the fail window) |
| 199635 / 199660 | 3 | on | 6 h / 45 min | NaN in `u` near 00:25, then SIGSEGV |
| 199659 | 3 | off | 6 h | same NaN |
| 199662 | 0 | on | 45 min | completed, 0 NaNs, bit-identical to all-fine |
| 199818 | 3 | off | 2 h | `d3fe6fb75` only: finite through ~1.5 h, first NaN in `w` at `2026-10-15_01:57:30`. all-fine critical 0. No SIGSEGV on this 2 h job |
| 199819 | 3 | on | 6 h | same first NaN on `d3fe6fb75`, then SIGSEGV. GWDO does not move the clock |
| 199900 | 3 | off | 2 h | `902a6084d` subcycle, rebuild 199893: `|u|` jumps from about −98 to −269 / +272 at 00:39:30 (level-3 step, `hts_tick` 79), stuck near ±269 until the next level-3 step at 00:43:30, then first `u` NaN at `2026-10-15_00:45:00` (`hts_tick` 90, level 0 only) and `w` after (wall 165.30 s, exit 11, SIGSEGV rank 3). Hotspot cell 22181 (31.84°N, 0.14°E, land, `meshDensity` 0.610). 0.5 h frame still finite (`u` relL2 1.68×10⁻², `w` 3.89, max `|Δw|` 18.85 m/s, land `skintemp` rms 0.27 K). 0 pressure inversions. all-fine 2 h OK (wall 373.60 s, critical 0) |

Onset on the runs that still used a fine `dt_pbl`: t = 0 radiation is fine; near 00:19, `pressure(1) < pressure(2)` on southern-ocean coarse cells; 00:25 `u` is NaN. The second radiation call is at 00:30, after the state is already corrupt. One-level integration was 147.4 s all-fine / 145.5 s HTS. Giving those schemes the level `dt_pbl` removed that 00:25 NaN. Jobs 199818 and 199819 are the pair after that correction.

`01:57:30` is a fine tick (`hts_tick` 235), not a radiation alarm. Level 3 last stepped at 01:55:30. The failing tick steps levels 0, 1, and 2. The log has 0 pressure-inversion warnings, and the last radiation call is 01:30. At 1 h the largest differences sit on land cells near 32°N, 0°E with `meshDensity` about 0.38–0.61 (level 0) beside a level-1 neighbor near 0.12. The level-0 step drives the interface-edge `|u|` to about 243 m/s (all-fine stays near 89). The level-1 impose then replaces that edge with the mean flux, about 47 m/s, and the level-2 step NaNs `w`.

The fine-dt PBL/Noah subcycle from `902a6084d` is removed in `cf1dc652d`. Job 199900 showed it moves the land-VR NaN earlier (00:45, first field `u`, after `|u|` reached about ±269 m/s on the 00:39:30 level-3 step) than the single level-dt call (01:57 on `d3fe6fb75`, first field `w`). The 0.5 h frame is still finite and the largest `|Δw|` (18.85 m/s) is on cell 22181, the same refined-land cell as that later blow-up. Skin temperature at 0.5 h is only 0.27 K rms on land. That subcycle is not how CPAS runs multi-rate physics. CPAS uses a WRF-style level step: PBL and Noah are called once with `Δt = config_dt * 2^L`, and the fine step reads a coarse state that moves across the interval. A coarse column is one `physics_driver` call with `dt_pbl = dt_level` again. Kessler and WSM6 are one call of that dt. Radiation stays a single K/s call. The 01:57 edge was integrating `tend_u_euler` against a neighbor whose `rho` and `theta` were frozen and whose acoustic `rtheta_pp` was zero. Held cells in the active level's diagnostic ring now follow `snap + alpha * previous_increment` in `rho` and `theta` only, `alpha = mod(hts_tick, 2^L) / 2^L`. Interior coarse columns stay on the snapshot. `u` and `w` stay on the snapshot, and the ramp is removed before the coarse RK. A coarser cell that shares an edge with the active level still updates `rho_pp`, `rtheta_pp`, and `rw_p` on that step. `config_hts_max_level = 0` is not the land-physics workaround. The ≥2 h rerun of this combination is still to be done on eds; this note does not record it and does not claim the land case is fixed. A 6 h multi-level land run waits until that pair is finite.

The land static and GFS IC are not on the cloud VM. A 45-minute extension of the local USGS-water RRTMG/YSU/Noah/Kessler smoke (the 8-minute case, not the eds land VR) finished on both the pre-`dt_pbl` binary and the level-step binary: critical 0, no NaNs, no `pressure(1) < pressure(2)`. At 45 minutes, level-3 lowest-layer `theta` differs from all-fine by a mean of −0.073 K with the level `dt_pbl` (maximum 0.28 K) and −0.046 K with the fine `dt_pbl` (maximum 0.23 K). Dry 6 h on the subcycled binary, before the interface-parent update, matches the previous HTS history bit for bit (`u`, `w`, `theta`, `rho`, `pressure`). Kessler 1 h does too on this IC. After the interface-parent update those two moved and stay inside the published all-fine table: dry `u` / `w` / `theta` / `rho` relative L2 6.95×10⁻⁴ / 0.404 / 3.05×10⁻⁵ / 3.94×10⁻⁶, Kessler 2.67×10⁻⁴ / 0.259 / 3.71×10⁻⁶ / 6.04×10⁻⁶. Integration on that sitting was 64.45 s and 11.74 s; quote the published 62.95 s and 11.17 s. A 2 h USGS-water multi-level run with the subcycle also finished, critical 0, no NaNs and no lowest-layer pressure inversions in the 30-minute frames. Integration on that pair was 132.9 s HTS against 74.1 s all-fine, and `u` relative L2 at 2 h was 8.0×10⁻⁴. That static has zero soil moisture and does not reproduce the land-VR NaN. The physics calls on coarse columns are back, so the 1.3×–1.57× land timing from the single long step will not carry over.
