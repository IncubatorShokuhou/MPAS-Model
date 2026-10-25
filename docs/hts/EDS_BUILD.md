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
