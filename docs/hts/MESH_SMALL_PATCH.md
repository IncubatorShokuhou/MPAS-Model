# Small fine-patch mesh (phase 2c)

The public variable-resolution meshes that fit this VM (`x5.30210`, `x4.10242`) put about half of their cells in the fine region. The large-step tendency ceiling on those meshes is 1.3–1.5×, and phase 2b did not beat all-fine. This mesh is the other shape: one circular fine patch, a few percent of the cells, and a spacing ratio near 10×.

The grid NetCDF is about 30 MB and is not in git. Regenerate it with the commands below. The 4-rank partition and a small metadata file are committed:

- `docs/hts/partitions/sp8.24755.graph.info.htsw.part.4`
- `docs/hts/meshes/sp8.24755.meta.json`

## How it was built

Generator: `mpas_tools.mesh.creation.build_mesh.build_spherical_mesh`, which is the JIGSAW path documented at <https://mpas-dev.github.io/MPAS-Tools/master/mesh_creation.html>. The driver in this repo is `docs/hts/scripts/build_small_patch_mesh.py`.

`mpas_tools` and `jigsawpy` are on conda-forge, not PyPI. On this VM:

```bash
# micromamba 2.x, linux-64, from https://micro.mamba.pm/api/micromamba/linux-64/latest
micromamba create -p "$HOME/micromamba/envs/mpas" -c conda-forge \
    python=3.11 mpas_tools jigsawpy netcdf4 numpy
export PATH="$HOME/micromamba/envs/mpas/bin:$PATH"
python docs/hts/scripts/build_small_patch_mesh.py --out-dir /tmp/sp8
```

The packages that produced the mesh used below were `mpas_tools` 2.1.0, `jigsawpy` 0.3.3, and JIGSAW 0.9.14 (`jigsaw` and `MpasMeshConverter.x` on that env's `PATH`). JIGSAW plus the converter finished in a few seconds. A full geographic static interpolation (`mpas_static`, case 7 with `config_static_interp = true`) was not run; this mesh is for the dry Jablonowski–Williamson init.

The cell-width field is a circular patch. Inside 180 km of (lat, lon) = (30°, 0°) the requested spacing is 20 km. It blends in log(dx), with a smoothstep, across the next 700 km, out to 160 km. The sample grid is 0.5°. `build_spherical_mesh` is called with `plot_cellWidth=False` (no cartopy coastline).

Two post-steps are required before `init_atmosphere` will accept the file.

1. `MpasMeshConverter.x` writes NetCDF-4 and sets `meshDensity = 1`. This build's SMIOL library cannot read NetCDF-4 (`NetCDF: Attempt to use feature that was not turned on when netCDF was built`). The script rewrites the file as `NETCDF3_64BIT_OFFSET`, the same format as `x1.40962.grid.nc`.
2. HTS levels are `nint(log2(max(meshDensity^(-1/4), 1)))`. With density left at 1, every cell is level 0. The script sets `meshDensity = (dx_min / dx)^4` with `dx = sqrt(areaCell)`, while the coordinates are still in meters. That is the NCAR convention: density is 1 on the finest cells.
3. The converter writes coordinates in meters and `sphere_radius = 6371229`. `init_atmosphere` multiplies stored coordinates by 6371229 m, which is correct for the unit-sphere NCAR grids and wrong for a mesh that is already in meters (areas grow by `R^2`, winds collapse). The script divides `x/y/z` of cells, edges, and vertices, and `dcEdge` and `dvEdge`, by `R`, divides areas by `R^2`, and sets `sphere_radius = 1`. `latCell`, `lonCell`, `meshDensity`, `angleEdge`, and `gridSpacing` are left as the converter wrote them.

The run that was integrated is `sp8.24755` (24755 cells, 74259 edges, 49506 vertices). Realized `sqrt(areaCell)` is 16.7–169 km. `(meshDensity_max / meshDensity_min)^0.25 = 10.08` against a requested ratio of 8. JIGSAW did not hit 20 km and 160 km exactly. There were 0 obtuse triangles.

## Cells per level

Unclamped levels, `L = nint(log2(max(dx/dx_min, 1)))`:

| level | cells | fraction | spacing, roughly |
| --- | ---: | ---: | --- |
| 0 | 885 | 3.58% | ~17–24 km |
| 1 | 452 | 1.83% | transition |
| 2 | 150 | 0.61% | transition |
| 3 | 23268 | 94.0% | ~150–170 km |

The fine patch is a disk, not a latitude band. Level 0 is 3.6% of the cells. The transition rings are smaller than the patch. Almost all of the mesh is the coarse level.

Large-step tendency work over one top-level interval, relative to taking every cell at the fine step. Acoustic work does not shrink, because every level keeps the fine acoustic `dts`.

| `config_hts_max_level` | cells by level | large-step work / all-fine | ceiling |
| --- | --- | ---: | ---: |
| 1 (2 levels) | 885 + 23870 | 0.518 | 1.93× |
| 2 (3 levels) | 885 + 452 + 23418 | 0.281 | 3.55× |
| 3 (4 levels) | 885 + 452 + 150 + 23268 | 0.164 | 6.10× |

The ceiling is only the masked slow tendency. Diagnostics stay full-domain, and the acoustic step stays the fine `dts` on every level (a 60 s interface acoustic step diverged in phase 2). The 6× figure is an upper bound on `atm_compute_dyn_tend`, not on `time integration`.

## Partition

```bash
docs/hts/scripts/hts_partition.py \
    sp8.24755.grid.nc sp8.24755.graph.info 4 --max-level 3
```

`gpmetis -ptype=kway` on a 4-constraint graph. The committed file is the result. Rank cell counts are 5902, 6367, 6122, 6364. Level 0 counts on those ranks are 226, 225, 207, 227. Edgecut 1305. One rank's partition has two components. The same file was used for the all-fine run and for every HTS clamp, so the comparison is not a change of decomposition.

## Init

JW case 2, 26 levels, one MPI rank (`init_atmosphere_model`). A 4-rank init on the double-scaled (meters stored as if they were the unit sphere) grid aborted in malloc; after the unit-sphere conversion, the 1-rank init exits 0. `u` is ±35 m/s, `areaCell` is 2.80×10^8 to 2.84×10^10 m², and `dcEdge` is 16.0–210 km. The grid has no terrain and no land use (`landmask` is 0, `ter` is 0). It is an aquaplanet JW mesh, not a real-data static file.
