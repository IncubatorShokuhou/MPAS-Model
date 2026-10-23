#!/usr/bin/env python3
"""Build a global MPAS mesh with one circular fine patch.

This is the Phase 2c speed-case mesh. It uses the MPAS-Tools / JIGSAW path
(`mpas_tools.mesh.creation.build_mesh.build_spherical_mesh`), which is the
mesh generator documented at
https://mpas-dev.github.io/MPAS-Tools/master/mesh_creation.html

JIGSAW's triangle file stores meshDensity = 1. HTS levels follow
meshDensity^(-1/4), so after the MPAS converter finishes this script rewrites

    meshDensity = (dx_min / dx)^4,   dx = sqrt(areaCell)

which is the same convention as the NCAR variable-resolution grids
(meshDensity is 1 on the finest cells).

Requires the conda-forge tools on PATH (jigsaw, MpasMeshConverter.x) and
mpas_tools + jigsawpy importable. See docs/hts/MESH_SMALL_PATCH.md.

Example:
  PATH=/home/ubuntu/micromamba/envs/mpas/bin:$PATH \\
  /home/ubuntu/micromamba/envs/mpas/bin/python \\
      docs/hts/scripts/build_small_patch_mesh.py --out-dir /tmp/sp8
"""

import argparse
import json
import os
import shutil
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def cell_width_field(lat, lon, dx_fine_km, dx_coarse_km, radius_km, width_km,
                     center_lat, center_lon):
    """Cell width (km) on a lat/lon grid. Log-blended with a smoothstep."""
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64)
    lat2, lon2 = np.meshgrid(lat, lon, indexing="ij")
    # Great-circle distance on a sphere of radius 6371 km.
    radius = 6371.0
    phi1 = np.deg2rad(lat2)
    phi2 = np.deg2rad(center_lat)
    dlam = np.deg2rad(lon2 - center_lon)
    hav = (np.sin((phi1 - phi2) / 2.0) ** 2
           + np.cos(phi1) * np.cos(phi2) * np.sin(dlam / 2.0) ** 2)
    dist = 2.0 * radius * np.arcsin(np.minimum(1.0, np.sqrt(hav)))
    s = np.clip((dist - radius_km) / width_km, 0.0, 1.0)
    s = s * s * (3.0 - 2.0 * s)
    log_dx = (1.0 - s) * np.log(dx_fine_km) + s * np.log(dx_coarse_km)
    return np.exp(log_dx)


def _to_classic(src_path, dst_path):
    src = Dataset(src_path, "r")
    dst = Dataset(dst_path, "w", format="NETCDF3_64BIT_OFFSET")
    for name, dim in src.dimensions.items():
        dst.createDimension(name, None if dim.isunlimited() else len(dim))
    for name, attr in src.__dict__.items():
        dst.setncattr(name, attr)
    for name, var in src.variables.items():
        fill = None
        if "_FillValue" in var.ncattrs():
            fill = var.getncattr("_FillValue")
        out = dst.createVariable(name, var.dtype, var.dimensions, fill_value=fill)
        for attr in var.ncattrs():
            if attr == "_FillValue":
                continue
            out.setncattr(attr, var.getncattr(attr))
        out[:] = var[:]
    src.close()
    dst.close()


def to_unit_sphere(path, earth_radius):
    """Match the NCAR grids: coordinates on the unit sphere, sphere_radius = 1.

    init_atmosphere multiplies the stored coordinates by 6371229 m. A mesh
    that is already in meters and carries sphere_radius = 6371229 gets scaled
    a second time.
    """
    length = (
        "xCell", "yCell", "zCell", "xEdge", "yEdge", "zEdge",
        "xVertex", "yVertex", "zVertex", "dcEdge", "dvEdge",
    )
    area = ("areaCell", "areaTriangle", "kiteAreasOnVertex")
    nc = Dataset(path, "a")
    radius = float(earth_radius)
    for name in length:
        if name in nc.variables:
            var = nc.variables[name]
            var[:] = np.asarray(var[:], dtype=np.float64) / radius
    for name in area:
        if name in nc.variables:
            var = nc.variables[name]
            var[:] = np.asarray(var[:], dtype=np.float64) / radius ** 2
    nc.sphere_radius = np.float64(1.0)
    nc.on_a_sphere = "YES"
    nc.close()


def rewrite_mesh_density(path):
    """Set meshDensity from cell area. Returns dx stats in km and the array."""
    nc = Dataset(path, "a")
    area = np.asarray(nc.variables["areaCell"][:], dtype=np.float64)
    # areaCell is m^2 when sphere_radius is in meters.
    dx = np.sqrt(np.maximum(area, 0.0))
    dx_min = float(dx.min())
    density = (dx_min / dx) ** 4
    var = nc.variables["meshDensity"]
    var[:] = density
    nc.setncattr(
        "hts_mesh_density",
        "meshDensity=(dx_min/dx)^4 with dx=sqrt(areaCell); finest cell is 1",
    )
    nc.close()
    return dx_min, dx, density


def level_counts(density, max_level):
    spacing = np.maximum(np.asarray(density, dtype=np.float64) ** (-0.25), 1.0)
    level = np.rint(np.log(spacing) / np.log(2.0)).astype(np.int32)
    level = np.clip(level, 0, max_level)
    nlev = int(level.max()) + 1
    counts = np.bincount(level, minlength=nlev)
    tend = 0.0
    n = float(level.size)
    for ell in range(nlev):
        tend += (counts[ell] / n) / (2 ** ell)
    return counts, tend


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--dx-fine-km", type=float, default=20.0)
    p.add_argument("--dx-coarse-km", type=float, default=160.0)
    p.add_argument("--radius-km", type=float, default=180.0)
    p.add_argument("--width-km", type=float, default=700.0)
    p.add_argument("--center-lat", type=float, default=30.0)
    p.add_argument("--center-lon", type=float, default=0.0)
    p.add_argument("--dlat", type=float, default=0.5)
    p.add_argument("--dlon", type=float, default=0.5)
    p.add_argument("--earth-radius", type=float, default=6371229.0)
    args = p.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    work = out / "jigsaw_work"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir()
    os.chdir(work)

    lat = np.linspace(-90.0, 90.0, int(round(180.0 / args.dlat)) + 1)
    lon = np.linspace(-180.0, 180.0, int(round(360.0 / args.dlon)) + 1)
    cell_width = cell_width_field(
        lat, lon, args.dx_fine_km, args.dx_coarse_km, args.radius_km,
        args.width_km, args.center_lat, args.center_lon,
    )

    from mpas_tools.mesh.creation.build_mesh import build_spherical_mesh

    build_spherical_mesh(
        cell_width, lon, lat, args.earth_radius,
        out_filename="base_mesh.nc", plot_cellWidth=False,
    )
    grid_src = work / "base_mesh.nc"
    graph_src = work / "graph.info"
    if not graph_src.is_file():
        raise SystemExit(f"MpasMeshConverter did not write {graph_src}")

    dx_min, dx, density = rewrite_mesh_density(grid_src)
    # MPAS SMIOL in this build reads classic 64-bit-offset files. The
    # converter writes NetCDF-4, which init_atmosphere rejects.
    classic = work / "base_mesh.classic.nc"
    _to_classic(grid_src, classic)
    to_unit_sphere(classic, args.earth_radius)
    grid_src = classic
    # Unclamped levels, then the ceilings the model will see at each cap.
    counts_raw, _ = level_counts(density, max_level=8)
    ceilings = {}
    for cap in range(0, len(counts_raw)):
        counts, tend = level_counts(density, max_level=cap)
        ceilings[str(cap)] = {
            "counts": counts.tolist(),
            "fine_fraction": float(counts[0] / counts.sum()),
            "large_step_work_over_all_fine": float(tend),
            "large_step_ceiling": float(1.0 / tend),
        }

    stem = f"sp8.{int(dx.size)}"
    grid_out = out / f"{stem}.grid.nc"
    graph_out = out / f"{stem}.graph.info"
    shutil.copy(grid_src, grid_out)
    shutil.copy(graph_src, graph_out)

    meta = {
        "generator": "mpas_tools.mesh.creation.build_mesh.build_spherical_mesh",
        "jigsaw": "conda-forge jigsawpy / jigsaw, via MpasMeshConverter.x",
        "dx_fine_km": args.dx_fine_km,
        "dx_coarse_km": args.dx_coarse_km,
        "requested_ratio": args.dx_coarse_km / args.dx_fine_km,
        "radius_km": args.radius_km,
        "transition_width_km": args.width_km,
        "center_lat": args.center_lat,
        "center_lon": args.center_lon,
        "earth_radius_m": args.earth_radius,
        "n_cells": int(dx.size),
        "dx_sqrt_area_km": {
            "min": float(dx.min() / 1000.0),
            "p50": float(np.median(dx) / 1000.0),
            "max": float(dx.max() / 1000.0),
        },
        "mesh_density_min": float(density.min()),
        "mesh_density_max": float(density.max()),
        "realized_ratio_(mdmax/mdmin)^0.25": float((density.max() / density.min()) ** 0.25),
        "unclamped_level_counts": counts_raw.tolist(),
        "ceilings_by_max_level": ceilings,
        "grid": str(grid_out),
        "graph": str(graph_out),
    }
    (out / f"{stem}.meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
