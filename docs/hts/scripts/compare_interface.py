#!/usr/bin/env python3
"""Interface-band comparison for hierarchical time stepping.

Cell levels follow mpas_atm_hts: L = nint(log2(meshDensity**(-1/4))),
clamped to --max-level. An interface cell has a neighbor at a different
level. Interior cells are everyone else, split into fine (level 0) and
coarse (level > 0). Distance is the cell-graph distance to the nearest
interface cell. Horizontal acoustic CFL uses dcEdge and
cs = sqrt(1.4 * pressure / rho) from the reference frame.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_histories import dry_mass, frame, linf, rel_l2


def cell_levels(mesh_density, max_level):
    md = np.asarray(mesh_density, dtype=np.float64)
    levels = np.zeros(md.shape, dtype=np.int32)
    if md.min() <= 0.0 or md.max() <= md.min() * (1.0 + 1.0e-5):
        return levels
    spacing = np.maximum(md ** (-0.25), 1.0)
    levels = np.rint(np.log(spacing) / np.log(2.0)).astype(np.int32)
    levels = np.clip(levels, 0, max_level)
    return levels


def main():
    p = argparse.ArgumentParser()
    p.add_argument("test")
    p.add_argument("reference")
    p.add_argument("grid", help="grid.nc with cellsOnCell and nEdgesOnCell")
    p.add_argument("--frame", type=int, default=-1)
    p.add_argument("--ref-frame", type=int, default=None)
    p.add_argument("--max-level", type=int, default=2)
    args = p.parse_args()

    test = Dataset(args.test)
    ref = Dataset(args.reference)
    grid = Dataset(args.grid)
    nt = len(test.dimensions["Time"])
    nr = len(ref.dimensions["Time"])
    it = args.frame if args.frame >= 0 else nt - 1
    ir = args.ref_frame if args.ref_frame is not None else (args.frame if args.frame >= 0 else nr - 1)

    md = np.asarray(test.variables["meshDensity"][:], dtype=np.float64)
    if md.ndim == 2:
        md = md[0]
    levels = cell_levels(md, args.max_level)
    n_cells = levels.size
    coc = np.asarray(grid.variables["cellsOnCell"][:], dtype=np.int32)
    nec = np.asarray(grid.variables["nEdgesOnCell"][:], dtype=np.int32)
    # NetCDF may store cellsOnCell as (nCells, maxEdges) or (maxEdges, nCells).
    if coc.shape[0] != n_cells and coc.shape[-1] == n_cells:
        coc = coc.T
    if coc.shape[0] != n_cells:
        sys.exit(f"cellsOnCell shape {coc.shape} does not match {n_cells} cells")

    interface = np.zeros(n_cells, dtype=bool)
    for i in range(n_cells):
        for j in range(int(nec[i])):
            nbr = int(coc[i, j]) - 1
            if nbr < 0 or nbr >= n_cells:
                continue
            if levels[nbr] != levels[i]:
                interface[i] = True
                break
    fine_interior = (levels == 0) & (~interface)
    coarse_interior = (levels > 0) & (~interface)
    dist = np.full(n_cells, -1, dtype=np.int32)
    queue = [int(i) for i in np.flatnonzero(interface)]
    dist[interface] = 0
    head = 0
    while head < len(queue):
        i = queue[head]
        head += 1
        for j in range(int(nec[i])):
            nbr = int(coc[i, j]) - 1
            if nbr < 0 or nbr >= n_cells or dist[nbr] >= 0:
                continue
            dist[nbr] = dist[i] + 1
            queue.append(nbr)
    print(f"frames test={it}/{nt} reference={ir}/{nr}  max_level={args.max_level}")
    for L in range(int(levels.max()) + 1):
        print(f"level {L} cells {int(np.sum(levels == L))}")
    print(f"interface cells {int(interface.sum())}  fine interior {int(fine_interior.sum())}  "
          f"coarse interior {int(coarse_interior.sum())}")
    for d in range(0, 4):
        print(f"distance {d} cells {int(np.sum(dist == d))}")
    print(f"distance >=4 cells {int(np.sum(dist >= 4))}")

    # u lives on edges. An interface edge has cellsOnEdge at two different levels.
    coe = np.asarray(grid.variables["cellsOnEdge"][:], dtype=np.int32)
    if coe.shape[-1] == 2 and coe.shape[0] != 2:
        c1 = coe[:, 0] - 1
        c2 = coe[:, 1] - 1
    else:
        c1 = coe[0, :] - 1
        c2 = coe[1, :] - 1
    valid = (c1 >= 0) & (c2 >= 0) & (c1 < n_cells) & (c2 < n_cells)
    edge_interface = np.zeros(c1.shape, dtype=bool)
    edge_interface[valid] = levels[c1[valid]] != levels[c2[valid]]

    def rms(x):
        x = np.asarray(x, dtype=np.float64)
        if x.size == 0:
            return float("nan")
        return float(np.sqrt(np.mean(x ** 2)))

    def report(name, a, b, label, mask):
        aa = a[mask]
        bb = b[mask]
        diff = aa - bb
        print(f"{name:7s} {label:22s} n {aa.size:8d}  relL2 {rel_l2(aa, bb):.6e}  "
              f"Linf {linf(aa, bb):.6e}  rms(diff) {rms(diff):.6e}  "
              f"rms(ref) {rms(bb):.6e}  rms(test) {rms(aa):.6e}")

    cell_regions = (
        ("all", np.ones(n_cells, dtype=bool)),
        ("interface", interface),
        ("fine interior", fine_interior),
        ("coarse interior", coarse_interior),
        ("dist 0", dist == 0),
        ("dist 1 fine", (dist == 1) & (levels == 0)),
        ("dist 1 coarse", (dist == 1) & (levels > 0)),
        ("dist>=2 fine", (dist >= 2) & (levels == 0)),
        ("dist>=2 coarse", (dist >= 2) & (levels > 0)),
    )
    edge_fine = valid & (levels[c1] == 0) & (levels[c2] == 0)
    edge_coarse = valid & (levels[c1] > 0) & (levels[c2] > 0)
    edge_regions = (
        ("all", np.ones(edge_interface.shape, dtype=bool)),
        ("interface", edge_interface),
        ("fine edges", edge_fine),
        ("coarse edges", edge_coarse),
    )

    for name in ("u", "w", "theta", "rho"):
        a = np.asarray(frame(test, name, it), dtype=np.float64)
        b = np.asarray(frame(ref, name, ir), dtype=np.float64)
        if a.shape[0] == n_cells:
            regions = cell_regions
        elif a.shape[0] == edge_interface.size:
            regions = edge_regions
        else:
            print(f"{name:7s} no region mask; leading dimension {a.shape[0]}")
            continue
        for label, mask in regions:
            report(name, a, b, label, mask)

    # Horizontal acoustic CFL. The acoustic step is horizontally explicit
    # and vertically implicit, so the limit is cs * dts / dcEdge.
    if "dcEdge" in grid.variables and "pressure" in ref.variables and "rho" in ref.variables:
        dc = np.asarray(grid.variables["dcEdge"][:], dtype=np.float64)
        if dc.ndim == 2:
            dc = dc[0]
        if np.median(dc) < 100.0:
            dc = dc * 6371229.0
        pres = np.asarray(frame(ref, "pressure", ir), dtype=np.float64)
        rho = np.asarray(frame(ref, "rho", ir), dtype=np.float64)
        cs = np.sqrt(1.4 * np.maximum(pres, 1.0) / np.maximum(rho, 1.0e-6))
        cs_max = np.max(cs, axis=1)
        cs_edge = np.zeros(dc.shape, dtype=np.float64)
        cs_edge[valid] = 0.5 * (cs_max[c1[valid]] + cs_max[c2[valid]])
        print("acoustic CFL = max(cs) * dts / dcEdge, gamma=1.4, reference frame")
        for label, mask in (
            ("fine edges", edge_fine),
            ("interface edges", edge_interface),
            ("coarse edges", edge_coarse),
        ):
            dcm = dc[mask]
            cse = cs_edge[mask]
            print(f"  {label:18s} n {int(mask.sum()):6d}  "
                  f"dc km p50 {np.percentile(dcm, 50)/1000:.1f}  "
                  f"p10 {np.percentile(dcm, 10)/1000:.1f}  "
                  f"min {np.min(dcm)/1000:.1f}")
            for dts, tag in ((30.0, "dts 30s (fine / matched)"), (60.0, "dts 60s (coarse N=6)")):
                cfl = cse * dts / dcm
                print(f"    {tag:28s} CFL p50 {np.percentile(cfl, 50):.3f}  "
                      f"p99 {np.percentile(cfl, 99):.3f}  max {np.max(cfl):.3f}")

    area = np.asarray(test.variables["areaCell"][:], dtype=np.float64)
    zgrid = np.asarray(test.variables["zgrid"][:], dtype=np.float64)
    if zgrid.ndim == 3:
        zgrid = zgrid[0]
    rho_t = np.asarray(test.variables["rho"][:], dtype=np.float64)
    mt = dry_mass(rho_t, area, zgrid)
    area_r = np.asarray(ref.variables["areaCell"][:], dtype=np.float64)
    zgrid_r = np.asarray(ref.variables["zgrid"][:], dtype=np.float64)
    if zgrid_r.ndim == 3:
        zgrid_r = zgrid_r[0]
    rho_r = np.asarray(ref.variables["rho"][:], dtype=np.float64)
    mr = dry_mass(rho_r, area_r, zgrid_r)
    print(f"mass test drift {(mt[it] - mt[0]) / mt[0]:.6e}")
    print(f"mass ref  drift {(mr[ir] - mr[0]) / mr[0]:.6e}")
    print(f"mass frame diff {(mt[it] - mr[ir]) / mr[ir]:.6e}")

    test.close()
    ref.close()
    grid.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
