#!/usr/bin/env python3
"""Level-balanced METIS partition for hierarchical time stepping.

MPAS graph.info is an unweighted cell graph. A default partition can put every
fine cell on one rank and every coarse cell on another. During a fine step the
coarse-only rank finishes immediately and the rest wait inside the halo
exchange. This script writes a multi-constraint graph whose vertex weights
balance each HTS level across ranks, then runs gpmetis.

Level of a cell matches the model: L = nint(log2(max(meshDensity^(-1/4), 1))),
clamped to --max-level. Edge level is not needed here; the weight is per cell.

Usage:
  docs/hts/scripts/hts_partition.py GRID.nc GRAPH.info NPARTS --max-level 1

Writes GRAPH.info.htsw.part.NPARTS (one 0-based rank id per cell, the file
MPAS reads through config_block_decomp_file_prefix). Also writes the weighted
graph next to it as GRAPH.info.htsw so the gpmetis command can be repeated:

  gpmetis -ptype=kway GRAPH.info.htsw NPARTS

The weighted graph header is "n m 10 ncon". Each line is ncon positive integer
weights followed by the 1-based neighbor ids from graph.info. The cell's own
level has weight 1000 and every other level has weight 1. METIS rejects a zero
weight. ncon is the number of levels present (max used level + 1).

`--physics-weight N` adds N to every constraint and adds it again, with the
dynamics weight, on the cell's own level. N = 0 leaves the graph above unchanged.
Use it when physics_driver is a large share of the fine step: that driver still
visits every column on every fine tick, so a uniform physics cost belongs on
every level constraint.
"""

import argparse
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def cell_levels(mesh_density, max_level):
    md = np.asarray(mesh_density, dtype=np.float64)
    spacing = np.maximum(md ** (-0.25), 1.0)
    level = np.rint(np.log(spacing) / np.log(2.0)).astype(np.int32)
    level = np.clip(level, 0, max_level)
    # A uniform mesh (or non-positive density) is one level, matching the model.
    if md.size == 0 or np.nanmax(md) <= np.nanmin(md) * (1.0 + 1.0e-5) or np.nanmin(md) <= 0.0:
        level[:] = 0
    return level


def read_graph(path):
    lines = Path(path).read_text().splitlines()
    header = lines[0].split()
    n = int(header[0])
    m = int(header[1])
    body = lines[1:]
    if len(body) != n:
        raise SystemExit(f"{path}: expected {n} adjacency lines, found {len(body)}")
    return n, m, body


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("grid")
    p.add_argument("graph")
    p.add_argument("nparts", type=int)
    p.add_argument("--max-level", type=int, default=1)
    p.add_argument("--weight", type=int, default=1000, help="weight on the cell's own level")
    p.add_argument("--physics-weight", type=int, default=0,
                   help="per-cell physics cost added on every level constraint; "
                        "0 keeps the dynamics-only weights")
    p.add_argument("--gpmetis", default="gpmetis")
    args = p.parse_args()

    grid = Dataset(args.grid)
    md = np.asarray(grid.variables["meshDensity"][:], dtype=np.float64)
    grid.close()
    n, m, body = read_graph(args.graph)
    if md.shape[0] != n:
        raise SystemExit(f"meshDensity length {md.shape[0]} != graph cells {n}")

    level = cell_levels(md, args.max_level)
    nlev = int(level.max()) + 1
    counts = np.bincount(level, minlength=nlev)
    print(f"cells {n}  levels {nlev}  counts {counts.tolist()}")
    fine = counts[0] / n
    # One coarse interval: fine cells take 2^L large steps per coarse step of the top level.
    # With matched acoustic dts the acoustic work per cell is flat, so the large-step
    # ceiling is the useful one.
    top = nlev - 1
    tend_frac = 0.0
    for L in range(nlev):
        tend_frac += (counts[L] / n) / (2 ** L)
    print(f"fine fraction {fine:.4f}")
    print(f"large-step work / all-fine over one top-level step: {tend_frac:.4f}")
    print(f"large-step ceiling (all-fine / that work): {1.0 / tend_frac:.3f}x")

    graph_path = Path(args.graph)
    weighted = graph_path.with_name(graph_path.name + ".htsw")
    with weighted.open("w") as out:
        # fmt 10: vertex weights present, no edge weights. ncon follows.
        out.write(f"{n} {m} 10 {nlev}\n")
        if args.physics_weight < 0:
            raise SystemExit("--physics-weight must be >= 0")
        # physics_weight == 0 is the phase 2b graph: own level `weight`, every
        # other level 1. A positive value is the cost of physics_driver, which
        # still runs on every column on every fine tick. It is added to every
        # constraint, and the masked dynamics cost stays on the cell's own level.
        # METIS rejects a zero, so the no-physics other-level weight stays 1.
        other = args.physics_weight if args.physics_weight > 0 else 1
        for i, adj in enumerate(body):
            weights = [other] * nlev
            weights[int(level[i])] = args.weight + args.physics_weight
            out.write(" ".join(str(w) for w in weights))
            adj = adj.strip()
            if adj:
                out.write(" " + adj)
            out.write("\n")
    print(f"wrote {weighted}")

    cmd = [args.gpmetis, "-ptype=kway", str(weighted), str(args.nparts)]
    print(" ".join(cmd))
    subprocess.check_call(cmd)
    produced = Path(str(weighted) + f".part.{args.nparts}")
    if not produced.is_file():
        raise SystemExit(f"gpmetis did not write {produced}")

    part = np.loadtxt(produced, dtype=np.int32)
    if part.shape[0] != n:
        raise SystemExit(f"partition length {part.shape[0]} != {n}")
    print(f"wrote {produced}")
    for rank in range(args.nparts):
        sel = part == rank
        lc = np.bincount(level[sel], minlength=nlev)
        print(f"  rank {rank}: cells {int(sel.sum())}  levels {lc.tolist()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
