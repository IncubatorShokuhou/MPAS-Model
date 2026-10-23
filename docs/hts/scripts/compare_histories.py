#!/usr/bin/env python3
"""Compare two MPAS-A history files on u, w, theta, rho and dry mass.

Relative L2 is sqrt(mean((a-b)^2)) / sqrt(mean(b^2)), with b the reference.
Dry mass is sum(rho * areaCell * dz) with dz from zgrid.
"""

import argparse
import sys

import numpy as np
from netCDF4 import Dataset


def dry_mass(rho, area, zgrid):
    if rho.ndim == 2:
        rho = rho[None, ...]
    dz = zgrid[:, 1:] - zgrid[:, :-1]
    return np.sum(rho * area[None, :, None] * dz[None, :, :], axis=(1, 2))


def rel_l2(test, ref):
    diff = np.asarray(test, dtype=np.float64) - np.asarray(ref, dtype=np.float64)
    denom = np.sqrt(np.mean(np.asarray(ref, dtype=np.float64) ** 2))
    num = np.sqrt(np.mean(diff ** 2))
    if denom == 0.0:
        return num
    return num / denom


def linf(test, ref):
    diff = np.abs(np.asarray(test, dtype=np.float64) - np.asarray(ref, dtype=np.float64))
    return float(np.max(diff))


def frame(nc, name, index):
    var = nc.variables[name]
    if "Time" in var.dimensions:
        return var[index, ...]
    return var[:]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("test")
    p.add_argument("reference")
    p.add_argument("--frame", type=int, default=-1, help="time index in both files (default last)")
    p.add_argument("--ref-frame", type=int, default=None, help="reference time index if it differs")
    args = p.parse_args()

    test = Dataset(args.test)
    ref = Dataset(args.reference)
    nt = len(test.dimensions["Time"])
    nr = len(ref.dimensions["Time"])
    it = args.frame if args.frame >= 0 else nt - 1
    ir = args.ref_frame if args.ref_frame is not None else (args.frame if args.frame >= 0 else nr - 1)
    print(f"frames test={it}/{nt} reference={ir}/{nr}")

    for name in ("u", "w", "theta", "rho"):
        a = frame(test, name, it)
        b = frame(ref, name, ir)
        print(f"{name:7s}  relL2 {rel_l2(a, b):.6e}  Linf {linf(a, b):.6e}")

    area = np.asarray(test.variables["areaCell"][:], dtype=np.float64)
    zgrid = np.asarray(test.variables["zgrid"][:], dtype=np.float64)
    if zgrid.ndim == 3:
        zgrid = zgrid[0]
    rho_t = np.asarray(test.variables["rho"][:], dtype=np.float64)
    rho_r = np.asarray(ref.variables["rho"][:], dtype=np.float64)
    mt = dry_mass(rho_t, area, zgrid)
    # Reference area/zgrid should match the same mesh. Use the reference file's own geometry.
    area_r = np.asarray(ref.variables["areaCell"][:], dtype=np.float64)
    zgrid_r = np.asarray(ref.variables["zgrid"][:], dtype=np.float64)
    if zgrid_r.ndim == 3:
        zgrid_r = zgrid_r[0]
    mr = dry_mass(rho_r, area_r, zgrid_r)
    print(f"mass test  t0 {mt[0]:.10e}  frame {mt[it]:.10e}  drift {(mt[it] - mt[0]) / mt[0]:.6e}")
    print(f"mass ref   t0 {mr[0]:.10e}  frame {mr[ir]:.10e}  drift {(mr[ir] - mr[0]) / mr[0]:.6e}")
    print(f"mass frame diff {(mt[it] - mr[ir]) / mr[ir]:.6e}")

    if "pressure" in test.variables and "pressure" in ref.variables:
        pt = np.asarray(frame(test, "pressure", it))
        pr = np.asarray(frame(ref, "pressure", ir))
        # Lowest model level is index 0 (near-surface).
        d = pt[:, 0] - pr[:, 0]
        print(f"psurf  relL2 {rel_l2(pt[:, 0], pr[:, 0]):.6e}  rms {np.sqrt(np.mean(d**2)):.6e} Pa  max {np.max(np.abs(d)):.6e} Pa")

    test.close()
    ref.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
