#!/usr/bin/env python3
"""Summarize a Jablonowski-Williamson baseline run directory.

Reads init.nc, history.nc, time.model.txt, and log.atmosphere.0000.out.
Prints a short text report: wall time, per-step time, dry-mass drift,
surface-pressure departure, and velocity extrema.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

G = 9.80616


def _arr(ds, name):
    if name not in ds.variables:
        return None
    return np.array(ds.variables[name][:], dtype=np.float64)


def dry_mass(rho, area, zgrid):
    """rho is (nCells, nVert) or (time, nCells, nVert). zgrid is (nCells, nVert+1)."""
    if rho.ndim == 2:
        rho = rho[None, ...]
    dz = zgrid[:, 1:] - zgrid[:, :-1]
    return np.sum(rho * area[None, :, None] * dz[None, :, :], axis=(1, 2))


def ps_stats(ps, ps0):
    # ps, ps0 shaped (nCells,)
    d = ps - ps0
    return float(np.sqrt(np.mean(d ** 2))), float(np.max(np.abs(d))), float(ps.min()), float(ps.max())


def parse_log(path: Path):
    step_times = []
    u_minmax = []
    integ = None
    calls = None
    for line in path.read_text(errors="replace").splitlines():
        if "Timing for integration step:" in line:
            step_times.append(float(line.split()[-2]))
        if "global min, max u" in line:
            parts = line.split()
            u_minmax.append((float(parts[-2]), float(parts[-1])))
        if "time integration" in line and "total time" not in line:
            m = re.search(r"time integration\s+([0-9.]+)\s+([0-9]+)", line)
            if m:
                integ = float(m.group(1))
                calls = int(m.group(2))
    return step_times, u_minmax, integ, calls


def parse_time(path: Path):
    out = {}
    if not path.exists():
        return out
    for line in path.read_text().splitlines():
        bits = line.split()
        if len(bits) >= 2:
            out[bits[0]] = bits[1]
    return out


def main(run_dir: str):
    d = Path(run_dir)
    init = Dataset(d / "init.nc")
    hist = Dataset(d / "history.nc")

    area = _arr(init, "areaCell")
    if area is None:
        area = _arr(hist, "areaCell")
        if area.ndim == 2:
            area = area[0]

    def pack(ds, it=None):
        rho = _arr(ds, "rho")
        zg = _arr(ds, "zgrid")
        ps = _arr(ds, "surface_pressure")
        if rho.ndim == 3:
            rho = rho[0 if it is None else it]
        if ps.ndim == 2:
            ps = ps[0 if it is None else it]
        if zg.ndim == 3:
            zg = zg[0 if it is None else it]
        return rho, zg, ps

    rho0, zg0, ps0 = pack(init)
    mass0 = dry_mass(rho0, area, zg0)[0]
    ps_mass0 = float(np.sum(ps0 * area) / G)

    xtime = hist.variables["xtime"][:] if "xtime" in hist.variables else None
    ntime = hist.variables["rho"].shape[0]
    print(f"run {d}")
    print(f"  init dry mass (rho dV) {mass0:.8e} kg")
    print(f"  init ps*area/g         {ps_mass0:.8e} kg")
    print(f"  history times {ntime}")
    for it in range(ntime):
        rho, zg, ps = pack(hist, it)
        mass = dry_mass(rho, area, zg)[0]
        rms, linf, pmin, pmax = ps_stats(ps, ps0)
        label = ""
        if xtime is not None:
            raw = xtime[it]
            if isinstance(raw, np.ndarray):
                label = "".join(c.decode() if isinstance(c, (bytes, np.bytes_)) else str(c) for c in raw).strip()
            else:
                label = str(raw).strip()
        print(
            f"  t={label or it}  mass {mass:.8e}  rel {(mass - mass0) / mass0:.3e}"
            f"  ps_rms_anom {rms:.4f} Pa  ps_absmax_anom {linf:.4f} Pa"
            f"  ps[{pmin:.2f},{pmax:.2f}]"
        )
    init.close()
    hist.close()

    steps, umm, integ, calls = parse_log(d / "log.atmosphere.0000.out")
    if steps:
        print(
            f"  integration steps {len(steps)}  mean {np.mean(steps):.4f} s"
            f"  min {np.min(steps):.4f}  max {np.max(steps):.4f}"
        )
    if integ is not None:
        print(f"  timer 'time integration' total {integ:.4f} s over {calls} calls")
    if umm:
        print(f"  first logged u min/max {umm[0][0]:.4f} {umm[0][1]:.4f}")
        print(f"  last  logged u min/max {umm[-1][0]:.4f} {umm[-1][1]:.4f}")
    print("  " + " ".join(f"{k}={v}" for k, v in parse_time(d / "time.model.txt").items()))
    print("  " + " ".join(f"{k}={v}" for k, v in parse_time(d / "time.init.txt").items()))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: summarize_baseline.py RUN_DIR [RUN_DIR...]")
    for rd in sys.argv[1:]:
        main(rd)
