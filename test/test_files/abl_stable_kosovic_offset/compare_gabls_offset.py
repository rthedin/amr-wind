#!/usr/bin/env python3
"""Compare the flat and terrain-offset GABLS1 runs with the Kosovic model.

The flat run (abl_stable_kosovic) uses the ABLWallFunction on the domain
floor. The offset run (abl_stable_kosovic_offset) raises the ground by
``--offset`` metres with TerrainDrag, and DragForcing/DragTempForcing set the
surface layer from ABL.surface_stability_type = heating_rate. Both cool the
surface at 0.25 K/hr from 265 K.

The script reads the ABLStats ASCII plane averages of both runs, averages them
over a time window, shifts the offset heights down by ``--offset`` and plots
wind speed, wind direction, potential temperature and the resolved vertical
heat flux. It prints the largest differences between the two runs.

Usage:
    python3 compare_gabls_offset.py FLAT_DIR OFFSET_DIR \
        [--offset 50] [--tstart 28800] [--tend 32400] [--out gabls_offset.png]

FLAT_DIR and OFFSET_DIR are the run directories that contain post_processing/.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_plane_average(path):
    """Return (time, z, values[t, z, comp]) from an ABLStats ASCII file."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = [p for p in line.replace(",", " ").split()]
            if len(parts) < 4:
                continue
            rows.append([float(p) for p in parts])
    data = np.array(rows)
    times = np.unique(data[:, 1])
    z = np.unique(data[:, 2])
    ncomp = data.shape[1] - 3
    vals = np.full((times.size, z.size, ncomp), np.nan)
    ti = np.searchsorted(times, data[:, 1])
    zi = np.searchsorted(z, data[:, 2])
    vals[ti, zi, :] = data[:, 3:]
    return times, z, vals


def time_mean(path, tstart, tend):
    times, z, vals = read_plane_average(path)
    sel = (times >= tstart) & (times <= tend)
    if not np.any(sel):
        raise SystemExit(
            f"{path}: no output between {tstart} and {tend} s "
            f"(last time {times[-1]:.1f} s)"
        )
    return z, np.nanmean(vals[sel], axis=0), times[sel]


def profiles(run_dir, tstart, tend, offset):
    pp = Path(run_dir) / "post_processing"
    z, vel, used = time_mean(pp / "plane_average_velocity.txt", tstart, tend)
    _, temp, _ = time_mean(pp / "plane_average_temperature.txt", tstart, tend)
    _, tu, _ = time_mean(pp / "second_moment_temperature_velocity.txt", tstart, tend)
    zag = z - offset
    keep = zag > 0.0
    speed = np.hypot(vel[:, 0], vel[:, 1])
    direction = np.degrees(np.arctan2(vel[:, 1], vel[:, 0]))
    return dict(
        z=zag[keep],
        speed=speed[keep],
        direction=direction[keep],
        theta=temp[keep, 0],
        wtheta=tu[keep, 2],
        n=used.size,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("flat_dir")
    ap.add_argument("offset_dir")
    ap.add_argument("--offset", type=float, default=50.0)
    ap.add_argument("--tstart", type=float, default=28800.0)
    ap.add_argument("--tend", type=float, default=32400.0)
    ap.add_argument("--out", default="gabls_offset.png")
    args = ap.parse_args()

    flat = profiles(args.flat_dir, args.tstart, args.tend, 0.0)
    off = profiles(args.offset_dir, args.tstart, args.tend, args.offset)

    print(f"averaged {flat['n']} flat and {off['n']} offset outputs "
          f"between {args.tstart:.0f} and {args.tend:.0f} s")
    print(f"{'quantity':>22s} {'max |diff|':>12s} {'at z (m)':>9s}")
    for key, label in (
        ("speed", "wind speed (m/s)"),
        ("direction", "direction (deg)"),
        ("theta", "theta (K)"),
        ("wtheta", "<w'theta'> (K m/s)"),
    ):
        other = np.interp(flat["z"], off["z"], off[key])
        diff = np.abs(other - flat[key])
        i = int(np.argmax(diff))
        print(f"{label:>22s} {diff[i]:12.4g} {flat['z'][i]:9.1f}")

    fig, axes = plt.subplots(1, 4, figsize=(13, 4.5), sharey=True)
    panels = (
        ("speed", "Wind speed (m/s)"),
        ("direction", "Wind direction (deg)"),
        ("theta", r"$\theta$ (K)"),
        ("wtheta", r"Resolved $\langle w'\theta'\rangle$ (K m/s)"),
    )
    for ax, (key, xlabel) in zip(axes, panels):
        ax.plot(flat[key], flat["z"], "k-", lw=1.8, label="Flat, ABLWallFunction")
        ax.plot(off[key], off["z"], "C3--", lw=1.8,
                label=f"Offset {args.offset:g} m, heating_rate")
        ax.set_xlabel(xlabel)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("Height above ground (m)")
    axes[0].set_ylim(0, 300)
    axes[0].legend(loc="upper left", fontsize=8)
    fig.suptitle(
        "GABLS1, Kosovic, surface cooling 0.25 K/hr; "
        f"mean over {args.tstart / 3600:g}-{args.tend / 3600:g} h"
    )
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
