#!/usr/bin/env python3
# sqcalc - structure factors from LAMMPS dump trajectories
# Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Pick --qmin, --qmax and --nq for a LAMMPS dump from its box.

The reciprocal methods sample q on the lattice of the box, so only the values
|q| = |h b1 + k b2 + l b3| exist.  Those values are not uniformly spaced: a
shell narrower than the gaps between them stays empty, and sqcalc writes a
flat 0 there, which is not a measurement.  A small box makes this worse,
because the lattice spacing 2*pi/L is large.

This script reads the box from the first frame of a dump, mirrors the mode
grid sqcalc builds, and reports

* the smallest |q| the box can represent, to use as --qmin (the default 0
  leaves a dead zone below it),
* the shell width dq, i.e. the --nq to use, that provably leaves no empty
  shell, next to what the usual --nq 500 would do,
* the grid size for --qmax, against the 4e8 mode limit, and
* the resulting option string, with --print-options for scripting.

Usage:
    choose_q.py traj.dump                    # report for qmax = 20 (default)
    choose_q.py traj.dump --qmax 15
    choose_q.py traj.dump --qmax 12 --print-options

Only numpy is needed.  Pair this with a frame count that gives enough
independent configurations: the error of a shell average falls as
1/sqrt(frames * modes in the shell), and the modes column printed here shows
where the statistics are thin.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

MAX_GRID = 4.0e8          # sqcalc refuses reciprocal grids above this size
MAX_NQ = 4000             # bound for the "finest shell width" search
ENUM_CAP = 4_000_000      # largest mode grid this script enumerates itself


def read_cell(path: Path) -> np.ndarray:
    """Cell vectors a1, a2, a3 of the first frame of a LAMMPS dump."""
    with path.open() as handle:
        for line in handle:
            if not line.startswith("ITEM: BOX BOUNDS"):
                continue
            triclinic = "xy" in line.split()
            rows = []
            for _ in range(3):
                rows.append([float(value) for value in handle.readline().split()])
            if triclinic:
                (xlo, xhi, xy), (ylo, yhi, xz), (zlo, zhi, yz) = rows
            else:
                (xlo, xhi), (ylo, yhi), (zlo, zhi) = rows
                xy = xz = yz = 0.0
            # LAMMPS box convention.
            return np.array([[xhi - xlo, 0.0, 0.0],
                             [xy, yhi - ylo, 0.0],
                             [xz, yz, zhi - zlo]])
    raise SystemExit(f"{path}: no 'ITEM: BOX BOUNDS' line found")


def modes_per_axis(qmax: float, cell: np.ndarray) -> np.ndarray:
    """Grid sizes sqcalc uses, one per cell vector."""
    lengths = np.linalg.norm(cell, axis=1)
    return 2*np.ceil(qmax*lengths/(2.0*np.pi)).astype(int) + 1


def lattice_qvec(cell: np.ndarray, modes: np.ndarray) -> np.ndarray:
    """Every |q| on the sqcalc mode grid, without the q = 0 mode."""
    volume = float(np.dot(cell[0], np.cross(cell[1], cell[2])))
    b = np.array([2.0*np.pi*np.cross(cell[1], cell[2])/volume,
                  2.0*np.pi*np.cross(cell[2], cell[0])/volume,
                  2.0*np.pi*np.cross(cell[0], cell[1])/volume])
    axes = [np.arange(1, int(n) + 1) - (int(n) + 1)//2 for n in modes]
    h1, h2, h3 = np.meshgrid(*axes, indexing="ij")
    h = np.stack([h1.ravel(), h2.ravel(), h3.ravel()], axis=1)
    h = h[~np.all(h == 0, axis=1)]
    return np.linalg.norm(h @ b, axis=1)


def shell_counts(qlen: np.ndarray, qmin: float, qmax: float, nq: int) -> np.ndarray:
    """Occupancy of each q shell, exactly as sqcalc bins it."""
    dq = (qmax - qmin)/nq
    index = np.clip(((qlen - qmin)/dq).astype(int), 0, nq - 1)
    return np.bincount(index, minlength=nq)


def enumeration_grid(cell: np.ndarray, qmax: float) -> tuple[float, np.ndarray]:
    """qmax and grid sizes for the enumeration, bounded by ENUM_CAP."""
    scan = qmax
    modes = modes_per_axis(scan, cell)
    while int(np.prod(modes)) > ENUM_CAP and scan > 0.5:
        scan *= 0.9
        modes = modes_per_axis(scan, cell)
    return scan, modes


def finest_nq(qlen: np.ndarray, qmin: float, qmax: float, upper: float,
              verify: bool) -> int:
    """Largest nq (finest dq) that leaves no empty shell.

    A shell that fits entirely between two neighbouring q values is empty, so
    requiring dq to cover the widest step between neighbours (with qmin and
    qmax as the outer ends) rules empty shells out.  The alignment of the
    shell edges can sometimes allow a finer dq, but that gain is small
    compared with the safety of this criterion; the occupancy is verified on
    the actual grid before returning.
    """
    inside = qlen[(qlen >= qmin) & (qlen <= upper)]
    edges = np.concatenate([[qmin], np.sort(inside), [upper]])
    widest = float(np.max(np.diff(edges)))
    nq = max(1, min(MAX_NQ, int((qmax - qmin)/widest)))
    if verify:
        while nq > 1 and not np.all(shell_counts(qlen, qmin, qmax, nq) > 0):
            nq -= 1
    return nq


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Choose --qmin/--qmax/--nq for a LAMMPS dump from its box."
    )
    parser.add_argument("dump", metavar="DUMP", help="LAMMPS dump file")
    parser.add_argument("--qmax", type=float, default=20.0,
                        help="largest |q| wanted [1/A] (default: 20, sqcalc's "
                             "own default)")
    parser.add_argument("--qmin", type=float, default=None,
                        help="smallest |q|; defaults to the first one the box "
                             "can represent")
    parser.add_argument("--print-options", action="store_true",
                        help="print only the option string")
    args = parser.parse_args()

    path = Path(args.dump)
    if not path.is_file():
        print(f"error: {path} not found", file=sys.stderr)
        return 2

    cell = read_cell(path)
    modes = modes_per_axis(args.qmax, cell)
    requested = int(np.prod(modes))
    if requested > MAX_GRID:
        # Back off qmax until the grid fits sqcalc's own limit.
        qmax = args.qmax
        while np.prod(modes_per_axis(qmax, cell)) > MAX_GRID and qmax > 0.1:
            qmax *= 0.98
        qmax = float(np.floor(qmax*100.0)/100.0)
        print(f"warning: qmax {args.qmax:g} needs a grid of {requested:.3g} modes "
              f"(limit {MAX_GRID:.0e}); using qmax {qmax:g} instead")
        args.qmax = qmax

    scan_qmax, modes = enumeration_grid(cell, args.qmax)
    truncated = scan_qmax < args.qmax - 1.0e-9
    qlen = lattice_qvec(cell, modes)
    inside = qlen[qlen <= args.qmax]
    if inside.size == 0:
        print(f"error: qmax {args.qmax:g} is below the first accessible |q| "
              f"{qlen.min():.4f} 1/A", file=sys.stderr)
        return 2

    first_q = float(qlen.min())
    if args.qmin is None:
        qmin = first_q
    else:
        qmin = args.qmin
        if qmin < first_q:
            print(f"warning: qmin {qmin:g} is below the first accessible |q| "
                  f"{first_q:.4f}; those shells can never be filled")
            qmin = first_q
    qlen = qlen[(qlen >= qmin) & (qlen <= args.qmax)]

    upper = scan_qmax if truncated else args.qmax
    nq = finest_nq(qlen, qmin, args.qmax, upper, verify=not truncated)
    dq = (args.qmax - qmin)/nq
    # Occupancy of the shells the enumeration actually covers, using the same
    # shell width the recommendation would produce.
    scanned = qlen[(qlen >= qmin) & (qlen <= upper)]
    n_scanned = min(nq, int((upper - qmin)/dq) + 1)
    counts = shell_counts(scanned, qmin, qmin + n_scanned*dq, n_scanned)
    default_nq = 500
    default_counts = None
    if not truncated:
        default_counts = shell_counts(scanned, qmin, args.qmax, default_nq)
    options = (f"--qmin {qmin:.4f} --qmax {args.qmax:g} --nq {nq}")

    if args.print_options:
        print(options)
        return 0

    lengths = np.linalg.norm(cell, axis=1)
    print(f"box                : {lengths[0]:.4f} x {lengths[1]:.4f} x "
          f"{lengths[2]:.4f}")
    print(f"grid               : {modes[0]} x {modes[1]} x {modes[2]} = "
          f"{int(np.prod(modes))} modes (limit {MAX_GRID:.0e})")
    print(f"first accessible q : {first_q:.4f} 1/A  (2*pi/L for a cubic box)")
    if truncated:
        print(f"modes scanned      : {qlen.size} up to q = {scan_qmax:.2f} "
              f"(larger grids are not enumerated; the |q| steps only get "
              f"narrower above that)")
    else:
        print(f"modes up to qmax   : {qlen.size}")
    print()
    print(f"recommended        : {options}")
    print(f"  shell width dq   : {dq:.4f} 1/A")
    where = (f" over the first {n_scanned} shells (q <= {upper:.2f})"
             if truncated else "")
    print(f"  modes per shell  : min {counts.min()}, 1% "
          f"{int(np.percentile(counts, 1))}, median {int(np.median(counts))}{where}")
    if truncated:
        print(f"  shells           : {nq}, all {n_scanned} enumerated shells are "
              f"occupied")
    else:
        print(f"  shells           : {nq}, none empty")
    if default_counts is not None:
        print(f"  at --nq {default_nq}     : "
              f"{int((default_counts == 0).sum())} of {default_nq} shells empty "
              f"(sqcalc writes 0 there, and the low-q shells hold only a few "
              f"modes)")
    print()
    print("Notes")
    print("  * qmin removes the shells below the smallest reciprocal vector,")
    print("    which cannot be filled by any amount of sampling.")
    print("  * dq may be coarsened further for smoother curves, or refined")
    print("    together with a bigger box; dq ~ 0.3-0.4 * 2*pi/L is the")
    print("    practical limit before shells start to empty out.")
    print("  * a small box cannot resolve small q at all: use a larger box,")
    print("    or --method debye, whose q grid is not tied to the box.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
