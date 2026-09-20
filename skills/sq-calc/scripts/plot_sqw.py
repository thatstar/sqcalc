#!/usr/bin/env python3
# sqcalc - structure factors from LAMMPS dump trajectories
# Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
#
# SPDX-License-Identifier: GPL-3.0-or-later

r"""Plot the dynamic tables written by ``sqcalc --dyn``.

Reads the ``# qx qy qz omega S(q,w) ...`` spectra table (or the ``tau F(q,t)``
table written by ``--fqt``, or the ``tau S4(q,t)`` table written by ``--s4``)
and draws the map of the quantity over $(|q|, \omega)$ next to a few selected
spectra, or single curves with ``--mode spectra``.

Usage:
    plot_sqw.py S_qw.dat                        # map and some spectra -> .png
    plot_sqw.py --mode spectra --q 1,4 S_qw.dat # only the spectra at |q| = 1, 4
    plot_sqw.py -o fig.pdf F_qt.dat             # vector output

Matplotlib is the only dependency.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


def read_table(path: Path) -> tuple[list[str], np.ndarray]:
    """Return the column names of the last header line and the numeric rows."""
    labels: list[str] = []
    rows: list[list[float]] = []
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                names = stripped.lstrip("#").split()
                # keep the last line that looks like the column list
                if len(names) >= 5:
                    labels = names
                continue
            rows.append([float(value) for value in stripped.split()])
    if not rows:
        raise SystemExit("no data rows in %s" % path)
    return labels, np.array(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("table", type=Path, metavar="TABLE")
    parser.add_argument("-o", "--output", metavar="FIGURE",
                        help="output figure (default: the table name .png)")
    parser.add_argument("--show", action="store_true", help="open a window")
    parser.add_argument("--mode", choices=["both", "map", "spectra"], default="both")
    parser.add_argument("--q", metavar="LIST",
                        help="comma separated |q| values to draw as spectra")
    parser.add_argument("--column", type=int, default=4,
                        help="column of the quantity to plot (default 4)")
    parser.add_argument("--vmax", type=float, default=None,
                        help="upper colour limit of the map")
    args = parser.parse_args()

    import matplotlib

    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels, data = read_table(args.table)
    qmag = np.hypot(np.hypot(data[:, 0], data[:, 1]), data[:, 2])
    axis = data[:, 3]
    values = data[:, args.column]
    q_values = np.unique(np.round(qmag, 8))
    axis_values = np.unique(np.round(axis, 10))
    grid = values.reshape(len(q_values), len(axis_values))
    axis_name = labels[3] if len(labels) > 4 else "axis"
    quantity = labels[args.column] if len(labels) > args.column else "S(q,w)"

    wanted = q_values
    if args.q:
        wanted = np.array([float(v) for v in args.q.split(",")])

    panels = 1 if args.mode in ("map", "spectra") else 2
    fig, axes = plt.subplots(1, panels, figsize=(5.0*panels, 4.0), dpi=140,
                             squeeze=False)
    column = 0
    if args.mode in ("map", "both"):
        ax = axes[0][column]
        mesh = ax.pcolormesh(axis_values, q_values, grid, shading="auto",
                             vmax=args.vmax if args.vmax else np.max(grid))
        fig.colorbar(mesh, ax=ax, label=quantity)
        ax.set_xlabel(axis_name)
        ax.set_ylabel(r"$|q|$ [1/A]")
        ax.set_title(quantity)
        column += 1
    if args.mode in ("spectra", "both"):
        ax = axes[0][column]
        for q in wanted:
            index = int(np.argmin(np.abs(q_values - q)))
            ax.plot(axis_values, grid[index], label=r"$|q| = %.2f$" % q_values[index])
        ax.set_xlabel(axis_name)
        ax.set_ylabel(quantity)
        ax.legend(fontsize="small")
    fig.tight_layout()

    output = args.output or str(args.table.with_suffix(".png"))
    fig.savefig(output)
    print("wrote %s" % output)
    if args.show:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
