#!/usr/bin/env python3
# sqcalc - structure factors from LAMMPS dump trajectories
# Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Plot the text tables written by sqcalc.

Reads the ``# q S(q) ...`` shell table, the ``# r g(r) ...`` Debye table or the
``# qx qy qz S(q)`` reciprocal grid table, takes the column names from the
leading comment line and draws the total curve plus the partials.  Several
tables can be passed at once to overlay their totals; a grid table is drawn as
a scatter against ``|q|``.

Usage:
    plot_sq.py S_q.dat                        # S(q) and partials -> S_q.png
    plot_sq.py -o compare.png a.dat b.dat     # overlay two totals
    plot_sq.py --show --partials rdf.dat      # interactive g(r) window
    plot_sq.py --refline 0 S_q.dat            # mark zero instead of the 1 line

Matplotlib is the only dependency.  The ``--output`` suffix picks the format
(``.png``, ``.pdf``, ``.svg``, ...); ``.pdf``/``.svg`` are the vector choices.
A dashed reference line is drawn at y = 1 by default, the value S(q) and g(r)
tend to at large q / r (and the Faber-Ziman partials tend to); ``--refline``
takes another value, or ``none`` to leave the figure clean.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


def read_table(path: Path) -> tuple[list[str] | None, list[list[float]]]:
    """Return the column names from the header line and the numeric rows."""
    labels: list[str] | None = None
    rows: list[list[float]] = []
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                # The first comment line is the column header, e.g.
                # "# q S(q) S(Si-Si) S(O-O)".
                if labels is None:
                    tokens = stripped.lstrip("#").split()
                    if len(tokens) >= 2:
                        labels = tokens
                continue
            try:
                rows.append([float(value) for value in stripped.split()])
            except ValueError:
                continue
    if not rows:
        raise SystemExit(f"{path}: no numeric rows found")
    return labels, rows


def kind_of(labels: list[str] | None, ncols: int) -> str:
    """Classify a table as a grid, a g(r) table or an S(q) table."""
    if labels:
        first = labels[0].lower()
        if first.startswith("q") and len(first) > 1:
            return "grid"
        if first.startswith("r"):
            return "rdf"
        if first.startswith("q"):
            return "sq"
    # No usable header: a grid has four columns, everything else is a curve.
    return "grid" if ncols >= 4 else "sq"


def column(labels: list[str] | None, index: int, default: str) -> str:
    if labels and index < len(labels):
        return labels[index]
    return default


def parse_refline(text: str) -> float | None:
    """Reference line value, or None when the line is switched off."""
    if text.strip().lower() in ("", "none", "off", "no"):
        return None
    try:
        return float(text)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"expects a number or 'none', not {text!r}")


def plot_curves(ax, tables, partials: bool) -> None:
    """Draw the S(q) or g(r) total of each table, plus the partials."""
    for path, labels, rows in tables:
        xs = [row[0] for row in rows]
        total = column(labels, 1, "y")
        name = total if len(tables) == 1 else f"{path.stem} ({total})"
        ax.plot(xs, [row[1] for row in rows], label=name)
        if not partials:
            continue
        for index in range(2, len(rows[0])):
            label = column(labels, index, f"column {index + 1}")
            ax.plot(xs, [row[index] for row in rows], linestyle=":", label=label)


def plot_grid(ax, tables) -> None:
    """Draw a reciprocal grid table as S(q) against |q|."""
    for path, labels, rows in tables:
        qs = [
            (row[0] ** 2 + row[1] ** 2 + row[2] ** 2) ** 0.5
            for row in rows
        ]
        total = column(labels, 3, "S(q)")
        name = total if len(tables) == 1 else f"{path.stem} ({total})"
        ax.scatter(qs, [row[3] for row in rows], s=4, label=name)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot the text output of sqcalc (S(q), g(r) or a q-grid)."
    )
    parser.add_argument("tables", nargs="+", metavar="TABLE",
                        help="text table written by sqcalc")
    parser.add_argument("-o", "--output", metavar="FIGURE",
                        help="figure file; defaults to the first table with a "
                             ".png suffix")
    parser.add_argument("--show", action="store_true",
                        help="open an interactive window instead of only saving")
    parser.add_argument("--partials", action="store_true",
                        help="also draw the partial columns of every table")
    parser.add_argument("--no-partials", action="store_true",
                        help="draw only the totals")
    parser.add_argument("--title", metavar="TITLE", help="figure title")
    parser.add_argument("--refline", metavar="Y", type=parse_refline,
                        default=1.0,
                        help="horizontal reference line, drawn at the 1 that "
                             "S(q), g(r) and the Faber-Ziman partials tend to "
                             "(default: 1.0; use 0 or 'none')")
    args = parser.parse_args()

    if args.partials and args.no_partials:
        parser.error("--partials and --no-partials are mutually exclusive")
    refline = args.refline

    # Read and validate the tables before pulling in matplotlib, so a typo or a
    # mixed set of tables fails without the import cost.
    tables = []
    kinds = set()
    for name in args.tables:
        path = Path(name)
        if not path.is_file():
            print(f"error: {path} not found", file=sys.stderr)
            return 2
        labels, rows = read_table(path)
        kinds.add(kind_of(labels, len(rows[0])))
        tables.append((path, labels, rows))
    if len(kinds) > 1:
        print("error: cannot mix grid tables with S(q)/g(r) tables",
              file=sys.stderr)
        return 2
    kind = kinds.pop()

    if not args.show:
        os.environ.setdefault("MPLBACKEND", "Agg")
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("error: plot_sq.py needs matplotlib (pip install matplotlib)",
              file=sys.stderr)
        return 2

    # Partials are only useful next to a single total.
    partials = args.partials or (len(tables) == 1 and not args.no_partials)

    figure, ax = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
    if kind == "grid":
        plot_grid(ax, tables)
        ax.set_xlabel("|q| [1/A]")
        ax.set_ylabel("S(q)")
    else:
        plot_curves(ax, tables, partials)
        ax.set_xlabel(column(tables[0][1], 0, "q [1/A]"))
        ax.set_ylabel(column(tables[0][1], 1, "S(q)"))
    if refline is not None:
        ax.axhline(refline, color="0.65", linewidth=0.8, linestyle="--",
                   zorder=0)
        ax.annotate(f"{refline:g}", xy=(1.0, refline),
                    xycoords=("axes fraction", "data"), xytext=(-3, 3),
                    textcoords="offset points", ha="right", va="bottom",
                    color="0.45", fontsize="small")
    ax.legend()
    ax.set_title(args.title or (tables[0][0].stem if len(tables) == 1
                                else "sqcalc"))

    output = Path(args.output) if args.output else \
        tables[0][0].with_suffix(".png")
    figure.savefig(output, dpi=150)
    print(f"wrote {output}")
    if args.show:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
