#!/usr/bin/env python3
"""Convert LAMMPS dump frames to AtomEye .cfg files (one file per frame).

Only a helper for cross-checking sqcalc against the external `debyer` program
(`test/compare_debyer.sh`); nothing in the sqcalc build or test suite depends on
debyer.  .cfg stores fractional coordinates and the cell in H0(i,j) lines, which
is what debyer's reader expects.
"""

import argparse

import numpy as np

import ref_sq


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--mapping", default="1:Si,2:O")
    parser.add_argument("--prefix", required=True, help="output prefix, frames get _000, _001, ...")
    args = parser.parse_args()

    mapping = ref_sq.parse_mapping(args.mapping)
    written = []
    for frame, (types, pos, cell, origin) in enumerate(ref_sq.read_dump(args.input)):
        frac = np.linalg.solve(cell.T, (pos - origin).T).T
        frac -= np.floor(frac)
        name = "%s_%03d.cfg" % (args.prefix, frame)
        with open(name, "w") as handle:
            handle.write("Number of particles = %d\n" % len(pos))
            handle.write("A = 1.0 Angstrom\n")
            for i in range(3):
                for j in range(3):
                    handle.write("H0(%d,%d) = %.12f A\n" % (i + 1, j + 1, cell[i][j]))
            handle.write("# sqcalc frame %d of %s\n" % (frame, args.input))
            for t, xyz in zip(types, frac):
                symbol = mapping.get(int(t), "X")
                handle.write("28.086 %s %.10f %.10f %.10f\n"
                             % (symbol, xyz[0], xyz[1], xyz[2]))
        written.append(name)
    print(" ".join(written))


if __name__ == "__main__":
    main()
