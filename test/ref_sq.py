#!/usr/bin/env python3
# sqcalc - structure factors from LAMMPS dump trajectories
# Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Independent numpy reference for the sqcalc shell averaged S(q).

This is deliberately a second implementation: it reads the dump with its own
parser and evaluates rho(q) with a plain numpy exponential sum, so that it can
cross check the Fortran/NUFFT results.  Weighted runs are supported for the
few elements tabulated below (the numbers are the same IT92 / NN92 data that
src/sqc_element_data.f90 carries).
"""

import argparse

import numpy as np

NEUTRON_B = {"Si": 4.1491, "O": 5.803, "H": -3.7390, "C": 6.646, "Na": 3.63, "Cl": 9.5770}

XRAY = {
    "Si": ([6.2915, 3.0353, 1.9891, 1.541], [2.4386, 32.3337, 0.6785, 81.6937], 1.1407),
    "O": ([3.0485, 2.2868, 1.5463, 0.867], [13.2771, 5.7011, 0.3239, 32.9089], 0.2508),
    "H": ([0.493002, 0.322912, 0.140191, 0.040810],
          [10.5109, 26.1257, 3.14236, 57.7997], 0.003038),
}


def xray_f(symbol, q):
    a, b, c = XRAY[symbol]
    stol2 = (np.asarray(q, dtype=float) / (4.0 * np.pi)) ** 2
    return sum(ai * np.exp(-bi * stol2) for ai, bi in zip(a, b)) + c


def amplitude_of(type_id, q, mapping, weight):
    """Scattering amplitude of one LAMMPS type id at momentum transfer q."""
    symbol = mapping[type_id]
    if weight == "unit":
        return np.ones_like(q)
    if weight == "neutron":
        return np.full_like(q, NEUTRON_B[symbol])
    return xray_f(symbol, q)


def parse_mapping(spec):
    mapping = {}
    for item in spec.split(","):
        if not item.strip():
            continue
        key, symbol = item.split(":")
        mapping[int(key)] = symbol.strip().capitalize()
    return mapping


def read_dump(path):
    """Yield (types, positions, cell) for every frame of a LAMMPS dump."""
    with open(path) as handle:
        lines = handle.readlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("ITEM: TIMESTEP"):
            i += 2
        elif line.startswith("ITEM: NUMBER OF ATOMS"):
            natoms = int(lines[i + 1])
            i += 2
        elif line.startswith("ITEM: BOX BOUNDS"):
            tilt = "xy" in line.split()
            bounds = [list(map(float, lines[i + 1 + k].split())) for k in range(3)]
            i += 4
            xlo, xhi, xy = bounds[0] + [0.0] * (3 - len(bounds[0]))
            ylo, yhi, xz = bounds[1] + [0.0] * (3 - len(bounds[1]))
            zlo, zhi, yz = bounds[2] + [0.0] * (3 - len(bounds[2]))
            if tilt:
                xlo, xhi = (xlo - min(0.0, xy, xz, xy + xz), xhi - max(0.0, xy, xz, xy + xz))
                ylo, yhi = (ylo - min(0.0, yz), yhi - max(0.0, yz))
            cell = np.array([[xhi - xlo, 0.0, 0.0], [xy, yhi - ylo, 0.0], [xz, yz, zhi - zlo]])
            origin = np.array([xlo, ylo, zlo])
        elif line.startswith("ITEM: ATOMS"):
            columns = line.split()[2:]
            body = np.array([list(map(float, lines[i + 1 + k].split())) for k in range(natoms)])
            i += 1 + natoms
            types = body[:, columns.index("type")].astype(int)
            ix = [columns.index(name) for name in ("x", "y", "z")] if "x" in columns else \
                 [columns.index(name) for name in ("xu", "yu", "zu")]
            pos = body[:, ix]
            yield types, pos, cell, origin
        else:
            i += 1


def modes_for(cell, qmin, qmax):
    # cell has the lattice vectors in its rows, so r = A s with A = cell.T and
    # the reciprocal vectors are the columns of 2 pi (A^-1)^T.
    reciprocal = 2.0 * np.pi * np.linalg.inv(cell.T).T
    dims = [
        2 * int(np.ceil(qmax * np.linalg.norm(cell[i]) / (2.0 * np.pi))) + 1 for i in range(3)
    ]
    # Odd grid sizes with CMCL ordering: indices -(m-1)/2 ... (m-1)/2.
    ranges = [np.arange(m) - (m - 1) // 2 for m in dims]
    grids = np.meshgrid(*ranges, indexing="ij")
    hkl = np.stack([g.ravel() for g in grids], axis=1)
    hkl = hkl[np.any(hkl != 0, axis=1)]
    qvec = hkl @ reciprocal.T
    qlen = np.linalg.norm(qvec, axis=1)
    keep = (qlen <= qmax) & (qlen >= qmin)
    return hkl[keep], qvec[keep], qlen[keep]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--mapping", default="1:Si,2:O")
    parser.add_argument("--weight", default="unit", choices=["unit", "neutron", "xray"])
    parser.add_argument("--norm", default="mean", choices=["mean", "self", "n"])
    parser.add_argument("--qmin", type=float, default=0.0)
    parser.add_argument("--qmax", type=float, default=20.0)
    parser.add_argument("--nq", type=int, default=500)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    mapping = parse_mapping(args.mapping)
    dq = (args.qmax - args.qmin) / args.nq
    num = np.zeros(args.nq)
    den = np.zeros(args.nq)
    first = True

    for types, pos, cell, origin in read_dump(args.input):
        if first:
            hkl, qvec, qlen = modes_for(cell, args.qmin, args.qmax)
            shell = np.clip(((qlen - args.qmin) / dq).astype(int), 0, args.nq - 1)
            atom_symbols = np.array([mapping[t] for t in types])
            natoms = len(types)
            first = False
        frac = np.linalg.solve(cell.T, (pos - origin).T).T
        frac = frac - np.floor(frac)
        phases = 2.0 * np.pi * (hkl @ frac.T)
        if args.weight == "unit":
            amplitudes = np.ones((len(hkl), natoms))
        elif args.weight == "neutron":
            weights = np.array([NEUTRON_B[s] for s in atom_symbols])
            amplitudes = np.tile(weights, (len(hkl), 1))
        else:
            symbols = sorted(set(atom_symbols))
            amplitudes = np.zeros((len(hkl), natoms))
            for symbol in symbols:
                fq = xray_f(symbol, qlen)
                amplitudes[:, atom_symbols == symbol] = fq[:, None]
        rho = np.sum(amplitudes * np.exp(1j * phases), axis=1)
        num += np.bincount(shell, weights=np.abs(rho) ** 2, minlength=args.nq)

        # Denominator of S(q): the per-frame normalization, accumulated once
        # per frame so that num/den is the trajectory average.
        if args.norm == "n":
            dmode = np.full(len(hkl), float(natoms))
        else:
            sum_w = np.zeros(len(hkl))
            sum_w2 = np.zeros(len(hkl))
            for type_id, count in {t: np.count_nonzero(types == t) for t in set(types)}.items():
                w = amplitude_of(type_id, qlen, mapping, args.weight)
                sum_w += count * w
                sum_w2 += count * w * w
            dmode = sum_w**2 / natoms if args.norm == "mean" else sum_w2
        den += np.bincount(shell, weights=dmode, minlength=args.nq)

    s_of_q = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    q_centers = args.qmin + (np.arange(args.nq) + 0.5) * dq
    with open(args.output, "w") as handle:
        handle.write("# q S(q)\n")
        for qc, value in zip(q_centers, s_of_q):
            handle.write("%14.6f  %20.12e\n" % (qc, value))


if __name__ == "__main__":
    main()
