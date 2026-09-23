#!/usr/bin/env python3
# sqcalc - structure factors from LAMMPS dump trajectories
# Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Independent numpy reference for the sqcalc dynamic structure factor (`sqcalc dyn`).

This is deliberately a second implementation: it parses the dump itself,
accumulates the multi-origin time correlation with its own loops and performs
the Fourier transform itself, so that the Fortran results can be cross checked
on the same trajectory.  The conventions are the ones of the dyn subcommand:

  C_ab(q,l) = <rho_a(q,t+l) rho_b*(q,t)>          with the count per lag
  F(q,tau)  = sum_ab w_a w_b C_ab(q,tau) / W(q)
  S(q,w)    = (dt/2pi) sum F exp(i w t)           (even extension)

so that sum_n S(q,w_n)*dw = F(q,0) = S(q) of the static table.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_sq import parse_mapping, amplitude_of  # noqa: E402  (sibling module)


def read_dump(path):
    """Yield (timestep, types, positions, cell) for every frame."""
    with open(path) as handle:
        lines = handle.readlines()
    i = 0
    while i < len(lines):
        if not lines[i].startswith("ITEM: TIMESTEP"):
            i += 1
            continue
        timestep = int(lines[i + 1])
        natoms = int(lines[i + 3])
        bounds = np.array([[float(v) for v in lines[i + 5 + k].split()[:2]]
                           for k in range(3)])
        tilt = np.zeros(3)
        if len(lines[i + 5].split()) > 2:
            tilt = np.array([[float(v) for v in lines[i + 5 + k].split()][2]
                             for k in range(3)])
        cell = np.diag(bounds[:, 1] - bounds[:, 0])
        cell[0, 1], cell[0, 2], cell[1, 2] = tilt
        header = lines[i + 8].split()[2:]
        start = i + 9
        values = np.array([[float(v) for v in lines[start + j].split()]
                           for j in range(natoms)])
        col = {name: k for k, name in enumerate(header)}
        types = values[:, col["type"]].astype(int)
        keys = ["xu" if "xu" in col else "x", "yu" if "yu" in col else "y",
                "zu" if "zu" in col else "z"]
        pos = values[:, [col[k] for k in keys]]
        yield timestep, types, pos, cell, bounds[:, 0]
        i = start + natoms


def unwrap_frames(frames):
    """Minimum image unwrapping, matching what the reader reconstructs."""
    out = []
    cell = frames[0][3]
    inv = np.linalg.inv(cell.T)
    prev = None
    unwrapped = None
    for timestep, types, pos, cell, origin in frames:
        frac = np.linalg.solve(cell.T, (pos - origin).T).T
        frac = frac - np.floor(frac)
        if prev is None:
            unwrapped = frac
        else:
            ds = frac - prev
            ds = ds - np.round(ds)
            unwrapped = unwrapped + ds
        prev = frac
        out.append((timestep, types, unwrapped @ cell, cell))
    return out


_FRAME_CACHE = {}


def load_frames(path):
    """Parsed, unwrapped frames of a dump, read at most once per process.

    The reference scripts are what makes the dynamic checks independent of the
    Fortran reader, but every one of them used to parse the whole trajectory
    again.  ref_batch.py runs several cases in one process, and this cache is
    what lets them share a single parse.
    """
    key = os.path.abspath(path)
    if key not in _FRAME_CACHE:
        _FRAME_CACHE[key] = unwrap_frames(list(read_dump(key)))
    return _FRAME_CACHE[key]


def correlation(frames, qvec, nspecies, species_of, maxframes, lag):
    """Multi-origin C_ab(q,l) with the count of origins per lag."""
    nq = len(qvec)
    nframes = len(frames)
    static = np.zeros((nq, nspecies, nspecies), dtype=complex)
    csum = np.zeros((nq, maxframes + 1, nspecies, nspecies), dtype=complex)
    ccnt = np.zeros((nq, maxframes + 1), dtype=int)
    rho_hist = []
    for frame_index, (timestep, types, pos, cell) in enumerate(frames):
        rho = np.zeros((nq, nspecies), dtype=complex)
        for isp in range(nspecies):
            sub = pos[species_of == isp]
            rho[:, isp] = np.sum(np.exp(1j * (sub @ qvec.T)), axis=0)
        static += rho[:, :, None] * np.conj(rho[:, None, :])
        rho_hist.append(rho)
        for l in range(min(frame_index, maxframes) + 1):
            if (frame_index - l) % lag:
                continue
            old = rho_hist[frame_index - l]
            csum[:, l, :, :] += rho[:, :, None] * np.conj(old[:, None, :])
            ccnt[:, l] += 1
    return static, csum, ccnt


def spectrum(f, dt):
    """One-sided DCT-I of a real even correlation function."""
    last = len(f) - 1
    n = np.arange(last + 1)
    out = np.zeros(last + 1)
    for k in range(last + 1):
        acc = f[0] + 2.0 * np.sum(f[1:last] * np.cos(np.pi * k * n[1:last] / last))
        acc += f[last] * (1.0 if k % 2 == 0 else -1.0)
        out[k] = acc * dt / (2.0 * np.pi)
    out[1:last] *= 2.0
    return out


def parse_dyn_q(spec):
    """q vectors of a --qpoints specification, with their quadrature weights.

    Returns (qvec, weights, radius): the weights are all 1 on a q line, and a
    shell returns its Lebedev weights plus the radius it belongs to.  The
    Lebedev grid comes from scipy, which is the reference the Fortran table
    was generated from.
    """
    spec = spec.strip()
    if spec.startswith("line:"):
        nint, s0, s1, dx, dy, dz = [float(v) for v in spec[5:].split(",")]
        nint = int(nint)
        direction = np.array([dx, dy, dz], dtype=float)
        direction /= np.linalg.norm(direction)
        scale = s0 + np.arange(nint + 1) * (s1 - s0) / nint
        return scale[:, None] * direction[None, :], np.ones(nint + 1), None
    if spec.startswith("shell:"):
        try:
            from scipy.integrate import lebedev_rule
        except ImportError as exc:  # pragma: no cover - depends on the host
            raise SystemExit("--qpoints shell needs scipy.integrate.lebedev_rule: %s" % exc)
        radius, accuracy = spec[6:].split(",")
        orders = {"low": 11, "medium": 17, "high": 23}
        x, w = lebedev_rule(orders[accuracy.strip().lower()])
        return float(radius) * x.T, w, float(radius)
    raise SystemExit("unsupported --qpoints specification %r" % spec)


def shell_average(values, weights):
    """Weighted average over the q modes of a shell (row 0 is the mode axis)."""
    return np.tensordot(weights, values, axes=(0, 0)) / weights.sum()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--qpoints", required=True, dest="dyn_q",
                        help='"-", "line:NINT,S0,S1,DX,DY,DZ" or "shell:Q,low|medium|high"')
    parser.add_argument("--maxframes", type=int, required=True)
    parser.add_argument("--lag", type=int, default=1)
    parser.add_argument("--dt", type=float, default=1.0)
    parser.add_argument("--weight", default="unit")
    parser.add_argument("--norm", default="mean")
    parser.add_argument("--mapping", default="")
    parser.add_argument("--output-fqt", required=True)
    parser.add_argument("--output-sqw", required=True)
    args = parser.parse_args()

    qvec, weights, radius = parse_dyn_q(args.dyn_q)

    frames = load_frames(args.input)
    types = frames[0][1]
    natoms = len(types)
    mapping = parse_mapping(args.mapping) if args.mapping else {}
    species = sorted(set(types.tolist()))
    species_of = np.array([species.index(t) for t in types])
    if not mapping:
        # unit weights do not need element data, but the amplitude helper looks
        # the type id up in the mapping first
        mapping = {t: "Si" for t in species}

    static, csum, ccnt = correlation(frames, qvec, len(species), species_of,
                                     args.maxframes, args.lag)
    nframes = len(frames)

    # weights and normalization, exactly like mode_denominator()
    qlen = np.linalg.norm(qvec, axis=1)
    w = np.zeros((len(qlen), len(species)))
    for isp, type_id in enumerate(species):
        w[:, isp] = amplitude_of(type_id, qlen, mapping, args.weight)
    den = np.zeros(len(qlen))
    for k, q in enumerate(qlen):
        if args.norm == "n":
            den[k] = natoms
            continue
        sum_w = sum(np.count_nonzero(types == t) * amplitude_of(t, q, mapping, args.weight)
                    for t in species)
        if args.norm == "mean":
            den[k] = sum_w ** 2 / natoms
        else:
            den[k] = sum(np.count_nonzero(types == t)
                         * amplitude_of(t, q, mapping, args.weight) ** 2 for t in species)

    ftau = np.zeros((len(qlen), args.maxframes + 1))
    for k in range(len(qlen)):
        for l in range(args.maxframes + 1):
            if l == 0:
                cab = static[k] / nframes
            else:
                cab = csum[k, l] / max(ccnt[k, l], 1)
            value = np.sum(w[k, :, None] * w[k, None, :] * np.real(cab))
            ftau[k, l] = value / den[k] if den[k] > 0 else 0.0

    dw = np.pi / (args.maxframes * args.dt)
    omega = np.arange(args.maxframes + 1) * dw
    tau = np.arange(args.maxframes + 1) * args.dt
    sqw = np.array([spectrum(ftau[k], args.dt) for k in range(len(qlen))])

    # partials in the plain (1/N) convention of the static columns
    present = set(types.tolist())
    pairs = [(a, b) for a in range(len(species)) for b in range(a, len(species))
             if species[a] in present and species[b] in present]
    labels = ["%d-%d" % (species[a], species[b]) for a, b in pairs]
    part_f = np.zeros((len(qlen), args.maxframes + 1, len(labels)))
    for p, (a, b) in enumerate(pairs):
        for k in range(len(qlen)):
            for l in range(args.maxframes + 1):
                if l == 0:
                    cab = static[k, a, b] / nframes
                else:
                    cab = csum[k, l, a, b] / max(ccnt[k, l], 1)
                    if a != b:
                        # symmetrized cross partial, as in the Fortran writer
                        cab = 0.5 * (cab + csum[k, l, b, a] / max(ccnt[k, l], 1))
                part_f[k, l, p] = np.real(cab) / natoms
    part_s = np.zeros((len(qlen), args.maxframes + 1, len(labels)))
    for k in range(len(qlen)):
        for p in range(len(labels)):
            part_s[k, :, p] = spectrum(part_f[k, :, p], args.dt)

    if radius is not None:
        # A shell is written as the single quadrature averaged row at |q| = Q.
        ftau = shell_average(ftau, weights)[None, :]
        part_f = shell_average(part_f, weights)[None, :, :]
        sqw = np.array([spectrum(ftau[0], args.dt)])
        part_s = np.zeros((1, args.maxframes + 1, len(labels)))
        for p in range(len(labels)):
            part_s[0, :, p] = spectrum(part_f[0, :, p], args.dt)
        qvec = np.array([[0.0, 0.0, radius]])

    def write(path, axis, axis_name, values, pvalues, total_name):
        with open(path, "w") as handle:
            handle.write("# %s\n" % total_name)
            handle.write("# qx qy qz %s %s %s\n"
                         % (axis_name, total_name, " ".join(labels)))
            for k in range(len(qvec)):
                for l in range(args.maxframes + 1):
                    handle.write("%14.8f  %14.8f  %14.8f  %16.8f  %20.12e"
                                 % (qvec[k, 0], qvec[k, 1], qvec[k, 2], axis[l], values[k, l]))
                    for i in range(len(labels)):
                        handle.write("  %20.12e" % pvalues[k, l, i])
                    handle.write("\n")

    write(args.output_fqt, tau, "tau", ftau, part_f, "F(q,t)")
    write(args.output_sqw, omega, "omega", sqw, part_s, "S(q,w)")


if __name__ == "__main__":
    main()
