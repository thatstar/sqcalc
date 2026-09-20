#!/usr/bin/env python3
# sqcalc - structure factors from LAMMPS dump trajectories
# Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Independent numpy reference for the self intermediate scattering function.

The conventions match the --fqt-self documentation:

  F_s(q,t)   = N^-1 <sum_i exp(i q.(r_i(t0+t)-r_i(t0)))>
  F_s^a(q,t) = N^-1 <sum_{i in a} exp(i q.(r_i(t0+t)-r_i(t0)))>

with unit weights and the same frame stride / time-origin schedule as the
Fortran implementation.  The total is the sum of the species columns.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_sqw import read_dump, unwrap_frames  # noqa: E402  (sibling module)


def self_function(frames, qvec, maxframes, lag, stride):
    """Return F_s(q,lag), the per-species parts and the origin counts."""
    types = frames[0][1]
    natoms = len(types)
    species = sorted(set(types.tolist()))
    species_of = np.array([species.index(t) for t in types])
    nq = len(qvec)
    nsteps = maxframes // stride
    fs_sum = np.zeros((nq, nsteps + 1, len(species)), dtype=complex)
    counts = np.zeros(nsteps + 1, dtype=int)

    for frame_index, (_, _, pos, _) in enumerate(frames):
        if frame_index % stride:
            continue
        sample = frame_index // stride
        for ell in range(min(sample, nsteps) + 1):
            if (frame_index - ell * stride) % lag:
                continue
            origin = frames[frame_index - ell * stride][2]
            disp = pos - origin
            counts[ell] += 1
            for isp in range(len(species)):
                sub = disp[species_of == isp]
                if len(sub):
                    fs_sum[:, ell, isp] += np.exp(1j * (sub @ qvec.T)).sum(axis=0)

    partial = np.zeros((nq, nsteps + 1, len(species)))
    for ell in range(nsteps + 1):
        n = counts[ell]
        if n > 0:
            partial[:, ell, :] = np.real(fs_sum[:, ell, :]) / n / natoms
    total = partial.sum(axis=2)
    return total, partial, species, counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--dyn", required=True, help="NINT,S0,S1,DX,DY,DZ")
    parser.add_argument("--maxframes", type=int, required=True)
    parser.add_argument("--lag", type=int, default=1)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--dt", type=float, default=1.0)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    nint, s0, s1, dx, dy, dz = [float(v) for v in args.dyn.split(",")]
    nint = int(nint)
    direction = np.array([dx, dy, dz], dtype=float)
    direction /= np.linalg.norm(direction)
    scale = s0 + np.arange(nint + 1) * (s1 - s0) / nint
    qvec = scale[:, None] * direction[None, :]

    frames = unwrap_frames(list(read_dump(args.input)))
    total, partial, species, _ = self_function(frames, qvec, args.maxframes, args.lag,
                                               args.stride)
    nsteps = args.maxframes // args.stride
    tau = np.arange(nsteps + 1) * args.stride * args.dt
    with open(args.output, "w") as handle:
        handle.write("# sqcalc 0.1.0 self intermediate scattering function F_s(q,t)\n")
        labels = " ".join("F_s(%d)" % s for s in species)
        handle.write("# qx qy qz tau F_s(q,t) %s\n" % labels)
        for k in range(len(qvec)):
            for ell, time in enumerate(tau):
                handle.write("%14.8f  %14.8f  %14.8f  %16.8f  %20.12e"
                             % (qvec[k, 0], qvec[k, 1], qvec[k, 2], time, total[k, ell]))
                for isp in range(len(species)):
                    handle.write("  %20.12e" % partial[k, ell, isp])
                handle.write("\n")


if __name__ == "__main__":
    main()
