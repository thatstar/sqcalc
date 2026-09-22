#!/usr/bin/env python3
# sqcalc - structure factors from LAMMPS dump trajectories
# Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Independent numpy reference for the mean squared displacement.

The conventions match the --msd documentation:

  MSD(t)   = N^-1 <sum_i |r_i(t0+t) - r_i(t0)|^2>_t0
  MSD^a(t) = N^-1 <sum_{i in a} |r_i(t0+t) - r_i(t0)|^2>_t0

with unit weights and the same frame stride / time-origin schedule as the
Fortran implementation.  The total is the sum of the species columns, so
MSD^a / x_a tends to 6 D_a t at long times.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_sqw import read_dump, unwrap_frames  # noqa: E402


def mean_squared_displacement(frames, maxframes, lag, stride):
    """Return MSD(t), the per-species parts and the origin counts."""
    types = frames[0][1]
    natoms = len(types)
    species = sorted(set(types.tolist()))
    species_of = np.array([species.index(t) for t in types])
    nsteps = maxframes // stride
    msd_sum = np.zeros((len(species), nsteps + 1))
    counts = np.zeros(nsteps + 1, dtype=int)

    for frame_index, (_, _, pos, _) in enumerate(frames):
        if frame_index % stride:
            continue
        sample = frame_index // stride
        for ell in range(min(sample, nsteps) + 1):
            if (frame_index - ell * stride) % lag:
                continue
            origin = frames[frame_index - ell * stride][2]
            disp2 = ((pos - origin) ** 2).sum(axis=1)
            counts[ell] += 1
            for isp in range(len(species)):
                msd_sum[isp, ell] += disp2[species_of == isp].sum()

    partial = np.zeros((len(species), nsteps + 1))
    for ell in range(nsteps + 1):
        if counts[ell] > 0:
            partial[:, ell] = msd_sum[:, ell] / counts[ell] / natoms
    total = partial.sum(axis=0)
    return total, partial, species, counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--maxframes", type=int, required=True)
    parser.add_argument("--lag", type=int, default=1)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--dt", type=float, default=1.0)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    frames = unwrap_frames(list(read_dump(args.input)))
    total, partial, species, _ = mean_squared_displacement(
        frames, args.maxframes, args.lag, args.stride)
    nsteps = args.maxframes // args.stride
    tau = np.arange(nsteps + 1) * args.stride * args.dt
    with open(args.output, "w") as handle:
        handle.write("# sqcalc 0.1.0 mean squared displacement MSD(t)\n")
        labels = " ".join("MSD(%d)" % s for s in species)
        handle.write("# tau MSD(t) %s\n" % labels)
        for ell, time in enumerate(tau):
            handle.write("%16.8f  %20.12e" % (time, total[ell]))
            for isp in range(len(species)):
                handle.write("  %20.12e" % partial[isp, ell])
            handle.write("\n")


if __name__ == "__main__":
    main()
