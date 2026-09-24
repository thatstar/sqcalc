#!/usr/bin/env python3
# sqcalc - structure factors from LAMMPS dump trajectories
# Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Independent numpy reference for the non-Gaussian parameter.

The conventions match the --ngp documentation:

  alpha2(t)   = 3 <r^4> / (5 <r^2>^2) - 1
  alpha2^a(t) = 3 <r^4>_a / (5 <r^2>_a^2) - 1

with unit weights, the moments averaged over the same frame stride / time
origins as the Fortran implementation, and <...>_a the species-restricted
moments (normalized by the species count N_a, so that this is the ratio of the
species moments and not the MSD partial rescaled).  This is the ratio of the
averaged moments, not the average of the per-origin ratios, and alpha2(0) is 0
because the displacement vanishes.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_sqw import read_dump, unwrap_frames  # noqa: E402


def non_gaussian(frames, maxframes, lag, stride):
    """Return alpha2(t), the per-species parts and the origin counts."""
    types = frames[0][1]
    natoms = len(types)
    species = sorted(set(types.tolist()))
    atoms_of = np.array([np.sum(types == s) for s in species])
    species_of = np.array([species.index(t) for t in types])
    nsteps = maxframes // stride
    m2_sum = np.zeros((len(species), nsteps + 1))
    m4_sum = np.zeros((len(species), nsteps + 1))
    counts = np.zeros(nsteps + 1, dtype=int)

    for frame_index, (_, _, pos, _) in enumerate(frames):
        if frame_index % stride:
            continue
        sample = frame_index // stride
        for ell in range(min(sample, nsteps) + 1):
            if (frame_index - ell * stride) % lag:
                continue
            origin = frames[frame_index - ell * stride][2]
            r2 = ((pos - origin) ** 2).sum(axis=1)
            counts[ell] += 1
            for isp in range(len(species)):
                mask = species_of == isp
                m2_sum[isp, ell] += r2[mask].sum()
                m4_sum[isp, ell] += (r2[mask] ** 2).sum()

    partial = np.zeros((len(species), nsteps + 1))
    total = np.zeros(nsteps + 1)
    for ell in range(nsteps + 1):
        n = counts[ell]
        if n == 0:
            continue
        m2 = m2_sum[:, ell].sum() / n / natoms
        m4 = m4_sum[:, ell].sum() / n / natoms
        if m2 > 0.0:
            total[ell] = 3.0 * m4 / (5.0 * m2 * m2) - 1.0
        for isp in range(len(species)):
            m2a = m2_sum[isp, ell] / n / atoms_of[isp]
            m4a = m4_sum[isp, ell] / n / atoms_of[isp]
            if m2a > 0.0:
                partial[isp, ell] = 3.0 * m4a / (5.0 * m2a * m2a) - 1.0
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
    total, partial, species, _ = non_gaussian(
        frames, args.maxframes, args.lag, args.stride)
    nsteps = args.maxframes // args.stride
    tau = np.arange(nsteps + 1) * args.stride * args.dt
    with open(args.output, "w") as handle:
        handle.write("# sqcalc 0.1.0 non-Gaussian parameter alpha2(t)\n")
        labels = " ".join("alpha2(%d)" % s for s in species)
        handle.write("# tau alpha2(t) %s\n" % labels)
        for ell, time in enumerate(tau):
            handle.write("%16.8f  %20.12e" % (time, total[ell]))
            for isp in range(len(species)):
                handle.write("  %20.12e" % partial[isp, ell])
            handle.write("\n")


if __name__ == "__main__":
    main()
