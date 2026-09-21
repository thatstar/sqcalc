#!/usr/bin/env python3
# sqcalc - structure factors from LAMMPS dump trajectories
# Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Independent numpy reference for the sqcalc four-point structure factor.

The conventions match the --s4/--chi4 documentation:

  w_i(t,t0) = theta(a - |r_i(t0+t) - r_i(t0)|)
  W(q,t)    = sum_i exp(i q.r_i(t0)) w_i(t,t0)
  S4(q,t)   = N^-1 (<W(q,t) W(-q,t)> - |<W(q,t)>|^2)
  Q(t)      = N^-1 <sum_i w_i(t,t0)>
  chi4(t)   = N^-1 (<W0(t)^2> - <W0(t)>^2),  W0(t) = sum_i w_i(t,t0)

The connected averages use the unbiased (n-1) denominator over the time
origins, exactly like the Fortran estimator.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_sqw import read_dump, unwrap_frames, parse_dyn_q, shell_average  # noqa: E402


def four_point(frames, qvec, cutoff, maxframes, lag, stride):
    """Return S4(q,lag), Q(lag), chi4(lag) and the origin counts."""
    nq = len(qvec)
    natoms = len(frames[0][1])
    nsteps = maxframes // stride
    s4_a = np.zeros((nq, nsteps + 1))
    s4_b = np.zeros((nq, nsteps + 1), dtype=complex)
    chi4_a = np.zeros(nsteps + 1)
    chi4_b = np.zeros(nsteps + 1)
    counts = np.zeros(nsteps + 1, dtype=int)
    cutoff2 = cutoff * cutoff

    for frame_index, (_, _, pos, _) in enumerate(frames):
        if frame_index % stride:
            continue
        sample = frame_index // stride
        for ell in range(min(sample, nsteps) + 1):
            if (frame_index - ell * stride) % lag:
                continue
            origin = frames[frame_index - ell * stride][2]
            displacement = pos - origin
            inside = np.einsum("ij,ij->i", displacement, displacement) <= cutoff2
            nover = int(np.count_nonzero(inside))
            counts[ell] += 1
            chi4_b[ell] += nover
            chi4_a[ell] += nover * nover
            if nover == 0:
                continue
            phase = origin[inside] @ qvec.T
            w = np.exp(1j * phase).sum(axis=0)
            s4_b[:, ell] += w
            s4_a[:, ell] += np.abs(w) ** 2

    s4 = np.zeros((nq, nsteps + 1))
    overlap = np.zeros(nsteps + 1)
    chi4 = np.zeros(nsteps + 1)
    for ell in range(nsteps + 1):
        n = counts[ell]
        if n > 0:
            overlap[ell] = (chi4_b[ell] / n) / natoms
        if n > 1:
            s4[:, ell] = (s4_a[:, ell] - np.abs(s4_b[:, ell]) ** 2 / n) / (n - 1) / natoms
            chi4[ell] = (chi4_a[ell] - chi4_b[ell] ** 2 / n) / (n - 1) / natoms
    return s4, overlap, chi4, counts


def write_s4(path, qvec, tau, s4):
    with open(path, "w") as handle:
        handle.write("# sqcalc 0.1.0 four-point structure factor S4(q,t)\n")
        handle.write("# qx qy qz tau S4(q,t)\n")
        for k in range(len(qvec)):
            for ell, time in enumerate(tau):
                handle.write("%14.8f  %14.8f  %14.8f  %16.8f  %20.12e\n"
                             % (qvec[k, 0], qvec[k, 1], qvec[k, 2], time, s4[k, ell]))


def write_chi4(path, tau, overlap, chi4):
    with open(path, "w") as handle:
        handle.write("# sqcalc 0.1.0 average overlap Q(t) and four-point "
                     "susceptibility chi4(t)\n")
        handle.write("# tau Q(t) chi4(t)\n")
        for ell, time in enumerate(tau):
            handle.write("%16.8f  %20.12e  %20.12e\n"
                         % (time, overlap[ell], chi4[ell]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--dyn-q", required=True, dest="dyn_q",
                        help='"-", "line:NINT,S0,S1,DX,DY,DZ" or "shell:Q,low|medium|high"')
    parser.add_argument("--maxframes", type=int, required=True)
    parser.add_argument("--lag", type=int, default=1)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--dt", type=float, default=1.0)
    parser.add_argument("--cutoff", type=float, required=True)
    parser.add_argument("--output-s4", default=None)
    parser.add_argument("--output-chi4", default=None)
    args = parser.parse_args()

    qvec, weights, radius = parse_dyn_q(args.dyn_q)

    frames = unwrap_frames(list(read_dump(args.input)))
    nsteps = args.maxframes // args.stride
    s4, overlap, chi4, _ = four_point(frames, qvec, args.cutoff, args.maxframes, args.lag,
                                      args.stride)
    if radius is not None:
        s4 = shell_average(s4, weights)[None, :]
        qvec = np.array([[0.0, 0.0, radius]])
    tau = np.arange(nsteps + 1) * args.stride * args.dt
    if args.output_s4:
        write_s4(args.output_s4, qvec, tau, s4)
    if args.output_chi4:
        write_chi4(args.output_chi4, tau, overlap, chi4)


if __name__ == "__main__":
    main()
