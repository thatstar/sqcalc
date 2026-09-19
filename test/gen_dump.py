#!/usr/bin/env python3
# sqcalc - structure factors from LAMMPS dump trajectories
# Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Generate synthetic LAMMPS dump files for testing sqcalc."""

import argparse

import numpy as np


def type_assignment(natoms, fractions, rng):
    """Assign type ids 1..len(fractions) with the requested fractions."""
    ntypes = len(fractions)
    counts = np.zeros(ntypes, dtype=int)
    remaining = natoms
    for i, frac in enumerate(fractions[:-1]):
        counts[i] = int(round(frac * natoms))
        remaining -= counts[i]
    counts[-1] = remaining
    types = np.repeat(np.arange(1, ntypes + 1), counts)
    rng.shuffle(types)
    return types


def ideal_gas(natoms, length, rng):
    return rng.random((natoms, 3)) * length


def lattice(natoms, length, rng, jitter):
    side = int(round(natoms ** (1.0 / 3.0)))
    if side ** 3 != natoms:
        raise SystemExit("--mode lattice needs a perfect cube number of atoms")
    grid = np.arange(side) * (length / side)
    pos = np.stack(np.meshgrid(grid, grid, grid, indexing="ij"), axis=-1).reshape(-1, 3)
    if jitter > 0:
        pos = pos + (rng.random(pos.shape) - 0.5) * jitter
    return pos


def ballistic(natoms, length, rng, nframes, temperature):
    """Free particles with Maxwell-Boltzmann velocities: an ideal gas in the
    ballistic regime.  The ensemble averaged coherent F(q,t) is
    exp(-q^2 <v^2> t^2 / 2), but a single finite box keeps the frozen in
    structure factor of the initial configuration as an amplitude, so the
    time averaged estimate is A exp(...) with A != 1; use --mode diffusive
    when an analytic decay has to be compared."""
    vel = rng.normal(0.0, np.sqrt(temperature), (natoms, 3))
    start = rng.random((natoms, 3)) * length
    return [start + vel * float(frame) for frame in range(nframes)]


def diffusive(natoms, length, rng, nframes, diffusivity):
    """Independent Brownian particles: F(q,t) = exp(-D q^2 t) for q L >> 1."""
    pos = rng.random((natoms, 3)) * length
    frames = [pos.copy()]
    for _ in range(nframes - 1):
        pos = pos + rng.normal(0.0, np.sqrt(2.0 * diffusivity), (natoms, 3))
        frames.append(pos.copy())
    return frames


def bounding_limits(length, tilt):
    xy, xz, yz = tilt
    xlo = 0.0 - min(0.0, xy, xz, xy + xz)
    xhi = length - max(0.0, xy, xz, xy + xz)
    ylo = 0.0 - min(0.0, yz)
    yhi = length - max(0.0, yz)
    return xlo, xhi, ylo, yhi, 0.0, length


def write_dump(path, frames, types, length, tilt, columns, step_increment=1,
               boundary="pp pp pp"):
    """Write a default style LAMMPS dump."""
    unwrapped = any(name in ("xu", "yu", "zu") for name in columns.split())
    periodic = [token.startswith("p") for token in boundary.split()]
    with open(path, "w") as handle:
        for frame, pos in enumerate(frames):
            handle.write("ITEM: TIMESTEP\n%d\n" % (frame * step_increment))
            handle.write("ITEM: NUMBER OF ATOMS\n%d\n" % len(types))
            if tilt is None:
                handle.write("ITEM: BOX BOUNDS %s\n" % boundary)
                for i in range(3):
                    handle.write("%.10f %.10f\n" % (0.0, length))
            else:
                handle.write("ITEM: BOX BOUNDS xy xz yz %s\n" % boundary)
                xlo, xhi, ylo, yhi, zlo, zhi = bounding_limits(length, tilt)
                handle.write("%.10f %.10f %.10f\n" % (xlo, xhi, tilt[0]))
                handle.write("%.10f %.10f %.10f\n" % (ylo, yhi, tilt[1]))
                handle.write("%.10f %.10f %.10f\n" % (zlo, zhi, tilt[2]))
            handle.write("ITEM: ATOMS %s\n" % columns)
            for i in range(len(types)):
                values = []
                for name in columns.split():
                    if name == "id":
                        values.append("%d" % (i + 1))
                    elif name == "type":
                        values.append("%d" % types[i])
                    elif name in ("x", "xu"):
                        values.append("%.10f" % (pos[i, 0] if unwrapped or not periodic[0]
                                                 else pos[i, 0] % length))
                    elif name in ("y", "yu"):
                        values.append("%.10f" % (pos[i, 1] if unwrapped or not periodic[1]
                                                 else pos[i, 1] % length))
                    elif name in ("z", "zu"):
                        values.append("%.10f" % (pos[i, 2] if unwrapped or not periodic[2]
                                                 else pos[i, 2] % length))
                    else:
                        raise SystemExit("unsupported column " + name)
                handle.write(" ".join(values) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--natoms", type=int, default=200)
    parser.add_argument("--length", type=float, default=20.0)
    parser.add_argument("--frames", type=int, default=5)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--mode",
                        choices=["random", "lattice", "ballistic", "diffusive"],
                        default="random")
    parser.add_argument("--fractions", default="0.5,0.5")
    parser.add_argument("--jitter", type=float, default=0.0)
    parser.add_argument("--tilt", type=float, nargs=3, default=None)
    parser.add_argument("--columns", default="id type x y z")
    parser.add_argument("--temperature", type=float, default=1.0,
                        help="--mode ballistic: Maxwell-Boltzmann temperature")
    parser.add_argument("--diffusivity", type=float, default=0.1,
                        help="--mode diffusive: diffusion coefficient")
    parser.add_argument("--step", type=int, default=1,
                        help="timestep increment written between frames")
    parser.add_argument("--boundary", default="pp pp pp",
                        help="LAMMPS periodicity tokens, e.g. 'pp pp ff' for a slab")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    fractions = [float(x) for x in args.fractions.split(",")]
    types = type_assignment(args.natoms, fractions, rng)

    if args.mode == "ballistic":
        frames = ballistic(args.natoms, args.length, rng, args.frames, args.temperature)
    elif args.mode == "diffusive":
        frames = diffusive(args.natoms, args.length, rng, args.frames, args.diffusivity)
    else:
        frames = []
        for _ in range(args.frames):
            if args.mode == "lattice":
                frames.append(lattice(args.natoms, args.length, rng, args.jitter))
            else:
                frames.append(ideal_gas(args.natoms, args.length, rng))
    write_dump(args.output, frames, types, args.length, args.tilt, args.columns,
               step_increment=args.step, boundary=args.boundary)


if __name__ == "__main__":
    main()
