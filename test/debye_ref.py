#!/usr/bin/env python3
"""Independent numpy reference for the sqcalc Debye method.

Builds the same partial pair histograms (minimum image, ordered pairs), the same
analytic r = 0 self term and the same normalization, so the Fortran result can be
checked column by column.  It also reports the exact (unbinned) Debye sum, which
quantifies the binning error.
"""

import argparse

import numpy as np

import ref_sq


def minimum_image(frac, rmax, cell):
    """Distances of all ordered, non-identical pairs within rmax."""
    natoms = len(frac)
    frac = frac - np.floor(frac)
    dr = frac[:, None, :] - frac[None, :, :]
    dr = dr - np.rint(dr)
    dist = np.linalg.norm(dr @ cell.T, axis=-1)
    idx = np.where(dist <= rmax)
    keep = idx[0] != idx[1]
    return idx[0][keep], idx[1][keep], dist[idx][keep]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--mapping", default="1:Si,2:O")
    parser.add_argument("--weight", default="unit", choices=["unit", "neutron", "xray"])
    parser.add_argument("--norm", default="mean", choices=["mean", "self", "n"])
    parser.add_argument("--qmin", type=float, default=0.0)
    parser.add_argument("--qmax", type=float, default=20.0)
    parser.add_argument("--nq", type=int, default=500)
    parser.add_argument("--rmax", type=float, required=True)
    parser.add_argument("--dr", type=float, default=0.01)
    parser.add_argument("--no-correction", action="store_true",
                        help="skip debyer's cut-off density correction")
    parser.add_argument("--output", required=True)
    parser.add_argument("--rdf", default=None, help="also write the total g(r)")
    args = parser.parse_args()

    mapping = ref_sq.parse_mapping(args.mapping)
    nbins = int(np.ceil(args.rmax / args.dr))
    nq = args.nq
    fq = args.qmin + (np.arange(nq) + 0.5) * (args.qmax - args.qmin) / nq

    frames = list(ref_sq.read_dump(args.input))
    nframes = len(frames)
    ntypes = max(max(mapping), max(int(t.max()) for t, _, _, _ in frames))
    hist = np.zeros((ntypes + 1, ntypes + 1, nbins))
    counts = np.zeros(ntypes + 1, dtype=int)
    volume = abs(np.linalg.det(frames[0][2]))
    exact = np.zeros((nq, ntypes + 1, ntypes + 1))

    for types, pos, cell, origin in frames:
        frac = np.linalg.solve(cell.T, (pos - origin).T).T
        i, j, dist = minimum_image(frac, args.rmax, cell)
        bin_idx = (dist / args.dr).astype(int)
        np.add.at(hist, (types[i], types[j], bin_idx), 1.0)
        for k in range(nq):
            exact[k] += np.bincount(types[i] * (ntypes + 1) + types[j],
                                    weights=np.sinc(fq[k] * dist / np.pi), minlength=(ntypes + 1) ** 2
                                    ).reshape(ntypes + 1, ntypes + 1)
        for t in range(1, ntypes + 1):
            counts[t] = np.count_nonzero(types == t)

    hist /= nframes
    exact /= nframes
    r_centers = (np.arange(nbins) + 0.5) * args.dr

    def weights(q):
        return np.array([0.0] + [ref_sq.amplitude_of(t, q, mapping, args.weight)
                                 for t in range(1, ntypes + 1)])

    s_of_q = np.zeros(nq)
    s_exact = np.zeros(nq)
    for k in range(nq):
        w = weights(fq[k])
        sinc = np.sinc(fq[k] * r_centers / np.pi)
        pair_sum = np.einsum("a,b,abm,m->", w, w, hist, sinc)
        self_sum = float(np.sum(counts * w**2))
        weight_sum = float(np.sum(counts * w))
        natoms = float(counts.sum())
        if args.norm == "n":
            den = natoms
        elif args.norm == "mean":
            den = weight_sum**2 / natoms
        else:
            den = float(np.sum(counts * w**2))
        correction = 0.0
        if not args.no_correction:
            # debyer's add_cutoff_correction (per-atom normalization, so the
            # numerator gains corr * N)
            avg = weight_sum / natoms
            rho0 = natoms / volume
            q = fq[k]
            correction = (avg**2 * 4.0 * np.pi * rho0 / q**2
                          * (args.rmax * np.cos(q * args.rmax) - np.sin(q * args.rmax) / q)
                          * natoms)
        s_of_q[k] = (pair_sum + self_sum + correction) / den
        s_exact[k] = (np.einsum("a,ab,b->", w, exact[k], w) + self_sum + correction) / den

    with open(args.output, "w") as handle:
        handle.write("# q S(q)\n")
        for qc, value in zip(fq, s_of_q):
            handle.write("%14.6f  %20.12e\n" % (qc, value))

    worst = np.max(np.abs(s_of_q - s_exact) / np.maximum(np.abs(s_exact), 1.0))
    print("debye_ref: binning error %.2e (dr=%g)" % (worst, args.dr))

    if args.rdf:
        # total g(r) with the same weighting as the Fortran writer (q = 0)
        w0 = weights(0.0)
        shell = 4.0 * np.pi * r_centers**2 * args.dr
        g = np.zeros((ntypes + 1, ntypes + 1, nbins))
        for a in range(1, ntypes + 1):
            for b in range(1, ntypes + 1):
                if counts[a] == 0 or counts[b] == 0:
                    continue
                rho_b = counts[b] / volume
                ideal = counts[a] * shell * rho_b
                good = ideal > 0
                g[a, b, good] = hist[a, b, good] / ideal[good]
        c = counts / counts.sum()
        f_mean = float(np.sum(c * w0))
        g_tot = np.einsum("a,b,abm->m", c * w0, c * w0, g) / f_mean**2 if f_mean else np.zeros(nbins)
        with open(args.rdf, "w") as handle:
            handle.write("# r g(r)\n")
            for rc, value in zip(r_centers, g_tot):
                handle.write("%12.6f  %18.10e\n" % (rc, value))


if __name__ == "__main__":
    main()
