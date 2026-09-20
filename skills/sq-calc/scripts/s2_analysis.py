#!/usr/bin/env python3
# sqcalc - structure factors from LAMMPS dump trajectories
# Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Analyse sqcalc pair-entropy outputs.

The Debye method writes the raw S2(r) accumulation curve, which is only a
finite-rmax estimate.  This script adds the model-dependent post-processing:

  * optional GCV smoothing of the partial g(r) before re-integrating,
  * Richardson extrapolation over several --dr values,
  * tail fits (power law, exponential, damped oscillation) of S2(r),
  * optional multi-box extrapolation.

No statistical error is estimated: S2 is treated as a qualitative/trend
quantity, and the dominant uncertainties are the rmax tail and finite-size
effects.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

import numpy as np


def read_lines(path: Path) -> list[str]:
    return path.read_text().splitlines()


def parse_metadata(headers: list[str]) -> dict:
    text = " ".join(headers)
    meta: dict = {}
    for key in ("natoms", "nframes"):
        match = re.search(r"\b%s\s+(\d+)" % key, text)
        if match:
            meta[key] = int(match.group(1))
    for key in ("volume", "rmax", "dr"):
        match = re.search(r"\b%s\s+([0-9eE.+-]+)" % key, text)
        if match:
            meta[key] = float(match.group(1))
    match = re.search(r"#\s*counts\s+([0-9 ]+)", text)
    if match:
        meta["counts"] = [int(v) for v in match.group(1).split()]
    return meta


def read_pair_entropy(path: Path) -> tuple[dict, dict]:
    lines = read_lines(path)
    headers = [line for line in lines if line.startswith("#")]
    meta = parse_metadata(headers)
    final: dict[str, float] = {}
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        if len(fields) >= 2:
            final[fields[0]] = float(fields[1])
    return meta, final


def read_accum(path: Path) -> tuple[dict, list[str], np.ndarray]:
    lines = read_lines(path)
    headers = [line for line in lines if line.startswith("#")]
    meta = parse_metadata(headers)
    labels: list[str] = []
    for line in headers:
        names = line.lstrip("#").split()
        if names and names[0] == "r" and len(names) >= 2:
            labels = names
    rows = []
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        rows.append([float(v) for v in stripped.split()])
    return meta, labels, np.array(rows)


def read_rdf(path: Path) -> tuple[list[str], np.ndarray]:
    lines = read_lines(path)
    labels: list[str] = []
    for line in lines:
        if line.startswith("#"):
            names = line.lstrip("#").split()
            if names and names[0] == "r" and len(names) >= 2:
                labels = names
    rows = []
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        rows.append([float(v) for v in stripped.split()])
    return labels, np.array(rows)


def triangle_pairs(ntypes: int) -> list[tuple[int, int]]:
    return [(ia, ib) for ia in range(ntypes) for ib in range(ia, ntypes)]


def ntypes_from_npairs(npair: int) -> int:
    return int((math.isqrt(1 + 8 * npair) - 1) // 2)


def shell_volume(r: np.ndarray, dr: float) -> np.ndarray:
    rlo = np.maximum(0.0, r - 0.5 * dr)
    rhi = r + 0.5 * dr
    return 4.0 * np.pi / 3.0 * (rhi**3 - rlo**3)


def integrate_s2(r, g_partial, meta, bias_correction=True):
    """Integrate S2 from partial g(r) with the exact shell volume."""
    counts = meta.get("counts")
    if counts is None:
        raise SystemExit("missing '# counts' metadata; use the sqcalc output")
    natoms = float(meta["natoms"])
    volume = float(meta["volume"])
    nframes = float(meta["nframes"])
    dr = float(meta["dr"])
    ntypes = len(counts)
    pairs = triangle_pairs(ntypes)
    npair = len(pairs)
    if g_partial.shape[1] < npair:
        raise SystemExit("the RDF file has fewer partial columns than expected")
    rho = natoms / volume
    dv = shell_volume(r, dr)
    curves = np.zeros((len(r), npair))
    final = np.zeros(npair)
    for p, (ia, ib) in enumerate(pairs):
        if (ia == ib and counts[ia] <= 1) or (ia != ib and (counts[ia] == 0 or counts[ib] == 0)):
            continue
        g = g_partial[:, p]
        integrand = np.where(g > 0.0, g * np.log(np.maximum(g, 1.0e-300)) - g + 1.0, 1.0)
        contrib = dv * integrand
        if bias_correction:
            bias = volume / (2.0 * nframes * counts[ia] * counts[ib])
            contrib = contrib - bias
        cum = np.cumsum(contrib)
        curves[:, p] = -2.0 * math.pi * rho * (counts[ia] / natoms) * (counts[ib] / natoms) * cum
        final[p] = curves[-1, p]
    total_curve = np.zeros(len(r))
    total = 0.0
    for p in range(npair):
        ia, ib = pairs[p]
        weight = 1.0 if ia == ib else 2.0
        total_curve += weight * curves[:, p]
        total += weight * final[p]
    return r, curves, total_curve, total, final


def smooth_gcv(r, y):
    """GCV smoothing spline (SciPy); fail with a clear message if missing."""
    try:
        from scipy.interpolate import make_smoothing_spline
    except Exception as exc:  # pragma: no cover - dependency error path
        raise SystemExit("--s2-smooth needs scipy (scipy.interpolate): %s" % exc)
    order = np.argsort(r)
    r = np.asarray(r)[order]
    y = np.asarray(y)[order]
    spline = make_smoothing_spline(r, y)
    return spline(r)


def fit_tail(r, y, kind, rmin, rmax):
    """Fit the tail models; return (S_inf, parameters, rms)."""
    try:
        from scipy.optimize import curve_fit
    except Exception as exc:  # pragma: no cover - dependency error path
        raise SystemExit("tail fits need scipy (scipy.optimize): %s" % exc)
    mask = (r >= rmin) & (r <= rmax)
    x = r[mask]
    z = y[mask]
    if len(x) < 4:
        raise SystemExit("the fit range has too few points")

    if kind == "power":
        def model(xx, sinf, amp, power):
            return sinf - amp * xx ** (-power)
        p0 = (z[0], max(1.0e-12, abs(z[0]) * x[0] ** 9), 9.0)
        bounds = ([-np.inf, 0.0, 0.0], [np.inf, np.inf, 30.0])
    elif kind == "exp":
        def model(xx, sinf, amp, xi):
            return sinf - amp * np.exp(-2.0 * xx / xi)
        p0 = (z[0], max(1.0e-12, abs(z[0])), max(1.0, x[-1]))
        bounds = ([-np.inf, 0.0, 1.0e-6], [np.inf, np.inf, np.inf])
    elif kind == "damped":
        def model(xx, sinf, amp, xi, k, phi):
            return sinf - amp * np.exp(-2.0 * xx / xi) * np.cos(2.0 * k * xx + phi)
        p0 = (z[0], max(1.0e-12, abs(z[0])), max(1.0, x[-1]), 1.0, 0.0)
        bounds = ([-np.inf, 0.0, 1.0e-6, 0.0, -np.pi], [np.inf, np.inf, np.inf, np.inf, np.pi])
    else:
        raise SystemExit("unknown fit %s" % kind)

    popt, _ = curve_fit(model, x, z, p0=p0, bounds=bounds, maxfev=20000)
    rms = float(np.sqrt(np.mean((model(x, *popt) - z) ** 2)))
    return float(popt[0]), popt, rms, model


def richardson(values):
    """values = [(dr, S2), ...]; return Richardson estimates."""
    values = sorted(values, key=lambda item: item[0], reverse=True)
    out = {}
    if len(values) >= 2:
        d1, s1 = values[0]
        d2, s2 = values[1]
        if abs(d2 / d1 - 0.5) > 0.05:
            print("warning: two-point Richardson needs dr/2; got dr = %.6g and %.6g"
                  % (d1, d2))
        else:
            out["two_point"] = (4.0 * s2 - s1) / 3.0
    if len(values) >= 3:
        d1, s1 = values[0]
        d2, s2 = values[1]
        d3, s3 = values[2]
        if abs(d2 / d1 - 0.5) > 0.05 or abs(d3 / d2 - 0.5) > 0.05:
            print("warning: three-point Richardson needs dr/2 and dr/4; got "
                  "dr = %.6g, %.6g, %.6g" % (d1, d2, d3))
        else:
            out["three_point"] = (64.0 * s3 - 20.0 * s2 + s1) / 45.0
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--s2", type=Path, help="sqcalc --pair-entropy text file")
    parser.add_argument("--rdf", type=Path, help="sqcalc --rdf text file (for --s2-smooth)")
    parser.add_argument("--accum", type=Path, nargs="+", required=True,
                        help="one or more --s2-accum files (multiple: --dr sweep)")
    parser.add_argument("--s2-smooth", action="store_true",
                        help="GCV-smooth the partial g(r) from --rdf and re-integrate")
    parser.add_argument("--fit", choices=["power", "exp", "damped"], default="power")
    parser.add_argument("--r-fit-min", type=float, default=None)
    parser.add_argument("--r-fit-max", type=float, default=None)
    parser.add_argument("--json", type=Path, help="write the results as JSON")
    parser.add_argument("--plot", type=Path, help="write a PNG of S2(r) and the fit")
    args = parser.parse_args()

    results: dict = {"accum": []}
    dr_values = []
    first_meta = None
    first_labels = None
    first_curve = None
    for path in args.accum:
        meta, labels, rows = read_accum(path)
        if first_meta is None:
            first_meta, first_labels, first_curve = meta, labels, rows
        r = rows[:, 0]
        total_curve = rows[:, 1]
        dr_values.append((float(meta["dr"]), float(total_curve[-1])))
        results["accum"].append({"path": str(path), "dr": meta["dr"],
                                 "S2_rmax": float(total_curve[-1])})

    r = first_curve[:, 0]
    total_curve = first_curve[:, 1]

    if args.s2_smooth:
        if args.s2 is None or args.rdf is None:
            raise SystemExit("--s2-smooth needs both --s2 and --rdf")
        meta, _ = read_pair_entropy(args.s2)
        rdf_labels, rdf_rows = read_rdf(args.rdf)
        r_rdf = rdf_rows[:, 0]
        dr = float(meta["dr"])
        dv = shell_volume(r_rdf, dr)
        # Convert the midpoint-shell g(r) of --rdf to the exact-shell g used
        # by the pair-entropy integral.
        shell_mid = 4.0*np.pi*r_rdf**2*dr
        g_partial = np.array(rdf_rows[:, 2:])*(shell_mid/dv)[:, None]
        for p in range(g_partial.shape[1]):
            g_partial[:, p] = smooth_gcv(r_rdf, g_partial[:, p])
            g_partial[:, p] = np.maximum(g_partial[:, p], 0.0)
        r_s, curves, total_curve, total, final = integrate_s2(
            r_rdf, g_partial, meta, bias_correction=False)
        results["gcv_smooth"] = {"S2_rmax": float(total)}
        r, total_curve = r_s, total_curve

    results["richardson"] = richardson(dr_values)

    rmin = args.r_fit_min if args.r_fit_min is not None else float(0.5 * r[-1])
    rmax = args.r_fit_max if args.r_fit_max is not None else float(r[-1])
    sinf, popt, rms, fit_model = fit_tail(r, total_curve, args.fit, rmin, rmax)
    results["tail_fit"] = {"model": args.fit, "S2_inf": sinf,
                           "parameters": [float(v) for v in popt], "rms": rms,
                           "r_fit_min": rmin, "r_fit_max": rmax}

    if args.json:
        args.json.write_text(json.dumps(results, indent=2))

    print("raw S2(rmax)      : %.6f" % float(total_curve[-1]))
    if "two_point" in results["richardson"]:
        print("Richardson (2 pt) : %.6f" % results["richardson"]["two_point"])
    if "three_point" in results["richardson"]:
        print("Richardson (3 pt) : %.6f" % results["richardson"]["three_point"])
    print("%-18s: %.6f (rms %.3g)" % ("tail fit " + args.fit, sinf, rms))

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fit_r = np.linspace(rmin, rmax, 200)
        fig, ax = plt.subplots(figsize=(5, 4), dpi=140)
        ax.plot(r, total_curve, label=r"$S_2(r)$")
        ax.plot(fit_r, fit_model(fit_r, *popt), "r--", label="%s fit" % args.fit)
        ax.axvline(rmin, color="k", ls=":", lw=0.8)
        ax.set_xlabel("r [A]")
        ax.set_ylabel(r"$S_2(r)/k_B$")
        ax.legend()
        fig.tight_layout()
        fig.savefig(args.plot)
        print("wrote %s" % args.plot)
    return 0


if __name__ == "__main__":
    sys.exit(main())
