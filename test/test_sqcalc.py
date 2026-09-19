#!/usr/bin/env python3
# sqcalc - structure factors from LAMMPS dump trajectories
# Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""End to end checks for the sqcalc executable.

The tests combine independent views of the same physics:

1. sqcalc with the NUFFT method (the production path),
2. sqcalc with the direct summation method,
3. sqcalc with the Debye pair-histogram method,
4. small numpy reference implementations (test/ref_sq.py, test/debye_ref.py),

plus qualitative checks (ideal gas -> S ~ 1, cubic lattice -> Bragg peaks),
the HDF5 and partial structure factor outputs, the g(r) output and a handful of
command line error cases.
"""

import argparse
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def run(cmd):
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise SystemExit("command failed (%d): %s\n%s\n%s"
                         % (result.returncode, " ".join(cmd), result.stdout, result.stderr))
    return result


def read_table(path):
    q, s = [], []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            # the S(q) table may carry partial structure factor columns after
            # the total; the comparison helpers only use q and S(q)
            values = line.split()
            q.append(float(values[0]))
            s.append(float(values[1]))
    return np.array(q), np.array(s)


def read_matrix(path):
    """Read a whitespace separated data table (any number of columns)."""
    rows = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            rows.append([float(value) for value in line.split()])
    return np.array(rows)


def compare(reference, candidate, label, rtol=1.0e-6, atol=1.0e-8):
    """Compare two S(q) tables and fail when they disagree."""
    scale = np.maximum(np.abs(reference), 1.0)
    worst = np.max(np.abs(candidate - reference) / scale)
    if worst > rtol:
        raise SystemExit("FAIL %s: max relative deviation %.3e exceeds %.1e"
                         % (label, worst, rtol))
    print("  ok   %-42s max deviation %.2e" % (label, worst))
    return worst


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--gpu", action="store_true",
                        help="also compare the CUDA backend against the CPU")
    parser.add_argument("--h5read", default=None,
                        help="path of the h5read helper (enables the HDF5 checks)")
    args = parser.parse_args()
    exe = os.path.abspath(args.executable)
    work = os.path.abspath(args.workdir)
    os.makedirs(work, exist_ok=True)

    def path(name):
        return os.path.join(work, name)

    def generate(out, **kwargs):
        cmd = [sys.executable, os.path.join(HERE, "gen_dump.py"), "--output", path(out)]
        for key, value in kwargs.items():
            flag = "--" + key.replace("_", "-")
            # --tilt takes three numbers, every other option a single token.
            cmd += [flag] + str(value).split() if key == "tilt" else [flag, str(value)]
        run(cmd)
        return path(out)

    def sqcalc(dump, out, *options):
        cmd = [exe, "-i", dump, "-m", "1:Si,2:O", *options, path(out)]
        return run(cmd)

    common = ["--qmax", "6", "--nq", "120", "--eps", "1e-10"]
    failures = 0

    # --- 1. ideal gas, unit weights: NUFFT vs direct vs numpy --------------
    print("ideal gas, unit weights")
    dump = generate("gas.dump", natoms=250, length=20, frames=4, seed=7)
    sqcalc(dump, "gas_nufft.dat", "-w", "unit", "--norm", "self", *common)
    sqcalc(dump, "gas_direct.dat", "-w", "unit", "--norm", "self", "--method", "direct", *common)
    run([sys.executable, os.path.join(HERE, "ref_sq.py"), "--input", dump,
         "--weight", "unit", "--norm", "self", "--qmax", "6", "--nq", "120",
         "--output", path("gas_ref.dat")])
    _, s_nufft = read_table(path("gas_nufft.dat"))
    _, s_direct = read_table(path("gas_direct.dat"))
    _, s_ref = read_table(path("gas_ref.dat"))
    compare(s_direct, s_nufft, "NUFFT vs direct (unit weights)")
    compare(s_direct, s_ref, "direct vs numpy reference (unit weights)")

    q, _ = read_table(path("gas_nufft.dat"))
    high_q = q > 3.0
    mean_s = float(np.mean(s_nufft[high_q]))
    if abs(mean_s - 1.0) > 0.35:
        raise SystemExit("FAIL ideal gas: <S(q)> = %.3f at large q, expected ~1" % mean_s)
    print("  ok   ideal gas <S(q)> = %.3f for q > 3 1/A" % mean_s)

    # --- 2. weighted schemes: NUFFT vs direct -----------------------------
    print("weighted schemes")
    sqcalc(dump, "gas_neutron.dat", "-w", "neutron", "--norm", "mean", *common)
    sqcalc(dump, "gas_neutron_direct.dat", "-w", "neutron", "--norm", "mean",
           "--method", "direct", *common)
    _, s_neutron = read_table(path("gas_neutron.dat"))
    _, s_neutron_direct = read_table(path("gas_neutron_direct.dat"))
    compare(s_neutron_direct, s_neutron, "NUFFT vs direct (neutron)")

    sqcalc(dump, "gas_xray.dat", "-w", "xray", "--norm", "mean", *common)
    sqcalc(dump, "gas_xray_direct.dat", "-w", "xray", "--norm", "mean",
           "--method", "direct", *common)
    run([sys.executable, os.path.join(HERE, "ref_sq.py"), "--input", dump,
         "--weight", "xray", "--norm", "mean", "--qmax", "6", "--nq", "120",
         "--output", path("gas_xray_ref.dat")])
    _, s_xray = read_table(path("gas_xray.dat"))
    _, s_xray_direct = read_table(path("gas_xray_direct.dat"))
    _, s_xray_ref = read_table(path("gas_xray_ref.dat"))
    compare(s_xray_direct, s_xray, "NUFFT vs direct (X-ray)")
    compare(s_xray_direct, s_xray_ref, "direct vs numpy reference (X-ray)")

    # --- 3. simple cubic lattice: Bragg peaks ------------------------------
    print("simple cubic lattice")
    dump_lattice = generate("lattice.dump", natoms=125, length=20, frames=3, seed=3,
                            mode="lattice", fractions="1.0")
    run([exe, "-i", dump_lattice, "-w", "unit", "--norm", "self", "--qmax", "6",
         "--nq", "120", path("lattice.dat")])
    q_lat, s_lat = read_table(path("lattice.dat"))
    peak = np.max(s_lat)
    median = np.median(s_lat)
    if peak < 20.0 or median > 5.0:
        raise SystemExit("FAIL lattice: peak %.1f median %.2f" % (peak, median))
    q_peaks = q_lat[s_lat > 0.5 * peak]
    print("  ok   lattice peak %.1f at q = %s (median %.2f)"
          % (peak, np.array2string(q_peaks[:3], precision=2), median))

    # --- 4. reciprocal grid output ----------------------------------------
    print("reciprocal grid output")
    sqcalc(dump, "gas_nufft_grid.dat", "-w", "unit", "--norm", "self",
           "--grid", path("gas_grid.dat"), *common)
    grid = np.loadtxt(path("gas_grid.dat"))
    if grid.ndim != 2 or grid.shape[1] != 4:
        raise SystemExit("FAIL grid output has unexpected shape %s" % (grid.shape,))
    if not np.all(np.isfinite(grid)):
        raise SystemExit("FAIL grid output contains non-finite values")
    # Shell average of the grid table must reproduce the 1D table.
    qlen = np.linalg.norm(grid[:, :3], axis=1)
    shell = np.clip(((qlen - 0.0) / (6.0 / 120)).astype(int), 0, 119)
    grid_sum = np.bincount(shell, weights=grid[:, 3], minlength=120)
    grid_cnt = np.bincount(shell, minlength=120)
    shell_avg = np.divide(grid_sum, grid_cnt, out=np.zeros(120), where=grid_cnt > 0)
    _, s_from_grid = read_table(path("gas_nufft.dat"))
    # The 1D table averages |rho|^2 per mode, so compare with a loose tolerance.
    worst = np.max(np.abs(shell_avg - s_from_grid)/(1.0 + np.abs(s_from_grid)))
    if worst > 0.05:
        raise SystemExit("FAIL grid and shell tables disagree (%.3f)" % worst)
    print("  ok   grid table (%d modes) is consistent with the 1D table" % len(grid))

    # --- 5. triclinic box and unusual column order -------------------------
    print("triclinic box and column order")
    dump_tri = generate("tri.dump", natoms=150, length=18, frames=3, seed=11,
                        tilt="1.5 0.5 -2.0", columns="id xu yu zu type")
    sqcalc(dump_tri, "tri_nufft.dat", "-w", "unit", "--norm", "self", *common)
    sqcalc(dump_tri, "tri_direct.dat", "-w", "unit", "--norm", "self", "--method",
           "direct", *common)
    run([sys.executable, os.path.join(HERE, "ref_sq.py"), "--input", dump_tri,
         "--weight", "unit", "--norm", "self", "--qmax", "6", "--nq", "120",
         "--output", path("tri_ref.dat")])
    _, s_tri_nufft = read_table(path("tri_nufft.dat"))
    _, s_tri_direct = read_table(path("tri_direct.dat"))
    _, s_tri_ref = read_table(path("tri_ref.dat"))
    compare(s_tri_direct, s_tri_nufft, "NUFFT vs direct (triclinic, xu yu zu)")
    compare(s_tri_direct, s_tri_ref, "direct vs numpy reference (triclinic)")

    # --- 5b. HDF5 grid output ---------------------------------------------
    if args.h5read:
        print("hdf5 output")
        sqcalc(dump, "h5_shell.dat", "-w", "unit", "--norm", "self",
               "--grid", path("grid.h5"), *common)
        sqcalc(dump, "txt_shell.dat", "-w", "unit", "--norm", "self",
               "--grid", path("grid.txt"), *common)
        # HDF5 grid table must reproduce the text grid table exactly.
        h5_grid = run([args.h5read, path("grid.h5"), "grid"]).stdout
        with open(path("grid_from_h5.txt"), "w") as handle:
            handle.write(h5_grid)
        s_h5 = read_matrix(path("grid_from_h5.txt"))[:, 3]
        s_txt = read_matrix(path("grid.txt"))[:, 3]
        if s_h5.shape != s_txt.shape:
            raise SystemExit("FAIL hdf5 grid has %d rows, text grid has %d"
                             % (len(s_h5), len(s_txt)))
        worst = np.max(np.abs(s_h5 - s_txt)/np.maximum(np.abs(s_txt), 1.0))
        if worst > 1.0e-12:
            raise SystemExit("FAIL hdf5 vs text grid: %.3e" % worst)
        print("  ok   %-42s max deviation %.2e" % ("HDF5 vs text grid table", worst))

        # Partials stored under /shell/S_partial must match the text columns.
        sqcalc(dump, "h5_fz_shell.dat", "-w", "unit", "--norm", "self", "-fz",
               "--grid", path("grid_fz.h5"), *common)
        h5_part = run([args.h5read, path("grid_fz.h5"), "partials"]).stdout
        with open(path("part_from_h5.txt"), "w") as handle:
            handle.write(h5_part)
        part_h5 = read_matrix(path("part_from_h5.txt"))
        part_txt = read_matrix(path("h5_fz_shell.dat"))
        if part_h5.shape != part_txt[:, :part_h5.shape[1]].shape:
            raise SystemExit("FAIL hdf5 partials shape %s vs text %s"
                             % (part_h5.shape, part_txt.shape))
        worst = float(np.max(np.abs(part_h5[:, 1:] - part_txt[:, 2:2+part_h5.shape[1]-1])))
        if worst > 1.0e-12:
            raise SystemExit("FAIL hdf5 vs text partials: %.3e" % worst)
        print("  ok   %-42s max deviation %.2e" % ("HDF5 vs text partial columns", worst))

        # The shell table stored in the same file must match the text one.
        h5_shell = run([args.h5read, path("grid.h5"), "shell"]).stdout
        with open(path("shell_from_h5.txt"), "w") as handle:
            handle.write(h5_shell)
        _, s_shell_h5 = read_table(path("shell_from_h5.txt"))
        _, s_shell_txt = read_table(path("h5_shell.dat"))
        compare(s_shell_txt, s_shell_h5, "HDF5 vs text shell table", rtol=1.0e-12)
    else:
        print("hdf5 output: skipped (no HDF5 in this build)")

    # --- 6. command line errors -------------------------------------------
    print("command line handling")
    cases = [
        ([exe, "-i", dump], "missing output argument"),
        ([exe, dump], "missing input"),
        ([exe, "-i", dump, "--bogus", "1", "-"], "unknown option"),
        ([exe, "-i", dump, "-w", "neutron", "-"], "weight without mapping"),
    ]
    for cmd, label in cases:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            failures += 1
            print("  FAIL %s was accepted" % label)
        else:
            print("  ok   %s rejected" % label)
    if failures:
        raise SystemExit("%d command line checks failed" % failures)

    # --- 7. GPU backend (optional) ----------------------------------------
    if args.gpu:
        print("gpu backend")
        probe = subprocess.run([exe, "-i", dump, "--device", "gpu", "--qmax", "3",
                                "--nq", "4", "-"], capture_output=True, text=True)
        skip_reasons = ("no GPU support", "no usable CUDA device", "CUDA device",
                        "incompatible cufinufft_opts")
        if probe.returncode != 0:
            if any(marker in probe.stderr for marker in skip_reasons):
                reason = next((line for line in probe.stderr.splitlines()
                               if any(marker in line for marker in skip_reasons)),
                              probe.stderr.strip().splitlines()[0])
                print("  skip %s" % reason.strip())
            else:
                raise SystemExit("FAIL gpu probe exited %d\n%s" % (probe.returncode, probe.stderr))
        else:
            gpu_common = ["--qmax", "6", "--nq", "40", "--eps", "1e-10"]
            for weight, norm in (("unit", "self"), ("neutron", "mean"), ("xray", "mean")):
                sqcalc(dump, "cpu_%s.dat" % weight, "-w", weight, "--norm", norm, *gpu_common)
                run([exe, "-i", dump, "-m", "1:Si,2:O", "-w", weight, "--norm", norm,
                     "--device", "gpu", *gpu_common, path("gpu_%s.dat" % weight)])
                _, s_cpu = read_table(path("cpu_%s.dat" % weight))
                _, s_gpu = read_table(path("gpu_%s.dat" % weight))
                compare(s_cpu, s_gpu, "GPU vs CPU (%s weights)" % weight)

            # Single precision GPU transform against the double precision CPU one.
            run([exe, "-i", dump, "-m", "1:Si,2:O", "-w", "unit", "--norm", "self",
                 "--device", "gpu", "--precision", "single", *gpu_common,
                 path("gpu_single.dat")])
            _, s_cpu_unit = read_table(path("cpu_unit.dat"))
            _, s_gpu_single = read_table(path("gpu_single.dat"))
            compare(s_cpu_unit, s_gpu_single, "GPU float32 vs CPU float64", rtol=1.0e-4)

            # Reciprocal grid output must agree as well.
            run([exe, "-i", dump, "-m", "1:Si,2:O", "-w", "unit", "--norm", "self",
                 "--grid", path("cpu_grid.dat"), *gpu_common, path("cpu_grid_sq.dat")])
            run([exe, "-i", dump, "-m", "1:Si,2:O", "-w", "unit", "--norm", "self",
                 "--device", "gpu", "--grid", path("gpu_grid.dat"), *gpu_common,
                 path("gpu_grid_sq.dat")])
            cpu_grid = np.loadtxt(path("cpu_grid.dat"))
            gpu_grid = np.loadtxt(path("gpu_grid.dat"))
            if cpu_grid.shape != gpu_grid.shape:
                raise SystemExit("FAIL gpu grid shape %s vs %s" % (gpu_grid.shape, cpu_grid.shape))
            worst = np.max(np.abs(gpu_grid[:, 3] - cpu_grid[:, 3])
                           / np.maximum(np.abs(cpu_grid[:, 3]), 1.0))
            if worst > 1.0e-6:
                raise SystemExit("FAIL gpu grid vs cpu grid: %.3e" % worst)
            print("  ok   %-42s max deviation %.2e"
                  % ("GPU vs CPU (reciprocal grid)", worst))

    # --- 8. Debye method --------------------------------------------------
    print("debye method")
    ref_debye = os.path.join(HERE, "debye_ref.py")
    debye_q = ["--qmax", "6", "--nq", "30", "--rmax", "6", "--dr", "0.01"]
    for weight, norm in (("unit", "self"), ("neutron", "mean"), ("xray", "mean")):
        sqcalc(dump, "deb_f_%s.dat" % weight, "-w", weight, "--norm", norm,
               "--method", "debye", *debye_q)
        run([sys.executable, ref_debye, "--input", dump, "--mapping", "1:Si,2:O",
             "--weight", weight, "--norm", norm, "--output", path("deb_r_%s.dat" % weight),
             *debye_q])
        _, s_fortran = read_table(path("deb_f_%s.dat" % weight))
        _, s_ref = read_table(path("deb_r_%s.dat" % weight))
        compare(s_ref, s_fortran, "Debye vs numpy reference (%s)" % weight, rtol=1.0e-9)

    # The default cut off is half of the smallest periodic side (20 A -> 10 A).
    sqcalc(dump, "deb_default.dat", "-w", "unit", "--norm", "self", "--method", "debye",
           "--qmax", "6", "--nq", "30", "--dr", "0.01")
    sqcalc(dump, "deb_rmax10.dat", "-w", "unit", "--norm", "self", "--method", "debye",
           "--qmax", "6", "--nq", "30", "--dr", "0.01", "--rmax", "10")
    _, s_default = read_table(path("deb_default.dat"))
    _, s_rmax = read_table(path("deb_rmax10.dat"))
    compare(s_rmax, s_default, "Debye default rmax = half the box side", rtol=1.0e-12)

    # After the cut-off density correction an ideal gas must give S(q) = 1
    # everywhere (debyer agrees to machine precision on this configuration).
    q_deb, s_deb = read_table(path("deb_f_unit.dat"))
    high = q_deb > 1.5
    rms = float(np.sqrt(np.mean((s_deb[high] - 1.0)**2)))
    if rms > 0.1:
        raise SystemExit("FAIL corrected ideal gas: rms(S-1) = %.3f for q > 1.5" % rms)
    print("  ok   %-42s rms(S-1) = %.4f" % ("ideal gas S(q) -> 1 (corrected)", rms))

    # The uncorrected variant must match the reference with --no-correction.
    sqcalc(dump, "deb_raw.dat", "-w", "unit", "--norm", "n", "--method", "debye",
           "--skin", "0", "--no-cutoff-correction", *debye_q)
    run([sys.executable, ref_debye, "--input", dump, "--mapping", "1:Si,2:O",
         "--weight", "unit", "--norm", "n", "--no-correction", "--output",
         path("deb_raw_ref.dat"), *debye_q])
    _, s_raw = read_table(path("deb_raw.dat"))
    _, s_raw_ref = read_table(path("deb_raw_ref.dat"))
    compare(s_raw_ref, s_raw, "Debye uncorrected vs reference", rtol=1.0e-9)

    # Skin reuse must not change a single number.
    sqcalc(dump, "deb_skin0.dat", "-w", "unit", "--norm", "self", "--method", "debye",
           "--skin", "0", *debye_q)
    _, s_skin0 = read_table(path("deb_skin0.dat"))
    _, s_skin1 = read_table(path("deb_f_unit.dat"))
    compare(s_skin0, s_skin1, "Debye skin 0 vs skin 1", rtol=1.0e-12)

    # Pair distribution function: total g(r) against the reference, and the
    # ideal gas limit g -> 1 away from the cut off.
    sqcalc(dump, "deb_shell.dat", "-w", "unit", "--norm", "self", "--method", "debye",
           "--qmax", "4", "--nq", "8", "--rmax", "8", "--dr", "0.05",
           "--rdf", path("deb.rdf"))
    run([sys.executable, ref_debye, "--input", dump, "--mapping", "1:Si,2:O",
         "--weight", "unit", "--qmax", "4", "--nq", "8", "--rmax", "8", "--dr", "0.05",
         "--output", path("deb_r_rdf.dat"), "--rdf", path("deb_ref.rdf")])
    rdf = read_matrix(path("deb.rdf"))
    rdf_ref = read_matrix(path("deb_ref.rdf"))
    # g_ab is symmetric in a <-> b for the same reason as S_ab: the ordered
    # pair counts satisfy n_ab = n_ba and the normalizations N_a rho_b and
    # N_b rho_a are equal, so swapping -m must not change the cross column.
    run([exe, "-i", dump, "-m", "1:O,2:Si", "-w", "unit", "--norm", "self",
         "--method", "debye", "--qmax", "4", "--nq", "8", "--rmax", "8", "--dr", "0.05",
         "--rdf", path("deb_swap.rdf"), path("deb_swap_sq.dat")])
    rdf_swap = read_matrix(path("deb_swap.rdf"))
    worst_swap = float(np.max(np.abs(rdf_swap[:, 3] - rdf[:, 3])))
    if worst_swap != 0.0:
        raise SystemExit("FAIL g_ab is not symmetric in a <-> b (%.3e)" % worst_swap)
    print("  ok   %-42s cross column identical" % "partial g(r) symmetry a<->b")
    if rdf.shape[0] != rdf_ref.shape[0]:
        raise SystemExit("FAIL rdf row count %d vs %d" % (rdf.shape[0], rdf_ref.shape[0]))
    worst = np.max(np.abs(rdf[:, 1] - rdf_ref[:, 1]))
    if worst > 1.0e-12:
        raise SystemExit("FAIL total g(r) vs reference: %.3e" % worst)
    print("  ok   %-42s max deviation %.2e" % ("total g(r) vs reference", worst))
    mid = (rdf[:, 0] > 4.0) & (rdf[:, 0] < 8.0)
    mean_g = float(np.mean(rdf[mid, 1]))
    if abs(mean_g - 1.0) > 0.25:
        raise SystemExit("FAIL ideal gas g(r) = %.3f for 4 < r < 8 A" % mean_g)
    print("  ok   %-42s mean g(r) = %.3f" % ("ideal gas g(r) -> 1", mean_g))

    # The HDF5 form of the same g(r) file must carry the same numbers (the text
    # file writes 10 significant digits, hence the relative comparison).
    if args.h5read:
        sqcalc(dump, "deb_rdf_h5.dat", "-w", "unit", "--norm", "self", "--method", "debye",
               "--qmax", "4", "--nq", "8", "--rmax", "8", "--dr", "0.05",
               "--rdf", path("deb_rdf.h5"))
        with open(path("deb_rdf_h5.txt"), "w") as handle:
            handle.write(run([args.h5read, path("deb_rdf.h5"), "rdf"]).stdout)
        rdf_h5 = read_matrix(path("deb_rdf_h5.txt"))
        if rdf_h5.shape != rdf.shape:
            raise SystemExit("FAIL HDF5 rdf shape %s vs text shape %s"
                             % (rdf_h5.shape, rdf.shape))
        compare(rdf[:, 1:], rdf_h5[:, 1:], "HDF5 rdf vs text rdf", rtol=1.0e-8)

    # Lattice: the RDF peaks at the neighbour shell distances of a cubic
    # lattice with a = 4 A (4, 5.66, 6.93, 8 A).
    lattice_dump = generate("lattice_debye.dump", natoms=125, length=20, frames=1,
                            seed=3, mode="lattice", fractions="1.0")
    run([exe, "-i", lattice_dump, "-m", "1:Si", "-w", "unit", "--norm", "self",
         "--method", "debye", "--qmax", "4", "--nq", "8", "--rmax", "9", "--dr", "0.02",
         "--rdf", path("lattice.rdf"), path("lattice_sq.dat")])
    rdf_lat = read_matrix(path("lattice.rdf"))
    found = []
    for expected in (4.0, 5.657, 6.928, 8.0):
        window = np.abs(rdf_lat[:, 0] - expected) < 0.08
        if not np.any(window) or np.max(rdf_lat[window, 1]) < 1.5:
            raise SystemExit("FAIL lattice g(r) has no peak near %.3f A" % expected)
        found.append(round(float(rdf_lat[window, 0][np.argmax(rdf_lat[window, 1])]), 3))
    print("  ok   %-42s peaks at %s A" % ("lattice g(r) neighbour shells", found))

    # --- 9. partial structure factors (OVITO convention + Faber-Ziman) -----
    print("partial structure factors")
    for method, extra in (("nufft", []), ("debye", ["--rmax", "6", "--dr", "0.01", "--skin", "0"])):
        sqcalc(dump, "part_%s.dat" % method, "-w", "unit", "--norm", "n", "--method", method,
               "--qmin", "0.5", "--qmax", "5.9", "--nq", "27", *extra)
        sqcalc(dump, "partfz_%s.dat" % method, "-w", "unit", "--norm", "n", "-fz",
               "--method", method, "--qmin", "0.5", "--qmax", "5.9", "--nq", "27", *extra)
        part = read_matrix(path("part_%s.dat" % method))
        fz = read_matrix(path("partfz_%s.dat" % method))
        if part.shape[1] != 5:
            raise SystemExit("FAIL %s: expected 5 columns, got %d" % (method, part.shape[1]))
        # OVITO sum rule: S(q) = S_aa + 2 S_ab + S_bb
        worst = float(np.max(np.abs(part[:, 1] - (part[:, 2] + 2.0*part[:, 3] + part[:, 4]))))
        if worst > 0.02:
            raise SystemExit("FAIL %s: partial sum rule is off by %.3e" % (method, worst))
        print("  ok   %-42s sum rule %.1e" % ("%s partials: S = Saa + 2Sab + Sbb" % method, worst))
        # high-q limits: S_aa -> x_a = 0.5, S_ab -> 0, A_ab -> 1
        high = part[:, 0] > 3.0
        saa = float(np.mean(part[high, 2]))
        sab = float(np.mean(np.abs(part[high, 3])))
        aab = float(np.mean(fz[high, 3]))
        if abs(saa - 0.5) > 0.05 or sab > 0.1 or abs(aab - 1.0) > 0.05:
            raise SystemExit("FAIL %s: high-q partial limits Saa=%.3f Sab=%.3f Aab=%.3f"
                             % (method, saa, sab, aab))
        print("  ok   %-42s Saa=%.3f Sab=%.3f Aab=%.3f"
              % ("%s high-q partial limits" % method, saa, sab, aab))

        # S_ab is symmetric under an exchange of the two types: it is a
        # same-time product, so Re(rho_a rho_b*) and Re(rho_b rho_a*) are equal
        # exactly (they are complex conjugates), and the Debye ordered pair
        # counts satisfy n_ab = n_ba identically.  Swapping -m must therefore
        # leave the cross column bit for bit unchanged.
        run([exe, "-i", dump, "-m", "1:O,2:Si", "-w", "unit", "--norm", "n",
             "--method", method, "--qmin", "0.5", "--qmax", "5.9", "--nq", "27", *extra,
             path("part_swap_%s.dat" % method)])
        swapped = read_matrix(path("part_swap_%s.dat" % method))
        worst = float(np.max(np.abs(swapped[:, 3] - part[:, 3])))
        if worst != 0.0:
            raise SystemExit("FAIL %s: S_ab is not symmetric in a <-> b (%.3e)"
                             % (method, worst))
        print("  ok   %-42s cross column identical" % ("%s partial symmetry a<->b" % method))

    # --- 9b. partials without an element mapping ---------------------------
    # A model system (Kob-Andersen style binary mixture) has no elements, so
    # the partials have to fall back to the LAMMPS type ids.  The numbers must
    # not depend on the labels, so a mapped run has to agree exactly.
    print("partials without a mapping")
    dump_ka = generate("ka.dump", natoms=400, length=16, frames=3, seed=5,
                       fractions="0.8,0.2")
    run([exe, "-i", dump_ka, "-w", "unit", "--qmax", "6", "--nq", "120",
         path("ka_unmapped.dat")])
    run([exe, "-i", dump_ka, "-m", "1:Si,2:O", "-w", "unit", "--qmax", "6",
         "--nq", "120", path("ka_mapped.dat")])
    with open(path("ka_unmapped.dat")) as handle:
        header = handle.readline().split()
    if header != ["#", "q", "S(q)", "S(1-1)", "S(1-2)", "S(2-2)"]:
        raise SystemExit("FAIL unmapped partial columns: %s" % " ".join(header))
    ka = read_matrix(path("ka_unmapped.dat"))
    compare(read_matrix(path("ka_mapped.dat")), ka,
            "unmapped vs mapped unit weights", rtol=1.0e-12)
    # OVITO sum rule and the high-q concentration limits.
    worst = float(np.max(np.abs(ka[:, 1] - (ka[:, 2] + 2.0*ka[:, 3] + ka[:, 4]))))
    if worst > 0.02:
        raise SystemExit("FAIL unmapped partials: sum rule is off by %.3e" % worst)
    high = ka[:, 0] > 4.0
    s_aa = float(np.mean(ka[high, 2]))
    s_bb = float(np.mean(ka[high, 4]))
    if abs(s_aa - 0.8) > 0.06 or abs(s_bb - 0.2) > 0.06:
        raise SystemExit("FAIL unmapped partials: high-q S11=%.3f S22=%.3f" % (s_aa, s_bb))
    print("  ok   %-42s sum rule %.1e, high-q S11=%.3f S22=%.3f"
          % ("unmapped partials are type labelled", worst, s_aa, s_bb))
    # The Debye g(r) partials carry the same type-id labels.
    run([exe, "-i", dump_ka, "-w", "unit", "--method", "debye", "--rmax", "6",
         "--rdf", path("ka.rdf"), path("ka_sq.dat")])
    with open(path("ka.rdf")) as handle:
        rdf_header = handle.readline().split()
    if rdf_header != ["#", "r", "g(r)", "g(1-1)", "g(1-2)", "g(2-2)"]:
        raise SystemExit("FAIL unmapped partial g(r) columns: %s"
                         % " ".join(rdf_header))

    # --- 9c. q sampling helper --------------------------------------------
    # The skill ships choose_q.py, which reads the box and returns a sampling
    # with no empty shell; a shell the box cannot fill is written as 0.
    helper = os.path.join(HERE, os.pardir, "skills", "sq-calc", "scripts",
                          "choose_q.py")
    if not os.path.isfile(helper):
        print("skip q sampling helper (skills/sq-calc/scripts/choose_q.py)")
    else:
        print("q sampling helper")
        for label, source in (("cubic", dump), ("triclinic", dump_tri)):
            options = run([sys.executable, helper, source, "--qmax", "6",
                           "--print-options"]).stdout.split()
            nq = int(options[options.index("--nq") + 1])
            run([exe, "-i", source, *options, path("chosen_q.dat")])
            chosen = np.loadtxt(path("chosen_q.dat"))
            if chosen.shape[0] != nq:
                raise SystemExit("FAIL choose_q (%s): %d rows for nq = %d"
                                 % (label, chosen.shape[0], nq))
            empty = int((chosen[:, 1] == 0.0).sum())
            if empty:
                raise SystemExit("FAIL choose_q (%s): %d of %d shells empty"
                                 % (label, empty, nq))
            print("  ok   %-42s %s" % ("choose_q keeps every shell filled (%s)"
                                       % label, " ".join(options)))

        # The Debye method evaluates S(q) straight from the pair histogram, so
        # its q grid is free: a finer grid samples the same curve, and any qmin
        # is allowed (no lattice to fill).
        options = run([sys.executable, helper, dump, "--method", "debye",
                       "--qmax", "6", "--print-options"]).stdout.split()
        if options[options.index("--qmin") + 1] != "0":
            raise SystemExit("FAIL choose_q debye qmin: %s" % " ".join(options))
        if int(options[options.index("--nq") + 1]) != 600:
            raise SystemExit("FAIL choose_q debye nq: %s" % " ".join(options))
        debye_opts = ["-w", "unit", "--method", "debye", "--rmax", "6",
                      "--dr", "0.01", "--qmax", "6"]
        run([exe, "-i", dump, *debye_opts, "--nq", "40", path("deb_coarse.dat")])
        run([exe, "-i", dump, *debye_opts, "--nq", "200", path("deb_fine.dat")])
        coarse = read_matrix(path("deb_coarse.dat"))
        fine = read_matrix(path("deb_fine.dat"))
        idx = 5*np.arange(1, 41) - 3          # coincident shell centres
        if not np.allclose(coarse[:, 0], fine[idx, 0], rtol=0.0, atol=1.0e-12):
            raise SystemExit("FAIL debye q grid: shell centres do not coincide")
        worst = float(np.max(np.abs(coarse[:, 1:] - fine[idx, 1:])))
        if worst > 1.0e-8:
            raise SystemExit("FAIL debye q grid: nq 40 vs 200 differ by %.3e"
                             % worst)
        print("  ok   %-42s max deviation %.1e"
              % ("Debye S(q) is independent of --nq", worst))
        print("  ok   %-42s %s" % ("choose_q debye advice", " ".join(options)))

    # --- 9d. plot script ---------------------------------------------------
    # plot_sq.py marks the asymptote (S(q), g(r) -> 1) by default.  It needs
    # matplotlib, so the check skips itself when that is not installed.
    plotter = os.path.join(HERE, os.pardir, "skills", "sq-calc", "scripts",
                           "plot_sq.py")
    try:
        import matplotlib  # noqa: F401
        have_matplotlib = True
    except ImportError:
        have_matplotlib = False
    if not os.path.isfile(plotter) or not have_matplotlib:
        print("skip plot script (needs skills/sq-calc/scripts/plot_sq.py and "
              "matplotlib)")
    else:
        print("plot script")
        for index, extra in enumerate(([], ["--refline", "0"],
                                       ["--refline", "none"])):
            figure = path("plot_refline_%d.png" % index)
            run([sys.executable, plotter, path("gas_nufft.dat"), *extra,
                 "-o", figure])
            if not os.path.isfile(figure) or os.path.getsize(figure) == 0:
                raise SystemExit("FAIL plot script wrote no figure for %s"
                                 % (" ".join(extra) or "the default refline"))
        print("  ok   %-42s %s" % ("plot_sq.py renders", "refline 1.0 / 0 / none"))

    # --- 10. Faber-Ziman cross-check and partial g(r) ---------------------
    print("faber-ziman cross-check and partial g(r)")
    for method, extra in (("nufft", []), ("debye", ["--rmax", "6", "--dr", "0.01", "--skin", "0"])):
        sqcalc(dump, "fzchk_%s.dat" % method, "-w", "unit", "--norm", "n", "-fz",
               "--method", method, "--qmin", "0.5", "--qmax", "5.9", "--nq", "27", *extra)
        part = read_matrix(path("part_%s.dat" % method))
        fz = read_matrix(path("partfz_%s.dat" % method))
        # A_ab = (S_ab - x_a delta_ab)/(x_a x_b) + 1 with x = 0.5 for both types
        x = 0.5
        expect = np.column_stack([(part[:, 2] - x)/(x*x) + 1.0,
                                  part[:, 3]/(x*x) + 1.0,
                                  (part[:, 4] - x)/(x*x) + 1.0])
        worst = float(np.max(np.abs(fz[:, 2:5] - expect)))
        if worst > 1.0e-9:
            raise SystemExit("FAIL %s: FZ transformation is off by %.3e" % (method, worst))
        print("  ok   %-42s max deviation %.1e"
              % ("%s: A = (S - x delta)/(x x) + 1" % method, worst))
        # FZ sum rule: sum_ab (2 - delta_ab) x_a x_b A_ab = S(q)
        fsum = fz[:, 2] + 2.0*fz[:, 3]*x*x/0.25 + fz[:, 4]
        fsum = x*x*fz[:, 2] + 2.0*x*x*fz[:, 3] + x*x*fz[:, 4]
        worst = float(np.max(np.abs(fsum - fz[:, 1])))
        if worst > 0.02:
            raise SystemExit("FAIL %s: FZ sum rule is off by %.3e" % (method, worst))
        print("  ok   %-42s sum rule %.1e"
              % ("%s: S = sum (2-d) x_a x_b A_ab" % method, worst))

    # partial g(r): ideal gas -> every partial approaches 1
    rdf = read_matrix(path("deb.rdf"))
    if rdf.shape[1] != 5:
        raise SystemExit("FAIL rdf should have r, total and 3 partial columns, got %d" % rdf.shape[1])
    mid = (rdf[:, 0] > 4.0) & (rdf[:, 0] < 8.0)
    for col, label in ((1, "total"), (2, "Si-Si"), (3, "Si-O"), (4, "O-O")):
        mean_g = float(np.mean(rdf[mid, col]))
        if abs(mean_g - 1.0) > 0.25:
            raise SystemExit("FAIL ideal gas partial g(%s) = %.3f" % (label, mean_g))
    print("  ok   %-42s total %.3f Si-Si %.3f Si-O %.3f O-O %.3f"
          % ("ideal gas partial g(r) -> 1", rdf[mid, 1].mean(), rdf[mid, 2].mean(),
             rdf[mid, 3].mean(), rdf[mid, 4].mean()))

    # --- 11. dynamic structure factor along a q line (--dyn) --------------
    print("dynamic structure factor (--dyn)")
    ref_dyn = os.path.join(HERE, "ref_sqw.py")
    dyn_spec = "4,1,4,1,0,0"          # 5 q points, |q| = 1..4 along x
    dyn_q = ["--dyn", dyn_spec, "--dt", "1", "--maxframes", "20", "--lag", "1"]
    nq_dyn = 5
    dyn_dump = generate("dyn_ball.dump", natoms=1000, length=20, frames=40, seed=11,
                        mode="ballistic", temperature=1.0, columns="id type xu yu zu")
    sqcalc(dyn_dump, "dyn_sq.dat", "-w", "unit", "--norm", "mean", *dyn_q,
           "--sqw", path("dyn_sqw.dat"), "--fsq", path("dyn_fsq.dat"))
    run([sys.executable, ref_dyn, "--input", dyn_dump, "--dyn", dyn_spec,
         "--maxframes", "20", "--lag", "1", "--dt", "1", "--weight", "unit",
         "--norm", "mean", "--output-fsq", path("dyn_ref_fsq.dat"),
         "--output-sqw", path("dyn_ref_sqw.dat")])
    fsq = read_matrix(path("dyn_fsq.dat"))
    sqw = read_matrix(path("dyn_sqw.dat"))
    compare(fsq, read_matrix(path("dyn_ref_fsq.dat")), "F(q,t) vs numpy reference", rtol=1.0e-9)
    compare(sqw, read_matrix(path("dyn_ref_sqw.dat")), "S(q,w) vs numpy reference", rtol=1.0e-9)

    # The spectrum integrates to F(q,0), which is the static OUTPUT column; the
    # one-sided folding makes this exact (dw = pi/(maxframes*dt)).
    naxis_dyn = 21
    _, s_static = read_table(path("dyn_sq.dat"))
    integral = sqw[:, 4].reshape(nq_dyn, naxis_dyn).sum(axis=1) * (np.pi/20.0)
    compare(s_static, integral, "sum rule int S(q,w) dw = static S(q)", rtol=1.0e-9)
    compare(s_static, fsq[:, 4].reshape(nq_dyn, naxis_dyn)[:, 0], "F(q,0) = static S(q)",
            rtol=1.0e-9)

    # Brownian particles: F(q,t) = exp(-D q^2 t) (the amplitude carries the
    # frozen-in structure factor of the configuration, which decays away over
    # the trajectory, so only the decay rate is a clean observable).
    diff_dump = generate("dyn_diff.dump", natoms=500, length=20, frames=1000, seed=3,
                         mode="diffusive", diffusivity=0.1, columns="id type xu yu zu")
    run([exe, "-i", diff_dump, "-m", "1:Si,2:O", "-w", "unit", "--norm", "mean",
         "--dyn", "1,2,3,1,0,0", "--dt", "1", "--maxframes", "8", "--lag", "1",
         "--fsq", path("dyn_diff_fsq.dat"), path("dyn_diff_sq.dat")])
    diff_f = read_matrix(path("dyn_diff_fsq.dat"))
    naxis_d = 9
    f_d = diff_f[:, 4].reshape(2, naxis_d)
    tau_d = diff_f[:, 3].reshape(2, naxis_d)[0]
    worst = 0.0
    for k, q_value in enumerate((2.0, 3.0)):
        analytic = np.exp(-0.1 * q_value**2 * tau_d)
        valid = analytic > 0.1
        worst = max(worst, float(np.max(np.abs(f_d[k, valid] - analytic[valid]))))
    if worst > 0.06:
        raise SystemExit("FAIL diffusive F(q,t) deviates from exp(-D q2 t) by %.3f" % worst)
    print("  ok   %-42s max deviation %.3f" % ("diffusive F(q,t) vs exp(-D q2 t)", worst))

    # The same trajectory with wrapped coordinates: the reader reconstructs the
    # unwrapped positions, so the result must not change.
    dyn_wrapped = generate("dyn_ball_wrapped.dump", natoms=1000, length=20, frames=40,
                           seed=11, mode="ballistic", temperature=1.0,
                           columns="id type x y z")
    sqcalc(dyn_wrapped, "dyn_sq_wrap.dat", "-w", "unit", "--norm", "mean", *dyn_q,
           "--sqw", path("dyn_sqw_wrap.dat"))
    compare(read_matrix(path("dyn_sqw_wrap.dat")), sqw,
            "wrapped dump (reader unwraps) vs xu", rtol=1.0e-9)

    # --dt is the time step of the trajectory: the same frames written every
    # fifth step with dt/5 have to give the same spectrum.
    dyn_step = generate("dyn_ball_step.dump", natoms=1000, length=20, frames=40, seed=11,
                        mode="ballistic", temperature=1.0, columns="id type xu yu zu",
                        step=5)
    run([exe, "-i", dyn_step, "-m", "1:Si,2:O", "-w", "unit", "--norm", "mean",
         "--dyn", dyn_spec, "--dt", "0.2", "--maxframes", "20", "--lag", "1",
         "--sqw", path("dyn_sqw_step.dat"), path("dyn_sq_step.dat")])
    compare(read_matrix(path("dyn_sqw_step.dat")), sqw,
            "dt = 0.2 over 5 steps equals dt = 1 over 1", rtol=1.0e-9)

    # Partial spectra: with unit weights and --norm n the OVITO sum rule
    # S(q,w) = S_aa + 2 S_ab + S_bb holds sample by sample.
    run([exe, "-i", dyn_dump, "-m", "1:Si,2:O", "-w", "unit", "--norm", "n", *dyn_q,
         "--sqw", path("dyn_part.dat"), path("dyn_sq_n.dat")])
    part = read_matrix(path("dyn_part.dat"))
    worst = float(np.max(np.abs(part[:, 4] - (part[:, 5] + 2.0*part[:, 6] + part[:, 7]))))
    if worst > 1.0e-9:
        raise SystemExit("FAIL dynamic partial sum rule is off by %.3e" % worst)
    print("  ok   %-42s sum rule %.1e"
          % ("dynamic partials: S = Saa + 2Sab + Sbb", worst))

    # Chemical weighting: the weighted total follows the same post-processing
    # convention as the static partials (neutron b, X-ray f(q)).
    for weight, norm in (("neutron", "mean"), ("xray", "self")):
        sqcalc(dyn_dump, "dyn_w_%s.dat" % weight, "-w", weight, "--norm", norm, *dyn_q,
               "--sqw", path("dyn_w_%s_sqw.dat" % weight))
        run([sys.executable, ref_dyn, "--input", dyn_dump, "--dyn", dyn_spec,
             "--maxframes", "20", "--lag", "1", "--dt", "1", "--weight", weight,
             "--norm", norm, "--mapping", "1:Si,2:O",
             "--output-fsq", path("dyn_w_%s_ref_fsq.dat" % weight),
             "--output-sqw", path("dyn_w_%s_ref_sqw.dat" % weight)])
        compare(read_matrix(path("dyn_w_%s_sqw.dat" % weight)),
                read_matrix(path("dyn_w_%s_ref_sqw.dat" % weight)),
                "%s weighted S(q,w) vs reference" % weight, rtol=1.0e-9)

    # --lag thins the time origins; the reference applies the same rule.
    run([exe, "-i", dyn_dump, "-m", "1:Si,2:O", "-w", "unit", "--norm", "mean",
         "--dyn", dyn_spec, "--dt", "1", "--maxframes", "20", "--lag", "3",
         "--fsq", path("dyn_lag_fsq.dat"), path("dyn_lag_sq.dat")])
    run([sys.executable, ref_dyn, "--input", dyn_dump, "--dyn", dyn_spec,
         "--maxframes", "20", "--lag", "3", "--dt", "1", "--weight", "unit",
         "--norm", "mean", "--output-fsq", path("dyn_lag_ref_fsq.dat"),
         "--output-sqw", path("dyn_lag_ref_sqw.dat")])
    compare(read_matrix(path("dyn_lag_fsq.dat")), read_matrix(path("dyn_lag_ref_fsq.dat")),
            "lag stride 3 vs reference", rtol=1.0e-9)

    # A non-periodic direction must not be wrapped: the dump writes the true
    # coordinate for it (x, not xu), so an atom that leaves the box has to be
    # kept outside it.  Both spellings of the same slab trajectory must agree.
    slab_kwargs = dict(natoms=400, length=10, frames=25, seed=9, mode="ballistic",
                       temperature=0.5, boundary="pp pp ff")
    slab_xu = generate("dyn_slab_xu.dump", columns="id type xu yu zu", **slab_kwargs)
    slab_x = generate("dyn_slab_x.dump", columns="id type x y z", **slab_kwargs)
    slab_q = ["--dyn", "2,1,3,0,0,1", "--dt", "1", "--maxframes", "5"]
    for dump_file, tag in ((slab_xu, "xu"), (slab_x, "x")):
        run([exe, "-i", dump_file, "-m", "1:Si,2:O", "-w", "unit", "--norm", "mean",
             *slab_q, "--fsq", path("dyn_slab_%s_fsq.dat" % tag),
             path("dyn_slab_%s_sq.dat" % tag)])
    compare(read_matrix(path("dyn_slab_x_fsq.dat")),
            read_matrix(path("dyn_slab_xu_fsq.dat")),
            "non-periodic axis: x y z (unwrapped) vs xu", rtol=1.0e-9)

    plot_dyn = os.path.join(HERE, os.pardir, "skills", "sq-calc", "scripts", "plot_sqw.py")
    if os.path.isfile(plot_dyn) and have_matplotlib:
        for mode in ("both", "map", "spectra"):
            figure = path("plot_sqw_%s.png" % mode)
            run([sys.executable, plot_dyn, "--mode", mode, path("dyn_sqw.dat"),
                 "-o", figure])
            if not os.path.isfile(figure) or os.path.getsize(figure) == 0:
                raise SystemExit("FAIL plot_sqw.py wrote no figure in mode %s" % mode)
        print("  ok   %-42s map / spectra / both" % "plot_sqw.py renders")

    if args.h5read:
        run([exe, "-i", dyn_dump, "-m", "1:Si,2:O", "-w", "unit", "--norm", "mean", *dyn_q,
             "--sqw", path("dyn_sqw.h5"), "--fsq", path("dyn_fsq.h5"), path("dyn_sq_h5.dat")])
        with open(path("dyn_sqw_h5.txt"), "w") as handle:
            handle.write(run([args.h5read, path("dyn_sqw.h5"), "sqw"]).stdout)
        with open(path("dyn_fsq_h5.txt"), "w") as handle:
            handle.write(run([args.h5read, path("dyn_fsq.h5"), "fsq"]).stdout)
        compare(read_matrix(path("dyn_sqw_h5.txt")), sqw, "HDF5 S(q,w) vs text", rtol=1.0e-9)
        compare(read_matrix(path("dyn_fsq_h5.txt")), fsq, "HDF5 F(q,t) vs text", rtol=1.0e-9)

    # option validation
    bad = [
        (["--dyn", dyn_spec, "--maxframes", "20"], "--dyn without --dt"),
        (["--dyn", dyn_spec, "--dt", "1"], "--dyn without --maxframes"),
        (["--dyn", dyn_spec, "--dt", "1", "--maxframes", "20", "--grid", "g.dat"],
         "--dyn with --grid"),
        (["--sqw", "x.dat"], "--sqw without --dyn"),
        (["--dyn", "0,1,4,1,0,0", "--dt", "1", "--maxframes", "20"], "--dyn without intervals"),
        (["--dyn", "4,4,1,1,0,0", "--dt", "1", "--maxframes", "20"], "--dyn with S1 <= S0"),
    ]
    rejected = 0
    for options, label in bad:
        result = subprocess.run([exe, "-i", dyn_dump] + options + ["-"],
                                capture_output=True, text=True)
        if result.returncode == 0:
            raise SystemExit("FAIL %s was accepted" % label)
        rejected += 1
    print("  ok   %-42s %d cases rejected" % ("--dyn option validation", rejected))

    # A trajectory that changes its atom count or its atom order cannot be
    # unwrapped index by index: the reader has to say so instead of running
    # past the end of its per atom arrays.
    for name, columns, label in (("grow", "id type x y z", "atom count changes"),
                                 ("reorder", "id type x y z", "atom order changes")):
        broken = path("dyn_%s.dump" % name)
        with open(broken, "w") as handle:
            for frame in (0, 1, 2):
                natoms = 20 if (name == "grow" and frame == 0) else 40
                handle.write("ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\n" % (frame, natoms))
                handle.write("ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n")
                handle.write("ITEM: ATOMS %s\n" % columns)
                order = range(natoms, 0, -1) if name == "reorder" and frame > 0 \
                    else range(1, natoms + 1)
                for index in order:
                    handle.write("%d 1 %.6f %.6f %.6f\n" % (index, (index % 10)*1.0,
                                                            (index % 7)*1.3, (index % 5)*1.7))
        result = subprocess.run([exe, "-i", broken, "-w", "unit", "--dyn", dyn_spec,
                                 "--dt", "1", "--maxframes", "20", "-"],
                                capture_output=True, text=True)
        if result.returncode == 0:
            raise SystemExit("FAIL %s was accepted" % label)
        print("  ok   %-42s rejected" % label)

    print("all sqcalc tests passed")


if __name__ == "__main__":
    main()
