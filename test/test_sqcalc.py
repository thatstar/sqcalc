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
import math
import os
import subprocess
import sys
from itertools import permutations, product

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def run(cmd, env=None):
    if env is not None:
        env = dict(os.environ, **env)
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
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


def dyn_line(spec):
    """The q sampling flags of a q line through the Gamma point.

    This uses the short alias, so every dynamic test exercises `-q` as well;
    the long `--qpoints` spelling is checked against it in the validation
    section at the end.
    """
    return ["-q", "line:" + spec]


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
        cmd = [exe, "static", "-i", dump, "-m", "1:Si,2:O", *options, path(out)]
        return run(cmd)

    def sqcalc_dyn(dump, *options):
        cmd = [exe, "dyn", "-i", dump, *options]
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

    # --- 2b. ion species: an ion label selects a different X-ray row -------
    print("ion species")
    run([exe, "static", "-i", dump, "-m", "1:Si4+,2:O2-", "-w", "xray", "--norm", "mean",
         *common, path("gas_ion.dat")])
    run([sys.executable, os.path.join(HERE, "ref_sq.py"), "--input", dump,
         "--weight", "xray", "--norm", "mean", "--mapping", "1:Si4+,2:O2-",
         "--qmax", "6", "--nq", "120", "--output", path("gas_ion_ref.dat")])
    _, s_ion = read_table(path("gas_ion.dat"))
    _, s_ion_ref = read_table(path("gas_ion_ref.dat"))
    compare(s_ion, s_ion_ref, "ionic X-ray vs numpy reference")
    q_ion, _ = read_table(path("gas_ion.dat"))
    low = q_ion < 1.0
    shift = np.max(np.abs(s_ion[low] - s_xray[low]) / (1.0 + np.abs(s_xray[low])))
    if shift < 1.0e-3:
        raise SystemExit("FAIL ion species: the charge state did not change S(q)")
    print("  ok   Si4+/O2- moves S(q) by %.3f at q < 1 1/A" % shift)

    # The neutron length is nuclear: an ion label must not change it.
    run([exe, "static", "-i", dump, "-m", "1:Si4+,2:O2-", "-w", "neutron", "--norm",
         "mean", *common, path("gas_ion_neutron.dat")])
    _, s_ion_neutron = read_table(path("gas_ion_neutron.dat"))
    compare(s_ion_neutron, s_neutron, "neutron ignores the charge state")

    # The pair columns carry the species label, and an unknown state is named.
    header = open(path("gas_ion.dat")).readline()
    if "S(Si4+-O2-)" not in header:
        raise SystemExit("FAIL ion species: pair label missing from %r" % header)
    print("  ok   pair columns are labelled with the species")
    bad = subprocess.run([exe, "static", "-i", dump, "-m", "1:O3-", "-w", "xray", "-"],
                         capture_output=True, text=True)
    if bad.returncode == 0 or "available: O, O1-, O2-" not in bad.stderr:
        raise SystemExit("FAIL ion species: O3- was accepted\n%s" % bad.stderr)
    print("  ok   an unknown charge state lists the available ones")

    # --- 3. simple cubic lattice: Bragg peaks ------------------------------
    print("simple cubic lattice")
    dump_lattice = generate("lattice.dump", natoms=125, length=20, frames=3, seed=3,
                            mode="lattice", fractions="1.0")
    run([exe, "static", "-i", dump_lattice, "-w", "unit", "--norm", "self", "--qmax", "6",
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
        ([exe, "static", "-i", dump], "missing output argument"),
        ([exe, "static", dump], "missing input"),
        ([exe, "static", "-i", dump, "--bogus", "1", "-"], "unknown option"),
        ([exe, "static", "-i", dump, "-w", "neutron", "-"], "weight without mapping"),
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
        probe = subprocess.run([exe, "static", "-i", dump, "--device", "gpu", "--qmax", "3",
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
                run([exe, "static", "-i", dump, "-m", "1:Si,2:O", "-w", weight, "--norm", norm,
                     "--device", "gpu", *gpu_common, path("gpu_%s.dat" % weight)])
                _, s_cpu = read_table(path("cpu_%s.dat" % weight))
                _, s_gpu = read_table(path("gpu_%s.dat" % weight))
                compare(s_cpu, s_gpu, "GPU vs CPU (%s weights)" % weight)

            # Single precision GPU transform against the double precision CPU one.
            run([exe, "static", "-i", dump, "-m", "1:Si,2:O", "-w", "unit", "--norm", "self",
                 "--device", "gpu", "--precision", "single", *gpu_common,
                 path("gpu_single.dat")])
            _, s_cpu_unit = read_table(path("cpu_unit.dat"))
            _, s_gpu_single = read_table(path("gpu_single.dat"))
            compare(s_cpu_unit, s_gpu_single, "GPU float32 vs CPU float64", rtol=1.0e-4)

            # Reciprocal grid output must agree as well.
            run([exe, "static", "-i", dump, "-m", "1:Si,2:O", "-w", "unit", "--norm", "self",
                 "--grid", path("cpu_grid.dat"), *gpu_common, path("cpu_grid_sq.dat")])
            run([exe, "static", "-i", dump, "-m", "1:Si,2:O", "-w", "unit", "--norm", "self",
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

            # Pair columns come from the per species amplitudes on the host
            # side of the download, so the S(q) partials have to match the CPU
            # run column by column.
            run([exe, "static", "-i", dump, "-m", "1:Si,2:O", "-w", "xray", *gpu_common,
                 path("cpu_p.dat")])
            run([exe, "static", "-i", dump, "-m", "1:Si,2:O", "-w", "xray", "--device", "gpu",
                 *gpu_common, path("gpu_p.dat")])
            cpu_part = read_matrix(path("cpu_p.dat"))
            gpu_part = read_matrix(path("gpu_p.dat"))
            if cpu_part.shape != gpu_part.shape:
                raise SystemExit("FAIL gpu partials shape %s vs %s"
                                 % (gpu_part.shape, cpu_part.shape))
            worst = float(np.max(np.abs(gpu_part[:, 1:] - cpu_part[:, 1:])
                                 / np.maximum(np.abs(cpu_part[:, 1:]), 1.0)))
            if worst > 1.0e-6:
                raise SystemExit("FAIL gpu vs cpu partials: %.3e" % worst)
            print("  ok   %-42s max deviation %.2e" % ("GPU vs CPU (S(q) partials)", worst))

            # The powder XRD pattern accumulates through the shared per-mode
            # reduction and carries the pair columns as well, so every column
            # has to match.
            xrd_gpu = ["--xrd-lambda", "1.5406", "--xrd-range", "10", "60",
                       "--xrd-step", "1"]
            run([exe, "static", "-i", dump, "-m", "1:Si,2:O", "-w", "xray", "--xrd",
                 path("cpu.xrd"), *xrd_gpu, *gpu_common, path("cpu_xrd_sq.dat")])
            run([exe, "static", "-i", dump, "-m", "1:Si,2:O", "-w", "xray", "--device", "gpu",
                 "--xrd", path("gpu.xrd"), *xrd_gpu, *gpu_common, path("gpu_xrd_sq.dat")])
            cpu_xrd = read_matrix(path("cpu.xrd"))
            gpu_xrd = read_matrix(path("gpu.xrd"))
            if cpu_xrd.shape != gpu_xrd.shape:
                raise SystemExit("FAIL gpu xrd shape %s vs %s"
                                 % (gpu_xrd.shape, cpu_xrd.shape))
            worst = float(np.max(np.abs(gpu_xrd[:, 1:] - cpu_xrd[:, 1:])
                                 / np.maximum(np.abs(cpu_xrd[:, 1:]), 1.0)))
            if worst > 1.0e-6:
                raise SystemExit("FAIL gpu vs cpu powder XRD: %.3e" % worst)
            print("  ok   %-42s max deviation %.2e"
                  % ("GPU vs CPU (powder XRD, pair columns)", worst))

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
    run([exe, "static", "-i", dump, "-m", "1:O,2:Si", "-w", "unit", "--norm", "self",
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

        # --format has to override the file name in both directions: HDF5 into
        # a plain name, text into a .h5 name.
        sqcalc(dump, "deb_rdf_fmt_sq.dat", "-w", "unit", "--norm", "self", "--method", "debye",
               "--qmax", "4", "--nq", "8", "--rmax", "8", "--dr", "0.05",
               "--format", "hdf5", "--rdf", path("deb_rdf_fmt.dat"))
        with open(path("deb_rdf_fmt.txt"), "w") as handle:
            handle.write(run([args.h5read, path("deb_rdf_fmt.dat"), "rdf"]).stdout)
        rdf_fmt = read_matrix(path("deb_rdf_fmt.txt"))
        if rdf_fmt.shape != rdf.shape:
            raise SystemExit("FAIL --format hdf5 rdf shape %s vs text shape %s"
                             % (rdf_fmt.shape, rdf.shape))
        compare(rdf[:, 1:], rdf_fmt[:, 1:],
                "rdf: --format hdf5 beats a .dat name", rtol=1.0e-8)

        sqcalc(dump, "deb_rdf_forced_sq.dat", "-w", "unit", "--norm", "self", "--method", "debye",
               "--qmax", "4", "--nq", "8", "--rmax", "8", "--dr", "0.05",
               "--format", "text", "--rdf", path("deb_rdf_forced.h5"))
        rdf_txt = read_matrix(path("deb_rdf_forced.h5"))
        compare(rdf[:, 1:], rdf_txt[:, 1:],
                "rdf: --format text beats a .h5 name", rtol=1.0e-8)

    # --- pair entropy from the Debye g(r) ---------------------------------
    print("pair entropy (--pair-entropy/--s2-accum)")
    s2_dump = generate("s2_gas.dump", natoms=1000, length=25, frames=200, seed=13)
    run([exe, "static", "-i", s2_dump, "-m", "1:Si,2:O", "-w", "unit", "--method", "debye",
         "--rmax", "8", "--dr", "0.05", "--rdf", path("s2.rdf"),
         "--pair-entropy", path("s2.dat"), "--s2-accum", path("s2acc.dat"),
         path("s2_sq.dat")])
    rdf = read_matrix(path("s2.rdf"))
    accum = read_matrix(path("s2acc.dat"))
    meta = {}
    with open(path("s2acc.dat")) as handle:
        for line in handle:
            if line.startswith("# counts"):
                meta["counts"] = [int(v) for v in line.split()[2:]]
            elif line.startswith("# input"):
                tokens = line.split()
                for i, token in enumerate(tokens):
                    if token in ("natoms", "nframes"):
                        meta[token] = int(tokens[i + 1])
                    if token in ("volume", "rmax", "dr"):
                        meta[token] = float(tokens[i + 1])
    counts = meta["counts"]
    natoms = float(meta["natoms"])
    volume = float(meta["volume"])
    nframes = float(meta["nframes"])
    dr = float(meta["dr"])
    npairs = rdf.shape[1] - 2
    ntypes = int((math.isqrt(1 + 8*npairs) - 1)//2)
    pairs = [(a, b) for a in range(ntypes) for b in range(a, ntypes)]
    r = rdf[:, 0]
    rho = natoms/volume
    dv = 4.0*math.pi/3.0*((r + 0.5*dr)**3 - np.maximum(0.0, r - 0.5*dr)**3)
    curves = np.zeros((len(r), npairs))
    for p, (a, b) in enumerate(pairs):
        if (a == b and counts[a] <= 1) or (a != b and (counts[a] == 0 or counts[b] == 0)):
            continue
        # The --rdf table uses the historical midpoint shell for g(r); the
        # pair-entropy integral uses the exact shell volume, so convert first.
        g = rdf[:, 2 + p]*(4.0*math.pi*r**2*dr)/dv
        integrand = np.where(g > 0.0,
                             g*np.log(np.maximum(g, 1.0e-300)) - g + 1.0, 1.0)
        contrib = dv*integrand
        bias = volume/(2.0*nframes*counts[a]*counts[b])
        curves[:, p] = -2.0*math.pi*rho*(counts[a]/natoms)*(counts[b]/natoms) \
            * np.cumsum(contrib - bias)
    total_curve = np.zeros(len(r))
    for p, (a, b) in enumerate(pairs):
        total_curve += (1.0 if a == b else 2.0)*curves[:, p]
    compare(accum[:, 1].reshape(-1, 1), total_curve.reshape(-1, 1),
            "S2 total curve vs reference", rtol=1.0e-9)
    compare(accum[:, 2:], curves, "S2 partial curves vs reference", rtol=1.0e-9)
    text_final = []
    with open(path("s2.dat")) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            text_final.append(float(line.split()[1]))
    text_final = np.array(text_final)
    if abs(text_final[0] - total_curve[-1]) > 1.0e-12:
        raise SystemExit("FAIL final S2 does not match the accumulation curve")
    if abs(text_final[0]) > 0.1:
        raise SystemExit("FAIL ideal gas S2 = %.3f is too large" % text_final[0])
    print("  ok   %-42s total %.4f, partials %d" %
          ("ideal gas S2 and partial curves", text_final[0], npairs))

    if args.h5read:
        run([exe, "static", "-i", s2_dump, "-m", "1:Si,2:O", "-w", "unit", "--method", "debye",
             "--rmax", "8", "--dr", "0.05", "--pair-entropy", path("s2.h5"),
             "--s2-accum", path("s2acc.h5"), path("s2_h5_sq.dat")])
        with open(path("s2_h5.txt"), "w") as handle:
            handle.write(run([args.h5read, path("s2.h5"), "pair_entropy"]).stdout)
        with open(path("s2acc_h5.txt"), "w") as handle:
            handle.write(run([args.h5read, path("s2acc.h5"), "s2_accum"]).stdout)
        compare(read_matrix(path("s2acc_h5.txt")), accum, "HDF5 S2 curve vs text", rtol=1.0e-8)
        compare(read_matrix(path("s2_h5.txt")).reshape(-1),
                text_final.reshape(-1), "HDF5 pair entropy vs text", rtol=1.0e-8)

        # --format decides the container of both entropy files, whatever the
        # names look like.
        run([exe, "static", "-i", s2_dump, "-m", "1:Si,2:O", "-w", "unit", "--method", "debye",
             "--rmax", "8", "--dr", "0.05", "--format", "hdf5",
             "--pair-entropy", path("s2_fmt.dat"),
             "--s2-accum", path("s2acc_fmt.dat"), path("s2_fmt_sq.dat")])
        with open(path("s2_fmt.txt"), "w") as handle:
            handle.write(run([args.h5read, path("s2_fmt.dat"), "pair_entropy"]).stdout)
        with open(path("s2acc_fmt.txt"), "w") as handle:
            handle.write(run([args.h5read, path("s2acc_fmt.dat"), "s2_accum"]).stdout)
        compare(read_matrix(path("s2acc_fmt.txt")), accum,
                "pair entropy: --format hdf5 beats a .dat name", rtol=1.0e-8)
        compare(read_matrix(path("s2_fmt.txt")).reshape(-1), text_final.reshape(-1),
                "S2 curve: --format hdf5 beats a .dat name", rtol=1.0e-8)

    # The optional Python analysis tool must run on the generated files.
    try:
        import scipy  # noqa: F401
        have_scipy = True
    except ImportError:
        have_scipy = False
    if have_scipy:
        tool = os.path.join(HERE, os.pardir, "skills", "sq-calc", "scripts", "s2_analysis.py")
        run([sys.executable, tool, "--s2", path("s2.dat"), "--rdf", path("s2.rdf"),
             "--accum", path("s2acc.dat"), "--s2-smooth", "--fit", "power",
             "--r-fit-min", "4", "--r-fit-max", "8", "--json", path("s2_tool.json")])
        import json
        with open(path("s2_tool.json")) as handle:
            data = json.load(handle)
        if "tail_fit" not in data or "gcv_smooth" not in data:
            raise SystemExit("FAIL s2_analysis.py JSON is missing results")
        print("  ok   %-42s tail fit + GCV" % "s2_analysis.py")

    # Lattice: the RDF peaks at the neighbour shell distances of a cubic
    # lattice with a = 4 A (4, 5.66, 6.93, 8 A).
    lattice_dump = generate("lattice_debye.dump", natoms=125, length=20, frames=1,
                            seed=3, mode="lattice", fractions="1.0")
    run([exe, "static", "-i", lattice_dump, "-m", "1:Si", "-w", "unit", "--norm", "self",
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
    for method, extra in (("nufft", []), ("direct", []),
                          ("debye", ["--rmax", "6", "--dr", "0.01", "--skin", "0"])):
        sqcalc(dump, "part_%s.dat" % method, "-w", "unit", "--norm", "n", "--method", method,
               "--qmin", "0.5", "--qmax", "5.9", "--nq", "27", *extra)
        sqcalc(dump, "partfz_%s.dat" % method, "-w", "unit", "--norm", "n", "-fz",
               "--method", method, "--qmin", "0.5", "--qmax", "5.9", "--nq", "27", *extra)
        part = read_matrix(path("part_%s.dat" % method))
        fz = read_matrix(path("partfz_%s.dat" % method))
        if part.shape[1] != 5:
            raise SystemExit("FAIL %s: expected 5 columns, got %d" % (method, part.shape[1]))
        # OVITO sum rule: S(q) = S_aa + 2 S_ab + S_bb.  This is an identity of
        # the accumulation (see "partials are accumulated exactly" below), so
        # it holds to roundoff for both methods, not to a few percent.
        worst = float(np.max(np.abs(part[:, 1] - (part[:, 2] + 2.0*part[:, 3] + part[:, 4]))))
        if worst > 1.0e-9:
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
        # leave the cross column unchanged; the comparison is to roundoff
        # rather than exactly, because the two runs reach the amplitudes
        # through a different order of operations inside the transform.
        run([exe, "static", "-i", dump, "-m", "1:O,2:Si", "-w", "unit", "--norm", "n",
             "--method", method, "--qmin", "0.5", "--qmax", "5.9", "--nq", "27", *extra,
             path("part_swap_%s.dat" % method)])
        swapped = read_matrix(path("part_swap_%s.dat" % method))
        worst = float(np.max(np.abs(swapped[:, 3] - part[:, 3])))
        if worst > 1.0e-12:
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
    run([exe, "static", "-i", dump_ka, "-w", "unit", "--qmax", "6", "--nq", "120",
         path("ka_unmapped.dat")])
    run([exe, "static", "-i", dump_ka, "-m", "1:Si,2:O", "-w", "unit", "--qmax", "6",
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
    if worst > 1.0e-9:
        raise SystemExit("FAIL unmapped partials: sum rule is off by %.3e" % worst)
    high = ka[:, 0] > 4.0
    s_aa = float(np.mean(ka[high, 2]))
    s_bb = float(np.mean(ka[high, 4]))
    if abs(s_aa - 0.8) > 0.06 or abs(s_bb - 0.2) > 0.06:
        raise SystemExit("FAIL unmapped partials: high-q S11=%.3f S22=%.3f" % (s_aa, s_bb))
    print("  ok   %-42s sum rule %.1e, high-q S11=%.3f S22=%.3f"
          % ("unmapped partials are type labelled", worst, s_aa, s_bb))
    # The Debye g(r) partials carry the same type-id labels.
    run([exe, "static", "-i", dump_ka, "-w", "unit", "--method", "debye", "--rmax", "6",
         "--rdf", path("ka.rdf"), path("ka_sq.dat")])
    with open(path("ka.rdf")) as handle:
        rdf_header = handle.readline().split()
    if rdf_header != ["#", "r", "g(r)", "g(1-1)", "g(1-2)", "g(2-2)"]:
        raise SystemExit("FAIL unmapped partial g(r) columns: %s"
                         % " ".join(rdf_header))

    # A sparse type id: an analysed group can leave one type without atoms, so
    # the pair sums have to be indexed by type id, not by the position of the
    # species in the list of types that are present.  Rewriting the type column
    # of the two component gas (2 -> 3) has to leave every number unchanged.
    with open(dump) as handle:
        lines = handle.readlines()
    rewritten, in_atoms = [], False
    for line in lines:
        if line.startswith("ITEM: ATOMS"):
            in_atoms = True
        elif line.startswith("ITEM"):
            in_atoms = False
        elif in_atoms and line.strip():
            fields = line.split()
            fields[1] = "3" if fields[1] == "2" else fields[1]
            line = " ".join(fields) + "\n"
        rewritten.append(line)
    with open(path("sparse.dump"), "w") as handle:
        handle.writelines(rewritten)
    sqcalc(dump, "dense_ids.dat", "-w", "unit", "--qmax", "6", "--nq", "30")
    run([exe, "static", "-i", path("sparse.dump"), "-m", "1:Si,3:O", "-w", "unit",
         "--qmax", "6", "--nq", "30", path("sparse_ids.dat")])
    dense = read_matrix(path("dense_ids.dat"))
    sparse = read_matrix(path("sparse_ids.dat"))
    if sparse.shape != dense.shape:
        raise SystemExit("FAIL sparse type ids: %s columns vs %s"
                         % (sparse.shape, dense.shape))
    worst = float(np.max(np.abs(sparse - dense)))
    # A tolerance rather than bit identity: the two runs number the types
    # differently, so the transform sees the same atoms in a different order
    # and the last bit of the amplitudes can differ.
    if worst > 1.0e-12:
        raise SystemExit("FAIL sparse type ids: columns differ by %.3e" % worst)
    if float(np.max(np.abs(sparse[:, 3]))) <= 0.0:
        raise SystemExit("FAIL sparse type ids: the cross column is empty")
    print("  ok   %-42s columns identical, cross term %.3f"
          % ("sparse type ids (1 and 3 of 3)", float(np.max(np.abs(sparse[:, 3])))))

    # --- 9c. partials are accumulated exactly ------------------------------
    # The partial sums are real numbers; they used to be accumulated in an
    # integer array, which truncated every shell sum to a whole number and
    # overflowed int32 once a shell held sharp Bragg peaks (|rho|^2 ~ N^2 per
    # mode, so a few 1e4 atoms already exceed 2^31).  Both effects show up
    # against quantities that are exact without a reference: for one species
    # S(q) and S(a-a) are the same number, and a commensurate lattice has
    # S(q) = N in a shell that holds only the Bragg family.
    print("partials are accumulated exactly")
    gas1 = generate("part_gas1.dump", natoms=400, length=20, frames=3, seed=5,
                    fractions=1.0)
    run([exe, "static", "-i", gas1, "-w", "unit", "--norm", "n", "--qmax", "6",
         "--nq", "60", path("part_gas1.dat")])
    one = read_matrix(path("part_gas1.dat"))
    if one.shape[1] != 3:
        raise SystemExit("FAIL partial accumulator: %d columns, expected 3"
                         % one.shape[1])
    worst = float(np.max(np.abs(one[:, 1] - one[:, 2])
                         / np.maximum(np.abs(one[:, 1]), 1.0)))
    if worst > 1.0e-6:
        raise SystemExit("FAIL partial accumulator: S(q) and S(a-a) differ by "
                         "%.3e for a single species (truncated?)" % worst)
    print("  ok   %-42s max deviation %.1e"
          % ("single species S(q) = S(a-a)", worst))

    # The cross term uses the same accumulator, so the OVITO sum rule has to
    # hold to roundoff rather than to the few percent the older check allowed.
    run([exe, "static", "-i", dump, "-m", "1:Si,2:O", "-w", "unit", "--norm", "n",
         "--qmax", "6", "--nq", "60", path("part_sum_rule.dat")])
    two = read_matrix(path("part_sum_rule.dat"))
    worst = float(np.max(np.abs(two[:, 1] - (two[:, 2] + 2.0*two[:, 3] + two[:, 4]))))
    if worst > 1.0e-6:
        raise SystemExit("FAIL partial accumulator: sum rule off by %.3e" % worst)
    print("  ok   %-42s off by %.1e" % ("S = Saa + 2Sab + Sbb", worst))

    # Sharp Bragg peaks: |rho|^2 = N^2 per mode, so summing the six <100> modes
    # of a 32768-atom lattice needs 6 N^2 ~ 6.4e9 and used to come back as a
    # fraction of a single mode (or a negative number).
    lat = generate("part_lat.dump", natoms=32768, length=20, frames=1,
                   mode="lattice", fractions=1.0, seed=3)
    run([exe, "static", "-i", lat, "-w", "unit", "--norm", "n", "--qmin", "10.05",
         "--qmax", "10.056", "--nq", "1", "--eps", "1e-10",
         path("part_lat.dat")])
    bragg = read_matrix(path("part_lat.dat"))
    peak = float(bragg[0, 1])
    if abs(peak - 32768.0) / 32768.0 > 1.0e-6:
        raise SystemExit("FAIL partial accumulator: Bragg shell S = %.6g, "
                         "expected N = 32768" % peak)
    worst = abs(bragg[0, 1] - bragg[0, 2]) / abs(bragg[0, 1])
    if worst > 1.0e-9:
        raise SystemExit("FAIL partial accumulator: Bragg shell S(q) and S(a-a) "
                         "differ by %.3e (integer overflow?)" % worst)
    print("  ok   %-42s S(q) = S(a-a) = N = %.6g"
          % ("Bragg shell, 32768 atoms", peak))

    # --- 9d. powder XRD pattern -------------------------------------------
    # LAMMPS's own example geometry: fcc Ni, a = 3.52 A, 20x20x20 cells, Cu
    # Kalpha.  The wavelength, range and bin width are the ones of
    # examples/PACKAGES/diffraction, so the pattern can be compared with the
    # reference histogram LAMMPS ships.  The intensity is per atom, as theirs
    # is (they divide the same sum by N), which is what makes the numbers
    # comparable at all.
    print("powder XRD pattern")
    ni = generate("xrd_ni.dump", natoms=32000, length=70.4, frames=1, mode="fcc",
                  fractions=1.0, seed=3)
    ni_opts = ["--xrd-lambda", "1.541838", "--xrd-range", "40", "80"]
    run([exe, "static", "-i", ni, "-m", "1:Ni", "-w", "xray", "--xrd", path("ni.xrd"),
         *ni_opts, "--xrd-step", "0.2"])
    xrd = read_matrix(path("ni.xrd"))
    if xrd.shape != (200, 3):
        raise SystemExit("FAIL xrd: got a %s table for 40-80 deg at 0.2 deg"
                         % (xrd.shape,))
    with open(path("ni.xrd")) as handle:
        header = [handle.readline() for _ in range(4)]
    if not all(line.startswith("#") for line in header) or "lambda" not in header[0]:
        raise SystemExit("FAIL xrd: the header does not describe the run")
    if "Lorentz-polarization" not in header[1]:
        raise SystemExit("FAIL xrd: the header does not mention the LP factor")
    if header[3].split() != ["#", "2theta[deg]", "I", "I(Ni-Ni)"]:
        raise SystemExit("FAIL xrd: unexpected columns %r" % header[3])
    # one species: the pair column is the total, since every pair is Ni-Ni
    worst = float(np.max(np.abs(xrd[:, 1] - xrd[:, 2])/np.maximum(np.abs(xrd[:, 1]), 1.0)))
    if worst > 1.0e-9:
        raise SystemExit("FAIL xrd: I(Ni-Ni) differs from the total by %.3e" % worst)

    # IT92 form factor of Ni, the row src/sqc_element_data.f90 carries.
    ni_a = [12.8376, 7.2920, 4.4438, 2.3800]
    ni_b = [3.8785, 0.2565, 12.1763, 66.3421]
    ni_c = 1.0341

    def form_factor(q):
        return ni_c + sum(ai*np.exp(-bi*(q/(4.0*np.pi))**2)
                          for ai, bi in zip(ni_a, ni_b))

    def lorentz_polarization(tth):
        theta = np.radians(tth)/2.0
        return ((1.0 + np.cos(np.radians(tth))**2)/(np.sin(theta)**2*np.cos(theta)))

    # fcc Ni: all-even or all-odd hkl, with the multiplicity of the cubic
    # lattice, so the relative intensities are m |F|^2 LP.
    lines = [("111", 1, 1, 1, 8), ("200", 2, 0, 0, 6), ("220", 2, 2, 0, 12)]
    expect = []
    for label, h, k, l, mult in lines:
        q = np.sqrt(h*h + k*k + l*l)*2.0*np.pi/3.52
        tth = 2.0*np.degrees(np.arcsin(q*1.541838/(4.0*np.pi)))
        expect.append((label, tth, mult*form_factor(q)**2*lorentz_polarization(tth)))
    peak = float(xrd[:, 1].max())
    ratios = []
    for label, tth, want in expect:
        j = int(np.argmin(np.abs(xrd[:, 0] - tth)))
        if abs(xrd[j, 0] - tth) > 0.15:
            raise SystemExit("FAIL xrd: no bin near %.2f deg for Ni(%s)" % (tth, label))
        got = xrd[j, 1]/peak
        rel = want/expect[0][2]
        ratios.append(got)
        if abs(got - rel) > 0.01*rel:
            raise SystemExit("FAIL xrd: Ni(%s) at %.2f deg is %.4f of the peak, "
                             "m|F|^2 LP predicts %.4f" % (label, tth, got, rel))
    # a perfect lattice puts every line in one bin and leaves the rest empty
    fourth = float(np.sort(xrd[:, 1])[-4])
    if fourth > 1.0e-6*peak:
        raise SystemExit("FAIL xrd: a fourth line holds %.3e of the peak" % (fourth/peak))
    print("  ok   %-42s 1 : %.4f : %.4f" % ("fcc Ni, Cu Ka, 40-80 deg",
                                            ratios[1], ratios[2]))

    # Two species: the pair columns have to add up to the total, which is what
    # makes them a decomposition rather than three unrelated curves.
    run([exe, "static", "-i", dump, "-m", "1:Si,2:O", "-w", "xray", "--xrd", path("mix.xrd"),
         "--xrd-lambda", "1.5418", "--xrd-range", "5", "60", "--xrd-step", "0.5"])
    mix = read_matrix(path("mix.xrd"))
    if mix.shape[1] != 5:
        raise SystemExit("FAIL xrd: expected 5 columns for two types, got %d"
                         % mix.shape[1])
    with open(path("mix.xrd")) as handle:
        for line in handle:
            if line.startswith("# 2theta"):
                columns = line.split()
                break
    if columns != ["#", "2theta[deg]", "I", "I(Si-Si)", "I(Si-O)", "I(O-O)"]:
        raise SystemExit("FAIL xrd: unexpected columns %s" % " ".join(columns))
    total = mix[:, 1]
    pair_sum = mix[:, 2] + 2.0*mix[:, 3] + mix[:, 4]
    worst = float(np.max(np.abs(total - pair_sum)/np.maximum(np.abs(total), 1.0)))
    if worst > 1.0e-9:
        raise SystemExit("FAIL xrd: pair columns are off the total by %.3e" % worst)
    # --no-partials drops them, like it does in the S(q) table
    run([exe, "static", "-i", dump, "-m", "1:Si,2:O", "-w", "xray", "--no-partials", "--xrd",
         path("mix_nopart.xrd"), "--xrd-lambda", "1.5418", "--xrd-range", "5", "60",
         "--xrd-step", "0.5"])
    if read_matrix(path("mix_nopart.xrd")).shape[1] != 2:
        raise SystemExit("FAIL xrd: --no-partials left pair columns in the table")
    print("  ok   %-42s off the total by %.1e" % ("Si-O pair columns and sum rule", worst))

    # Unit weights remove the form factor, so the bins are multiplicity times
    # LP alone - a second, independent check of the binning, and of the bin
    # width the box derives when --xrd-step is left out.
    sc = generate("xrd_sc.dump", natoms=64, length=8.0, frames=2, mode="lattice",
                  fractions=1.0, seed=3)
    run([exe, "static", "-i", sc, "-w", "unit", "--method", "direct", "--xrd", path("sc.xrd"),
         "--xrd-lambda", "1.5418", "--xrd-range", "40", "100"])
    cubic = read_matrix(path("sc.xrd"))
    if cubic.shape != (4, 3):
        raise SystemExit("FAIL xrd (unit weights): %s table, expected 4 bins "
                         "and one pair column"
                         % (cubic.shape,))
    # one species: the pair column carries the total, direct or not
    identity = float(np.max(np.abs(cubic[:, 1] - cubic[:, 2])
                            / np.maximum(np.abs(cubic[:, 1]), 1.0)))
    if identity > 1.0e-12:
        raise SystemExit("FAIL xrd (unit weights): I(1-1) differs from the total "
                         "by %.3e" % identity)
    for index, (m, mult) in enumerate(((1, 6), (2, 12), (3, 8))):
        tth = 2.0*np.degrees(np.arcsin(np.pi*np.sqrt(m)*1.5418/(4.0*np.pi)))
        want = 64.0*mult*lorentz_polarization(tth)
        got = cubic[index, 1]
        if abs(got - want) > 0.01*want:
            raise SystemExit("FAIL xrd (unit weights): bin %d is %.6g, "
                             "multiplicity x LP predicts %.6g" % (m, got, want))
    if cubic[3, 1] > 1.0e-6*cubic[0, 1]:
        raise SystemExit("FAIL xrd (unit weights): the empty bin holds %.3e"
                         % cubic[3, 1])
    print("  ok   %-42s cubic 100/110/111 in 4 box sized bins"
          % "unit weights, direct method")

    # --no-lp drops exactly the Lorentz-polarization factor.
    run([exe, "static", "-i", ni, "-m", "1:Ni", "-w", "xray", "--no-lp", "--xrd", path("ni_nolp.xrd"),
         *ni_opts, "--xrd-step", "0.2"])
    with open(path("ni_nolp.xrd")) as handle:
        handle.readline()
        if "without the Lorentz-polarization factor" not in handle.readline():
            raise SystemExit("FAIL xrd: --no-lp is not reported in the header")
    plain = read_matrix(path("ni_nolp.xrd"))
    for label, tth, _ in expect:
        j = int(np.argmin(np.abs(xrd[:, 0] - tth)))
        ratio = xrd[j, 1]/plain[j, 1]
        want = lorentz_polarization(tth)
        if abs(ratio - want) > 0.01*want:
            raise SystemExit("FAIL xrd: --no-lp ratio of Ni(%s) is %.4f, LP is %.4f"
                             % (label, ratio, want))
    print("  ok   %-42s LP within 1 per cent" % "--no-lp drops the LP factor")

    # The derived bin width must stay usable over the default range, where the
    # box spacing in two-theta diverges towards 180 degrees: without the cap on
    # its reference angle the whole 1-179 degree range collapses to two bins.
    wide = generate("xrd_wide.dump", natoms=4000, length=35.2, frames=1, mode="fcc",
                    fractions=1.0, seed=3)
    run([exe, "static", "-i", wide, "-m", "1:Ni", "-w", "xray", "--xrd", path("wide.xrd"),
         "--xrd-lambda", "1.541838"])
    wide_xrd = read_matrix(path("wide.xrd"))
    if wide_xrd.shape[0] < 20:
        raise SystemExit("FAIL xrd: the default range gave %d bins at %s"
                         % (wide_xrd.shape[0], "the derived step"))
    print("  ok   %-42s %d bins over the default 1-179 deg range"
          % ("default step stays usable", wide_xrd.shape[0]))

    # `direct` is a second implementation of the same sums, so it has to
    # reproduce every column of the pattern - the pair columns included, which
    # no other check pins against an independent computation.
    cross = ["--xrd-lambda", "1.541838", "--xrd-range", "40", "80", "--xrd-step", "0.5"]
    run([exe, "static", "-i", wide, "-m", "1:Ni", "-w", "xray", "--xrd", path("wide_direct.xrd"),
         *cross])
    run([exe, "static", "-i", wide, "-m", "1:Ni", "-w", "xray", "--method", "direct",
         "--xrd", path("wide_nufft.xrd"), *cross])
    xrd_direct = read_matrix(path("wide_direct.xrd"))
    xrd_nufft = read_matrix(path("wide_nufft.xrd"))
    if xrd_direct.shape != xrd_nufft.shape:
        raise SystemExit("FAIL xrd: direct gives %s columns, nufft %s"
                         % (xrd_direct.shape, xrd_nufft.shape))
    worst = float(np.max(np.abs(xrd_direct[:, 1:] - xrd_nufft[:, 1:])
                         / np.maximum(np.abs(xrd_nufft[:, 1:]), 1.0)))
    if worst > 1.0e-8:
        raise SystemExit("FAIL xrd: direct and nufft differ by %.3e" % worst)
    # one species: the pair column is the total, for direct as well
    identity = float(np.max(np.abs(xrd_direct[:, 1] - xrd_direct[:, 2])
                            / np.maximum(np.abs(xrd_direct[:, 1]), 1.0)))
    if identity > 1.0e-12:
        raise SystemExit("FAIL xrd: direct I(Ni-Ni) differs from the total by %.3e"
                         % identity)
    print("  ok   %-42s max deviation %.2e" % ("direct vs nufft, pattern + pairs", worst))

    # --no-partials still drops them
    run([exe, "static", "-i", wide, "-m", "1:Ni", "-w", "xray", "--method", "direct",
         "--no-partials", "--xrd", path("wide_plain.xrd"), *cross])
    if read_matrix(path("wide_plain.xrd")).shape[1] != 2:
        raise SystemExit("FAIL xrd: direct --no-partials left pair columns")
    print("  ok   %-42s direct writes the total only" % "--no-partials")

    if args.h5read:
        run([exe, "static", "-i", ni, "-m", "1:Ni", "-w", "xray", "--xrd", path("ni.h5"),
             *ni_opts, "--xrd-step", "0.2"])
        with open(path("ni_from_h5.txt"), "w") as handle:
            handle.write(run([args.h5read, path("ni.h5"), "xrd"]).stdout)
        h5_xrd = read_matrix(path("ni_from_h5.txt"))
        if h5_xrd.shape != xrd.shape:
            raise SystemExit("FAIL hdf5 xrd shape %s vs text %s"
                             % (h5_xrd.shape, xrd.shape))
        worst = float(np.max(np.abs(h5_xrd - xrd)/np.maximum(np.abs(xrd), 1.0)))
        if worst > 1.0e-12:
            raise SystemExit("FAIL hdf5 vs text xrd: %.3e" % worst)
        print("  ok   %-42s max deviation %.2e" % ("HDF5 vs text XRD table", worst))

    # The combinations that cannot mean anything have to be refused.
    def expect_error(argv, needle):
        result = subprocess.run(argv, capture_output=True, text=True)
        tail = " ".join(argv[1:])
        if result.returncode == 0:
            raise SystemExit("FAIL xrd: %s was accepted" % tail)
        if needle not in result.stdout + result.stderr:
            raise SystemExit("FAIL xrd: %s did not report %r" % (tail, needle))

    base = [exe, "static", "-i", ni, "-m", "1:Ni", "-w", "xray"]
    expect_error(base + ["--xrd", path("bad.xrd"), "--xrd-lambda", "1.5418",
                         "--method", "debye"], "--method debye sums")
    expect_error(base + ["--xrd", path("bad.xrd"), "--xrd-lambda", "1.5418",
                         "--qmin", "2"], "--qmin")
    expect_error(base + ["--xrd", path("bad.xrd")], "--xrd-lambda")
    expect_error(base + ["--xrd", path("bad.xrd"), "--xrd-lambda", "1.5418",
                         "--xrd-range", "0", "80"], "MIN2TH")
    expect_error(base + ["--xrd-step", "0.1", path("bad.xrd")], "belong to --xrd")
    expect_error(base + ["--no-lp", path("bad.dat")], "belong to --xrd")
    # the 32000 atom box needs 2.3e10 phase evaluations, over the direct limit
    expect_error([exe, "static", "-i", ni, "-m", "1:Ni", "-w", "xray", "--method", "direct",
                  "--xrd", path("bad.xrd"), "--xrd-lambda", "1.541838",
                  "--xrd-range", "40", "80"], "phase evaluations")
    print("  ok   %-42s 7 refusals" % "debye, qmin, lambda, range, stray flags")

    # --- 9e. B2 NiAl: the pair columns of a superlattice -------------------
    # B2 (CsCl type) puts Ni on the cell corners and Al in the cell centres.
    # With h + k + l odd the two sublattices scatter in antiphase, so the Ni-Al
    # cross term is negative there and - with equal sublattice populations and
    # unit weights - the total cancels to zero while both diagonal terms stay
    # large.  That is a pair signature no gas or single species test can show:
    # a missing or sign-flipped cross term leaves the total at N/2 instead of 0.
    print("B2 NiAl superlattice")
    b2_side, b2_a = 8, 2.88
    b2_natoms = 2*b2_side**3
    b2 = generate("b2.dump", natoms=b2_natoms, length=b2_side*b2_a, frames=1,
                  mode="b2", seed=3)
    q_100 = 2.0*np.pi/b2_a
    q_110 = q_100*np.sqrt(2.0)

    def b2_shell(name, q, weight):
        # a window of +-5e-5 1/A holds the six (or twelve) lattice vectors of
        # one family and nothing else, so the shell average is the family mean
        run([exe, "static", "-i", b2, "-m", "1:Ni,2:Al", "-w", weight, "--norm", "mean",
             "--qmin", "%.6f" % (q - 5.0e-5), "--qmax", "%.6f" % (q + 5.0e-5),
             "--nq", "1", path(name)])
        return read_matrix(path(name))[0]

    # IT92 rows of the element table, for the analytic contrast below.
    it92 = {"Ni": ([12.8376, 7.2920, 4.4438, 2.3800],
                   [3.8785, 0.2565, 12.1763, 66.3421], 1.0341),
            "Al": ([6.4202, 1.9002, 1.5936, 1.9646],
                   [3.0387, 0.7426, 31.5472, 85.0886], 1.1151)}

    def form_factor(symbol, q):
        a, b, c = it92[symbol]
        return c + sum(ai*np.exp(-bi*(q/(4.0*np.pi))**2) for ai, bi in zip(a, b))

    shell_aa = b2_natoms/4.0
    shells = {}
    for label, q, cross in (("odd", q_100, -shell_aa), ("even", q_110, shell_aa)):
        row = b2_shell("b2_%s.dat" % label, q, "unit")
        shells[label] = row
        if abs(row[2] - shell_aa) > 1.0e-6*shell_aa or abs(row[4] - shell_aa) > 1.0e-6*shell_aa:
            raise SystemExit("FAIL B2 %s: diagonal pairs are %.6g and %.6g, expected %.6g"
                             % (label, row[2], row[4], shell_aa))
        if abs(row[3] - cross) > 1.0e-6*shell_aa:
            raise SystemExit("FAIL B2 %s: cross term %.6g, expected %.6g"
                             % (label, row[3], cross))
        worst = abs(row[1] - (row[2] + 2.0*row[3] + row[4]))/max(abs(row[1]), 1.0)
        if worst > 1.0e-12:
            raise SystemExit("FAIL B2 %s: sum rule off by %.3e" % (label, worst))
    odd_total = shells["odd"][1]
    even_total = shells["even"][1]
    if abs(odd_total) > 1.0e-6*shell_aa:
        raise SystemExit("FAIL B2: the antiphase total is %.6g, expected 0" % odd_total)
    if abs(even_total - b2_natoms) > 1.0e-6*b2_natoms:
        raise SystemExit("FAIL B2: the in-phase total is %.6g, expected N = %d"
                         % (even_total, b2_natoms))
    print("  ok   %-42s total 0 vs N, pairs +-%.0f" % ("antiphase (100) vs (110)", shell_aa))

    # Chemical weighting: the superlattice intensity is the contrast factor of
    # the two form factors, relative to the in-phase line which is 1 (up to the
    # N of a coherent crystal).
    contrast = ((form_factor("Ni", q_100) - form_factor("Al", q_100))
                / (form_factor("Ni", q_100) + form_factor("Al", q_100)))**2
    odd_x = b2_shell("b2_odd_x.dat", q_100, "xray")
    if abs(odd_x[1] - b2_natoms*contrast) > 1.0e-6*b2_natoms*contrast:
        raise SystemExit("FAIL B2: xray superlattice is %.6g, contrast predicts %.6g"
                         % (odd_x[1], b2_natoms*contrast))
    print("  ok   %-42s S(100) = N (f_Ni-f_Al)^2/(f_Ni+f_Al)^2 = %.4f"
          % ("xray contrast", odd_x[1]/b2_natoms))

    # The pattern keeps the same signature: the cross column is negative at the
    # odd lines and positive at the even ones, and direct reproduces it all.
    xrd_b2 = ["--xrd-lambda", "1.541838", "--xrd-range", "25", "80", "--xrd-step", "0.2"]
    run([exe, "static", "-i", b2, "-m", "1:Ni,2:Al", "-w", "xray", "--xrd", path("b2.xrd"), *xrd_b2])
    run([exe, "static", "-i", b2, "-m", "1:Ni,2:Al", "-w", "xray", "--method", "direct",
         "--xrd", path("b2_direct.xrd"), *xrd_b2])
    b2_xrd = read_matrix(path("b2.xrd"))
    b2_dir = read_matrix(path("b2_direct.xrd"))
    if b2_xrd.shape != b2_dir.shape:
        raise SystemExit("FAIL B2 xrd: %s columns vs %s" % (b2_dir.shape, b2_xrd.shape))
    worst = float(np.max(np.abs(b2_dir[:, 1:] - b2_xrd[:, 1:])
                         / np.maximum(np.abs(b2_xrd[:, 1:]), 1.0)))
    if worst > 1.0e-8:
        raise SystemExit("FAIL B2 xrd: direct and nufft differ by %.3e" % worst)
    for h, k, l, sign in ((1, 0, 0, -1), (1, 1, 0, 1), (1, 1, 1, -1), (2, 0, 0, 1)):
        q = q_100*np.sqrt(h*h + k*k + l*l)
        tth = 2.0*np.degrees(np.arcsin(q*1.541838/(4.0*np.pi)))
        j = int(np.argmin(np.abs(b2_xrd[:, 0] - tth)))
        if abs(b2_xrd[j, 0] - tth) > 0.2:
            raise SystemExit("FAIL B2 xrd: no bin near %.2f deg for (%d%d%d)"
                             % (tth, h, k, l))
        if np.sign(b2_xrd[j, 3]) != sign:
            raise SystemExit("FAIL B2 xrd: (%d%d%d) cross column is %.4g, expected sign %d"
                             % (h, k, l, b2_xrd[j, 3], sign))
        worst = abs(b2_xrd[j, 1] - (b2_xrd[j, 2] + 2.0*b2_xrd[j, 3] + b2_xrd[j, 4]))
        if worst > 1.0e-6*abs(b2_xrd[j, 1]):
            raise SystemExit("FAIL B2 xrd: (%d%d%d) sum rule off by %.3e" % (h, k, l, worst))
    print("  ok   %-42s cross term -/+/-, direct matches to %.1e"
          % ("pattern: (100), (110), (111), (200)", worst))

    # --- 9f. q sampling helper --------------------------------------------
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
            run([exe, "static", "-i", source, *options, path("chosen_q.dat")])
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
        run([exe, "static", "-i", dump, *debye_opts, "--nq", "40", path("deb_coarse.dat")])
        run([exe, "static", "-i", dump, *debye_opts, "--nq", "200", path("deb_fine.dat")])
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
    for method, extra in (("nufft", []), ("direct", []),
                          ("debye", ["--rmax", "6", "--dr", "0.01", "--skin", "0"])):
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
        fsum = x*x*fz[:, 2] + 2.0*x*x*fz[:, 3] + x*x*fz[:, 4]
        worst = float(np.max(np.abs(fsum - fz[:, 1])))
        if worst > 1.0e-9:
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

    # --- 11. dynamic structure factor along a q line (sqcalc dyn) ----------
    print("dynamic structure factor along a q line")
    ref_dyn = os.path.join(HERE, "ref_sqw.py")
    dyn_spec = "4,1,4,1,0,0"          # 5 q points, |q| = 1..4 along x
    dyn_q = ["--qpoints", "line:" + dyn_spec, "--dt", "1", "--maxframes", "20",
             "--lag", "1"]
    nq_dyn = 5
    dyn_dump = generate("dyn_ball.dump", natoms=1000, length=20, frames=40, seed=11,
                        mode="ballistic", temperature=1.0, columns="id type xu yu zu")
    sqcalc_dyn(dyn_dump, "-w", "unit", "--norm", "mean", *dyn_q,
               "--sq", path("dyn_sq.dat"), "--sqw", path("dyn_sqw.dat"),
               "--fqt", path("dyn_fqt.dat"))
    run([sys.executable, ref_dyn, "--input", dyn_dump,
         "--qpoints", "line:" + dyn_spec,
         "--maxframes", "20", "--lag", "1", "--dt", "1", "--weight", "unit",
         "--norm", "mean", "--output-fqt", path("dyn_ref_fqt.dat"),
         "--output-sqw", path("dyn_ref_sqw.dat")])
    fqt = read_matrix(path("dyn_fqt.dat"))
    sqw = read_matrix(path("dyn_sqw.dat"))
    compare(fqt, read_matrix(path("dyn_ref_fqt.dat")), "F(q,t) vs numpy reference", rtol=1.0e-9)
    compare(sqw, read_matrix(path("dyn_ref_sqw.dat")), "S(q,w) vs numpy reference", rtol=1.0e-9)

    # The spectrum integrates to F(q,0), which is the static OUTPUT column; the
    # one-sided folding makes this exact (dw = pi/(maxframes*dt)).
    naxis_dyn = 21
    _, s_static = read_table(path("dyn_sq.dat"))
    integral = sqw[:, 4].reshape(nq_dyn, naxis_dyn).sum(axis=1) * (np.pi/20.0)
    compare(s_static, integral, "sum rule int S(q,w) dw = static S(q)", rtol=1.0e-9)
    compare(s_static, fqt[:, 4].reshape(nq_dyn, naxis_dyn)[:, 0], "F(q,0) = static S(q)",
            rtol=1.0e-9)

    # Brownian particles: F(q,t) = exp(-D q^2 t) (the amplitude carries the
    # frozen-in structure factor of the configuration, which decays away over
    # the trajectory, so only the decay rate is a clean observable).
    diff_dump = generate("dyn_diff.dump", natoms=500, length=20, frames=1000, seed=3,
                         mode="diffusive", diffusivity=0.1, columns="id type xu yu zu")
    run([exe, "dyn", "-i", diff_dump, "-m", "1:Si,2:O", "-w", "unit", "--norm", "mean",
         "--qpoints", "line:1,2,3,1,0,0", "--dt", "1", "--maxframes", "8",
         "--lag", "1",
         "--fqt", path("dyn_diff_fqt.dat"), "--sq", path("dyn_diff_sq.dat")])
    diff_f = read_matrix(path("dyn_diff_fqt.dat"))
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
    sqcalc_dyn(dyn_wrapped, "-w", "unit", "--norm", "mean", *dyn_q,
               "--sq", path("dyn_sq_wrap.dat"), "--sqw", path("dyn_sqw_wrap.dat"))
    compare(read_matrix(path("dyn_sqw_wrap.dat")), sqw,
            "wrapped dump (reader unwraps) vs xu", rtol=1.0e-9)

    # --dt is the time step of the trajectory: the same frames written every
    # fifth step with dt/5 have to give the same spectrum.
    dyn_step = generate("dyn_ball_step.dump", natoms=1000, length=20, frames=40, seed=11,
                        mode="ballistic", temperature=1.0, columns="id type xu yu zu",
                        step=5)
    run([exe, "dyn", "-i", dyn_step, "-m", "1:Si,2:O", "-w", "unit", "--norm", "mean",
         "--qpoints", "line:" + dyn_spec, "--dt", "0.2", "--maxframes", "20",
         "--lag", "1",
         "--sqw", path("dyn_sqw_step.dat"), "--sq", path("dyn_sq_step.dat")])
    compare(read_matrix(path("dyn_sqw_step.dat")), sqw,
            "dt = 0.2 over 5 steps equals dt = 1 over 1", rtol=1.0e-9)

    # Partial spectra: with unit weights and --norm n the OVITO sum rule
    # S(q,w) = S_aa + 2 S_ab + S_bb holds sample by sample.
    run([exe, "dyn", "-i", dyn_dump, "-m", "1:Si,2:O", "-w", "unit", "--norm", "n", *dyn_q,
         "--sqw", path("dyn_part.dat"), "--sq", path("dyn_sq_n.dat")])
    part = read_matrix(path("dyn_part.dat"))
    worst = float(np.max(np.abs(part[:, 4] - (part[:, 5] + 2.0*part[:, 6] + part[:, 7]))))
    if worst > 1.0e-9:
        raise SystemExit("FAIL dynamic partial sum rule is off by %.3e" % worst)
    print("  ok   %-42s sum rule %.1e"
          % ("dynamic partials: S = Saa + 2Sab + Sbb", worst))

    # Chemical weighting: the weighted total follows the same post-processing
    # convention as the static partials (neutron b, X-ray f(q)).
    for weight, norm in (("neutron", "mean"), ("xray", "self")):
        sqcalc_dyn(dyn_dump, "-w", weight, "--norm", norm, *dyn_q,
                   "-m", "1:Si,2:O",
                   "--sq", path("dyn_w_%s.dat" % weight),
                   "--sqw", path("dyn_w_%s_sqw.dat" % weight))
        run([sys.executable, ref_dyn, "--input", dyn_dump,
             "--qpoints", "line:" + dyn_spec,
             "--maxframes", "20", "--lag", "1", "--dt", "1", "--weight", weight,
             "--norm", norm, "--mapping", "1:Si,2:O",
             "--output-fqt", path("dyn_w_%s_ref_fqt.dat" % weight),
             "--output-sqw", path("dyn_w_%s_ref_sqw.dat" % weight)])
        compare(read_matrix(path("dyn_w_%s_sqw.dat" % weight)),
                read_matrix(path("dyn_w_%s_ref_sqw.dat" % weight)),
                "%s weighted S(q,w) vs reference" % weight, rtol=1.0e-9)

    # --lag thins the time origins; the reference applies the same rule.
    run([exe, "dyn", "-i", dyn_dump, "-m", "1:Si,2:O", "-w", "unit", "--norm", "mean",
         "--qpoints", "line:" + dyn_spec, "--dt", "1", "--maxframes", "20",
         "--lag", "3",
         "--fqt", path("dyn_lag_fqt.dat"), "--sq", path("dyn_lag_sq.dat")])
    run([sys.executable, ref_dyn, "--input", dyn_dump,
         "--qpoints", "line:" + dyn_spec,
         "--maxframes", "20", "--lag", "3", "--dt", "1", "--weight", "unit",
         "--norm", "mean", "--output-fqt", path("dyn_lag_ref_fqt.dat"),
         "--output-sqw", path("dyn_lag_ref_sqw.dat")])
    compare(read_matrix(path("dyn_lag_fqt.dat")), read_matrix(path("dyn_lag_ref_fqt.dat")),
            "lag stride 3 vs reference", rtol=1.0e-9)

    # --- four point structure factor S4(q,t), Q(t) and chi4(t) ------------
    print("four-point structure factor (--s4/--chi4)")
    ref_four = os.path.join(HERE, "ref_s4.py")
    s4_spec = "1,2,3,1,0,0"            # 2 q points, |q| = 2 and 3 along x
    s4_cutoff = "0.5"
    s4_base = ["--qpoints", "line:" + s4_spec, "--dt", "1", "--maxframes", "8"]
    s4_q = s4_base + ["--lag", "1"]
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", "--norm", "mean", *s4_q,
         "--s4-cutoff", s4_cutoff, "--s4", path("s4.dat"),
         "--chi4", path("chi4.dat"), "--sq", path("s4_sq.dat")])
    run([sys.executable, ref_four, "--input", diff_dump, "--qpoints", "line:" + s4_spec,
         "--maxframes", "8", "--lag", "1", "--dt", "1", "--cutoff", s4_cutoff,
         "--output-s4", path("s4_ref.dat"), "--output-chi4", path("chi4_ref.dat")])
    s4 = read_matrix(path("s4.dat"))
    chi4 = read_matrix(path("chi4.dat"))
    compare(s4, read_matrix(path("s4_ref.dat")), "S4(q,t) vs numpy reference", rtol=1.0e-9)
    compare(chi4, read_matrix(path("chi4_ref.dat")), "Q(t)/chi4(t) vs numpy reference",
            rtol=1.0e-9)
    if abs(chi4[0, 1] - 1.0) > 1.0e-12 or abs(chi4[0, 2]) > 1.0e-12:
        raise SystemExit("FAIL Q(0)=1 and chi4(0)=0")
    print("  ok   %-42s Q(0)=%.6f chi4(0)=%.1e"
          % ("Q(0) and chi4(0)", chi4[0, 1], chi4[0, 2]))

    # The overlap scan compacts one block of atoms per thread and joins the
    # segments in thread order, so an S4 run cannot depend on the thread
    # count.  The two tables have to be identical, not merely close.
    for threads in ("1", "2"):
        run([exe, "dyn", "-i", diff_dump, "-w", "unit", *s4_base, "--lag", "1",
             "--s4-cutoff", s4_cutoff, "--s4", path("s4_scan_t%s.dat" % threads),
             "--sq", path("s4_scan_t%s_sq.dat" % threads)],
            env={"OMP_NUM_THREADS": threads})
    compare(read_matrix(path("s4_scan_t1.dat")), read_matrix(path("s4_scan_t2.dat")),
            "S4 is independent of the thread count", rtol=0.0)

    # The q = 0 row of the S4 table is the scalar chi4(t).
    s4_q0 = ["--qpoints", "line:2,0,2,1,0,0", "--dt", "1", "--maxframes", "8",
             "--lag", "1"]
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", *s4_q0, "--s4-cutoff", s4_cutoff,
         "--s4", path("s4_q0.dat"), "--sq", path("s4_q0_sq.dat")])
    compare(read_matrix(path("s4_q0.dat"))[:9, 4].reshape(-1, 1),
            chi4[:, 2].reshape(-1, 1), "S4(q=0,t) = chi4(t)", rtol=1.0e-12)

    # --partials must not add S4 columns.
    run([exe, "dyn", "-i", diff_dump, "-m", "1:Si,2:O", "--partials", *s4_q,
         "--s4-cutoff", s4_cutoff, "--s4", path("s4_part.dat"), "--sq", path("s4_part_sq.dat")])
    if read_matrix(path("s4_part.dat")).shape != s4.shape:
        raise SystemExit("FAIL --partials changed the S4 table layout")
    print("  ok   %-42s no partial columns" % "--partials is ignored by S4")

    # --chi4 alone must reproduce the chi4 of the combined run.
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", *s4_q, "--s4-cutoff", s4_cutoff,
         "--chi4", path("chi4_only.dat"), "--sq", path("chi4_only_sq.dat")])
    compare(read_matrix(path("chi4_only.dat")), chi4, "chi4 without --s4", rtol=1.0e-9)

    # --format overrides the suffix for the dynamic outputs.
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", *s4_q, "--s4-cutoff", s4_cutoff,
         "--format", "text", "--s4", path("s4_forced.h5"),
         "--chi4", path("chi4_forced.h5"), "--sq", path("s4_forced_sq.dat")])
    compare(read_matrix(path("s4_forced.h5")), s4, "S4 --format text override", rtol=1.0e-9)
    compare(read_matrix(path("chi4_forced.h5")), chi4, "chi4 --format text override",
            rtol=1.0e-9)

    # A run with only --s4/--chi4 must not need the coherent ring buffers;
    # adding --sqw has to leave S4, chi4 and the static S(q) table unchanged.
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", *s4_q, "--s4-cutoff", s4_cutoff,
         "--s4", path("s4_coherent.dat"), "--chi4", path("chi4_coherent.dat"),
         "--sqw", path("s4_coherent_sqw.dat"), "--sq", path("s4_coherent_sq.dat")])
    compare(read_matrix(path("s4_coherent.dat")), s4, "S4-only == S4 with --sqw", rtol=1.0e-9)
    compare(read_matrix(path("chi4_coherent.dat")), chi4, "chi4-only == chi4 with --sqw",
            rtol=1.0e-9)
    compare(read_table(path("s4_sq.dat"))[1].reshape(-1, 1),
            read_table(path("s4_coherent_sq.dat"))[1].reshape(-1, 1),
            "S4-only static S(q) == full", rtol=1.0e-9)

    # --- reciprocal-lattice grid sampling (--qpoints grid) ------------------
    print("reciprocal-lattice grid sampling (--qpoints grid)")
    q1 = 2.0*math.pi/20.0              # smallest reciprocal vector of the box
    grid_base = ["--qpoints", "grid:%.6f" % (2.0*q1), "--dt", "1",
                 "--maxframes", "8", "--lag", "1"]
    grid_run = run([exe, "dyn", "-i", diff_dump, "-w", "unit", *grid_base, "--s4-cutoff", s4_cutoff,
                    "--s4", path("s4_grid.dat"), "--chi4", path("chi4_grid.dat"), "--sq", 
                    path("s4_grid_sq.dat")])
    grid = read_matrix(path("s4_grid.dat"))
    grid_chi4 = read_matrix(path("chi4_grid.dat"))
    # By default the rows are the |q| shells: the shells |n| = 1, sqrt(2),
    # sqrt(3), 2 plus the Gamma point.
    shells = np.unique(np.round(grid[:, :3], 6), axis=0)
    if shells.shape[0] != 5:
        raise SystemExit("FAIL grid: expected 5 shell rows, found %d" % shells.shape[0])
    if np.any(np.abs(shells[:, 0]) > 1.0e-9) or np.any(np.abs(shells[:, 1]) > 1.0e-9):
        raise SystemExit("FAIL grid: the shell rows are not keyed by |q|")
    # The per-mode rows are still available on request.
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", *grid_base, "--keep-modes",
         "--s4-cutoff", s4_cutoff, "--s4", path("s4_grid_modes.dat"), "--sq", 
         path("s4_grid_modes_sq.dat")])
    mode_table = read_matrix(path("s4_grid_modes.dat"))
    modes = np.unique(np.round(mode_table[:, :3], 6), axis=0)
    # The shells hold 6+12+8+6 = 32 vectors; +- keeps 16 and Gamma adds one.
    if modes.shape[0] != 17:
        raise SystemExit("FAIL grid: expected 17 modes, found %d" % modes.shape[0])
    # The run has to report the isotropy probe next to the mode counts.
    if "isotropy" not in grid_run.stderr:
        raise SystemExit("FAIL grid: the isotropy probe is missing from the summary")

    def rows_at(table, qvec):
        sel = table[np.all(np.abs(table[:, :3] - np.array(qvec)) < 1.0e-6, axis=1)]
        return sel[np.argsort(sel[:, 3])]

    # Gamma is the q -> 0 limit, i.e. the scalar chi4(t).
    compare(rows_at(grid, [0.0, 0.0, 0.0])[:, 4].reshape(-1, 1),
            grid_chi4[:, 2].reshape(-1, 1), "grid: S4(q=0,t) = chi4(t)", rtol=1.0e-12)

    # A lattice mode has to agree with a q line through the same vector.
    run([exe, "dyn", "-i", diff_dump, "-w", "unit",
         *dyn_line("1,0,%.12f,1,0,0" % q1), "--dt", "1", "--maxframes", "8",
         "--lag", "1", "--s4-cutoff", s4_cutoff, "--s4", path("s4_grid_line.dat"), "--sq", 
         path("s4_grid_line_sq.dat")])
    compare(rows_at(mode_table, [q1, 0.0, 0.0])[:, 4].reshape(-1, 1),
            rows_at(read_matrix(path("s4_grid_line.dat")), [q1, 0.0, 0.0])[:, 4].reshape(-1, 1),
            "grid: a lattice mode matches the q line", rtol=1.0e-12)

    # The Gamma row is the q -> 0 limit in either flavour.
    compare(rows_at(mode_table, [0.0, 0.0, 0.0])[:, 4].reshape(-1, 1),
            grid_chi4[:, 2].reshape(-1, 1), "grid: per-mode Gamma row = chi4(t)",
            rtol=1.0e-12)

    # The mode budget drops whole shells, keeps the two ends of the q range and
    # is reported in the table header.
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", *grid_base, "--modes", "12",
         "--s4-cutoff", s4_cutoff, "--s4", path("s4_grid_thin.dat"), "--sq", 
         path("s4_grid_thin_sq.dat")])
    thin_modes = np.unique(np.round(read_matrix(path("s4_grid_thin.dat"))[:, :3], 6), axis=0)
    if not 1 < thin_modes.shape[0] <= 12:
        raise SystemExit("FAIL grid: --modes 12 gave %d modes" % thin_modes.shape[0])
    if not any("thinned" in line for line in open(path("s4_grid_thin.dat"))):
        raise SystemExit("FAIL grid: the thinning is not reported in the header")
    print("  ok   %-42s %d modes, thinned to %d"
          % ("grid modes and the mode budget", modes.shape[0], thin_modes.shape[0]))

    # The coherent and self amplitudes use the same separable evaluation, so
    # compare them against a q line, which still uses the direct loop.  One
    # vector per axis count: at (1,0,0) only the first factor matters, at
    # (1,1,0) the first two do.
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", *grid_base, "--keep-modes",
         "--fqt", path("fq_grid.dat"), "--fqt-self", path("fs_grid.dat"), "--sq", 
         path("fq_grid_sq.dat")])
    for direction, qvec in (("1,0,0", [q1, 0.0, 0.0]), ("1,1,0", [q1, q1, 0.0])):
        run([exe, "dyn", "-i", diff_dump, "-w", "unit",
             *dyn_line("1,0,%.12f,%s" % (math.sqrt(sum(v*v for v in qvec)), direction)),
             "--dt", "1", "--maxframes", "8", "--lag", "1", "--fqt", path("fq_line.dat"),
             "--fqt-self", path("fs_line.dat"), "--sq", path("fq_line_sq.dat")])
        for table, label in (("fq", "F(q,t)"), ("fs", "F_s(q,t)")):
            # The two evaluations differ only in the order the atoms are
            # summed, which the 13 significant digits of the table cannot
            # resolve; a wrong phase or a wrong atom set would be orders of
            # magnitude larger than this.
            compare(rows_at(read_matrix(path("%s_grid.dat" % table)), qvec)[:, 4].reshape(-1, 1),
                    rows_at(read_matrix(path("%s_line.dat" % table)), qvec)[:, 4].reshape(-1, 1),
                    "grid: %s matches the line at (%s)" % (label, direction), rtol=1.0e-11)
    print("  ok   %-42s F and F_s, one and two axes" % "grid: separable rho and F_s")

    # A budget that only the orbit route can meet must be met, not reported as
    # impossible: three modes are the Gamma point plus two shell
    # representatives, while one whole shell already costs four.
    tight = run([exe, "dyn", "-i", diff_dump, "-w", "unit", *grid_base, "--modes", "3",
                 "--s4-cutoff", s4_cutoff, "--s4", path("s4_grid_tight.dat"), "--sq", 
                 path("s4_grid_tight_sq.dat")])
    tight_rows = np.unique(np.round(read_matrix(path("s4_grid_tight.dat"))[:, :3], 6), axis=0)
    if tight_rows.shape[0] != 3:
        raise SystemExit("FAIL grid: --modes 3 kept %d rows" % tight_rows.shape[0])
    if "could not be met" in tight.stderr:
        raise SystemExit("FAIL grid: --modes 3 was reported as impossible")
    print("  ok   %-42s %d rows" % ("grid: the budget ladder reaches the orbit route",
                                    tight_rows.shape[0]))

    # A shell with two orbits: |n|^2 = 9 holds six (3,0,0) vectors and twenty
    # four (2,2,1) ones, so a shell average that ignores the multiplicities is
    # wrong.  Under --thin orbits one representative per orbit survives and
    # the collapsed row must be their 6:24 weighted mean, which is an algebraic
    # identity and therefore comparable to machine precision.
    def rows_at_q(table, q):
        sel = table[np.abs(np.linalg.norm(table[:, :3], axis=1) - q) < 1.0e-4]
        return sel[np.argsort(sel[:, 3])]

    def cubic_orbit_size(miller):
        return len({tuple(np.array(perm)*np.array(sign))
                    for perm in permutations(miller)
                    for sign in product([1, -1], repeat=3)})

    orbit_run = ["--qpoints", "grid:%.6f" % (3.1*q1), "--dt", "1",
                 "--maxframes", "8", "--lag", "1", "--thin", "orbits",
                 "--modes", "10"]
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", *orbit_run, "--keep-modes",
         "--s4-cutoff", s4_cutoff, "--s4", path("s4_orbit_modes.dat"), "--sq", 
         path("s4_orbit_modes_sq.dat")])
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", *orbit_run, "--s4-cutoff", s4_cutoff,
         "--s4", path("s4_orbit_avg.dat"), "--sq", path("s4_orbit_avg_sq.dat")])
    shell = rows_at_q(read_matrix(path("s4_orbit_modes.dat")), 3.0*q1)
    vectors = np.unique(np.round(shell[:, :3], 6), axis=0)
    if vectors.shape[0] != 2:
        raise SystemExit("FAIL grid: expected two orbit representatives at |n| = 3, got %d"
                         % vectors.shape[0])
    curves = []
    for vector in vectors:
        rows = shell[np.all(np.abs(shell[:, :3] - vector) < 1.0e-6, axis=1)]
        curves.append(rows[np.argsort(rows[:, 3])][:, 4])
    curves = np.array(curves)
    weights = np.array([cubic_orbit_size(np.round(vector/q1).astype(int)) for vector in vectors],
                       dtype=float)
    if sorted(weights) != [6.0, 24.0]:
        raise SystemExit("FAIL grid: unexpected orbit multiplicities %s" % weights)
    weighted = (weights[:, None]*curves).sum(axis=0)/weights.sum()
    collapsed = rows_at_q(read_matrix(path("s4_orbit_avg.dat")), 3.0*q1)[:, 4]
    # Same reason as above: the identity is exact in arithmetic, so what is
    # left is the print resolution of the tables.
    compare(collapsed.reshape(-1, 1), weighted.reshape(-1, 1),
            "grid: the shell row is the multiplicity weighted mean", rtol=1.0e-11)
    if np.allclose(collapsed, curves.mean(axis=0), rtol=1.0e-9):
        raise SystemExit("FAIL grid: the weighted and plain means coincide, test is vacuous")
    print("  ok   %-42s weights %s" % ("grid: two-orbit shell weighting",
                                       weights.astype(int)))

    # The partial columns have to add up to the total and to match the static
    # bin of the same shell.  The grid row divides the weighted partial
    # numerator by the number of lattice vectors the shell averages; dividing by
    # the number of modes instead leaves it larger by the mean mode weight (two
    # after the +- reduction), which is what this checks.
    run([exe, "dyn", "-i", diff_dump, "-m", "1:Si,2:O", "-w", "unit", "--qpoints",
         "grid:%.6f" % (2.0*q1), "--dt", "1", "--maxframes", "8", "--lag", "1",
         "--sq", path("grid_partials.dat")])
    grid_part = read_matrix(path("grid_partials.dat"))
    used = np.abs(grid_part[:, 1]) > 1.0e-12
    ratio = ((grid_part[used, 2] + 2.0*grid_part[used, 3] + grid_part[used, 4])
             /grid_part[used, 1])
    worst = float(np.max(np.abs(ratio - 1.0)))
    if worst > 1.0e-12:
        raise SystemExit("FAIL grid: the partial columns are off by a factor of "
                         "up to %.6f" % (1.0 + worst))
    run([exe, "static", "-i", diff_dump, "-m", "1:Si,2:O", "-w", "unit",
         "--qmin", "%.6f" % (q1 - 1.0e-3), "--qmax", "%.6f" % (q1 + 1.0e-3),
         "--nq", "1", path("grid_partials_bin.dat")])
    shell_row = grid_part[np.abs(grid_part[:, 0] - q1) < 1.0e-6]
    if shell_row.shape[0] != 1:
        raise SystemExit("FAIL grid: the |n| = 1 shell is missing from the table")
    compare(shell_row[0:1, 1:5], read_matrix(path("grid_partials_bin.dat"))[0:1, 1:5],
            "grid: partials match the static bin", rtol=1.0e-9)
    print("  ok   %-42s %d rows, max |ratio-1| = %.1e"
          % ("grid partial sum rule", int(np.sum(used)), worst))

    # Under orbit thinning the weights are no longer uniform (one representative
    # per orbit carries the full multiplicity), so the divisor has to be the
    # weight sum there too; counting modes instead left a factor of the mean
    # orbit multiplicity rather than the factor of two seen above.
    run([exe, "dyn", "-i", diff_dump, "-m", "1:Si,2:O", "-w", "unit", "--qpoints",
         "grid:%.6f" % (2.0*q1), "--modes", "10", "--thin", "orbits", "--dt", "1",
         "--maxframes", "8", "--lag", "1", "--sq", path("grid_partials_thin.dat")])
    thin = read_matrix(path("grid_partials_thin.dat"))
    used_thin = np.abs(thin[:, 1]) > 1.0e-12
    worst_thin = float(np.max(np.abs((thin[used_thin, 2] + 2.0*thin[used_thin, 3]
                                      + thin[used_thin, 4])/thin[used_thin, 1] - 1.0)))
    if worst_thin > 1.0e-12:
        raise SystemExit("FAIL grid: the thinned partial columns are off by a factor of "
                         "up to %.6f" % (1.0 + worst_thin))
    print("  ok   %-42s %d rows, max |ratio-1| = %.1e"
          % ("grid partial sum rule, orbit thinning", int(np.sum(used_thin)), worst_thin))

    # --- skewed triclinic cells (point group search) ------------------------
    print("skewed triclinic cells")
    # The point-group search is a bootstrap over Miller indices, and a strongly
    # sheared cell used to abort the run.  Two cases: one whose operations fit
    # inside the raised bound (group found, no note) and one that does not
    # (trivial group, note reported, sum rules unchanged).
    skew_q1 = 2.0*math.pi/18.0
    skewed = generate("skew.dump", natoms=400, length=18, frames=12, seed=17,
                      mode="diffusive", tilt="3 4 2", columns="id type xu yu zu")
    skew_run = run([exe, "dyn", "-i", skewed, "-w", "unit", "--qpoints",
                    "grid:%.6f" % (2.0*skew_q1), "--dt", "1", "--maxframes", "8",
                    "--lag", "1", "--sq", path("skew_sq.dat"),
                    "--fqt", path("skew_fqt.dat"), "--fqt-self", path("skew_fs.dat")])
    if "point group order 2" not in skew_run.stderr:
        raise SystemExit("FAIL skew: the point group of the sheared cell is not reported")
    if "orbit thinning is off" in skew_run.stderr:
        raise SystemExit("FAIL skew: fallback reported for a cell inside the bound")
    skew_q = read_matrix(path("skew_sq.dat"))
    skew_f = read_matrix(path("skew_fqt.dat"))
    skew_s = read_matrix(path("skew_fs.dat"))
    compare(skew_f[skew_f[:, 3] == 0.0, 4].reshape(-1, 1), skew_q[:, 1].reshape(-1, 1),
            "skew cell: F(q,0) = S(q)", rtol=0.0)
    compare(skew_s[skew_s[:, 3] == 0.0, 4].reshape(-1, 1),
            np.ones((skew_q.shape[0], 1)), "skew cell: F_s(q,0) = 1", rtol=0.0)

    # Past the bound the build continues with the trivial group.  This cell has
    # a 2 A edge, which needs Miller indices beyond the search box.
    extreme = path("skew_extreme.dump")
    with open(extreme, "w") as handle:
        for frame in range(3):
            handle.write("ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n4\n" % frame)
            handle.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n"
                         "0.0 9.0 3.0\n0.0 17.0 4.0\n0.0 18.0 2.0\n")
            handle.write("ITEM: ATOMS id type xu yu zu\n")
            for index, kind, x, y, z in ((1, 1, 0.5, 0.5, 0.5), (2, 1, 1.2, 3.0, 4.0),
                                         (3, 2, 1.5, 7.0, 9.0), (4, 2, 0.8, 11.0, 14.0)):
                handle.write("%d %d %.6f %.6f %.6f\n" % (index, kind, x, y, z))
    extreme_run = run([exe, "dyn", "-i", extreme, "-w", "unit", "--qpoints", "grid:1.0",
                       "--dt", "1", "--maxframes", "2", "--lag", "1",
                       "--sq", path("skew_ext_sq.dat"), "--fqt", path("skew_ext_fqt.dat"),
                       "--fqt-self", path("skew_ext_fs.dat")])
    if "orbit thinning is off" not in extreme_run.stderr:
        raise SystemExit("FAIL skew: the point group fallback is not reported")
    ext_q = read_matrix(path("skew_ext_sq.dat"))
    ext_f = read_matrix(path("skew_ext_fqt.dat"))
    ext_s = read_matrix(path("skew_ext_fs.dat"))
    # Brute force: ten +- reduced vectors below |q| = 1, one per shell, plus the
    # Gamma row the grid table carries (S4 at q -> 0 is chi4).
    if ext_q.shape[0] != 11:
        raise SystemExit("FAIL skew: expected 10 shells plus Gamma, found %d"
                         % ext_q.shape[0])
    compare(ext_f[ext_f[:, 3] == 0.0, 4].reshape(-1, 1), ext_q[:, 1].reshape(-1, 1),
            "skew cell: fallback keeps F(q,0) = S(q)", rtol=0.0)
    compare(ext_s[ext_s[:, 3] == 0.0, 4].reshape(-1, 1),
            np.ones((ext_q.shape[0], 1)), "skew cell: fallback keeps F_s(q,0) = 1",
            rtol=0.0)
    print("  ok   %-42s %s" % ("point group fallback",
                               "10 shells, trivial group"))

    # --- lattice-shell window (--qpoints powder) ----------------------------
    print("lattice-shell window (--qpoints powder)")

    def powder_counts(text):
        """Lattice vectors and shells of a powder summary, for the checks."""
        head, _, tail = text.replace("\n", " ").partition(" lattice vectors in ")
        return int(head.split()[-1]), int(tail.split()[0])

    # The window average uses the same estimator as the static table's bins, so
    # the two agree in value while their row labels differ by design: the
    # static table writes the bin centre, the powder row writes <|q|>.
    powder_dq = 0.125
    powder_run = run([exe, "dyn", "-i", diff_dump, "-w", "unit", "--qpoints",
                      "powder:%.6f,dq=%.6f" % (2.5, powder_dq), "--dt", "1",
                      "--maxframes", "8", "--lag", "1",
                      "--sq", path("powder_sq.dat"), "--fqt", path("powder_fqt.dat"),
                      "--fqt-self", path("powder_fs.dat")])
    if "window held no shell of its own" in powder_run.stderr:
        raise SystemExit("FAIL powder: unexpected nearest-shell fallback")
    powder_q = read_matrix(path("powder_sq.dat"))
    powder_f = read_matrix(path("powder_fqt.dat"))
    powder_s = read_matrix(path("powder_fs.dat"))
    if powder_q.shape[0] != 1:
        raise SystemExit("FAIL powder: expected one row, found %d" % powder_q.shape[0])
    if not 2.5 - powder_dq - 1.0e-9 <= powder_q[0, 0] <= 2.5 + powder_dq + 1.0e-9:
        raise SystemExit("FAIL powder: the label %.6f is outside the window" % powder_q[0, 0])
    if abs(powder_q[0, 0] - 2.5) < 1.0e-6:
        raise SystemExit("FAIL powder: the row is labelled with the requested |q|")
    compare(powder_f[powder_f[:, 3] == 0.0, 4].reshape(-1, 1),
            powder_q[:, 1].reshape(-1, 1), "powder: F(q,0) = S(q)", rtol=0.0)
    compare(powder_s[powder_s[:, 3] == 0.0, 4].reshape(-1, 1), np.ones((1, 1)),
            "powder: F_s(q,0) = 1", rtol=0.0)
    run([exe, "static", "-i", diff_dump, "-w", "unit",
         "--qmin", "%.6f" % (2.5 - powder_dq), "--qmax", "%.6f" % (2.5 + powder_dq),
         "--nq", "1", path("powder_bin.dat")])
    compare(powder_q[:, 1].reshape(-1, 1),
            read_matrix(path("powder_bin.dat"))[:, 1].reshape(-1, 1),
            "powder: equals the static bin of the same window", rtol=1.0e-9)
    with open(path("powder_fqt.dat")) as handle:
        header = "".join(line for line in handle if line.startswith("#"))
    for needle in ("powder shell", "requested |q|", "mean |q|", "dq"):
        if needle not in header:
            raise SystemExit("FAIL powder: %r missing from the table header" % needle)
    print("  ok   %-42s %s" % ("powder row, header and static bin",
                               "%.6f, %d vectors in %d shells"
                               % (powder_q[0, 0], *powder_counts(powder_run.stderr))))

    # M is the number of lattice vectors the window should hold.  The width it
    # implies is capped at Q/20, so targets past the cap give exactly the window
    # of an explicit dq at the cap, while the default (M = 50) stays narrower.
    default_run = run([exe, "dyn", "-i", diff_dump, "-w", "unit", "--qpoints", "powder:2.5",
                       "--dt", "1", "--maxframes", "8", "--lag", "1",
                       "--sq", path("powder_def.dat")])
    wide = run([exe, "dyn", "-i", diff_dump, "-w", "unit", "--qpoints", "powder:2.5,500",
                "--dt", "1", "--maxframes", "8", "--lag", "1", "--sq", path("powder_m.dat")])
    for count, name in ((2000, "powder_cap_a.dat"), (20000, "powder_cap_b.dat")):
        run([exe, "dyn", "-i", diff_dump, "-w", "unit", "--qpoints",
             "powder:2.5,%d" % count, "--dt", "1", "--maxframes", "8", "--lag", "1",
             "--sq", path(name)])
    compare(read_matrix(path("powder_cap_a.dat")), read_matrix(path("powder_cap_b.dat")),
            "powder: the window width is capped at Q/20", rtol=0.0)
    compare(read_matrix(path("powder_cap_a.dat")), powder_q,
            "powder: a capped M equals dq at the cap", rtol=0.0)
    n_big, s_big = powder_counts(wide.stderr)
    n_def, s_def = powder_counts(default_run.stderr)
    if n_big <= n_def or s_big < s_def:
        raise SystemExit("FAIL powder: M = 500 did not widen the window (%d <= %d)"
                         % (n_big, n_def))
    if abs(read_matrix(path("powder_def.dat"))[0, 1]
           - read_matrix(path("powder_m.dat"))[0, 1]) < 1.0e-9:
        raise SystemExit("FAIL powder: the default and M = 500 windows agree, test is vacuous")
    if abs(read_matrix(path("powder_m.dat"))[0, 1]
           - read_matrix(path("powder_cap_a.dat"))[0, 1]) < 1.0e-9:
        raise SystemExit("FAIL powder: the capped and uncapped windows agree, test is vacuous")
    print("  ok   %-42s %d vs %d lattice vectors"
          % ("powder: M sets the window", n_big, n_def))

    # A window that falls between two shells of a coarse lattice is widened to
    # the nearest one, and the summary says so.
    nearest_run = run([exe, "dyn", "-i", extreme, "-w", "unit", "--qpoints",
                       "powder:0.5465", "--dt", "1", "--maxframes", "2", "--lag", "1",
                       "--sq", path("powder_near.dat"),
                       "--fqt-self", path("powder_near_fs.dat")])
    if "window held no shell of its own" not in nearest_run.stderr:
        raise SystemExit("FAIL powder: the nearest-shell fallback is not reported")
    near_q = read_matrix(path("powder_near.dat"))
    if near_q.shape[0] != 1:
        raise SystemExit("FAIL powder: the fallback wrote %d rows" % near_q.shape[0])
    # The two bracketing shells are 0.516701 and 0.576164, so the nearest one
    # is the upper of the pair.
    if abs(near_q[0, 0] - 0.576164) > 1.0e-5:
        raise SystemExit("FAIL powder: the fallback label %.6f is not the nearest shell"
                         % near_q[0, 0])
    compare(read_matrix(path("powder_near_fs.dat"))[0:1, 4].reshape(-1, 1), np.ones((1, 1)),
            "powder: fallback keeps F_s(q,0) = 1", rtol=0.0)
    print("  ok   %-42s %s" % ("powder: nearest-shell fallback",
                               "|q| = %.6f" % near_q[0, 0]))

    # A box whose shells are coarser than the window, with the nearest shell
    # *above* the request: finding it means looking past Q + dq, which is what
    # the enumeration of the first attempt stops at.
    coarse = path("powder_coarse.dump")
    with open(coarse, "w") as handle:
        for frame in range(3):
            handle.write("ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n4\n" % frame)
            handle.write("ITEM: BOX BOUNDS pp pp pp\n0.0 4.0\n0.0 4.0\n0.0 4.0\n")
            handle.write("ITEM: ATOMS id type xu yu zu\n")
            for index, kind, x, y, z in ((1, 1, 0.5, 0.5, 0.5), (2, 1, 1.5, 1.5, 1.5),
                                         (3, 2, 2.5, 2.5, 2.5), (4, 2, 3.0, 1.0, 2.0)):
                handle.write("%d %d %.6f %.6f %.6f\n" % (index, kind, x, y, z))
    coarse_run = run([exe, "dyn", "-i", coarse, "-w", "unit", "--qpoints", "powder:2.0",
                      "--dt", "1", "--maxframes", "2", "--lag", "1",
                      "--sq", path("powder_coarse_sq.dat")])
    if "window held no shell of its own" not in coarse_run.stderr:
        raise SystemExit("FAIL powder: the coarse-box fallback is not reported")
    coarse_q = read_matrix(path("powder_coarse_sq.dat"))
    # The box has |n|^2 = 1 at 1.5708 and |n|^2 = 2 at 2.2214, so the nearest
    # shell to 2.0 is the one above it.
    if abs(coarse_q[0, 0] - 2.221441) > 1.0e-5:
        raise SystemExit("FAIL powder: the fallback took |q| = %.6f instead of the nearest "
                         "shell at 2.221441" % coarse_q[0, 0])
    run([exe, "dyn", "-i", coarse, "-w", "unit", "--qpoints", "grid:3.0", "--dt", "1",
         "--maxframes", "2", "--lag", "1", "--sq", path("powder_coarse_grid.dat")])
    coarse_grid = read_matrix(path("powder_coarse_grid.dat"))
    row = coarse_grid[np.abs(coarse_grid[:, 0] - 2.221441) < 1.0e-5]
    if row.shape[0] != 1:
        raise SystemExit("FAIL powder: the coarse grid table has no shell at 2.221441")
    compare(coarse_q[0:1, 1].reshape(-1, 1), row[0:1, 1].reshape(-1, 1),
            "powder: the fallback averages the nearest shell", rtol=1.0e-9)
    print("  ok   %-42s %s" % ("powder: nearest shell above the window",
                               "|q| = %.6f" % coarse_q[0, 0]))

    # The window is not derivable from the q dataset, so HDF5 records it.
    if args.h5read:
        run([exe, "dyn", "-i", diff_dump, "-w", "unit", "--qpoints", "powder:2.5,dq=0.05",
             "--dt", "1", "--maxframes", "8", "--lag", "1", "--sqw", path("powder_sqw.h5")])
        for name in ("q_requested", "dq", "n_lattice_vectors"):
            value = run([args.h5read, path("powder_sqw.h5"), "attr", name]).stdout
            if not value.strip():
                raise SystemExit("FAIL powder: HDF5 attribute %s is missing" % name)
        sampling = run([args.h5read, path("powder_sqw.h5"), "attr", "sampling"]).stdout
        if "powder" not in sampling:
            raise SystemExit("FAIL powder: HDF5 does not record the sampling")
        print("  ok   %-42s %s" % ("powder: HDF5 window metadata",
                                   "q_requested, dq, n_lattice_vectors"))

    # --- one reciprocal-lattice vector (--qpoints single) --------------------
    print("single reciprocal-lattice vector (--qpoints single)")
    single_base = ["--qpoints", "single:1,0,0", "--dt", "1", "--maxframes", "8",
                   "--lag", "1"]
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", *single_base, "--s4-cutoff", s4_cutoff,
         "--s4", path("s4_single.dat"), "--chi4", path("chi4_single.dat"), "--sq", 
         path("s4_single_sq.dat")])
    single_table = read_matrix(path("s4_single.dat"))
    # The row is keyed by |q|, the indices and the vector go into the header.
    if np.any(np.abs(single_table[:, :2]) > 1.0e-9) or \
            np.max(np.abs(single_table[:, 2] - q1)) > 1.0e-6:
        raise SystemExit("FAIL single: the rows are not labelled by |q|")
    single_header = "".join(line for line in open(path("s4_single.dat")) if line.startswith("#"))
    if "n = (1 0 0" not in single_header or "q = (" not in single_header:
        raise SystemExit("FAIL single: the header does not carry n and the q vector")
    # It has to agree, bit for bit, with both existing routes to that vector.
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", *dyn_line("1,0,%.12f,1,0,0" % q1),
         "--dt", "1", "--maxframes", "8", "--lag", "1", "--s4-cutoff", s4_cutoff,
         "--s4", path("s4_single_line.dat"), "--sq", path("s4_single_line_sq.dat")])
    compare(rows_at_q(single_table, q1)[:, 4].reshape(-1, 1),
            rows_at(read_matrix(path("s4_single_line.dat")), [q1, 0.0, 0.0])[:, 4].reshape(-1, 1),
            "single: matches the q line at that vector", rtol=1.0e-12)
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", *grid_base, "--keep-modes",
         "--s4-cutoff", s4_cutoff, "--s4", path("s4_single_grid.dat"), "--sq", 
         path("s4_single_grid_sq.dat")])
    grid_row = rows_at(read_matrix(path("s4_single_grid.dat")), [q1, 0.0, 0.0])
    compare(rows_at_q(single_table, q1)[:, 4].reshape(-1, 1),
            grid_row[:, 4].reshape(-1, 1),
            "single: matches the grid at that vector", rtol=1.0e-12)
    if np.max(np.abs(rows_at_q(single_table, q1)[:, 4] - grid_row[:, 4])) == 0.0 and \
            np.allclose(grid_row[:, 4], 0.0):
        raise SystemExit("FAIL single: the comparison is vacuous")
    # A vector with more than one nonzero index: at (1,0,0) the phase is
    # t_1(1) alone, because every t_j(0) is 1, so that comparison cannot see a
    # mix-up in the other two axes.
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", "--qpoints", "single:1,1,0", "--dt", "1",
         "--maxframes", "8", "--lag", "1", "--s4-cutoff", s4_cutoff,
         "--s4", path("s4_single_diag.dat"), "--sq", path("s4_single_diag_sq.dat")])
    compare(rows_at_q(read_matrix(path("s4_single_diag.dat")),
                      math.sqrt(2.0)*q1)[:, 4].reshape(-1, 1),
            rows_at(read_matrix(path("s4_single_grid.dat")), [q1, q1, 0.0])[:, 4].reshape(-1, 1),
            "single: a two-axis vector matches the grid", rtol=1.0e-12)
    # The mirrored indices are the same vector up to a sign, and S4 is even.
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", "--qpoints", "single:-1,0,0",
         "--dt", "1", "--maxframes", "8", "--lag", "1", "--s4-cutoff", s4_cutoff,
         "--s4", path("s4_single_neg.dat"), "--sq", path("s4_single_neg_sq.dat")])
    compare(rows_at_q(read_matrix(path("s4_single_neg.dat")), q1)[:, 4].reshape(-1, 1),
            rows_at_q(single_table, q1)[:, 4].reshape(-1, 1),
            "single: (-1,0,0) agrees with (1,0,0)", rtol=1.0e-12)
    # Gamma is a legal lattice vector, and its row is the scalar chi4(t).
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", "--qpoints", "single:0,0,0",
         "--dt", "1", "--maxframes", "8", "--lag", "1", "--s4-cutoff", s4_cutoff,
         "--s4", path("s4_single_gamma.dat"), "--chi4", path("chi4_single_gamma.dat"), "--sq", 
         path("s4_single_gamma_sq.dat")])
    compare(read_matrix(path("s4_single_gamma.dat"))[:, 4].reshape(-1, 1),
            read_matrix(path("chi4_single_gamma.dat"))[:, 2].reshape(-1, 1),
            "single: the Gamma row is chi4(t)", rtol=1.0e-12)
    # The HDF5 flavour has no header, so the indices have to be attributes:
    # |q| alone cannot tell single:1,0,0 from single:0,1,0.
    if args.h5read:
        run([exe, "dyn", "-i", diff_dump, "-w", "unit", *single_base, "--format", "hdf5",
             "--s4-cutoff", s4_cutoff, "--s4", path("s4_single.h5"), "--sq", 
             path("s4_single_h5_sq.dat")])
        attrs = {name: run([args.h5read, path("s4_single.h5"), "attr", name]).stdout.strip()
                 for name in ("sampling", "single_n1", "single_n2", "single_n3")}
        if attrs["sampling"] != "single" or \
                (attrs["single_n1"], attrs["single_n2"], attrs["single_n3"]) != ("1", "0", "0"):
            raise SystemExit("FAIL single: HDF5 attributes are %s" % attrs)
        run([exe, "dyn", "-i", diff_dump, "-w", "unit", "--qpoints", "single:0,1,0",
             "--dt", "1", "--maxframes", "8", "--lag", "1", "--format", "hdf5",
             "--s4-cutoff", s4_cutoff, "--s4", path("s4_single_y.h5"), "--sq", 
             path("s4_single_y_h5_sq.dat")])
        if run([args.h5read, path("s4_single_y.h5"), "attr", "single_n2"]).stdout.strip() != "1":
            raise SystemExit("FAIL single: the HDF5 indices do not follow the request")
        print("  ok   %-42s n = (%s %s %s)"
              % ("single: HDF5 records the indices", attrs["single_n1"], attrs["single_n2"],
                 attrs["single_n3"]))
    print("  ok   %-42s |q| = %.4f" % ("single lattice vector", q1))

    # --buffer-limit is in GB and can be lowered for this small buffer.
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", *s4_q, "--s4-cutoff", s4_cutoff,
         "--buffer-limit", "0.001", "--s4", path("s4_limit.dat"), "--sq", 
         path("s4_limit_sq.dat")])
    compare(read_matrix(path("s4_limit.dat"))[:, 4].reshape(-1, 1),
            s4[:, 4].reshape(-1, 1), "S4 with --buffer-limit 0.001", rtol=1.0e-9)

    # --stride subsamples the S4/chi4 trajectory.  maxframes is rounded
    # down to a multiple of the stride; the coherent window is unchanged.
    for stride in (2, 3):
        run([exe, "dyn", "-i", diff_dump, "-w", "unit", *s4_base, "--lag", "1",
             "--s4-cutoff", s4_cutoff, "--stride", str(stride),
             "--s4", path("stride.dat"), "--chi4", path("chi4_stride.dat"), "--sq", 
             path("stride_sq.dat")])
        run([sys.executable, ref_four, "--input", diff_dump, "--qpoints", "line:" + s4_spec,
             "--maxframes", "8", "--lag", "1", "--stride", str(stride), "--dt", "1",
             "--cutoff", s4_cutoff, "--output-s4", path("stride_ref.dat"),
             "--output-chi4", path("chi4_stride_ref.dat")])
        compare(read_matrix(path("stride.dat")), read_matrix(path("stride_ref.dat")),
                "S4 stride %d vs reference" % stride, rtol=1.0e-9)
        compare(read_matrix(path("chi4_stride.dat")), read_matrix(path("chi4_stride_ref.dat")),
                "chi4 stride %d vs reference" % stride, rtol=1.0e-9)

    s4_round = ["--qpoints", "line:" + s4_spec, "--dt", "1", "--maxframes", "9"]
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", *s4_round, "--lag", "1",
         "--s4-cutoff", s4_cutoff, "--stride", "2", "--s4", path("s4_round.dat"),
         "--chi4", path("chi4_round.dat"), "--sq", path("s4_round_sq.dat")])
    run([sys.executable, ref_four, "--input", diff_dump, "--qpoints", "line:" + s4_spec,
         "--maxframes", "9", "--lag", "1", "--stride", "2", "--dt", "1",
         "--cutoff", s4_cutoff, "--output-s4", path("s4_round_ref.dat"),
         "--output-chi4", path("chi4_round_ref.dat")])
    compare(read_matrix(path("s4_round.dat")), read_matrix(path("s4_round_ref.dat")),
            "S4 stride rounding vs reference", rtol=1.0e-9)
    compare(read_matrix(path("chi4_round.dat")), read_matrix(path("chi4_round_ref.dat")),
            "chi4 stride rounding vs reference", rtol=1.0e-9)
    if read_matrix(path("s4_round.dat"))[:, 3].max() != 8.0:
        raise SystemExit("FAIL maxframes 9 with stride 2 should use an effective window of 8")
    print("  ok   %-42s effective window 8" % "S4 stride rounding")

    # --- self intermediate scattering function (--fqt-self) ----------------
    print("self intermediate scattering function (--fqt-self)")
    ref_self = os.path.join(HERE, "ref_fsqt.py")
    fs_base = ["--qpoints", "line:" + s4_spec, "--dt", "1", "--maxframes", "8"]
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", *fs_base, "--lag", "1",
         "--fqt-self", path("fs.dat"), "--sq", path("fs_sq.dat")])
    run([sys.executable, ref_self, "--input", diff_dump, "--qpoints", "line:" + s4_spec,
         "--maxframes", "8", "--lag", "1", "--stride", "1", "--dt", "1",
         "--output", path("fs_ref.dat")])
    fs = read_matrix(path("fs.dat"))
    compare(fs, read_matrix(path("fs_ref.dat")), "F_s(q,t) vs numpy reference", rtol=1.0e-9)
    if abs(fs[0, 4] - 1.0) > 1.0e-12:
        raise SystemExit("FAIL total F_s(q,0) != 1")
    if float(np.max(np.abs(fs[:, 4] - fs[:, 5:].sum(axis=1)))) > 1.0e-12:
        raise SystemExit("FAIL species F_s columns do not sum to the total")
    print("  ok   %-42s F_s(q,0)=%.6f" % ("F_s normalization", fs[0, 4]))

    # F_s uses every atom and never the overlap cutoff.
    for cut in ("0.3", "1.5"):
        run([exe, "dyn", "-i", diff_dump, "-w", "unit", *fs_base, "--lag", "1",
             "--s4-cutoff", cut, "--s4", path("fs_cut_s4.dat"),
             "--fqt-self", path("fs_cut.dat"), "--sq", path("fs_cut_sq.dat")])
        compare(read_matrix(path("fs_cut.dat")), fs, "F_s cutoff %s invariance" % cut,
                rtol=1.0e-9)

    # F_s with stride and origin lag against the reference.
    for stride in (2, 3):
        run([exe, "dyn", "-i", diff_dump, "-w", "unit", *fs_base, "--lag", "3",
             "--stride", str(stride), "--fqt-self", path("fs_stride.dat"), "--sq", 
             path("fs_stride_sq.dat")])
        run([sys.executable, ref_self, "--input", diff_dump, "--qpoints", "line:" + s4_spec,
             "--maxframes", "8", "--lag", "3", "--stride", str(stride), "--dt", "1",
             "--output", path("fs_stride_ref.dat")])
        compare(read_matrix(path("fs_stride.dat")), read_matrix(path("fs_stride_ref.dat")),
                "F_s stride %d lag 3 vs reference" % stride, rtol=1.0e-9)

    # Diffusive ideal gas: F_s(q,t) = exp(-D q^2 t).
    tau_self = fs[:, 3].reshape(2, 9)[0]
    worst_self = 0.0
    for iq, qval in enumerate((2.0, 3.0)):
        block = fs[iq*9:(iq+1)*9, 4]
        analytic = np.exp(-0.1 * qval**2 * tau_self)
        worst_self = max(worst_self, float(np.max(np.abs(block - analytic))))
    if worst_self > 0.01:
        raise SystemExit("FAIL self F_s deviates from exp(-D q2 t) by %.3f" % worst_self)
    print("  ok   %-42s max deviation %.3f" % ("diffusive F_s vs exp(-D q2 t)", worst_self))

    # --- mean squared displacement (--msd) ---------------------------------
    print("mean squared displacement (--msd)")
    ref_msd = os.path.join(HERE, "ref_msd.py")
    msd_base = ["--qpoints", "line:" + s4_spec, "--dt", "1", "--maxframes", "8"]
    run([exe, "dyn", "-i", diff_dump, "-m", "1:Si,2:O", *msd_base, "--lag", "1",
         "--msd", path("msd.dat"), "--sq", path("msd_sq.dat")])
    run([sys.executable, ref_msd, "--input", diff_dump, "--maxframes", "8",
         "--lag", "1", "--stride", "1", "--dt", "1", "--output", path("msd_ref.dat")])
    msd = read_matrix(path("msd.dat"))
    compare(msd, read_matrix(path("msd_ref.dat")), "MSD(t) vs numpy reference", rtol=1.0e-9)
    if float(np.max(np.abs(msd[0, 1:]))) > 1.0e-12:
        raise SystemExit("FAIL MSD(0) != 0")
    # the text columns carry 13 significant digits, so the sum rule holds to a
    # scale-following tolerance (MSD grows to ~5 in this window)
    if float(np.max(np.abs(msd[:, 1] - msd[:, 2:].sum(axis=1)))) > 1.0e-10:
        raise SystemExit("FAIL species MSD columns do not sum to the total")
    print("  ok   %-42s MSD(0)=0" % "MSD normalization")

    # The diffusive ideal gas has MSD = 6 D t, so the slope is 6*0.1 = 0.6 and
    # each species column carries the same slope in this equilibrated mixture.
    tau_msd = msd[:, 0]
    slope_msd = (msd[-1, 1] - msd[0, 1])/(tau_msd[-1] - tau_msd[0])
    if abs(slope_msd - 0.6) > 0.02:
        raise SystemExit("FAIL diffusive MSD slope %.4f != 6D = 0.6" % slope_msd)
    print("  ok   %-42s slope %.4f (6D = 0.6)" % ("diffusive MSD slope", slope_msd))

    # MSD is unit weighted: -w/--norm, the overlap cutoff and every other
    # output that shares the position buffer leave it untouched.
    run([exe, "dyn", "-i", diff_dump, "-m", "1:Si,2:O", *msd_base, "--lag", "1",
         "--s4-cutoff", "0.5", "--s4", path("msd_s4.dat"), "--chi4", path("msd_chi4.dat"),
         "--fqt-self", path("msd_fs.dat"), "--sqw", path("msd_sqw.dat"),
         "--msd", path("msd_all.dat"), "--sq", path("msd_all_sq.dat")])
    compare(read_matrix(path("msd_all.dat")), msd, "MSD with S4/chi4/F_s/sqw", rtol=1.0e-9)
    run([exe, "dyn", "-i", diff_dump, "-m", "1:Si,2:O", "-w", "unit", "--norm", "n",
         *msd_base, "--lag", "1", "--msd", path("msd_n.dat"), "--sq", path("msd_n_sq.dat")])
    compare(read_matrix(path("msd_n.dat")), msd, "MSD --norm n invariance", rtol=1.0e-9)

    # stride and origin lag against the reference.
    for stride in (2, 3):
        run([exe, "dyn", "-i", diff_dump, "-m", "1:Si,2:O", *msd_base, "--lag", "3",
             "--stride", str(stride), "--msd", path("msd_stride.dat"), "--sq", 
             path("msd_stride_sq.dat")])
        run([sys.executable, ref_msd, "--input", diff_dump, "--maxframes", "8",
             "--lag", "3", "--stride", str(stride), "--dt", "1",
             "--output", path("msd_stride_ref.dat")])
        compare(read_matrix(path("msd_stride.dat")), read_matrix(path("msd_stride_ref.dat")),
                "MSD stride %d lag 3 vs reference" % stride, rtol=1.0e-9)

    # Without --qpoints nothing is sampled, so there is no S(q) table at all and the
    # MSD is the only output; --no-partials drops the species columns.
    run([exe, "dyn", "-i", diff_dump, "-m", "1:Si,2:O", "--dt", "1",
         "--maxframes", "8", "--msd", path("msd_none.dat")])
    compare(read_matrix(path("msd_none.dat")), msd, "MSD without q points", rtol=1.0e-9)
    run([exe, "dyn", "-i", diff_dump, "-m", "1:Si,2:O", "--no-partials", *msd_base, "--lag", "1",
         "--msd", path("msd_total.dat"), "--sq", path("msd_total_sq.dat")])
    if read_matrix(path("msd_total.dat")).shape[1] != 2:
        raise SystemExit("FAIL --no-partials still wrote MSD species columns")
    print("  ok   %-42s total only" % "--no-partials MSD columns")

    # The origin stride is shared with the coherent correlations.
    s4_lag = s4_base + ["--lag", "3"]
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", *s4_lag, "--s4-cutoff", s4_cutoff,
         "--s4", path("s4_lag.dat"), "--chi4", path("chi4_lag.dat"), "--sq", 
         path("s4_lag_sq.dat")])
    run([sys.executable, ref_four, "--input", diff_dump, "--qpoints", "line:" + s4_spec,
         "--maxframes", "8", "--lag", "3", "--dt", "1", "--cutoff", s4_cutoff,
         "--output-s4", path("s4_lag_ref.dat"), "--output-chi4", path("chi4_lag_ref.dat")])
    compare(read_matrix(path("s4_lag.dat")), read_matrix(path("s4_lag_ref.dat")),
            "S4 lag stride 3 vs reference", rtol=1.0e-9)
    compare(read_matrix(path("chi4_lag.dat")), read_matrix(path("chi4_lag_ref.dat")),
            "chi4 lag stride 3 vs reference", rtol=1.0e-9)

    # Ideal-gas diffusive limit: for independent Brownian particles and q
    # away from the reciprocal lattice, S4(q,t) -> p(t), Q(t) -> p(t) and
    # chi4(t) -> p(t)(1-p(t)), with p the 3D Gaussian overlap probability.
    diffusivity = 0.1
    cutoff = float(s4_cutoff)
    tau = chi4[:, 0]

    def overlap_probability(time):
        if time <= 0.0:
            return 1.0
        sigma = math.sqrt(diffusivity*time)
        return (math.erf(cutoff/(2.0*sigma))
                - cutoff/(math.sqrt(math.pi)*sigma)
                * math.exp(-cutoff*cutoff/(4.0*diffusivity*time)))

    p = np.array([overlap_probability(t) for t in tau])
    worst_q = max(float(np.max(np.abs(s4[1:9, 4] - p[1:]))),
                  float(np.max(np.abs(s4[10:18, 4] - p[1:]))))
    worst_q0 = float(np.max(np.abs(chi4[1:, 1] - p[1:])))
    worst_c = float(np.max(np.abs(chi4[1:, 2] - p[1:]*(1.0 - p[1:]))))
    if worst_q > 0.02 or worst_q0 > 0.01 or worst_c > 0.015:
        raise SystemExit("FAIL ideal gas diffusive limit: S4 %.3f Q %.3f chi4 %.3f"
                         % (worst_q, worst_q0, worst_c))
    print("  ok   %-42s max deviations S4 %.3f Q %.3f chi4 %.3f"
          % ("ideal gas diffusive limit", worst_q, worst_q0, worst_c))

    # A non-periodic direction must not be wrapped: the dump writes the true
    # coordinate for it (x, not xu), so an atom that leaves the box has to be
    # kept outside it.  Both spellings of the same slab trajectory must agree.
    slab_kwargs = dict(natoms=400, length=10, frames=25, seed=9, mode="ballistic",
                       temperature=0.5, boundary="pp pp ff")
    slab_xu = generate("dyn_slab_xu.dump", columns="id type xu yu zu", **slab_kwargs)
    slab_x = generate("dyn_slab_x.dump", columns="id type x y z", **slab_kwargs)
    slab_q = dyn_line("2,1,3,0,0,1") + ["--dt", "1", "--maxframes", "5"]
    for dump_file, tag in ((slab_xu, "xu"), (slab_x, "x")):
        run([exe, "dyn", "-i", dump_file, "-m", "1:Si,2:O", "-w", "unit", "--norm", "mean",
             *slab_q, "--fqt", path("dyn_slab_%s_fqt.dat" % tag),
             "--sq", path("dyn_slab_%s_sq.dat" % tag)])
    compare(read_matrix(path("dyn_slab_x_fqt.dat")),
            read_matrix(path("dyn_slab_xu_fqt.dat")),
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
        run([exe, "dyn", "-i", dyn_dump, "-m", "1:Si,2:O", "-w", "unit", "--norm", "mean", *dyn_q,
             "--sqw", path("dyn_sqw.h5"), "--fqt", path("dyn_fqt.h5"), "--sq", path("dyn_sq_h5.dat")])
        with open(path("dyn_sqw_h5.txt"), "w") as handle:
            handle.write(run([args.h5read, path("dyn_sqw.h5"), "sqw"]).stdout)
        with open(path("dyn_fqt_h5.txt"), "w") as handle:
            handle.write(run([args.h5read, path("dyn_fqt.h5"), "fqt"]).stdout)
        compare(read_matrix(path("dyn_sqw_h5.txt")), sqw, "HDF5 S(q,w) vs text", rtol=1.0e-9)
        compare(read_matrix(path("dyn_fqt_h5.txt")), fqt, "HDF5 F(q,t) vs text", rtol=1.0e-9)
        run([exe, "dyn", "-i", diff_dump, "-w", "unit", *s4_q, "--s4-cutoff", s4_cutoff,
             "--s4", path("s4.h5"), "--chi4", path("chi4.h5"), "--sq", path("s4_h5_sq.dat")])
        with open(path("s4_h5.txt"), "w") as handle:
            handle.write(run([args.h5read, path("s4.h5"), "s4"]).stdout)
        with open(path("chi4_h5.txt"), "w") as handle:
            handle.write(run([args.h5read, path("chi4.h5"), "chi4"]).stdout)
        compare(read_matrix(path("s4_h5.txt")), s4, "HDF5 S4(q,t) vs text", rtol=1.0e-9)
        compare(read_matrix(path("chi4_h5.txt")), chi4, "HDF5 chi4 vs text", rtol=1.0e-9)
        run([exe, "dyn", "-i", diff_dump, "-w", "unit", *fs_base, "--lag", "1",
             "--fqt-self", path("fs.h5"), "--sq", path("fs_h5_sq.dat")])
        with open(path("fs_h5.txt"), "w") as handle:
            handle.write(run([args.h5read, path("fs.h5"), "fqt_self"]).stdout)
        compare(read_matrix(path("fs_h5.txt")), fs, "HDF5 F_s(q,t) vs text", rtol=1.0e-9)
        run([exe, "dyn", "-i", diff_dump, "-m", "1:Si,2:O", *fs_base, "--lag", "1",
             "--msd", path("msd.h5"), "--sq", path("msd_h5_sq.dat")])
        with open(path("msd_h5.txt"), "w") as handle:
            handle.write(run([args.h5read, path("msd.h5"), "msd"]).stdout)
        compare(read_matrix(path("msd_h5.txt")), msd, "HDF5 MSD vs text", rtol=1.0e-9)

    # --- the q sampling modes of --qpoints --------------------------------------
    print("--qpoints sampling modes")

    # Without q points only the scalar overlap is left: chi4 must come out
    # exactly as it does on a q line, and there is no S(q) table at all.
    no_q = ["--dt", "1", "--maxframes", "8", "--lag", "1"]
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", *no_q, "--s4-cutoff", s4_cutoff,
         "--chi4", path("noq_chi4.dat")])
    compare(read_matrix(path("noq_chi4.dat")), chi4, "chi4 without q points", rtol=1.0e-12)
    result = subprocess.run([exe, "dyn", "-i", diff_dump, "-w", "unit", *no_q,
                             "--s4-cutoff", s4_cutoff, "--chi4", path("noq_c.dat"), "-"],
                            capture_output=True, text=True)
    if result.returncode == 0:
        raise SystemExit("FAIL dyn accepted a positional S(q) table")
    print("  ok   %-42s rejected" % "dyn with a positional argument")

    # A shell is one |q| averaged over every direction: all the tables hold a
    # single q row, whose value is the Lebedev average.
    shell_spec = "shell:2.5,medium"
    shell_base = ["--qpoints", shell_spec, "--dt", "1", "--maxframes", "8"]
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", *shell_base, "--lag", "1",
         "--s4-cutoff", s4_cutoff, "--sqw", path("shell_sqw.dat"),
         "--fqt", path("shell_fqt.dat"), "--fqt-self", path("shell_fs.dat"),
         "--s4", path("shell_s4.dat"), "--chi4", path("shell_chi4.dat"), "--sq", 
         path("shell_sq.dat")])
    shell_q, shell_s = read_table(path("shell_sq.dat"))
    shell_fqt = read_matrix(path("shell_fqt.dat"))
    shell_sqw = read_matrix(path("shell_sqw.dat"))
    shell_fs = read_matrix(path("shell_fs.dat"))
    if shell_q.shape != (1,) or abs(shell_q[0] - 2.5) > 1.0e-9:
        raise SystemExit("FAIL a shell must write one S(q) row at |q| = 2.5")
    if shell_fqt.shape[0] != 9 or shell_sqw.shape[0] != 9 or shell_fs.shape[0] != 9:
        raise SystemExit("FAIL the shell tables must hold a single q row")
    compare(shell_s.reshape(-1, 1), shell_fqt[0:1, 4], "shell F(q,0) = static S(q)",
            rtol=1.0e-9)
    compare(shell_s.reshape(-1, 1),
            shell_sqw[:, 4].reshape(1, 9).sum(axis=1, keepdims=True)*(np.pi/8.0),
            "shell sum rule int S(q,w) dw = S(q)", rtol=1.0e-9)
    worst = float(np.max(np.abs(shell_fqt[:, 4]
                                - (shell_fqt[:, 5] + 2.0*shell_fqt[:, 6] + shell_fqt[:, 7]))))
    if worst > 1.0e-9:
        raise SystemExit("FAIL shell partial sum rule is off by %.3e" % worst)
    # chi4 is the q -> 0 scalar: the shell must not touch it.
    compare(read_matrix(path("shell_chi4.dat")), chi4, "chi4 is q independent",
            rtol=1.0e-12)
    if abs(shell_fs[0, 4] - 1.0) > 1.0e-12:
        raise SystemExit("FAIL shell F_s(q,0) != 1")
    print("  ok   %-42s S(q)=%.6f" % ("shell row and sum rules", shell_s[0]))
    if args.h5read:
        run([exe, "dyn", "-i", diff_dump, "-w", "unit", *shell_base, "--lag", "1",
             "--s4-cutoff", s4_cutoff, "--sqw", path("shell_sqw.h5"),
             "--fqt-self", path("shell_fs.h5"), "--s4", path("shell_s4.h5"), "--sq", 
             path("shell_h5_sq.dat")])
        with open(path("shell_sqw_h5.txt"), "w") as handle:
            handle.write(run([args.h5read, path("shell_sqw.h5"), "sqw"]).stdout)
        with open(path("shell_fs_h5.txt"), "w") as handle:
            handle.write(run([args.h5read, path("shell_fs.h5"), "fqt_self"]).stdout)
        compare(read_matrix(path("shell_sqw_h5.txt")), shell_sqw,
                "HDF5 shell S(q,w) vs text", rtol=1.0e-9)
        compare(read_matrix(path("shell_fs_h5.txt")), shell_fs,
                "HDF5 shell F_s(q,t) vs text", rtol=1.0e-9)
        # The per-lag origin counts have to survive the shell collapse: the
        # 1 x nlag /sqw/count must equal the /s4/count of the same run, not
        # the lag-0 count repeated.
        with open(path("shell_sqw_count.txt"), "w") as handle:
            handle.write(run([args.h5read, path("shell_sqw.h5"), "count", "sqw"]).stdout)
        with open(path("shell_s4_count.txt"), "w") as handle:
            handle.write(run([args.h5read, path("shell_s4.h5"), "count", "s4"]).stdout)
        count_sqw = read_matrix(path("shell_sqw_count.txt"))
        count_s4 = read_matrix(path("shell_s4_count.txt"))
        if count_sqw.shape != count_s4.shape or count_sqw.shape[0] != 1:
            raise SystemExit("FAIL shell origin counts have shapes %s and %s"
                             % (count_sqw.shape, count_s4.shape))
        compare(count_sqw, count_s4, "HDF5 shell per-lag origin counts", rtol=0.0)
    if os.path.isfile(plot_dyn) and have_matplotlib:
        figure = path("plot_shell.png")
        run([sys.executable, plot_dyn, "--mode", "spectra", path("shell_sqw.dat"),
             "-o", figure])
        if not os.path.isfile(figure) or os.path.getsize(figure) == 0:
            raise SystemExit("FAIL plot_sqw.py wrote no figure for a shell table")
        print("  ok   %-42s single q row" % "plot_sqw.py renders a shell")

    try:
        from scipy.integrate import lebedev_rule  # noqa: F401
        have_scipy = True
    except ImportError:
        have_scipy = False
    if have_scipy:
        run([sys.executable, ref_dyn, "--input", diff_dump, "--qpoints", shell_spec,
             "--maxframes", "8", "--lag", "1", "--dt", "1", "--weight", "unit",
             "--norm", "mean", "--output-fqt", path("shell_ref_fqt.dat"),
             "--output-sqw", path("shell_ref_sqw.dat")])
        compare(shell_fqt, read_matrix(path("shell_ref_fqt.dat")),
                "shell F(q,t) vs scipy reference", rtol=1.0e-9)
        compare(shell_sqw, read_matrix(path("shell_ref_sqw.dat")),
                "shell S(q,w) vs scipy reference", rtol=1.0e-9)
        run([sys.executable, ref_four, "--input", diff_dump, "--qpoints", shell_spec,
             "--maxframes", "8", "--lag", "1", "--dt", "1", "--cutoff", s4_cutoff,
             "--output-s4", path("shell_ref_s4.dat")])
        compare(read_matrix(path("shell_s4.dat")), read_matrix(path("shell_ref_s4.dat")),
                "shell S4(q,t) vs scipy reference", rtol=1.0e-9)
        run([sys.executable, ref_self, "--input", diff_dump, "--qpoints", shell_spec,
             "--maxframes", "8", "--lag", "1", "--stride", "1", "--dt", "1",
             "--output", path("shell_ref_fs.dat")])
        compare(shell_fs, read_matrix(path("shell_ref_fs.dat")),
                "shell F_s(q,t) vs scipy reference", rtol=1.0e-9)
    else:
        print("  skip %-42s scipy is not installed" % "shell sampling vs reference")

    # --- a q line long enough for the separable phase evaluation ----------
    # The scale of a line steps by a constant, so its phase is a
    # mode-independent start s0 times the n-th power of one factor table - the
    # same identity the grid uses, but on a single axis.  A line switches to
    # it from phase_min_modes = 16 points on, which is longer than every other
    # line in this suite, so this is the only coverage of the line tables.
    print("separable q line (--qpoints line, 16 points or more)")
    long_spec = "20,2,3,1,0,0"         # 21 q points from 2 to 3 along x
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", "--qpoints", "line:" + long_spec,
         "--dt", "1", "--maxframes", "8", "--lag", "1", "--s4-cutoff", s4_cutoff,
         "--s4", path("s4_long.dat"), "--chi4", path("chi4_long.dat"),
         "--fqt", path("fq_long.dat"), "--fqt-self", path("fs_long.dat"),
         "--sqw", path("sqw_long.dat"), "--sq", path("s4_long_sq.dat")])
    # The references sum the line directly, so they pin the factor tables and
    # the base term of the first axis independently of the recurrence.
    run([sys.executable, ref_four, "--input", diff_dump, "--qpoints", "line:" + long_spec,
         "--maxframes", "8", "--lag", "1", "--dt", "1", "--cutoff", s4_cutoff,
         "--output-s4", path("s4_long_ref.dat")])
    run([sys.executable, ref_self, "--input", diff_dump, "--qpoints", "line:" + long_spec,
         "--maxframes", "8", "--lag", "1", "--stride", "1", "--dt", "1",
         "--output", path("fs_long_ref.dat")])
    run([sys.executable, ref_dyn, "--input", diff_dump, "--qpoints", "line:" + long_spec,
         "--maxframes", "8", "--lag", "1", "--dt", "1", "--weight", "unit",
         "--norm", "mean", "--output-fqt", path("fq_long_ref.dat"),
         "--output-sqw", path("sqw_long_ref.dat")])
    compare(read_matrix(path("s4_long.dat")), read_matrix(path("s4_long_ref.dat")),
            "separable line: S4(q,t) vs reference", rtol=1.0e-9)
    compare(read_matrix(path("fs_long.dat")), read_matrix(path("fs_long_ref.dat")),
            "separable line: F_s(q,t) vs reference", rtol=1.0e-9)
    compare(read_matrix(path("fq_long.dat")), read_matrix(path("fq_long_ref.dat")),
            "separable line: F(q,t) vs reference", rtol=1.0e-9)
    compare(read_matrix(path("sqw_long.dat")), read_matrix(path("sqw_long_ref.dat")),
            "separable line: S(q,w) vs reference", rtol=1.0e-9)
    # The two ends of the long line are the two points of the short line the
    # rest of the suite uses, and that one stays on the direct branch: the two
    # evaluations have to meet there.
    long_s4 = read_matrix(path("s4_long.dat"))
    for q_value in (2.0, 3.0):
        compare(rows_at(read_matrix(path("s4.dat")), [q_value, 0.0, 0.0])[:, 4].reshape(-1, 1),
                rows_at(long_s4, [q_value, 0.0, 0.0])[:, 4].reshape(-1, 1),
                "separable line: direct vs separable at |q| = %g" % q_value, rtol=1.0e-11)
    # A line through Gamma turns the base term off; it has to agree too.
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", "--qpoints", "line:20,0,4,1,0,0",
         "--dt", "1", "--maxframes", "8", "--lag", "1", "--s4-cutoff", s4_cutoff,
         "--s4", path("s4_long0.dat"), "--sq", path("s4_long0_sq.dat")])
    run([sys.executable, ref_four, "--input", diff_dump, "--qpoints", "line:20,0,4,1,0,0",
         "--maxframes", "8", "--lag", "1", "--dt", "1", "--cutoff", s4_cutoff,
         "--output-s4", path("s4_long0_ref.dat")])
    compare(read_matrix(path("s4_long0.dat")), read_matrix(path("s4_long0_ref.dat")),
            "separable line: s0 = 0 vs reference", rtol=1.0e-9)
    print("  ok   %-42s 21 points, base on and off" % "separable q line")

    # The short alias and the long spelling name the same option.
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", "-q", "line:2,0.5,2,1,0,0",
         "--dt", "1", "--maxframes", "8", "--lag", "1", "--sq", path("alias_short.dat")])
    run([exe, "dyn", "-i", diff_dump, "-w", "unit", "--qpoints", "line:2,0.5,2,1,0,0",
         "--dt", "1", "--maxframes", "8", "--lag", "1", "--sq", path("alias_long.dat")])
    compare(read_matrix(path("alias_short.dat")), read_matrix(path("alias_long.dat")),
            "-q is an alias of --qpoints", rtol=0.0)

    # option validation: a flag of the other subcommand, and a dyn run that
    # cannot mean anything, both have to be refused.
    bad_static = [
        (["--qpoints", "line:4,1,4,1,0,0", "--dt", "1", "--maxframes", "8"],
         "--qpoints in static"),
        (["-q", "line:4,1,4,1,0,0"], "-q in static"),
        (["--sqw", "w.dat"], "--sqw in static"),
        (["--fqt-self", "f.dat"], "--fqt-self in static"),
        (["--msd", "m.dat"], "--msd in static"),
        (["--s4", "s.dat", "--s4-cutoff", "0.5"], "--s4 in static"),
        (["--pair-entropy", "s.dat"], "--pair-entropy without --method debye"),
        (["--method", "debye", "--s2-accum", "a.dat"],
         "--s2-accum without --pair-entropy"),
    ]
    bad_dyn = [
        (dyn_line(dyn_spec) + ["--maxframes", "20"], "dyn without --dt"),
        (dyn_line(dyn_spec) + ["--dt", "1"], "dyn without --maxframes"),
        (dyn_line(dyn_spec) + ["--dt", "1", "--maxframes", "20", "--grid", "g.dat"],
         "--grid in dyn"),
        (dyn_line(dyn_spec) + ["--dt", "1", "--maxframes", "20", "--rmax", "4"],
         "--rmax in dyn"),
        (dyn_line("0,1,4,1,0,0") + ["--dt", "1", "--maxframes", "20"],
         "--qpoints line without intervals"),
        (dyn_line("4,4,1,1,0,0") + ["--dt", "1", "--maxframes", "20"],
         "--qpoints line with S1 <= S0"),
        (dyn_line(dyn_spec) + ["--dt", "1", "--maxframes", "20", "--s4", "s.dat"],
         "--s4 without --s4-cutoff"),
        (dyn_line(dyn_spec) + ["--dt", "1", "--maxframes", "20", "--s4-cutoff", "0.5"],
         "--s4-cutoff without --s4/--chi4"),
        (dyn_line(dyn_spec) + ["--dt", "1", "--maxframes", "20", "--s4-cutoff", "0.5",
          "--format", "bogus", "--s4", "s.dat"], "--format bogus"),
        (dyn_line(dyn_spec) + ["--dt", "1", "--maxframes", "100000000",
          "--s4-cutoff", "0.5", "--s4", "s.dat"], "S4 position buffer too large"),
        (dyn_line(dyn_spec) + ["--dt", "1", "--maxframes", "8", "--s4-cutoff", "0.5",
          "--buffer-limit", "0.0001", "--s4", "s.dat"], "S4 buffer limit too small"),
        (dyn_line(dyn_spec) + ["--dt", "1", "--maxframes", "8", "--buffer-limit", "2.0"],
         "--buffer-limit without S4/chi4/F_s"),
        (dyn_line(dyn_spec) + ["--dt", "1", "--maxframes", "8", "--s4-cutoff", "0.5",
          "--buffer-limit", "0", "--s4", "s.dat"], "--buffer-limit zero"),
        (dyn_line(dyn_spec) + ["--dt", "1", "--maxframes", "8", "--s4-cutoff", "0.5",
          "--stride", "9", "--s4", "s.dat"], "--stride exceeds --maxframes"),
        (dyn_line(dyn_spec) + ["--dt", "1", "--maxframes", "8", "--s4-cutoff", "0.5",
          "--stride", "0", "--s4", "s.dat"], "--stride zero"),
        (dyn_line(dyn_spec) + ["--dt", "1", "--maxframes", "8", "--stride", "2"],
         "--stride without S4/chi4/F_s"),
        (dyn_line(dyn_spec) + ["--dt", "1", "--maxframes", "8", "--s4-cutoff", "0.5",
          "--fqt-self", "f.dat"], "--s4-cutoff with --fqt-self but no S4/chi4"),
        (["--qpoints", "bogus", "--dt", "1", "--maxframes", "8"],
         "--qpoints with an unknown spec"),
        (["--qpoints", "shell:2.5,bogus", "--dt", "1", "--maxframes", "8"],
         "--qpoints shell with an unknown accuracy"),
        (["--qpoints", "shell:0,low", "--dt", "1", "--maxframes", "8"],
         "--qpoints shell with a zero radius"),
        (["--qpoints", "line:4,1,4,1,0,0,7", "--dt", "1", "--maxframes", "8"],
         "--qpoints line with too many fields"),
        (["--dt", "1", "--maxframes", "8", "--sqw", "w.dat"],
         "--sqw without a q sampling"),
        (["--dt", "1", "--maxframes", "8", "--s4-cutoff", "0.5", "--s4", "s.dat"],
         "--s4 without a q sampling"),
        (["--qpoints", "-", "--dt", "1", "--maxframes", "8"],
         "--qpoints - is not a spec"),
        (["--dt", "1", "--maxframes", "8"],
         "dyn without a q sampling and without --chi4/--msd"),
        (["--qpoints", "grid:0", "--dt", "1", "--maxframes", "8"],
         "--qpoints grid with a zero qmax"),
        (["--qpoints", "powder:2.5,80,dq=0.03", "--dt", "1", "--maxframes", "8"],
         "powder with both M and dq"),
        (["--qpoints", "powder:0", "--dt", "1", "--maxframes", "8"],
         "powder with a zero radius"),
        (["--qpoints", "line:4,1,4,1,0,0", "--modes", "4",
          "--dt", "1", "--maxframes", "8"],
         "--modes with a q line"),
        (["--qpoints", "line:4,1,4,1,0,0", "--thin", "orbits",
          "--dt", "1", "--maxframes", "8"],
         "--thin with a q line"),
        (["--qpoints", "grid:1", "--thin", "diagonal",
          "--dt", "1", "--maxframes", "8"],
         "--thin with an unknown policy"),
        (["--qpoints", "single:1,0", "--dt", "1", "--maxframes", "8"],
         "--qpoints single with two indices"),
        (["--qpoints", "single:1,0,0,0", "--dt", "1", "--maxframes", "8"],
         "--qpoints single with four indices"),
        (["--qpoints", "single:1.5,0,0", "--dt", "1", "--maxframes", "8"],
         "--qpoints single with a non-integer index"),
        (["--qpoints", "single:1,0,0", "--modes", "4",
          "--dt", "1", "--maxframes", "8"],
         "--modes with a single q"),
    ]
    rejected = 0
    for options, label in bad_static:
        result = subprocess.run([exe, "static", "-i", dyn_dump] + options + ["-"],
                                capture_output=True, text=True)
        if result.returncode == 0:
            raise SystemExit("FAIL %s was accepted" % label)
        rejected += 1
    for options, label in bad_dyn:
        result = subprocess.run([exe, "dyn", "-i", dyn_dump] + options,
                                capture_output=True, text=True)
        if result.returncode == 0:
            raise SystemExit("FAIL %s was accepted" % label)
        rejected += 1
    print("  ok   %-42s %d cases rejected"
          % ("subcommand and option validation", rejected))

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
        result = subprocess.run([exe, "dyn", "-i", broken, "-w", "unit", *dyn_line(dyn_spec),
                                 "--dt", "1", "--maxframes", "20"],
                                capture_output=True, text=True)
        if result.returncode == 0:
            raise SystemExit("FAIL %s was accepted" % label)
        print("  ok   %-42s rejected" % label)

    print("all sqcalc tests passed")


if __name__ == "__main__":
    main()
