#!/usr/bin/env python3
"""End to end checks for the sqcalc executable.

The tests combine three independent views of the same physics:

1. sqcalc with the NUFFT method (the production path),
2. sqcalc with the direct summation method,
3. a small numpy reference implementation (test/ref_sq.py),

plus two qualitative checks (ideal gas -> S ~ 1, cubic lattice -> Bragg peaks)
and a handful of command line error cases.
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

    print("all sqcalc tests passed")


if __name__ == "__main__":
    main()
