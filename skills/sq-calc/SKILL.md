---
name: sq-calc
description: Build, run, and troubleshoot sqcalc, the Fortran/CMake structure factor program that computes S(q) and g(r) from LAMMPS dump trajectories. Use for its CMake configuration, command line options, output formats, and tests.
metadata:
  short-description: Build and use the sqcalc structure factor code
---

# sqcalc

`sqcalc` reads a LAMMPS `atoms` dump trajectory, averages over the frames, and
writes the total structure factor $S(q)$ plus one partial column per LAMMPS
type pair as a `# q S(q)` table.  Two independent routes are available: a
reciprocal space non-uniform FFT (FINUFFT on the CPU, cufinufft/cuFFT on a GPU)
and a real space Debye pair histogram, which can also write $g(r)$.

With `--dyn` it instead keeps the time axis and writes the dynamic structure
factor $S(q,\omega)$ along a line in reciprocal space (`--sqw`), optionally
with the intermediate scattering function $F(q,t)$ (`--fsq`).  The same run
can write the total four-point structure factor $S_4(q,t)$ (`--s4`) and the
average overlap and dynamic susceptibility $Q(t)$, $\chi_4(t)$ (`--chi4`);
both need the overlap cutoff `--s4-cutoff`.

It is one self-contained binary with no data files, plugins or environment
variables, so it belongs on `PATH` and is called as `sqcalc`; `command -v
sqcalc` is the check that it resolves.  Nothing else has to be located to run a
calculation:

```sh
sqcalc -i traj.dump -m 1:Si,2:O -w neutron -t 8 S_q.dat
```

## Documents

Both ship with the source and are checked against it, so read the relevant one
instead of reconstructing the details from `--help` or a build tree.

* [references/usage.md](references/usage.md) - every command line option, what
  the weighting schemes, normalizations, methods and partials mean, the output
  formats (text, HDF5, `g(r)`) and worked examples.  Read it before choosing
  flags or explaining a result.
* [references/building.md](references/building.md) - requirements, the CMake
  configure/build/test commands, the cache options (FFTW, CUDA, HDF5, vendored
  versus system FINUFFT) and the test suite.  Read it when `sqcalc` is missing
  or the task is to change the code.

`scripts/plot_sq.py` turns the text tables into a figure.  Pass the `S(q)` table,
a Debye `g(r)` table or a reciprocal grid table and it labels the columns from
the header line and draws the total plus the partials:

```sh
python3 scripts/plot_sq.py S_q.dat                     # -> S_q.png
python3 scripts/plot_sq.py --output g.png --show g.dat # interactive window
python3 scripts/plot_sq.py -o compare.png a.dat b.dat  # overlay two totals
```

It needs matplotlib and nothing else; `--output` picks the format from the
suffix, so `.pdf` or `.svg` gives vector output.  A dashed line marks the 1 that
S(q), g(r) and the Faber-Ziman partials tend to; `--refline 0` or
`--refline none` changes or removes it.

`scripts/plot_sqw.py` plots a `--sqw`/`--fsq`/`--s4` table instead: the map of
the quantity over $(|q|,\omega)$ (or $(|q|,\tau)$) next to a few selected
spectra, or single curves with `--mode spectra --q 1,4`:

```sh
python3 scripts/plot_sqw.py S_qw.dat
python3 scripts/plot_sqw.py --mode spectra --q 1,4 -o spectra.pdf F_qt.dat
```

`scripts/choose_q.py DUMP` picks the q sampling from the box, which the
reciprocal methods need because their q grid is the lattice of that box:

```sh
python3 scripts/choose_q.py traj.dump --qmax 15            # report
python3 scripts/choose_q.py traj.dump --qmax 15 --print-options
sqcalc -i traj.dump $(python3 scripts/choose_q.py traj.dump --qmax 15 --print-options) S_q.dat
```

It needs numpy, reads only the first frame, and prints `--qmin/--qmax/--nq`
such that no shell is empty, with the mode counts that show where the
statistics are thin.  Pass `--method debye` for the real-space method, whose q
grid is not tied to the box: there a fine dq costs nothing and adds no noise,
and `--dr` instead of the box sets the usable q range.

## Reminders

* `-m` (LAMMPS type id to element) is required for the `neutron`/`xray` weights;
  without it the partials are still written, labelled by LAMMPS type id
  (`S(1-1)`, `S(1-2)`, ...).
* The reciprocal methods sample $q$ on the lattice of the dump box, so the box
  must be constant and the grid costs $(q_{\max} L)^3$; the Debye method needs
  no periodic box and wins for large, sparse systems.
* Choose that sampling from the box: `--qmin` is the smallest accessible $|q|$
  ($2\pi/L$ for a cubic box) and $\Delta q = (q_{\max}-q_{\min})/n_q$ must not
  undercut the spacing of the lattice, or shells come out empty and sqcalc
  writes them as 0.  `scripts/choose_q.py` prints a safe `--qmin`, `--qmax`
  and `--nq`; for `--method debye` it reports the opposite advice, since that
  q grid is free and `--dr` limits the range instead.
* Report partials with `-fz` (Faber-Ziman) in preference to the default OVITO
  columns: they tend to 1 for every pair, so they do not hide the structure
  behind the concentrations and compare directly between systems.
* `--s4` and `--chi4` are total overlap quantities and use unit weights, so
  `-w`/`--norm` and `--partials` do not change them.  They need a physically
  meaningful `--s4-cutoff` (a fraction of the particle diameter), and the
  full $S_4$ calculation is the most expensive dynamic output; use a large
  `--lag` and a short low-$q$ line when needed.  The position buffer is
  limited to 2 GB by default (`--s4-buffer-limit GB` raises it), and a run
  with only `--s4`/`--chi4` skips the coherent $F(q,t)/S(q,\omega)$ buffers.
* The table goes to `OUTPUT` and progress to stderr, so `-q` keeps logs clean.
* The tables are plain text, so `scripts/plot_sq.py` (or any column reader)
  plots them without further tooling.

The formulas behind all of this are in the repository `README.md`.
