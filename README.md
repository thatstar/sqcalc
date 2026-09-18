# sqcalc

`sqcalc` computes the total structure factor $S(q)$ of a LAMMPS trajectory.  It
reads the default LAMMPS `atoms` dump style, averages over all snapshots and
writes a `# q S(q)` table.  Two independent evaluations are available: the
reciprocal space transform with a non-uniform FFT (FINUFFT) on the CPU, or on a
CUDA GPU with cuFFT, which averages $|\rho(q)|^2$, and the real space Debye
pair-histogram method.  With an element mapping the table also carries the
partial structure factors $S_{ab}(q)$.

```
sqcalc -i traj.dump -m 1:Si,2:O -w neutron -t 8 S_q.dat
```

## Building

Requirements

* a Fortran 2008 compiler (gfortran >= 10), CMake >= 3.20 and OpenMP
* an FFTW3 installation (library and headers); CMake looks in `FFTW_ROOT`,
  `$CONDA_PREFIX`, any `~/.pixi/envs/*` and the usual system prefixes
* FINUFFT, which is vendored in `external/finufft` (git subtree, v2.5.1) and
  built together with sqcalc, so nothing has to be downloaded by hand

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build          # optional, runs the test suite
```

If FFTW lives in a non standard prefix pass `-DFFTW_ROOT=/path/to/fftw`.  The
first configure of the vendored FINUFFT obtains its build helpers
(`CPM.cmake`, `findFFTW`, `xsimd` and, unless a system FFTW is found, FFTW3)
over the network; they are cached in `.cpm-cache/` next to the sources, so
later configures work offline.  Set `-DCPM_SOURCE_CACHE=/path` to keep that
cache elsewhere.  With `-DSQC_USE_SYSTEM_FINUFFT=ON -DSQC_FINUFFT_ROOT=/path`
an installed FINUFFT can be used instead.

## Usage

```
sqcalc -i DUMP [options] OUTPUT
```

`OUTPUT` is the shell averaged $S(q)$ table; use `-` to write it to stdout.

| option | meaning |
| --- | --- |
| `-i, --input FILE` | LAMMPS dump trajectory (required) |
| `-m, --mapping LIST` | LAMMPS type id to element symbol, e.g. `1:Si,2:O` |
| `-w, --weight SCHEME` | `unit` (default), `neutron` or `xray` |
| `-t, --threads N` | OpenMP threads (default: all available) |
| `--qmin`, `--qmax`, `--nq` | $q$ range and number of shells (defaults 0, 20 1/A, 500) |
| `--grid FILE` | also write $S(q)$ on every reciprocal lattice point |
| `--grid-format NAME` | `text` (default) or `hdf5`; a `.h5`/`.hdf5` name implies hdf5 |
| `--method NAME` | `nufft` (default), `direct` ($O(N \cdot N_{\text{modes}})$ reference) or `debye` (real space pair histograms, see below) |
| `--device NAME` | `cpu` (default) or `gpu` (needs `-DSQC_ENABLE_CUDA=ON`) |
| `--gpu-id N` | CUDA device to use when `--device gpu` (default 0) |
| `--precision NAME` | `double` (default) or `single` (float32, GPU only) |
| `--norm NAME` | `mean` (default), `self` or `n` |
| `--eps VALUE` | NUFFT tolerance (default 1e-9; 1e-5 for the float32 GPU path) |
| `--partials` | append the partial structure factor columns (default: on when `-m` is given) |
| `--no-partials` | do not append the partial columns |
| `-fz, --faber-ziman` | report the partials in the Faber-Ziman normalization |
| `-q, --quiet` | suppress progress output on stderr |
| `-h, --help` | show the option summary |
| `-v, --version` | print the program version |

Examples

```sh
# unit weights, stdout, 4 threads
sqcalc -i traj.dump -t 4 - > S_q.dat

# neutron weighting of a two component glass, plus the reciprocal grid table
sqcalc -i traj.dump -m 1:Si,2:O -w neutron --qmax 25 --nq 1000 \
       --grid S_q_grid.dat S_q.dat

# X-ray weighting (q dependent form factors) and a cross check with the
# brute force implementation
sqcalc -i traj.dump -m 1:Si,2:O -w xray --method direct --qmax 6 S_q_direct.dat

# same calculation on the GPU (cufinufft / cuFFT)
sqcalc -i traj.dump -m 1:Si,2:O -w unit --device gpu --qmax 20 -t 4 S_q_gpu.dat

# float32 GPU transform: about 2.5x faster than float64 on a consumer GPU,
# accurate to a few 1e-6 in S(q)
sqcalc -i traj.dump -m 1:Si,2:O -w unit --device gpu --precision single S_q_gpu32.dat

# reciprocal grid into HDF5 (much faster to write than the text table)
sqcalc -i traj.dump -m 1:Si,2:O -w unit --grid S_q_grid.h5 S_q.dat
```

## GPU execution (cuFFT)

Configure with the CUDA backend to get `--device gpu`, which runs the type-1
transforms with cufinufft (cuFFT) while the cheap combine/binning stays on the
host:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSQC_ENABLE_CUDA=ON
cmake --build build -j
ctest --test-dir build -R sqcalc_gpu --output-on-failure
```

Notes

* The CUDA toolkit is required; the compute capability is taken from
  `nvidia-smi`, or pass `-DSQC_CUDA_ARCHITECTURES=<cc>` (for example `75`) when
  configuring on a machine without a visible GPU.
* FINUFFT's CUDA build downloads CCCL from GitHub once (even for CUDA 12, where
  the toolkit ships its own copy).  It is cached in `.cpm-cache/`, so later
  reconfigures work offline.  Point `-DCPM_CCCL_SOURCE=<dir>` at a local CCCL to
  avoid that download entirely.
* cufinufft is linked as a shared library, so the CPU `libfinufft` is built
  shared as well in a CUDA build.
* Accuracy is the same as on the CPU (`--eps` controls the tolerance).
  Measured on an RTX 2060 (10 000 atoms, 10 frames, $193^3$ grid, 4 threads):
  CPU 17.3 s, GPU float64 5.6 s.  The GPU wins for large grids (roughly
  $q_{\max} L > 300$, i.e. from a $97^3$ grid upwards here); for small grids the
  fixed plan cost and per-frame overhead make the CPU faster ($49^3$ grid:
  CPU 0.6 s, GPU 2.1 s).
* The `--method direct` reference path is CPU only, and `--grid` accumulation is
  done on the host after downloading the transform.
* Consumer GPUs run double precision at a much lower rate (the RTX 2060 at
  1/32 of fp32), so the float64 GPU win is limited by the FP64 FFT and
  spreading; float32 is substantially faster here and accurate enough for $S(q)$.
* `--precision single` selects the float32 GPU transform (cufinufftf).  It only
  affects the transform: coordinates and strengths are converted down before the
  upload and the grid is widened back to double for the combine/binning, so the
  accumulators stay in double.  Because float32 cannot reach the double
  precision default, `--eps` defaults to `1e-5` in that mode (pass `--eps`
  explicitly to override).  Measured on the RTX 2060 (10k atoms, 10 frames,
  $193^3$ grid): CPU 17.3 s, GPU float64 5.6 s, GPU float32 2.1 s; float32
  differs from float64 by ~2e-6 in $S(q)$.

## Debye method (real space pair histograms)

`--method debye` replaces the reciprocal space transform by counting atom pairs
as a function of distance, for every pair of chemical types, and averaging those
histograms over the trajectory - the recipe used by the `debyer` program:

$$
S(q) = \frac{1}{W(q)} \left[
  \sum_{ab} f_a(q) f_b(q) \sum_k n_{ab}[k] \frac{\sin(q r_k)}{q r_k} +
  \sum_a N_a f_a(q)^2 \right]
$$

with $r_k = (k - \tfrac{1}{2})\,dr$ the bin centers, the second sum the $r = 0$
self term (an atom with itself, added analytically because no histogram can hold
it) and $W(q)$ the same normalization as the grid method, so `--weight`, `--norm`
and the output format behave identically:

```sh
sqcalc -i traj.dump -m 1:Si,2:O -w neutron --method debye \
       --rmax 12 --dr 0.01 --skin 1.0 --rdf g_of_r.dat S_q.dat
```

| option | meaning |
| --- | --- |
| `--rmax VALUE` | pair cutoff [A]; default: half of the smallest periodic box side (minimum image convention, so no self images), or all pairs when the box is not periodic |
| `--dr VALUE` | radial bin width (default 0.01 A, reliable up to $q \sim \pi/(2\,dr)$; ~150 1/A at the default) |
| `--skin VALUE` | Verlet skin for reusing the pair list (default 1.0 A, `0` rebuilds every frame) |
| `--rdf FILE` | total and all partial $g(r)$ in one file (text, or HDF5 for `.h5`/`.hdf5`) |
| `--no-cutoff-correction` | disable the cut-off density correction (for comparison; applied by default when at least one direction is periodic) |

### Partial structure factors

Whenever an element mapping (`-m`) is given, the $S(q)$ table gets one extra
column per type pair, using the convention of OVITO's structure factor modifier:

```
# q S(q) S(Si-Si) S(Si-O) S(O-O)
```

with

$$
S_{ab}(q) = \frac{1}{N} \left\langle
  \sum_{j \in a} \sum_{k \in b} \frac{\sin(q r_{jk})}{q r_{jk}} \right\rangle
$$

including the $r = 0$ self term for $a = b$.  The columns satisfy the sum rule

$$
S(q) = \sum_{ab} (2 - \delta_{ab})\, S_{ab}(q),
$$

and at large $q$, $S_{aa} \to x_a$ (the concentration) while $S_{ab} \to 0$ for
$a \neq b$.  `-fz` switches to the Faber-Ziman form used in PDF work,

$$
A_{ab}(q) = \frac{S_{ab}(q) - x_a \delta_{ab}}{x_a x_b} + 1,
$$

which tends to 1 for every pair (and implies the equivalent weighted sum rule
$S(q) = \sum_{ab} (2 - \delta_{ab})\, x_a x_b A_{ab}(q)\, w_a w_b / \langle w \rangle^2$).
Partials are available for both the reciprocal and the Debye method;
`--no-partials` turns them off, and they are skipped automatically when no
mapping is given.

Verified: the sum rule holds to 4.5e-13 for the Debye method and 7e-3 for the
reciprocal method (the latter residual is the per-shell mode weighting of the
total), and the high-q limits are reproduced by both methods
($S_{aa} = 0.50 = x_a$, $S_{ab} < 0.01$, $A_{ab} = 0.99$).

The `-fz` columns are checked against the algebraic transformation of the plain
partial columns ($A_{ab} = (S_{ab} - x_a \delta_{ab})/(x_a x_b) + 1$, agreement
6e-13) and against the equivalent weighted sum rule
$S(q) = \sum_{ab} (2 - \delta_{ab})\, x_a x_b A_{ab}(q)$ (Debye 5e-13).

The reciprocal grid table (`--grid`, text or HDF5) deliberately keeps writing
the total $S(q)$ only - a per-partial decomposition there would multiply the
storage for output that is rarely used.  Partial $g(r)$ are available in the
`--rdf` file (Debye method, text and HDF5).

In the HDF5 output the partials *are* stored for the shell-averaged table
(small: `npairs x nq` values), as one dataset per pair under `/shell/S_partial/`
(`Si-Si`, `Si-O`, ...), with the label list in `/shell/pairs` and an attribute
`partial_normalization` recording whether they are `ovito` or `faber-ziman`
(the same distinction the text header makes with `S(...)` versus `A(...)`).
The grid group stays total-only, as decided.

Cross-check of the partial $g(r)$ against debyer (`-g -p`, which writes one
column per pair plus a sum): debyer's partial columns are scaled by $x_a x_b$ and
use its half-list pair counting, so dividing by $x_a x_b$ brings them onto our
definition - then $g_{\mathrm{SiSi}} = 0.988$ (debyer) versus $1.029$ (ours) and
$g_{\mathrm{OO}} = 0.979$ versus $1.027$, i.e. agreement within the statistics of
the 250-atom, 4-frame test trajectory.  Our partials use the standard
normalization $g_{ab} \to 1$.

`--grid` is not available with this method (it evaluates $S(q)$ directly) and it
runs on the CPU.

The RDF file holds one row per r bin, `# r g(r) g(a-a) g(a-b) ...` in the order
given by `-m`; the partials are normalized so that $g_{ab} \to 1$ at large $r$
(HDF5: `/rdf/r`, the total in `/rdf/g`, one dataset per pair under
`/rdf/g_partial/<label>` and the labels in `/rdf/pairs`).  The weighted total
uses the $q \to 0$ amplitudes, the usual PDF convention.

When to use it (measured here, 10k atoms, 3-10 frames, $q_{\max} = 20$ 1/A):

| system | Debye | reciprocal grid (CPU, 4 threads) |
| --- | --- | --- |
| 30 A box ($r_{\max} = L/2$, essentially all pairs) | 44.5 s | 8.1 s |
| 60 A box ($r_{\max} = 15$ A, ~650 neighbours per atom) | 5.7 s | 83.9 s |

The pair count grows like $N \cdot n_{\text{neigh}}(r_{\max})$ while the grid
method grows like $(q_{\max} L)^3$, so Debye wins for large boxes and high $q$,
and the grid method wins for small, dense, periodically replicated cells.

Both methods accept non-periodic boxes (`ITEM: BOX BOUNDS ff ff ff`).  The
reciprocal method wraps each atom into the dump box, which leaves $\rho(q)$
exactly unchanged for that box's reciprocal lattice vectors, so a non-periodic
configuration is transformed correctly - but $q$ is only sampled on that
lattice, so a cluster needs enough vacuum padding that $2\pi/L$ resolves the
structure, and the $(q_{\max} L)^3$ grid cost grows with the padding.  Debye
needs no padding: it reads the `pp`/`ff` flags per axis, evaluates $S(q)$ on any
requested $q$ grid and defaults to all pairs (`--rmax` = diagonal of the atom
cloud) when no direction is periodic, which usually makes it the cheaper choice
for clusters,
nanoparticles and surfaces.

The Verlet skin only pays off when consecutive frames are *correlated*: on a
rattled trajectory the list survives many frames (the test suite checks that 20
rattled frames cost a single rebuild), while on a series of independent random
configurations every frame rebuilds the list and the larger candidate list makes
`--skin 0` about 4x faster.  Use `--skin 0` for uncorrelated frames.

### Cross-check against the debyer program

`test/dump_to_cfg.py` converts dump frames to AtomEye `.cfg` files so the
external `debyer` binary can be run on exactly the same configurations (one file
per frame, then averaging its $S(q)$ curves).  debyer is *not* a build or test
dependency; the script is a manual validation helper.

On an ideal gas (250 atoms, 20 A box, $r_{\max} = 6$ A, 4 frames, `-c sf` versus
`--weight unit --norm n`) the raw pair histograms agree in structure but not in
detail - the uncorrected sqcalc curve still carries the finite cut-off artifacts
(the low-q rise and the residual few-percent bias at high q):

| q [1/A] | 0.6 | 1.6 | 2.6 | 3.6 | 4.6 | 5.6 |
| --- | --- | --- | --- | --- | --- | --- |
| debyer | 0.91 | 1.02 | 1.02 | 0.98 | 1.01 | 1.00 |
| sqcalc (uncorrected) | 5.98 | 1.91 | 1.37 | 1.15 | 1.10 | 1.05 |

debyer removes those artifacts with the cut-off density correction
(`add_cutoff_correction`), which it applies automatically whenever a cut-off is
given; sqcalc implements the same correction and enables it by default
(`--no-cutoff-correction` disables it).  With it, sqcalc reproduces debyer's
curve for this configuration **exactly** (max |difference| = 0.0000 over the 27
q points), and the ideal gas RMS deviation from $S(q) = 1$ drops from 0.267
(uncorrected) to 0.0146, identical to debyer's own value.  The pair counting
itself is validated independently: `test/debye_ref.py` reproduces sqcalc's
$S(q)$ and $g(r)$ exactly, and both programs reach the correct high-q plateau
(1.00).

The same configuration against our reciprocal (NUFFT) method, all three
evaluated on the same $q$ grid.  An ideal gas has $S(q) = 1$ everywhere, so any
deviation is an artifact:

| q [1/A] | 0.6 | 1.2 | 1.8 | 2.4 | 3.0 | 3.6 | 4.2 | 4.8 | 5.4 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| debyer | 0.91 | 1.03 | 1.01 | 1.01 | 1.02 | 0.98 | 0.99 | 1.01 | 1.00 |
| sqcalc Debye | 5.98 | 0.21 | 1.09 | 1.15 | 0.84 | 1.15 | 0.86 | 1.10 | 0.96 |
| sqcalc Debye (corrected) | 0.91 | 1.03 | 1.01 | 1.01 | 1.02 | 0.98 | 0.99 | 1.01 | 1.00 |
| sqcalc NUFFT | 1.28 | 1.04 | 0.99 | 1.01 | 0.99 | 1.02 | 1.00 | 1.01 | 1.00 |

RMS deviation from 1 for $q > 1.5$ 1/A: debyer 0.015, our corrected Debye 0.015,
our NUFFT 0.024, our uncorrected Debye 0.267.  The reciprocal method and the
corrected Debye method agree with debyer and with the exact ideal-gas value; the
small NUFFT scatter is the finite shell resolution plus 4-frame statistics, not
physics.

## Output

The first table is the isotropic average over q shells:

```
# q S(q)
      0.100000    1.234567890123E+00
      ...
```

`--grid FILE` writes every reciprocal lattice vector of the requested range,
averaged over the trajectory:

```
# qx qy qz S(q)
```

Both tables exclude the trivial $q = 0$ mode.

### HDF5 grid output

When `--grid` names a `.h5`/`.hdf5` file (or `--grid-format hdf5` is given) the
results are written as HDF5 instead of text.  This needs HDF5 to have been found
at configure time (`-DSQC_ENABLE_HDF5=OFF` disables it).  Every dataset is a
plain one dimensional array, so no reader has to care about dimension ordering:

```
/shell/q        nq      doubles   shell centers [1/A]
/shell/S        nq      doubles   S(q) averaged over the shell
/shell/count    nq      int64     reciprocal lattice points per shell
/grid/qx,qy,qz  nmodes  doubles   reciprocal lattice vector [1/A]
/grid/h,k,l     nmodes  int32     Miller indices
/grid/S         nmodes  doubles   S(q) at that reciprocal lattice point
```

When partial structure factors are written they are added as one dataset per
pair under `/shell/S_partial/<pair>`, with the label list in `/shell/pairs` (see
[Partial structure factors](#partial-structure-factors)).

File attributes carry the run metadata: `nframes`, `natoms`, `nmodes`, `nq`,
`qmin`, `qmax`, `eps`, `cell` (3x3), `weight`, `norm`, `device`, `precision`
and `mapping`.

```python
import h5py, numpy as np
with h5py.File("S_q_grid.h5") as f:
    q = np.column_stack([f["grid/qx"][:], f["grid/qy"][:], f["grid/qz"][:]])
    s = f["grid/S"][:]
    cell = f.attrs["cell"]          # 3x3 matrix
```

Writing 3.65M grid points costs ~0.05 s as HDF5 versus ~5 s as text (the text
file is 252 MB, the HDF5 file 160 MB).

## Conventions

For every frame the reciprocal methods (NUFFT and direct) evaluate the scattering
amplitude on the reciprocal lattice of the dump box,
$q = h b_1 + k b_2 + l b_3$ with $a_i \cdot b_j = 2\pi \delta_{ij}$ (orthogonal
and restricted triclinic boxes are supported):

$$
\rho(q) = \sum_j w_j \exp(i\, q \cdot r_j)
$$

with the per-atom weight $w_j$

* `unit`   : $w = 1$ for every atom,
* `neutron`: $w = b(\text{element})$, the bound coherent neutron scattering length,
* `xray`   : $w = f(\text{element}, q) = \sum_i a_i \exp(-b_i (q/4\pi)^2) + c$ (IT92).

$S(q)$ is the trajectory average

$$
S(q) = \frac{\langle |\rho(q)|^2 \rangle}{W(q)}
$$

where the normalization is

| `--norm` | $W(q)$ |
| --- | --- |
| `mean` (default) | $N \langle w \rangle^2$, the Faber-Ziman total $S(q)$ |
| `self` | $\sum_j w_j^2$, so that $S(q) \to 1$ at large $q$ |
| `n` | $N$ (the convention used by the debyer program) |

The weights are looked up from a table of 104 elements covering the periodic
table.  Without `-m` no element is known, so every atom has weight 1.0; a
weighted scheme therefore requires a type id mapping.

Element data (masses, IT92 X-ray coefficients, NN92 neutron scattering lengths)
was converted from the debyer program (`debyer/debyer/atomtables.c`, GPL-2 or
later for the code, data taken from the International Tables for Crystallography
Vol. C (1992) table 6.1.1.4 and from Neutron News 3 (1992) 29-37).  Only the
tabulated numbers are reused here, via `tools/gen_element_data.py`; debyer's own
notice covers the code and not the data:

> Copyright 2009 Marcin Wojdyr (only code, not the tabular data)

More about the debyer program: <https://github.com/wojdyr/debyer>.

## Assumptions and limits

* The simulation box must not change along the trajectory (the reciprocal grid
  is built once from the first frame); the program stops with an error if it
  does.
* Non-periodic boxes (`ITEM: BOX BOUNDS ff ff ff`, for example a cluster in a
  bounding box) are accepted by every method.  The reciprocal methods still
  sample $q$ on the reciprocal lattice of that box, so the shell resolution and
  the grid cost stay tied to its size; Debye reads the `pp`/`ff` flags per
  direction, counts the finite pair list, and applies its isotropic cut-off
  density correction whenever at least one direction is periodic (it is skipped
  entirely for a fully non-periodic box).
* Atom coordinates are taken as dumped.  A wrapped coordinate set (`x y z`,
  the LAMMPS default) is what the reciprocal grid method expects; `xu yu zu`
  columns are accepted and wrapped internally.
* Atoms with a type that has no entry in `-m`, or an element without tabulated
  X-ray/neutron data, are rejected with a clear message.
* The grid cost grows like $(q_{\max} L)^3$; the program refuses grids larger
  than 4e8 points (about 6 GB) and asks for a smaller `--qmax`.
* `--grid` is written single threaded; the shell table is not affected.

## Tests

`ctest` runs

* `finufft_opts`: verifies at run time that the Fortran mirror of
  `finufft_opts` matches the linked FINUFFT library,
* `neighbour_list`: the Debye cell list (pair enumeration, periodic images and
  skin reuse) against brute force pair counting,
* `sqcalc_physics`: NUFFT versus direct summation versus an independent numpy
  implementation (`test/ref_sq.py`, `test/debye_ref.py`) for unit, neutron and
  X-ray weights, an ideal gas ($S \to 1$), a simple cubic lattice (Bragg peaks,
  elsewhere zero), the reciprocal grid output, the Debye method with the cut-off
  correction, $g(r)$, the partial structure factors, a triclinic box with
  `xu yu zu` columns, the HDF5 output (when HDF5 was found) and command line
  error handling,
* `sqcalc_gpu` (CUDA builds only): compares the GPU backend against the CPU for
  unit, neutron and X-ray weights, the float32 transform and the reciprocal grid
  output; it skips itself when no CUDA device is visible.

## License

Copyright (C) 2026 Rui Su, Hangzhou Dianzi University

`sqcalc` is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the GNU General Public License (the `LICENSE` file, or
<https://www.gnu.org/licenses/>) for more details.

### Third party code and dependencies

Version 3 is what covers the whole program, because the default build links
components under terms that do not otherwise combine:

* **FINUFFT** is vendored with `git subtree` in `external/finufft` (release
  v2.5.1, upstream commit `679d9ae5`).  It is Apache-2.0, Copyright (C)
  2017-2026 The Simons Foundation, Inc.; its `LICENSE` and `NOTICE` files are
  kept as upstream ships them.  No file under `external/finufft` is modified -
  the top level `CMakeLists.txt` only sets the options that build it.
  Apache-2.0 can be combined with GPLv3, but not with GPLv2, which is why the
  license above is version 3.
* **FFTW3** (GPLv2 or later) is what FINUFFT transforms with, unless FINUFFT is
  built against DUCC0 instead, so a normally built `sqcalc` binary links a GPL
  library.
* **HDF5** (BSD-3-Clause) backs `--grid-format hdf5`, and the OpenMP runtime
  does the threading.
* **CUDA / cuFFT**, used only by the optional `-DSQC_ENABLE_CUDA=ON` build (and
  `--device gpu`), is under NVIDIA's terms.

The element tables in `src/sqc_element_data.f90` are generated by
`tools/gen_element_data.py`; their provenance is described under
[Conventions](#conventions).
