# sqcalc

`sqcalc` computes the total structure factor $S(q)$ of a LAMMPS trajectory.  It
reads the default LAMMPS `atoms` dump style, averages over all snapshots and
writes a `# q S(q)` table.  Two independent evaluations are available: the
reciprocal space transform with a non-uniform FFT (FINUFFT) on the CPU, or on a
CUDA GPU with cuFFT, which averages $|\rho(q)|^2$, and the real space Debye
pair-histogram method.  The table also carries the partial structure factors
$S_{ab}(q)$, one column per pair of LAMMPS types.

With `--dyn` the time axis is kept instead of being averaged away, and the
run writes the dynamic structure factor $S(q,\omega)$ along a line in
reciprocal space, plus optionally the intermediate scattering function
$F(q,t)$.  That method correlates the density amplitudes in time rather than
averaging their squares.

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
first configure of the vendored FINUFFT obtains its build helpers (`CPM.cmake`,
`findFFTW`, `xsimd` and, unless a system FFTW is found, FFTW3) over the network;
they are cached in `.cpm-cache/` next to the sources, so later configures work
offline.  Set `-DCPM_SOURCE_CACHE=/path` to keep that cache elsewhere.  With
`-DSQC_USE_SYSTEM_FINUFFT=ON -DSQC_FINUFFT_ROOT=/path` an installed FINUFFT can
be used instead.

The optional CUDA backend, which adds `--device gpu`, is configured with

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSQC_ENABLE_CUDA=ON
cmake --build build -j
```

The CUDA toolkit is required; the compute capability is taken from
`nvidia-smi`, or pass `-DSQC_CUDA_ARCHITECTURES=<cc>` (for example `75`) when
configuring on a machine without a visible GPU.  Measurements on an RTX 2060
show the GPU path winning on large grids, while `--precision single` (float32,
faster but accurate to only a few 1e-6 in $S(q)$) is the better choice on
consumer cards whose FP64 rate is low.

## Usage

```
sqcalc -i DUMP [options] OUTPUT
```

`OUTPUT` is the shell averaged $S(q)$ table; use `-` to write it to stdout.  The
shell table is text, `# q S(q)`; `--grid FILE` additionally writes every
reciprocal lattice vector, in text (`# qx qy qz S(q)`) or as a one dimensional
array per quantity in HDF5 when the file name ends in `.h5`/`.hdf5`.  Both
tables exclude the trivial $q = 0$ mode.  When partials are written they appear
as one extra column per type pair, labelled with the element symbols of `-m`
(e.g. `# q S(q) S(Si-Si) S(Si-O) S(O-O)`) or, without a mapping, with the
LAMMPS type ids (`# q S(q) S(1-1) S(1-2) S(2-2)`).

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
| `--partials` | append the partial structure factor columns (default: on) |
| `--no-partials` | do not append the partial columns |
| `-fz, --faber-ziman` | report the partials in the Faber-Ziman normalization |
| `-q, --quiet` | suppress progress output on stderr |
| `-h, --help` | show the option summary |
| `-v, --version` | print the program version |

With `--method debye` these options select the pair histogram:

| Debye option | meaning |
| --- | --- |
| `--rmax VALUE` | pair cutoff [A]; default: half of the smallest periodic box side (minimum image convention, so no self images), or all pairs when the box is not periodic |
| `--dr VALUE` | radial bin width (default 0.01 A, reliable up to $q \sim \pi/(2\,dr)$; ~150 1/A at the default) |
| `--skin VALUE` | Verlet skin for reusing the pair list (default 1.0 A, `0` rebuilds every frame) |
| `--rdf FILE` | total and all partial $g(r)$ in one file (text, or HDF5 for `.h5`/`.hdf5`) |
| `--no-cutoff-correction` | disable the cut-off density correction (for comparison; applied by default when at least one direction is periodic) |
| `--dyn NINT,S0,S1,DX,DY,DZ` | dynamic structure factor $S(q,\omega)$ along a $q$ line (see below) |
| `--dt VALUE` | time step of the trajectory (the LAMMPS `timestep`), required by `--dyn` |
| `--maxframes L` | correlation window in frames (largest lag kept), required by `--dyn` |
| `--lag N` | frames between consecutive time origins (default 1) |
| `--sqw FILE` | the $S(q,\omega)$ spectra (text, or HDF5 for `.h5`/`.hdf5`) |
| `--fsq FILE` | the intermediate scattering function $F(q,t)$ |
| `--sqw-format NAME` | `text` (default) or `hdf5` for the two files above |

The HDF5 output needs HDF5 to have been found at configure time
(`-DSQC_ENABLE_HDF5=OFF` disables it); its attributes record the run metadata
(`nframes`, `natoms`, `nmodes`, `nq`, `qmin`, `qmax`, `eps`, `cell`, `weight`,
`norm`, `device`, `precision`, `mapping`), and the partials are stored as one
dataset per pair under `/shell/S_partial/`.

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

# the same calculation on the GPU, and the real space Debye method with g(r)
sqcalc -i traj.dump -m 1:Si,2:O -w unit --device gpu --qmax 20 -t 4 S_q_gpu.dat
sqcalc -i traj.dump -m 1:Si,2:O -w neutron --method debye \
       --rmax 12 --dr 0.01 --skin 1.0 --rdf g_of_r.dat S_q.dat

# dynamic structure factor along (1,1,0): 101 q points from 0.5 to 20 1/A,
# frames every 10 steps of dt = 0.005, a 400 frame correlation window
sqcalc -i traj.dump -m 1:Si,2:O -w neutron \
       --dyn 100,0.5,20,1,1,0 --dt 0.005 --maxframes 400 \
       --sqw S_qw.dat S_q.dat
```

## Theoretical background

For every frame the reciprocal methods (NUFFT and direct) evaluate the
scattering amplitude on the reciprocal lattice of the dump box,
$q = h b_1 + k b_2 + l b_3$ with $a_i \cdot b_j = 2\pi \delta_{ij}$
(orthogonal and restricted triclinic boxes are supported):

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

The Debye method replaces the reciprocal space transform by counting atom pairs
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
and the output format behave identically.  A finite pair cutoff biases the low
$q$ region; sqcalc applies the same cut-off density correction as `debyer`
whenever at least one direction is periodic (`--no-cutoff-correction` disables
it).  `--grid` is not available with this method, and it runs on the CPU.

### Dynamic structure factor

`--dyn` follows the recipe used by the dynasor and MDANSE packages: pick a
line in reciprocal space, evaluate the density amplitudes there, correlate
them in time and Fourier transform the correlation.  With
`q_i = s_i \hat{u}`, `s_i = S_0 + i (S_1-S_0)/N_{\rm int}` (the line always
passes through $\Gamma$, the origin of reciprocal space),

$$
F(q,t) = \frac{1}{W(q)} \sum_{ab} w_a(q) w_b(q)
         \frac{\langle \rho_a(q,t'+\tau)\,\rho_b^*(q,t')\rangle}{N_{\rm origins}(\tau)},
\qquad
S(q,\omega) = \frac{1}{2\pi}\int_{-\infty}^{\infty} e^{i\omega t} F(q,t)\,dt
$$

with the same $\rho$, $w$ and $W(q)$ as the static calculation, so the
normalization conventions (`--weight`, `--norm`) carry over unchanged and the
zeroth moment is exact: $\int S(q,\omega)\,d\omega = F(q,0) = S(q)$ of the
`OUTPUT` table.

The average over time origins accumulates a sum *and* a count for every lag,
which is what makes the estimator exact for a trajectory that is not a whole
number of windows long, and the ring buffer means the memory does not grow
with the trajectory.  `--maxframes L` sets the largest lag kept (and hence the
resolution $\Delta\omega = \pi/(L\Delta t)$), while `--lag N` thins the time
origins for cheaper or more nearly independent averages.  The frame interval
$\Delta t$ is `--dt` times the step increment read from the dump, which must
be constant; the highest frequency the spectrum can represent is
$\pi/\Delta t$, so the dump cadence - not `--maxframes` - decides what is
resolvable.

The $q$ line may be off the reciprocal lattice of the box.  That is the
experimental situation (the $q$ of an inelastic scattering experiment is set by
the scattering angle, not by the simulation box) and it is why the amplitudes
are summed directly instead of being transformed on a grid: for the ~100 q
points of a line scan the direct sum costs $O(N\,n_q)$ per frame, far less
than transforming the $(q_{\max}L)^3$ grid.  Off-lattice $q$ requires
unwrapped coordinates, which `--dyn` reconstructs while reading if the dump
carries only wrapped `x y z` (preferring `xu yu zu` when available).

### Partial structure factors

The $S(q)$ table gets one extra column per type pair, using the convention of
OVITO's structure factor modifier:

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

which tends to 1 for every pair.  Partials are available for both the
reciprocal and the Debye method; `--no-partials` turns them off.  They only
need the LAMMPS type ids, so a model system without elements gets them too,
labelled `S(1-1)`, `S(1-2)`, ...; `-m` replaces those labels with element
symbols and is required for the weighted schemes.  Report them with `-fz`:
every Faber-Ziman pair tends to 1, so the curves can be read and compared
across concentrations, which the OVITO pairs (tending to $x_a$) do not allow.

A few limits are worth keeping in mind: the simulation box must not change
along the trajectory (the reciprocal grid is built once from the first frame),
both methods accept a non-periodic box (the reciprocal one still samples $q$ on
the lattice of that box, so a cluster needs enough vacuum padding, whereas
Debye needs none).  That lattice also fixes which $q$ are measurable: set
`--qmin` to the smallest reciprocal vector and keep the shell width
`(qmax - qmin)/nq` comparable to the lattice spacing, or the table comes back
with zero-filled shells the box can never fill; the skill's
`scripts/choose_q.py` prints a safe set of options for a given dump.  Atoms of a
type missing from `-m` under `-w neutron`/`xray`, or with an element that has
no tabulated data, are rejected with a clear message.

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
* The element data has a different origin from the code: masses, IT92 X-ray
  coefficients and NN92 neutron scattering lengths were converted from the
  debyer program (`debyer/debyer/atomtables.c`) by `tools/gen_element_data.py`.
  Only the tabulated numbers are reused, taken from the International Tables
  for Crystallography Vol. C (1992) table 6.1.1.4 and from Neutron News 3
  (1992) 29-37; debyer's GPL notice covers its code and not the data
  ("Copyright 2009 Marcin Wojdyr (only code, not the tabular data)").  More
  about debyer: <https://github.com/wojdyr/debyer>.
