# sqcalc

`sqcalc` computes the total structure factor S(q) of a LAMMPS trajectory.  It
reads the default LAMMPS `atoms` dump style, evaluates the scattering amplitude
on the reciprocal lattice of the dump box with a non-uniform FFT (FINUFFT),
averages |rho(q)|^2 over all snapshots and writes a `# q S(q)` table.

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
ctest --test-dir build          # optional, runs the physics test suite
```

If FFTW lives in a non standard prefix pass `-DFFTW_ROOT=/path/to/fftw`.  The
first configure of the vendored FINUFFT fetches `CPM.cmake`, `findFFTW` and
`xsimd` (network access needed once, cached in the build directory afterwards);
set `-DCPM_SOURCE_CACHE=/path` to keep that cache outside the build tree.  With
`-DSQC_USE_SYSTEM_FINUFFT=ON -DSQC_FINUFFT_ROOT=/path` an installed FINUFFT can
be used instead.

## Usage

```
sqcalc -i DUMP [options] OUTPUT
```

`OUTPUT` is the shell averaged S(q) table; use `-` to write it to stdout.

| option | meaning |
| --- | --- |
| `-i, --input FILE` | LAMMPS dump trajectory (required) |
| `-m, --mapping LIST` | LAMMPS type id to element symbol, e.g. `1:Si,2:O` |
| `-w, --weight SCHEME` | `unit` (default), `neutron` or `xray` |
| `-t, --threads N` | OpenMP threads (default: all available) |
| `--qmin`, `--qmax`, `--nq` | q range and number of shells (defaults 0, 20 1/A, 500) |
| `--grid FILE` | also write S(q) on every reciprocal lattice point |
| `--grid-format NAME` | `text` (default) or `hdf5`; a `.h5`/`.hdf5` name implies hdf5 |
| `--method NAME` | `nufft` (default) or `direct` (O(N * Nmodes) reference) |
| `--device NAME` | `cpu` (default) or `gpu` (needs `-DSQC_ENABLE_CUDA=ON`) |
| `--gpu-id N` | CUDA device to use when `--device gpu` (default 0) |
| `--precision NAME` | `double` (default) or `single` (float32, GPU only) |
| `--norm NAME` | `mean` (default), `self` or `n` |
| `--eps VALUE` | FINUFFT tolerance (default 1e-9) |
| `-q, --quiet` | suppress progress output on stderr |

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
  Measured on an RTX 2060 (10 000 atoms, 10 frames, 193^3 grid, 4 threads):
  CPU 17.3 s, GPU 5.5 s.  The GPU wins for large grids (roughly qmax * L > 300,
  i.e. from a 97^3 grid upwards here); for small grids the fixed plan cost and
  per-frame overhead make the CPU faster (49^3 grid: CPU 0.6 s, GPU 2.1 s).
* The `--method direct` reference path is CPU only, and `--grid` accumulation is
  done on the host after downloading the transform.
* Consumer GPUs run double precision at a much lower rate (the RTX 2060 at
  1/32 of fp32), so the GPU win is limited by the FP64 FFT and spreading; a
  single precision GPU transform would be substantially faster and is accurate
  enough for S(q) - a possible follow-up.
* `--precision single` selects the float32 GPU transform (cufinufftf).  It only
  affects the transform: coordinates and strengths are converted down before the
  upload and the grid is widened back to double for the combine/binning, so the
  accumulators stay in double.  Because float32 cannot reach the double
  precision default, `--eps` defaults to `1e-5` in that mode (pass `--eps`
  explicitly to override).  Measured on the RTX 2060 (10k atoms, 10 frames,
  193^3 grid): CPU 17.3 s, GPU float64 5.6 s, GPU float32 2.1 s; float32 differs
  from float64 by ~2e-6 in S(q).

## Debye method (real space pair histograms)

`--method debye` replaces the reciprocal space transform by counting atom pairs
as a function of distance, for every pair of chemical types, and averaging those
histograms over the trajectory - the recipe used by the `debyer` program:

```
S(q) = [ sum_ab f_a(q) f_b(q) sum_k n_ab[k] sin(q r_k)/(q r_k)
         + sum_a N_a f_a(q)^2 ] / W(q)
```

with `r_k = (k-1/2) dr` the bin centers, the second sum the r = 0 self term (an
atom with itself, added analytically because no histogram can hold it) and `W(q)`
the same normalization as the grid method, so `--weight`, `--norm` and the output
format behave identically:

```sh
sqcalc -i traj.dump -m 1:Si,2:O -w neutron --method debye \
       --rmax 12 --dr 0.01 --skin 1.0 --rdf g_of_r.dat S_q.dat
```

| option | meaning |
| --- | --- |
| `--rmax VALUE` | pair cutoff [A]; default: half of the smallest periodic box side (minimum image convention, so no self images), or all pairs when the box is not periodic |
| `--dr VALUE` | radial bin width (default 0.01 A, reliable up to q ~ pi/(2 dr) ~ 150 1/A) |
| `--skin VALUE` | Verlet skin for reusing the pair list (default 1.0 A, `0` rebuilds every frame) |
| `--rdf FILE` | total and all partial g(r) in one file (text, or HDF5 for `.h5`) |

`--grid` is not available with this method (it evaluates S(q) directly) and it
runs on the CPU.

The RDF file holds one row per r bin, `# r g(r) g(a-a) g(a-b) ...` in the order
given by `-m`; the partials are normalized so that `g_ab -> 1` at large r (HDF5:
`/rdf/r`, `/rdf/g`, one dataset per pair under `/rdf/g/<label>` and the labels in
`/rdf/pairs`).  The weighted total uses the q -> 0 amplitudes, the usual PDF
convention.

When to use it (measured here, 10k atoms, 3-10 frames, qmax 20 1/A):

| system | Debye | reciprocal grid (CPU, 4 threads) |
| --- | --- | --- |
| 30 A box (rmax = L/2, essentially all pairs) | 44.5 s | 8.1 s |
| 60 A box (rmax = 15 A, ~650 neighbours per atom) | 5.7 s | 83.9 s |

The pair count grows like `N * neighbours(rmax)` while the grid method grows like
`(qmax*L)^3`, so Debye wins for large boxes and high q, and the grid method wins
for small, dense, periodically replicated cells.  Debye is also the only option
for non-periodic systems (clusters, nanoparticles, surfaces).

The Verlet skin only pays off when consecutive frames are *correlated*: on a
rattled trajectory the list survives many frames (the test suite checks that 20
rattled frames cost a single rebuild), while on a series of independent random
configurations every frame rebuilds the list and the larger candidate list makes
`--skin 0` about 4x faster.  Use `--skin 0` for uncorrelated frames.

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

Both tables exclude the trivial q = 0 mode.

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

For every frame the scattering amplitude is evaluated on the reciprocal lattice
of the dump box, q = h b1 + k b2 + l b3 with a_i . b_j = 2 pi delta_ij
(orthogonal and restricted triclinic boxes are supported):

```
rho(q) = sum_j w_j exp(i q . r_j)
```

with the per-atom weight `w_j`

* `unit`   : w = 1 for every atom,
* `neutron`: w = b(element), the bound coherent neutron scattering length,
* `xray`   : w = f(element, q) = sum_i a_i exp(-b_i (q/4pi)^2) + c (IT92).

S(q) is the trajectory average

```
S(q) = < |rho(q)|^2 > / W(q)
```

where the normalization is

| `--norm` | W(q) |
| --- | --- |
| `mean` (default) | N <w>^2, the Faber-Ziman total S(q) |
| `self` | sum_j w_j^2, so that S(q) -> 1 at large q |
| `n` | N (the convention used by the debyer program) |

The weights are looked up from a table of 104 elements covering the periodic
table.  Without `-m` no element is known, so every atom has weight 1.0; a
weighted scheme therefore requires a type id mapping.

Element data (masses, IT92 X-ray coefficients, NN92 neutron scattering lengths)
was converted from the debyer program (`debyer/debyer/atomtables.c`, GPL-2 for
the code, data taken from the International Tables for Crystallography Vol. C
(1992) table 6.1.1.4 and from Neutron News 3 (1992) 29-37).  Only the tabulated
numbers are reused here, via `tools/gen_element_data.py`.

More about the debyer program: <https://github.com/wojdyr/debyer>.

## Assumptions and limits

* The simulation box must not change along the trajectory (the reciprocal grid
  is built once from the first frame); the program stops with an error if it
  does.
* Atom coordinates are taken as dumped.  A wrapped coordinate set (`x y z`,
  the LAMMPS default) is what the reciprocal grid method expects; `xu yu zu`
  columns are accepted and wrapped internally.
* Atoms with a type that has no entry in `-m`, or an element without tabulated
  X-ray/neutron data, are rejected with a clear message.
* The grid cost grows like (qmax * L)^3; the program refuses grids larger than
  4e8 points (about 6 GB) and asks for a smaller `--qmax`.
* `--grid` is written single threaded; the shell table is not affected.

## Tests

`ctest` runs

* `finufft_opts`: verifies at run time that the Fortran mirror of
  `finufft_opts` matches the linked FINUFFT library,
* `sqcalc_physics`: NUFFT versus direct summation versus an independent numpy
  implementation (`test/ref_sq.py`) for unit, neutron and X-ray weights,
  an ideal gas (`S -> 1`), a simple cubic lattice (Bragg peaks, elsewhere zero),
  the reciprocal grid output and a triclinic box with `xu yu zu` columns,
  plus command line error handling.
* `sqcalc_gpu` (CUDA builds only): compares the GPU backend against the CPU for
  unit, neutron and X-ray weights and for the reciprocal grid output; it skips
  itself when no CUDA device is visible.
