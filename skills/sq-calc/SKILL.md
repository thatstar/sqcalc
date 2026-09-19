---
name: sq-calc
description: Build, run, and troubleshoot sqcalc, the Fortran/CMake structure factor program that computes S(q) and g(r) from LAMMPS dump trajectories. Use for its CMake configuration, command line options, output formats, and tests.
metadata:
  short-description: Build and use the sqcalc structure factor code
---

# sqcalc

`sqcalc` reads a LAMMPS `atoms` dump trajectory, averages over the frames, and
writes the total structure factor $S(q)$ (plus partials when an element mapping
is given) as a `# q S(q)` table.  Two independent routes are available: a
reciprocal space non-uniform FFT (FINUFFT on the CPU, cufinufft/cuFFT on a GPU)
and a real space Debye pair histogram, which can also write $g(r)$.

## Locating the program

Running sqcalc needs nothing but the `sqcalc` executable: a single binary with
no data files, plugins or environment variables.  Call it as `sqcalc` and let
the shell find it; when in doubt, `command -v sqcalc` proves it is on `PATH`.
A build tree is not required at all, so do not go looking for a source checkout
unless the binary is missing, plainly out of date, or the task is to change the
code.

This skill is versioned inside the sqcalc repository at `skills/sq-calc/`, so if
a build does turn out to be necessary, the checkout root is the parent of this
skill's own directory (upstream <https://github.com/thatstar/sqcalc>) - resolve
it from there instead of assuming a remembered location.

## Building

Needs a Fortran 2008 compiler (gfortran >= 10), CMake >= 3.20, OpenMP and an
FFTW3 installation with headers.  FINUFFT is vendored in `external/finufft`, so
only FFTW3 has to exist locally.  Run the commands from the checkout root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The freshly built binary is `build/sqcalc`; `cmake --install build --prefix
PREFIX` copies it to `PREFIX/bin/sqcalc` if it should end up on `PATH`.

Relevant CMake options:

| option | effect |
| --- | --- |
| `-DFFTW_ROOT=PREFIX` | FFTW3 prefix; otherwise CMake probes `$FFTW_ROOT`, `$CONDA_PREFIX`, `~/.pixi/envs/*`, `/usr/local` and `/usr` |
| `-DSQC_ENABLE_CUDA=ON` | adds the cufinufft/cuFFT backend (`--device gpu`); also pass `-DSQC_CUDA_ARCHITECTURES=<cc>` (e.g. `75`) when no GPU is visible at configure time |
| `-DSQC_ENABLE_HDF5=OFF` | drops the HDF5 grid/`g(r)` output (on by default when HDF5 is found) |
| `-DSQC_USE_SYSTEM_FINUFFT=ON -DSQC_FINUFFT_ROOT=PREFIX` | link an installed FINUFFT instead of the vendored subtree |
| `-DCPM_SOURCE_CACHE=DIR` | where FINUFFT's build helpers (`CPM.cmake`, `findFFTW`, `xsimd`, CCCL, possibly FFTW3) are cached; defaults to `.cpm-cache/` in the source tree |
| `-DCMAKE_BUILD_TYPE=Debug` | adds `-fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow` |

The first configure of the vendored FINUFFT downloads those helpers, and the
CUDA build additionally downloads CCCL; both are cached, so later configures
work offline.  A CUDA build links cufinufft shared, which also makes the CPU
`libfinufft` shared.  Network is needed only for that first configure, so on a
sandboxed machine expect to request approval for it, or point
`-DCPM_SOURCE_CACHE` / `-DFFTW_ROOT` at an existing cache first.

## Running

```sh
sqcalc -i DUMP [options] OUTPUT
```

`OUTPUT` is the shell averaged $S(q)$ table, `-` writes it to stdout.  `-i` is
the only required option; `sqcalc -h` prints the built-in summary.

| option | meaning |
| --- | --- |
| `-i, --input FILE` | LAMMPS dump trajectory (required) |
| `-m, --mapping LIST` | LAMMPS type id to element symbol, e.g. `1:Si,2:O`; required for the `neutron`/`xray` weights and for partials |
| `-w, --weight SCHEME` | `unit` (default), `neutron` or `xray` (q dependent IT92 form factors) |
| `-t, --threads N` | OpenMP threads (default: all available) |
| `--qmin`, `--qmax`, `--nq` | $q$ range in 1/A and number of shells (defaults 0, 20, 500) |
| `--grid FILE` | also write $S(q)$ on every reciprocal lattice point |
| `--grid-format NAME` | `text` (default) or `hdf5`; a `.h5`/`.hdf5` name implies hdf5 |
| `--method NAME` | `nufft` (default), `direct` (slow $O(N N_{\text{modes}})$ reference) or `debye` |
| `--device NAME`, `--gpu-id N` | `cpu` (default) or `gpu` on CUDA builds, and which device |
| `--precision NAME` | `double` (default) or `single` (float32, GPU only; `--eps` then defaults to `1e-5`) |
| `--norm NAME` | `mean` (default, Faber-Ziman), `self` ($S \to 1$ at large $q$) or `n` (debyer convention) |
| `--eps VALUE` | NUFFT tolerance (default `1e-9`) |
| `--partials`, `--no-partials`, `-fz` | partial $S_{ab}(q)$ columns (on by default with `-m`), and the Faber-Ziman $A_{ab}(q)$ form |
| `--rmax`, `--dr`, `--skin`, `--rdf`, `--no-cutoff-correction` | Debye method: pair cutoff, radial bin width, Verlet skin, $g(r)$ output file, and switching off the cut-off density correction |
| `-q, --quiet`, `-h, --help`, `-v, --version` | progress output, help, version |

Typical invocations:

```sh
# neutron weighted two component glass, shell table plus reciprocal grid
sqcalc -i traj.dump -m 1:Si,2:O -w neutron --qmax 25 --nq 1000 \
       --grid S_q_grid.dat S_q.dat

# real space Debye method with the partial g(r)
sqcalc -i traj.dump -m 1:Si,2:O --method debye --rmax 12 --rdf g_of_r.dat S_q.dat

# GPU transform, single precision
sqcalc -i traj.dump -m 1:Si,2:O --device gpu --precision single S_q_gpu.dat
```

Things worth knowing before suggesting flags:

* The reciprocal methods sample $q$ on the reciprocal lattice of the dump box,
  so the box must not change between frames and a non-periodic cluster needs
  enough vacuum padding; the grid cost grows like $(q_{\max} L)^3$ and grids
  above 4e8 points are refused.  Debye needs no padding and wins for large,
  low-density boxes and high $q$; the reciprocal method wins for small dense
  periodic cells.
* `--weight neutron`/`xray` and the partials require `-m`; unknown type ids, or
  elements without tabulated data, are rejected with a clear message.
* Coordinates are taken as dumped: wrapped `x y z` is the native input, `xu yu
  zu` is wrapped internally.
* `--grid` is total-only, written single threaded, and unavailable with
  `--method debye`.  HDF5 output stores one one-dimensional dataset per
  quantity (e.g. `/shell/q`, `/shell/S`, `/grid/qx`, `/grid/S`), the partials
  under `/shell/S_partial/<pair>`, and the run metadata as file attributes.
* Progress goes to stderr and the table to `OUTPUT`, so `-q` keeps logs clean.

## Tests

`ctest --test-dir build --output-on-failure` runs `finufft_opts` (Fortran mirror
of the FINUFFT options struct), `neighbour_list` (Debye cell list against brute
force), `sqcalc_physics` (NUFFT versus direct summation and independent numpy
references, plus CLI error handling) and, on CUDA builds, `sqcalc_gpu` (GPU
against the CPU reference, skipping itself without a visible device).  A single
test: `ctest --test-dir build -R sqcalc_physics`.

## Reference

`README.md` in the checkout carries the theory (weighting, normalization, Debye
and partial structure factor formulas) and the build details; `AGENTS.md`
documents the module layout, coding style and test conventions to follow when
changing the code.

This skill is versioned with the code in `skills/sq-calc/`, so a fresh clone
picks it up from there.  To install it for another user or machine, copy the
folder into that machine's skills directory (`~/.agents/skills` or
`$CODEX_HOME/skills`), or let the skill installer fetch it:

```sh
install-skill-from-github.py --url \
  https://github.com/thatstar/sqcalc/tree/master/skills/sq-calc
```
