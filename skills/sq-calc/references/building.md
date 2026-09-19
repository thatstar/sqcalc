# Building sqcalc

## Requirements

* a Fortran 2008 compiler (gfortran >= 10), CMake >= 3.20 and OpenMP
* an FFTW3 installation (library and headers); CMake looks in `FFTW_ROOT`,
  `$CONDA_PREFIX`, any `~/.pixi/envs/*` and the usual system prefixes
* FINUFFT, vendored in `external/finufft` (git subtree, v2.5.1) and built
  together with sqcalc, so nothing has to be downloaded by hand

## Configure, build, test

From the checkout root (the directory containing `CMakeLists.txt`):

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The build produces `build/sqcalc`.

## CMake options

| option | effect |
| --- | --- |
| `-DFFTW_ROOT=PREFIX` | FFTW3 prefix; otherwise CMake probes `$FFTW_ROOT`, `$CONDA_PREFIX`, `~/.pixi/envs/*`, `/usr/local` and `/usr` |
| `-DSQC_ENABLE_CUDA=ON` | adds the cufinufft/cuFFT backend (`--device gpu`); also pass `-DSQC_CUDA_ARCHITECTURES=<cc>` (e.g. `75`) when no GPU is visible at configure time |
| `-DSQC_ENABLE_HDF5=OFF` | drops the HDF5 grid/`g(r)` output (on by default when HDF5 is found) |
| `-DSQC_USE_SYSTEM_FINUFFT=ON` | link an installed FINUFFT instead of the vendored subtree |
| `-DSQC_FINUFFT_ROOT=PREFIX` | prefix of that installed FINUFFT |
| `-DCPM_SOURCE_CACHE=DIR` | where FINUFFT's build helpers (`CPM.cmake`, `findFFTW`, `xsimd`, CCCL, possibly FFTW3) are cached; defaults to `.cpm-cache/` in the source tree |
| `-DSQC_BUILD_TESTS=OFF` | skip the test targets |
| `-DSQC_ENABLE_RUNTIME_CHECKS=ON` | add `-fcheck=bounds -fbacktrace` to non Debug builds; off by default because the checks are a development aid (they catch out of bounds accesses but cost a few percent and abort the run) |
| `-DCMAKE_BUILD_TYPE=Debug` | adds `-fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow` |

## Network and caching

The first configure of the vendored FINUFFT downloads its build helpers
(`CPM.cmake`, `findFFTW`, `xsimd` and, unless a system FFTW is found, FFTW3),
and a CUDA build additionally downloads CCCL, even though CUDA 12 ships its own
copy.  Everything is cached, so later configures work offline; point
`-DCPM_SOURCE_CACHE` at a shared cache, or `-DCPM_CCCL_SOURCE` at a local CCCL,
to avoid the download.  On a machine without network access, configure with a
populated cache, with `-DFFTW_ROOT`, or pass an existing build tree around.

A CUDA build links cufinufft as a shared library, which also makes the CPU
`libfinufft` shared.

## Tests

`ctest --test-dir build --output-on-failure` runs

* `finufft_opts` - the Fortran mirror of `finufft_opts` against the linked
  FINUFFT library,
* `neighbour_list` - the Debye cell list (pair enumeration, periodic images,
  skin reuse) against brute force pair counting,
* `sqcalc_physics` - NUFFT versus direct summation and independent numpy
  references for the weights, ideal gas, a simple cubic lattice, the grid
  output, the Debye method, $g(r)$, the partials, a triclinic box with `xu yu
  zu` columns, HDF5 and command line error handling,
* `sqcalc_gpu` (CUDA builds) - the GPU backend against the CPU, skipping itself
  when no device is visible.
* `skills_doc_sync` - the documents in `skills/sq-calc/references/` against
  `src/sqc_options.f90` and the CMake options above, via
  `tools/check_doc_sync.py` (registered only when Python is present).  The
  checker is a development tool and does not ship with the skill.

A single test: `ctest --test-dir build -R sqcalc_physics`.  Test data is
generated at run time by `test/gen_dump.py`, so no trajectories are committed.

## Repository conventions

`AGENTS.md` in the checkout documents the module layout, coding style and the
conventions for adding new modules and tests.
