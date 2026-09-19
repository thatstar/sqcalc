# Repository Guidelines

## Project Structure & Module Organization

`sqcalc` is a Fortran 2008 / CMake project that computes structure factors from
LAMMPS dump trajectories.

- `app/main.f90` — command line driver; option parsing lives in
  `src/sqc_options.f90`.
- `src/` — one module per concern, prefix `sqc_`: `sqc_cell`, `sqc_dump`,
  `sqc_weights`, `sqc_elements`, `sqc_finufft`, `sqc_structure_factor`,
  `sqc_debye`, `sqc_hdf5`, `sqc_gpu`/`sqc_cufinufft`. New modules follow the
  `sqc_<topic>.f90` pattern and are added to `sqcalc_lib` in `CMakeLists.txt`.
- `test/` — Fortran unit tests and `test_sqcalc.py`, the end-to-end suite.
- `external/finufft/` — vendored FINUFFT (git subtree, do not edit).
- `skills/sq-calc/` — the agent skill shipped with the code: its
  `references/usage.md` and `references/building.md` document the CLI and the
  build, and `scripts/check_doc_sync.py` (run by the `skills_doc_sync` test)
  keeps them in step with `src/sqc_options.f90` and `CMakeLists.txt`.
- `.devdocs/` — local design notes (gitignored).

## Build, Test, and Development Commands

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release   # configure
cmake --build build -j                           # compile sqcalc
ctest --test-dir build --output-on-failure       # run the test suite
ctest --test-dir build -R sqcalc_physics         # single test
```

Useful options: `-DFFTW_ROOT=<prefix>` to point at FFTW3, and
`-DSQC_ENABLE_CUDA=ON` plus `-DSQC_CUDA_ARCHITECTURES=<cc>` for the GPU
backend (tested by `ctest -R sqcalc_gpu`). Debug builds add
`-fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`.

## Coding Style & Naming Conventions

- Fortran free form, three-space indent, `implicit none`, modules default to
  `private` with explicit `public` lists.
- Types end in `_t`, procedures use `module_name_verb` (e.g.
  `cell_init_from_bounds`), and kinds come from `sqc_kinds` (`rk`, `ik`).
- Start every source with the project banner, copyright line and
  `! SPDX-License-Identifier: GPL-3.0-or-later`; document exported routines
  with `!>` doxygen comments.
- Python tests follow PEP 8, four-space indent, `snake_case`.

## Testing Guidelines

New physics or I/O behaviour needs a test. Add Fortran unit tests as
`test/test_<subject>.f90` and register them with `add_test` in
`test/CMakeLists.txt`; extend `test/test_sqcalc.py` for end-to-end coverage,
cross-checking against the numpy references in `test/ref_sq.py` and
`test/debye_ref.py`. Test data is generated at run time by `test/gen_dump.py`;
do not commit `.dump` files.

## Commit & Pull Request Guidelines

Commit messages are short, imperative and capitalized, often with a module
scope prefix: `Debye: cross-check against the external debyer program`,
`Fix the HDF5 g(r) layout`. Keep each commit focused and describe physics or
output changes in the body. Pull requests should state the motivation, list the
commands run (`ctest --test-dir build`), note any result changes for existing
inputs, and reference related issues.
