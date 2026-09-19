# sqcalc command line usage

```
sqcalc -i DUMP [options] OUTPUT
```

`OUTPUT` is the shell averaged $S(q)$ table and `-` writes it to stdout.  `-i`
is the only required option.  The table goes to `OUTPUT` and progress goes to
stderr, so `-q`/`--quiet` keeps logs clean.

## Options

| option | meaning |
| --- | --- |
| `-i, --input FILE` | LAMMPS dump trajectory (required) |
| `-m, --mapping LIST` | LAMMPS type id to element symbol, e.g. `1:Si,2:O` |
| `-w, --weight SCHEME` | `unit` (default), `neutron` or `xray` |
| `-t, --threads N` | OpenMP threads (default: all available) |
| `--qmin`, `--qmax`, `--nq` | $q$ range [1/A] and number of shells (defaults 0, 20, 500) |
| `--grid FILE` | also write $S(q)$ on every reciprocal lattice point |
| `--grid-format NAME` | `text` (default) or `hdf5`; a `.h5`/`.hdf5` name implies hdf5 |
| `--method NAME` | `nufft` (default), `direct` ($O(N N_{\text{modes}})$ reference) or `debye` (real space pair histograms) |
| `--device NAME` | `cpu` (default) or `gpu` (CUDA builds only) |
| `--gpu-id N` | CUDA device to use when `--device gpu` (default 0) |
| `--precision NAME` | `double` (default) or `single` (float32, GPU only) |
| `--norm NAME` | `mean` (default), `self` or `n` (see below) |
| `--eps VALUE` | NUFFT tolerance (default 1e-9; 1e-5 for the float32 GPU path) |
| `--partials`, `--no-partials` | write / omit the partial structure factor columns (default: on) |
| `-fz, --faber-ziman` | report the partials in the Faber-Ziman normalization |
| `--rmax VALUE` | Debye pair cutoff [A] (default: half the smallest periodic box side, or all pairs without periodicity) |
| `--dr VALUE` | Debye radial bin width [A] (default 0.01, reliable up to $q \sim \pi/(2\,dr)$) |
| `--skin VALUE` | Debye Verlet skin for reusing the pair list [A] (default 1.0, `0` rebuilds every frame) |
| `--rdf FILE` | total and partial $g(r)$ in one file (text, or HDF5 for `.h5`/`.hdf5`) |
| `--no-cutoff-correction` | disable the Debye cut-off density correction (applied by default when at least one direction is periodic) |
| `-q, --quiet` | suppress the progress output on stderr |
| `-h, --help` | option summary |
| `-v, --version` | program version |

## What the options mean

### Weighting

The per-atom weight $w_j$ is `1`, the bound coherent neutron scattering length
$b$ (`neutron`) or the IT92 form factor $f(q)$ (`xray`).  The weights come from
a table of 104 elements, so `-m` is required for anything but `unit`; without it
every atom has weight 1.0.

`S(q) = \langle |\rho(q)|^2 \rangle / W(q)` with `--norm` selecting $W(q)$:

| `--norm` | $W(q)$ | use |
| --- | --- | --- |
| `mean` (default) | $N \langle w \rangle^2$ | Faber-Ziman total $S(q)$ |
| `self` | $\sum_j w_j^2$ | $S(q) \to 1$ at large $q$ |
| `n` | $N$ | the debyer convention |

### Methods

* `nufft` - the production path: reciprocal space transform with FINUFFT (CPU)
  or cufinufft/cuFFT (`--device gpu`, CUDA builds), averaged over frames.
* `direct` - brute force $O(N N_{\text{modes}})$ summation over the same
  reciprocal lattice; a slow reference to cross-check results against.
* `debye` - real space pair histograms in the style of the debyer program.
  `--rmax`, `--dr`, `--skin` and `--rdf` apply here, `--grid` does not, and it
  runs on the CPU.  It is the cheaper choice for large, sparse or non-periodic
  systems, where the grid methods pay for $(q_{\max} L)^3$ lattice points.

`--precision single` switches only the GPU transform to float32 (coordinates and
strengths are converted down, the grid is widened back for the double precision
combine), which is worthwhile on consumer cards with a low FP64 rate; it differs
from float64 by roughly 1e-6 in $S(q)$.

### Choosing the q sampling

The reciprocal methods (`nufft`, `direct`) can only evaluate the reciprocal
vectors of the box, $|q| = |h b_1 + k b_2 + l b_3|$, so the sampling has to
follow the box.  `scripts/choose_q.py DUMP` reads the box from the first frame
and prints the option string to use; the reasoning is

* `--qmin` should be the smallest accessible $|q|$ - `2 pi / L` for a cubic
  box.  With the default `--qmin 0` every shell below it is empty, and those
  rows are not measurements.
* `--nq` sets the shell width, `dq = (qmax - qmin) / nq`.  The $|q|$ values of
  a lattice are not evenly spaced, and a shell that falls between two
  neighbours stays empty, in which case sqcalc writes a flat 0.  For the
  default `--nq 500` on a 1000-atom box at $rho = 1.2$ this hits a third of
  the shells.  A width of about $0.3$-$0.4 \times 2 pi / L$ is the practical
  limit before shells start to empty out, and the script reports the widest
  width that still populates every shell.
* A shell that holds only a few modes is noisy: the spread of a shell value
  falls roughly as $1/sqrt("modes" \times "frames")$, so the low-$q$ end is the
  first place to look when a curve looks ragged (the script prints the mode
  counts).
* `--qmax` costs $(q_{\max} L)^3$ lattice points; grids above 4e8 are refused,
  which caps `qmax` for large boxes.  The script reports the grid it would
  build.
* A small box cannot resolve small $q$ at all: use a larger box, or the Debye
  method, whose $q$ grid is not tied to the box.

### Partials

Each type pair gets a column
($S_{ab}(q) = \langle \sum_{j \in a} \sum_{k \in b} \sin(q r_{jk})/q r_{jk}
\rangle / N$, OVITO convention, including the $r = 0$ self term).  They satisfy
$S(q) = \sum_{ab} (2 - \delta_{ab}) S_{ab}(q)$ and tend to $S_{aa} \to x_a$,
$S_{ab} \to 0$ at large $q$.  `-fz` reports instead
$A_{ab}(q) = (S_{ab}(q) - x_a \delta_{ab})/(x_a x_b) + 1$, which tends to 1.
Partials work with both the reciprocal and the Debye method and are written by
default, so a model system without elements (`-w unit`, no `-m`) still gets
them: the pair is labelled with the LAMMPS type ids, `S(1-1) S(1-2) S(2-2)`,
and `-m` only replaces those labels with element symbols.  The derivations are
in the repository `README.md`.  **Prefer `-fz` when reporting partials**: the
Faber-Ziman form tends to 1 for every pair, so the curves share a common
asymptote and can be compared directly with the partials of other systems,
while the default (OVITO) partials tend to the concentrations $x_a$ and hide
the structure behind the composition.

## Output

The shell table is text:

```
# q S(q)
      0.100000    1.234567890123E+00
```

with one extra column per type pair when partials are written, e.g.
`# q S(q) S(Si-Si) S(Si-O) S(O-O)` with `-m 1:Si,2:O` or
`# q S(q) S(1-1) S(1-2) S(2-2)` without a mapping.

`--grid FILE` writes every reciprocal lattice vector of the requested range as
`# qx qy qz S(q)` (text), excluding the trivial $q = 0$ mode in both tables.
With a `.h5`/`.hdf5` name or `--grid-format hdf5` the result is HDF5 instead:
one one-dimensional dataset per quantity (`/shell/q`, `/shell/S`,
`/shell/count`, `/grid/qx`, `/grid/qy`, `/grid/qz`, `/grid/h`, `/grid/k`,
`/grid/l`, `/grid/S`), partials under `/shell/S_partial/<pair>` with the labels
in `/shell/pairs`, and the run metadata as file attributes (`nframes`, `natoms`,
`nmodes`, `nq`, `qmin`, `qmax`, `eps`, `cell`, `weight`, `norm`, `device`,
`precision`, `mapping`).  Writing large grids is far cheaper in HDF5 than in
text.

`--rdf FILE` (Debye method) holds one row per $r$ bin, `# r g(r) g(a-a) ...` in
the order of `-m`, with the partials normalized to $g_{ab} \to 1$ at large $r$
(HDF5: `/rdf/r`, `/rdf/g`, `/rdf/g_partial/<label>`, `/rdf/pairs`).  The
weighted total uses the $q \to 0$ amplitudes.

The shell, `g(r)` and grid tables are plain text, with the column names in the
header line, so any plotting tool reads them as they are.  The skill ships
`scripts/plot_sq.py` for this: it takes the axis labels and the partials from
the header and saves a figure (matplotlib is its only dependency).  HDF5 output
is meant for further analysis rather than plotting and needs h5py or the test
helper `h5read`.

## Examples

```sh
# unit weights, stdout, 4 threads
sqcalc -i traj.dump -t 4 - > S_q.dat

# neutron weighting of a two component glass, plus the reciprocal grid table
sqcalc -i traj.dump -m 1:Si,2:O -w neutron --qmax 25 --nq 1000 \
       --grid S_q_grid.dat S_q.dat

# X-ray weighting (q dependent form factors), cross-checked by direct summation
sqcalc -i traj.dump -m 1:Si,2:O -w xray --method direct --qmax 6 S_q_direct.dat

# Debye method with partial g(r), and a GPU transform
sqcalc -i traj.dump -m 1:Si,2:O --method debye --rmax 12 --rdf g_of_r.dat S_q.dat
sqcalc -i traj.dump -m 1:Si,2:O --device gpu --precision single S_q_gpu.dat
```

## Limits and pitfalls

* The simulation box must not change along the trajectory; the reciprocal
  methods build the lattice once from the first frame.
* The reciprocal methods sample $q$ on the lattice of the dump box, so a
  non-periodic cluster needs enough vacuum padding for $2\pi/L$ to resolve the
  structure, and the grid cost grows with the padding.  Debye removes the
  padding requirement.
* A finite Debye cutoff biases the low $q$ region, which is what the cut-off
  density correction removes; on uncorrelated frames `--skin 0` is faster.
* Coordinates are taken as dumped: `x y z` is native, `xu yu zu` is wrapped
  internally.
* With `-w neutron`/`xray` a type missing from the mapping, or an element
  without tabulated data, is rejected with a clear message; `-w unit` needs no
  mapping and labels the partials by type id instead.
* The shell table reports 0 for shells the box cannot sample.  Pick `--qmin`
  and `--nq` with `scripts/choose_q.py` before trusting the low-$q$ end, and
  read the Faber-Ziman partials (`-fz`) when comparing partials across
  concentrations.
* Grids above 4e8 points are refused (about 6 GB); `--grid` is written single
  threaded and does not affect the shell table.
