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
| `--partials`, `--no-partials` | write / omit the partial structure factor columns (default: on when `-m` is given) |
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

### Partials

With a mapping, each type pair gets a column
($S_{ab}(q) = \langle \sum_{j \in a} \sum_{k \in b} \sin(q r_{jk})/q r_{jk}
\rangle / N$, OVITO convention, including the $r = 0$ self term).  They satisfy
$S(q) = \sum_{ab} (2 - \delta_{ab}) S_{ab}(q)$ and tend to $S_{aa} \to x_a$,
$S_{ab} \to 0$ at large $q$.  `-fz` reports instead
$A_{ab}(q) = (S_{ab}(q) - x_a \delta_{ab})/(x_a x_b) + 1$, which tends to 1.
Partials work with both the reciprocal and the Debye method.  The derivations
are in the repository `README.md`.

## Output

The shell table is text:

```
# q S(q)
      0.100000    1.234567890123E+00
```

with one extra column per type pair when partials are written, e.g.
`# q S(q) S(Si-Si) S(Si-O) S(O-O)`.

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
* Unknown type ids, or elements without tabulated X-ray/neutron data, are
  rejected with a clear message.
* Grids above 4e8 points are refused (about 6 GB); `--grid` is written single
  threaded and does not affect the shell table.
