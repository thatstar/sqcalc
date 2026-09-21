# sqcalc command line usage

```
sqcalc -i DUMP [options] OUTPUT
```

`OUTPUT` is the shell averaged $S(q)$ table and `-` writes it to stdout.  `-i`
is the only required option, and a `--dyn` run with `--dyn-q -` (the default
there) takes no `OUTPUT` argument at all, because it computes no $S(q)$.  The
table goes to `OUTPUT` and progress goes to stderr, so `-q`/`--quiet` keeps
logs clean.

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
| `--pair-entropy FILE` | total and partial pair entropy $S_2/k_B$ from the Debye $g(r)$ |
| `--s2-accum FILE` | $S_2(r)$ accumulation curve for tail extrapolation |
| `--no-cutoff-correction` | disable the Debye cut-off density correction (applied by default when at least one direction is periodic) |
| `--dyn` | keep the time axis: dynamic structure factor, four-point structure factor and overlap (see below) |
| `--dyn-q SPEC` | $q$ sampling of `--dyn`: `-` (default, no $q$ points), `line:NINT,S0,S1,DX,DY,DZ`, `shell:Q,ACC`, `grid:QMAX` or `single:N1,N2,N3` |
| `--dyn-modes N` | mode budget of `--dyn-q grid`: 0 = unlimited (default) |
| `--dyn-thin KIND` | order of the mode thinning: `shells` (default) or `orbits` |
| `--dyn-keep-modes` | write the grid rows per lattice vector instead of the `\|q\|` shell average |
| `--dt VALUE` | time step of the trajectory (the LAMMPS `timestep`), required by `--dyn` |
| `--maxframes L` | correlation window in frames (the largest lag kept), required by `--dyn` |
| `--lag N` | frames between two consecutive time origins (default 1) |
| `--sqw FILE` | the $S(q,\omega)$ spectra (text, or HDF5 for `.h5`/`.hdf5`) |
| `--fqt FILE` | the coherent intermediate scattering function $F(q,t)$ |
| `--fqt-self FILE` | the self intermediate scattering function $F_s(q,t)$ |
| `--dyn-format NAME` | `text` (default) or `hdf5` for all dynamic outputs |
| `--s4 FILE` | the total four-point structure factor $S_4(q,t)$ |
| `--chi4 FILE` | the average overlap $Q(t)$ and dynamic susceptibility $\chi_4(t)$ |
| `--s4-cutoff A` | overlap cutoff $a$ for `--s4` and `--chi4`, in dump length units |
| `--buffer-limit GB` | position buffer limit for S4/chi4/F_s (default 2.0; GB = $10^9$ bytes) |
| `--stride N` | use every N-th dump frame for S4/chi4/F_s (default 1) |
| `-q, --quiet` | suppress the progress output on stderr |
| `-h, --help` | option summary |
| `-v, --version` | program version |

## What the options mean

### Weighting

The per-atom weight $w_j$ is `1`, the bound coherent neutron scattering length
$b$ (`neutron`) or the IT92 form factor $f(q)$ (`xray`).  The weights come from
a table of 104 elements, so `-m` is required for anything but `unit`; without it
every atom has weight 1.0.

$$S(q) = \frac{\langle |\rho(q)|^2 \rangle}{W(q)}$$

with `--norm` selecting $W(q)$:

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

The Debye method is the exception to all of the above: it evaluates $S(q)$
directly from the pair histogram, not from a transform, so every $q$ is
available.  There are no empty shells and no lattice-imposed shell width - a
small $\Delta q$ just samples the same smooth curve at more points, at the cost
of one pass over the histogram bins per $q$ (negligible next to the pair list).
The limits are in real space instead: `--dr` sets the $q$ range, since the
histogram is only reliable up to $q \sim \pi/(2\,dr)$ (sqcalc prints a note
when `--qmax` exceeds it), and `--rmax` together with the cut-off correction
fixes the low $q$ end.  So `--nq` can be pushed as far as the table length is
convenient: `scripts/choose_q.py DUMP --method debye` suggests
`--qmin 0 --qmax 20 --nq 2000` ($\Delta q = 0.01$) and checks `--dr` against
the requested `--qmax`.

### Pair entropy (`--pair-entropy`, `--s2-accum`)

The Debye histograms also give the two-body excess entropy (pair entropy) per
particle in units of $k_B$.  The partial contribution is

$$
S_2^{ab} = -2\pi\rho\,x_a x_b \int_0^{r_{\max}} r^2
\left[g_{ab}\ln g_{ab} - g_{ab} + 1\right]dr
$$

and the total follows the usual partial sum rule,

$$
S_2 = \sum_a S_2^{aa} + 2\sum_{a<b}S_2^{ab}.
$$

The integration uses the exact shell volume of each radial bin (not
$4\pi r^2 dr$), takes the $g\to0$ limit of the integrand as 1, and symmetrizes
the cross partials before taking the logarithm.  A leading-order Poisson bias
correction removes the noise-induced offset of the nonlinear integrand, so
even a noisy ideal gas stays close to $S_2=0$; no smoothing is applied by
default.  No statistical error is estimated: $S_2$ is used as a
qualitative/trend quantity, and the dominant uncertainties are the
$r_{\max}$ tail, the finite box and the multiparticle corrections.

`--s2-accum` writes the $S_2(r)$ accumulation curve instead of only the final
value at $r_{\max}$.  A single Debye run cannot make a rigorous tail
correction (the $g(r)$ data stop at half the periodic box, and the Debye
$S(q)$ is not suitable for the reciprocal-space route), so the model-dependent
analysis is left to `scripts/s2_analysis.py`: it can GCV-smooth the partial
$g(r)$ (with `--rdf`), apply the `--dr` Richardson extrapolation, and fit the
tail with power-law, exponential or damped-oscillation forms.  Do not use the
Debye $S(q)$ to estimate the reciprocal-space entropy; use the NUFFT/direct
method for that.

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

### The dynamic run (`--dyn`, `--dyn-q`)

The other methods average the snapshots as an unordered ensemble and return
$S(q)$; `--dyn` keeps the time axis instead.  It is a plain flag: `--dt` and
`--maxframes` set the time axis, and `--dyn-q SPEC` chooses how reciprocal
space is sampled,

| `--dyn-q` | q points | outputs |
| --- | --- | --- |
| `-` (default) | none | only `--chi4` ($Q(t)$ and $\chi_4(t)$) |
| `line:NINT,S0,S1,DX,DY,DZ` | `NINT+1` points on a line through Gamma | all |
| `shell:Q,ACC` | one shell $\|q\| = Q$, averaged over every direction | all |
| `grid:QMAX` | every reciprocal-lattice vector with $\|q\| \le Q_{\max}$ | all |
| `single:N1,N2,N3` | one reciprocal-lattice vector of the box | all |

Without q points there is no $S(q)$ at all: the run takes no `OUTPUT` table,
and `--sqw`, `--fqt`, `--fqt-self` and `--s4` are rejected.  That is the
cheapest way to get $Q(t)$ and $\chi_4(t)$, which only need the overlap of the
displaced positions.

#### A q line

`line:NINT,S0,S1,DX,DY,DZ` are `NINT` intervals (so `NINT+1` q points) with the
scale running from `S0` to `S1` in 1/A along the direction `(DX,DY,DZ)`; every
line passes through the Gamma point.  `--dyn-q line:100,0.5,20,1,1,0` is
therefore 101 q points from 0.5 to 20 1/A along (1,1,0).  Off-lattice $q$ is
the point of the method: it is what an experiment at that $q$ measures, and
the reciprocal grid methods cannot reach it.

#### A q shell

`shell:Q,ACC` puts every direction of the sphere $|q| = Q$ on a Lebedev grid
and averages the results,

$$
\bar X(Q) = \frac{\sum_k w_k X(q_k)}{\sum_k w_k},
$$

with the standard Lebedev weights $w_k$ (they sum to $4\pi$).  `ACC` selects
the rule: `low` is the order 11 rule with 50 directions, `medium` order 17 with
110, and `high` order 23 with 194.  That is enough to integrate a smooth
$S(q,\omega)$ over the sphere to better than the statistical noise, and it
gives the isotropic average that a powder or a simulation with no preferred
direction actually measures.  The static table, $F(q,t)$, $S(q,\omega)$,
$S_4(q,t)$ and $F_s(q,t)$ are all shell averages; each is written as a single
row at $q = (0,0,Q)$, and the static table is labelled by $|q| = Q$.  The
average is taken over the *normalized* $S(q)$ of every direction, which matters
for `-w xray`, where the form factor and hence the denominator depend on $q$
even at fixed $|q|$.

#### A reciprocal-lattice grid

`grid:QMAX` samples the reciprocal lattice of the dump box instead of an
arbitrary line or sphere,

$$
\mathbf q = n_1 \mathbf b_1 + n_2 \mathbf b_2 + n_3 \mathbf b_3,
\qquad 0 < |\mathbf q| \le Q_{\max},
$$

plus the Gamma point.  This is the sampling to use when $S_4(q,t)$ is fitted to
the Ornstein-Zernike form: only a reciprocal-lattice vector carries a density
amplitude that does not depend on how the periodic images are chosen, so the
grid is free of the box-form-factor contamination that makes the off-lattice
Lebedev shell unusable for $S_4$.  The tables hold one row per lattice shell
and the Gamma row is the $q \to 0$ limit, i.e. $\chi_4(\tau)$.

The modes are reduced in three steps, in this order:

* `W(-q) = conjg(W(q))` because the overlap weights are real, so one vector of
  every $\pm$ pair is kept; this halves the work and loses nothing;
* vectors that are related by the point group of the periodic lattice have the
  same expectation for an isotropic system, so `--dyn-thin orbits` may keep one
  representative per orbit.  The orbit multiplicity is kept, because a shell
  can hold several orbits (for $|\mathbf n|^2 = 9$, a 6-fold $(3,0,0)$ orbit and
  a 24-fold $(2,2,1)$ one);
* `--dyn-thin shells`, the default, drops whole shells instead, which leaves the
  remaining shells complete and therefore costs no accuracy in the points that
  survive.  The kept shells are spread evenly over the q range and always
  include its two ends.

The tables are already the isotropic average: one row per lattice shell, at
$q = (0,0,|q|)$, taken over the reciprocal-lattice vectors of that $|q|$ and
weighted by their orbit multiplicities.  That is the object an
Ornstein-Zernike fit and a powder average want, and the code does the weighting
because a reader of the table cannot: the multiplicity of a vector such as
$(1,2,2)$ follows from the point group of the cell, not from the row.
`--dyn-keep-modes` writes one row per lattice vector instead, which is what a
directional analysis and the verification below need.

These shells are the exact ones, i.e. all lattice vectors of one $|\mathbf q|$.
The static grid under `--method nufft` bins $|\mathbf q|$ uniformly instead,
with a width the user sets through `--nq`, because it reads the whole transform
at once and has no per-mode cost to control.

`--dyn-modes N` caps how many modes the run pays for; without it the cell sets
the count, which grows as $|q_{\max} L|^3$.  The budget is a target rather than
a hard limit: when even the two end shells need more modes than it allows, they
are kept anyway and the run summary says so.  Both the shell and the orbit
reduction assume that the system is isotropic and in equilibrium.  The cheapest
way to check that is built in: a `--dyn-q grid` run with `--s4` prints an
`isotropy` line in its summary with the spread of the modes of the smallest
shell at the lag of the $\chi_4$ peak, next to their mean value.  Compare that
spread with the run-to-run scatter of a single mode (a second run with a
different `--lag`, or the two halves of the trajectory): a spread that stays
far above it, and that does not shrink when the box is enlarged, means the
trajectory is not equilibrated or the system is not isotropic, and then
neither the orbit reduction nor an isotropic correlation length may be
trusted.  The per-mode numbers behind it are in the `S4` table of a
`--dyn-keep-modes` run.

#### One lattice vector

`single:N1,N2,N3` samples exactly one reciprocal-lattice vector of the box,

$$
\mathbf q = N_1 \mathbf b_1 + N_2 \mathbf b_2 + N_3 \mathbf b_3,
$$

with integer indices that may be negative.  This is the exact way to ask for a
single lattice mode: the vector is built from the integers, so it sits on the
lattice however the box is oriented, whereas a hand written $|q|$ on a `line:`
would be a few ulps away, and a density amplitude that is even slightly off the
lattice stops being independent of the choice of the periodic images.  The
table has one row, labelled by $|\mathbf q|$, while the header carries the
indices and the vector.  `single:0,0,0` is the Gamma point, whose $S_4$ row is
$\chi_4(t)$; a run that wants only $Q(t)$ and $\chi_4(t)$ is cheaper with
`--dyn-q -` and `--chi4`.

The dynamics come from the same density amplitudes $\rho_a(q,t)$, evaluated by
direct summation at the sampled $q$ and correlated in time with a multi-origin
estimator over a ring buffer,

$$C_{ab}(q,\tau) = \langle \rho_a(q,t+\tau)\rho_b^*(q,t)\rangle,$$

accumulating a sum *and* a count per lag, so the memory does not grow with the
trajectory length and a trajectory whose length is not a multiple of the window
needs no padding.  `--maxframes L` is the window: the largest lag kept, and
hence the resolution $\Delta\omega = \pi/(L\,\Delta t_{\text{frame}})$.  `--lag N`
is the distance in frames between consecutive time origins (default 1, the most
overlapping averages; larger values are cheaper and more nearly independent).

`--dt` is the time step of the trajectory, i.e. the LAMMPS `timestep` value:
the frame interval is `dt` times the step increment found in the dump, which
must be constant (a dump written at an irregular cadence is rejected).  The
unit of $\omega$ is the inverse of the unit of `--dt` (a `dt` in ps gives
$\omega$ in rad/ps); the trajectory carries no unit, so it is up to the user to
keep `--dt` and the interpretation consistent.

The transform uses the one-sided folding $S(q,-\omega) = S(q,\omega)$ and the
normalization $W(q)$ of the static table, so the zeroth moment is exact:
$\int S(q,\omega)\,d\omega = F(q,0) = S(q)$, the value in `OUTPUT`.  Chemical
weighting (`-w neutron`/`xray`, `--norm`) is applied to the partial correlation
functions exactly as for the static partials, so the same conventions hold.

### Four-point structure factor (`--s4`, `--chi4`)

`--s4` adds the total four-point structure factor to the same dynamic run,
and `--chi4` writes the average overlap and the dynamic susceptibility.  Both
need the overlap cutoff,

```sh
--s4-cutoff A
```

in the length unit of the dump.  The cutoff is a fraction of the particle
diameter; a common choice is a value near the first minimum of $g(r)$, for
example $0.3\,\sigma$ for a Lennard-Jones glass.  There is no default: the
dump does not carry the particle size.

The overlap of atom $i$ between the origin $t_0$ and the later time
$t_0+t$ is

$$
w_i(t,t_0) = \Theta\!\left(a - |\mathbf r_i(t_0+t) - \mathbf r_i(t_0)|\right),
$$

and the four-point structure factor is the connected fluctuation of its
origin-position Fourier transform,

$$
W(\mathbf q;t_0,t) = \sum_i e^{i\mathbf q\cdot\mathbf r_i(t_0)}\, w_i(t,t_0),
$$

$$
\begin{aligned}
S_4(\mathbf q,t) = \frac{1}{N}\Bigl[
  &\langle W(\mathbf q;t_0,t) W(-\mathbf q;t_0,t)\rangle_{t_0} \\
  &- \left|\langle W(\mathbf q;t_0,t)\rangle_{t_0}\right|^2
\Bigr].
\end{aligned}
$$

`--chi4` writes the scalar $q = 0$ case,

$$
Q(t) = \frac{1}{N}\left\langle \sum_i w_i(t,t_0)\right\rangle_{t_0},
$$

$$
\chi_4(t) = \frac{1}{N}\left(\langle W_0(t)^2\rangle-\langle W_0(t)\rangle^2\right),
$$

with $W_0(t) = \sum_i w_i(t,t_0)$.  This is the same quantity as the $q \to 0$
limit of $S_4(q,t)$,

$$
\lim_{q\to0} S_4(q,t) = \chi_4(t),
$$

and a small-$q$ Ornstein-Zernike fit,

$$
S_4(q,t) \simeq \frac{\chi_4(t)}{1+[q\,\xi(t)]^2},
$$

gives the dynamic correlation length $\xi(t)$.  $Q(0) = 1$ and
$\chi_4(0) = 0$; $\chi_4(t)$ peaks near the structural relaxation time and
both functions decay at long times.

The estimator is the unbiased connected covariance over time origins, with
the same origin count per lag as the coherent correlation.  A lag with a
single time origin is written as zero; the HDF5 `count` dataset shows how
many origins contributed.  The overlap uses unit weights, so `-w` and
`--norm` affect $S(q)$, $F(q,t)$ and $S(q,\omega)$ but not $S_4$ or
$\chi_4$.  Both are total quantities: `--partials` does not add partial
four-point columns.

`--chi4` can be used without `--s4`.  It still needs the position ring buffer
to evaluate the overlap, but it skips the q-resolved sums, so it is much
cheaper.  The full $S_4$ calculation is dominated by the per-atom displacement
check and the phase sum over the atoms inside the cutoff; a large `--lag`
(more widely spaced time origins) and a short, low-$q$ line are the practical
ways to control the cost.  The position buffer grows with the effective S4
window (`--maxframes` divided by `--stride`, rounded down);
`--buffer-limit GB` sets its limit (default 2.0, GB = $10^9$ bytes) and an
oversized request is rejected before any large allocation.

`--stride N` subsamples the S4/chi4 trajectory: only every N-th dump frame
contributes, the lag axis is $0, N, 2N, \ldots$ dump frames, and the position
buffer stores only those frames.  This reduces the buffer and the S4/chi4
work by roughly a factor $N$, at the cost of a coarser time axis.  The
coherent $F(q,t)/S(q,\omega)$ outputs still use every dump frame and the full
`--maxframes`; the S4/chi4 window is reduced to the largest multiple of $N$
that does not exceed `--maxframes`, and the run summary prints a note when
that happens.  `--lag` is still counted in dump frames, so with stride $N$,
`--lag m` selects origins every $m$ dump frames, not every $m$ S4 samples.

The dynamic method allocates only what the requested outputs need.  A run with
only `--s4`/`--chi4` does not allocate the coherent ring buffer, the
multi-origin correlation, or the $S(q,\omega)$/Fourier transform buffers; the
static `OUTPUT` table is still written.  A run with only `--sqw`/`--fqt` does
not allocate the S4 position buffer.  `--fqt-self` allocates that shared
position buffer but not the coherent $F(q,t)/S(q,\omega)$ buffers.  With
`--dyn-q -` there are no density amplitudes at all, so only `--chi4` is left
and no `OUTPUT` table is taken.

```sh
# low-q S4(q,t) and chi4(t), overlap cutoff a = 1.0 in dump length units
sqcalc -i traj.dump -w unit --dyn --dyn-q line:20,0.5,4,1,1,0 \
       --dt 0.005 --maxframes 400 \
       --s4-cutoff 1.0 --s4 S4.dat --chi4 chi4.dat S_q.dat

# the same overlap dynamics without any q points, and no S(q) table
sqcalc -i traj.dump -w unit --dyn --dt 0.005 --maxframes 400 \
       --s4-cutoff 1.0 --chi4 chi4.dat

# the isotropic average on one shell: 110 Lebedev directions at |q| = 2.5 1/A
sqcalc -i traj.dump -w unit --dyn --dyn-q shell:2.5,medium \
       --dt 0.005 --maxframes 400 --sqw S_qw.dat S_q.dat
```

### Self intermediate scattering function (`--fqt-self`)

`--fqt` is the coherent intermediate scattering function; `--fqt-self` writes
the self/incoherent counterpart

$$
F_s(q,\tau)=\frac{1}{N}\left\langle \sum_i
  e^{i q\cdot[r_i(t+\tau) - r_i(t)]}\right\rangle,
$$

with unit weights.  The per-species columns are

$$
F_s^{a}(q,\tau)=\frac{1}{N}\left\langle \sum_{i\in a}
  e^{i q\cdot[r_i(t+\tau) - r_i(t)]}\right\rangle,
$$

and the total is their sum,

$$
F_s(q,\tau)=\sum_a F_s^{a}(q,\tau).
$$

At $\tau = 0$ the total is 1 and each species column is $N_a/N = x_a$; the
columns are labelled by element symbol when `-m` is given, otherwise by the
LAMMPS type id (`F_s(Si)`, `F_s(O)` or `F_s(1)`, `F_s(2)`).

`--fqt-self` always includes every atom and **never** applies `--s4-cutoff`,
even when `--s4` is also requested.  It reuses the shared position buffer and
the `--stride`/`--lag` schedule: only every N-th dump frame contributes, lag
$j$ is $jN$ dump frames, and the effective window is the largest multiple of
$N$ that does not exceed `--maxframes`.  The position buffer limit
`--buffer-limit GB` applies to the same buffer.  `-w`/`--norm` do not affect
$F_s$.

The text columns are

```
# qx qy qz tau F_s(q,t) F_s(Si) F_s(O)
```

and HDF5 writes `/fqt_self/q`, `/fqt_self/tau`, `/fqt_self/F_s`,
`/fqt_self/count`, `/fqt_self/pairs` and
`/fqt_self/F_s_partial/<species>`, with `stride` and `effective_maxframes`
attributes.  The coherent `--fqt` output uses the `/fqt` group.

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

`--pair-entropy FILE` (Debye method) holds the final total and partial
$S_2/k_B$, one `total` row followed by one row per type pair, with the
`natoms`/`volume`/`nframes`/`rmax`/`dr`/`counts` metadata in the header.
`--s2-accum FILE` holds the accumulation curve
`# r S2_total(r)/kB S2(a-a)(r)/kB ...`.  Both accept a `.h5`/`.hdf5` name;
the HDF5 group is `/pair_entropy` with `S2`, `S2_total`, `pairs`, and, for the
accumulation file, `r`, `S2_curve` and `S2_total_curve`.

`--sqw FILE` (dynamic method) holds one row per $(q,\omega)$ sample,
`# qx qy qz omega S(q,w) S(a-a) ...`, and `--fqt FILE` the intermediate
scattering function with `tau` in place of `omega` and the same columns.  Both
accept a `.h5`/`.hdf5` name (or `--dyn-format hdf5`) for HDF5, where the
datasets are `/sqw/q` (one row per q: `qx qy qz |q|`), `/sqw/omega`,
`/sqw/S`, `/sqw/count` (the time origins actually accumulated per lag),
`/sqw/pairs` and `/sqw/S_partial/<pair>`, with the same layout under `/fqt`.
The partial spectra are always in the plain (OVITO) convention of the static
partial columns: transforming to the Faber-Ziman form would need the self part,
which the method does not separate, so `-fz` applies to the static table only.

`--s4 FILE` (dynamic method) holds one row per $(q,\tau)$ sample,
`# qx qy qz tau S4(q,t)`, with no partial columns.  `--chi4 FILE` is the
scalar time series `# tau Q(t) chi4(t)`.  Both accept a `.h5`/`.hdf5` name (or
`--dyn-format hdf5`) for HDF5, where the datasets are `/s4/q`, `/s4/tau`,
`/s4/S4`, `/s4/count` and `/chi4/tau`, `/chi4/Q`, `/chi4/chi4`,
`/chi4/count`, with the overlap cutoff in the `overlap` attribute.  The count
is the number of time origins at each lag; a lag with a single origin has no
measurable fluctuation and is written as zero.

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

# dynamic structure factor along (1,1,0), 101 q points, dumped every 10 steps
# of dt = 0.005 (frame interval 0.05), window 400 frames
sqcalc -i traj.dump -m 1:Si,2:O -w neutron \
       --dyn --dyn-q line:100,0.5,20,1,1,0 --dt 0.005 --maxframes 400 \
       --sqw S_qw.dat S_q.dat
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
  internally.  The dynamic method always needs *unwrapped* coordinates: it
  prefers `xu yu zu` and otherwise reconstructs them from `x y z` while
  reading, which requires the atom order (`id`) to be stable between frames.
* `--dyn` needs a dump written at a constant timestep interval, and enough of
  them: the frame interval sets the highest frequency the spectrum can see
  ($\omega_{\max} = \pi/\Delta t_{\text{frame}}$), so a trajectory dumped every
  100 MD steps resolves far less than one dumped every 10.  sqcalc prints the
  derived frame interval and warns when the frame-to-frame displacements
  suggest the dynamics are undersampled.
* With `-w neutron`/`xray` a type missing from the mapping, or an element
  without tabulated data, is rejected with a clear message; `-w unit` needs no
  mapping and labels the partials by type id instead.
* The `OUTPUT` table reports 0 for the grid shells the box cannot sample (the
  `--dyn-q shell:` average is a single shell and is never empty).  Pick
  `--qmin` and `--nq` with `scripts/choose_q.py` before trusting the low-$q$
  end, and read the Faber-Ziman partials (`-fz`) when comparing partials across
  concentrations.
* Grids above 4e8 points are refused (about 6 GB); `--grid` is written single
  threaded and does not affect the shell table.
