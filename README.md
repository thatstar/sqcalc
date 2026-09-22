# sqcalc

`sqcalc` computes the total structure factor $S(q)$ of a LAMMPS trajectory.  It
reads the default LAMMPS `atoms` dump style, averages over all snapshots and
writes a `# q S(q)` table.  Two independent evaluations are available: the
reciprocal space transform with a non-uniform FFT (FINUFFT) on the CPU, or on a
CUDA GPU with cuFFT, which averages $|\rho(q)|^2$, and the real space Debye
pair-histogram method.  The table also carries the partial structure factors
$S_{ab}(q)$, one column per pair of LAMMPS types.

With `--dyn` the time axis is kept instead of being averaged away.  `--dyn-q`
chooses what is sampled in reciprocal space - a line, a single spherical shell,
or nothing at all - and the run writes the dynamic structure factor
$S(q,\omega)$, optionally with the intermediate scattering function $F(q,t)$.
The same dynamic run can also write the four-point structure factor $S_4(q,t)$,
the average overlap $Q(t)$ and the dynamic susceptibility $\chi_4(t)$ with
`--s4` and `--chi4`, and the mean squared displacement $\text{MSD}(t)$ with
`--msd`.  That method correlates the density amplitudes (or the overlap field)
in time rather than averaging their squares.

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
as one extra column per type pair, labelled with the species of `-m`
(e.g. `# q S(q) S(Si-Si) S(Si-O) S(O-O)`) or, without a mapping, with the
LAMMPS type ids (`# q S(q) S(1-1) S(1-2) S(2-2)`).

| option                           | meaning                                                                                                                          |
| -------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| `-i, --input FILE`             | LAMMPS dump trajectory (required)                                                                                                |
| `-m, --mapping LIST`           | LAMMPS type id to element or species, e.g.`1:Si,2:O` or `1:Si4+,2:O2-`                                                    |
| `-w, --weight SCHEME`          | `unit` (default), `neutron` or `xray`                                                                                      |
| `-t, --threads N`              | OpenMP threads (default: all available)                                                                                          |
| `--qmin`, `--qmax`, `--nq` | $q$ range and number of shells (defaults 0, 20 1/A, 500)                                                                       |
| `--grid FILE`                  | also write$S(q)$ on every reciprocal lattice point                                                                             |
| `--grid-format NAME`           | `text` (default) or `hdf5`; a `.h5`/`.hdf5` name implies hdf5                                                            |
| `--xrd FILE`                   | also write a powder XRD pattern;`.h5`/`.hdf5` implies HDF5                                                              |
| `--xrd-lambda VALUE`           | incident wavelength of`--xrd`, in the dump length unit                                                                    |
| `--xrd-range MIN MAX`          | two-theta range of`--xrd` in degrees (default `1 179`)                                                                    |
| `--xrd-step DEG`               | two-theta bin width in degrees (default: derived from the box)                                                            |
| `--lp`                         | apply the Lorentz-polarization factor of`--xrd` (default)                                                                 |
| `--no-lp`                      | drop the Lorentz-polarization factor of`--xrd`                                                                            |
| `--method NAME`                | `nufft` (default), `direct` ($O(N \cdot N_{\text{modes}})$ reference) or `debye` (real space pair histograms, see below) |
| `--device NAME`                | `cpu` (default) or `gpu` (needs `-DSQC_ENABLE_CUDA=ON`)                                                                    |
| `--gpu-id N`                   | CUDA device to use when`--device gpu` (default 0)                                                                              |
| `--precision NAME`             | `double` (default) or `single` (float32, GPU only)                                                                           |
| `--norm NAME`                  | `mean` (default), `self` or `n`                                                                                            |
| `--eps VALUE`                  | NUFFT tolerance (default 1e-9; 1e-5 for the float32 GPU path)                                                                    |
| `--partials`                   | append the partial structure factor columns (default: on)                                                                        |
| `--no-partials`                | do not append the partial columns                                                                                                |
| `-fz, --faber-ziman`           | report the partials in the Faber-Ziman normalization                                                                             |
| `-q, --quiet`                  | suppress progress output on stderr                                                                                               |
| `-h, --help`                   | show the option summary                                                                                                          |
| `-v, --version`                | print the program version                                                                                                        |

With `--method debye` these options select the pair histogram:

| Debye option               | meaning                                                                                                                                                   |
| -------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--rmax VALUE`           | pair cutoff [A]; default: half of the smallest periodic box side (minimum image convention, so no self images), or all pairs when the box is not periodic |
| `--dr VALUE`             | radial bin width (default 0.01 A, reliable up to$q \sim \pi/(2\,dr)$; ~150 1/A at the default)                                                          |
| `--skin VALUE`           | Verlet skin for reusing the pair list (default 1.0 A,`0` rebuilds every frame)                                                                          |
| `--rdf FILE`             | total and all partial$g(r)$ in one file (text, or HDF5 for `.h5`/`.hdf5`)                                                                           |
| `--pair-entropy FILE`    | total and partial pair entropy $S_2/k_B$ from the Debye $g(r)$                                                                                     |
| `--s2-accum FILE`        | the $S_2(r)$ accumulation curve for tail extrapolation                                                                                              |
| `--no-cutoff-correction` | disable the cut-off density correction (for comparison; applied by default when at least one direction is periodic)                                       |

With `--dyn` these options select the $q$ sampling and the correlation window:

| dynamic option                | meaning                                                                     |
| ----------------------------- | --------------------------------------------------------------------------- |
| `--dyn`                     | keep the time axis: dynamic structure factor, S4, overlap and MSD (see below) |
| `--dyn-q SPEC`              | q sampling: `-` (default, no q points), `line:NINT,S0,S1,DX,DY,DZ`, `shell:Q,ACC` or `grid:QMAX` |
| `--dyn-modes N`             | mode budget of `--dyn-q grid` (0 = unlimited, the default)                  |
| `--dyn-thin KIND`           | order of the grid thinning: `shells` (default) or `orbits`                  |
| `--dyn-keep-modes`          | write the grid rows per lattice vector instead of the $\|q\|$ shell average  |
| `--dt VALUE`                | time step of the trajectory (the LAMMPS`timestep`), required by `--dyn` |
| `--maxframes L`             | correlation window in frames (largest lag kept), required by`--dyn`       |
| `--lag N`                   | frames between consecutive time origins (default 1)                         |
| `--sqw FILE`                | the$S(q,\omega)$ spectra (text, or HDF5 for `.h5`/`.hdf5`)            |
| `--fqt FILE`                | the coherent intermediate scattering function$F(q,t)$                     |
| `--fqt-self FILE`           | the self intermediate scattering function$F_s(q,t)$                       |
| `--dyn-format NAME`         | `text` (default) or `hdf5` for all dynamic outputs                      |
| `--s4 FILE`                 | the total four-point structure factor$S_4(q,t)$                           |
| `--chi4 FILE`               | the average overlap$Q(t)$ and the susceptibility $\chi_4(t)$            |
| `--msd FILE`                | the mean squared displacement $\text{MSD}(t)$ and its species columns      |
| `--s4-cutoff A`             | overlap cutoff$a$, required by `--s4` and `--chi4`                    |
| `--buffer-limit GB`         | position buffer limit for S4/chi4/F_s/MSD (default 2.0 GB)                  |
| `--stride N`                | use every N-th dump frame for S4/chi4/F_s/MSD (default 1)                   |

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
       --dyn --dyn-q line:100,0.5,20,1,1,0 --dt 0.005 --maxframes 400 \
       --sqw S_qw.dat S_q.dat

# the same run averaged over a spherical shell of |q| = 2.5 1/A
sqcalc -i traj.dump -m 1:Si,2:O -w neutron \
       --dyn --dyn-q shell:2.5,medium --dt 0.005 --maxframes 400 \
       --sqw S_qw.dat S_q.dat

# and the overlap dynamics alone: no q points, no S(q) table
sqcalc -i traj.dump -w unit --dyn --dt 0.005 --maxframes 400 \
       --s4-cutoff 1.0 --chi4 chi4.dat

# the mean squared displacement, also without q points
sqcalc -i traj.dump -m 1:Si,2:O --dyn --dyn-q - --dt 0.005 --maxframes 400 \
       --msd msd.dat
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
* `xray`   : $w = f(\text{species}, q) = \sum_i a_i \exp(-b_i (q/4\pi)^2) + c$ (IT92).

The X-ray weights are tabulated per species, so the mapping of a type may name
an ion or a valence state instead of the neutral atom (`-m 1:Si4+,2:O2-`); the
neutron length is nuclear and uses the element of the species either way.  The
available labels are listed in `skills/sq-calc/references/usage.md`.

$S(q)$ is the trajectory average

$$
S(q) = \frac{\langle |\rho(q)|^2 \rangle}{W(q)}
$$

where the normalization is

| `--norm`         | $W(q)$                                                  |
| ------------------ | --------------------------------------------------------- |
| `mean` (default) | $N \langle w \rangle^2$, the Faber-Ziman total $S(q)$ |
| `self`           | $\sum_j w_j^2$, so that $S(q) \to 1$ at large $q$   |
| `n`              | $N$ (the convention used by the debyer program)         |

The weights are looked up from a table of 104 elements covering the periodic
table.  Without `-m` no element is known, so every atom has weight 1.0; a
weighted scheme therefore requires a type id mapping.

The three routes trade memory for time, and every run prints the memory of the
one it is about to use - the grid estimate, or the mode table of `direct` -
before it reads the trajectory:

| route | memory | work per frame |
| --- | --- | --- |
| `nufft` | a grid of $(q_{\max} L)^3$ lattice points, about 0.16 kB per point at the default tolerance | $O(N + M \log M)$ |
| `direct` | only the lattice points inside the $q$ window, 60 B per point | $O(N M)$ |
| `debye` | the pair histograms, independent of $N$ and of $q_{\max}$ | $O(N \times \text{neighbours})$ |

with $M$ the number of kept lattice points, which grows as $L^3$ in both grid
routes and therefore linearly with the atom count at a fixed $q_{\max}$.  The
NUFFT's internally upsampled fine grid is its largest single allocation:
FINUFFT upsamples by 1.25 while the requested tolerance is 1e-9 or looser
(about twice the mode grid) and by 2.0 for a tighter `--eps` (about eight
times), so asking for more than nine digits roughly doubles the memory - and
`cufinufft` implements only 2.0, so the GPU path stays there.  The direct sum
is the light one but pays $N$ times the work; the Debye histograms are the only
route that does not grow with the system size.

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

The same histograms give the two-body excess entropy (pair entropy) per
particle in units of $k_B$,

$$
S_2^{ab} = -2\pi\rho\,x_a x_b \int_0^{r_{\max}} r^2
\left[g_{ab}\ln g_{ab}-g_{ab}+1\right]dr,
$$

with the total following the partial sum rule
$S_2 = \sum_a S_2^{aa} + 2\sum_{a \lt b} S_2^{ab}$.  `--pair-entropy FILE` writes
the final total and partial values; `--s2-accum FILE` writes the $S_2(r)$
accumulation curve.  The integral uses the exact shell volume of each bin, the
$g\to0$ limit of the integrand and a leading-order Poisson bias correction;
no smoothing is applied by default and no statistical error is estimated.  A
single Debye run cannot make a rigorous tail correction (the $g(r)$ stops at
half the periodic box, and the Debye $S(q)$ has truncation ripples), so the
model-dependent part is left to `scripts/s2_analysis.py`, which can GCV-smooth
the partial $g(r)$, apply the `--dr` Richardson extrapolation and fit the tail
or multiple box sizes.  Do not use the Debye $S(q)$ for the reciprocal-space
entropy; use the NUFFT/direct $S(q)$ for that.

### Dynamic structure factor

`--dyn` keeps the time axis instead of averaging the snapshots away.  The
reciprocal space sampling is a separate choice, `--dyn-q`:

| `--dyn-q` | q points | outputs |
| --- | --- | --- |
| `-` (default) | none | `--chi4`: the overlap $Q(t)$, $\chi_4(t)$; `--msd` |
| `line:NINT,S0,S1,DX,DY,DZ` | `NINT+1` points on a line through $\Gamma$ | all |
| `shell:Q,ACC` | every direction of the shell $\|q\| = Q$ | all |
| `grid:QMAX` | every reciprocal-lattice vector with $\|q\| \le Q_{\max}$ | all |
| `single:N1,N2,N3` | one reciprocal-lattice vector of the box | all |

With `grid:QMAX` the q points are the reciprocal-lattice vectors of the box with
$|\mathbf q| \le Q_{\max}$, plus the $\Gamma$ point.  Only a lattice vector
carries a density amplitude that does not depend on the choice of the periodic
images, so this is the sampling an Ornstein-Zernike fit of $S_4(q,t)$ wants:
the tables hold one row per lattice shell at $q = (0,0,|q|)$, averaged over the
lattice vectors of that $|q|$ with their orbit multiplicities, the $\Gamma$ row
is $\chi_4(t)$, and a run with `--s4` also reports the spread of the symmetry
related modes of the smallest shell, which is the isotropy assumption behind
the reduction.  Add `--dyn-keep-modes` for the per-vector rows.

With `single:N1,N2,N3` the run samples exactly one reciprocal-lattice vector,
$\mathbf q = N_1\mathbf b_1 + N_2\mathbf b_2 + N_3\mathbf b_3$, with integer
indices that may be negative.  Building the vector from the integers keeps it
exactly on the lattice, which is what $S_4$ needs: a hand written $|q|$ would
sit a few ulps away and the box form factor would come back.  The single table
row is labelled by $|\mathbf q|$ and the header carries the indices and the
vector, so nothing is lost.

A run without q points computes no $S(q)$ at all: it takes no `OUTPUT` table,
and `--sqw`, `--fqt`, `--fqt-self` and `--s4` are rejected, which leaves the
cheapest route to $Q(t)$ and $\chi_4(t)$.

With a q line, `--dyn` follows the recipe used by the dynasor and MDANSE
packages: pick a line in reciprocal space, evaluate the density amplitudes
there, correlate them in time and Fourier transform the correlation.  With
$q_i = s_i\,\hat{u}$ and

$$
s_i = S_0 + i\,\frac{S_1 - S_0}{N_{\mathrm{int}}}
$$

(the line always passes through $\Gamma$, the origin of reciprocal space),

$$
F(q,t) = \frac{1}{W(q)} \sum_{ab} w_a(q)\, w_b(q)\,
         \frac{\langle \rho_a(q,t'+\tau)\,\rho_b^*(q,t')\rangle}
              {N_{\text{origins}}(\tau)}
$$

and

$$
S(q,\omega) = \frac{1}{2\pi}\int_{-\infty}^{\infty} e^{i\omega t} F(q,t)\,dt,
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

Because the scale steps by a constant, a line is also separable: writing
$q_i \cdot \mathbf r = S_0 (\hat u \cdot \mathbf r) + i\,\Delta s\,
(\hat u \cdot \mathbf r)$, the phase of mode $i$ is a mode independent start
times the $i$-th power of one factor, exactly as for a lattice mode.  A line of
16 points or more is therefore evaluated with one factor table per atom, and a
shorter one - where the table build would cost more than it saves - with a
sine and cosine per (mode, atom) pair.

With `--dyn-q shell:Q,ACC` the quantity is instead the isotropic average over
the sphere of radius $Q$,

$$
\bar X(Q) = \frac{\sum_k w_k X(q_k)}{\sum_k w_k},
$$

where the $q_k$ are the points of a Lebedev rule and the $w_k$ its weights
(which sum to $4\pi$).  `ACC` selects the rule: `low` is order 11 with 50
directions, `medium` order 17 with 110 and `high` order 23 with 194.  Every
output becomes a single row at $q = (0,0,Q)$, and the average is taken over the
normalized $S(q)$ of every direction, so that the $q$ dependent X-ray form
factors are averaged correctly even at a fixed $|q|$.

The rule is invariant under $q \to -q$ with the same weights, and every
quantity averaged on the shell is even in $q$ ($S_4$ and $F_s$ are squared
amplitudes, $F$ and $S$ enter through their real parts), so one vector of each
$\pm$ pair carries both with twice the weight.  A shell therefore accumulates
half of the rule's directions - 25, 55 or 97 - while its average is unchanged.

### Four-point structure factor and dynamic susceptibility

`--s4` and `--chi4` extend the same dynamic run to the four-point quantities
used to characterize dynamical heterogeneity.  They need the overlap cutoff
$a$ in the length unit of the dump,

$$
w_i(t,t_0) = \Theta\!\left(a - |\mathbf r_i(t_0+t) - \mathbf r_i(t_0)|\right),
$$

and the Fourier transform of the overlap field at the time origin,

$$
W(\mathbf q;t_0,t) = \sum_i e^{i\mathbf q\cdot\mathbf r_i(t_0)}\, w_i(t,t_0).
$$

The total four-point structure factor is the connected fluctuation of that
field,

$$
\begin{aligned}
S_4(\mathbf q,t) = \frac{1}{N}\Bigl[
  &\langle W(\mathbf q;t_0,t) W(-\mathbf q;t_0,t)\rangle_{t_0} \\
  &- \left|\langle W(\mathbf q;t_0,t)\rangle_{t_0}\right|^2
\Bigr],
\end{aligned}
$$

and `--chi4` writes the scalar $q = 0$ case,

$$
Q(t) = \frac{1}{N}\left\langle \sum_i w_i(t,t_0)\right\rangle_{t_0},
$$

$$
\chi_4(t) = \frac{1}{N}\left(\langle W_0(t)^2\rangle-\langle W_0(t)\rangle^2\right),
$$

with $W_0(t) = \sum_i w_i(t,t_0)$.  Equivalently,

$$
\chi_4(t) = N\left(\langle Q(t)^2\rangle-\langle Q(t)\rangle^2\right),
$$

so $Q(0) = 1$, $\chi_4(0) = 0$, and both decay at long times.  The $q \to 0$
limit of the four-point structure factor is the susceptibility,

$$
\lim_{q\to0} S_4(q,t) = \chi_4(t),
$$

and a small-$q$ Ornstein-Zernike fit,

$$
S_4(q,t) \simeq \frac{\chi_4(t)}{1+[q\,\xi(t)]^2},
$$

gives the dynamic correlation length $\xi(t)$.

The estimator is the unbiased connected covariance over time origins, with
the same origin count per lag as the coherent correlation.  `--s4-cutoff` is
required and has no default: the relevant $a$ is a fraction of the particle
diameter (or of the first-neighbour distance from $g(r)$), which the dump
does not carry.  The four-point quantities use unit overlap weights, so
`-w`/`--norm` continue to affect $S(q)$, $F(q,t)$ and $S(q,\omega)$ but not
$S_4$ or $\chi_4$.  They are total quantities: `--partials` does not add
partial four-point columns.

The $S_4$ table is written as one row per $(q,t)$ sample,

```
# qx qy qz tau S4(q,t)
```

and the chi4 table is a single time series,

```
# tau Q(t) chi4(t)
```

Both accept a `.h5`/`.hdf5` name or `--dyn-format hdf5`.  The HDF5 layout is
`/s4/q`, `/s4/tau`, `/s4/S4`, `/s4/count` and `/chi4/tau`, `/chi4/Q`,
`/chi4/chi4`, `/chi4/count`, with the overlap cutoff in the `overlap`
attribute.  The count is the number of time origins at each lag; a lag with a
single origin has no measurable fluctuation and is written as zero.

The calculation is more expensive than $F(q,t)$ because every atom and every
lag has to be compared with its origin position, and every atom inside the
cutoff contributes to every $q$ mode.  A large `--lag` thins the time origins,
and a short, low-$q$ line is usually enough for the Ornstein-Zernike fit.  The
position ring buffer grows with the effective window; `--buffer-limit GB`
sets its limit (default 2.0 GB, where GB is $10^9$ bytes).  The estimated size
is printed in the run summary, and an oversized request is rejected before any
large allocation.

`--stride N` subsamples the S4/chi4/F_s/MSD trajectory: only every N-th dump frame
contributes, the lag axis becomes $0, N, 2N, \ldots$ in dump frames, and the
position buffer stores only those frames.  This reduces both the buffer and
that work by roughly a factor $N$, at the cost of a coarser time axis.
The coherent $F(q,t)/S(q,\omega)$ outputs still use every dump frame and the
full `--maxframes`; only the S4/chi4/F_s/MSD window is reduced to the largest
multiple of $N$ that does not exceed `--maxframes`.  The run summary prints
the effective window and a note when it is shorter than the requested one.
`--lag` keeps its meaning in dump frames: an origin is used when its original
frame index is a multiple of `--lag`, so with stride $N$, `--lag m` gives an
origin every $m$ dump frames, not every $m$ S4 samples.

`--fqt-self` writes the self intermediate scattering function

$$
F_s(q,\tau) = \frac{1}{N}\left\langle \sum_i
  e^{i q\cdot[r_i(t+\tau)-r_i(t)]}\right\rangle,
$$

with unit weights.  The total is 1 at $\tau = 0$, and the per-species columns
$F_s(a)$ sum to it.  Like S4, it reuses the shared position buffer and the
`--stride`/`--lag` schedule, but it always includes every atom and never
applies `--s4-cutoff`, even when `--s4` is also requested.

### Mean squared displacement

`--msd` writes the mean squared displacement of the same trajectory,

$$
\text{MSD}(t) = \frac{1}{N}\left\langle \sum_i
  \left|\mathbf r_i(t_0+t) - \mathbf r_i(t_0)\right|^2\right\rangle_{t_0},
$$

with unit weights, so `-w`/`--norm` do not affect it and no overlap cutoff is
needed.  It is the self quantity that parallels $F_s$,

$$
\text{MSD}(t) = \lim_{q\to0} \frac{6\left[1 - F_s(q,t)\right]}{q^2},
$$

and its long-time slope gives the diffusion coefficient, $D =
\text{slope}/6$ in three dimensions.  Like S4 and $F_s$ it reuses the shared
position buffer and the `--stride`/`--lag` schedule, and it is available
without any q sampling (`--dyn-q -`), where a run writes no `OUTPUT` table.

The table carries one species column per LAMMPS type when `--partials` is on
(the default), with the same $1/N$ normalization as the $F_s$ columns, so the
species columns sum to the total and a partial divided by its concentration
tends to $6 D_a t$,

```
# tau MSD(t) MSD(Si) MSD(O)
      0.00000000    0.000000000000E+00    0.000000000000E+00    0.000000000000E+00
      1.00000000    5.961944468565E-01    2.996201670279E-01    2.965742798286E-01
```

`--no-partials` leaves only the total column.  A `.h5`/`.hdf5` name or
`--dyn-format hdf5` selects HDF5, with `/msd/tau`, `/msd/MSD`, `/msd/count`
and, with `--partials`, `/msd/pairs` and `/msd/MSD_partial/<species>`.

The dynamic method allocates only what the requested outputs need.  A run with
only `--s4`/`--chi4` does not allocate the coherent ring buffer, the
multi-origin correlation, or the $S(q,\omega)$/Fourier transform buffers; the
static `OUTPUT` table is still written.  Conversely, a run with only
`--sqw`/`--fqt` does not allocate the S4 position buffer, while `--fqt-self`
and `--msd` allocate that shared buffer but not the coherent
$F(q,t)/S(q,\omega)$ buffers.  With `--dyn-q -` there are no density
amplitudes at all, so `--chi4` and `--msd` are the only outputs and no
`OUTPUT` table is taken.

```sh
# low-q S4(q,t) and chi4(t), overlap cutoff a = 1.0 in dump length units
sqcalc -i traj.dump -w unit --dyn --dyn-q line:20,0.5,4,1,1,0 \
       --dt 0.005 --maxframes 400 \
       --s4-cutoff 1.0 --s4 S4.dat --chi4 chi4.dat S_q.dat
```

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
labelled `S(1-1)`, `S(1-2)`, ...; `-m` replaces those labels with the species
of the mapping and is required for the weighted schemes.  Report them with `-fz`:
every Faber-Ziman pair tends to 1, so the curves can be read and compared
across concentrations, which the OVITO pairs (tending to $x_a$) do not allow.
Every method accumulates the columns: `nufft` (CPU or CUDA), `direct` and
`debye`.  `direct` computes them in the same pass as the total, so the
reference implementation covers the partials as well.

A few limits are worth keeping in mind:

- the simulation box must not change along the trajectory (the reciprocal grid is built once from the first frame), both of nufft and Debye methods accept a non-periodic box (the reciprocal one still samples $q$ on the lattice of that box, so a cluster needs enough vacuum padding, whereas Debye needs none).

- The lattice fixes which $q$ are measurable: set `--qmin` to the smallest reciprocal vector and keep the shell width `(qmax - qmin)/nq` comparable to the lattice spacing, or the table comes back with zero-filled shells the box can never fill. The skill's `scripts/choose_q.py` prints a safe set of options for a given dump.

- Atoms of a type missing from `-m` under `-w neutron`/`xray`, or with an element that has no tabulated data, are rejected with a clear message.

### Powder XRD pattern

`--xrd FILE` adds a powder diffraction pattern in the convention of LAMMPS's
`compute xrd`:

$$
I(2\theta_b) = \frac{1}{N\,n_{\text{frames}}}
\sum_{q \in b} \lvert \rho(q) \rvert^{2} \, \mathrm{LP}(2\theta_q)
$$

The sum runs over the reciprocal lattice vectors $q$ of the box that fall into
the two-theta bin $b$ and $\rho(q)$ is the amplitude of the chosen weighting
(so the same form factors as the $S(q)$ table).  The Lorentz-polarization
factor is

$$
\mathrm{LP}(2\theta) = \frac{1 + \cos^2 2\theta}{\sin^2\theta \cos\theta},
$$

which `--no-lp` drops.  The wavelength fixes the $q$ window through

$$
q = \frac{4\pi \sin\theta}{\lambda},
$$

so `--xrd-lambda` and `--xrd-range` replace `--qmin` and `--qmax` for such a
run, and `--xrd` implies `-w xray` when no `-w` was given; `-w neutron` and
`-w unit` are accepted as well (the latter is a structure-only pattern).  The
intensity is per atom and per frame, so it is directly comparable with the
histogram `compute xrd` produces.  The self term is included, and the pattern
is the average over the trajectory of the intensity of each frame.

Bins are uniform in two-theta.  Without `--xrd-step` the width comes from the
box: the reciprocal lattice spacing $2\pi/L$ maps to
$\lambda/(L\cos\theta)$ in two-theta, which is the width a box of size $L$ can
give, so the default puts every line of the model into a single bin.  That
width diverges towards $2\theta = 180^\circ$, so the reference angle is capped
at $120^\circ$ of two-theta: a range that reaches beyond it is sampled more
coarsely than the box could, and the run counts the bins that then hold no
lattice point.  A finer
step is allowed and then shows the individual reciprocal lattice points of the
box (the run reports how many bins hold none), and a coarser step merges
neighbouring lines.  The lines are as sharp as the box permits: that is the
size broadening of the periodic model, not an instrument profile.

With `--partials` (the default) the pattern carries one column per type pair.
The pair intensity is the same cross term the $S(q)$ partials use, weighted by
the form factors of the pair:

$$
I_{ab}(2\theta_b) = \frac{1}{N\,n_{\text{frames}}}
\sum_{q \in b} f_a(q) f_b(q)\, \mathrm{Re}\left[\rho_a(q) \rho_b^*(q)\right]
\mathrm{LP}(2\theta_q),
$$

so the columns are a decomposition of the total:

$$
I(2\theta) = \sum_a I_{aa}(2\theta) + 2\sum_{a \lt b} I_{ab}(2\theta).
$$

They are labelled `I(Si-Si)`, `I(Si-O)`, ... from `-m`, or `I(1-1)`, ... when
no mapping was given; the self term sits in the diagonal ones and `-fz` does
not apply, because the table is in intensity units.  This is the simulation
side of what isotope substitution or anomalous scattering provides in an
experiment - in a two component pattern the pair columns say which sub-lattice
a feature comes from.  Every method accumulates them, `direct` included;
`--no-partials` drops them, as it does in the $S(q)$ table.

`--method debye` evaluates a different quantity - the orientation average
$\sum_{ij} f_i f_j \sin(Q r_{ij})/(Q r_{ij})$, which carries no multiplicity -
so `--xrd` refuses that combination instead of mixing the two conventions.

The table is text (`# 2theta[deg] I`, with the wavelength, range, step and the
LP switch in the comment lines above it) or HDF5 when the name ends in
`.h5`/`.hdf5` (`/xrd/two_theta`, `/xrd/I`, `/xrd/count`, with the settings as
file attributes).

## License

Copyright (C) 2026 Rui Su, Hangzhou Dianzi University

`sqcalc` is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the GNU General Public License (the `LICENSE` file, or
[https://www.gnu.org/licenses/](https://www.gnu.org/licenses/)) for more details.

### Third party code and dependencies

Version 3 is what covers the whole program, because the default build links
components under terms that do not otherwise combine:

* **FINUFFT** is vendored with `git subtree` in `external/finufft` (release
  v2.5.1, upstream commit `679d9ae5`).  It is Apache-2.0, Copyright (C)
  2017-2026 The Simons Foundation, Inc.; its `LICENSE` and `NOTICE` files are
  kept as upstream ships them.
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
  about debyer: [https://github.com/wojdyr/debyer](https://github.com/wojdyr/debyer).
  The X-ray table keeps the ion and valence rows of that source next to the
  neutral atoms (213 species), and the generator checks every row against the
  electron-count sum rule $f(0) = Z - q$ before writing it out.
