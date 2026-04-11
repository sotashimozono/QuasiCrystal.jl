# Cut-and-project construction

This page is a short mathematical primer on cut-and-project
quasicrystals — enough background to read the rest of the docs
and judge whether the numerical results make sense.

## Host lattice, parallel, and perpendicular subspaces

A cut-and-project quasicrystal starts from a *host* lattice
$\mathbb{Z}^{D_\mathrm{hyper}}$ embedded in
$\mathbb{R}^{D_\mathrm{hyper}}$. Write the ambient space as a
direct sum

```math
\mathbb{R}^{D_\mathrm{hyper}} \;=\; E_\parallel \oplus E_\perp,
```

where $E_\parallel$ has dimension $D_\mathrm{phys}$ (the "physical"
space the quasicrystal lives in) and $E_\perp$ has dimension
$D_\mathrm{perp} = D_\mathrm{hyper} - D_\mathrm{phys}$. Let
$\pi_\parallel$ and $\pi_\perp$ denote the projections onto those
subspaces.

For each of the quasicrystals this package ships:

| Model            | $D_\mathrm{hyper}$ | $D_\mathrm{phys}$ | $D_\mathrm{perp}$ |
|------------------|:------------------:|:-----------------:|:-----------------:|
| Fibonacci chain  | 2                  | 1                 | 1                 |
| Ammann–Beenker   | 4                  | 2                 | 2                 |
| Penrose P3       | 5                  | 2                 | 3                 |

The two projection matrices live on the returned
`HyperReciprocalLattice` object as `parallel_proj`
(`DPhys × DHyper`) and `perp_proj` (`DPerp × DHyper`).

## Acceptance window and the point set

An acceptance window $W \subset E_\perp$ selects which integer
host points survive the projection:

```math
\Lambda \;=\; \bigl\{\,\pi_\parallel(n) \;:\; n \in \mathbb{Z}^{D_\mathrm{hyper}},\;
\pi_\perp(n) \in W\,\bigr\}.
```

The window shape is a property of the model:

- Fibonacci: `IntervalWindow` (a symmetric 1D interval).
- Ammann–Beenker: `BoxWindow{2}` (a square).
- Penrose P3: `BoxWindow{3}` (a cube in the 3D internal space).

The point set $\Lambda$ is the output of
`generate_fibonacci_projection`,
`generate_ammann_beenker_projection`, and
`generate_penrose_projection`.

## Bragg peaks from the window's Fourier transform

The *key* identity of the cut-and-project construction is that
the diffraction pattern of $\Lambda$ is a pure point measure
supported on the projection of the host reciprocal lattice:

```math
\hat{\rho}(k) \;\propto\; \sum_{g \in (\mathbb{Z}^{D_\mathrm{hyper}})^{*}}
\hat{W}\!\bigl(\pi_\perp(g)\bigr)\, \delta\!\bigl(k - \pi_\parallel(g)\bigr).
```

In words: **every Bragg peak sits at a projected reciprocal
lattice vector, and its amplitude is the Fourier transform of the
acceptance window evaluated at the *perpendicular* component of
that reciprocal vector.** Bright peaks correspond to reciprocal
vectors whose perpendicular projection is small (so $\hat{W}$ is
near its maximum); dim peaks correspond to large perpendicular
components (where $\hat{W}$ has decayed).

Intensities are the squared magnitude:

```math
I(g) \;\propto\; \bigl|\hat{W}\!\bigl(\pi_\perp(g)\bigr)\bigr|^2.
```

This is what the generic [`bragg_peaks`](@ref) enumeration
computes: walk every nonzero integer point of the host reciprocal
lattice inside a cutoff ball, evaluate $\hat{W}$ at the
perpendicular projection, and keep the ones whose intensity
exceeds `intensity_cutoff`.

## Practical notes on the shipped models

- The current Fibonacci projection uses the canonical physical
  direction $(1, 1/\varphi)$ with $\varphi = (1 + \sqrt{5})/2$ and
  its orthogonal complement as the perpendicular line.
- The Ammann–Beenker and Penrose P3 generators use a non-canonical
  shifted perpendicular star (matching the generator's current
  window definition). The resulting peak set is still closed
  under the physical-space rotational symmetry (C₄ and C₅
  respectively), but the "sharpest" canonical variant — which
  needs an octagonal or pentagonal window in an orthogonal
  Galois-conjugate perpendicular plane — is still future work.
  See the note at the end of `src/core/fourier/fourier.jl`.
