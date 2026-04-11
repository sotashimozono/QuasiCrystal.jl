# Fourier analysis pipeline

The Fourier pipeline is a narrow, composable stack of four pieces
that together deliver a concrete `BraggPeakSet` from a
`QuasicrystalData`:

```text
generate_*_projection(...)          # QuasicrystalData
        │
        ▼
hyper_reciprocal_lattice(qc)        # HyperReciprocalLattice
        │
        ▼
bragg_peaks(qc; kmax, intensity_cutoff)  # BraggPeakSet
        │
        ▼
momentum_lattice(qc) / structure_factor(qc, state, peaks)
```

All three shipped quasicrystals route through the same generic
enumeration — the only per-model piece is
`hyper_reciprocal_lattice`.

## Step 1 — `HyperReciprocalLattice`

[`hyper_reciprocal_lattice`](@ref) constructs the
higher-dimensional object that describes the quasicrystal's
Fourier structure. It holds four fields:

| Field            | Type                           | Meaning                          |
|------------------|--------------------------------|----------------------------------|
| `hyper_basis`    | `SMatrix{DHyper, DHyper, T}`   | Host reciprocal basis            |
| `parallel_proj`  | `SMatrix{DPhys, DHyper, T}`    | $\pi_\parallel$                  |
| `perp_proj`      | `SMatrix{DPerp, DHyper, T}`    | $\pi_\perp$                      |
| `window`         | `<:AcceptanceWindow`           | Shape of the window in $E_\perp$ |

For the three shipped models the host lattice is
$\mathbb{Z}^{D_\mathrm{hyper}}$, so `hyper_basis = 2π · I`.

## Step 2 — Window Fourier transforms

An acceptance window's brightness law
$I(g) \propto |\hat{W}(\pi_\perp g)|^2$ needs $\hat{W}$ to be
cheap to evaluate at arbitrary perpendicular momenta. The two
windows shipped in QuasiCrystal both have closed-form transforms:

- [`IntervalWindow`](@ref) (1D) has
  $\hat{W}(q) = 2 \sin(qa)/q$ where $a$ is the half-width.
- [`BoxWindow`](@ref) (axis-aligned hyper-rectangle) has a
  separable product of 1D sincs:
  $\hat{W}(q) = \prod_i 2 \sin(q_i a_i) / q_i$.

Both implementations handle the $q \to 0$ limit exactly
(returning $2a$ or $\prod 2a_i$) to avoid a $0/0$ at Γ.

## Step 3 — Bragg peak enumeration

[`bragg_peaks`](@ref) walks every integer hyper index inside a
conservative bounding box, evaluates the peak location
$k_\parallel = \pi_\parallel (2\pi n)$ and window amplitude
$\hat{W}(\pi_\perp (2\pi n))$, and retains peaks with
$\|k_\parallel\| \le k_\mathrm{max}$ and intensity above
`intensity_cutoff`.

The result is a `BraggPeakSet{DPhys, DHyper, T}` from
LatticeCore, which is a subtype of `AbstractMomentumLattice` and
therefore usable as a drop-in momentum lattice.

## Step 4 — Trait-based dispatch

Because `QuasicrystalData` reports
`reciprocal_support(qc) = HasFourierModule()`, the generic
`LatticeCore.momentum_lattice(qc)` helper routes through
`fourier_module(qc)`, which in turn calls `bragg_peaks` with a
default cutoff. The upshot is that any code written against
`AbstractMomentumLattice` — including
`LatticeCore.structure_factor` — works on a quasicrystal exactly
the same way it works on a periodic lattice. No branches,
no specialised code paths.
