"""
Cut-and-project Fourier analysis: build a `HyperReciprocalLattice`,
enumerate Bragg peaks up to a cutoff, and hand the result back as a
`BraggPeakSet` that plugs into `LatticeCore.structure_factor` and
`LatticeCore.momentum_lattice`.

This file ships the generic peak-enumeration routine
(`bragg_peaks`) plus the topology-specific constructor for
[`FibonacciLattice`](@ref)'s `HyperReciprocalLattice`. Penrose and
Ammann–Beenker land in follow-up PRs, once their polygonal
acceptance windows are coded.
"""

# ---- hyper_reciprocal_lattice ---------------------------------------

"""
    hyper_reciprocal_lattice(qc::QuasicrystalData) → HyperReciprocalLattice

Construct the higher-dimensional reciprocal lattice description
for a cut-and-project quasicrystal `qc`. The returned object
carries

- `hyper_basis::SMatrix{DHyper, DHyper, T}` — basis of the hyper
  reciprocal lattice (`(Z^DHyper)*` = `2π · Z^DHyper` for the
  integer host lattices used by Fibonacci / Ammann–Beenker /
  Penrose)
- `parallel_proj::SMatrix{DPhys, DHyper, T}` — `π_∥`, the
  projection onto physical k-space
- `perp_proj::SMatrix{DPerp, DHyper, T}` — `π_⊥`, the projection
  onto internal (perpendicular) space
- `window::AcceptanceWindow` — the shape of the acceptance window
  in `E_⊥`; its Fourier transform governs Bragg peak intensities

Concrete lattice methods dispatch on the `topology` type parameter
of `qc`.
"""
function hyper_reciprocal_lattice end

# ---- Fibonacci: 2D host, 1D physical, 1D perp -----------------------

"""
    hyper_reciprocal_lattice(qc::QuasicrystalData{1, Float64, FibonacciLattice})

Build the Fibonacci chain's cut-and-project reciprocal structure.

Geometry:

- **Host lattice**: `Z²`, reciprocal `2π · Z²`, so
  `hyper_basis = 2π · I₂`.
- **Physical line `E_∥`**: direction `(1, 1/ϕ) / ‖·‖`. The parallel
  projection matrix is therefore

    `parallel_proj = (1, 1/ϕ) / √(1 + 1/ϕ²)` (1×2).

- **Perpendicular line `E_⊥`**: orthogonal direction `(-1/ϕ, 1) /
  ‖·‖`, so

    `perp_proj = (-1/ϕ, 1) / √(1 + 1/ϕ²)` (1×2).

- **Acceptance window**: a symmetric 1D interval. The Fibonacci
  generator's default choice corresponds to half-width `1/2` in
  the `(1, 1/ϕ)` frame, scaled into the orthogonal `(e_∥, e_⊥)`
  frame by the same `1 / √(1 + 1/ϕ²)` factor. We expose the
  projected half-width directly.

`ϕ = (1 + √5)/2` is the golden ratio (`LatticeCore.ϕ` is also
re-exported through QuasiCrystal).
"""
function hyper_reciprocal_lattice(
    ::QuasicrystalData{1,Float64,FibonacciLattice}
)
    T = Float64

    # Host lattice is Z²; its reciprocal is 2π·Z².
    hyper_basis = SMatrix{2,2,T}(2π, 0.0, 0.0, 2π)

    # Unit vectors for the physical and perpendicular 1D subspaces.
    norm_par = sqrt(1 + inv(ϕ)^2)
    parallel_proj = SMatrix{1,2,T}(1 / norm_par, inv(ϕ) / norm_par)
    perp_proj = SMatrix{1,2,T}(-inv(ϕ) / norm_par, 1 / norm_par)

    # Acceptance window in the perpendicular subspace. The default
    # projection builder (`generate_fibonacci_projection`) accepts
    # points inside a unit-width strip in the (1, 1/ϕ) frame, which
    # becomes a 1 / √(1+1/ϕ²) half-width = 1 / (2·norm_par) on
    # E_⊥ after the orthonormalisation above.
    window = IntervalWindow(T(1 / (2 * norm_par)))

    return HyperReciprocalLattice{1,2,1,T,typeof(window)}(
        hyper_basis, parallel_proj, perp_proj, window
    )
end

# ---- Bragg peak enumeration -----------------------------------------

"""
    bragg_peaks(qc::QuasicrystalData;
                kmax::Real,
                intensity_cutoff::Real = 0.0) → BraggPeakSet

Enumerate the finite set of Bragg peaks of a cut-and-project
quasicrystal `qc` whose projected physical momentum has
`‖k‖ ≤ kmax` and whose intensity `|Ŵ(g_⊥)|²` exceeds
`intensity_cutoff`. Returns a `BraggPeakSet{DPhys, DHyper, T}`
that inherits from `LatticeCore.AbstractMomentumLattice` and is
therefore directly usable as a momentum lattice for
`LatticeCore.structure_factor` and the MC observer stack.

The enumeration walks every integer point in the hyper reciprocal
lattice whose parallel projection is inside the `kmax` ball. The
bound `n_max` on the integer coordinates is chosen conservatively
from the operator norm of `parallel_proj` so that nothing bright
slips through unnoticed.
"""
function bragg_peaks(
    qc::QuasicrystalData; kmax::Real, intensity_cutoff::Real=0.0
)
    return _bragg_peaks_impl(hyper_reciprocal_lattice(qc), kmax, intensity_cutoff)
end

function _bragg_peaks_impl(
    hrl::HyperReciprocalLattice{DPhys,DHyper,DPerp,T,W},
    kmax::Real,
    intensity_cutoff::Real,
) where {DPhys,DHyper,DPerp,T,W}
    par = hrl.parallel_proj
    perp = hrl.perp_proj
    A = hrl.hyper_basis

    # Bound the enumeration box. The worst-case contribution to
    # ‖k_∥‖ from a hyper index `n` with `|nᵢ| ≤ n_max` is
    # `opnorm(par · A) * n_max * √DHyper`. Invert for n_max.
    par_A = par * A
    row_norms = [sqrt(sum(abs2, par_A[i, :])) for i in 1:DPhys]
    bound = maximum(row_norms) * sqrt(DHyper)
    n_max = bound > 0 ? ceil(Int, T(kmax) / bound) + 1 : 1

    peaks = SVector{DPhys,T}[]
    intensities = T[]
    hyper_indices = NTuple{DHyper,Int}[]

    for idx in Iterators.product(ntuple(_ -> (-n_max):n_max, DHyper)...)
        idx_int = SVector{DHyper,Int}(idx)
        # Zero hyper index ⇒ trivial Γ peak; keep it so users see a
        # non-empty result for "ferromagnet on a quasicrystal".
        g_hyper = A * SVector{DHyper,T}(T.(Tuple(idx_int)))
        k_par = par * g_hyper
        # StaticArrays returns SVector for matrix-vector multiply.
        k_perp = perp * g_hyper

        # Radius filter in physical k-space.
        kpar_norm = sqrt(sum(abs2, k_par))
        kpar_norm > T(kmax) && continue

        # Window Fourier transform at the perpendicular momentum.
        amp = window_fourier(hrl.window, k_perp)
        I = T(abs2(amp))

        I < T(intensity_cutoff) && continue

        push!(peaks, SVector{DPhys,T}(k_par))
        push!(intensities, I)
        push!(hyper_indices, Tuple(idx_int))
    end

    return BraggPeakSet{DPhys,DHyper,T}(peaks, intensities, hyper_indices)
end

# ---- fourier_module dispatch ----------------------------------------

"""
    LatticeCore.fourier_module(qc::QuasicrystalData; kmax, intensity_cutoff) → BraggPeakSet

QuasiCrystal-side override of the generic
[`LatticeCore.fourier_module`](@ref LatticeCore.fourier_module)
entry point. Calls [`bragg_peaks`](@ref) with the requested cutoff
and returns a `BraggPeakSet`.

Currently implemented for `FibonacciLattice` only. Passing a
Penrose / Ammann–Beenker lattice raises a `MethodError` on
`hyper_reciprocal_lattice`.
"""
function LatticeCore.fourier_module(
    qc::QuasicrystalData; kmax::Real=20.0, intensity_cutoff::Real=0.0
)
    return bragg_peaks(qc; kmax=kmax, intensity_cutoff=intensity_cutoff)
end
