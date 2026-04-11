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
function hyper_reciprocal_lattice(::QuasicrystalData{1,Float64,FibonacciLattice})
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

# ---- Ammann–Beenker: 4D host, 2D physical, 2D perp ------------------

"""
    hyper_reciprocal_lattice(qc::QuasicrystalData{2, Float64, AmmannBeenker})

Build the Ammann–Beenker cut-and-project reciprocal structure.

Geometry (matching [`generate_ammann_beenker_projection`](@ref)):

- **Host lattice**: `Z⁴`, so `hyper_basis = 2π · I₄`.
- **Physical plane `E_∥`**: the four host basis vectors `eₖ`
  (k = 1..4) project to the star vectors
  `(cos((k-1)π/4), sin((k-1)π/4))`, giving a 2×4
  `parallel_proj` whose columns are those star vectors.
- **Perpendicular plane `E_⊥`**: shifted star vectors with angles
  `(k-1)π/4 + π/4`. (Matches the generator's current construction
  — note this is not the canonical Galois-conjugate embedding, but
  it is the one that is consistent with the points actually
  produced by the generator.)
- **Acceptance window**: a square `|y₁|, |y₂| ≤ ½`, modelled as a
  [`BoxWindow{2}`](@ref).

The resulting `BraggPeakSet` plugs into
`LatticeCore.momentum_lattice` and composes with
`structure_factor` just like the Fibonacci version.
"""
function hyper_reciprocal_lattice(::QuasicrystalData{2,Float64,AmmannBeenker})
    T = Float64
    hyper_basis = SMatrix{4,4,T}(2π, 0, 0, 0, 0, 2π, 0, 0, 0, 0, 2π, 0, 0, 0, 0, 2π)

    θ = T(π / 4)
    c0, s0 = cos(zero(T)), sin(zero(T))
    c1, s1 = cos(θ), sin(θ)
    c2, s2 = cos(2θ), sin(2θ)
    c3, s3 = cos(3θ), sin(3θ)
    c4, s4 = cos(4θ), sin(4θ)

    # Physical star (k = 1..4): angles 0, π/4, π/2, 3π/4.
    # SMatrix is column-major, so each (col_1, col_2) pair is one
    # star vector.
    parallel_proj = SMatrix{2,4,T}(c0, s0, c1, s1, c2, s2, c3, s3)

    # Perpendicular star: angles π/4, π/2, 3π/4, π — the shifted
    # choice the generator currently filters with.
    perp_proj = SMatrix{2,4,T}(c1, s1, c2, s2, c3, s3, c4, s4)

    # Square acceptance window of half-width 0.5 on each perp axis,
    # matching `generate_ammann_beenker_projection`'s
    # `all(abs.(pos_perp) .<= 0.5)` filter.
    window = BoxWindow(SVector{2,T}(0.5, 0.5))

    return HyperReciprocalLattice{2,4,2,T,typeof(window)}(
        hyper_basis, parallel_proj, perp_proj, window
    )
end

# ---- Penrose P3: 5D host, 2D physical, 3D perp ----------------------

"""
    hyper_reciprocal_lattice(qc::QuasicrystalData{2, Float64, PenroseP3})

Build the Penrose P3 cut-and-project reciprocal structure.

Geometry (matching [`generate_penrose_projection`](@ref)):

- **Host lattice**: `Z⁵`, so `hyper_basis = 2π · I₅`.
- **Physical plane `E_∥`**: the five host basis vectors `eₖ`
  (k = 1..5) project to the 5-fold star
  `(cos((k-1)·2π/5), sin((k-1)·2π/5))`. `parallel_proj` is the
  2×5 matrix whose columns are those star vectors.
- **Perpendicular space `E_⊥ ⊂ R³`**: coordinates are
  `(cos(2(k-1)·2π/5), sin(2(k-1)·2π/5), cos(3(k-1)·2π/5))`,
  matching the generator's current (non-canonical) choice. The
  canonical 5-fold Penrose uses a 2D Galois-conjugate
  perpendicular plane with a pentagonal acceptance window — see
  the follow-up note at the bottom of this file.
- **Acceptance window**: a cube `|yᵢ| ≤ ½` in `R³`, modelled as a
  [`BoxWindow{3}`](@ref).
"""
function hyper_reciprocal_lattice(::QuasicrystalData{2,Float64,PenroseP3})
    T = Float64
    hyper_basis = SMatrix{5,5,T}(
        2π, 0, 0, 0, 0,
        0, 2π, 0, 0, 0,
        0, 0, 2π, 0, 0,
        0, 0, 0, 2π, 0,
        0, 0, 0, 0, 2π,
    )

    θ = T(2π / 5)

    # parallel_proj columns: (cos((k-1)θ), sin((k-1)θ)), k=1..5
    parallel_proj = SMatrix{2,5,T}(
        cos(0 * θ), sin(0 * θ),
        cos(1 * θ), sin(1 * θ),
        cos(2 * θ), sin(2 * θ),
        cos(3 * θ), sin(3 * θ),
        cos(4 * θ), sin(4 * θ),
    )

    # perp_proj columns: (cos(2(k-1)θ), sin(2(k-1)θ), cos(3(k-1)θ))
    perp_proj = SMatrix{3,5,T}(
        cos(2 * 0 * θ), sin(2 * 0 * θ), cos(3 * 0 * θ),
        cos(2 * 1 * θ), sin(2 * 1 * θ), cos(3 * 1 * θ),
        cos(2 * 2 * θ), sin(2 * 2 * θ), cos(3 * 2 * θ),
        cos(2 * 3 * θ), sin(2 * 3 * θ), cos(3 * 3 * θ),
        cos(2 * 4 * θ), sin(2 * 4 * θ), cos(3 * 4 * θ),
    )

    window = BoxWindow(SVector{3,T}(0.5, 0.5, 0.5))

    return HyperReciprocalLattice{2,5,3,T,typeof(window)}(
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
function bragg_peaks(qc::QuasicrystalData; kmax::Real, intensity_cutoff::Real=0.0)
    return _bragg_peaks_impl(hyper_reciprocal_lattice(qc), kmax, intensity_cutoff)
end

function _bragg_peaks_impl(
    hrl::HyperReciprocalLattice{DPhys,DHyper,DPerp,T,W}, kmax::Real, intensity_cutoff::Real
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
