"""
    QuasiCrystal single-particle spectral function A(k, ω)
    ======================================================

The dynamic / momentum-resolved spectral function of the tight-binding model
on a quasicrystal,

```math
A(k, ω) = \\sum_n |\\langle k | ψ_n \\rangle|^2 \\, δ(ω - E_n),
\\qquad
\\langle k | ψ_n \\rangle = \\frac{1}{\\sqrt{N}} \\sum_l e^{-i k · r_l} ψ_n(l),
```

where `(E_n, ψ_n)` are the [`eigenstates`](@ref) of the tight-binding
Hamiltonian and `|k⟩` is the plane wave sampled on the (aperiodic) site set.
On a periodic crystal this collapses onto the Bloch bands `E(k)`; on a
quasicrystal it resolves the fragmented "band structure" and the pseudo-gaps
of the singular-continuous spectrum (issue #27).

The δ-functions are replaced by a normalized Gaussian of width `broadening`,
so the spectral weight is conserved:

```math
\\int A(k, ω)\\, dω = \\sum_n |\\langle k | ψ_n \\rangle|^2 = \\langle k | k \\rangle = 1 .
```
"""

# ---- spectral function ----------------------------------------------

"""
    dynamic_structure_factor(qc::QuasicrystalData,
                             kpoints::AbstractVector{<:AbstractVector};
                             t=1.0, onsite=0.0,
                             omega=nothing, nomega::Int=200, broadening::Real=0.05)
        -> (omegas::Vector{Float64}, A::Matrix{Float64})

Momentum-resolved single-particle spectral function `A(k, ω)` of the
tight-binding model on `qc`, evaluated at every k-vector in `kpoints`.
Returns the shared frequency grid `omegas` (length `nomega`) and a
`length(kpoints) × nomega` matrix `A`, with `A[j, m] = A(kpoints[j], omegas[m])`.

Each δ-peak is broadened by a normalized Gaussian of width `broadening`, so
`sum(A[j, :]) * Δω ≈ 1` for every `j` (spectral-weight sum rule). Pass an
explicit `omega` grid (a range or vector) to override the auto-range, which
otherwise spans the spectrum padded by `3·broadening`.

See also [`spectrum`](@ref), [`structure_factor`](@ref).
"""
function dynamic_structure_factor(
    qc::QuasicrystalData,
    kpoints::AbstractVector{<:AbstractVector};
    t::Real=1.0,
    onsite=0.0,
    omega=nothing,
    nomega::Int=200,
    broadening::Real=0.05,
)
    nomega ≥ 1 || throw(ArgumentError("nomega must be ≥ 1, got $nomega"))
    broadening > 0 || throw(ArgumentError("broadening must be > 0, got $broadening"))

    E, V = eigenstates(qc; t=t, onsite=onsite)
    N = length(E)
    R = _positions_matrix(qc, N)                     # D × N (reused from FINUFFT ext helper? no)

    omegas = _omega_grid(omega, E, Float64(broadening), nomega)
    σ = Float64(broadening)
    gnorm = 1 / (σ * sqrt(2π))

    nk = length(kpoints)
    A = zeros(Float64, nk, length(omegas))

    # weight w_n(k) = |⟨k|ψ_n⟩|²  with ⟨k|ψ_n⟩ = (1/√N) Σ_l e^{-i k·r_l} ψ_n(l)
    phase = Vector{ComplexF64}(undef, N)
    @inbounds for j in 1:nk
        k = kpoints[j]
        for l in 1:N
            kr = 0.0
            for d in eachindex(k)
                kr += k[d] * R[d, l]
            end
            phase[l] = cis(-kr)
        end
        for n in 1:N
            c = zero(ComplexF64)
            for l in 1:N
                c += phase[l] * V[l, n]
            end
            w = abs2(c) / N
            En = E[n]
            for m in eachindex(omegas)
                A[j, m] += w * gnorm * exp(-((omegas[m] - En)^2) / (2σ^2))
            end
        end
    end
    return omegas, A
end

"""
    dynamic_structure_factor(qc::QuasicrystalData, ml::AbstractMomentumLattice; kwargs...)

Convenience method taking the k-points from a momentum lattice (e.g. a
`BraggPeakSet`), so the spectral function can be evaluated on the same
k-grid used for the static [`structure_factor`](@ref).
"""
function dynamic_structure_factor(
    qc::QuasicrystalData, ml::LatticeCore.AbstractMomentumLattice; kwargs...
)
    kpoints = [k_point(ml, j) for j in 1:num_k_points(ml)]
    return dynamic_structure_factor(qc, kpoints; kwargs...)
end

# ---- helpers ---------------------------------------------------------

function _omega_grid(omega, E::AbstractVector, σ::Float64, nomega::Int)
    omega === nothing || return collect(Float64.(omega))
    Emin, Emax = extrema(E)
    if Emax ≈ Emin
        Emin -= 0.5
        Emax += 0.5
    end
    return collect(range(Emin - 3σ, Emax + 3σ; length=nomega))
end

# Local D × N Float64 position matrix (self-contained; the FINUFFT extension
# has its own copy that is only loaded with FINUFFT).
function _positions_matrix(qc::QuasicrystalData, N::Int)
    D = _spatial_dim(qc)
    R = Matrix{Float64}(undef, D, N)
    @inbounds for l in 1:N
        p = position(qc, l)
        for d in 1:D
            R[d, l] = Float64(p[d])
        end
    end
    return R
end

_spatial_dim(::LatticeCore.AbstractLattice{D}) where {D} = D
