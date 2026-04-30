module QuasiCrystalNFFTExt

using QuasiCrystal
using LatticeCore
using NFFT
using StaticArrays
using LinearAlgebra

"""
    QuasiCrystalNFFTExt

Optional extension that installs an NUFFT-backed fast path for
`structure_factor(qc, state, ml)` when `qc::QuasicrystalData` carries
the `HasFourierModule` trait.

Loaded automatically once the user issues `using NFFT`. The
implementation routes through `NFFT.NNDFTPlan` (non-uniform → non-uniform
direct DFT) so that both site positions and target k-points are treated
as arbitrary point clouds, which matches the structure-factor problem
on a quasicrystal: site positions live on a quasiperiodic point set and
the target k-grid is typically a `BraggPeakSet`.

`NNDFTPlan` is mathematically equivalent to the naive `O(N · M)` double
loop — its job is to provide a structured, library-backed entry point
that other extensions / tooling can specialise. Numerical results match
the naive path to floating-point round-off; the test-suite asserts this
explicitly. A genuine `O((N + M) log)` Type-3 NUFFT (gridding + FFT +
interpolation) is the natural next step but is not yet implemented in
NFFT.jl's public API; this extension is the dispatch site for that
follow-up.
"""
QuasiCrystalNFFTExt

# ---- Trait override --------------------------------------------------

# Override the LatticeCore stub installed by `LatticeCoreNFFTExt`.
# Both extensions specialise on `HasFourierModule`; the QC override is
# strictly more specific (lat type) and so wins dispatch whenever
# QuasiCrystal + NFFT are both loaded.
function LatticeCore._structure_factor_fast(
    ::LatticeCore.HasFourierModule,
    lat::QuasiCrystal.QuasicrystalData,
    state::AbstractVector,
    ml::LatticeCore.AbstractMomentumLattice,
)
    return _structure_factor_nufft(lat, state, ml)
end

# ---- NUFFT path ------------------------------------------------------

# Compute S(k_j) = |Σ_n state_n exp(-i k_j · r_n)|² / N for every k-point
# of `ml`, using NFFT.jl's `NNDFTPlan`.
#
# `NNDFTPlan(k, y)` evaluates
#
#     g[j] = Σ_l f[l] · exp(-2π i Σ_d k[d, j] · y[d, l])
#
# i.e. its "k" matrix carries an implicit factor of 2π relative to the
# physical k-vectors used in `structure_factor`. We absorb this by
# passing `k / (2π)`.
function _structure_factor_nufft(
    lat::QuasiCrystal.QuasicrystalData{D,T},
    state::AbstractVector,
    ml::LatticeCore.AbstractMomentumLattice,
) where {D,T}
    N = LatticeCore.num_sites(lat)
    M = LatticeCore.num_k_points(ml)

    # Stack site positions into a D × N matrix. We promote to Float64
    # to match NFFT's preferred working type and to avoid mixed-type
    # plan construction.
    R = _positions_matrix(lat, N)

    # Stack target k-points and rescale by 1/(2π).
    K = _k_matrix(ml, M)
    K ./= 2π

    plan = NFFT.NNDFTPlan(K, R)
    s_complex = ComplexF64.(state)
    g = Vector{ComplexF64}(undef, M)
    mul!(g, plan, s_complex)

    out = Vector{Float64}(undef, M)
    invN = 1.0 / N
    @inbounds for j in 1:M
        out[j] = abs2(g[j]) * invN
    end
    return out
end

@inline function _positions_matrix(lat::LatticeCore.AbstractLattice, N::Int)
    D = _spatial_dim(lat)
    R = Matrix{Float64}(undef, D, N)
    @inbounds for n in 1:N
        p = LatticeCore.position(lat, n)
        for d in 1:D
            R[d, n] = Float64(p[d])
        end
    end
    return R
end

@inline function _k_matrix(ml::LatticeCore.AbstractMomentumLattice, M::Int)
    D = _spatial_dim(ml)
    K = Matrix{Float64}(undef, D, M)
    @inbounds for j in 1:M
        k = LatticeCore.k_point(ml, j)
        for d in 1:D
            K[d, j] = Float64(k[d])
        end
    end
    return K
end

# Helper: extract D from an AbstractLattice{D, T} subtype.
_spatial_dim(::LatticeCore.AbstractLattice{D}) where {D} = D

end # module QuasiCrystalNFFTExt
