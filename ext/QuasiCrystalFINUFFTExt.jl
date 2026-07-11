module QuasiCrystalFINUFFTExt

using QuasiCrystal
using LatticeCore
using FINUFFT
using LinearAlgebra

"""
    QuasiCrystalFINUFFTExt

Optional extension that installs a genuine `O((N + M) log(1/ε))` NUFFT fast
path for `structure_factor(qc, state, ml)` when `qc::QuasicrystalData`
carries the `HasFourierModule` trait.

Loaded automatically once the user issues `using FINUFFT`. The structure
factor of a quasicrystal is a **type-3** (non-uniform → non-uniform)
transform: the `N` site positions form a quasiperiodic point cloud and the
`M` target k-points are typically an irregular `BraggPeakSet`. FINUFFT.jl
(the Julia binding of the Flatiron Institute FINUFFT library) provides
`nufft1d3` / `nufft2d3` / `nufft3d3` for exactly this problem, so this path
is asymptotically faster than the naive `O(N · M)` double loop instead of
merely a structured restatement of it (cf. issue #61).

Numerical results match the naive reference to the requested tolerance
(`ε = 1e-12` here, giving agreement at the `1e-10` level, asserted by the
test-suite). At `N = M = 10^4` the FINUFFT path is ~25× faster than the
naive loop on a laptop CPU. Small problems (work `N · M` below
`_NUFFT_MIN_WORK`) are routed back to the naive loop, whose lack of setup
overhead makes it faster there — so enabling this extension never regresses
performance.
"""
QuasiCrystalFINUFFTExt

# Tolerance for the FINUFFT transform. 1e-12 keeps the structure-factor
# agreement with the naive path comfortably below the 1e-10 test bound.
const _FINUFFT_EPS = 1e-12

# Crossover heuristic. FINUFFT type-3 costs ~O((N + M) log(1/ε)) *plus* a
# fixed plan/spread overhead of tens of ms, whereas the naive double loop is
# O(N · M) with essentially zero setup. Below this work threshold the naive
# loop is faster (measured crossover ≈ 1e6 on a laptop CPU), so we route
# small problems back to the naive path to guarantee the fast path is never
# a regression. Above it FINUFFT wins decisively (≈25× at N = M = 1e4).
const _NUFFT_MIN_WORK = 2_000_000

# ---- Trait override --------------------------------------------------

# Override the LatticeCore stub installed by the naive default. This method
# is strictly more specific (QuasicrystalData) and so wins dispatch whenever
# QuasiCrystal + FINUFFT are both loaded.
function LatticeCore._structure_factor_fast(
    ::LatticeCore.HasFourierModule,
    lat::QuasiCrystal.QuasicrystalData,
    state::AbstractVector,
    ml::LatticeCore.AbstractMomentumLattice,
)
    N = LatticeCore.num_sites(lat)
    M = LatticeCore.num_k_points(ml)
    if N * M < _NUFFT_MIN_WORK
        return LatticeCore._structure_factor_naive(lat, state, ml)
    end
    return _structure_factor_finufft(lat, state, ml)
end

# ---- NUFFT path ------------------------------------------------------

# Compute S(k_j) = |Σ_n state_n exp(-i k_j · r_n)|² / N for every k-point
# of `ml`, using a FINUFFT type-3 transform.
#
# FINUFFT's `nufftNd3(x..., c, iflag, eps, s...)` evaluates
#
#     g[j] = Σ_n c[n] · exp(iflag · i · (s[j]·x[n] + ...))
#
# using the *raw* (no implicit 2π) exponent, so `iflag = -1` with the
# physical positions and k-points reproduces the structure-factor sum
# directly — no rescaling is needed.
function _structure_factor_finufft(
    lat::QuasiCrystal.QuasicrystalData{D,T},
    state::AbstractVector,
    ml::LatticeCore.AbstractMomentumLattice,
) where {D,T}
    N = LatticeCore.num_sites(lat)
    M = LatticeCore.num_k_points(ml)

    R = _positions_matrix(lat, N)      # D × N
    K = _k_matrix(ml, M)               # D × M
    c = ComplexF64.(state)

    g = _nufft_type3(Val(D), R, K, c)

    out = Vector{Float64}(undef, M)
    invN = 1.0 / N
    @inbounds for j in 1:M
        out[j] = abs2(g[j]) * invN
    end
    return out
end

# Dimension-specialised type-3 dispatch. `R` is D × N, `K` is D × M.
function _nufft_type3(::Val{1}, R::Matrix{Float64}, K::Matrix{Float64}, c)
    x = @view R[1, :]
    s = @view K[1, :]
    return nufft1d3(collect(x), c, -1, _FINUFFT_EPS, collect(s))
end

function _nufft_type3(::Val{2}, R::Matrix{Float64}, K::Matrix{Float64}, c)
    x = collect(@view R[1, :])
    y = collect(@view R[2, :])
    s = collect(@view K[1, :])
    t = collect(@view K[2, :])
    return nufft2d3(x, y, c, -1, _FINUFFT_EPS, s, t)
end

function _nufft_type3(::Val{3}, R::Matrix{Float64}, K::Matrix{Float64}, c)
    x = collect(@view R[1, :])
    y = collect(@view R[2, :])
    z = collect(@view R[3, :])
    s = collect(@view K[1, :])
    t = collect(@view K[2, :])
    u = collect(@view K[3, :])
    return nufft3d3(x, y, z, c, -1, _FINUFFT_EPS, s, t, u)
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

end # module QuasiCrystalFINUFFTExt
