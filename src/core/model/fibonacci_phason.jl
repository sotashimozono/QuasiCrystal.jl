"""
Phason-space utilities for the Fibonacci chain.

These are the **geometric** primitives that any quasi-1D Fibonacci
model needs but that do not belong inside a specific Hamiltonian
or variational ansatz: the binary L/S word produced by the
substitution `L → LS, S → L`; the phason orbit on the torus
`R / Z`; the integer permutation that an irrational shift induces
on a uniform `M`-bin θ-grid; the inflation matrix; and the Fourier
expansion of an L/S step-function coupling `J(θ)` keyed to that
torus.

They are dispatched on [`FibonacciLattice`](@ref) so that the
analogous internal-space utilities for `PenroseP3` and
`AmmannBeenker` can be added later without reusing names.
"""

# ---- Default phason intercept --------------------------------------

"""
    PHASON_INTERCEPT_FIBONACCI :: Float64

Default phason intercept α = 1/ϕ used by the Fibonacci chain. Every
function in this file takes `α` as a keyword argument and defaults
to this constant, so callers can sweep other slopes for
experiments without code changes.
"""
const PHASON_INTERCEPT_FIBONACCI = inv(ϕ)

# ---- Substitution word ---------------------------------------------

"""
    fibonacci_word(::FibonacciLattice, gen::Int) → Vector{Int}

Substitution word `L → LS, S → L` after `gen` iterations starting
from the axiom `[L]`, returned as a binary vector with `0 = L`
and `1 = S`. Length is `F_{gen + 2}` (Fibonacci numbers indexed by
`F_1 = F_2 = 1`).

```jldoctest
julia> using QuasiCrystal

julia> fibonacci_word(FibonacciLattice(), 0)
1-element Vector{Int64}:
 0

julia> fibonacci_word(FibonacciLattice(), 3)
5-element Vector{Int64}:
 0
 1
 0
 0
 1
```
"""
function fibonacci_word(::FibonacciLattice, gen::Int)
    seq = [0]
    for _ in 1:gen
        new_seq = Int[]
        for s in seq
            if s == 0
                push!(new_seq, 0, 1)
            else
                push!(new_seq, 0)
            end
        end
        seq = new_seq
    end
    return seq
end

# ---- Phason orbit on the torus -------------------------------------

"""
    phason_orbit_at(::FibonacciLattice, i::Int;
                    θ0::Real = 0.0,
                    α::Real  = PHASON_INTERCEPT_FIBONACCI) → Float64

Phason coordinate of site `i` on the unit torus `R / Z`:

```math
θ_i \\;=\\; \\{ θ_0 + i\\,α \\} \\;∈\\; [0, 1).
```

For α irrational (the Fibonacci case α = 1/ϕ) the orbit `{θ_i}` is
equidistributed in `[0,1)`.
"""
function phason_orbit_at(
    ::FibonacciLattice, i::Int; θ0::Real=0.0, α::Real=PHASON_INTERCEPT_FIBONACCI
)
    return mod(float(θ0) + i * float(α), 1.0)
end

# ---- Bond coupling pattern -----------------------------------------

"""
    bond_couplings(::FibonacciLattice, JL::Real, JS::Real;
                   gen::Int) → Vector{Float64}

Per-bond coupling array of length `length(word) - 1` for the
substitution word at the requested generation `gen`. Bond `i`
takes the value `JL` when `word[i] == 0` (L bond) and `JS`
otherwise.

This is the discrete-site version of the L/S step function
`J(θ)` whose Fourier expansion is exposed by
[`J_fourier_coeffs`](@ref).
"""
function bond_couplings(lat::FibonacciLattice, JL::Real, JS::Real; gen::Int)
    seq = fibonacci_word(lat, gen)
    return Float64[seq[i] == 0 ? Float64(JL) : Float64(JS) for i in 1:(length(seq) - 1)]
end

# ---- Discrete phason-grid shift ------------------------------------

"""
    phason_grid_shift(::FibonacciLattice, M::Int;
                      α::Real = PHASON_INTERCEPT_FIBONACCI) → Vector{Int}

Permutation `σ : 1:M → 1:M` such that
`θ_{σ(k)} ≈ θ_k + α  (mod 1)` on the uniform M-bin grid
`θ_k = (k - 1) / M`. Each entry is computed by rounding the
shifted continuous coordinate to the nearest grid index and
wrapping with `mod1`.

For irrational `α` this is the discretisation that an `M`-binned
phason-space cocycle uses to advance one site.
"""
function phason_grid_shift(::FibonacciLattice, M::Int; α::Real=PHASON_INTERCEPT_FIBONACCI)
    return map(1:M) do k
        mod1(round(Int, mod((k - 1) / M + α, 1.0) * M) + 1, M)
    end
end

# ---- Inflation matrix ----------------------------------------------

"""
    inflation_matrix(::FibonacciLattice) → SMatrix{2,2,Int,4}

`M = ((1,1),(1,0))`, the abelianisation of the Fibonacci
substitution. Eigenvalues are `(ϕ, -1/ϕ)`; the parallel /
perpendicular projectors of the cut-and-project description are
the eigenvectors of `M`.
"""
inflation_matrix(::FibonacciLattice) = SMatrix{2,2,Int}(1, 1, 1, 0)

# ---- Step-function Fourier coefficients ----------------------------

"""
    _J_fourier_intervals(::FibonacciLattice, α::Real)
        → (L_intervals, S_intervals)

L/S partition of the phason torus that `J_fourier_coeffs` integrates
over. For α ∈ (0, 1) the partition is

  - L: `[0, 1 - α) ∪ [2 - 2α, 1)`
  - S: `[1 - α, 2 - 2α)`

This matches the 2-step orbital labeling used downstream by the
Fibonacci VUMPS cocycle (the bond between sites `i` and `i+1`
carries the phason midpoint `θ_i + α / 2 mod 1`).
"""
function _J_fourier_intervals(::FibonacciLattice, α::Real)
    αf = float(α)
    L_intervals = ((0.0, 1.0 - αf), (2.0 - 2αf, 1.0))
    S_intervals = ((1.0 - αf, 2.0 - 2αf),)
    return L_intervals, S_intervals
end

"""
    _indicator_fourier_coefficient(a::Real, b::Real, k::Int) → ComplexF64

`∫_a^b e^{-2πi k θ} dθ` for the indicator of `[a, b)`. At `k = 0`
this is `b - a` (mean contribution); otherwise it is the closed
form `(e^{-2πi k a} - e^{-2πi k b}) / (2πi k)`.
"""
function _indicator_fourier_coefficient(a::Real, b::Real, k::Int)
    if k == 0
        return ComplexF64(b - a)
    end
    kk = 2π * im * k
    return (exp(-kk * a) - exp(-kk * b)) / kk
end

"""
    J_fourier_coeffs(::FibonacciLattice, JL::Real, JS::Real, K_J::Int;
                     α::Real = PHASON_INTERCEPT_FIBONACCI)
        → Vector{ComplexF64}

Fourier coefficients `Ĵ_k` for `k ∈ -K_J : K_J` of the L/S step
function `J(θ)` on the unit phason torus, under the convention

```math
J(θ) \\;=\\; \\sum_{k} \\hat J_k\\, e^{2\\pi i k θ}.
```

The partition is taken from [`_J_fourier_intervals`](@ref), so

```math
\\hat J_k \\;=\\; J_L \\sum_{(a,b) \\in L} \\hat χ_{[a,b)}(k)
            \\;+\\; J_S \\sum_{(a,b) \\in S} \\hat χ_{[a,b)}(k),
```

where `χ̂_{[a,b)}(k)` is the per-interval coefficient from
[`_indicator_fourier_coefficient`](@ref). The DC term reproduces
the orbit average `⟨J⟩ = α·J_L + (1-α)·J_S` (mean Birkhoff sum).
"""
function J_fourier_coeffs(
    lat::FibonacciLattice, JL::Real, JS::Real, K_J::Int; α::Real=PHASON_INTERCEPT_FIBONACCI
)
    L_intervals, S_intervals = _J_fourier_intervals(lat, α)
    JL_f = Float64(JL)
    JS_f = Float64(JS)
    return ComplexF64[
        JL_f * sum(_indicator_fourier_coefficient(a, b, k) for (a, b) in L_intervals) +
        JS_f * sum(_indicator_fourier_coefficient(a, b, k) for (a, b) in S_intervals) for
        k in (-K_J):K_J
    ]
end
