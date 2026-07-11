"""
    QuasiCrystal localization diagnostics
    =====================================

Finite-size scaling of the inverse participation ratio (IPR), the standard
numerical probe of Anderson / quasiperiodic localization used in issue #26
(phason-strain localization of the tight-binding spectrum).

The mean IPR over the spectrum scales with system size `N` as

```math
\\overline{\\mathrm{IPR}}(N) \\sim N^{-τ},
```

and the exponent `τ` classifies the states:

| regime      | `τ`        | meaning                                   |
|:------------|:-----------|:------------------------------------------|
| extended    | `≈ 1`      | weight spread over `∝ N` sites            |
| critical    | `0 < τ < 1`| multifractal (perfect quasicrystal)       |
| localized   | `≈ 0`      | weight on `O(1)` sites, `N`-independent   |

Build a size-ordered family of Hamiltonians (uniform, phason-modulated
hoppings, or with phason-disordered on-site energies), take the mean IPR of
each, and read off `τ` with [`ipr_scaling_exponent`](@ref). These reuse the
tight-binding core ([`tight_binding_hamiltonian`](@ref),
[`inverse_participation_ratio`](@ref)).
"""

# ---- mean IPR of a spectrum -----------------------------------------

"""
    mean_inverse_participation_ratio(H::AbstractMatrix; energy_window=nothing)
        -> Float64

Mean inverse participation ratio over the eigenstates of the (real symmetric)
tight-binding Hamiltonian `H`. When `energy_window = (Emin, Emax)` is given,
only eigenstates with energy in that closed interval are averaged (useful for
probing a mobility edge); it is an error for the window to be empty.
"""
function mean_inverse_participation_ratio(H::AbstractMatrix; energy_window=nothing)
    F = eigen(Symmetric(Matrix(H)))
    sel = _window_indices(F.values, energy_window)
    s = 0.0
    @inbounds for n in sel
        s += inverse_participation_ratio(@view F.vectors[:, n])
    end
    return s / length(sel)
end

"""
    mean_inverse_participation_ratio(qc::QuasicrystalData; t=1.0, onsite=0.0,
                                     energy_window=nothing) -> Float64

Convenience method that builds the uniform-hopping (or on-site-modulated)
tight-binding Hamiltonian of `qc` and returns its mean IPR. For
phason-modulated *hoppings*, build the Hamiltonian with the per-bond
[`tight_binding_hamiltonian`](@ref) method and pass it to the matrix form.
"""
function mean_inverse_participation_ratio(
    qc::QuasicrystalData; t::Real=1.0, onsite=0.0, energy_window=nothing
)
    H = tight_binding_hamiltonian(qc; t=t, onsite=onsite)
    return mean_inverse_participation_ratio(H; energy_window=energy_window)
end

function _window_indices(vals::AbstractVector, ::Nothing)
    return eachindex(vals)
end
function _window_indices(vals::AbstractVector, window)
    lo, hi = window
    sel = findall(e -> lo ≤ e ≤ hi, vals)
    isempty(sel) && throw(ArgumentError("no eigenstates in energy_window $((lo, hi))"))
    return sel
end

# ---- scaling sweep & exponent ---------------------------------------

"""
    ipr_scaling(qcs; t=1.0, onsite=0.0, energy_window=nothing)
        -> (sizes::Vector{Int}, mean_iprs::Vector{Float64})

Mean IPR of each lattice in the size-ordered collection `qcs` (uniform-hopping
model), paired with `num_sites`. Feed the result to
[`ipr_scaling_exponent`](@ref) to obtain the localization exponent `τ`. For
phason-modulated or disordered models, compute `mean_inverse_participation_ratio`
per Hamiltonian yourself and call `ipr_scaling_exponent` directly.
"""
function ipr_scaling(qcs; t::Real=1.0, onsite=0.0, energy_window=nothing)
    sizes = Int[num_sites(qc) for qc in qcs]
    mis = Float64[
        mean_inverse_participation_ratio(
            qc; t=t, onsite=onsite, energy_window=energy_window
        ) for qc in qcs
    ]
    return sizes, mis
end

"""
    ipr_scaling_exponent(sizes, mean_iprs) -> Float64

Localization exponent `τ` from a least-squares fit of
`log(mean_iprs) = a - τ · log(sizes)`. `τ ≈ 1` for extended states, `τ ≈ 0`
for localized states, and `0 < τ < 1` for critical (multifractal) states.

Requires at least two size points, and all `mean_iprs` must be positive.
"""
function ipr_scaling_exponent(sizes, mean_iprs)
    length(sizes) == length(mean_iprs) ||
        throw(DimensionMismatch("sizes and mean_iprs have different lengths"))
    n = length(sizes)
    n ≥ 2 || throw(ArgumentError("need at least two size points, got $n"))
    all(>(0), mean_iprs) ||
        throw(ArgumentError("all mean_iprs must be positive to take a log"))
    x = log.(float.(collect(sizes)))
    y = log.(float.(collect(mean_iprs)))
    x̄ = sum(x) / n
    ȳ = sum(y) / n
    sxx = sum((xi - x̄)^2 for xi in x)
    sxy = sum((x[i] - x̄) * (y[i] - ȳ) for i in 1:n)
    slope = sxy / sxx
    return -slope
end
