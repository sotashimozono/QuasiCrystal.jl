"""
    QuasiCrystal tight-binding & spectral analysis
    ==============================================

A minimal single-particle tight-binding layer on top of any
`QuasicrystalData` whose nearest-neighbour bonds have been populated (via
[`build_nearest_neighbor_bonds!`](@ref)). It provides the shared numerical
core for electronic-structure studies on quasicrystals — density of states,
and localization diagnostics through the inverse participation ratio — that
issues #26 (phason-strain localization) and #27 (QAtlas tight-binding /
dynamic structure factor) both build on.

The Hamiltonian is

```math
H = ε \\, \\mathbb{1} - \\sum_{\\langle i,j\\rangle} t_{ij}
        \\bigl( |i\\rangle\\langle j| + |j\\rangle\\langle i| \\bigr),
```

with a uniform hopping `t` by default, or per-bond hoppings `t_{ij}` (e.g.
the `bond_couplings` of a phason model) supplied in `bonds(qc)` order. The
on-site energy `ε` may be a scalar or a per-site vector.
"""

# ---- Hamiltonian -----------------------------------------------------

"""
    tight_binding_hamiltonian(qc::QuasicrystalData; t=1.0, onsite=0.0)
        -> SparseMatrixCSC{Float64,Int}

Build the real symmetric single-particle tight-binding Hamiltonian of `qc`
with uniform hopping `t` on every nearest-neighbour bond. `onsite` is either
a scalar on-site energy applied to every site or a length-`num_sites(qc)`
vector of site energies.

Requires that the bonds of `qc` have been populated
([`build_nearest_neighbor_bonds!`](@ref)); an all-zero (bond-less)
Hamiltonian otherwise.

See also [`spectrum`](@ref), [`inverse_participation_ratio`](@ref),
[`density_of_states`](@ref).
"""
function tight_binding_hamiltonian(qc::QuasicrystalData; t::Real=1.0, onsite=0.0)
    bs = bonds(qc)
    hoppings = fill(Float64(t), length(bs))
    return _assemble_hamiltonian(qc, hoppings, onsite)
end

"""
    tight_binding_hamiltonian(qc::QuasicrystalData, hoppings::AbstractVector;
                              onsite=0.0) -> SparseMatrixCSC{Float64,Int}

Per-bond hopping variant: `hoppings[k]` is the hopping amplitude on the
`k`-th bond of `bonds(qc)`. Use this to feed a bond-resolved coupling list,
e.g. the `bond_couplings(FibonacciLattice(), JL, JS; gen)` of a phason model.

`length(hoppings)` must equal `num_bonds(qc)`.
"""
function tight_binding_hamiltonian(
    qc::QuasicrystalData, hoppings::AbstractVector; onsite=0.0
)
    length(hoppings) == length(bonds(qc)) || throw(
        DimensionMismatch(
            "hoppings has length $(length(hoppings)) but qc has $(length(bonds(qc))) bonds",
        ),
    )
    return _assemble_hamiltonian(qc, hoppings, onsite)
end

function _assemble_hamiltonian(qc::QuasicrystalData, hoppings::AbstractVector, onsite)
    N = num_sites(qc)
    bs = bonds(qc)
    onsite_vec = _onsite_vector(onsite, N)

    I = Int[]
    J = Int[]
    V = Float64[]
    sizehint!(I, 2 * length(bs) + N)
    sizehint!(J, 2 * length(bs) + N)
    sizehint!(V, 2 * length(bs) + N)

    @inbounds for k in eachindex(bs)
        b = bs[k]
        tij = Float64(hoppings[k])
        push!(I, b.i);
        push!(J, b.j);
        push!(V, -tij)
        push!(I, b.j);
        push!(J, b.i);
        push!(V, -tij)
    end
    @inbounds for i in 1:N
        if onsite_vec[i] != 0
            push!(I, i);
            push!(J, i);
            push!(V, onsite_vec[i])
        end
    end
    return sparse(I, J, V, N, N)
end

_onsite_vector(onsite::Real, N::Int) = fill(Float64(onsite), N)
function _onsite_vector(onsite::AbstractVector, N::Int)
    length(onsite) == N ||
        throw(DimensionMismatch("onsite vector has length $(length(onsite)), expected $N"))
    return Float64.(onsite)
end

# ---- Spectrum & eigenstates -----------------------------------------

"""
    spectrum(qc::QuasicrystalData; t=1.0, onsite=0.0) -> Vector{Float64}

Sorted eigenvalues of the tight-binding Hamiltonian of `qc`. Convenience
wrapper over `eigvals` of [`tight_binding_hamiltonian`](@ref).
"""
function spectrum(qc::QuasicrystalData; t::Real=1.0, onsite=0.0)
    H = tight_binding_hamiltonian(qc; t=t, onsite=onsite)
    return eigvals(Symmetric(Matrix(H)))
end

"""
    eigenstates(qc::QuasicrystalData; t=1.0, onsite=0.0)
        -> (values::Vector{Float64}, vectors::Matrix{Float64})

Full eigendecomposition of the tight-binding Hamiltonian: `values` sorted
ascending, `vectors[:, n]` the eigenvector for `values[n]`.
"""
function eigenstates(qc::QuasicrystalData; t::Real=1.0, onsite=0.0)
    H = tight_binding_hamiltonian(qc; t=t, onsite=onsite)
    F = eigen(Symmetric(Matrix(H)))
    return F.values, F.vectors
end

# ---- Localization: inverse participation ratio ----------------------

"""
    inverse_participation_ratio(ψ::AbstractVector) -> Float64

Inverse participation ratio of a single state,

```math
\\mathrm{IPR}(ψ) = \\frac{\\sum_i |ψ_i|^4}{\\left(\\sum_i |ψ_i|^2\\right)^2}.
```

It ranges in `[1/N, 1]`: a fully extended state (uniform amplitude) gives
`1/N`, a state localized on a single site gives `1`. The *participation
number* `1 / IPR` estimates how many sites the state occupies.
"""
function inverse_participation_ratio(ψ::AbstractVector)
    n2 = zero(real(eltype(ψ)))
    n4 = zero(real(eltype(ψ)))
    @inbounds for a in ψ
        p = abs2(a)
        n2 += p
        n4 += p * p
    end
    n2 == 0 && throw(ArgumentError("inverse_participation_ratio of the zero vector"))
    return n4 / (n2 * n2)
end

"""
    inverse_participation_ratios(qc::QuasicrystalData; t=1.0, onsite=0.0)
        -> (energies::Vector{Float64}, iprs::Vector{Float64})

Per-eigenstate IPR of the tight-binding spectrum of `qc`, paired with the
corresponding energies (ascending). Rising IPR (toward `O(1)`) as system
size grows signals localized states.
"""
function inverse_participation_ratios(qc::QuasicrystalData; t::Real=1.0, onsite=0.0)
    vals, vecs = eigenstates(qc; t=t, onsite=onsite)
    iprs = [inverse_participation_ratio(@view vecs[:, n]) for n in axes(vecs, 2)]
    return vals, iprs
end

# ---- Density of states ----------------------------------------------

"""
    density_of_states(qc::QuasicrystalData; t=1.0, onsite=0.0,
                      nbins::Int=100, broadening::Real=0.0)
        -> (centers::Vector{Float64}, dos::Vector{Float64})

Density of states of the tight-binding spectrum of `qc`, normalized so that
`sum(dos) * step ≈ num_sites(qc)` (i.e. it integrates to the total number of
states). With `broadening = 0` a plain histogram over `nbins` bins is
returned; with `broadening = σ > 0` each eigenvalue is smeared by a Gaussian
of width `σ` and sampled at the bin centres.
"""
function density_of_states(
    qc::QuasicrystalData; t::Real=1.0, onsite=0.0, nbins::Int=100, broadening::Real=0.0
)
    nbins ≥ 1 || throw(ArgumentError("nbins must be ≥ 1, got $nbins"))
    E = spectrum(qc; t=t, onsite=onsite)
    return _density_of_states(E, nbins, Float64(broadening))
end

function _density_of_states(E::AbstractVector, nbins::Int, σ::Float64)
    Emin, Emax = extrema(E)
    # Pad a degenerate spectrum so the range is non-zero.
    if Emax ≈ Emin
        Emin -= 0.5
        Emax += 0.5
    end
    # When broadening, widen the window by a few σ on each side so the
    # Gaussian weight of edge eigenvalues is captured; otherwise the DOS
    # would integrate to less than the number of states.
    if σ > 0
        Emin -= 3σ
        Emax += 3σ
    end
    edges = range(Emin, Emax; length=nbins + 1)
    dE = step(edges)
    centers = collect(edges[1:(end - 1)] .+ dE / 2)
    dos = zeros(Float64, nbins)

    if σ > 0
        norm = 1 / (σ * sqrt(2π))
        @inbounds for e in E, b in 1:nbins
            dos[b] += norm * exp(-((centers[b] - e)^2) / (2σ^2))
        end
    else
        @inbounds for e in E
            b = clamp(floor(Int, (e - Emin) / dE) + 1, 1, nbins)
            dos[b] += 1
        end
        dos ./= dE                       # counts -> density
    end
    return centers, dos
end
