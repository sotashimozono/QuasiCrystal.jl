# Infinite-abstract quasicrystal + lazy `materialize` bridge.
#
# `QuasicrystalData` is eager: it stores every site position densely and
# is produced by running a generator to a finite cutoff. This file adds
# the *ideal infinite* counterpart for the cut-and-project families
# (Fibonacci, Penrose, Ammann–Beenker): `InfiniteQuasicrystal` carries
# only the topology marker and a site layout — no dense arrays — and
# reports `InfiniteSize`, so a caller can hold the thermodynamic-limit
# object and generate sites only on demand via
# `materialize(inf; cutoff)`.
#
# Scope: cut-and-project only. The substitution generators are keyed by
# inflation depth, `rescale` regenerates from scratch, and containment
# across generations is not guaranteed — so an "infinite substitution
# quasicrystal materialised to a cutoff" is not well defined and is
# deliberately out of scope here. `InfiniteQuasicrystal` therefore always
# materialises through the projection generator, which has a genuine
# window/radius cutoff on an ideal infinite object.

"""
    InfiniteQuasicrystal{D, T, Q, L}(topology; layout = UniformLayout(IsingSite()))

The ideal infinite cut-and-project quasicrystal of family `Q` (a
[`FibonacciLattice`](@ref), [`PenroseP3`](@ref), or
[`AmmannBeenker`](@ref) marker), in `D` physical dimensions. It carries
no site data — `size_trait` is `InfiniteSize` and
[`num_sites`](@ref)/[`position`](@ref)/[`neighbors`](@ref) throw — and
exists so the thermodynamic-limit object can be passed around and
generated lazily:

```julia
inf = InfiniteQuasicrystal(FibonacciLattice())
d   = materialize(inf; n_points = 200)     # finite QuasicrystalData
```

The cutoff keyword is family-specific, matching the underlying
generator: `n_points` for Fibonacci, `radius` for Penrose and
Ammann–Beenker. Positions use `Float64`, matching the generators.

Being the ideal infinite object, it is also the natural home for the
self-similar RG maps [`inflate`](@ref) / [`deflate`](@ref), stored
lazily as an integer `inflation` power and applied only on
`materialize`:

```julia
big = materialize(inflate(inf); n_points = 100)   # points scaled by ϕ
```

Only cut-and-project generation is modelled; substitution/inflation is
a finite construction keyed by depth, not an infinite-abstract object,
so it is out of scope for this type.
"""
struct InfiniteQuasicrystal{
    D,T<:AbstractFloat,Q<:AbstractQuasicrystal{D},L<:AbstractSiteLayout
} <: AbstractLattice{D,T}
    topology::Q
    layout::L
    inflation::Int
end

function InfiniteQuasicrystal(
    topology::Q; layout::AbstractSiteLayout=UniformLayout(IsingSite()), inflation::Integer=0
) where {D,Q<:AbstractQuasicrystal{D}}
    return InfiniteQuasicrystal{D,Float64,Q,typeof(layout)}(
        topology, layout, Int(inflation)
    )
end

# ---- RG inflate / deflate -------------------------------------------
#
# A cut-and-project quasicrystal is exactly self-similar: with λ the
# Perron eigenvalue of the family's substitution, `λ·S ⊆ S` (verified to
# ~1e-14 for the Fibonacci and Penrose generators — every point of the
# scaled set is a genuine point of the base set). Inflation is therefore
# an exact self-map on the *infinite* object — no finite-patch
# containment caveat, unlike `rescale` on a materialised
# `QuasicrystalData`. It is represented lazily as an integer power on the
# abstract lattice and applied only at `materialize` time by scaling
# positions by `λ^inflation`.
#
# Ammann–Beenker is excluded on purpose: its shipped projection
# generator is not exactly self-similar (`λ·S ⊄ S` numerically for every
# candidate λ), so an exact inflation self-map is not well defined.

"""
    inflation_factor(topo::AbstractQuasicrystal) → Float64

The linear inflation factor λ (Perron eigenvalue of the family's
substitution): the golden ratio ϕ for [`FibonacciLattice`](@ref) and
[`PenroseP3`](@ref), and the silver ratio `1 + √2` for
[`AmmannBeenker`](@ref).

For each shipped cut-and-project generator `λ·S ⊆ S` holds exactly, so
[`inflate`](@ref) / [`deflate`](@ref) are well-posed.
"""
inflation_factor(::FibonacciLattice) = (1 + sqrt(5)) / 2
inflation_factor(::PenroseP3) = (1 + sqrt(5)) / 2
inflation_factor(::AmmannBeenker) = 1 + sqrt(2)

"""
    inflate(inf::InfiniteQuasicrystal, k::Integer = 1) → InfiniteQuasicrystal

Apply `k` RG inflations: return the same infinite quasicrystal viewed at
a length scale `λ^k` larger (λ = [`inflation_factor`](@ref)). Because the
cut-and-project set is exactly self-similar (`λ·S ⊆ S`), the inflated
lattice's materialised points are genuine points of the same set scaled
by `λ^k`. [`deflate`](@ref) is the inverse. Topology and layout carry
over.

Supported for Fibonacci and Penrose; throws for Ammann–Beenker.
"""
function inflate(inf::InfiniteQuasicrystal, k::Integer=1)
    inflation_factor(inf.topology)   # validate the family supports inflation
    return InfiniteQuasicrystal(
        inf.topology; layout=inf.layout, inflation=inf.inflation + k
    )
end

"""
    deflate(inf::InfiniteQuasicrystal, k::Integer = 1) → InfiniteQuasicrystal

Apply `k` RG deflations, the inverse of [`inflate`](@ref): the same
infinite quasicrystal at a length scale `λ^k` finer.
"""
deflate(inf::InfiniteQuasicrystal, k::Integer=1) = inflate(inf, -k)

# Scale a freshly generated patch by `λ^inflation`. `inflation == 0`
# (the common case, and the only one reachable for Ammann–Beenker) is a
# no-op that never touches `inflation_factor`.
function _apply_inflation(d::QuasicrystalData{D,T}, inf::InfiniteQuasicrystal) where {D,T}
    inf.inflation == 0 && return d
    s = T(inflation_factor(inf.topology))^inf.inflation
    scaled = [s .* p for p in d.positions]
    return QuasicrystalData{D,T}(
        d.topology, scaled, d.tiles, d.generation_method, d.parameters; layout=d.layout
    )
end

# ---- Traits ---------------------------------------------------------

LatticeCore.size_trait(::InfiniteQuasicrystal) = InfiniteSize()
LatticeCore.periodicity(::InfiniteQuasicrystal) = Aperiodic()
LatticeCore.reciprocal_support(::InfiniteQuasicrystal) = HasFourierModule()
LatticeCore.topology(::InfiniteQuasicrystal) = TopologyTrait{:quasiperiodic}()
LatticeCore.site_layout(l::InfiniteQuasicrystal) = l.layout

function LatticeCore.boundary(::InfiniteQuasicrystal{D}) where {D}
    return LatticeBoundary(ntuple(_ -> OpenAxis(), D), NoModifier())
end

# ---- Linear-index API is undefined for the infinite object ----------

function LatticeCore.num_sites(::InfiniteQuasicrystal)
    return throw(
        DomainError(
            InfiniteSize(),
            "InfiniteQuasicrystal has no finite site count; call " *
            "`materialize(inf; cutoff)` (n_points for Fibonacci, radius otherwise)",
        ),
    )
end

function LatticeCore.position(::InfiniteQuasicrystal, ::Int)
    return throw(
        DomainError(
            InfiniteSize(),
            "InfiniteQuasicrystal has no site positions until materialised; call " *
            "`materialize(inf; cutoff)`",
        ),
    )
end

function LatticeCore.neighbors(::InfiniteQuasicrystal, ::Int)
    return throw(
        DomainError(
            InfiniteSize(),
            "InfiniteQuasicrystal has no neighbours until materialised; call " *
            "`materialize(inf; cutoff)`",
        ),
    )
end

# ---- Lazy materialisation through the projection generators ---------

# Re-attach the infinite lattice's layout to a freshly generated
# `QuasicrystalData` (the generators build with the default layout). The
# position and tile vectors are shared, not copied.
function _relayout(d::QuasicrystalData{D,T}, layout::AbstractSiteLayout) where {D,T}
    return QuasicrystalData{D,T}(
        d.topology, d.positions, d.tiles, d.generation_method, d.parameters; layout=layout
    )
end

"""
    materialize(inf::InfiniteQuasicrystal{1, T, FibonacciLattice}; n_points = 100)
        → QuasicrystalData

Generate a finite Fibonacci chain of `n_points` sites by cut-and-project
([`generate_fibonacci_projection`](@ref)), carrying `inf`'s site layout.
The returned lattice has empty bonds; run
[`build_nearest_neighbor_bonds!`](@ref) before walking the graph, as for
any generated `QuasicrystalData`.
"""
function LatticeCore.materialize(
    inf::InfiniteQuasicrystal{1,T,FibonacciLattice}; n_points::Int=100
) where {T}
    return _apply_inflation(
        _relayout(generate_fibonacci_projection(n_points), inf.layout), inf
    )
end

"""
    materialize(inf::InfiniteQuasicrystal{2, T, PenroseP3}; radius = 5.0)
        → QuasicrystalData

Generate a finite Penrose (P3) point set within `radius` by
cut-and-project ([`generate_penrose_projection`](@ref)), carrying
`inf`'s site layout.
"""
function LatticeCore.materialize(
    inf::InfiniteQuasicrystal{2,T,PenroseP3}; radius::Real=5.0
) where {T}
    return _apply_inflation(_relayout(generate_penrose_projection(radius), inf.layout), inf)
end

"""
    materialize(inf::InfiniteQuasicrystal{2, T, AmmannBeenker}; radius = 5.0)
        → QuasicrystalData

Generate a finite Ammann–Beenker point set within `radius` by
cut-and-project ([`generate_ammann_beenker_projection`](@ref)), carrying
`inf`'s site layout.
"""
function LatticeCore.materialize(
    inf::InfiniteQuasicrystal{2,T,AmmannBeenker}; radius::Real=5.0
) where {T}
    return _apply_inflation(
        _relayout(generate_ammann_beenker_projection(radius), inf.layout), inf
    )
end
