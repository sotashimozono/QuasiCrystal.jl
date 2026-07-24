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

Only cut-and-project generation is modelled; substitution/inflation is
a finite construction keyed by depth, not an infinite-abstract object,
so it is out of scope for this type.
"""
struct InfiniteQuasicrystal{
    D,T<:AbstractFloat,Q<:AbstractQuasicrystal{D},L<:AbstractSiteLayout
} <: AbstractLattice{D,T}
    topology::Q
    layout::L
end

function InfiniteQuasicrystal(
    topology::Q; layout::AbstractSiteLayout=UniformLayout(IsingSite())
) where {D,Q<:AbstractQuasicrystal{D}}
    return InfiniteQuasicrystal{D,Float64,Q,typeof(layout)}(topology, layout)
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
    return _relayout(generate_fibonacci_projection(n_points), inf.layout)
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
    return _relayout(generate_penrose_projection(radius), inf.layout)
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
    return _relayout(generate_ammann_beenker_projection(radius), inf.layout)
end
