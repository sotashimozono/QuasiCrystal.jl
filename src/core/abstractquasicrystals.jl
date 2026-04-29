"""
Quasicrystal topology markers, generation method tags, and the
`QuasicrystalData` concrete lattice type.

After the LatticeCore migration, `AbstractLattice`, `Bond`, and
`UnitCell` are imported from LatticeCore; only QuasiCrystal-specific
abstractions live here.
"""

const GOLDEN_RATIO = (1 + sqrt(5)) / 2
const ϕ = GOLDEN_RATIO

"""
    AbstractQuasicrystal{D}

Topology marker for cut-and-project / substitution quasicrystals in
`D` physical dimensions. Unlike `LatticeCore.AbstractLattice`, this
is **not** a lattice type — it is a *dispatch key* identifying the
family of algorithms used to generate a concrete
[`QuasicrystalData`](@ref) instance. The shipped markers are
[`FibonacciLattice`](@ref), [`PenroseP3`](@ref), and
[`AmmannBeenker`](@ref).
"""
abstract type AbstractQuasicrystal{D} end

"""
    AbstractGenerationMethod

Tag type for the algorithmic strategy used to build a quasicrystal:
[`ProjectionMethod`](@ref) (cut-and-project from a higher-dimensional
periodic lattice) or [`SubstitutionMethod`](@ref) (inflation /
substitution rules).
"""
abstract type AbstractGenerationMethod end

"""Cut-and-project generation tag."""
struct ProjectionMethod <: AbstractGenerationMethod end

"""
    AbstractSubstitutionAlgorithm

Dispatch trait to select the specific deflation/inflation algorithm
during substitution generation.
"""
abstract type AbstractSubstitutionAlgorithm end

"""Fallback to the fastest / standard algorithm available."""
struct DefaultSubstitution <: AbstractSubstitutionAlgorithm end

"""Subdivide by passing through Robinson triangles (Half-kite / Half-dart)."""
struct RobinsonTriangleInflation <: AbstractSubstitutionAlgorithm end

"""Subdivide tiles directly based on vertex/edge matching."""
struct DirectTileInflation <: AbstractSubstitutionAlgorithm end

"""Substitution / inflation generation tag."""
struct SubstitutionMethod{A<:AbstractSubstitutionAlgorithm} <: AbstractGenerationMethod
    algorithm::A
end
SubstitutionMethod() = SubstitutionMethod(DefaultSubstitution())

"""
    Tile{D, T}

A single tile in a quasicrystalline tiling. Carries the vertices of
the tile, a semantic [`TileType`](@ref) tag (e.g. [`FatRhombus`](@ref)
vs [`ThinRhombus`](@ref)), and the tile centre.

!!! note "Breaking change (v0.5)"
    `type::Int` was replaced by `type::TileType`. Migrate `1/2`
    integer tags to the geometry-specific singletons defined in
    `src/core/tile_types.jl`.
"""
struct Tile{D,T}
    vertices::Vector{SVector{D,T}}
    type::TileType
    center::SVector{D,T}
end

"""
    QuasicrystalData{D, T, Topo, TileType, L} <: LatticeCore.AbstractLattice{D, T}

The concrete lattice instance returned by the `generate_*` family
of functions. Subtype of
`LatticeCore.AbstractLattice{D, T}` so that every LatticeCore
observer, trait, and test-suite helper works on quasicrystals as
well as on periodic lattices.

The `topology::Topo` field carries the topology marker singleton
(e.g. `FibonacciLattice()`, `PenroseP3()`, `AmmannBeenker()`).
Fourier-analysis entry points like
[`hyper_reciprocal_lattice`](@ref) and
[`LatticeCore.fourier_module`](@ref LatticeCore.fourier_module)
dispatch on `topology`'s concrete type to pick the right
projection matrices and acceptance window.

# Fields

- `topology::Topo` — topology marker singleton
- `positions::Vector{SVector{D, T}}` — physical positions of every
  site
- `tiles::Vector{TileType}` — list of tiles in the tiling (may be
  empty for 1D lattices or ungenerated tilings)
- `generation_method::AbstractGenerationMethod` — which algorithm
  built this instance
- `parameters::Dict{Symbol, Any}` — free-form parameter bag
- `bonds::Vector{Bond{D, T}}` — nearest-neighbour bonds, populated
  by [`build_nearest_neighbor_bonds!`](@ref)
- `nearest_neighbors::Vector{Vector{Int}}` — neighbour adjacency
  lists, one per site
- `layout::AbstractSiteLayout` — LatticeCore site layout (defaults
  to `UniformLayout(IsingSite())`)

`bonds` and `nearest_neighbors` start empty and must be filled in
with `build_nearest_neighbor_bonds!(data; cutoff)` before Monte
Carlo observers walk the graph.
"""
struct QuasicrystalData{
    D,T<:AbstractFloat,Topo<:AbstractQuasicrystal{D},TileType,L<:AbstractSiteLayout
} <: AbstractLattice{D,T}
    topology::Topo
    positions::Vector{SVector{D,T}}
    tiles::Vector{TileType}
    generation_method::AbstractGenerationMethod
    parameters::Dict{Symbol,Any}
    bonds::Vector{Bond{D,T}}
    nearest_neighbors::Vector{Vector{Int}}
    layout::L
end

# ---- Convenience constructors ---------------------------------------

function QuasicrystalData{D,T}(
    topology::Topo,
    positions::Vector{SVector{D,T}},
    tiles::Vector{TT},
    method::AbstractGenerationMethod,
    params::Dict{Symbol,Any};
    layout::AbstractSiteLayout=UniformLayout(IsingSite()),
) where {D,T<:AbstractFloat,Topo<:AbstractQuasicrystal{D},TT}
    n = length(positions)
    bonds = Bond{D,T}[]
    nn = [Int[] for _ in 1:n]
    return QuasicrystalData{D,T,Topo,TT,typeof(layout)}(
        topology, positions, tiles, method, params, bonds, nn, layout
    )
end

# ---- LatticeCore required interface ---------------------------------

LatticeCore.num_sites(data::QuasicrystalData) = length(data.positions)

LatticeCore.position(data::QuasicrystalData, i::Int) = data.positions[i]

function LatticeCore.boundary(::QuasicrystalData{D}) where {D}
    LatticeBoundary(ntuple(_ -> OpenAxis(), D), NoModifier())
end

LatticeCore.site_layout(data::QuasicrystalData) = data.layout

function LatticeCore.size_trait(data::QuasicrystalData)
    return FiniteSize((length(data.positions),))
end

LatticeCore.bonds(data::QuasicrystalData) = data.bonds

# Quasicrystals are aperiodic by construction: they admit no
# Bravais reciprocal lattice (they have a dense Fourier module
# instead, which is handled by `LatticeCore.fourier_module`).
LatticeCore.periodicity(::QuasicrystalData) = Aperiodic()
LatticeCore.reciprocal_support(::QuasicrystalData) = HasFourierModule()
