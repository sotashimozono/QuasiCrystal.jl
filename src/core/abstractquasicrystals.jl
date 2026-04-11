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

"""Substitution / inflation generation tag."""
struct SubstitutionMethod <: AbstractGenerationMethod end

"""
    Tile{D, T}

A single tile in a quasicrystalline tiling. Carries the vertices of
the tile, an integer type id (e.g. fat vs thin rhombus), and the
tile centre.
"""
struct Tile{D,T}
    vertices::Vector{SVector{D,T}}
    type::Int
    center::SVector{D,T}
end

"""
    QuasicrystalData{D, T, TileType} <: LatticeCore.AbstractLattice{D, T}

The concrete lattice instance returned by the `generate_*` family
of functions. Subtype of
[`LatticeCore.AbstractLattice`](@ref LatticeCore.AbstractLattice)`{D, T}`
so that every LatticeCore observer, trait, and test-suite helper
works on quasicrystals as well as on periodic lattices.

# Fields

- `positions::Vector{SVector{D, T}}` — physical positions of every
  site
- `tiles::Vector{TileType}` — list of tiles in the tiling (may be
  empty for 1D lattices or ungenerated tilings)
- `generation_method::AbstractGenerationMethod` — which algorithm
  built this instance
- `parameters::Dict{Symbol, Any}` — free-form parameter bag kept
  for backwards compatibility with the pre-migration API
- `bonds::Vector{Bond{D, T}}` — nearest-neighbour bonds, populated
  by [`build_nearest_neighbor_bonds!`](@ref) or by generation
  functions that know the local connectivity
- `nearest_neighbors::Vector{Vector{Int}}` — neighbour adjacency
  lists, one per site
- `layout::AbstractSiteLayout` — LatticeCore site layout (defaults
  to `UniformLayout(IsingSite())`)

`bonds` and `nearest_neighbors` start empty and must be filled in
with `build_nearest_neighbor_bonds!(data; cutoff)` before Monte
Carlo observers walk the graph.
"""
struct QuasicrystalData{D,T<:AbstractFloat,TileType,L<:AbstractSiteLayout} <:
       AbstractLattice{D,T}
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
    positions::Vector{SVector{D,T}},
    tiles::Vector{TT},
    method::AbstractGenerationMethod,
    params::Dict{Symbol,Any};
    layout::AbstractSiteLayout=UniformLayout(IsingSite()),
) where {D,T<:AbstractFloat,TT}
    n = length(positions)
    bonds = Bond{D,T}[]
    nn = [Int[] for _ in 1:n]
    return QuasicrystalData{D,T,TT,typeof(layout)}(
        positions, tiles, method, params, bonds, nn, layout
    )
end

# ---- LatticeCore required interface ---------------------------------

LatticeCore.num_sites(data::QuasicrystalData) = length(data.positions)

LatticeCore.position(data::QuasicrystalData, i::Int) = data.positions[i]

LatticeCore.neighbors(data::QuasicrystalData, i::Int) = data.nearest_neighbors[i]

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
