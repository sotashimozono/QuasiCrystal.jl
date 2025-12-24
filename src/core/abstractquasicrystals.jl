"""
Quasicrystal structures and generation methods.
This file contains the infrastructure for generating quasicrystalline patterns.
"""

const GOLDEN_RATIO = (1 + sqrt(5)) / 2
const ϕ = GOLDEN_RATIO

"""
    AbstractQuasicrystal{D}
Abstract type for quasicrystal topologies in D dimensions.
Unlike periodic lattices, quasicrystals lack translational symmetry but have long-range order.
Inherits from AbstractLattice to provide a unified interface with periodic lattices.
"""
abstract type AbstractQuasicrystal{D} <: AbstractLattice{D} end

"""
    AbstractGenerationMethod
Abstract type for quasicrystal generation methods.
Different methods can be used to generate the same quasicrystal pattern.
"""
abstract type AbstractGenerationMethod end

"""
    ProjectionMethod <: AbstractGenerationMethod
Generate quasicrystals via projection from higher-dimensional periodic lattices.
This is the most common theoretical approach.
"""
struct ProjectionMethod <: AbstractGenerationMethod end

"""
    SubstitutionMethod <: AbstractGenerationMethod
Generate quasicrystals via substitution (inflation) rules.
This is an alternative algorithmic approach.
"""
struct SubstitutionMethod <: AbstractGenerationMethod end

"""
    QuasicrystalData{D,T,TileType}
Data structure holding the generated quasicrystal pattern.
- `positions::Vector{Vector{T}}`: positions of vertices/sites
- `tiles::Vector{TileType}`: list of tiles in the pattern
- `generation_method::AbstractGenerationMethod`: method used to generate
- `parameters::Dict{Symbol,Any}`: generation parameters
- `bonds::Vector{Bond}`: list of bonds connecting sites (optional, empty by default)
- `nearest_neighbors::Vector{Vector{Int}}`: nearest neighbor indices for each site (optional, empty by default)
"""
struct QuasicrystalData{D,T,TileType}
    positions::Vector{Vector{T}}
    tiles::Vector{TileType}
    generation_method::AbstractGenerationMethod
    parameters::Dict{Symbol,Any}
    bonds::Vector{Bond}
    nearest_neighbors::Vector{Vector{Int}}
end

# Helper function to initialize empty bonds and nearest neighbors
function _init_empty_connectivity(n_positions::Int)
    bonds = Bond[]
    nearest_neighbors = Vector{Int}[Int[] for _ in 1:n_positions]
    return bonds, nearest_neighbors
end

# Convenience constructor that infers TileType
function QuasicrystalData{D,T}(
    positions::Vector{Vector{T}},
    tiles::Vector{TT},
    method::AbstractGenerationMethod,
    params::Dict{Symbol,Any},
) where {D,T,TT}
    # Default: no bonds or nearest neighbors
    bonds, nearest_neighbors = _init_empty_connectivity(length(positions))
    return QuasicrystalData{D,T,TT}(
        positions, tiles, method, params, bonds, nearest_neighbors
    )
end

# Constructor with bonds and nearest neighbors
function QuasicrystalData{D,T}(
    positions::Vector{Vector{T}},
    tiles::Vector{TT},
    method::AbstractGenerationMethod,
    params::Dict{Symbol,Any},
    bonds::Vector{Bond},
    nearest_neighbors::Vector{Vector{Int}},
) where {D,T,TT}
    return QuasicrystalData{D,T,TT}(
        positions, tiles, method, params, bonds, nearest_neighbors
    )
end

"""
    Tile{D,T}
Represents a single tile in the quasicrystal pattern.
- `vertices::Vector{Vector{T}}`: corner positions of the tile
- `type::Int`: tile type identifier (e.g., fat vs thin rhombus)
- `center::Vector{T}`: center position of the tile
"""
struct Tile{D,T}
    vertices::Vector{Vector{T}}
    type::Int
    center::Vector{T}
end

export AbstractQuasicrystal, AbstractGenerationMethod
export ProjectionMethod, SubstitutionMethod
export QuasicrystalData, Tile
export GOLDEN_RATIO, ϕ
