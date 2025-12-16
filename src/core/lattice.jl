"""
Periodic lattice implementation with unit cell structure.
"""

"""
    AbstractBoundaryCondition
Abstract type for boundary conditions.
"""
abstract type AbstractBoundaryCondition end

"""
    PeriodicBoundary <: AbstractBoundaryCondition
Periodic boundary condition.
"""
struct PeriodicBoundary <: AbstractBoundaryCondition end

"""
    OpenBoundary <: AbstractBoundaryCondition
Open boundary condition (no wrapping).
"""
struct OpenBoundary <: AbstractBoundaryCondition end

"""
    AbstractIndexing
Abstract type for indexing methods.
"""
abstract type AbstractIndexing end

"""
    LinearIndexing <: AbstractIndexing
Linear indexing method (sites numbered 1, 2, 3, ..., N).
"""
struct LinearIndexing <: AbstractIndexing end

"""
    CartesianIndexing <: AbstractIndexing
Cartesian indexing method (sites indexed by (i, j) coordinates).
"""
struct CartesianIndexing <: AbstractIndexing end

"""
    Lattice{Topology<:AbstractTopology, T, B<:AbstractBoundaryCondition, I<:AbstractIndexing}
Represents a periodic lattice constructed by repeating a unit cell.
Mainly for 2-dimensional lattices, but can be used for 1-dimensional lattices as well.

# Fields
- `Lx::Int`: x direction lattice size (number of unit cells)
- `Ly::Int`: y direction lattice size (number of unit cells)
- `N::Int`: total number of sites
- `positions::Vector{Vector{T}}`: position vectors of each site
- `nearest_neighbors::Vector{Vector{Int}}`: nearest neighbor indices for each site
- `bonds::Vector{Bond}`: list of bonds (edges) in the lattice
- `basis_vectors::Vector{Vector{T}}`: lattice basis vectors
- `reciprocal_vectors::Union{Vector{Vector{T}}, Nothing}`: reciprocal lattice vectors
- `sublattice_ids::Vector{Int}`: sublattice IDs of each site
- `is_bipartite::Bool`: whether the lattice is bipartite
- `site_map::Union{Matrix{Int}, Nothing}`: mapping of site indices on the lattice
- `translation_x::Vector{Int}`: x direction translation vector
- `translation_y::Vector{Int}`: y direction translation vector
- `boundary::B`: boundary condition
- `index_method::I`: indexing method
"""
struct Lattice{
  Topology<:AbstractTopology,T,B<:AbstractBoundaryCondition,I<:AbstractIndexing
} <: AbstractLattice{2}
  Lx::Int
  Ly::Int
  N::Int
  positions::Vector{Vector{T}}
  # Graph representation
  nearest_neighbors::Vector{Vector{Int}}
  bonds::Vector{Bond}
  # Topology information
  basis_vectors::Vector{Vector{T}}
  reciprocal_vectors::Union{Vector{Vector{T}},Nothing}
  sublattice_ids::Vector{Int}
  is_bipartite::Bool
  site_map::Union{Matrix{Int},Nothing}
  translation_x::Vector{Int}
  translation_y::Vector{Int}
  boundary::B
  index_method::I
end

export AbstractBoundaryCondition, PeriodicBoundary, OpenBoundary
export AbstractIndexing, LinearIndexing, CartesianIndexing
export Lattice
