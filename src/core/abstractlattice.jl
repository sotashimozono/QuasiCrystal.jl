"""
Abstract lattice interface for unified access to periodic lattices and quasicrystals.

This module provides a common interface for both periodic lattices (with unit cell structure)
and aperiodic quasicrystals, enabling applications like Lattice2DMonteCarlo to work uniformly
with either structure type through consistent access to bonds, sites, and indexing.
"""

"""
    AbstractLattice{D}
Abstract type for lattices in D dimensions.
This is the base type for all lattice-like structures including both periodic lattices
and aperiodic quasicrystals.
"""
abstract type AbstractLattice{D} end

"""
    AbstractLatticeConnection
Abstract type for lattice connections (edges, bonds).
"""
abstract type AbstractLatticeConnection end

"""
    Bond <: AbstractLatticeConnection
Struct for a bond (edge) in the lattice.
- `src::Int`: start site index
- `dst::Int`: destination site index
- `type::Int`: type of the bond for categorization
- `vector::Vector{Float64}`: expresses the bond vector from src to dst
"""
struct Bond <: AbstractLatticeConnection
    src::Int
    dst::Int
    type::Int
    vector::Vector{Float64}
end

"""
    Connection <: AbstractLatticeConnection
Connection rules within or between unit cells (used for periodic lattices).
- `src_sub::Int`: sublattice index of the start point (1, 2, ...)
- `dst_sub::Int`: sublattice index of the end point
- `dx::Int`: relative cell position in x direction (0 means within the same unit cell)
- `dy::Int`: relative cell position in y direction (0 means within the same unit cell)
- `type::Int`: type of the connection
"""
struct Connection <: AbstractLatticeConnection
    src_sub::Int
    dst_sub::Int
    dx::Int
    dy::Int
    type::Int
end

"""
    UnitCell{D,T}
Geometric definition data of periodic lattices.
The lattice is constructed based on repeating this unit cell.
- `basis::Vector{Vector{T}}`: lattice basis vectors
- `sublattice_positions::Vector{Vector{T}}`: positions of sublattice sites within the unit cell
- `connections::Vector{Connection}`: connection rules between sites
"""
struct UnitCell{D,T} <: AbstractLattice{D}
    basis::Vector{Vector{T}}
    sublattice_positions::Vector{Vector{T}}
    connections::Vector{Connection}
end

export AbstractLattice, AbstractLatticeConnection
export Bond, Connection, UnitCell
