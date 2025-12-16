"""
Topology definitions for periodic lattices.
"""

"""
    AbstractTopology{D}
Abstract type for lattice topologies in D dimensions.
Each concrete topology defines a specific lattice structure (e.g., Square, Triangular, Honeycomb).
"""
abstract type AbstractTopology{D} <: AbstractLattice{D} end

"""
    get_unit_cell(::Type{T}) where T <: AbstractTopology
Returns the UnitCell associated with the given Topology type.
If the Topology type is not recognized, it throws an error.

# Returns
- `UnitCell{D,Float64}`: Unit cell definition for the topology
"""
function get_unit_cell(::Type{T}) where {T<:AbstractTopology}
  return error("UnitCell not defined for $T")
end

"""
    Square <: AbstractTopology{2}
Square lattice topology.
A simple 2D lattice with one site per unit cell and square geometry.
"""
struct Square <: AbstractTopology{2} end

"""
    get_unit_cell(::Type{Square})
Returns the unit cell for a square lattice.

# Returns
- Unit cell with basis vectors a1=[1,0], a2=[0,1], one sublattice site at origin,
  and connections to nearest neighbors in x and y directions.
"""
function get_unit_cell(::Type{Square})
  a1 = [1.0, 0.0]
  a2 = [0.0, 1.0]
  # Connection(src_sub, dst_sub, dx, dy, type)
  # type 1: horizontal bond, type 2: vertical bond
  conns = [Connection(1, 1, 1, 0, 1), Connection(1, 1, 0, 1, 2)]
  return UnitCell{2,Float64}([a1, a2], [[0.0, 0.0]], conns)
end

export AbstractTopology, get_unit_cell, Square
