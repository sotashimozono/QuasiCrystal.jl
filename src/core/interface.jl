"""
Common interface methods for AbstractLattice types.

These methods provide unified access to lattice properties regardless of whether
the structure is a periodic lattice or an aperiodic quasicrystal.
"""

using LinearAlgebra

# Position tolerance for duplicate detection
const POSITION_TOLERANCE = 1e-10

"""
    get_positions(lattice::AbstractLattice)
Get the positions of all sites in the lattice.

# Returns
- `Vector{Vector{T}}`: positions of all sites
"""
function get_positions(data::QuasicrystalData{D,T,TT}) where {D,T,TT}
  return data.positions
end

function get_positions(lattice::Lattice{Topology,T,B,I}) where {Topology,T,B,I}
  return lattice.positions
end

"""
    get_bonds(lattice::AbstractLattice)
Get the bonds (edges) in the lattice.

# Returns
- `Vector{Bond}`: list of bonds
"""
function get_bonds(data::QuasicrystalData{D,T,TT}) where {D,T,TT}
  return data.bonds
end

function get_bonds(lattice::Lattice{Topology,T,B,I}) where {Topology,T,B,I}
  return lattice.bonds
end

"""
    get_nearest_neighbors(lattice::AbstractLattice)
Get the nearest neighbor indices for each site.

# Returns
- `Vector{Vector{Int}}`: nearest neighbor indices for each site
"""
function get_nearest_neighbors(data::QuasicrystalData{D,T,TT}) where {D,T,TT}
  return data.nearest_neighbors
end

function get_nearest_neighbors(
  lattice::Lattice{Topology,T,B,I}
) where {Topology,T,B,I}
  return lattice.nearest_neighbors
end

"""
    num_sites(lattice::AbstractLattice)
Get the total number of sites in the lattice.

# Returns
- `Int`: number of sites
"""
function num_sites(data::QuasicrystalData{D,T,TT}) where {D,T,TT}
  return length(data.positions)
end

function num_sites(lattice::Lattice{Topology,T,B,I}) where {Topology,T,B,I}
  return lattice.N
end

"""
    num_bonds(lattice::AbstractLattice)
Get the total number of bonds in the lattice.

# Returns
- `Int`: number of bonds
"""
function num_bonds(data::QuasicrystalData{D,T,TT}) where {D,T,TT}
  return length(data.bonds)
end

function num_bonds(lattice::Lattice{Topology,T,B,I}) where {Topology,T,B,I}
  return length(lattice.bonds)
end

"""
    build_nearest_neighbor_bonds!(data::QuasicrystalData{D,T,TT}; cutoff::Real) where {D,T,TT}
Build nearest neighbor bonds for a quasicrystal based on distance cutoff.
Updates the `bonds` and `nearest_neighbors` fields in place.

# Arguments
- `data`: Quasicrystal data structure
- `cutoff`: Maximum distance for nearest neighbors

# Returns
- Modified `data` with bonds and nearest neighbors populated
"""
function build_nearest_neighbor_bonds!(
  data::QuasicrystalData{D,T,TT}; cutoff::Real
) where {D,T,TT}
  n = length(data.positions)

  # Initialize nearest neighbors if empty
  if isempty(data.nearest_neighbors)
    data.nearest_neighbors = [Int[] for _ in 1:n]
  end

  # Clear existing bonds
  empty!(data.bonds)

  # Build bonds based on distance
  for i in 1:n
    for j in (i + 1):n
      pos_i = data.positions[i]
      pos_j = data.positions[j]
      dist = norm(pos_j - pos_i)

      if dist < cutoff && dist > POSITION_TOLERANCE  # Avoid duplicate positions
        # Add bond
        bond_vector = pos_j - pos_i
        push!(data.bonds, Bond(i, j, 1, bond_vector))

        # Update nearest neighbors
        push!(data.nearest_neighbors[i], j)
        push!(data.nearest_neighbors[j], i)
      end
    end
  end

  return data
end

export get_positions, get_bonds, get_nearest_neighbors
export num_sites, num_bonds
export build_nearest_neighbor_bonds!
