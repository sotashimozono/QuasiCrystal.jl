"""
Backwards-compat accessors for `QuasicrystalData` plus the
distance-based nearest-neighbour builder.

After the LatticeCore migration the canonical accessors are
`num_sites`, `position`, `neighbors`, and `bonds` — they are all
exported by LatticeCore and work on any `AbstractLattice`. The
names in this file remain for legacy code that spells them out.
"""

"""Position tolerance for duplicate-detection in `build_nearest_neighbor_bonds!`."""
const POSITION_TOLERANCE = 1e-10

"""
    get_positions(data::QuasicrystalData) → Vector{SVector{D, T}}

Return the list of site positions. Equivalent to the LatticeCore
idiom `collect(positions(data))`.
"""
get_positions(data::QuasicrystalData) = data.positions

"""
    get_bonds(data::QuasicrystalData) → Vector{Bond{D, T}}

Return the list of `LatticeCore.Bond` objects populated by
[`build_nearest_neighbor_bonds!`](@ref).
"""
get_bonds(data::QuasicrystalData) = data.bonds

"""
    get_nearest_neighbors(data::QuasicrystalData) → Vector{Vector{Int}}

Return the neighbour adjacency lists. Equivalent to
`[neighbors(data, i) for i in 1:num_sites(data)]`.
"""
get_nearest_neighbors(data::QuasicrystalData) = data.nearest_neighbors

"""
    num_bonds(data::QuasicrystalData) → Int

Number of bonds currently stored. Zero until
[`build_nearest_neighbor_bonds!`](@ref) populates the bond list.
"""
num_bonds(data::QuasicrystalData) = length(data.bonds)

"""
    build_nearest_neighbor_bonds!(data::QuasicrystalData{D, T}; cutoff::Real)

Populate `data.bonds` and `data.nearest_neighbors` with all pairs
of sites whose Euclidean distance is strictly less than `cutoff`
and strictly greater than `POSITION_TOLERANCE`. Existing bonds are
cleared first. Mutates `data.bonds` and the inner vectors of
`data.nearest_neighbors` in place.

Returns `data` for chaining.
"""
function build_nearest_neighbor_bonds!(
    data::QuasicrystalData{D,T}; cutoff::Real
) where {D,T}
    n = num_sites(data)

    # Clear any previously populated connectivity.
    empty!(data.bonds)
    for nb in data.nearest_neighbors
        empty!(nb)
    end

    for i in 1:n
        pos_i = data.positions[i]
        for j in (i + 1):n
            pos_j = data.positions[j]
            bond_vec = pos_j - pos_i
            dist = norm(bond_vec)
            if dist < cutoff && dist > POSITION_TOLERANCE
                push!(data.bonds, Bond{D,T}(i, j, bond_vec, :nearest))
                push!(data.nearest_neighbors[i], j)
                push!(data.nearest_neighbors[j], i)
            end
        end
    end
    return data
end

"""
    build_quasicrystal(type::Type{<:AbstractQuasicrystal};
                       generator::Symbol = :projection,
                       radius = 3.0,
                       generations::Int = 4,
                       n_points::Int = 200)

High-level dispatch wrapper: selects the right generator for a
topology marker.
"""
function build_quasicrystal(
    type::Type{<:AbstractQuasicrystal};
    generator::Symbol=:projection,
    radius=3.0,
    generations::Int=4,
    n_points::Int=200,
)
    if type == PenroseP3
        return if generator == :projection
            generate_penrose_projection(radius)
        else
            generate_penrose_substitution(generations)
        end
    elseif type == AmmannBeenker
        return if generator == :projection
            generate_ammann_beenker_projection(radius)
        else
            generate_ammann_beenker_substitution(generations)
        end
    elseif type == FibonacciLattice
        return if generator == :projection
            generate_fibonacci_projection(n_points)
        else
            generate_fibonacci_substitution(generations)
        end
    else
        error(
            "Unsupported type: $(type). Choose PenroseP3, AmmannBeenker, or FibonacciLattice.",
        )
    end
end
