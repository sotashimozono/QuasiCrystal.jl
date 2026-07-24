"""
Backwards-compat accessors for `QuasicrystalData` plus the
distance-based nearest-neighbour builder.

After the LatticeCore migration the canonical accessors are
`num_sites`, `position`, `neighbors`, and `bonds` — they are all
exported by LatticeCore and work on any `AbstractLattice`. The
names in this file remain for legacy code that spells them out.
"""

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
and strictly greater than `VERTEX_MERGE_TOL`. Existing bonds are
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
            if dist < cutoff && dist > VERTEX_MERGE_TOL
                push!(data.bonds, Bond{D,T}(i, j, bond_vec, :nearest))
                push!(data.nearest_neighbors[i], j)
                push!(data.nearest_neighbors[j], i)
            end
        end
    end
    return data
end

"""
    build_chain_bonds!(data::QuasicrystalData{1, T}) → data

Populate `data.bonds` / `data.nearest_neighbors` with the edges of the
1D chain: each site is bonded to its neighbours in position order (the
tile edges of the quasiperiodic chain). Existing bonds are cleared
first.

Unlike [`build_nearest_neighbor_bonds!`](@ref), whose absolute distance
`cutoff` selects a *fixed length scale*, consecutive-adjacency is
**scale invariant**: it depends only on the ordering of the sites, not
on their absolute spacing. It is therefore covariant under inflation —
`build_chain_bonds!` on `materialize(inflate(inf))` yields the same
bond graph as on `materialize(inf)`, with every bond vector scaled by
the inflation factor λ — whereas a fixed `cutoff` would silently drop or
add bonds once the chain is rescaled.
"""
function build_chain_bonds!(data::QuasicrystalData{1,T}) where {T}
    empty!(data.bonds)
    for nb in data.nearest_neighbors
        empty!(nb)
    end

    order = sortperm(data.positions; by=p -> p[1])
    for k in 1:(length(order) - 1)
        i, j = order[k], order[k + 1]
        bond_vec = data.positions[j] - data.positions[i]
        push!(data.bonds, Bond{1,T}(i, j, bond_vec, :nearest))
        push!(data.nearest_neighbors[i], j)
        push!(data.nearest_neighbors[j], i)
    end
    return data
end

"""
    build_tile_bonds!(data::QuasicrystalData{D, T}) → data

Populate `data.bonds` / `data.nearest_neighbors` with the edges of the
tiling: every boundary edge of every tile in `data.tiles`, deduplicated
across shared edges. Existing bonds are cleared first. The tiling's
vertices are resolved to site indices through the same path as
`plaquettes(data)`.

This is the 2D counterpart of [`build_chain_bonds!`](@ref): the bond set
is the *tiling-intrinsic* adjacency (which sites share a tile edge), so
it carries **no length-scale parameter** — unlike
[`build_nearest_neighbor_bonds!`](@ref), whose absolute `cutoff` fixes a
scale and can catch short tile diagonals or miss edges once the patch is
rescaled. For the Penrose and Ammann–Beenker rhombus tilings every such
edge is a unit edge, so the graph is exactly the unit-edge network at any
inflation generation.

Throws if `data` has no tiling (e.g. a projection-generated patch, which
records no tiles): use [`build_chain_bonds!`](@ref) for a 1D chain or
[`build_nearest_neighbor_bonds!`](@ref) otherwise.
"""
function build_tile_bonds!(data::QuasicrystalData{D,T}) where {D,T}
    empty!(data.bonds)
    for nb in data.nearest_neighbors
        empty!(nb)
    end
    isempty(data.tiles) && throw(
        ArgumentError(
            "build_tile_bonds! needs a tiling, but this patch records no tiles " *
            "(projection-generated?). Use `build_chain_bonds!` (1D chain) or " *
            "`build_nearest_neighbor_bonds!`.",
        ),
    )

    seen = Set{Tuple{Int,Int}}()
    for p in plaquettes(data)
        vs = p.vertices
        n = length(vs)
        for k in 1:n
            i, j = vs[k], vs[mod1(k + 1, n)]
            i == j && continue
            key = (min(i, j), max(i, j))
            key in seen && continue
            push!(seen, key)
            bond_vec = data.positions[j] - data.positions[i]
            push!(data.bonds, Bond{D,T}(i, j, bond_vec, :nearest))
            push!(data.nearest_neighbors[i], j)
            push!(data.nearest_neighbors[j], i)
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
