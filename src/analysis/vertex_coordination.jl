"""
Vertex coordination and vertex-type analysis for `QuasicrystalData`.

`coordination` works on any `QuasicrystalData` and reads from the
populated `data.nearest_neighbors` adjacency lists, so the bond list
must already have been built via
[`build_nearest_neighbor_bonds!`](@ref).

`vertex_type` and `vertex_type_counts` are Penrose-specific (the
canonical 8 vertex configurations Sun / Star / Ace / Deuce / Jack /
Queen / King are well-defined only on the Penrose P3 tiling). The
classifier delegates to the existing
[`vertex_configuration`](@ref) helper, which reads off interior
angles around the vertex from `plaquettes(data)`. Vertices not
matching any of the 8 standard signatures (e.g. boundary vertices on
a finite patch) are reported as `:Other`.
"""

# ---- coordination -------------------------------------------------

"""
    coordination(data::QuasicrystalData, i::Int) → Int

Number of bonds incident to vertex `i`. Reads from
`data.nearest_neighbors`, so [`build_nearest_neighbor_bonds!`](@ref)
must have populated the adjacency lists first; otherwise the result
is `0`.
"""
function coordination(data::QuasicrystalData, i::Int)
    1 ≤ i ≤ num_sites(data) || throw(BoundsError(data.positions, i))
    return length(data.nearest_neighbors[i])
end

"""
    coordination(data::QuasicrystalData) → Vector{Int}

Coordination number of every vertex, as a `Vector{Int}` of length
`num_sites(data)`.
"""
function coordination(data::QuasicrystalData)
    return [length(nb) for nb in data.nearest_neighbors]
end

# ---- vertex_type (Penrose) ----------------------------------------

"""
    vertex_type(data::QuasicrystalData{2, Float64, PenroseP3}, i::Int) → Symbol

Classify the local environment of vertex `i` on a Penrose P3 tiling
into one of the 8 canonical vertex configurations
(`:Sun`, `:Star`, `:Ace`, `:Deuce`, `:Jack`, `:Queen`, `:King`) or
`:Other` if the angle signature does not match a standard interior
configuration (typically a boundary vertex of a finite patch).

Implementation: delegates to the existing
[`vertex_configuration`](@ref) helper for the angle-signature lookup
and re-maps any non-standard return (`:Boundary`, `Symbol("V_…")`)
to `:Other`. Hence the result is always one of the eight symbols
`{:Sun, :Star, :Ace, :Deuce, :Jack, :Queen, :King, :Other}`.

The classifier requires `data.tiles` to be populated — the
projection-method generator currently does not emit tiles, so call
[`generate_penrose_substitution`](@ref) (or build tiles separately)
before using this API.
"""
function vertex_type(data::QuasicrystalData{2,Float64,PenroseP3}, i::Int)
    1 ≤ i ≤ num_sites(data) || throw(BoundsError(data.positions, i))
    cfg = vertex_configuration(data, i)
    if cfg in (:Sun, :Star, :Ace, :Deuce, :Jack, :Queen, :King)
        return cfg
    else
        return :Other
    end
end

"""
    vertex_type_counts(data::QuasicrystalData{2, Float64, PenroseP3})
        → Dict{Symbol, Int}

Histogram of `vertex_type(data, i)` over `i = 1 … num_sites(data)`.
The sum of values equals `num_sites(data)`. Keys are a subset of
`{:Sun, :Star, :Ace, :Deuce, :Jack, :Queen, :King, :Other}`; types
that do not appear in the patch are simply absent from the
returned dictionary.
"""
function vertex_type_counts(data::QuasicrystalData{2,Float64,PenroseP3})
    counts = Dict{Symbol,Int}()
    for i in 1:num_sites(data)
        t = vertex_type(data, i)
        counts[t] = get(counts, t, 0) + 1
    end
    return counts
end
