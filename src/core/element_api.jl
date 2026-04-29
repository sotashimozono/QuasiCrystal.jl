"""
    element_api.jl

QuasiCrystal-side implementation of the LatticeCore element-center
and incidence API (`plaquettes`, `num_plaquettes`, `bond_type`,
`neighbors(..., shell=k)`, and the full element-center surface).

QuasiCrystal's own `Tile` type already describes the plaquettes of a
cut-and-project tiling — the per-tile `vertices::Vector{SVector{D,T}}`
and `center::SVector{D,T}` are exactly what
`LatticeCore.Plaquette{D,T}` wants, except that the LatticeCore type
identifies vertices by **integer site indices** rather than positions.
Here we lazily convert `data.tiles` into a `Vector{Plaquette{D,T}}`
on first access and stash the result in `data.parameters[:plaquettes]`
so subsequent calls are O(1).

The conversion looks up each tile vertex in `data.positions` via a
tolerant linear search — tile vertices are always a subset of the
lattice positions because both come from the same projection
pipeline, so the search always succeeds.
"""

# ---- Plaquettes: tile → Plaquette{D, T} conversion ---------------

"""
    plaquettes(data::QuasicrystalData{D, T}) → Vector{Plaquette{D, T}}

Return the quasicrystal's plaquette list, derived from `data.tiles`
on first access and cached in `data.parameters[:plaquettes]`. The
plaquette type tag is `Symbol("tile_type_\$n")` where `n` is the
integer `Tile.type` — downstream code that needs a richer tag can
override the conversion after construction.
"""
function LatticeCore.plaquettes(data::QuasicrystalData{D,T}) where {D,T}
    _ensure_plaquettes!(data)
    return data.parameters[:plaquettes]::Vector{Plaquette{D,T}}
end

"""
    _ensure_plaquettes!(data::QuasicrystalData)

Populate `data.parameters[:plaquettes]` on first access. The stash
lives inside the already-mutable `parameters::Dict{Symbol,Any}` bag,
so the struct itself stays immutable.
"""
function _ensure_plaquettes!(data::QuasicrystalData{D,T}) where {D,T}
    if haskey(data.parameters, :plaquettes)
        return nothing
    end
    data.parameters[:plaquettes] = _materialise_plaquettes(data)
    return nothing
end

function _materialise_plaquettes(data::QuasicrystalData{D,T}) where {D,T}
    out = Plaquette{D,T}[]
    isempty(data.tiles) && return out
    pidx = _ensure_position_index!(data)
    for tile in data.tiles
        vertex_ids = _resolve_tile_vertices(data, tile, pidx)
        push!(out, Plaquette{D,T}(vertex_ids, tile.center, Symbol("tile_type_", tile.type)))
    end
    return out
end

"""
    _ensure_position_index!(data) → PositionIndex

Lazily build and cache a KDTree over `data.positions` in
`data.parameters[:position_index]`. The cache is keyed on the
positions vector identity so it survives multiple plaquette / lookup
calls without rebuilding.
"""
function _ensure_position_index!(data::QuasicrystalData{D,T}) where {D,T}
    cached = get(data.parameters, :position_index, nothing)
    if cached !== nothing
        return cached::PositionIndex
    end
    pidx = build_position_index(data.positions)
    data.parameters[:position_index] = pidx
    return pidx
end

"""
    _resolve_tile_vertices(data, tile, pidx) → Vector{Int}

Map each of `tile.vertices` (a position in real space) back to its
integer index in `data.positions` via a tolerant KDTree query.
Throws if the tile vertex doesn't match any lattice site.
"""
function _resolve_tile_vertices(
    data::QuasicrystalData{D,T}, tile::Tile{D,T}, pidx::PositionIndex
) where {D,T}
    out = Int[]
    sizehint!(out, length(tile.vertices))
    for tv in tile.vertices
        idx = find_position_index(pidx, tv, POSITION_TOLERANCE)
        idx == 0 && throw(
            ArgumentError(
                "tile vertex at $(tv) does not match any site position on $(typeof(data).name.name)",
            ),
        )
        push!(out, idx)
    end
    return out
end

"""
    num_plaquettes(data::QuasicrystalData) → Int

Number of plaquettes (tiles) on `data`. O(1) after first access.
"""
num_plaquettes(data::QuasicrystalData) = length(plaquettes(data))

# ---- Element-center API specialisations -------------------------

LatticeCore.num_elements(data::QuasicrystalData, ::VertexCenter) = num_sites(data)
LatticeCore.num_elements(data::QuasicrystalData, ::BondCenter) = length(data.bonds)
function LatticeCore.num_elements(data::QuasicrystalData, ::PlaquetteCenter)
    return length(plaquettes(data))
end

LatticeCore.elements(data::QuasicrystalData, ::VertexCenter) = 1:num_sites(data)
LatticeCore.elements(data::QuasicrystalData, ::BondCenter) = data.bonds
LatticeCore.elements(data::QuasicrystalData, ::PlaquetteCenter) = plaquettes(data)

function LatticeCore.element_position(data::QuasicrystalData, ::BondCenter, i::Int)
    return bond_center(data, data.bonds[i])
end

function LatticeCore.element_position(data::QuasicrystalData, ::PlaquetteCenter, i::Int)
    return plaquettes(data)[i].center
end

# ---- Multi-shell neighbours (distance-based) ---------------------

"""
    neighbors(data::QuasicrystalData, i::Int; shell::Int = nothing)

Without `shell`: same as the existing adjacency list
(`data.nearest_neighbors[i]`).

With `shell = k ≥ 1`: the `k`-th geometric neighbour shell around
site `i`, ranked by Euclidean distance with a float tolerance.
Unlike Lattice2D's periodic multi-shell walk, quasicrystals are
genuinely aperiodic so the search is over every other site — this is
O(N) per call and is intended for inspection and observables, not
MC hot paths.
"""
function LatticeCore.neighbors(
    data::QuasicrystalData, i::Int; shell::Union{Nothing,Int}=nothing
)
    if shell === nothing
        return data.nearest_neighbors[i]
    else
        shell ≥ 1 || throw(ArgumentError("shell must be ≥ 1, got $shell"))
        return _neighbors_by_shell(data, i, shell)
    end
end

function _neighbors_by_shell(data::QuasicrystalData{D,T}, i::Int, k::Int) where {D,T}
    N = num_sites(data)
    p_i = position(data, i)
    dist_map = Dict{Int,T}()
    for j in 1:N
        j == i && continue
        d = norm(position(data, j) - p_i)
        dist_map[j] = d
    end

    isempty(dist_map) && return Int[]

    sorted_d = sort(collect(values(dist_map)))
    shells = T[sorted_d[1]]
    for d in sorted_d
        ref = shells[end]
        if d - ref > 10 * eps(T) * max(one(T), ref)
            push!(shells, d)
        end
    end

    k > length(shells) && return Int[]
    target = shells[k]
    tol = 10 * eps(T) * max(one(T), target)
    return sort!([j for (j, d) in dist_map if abs(d - target) ≤ tol])
end

# ---- bond_type query ---------------------------------------------

"""
    bond_type(data::QuasicrystalData, i::Int, j::Int) → Symbol

Return the bond-type tag of the edge connecting sites `i` and `j`.
Currently the `build_nearest_neighbor_bonds!` populator emits only
`:nearest` bonds, so the query returns `:nearest` on a successful
lookup; more exotic bond-type populators can override without
changing the API.

Throws `ArgumentError` if there is no bond between `i` and `j` on
`data`.
"""
function bond_type(data::QuasicrystalData, i::Int, j::Int)
    for b in data.bonds
        if (b.i == i && b.j == j) || (b.i == j && b.j == i)
            return b.type
        end
    end
    throw(ArgumentError("no bond between sites $i and $j on $(typeof(data).name.name)"))
end
