"""
    analysis/tile_statistics.jl

Tile-level statistics for `QuasicrystalData`. The functions here treat
the tiling itself as the object of interest (rather than the vertex
graph), so they all key off `data.tiles` / `plaquettes(data)`.

Provided API:

- [`tile_counts`](@ref) — per-`TileType` tile multiplicity.
- [`tile_density`](@ref) — per-`TileType` fraction of all tiles.
- [`tile_area`](@ref) — geometric area of a unit-edge tile, dispatched
  on `TileType`.
- [`golden_ratio_check`](@ref) — Penrose-specific check that the
  observed `#FatRhombus / #ThinRhombus` ratio approaches `ϕ`.
- [`tile_perimeter`](@ref) — total length of the union of tile
  boundary edges (each shared edge counted once).

These are pure analysis helpers: they do not mutate `data` (apart from
the existing lazy `parameters[:plaquettes]` / `parameters[:position_index]`
caches used by [`plaquettes`](@ref)).
"""

# ---- Tile counts and densities -----------------------------------

"""
    tile_counts(data::QuasicrystalData) → Dict{TileType, Int}

Return a dictionary mapping each [`TileType`](@ref) singleton present
in `data.tiles` to the number of tiles of that type. Returns an empty
dict for tilings without tiles (e.g. 1D Fibonacci, where the chain is
described purely by site positions).

```jldoctest; setup = :(using QuasiCrystal)
julia> qc = generate_penrose_substitution(3);

julia> counts = tile_counts(qc);

julia> sum(values(counts)) == length(qc.tiles)
true
```
"""
function tile_counts(data::QuasicrystalData)
    counts = Dict{TileType,Int}()
    for tile in data.tiles
        counts[tile.type] = get(counts, tile.type, 0) + 1
    end
    return counts
end

"""
    tile_density(data::QuasicrystalData) → Dict{TileType, Float64}

Per-`TileType` fraction of all tiles. The values sum to `1.0` (within
rounding) for tilings with at least one tile, and the dict is empty
when `data.tiles` is empty.
"""
function tile_density(data::QuasicrystalData)
    counts = tile_counts(data)
    total = length(data.tiles)
    total == 0 && return Dict{TileType,Float64}()
    return Dict{TileType,Float64}(t => c / total for (t, c) in counts)
end

# ---- Per-tile-type area constants --------------------------------

"""
    tile_area(::TileType) → Float64

Geometric area of a unit-edge tile. Per-tile-type constants:

- [`FatRhombus`](@ref):   `sin(2π/5) ≈ 0.9511` (Penrose, 72°/108°)
- [`ThinRhombus`](@ref):  `sin(π/5)  ≈ 0.5878` (Penrose, 36°/144°)
- [`Square`](@ref):       `1.0` (Ammann–Beenker)
- [`RhombusAB`](@ref):    `sin(π/4)  ≈ 0.7071` (Ammann–Beenker, 45°/135°)

The values assume the canonical unit edge length used by the
`generate_*` functions; tiles produced by inflation algorithms that
rescale lengths should be normalised before applying these constants.
"""
tile_area(::FatRhombus) = sin(2π / 5)
tile_area(::ThinRhombus) = sin(π / 5)
tile_area(::Square) = 1.0
tile_area(::RhombusAB) = sin(π / 4)

# ---- Penrose golden-ratio check ----------------------------------

"""
    golden_ratio_check(data::QuasicrystalData; tol = 0.05)
        → NamedTuple{(:observed, :expected, :ok)}

For Penrose P3 tilings, the asymptotic ratio
`#FatRhombus / #ThinRhombus → ϕ`. This helper measures the observed
ratio on `data` and compares it to `ϕ` with relative tolerance `tol`.

The returned NamedTuple has fields

- `observed::Float64` — `#FatRhombus / #ThinRhombus` on `data`
- `expected::Float64` — the golden ratio [`ϕ`](@ref)
- `ok::Bool` — `abs(observed - expected) / expected ≤ tol`

Throws `ArgumentError` when `data` has no `ThinRhombus` tiles, since
the ratio is undefined in that case.
"""
function golden_ratio_check(data::QuasicrystalData; tol::Real=0.05)
    counts = tile_counts(data)
    n_fat = get(counts, FatRhombus(), 0)
    n_thin = get(counts, ThinRhombus(), 0)
    n_thin == 0 && throw(
        ArgumentError(
            "golden_ratio_check needs at least one ThinRhombus tile; got $(counts)"
        ),
    )
    observed = n_fat / n_thin
    expected = float(ϕ)
    ok = abs(observed - expected) / expected ≤ tol
    return (observed=observed, expected=expected, ok=ok)
end

# ---- Tile perimeter (deduplicated boundary edges) ----------------

"""
    tile_perimeter(data::QuasicrystalData) → Float64

Total length of the union of tile boundary edges. Edges shared by two
tiles (interior bonds of the tiling) are counted **once**, so the
result equals `sum_of_distinct_edge_lengths`.

Returns `0.0` for tilings without tiles. The implementation reuses
[`plaquettes`](@ref)' integer vertex indices to deduplicate edges
exactly (no float-tolerance heuristic) and reads the geometric edge
length from `data.positions`.
"""
function tile_perimeter(data::QuasicrystalData{D,T}) where {D,T}
    isempty(data.tiles) && return zero(float(T))
    plaqs = plaquettes(data)
    seen = Set{Tuple{Int,Int}}()
    total = zero(float(T))
    for p in plaqs
        vs = p.vertices
        n = length(vs)
        for k in 1:n
            a = vs[k]
            b = vs[mod1(k + 1, n)]
            edge = a < b ? (a, b) : (b, a)
            if !(edge in seen)
                push!(seen, edge)
                total += norm(data.positions[edge[1]] - data.positions[edge[2]])
            end
        end
    end
    return total
end
