"""
    numerics.jl

Shared numerical utilities for stable floating-point dedup and lookup
in the substitution / projection pipelines.

The two routines here are used by the tile generators (Penrose,
Ammann–Beenker) and by the tile→plaquette conversion in
`element_api.jl`.

* `snap_to_grid(pos, eps)` produces a stable `NTuple{D, Int}` hash key
  by snapping each coordinate to a multiple of `eps`. Replaces ad-hoc
  `round(Int, c[i] * 1e5)` patterns that were sensitive to drift in
  the last digits.
* `PositionIndex` wraps a `KDTree` (from NearestNeighbors.jl) over the
  unique-vertex list so that the tile→plaquette conversion can do
  O(log N) tolerant lookups instead of the previous O(N) linear scan.
"""

using NearestNeighbors

# ---- snap_to_grid -----------------------------------------------------

"""
    snap_to_grid(pos::SVector{D,T}, eps::Real) -> NTuple{D, Int}

Snap each component of `pos` to the nearest multiple of `eps` and
return the resulting integer tuple. Used as a stable hash key for
deduplicating tile centres / vertices that should be equal up to
floating-point round-off.

`eps` should be small relative to the minimum spacing of distinct
sites but comfortably larger than the accumulated round-off in the
substitution recursion. The defaults in callers use `1e-5` (matches
the legacy `round(Int, c*1e5)` resolution).
"""
@inline function snap_to_grid(pos::SVector{D,T}, eps::Real) where {D,T}
    inv_eps = inv(T(eps))
    return ntuple(i -> round(Int, pos[i] * inv_eps), D)
end

@inline snap_to_grid(pos::SVector{D,T}) where {D,T} = snap_to_grid(pos, T(1e-5))

# ---- PositionIndex ---------------------------------------------------

"""
    PositionIndex{D,T}

KDTree-backed accelerator for tolerant lookup of a position in a
`Vector{SVector{D,T}}`. Build once via `build_position_index` and
query with `find_position_index`.
"""
struct PositionIndex{D,T,Tree}
    tree::Tree
    npts::Int
end

"""
    build_position_index(positions::Vector{SVector{D,T}}) -> PositionIndex

Construct a KDTree over `positions`. Each query is O(log N).
"""
function build_position_index(positions::Vector{SVector{D,T}}) where {D,T}
    if isempty(positions)
        # KDTree requires at least one point; return a sentinel
        return PositionIndex{D,T,Nothing}(nothing, 0)
    end
    # NearestNeighbors accepts a D×N matrix; build one explicitly.
    N = length(positions)
    data = Matrix{T}(undef, D, N)
    @inbounds for j in 1:N
        for k in 1:D
            data[k, j] = positions[j][k]
        end
    end
    tree = KDTree(data)
    return PositionIndex{D,T,typeof(tree)}(tree, length(positions))
end

"""
    find_position_index(idx::PositionIndex, target, tol) -> Int

Return the 1-based index of the closest point to `target` if its
distance is `< tol`, else `0`.
"""
function find_position_index(
    idx::PositionIndex{D,T,Nothing}, target::SVector{D,T}, tol::Real
) where {D,T}
    return 0
end

function find_position_index(
    idx::PositionIndex{D,T}, target::SVector{D,T}, tol::Real
) where {D,T}
    i, d = nn(idx.tree, target)
    return d < tol ? i : 0
end
