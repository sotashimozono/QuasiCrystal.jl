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

# ---- Tolerance constants ---------------------------------------------
#
# Centralised numerical tolerances for the substitution / projection
# pipelines. Each call site that previously hardcoded a magic number
# (`1e-10`, `1e-5`, `1e-4`) now references one of these named
# constants so the trade-offs are documented in one place and any
# future tightening / loosening can be done in a single edit.

"""
    VERTEX_MERGE_TOL

Distance tolerance below which two real-space site / vertex positions
are considered the same point. Used by the bond builder
(`build_nearest_neighbor_bonds!`) to skip self-pairs and by the
tile→plaquette KDTree lookup in `_resolve_tile_vertices`.

Value: `1e-10`. This is several orders of magnitude smaller than the
minimum inter-site spacing produced by the shipped projection /
substitution generators (`O(1)` in the natural units of the
construction), and comfortably larger than the round-off accumulated
by a few generations of recursion.
"""
const VERTEX_MERGE_TOL = 1e-10

"""
    POSITION_TOLERANCE

Legacy alias for [`VERTEX_MERGE_TOL`](@ref). Retained for
backwards-compatibility with callers that imported the old name.
"""
const POSITION_TOLERANCE = VERTEX_MERGE_TOL

"""
    SNAP_GRID_EPS

Resolution of the integer-tuple hash key produced by
[`snap_to_grid`](@ref). Distinct sites produced by the substitution
pipeline are well separated at this scale; round-off in deep
recursion is comfortably below it.

Value: `1e-5` — matches the legacy `round(Int, c * 1e5)` resolution.
"""
const SNAP_GRID_EPS = 1e-5

"""
    STAR_DIRECTION_TOL

Tolerance for matching an inflated tile edge to one of the unit
star vectors `e_k = (cos(2πk/n), sin(2πk/n))` in the substitution
inflation rules (Penrose, Ammann–Beenker). Looser than
[`VERTEX_MERGE_TOL`](@ref) because the comparison is between unit
vectors after a long chain of multiplications by `ϕ`, not between
sites of the lattice.

Value: `1e-4`.
"""
const STAR_DIRECTION_TOL = 1e-4

"""
    positions_equal(a, b; tol = VERTEX_MERGE_TOL) -> Bool

Tolerant equality test between two real-space positions.
Returns `true` iff `norm(a - b) < tol`.
"""
@inline function positions_equal(
    a::SVector{D,T}, b::SVector{D,T}; tol::Real=VERTEX_MERGE_TOL
) where {D,T}
    return norm(a - b) < tol
end

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

@inline snap_to_grid(pos::SVector{D,T}) where {D,T} = snap_to_grid(pos, T(SNAP_GRID_EPS))

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
