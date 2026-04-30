"""
Penrose P3 (rhombus) tiling implementation: cut-and-project
quasicrystal with 10-fold (5-fold) rotational symmetry from a 5D
hypercubic host lattice, plus a Robinson-triangle inflation
substitution generator.
"""

"""
    PenroseP3 <: AbstractQuasicrystal{2}

Topology marker for the Penrose P3 rhombus tiling. Pass to
[`build_quasicrystal`](@ref) or call
[`generate_penrose_projection`](@ref) /
[`generate_penrose_substitution`](@ref) directly.
"""
struct PenroseP3 <: AbstractQuasicrystal{2} end

"""
    generate_penrose_projection(radius::Real;
                                method::ProjectionMethod = ProjectionMethod())
        → QuasicrystalData{2, Float64}

Generate a Penrose P3 point set by projecting integer lattice
points of `Z^5` onto the 2D physical subspace of the cut-and-
project construction. The perpendicular acceptance window is a
5-cube of half-width `window_size`.
"""
function generate_penrose_projection(
    radius::Real; method::ProjectionMethod=ProjectionMethod()
)
    theta = 2π / 5

    # Parallel projection: 5D → 2D (physical)
    E_par = zeros(5, 2)
    for i in 1:5
        E_par[i, 1] = cos((i - 1) * theta)
        E_par[i, 2] = sin((i - 1) * theta)
    end

    # Perpendicular projection: 5D → 3D (internal, acceptance window)
    E_perp = zeros(5, 3)
    for i in 1:5
        E_perp[i, 1] = cos(2 * (i - 1) * theta)
        E_perp[i, 2] = sin(2 * (i - 1) * theta)
        E_perp[i, 3] = cos(3 * (i - 1) * theta)
    end

    window_size = 0.5
    n_max = ceil(Int, radius * 1.5)

    # Pre-typed SVector buffer: avoids the intermediate `Vector{Float64}[]`
    # and the second pass that converted each entry to `SVector{2,Float64}`.
    positions = SVector{2,Float64}[]
    P = SMatrix{2,5,Float64}(E_par')
    Q = SMatrix{3,5,Float64}(E_perp')

    for n1 in (-n_max):n_max,
        n2 in (-n_max):n_max, n3 in (-n_max):n_max, n4 in (-n_max):n_max,
        n5 in (-n_max):n_max

        lp = SVector{5,Float64}(float(n1), float(n2), float(n3), float(n4), float(n5))
        pos_par = P * lp
        norm(pos_par) > radius && continue

        pos_perp = Q * lp
        if abs(pos_perp[1]) <= window_size &&
            abs(pos_perp[2]) <= window_size &&
            abs(pos_perp[3]) <= window_size
            push!(positions, pos_par)
        end
    end
    tiles = Tile{2,Float64}[]   # Tile construction is a future task.

    params = Dict{Symbol,Any}(
        :radius => radius,
        :n_max => n_max,
        :window_size => window_size,
        :window_shape => :box_3d,
        :n_vertices => length(positions),
    )
    return QuasicrystalData{2,Float64}(PenroseP3(), positions, tiles, method, params)
end

# ---- Robinson-triangle representation ---------------------------------
#
# A `RobTri` encodes one of the two Penrose Robinson triangles:
#
#   - `:red`  = acute "golden triangle" = half of a thin rhombus.
#               Sides (1, 1, 1/φ); apex angle 36° at vertex `A`;
#               base angles 72° at `B` and `C`.
#   - `:blue` = obtuse "golden gnomon" = half of a fat rhombus.
#               Sides (1, 1, φ); apex angle 108° at vertex `A`;
#               base angles 36° at `B` and `C`.
#
# `A` is always the apex; `B` and `C` are the base. The (B, C)
# ordering encodes chirality, which the subdivision rules respect
# so that pairs of half-tiles can be reassembled into a rhombus by
# matching the appropriate edge.

"""
    RobTri

A single Robinson triangle used by the Penrose P3 substitution
pipeline.

# Fields
- `kind::Symbol` — `:red` (acute / half-thin) or `:blue` (obtuse /
  half-fat).
- `A::SVector{2,Float64}` — apex vertex.
- `B::SVector{2,Float64}` — first base vertex.
- `C::SVector{2,Float64}` — second base vertex.

`A`, `B`, `C` are oriented so that `(B, C)` traces the base from
`B` to `C` consistently within an inflation generation; sub-tiles
emitted by [`subdivide_red`](@ref) / [`subdivide_blue`](@ref)
preserve this orientation so that pairs can be re-merged into
rhombi via [`merge_robinson_to_rhombi`](@ref).
"""
struct RobTri
    kind::Symbol
    A::SVector{2,Float64}
    B::SVector{2,Float64}
    C::SVector{2,Float64}
end

"""
    subdivide_red(t::RobTri) → Vector{RobTri}

Apply the Penrose deflation rule to an acute Robinson triangle
(`t.kind == :red`). Each red triangle decomposes into one smaller
red plus one smaller blue, with all sub-edges scaled by `1/ϕ`:

    red(A, B, C) → red(C, P, B), blue(P, C, A)

where `P = A + (B − A) / ϕ`.
"""
function subdivide_red(t::RobTri)
    A, B, C = t.A, t.B, t.C
    P = A + (B - A) / ϕ
    return RobTri[RobTri(:red, C, P, B), RobTri(:blue, P, C, A)]
end

"""
    subdivide_blue(t::RobTri) → Vector{RobTri}

Apply the Penrose deflation rule to an obtuse Robinson triangle
(`t.kind == :blue`). Each blue triangle decomposes into one red
plus two blues, with all sub-edges scaled by `1/ϕ`:

    blue(A, B, C) → blue(Q, B, R), blue(R, C, A), red(R, A, Q)

where `Q = B + (A − B) / ϕ` and `R = B + (C − B) / ϕ`.

Combined with [`subdivide_red`](@ref), the rhombus-level counts
satisfy `(#fat, #thin) ↦ (2·#fat + #thin, #fat + #thin)` per
generation, so `#fat / #thin → ϕ` asymptotically.

Each sub-triangle is labelled with its apex as the first argument
(`A_new`); the `(B_new, C_new)` order preserves the parent's
chirality so the rule can be re-applied recursively without an
explicit mirror correction.
"""
function subdivide_blue(t::RobTri)
    A, B, C = t.A, t.B, t.C
    Q = B + (A - B) / ϕ
    R = B + (C - B) / ϕ
    return RobTri[RobTri(:blue, Q, B, R), RobTri(:blue, R, C, A), RobTri(:red, R, A, Q)]
end

"""
    subdivide(t::RobTri) → Vector{RobTri}

Dispatch on `t.kind`.
"""
subdivide(t::RobTri) = t.kind === :red ? subdivide_red(t) : subdivide_blue(t)

"""
    deflate_robinson(tris::Vector{RobTri}) → Vector{RobTri}

Apply one generation of the Robinson triangle substitution rule to
every triangle in `tris`. The output edge length is the input edge
length divided by `ϕ`.
"""
function deflate_robinson(tris::Vector{RobTri})
    out = Vector{RobTri}(undef, 0)
    sizehint!(out, 3 * length(tris))
    for t in tris
        append!(out, subdivide(t))
    end
    return out
end

# ---- Initial seeds ----------------------------------------------------

"""
    sun_seed_blue() → Vector{RobTri}

Initial "Sun" patch for the Penrose substitution: 10 obtuse
(`:blue`) Robinson triangles whose 36° base vertices share the
origin, forming a wheel of 5 fat rhombi (each fat rhombus = 2 blue
half-tiles glued along the long diagonal of length `ϕ`). The
apex (108°) of each blue triangle lies on a circle of radius `1/ϕ`
... actually on the unit circle; see implementation.
"""
function sun_seed_blue()
    # We construct the Sun = 5 fat rhombi sharing the origin at one
    # 72° corner each. For fat rhombus k (k = 0..4), the 72° corner
    # at the origin spans the angular sector
    #   [k·72° − 36°, k·72° + 36°]  (in 1-based ϕ-rotation about z).
    # Adjacent edges go in directions (k·72° ± 36°), each unit length;
    # the two 108° corners sit at those unit-vector endpoints, and
    # the opposite 72° corner is at their sum (2·cos(36°) = ϕ along
    # the bisector).  The fat rhombus is then split along its long
    # diagonal (origin → opposite-72° corner) into two blue
    # triangles, with apex (108°) at each of the two 108° corners.
    #
    # Both halves of one rhombus get the **same** `(B, C)` labelling
    # (`B = O`, `C = Q`) so that the deflation rule
    # [`subdivide_blue`](@ref) places its sub-tile cuts at the same
    # spatial point on both sides of the shared diagonal — without
    # this, the two parents subdivide in mismatched ways and
    # half-tiles fail to pair into rhombi after one generation.
    # Geometrically the two halves of a fat rhombus are mirror
    # images of each other; using the same `(B, C)` orientation for
    # both means the apex parity (`B → C` traversed CCW vs CW around
    # `A`) flips between the two halves, and the deflation rule is
    # designed to be parity-agnostic.
    tris = RobTri[]
    sizehint!(tris, 10)
    for k in 0:4
        θ_minus = (k * 2π / 5) - π / 5
        θ_plus = (k * 2π / 5) + π / 5
        P_minus = SVector{2,Float64}(cos(θ_minus), sin(θ_minus))
        P_plus = SVector{2,Float64}(cos(θ_plus), sin(θ_plus))
        O = SVector{2,Float64}(0.0, 0.0)
        Q = P_minus + P_plus  # opposite 72° corner, |Q| = ϕ

        # Blue triangle 1: apex at P_plus (108°, "above" the diagonal),
        # base (B = O, C = Q).
        # Blue triangle 2: apex at P_minus (108°, "below" the diagonal),
        # base (B = O, C = Q) — same labelling, mirror parity.
        push!(tris, RobTri(:blue, P_plus, O, Q))
        push!(tris, RobTri(:blue, P_minus, O, Q))
    end
    return tris
end

# ---- Half-tile → rhombus reassembly ----------------------------------

"""
    merge_robinson_to_rhombi(tris::Vector{RobTri}) → Vector{Tile{2,Float64}}

Reassemble a list of Robinson half-tiles into Penrose rhombi.

Each pair of half-tiles that share their long base (for `:blue` =
half-fat: the side of length `ϕ`·s) merges into a fat rhombus;
each pair of `:red` half-tiles that share their short base (length
`s/ϕ`) merges into a thin rhombus, where `s` is the current
edge-length scale of the tiling.

The implementation buckets half-tiles by their (rounded) shared
edge midpoint and pairs them up. Unpaired half-tiles (typically on
the patch boundary) are discarded — they correspond to truncated
rhombi at the convex hull of the patch.
"""
function merge_robinson_to_rhombi(tris::Vector{RobTri})
    isempty(tris) && return Tile{2,Float64}[]

    # Bucket half-tiles by their "outer base" midpoint:
    # - blue (half-fat): the long base is |B − C| = ϕ·s; this is
    #   the cut of the original fat rhombus' long diagonal. Two
    #   blues that came from the same fat rhombus share this base.
    # - red (half-thin): the short base is |B − C| = s/ϕ; two reds
    #   from the same thin rhombus share this base.
    # We use the midpoint of the (B, C) base as the bucket key
    # (snapped to a stable integer grid) and merge buckets of size
    # 2.

    # Use a tighter grid eps that scales with s (the half-tile edge
    # length): SNAP_GRID_EPS = 1e-5 is plenty for s ≥ 1e-3 and
    # comfortably above the deflation round-off (a few × ulp(s)
    # per generation).
    eps_snap = SNAP_GRID_EPS

    blue_groups = Dict{NTuple{2,Int},Vector{Int}}()
    red_groups = Dict{NTuple{2,Int},Vector{Int}}()
    for (i, t) in enumerate(tris)
        midpoint = (t.B + t.C) / 2
        key = snap_to_grid(midpoint, eps_snap)
        if t.kind === :blue
            push!(get!(() -> Int[], blue_groups, key), i)
        else
            push!(get!(() -> Int[], red_groups, key), i)
        end
    end

    out = Tile{2,Float64}[]
    sizehint!(out, length(tris) ÷ 2)

    # --- Fat rhombi from blue pairs ---
    for (_, idxs) in blue_groups
        length(idxs) == 2 || continue
        t1, t2 = tris[idxs[1]], tris[idxs[2]]
        # The two apices A1, A2 are the two 108° corners of the
        # rhombus; the base endpoints B, C (shared between the
        # halves up to permutation) are the two 72° corners.
        # Rhombus vertices (CCW): B, A1, C, A2.
        A1, A2 = t1.A, t2.A
        # Make sure (B, C) corresponds to t1's base.
        B, C = t1.B, t1.C
        # CCW order check: the rhombus ordering (B, A1, C, A2) is
        # convex; pick the orientation such that A1 lies on the
        # left of B→C.
        if _signed_orientation(B, A1, C) < 0
            A1, A2 = A2, A1
        end
        v1, v2, v3, v4 = B, A1, C, A2
        center = (B + C) / 2
        push!(out, Tile{2,Float64}([v1, v2, v3, v4], FatRhombus(), center))
    end

    # --- Thin rhombi from red pairs ---
    for (_, idxs) in red_groups
        length(idxs) == 2 || continue
        t1, t2 = tris[idxs[1]], tris[idxs[2]]
        A1, A2 = t1.A, t2.A  # the two 72° corners (apex of acute = 36°? wait)
        # NB: red is acute, apex angle 36°; the two reds that share
        # the short base BC together form a thin rhombus. The
        # apices A1, A2 become the two 36° corners of the thin
        # rhombus, while B and C become the two 144° corners.
        B, C = t1.B, t1.C
        if _signed_orientation(B, A1, C) < 0
            A1, A2 = A2, A1
        end
        v1, v2, v3, v4 = B, A1, C, A2
        center = (B + C) / 2
        push!(out, Tile{2,Float64}([v1, v2, v3, v4], ThinRhombus(), center))
    end

    return out
end

"""
    _signed_orientation(p, q, r) → Float64

Sign of the cross product `(q − p) × (r − p)`. Positive when
`(p, q, r)` is counter-clockwise.
"""
@inline function _signed_orientation(
    p::SVector{2,Float64}, q::SVector{2,Float64}, r::SVector{2,Float64}
)
    return (q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1])
end

# ---- Public substitution generator -----------------------------------

"""
    generate_penrose_substitution(generations::Int;
                                  method::SubstitutionMethod = SubstitutionMethod())
        → QuasicrystalData{2, Float64}

Generate a Penrose P3 patch by `generations` rounds of Robinson
triangle deflation starting from a 5-fold "Sun" seed (5 fat
rhombi = 10 obtuse half-tiles).

Each generation applies the standard P3 substitution rule
(see [`subdivide_red`](@ref) / [`subdivide_blue`](@ref)):

- each acute Robinson triangle (red, half-thin) → 1 red + 1 blue;
- each obtuse Robinson triangle (blue, half-fat) → 1 red + 2 blue.

Edge lengths are divided by `ϕ` per generation. The patch is
re-scaled at the end so the final rhombus edge has unit length, so
the returned tiling has edge 1 regardless of `generations`.

After deflation, half-tiles are paired up into rhombi via
[`merge_robinson_to_rhombi`](@ref). The asymptotic tile-count
ratio satisfies `#FatRhombus / #ThinRhombus → ϕ`; see
[`golden_ratio_check`](@ref).
"""
function generate_penrose_substitution(
    generations::Int; method::SubstitutionMethod=SubstitutionMethod()
)
    generations ≥ 0 || throw(ArgumentError("generations must be ≥ 0, got $generations"))

    tris = sun_seed_blue()
    for _ in 1:generations
        tris = deflate_robinson(tris)
    end

    # Rescale so that the final rhombus edge is 1.
    s = isempty(tris) ? 1.0 : norm(tris[1].A - tris[1].B)
    if !(isapprox(s, 1.0))
        scale = 1.0 / s
        tris = RobTri[RobTri(t.kind, t.A * scale, t.B * scale, t.C * scale) for t in tris]
    end

    tiles_raw = merge_robinson_to_rhombi(tris)

    # Final dedup on tile centres (typically a no-op after the
    # half-tile pairing, but cheap insurance against accumulated
    # round-off near the patch boundary).
    tile_dict = Dict{NTuple{2,Int},Tile{2,Float64}}()
    for tile in tiles_raw
        key = snap_to_grid(tile.center, SNAP_GRID_EPS)
        get!(tile_dict, key, tile)
    end
    tiles = collect(values(tile_dict))

    # Collect unique vertices via stable grid snap.
    pos_dict = Dict{NTuple{2,Int},SVector{2,Float64}}()
    for tile in tiles
        for v in tile.vertices
            k = snap_to_grid(v, SNAP_GRID_EPS)
            get!(pos_dict, k, v)
        end
    end
    positions = collect(values(pos_dict))

    params = Dict{Symbol,Any}(
        :generations => generations,
        :n_tiles => length(tiles),
        :n_vertices => length(positions),
        :window_shape => :none,
    )
    return QuasicrystalData{2,Float64}(PenroseP3(), positions, tiles, method, params)
end

# ---- inflate_penrose_tiles: substitution on existing tile lists -------

"""
    inflate_penrose_tiles(tiles::Vector{Tile{2,Float64}},
                          alg::AbstractSubstitutionAlgorithm)
        → Vector{Tile{2,Float64}}

Apply one generation of Penrose P3 substitution to an existing
tile list. The implementation:

1. Decomposes every input rhombus into its two Robinson
   half-tiles ([`_split_rhombus_to_robinson`](@ref)).
2. Deflates the half-tile list ([`deflate_robinson`](@ref)).
3. Re-pairs deflated half-tiles into rhombi
   ([`merge_robinson_to_rhombi`](@ref)).

`DefaultSubstitution` and [`RobinsonTriangleInflation`](@ref) are
both fully wired through this routine; [`DirectTileInflation`](@ref)
is currently a not-yet-implemented placeholder.
"""
function inflate_penrose_tiles(tiles::Vector{Tile{2,Float64}}, ::RobinsonTriangleInflation)
    tris = RobTri[]
    sizehint!(tris, 2 * length(tiles))
    for tile in tiles
        append!(tris, _split_rhombus_to_robinson(tile))
    end
    deflated = deflate_robinson(tris)
    return merge_robinson_to_rhombi(deflated)
end

function inflate_penrose_tiles(tiles::Vector{Tile{2,Float64}}, ::DefaultSubstitution)
    # `DefaultSubstitution` is contractually the fastest fully-working
    # algorithm for the family — for Penrose that is the Robinson
    # triangle inflation route.
    return inflate_penrose_tiles(tiles, RobinsonTriangleInflation())
end

# `DirectTileInflation` for Penrose is a placeholder shell — the
# direct tile-by-tile substitution rules have not been ported yet.
# Make the unimplemented status explicit instead of silently
# delegating to a different algorithm.
function inflate_penrose_tiles(::Vector{Tile{2,Float64}}, ::DirectTileInflation)
    error(
        "DirectTileInflation is not yet implemented for PenroseP3. " *
        "Use DefaultSubstitution() or RobinsonTriangleInflation() instead.",
    )
end

# Single-dispatch on the algorithm: `RobinsonTriangleInflation` is
# Penrose-specific, so this overload is unambiguous.
function inflate_tiles(tiles::Vector{Tile{2,Float64}}, alg::RobinsonTriangleInflation)
    return inflate_penrose_tiles(tiles, alg)
end

"""
    _split_rhombus_to_robinson(tile::Tile{2,Float64}) → Vector{RobTri}

Split a Penrose rhombus into its two Robinson half-tiles. Splits a
fat rhombus along its long diagonal (between the two 72° corners)
into two `:blue` (obtuse) Robinson triangles, and a thin rhombus
along its short diagonal (between the two 144° corners) into two
`:red` (acute) Robinson triangles.

The tile's vertices are assumed to be ordered cyclically
(`v1, v2, v3, v4`). For a fat rhombus, the diagonal (`v1, v3`)
corresponds to the long diagonal because the canonical generators
emit fat rhombi with `v1, v3` at the two 72° corners (see
[`merge_robinson_to_rhombi`](@ref)). For a thin rhombus, the same
diagonal `(v1, v3)` corresponds to the short diagonal (between the
two 144° corners) for the same reason.
"""
function _split_rhombus_to_robinson(tile::Tile{2,Float64})
    v = tile.vertices
    length(v) == 4 || throw(ArgumentError("Penrose rhombus tiles must have 4 vertices"))
    v1, v2, v3, v4 = v[1], v[2], v[3], v[4]
    if tile.type isa FatRhombus
        # v1, v3 are the 72° corners (long-diagonal endpoints); v2, v4
        # are the 108° corners (apex of each blue half-tile). Both
        # halves get the **same** `(B, C) = (v1, v3)` labelling so
        # that [`subdivide_blue`](@ref) places its sub-tile cuts at
        # the same spatial point on both sides of the shared
        # diagonal — the two half-tiles still differ in apex parity
        # (`B → C` traversed CCW vs CW around `A`), but the deflation
        # rule is parity-agnostic.
        return RobTri[RobTri(:blue, v2, v1, v3), RobTri(:blue, v4, v1, v3)]
    elseif tile.type isa ThinRhombus
        # v1, v3 are the 144° corners (short-diagonal endpoints);
        # v2, v4 are the 36° corners (apex of each red half-tile).
        # Same `(B, C) = (v1, v3)` labelling rationale as above.
        return RobTri[RobTri(:red, v2, v1, v3), RobTri(:red, v4, v1, v3)]
    else
        throw(ArgumentError("unsupported tile type $(tile.type) for Penrose split"))
    end
end

export vertex_angles, vertex_configuration

"""
    vertex_angles(data::QuasicrystalData{2, Float64, PenroseP3}, v_idx::Int) → Vector{Int}

Compute the internal angles of the tiles surrounding the vertex `v_idx` in degrees.
Returns a sorted vector of integers representing the angles (e.g., `[72, 72, 72, 72, 72]`).
"""
function vertex_angles(data::QuasicrystalData{2,Float64,PenroseP3}, v_idx::Int)
    v_pos = data.positions[v_idx]

    # We find all tiles that contain v_pos
    _ensure_plaquettes!(data)
    plaqs = data.parameters[:plaquettes]

    angles = Int[]
    for p in plaqs
        if v_idx in p.vertices
            # Find its position in the cycle
            idx_in_p = findfirst(==(v_idx), p.vertices)
            N = length(p.vertices)
            prev_idx = p.vertices[mod1(idx_in_p - 1, N)]
            next_idx = p.vertices[mod1(idx_in_p + 1, N)]

            p_prev = data.positions[prev_idx]
            p_next = data.positions[next_idx]

            # Vectors from v_pos
            vec1 = p_prev - v_pos
            vec2 = p_next - v_pos

            # Compute angle
            cos_theta = dot(vec1, vec2) / (norm(vec1) * norm(vec2))
            cos_theta = clamp(cos_theta, -1.0, 1.0)
            angle_deg = round(Int, acos(cos_theta) * 180 / π)
            push!(angles, angle_deg)
        end
    end
    return sort(angles)
end

"""
    vertex_configuration(data::QuasicrystalData{2, Float64, PenroseP3}, v_idx::Int) → Symbol

Identify the vertex configuration around `v_idx` using the 8 standard P3 vertex types.
Standard combinations (sorted degrees summing to 360):
Returns the signature or standard name if known.
"""
function vertex_configuration(data::QuasicrystalData{2,Float64,PenroseP3}, v_idx::Int)
    angles = vertex_angles(data, v_idx)

    # The sum of angles must be 360 for interior vertices.
    if sum(angles) != 360
        return :Boundary  # Not fully surrounded
    end

    # Map from sorted angle lists to Conway/De Bruijn P3 identifiers
    if angles == [72, 72, 72, 72, 72]
        return :Sun
    elseif angles == [72, 144, 144]
        return :Star
    elseif angles == [36, 108, 108, 108]
        return :King
    elseif angles == [72, 72, 108, 108]
        return :Queen
    elseif angles == [36, 72, 108, 144]
        return :Jack
    elseif angles == [36, 36, 144, 144]
        return :Deuce
    elseif angles == [36, 36, 72, 108, 108]
        return :Ace
    else
        return Symbol("V_", join(angles, "_"))
    end
end
