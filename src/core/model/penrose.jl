"""
Penrose P3 (rhombus) tiling implementation: cut-and-project
quasicrystal with 10-fold (5-fold) rotational symmetry from a 5D
hypercubic host lattice.
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
    raw_positions = Vector{Float64}[]
    n_max = ceil(Int, radius * 1.5)

    for n1 in (-n_max):n_max,
        n2 in (-n_max):n_max, n3 in (-n_max):n_max, n4 in (-n_max):n_max,
        n5 in (-n_max):n_max

        lattice_point = [float(n1), float(n2), float(n3), float(n4), float(n5)]
        pos_par = E_par' * lattice_point

        norm(pos_par) > radius && continue

        pos_perp = E_perp' * lattice_point
        if all(abs.(pos_perp) .<= window_size)
            push!(raw_positions, pos_par)
        end
    end

    positions = [SVector{2,Float64}(p[1], p[2]) for p in raw_positions]
    tiles = Tile{2,Float64}[]   # Tile construction is a future task.

    params = Dict{Symbol,Any}(
        :radius => radius,
        :n_max => n_max,
        :window_size => window_size,
        :n_vertices => length(positions),
    )
    return QuasicrystalData{2,Float64}(PenroseP3(), positions, tiles, method, params)
end

"""
    generate_penrose_substitution(generations::Int;
                                  method::SubstitutionMethod = SubstitutionMethod())
        → QuasicrystalData{2, Float64}

Generate a Penrose P3 point set by substitution (inflation) rules.

The inflation routine currently implemented is a **placeholder**
that scales tiles by the golden ratio rather than subdividing them
properly into fat + thin rhombi. It reproduces the correct point
set only for zero generations; higher `generations` return a
self-consistently scaled tiling whose positions are *not* the true
Penrose tiling. This is a known pre-migration limitation and will
be fixed in a dedicated follow-up PR alongside the Bragg peak
enumeration work.
"""
function generate_penrose_substitution(
    generations::Int; method::SubstitutionMethod=SubstitutionMethod()
)
    # Start with a "Sun" of 5 fat rhombi
    # Each fat rhombus is spanned by (e_i, e_{i+1})
    star = [SVector(cos(i * 2π / 5), sin(i * 2π / 5)) for i in 0:4]

    current_rhombi = []
    for i in 0:4
        # (index1, index2, offset)
        push!(current_rhombi, (i, mod(i+1, 5), SVector(0.0, 0.0)))
    end

    for _ in 1:generations
        new_rhombi = []
        for (i, j, w) in current_rhombi
            # Grid substitution: e_i -> e_{i-1} + e_i + e_{i+1}
            # The tile (e_i, e_j) is replaced by 9 tiles (e_a, e_b)
            # a in {i-1, i, i+1}, b in {j-1, j, j+1}

            # Sub-vectors for i
            U = [star[mod(i - 1, 5) + 1], star[i + 1], star[mod(i + 1, 5) + 1]]
            # Sub-vectors for j
            V = [star[mod(j - 1, 5) + 1], star[j + 1], star[mod(j + 1, 5) + 1]]

            # Scale old offset
            w_scaled = w * ϕ

            for (ai, u) in enumerate(U), (bi, v) in enumerate(V)
                # Position of sub-tile (ai, bi)
                # offset is sum of previous vectors in the expansion
                pos = w_scaled
                for ak in 1:(ai - 1)
                    pos += U[ak]
                end
                for bk in 1:(bi - 1)
                    pos += V[bk]
                end

                # New indices
                idx_a = mod(i + (ai-2), 5)
                idx_b = mod(j + (bi-2), 5)

                if idx_a != idx_b
                    push!(new_rhombi, (idx_a, idx_b, pos))
                end
            end
        end
        current_rhombi = new_rhombi
    end

    # Convert to Tiles and deduplicate
    tile_dict = Dict{NTuple{2,Int},Tile{2,Float64}}()
    star_vectors = star
    for (i, j, w) in current_rhombi
        v1 = w
        v2 = w + star_vectors[i + 1]
        v3 = w + star_vectors[i + 1] + star_vectors[j + 1]
        v4 = w + star_vectors[j + 1]

        # Canonical key for deduplication: snap centre to a fixed grid
        # (stable hash key under floating-point round-off).
        center = (v1 + v3) / 2
        key = snap_to_grid(center, 1e-5)

        # Determine type: Fat if |i-j| == 1 or 4, Thin if |i-j| == 2 or 3
        diff = mod(abs(i - j), 5)
        type = (diff == 1 || diff == 4) ? 1 : 2

        if !haskey(tile_dict, key)
            tile_dict[key] = Tile{2,Float64}([v1, v2, v3, v4], type, center)
        end
    end

    tiles = collect(values(tile_dict))

    # Collect unique vertices via stable grid snap
    pos_dict = Dict{NTuple{2,Int},SVector{2,Float64}}()
    for tile in tiles
        for v in tile.vertices
            k = snap_to_grid(v, 1e-5)
            get!(pos_dict, k, v)
        end
    end
    positions = collect(values(pos_dict))

    params = Dict{Symbol,Any}(
        :generations => generations,
        :n_tiles => length(tiles),
        :n_vertices => length(positions),
    )
    return QuasicrystalData{2,Float64}(PenroseP3(), positions, tiles, method, params)
end

"""
    inflate_penrose_tiles(tiles::Vector{Tile{2, Float64}}, alg::AbstractSubstitutionAlgorithm)

Apply the Penrose substitution rules tile-by-tile.
"""
struct RobTri
    type::Int
    parity::Int # 1 for Left, -1 for Right
    a::SVector{2,Float64}
    b::SVector{2,Float64}
    c::SVector{2,Float64}
end

function inflate_penrose_tiles(
    tiles::Vector{Tile{2,Float64}}, alg::RobinsonTriangleInflation
)
    # For now, we reuse the robust logic by identifying star indices from tiles
    star = [SVector(cos(i * 2π / 5), sin(i * 2π / 5)) for i in 0:4]

    current_rhombi = []
    for tile in tiles
        v = tile.vertices
        # Identify spanning vectors from origin (v1)
        v1 = v[1]
        e1 = v[2] - v1
        e2 = v[4] - v1

        # Find closest star indices
        i = -1
        j = -1
        for k in 0:4
            if norm(e1 - star[k + 1]) < 1e-4
                i = k
            elseif norm(e1 + star[k + 1]) < 1e-4
                # Handle mirrored/inverted tiles if necessary
                # In P3 they are usually aligned to star
            end
            if norm(e2 - star[k + 1]) < 1e-4
                j = k
            end
        end

        if i != -1 && j != -1
            push!(current_rhombi, (i, j, v1))
        end
    end

    # Apply one generation of vector inflation
    new_rhombi = []
    for (i, j, w) in current_rhombi
        U = [star[mod(i - 1, 5) + 1], star[i + 1], star[mod(i + 1, 5) + 1]]
        V = [star[mod(j - 1, 5) + 1], star[j + 1], star[mod(j + 1, 5) + 1]]
        w_scaled = w * ϕ
        for (ai, u) in enumerate(U), (bi, v) in enumerate(V)
            pos = w_scaled
            for ak in 1:(ai - 1)
                pos += U[ak]
            end
            for bk in 1:(bi - 1)
                pos += V[bk]
            end
            idx_a = mod(i + (ai-2), 5)
            idx_b = mod(j + (bi-2), 5)
            if idx_a != idx_b
                push!(new_rhombi, (idx_a, idx_b, pos))
            end
        end
    end

    # Deduplicate and return Tiles
    tile_dict = Dict{Tuple{Int,Int},Tile{2,Float64}}()
    for (i, j, w) in new_rhombi
        v1, v2, v3, v4 = w, w + star[i + 1], w + star[i + 1] + star[j + 1], w + star[j + 1]
        c = (v1 + v3) / 2
        key = (round(Int, c[1]*1e5), round(Int, c[2]*1e5))
        diff = mod(abs(i - j), 5)
        type = (diff == 1 || diff == 4) ? 1 : 2
        if !haskey(tile_dict, key)
            tile_dict[key] = Tile{2,Float64}([v1, v2, v3, v4], type, c)
        end
    end
    return collect(values(tile_dict))
end

# Keep internal helpers for backward compatibility if needed, but they are no longer used by main loop

function inflate_penrose_tiles(tiles::Vector{Tile{2,Float64}}, alg::DefaultSubstitution)
    inflate_penrose_tiles(tiles, RobinsonTriangleInflation())
end

# Placeholder for direct tile inflation
function inflate_penrose_tiles(tiles::Vector{Tile{2,Float64}}, alg::DirectTileInflation)
    # Fallback to Robinson for now
    return inflate_penrose_tiles(tiles, RobinsonTriangleInflation())
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
