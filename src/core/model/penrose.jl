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
    angle_fat = deg2rad(72)

    initial_tris = RobTri[]
    for i in 0:4
        angle = i * 2π / 5
        v1 = SVector(0.0, 0.0)
        v2 = SVector(cos(angle), sin(angle))
        v3 = v2 + SVector(cos(angle + angle_fat), sin(angle + angle_fat))
        v4 = SVector(cos(angle + angle_fat), sin(angle + angle_fat))
        
        # A Fat rhombus (72, 108) consists of two Obtuse (type 2) triangles
        # sharing their long edge (v1-v3 in our setup).
        push!(initial_tris, RobTri(2, 1, v2, v3, v1))
        push!(initial_tris, RobTri(2, -1, v4, v1, v3))
    end

    triangles = initial_tris
    for _ in 1:generations
        triangles = _inflate_rob_tris(triangles)
    end
    
    # Pair triangles into rhombi at the very end to avoid erosion
    tiles = _pair_rob_tris(triangles)

    # Collect unique vertices from every tile.
    position_set = Set{SVector{2,Float64}}()
    for tile in tiles
        for v in tile.vertices
            push!(position_set, v)
        end
    end
    positions = collect(position_set)

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
    a::SVector{2, Float64}
    b::SVector{2, Float64}
    c::SVector{2, Float64}
end

function inflate_penrose_tiles(tiles::Vector{Tile{2,Float64}}, alg::RobinsonTriangleInflation)
    # Decompose into Robinson Triangles
    triangles = RobTri[]
    for tile in tiles
        v = tile.vertices
        if tile.type == 1 # Fat rhombus (108, 72)
            # Find the acute angles (72). They are at opposite vertices.
            # Long diagonal connects the 72-degree vertices.
            diagonal_sq1 = sum(abs2, v[1] - v[3])
            diagonal_sq2 = sum(abs2, v[2] - v[4])
            if diagonal_sq1 > diagonal_sq2
                push!(triangles, RobTri(2, 1, v[2], v[3], v[1]))
                push!(triangles, RobTri(2, -1, v[4], v[1], v[3]))
            else
                push!(triangles, RobTri(2, 1, v[1], v[2], v[4]))
                push!(triangles, RobTri(2, -1, v[3], v[4], v[2]))
            end
        else # Thin rhombus (144, 36)
            diagonal_sq1 = sum(abs2, v[1] - v[3])
            diagonal_sq2 = sum(abs2, v[2] - v[4])
            if diagonal_sq1 < diagonal_sq2
                push!(triangles, RobTri(1, 1, v[2], v[3], v[1]))
                push!(triangles, RobTri(1, -1, v[4], v[1], v[3]))
            else
                push!(triangles, RobTri(1, 1, v[1], v[2], v[4]))
                push!(triangles, RobTri(1, -1, v[3], v[4], v[2]))
            end
        end
    end
    
    inflated = _inflate_rob_tris(triangles)
    return _pair_rob_tris(inflated)
end

"""
    _inflate_rob_tris(triangles::Vector{RobTri}) → Vector{RobTri}

Internal helper to deflate one set of Robinson triangles and scale up by ϕ.
Preserves all fragments (erosion-free).
"""
function _inflate_rob_tris(triangles::Vector{RobTri})
    new_tris = RobTri[]
    for t in triangles
        a, b, c = t.a, t.b, t.c
        if t.type == 1 # Acute
            if t.parity == 1
                p = a + (b - a) / ϕ
                push!(new_tris, RobTri(1, -1, c, p, b))
                push!(new_tris, RobTri(2, 1, p, c, a))
            else
                p = a + (c - a) / ϕ
                push!(new_tris, RobTri(1, 1, b, p, c))
                push!(new_tris, RobTri(2, -1, p, b, a))
            end
        else # Obtuse
            if t.parity == 1
                p = b + (c - b) / ϕ
                push!(new_tris, RobTri(1, 1, b, p, a))
                push!(new_tris, RobTri(2, -1, p, a, c))
            else
                p = c + (b - c) / ϕ
                push!(new_tris, RobTri(1, -1, c, p, a))
                push!(new_tris, RobTri(2, 1, p, a, b))
            end
        end
    end
    
    # Scale up by phi
    scaled = RobTri[]
    for t in new_tris
        push!(scaled, RobTri(t.type, t.parity, t.a * ϕ, t.b * ϕ, t.c * ϕ))
    end
    return scaled
end

"""
    _pair_rob_tris(triangles::Vector{RobTri}) → Vector{Tile}

Internal helper to pair compatible Robinson triangles into fat/thin rhombi.
Triangles without a neighbor on their base (BC) are dropped.
"""
function _pair_rob_tris(triangles::Vector{RobTri})
    new_tiles = Tile{2, Float64}[]
    edge_dict = Dict{Tuple{Float64, Float64}, Vector{RobTri}}()
    
    for t in triangles
        mid = (t.b + t.c) / 2
        # Use rounding for stable dictionary keys
        key = (round(mid[1], digits=5), round(mid[2], digits=5))
        if !haskey(edge_dict, key)
            edge_dict[key] = [t]
        else
            push!(edge_dict[key], t)
        end
    end
    
    for (key, list) in edge_dict
        if length(list) == 2
            t1, t2 = list[1], list[2]
            if t1.type == 1 && t2.type == 1
                # Thin rhombus (two Acute triangles)
                v = [t1.a, t1.b, t2.a, t1.c]
                center = sum(v) / 4.0
                push!(new_tiles, Tile{2, Float64}(v, 2, center))  # type 2 is Thin
            elseif t1.type == 2 && t2.type == 2
                # Fat rhombus (two Obtuse triangles)
                v = [t1.a, t1.b, t2.a, t1.c]
                center = sum(v) / 4.0
                push!(new_tiles, Tile{2, Float64}(v, 1, center))  # type 1 is Fat
            end
        end
    end
    return new_tiles
end

inflate_penrose_tiles(tiles::Vector{Tile{2,Float64}}, alg::DefaultSubstitution) = inflate_penrose_tiles(tiles, RobinsonTriangleInflation())

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
function vertex_angles(data::QuasicrystalData{2, Float64, PenroseP3}, v_idx::Int)
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
function vertex_configuration(data::QuasicrystalData{2, Float64, PenroseP3}, v_idx::Int)
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


