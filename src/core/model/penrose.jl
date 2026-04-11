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

    initial_tiles = Tile{2,Float64}[]
    for i in 0:4
        angle = i * 2π / 5
        v1 = SVector(0.0, 0.0)
        v2 = SVector(cos(angle), sin(angle))
        v3 = v2 + SVector(cos(angle + angle_fat), sin(angle + angle_fat))
        v4 = SVector(cos(angle + angle_fat), sin(angle + angle_fat))
        center = (v1 + v2 + v3 + v4) / 4
        push!(initial_tiles, Tile{2,Float64}([v1, v2, v3, v4], 1, center))
    end

    tiles = initial_tiles
    for _ in 1:generations
        tiles = inflate_penrose_tiles(tiles)
    end

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
    inflate_penrose_tiles(tiles::Vector{Tile{2, Float64}})

Apply the (placeholder) Penrose substitution rules tile-by-tile.
See the note in [`generate_penrose_substitution`](@ref) — the
current implementation scales rather than subdivides.
"""
function inflate_penrose_tiles(tiles::Vector{Tile{2,Float64}})
    new_tiles = Tile{2,Float64}[]
    for tile in tiles
        if tile.type == 1
            append!(new_tiles, inflate_fat_rhombus(tile))
        else
            append!(new_tiles, inflate_thin_rhombus(tile))
        end
    end
    return new_tiles
end

function inflate_fat_rhombus(tile::Tile{2,Float64})
    # TODO: proper inflation (1 fat → 1 fat + 2 thin).
    v = tile.vertices
    return [Tile{2,Float64}([ϕ * vertex for vertex in v], 1, ϕ * tile.center)]
end

function inflate_thin_rhombus(tile::Tile{2,Float64})
    # TODO: proper inflation (1 thin → 1 fat).
    v = tile.vertices
    return [Tile{2,Float64}([ϕ * vertex for vertex in v], 2, ϕ * tile.center)]
end
