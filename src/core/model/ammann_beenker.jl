"""
Ammann–Beenker tiling: 8-fold rotationally symmetric cut-and-project
quasicrystal built from squares and 45° rhombi, from a 4D hypercubic
host lattice.
"""

"""
    AmmannBeenker <: AbstractQuasicrystal{2}

Topology marker for the Ammann–Beenker tiling.
"""
struct AmmannBeenker <: AbstractQuasicrystal{2} end

const SQRT2 = sqrt(2.0)

"""
    generate_ammann_beenker_projection(radius::Real;
                                       method::ProjectionMethod = ProjectionMethod())
        → QuasicrystalData{2, Float64}

Generate an Ammann–Beenker point set by projecting `Z^4` onto the
2D physical subspace of the cut-and-project construction.
Accepts lattice points whose perpendicular projection falls inside
a 2D window of half-width `window_size`.
"""
function generate_ammann_beenker_projection(
    radius::Real; method::ProjectionMethod=ProjectionMethod()
)
    theta = π / 4

    E_par = zeros(4, 2)
    for i in 1:4
        E_par[i, 1] = cos((i - 1) * theta)
        E_par[i, 2] = sin((i - 1) * theta)
    end

    E_perp = zeros(4, 2)
    for i in 1:4
        E_perp[i, 1] = cos((i - 1) * theta + π / 4)
        E_perp[i, 2] = sin((i - 1) * theta + π / 4)
    end

    window_size = 0.5
    raw_positions = Vector{Float64}[]
    n_max = ceil(Int, radius * 1.5)

    for n1 in (-n_max):n_max,
        n2 in (-n_max):n_max, n3 in (-n_max):n_max,
        n4 in (-n_max):n_max

        lattice_point = [float(n1), float(n2), float(n3), float(n4)]
        pos_par = E_par' * lattice_point
        norm(pos_par) > radius && continue
        pos_perp = E_perp' * lattice_point
        if all(abs.(pos_perp) .<= window_size)
            push!(raw_positions, pos_par)
        end
    end

    positions = [SVector{2,Float64}(p[1], p[2]) for p in raw_positions]
    tiles = Tile{2,Float64}[]

    params = Dict{Symbol,Any}(
        :radius => radius,
        :n_max => n_max,
        :window_size => window_size,
        :n_vertices => length(positions),
        :symmetry => 8,
    )
    return QuasicrystalData{2,Float64}(AmmannBeenker(), positions, tiles, method, params)
end

"""
    generate_ammann_beenker_substitution(generations::Int;
                                         method::SubstitutionMethod = SubstitutionMethod())
        → QuasicrystalData{2, Float64}

Generate an Ammann–Beenker point set by iterating the substitution
(inflation) rules with an inflation factor of `1 + √2`.

The current inflation routine is a **placeholder** that scales
tiles rather than subdividing them properly into smaller squares
and rhombi. A full implementation is a future task.
"""
function generate_ammann_beenker_substitution(
    generations::Int; method::SubstitutionMethod=SubstitutionMethod()
)
    initial_tiles = Tile{2,Float64}[]

    # Seed: a unit square.
    square_vertices = [
        SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0), SVector(0.0, 1.0)
    ]
    square_center = SVector(0.5, 0.5)
    push!(initial_tiles, Tile{2,Float64}(square_vertices, 1, square_center))

    # Seed: 8 rhombi in an octagonal wreath.
    for i in 0:7
        angle = i * π / 4
        v1 = SVector(cos(angle), sin(angle))
        v2 = v1 + SVector(cos(angle + π / 4), sin(angle + π / 4))
        v3 = v2 + SVector(cos(angle + π), sin(angle + π))
        v4 = v1 + SVector(cos(angle + π), sin(angle + π))
        center = (v1 + v2 + v3 + v4) / 4
        push!(initial_tiles, Tile{2,Float64}([v1, v2, v3, v4], 2, center))
    end

    tiles = initial_tiles
    for _ in 1:generations
        tiles = inflate_ammann_beenker_tiles(tiles)
    end

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
        :symmetry => 8,
    )
    return QuasicrystalData{2,Float64}(AmmannBeenker(), positions, tiles, method, params)
end

"""
    inflate_ammann_beenker_tiles(tiles::Vector{Tile{2, Float64}})

Placeholder inflation dispatch: see `inflate_ab_square` and
`inflate_ab_rhombus`.
"""
function inflate_ammann_beenker_tiles(tiles::Vector{Tile{2,Float64}})
    new_tiles = Tile{2,Float64}[]
    inflation_factor = 1 + SQRT2
    for tile in tiles
        if tile.type == 1
            append!(new_tiles, inflate_ab_square(tile, inflation_factor))
        else
            append!(new_tiles, inflate_ab_rhombus(tile, inflation_factor))
        end
    end
    return new_tiles
end

function inflate_ab_square(tile::Tile{2,Float64}, factor::Float64)
    v = tile.vertices
    return [Tile{2,Float64}([factor * vertex for vertex in v], 1, factor * tile.center)]
end

function inflate_ab_rhombus(tile::Tile{2,Float64}, factor::Float64)
    v = tile.vertices
    return [Tile{2,Float64}([factor * vertex for vertex in v], 2, factor * tile.center)]
end
