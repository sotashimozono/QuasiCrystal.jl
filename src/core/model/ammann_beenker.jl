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
    # Start with an octagonal wreath of 8 rhombi
    star = [SVector(cos(i * π / 4), sin(i * π / 4)) for i in 0:7]

    current_rhombi = []
    for i in 0:7
        # Tile spanned by (e_i, e_{i+1})
        push!(current_rhombi, (i, mod(i+1, 8), SVector(0.0, 0.0)))
    end

    λ = 1 + sqrt(2)

    for _ in 1:generations
        new_rhombi = []
        for (i, j, w) in current_rhombi
            # Grid substitution: e_i -> e_{i-1} + e_i + e_{i+1}
            # The tile (e_i, e_j) is replaced by 9 tiles (e_a, e_b)
            # a in {i-1, i, i+1}, b in {j-1, j, j+1}

            U = [star[mod(i - 1, 8) + 1], star[i + 1], star[mod(i + 1, 8) + 1]]
            V = [star[mod(j - 1, 8) + 1], star[j + 1], star[mod(j + 1, 8) + 1]]

            w_scaled = w * λ

            for (ai, u) in enumerate(U), (bi, v) in enumerate(V)
                pos = w_scaled
                for ak in 1:(ai - 1)
                    pos += U[ak]
                end
                for bk in 1:(bi - 1)
                    pos += V[bk]
                end

                idx_a = mod(i + (ai-2), 8)
                idx_b = mod(j + (bi-2), 8)

                if idx_a != idx_b && mod(abs(idx_a - idx_b), 8) != 4
                    push!(new_rhombi, (idx_a, idx_b, pos))
                end
            end
        end
        current_rhombi = new_rhombi
    end

    # Convert to Tiles and deduplicate
    tile_dict = Dict{NTuple{2,Int},Tile{2,Float64}}()
    for (i, j, w) in current_rhombi
        v1, v2, v3, v4 = w, w + star[i + 1], w + star[i + 1] + star[j + 1], w + star[j + 1]
        c = (v1 + v3) / 2
        key = snap_to_grid(c, 1e-5)

        # Determine type: Square if |i-j| == 2 or 6, Rhombus otherwise.
        diff = mod(abs(i - j), 8)
        type = (diff == 2 || diff == 6) ? Square() : RhombusAB()

        if !haskey(tile_dict, key)
            tile_dict[key] = Tile{2,Float64}([v1, v2, v3, v4], type, c)
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
        :symmetry => 8,
    )
    return QuasicrystalData{2,Float64}(AmmannBeenker(), positions, tiles, method, params)
end

"""
    inflate_ammann_beenker_tiles(tiles::Vector{Tile{2, Float64}})

Inflation rule for AB: λe₁ = e₁₋₁ + e₁ + e₁₊₁ where λ = 1+√2
"""
function inflate_ammann_beenker_tiles(tiles::Vector{Tile{2,Float64}})
    # Robust rule for AB: λe₁ = e₁₋₁ + e₁ + e₁₊₁ where λ = 1+√2
    star = [SVector(cos(i * π / 4), sin(i * π / 4)) for i in 0:7]
    λ = 1 + sqrt(2)

    current_rhombi = []
    for tile in tiles
        v = tile.vertices
        v1 = v[1]
        e1, e2 = v[2] - v1, v[4] - v1
        i = j = -1
        for k in 0:7
            if norm(e1 - star[k + 1]) < 1e-4
                i = k
            end
            if norm(e2 - star[k + 1]) < 1e-4
                j = k
            end
        end
        if i != -1 && j != -1
            push!(current_rhombi, (i, j, v1))
        end
    end

    new_rhombi = []
    for (i, j, w) in current_rhombi
        U = [star[mod(i - 1, 8) + 1], star[i + 1], star[mod(i + 1, 8) + 1]]
        V = [star[mod(j - 1, 8) + 1], star[j + 1], star[mod(j + 1, 8) + 1]]
        w_scaled = w * λ
        for (ai, u) in enumerate(U), (bi, v) in enumerate(V)
            pos = w_scaled
            for ak in 1:(ai - 1)
                pos += U[ak]
            end
            for bk in 1:(bi - 1)
                pos += V[bk]
            end
            idx_a, idx_b = mod(i + (ai-2), 8), mod(j + (bi-2), 8)
            if idx_a != idx_b && mod(abs(idx_a - idx_b), 8) != 4
                push!(new_rhombi, (idx_a, idx_b, pos))
            end
        end
    end

    tile_dict = Dict{Tuple{Int,Int},Tile{2,Float64}}()
    for (i, j, w) in new_rhombi
        v1, v2, v3, v4 = w, w + star[i + 1], w + star[i + 1] + star[j + 1], w + star[j + 1]
        c = (v1 + v3) / 2
        key = (round(Int, c[1]*1e5), round(Int, c[2]*1e5))
        diff = mod(abs(i - j), 8)
        type = (diff == 2 || diff == 6) ? Square() : RhombusAB()
        if !haskey(tile_dict, key)
            tile_dict[key] = Tile{2,Float64}([v1, v2, v3, v4], type, c)
        end
    end
    return collect(values(tile_dict))
end
